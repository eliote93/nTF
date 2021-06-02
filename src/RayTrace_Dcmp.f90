#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTrace_Dcmp(RayInfo, CoreInfo, iz, gb, ge, lJout, lSubGrp, lScat1, lLinSrcCASMO, lHybrid)

USE OMP_LIB
USE PARAM,       ONLY : ZERO, RED, BLACK, GREEN
USE TYPEDEF,     ONLY : RayInfo_Type, Coreinfo_type, Asy_Type, AsyInfo_Type
USE Moc_Mod,     ONLY : TrackingDat, phisnm, phimnm, srcnm, xstnm, MocJoutnm, PhiAngInnm, DcmpPhiAngIn, DcmpPhiAngOut, &
                        RayTraceDcmp_NM, RayTraceDcmp_Pn, RayTraceDcmp_LSCASMO, DcmpGatherBoundaryFlux, DcmpScatterBoundaryFlux, DcmpLinkBoundaryFlux
USE Core_mod,    ONLY : phisSlope, srcSlope
USE PE_MOD,      ONLY : PE
USE CNTL,        ONLY : nTracerCntl
USE itrcntl_mod, ONLY : itrcntl

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

INTEGER :: iz, gb, ge
LOGICAL :: lJout, lSubGrp, lScat1, lLinSrcCASMO, lHybrid
! ----------------------------------------------------
TYPE (AsyInfo_Type), POINTER, DIMENSION(:) :: AsyInfo
TYPE (Asy_Type),     POINTER, DIMENSION(:) :: Asy

INTEGER :: color, ithr, nThr, nPolarAngle, nAziAngle, ScatOd, AsyType, iAsy, startColor, endColor, colorInc
! ----------------------------------------------------

nPolarAngle = RayInfo%nPolarAngle
nAziAngle   = RayInfo%nAziAngle

AsyInfo => CoreInfo%AsyInfo
Asy     => CoreInfo%Asy

ScatOd = nTracerCntl%ScatOd

IF (.NOT. nTracerCntl%lHex) THEN
  IF (mod(itrcntl%mocit, 2) .EQ. 0) THEN
    startColor = BLACK
    endColor   = GREEN
    colorInc   = (GREEN - BLACK)/2
  ELSE
    startColor = GREEN
    endColor   = BLACK
    colorInc   = (BLACK - GREEN)/2
  END IF
ELSE
  IF (mod(itrcntl%mocit, 2) .EQ. 0) THEN
    startColor = BLACK
    endColor   = RED
    colorInc   = RED - BLACK
  ELSE
    startColor = RED
    endColor   = BLACK
    colorInc   = BLACK - RED
  END IF
END IF

CALL OMP_SET_NUM_THREADS(nThr)
! ----------------------------------------------------
DO ithr = 1, nThr
  TrackingDat(ithr)%srcnm => srcnm
  TrackingDat(ithr)%xstnm => xstnm
    
  IF (lLinSrcCASMO) TrackingDat(ithr)%srcSlope => srcSlope(:, :, :, iz)
END DO

DcmpPhiAngOut(:, gb:ge, :, :, :) = ZERO
! ----------------------------------------------------
DO color = startColor, endColor, colorInc
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) CALL DcmpScatterBoundaryFlux(RayInfo, PhiAngInnm, DcmpPhiAngIn)
#endif

  !$OMP PARALLEL PRIVATE(iAsy, AsyType)
  ithr = 1
  !$ ithr = omp_get_thread_num()+1
  
  TrackingDat(ithr)%PhiAngInnm    => PhiAngInnm
  TrackingDat(ithr)%DcmpPhiAngIn  => DcmpPhiAngIn
  TrackingDat(ithr)%DcmpPhiAngOut => DcmpPhiAngOut
  !$OMP BARRIER
  !$OMP DO SCHEDULE(DYNAMIC)
  DO iAsy = PE%myAsyBeg, PE%myAsyEnd
    IF (Asy(iAsy)%color .NE. color) CYCLE
    
    AsyType = Asy(iAsy)%AsyType
    
    IF (lScat1 .AND. .NOT.lSubGrp) THEN
      CALL RayTraceDcmp_Pn(RayInfo, CoreInfo, phisnm, phimnm, MocJoutnm, iz, iAsy, gb, ge, ScatOd, lJout)
    ELSE
      IF (.NOT.lLinSrcCASMO .OR. (AsyInfo(AsyType)%lFuel.AND.lHybrid)) THEN
        CALL RayTraceDcmp_NM(RayInfo, CoreInfo, phisnm, MocJoutnm, iz, iAsy, gb, ge, ljout)
      ELSE
        CALL RayTraceDcmp_LSCASMO(RayInfo, CoreInfo, phisnm, phisSlope, MocJoutnm, iz, iAsy, gb, ge, ljout)
      END IF
    END IF
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
  
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) CALL DcmpGatherBoundaryFlux(RayInfo, DcmpPhiAngOut)
#endif
  
  IF (PE%RTMASTER) CALL DcmpLinkBoundaryFlux(CoreInfo, RayInfo, PhiAngInnm, DcmpPhiAngIn, DcmpPhiAngOut, gb, ge, color)
END DO

NULLIFY (Asy)
NULLIFY (AsyInfo)
! ----------------------------------------------------

END SUBROUTINE RayTrace_Dcmp
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmp_NM(RayInfo, CoreInfo, phisnm, joutnm, iz, iasy, gb, ge, ljout)

USE OMP_LIB
USE PARAM,   ONLY : ZERO
USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, AsyInfo_Type, Asy_Type, Pin_Type, Cell_Type, DcmpAsyRayInfo_Type
USE Moc_Mod, ONLY : RecTrackRotRayNM_Dcmp, HexTrackRotRayNM_Dcmp, TrackingDat
USE CNTL,    ONLY : nTracerCntl

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:,:)     :: phisnm
REAL, POINTER, DIMENSION(:,:,:,:) :: joutnm

INTEGER :: iz, iAsy, gb, ge
LOGICAL :: ljout
! ----------------------------------------------------
TYPE (AsyInfo_Type),        POINTER, DIMENSION(:)   :: AsyInfo
TYPE (Asy_Type),            POINTER, DIMENSION(:)   :: Asy
TYPE (Cell_Type),           POINTER, DIMENSION(:)   :: Cell
TYPE (Pin_Type),            POINTER, DIMENSION(:)   :: Pin
TYPE (DcmpAsyRayInfo_Type), POINTER, DIMENSION(:,:) :: DcmpAsyRay

INTEGER, POINTER, DIMENSION(:) :: DcmpAsyRayCount

REAL, POINTER, DIMENSION(:,:) :: xstnm, srcnm

INTEGER :: nxy, iAsyRay, ithr, icel, ig, ipin, ifsr, jfsr, krot, AsyType, FsrIdxSt, FsrIdxEnd, PinIdxSt, PinIdxEnd
! ----------------------------------------------------

AsyInfo => CoreInfo%AsyInfo
Asy     => CoreInfo%Asy
Cell    => CoreInfo%CellInfo
Pin     => CoreInfo%Pin

AsyType   = Asy(iAsy)%AsyType
PinIdxSt  = Asy(iAsy)%GlobalPinIdx(1)
PinIdxEnd = Asy(iAsy)%GlobalPinIdx(nxy)

nxy = AsyInfo(AsyType)%nxy

FsrIdxSt  = Pin(PinIdxSt)%FsrIdxSt
FsrIdxEnd = Pin(PinIdxEnd)%FsrIdxSt + Cell(Pin(PinIdxEnd)%Cell(iz))%nFsr - 1

DcmpAsyRay      => RayInfo%DcmpAsyRay
DcmpAsyRayCount => RayInfo%DcmpAsyRayCount

ithr = omp_get_thread_num() + 1

xstnm => TrackingDat(ithr)%xstnm
srcnm => TrackingDat(ithr)%srcnm

phisnm(gb:ge, FsrIdxSt:FsrIdxEnd) = ZERO

IF (ljout) joutnm(:, gb:ge, :, PinIdxSt:PinIdxEnd) = ZERO
! ----------------------------------------------------
IF (nTracerCntl%lHex) THEN
  DO iAsyRay = 1, DcmpAsyRayCount(iAsy)
    DO krot = 1, 2
      CALL HexTrackRotRayNM_Dcmp(RayInfo, CoreInfo, TrackingDat(ithr), phisnm, joutnm, ljout, DcmpAsyRay(iAsyRay, iAsy), iz, gb, ge, krot)
    END DO
  END DO
ELSE
  DO iAsyRay = 1, DcmpAsyRayCount(iAsy)
    DO krot = 1, 2
      CALL RecTrackRotRayNM_Dcmp(RayInfo, CoreInfo, TrackingDat(ithr), phisnm, joutnm, ljout, DcmpAsyRay(iAsyRay, iAsy), iz, gb, ge, krot)
    END DO
  END DO
END IF

DO ipin = PinIdxSt, PinIdxEnd
  FsrIdxSt = Pin(ipin)%FsrIdxSt
  icel     = Pin(ipin)%Cell(iz)
  
  DO ifsr = 1, Cell(icel)%nFsr
    jfsr = FsrIdxSt + ifsr - 1
    
    DO ig = gb, ge
      phisnm(ig, jfsr) = phisnm(ig, jfsr) / (xstnm(ig, jfsr) * Cell(icel)%vol(ifsr)) + srcnm(ig, jfsr)
    END DO
  END DO
END DO
! ----------------------------------------------------
NULLIFY (AsyInfo)
NULLIFY (Asy)
NULLIFY (Cell)
NULLIFY (Pin)
NULLIFY (DcmpAsyRay)
NULLIFY (DcmpAsyRayCount)
NULLIFY (xstnm)
NULLIFY (srcnm)
! ----------------------------------------------------

END SUBROUTINE RayTraceDcmp_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayNM_Dcmp(RayInfo, CoreInfo, TrackingDat, phis, jout, ljout, DcmpAsyRay, iz, gb, ge, krot)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, Asy_Type, Cell_Type, AsyRayInfo_type, CellRayInfo_type, TrackingDat_Type, DcmpAsyRayInfo_Type

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

REAL, POINTER, DIMENSION(:,:)     :: phis
REAL, POINTER, DIMENSION(:,:,:,:) :: jout

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: iz, gb, ge, krot
! ----------------------------------------------------
TYPE (Pin_Type),        POINTER, DIMENSION(:) :: Pin
TYPE (Asy_Type),        POINTER, DIMENSION(:) :: Asy
TYPE (Cell_Type),       POINTER, DIMENSION(:) :: Cell
TYPE (AsyRayInfo_type), POINTER, DIMENSION(:) :: AsyRay

TYPE (CellRayInfo_Type),  POINTER :: CellRay

REAL, POINTER, DIMENSION(:)         :: LenSeg
REAL, POINTER, DIMENSION(:,:)       :: src, xst, EXPA, EXPB, wtang
REAL, POINTER, DIMENSION(:,:,:)     :: PhiAngIn, wtsurf
REAL, POINTER, DIMENSION(:,:,:,:,:) :: DcmpPhiAngIn, DcmpPhiAngOut

INTEGER, POINTER, DIMENSION(:) :: LocalFsrIdx, AsyRayList, DirList, AziList

INTEGER :: mp(2)
INTEGER :: iazi, ipol, irotray, iasyray, iceray, iray, irayseg, idir, ipin, icel, ibcel, iasy, ireg, isurf, ig, jbeg, jend, jinc, imray, ipray
INTEGER :: ipst, iped, ipinc, isfst, isfed, isgst, isged, isginc
INTEGER :: nAsyRay, nPinRay, nRaySeg, FsrIdxSt, nAziAng, nPolarAng, PhiAnginSvIdx, ExpAppIdx

REAL :: wt(10), wt2(10, 4)

REAL, DIMENSION(RayInfo%nPolarAngle, gb:ge) :: PhiAngOut
REAL :: phid, tau, ExpApp

DATA mp /2, 1/
! ----------------------------------------------------

! Ray Info.
nPolarAng = RayInfo%nPolarAngle
nAziAng   = RayInfo%nAziAngle
AsyRay   => RayInfo%AsyRay

! Geo.
Asy  => CoreInfo%Asy
Pin  => CoreInfo%Pin
Cell => CoreInfo%CellInfo

! Tracking Dat Pointing
src           => TrackingDat%srcnm
xst           => TrackingDat%xstnm
PhiAngIn      => TrackingDat%phiAngInnm
DcmpPhiAngIn  => TrackingDat%DcmpPhiAngIn
DcmpPhiAngOut => TrackingDat%DcmpPhiAngOut
EXPA          => TrackingDat%EXPA
EXPB          => TrackingDat%EXPB
wtang         => TrackingDat%wtang
wtsurf        => TrackingDat%wtsurf

! Dcmp. Ray
nAsyRay     = DcmpAsyRay%nAsyRay
iRotRay     = DcmpAsyRay%iRotRay
iAsy        = DcmpAsyRay%iAsy
iRay        = DcmpAsyRay%iRay
AsyRayList => DcmpAsyRay%AsyRayList
DirList    => DcmpAsyRay%DirList
AziList    => DcmpAsyRay%AziList

! Ray. B.C.
PhiAngInSvIdx = RayInfo%PhiAngInSvIdx(iRotRay, krot)

IF (DcmpAsyRay%lRotRayBeg(krot)) THEN
  PhiAngOut = PhiAngIn(1:nPolarAng, gb:ge, PhiAnginSvIdx)
ELSE
  PhiAngOut = DcmpPhiAngIn(1:nPolarAng, gb:ge, krot, iRay, iAsy)
END IF
  
IF (krot .EQ. 2) THEN
  jbeg = nAsyRay; jend = 1; jinc = -1
ELSE
  jbeg = 1; jend = nAsyRay; jinc = 1
END IF
! ----------------------------------------------------  
DO imray = jbeg, jend, jinc
  iAsyRay = AsyRayList(imray)
  nPinRay = AsyRay(iAsyRay)%nCellRay
  iazi    = AziList(imray)
  idir    = DirList(imray)
  IF (krot .EQ. 2) idir = mp(idir)
  
  DO ipol = 1, nPolarAng
    wt(ipol) = wtang(ipol, iazi)
    
    IF (lJout) wt2(ipol, 1:4) = wtsurf(ipol, iazi, 1:4)
  END DO
  
  IF (idir .EQ. 1) THEN
    ipst = 1; iped = nPinRay; ipinc = 1;  isfst = 1; isfed = 2;
  ELSE
    ipst = nPinRay; iped = 1; ipinc = -1; isfst = 2; isfed = 1;
  END IF
  
  DO ipray = ipst, iped, ipinc
    ipin   = AsyRay(iAsyRay)%PinIdx(ipray)
    iceray = AsyRay(iAsyRay)%PinRayIdx(ipray)
    
    ipin = Asy(iAsy)%GlobalPinIdx(ipin)   
    
    icel     = Pin(ipin)%Cell(iz)
    FsrIdxSt = Pin(ipin)%FsrIdxSt
    
    ibcel    = Cell(icel)%basecellstr
    CellRay => Cell(ibcel)%CellRay(iceray)
    
    nRaySeg      = CellRay%nSeg
    LocalFsrIdx => CellRay%LocalFsrIdx
    LenSeg      => CellRay%LenSeg
    
    IF (lJout) THEN
      isurf = AsyRay(iAsyRay)%PinRaySurf(isfst, ipray)
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          Jout(1, ig, isurf, ipin) = Jout(1, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
          Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
        END DO
      END DO
    END IF
    
    IF (idir .EQ. 1) THEN
      isgst = 1; isged = nRaySeg; isginc = 1
    ELSE
      isgst = nRaySeg; isged = 1; isginc = -1
    END IF
    
    DO iRaySeg = isgst, isged, isginc
      ireg = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
      
      DO ig = gb, ge
        tau = - LenSeg(iRaySeg) * xst(ig, ireg)
        
        ExpAppIdx = max(INT(tau), -40000)
        ExpAppIdx = min(0, ExpAppIdx)
        
        DO ipol = 1, nPolarAng
          ExpApp = EXPA(ipol, ExpAppIdx) * tau + EXPB(ipol, ExpAppIdx)
          
          phid = (PhiAngOut(ipol, ig) - src(ig, ireg)) * ExpApp
          
          PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
          
          phis(ig, ireg) = phis(ig, ireg) + wt(ipol) * phid
        END DO
      END DO
    END DO
    
    IF (lJout) THEN
      isurf = AsyRay(iAsyRay)%PinRaySurf(isfed, ipray)
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          Jout(2, ig, isurf, ipin) = Jout(2, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
          Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
        END DO
      END DO
    END IF
  END DO
END DO

DcmpPhiAngOut(1:nPolarAng, gb:ge, krot, iRay, iAsy) = PhiAngOut
! ----------------------------------------------------
! Geo.
NULLIFY (Pin)
NULLIFY (Asy)
NULLIFY (Cell)

! Ray
NULLIFY (AsyRay)
NULLIFY (CellRay)

! Loc.
NULLIFY (LenSeg)
NULLIFY (src)
NULLIFY (xst)
NULLIFY (EXPA)
NULLIFY (EXPB)
NULLIFY (wtang)
NULLIFY (wtsurf)
NULLIFY (PhiAngIn)
NULLIFY (LocalFsrIdx)

! Dcmp.
NULLIFY (DcmpPhiAngIn)
NULLIFY (DcmpPhiAngOut)
NULLIFY (AsyRayList)
NULLIFY (DirList)
NULLIFY (AziList)
! ----------------------------------------------------

END SUBROUTINE RecTrackRotRayNM_Dcmp
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayNM_Dcmp(RayInfo, CoreInfo, TrackingDat, phis, jout, ljout, DcmpAsyRay, iz, gb, ge, krot)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, TrackingDat_Type, DcmpAsyRayInfo_Type
USE HexData, ONLY : hAsy
USE HexType, ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexCoreRay, Type_HexRotRay
USE HexData, ONLY : haRay, hAsyTypInfo

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

REAL, POINTER, DIMENSION(:,:)     :: phis
REAL, POINTER, DIMENSION(:,:,:,:) :: jout

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: iz, gb, ge, krot
! ----------------------------------------------------
TYPE (Pin_Type), POINTER, DIMENSION(:) :: Pin

TYPE (Type_HexAsyRay), POINTER :: haRay_Loc
TYPE (Type_HexCelRay), POINTER :: CelRay_Loc

REAL, POINTER, DIMENSION(:,:)       :: src, xst, EXPA, EXPB, wtang, wthcs, wthsn
REAL, POINTER, DIMENSION(:,:,:)     :: PhiAngIn
REAL, POINTER, DIMENSION(:,:,:,:,:) :: DcmpPhiAngIn, DcmpPhiAngOut

INTEGER, POINTER, DIMENSION(:) :: AsyRayList, DirList, AziList

INTEGER :: mp(2)
INTEGER :: iazi, ipol, irotray, iasyray, iray, irayseg, idir, ipin, icel, iasy, ireg, isurf, ig, jbeg, jend, jinc, imray, ipray
INTEGER :: ipst, iped, ipinc, isfst, isfed, isgst, isged, isginc, iAsyTyp, iGeoTyp, icBss, jhPin, jcBss
INTEGER :: nAsyRay, nPinRay, nRaySeg, FsrIdxSt, nAziAng, nPolarAng, PhiAnginSvIdx, ExpAppIdx

REAL :: wtazi(10), wt2cs(10), wt2sn(10)

REAL, DIMENSION(RayInfo%nPolarAngle, gb:ge) :: PhiAngOut
REAL :: phid, tau, ExpApp, wtsurf

DATA mp /2, 1/
! ----------------------------------------------------

! Ray Info.
nPolarAng = RayInfo%nPolarAngle
nAziAng   = RayInfo%nAziAngle

! Geo.
Pin => CoreInfo%Pin

! Tracking Dat Pointing
src           => TrackingDat%srcnm
xst           => TrackingDat%xstnm
PhiAngIn      => TrackingDat%phiAngInnm
DcmpPhiAngIn  => TrackingDat%DcmpPhiAngIn
DcmpPhiAngOut => TrackingDat%DcmpPhiAngOut
EXPA          => TrackingDat%EXPA
EXPB          => TrackingDat%EXPB
wtang         => TrackingDat%wtang
wthcs         => TrackingDat%wthcs
wthsn         => TrackingDat%wthsn

! Dcmp. Ray
iRotRay     = DcmpAsyRay%iRotRay
iAsy        = DcmpAsyRay%iAsy
iRay        = DcmpAsyRay%iRay
AsyRayList => DcmpAsyRay%AsyRayList
DirList    => DcmpAsyRay%DirList
AziList    => DcmpAsyRay%AziList

! Ray. B.C.
PhiAngInSvIdx = RayInfo%PhiAngInSvIdx(iRotRay, krot)

IF (DcmpAsyRay%lRotRayBeg(krot)) THEN
  PhiAngOut = PhiAngIn(1:nPolarAng, gb:ge, PhiAnginSvIdx)
ELSE
  PhiAngOut = DcmpPhiAngIn(1:nPolarAng, gb:ge, krot, iRay, iAsy)
END IF
  
IF (krot .EQ. 2) THEN
  jbeg = nAsyRay; jend = 1; jinc = -1
ELSE
  jbeg = 1; jend = nAsyRay; jinc = 1
END IF
! ----------------------------------------------------
DO imray = jbeg, jend, jinc
  iAsyRay = AsyRayList(imray)
  iazi    = AziList(imray)
  idir    = DirList(imray)
  IF (krot .EQ. 2) idir = mp(idir)
    
  DO ipol = 1, nPolarAng
    wtazi(ipol) = wtang(ipol, iazi)
    wt2cs(ipol) = wthcs(ipol, iazi)
    wt2sn(ipol) = wthsn(ipol, iazi)
  END DO
  
  iAsyTyp = hAsy(iAsy)%AsyTyp
  iGeoTyp = hAsy(iAsy)%GeoTyp
  icBss   = hAsyTypInfo(iAsyTyp)%iBss
  
  haRay_Loc => haRay(iGeoTyp, icBss, iAsyRay)
  
  nPinRay = haRay_Loc%nhpRay
  
  IF (idir .EQ. 1) THEN
    ipst = 1; iped = nPinRay; ipinc = 1;  isfst = 1; isfed = 2;
  ELSE
    ipst = nPinRay; iped = 1; ipinc = -1; isfst = 2; isfed = 1;
  END IF
  ! --------------------------------------------------
  DO ipray = ipst, iped, ipinc
    jhPin = haRay_Loc%CelRay(ipRay)%hPinIdx
    jhPin = hAsy(iAsy)%PinIdxSt + hAsyTypInfo(iAsyTyp)%PinLocIdx(iGeoTyp, jhPin) - 1
    jcBss = Pin(jhPin)%hCelGeo(iz)
    
    CelRay_Loc => haRay(iGeoTyp, jcBss, iAsyRay)%CelRay(ipRay)
    
    ! Surface : In-coming
    IF (lJout) THEN
      iSurf = CelRay_Loc%hSufIdx(isfst)
      
      DO ipol = 1, nPolarAng
        wtsurf = abs(wt2cs(ipol) * CelRay_Loc%hsn(isfst) - wt2sn(ipol) * CelRay_Loc%hcs(isfst))
        
        DO ig = gb, ge
          Jout(1, ig, isurf, ipin) = Jout(1, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wtazi(ipol)
          Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wtsurf
        END DO
      END DO
    END IF
    
    ! Iter. : FSR
    nRaySeg = CelRay_Loc%nSegRay
    
    IF (idir .EQ. 1) THEN
      isgst = 1; isged = nRaySeg; isginc = 1
    ELSE
      isgst = nRaySeg; isged = 1; isginc = -1
    END IF
    
    DO iRaySeg = isgst, isged, isginc
      ireg = CelRay_Loc%MshIdx(iRaySeg) + Pin(jhPin)%FsrIdxSt - 1
      
      DO ig = gb, ge
        tau = -CelRay_Loc%SegLgh(iRaySeg) * xst(ig, ireg) ! Optimum Length
        
        ExpAppIdx = max(INT(tau), -40000)
        ExpAppIdx = min(0, ExpAppIdx)
        
        DO ipol = 1, nPolarAng
          ExpApp = EXPA(ipol, ExpAppIdx) * tau + EXPB(ipol, ExpAppIdx)
          
          phid = (PhiAngOut(ipol, ig) - src(ig, ireg)) * ExpApp
          
          PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
          
          phis(ig, ireg) = phis(ig, ireg) + wtazi(ipol) * phid
        END DO
      END DO
    END DO
    
    ! Surface : Out-going
    IF (lJout) THEN
      isurf = CelRay_Loc%hSufIdx(isfed) ! y : Big
      
      DO ipol = 1, nPolarAng
        wtsurf = abs(wt2cs(ipol) * CelRay_Loc%hsn(isfed) - wt2sn(ipol) * CelRay_Loc%hcs(isfed))
        
        DO ig = gb, ge
          Jout(2, ig, isurf, ipin) = Jout(2, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wtazi(ipol)
          Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wtsurf
        END DO
      END DO
    END IF
  END DO
END DO

DcmpPhiAngOut(1:nPolarAng, gb:ge, krot, iRay, iAsy) = PhiAngOut
! ----------------------------------------------------
! Geo.
NULLIFY (Pin)

! Hex.
NULLIFY (haRay_Loc)
NULLIFY (CelRay_Loc)

! Loc.
NULLIFY (src)
NULLIFY (xst)
NULLIFY (EXPA)
NULLIFY (EXPB)
NULLIFY (wtang)
NULLIFY (wthcs)
NULLIFY (wthsn)
NULLIFY (PhiAngIn)

! Dcmp.
NULLIFY (DcmpPhiAngIn)
NULLIFY (DcmpPhiAngOut)
NULLIFY (AsyRayList)
NULLIFY (DirList)
NULLIFY (AziList)
! ----------------------------------------------------

END SUBROUTINE HexTrackRotRayNM_Dcmp
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmp_LSCASMO(RayInfo, CoreInfo, phisnm, phisSlope, joutnm, iz, iasy, gb, ge, ljout)

USE PARAM
USE TYPEDEF,    ONLY :  RayInfo_Type,       Coreinfo_type,      Pin_Type,           Asy_Type,               &
                        AsyInfo_Type,       PinInfo_Type,       Cell_Type,          DcmpAsyRayInfo_Type
USE Moc_Mod,    ONLY :  nMaxRaySeg,         nMaxCellRay,        nMaxAsyRay,         nMaxCoreRay,            &
                        Expa,               Expb,               ApproxExp,          TrackingDat,            &
                        DcmpPhiAngIn,       DcmpPhiAngOut,      TrackRotRayLSDcmp_CASMO
USE geom,       ONLY :  ng
USE PE_Mod,     ONLY :  PE
USE ALLOCS
USE OMP_LIB
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phisnm(:, :)
REAL, POINTER :: joutnm(:, :, :, :)
REAL, POINTER :: phisSlope(:, :, :, :)
INTEGER :: iz, iasy, gb, ge
LOGICAL :: ljout

TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(DcmpAsyRayInfo_Type), POINTER :: DcmpAsyRay(:, :)
INTEGER, POINTER :: DcmpAsyRayCount(:)
REAL, POINTER :: FsrMxx(:), FsrMyy(:), FsrMxy(:)

INTEGER :: AsyType

INTEGER :: nxy, nFsr
INTEGER :: iAsyRay
INTEGER :: tid
INTEGER :: FsrIdxSt, FsrIdxEnd
INTEGER :: PinIdxSt, PinIdxEnd
INTEGER :: icel, ireg
INTEGER :: i, j, l, g
REAL :: detinv

REAL, POINTER :: phimx(:, :, :), phimy(:, :, :)
REAL, POINTER :: srcnm(:, :), xstnm(:, :)

AsyInfo => CoreInfo%AsyInfo; Asy => CoreInfo%Asy
Cell => CoreInfo%CellInfo; Pin => CoreInfo%Pin
FsrMxx => CoreInfo%CoreFsrMxx(:, iz)
FsrMyy => CoreInfo%CoreFsrMyy(:, iz)
FsrMxy => CoreInfo%CoreFsrMxy(:, iz)

AsyType = Asy(iAsy)%AsyType
nxy = AsyInfo(AsyType)%nxy

PinIdxSt = Asy(iAsy)%GlobalPinIdx(1)
PinIdxEnd = Asy(iAsy)%GlobalPinIdx(nxy)
FsrIdxSt = Pin(PinIdxSt)%FsrIdxSt
FsrIdxEnd = Pin(PinIdxEnd)%FsrIdxSt + Cell(Pin(PinIdxEnd)%Cell(iz))%nFsr - 1

DcmpAsyRay => RayInfo%DcmpAsyRay
DcmpAsyRayCount => RayInfo%DcmpAsyRayCount

ALLOCATE(phimx(2, gb : ge, FsrIdxSt : FsrIdxEnd), phimy(2, gb : ge, FsrIdxSt : FsrIdxEnd))
phimx(:, :, :) = zero; phimy(:, :, :) = zero

phisnm(gb : ge, FsrIdxSt : FsrIdxEnd) = zero
IF(ljout) joutnm(:, gb : ge, :, PinIdxSt : PinIdxEnd) = zero

tid = omp_get_thread_num() + 1

xstnm => TrackingDat(tid)%xstnm
srcnm => TrackingDat(tid)%srcnm

DO iAsyRay = 1, DcmpAsyRayCount(iAsy)
  CALL TrackRotRayLSDcmp_CASMO(RayInfo, CoreInfo, TrackingDat(tid), phisnm, phimx, phimy, joutnm,   &
                               DcmpAsyRay(iAsyRay, iAsy), ljout, iz, gb, ge)
END DO
    
DO i = FsrIdxSt, FsrIdxEnd
  DO g = gb, ge
    phimx(1, g, i) = phimx(1, g, i) / xstnm(g, i)
    phimy(1, g, i) = phimy(1, g, i) / xstnm(g, i)
  END DO
END DO

DO l = PinIdxSt, PinIdxEnd
  icel = Pin(l)%Cell(iz);
  DO j = 1, Cell(icel)%nFsr
    ireg = Pin(l)%FsrIdxSt + j - 1
    DO g = gb, ge
      phisnm(g, ireg) = phisnm(g, ireg) / (xstnm(g, ireg) * Cell(icel)%vol(j)) + srcnm(g, ireg)
      phimx(:, g, ireg) = phimx(:, g, ireg) / (xstnm(g, ireg) * Cell(icel)%vol(j))
      phimy(:, g, ireg) = phimy(:, g, ireg) / (xstnm(g, ireg) * Cell(icel)%vol(j))
    END DO
  END DO  
END DO

DO i = FsrIdxSt, FsrIdxEnd
  detinv = 1.0 / (FsrMxx(i) * FsrMyy(i) - FsrMxy(i) * FsrMxy(i))
  DO g = gb, ge
    phisSlope(1, g, i, iz) = detinv * (FsrMyy(i) * SUM(phimx(:, g, i)) - FsrMxy(i) * SUM(phimy(:, g, i)))
    phisSlope(2, g, i, iz) = detinv * (FsrMxx(i) * SUM(phimy(:, g, i)) - FsrMxy(i) * SUM(phimx(:, g, i)))
  END DO
END DO

DEALLOCATE(phimx, phimy)

NULLIFY(Cell); NULLIFY(Pin)
NULLIFY(FsrMxx); NULLIFY(FsrMyy); NULLIFY(FsrMxy)
NULLIFY(xstnm); NULLIFY(srcnm);

END SUBROUTINE RayTraceDcmp_LSCASMO
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE TrackRotRayLSDcmp_CASMO(RayInfo, CoreInfo, TrackingDat, phisC, phimx, phimy, jout, DcmpAsyRay, ljout, iz, gb, ge)
USE PARAM
USE TYPEDEF,    ONLY :  RayInfo_Type,       Coreinfo_type,      Pin_Type,           Asy_Type,               &
                        AsyInfo_Type,       PinInfo_Type,       Cell_Type,          AziAngleInfo_Type,      &
                        PolarAngle_Type,    AsyRayInfo_type,    CoreRayInfo_Type,   DcmpAsyRayInfo_Type,    &
                        RotRayInfo_Type,    CellRayInfo_type,   TrackingDat_Type
USE Moc_Mod,    ONLY :  nMaxCellRay,        nMaxCoreRay

IMPLICIT NONE

TYPE(RayInfo_Type), INTENT(INOUT) :: RayInfo
TYPE(CoreInfo_Type), INTENT(INOUT) :: CoreInfo
TYPE(TrackingDat_Type), INTENT(INOUT) :: TrackingDat
TYPE(DcmpAsyRayInfo_Type) :: DcmpAsyRay
REAL, POINTER :: phisC(:, :), phimx(:, :, :), phimy(:, :, :), jout(:, :, :, :)
LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: iz, gb, ge

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(PinInfo_Type), POINTER :: PinInfo(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(AziAngleInfo_Type), POINTER :: AziAng(:)
TYPE(PolarAngle_Type), POINTER :: PolarAng(:)
TYPE(AsyRayInfo_type), POINTER :: AsyRay(:)
TYPE(CoreRayInfo_Type), POINTER :: CoreRay(:)
TYPE(RotRayInfo_Type), POINTER :: RotRay(:)
TYPE(CellRayInfo_Type), POINTER :: CellRay

REAL, POINTER :: LenSeg(:)
INTEGER, POINTER :: LocalFsrIdx(:)
INTEGER, POINTER :: FsrIdx(:, :),  ExpAppIdx(:, :, :)
REAL, POINTER :: OptLenList(:, :, :)
REAL, POINTER :: srcC(:, :), srcSlope(:, :, :), xst(:, :)
REAL, POINTER :: PhiAngOut(:, :, :), PhiAngIn(:, :, :)
REAL, POINTER :: DcmpPhiAngOut(:, :, :, :, :), DcmpPhiAngIn(:, :, :, :, :)
REAL, POINTER :: EXPA(:, :), EXPB(:, :)
REAL, POINTER :: E1(:, :, :, :), E3(:, :, :, :), R1(:, :, :, :), R3(:, :, :, :)
REAL, POINTER :: cmOptLen(:, :, :, :), cmOptLenInv(:, :, :, :)
REAL, POINTER :: q0(:, :, :), q1(:, :, :, :)
REAL, POINTER :: x0(:, :, :), y0(:, :, :)
REAL, POINTER :: FsrCentroid(:, :)
INTEGER, POINTER :: AsyRayList(:), DirList(:), AziList(:)

INTEGER :: mp(2)
INTEGER :: iazi, ipol, ig, iRotRay, iasyray, iceray, irayseg, irot, itype, idir, iray
INTEGER :: nRotRay, nCoreRay, nAsyRay, nPinRay, nRaySeg, FsrIdxSt
INTEGER :: nPolarAng, nAziAng, nPhiAngSv
INTEGER :: ipin, icel, iasy, ireg, isurf
INTEGER :: irsegidx, icellrayidx
INTEGER :: nFsr, nxy
INTEGER :: i, j, k, l, m, jbeg, jend, jinc, ir, ir1
INTEGER :: ibcel

INTEGER :: CellRayIdxSt(DcmpAsyRay%nMaxCellRay, DcmpAsyRay%nAsyRay, 2)
INTEGER :: PinIdx(DcmpAsyRay%nMaxCellRay, DcmpAsyRay%nAsyRay)
INTEGER :: SurfIdx(DcmpAsyRay%nMaxCellRay, DcmpAsyRay%nAsyRay, 2)
INTEGER :: PhiAnginSvIdx(2)
INTEGER :: nTotRaySeg(DcmpAsyRay%nAsyRay), nTotCellRay(DcmpAsyRay%nAsyRay)

REAL :: wtang(10, 100), wttemp, wtang2(10, 100, 4)
REAL :: tau, phiobd(10, gb : ge), phid, phim, wt(10), wt2(10, 4)
REAL :: ax(10), ay(10), polsininv(10)
REAL :: segCenter(2), segSlope(2, 10)

DATA mp /2, 1/

!Ray Info Pointing
AziAng => RayInfo%AziAngle;
PolarAng => RayInfo%PolarAngle;
AsyRay => RayInfo%AsyRay
CoreRay => RayInfo%CoreRay
RotRay => RayInfo%RotRay

!Geometry Info Pointing
Asy => CoreInfo%Asy
Pin => CoreInfo%Pin
PinInfo => CoreInfo%Pininfo
Cell => CoreInfo%CellInfo
FsrCentroid => CoreInfo%CoreFsrCentroid(:, :, iz)

!Tracking Dat Pointing
FsrIdx => TrackingDat%FsrIdx
ExpAppIdx => TrackingDat%ExpAppIdxnm
OptLenList => TrackingDat%OptLenListnm
cmOptLen => TrackingDat%cmOptLen; cmOptLenInv => TrackingDat%cmOptLenInv
srcC => TrackingDat%srcnm; srcSlope => TrackingDat%srcSlope
xst => TrackingDat%xstnm
PhiAngOut => TrackingDat%PhiAngOutnm
PhiAngIn => TrackingDat%phiAngInnm
DcmpPhiAngIn => TrackingDat%DcmpPhiAngIn
DcmpPhiAngOut => TrackingDat%DcmpPhiAngOut
EXPA => TrackingDat%EXPA
EXPB => TrackingDat%EXPB
E1 => TrackingDat%E1; E3 => TrackingDat%E3
R1 => TrackingDat%R1; R3 => TrackingDat%R3
q0 => TrackingDat%q0; q1 => TrackingDat%q1
x0 => TrackingDat%x0; y0 => TrackingDat%y0

nAziAng = RayInfo%nAziAngle; nPolarAng = RayInfo%nPolarAngle
nPhiAngSv = RayInfo%nPhiAngSv
nFsr = CoreInfo%nCoreFsr; nxy = CoreInfo%nxy

DO ipol = 1, nPolarAng
  wttemp = PolarAng(ipol)%weight * PolarAng(ipol)%sinv
  polsininv(ipol) = 1._8 / PolarAng(ipol)%sinv
  DO iazi = 1, nAziAng
    wtang(ipol, iazi) = wttemp * AziAng(iazi)%weight * AziAng(iazi)%del
  END DO
END DO
IF(lJout)THEN
  DO ipol = 1, nPolarAng
    DO iazi = 1, nAziAng
      wtang2(ipol, iazi, 1) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv)
      wtang2(ipol, iazi, 3) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv)
      wtang2(ipol, iazi, 2) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv)
      wtang2(ipol, iazi, 4) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv)
    END DO
  END DO
END IF

nAsyRay = DcmpAsyRay%nAsyRay
iRotRay = DcmpAsyRay%iRotRay
AsyRayList => DcmpAsyRay%AsyRayList
DirList => DcmpAsyRay%DirList
AziList => DcmpAsyRay%AziList
iAsy = DcmpAsyRay%iAsy; iRay = DcmpAsyRay%iRay
PhiAngInSvIdx = RayInfo%PhiAngInSvIdx(iRotRay, 1 : 2)

DO j = 1, nAsyRay
  irsegidx = 0
  iAsyRay = AsyRayList(j); iAzi = AziList(j)
  nPinRay = AsyRay(iAsyRay)%nCellRay
  DO ipol = 1, nPolarAng
    segSlope(:, ipol) = (/ AziAng(iazi)%cosv * PolarAng(ipol)%sinv, AziAng(iazi)%sinv * PolarAng(ipol)%sinv /)
  END DO
  DO l = 1, nPinRay
    ipin = AsyRay(iAsyRay)%PinIdx(l)
    iceray = AsyRay(iAsyRay)%PinRayIdx(l)
    ipin = Asy(iAsy)%GlobalPinIdx(ipin)
    icel = Pin(ipin)%Cell(iz)
    FsrIdxSt = Pin(ipin)%FsrIdxSt
    ibcel = Cell(icel)%basecellstr
    CellRay => Cell(ibcel)%CellRay(iceray)
    PinIdx(l, j) = ipin
    CellRayIdxSt(l, j, 2) = irsegidx + 1
    nRaySeg = CellRay%nSeg
    LocalFsrIdx => CellRay%LocalFsrIdx
    LenSeg => CellRay%LenSeg
    DO iRaySeg = 1, nRaySeg
      ireg = FsrIdxSt + CellRay%LocalFsrIdx(iRaySeg) - 1
      irsegidx = irsegidx + 1
      FsrIdx(irsegidx, j) = ireg
      segCenter(:) = half * (CellRay%pts(:, iRaySeg) + CellRay%pts(:, iRaySeg + 1)) - FsrCentroid(:, ireg)
      x0(1, irsegidx, j) = segCenter(1) - LenSeg(iRaySeg) * 0.0005_8 * AziAng(iazi)%cosv
      x0(2, irsegidx, j) = segCenter(1) + LenSeg(iRaySeg) * 0.0005_8 * AziAng(iazi)%cosv
      y0(1, irsegidx, j) = segCenter(2) - LenSeg(iRaySeg) * 0.0005_8 * AziAng(iazi)%sinv
      y0(2, irsegidx, j) = segCenter(2) + LenSeg(iRaySeg) * 0.0005_8 * AziAng(iazi)%sinv
      DO ig = gb, ge
        tau = - LenSeg(iRaySeg) * xst(ig, ireg)
        OptLenList(ig, irsegidx, j) = tau
        q0(ig, irsegidx, j) = srcC(ig, ireg) + sum(srcSlope(1 : 2, ig, ireg) * segCenter(:))
        DO ipol = 1, nPolarAng
          cmOptLen(ipol, ig, irsegidx, j) = - tau * polsininv(ipol) * 0.001_8
          cmOptLenInv(ipol, ig, irsegidx, j) = 1._8 / cmOptLen(ipol, ig, irsegidx, j)
          q1(ipol, ig, irsegidx, j) = half * sum(srcSlope(3 : 4, ig, ireg) * segSlope(:, ipol))
        END DO
        ExpAppIdx(ig, irsegidx, j) = max(INT(tau), -40000)
        ExpAppIdx(ig, irsegidx, j) = min(0, ExpAppIdx(ig, irsegidx, j))
      END DO
    END DO
    CellRayIdxSt(l, j, 1) = irsegidx
    SurfIdx(l, j, 1) = AsyRay(iAsyRay)%PinRaySurf(2, l) !OutSurface
    SurfIdx(l, j, 2) = AsyRay(iAsyRay)%PinRaySurf(1, l) !Insurface
  END DO
  nTotRaySeg(j) = irsegidx
  nTotCellRay(j) = nPinRay
END DO

DO ipol = 1, nPolarAng
  DO j = 1, nAsyRay
    DO l = 1, nTotRaySeg(j)
      DO ig = gb, ge
        E1(ipol, ig, l, j) = EXPA(ExpAppIdx(ig, l, j), ipol) * OptLenList(ig, l, j) + EXPB(ExpAppIdx(ig, l, j), ipol)
        E3(ipol, ig, l, j) = 2._8 * (cmOptLen(ipol, ig, l, j) - E1(ipol, ig, l, j)) - cmOptLen(ipol, ig, l, j) * E1(ipol, ig, l, j)
        R1(ipol, ig, l, j) = 1._8 + half * cmOptLen(ipol, ig, l, j) - (1._8 + cmOptLenInv(ipol, ig, l, j)) * E1(ipol, ig, l, j)
        R3(ipol, ig, l, j) = rsix * cmOptLen(ipol, ig, l, j) - 2._8 * cmOptLenInv(ipol, ig, l, j) - 2._8  &
                             + (1._8 + cmOptLenInv(ipol, ig, l, j)) * (1._8 + 2._8 * cmOptLenInv(ipol, ig, l, j)) * E1(ipol, ig, l, j)
      END DO
    END DO
  END DO
END DO

DO irot = 1, 2
  IF (DcmpAsyRay%lRotRayBeg(irot)) THEN
    phiobd(1 : nPolarAng, gb : ge) = PhiAngIn(1 : nPolarAng, gb : ge, PhiAnginSvIdx(irot))
  ELSE
    phiobd(1 : nPolarAng, gb : ge) = DcmpPhiAngIn(1 : nPolarAng, gb : ge, irot, iRay, iAsy)
  END IF
  jinc = 1; jbeg = 1; jend = nAsyRay
  IF(irot .eq. 2) THEN
    jinc = -1; jbeg = nAsyRay; jend = 1
  END IF
  DO j = jbeg, jend, jinc
    idir = DirList(j); iazi = AziList(j)
    nRaySeg = nTotRaySeg(j)
    IF(irot .eq. 2) idir = mp(idir)
    DO ipol = 1, nPolarAng
      wt(ipol) = wtang(ipol, iazi)
      ax(ipol) = wt(ipol) * AziAng(iazi)%cosv * PolarAng(ipol)%sinv
      ay(ipol) = wt(ipol) * AziAng(iazi)%sinv * PolarAng(ipol)%sinv
    END DO
    IF(lJout) THEN
      wt2(1 : nPolarAng, :) = wtang2(1 : nPolarAng, iazi, :)
    END IF
    nRaySeg = nTotRaySeg(j)
    IF(idir .eq. 1) THEN
      PhiAngOut(1 : nPolarAng, gb : ge, 1) = phiobd(1 : nPolarAng, gb : ge)
      DO ir = 1, nRaySeg
        ireg = FsrIdx(ir, j)
        DO ig = gb, ge
          DO ipol = 1, nPolarAng          
            phid = (q0(ig, ir, j) - PhiAngOut(ipol, ig, ir)) * E1(ipol, ig, ir, j) + q1(ipol, ig, ir, j) * E3(ipol, ig, ir, j)
            PhiAngOut(ipol, ig, ir + 1) = PhiAngOut(ipol, ig, ir) + phid
            phisC(ig, ireg) = phisC(ig, ireg) - wt(ipol) * phid
            phim = PhiAngOut(ipol, ig, ir) * (cmOptLen(ipol, ig, ir, j) * half) + (q0(ig, ir, j) - PhiAngOut(ipol, ig, ir)) * R1(ipol, ig, ir, j)  &
                   + q1(ipol, ig, ir, j) * cmOptLen(ipol, ig, ir, j) * R3(ipol, ig, ir, j)
            phimx(1, ig, ireg) = phimx(1, ig, ireg) + ax(ipol) * phim * cmOptLen(ipol, ig, ir, j)
            phimy(1, ig, ireg) = phimy(1, ig, ireg) + ay(ipol) * phim * cmOptLen(ipol, ig, ir, j)
            phimx(2, ig, ireg) = phimx(2, ig, ireg) + wt(ipol) * x0(idir, ir, j) * (-phid + q0(ig, ir, j) * cmOptLen(ipol, ig, ir, j))
            phimy(2, ig, ireg) = phimy(2, ig, ireg) + wt(ipol) * y0(idir, ir, j) * (-phid + q0(ig, ir, j) * cmOptLen(ipol, ig, ir, j))
          END DO
        END DO
      END DO
      phiobd(1 : nPolarAng, gb : ge) = PhiAngOut(1 : nPolarAng, gb : ge, nRaySeg + 1)
      IF(ljout) THEN
        DO ir = 1, nTotCellRay(j)
          icel = PinIdx(ir, j); isurf = SurfIdx(ir, j, 1)
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(2, ig, isurf, icel) = Jout(2, ig, isurf, icel) + wt(ipol) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 1) + 1)
              Jout(3, ig, isurf, icel) = Jout(3, ig, isurf, icel) + wt2(ipol, isurf) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 1) + 1)
            END DO
          END DO
          isurf = SurfIdx(ir, j, 2)
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(1, ig, isurf, icel) = Jout(1, ig, isurf, icel) + wt(ipol) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 2))
              Jout(3, ig, isurf, icel) = Jout(3, ig, isurf, icel) + wt2(ipol, isurf) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 2))
            END DO
          END DO
        END DO
      END IF
    ELSE
      PhiAngOut(1 : nPolarAng, gb : ge, nRaySeg + 2) = phiobd(1 : nPolarAng, gb : ge)
      ir = nRaySeg + 1
      DO ir1 = 1, nRaySeg
        ir = ir - 1
        ireg = FsrIdx(ir, j)
        DO ig = gb, ge
          DO ipol = 1, nPolarAng
            phid = (q0(ig, ir, j) - PhiAngOut(ipol, ig, ir + 2)) * E1(ipol, ig, ir, j) - q1(ipol, ig, ir, j) * E3(ipol, ig, ir, j)
            PhiAngOut(ipol, ig, ir + 1) = PhiAngOut(ipol, ig, ir + 2) + phid
            phisC(ig, ireg) = phisC(ig, ireg) - wt(ipol) * phid
            phim = PhiAngOut(ipol, ig, ir + 2) * (cmOptLen(ipol, ig, ir, j) * half) + (q0(ig, ir, j) - PhiAngOut(ipol, ig, ir + 2)) * R1(ipol, ig, ir, j)  &
                   - q1(ipol, ig, ir, j) * cmOptLen(ipol, ig, ir, j) * R3(ipol, ig, ir, j)
            phimx(1, ig, ireg) = phimx(1, ig, ireg) - ax(ipol) * phim * cmOptLen(ipol, ig, ir, j)
            phimy(1, ig, ireg) = phimy(1, ig, ireg) - ay(ipol) * phim * cmOptLen(ipol, ig, ir, j)
            phimx(2, ig, ireg) = phimx(2, ig, ireg) + wt(ipol) * x0(idir, ir, j) * (-phid + q0(ig, ir, j) * cmOptLen(ipol, ig, ir, j))
            phimy(2, ig, ireg) = phimy(2, ig, ireg) + wt(ipol) * y0(idir, ir, j) * (-phid + q0(ig, ir, j) * cmOptLen(ipol, ig, ir, j))
          END DO
        END DO
      END DO
      phiobd(1 : nPolarAng, gb : ge) = PhiAngOut(1 : nPolarAng, gb : ge, 2)
      IF(lJout) THEN
        DO ir = 1, nTotCellRay(j)
          icel = PinIdx(ir, j); isurf = SurfIdx(ir, j, 2)
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(2, ig, isurf, icel) = Jout(2, ig, isurf, icel) + wt(ipol) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 2) + 1)
              Jout(3, ig, isurf, icel) = Jout(3, ig, isurf, icel) + wt2(ipol, isurf) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 2) + 1)
            END DO
          END DO
          isurf = SurfIdx(ir, j, 1)
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(1, ig, isurf, icel) = Jout(1, ig, isurf, icel) + wt(ipol) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 1) + 2)
              Jout(3, ig, isurf, icel) = Jout(3, ig, isurf, icel) + wt2(ipol, isurf) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 1) + 2)
            END DO
          END DO
        END DO
      END IF
    END IF
  END DO
  DcmpPhiAngOut(1 : nPolarAng, gb : ge, irot, iRay, iAsy) = phiobd(1 : nPolarAng, gb : ge)
END DO

NULLIFY(AziAng); NULLIFY(PolarAng)
NULLIFY(AsyRay); NULLIFY(CoreRay)
NULLIFY(RotRay); NULLIFY(CellRay)

!Geometry Info Pointing
NULLIFY(Asy); NULLIFY(Pin)
NULLIFY(PinInfo); NULLIFY(Cell)
NULLIFY(FsrCentroid)

!Tracking Dat Pointing
NULLIFY(FsrIdx); NULLIFY(ExpAppIdx)
NULLIFY(OptLenList); NULLIFY(cmOptLen); NULLIFY(cmOptLenInv)
NULLIFY(LenSeg); NULLIFY(LocalFsrIdx)
NULLIFY(srcC); NULLIFY(srcSlope)
NULLIFY(xst)
NULLIFY(PhiAngOut); NULLIFY(PhiAngIn)
NULLIFY(EXPA); NULLIFY(EXPB)
NULLIFY(E1); NULLIFY(E3); NULLIFY(R1); NULLIFY(R3)
NULLIFY(q0); NULLIFY(q1); NULLIFY(x0); NULLIFY(y0)

END SUBROUTINE TrackRotRayLSDcmp_CASMO
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmp_Pn(RayInfo, CoreInfo, phisnm, phimnm, joutnm, iz, iAsy, gb, ge, ScatOd, lJout)

USE ALLOCS
USE OMP_LIB
USE PARAM,   ONLY : TRUE, ZERO, RTHREE, RFIVE, RSEVEN
USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type, Pin_Type, Cell_Type, Asy_Type, AsyInfo_Type, DcmpAsyRayInfo_Type
USE Moc_Mod, ONLY : TrackingDat, Comp, mwt, AziMap, DcmpPhiAngIn, DcmpPhiAngOut, TrackRotRayPn_Dcmp, wtang, srcmnm
USE PE_MOD,  ONLY : PE

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:,:)     :: phisnm
REAL, POINTER, DIMENSION(:,:,:)   :: phimnm
REAL, POINTER, DIMENSION(:,:,:,:) :: joutnm

INTEGER :: iz, iAsy, gb, ge, ScatOd
LOGICAL :: ljout
! ----------------------------------------------------
TYPE (Asy_Type),            POINTER, DIMENSION(:)   :: Asy
TYPE (AsyInfo_Type),        POINTER, DIMENSION(:)   :: AsyInfo
TYPE (Cell_Type),           POINTER, DIMENSION(:)   :: Cell
TYPE (Pin_Type),            POINTER, DIMENSION(:)   :: Pin
TYPE (DcmpAsyRayInfo_Type), POINTER, DIMENSION(:,:) :: DcmpAsyRay

REAL, POINTER, DIMENSION(:,:)     :: xstnm, srcnm
REAL, POINTER, DIMENSION(:,:,:,:) :: phianm, SrcAngnm

INTEGER, POINTER, DIMENSION(:,:,:) :: DcmpAsyAziList

INTEGER :: nAziAng, nPolarAng, nAsyRay, nxy, ithr, FsrIdxSt, FsrIdxEnd, PinIdxSt, PinIdxEnd, AziIdx
INTEGER :: icel, ipin, iazi, ipol, iDcmpAsyRay, imRay, ig, iFSR, jFSR
REAL :: wttemp, tempsrc
! ----------------------------------------------------

Asy     => CoreInfo%Asy
AsyInfo => CoreInfo%AsyInfo
Cell    => CoreInfo%CellInfo
Pin     => CoreInfo%Pin

DcmpAsyRay     => RayInfo%DcmpAsyRay
DcmpAsyAziList => RayInfo%DcmpAsyAziList
nAziAng         = RayInfo%nAziAngle
nPolarAng       = RayInfo%nPolarAngle

nxy = AsyInfo(Asy(iAsy)%AsyType)%nxy

PinIdxSt  = Asy(iAsy)%GlobalPinIdx(1)
PinIdxEnd = Asy(iAsy)%GlobalPinIdx(nxy)

FsrIdxSt  = Pin(Asy(iAsy)%GlobalPinIdx(1))%FsrIdxSt
FsrIdxEnd = Pin(Asy(iAsy)%GlobalPinIdx(nxy))%FsrIdxSt + Cell(Pin(Asy(iAsy)%GlobalPinIdx(nxy))%Cell(iz))%nFsr - 1
! ----------------------------------------------------
ithr = omp_get_thread_num() + 1

CALL dmalloc0(TrackingDat(ithr)%phianm,   1, nPolarAng, gb, ge, FsrIdxSt, FsrIdxEnd, 1, 4)
CALL dmalloc0(TrackingDat(ithr)%SrcAngnm, 1, nPolarAng, gb, ge, FsrIdxSt, FsrIdxEnd, 1, 4)

xstnm    => TrackingDat(ithr)%xstnm
srcnm    => TrackingDat(ithr)%srcnm
phianm   => TrackingDat(ithr)%phianm
SrcAngnm => TrackingDat(ithr)%SrcAngnm

phisnm(   gb:ge, FsrIdxSt:FsrIdxEnd) = ZERO
phimnm(:, gb:ge, FsrIdxSt:FsrIdxEnd) = ZERO

IF (lJout) joutnm(:, gb:ge, :, PinIdxSt:PinIdxEnd) = ZERO
! ----------------------------------------------------
DO iAzi = 1, nAziAng / 2
  ! SET : Src Ang
  IF (ScatOd .EQ. 1) THEN
    DO iFSR = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = srcnm(ig, iFSR)
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = srcnm(ig, iFSR)
          
          tempsrc = Comp(1, ipol, AziIdx) * srcmnm(1, ig, iFSR) + Comp(2, ipol, AziIdx) * srcmnm(2, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) - tempsrc
        END DO
      END DO
      
      AziIdx = nAziAng - iAzi + 1
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = srcnm(ig, iFSR)
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = srcnm(ig, iFSR)
          
          tempsrc = Comp(1, ipol, AziIdx) * srcmnm(1, ig, iFSR) + Comp(2, ipol, AziIdx) * srcmnm(2, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) - tempsrc
        END DO
      END DO
    END DO
  ELSE IF (ScatOd .EQ. 2) THEN 
    DO iFSR = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = srcnm(ig, iFSR)
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = srcnm(ig, iFSR)
          
          tempsrc = Comp(1, ipol, AziIdx) * srcmnm(1, ig, iFSR) + Comp(2, ipol, AziIdx) * srcmnm(2, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) - tempsrc
          
          tempsrc = Comp(3, ipol, AziIdx) * srcmnm(3, ig, iFSR) + Comp(4, ipol, AziIdx) * srcmnm(4, ig, iFSR) + Comp(5, ipol, AziIdx) * srcmnm(5, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) + tempsrc
        END DO
      END DO
      
      AziIdx = nAziAng - iAzi + 1
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = srcnm(ig, iFSR)
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = srcnm(ig, iFSR)
          
          tempsrc = Comp(1, ipol, AziIdx) * srcmnm(1, ig, iFSR) + Comp(2, ipol, AziIdx) * srcmnm(2, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) - tempsrc
          
          tempsrc = Comp(3, ipol, AziIdx) * srcmnm(3, ig, iFSR) + Comp(4, ipol, AziIdx) * srcmnm(4, ig, iFSR) + Comp(5, ipol, AziIdx) * srcmnm(5, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) + tempsrc
        END DO
      END DO
    END DO
  ELSEIF (ScatOd .EQ. 3) THEN
    DO iFSR = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = srcnm(ig, iFSR)
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = srcnm(ig, iFSR)
          
          tempsrc = Comp(1, ipol, AziIdx) * srcmnm(1, ig, iFSR) + Comp(2, ipol, AziIdx) * srcmnm(2, ig, iFSR) &
                  + Comp(6, ipol, AziIdx) * srcmnm(6, ig, iFSR) + Comp(7, ipol, AziIdx) * srcmnm(7, ig, iFSR) &
                  + Comp(8, ipol, AziIdx) * srcmnm(8, ig, iFSR) + Comp(9, ipol, AziIdx) * srcmnm(9, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) - tempsrc
          
          tempsrc = Comp(3, ipol, AziIdx) * srcmnm(3, ig, iFSR) + Comp(4, ipol, AziIdx) * srcmnm(4, ig, iFSR) + Comp(5, ipol, AziIdx) * srcmnm(5, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) + tempsrc
        END DO
      END DO
      
      AziIdx = nAziAng - iAzi + 1
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = srcnm(ig, iFSR)
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = srcnm(ig, iFSR)
          
          tempsrc = Comp(1, ipol, AziIdx) * srcmnm(1, ig, iFSR) + Comp(2, ipol, AziIdx) * srcmnm(2, ig, iFSR) &
                  + Comp(6, ipol, AziIdx) * srcmnm(6, ig, iFSR) + Comp(7, ipol, AziIdx) * srcmnm(7, ig, iFSR) &
                  + Comp(8, ipol, AziIdx) * srcmnm(8, ig, iFSR) + Comp(9, ipol, AziIdx) * srcmnm(9, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) - tempsrc
          
          tempsrc = Comp(3, ipol, AziIdx) * srcmnm(3, ig, iFSR) + Comp(4, ipol, AziIdx) * srcmnm(4, ig, iFSR) + Comp(5, ipol, AziIdx) * srcmnm(5, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) + tempsrc
        END DO
      END DO
    END DO
  END IF
  
  ! TRACK
  phianm = ZERO
  
  DO imRay = 1, DcmpAsyAziList(0, iAzi, iAsy)
    iDcmpAsyRay = DcmpAsyAziList(imRay, iAzi, iAsy)
    
    CALL TrackRotRayPn_Dcmp(RayInfo, CoreInfo, TrackingDat(ithr), DcmpAsyRay(iDcmpAsyRay, iAsy), Joutnm, lJout, iz, gb, ge)
  END DO
  
  ! CnP
  IF (ScatOd .EQ. 1) THEN
    DO iFSR = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          phisnm(ig, iFSR) = phisnm(ig, iFSR) + wtang(ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          
          phimnm(1:2, ig, iFSR) = phimnm(1:2, ig, iFSR) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
        END DO
      END DO
      
      AziIdx = nAziAng - iAzi + 1
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          phisnm(ig, iFSR) = phisnm(ig, iFSR) + wtang(ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          
          phimnm(1:2, ig, iFSR) = phimnm(1:2, ig, iFSR) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
        END DO
      END DO
    END DO
  ELSEIF (ScatOd .EQ. 2) THEN
    DO iFSR = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          phisnm(ig, iFSR) = phisnm(ig, iFSR) + wtang(ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          
          phimnm(1:2, ig, iFSR) = phimnm(1:2, ig, iFSR) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          phimnm(3:5, ig, iFSR) = phimnm(3:5, ig, iFSR) + mwt(3:5, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
        END DO
      END DO
      
      AziIdx = nAziAng - iAzi + 1
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          phisnm(ig, iFSR) = phisnm(ig, iFSR) + wtang(ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          
          phimnm(1:2, ig, iFSR) = phimnm(1:2, ig, iFSR) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          phimnm(3:5, ig, iFSR) = phimnm(3:5, ig, iFSR) + mwt(3:5, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
        END DO
      END DO
    END DO
  ELSEIF (ScatOd .EQ. 3) THEN
    DO iFSR = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          phisnm(ig, iFSR) = phisnm(ig, iFSR) + wtang(ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          
          phimnm(1:2, ig, iFSR) = phimnm(1:2, ig, iFSR) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          phimnm(3:5, ig, iFSR) = phimnm(3:5, ig, iFSR) + mwt(3:5, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          phimnm(6:9, ig, iFSR) = phimnm(6:9, ig, iFSR) + mwt(6:9, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
        END DO
      END DO
      
      AziIdx = nAziAng - iAzi + 1
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          phisnm(ig, iFSR) = phisnm(ig, iFSR) + wtang(ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          
          phimnm(1:2, ig, iFSR) = phimnm(1:2, ig, iFSR) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          phimnm(3:5, ig, iFSR) = phimnm(3:5, ig, iFSR) + mwt(3:5, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          phimnm(6:9, ig, iFSR) = phimnm(6:9, ig, iFSR) + mwt(6:9, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
        END DO
      END DO
    END DO
  END IF
END DO

DEALLOCATE (TrackingDat(ithr)%phianm)
DEALLOCATE (TrackingDat(ithr)%SrcAngnm)
! ----------------------------------------------------
IF (ScatOd .EQ. 1) THEN
  !$OMP PARALLEL DO PRIVATE(iPin, FsrIdxSt, iFSR, jFSR, icel, wttemp) SCHEDULE(DYNAMIC)
  DO iPin = PinIdxSt, PinIdxEnd
    FsrIdxSt = Pin(iPin)%FsrIdxSt
    icel     = Pin(iPin)%Cell(iz)
    
    DO iFSR = 1, Cell(icel)%nFsr
      jFSR = FsrIdxSt + iFSR - 1
      
      DO ig = gb, ge
        wttemp = 1._8 / (xstnm(ig, jFSR) * Cell(icel)%vol(iFSR))
        
        phisnm(ig, jFSR) = phisnm(ig, jFSR) * wttemp + srcnm(ig, jFSR)
        
        phimnm(1:2, ig, jFSR) = phimnm(1:2, ig, jFSR) * wttemp + srcmnm(1:2, ig, jFSR) * RTHREE
      END DO
    END DO  
  END DO
  !$OMP END PARALLEL DO
ELSEIF (ScatOd .EQ. 2) THEN
  !$OMP PARALLEL DO PRIVATE(iPin, FsrIdxSt, iFSR, jFSR, icel, wttemp) SCHEDULE(DYNAMIC)
  DO iPin = PinIdxSt, PinIdxEnd
    FsrIdxSt = Pin(iPin)%FsrIdxSt
    icel     = Pin(iPin)%Cell(iz)
    
    DO iFSR = 1, Cell(icel)%nFsr
      jFSR = FsrIdxSt + iFSR - 1
      
      DO ig = gb, ge
        wttemp = 1._8 / (xstnm(ig, jFSR) * Cell(icel)%vol(iFSR))
        
        phisnm(ig, jFSR) = phisnm(ig, jFSR) * wttemp + srcnm(ig, jFSR)
        
        phimnm(1:2, ig, jFSR) = phimnm(1:2, ig, jFSR) * wttemp + srcmnm(1:2, ig, jFSR) * rthree
        phimnm(3:5, ig, jFSR) = phimnm(3:5, ig, jFSR) * wttemp + srcmnm(3:5, ig, jFSR) * rfive
      END DO
    END DO  
  END DO
  !$OMP END PARALLEL DO
ELSEIF (ScatOd .EQ. 3) THEN
  !$OMP PARALLEL DO PRIVATE(iPin, FsrIdxSt, iFSR, jFSR, icel, wttemp) SCHEDULE(DYNAMIC)
  DO iPin = PinIdxSt, PinIdxEnd
    FsrIdxSt = Pin(iPin)%FsrIdxSt
    icel     = Pin(iPin)%Cell(iz)
    
    DO iFSR = 1, Cell(icel)%nFsr
      jFSR = FsrIdxSt + iFSR - 1
      
      DO ig = gb, ge
        wttemp = 1._8 / (xstnm(ig, jFSR) * Cell(icel)%vol(iFSR))
        
        phisnm(ig, jFSR) = phisnm(ig, jFSR) * wttemp + srcnm(ig, jFSR)
        
        phimnm(1:2, ig, jFSR) = phimnm(1:2, ig, jFSR) * wttemp + srcmnm(1:2, ig, jFSR) * rthree
        phimnm(3:5, ig, jFSR) = phimnm(3:5, ig, jFSR) * wttemp + srcmnm(3:5, ig, jFSR) * rfive
        phimnm(6:9, ig, jFSR) = phimnm(6:9, ig, jFSR) * wttemp + srcmnm(6:9, ig, jFSR) * rseven
      END DO
    END DO  
  END DO
  !$OMP END PARALLEL DO
END IF
! ----------------------------------------------------
NULLIFY (Asy)
NULLIFY (AsyInfo)
NULLIFY (Cell)
NULLIFY (Pin)
NULLIFY (DcmpAsyRay)
NULLIFY (DcmpAsyAziList)
NULLIFY (xstnm)
NULLIFY (srcnm)
NULLIFY (phianm)
NULLIFY (SrcAngnm)
! ----------------------------------------------------

END SUBROUTINE RayTraceDcmp_Pn
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE TrackRotRayPn_Dcmp(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, Jout, lJout, iz, gb, ge)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, Asy_Type, Cell_Type, AsyRayInfo_type, CoreRayInfo_Type, CellRayInfo_type, TrackingDat_Type, DcmpAsyRayInfo_Type

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

REAL, POINTER, DIMENSION(:,:,:,:) :: Jout

LOGICAL :: lJout
INTEGER :: iz, gb, ge, ScatOd

TYPE (Pin_Type),          POINTER, DIMENSION(:) :: Pin
TYPE (Asy_Type),          POINTER, DIMENSION(:) :: Asy
TYPE (Cell_Type),         POINTER, DIMENSION(:) :: Cell
TYPE (AsyRayInfo_type),   POINTER, DIMENSION(:) :: AsyRay
TYPE (CoreRayInfo_Type),  POINTER, DIMENSION(:) :: CoreRay

TYPE (CellRayInfo_Type), POINTER :: CellRay

REAL, POINTER, DIMENSION(:)         :: LenSeg
REAL, POINTER, DIMENSION(:,:)       :: xst, EXPA, EXPB, wtang
REAL, POINTER, DIMENSION(:,:,:)     :: PhiAngIn, wtsurf
REAL, POINTER, DIMENSION(:,:,:,:)   :: SrcAng, phia
REAL, POINTER, DIMENSION(:,:,:,:,:) :: DcmpPhiAngIn, DcmpPhiAngOut

INTEGER, POINTER, DIMENSION(:)   :: LocalFsrIdx, AsyRayList, DirList, AziList
INTEGER, POINTER, DIMENSION(:,:) :: AziMap

INTEGER :: mp(2), AziSvIdx(2)
INTEGER :: iAzi, ipol, iRotRay, iAsyRay, jAsyRay, iRay, iRaySeg, irot, idir, jbeg, jend, jinc
INTEGER :: ipin, icel, iasy, ireg, isurf, ibcel, iceray, ig, ir, iPinRay
INTEGER :: nCoreRay, nAsyRay, nDcmpRay, nPinRay, nRaySeg, FsrIdxSt, nPolarAng, nAziAng, PhiAnginSvIdx, ExpAppIdx

REAL :: wt(10), wt2(10, 4)
REAL :: PhiAngOut(RayInfo%nPolarAngle, gb : ge)
REAL :: phid, tau, ExpApp

DATA mp /2, 1/
! ----------------------------------------------------

! Ray Info Pointing
nAziAng   = RayInfo%nAziAngle
nPolarAng = RayInfo%nPolarAngle
AsyRay   => RayInfo%AsyRay
CoreRay  => RayInfo%CoreRay

! Geometry Info Pointing
Asy  => CoreInfo%Asy
Pin  => CoreInfo%Pin
Cell => CoreInfo%CellInfo

! Tracking Dat Pointing
xst           => TrackingDat%xstnm
PhiAngIn      => TrackingDat%PhiAngInnm
DcmpPhiAngIn  => TrackingDat%DcmpPhiAngIn
DcmpPhiAngOut => TrackingDat%DcmpPhiAngOut
EXPA          => TrackingDat%EXPA
EXPB          => TrackingDat%EXPB
wtang         => TrackingDat%wtang
wtsurf        => TrackingDat%wtsurf
AziMap        => TrackingDat%AziMap
phia          => TrackingDat%phianm
SrcAng        => TrackingDat%SrcAngnm

nAsyRay     = DcmpAsyRay%nAsyRay
iRotRay     = DcmpAsyRay%iRotRay
iAsy        = DcmpAsyRay%iAsy
iRay        = DcmpAsyRay%iRay
AsyRayList => DcmpAsyRay%AsyRayList
DirList    => DcmpAsyRay%DirList
AziList    => DcmpAsyRay%AziList
! ----------------------------------------------------
DO irot = 1, 2  
  PhiAngInSvIdx = RayInfo%PhiAngInSvIdx(iRotRay, irot)
  
  IF (DcmpAsyRay%lRotRayBeg(irot)) THEN
    PhiAngOut = PhiAngIn(1:nPolarAng, gb:ge, PhiAnginSvIdx)
  ELSE
    PhiAngOut = DcmpPhiAngIn(1:nPolarAng, gb:ge, irot, iRay, iAsy)
  END IF
  
  jbeg = 1; jend = nAsyRay; jinc = 1
  IF (irot .EQ. 2) THEN
    jbeg = nAsyRay; jend = 1; jinc = -1
  END IF
  
  DO iAsyRay = jbeg, jend, jinc
    jAsyRay = AsyRayList(iAsyRay)
    nPinRay = AsyRay(jAsyRay)%nCellRay
    idir    = DirList(iAsyRay)
    iazi    = AziList(iAsyRay)
    IF (irot .eq. 2) idir = mp(idir)
    
    AziSvIdx = AziMap(iAzi, :)
    
    DO ipol = 1, nPolarAng
      wt(ipol) = wtang(ipol, iazi)
      
      IF (lJout) wt2(ipol, 1:4) = wtsurf(ipol, iazi, 1:4)
    END DO
    
    IF (idir .EQ. 1) THEN
      DO iPinRay = 1, nPinRay
        ipin   = AsyRay(jAsyRay)%PinIdx   (iPinRay)
        iceray = AsyRay(jAsyRay)%PinRayIdx(iPinRay)
        
        ipin = Asy(iAsy)%GlobalPinIdx(ipin)   
        
        icel     = Pin(ipin)%Cell(iz)
        FsrIdxSt = Pin(ipin)%FsrIdxSt
        
        ibcel    = Cell(icel)%basecellstr
        CellRay => Cell(ibcel)%CellRay(iceray)
        nRaySeg  = CellRay%nSeg
        
        LocalFsrIdx => CellRay%LocalFsrIdx
        LenSeg      => CellRay%LenSeg
        
        IF (lJout) THEN
          isurf = AsyRay(jAsyRay)%PinRaySurf(1, iPinRay)
          
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(1, ig, isurf, ipin) = Jout(1, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
              Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
            END DO
          END DO
        END IF
        
        DO iRaySeg = 1, nRaySeg
          ireg = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
          
          DO ig = gb, ge
            tau = - LenSeg(iRaySeg) * xst(ig, ireg)
            
            ExpAppIdx = max(INT(tau), -40000)
            ExpAppIdx = min(0, ExpAppIdx)
            
            DO ipol = 1, nPolarAng
              ExpApp = EXPA(ipol, ExpAppIdx) * tau + EXPB(ipol, ExpAppIdx)
              
              phid = (PhiAngOut(ipol, ig) - SrcAng(ipol, ig, ireg, AziSvIdx(idir))) * ExpApp
              
              PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
              
              phia(ipol, ig, ireg, AziSvIdx(idir)) = phia(ipol, ig, ireg, AziSvIdx(idir)) + phid
            END DO
          END DO
        END DO
        
        IF (lJout) THEN
          isurf = AsyRay(jAsyRay)%PinRaySurf(2, iPinRay)
          
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(2, ig, isurf, ipin) = Jout(2, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
              Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
            END DO
          END DO
        END IF
      END DO
    ELSE
      DO iPinRay = nPinRay, 1, -1
        ipin   = AsyRay(jAsyRay)%PinIdx   (iPinRay)
        iceray = AsyRay(jAsyRay)%PinRayIdx(iPinRay)
        
        ipin = Asy(iAsy)%GlobalPinIdx(ipin)   
        
        icel     = Pin(ipin)%Cell(iz)             
        FsrIdxSt = Pin(ipin)%FsrIdxSt
        
        ibcel    = Cell(icel)%basecellstr
        CellRay => Cell(ibcel)%CellRay(iceray)
        nRaySeg  = CellRay%nSeg
        
        LocalFsrIdx => CellRay%LocalFsrIdx
        LenSeg      => CellRay%LenSeg
        
        IF (lJout) THEN
          isurf = AsyRay(jAsyRay)%PinRaySurf(2, iPinRay)
          
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(1, ig, isurf, ipin) = Jout(1, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
              Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
            END DO
          END DO
        END IF
        
        DO iRaySeg = nRaySeg, 1, -1
          ireg = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
          
          DO ig = gb, ge
            tau = - LenSeg(iRaySeg) * xst(ig, ireg)
            
            ExpAppIdx = max(INT(tau), -40000)
            ExpAppIdx = min(0, ExpAppIdx)
            
            DO ipol = 1, nPolarAng
              ExpApp = EXPA(ipol, ExpAppIdx) * tau + EXPB(ipol, ExpAppIdx)
              
              phid = (PhiAngOut(ipol, ig) - SrcAng(ipol, ig, ireg, AziSvIdx(idir))) * ExpApp
              
              PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
              
              phia(ipol, ig, ireg, AziSvIdx(idir)) = phia(ipol, ig, ireg, AziSvIdx(idir)) + phid
            END DO
          END DO
        END DO
        
        IF (lJout) THEN
          isurf = AsyRay(jAsyRay)%PinRaySurf(1, iPinRay)
          
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(2, ig, isurf, ipin) = Jout(2, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
              Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
            END DO
          END DO
        END IF
      END DO
    END IF
  END DO
  
  DcmpPhiAngOut(1:nPolarAng, gb:ge, irot, iRay, iAsy) = PhiAngOut
END DO
! ----------------------------------------------------
! Geom
NULLIFY (Pin)
NULLIFY (Asy)
NULLIFY (Cell)

! Ray
NULLIFY (AsyRay)
NULLIFY (CellRay)

! Loc.
NULLIFY (LenSeg)
NULLIFY (srcAng)
NULLIFY (phia)
NULLIFY (xst)
NULLIFY (EXPA)
NULLIFY (EXPB)
NULLIFY (wtang)
NULLIFY (wtsurf)
NULLIFY (PhiAngIn)
NULLIFY (LocalFsrIdx)

! Dcmp.
NULLIFY (DcmpPhiAngIn)
NULLIFY (DcmpPhiAngOut)
NULLIFY (AsyRayList)
NULLIFY (DirList)
NULLIFY (AziList)
NULLIFY (AziMap)
! ----------------------------------------------------

END SUBROUTINE TrackRotRayPn_Dcmp
! ------------------------------------------------------------------------------------------------------------