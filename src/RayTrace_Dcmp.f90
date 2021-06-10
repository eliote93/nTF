#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTrace_Dcmp(RayInfo, CoreInfo, iz, gb, ge, lJout, lSubGrp, lScat1, lLinSrcCASMO, lHybrid)

USE OMP_LIB
USE PARAM,       ONLY : ZERO, RED, BLACK, GREEN
USE TYPEDEF,     ONLY : RayInfo_Type, Coreinfo_type, Asy_Type, AsyInfo_Type
USE Moc_Mod,     ONLY : TrackingDat, phisNM, phimNM, srcNM, xstNM, MocjoutNM, PhiAngInNM, DcmpPhiAngIn, DcmpPhiAngOut, &
                        RayTraceDcmp_OMP, RayTraceDcmp_Pn, RayTraceDcmp_LSCASMO, DcmpGatherBoundaryFlux, DcmpScatterBoundaryFlux, DcmpLinkBoundaryFlux
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

INTEGER :: color, startColor, endColor, colorInc, ithr, nThr, ScatOd, AsyType, iAsy
! ----------------------------------------------------

AsyInfo => CoreInfo%AsyInfo
Asy     => CoreInfo%Asy

ScatOd = nTracerCntl%ScatOd

IF (.NOT. nTracerCntl%lHex) THEN
  IF (mod(itrcntl%mocit, 2) .EQ. 0) THEN
    startColor = BLACK
    endColor   = RED
    colorInc   = RED - BLACK
  ELSE
    startColor = RED
    endColor   = BLACK
    colorInc   = BLACK - RED
  END IF
ELSE
  IF (mod(itrcntl%mocit, 2) .EQ. 0) THEN
    startColor = GREEN
    endColor   = RED
    colorInc   = (RED - GREEN)/2
  ELSE
    startColor = RED
    endColor   = GREEN
    colorInc   = (GREEN - RED)/2
  END IF
END IF

nthr = PE%nthread
CALL OMP_SET_NUM_THREADS(nThr)
! ----------------------------------------------------
DO ithr = 1, nThr
  TrackingDat(ithr)%srcNM => srcNM
  TrackingDat(ithr)%xstNM => xstNM
  
  IF (lLinSrcCASMO) TrackingDat(ithr)%srcSlope => srcSlope(:, :, :, iz)
END DO

DcmpPhiAngOut(:, gb:ge, :, :, :) = ZERO
! ----------------------------------------------------
DO color = startColor, endColor, colorInc
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) CALL DcmpScatterBoundaryFlux(RayInfo, PhiAngInNM, DcmpPhiAngIn)
#endif

  !$OMP PARALLEL PRIVATE(ithr, iAsy, AsyType)
  ithr = 1
  !$ ithr = omp_get_thread_num()+1
  
  TrackingDat(ithr)%PhiAngInNM    => PhiAngInNM
  TrackingDat(ithr)%DcmpPhiAngIn  => DcmpPhiAngIn
  TrackingDat(ithr)%DcmpPhiAngOut => DcmpPhiAngOut
  !$OMP BARRIER
  !$OMP DO SCHEDULE(GUIDED)
  DO iAsy = PE%myAsyBeg, PE%myAsyEnd
    IF (Asy(iAsy)%color .NE. color) CYCLE
    
    AsyType = Asy(iAsy)%AsyType
    
    IF (lScat1 .AND. .NOT.lSubGrp) THEN
      CALL RayTraceDcmp_Pn(RayInfo, CoreInfo, phisNM, phimNM, MocjoutNM, iz, iAsy, gb, ge, ScatOd, lJout)
    ELSE
      IF (.NOT.lLinSrcCASMO .OR. (AsyInfo(AsyType)%lFuel.AND.lHybrid)) THEN
        CALL RayTraceDcmp_OMP(RayInfo, CoreInfo, phisNM, MocjoutNM, iz, iAsy, gb, ge, ljout)
      ELSE
        CALL RayTraceDcmp_LSCASMO(RayInfo, CoreInfo, phisNM, phisSlope, MocjoutNM, iz, iAsy, gb, ge, ljout)
      END IF
    END IF
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
  
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) CALL DcmpGatherBoundaryFlux(RayInfo, DcmpPhiAngOut)
#endif
  
  IF (PE%RTMASTER) CALL DcmpLinkBoundaryFlux(CoreInfo, RayInfo, PhiAngInNM, DcmpPhiAngIn, DcmpPhiAngOut, gb, ge, color)
END DO

NULLIFY (Asy)
NULLIFY (AsyInfo)
! ----------------------------------------------------

END SUBROUTINE RayTrace_Dcmp
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmp_OMP(RayInfo, CoreInfo, phisNM, joutNM, iz, iasy, gb, ge, ljout)

USE OMP_LIB
USE PARAM,   ONLY : ZERO
USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, AsyInfo_Type, Asy_Type, Pin_Type, Cell_Type, DcmpAsyRayInfo_Type
USE Moc_Mod, ONLY : RecTrackRotRayOMP_Dcmp, HexTrackRotRayOMP_Dcmp, TrackingDat
USE CNTL,    ONLY : nTracerCntl
USE HexData, ONLY : hAsy

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:,:)     :: phisNM
REAL, POINTER, DIMENSION(:,:,:,:) :: joutNM

INTEGER :: iz, iAsy, gb, ge
LOGICAL :: ljout
! ----------------------------------------------------
TYPE (AsyInfo_Type),        POINTER, DIMENSION(:)   :: AsyInfo
TYPE (Asy_Type),            POINTER, DIMENSION(:)   :: Asy
TYPE (Cell_Type),           POINTER, DIMENSION(:)   :: Cell
TYPE (Pin_Type),            POINTER, DIMENSION(:)   :: Pin
TYPE (DcmpAsyRayInfo_Type), POINTER, DIMENSION(:,:) :: DcmpAsyRay

INTEGER, POINTER, DIMENSION(:) :: DcmpAsyRayCount

REAL, POINTER, DIMENSION(:,:) :: xstNM, srcNM

INTEGER :: nxy, iAsyRay, ithr, icel, ig, ipin, ifsr, jfsr, krot, FsrIdxSt, FsrIdxEnd, PinIdxSt, PinIdxEd
! ----------------------------------------------------

AsyInfo => CoreInfo%AsyInfo
Asy     => CoreInfo%Asy
Cell    => CoreInfo%CellInfo
Pin     => CoreInfo%Pin

DcmpAsyRay      => RayInfo%DcmpAsyRay
DcmpAsyRayCount => RayInfo%DcmpAsyRayCount
! ----------------------------------------------------
IF (.NOT. nTracerCntl%lHex) THEN
  nxy = AsyInfo(Asy(iAsy)%AsyType)%nxy
  
  PinIdxSt = Asy(iAsy)%GlobalPinIdx(1)
  PinIdxEd = Asy(iAsy)%GlobalPinIdx(nxy)
ELSE
  PinIdxSt = hAsy(iAsy)%PinIdxSt
  PinIdxEd = hAsy(iAsy)%PinIdxSt + hAsy(iAsy)%nTotPin - 1
END IF

FsrIdxSt  = Pin(PinIdxSt)%FsrIdxSt
FsrIdxEnd = Pin(PinIdxEd)%FsrIdxSt + Cell(Pin(PinIdxEd)%Cell(iz))%nFsr - 1

phisNM(gb:ge, FsrIdxSt:FsrIdxEnd) = ZERO

IF (ljout) joutNM(:, gb:ge, :, PinIdxSt:PinIdxEd) = ZERO
! ----------------------------------------------------
ithr = omp_get_thread_num() + 1

IF (nTracerCntl%lHex) THEN
  DO iAsyRay = 1, DcmpAsyRayCount(iAsy)
    DO krot = 1, 2
      CALL HexTrackRotRayOMP_Dcmp(RayInfo, CoreInfo, TrackingDat(ithr), DcmpAsyRay(iAsyRay, iAsy), phisNM, joutNM, ljout, iz, gb, ge, krot)
    END DO
  END DO
ELSE
  DO iAsyRay = 1, DcmpAsyRayCount(iAsy)
    DO krot = 1, 2
      CALL RecTrackRotRayOMP_Dcmp(RayInfo, CoreInfo, TrackingDat(ithr), DcmpAsyRay(iAsyRay, iAsy), phisNM, joutNM, ljout, iz, gb, ge, krot)
    END DO
  END DO
END IF
! ----------------------------------------------------
xstNM => TrackingDat(ithr)%xstNM
srcNM => TrackingDat(ithr)%srcNM

DO ipin = PinIdxSt, PinIdxEd
  FsrIdxSt = Pin(ipin)%FsrIdxSt
  icel     = Pin(ipin)%Cell(iz)
  
  DO ifsr = 1, Cell(icel)%nFsr
    jfsr = FsrIdxSt + ifsr - 1
    
    DO ig = gb, ge
      phisNM(ig, jfsr) = phisNM(ig, jfsr) / (xstNM(ig, jfsr) * Cell(icel)%vol(ifsr)) + srcNM(ig, jfsr)
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
NULLIFY (xstNM)
NULLIFY (srcNM)
! ----------------------------------------------------

END SUBROUTINE RayTraceDcmp_OMP
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayOMP_Dcmp(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, phisNM, joutNM, ljout, iz, gb, ge, krot)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, Asy_Type, Cell_Type, AsyRayInfo_type, CellRayInfo_type, TrackingDat_Type, DcmpAsyRayInfo_Type

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

REAL, POINTER, DIMENSION(:,:)     :: phisNM
REAL, POINTER, DIMENSION(:,:,:,:) :: joutNM

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: iz, gb, ge, krot
! ----------------------------------------------------
TYPE (Pin_Type),        POINTER, DIMENSION(:) :: Pin
TYPE (Asy_Type),        POINTER, DIMENSION(:) :: Asy
TYPE (Cell_Type),       POINTER, DIMENSION(:) :: Cell
TYPE (AsyRayInfo_type), POINTER, DIMENSION(:) :: AsyRay

TYPE (CellRayInfo_Type),  POINTER :: CellRay

REAL, POINTER, DIMENSION(:)         :: LenSeg
REAL, POINTER, DIMENSION(:,:)       :: srcNM, xstNM, EXPA, EXPB, wtang
REAL, POINTER, DIMENSION(:,:,:)     :: PhiAngInNM, wtsurf
REAL, POINTER, DIMENSION(:,:,:,:,:) :: DcmpPhiAngIn, DcmpPhiAngOut

INTEGER, POINTER, DIMENSION(:) :: LocalFsrIdx, AsyRayList, DirList, AziList

INTEGER :: mp(2)
INTEGER :: iazi, ipol, irotray, iasyray, iceray, iray, irayseg, idir, ipin, icel, ibcel, iasy, iFSR, isurf, ig, jbeg, jend, jinc, imray, ipray
INTEGER :: ipst, iped, ipinc, isfst, isfed, isgst, isged, isginc
INTEGER :: nAsyRay, nPinRay, nRaySeg, FsrIdxSt, nPolarAng, PhiAnginSvIdx, ExpAppIdx

REAL :: wt(10), wt2(10, 4)

REAL, DIMENSION(RayInfo%nPolarAngle, gb:ge) :: PhiAngOut
REAL :: phid, tau, ExpApp

DATA mp /2, 1/
! ----------------------------------------------------

! Ray Info.
nPolarAng = RayInfo%nPolarAngle
AsyRay   => RayInfo%AsyRay

! Geo.
Asy  => CoreInfo%Asy
Pin  => CoreInfo%Pin
Cell => CoreInfo%CellInfo

! Tracking Dat Pointing
srcNM         => TrackingDat%srcNM
xstNM         => TrackingDat%xstNM
PhiAngInNM    => TrackingDat%PhiAngInNM
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
  PhiAngOut(1:nPolarAng, gb:ge) = PhiAngInNM  (1:nPolarAng, gb:ge, PhiAnginSvIdx)
ELSE
  PhiAngOut(1:nPolarAng, gb:ge) = DcmpPhiAngIn(1:nPolarAng, gb:ge, krot, iRay, iAsy)
END IF
  
IF (krot .EQ. 2) THEN
  jbeg = nAsyRay; jend = 1; jinc = -1
ELSE
  jbeg = 1; jend = nAsyRay; jinc = 1
END IF
! ----------------------------------------------------  
DO imray = jbeg, jend, jinc
  iAsyRay = AsyRayList(imray)
  iazi    = AziList   (imray)
  idir    = DirList   (imray)
  IF (krot .EQ. 2) idir = mp(idir)
  
  DO ipol = 1, nPolarAng
    wt(ipol) = wtang(ipol, iazi)
    
    IF (lJout) wt2(ipol, 1:4) = wtsurf(ipol, iazi, 1:4)
  END DO
  
  nPinRay = AsyRay(iAsyRay)%nCellRay
  
  IF (idir .EQ. 1) THEN
    ipst = 1; iped = nPinRay; ipinc = 1;  isfst = 1; isfed = 2
  ELSE
    ipst = nPinRay; iped = 1; ipinc = -1; isfst = 2; isfed = 1
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
          joutNM(1, ig, isurf, ipin) = joutNM(1, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
          joutNM(3, ig, isurf, ipin) = joutNM(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
        END DO
      END DO
    END IF
    
    IF (idir .EQ. 1) THEN
      isgst = 1; isged = nRaySeg; isginc = 1
    ELSE
      isgst = nRaySeg; isged = 1; isginc = -1
    END IF
    
    DO iRaySeg = isgst, isged, isginc
      iFSR = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
      
      DO ig = gb, ge
        tau = -LenSeg(iRaySeg) * xstNM(ig, iFSR)
        
        ExpAppIdx = max(INT(tau), -40000)
        ExpAppIdx = min(0, ExpAppIdx)
        
        DO ipol = 1, nPolarAng
          ExpApp = ExpA(ExpAppIdx, ipol) * tau + ExpB(ExpAppIdx, ipol)
          
          phid = (PhiAngOut(ipol, ig) - srcNM(ig, iFSR)) * ExpApp
          
          PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
          
          phisNM(ig, iFSR) = phisNM(ig, iFSR) + wt(ipol) * phid
        END DO
      END DO
    END DO
    
    IF (lJout) THEN
      isurf = AsyRay(iAsyRay)%PinRaySurf(isfed, ipray)
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          joutNM(2, ig, isurf, ipin) = joutNM(2, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
          joutNM(3, ig, isurf, ipin) = joutNM(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
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
NULLIFY (srcNM)
NULLIFY (xstNM)
NULLIFY (EXPA)
NULLIFY (EXPB)
NULLIFY (wtang)
NULLIFY (wtsurf)
NULLIFY (PhiAngInNM)
NULLIFY (LocalFsrIdx)

! Dcmp.
NULLIFY (DcmpPhiAngIn)
NULLIFY (DcmpPhiAngOut)
NULLIFY (AsyRayList)
NULLIFY (DirList)
NULLIFY (AziList)
! ----------------------------------------------------

END SUBROUTINE RecTrackRotRayOMP_Dcmp
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayOMP_Dcmp(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, phisNM, joutNM, ljout, iz, gb, ge, krot)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, TrackingDat_Type, DcmpAsyRayInfo_Type, AziAngleInfo_Type
USE HexData, ONLY : hAsy
USE HexType, ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexCoreRay, Type_HexRotRay
USE HexData, ONLY : haRay, hAsyTypInfo

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

REAL, POINTER, DIMENSION(:,:)     :: phisNM
REAL, POINTER, DIMENSION(:,:,:,:) :: joutNM

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: iz, gb, ge, krot
! ----------------------------------------------------
TYPE (Pin_Type),          POINTER, DIMENSION(:) :: Pin
TYPE (AziAngleInfo_Type), POINTER, DIMENSION(:) :: AziAng
  
TYPE (Type_HexAsyRay), POINTER :: haRay_Loc
TYPE (Type_HexCelRay), POINTER :: CelRay_Loc

REAL, POINTER, DIMENSION(:,:)       :: srcNM, xstNM, EXPA, EXPB, wtang, hwt
REAL, POINTER, DIMENSION(:,:,:)     :: PhiAngInNM
REAL, POINTER, DIMENSION(:,:,:,:,:) :: DcmpPhiAngIn, DcmpPhiAngOut

INTEGER, POINTER, DIMENSION(:) :: AsyRayList, DirList, AziList

INTEGER :: mp(2)
INTEGER :: iazi, ipol, irotray, iasyray, iray, irayseg, idir, icel, iasy, iFSR, isurf, ig, jbeg, jend, jinc, imray, ipray
INTEGER :: ipst, iped, ipinc, isfst, isfed, isgst, isged, isginc, iAsyTyp, iGeoTyp, icBss, jhPin, jcBss
INTEGER :: nAsyRay, nPinRay, nRaySeg, FsrIdxSt, nPolarAng, PhiAnginSvIdx, ExpAppIdx

REAL :: wtazi(10), locwt(10), loccs, locsn

REAL, DIMENSION(RayInfo%nPolarAngle, gb:ge) :: PhiAngOut
REAL :: phid, tau, ExpApp, wtsurf

DATA mp /2, 1/
! ----------------------------------------------------

! Ray Info.
nPolarAng = RayInfo%nPolarAngle
AziAng   => RayInfo%AziAngle

! Geo.
Pin => CoreInfo%Pin

! Tracking Dat
srcNM         => TrackingDat%srcNM
xstNM         => TrackingDat%xstNM
PhiAngInNM    => TrackingDat%PhiAngInNM
DcmpPhiAngIn  => TrackingDat%DcmpPhiAngIn
DcmpPhiAngOut => TrackingDat%DcmpPhiAngOut
EXPA          => TrackingDat%EXPA
EXPB          => TrackingDat%EXPB
wtang         => TrackingDat%wtang
hwt           => TrackingDat%hwt

! Dcmp.
nAsyRay     = DcmpAsyRay%nAsyRay
iRotRay     = DcmpAsyRay%iRotRay
iAsy        = DcmpAsyRay%iAsy
iRay        = DcmpAsyRay%iRay
AsyRayList => DcmpAsyRay%AsyRayList
DirList    => DcmpAsyRay%DirList
AziList    => DcmpAsyRay%AziList

iAsyTyp = hAsy(iAsy)%AsyTyp
iGeoTyp = hAsy(iAsy)%GeoTyp
icBss   = hAsyTypInfo(iAsyTyp)%iBss

! Ray. B.C.
PhiAngInSvIdx = RayInfo%PhiAngInSvIdx(iRotRay, krot)

IF (DcmpAsyRay%lRotRayBeg(krot)) THEN
  PhiAngOut(1:nPolarAng, gb:ge) = PhiAngInNM  (1:nPolarAng, gb:ge, PhiAnginSvIdx)
ELSE
  PhiAngOut(1:nPolarAng, gb:ge) = DcmpPhiAngIn(1:nPolarAng, gb:ge, krot, iRay, iAsy)
END IF

IF (krot .EQ. 1) THEN
  jbeg = 1; jend = nAsyRay; jinc = 1
ELSE
  jbeg = nAsyRay; jend = 1; jinc = -1
END IF
! ----------------------------------------------------
DO imray = jbeg, jend, jinc
  iAsyRay = AsyRayList(imray)
  iazi    = AziList   (imray)
  idir    = DirList   (imray)
  IF (krot .EQ. 2) idir = mp(idir)
  
  DO ipol = 1, nPolarAng
    wtazi(ipol) = wtang(ipol, iazi)
    locwt(ipol) = hwt  (ipol, iazi)
  END DO
  
  loccs = AziAng(iazi)%cosv
  locsn = AziAng(iazi)%sinv
  
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
        wtsurf = locwt(ipol) / abs(loccs * CelRay_Loc%hsn(isfst) - locsn * CelRay_Loc%hcs(isfst))
        
        DO ig = gb, ge
          joutNM(1, ig, isurf, jhpin) = joutNM(1, ig, isurf, jhpin) + PhiAngOut(ipol, ig) * wtazi(ipol)
          joutNM(3, ig, isurf, jhpin) = joutNM(3, ig, isurf, jhpin) + PhiAngOut(ipol, ig) * wtsurf
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
      iFSR = CelRay_Loc%MshIdx(iRaySeg) + Pin(jhPin)%FsrIdxSt - 1
      
      DO ig = gb, ge
        tau = -CelRay_Loc%SegLgh(iRaySeg) * xstNM(ig, iFSR) ! Optimum Length
        
        ExpAppIdx = max(INT(tau), -40000)
        ExpAppIdx = min(0, ExpAppIdx)
        
        DO ipol = 1, nPolarAng
          ExpApp = ExpA(ExpAppIdx, ipol) * tau + ExpB(ExpAppIdx, ipol)
          
          phid = (PhiAngOut(ipol, ig) - srcNM(ig, iFSR)) * ExpApp
          
          PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
          
          phisNM(ig, iFSR) = phisNM(ig, iFSR) + wtazi(ipol) * phid
        END DO
      END DO
    END DO
    
    ! Surface : Out-going
    IF (lJout) THEN
      isurf = CelRay_Loc%hSufIdx(isfed) ! y : Big
      
      DO ipol = 1, nPolarAng
        wtsurf = locwt(ipol) / abs(loccs * CelRay_Loc%hsn(isfed) - locsn * CelRay_Loc%hcs(isfed))
        
        DO ig = gb, ge
          joutNM(2, ig, isurf, jhpin) = joutNM(2, ig, isurf, jhpin) + PhiAngOut(ipol, ig) * wtazi(ipol)
          joutNM(3, ig, isurf, jhpin) = joutNM(3, ig, isurf, jhpin) + PhiAngOut(ipol, ig) * wtsurf
        END DO
      END DO
    END IF
  END DO
END DO

DcmpPhiAngOut(1:nPolarAng, gb:ge, krot, iRay, iAsy) = PhiAngOut(1:nPolarAng, gb:ge)
! ----------------------------------------------------
! Geo.
NULLIFY (Pin)

! Ray
NULLIFY (AziAng)
NULLIFY (haRay_Loc)
NULLIFY (CelRay_Loc)

! Loc.
NULLIFY (srcNM)
NULLIFY (xstNM)
NULLIFY (EXPA)
NULLIFY (EXPB)
NULLIFY (wtang)
NULLIFY (hwt)
NULLIFY (PhiAngInNM)

! Dcmp.
NULLIFY (DcmpPhiAngIn)
NULLIFY (DcmpPhiAngOut)
NULLIFY (AsyRayList)
NULLIFY (DirList)
NULLIFY (AziList)
! ----------------------------------------------------

END SUBROUTINE HexTrackRotRayOMP_Dcmp
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmp_LSCASMO(RayInfo, CoreInfo, phisNM, phisSlope, joutNM, iz, iasy, gb, ge, ljout)

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
REAL, POINTER :: phisNM(:, :)
REAL, POINTER :: joutNM(:, :, :, :)
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
INTEGER :: PinIdxSt, PinIdxEd
INTEGER :: icel, ireg
INTEGER :: i, j, l, g
REAL :: detinv

REAL, POINTER :: phimx(:, :, :), phimy(:, :, :)
REAL, POINTER :: srcNM(:, :), xstNM(:, :)

AsyInfo => CoreInfo%AsyInfo; Asy => CoreInfo%Asy
Cell => CoreInfo%CellInfo; Pin => CoreInfo%Pin
FsrMxx => CoreInfo%CoreFsrMxx(:, iz)
FsrMyy => CoreInfo%CoreFsrMyy(:, iz)
FsrMxy => CoreInfo%CoreFsrMxy(:, iz)

AsyType = Asy(iAsy)%AsyType
nxy = AsyInfo(AsyType)%nxy

PinIdxSt = Asy(iAsy)%GlobalPinIdx(1)
PinIdxEd = Asy(iAsy)%GlobalPinIdx(nxy)
FsrIdxSt = Pin(PinIdxSt)%FsrIdxSt
FsrIdxEnd = Pin(PinIdxEd)%FsrIdxSt + Cell(Pin(PinIdxEd)%Cell(iz))%nFsr - 1

DcmpAsyRay => RayInfo%DcmpAsyRay
DcmpAsyRayCount => RayInfo%DcmpAsyRayCount

ALLOCATE(phimx(2, gb : ge, FsrIdxSt : FsrIdxEnd), phimy(2, gb : ge, FsrIdxSt : FsrIdxEnd))
phimx(:, :, :) = zero; phimy(:, :, :) = zero

phisNM(gb : ge, FsrIdxSt : FsrIdxEnd) = zero
IF(ljout) joutNM(:, gb : ge, :, PinIdxSt : PinIdxEd) = zero

tid = omp_get_thread_num() + 1

xstNM => TrackingDat(tid)%xstNM
srcNM => TrackingDat(tid)%srcNM

DO iAsyRay = 1, DcmpAsyRayCount(iAsy)
  CALL TrackRotRayLSDcmp_CASMO(RayInfo, CoreInfo, TrackingDat(tid), phisNM, phimx, phimy, joutNM,   &
                               DcmpAsyRay(iAsyRay, iAsy), ljout, iz, gb, ge)
END DO
    
DO i = FsrIdxSt, FsrIdxEnd
  DO g = gb, ge
    phimx(1, g, i) = phimx(1, g, i) / xstNM(g, i)
    phimy(1, g, i) = phimy(1, g, i) / xstNM(g, i)
  END DO
END DO

DO l = PinIdxSt, PinIdxEd
  icel = Pin(l)%Cell(iz);
  DO j = 1, Cell(icel)%nFsr
    ireg = Pin(l)%FsrIdxSt + j - 1
    DO g = gb, ge
      phisNM(g, ireg) = phisNM(g, ireg) / (xstNM(g, ireg) * Cell(icel)%vol(j)) + srcNM(g, ireg)
      phimx(:, g, ireg) = phimx(:, g, ireg) / (xstNM(g, ireg) * Cell(icel)%vol(j))
      phimy(:, g, ireg) = phimy(:, g, ireg) / (xstNM(g, ireg) * Cell(icel)%vol(j))
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
NULLIFY(xstNM); NULLIFY(srcNM);

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
REAL, POINTER :: srcC(:, :), srcSlope(:, :, :), xstNM(:, :)
REAL, POINTER :: PhiAngOut(:, :, :), PhiAngInNM(:, :, :)
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
srcC => TrackingDat%srcNM; srcSlope => TrackingDat%srcSlope
xstNM => TrackingDat%xstNM
PhiAngOut => TrackingDat%PhiAngOutnm
PhiAngInNM => TrackingDat%PhiAngInNM
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
        tau = - LenSeg(iRaySeg) * xstNM(ig, ireg)
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
    phiobd(1 : nPolarAng, gb : ge) = PhiAngInNM(1 : nPolarAng, gb : ge, PhiAnginSvIdx(irot))
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
NULLIFY(xstNM)
NULLIFY(PhiAngOut); NULLIFY(PhiAngInNM)
NULLIFY(EXPA); NULLIFY(EXPB)
NULLIFY(E1); NULLIFY(E3); NULLIFY(R1); NULLIFY(R3)
NULLIFY(q0); NULLIFY(q1); NULLIFY(x0); NULLIFY(y0)

END SUBROUTINE TrackRotRayLSDcmp_CASMO
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmp_Pn(RayInfo, CoreInfo, phisNM, phimNM, joutNM, iz, iAsy, gb, ge, ScatOd, lJout)

USE ALLOCS
USE OMP_LIB
USE PARAM,   ONLY : TRUE, ZERO, RTHREE, RFIVE, RSEVEN
USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type, Pin_Type, Cell_Type, Asy_Type, AsyInfo_Type, DcmpAsyRayInfo_Type
USE Moc_Mod, ONLY : TrackingDat, Comp, mwt, AziMap, DcmpPhiAngIn, DcmpPhiAngOut, RecTrackRotRayPn_Dcmp, wtang, srcmNM
USE PE_MOD,  ONLY : PE
USE CNTL,    ONLY : nTracerCntl
USE HexData, ONLY : hAsy

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:,:)     :: phisNM
REAL, POINTER, DIMENSION(:,:,:)   :: phimNM
REAL, POINTER, DIMENSION(:,:,:,:) :: joutNM

INTEGER :: iz, iAsy, gb, ge, ScatOd
LOGICAL :: ljout
! ----------------------------------------------------
TYPE (Asy_Type),            POINTER, DIMENSION(:)   :: Asy
TYPE (AsyInfo_Type),        POINTER, DIMENSION(:)   :: AsyInfo
TYPE (Cell_Type),           POINTER, DIMENSION(:)   :: Cell
TYPE (Pin_Type),            POINTER, DIMENSION(:)   :: Pin
TYPE (DcmpAsyRayInfo_Type), POINTER, DIMENSION(:,:) :: DcmpAsyRay

REAL, POINTER, DIMENSION(:,:)     :: xstNM, srcNM
REAL, POINTER, DIMENSION(:,:,:,:) :: phianm, SrcAngnm

INTEGER, POINTER, DIMENSION(:) :: DcmpAsyRayCount

INTEGER :: nAziAng, nPolarAng, nAsyRay, nxy, ithr, FsrIdxSt, FsrIdxEnd, PinIdxSt, PinIdxEd, AziIdx
INTEGER :: icel, ipin, iazi, ipol, iDcmpAsyRay, imRay, ig, iFSR, jFSR, iAsyRay, krot
REAL :: wttemp, tempsrc
! ----------------------------------------------------

Asy     => CoreInfo%Asy
AsyInfo => CoreInfo%AsyInfo
Cell    => CoreInfo%CellInfo
Pin     => CoreInfo%Pin

nAziAng          = RayInfo%nAziAngle
nPolarAng        = RayInfo%nPolarAngle
DcmpAsyRay      => RayInfo%DcmpAsyRay
DcmpAsyRayCount => RayInfo%DcmpAsyRayCount

IF (.NOT. nTracerCntl%lHex) THEN
  nxy = AsyInfo(Asy(iAsy)%AsyType)%nxy
  
  PinIdxSt = Asy(iAsy)%GlobalPinIdx(1)
  PinIdxEd = Asy(iAsy)%GlobalPinIdx(nxy)
ELSE
  PinIdxSt = hAsy(iAsy)%PinIdxSt
  PinIdxEd = hAsy(iAsy)%PinIdxSt + hAsy(iAsy)%nTotPin - 1
END IF

FsrIdxSt  = Pin(PinIdxSt)%FsrIdxSt
FsrIdxEnd = Pin(PinIdxEd)%FsrIdxSt + Cell(Pin(PinIdxEd)%Cell(iz))%nFsr - 1
! ----------------------------------------------------
ithr = omp_get_thread_num() + 1

CALL dmalloc0(TrackingDat(ithr)%phianm,   1, nPolarAng, gb, ge, FsrIdxSt, FsrIdxEnd, 1, 4) ! Can be Huge Time-consuming
CALL dmalloc0(TrackingDat(ithr)%SrcAngnm, 1, nPolarAng, gb, ge, FsrIdxSt, FsrIdxEnd, 1, 4)

xstNM    => TrackingDat(ithr)%xstNM
srcNM    => TrackingDat(ithr)%srcNM
SrcAngnm => TrackingDat(ithr)%SrcAngnm

phisNM(   gb:ge, FsrIdxSt:FsrIdxEnd) = ZERO
phimNM(:, gb:ge, FsrIdxSt:FsrIdxEnd) = ZERO

IF (lJout) joutNM(:, gb:ge, :, PinIdxSt:PinIdxEd) = ZERO
! ----------------------------------------------------
DO iAzi = 1, nAziAng / 2
  IF (ScatOd .EQ. 1) THEN
    DO iFSR = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = srcNM(ig, iFSR)
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = srcNM(ig, iFSR)
          
          tempsrc = Comp(1, ipol, AziIdx) * srcmNM(1, ig, iFSR) + Comp(2, ipol, AziIdx) * srcmNM(2, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) - tempsrc
        END DO
      END DO
      
      AziIdx = nAziAng - iAzi + 1
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = srcNM(ig, iFSR)
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = srcNM(ig, iFSR)
          
          tempsrc = Comp(1, ipol, AziIdx) * srcmNM(1, ig, iFSR) + Comp(2, ipol, AziIdx) * srcmNM(2, ig, iFSR)
          
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
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = srcNM(ig, iFSR)
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = srcNM(ig, iFSR)
          
          tempsrc = Comp(1, ipol, AziIdx) * srcmNM(1, ig, iFSR) + Comp(2, ipol, AziIdx) * srcmNM(2, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) - tempsrc
          
          tempsrc = Comp(3, ipol, AziIdx) * srcmNM(3, ig, iFSR) + Comp(4, ipol, AziIdx) * srcmNM(4, ig, iFSR) + Comp(5, ipol, AziIdx) * srcmNM(5, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) + tempsrc
        END DO
      END DO
      
      AziIdx = nAziAng - iAzi + 1
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = srcNM(ig, iFSR)
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = srcNM(ig, iFSR)
          
          tempsrc = Comp(1, ipol, AziIdx) * srcmNM(1, ig, iFSR) + Comp(2, ipol, AziIdx) * srcmNM(2, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) - tempsrc
          
          tempsrc = Comp(3, ipol, AziIdx) * srcmNM(3, ig, iFSR) + Comp(4, ipol, AziIdx) * srcmNM(4, ig, iFSR) + Comp(5, ipol, AziIdx) * srcmNM(5, ig, iFSR)
          
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
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = srcNM(ig, iFSR)
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = srcNM(ig, iFSR)
          
          tempsrc = Comp(1, ipol, AziIdx) * srcmNM(1, ig, iFSR) + Comp(2, ipol, AziIdx) * srcmNM(2, ig, iFSR) &
                  + Comp(6, ipol, AziIdx) * srcmNM(6, ig, iFSR) + Comp(7, ipol, AziIdx) * srcmNM(7, ig, iFSR) &
                  + Comp(8, ipol, AziIdx) * srcmNM(8, ig, iFSR) + Comp(9, ipol, AziIdx) * srcmNM(9, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) - tempsrc
          
          tempsrc = Comp(3, ipol, AziIdx) * srcmNM(3, ig, iFSR) + Comp(4, ipol, AziIdx) * srcmNM(4, ig, iFSR) + Comp(5, ipol, AziIdx) * srcmNM(5, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) + tempsrc
        END DO
      END DO
      
      AziIdx = nAziAng - iAzi + 1
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = srcNM(ig, iFSR)
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = srcNM(ig, iFSR)
          
          tempsrc = Comp(1, ipol, AziIdx) * srcmNM(1, ig, iFSR) + Comp(2, ipol, AziIdx) * srcmNM(2, ig, iFSR) &
                  + Comp(6, ipol, AziIdx) * srcmNM(6, ig, iFSR) + Comp(7, ipol, AziIdx) * srcmNM(7, ig, iFSR) &
                  + Comp(8, ipol, AziIdx) * srcmNM(8, ig, iFSR) + Comp(9, ipol, AziIdx) * srcmNM(9, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) - tempsrc
          
          tempsrc = Comp(3, ipol, AziIdx) * srcmNM(3, ig, iFSR) + Comp(4, ipol, AziIdx) * srcmNM(4, ig, iFSR) + Comp(5, ipol, AziIdx) * srcmNM(5, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) + tempsrc
        END DO
      END DO
    END DO
  END IF
END DO
! ----------------------------------------------------
IF (nTracerCntl%lHex) THEN
  DO iAsyRay = 1, DcmpAsyRayCount(iAsy)
    DO krot = 1, 2
      !CALL HexTrackRotRayPn_Dcmp(RayInfo, CoreInfo, TrackingDat(ithr), DcmpAsyRay(iAsyRay, iAsy), joutNM, ljout, iz, gb, ge, krot)
    END DO
  END DO
ELSE
  DO iAsyRay = 1, DcmpAsyRayCount(iAsy)
    DO krot = 1, 2
      CALL RecTrackRotRayPn_Dcmp(RayInfo, CoreInfo, TrackingDat(ithr), DcmpAsyRay(iAsyRay, iAsy), joutNM, ljout, iz, gb, ge, krot)
    END DO
  END DO
END IF

!DO iAzi = 1, nAziAng / 2
!  DO imRay = 1, DcmpAsyAziList(0, iAzi, iAsy)
!    iDcmpAsyRay = DcmpAsyAziList(imRay, iAzi, iAsy)
!    
!    CALL RecTrackRotRayPn_Dcmp(RayInfo, CoreInfo, TrackingDat(ithr), DcmpAsyRay(iDcmpAsyRay, iAsy), joutNM, lJout, iz, gb, ge)
!  END DO
!END DO
! ----------------------------------------------------
phianm => TrackingDat(ithr)%phianm

DO iAzi = 1, nAziAng / 2
  IF (ScatOd .EQ. 1) THEN
    DO iFSR = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          phisNM(ig, iFSR) = phisNM(ig, iFSR) + wtang(ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          
          phimNM(1:2, ig, iFSR) = phimNM(1:2, ig, iFSR) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
        END DO
      END DO
      
      AziIdx = nAziAng - iAzi + 1
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          phisNM(ig, iFSR) = phisNM(ig, iFSR) + wtang(ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          
          phimNM(1:2, ig, iFSR) = phimNM(1:2, ig, iFSR) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
        END DO
      END DO
    END DO
  ELSEIF (ScatOd .EQ. 2) THEN
    DO iFSR = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          phisNM(ig, iFSR) = phisNM(ig, iFSR) + wtang(ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          
          phimNM(1:2, ig, iFSR) = phimNM(1:2, ig, iFSR) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          phimNM(3:5, ig, iFSR) = phimNM(3:5, ig, iFSR) + mwt(3:5, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
        END DO
      END DO
      
      AziIdx = nAziAng - iAzi + 1
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          phisNM(ig, iFSR) = phisNM(ig, iFSR) + wtang(ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          
          phimNM(1:2, ig, iFSR) = phimNM(1:2, ig, iFSR) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          phimNM(3:5, ig, iFSR) = phimNM(3:5, ig, iFSR) + mwt(3:5, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
        END DO
      END DO
    END DO
  ELSEIF (ScatOd .EQ. 3) THEN
    DO iFSR = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          phisNM(ig, iFSR) = phisNM(ig, iFSR) + wtang(ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          
          phimNM(1:2, ig, iFSR) = phimNM(1:2, ig, iFSR) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          phimNM(3:5, ig, iFSR) = phimNM(3:5, ig, iFSR) + mwt(3:5, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          phimNM(6:9, ig, iFSR) = phimNM(6:9, ig, iFSR) + mwt(6:9, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
        END DO
      END DO
      
      AziIdx = nAziAng - iAzi + 1
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          phisNM(ig, iFSR) = phisNM(ig, iFSR) + wtang(ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          
          phimNM(1:2, ig, iFSR) = phimNM(1:2, ig, iFSR) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          phimNM(3:5, ig, iFSR) = phimNM(3:5, ig, iFSR) + mwt(3:5, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          phimNM(6:9, ig, iFSR) = phimNM(6:9, ig, iFSR) + mwt(6:9, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
        END DO
      END DO
    END DO
  END IF
END DO

DEALLOCATE (TrackingDat(ithr)%phianm)
DEALLOCATE (TrackingDat(ithr)%SrcAngnm)
! ----------------------------------------------------
IF (ScatOd .EQ. 1) THEN
  DO iPin = PinIdxSt, PinIdxEd
    FsrIdxSt = Pin(iPin)%FsrIdxSt
    icel     = Pin(iPin)%Cell(iz)
    
    DO iFSR = 1, Cell(icel)%nFsr
      jFSR = FsrIdxSt + iFSR - 1
      
      DO ig = gb, ge
        wttemp = 1._8 / (xstNM(ig, jFSR) * Cell(icel)%vol(iFSR))
        
        phisNM(ig, jFSR) = phisNM(ig, jFSR) * wttemp + srcNM(ig, jFSR)
        
        phimNM(1:2, ig, jFSR) = phimNM(1:2, ig, jFSR) * wttemp + srcmNM(1:2, ig, jFSR) * RTHREE
      END DO
    END DO  
  END DO
ELSE IF (ScatOd .EQ. 2) THEN
  DO iPin = PinIdxSt, PinIdxEd
    FsrIdxSt = Pin(iPin)%FsrIdxSt
    icel     = Pin(iPin)%Cell(iz)
    
    DO iFSR = 1, Cell(icel)%nFsr
      jFSR = FsrIdxSt + iFSR - 1
      
      DO ig = gb, ge
        wttemp = 1._8 / (xstNM(ig, jFSR) * Cell(icel)%vol(iFSR))
        
        phisNM(ig, jFSR) = phisNM(ig, jFSR) * wttemp + srcNM(ig, jFSR)
        
        phimNM(1:2, ig, jFSR) = phimNM(1:2, ig, jFSR) * wttemp + srcmNM(1:2, ig, jFSR) * rthree
        phimNM(3:5, ig, jFSR) = phimNM(3:5, ig, jFSR) * wttemp + srcmNM(3:5, ig, jFSR) * rfive
      END DO
    END DO  
  END DO
ELSE IF (ScatOd .EQ. 3) THEN
  DO iPin = PinIdxSt, PinIdxEd
    FsrIdxSt = Pin(iPin)%FsrIdxSt
    icel     = Pin(iPin)%Cell(iz)
    
    DO iFSR = 1, Cell(icel)%nFsr
      jFSR = FsrIdxSt + iFSR - 1
      
      DO ig = gb, ge
        wttemp = 1._8 / (xstNM(ig, jFSR) * Cell(icel)%vol(iFSR))
        
        phisNM(ig, jFSR) = phisNM(ig, jFSR) * wttemp + srcNM(ig, jFSR)
        
        phimNM(1:2, ig, jFSR) = phimNM(1:2, ig, jFSR) * wttemp + srcmNM(1:2, ig, jFSR) * rthree
        phimNM(3:5, ig, jFSR) = phimNM(3:5, ig, jFSR) * wttemp + srcmNM(3:5, ig, jFSR) * rfive
        phimNM(6:9, ig, jFSR) = phimNM(6:9, ig, jFSR) * wttemp + srcmNM(6:9, ig, jFSR) * rseven
      END DO
    END DO  
  END DO
END IF
! ----------------------------------------------------
! Geo.
NULLIFY (Asy)
NULLIFY (AsyInfo)
NULLIFY (Cell)
NULLIFY (Pin)

! Loc.
NULLIFY (xstNM)
NULLIFY (srcNM)
NULLIFY (phianm)
NULLIFY (SrcAngnm)

! Dcmp.
NULLIFY (DcmpAsyRay)
NULLIFY (DcmpAsyRayCount)
! ----------------------------------------------------

END SUBROUTINE RayTraceDcmp_Pn
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayPn_Dcmp(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, joutNM, lJout, iz, gb, ge, krot)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, Asy_Type, Cell_Type, AsyRayInfo_type, CoreRayInfo_Type, CellRayInfo_type, TrackingDat_Type, DcmpAsyRayInfo_Type

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

REAL, POINTER, DIMENSION(:,:,:,:) :: joutNM

LOGICAL :: lJout
INTEGER :: iz, gb, ge, krot
! ----------------------------------------------------
TYPE (Pin_Type),          POINTER, DIMENSION(:) :: Pin
TYPE (Asy_Type),          POINTER, DIMENSION(:) :: Asy
TYPE (Cell_Type),         POINTER, DIMENSION(:) :: Cell
TYPE (AsyRayInfo_type),   POINTER, DIMENSION(:) :: AsyRay
TYPE (CoreRayInfo_Type),  POINTER, DIMENSION(:) :: CoreRay

TYPE (CellRayInfo_Type), POINTER :: CellRay

REAL, POINTER, DIMENSION(:)         :: LenSeg
REAL, POINTER, DIMENSION(:,:)       :: xstNM, EXPA, EXPB, wtang
REAL, POINTER, DIMENSION(:,:,:)     :: PhiAngInNM, wtsurf
REAL, POINTER, DIMENSION(:,:,:,:)   :: SrcAngNM, phiaNM
REAL, POINTER, DIMENSION(:,:,:,:,:) :: DcmpPhiAngIn, DcmpPhiAngOut

INTEGER, POINTER, DIMENSION(:)   :: LocalFsrIdx, AsyRayList, DirList, AziList
INTEGER, POINTER, DIMENSION(:,:) :: AziMap

INTEGER :: mp(2), AziSvIdx(2)
INTEGER :: iAzi, ipol, iRotRay, iAsyRay, jAsyRay, iRay, iRaySeg, idir, jbeg, jend, jinc
INTEGER :: ipin, icel, iasy, ireg, isurf, ibcel, iceray, ig, ir, iPinRay
INTEGER :: nCoreRay, nAsyRay, nDcmpRay, nPinRay, nRaySeg, FsrIdxSt, nPolarAng, PhiAnginSvIdx, ExpAppIdx

REAL :: wt(10), wt2(10, 4)
REAL :: PhiAngOut(RayInfo%nPolarAngle, gb : ge)
REAL :: phid, tau, ExpApp

DATA mp /2, 1/
! ----------------------------------------------------

! Ray
nPolarAng = RayInfo%nPolarAngle
AsyRay   => RayInfo%AsyRay
CoreRay  => RayInfo%CoreRay

! Geo.
Asy  => CoreInfo%Asy
Pin  => CoreInfo%Pin
Cell => CoreInfo%CellInfo

! Tracking Dat
xstNM         => TrackingDat%xstNM
PhiAngInNM    => TrackingDat%PhiAngInNM
DcmpPhiAngIn  => TrackingDat%DcmpPhiAngIn
DcmpPhiAngOut => TrackingDat%DcmpPhiAngOut
EXPA          => TrackingDat%EXPA
EXPB          => TrackingDat%EXPB
wtang         => TrackingDat%wtang
wtsurf        => TrackingDat%wtsurf
AziMap        => TrackingDat%AziMap
phiaNM        => TrackingDat%phianm
SrcAngNM      => TrackingDat%SrcAngnm

! Dcmp.
nAsyRay     = DcmpAsyRay%nAsyRay
iRotRay     = DcmpAsyRay%iRotRay
iAsy        = DcmpAsyRay%iAsy
iRay        = DcmpAsyRay%iRay
AsyRayList => DcmpAsyRay%AsyRayList
DirList    => DcmpAsyRay%DirList
AziList    => DcmpAsyRay%AziList
! ----------------------------------------------------
PhiAngInSvIdx = RayInfo%PhiAngInSvIdx(iRotRay, krot)

IF (DcmpAsyRay%lRotRayBeg(krot)) THEN
  PhiAngOut = PhiAngInNM(1:nPolarAng, gb:ge, PhiAnginSvIdx)
ELSE
  PhiAngOut = DcmpPhiAngIn(1:nPolarAng, gb:ge, krot, iRay, iAsy)
END IF

IF (krot .EQ. 2) THEN
  jbeg = nAsyRay; jend = 1; jinc = -1
ELSE
  jbeg = 1; jend = nAsyRay; jinc = 1
END IF

DO iAsyRay = jbeg, jend, jinc
  jAsyRay = AsyRayList(iAsyRay)
  nPinRay = AsyRay(jAsyRay)%nCellRay
  iazi    = AziList(iAsyRay)
  idir    = DirList(iAsyRay)
  IF (krot .eq. 2) idir = mp(idir)
  
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
            joutNM(1, ig, isurf, ipin) = joutNM(1, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
            joutNM(3, ig, isurf, ipin) = joutNM(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
          END DO
        END DO
      END IF
      
      DO iRaySeg = 1, nRaySeg
        ireg = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
        
        DO ig = gb, ge
          tau = - LenSeg(iRaySeg) * xstNM(ig, ireg)
          
          ExpAppIdx = max(INT(tau), -40000)
          ExpAppIdx = min(0, ExpAppIdx)
          
          DO ipol = 1, nPolarAng
            ExpApp = ExpA(ExpAppIdx, ipol) * tau + EXPB(ExpAppIdx, ipol)
            
            phid = (PhiAngOut(ipol, ig) - SrcAngNM(ipol, ig, ireg, AziSvIdx(idir))) * ExpApp
            
            PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
            
            phiaNM(ipol, ig, ireg, AziSvIdx(idir)) = phiaNM(ipol, ig, ireg, AziSvIdx(idir)) + phid
          END DO
        END DO
      END DO
      
      IF (lJout) THEN
        isurf = AsyRay(jAsyRay)%PinRaySurf(2, iPinRay)
        
        DO ig = gb, ge
          DO ipol = 1, nPolarAng
            joutNM(2, ig, isurf, ipin) = joutNM(2, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
            joutNM(3, ig, isurf, ipin) = joutNM(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
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
            joutNM(1, ig, isurf, ipin) = joutNM(1, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
            joutNM(3, ig, isurf, ipin) = joutNM(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
          END DO
        END DO
      END IF
      
      DO iRaySeg = nRaySeg, 1, -1
        ireg = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
        
        DO ig = gb, ge
          tau = - LenSeg(iRaySeg) * xstNM(ig, ireg)
          
          ExpAppIdx = max(INT(tau), -40000)
          ExpAppIdx = min(0, ExpAppIdx)
          
          DO ipol = 1, nPolarAng
            ExpApp = EXPA(ExpAppIdx, ipol) * tau + EXPB(ExpAppIdx, ipol)
            
            phid = (PhiAngOut(ipol, ig) - SrcAngNM(ipol, ig, ireg, AziSvIdx(idir))) * ExpApp
            
            PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
            
            phiaNM(ipol, ig, ireg, AziSvIdx(idir)) = phiaNM(ipol, ig, ireg, AziSvIdx(idir)) + phid
          END DO
        END DO
      END DO
      
      IF (lJout) THEN
        isurf = AsyRay(jAsyRay)%PinRaySurf(1, iPinRay)
        
        DO ig = gb, ge
          DO ipol = 1, nPolarAng
            joutNM(2, ig, isurf, ipin) = joutNM(2, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
            joutNM(3, ig, isurf, ipin) = joutNM(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
          END DO
        END DO
      END IF
    END DO
  END IF
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
NULLIFY (srcAngNM)
NULLIFY (phiaNM)
NULLIFY (xstNM)
NULLIFY (EXPA)
NULLIFY (EXPB)
NULLIFY (wtang)
NULLIFY (wtsurf)
NULLIFY (PhiAngInNM)
NULLIFY (LocalFsrIdx)

! Dcmp.
NULLIFY (DcmpPhiAngIn)
NULLIFY (DcmpPhiAngOut)
NULLIFY (AsyRayList)
NULLIFY (DirList)
NULLIFY (AziList)
NULLIFY (AziMap)
! ----------------------------------------------------

END SUBROUTINE RecTrackRotRayPn_Dcmp
! ------------------------------------------------------------------------------------------------------------