#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmp_NM(RayInfo, CoreInfo, phisNM, PhiAngInNM, xstNM, srcNM, MocJoutNM, iz, gb, ge, lJout)

USE OMP_LIB
USE allocs
USE PARAM,       ONLY : ZERO, RED, BLACK, GREEN
USE TYPEDEF,     ONLY : RayInfo_Type, Coreinfo_type, Asy_Type, AsyInfo_Type, Pin_Type, Cell_Type, DcmpAsyRayInfo_Type
USE Moc_Mod,     ONLY : TrackingDat, DcmpPhiAngIn, DcmpPhiAngOut, DcmpColorAsy, RayTraceDcmp_OMP, DcmpGatherBoundaryFlux, DcmpScatterBoundaryFlux, DcmpLinkBoundaryFlux
USE PE_MOD,      ONLY : PE
USE CNTL,        ONLY : nTracerCntl
USE geom,        ONLY : nbd
USE itrcntl_mod, ONLY : itrcntl
USE HexData,     ONLY : hAsy, hLgc

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:,:)     :: phisNM, xstNM, srcNM
REAL, POINTER, DIMENSION(:,:,:)   :: PhiAngInNM
REAL, POINTER, DIMENSION(:,:,:,:) :: MocJoutNM

INTEGER :: iz, gb, ge
LOGICAL :: lJout
! ----------------------------------------------------
TYPE (AsyInfo_Type), POINTER, DIMENSION(:) :: AsyInfo
TYPE (Asy_Type),     POINTER, DIMENSION(:) :: Asy
TYPE (Cell_Type),    POINTER, DIMENSION(:) :: Cell
TYPE (Pin_Type),     POINTER, DIMENSION(:) :: Pin

TYPE (DcmpAsyRayInfo_Type), POINTER, DIMENSION(:,:) :: DcmpAsyRay

INTEGER, POINTER, DIMENSION(:) :: DcmpAsyRayCount

INTEGER :: ithr, nThr, iAsy, jAsy, ifsr, ibd, ixy, nxy, FsrIdxSt, icel, jfsr, ig, iAsyRay, krot, icolor, jcolor, ncolor, iit, PinSt, PinEd, FsrSt, FsrEd
LOGICAL :: lHex

INTEGER, PARAMETER :: AuxRec(2, 0:1) = [2, 1,  1, 2]
INTEGER, PARAMETER :: AuxHex(3, 0:2) = [3, 1, 2,  1, 2, 3,  2, 3, 1]
! ----------------------------------------------------

nxy      = CoreInfo%nxy
Asy     => CoreInfo%Asy
AsyInfo => CoreInfo%AsyInfo
Cell    => CoreInfo%CellInfo
Pin     => CoreInfo%Pin

DcmpAsyRay      => RayInfo%DcmpAsyRay
DcmpAsyRayCount => RayInfo%DcmpAsyRayCount

lHex = nTracerCntl%lHex

IF (lHex) THEN
  ncolor = 3; iit = mod(itrcntl%mocit, 3)
  
  IF (hLgc%l060) ncolor = 1
ELSE
  ncolor = 2; iit = mod(itrcntl%mocit, 2)
END IF

nthr = PE%nthread
CALL OMP_SET_NUM_THREADS(nThr)
! ----------------------------------------------------
DO ithr = 1, nThr
  TrackingDat(ithr)%srcNM => srcNM
  TrackingDat(ithr)%xstNM => xstNM
END DO

DcmpPhiAngOut(:, gb:ge, :, :, :) = ZERO
phisNM(gb:ge, :) = ZERO
IF (ljout) MocjoutNM(:, gb:ge, :, :) = ZERO
! ----------------------------------------------------
DO icolor = 1, ncolor
  IF (lHex) THEN
    jcolor = AuxHex(icolor, iit)
    
    IF (hLgc%l060) jcolor = icolor
  ELSE
    jcolor = AuxRec(icolor, iit)
  END IF
  
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) CALL DcmpScatterBoundaryFlux(RayInfo, PhiAngInNM, DcmpPhiAngIn)
#endif
  
  DO ithr = 1, nthr
    TrackingDat(ithr)%PhiAngInNM    => PhiAngInNM
    TrackingDat(ithr)%DcmpPhiAngIn  => DcmpPhiAngIn
    TrackingDat(ithr)%DcmpPhiAngOut => DcmpPhiAngOut
  END DO
  
  !$OMP PARALLEL PRIVATE(ithr, iAsy, jAsy, PinSt, PinEd, FsrSt, FsrEd, krot, iAsyRay, ig, ifsr, ixy, ibd)
  ithr = 1
  !$ ithr = omp_get_thread_num()+1
  !$OMP DO SCHEDULE(GUIDED)
  DO iAsy = 1, DcmpColorAsy(0, jcolor)
    jAsy = DcmpColorAsy(iAsy, jcolor)
    
    ! SET : Range
    IF (lHex) THEN
      PinSt = hAsy(jAsy)%PinIdxSt
      PinEd = hAsy(jAsy)%PinIdxSt + hAsy(jAsy)%nTotPin - 1
    ELSE
      PinSt = Asy(jAsy)%GlobalPinIdx(1)
      PinEd = Asy(jAsy)%GlobalPinIdx(AsyInfo(Asy(jAsy)%AsyType)%nxy)
    END IF
    
    FsrSt = Pin(PinSt)%FsrIdxSt
    FsrEd = Pin(PinEd)%FsrIdxSt + Cell(Pin(PinEd)%Cell(iz))%nFsr - 1
    
    ! ALLOC
    CALL dmalloc0(TrackingDat(ithr)%phisNM, gb, ge, FsrSt, FsrEd)
    IF (ljout) CALL dmalloc0(TrackingDat(ithr)%JoutNM, 1, 3, gb, ge, 1, nbd, PinSt, PinEd)
    
    ! RT
    DO krot = 1, 2
      DO iAsyRay = 1, DcmpAsyRayCount(jAsy)
        IF (lHex) THEN
          CALL HexTrackRotRayDcmp_NM(RayInfo, CoreInfo, TrackingDat(ithr), DcmpAsyRay(iAsyRay, jAsy), ljout, iz, gb, ge, krot)
        ELSE
          CALL RecTrackRotRayDcmp_NM(RayInfo, CoreInfo, TrackingDat(ithr), DcmpAsyRay(iAsyRay, jAsy), ljout, iz, gb, ge, krot)
        END IF
      END DO
    END DO
    
    ! GATHER
    DO ifsr = FsrSt, FsrEd
      DO ig = gb, ge
        phisNM(ig, ifsr) = phisNM(ig, ifsr) + TrackingDat(ithr)%phisNM(ig, ifsr)
      END DO
    END DO
    
    IF (ljout) THEN
      DO ixy = PinSt, PinEd
        DO ibd = 1, nbd
          DO ig = gb, ge
            MocjoutNM(:, ig, ibd, ixy) = MocjoutNM(:, ig, ibd, ixy) + TrackingDat(ithr)%joutNM(:, ig, ibd, ixy)
          END DO
        END DO
      END DO
    END IF
    
    ! FREE
    DEALLOCATE (TrackingDat(ithr)%phisNM)
    IF (ljout) DEALLOCATE (TrackingDat(ithr)%JoutNM)
  END DO
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL
  
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) CALL DcmpGatherBoundaryFlux(RayInfo, DcmpPhiAngOut)
#endif
  
  IF (PE%RTMASTER) CALL DcmpLinkBoundaryFlux(CoreInfo, RayInfo, PhiAngInNM, DcmpPhiAngIn, DcmpPhiAngOut, gb, ge, jcolor)
END DO
! ----------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ixy, FsrIdxSt, icel, ifsr, jfsr, ig)
!$OMP DO SCHEDULE(GUIDED)
DO ixy = 1, nxy
  FsrIdxSt = Pin(ixy)%FsrIdxSt
  icel     = Pin(ixy)%Cell(iz)
  
  DO ifsr = 1, Cell(icel)%nFsr
    jfsr = FsrIdxSt + ifsr - 1
    
    DO ig = gb, ge
      phisNM(ig, jfsr) = phisNM(ig, jfsr) / (xstNM(ig, jfsr) * Cell(icel)%vol(ifsr)) + srcNM(ig, jfsr)
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
! ----------------------------------------------------
NULLIFY (Cell)
NULLIFY (Pin)
NULLIFY (Asy)
NULLIFY (AsyInfo)
NULLIFY (DcmpAsyRay)
NULLIFY (DcmpAsyRayCount)
! ----------------------------------------------------

END SUBROUTINE RayTraceDcmp_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayDcmp_NM(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, ljout, iz, gb, ge, krot)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, Asy_Type, Cell_Type, AsyRayInfo_type, CellRayInfo_type, TrackingDat_Type, DcmpAsyRayInfo_Type

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: iz, gb, ge, krot
! ----------------------------------------------------
TYPE (Pin_Type),        POINTER, DIMENSION(:) :: Pin
TYPE (Asy_Type),        POINTER, DIMENSION(:) :: Asy
TYPE (Cell_Type),       POINTER, DIMENSION(:) :: Cell
TYPE (AsyRayInfo_type), POINTER, DIMENSION(:) :: AsyRay

TYPE (CellRayInfo_Type),  POINTER :: CellRay

REAL, POINTER, DIMENSION(:)         :: LenSeg
REAL, POINTER, DIMENSION(:,:)       :: phisNM, srcNM, xstNM, EXPA, EXPB, wtang
REAL, POINTER, DIMENSION(:,:,:)     :: PhiAngInNM, wtsurf
REAL, POINTER, DIMENSION(:,:,:,:)   :: joutNM
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
phisNM        => TrackingDat%phisNM
srcNM         => TrackingDat%srcNM
xstNM         => TrackingDat%xstNM
PhiAngInNM    => TrackingDat%PhiAngInNM
joutNM        => TrackingDat%joutNM
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
NULLIFY (phisNM)
NULLIFY (srcNM)
NULLIFY (xstNM)
NULLIFY (joutNM)
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

END SUBROUTINE RecTrackRotRayDcmp_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayDcmp_NM(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, ljout, iz, gb, ge, krot)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, TrackingDat_Type, DcmpAsyRayInfo_Type, AziAngleInfo_Type
USE HexType, ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexCoreRay, Type_HexRotRay
USE HexData, ONLY : hAsy, haRay, hAsyTypInfo

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: iz, gb, ge, krot
! ----------------------------------------------------
TYPE (Pin_Type),          POINTER, DIMENSION(:) :: Pin
TYPE (AziAngleInfo_Type), POINTER, DIMENSION(:) :: AziAng
  
TYPE (Type_HexAsyRay), POINTER :: haRay_Loc
TYPE (Type_HexCelRay), POINTER :: CelRay_Loc

REAL, POINTER, DIMENSION(:,:)       :: phisNM, srcNM, xstNM, EXPA, EXPB, wtang, hwt
REAL, POINTER, DIMENSION(:,:,:)     :: PhiAngInNM
REAL, POINTER, DIMENSION(:,:,:,:)   :: joutNM
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
phisNM        => TrackingDat%phisNM
srcNM         => TrackingDat%srcNM
xstNM         => TrackingDat%xstNM
PhiAngInNM    => TrackingDat%PhiAngInNM
joutNM        => TrackingDat%joutNM
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
NULLIFY (EXPA)
NULLIFY (EXPB)
NULLIFY (wtang)
NULLIFY (hwt)
NULLIFY (phisNM)
NULLIFY (srcNM)
NULLIFY (xstNM)
NULLIFY (PhiAngInNM)
NULLIFY (joutNM)

! Dcmp.
NULLIFY (AsyRayList)
NULLIFY (DirList)
NULLIFY (AziList)
NULLIFY (DcmpPhiAngIn)
NULLIFY (DcmpPhiAngOut)
! ----------------------------------------------------

END SUBROUTINE HexTrackRotRayDcmp_NM
! ------------------------------------------------------------------------------------------------------------