#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmp_NM(RayInfo, CoreInfo, phisNg, PhiAngInNg, xstNg, srcNg, MocJoutNg, iz, gb, ge, lJout)

USE OMP_LIB
USE PARAM,       ONLY : ZERO
USE TYPEDEF,     ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, Cell_Type
USE Moc_Mod,     ONLY : TrackingDat, DcmpPhiAngInNg, DcmpPhiAngOutNg, DcmpAsyClr, DcmpGatherBndyFluxNg, DcmpScatterBndyFluxNg, DcmpLinkBndyFluxNg, RtDcmpThr_NM, nClr, setDcmpClr
USE PE_MOD,      ONLY : PE
USE CNTL,        ONLY : nTracerCntl
USE itrcntl_mod, ONLY : itrcntl

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:,:)     :: phisNg, xstNg, srcNg
REAL, POINTER, DIMENSION(:,:,:)   :: PhiAngInNg
REAL, POINTER, DIMENSION(:,:,:,:) :: MocJoutNg

INTEGER :: iz, gb, ge
LOGICAL :: lJout
! ----------------------------------------------------
TYPE (Pin_Type),  POINTER, DIMENSION(:) :: Pin
TYPE (Cell_Type), POINTER, DIMENSION(:) :: Cell

INTEGER :: ithr, nThr, iAsy, jAsy, ixy, nxy, icel, ifsr, jfsr, FsrIdxSt, ig, iClr, jClr
LOGICAL :: lHex, lRGB, lAFSS
! ----------------------------------------------------

nxy   = CoreInfo%nxy
Cell => CoreInfo%CellInfo
Pin  => CoreInfo%Pin

lHex  = nTracerCntl%lHex
lRGB  = nTracerCntl%lRGB
lAFSS = nTracerCntl%lAFSS

nthr = PE%nthread
CALL OMP_SET_NUM_THREADS(nThr)
! ----------------------------------------------------
DO ithr = 1, nThr
  TrackingDat(ithr)%srcNg => srcNg
  TrackingDat(ithr)%xstNg => xstNg
END DO

DcmpPhiAngOutNg(:, gb:ge, :, :, :) = ZERO

phisNg(gb:ge, :) = ZERO
IF (ljout) MocJoutNg(:, gb:ge, :, :) = ZERO
! ----------------------------------------------------
DO iClr = 1, nClr
  jClr = setDcmpClr(lHex, lRGB, iClr, itrcntl%mocit)
  
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) CALL DcmpScatterBndyFluxNg(RayInfo, PhiAngInNg, DcmpPhiAngInNg)
#endif
  
  DO ithr = 1, nthr
    TrackingDat(ithr)%PhiAngInNg      => PhiAngInNg
    TrackingDat(ithr)%DcmpPhiAngInNg  => DcmpPhiAngInNg
    TrackingDat(ithr)%DcmpPhiAngOutNg => DcmpPhiAngOutNg
  END DO
  
  !$OMP PARALLEL PRIVATE(ithr, iAsy, jAsy)
  ithr = 1
  !$ ithr = omp_get_thread_num()+1
  !$OMP DO SCHEDULE(GUIDED)
  DO iAsy = 1, DcmpAsyClr(0, jClr)
    jAsy = DcmpAsyClr(iAsy, jClr)
    
    CALL RtDcmpThr_NM(RayInfo, CoreInfo, TrackingDat(ithr), phisNg, MocJoutNg, jAsy, iz, gb, ge, lJout, lHex, lAFSS)
  END DO
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL
  
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) CALL DcmpGatherBndyFluxNg(RayInfo, DcmpPhiAngOutNg)
#endif
  
  IF (PE%RTMASTER) CALL DcmpLinkBndyFluxNg(CoreInfo, RayInfo, PhiAngInNg, DcmpPhiAngInNg, DcmpPhiAngOutNg, gb, ge, jClr)
END DO
! ----------------------------------------------------
!$OMP PARALLEL PRIVATE(ixy, FsrIdxSt, icel, ifsr, jfsr, ig)
!$OMP DO SCHEDULE(GUIDED)
DO ixy = 1, nxy
  FsrIdxSt = Pin(ixy)%FsrIdxSt
  icel     = Pin(ixy)%Cell(iz)
  
  DO ifsr = 1, Cell(icel)%nFsr
    jfsr = FsrIdxSt + ifsr - 1
    
    DO ig = gb, ge
      phisNg(ig, jfsr) = phisNg(ig, jfsr) / (xstNg(ig, jfsr) * Cell(icel)%vol(ifsr)) + srcNg(ig, jfsr)
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

NULLIFY (Pin)
NULLIFY (Cell)
! ----------------------------------------------------

END SUBROUTINE RayTraceDcmp_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RtDcmpThr_NM(RayInfo, CoreInfo, TrackingLoc, phisNg, MocJoutNg, jAsy, iz, gb, ge, lJout, lHex, lAFSS)

USE allocs
USE PARAM,   ONLY : ZERO
USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type, Pin_Type, Cell_Type, DcmpAsyRayInfo_Type
USE geom,    ONLY : nbd
USE MOC_MOD, ONLY : RecTrackRotRayDcmp_NM, HexTrackRotRayDcmp_NM, DcmpAziRay, wtang
USE HexData, ONLY : hAsy, hLgc

IMPLICIT NONE

TYPE (RayInfo_Type)     :: RayInfo
TYPE (CoreInfo_Type)    :: CoreInfo
TYPE (TrackingDat_Type) :: TrackingLoc

REAL, POINTER, DIMENSION(:,:)     :: phisNg
REAL, POINTER, DIMENSION(:,:,:,:) :: MocJoutNg

INTEGER :: jAsy, iz, gb, ge
LOGICAL :: lJout, lHex, lAFSS
! ----------------------------------------------------
TYPE (Pin_Type),  POINTER, DIMENSION(:) :: Pin
TYPE (Cell_Type), POINTER, DIMENSION(:) :: Cell

TYPE (DcmpAsyRayInfo_Type), POINTER, DIMENSION(:,:) :: DcmpAsyRay

INTEGER, POINTER, DIMENSION(:) :: DcmpAsyRayCount

INTEGER :: kRot, iAsyRay, jAsyRay, ifsr, ig, ixy, ibd, iazi, ipol, PinSt, PinEd, FsrSt, FsrEd, nAzi, nPol
! ----------------------------------------------------

Pin  => CoreInfo%Pin
Cell => CoreInfo%Cellinfo

nAzi             = RayInfo%nAziAngle
nPol             = RayInfo%nPolarAngle
DcmpAsyRay      => RayInfo%DcmpAsyRay
DcmpAsyRayCount => RayInfo%DcmpAsyRayCount

IF (lHex) THEN
  PinSt = hAsy(jAsy)%PinIdxSt
  PinEd = hAsy(jAsy)%PinIdxSt + hAsy(jAsy)%nTotPin - 1
ELSE
  PinSt = CoreInfo%Asy(jAsy)%GlobalPinIdx(1)
  PinEd = CoreInfo%Asy(jAsy)%GlobalPinIdx(CoreInfo%AsyInfo(CoreInfo%Asy(jAsy)%AsyType)%nxy)
END IF

FsrSt = Pin(PinSt)%FsrIdxSt
FsrEd = Pin(PinEd)%FsrIdxSt + Cell(Pin(PinEd)%Cell(iz))%nFsr - 1

CALL dmalloc0(TrackingLoc%phisNg, gb, ge, FsrSt, FsrEd)
IF (ljout) CALL dmalloc0(TrackingLoc%JoutNg, 1, 3, gb, ge, 1, nbd, PinSt, PinEd)
! ----------------------------------------------------
IF (hLgc%lNoRef .AND. lAFSS) THEN
  CALL dmalloc0(TrackingLoc%phiaNg1, 1, nPol, gb, ge, Fsrst, FsrEd, 1, 1)
  CALL dmalloc0(TrackingLoc%phiaNg2, 1, nPol, gb, ge, Fsrst, FsrEd, 1, 1)
  
  DO iazi = 1, nAzi
    TrackingLoc%phiaNg1 = ZERO
    TrackingLoc%phiaNg2 = ZERO
    
    IF (lHex) THEN
      DO krot = 1, 2
        DO iAsyRay = 1, DcmpAziRay(0, iazi, jAsy)
          jAsyRay = DcmpAziRay(iAsyRay, iazi, jAsy)
          
          CALL HexTrackRotRayDcmp_NM(RayInfo, CoreInfo, TrackingLoc, DcmpAsyRay(jAsyRay, jAsy), lJout, iz, gb, ge, krot, lAFSS)
        END DO
      END DO
    ELSE
      DO krot = 1, 2
        DO iAsyRay = 1, DcmpAziRay(0, iazi, jAsy)
          jAsyRay = DcmpAziRay(iAsyRay, iazi, jAsy)
          
          !CALL RecTrackRotRayDcmp_NM(RayInfo, CoreInfo, TrackingLoc, DcmpAsyRay(jAsyRay, jAsy), lJout, iz, gb, ge, krot, lAFSS)
        END DO
      END DO
    END IF
    
    DO ifsr = FsrSt, FsrEd
      DO ig = gb, ge
        DO ipol = 1, nPol
          phisNg(ig, ifsr) = phisNg(ig, ifsr) + wtang(ipol, iazi) * (TrackingLoc%phiaNg1(ipol, ig, ifsr, 1) + TrackingLoc%phiaNg2(ipol, ig, ifsr, 1)) ! NOTICE
        END DO
      END DO
    END DO
  END DO
ELSE
  IF (lAFSS) CALL dmalloc0(TrackingLoc%phiaNg1, 1, nPol, gb, ge, Fsrst, FsrEd, 1, nAzi)
  IF (lAFSS) CALL dmalloc0(TrackingLoc%phiaNg2, 1, nPol, gb, ge, Fsrst, FsrEd, 1, nAzi)
  
  IF (lHex) THEN
    DO krot = 1, 2
      DO iAsyRay = 1, DcmpAsyRayCount(jAsy)
        CALL HexTrackRotRayDcmp_NM(RayInfo, CoreInfo, TrackingLoc, DcmpAsyRay(iAsyRay, jAsy), lJout, iz, gb, ge, krot, lAFSS)
      END DO
    END DO
  ELSE
    DO krot = 1, 2
      DO iAsyRay = 1, DcmpAsyRayCount(jAsy)
        !CALL RecTrackRotRayDcmp_NM(RayInfo, CoreInfo, TrackingLoc, DcmpAsyRay(iAsyRay, jAsy), lJout, iz, gb, ge, krot, lAFSS)
      END DO
    END DO
  END IF
  
  IF (lAFSS) THEN
    DO iazi = 1, nAzi
      DO ifsr = FsrSt, FsrEd
        DO ig = gb, ge
          DO ipol = 1, nPol
            phisNg(ig, ifsr) = phisNg(ig, ifsr) + wtang(ipol, iazi) * (TrackingLoc%phiaNg1(ipol, ig, ifsr, iazi) + TrackingLoc%phiaNg2(ipol, ig, ifsr, iazi))
          END DO
        END DO
      END DO
    END DO
  END IF
END IF
! ----------------------------------------------------
IF (.NOT. lAFSS) THEN
  DO ifsr = FsrSt, FsrEd
    DO ig = gb, ge
      phisNg(ig, ifsr) = phisNg(ig, ifsr) + TrackingLoc%phisNg(ig, ifsr)
    END DO
  END DO
END IF

IF (lJout) THEN
  DO ixy = PinSt, PinEd
    DO ibd = 1, nbd
      DO ig = gb, ge
        MocJoutNg(:, ig, ibd, ixy) = MocJoutNg(:, ig, ibd, ixy) + TrackingLoc%JoutNg(:, ig, ibd, ixy)
      END DO
    END DO
  END DO
END IF
! ----------------------------------------------------
DEALLOCATE (TrackingLoc%phisNg)
IF (lJout) DEALLOCATE (TrackingLoc%JoutNg)
IF (lAFSS) DEALLOCATE (TrackingLoc%phiaNg1)
IF (lAFSS) DEALLOCATE (TrackingLoc%phiaNg2)

NULLIFY (Pin)
NULLIFY (Cell)
NULLIFY (DcmpAsyRay)
NULLIFY (DcmpAsyRayCount)
! ----------------------------------------------------

END SUBROUTINE RtDcmpThr_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayDcmp_NM(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, ljout, iz, gb, ge, krot, lAFSS)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, Asy_Type, Cell_Type, AsyRayInfo_type, CellRayInfo_type, TrackingDat_Type, DcmpAsyRayInfo_Type

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

LOGICAL, INTENT(IN) :: ljout, lAFSS
INTEGER, INTENT(IN) :: iz, gb, ge, krot
! ----------------------------------------------------
TYPE (Pin_Type),        POINTER, DIMENSION(:) :: Pin
TYPE (Asy_Type),        POINTER, DIMENSION(:) :: Asy
TYPE (Cell_Type),       POINTER, DIMENSION(:) :: Cell
TYPE (AsyRayInfo_type), POINTER, DIMENSION(:) :: AsyRay

TYPE (CellRayInfo_Type),  POINTER :: CellRay

REAL, POINTER, DIMENSION(:)         :: LenSeg
REAL, POINTER, DIMENSION(:,:)       :: phisNg, srcNg, xstNg, EXPA, EXPB, wtang
REAL, POINTER, DIMENSION(:,:,:)     :: PhiAngInNg, wtsurf
REAL, POINTER, DIMENSION(:,:,:,:)   :: JoutNg, phiaNg
REAL, POINTER, DIMENSION(:,:,:,:,:) :: DcmpPhiAngInNg, DcmpPhiAngOutNg

INTEGER, POINTER, DIMENSION(:) :: LocalFsrIdx, AsyRayList, DirList, AziList

INTEGER :: mp(2)
INTEGER :: iazi, ipol, irotray, iasyray, iceray, iray, irayseg, idir, ipin, icel, ibcel, iasy, ifsr, isurf, ig, jbeg, jend, jinc, imray, ipray
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
phisNg          => TrackingDat%phisNg
srcNg           => TrackingDat%srcNg
xstNg           => TrackingDat%xstNg
PhiAngInNg      => TrackingDat%PhiAngInNg
JoutNg          => TrackingDat%JoutNg
DcmpPhiAngInNg  => TrackingDat%DcmpPhiAngInNg
DcmpPhiAngOutNg => TrackingDat%DcmpPhiAngOutNg
EXPA            => TrackingDat%EXPA
EXPB            => TrackingDat%EXPB
wtang           => TrackingDat%wtang
wtsurf          => TrackingDat%wtsurf

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
  DO ig = gb, ge
    PhiAngOut(1:nPolarAng, ig) = PhiAngInNg(1:nPolarAng, PhiAnginSvIdx, ig)
  END DO
ELSE
  PhiAngOut(1:nPolarAng, gb:ge) = DcmpPhiAngInNg(1:nPolarAng, gb:ge, krot, iRay, iAsy)
END IF
  
IF (krot .EQ. 2) THEN
  jbeg = nAsyRay; jend = 1; jinc = -1
ELSE
  jend = nAsyRay; jbeg = 1; jinc = 1
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
    ipst = 1; iped = nPinRay; ipinc = 1;  isfst = 1; isfed = 2; IF (lAFSS) phiaNg => TrackingDat%phiaNg1
  ELSE
    iped = 1; ipst = nPinRay; ipinc = -1; isfed = 1; isfst = 2; IF (lAFSS) phiaNg => TrackingDat%phiaNg2
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
          JoutNg(1, ig, isurf, ipin) = JoutNg(1, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
          JoutNg(3, ig, isurf, ipin) = JoutNg(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
        END DO
      END DO
    END IF
    
    IF (idir .EQ. 1) THEN
      isgst = 1; isged = nRaySeg; isginc = 1
    ELSE
      isged = 1; isgst = nRaySeg; isginc = -1
    END IF
    
    DO iRaySeg = isgst, isged, isginc
      ifsr = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
      
      DO ig = gb, ge
        tau = -LenSeg(iRaySeg) * xstNg(ig, ifsr)
        
        ExpAppIdx = max(INT(tau), -40000)
        ExpAppIdx = min(0, ExpAppIdx)
        
        DO ipol = 1, nPolarAng
          ExpApp = ExpA(ExpAppIdx, ipol) * tau + ExpB(ExpAppIdx, ipol)
          
          phid = (PhiAngOut(ipol, ig) - srcNg(ig, ifsr)) * ExpApp
          
          PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
          
          IF (lAFSS) THEN
            phiaNg(ipol, ig, ifsr, iazi) = phiaNg(ipol, ig, ifsr, iazi)  + phid
          ELSE
            phisNg(ig, ifsr) = phisNg(ig, ifsr) + wt(ipol) * phid
          END IF
        END DO
      END DO
    END DO
    
    IF (lJout) THEN
      isurf = AsyRay(iAsyRay)%PinRaySurf(isfed, ipray)
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          JoutNg(2, ig, isurf, ipin) = JoutNg(2, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
          JoutNg(3, ig, isurf, ipin) = JoutNg(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
        END DO
      END DO
    END IF
  END DO
END DO

DcmpPhiAngOutNg(1:nPolarAng, gb:ge, krot, iRay, iAsy) = PhiAngOut
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
NULLIFY (phisNg)
NULLIFY (srcNg)
NULLIFY (xstNg)
NULLIFY (JoutNg)
NULLIFY (EXPA)
NULLIFY (EXPB)
NULLIFY (wtang)
NULLIFY (wtsurf)
NULLIFY (PhiAngInNg)
NULLIFY (LocalFsrIdx)

IF (lAFSS) NULLIFY (phiaNg)

! Dcmp.
NULLIFY (DcmpPhiAngInNg)
NULLIFY (DcmpPhiAngOutNg)
NULLIFY (AsyRayList)
NULLIFY (DirList)
NULLIFY (AziList)
! ----------------------------------------------------

END SUBROUTINE RecTrackRotRayDcmp_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayDcmp_NM(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, ljout, iz, gb, ge, krot, lAFSS)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type, DcmpAsyRayInfo_Type, Pin_Type, AziAngleInfo_Type
USE HexType, ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexCoreRay, Type_HexRotRay
USE HexData, ONLY : hAsy, haRay, hAsyTypInfo, hLgc

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

LOGICAL, INTENT(IN) :: ljout, lAFSS
INTEGER, INTENT(IN) :: iz, gb, ge, krot
! ----------------------------------------------------
TYPE (Pin_Type),          POINTER, DIMENSION(:) :: Pin
TYPE (AziAngleInfo_Type), POINTER, DIMENSION(:) :: AziAng
  
TYPE (Type_HexAsyRay), POINTER :: haRay_Loc
TYPE (Type_HexCelRay), POINTER :: CelRay_Loc

REAL, POINTER, DIMENSION(:,:)       :: phisNg, srcNg, xstNg, EXPA, EXPB, wtang, hwt
REAL, POINTER, DIMENSION(:,:,:)     :: PhiAngInNg
REAL, POINTER, DIMENSION(:,:,:,:)   :: JoutNg, phiaNg
REAL, POINTER, DIMENSION(:,:,:,:,:) :: DcmpPhiAngInNg, DcmpPhiAngOutNg

INTEGER, POINTER, DIMENSION(:) :: AsyRayList, DirList, AziList

INTEGER :: mp(2)
INTEGER :: iazi, jazi, ipol, irotray, iasyray, iray, irayseg, idir, icel, iasy, ifsr, isurf, ig, jbeg, jend, jinc, imray, ipray
INTEGER :: ipst, iped, ipinc, isfst, isfed, isgst, isged, isginc, iAsyTyp, iGeoTyp, icBss, jhPin, jcBss
INTEGER :: nAsyRay, nPinRay, nRaySeg, FsrIdxSt, nPolarAng, PhiAnginSvIdx, ExpAppIdx

REAL :: wtazi(10), locwt(10), loccs, locsn

REAL, DIMENSION(RayInfo%nPolarAngle, gb:ge) :: PhiAngOut
REAL :: phid, tau, ExpApp, wtsurf

DATA mp /2, 1/
! ----------------------------------------------------

! Dcmp.
nAsyRay     = DcmpAsyRay%nAsyRay
iRotRay     = DcmpAsyRay%iRotRay
iAsy        = DcmpAsyRay%iAsy
iRay        = DcmpAsyRay%iRay
AsyRayList => DcmpAsyRay%AsyRayList
DirList    => DcmpAsyRay%DirList
AziList    => DcmpAsyRay%AziList

DcmpPhiAngInNg  => TrackingDat%DcmpPhiAngInNg
DcmpPhiAngOutNg => TrackingDat%DcmpPhiAngOutNg

! Ray
nPolarAng     = RayInfo%nPolarAngle
AziAng       => RayInfo%AziAngle
PhiAngInSvIdx = RayInfo%PhiAngInSvIdx(iRotRay, krot)

! Geo.
Pin => CoreInfo%Pin

iAsyTyp = hAsy(iAsy)%AsyTyp
iGeoTyp = hAsy(iAsy)%GeoTyp
icBss   = hAsyTypInfo(iAsyTyp)%iBss

! Loc.
phisNg     => TrackingDat%phisNg
srcNg      => TrackingDat%srcNg
xstNg      => TrackingDat%xstNg
JoutNg     => TrackingDat%JoutNg
PhiAngInNg => TrackingDat%PhiAngInNg
wtang      => TrackingDat%wtang
EXPA       => TrackingDat%EXPA
EXPB       => TrackingDat%EXPB
hwt        => TrackingDat%hwt

! Iter.
IF (DcmpAsyRay%lRotRayBeg(krot)) THEN
  DO ig = gb, ge
    PhiAngOut(1:nPolarAng, ig) = PhiAngInNg(1:nPolarAng, PhiAnginSvIdx, ig)
  END DO
ELSE
  PhiAngOut(1:nPolarAng, gb:ge) = DcmpPhiAngInNg(1:nPolarAng, gb:ge, krot, iRay, iAsy)
END IF

IF (krot .EQ. 1) THEN
  jbeg = 1; jend = nAsyRay; jinc = 1
ELSE
  jend = 1; jbeg = nAsyRay; jinc = -1
END IF
! ----------------------------------------------------
DO imray = jbeg, jend, jinc
  iAsyRay = AsyRayList(imray)
  iazi    = AziList   (imray)
  jazi    = iazi
  idir    = DirList   (imray)
  IF (hLgc%lNoRef) jazi = 1
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
    ipst = 1; iped = nPinRay; ipinc = 1;  isfst = 1; isfed = 2; phiaNg => TrackingDat%phiaNg1
  ELSE
    iped = 1; ipst = nPinRay; ipinc = -1; isfed = 1; isfst = 2; phiaNg => TrackingDat%phiaNg2 
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
          JoutNg(1, ig, isurf, jhpin) = JoutNg(1, ig, isurf, jhpin) + PhiAngOut(ipol, ig) * wtazi(ipol)
          JoutNg(3, ig, isurf, jhpin) = JoutNg(3, ig, isurf, jhpin) + PhiAngOut(ipol, ig) * wtsurf
        END DO
      END DO
    END IF
    
    ! Iter. : FSR
    nRaySeg = CelRay_Loc%nSegRay
    
    IF (idir .EQ. 1) THEN
      isgst = 1; isged = nRaySeg; isginc = 1
    ELSE
      isged = 1; isgst = nRaySeg; isginc = -1
    END IF
    
    DO iRaySeg = isgst, isged, isginc
      ifsr = CelRay_Loc%MshIdx(iRaySeg) + Pin(jhPin)%FsrIdxSt - 1
      
      DO ig = gb, ge
        tau = -CelRay_Loc%SegLgh(iRaySeg) * xstNg(ig, ifsr) ! Optimum Length
        
        ExpAppIdx = max(INT(tau), -40000)
        ExpAppIdx = min(0, ExpAppIdx)
        
        DO ipol = 1, nPolarAng
          ExpApp = ExpA(ExpAppIdx, ipol) * tau + ExpB(ExpAppIdx, ipol)
          
          phid = (PhiAngOut(ipol, ig) - srcNg(ig, ifsr)) * ExpApp
          
          PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
          
          IF (lAFSS) THEN
            phiaNg(ipol, ig, ifsr, jazi) = phiaNg(ipol, ig, ifsr, jazi) + phid
          ELSE
            phisNg(ig, ifsr) = phisNg(ig, ifsr) + wtazi(ipol) * phid
          END IF
        END DO
      END DO
    END DO
    
    ! Surface : Out-going
    IF (lJout) THEN
      isurf = CelRay_Loc%hSufIdx(isfed) ! y : Big
      
      DO ipol = 1, nPolarAng
        wtsurf = locwt(ipol) / abs(loccs * CelRay_Loc%hsn(isfed) - locsn * CelRay_Loc%hcs(isfed))
        
        DO ig = gb, ge
          JoutNg(2, ig, isurf, jhpin) = JoutNg(2, ig, isurf, jhpin) + PhiAngOut(ipol, ig) * wtazi(ipol)
          JoutNg(3, ig, isurf, jhpin) = JoutNg(3, ig, isurf, jhpin) + PhiAngOut(ipol, ig) * wtsurf
        END DO
      END DO
    END IF
  END DO
END DO

DcmpPhiAngOutNg(1:nPolarAng, gb:ge, krot, iRay, iAsy) = PhiAngOut(1:nPolarAng, gb:ge)
! ----------------------------------------------------
! Ray
NULLIFY (AziAng)
NULLIFY (haRay_Loc)
NULLIFY (CelRay_Loc)

! Geo.
NULLIFY (Pin)

! Loc.
NULLIFY (phisNg)
NULLIFY (srcNg)
NULLIFY (xstNg)
NULLIFY (PhiAngInNg)
NULLIFY (JoutNg)
NULLIFY (wtang)
NULLIFY (EXPA)
NULLIFY (EXPB)
NULLIFY (hwt)

IF (lAFSS) NULLIFY (phiaNg)

! Dcmp.
NULLIFY (AsyRayList)
NULLIFY (DirList)
NULLIFY (AziList)
NULLIFY (DcmpPhiAngInNg)
NULLIFY (DcmpPhiAngOutNg)
! ----------------------------------------------------

END SUBROUTINE HexTrackRotRayDcmp_NM
! ------------------------------------------------------------------------------------------------------------