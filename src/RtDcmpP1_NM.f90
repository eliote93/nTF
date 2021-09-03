#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmpP1_NM(RayInfo, CoreInfo, phisNg, phimNg, PhiAngInNg, xstNg, srcNg, srcmNg, MocJoutNg, iz, gb, ge, lJout)

USE OMP_LIB
USE PARAM,       ONLY : ZERO, ONE
USE TYPEDEF,     ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, Cell_Type
USE Moc_Mod,     ONLY : TrackingDat, DcmpPhiAngInNg, DcmpPhiAngOutNg, DcmpAsyClr, DcmpGatherBndyFluxNg, DcmpScatterBndyFluxNg, DcmpLinkBndyFluxNg, RtDcmpP1Thr_NM, nClr, setDcmpClr
USE PE_MOD,      ONLY : PE
USE CNTL,        ONLY : nTracerCntl
USE itrcntl_mod, ONLY : itrcntl
USE HexData,     ONLY : hLgc

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:,:)     :: phisNg, xstNg, srcNg
REAL, POINTER, DIMENSION(:,:,:)   :: PhiAngInNg, phimNg, srcmNg
REAL, POINTER, DIMENSION(:,:,:,:) :: MocJoutNg

INTEGER :: iz, gb, ge
LOGICAL :: lJout
! ----------------------------------------------------
TYPE (Pin_Type),  POINTER, DIMENSION(:) :: Pin
TYPE (Cell_Type), POINTER, DIMENSION(:) :: Cell

INTEGER :: ithr, nThr, iAsy, jAsy, ixy, nxy, icel, ifsr, jfsr, FsrIdxSt, ig, iClr, jClr, ScatOd
LOGICAL :: lHex
REAL :: wttmp

INTEGER, PARAMETER :: AuxRec(2, 0:1) = [2, 1,  1, 2]
INTEGER, PARAMETER :: AuxHex(3, 0:2) = [3, 1, 2,  1, 2, 3,  2, 3, 1]
! ----------------------------------------------------

nxy   = CoreInfo%nxy
Cell => CoreInfo%CellInfo
Pin  => CoreInfo%Pin

lHex   = nTracerCntl%lHex
ScatOd = nTracerCntl%ScatOd

nthr = PE%nthread
CALL OMP_SET_NUM_THREADS(nThr)
! ----------------------------------------------------
DO ithr = 1, nThr
  TrackingDat(ithr)%srcNg => srcNg
  TrackingDat(ithr)%xstNg => xstNg
END DO

DcmpPhiAngOutNg(:, gb:ge, :, :, :) = ZERO

phisNg(gb:ge, :) = ZERO
IF (lJout) MocJoutNg(:, gb:ge, :, :) = ZERO
phimNg(:, :, gb:ge) = ZERO
! ----------------------------------------------------
DO iClr = 1, nClr
  jClr = setDcmpClr(lHex, hLgc%l060, iClr, itrcntl%mocit)
  
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
    
    CALL RtDcmpP1Thr_NM(RayInfo, CoreInfo, TrackingDat(ithr), phisNg, phimNg, srcmNg, MocJoutNg, jAsy, iz, gb, ge, ScatOd, lJout, lHex)
  END DO
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL
  
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) CALL DcmpGatherBndyFluxNg(RayInfo, DcmpPhiAngOutNg)
#endif
  
  IF (PE%RTMASTER) CALL DcmpLinkBndyFluxNg(CoreInfo, RayInfo, PhiAngInNg, DcmpPhiAngInNg, DcmpPhiAngOutNg, gb, ge, jClr)
END DO
! ----------------------------------------------------
!$OMP PARALLEL PRIVATE(ixy, FsrIdxSt, icel, ifsr, jfsr, ig, wttmp)
!$OMP DO SCHEDULE(GUIDED)
DO ixy = 1, nxy
  FsrIdxSt = Pin(ixy)%FsrIdxSt
  icel     = Pin(ixy)%Cell(iz)
  
  DO ifsr = 1, Cell(icel)%nFsr
    jfsr = FsrIdxSt + ifsr - 1
    
    DO ig = gb, ge
      wttmp = ONE / (xstNg(ig, jfsr) * Cell(icel)%vol(ifsr))
      
      phisNg(ig, jfsr) = phisNg(ig, jfsr) * wttmp + srcNg(ig, jfsr)
      
      phimNg(1:2, jfsr, ig) = phimNg(1:2, jfsr, ig) * wttmp + srcmNg(1:2, jfsr, ig)
      
      IF (ScatOd .GE. 2) phimNg(3:5, jfsr, ig) = phimNg(3:5, jfsr, ig) * wttmp + srcmNg(3:5, jfsr, ig)
      IF (ScatOd .EQ. 3) phimNg(6:9, jfsr, ig) = phimNg(6:9, jfsr, ig) * wttmp + srcmNg(6:9, jfsr, ig)
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

NULLIFY (Pin)
NULLIFY (Cell)
! ----------------------------------------------------

END SUBROUTINE RayTraceDcmpP1_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RtDcmpP1Thr_NM(RayInfo, CoreInfo, TrackingLoc, phisNg, phimNg, srcmNg, MocJoutNg, jAsy, iz, gb, ge, ScatOd, lJout, lHex)

USE ALLOCS
USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type, Pin_Type, Cell_Type, DcmpAsyRayInfo_Type
USE geom,    ONLY : nbd
USE Moc_Mod, ONLY : RecTrackRotRayDcmpP1_NM, HexTrackRotRayDcmpP1_NM, Comp, DcmpAziRay
USE HexData, ONLY : hAsy, hLgc

IMPLICIT NONE

TYPE (RayInfo_Type)     :: RayInfo
TYPE (CoreInfo_Type)    :: CoreInfo
TYPE (TrackingDat_Type) :: TrackingLoc

REAL, POINTER, DIMENSION(:,:)     :: phisNg
REAL, POINTER, DIMENSION(:,:,:)   :: phimNg, srcmNg
REAL, POINTER, DIMENSION(:,:,:,:) :: MocJoutNg

INTEGER :: jAsy, iz, gb, ge, ScatOd
LOGICAL :: lJout, lHex
! ----------------------------------------------------
TYPE (Cell_Type), POINTER, DIMENSION(:) :: Cell
TYPE (Pin_Type),  POINTER, DIMENSION(:) :: Pin

TYPE (DcmpAsyRayInfo_Type), POINTER, DIMENSION(:,:) :: DcmpAsyRay

REAL, POINTER, DIMENSION(:,:)     :: srcNg
REAL, POINTER, DIMENSION(:,:,:,:) :: SrcAngNg1, SrcAngNg2

INTEGER, POINTER, DIMENSION(:) :: DcmpAsyRayCount

INTEGER :: PinSt, PinEd, FsrSt, FsrEd, kRot, iAsyRay, ifsr, ig, ixy, ibd
INTEGER :: nAziAng, nPolarAng, nOd, iOd, iazi, ipol, jAsyRay
REAL :: srctmp
! ----------------------------------------------------

Cell => CoreInfo%Cellinfo
Pin  => CoreInfo%Pin

nAziAng          = RayInfo%nAziAngle
nPolarAng        = RayInfo%nPolarAngle
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

SELECT CASE (ScatOd)
CASE (1); nod = 2
CASE (2); nod = 5
CASE (3); nod = 9
END SELECT
! ----------------------------------------------------
CALL dmalloc0(TrackingLoc%phisNg, gb, ge, FsrSt, FsrEd)
IF (ljout) CALL dmalloc0(TrackingLoc%JoutNg, 1, 3, gb, ge, 1, nbd, PinSt, PinEd)

CALL dmalloc0(TrackingLoc%phimNg, 1, nOd, FsrSt, FsrEd, gb, ge)

srcNg => TrackingLoc%srcNg
! ----------------------------------------------------
IF (hLgc%lNoRef) THEN
  CALL dmalloc0(TrackingLoc%SrcAngNg1, gb, ge, 1, nPolarAng, FsrSt, FsrEd, 1, 1)
  CALL dmalloc0(TrackingLoc%SrcAngNg2, gb, ge, 1, nPolarAng, FsrSt, FsrEd, 1, 1)
  
  DO iazi = 1, nAziAng
    SrcAngNg1 => TrackingLoc%SrcAngNg1
    SrcAngNg2 => TrackingLoc%SrcAngNg2
    
    DO ifsr = FsrSt, FsrEd
      DO ipol = 1, nPolarAng
        DO ig = gb, ge
          SrcAngNg1(ig, ipol, ifsr, 1) = srcNg(ig, ifsr)
          SrcAngNg2(ig, ipol, ifsr, 1) = srcNg(ig, ifsr)
          
          srctmp = comp(1, ipol, iazi) * srcmNg(1, ifsr, ig) + comp(2, ipol, iazi) * srcmNg(2, ifsr, ig)
          
          SrcAngNg1(ig, ipol, ifsr, 1) = SrcAngNg1(ig, ipol, ifsr, 1) + srctmp
          SrcAngNg2(ig, ipol, ifsr, 1) = SrcAngNg2(ig, ipol, ifsr, 1) - srctmp
          
          IF (ScatOd .LT. 2) CYCLE
          
          srctmp = comp(3, ipol, iazi) * srcmNg(3, ifsr, ig) + comp(4, ipol, iazi) * srcmNg(4, ifsr, ig) + comp(5, ipol, iazi) * srcmNg(5, ifsr, ig)
          
          SrcAngNg1(ig, ipol, ifsr, 1) = SrcAngNg1(ig, ipol, ifsr, 1) + srctmp
          SrcAngNg2(ig, ipol, ifsr, 1) = SrcAngNg2(ig, ipol, ifsr, 1) + srctmp
          
          IF (ScatOd .LT. 3) CYCLE
          
          srctmp = comp(6, ipol, iazi) * srcmNg(6, ifsr, ig) + comp(7, ipol, iazi) * srcmNg(7, ifsr, ig) + comp(8, ipol, iazi) * srcmNg(8, ifsr, ig) + comp(9, ipol, iazi) * srcmNg(9, ifsr, ig)
          
          SrcAngNg1(ig, ipol, ifsr, 1) = SrcAngNg1(ig, ipol, ifsr, 1) + srctmp
          SrcAngNg2(ig, ipol, ifsr, 1) = SrcAngNg2(ig, ipol, ifsr, 1) - srctmp
        END DO
      END DO
    END DO
    
    IF (lHex) THEN
      DO krot = 1, 2
        DO iAsyRay = 1, DcmpAziRay(0, iazi, jAsy)
          jAsyRay = DcmpAziRay(iAsyRay, iazi, jAsy)
          
          CALL HexTrackRotRayDcmpP1_NM(RayInfo, CoreInfo, TrackingLoc, DcmpAsyRay(jAsyRay, jAsy), lJout, iz, gb, ge, krot, ScatOd)
        END DO
      END DO
    ELSE
      DO krot = 1, 2
        DO iAsyRay = 1, DcmpAziRay(0, iazi, jAsy)
          jAsyRay = DcmpAziRay(iAsyRay, iazi, jAsy)
          
          !CALL RecTrackRotRayDcmpP1_NM(RayInfo, CoreInfo, TrackingLoc, DcmpAsyRay(jAsyRay, jAsy), lJout, iz, gb, ge, krot, ScatOd)
        END DO
      END DO
    END IF
  END DO
ELSE
  CALL dmalloc0(TrackingLoc%SrcAngNg1, gb, ge, 1, nPolarAng, FsrSt, FsrEd, 1, nAziAng)
  CALL dmalloc0(TrackingLoc%SrcAngNg2, gb, ge, 1, nPolarAng, FsrSt, FsrEd, 1, nAziAng)
  
  SrcAngNg1 => TrackingLoc%SrcAngNg1
  SrcAngNg2 => TrackingLoc%SrcAngNg2
  
  DO iazi = 1, nAziAng
    DO ifsr = FsrSt, FsrEd
      DO ipol = 1, nPolarAng
        DO ig = gb, ge
          SrcAngNg1(ig, ipol, ifsr, iazi) = srcNg(ig, ifsr)
          SrcAngNg2(ig, ipol, ifsr, iazi) = srcNg(ig, ifsr)
          
          srctmp = comp(1, ipol, iazi) * srcmNg(1, ifsr, ig) + comp(2, ipol, iazi) * srcmNg(2, ifsr, ig)
          
          SrcAngNg1(ig, ipol, ifsr, iazi) = SrcAngNg1(ig, ipol, ifsr, iazi) + srctmp
          SrcAngNg2(ig, ipol, ifsr, iazi) = SrcAngNg2(ig, ipol, ifsr, iazi) - srctmp
          
          IF (ScatOd .LT. 2) CYCLE
          
          srctmp = comp(3, ipol, iazi) * srcmNg(3, ifsr, ig) + comp(4, ipol, iazi) * srcmNg(4, ifsr, ig) + comp(5, ipol, iazi) * srcmNg(5, ifsr, ig)
          
          SrcAngNg1(ig, ipol, ifsr, iazi) = SrcAngNg1(ig, ipol, ifsr, iazi) + srctmp
          SrcAngNg2(ig, ipol, ifsr, iazi) = SrcAngNg2(ig, ipol, ifsr, iazi) + srctmp
          
          IF (ScatOd .LT. 3) CYCLE
          
          srctmp = comp(6, ipol, iazi) * srcmNg(6, ifsr, ig) + comp(7, ipol, iazi) * srcmNg(7, ifsr, ig) + comp(8, ipol, iazi) * srcmNg(8, ifsr, ig) + comp(9, ipol, iazi) * srcmNg(9, ifsr, ig)
          
          SrcAngNg1(ig, ipol, ifsr, iazi) = SrcAngNg1(ig, ipol, ifsr, iazi) + srctmp
          SrcAngNg2(ig, ipol, ifsr, iazi) = SrcAngNg2(ig, ipol, ifsr, iazi) - srctmp
        END DO
      END DO
    END DO
  END DO
  
  IF (lHex) THEN
    DO krot = 1, 2
      DO iAsyRay = 1, DcmpAsyRayCount(jAsy)
        CALL HexTrackRotRayDcmpP1_NM(RayInfo, CoreInfo, TrackingLoc, DcmpAsyRay(iAsyRay, jAsy), lJout, iz, gb, ge, krot, ScatOd)
      END DO
    END DO
  ELSE
    DO krot = 1, 2
      DO iAsyRay = 1, DcmpAsyRayCount(jAsy)
        !CALL RecTrackRotRayDcmpP1_NM(RayInfo, CoreInfo, TrackingLoc, DcmpAsyRay(iAsyRay, jAsy), lJout, iz, gb, ge, krot, ScatOd)
      END DO
    END DO
  END IF
END IF
! ----------------------------------------------------
DO ifsr = FsrSt, FsrEd
  DO ig = gb, ge
    phisNg(ig, ifsr) = phisNg(ig, ifsr) + TrackingLoc%phisNg(ig, ifsr)
    
    DO iod = 1, nod
      phimNg(iod, ifsr, ig) = phimNg(iod, ifsr, ig) + TrackingLoc%phimNg(iod, ifsr, ig)
    END DO
  END DO
END DO

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

DEALLOCATE (TrackingLoc%phimNg)
DEALLOCATE (TrackingLoc%SrcAngNg1)
DEALLOCATE (TrackingLoc%SrcAngNg2)

NULLIFY (Cell)
NULLIFY (Pin)
NULLIFY (DcmpAsyRay)
NULLIFY (DcmpAsyRayCount)

NULLIFY (srcNg)
NULLIFY (SrcAngNg1)
NULLIFY (SrcAngNg2)
! ----------------------------------------------------

END SUBROUTINE RtDcmpP1Thr_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayDcmpP1_NM(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, ljout, iz, gb, ge, krot, ScatOd)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type, DcmpAsyRayInfo_Type, Pin_Type, AziAngleInfo_Type
USE HexType, ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexCoreRay, Type_HexRotRay
USE HexData, ONLY : hAsy, haRay, hAsyTypInfo, hLgc

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: iz, gb, ge, krot, ScatOd
! ----------------------------------------------------
TYPE (Pin_Type),          POINTER, DIMENSION(:) :: Pin
TYPE (AziAngleInfo_Type), POINTER, DIMENSION(:) :: AziAng
  
TYPE (Type_HexAsyRay), POINTER :: haRay_Loc
TYPE (Type_HexCelRay), POINTER :: CelRay_Loc

REAL, POINTER, DIMENSION(:,:)       :: phisNg, srcNg, xstNg, EXPA, EXPB, wtang, hwt
REAL, POINTER, DIMENSION(:,:,:)     :: PhiAngInNg, phimNg, LocMwt
REAL, POINTER, DIMENSION(:,:,:,:)   :: JoutNg, LocSrc
REAL, POINTER, DIMENSION(:,:,:,:,:) :: DcmpPhiAngInNg, DcmpPhiAngOutNg

INTEGER, POINTER, DIMENSION(:) :: AsyRayList, DirList, AziList

INTEGER :: mp(2)
INTEGER :: iazi, jazi, ipol, irotray, iasyray, iray, irayseg, idir, icel, iasy, iFSR, isurf, ig, jbeg, jend, jinc, imray, ipray
INTEGER :: ipst, iped, ipinc, isfst, isfed, isgst, isged, isginc, iAsyTyp, iGeoTyp, icBss, jhPin, jcBss, iOd, nOd
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

! P1
phimNg => TrackingDat%phimNg

SELECT CASE (ScatOd)
CASE (1); nOd = 2
CASE (2); nOd = 5
CASE (3); nOd = 9
END SELECT
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
    ipst = 1; iped = nPinRay; ipinc = 1;  isfst = 1; isfed = 2; LocMwt => TrackingDat%mwt
  ELSE
    iped = 1; ipst = nPinRay; ipinc = -1; isfed = 1; isfst = 2; LocMwt => TrackingDat%mwt2
  END IF
  
  IF (idir .EQ. 1) THEN
    LocSrc => TrackingDat%SrcAngNg1
  ELSE
    LocSrc => TrackingDat%SrcAngNg2
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
      iFSR = CelRay_Loc%MshIdx(iRaySeg) + Pin(jhPin)%FsrIdxSt - 1
      
      DO ig = gb, ge
        tau = -CelRay_Loc%SegLgh(iRaySeg) * xstNg(ig, iFSR) ! Optimum Length
        
        ExpAppIdx = max(INT(tau), -40000)
        ExpAppIdx = min(0, ExpAppIdx)
        
        DO ipol = 1, nPolarAng
          ExpApp = ExpA(ExpAppIdx, ipol) * tau + ExpB(ExpAppIdx, ipol)
          
          phid = (PhiAngOut(ipol, ig) - LocSrc(ig, ipol, ifsr, jazi)) * ExpApp
          
          PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
          
          phisNg(ig, iFSR) = phisNg(ig, iFSR) + wtazi(ipol) * phid
          
          DO iod = 1, nod
            phimNg(iod, ifsr, ig) = phimNg(iod, ifsr, ig) + LocMwt(iod, ipol, iazi) * phid ! NOTICE
          END DO
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

! Dcmp.
NULLIFY (AsyRayList)
NULLIFY (DirList)
NULLIFY (AziList)
NULLIFY (DcmpPhiAngInNg)
NULLIFY (DcmpPhiAngOutNg)

! P1
NULLIFY (phimNg)
NULLIFY (LocMwt)
NULLIFY (LocSrc)
! ----------------------------------------------------

END SUBROUTINE HexTrackRotRayDcmpP1_NM
! ------------------------------------------------------------------------------------------------------------