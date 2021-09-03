#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmp_GM(RayInfo, CoreInfo, phis1g, PhiAngIn1g, xst1g, src1g, MocJout1g, iz, lJout)

USE OMP_LIB
USE PARAM,       ONLY : ZERO
USE TYPEDEF,     ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, Cell_Type
USE Moc_Mod,     ONLY : TrackingDat, DcmpPhiAngIn1g, DcmpPhiAngOut1g, DcmpAsyClr, DcmpGatherBndyFlux1g, DcmpScatterBndyFlux1g, DcmpLinkBndyFlux1g, RtDcmpThr_GM, nClr, setDcmpClr
USE PE_MOD,      ONLY : PE
USE CNTL,        ONLY : nTracerCntl
USE itrcntl_mod, ONLY : itrcntl
USE HexData,     ONLY : hLgc

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:)     :: phis1g, xst1g, src1g
REAL, POINTER, DIMENSION(:,:)   :: PhiAngIn1g
REAL, POINTER, DIMENSION(:,:,:) :: MocJout1g

INTEGER :: iz
LOGICAL :: lJout
! ----------------------------------------------------
TYPE (Pin_Type),  POINTER, DIMENSION(:) :: Pin
TYPE (Cell_Type), POINTER, DIMENSION(:) :: Cell

INTEGER :: ithr, nThr, iAsy, jAsy, ixy, nxy, icel, ifsr, jfsr, FsrIdxSt, iClr, jClr
LOGICAL :: lHex, lAFSS

INTEGER, PARAMETER :: AuxRec(2, 0:1) = [2, 1,  1, 2]
INTEGER, PARAMETER :: AuxHex(3, 0:2) = [3, 1, 2,  1, 2, 3,  2, 3, 1]
! ----------------------------------------------------

nxy   = CoreInfo%nxy
Cell => CoreInfo%CellInfo
Pin  => CoreInfo%Pin

lHex  = nTracerCntl%lHex
lAFSS = nTracerCntl%lAFSS

nthr = PE%nthread
CALL OMP_SET_NUM_THREADS(nThr)
! ----------------------------------------------------
DO ithr = 1, nThr
  TrackingDat(ithr)%src1g => src1g
  TrackingDat(ithr)%xst1g => xst1g
END DO

DcmpPhiAngOut1g = ZERO

phis1g = ZERO
IF (ljout) Mocjout1g = ZERO
! ----------------------------------------------------
DO iClr = 1, nClr
  jClr = setDcmpClr(lHex, hLgc%l060, iClr, itrcntl%mocit)
  
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) CALL DcmpScatterBndyFlux1g(RayInfo, PhiAngIn1g, DcmpPhiAngIn1g)
#endif
  
  DO ithr = 1, nthr
    TrackingDat(ithr)%PhiAngIn1g      => PhiAngIn1g
    TrackingDat(ithr)%DcmpPhiAngIn1g  => DcmpPhiAngIn1g
    TrackingDat(ithr)%DcmpPhiAngOut1g => DcmpPhiAngOut1g
  END DO
  
  !$OMP PARALLEL PRIVATE(ithr, iAsy, jAsy)
  ithr = 1
  !$ ithr = omp_get_thread_num()+1
  !$OMP DO SCHEDULE(GUIDED)
  DO iAsy = 1, DcmpAsyClr(0, jClr)
    jAsy = DcmpAsyClr(iAsy, jClr)
    
    CALL RtDcmpThr_GM(RayInfo, CoreInfo, TrackingDat(ithr), phis1g, MocJout1g, jAsy, iz, lJout, lHex, lAFSS)
  END DO
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL
  
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) CALL DcmpGatherBndyFlux1g(RayInfo, DcmpPhiAngOut1g)
#endif
  
  IF (PE%RTMASTER) CALL DcmpLinkBndyFlux1g(CoreInfo, RayInfo, PhiAngIn1g, DcmpPhiAngIn1g, DcmpPhiAngOut1g, jClr)
END DO
! ----------------------------------------------------
!$OMP PARALLEL PRIVATE(ixy, FsrIdxSt, icel, ifsr, jfsr)
!$OMP DO SCHEDULE(GUIDED)
DO ixy = 1, nxy
  FsrIdxSt = Pin(ixy)%FsrIdxSt
  icel     = Pin(ixy)%Cell(iz)
  
  DO ifsr = 1, Cell(icel)%nFsr
    jfsr = FsrIdxSt + ifsr - 1
    
    phis1g(jfsr) = phis1g(jfsr) / (xst1g(jfsr) * Cell(icel)%vol(ifsr)) + src1g(jfsr)
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

NULLIFY (Pin)
NULLIFY (Cell)
! ----------------------------------------------------

END SUBROUTINE RayTraceDcmp_GM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RtDcmpThr_GM(RayInfo, CoreInfo, TrackingLoc, phis1g, MocJout1g, jAsy, iz, lJout, lHex, lAFSS)

USE allocs
USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type, Cell_Type, Pin_Type, DcmpAsyRayInfo_Type
USE geom,    ONLY : nbd
USE MOC_MOD, ONLY : HexTrackRotRayDcmp_GM, wtang
USE HexData, ONLY : hAsy

IMPLICIT NONE

TYPE (RayInfo_Type)     :: RayInfo
TYPE (CoreInfo_Type)    :: CoreInfo
TYPE (TrackingDat_Type) :: TrackingLoc

REAL, POINTER, DIMENSION(:)     :: phis1g
REAL, POINTER, DIMENSION(:,:,:) :: MocJout1g

INTEGER :: jAsy, iz
LOGICAL :: lJout, lHex, lAFSS
! ----------------------------------------------------
TYPE (Cell_Type), POINTER, DIMENSION(:) :: Cell
TYPE (Pin_Type),  POINTER, DIMENSION(:) :: Pin

TYPE (DcmpAsyRayInfo_Type), POINTER, DIMENSION(:,:) :: DcmpAsyRay

INTEGER, POINTER, DIMENSION(:) :: DcmpAsyRayCount

INTEGER :: kRot, iAsyRay, ifsr, ixy, ibd, iazi, ipol, PinSt, PinEd, FsrSt, FsrEd, nPol, nAzi
! ----------------------------------------------------

Cell => CoreInfo%Cellinfo
Pin  => CoreInfo%Pin

nPol             = RayInfo%nPolarAngle
nAzi             = RayInfo%nAziAngle
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

! ALLOC
CALL dmalloc0(TrackingLoc%phis1g, FsrSt, FsrEd)
IF (ljout) CALL dmalloc0(TrackingLoc%Jout1g, 1, 3, 1, nbd, PinSt, PinEd)
IF (lAFSS) CALL dmalloc0(TrackingLoc%phia1g1, 1, nPol, FsrSt, FsrEd, 1, nAzi)
IF (lAFSS) CALL dmalloc0(TrackingLoc%phia1g2, 1, nPol, FsrSt, FsrEd, 1, nAzi)

! RT
IF (lHex) THEN
  DO krot = 1, 2
    DO iAsyRay = 1, DcmpAsyRayCount(jAsy)
      CALL HexTrackRotRayDcmp_GM(RayInfo, CoreInfo, TrackingLoc, DcmpAsyRay(iAsyRay, jAsy), lJout, iz, krot, lAFSS)
    END DO
  END DO
ELSE
  DO krot = 1, 2
    DO iAsyRay = 1, DcmpAsyRayCount(jAsy)
      !CALL RecTrackRotRayDcmp_GM(RayInfo, CoreInfo, TrackingLoc, DcmpAsyRay(iAsyRay, jAsy), lJout, iz, krot, lAFSS)
    END DO
  END DO
END IF

! GATHER
IF (.NOT. lAFSS) THEN
  DO ifsr = FsrSt, FsrEd
    phis1g(ifsr) = phis1g(ifsr) + TrackingLoc%phis1g(ifsr)
  END DO
ELSE
  DO iazi = 1, nAzi
    DO ifsr = FsrSt, FsrEd
      DO ipol = 1, nPol
        phis1g(ifsr) = phis1g(ifsr) + wtang(ipol, iazi) * (TrackingLoc%phia1g1(ipol, ifsr, iazi) + TrackingLoc%phia1g2(ipol, ifsr, iazi))
      END DO
    END DO
  END DO
END IF

IF (lJout) THEN
  DO ixy = PinSt, PinEd
    DO ibd = 1, nbd
      Mocjout1g(:, ibd, ixy) = Mocjout1g(:, ibd, ixy) + TrackingLoc%jout1g(:, ibd, ixy)
    END DO
  END DO
END IF

! FREE
DEALLOCATE (TrackingLoc%phis1g)
IF (lJout) DEALLOCATE (TrackingLoc%Jout1g)
IF (lAFSS) DEALLOCATE (TrackingLoc%phia1g1)
IF (lAFSS) DEALLOCATE (TrackingLoc%phia1g2)

NULLIFY (Cell)
NULLIFY (Pin)
NULLIFY (DcmpAsyRay)
NULLIFY (DcmpAsyRayCount)
! ----------------------------------------------------

END SUBROUTINE RtDcmpThr_GM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayDcmp_GM(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, ljout, iz, krot, lAFSS)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type, DcmpAsyRayInfo_Type, Pin_Type, AziAngleInfo_Type
USE HexType, ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexCoreRay, Type_HexRotRay
USE HexData, ONLY : hAsy, haRay, hAsyTypInfo

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

LOGICAL, INTENT(IN) :: ljout, lAFSS
INTEGER, INTENT(IN) :: iz, krot
! ----------------------------------------------------
TYPE (Pin_Type),          POINTER, DIMENSION(:) :: Pin
TYPE (AziAngleInfo_Type), POINTER, DIMENSION(:) :: AziAng
  
TYPE (Type_HexAsyRay), POINTER :: haRay_Loc
TYPE (Type_HexCelRay), POINTER :: CelRay_Loc

REAL, POINTER, DIMENSION(:)       :: phis1g, src1g, xst1g
REAL, POINTER, DIMENSION(:,:)     :: PhiAngIn1g, EXPA, EXPB, wtang, hwt
REAL, POINTER, DIMENSION(:,:,:)   :: jout1g, phia1g
REAL, POINTER, DIMENSION(:,:,:,:) :: DcmpPhiAngIn1g, DcmpPhiAngOut1g

INTEGER, POINTER, DIMENSION(:) :: AsyRayList, DirList, AziList

INTEGER :: mp(2)
INTEGER :: iazi, ipol, irotray, iasyray, iray, irayseg, idir, icel, iasy, ifsr, isurf, jbeg, jend, jinc, imray, ipray
INTEGER :: ipst, iped, ipinc, isfst, isfed, isgst, isged, isginc, iAsyTyp, iGeoTyp, icBss, jhPin, jcBss
INTEGER :: nAsyRay, nPinRay, nRaySeg, FsrIdxSt, nPolarAng, PhiAnginSvIdx, ExpAppIdx

REAL :: wtazi(10), locwt(10), loccs, locsn

REAL, DIMENSION(RayInfo%nPolarAngle) :: PhiAngOut1g
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

DcmpPhiAngIn1g  => TrackingDat%DcmpPhiAngIn1g
DcmpPhiAngOut1g => TrackingDat%DcmpPhiAngOut1g

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
phis1g     => TrackingDat%phis1g
src1g      => TrackingDat%src1g
xst1g      => TrackingDat%xst1g
PhiAngIn1g => TrackingDat%PhiAngIn1g
jout1g     => TrackingDat%jout1g
wtang      => TrackingDat%wtang
EXPA       => TrackingDat%EXPA
EXPB       => TrackingDat%EXPB
hwt        => TrackingDat%hwt

! Iter.
IF (DcmpAsyRay%lRotRayBeg(krot)) THEN
  PhiAngOut1g(1:nPolarAng) = PhiAngIn1g    (1:nPolarAng, PhiAnginSvIdx)
ELSE
  PhiAngOut1g(1:nPolarAng) = DcmpPhiAngIn1g(1:nPolarAng, krot, iRay, iAsy)
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
    ipst = 1; iped = nPinRay; ipinc = 1;  isfst = 1; isfed = 2; IF (lAFSS) phia1g => TrackingDat%phia1g1
  ELSE
    iped = 1; ipst = nPinRay; ipinc = -1; isfed = 1; isfst = 2; IF (lAFSS) phia1g => TrackingDat%phia1g2
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
        
        jout1g(1, isurf, jhpin) = jout1g(1, isurf, jhpin) + PhiAngOut1g(ipol) * wtazi(ipol)
        jout1g(3, isurf, jhpin) = jout1g(3, isurf, jhpin) + PhiAngOut1g(ipol) * wtsurf
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
      
      tau = -CelRay_Loc%SegLgh(iRaySeg) * xst1g(ifsr) ! Optimum Length
      
      ExpAppIdx = max(INT(tau), -40000)
      ExpAppIdx = min(0, ExpAppIdx)
      
      DO ipol = 1, nPolarAng
        ExpApp = ExpA(ExpAppIdx, ipol) * tau + ExpB(ExpAppIdx, ipol)
        
        phid = (PhiAngOut1g(ipol) - src1g(ifsr)) * ExpApp
        
        PhiAngOut1g(ipol) = PhiAngOut1g(ipol) - phid
        
        IF (lAFSS) THEN
          phia1g(ipol, ifsr, iazi) = phia1g(ipol, ifsr, iazi) + phid
        ELSE
          phis1g(ifsr) = phis1g(ifsr) + wtazi(ipol) * phid
        END IF
      END DO
    END DO
    
    ! Surface : Out-going
    IF (lJout) THEN
      isurf = CelRay_Loc%hSufIdx(isfed) ! y : Big
      
      DO ipol = 1, nPolarAng
        wtsurf = locwt(ipol) / abs(loccs * CelRay_Loc%hsn(isfed) - locsn * CelRay_Loc%hcs(isfed))
        
        jout1g(2, isurf, jhpin) = jout1g(2, isurf, jhpin) + PhiAngOut1g(ipol) * wtazi(ipol)
        jout1g(3, isurf, jhpin) = jout1g(3, isurf, jhpin) + PhiAngOut1g(ipol) * wtsurf
      END DO
    END IF
  END DO
END DO

DcmpPhiAngOut1g(1:nPolarAng, krot, iRay, iAsy) = PhiAngOut1g(1:nPolarAng)
! ----------------------------------------------------
! Ray
NULLIFY (AziAng)
NULLIFY (haRay_Loc)
NULLIFY (CelRay_Loc)

! Geo.
NULLIFY (Pin)

! Loc.
NULLIFY (phis1g)
NULLIFY (src1g)
NULLIFY (xst1g)
NULLIFY (PhiAngIn1g)
NULLIFY (jout1g)
NULLIFY (wtang)
NULLIFY (EXPA)
NULLIFY (EXPB)
NULLIFY (hwt)

IF (lAFSS) NULLIFY (phia1g)

! Dcmp.
NULLIFY (AsyRayList)
NULLIFY (DirList)
NULLIFY (AziList)
NULLIFY (DcmpPhiAngIn1g)
NULLIFY (DcmpPhiAngOut1g)
! ----------------------------------------------------

END SUBROUTINE HexTrackRotRayDcmp_GM
! ------------------------------------------------------------------------------------------------------------