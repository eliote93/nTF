#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmpP1_GM(RayInfo, CoreInfo, phis1g, phim1g, PhiAngIn1g, xst1g, src1g, srcm1g, MocJout1g, iz, lJout)

USE OMP_LIB
USE PARAM,       ONLY : ZERO, ONE, RTHREE, RFIVE, RSEVEN
USE TYPEDEF,     ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, Cell_Type
USE Moc_Mod,     ONLY : TrackingDat, DcmpPhiAngIn1g, DcmpPhiAngOut1g, DcmpAsyClr, DcmpGatherBndyFlux1g, DcmpScatterBndyFlux1g, DcmpLinkBndyFlux1g, RtDcmpP1Thr_GM
USE PE_MOD,      ONLY : PE
USE CNTL,        ONLY : nTracerCntl
USE itrcntl_mod, ONLY : itrcntl
USE HexData,     ONLY : hLgc

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:)     :: phis1g, xst1g, src1g
REAL, POINTER, DIMENSION(:,:)   :: PhiAngIn1g, phim1g, srcm1g
REAL, POINTER, DIMENSION(:,:,:) :: MocJout1g

INTEGER :: iz
LOGICAL :: lJout
! ----------------------------------------------------
TYPE (Pin_Type),  POINTER, DIMENSION(:) :: Pin
TYPE (Cell_Type), POINTER, DIMENSION(:) :: Cell

INTEGER :: ithr, nThr, iAsy, jAsy, ixy, nxy, icel, ifsr, jfsr, FsrIdxSt, iClr, jClr, nClr, iit, ScatOd
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

IF (lHex) THEN
  nClr = 3; iit = mod(itrcntl%mocit, 3)
  
  IF (hLgc%l060) nClr = 1
ELSE
  nClr = 2; iit = mod(itrcntl%mocit, 2)
END IF

nthr = PE%nthread
CALL OMP_SET_NUM_THREADS(nThr)
! ----------------------------------------------------
DO ithr = 1, nThr
  TrackingDat(ithr)%src => src1g
  TrackingDat(ithr)%xst => xst1g
END DO

DcmpPhiAngOut1g = ZERO

phis1g = ZERO
IF (ljout) Mocjout1g = ZERO
phim1g = ZERO
! ----------------------------------------------------
DO iClr = 1, nClr
  IF (lHex) THEN
    jClr = AuxHex(iClr, iit)
    
    IF (hLgc%l060) jClr = iClr
  ELSE
    jClr = AuxRec(iClr, iit)
  END IF
  
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) CALL DcmpScatterBndyFlux1g(RayInfo, PhiAngIn1g, DcmpPhiAngIn1g)
#endif
  
  DO ithr = 1, nthr
    TrackingDat(ithr)%PhiAngIn        => PhiAngIn1g
    TrackingDat(ithr)%DcmpPhiAngIn1g  => DcmpPhiAngIn1g
    TrackingDat(ithr)%DcmpPhiAngOut1g => DcmpPhiAngOut1g
  END DO
  
  !$OMP PARALLEL PRIVATE(ithr, iAsy, jAsy)
  ithr = 1
  !$ ithr = omp_get_thread_num()+1
  !$OMP DO SCHEDULE(GUIDED)
  DO iAsy = 1, DcmpAsyClr(0, jClr)
    jAsy = DcmpAsyClr(iAsy, jClr)
    
    CALL RtDcmpP1Thr_GM(RayInfo, CoreInfo, TrackingDat(ithr), phis1g, phim1g, srcm1g, MocJout1g, jAsy, iz, lJout, lHex, ScatOd)
  END DO
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL
  
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) CALL DcmpGatherBndyFlux1g(RayInfo, DcmpPhiAngOut1g)
#endif
  
  IF (PE%RTMASTER) CALL DcmpLinkBndyFlux1g(CoreInfo, RayInfo, PhiAngIn1g, DcmpPhiAngIn1g, DcmpPhiAngOut1g, jClr)
END DO
! ----------------------------------------------------
!$OMP PARALLEL PRIVATE(ixy, FsrIdxSt, icel, ifsr, jfsr, wttmp)
!$OMP DO SCHEDULE(GUIDED)
DO ixy = 1, nxy
  FsrIdxSt = Pin(ixy)%FsrIdxSt
  icel     = Pin(ixy)%Cell(iz)
  
  DO ifsr = 1, Cell(icel)%nFsr
    jfsr = FsrIdxSt + ifsr - 1
    
    wttmp = ONE / (xst1g(jfsr) * Cell(icel)%vol(ifsr))
    
    phis1g(jfsr) = phis1g(jfsr) * wttmp + src1g(jfsr)
    
    phim1g(1:2, jfsr) = phim1g(1:2, jfsr) * wttmp + srcm1g(1:2, jfsr) * RTHREE
    
    IF (ScatOd .LT. 2) CYCLE
    
    phim1g(3:5, jfsr) = phim1g(3:5, jfsr) * wttmp + srcm1g(3:5, jfsr) * RFIVE
    
    IF (ScatOd .LT. 3) CYCLE
    
    phim1g(6:9, jfsr) = phim1g(6:9, jfsr) * wttmp + srcm1g(6:9, jfsr) * RSEVEN
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

NULLIFY (Pin)
NULLIFY (Cell)
! ----------------------------------------------------

END SUBROUTINE RayTraceDcmpP1_GM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RtDcmpP1Thr_GM(RayInfo, CoreInfo, TrackingLoc, phis1g, phim1g, srcm1g, MocJout1g, jAsy, iz, lJout, lHex, ScatOd)

USE allocs
USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type, Cell_Type, Pin_Type, DcmpAsyRayInfo_Type
USE geom,    ONLY : nbd
USE MOC_MOD, ONLY : HexTrackRotRayDcmpP1_GM, Comp, DcmpPhiAngIn1g, DcmpPhiAngOut1g, DcmpAziRay
USE HexData, ONLY : hAsy, hLgc

IMPLICIT NONE

TYPE (RayInfo_Type)     :: RayInfo
TYPE (CoreInfo_Type)    :: CoreInfo
TYPE (TrackingDat_Type) :: TrackingLoc

REAL, POINTER, DIMENSION(:)     :: phis1g
REAL, POINTER, DIMENSION(:,:)   :: phim1g, srcm1g
REAL, POINTER, DIMENSION(:,:,:) :: MocJout1g

INTEGER :: jAsy, iz, ScatOd
LOGICAL :: lJout, lHex
! ----------------------------------------------------
TYPE (Cell_Type), POINTER, DIMENSION(:) :: Cell
TYPE (Pin_Type),  POINTER, DIMENSION(:) :: Pin

TYPE (DcmpAsyRayInfo_Type), POINTER, DIMENSION(:,:) :: DcmpAsyRay

REAL, POINTER, DIMENSION(:)     :: src1g
REAL, POINTER, DIMENSION(:,:,:) :: SrcAng1g1, SrcAng1g2

INTEGER, POINTER, DIMENSION(:) :: DcmpAsyRayCount

INTEGER :: PinSt, PinEd, FsrSt, FsrEd, kRot, iAsyRay, ifsr, ixy, ibd
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
CALL dmalloc0(TrackingLoc%phis, FsrSt, FsrEd)
IF (ljout) CALL dmalloc0(TrackingLoc%Jout, 1, 3, 1, nbd, PinSt, PinEd)

CALL dmalloc0(TrackingLoc%phim, 1, nOd, FsrSt, FsrEd)

src1g => TrackingLoc%src
! ----------------------------------------------------
IF (lHex .AND. hLgc%l360) THEN
  CALL dmalloc0(TrackingLoc%SrcAng1, 1, nPolarAng, FsrSt, FsrEd, 1, 1)
  CALL dmalloc0(TrackingLoc%SrcAng2, 1, nPolarAng, FsrSt, FsrEd, 1, 1)
  
  DO iazi = 1, nAziAng
    SrcAng1g1 => TrackingLoc%SrcAng1
    SrcAng1g2 => TrackingLoc%SrcAng2
    
    DO ifsr = FsrSt, FsrEd
      DO ipol = 1, nPolarAng
        SrcAng1g1(ipol, ifsr, 1) = src1g(ifsr)
        SrcAng1g2(ipol, ifsr, 1) = src1g(ifsr)
        
        srctmp = comp(1, ipol, iazi) * srcm1g(1, ifsr) + comp(2, ipol, iazi) * srcm1g(2, ifsr)
        
        SrcAng1g1(ipol, ifsr, 1) = SrcAng1g1(ipol, ifsr, 1) + srctmp
        SrcAng1g2(ipol, ifsr, 1) = SrcAng1g2(ipol, ifsr, 1) - srctmp
        
        IF (ScatOd .LT. 2) CYCLE
        
        srctmp = comp(3, ipol, iazi) * srcm1g(3, ifsr) + comp(4, ipol, iazi) * srcm1g(4, ifsr) + comp(5, ipol, iazi) * srcm1g(5, ifsr)
        
        SrcAng1g1(ipol, ifsr, 1) = SrcAng1g1(ipol, ifsr, 1) + srctmp
        SrcAng1g2(ipol, ifsr, 1) = SrcAng1g2(ipol, ifsr, 1) + srctmp
        
        IF (ScatOd .LT. 3) CYCLE
        
        srctmp = comp(6, ipol, iazi) * srcm1g(6, ifsr) + comp(7, ipol, iazi) * srcm1g(7, ifsr) + comp(8, ipol, iazi) * srcm1g(8, ifsr) + comp(9, ipol, iazi) * srcm1g(9, ifsr)
        
        SrcAng1g1(ipol, ifsr, 1) = SrcAng1g1(ipol, ifsr, 1) + srctmp
        SrcAng1g2(ipol, ifsr, 1) = SrcAng1g2(ipol, ifsr, 1) - srctmp
      END DO
    END DO
    
    DO krot = 1, 2
      DO iAsyRay = 1, DcmpAziRay(0, iazi, jAsy)
        jAsyRay = DcmpAziRay(iAsyRay, iazi, jAsy)
        
        IF (lHex) THEN
          CALL HexTrackRotRayDcmpP1_GM(RayInfo, CoreInfo, TrackingLoc, DcmpAsyRay(jAsyRay, jAsy), ljout, iz, krot, ScatOd)
        ELSE
          !CALL RecTrackRotRayDcmpP1_GM(RayInfo, CoreInfo, TrackingLoc, DcmpAsyRay(jAsyRay, jAsy), ljout, iz, krot, ScatOd)
        END IF
      END DO
    END DO
  END DO
ELSE
  CALL dmalloc0(TrackingLoc%SrcAng1, 1, nPolarAng, FsrSt, FsrEd, 1, nAziAng)
  CALL dmalloc0(TrackingLoc%SrcAng2, 1, nPolarAng, FsrSt, FsrEd, 1, nAziAng)
  
  SrcAng1g1 => TrackingLoc%SrcAng1
  SrcAng1g2 => TrackingLoc%SrcAng2
  
  DO iazi = 1, nAziAng
    DO ifsr = FsrSt, FsrEd
      DO ipol = 1, nPolarAng
        SrcAng1g1(ipol, ifsr, iazi) = src1g(ifsr)
        SrcAng1g2(ipol, ifsr, iazi) = src1g(ifsr)
        
        srctmp = comp(1, ipol, iazi) * srcm1g(1, ifsr) + comp(2, ipol, iazi) * srcm1g(2, ifsr)
        
        SrcAng1g1(ipol, ifsr, iazi) = SrcAng1g1(ipol, ifsr, iazi) + srctmp
        SrcAng1g2(ipol, ifsr, iazi) = SrcAng1g2(ipol, ifsr, iazi) - srctmp
        
        IF (ScatOd .LT. 2) CYCLE
        
        srctmp = comp(3, ipol, iazi) * srcm1g(3, ifsr) + comp(4, ipol, iazi) * srcm1g(4, ifsr) + comp(5, ipol, iazi) * srcm1g(5, ifsr)
        
        SrcAng1g1(ipol, ifsr, iazi) = SrcAng1g1(ipol, ifsr, iazi) + srctmp
        SrcAng1g2(ipol, ifsr, iazi) = SrcAng1g2(ipol, ifsr, iazi) + srctmp
        
        IF (ScatOd .LT. 3) CYCLE
        
        srctmp = comp(6, ipol, iazi) * srcm1g(6, ifsr) + comp(7, ipol, iazi) * srcm1g(7, ifsr) + comp(8, ipol, iazi) * srcm1g(8, ifsr) + comp(9, ipol, iazi) * srcm1g(9, ifsr)
        
        SrcAng1g1(ipol, ifsr, iazi) = SrcAng1g1(ipol, ifsr, iazi) + srctmp
        SrcAng1g2(ipol, ifsr, iazi) = SrcAng1g2(ipol, ifsr, iazi) - srctmp
      END DO
    END DO
  END DO
  
  DO krot = 1, 2
    DO iAsyRay = 1, DcmpAsyRayCount(jAsy)
      IF (lHex) THEN
        CALL HexTrackRotRayDcmpP1_GM(RayInfo, CoreInfo, TrackingLoc, DcmpAsyRay(iAsyRay, jAsy), ljout, iz, krot, ScatOd)
      ELSE
        !CALL RecTrackRotRayDcmpP1_GM(RayInfo, CoreInfo, TrackingLoc, DcmpAsyRay(iAsyRay, jAsy), ljout, iz, krot, ScatOd)
      END IF
    END DO
  END DO
END IF
! ----------------------------------------------------
DO ifsr = FsrSt, FsrEd
  phis1g(ifsr) = phis1g(ifsr) + TrackingLoc%phis(ifsr)
  
  DO iod = 1, nod
    phim1g(iod, ifsr) = phim1g(iod, ifsr) + TrackingLoc%phim(iod, ifsr)
  END DO
END DO

IF (lJout) THEN
  DO ixy = PinSt, PinEd
    DO ibd = 1, nbd
      Mocjout1g(:, ibd, ixy) = Mocjout1g(:, ibd, ixy) + TrackingLoc%jout(:, ibd, ixy)
    END DO
  END DO
END IF
! ----------------------------------------------------
DEALLOCATE (TrackingLoc%phis)
IF (lJout) DEALLOCATE (TrackingLoc%Jout)

DEALLOCATE (TrackingLoc%phim)
DEALLOCATE (TrackingLoc%SrcAng1)
DEALLOCATE (TrackingLoc%SrcAng2)

NULLIFY (Cell)
NULLIFY (Pin)
NULLIFY (DcmpAsyRay)
NULLIFY (DcmpAsyRayCount)

NULLIFY (src1g)
NULLIFY (SrcAng1g1)
NULLIFY (SrcAng1g2)
! ----------------------------------------------------

END SUBROUTINE RtDcmpP1Thr_GM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayDcmpP1_GM(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, ljout, iz, krot, ScatOd)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type, DcmpAsyRayInfo_Type, Pin_Type, AziAngleInfo_Type
USE HexType, ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexCoreRay, Type_HexRotRay
USE HexData, ONLY : hAsy, haRay, hAsyTypInfo, hLgc

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: iz, krot, ScatOd
! ----------------------------------------------------
TYPE (Pin_Type),          POINTER, DIMENSION(:) :: Pin
TYPE (AziAngleInfo_Type), POINTER, DIMENSION(:) :: AziAng
  
TYPE (Type_HexAsyRay), POINTER :: haRay_Loc
TYPE (Type_HexCelRay), POINTER :: CelRay_Loc

REAL, POINTER, DIMENSION(:)       :: phis1g, src1g, xst1g
REAL, POINTER, DIMENSION(:,:)     :: PhiAngIn1g, EXPA, EXPB, wtang, hwt, phim1g
REAL, POINTER, DIMENSION(:,:,:)   :: jout1g, LocMwt, LocSrc
REAL, POINTER, DIMENSION(:,:,:,:) :: DcmpPhiAngIn1g, DcmpPhiAngOut1g

INTEGER, POINTER, DIMENSION(:) :: AsyRayList, DirList, AziList

INTEGER :: mp(2)
INTEGER :: iazi, jazi, ipol, irotray, iasyray, iray, irayseg, idir, icel, iasy, ifsr, isurf, jbeg, jend, jinc, imray, ipray
INTEGER :: ipst, iped, ipinc, isfst, isfed, isgst, isged, isginc, iAsyTyp, iGeoTyp, icBss, jhPin, jcBss, iOd, nOd
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
phis1g     => TrackingDat%phis
src1g      => TrackingDat%src
xst1g      => TrackingDat%xst
PhiAngIn1g => TrackingDat%PhiAngIn
jout1g     => TrackingDat%jout
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

! P1
phim1g => TrackingDat%phim

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
  IF (hLgc%l360)   jazi = 1
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
    LocSrc => TrackingDat%SrcAng1
  ELSE
    LocSrc => TrackingDat%SrcAng2
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
        
        phid = (PhiAngOut1g(ipol) - LocSrc(ipol, ifsr, jazi)) * ExpApp
        
        PhiAngOut1g(ipol) = PhiAngOut1g(ipol) - phid
        
        phis1g(ifsr) = phis1g(ifsr) + wtazi(ipol) * phid
        
        DO iod = 1, nod
          phim1g(iod, ifsr) = phim1g(iod, ifsr) + LocMwt(iod, ipol, iazi) * phid ! NOTICE
        END DO
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

! Dcmp.
NULLIFY (AsyRayList)
NULLIFY (DirList)
NULLIFY (AziList)
NULLIFY (DcmpPhiAngIn1g)
NULLIFY (DcmpPhiAngOut1g)

! P1
NULLIFY (phim1g)
NULLIFY (LocMwt)
NULLIFY (LocSrc)
! ----------------------------------------------------

END SUBROUTINE HexTrackRotRayDcmpP1_GM
! ------------------------------------------------------------------------------------------------------------