#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmp_GM(RayInfo, CoreInfo, phis, PhiAngIn, xst, src, MocJout, iz, lJout)

USE OMP_LIB
USE allocs
USE PARAM,       ONLY : ZERO, RED, BLACK, GREEN
USE TYPEDEF,     ONLY : RayInfo_Type, Coreinfo_type, Asy_Type, AsyInfo_Type, Pin_Type, Cell_Type, DcmpAsyRayInfo_Type
USE Moc_Mod,     ONLY : TrackingDat, DcmpPhiAngIn1g, DcmpPhiAngOut1g, DcmpColorAsy, RayTraceDcmp_OMP, DcmpGatherBndyFlux1g, DcmpScatterBndyFlux1g, DcmpLinkBndyFlux1g
USE PE_MOD,      ONLY : PE
USE CNTL,        ONLY : nTracerCntl
USE geom,        ONLY : nbd
USE itrcntl_mod, ONLY : itrcntl
USE HexData,     ONLY : hAsy, hLgc

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:)     :: phis, xst, src
REAL, POINTER, DIMENSION(:,:)   :: PhiAngIn
REAL, POINTER, DIMENSION(:,:,:) :: MocJout

INTEGER :: iz
LOGICAL :: lJout
! ----------------------------------------------------
TYPE (AsyInfo_Type), POINTER, DIMENSION(:) :: AsyInfo
TYPE (Asy_Type),     POINTER, DIMENSION(:) :: Asy
TYPE (Cell_Type),    POINTER, DIMENSION(:) :: Cell
TYPE (Pin_Type),     POINTER, DIMENSION(:) :: Pin

TYPE (DcmpAsyRayInfo_Type), POINTER, DIMENSION(:,:) :: DcmpAsyRay

INTEGER, POINTER, DIMENSION(:) :: DcmpAsyRayCount

INTEGER :: ithr, nThr, iAsy, jAsy, ifsr, ibd, ixy, nxy, FsrIdxSt, icel, jfsr, iAsyRay, krot, icolor, jcolor, ncolor, iit, PinSt, PinEd, FsrSt, FsrEd
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
  TrackingDat(ithr)%src => src
  TrackingDat(ithr)%xst => xst
END DO

DcmpPhiAngOut1g = ZERO
phis = ZERO
IF (ljout) Mocjout = ZERO
! ----------------------------------------------------
DO icolor = 1, ncolor
  IF (lHex) THEN
    jcolor = AuxHex(icolor, iit)
    
    IF (hLgc%l060) jcolor = icolor
  ELSE
    jcolor = AuxRec(icolor, iit)
  END IF
  
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) CALL DcmpScatterBndyFlux1g(RayInfo, PhiAngIn, DcmpPhiAngIn1g)
#endif
  
  DO ithr = 1, nthr
    TrackingDat(ithr)%PhiAngIn        => PhiAngIn
    TrackingDat(ithr)%DcmpPhiAngIn1g  => DcmpPhiAngIn1g
    TrackingDat(ithr)%DcmpPhiAngOut1g => DcmpPhiAngOut1g
  END DO
  
  !$OMP PARALLEL PRIVATE(ithr, iAsy, jAsy, PinSt, PinEd, FsrSt, FsrEd, krot, iAsyRay, ifsr, ixy, ibd)
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
    CALL dmalloc0(TrackingDat(ithr)%phis, FsrSt, FsrEd)
    IF (ljout) CALL dmalloc0(TrackingDat(ithr)%Jout, 1, 3, 1, nbd, PinSt, PinEd)
    
    ! RT
    DO krot = 1, 2
      DO iAsyRay = 1, DcmpAsyRayCount(jAsy)
        IF (lHex) THEN
          CALL HexTrackRotRayDcmp_GM(RayInfo, CoreInfo, TrackingDat(ithr), DcmpAsyRay(iAsyRay, jAsy), ljout, iz, krot)
        ELSE
          !CALL RecTrackRotRayDcmp_GM(RayInfo, CoreInfo, TrackingDat(ithr), DcmpAsyRay(iAsyRay, jAsy), ljout, iz,  krot)
        END IF
      END DO
    END DO
    
    ! GATHER
    DO ifsr = FsrSt, FsrEd
      phis(ifsr) = phis(ifsr) + TrackingDat(ithr)%phis(ifsr)
    END DO
    
    IF (ljout) THEN
      DO ixy = PinSt, PinEd
        DO ibd = 1, nbd
          Mocjout(:, ibd, ixy) = Mocjout(:, ibd, ixy) + TrackingDat(ithr)%jout(:, ibd, ixy)
        END DO
      END DO
    END IF
    
    ! FREE
    DEALLOCATE (TrackingDat(ithr)%phis)
    IF (ljout) DEALLOCATE (TrackingDat(ithr)%Jout)
  END DO
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL
  
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) CALL DcmpGatherBndyFlux1g(RayInfo, DcmpPhiAngOut1g)
#endif
  
  IF (PE%RTMASTER) CALL DcmpLinkBndyFlux1g(CoreInfo, RayInfo, PhiAngIn, DcmpPhiAngIn1g, DcmpPhiAngOut1g, jcolor)
END DO
! ----------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ixy, FsrIdxSt, icel, ifsr, jfsr)
!$OMP DO SCHEDULE(GUIDED)
DO ixy = 1, nxy
  FsrIdxSt = Pin(ixy)%FsrIdxSt
  icel     = Pin(ixy)%Cell(iz)
  
  DO ifsr = 1, Cell(icel)%nFsr
    jfsr = FsrIdxSt + ifsr - 1
    
    phis(jfsr) = phis(jfsr) / (xst(jfsr) * Cell(icel)%vol(ifsr)) + src(jfsr)
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

END SUBROUTINE RayTraceDcmp_GM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayDcmp_GM(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, ljout, iz, krot)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, TrackingDat_Type, DcmpAsyRayInfo_Type, AziAngleInfo_Type
USE HexType, ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexCoreRay, Type_HexRotRay
USE HexData, ONLY : hAsy, haRay, hAsyTypInfo

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: iz, krot
! ----------------------------------------------------
TYPE (Pin_Type),          POINTER, DIMENSION(:) :: Pin
TYPE (AziAngleInfo_Type), POINTER, DIMENSION(:) :: AziAng
  
TYPE (Type_HexAsyRay), POINTER :: haRay_Loc
TYPE (Type_HexCelRay), POINTER :: CelRay_Loc

REAL, POINTER, DIMENSION(:)       :: phis1g, src1g, xst1g
REAL, POINTER, DIMENSION(:,:)     :: PhiAngIn1g, EXPA, EXPB, wtang, hwt
REAL, POINTER, DIMENSION(:,:,:)   :: jout1g
REAL, POINTER, DIMENSION(:,:,:,:) :: DcmpPhiAngIn1g, DcmpPhiAngOut1g

INTEGER, POINTER, DIMENSION(:) :: AsyRayList, DirList, AziList

INTEGER :: mp(2)
INTEGER :: iazi, ipol, irotray, iasyray, iray, irayseg, idir, icel, iasy, iFSR, isurf, jbeg, jend, jinc, imray, ipray
INTEGER :: ipst, iped, ipinc, isfst, isfed, isgst, isged, isginc, iAsyTyp, iGeoTyp, icBss, jhPin, jcBss
INTEGER :: nAsyRay, nPinRay, nRaySeg, FsrIdxSt, nPolarAng, PhiAnginSvIdx, ExpAppIdx

REAL :: wtazi(10), locwt(10), loccs, locsn

REAL, DIMENSION(RayInfo%nPolarAngle) :: PhiAngOut1g
REAL :: phid, tau, ExpApp, wtsurf

DATA mp /2, 1/
! ----------------------------------------------------

! Ray Info.
nPolarAng = RayInfo%nPolarAngle
AziAng   => RayInfo%AziAngle

! Geo.
Pin => CoreInfo%Pin

! Tracking Dat
phis1g          => TrackingDat%phis
src1g           => TrackingDat%src
xst1g           => TrackingDat%xst
PhiAngIn1g      => TrackingDat%PhiAngIn
jout1g          => TrackingDat%jout
DcmpPhiAngIn1g  => TrackingDat%DcmpPhiAngIn1g
DcmpPhiAngOut1g => TrackingDat%DcmpPhiAngOut1g
EXPA            => TrackingDat%EXPA
EXPB            => TrackingDat%EXPB
wtang           => TrackingDat%wtang
hwt             => TrackingDat%hwt

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
    ipst = 1; iped = nPinRay; ipinc = 1;  isfst = 1; isfed = 2
  ELSE
    iped = 1; ipst = nPinRay; ipinc = -1; isfed = 1; isfst = 2
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
      iFSR = CelRay_Loc%MshIdx(iRaySeg) + Pin(jhPin)%FsrIdxSt - 1
      
      tau = -CelRay_Loc%SegLgh(iRaySeg) * xst1g(iFSR) ! Optimum Length
      
      ExpAppIdx = max(INT(tau), -40000)
      ExpAppIdx = min(0, ExpAppIdx)
      
      DO ipol = 1, nPolarAng
        ExpApp = ExpA(ExpAppIdx, ipol) * tau + ExpB(ExpAppIdx, ipol)
        
        phid = (PhiAngOut1g(ipol) - src1g(iFSR)) * ExpApp
        
        PhiAngOut1g(ipol) = PhiAngOut1g(ipol) - phid
        
        phis1g(iFSR) = phis1g(iFSR) + wtazi(ipol) * phid
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
NULLIFY (phis1g)
NULLIFY (src1g)
NULLIFY (xst1g)
NULLIFY (PhiAngIn1g)
NULLIFY (jout1g)

! Dcmp.
NULLIFY (AsyRayList)
NULLIFY (DirList)
NULLIFY (AziList)
NULLIFY (DcmpPhiAngIn1g)
NULLIFY (DcmpPhiAngOut1g)
! ----------------------------------------------------

END SUBROUTINE HexTrackRotRayDcmp_GM
! ------------------------------------------------------------------------------------------------------------