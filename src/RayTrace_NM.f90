#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceNM_OMP(RayInfo, CoreInfo, phisnm, PhiAngInnm, xstnm, srcnm, joutnm, iz, gb, ge, ljout, lDomainDcmp, FastMocLv)

USE PARAM
USE TIMER
USE ALLOCS
USE OMP_LIB
USE TYPEDEF,     ONLY : RayInfo_Type, CoreInfo_type, Pin_Type, Cell_Type, MultigridInfo_Type
USE ITRCNTL_MOD, ONLY : ItrCntl
USE Moc_Mod,     ONLY : TrackingDat
USE geom,        ONLY : ng, nbd
USE PE_MOD,      ONLY : PE
USE CNTL,        ONLY : nTracerCntl

IMPLICIT NONE

TYPE(RayInfo_Type)  :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo

INTEGER :: iz, gb, ge, FastMocLv
LOGICAL :: ljout, lDomainDcmp

REAL, POINTER, DIMENSION(:,:)     :: phisnm, xstnm, srcnm
REAL, POINTER, DIMENSION(:,:,:)   :: PhiAngInnm
REAL, POINTER, DIMENSION(:,:,:,:) :: joutnm

LOGICAL, SAVE :: lfirst
DATA lfirst /TRUE/

TYPE(Cell_Type),          POINTER, DIMENSION(:) :: Cell
TYPE(Pin_Type),           POINTER, DIMENSION(:) :: Pin
TYPE(MultigridInfo_Type), POINTER, DIMENSION(:) :: MultigridInfo

INTEGER :: nAziAng, nPolarAng, nPhiAngSv, nRotray, nFsr, nAsy, nxy, nThread
INTEGER :: iRotRay, ipol, iazi, krot, AziIdx, ithr, FsrIdxSt, ipin, icel, ibd, ilv, iray, ifsr, jfsr, ig
! ----------------------------------------------------

nAsy = CoreInfo%nxya
nFsr = CoreInfo%nCoreFsr
nxy  = CoreInfo%nxy

MultigridInfo => RayInfo%MultigridInfo
nPhiAngSv      = RayInfo%nPhiAngSv
nRotRay        = RayInfo%nRotRay

nThread = PE%nThread

ilv = nTracerCntl%gridStr(ItrCntl%mocit)
IF (ItrCntl%mocit .GT. nTracerCntl%gridNum) ilv = nTracerCntl%MultigridLV

nAziAng   = MultigridInfo(ilv)%nAzi
nPolarAng = MultigridInfo(ilv)%nPolar

CALL omp_set_num_threads(nThread)
! ----------------------------------------------------
IF (lfirst) THEN
  lFirst = FALSE

  DO ithr = 1, nThread
    IF (TrackingDat(ithr)%lAllocNM) CYCLE
    
    CALL Dmalloc(TrackingDat(ithr)%phisnm, ng, nFsr)
    
    IF (.NOT. lDomainDcmp) CALL Dmalloc(TrackingDat(ithr)%Joutnm, 3, ng, nbd, nxy)
    
    TrackingDat(ithr)%lAllocNM = TRUE
  END DO
END IF
! ----------------------------------------------------
!$OMP PARALLEL PRIVATE(ithr, krot, iazi, iRotRay, iray, AziIdx)
ithr = omp_get_thread_num() + 1

TrackingDat(ithr)%phisnm = ZERO

TrackingDat(ithr)%EXPA       => MultigridInfo(ilv)%EXPA
TrackingDat(ithr)%EXPB       => MultigridInfo(ilv)%EXPB
TrackingDat(ithr)%srcnm      => srcnm
TrackingDat(ithr)%xstnm      => xstnm
TrackingDat(ithr)%PhiAngInnm => PhiAngInnm

IF (ljout) TrackingDat(ithr)%joutnm = ZERO

DO krot = 1, 2
  DO iazi = 1, nAziAng / 2
    AziIdx = MultigridInfo(ilv)%AziList(iazi)
    
    IF (nTracerCntl%lHex) THEN
      !$OMP DO SCHEDULE(GUIDED)
      DO iray = 1, RayInfo%RotRayAziList(0, Aziidx)
        iRotRay = RayInfo%RotRayAziList(iray, AziIdx)
        
        CALL HexTrackRotRayNM_OMP(RayInfo, CoreInfo, TrackingDat(ithr), ljout, iRotRay, iz, ilv, krot, gb, ge)
      END DO
      !$OMP END DO NOWAIT
    ELSE
      !$OMP DO SCHEDULE(GUIDED)
      DO iray = 1, RayInfo%RotRayAziList(0, Aziidx)
        iRotRay = RayInfo%RotRayAziList(iray, AziIdx)
        
        CALL RecTrackRotRayNM_OMP(RayInfo, CoreInfo, TrackingDat(ithr), ljout, iRotRay, iz, ilv, krot, gb, ge)
      END DO
      !$OMP END DO NOWAIT
    END IF
  END DO
END DO
!$OMP END PARALLEL
! ----------------------------------------------------
phisnm(gb:ge, :) = ZERO
IF(ljout) joutnm(:, gb:ge, :, :) = ZERO

DO ithr = 1, nThread
  phisnm(gb:ge, 1:nfsr) = phisnm(gb:ge, 1:nfsr) + TrackingDat(ithr)%phisnm(gb:ge, 1:nfsr)
  
  IF (.NOT. ljout) CYCLE
    
  joutnm(:, gb:ge, 1:nbd, 1:nxy) = joutnm(:, gb:ge, 1:nbd, 1:nxy) + TrackingDat(ithr)%joutnm(:, gb:ge, 1:nbd, 1:nxy)
END DO
! ----------------------------------------------------
Cell => CoreInfo%CellInfo
Pin  => CoreInfo%Pin

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ipin, FsrIdxSt, icel, ifsr, jfsr, ig)
!$OMP DO SCHEDULE(GUIDED)
DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt
  icel     = Pin(ipin)%Cell(iz)
  
  DO ifsr = 1, Cell(icel)%nFsr
    jfsr = FsrIdxSt + ifsr - 1
    
    DO ig = gb, ge
      phisnm(ig, jfsr) = phisnm(ig, jfsr) / (xstnm(ig, jfsr) * Cell(icel)%vol(ifsr)) + srcnm(ig, jfsr)
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
! ----------------------------------------------------
NULLIFY (Cell)
NULLIFY (Pin)
NULLIFY (MultigridInfo)
! ----------------------------------------------------

END SUBROUTINE RayTraceNM_OMP
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayNM_OMP(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, ilv, krot, gb, ge)

USE PARAM
USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, Asy_Type, PinInfo_Type, Cell_Type, AsyRayInfo_type, CoreRayInfo_Type, RotRayInfo_Type, CellRayInfo_type, TrackingDat_Type

IMPLICIT NONE

TYPE(RayInfo_Type)     :: RayInfo
TYPE(CoreInfo_Type)    :: CoreInfo
TYPE(TrackingDat_Type) :: TrackingDat

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz, ilv, krot, gb, ge

TYPE(Pin_Type),         POINTER, DIMENSION(:) :: Pin
TYPE(Asy_Type),         POINTER, DIMENSION(:) :: Asy
TYPE(PinInfo_Type),     POINTER, DIMENSION(:) :: PinInfo
TYPE(Cell_Type),        POINTER, DIMENSION(:) :: Cell
TYPE(AsyRayInfo_type),  POINTER, DIMENSION(:) :: AsyRay
TYPE(CoreRayInfo_Type), POINTER, DIMENSION(:) :: CoreRay
TYPE(RotRayInfo_Type),  POINTER, DIMENSION(:) :: RotRay

TYPE(CellRayInfo_Type), POINTER :: CellRay

INTEGER :: iazi, ipol, icray, jcray, iaray, jaray, ipray, iceray, irayseg, idir, ipin, icel, iasy, ifsr, isurf, jbeg, jend, jinc, ig, ibcel
INTEGER :: nRotRay, nCoreRay, nAsyRay, nPinRay, nRaySeg, nPolarAng, FsrIdxSt, PhiAnginSvIdx, PhiAngOutSvIdx, ExpAppIdx

INTEGER :: mp(2)

INTEGER, POINTER, DIMENSION(:) :: LocalFsrIdx

REAL :: PhiAngOut(RayInfo%nPolarAngle, gb:ge)
REAL :: phid, tau, locsrc, ExpApp
REAL :: wt(10), wt2(10, 4)

REAL, POINTER, DIMENSION(:)       :: LenSeg
REAL, POINTER, DIMENSION(:,:)     :: phis, src, xst, expa, expb, wtang, wtazi
REAL, POINTER, DIMENSION(:,:,:)   :: PhiAngIn, wtsurf
REAL, POINTER, DIMENSION(:,:,:,:) :: jout

DATA mp /2, 1/
! ----------------------------------------------------

AsyRay  => RayInfo%AsyRay
CoreRay => RayInfo%CoreRay
RotRay  => RayInfo%RotRay

wtang    => RayInfo%MultigridInfo(ilv)%wtang
wtsurf   => RayInfo%MultigridInfo(ilv)%wtsurf
nPolarAng = RayInfo%MultigridInfo(ilv)%nPolar

Asy     => CoreInfo%Asy
Pin     => CoreInfo%Pin
PinInfo => CoreInfo%Pininfo
Cell    => CoreInfo%CellInfo

phis     => TrackingDat%phisnm
src      => TrackingDat%srcnm
xst      => TrackingDat%xstnm
jout     => TrackingDat%joutnm
PhiAngIn => TrackingDat%phiAngInnm
EXPA     => TrackingDat%EXPA
EXPB     => TrackingDat%EXPB

nCoreRay = RotRay(irotRay)%nRay

PhiAngInSvIdx  = RayInfo%PhiAngInSvIdx (iRotRay, krot)
PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay, krot)

PhiAngOut = PhiAngIn(:, gb:ge, PhiAnginSvIdx)

IF (krot .EQ. 1) THEN
  jbeg = 1; jend = nCoreRay; jinc = 1
ELSE
  jbeg = nCoreRay; jend = 1; jinc = -1
END IF
! ----------------------------------------------------
DO icray = jbeg, jend, jinc
  jcray   = RotRay(iRotRay)%RayIdx(icray)
  nAsyRay = CoreRay(jcray)%nRay
  iazi    = CoreRay(jcray)%iang
  
  idir = RotRay(iRotRay)%dir(icray)
  IF (krot .eq. 2) idir = mp(idir)
  
  DO ipol = 1, nPolarAng
    wt(ipol) = wtang(ipol, iazi)
    
    IF (lJout) wt2(ipol, 1:4) = wtsurf(ipol, iazi, 1:4)
  END DO
  ! --------------------------------------------------
  IF (idir .EQ. 1) THEN
    DO iaray = 1, nAsyRay
      jaray = CoreRay(jcray)%AsyRayIdx(iaray)
      iasy  = CoreRay(jcray)%AsyIdx(iaray)
      
      IF (iAsy .EQ. 0) CYCLE
      
      nPinRay = AsyRay(jaray)%nCellRay
      
      DO ipray = 1, nPinRay
        ipin   = AsyRay(jaray)%PinIdx(ipray)
        iceray = AsyRay(jaray)%PinRayIdx(ipray)
        ipin   = Asy(iAsy)%GlobalPinIdx(ipin)
        
        icel     = Pin(ipin)%Cell(iz)
        FsrIdxSt = Pin(ipin)%FsrIdxSt
        
        ibcel    = Cell(icel)%basecellstr
        CellRay => Cell(ibcel)%CellRay(iceray)
        
        nRaySeg      = CellRay%nSeg
        LocalFsrIdx => CellRay%LocalFsrIdx
        LenSeg      => CellRay%LenSeg
        
        IF (lJout) THEN
          isurf = AsyRay(jaray)%PinRaySurf(1, ipray)
          
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(1, ig, isurf, ipin) = Jout(1, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
              Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
            END DO
          END DO
        END IF
        
        DO iRaySeg = 1, nRaySeg
          ifsr = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
          
          DO ig = gb, ge
            tau = -LenSeg(iRaySeg) * xst(ig, ifsr)
            
            ExpAppIdx = max(INT(tau), -40000)
            ExpAppIdx = min(0, ExpAppIdx)
            
            locsrc = src(ig, ifsr)
            
            DO ipol = 1, nPolarAng
              ExpApp = EXPA(ipol, ExpAppIdx) * tau + EXPB(ipol, ExpAppIdx)
              
              phid = (PhiAngOut(ipol, ig) - locsrc) * ExpApp
              
              PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
              
              phis(ig, ifsr) = phis(ig, ifsr) + wt(ipol) * phid
            END DO
          END DO
        END DO
        
        IF (lJout) THEN
          isurf = AsyRay(jaray)%PinRaySurf(2, ipray)
          
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(2, ig, isurf, ipin) = Jout(2, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
              Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
            END DO
          END DO
        END IF
      END DO
    END DO
  ! --------------------------------------------------
  ELSE
    DO iaray = nAsyRay, 1, -1
      jaray = CoreRay(jcray)%AsyRayIdx(iaray)
      iasy  = CoreRay(jcray)%AsyIdx(iaray)
      
      IF (iAsy .EQ. 0) CYCLE
      
      nPinRay = AsyRay(jaray)%nCellRay
      
      DO ipray = nPinRay, 1, -1
        ipin   = AsyRay(jaray)%PinIdx(ipray)
        iceray = AsyRay(jaray)%PinRayIdx(ipray)
        ipin   = Asy(iAsy)%GlobalPinIdx(ipin)
        
        icel     = Pin(ipin)%Cell(iz)
        FsrIdxSt = Pin(ipin)%FsrIdxSt
        
        ibcel    = Cell(icel)%basecellstr
        CellRay => Cell(ibcel)%CellRay(iceray)
        
        nRaySeg      = CellRay%nSeg
        LocalFsrIdx => CellRay%LocalFsrIdx
        LenSeg      => CellRay%LenSeg
        
        IF (lJout) THEN
          isurf = AsyRay(jaray)%PinRaySurf(2, ipray)
          
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(1, ig, isurf, ipin) = Jout(1, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
              Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
            END DO
          END DO
        END IF
        
        DO iRaySeg = nRaySeg, 1, -1
          ifsr = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
          
          DO ig = gb, ge
            tau = -LenSeg(iRaySeg) * xst(ig, ifsr)
            
            ExpAppIdx = max(INT(tau), -40000)
            ExpAppIdx = min(0, ExpAppIdx)
            
            locsrc = src(ig, ifsr)
            
            DO ipol = 1, nPolarAng
              ExpApp = EXPA(ipol, ExpAppIdx) * tau + EXPB(ipol, ExpAppIdx)
              
              phid = (PhiAngOut(ipol, ig) - locsrc) * ExpApp
              
              PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
              
              phis(ig, ifsr) = phis(ig, ifsr) + wt(ipol) * phid
            END DO
          END DO
        END DO
        
        IF (lJout) THEN
          isurf = AsyRay(jaray)%PinRaySurf(1, ipray)
          
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(2, ig, isurf, ipin) = Jout(2, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
              Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
            END DO
          END DO
        END IF
      END DO
    END DO
  END IF
END DO

PhiAngIn(:, gb:ge, PhiAngOutSvIdx) = PhiAngOut
! ----------------------------------------------------

END SUBROUTINE RecTrackRotRayNM_OMP
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayNM_OMP(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, ilv, krot, gb, ge)

USE PARAM
USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, TrackingDat_Type, Pin_Type
USE HexData, ONLY : hAsy
USE HexType, ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexCoreRay, Type_HexRotRay
USE HexData, ONLY : haRay, hcRay, hRotRay, hAsyTypInfo

IMPLICIT NONE

TYPE(RayInfo_Type)     :: RayInfo
TYPE(CoreInfo_Type)    :: CoreInfo
TYPE(TrackingDat_Type) :: TrackingDat

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz, ilv, krot, gb, ge
! ----------------------------------------------------

TYPE(Pin_Type), POINTER, DIMENSION(:) :: Pin

TYPE(Type_HexAsyRay), POINTER :: haRay_Loc
TYPE(Type_HexCelRay), POINTER :: CelRay_Loc
TYPE(Type_HexRotRay), POINTER :: hRotRay_Loc

INTEGER :: iazi, ipol, icRay, jcRay, iAsyRay, iRaySeg, imRay, ihpRay, iAsy, ifsr, iSurf, jbeg, jend, jinc, ig, iGeoTyp, iAsyTyp, jhPin, icBss, jcBss
INTEGER :: nCoreRay, nAsyRay, nPolarAng, PhiAnginSvIdx, PhiAngOutSvIdx, ExpAppIdx

REAL :: phid, tau, locsrc, ExpApp

REAL :: wtazi(10)
REAL :: PhiAngOut(RayInfo%nPolarAngle, gb:ge)

REAL, POINTER, DIMENSION(:,:)     :: phis, src, xst, EXPA, EXPB, wtang
REAL, POINTER, DIMENSION(:,:,:,:) :: jout
! ----------------------------------------------------

wtang    => RayInfo%MultigridInfo(ilv)%wtang
!wtsurf  => RayInfo%MultigridInfo(ilv)%wtsurf
nPolarAng = RayInfo%MultigridInfo(ilv)%nPolar

PhiAngInSvIdx  = RayInfo%PhiAngInSvIdx (iRotRay, krot)
PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay, krot)

Pin => CoreInfo%Pin

Phis     => TrackingDat%phisnm
src      => TrackingDat%srcnm
xst      => TrackingDat%xstnm
jout     => TrackingDat%joutnm
PhiAngOut = TrackingDat%phiAngInnm(:, gb:ge, PhiAnginSvIdx)
EXPA     => TrackingDat%EXPA
EXPB     => TrackingDat%EXPB

hRotRay_Loc => hRotRay(iRotRay)
nCoreRay     = hRotRay_Loc%ncRay

IF (krot .EQ. 1) THEN
  jbeg = 1; jend = nCoreRay; jinc = 1
ELSE
  jbeg = nCoreRay; jend = 1; jinc = -1
END IF
! ----------------------------------------------------
DO icRay = jbeg, jend, jinc
  jcRay   = hRotRay_Loc%cRayIdx(icRay)
  nAsyRay = hcRay(abs(jcRay))%nmRay
  iazi    = hcRay(abs(jcRay))%AzmIdx
  
  IF (krot .EQ. 2) jcRay = -jcRay ! Reverse the Sweep Direction
  
  DO ipol = 1, nPolarAng
    wtazi(ipol) = wtang(ipol, iazi)
    
    !IF (lJout) wtazi2(ipol, 1:4) = wtsurf(ipol, iazi, 1:4)
  END DO
  ! --------------------------------------------------
  IF (jcRay .GT. 0) THEN
    DO imRay = 1, nAsyRay
      iAsyRay = hcRay(abs(jcRay))%mRayIdx(imRay)
      iAsy    = hcRay(abs(jcRay))%AsyIdx(imRay)
      
      IF (iAsy .EQ. 0) CYCLE
      
      iAsyTyp = hAsy(iAsy)%AsyTyp
      iGeoTyp = hAsy(iAsy)%GeoTyp
      icBss   = hAsyTypInfo(iAsyTyp)%iBss
      
      haRay_Loc => haRay(iGeoTyp, icBss, iAsyRay)
      
      DO ihpRay = 1, haRay_Loc%nhpRay
        jhPin = haRay_Loc%CelRay(ihpRay)%hPinIdx
        jhPin = hAsy(iAsy)%PinIdxSt + hAsyTypInfo(iAsyTyp)%PinLocIdx(iGeoTyp, jhPin) - 1
        jcBss = Pin(jhPin)%hCelGeo(iz)
        
        CelRay_Loc => haRay(iGeoTyp, jcBss, iAsyRay)%CelRay(ihpRay)
        
        IF (lJout) THEN
          iSurf = CelRay_Loc%SurfIdx(1) ! y : small
          
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(1, ig, iSurf, jhPin) = Jout(1, ig, isurf, jhPin) + wtazi(ipol) * PhiAngOut(ipol, ig)
              !Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + wtazi2(ipol, isurf) * PhiAngOut(ipol, ig)
            END DO
          END DO
        END IF
        
        DO iRaySeg = 1, CelRay_Loc%nSegRay
          ifsr = CelRay_Loc%MshIdx(iRaySeg) + Pin(jhPin)%FsrIdxSt - 1
          
          DO ig = gb, ge
            tau = -CelRay_Loc%SegLgh(iRaySeg) * xst(ig, ifsr) ! Optimum Length
            
            ExpAppIdx = max(INT(tau), -40000)
            ExpAppIdx = min(0, ExpAppIdx)
            
            locsrc = src(ig, ifsr)
            
            DO ipol = 1, nPolarAng
              ExpApp = EXPA(ipol, ExpAppIdx) * tau + EXPB(ipol, ExpAppIdx)
              
              phid = (PhiAngOut(ipol, ig) - locsrc) * ExpApp
              
              PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
              
              phis(ig, ifsr) = phis(ig, ifsr) + wtazi(ipol) * phid
            END DO
          END DO
        END DO
        
        IF (lJout) THEN
          isurf = CelRay_Loc%SurfIdx(2) ! y : Big
          
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(2, ig, iSurf, jhPin) = Jout(2, ig, iSurf, jhPin) + wtazi(ipol) * PhiAngOut(ipol, ig)
              !Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + wtazi2(ipol, isurf) * PhiAngOut(ipol, ig)
            END DO
          END DO
        END IF
      END DO
    END DO
  ! --------------------------------------------------
  ELSE
    DO imRay = nAsyRay, 1, -1
      iAsyRay = hcRay(abs(jcRay))%mRayIdx(imRay)
      iAsy    = hcRay(abs(jcRay))%AsyIdx(imRay)
      
      IF (iAsy .EQ. 0) CYCLE
      
      iAsyTyp = hAsy(iAsy)%AsyTyp
      iGeoTyp = hAsy(iAsy)%GeoTyp
      icBss   = hAsyTypInfo(iAsyTyp)%iBss
      
      haRay_Loc => haRay(iGeoTyp, icBss, iAsyRay)
      
      DO ihpRay = haRay_Loc%nhpRay, 1, -1
        jhPin = haRay_Loc%CelRay(ihpRay)%hPinIdx
        jhPin = hAsy(iAsy)%PinIdxSt + hAsyTypInfo(iAsyTyp)%PinLocIdx(iGeoTyp, jhPin) - 1
        jcBss = Pin(jhPin)%hCelGeo(iz)
        
        CelRay_Loc => haRay(iGeoTyp, jcBss, iAsyRay)%CelRay(ihpRay)
        
        IF (lJout) THEN
          iSurf = CelRay_Loc%SurfIdx(2) ! y : Big
          
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(1, ig, iSurf, jhPin) = Jout(1, ig, isurf, jhPin) + wtazi(ipol) * PhiAngOut(ipol, ig)
              !Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + wtazi2(ipol, isurf) * PhiAngOut(ipol, ig)
            END DO
          END DO
        END IF
        
        DO iRaySeg = CelRay_Loc%nSegRay, 1, -1
          ifsr = CelRay_Loc%MshIdx(iRaySeg) + Pin(jhPin)%FsrIdxSt - 1
          
          DO ig = gb, ge
            tau = -CelRay_Loc%SegLgh(iRaySeg) * XsT(ig, ifsr) ! Optimum Length
            
            ExpAppIdx = max(INT(tau), -40000)
            ExpAppIdx = min(0, ExpAppIdx)
            
            locsrc = src(ig, ifsr)
            
            DO ipol = 1, nPolarAng
              ExpApp = EXPA(ipol, ExpAppIdx) * tau + EXPB(ipol, ExpAppIdx)
              
              phid = (PhiAngOut(ipol, ig) - locsrc) * ExpApp
              
              PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
              
              phis(ig, ifsr) = phis(ig, ifsr) + wtazi(ipol) * phid
            END DO
          END DO
        END DO
        
        IF (lJout) THEN
          isurf = CelRay_Loc%SurfIdx(1) ! y : small
          
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(2, ig, iSurf, jhPin) = Jout(2, ig, iSurf, jhPin) + wtazi(ipol) * PhiAngOut(ipol, ig)
              !Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + wtazi2(ipol, isurf) * PhiAngOut(ipol, ig)
            END DO
          END DO
        END IF
      END DO
    END DO
  END IF
END DO

TrackingDat%PhiAngInNM(:, gb:ge, PhiAngOutSvIdx) = PhiAngOut
! ----------------------------------------------------
NULLIFY (wtang)
NULLIFY (phis)
NULLIFY (src)
NULLIFY (xst)
NULLIFY (jout)
NULLIFY (EXPA)
NULLIFY (EXPB)
NULLIFY (Pin)
NULLIFY (hRotRay_Loc)
NULLIFY (CelRay_Loc)
NULLIFY (hRotRay_Loc)
! ----------------------------------------------------

END SUBROUTINE HexTrackRotRayNM_OMP
! ------------------------------------------------------------------------------------------------------------