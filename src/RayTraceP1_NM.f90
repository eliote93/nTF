#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceP1NM_OMP(RayInfo, CoreInfo, phisnm, phimnm, PhiAngInnm, xstnm, srcnm, srcmnm, joutnm, iz, gb, ge, ljout, lDomainDcmp)

USE PARAM
USE TIMER
USE ALLOCS
USE OMP_LIB
USE TYPEDEF,     ONLY : RayInfo_Type, CoreInfo_type, Pin_Type, Cell_Type, MultigridInfo_Type
USE ITRCNTL_MOD, ONLY : ItrCntl
USE Moc_Mod,     ONLY : TrackingDat, SrcAngnm1, SrcAngnm2
USE geom,        ONLY : ng, nbd
USE PE_MOD,      ONLY : PE
USE CNTL,        ONLY : nTracerCntl

IMPLICIT NONE

TYPE(RayInfo_Type)  :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo

INTEGER :: iz, gb, ge
LOGICAL :: ljout, lDomainDcmp

REAL, POINTER, DIMENSION(:,:)     :: phisnm, xstnm, srcnm
REAL, POINTER, DIMENSION(:,:,:)   :: phimnm, PhiAngInnm, srcmnm
REAL, POINTER, DIMENSION(:,:,:,:) :: joutnm

LOGICAL, SAVE :: lfirst
DATA lfirst /TRUE/

TYPE(Cell_Type),          POINTER, DIMENSION(:) :: Cell
TYPE(Pin_Type),           POINTER, DIMENSION(:) :: Pin
TYPE(MultigridInfo_Type), POINTER, DIMENSION(:) :: MultigridInfo

INTEGER :: nAziAng, nPolarAng, nFsr, nxy, nThread, scatod, nod
INTEGER :: iRotRay, ipol, iazi, irot, AziIdx, ithr, FsrIdxSt, ipin, icel, ibd, ilv, ig, ifsr, jfsr, iray, iod

REAL :: wttmp, wtsin2, wtcos, wtpolar, srctmp

REAL, POINTER, DIMENSION(:,:,:) :: Comp
! ----------------------------------------------------

nFsr    = CoreInfo%nCoreFsr
nxy     = CoreInfo%nxy
nThread = PE%nThread
ScatOd  = nTracerCntl%ScatOd

ilv = nTracerCntl%gridStr(ItrCntl%mocit)
IF (ItrCntl%mocit .GT. nTracerCntl%gridNum) ilv = nTracerCntl%MultigridLV

MultigridInfo => RayInfo%MultigridInfo

Comp     => MultigridInfo(ilv)%Comp
nAziAng   = MultigridInfo(ilv)%nAzi
nPolarAng = MultigridInfo(ilv)%nPolar

nod = 2
IF (ScatOd .EQ. 2) nod = 5
IF (ScatOd .EQ. 3) nod = 9

CALL omp_set_num_threads(nThread)
! ----------------------------------------------------
IF (lfirst) THEN
  lFirst = FALSE
  
  ! Tracking Dat.
  DO ithr = 1, nThread
    IF (TrackingDat(ithr)%lAllocNM) CYCLE
    
    CALL dmalloc(TrackingDat(ithr)%phisnm,      ng, nFsr)
    CALL dmalloc(TrackingDat(ithr)%phimnm, nod, ng, nFsr)
    
    IF (.NOT. lDomainDcmp) CALL dmalloc(TrackingDat(ithr)%Joutnm, 3, ng, nbd, nxy)
    
    TrackingDat(ithr)%lAllocNM = TRUE
  END DO
  
  ! Angular Src.
  CALL dmalloc(SrcAngnm1, ng, nPolarAng, nFsr, nAziAng)
  CALL dmalloc(SrcAngnm2, ng, nPolarAng, nFsr, nAziAng)
END IF
! ----------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ithr, iazi, Ifsr, ipol, ig, srctmp)

ithr = 1

!$ ithr = omp_get_thread_num()+1

!$OMP DO SCHEDULE(GUIDED) COLLAPSE(4)
DO iazi = 1, nAziAng
  DO ifsr = PE%myOmpFsrBeg(ithr), PE%myOmpFsrEnd(ithr)
    DO ipol = 1, nPolarAng
      DO ig = gb, ge
        SrcAngnm1(ig, ipol, ifsr, iazi) = srcnm(ig, ifsr)
        SrcAngnm2(ig, ipol, ifsr, iazi) = srcnm(ig, ifsr)
        
        srctmp = comp(1, ipol, iazi) * srcmnm(1, ig, ifsr) + comp(2, ipol, iazi) * srcmnm(2, ig, ifsr)
        
        SrcAngnm1(ig, ipol, ifsr, iazi) = SrcAngnm1(ig, ipol, ifsr, iazi) + srctmp
        SrcAngnm2(ig, ipol, ifsr, iazi) = SrcAngnm2(ig, ipol, ifsr, iazi) - srctmp
        
        IF (ScatOd .LT. 2) CYCLE
        
        srctmp = comp(3, ipol, iazi) * srcmnm(3, ig, ifsr) + comp(4, ipol, iazi) * srcmnm(4, ig, ifsr) + comp(5, ipol, iazi) * srcmnm(5, ig, ifsr)
        
        SrcAngnm1(ig, ipol, ifsr, iazi) = SrcAngnm1(ig, ipol, ifsr, iazi) + srctmp
        SrcAngnm2(ig, ipol, ifsr, iazi) = SrcAngnm2(ig, ipol, ifsr, iazi) + srctmp
        
        IF (ScatOd .LT. 3) CYCLE
        
        srctmp = comp(6, ipol, iazi) * srcmnm(6, ig, ifsr) + comp(7, ipol, iazi) * srcmnm(7, ig, ifsr) + comp(8, ipol, iazi) * srcmnm(8, ig, ifsr) + comp(9, ipol, iazi) * srcmnm(9, ig, ifsr)
        
        SrcAngnm1(ig, ipol, ifsr, iazi) = SrcAngnm1(ig, ipol, ifsr, iazi) + srctmp
        SrcAngnm2(ig, ipol, ifsr, iazi) = SrcAngnm2(ig, ipol, ifsr, iazi) - srctmp
      END DO
     END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
! ----------------------------------------------------
!$OMP PARALLEL PRIVATE(ithr, iRotRay, AziIdx)
ithr = omp_get_thread_num() + 1

TrackingDat(ithr)%phisnm = ZERO
TrackingDat(ithr)%phimnm = ZERO

TrackingDat(ithr)%EXPA       => MultigridInfo(ilv)%EXPA
TrackingDat(ithr)%EXPB       => MultigridInfo(ilv)%EXPB
TrackingDat(ithr)%srcnm      => srcnm
TrackingDat(ithr)%xstnm      => xstnm
TrackingDat(ithr)%PhiAngInnm => PhiAngInnm
TrackingDat(ithr)%SrcAngnm1  => SrcAngnm1
TrackingDat(ithr)%SrcAngnm2  => SrcAngnm2

IF (ljout) TrackingDat(ithr)%joutnm = ZERO

DO irot = 1, 2
  DO iazi = 1, nAziAng / 2
    AziIdx = MultigridInfo(ilv)%AziList(iazi)
    
    IF (nTracerCntl%lHex) THEN
      !$OMP DO SCHEDULE(GUIDED)
      DO iray = 1, RayInfo%RotRayAziList(0, Aziidx)
        iRotRay = RayInfo%RotRayAziList(iray, AziIdx)
        
        CALL HexTrackRotRayP1NM_OMP(RayInfo, CoreInfo, TrackingDat(ithr), ljout, iRotRay, iz, ilv, irot, gb, ge, scatod)
      END DO
      !$OMP END DO NOWAIT
    ELSE
      !$OMP DO SCHEDULE(GUIDED)
      DO iray = 1, RayInfo%RotRayAziList(0, Aziidx)
        iRotRay = RayInfo%RotRayAziList(iray, AziIdx)
        
        CALL RecTrackRotRayP1NM_OMP(RayInfo, CoreInfo, TrackingDat(ithr), ljout, iRotRay, iz, ilv, irot, gb, ge, scatod)
      END DO
      !$OMP END DO NOWAIT
    END IF
  END DO
END DO
!$OMP END PARALLEL
! ----------------------------------------------------
phisnm(gb:ge, :) = ZERO
phimnm(:, gb:ge, :) = ZERO
IF (ljout) joutnm(:, gb:ge, :, :) = ZERO

DO ithr = 1, nThread
  DO ifsr = 1, nFsr
    DO ig = gb, ge
      phisnm(ig, ifsr) = phisnm(ig, ifsr) + TrackingDat(ithr)%phisnm(ig, ifsr)
        
      DO iod = 1, nod
        phimnm(iod, ig, ifsr) = phimnm(iod, ig, ifsr) + TrackingDat(ithr)%phimnm(iod, ig, ifsr)
      END DO
    END DO
  END DO
  
  IF (.NOT. ljout) CYCLE
  
  DO ipin = 1, nxy
    DO ibd = 1, nbd
      DO ig = gb, ge
        joutnm(:, ig, ibd, ipin) = joutnm(:, ig, ibd, ipin) + TrackingDat(ithr)%joutnm(:, ig, ibd, ipin)
      END DO
    END DO
  END DO
END DO
! ----------------------------------------------------
Cell => CoreInfo%CellInfo
Pin  => CoreInfo%Pin

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(FsrIdxSt, icel, jfsr)
!$OMP DO SCHEDULE(GUIDED)
DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt
  icel     = Pin(ipin)%Cell(iz)
  
  DO ifsr = 1, Cell(icel)%nFsr
    jfsr = FsrIdxSt + ifsr - 1
    
    DO ig = gb, ge
      wttmp = ONE / (xstnm(ig, jfsr) * Cell(icel)%vol(ifsr))
      
      phisnm(ig, jfsr) = phisnm(ig, jfsr) * wttmp + srcnm(ig, jfsr)
      
      phimnm(1:2, ig, jfsr) = phimnm(1:2, ig, jfsr) * wttmp + srcmnm(1:2, ig, jfsr) / 3._8
      
      IF (ScatOd .LT. 2) CYCLE
      
      phimnm(3:5, ig, jfsr) = phimnm(3:5, ig, jfsr) * wttmp + srcmnm(3:5, ig, jfsr) / 5._8
      
      IF (ScatOd .LT. 3) CYCLE
      
      phimnm(6:9, ig, jfsr) = phimnm(6:9, ig, jfsr) * wttmp + srcmnm(6:9, ig, jfsr) / 7._8
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
! ----------------------------------------------------
NULLIFY (Cell)
NULLIFY (Pin)
NULLIFY (MultigridInfo)
NULLIFY (Comp)
! ----------------------------------------------------

END SUBROUTINE RayTraceP1NM_OMP
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayP1NM_OMP(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, ilv, irot, gb, ge, scatod)

USE PARAM
USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, Asy_Type, PinInfo_Type, Cell_Type, AsyRayInfo_type, CoreRayInfo_Type, RotRayInfo_Type, CellRayInfo_type, TrackingDat_Type

IMPLICIT NONE

TYPE(RayInfo_Type)     :: RayInfo
TYPE(CoreInfo_Type)    :: CoreInfo
TYPE(TrackingDat_Type) :: TrackingDat

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz, ilv, irot, gb, ge, scatod

TYPE(Pin_Type),         POINTER, DIMENSION(:) :: Pin
TYPE(Asy_Type),         POINTER, DIMENSION(:) :: Asy
TYPE(PinInfo_Type),     POINTER, DIMENSION(:) :: PinInfo
TYPE(Cell_Type),        POINTER, DIMENSION(:) :: Cell
TYPE(AsyRayInfo_type),  POINTER, DIMENSION(:) :: AsyRay
TYPE(CoreRayInfo_Type), POINTER, DIMENSION(:) :: CoreRay
TYPE(RotRayInfo_Type),  POINTER, DIMENSION(:) :: RotRay

TYPE(CellRayInfo_Type), POINTER :: CellRay

INTEGER :: iazi, ipol, icray, jcray, iaray, jaray, ipray, iceray, irayseg, idir, ipin, icel, iasy, ifsr, isurf, jbeg, jend, jinc, ig, ibcel, iod, nod
INTEGER :: nRotRay, nCoreRay, nAsyRay, nPinRay, nRaySeg, nPolarAng, FsrIdxSt, PhiAnginSvIdx, PhiAngOutSvIdx, ExpAppIdx

INTEGER :: mp(2)

INTEGER, POINTER, DIMENSION(:) :: LocalFsrIdx

REAL :: PhiAngOut(RayInfo%nPolarAngle, gb : ge)
REAL :: phid, tau, local_src, ExpApp

REAL :: wt(10), wt2(10, 4)

REAL, POINTER, DIMENSION(:)       :: LenSeg
REAL, POINTER, DIMENSION(:,:)     :: phis, src, xst, expa, expb, wtang, wtazi
REAL, POINTER, DIMENSION(:,:,:)   :: phim, PhiAngIn, wtsurf, mwt, mwt2
REAL, POINTER, DIMENSION(:,:,:,:) :: jout, SrcAngnm1, SrcAngnm2

DATA mp /2, 1/
! ----------------------------------------------------

AsyRay  => RayInfo%AsyRay
CoreRay => RayInfo%CoreRay
RotRay  => RayInfo%RotRay

wtang    => RayInfo%MultigridInfo(ilv)%wtang
wtsurf   => RayInfo%MultigridInfo(ilv)%wtsurf
nPolarAng = RayInfo%MultigridInfo(ilv)%nPolar
mwt      => RayInfo%MultigridInfo(ilv)%mwt
mwt2     => RayInfo%MultigridInfo(ilv)%mwt2

Asy     => CoreInfo%Asy
Pin     => CoreInfo%Pin
PinInfo => CoreInfo%Pininfo
Cell    => CoreInfo%CellInfo

phis      => TrackingDat%phisnm
phim      => TrackingDat%phimnm
src       => TrackingDat%srcnm
xst       => TrackingDat%xstnm
jout      => TrackingDat%joutnm
PhiAngIn  => TrackingDat%phiAngInnm
EXPA      => TrackingDat%EXPA
EXPB      => TrackingDat%EXPB
SrcAngnm1 => TrackingDat%SrcAngnm1
SrcAngnm2 => TrackingDat%SrcAngnm2

nCoreRay = RotRay(irotRay)%nRay

PhiAngInSvIdx  = RayInfo%PhiAngInSvIdx (iRotRay, irot)
PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay, irot)

PhiAngOut = PhiAngIn(:, gb:ge, PhiAnginSvIdx)

nod = 2
IF(ScatOd .EQ. 2) nod = 5
IF(ScatOd .EQ. 3) nod = 9

IF (irot .EQ. 1) THEN
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
  IF (irot .eq. 2) idir = mp(idir)
  
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
            
            DO ipol = 1, nPolarAng
              ExpApp = EXPA(ipol, ExpAppIdx) * tau + EXPB(ipol, ExpAppIdx)
              phid   = (PhiAngOut(ipol, ig) - SrcAngnm1(ig, ipol, ifsr, iazi)) * ExpApp ! NOTICE : 1
              
              PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
              
              phis(ig, ifsr) = phis(ig, ifsr) + wt(ipol) * phid
              
              DO iod = 1, nod
                phim(iod, ig, ifsr) = phim(iod, ig, ifsr) + mwt(iod, ipol, iazi) * phid ! NOTICE : 1
              END DO
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
            
            DO ipol = 1, nPolarAng
              ExpApp = EXPA(ipol, ExpAppIdx) * tau + EXPB(ipol, ExpAppIdx)
              
              phid = (PhiAngOut(ipol, ig) - SrcAngnm2(ig, ipol, ifsr, iazi)) * ExpApp ! NOTICE : 2
              
              PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
              
              phis(ig, ifsr) = phis(ig, ifsr) + wt(ipol) * phid
              
              DO iod = 1, nod
                phim(iod, ig, ifsr) = phim(iod, ig, ifsr) + mwt2(iod, ipol, iazi) * phid ! NOTICE : 2
              END DO
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

END SUBROUTINE RecTrackRotRayP1NM_OMP
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayP1NM_OMP(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, ilv, irot, gb, ge, ScatOd)

USE PARAM
USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, PolarAngle_Type, TrackingDat_Type, Pin_Type
USE Moc_Mod, ONLY : nMaxCellRay, nMaxCoreRay
USE HexData, ONLY : hAsy
USE HexType, ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexCoreRay, Type_HexRotRay
USE HexData, ONLY : haRay, hCelBss, hCel, hLgc, hcRay, hRotRay, hAsyTypInfo

IMPLICIT NONE

TYPE(RayInfo_Type)     :: RayInfo
TYPE(CoreInfo_Type)    :: CoreInfo
TYPE(TrackingDat_Type) :: TrackingDat
LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz, ilv, irot, gb, ge, ScatOd

TYPE(Pin_Type), POINTER, DIMENSION(:) :: Pin

TYPE(Type_HexAsyRay), POINTER :: haRay_Loc
TYPE(Type_HexCelRay), POINTER :: CelRay_Loc
TYPE(Type_HexRotRay), POINTER :: hRotRay_Loc

INTEGER :: iAzi, iPol, icRay, jcRay, iAsyRay, iRaySeg, imRay, ihpRay, iAsy, ifsr, iSurf, jbeg, jend, jinc, ig, iGeoTyp, iAsyTyp, jhPin, icBss, jcBss, iod
INTEGER :: nCoreRay, nAsyRay, nPolarAng, PhiAnginSvIdx, PhiAngOutSvIdx, ExpAppIdx, nod

REAL :: phid, tau, ExpApp

REAL :: wt(10)
REAL :: PhiAngOut(RayInfo%nPolarAngle, gb:ge)

REAL, POINTER, DIMENSION(:,:)     :: phis, src, xst, expa, expb, wtang
REAL, POINTER, DIMENSION(:,:,:)   :: phim, PhiAngIn, mwt, mwt2
REAL, POINTER, DIMENSION(:,:,:,:) :: jout, SrcAngnm1, SrcAngnm2
! ----------------------------------------------------

wtang    => RayInfo%MultigridInfo(ilv)%wtang
!wtsurf  => RayInfo%MultigridInfo(ilv)%wtsurf
nPolarAng = RayInfo%MultigridInfo(ilv)%nPolar
mwt      => RayInfo%MultigridInfo(ilv)%mwt
mwt2     => RayInfo%MultigridInfo(ilv)%mwt2

hRotRay_Loc => hRotRay(iRotRay)
nCoreRay     = hRotRay_Loc%ncRay

Pin => CoreInfo%Pin

phis      => TrackingDat%phisnm
phim      => TrackingDat%phimnm
src       => TrackingDat%srcnm
xst       => TrackingDat%xstnm
jout      => TrackingDat%joutnm
PhiAngIn  => TrackingDat%phiAngInnm
EXPA      => TrackingDat%EXPA
EXPB      => TrackingDat%EXPB
SrcAngnm1 => TrackingDat%SrcAngnm1
SrcAngnm2 => TrackingDat%SrcAngnm2

PhiAngInSvIdx  = RayInfo%PhiAngInSvIdx (iRotRay, irot)
PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay, irot)

PhiAngOut = PhiAngIn(:, gb:ge, PhiAnginSvIdx)

nod = 2
IF(ScatOd .EQ. 2) nod = 5
IF(ScatOd .EQ. 3) nod = 9

IF (iRot .EQ. 1) THEN
  jbeg = 1; jend = nCoreRay; jinc = 1
ELSE
  jbeg = nCoreRay; jend = 1; jinc = -1
END IF
! ----------------------------------------------------
DO icRay = jbeg, jend, jinc
  jcRay   = hRotRay_Loc%cRayIdx(icRay)
  nAsyRay = hcRay(abs(jcRay))%nmRay
  iAzi    = hcRay(abs(jcRay))%AzmIdx
  
  IF (iRot .eq. 2) jcRay = -jcRay !Reverse the Sweep Direction
  
  DO iPol = 1, nPolarAng
    wt(ipol) = wtang(iPol, iAzi)
    
    !IF (lJout) wt2(ipol, 1 : 4) = wtsurf(ipol, iazi, 1 : 4)
  END DO
  ! --------------------------------------------------
  IF(jcRay > 0) THEN
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
            DO iPol = 1, nPolarAng
              Jout(1, ig, iSurf, jhPin) = Jout(1, ig, isurf, jhPin) + wt(ipol) * PhiAngOut(iPol, ig)
              !Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol, ig)
            END DO
          END DO
        END IF
        
        DO iRaySeg = 1, CelRay_Loc%nSegRay
          ifsr = CelRay_Loc%MshIdx(iRaySeg) + Pin(jhPin)%FsrIdxSt - 1
          
          DO ig = gb, ge
            tau = -CelRay_Loc%SegLgh(iRaySeg) * xst(ig, ifsr) ! Optimum Length
            
            ExpAppIdx = max(INT(tau), -40000)
            ExpAppIdx = min(0, ExpAppIdx)
            
            DO iPol = 1, nPolarAng
              ExpApp = EXPA(iPol, ExpAppIdx) * tau + EXPB(iPol, ExpAppIdx)
              
              phid = (PhiAngOut(iPol, ig) - SrcAngnm1(ig, ipol, ifsr, iazi)) * ExpApp ! NOTICE
              
              PhiAngOut(iPol, ig) = PhiAngOut(iPol, ig) - phid
              
              phis(ig, ifsr) = phis(ig, ifsr) + wt(iPol) * phid
              
              DO iod = 1, nod
                phim(iod, ig, ifsr) = phim(iod, ig, ifsr) + mwt(iod, ipol, iazi) * phid ! NOTICE
              END DO
            END DO
          END DO
        END DO
        
        IF (lJout) THEN
          isurf = CelRay_Loc%SurfIdx(2) ! y : Big
          
          DO ig = gb, ge
            DO iPol = 1, nPolarAng
              Jout(2, ig, iSurf, jhPin) = Jout(2, ig, iSurf, jhPin) + wt(iPol) * PhiAngOut(iPol, ig)
              
              !Jout(3, ig, isurf, jhPin) = Jout(3, ig, isurf, jhPin) + wt2(ipol, isurf) * PhiAngOut(ipol, ig)
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
            DO iPol = 1, nPolarAng
              Jout(1, ig, iSurf, jhPin) = Jout(1, ig, isurf, jhPin) + wt(ipol) * PhiAngOut(iPol, ig)
              
              !Jout(3, ig, isurf, jhPin) = Jout(3, ig, isurf, jhpin) + wt2(ipol, isurf) * PhiAngOut(ipol, ig)
            END DO
          END DO
        END IF
        
        DO iRaySeg = CelRay_Loc%nSegRay, 1, -1
          ifsr = CelRay_Loc%MshIdx(iRaySeg) + Pin(jhPin)%FsrIdxSt - 1
          
          DO ig = gb, ge
            tau = -CelRay_Loc%SegLgh(iRaySeg) * XsT(ig, ifsr) ! Optimum Length
            
            ExpAppIdx = max(INT(tau), -40000)
            ExpAppIdx = min(0, ExpAppIdx)
            
            DO iPol = 1, nPolarAng
              ExpApp = EXPA(iPol, ExpAppIdx) * tau + EXPB(iPol, ExpAppIdx)
              
              phid = (PhiAngOut(iPol, ig) - SrcAngnm2(ig, ipol, ifsr, iazi)) * ExpApp
              
              PhiAngOut(iPol, ig) = PhiAngOut(iPol, ig) - phid
              
              phis(ig, ifsr) = phis(ig, ifsr) + wt(iPol) * phid
              
              DO iod = 1, nod
                phim(iod, ig, ifsr) = phim(iod, ig, ifsr) + mwt2(iod, ipol, iazi) * phid ! NOTICE : 2
              END DO
            END DO
          END DO
        END DO
        
        IF (lJout) THEN
          isurf = CelRay_Loc%SurfIdx(1) ! y : small
          
          DO ig = gb, ge
            DO iPol = 1, nPolarAng
              Jout(2, ig, iSurf, jhPin) = Jout(2, ig, iSurf, jhPin) + wt(iPol) * PhiAngOut(iPol, ig)
              !Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol, ig)
            END DO
          END DO
        END IF
      END DO
    END DO
  END IF
END DO

PhiAngIn(:, gb:ge, PhiAngOutSvIdx) = PhiAngOut
! ----------------------------------------------------

END SUBROUTINE HexTrackRotRayP1NM_OMP
! ------------------------------------------------------------------------------------------------------------