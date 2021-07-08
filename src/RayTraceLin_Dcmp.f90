#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceLin_Dcmp(RayInfo, CoreInfo, iz, gb, ge, lJout, lHybrid)

USE OMP_LIB
USE PARAM,       ONLY : ZERO, RED, BLACK, GREEN
USE TYPEDEF,     ONLY : RayInfo_Type, Coreinfo_type, Asy_Type, AsyInfo_Type
USE Moc_Mod,     ONLY : TrackingDat, phisNM, phimNM, srcNM, xstNM, MocjoutNM, PhiAngInNM, DcmpPhiAngInNg, DcmpPhiAngOutNg, &
                        RayTraceDcmp_OMP, RayTraceDcmp_Pn, RayTraceDcmp_LSCASMO, DcmpGatherBndyFlux, DcmpScatterBndyFlux, DcmpLinkBndyFlux
USE Core_mod,    ONLY : phisSlope, srcSlope
USE PE_MOD,      ONLY : PE
USE CNTL,        ONLY : nTracerCntl
USE itrcntl_mod, ONLY : itrcntl

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

INTEGER :: iz, gb, ge
LOGICAL :: lJout, lHybrid
! ----------------------------------------------------
TYPE (AsyInfo_Type), POINTER, DIMENSION(:) :: AsyInfo
TYPE (Asy_Type),     POINTER, DIMENSION(:) :: Asy

INTEGER :: color, startColor, endColor, colorInc, ithr, nThr, AsyType, iAsy
! ----------------------------------------------------

AsyInfo => CoreInfo%AsyInfo
Asy     => CoreInfo%Asy

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
  
  TrackingDat(ithr)%phisNM(gb:ge, :) = ZERO
  IF (ljout) TrackingDat(ithr)%JoutNM(:, gb:ge, :, :) = ZERO
  
  TrackingDat(ithr)%srcSlope => srcSlope(:, :, :, iz) ! Lin. CASMO
END DO

DcmpPhiAngOutNg(:, gb:ge, :, :, :) = ZERO
! ----------------------------------------------------
DO color = startColor, endColor, colorInc
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) CALL DcmpScatterBndyFlux(RayInfo, PhiAngInNM, DcmpPhiAngInNg)
#endif

  !$OMP PARALLEL PRIVATE(ithr, iAsy, AsyType)
  ithr = 1
  !$ ithr = omp_get_thread_num()+1
  
  TrackingDat(ithr)%PhiAngInNM    => PhiAngInNM
  TrackingDat(ithr)%DcmpPhiAngInNg  => DcmpPhiAngInNg
  TrackingDat(ithr)%DcmpPhiAngOutNg => DcmpPhiAngOutNg
  !$OMP BARRIER
  !$OMP DO SCHEDULE(GUIDED)
  DO iAsy = PE%myAsyBeg, PE%myAsyEnd
    IF (Asy(iAsy)%color .NE. color) CYCLE
    
    AsyType = Asy(iAsy)%AsyType
    
    IF (lHybrid .AND. AsyInfo(AsyType)%lFuel) THEN
      !CALL RayTraceDcmp_OMP(RayInfo, CoreInfo, iz, iAsy, gb, ge, ljout)
    ELSE
      CALL RayTraceDcmp_LSCASMO(RayInfo, CoreInfo, phisNM, phisSlope, MocjoutNM, iz, iAsy, gb, ge, ljout)
    END IF
  END DO
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL
  
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) CALL DcmpGatherBndyFlux(RayInfo, DcmpPhiAngOutNg)
#endif
  
  IF (PE%RTMASTER) CALL DcmpLinkBndyFlux(CoreInfo, RayInfo, PhiAngInNM, DcmpPhiAngInNg, DcmpPhiAngOutNg, gb, ge, color)
END DO

NULLIFY (Asy)
NULLIFY (AsyInfo)
! ----------------------------------------------------

END SUBROUTINE RayTraceLin_Dcmp
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmp_LSCASMO(RayInfo, CoreInfo, phisNM, phisSlope, joutNM, iz, iasy, gb, ge, ljout)

USE PARAM
USE ALLOCS
USE OMP_LIB
USE TYPEDEF,    ONLY :  RayInfo_Type,       Coreinfo_type,      Pin_Type,           Asy_Type,               &
                        AsyInfo_Type,       PinInfo_Type,       Cell_Type,          DcmpAsyRayInfo_Type
USE Moc_Mod,    ONLY :  nMaxRaySeg,         nMaxCellRay,        nMaxAsyRay,         nMaxCoreRay,            &
                        Expa,               Expb,               ApproxExp,          TrackingDat,            &
                        DcmpPhiAngInNg,       DcmpPhiAngOutNg,      TrackRotRayLSDcmp_CASMO
USE geom,       ONLY :  ng
USE PE_Mod,     ONLY :  PE

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
REAL, POINTER :: DcmpPhiAngOutNg(:, :, :, :, :), DcmpPhiAngInNg(:, :, :, :, :)
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
DcmpPhiAngInNg => TrackingDat%DcmpPhiAngInNg
DcmpPhiAngOutNg => TrackingDat%DcmpPhiAngOutNg
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
    phiobd(1 : nPolarAng, gb : ge) = DcmpPhiAngInNg(1 : nPolarAng, gb : ge, irot, iRay, iAsy)
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
  DcmpPhiAngOutNg(1 : nPolarAng, gb : ge, irot, iRay, iAsy) = phiobd(1 : nPolarAng, gb : ge)
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