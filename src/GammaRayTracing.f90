#include <defines.h>
#ifdef __GAMMA_TRANSPORT
!--- CNJ Edit : Ray Tracing Routines for Gamma MOC Calculation
SUBROUTINE RayTraceGamma(RayInfo, CoreInfo, phisnm, PhiAngInnm, xstnm, srcnm,                 &
                         joutnm, iz, gb, ge, ljout)
USE PARAM
USE TYPEDEF,            ONLY : RayInfo_Type,    CoreInfo_type,      Pin_Type,        Cell_Type
USE Moc_Mod,            ONLY : Expa_p,          Expb_p
USE GamMOC_MOD,         ONLY : TrackingDat
!USE Core_mod,           ONLY : GroupInfo
USE TIMER
USE ALLOCS
USE PE_MOD,             ONLY : PE
USE OMP_LIB
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phisnm(:, :), PhiAngInnm(:, :, :), xstnm(:, :), srcnm(:, :), joutnm(:, :, :, :)
INTEGER :: iz, gb, ge
LOGICAL :: ljout

TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(Pin_Type), POINTER :: Pin(:)

INTEGER :: ng, nAziAng, nPolarAng, nPhiAngSv, nRotray, nFsr, nxy, nThread
INTEGER :: i, j, l, g, icel, ireg, iRotRay, iAsy, iAsyRay, iazi, ipol, irot, tid
INTEGER :: FsrIdxSt

nAziAng = RayInfo%nAziAngle; nPolarAng = RayInfo%nPolarAngle
nPhiAngSv = RayInfo%nPhiAngSv
nRotRay = RayInfo%nRotRay
nFsr = CoreInfo%nCoreFsr; nxy = CoreInfo%nxy
!ng = GroupInfo%ngg
nThread = PE%nThread

DO tid = 1, nThread
  TrackingDat(tid)%Expa => Expa_p; TrackingDat(tid)%Expb => Expb_p
  TrackingDat(tid)%srcnm => srcnm; TrackingDat(tid)%xstnm => xstnm
  TrackingDat(tid)%PhiAngInnm => PhiAngInnm
ENDDO

phisnm(gb : ge, :) = zero
IF(ljout) joutnm(:, gb : ge, :, :) = zero

CALL omp_set_num_threads(nThread)

!!$OMP PARALLEL PRIVATE(tid, iRotRay)
tid = omp_get_thread_num() + 1

TrackingDat(tid)%phisnm = zero
IF (ljout) THEN
  TrackingDat(tid)%joutnm = zero
ENDIF

DO irot = 1, 2
  DO iazi = 1, nAziAng / 2
    ! !$OMP DO SCHEDULE(GUIDED)
    DO i = 1, RayInfo%RotRayAziList(0, iazi)
      iRotRay = RayInfo%RotRayAziList(i, iazi)
      CALL TrackRotRayGamma(RayInfo, CoreInfo, TrackingDat(tid), ljout, iRotRay, iz, irot, gb, ge)
    ENDDO
    ! !$OMP END DO NOWAIT
  ENDDO
ENDDO

!!$OMP END PARALLEL

DO j = 1, nThread
  phisnm = phisnm + TrackingDat(j)%phisnm
  IF (ljout) THEN
    joutnm = joutnm + TrackingDat(j)%joutnm
  ENDIF
ENDDO

Cell => CoreInfo%CellInfo; Pin => CoreInfo%Pin
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(FsrIdxSt, icel, ireg)
!$OMP DO SCHEDULE(GUIDED)
DO l = 1, nxy
  FsrIdxSt = Pin(l)%FsrIdxSt; icel = Pin(l)%Cell(iz);
  DO j = 1, Cell(icel)%nFsr
    ireg = FsrIdxSt + j - 1
    DO g = gb, ge
      phisnm(g, ireg) = phisnm(g, ireg) / (xstnm(g, ireg) * Cell(icel)%vol(j)) + srcnm(g, ireg)
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE TrackRotRayGamma(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, irot, gb, ge)
USE PARAM
USE TYPEDEF,        ONLY : RayInfo_Type,        Coreinfo_type,      Pin_Type,           Asy_Type,                   &
                           PinInfo_Type,        Cell_Type,          AziAngleInfo_Type,  PolarAngle_Type,            &
                           AsyRayInfo_type,     CoreRayInfo_Type,   RotRayInfo_Type,    CellRayInfo_type
USE GammaTYPEDEF,   ONLY : GammaTrackingDat_Type
USE cntl,           ONLY : nTracerCntl
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
TYPE(GammaTrackingDat_Type) :: TrackingDat
LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz, irot, gb, ge

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

INTEGER, POINTER :: LocalFsrIdx(:)
REAL, POINTER :: LenSeg(:)
REAL, POINTER :: phis(:, :), src(:, :), xst(:, :), jout(:, :, :, :)
REAL, POINTER :: PhiAngIn(:, :, :)
REAL, POINTER :: EXPA(:, :), EXPB(:, :), wtang(:, :), wtsurf(:, :, :)

INTEGER :: mp(2)
INTEGER :: iazi, ipol, iCoreRay, iasyray, iceray, irayseg, idir
INTEGER :: nRotRay, nCoreRay, nAsyRay, nPinRay, nRaySeg, FsrIdxSt
INTEGER :: nPolarAng, nAziAng
INTEGER :: ipin, icel, iasy, ireg, isurf
INTEGER :: PhiAnginSvIdx, PhiAngOutSvIdx
INTEGER :: i, j, k, l, m, jbeg, jend, jinc, ir, ir1, ig
INTEGER :: ibcel

REAL :: wttemp
REAL :: wt(10), wt2(10, 4)

REAL :: PhiAngOut(RayInfo%nPolarAngle, gb : ge)
REAL :: phid, tau, local_src, ExpApp
INTEGER :: ExpAppIdx

LOGICAL :: lFast

DATA mp /2, 1/

IF(nTracerCntl%lHex) THEN 
    CALL HexTrackRotRayGamma(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, irot, gb, ge)
    RETURN
END IF

AziAng => RayInfo%AziAngle;
PolarAng => RayInfo%PolarAngle;
AsyRay => RayInfo%AsyRay
CoreRay => RayInfo%CoreRay
RotRay => RayInfo%RotRay
nPolarAng = RayInfo%nPolarAngle
nAziAng = RayInfo%nAziAngle
Asy => CoreInfo%Asy
Pin => CoreInfo%Pin
PinInfo => CoreInfo%Pininfo
Cell => CoreInfo%CellInfo

phis => TrackingDat%phisnm; src => TrackingDat%srcnm
xst => TrackingDat%xstnm; jout => TrackingDat%joutnm
PhiAngIn => TrackingDat%PhiAngInnm
EXPA => TrackingDat%EXPA
EXPB => TrackingDat%EXPB
wtang => TrackingDat%wtang
wtsurf => TrackingDat%wtsurf

nCoreRay = RotRay(irotRay)%nRay

! DO irot = 1, 2
  PhiAngInSvIdx = RayInfo%PhiAngInSvIdx(iRotRay, irot)
  PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay, irot)
  PhiAngOut = PhiAngIn(:, gb : ge, PhiAnginSvIdx)
  jbeg = 1; jend = nCoreRay; jinc = 1
  IF (irot .EQ. 2) THEN
    jbeg = nCoreRay; jend = 1; jinc = -1
  ENDIF
  DO j = jbeg, jend, jinc
    iCoreRay = RotRay(iRotRay)%RayIdx(j)
    nAsyRay = CoreRay(iCoreRay)%nRay
    idir = RotRay(iRotRay)%dir(j); iazi = CoreRay(iCoreRay)%iang
    IF (irot .eq. 2) idir = mp(idir)
    DO ipol = 1, nPolarAng
      wt(ipol) = wtang(ipol, iazi)
      IF (lJout) wt2(ipol, 1 : 4) = wtsurf(ipol, iazi, 1 : 4)
    ENDDO
    IF (idir .EQ. 1) THEN
      DO k = 1, nAsyRay
        iasyray = CoreRay(iCoreRay)%AsyRayIdx(k)
        iasy = CoreRay(iCoreRay)%AsyIdx(k)
        IF (iAsy .EQ. 0) CYCLE
        nPinRay = AsyRay(iAsyRay)%nCellRay
        DO l = 1, nPinRay
          ipin = AsyRay(iAsyRay)%PinIdx(l)
          iceray = AsyRay(iAsyRay)%PinRayIdx(l)
          ipin = Asy(iAsy)%GlobalPinIdx(ipin)
          icel = Pin(ipin)%Cell(iz)
          FsrIdxSt = Pin(ipin)%FsrIdxSt
          ibcel = Cell(icel)%basecellstr
          CellRay => Cell(ibcel)%CellRay(iceray)
          nRaySeg = CellRay%nSeg
          LocalFsrIdx => CellRay%LocalFsrIdx
          LenSeg => CellRay%LenSeg
          IF (lJout) THEN
            isurf = AsyRay(iAsyRay)%PinRaySurf(1, l)
            DO ig = gb, ge
              DO ipol = 1, nPolarAng
                Jout(1, ig, isurf, ipin) = Jout(1, ig, isurf, ipin) + wt(ipol) * PhiAngOut(ipol, ig)
                Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol, ig)
              ENDDO
            ENDDO
          ENDIF
          DO iRaySeg = 1, nRaySeg
            ireg = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
            DO ig = gb, ge
              tau = - LenSeg(iRaySeg) * xst(ig, ireg)
              ExpAppIdx = max(INT(tau), -40000)
              ExpAppIdx = min(0, ExpAppIdx)
              local_src = src(ig, ireg)
              DO ipol = 1, nPolarAng
                ExpApp = EXPA(ipol, ExpAppIdx) * tau + EXPB(ipol, ExpAppIdx)
                phid = (PhiAngOut(ipol, ig) - local_src) * ExpApp
                PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
                phis(ig, ireg) = phis(ig, ireg) + wt(ipol) * phid
              ENDDO
            ENDDO
          ENDDO
          IF (lJout) THEN
            isurf = AsyRay(iAsyRay)%PinRaySurf(2, l)
            DO ig = gb, ge
              DO ipol = 1, nPolarAng
                Jout(2, ig, isurf, ipin) = Jout(2, ig, isurf, ipin) + wt(ipol) * PhiAngOut(ipol, ig)
                Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol, ig)
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    ELSE
      DO k = nAsyRay, 1, -1
        iasyray = CoreRay(iCoreRay)%AsyRayIdx(k)
        iasy = CoreRay(iCoreRay)%AsyIdx(k)
        IF (iAsy .EQ. 0) CYCLE
        nPinRay = AsyRay(iAsyRay)%nCellRay
        DO l = nPinRay, 1, -1
          ipin = AsyRay(iAsyRay)%PinIdx(l)
          iceray = AsyRay(iAsyRay)%PinRayIdx(l)
          ipin = Asy(iAsy)%GlobalPinIdx(ipin)
          icel = Pin(ipin)%Cell(iz)
          FsrIdxSt = Pin(ipin)%FsrIdxSt
          ibcel = Cell(icel)%basecellstr
          CellRay => Cell(ibcel)%CellRay(iceray)
          nRaySeg = CellRay%nSeg
          LocalFsrIdx => CellRay%LocalFsrIdx
          LenSeg => CellRay%LenSeg
          IF (lJout) THEN
            isurf = AsyRay(iAsyRay)%PinRaySurf(2, l)
            DO ig = gb, ge
              DO ipol = 1, nPolarAng
                Jout(1, ig, isurf, ipin) = Jout(1, ig, isurf, ipin) + wt(ipol) * PhiAngOut(ipol, ig)
                Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol, ig)
              ENDDO
            ENDDO
          ENDIF
          DO iRaySeg = nRaySeg, 1, -1
            ireg = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
            DO ig = gb, ge
              tau = - LenSeg(iRaySeg) * xst(ig, ireg)
              ExpAppIdx = max(INT(tau), -40000)
              ExpAppIdx = min(0, ExpAppIdx)
              local_src = src(ig, ireg)
              DO ipol = 1, nPolarAng
                ExpApp = EXPA(ipol, ExpAppIdx) * tau + EXPB(ipol, ExpAppIdx)
                phid = (PhiAngOut(ipol, ig) - local_src) * ExpApp
                PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
                phis(ig, ireg) = phis(ig, ireg) + wt(ipol) * phid
              ENDDO
            ENDDO
          ENDDO
          IF (lJout) THEN
            isurf = AsyRay(iAsyRay)%PinRaySurf(1, l)
            DO ig = gb, ge
              DO ipol = 1, nPolarAng
                Jout(2, ig, isurf, ipin) = Jout(2, ig, isurf, ipin) + wt(ipol) * PhiAngOut(ipol, ig)
                Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol, ig)
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDDO
  PhiAngIn(:, gb : ge, PhiAngOutSvIdx) = PhiAngOut
! ENDDO

END SUBROUTINE

SUBROUTINE RayTraceGamma_Pn(RayInfo, CoreInfo, phis, phim, PhiAngIn, xst, src, srcm, jout, iz, ScatOd, lJout)
USE PARAM
USE TYPEDEF,        ONLY : RayInfo_Type,    CoreInfo_type,      Pin_Type,           Cell_Type
USE Moc_Mod,        ONLY : Expa_p,          Expb_p,             AziMap
USE GamMOC_MOD,     ONLY : TrackingDat,     Comp,               mwt,                wtang
USE ALLOCS
USE PE_MOD,         ONLY : PE
USE OMP_LIB
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phis(:), phim(:, :), PhiAngIn(:, :), xst(:)
REAL, POINTER :: src(:), srcm(:, :), jout(:, :, :)
INTEGER :: iz
LOGICAL :: ljout
INTEGER :: ScatOd

TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(Pin_Type), POINTER :: Pin(:)
REAL, POINTER :: phia(:, :, :), SrcAng(:, :, :)

INTEGER :: nAziAng, nPolarAng
INTEGER :: nxy, nFsr, nThread, nAzipThread
INTEGER :: tid, od, icel, ireg, ipin, iazi, ipol, irot, iRotRay
INTEGER :: AziIdx, FsrIdxSt
INTEGER :: i, j, l, k, m
INTEGER, POINTER :: AziList(:, :)
REAL :: wttemp, tempsrc

Cell => CoreInfo%CellInfo
Pin => CoreInfo%Pin
nxy = CoreInfo%nxy
nFsr = CoreInfo%nCoreFsr

nAziAng = RayInfo%nAziAngle
nPolarAng = RayInfo%nPolarAngle
nThread = PE%nThread

IF (ScatOd .EQ. 1) od = 2
IF (ScatOd .EQ. 2) od = 5
IF (ScatOd .EQ. 3) od = 9

phis = zero
phim = zero
IF (lJout) jout = zero

DO tid = 1, nThread
  TrackingDat(tid)%Expa => Expa_p; TrackingDat(tid)%Expb => Expb_p
  TrackingDat(tid)%xst1g => xst
  TrackingDat(tid)%PhiAngIn1g => PhiAngIn
ENDDO

CALL omp_set_num_threads(nThread)

IF (mod(nAziAng / 2, nThread) .NE. 0) THEN
  PRINT *, '# Azimuthal Angles Should be a Multiple of # Threads'
  STOP
ENDIF

nAzipThread = nAziAng / 2 / nThread

ALLOCATE(AziList(nThread, nAzipThread))

iAzi = 0
DO i = 1, nThread
  DO j = 1, nAzipThread
    iAzi = iAzi + 1
    AziList(i, j) = iAzi
  ENDDO
ENDDO

DO i = 1, nAzipThread

  !$OMP PARALLEL PRIVATE(iRotRay, iAzi, AziIdx, SrcAng, tempsrc, tid)
  tid = omp_get_thread_num() + 1; iAzi = AziList(tid, i)
  SrcAng => TrackingDat(tid)%SrcAng1g

  IF (ScatOd .EQ. 1) THEN
    DO ireg = 1, nFsr
      AziIdx = iAzi
      DO ipol = 1, nPolarAng
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = src(ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = src(ireg)
        tempsrc = Comp(1, ipol, AziIdx) * srcm(1, ireg) + Comp(2, ipol, AziIdx) * srcm(2, ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = SrcAng(ipol, ireg, AziMap(AziIdx, 1)) + tempsrc
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = SrcAng(ipol, ireg, AziMap(AziIdx, 2)) - tempsrc
      ENDDO
      AziIdx = nAziAng - iAzi + 1
      DO ipol = 1, nPolarAng
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = src(ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = src(ireg)
        tempsrc = Comp(1, ipol, AziIdx) * srcm(1, ireg) + Comp(2, ipol, AziIdx) * srcm(2, ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = SrcAng(ipol, ireg, AziMap(AziIdx, 1)) + tempsrc
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = SrcAng(ipol, ireg, AziMap(AziIdx, 2)) - tempsrc
      ENDDO
    ENDDO
  ELSEIF (ScatOd .EQ. 2) THEN
    DO ireg = 1, nFsr
      AziIdx = iAzi
      DO ipol = 1, nPolarAng
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = src(ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = src(ireg)
        tempsrc = Comp(1, ipol, AziIdx) * srcm(1, ireg) + Comp(2, ipol, AziIdx) * srcm(2, ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = SrcAng(ipol, ireg, AziMap(AziIdx, 1)) + tempsrc
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = SrcAng(ipol, ireg, AziMap(AziIdx, 2)) - tempsrc
        tempsrc = Comp(3, ipol, AziIdx) * srcm(3, ireg) + Comp(4, ipol, AziIdx) * srcm(4, ireg)                     &
                  + Comp(5, ipol, AziIdx) * srcm(5, ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = SrcAng(ipol, ireg, AziMap(AziIdx, 1)) + tempsrc
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = SrcAng(ipol, ireg, AziMap(AziIdx, 2)) + tempsrc
      ENDDO
      AziIdx = nAziAng - iAzi + 1
      DO ipol = 1, nPolarAng
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = src(ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = src(ireg)
        tempsrc = Comp(1, ipol, AziIdx) * srcm(1, ireg) + Comp(2, ipol, AziIdx) * srcm(2, ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = SrcAng(ipol, ireg, AziMap(AziIdx, 1)) + tempsrc
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = SrcAng(ipol, ireg, AziMap(AziIdx, 2)) - tempsrc
        tempsrc = Comp(3, ipol, AziIdx) * srcm(3, ireg) + Comp(4, ipol, AziIdx) * srcm(4, ireg)                     &
                  + Comp(5, ipol, AziIdx) * srcm(5, ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = SrcAng(ipol, ireg, AziMap(AziIdx, 1)) + tempsrc
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = SrcAng(ipol, ireg, AziMap(AziIdx, 2)) + tempsrc
      ENDDO
    ENDDO
  ELSEIF (ScatOd .EQ. 3) THEN
    DO ireg = 1, nFsr
      AziIdx = iAzi
      DO ipol = 1, nPolarAng
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = src(ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = src(ireg)
        tempsrc = Comp(1, ipol, AziIdx) * srcm(1, ireg) + Comp(2, ipol, AziIdx) * srcm(2, ireg)                     &
                  + Comp(6, ipol, AziIdx) * srcm(6, ireg) + Comp(7, ipol, AziIdx) * srcm(7, ireg)                   &
                  + Comp(8, ipol, AziIdx) * srcm(8, ireg) + Comp(9, ipol, AziIdx) * srcm(9, ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = SrcAng(ipol, ireg, AziMap(AziIdx, 1)) + tempsrc
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = SrcAng(ipol, ireg, AziMap(AziIdx, 2)) - tempsrc
        tempsrc = Comp(3, ipol, AziIdx) * srcm(3, ireg) + Comp(4, ipol, AziIdx) * srcm(4, ireg)                     &
                  + Comp(5, ipol, AziIdx) * srcm(5, ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = SrcAng(ipol, ireg, AziMap(AziIdx, 1)) + tempsrc
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = SrcAng(ipol, ireg, AziMap(AziIdx, 2)) + tempsrc
      ENDDO
      AziIdx = nAziAng - iAzi + 1
      DO ipol = 1, nPolarAng
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = src(ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = src(ireg)
        tempsrc = Comp(1, ipol, AziIdx) * srcm(1, ireg) + Comp(2, ipol, AziIdx) * srcm(2, ireg)                     &
                  + Comp(6, ipol, AziIdx) * srcm(6, ireg) + Comp(7, ipol, AziIdx) * srcm(7, ireg)                   &
                  + Comp(8, ipol, AziIdx) * srcm(8, ireg) + Comp(9, ipol, AziIdx) * srcm(9, ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = SrcAng(ipol, ireg, AziMap(AziIdx, 1)) + tempsrc
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = SrcAng(ipol, ireg, AziMap(AziIdx, 2)) - tempsrc
        tempsrc = Comp(3, ipol, AziIdx) * srcm(3, ireg) + Comp(4, ipol, AziIdx) * srcm(4, ireg)                     &
                  + Comp(5, ipol, AziIdx) * srcm(5, ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = SrcAng(ipol, ireg, AziMap(AziIdx, 1)) + tempsrc
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = SrcAng(ipol, ireg, AziMap(AziIdx, 2)) + tempsrc
      ENDDO
    ENDDO
  ENDIF

  TrackingDat(tid)%phia1g = zero
  IF (lJout) TrackingDat(tid)%Jout1g = zero

  DO irot = 1, 2
    DO j = 1, RayInfo%RotRayAziList(0, iAzi)
      iRotRay = RayInfo%RotRayAziList(j, iAzi)
      CALL TrackRotRayGamma_Pn(RayInfo, CoreInfo, TrackingDat(tid), ljout, iRotRay, iz, irot)
    ENDDO
  ENDDO

  !$OMP END PARALLEL

  DO tid = 1, nThread

    iAzi = AziList(tid, i)
    phia => TrackingDat(tid)%phia1g

    IF (ScatOd .EQ. 1) THEN
    !$OMP PARALLEL DO PRIVATE(AziIdx) SCHEDULE(GUIDED)
      DO ireg = 1, nFsr
        AziIdx = iAzi
        DO ipol = 1, nPolarAng
          phis(ireg) = phis(ireg) + wtang(ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) + phia(ipol, ireg, AziMap(AziIdx, 2)))
          phim(1:2, ireg) = phim(1:2, ireg) + mwt(1:2, ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) - phia(ipol, ireg, AziMap(AziIdx, 2)))
        ENDDO
        AziIdx = nAziAng - iAzi + 1
        DO ipol = 1, nPolarAng
          phis(ireg) = phis(ireg) + wtang(ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) + phia(ipol, ireg, AziMap(AziIdx, 2)))
          phim(1:2, ireg) = phim(1:2, ireg) + mwt(1:2, ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) - phia(ipol, ireg, AziMap(AziIdx, 2)))
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
    ELSEIF (ScatOd .EQ. 2) THEN
      !$OMP PARALLEL DO PRIVATE(AziIdx) SCHEDULE(GUIDED)
      DO ireg = 1, nFsr
        AziIdx = iAzi
        DO ipol = 1, nPolarAng
          phis(ireg) = phis(ireg) + wtang(ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) + phia(ipol, ireg, AziMap(AziIdx, 2)))
          phim(1:2, ireg) = phim(1:2, ireg) + mwt(1:2, ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) - phia(ipol, ireg, AziMap(AziIdx, 2)))
          phim(3:5, ireg) = phim(3:5, ireg) + mwt(3:5, ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) + phia(ipol, ireg, AziMap(AziIdx, 2)))
        ENDDO
        AziIdx = nAziAng - iAzi + 1
        DO ipol = 1, nPolarAng
          phis(ireg) = phis(ireg) + wtang(ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) + phia(ipol, ireg, AziMap(AziIdx, 2)))
          phim(1:2, ireg) = phim(1:2, ireg) + mwt(1:2, ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) - phia(ipol, ireg, AziMap(AziIdx, 2)))
          phim(3:5, ireg) = phim(3:5, ireg) + mwt(3:5, ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) + phia(ipol, ireg, AziMap(AziIdx, 2)))
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
    ELSEIF (ScatOd .EQ. 3) THEN
      !$OMP PARALLEL DO PRIVATE(AziIdx) SCHEDULE(GUIDED)
      DO ireg = 1, nFsr
        AziIdx = iAzi
        DO ipol = 1, nPolarAng
          phis(ireg) = phis(ireg) + wtang(ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) + phia(ipol, ireg, AziMap(AziIdx, 2)))
          phim(1:2, ireg) = phim(1:2, ireg) + mwt(1:2, ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) - phia(ipol, ireg, AziMap(AziIdx, 2)))
          phim(3:5, ireg) = phim(3:5, ireg) + mwt(3:5, ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) + phia(ipol, ireg, AziMap(AziIdx, 2)))
          phim(6:9, ireg) = phim(6:9, ireg) + mwt(6:9, ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) - phia(ipol, ireg, AziMap(AziIdx, 2)))
        ENDDO
        AziIdx = nAziAng - iAzi + 1
        DO ipol = 1, nPolarAng
          phis(ireg) = phis(ireg) + wtang(ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) + phia(ipol, ireg, AziMap(AziIdx, 2)))
          phim(1:2, ireg) = phim(1:2, ireg) + mwt(1:2, ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) - phia(ipol, ireg, AziMap(AziIdx, 2)))
          phim(3:5, ireg) = phim(3:5, ireg) + mwt(3:5, ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) + phia(ipol, ireg, AziMap(AziIdx, 2)))
          phim(6:9, ireg) = phim(6:9, ireg) + mwt(6:9, ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) - phia(ipol, ireg, AziMap(AziIdx, 2)))
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
    ENDIF

    IF (lJout) Jout = Jout + TrackingDat(tid)%Jout1g

  ENDDO

ENDDO

IF (ScatOd .EQ. 1) THEN
  !$OMP PARALLEL DO PRIVATE(FsrIdxSt, ireg, icel, wttemp) SCHEDULE(GUIDED)
  DO l = 1, nxy
    FsrIdxSt = Pin(l)%FsrIdxSt; icel = Pin(l)%Cell(iz);
    DO j = 1, Cell(icel)%nFsr
      ireg = FsrIdxSt + j - 1
      wttemp = 1._8 / (xst(ireg) * Cell(icel)%vol(j))
      phis(ireg) = phis(ireg) * wttemp + src(ireg)
      phim(1:2, ireg) = phim(1:2, ireg) * wttemp + srcm(1:2, ireg) * rthree
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
ELSEIF (ScatOd .EQ. 2) THEN
  !$OMP PARALLEL DO PRIVATE(FsrIdxSt, ireg, icel, wttemp) SCHEDULE(GUIDED)
  DO l = 1, nxy
    FsrIdxSt = Pin(l)%FsrIdxSt; icel = Pin(l)%Cell(iz);
    DO j = 1, Cell(icel)%nFsr
      ireg = FsrIdxSt + j - 1
      wttemp = 1._8 / (xst(ireg) * Cell(icel)%vol(j))
      phis(ireg) = phis(ireg) * wttemp + src(ireg)
      phim(1:2, ireg) = phim(1:2, ireg) * wttemp + srcm(1:2, ireg) * rthree
      phim(3:5, ireg) = phim(3:5, ireg) * wttemp + srcm(3:5, ireg) * rfive
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
ELSEIF (ScatOd .EQ. 3) THEN
  !$OMP PARALLEL DO PRIVATE(FsrIdxSt, ireg, icel, wttemp) SCHEDULE(GUIDED)
  DO l = 1, nxy
    FsrIdxSt = Pin(l)%FsrIdxSt; icel = Pin(l)%Cell(iz);
    DO j = 1, Cell(icel)%nFsr
      ireg = FsrIdxSt + j - 1
      wttemp = 1._8 / (xst(ireg) * Cell(icel)%vol(j))
      phis(ireg) = phis(ireg) * wttemp + src(ireg)
      phim(1:2, ireg) = phim(1:2, ireg) * wttemp + srcm(1:2, ireg) * rthree
      phim(3:5, ireg) = phim(3:5, ireg) * wttemp + srcm(3:5, ireg) * rfive
      phim(6:9, ireg) = phim(6:9, ireg) * wttemp + srcm(6:9, ireg) * rseven
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
ENDIF

END SUBROUTINE

SUBROUTINE TrackRotRayGamma_Pn(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, irot)
USE PARAM
USE TYPEDEF,        ONLY : RayInfo_Type,        Coreinfo_type,      Pin_Type,           Asy_Type,                   &
                           PinInfo_Type,        Cell_Type,          AziAngleInfo_Type,  PolarAngle_Type,            &
                           AsyRayInfo_type,     CoreRayInfo_Type,   RotRayInfo_Type,    CellRayInfo_type
USE GammaTYPEDEF,   ONLY : GammaTrackingDat_Type
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
TYPE(GammaTrackingDat_Type) :: TrackingDat
LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz, irot

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

INTEGER, POINTER :: LocalFsrIdx(:), AziMap(:, :)
REAL, POINTER :: LenSeg(:)
REAL, POINTER :: phia(:, :, :), SrcAng(:, :, :), xst(:), jout(:, :, :)
REAL, POINTER :: PhiAngIn(:, :)
REAL, POINTER :: EXPA(:, :), EXPB(:, :)
REAL, POINTER :: wtang(:, :), wtsurf(:, :, :)

INTEGER :: mp(2)
INTEGER :: iazi, ipol, iCoreRay, iasyray, iceray, irayseg, idir
INTEGER :: nRotRay, nCoreRay, nAsyRay, nPinRay, nRaySeg, FsrIdxSt
INTEGER :: nPolarAng, nAziAng
INTEGER :: ipin, icel, iasy, ireg, isurf
INTEGER :: PhiAnginSvIdx, PhiAngOutSvIdx, AziSvIdx(2)
INTEGER :: i, j, k, l, m, jbeg, jend, jinc, ir, ir1
INTEGER :: ibcel

REAL :: wttemp
REAL :: wt(10), wt2(10, 4)

REAL :: PhiAngOut(RayInfo%nPolarAngle)
REAL :: phid, tau, ExpApp
INTEGER :: ExpAppIdx

LOGICAL :: lFast

DATA mp /2, 1/

AziAng => RayInfo%AziAngle;
PolarAng => RayInfo%PolarAngle;
AsyRay => RayInfo%AsyRay
CoreRay => RayInfo%CoreRay
RotRay => RayInfo%RotRay

Asy => CoreInfo%Asy
Pin => CoreInfo%Pin
PinInfo => CoreInfo%Pininfo
Cell => CoreInfo%CellInfo

phia => TrackingDat%phia1g; SrcAng => TrackingDat%SrcAng1g
xst => TrackingDat%xst1g; jout => TrackingDat%Jout1g
PhiAngIn => TrackingDat%PhiAngIn1g
EXPA => TrackingDat%EXPA
EXPB => TrackingDat%EXPB
AziMap => TrackingDat%AziMap
wtang => TrackingDat%wtang
wtsurf => TrackingDat%wtsurf

nPolarAng = RayInfo%nPolarAngle
nAziAng = RayInfo%nAziAngle
nCoreRay = RotRay(irotRay)%nRay

! DO irot = 1, 2
  PhiAngInSvIdx = RayInfo%PhiAngInSvIdx(iRotRay, irot)
  PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay, irot)
  PhiAngOut = PhiAngIn(:, PhiAnginSvIdx)
  jbeg = 1; jend = nCoreRay; jinc = 1
  IF (irot .EQ. 2) THEN
    jbeg = nCoreRay; jend = 1; jinc = -1
  ENDIF
  DO j = jbeg, jend, jinc
    iCoreRay = RotRay(iRotRay)%RayIdx(j)
    nAsyRay = CoreRay(iCoreRay)%nRay
    idir = RotRay(iRotRay)%dir(j); iazi = CoreRay(iCoreRay)%iang
    AziSvIdx = AziMap(iazi, :)
    IF (irot .eq. 2) idir = mp(idir)
    DO ipol = 1, nPolarAng
      wt(ipol) = wtang(ipol, iazi)
      IF (lJout) wt2(ipol, 1 : 4) = wtsurf(ipol, iazi, 1 : 4)
    ENDDO
    IF (idir .EQ. 1) THEN
      DO k = 1, nAsyRay
        iasyray = CoreRay(iCoreRay)%AsyRayIdx(k)
        iasy = CoreRay(iCoreRay)%AsyIdx(k)
        IF (iAsy .EQ. 0) CYCLE
        nPinRay = AsyRay(iAsyRay)%nCellRay
        DO l = 1, nPinRay
          ipin = AsyRay(iAsyRay)%PinIdx(l)
          iceray = AsyRay(iAsyRay)%PinRayIdx(l)
          ipin = Asy(iAsy)%GlobalPinIdx(ipin)
          icel = Pin(ipin)%Cell(iz)
          FsrIdxSt = Pin(ipin)%FsrIdxSt
          ibcel = Cell(icel)%basecellstr
          CellRay => Cell(ibcel)%CellRay(iceray)
          nRaySeg = CellRay%nSeg
          LocalFsrIdx => CellRay%LocalFsrIdx
          LenSeg => CellRay%LenSeg
          IF (lJout) THEN
            isurf = AsyRay(iAsyRay)%PinRaySurf(1, l)
            DO ipol = 1, nPolarAng
              Jout(1, isurf, ipin) = Jout(1, isurf, ipin) + wt(ipol) * PhiAngOut(ipol)
              Jout(3, isurf, ipin) = Jout(3, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol)
            ENDDO
          ENDIF
          DO iRaySeg = 1, nRaySeg
            ireg = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
            tau = - LenSeg(iRaySeg) * xst(ireg)
            ExpAppIdx = max(INT(tau), -40000)
            ExpAppIdx = min(0, ExpAppIdx)
            DO ipol = 1, nPolarAng
              ExpApp = EXPA(ipol, ExpAppIdx) * tau + EXPB(ipol, ExpAppIdx)
              phid = (PhiAngOut(ipol) - SrcAng(ipol, ireg, AziSvIdx(idir))) * ExpApp
              PhiAngOut(ipol) = PhiAngOut(ipol) - phid
              phia(ipol, ireg, AziSvIdx(idir)) = phia(ipol, ireg, AziSvIdx(idir)) + phid
            ENDDO
          ENDDO
          IF (lJout) THEN
            isurf = AsyRay(iAsyRay)%PinRaySurf(2, l)
            DO ipol = 1, nPolarAng
              Jout(2, isurf, ipin) = Jout(2, isurf, ipin) + wt(ipol) * PhiAngOut(ipol)
              Jout(3, isurf, ipin) = Jout(3, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    ELSE
      DO k = nAsyRay, 1, -1
        iasyray = CoreRay(iCoreRay)%AsyRayIdx(k)
        iasy = CoreRay(iCoreRay)%AsyIdx(k)
        IF (iAsy .EQ. 0) CYCLE
        nPinRay = AsyRay(iAsyRay)%nCellRay
        DO l = nPinRay, 1, -1
          ipin = AsyRay(iAsyRay)%PinIdx(l)
          iceray = AsyRay(iAsyRay)%PinRayIdx(l)
          ipin = Asy(iAsy)%GlobalPinIdx(ipin)
          icel = Pin(ipin)%Cell(iz)
          FsrIdxSt = Pin(ipin)%FsrIdxSt
          ibcel = Cell(icel)%basecellstr
          CellRay => Cell(ibcel)%CellRay(iceray)
          nRaySeg = CellRay%nSeg
          LocalFsrIdx => CellRay%LocalFsrIdx
          LenSeg => CellRay%LenSeg
          IF (lJout) THEN
            isurf = AsyRay(iAsyRay)%PinRaySurf(2, l)
            DO ipol = 1, nPolarAng
              Jout(1, isurf, ipin) = Jout(1, isurf, ipin) + wt(ipol) * PhiAngOut(ipol)
              Jout(3, isurf, ipin) = Jout(3, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol)
            ENDDO
          ENDIF
          DO iRaySeg = nRaySeg, 1, -1
            ireg = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
            tau = - LenSeg(iRaySeg) * xst(ireg)
            ExpAppIdx = max(INT(tau), -40000)
            ExpAppIdx = min(0, ExpAppIdx)
            DO ipol = 1, nPolarAng
              ExpApp = EXPA(ipol, ExpAppIdx) * tau + EXPB(ipol, ExpAppIdx)
              phid = (PhiAngOut(ipol) - SrcAng(ipol, ireg, AziSvIdx(idir))) * ExpApp
              PhiAngOut(ipol) = PhiAngOut(ipol) - phid
              phia(ipol, ireg, AziSvIdx(idir)) = phia(ipol, ireg, AziSvIdx(idir)) + phid
            ENDDO
          ENDDO
          IF (lJout) THEN
            isurf = AsyRay(iAsyRay)%PinRaySurf(1, l)
            DO ipol = 1, nPolarAng
              Jout(2, isurf, ipin) = Jout(2, isurf, ipin) + wt(ipol) * PhiAngOut(ipol)
              Jout(3, isurf, ipin) = Jout(3, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDDO
  PhiAngIn(:, PhiAngOutSvIdx) = PhiAngOut
! ENDDO

END SUBROUTINE

SUBROUTINE HexTrackRotRayGamma(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, irot, gb, ge)

USE PARAM
USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, PolarAngle_Type, Pin_Type
USE GammaTYPEDEF,  ONLY : GammaTrackingDat_Type
USE Moc_Mod,  ONLY : nMaxCellRay, nMaxCoreRay
USE HexData, ONLY : hAsy
USE HexType, ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexCoreRay, Type_HexRotRay
USE HexData, ONLY : haRay, hCelBss, hCel, hLgc, hcRay, hRotRay, hAsyTypInfo

IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
TYPE(GammaTrackingDat_Type) :: TrackingDat
LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz, irot, gb, ge

TYPE(Pin_Type),       POINTER :: Pin(:)
TYPE(Type_HexAsyRay), POINTER :: haRay_Loc
TYPE(Type_HexCelRay), POINTER :: CelRay_Loc
TYPE(Type_HexRotRay), POINTER :: hRotRay_Loc

REAL, POINTER :: phis(:, :), src(:, :), xst(:, :), jout(:, :, :, :)
REAL, POINTER :: PhiAngIn(:, :, :)
REAL, POINTER :: EXPA(:, :), EXPB(:, :)

INTEGER :: iAzi, iPol, icRay, jcRay, iAsyRay, iRaySeg, imRay, ihpRay
INTEGER :: nCoreRay, nAsyRay, nPolarAng
INTEGER :: iAsy, iReg, iSurf
INTEGER :: PhiAnginSvIdx, PhiAngOutSvIdx
INTEGER :: jbeg, jend, jinc, ig
INTEGER :: iGeoTyp, iAsyTyp, jhPin, icBss, jcBss

REAL, POINTER :: wtang(:, :)!, wtsurf(:, :, :)
REAL :: wt(10)!, wt2(10, 4)

REAL :: PhiAngOut(RayInfo%nPolarAngle, gb:ge)
REAL :: phid, tau, local_src, ExpApp
INTEGER :: ExpAppIdx
! ----------------------------------------------------

!Ray Info Pointing
wtang => TrackingDat%wtang
nPolarAng = RayInfo%nPolarAngle
hRotRay_Loc => hRotRay(iRotRay)

! Geometry Info Pointing
Pin => CoreInfo%Pin

!Tracking Dat Pointing
Phis     => TrackingDat%phisnm
src      => TrackingDat%srcnm
xst      => TrackingDat%xstnm
jout     => TrackingDat%joutnm
PhiAngIn => TrackingDat%phiAngInnm
EXPA     => TrackingDat%EXPA
EXPB     => TrackingDat%EXPB

nCoreRay       = hRotRay_Loc%ncRay
PhiAngInSvIdx  = RayInfo%PhiAngInSvIdx (iRotRay, irot)
PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay, irot)
PhiAngOut      = PhiAngIn(:, gb:ge, PhiAnginSvIdx)

jbeg = 1
jend = nCoreRay
jinc = 1

IF (iRot .EQ. 2) THEN
  jbeg = nCoreRay
  jend = 1
  jinc = -1
END IF
! ----------------------------------------------------
!                TRACK
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
  ! ----------------------------
  !      1. Forward
  ! ----------------------------
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
          iReg = CelRay_Loc%MshIdx(iRaySeg) + Pin(jhPin)%FsrIdxSt - 1
          
          DO ig = gb, ge
            tau       = -CelRay_Loc%SegLgh(iRaySeg) * XsT(ig, iReg) ! Optimum Length
            ExpAppIdx = max(INT(tau), -40000)
            ExpAppIdx = min(0, ExpAppIdx)
            local_src = src(ig, iReg)
            
            DO iPol = 1, nPolarAng
              ExpApp = EXPA(iPol, ExpAppIdx) * tau + EXPB(iPol, ExpAppIdx)
              phid   = (PhiAngOut(iPol, ig) - local_src) * ExpApp
              
              PhiAngOut(iPol, ig) = PhiAngOut(iPol, ig) - phid
              phis(ig, iReg)      = phis(ig, iReg) + wt(iPol) * phid
            END DO
          END DO
        END DO
        
        IF (lJout) THEN
          isurf = CelRay_Loc%SurfIdx(2) ! y : Big
          
          DO ig = gb, ge
            DO iPol = 1, nPolarAng
              Jout(2, ig, iSurf, jhPin) = Jout(2, ig, iSurf, jhPin) + wt(iPol) * PhiAngOut(iPol, ig)
              !Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol, ig)
            END DO
          END DO
        END IF
      END DO
    END DO
  ! ----------------------------
  !      2. Backward
  ! ----------------------------
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
              !Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol, ig)
            END DO
          END DO
        END IF
        
        DO iRaySeg = CelRay_Loc%nSegRay, 1, -1
          iReg = CelRay_Loc%MshIdx(iRaySeg) + Pin(jhPin)%FsrIdxSt - 1
          
          DO ig = gb, ge
            tau       = -CelRay_Loc%SegLgh(iRaySeg) * XsT(ig, iReg) ! Optimum Length
            ExpAppIdx = max(INT(tau), -40000)
            ExpAppIdx = min(0, ExpAppIdx)
            local_src = src(ig, iReg)
            
            DO iPol = 1, nPolarAng
              ExpApp = EXPA(iPol, ExpAppIdx) * tau + EXPB(iPol, ExpAppIdx)
              phid   = (PhiAngOut(iPol, ig) - local_src) * ExpApp
              
              PhiAngOut(iPol, ig) = PhiAngOut(iPol, ig) - phid
              phis(ig, iReg)      = phis(ig, iReg) + wt(iPol) * phid
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

END SUBROUTINE

#endif
