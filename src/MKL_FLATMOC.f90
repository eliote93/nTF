#include <defines.h>
!--- CNJ Edit : 1D Axial Flat MOC Modules with Intel MKL
#ifdef __INTEL_MKL

#define FLUX_SHAPE

MODULE MKL_FLATMOC

USE MKL_3D
IMPLICIT NONE

REAL, POINTER, PRIVATE :: phisShape(:, :, :), phis(:, :, :), phim(:, :, :, :)
REAL, POINTER, PRIVATE :: psi(:, :), psiShape(:, :), src(:, :, :), srca(:, :, :, :, :), srcm(:, :, :, :), lkg(:, :, :)
REAL, POINTER, PRIVATE :: xst(:, :, :), pxs(:, :, :)
REAL, POINTER, PRIVATE :: S(:, :, :, :), F(:, :, :), Chi(:, :, :)

REAL, POINTER, PRIVATE :: wtExp(:, :, :, :), Comp(:, :), mwt(:, :)

REAL, POINTER, PRIVATE :: hzMOC(:)
INTEGER, POINTER, PRIVATE :: cmRange(:, :), fmRange(:, :)
INTEGER, POINTER, PRIVATE :: cmMap(:), fmMap(:)
INTEGER, PRIVATE :: nzMOC, nDiv(100)

LOGICAL, PRIVATE :: lFirst = .TRUE.

PRIVATE
PUBLIC :: AllocFlatMOC, FlatMOCDriver

CONTAINS

!--- Public Routines ------------------------------------------------------------------------------

SUBROUTINE AllocFlatMOC()

IMPLICIT NONE

INTEGER :: ng, nxy, nzCMFD, nPolarAngle

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
nPolarAngle = mklGeom%nPolar1D

CALL SetSphericalHarmonics()
CALL SetSubmesh()

ALLOCATE(phis(ng, nzMOC, nxy))
ALLOCATE(phim(3, ng, nzMOC, nxy)); phim = 0.0
#ifdef FLUX_SHAPE
ALLOCATE(phisShape(ng, nzMOC, nxy)); phisShape = 1.0
#endif
ALLOCATE(psi(nzMOC, nxy))
#ifndef FLUX_SHAPE
ALLOCATE(psiShape(nzMOC, nxy)); psiShape = 1.0
#endif
ALLOCATE(src(ng, nzMOC, nxy))
ALLOCATE(srca(2, nPolarAngle, ng, nzMOC, nxy)); srca = 0.0
ALLOCATE(srcm(3, ng, nzMOC, nxy)); srcm = 0.0
ALLOCATE(lkg(ng, nzMOC, nxy))
ALLOCATE(xst(ng, nzMOC, nxy))
ALLOCATE(pxs(ng, nzMOC, nxy))

ALLOCATE(wtExp(nPolarAngle, ng, nzMOC, nxy))

ALLOCATE(S(ng, ng, nzCMFD, nxy)); S = 0.0
ALLOCATE(F(ng, nzCMFD, nxy)); F = 0.0
ALLOCATE(Chi(ng, nzCMFD, nxy)); Chi = 0.0

END SUBROUTINE

SUBROUTINE FlatMOCDriver(PinXS, eigv)
USE TYPEDEF,        ONLY : PinXS_Type
USE PE_MOD,         ONLY : PE
USE TIMER,          ONLY : nTracer_dclock,      TimeChk
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)
REAL :: eigv

INTEGER :: i, ig, ixy, iz, izf, ng, nxy, nzCMFD, nNeg(1000)
INTEGER :: ierr, iter, itermax = 10
REAL :: AxNSolverBeg, AxNSolverEnd

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

CALL SetSourceOperator(PinXS)
CALL SetFlux()
CALL SetPsi()
CALL SetLeakage(PinXS)
CALL SetPseudoAbsorption()
CALL SetCrossSection(PinXS)
CALL SetConstant(PinXS)

! PRINT *, 'negative feed flux', mklGeom%myzb, mklGeom%myze, COUNT(phis .LT. 0.0)
! PRINT *, 'negative incoming', mklGeom%myzb, mklGeom%myze, COUNT(mklAxial%PhiAngIn .LT. 0.0)

DO iter = 1, itermax
  AxNSolverBeg = nTracer_dclock(.FALSE., .FALSE.)
  CALL SetSource(eigv)
  CALL SetSourceMoment(); CALL SetP1Source()
  phis = 0.0; phim = 0.0
  !$OMP PARALLEL DO SCHEDULE(GUIDED) NUM_THREADS(6)
  DO ixy = 1, nxy
    CALL FlatRayTraceP0(ixy, FALSE)
  ENDDO
  !$OMP END PARALLEL DO
  CALL SetFluxMoment()
  AxNSolverEnd = nTracer_dclock(.FALSE., .FALSE.)
  TimeChk%AxNSolverTime = TimeChk%AxNSolverTime + (AxNSolverEnd - AxNSolverBeg)
  CALL SetBoundaryFlux()
ENDDO

phis = 0.0; phim = 0.0; mklAxial%Jout = 0.0
AxNSolverBeg = nTracer_dclock(.FALSE., .FALSE.)
!$OMP PARALLEL DO SCHEDULE(GUIDED)
DO ixy = 1, nxy
  CALL FlatRayTraceP0(ixy, TRUE)
ENDDO
!$OMP END PARALLEL DO
CALL SetFluxMoment()
AxNSolverEnd = nTracer_dclock(.FALSE., .FALSE.)
TimeChk%AxNSolverTime = TimeChk%AxNSolverTime + (AxNSolverEnd - AxNSolverBeg)
CALL SetBoundaryFlux()

#ifdef FLUX_SHAPE
CALL SetFluxShape()
#else
CALL SetPsiShape()
CALL SetCellFlux()
#endif

! DO iz = 1, nzMOC
!   nNeg(iz) = COUNT(phis(:, iz, :) .LT. 0.0)
! ENDDO
! 
! PRINT *, 'negative flux', mklGeom%myzb, mklGeom%myze
! PRINT *, nNeg(1 : nzMOC)
! PRINT *, 'negative outgoing', mklGeom%myzb, mklGeom%myze, COUNT(mklAxial%PhiAngOut .LT. 0.0)

lFirst = .FALSE.

END SUBROUTINE

!--- Private Routines -----------------------------------------------------------------------------

SUBROUTINE SetSubmesh()

IMPLICIT NONE

INTEGER :: nzCMFD
INTEGER :: iz, izc, izf
INTEGER :: myzb, myze
REAL, POINTER :: hzfm(:)

nzCMFD = mklGeom%nzCMFD
myzb = mklGeom%myzb
myze = mklGeom%myze
hzfm => mklGeom%hzfm

nzMOC = 0

DO iz = 1, nzCMFD
  nDiv(iz) = INT(mklGeom%hzfm(iz) / mklCntl%MOCHeight) + 1
  nzMOC = nzMOC + nDiv(iz)
ENDDO

ALLOCATE(hzMOC(nzMOC))
ALLOCATE(cmRange(myzb : myze, 2), fmRange(nzCMFD, 2))
ALLOCATE(cmMap(nzMOC), fmMap(nzMOC))

izf = 0
DO izc = myzb, myze
  cmRange(izc, 1) = izf + 1
  DO iz = mklGeom%fmRange(izc, 1), mklGeom%fmRange(izc, 2)
    fmRange(iz, 1) = izf + 1
    fmRange(iz, 2) = izf + nDiv(iz)
    fmMap(fmRange(iz, 1) : fmRange(iz, 2)) = iz
    hzMOC(izf + 1 : izf + nDiv(iz)) = hzfm(iz) / nDiv(iz)
    izf = izf + nDiv(iz)
  ENDDO
  cmRange(izc, 2) = izf
  cmMap(cmRange(izc, 1) : cmRange(izc, 2)) = izc
ENDDO

END SUBROUTINE

SUBROUTINE SetSphericalHarmonics()

IMPLICIT NONE

TYPE(mklAngle_Type), POINTER :: Angle(:)
INTEGER :: nPolarAngle
INTEGER :: ipol

Angle => mklAxial%Angle
nPolarAngle = mklGeom%nPolar1D

ALLOCATE(Comp(3, nPolarAngle), mwt(3, nPolarAngle))

DO ipol = 1, nPolarAngle
  Comp(1, ipol) = Angle(ipol)%cosv
  mwt(1, ipol) = Comp(1, ipol) * Angle(ipol)%wtsurf
ENDDO

DO ipol = 1, nPolarAngle
  Comp(2, ipol) = 0.5 * (3.0 * Angle(ipol)%cosv ** 2 - 1.0)
  mwt(2, ipol) = Comp(2, ipol) * Angle(ipol)%wtsurf
ENDDO

DO ipol = 1, nPolarAngle
  Comp(3, ipol) = 0.5 * (5.0 * Angle(ipol)%cosv ** 3 - 3.0 * Angle(ipol)%cosv)
  mwt(3, ipol) = Comp(3, ipol) * Angle(ipol)%wtsurf
ENDDO

END SUBROUTINE

SUBROUTINE SetSourceOperator(PinXS)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)

INTEGER :: ig, igs, ipin, ipin_map, iz, izf
INTEGER :: ng, nxy, nzCMFD
INTEGER, POINTER :: pinMap(:), planeMap(:)

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
pinMap => mklGeom%pinMap
planeMap => mklGeom%planeMap

!$OMP PARALLEL PRIVATE(iz, ipin_map)
!$OMP DO SCHEDULE(GUIDED)
DO ipin = 1, nxy
  ipin_map = pinMap(ipin)
  DO izf = 1, nzCMFD
    iz = planeMap(izf)
    DO ig = 1, ng
      DO igs = 1, ng
        IF (PinXS(ipin_map, iz)%XSs(ig)%ib .GT. igs) CYCLE
        IF (PinXS(ipin_map, iz)%XSs(ig)%ie .LT. igs) CYCLE
        IF (igs .EQ. ig) THEN
          S(igs, ig, izf, ipin) = PinXS(ipin_map, iz)%XSs(ig)%self
        ELSE
          S(igs, ig, izf, ipin) = PinXS(ipin_map, iz)%XSs(ig)%from(igs)
        ENDIF
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(iz, ipin_map)
!$OMP DO SCHEDULE(GUIDED)
DO ipin = 1, nxy
  ipin_map = pinMap(ipin)
  DO izf = 1, nzCMFD
    iz = planeMap(izf)
    DO ig = 1, ng
      F(ig, izf, ipin) = PinXS(ipin_map, iz)%XSnf(ig)
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(iz, ipin_map)
!$OMP DO SCHEDULE(GUIDED)
DO ipin = 1, nxy
  ipin_map = pinMap(ipin)
  DO izf = 1, nzCMFD
    iz = planeMap(izf)
    DO ig = 1, ng
      Chi(ig, izf, ipin) = PinXS(ipin_map, iz)%Chi(ig)
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

#ifdef FLUX_SHAPE

SUBROUTINE SetFlux()

IMPLICIT NONE

INTEGER :: ig, ipin, iz, izf
INTEGER :: ng, nxy, nzCMFD

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD
    DO ig = 1, ng
      IF (mklCMFD%phis(ipin, iz, ig) .LT. 0.0) CYCLE
      DO izf = fmRange(iz, 1), fmRange(iz, 2)
        phis(ig, izf, ipin) = mklCMFD%phis(ipin, iz, ig)
        phis(ig, izf, ipin) = phis(ig, izf, ipin) * phisShape(ig, izf, ipin)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetFluxShape()

IMPLICIT NONE

INTEGER :: ig, ipin, iz, izf
INTEGER :: ng, nxy, nzCMFD
LOGICAL :: lNegative

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

IF (mklCntl%lRefPinFDM) THEN
  !$OMP PARALLEL
  !$OMP DO SCHEDULE(GUIDED)
  DO ipin = 1, nxy
    IF (.NOT. mklGeom%lRefPin(ipin)) CYCLE
    DO izf = 1, nzMOC
      phis(:, izf, ipin) = 1.0
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDIF

mklAxial%phic = 0.0

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      DO ig = 1, ng
        mklAxial%phic(ig, ipin, iz) = mklAxial%phic(ig, ipin, iz) + phis(ig, izf, ipin)
      ENDDO
    ENDDO
    mklAxial%phic(:, ipin, iz) = mklAxial%phic(:, ipin, iz) / nDiv(iz)
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      DO ig = 1, ng
        phisShape(ig, izf, ipin) = phis(ig, izf, ipin) / mklAxial%phic(ig, ipin, iz)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(lNegative)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD
    DO ig = 1, ng
      lNegative = .FALSE.
      DO izf = fmRange(iz, 1), fmRange(iz, 2)
        IF (phis(ig, izf, ipin) .LT. 0.0) lNegative = .TRUE.
      ENDDO
      IF (lNegative) THEN
        DO izf = fmRange(iz, 1), fmRange(iz, 2)
          phisShape(ig, izf, ipin) = 1.0
        ENDDO
      ENDIF
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

#else

SUBROUTINE SetFlux()

IMPLICIT NONE

INTEGER :: ig, ipin, iz, izf
INTEGER :: ng, nxy, nzCMFD
REAL :: fmult

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

!$OMP PARALLEL PRIVATE(fmult)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD
    DO ig = 1, ng
      IF (lFirst) THEN
        DO izf = fmRange(iz, 1), fmRange(iz, 2)
          phis(ig, izf, ipin) = mklCMFD%phis(ipin, iz, ig)
        ENDDO
      ELSE
!        fmult = mklCMFD%phis(ipin, iz, ig) / mklAxial%phic(ig, ipin, iz)
!        DO izf = fmRange(iz, 1), fmRange(iz, 2)
!          phis(ig, izf, ipin) = fmult * phis(ig, izf, ipin)
!        ENDDO
      ENDIF
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetCellFlux()

IMPLICIT NONE

INTEGER :: ig, ipin, iz, izf
INTEGER :: ng, nxy, nzCMFD
LOGICAL :: lNegative

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

mklAxial%phic = 0.0

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      DO ig = 1, ng
        mklAxial%phic(ig, ipin, iz) = mklAxial%phic(ig, ipin, iz) + phis(ig, izf, ipin)
      ENDDO
    ENDDO
    mklAxial%phic(:, ipin, iz) = mklAxial%phic(:, ipin, iz) / nDiv(iz)
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

#endif

SUBROUTINE SetFluxMoment()

IMPLICIT NONE

TYPE(mklAngle_Type), POINTER :: Angle(:)
INTEGER :: ipol, ig, ipin, izf
INTEGER :: ng, nxy, nPolarAngle
REAL :: rsigv

Angle => mklAxial%Angle
ng = mklGeom%ng
nxy = mklGeom%nxy
nPolarAngle = mklGeom%nPolar1D

!$OMP PARALLEL PRIVATE(rsigv) NUM_THREADS(6)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ipin = 1, nxy
  DO izf = 1, nzMOC
    DO ig = 1, ng
      rsigv = 1.0 / hzMOC(izf) / xst(ig, izf, ipin)
      phis(ig, izf, ipin) = phis(ig, izf, ipin) * rsigv + src(ig, izf, ipin)
      phim(1, ig, izf, ipin) = phim(1, ig, izf, ipin) * rsigv + 1.0 / 3.0 * srcm(1, ig, izf, ipin)
!      phim(2, ig, izf, ipin) = phim(2, ig, izf, ipin) * rsigv + 1.0 / 5.0 * srcm(2, ig, izf, ipin)
!      phim(3, ig, izf, ipin) = phim(3, ig, izf, ipin) * rsigv + 1.0 / 7.0 * srcm(3, ig, izf, ipin)
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetCrossSection(PinXS)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)

INTEGER :: ig, ipin, ipin_map, iz, izc, izf
INTEGER :: ng, nxy, nzCMFD
INTEGER, POINTER :: pinMap(:), planeMap(:)

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
pinMap => mklGeom%pinMap
planeMap => mklGeom%planeMap

!$OMP PARALLEL PRIVATE(ipin_map, izc)
!$OMP DO SCHEDULE(GUIDED)
DO ipin = 1, nxy
  ipin_map = pinMap(ipin)
  DO iz = 1, nzCMFD
    izc = planeMap(iz)
    DO ig = 1, ng
      IF (mklGeom%lH2OCell(iz, ipin)) THEN
        S(ig, ig, iz, ipin) = S(ig, ig, iz, ipin) + (PinXS(ipin_map, izc)%XSt(ig) - PinXS(ipin_map, izc)%XStr(ig))
        DO izf = fmRange(iz, 1), fmRange(iz, 2)
          xst(ig, izf, ipin) = pxs(ig, izf, ipin) + PinXS(ipin_map, izc)%XSt(ig)
        ENDDO
      ELSE
        DO izf = fmRange(iz, 1), fmRange(iz, 2)
          xst(ig, izf, ipin) = pxs(ig, izf, ipin) + PinXS(ipin_map, izc)%XStr(ig)
        ENDDO
      ENDIF
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL
  
END SUBROUTINE

#ifdef FLUX_SHAPE

SUBROUTINE SetPsi()

IMPLICIT NONE

INTEGER :: ig, igs, ipin, iz, izf
INTEGER :: ng, nxy, nzCMFD

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

psi = 0.0

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      DO ig = 1, ng
        psi(izf, ipin) = psi(izf, ipin) + F(ig, iz, ipin) * phis(ig, izf, ipin)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

#else

SUBROUTINE SetPsi()

IMPLICIT NONE

INTEGER :: ig, igs, ipin, iz, izf
INTEGER :: ng, nxy, nzCMFD

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      psi(izf, ipin) = mklCMFD%psi(ipin, iz)
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

CALL vdmul(nxy * nzMOC, psiShape, psi, psi)

END SUBROUTINE

SUBROUTINE SetPsiShape()

IMPLICIT NONE

INTEGER :: ig, ipin, iz, izf
INTEGER :: ng, nxy, nzCMFD
REAL :: psic

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

psi = 0.0

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      DO ig = 1, ng
        psi(izf, ipin) = psi(izf, ipin) + F(ig, iz, ipin) * phis(ig, izf, ipin)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

psiShape = 0.0

!$OMP PARALLEL PRIVATE(psic)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD
    psic = 0.0
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      psic = psic + psi(izf, ipin)
    ENDDO
    IF (psic .NE. 0.0) THEN
      DO izf = fmRange(iz, 1), fmRange(iz, 2)
        psiShape(izf, ipin) = psi(izf, ipin) / psic
      ENDDO
    ENDIF
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

#endif

SUBROUTINE SetLeakage(PinXS)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)

TYPE(superPin_Type), POINTER :: Pin(:)
INTEGER :: ig, ibd, ipin, ineighpin, ipin_map, iz, izc, izf
INTEGER :: ng, nxy, nzCMFD
INTEGER, POINTER :: pinMap(:), pinMapRev(:), planeMap(:)
REAL :: Dtil, Dhat, myphi, neighphi
REAL, POINTER :: radLkg(:, :, :)

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
Pin => mklGeom%superPin
pinMap => mklGeom%pinMap
pinMapRev => mklGeom%pinMapRev
planeMap => mklGeom%planeMap

ALLOCATE(radLkg(ng, nxy, 0 : nzCMFD + 1)); radLkg = 0.0

DO iz = 1, nzCMFD
  izc = planeMap(iz)
  !$OMP PARALLEL PRIVATE(ipin_map, ineighpin, myphi, neighphi, Dtil, Dhat)
  !$OMP DO SCHEDULE(GUIDED)
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    DO ig = 1, ng
      myphi = mklCMFD%phis(ipin, iz, ig)
      DO ibd = 1, Pin(ipin)%nNgh
        ineighpin = Pin(ipin_map)%NeighIdx(ibd)
        ineighpin = pinMapRev(ineighpin)
        IF (ineighpin .EQ. VoidCell) THEN
          neighphi = 0.0
        ELSEIF (ineighpin .EQ. RefCell) THEN
          neighphi = myphi
        ELSE
          neighphi = mklCMFD%phis(ineighpin, iz, ig)
        ENDIF
        Dtil = PinXS(ipin_map, izc)%Dtil(ibd, ig)
        Dhat = PinXS(ipin_map, izc)%Dhat(ibd, ig)
        radLkg(ig, ipin, iz) = radLkg(ig, ipin, iz) - Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
      ENDDO
      radLkg(ig, ipin, iz) = radLkg(ig, ipin, iz) * mklGeom%hzfm(iz) / mklGeom%PinVolFm(ipin, iz)
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDDO

CALL InitFastComm()
CALL GetNeighborFast(ng * nxy, radLkg(:, :, 1), radLkg(:, :, nzCMFD + 1), bottom)
CALL GetNeighborFast(ng * nxy, radLkg(:, :, nzCMFD), radLkg(:, :, 0), top)
CALL FinalizeFastComm()

IF (mklGeom%lBottom) THEN
  IF (mklGeom%AxBC(bottom) .EQ. VoidCell) radLkg(:, :, 0) = 0.0
  IF (mklGeom%AxBC(bottom) .EQ. RefCell) CALL dcopy(ng * nxy, radLkg(:, :, 1), 1, radLkg(:, :, 0), 1)
ENDIF

IF (mklGeom%lTop) THEN
  IF (mklGeom%AxBC(top) .EQ. VoidCell) radLkg(:, :, nzCMFD + 1) = 0.0
  IF (mklGeom%AxBC(top) .EQ. RefCell) CALL dcopy(ng * nxy, radLkg(:, :, nzCMFD), 1, radLkg(:, :, nzCMFD + 1), 1)
ENDIF

DO iz = 1, nzCMFD
  CALL LeakageExpansion(radLkg(:, :, iz - 1), radLkg(:, :, iz), radLkg(:, :, iz + 1), iz)
ENDDO

DEALLOCATE(radLkg)

END SUBROUTINE

SUBROUTINE LeakageExpansion(L0, L1, L2, iz)

IMPLICIT NONE

REAL :: L0(:, :), L1(:, :), L2(:, :)
INTEGER :: iz

REAL :: n0(3), n1(3), n2(3)
REAL :: d0, d1, d2
REAL :: h0, h1, h2, dh
REAL :: x0, x1
REAL :: a, b, c
INTEGER :: ig, ipin, izf
INTEGER :: ng, nxy

ng = mklGeom%ng
nxy = mklGeom%nxy

h0 = mklGeom%hzfm(iz - 1)
h1 = mklGeom%hzfm(iz)
h2 = mklGeom%hzfm(iz + 1)
dh = h1 / nDiv(iz)

n0(1) = h1 ** 3 + 2.0 * h1 ** 2 * h2 + h1 * h2 ** 2
n0(2) = 2.0 * h0 ** 2 * h1 + 3.0 * h0 * h1 ** 2 + h0 ** 2 * h2 + 3.0 * h0 * h1 * h2 + h0 * h2 ** 2
n0(3) = - h0 ** 2 * h1 - h0 * h1 ** 2

n1(1) = 2.0 * h1 ** 2 + 3.0 * h1 * h2 + h2 ** 2
n1(2) = h0 ** 2 - 3.0 * h1 ** 2 - 3.0 * h1 * h2 - h2 ** 2
n1(3) = - h0 ** 2 + h1 ** 2

n2(1) = h1 + h2
n2(2) = - h0 - 2.0 * h1 - h2
n2(3) = h0 + h1

d0 = (h1 + h2) * (h0 ** 2 + 2.0 * h0 * h1 + h1 ** 2 + h0 * h2 + h1 * h2)
d1 = (h0 + h1) * (h1 + h2) * (h0 + h1 + h2)
d2 = (h0 + h1) * (h1 + h2) * (h0 + h1 + h2)

!$OMP PARALLEL PRIVATE(a, b, c, x0, x1)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO ig = 1, ng
    a = 3.0 * (L0(ig, ipin) * n2(1) + L1(ig, ipin) * n2(2) + L2(ig, ipin) * n2(3)) / d2
    b = - 2.0 * (L0(ig, ipin) * n1(1) + L1(ig, ipin) * n1(2) + L2(ig, ipin) * n1(3)) / d1
    c = (L0(ig, ipin) * n0(1) + L1(ig, ipin) * n0(2) + L2(ig, ipin) * n0(3)) / d0
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      x0 = (izf - fmRange(iz, 1)) * dh; x1 = x0 + dh
      lkg(ig, izf, ipin) = Integral(a, b, c, x0, x1) / dh
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

CONTAINS

FUNCTION Integral(a, b, c, x0, x1) RESULT(val)

IMPLICIT NONE

REAL :: a, b, c, x0, x1
REAL :: val

val = a * (x1 ** 3 - x0 ** 3) / 3.0 + b * (x1 ** 2 - x0 ** 2) / 2.0 + c * (x1 - x0)

END FUNCTION

END SUBROUTINE

SUBROUTINE SetPseudoAbsorption()

IMPLICIT NONE

INTEGER :: ig, ipin, iz, izf
INTEGER :: ng, nxy, nzCMFD

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

pxs = 0.0

!$OMP PARALLEL PRIVATE(iz)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO izf = 1, nzMOC
    iz = fmMap(izf); IF (.NOT. mklGeom%lRefCell(iz, ipin)) CYCLE
    DO ig = 1, ng
      IF (lkg(ig, izf, ipin) .GT. 0.0 .AND. phis(ig, izf, ipin) .GT. 0.0) THEN
        pxs(ig, izf, ipin) = lkg(ig, izf, ipin) / phis(ig, izf, ipin)
        lkg(ig, izf, ipin) = 0.0
      ENDIF
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD
    IF (.NOT. mklGeom%lRefCell(iz, ipin)) CYCLE
    IF (mklGeom%lH2OCell(iz, ipin)) CYCLE
    DO ig = 1, ng
      IF (S(ig, ig, iz, ipin) .LT. 0.0) THEN
        DO izf = fmRange(iz, 1), fmRange(iz, 2)
          pxs(ig, izf, ipin) = pxs(ig, izf, ipin) - S(ig, ig, iz, ipin)
        ENDDO
        S(ig, ig, iz, ipin) = 0.0
      ENDIF
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetConstant(PinXS)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)

TYPE(mklAngle_Type), POINTER :: Angle(:)
INTEGER :: ipol, ig, ipin, iz, izf
INTEGER :: ng, nxy, nPolarAngle
INTEGER, POINTER :: pinMap(:)
REAL :: tau

ng = mklGeom%ng
nxy = mklGeom%nxy
nPolarAngle = mklGeom%nPolar1D
pinMap => mklGeom%pinMap
Angle => mklAxial%Angle

!$OMP PARALLEL PRIVATE(tau)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(4)
DO ipin = 1, nxy
  DO izf = 1, nzMOC
    DO ig = 1, ng
      DO ipol = 1, nPolarAngle
        tau = xst(ig, izf, ipin) * hzMOC(izf) * Angle(ipol)%rcosv
        wtExp(ipol, ig, izf, ipin) = 1.0 - EXP(- tau)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetSourceMoment()

INTEGER :: ig, igs, ipin, iz, izf
INTEGER :: ng, nxy, nzCMFD
INTEGER :: gb, ge

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

srcm = 0.0

!$OMP PARALLEL PRIVATE(gb, ge) NUM_THREADS(6)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD
    IF (.NOT. mklGeom%lH2OCell(iz, ipin)) CYCLE
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      DO ig = 1, ng
        gb = mklGeom%InScatRange(1, ig)
        ge = mklGeom%InScatRange(2, ig)
        DO igs = gb, ge
          srcm(1, ig, izf, ipin) = srcm(1, ig, izf, ipin) + mklAxial%SmP1(igs, ig, iz, ipin) * phim(1, igs, izf, ipin)
!          srcm(2, ig, izf, ipin) = srcm(2, ig, izf, ipin) + mklAxial%SmP2(igs, ig, iz, ipin) * phim(2, igs, izf, ipin)
!          srcm(3, ig, izf, ipin) = srcm(3, ig, izf, ipin) + mklAxial%SmP3(igs, ig, iz, ipin) * phim(3, igs, izf, ipin)
        ENDDO
        srcm(1, ig, izf, ipin) = 3.0 * srcm(1, ig, izf, ipin) / xst(ig, izf, ipin)
!        srcm(2, ig, izf, ipin) = 5.0 *srcm(2, ig, izf, ipin) / xst(ig, izf, ipin)
!        srcm(3, ig, izf, ipin) = 7.0 * srcm(3, ig, izf, ipin) / xst(ig, izf, ipin)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetSource(eigv)

IMPLICIT NONE

REAL :: eigv, reigv

TYPE(mklAngle_Type), POINTER :: Angle(:)
INTEGER :: ipol, ig, igs, ipin, iz, izf
INTEGER :: ng, nxy, nzCMFD, nPolarAngle
INTEGER :: gb, ge
REAL :: tau

Angle => mklAxial%Angle
ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
nPolarAngle = mklGeom%nPolar1D

reigv = 1.0 / eigv

src = 0.0

!$OMP PARALLEL PRIVATE(gb, ge) NUM_THREADS(6)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      DO ig = 1, ng
        gb = mklGeom%InScatRange(1, ig)
        ge = mklGeom%InScatRange(2, ig)
        DO igs = gb, ge
          src(ig, izf, ipin) = src(ig, izf, ipin) + S(igs, ig, iz, ipin) * phis(igs, izf, ipin)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL 
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      DO ig = 1, ng
        src(ig, izf, ipin) = src(ig, izf, ipin) + reigv * Chi(ig, iz, ipin) * psi(izf, ipin)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

CALL vdsub(ng * nxy * nzMOC, src, lkg, src)
CALL vddiv(ng * nxy * nzMOC, src, xst, src)

END SUBROUTINE

SUBROUTINE SetP1Source()

IMPLICIT NONE

INTEGER :: ipol, ig, ipin, izf
INTEGER :: ng, nxy, nzCMFD, nPolarAngle

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
nPolarAngle = mklGeom%nPolar1D

!$OMP PARALLEL NUM_THREADS(6)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(4)
DO ipin = 1, nxy
  DO izf = 1, nzMOC
    DO ig = 1, ng
      DO ipol = 1, nPolarAngle
        srca(1, ipol, ig, izf, ipin) = src(ig, izf, ipin) + Comp(1, ipol) * srcm(1, ig, izf, ipin)
!        srca(1, ipol, ig, izf, ipin) = srca(1, ipol, ig, izf, ipin) + Comp(2, ipol) * srcm(2, ig, izf, ipin)
!        srca(1, ipol, ig, izf, ipin) = srca(1, ipol, ig, izf, ipin) + Comp(3, ipol) * srcm(3, ig, izf, ipin)
        srca(2, ipol, ig, izf, ipin) = src(ig, izf, ipin) - Comp(1, ipol) * srcm(1, ig, izf, ipin)
!        srca(2, ipol, ig, izf, ipin) = srca(2, ipol, ig, izf, ipin) + Comp(2, ipol) * srcm(2, ig, izf, ipin)
!        srca(2, ipol, ig, izf, ipin) = srca(2, ipol, ig, izf, ipin) - Comp(3, ipol) * srcm(3, ig, izf, ipin)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetBoundaryFlux()

IMPLICIT NONE

INTEGER :: n, ng, nxy, nPolarAngle
INTEGER :: ipol, ig, ipin, iz
REAL, POINTER :: PhiAngIn(:, :, :, :)

ng = mklGeom%ng
nxy = mklGeom%nxy
nPolarAngle = mklGeom%nPolar1D
n = nPolarAngle * ng * nxy

ALLOCATE(PhiAngIn(nPolarAngle, ng, nxy, 2))

CALL InitFastComm()
CALL GetNeighborFast(n, mklAxial%PhiAngOut(:, :, :, bottom), PhiAngIn(:, :, :, top), bottom)
CALL GetNeighborFast(n, mklAxial%PhiAngOut(:, :, :, top), PhiAngIn(:, :, :, bottom), top)
CALL FinalizeFastComm()

IF (mklGeom%lBottom) THEN
  IF (mklGeom%AxBC(bottom) .EQ. VoidCell) THEN
    PhiAngIn(:, :, :, bottom) = 0.0
  ELSEIF (mklGeom%AxBC(bottom) .EQ. RefCell) THEN
    CALL dcopy(n, mklAxial%PhiAngOut(:, :, :, bottom), 1, PhiAngIn(:, :, :, bottom), 1)
  ENDIF
ENDIF

IF (mklGeom%lTop) THEN
  IF (mklGeom%AxBC(top) .EQ. VoidCell) THEN
    PhiAngIn(:, :, :, top) = 0.0
  ELSEIF (mklGeom%AxBC(top) .EQ. RefCell) THEN
    CALL dcopy(n, mklAxial%PhiAngOut(:, :, :, top), 1, PhiAngIn(:, :, :, top), 1)
  ENDIF
ENDIF

CALL dcopy(n * 2, PhiAngIn, 1, mklAxial%PhiAngIn, 1)

DEALLOCATE(PhiAngIn)

END SUBROUTINE

SUBROUTINE FlatRayTraceP0(ipin, lJout)

IMPLICIT NONE

TYPE(mklAngle_Type), POINTER :: Angle(:)
REAL, POINTER :: Jout(:, :, :, :), PhiAngIn(:, :, :, :), PhiAngOut(:, :, :, :)
REAL :: del_phi, track_phi(mklGeom%nPolar1D, mklGeom%ng)
INTEGER :: ig, ipin, ipol, iz, izf
INTEGER :: ng, nxy, nzCMFD, nPolarAngle
LOGICAL :: lJout

IF (mklCntl%lRefPinFDM) THEN
  IF (mklGeom%lRefPin(ipin)) RETURN
ENDIF

PhiAngIn => mklAxial%PhiAngIn
PhiAngOut => mklAxial%PhiAngOut
Jout => mklAxial%Jout(:, :, :, :, ipin)
Angle => mklAxial%Angle

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
nPolarAngle = mklGeom%nPolar1D

!--- Upward Sweep

track_phi = PhiAngIn(:, :, ipin, bottom)

DO iz = 1, nzCMFD
  IF (lJout) THEN
    DO ig = 1, ng
      DO ipol = 1, nPolarAngle
        Jout(in, ig, bottom, iz) = Jout(in, ig, bottom, iz) + Angle(ipol)%wtsurf * track_phi(ipol, ig)
      ENDDO
    ENDDO
  ENDIF
  IF (mklGeom%lH2OCell(iz, ipin)) THEN
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      DO ig = 1, ng
        DO ipol = 1, nPolarAngle
          del_phi = (srca(1, ipol, ig, izf, ipin) - track_phi(ipol, ig)) * wtExp(ipol, ig, izf, ipin)
          phis(ig, izf, ipin) = phis(ig, izf, ipin) - Angle(ipol)%wtsurf * del_phi
          phim(1, ig, izf, ipin) = phim(1, ig, izf, ipin) - mwt(1, ipol) * del_phi
!          phim(2, ig, izf, ipin) = phim(2, ig, izf, ipin) - mwt(2, ipol) * del_phi
!          phim(3, ig, izf, ipin) = phim(3, ig, izf, ipin) - mwt(3, ipol) * del_phi
          track_phi(ipol, ig) = track_phi(ipol, ig) + del_phi
        ENDDO
      ENDDO
    ENDDO
  ELSE
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      DO ig = 1, ng
        DO ipol = 1, nPolarAngle
          del_phi = (srca(1, ipol, ig, izf, ipin) - track_phi(ipol, ig)) * wtExp(ipol, ig, izf, ipin)
          phis(ig, izf, ipin) = phis(ig, izf, ipin) - Angle(ipol)%wtsurf * del_phi
          track_phi(ipol, ig) = track_phi(ipol, ig) + del_phi
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  IF (lJout) THEN
    DO ig = 1, ng
      DO ipol = 1, nPolarAngle
        Jout(out, ig, top, iz) = Jout(out, ig, top, iz) + Angle(ipol)%wtsurf * track_phi(ipol, ig)
      ENDDO
    ENDDO
  ENDIF
ENDDO

PhiAngOut(:, :, ipin, top) = track_phi

!--- Downward Sweep

track_phi = PhiAngIn(:, :, ipin, top)

DO iz = nzCMFD, 1, -1
  IF (lJout) THEN
    DO ig = 1, ng
      DO ipol = 1, nPolarAngle
        Jout(in, ig, top, iz) = Jout(in, ig, top, iz) + Angle(ipol)%wtsurf * track_phi(ipol, ig)
      ENDDO
    ENDDO
  ENDIF
  IF (mklGeom%lH2OCell(iz, ipin)) THEN
    DO izf = fmRange(iz, 2), fmRange(iz, 1), -1
      DO ig = 1, ng
        DO ipol = 1, nPolarAngle
          del_phi = (srca(2, ipol, ig, izf, ipin) - track_phi(ipol, ig)) * wtExp(ipol, ig, izf, ipin)
          phis(ig, izf, ipin) = phis(ig, izf, ipin) - Angle(ipol)%wtsurf * del_phi
          phim(1, ig, izf, ipin) = phim(1, ig, izf, ipin) + mwt(1, ipol) * del_phi
!          phim(2, ig, izf, ipin) = phim(2, ig, izf, ipin) - mwt(2, ipol) * del_phi
!          phim(3, ig, izf, ipin) = phim(3, ig, izf, ipin) + mwt(3, ipol) * del_phi
          track_phi(ipol, ig) = track_phi(ipol, ig) + del_phi
        ENDDO
      ENDDO
    ENDDO
  ELSE
    DO izf = fmRange(iz, 2), fmRange(iz, 1), -1
      DO ig = 1, ng
        DO ipol = 1, nPolarAngle
          del_phi = (srca(2, ipol, ig, izf, ipin) - track_phi(ipol, ig)) * wtExp(ipol, ig, izf, ipin)
          phis(ig, izf, ipin) = phis(ig, izf, ipin) - Angle(ipol)%wtsurf * del_phi
          track_phi(ipol, ig) = track_phi(ipol, ig) + del_phi
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  IF (lJout) THEN
    DO ig = 1, ng
      DO ipol = 1, nPolarAngle
        Jout(out, ig, bottom, iz) = Jout(out, ig, bottom, iz) + Angle(ipol)%wtsurf * track_phi(ipol, ig)
      ENDDO
    ENDDO
  ENDIF
ENDDO

PhiAngOut(:, :, ipin, bottom) = track_phi

END SUBROUTINE

END MODULE

#endif