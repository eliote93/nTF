#include <defines.h>
!--- CNJ Edit : Solver Routine for Gamma Axial Calculation
#ifdef __INTEL_MKL
#ifdef __GAMMA_TRANSPORT

MODULE GAMMA_FLATMOC

USE MKL_3D
IMPLICIT NONE

REAL, POINTER, PRIVATE :: phis(:, :, :), phim(:, :, :, :)
REAL, POINTER, PRIVATE :: prod(:, :, :), src(:, :, :), srca(:, :, :, :, :), srcm(:, :, :, :), lkg(:, :, :)
REAL, POINTER, PRIVATE :: xst(:, :, :), pxs(:, :, :)
REAL, POINTER, PRIVATE :: S(:, :, :, :), P(:, :, :, :)

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

SUBROUTINE AllocFlatMOC(Axial)

IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial

INTEGER :: ng, ngg, nxy, nzCMFD, nPolarAngle

ng = mklGeom%ng
ngg = mklGeom%ngg
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
nPolarAngle = mklGeom%nPolar1D

CALL SetSphericalHarmonics(Axial)
CALL SetSubmesh()

ALLOCATE(phis(ngg, nzMOC, nxy))
ALLOCATE(phim(3, ngg, nzMOC, nxy)); phim = 0.0
ALLOCATE(prod(ngg, nzMOC, nxy))
ALLOCATE(src(ngg, nzMOC, nxy))
ALLOCATE(srca(2, nPolarAngle, ngg, nzMOC, nxy)); srca = 0.0
ALLOCATE(srcm(3, ngg, nzMOC, nxy)); srcm = 0.0
ALLOCATE(lkg(ngg, nzMOC, nxy))
ALLOCATE(xst(ngg, nzMOC, nxy))
ALLOCATE(pxs(ngg, nzMOC, nxy))

ALLOCATE(wtExp(nPolarAngle, ngg, nzMOC, nxy))

ALLOCATE(S(ngg, ngg, nzCMFD, nxy)); S = 0.0
ALLOCATE(P(ng, ngg, nzCMFD, nxy)); P = 0.0

END SUBROUTINE

SUBROUTINE FlatMOCDriver(GammaCMFD, Axial, GammaAxial)
USE PE_MOD,         ONLY : PE
USE TIMER,          ONLY : nTracer_dclock,      TimeChk
IMPLICIT NONE

TYPE(mklCMFD_Type) :: GammaCMFD
TYPE(mklAxial_Type) :: Axial, GammaAxial

TYPE(GPinXS_Type), POINTER :: GPinXS(:, :)
INTEGER :: i, ig, ixy, iz, izf, nxy, nzCMFD, nNeg(1000)
INTEGER :: ierr, iter, itermax = 10
REAL :: AxNSolverBeg, AxNSolverEnd

GPinXS => GammaCMFD%GPinXS
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

CALL SetSourceOperator(GPinXS)
CALL SetFlux(GammaCMFD, GammaAxial)
CALL SetProd(Axial)
CALL SetLeakage(GammaCMFD, GammaAxial, GPinXS)
CALL SetPseudoAbsorption(GammaAxial)
CALL SetCrossSection(GammaAxial, GPinXS)
CALL SetConstant(GammaAxial)

!PRINT *, 'negative feed flux', mklGeom%myzb, mklGeom%myze, COUNT(phis .LT. 0.0)
!PRINT *, 'negative incoming', mklGeom%myzb, mklGeom%myze, COUNT(GammaAxial%PhiAngIn .LT. 0.0)

DO iter = 1, itermax
  AxNSolverBeg = nTracer_dclock(.FALSE., .FALSE.)
  CALL SetSource(GammaCMFD, GammaAxial)
  CALL SetSourceMoment(GammaCMFD, GammaAxial); CALL SetP1Source(GammaAxial)
  phis = 0.0; phim = 0.0
  !$OMP PARALLEL DO SCHEDULE(GUIDED)
  DO ixy = 1, nxy
    CALL FlatRayTraceP0(GammaAxial, ixy, FALSE)
  ENDDO
  !$OMP END PARALLEL DO
  CALL SetFluxMoment(GammaAxial)
  AxNSolverEnd = nTracer_dclock(.FALSE., .FALSE.)
  TimeChk%AxNSolverTime = TimeChk%AxNSolverTime + (AxNSolverEnd - AxNSolverBeg)
  CALL SetBoundaryFlux(GammaAxial)
ENDDO

phis = 0.0; phim = 0.0; GammaAxial%Jout = 0.0
AxNSolverBeg = nTracer_dclock(.FALSE., .FALSE.)
!$OMP PARALLEL DO SCHEDULE(GUIDED)
DO ixy = 1, nxy
  CALL FlatRayTraceP0(GammaAxial, ixy, TRUE)
ENDDO
!$OMP END PARALLEL DO
CALL SetFluxMoment(GammaAxial)
AxNSolverEnd = nTracer_dclock(.FALSE., .FALSE.)
TimeChk%AxNSolverTime = TimeChk%AxNSolverTime + (AxNSolverEnd - AxNSolverBeg)
CALL SetBoundaryFlux(GammaAxial)
CALL SetCellFlux(GammaAxial)

DO iz = 1, nzMOC
  nNeg(iz) = COUNT(phis(:, iz, :) .LT. 0.0)
ENDDO

!PRINT *, 'negative flux', mklGeom%myzb, mklGeom%myze
!PRINT *, nNeg(1 : nzMOC)
!PRINT *, 'negative outgoing', mklGeom%myzb, mklGeom%myze, COUNT(GammaAxial%PhiAngOut .LT. 0.0)

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

SUBROUTINE SetSphericalHarmonics(Axial)

IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial

TYPE(mklAngle_Type), POINTER :: Angle(:)
INTEGER :: nPolarAngle
INTEGER :: ipol

Angle => Axial%Angle
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

SUBROUTINE SetSourceOperator(GPinXS)

IMPLICIT NONE

TYPE(GPinXS_Type), POINTER :: GPinXS(:, :)

INTEGER :: ig, igg, igs, ipin, ipin_map, iz, izf
INTEGER :: ng, ngg, nxy, nzCMFD
INTEGER, POINTER :: pinMap(:), planeMap(:)

ng = mklGeom%ng
ngg = mklGeom%ngg
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
    DO igg = 1, ngg
      DO igs = 1, ngg
        IF (GPinXS(ipin_map, iz)%XSs(igg)%ib .GT. igs) CYCLE
        IF (GPinXS(ipin_map, iz)%XSs(igg)%ie .LT. igs) CYCLE
        IF (igs .EQ. igg) THEN
          S(igs, igg, izf, ipin) = GPinXS(ipin_map, iz)%XSs(igg)%self
        ELSE
          S(igs, igg, izf, ipin) = GPinXS(ipin_map, iz)%XSs(igg)%from(igs)
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
    P(:, :, izf, ipin) = GPinXS(ipin_map, iz)%XSP
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetFlux(CMFD, Axial)

IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
TYPE(mklAxial_Type) :: Axial

INTEGER :: ig, ipin, iz, izf
INTEGER :: ng, nxy, nzCMFD
REAL :: fmult

ng = Axial%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

!$OMP PARALLEL PRIVATE(fmult)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD
    DO ig = 1, ng
      IF (lFirst) THEN
        DO izf = fmRange(iz, 1), fmRange(iz, 2)
          phis(ig, izf, ipin) = CMFD%phis(ipin, iz, ig)
        ENDDO
      ELSE
        fmult = CMFD%phis(ipin, iz, ig) / Axial%phic(ig, ipin, iz)
        DO izf = fmRange(iz, 1), fmRange(iz, 2)
          phis(ig, izf, ipin) = fmult * phis(ig, izf, ipin)
        ENDDO
      ENDIF
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetCellFlux(Axial)

IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial

INTEGER :: ig, ipin, iz, izf
INTEGER :: ng, nxy, nzCMFD
LOGICAL :: lNegative

ng = Axial%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

Axial%phic = 0.0

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      DO ig = 1, ng
        Axial%phic(ig, ipin, iz) = Axial%phic(ig, ipin, iz) + phis(ig, izf, ipin)
      ENDDO
    ENDDO
    Axial%phic(:, ipin, iz) = Axial%phic(:, ipin, iz) / nDiv(iz)
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetFluxMoment(Axial)

IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial

TYPE(mklAngle_Type), POINTER :: Angle(:)
INTEGER :: ipol, ig, ipin, izf
INTEGER :: ng, nxy, nPolarAngle
REAL :: rsigv

Angle => Axial%Angle
ng = Axial%ng
nxy = mklGeom%nxy
nPolarAngle = mklGeom%nPolar1D

!$OMP PARALLEL PRIVATE(rsigv)
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

SUBROUTINE SetCrossSection(Axial, GPinXS)

IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial
TYPE(GPinXS_Type), POINTER :: GPinXS(:, :)

INTEGER :: ig, ipin, ipin_map, iz, izc, izf
INTEGER :: ng, nxy, nzCMFD
INTEGER, POINTER :: pinMap(:), planeMap(:)

ng = Axial%ng
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
        S(ig, ig, iz, ipin) = S(ig, ig, iz, ipin) + (GPinXS(ipin_map, izc)%XSt(ig) - GPinXS(ipin_map, izc)%XStr(ig))
        DO izf = fmRange(iz, 1), fmRange(iz, 2)
          xst(ig, izf, ipin) = pxs(ig, izf, ipin) + GPinXS(ipin_map, izc)%XSt(ig)
        ENDDO
      ELSE
        DO izf = fmRange(iz, 1), fmRange(iz, 2)
          xst(ig, izf, ipin) = pxs(ig, izf, ipin) + GPinXS(ipin_map, izc)%XStr(ig)
        ENDDO
      ENDIF
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetProd(Axial)

IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial

INTEGER :: ig, igg, ipin, iz, izf
INTEGER :: ng, ngg, nxy, nzCMFD

ng = mklGeom%ng
ngg = mklGeom%ngg
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

prod = 0.0

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      DO igg = 1, ngg
        DO ig = 1, ng
          prod(igg, izf, ipin) = prod(igg, izf, ipin) + P(ig, igg, iz, ipin) * phis(ig, izf, ipin)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetLeakage(CMFD, Axial, GPinXS)

IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
TYPE(mklAxial_Type) :: Axial
TYPE(GPinXS_Type), POINTER :: GPinXS(:, :)

TYPE(superPin_Type), POINTER :: Pin(:)
INTEGER :: ig, ibd, ipin, ineighpin, ipin_map, iz, izc, izf
INTEGER :: ng, nxy, nzCMFD
INTEGER, POINTER :: pinMap(:), pinMapRev(:), planeMap(:)
REAL :: Dtil, Dhat, myphi, neighphi
REAL, POINTER :: radLkg(:, :, :)

ng = Axial%ng
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
      myphi = CMFD%phis(ipin, iz, ig)
      DO ibd = 1, 4
        ineighpin = Pin(ipin_map)%NeighIdx(ibd)
        ineighpin = pinMapRev(ineighpin)
        IF (ineighpin .EQ. VoidCell) THEN
          neighphi = 0.0
        ELSEIF (ineighpin .EQ. RefCell) THEN
          neighphi = myphi
        ELSE
          neighphi = CMFD%phis(ineighpin, iz, ig)
        ENDIF
        Dtil = GPinXS(ipin_map, izc)%Dtil(ibd, ig)
        Dhat = GPinXS(ipin_map, izc)%Dhat(ibd, ig)
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
  CALL LeakageExpansion(radLkg(:, :, iz - 1), radLkg(:, :, iz), radLkg(:, :, iz + 1), iz, ng, nxy)
ENDDO

DEALLOCATE(radLkg)

END SUBROUTINE

SUBROUTINE LeakageExpansion(L0, L1, L2, iz, ng, nxy)

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
      lkg(ig, izf, ipin) = L1(ig, ipin) ! / dh
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

SUBROUTINE SetPseudoAbsorption(Axial)

IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial

INTEGER :: ig, ipin, iz, izf
INTEGER :: ng, nxy, nzCMFD

ng = Axial%ng
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

SUBROUTINE SetConstant(Axial)

IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial

TYPE(mklAngle_Type), POINTER :: Angle(:)
INTEGER :: ipol, ig, ipin, iz, izf
INTEGER :: ng, nxy, nPolarAngle
INTEGER, POINTER :: pinMap(:)
REAL :: tau

ng = Axial%ng
nxy = mklGeom%nxy
nPolarAngle = mklGeom%nPolar1D
pinMap => mklGeom%pinMap
Angle => Axial%Angle

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

SUBROUTINE SetSourceMoment(CMFD, Axial)

IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
TYPE(mklAxial_Type) :: Axial

INTEGER :: ig, igs, ipin, iz, izf
INTEGER :: ng, nxy, nzCMFD
INTEGER :: gb, ge

ng = Axial%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

srcm = 0.0

!$OMP PARALLEL PRIVATE(gb, ge)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD
    IF (.NOT. mklGeom%lH2OCell(iz, ipin)) CYCLE
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      DO ig = 1, ng
        gb = CMFD%InScatRange(1, ig)
        ge = CMFD%InScatRange(2, ig)
        DO igs = gb, ge
          srcm(1, ig, izf, ipin) = srcm(1, ig, izf, ipin) + Axial%SmP1(igs, ig, iz, ipin) * phim(1, igs, izf, ipin)
!          srcm(2, ig, izf, ipin) = srcm(2, ig, izf, ipin) + Axial%SmP2(igs, ig, iz, ipin) * phim(2, igs, izf, ipin)
!          srcm(3, ig, izf, ipin) = srcm(3, ig, izf, ipin) + Axial%SmP3(igs, ig, iz, ipin) * phim(3, igs, izf, ipin)
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

SUBROUTINE SetSource(CMFD, Axial)

IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
TYPE(mklAxial_Type) :: Axial

INTEGER :: ipol, ig, igs, ipin, iz, izf
INTEGER :: ng, nxy, nzCMFD, nPolarAngle
INTEGER :: gb, ge
REAL :: tau

ng = Axial%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
nPolarAngle = mklGeom%nPolar1D

src = 0.0

!$OMP PARALLEL PRIVATE(gb, ge)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      DO ig = 1, ng
        gb = CMFD%InScatRange(1, ig)
        ge = CMFD%InScatRange(2, ig)
        DO igs = gb, ge
          src(ig, izf, ipin) = src(ig, izf, ipin) + S(igs, ig, iz, ipin) * phis(igs, izf, ipin)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

CALL vdadd(ng * nxy * nzMOC, src, prod, src)
CALL vdsub(ng * nxy * nzMOC, src, lkg, src)
CALL vddiv(ng * nxy * nzMOC, src, xst, src)

END SUBROUTINE

SUBROUTINE SetP1Source(Axial)

IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial

INTEGER :: ipol, ig, ipin, izf
INTEGER :: ng, nxy, nzCMFD, nPolarAngle

ng = Axial%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
nPolarAngle = mklGeom%nPolar1D

!$OMP PARALLEL
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

SUBROUTINE SetBoundaryFlux(Axial)

IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial

INTEGER :: n, ng, nxy, nPolarAngle
INTEGER :: ipol, ig, ipin, iz
REAL, POINTER :: PhiAngIn(:, :, :, :)

ng = Axial%ng
nxy = mklGeom%nxy
nPolarAngle = mklGeom%nPolar1D
n = nPolarAngle * ng * nxy

ALLOCATE(PhiAngIn(nPolarAngle, ng, nxy, 2))

CALL InitFastComm()
CALL GetNeighborFast(n, Axial%PhiAngOut(:, :, :, bottom), PhiAngIn(:, :, :, top), bottom)
CALL GetNeighborFast(n, Axial%PhiAngOut(:, :, :, top), PhiAngIn(:, :, :, bottom), top)
CALL FinalizeFastComm()

IF (mklGeom%lBottom) THEN
  IF (mklGeom%AxBC(bottom) .EQ. VoidCell) THEN
    PhiAngIn(:, :, :, bottom) = 0.0
  ELSEIF (mklGeom%AxBC(bottom) .EQ. RefCell) THEN
    CALL dcopy(n, Axial%PhiAngOut(:, :, :, bottom), 1, PhiAngIn(:, :, :, bottom), 1)
  ENDIF
ENDIF

IF (mklGeom%lTop) THEN
  IF (mklGeom%AxBC(top) .EQ. VoidCell) THEN
    PhiAngIn(:, :, :, top) = 0.0
  ELSEIF (mklGeom%AxBC(top) .EQ. RefCell) THEN
    CALL dcopy(n, Axial%PhiAngOut(:, :, :, top), 1, PhiAngIn(:, :, :, top), 1)
  ENDIF
ENDIF

CALL dcopy(n * 2, PhiAngIn, 1, Axial%PhiAngIn, 1)

DEALLOCATE(PhiAngIn)

END SUBROUTINE

SUBROUTINE FlatRayTraceP0(Axial, ipin, lJout)

IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial

TYPE(mklAngle_Type), POINTER :: Angle(:)
REAL, POINTER :: Jout(:, :, :, :), PhiAngIn(:, :, :, :), PhiAngOut(:, :, :, :)
REAL :: del_phi, track_phi(mklGeom%nPolar1D, Axial%ng)
INTEGER :: ig, ipin, ipol, iz, izf
INTEGER :: ng, nxy, nzCMFD, nPolarAngle
LOGICAL :: lJout

IF (mklCntl%lRefPinFDM) THEN
  IF (mklGeom%lRefPin(ipin)) RETURN
ENDIF

PhiAngIn => Axial%PhiAngIn
PhiAngOut => Axial%PhiAngOut
Jout => Axial%Jout(:, :, :, :, ipin)
Angle => Axial%Angle

ng = Axial%ng
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
#endif
