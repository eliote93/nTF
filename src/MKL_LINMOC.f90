#include <defines.h>
!--- CNJ Edit : 1D Axial Linear MOC Modules with Intel MKL
#ifdef __INTEL_MKL

MODULE MKL_LINMOC

USE MKL_3D
IMPLICIT NONE

REAL, POINTER, PRIVATE :: Qc(:, :, :, :, :), Qz(:, :, :, :, :)
REAL, POINTER, PRIVATE :: E1(:, :, :, :), E3(:, :, :, :), R1(:, :, :, :), R3(:, :, :, :)
REAL, POINTER, PRIVATE :: Comp(:, :), mwt(:, :)

REAL, POINTER, PRIVATE :: phisShape(:, :, :), phisCoeff(:, :, :, :), phimCoeff(:, :, :, :)
REAL, POINTER, PRIVATE :: src(:, :, :, :), srcm(:, :, :, :), psi(:, :, :), lkg(:, :, :)
REAL, POINTER, PRIVATE :: xst(:, :, :), pxs(:, :, :)
REAL, POINTER, PRIVATE :: phi(:, :, :, :, :), phimx(:, :, :)
REAL, POINTER, PRIVATE :: S(:, :, :, :), F(:, :, :), Chi(:, :, :)

REAL, POINTER, PRIVATE :: hzMOC(:)
INTEGER, POINTER, PRIVATE :: cmRange(:, :), fmRange(:, :)
INTEGER, POINTER, PRIVATE :: cmMap(:), fmMap(:)
INTEGER, PRIVATE :: nzMOC, nDiv(100)

PRIVATE
PUBLIC :: AllocLinearMOC, LinearMOCDriver

CONTAINS

!--- Public Routines ------------------------------------------------------------------------------

SUBROUTINE AllocLinearMOC()

IMPLICIT NONE

INTEGER :: ng, nxy, nzCMFD, nPolarAngle

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
nPolarAngle = mklGeom%nPolar1D

CALL SetSphericalHarmonics()
CALL SetSubmesh()

ALLOCATE(Qc(nPolarAngle, ng, nzMOC, 2, nxy))
ALLOCATE(Qz(nPolarAngle, ng, nzMOC, 2, nxy))
ALLOCATE(E1(nPolarAngle, ng, nzMOC, nxy))
ALLOCATE(E3(nPolarAngle, ng, nzMOC, nxy))
ALLOCATE(R1(nPolarAngle, ng, nzMOC, nxy))
ALLOCATE(R3(nPolarAngle, ng, nzMOC, nxy))

ALLOCATE(phisShape(ng, nzMOC, nxy)); phisShape = 1.0
ALLOCATE(phisCoeff(ng, nzMOC, nxy, 0 : 1)); phisCoeff = 0.0
ALLOCATE(phimCoeff(1, ng, nzMOC, nxy)); phimCoeff = 0.0
ALLOCATE(phi(0 : 1, ng, nzMOC, 2, nxy))
ALLOCATE(phimx(ng, nzMOC, nxy))

ALLOCATE(psi(nzMOC, nxy, 0 : 1))
ALLOCATE(src(ng, nzMOC, nxy, 0 : 1))
ALLOCATE(srcm(1, ng, nzMOC, nxy)); srcm = 0.0
ALLOCATE(xst(ng, nzMOC, nxy))
ALLOCATE(pxs(ng, nzMOC, nxy)); pxs = 0.0
ALLOCATE(lkg(ng, nzMOC, nxy)); lkg = 0.0

ALLOCATE(S(ng, ng, nzCMFD, nxy)); S = 0.0
ALLOCATE(F(ng, nzCMFD, nxy)); F = 0.0
ALLOCATE(Chi(ng, nzCMFD, nxy)); Chi = 0.0

END SUBROUTINE

SUBROUTINE LinearMOCDriver(PinXS, eigv)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)
REAL :: eigv

INTEGER :: ixy, nxy
INTEGER :: iter, itermax = 20

nxy = mklGeom%nxy

CALL SetSourceOperator(PinXS)
CALL SetFlux()
CALL SetPsi(0)
CALL SetLeakage(PinXS)
CALL SetPseudoAbsorption()
! CALL LeakageSplit(PinXS)
CALL SetCrossSection(PinXS)
CALL SetLinearRTCoeff(PinXS, TRUE)

!PRINT *, 'negative feed flux', mklGeom%myzb, mklGeom%myze, COUNT(phisCoeff(:, :, :, 0) .LT. 0.0)
!PRINT *, 'negative incoming', mklGeom%myzb, mklGeom%myze, COUNT(mklAxial%PhiAngIn .LT. 0.0)

DO iter = 1, itermax
  CALL SetPsi(1)
  CALL SetSource(eigv, 0); CALL SetSource(eigv, 1)
  CALL SetSourceMoment()
  CALL SetLinearRTCoeff(PinXS, FALSE)
  phi = 0.0; phimx = 0.0
  !$OMP PARALLEL DO SCHEDULE(GUIDED)
  DO ixy = 1, nxy
    CALL LinearRayTrace(ixy, FALSE)
  ENDDO
  !$OMP END PARALLEL DO
  CALL SetBoundaryFlux()
  CALL SetFluxMoment()
ENDDO

phi = 0.0; phimx = 0.0; mklAxial%Jout = 0.0
CALL SetLinearRTCoeff(PinXS, FALSE)
!$OMP PARALLEL DO SCHEDULE(GUIDED)
DO ixy = 1, nxy
  CALL LinearRayTrace(ixy, TRUE)
ENDDO
!$OMP END PARALLEL DO
CALL SetBoundaryFlux()
CALL SetFluxMoment()
CALL SetFluxShape()

!PRINT *, 'negative flux', mklGeom%myzb, mklGeom%myze, COUNT(phisCoeff(:, :, :, 0) .LT. 0.0)
!PRINT *, 'negative outgoing', mklGeom%myzb, mklGeom%myze, COUNT(mklAxial%PhiAngOut .LT. 0.0)

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

SUBROUTINE SetFlux()

IMPLICIT NONE

INTEGER :: ig, ipin, iz, izf
INTEGER :: ng, nxy, nzCMFD

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      DO ig = 1, ng
        phisCoeff(ig, izf, ipin, 0) = mklCMFD%phis(ipin, iz, ig)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

CALL vdmul(ng * nxy * nzMOC, phisShape, phisCoeff(:, :, :, 0), phisCoeff(:, :, :, 0))

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

SUBROUTINE SetPsi(Order)

IMPLICIT NONE

INTEGER :: Order

INTEGER :: ig, igs, ipin, iz, izf
INTEGER :: ng, nxy, nzCMFD

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

psi(:, :, Order) = 0.0

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      DO ig = 1, ng
        psi(izf, ipin, Order) = psi(izf, ipin, Order) + F(ig, iz, ipin) * phisCoeff(ig, izf, ipin, Order)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

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
      DO ibd = 1, 4
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

! SUBROUTINE LeakageSplit(PinXS)
! USE TYPEDEF,        ONLY : PinXS_Type
! IMPLICIT NONE
! 
! TYPE(PinXS_Type), POINTER :: PinXS(:, :)
! 
! TYPE(superPin_Type), POINTER :: Pin(:)
! INTEGER :: ig, ibd, ipin, ineighpin, ipin_map, iz, izc, izf
! INTEGER :: ng, nxy, nzCMFD
! INTEGER, POINTER :: pinMap(:), pinMapRev(:), planeMap(:)
! REAL :: Dtil, pDhat(2), myphi, neighphi
! REAL, POINTER :: pLkg(:, :, :, :)
! 
! ng = mklGeom%ng
! nxy = mklGeom%nxy
! nzCMFD = mklGeom%nzCMFD
! Pin => mklGeom%superPin
! pinMap => mklGeom%pinMap
! pinMapRev => mklGeom%pinMapRev
! planeMap => mklGeom%planeMap
! 
! ALLOCATE(pLkg(2, ng, nxy, 0 : nzCMFD + 1)); pLkg = 0.0
! 
! DO iz = 1, nzCMFD
!   izc = planeMap(iz)
!   !$OMP PARALLEL PRIVATE(ipin_map, ineighpin, myphi, neighphi, Dtil, pDhat)
!   !$OMP DO SCHEDULE(GUIDED)
!   DO ipin = 1, nxy
!     ipin_map = pinMap(ipin)
!     DO ig = 1, ng
!       myphi = mklCMFD%phis(ipin, iz, ig)
!       DO ibd = 1, 4
!         ineighpin = Pin(ipin_map)%NeighIdx(ibd)
!         ineighpin = pinMapRev(ineighpin)
!         IF (ineighpin .LE. 0) THEN
!           neighphi = myphi
!         ELSE
!           neighphi = mklCMFD%phis(ineighpin, iz, ig)
!         ENDIF
!         Dtil = PinXS(ipin_map, izc)%Dtil(ibd, ig)
!         pDhat = PinXS(ipin_map, izc)%partialDhat(:, ibd, ig)
!         pLkg(in, ig, ipin, iz) = pLkg(in, ig, ipin, iz) - 0.5 * Dtil * (neighphi - myphi) - pDhat(in) * neighphi
!         pLkg(out, ig, ipin, iz) = pLkg(out, ig, ipin, iz) - 0.5 * Dtil * (neighphi - myphi) + pDhat(out) * myphi
!       ENDDO
!       pLkg(:, ig, ipin, iz) = pLkg(:, ig, ipin, iz) * mklGeom%hzfm(iz) / mklGeom%PinVolFm(ipin, iz)
!     ENDDO
!   ENDDO
!   !$OMP END DO
!   !$OMP END PARALLEL
! ENDDO
! 
! pxs = 0.0
! 
! !$OMP PARALLEL
! !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
! DO ipin = 1, nxy
!   DO iz = 1, nzCMFD
!     IF (.NOT. mklGeom%lRefCell(iz, ipin)) CYCLE
!     DO izf = fmRange(iz, 1), fmRange(iz, 2)
!       DO ig = 1, ng
!         lkg(ig, izf, ipin) = pLkg(in, ig, ipin, iz)
!         IF (mklCMFD%phis(ipin, iz, ig) .GT. 0.0) THEN
!           pxs(ig, izf, ipin) = pLkg(out, ig, ipin, iz) / mklCMFD%phis(ipin, iz, ig)
!         ELSE
!           lkg(ig, izf, ipin) = lkg(ig, izf, ipin) + pLkg(out, ig, ipin, iz)
!         ENDIF
!       ENDDO
!     ENDDO
!   ENDDO
! ENDDO
! !$OMP END DO
! !$OMP END PARALLEL
! 
! DEALLOCATE(pLkg)
! 
! END SUBROUTINE

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
      IF (lkg(ig, izf, ipin) .GT. 0.0 .AND. phisCoeff(ig, izf, ipin, 0) .GT. 0.0) THEN
        pxs(ig, izf, ipin) = lkg(ig, izf, ipin) / phisCoeff(ig, izf, ipin, 0)
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

SUBROUTINE SetSourceMoment()

INTEGER :: ig, igs, ipin, iz, izf
INTEGER :: ng, nxy, nzCMFD
INTEGER :: gb, ge

ng = mklGeom%ng
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
        gb = mklGeom%InScatRange(1, ig)
        ge = mklGeom%InScatRange(2, ig)
        DO igs = gb, ge
          srcm(1, ig, izf, ipin) = srcm(1, ig, izf, ipin) + 3.0 * mklAxial%SmP1(igs, ig, iz, ipin) * phimCoeff(1, igs, izf, ipin)
        ENDDO
        srcm(:, ig, izf, ipin) = srcm(:, ig, izf, ipin) / xst(ig, izf, ipin)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetSource(eigv, Order)

IMPLICIT NONE

REAL :: eigv, reigv
INTEGER :: Order

INTEGER :: ig, igs, ipin, iz, izf
INTEGER :: ng, nxy, nzCMFD
INTEGER :: gb, ge

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

reigv = 1.0 / eigv

src(:, :, :, Order) = 0.0

!$OMP PARALLEL PRIVATE(gb, ge)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      DO ig = 1, ng
        gb = mklGeom%InScatRange(1, ig)
        ge = mklGeom%InScatRange(2, ig)
        DO igs = gb, ge
          src(ig, izf, ipin, Order) = src(ig, izf, ipin, Order) + S(igs, ig, iz, ipin) * phisCoeff(igs, izf, ipin, Order)
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
        src(ig, izf, ipin, Order) = src(ig, izf, ipin, Order) + reigv * Chi(ig, iz, ipin) * psi(izf, ipin, Order)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

IF (Order .EQ. 0) THEN
  CALL vdsub(ng * nxy * nzMOC, src(:, :, :, Order), lkg, src(:, :, :, Order))
  CALL vddiv(ng * nxy * nzMOC, src(:, :, :, Order), xst, src(:, :, :, Order))
ENDIF

END SUBROUTINE

SUBROUTINE SetLinearRTCoeff(PinXS, lFirst)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)
LOGICAL :: lFirst

TYPE(mklAngle_Type), POINTER :: Angle(:)
INTEGER :: ng, nxy, nPolarAngle
INTEGER :: ig, ipol, ipin, ipin_map, iz, izf
INTEGER :: myzb, myze
INTEGER, POINTER :: pinMap(:)
REAL :: ex, s, sigt, tau

ng = mklGeom%ng
nxy = mklGeom%nxy
nPolarAngle = mklGeom%nPolar1D
myzb = mklGeom%myzb
myze = mklGeom%myze
pinMap => mklGeom%pinMap
Angle => mklAxial%Angle

IF (lFirst) THEN
  !$OMP PARALLEL PRIVATE(ex, s, tau, ipin_map)
  !$OMP DO SCHEDULE(GUIDED)
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    DO iz = myzb, myze
      DO izf = cmRange(iz, 1), cmRange(iz, 2)
        DO ig = 1, ng
          DO ipol = 1, nPolarAngle
            s = hzMOC(izf) * Angle(ipol)%rcosv
            tau = xst(ig, izf, ipin) * s
            ex = 1.0 - EXP(- tau)
            E1(ipol, ig, izf, ipin) = ex
            E3(ipol, ig, izf, ipin) = 2.0 * (tau - ex) - tau * ex
            R1(ipol, ig, izf, ipin) = (1.0 + tau / 2.0 - (1.0 + 1.0 / tau) * ex) / tau
            R3(ipol, ig, izf, ipin) = tau / 6.0 - 2.0 / tau - 2.0 + (1.0 + 1.0 / tau) * (1.0 + 2.0 / tau) * ex
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ELSE
  !$OMP PARALLEL PRIVATE(sigt, ipin_map)
  !$OMP DO SCHEDULE(GUIDED)
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    DO iz = myzb, myze
      DO izf = cmRange(iz, 1), cmRange(iz, 2)
        DO ig = 1, ng
          sigt = xst(ig, izf, ipin)
          DO ipol = 1, nPolarAngle
            Qc(ipol, ig, izf, 1, ipin) = src(ig, izf, ipin, 0) + Comp(1, ipol) * srcm(1, ig, izf, ipin)
            Qc(ipol, ig, izf, 2, ipin) = src(ig, izf, ipin, 0) - Comp(1, ipol) * srcm(1, ig, izf, ipin)
            Qz(ipol, ig, izf, :, ipin) = Angle(ipol)%cosv * src(ig, izf, ipin, 1) / 2.0 / sigt ** 2
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDIF

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

SUBROUTINE LinearRayTrace(ipin, lJout)

IMPLICIT NONE

TYPE(mklAngle_Type), POINTER :: Angle(:)
REAL, POINTER :: Jout(:, :, :, :), PhiAngIn(:, :, :, :), PhiAngOut(:, :, :, :)
REAL, POINTER :: myQc(:, :, :), myQz(:, :, :)
REAL, POINTER :: myE1(:, :, :), myE3(:, :, :), myR1(:, :, :), myR3(:, :, :)
REAL :: del_phi, del_phim, track_phi(mklGeom%nPolar1D, mklGeom%ng)
INTEGER :: ig, ipin, ipol, iz, izf
INTEGER :: ng, nPolarAngle, nzCMFD
LOGICAL :: lJout

IF (mklCntl%lRefPinFDM) THEN
  IF (mklGeom%lRefPin(ipin)) RETURN
ENDIF

Jout => mklAxial%Jout(:, :, :, :, ipin)
PhiAngIn => mklAxial%PhiAngIn
PhiAngOut => mklAxial%PhiAngOut
myE1 => E1(:, :, :, ipin)
myE3 => E3(:, :, :, ipin)
myR1 => R1(:, :, :, ipin)
myR3 => R3(:, :, :, ipin)
Angle => mklAxial%Angle

ng = mklGeom%ng
nPolarAngle = mklGeom%nPolar1D
nzCMFD = mklGeom%nzCMFD

!--- Upward Sweep

track_phi = PhiAngIn(:, :, ipin, bottom)

myQc => Qc(:, :, :, 1, ipin)
myQz => Qz(:, :, :, 1, ipin)

DO iz = 1, nzCMFD
  IF (lJout) THEN
    DO ig = 1, ng
      DO ipol = 1, nPolarAngle
        Jout(in, ig, bottom, iz) = Jout(in, ig, bottom, iz) + Angle(ipol)%wtsurf * track_phi(ipol, ig)
      ENDDO
    ENDDO
  ENDIF
  DO izf = fmRange(iz, 1), fmRange(iz, 2)
    DO ig = 1, ng
      DO ipol = 1, nPolarAngle
        del_phi = (myQc(ipol, ig, izf) - track_phi(ipol, ig)) * myE1(ipol, ig, izf)                                 &
                  + myQz(ipol, ig, izf) * myE3(ipol, ig, izf)
        del_phim = track_phi(ipol, ig) * 0.5 + (myQc(ipol, ig, izf) - track_phi(ipol, ig)) * myR1(ipol, ig, izf)    &
                   + myQz(ipol, ig, izf) * myR3(ipol, ig, izf)
        phi(0, ig, izf, ipin, 1) = phi(0, ig, izf, ipin, 1) - Angle(ipol)%wtsurf * del_phi
        phi(1, ig, izf, ipin, 1) = phi(1, ig, izf, ipin, 1) - mwt(1, ipol) * del_phi
        phimx(ig, izf, ipin) = phimx(ig, izf, ipin) + Angle(ipol)%wt * del_phim
        track_phi(ipol, ig) = track_phi(ipol, ig) + del_phi
      ENDDO
    ENDDO
  ENDDO
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

myQc => Qc(:, :, :, 2, ipin)
myQz => Qz(:, :, :, 2, ipin)

DO iz = nzCMFD, 1, -1
  IF (lJout) THEN
    DO ig = 1, ng
      DO ipol = 1, nPolarAngle
        Jout(in, ig, top, iz) = Jout(in, ig, top, iz) + Angle(ipol)%wtsurf * track_phi(ipol, ig)
      ENDDO
    ENDDO
  ENDIF
  DO izf = fmRange(iz, 2), fmRange(iz, 1), -1
    DO ig = 1, ng
      DO ipol = 1, nPolarAngle
        del_phi = (myQc(ipol, ig, izf) - track_phi(ipol, ig)) * myE1(ipol, ig, izf)                                 &
                  - myQz(ipol, ig, izf) * myE3(ipol, ig, izf)
        del_phim = track_phi(ipol, ig) * 0.5 + (myQc(ipol, ig, izf) - track_phi(ipol, ig)) * myR1(ipol, ig, izf)    &
                   - myQz(ipol, ig, izf) * myR3(ipol, ig, izf)
        phi(0, ig, izf, ipin, 2) = phi(0, ig, izf, ipin, 2) - Angle(ipol)%wtsurf * del_phi
        phi(1, ig, izf, ipin, 2) = phi(1, ig, izf, ipin, 2) + mwt(1, ipol) * del_phi
        phimx(ig, izf, ipin) = phimx(ig, izf, ipin) - Angle(ipol)%wt * del_phim
        track_phi(ipol, ig) = track_phi(ipol, ig) + del_phi
      ENDDO
    ENDDO
  ENDDO
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

SUBROUTINE SetFluxMoment()

IMPLICIT NONE

TYPE(mklAngle_Type), POINTER :: Angle(:)
INTEGER :: ipol, ig, ipin, izf
INTEGER :: ng, nxy, nPolarAngle
REAL :: myphis(2), myphim(2)

ng = mklGeom%ng
nxy = mklGeom%nxy
nPolarAngle = mklGeom%nPolar1D
Angle => mklAxial%Angle

!$OMP PARALLEL PRIVATE(myphis, myphim)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ipin = 1, nxy
  DO izf = 1, nzMOC
    DO ig = 1, ng
      myphis = phi(0, ig, izf, ipin, :) / (hzMOC(izf) * xst(ig, izf, ipin))
      myphim = phi(1, ig, izf, ipin, :) / (hzMOC(izf) * xst(ig, izf, ipin))
      phisCoeff(ig, izf, ipin, 0) = myphis(1) + myphis(2) + src(ig, izf, ipin, 0)
      phisCoeff(ig, izf, ipin, 1) = (- 6.0 * (myphis(1) - myphis(2)) + 12.0 * phimx(ig, izf, ipin)) / hzMOC(izf)
      phimCoeff(1, ig, izf, ipin) = myphim(1) + myphim(2) + 1.0 / 3.0 * srcm(1, ig, izf, ipin)
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
      phisCoeff(:, izf, ipin, 0) = 1.0
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDIF

!!$OMP PARALLEL
!!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      DO ig = 1, ng
        mklAxial%phic(ig, iz, ipin) = mklAxial%phic(ig, iz, ipin) + phisCoeff(ig, izf, ipin, 0)
      ENDDO
    ENDDO
    mklAxial%phic(:, iz, ipin) = mklAxial%phic(:, iz, ipin) / nDiv(iz)
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      DO ig = 1, ng
        phisShape(ig, izf, ipin) = phisCoeff(ig, izf, ipin, 0) / mklAxial%phic(ig, iz, ipin)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!!$OMP END DO
!!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(lNegative)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD
    DO ig = 1, ng
      lNegative = .FALSE.
      DO izf = fmRange(iz, 1), fmRange(iz, 2)
        IF (phisCoeff(ig, izf, ipin, 0) .LT. 0.0) lNegative = .TRUE.
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

END MODULE

#endif    