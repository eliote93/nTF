#include <defines.h>
!--- CNJ Edit : 1D Axial Nodal Modules with Intel MKL
#ifdef __INTEL_MKL

MODULE MKL_NODAL

USE MKL_3D
IMPLICIT NONE

REAL, POINTER, PRIVATE :: phisCoeff(:, :, :, :), srcCoeff(:, :, :, :), lkgCoeff(:, :, :, :), psi(:, :), pxs(:, :, :)

PRIVATE
PUBLIC :: AllocNodal, NodalDriver

CONTAINS

!--- Public Routines ------------------------------------------------------------------------------

SUBROUTINE AllocNodal()

IMPLICIT NONE

INTEGER :: ng, nxy, nzCMFD

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

ALLOCATE(phisCoeff(0 : 4, ng, nzCMFD, nxy)); phisCoeff = 0.0
ALLOCATE(srcCoeff(0 : 4, ng, nzCMFD, nxy)); srcCoeff = 0.0
ALLOCATE(lkgCoeff(0 : 2, ng, nzCMFD, nxy)); lkgCoeff = 0.0
ALLOCATE(psi(nzCMFD, nxy)); psi = 0.0
ALLOCATE(pxs(ng, nzCMFD, nxy)); pxs = 0.0

END SUBROUTINE

SUBROUTINE NodalDriver(PinXS, eigv)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)
REAL :: eigv

INTEGER :: ng, nxy, nzCMFD
INTEGER :: ipin, ipin_map, ig, iz, izf, iter
INTEGER, POINTER :: pinMap(:), planeMap(:)
REAL, POINTER :: a(:, :), q(:, :), l(:, :), jL(:, :), jR(:, :)
REAL :: h, err, errmax, tol = 1.0D-03
LOGICAL :: lConv

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
pinMap => mklGeom%pinMap
planeMap => mklGeom%planeMap

CALL SetAverageFlux()
CALL SetTransverseLeakage(PinXS)
! CALL SetPseudoAbsorption()
CALL SetAveragePsi(PinXS)

iter = 0; lConv = .FALSE.
DO iter = 1, 50 ! WHILE (.NOT. lConv)
  errmax = 0.0
  !$OMP PARALLEL PRIVATE(iz, ipin_map, a, q, l, jL, jR, h, err) REDUCTION(MAX : errmax)
  !$OMP DO SCHEDULE(GUIDED)
  DO ipin = 1, nxy
    IF (mklCntl%lRefPinFDM) THEN
      IF (mklGeom%lRefPin(ipin)) CYCLE
    ENDIF
    ipin_map = pinMap(ipin)
    DO izf = 1, nzCMFD
      iz = planeMap(izf)
      h = mklGeom%hzfm(izf)
      a => phisCoeff(:, :, izf, ipin); q => srcCoeff(:, :, izf, ipin); l => lkgCoeff(:, :, izf, ipin)
      jL => mklAxial%Jout(1 : 2, :, bottom, izf, ipin); jR => mklAxial%Jout(1 : 2, :, top, izf, ipin)
      CALL SetSrcCoeff(PinXS(ipin_map, iz), a, q, l, psi(izf, ipin), eigv)
      IF (mklCntl%lSENM) THEN
        CALL SolveSENM(PinXS(ipin_map, iz), pxs(:, izf, ipin), a, q, jL, jR, h, err)
      ELSE
        CALL SolveNEM(PinXS(ipin_map, iz), pxs(:, izf, ipin), a, q, jL, jR, h, err)
      ENDIF
      errmax = max(err, errmax)
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
  CALL GetIncomingCurrent()
  lConv = errmax .LE. tol
!  iter = iter + 1
ENDDO

PRINT *, COUNT(phisCoeff(0, :, :, :) .LE. 0.0)

END SUBROUTINE

!--- Private Rotines ------------------------------------------------------------------------------

SUBROUTINE SetAverageFlux()

IMPLICIT NONE

INTEGER :: ng, nxy, nzCMFD
INTEGER :: ig, ipin, iz

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO iz = 1, nzCMFD
    DO ipin = 1, nxy
      phisCoeff(0, ig, iz, ipin) = mklCMFD%phis(ipin, iz, ig)
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetAveragePsi(PinXS)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)

INTEGER :: ig, ipin, ipin_map, iz, izf
INTEGER :: ng, nxy, nzCMFD
INTEGER, POINTER :: pinMap(:), planeMap(:)

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
pinMap => mklGeom%pinMap
planeMap => mklGeom%planeMap

psi = 0.0

!$OMP PARALLEL PRIVATE(ipin_map, iz)
!$OMP DO SCHEDULE(GUIDED)
DO ipin = 1, nxy
  ipin_map = pinMap(ipin)
  DO izf = 1, nzCMFD
    iz = planeMap(izf)
    DO ig = 1, ng
      psi(izf, ipin) = psi(izf, ipin) + phisCoeff(0, ig, izf, ipin) * PinXS(ipin_map, iz)%XSnf(ig)
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetSrcCoeff(PinXS, a, q, l, psi, eigv)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type) :: PinXS
REAL :: a(0 : 4, mklGeom%ng), q(0 : 4, mklGeom%ng), l(0 : 2, mklGeom%ng)
REAL :: eigv, reigv, psi, psiCoeff(0 : 4)
INTEGER :: ig, igs, gb, ge, ng

ng = mklGeom%ng
reigv = 1.0 / eigv

psiCoeff = 0.0
psiCoeff(0) = psi
DO ig = 1, ng
  psiCoeff(1 : 4) = psiCoeff(1 : 4) + a(1 : 4, ig) * PinXS%XSnf(ig)
ENDDO

q = 0.0
DO ig = 1, ng
  gb = PinXS%XSs(ig)%ib
  ge = PinXS%XSs(ig)%ie
  DO igs = gb, ge
    q(0 : 4, ig) = q(0 : 4, ig) + a(0 : 4, igs) * PinXS%XSs(ig)%from(igs)
  ENDDO
  q(0 : 4, ig) = q(0 : 4, ig) + reigv * psiCoeff(0 : 4) * PinXS%Chi(ig)
  q(0 : 2, ig) = q(0 : 2, ig) - l(0 : 2, ig)
ENDDO

END SUBROUTINE

! SUBROUTINE SetIntraNodalLeakage(PinXS)
! USE TYPEDEF,        ONLY : PinXS_Type
! IMPLICIT NONE
! 
! TYPE(PinXS_Type), POINTER :: PinXS(:, :)
! 
! TYPE(superPin_Type), POINTER :: Pin(:)
! INTEGER :: ng, nxy, nzCMFD
! INTEGER :: ig, ibd, ipin, ipin_map, ineighpin, iz, izf
! INTEGER, POINTER :: pinMap(:), pinMapRev(:), planeMap(:)
! REAL :: myphi(0 : 2), neighphi(0 : 2), lkg(0 : 2)
! REAL :: Dtil, Dhat, pDhat(2)
! REAL, POINTER :: intraFlux(:, :, :, :)
! 
! Pin => mklGeom%superPin
! ng = mklGeom%ng
! nxy = mklGeom%nxy
! nzCMFD = mklGeom%nzCMFD
! pinMap => mklGeom%pinMap
! pinMapRev => mklGeom%pinMapRev
! planeMap => mklGeom%planeMap
! 
! ALLOCATE(intraFlux(0 : 2, ng, nzCMFD, nxy))
! 
! !$OMP PARALLEL
! !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
! DO ipin = 1, nxy
!   DO izf = 1, nzCMFD
!     DO ig = 1, ng
!       IF (mklCntl%lSENM) THEN
!         CALL IntegrateSENMFlux(phisCoeff(:, ig, izf, ipin), intraFlux(:, ig, izf, ipin))
!       ELSE
!         CALL IntegrateNEMFlux(phisCoeff(:, ig, izf, ipin), intraFlux(:, ig, izf, ipin))
!       ENDIF
!     ENDDO
!   ENDDO
! ENDDO
! !$OMP END DO
! !$OMP END PARALLEL
! 
! DO izf = 1, nzCMFD
!   iz = planeMap(izf)
!   !$OMP PARALLEL PRIVATE(ipin_map, ineighpin, myphi, neighphi, lkg, Dtil, Dhat, pDhat)
!   !$OMP DO SCHEDULE(GUIDED)
!   DO ipin = 1, nxy
!     ipin_map = pinMap(ipin)
!     DO ig = 1, ng
!       myphi = intraFlux(:, ig, izf, ipin)
!       lkg = 0.0
!       DO ibd = 1, 4
!         ineighpin = Pin(ipin_map)%NeighIdx(ibd)
!         ineighpin = pinMapRev(ineighpin)
!         IF (ineighpin .EQ. VoidCell) THEN
!           neighphi = 0.0
!         ELSEIF (ineighpin .EQ. RefCell) THEN
!           neighphi = myphi
!         ELSE
!           neighphi = intraFlux(:, ig, izf, ineighpin)
!         ENDIF
!         IF (mklCntl%pCMFD) THEN
!           Dtil = PinXS(ipin_map, iz)%Dtil(ibd, ig)
!           pDhat = PinXS(ipin_map, iz)%partialDhat(:, ibd, ig)
!           lkg = lkg - Dtil * (neighphi - myphi) - (pDhat(in) * neighphi - pDhat(out) * myphi)
!         ELSE
!           Dtil = PinXS(ipin_map, iz)%Dtil(ibd, ig)
!           Dhat = PinXS(ipin_map, iz)%Dhat(ibd, ig)
!           lkg = lkg - Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
!         ENDIF
!       ENDDO
!       lkg = lkg * mklGeom%hzfm(izf) / mklGeom%PinVolFm(ipin, izf)
!       IF (mklCntl%lSENM) THEN
!         CALL ExpandSENMLeakage(lkg, lkgCoeff(:, ig, izf, ipin))
!       ELSE
!         CALL ExpandNEMLeakage(lkg, lkgCoeff(:, ig, izf, ipin))
!       ENDIF
!     ENDDO
!   ENDDO
!   !$OMP END DO
!   !$OMP END PARALLEL
! ENDDO
! 
! DEALLOCATE(intraFlux)
! 
! CONTAINS
! 
! SUBROUTINE IntegrateNEMFlux(phi, intraFlux)
! 
! IMPLICIT NONE
! 
! REAL :: phi(0 : 4), intraFlux(0 : 2)
! 
! intraFlux(0) = 4.0 / 81.0 * phi(4) - 4.0 / 27.0 * phi(3) - 2.0 / 27.0 * phi(2) - 2.0 / 9.0 * phi(1) + 1.0 / 3.0 * phi(0)
! intraFlux(1) = - 8.0 / 81.0 * phi(4) + 4.0 / 27.0 * phi(2) + 1.0 / 3.0 * phi(0)
! intraFlux(2) = 4.0 / 81.0 * phi(4) + 4.0 / 27.0 * phi(3) - 2.0 / 27.0 * phi(2) + 2.0 / 9.0 * phi(1) + 1.0 / 3.0 * phi(0)
! 
! intraFlux = intraFlux * 3.0
! 
! END SUBROUTINE
! 
! SUBROUTINE ExpandNEMLeakage(lkg, intraLkg)
! 
! IMPLICIT NONE
! 
! REAL :: lkg(0 : 2), intraLkg(0 : 2)
! 
! intraLkg(0) = 1.0 / 3.0 * (lkg(0) + lkg(1) + lkg(2))
! intraLkg(1) = - 3.0 / 4.0 * (lkg(0) - lkg(2))
! intraLkg(2) = - 3.0 / 4.0 * (lkg(0) - 2.0 * lkg(1) + lkg(2))
! 
! END SUBROUTINE
! 
! SUBROUTINE IntegrateSENMFlux(phi, intraFlux)
! 
! IMPLICIT NONE
! 
! REAL :: phi(0 : 4), intraFlux(0 : 2)
! 
! intraFlux(0) = - 20.0 / 243.0 * phi(4) + 4.0 / 81.0 * phi(3) + 4.0 / 27.0 * phi(2) - 4.0 / 9.0 * phi(1) + 2.0 / 3.0 * phi(0)
! intraFlux(1) = 40.0 / 243.0 * phi(4) - 8.0 / 27.0 * phi(2) + 2.0 / 3.0 * phi(0)
! intraFlux(2) = - 20.0 / 243.0 * phi(4) - 4.0 / 81.0 * phi(3) + 4.0 / 27.0 * phi(2) + 4.0 / 9.0 * phi(1) + 2.0 / 3.0 * phi(0)
! 
! intraFlux = intraFlux * 3.0 / 2.0
! 
! END SUBROUTINE
! 
! SUBROUTINE ExpandSENMLeakage(lkg, intraLkg)
! 
! IMPLICIT NONE
! 
! REAL :: lkg(0 : 2), intraLkg(0 : 2)
! 
! intraLkg(0) = 1.0 / 3.0 * (lkg(0) + lkg(1) + lkg(2))
! intraLkg(1) = - 3.0 / 4.0 * (lkg(0) - lkg(2))
! intraLkg(2) = 3.0 / 4.0 * (lkg(0) - 2.0 * lkg(1) + lkg(2))
! 
! END SUBROUTINE
! 
! END SUBROUTINE

SUBROUTINE SetTransverseLeakage(PinXS)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)

TYPE(superPin_Type), POINTER :: Pin(:)
INTEGER :: ng, nxy, nzCMFD
INTEGER :: ig, ibd, ipin, ipin_map, ineighpin, iz, izf
INTEGER, POINTER :: pinMap(:), pinMapRev(:), planeMap(:)
REAL :: myphi, neighphi, lkg
REAL :: Dtil, Dhat
REAL :: h0, h2
REAL, POINTER :: l(:, :, :), L0(:, :), L1(:, :), L2(:, :)
REAL, POINTER :: radLkg(:, :, :)

Pin => mklGeom%superPin
ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
pinMap => mklGeom%pinMap
pinMapRev => mklGeom%pinMapRev
planeMap => mklGeom%planeMap

ALLOCATE(radLkg(ng, nxy, 0 : nzCMFD + 1)); radLkg = 0.0

DO izf = 1, nzCMFD
  iz = planeMap(izf)
  !$OMP PARALLEL PRIVATE(ipin_map, ineighpin, myphi, neighphi, lkg, Dtil, Dhat)
  !$OMP DO SCHEDULE(GUIDED)
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    DO ig = 1, ng
      myphi = mklCMFD%phis(ipin, izf, ig)
      lkg = 0.0
      DO ibd = 1, 4
        ineighpin = Pin(ipin_map)%NeighIdx(ibd)
        ineighpin = pinMapRev(ineighpin)
        IF (ineighpin .EQ. VoidCell) THEN
          neighphi = 0.0
        ELSEIF (ineighpin .EQ. RefCell) THEN
          neighphi = myphi
        ELSE
          neighphi = mklCMFD%phis(ineighpin, izf, ig)
        ENDIF
        Dtil = PinXS(ipin_map, iz)%Dtil(ibd, ig)
        Dhat = PinXS(ipin_map, iz)%Dhat(ibd, ig)
        lkg = lkg - Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
      ENDDO
      radLkg(ig, ipin, izf) = lkg * mklGeom%hzfm(izf) / mklGeom%PinVolFm(ipin, izf)
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
  l => lkgCoeff(:, :, iz, :)
  L0 => radLkg(:, :, iz - 1); L1 => radLkg(:, :, iz); L2 => radLkg(:, :, iz + 1)
  IF (mklCntl%lSENM) THEN
    h0 = mklGeom%hzfm(iz - 1) / mklGeom%hzfm(iz) * 2.0
    h2 = mklGeom%hzfm(iz + 1) / mklGeom%hzfm(iz) * 2.0
    CALL SetSENMLkgCoeff(l, L0, L1, L2, h0, h2)
  ELSE
    h0 = mklGeom%hzfm(iz - 1) / mklGeom%hzfm(iz)
    h2 = mklGeom%hzfm(iz + 1) / mklGeom%hzfm(iz)
    CALL SetNEMLkgCoeff(l, L0, L1, L2, h0, h2)
  ENDIF
ENDDO

DEALLOCATE(radLkg)

CONTAINS

SUBROUTINE SetNEMLkgCoeff(l, L0, L1, L2, h0, h2)

IMPLICIT NONE

REAL :: l(0 : 2, mklGeom%ng, mklGeom%nxy)
REAL :: L0(:, :), L1(:, :), L2(:, :)
REAL :: h0, h2, n1(2), n2(2), d1, d2
INTEGER :: ng, nxy
INTEGER :: ig, ipin

ng = mklGeom%ng
nxy = mklGeom%nxy

n1(1) = 1.0 + 3.0 * h2 + 2.0 * h2 ** 2
n1(2) = 1.0 + 3.0 * h0 + 2.0 * h0 ** 2
d1 = 2.0 * (1.0 + h0) * (1.0 + h2) * (1.0 + h0 + h2)

n2(1) = 1.0 + h2
n2(2) = 1.0 + h0
d2 = 2.0 * (1.0 + h2) * (1.0 + 2.0 * h0 + h0 ** 2 + h2 + h0 * h2)

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO ig = 1, ng
    l(0, ig, ipin) = L1(ig, ipin)
    l(1, ig, ipin) = - (L0(ig, ipin) * n1(1) + L1(ig, ipin) * (n1(2) - n1(1)) - L2(ig, ipin) * n1(2)) / d1
    l(2, ig, ipin) = - (L0(ig, ipin) * n2(1) - L1(ig, ipin) * (n2(2) + n2(1)) + L2(ig, ipin) * n2(2)) / d2
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetSENMLkgCoeff(l, L0, L1, L2, h0, h2)

IMPLICIT NONE

REAL :: l(0 : 2, mklGeom%ng, mklGeom%nxy)
REAL :: L0(:, :), L1(:, :), L2(:, :)
REAL :: h0, h2, n1(2), n2(2), d1, d2
INTEGER :: ng, nxy
INTEGER :: ig, ipin

ng = mklGeom%ng
nxy = mklGeom%nxy

n1(1) = 2.0 + 3.0 * h2 + h2 ** 2
n1(2) = 2.0 + 3.0 * h0 + h0 ** 2
d1 = (2.0 + h0) * (2.0 + h2) * (2.0 + h0 + h2)

n2(1) = 2.0 + h2
n2(2) = 2.0 + h0
d2 = (2.0 + h2) * (4.0 + 4.0 * h0 + h0 ** 2 + 2.0 * h2 + h0 * h2)

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO ig = 1, ng
    l(0, ig, ipin) = L1(ig, ipin)
    l(1, ig, ipin) = - 2.0 * (L0(ig, ipin) * n1(1) + L1(ig, ipin) * (n1(2) - n1(1)) - L2(ig, ipin) * n1(2)) / d1
    l(2, ig, ipin) = 2.0 * (L0(ig, ipin) * n2(1) - L1(ig, ipin) * (n2(2) + n2(1)) + L2(ig, ipin) * n2(2)) / d2
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

END SUBROUTINE

SUBROUTINE SetPseudoAbsorption()

IMPLICIT NONE

INTEGER :: ig, ipin, iz
INTEGER :: ng, nxy, nzCMFD

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

pxs = 0.0

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD
    IF (.NOT. mklGeom%lRefCell(iz, ipin)) CYCLE
    DO ig = 1, ng
      IF (lkgCoeff(0, ig, iz, ipin) .GT. 0.0 .AND. phisCoeff(0, ig, iz, ipin) .GT. 0.0) THEN
        pxs(ig, iz, ipin) = lkgCoeff(0, ig, iz, ipin) / phisCoeff(0, ig, iz, ipin)
        lkgCoeff(:, ig, iz, ipin) = 0.0
      ENDIF
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetSurfaceFlux()

IMPLICIT NONE

INTEGER :: ng, nxy, nzCMFD
INTEGER :: ig, ipin, iz

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD
    DO ig = 1, ng
      mklAxial%Jout(surf, ig, bottom, iz, ipin) = 2.0 * sum(mklAxial%Jout(in : out, ig, bottom, iz, ipin))
      mklAxial%Jout(surf, ig, top, iz, ipin) = 2.0 * sum(mklAxial%Jout(in : out, ig, top, iz, ipin))
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE GetIncomingCurrent()

IMPLICIT NONE

INTEGER :: ng, nxy, nzCMFD
INTEGER :: ig, ipin, iz
REAL, POINTER :: myJout(:, :, :), neighJin(:, :, :)

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

ALLOCATE(myJout(ng, nxy, 2))
ALLOCATE(neighJin(ng, nxy, 2))

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO ig = 1, ng
    myJout(ig, ipin, bottom) = mklAxial%Jout(out, ig, bottom, 1, ipin)
    myJout(ig, ipin, top) = mklAxial%Jout(out, ig, top, nzCMFD, ipin)
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

CALL InitFastComm()

CALL GetNeighborFast(ng * nxy, myJout(:, :, bottom), neighJin(:, :, top), bottom)
CALL GetNeighborFast(ng * nxy, myJout(:, :, top), neighJin(:, :, bottom), top)

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO iz = 1, nzCMFD
  DO ipin = 1, nxy
    DO ig = 1, ng
      IF (iz .NE. 1) mklAxial%Jout(in, ig, bottom, iz, ipin) = mklAxial%Jout(out, ig, top, iz - 1, ipin)
      IF (iz .NE. nzCMFD) mklAxial%Jout(in, ig, top, iz, ipin) = mklAxial%Jout(out, ig, bottom, iz + 1, ipin)
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

CALL FinalizeFastComm()

IF (mklGeom%lBottom) THEN
  IF (mklGeom%AxBC(bottom) .EQ. VoidCell) neighJin(:, :, bottom) = 0.0
  IF (mklGeom%AxBC(bottom) .EQ. RefCell) CALL dcopy(ng * nxy, myJout(:, :, bottom), 1, neighJin(:, :, bottom), 1)
ENDIF

IF (mklGeom%lTop) THEN
  IF (mklGeom%AxBC(top) .EQ. VoidCell) neighJin(:, :, top) = 0.0
  IF (mklGeom%AxBC(top) .EQ. RefCell) CALL dcopy(ng * nxy, myJout(:, :, top), 1, neighJin(:, :, top), 1)
ENDIF

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO ig = 1, ng
    mklAxial%Jout(in, ig, bottom, 1, ipin) = neighJin(ig, ipin, bottom)
    mklAxial%Jout(in, ig, top, nzCMFD, ipin) = neighJin(ig, ipin, top)
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

DEALLOCATE(myJout)
DEALLOCATE(neighJin)

END SUBROUTINE

SUBROUTINE SolveNEM(PinXS, pxs, a, q, jL, jR, h, errmax)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type) :: PinXS
REAL :: pxs(mklGeom%ng), a(0 : 4, mklGeom%ng), q(0 : 4, mklGeom%ng), c(4)
REAL :: jL(2, mklGeom%ng), jR(2, mklGeom%ng)
REAL :: h, beta(mklGeom%ng), sigD(mklGeom%ng), sigR(mklGeom%ng)
REAL :: err, errmax
INTEGER :: ig, ng

ng = mklGeom%ng
beta = PinXS%XSD / h
sigD = PinXS%XSD / h ** 2
sigR = PinXS%XSr + pxs

errmax = 0.0

DO ig = 1, ng

  !--- Constants
  c(1) = 6.0 * beta(ig) / (1.0 + 12.0 * beta(ig))
  c(2) = - 8.0 * beta(ig) / ((1.0 + 4.0 * beta(ig)) * (1.0 + 12.0 * beta(ig)))
  c(3) = (1.0 - 48.0 * beta(ig) ** 2) / ((1.0 + 4.0 * beta(ig)) * (1.0 + 12.0 * beta(ig)))
  c(4) = 6.0 * beta(ig) / (1.0 + 4.0 * beta(ig))

  !--- Flux Coefficients
  a(1, ig) = jR(out, ig) + jR(in, ig) - jL(out, ig) - jL(in, ig)
  a(2, ig) = a(0, ig) - (jR(out, ig) + jR(in, ig) + jL(out, ig) + jL(in, ig))
  a(3, ig) = (5.0 * q(1, ig) + 3.0 * q(3, ig) - 5.0 * a(1, ig) * sigR(ig)) / (3.0 * (60.0 * sigD(ig) + sigR(ig)))
  a(4, ig) = (- 7.0 * q(2, ig) + 3.0 * q(4, ig) + 7.0 * a(2, ig) * sigR(ig)) / (3.0 * (140.0 * sigD(ig) + sigR(ig)))
  !--- Average Flux
  err = a(0, ig)
  a(0, ig) = q(0, ig) - (2.0 * a(4, ig) * c(1) - (1.0 - c(2) - c(3)) * (jL(in, ig) + jR(in, ig))) / h
  a(0, ig) = a(0, ig) / (sigR(ig) + 2.0 * c(1) / h)
  err = abs(err - a(0, ig)) / a(0, ig); errmax = max(err, errmax)
  !--- Outgoing Current
  jL(out, ig) = c(1) * (a(0, ig) + a(4, ig)) + c(3) * jL(in, ig) + c(2) * jR(in, ig) - c(4) * a(3, ig)
  jR(out, ig) = c(1) * (a(0, ig) + a(4, ig)) + c(2) * jL(in, ig) + c(3) * jR(in, ig) + c(4) * a(3, ig)

ENDDO

END SUBROUTINE

SUBROUTINE SolveSENM(PinXS, pxs, a, q, jL, jR, h, errmax)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type) :: PinXS
REAL :: pxs(mklGeom%ng), a(0 : 4, mklGeom%ng), q(0 : 4, mklGeom%ng), c(6)
REAL :: jL(2, mklGeom%ng), jR(2, mklGeom%ng), sigD(mklGeom%ng), sigR(mklGeom%ng)
REAL :: h, k, alpha(2), beta, tau, phi(2), j(2)
REAL :: coshk, sinhk
REAL :: err, errmax
INTEGER :: ig, ng

ng = mklGeom%ng
sigD = 4.0 * PinXS%XSD / h ** 2
sigR = PinXS%XSr + pxs

errmax = 0.0

DO ig = 1, ng

  !--- Constants
  beta = PinXS%XSD(ig) / h
  k = h / 2.0 * sqrt(sigR(ig) / PinXS%XSD(ig))
  coshk = dcosh(k); sinhk = dsinh(k)
  tau = k * sinhk / (4.0 * beta * k * sinhk + coshk - sinhk / k)

  !--- Legendre Expansion Coefficients
  c(1) = (1.0 / sigR(ig)) * (q(1, ig) + 15.0 / k ** 2 * q(3, ig))
  c(2) = (1.0 / sigR(ig)) * (q(2, ig) + 35.0 / k ** 2 * q(4, ig))
  c(3) = q(3, ig) / sigR(ig)
  c(4) = q(4, ig) / sigR(ig)
  !--- Average Flux
  err = a(0, ig)
  a(0, ig) = q(0, ig) + sigD(ig) * (2.0 * tau * (jL(in, ig) + jR(in, ig)) + (3.0 - tau * (1.0 + 12.0 * beta)) * c(2) + (10.0 - tau * (1.0 + 40.0 * beta)) * c(4))
  a(0, ig) = a(0, ig) / (sigR(ig) + sigD(ig) * tau)
  err = abs(err - a(0, ig)) / a(0, ig); errmax = max(err, errmax)
  !--- Hyperbolic Term Coefficients
  c(5) = 2.0 * (jR(in, ig) - jL(in, ig)) - (1.0 + 4.0 * beta) * c(1) - (1.0 + 24.0 * beta) * c(3)
  c(5) = c(5) / (4.0 * beta * k * coshk + sinhk)
  c(6) = - a(0, ig) + 2.0 * (jR(in, ig) + jL(in, ig)) - (1.0 + 12.0 * beta) * c(2) - (1.0 + 40.0 * beta) * c(4)
  c(6) = c(6) / (4.0 * beta * k * sinhk + coshk - sinhk / k)
  !--- Approximated Flux Coefficients
  a(1, ig) = c(1) + (3.0 / k) * (coshk - sinhk / k) * c(5)
  a(2, ig) = c(2) - 5.0 * (3.0 * coshk / k ** 2 - (1.0 + 3.0 / k ** 2) * sinhk / k) * c(6)
  a(3, ig) = c(3) + (7.0 / k) * ((1.0 + 15.0 / k ** 2) * coshk - (6.0 + 15.0 / k ** 2) * sinhk / k) * c(5)
  a(4, ig) = c(4) - 9.0 * (5.0 / k ** 2 * (2.0 + 21.0 / k ** 2) * coshk - (1.0 + 45.0 / k ** 2 + 105.0 / k ** 4) * sinhk / k) * c(6)
  !--- Outgoing Current
  alpha(1) = - beta * k * coshk + 0.25 * sinhk
  alpha(2) = - beta * k * sinhk + 0.25 * (coshk - sinhk / k)
  phi(1) = - c(1) + c(2) - c(3) + c(4)
  phi(2) = c(1) + c(2) + c(3) + c(4)
  j(1) = - 2.0 * beta * (c(1) - 3.0 * c(2) + 6.0 * c(3) - 10.0 * c(4))
  j(2) = - 2.0 * beta * (c(1) + 3.0 * c(2) + 6.0 * c(3) + 10.0 * c(4))
  jL(out, ig) = - alpha(1) * c(5) + alpha(2) * c(6) + 0.25 * (a(0, ig) + phi(1)) - 0.5 * j(1)
  jR(out, ig) = alpha(1) * c(5) + alpha(2) * c(6) + 0.25 * (a(0, ig) + phi(2)) + 0.5 * j(2)

ENDDO

END SUBROUTINE

END MODULE

MODULE MKL_SP3NODAL

USE MKL_3D
IMPLICIT NONE

CONTAINS

!--- Temporary Data Transfer to Use nTRACER Axial Solvers

SUBROUTINE CopyFlux(CoreInfo, CmInfo)
USE TYPEDEF,        ONLY : CoreInfo_Type,       CmInfo_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(CmInfo_Type) :: CmInfo

INTEGER :: ng, nxy, nzCMFD, nSubplane
INTEGER :: ig, ipin, ipin_map, iz, izf, izf_map(2)
INTEGER :: myzb, myze
INTEGER, POINTER :: pinMap(:), fmRange(:, :), SubplaneRange(:, :)

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
myzb = mklGeom%myzb
myze = mklGeom%myze
pinMap => mklGeom%pinMap
fmRange => mklGeom%fmRange
nSubplane = CoreInfo%nSubplane
SubplaneRange => CoreInfo%SubplaneRange

DO ig = 1, ng
  DO iz = myzb, myze
    DO izf = 1, nSubplane
      DO ipin = 1, nxy
        ipin_map = pinMap(ipin)
        izf_map(1) = SubplaneRange(1, iz) + izf - 1
        izf_map(2) = fmRange(iz, 1) + izf - 1
        CmInfo%PhiFm(ipin_map, izf_map(1), ig) = mklCMFD%phis(ipin, izf_map(2), ig)
      ENDDO
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE CopyDhat(CoreInfo, CmInfo)
USE TYPEDEF,        ONLY : CoreInfo_Type,       CmInfo_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(CmInfo_Type) :: CmInfo

INTEGER :: ng, nxy, nzCMFD, nSubplane
INTEGER :: ig, ipin, ipin_map, iz, izf, izf_map(2)
INTEGER :: myzb, myze
INTEGER, POINTER :: pinMap(:), fmRange(:, :), SubplaneRange(:, :)

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
myzb = mklGeom%myzb
myze = mklGeom%myze
pinMap => mklGeom%pinMap
fmRange => mklGeom%fmRange
nSubplane = CoreInfo%nSubplane
SubplaneRange => CoreInfo%SubplaneRange

DO ig = 1, ng
  DO iz = myzb, myze
    DO izf = 1, nSubplane
      DO ipin = 1, nxy
        ipin_map = pinMap(ipin)
        izf_map(1) = SubplaneRange(1, iz) + izf - 1
        izf_map(2) = fmRange(iz, 1) + izf - 1
        mklCMFD%AxDhat(:, ipin, izf_map(2), ig) = CmInfo%AxDhat(:, ipin_map, izf_map(1), ig)
      ENDDO
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE

END MODULE

#endif
