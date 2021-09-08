#include <defines.h>
!--- CNJ Edit : 1D Axial Flat MOC Modules with Intel MKL
#ifdef __INTEL_MKL

MODULE MKL_FLATMOC

USE MKL_3D

IMPLICIT NONE

REAL, POINTER, PRIVATE, DIMENSION(:)         :: hzMOC
REAL, POINTER, PRIVATE, DIMENSION(:,:)       :: psi, Comp, mwt
REAL, POINTER, PRIVATE, DIMENSION(:,:,:)     :: phis, src, lkg, xst, pxs, F, chi
REAL, POINTER, PRIVATE, DIMENSION(:,:,:,:)   :: phim, srcm, S, wtExp
REAL, POINTER, PRIVATE, DIMENSION(:,:,:,:,:) :: srca

INTEGER, POINTER, PRIVATE, DIMENSION(:)   :: cmMap, fmMap
INTEGER, POINTER, PRIVATE, DIMENSION(:,:) :: cmRange, fmRange
INTEGER, PRIVATE :: nzMOC, nDiv(100)

LOGICAL, PRIVATE :: lFirst = .TRUE.

PRIVATE

PUBLIC :: AllocFlatMOC, FlatMOCDriver

CONTAINS
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE AllocFlatMOC(Axial)

USE allocs

IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial

INTEGER :: ng, nxy, nzCMFD, nPolarAngle
! ----------------------------------------------------

ng          = Axial%ng
nxy         = mklGeom%nxy
nzCMFD      = mklGeom%nzCMFD
nPolarAngle = mklGeom%nPolar1D

CALL SetSphericalHarmonics(Axial)
CALL SetSubmesh

CALL dmalloc(psi,         nzMOC, nxy)
CALL dmalloc(phis,    ng, nzMOC, nxy)
CALL dmalloc(src,     ng, nzMOC, nxy)
CALL dmalloc(lkg,     ng, nzMOC, nxy)
CALL dmalloc(xst,     ng, nzMOC, nxy)
CALL dmalloc(pxs,     ng, nzMOC, nxy)
CALL dmalloc(phim, 3, ng, nzMOC, nxy)
CALL dmalloc(srcm, 3, ng, nzMOC, nxy)

CALL dmalloc(S,   ng, ng, nzCMFD, nxy)
CALL dmalloc(F,       ng, nzCMFD, nxy)
CALL dmalloc(Chi,     ng, nzCMFD, nxy)

CALL dmalloc(srca, 2, nPolarAngle, ng, nzMOC, nxy)
CALL dmalloc(wtExp,   nPolarAngle, ng, nzMOC, nxy)
! ----------------------------------------------------

END SUBROUTINE AllocFlatMOC
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE FlatMOCDriver(CMFD, Axial, eigv)

USE param,   ONLY : ZERO, FALSE, TRUE
USE TYPEDEF, ONLY : PinXS_Type
USE PE_MOD,  ONLY : PE
USE TIMER,   ONLY : nTracer_dclock, TimeChk

IMPLICIT NONE

TYPE (mklCMFD_Type)  :: CMFD
TYPE (mklAxial_Type) :: Axial
REAL :: eigv
! ----------------------------------------------------
TYPE (PinXS_Type), POINTER, DIMENSION(:,:) :: PinXS

INTEGER :: i, ig, ixy, iz, izf, ng, nxy, nzCMFD, ierr, iter, nNeg(1000)
INTEGER :: itermax = 10
REAL :: AxNSolverBeg, AxNSolverEnd
! ----------------------------------------------------

PinXS => CMFD%PinXS
ng     = CMFD%ng

nxy    = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

CALL SetSourceOperator(Axial, PinXS)
CALL SetFlux(CMFD, Axial)
CALL SetPsi(Axial)
CALL SetLeakage(CMFD, Axial, PinXS)
CALL SetPseudoAbsorption(Axial)
CALL SetCrossSection(Axial, PinXS)
CALL SetConstant(Axial, PinXS)
! ----------------------------------------------------
DO iter = 1, itermax
  AxNSolverBeg = nTracer_dclock(FALSE, FALSE)
  
  CALL SetSource(Axial, eigv)
  CALL SetSourceMoment(Axial)
  CALL SetP1Source(Axial)
  
  phis = ZERO
  phim = ZERO
  
  !$OMP PARALLEL DO SCHEDULE(GUIDED)
  DO ixy = 1, nxy
    CALL FlatRayTraceP0(Axial, ixy, FALSE)
  ENDDO
  !$OMP END PARALLEL DO
  
  CALL SetFluxMoment(Axial)
  
  AxNSolverEnd = nTracer_dclock(FALSE, FALSE)
  TimeChk%AxNSolverTime = TimeChk%AxNSolverTime + (AxNSolverEnd - AxNSolverBeg)
  
  CALL SetBoundaryFlux(Axial)
END DO
! ----------------------------------------------------
phis = ZERO
phim = ZERO
Axial%Jout = ZERO

AxNSolverBeg = nTracer_dclock(FALSE, FALSE)

!$OMP PARALLEL DO SCHEDULE(GUIDED)
DO ixy = 1, nxy
  CALL FlatRayTraceP0(Axial, ixy, TRUE)
ENDDO
!$OMP END PARALLEL DO

CALL SetFluxMoment(Axial)

AxNSolverEnd = nTracer_dclock(FALSE, FALSE)
TimeChk%AxNSolverTime = TimeChk%AxNSolverTime + (AxNSolverEnd - AxNSolverBeg)
! ----------------------------------------------------
CALL SetBoundaryFlux(Axial)
CALL SetCellFlux(Axial)

DO iz = 1, nzMOC
  nNeg(iz) = COUNT(phis(:, iz, :) .LT. ZERO)
END DO

lFirst = FALSE
! ----------------------------------------------------

END SUBROUTINE FlatMOCDriver
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetSubmesh()
! MoC Pln. <= Sub-Pln. <= Pln.

USE allocs

IMPLICIT NONE

INTEGER :: nzCMFD, iz, izc, izf, myzb, myze
REAL, POINTER, DIMENSION(:) :: hzfm
! ----------------------------------------------------

nzCMFD = mklGeom%nzCMFD
myzb   = mklGeom%myzb
myze   = mklGeom%myze
hzfm  => mklGeom%hzfm

! SET : # of MoC Pln.
nzMOC = 0

DO iz = 1, nzCMFD
  nDiv(iz) = INT(mklGeom%hzfm(iz) / mklCntl%MOCHeight) + 1
  nzMOC = nzMOC + nDiv(iz)
END DO

CALL dmalloc(hzMOC, nzMOC)
CALL dmalloc(cmMap, nzMOC)
CALL dmalloc(fmMap, nzMOC)
CALL dmalloc(fmRange, nzCMFD, 2)
CALL dmalloc0(cmRange, myzb, myze, 1, 2)

! SET : Map of MoC Pln.
izf = 0

DO izc = myzb, myze ! Pln.
  cmRange(izc, 1) = izf + 1 ! MoC Pln St. at Pln.
  
  DO iz = mklGeom%fmRange(izc, 1), mklGeom%fmRange(izc, 2) ! Sub-Pln.
    fmRange(iz, 1) = izf + 1        ! MoC Pln St. at Sub-Pln.
    fmRange(iz, 2) = izf + nDiv(iz) ! MoC Pln Ed. at Sub-Pln.
    
    fmMap(fmRange(iz, 1):fmRange(iz, 2)) = iz ! MoC Pln. to Sub-Pln.
    hzMOC(izf + 1:izf + nDiv(iz)) = hzfm(iz) / nDiv(iz) ! Hgt. of MoC Pln.
    izf = izf + nDiv(iz)
  END DO
  
  cmRange(izc, 2) = izf ! MoC Pln Ed.
  cmMap(cmRange(izc, 1):cmRange(izc, 2)) = izc ! MoC Pln. to Pln.
END DO

NULLIFY (hzfm)
! ----------------------------------------------------

END SUBROUTINE SetSubmesh
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetSphericalHarmonics(Axial)

USE allocs

IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial
TYPE(mklAngle_Type), POINTER, DIMENSION(:) :: Angle

INTEGER :: nPolarAngle, ipol
REAL :: cosv
! ----------------------------------------------------

Angle => Axial%Angle
nPolarAngle = mklGeom%nPolar1D

CALL dmalloc(Comp, 3, nPolarAngle)
CALL dmalloc(mwt,  3, nPolarAngle)

DO ipol = 1, nPolarAngle
  cosv = Angle(ipol)%cosv
  
  Comp(1, ipol) = cosv
  Comp(2, ipol) = 0.5 * (3.0 * cosv ** 2 - 1.0)
  Comp(3, ipol) = 0.5 * (5.0 * cosv ** 3 - 3.0 * cosv)
  
  mwt (:, ipol) = Comp(:, ipol) * Angle(ipol)%wtsurf
END DO
! ----------------------------------------------------

END SUBROUTINE SetSphericalHarmonics
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetSourceOperator(Axial, PinXS)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial
TYPE(PinXS_Type), POINTER :: PinXS(:, :)

INTEGER :: ig, igs, ipin, ipin_map, iz, izf
INTEGER :: ng, nxy, nzCMFD
INTEGER, POINTER :: pinMap(:), planeMap(:)

ng = Axial%ng
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

! SUBROUTINE SetFlux()
! 
! IMPLICIT NONE
! 
! INTEGER :: ig, ipin, iz, izf
! INTEGER :: ng, nxy, nzCMFD
! 
! ng = mklGeom%ng
! nxy = mklGeom%nxy
! nzCMFD = mklGeom%nzCMFD
! 
! !$OMP PARALLEL
! !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
! DO ipin = 1, nxy
!   DO iz = 1, nzCMFD
!     DO ig = 1, ng
!       IF (mklCMFD%phis(ipin, iz, ig) .LT. 0.0) CYCLE
!       DO izf = fmRange(iz, 1), fmRange(iz, 2)
!         phis(ig, izf, ipin) = mklCMFD%phis(ipin, iz, ig)
!         phis(ig, izf, ipin) = phis(ig, izf, ipin) * phisShape(ig, izf, ipin)
!       ENDDO
!     ENDDO
!   ENDDO
! ENDDO
! !$OMP END DO
! !$OMP END PARALLEL
! 
! END SUBROUTINE
!
! SUBROUTINE SetFluxShape()
! 
! IMPLICIT NONE
! 
! INTEGER :: ig, ipin, iz, izf
! INTEGER :: ng, nxy, nzCMFD
! LOGICAL :: lNegative
! 
! ng = mklGeom%ng
! nxy = mklGeom%nxy
! nzCMFD = mklGeom%nzCMFD
! 
! IF (mklCntl%lRefPinFDM) THEN
!   !$OMP PARALLEL
!   !$OMP DO SCHEDULE(GUIDED)
!   DO ipin = 1, nxy
!     IF (.NOT. mklGeom%lRefPin(ipin)) CYCLE
!     DO izf = 1, nzMOC
!       phis(:, izf, ipin) = 1.0
!     ENDDO
!   ENDDO
!   !$OMP END DO
!   !$OMP END PARALLEL
! ENDIF
! 
! mklAxial%phic = 0.0
! 
! !$OMP PARALLEL
! !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
! DO ipin = 1, nxy
!   DO iz = 1, nzCMFD
!     DO izf = fmRange(iz, 1), fmRange(iz, 2)
!       DO ig = 1, ng
!         mklAxial%phic(ig, ipin, iz) = mklAxial%phic(ig, ipin, iz) + phis(ig, izf, ipin)
!       ENDDO
!     ENDDO
!     mklAxial%phic(:, ipin, iz) = mklAxial%phic(:, ipin, iz) / nDiv(iz)
!     DO izf = fmRange(iz, 1), fmRange(iz, 2)
!       DO ig = 1, ng
!         phisShape(ig, izf, ipin) = phis(ig, izf, ipin) / mklAxial%phic(ig, ipin, iz)
!       ENDDO
!     ENDDO
!   ENDDO
! ENDDO
! !$OMP END DO
! !$OMP END PARALLEL
! 
! !$OMP PARALLEL PRIVATE(lNegative)
! !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
! DO ipin = 1, nxy
!   DO iz = 1, nzCMFD
!     DO ig = 1, ng
!       lNegative = .FALSE.
!       DO izf = fmRange(iz, 1), fmRange(iz, 2)
!         IF (phis(ig, izf, ipin) .LT. 0.0) lNegative = .TRUE.
!       ENDDO
!       IF (lNegative) THEN
!         DO izf = fmRange(iz, 1), fmRange(iz, 2)
!           phisShape(ig, izf, ipin) = 1.0
!         ENDDO
!       ENDIF
!     ENDDO
!   ENDDO
! ENDDO
! !$OMP END DO
! !$OMP END PARALLEL
! 
! END SUBROUTINE

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

SUBROUTINE SetCrossSection(Axial, PinXS)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial
TYPE(PinXS_Type), POINTER :: PinXS(:, :)

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

SUBROUTINE SetPsi(Axial)

IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial

INTEGER :: ig, igs, ipin, iz, izf
INTEGER :: ng, nxy, nzCMFD

ng = Axial%ng
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

SUBROUTINE SetLeakage(CMFD, Axial, PinXS)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
TYPE(mklAxial_Type) :: Axial
TYPE(PinXS_Type), POINTER :: PinXS(:, :)

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
      DO ibd = 1, Pin(ipin_map)%nNgh
        ineighpin = Pin(ipin_map)%NeighIdx(ibd)
        ineighpin = pinMapRev(ineighpin)
        IF (ineighpin .EQ. VoidCell) THEN
          neighphi = 0.0
        ELSEIF (ineighpin .EQ. RefCell) THEN
          neighphi = myphi
        ELSE
          neighphi = CMFD%phis(ineighpin, iz, ig)
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

SUBROUTINE SetConstant(Axial, PinXS)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial
TYPE(PinXS_Type), POINTER :: PinXS(:, :)

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

SUBROUTINE SetSourceMoment(Axial)

IMPLICIT NONE

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
        gb = mklGeom%InScatRange(1, ig)
        ge = mklGeom%InScatRange(2, ig)
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

SUBROUTINE SetSource(Axial, eigv)

IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial
REAL :: eigv, reigv

TYPE(mklAngle_Type), POINTER :: Angle(:)
INTEGER :: ipol, ig, igs, ipin, iz, izf
INTEGER :: ng, nxy, nzCMFD, nPolarAngle
INTEGER :: gb, ge
REAL :: tau

Angle => Axial%Angle
ng = Axial%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
nPolarAngle = mklGeom%nPolar1D

reigv = 1.0 / eigv

src = 0.0

!$OMP PARALLEL PRIVATE(gb, ge)
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