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
SUBROUTINE FlatMOCDriver(CMFD, Axial, eigv, itermax)

USE param,   ONLY : ZERO, FALSE, TRUE
USE TYPEDEF, ONLY : PinXS_Type
USE PE_MOD,  ONLY : PE
USE TIMER,   ONLY : nTracer_dclock, TimeChk

IMPLICIT NONE

TYPE (mklCMFD_Type)  :: CMFD
TYPE (mklAxial_Type) :: Axial
REAL :: eigv
INTEGER :: itermax
! ----------------------------------------------------
TYPE (PinXS_Type), POINTER, DIMENSION(:,:) :: PinXS

INTEGER :: i, ig, ixy, iz, izf, ng, nxy, nzCMFD, ierr, iter, nNgt(1000)
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
  
  CALL SetBndyFlux(Axial)
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
CALL SetBndyFlux(Axial)
CALL SetCellFlux(Axial)

DO iz = 1, nzMOC
  nNgt(iz) = COUNT(phis(:, iz, :) .LT. ZERO) ! DEBUG
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

USE TYPEDEF, ONLY : PinXS_Type

IMPLICIT NONE

TYPE (mklAxial_Type) :: Axial
TYPE (PinXS_Type), POINTER, DIMENSION(:,:) :: PinXS

INTEGER :: ig, jg, ipin, ipin_map, iz, jz, ng, nxy, nzCMFD
INTEGER, POINTER, DIMENSION(:) :: pinMap, planeMap
! ----------------------------------------------------

ng = Axial%ng

nxy       = mklGeom%nxy
nzCMFD    = mklGeom%nzCMFD
pinMap   => mklGeom%pinMap
planeMap => mklGeom%planeMap

!$OMP PARALLEL PRIVATE(jz, ipin_map)
!$OMP DO SCHEDULE(GUIDED)
DO ipin = 1, nxy
  ipin_map = pinMap(ipin)
  
  DO iz = 1, nzCMFD ! Coarse Pln.
    jz = planeMap(iz)
    
    DO ig = 1, ng
      F  (ig, iz, ipin) = PinXS(ipin_map, jz)%XSnf(ig)
      Chi(ig, iz, ipin) = PinXS(ipin_map, jz)%Chi (ig)
      
      DO jg = 1, ng ! From
        IF (PinXS(ipin_map, jz)%XSs(ig)%ib .GT. jg) CYCLE
        IF (PinXS(ipin_map, jz)%XSs(ig)%ie .LT. jg) CYCLE
        
        IF (jg .EQ. ig) THEN
          S(jg, ig, iz, ipin) = PinXS(ipin_map, jz)%XSs(ig)%self
        ELSE
          S(jg, ig, iz, ipin) = PinXS(ipin_map, jz)%XSs(ig)%from(jg)
        END IF
      END DO
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
! ----------------------------------------------------

END SUBROUTINE SetSourceOperator
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetFlux(CMFD, Axial)

IMPLICIT NONE

TYPE (mklCMFD_Type)  :: CMFD
TYPE (mklAxial_Type) :: Axial

INTEGER :: ig, ipin, iz, izf, ng, nxy, nzCMFD
REAL :: fmult
! ----------------------------------------------------

ng = Axial%ng

nxy    = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

IF (lFirst) THEN
  !$OMP PARALLEL
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
  DO ipin = 1, nxy
    DO iz = 1, nzCMFD ! Coarse Pln.
      DO ig = 1, ng
        DO izf = fmRange(iz, 1), fmRange(iz, 2)
          phis(ig, izf, ipin) = CMFD%phis(ipin, iz, ig)
        END DO
      END DO
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
ELSE
  !$OMP PARALLEL PRIVATE(fmult)
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
  DO ipin = 1, nxy
    DO iz = 1, nzCMFD ! Coarse Pln.
      DO ig = 1, ng
        fmult = CMFD%phis(ipin, iz, ig) / Axial%phic(ig, ipin, iz)
        
        DO izf = fmRange(iz, 1), fmRange(iz, 2) ! Fine Pln.
          phis(ig, izf, ipin) = fmult * phis(ig, izf, ipin)
        END DO
      END DO
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
END IF
! ----------------------------------------------------

END SUBROUTINE SetFlux
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetCellFlux(Axial)

USE param, ONLY : ZERO

IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial

INTEGER :: ig, ipin, iz, izf, ng, nxy, nzCMFD
LOGICAL :: lNegative
! ----------------------------------------------------

ng = Axial%ng

nxy    = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

Axial%phic = ZERO

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD ! Coarse Pln.
    DO izf = fmRange(iz, 1), fmRange(iz, 2) ! Fine Pln.
      DO ig = 1, ng
        Axial%phic(ig, ipin, iz) = Axial%phic(ig, ipin, iz) + phis(ig, izf, ipin)
      END DO
    END DO
    
    Axial%phic(:, ipin, iz) = Axial%phic(:, ipin, iz) / nDiv(iz)
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
! ----------------------------------------------------

END SUBROUTINE SetCellFlux
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetFluxMoment(Axial)

USE param, ONLY : RTHREE, RFIVE, RSEVEN

IMPLICIT NONE

TYPE (mklAxial_Type) :: Axial
! ----------------------------------------------------
TYPE (mklAngle_Type), POINTER, DIMENSION(:) :: Angle

INTEGER :: ipol, ig, ipin, izf, ng, nxy, nPolarAngle
REAL :: rsigv
! ----------------------------------------------------

Angle => Axial%Angle
ng     = Axial%ng

nxy         = mklGeom%nxy
nPolarAngle = mklGeom%nPolar1D

!$OMP PARALLEL PRIVATE(rsigv)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ipin = 1, nxy
  DO izf = 1, nzMOC ! Fine Pln.
    DO ig = 1, ng
      rsigv = 1.0 / hzMOC(izf) / xst(ig, izf, ipin)
      
      phis(   ig, izf, ipin) = phis(   ig, izf, ipin) * rsigv + src(    ig, izf, ipin)
      phim(1, ig, izf, ipin) = phim(1, ig, izf, ipin) * rsigv + srcm(1, ig, izf, ipin) * RTHREE
!      phim(2, ig, izf, ipin) = phim(2, ig, izf, ipin) * rsigv + srcm(2, ig, izf, ipin) * RFIVE
!      phim(3, ig, izf, ipin) = phim(3, ig, izf, ipin) * rsigv + srcm(3, ig, izf, ipin) * RSEVEN
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
! ----------------------------------------------------

END SUBROUTINE SetFluxMoment
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetCrossSection(Axial, PinXS)

USE TYPEDEF, ONLY : PinXS_Type

IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial
TYPE(PinXS_Type), POINTER, DIMENSION(:,:) :: PinXS

INTEGER :: ig, ipin, ipin_map, iz, jz, izf, ng, nxy, nzCMFD
INTEGER, POINTER, DIMENSION(:) :: pinMap, planeMap
! ----------------------------------------------------

ng = Axial%ng

nxy       = mklGeom%nxy
nzCMFD    = mklGeom%nzCMFD
pinMap   => mklGeom%pinMap
planeMap => mklGeom%planeMap

!$OMP PARALLEL PRIVATE(ipin_map, jz)
!$OMP DO SCHEDULE(GUIDED)
DO ipin = 1, nxy
  ipin_map = pinMap(ipin)
  
  DO iz = 1, nzCMFD ! Coarse Pln.
    jz = planeMap(iz)
    
    IF (mklGeom%lH2OCell(iz, ipin)) THEN
      DO ig = 1, ng
        S(ig, ig, iz, ipin) = S(ig, ig, iz, ipin) + (PinXS(ipin_map, jz)%XSt(ig) - PinXS(ipin_map, jz)%XStr(ig))
        
        DO izf = fmRange(iz, 1), fmRange(iz, 2)
          xst(ig, izf, ipin) = pxs(ig, izf, ipin) + PinXS(ipin_map, jz)%XSt(ig)
        END DO
      END DO
    ELSE
      DO ig = 1, ng
        DO izf = fmRange(iz, 1), fmRange(iz, 2)
          xst(ig, izf, ipin) = pxs(ig, izf, ipin) + PinXS(ipin_map, jz)%XStr(ig)
        END DO
      END DO
    END IF
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
! ----------------------------------------------------

END SUBROUTINE SetCrossSection
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetPsi(Axial)

USE param, ONLY : ZERO

IMPLICIT NONE

TYPE (mklAxial_Type) :: Axial

INTEGER :: ig, igs, ipin, iz, izf, ng, nxy, nzCMFD
! ----------------------------------------------------

ng = Axial%ng

nxy    = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

psi = ZERO

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD ! Coarse Pln.
    DO izf = fmRange(iz, 1), fmRange(iz, 2) ! Fine Pln.
      DO ig = 1, ng
        psi(izf, ipin) = psi(izf, ipin) + F(ig, iz, ipin) * phis(ig, izf, ipin)
      END DO
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
! ----------------------------------------------------

END SUBROUTINE SetPsi
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetLeakage(CMFD, Axial, PinXS)

USE allocs
USE param,   ONLY : ZERO
USE TYPEDEF, ONLY : PinXS_Type

IMPLICIT NONE

TYPE (mklCMFD_Type) :: CMFD
TYPE (mklAxial_Type) :: Axial
TYPE (PinXS_Type), POINTER, DIMENSION(:,:) :: PinXS
! ----------------------------------------------------
TYPE(superPin_Type), POINTER, DIMENSION(:) :: Pin

INTEGER :: ig, ibd, ipin, inghpin, ipin_map, iz, izc, izf, ng, nxy, nzCMFD
INTEGER, POINTER, DIMENSION(:) :: pinMap, pinMapRev, planeMap

REAL :: Dtil, Dhat, myphi, nghphi
REAL, POINTER, DIMENSION(:,:,:) :: radLkg
! ----------------------------------------------------

ng = Axial%ng

nxy        = mklGeom%nxy
nzCMFD     = mklGeom%nzCMFD
Pin       => mklGeom%superPin
pinMap    => mklGeom%pinMap
pinMapRev => mklGeom%pinMapRev
planeMap  => mklGeom%planeMap

CALL dmalloc0(radLkg, 1, ng, 1, nxy, 0, nzCMFD+1)

DO iz = 1, nzCMFD ! Coarse Pln.
  izc = planeMap(iz)
  
  !$OMP PARALLEL PRIVATE(ipin_map, inghpin, myphi, nghphi, Dtil, Dhat)
  !$OMP DO SCHEDULE(GUIDED)
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    
    DO ig = 1, ng
      myphi = CMFD%phis(ipin, iz, ig)
      
      DO ibd = 1, Pin(ipin_map)%nNgh
        inghpin = Pin(ipin_map)%NeighIdx(ibd)
        inghpin = pinMapRev(inghpin)
        
        IF (inghpin .EQ. VOIDCELL) THEN
          nghphi = ZERO
        ELSE IF (inghpin .EQ. REFCELL) THEN
          nghphi = myphi
        ELSE
          nghphi = CMFD%phis(inghpin, iz, ig)
        END IF
        
        Dtil = PinXS(ipin_map, izc)%Dtil(ibd, ig)
        Dhat = PinXS(ipin_map, izc)%Dhat(ibd, ig)
        
        radLkg(ig, ipin, iz) = radLkg(ig, ipin, iz) - Dtil * (nghphi - myphi) - Dhat * (nghphi + myphi)
      END DO
      
      radLkg(ig, ipin, iz) = radLkg(ig, ipin, iz) * mklGeom%hzfm(iz) / mklGeom%PinVolFm(ipin, iz)
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
END DO
! ----------------------------------------------------
CALL InitFastComm()
CALL GetNeighborFast(ng * nxy, radLkg(:, :,      1), radLkg(:, :, nzCMFD+1), BOTTOM)
CALL GetNeighborFast(ng * nxy, radLkg(:, :, nzCMFD), radLkg(:, :,        0), TOP)
CALL FinalizeFastComm()

IF (mklGeom%lBottom) THEN
  IF (mklGeom%AxBC(BOTTOM) .EQ. VOIDCELL) radLkg(:, :, 0) = ZERO
  IF (mklGeom%AxBC(BOTTOM) .EQ. REFCELL) CALL dcopy(ng * nxy, radLkg(:, :, 1), 1, radLkg(:, :, 0), 1)
ENDIF

IF (mklGeom%lTop) THEN
  IF (mklGeom%AxBC(TOP) .EQ. VOIDCELL) radLkg(:, :, nzCMFD+1) = ZERO
  IF (mklGeom%AxBC(TOP) .EQ. REFCELL) CALL dcopy(ng * nxy, radLkg(:, :, nzCMFD), 1, radLkg(:, :, nzCMFD+1), 1)
ENDIF
! ----------------------------------------------------
DO iz = 1, nzCMFD
  CALL LeakageExpansion(radLkg(:, :, iz-1), radLkg(:, :, iz), radLkg(:, :, iz+1), iz, ng, nxy)
ENDDO

DEALLOCATE (radLkg)
! ----------------------------------------------------

END SUBROUTINE SetLeakage
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE LeakageExpansion(L0, L1, L2, iz, ng, nxy)

IMPLICIT NONE

REAL, DIMENSION(:,:) :: L0, L1, L2
INTEGER :: iz, ng, nxy
! ----------------------------------------------------
REAL, DIMENSION(3) :: n0, n1, n2
REAL :: d0, d1, d2, h0, h1, h2, dh, x0, x1, a, b, c
INTEGER :: ig, ipin, izf
! ----------------------------------------------------

h0 = mklGeom%hzfm(iz - 1)
h1 = mklGeom%hzfm(iz)
h2 = mklGeom%hzfm(iz + 1)
dh = h1 / nDiv(iz)

n0(1) = h1 ** 3 + 2.0 * h1 ** 2 * h2 + h1 * h2 ** 2
n0(2) = 2.0 * h0 ** 2 * h1 + 3.0 * h0 * h1 ** 2 + h0 ** 2 * h2 + 3.0 * h0 * h1 * h2 + h0 * h2 ** 2
n0(3) = -h0 ** 2 * h1 - h0 * h1 ** 2

n1(1) = 2.0 * h1 ** 2 + 3.0 * h1 * h2 + h2 ** 2
n1(2) = h0 ** 2 - 3.0 * h1 ** 2 - 3.0 * h1 * h2 - h2 ** 2
n1(3) = -h0 ** 2 + h1 ** 2

n2(1) = h1 + h2
n2(2) = -h0 - 2.0 * h1 - h2
n2(3) = h0 + h1

d0 = (h1 + h2) * (h0 ** 2 + 2.0 * h0 * h1 + h1 ** 2 + h0 * h2 + h1 * h2)
d1 = (h0 + h1) * (h1 + h2) * (h0 + h1 + h2)
d2 = (h0 + h1) * (h1 + h2) * (h0 + h1 + h2)

!$OMP PARALLEL PRIVATE(a, b, c, x0, x1)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO ig = 1, ng
    a =   3.0 * (L0(ig, ipin) * n2(1) + L1(ig, ipin) * n2(2) + L2(ig, ipin) * n2(3)) / d2
    b = - 2.0 * (L0(ig, ipin) * n1(1) + L1(ig, ipin) * n1(2) + L2(ig, ipin) * n1(3)) / d1
    c = (L0(ig, ipin) * n0(1) + L1(ig, ipin) * n0(2) + L2(ig, ipin) * n0(3)) / d0
    
    DO izf = fmRange(iz, 1), fmRange(iz, 2) ! Fine Pln.
      x0 = (izf - fmRange(iz, 1)) * dh
      x1 = x0 + dh
      
      lkg(ig, izf, ipin) = Integral(a, b, c, x0, x1) / dh
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

CONTAINS
! ----------------------------------------------------
FUNCTION Integral(a, b, c, x0, x1) RESULT(val)

IMPLICIT NONE

REAL :: a, b, c, x0, x1
REAL :: val

val = a * (x1 ** 3 - x0 ** 3) / 3.0 + b * (x1 ** 2 - x0 ** 2) / 2.0 + c * (x1 - x0)

END FUNCTION Integral
! ----------------------------------------------------

END SUBROUTINE LeakageExpansion
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetPseudoAbsorption(Axial)

USE param, ONLY : ZERO

IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial

INTEGER :: ig, ipin, iz, izf, ng, nxy, nzCMFD
! ----------------------------------------------------

ng = Axial%ng

nxy    = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

pxs = ZERO

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD ! Coarse Pln.
    IF (.NOT. mklGeom%lREFCELL(iz, ipin)) CYCLE
    
    DO izf = fmRange(iz, 1), fmRange(iz, 2) ! Fine Pln.
      DO ig = 1, ng
        IF (lkg (ig, izf, ipin) .LT. ZERO) CYCLE
        IF (phis(ig, izf, ipin) .LT. ZERO) CYCLE
        
        pxs(ig, izf, ipin) = lkg(ig, izf, ipin) / phis(ig, izf, ipin)
        lkg(ig, izf, ipin) = ZERO
      END DO
    END DO
    
    IF (mklGeom%lH2OCell(iz, ipin)) CYCLE
    
    DO ig = 1, ng
      IF (S(ig, ig, iz, ipin) .GT. ZERO) CYCLE
      
      DO izf = fmRange(iz, 1), fmRange(iz, 2) ! Fine Pln.
        pxs(ig, izf, ipin) = pxs(ig, izf, ipin) - S(ig, ig, iz, ipin)
      END DO
      
      S(ig, ig, iz, ipin) = ZERO
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
! ----------------------------------------------------

END SUBROUTINE SetPseudoAbsorption
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetConstant(Axial, PinXS)

USE PARAM,   ONLY : ONE
USE TYPEDEF, ONLY : PinXS_Type

IMPLICIT NONE

TYPE (mklAxial_Type) :: Axial
TYPE (PinXS_Type), POINTER, DIMENSION(:,:) :: PinXS
! ----------------------------------------------------
TYPE (mklAngle_Type), POINTER, DIMENSION(:) :: Angle

INTEGER :: ipol, ig, ipin, iz, izf, ng, nxy, nPol
INTEGER, POINTER, DIMENSION(:) :: pinMap
REAL :: tau
! ----------------------------------------------------

ng     = Axial%ng
Angle => Axial%Angle

nxy     = mklGeom%nxy
nPol    = mklGeom%nPolar1D
pinMap => mklGeom%pinMap


!$OMP PARALLEL PRIVATE(tau)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(4)
DO ipin = 1, nxy
  DO izf = 1, nzMOC ! Fine Pln.
    DO ig = 1, ng
      DO ipol = 1, nPol
        tau = xst(ig, izf, ipin) * hzMOC(izf) * Angle(ipol)%rcosv
        
        wtExp(ipol, ig, izf, ipin) = ONE - EXP(-tau)
      END DO
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
! ----------------------------------------------------

END SUBROUTINE SetConstant
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetSourceMoment(Axial)

USE param, ONLY : ZERO

IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial

INTEGER :: ig, jg, ipin, iz, izf, ng, nxy, nzCMFD, gb, ge
! ----------------------------------------------------

ng = Axial%ng

nxy    = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

srcm = ZERO

!$OMP PARALLEL PRIVATE(gb, ge)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD ! Coarse Pln.
    IF (.NOT. mklGeom%lH2OCell(iz, ipin)) CYCLE
    
    DO izf = fmRange(iz, 1), fmRange(iz, 2) ! Fine Pln.
      DO ig = 1, ng
        gb = mklGeom%InScatRange(1, ig)
        ge = mklGeom%InScatRange(2, ig)
        
        DO jg = gb, ge
          srcm(1, ig, izf, ipin) = srcm(1, ig, izf, ipin) + Axial%SmP1(jg, ig, iz, ipin) * phim(1, jg, izf, ipin)
          srcm(2, ig, izf, ipin) = srcm(2, ig, izf, ipin) + Axial%SmP2(jg, ig, iz, ipin) * phim(2, jg, izf, ipin)
          srcm(3, ig, izf, ipin) = srcm(3, ig, izf, ipin) + Axial%SmP3(jg, ig, iz, ipin) * phim(3, jg, izf, ipin)
        END DO
        
        srcm(1, ig, izf, ipin) = 3.0 * srcm(1, ig, izf, ipin) / xst(ig, izf, ipin)
        srcm(2, ig, izf, ipin) = 5.0 *srcm(2, ig, izf, ipin) / xst(ig, izf, ipin)
        srcm(3, ig, izf, ipin) = 7.0 * srcm(3, ig, izf, ipin) / xst(ig, izf, ipin)
      END DO
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
! ----------------------------------------------------

END SUBROUTINE SetSourceMoment
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetSource(Axial, eigv)

USE param, ONLY : ZERO, ONE

IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial
REAL :: eigv
! ----------------------------------------------------
TYPE(mklAngle_Type), POINTER, DIMENSION(:) :: Angle
INTEGER :: ipol, ig, jg, ipin, iz, izf, ng, nxy, nzCMFD, nPolarAngle, gb, ge
REAL :: tau, reigv
! ----------------------------------------------------

Angle => Axial%Angle
ng     = Axial%ng

nxy         = mklGeom%nxy
nzCMFD      = mklGeom%nzCMFD
nPolarAngle = mklGeom%nPolar1D

reigv = ONE / eigv
src   = ZERO

!$OMP PARALLEL PRIVATE(gb, ge)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO iz = 1, nzCMFD ! Coarse Pln.
    DO izf = fmRange(iz, 1), fmRange(iz, 2) ! Fine Pln.
      DO ig = 1, ng
        ! Sct.
        gb = mklGeom%InScatRange(1, ig)
        ge = mklGeom%InScatRange(2, ig)
        
        DO jg = gb, ge
          src(ig, izf, ipin) = src(ig, izf, ipin) + S(jg, ig, iz, ipin) * phis(jg, izf, ipin)
        END DO
        
        ! Fis.
        src(ig, izf, ipin) = src(ig, izf, ipin) + reigv * Chi(ig, iz, ipin) * psi(izf, ipin)
      END DO
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

CALL vdsub(ng * nxy * nzMOC, src, lkg, src)
CALL vddiv(ng * nxy * nzMOC, src, xst, src)
! ----------------------------------------------------

END SUBROUTINE SetSource
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetP1Source(Axial)

IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial

INTEGER :: ipol, ig, ipin, izf, ng, nxy, nzCMFD, nPol
! ----------------------------------------------------

ng = Axial%ng

nxy    = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
nPol   = mklGeom%nPolar1D

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(4)
DO ipin = 1, nxy
  DO izf = 1, nzMOC ! Fine Pln.
    DO ig = 1, ng
      DO ipol = 1, nPol
        srca(1, ipol, ig, izf, ipin) = src(ig, izf, ipin) + Comp(1, ipol) * srcm(1, ig, izf, ipin)
        srca(1, ipol, ig, izf, ipin) = srca(1, ipol, ig, izf, ipin) + Comp(2, ipol) * srcm(2, ig, izf, ipin)
        srca(1, ipol, ig, izf, ipin) = srca(1, ipol, ig, izf, ipin) + Comp(3, ipol) * srcm(3, ig, izf, ipin)
        
        srca(2, ipol, ig, izf, ipin) = src(ig, izf, ipin) - Comp(1, ipol) * srcm(1, ig, izf, ipin)
        srca(2, ipol, ig, izf, ipin) = srca(2, ipol, ig, izf, ipin) + Comp(2, ipol) * srcm(2, ig, izf, ipin)
        srca(2, ipol, ig, izf, ipin) = srca(2, ipol, ig, izf, ipin) - Comp(3, ipol) * srcm(3, ig, izf, ipin)
      END DO
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
! ----------------------------------------------------

END SUBROUTINE SetP1Source
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetBndyFlux(Axial)

USE allocs
USE param, ONLY : ZERO

IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial

INTEGER :: n, ng, nxy, nPol, ipol, ig, ipin, iz
REAL, POINTER, DIMENSION(:,:,:,:) :: PhiAngIn
! ----------------------------------------------------

ng = Axial%ng
nxy  = mklGeom%nxy
nPol = mklGeom%nPolar1D
n = nPol * ng * nxy

CALL dmalloc(PhiAngIn, nPol, ng, nxy, 2)

CALL InitFastComm
CALL GetNeighborFast(n, Axial%PhiAngOut(:, :, :, BOTTOM), PhiAngIn(:, :, :,    TOP), BOTTOM)
CALL GetNeighborFast(n, Axial%PhiAngOut(:, :, :,    TOP), PhiAngIn(:, :, :, BOTTOM),    TOP)
CALL FinalizeFastComm

IF (mklGeom%lBottom) THEN
  IF (mklGeom%AxBC(BOTTOM) .EQ. VOIDCELL) THEN
    PhiAngIn(:, :, :, BOTTOM) = ZERO
  ELSE IF (mklGeom%AxBC(BOTTOM) .EQ. REFCELL) THEN
    CALL dcopy(n, Axial%PhiAngOut(:, :, :, BOTTOM), 1, PhiAngIn(:, :, :, BOTTOM), 1)
  END IF
END IF

IF (mklGeom%lTop) THEN
  IF (mklGeom%AxBC(TOP) .EQ. VOIDCELL) THEN
    PhiAngIn(:, :, :, TOP) = ZERO
  ELSE IF (mklGeom%AxBC(TOP) .EQ. REFCELL) THEN
    CALL dcopy(n, Axial%PhiAngOut(:, :, :, TOP), 1, PhiAngIn(:, :, :, TOP), 1)
  END IF
END IF

CALL dcopy(n * 2, PhiAngIn, 1, Axial%PhiAngIn, 1)

DEALLOCATE (PhiAngIn)
! ----------------------------------------------------

END SUBROUTINE SetBndyFlux
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE FlatRayTraceP0(Axial, ipin, lJout)

IMPLICIT NONE

TYPE(mklAxial_Type) :: Axial
INTEGER :: ipin
LOGICAL :: lJout
! ----------------------------------------------------
TYPE(mklAngle_Type), POINTER :: Angle(:)

INTEGER :: ig, ipol, iz, izf, ng, nxy, nzCMFD, nPol

REAL, POINTER, DIMENSION(:,:,:,:) :: Jout, PhiAngIn, PhiAngOut
REAL, DIMENSION(mklGeom%nPolar1D, Axial%ng) :: track_phi
REAL :: del_phi
! ----------------------------------------------------

IF (mklCntl%lRefPinFDM) THEN
  IF (mklGeom%lRefPin(ipin)) RETURN
END IF

ng         = Axial%ng
Angle     => Axial%Angle
PhiAngIn  => Axial%PhiAngIn
PhiAngOut => Axial%PhiAngOut
Jout      => Axial%Jout(:, :, :, :, ipin)

nxy    = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
nPol   = mklGeom%nPolar1D

! Upward
track_phi = PhiAngIn(:, :, ipin, BOTTOM)

DO iz = 1, nzCMFD ! Coarse Pln.
  IF (lJout) THEN
    DO ig = 1, ng
      DO ipol = 1, nPol
        Jout(in, ig, BOTTOM, iz) = Jout(in, ig, BOTTOM, iz) + Angle(ipol)%wtsurf * track_phi(ipol, ig)
      END DO
    END DO
  END IF
  
  IF (mklGeom%lH2OCell(iz, ipin)) THEN
    DO izf = fmRange(iz, 1), fmRange(iz, 2) ! Fine Pln.
      DO ig = 1, ng
        DO ipol = 1, nPol
          del_phi = (srca(1, ipol, ig, izf, ipin) - track_phi(ipol, ig)) * wtExp(ipol, ig, izf, ipin)
          
          phis(   ig, izf, ipin) = phis(   ig, izf, ipin) - del_phi * Angle(ipol)%wtsurf
          phim(1, ig, izf, ipin) = phim(1, ig, izf, ipin) - del_phi * mwt(1, ipol)
          phim(2, ig, izf, ipin) = phim(2, ig, izf, ipin) - del_phi * mwt(2, ipol)
          phim(3, ig, izf, ipin) = phim(3, ig, izf, ipin) - del_phi * mwt(3, ipol)
          
          track_phi(ipol, ig) = track_phi(ipol, ig) + del_phi
        END DO
      END DO
    END DO
  ELSE
    DO izf = fmRange(iz, 1), fmRange(iz, 2) ! Fine Pln.
      DO ig = 1, ng
        DO ipol = 1, nPol
          del_phi = (srca(1, ipol, ig, izf, ipin) - track_phi(ipol, ig)) * wtExp(ipol, ig, izf, ipin)
          
          phis(ig, izf, ipin) = phis(ig, izf, ipin) - del_phi * Angle(ipol)%wtsurf
          
          track_phi(ipol, ig) = track_phi(ipol, ig) + del_phi
        END DO
      END DO
    END DO
  END IF
  
  IF (lJout) THEN
    DO ig = 1, ng
      DO ipol = 1, nPol
        Jout(out, ig, TOP, iz) = Jout(out, ig, TOP, iz) + Angle(ipol)%wtsurf * track_phi(ipol, ig)
      END DO
    END DO
  END IF
END DO

PhiAngOut(:, :, ipin, TOP) = track_phi

! Downward
track_phi = PhiAngIn(:, :, ipin, TOP)

DO iz = nzCMFD, 1, -1 ! Coarse Pln.
  IF (lJout) THEN
    DO ig = 1, ng
      DO ipol = 1, nPol
        Jout(in, ig, TOP, iz) = Jout(in, ig, TOP, iz) + Angle(ipol)%wtsurf * track_phi(ipol, ig)
      END DO
    END DO
  END IF
  
  IF (mklGeom%lH2OCell(iz, ipin)) THEN
    DO izf = fmRange(iz, 2), fmRange(iz, 1), -1 ! Fine Pln.
      DO ig = 1, ng
        DO ipol = 1, nPol
          del_phi = (srca(2, ipol, ig, izf, ipin) - track_phi(ipol, ig)) * wtExp(ipol, ig, izf, ipin)
          
          phis(   ig, izf, ipin) = phis(   ig, izf, ipin) - del_phi * Angle(ipol)%wtsurf
          phim(1, ig, izf, ipin) = phim(1, ig, izf, ipin) + del_phi * mwt(1, ipol)
          phim(2, ig, izf, ipin) = phim(2, ig, izf, ipin) - del_phi * mwt(2, ipol)
          phim(3, ig, izf, ipin) = phim(3, ig, izf, ipin) + del_phi * mwt(3, ipol)
          
          track_phi(ipol, ig) = track_phi(ipol, ig) + del_phi
        END DO
      END DO
    END DO
  ELSE
    DO izf = fmRange(iz, 2), fmRange(iz, 1), -1 ! Fine Pln.
      DO ig = 1, ng
        DO ipol = 1, nPol
          del_phi = (srca(2, ipol, ig, izf, ipin) - track_phi(ipol, ig)) * wtExp(ipol, ig, izf, ipin)
          
          phis(ig, izf, ipin) = phis(ig, izf, ipin) - del_phi * Angle(ipol)%wtsurf
          
          track_phi(ipol, ig) = track_phi(ipol, ig) + del_phi
        END DO
      END DO
    END DO
  END IF
  
  IF (lJout) THEN
    DO ig = 1, ng
      DO ipol = 1, nPol
        Jout(out, ig, BOTTOM, iz) = Jout(out, ig, BOTTOM, iz) + Angle(ipol)%wtsurf * track_phi(ipol, ig)
      END DO
    END DO
  END IF
END DO

PhiAngOut(:, :, ipin, BOTTOM) = track_phi

END SUBROUTINE FlatRayTraceP0
! ------------------------------------------------------------------------------------------------------------

END MODULE MKL_FLATMOC
#endif