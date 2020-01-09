#include <defines.h>
!--- CNJ Edit : 1D Axial Polynomial MOC Modules with Intel MKL
#ifdef __INTEL_MKL

MODULE MKL_POLYMOC

USE MKL_3D
IMPLICIT NONE

REAL, POINTER, PRIVATE :: C0(:, :, :, :, :), C1(:, :, :, :, :), EX(:, :, :, :)
REAL, POINTER, PRIVATE :: phisCoeff(:, :, :, :), srcCoeff(:, :, :, :)
REAL, POINTER, PRIVATE :: phis(:, :, :, :), phia(:, :, :, :, :), phim(:, :, :, :)
REAL, POINTER, PRIVATE :: src(:, :, :, :), srcm(:, :, :, :), srcAng(:, :, :, :, :)
REAL, POINTER, PRIVATE :: Comp(:, :), mwt(:, :)

PRIVATE
PUBLIC :: AllocPolynomialMOC, PolynomialMOCDriver

CONTAINS

!--- Public Routines ------------------------------------------------------------------------------

SUBROUTINE AllocPolynomialMOC()

IMPLICIT NONE

INTEGER :: ng, nxy, nzCMFD, nPolar1D
INTEGER :: ScatOd, PolyOd

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
nPolar1D = mklGeom%nPolar1D
ScatOd = mklCntl%scatOrder
PolyOd = mklCntl%polyOrder + 1

ALLOCATE(C0(nPolar1D, ng, nzCMFD, 2, nxy))
ALLOCATE(C1(nPolar1D, ng, nzCMFD, 2, nxy))
ALLOCATE(EX(nPolar1D, ng, nzCMFD, nxy))
ALLOCATE(phisCoeff(nxy, nzCMFD, ng, PolyOd))
ALLOCATE(srcCoeff(nxy, nzCMFD, ng, PolyOd))
  
IF (mklCntl%lScat1) THEN
  ALLOCATE(phia(nPolar1D, ng, nzCMFD, 2, nxy))
  ALLOCATE(phim(nxy, nzCMFD, ng, ScatOd))
  ALLOCATE(srcm(nxy, nzCMFD, ng, ScatOd))
  ALLOCATE(srcAng(nPolar1D, ng, nzCMFD, 2, nxy))
  ALLOCATE(Comp(nPolar1D, ScatOd))
  ALLOCATE(mwt(nPolar1D, ScatOd))
  CALL SetSphericalHarmonics()
ENDIF

END SUBROUTINE

SUBROUTINE PolynomialMOCDriver(PinXS, eigv)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)
REAL :: eigv

REAL :: lkg, rmv, src
INTEGER :: ig, ipin, ipin_map, iz, izf
INTEGER :: ng, nxy, nzCMFD
INTEGER :: iter, itermax = 20
LOGICAL :: lConv

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

CALL SetPolynomialFlux(PinXS)
CALL SetPolynomialSrc(PinXS, eigv)
CALL SetPolynomialRTCoeff(PinXS)
    
lConv = FALSE; iter = 0
IF (mklCntl%lScat1) THEN

  DO WHILE (.NOT. lConv)
    iter = iter + 1
    CALL SetSrcMoment()
    CALL SetPNScatSrc()
    phia = 0.0
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    DO ipin = 1, nxy
      CALL PolynomialRayTracePN(ipin, FALSE)
    ENDDO
    !$OMP END PARALLEL DO
    CALL SetFluxMoment()
    CALL CommBoundaryFluxConv()
    lConv = lConv .OR. (iter .GE. itermax)
  ENDDO
  
  CALL SetSrcMoment()
  CALL SetPNScatSrc()
  phia = 0.0
  mklAxial%Jout = 0.0
  !$OMP PARALLEL DO SCHEDULE(GUIDED)
  DO ipin = 1, nxy
    CALL PolynomialRayTracePN(ipin, TRUE)
  ENDDO
  !$OMP END PARALLEL DO
  CALL SetFluxMoment()
  CALL CommBoundaryFluxConv()
  
ELSE

  DO WHILE (.NOT. lConv)
    iter = iter + 1
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    DO ipin = 1, nxy
      CALL PolynomialRayTraceP0(ipin, FALSE)
    ENDDO
    !$OMP END PARALLEL DO
    CALL CommBoundaryFluxConv(lConv)
    lConv = lConv .OR. (iter .GE. itermax)
  ENDDO
  
  mklAxial%Jout = 0.0
  !$OMP PARALLEL DO SCHEDULE(GUIDED)
  DO ipin = 1, nxy
    CALL PolynomialRayTraceP0(ipin, TRUE)
  ENDDO
  !$OMP END PARALLEL DO
  CALL CommBoundaryFluxConv()
  
ENDIF

END SUBROUTINE

!--- Private Routines -----------------------------------------------------------------------------

SUBROUTINE SetSphericalHarmonics()

IMPLICIT NONE

TYPE(mklAngle_Type), POINTER :: Angle(:)
INTEGER :: nPolar1D
INTEGER :: ipol

Angle => mklAxial%Angle
nPolar1D = mklGeom%nPolar1D

DO ipol = 1, nPolar1D
  Comp(ipol, 1) = Angle(ipol)%cosv
  mwt(ipol, 1) = Comp(ipol, 1) * Angle(ipol)%wt
ENDDO
IF (mklCntl%scatOrder .GE. 2) THEN
  DO ipol = 1, nPolar1D
    Comp(ipol, 2) = 1.5 * Angle(ipol)%cosv ** 2 - 1.0
    mwt(ipol, 2) = Comp(ipol, 2) * Angle(ipol)%wt
  ENDDO
ENDIF
IF (mklCntl%scatOrder .EQ. 3) THEN
  DO ipol = 1, nPolar1D
    Comp(ipol, 3) = 2.5 * Angle(ipol)%cosv ** 3 - 1.5 * Angle(ipol)%cosv
    mwt(ipol, 3) = Comp(ipol, 3) * Angle(ipol)%wt
  ENDDO
ENDIF

END SUBROUTINE

SUBROUTINE QuadraticExpansion(a, b, c, h, y0, y1, yavg)

IMPLICIT NONE

REAL :: a, b, c
REAL :: h, y0, y1, yavg

a = 3.0 * (y0 + y1) / h ** 2 - 6.0 * yavg / h ** 2
b = - (4.0 * y0 + 2.0 * y1) / h + 6.0 * yavg / h
c = y0

END SUBROUTINE

SUBROUTINE QuarticExpansion(Coeff, h, phiL, phiR, phiavg, jL, jR)

IMPLICIT NONE

REAL :: Coeff(5)
REAL :: h, phiL, phiR, phiavg, jL, jR

Coeff(1) = - 5.0 * (h * (jL - jR) - 12.0 * phiavg + 6.0 * (phiL + phiR)) / (2.0 * h ** 4)
Coeff(2) = 2.0 * (h * (3.0 * jL - 2.0 * jR) - 30.0 * phiavg + 16.0 * phiL + 14.0 * phiR) / h ** 3
Coeff(3) = - 3.0 * (h * (3.0 * jL - jR) - 20.0 * phiavg + 12.0 * phiL + 8.0 * phiR) / (2.0 * h ** 2)
Coeff(4) = jL
Coeff(5) = phiL

END SUBROUTINE

SUBROUTINE HyperbolicExpansion(Coeff, h, phiL, phiR, phiavg, jL, jR)

IMPLICIT NONE

REAL :: Coeff(5)
REAL :: h, phiL, phiR, phiavg, jL, jR
REAL :: cosh1, sinh1

cosh1 = dcosh(1.0D0)
sinh1 = dsinh(1.0D0)

Coeff(1) = - (cosh1 * (h * (- jL + jR) + 6.0 * phiavg - 3.0 * (phiL + phiR))                                        &
              + sinh1 * (h * (jL - jR) - 6.0 * phiavg + 3.0 * (phiL + phiR)))                                       &
           / (2.0 * (3.0 * cosh1 - 4.0 * sinh1) * (cosh1 - sinh1))
Coeff(2) = - (h * (jL + jR) + (phiL - phiR)) / (2.0 * (- cosh1 + sinh1))
Coeff(3) = - 3.0 * (cosh1 * h * (jL - jR) + sinh1 * (h * (- jL + jR) - 2.0 * phiavg + phiL + phiR))                 &
           / (4.0 * (3.0 * cosh1 - 4.0 * sinh1))
Coeff(4) = - (cosh1 * (phiL - phiR) + sinh1 * h * (jL + jR)) / (2.0 * (cosh1 - sinh1))
Coeff(5) = - (cosh1 * (h * (- jL + jR) - 12.0 * phiavg)                                                             &
              + sinh1 * (3.0 * h * (jL - jR) + 6.0 * phiavg + 5.0 * (phiL + phiR)))                                 &
           / (4.0 * (3.0 * cosh1 - 4.0 * sinh1))
              
END SUBROUTINE

SUBROUTINE SetPolynomialFlux(PinXS)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)

IF (mklCntl%polyOrder .EQ. 2) THEN
  CALL SetQuadraticFlux()
ELSEIF (mklCntl%polyOrder .EQ. 4) THEN
  IF (mklCntl%lHyperbolic) THEN
    CALL SetHyperbolicFlux(PinXS)
  ELSE
    CALL SetQuarticFlux(PinXS)
  ENDIF
ENDIF

END SUBROUTINE

SUBROUTINE SetQuadraticFlux()

IMPLICIT NONE

INTEGER :: ig, ipin, izf
INTEGER :: ng, nxy, nzCMFD
REAL :: myphi, neighphi, surfphi(2)
REAL :: atil, ahat
REAL :: a, b, c

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

!$OMP PARALLEL PRIVATE(myphi, neighphi, surfphi, atil, ahat, a, b, c)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO izf = 1, nzCMFD
    DO ipin = 1, nxy
      myphi = mklCMFD%phis(ipin, izf, ig)
      !--- Bottom Surface Flux
      IF (izf .EQ. 1) THEN
        neighphi = mklCMFD%neighphis(ipin, ig, bottom)
      ELSE
        neighphi = mklCMFD%phis(ipin, izf - 1, ig)
      ENDIF
      atil = mklAxial%atil(bottom, ipin, izf, ig)
      surfphi(bottom) = atil * myphi + (1.0 - atil) * neighphi + ahat * (myphi + neighphi)
      !--- Top Surface Flux
      IF (izf .EQ. nzCMFD) THEN
        neighphi = mklCMFD%neighphis(ipin, ig, top)
      ELSE
        neighphi = mklCMFD%phis(ipin, izf + 1, ig)
      ENDIF
      atil = mklAxial%atil(top, ipin, izf, ig)
      surfphi(top) = atil * myphi + (1.0 - atil) * neighphi + ahat * (myphi + neighphi)
      !--- Expand Cell Flux
      CALL QuadraticExpansion(a, b, c, mklGeom%hzfm(izf), surfphi(bottom), surfphi(top), myphi)
      phisCoeff(ipin, izf, ig, 1) = a
      phisCoeff(ipin, izf, ig, 2) = b
      phisCoeff(ipin, izf, ig, 3) = c
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetQuarticFlux(PinXS)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)

INTEGER :: ig, ipin, izf
INTEGER :: ng, nxy, nzCMFD
INTEGER, POINTER :: pinMap(:), planeMap(:)
REAL :: myphi, neighphi, jnet(2), surfphi(2)
REAL :: D, Dtil, Dhat, pDhat(2), atil, ahat
REAL :: Coeff(5)

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
pinMap => mklGeom%pinMap
planeMap => mklGeom%planeMap

!$OMP PARALLEL PRIVATE(myphi, neighphi, jnet, surfphi, D, Dtil, Dhat, pDhat, atil, ahat, Coeff)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO izf = 1, nzCMFD
    DO ipin = 1, nxy
      myphi = mklCMFD%phis(ipin, izf, ig)
      D = PinXS(pinMap(ipin), planeMap(izf))%XSD(ig)
      !--- Bottom Surface Flux
      IF (izf .EQ. 1) THEN
        neighphi = mklCMFD%neighphis(ipin, ig, bottom)
      ELSE
        neighphi = mklCMFD%phis(ipin, izf - 1, ig)
      ENDIF
      IF (mklCntl%pCMFD) THEN
        Dtil = mklCMFD%AxDtil(bottom, ipin, izf, ig)
        pDhat = mklCMFD%partialAxDhat(:, bottom, ipin, izf, ig)
        jnet(bottom) = - Dtil * (neighphi - myphi) - (pDhat(in) * neighphi - pDhat(out) * myphi)
      ELSE
        Dtil = mklCMFD%AxDtil(bottom, ipin, izf, ig)
        Dhat = mklCMFD%AxDhat(bottom, ipin, izf, ig)
        jnet(bottom) = - Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
      ENDIF
      atil = mklAxial%atil(bottom, ipin, izf, ig)
      surfphi(bottom) = atil * myphi + (1.0 - atil) * neighphi + ahat * (myphi + neighphi)
      !--- Top Surface Flux
      IF (izf .EQ. nzCMFD) THEN
        neighphi = mklCMFD%neighphis(ipin, ig, top)
      ELSE
        neighphi = mklCMFD%phis(ipin, izf + 1, ig)
      ENDIF
      IF (mklCntl%pCMFD) THEN
        Dtil = mklCMFD%AxDtil(top, ipin, izf, ig)
        pDhat = mklCMFD%partialAxDhat(:, top, ipin, izf, ig)
        jnet(top) = Dtil * (neighphi - myphi) + (pDhat(in) * neighphi - pDhat(out) * myphi)
      ELSE
        Dtil = mklCMFD%AxDtil(top, ipin, izf, ig)
        Dhat = mklCMFD%AxDhat(top, ipin, izf, ig)
        jnet(top) = Dtil * (neighphi - myphi) + Dhat * (neighphi + myphi)
      ENDIF
      atil = mklAxial%atil(top, ipin, izf, ig)
      surfphi(top) = atil * myphi + (1.0 - atil) * neighphi + ahat * (myphi + neighphi)
      !--- Expand Cell Flux
      jnet = jnet / D
      CALL QuarticExpansion(Coeff, mklGeom%hzfm(izf), surfphi(bottom), surfphi(top), myphi, jnet(bottom), jnet(top))
      phisCoeff(ipin, izf, ig, :) = Coeff
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetHyperbolicFlux(PinXS)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)

INTEGER :: ig, ipin, izf
INTEGER :: ng, nxy, nzCMFD
INTEGER, POINTER :: pinMap(:), planeMap(:)
REAL :: myphi, neighphi, jnet(2), surfphi(2)
REAL :: D, Dtil, Dhat, atil, ahat
REAL :: Coeff(5)

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
pinMap => mklGeom%pinMap
planeMap => mklGeom%planeMap

!$OMP PARALLEL PRIVATE(myphi, neighphi, jnet, surfphi, D, Dtil, Dhat, atil, ahat, Coeff)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO izf = 1, nzCMFD
    DO ipin = 1, nxy
      myphi = mklCMFD%phis(ipin, izf, ig)
      D = PinXS(pinMap(ipin), planeMap(izf))%XSD(ig)
      !--- Bottom Surface Flux
      IF (izf .EQ. 1) THEN
        neighphi = mklCMFD%neighphis(ipin, ig, bottom)
      ELSE
        neighphi = mklCMFD%phis(ipin, izf - 1, ig)
      ENDIF
      Dtil = mklCMFD%AxDtil(bottom, ipin, izf, ig)
      Dhat = mklCMFD%AxDhat(bottom, ipin, izf, ig)
      atil = mklAxial%atil(bottom, ipin, izf, ig)
      jnet(bottom) = - Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
      surfphi(bottom) = atil * myphi + (1.0 - atil) * neighphi + ahat * (myphi + neighphi)
      !--- Top Surface Flux
      IF (izf .EQ. nzCMFD) THEN
        neighphi = mklCMFD%neighphis(ipin, ig, top)
      ELSE
        neighphi = mklCMFD%phis(ipin, izf + 1, ig)
      ENDIF
      Dtil = mklCMFD%AxDtil(top, ipin, izf, ig)
      Dhat = mklCMFD%AxDhat(top, ipin, izf, ig)
      atil = mklAxial%atil(top, ipin, izf, ig)
      jnet(top) = Dtil * (neighphi - myphi) + Dhat * (neighphi + myphi)
      surfphi(top) = atil * myphi + (1.0 - atil) * neighphi + ahat * (myphi + neighphi)
      !--- Expand Cell Flux
      jnet = jnet / D
      CALL HyperbolicExpansion(Coeff, 0.5 * mklGeom%hzfm(izf), surfphi(bottom), surfphi(top), myphi, jnet(bottom), jnet(top))
      phisCoeff(ipin, izf, ig, :) = Coeff
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetPolynomialSrc(PinXS, eigv)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)
REAL :: eigv

TYPE(superPin_Type), POINTER :: Pin(:)
INTEGER :: Order
INTEGER :: i, ig, igf, iz, izf, ibd, ipin, ipin_map, ineighpin
INTEGER :: n, ng, nxy, nzCMFD
INTEGER, POINTER :: pinMap(:), pinMapRev(:), planeMap(:)
REAL, POINTER :: scatSrc(:), fisSrc(:), src(:, :), psi(:)
REAL :: Dtil, Dhat, pDhat(2), radLkg, myphi, neighphi

Pin => mklGeom%superPin
ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
pinMap => mklGeom%pinMap
pinMapRev => mklGeom%pinMapRev
planeMap => mklGeom%planeMap
Order = mklCntl%polyOrder + 1

!--- Compute Source Expansion Coefficients

ALLOCATE(scatSrc(nxy * nzCMFD))
ALLOCATE(fisSrc(nxy * nzCMFD))
ALLOCATE(src(nxy * nzCMFD, ng))
ALLOCATE(psi(nxy * nzCMFD))

DO i = 1, Order

  src = 0.0
  psi = 0.0
  
  DO ig = 1, ng
    CALL vdmul(nxy * nzCMFD, mklCMFD%F(:, ig), phisCoeff(:, :, ig, i), fisSrc)
    CALL vdadd(nxy * nzCMFD, fisSrc, psi, psi)
    DO igf = 1, ng
      CALL vdmul(nxy * nzCMFD, mklCMFD%S(:, igf, ig), phisCoeff(:, :, igf, i), scatSrc)
      CALL vdadd(nxy * nzCMFD, scatSrc, src(:, ig), src(:, ig))
    ENDDO
  ENDDO
  
  CALL dscal(nxy * nzCMFD, 1.0 / eigv, psi, 1)
  
  DO ig = 1, ng
    CALL vdmul(nxy * nzCMFD, mklCMFD%Chi(:, ig), psi, fisSrc)
    CALL vdadd(nxy * nzCMFD, src(:, ig), fisSrc, srcCoeff(:, :, ig, i))
    CALL vddiv(nxy * nzCMFD, srcCoeff(:, :, ig, i), mklGeom%PinVolFm, srcCoeff(:, :, ig, i))
  ENDDO
  
ENDDO

!--- Compute Radial Leakage

DO ig = 1, ng
  DO izf = 1, nzCMFD
    iz = planeMap(izf)
    !$OMP PARALLEL PRIVATE(ipin_map, ineighpin, myphi, neighphi, radLkg, Dtil, Dhat, pDhat)
    !$OMP DO SCHEDULE(GUIDED)
    DO ipin = 1, nxy
      ipin_map = pinMap(ipin)
      myphi = mklCMFD%phis(ipin, izf, ig)
      radLkg = 0.0
      DO ibd = 1, 4
        ineighpin = Pin(ipin_map)%NeighIdx(ibd)
        IF (ineighpin .EQ. VoidCell) THEN
          neighphi = 0.0
        ELSEIF (ineighpin .EQ. RefCell) THEN
          neighphi = myphi
        ELSE
          ineighpin = pinMapRev(ineighpin)
          neighphi = mklCMFD%phis(ineighpin, izf, ig)
        ENDIF
        IF (mklCntl%pCMFD) THEN
          Dtil = PinXS(ipin_map, iz)%Dtil(ibd, ig)
          pDhat = PinXS(ipin_map, iz)%partialDhat(:, ibd, ig)
          radLkg = radLkg - Dtil * (neighphi - myphi) - (pDhat(in) * neighphi - pDhat(out) * myphi)
        ELSE
          Dtil = PinXS(ipin_map, iz)%Dtil(ibd, ig)
          Dhat = PinXS(ipin_map, iz)%Dhat(ibd, ig)
          radLkg = radLkg - Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
        ENDIF
      ENDDO
      radLkg = radLkg * mklGeom%hzfm(izf) / mklGeom%PinVolFm(ipin, izf)
      srcCoeff(ipin, izf, ig, Order) = srcCoeff(ipin, izf, ig, Order) - radLkg
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL
  ENDDO
ENDDO

DEALLOCATE(scatSrc)
DEALLOCATE(fisSrc)
DEALLOCATE(src)
DEALLOCATE(psi)

END SUBROUTINE

SUBROUTINE SetFluxMoment()

IMPLICIT NONE

INTEGER :: ng, nxy, nzCMFD, nPolarAngle
INTEGER :: i, ig, ipin, iz, ipol
REAL :: phi_moment, phia_diff, phia_sum

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
nPolarAngle = mklGeom%nPolar1D

!$OMP PARALLEL PRIVATE(phim, phia_diff, phia_sum)
DO i = 1, mklCntl%scatOrder
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
  DO ipin = 1, nxy
    DO iz = 1, nzCMFD
      DO ig = 1, ng
        phi_moment = 0.0
        IF (i .EQ. 2) THEN
          DO ipol = 1, nPolarAngle
            phia_sum = phia(ipol, ig, iz, 1, ipin) + phia(ipol, ig, iz, 2, ipin)
            phi_moment = phi_moment + mwt(ipol, i) * phia_sum
          ENDDO
        ELSE
          DO ipol = 1, nPolarAngle
            phia_diff = phia(ipol, ig, iz, 1, ipin) - phia(ipol, ig, iz, 2, ipin)
            phi_moment = phi_moment + mwt(ipol, i) * phia_diff
          ENDDO
        ENDIF
        phim(ipin, iz, ig, i) = phi_moment
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
ENDDO
!$OMP END PARALLEL

CALL daxpy(ng * nxy * nzCMFD, 1.0D0 / 3.0D0, srcm(:, :, :, 1), 1, phim(:, :, :, 1), 1)

IF (mklCntl%scatOrder .GE. 2) THEN
  CALL daxpy(ng * nxy * nzCMFD, 1.0D0 / 5.0D0, srcm(:, :, :, 2), 1, phim(:, :, :, 2), 1)
ENDIF

IF (mklCntl%scatOrder .EQ. 3) THEN
  CALL daxpy(ng * nxy * nzCMFD, 1.0D0 / 7.0D0, srcm(:, :, :, 3), 1, phim(:, :, :, 3), 1)
ENDIF

END SUBROUTINE

SUBROUTINE SetSrcMoment()

IMPLICIT NONE

INTEGER :: ng, nxy, nzCMFD
INTEGER :: i, ig, igf
REAL, POINTER :: scatSrc(:)

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

ALLOCATE(scatSrc(nxy * nzCMFD))

DO i = 1, mklCntl%scatOrder
  DO ig = 1, ng
    DO igf = 1, ng

      CALL vdadd(nxy * nzCMFD, srcm(:, :, ig, i), scatSrc, srcm(:, :, ig, i))
    ENDDO
  ENDDO
ENDDO

CALL dscal(ng * nxy * nzCMFD, 3.0D0, srcm(:, :, :, 1), 1)

IF (mklCntl%scatOrder .GE. 2) THEN
  CALL dscal(ng * nxy * nzCMFD, 5.0D0, srcm(:, :, :, 2), 1)
ENDIF

IF (mklCntl%scatOrder .EQ. 3) THEN
  CALL dscal(ng * nxy * nzCMFD, 7.0D0, srcm(:, :, :, 3), 1)
ENDIF

DEALLOCATE(scatSrc)

END SUBROUTINE

SUBROUTINE SetPNScatSrc()
  
IMPLICIT NONE

INTEGER :: ng, nxy, nzCMFD, nPolarAngle
INTEGER :: i, ig, iz, ipin, ipol
REAL :: src(2)

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
nPolarAngle = mklGeom%nPolar1D

!$OMP PARALLEL PRIVATE(src)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO iz = 1, nzCMFD
    DO ipin = 1, nxy
      SrcAng = 0.0
      DO i = 1, mklCntl%scatOrder
        IF (i .EQ. 2) THEN
          DO ipol = 1, nPolarAngle
            src(1) = src(1) + Comp(ipol, i) * srcm(ipin, iz, ig, i)
            src(2) = src(2) + Comp(ipol, i) * srcm(ipin, iz, ig, i)
          ENDDO
        ELSE
          DO ipol = 1, nPolarAngle
            src(1) = src(1) + Comp(ipol, i) * srcm(ipin, iz, ig, i)
            src(2) = src(2) - Comp(ipol, i) * srcm(ipin, iz, ig, i)
          ENDDO
        ENDIF
      ENDDO
      SrcAng(ipol, ig, iz, 1, ipin) = src(1)
      SrcAng(ipol, ig, iz, 2, ipin) = src(2)
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetPolynomialRTCoeff(PinXS)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)

IF (mklCntl%polyOrder .EQ. 2) THEN
  CALL SetQuadraticRTCoeff(PinXS)
ELSEIF (mklCntl%polyOrder .EQ. 4) THEN
  IF (mklCntl%lHyperbolic) THEN
    CALL SetHyperbolicRTCoeff(PinXS)
  ELSE
    CALL SetQuarticRTCoeff(PinXS)
  ENDIF
ENDIF

END SUBROUTINE

SUBROUTINE SetQuadraticRTCoeff(PinXS)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)

TYPE(mklAngle_Type), POINTER :: Angle(:)
INTEGER :: ng, nxy, nzCMFD, nPolarAngle
INTEGER :: ig, ipol, ipin, iz
INTEGER, POINTER :: pinMap(:), planeMap(:)
REAL :: a, b, c, sigt(mklGeom%nPolar1D), cosv, tau, h
REAL, POINTER :: hzfm(:)

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
nPolarAngle = mklGeom%nPolar1D
hzfm => mklGeom%hzfm
pinMap => mklGeom%pinMap
planeMap => mklGeom%planeMap
Angle => mklAxial%Angle

!$OMP PARALLEL PRIVATE(a, b, c, sigt, cosv, tau, h)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO iz = 1, nzCMFD
    DO ipin = 1, nxy
      h = hzfm(iz)
      IF (mklCntl%lPolarXS) THEN

      ELSE
        sigt = PinXS(pinMap(ipin), planeMap(iz))%XStr(ig)
      ENDIF
      !--- Forward Coefficients
      a = srcCoeff(ipin, iz, ig, 1)
      b = srcCoeff(ipin, iz, ig, 2)
      c = srcCoeff(ipin, iz, ig, 3)
      DO ipol = 1, nPolarAngle
        cosv = Angle(ipol)%cosv
        C0(ipol, ig, iz, 1, ipin) = (2.0 * cosv ** 2 * a / sigt(ipol) ** 3 - cosv * b / sigt(ipol) ** 2             &
                                     + c / sigt(ipol))
        C1(ipol, ig, iz, 1, ipin) = (- 2.0 * cosv * a / sigt(ipol) ** 2 + b / sigt(ipol)) * h
        C1(ipol, ig, iz, 1, ipin) = C1(ipol, ig, iz, 1, ipin) + (a / sigt(ipol)) * h ** 2
      ENDDO
      !--- Backward Coefficients
      c = a * h ** 2 + b * h + c
      b = - (2.0 * a * h + b)
      DO ipol = 1, nPolarAngle
        cosv = Angle(ipol)%cosv
        C0(ipol, ig, iz, 2, ipin) = (2.0 * cosv ** 2 * a / sigt(ipol) ** 3 - cosv * b / sigt(ipol) ** 2             &
                                     + c / sigt(ipol))
        C1(ipol, ig, iz, 2, ipin) = (- 2.0 * cosv * a / sigt(ipol) ** 2 + b / sigt(ipol)) * h
        C1(ipol, ig, iz, 2, ipin) = C1(ipol, ig, iz, 2, ipin) + (a / sigt(ipol)) * h ** 2
      ENDDO
      DO ipol = 1, nPolarAngle
        tau = sigt(ipol) * h * Angle(ipol)%rcosv
        EX(ipol, ig, iz, ipin) = 1.0 - EXP(- tau)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetQuarticRTCoeff(PinXS)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)

TYPE(mklAngle_Type), POINTER :: Angle(:)
INTEGER :: ng, nxy, nzCMFD, nPolarAngle
INTEGER :: ig, ipol, ipin, iz
INTEGER, POINTER :: pinMap(:), planeMap(:)
REAL :: a, b, c, d, e, Coeff(5), sigt(mklGeom%nPolar1D), cosv, tau, h
REAL, POINTER :: hzfm(:)

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
nPolarAngle = mklGeom%nPolar1D
hzfm => mklGeom%hzfm
pinMap => mklGeom%pinMap
planeMap => mklGeom%planeMap
Angle => mklAxial%Angle

!$OMP PARALLEL PRIVATE(a, b, c, d, e, Coeff, sigt, cosv, tau, h)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO iz = 1, nzCMFD
    DO ipin = 1, nxy
      h = hzfm(iz)
      IF (mklCntl%lPolarXS) THEN

      ELSE
        sigt = PinXS(pinMap(ipin), planeMap(iz))%XStr(ig)
      ENDIF
      Coeff = srcCoeff(ipin, iz, ig, :)
      !--- Forward Coefficients
      a = Coeff(1); b = Coeff(2); c = Coeff(3); d = Coeff(4); e = Coeff(5)
      DO ipol = 1, nPolarAngle
        cosv = Angle(ipol)%cosv
        C0(ipol, ig, iz, 1, ipin) = (24.0 * cosv ** 4 * a / sigt(ipol) ** 5                                         &
                                     - 6.0 * cosv ** 3 * b / sigt(ipol) ** 4                                        &
                                     + 2.0 * cosv ** 2 * c / sigt(ipol) ** 3                                        &
                                     - cosv * d / sigt(ipol) ** 2                                                   &
                                     + e / sigt(ipol))
        C1(ipol, ig, iz, 1, ipin) = (- 24.0 * cosv ** 3 * a / sigt(ipol) ** 4                                       &
                                     + 6.0 * cosv ** 2 * b / sigt(ipol) ** 3                                        &
                                     - 2.0 * cosv * c / sigt(ipol) ** 2                                             &
                                     + d / sigt(ipol)) * h
        C1(ipol, ig, iz, 1, ipin) = C1(ipol, ig, iz, 1, ipin) +                                                     &
                                    (12.0 * cosv ** 2 * a / sigt(ipol) ** 3                                         &
                                     - 3.0 * cosv * b / sigt(ipol) ** 2                                             &
                                     + c / sigt(ipol)) * h ** 2
        C1(ipol, ig, iz, 1, ipin) = C1(ipol, ig, iz, 1, ipin) +                                                     &
                                    (- 4.0 * cosv * a / sigt(ipol) ** 2                                             &
                                     + b / sigt(ipol)) * h ** 3
        C1(ipol, ig, iz, 1, ipin) = C1(ipol, ig, iz, 1, ipin) +                                                     &
                                    (a / sigt(ipol)) * h ** 4
      ENDDO
      !--- Backward Coefficients
      a = Coeff(1)
      b = - Coeff(2) - 4.0 * Coeff(1) * h
      c = Coeff(3) + 3.0 * Coeff(2) * h + 6.0 * Coeff(1) * h ** 2
      d = - Coeff(4) - 2.0 * Coeff(3) * h - 3.0 * Coeff(2) * h ** 2 - 4.0 * Coeff(1) * h ** 3
      e = Coeff(5) + Coeff(4) * h + Coeff(3) * h ** 2 + Coeff(2) * h ** 3 + Coeff(1) * h ** 4
      DO ipol = 1, nPolarAngle
        cosv = Angle(ipol)%cosv
        C0(ipol, ig, iz, 2, ipin) = (24.0 * cosv ** 4 * a / sigt(ipol) ** 5                                         &
                                     - 6.0 * cosv ** 3 * b / sigt(ipol) ** 4                                        &
                                     + 2.0 * cosv ** 2 * c / sigt(ipol) ** 3                                        &
                                     - cosv * d / sigt(ipol) ** 2                                                   &
                                     + e / sigt(ipol))
        C1(ipol, ig, iz, 2, ipin) = (- 24.0 * cosv ** 3 * a / sigt(ipol) ** 4                                       &
                                     + 6.0 * cosv ** 2 * b / sigt(ipol) ** 3                                        &
                                     - 2.0 * cosv * c / sigt(ipol) ** 2                                             &
                                     + d / sigt(ipol)) * h
        C1(ipol, ig, iz, 2, ipin) = C1(ipol, ig, iz, 2, ipin) +                                                     &
                                    (12.0 * cosv ** 2 * a / sigt(ipol) ** 3                                         &
                                     - 3.0 * cosv * b / sigt(ipol) ** 2                                             &
                                     + c / sigt(ipol)) * h ** 2
        C1(ipol, ig, iz, 2, ipin) = C1(ipol, ig, iz, 2, ipin) +                                                     &
                                    (- 4.0 * cosv * a / sigt(ipol) ** 2                                             &
                                     + b / sigt(ipol)) * h ** 3
        C1(ipol, ig, iz, 2, ipin) = C1(ipol, ig, iz, 2, ipin) +                                                     &
                                    (a / sigt(ipol)) * h ** 4
      ENDDO
      DO ipol = 1, nPolarAngle
        tau = sigt(ipol) * h * Angle(ipol)%rcosv
        EX(ipol, ig, iz, ipin) = 1.0 - EXP(- tau)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetHyperbolicRTCoeff(PinXS)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)

TYPE(mklAngle_Type), POINTER :: Angle(:)
INTEGER :: ng, nxy, nzCMFD, nPolarAngle
INTEGER :: ig, ipol, ipin, iz
INTEGER, POINTER :: pinMap(:), planeMap(:)
REAL :: cosh1, sinh1
REAL :: a, b, c, d, e, Coeff(5), sigt(mklGeom%nPolar1D), cosv, tau, h, s, l
REAL, POINTER :: hzfm(:)

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
nPolarAngle = mklGeom%nPolar1D
hzfm => mklGeom%hzfm
pinMap => mklGeom%pinMap
planeMap => mklGeom%planeMap
Angle => mklAxial%Angle

cosh1 = dcosh(1.0D0)
sinh1 = dsinh(1.0D0)

!$OMP PARALLEL PRIVATE(a, b, c, d, e, Coeff, sigt, cosv, tau, h, s, l)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO iz = 1, nzCMFD
    DO ipin = 1, nxy
      h = 0.5 * hzfm(iz)
      IF (mklCntl%lPolarXS) THEN

      ELSE
        sigt = PinXS(pinMap(ipin), planeMap(iz))%XStr(ig)
      ENDIF
      Coeff = srcCoeff(ipin, iz, ig, :)
      !--- Forward Coefficients
      a = Coeff(1); b = Coeff(2); c = Coeff(3) / h ** 2; d = Coeff(4) / h; e = Coeff(5)
      DO ipol = 1, nPolarAngle
        cosv = Angle(ipol)%cosv; s = sigt(ipol); l = h / cosv

      ENDDO
      !--- Backward Coefficients
      d = - (Coeff(4) + 2.0 * Coeff(3)) / hzfm(iz)
      e = Coeff(5) + Coeff(4) + Coeff(3)
      DO ipol = 1, nPolarAngle
        cosv = Angle(ipol)%cosv; s = sigt(ipol); h = 0.5 * hzfm(iz)

      ENDDO
      DO ipol = 1, nPolarAngle
        tau = sigt(ipol) * hzfm(iz) * Angle(ipol)%rcosv
        EX(ipol, ig, iz, ipin) = 1.0 - EXP(- tau)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE CommBoundaryFluxConv(lConv)
USE PE_MOD,         ONLY : PE
IMPLICIT NONE

INTEGER :: n, ng, nxy, nzCMFD, nPolarAngle
INTEGER :: ipol, ig, ipin, iz, ierr
REAL :: err, errmax, tol = 1.0D-06
REAL, POINTER :: PhiAngIn(:, :, :), PhiAngOut(:, :, :)
LOGICAL, OPTIONAL :: lConv
LOGICAL :: lConvCheck

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
nPolarAngle = mklGeom%nPolar1D
lConvCheck = PRESENT(lConv)

n = nPolarAngle * ng * nxy
errmax = 0.0

ALLOCATE(PhiAngIn(nPolarAngle, ng, nxy))

PhiAngOut => mklAxial%PhiAngOut(:, :, :, top)
CALL GetNeighbor(n, PhiAngOut, PhiAngIn, top)
IF (mklGeom%lBottom) THEN
  IF (mklGeom%AxBC(bottom) .EQ. VoidCell) THEN
    PhiAngIn = 0.0
  ELSEIF (mklGeom%AxBC(bottom) .EQ. RefCell) THEN
    CALL dcopy(n, mklAxial%PhiAngOut(:, :, :, bottom), 1, PhiAngIn, 1)
  ENDIF
ENDIF

IF (lConvCheck) THEN
  !$OMP PARALLEL PRIVATE(err)
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3) REDUCTION(MAX : errmax)
  DO ipin = 1, nxy
    DO ig = 1, ng
      DO ipol = 1, nPolarAngle
        IF (PhiAngIn(ipol, ig, ipin) .EQ. 0.0) CYCLE
        err = abs(PhiAngIn(ipol, ig, ipin) - mklAxial%PhiAngIn(ipol, ig, ipin, bottom)) / PhiAngIn(ipol, ig, ipin)
        errmax = max(err, errmax)
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDIF

CALL dcopy(n, PhiAngIn, 1, mklAxial%PhiAngIn(:, :, :, bottom), 1)

PhiAngOut => mklAxial%PhiAngOut(:, :, :, bottom)
CALL GetNeighbor(n, PhiAngOut, PhiAngIn, bottom)
IF (mklGeom%lTop) THEN
  IF (mklGeom%AxBC(top) .EQ. VoidCell) THEN
    PhiAngIn = 0.0
  ELSEIF (mklGeom%AxBC(top) .EQ. RefCell) THEN
    CALL dcopy(n, mklAxial%PhiAngOut(:, :, :, top), 1, PhiAngIn, 1)
  ENDIF
ENDIF

IF (lConvCheck) THEN
  !$OMP PARALLEL PRIVATE(err)
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3) REDUCTION(MAX : errmax)
  DO ipin = 1, nxy
    DO ig = 1, ng
      DO ipol = 1, nPolarAngle
        IF (PhiAngIn(ipol, ig, ipin) .EQ. 0.0) CYCLE
        err = abs(PhiAngIn(ipol, ig, ipin) - mklAxial%PhiAngIn(ipol, ig, ipin, top)) / PhiAngIn(ipol, ig, ipin)
        errmax = max(err, errmax)
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDIF

CALL dcopy(n, PhiAngIn, 1, mklAxial%PhiAngIn(:, :, :, top), 1)

DEALLOCATE(PhiAngIn)

IF (lConvCheck) THEN
  err = errmax
  CALL MPI_ALLREDUCE(err, errmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, PE%MPI_CMFD_COMM, ierr)
  lConv = errmax .LE. tol
ENDIF

END SUBROUTINE

SUBROUTINE PolynomialRayTraceP0(ipin, lJout)

IMPLICIT NONE

TYPE(mklAngle_Type), POINTER :: Angle(:)
REAL, POINTER :: Jout(:, :, :, :), PhiAngIn(:, :, :, :), PhiAngOut(:, :, :, :)
REAL, POINTER :: myC0(:, :, :, :), myC1(:, :, :, :), myEX(:, :, :)
REAL :: phid, track_phi(mklGeom%nPolar1D, mklGeom%ng)
INTEGER :: ig, ipin, ipol, izf
INTEGER :: ng, nPolarAngle, nzCMFD
LOGICAL :: lJout

Jout => mklAxial%Jout(:, :, :, :, ipin)
PhiAngIn => mklAxial%PhiAngIn
PhiAngOut => mklAxial%PhiAngOut
myC0 => C0(:, :, :, :, ipin)
myC1 => C1(:, :, :, :, ipin)
myEX => EX(:, :, :, ipin)
Angle => mklAxial%Angle

ng = mklGeom%ng
nPolarAngle = mklGeom%nPolar1D
nzCMFD = mklGeom%nzCMFD

!--- Upward Sweep

track_phi = PhiAngIn(:, :, ipin, bottom)

DO izf = 1, nzCMFD
  IF (lJout) THEN
    DO ig = 1, ng
      DO ipol = 1, nPolarAngle
        Jout(in, ig, bottom, izf) = Jout(in, ig, bottom, izf) + Angle(ipol)%wtsurf * track_phi(ipol, ig)
        Jout(surf, ig, bottom, izf) = Jout(surf, ig, bottom, izf) + Angle(ipol)%wt * track_phi(ipol, ig)
      ENDDO
    ENDDO
  ENDIF
  DO ig = 1, ng
    DO ipol = 1, nPolarAngle
      phid = myC1(ipol, ig, izf, 1) + (myC0(ipol, ig, izf, 1) - track_phi(ipol, ig)) * myEX(ipol, ig, izf)
      track_phi(ipol, ig) = track_phi(ipol, ig) + phid
    ENDDO
  ENDDO
  IF (lJout) THEN
    DO ig = 1, ng
      DO ipol = 1, nPolarAngle
        Jout(out, ig, top, izf) = Jout(out, ig, top, izf) + Angle(ipol)%wtsurf * track_phi(ipol, ig)
        Jout(surf, ig, top, izf) = Jout(surf, ig, top, izf) + Angle(ipol)%wt * track_phi(ipol, ig)
      ENDDO
    ENDDO
  ENDIF
ENDDO

PhiAngOut(:, :, ipin, top) = track_phi

!--- Downward Sweep

track_phi = PhiAngIn(:, :, ipin, top)

DO izf = nzCMFD, 1, -1
  IF (lJout) THEN
    DO ig = 1, ng
      DO ipol = 1, nPolarAngle
        Jout(in, ig, top, izf) = Jout(in, ig, top, izf) + Angle(ipol)%wtsurf * track_phi(ipol, ig)
        Jout(surf, ig, top, izf) = Jout(surf, ig, top, izf) + Angle(ipol)%wt * track_phi(ipol, ig)
      ENDDO
    ENDDO
  ENDIF
  DO ig = 1, ng
    DO ipol = 1, nPolarAngle
      phid = myC1(ipol, ig, izf, 2) + (myC0(ipol, ig, izf, 2) - track_phi(ipol, ig)) * myEX(ipol, ig, izf)
      track_phi(ipol, ig) = track_phi(ipol, ig) + phid
    ENDDO
  ENDDO
  IF (lJout) THEN
    DO ig = 1, ng
      DO ipol = 1, nPolarAngle
        Jout(out, ig, bottom, izf) = Jout(out, ig, bottom, izf) + Angle(ipol)%wtsurf * track_phi(ipol, ig)
        Jout(surf, ig, bottom, izf) = Jout(surf, ig, bottom, izf) + Angle(ipol)%wt * track_phi(ipol, ig)
      ENDDO
    ENDDO
  ENDIF
ENDDO

PhiAngOut(:, :, ipin, bottom) = track_phi

END SUBROUTINE

SUBROUTINE PolynomialRayTracePN(ipin, lJout)

IMPLICIT NONE

TYPE(mklAngle_Type), POINTER :: Angle(:)
REAL, POINTER :: Jout(:, :, :, :), PhiAngIn(:, :, :, :), PhiAngOut(:, :, :, :)
REAL, POINTER :: tauInv(:, :, :)
REAL :: src, phid, track_phi(mklGeom%nPolar1D, mklGeom%ng)
INTEGER :: ig, ipin, ipol, izf
INTEGER :: ng, nPolarAngle, nzCMFD
LOGICAL :: lJout

Jout => mklAxial%Jout(:, :, :, :, ipin)
PhiAngIn => mklAxial%PhiAngIn
PhiAngOut => mklAxial%PhiAngOut
Angle => mklAxial%Angle

ng = mklGeom%ng
nPolarAngle = mklGeom%nPolar1D
nzCMFD = mklGeom%nzCMFD

!--- Upward Sweep

track_phi = PhiAngIn(:, :, ipin, bottom)

DO izf = 1, nzCMFD
  IF (lJout) THEN
    DO ig = 1, ng
      DO ipol = 1, nPolarAngle
        Jout(in, ig, bottom, izf) = Jout(in, ig, bottom, izf) + Angle(ipol)%wtsurf * track_phi(ipol, ig)
        Jout(surf, ig, bottom, izf) = Jout(surf, ig, bottom, izf) + Angle(ipol)%wt * track_phi(ipol, ig)
      ENDDO
    ENDDO
  ENDIF
  DO ig = 1, ng
    DO ipol = 1, nPolarAngle
      src = C0(ipol, ig, izf, 1, ipin) + SrcAng(ipol, ig, izf, 1, ipin)
      phid = C1(ipol, ig, izf, 1, ipin) + (src - track_phi(ipol, ig)) * EX(ipol, ig, izf, ipin)
      track_phi(ipol, ig) = track_phi(ipol, ig) + phid
      phia(ipol, ig, izf, 1, ipin) = - phid * tauInv(ipol, ig, izf)
    ENDDO
  ENDDO
  IF (lJout) THEN
    DO ig = 1, ng
      DO ipol = 1, nPolarAngle
        Jout(out, ig, top, izf) = Jout(out, ig, top, izf) + Angle(ipol)%wtsurf * track_phi(ipol, ig)
        Jout(surf, ig, top, izf) = Jout(surf, ig, top, izf) + Angle(ipol)%wt * track_phi(ipol, ig)
      ENDDO
    ENDDO
  ENDIF
ENDDO

PhiAngOut(:, :, ipin, top) = track_phi

!--- Downward Sweep

track_phi = PhiAngIn(:, :, ipin, top)

DO izf = nzCMFD, 1, -1
  IF (lJout) THEN
    DO ig = 1, ng
      DO ipol = 1, nPolarAngle
        Jout(in, ig, top, izf) = Jout(in, ig, top, izf) + Angle(ipol)%wtsurf * track_phi(ipol, ig)
        Jout(surf, ig, top, izf) = Jout(surf, ig, top, izf) + Angle(ipol)%wt * track_phi(ipol, ig)
      ENDDO
    ENDDO
  ENDIF
  DO ig = 1, ng
    DO ipol = 1, nPolarAngle
      src = C0(ipol, ig, izf, 1, ipin) + SrcAng(ipol, ig, izf, 1, ipin)
      phid = C1(ipol, ig, izf, 2, ipin) + (src - track_phi(ipol, ig)) * EX(ipol, ig, izf, ipin)
      track_phi(ipol, ig) = track_phi(ipol, ig) + phid
      phia(ipol, ig, izf, 2, ipin) = - phid * tauInv(ipol, ig, izf)
    ENDDO
  ENDDO
  IF (lJout) THEN
    DO ig = 1, ng
      DO ipol = 1, nPolarAngle
        Jout(out, ig, bottom, izf) = Jout(out, ig, bottom, izf) + Angle(ipol)%wtsurf * track_phi(ipol, ig)
        Jout(surf, ig, bottom, izf) = Jout(surf, ig, bottom, izf) + Angle(ipol)%wt * track_phi(ipol, ig)
      ENDDO
    ENDDO
  ENDIF
ENDDO

PhiAngOut(:, :, ipin, bottom) = track_phi

END SUBROUTINE

END MODULE

#endif    