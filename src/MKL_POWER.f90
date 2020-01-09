#include <defines.h>
!--- CNJ Edit : 3D CMFD Acceleration Modules with Intel MKL 
#ifdef __INTEL_MKL

!--- Chebyshev Acceleration Routines --------------------------------------------------------------

MODULE MKL_CHEBYSHEV

USE MKL_3D
IMPLICIT NONE

REAL :: sigma, xsi, err(2)
REAL, POINTER :: chebyCoeff(:, :)
INTEGER :: Order

CONTAINS

SUBROUTINE ChebyshevInit()

IMPLICIT NONE

INTEGER :: m, l

ALLOCATE(chebyCoeff(0 : mklCntl%chebyOrder, 0 : mklCntl%chebyOrder))
chebyCoeff = 0.0

!--- Zeroth Order
chebyCoeff(0, 0) = 1.0

!--- First Order
chebyCoeff(1, 1) = 1.0

!--- Second Order
chebyCoeff(0, 2) = - 1.0
chebyCoeff(2, 2) = 2.0

!--- Higher Order
DO m = 3, mklCntl%chebyOrder
  chebyCoeff(0, m) = - chebyCoeff(0, m - 2)
  DO l = 1, m
    chebyCoeff(l, m) = chebyCoeff(l - 1, m - 1) * 2.0 - chebyCoeff(l, m - 2)
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE ChebyshevReset()
USE IOUTIL,         ONLY : message
USE FILES,          ONLY : io8
USE PE_MOD,         ONLY : PE
IMPLICIT NONE

Order = 0

WRITE(mesg, '(a)') 'Restarting Chebyshev Acceleration Cycle...'
IF (PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)

WRITE(mesg, '(a, f6.4)') 'Estimated Dominance Ratio = ', sigma
IF (PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)

END SUBROUTINE

SUBROUTINE ChebyshevAcc(CMFD, lAcc)
USE PE_MOD,         ONLY : PE
IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
LOGICAL :: lAcc

INTEGER :: n, ng, nxy, nzCMFD
REAL :: alpha, beta, gamma, norm(2)
REAL, POINTER :: errVec(:, :)

ng = CMFD%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
n = ng * nxy * nzCMFD

IF (lAcc) THEN
  
  !--- Determine Extrapolation Parameters
  
  IF (Order .EQ. 0) THEN
    alpha = 1.0
    beta = 0.0
  ELSEIF (Order .EQ. 1) THEN
    alpha = 2.0 / (2.0 - sigma)
    beta = 0.0
  ELSE
    gamma = 2.0 / sigma - 1.0
    alpha = 4.0 / sigma * Chebyshev(gamma, Order - 1) / Chebyshev(gamma, Order)
    beta = Chebyshev(gamma, Order - 2) / Chebyshev(gamma, Order)
  ENDIF

  !--- Error Monitering
  
  ALLOCATE(errVec(n, 1))
  
  CALL vdsub(n, CMFD%phis(:, :, :), CMFD%phisd(:, :, :, 1), errVec(:, 1))
  
  norm(1) = normMPI(errVec(:, 1), n, PE%MPI_CMFD_COMM)
  
  err(2) = norm(1)
  IF (Order .EQ. 1) err(1) = norm(1)
  
  DEALLOCATE(errVec)
  
  !--- Extrapolation
  
  CALL daxpby(n, 1.0 - alpha + beta, CMFD%phisd(:, :, :, 1), 1, alpha, CMFD%phis, 1)
  CALL daxpy(n, -beta, CMFD%phisd(:, :, :, 2), 1, CMFD%phis, 1)
  
  !--- Increase Polynomial Order
  
  Order = Order + 1
  
  !--- Update Dominance Ratio and Check Cycle Status
  
  IF (Order .GT. 2) CALL ChebyshevCheckStatus()
  
ELSE

  !--- Estimate Dominance Ratio

  ALLOCATE(errVec(n, 2))

  CALL vdsub(n, CMFD%phis(:, :, :), CMFD%phisd(:, :, :, 1), errVec(:, 1))
  CALL vdsub(n, CMFD%phisd(:, :, :, 1), CMFD%phisd(:, :, :, 2), errVec(:, 2))

  norm(1) = normMPI(errVec(:, 1), n, PE%MPI_CMFD_COMM)
  norm(2) = normMPI(errVec(:, 2), n, PE%MPI_CMFD_COMM)
  
  sigma = norm(1) / norm(2)
  
  DEALLOCATE(errVec)
  
ENDIF

CALL dcopy(n, CMFD%phisd(:, :, :, 1), 1, CMFD%phisd(:, :, :, 2), 1)
CALL dcopy(n, CMFD%phis(:, :, :), 1, CMFD%phisd(:, :, :, 1), 1)
  
END SUBROUTINE

SUBROUTINE ChebyshevCheckStatus()

IMPLICIT NONE

REAL :: gamma

IF (Order .GT. mklCntl%chebyOrder) THEN
  CALL ChebyshevReset(); RETURN
ENDIF

IF (mklCntl%lDcpl) RETURN

! gamma = 2.0 / sigma - 1.0
! xsi = err(2) / err(1) * Chebyshev(gamma, Order - 1)
!  
! IF (xsi .LE. 1.0) gamma = dcos(dacos(xsi) / (Order - 1))
! IF (xsi .GT. 1.0) gamma = dcosh(dacosh(xsi) / (Order - 1))
! 
! IF (xsi .GT. 1.0) THEN
!   sigma = sigma * (1.0 + gamma) / 2.0  
!   CALL ChebyshevReset()
! ENDIF

END SUBROUTINE

FUNCTION Chebyshev(x, m) RESULT(val)

IMPLICIT NONE

REAL :: x, val
INTEGER :: m, l

val = 0.0
DO l = 0, m
  val = val + chebyCoeff(l, m) * x ** l
ENDDO

END FUNCTION

END MODULE
    
!--- Power Method Routines ------------------------------------------------------------------------
    
MODULE MKL_POWER

USE MKL_3D
IMPLICIT NONE

CONTAINS

SUBROUTINE SetCsrBiCGSystem(CMFD, PinXS, l3dim, lPrecond, seigv)
USE PARAM
USE geom,           ONLY : ncbd
USE TYPEDEF,        ONLY : PinXS_Type
USE MKL_BILU,       ONLY : MKL_PrepareILU,     MKL_PrepareSPAI
IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
LOGICAL :: l3dim, lPrecond
REAL :: seigv

TYPE(superPin_Type), POINTER :: Pin(:)
INTEGER, POINTER :: pinMap(:), pinMapRev(:), planeMap(:)
INTEGER :: DOWN, UP, SELF
INTEGER :: ir, ic, iz, izf, ig, igf, igt, ibd, isurf, ipin, ipin_map, ineighpin
INTEGER :: dz, gb, ge
INTEGER :: ng, nxy, nzCMFD, nbd
REAL, POINTER :: PinVolFm(:, :), hzfm(:)
REAL :: diagVal(ncbd+3), Dtil, Dhat, val

Pin => mklGeom%superPin
ng = CMFD%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
hzfm => mklGeom%hzfm
pinMap => mklGeom%pinMap
pinMapRev => mklGeom%pinMapRev
planeMap => CMFD%planeMap
PinVolFm => mklGeom%PinVolFm

DOWN = ncbd+1
UP   = ncbd+2
SELF = ncbd+3

!--- Set Group Major Diffusion Operator

!$OMP PARALLEL PRIVATE(diagVal, iz, ipin_map, Dtil, Dhat, isurf, ineighpin, dz, ir, ic, val)
!$OMP DO SCHEDULE(DYNAMIC)
DO ig = 1, ng
  CALL createCsr(CMFD%M(ig), (ncbd+3) * nxy * nzCMFD, nxy * nzCMFD, nxy * nzCMFD)
  DO izf = 1, nzCMFD
    iz = planeMap(izf)
    DO ipin = 1, nxy
      ir = ipin + (izf - 1) * nxy
      ipin_map = pinMap(ipin)
      diagVal = 0.0
      DO ibd = 1, pin(ipin)%nNgh
        Dtil = PinXS(ipin_map, iz)%Dtil(ibd, ig)
        Dhat = PinXS(ipin_map, iz)%Dhat(ibd, ig)
        diagVal(ibd) = - (Dtil + Dhat) * hzfm(izf)
        diagVal(SELF) = diagVal(SELF) + (Dtil - Dhat) * hzfm(izf)
      ENDDO
      IF (l3dim) THEN
        DO ibd = DOWN, UP
          Dtil = CMFD%AxDtil(ibd - ncbd, ipin, izf, ig)
          Dhat = CMFD%AxDhat(ibd - ncbd, ipin, izf, ig)
          diagVal(ibd) = - (Dtil + Dhat) * PinVolFm(ipin, izf) / hzfm(izf)
          diagVal(SELF) = diagVal(SELF) + (Dtil - Dhat) * PinVolFm(ipin, izf) / hzfm(izf)
        ENDDO
      ENDIF
      diagVal(SELF) = diagVal(SELF) + PinVolFm(ipin, izf) * PinXS(ipin_map, iz)%XSr(ig)
      IF (mklCntl%lShift) diagVal(SELF) = diagVal(SELF) - CMFD%F(ir, ig) * CMFD%Chi(ir, ig) * seigv
      DO ibd = 1, ncbd+3
        isurf = ibd
        
        IF (isurf .EQ. SELF) THEN
          ineighpin = ipin
          dz = 0
        ELSE IF (isurf .EQ. UP) THEN
          ineighpin = ipin
          IF (mklCntl%DcplLv .EQ. 2) diagVal(isurf) = 0.0
          dz = 1
        ELSE IF (isurf .EQ. DOWN) THEN
          ineighpin = ipin
          IF (mklCntl%DcplLv .EQ. 2) diagVal(isurf) = 0.0
          dz = -1
        ELSE
          ineighpin = Pin(ipin_map)%NeighIdx(isurf)
          ineighpin = pinMapRev(ineighpin)
          IF (ineighpin .LE. 0) diagVal(isurf) = 0.0
          IF (ineighpin .EQ. ipin) THEN
            Dtil = PinXS(ipin_map, iz)%Dtil(isurf, ig)
            Dhat = PinXS(ipin_map, iz)%Dhat(isurf, ig)
            diagVal(SELF) = diagVal(SELF) - (Dtil + Dhat) * hzfm(izf)
            diagVal(isurf) = 0.0
          ENDIF
          dz = 0
        END IF
        ic = ir + (ineighpin - ipin) + dz * nxy
        val = diagVal(isurf)
        CALL pushCsr(CMFD%M(ig), val, ir, ic)
      ENDDO
    ENDDO
  ENDDO
  CALL finalizeSortCsr(CMFD%M(ig), FALSE)
ENDDO
!$OMP END DO
!$OMP END PARALLEL

!--- Factorization

IF (mklCntl%lPardiso) THEN
  DO ig = 1, ng
    CALL pardisoCreate(CMFD%M(ig))
  ENDDO
ELSE
  IF (lPrecond) THEN
    IF (mklCntl%lSPAI) THEN
      !$OMP PARALLEL DO SCHEDULE(DYNAMIC)
      DO ig = 1, ng
        CALL MKL_PrepareSPAI(CMFD%M(ig), CMFD%SPAI(ig))
      ENDDO
      !$OMP END PARALLEL DO
    ELSE
      !$OMP PARALLEL DO SCHEDULE(DYNAMIC)
      DO ig = 1, ng
        CALL MKL_PrepareILU(CMFD%M(ig), CMFD%ILU(ig))
      ENDDO
      !$OMP END PARALLEL DO
    ENDIF
  ENDIF
ENDIF

!--- Set Axial Off-diagonals for MPI

IF (l3dim) THEN
  IF (.NOT. mklGeom%lBottom) THEN
    !$OMP PARALLEL PRIVATE(Dtil, Dhat, val)
    !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
    DO ig = 1, ng
      DO ipin = 1, nxy
        Dtil = CMFD%AxDtil(bottom, ipin, 1, ig)
        Dhat = CMFD%AxDhat(bottom, ipin, 1, ig)
        val = - (Dtil + Dhat) * PinVolFm(ipin, 1) / hzfm(1)
        CMFD%AxOffDiag(ipin, bottom, ig) = val
      ENDDO
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL
  ENDIF
  IF (.NOT. mklGeom%lTop) THEN
    !$OMP PARALLEL PRIVATE(Dtil, Dhat, val)
    !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
    DO ig = 1, ng
      DO ipin = 1, nxy
        Dtil = CMFD%AxDtil(top, ipin, nzCMFD, ig)
        Dhat = CMFD%AxDhat(top, ipin, nzCMFD, ig)
        val = - (Dtil + Dhat) * PinVolFm(ipin, nzCMFD) / hzfm(nzCMFD)
        CMFD%AxOffDiag(ipin, top, ig) = val
      ENDDO
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL
  ENDIF
ENDIF

!--- Set Axial Off-diagonals for Level 2 Decoupling

IF (l3dim .AND. mklCntl%DcplLv .EQ. 2) THEN
  !$OMP PARALLEL PRIVATE(Dtil, Dhat)
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
  DO ig = 1, ng
    DO izf = 1, nzCMFD
      DO ipin = 1, nxy
        IF (izf .NE. 1) THEN
          Dtil = CMFD%AxDtil(bottom, ipin, izf, ig)
          Dhat = CMFD%AxDhat(bottom, ipin, izf, ig)
          CMFD%dcplAxOffDiag(ipin, izf, bottom, ig) = - (Dtil + Dhat) * PinVolFm(ipin, izf) / hzfm(izf)
        ENDIF
        IF (izf .NE. nzCMFD) THEN
          Dtil = CMFD%AxDtil(top, ipin, izf, ig)
          Dhat = CMFD%AxDhat(top, ipin, izf, ig)
          CMFD%dcplAxOffDiag(ipin, izf, top, ig) = - (Dtil + Dhat) * PinVolFm(ipin, izf) / hzfm(izf)
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDIF

END SUBROUTINE

SUBROUTINE SetCsrJacobiSystem(CMFD, PinXS, l3dim, seigv)
USE PARAM
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
LOGICAL :: l3dim
REAL :: seigv

TYPE(superPin_Type), POINTER :: Pin(:)
INTEGER, POINTER :: pinMap(:), pinMapRev(:), planeMap(:)
INTEGER, PARAMETER :: DOWN = 5, UP = 6, SELF = 7
INTEGER :: surf(7) = (/ UP, NORTH, WEST, SELF, EAST, SOUTH, DOWN /)
INTEGER :: ir, ic, iz, izf, ig, igf, igt, ibd, isurf, ipin, ipin_map, ineighpin
INTEGER :: dz, gb, ge
INTEGER :: ng, nxy, nzCMFD, nbd
REAL, POINTER :: PinVolFm(:, :), hzfm(:)
REAL :: diagVal(7), Dtil, Dhat, val

Pin => mklGeom%superPin
ng = CMFD%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
hzfm => mklGeom%hzfm
planeMap => mklGeom%planeMap
pinMap => mklGeom%pinMap
pinMapRev => mklGeom%pinMapRev
PinVolFm => mklGeom%PinVolFm

!--- Set Group Major Diffusion Operator

!$OMP PARALLEL PRIVATE(diagVal, iz, ipin_map, Dtil, Dhat, isurf, ineighpin, dz, ir, ic, val)
!$OMP DO SCHEDULE(DYNAMIC)
DO ig = 1, ng
  CALL createCsr(CMFD%M(ig), 7 * nxy * nzCMFD, nxy * nzCMFD, nxy * nzCMFD)
  CALL createCsr(CMFD%Jacobi(ig)%M, 6 * nxy * nzCMFD, nxy * nzCMFD, nxy * nzCMFD)
  DO izf = 1, nzCMFD
    iz = planeMap(izf)
    DO ipin = 1, nxy
      ir = ipin + (izf - 1) * nxy
      ipin_map = pinMap(ipin)
      diagVal = 0.0
      DO ibd = 1, 4
        Dtil = PinXS(ipin_map, iz)%Dtil(ibd, ig)
        Dhat = PinXS(ipin_map, iz)%Dhat(ibd, ig)
        diagVal(ibd) = - (Dtil + Dhat) * hzfm(izf)
        diagVal(SELF) = diagVal(SELF) + (Dtil - Dhat) * hzfm(izf)
      ENDDO
      IF (l3dim) THEN
        DO ibd = 5, 6
          Dtil = CMFD%AxDtil(ibd - 4, ipin, izf, ig)
          Dhat = CMFD%AxDhat(ibd - 4, ipin, izf, ig)
          diagVal(ibd) = - (Dtil + Dhat) * PinVolFm(ipin, izf) / hzfm(izf)
          diagVal(SELF) = diagVal(SELF) + (Dtil - Dhat) * PinVolFm(ipin, izf) / hzfm(izf)
        ENDDO
      ENDIF
      diagVal(SELF) = diagVal(SELF) + PinVolFm(ipin, izf) * PinXS(ipin_map, iz)%XSr(ig)
      IF (mklCntl%lShift) diagVal(SELF) = diagVal(SELF) - CMFD%F(ir, ig) * CMFD%Chi(ir, ig) * seigv
      DO ibd = 1, 7
        isurf = surf(ibd)
        SELECT CASE (isurf)
        CASE (SELF)
          ineighpin = ipin
          dz = 0
        CASE (UP)
          ineighpin = ipin
          dz = 1
        CASE (DOWN)
          ineighpin = ipin
          dz = -1
        CASE (NORTH, WEST, EAST, SOUTH)
          ineighpin = Pin(ipin_map)%NeighIdx(isurf)
          ineighpin = pinMapRev(ineighpin)
          IF (ineighpin .LE. 0) diagVal(isurf) = 0.0
          dz = 0
        END SELECT
        ic = ir + (ineighpin - ipin) + dz * nxy
        val = diagVal(isurf)
        CALL pushCsr(CMFD%M(ig), val, ir, ic)
        IF (isurf .EQ. SELF) THEN
          CMFD%Jacobi(ig)%invDiag(ir) = 1.0 / val; CYCLE
        ENDIF
        CALL pushCsr(CMFD%Jacobi(ig)%M, val, ir, ic)
      ENDDO
    ENDDO
  ENDDO
  CALL finalizeSortCsr(CMFD%M(ig), FALSE)
  CALL finalizeSortCsr(CMFD%Jacobi(ig)%M, FALSE)
ENDDO
!$OMP END DO
!$OMP END PARALLEL

!--- Set Axial Off-diagonals for MPI

IF (l3dim) THEN
  IF (.NOT. mklGeom%lBottom) THEN
    !$OMP PARALLEL PRIVATE(Dtil, Dhat, val)
    !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
    DO ig = 1, ng
      DO ipin = 1, nxy
        Dtil = CMFD%AxDtil(bottom, ipin, 1, ig)
        Dhat = CMFD%AxDhat(bottom, ipin, 1, ig)
        val = - (Dtil + Dhat) * PinVolFm(ipin, 1) / hzfm(1)
        CMFD%Jacobi(ig)%offDiag(ipin, bottom) = val
      ENDDO
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL
  ENDIF
  IF (.NOT. mklGeom%lTop) THEN
    !$OMP PARALLEL PRIVATE(Dtil, Dhat, val)
    !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
    DO ig = 1, ng
      DO ipin = 1, nxy
        Dtil = CMFD%AxDtil(top, ipin, nzCMFD, ig)
        Dhat = CMFD%AxDhat(top, ipin, nzCMFD, ig)
        val = - (Dtil + Dhat) * PinVolFm(ipin, nzCMFD) / hzfm(nzCMFD)
        CMFD%Jacobi(ig)%offDiag(ipin, top) = val
      ENDDO
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL
  ENDIF
ENDIF

END SUBROUTINE

SUBROUTINE CMFDSrcUpdt(CMFD, ig, eigv)
USE TIMER,          ONLY : nTracer_dclock,      TimeChk
IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
INTEGER :: iz, ig, igf, gb, ge
INTEGER :: n, ng, nxy, nzCMFD
REAL :: Tbeg, Tend
REAL :: eigv
REAL, POINTER :: scatSrc(:), fisSrc(:)
REAL, POINTER :: bottomSrc(:), topSrc(:), dcplSrc(:, :)

ng = CMFD%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
n = nxy * nzCMFD

ALLOCATE(scatSrc(n), fisSrc(n))
IF (mklCntl%lDcpl) THEN
  ALLOCATE(dcplSrc(nxy, nzCMFD)); dcplSrc = 0.0
  ALLOCATE(bottomSrc(nxy))
  ALLOCATE(topSrc(nxy))
ENDIF

CMFD%src(:, ig) = 0.0

gb = CMFD%InScatRange(1, ig)
ge = CMFD%InScatRange(2, ig)

Tbeg = nTracer_dclock(FALSE, FALSE)

DO igf = gb, ge
  IF (igf .EQ. ig) CYCLE
  CALL vdmul(n, CMFD%S(:, igf, ig), CMFD%phis(:, :, igf), scatSrc)
  CALL vdadd(n, scatSrc, CMFD%src(:, ig), CMFD%src(:, ig))
ENDDO

Tend = nTracer_dclock(FALSE, FALSE)
TimeChk%AxBTime = TimeChk%AxBTime + (Tend - Tbeg)

CALL vdmul(n, CMFD%Chi(:, ig), CMFD%psi, fisSrc)
CALL dscal(n, 1.0D0 / eigv, fisSrc, 1)

CALL vdadd(n, fisSrc, CMFD%src(:, ig), CMFD%src(:, ig))

IF (mklCntl%lDcpl) THEN
  IF (mklCntl%DcplLv .EQ. 2) THEN
    DO iz = 1, nzCMFD
      bottomSrc = 0.0; topSrc = 0.0
      IF (iz .NE. 1) CALL vdmul(nxy, CMFD%dcplAxOffDiag(:, iz, bottom, ig), CMFD%phis(:, iz - 1, ig), bottomSrc)
      IF (iz .NE. nzCMFD) CALL vdmul(nxy, CMFD%dcplAxOffDiag(:, iz, top, ig), CMFD%phis(:, iz + 1, ig), topSrc)
      CALL vdadd(nxy, bottomSrc, topSrc, dcplSrc(:, iz))
    ENDDO
  ENDIF
  IF (.NOT. mklGeom%lBottom) THEN
    CALL vdmul(nxy, CMFD%AxOffDiag(:, bottom, ig), CMFD%neighphis(:, ig, bottom), bottomSrc)
    CALL vdadd(nxy, bottomSrc, dcplSrc(:, 1), dcplSrc(:, 1))
  ENDIF
  IF (.NOT. mklGeom%lTop) THEN
    CALL vdmul(nxy, CMFD%AxOffDiag(:, top, ig), CMFD%neighphis(:, ig, top), topSrc)
    CALL vdadd(nxy, topSrc, dcplSrc(:, nzCMFD), dcplSrc(:, nzCMFD))
  ENDIF
  CALL vdsub(nxy * nzCMFD, CMFD%src(:, ig), dcplSrc, CMFD%src(:, ig))
ENDIF

DEALLOCATE(scatSrc, fisSrc)
IF (mklCntl%lDcpl) THEN
  DEALLOCATE(dcplSrc)
  DEALLOCATE(bottomSrc)
  DEALLOCATE(topSrc)
ENDIF

END SUBROUTINE

SUBROUTINE CMFDPsiUpdt(CMFD)

IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
INTEGER :: ig
INTEGER :: n, ng
REAL, POINTER :: fisSrc(:)

ng = CMFD%ng
n = mklGeom%nxy * mklGeom%nzCMFD

ALLOCATE(fisSrc(n))
fisSrc = 0.0

CALL dcopy(n, CMFD%psi, 1, CMFD%psid, 1)
CMFD%psi = 0.0

DO ig = 1, ng
  CALL vdmul(n, CMFD%F(:, ig), CMFD%phis(:, :, ig), fisSrc)
  CALL vdadd(n, fisSrc, CMFD%psi, CMFD%psi)
ENDDO

DEALLOCATE(fisSrc)

END SUBROUTINE

SUBROUTINE CMFDEigUpdt(CMFD, eigv)
USE PE_MOD,         ONLY : PE
USE IEEE_ARITHMETIC
IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
REAL :: eigv
REAL :: psipsi(2)
INTEGER :: n, ierr

n = mklGeom%nxy * mklGeom%nzCMFD

psipsi(1) = dotMPI(CMFD%psi, CMFD%psi, n, PE%MPI_CMFD_COMM)
psipsi(2) = dotMPI(CMFD%psi, CMFD%psid, n, PE%MPI_CMFD_COMM)

eigv = eigv * psipsi(1) / psipsi(2)

END SUBROUTINE

FUNCTION CMFDResidual(CMFD, eigv) RESULT(res)
USE PE_MOD,         ONLY : PE
IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
REAL :: eigv
REAL :: res, tsrc, buf_res, buf_tsrc
REAL, POINTER :: x(:, :), r(:)
REAL, POINTER :: csrVal(:)
INTEGER, POINTER :: csrRowPtr(:), csrColIdx(:)
INTEGER :: ig, ierr
INTEGER :: ng, nxy, nzCMFD, nr, nc
INTEGER :: bottomRange(2), topRange(2)

ng = CMFD%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
bottomRange = CMFD%bottomRange
topRange = CMFD%topRange

res = 0.0
tsrc = 0.0

DO ig = 1, ng
  
  CALL CMFDSrcUpdt(CMFD, ig, eigv)
  
  x => CMFD%phis(:, :, ig)
  csrVal => CMFD%M(ig)%csrVal
  csrRowPtr => CMFD%M(ig)%csrRowPtr
  csrColIdx => CMFD%M(ig)%csrColIdx
  nr = CMFD%M(ig)%nr
  nc = CMFD%M(ig)%nc
  
  ALLOCATE(r(nc))
  
  CALL MatMulMPIFast(CMFD%M(ig), CMFD%AxOffDiag(:, :, ig), x, r, bottomRange, topRange)
  CALL vdsub(nc, CMFD%src(:, ig), r, r)
  res = res + ddot(nc, r, 1, r, 1)
  tsrc = tsrc + ddot(nc, CMFD%src(:, ig), 1, CMFD%src(:, ig), 1)
  
  DEALLOCATE(r)
  
ENDDO

CALL MPI_ALLREDUCE(res, buf_res, 1, MPI_DOUBLE_PRECISION, MPI_SUM, PE%MPI_CMFD_COMM, ierr)
CALL MPI_ALLREDUCE(tsrc, buf_tsrc, 1, MPI_DOUBLE_PRECISION, MPI_SUM, PE%MPI_CMFD_COMM, ierr)

res = sqrt(buf_res / buf_tsrc)

END FUNCTION

FUNCTION CMFDPsiErr(CMFD) RESULT(err)
USE PE_MOD,         ONLY : PE
IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD

INTEGER :: ierr
INTEGER :: n, nxy, nzCMFD
REAL, POINTER :: e(:)
REAL :: err

nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
n = nxy * nzCMFD

ALLOCATE(e(n))

CALL vdsub(n, CMFD%psi, CMFD%psid, e)
err = normMPI(e, n, PE%MPI_CMFD_COMM) / normMPI(CMFD%psi, n, PE%MPI_CMFD_COMM)

DEALLOCATE(e)

END FUNCTION

! ===================================================================== !
!    Preconditioned Bi-Conjugate Gradient Stablized(BiCGSTAB) Solver    !
!  http://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method  !
! ===================================================================== !

SUBROUTINE BiCGSTAB(CMFD, ig)
USE PE_MOD,         ONLY : PE
USE TIMER,          ONLY : nTracer_dclock,      TimeChk
USE MKL_BILU,       ONLY : MKL_SolveILU,        MKL_ApplySPAI
IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
INTEGER :: ig

TYPE(CSR_DOUBLE), POINTER :: M, ILU, SPAI
INTEGER :: i, iter, itermin, itermax
INTEGER :: nr, nc, nnz
INTEGER :: nxy, nzCMFD
REAL :: Tbeg, Tend
REAL :: err, tol
REAL :: rho0, rho1, w0, w1, norm0, norm1, dot0, dot1, alpha, beta
REAL, POINTER :: x(:, :), b(:)
REAL, POINTER :: h(:), s(:), t(:), y(:), z(:)
REAL, POINTER :: r0(:), r1(:), rhat(:), v0(:), v1(:), p0(:), p1(:)
REAL, POINTER :: csrVal(:)
INTEGER, POINTER :: csrRowPtr(:), csrColIdx(:)
INTEGER :: bottomRange(2), topRange(2)
LOGICAL :: lConv

Tbeg = nTracer_dclock(FALSE, FALSE)

M => CMFD%M(ig)
ILU => CMFD%ILU(ig)
SPAI => CMFD%SPAI(ig)
x => CMFD%phis(:, :, ig)
b => CMFD%src(:, ig)
csrVal => M%csrVal
csrRowPtr => M%csrRowPtr
csrColIdx => M%csrColIdx

nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
bottomRange = CMFD%bottomRange
topRange = CMFD%topRange

itermin = mklCntl%minInner
itermax = mklCntl%maxInner
tol = mklCntl%innerConv

nr = M%nr
nc = M%nc
nnz = M%nnz

ALLOCATE(h(nc), s(nc), t(nc), y(nc), z(nc))
ALLOCATE(r0(nc), r1(nc), rhat(nc))
ALLOCATE(v0(nc), v1(nc))
ALLOCATE(p0(nc), p1(nc))

!--- 1 ----------------------
CALL MatMulMPIFast(M, CMFD%AxOffDiag(:, :, ig), x, r1, bottomRange, topRange)
CALL vdsub(nc, b, r1, r1)
!--- 2 ----------------------
CALL dcopy(nc, r1, 1, rhat, 1)
!--- 3 ----------------------
IF (mklCntl%lDcpl) THEN
  norm0 = dnrm2(nc, r1, 1)
ELSE
  norm0 = normMPI(r1, nc, PE%MPI_CMFD_COMM)
ENDIF
rho1 = 1.0; w1 = 1.0; alpha = 1.0
!--- 4 ----------------------
v1 = 0.0; p1 = 0.0
!--- 5 ----------------------
lConv = FALSE; iter = 0
DO WHILE (.NOT. lConv)
  !--- 5.0 ------------------
  iter = iter + 1
  rho0 = rho1; w0 = w1
  CALL dcopy(nc, r1, 1, r0, 1)
  CALL dcopy(nc, p1, 1, p0, 1)
  CALL dcopy(nc, v1, 1, v0, 1)
  !--- 5.1 ------------------
  IF (mklCntl%lDcpl) THEN
    rho1 = ddot(nc, rhat, 1, r0, 1)
  ELSE
    rho1 = dotMPI(rhat, r0, nc, PE%MPI_CMFD_COMM)
  ENDIF
  !--- 5.2 ------------------
  beta = (rho1 / rho0) * (alpha / w0)
  !--- 5.3 ------------------
  CALL dcopy(nc, r0, 1, p1, 1)
  CALL daxpy(nc, -w0, v0, 1, p0, 1)
  CALL daxpy(nc, beta, p0, 1, p1, 1)
  !--- 5.4 ------------------
  IF (mklCntl%lSPAI) THEN
    CALL MKL_ApplySPAI(SPAI, p1, y)
  ELSE
    CALL MKL_SolveILU(ILU, y, p1)
  ENDIF
  !--- 5.5 ------------------
  CALL MatMulMPIFast(M, CMFD%AxOffDiag(:, :, ig), y, v1, bottomRange, topRange)
  !--- 5.6 ------------------
  IF (mklCntl%lDcpl) THEN
    dot1 = ddot(nc, rhat, 1, v1, 1)
  ELSE
    dot1 = dotMPI(rhat, v1, nc, PE%MPI_CMFD_COMM)
  ENDIF
  alpha = rho1 / dot1
  !--- 5.7 ------------------
  CALL dcopy(nc, x, 1, h, 1)
  CALL daxpy(nc, alpha, y, 1, h, 1)
  !--- 5.9 ------------------
  CALL dcopy(nc, r0, 1, s, 1)
  CALL daxpy(nc, -alpha, v1, 1, s, 1)
  !--- 5.10 -----------------
  IF (mklCntl%lSPAI) THEN
    CALL MKL_ApplySPAI(SPAI, s, z)
  ELSE
    CALL MKL_SolveILU(ILU, z, s)
  ENDIF
  !--- 5.11 -----------------
  CALL MatMulMPIFast(M, CMFD%AxOffDiag(:, :, ig), z, t, bottomRange, topRange)
  !--- 5.12 -----------------
  IF (mklCntl%lDcpl) THEN
    dot0 = ddot(nc, t, 1, s, 1)
    dot1 = ddot(nc, t, 1, t, 1)
  ELSE
    dot0 = dotMPI(t, s, nc, PE%MPI_CMFD_COMM)
    dot1 = dotMPI(t, t, nc, PE%MPI_CMFD_COMM)
  ENDIF
  w1 = dot0 / dot1
  !--- 5.13 -----------------
  CALL dcopy(nc, h, 1, x, 1)
  CALL daxpy(nc, w1, z, 1, x, 1)
  !--- 5.15 -----------------
  CALL dcopy(nc, s, 1, r1, 1)
  CALL daxpy(nc, -w1, t, 1, r1, 1)
  !--- Convergence ----------
  IF (mklCntl%lDcpl) THEN
    norm1 = dnrm2(nc, r1, 1)
  ELSE
    norm1 = normMPI(r1, nc, PE%MPI_CMFD_COMM)
  ENDIF
  err = norm1 / norm0; lConv = (err .LE. tol) .AND. (iter .GE. itermin)
  lConv = lConv .OR. (iter .GE. itermax)
ENDDO

DEALLOCATE(h, s, t, y, z)
DEALLOCATE(r0, r1, rhat)
DEALLOCATE(v0, v1)
DEALLOCATE(p0, p1)

Tend = nTracer_dclock(FALSE, FALSE)
TimeChk%AxBTime = TimeChk%AxBTime + (Tend - Tbeg)

END SUBROUTINE

! ============================================================ !
!            Intel MKL Sparse Direct Solver PARDISO            !
!  http://software.intel.com/en-us/articles/intel-mkl-pardiso  !
! ============================================================ !

SUBROUTINE pardisoCreate(M)

IMPLICIT NONE

TYPE(CSR_DOUBLE) :: M
INTEGER :: ierr
INTEGER :: nr, nc, nnz
INTEGER, PARAMETER :: mtype = 1, phase = 12
INTEGER, POINTER :: csrRowPtr(:), csrColIdx(:)
REAL, POINTER :: csrVal(:)
REAL, POINTER :: x(:), b(:)   !--- Dummy

csrVal => M%csrVal
csrRowPtr => M%csrRowPtr
csrColIdx => M%csrColIdx
nr = M%nr
nc = M%nc
nnz = M%nnz

CALL pardisoinit(M%pardisoPtr, mtype, M%pardisoParam)

M%pardisoParam(5) = 2
M%pardisoParam(6) = 0
M%pardisoParam(28) = 0
M%pardisoParam(31) = 0
M%pardisoParam(36) = 0

CALL pardiso(M%pardisoPtr, 1, 1, mtype, phase, nr, csrVal, csrRowPtr, csrColIdx, M%pardisoPermute, 1,               &
             M%pardisoParam, 0, b, x, ierr)

END SUBROUTINE

SUBROUTINE pardisoSolve(M, x, b)
USE PE_MOD,         ONLY : PE
USE TIMER,          ONLY : nTracer_dclock,      TimeChk
IMPLICIT NONE

TYPE(CSR_DOUBLE) :: M
INTEGER :: ierr
INTEGER :: nr, nc, nnz
INTEGER, PARAMETER :: mtype = 1, phase = 33
INTEGER, POINTER :: csrRowPtr(:), csrColIdx(:)
REAL, POINTER :: csrVal(:)
REAL :: x(*), b(:)
REAL :: Tbeg, Tend

Tbeg = nTracer_dclock(FALSE, FALSE)

csrVal => M%csrVal
csrRowPtr => M%csrRowPtr
csrColIdx => M%csrColIdx
nr = M%nr
nc = M%nc
nnz = M%nnz

CALL pardiso(M%pardisoPtr, 1, 1, mtype, phase, nr, csrVal, csrRowPtr, csrColIdx, M%pardisoPermute, 1,               &
             M%pardisoParam, 0, b, x, ierr)

CALL MPI_BARRIER(PE%MPI_CMFD_COMM, ierr)

Tend = nTracer_dclock(FALSE, FALSE)
TimeChk%AxBTime = TimeChk%AxBTime + (Tend - Tbeg)

END SUBROUTINE

SUBROUTINE pardisoDelete(M)

IMPLICIT NONE

TYPE(CSR_DOUBLE) :: M
INTEGER :: ierr
INTEGER :: nr, nc, nnz
INTEGER, PARAMETER :: mtype = 1, phase = -1
INTEGER, POINTER :: csrRowPtr(:), csrColIdx(:)
REAL, POINTER :: csrVal(:)
REAL, POINTER :: x(:), b(:)   !--- Dummy

csrVal => M%csrVal
csrRowPtr => M%csrRowPtr
csrColIdx => M%csrColIdx
nr = M%nr
nc = M%nc
nnz = M%nnz

CALL pardiso(M%pardisoPtr, 1, 1, mtype, phase, nr, csrVal, csrRowPtr, csrColIdx, M%pardisoPermute, 1,               &
             M%pardisoParam, 0, b, x, ierr)

END SUBROUTINE

! ============================================================ !
!     Anderson Acceleration of the Jacobi Iterative Method     !
!   Phanisri P. Pratapa, Phanish Suryanarayana, John E. Pask   !
! ============================================================ !

SUBROUTINE AAJ(CMFD, Jacobi, x, b)
USE PE_MOD,         ONLY : PE
USE TIMER,          ONLY : nTracer_dclock,      TimeChk
IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
TYPE(mklJacobi_Type) :: Jacobi
REAL :: x(*), b(*)

INTEGER :: p = 3, k, i, iter, itermin = 6, ierr
INTEGER :: n, nxy, nzCMFD
INTEGER :: bottomRange(2), topRange(2)
REAL, POINTER :: Xk(:, :), Fk(:, :), Fr(:, :)
REAL, POINTER :: Xkp(:, :), Fkp(:, :), Frp(:, :), invFr(:, :)
REAL, POINTER :: M(:, :), gamma(:), v(:), f(:), fp(:), xp(:)
REAL :: w = 0.5, beta = 0.5, err, tol = 0.1
REAL :: r0, r1
REAL :: Tbeg, Tend

Tbeg = nTracer_dclock(FALSE, FALSE)

nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
bottomRange = CMFD%bottomRange
topRange = CMFD%topRange
n = nxy * nzCMFD

ALLOCATE(f(n), fp(n), xp(n))

k = 1; iter = 0   !--- Initial Step

CALL dcopy(n, x, 1, xp, 1)

CALL MatMulMPIFast(Jacobi%M, Jacobi%offDiag, xp, f, bottomRange, topRange)
CALL vdsub(n, b, f, f); CALL vdmul(n, Jacobi%invDiag, f, x); CALL vdsub(n, x, xp, f)

r0 = dnrm2(n, f, 1)

ALLOCATE(Xk(n, k))
CALL vdsub(n, x, xp, xp)
CALL dcopy(n, xp, 1, Xk(:, k), 1)

DO WHILE (TRUE)
  
  CALL dcopy(n, x, 1, xp, 1)
  CALL dcopy(n, f, 1, fp, 1)

  !--- Compute f(x)

  CALL MatMulMPIFast(Jacobi%M, Jacobi%offDiag, xp, f, bottomRange, topRange)
  CALL vdsub(n, b, f, f); CALL vdmul(n, Jacobi%invDiag, f, f); CALL vdsub(n, f, xp, f)
  
  !--- Check Convergence
  
  r1 = dnrm2(n, f, 1); err = r1 / r0
  IF (err .LE. tol .AND. iter .GT. itermin) EXIT
  
  !--- Construct F
  
  IF (k .GT. 1) THEN
    Fkp => Fk; ALLOCATE(Fk(n, k))
    CALL dcopy(n * (k - 1), Fkp, 1, Fk, 1); DEALLOCATE(Fkp)
  ELSE
    ALLOCATE(Fk(n, k))
  ENDIF
  CALL vdsub(n, f, fp, fp)
  CALL dcopy(n, fp, 1, Fk(:, k), 1)
  
  !--- Compute F'F
  
  IF (k .GT. 1) THEN
    Frp => Fr; ALLOCATE(Fr(k, k))
    Fr(1 : k - 1, 1 : k - 1) = Frp; DEALLOCATE(Frp)
  ELSE
    ALLOCATE(Fr(k, k))
  ENDIF
  
  DO i = 1, k
    Fr(i, k) = ddot(n, Fk(:, i), 1, Fk(:, k), 1)
    Fr(k, i) = Fr(i, k)
  ENDDO

  !--- Alternating Anderson-Jacobi
  
  IF (mod(k, p) .NE. 0) THEN
    CALL WeightedJacobi()
  ELSE
    CALL AndersonJacobi()
  ENDIF
  
  !--- Iteration Status Update

  k = k + 1; iter = iter + 1
  
  !--- Construct X
  
  Xkp => Xk; ALLOCATE(Xk(n, k))
  CALL dcopy(n * (k - 1), Xkp, 1, Xk, 1); DEALLOCATE(Xkp)
  CALL vdsub(n, x, xp, xp)
  CALL dcopy(n, xp, 1, Xk(:, k), 1)
  
ENDDO

DEALLOCATE(Xk, Fk, Fr)
DEALLOCATE(f, fp, xp)

CALL MPI_BARRIER(PE%MPI_CMFD_COMM, ierr)

Tend = nTracer_dclock(FALSE, FALSE)
TimeChk%AxBTime = TimeChk%AxBTime + (Tend - Tbeg)

CONTAINS

SUBROUTINE WeightedJacobi()

!--- Solve for x

CALL daxpy(n, w, f, 1, x, 1)
  
END SUBROUTINE

SUBROUTINE AndersonJacobi()

!--- Compute inv(F'F)

ALLOCATE(invFr(k, k))

IF (mklCntl%lDcpl) THEN
  CALL dcopy(k ** 2, Fr, 1, invFr, 1)
ELSE
  CALL MPI_ALLREDUCE(Fr, invFr, k ** 2, MPI_DOUBLE_PRECISION, MPI_SUM, PE%MPI_CMFD_COMM, ierr)
ENDIF

CALL MatInv(invFr, k)

!--- Anderson Extrapolation

ALLOCATE(M(n, k), gamma(k), v(k))

CALL dcopy(n * k, Xk, 1, M, 1)
CALL daxpy(n * k, beta, Fk, 1, M, 1)
CALL dgemv('t', n, k, 1.0, Fk, n, f, 1, 0.0, v, 1)
CALL dgemv('n', k, k, 1.0, invFr, k, v, 1, 0.0, gamma, 1)
CALL dgemv('n', n, k, -1.0, M, n, gamma, 1, 1.0, x, 1)
CALL daxpy(n, beta, f, 1, x, 1)

DEALLOCATE(M, gamma, v)

DEALLOCATE(invFr)

END SUBROUTINE

SUBROUTINE MatInv(A, n)
            
IMPLICIT NONE

REAL :: A(*)
INTEGER :: n

INTEGER :: ipiv(n), ierr, worksize
REAL, POINTER :: workspace(:)
    
CALL dgetrf(n, n, A, n, ipiv, ierr)

worksize = n * 64; ALLOCATE(workspace(worksize))
CALL dgetri(n, A, n, ipiv, workspace, worksize, ierr)
DEALLOCATE(workspace)
        
END SUBROUTINE

END SUBROUTINE

END MODULE
    
#endif