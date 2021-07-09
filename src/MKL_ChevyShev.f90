#include <defines.h>
! 3D CMFD Acceleration Modules with Intel MKL
! Chebyshev Acceleration Routines
#ifdef __INTEL_MKL
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

END MODULE MKL_CHEBYSHEV
#endif