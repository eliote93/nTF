#include <CUDADEFINES.h>

#ifdef __PGI

MODULE CUDA_SOLVER

USE CUDA_MASTER
USE CUDA_UTIL
IMPLICIT NONE

REAL(GPU_CMFD_PRECISION), PARAMETER, PRIVATE :: zero = 0.0, one = 1.0

CONTAINS

SUBROUTINE cuInnerSolve(cuCMFD, cuDevice, tol)
USE TIMER,          ONLY : nTracer_dclock,      TimeChk
IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
REAL :: Tbeg, Tend
REAL :: tol
LOGICAL :: lConv

Tbeg = nTracer_dclock(.FALSE., .FALSE.)

#if (GPU_CMFD_PRECISION == 4)
CALL cuIterativeRefinement(cuCMFD, cuDevice, tol)
#elif (GPU_CMFD_PRECISION == 8)
IF (cuCntl%lNatural) THEN
  CALL cuBiCGSTAB(cuCMFD, cuDevice, cuCMFD%phis8, cuCMFD%src8, tol)
ELSE
  CALL cuRedblackSOR(cuCMFD, cuDevice, cuCMFD%phis8, cuCMFD%src8, tol)
ENDIF
#endif

Tend = nTracer_dclock(.FALSE., .FALSE.)
TimeChk%AxBTime = TimeChk%AxBTime + (Tend - Tbeg)

END SUBROUTINE

SUBROUTINE cuInnerSolve_outInfo(cuCMFD, cuDevice, tol, finRes, Res0, initer)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
REAL :: tol

REAL :: finRes, Res0
INTEGER :: initer

LOGICAL :: lConv
#if (GPU_CMFD_PRECISION == 4)
CALL cuIterativeRefinement_outInfo(cuCMFD, cuDevice, tol, finRes, Res0, initer)
#elif (GPU_CMFD_PRECISION == 8)
IF (cuCntl%lNatural) THEN
  CALL cuBiCGSTAB_outinfo(cuCMFD, cuDevice, cuCMFD%phis8, cuCMFD%src8, tol, finRes, Res0, initer)
  !IF (cuCntl%CMFDSolver .EQ. 2) CALL cuAAJ(cuCMFD, cuDevice, cuCMFD%phis8, cuCMFD%src8, tol, lConv)
ELSE
  CALL cuRedblackSOR(cuCMFD, cuDevice, cuCMFD%phis8, cuCMFD%src8, tol)
ENDIF
#endif

END SUBROUTINE

#if (GPU_CMFD_PRECISION == 4)

SUBROUTINE cuIterativeRefinement_outinfo(cuCMFD, cuDevice, tol, finRes, Res0, initer)
USE TIMER,          ONLY : nTracer_dclock,      TimeChk
IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
REAL :: tol
REAL :: finRes, Res0
INTEGER :: initer

REAL(4), POINTER, DEVICE :: r4(:, :), d4(:, :)
REAL(8), POINTER, DEVICE :: r8(:, :), d8(:, :)
INTEGER :: ierr
INTEGER :: n, ng, nxy, nzCMFD
INTEGER :: bottomRange(2), topRange(2), bottomRangeRB(2, 2), topRangeRB(2, 2)
REAL :: Tbeg, Tend

Tbeg = nTracer_dclock(.FALSE., .FALSE.)

ng = cuCMFD%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
bottomRange = cuCMFD%bottomRange
bottomRangeRB = cuCMFD%bottomRangeRB
topRange = cuCMFD%topRange
topRangeRB = cuCMFD%topRangeRB
n = ng * nxy * nzCMFD

ALLOCATE(r4(ng, nxy * nzCMFD), r8(ng, nxy * nzCMFD)); r4 = 0.0
ALLOCATE(d4(ng, nxy * nzCMFD), d8(ng, nxy * nzCMFD)); d4 = 0.0

!--- Compute Residual
IF (cuCntl%lNatural) THEN
  CALL cuMatMul3D(cuCMFD%M, cuCMFD%offDiag8, cuCMFD%phis8, r8, bottomRange, topRange,                               &
                  cuDevice%mySparseHandle, cuDevice%myStream, MPI_CUDA_COMM)
ELSE
  CALL cuMatMul3D(cuCMFD%M, cuCMFD%offDiag8, cuCMFD%phis8, r8, bottomRangeRB, topRangeRB,                           &
                  cuDevice%mySparseHandle, cuDevice%myStream, MPI_CUDA_COMM)
ENDIF
CALL cuVectorOp('-', n, cuCMFD%src8, r8, r8, cuDevice%myStream)

!--- Find Approximate Correction
CALL cuTypecast(n, r8, r4, cuDevice%myStream)
IF (cuCntl%lNatural) THEN
  CALL cuBiCGSTAB_outinfo(cuCMFD, cuDevice, d4, r4, tol, finRes, Res0, initer)
!  IF (cuCntl%CMFDSolver .EQ. 2) CALL cuAAJ(cuCMFD, cuDevice, d4, r4, tol)
ELSE
  CALL cuRedblackSOR(cuCMFD, cuDevice, d4, r4, tol)
ENDIF
CALL cuTypecast(n, d4, d8, cuDevice%myStream)

!--- Add Correction
CALL cuVectorOp('+', n, cuCMFD%phis8, d8, cuCMFD%phis8, cuDevice%myStream)

DEALLOCATE(r4, r8)
DEALLOCATE(d4, d8)

ierr = cudaStreamSynchronize(cuDevice%myStream)

Tend = nTracer_dclock(.FALSE., .FALSE.)
TimeChk%AxBTime = TimeChk%AxBTime + (Tend - Tbeg)

END SUBROUTINE

SUBROUTINE cuIterativeRefinement(cuCMFD, cuDevice, tol)
IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
REAL :: tol

REAL(4), POINTER, DEVICE :: r4(:, :), d4(:, :)
REAL(8), POINTER, DEVICE :: r8(:, :), d8(:, :)
INTEGER :: ierr
INTEGER :: n, ng, nxy, nzCMFD
INTEGER :: bottomRange(2), topRange(2), bottomRangeRB(2, 2), topRangeRB(2, 2)

ng = cuCMFD%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
bottomRange = cuCMFD%bottomRange
bottomRangeRB = cuCMFD%bottomRangeRB
topRange = cuCMFD%topRange
topRangeRB = cuCMFD%topRangeRB
n = ng * nxy * nzCMFD

ALLOCATE(r4(ng, nxy * nzCMFD), r8(ng, nxy * nzCMFD)); r4 = 0.0
ALLOCATE(d4(ng, nxy * nzCMFD), d8(ng, nxy * nzCMFD)); d4 = 0.0

!--- Compute Residual
IF (cuCntl%lNatural) THEN
  CALL cuMatMul3D(cuCMFD%M, cuCMFD%offDiag8, cuCMFD%phis8, r8, bottomRange, topRange,                               &
                  cuDevice%mySparseHandle, cuDevice%myStream, MPI_CUDA_COMM)
ELSE
  CALL cuMatMul3D(cuCMFD%M, cuCMFD%offDiag8, cuCMFD%phis8, r8, bottomRangeRB, topRangeRB,                           &
                  cuDevice%mySparseHandle, cuDevice%myStream, MPI_CUDA_COMM)
ENDIF
CALL cuVectorOp('-', n, cuCMFD%src8, r8, r8, cuDevice%myStream)

!--- Find Approximate Correction
CALL cuTypecast(n, r8, r4, cuDevice%myStream)
IF (cuCntl%lNatural) THEN
  CALL cuBiCGSTAB(cuCMFD, cuDevice, d4, r4, tol)
ELSE
  CALL cuRedblackSOR(cuCMFD, cuDevice, d4, r4, tol)
ENDIF
CALL cuTypecast(n, d4, d8, cuDevice%myStream)

!--- Add Correction
CALL cuVectorOp('+', n, cuCMFD%phis8, d8, cuCMFD%phis8, cuDevice%myStream)

DEALLOCATE(r4, r8)
DEALLOCATE(d4, d8)

ierr = cudaStreamSynchronize(cuDevice%myStream)

END SUBROUTINE

#endif

SUBROUTINE cuBiCGSTAB_outinfo(cuCMFD, cuDevice, x, b, tol, finRes, Res0, initer)
USE CUDA_PCR
IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
REAL(GPU_CMFD_PRECISION), DEVICE :: x(*), b(*)
REAL :: err, tol
REAL :: finRes, Res0
INTEGER :: initer

INTEGER :: i, ierr, iter, itermin = 3, itermax = 200 ! 50
INTEGER :: nr, nc, nnz
INTEGER :: bottomRange(2), topRange(2)
REAL(GPU_CMFD_PRECISION) :: rho0, rho1, w0, w1, norm0, norm1, dot0, dot1, alpha, beta
REAL(GPU_CMFD_PRECISION), ALLOCATABLE, DEVICE :: h(:), s(:), t(:), y(:), z(:)
REAL(GPU_CMFD_PRECISION), ALLOCATABLE, DEVICE :: r0(:), r1(:), rhat(:), v0(:), v1(:), p0(:), p1(:)
LOGICAL :: lConv

REAL(GPU_CMFD_PRECISION), POINTER :: tester(:)

nr = cuCMFD%M%nr
nc = cuCMFD%M%nc
nnz = cuCMFD%M%nnz
bottomRange = cuCMFD%bottomRange
topRange = cuCMFD%topRange

! ======================================================================== !
!   CUDA Preconditioned Bi-Conjugate Gradient Stablized(BiCGSTAB) Solver   !
!   https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method   !
! ======================================================================== !

ALLOCATE(h(nc), s(nc), t(nc), y(nc), z(nc))
ALLOCATE(r0(nc), r1(nc), rhat(nc))
ALLOCATE(v0(nc), v1(nc))
ALLOCATE(p0(nc), p1(nc))


!--- 1 ----------------------
CALL cuMatMul3D(cuCMFD%M, cuCMFD%offDiag, x, r1, bottomRange, topRange, cuDevice%mySparseHandle,                    &
                cuDevice%myStream, MPI_CUDA_COMM)
CALL cuVectorOp('-', nc, b, r1, r1, cuDevice%myStream)
!--- 2 ----------------------
ierr = cublasPcopy_v2(cuDevice%myblasHandle, nc, r1, 1, rhat, 1)
!--- 3 ----------------------
IF (cuCntl%lDcpl) THEN
  IF (cuCntl%DcplLv .EQ. 1) THEN
    norm0 = normMulti(r1, nc, MPI_CUDA_COMM, cuDevice%myblasHandle, .FALSE.)
  ELSEIF (cuCntl%DcplLv .EQ. 2) THEN
    ierr = cublasPnrm2_v2(cuDevice%myblasHandle, nc, r1, 1, norm0)
  ENDIF
ELSE
  norm0 = normMulti(r1, nc, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
ENDIF
rho1 = 1.0; w1 = 1.0; alpha = 1.0
!--- 4 ----------------------
v1 = 0.0; p1 = 0.0
!--- 5 ----------------------
lConv = .FALSE.; iter = 0
DO WHILE (.NOT. lConv)
  !--- 5.0 ------------------
  iter = iter + 1
  rho0 = rho1
  w0 = w1
  ierr = cublasPcopy_v2(cuDevice%myblasHandle, nc, r1, 1, r0, 1)
  ierr = cublasPcopy_v2(cuDevice%myblasHandle, nc, p1, 1, p0, 1)
  ierr = cublasPcopy_v2(cuDevice%myblasHandle, nc, v1, 1, v0, 1)
  !--- 5.1 ------------------
  IF (cuCntl%lDcpl) THEN
    IF (cuCntl%DcplLv .EQ. 1) THEN
      rho1 = dotMulti(rhat, r0, nc, MPI_CUDA_COMM, cuDevice%myblasHandle, .FALSE.)
    ELSEIF (cuCntl%DcplLv .EQ. 2) THEN
      ierr = cublasPdot_v2(cuDevice%myblasHandle, nc, rhat, 1, r0, 1, rho1)
      ierr = cudaStreamSynchronize(cuDevice%myStream)
    ENDIF
  ELSE
    rho1 = dotMulti(rhat, r0, nc, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
  ENDIF

  !--- 5.2 ------------------
  beta = (rho1 / rho0) * (alpha / w0)
  !--- 5.3 ------------------
  ierr = cublasPcopy_v2(cuDevice%myblasHandle, nc, r0, 1, p1, 1)
  ierr = cublasPaxpy_v2(cuDevice%myblasHandle, nc, -w0, v0, 1, p0, 1)
  ierr = cublasPaxpy_v2(cuDevice%myblasHandle, nc, beta, p0, 1, p1, 1)
  !--- 5.4 ------------------
  CALL cuApplyPreconditioner(cuCMFD, cuDevice, p1, y)
  !--- 5.5 ------------------
  CALL cuMatMul3D(cuCMFD%M, cuCMFD%offDiag, y, v1, bottomRange, topRange, cuDevice%mySparseHandle,                  &
                  cuDevice%myStream, MPI_CUDA_COMM)
  !--- 5.6 ------------------
  IF (cuCntl%lDcpl) THEN
    IF (cuCntl%DcplLv .EQ. 1) THEN
      dot1 = dotMulti(rhat, v1, nc, MPI_CUDA_COMM, cuDevice%myblasHandle, .FALSE.)
    ELSEIF (cuCntl%DcplLv .EQ. 2) THEN
      ierr = cublasPdot_v2(cuDevice%myblasHandle, nc, rhat, 1, v1, 1, dot1)
      ierr = cudaStreamSynchronize(cuDevice%myStream)
    ENDIF
  ELSE
    dot1 = dotMulti(rhat, v1, nc, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
  ENDIF
  alpha = rho1 / dot1
  !--- 5.7 ------------------
  ierr = cublasPcopy_v2(cuDevice%myblasHandle, nc, x, 1, h, 1)
  ierr = cublasPaxpy_v2(cuDevice%myblasHandle, nc, alpha, y, 1, h, 1)
  !--- 5.9 ------------------
  ierr = cublasPcopy_v2(cuDevice%myblasHandle, nc, r0, 1, s, 1)
  ierr = cublasPaxpy_v2(cuDevice%myblasHandle, nc, -alpha, v1, 1, s, 1)
  !--- 5.10 -----------------
  CALL cuApplyPreconditioner(cuCMFD, cuDevice, s, z)
  !--- 5.11 -----------------
  CALL cuMatMul3D(cuCMFD%M, cuCMFD%offDiag, z, t, bottomRange, topRange, cuDevice%mySparseHandle,                   &
                  cuDevice%myStream, MPI_CUDA_COMM)
  !--- 5.12 -----------------
  IF (cuCntl%lDcpl) THEN
    IF (cuCntl%DcplLv .EQ. 1) THEN
      dot0 = dotMulti(t, s, nc, MPI_CUDA_COMM, cuDevice%myblasHandle, .FALSE.)
      dot1 = dotMulti(t, t, nc, MPI_CUDA_COMM, cuDevice%myblasHandle, .FALSE.)
    ELSEIF (cuCntl%DcplLv .EQ. 2) THEN
      ierr = cublasPdot_v2(cuDevice%myblasHandle, nc, t, 1, s, 1, dot0)
      ierr = cublasPdot_v2(cuDevice%myblasHandle, nc, t, 1, t, 1, dot1)
      ierr = cudaStreamSynchronize(cuDevice%myStream)
    ENDIF
  ELSE
    dot0 = dotMulti(t, s, nc, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
    dot1 = dotMulti(t, t, nc, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
  ENDIF
  w1 = dot0 / dot1
  !--- 5.13 -----------------
  ierr = cublasPcopy_v2(cuDevice%myblasHandle, nc, h, 1, x, 1)
  ierr = cublasPaxpy_v2(cuDevice%myblasHandle, nc, w1, z, 1, x, 1)
  !--- 5.15 -----------------
  ierr = cublasPcopy_v2(cuDevice%myblasHandle, nc, s, 1, r1, 1)
  ierr = cublasPaxpy_v2(cuDevice%myblasHandle, nc, -w1, t, 1, r1, 1)
  !--- Convergence ----------
  IF (cuCntl%lDcpl) THEN
    IF (cuCntl%DcplLv .EQ. 1) THEN
      norm1 = normMulti(r1, nc, MPI_CUDA_COMM, cuDevice%myblasHandle, .FALSE.)
    ELSEIF (cuCntl%DcplLv .EQ. 2) THEN
      ierr = cublasPnrm2_v2(cuDevice%myblasHandle, nc, r1, 1, norm1)
      ierr = cudaStreamSynchronize(cuDevice%myStream)
    ENDIF
  ELSE
    norm1 = normMulti(r1, nc, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
  ENDIF
  err = norm1 / norm0; lConv = (err .LE. tol) .AND. (iter .GE. itermin)
  lConv = lConv .OR. (iter .GE. itermax)
  !IF(ABS(rho1) .LT. 1.e-10) THEN
  !  lConv = .TRUE.
  !END IF
  !IF(cuDevice%myzb .EQ. 1) WRITE(108,'(i3,4es14.6)') iter, norm1, rho1, dot1, norm0
ENDDO

initer = iter
finRes = norm1
Res0 = norm0

DEALLOCATE(h, s, t, y, z)
DEALLOCATE(r0, r1, rhat)
DEALLOCATE(v0, v1)
DEALLOCATE(p0, p1)

END SUBROUTINE

SUBROUTINE cuBiCGSTAB(cuCMFD, cuDevice, x, b, tol)
USE CUDA_PCR
IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
REAL(GPU_CMFD_PRECISION), DEVICE :: x(*), b(*)
REAL :: err, tol

INTEGER :: i, ierr, iter, itermin = 3, itermax = 100
INTEGER :: nr, nc, nnz
INTEGER :: bottomRange(2), topRange(2)
REAL(GPU_CMFD_PRECISION) :: rho0, rho1, w0, w1, norm0, norm1, dot0, dot1, alpha, beta
REAL(GPU_CMFD_PRECISION), ALLOCATABLE, DEVICE :: h(:), s(:), t(:), y(:), z(:)
REAL(GPU_CMFD_PRECISION), ALLOCATABLE, DEVICE :: r0(:), r1(:), rhat(:), v0(:), v1(:), p0(:), p1(:)
LOGICAL :: lConv

nr = cuCMFD%M%nr
nc = cuCMFD%M%nc
nnz = cuCMFD%M%nnz
bottomRange = cuCMFD%bottomRange
topRange = cuCMFD%topRange

! ======================================================================== !
!   CUDA Preconditioned Bi-Conjugate Gradient Stablized(BiCGSTAB) Solver   !
!   https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method   !
! ======================================================================== !

ALLOCATE(h(nc), s(nc), t(nc), y(nc), z(nc))
ALLOCATE(r0(nc), r1(nc), rhat(nc))
ALLOCATE(v0(nc), v1(nc))
ALLOCATE(p0(nc), p1(nc))

!--- 1 ----------------------
CALL cuMatMul3D(cuCMFD%M, cuCMFD%offDiag, x, r1, bottomRange, topRange, cuDevice%mySparseHandle,                    &
                cuDevice%myStream, MPI_CUDA_COMM)
CALL cuVectorOp('-', nc, b, r1, r1, cuDevice%myStream)
!--- 2 ----------------------
ierr = cublasPcopy_v2(cuDevice%myblasHandle, nc, r1, 1, rhat, 1)
!--- 3 ----------------------
IF (cuCntl%lDcpl) THEN
  IF (cuCntl%DcplLv .EQ. 1) THEN
    norm0 = normMulti(r1, nc, MPI_CUDA_COMM, cuDevice%myblasHandle, .FALSE.)
  ELSEIF (cuCntl%DcplLv .EQ. 2) THEN
    ierr = cublasPnrm2_v2(cuDevice%myblasHandle, nc, r1, 1, norm0)
  ENDIF
ELSE
  norm0 = normMulti(r1, nc, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
ENDIF
rho1 = 1.0; w1 = 1.0; alpha = 1.0
!--- 4 ----------------------
v1 = 0.0; p1 = 0.0
!--- 5 ----------------------
lConv = .FALSE.; iter = 0
DO WHILE (.NOT. lConv)
  !--- 5.0 ------------------
  iter = iter + 1
  rho0 = rho1
  w0 = w1
  ierr = cublasPcopy_v2(cuDevice%myblasHandle, nc, r1, 1, r0, 1)
  ierr = cublasPcopy_v2(cuDevice%myblasHandle, nc, p1, 1, p0, 1)
  ierr = cublasPcopy_v2(cuDevice%myblasHandle, nc, v1, 1, v0, 1)
  !--- 5.1 ------------------
  IF (cuCntl%lDcpl) THEN
    IF (cuCntl%DcplLv .EQ. 1) THEN
      rho1 = dotMulti(rhat, r0, nc, MPI_CUDA_COMM, cuDevice%myblasHandle, .FALSE.)
    ELSEIF (cuCntl%DcplLv .EQ. 2) THEN
      ierr = cublasPdot_v2(cuDevice%myblasHandle, nc, rhat, 1, r0, 1, rho1)
      ierr = cudaStreamSynchronize(cuDevice%myStream)
    ENDIF
  ELSE
    rho1 = dotMulti(rhat, r0, nc, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
  ENDIF
  !--- 5.2 ------------------
  beta = (rho1 / rho0) * (alpha / w0)
  !--- 5.3 ------------------
  ierr = cublasPcopy_v2(cuDevice%myblasHandle, nc, r0, 1, p1, 1)
  ierr = cublasPaxpy_v2(cuDevice%myblasHandle, nc, -w0, v0, 1, p0, 1)
  ierr = cublasPaxpy_v2(cuDevice%myblasHandle, nc, beta, p0, 1, p1, 1)
  !--- 5.4 ------------------
  CALL cuApplyPreconditioner(cuCMFD, cuDevice, p1, y)
  !--- 5.5 ------------------
  CALL cuMatMul3D(cuCMFD%M, cuCMFD%offDiag, y, v1, bottomRange, topRange, cuDevice%mySparseHandle,                  &
                  cuDevice%myStream, MPI_CUDA_COMM)
  !--- 5.6 ------------------
  IF (cuCntl%lDcpl) THEN
    IF (cuCntl%DcplLv .EQ. 1) THEN
      dot1 = dotMulti(rhat, v1, nc, MPI_CUDA_COMM, cuDevice%myblasHandle, .FALSE.)
    ELSEIF (cuCntl%DcplLv .EQ. 2) THEN
      ierr = cublasPdot_v2(cuDevice%myblasHandle, nc, rhat, 1, v1, 1, dot1)
      ierr = cudaStreamSynchronize(cuDevice%myStream)
    ENDIF
  ELSE
    dot1 = dotMulti(rhat, v1, nc, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
  ENDIF
  alpha = rho1 / dot1
  !--- 5.7 ------------------
  ierr = cublasPcopy_v2(cuDevice%myblasHandle, nc, x, 1, h, 1)
  ierr = cublasPaxpy_v2(cuDevice%myblasHandle, nc, alpha, y, 1, h, 1)
  !--- 5.9 ------------------
  ierr = cublasPcopy_v2(cuDevice%myblasHandle, nc, r0, 1, s, 1)
  ierr = cublasPaxpy_v2(cuDevice%myblasHandle, nc, -alpha, v1, 1, s, 1)
  !--- 5.10 -----------------
  CALL cuApplyPreconditioner(cuCMFD, cuDevice, s, z)
  !--- 5.11 -----------------
  CALL cuMatMul3D(cuCMFD%M, cuCMFD%offDiag, z, t, bottomRange, topRange, cuDevice%mySparseHandle,                   &
                  cuDevice%myStream, MPI_CUDA_COMM)
  !--- 5.12 -----------------
  IF (cuCntl%lDcpl) THEN
    IF (cuCntl%DcplLv .EQ. 1) THEN
      dot0 = dotMulti(t, s, nc, MPI_CUDA_COMM, cuDevice%myblasHandle, .FALSE.)
      dot1 = dotMulti(t, t, nc, MPI_CUDA_COMM, cuDevice%myblasHandle, .FALSE.)
    ELSEIF (cuCntl%DcplLv .EQ. 2) THEN
      ierr = cublasPdot_v2(cuDevice%myblasHandle, nc, t, 1, s, 1, dot0)
      ierr = cublasPdot_v2(cuDevice%myblasHandle, nc, t, 1, t, 1, dot1)
      ierr = cudaStreamSynchronize(cuDevice%myStream)
    ENDIF
  ELSE
    dot0 = dotMulti(t, s, nc, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
    dot1 = dotMulti(t, t, nc, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
  ENDIF
  w1 = dot0 / dot1
  !--- 5.13 -----------------
  ierr = cublasPcopy_v2(cuDevice%myblasHandle, nc, h, 1, x, 1)
  ierr = cublasPaxpy_v2(cuDevice%myblasHandle, nc, w1, z, 1, x, 1)
  !--- 5.15 -----------------
  ierr = cublasPcopy_v2(cuDevice%myblasHandle, nc, s, 1, r1, 1)
  ierr = cublasPaxpy_v2(cuDevice%myblasHandle, nc, -w1, t, 1, r1, 1)
  !--- Convergence ----------
  IF (cuCntl%lDcpl) THEN
    IF (cuCntl%DcplLv .EQ. 1) THEN
      norm1 = normMulti(r1, nc, MPI_CUDA_COMM, cuDevice%myblasHandle, .FALSE.)
    ELSEIF (cuCntl%DcplLv .EQ. 2) THEN
      ierr = cublasPnrm2_v2(cuDevice%myblasHandle, nc, r1, 1, norm1)
      ierr = cudaStreamSynchronize(cuDevice%myStream)
    ENDIF
  ELSE
    norm1 = normMulti(r1, nc, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
  ENDIF
  err = norm1 / norm0; lConv = (err .LE. tol) .AND. (iter .GE. itermin)
  lConv = lConv .OR. (iter .GE. itermax)
ENDDO

DEALLOCATE(h, s, t, y, z)
DEALLOCATE(r0, r1, rhat)
DEALLOCATE(v0, v1)
DEALLOCATE(p0, p1)

END SUBROUTINE

SUBROUTINE cuRedblackGaussSeidel(cuCMFD, cuDevice, x, b, tol)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
REAL :: tol

INTEGER(KIND = cuda_stream_kind) :: myStream
INTEGER :: i, ierr, nIter
INTEGER :: nr, nc, nnz, ng, nxy, nzCMFD
INTEGER :: rbBeg(2), rbEnd(2), bottomRange(2, 2), topRange(2, 2)
INTEGER, POINTER, DEVICE :: csrRowPtr(:), csrColIdx(:)
REAL(GPU_CMFD_PRECISION), POINTER, DEVICE :: csrVal(:)
REAL(GPU_CMFD_PRECISION), POINTER, DEVICE :: b(:, :), xr(:, :), xb(:, :)
REAL(GPU_CMFD_PRECISION), POINTER, DEVICE :: x(:, :), br(:, :), bb(:, :)
REAL(GPU_CMFD_PRECISION), POINTER, DEVICE :: t(:, :), tr(:, :), tb(:, :)
REAL(GPU_CMFD_PRECISION), POINTER, DEVICE :: e(:)
REAL :: rho, rho_prev, rho_err, err(2)

ng = cuCMFD%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
rbBeg = cuDevice%rbBeg
rbEnd = cuDevice%rbEnd
bottomRange = cuCMFD%bottomRangeRB
topRange = cuCMFD%topRangeRB

ALLOCATE(t(ng, nxy * nzCMFD))
ALLOCATE(e(ng * maxval(cuDevice%nxyzRB)))

xr => x(:, rbBeg(red) : rbEnd(red))
xb => x(:, rbBeg(black) : rbEnd(black))
br => b(:, rbBeg(red) : rbEnd(red))
bb => b(:, rbBeg(black) : rbEnd(black))
tr => t(:, rbBeg(red) : rbEnd(red))
tb => t(:, rbBeg(black) : rbEnd(black))

rho_err = 1.0
DO WHILE (rho_err .GT. 0.0001)

  !--- Solve for Red
  csrVal => cuCMFD%rbDiag(red)%d_csrVal
  csrRowPtr => cuCMFD%rbDiag(red)%d_csrRowPtr
  csrColIdx => cuCMFD%rbDiag(red)%d_csrColIdx
  nr = cuCMFD%rbM(red)%nr; nc = cuCMFD%rbM(red)%nc; nnz = cuCMFD%rbM(red)%nnz
  CALL cuMatMul3D(cuCMFD%rbM(red), cuCMFD%rOffDiag, xb, tr, bottomRange(:, red), topRange(:, red),                  &
                  cuDevice%mySparseHandle, cuDevice%myStream, MPI_CUDA_COMM)
  CALL cuVectorOp('-', nr, br, tr, tr, cuDevice%myStream)
  ierr = cusparsePcsrsv_solve(cuDevice%mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, one,                   &
							  cuCMFD%rbDiag(red)%descrLU(1), csrVal, csrRowPtr, csrColIdx,                          &
                              cuCMFD%rbDiag(red)%infoLU(1), tr, tr)
  CALL cuVectorOp('-', nr, xr, tr, e, cuDevice%myStream)
  ierr = cublasPcopy_v2(cuDevice%myblasHandle, nr, tr, 1, xr, 1)
  err(2) = dotMulti(e, e, nr, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)

  !--- Solve for Black
  csrVal => cuCMFD%rbDiag(black)%d_csrVal
  csrRowPtr => cuCMFD%rbDiag(black)%d_csrRowPtr
  csrColIdx => cuCMFD%rbDiag(black)%d_csrColIdx
  nr = cuCMFD%rbM(black)%nr; nc = cuCMFD%rbM(black)%nc; nnz = cuCMFD%rbM(black)%nnz
  CALL cuMatMul3D(cuCMFD%rbM(black), cuCMFD%bOffDiag, xr, tb, bottomRange(:, black), topRange(:, black),            &
                  cuDevice%mySparseHandle, cuDevice%myStream, MPI_CUDA_COMM)
  CALL cuVectorOp('-', nr, bb, tb, tb, cuDevice%myStream)
  ierr = cusparsePcsrsv_solve(cuDevice%mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, one,                   &
							  cuCMFD%rbDiag(black)%descrLU(1), csrVal, csrRowPtr, csrColIdx,                        &
                              cuCMFD%rbDiag(black)%infoLU(1), tb, tb)
  CALL cuVectorOp('-', nr, xb, tb, e, cuDevice%myStream)
  ierr = cublasPcopy_v2(cuDevice%myblasHandle, nr, tb, 1, xb, 1)
  err(2) = err(2) + dotMulti(e, e, nr, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)

  !--- Estimate Spectral Radius
  rho_prev = rho
  rho = sqrt(err(2)) / sqrt(err(1))
  rho_err = abs(rho - rho_prev)
  err(1) = err(2)

ENDDO

cuCMFD%wSOR = 2.0 / (1.0 + sqrt(1.0 - rho))
cuCMFD%nIter = INT(log(tol) / log(cuCMFD%wSOR - 1.0)) + 1

DEALLOCATE(t, e)

END SUBROUTINE

SUBROUTINE cuRedblackSOR(cuCMFD, cuDevice, x, b, tol)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
REAL :: tol

INTEGER(KIND = cuda_stream_kind) :: myStream
INTEGER :: i, ierr, iter, nIter
INTEGER :: nr, nc, nnz, ng, nxy, nzCMFD
INTEGER :: rbBeg(2), rbEnd(2), bottomRange(2, 2), topRange(2, 2)
INTEGER, POINTER, DEVICE :: csrRowPtr(:), csrColIdx(:)
REAL(GPU_CMFD_PRECISION), POINTER, DEVICE :: csrVal(:)
REAL(GPU_CMFD_PRECISION), POINTER, DEVICE :: b(:, :), xr(:, :), xb(:, :)
REAL(GPU_CMFD_PRECISION), POINTER, DEVICE :: x(:, :), br(:, :), bb(:, :)
REAL(GPU_CMFD_PRECISION), POINTER, DEVICE :: t(:, :), tr(:, :), tb(:, :)
REAL(GPU_CMFD_PRECISION) :: wSOR

IF (cuCMFD%lFirst) THEN
  CALL cuRedblackGaussSeidel(cuCMFD, cuDevice, x, b, tol)
  cuCMFD%lFirst = .FALSE.; RETURN
ENDIF

ng = cuCMFD%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
rbBeg = cuDevice%rbBeg
rbEnd = cuDevice%rbEnd
bottomRange = cuCMFD%bottomRangeRB
topRange = cuCMFD%topRangeRB

nIter = cuCMFD%nIter
wSOR = cuCMFD%wSOR

ALLOCATE(t(ng, nxy * nzCMFD))

xr => x(:, rbBeg(red) : rbEnd(red))
xb => x(:, rbBeg(black) : rbEnd(black))
br => b(:, rbBeg(red) : rbEnd(red))
bb => b(:, rbBeg(black) : rbEnd(black))
tr => t(:, rbBeg(red) : rbEnd(red))
tb => t(:, rbBeg(black) : rbEnd(black))

DO iter = 1, nIter

  !--- Solve for Red
  csrVal => cuCMFD%rbDiag(red)%d_csrVal
  csrRowPtr => cuCMFD%rbDiag(red)%d_csrRowPtr
  csrColIdx => cuCMFD%rbDiag(red)%d_csrColIdx
  nr = cuCMFD%rbM(red)%nr; nc = cuCMFD%rbM(red)%nc; nnz = cuCMFD%rbM(red)%nnz
  CALL cuMatMul3D(cuCMFD%rbM(red), cuCMFD%rOffDiag, xb, tr, bottomRange(:, red), topRange(:, red),                  &
                  cuDevice%mySparseHandle, cuDevice%myStream, MPI_CUDA_COMM)
  CALL cuVectorOp('-', nr, br, tr, tr, cuDevice%myStream)
  ierr = cusparsePcsrsv_solve(cuDevice%mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, wSOR,                  &
							  cuCMFD%rbDiag(red)%descrLU(1), csrVal, csrRowPtr, csrColIdx,                          &
                              cuCMFD%rbDiag(red)%infoLU(1), tr, tr)
  ierr = cublasPaxpy_v2(cuDevice%myblasHandle, nr, one - wSOR, xr, 1, tr, 1)
  ierr = cublasPcopy_v2(cuDevice%myblasHandle, nr, tr, 1, xr, 1)

  !--- Solve for Black
  csrVal => cuCMFD%rbDiag(black)%d_csrVal
  csrRowPtr => cuCMFD%rbDiag(black)%d_csrRowPtr
  csrColIdx => cuCMFD%rbDiag(black)%d_csrColIdx
  nr = cuCMFD%rbM(black)%nr; nc = cuCMFD%rbM(black)%nc; nnz = cuCMFD%rbM(black)%nnz
  CALL cuMatMul3D(cuCMFD%rbM(black), cuCMFD%bOffDiag, xr, tb, bottomRange(:, black), topRange(:, black),            &
                  cuDevice%mySparseHandle, cuDevice%myStream, MPI_CUDA_COMM)
  CALL cuVectorOp('-', nr, bb, tb, tb, cuDevice%myStream)
  ierr = cusparsePcsrsv_solve(cuDevice%mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, wSOR,                  &
							  cuCMFD%rbDiag(black)%descrLU(1), csrVal, csrRowPtr, csrColIdx,                        &
                              cuCMFD%rbDiag(black)%infoLU(1), tb, tb)
  ierr = cublasPaxpy_v2(cuDevice%myblasHandle, nr, one - wSOR, xb, 1, tb, 1)
  ierr = cublasPcopy_v2(cuDevice%myblasHandle, nr, tb, 1, xb, 1)

ENDDO

DEALLOCATE(t)

END SUBROUTINE

END MODULE

#endif
