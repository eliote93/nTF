#include <CUDADEFINES.h>

#ifdef __PGI

MODULE CUDA_PCR

USE CUDA_MASTER
USE CUDA_UTIL
IMPLICIT NONE

REAL(GPU_CMFD_PRECISION), PARAMETER :: zero = 0.0, one = 1.0

PRIVATE
public :: cuPrepareSPAI
PUBLIC :: cuSetPreconditioner, cuApplyPreconditioner

INTERFACE cuPrepareSPAI
  MODULE PROCEDURE cuPrepareDSPAI
  MODULE PROCEDURE cuPrepareSSPAI
END INTERFACE

CONTAINS

SUBROUTINE cuSetPreconditioner(cuCMFD, cuDevice)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice

IF (cuCntl%lSPAI) THEN
  CALL cuPrepareSPAI(cuCMFD%D, cuCMFD%SPAI)
ELSE
  CALL cuPrepareILU(cuCMFD%M, cuCMFD%ILU, cuDevice%mySparseHandle)
ENDIF

END SUBROUTINE

SUBROUTINE cuApplyPreconditioner(cuCMFD, cuDevice, x, y)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
REAL(GPU_CMFD_PRECISION), DEVICE :: x(*), y(*)

IF (cuCntl%lSPAI) THEN
  CALL cuApplySPAI(cuCMFD%SPAI, cuDevice%mySparseHandle, x, y)
ELSE
  CALL cuSolveILU(cuCMFD%ILU, cuDevice%mySparseHandle, x, y)
ENDIF

END SUBROUTINE

SUBROUTINE cuPrepareILU(M, ILU, mySparseHandle)

IMPLICIT NONE

TYPE(CSR_CMFD_PRECISION) :: M, ILU
TYPE(cusparseHandle) :: mySparseHandle
LOGICAL :: lFirst

REAL(GPU_CMFD_PRECISION), POINTER, DEVICE :: csrVal(:)
INTEGER, POINTER, DEVICE :: csrRowPtr(:), csrColIdx(:)
INTEGER :: ierr
INTEGER :: nr, nc, nnz

ILU = M
CALL finalizeCsr(ILU, .TRUE.)

nr = ILU%nr
nc = ILU%nc
nnz = ILU%nnz
csrVal => ILU%d_csrVal
csrRowPtr => ILU%d_csrRowPtr
csrColIdx => ILU%d_csrColIdx

ierr = cusparseDestroySolveAnalysisInfo(ILU%info)
ierr = cusparseDestroySolveAnalysisInfo(ILU%infoLU(1))
ierr = cusparseDestroySolveAnalysisInfo(ILU%infoLU(2))
ierr = cusparseCreateSolveAnalysisInfo(ILU%info)
ierr = cusparseCreateSolveAnalysisInfo(ILU%infoLU(1))
ierr = cusparseCreateSolveAnalysisInfo(ILU%infoLU(2))
ierr = cusparseSetMatDiagType(ILU%descrLU(1), CUSPARSE_DIAG_TYPE_UNIT)

ierr = cusparsePcsrsv_analysis(mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nnz,							&
  						       ILU%descr, csrVal, csrRowPtr, csrColIdx, ILU%info)

ierr = cusparsePcsrilu0(mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr,										&
						ILU%descr, csrVal, csrRowPtr, csrColIdx, ILU%info)

ierr = cusparsePcsrsv_analysis(mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nnz,							&
 						       ILU%descrLU(1), csrVal, csrRowPtr, csrColIdx, ILU%infoLU(1))
ierr = cusparsePcsrsv_analysis(mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nnz,							&
						       ILU%descrLU(2), csrVal, csrRowPtr, csrColIdx, ILU%infoLU(2))

END SUBROUTINE

SUBROUTINE cuSolveILU(ILU, mySparseHandle, x, y)

IMPLICIT NONE

TYPE(CSR_CMFD_PRECISION) :: ILU
TYPE(cusparseHandle) :: mySparseHandle
REAL(GPU_CMFD_PRECISION), DEVICE :: x(*), y(*)

REAL(GPU_CMFD_PRECISION), ALLOCATABLE, DEVICE :: t(:)
REAL(GPU_CMFD_PRECISION), POINTER, DEVICE :: csrVal(:)
INTEGER, POINTER, DEVICE :: csrRowPtr(:), csrColIdx(:)
INTEGER :: ierr
INTEGER :: nr, nc, nnz

nr = ILU%nr
nc = ILU%nc
nnz = ILU%nnz
csrVal => ILU%d_csrVal
csrRowPtr => ILU%d_csrRowPtr
csrColIdx => ILU%d_csrColIdx

ALLOCATE(t(nc))

ierr = cusparsePcsrsv_solve(mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, one,                              &
                            ILU%descrLU(1), csrVal, csrRowPtr, csrColIdx, ILU%infoLU(1), x, t)
ierr = cusparsePcsrsv_solve(mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, one,                              &
                            ILU%descrLU(2), csrVal, csrRowPtr, csrColIdx, ILU%infoLU(2), t, y)

DEALLOCATE(t)

END SUBROUTINE

! ============================================================================ !
!                  Sparse Approximate Inverse Preconditioner                   !
!    "Parallel Preconditioning with Sparse Approximate Inverses", MJ Grote     !
! ============================================================================ !

SUBROUTINE cuPrepareSSPAI(M, SPAI)
USE PE_MOD,         ONLY : PE
IMPLICIT NONE

TYPE(CSR_MIXED) :: M, SPAI

TYPE(CSR_MIXED) :: trM, trSPAI
TYPE(CSR_MIXED), POINTER :: sub_trM(:)
REAL, POINTER :: subM(:, :), tau(:), work(:), ek(:, :), Pk(:)
INTEGER, POINTER :: sub_colIdx(:), sub_csrRowIdx(:), sub_rowIdx(:)

REAL :: tmpVal
INTEGER :: tmpInt
INTEGER :: lwork, info, lda

INTEGER :: nThread
INTEGER :: nc, nr, nnz
INTEGER :: ic, ibeg, iend
INTEGER :: tid

nnz = M%nnz
nr = M%nr
nc = M%nc

nThread = PE%nThread

CALL transposeCsr(M, trM)
CALL createCsr(trSPAI, nnz, nc, nr, nThread)

!$OMP PARALLEL PRIVATE(tid)
tid = omp_get_thread_num() + 1
!$OMP DO SCHEDULE(STATIC)
DO ic = 1, nc
  CALL cuPrepareSSPAI_1Col(trM, trSPAI, ic, nc, tid)
END DO
!$OMP END DO
!$OMP END PARALLEL

CALL finalizeCsr(trSPAI, .FALSE.)
CALL transposeCsr(trSPAI, SPAI)
CALL finalizeCsr(SPAI, .TRUE.)

CALL destroyCsr(trSPAI)
CALL destroyCsr(trM)

END SUBROUTINE

SUBROUTINE cuPrepareDSPAI(M, SPAI)
USE PE_MOD,         ONLY : PE
IMPLICIT NONE

TYPE(CSR_DOUBLE) :: M, SPAI

TYPE(CSR_DOUBLE) :: trM, trSPAI
TYPE(CSR_DOUBLE), POINTER :: sub_trM(:)
REAL, POINTER :: subM(:, :), tau(:), work(:), ek(:, :), Pk(:)
INTEGER, POINTER :: sub_colIdx(:), sub_csrRowIdx(:), sub_rowIdx(:)

REAL :: tmpVal
INTEGER :: tmpInt
INTEGER :: lwork, info, lda

INTEGER :: nThread
INTEGER :: nc, nr, nnz
INTEGER :: ic, ibeg, iend
INTEGER :: tid

nnz = M%nnz
nr = M%nr
nc = M%nc

nThread = PE%nThread

CALL transposeCsr(M, trM)
CALL createCsr(trSPAI, nnz, nc, nr, nThread)

!$OMP PARALLEL PRIVATE(tid)
tid = omp_get_thread_num() + 1
!$OMP DO SCHEDULE(STATIC)
DO ic = 1, nc
  CALL cuPrepareDSPAI_1Col(trM, trSPAI, ic, nc, tid)
END DO
!$OMP END DO
!$OMP END PARALLEL

CALL finalizeCsr(trSPAI, .FALSE.)
CALL transposeCsr(trSPAI, SPAI)
CALL finalizeCsr(SPAI, .TRUE.)

CALL destroyCsr(trSPAI)
CALL destroyCsr(trM)

END SUBROUTINE

SUBROUTINE cuPrepareSSPAI_1Col(trM, trSPAI, i, nc, tid)

IMPLICIT NONE

TYPE(CSR_MIXED) :: trM, trSPAI
INTEGER :: i, nc, tid

TYPE(CSR_MIXED) :: sub_trM
REAL, POINTER :: subM(:,:), tau(:), work(:), ek(:,:), Pk(:)
INTEGER, POINTER :: sub_colIdx(:), sub_csrRowIdx(:), sub_rowIdx(:)

REAL :: tmpVal
INTEGER :: tmpInt

INTEGER :: ntau, kIdx
INTEGER :: sub_nc, sub_nr, max_subr, sub_nnz
INTEGER :: subr, subc
INTEGER :: icol, tr_col
INTEGER :: offset, info
INTEGER :: j, k

sub_nc = trM%csrRowPtr(i + 1) - trM%csrRowPtr(i)
ALLOCATE(sub_colIdx(sub_nc))
offset = trM%csrRowPtr(i) - 1
DO j = trM%csrRowPtr(i), trM%csrRowPtr(i + 1) - 1
  sub_colIdx(j - offset) = trM%csrColIdx(j)
ENDDO

max_subr = 0
DO j = 1, sub_nc
  icol = sub_colIdx(j)
  max_subr = max_subr + trM%csrRowPtr(icol + 1) - trM%csrRowPtr(icol)
ENDDO

!--- Save Transposed Submatrix in CSR format
CALL createCsr(sub_trM, max_subr, sub_nc, nc)
DO j = 1, sub_nc
  icol = sub_colIdx(j)
  DO k = trM%csrRowPtr(icol), trM%csrRowPtr(icol + 1) - 1
    tr_col = trM%csrColIdx(k)
    tmpVal = trM%csrVal(k)
    CALL pushCsr(sub_trM, tmpVal, j, tr_col)
  ENDDO
ENDDO
CALL finalizeCsr(sub_trM, .FALSE.)


!--- Sorting
sub_nnz = sub_trM%nnz
ALLOCATE(sub_csrRowIdx(sub_nnz), sub_rowIdx(sub_nnz))
DO j = 1, sub_nc
  DO k = sub_trM%csrRowPtr(j), sub_trM%csrRowPtr(j + 1) - 1
    sub_csrRowIdx(k) = j
  ENDDO
ENDDO

DO j = sub_nnz, 1, -1
  DO k = 1, j - 1
    IF(sub_trM%csrColIdx(k) .GT. sub_trM%csrColIdx(k + 1)) THEN
      tmpInt = sub_trM%csrColIdx(k)
      sub_trM%csrColIdx(k) = sub_trM%csrColIdx(k + 1)
      sub_trM%csrColIdx(k + 1) = tmpInt

      tmpInt = sub_csrRowIdx(k)
      sub_csrRowIdx(k) = sub_csrRowIdx(k + 1)
      sub_csrRowIdx(k + 1) = tmpInt

      tmpVal = sub_trM%csrVal(k)
      sub_trM%csrVal(k) = sub_trM%csrVal(k + 1)
      sub_trM%csrVal(k + 1) = tmpVal
    ENDIF
  ENDDO
ENDDO

sub_rowIdx(1) = 1
IF(sub_trM%csrColIdx(1) .EQ. i) kIdx =  1
DO j = 2, sub_nnz
  IF(sub_trM%csrColIdx(j) .EQ. sub_trM%csrColIdx(j - 1)) THEN
    sub_rowIdx(j) = sub_rowIdx(j - 1)
  ELSE
    sub_rowIdx(j) = sub_rowIdx(j - 1) + 1
  ENDIF
  IF(sub_trM%csrColIdx(j) .EQ. i) kIdx =  sub_rowIdx(j)
ENDDO
sub_nr = sub_rowIdx(sub_nnz)
ntau = min(sub_nr, sub_nc)

ALLOCATE(subM(sub_nr, sub_nc), ek(sub_nr, 1), tau(ntau), work(sub_nc * 80), Pk(sub_nc))
subM = 0.
ek = 0.
ek(kIdx, 1) = 1.

DO j = 1, sub_nnz
  subr = sub_rowIdx(j)
  subc = sub_csrRowIdx(j)
  subM(subr, subc) = sub_trM%csrVal(j)
ENDDO

CALL DGEQRF(sub_nr, sub_nc, subM, sub_nr, tau, work, sub_nc * 80, info)
CALL DORMQR('L', 't', sub_nr, 1, ntau, subM, sub_nr, tau, ek, sub_nr, work, sub_nc * 80, info)

DO j = sub_nc, 1, -1
  tmpVal = ek(j, 1)
  DO k = sub_nc, j + 1, -1
    tmpVal = tmpVal - subM(j, k) * Pk(k)
  ENDDO
  Pk(j) = tmpVal / subM(j, j)
ENDDO

DO j = 1, sub_nc
  icol = sub_colIdx(j)
  CALL pushCsr(trSPAI, Pk(j), i, icol, tid)
ENDDO

CALL destroyCsr(sub_trM)
DEALLOCATE(subM, ek, tau, work, Pk)
DEALLOCATE(sub_colIdx, sub_csrRowIdx, sub_rowIdx)

END SUBROUTINE

SUBROUTINE cuPrepareDSPAI_1Col(trM, trSPAI, i, nc, tid)

IMPLICIT NONE

TYPE(CSR_DOUBLE) :: trM, trSPAI
INTEGER :: i, nc, tid

TYPE(CSR_DOUBLE) :: sub_trM
REAL, POINTER :: subM(:,:), tau(:), work(:), ek(:,:), Pk(:)
INTEGER, POINTER :: sub_colIdx(:), sub_csrRowIdx(:), sub_rowIdx(:)

REAL :: tmpVal
INTEGER :: tmpInt

INTEGER :: ntau, kIdx
INTEGER :: sub_nc, sub_nr, max_subr, sub_nnz
INTEGER :: subr, subc
INTEGER :: icol, tr_col
INTEGER :: offset, info
INTEGER :: j, k

sub_nc = trM%csrRowPtr(i + 1) - trM%csrRowPtr(i)
ALLOCATE(sub_colIdx(sub_nc))
offset = trM%csrRowPtr(i) - 1
DO j = trM%csrRowPtr(i), trM%csrRowPtr(i + 1) - 1
  sub_colIdx(j - offset) = trM%csrColIdx(j)
ENDDO

max_subr = 0
DO j = 1, sub_nc
  icol = sub_colIdx(j)
  max_subr = max_subr + trM%csrRowPtr(icol + 1) - trM%csrRowPtr(icol)
ENDDO

!--- Save Transposed Submatrix in CSR format
CALL createCsr(sub_trM, max_subr, sub_nc, nc)
DO j = 1, sub_nc
  icol = sub_colIdx(j)
  DO k = trM%csrRowPtr(icol), trM%csrRowPtr(icol + 1) - 1
    tr_col = trM%csrColIdx(k)
    tmpVal = trM%csrVal(k)
    CALL pushCsr(sub_trM, tmpVal, j, tr_col)
  ENDDO
ENDDO
CALL finalizeCsr(sub_trM, .FALSE.)


!--- Sorting
sub_nnz = sub_trM%nnz
ALLOCATE(sub_csrRowIdx(sub_nnz), sub_rowIdx(sub_nnz))
DO j = 1, sub_nc
  DO k = sub_trM%csrRowPtr(j), sub_trM%csrRowPtr(j + 1) - 1
    sub_csrRowIdx(k) = j
  ENDDO
ENDDO

DO j = sub_nnz, 1, -1
  DO k = 1, j - 1
    IF(sub_trM%csrColIdx(k) .GT. sub_trM%csrColIdx(k + 1)) THEN
      tmpInt = sub_trM%csrColIdx(k)
      sub_trM%csrColIdx(k) = sub_trM%csrColIdx(k + 1)
      sub_trM%csrColIdx(k + 1) = tmpInt

      tmpInt = sub_csrRowIdx(k)
      sub_csrRowIdx(k) = sub_csrRowIdx(k + 1)
      sub_csrRowIdx(k + 1) = tmpInt

      tmpVal = sub_trM%csrVal(k)
      sub_trM%csrVal(k) = sub_trM%csrVal(k + 1)
      sub_trM%csrVal(k + 1) = tmpVal
    ENDIF
  ENDDO
ENDDO

sub_rowIdx(1) = 1
IF(sub_trM%csrColIdx(1) .EQ. i) kIdx =  1
DO j = 2, sub_nnz
  IF(sub_trM%csrColIdx(j) .EQ. sub_trM%csrColIdx(j - 1)) THEN
    sub_rowIdx(j) = sub_rowIdx(j - 1)
  ELSE
    sub_rowIdx(j) = sub_rowIdx(j - 1) + 1
  ENDIF
  IF(sub_trM%csrColIdx(j) .EQ. i) kIdx =  sub_rowIdx(j)
ENDDO
sub_nr = sub_rowIdx(sub_nnz)
ntau = min(sub_nr, sub_nc)

ALLOCATE(subM(sub_nr, sub_nc), ek(sub_nr, 1), tau(ntau), work(sub_nc * 80), Pk(sub_nc))
subM = 0.
ek = 0.
ek(kIdx, 1) = 1.

DO j = 1, sub_nnz
  subr = sub_rowIdx(j)
  subc = sub_csrRowIdx(j)
  subM(subr, subc) = sub_trM%csrVal(j)
ENDDO

CALL DGEQRF(sub_nr, sub_nc, subM, sub_nr, tau, work, sub_nc * 80, info)
CALL DORMQR('L', 't', sub_nr, 1, ntau, subM, sub_nr, tau, ek, sub_nr, work, sub_nc * 80, info)

DO j = sub_nc, 1, -1
  tmpVal = ek(j, 1)
  DO k = sub_nc, j + 1, -1
    tmpVal = tmpVal - subM(j, k) * Pk(k)
  ENDDO
  Pk(j) = tmpVal / subM(j, j)
ENDDO

DO j = 1, sub_nc
  icol = sub_colIdx(j)
  CALL pushCsr(trSPAI, Pk(j), i, icol, tid)
ENDDO

CALL destroyCsr(sub_trM)
DEALLOCATE(subM, ek, tau, work, Pk)
DEALLOCATE(sub_colIdx, sub_csrRowIdx, sub_rowIdx)

END SUBROUTINE

SUBROUTINE cuApplySPAI(SPAI, mySparseHandle, x, y)

IMPLICIT NONE

TYPE(CSR_CMFD_PRECISION) :: SPAI
TYPE(cusparseHandle) :: mySparseHandle
REAL(GPU_CMFD_PRECISION), DEVICE :: x(*), y(*)

REAL(GPU_CMFD_PRECISION), POINTER, DEVICE :: csrVal(:)
INTEGER, POINTER, DEVICE :: csrRowPtr(:), csrColIdx(:)
INTEGER :: nr, nc, nnz
INTEGER :: ierr

nr = SPAI%nr
nc = SPAI%nc
nnz = SPAI%nnz
csrVal => SPAI%d_csrVal
csrRowPtr => SPAI%d_csrRowPtr
csrColIdx => SPAI%d_csrColIdx

ierr = cusparsePcsrmv(mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, one, SPAI%descr,             &
                      csrVal, csrRowPtr, csrColIdx, x, zero, y)

END SUBROUTINE

END MODULE

#endif
