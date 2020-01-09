#ifdef __INTEL_MKL
INCLUDE 'mkl_pardiso.f90'
#endif

MODULE CSRMATRIX

#ifdef __INTEL_MKL
USE MKL_PARDISO
#endif

#ifdef __PGI
USE CUSPARSE
#endif

IMPLICIT NONE

INTEGER, PARAMETER, PRIVATE :: buffer = 2

TYPE CSR_FLOAT
  REAL(4), POINTER :: csrVal(:)
  INTEGER, POINTER :: csrRowPtr(:), csrColIdx(:)
  !--- Parallel Generation Variables
  REAL(4), POINTER :: bufVal(:, :)
  INTEGER, POINTER :: bufColIdx(:, :)
#ifdef __PGI
  REAL(4), POINTER, DEVICE :: d_csrVal(:)
  INTEGER, POINTER, DEVICE :: d_csrRowPtr(:), d_csrColIdx(:)
  TYPE(cusparseSolveAnalysisInfo) :: info, infoLU(2)
  TYPE(cusparseMatDescr) :: descr, descrLU(2)
#endif
  INTEGER :: nr, nc, nnz, nnzOmp(64), nThread
  LOGICAL :: lAlloc = .FALSE., lFinalized = .FALSE., lDevice = .FALSE., lParallel = .FALSE.
#ifdef __INTEL_MKL
  !--- PARDISO Variables
  TYPE(MKL_PARDISO_HANDLE) :: pardisoPtr(64)
  INTEGER :: pardisoParam(64)
  INTEGER, POINTER :: pardisoPermute(:)
#endif
END TYPE

TYPE CSR_DOUBLE
  REAL(8), POINTER :: csrVal(:)
  INTEGER, POINTER :: csrRowPtr(:), csrColIdx(:)
  !--- Parallel Generation Variables
  REAL(8), POINTER :: bufVal(:, :)
  INTEGER, POINTER :: bufColIdx(:, :)
#ifdef __PGI
  REAL(8), POINTER, DEVICE :: d_csrVal(:)
  INTEGER, POINTER, DEVICE :: d_csrRowPtr(:), d_csrColIdx(:)
  TYPE(cusparseSolveAnalysisInfo) :: info, infoLU(2)
  TYPE(cusparseMatDescr) :: descr, descrLU(2)
#endif
  INTEGER :: nr, nc, nnz, nnzOmp(64), nThread
  LOGICAL :: lAlloc = .FALSE., lFinalized = .FALSE., lDevice = .FALSE., lParallel = .FALSE.
#ifdef __INTEL_MKL
  !--- PARDISO Variables
  TYPE(MKL_PARDISO_HANDLE) :: pardisoPtr(64)
  INTEGER :: pardisoParam(64)
  INTEGER, POINTER :: pardisoPermute(:)
#endif
END TYPE

TYPE CSR_MIXED
  REAL(8), POINTER :: csrVal(:)
  INTEGER, POINTER :: csrRowPtr(:), csrColIdx(:)
  !--- Parallel Generation Variables
  REAL(8), POINTER :: bufVal(:, :)
  INTEGER, POINTER :: bufColIdx(:, :)
#ifdef __PGI
  REAL(4), POINTER, DEVICE :: d_csrVal(:)
  REAL(8), POINTER, DEVICE :: d_csrVal8(:)
  INTEGER, POINTER, DEVICE :: d_csrRowPtr(:), d_csrColIdx(:)
  TYPE(cusparseSolveAnalysisInfo) :: info, infoLU(2)
  TYPE(cusparseMatDescr) :: descr, descrLU(2)
#endif
  INTEGER :: nr, nc, nnz, nnzOmp(64), nThread
  LOGICAL :: lAlloc = .FALSE., lFinalized = .FALSE., lDevice = .FALSE., lParallel = .FALSE.
#ifdef __INTEL_MKL
  !--- PARDISO Variables
  TYPE(MKL_PARDISO_HANDLE) :: pardisoPtr(64)
  INTEGER :: pardisoParam(64)
  INTEGER, POINTER :: pardisoPermute(:)
#endif
END TYPE

INTERFACE createCsr
  MODULE PROCEDURE createCsrFloat
  MODULE PROCEDURE createCsrDouble
  MODULE PROCEDURE createCsrMixed
  MODULE PROCEDURE createCsrFloatParallel
  MODULE PROCEDURE createCsrDoubleParallel
  MODULE PROCEDURE createCsrMixedParallel
END INTERFACE

INTERFACE clearCsr
  MODULE PROCEDURE clearCsrFloat
  MODULE PROCEDURE clearCsrDouble
  MODULE PROCEDURE clearCsrMixed
END INTERFACE

INTERFACE destroyCsr
  MODULE PROCEDURE destroyCsrFloat
  MODULE PROCEDURE destroyCsrDouble
  MODULE PROCEDURE destroyCsrMixed
END INTERFACE

INTERFACE pushCsr
  MODULE PROCEDURE pushCsrFloat4
  MODULE PROCEDURE pushCsrFloat8
  MODULE PROCEDURE pushCsrDouble4
  MODULE PROCEDURE pushCsrDouble8
  MODULE PROCEDURE pushCsrMixed4
  MODULE PROCEDURE pushCsrMixed8
  MODULE PROCEDURE pushCsrParallelFloat4
  MODULE PROCEDURE pushCsrParallelFloat8
  MODULE PROCEDURE pushCsrParallelDouble4
  MODULE PROCEDURE pushCsrParallelDouble8
  MODULE PROCEDURE pushCsrParallelMixed4
  MODULE PROCEDURE pushCsrParallelMixed8
END INTERFACE

#ifdef __INTEL_MKL
INTERFACE transposeCsr
  MODULE PROCEDURE transposeCsrFloat
  MODULE PROCEDURE transposeCsrDouble
  MODULE PROCEDURE transposeCsrMixed
END INTERFACE
#endif

INTERFACE finalizeCsr
  MODULE PROCEDURE finalizeCsrFloat
  MODULE PROCEDURE finalizeCsrDouble
  MODULE PROCEDURE finalizeCsrMixed
END INTERFACE

INTERFACE finalizeSortCsr
  MODULE PROCEDURE finalizeSortCsrFloat
  MODULE PROCEDURE finalizeSortCsrDouble
  MODULE PROCEDURE finalizeSortCsrMixed
END INTERFACE

INTERFACE printCsr
  MODULE PROCEDURE printCsrFloat
  MODULE PROCEDURE printCsrDouble
  MODULE PROCEDURE printCsrMixed
END INTERFACE

INTERFACE ASSIGNMENT (=)
  MODULE PROCEDURE copyCsrFloat
  MODULE PROCEDURE copyCsrDouble
  MODULE PROCEDURE copyCsrMixed
  MODULE PROCEDURE copyCsrFloat2Double
  MODULE PROCEDURE copyCsrFloat2Mixed
  MODULE PROCEDURE copyCsrDouble2Float
  MODULE PROCEDURE copyCsrDouble2Mixed
  MODULE PROCEDURE copyCsrMixed2Float
  MODULE PROCEDURE copyCsrMixed2Double
END INTERFACE

INTERFACE OPERATOR (-)
  MODULE PROCEDURE subCsrCsrFloat
  MODULE PROCEDURE subCsrCsrDouble
END INTERFACE

PRIVATE :: allocCsrParallelBufferFloat, allocCsrParallelBufferDouble, allocCsrParallelBufferMixed
PRIVATE :: collapseCsrParallelBufferFloat, collapseCsrParallelBufferDouble, collapseCsrParallelBufferMixed

CONTAINS

SUBROUTINE createCsrFloat(csrFloat, nnz, nr, nc)

IMPLICIT NONE

TYPE(CSR_FLOAT) :: csrFloat
INTEGER :: nnz, nr, nc
INTEGER :: ierr

IF (csrFloat%lAlloc) CALL destroyCsr(csrFloat)

csrFloat%nnz = 0
csrFloat%nr = nr
csrFloat%nc = nc

ALLOCATE(csrFloat%csrVal(nnz))
ALLOCATE(csrFloat%csrColIdx(nnz))
ALLOCATE(csrFloat%csrRowPtr(nr + 1))
csrFloat%csrVal = 0.0
csrFloat%csrColIdx = 0
csrFloat%csrRowPtr = 0

#ifdef __INTEL_MKL
ALLOCATE(csrFloat%pardisoPermute(nr))
#endif

#ifdef __PGI
ierr = cusparseCreateMatDescr(csrFloat%descr)
ierr = cusparseCreateMatDescr(csrFloat%descrLU(1))
ierr = cusparseCreateMatDescr(csrFloat%descrLU(2))
#endif

csrFloat%lAlloc = .TRUE.

END SUBROUTINE

SUBROUTINE createCsrDouble(csrDouble, nnz, nr, nc)

IMPLICIT NONE

TYPE(CSR_DOUBLE) :: csrDouble
INTEGER :: nnz, nr, nc
INTEGER :: ierr

IF (csrDouble%lAlloc) CALL destroyCsr(csrDouble)

csrDouble%nnz = 0
csrDouble%nr = nr
csrDouble%nc = nc

ALLOCATE(csrDouble%csrVal(nnz))
ALLOCATE(csrDouble%csrColIdx(nnz))
ALLOCATE(csrDouble%csrRowPtr(nr + 1))
csrDouble%csrVal = 0.0
csrDouble%csrColIdx = 0
csrDouble%csrRowPtr = 0

#ifdef __INTEL_MKL
ALLOCATE(csrDouble%pardisoPermute(nr))
#endif

#ifdef __PGI
ierr = cusparseCreateMatDescr(csrDouble%descr)
ierr = cusparseCreateMatDescr(csrDouble%descrLU(1))
ierr = cusparseCreateMatDescr(csrDouble%descrLU(2))
#endif

csrDouble%lAlloc = .TRUE.

END SUBROUTINE

SUBROUTINE createCsrMixed(csrMixed, nnz, nr, nc)

IMPLICIT NONE

TYPE(CSR_MIXED) :: csrMixed
INTEGER :: nnz, nr, nc
INTEGER :: ierr

IF (csrMixed%lAlloc) CALL destroyCsr(csrMixed)

csrMixed%nnz = 0
csrMixed%nr = nr
csrMixed%nc = nc

ALLOCATE(csrMixed%csrVal(nnz))
ALLOCATE(csrMixed%csrColIdx(nnz))
ALLOCATE(csrMixed%csrRowPtr(nr + 1))
csrMixed%csrVal = 0.0
csrMixed%csrColIdx = 0
csrMixed%csrRowPtr = 0

#ifdef __INTEL_MKL
ALLOCATE(csrMixed%pardisoPermute(nr))
#endif

#ifdef __PGI
ierr = cusparseCreateMatDescr(csrMixed%descr)
ierr = cusparseCreateMatDescr(csrMixed%descrLU(1))
ierr = cusparseCreateMatDescr(csrMixed%descrLU(2))
#endif

csrMixed%lAlloc = .TRUE.

END SUBROUTINE

SUBROUTINE createCsrFloatParallel(csrFloat, nnz, nr, nc, nThread)

IMPLICIT NONE

TYPE(CSR_FLOAT) :: csrFloat
INTEGER :: nnz, nr, nc, nThread

CALL createCsr(csrFloat, nnz, nr, nc)
CALL allocCsrParallelBufferFloat(csrFloat, nnz, nThread)

END SUBROUTINE

SUBROUTINE createCsrDoubleParallel(csrDouble, nnz, nr, nc, nThread)

IMPLICIT NONE

TYPE(CSR_DOUBLE) :: csrDouble
INTEGER :: nnz, nr, nc, nThread

CALL createCsr(csrDouble, nnz, nr, nc)
CALL allocCsrParallelBufferDouble(csrDouble, nnz, nThread)

END SUBROUTINE

SUBROUTINE createCsrMixedParallel(csrMixed, nnz, nr, nc, nThread)

IMPLICIT NONE

TYPE(CSR_MIXED) :: csrMixed
INTEGER :: nnz, nr, nc, nThread

CALL createCsr(csrMixed, nnz, nr, nc)
CALL allocCsrParallelBufferMixed(csrMixed, nnz, nThread)

END SUBROUTINE

SUBROUTINE clearCsrFloat(csrFloat)

IMPLICIT NONE

TYPE(CSR_FLOAT) :: csrFloat

csrFloat%csrVal = 0.0
csrFloat%csrRowPtr = 0
csrFloat%csrColIdx = 0

csrFloat%nnz = 0
csrFloat%lFinalized = .FALSE.

END SUBROUTINE

SUBROUTINE clearCsrDouble(csrDouble)

IMPLICIT NONE

TYPE(CSR_DOUBLE) :: csrDouble

csrDouble%csrVal = 0.0
csrDouble%csrRowPtr = 0
csrDouble%csrColIdx = 0

csrDouble%nnz = 0
csrDouble%lFinalized = .FALSE.

END SUBROUTINE

SUBROUTINE clearCsrMixed(csrMixed)

IMPLICIT NONE

TYPE(CSR_MIXED) :: csrMixed

csrMixed%csrVal = 0.0
csrMixed%csrRowPtr = 0
csrMixed%csrColIdx = 0

csrMixed%nnz = 0
csrMixed%lFinalized = .FALSE.

END SUBROUTINE

SUBROUTINE destroyCsrFloat(csrFloat)

IMPLICIT NONE

TYPE(CSR_FLOAT) :: csrFloat
INTEGER :: ierr

IF (.NOT. csrFloat%lAlloc) RETURN

DEALLOCATE(csrFloat%csrVal)
DEALLOCATE(csrFloat%csrColIdx)
DEALLOCATE(csrFloat%csrRowPtr)

#ifdef __PGI
IF (csrFloat%lDevice) THEN
  DEALLOCATE(csrFloat%d_csrVal)
  DEALLOCATE(csrFloat%d_csrColIdx)
  DEALLOCATE(csrFloat%d_csrRowPtr)
ENDIF

ierr = cusparseDestroyMatDescr(csrFloat%descr)
ierr = cusparseDestroyMatDescr(csrFloat%descrLU(1))
ierr = cusparseDestroyMatDescr(csrFloat%descrLU(2))
#endif

#ifdef __INTEL_MKL
DEALLOCATE(csrFloat%pardisoPermute)
#endif

csrFloat%lAlloc = .FALSE.
csrFloat%lFinalized = .FALSE.
csrFloat%lDevice = .FALSE.

END SUBROUTINE

SUBROUTINE destroyCsrDouble(csrDouble)

IMPLICIT NONE

TYPE(CSR_DOUBLE) :: csrDouble
INTEGER :: ierr

IF (.NOT. csrDouble%lAlloc) RETURN

DEALLOCATE(csrDouble%csrVal)
DEALLOCATE(csrDouble%csrColIdx)
DEALLOCATE(csrDouble%csrRowPtr)

#ifdef __PGI
IF (csrDouble%lDevice) THEN
  DEALLOCATE(csrDouble%d_csrVal)
  DEALLOCATE(csrDouble%d_csrColIdx)
  DEALLOCATE(csrDouble%d_csrRowPtr)
ENDIF

ierr = cusparseDestroyMatDescr(csrDouble%descr)
ierr = cusparseDestroyMatDescr(csrDouble%descrLU(1))
ierr = cusparseDestroyMatDescr(csrDouble%descrLU(2))
#endif

#ifdef __INTEL_MKL
DEALLOCATE(csrDouble%pardisoPermute)
#endif

csrDouble%lAlloc = .FALSE.
csrDouble%lFinalized = .FALSE.
csrDouble%lDevice = .FALSE.

END SUBROUTINE

SUBROUTINE destroyCsrMixed(csrMixed)

IMPLICIT NONE

TYPE(CSR_MIXED) :: csrMixed
INTEGER :: ierr

IF (.NOT. csrMixed%lAlloc) RETURN

DEALLOCATE(csrMixed%csrVal)
DEALLOCATE(csrMixed%csrColIdx)
DEALLOCATE(csrMixed%csrRowPtr)

#ifdef __PGI
IF (csrMixed%lDevice) THEN
  DEALLOCATE(csrMixed%d_csrVal)
  DEALLOCATE(csrMixed%d_csrVal8)
  DEALLOCATE(csrMixed%d_csrColIdx)
  DEALLOCATE(csrMixed%d_csrRowPtr)
ENDIF

ierr = cusparseDestroyMatDescr(csrMixed%descr)
ierr = cusparseDestroyMatDescr(csrMixed%descrLU(1))
ierr = cusparseDestroyMatDescr(csrMixed%descrLU(2))
#endif

#ifdef __INTEL_MKL
DEALLOCATE(csrMixed%pardisoPermute)
#endif

csrMixed%lAlloc = .FALSE.
csrMixed%lFinalized = .FALSE.
csrMixed%lDevice = .FALSE.

END SUBROUTINE

SUBROUTINE allocCsrParallelBufferFloat(csrFloat, nnz, nThread)

IMPLICIT NONE

TYPE(CSR_FLOAT) :: csrFloat
INTEGER :: nr, nc, nnz, nThread

IF (csrFloat%lParallel) RETURN

nr = csrFloat%nr
nc = csrFloat%nc

csrFloat%nnzOmp = 0

ALLOCATE(csrFloat%bufVal(buffer * nnz / nThread, nThread))
ALLOCATE(csrFloat%bufColIdx(buffer * nnz / nThread, nThread))
csrFloat%bufVal = 0.0
csrFloat%bufColIdx = 0

csrFloat%nThread = nThread
csrFloat%lParallel = .TRUE.

END SUBROUTINE

SUBROUTINE allocCsrParallelBufferDouble(csrDouble, nnz, nThread)

IMPLICIT NONE

TYPE(CSR_DOUBLE) :: csrDouble
INTEGER :: nr, nc, nnz, nThread

IF (csrDouble%lParallel) RETURN

nr = csrDouble%nr
nc = csrDouble%nc

csrDouble%nnzOmp = 0

ALLOCATE(csrDouble%bufVal(buffer * nnz / nThread, nThread))
ALLOCATE(csrDouble%bufColIdx(buffer * nnz / nThread, nThread))
csrDouble%bufVal = 0.0
csrDouble%bufColIdx = 0

csrDouble%nThread = nThread
csrDouble%lParallel = .TRUE.

END SUBROUTINE

SUBROUTINE allocCsrParallelBufferMixed(csrMixed, nnz, nThread)

IMPLICIT NONE

TYPE(CSR_MIXED) :: csrMixed
INTEGER :: nr, nc, nnz, nThread

IF (csrMixed%lParallel) RETURN

nr = csrMixed%nr
nc = csrMixed%nc

csrMixed%nnzOmp = 0

ALLOCATE(csrMixed%bufVal(buffer * nnz / nThread + 1, nThread))
ALLOCATE(csrMixed%bufColIdx(buffer * nnz / nThread + 1, nThread))
csrMixed%bufVal = 0.0
csrMixed%bufColIdx = 0

csrMixed%nThread = nThread
csrMixed%lParallel = .TRUE.

END SUBROUTINE

SUBROUTINE collapseCsrParallelBufferFloat(csrFloat)

IMPLICIT NONE

TYPE(CSR_FLOAT) :: csrFloat
INTEGER :: nThread
INTEGER :: i, ir, inz, tid

IF (.NOT. csrFloat%lParallel) RETURN

nThread = csrFloat%nThread

inz = 0
DO tid = 1, nThread
  csrFloat%nnz = csrFloat%nnz + csrFloat%nnzOmp(tid)
  DO i = 1, csrFloat%nnzOmp(tid)
    inz = inz + 1
    csrFloat%csrColIdx(inz) = csrFloat%bufColIdx(i, tid)
    csrFloat%csrVal(inz) = csrFloat%bufVal(i, tid)
  ENDDO
ENDDO

DEALLOCATE(csrFloat%bufVal)
DEALLOCATE(csrFloat%bufColIdx)

csrFloat%lParallel = .FALSE.

END SUBROUTINE

SUBROUTINE collapseCsrParallelBufferDouble(csrDouble)

IMPLICIT NONE

TYPE(CSR_DOUBLE) :: csrDouble
INTEGER :: nThread
INTEGER :: i, ir, inz, tid

IF (.NOT. csrDouble%lParallel) RETURN

nThread = csrDouble%nThread

inz = 0
DO tid = 1, nThread
  csrDouble%nnz = csrDouble%nnz + csrDouble%nnzOmp(tid)
  DO i = 1, csrDouble%nnzOmp(tid)
    inz = inz + 1
    csrDouble%csrColIdx(inz) = csrDouble%bufColIdx(i, tid)
    csrDouble%csrVal(inz) = csrDouble%bufVal(i, tid)
  ENDDO
ENDDO

DEALLOCATE(csrDouble%bufVal)
DEALLOCATE(csrDouble%bufColIdx)

csrDouble%lParallel = .FALSE.

END SUBROUTINE

SUBROUTINE collapseCsrParallelBufferMixed(csrMixed)

IMPLICIT NONE

TYPE(CSR_MIXED) :: csrMixed
INTEGER :: nThread
INTEGER :: i, ir, inz, tid

IF (.NOT. csrMixed%lParallel) RETURN

nThread = csrMixed%nThread

inz = 0
DO tid = 1, nThread
  csrMixed%nnz = csrMixed%nnz + csrMixed%nnzOmp(tid)
  DO i = 1, csrMixed%nnzOmp(tid)
    inz = inz + 1
    csrMixed%csrColIdx(inz) = csrMixed%bufColIdx(i, tid)
    csrMixed%csrVal(inz) = csrMixed%bufVal(i, tid)
  ENDDO
ENDDO

DEALLOCATE(csrMixed%bufVal)
DEALLOCATE(csrMixed%bufColIdx)

csrMixed%lParallel = .FALSE.

END SUBROUTINE

SUBROUTINE pushCsrFloat4(csrFloat, val, ir, ic)

IMPLICIT NONE

TYPE(CSR_FLOAT) :: csrFloat
REAL(4) :: val
INTEGER :: ir, ic

IF (abs(val) .LE. 1.0E-20)      RETURN
IF (ir .LE. 0)                  RETURN
IF (ir .GT. csrFloat%nr)        RETURN
IF (ic .LE. 0)                  RETURN
IF (ic .GT. csrFloat%nc)        RETURN

csrFloat%nnz = csrFloat%nnz + 1
csrFloat%csrVal(csrFloat%nnz) = val
csrFloat%csrColIdx(csrFloat%nnz) = ic
csrFloat%csrRowPtr(ir) = csrFloat%csrRowPtr(ir) + 1

END SUBROUTINE

SUBROUTINE pushCsrFloat8(csrFloat, val, ir, ic)

IMPLICIT NONE

TYPE(CSR_FLOAT) :: csrFloat
REAL(8) :: val
INTEGER :: ir, ic

IF (abs(val) .LE. 1.0E-20)      RETURN
IF (ir .LE. 0)                  RETURN
IF (ir .GT. csrFloat%nr)        RETURN
IF (ic .LE. 0)                  RETURN
IF (ic .GT. csrFloat%nc)        RETURN

csrFloat%nnz = csrFloat%nnz + 1
csrFloat%csrVal(csrFloat%nnz) = val
csrFloat%csrColIdx(csrFloat%nnz) = ic
csrFloat%csrRowPtr(ir) = csrFloat%csrRowPtr(ir) + 1

END SUBROUTINE

SUBROUTINE pushCsrDouble4(csrDouble, val, ir, ic)

IMPLICIT NONE

TYPE(CSR_DOUBLE) :: csrDouble
REAL(4) :: val
INTEGER :: ir, ic

IF (abs(val) .LE. 1.0E-20)      RETURN
IF (ir .LE. 0)                  RETURN
IF (ir .GT. csrDouble%nr)       RETURN
IF (ic .LE. 0)                  RETURN
IF (ic .GT. csrDouble%nc)       RETURN

csrDouble%nnz = csrDouble%nnz + 1
csrDouble%csrVal(csrDouble%nnz) = val
csrDouble%csrColIdx(csrDouble%nnz) = ic
csrDouble%csrRowPtr(ir) = csrDouble%csrRowPtr(ir) + 1

END SUBROUTINE

SUBROUTINE pushCsrDouble8(csrDouble, val, ir, ic)

IMPLICIT NONE

TYPE(CSR_DOUBLE) :: csrDouble
REAL(8) :: val
INTEGER :: ir, ic

IF (abs(val) .LE. 1.0E-20)      RETURN
IF (ir .LE. 0)                  RETURN
IF (ir .GT. csrDouble%nr)       RETURN
IF (ic .LE. 0)                  RETURN
IF (ic .GT. csrDouble%nc)       RETURN

csrDouble%nnz = csrDouble%nnz + 1
csrDouble%csrVal(csrDouble%nnz) = val
csrDouble%csrColIdx(csrDouble%nnz) = ic
csrDouble%csrRowPtr(ir) = csrDouble%csrRowPtr(ir) + 1

END SUBROUTINE

SUBROUTINE pushCsrMixed4(csrMixed, val, ir, ic)

IMPLICIT NONE

TYPE(CSR_MIXED) :: csrMixed
REAL(4) :: val
INTEGER :: ir, ic

IF (abs(val) .LE. 1.0E-20)      RETURN
IF (ir .LE. 0)                  RETURN
IF (ir .GT. csrMixed%nr)        RETURN
IF (ic .LE. 0)                  RETURN
IF (ic .GT. csrMixed%nc)        RETURN

csrMixed%nnz = csrMixed%nnz + 1
csrMixed%csrVal(csrMixed%nnz) = val
csrMixed%csrColIdx(csrMixed%nnz) = ic
csrMixed%csrRowPtr(ir) = csrMixed%csrRowPtr(ir) + 1

END SUBROUTINE

SUBROUTINE pushCsrMixed8(csrMixed, val, ir, ic)

IMPLICIT NONE

TYPE(CSR_MIXED) :: csrMixed
REAL(8) :: val
INTEGER :: ir, ic

IF (abs(val) .LE. 1.0E-20)      RETURN
IF (ir .LE. 0)                  RETURN
IF (ir .GT. csrMixed%nr)        RETURN
IF (ic .LE. 0)                  RETURN
IF (ic .GT. csrMixed%nc)        RETURN

csrMixed%nnz = csrMixed%nnz + 1
csrMixed%csrVal(csrMixed%nnz) = val
csrMixed%csrColIdx(csrMixed%nnz) = ic
csrMixed%csrRowPtr(ir) = csrMixed%csrRowPtr(ir) + 1

END SUBROUTINE

SUBROUTINE pushCsrParallelFloat4(csrFloat, val, ir, ic, thread)

IMPLICIT NONE

TYPE(CSR_FLOAT) :: csrFloat
REAL(4) :: val
INTEGER :: ir, ic
INTEGER :: thread

IF (abs(val) .LE. 1.0E-20)      RETURN
IF (ir .LE. 0)                  RETURN
IF (ir .GT. csrFloat%nr)        RETURN
IF (ic .LE. 0)                  RETURN
IF (ic .GT. csrFloat%nc)        RETURN

csrFloat%nnzOmp(thread) = csrFloat%nnzOmp(thread) + 1
csrFloat%bufVal(csrFloat%nnzOmp(thread), thread) = val
csrFloat%bufColIdx(csrFloat%nnzOmp(thread), thread) = ic
csrFloat%csrRowPtr(ir) = csrFloat%csrRowPtr(ir) + 1

END SUBROUTINE

SUBROUTINE pushCsrParallelFloat8(csrFloat, val, ir, ic, thread)

IMPLICIT NONE

TYPE(CSR_FLOAT) :: csrFloat
REAL(8) :: val
INTEGER :: ir, ic
INTEGER :: thread

IF (abs(val) .LE. 1.0E-20)      RETURN
IF (ir .LE. 0)                  RETURN
IF (ir .GT. csrFloat%nr)        RETURN
IF (ic .LE. 0)                  RETURN
IF (ic .GT. csrFloat%nc)        RETURN

csrFloat%nnzOmp(thread) = csrFloat%nnzOmp(thread) + 1
csrFloat%bufVal(csrFloat%nnzOmp(thread), thread) = val
csrFloat%bufColIdx(csrFloat%nnzOmp(thread), thread) = ic
csrFloat%csrRowPtr(ir) = csrFloat%csrRowPtr(ir) + 1

END SUBROUTINE

SUBROUTINE pushCsrParallelDouble4(csrDouble, val, ir, ic, thread)

IMPLICIT NONE

TYPE(CSR_DOUBLE) :: csrDouble
REAL(4) :: val
INTEGER :: ir, ic
INTEGER :: thread

IF (abs(val) .LE. 1.0E-20)      RETURN
IF (ir .LE. 0)                  RETURN
IF (ir .GT. csrDouble%nr)       RETURN
IF (ic .LE. 0)                  RETURN
IF (ic .GT. csrDouble%nc)       RETURN

csrDouble%nnzOmp(thread) = csrDouble%nnzOmp(thread) + 1
csrDouble%bufVal(csrDouble%nnzOmp(thread), thread) = val
csrDouble%bufColIdx(csrDouble%nnzOmp(thread), thread) = ic
csrDouble%csrRowPtr(ir) = csrDouble%csrRowPtr(ir) + 1

END SUBROUTINE

SUBROUTINE pushCsrParallelDouble8(csrDouble, val, ir, ic, thread)

IMPLICIT NONE

TYPE(CSR_DOUBLE) :: csrDouble
REAL(8) :: val
INTEGER :: ir, ic
INTEGER :: thread

IF (abs(val) .LE. 1.0E-20)      RETURN
IF (ir .LE. 0)                  RETURN
IF (ir .GT. csrDouble%nr)       RETURN
IF (ic .LE. 0)                  RETURN
IF (ic .GT. csrDouble%nc)       RETURN

csrDouble%nnzOmp(thread) = csrDouble%nnzOmp(thread) + 1
csrDouble%bufVal(csrDouble%nnzOmp(thread), thread) = val
csrDouble%bufColIdx(csrDouble%nnzOmp(thread), thread) = ic
csrDouble%csrRowPtr(ir) = csrDouble%csrRowPtr(ir) + 1

END SUBROUTINE

SUBROUTINE pushCsrParallelMixed4(csrMixed, val, ir, ic, thread)

IMPLICIT NONE

TYPE(CSR_MIXED) :: csrMixed
REAL(4) :: val
INTEGER :: ir, ic
INTEGER :: thread

IF (abs(val) .LE. 1.0E-20)      RETURN
IF (ir .LE. 0)                  RETURN
IF (ir .GT. csrMixed%nr)        RETURN
IF (ic .LE. 0)                  RETURN
IF (ic .GT. csrMixed%nc)        RETURN

csrMixed%nnzOmp(thread) = csrMixed%nnzOmp(thread) + 1
csrMixed%bufVal(csrMixed%nnzOmp(thread), thread) = val
csrMixed%bufColIdx(csrMixed%nnzOmp(thread), thread) = ic
csrMixed%csrRowPtr(ir) = csrMixed%csrRowPtr(ir) + 1

END SUBROUTINE

SUBROUTINE pushCsrParallelMixed8(csrMixed, val, ir, ic, thread)

IMPLICIT NONE

TYPE(CSR_MIXED) :: csrMixed
REAL(8) :: val
INTEGER :: ir, ic
INTEGER :: thread

IF (abs(val) .LE. 1.0E-20)      RETURN
IF (ir .LE. 0)                  RETURN
IF (ir .GT. csrMixed%nr)        RETURN
IF (ic .LE. 0)                  RETURN
IF (ic .GT. csrMixed%nc)        RETURN

csrMixed%nnzOmp(thread) = csrMixed%nnzOmp(thread) + 1
csrMixed%bufVal(csrMixed%nnzOmp(thread), thread) = val
csrMixed%bufColIdx(csrMixed%nnzOmp(thread), thread) = ic
csrMixed%csrRowPtr(ir) = csrMixed%csrRowPtr(ir) + 1

END SUBROUTINE

#ifdef __INTEL_MKL

SUBROUTINE transposeCsrFloat(inputCsr, outputCsr)

IMPLICIT NONE

TYPE(CSR_FLOAT) :: inputCsr, outputCsr
INTEGER :: nr, nc, nnz
INTEGER :: job(6) = (/ 0, 1, 1, 0, 0, 1 /)
INTEGER :: ierr

nr = inputCsr%nr
nc = inputCsr%nc
nnz = inputCsr%nnz

CALL createCsr(outputCsr, nnz, nr, nc)

CALL mkl_scsrcsc(job, nr, inputCsr%csrVal, inputCsr%csrColIdx, inputCsr%csrRowPtr,                                  &
                 outputCsr%csrVal, outputCsr%csrColIdx, outputCsr%csrRowPtr, ierr)

outputCsr%nnz = nnz
outputCsr%lFinalized = .TRUE.

END SUBROUTINE

SUBROUTINE transposeCsrDouble(inputCsr, outputCsr)

IMPLICIT NONE

TYPE(CSR_DOUBLE) :: inputCsr, outputCsr
INTEGER :: nr, nc, nnz
INTEGER :: job(6) = (/ 0, 1, 1, 0, 0, 1 /)
INTEGER :: ierr

nr = inputCsr%nr
nc = inputCsr%nc
nnz = inputCsr%nnz

CALL createCsr(outputCsr, nnz, nr, nc)

CALL mkl_dcsrcsc(job, nr, inputCsr%csrVal, inputCsr%csrColIdx, inputCsr%csrRowPtr,                                  &
                 outputCsr%csrVal, outputCsr%csrColIdx, outputCsr%csrRowPtr, ierr)

outputCsr%nnz = nnz
outputCsr%lFinalized = .TRUE.

END SUBROUTINE

SUBROUTINE transposeCsrMixed(inputCsr, outputCsr)

IMPLICIT NONE

TYPE(CSR_MIXED) :: inputCsr, outputCsr
INTEGER :: nr, nc, nnz
INTEGER :: job(6) = (/ 0, 1, 1, 0, 0, 1 /)
INTEGER :: ierr

nr = inputCsr%nr
nc = inputCsr%nc
nnz = inputCsr%nnz

CALL createCsr(outputCsr, nnz, nr, nc)

CALL mkl_dcsrcsc(job, nr, inputCsr%csrVal, inputCsr%csrColIdx, inputCsr%csrRowPtr,                                  &
                 outputCsr%csrVal, outputCsr%csrColIdx, outputCsr%csrRowPtr, ierr)

outputCsr%nnz = nnz
outputCsr%lFinalized = .TRUE.

END SUBROUTINE

#endif

SUBROUTINE finalizeCsrFloat(csrFloat, lDevice)

IMPLICIT NONE

TYPE(CSR_FLOAT) :: csrFloat
LOGICAL :: lDevice
REAL(4), POINTER :: csrVal(:)
INTEGER, POINTER :: csrColIdx(:)
INTEGER, ALLOCATABLE :: csrRowPtr(:)
INTEGER :: i, ierr

IF (csrFloat%lParallel) CALL collapseCsrParallelBufferFloat(csrFloat)

IF (.NOT. csrFloat%lFinalized) THEN

  csrVal => csrFloat%csrVal
  csrColIdx => csrFloat%csrColIdx

  ALLOCATE(csrFloat%csrVal(csrFloat%nnz))
  ALLOCATE(csrFloat%csrColIdx(csrFloat%nnz))

  DO i = 1, csrFloat%nnz
    csrFloat%csrVal(i) = csrVal(i)
    csrFloat%csrColIdx(i) = csrColIdx(i)
  ENDDO

  DEALLOCATE(csrVal)
  DEALLOCATE(csrColIdx)

  ALLOCATE(csrRowPtr(csrFloat%nr + 1))

  csrRowPtr = csrFloat%csrRowPtr
  csrFloat%csrRowPtr(1) = 1
  DO i = 2, csrFloat%nr + 1
    csrFloat%csrRowPtr(i) = csrFloat%csrRowPtr(i - 1) + csrRowPtr(i - 1)
  ENDDO

  DEALLOCATE(csrRowPtr)

  csrFloat%lFinalized = .TRUE.

ENDIF

#ifdef __PGI
IF (lDevice) THEN
  IF (.NOT. csrFloat%lDevice) THEN
    ALLOCATE(csrFloat%d_csrVal(csrFloat%nnz))
    ALLOCATE(csrFloat%d_csrColIdx(csrFloat%nnz))
    ALLOCATE(csrFloat%d_csrRowPtr(csrFloat%nr + 1))
    ierr = cusparseSetMatIndexBase(csrFloat%descr, CUSPARSE_INDEX_BASE_ONE)
    ierr = cusparseSetMatIndexBase(csrFloat%descrLU(1), CUSPARSE_INDEX_BASE_ONE)
    ierr = cusparseSetMatIndexBase(csrFloat%descrLU(2), CUSPARSE_INDEX_BASE_ONE)
    ierr = cusparseSetMatType(csrFloat%descr, CUSPARSE_MATRIX_TYPE_GENERAL)
    ierr = cusparseSetMatType(csrFloat%descrLU(1), CUSPARSE_MATRIX_TYPE_GENERAL)
    ierr = cusparseSetMatType(csrFloat%descrLU(2), CUSPARSE_MATRIX_TYPE_GENERAL)
    ierr = cusparseSetMatFillMode(csrFloat%descrLU(1), CUSPARSE_FILL_MODE_LOWER)
    ierr = cusparseSetMatFillMode(csrFloat%descrLU(2), CUSPARSE_FILL_MODE_UPPER)
    ierr = cusparseSetMatDiagType(csrFloat%descrLU(1), CUSPARSE_DIAG_TYPE_NON_UNIT)
    ierr = cusparseSetMatDiagType(csrFloat%descrLU(2), CUSPARSE_DIAG_TYPE_NON_UNIT)
    csrFloat%lDevice = .TRUE.
  ENDIF
  csrFloat%d_csrVal = csrFloat%csrVal
  csrFloat%d_csrColIdx = csrFloat%csrColIdx
  csrFloat%d_csrRowPtr = csrFloat%csrRowPtr
ENDIF
#endif

END SUBROUTINE

SUBROUTINE finalizeCsrDouble(csrDouble, lDevice)

IMPLICIT NONE

TYPE(CSR_DOUBLE) :: csrDouble
LOGICAL :: lDevice
REAL(8), POINTER :: csrVal(:)
INTEGER, POINTER :: csrColIdx(:)
INTEGER, ALLOCATABLE :: csrRowPtr(:)
INTEGER :: i, ierr

IF (csrDouble%lParallel) CALL collapseCsrParallelBufferDouble(csrDouble)

IF (.NOT. csrDouble%lFinalized) THEN

  csrVal => csrDouble%csrVal
  csrColIdx => csrDouble%csrColIdx

  ALLOCATE(csrDouble%csrVal(csrDouble%nnz))
  ALLOCATE(csrDouble%csrColIdx(csrDouble%nnz))

  DO i = 1, csrDouble%nnz
    csrDouble%csrVal(i) = csrVal(i)
    csrDouble%csrColIdx(i) = csrColIdx(i)
  ENDDO

  DEALLOCATE(csrVal)
  DEALLOCATE(csrColIdx)

  ALLOCATE(csrRowPtr(csrDouble%nr + 1))

  csrRowPtr = csrDouble%csrRowPtr
  csrDouble%csrRowPtr(1) = 1
  DO i = 2, csrDouble%nr + 1
    csrDouble%csrRowPtr(i) = csrDouble%csrRowPtr(i - 1) + csrRowPtr(i - 1)
  ENDDO

  DEALLOCATE(csrRowPtr)

  csrDouble%lFinalized = .TRUE.

ENDIF

#ifdef __PGI
IF (lDevice) THEN
  IF (.NOT. csrDouble%lDevice) THEN
    ALLOCATE(csrDouble%d_csrVal(csrDouble%nnz))
    ALLOCATE(csrDouble%d_csrColIdx(csrDouble%nnz))
    ALLOCATE(csrDouble%d_csrRowPtr(csrDouble%nr + 1))
    ierr = cusparseSetMatIndexBase(csrDouble%descr, CUSPARSE_INDEX_BASE_ONE)
    ierr = cusparseSetMatIndexBase(csrDouble%descrLU(1), CUSPARSE_INDEX_BASE_ONE)
    ierr = cusparseSetMatIndexBase(csrDouble%descrLU(2), CUSPARSE_INDEX_BASE_ONE)
    ierr = cusparseSetMatType(csrDouble%descr, CUSPARSE_MATRIX_TYPE_GENERAL)
    ierr = cusparseSetMatType(csrDouble%descrLU(1), CUSPARSE_MATRIX_TYPE_GENERAL)
    ierr = cusparseSetMatType(csrDouble%descrLU(2), CUSPARSE_MATRIX_TYPE_GENERAL)
    ierr = cusparseSetMatFillMode(csrDouble%descrLU(1), CUSPARSE_FILL_MODE_LOWER)
    ierr = cusparseSetMatFillMode(csrDouble%descrLU(2), CUSPARSE_FILL_MODE_UPPER)
    ierr = cusparseSetMatDiagType(csrDouble%descrLU(1), CUSPARSE_DIAG_TYPE_NON_UNIT)
    ierr = cusparseSetMatDiagType(csrDouble%descrLU(2), CUSPARSE_DIAG_TYPE_NON_UNIT)
    csrDouble%lDevice = .TRUE.
  ENDIF
  csrDouble%d_csrVal = csrDouble%csrVal
  csrDouble%d_csrColIdx = csrDouble%csrColIdx
  csrDouble%d_csrRowPtr = csrDouble%csrRowPtr
ENDIF
#endif

END SUBROUTINE

SUBROUTINE finalizeCsrMixed(csrMixed, lDevice)

IMPLICIT NONE

TYPE(CSR_MIXED) :: csrMixed
LOGICAL :: lDevice
REAL(8), POINTER :: csrVal(:)
INTEGER, POINTER :: csrColIdx(:)
INTEGER, ALLOCATABLE :: csrRowPtr(:)
INTEGER :: i, ierr

IF (csrMixed%lParallel) CALL collapseCsrParallelBufferMixed(csrMixed)

IF (.NOT. csrMixed%lFinalized) THEN

  csrVal => csrMixed%csrVal
  csrColIdx => csrMixed%csrColIdx

  ALLOCATE(csrMixed%csrVal(csrMixed%nnz))
  ALLOCATE(csrMixed%csrColIdx(csrMixed%nnz))

  DO i = 1, csrMixed%nnz
    csrMixed%csrVal(i) = csrVal(i)
    csrMixed%csrColIdx(i) = csrColIdx(i)
  ENDDO

  DEALLOCATE(csrVal)
  DEALLOCATE(csrColIdx)

  ALLOCATE(csrRowPtr(csrMixed%nr + 1))

  csrRowPtr = csrMixed%csrRowPtr
  csrMixed%csrRowPtr(1) = 1
  DO i = 2, csrMixed%nr + 1
    csrMixed%csrRowPtr(i) = csrMixed%csrRowPtr(i - 1) + csrRowPtr(i - 1)
  ENDDO

  DEALLOCATE(csrRowPtr)

  csrMixed%lFinalized = .TRUE.

ENDIF

#ifdef __PGI
IF (lDevice) THEN
  IF (.NOT. csrMixed%lDevice) THEN
    ALLOCATE(csrMixed%d_csrVal(csrMixed%nnz))
    ALLOCATE(csrMixed%d_csrVal8(csrMixed%nnz))
    ALLOCATE(csrMixed%d_csrColIdx(csrMixed%nnz))
    ALLOCATE(csrMixed%d_csrRowPtr(csrMixed%nr + 1))
    ierr = cusparseSetMatIndexBase(csrMixed%descr, CUSPARSE_INDEX_BASE_ONE)
    ierr = cusparseSetMatIndexBase(csrMixed%descrLU(1), CUSPARSE_INDEX_BASE_ONE)
    ierr = cusparseSetMatIndexBase(csrMixed%descrLU(2), CUSPARSE_INDEX_BASE_ONE)
    ierr = cusparseSetMatType(csrMixed%descr, CUSPARSE_MATRIX_TYPE_GENERAL)
    ierr = cusparseSetMatType(csrMixed%descrLU(1), CUSPARSE_MATRIX_TYPE_GENERAL)
    ierr = cusparseSetMatType(csrMixed%descrLU(2), CUSPARSE_MATRIX_TYPE_GENERAL)
    ierr = cusparseSetMatFillMode(csrMixed%descrLU(1), CUSPARSE_FILL_MODE_LOWER)
    ierr = cusparseSetMatFillMode(csrMixed%descrLU(2), CUSPARSE_FILL_MODE_UPPER)
    ierr = cusparseSetMatDiagType(csrMixed%descrLU(1), CUSPARSE_DIAG_TYPE_NON_UNIT)
    ierr = cusparseSetMatDiagType(csrMixed%descrLU(2), CUSPARSE_DIAG_TYPE_NON_UNIT)
    csrMixed%lDevice = .TRUE.
  ENDIF
  csrMixed%d_csrVal = csrMixed%csrVal
  csrMixed%d_csrVal8 = csrMixed%csrVal
  csrMixed%d_csrColIdx = csrMixed%csrColIdx
  csrMixed%d_csrRowPtr = csrMixed%csrRowPtr
ENDIF
#endif

END SUBROUTINE

SUBROUTINE finalizeSortCsrFloat(csrFloat, lDevice)

IMPLICIT NONE

TYPE(CSR_FLOAT) :: csrFloat
LOGICAL :: lDevice
REAL(4), ALLOCATABLE :: rowVal(:, :)
REAL(4) :: val
INTEGER, ALLOCATABLE :: rowCol(:, :)
INTEGER :: maxRowEntry
INTEGER :: idx, ic, ir, i, j

IF (csrFloat%lFinalized) RETURN
IF (csrFloat%lParallel) CALL collapseCsrParallelBufferFloat(csrFloat)

maxRowEntry = maxval(csrFloat%csrRowPtr)

ALLOCATE(rowVal(maxRowEntry, csrFloat%nr))
ALLOCATE(rowCol(maxRowEntry, csrFloat%nr))

idx = 0

DO ir = 1, csrFloat%nr
  DO ic = 1, csrFloat%csrRowPtr(ir)
    idx = idx + 1
    rowVal(ic, ir) = csrFloat%csrVal(idx)
    rowCol(ic, ir) = csrFloat%csrColIdx(idx)
  ENDDO
ENDDO

!$OMP PARALLEL PRIVATE(val, ic)
!$OMP DO SCHEDULE(GUIDED)
DO ir = 1, csrFloat%nr
  DO j = csrFloat%csrRowPtr(ir), 1, -1
    DO i = 1, j - 1
      IF (rowCol(i, ir) .GT. rowCol(i + 1, ir)) THEN
        val = rowVal(i, ir)
        rowVal(i, ir) = rowVal(i + 1, ir); rowVal(i + 1, ir) = val
        ic = rowCol(i, ir)
        rowCol(i, ir) = rowCol(i + 1, ir); rowCol(i + 1, ir) = ic
      ENDIF
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

idx = 0

DO ir = 1, csrFloat%nr
  DO ic = 1, csrFloat%csrRowPtr(ir)
    idx = idx + 1
    csrFloat%csrVal(idx) = rowVal(ic, ir)
    csrFloat%csrColIdx(idx) = rowCol(ic, ir)
  ENDDO
ENDDO

CALL finalizeCsrFloat(csrFloat, lDevice)

DEALLOCATE(rowVal)
DEALLOCATE(rowCol)

END SUBROUTINE

SUBROUTINE finalizeSortCsrDouble(csrDouble, lDevice)

IMPLICIT NONE

TYPE(CSR_DOUBLE) :: csrDouble
LOGICAL :: lDevice
REAL(8), ALLOCATABLE :: rowVal(:, :)
REAL(8) :: val
INTEGER, ALLOCATABLE :: rowCol(:, :)
INTEGER :: maxRowEntry
INTEGER :: idx, ic, ir, i, j

IF (csrDouble%lFinalized) RETURN
IF (csrDouble%lParallel) CALL collapseCsrParallelBufferDouble(csrDouble)

maxRowEntry = maxval(csrDouble%csrRowPtr)

ALLOCATE(rowVal(maxRowEntry, csrDouble%nr))
ALLOCATE(rowCol(maxRowEntry, csrDouble%nr))

idx = 0

DO ir = 1, csrDouble%nr
  DO ic = 1, csrDouble%csrRowPtr(ir)
    idx = idx + 1
    rowVal(ic, ir) = csrDouble%csrVal(idx)
    rowCol(ic, ir) = csrDouble%csrColIdx(idx)
  ENDDO
ENDDO

!$OMP PARALLEL PRIVATE(val, ic)
!$OMP DO SCHEDULE(GUIDED)
DO ir = 1, csrDouble%nr
  DO j = csrDouble%csrRowPtr(ir), 1, -1
    DO i = 1, j - 1
      IF (rowCol(i, ir) .GT. rowCol(i + 1, ir)) THEN
        val = rowVal(i, ir)
        rowVal(i, ir) = rowVal(i + 1, ir); rowVal(i + 1, ir) = val
        ic = rowCol(i, ir)
        rowCol(i, ir) = rowCol(i + 1, ir); rowCol(i + 1, ir) = ic
      ENDIF
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

idx = 0

DO ir = 1, csrDouble%nr
  DO ic = 1, csrDouble%csrRowPtr(ir)
    idx = idx + 1
    csrDouble%csrVal(idx) = rowVal(ic, ir)
    csrDouble%csrColIdx(idx) = rowCol(ic, ir)
  ENDDO
ENDDO

CALL finalizeCsrDouble(csrDouble, lDevice)

DEALLOCATE(rowVal)
DEALLOCATE(rowCol)

END SUBROUTINE

SUBROUTINE finalizeSortCsrMixed(csrMixed, lDevice)

IMPLICIT NONE

TYPE(CSR_MIXED) :: csrMixed
LOGICAL :: lDevice
REAL(8), ALLOCATABLE :: rowVal(:, :)
REAL(8) :: val
INTEGER, ALLOCATABLE :: rowCol(:, :)
INTEGER :: maxRowEntry
INTEGER :: idx, ic, ir, i, j

IF (csrMixed%lFinalized) RETURN
IF (csrMixed%lParallel) CALL collapseCsrParallelBufferMixed(csrMixed)

maxRowEntry = maxval(csrMixed%csrRowPtr)

ALLOCATE(rowVal(maxRowEntry, csrMixed%nr))
ALLOCATE(rowCol(maxRowEntry, csrMixed%nr))

idx = 0

DO ir = 1, csrMixed%nr
  DO ic = 1, csrMixed%csrRowPtr(ir)
    idx = idx + 1
    rowVal(ic, ir) = csrMixed%csrVal(idx)
    rowCol(ic, ir) = csrMixed%csrColIdx(idx)
  ENDDO
ENDDO

!$OMP PARALLEL PRIVATE(val, ic)
!$OMP DO SCHEDULE(GUIDED)
DO ir = 1, csrMixed%nr
  DO j = csrMixed%csrRowPtr(ir), 1, -1
    DO i = 1, j - 1
      IF (rowCol(i, ir) .GT. rowCol(i + 1, ir)) THEN
        val = rowVal(i, ir)
        rowVal(i, ir) = rowVal(i + 1, ir); rowVal(i + 1, ir) = val
        ic = rowCol(i, ir)
        rowCol(i, ir) = rowCol(i + 1, ir); rowCol(i + 1, ir) = ic
      ENDIF
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

idx = 0

DO ir = 1, csrMixed%nr
  DO ic = 1, csrMixed%csrRowPtr(ir)
    idx = idx + 1
    csrMixed%csrVal(idx) = rowVal(ic, ir)
    csrMixed%csrColIdx(idx) = rowCol(ic, ir)
  ENDDO
ENDDO

CALL finalizeCsrMixed(csrMixed, lDevice)

DEALLOCATE(rowVal)
DEALLOCATE(rowCol)

END SUBROUTINE

SUBROUTINE printCsrFloat(csrFloat, filename, io)

IMPLICIT NONE

TYPE(CSR_FLOAT) :: csrFloat
CHARACTER(*) :: filename
INTEGER :: io
INTEGER :: i

OPEN(io, FILE = filename)
WRITE(io, *), csrFloat%nr, csrFloat%nc, csrFloat%nnz
DO i = 1, csrFloat%nnz
  IF (i .LE. csrFloat%nr + 1) THEN
    WRITE(io, *), csrFloat%csrVal(i), csrFloat%csrColIdx(i), csrFloat%csrRowPtr(i)
  ELSE
    WRITE(io, *), csrFloat%csrVal(i), csrFloat%csrColIdx(i), 0
  ENDIF
ENDDO
CLOSE(io)

END SUBROUTINE

SUBROUTINE printCsrDouble(csrDouble, filename, io)

IMPLICIT NONE

TYPE(CSR_DOUBLE) :: csrDouble
CHARACTER(*) :: filename
INTEGER :: io
INTEGER :: i

OPEN(io, FILE = filename)
WRITE(io, *), csrDouble%nr, csrDouble%nc, csrDouble%nnz
DO i = 1, csrDouble%nnz
  IF (i .LE. csrDouble%nr + 1) THEN
    WRITE(io, *), csrDouble%csrVal(i), csrDouble%csrColIdx(i), csrDouble%csrRowPtr(i)
  ELSE
    WRITE(io, *), csrDouble%csrVal(i), csrDouble%csrColIdx(i), 0
  ENDIF
ENDDO
CLOSE(io)

END SUBROUTINE

SUBROUTINE printCsrMixed(csrMixed, filename, io)

IMPLICIT NONE

TYPE(CSR_MIXED) :: csrMixed
CHARACTER(*) :: filename
INTEGER :: io
INTEGER :: i

OPEN(io, FILE = filename)
WRITE(io, *), csrMixed%nr, csrMixed%nc, csrMixed%nnz
DO i = 1, csrMixed%nnz
  IF (i .LE. csrMixed%nr + 1) THEN
    WRITE(io, *), csrMixed%csrVal(i), csrMixed%csrColIdx(i), csrMixed%csrRowPtr(i)
  ELSE
    WRITE(io, *), csrMixed%csrVal(i), csrMixed%csrColIdx(i), 0
  ENDIF
ENDDO
CLOSE(io)

END SUBROUTINE

SUBROUTINE copyCsrFloat(outputCsr, inputCsr)

IMPLICIT NONE

TYPE(CSR_FLOAT), INTENT(INOUT) :: outputCsr
TYPE(CSR_FLOAT), INTENT(IN) :: inputCsr

CALL createCsr(outputCsr, inputCsr%nnz, inputCsr%nr, inputCsr%nc)

outputCsr%csrVal = inputCsr%csrVal
outputCsr%csrRowPtr = inputCsr%csrRowPtr
outputCsr%csrColIdx = inputCsr%csrColIdx
outputCsr%nr = inputCsr%nr
outputCsr%nc = inputCsr%nc
outputCsr%nnz = inputCsr%nnz
outputCsr%lFinalized = inputCsr%lFinalized

#ifdef __INTEL_MKL
outputCsr%pardisoPermute = inputCsr%pardisoPermute
#endif

END SUBROUTINE

SUBROUTINE copyCsrDouble(outputCsr, inputCsr)

IMPLICIT NONE

TYPE(CSR_DOUBLE), INTENT(INOUT) :: outputCsr
TYPE(CSR_DOUBLE), INTENT(IN) :: inputCsr

CALL createCsr(outputCsr, inputCsr%nnz, inputCsr%nr, inputCsr%nc)

outputCsr%csrVal = inputCsr%csrVal
outputCsr%csrRowPtr = inputCsr%csrRowPtr
outputCsr%csrColIdx = inputCsr%csrColIdx
outputCsr%nr = inputCsr%nr
outputCsr%nc = inputCsr%nc
outputCsr%nnz = inputCsr%nnz
outputCsr%lFinalized = inputCsr%lFinalized

#ifdef __INTEL_MKL
outputCsr%pardisoPermute = inputCsr%pardisoPermute
#endif

END SUBROUTINE

SUBROUTINE copyCsrMixed(outputCsr, inputCsr)

IMPLICIT NONE

TYPE(CSR_MIXED), INTENT(INOUT) :: outputCsr
TYPE(CSR_MIXED), INTENT(IN) :: inputCsr

CALL createCsr(outputCsr, inputCsr%nnz, inputCsr%nr, inputCsr%nc)

outputCsr%csrVal = inputCsr%csrVal
outputCsr%csrRowPtr = inputCsr%csrRowPtr
outputCsr%csrColIdx = inputCsr%csrColIdx
outputCsr%nr = inputCsr%nr
outputCsr%nc = inputCsr%nc
outputCsr%nnz = inputCsr%nnz
outputCsr%lFinalized = inputCsr%lFinalized

#ifdef __INTEL_MKL
outputCsr%pardisoPermute = inputCsr%pardisoPermute
#endif

END SUBROUTINE

SUBROUTINE copyCsrFloat2Double(outputCsr, inputCsr)

IMPLICIT NONE

TYPE(CSR_DOUBLE), INTENT(INOUT) :: outputCsr
TYPE(CSR_FLOAT), INTENT(IN) :: inputCsr

CALL createCsr(outputCsr, inputCsr%nnz, inputCsr%nr, inputCsr%nc)

outputCsr%csrVal = inputCsr%csrVal
outputCsr%csrRowPtr = inputCsr%csrRowPtr
outputCsr%csrColIdx = inputCsr%csrColIdx
outputCsr%nr = inputCsr%nr
outputCsr%nc = inputCsr%nc
outputCsr%nnz = inputCsr%nnz
outputCsr%lFinalized = inputCsr%lFinalized

#ifdef __INTEL_MKL
outputCsr%pardisoPermute = inputCsr%pardisoPermute
#endif

END SUBROUTINE

SUBROUTINE copyCsrFloat2Mixed(outputCsr, inputCsr)

IMPLICIT NONE

TYPE(CSR_MIXED), INTENT(INOUT) :: outputCsr
TYPE(CSR_FLOAT), INTENT(IN) :: inputCsr

CALL createCsr(outputCsr, inputCsr%nnz, inputCsr%nr, inputCsr%nc)

outputCsr%csrVal = inputCsr%csrVal
outputCsr%csrRowPtr = inputCsr%csrRowPtr
outputCsr%csrColIdx = inputCsr%csrColIdx
outputCsr%nr = inputCsr%nr
outputCsr%nc = inputCsr%nc
outputCsr%nnz = inputCsr%nnz
outputCsr%lFinalized = inputCsr%lFinalized

#ifdef __INTEL_MKL
outputCsr%pardisoPermute = inputCsr%pardisoPermute
#endif

END SUBROUTINE

SUBROUTINE copyCsrDouble2Float(outputCsr, inputCsr)

IMPLICIT NONE

TYPE(CSR_FLOAT), INTENT(INOUT) :: outputCsr
TYPE(CSR_DOUBLE), INTENT(IN) :: inputCsr

CALL createCsr(outputCsr, inputCsr%nnz, inputCsr%nr, inputCsr%nc)

outputCsr%csrVal = inputCsr%csrVal
outputCsr%csrRowPtr = inputCsr%csrRowPtr
outputCsr%csrColIdx = inputCsr%csrColIdx
outputCsr%nr = inputCsr%nr
outputCsr%nc = inputCsr%nc
outputCsr%nnz = inputCsr%nnz
outputCsr%lFinalized = inputCsr%lFinalized

#ifdef __INTEL_MKL
outputCsr%pardisoPermute = inputCsr%pardisoPermute
#endif

END SUBROUTINE

SUBROUTINE copyCsrDouble2Mixed(outputCsr, inputCsr)

IMPLICIT NONE

TYPE(CSR_MIXED), INTENT(INOUT) :: outputCsr
TYPE(CSR_DOUBLE), INTENT(IN) :: inputCsr

CALL createCsr(outputCsr, inputCsr%nnz, inputCsr%nr, inputCsr%nc)

outputCsr%csrVal = inputCsr%csrVal
outputCsr%csrRowPtr = inputCsr%csrRowPtr
outputCsr%csrColIdx = inputCsr%csrColIdx
outputCsr%nr = inputCsr%nr
outputCsr%nc = inputCsr%nc
outputCsr%nnz = inputCsr%nnz
outputCsr%lFinalized = inputCsr%lFinalized

#ifdef __INTEL_MKL
outputCsr%pardisoPermute = inputCsr%pardisoPermute
#endif

END SUBROUTINE

SUBROUTINE copyCsrMixed2Float(outputCsr, inputCsr)

IMPLICIT NONE

TYPE(CSR_FLOAT), INTENT(INOUT) :: outputCsr
TYPE(CSR_MIXED), INTENT(IN) :: inputCsr

CALL createCsr(outputCsr, inputCsr%nnz, inputCsr%nr, inputCsr%nc)

outputCsr%csrVal = inputCsr%csrVal
outputCsr%csrRowPtr = inputCsr%csrRowPtr
outputCsr%csrColIdx = inputCsr%csrColIdx
outputCsr%nr = inputCsr%nr
outputCsr%nc = inputCsr%nc
outputCsr%nnz = inputCsr%nnz
outputCsr%lFinalized = inputCsr%lFinalized

#ifdef __INTEL_MKL
outputCsr%pardisoPermute = inputCsr%pardisoPermute
#endif

END SUBROUTINE

SUBROUTINE copyCsrMixed2Double(outputCsr, inputCsr)

IMPLICIT NONE

TYPE(CSR_DOUBLE), INTENT(INOUT) :: outputCsr
TYPE(CSR_MIXED), INTENT(IN) :: inputCsr

CALL createCsr(outputCsr, inputCsr%nnz, inputCsr%nr, inputCsr%nc)

outputCsr%csrVal = inputCsr%csrVal
outputCsr%csrRowPtr = inputCsr%csrRowPtr
outputCsr%csrColIdx = inputCsr%csrColIdx
outputCsr%nr = inputCsr%nr
outputCsr%nc = inputCsr%nc
outputCsr%nnz = inputCsr%nnz
outputCsr%lFinalized = inputCsr%lFinalized

#ifdef __INTEL_MKL
outputCsr%pardisoPermute = inputCsr%pardisoPermute
#endif

END SUBROUTINE

FUNCTION subCsrCsrFloat(leftCsr, rightCsr) RESULT(resultCsr)

IMPLICIT NONE

TYPE(CSR_FLOAT), INTENT(IN) :: leftCsr, rightCsr
TYPE(CSR_FLOAT) :: resultCsr
REAL(4) :: diff
INTEGER :: i, left_ptr, right_ptr

CALL createCsr(resultCsr, leftCsr%nnz + rightCsr%nnz, leftCsr%nr, leftCsr%nc)

left_ptr = 1; right_ptr = 1
DO i = 1, leftCsr%nr

  DO WHILE (left_ptr .LT. leftCsr%csrRowPtr(i + 1) .AND. right_ptr .LT. rightCsr%csrRowPtr(i + 1))
    IF (leftCsr%csrColIdx(left_ptr) .LT. rightCsr%csrColIdx(right_ptr)) THEN
      CALL pushCsr(resultCsr, leftCsr%csrVal(left_ptr), i, leftCsr%csrColIdx(left_ptr))
      left_ptr = left_ptr + 1
    ELSEIF (leftCsr%csrColIdx(left_ptr) .EQ. rightCsr%csrColIdx(right_ptr)) THEN
      diff = leftCsr%csrVal(left_ptr) - rightCsr%csrVal(right_ptr)
      IF (abs(diff) .GT. 1.0D-10) CALL pushCsr(resultCsr, diff, i, leftCsr%csrColIdx(left_ptr))
      left_ptr = left_ptr + 1
      right_ptr = right_ptr + 1
    ELSE
      CALL pushCsr(resultCsr, - rightCsr%csrVal(right_ptr), i, rightCsr%csrColIdx(right_ptr))
      right_ptr = right_ptr + 1
    ENDIF
  ENDDO

  DO WHILE (left_ptr .LT. leftCsr%csrRowPtr(i + 1))
    CALL pushCsr(resultCsr, leftCsr%csrVal(left_ptr), i, leftCsr%csrColIdx(left_ptr))
    left_ptr = left_ptr + 1
  ENDDO

  DO WHILE (right_ptr .LT. rightCsr%csrRowPtr(i + 1))
    CALL pushCsr(resultCsr, - rightCsr%csrVal(right_ptr), i, rightCsr%csrColIdx(right_ptr))
    right_ptr = right_ptr + 1
  ENDDO

ENDDO

CALL finalizeCsr(resultCsr, .FALSE.)

END FUNCTION

FUNCTION subCsrCsrDouble(leftCsr, rightCsr) RESULT(resultCsr)

IMPLICIT NONE

TYPE(CSR_DOUBLE), INTENT(IN) :: leftCsr, rightCsr
TYPE(CSR_DOUBLE) :: resultCsr
REAL(8) :: diff
INTEGER :: i, left_ptr, right_ptr

CALL createCsr(resultCsr, leftCsr%nnz + rightCsr%nnz, leftCsr%nr, leftCsr%nc)

left_ptr = 1; right_ptr = 1
DO i = 1, leftCsr%nr

  DO WHILE (left_ptr .LT. leftCsr%csrRowPtr(i + 1) .AND. right_ptr .LT. rightCsr%csrRowPtr(i + 1))
    IF (leftCsr%csrColIdx(left_ptr) .LT. rightCsr%csrColIdx(right_ptr)) THEN
      CALL pushCsr(resultCsr, leftCsr%csrVal(left_ptr), i, leftCsr%csrColIdx(left_ptr))
      left_ptr = left_ptr + 1
    ELSEIF (leftCsr%csrColIdx(left_ptr) .EQ. rightCsr%csrColIdx(right_ptr)) THEN
      diff = leftCsr%csrVal(left_ptr) - rightCsr%csrVal(right_ptr)
      IF (abs(diff) .GT. 1.0D-10) CALL pushCsr(resultCsr, diff, i, leftCsr%csrColIdx(left_ptr))
      left_ptr = left_ptr + 1
      right_ptr = right_ptr + 1
    ELSE
      CALL pushCsr(resultCsr, - rightCsr%csrVal(right_ptr), i, rightCsr%csrColIdx(right_ptr))
      right_ptr = right_ptr + 1
    ENDIF
  ENDDO

  DO WHILE (left_ptr .LT. leftCsr%csrRowPtr(i + 1))
    CALL pushCsr(resultCsr, leftCsr%csrVal(left_ptr), i, leftCsr%csrColIdx(left_ptr))
    left_ptr = left_ptr + 1
  ENDDO

  DO WHILE (right_ptr .LT. rightCsr%csrRowPtr(i + 1))
    CALL pushCsr(resultCsr, - rightCsr%csrVal(right_ptr), i, rightCsr%csrColIdx(right_ptr))
    right_ptr = right_ptr + 1
  ENDDO

ENDDO

CALL finalizeCsr(resultCsr, .FALSE.)

END FUNCTION

END MODULE

MODULE BSRMATRIX

#ifdef __PGI
USE CUSPARSE
#endif

IMPLICIT NONE

TYPE BSR_FLOAT
  REAL(4), POINTER :: bsrVal(:, :, :)
  INTEGER, POINTER :: bsrRowPtr(:), bsrColIdx(:)
#ifdef __PGI
  REAL(4), POINTER, DEVICE :: d_bsrVal(:, :, :)
  INTEGER, POINTER, DEVICE :: d_bsrRowPtr(:), d_bsrColIdx(:)
  TYPE(cusparseBsrsv2Info) :: info, infoLU(2)
  TYPE(cusparseMatDescr) :: descr, descrLU(2)
  CHARACTER(c_char), POINTER, DEVICE :: pBuffer(:)
#endif
  INTEGER :: nbr, nbc, nBlock, blockSize
  LOGICAL :: lAlloc = .FALSE., lDevice = .FALSE.
END TYPE

TYPE BSR_DOUBLE
  REAL(8), POINTER :: bsrVal(:, :, :)
  INTEGER, POINTER :: bsrRowPtr(:), bsrColIdx(:)
#ifdef __PGI
  REAL(8), POINTER, DEVICE :: d_bsrVal(:, :, :)
  INTEGER, POINTER, DEVICE :: d_bsrRowPtr(:), d_bsrColIdx(:)
  TYPE(cusparseBsrsv2Info) :: info, infoLU(2)
  TYPE(cusparseMatDescr) :: descr, descrLU(2)
  CHARACTER(c_char), POINTER, DEVICE :: pBuffer(:)
#endif
  INTEGER :: nbr, nbc, nBlock, blockSize
  LOGICAL :: lAlloc = .FALSE., lDevice = .FALSE.
END TYPE

TYPE BSR_MIXED
  REAL(8), POINTER :: bsrVal(:, :, :)
  INTEGER, POINTER :: bsrRowPtr(:), bsrColIdx(:)
#ifdef __PGI
  REAL(4), POINTER, DEVICE :: d_bsrVal(:, :, :)
  REAL(8), POINTER, DEVICE :: d_bsrVal8(:, :, :)
  INTEGER, POINTER, DEVICE :: d_bsrRowPtr(:), d_bsrColIdx(:)
  TYPE(cusparseBsrsv2Info) :: info, infoLU(2)
  TYPE(cusparseMatDescr) :: descr, descrLU(2)
  CHARACTER(c_char), POINTER, DEVICE :: pBuffer(:)
#endif
  INTEGER :: nbr, nbc, nBlock, blockSize
  LOGICAL :: lAlloc = .FALSE., lDevice = .FALSE.
END TYPE

INTERFACE createBsr
  MODULE PROCEDURE createBsrFloat
  MODULE PROCEDURE createBsrDouble
  MODULE PROCEDURE createBsrMixed
END INTERFACE

INTERFACE clearBsr
  MODULE PROCEDURE clearBsrFloat
  MODULE PROCEDURE clearBsrDouble
  MODULE PROCEDURE clearBsrMixed
END INTERFACE

INTERFACE destroyBsr
  MODULE PROCEDURE destroyBsrFloat
  MODULE PROCEDURE destroyBsrDouble
  MODULE PROCEDURE destroyBsrMixed
END INTERFACE

#ifdef __PGI
INTERFACE copyDeviceBsr
  MODULE PROCEDURE copyDeviceBsrFloat
  MODULE PROCEDURE copyDeviceBsrDouble
  MODULE PROCEDURE copyDeviceBsrMixed
END INTERFACE
#endif

! INTERFACE ASSIGNMENT (=)
!   MODULE PROCEDURE copyBsrFloat
!   MODULE PROCEDURE copyBsrDouble
!   MODULE PROCEDURE copyBsrMixed
!   MODULE PROCEDURE copyBsrFloat2Double
!   MODULE PROCEDURE copyBsrFloat2Mixed
!   MODULE PROCEDURE copyBsrDouble2Float
!   MODULE PROCEDURE copyBsrDouble2Mixed
!   MODULE PROCEDURE copyBsrMixed2Float
!   MODULE PROCEDURE copyBsrMixed2Double
! END INTERFACE

CONTAINS

SUBROUTINE createBsrFloat(bsrFloat, nbr, nbc, nBlock, blockSize)

IMPLICIT NONE

TYPE(BSR_FLOAT) :: bsrFloat
INTEGER :: nbr, nbc, nBlock, blockSize
INTEGER :: ierr

IF (bsrFloat%lAlloc) CALL destroyBsr(bsrFloat)

bsrFloat%nbr = nbr
bsrFloat%nbc = nbc
bsrFloat%nBlock = nBlock
bsrFloat%blockSize = blockSize

ALLOCATE(bsrFloat%bsrVal(blockSize, blockSize, nBlock))
ALLOCATE(bsrFloat%bsrColIdx(nBlock))
ALLOCATE(bsrFloat%bsrRowPtr(nbr + 1))
bsrFloat%bsrVal = 0.0
bsrFloat%bsrColIdx = 0
bsrFloat%bsrRowPtr = 0

#ifdef __PGI
ierr = cusparseCreateBsrsv2Info(bsrFloat%info)
ierr = cusparseCreateBsrsv2Info(bsrFloat%infoLU(1))
ierr = cusparseCreateBsrsv2Info(bsrFloat%infoLU(2))
ierr = cusparseCreateMatDescr(bsrFloat%descr)
ierr = cusparseCreateMatDescr(bsrFloat%descrLU(1))
ierr = cusparseCreateMatDescr(bsrFloat%descrLU(2))
ierr = cusparseSetMatIndexBase(bsrFloat%descr, CUSPARSE_INDEX_BASE_ONE)
ierr = cusparseSetMatIndexBase(bsrFloat%descrLU(1), CUSPARSE_INDEX_BASE_ONE)
ierr = cusparseSetMatIndexBase(bsrFloat%descrLU(2), CUSPARSE_INDEX_BASE_ONE)
ierr = cusparseSetMatType(bsrFloat%descr, CUSPARSE_MATRIX_TYPE_GENERAL)
ierr = cusparseSetMatType(bsrFloat%descrLU(1), CUSPARSE_MATRIX_TYPE_GENERAL)
ierr = cusparseSetMatType(bsrFloat%descrLU(2), CUSPARSE_MATRIX_TYPE_GENERAL)
ierr = cusparseSetMatFillMode(bsrFloat%descrLU(1), CUSPARSE_FILL_MODE_LOWER)
ierr = cusparseSetMatFillMode(bsrFloat%descrLU(2), CUSPARSE_FILL_MODE_UPPER)
ierr = cusparseSetMatDiagType(bsrFloat%descrLU(1), CUSPARSE_DIAG_TYPE_NON_UNIT)
ierr = cusparseSetMatDiagType(bsrFloat%descrLU(2), CUSPARSE_DIAG_TYPE_NON_UNIT)
#endif

bsrFloat%lAlloc = .TRUE.

END SUBROUTINE

SUBROUTINE createBsrDouble(bsrDouble, nbr, nbc, nBlock, blockSize)

IMPLICIT NONE

TYPE(BSR_DOUBLE) :: bsrDouble
INTEGER :: nbr, nbc, nBlock, blockSize
INTEGER :: ierr

IF (bsrDouble%lAlloc) CALL destroyBsr(bsrDouble)

bsrDouble%nbr = nbr
bsrDouble%nbc = nbc
bsrDouble%nBlock = nBlock
bsrDouble%blockSize = blockSize

ALLOCATE(bsrDouble%bsrVal(blockSize, blockSize, nBlock))
ALLOCATE(bsrDouble%bsrColIdx(nBlock))
ALLOCATE(bsrDouble%bsrRowPtr(nbr + 1))
bsrDouble%bsrVal = 0.0
bsrDouble%bsrColIdx = 0
bsrDouble%bsrRowPtr = 0

#ifdef __PGI
ierr = cusparseCreateBsrsv2Info(bsrDouble%info)
ierr = cusparseCreateBsrsv2Info(bsrDouble%infoLU(1))
ierr = cusparseCreateBsrsv2Info(bsrDouble%infoLU(2))
ierr = cusparseCreateMatDescr(bsrDouble%descr)
ierr = cusparseCreateMatDescr(bsrDouble%descrLU(1))
ierr = cusparseCreateMatDescr(bsrDouble%descrLU(2))
ierr = cusparseSetMatIndexBase(bsrDouble%descr, CUSPARSE_INDEX_BASE_ONE)
ierr = cusparseSetMatIndexBase(bsrDouble%descrLU(1), CUSPARSE_INDEX_BASE_ONE)
ierr = cusparseSetMatIndexBase(bsrDouble%descrLU(2), CUSPARSE_INDEX_BASE_ONE)
ierr = cusparseSetMatType(bsrDouble%descr, CUSPARSE_MATRIX_TYPE_GENERAL)
ierr = cusparseSetMatType(bsrDouble%descrLU(1), CUSPARSE_MATRIX_TYPE_GENERAL)
ierr = cusparseSetMatType(bsrDouble%descrLU(2), CUSPARSE_MATRIX_TYPE_GENERAL)
ierr = cusparseSetMatFillMode(bsrDouble%descrLU(1), CUSPARSE_FILL_MODE_LOWER)
ierr = cusparseSetMatFillMode(bsrDouble%descrLU(2), CUSPARSE_FILL_MODE_UPPER)
ierr = cusparseSetMatDiagType(bsrDouble%descrLU(1), CUSPARSE_DIAG_TYPE_NON_UNIT)
ierr = cusparseSetMatDiagType(bsrDouble%descrLU(2), CUSPARSE_DIAG_TYPE_NON_UNIT)
#endif

bsrDouble%lAlloc = .TRUE.

END SUBROUTINE

SUBROUTINE createBsrMixed(bsrMixed, nbr, nbc, nBlock, blockSize)

IMPLICIT NONE

TYPE(BSR_MIXED) :: bsrMixed
INTEGER :: nbr, nbc, nBlock, blockSize
INTEGER :: ierr

IF (bsrMixed%lAlloc) CALL destroyBsr(bsrMixed)

bsrMixed%nbr = nbr
bsrMixed%nbc = nbc
bsrMixed%nBlock = nBlock
bsrMixed%blockSize = blockSize

ALLOCATE(bsrMixed%bsrVal(blockSize, blockSize, nBlock))
ALLOCATE(bsrMixed%bsrColIdx(nBlock))
ALLOCATE(bsrMixed%bsrRowPtr(nbr + 1))
bsrMixed%bsrVal = 0.0
bsrMixed%bsrColIdx = 0
bsrMixed%bsrRowPtr = 0

#ifdef __PGI
ierr = cusparseCreateBsrsv2Info(bsrMixed%info)
ierr = cusparseCreateBsrsv2Info(bsrMixed%infoLU(1))
ierr = cusparseCreateBsrsv2Info(bsrMixed%infoLU(2))
ierr = cusparseCreateMatDescr(bsrMixed%descr)
ierr = cusparseCreateMatDescr(bsrMixed%descrLU(1))
ierr = cusparseCreateMatDescr(bsrMixed%descrLU(2))
ierr = cusparseSetMatIndexBase(bsrMixed%descr, CUSPARSE_INDEX_BASE_ONE)
ierr = cusparseSetMatIndexBase(bsrMixed%descrLU(1), CUSPARSE_INDEX_BASE_ONE)
ierr = cusparseSetMatIndexBase(bsrMixed%descrLU(2), CUSPARSE_INDEX_BASE_ONE)
ierr = cusparseSetMatType(bsrMixed%descr, CUSPARSE_MATRIX_TYPE_GENERAL)
ierr = cusparseSetMatType(bsrMixed%descrLU(1), CUSPARSE_MATRIX_TYPE_GENERAL)
ierr = cusparseSetMatType(bsrMixed%descrLU(2), CUSPARSE_MATRIX_TYPE_GENERAL)
ierr = cusparseSetMatFillMode(bsrMixed%descrLU(1), CUSPARSE_FILL_MODE_LOWER)
ierr = cusparseSetMatFillMode(bsrMixed%descrLU(2), CUSPARSE_FILL_MODE_UPPER)
ierr = cusparseSetMatDiagType(bsrMixed%descrLU(1), CUSPARSE_DIAG_TYPE_NON_UNIT)
ierr = cusparseSetMatDiagType(bsrMixed%descrLU(2), CUSPARSE_DIAG_TYPE_NON_UNIT)
#endif

bsrMixed%lAlloc = .TRUE.

END SUBROUTINE

SUBROUTINE clearBsrFloat(bsrFloat)

IMPLICIT NONE

TYPE(BSR_FLOAT) :: bsrFloat

bsrFloat%bsrVal = 0.0
bsrFloat%bsrColIdx = 0
bsrFloat%bsrRowPtr = 0

END SUBROUTINE

SUBROUTINE clearBsrDouble(bsrDouble)

IMPLICIT NONE

TYPE(BSR_DOUBLE) :: bsrDouble

bsrDouble%bsrVal = 0.0
bsrDouble%bsrColIdx = 0
bsrDouble%bsrRowPtr = 0

END SUBROUTINE

SUBROUTINE clearBsrMixed(bsrMixed)

IMPLICIT NONE

TYPE(BSR_MIXED) :: bsrMixed

bsrMixed%bsrVal = 0.0
bsrMixed%bsrColIdx = 0
bsrMixed%bsrRowPtr = 0

END SUBROUTINE

SUBROUTINE destroyBsrFloat(bsrFloat)

IMPLICIT NONE

TYPE(BSR_FLOAT) :: bsrFloat
INTEGER :: ierr

IF (.NOT. bsrFloat%lAlloc) RETURN

DEALLOCATE(bsrFloat%bsrVal)
DEALLOCATE(bsrFloat%bsrColIdx)
DEALLOCATE(bsrFloat%bsrRowPtr)

#ifdef __PGI
IF (bsrFloat%lDevice) THEN
  DEALLOCATE(bsrFloat%d_bsrVal)
  DEALLOCATE(bsrFloat%d_bsrColIdx)
  DEALLOCATE(bsrFloat%d_bsrRowPtr)
ENDIF

ierr = cusparseDestroyMatDescr(bsrFloat%descr)
ierr = cusparseDestroyMatDescr(bsrFloat%descrLU(1))
ierr = cusparseDestroyMatDescr(bsrFloat%descrLU(2))
#endif

bsrFloat%lAlloc = .FALSE.
bsrFloat%lDevice = .FALSE.

END SUBROUTINE

SUBROUTINE destroyBsrDouble(bsrDouble)

IMPLICIT NONE

TYPE(BSR_DOUBLE) :: bsrDouble
INTEGER :: ierr

IF (.NOT. bsrDouble%lAlloc) RETURN

DEALLOCATE(bsrDouble%bsrVal)
DEALLOCATE(bsrDouble%bsrColIdx)
DEALLOCATE(bsrDouble%bsrRowPtr)

#ifdef __PGI
IF (bsrDouble%lDevice) THEN
  DEALLOCATE(bsrDouble%d_bsrVal)
  DEALLOCATE(bsrDouble%d_bsrColIdx)
  DEALLOCATE(bsrDouble%d_bsrRowPtr)
ENDIF

ierr = cusparseDestroyMatDescr(bsrDouble%descr)
ierr = cusparseDestroyMatDescr(bsrDouble%descrLU(1))
ierr = cusparseDestroyMatDescr(bsrDouble%descrLU(2))
#endif

bsrDouble%lAlloc = .FALSE.
bsrDouble%lDevice = .FALSE.

END SUBROUTINE

SUBROUTINE destroyBsrMixed(bsrMixed)

IMPLICIT NONE

TYPE(BSR_MIXED) :: bsrMixed
INTEGER :: ierr

IF (.NOT. bsrMixed%lAlloc) RETURN

DEALLOCATE(bsrMixed%bsrVal)
DEALLOCATE(bsrMixed%bsrColIdx)
DEALLOCATE(bsrMixed%bsrRowPtr)

#ifdef __PGI
IF (bsrMixed%lDevice) THEN
  DEALLOCATE(bsrMixed%d_bsrVal)
  DEALLOCATE(bsrMixed%d_bsrVal8)
  DEALLOCATE(bsrMixed%d_bsrColIdx)
  DEALLOCATE(bsrMixed%d_bsrRowPtr)
ENDIF

ierr = cusparseDestroyMatDescr(bsrMixed%descr)
ierr = cusparseDestroyMatDescr(bsrMixed%descrLU(1))
ierr = cusparseDestroyMatDescr(bsrMixed%descrLU(2))
#endif

bsrMixed%lAlloc = .FALSE.
bsrMixed%lDevice = .FALSE.

END SUBROUTINE

#ifdef __PGI

SUBROUTINE copyDeviceBsrFloat(bsrFloat)

IMPLICIT NONE

TYPE(BSR_FLOAT) :: bsrFloat
INTEGER :: nbr, nbc, nBlock, blockSize

nbr = bsrFloat%nbr
nbc = bsrFloat%nbc
nBlock = bsrFloat%nBlock
blockSize = bsrFloat%blockSize

IF (.NOT. bsrFloat%lDevice) THEN
  ALLOCATE(bsrFloat%d_bsrVal(blockSize, blockSize, nBlock))
  ALLOCATE(bsrFloat%d_bsrColIdx(nBlock))
  ALLOCATE(bsrFloat%d_bsrRowPtr(nbr + 1))
  bsrFloat%lDevice = .TRUE.
ENDIF

bsrFloat%d_bsrVal = bsrFloat%bsrVal
bsrFloat%d_bsrColIdx = bsrFloat%bsrColIdx
bsrFloat%d_bsrRowPtr = bsrFloat%bsrRowPtr

END SUBROUTINE

SUBROUTINE copyDeviceBsrDouble(bsrDouble)

IMPLICIT NONE

TYPE(BSR_DOUBLE) :: bsrDouble
INTEGER :: nbr, nbc, nBlock, blockSize


nbr = bsrDouble%nbr
nbc = bsrDouble%nbc
nBlock = bsrDouble%nBlock
blockSize = bsrDouble%blockSize

IF (bsrDouble%lDevice) THEN
  ALLOCATE(bsrDouble%d_bsrVal(blockSize, blockSize, nBlock))
  ALLOCATE(bsrDouble%d_bsrColIdx(nBlock))
  ALLOCATE(bsrDouble%d_bsrRowPtr(nbr + 1))
  bsrDouble%lDevice = .TRUE.
ENDIF

bsrDouble%d_bsrVal = bsrDouble%bsrVal
bsrDouble%d_bsrColIdx = bsrDouble%bsrColIdx
bsrDouble%d_bsrRowPtr = bsrDouble%bsrRowPtr

END SUBROUTINE

SUBROUTINE copyDeviceBsrMixed(bsrMixed)

IMPLICIT NONE

TYPE(BSR_MIXED) :: bsrMixed
INTEGER :: nbr, nbc, nBlock, blockSize

nbr = bsrMixed%nbr
nbc = bsrMixed%nbc
nBlock = bsrMixed%nBlock
blockSize = bsrMixed%blockSize

IF (bsrMixed%lDevice) THEN
  ALLOCATE(bsrMixed%d_bsrVal(blockSize, blockSize, nBlock))
  ALLOCATE(bsrMixed%d_bsrVal8(blockSize, blockSize, nBlock))
  ALLOCATE(bsrMixed%d_bsrColIdx(nBlock))
  ALLOCATE(bsrMixed%d_bsrRowPtr(nbr + 1))
  bsrMixed%lDevice = .TRUE.
ENDIF

bsrMixed%d_bsrVal = bsrMixed%bsrVal
bsrMixed%d_bsrVal8 = bsrMixed%bsrVal
bsrMixed%d_bsrColIdx = bsrMixed%bsrColIdx
bsrMixed%d_bsrRowPtr = bsrMixed%bsrRowPtr

END SUBROUTINE

#endif

END MODULE
