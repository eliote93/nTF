MODULE AuxilCSR
  USE CSRMATRIX
  IMPLICIT NONE
  
  CONTAINS
  
  SUBROUTINE Full2CSR(A, A_csr)
  IMPLICIT NONE
  REAL(8), INTENT(IN) :: A(:,:)
  TYPE(CSR_DOUBLE), INTENT(INOUT) :: A_csr
  
  INTEGER :: nr, nc, nnz
  INTEGER :: i, j
  
  REAL(8), ALLOCATABLE :: val(:)
  INTEGER, ALLOCATABLE :: rowptr(:), colIdx(:)
  
  nr = SIZE(A, 1); nc = SIZE(A, 2)
  nnz = 0
  DO i = 1, nr
    DO j = 1, nc
      IF (A(i,j) .NE. 0._8) nnz = nnz+1
    END DO
  END DO
  
  CALL createCSR(A_csr, nnz, nr, nc)
  
  nnz = 0;
  rowptr(:) = 0
  DO i = 1, nr
    DO j = 1, nc
      IF (A(i,j) .NE. 0._8) CALL pushCsr(A_Csr, A(i,j), i, j)
    END DO
  END DO
  
  CALL finalizeCSR(A_Csr, .FALSE.)
  
  END SUBROUTINE
  
  SUBROUTINE Full2CSRZ(A, A_csr)
  IMPLICIT NONE
  REAL(8), INTENT(IN) :: A(:,:)
  TYPE(CSR_DOUBLE_COMPLEX), INTENT(INOUT) :: A_csr
  
  INTEGER :: nr, nc, nnz
  INTEGER :: i, j
  
  COMPLEX(8), ALLOCATABLE :: val(:)
  INTEGER, ALLOCATABLE :: rowptr(:), colIdx(:)
  
  nr = SIZE(A, 1); nc = SIZE(A, 2)
  nnz = 0
  DO i = 1, nr
    DO j = 1, nc
      IF (A(i,j) .NE. 0._8) nnz = nnz+1
    END DO
  END DO
  
  CALL createCSR(A_csr, nnz, nr, nc)
  
  nnz = 0;
  rowptr(:) = 0
  DO i = 1, nr
    DO j = 1, nc
      IF (A(i,j) .NE. 0._8) CALL pushCsr(A_Csr, A(i,j), i, j)
    END DO
  END DO
  
  CALL finalizeCSR(A_Csr, .FALSE.)
  
  END SUBROUTINE
  
  SUBROUTINE FullZ2CSRZ(A, A_csr)
  IMPLICIT NONE
  COMPLEX(8), INTENT(IN) :: A(:,:)
  TYPE(CSR_DOUBLE_COMPLEX), INTENT(INOUT) :: A_csr
  
  INTEGER :: nr, nc, nnz
  INTEGER :: i, j
  
  COMPLEX(8), ALLOCATABLE :: val(:)
  INTEGER, ALLOCATABLE :: rowptr(:), colIdx(:)
  
  nr = SIZE(A, 1); nc = SIZE(A, 2)
  nnz = 0
  DO i = 1, nr
    DO j = 1, nc
      IF (A(i,j) .NE. 0._8) nnz = nnz+1
    END DO
  END DO
  
  CALL createCSR(A_csr, nnz, nr, nc)
  
  nnz = 0;
  rowptr(:) = 0
  DO i = 1, nr
    DO j = 1, nc
      IF (A(i,j) .NE. 0._8) CALL pushCsr(A_Csr, A(i,j), i, j)
    END DO
  END DO
  
  CALL finalizeCSR(A_Csr, .FALSE.)
  
  END SUBROUTINE
  
  SUBROUTINE CSR2Full(A, A_csr)
  IMPLICIT NONE
  REAL(8), INTENT(INOUT) :: A(:,:)
  TYPE(CSR_DOUBLE), INTENT(IN) :: A_csr
  
  INTEGER :: nr, nc, nnz
  INTEGER :: i, j
  
  INTEGER, ALLOCATABLE :: rowptr(:), colIdx(:)
  REAL(8), ALLOCATABLE :: val(:)
  
  nr = A_csr%nr; nc = A_csr%nc; nnz= A_csr%nnz
  
  ALLOCATE(rowptr(nr+1), colIdx(nnz), val(nnz))
  
  rowptr(:) = A_csr%csrRowPtr(:); colIdx(:) = A_csr%csrColIdx(:)
  val(:) = A_csr%csrVal(:)
  
  A(:,:) = 0._8
  
  nr = A_csr%nr; nc = A_csr%nc; nnz = A_csr%nnz;
  DO i = 1, nr
    DO j = rowptr(i),rowptr(i+1)-1
      A(i, colIdx(j)) = val(j)
    END DO
  END DO
  DEALLOCATE(rowptr, colIdx, val)
  END SUBROUTINE
  
  SUBROUTINE CSRZ2Full(A, A_csr)
  IMPLICIT NONE
  REAL(8), INTENT(INOUT) :: A(:,:)
  TYPE(CSR_DOUBLE_COMPLEX), INTENT(IN) :: A_csr
  
  INTEGER :: nr, nc, nnz
  INTEGER :: i, j
  
  INTEGER, ALLOCATABLE :: rowptr(:), colIdx(:)
  COMPLEX(8), ALLOCATABLE :: val(:)
  
  nr = A_csr%nr; nc = A_csr%nc; nnz= A_csr%nnz
  
  ALLOCATE(rowptr(nr+1), colIdx(nnz), val(nnz))
  
  rowptr(:) = A_csr%csrRowPtr(:); colIdx(:) = A_csr%csrColIdx(:)
  val(:) = A_csr%csrVal(:)
  
  A(:,:) = 0._8
  
  nr = A_csr%nr; nc = A_csr%nc; nnz = A_csr%nnz;
  DO i = 1, nr
    DO j = rowptr(i),rowptr(i+1)-1
      A(i, colIdx(j)) = dble(val(j))
    END DO
  END DO
  DEALLOCATE(rowptr, colIdx, val)
  END SUBROUTINE
  
  SUBROUTINE CSR2FullZ(A, A_csr)
  IMPLICIT NONE
  COMPLEX(8), INTENT(INOUT) :: A(:,:)
  TYPE(CSR_DOUBLE), INTENT(IN) :: A_csr
  
  INTEGER :: nr, nc, nnz
  INTEGER :: i, j
  
  INTEGER, ALLOCATABLE :: rowptr(:), colIdx(:)
  REAL(8), ALLOCATABLE :: val(:)
  
  nr = A_csr%nr; nc = A_csr%nc; nnz= A_csr%nnz
  
  ALLOCATE(rowptr(nr+1), colIdx(nnz), val(nnz))
  
  rowptr(:) = A_csr%csrRowPtr(:); colIdx(:) = A_csr%csrColIdx(:)
  val(:) = A_csr%csrVal(:)
  
  A(:,:) = 0._8
  
  nr = A_csr%nr; nc = A_csr%nc; nnz = A_csr%nnz;
  DO i = 1, nr
    DO j = rowptr(i),rowptr(i+1)-1
      A(i, colIdx(j)) = val(j)
    END DO
  END DO
  DEALLOCATE(rowptr, colIdx, val)
  END SUBROUTINE
  
  SUBROUTINE CSRZ2FullZ(A, A_csr)
  IMPLICIT NONE
  COMPLEX(8), INTENT(INOUT) :: A(:,:)
  TYPE(CSR_DOUBLE_COMPLEX), INTENT(IN) :: A_csr
  
  INTEGER :: nr, nc, nnz
  INTEGER :: i, j
  
  INTEGER, ALLOCATABLE :: rowptr(:), colIdx(:)
  COMPLEX(8), ALLOCATABLE :: val(:)
  
  nr = A_csr%nr; nc = A_csr%nc; nnz= A_csr%nnz
  
  ALLOCATE(rowptr(nr+1), colIdx(nnz), val(nnz))
  
  rowptr(:) = A_csr%csrRowPtr(:); colIdx(:) = A_csr%csrColIdx(:)
  val(:) = A_csr%csrVal(:)
  
  A(:,:) = 0._8
  
  nr = A_csr%nr; nc = A_csr%nc; nnz = A_csr%nnz;
  DO i = 1, nr
    DO j = rowptr(i),rowptr(i+1)-1
      A(i, colIdx(j)) = val(j)
    END DO
  END DO
  DEALLOCATE(rowptr, colIdx, val)
  END SUBROUTINE
  
END MODULE