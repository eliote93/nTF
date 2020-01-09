MODULE SenmMatSolver_mod
USE SP3DEF,   ONLY : SENMMAT_TYPE
USE MAT2x2OP, ONLY : SUBMATVECOP,  SubMatOP, MatInv
IMPLICIT NONE

CONTAINS
SUBROUTINE SenmMat_LU(Mat, N)
IMPLICIT NONE
TYPE(SENMMAT_TYPE) :: Mat
INTEGER :: n

INTEGER, POINTER :: nElmt(:), ElmtLoc(:, :)
REAL, POINTER :: Elmt(:, :, :)
INTEGER, POINTER :: LowerElmtIdx(:, :), nLowerElmt(:),DiagIdx(:)
!INTEGER, POINTER :: nLElmt(:), LElmtLoc(:, :)
!INTEGER, POINTER :: nUElmt(:), UElmtLoc(:, :)
REAL, POINTER :: LU(:, :, :)

INTEGER :: I, J, K
INTEGER :: m, irow, icol, irow0, icol0, IDX
INTEGER :: RowElmtIdx(n)
REAL :: Mult(4), DiagElmt(4), InvDiag(4), TEMP(4)

nElmt => Mat%nElmt; ElmtLoc => Mat%ElmtLoc
DiagIdx => Mat%DiagIdx; Elmt => Mat%Elmt
LowerElmtIdx => Mat%LowerElmtIdx; nLowerElmt => Mat%nLowerElmt
LU => Mat%LU

DO I = 1, N
  Mat%LU(:, :, I) = Mat%Elmt(:, :, I)
ENDDO
DO I = 1, N-1
  irow = i
  m = nLowerElmt(irow); idx = DiagIdx(irow)
  InvDiag = MatInv(LU(:, idx, irow))
  DO j = 1, m
    RowElmtIdx = 0
    irow0 = LowerElmtIdx(j, irow)
    DO k = 1, nElmt(irow0)
      icol0 = ElmtLoc(k, irow0)
      !IF(ICOL .LT. IROW) CYCLE
      RowElmtIdx(icol0) = k
    ENDDO
    !Finding Pivot
    Idx = RowElmtIdx(irow)
    Mult = SubMatOP(LU(:, Idx ,irow0), InvDiag)
    !Save Pivot Element
    LU(:, Idx, irow0) = Mult
    continue
    !Reduction
    DO k = DiagIdx(irow) + 1, nElmt(irow)
      ICOL = ElmtLoc(k, irow)
      !TEMP = LU(:, k, irow)
      IF(RowElmtIdx(ICOL) .EQ. 0) CYCLE
      TEMP = SubMatOP(MULT, LU(:, k, irow))
      IDX = RowElmtIdx(ICOL)
      LU(:, IDX, IROW0) = LU(:, IDX, IROW0) - TEMP 
    ENDDO
  ENDDO
ENDDO
!CALL Matrix_PRINT(MAT, N)
DO I = 1, N
  irow = i; idx = DiagIdx(irow)
  InvDiag = MatInv(LU(:, idx, irow))
  LU(:, idx, irow) = InvDiag
ENDDO
NULLIFY(nElmt, ElmtLoc, Elmt, DiagIdx)
NULLIFY(LowerElmtIdx, nLowerElmt, LU) 

END SUBROUTINE

SUBROUTINE SenmMatSolve_LU(A,  SOL, RHS, N)
IMPLICIT NONE
TYPE(SenmMat_TYPE) :: A
REAL :: RHS(2, n), SOL(2, n)
INTEGER :: N

INTEGER, POINTER :: nElmt(:), ElmtLoc(:, :), DiagIdx(:)
REAL, POINTER :: LU(:, :, :) 

INTEGER :: I, J, K, L, M
REAL :: Y(2,N),LMNT(2)

NELMT=> A%NELMT; ELMTLOC => A%ELMTLOC
DIAGIDX => A%DIAGIDX
LU => A%LU
!FORWARD
Y=0
Y(:, 1) = RHS(:, 1)
DO I = 2, N
  Y(:, I) = RHS(:, I)
  K = DIAGIDX(I)
  DO J = 1, K-1
    L = ELMTLOC(J, I)
    Y(:, I) = Y(:, I) - SUBMATVECOP(LU(:, J, I), Y(:, L))
  ENDDO
ENDDO
!BACKWORLD
K = DIAGIDX(N)
SOL(:, N) = SUBMATVECOP(LU(:, K, N), Y(:,N))
DO I = N-1, 1, -1
  LMNT = Y(:, I)
  K = DIAGIDX(I)
  M= NELMT(I)
  DO J = K+1, M
    L = ELMTLOC(J, I)
    LMNT = LMNT - SUBMATVECOP(LU(:, J, I), SOL(:, L))
  ENDDO
  SOL(:, I) = SUBMATVECOP(LU(:,K,I), LMNT)
ENDDO
!CHECK SOLUTION
DO i = 1, N
  M= NELMT(I)
  LMNT = 0
  DO J = 1, M
    L = ELMTLOC(J, I)
    LMNT = LMNT + SUBMATVECOP(A%ELMT(:, J, I), SOL(:, L))
  ENDDO
  Y(:, I) = LMNT
ENDDO
DO I = 1, N
  Y(:, I) = Y(:, I) - RHS(:, I)
ENDDO
!CONTINUE
!CALL Matrix_PRINT(A, N)
!write(43, '(A)') 'B=['
!DO i = 1,N
!  write(43, '(200(e40.12))') RHS(1, i) 
!  write(43, '(200(e40.12))') RHS(2, i) 
!ENDDO
!write(43, '(A)') '];'
!write(43, '(A)') 'X=['
!DO i = 1,N
!  write(43, '(200(e40.12))') SOL(1, i) 
!  write(43, '(200(e40.12))') SOL(2, i) 
!ENDDO
!write(43, '(A)') '];'
!STOP
END SUBROUTINE

SUBROUTINE Matrix_PRINT(Mat, N)
IMPLICIT NONE
TYPE(SenmMat_TYPE) :: Mat
INTEGER :: N
INTEGER :: I, J, K
INTEGER :: icol, irow, icol0, irow0
REAL :: A(2*N, 2*N), TEMP(4)
REAL :: L(2*N, 2*N), U(2*n, 2*n)
A = 0
L = 0
U = 0
DO i = 1, 2*N
  L(i, i) = 1
ENDDO
DO I = 1, N
  irow0 = i
  DO j = 1, Mat%nElmt(i)
    icol0 = Mat%ElmtLoc(j, i)
    TEMP = Mat%Elmt(:, j, i)
    icol = 2*(icol0) - 1; irow = 2*(irow0) - 1
    A(irow, icol) = Temp(1)
    icol = 2*(icol0); irow = 2*(irow0) -1
    A(irow, icol) = Temp(2)
    icol = 2*(icol0) - 1; irow = 2*(irow0) 
    A(irow, icol) = Temp(3)
    icol = 2*(icol0); irow = 2*(irow0)
    A(irow, icol) = Temp(4)   
  ENDDO
  k = MAT%DIagIdx(i)
  DO j = 1, k-1
    icol0 = Mat%ElmtLoc(j, i)
    TEMP = Mat%LU(:, j, i)
    icol = 2*(icol0) - 1; irow = 2*(irow0) - 1
    L(irow, icol) = Temp(1)
    icol = 2*(icol0); irow = 2*(irow0) -1
    L(irow, icol) = Temp(2)
    icol = 2*(icol0) - 1; irow = 2*(irow0) 
    L(irow, icol) = Temp(3)
    icol = 2*(icol0); irow = 2*(irow0)
    L(irow, icol) = Temp(4)   
  ENDDO
  
  DO j = k, Mat%nElmt(i)
    icol0 = Mat%ElmtLoc(j, i)
    TEMP = Mat%LU(:, j, i)
    icol = 2*(icol0) - 1; irow = 2*(irow0) - 1
    U(irow, icol) = Temp(1)
    icol = 2*(icol0); irow = 2*(irow0) -1
    U(irow, icol) = Temp(2)
    icol = 2*(icol0) - 1; irow = 2*(irow0) 
    U(irow, icol) = Temp(3)
    icol = 2*(icol0); irow = 2*(irow0)
    U(irow, icol) = Temp(4)   
  ENDDO  
ENDDO
write(43, '(A)') 'A=['
DO i = 1,2*N
  write(43, '(200(e40.12))') A(i, 1:2*N) 
ENDDO
write(43, '(A)') '];'
write(43, '(200(e40.12))')
write(43, '(A)') 'L=['
DO i = 1,2*N
  write(43, '(200(e40.12))') L(i, 1:2*N) 
ENDDO
write(43, '(A)') '];'
write(43, '(200(e40.12))')
write(43, '(A)') 'U=['
DO i = 1,2*N
  write(43, '(200(e40.12))') U(i, 1:2*N) 
ENDDO
write(43, '(A)') '];'
END SUBROUTINE


SUBROUTINE UpdtMatElmt(MAT, Elmt, irow, icol, ng)
IMPLICIT NONE
TYPE(SENMMAT_TYPE) :: MAT(NG)
REAL :: Elmt(4, ng)
INTEGER :: irow, icol, ng
INTEGER :: ig, j

DO IG = 1, NG
  Mat(IG)%nElmt(irow) = Mat(IG)%nElmt(irow) + 1
  j = Mat(IG)%nElmt(irow)
  Mat(IG)%ElmtLoc(j, irow) = icol
  Mat(IG)%Elmt(:, j, irow) = Elmt(:, ig)
  IF(IROW .EQ. ICOL) Mat(IG)%DiagIdx(irow) = j
  IF(IROW .GT. ICOL) THEN
    MAT(IG)%nLowerElmt(ICOL) = MAT(IG)%nLowerElmt(ICOL) + 1
    J = MAT(IG)%nLowerElmt(ICOL)
    MAT(IG)%LowerElmtIdx(J, ICOL) = IROW
  ENDIF
ENDDO
END SUBROUTINE

SUBROUTINE Init_SenmMAT(SENMMAT, NMESH, NG)
USE SP3DEF,      ONLY : SENMMAT_TYPE
IMPLICIT NONE
TYPE(SENMMAT_TYPE) :: SENMMAT(NG)
INTEGER :: NG, NMESH

INTEGER :: IG, N
n = 2*NMESH
DO IG = 1, NG
   SenmMat(IG)%nElmt = 0
   SENMMAT(IG)%ElmtLoc = 0
   SENMMAT(IG)%Elmt = 0
   SENMMAT(IG)%LowerElmtIdx = 0
   SENMMAT(IG)%nLowerElmt = 0
   SENMMAT(IG)%LU = 0   
ENDDO
END SUBROUTINE

SUBROUTINE Alloc_SENMMAT(SENMMAT, NMESH, NG)
USE SP3DEF,      ONLY : SENMMAT_TYPE
IMPLICIT NONE
TYPE(SENMMAT_TYPE) :: SENMMAT(NG)
INTEGER :: NG, NMESH

INTEGER :: IG, N
n = 2*NMESH
DO IG = 1, NG
   ALLOCATE(SenmMat(IG)%nElmt(N))
   ALLOCATE(SENMMAT(IG)%ElmtLoc(4, N))
   ALLOCATE(SENMMAT(IG)%DiagIdx(N))
   ALLOCATE(SENMMAT(IG)%Elmt(4, 4, N))
   ALLOCATE(SENMMAT(IG)%LowerElmtIdx(2, N))
   ALLOCATE(SENMMAT(IG)%nLowerElmt(N))
   ALLOCATE(SENMMAT(IG)%LU(4, 4, N))  
   SenmMat(IG)%nElmt = 0
   SENMMAT(IG)%ElmtLoc = 0
   SENMMAT(IG)%Elmt = 0
   SENMMAT(IG)%LowerElmtIdx = 0
   SENMMAT(IG)%nLowerElmt = 0
   SENMMAT(IG)%LU = 0   
ENDDO
END SUBROUTINE

SUBROUTINE FreeSENMMAT(SENMMAT, NG, NMESH)
USE SP3DEF,      ONLY : SENMMAT_TYPE
IMPLICIT NONE
TYPE(SENMMAT_TYPE) :: SENMMAT(NG)
INTEGER :: NG, NMESH
INTEGER :: IG

DO IG = 1, NG
   DEALLOCATE(SenmMat(IG)%nElmt)
   DEALLOCATE(SENMMAT(IG)%ElmtLoc)
   DEALLOCATE(SENMMAT(IG)%DiagIdx)
   DEALLOCATE(SENMMAT(IG)%Elmt)
   DEALLOCATE(SENMMAT(IG)%LowerElmtIdx)
   DEALLOCATE(SENMMAT(IG)%nLowerElmt) 
   DEALLOCATE(SENMMAT(IG)%LU)   
ENDDO

END SUBROUTINE


END MODULE