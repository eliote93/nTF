MODULE BasicOperation
USE PARAM
USE OMP_LIB
IMPLICIT NONE
!INTEGER :: i, j, k, l


INTERFACE CP_CA
  MODULE PROCEDURE CP_CA_1R
  MODULE PROCEDURE CP_CA_1I
  MODULE PROCEDURE CP_CA_2R
  MODULE PROCEDURE CP_CA_2I
  MODULE PROCEDURE CP_CA_3R
  MODULE PROCEDURE CP_CA_3I
  MODULE PROCEDURE CP_CA_4R
  MODULE PROCEDURE CP_CA_4I
END INTERFACE

INTERFACE CP_VA
  MODULE PROCEDURE CP_VA_1R
  MODULE PROCEDURE CP_VA_1I
  MODULE PROCEDURE CP_VA_2R
  MODULE PROCEDURE CP_VA_2I
  MODULE PROCEDURE CP_VA_3R
  MODULE PROCEDURE CP_VA_3I
  MODULE PROCEDURE CP_VA_4R
  MODULE PROCEDURE CP_VA_4I
END INTERFACE

INTERFACE MULTI_VA
  MODULE PROCEDURE MULTI_VA_1R
  MODULE PROCEDURE MULTI_VA_1I
  MODULE PROCEDURE MULTI_VA_2R
  MODULE PROCEDURE MULTI_VA_2I
  MODULE PROCEDURE MULTI_VA_3R
  MODULE PROCEDURE MULTI_VA_3I
  MODULE PROCEDURE MULTI_VA_4R
END INTERFACE

INTERFACE MULTI_CA
  MODULE PROCEDURE MULTI_CA_1R
  MODULE PROCEDURE MULTI_CA_1I
  MODULE PROCEDURE MULTI_CA_2R
  MODULE PROCEDURE MULTI_CA_2I
  MODULE PROCEDURE MULTI_CA_3R
  MODULE PROCEDURE MULTI_CA_3I
  MODULE PROCEDURE MULTI_CA_4R
END INTERFACE

INTERFACE DIV_VA
  MODULE PROCEDURE DIV_VA_1R
  !MODULE PROCEDURE MULTI_CA_1I
  MODULE PROCEDURE DIV_VA_2R
  !MODULE PROCEDURE MULTI_CA_2I
  MODULE PROCEDURE DIV_VA_3R
  !MODULE PROCEDURE MULTI_CA_3I
END INTERFACE

INTERFACE AD_VA
  MODULE PROCEDURE AD_VA_1R
  MODULE PROCEDURE AD_VA_1I
  MODULE PROCEDURE AD_VA_2R
  MODULE PROCEDURE AD_VA_2I
  MODULE PROCEDURE AD_VA_3R
  MODULE PROCEDURE AD_VA_3I
  MODULE PROCEDURE AD_VA_4R
  MODULE PROCEDURE AD_VA_4I
END INTERFACE

INTERFACE SUB_VA
  MODULE PROCEDURE SUB_VA_1R
  MODULE PROCEDURE SUB_VA_1I
  MODULE PROCEDURE SUB_VA_2R
  MODULE PROCEDURE SUB_VA_2I
  MODULE PROCEDURE SUB_VA_3R
  MODULE PROCEDURE SUB_VA_3I
  MODULE PROCEDURE SUB_VA_4R
  MODULE PROCEDURE SUB_VA_4I
END INTERFACE

INTERFACE CP_CAVB
  MODULE PROCEDURE CP_CAVB_1R
  MODULE PROCEDURE CP_CAVB_1I
  MODULE PROCEDURE CP_CAVB_2R
  MODULE PROCEDURE CP_CAVB_2I
  MODULE PROCEDURE CP_CAVB_3R
  MODULE PROCEDURE CP_CAVB_3I
END INTERFACE

INTERFACE DotProduct
  MODULE PROCEDURE DotProduct_1I
  MODULE PROCEDURE DotProduct_1R
  MODULE PROCEDURE DotProduct_2I
  MODULE PROCEDURE DotProduct_2R
  MODULE PROCEDURE DotProduct_3I
  MODULE PROCEDURE DotProduct_3R
END INTERFACE

INTERFACE CP_VA_OMP
  MODULE PROCEDURE CP_VA_OMP_2R
  MODULE PROCEDURE CP_VA_OMP_2I
END INTERFACE

INTERFACE CP_CA_OMP
  MODULE PROCEDURE CP_CA_OMP_2R
  MODULE PROCEDURE CP_CA_OMP_2I
END INTERFACE

INTERFACE AD_VA_OMP
  MODULE PROCEDURE AD_VA_OMP_2R
  MODULE PROCEDURE AD_VA_OMP_2I
END INTERFACE

INTERFACE MULTI_CA_OMP
  MODULE PROCEDURE MULTI_CA_OMP_2R
  MODULE PROCEDURE MULTI_CA_OMP_2I
END INTERFACE

CONTAINS

SUBROUTINE CP_CA_1R(A, r, n1)
REAL :: A(n1)
REAL :: r
INTEGER :: n1, n2
INTEGER :: i
DO i=1, n1
  A(i) = r
ENDDO
END SUBROUTINE

SUBROUTINE CP_CA_1I(A, r, n1)
INTEGER :: A(n1)
INTEGER :: r
INTEGER :: n1, n2
INTEGER :: i
DO i=1, n1
  A(i) = r
ENDDO
END SUBROUTINE

SUBROUTINE CP_CA_2R(A, r, n1, n2)
REAL :: A(n1, n2)
REAL :: r
INTEGER :: n1, n2
INTEGER :: i, j
DO i = 1, n2
  DO j = 1, n1
    A(j, i) = r
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE CP_CA_2I(A, r, n1, n2)
INTEGER :: A(n1, n2)
INTEGER :: r
INTEGER :: n1, n2
INTEGER :: i, j
DO i = 1, n2
  DO j = 1, n1
    A(j, i) = r
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE CP_CA_3R(A, r, n1, n2, n3)
REAL ::  A(n1, n2, n3)
REAL :: r
INTEGER :: n1, n2, n3
INTEGER :: i, j, k
DO i = 1, n3
  DO j = 1, n2
    DO k = 1, n1
      A(k, j, i) = r
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE CP_CA_3I(A, r, n1, n2, n3)
INTEGER :: A(n1, n2, n3)
INTEGER :: r
INTEGER :: n1, n2, n3
INTEGER :: i, j, k
DO i = 1, n3
  DO j = 1, n2
    DO k = 1, n1
      A(k, j, i) = r
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE CP_CA_4R(A, r, n1, n2, n3, n4)
REAL :: A(n1, n2, n3, n4)
REAL :: r
INTEGER :: n1, n2, n3, n4
INTEGER :: i, j, k, l
DO l = 1, n4
  DO k = 1, n3
    DO j = 1, n2
      DO i = 1, n1
        A(i, j, k, l) = r
      ENDDO
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE CP_CA_4I(A, r, n1, n2, n3, n4)
INTEGER :: A(n1, n2, n3, n4)
INTEGER :: r
INTEGER :: n1, n2, n3, n4
INTEGER :: i, j, k, l
DO l = 1, n4
  DO k = 1, n3
    DO j = 1, n2
      DO i = 1, n1
        A(i, j, k, l) = r
      ENDDO
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CP_VA_1R(A, B, n1)
REAL :: A(n1), B(n1)
INTEGER :: n1, n2
INTEGER :: i
CONTINUE

DO i=1, n1
  A(i) = B(i)
ENDDO
END SUBROUTINE

SUBROUTINE CP_VA_1I(A, B, n1)
INTEGER :: A(n1), B(n1)
INTEGER :: n1, n2
INTEGER :: i
DO i=1, n1
  A(i) = B(i)
ENDDO
END SUBROUTINE

SUBROUTINE CP_VA_2R(A, B, n1, n2)
REAL :: A(n1, n2), B(n1, n2)
INTEGER :: n1, n2
INTEGER :: i, j
CONTINUE
DO i = 1, n2
  DO j = 1, n1
    A(j, i) = B(j, i)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE CP_VA_2I(A, B, n1, n2)
INTEGER :: A(n1, n2), B(n1, n2)
INTEGER :: n1, n2
INTEGER :: i, j
CONTINUE
DO i = 1, n2
  DO j = 1, n1
    A(j, i) = B(j, i)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE CP_VA_3R(A, B, n1, n2, n3)
REAL :: A(n1, n2, n3), B(n1, n2, n3)
INTEGER :: n1, n2, n3
INTEGER :: i, j, k
DO k = 1, n3
  DO j = 1, n2
    DO i = 1, n1
      A(i, j, k) = B(i, j, k)
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE CP_VA_3I(A, B, n1, n2, n3)
INTEGER :: A(n1, n2, n3), B(n1, n2, n3)
INTEGER :: n1, n2, n3
INTEGER :: i, j, k
DO k = 1, n3
  DO j = 1, n2
    DO i = 1, n1
      A(i, j, k) = B(i, j, k)
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE CP_VA_4R(A, B, n1, n2, n3, n4)
REAL :: A(n1, n2, n3, n4), B(n1, n2, n3, n4)
INTEGER :: n1, n2, n3, n4
INTEGER :: i, j, k, l
DO l = 1, n4
  DO k = 1, n3
    DO j = 1, n2
      DO i = 1, n1
        A(i, j, k, l) = B(i, j, k, l)
      ENDDO
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE CP_VA_4I(A, B, n1, n2, n3, n4)
INTEGER :: A(n1, n2, n3, n4), B(n1, n2, n3, n4)
INTEGER :: n1, n2, n3, n4
INTEGER :: i, j, k, l
DO l = 1, n4
  DO k = 1, n3
    DO j = 1, n2
      DO i = 1, n1
        A(i, j, k, l) = B(i, j, k, l)
      ENDDO
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE MULTI_VA_1R(A, B, n1)
REAL :: A(n1), B(n1)
INTEGER :: n1, n2
INTEGER :: i, j
CONTINUE
DO i=1, n1
  B(i) = A(i) * B(i)
ENDDO
END SUBROUTINE

SUBROUTINE MULTI_VA_2R(A, B, n1, n2)
REAL :: A(n1, n2), B(n1, n2)
INTEGER :: n1, n2
INTEGER :: i, j
CONTINUE
DO i = 1, n2
  DO j = 1, n1
    B(j, i) = A(j, i) * B(j, i)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE MULTI_VA_3R(A, B, n1, n2, n3)
REAL :: A(n1, n2, n3), B(n1, n2, n3)
INTEGER :: n1, n2, n3
INTEGER :: i, j, k
DO i = 1, n3
  DO j = 1, n2
    DO k = 1, n1
      B(k, j, i) = A(k, j, i) * B(k, j, i)
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE


SUBROUTINE MULTI_VA_4R(A, B, n1, n2, n3, n4)
REAL :: A(n1, n2, n3, n4), B(n1, n2, n3, n4)
INTEGER :: n1, n2, n3, n4
INTEGER :: i, j, k, l
DO l = 1, n4
  DO i = 1, n3
    DO j = 1, n2
      DO k = 1, n1
        B(k, j, i, l) = A(k, j, i, l) * B(k, j, i, l)
      ENDDO
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE MULTI_VA_1I(A, B, n1)
INTEGER :: A(n1), B(n1)
INTEGER :: n1, n2
INTEGER :: i, j
CONTINUE
DO i=1, n1
  B(i) = A(i) * B(i)
ENDDO
END SUBROUTINE

SUBROUTINE MULTI_VA_2I(A, B, n1, n2)
INTEGER :: A(n1, n2), B(n1, n2)
INTEGER :: n1, n2
INTEGER :: i, j
CONTINUE
DO i = 1, n2
  DO j = 1, n1
    B(j, i) = A(j, i) * B(j, i)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE MULTI_VA_3I(A, B, n1, n2, n3)
INTEGER :: A(n1, n2, n3), B(n1, n2, n3)
INTEGER :: n1, n2, n3
INTEGER :: i, j, k
DO i = 1, n3
  DO j = 1, n2
    DO k = 1, n1
      B(k, j, i) = A(k, j, i) * B(k, j, i)
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE MULTI_CA_1R(R, A, n1)
REAL :: R, A(n1)
INTEGER :: n1
INTEGER :: i
CONTINUE
DO i=1, n1
  A(i) = A(i) * R
ENDDO
END SUBROUTINE

SUBROUTINE MULTI_CA_2R(R, A, n1, n2)
REAL :: R, A(n1, n2)
INTEGER :: n1, n2
INTEGER :: i, j
CONTINUE
DO i = 1, n2
  DO j = 1, n1
    A(j, i) = A(j, i) * R
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE MULTI_CA_3R(R, A, n1, n2, n3)
REAL :: A(n1, n2, n3), R
INTEGER :: n1, n2, n3
INTEGER :: i, j, k
DO i = 1, n3
  DO j = 1, n2
    DO k = 1, n1
      A(k, j, i) = A(k, j, i) * R
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE MULTI_CA_4R(R, A, n1, n2, n3, n4)
REAL :: A(n1, n2, n3, n4), R
INTEGER :: n1, n2, n3, n4
INTEGER :: i, j, k, l
DO l = 1, n4
  DO i = 1, n3
    DO j = 1, n2
      DO k = 1, n1
        A(k, j, i, l) = A(k, j, i, l) * R
      ENDDO
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE
!
SUBROUTINE MULTI_CA_1I(R, A, n1)
INTEGER :: A(n1), R
INTEGER :: n1, n2
INTEGER :: i, j
CONTINUE
DO i=1, n1
  A(i) = A(i) * R
ENDDO
END SUBROUTINE

SUBROUTINE MULTI_CA_2I(R, A, n1, n2)
INTEGER :: A(n1, n2), R
INTEGER :: n1, n2
INTEGER :: i, j
CONTINUE
DO i = 1, n2
  DO j = 1, n1
    A(j, i) = A(j, i) * R
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE MULTI_CA_3I(R, A, n1, n2, n3)
INTEGER :: R, A(n1, n2, n3)
INTEGER :: n1, n2, n3
INTEGER :: i, j, k
DO i = 1, n3
  DO j = 1, n2
    DO k = 1, n1
      A(k, j, i) = A(k, j, i) * R
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE DIV_VA_1R(A, B, n1)
REAL :: A(n1), B(n1)
INTEGER :: n1
INTEGER :: i
CONTINUE
DO i=1, n1
  A(i) = A(i) / B(i)
ENDDO
END SUBROUTINE

SUBROUTINE DIV_VA_2R(A, B, n1, n2)
REAL :: A(n1, n2), B(n1, n2)
INTEGER :: n1, n2
INTEGER :: i, j
CONTINUE
DO i = 1, n2
  DO j = 1, n1
    A(j, i) = A(j, i) / B(j, i)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE DIV_VA_3R(A, B, n1, n2, n3)
REAL :: A(n1, n2, n3), B(n1, n2, n3)
INTEGER :: n1, n2, n3
INTEGER :: i, j, k
DO i = 1, n3
  DO j = 1, n2
    DO k = 1, n1
      B(k, j, i) = A(k, j, i) / B(k, j, i)
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE AD_VA_1R(C, A, B, n1)
REAL :: C(n1), A(n1), B(n1)
INTEGER :: n1, n2
INTEGER :: i, j
CONTINUE
DO i=1, n1
  C(i) = A(i) + B(i)
ENDDO
END SUBROUTINE

SUBROUTINE AD_VA_1I(C, A, B, n1)
INTEGER :: C(n1), A(n1), B(n1)
INTEGER :: n1, n2
INTEGER :: i, j
DO i=1, n1
  C(i) = A(i) + B(i)
ENDDO
END SUBROUTINE

SUBROUTINE AD_VA_2R(C, A, B, n1, n2)
REAL :: C(n1, n2), A(n1, n2), B(n1, n2)
INTEGER :: n1, n2
INTEGER :: i, j
CONTINUE
DO i = 1, n2
  DO j = 1, n1
    C(j, i) = A(j, i) + B(j, i)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE AD_VA_2I(C, A, B, n1, n2)
INTEGER :: C(n1, n2), A(n1, n2), B(n1, n2)
INTEGER :: n1, n2
INTEGER :: i, j
CONTINUE
DO i = 1, n2
  DO j = 1, n1
    C(j, i) = A(j, i) + B(j, i)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE AD_VA_3R(C, A, B, n1, n2, n3)
REAL :: C(n1, n2, n3), A(n1, n2, n3), B(n1, n2, n3)
INTEGER :: n1, n2, n3
INTEGER :: i, j, k
DO k = 1, n3
  DO j = 1, n2
    DO i = 1, n1
      C(i, j, k) = A(i, j, k) + B(i, j, k)
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE AD_VA_3I(C, A, B, n1, n2, n3)
INTEGER :: C(n1, n2, n3), A(n1, n2, n3), B(n1, n2, n3)
INTEGER :: n1, n2, n3
INTEGER :: i, j, k
DO k = 1, n3
  DO j = 1, n2
    DO i = 1, n1
      C(i, j, k) = A(i, j, k) + B(i, j, k)
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE AD_VA_4R(C, A, B, n1, n2, n3, n4)
REAL :: C(n1, n2, n3, n4), A(n1, n2, n3, n4), B(n1, n2, n3, n4)
INTEGER :: n1, n2, n3, n4
INTEGER :: i, j, k, l
DO l = 1, n4
  DO k = 1, n3
    DO j = 1, n2
      DO i = 1, n1
        C(i, j, k, l) = A(i, j, k, l) + B(i, j, k, l)
      ENDDO
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE AD_VA_4I(C, A, B, n1, n2, n3, n4)
INTEGER :: C(n1, n2, n3, n4), A(n1, n2, n3, n4), B(n1, n2, n3, n4)
INTEGER :: n1, n2, n3, n4
INTEGER :: i, j, k, l
DO l = 1, n4
  DO k = 1, n3
    DO j = 1, n2
      DO i = 1, n1
        C(i, j, k, l) = A(i, j, k, l) + B(i, j, k, l)
      ENDDO
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SUB_VA_1R(C, A, B, n1)
REAL :: C(n1), A(n1), B(n1)
INTEGER :: n1, n2
INTEGER :: i
CONTINUE
DO i=1, n1
  C(i) = A(i) - B(i)
ENDDO
END SUBROUTINE

SUBROUTINE SUB_VA_1I(C, A, B, n1)
INTEGER :: C(n1), A(n1), B(n1)
INTEGER :: n1, n2
INTEGER :: i
DO i=1, n1
  C(i) = A(i) - B(i)
ENDDO
END SUBROUTINE

SUBROUTINE SUB_VA_2R(C, A, B, n1, n2)
REAL :: C(n1, n2), A(n1, n2), B(n1, n2)
INTEGER :: n1, n2
INTEGER :: i, j
CONTINUE
DO i = 1, n2
  DO j = 1, n1
    C(j, i) = A(j, i) - B(j, i)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE SUB_VA_2I(C, A, B, n1, n2)
INTEGER :: C(n1, n2), A(n1, n2), B(n1, n2)
INTEGER :: n1, n2
INTEGER :: i, j
CONTINUE
DO i = 1, n2
  DO j = 1, n1
    C(j, i) = A(j, i) - B(j, i)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE SUB_VA_3R(C, A, B, n1, n2, n3)
REAL :: C(n1, n2, n3), A(n1, n2, n3), B(n1, n2, n3)
INTEGER :: n1, n2, n3
INTEGER :: i, j, k
DO k = 1, n3
  DO j = 1, n2
    DO i = 1, n1
      C(i, j, k) = A(i, j, k) - B(i, j, k)
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE SUB_VA_3I(C, A, B, n1, n2, n3)
INTEGER :: C(n1, n2, n3), A(n1, n2, n3), B(n1, n2, n3)
INTEGER :: n1, n2, n3
INTEGER :: i, j, k
DO k = 1, n3
  DO j = 1, n2
    DO i = 1, n1
      C(i, j, k) = A(i, j, k) - B(i, j, k)
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE SUB_VA_4R(C, A, B, n1, n2, n3, n4)
REAL :: C(n1, n2, n3, n4), A(n1, n2, n3, n4), B(n1, n2, n3, n4)
INTEGER :: n1, n2, n3, n4
INTEGER :: i, j, k, l
DO l = 1, n4
  DO k = 1, n3
    DO j = 1, n2
      DO i = 1, n1
        C(i, j, k, l) = A(i, j, k, l) - B(i, j, k, l)
      ENDDO
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE SUB_VA_4I(C, A, B, n1, n2, n3, n4)
INTEGER :: C(n1, n2, n3, n4), A(n1, n2, n3, n4), B(n1, n2, n3, n4)
INTEGER :: n1, n2, n3, n4
INTEGER :: i, j, k, l
DO l = 1, n4
  DO k = 1, n3
    DO j = 1, n2
      DO i = 1, n1
        C(i, j, k, l) = A(i, j, k, l) - B(i, j, k, l)
      ENDDO
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CP_CAVB_1I(r, A, B, n1)
IMPLICIT NONE
INTEGER :: r
INTEGER :: A(n1), B(n1)
INTEGER :: n1
INTEGER :: i
DO i = 1, n1
  B(i) =r * A(i)
ENDDO
END SUBROUTINE

SUBROUTINE CP_CAVB_2I(r, A, B, n1, n2)
IMPLICIT NONE
INTEGER :: r
INTEGER :: A(n1, n2), B(n1, n2)
INTEGER :: n1, n2
INTEGER :: i, j
DO j = 1 ,n2
  DO i = 1, n1
    B(i, j) = r * A(i, j)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE CP_CAVB_3I(r, A, B, n1, n2, n3)
IMPLICIT NONE
INTEGER :: r
INTEGER :: A(n1, n2, n3), B(n1, n2, n3)
INTEGER :: n1, n2, n3
INTEGER :: i, j, k
DO k = 1, n3
  DO j = 1 ,n2
    DO i = 1, n1
      B(i, j, k) = r * A(i, j, k)
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE CP_CAVB_1R(r, A, B, n1)
IMPLICIT NONE
REAL :: r
REAL :: A(n1), B(n1)
INTEGER :: n1
INTEGER :: i
DO i = 1, n1
  B(i) = r * A(i)
ENDDO
END SUBROUTINE

SUBROUTINE CP_CAVB_2R(r, A, B, n1, n2)
IMPLICIT NONE
REAL :: r
REAL :: A(n1, n2), B(n1, n2)
INTEGER :: n1, n2
INTEGER :: i, j
DO j = 1 ,n2
  DO i = 1, n1
    B(i, j) = r * A(i, j)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE CP_CAVB_3R(r, A, B, n1, n2, n3)
IMPLICIT NONE
REAL :: r
REAL :: A(n1, n2, n3), B(n1, n2, n3)
INTEGER :: n1, n2, n3
INTEGER :: i, j, k
DO k = 1, n3
  DO j = 1 ,n2
    DO i = 1, n1
      B(i, j, k) = r * A(i, j, k)
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION DotProduct_1I(x, y, n1)
IMPLICIT NONE
INTEGER :: x(n1), y(n1)
INTEGER :: n1
INTEGER :: i, j, k
INTEGER :: DotProduct_1I
DotProduct_1I = 0
DO i = 1, n1
  DotProduct_1I = DotProduct_1I + x(i) * y(i)
ENDDO
END FUNCTION

FUNCTION DotProduct_1R(x, y, n1)
IMPLICIT NONE
REAL :: x(n1), y(n1)
INTEGER :: n1
INTEGER :: i, j, k
REAL :: DotProduct_1R
DotProduct_1R = 0
DO i = 1, n1
  DotProduct_1R = DotProduct_1R + x(i) * y(i)
ENDDO

END FUNCTION

FUNCTION DotProduct_2I(x, y, n1, n2)
IMPLICIT NONE
INTEGER :: x(n1, n2), y(n1, n2)
INTEGER :: n1, n2
INTEGER :: i, j, k
INTEGER :: DotProduct_2I
DotProduct_2I = 0
DO j = 1, n2
  DO i = 1, n1
    DotProduct_2I = DotProduct_2I + x(i, j) * y(i, j)
  ENDDO
ENDDO
END FUNCTION

FUNCTION DotProduct_2R(x, y, n1, n2)
IMPLICIT NONE
REAL :: x(n1, n2), y(n1, n2)
INTEGER :: n1, n2
INTEGER :: i, j, k
REAL :: DotProduct_2R
DotProduct_2R = 0
DO j = 1, n2
  DO i = 1, n1
    DotProduct_2R = DotProduct_2R + x(i, j) * y(i, j)
  ENDDO
ENDDO

END FUNCTION

FUNCTION DotProduct_3I(x, y, n1, n2, n3)
IMPLICIT NONE
INTEGER :: x(n1, n2, n3), y(n1, n2, n3)
INTEGER :: n1, n2, n3
INTEGER :: i, j, k
INTEGER :: DotProduct_3I
DotProduct_3I = 0
DO k = 1, n3
  DO j = 1, n2
    DO i = 1, n1
      DotProduct_3I = DotProduct_3I + x(i, j, k) * y(i, j, k)
    ENDDO
  ENDDO
ENDDO
END FUNCTION

FUNCTION DotProduct_3R(x, y, n1, n2, n3)
IMPLICIT NONE
REAL :: x(n1, n2, n3), y(n1, n2, n3)
INTEGER :: n1, n2, n3
INTEGER :: i, j, k
REAL :: DotProduct_3R
DotProduct_3R = 0
DO k = 1, n3
  DO j = 1, n2
    DO i = 1, n1
      DotProduct_3R = DotProduct_3R + x(i, j, k) * y(i, j, k)
    ENDDO
  ENDDO
ENDDO
END FUNCTION


SUBROUTINE AD_VA_OMP_2R(C, A, B, n1, n2, nThread)
REAL :: C(n1, n2), A(n1, n2), B(n1, n2)
INTEGER :: n1, n2, nThread
INTEGER :: i, j
CONTINUE
!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(nThread)
!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(j)
DO i = 1, n2
!$OMP DO
  DO j = 1, n1
    C(j, i) = A(j, i) + B(j, i)
  ENDDO
!$OMP END DO
ENDDO
!$OMP END PARALLEL
END SUBROUTINE

SUBROUTINE AD_VA_OMP_2I(C, A, B, n1, n2, nThread)
INTEGER :: C(n1, n2), A(n1, n2), B(n1, n2)
INTEGER :: n1, n2, nThread
INTEGER :: i, j
CONTINUE
!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(nThread)
!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(j)
DO i = 1, n2
!$OMP DO
  DO j = 1, n1
    C(j, i) = A(j, i) + B(j, i)
  ENDDO
!$OMP END DO
ENDDO
!$OMP END PARALLEL
END SUBROUTINE

SUBROUTINE CP_CA_OMP_2R(A, r, n1, n2, nThread)
REAL :: A(n1, n2)
REAL :: r
INTEGER :: n1, n2, nThread
INTEGER :: i, j
!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(nThread)
!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(j)
DO i = 1, n2
!$OMP DO
  DO j = 1, n1
    A(j, i) = r
  ENDDO
!$OMP END DO
ENDDO
!$OMP END PARALLEL
END SUBROUTINE

SUBROUTINE CP_CA_OMP_2I(A, r, n1, n2, nThread)
INTEGER :: A(n1, n2)
INTEGER :: r
INTEGER :: n1, n2, nThread
INTEGER :: i, j
!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(nThread)
!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(j)
DO i = 1, n2
!$OMP DO
  DO j = 1, n1
    A(j, i) = r
  ENDDO
!$OMP END DO
ENDDO
!$OMP END PARALLEL
END SUBROUTINE

SUBROUTINE MULTI_CA_OMP_2R(R, A, n1, n2, nThread)
REAL :: R, A(n1, n2)
INTEGER :: n1, n2, nThread
INTEGER :: i, j
CONTINUE
!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(nThread)
!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(j)
DO i = 1, n2
!$OMP DO
  DO j = 1, n1
    A(j, i) = A(j, i) * R
  ENDDO
!$OMP END DO
ENDDO
!$OMP END PARALLEL
END SUBROUTINE

SUBROUTINE CP_VA_OMP_2R(A, B, n1, n2, nThread)
REAL :: A(n1, n2), B(n1, n2)
INTEGER :: n1, n2, nThread
INTEGER :: i, j
CONTINUE
!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(nThread)
!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(j)
DO i = 1, n2
!$OMP DO
  DO j = 1, n1
    A(j, i) = B(j, i)
  ENDDO
!$OMP END DO
ENDDO
!$OMP END PARALLEL
END SUBROUTINE

SUBROUTINE CP_VA_OMP_2I(A, B, n1, n2, nThread)
INTEGER :: A(n1, n2), B(n1, n2)
INTEGER :: n1, n2, nThread
INTEGER :: i, j
CONTINUE
!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(nThread)
!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(j)
DO i = 1, n2
!$OMP DO
  DO j = 1, n1
    A(j, i) = B(j, i)
  ENDDO
!$OMP END DO
ENDDO
!$OMP END PARALLEL
END SUBROUTINE

SUBROUTINE MULTI_CA_OMP_2I(R, A, n1, n2, nThread)
INTEGER :: R, A(n1, n2)
INTEGER :: n1, n2, nThread
INTEGER :: i, j
CONTINUE
!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(nThread)
!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(j)
DO i = 1, n2
!$OMP DO
  DO j = 1, n1
    A(j, i) = A(j, i) * R
  ENDDO
!$OMP END DO
ENDDO
!$OMP END PARALLEL
END SUBROUTINE

END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
