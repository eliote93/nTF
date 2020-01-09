MODULE LU_MOD
USE PARAM
IMPLICIT NONE

INTEGER, PARAMETER :: DB=8

CONTAINS

SUBROUTINE DirLU1Dsolver(Diag, L, U, x, b, n)
REAL, INTENT(InOUT) :: Diag(n), L(n), U(n), b(n)
REAL, INTENT(InOut) :: x(n)
INTEGER :: n

REAL :: pivot
INTEGER :: i, j

j = 1;
DO i = 2, n
  pivot = L(i)/Diag(j)
  Diag(i) = Diag(i) - pivot * U(j)
  B(i) = B(i) - pivot * b(j)  
  j = i
ENDDO
j = n
x(n) = b(n)/ Diag(n)
DO i = n-1, 1, -1
  x(i) = (b(i) - U(i) * x(j)) / Diag(i)
  j = i
ENDDO
END SUBROUTINE

subroutine LU_factorize(A,LU,pvec,n)
implicit none
real(DB),intent(in) :: A(n,n)
real(DB),intent(out) :: LU(n,n)
integer,intent(Out) :: pvec(n)
integer,intent(in) :: n
integer :: i,j,k,iipvt(1),ipvt,itp
real(DB) :: pvt,rowtp(n),lmnt
equivalence(ipvt,iipvt(1))
LU=A
DO i =1,n
  pvec(i) = i
ENDDO
DO i = 1,n-1
  pvt = LU(i,i)
  iipvt = MAXLOC(abs(LU(i:n,i))) + i-1
  pvt = LU(ipvt,i)
  itp = pvec(i); pvec(i) = ipvt; pvec(ipvt)=itp
  rowtp(:) = LU(i,:);LU(i,:) = LU(ipvt,:)
  LU(ipvt,:) = rowtp(:)
  continue
  pvt=1._DB/pvt
  DO j = i+1,n
    lmnt = pvt*LU(j,i)
    LU(j,i) = lmnt
    LU(j,i+1:n)=LU(j,i+1:n) - lmnt*(LU(i,i+1:n))
  ENDDO
  continue
ENDDO
DO i=1,n
  LU(i,i) = 1._DB/LU(i,i)
ENDDO
continue
end subroutine      

subroutine LU_solver(LU,pvec,b,sol,n)
implicit none
REAL(DB),intent(in) :: LU(n,n),b(n)
integer,intent(in) :: pvec(n),n
REAL(DB),intent(out) :: sol(n)
REAL(DB) :: y(n),lmnt
integer :: i,j
!Forward
sol = 0.0_DB
DO i=1,n
  lmnt = 0.0_DB
  DO j=1,i-1,1
    lmnt = lmnt - LU(i,j)*y(j)
  ENDDO 
  y(i) = lmnt + b(pvec(i))
ENDDO

DO i=n,1,-1
  lmnt = 0.0_DB
  DO j= i+1,n,1
    lmnt = lmnt - LU(i,j)*sol(j)
  ENDDO
  sol(i) = (lmnt + y(i))*LU(i,i)
ENDDO
end subroutine

END MODULE