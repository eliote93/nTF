MODULE UtilFunction
USE PARAM
IMPLICIT NONE

PUBLIC
PRIVATE :: ObjectiveF, NewLineDcmp, GetIntRand

INTERFACE MDArrayTo1D
  MODULE PROCEDURE MDArrayTo1D_2DR
  MODULE PROCEDURE MDArrayTo1D_2DI
  MODULE PROCEDURE MDArrayTo1D_3DR
  MODULE PROCEDURE MDArrayTo1D_3DI
  MODULE PROCEDURE MDArrayTo1D_4DR
  MODULE PROCEDURE MDArrayTo1D_4DI
END INTERFACE

INTERFACE Array1DToMD
  MODULE PROCEDURE Array1DToMD_2DR
  MODULE PROCEDURE Array1DToMD_2DI
  MODULE PROCEDURE Array1DToMD_3DR
  MODULE PROCEDURE Array1DToMD_3DI
  MODULE PROCEDURE Array1DToMD_4DR
  MODULE PROCEDURE Array1DToMD_4DI
END INTERFACE

CONTAINS


SUBROUTINE ArrayLoc(idx, n, Dims, Locs)
IMPLICIT NONE
INTEGER :: n, idx
INTEGER :: Dims(n), LOCS(n)
INTEGER :: m, i, j, k, nn
m = 1
DO i = 1, n
  m = m * Dims(i)
ENDDO
nn = idx
DO i = n, 1, -1
  m = m / Dims(i)
  IF(i .EQ. 1) THEN
    LOCS(1) = nn
    EXIT
  ENDIF
  DO j = 1, DIMS(i)
    IF(nn .LE. j * m) EXIT
  ENDDO
  j = nn / m;   k = MOD(nn, m)
  IF(k .EQ. 0) THEN
    j = j
  ELSE
    j = j + 1
  ENDIF
  nn  = nn - m * (j-1)
  LOCS(i) = j
  IF(nn .eq. 0 .AND. i .NE. 1) THEN
    LOCS(1:i-1) = DIMS(1:i-1)
    EXIT
  ENDIF
ENDDO

END SUBROUTINE

SUBROUTINE Array1DToMD_2DR(x,y,n1,n2,n)
IMPLICIT NONE
REAL :: X(n1, n2)
REAL :: Y(n)
INTEGER :: n1, n2, n
INTEGER :: i, j, k, l, m
m = 0
DO j = 1, n2
  DO i = 1, n1
    m = m + 1
    Y(m) = X(i, j)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE Array1DToMD_2DI(x,y,n1,n2,n)
IMPLICIT NONE
INTEGER :: X(n1, n2)
INTEGER :: Y(n)
INTEGER :: n1, n2, n
INTEGER :: i, j, k, l, m
m = 0
DO j = 1, n2
  DO i = 1, n1
    m = m + 1
    Y(m) = X(i, j)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE Array1DToMD_3DR(x,y,n1,n2,n3,n)
IMPLICIT NONE
REAL :: X(n1, n2,n3)
REAL :: Y(n)
INTEGER :: n1, n2, n3, n
INTEGER :: i, j, k, l, m
m = 0
DO k = 1, n3
  DO j = 1, n2
    DO i = 1, n1
      m = m + 1
      Y(m) = X(i, j, k)
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE Array1DToMD_3DI(x,y,n1,n2,n3,n)
IMPLICIT NONE
INTEGER :: X(n1, n2,n3)
INTEGER :: Y(n)
INTEGER :: n1, n2, n3, n
INTEGER :: i, j, k, l, m
m = 0
DO k = 1, n3
  DO j = 1, n2
    DO i = 1, n1
      m = m + 1
      Y(m) = X(i, j, k)
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE Array1DToMD_4DR(x,y,n1,n2,n3,n4,n)
IMPLICIT NONE
REAL :: X(n1, n2, n3, n4)
REAL :: Y(n)
INTEGER :: n1, n2, n3, n4, n
INTEGER :: i, j, k, l, m
m = 0
DO l = 1, n4
  DO k = 1, n3
    DO j = 1, n2
      DO i = 1, n1
        m = m + 1
        Y(m) = X(i, j, k, l)
      ENDDO
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE Array1DToMD_4DI(x,y,n1,n2,n3,n4,n)
IMPLICIT NONE
INTEGER :: X(n1, n2, n3, n4)
INTEGER :: Y(n)
INTEGER :: n1, n2, n3, n4, n
INTEGER :: i, j, k, l, m
m = 0
DO l = 1, n4
  DO k = 1, n3
    DO j = 1, n2
      DO i = 1, n1
        m = m + 1
        Y(m) = X(i, j, k, l)
      ENDDO
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE


SUBROUTINE MDArrayTo1D_2DR(x,y,n1,n2,n)
IMPLICIT NONE
REAL :: X(n1, n2)
REAL :: Y(n)
INTEGER :: n1, n2, n
INTEGER :: i, j, k, l, m
m = 0
DO j = 1, n2
  DO i = 1, n1
    m = m + 1
    Y(m) = X(i, j)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE MDArrayTo1D_2DI(x,y,n1,n2,n)
IMPLICIT NONE
INTEGER :: X(n1, n2)
INTEGER :: Y(n)
INTEGER :: n1, n2, n
INTEGER :: i, j, k, l, m
m = 0
DO j = 1, n2
  DO i = 1, n1
    m = m + 1
    Y(m) = X(i, j)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE MDArrayTo1D_3DR(x,y,n1,n2,n3,n)
IMPLICIT NONE
REAL :: X(n1, n2,n3)
REAL :: Y(n)
INTEGER :: n1, n2, n3, n
INTEGER :: i, j, k, l, m
m = 0
DO k = 1, n3
  DO j = 1, n2
    DO i = 1, n1
      m = m + 1
      Y(m) = X(i, j, k)
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE MDArrayTo1D_3DI(x,y,n1,n2,n3,n)
IMPLICIT NONE
INTEGER :: X(n1, n2,n3)
INTEGER :: Y(n)
INTEGER :: n1, n2, n3, n
INTEGER :: i, j, k, l, m
m = 0
DO k = 1, n3
  DO j = 1, n2
    DO i = 1, n1
      m = m + 1
      Y(m) = X(i, j, k)
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE MDArrayTo1D_4DR(x,y,n1,n2,n3,n4,n)
IMPLICIT NONE
REAL :: X(n1, n2, n3, n4)
REAL :: Y(n)
INTEGER :: n1, n2, n3, n4, n
INTEGER :: i, j, k, l, m
m = 0
DO l = 1, n4
  DO k = 1, n3
    DO j = 1, n2
      DO i = 1, n1
        m = m + 1
        Y(m) = X(i, j, k, l)
      ENDDO
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE MDArrayTo1D_4DI(x,y,n1,n2,n3,n4,n)
IMPLICIT NONE
INTEGER :: X(n1, n2, n3, n4)
INTEGER :: Y(n)
INTEGER :: n1, n2, n3, n4, n
INTEGER :: i, j, k, l, m
m = 0
DO l = 1, n4
  DO k = 1, n3
    DO j = 1, n2
      DO i = 1, n1
        m = m + 1
        Y(m) = X(i, j, k, l)
      ENDDO
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Array2DSORT(X, Y, n, npt, lAcending, lRmvRpt, iPriority)
REAL :: X(n), Y(n)
INTEGER :: n                   !array Size
INTEGER :: npt
INTEGER :: iPriority           !1 : X sorting, 2: Y sorting
LOGICAL :: lAcending           !Ascending Order Flag
LOGICAL :: lRmvRpt             !Remove Repetition Points
REAL, POINTER :: A(:), B(:)
!REAL :: A(n), B(n)

REAL :: temp
INTEGER :: i, j
ALLOCATE(A(n), B(n))
IF(iPriority .EQ. 1) THEN
  A = X; B = Y
ELSE
  A = Y; B = X;
ENDIF

IF(lAcending) THEN
  DO j = n, 1, -1
    DO i = 1, j - 1
      IF(A(i) > A(i+1)) THEN
        temp = A(i)
        A(i) = A(i+1); A(i+1) = temp
        temp = B(i)
        B(i) = B(i+1); B(i+1) = temp
      ENDIF
    ENDDO
  ENDDO
ELSE
  DO j = n, 1, -1
    DO i = 1, j - 1
      IF(A(i) < A(i+1)) THEN
        temp = A(i)
        A(i) = A(i+1); A(i+1) = temp
        temp = B(i)
        B(i) = B(i+1); B(i+1) = temp
      ENDIF
    ENDDO
  ENDDO
ENDIF

IF(lRmvRpt) THEN
  CALL Array2DRmvRpt(A, B, n, npt, 1)
  !CALL Array2DRmvRpt(A, B, n, npt, 2)
ELSE
  npt = n
ENDIF

IF(iPriority .EQ. 1) THEN
  X = A; Y = B
ELSE
  X = B; Y = A
ENDIF
DEALLOCATE(A, B)
END SUBROUTINE

SUBROUTINE Array1DSORT(X, n, npt, lAcending, lRmvRpt)
REAL :: X(n)
INTEGER :: n                   !array Size
INTEGER :: npt
LOGICAL :: lAcending           !Ascending Order Flag
LOGICAL :: lRmvRpt             !Remove Repetition Points
REAL, POINTER :: A(:), B(:)

REAL :: temp
INTEGER :: i, j
ALLOCATE(A(n))
A = X

IF(lAcending) THEN
  DO j = n, 1, -1
    DO i = 1, j - 1
      IF(A(i) > A(i+1)) THEN
        temp = A(i)
        A(i) = A(i+1); A(i+1) = temp
      ENDIF
    ENDDO
  ENDDO
ELSE
  DO j = n, 1, -1
    DO i = 1, j - 1
      IF(A(i) < A(i+1)) THEN
        temp = A(i)
        A(i) = A(i+1); A(i+1) = temp
      ENDIF
    ENDDO
  ENDDO
ENDIF

IF(lRmvRpt) THEN
  CALL Array1DRmvRpt(A, n, npt)
ELSE
  npt = n
ENDIF
X = A
DEALLOCATE(A)
END SUBROUTINE

SUBROUTINE Array2DRmvRpt(X, Y, n, npt,  iPriority)
IMPLICIT NONE
REAL :: X(n), Y(n)
INTEGER :: n, npt, iPriority
!REAL :: tA(n), tB(n), A(n), B(n)
REAL, POINTER :: tA(:), tB(:), A(:), B(:)
INTEGER :: i

ALLOCATE(tA(n), tB(n), A(n), B(n))
tA = X; tB = Y
IF(iPriority .EQ. 2) THEN
  tA = Y; tB = X
ENDIF

A = 0; B = 0
npt = 1
A(1) = tA(1); B(1) = tB(1)
DO i = 2, n
  IF(abs(A(npt) - tA(i)) > epsm5) THEN
    npt = npt + 1
    A(npt) = tA(i); B(npt) = tB(i)
  ENDIF
ENDDO

X = A; Y = B
IF(iPriority .EQ. 2) THEN
  Y = A; X = B
ENDIF

DEALLOCATE(tA, tB, A, B)
END SUBROUTINE

SUBROUTINE Array1DRmvRpt(X, n, npt)
IMPLICIT NONE
REAL :: X(n)
INTEGER :: n, npt
REAL :: tA(n), tB(n), A(n), B(n)

INTEGER :: i

tA = X

A = 0;
npt = 1; A(1) = tA(1)
DO i = 2, n
  IF(abs(tA(i-1) - tA(i)) > epsm10) THEN
    npt = npt + 1
    A(npt) = tA(i)
  ENDIF
ENDDO
X = A
END SUBROUTINE

FUNCTION Strlen(a, nsize)
IMPLICIT NONE
INTEGER :: nsize
CHARACTER(nsize) :: a
INTEGER :: StrLen
INTEGER :: i
DO i = 1, nsize
  IF(a(i:i) .eq. '') EXIT
ENDDO
IF(i .GT. nsize) i = nsize
Strlen = i
END FUNCTION

SUBROUTINE CopyP1AxFlx(Flx1, Flx2, lTran)
USE PARAM
USE TYPEDEF, ONLY : AxFlx_Type
IMPLICIT NONE
TYPE(AxFlx_Type) :: Flx1, Flx2
LOGICAL :: lTran
Flx1%Phi = Flx2%Phi; Flx1%Psi = Flx2%Psi
Flx1%Tlkg = Flx2%Tlkg; Flx1%Jout = Flx2%Jout
Flx1%LkgShape = Flx2%LkgShape;
Flx1%Dhat = Flx2%Dhat; Flx1%Dtil = Flx2%Dtil
Flx1%PDhat = Flx2%PDhat
IF(lTran) THEN
  Flx1%TranSrc = Flx2%TranSrc
ENDIF

END SUBROUTINE

SUBROUTINE CopySP3AxFlx(Flx1, Flx2, lTran)
USE PARAM
USE TYPEDEF, ONLY : AxFlx_Type
IMPLICIT NONE
TYPE(AxFlx_Type) :: Flx1, Flx2
LOGICAL :: lTran
Flx1%Phi = Flx2%Phi; Flx1%Psi = Flx2%Psi
Flx1%Tlkg = Flx2%Tlkg; Flx1%Jout = Flx2%Jout
Flx1%LkgShape = Flx2%LkgShape;

Flx1%S = Flx2%S; Flx1%KSQ = Flx2%KSQ
Flx1%QT = Flx2%QT;

Flx1%Dhat = Flx2%Dhat; Flx1%Dtil = Flx2%Dtil
Flx1%PDhat = Flx2%PDhat

ENDSUBROUTINE

SUBROUTINE CopyAxXs(XS1, XS2, ng)
USE PARAM
USE TYPEDEF, ONLY : PinXS_TYPE
TYPE(PinXS_TYPE) :: XS1, XS2
INTEGER :: ng
INTEGER :: ig
XS1%Dtil = XS2%Dtil;  XS1%Dhat = XS2%Dhat
XS1%XSD = XS2%XSD;    XS1%XSD2 = XS2%XSD2
XS1%XST = XS2%XST;    XS1%XSTR = XS2%XSTR
XS1%XSR = XS2%XSR;    XS1%XSNF = XS2%XSNF
XS1%XSKF = XS2%XSKF;  XS1%CHI = XS2%CHI
DO ig = 1, ng
  XS1%xss(ig)%WithInGroupScat = XS2%xss(ig)%WithInGroupScat
  XS1%xss(ig)%from = XS2%xss(ig)%from
ENDDO
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine segev(jt1,rsinp,youtm,xdat,ydat)
!     --------------------------------------------------------------------------------
!     Function   : 'SEGEV' interpolation for the resonance integral
!     Institute  : Oak Ridge National Laboratory
!     Authors    : Kang-Seog Kim, Mark L. Williams
!     Revision   : 02/15/2009  Initial version
!     --------------------------------------------------------------------------------
USE PARAM
IMPLICIT NONE
!
INTEGER :: jt1
REAL :: rsinp, youtm
REAL xdat(20),ydat(20), g1, g2
REAL :: err, epm9, epm3, pval, eta, dif1, dif2, rho, w1, w2, w3, gw
INTEGER :: i, max, min, krt
data err/1.0e-5/,epm9/1.0e-9/,epm3/1.0e-3/
!
!     --------------------------------------------------------------------------------
!     @ jt1                    (r)    the number of sigma-b                          +
!     @ xdat(i)                (r)    background cross sections (i=1,jt1)            +
!     @ ydat(i)                (r)    ri/du (i=1,jt1)                                +
!     --------------------------------------------------------------------------------
!
!     [when the target sigma-b is less than sigma-b in the table]
if (rsinp.le.xdat(1)) then
  pval=1.0
  eta=xdat(1)*(ydat(jt1)-ydat(1))/ydat(1)
  youtm=ydat(jt1)*(rsinp/(rsinp+eta))**pval
!
!     [when the target sigma-b is larger than sigma-b in the table]
elseif (rsinp.ge.xdat(jt1)) then
  youtm=ydat(jt1)
!
!     [when the target sigma-b is between two entries]
else
!       +decide the position of 'rsinp'+
  do i=1,jt1-1
    dif1=rsinp-xdat(i)
    dif2=rsinp-xdat(i+1)
    if (dif1*dif2.le.0.0) then
      min=i
      max=i+1
      exit
    endif
  enddo
!       +when the target sigma-b is between the maximum and the next+
  if (max.eq.jt1) then
    pval=1.0
    eta=xdat(min)*(ydat(jt1)-ydat(min))/ydat(min)
    youtm=ydat(jt1)*(rsinp/(rsinp+eta))**pval
!
!       +when the target sigma-b is between two values+
  else
    rho=xdat(min)/xdat(max)
    g1=ydat(jt1)/ydat(min)
    g2=ydat(jt1)/ydat(max)
!         this resolve the situation where ri at any background xs is larger
!         than ri(infinite) or too close to ri(infinite).
    if (dabs(ydat(min)/ydat(max)-1.0d0).le.epm3) then
       youtm=(ydat(min)+ydat(max))/2.0d0
    elseif ((ydat(min).ge.ydat(jt1)).or.(ydat(max).ge.ydat(jt1))) then
       youtm=(ydat(min)+ydat(max))/2.0d0
    else
!           +data check+
!      if (g1.le.0.0.or.g2.le.0.0.or.rho.le.0.0.or.dlog(g1)/dlog(g2).le.0.0) then
!!YS MODI - 03/19/2009
!         call decart_terminate(io8,'Error: Cannot perform SEGEV interpolation')
!!              call error(" ",8)
!!YS MODI END
!      endif
!
      w1=dlog(rho*dlog(g1)/dlog(g2))/dlog(g2/g1)
!           +when 'w' is negative & positive+
      if (w1.lt.0.0) then
        !if (rho.ge.1.0) call decart_terminate(io8,'Error: Cannot perform SEGEV interpolation')
        w2=dlog(1.0d0-rho)/dlog(g2)
!             +when 'w' is positive+
      elseif (w1.gt.0.0) then
        w2=dlog(rho)/dlog(g2/g1)
      endif
!           +search for the solution for 'g(w)'+
      krt=0
      do i=1,500
        if (krt.ne.0) exit
        w3=w1+(w2-w1)/2.0d0
        gw=g2**w3-rho*g1**w3-(1.0d0-rho)
        if (dabs(gw).lt.err) then
          krt=1
        else
          if (gw.gt.0.0) then
            w1=w3
          elseif (gw.lt.0.0) then
            w2=w3
          endif
        endif
      enddo
      pval=1.0d0/w3
      eta=xdat(min)*(g1**w3-1.0d0)
      youtm=ydat(jt1)*(rsinp/(rsinp+eta))**pval
    endif
  endif
endif


return
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GetRangeDecomp(n1, n2, nproc, myrank, ibeg, iend)
IMPLICIT NONE
INTEGER :: n1, n2, nproc, myrank, ibeg, iend
INTEGER :: iwork1, iwork2
iwork1 = ( n2 - n1 + 1 ) / nproc
iwork2 = MOD( n2 - n1 + 1 , nproc)
ibeg = myrank * iwork1 + n1 + MIN( myrank , iwork2 )
iend = ibeg + iwork1 - 1
IF(iwork2 .gt. myrank) iend=iend+1
IF(ibeg .gt. iend) iend=ibeg
END SUBROUTINE

!!!!!!!!!!!!!!!!!!Radial Domain Decomposition with Simulated Annealing!!!!!!!!!!!!!!!!!!!

SUBROUTINE GetRadDecomp_SA(Dcmp, LineInfo, nline, nproc)
INTEGER :: nline, nproc
INTEGER :: Dcmp(2, nline), LineInfo(nline)


INTEGER :: NewDcmp(nproc), RefDcmp(nproc), BaseDcmp(nproc), BestDcmp(nproc)
INTEGER :: i, j, k, ibeg, iend, m1, m2, m
INTEGER :: nstag, ntemp

REAL :: avg, std, std0, stdbase, stdbest
REAL :: C, Cbase, jpavg, p, p0

ntemp = 50; nstag = 150
CALL RANDOM_SEED(); CALL RANDOM_NUMBER(p0)
DO i = 1, nproc
  RefDcmp(i) = Dcmp(2, i) - Dcmp(1, i) + 1
ENDDO
avg = sum(LineInfo(:))/nproc
std0 = ObjectiveF(RefDcmp, LineInfo, nLine, nproc, avg)
stdbest = std0; StdBase = std0
BestDcmp = RefDcmp; BaseDcmp=RefDcmp
m = 50; m2 = 0
NewDcmp = RefDcmp
DO i = 1, m
  NewDcmp = NewLineDcmp(RefDcmp, nline, nproc)
  std = ObjectiveF(NewDcmp, LineInfo, nLine, nproc, avg)
  IF(std .GT. stdbase) THEN
    m2 = m2 + 1
    jpavg = jpavg + (std-stdbase)
  ENDIF
ENDDO
m1 = m - m2; jpavg = jpavg / m2
Cbase = jpavg / log(m2/(m2*0.95-m1*(1-0.95)));
C = Cbase

DO i = 1, ntemp
  DO j = 1, nstag
    NewDcmp = NewLineDcmp(BaseDcmp, nline, nproc)
    std = ObjectiveF(NewDcmp, LineInfo, nLine, nproc, avg)
    IF(std .LT. stdbase) THEN
      p=1;
    ELSE
      p = exp(-(std-stdbase)/C)
    ENDIF
    CALL RANDOM_NUMBER(p0)
    IF(p0 .LT. p) THEN
      BaseDcmp = NewDcmp;
      stdbase = std
      IF(std .LT. stdbest) THEN
        StdBest = Std; BestDcmp = BaseDcmp
      ENDIF
    ENDIF
  ENDDO
  !Temperature
  C=C*0.85_8
ENDDO
iend = 0
DO i = 1, nproc
  ibeg = iend + 1; iend = ibeg + BestDcmp(i) - 1
  Dcmp(1, i) = ibeg; Dcmp(2, i) = iend;
ENDDO
END SUBROUTINE

FUNCTION ObjectiveF(LineDcmp, LineInfo, nLine, nproc, avg)
IMPLICIT NONE
INTEGER :: nLine, nProc
INTEGER :: LineInfo(nLine), LineDcmp(nproc)
INTEGER :: ndat(nproc)
INTEGER :: ibeg, iend, i
REAL :: ObjectiveF, avg
ObjectiveF = 0; iend = 0
DO i=1, nproc
  ibeg = iend + 1
  iend = ibeg + LineDcmp(i) - 1
  nDat(i) = sum(LineInfo(ibeg:iend))
  ObjectiveF = ObjectiveF + (nDat(i)-avg)**2
ENDDO

ObjectiveF =sqrt(ObjectiveF/nproc)
END FUNCTION

FUNCTION NewLineDcmp(OldLineDcmp, nline, nproc)
IMPLICIT NONE
INTEGER :: nline, nproc
INTEGER :: OldLineDcmp(nproc), NewLineDcmp(nproc)
INTEGER :: m1, m2, m3, j
1000 CONTINUE
m1 = GetIntRand(nproc)
j=0
DO
  j = j + 1
  IF(j .GT. 20) GOTO 1000
  m2 = GetIntRand(nproc)
  IF(m2 .EQ. m1) CYCLE
  IF(OldLineDcmp(m2) .EQ. 1) CYCLE
  EXIT
ENDDO
m3 = GetIntRand(OldLineDcmp(m2)-1)
NewLineDcmp = OldLineDcmp
NewLineDcmp(m1) = NewLineDcmp(m1) + m3
NewLineDcmp(m2) = NewLineDcmp(m2) - m3
END FUNCTION


FUNCTION GetIntRand(nmax)
INTEGER:: GetIntRand, nmax
REAL :: tmp
CALL RANDOM_NUMBER(tmp)
tmp = tmp * nmax + 1
GetIntRand = int(tmp)
END FUNCTION

FUNCTION LineIntPol(x, ndat, xdat, ydat)
USE PARAM
IMPLICIT NONE
REAL :: x, LineIntPol
INTEGER :: ndat
REAL :: xdat(ndat), ydat(ndat)
REAL :: wt1, wt2
INTEGER :: i, n1, n2

IF(ndat .EQ. 1) THEN
  LineIntPol = ydat(1); RETURN
ENDIF

DO i = 2, ndat
  if(x .le. xdat(i)) EXIT
ENDDO
IF(i .GT. ndat) i = ndat
n2 = i; n1 = n2 - 1
wt2 = (x - xdat(n1))/(xdat(n2) - xdat(n1))
wt1 = 1 - wt2
LineIntPol = wt2 * ydat(n2) + wt1 * ydat(n1)
END FUNCTION

!FUNCTION PolyIntegral(a, n, x1, x2)
!USE PARAM
!IMPLICIT NONE
!INTEGER :: n
!REAL :: a(0:n), x1, x2
!REAL :: PolyIntegral
!
!REAL :: tmp, sol
!INTEGER :: i
!
!PolyIntegral = 0
!DO i = 0, n
!  tmp = (x2**(n+1) - x1**(n+1))/REAL(n+1, 8)
!  tmp = a(i) * tmp / REAL(n+1, 8)
!  PolyIntegral = PolyIntegral + tmp
!ENDDO
!
!END FUNCTION

FUNCTION Lp4thIntegral(C, x1, x2)
USE PARAM
IMPLICIT NONE
REAL :: Lp4thIntegral
REAL :: C(0:4), x1, x2
Lp4thIntegral = IndefLpIntegral(C, x1) - IndefLpIntegral(C, x2)
CONTAINS

FUNCTION IndefLpIntegral(C, x)
IMPLICIT NONE
REAL :: IndefLpIntegral
REAL :: A, B, C(0:4), x
REAL :: y

y = x * c(0);
y = y + c(1) * (0.5_8 * x *x)
y = y + c(2) * 0.5_8 * (-x + x**3)
y = y + c(3) * (-0.75_8 * (x**2) + 0.625 * (x**4))
y = y + c(4) * (0.375 * x - 1.25*(x**3) + 0.875 * (x**5))
IndefLpIntegral = y
END FUNCTION
END FUNCTION

SUBROUTINE CombineIntList(List, n, List1, n1, List2, n2, ndim)
USE BasicOperation,         ONLY : CP_CA
IMPLICIT NONE
INTEGER :: n, n1, n2, ndim
INTEGER :: List(ndim), List1(ndim), List2(ndim)

INTEGER :: i, j, k
INTEGER :: m, m1, m2
INTEGER :: v1, v2
m = 0; m1 = 0; m2 = 0
CALL CP_CA(List, 0, ndim)
v1 = 1000000000; v2 = 1000000000
IF(n1 .NE. 0) v1 = List1(1)
IF(n2 .NE. 0) v2 = List2(1)
DO
  IF((m1 .EQ. n1) .AND. (m2 .EQ. n2)) EXIT
  IF(v1 .EQ. v2) THEN
    m = m + 1; m1= m1 + 1; m2 = m2 + 1
    List(m) = v1
    v1 = 1000000000; v2 = 1000000000
    IF(m1 .LT. n1) v1 = List1(m1 +1)
    IF(m2 .LT. n2) v2 = List2(m2 + 1)
  ELSEIF(v1 .LT. v2) THEN
    m = m + 1; m1 = m1 + 1
    List(m) = v1
    v1 = 1000000000
    IF(m1 .LT. n1) v1 = List1(m1 +1)
  ELSE
    m = m + 1; m2 = m2 + 1
    List(m) = v2
    v2 = 1000000000
    IF(m2 .LT. n2) v2 = List2(m2 +1)
  ENDIF
ENDDO
n = m
END SUBROUTINE

END MODULE
