MODULE GeomTreatment 
!Set of procedure for finding interception Points
USE PARAM
USE typedef, ONLY : BasicGeom
USE IOUTIL, ONLY : Terminate
IMPLICIT NONE

CONTAINS

SUBROUTINE Cell_RayIntersection(GeomInfo, ray, pts, npt, tor)
USE UtilFunction, ONLY : Array2DSORT
!
IMPLICIT NONE
TYPE(BasicGeom), INTENT(IN) :: GeomInfo
REAL :: RAY(3)
REAL :: pts(2, 100)                      !Intersection Points
INTEGER :: npt
REAL :: tor
!
REAL :: xlim(2), ylim(2)                 !Pin Domain
REAL :: Lines(3, 500)                    !Lines which are boundarys of FSR
REAL :: SOL(2, 100)
REAL :: Circle(3)                  !Circles which are boundarays of FSR
REAL :: localsol(2,2), x, y
REAL :: CRITX, CRITY
INTEGER :: nline, ncircle, npt0
INTEGER :: i,j, n
!EQUIVALENCE(Localsol(1), x)
!EQUIVALENCE(Localsol(2), y)
!Rectangular Cell Domain
xlim(1) = GeomInfo%x(1); xlim(2) = GeomInfo%x(2); 
ylim(1) = GeomInfo%y(1); ylim(2) = GeomInfo%y(2); 

nline = GeomInfo%nline
ncircle = GeomInfo%nCircle
lines(:, 1:nline) = GeomInfo%line(:, 1:nline)
npt = 0
DO i = 1, nline
  localsol(:,1) = Line_Line_Intersection(lines(:, i), RAY)
  !interection points adaptiveness check
  x= localsol(1,1); y=localsol(2,1)
  critx=(x-xlim(1))*(x-xlim(2))
  crity=(y-ylim(1))*(y-ylim(2))
  !IF((x-xlim(1))*(x-xlim(2)) .LT. tor) THEN
  IF((x-xlim(1)) .GT. tor .AND. (x-xlim(2)) .LT. -tor) THEN
    IF((y-ylim(1)) .LT. tor .OR. (y-ylim(2)) .GT. -tor) CYCLE
    !IF((y-ylim(1))*(y-ylim(2)) .GT. tor) CYCLE
    npt = npt + 1  
    sol(:,npt) = LocalSol(:,1)
  ENDIF
ENDDO

IF(GeomInfo%lcircle) THEN
  ncircle = GeomInfo%nCircle
  DO i = 1, nCircle
    Circle = GeomInfo%circle(:, i)
    localsol = Line_Circle_Intersection(RAY, circle, n)
    IF(n .GE. 1) THEN
      DO j = 1, n
        IF((localsol(1, j)-xlim(1))*(localsol(1, j)-xlim(2)) .LT. epsm6) THEN
          IF((localsol(2, j)-ylim(1))*(localsol(2, j)-ylim(2)) .GT. epsm6) CYCLE
          npt = npt + 1  
          sol(:,npt) = LocalSol(:,j)
        ENDIF       
      ENDDO
      !SOL(:, (npt + 1):(npt + n)) = LocalSol(:, 1:npt)
      !npt = npt + n
    ENDIF
  ENDDO
ENDIF

IF(npt .NE. 0) THEN
  !Solution Sorting and 
  CALL Array2DSort(sol(1,1:npt), sol(2,1:npt), npt, n, TRUE, TRUE, 2)
  npt = n
  pts(:,1:n) = Sol(:,1:n)
ELSE
  pts = 0; npt = 0
ENDIF
END SUBROUTINE

SUBROUTINE GeomBoundary_RayIntersection(GeomInfo, InOut, InOutSurf, ray, npt)
USE UtilFunction, ONLY : Array2DSORT
IMPLICIT NONE
TYPE(BasicGeom), INTENT(IN) ::  GeomInfo
REAL, INTENT(IN) :: ray(3)
REAL, INTENT(OUT) :: InOut(2,2)
INTEGER, INTENT(OUT) :: InOutSurf(2)
INTEGER :: NbdLine, npt, n                      !Number of boundary line 
REAL :: BdLine(3, 20), SOL(2,20), LocalSol(2)
REAL :: xlim(2), ylim(2), x, y, temp, temp1, temp2
INTEGER :: i, j, k
EQUIVALENCE(Localsol(1), x)
EQUIVALENCE(Localsol(2), y)
NbdLine = GeomInfo%nbd                                    !Number of Boundary Line
BdLine(1:3, 1:NbdLine) = GeomInfo%Bdline(1:3, 1:NbdLine)
xlim = (/GeomInfo%x(1), GeomInfo%x(2)/)
ylim = (/GeomInfo%y(1), GeomInfo%y(2)/)
npt = 0
DO i = 1, NbdLine
  LocalSol = Line_Line_Intersection(BdLine(:,i), ray)
  !temp = (x-xlim(1))*(x-xlim(2))
  IF((x-xlim(1))*(x-xlim(2)) .LT. epsm6) THEN
    IF((y-ylim(1))*(y-ylim(2)) .GT. epsm6) CYCLE
    npt = npt + 1  
    Sol(:,npt) = LocalSol
  ENDIF
ENDDO

IF(npt .EQ. 0) RETURN

!CALL Array2DSort(Sol(1,1:npt), Sol(2,1:npt), npt, n, TRUE, FALSE, 2)
CALL Array2DSort(Sol(1,1:npt), Sol(2,1:npt), npt, n, TRUE, TRUE, 2)
npt = n
InOut = Sol(:,1:2)

!DO i = 1, 2
!  write(82, '(2e20.5)') Sol(1, i), Sol(2, i)
!ENDDO
!Determine In 
DO j = 1, 2
  DO i = 1, NbdLine
    IF(abs(Bdline(1, i)*Sol(1, j)+Bdline(2, i)*Sol(2, j)+Bdline(3, i)) .LT. epsm8) EXIT
  ENDDO
  InOutSurf(j) = i
ENDDO
!In & Out Surface 
END SUBROUTINE



!Basic FUNCTION UNIT
FUNCTION Line_Circle_Intersection(line, circle, n)
IMPLICIT NONE
REAL :: line(3), circle(3)
INTEGER, INTENT(INOUT) :: n
REAL :: Line_Circle_Intersection(2,2), Sol(2, 2), Sol0(2, 2)
REAL :: a, b, c, cx, cy, r, det

REAL :: det1,det2,det3,det4

LOGICAL :: la0


cx = circle(1); cy = circle(2); r = circle(3)
a = line(1); b = line(2); c = line(3)+ a*cx + b*cy

la0 = .FALSE.
IF(abs(a) .LT.  epsm5) THEN
  la0 = .TRUE.
  a = line(2); b = line(1); c = line(3)+ a*cx + b*cy
ENDIF

det = a**2*(-c**2 + (a**2 + b**2)*r**2)
n = 2
IF(abs(det) .lt. epsm5) THEN
 n = 1; det =0
ELSEIF(det .lt. 0) THEN
 n = 0;  Line_Circle_Intersection = 0
 return
ENDIF
det = SQRT(det)
Sol(1, 1) = -((a**2*c + b*Det)/(a*(a**2 + b**2))) + cx
Sol(2, 1) = (-(b*c) + Det)/(a**2 + b**2)  + cy
Sol(1, 2) = (-(a**2*c) + b*Det)/(a*(a**2 + b**2)) + cx
Sol(2, 2) = -((b*c + Det)/(a**2 + b**2)) + cy
IF(la0) THEN
  Sol0(2, 1) =  Sol(1, 1); Sol0(1, 1) = Sol(2, 1)
  Sol0(2, 2) =  Sol(1, 2); Sol0(1, 2) = Sol(2, 2)
  Sol = Sol0
ENDIF

Line_Circle_Intersection = SOL

END FUNCTION
!FUNCTION Length(x1,x2)
!
!END FUNCTION

FUNCTION FindCellRayBase(SearchTree, x, SurfNo)
USE TYPEDEF, ONLY : CellRayBase_SearchTree
USE IOUTIL, ONLY : Terminate
IMPLICIT NONE

TYPE(CellRayBase_SearchTree) :: SearchTree
REAL :: x(2)
INTEGER :: SurfNo
INTEGER :: FindCellRayBase

INTEGER :: Xint(2)
INTEGER :: i, j, jbeg, jend, k, SearchDir
INTEGER :: k1, k2
INTEGER :: jmax, jmin
INTEGER :: maxv, minv
INTEGER :: imod0, idiv
INTEGER, POINTER :: pts(:,:,:)
REAL :: MULT
LOGICAL :: lFlag


MULT = SearchTree%MULT  !Real -> Integer Factor
!REAL -> Integer Converting
xint(1) = INT(MULT*X(1)); xint(2) = INT(MULT*X(2))

jbeg = SearchTree%InSurfIdxSt(SurfNo)
jend = jbeg + SearchTree%nPtsInSruf(SurfNo) - 1
!jbeg = 1
!jend = SearchTree%nline
pts => SearchTree%pts
IF(SurfNo .EQ. SOUTH .OR. SurfNo .EQ. NORTH) THEN  !X-coordinates comparision
  IF(pts(1, jbeg, 1) .GT. pts(1, jend, 1)) THEN
    maxv = pts(1, jbeg, 1); minv = pts(1, jend, 1);
    jmax = jbeg; jmin = jend
  ELSE
    maxv = pts(1, jend, 1); minv = pts(1, jbeg, 1); 
    jmax = jend; jmin = jbeg
  ENDIF
  
  IF(xint(1) .GE. maxv) THEN
    k = jmax
  ELSEIF(xint(1) .LE. minv) THEN
    k = jmin
  ELSE
    !Point to Point Comparsion
    DO j = jbeg, jend -1
      !Flag = 1.*(xint(1) - pts(1, j, 1))*(xint(1) - pts(1, j+1, 1))
      lFlag = FALSE
      IF( (xint(1) - pts(1, j, 1)) .GE. 0 .AND. (xint(1) - pts(1, j+1, 1)) .LT. 0) lFlag = TRUE
      IF( (xint(1) - pts(1, j, 1)) .LE. 0 .AND. (xint(1) - pts(1, j+1, 1)) .GT. 0) lFlag = TRUE
      IF(lFlag) THEN
        IF( abs(xint(1) - pts(1, j, 1)) .LT. abs(pts(1, j+1, 1) - xint(1)) ) THEN
          k = j
        ELSE
          k = j +1
        ENDIF
        EXIT
      ENDIF
    ENDDO
    IF(J .GT. jend - 1) THEN
      CALL Terminate('CAN NOT FIND Cell Ray Information')
    ENDIF    
  ENDIF
ELSE  !Y-coordinate Comparision
  IF(pts(2, jbeg, 1) .GT. pts(2, jend, 1)) THEN
    maxv = pts(2, jbeg, 1); minv = pts(2, jend, 1);
    jmax = jbeg; jmin = jend
  ELSE
    maxv = pts(2, jend, 1); minv = pts(2, jbeg, 1); 
    jmax = jend; jmin = jbeg
  ENDIF

  IF(xint(2) .GE. maxv) THEN
    k = jmax
  ELSEIF(xint(2) .LE. minv) THEN
    k = jmin
  ELSE
    DO j = jbeg, jend -1
      lFlag = FALSE
      IF((xint(2) - pts(2, j, 1)) .GE. 0 .and. (xint(2) - pts(2, j+1, 1)) .LT.0) lFlag = TRUE
      IF((xint(2) - pts(2, j, 1)) .LE. 0 .and. (xint(2) - pts(2, j+1, 1)) .GT. 0) lFlag = TRUE
      IF(lFlag) THEN
        IF( abs(xint(2) - pts(2, j, 1)) .LT. abs(pts(2, j+1, 1) - xint(2)) ) THEN
          k = j
        ELSE
          k = j +1
        ENDIF
        EXIT
      ENDIF
    ENDDO
    
    IF(J .GT. jend - 1) THEN
      CALL Terminate('CAN NOT FIND Cell Ray Information')
    ENDIF    
  ENDIF 
ENDIF
FindCellRayBase = k

END FUNCTION

FUNCTION RectAsyPinIdx(AsyGeom, x)
USE TYPEDEF, ONLY : BasicGeom
IMPLICIT NONE
TYPE(BasicGeom), INTENT(IN) :: AsyGeom
REAL :: x(2)
INTEGER :: RectAsyPinIdx
INTEGER :: nx, ny, ix, iy, ireg
REAL :: hx(0:1000), hy(0:1000)

nx = AsyGeom%nx; ny = AsyGeom%ny
hx(0) = AsyGeom%x(1); hx(nx) = AsyGeom%x(2)
hx(1:nx-1) = -AsyGeom%Line(3,1:nx-1)
hy(0) = AsyGeom%y(1); hy(ny) = AsyGeom%y(2)
hy(1:ny-1) = -AsyGeom%Line(3,nx:nx+ny-2)

DO ix = 1, nx !Serching X Cell Coordinate
  IF((hx(ix)-x(1))*(hx(ix-1)-x(1)) .LT. 0.) EXIT
ENDDO

DO iy = 1, ny 
  IF((hy(iy)-x(2))*(hy(iy-1)-x(2)) .LT. 0.) EXIT
ENDDO
iy = ny - iy + 1
ireg = nx * (iy - 1) + ix
RectAsyPinIdx = ireg 
END FUNCTION

FUNCTION RectSquareCellFsrIdx(CellInfo, x)
USE TYPEDEF, ONLY : Cell_Type
IMPLICIT NONE
TYPE(Cell_type), INTENT(IN) :: CellInfo
REAL :: x(2)
INTEGER :: RectSquareCellFsrIdx
INTEGER :: nx, ny, ix, iy, ireg
REAL :: hx(0:1000), hy(0:1000)

nx = CellInfo%geom%nx; ny = CellInfo%geom%ny
hx(0) = Cellinfo%geom%x(1); hx(nx) = Cellinfo%geom%x(2)
hx(1:nx-1) = -CellInfo%Geom%Line(3,1:nx-1)
hy(0) = Cellinfo%geom%y(1); hy(ny) = Cellinfo%geom%y(2)
hy(1:ny-1) = -CellInfo%Geom%Line(3,nx:nx+ny-2)

DO ix = 1, nx !Serching X Cell Coordinate
  IF((hx(ix)-x(1))*(hx(ix-1)-x(1)) .LT. epsm6) EXIT
ENDDO

DO iy = 1, ny 
  IF((hy(iy)-x(2))*(hy(iy-1)-x(2)) .LT. epsm6) EXIT
ENDDO
iy = ny - iy + 1
!y = iy;
ireg = nx * (iy - 1) + ix 
RectSquareCellFsrIdx = ireg
END FUNCTION

FUNCTION RectAnnularCellFsrIdx(CellInfo, x)
USE TYPEDEF, ONLY : Cell_Type
IMPLICIT NONE
TYPE(Cell_type), INTENT(IN) :: CellInfo
REAL :: x(2)
INTEGER :: RectAnnularCellFsrIdx

INTEGER :: CCentIdx(2, 4)
REAL :: sinv, cosv, radius, del_ang, ang
REAL :: cx, cy, x0, y0, r(0:100)
TYPE(BasicGeom), POINTER :: Geom
INTEGER :: nCircle, nDivAzi
INTEGER :: ir , iang, ireg

DATA CCentIdx / 1, 1,  1, 2,  2, 2,  2, 1 /

!Geom => CellInfo%Geom

IF(.NOT. CellInfo%Geom%lCCent) THEN
  R = 0
  nDivAzi = CellInfo%nDivAzi
  Del_ang = 2 * PI / nDivAzi
  nCircle = CellInfo%geom%nCircle;
  R(0) = CellInfo%Geom%lx*2
  R(1:nCircle) = CellInfo%Geom%Circle(3,:)
  !r(0:100)
  radius = SQRT(x(1) * x(1) + x(2) * x(2))
  !Ring Number
  DO ir = 1, nCircle
    IF((R(ir)-Radius)*(R(ir-1)-Radius) .LT. 0._8) EXIT
  ENDDO
  
  IF(x(1) .LT. 0.) THEN
    cosv = x(2) / radius
    ang = 0
  ELSE
    cosv = -x(2) / radius
    ang = PI
  ENDIF
  
  IF(CellInfo%lCentXY) THEN
    !ang = -PI / 2.
    ang = 0
    Del_ang = PI / nDivAzi/2
  ELSEIF(CellInfo%lCentX) THEN
    Del_ang = PI / nDivAzi
    IF(x(1) .LT. 0.) THEN
      ang = -PI / 2.
    ELSE  
      ang = PI / 2.
    ENDIF
  ELSEIF(CellInfo%lCentY) THEN
     Del_ang = PI / nDivAzi
     IF(x(1) .GT. 0.) ang = 0
  ENDIF
  ang = ang + ACOS(COSV)
  iang = INT(ang/Del_ang)
  ireg = nDivAzi*(ir-1) + iang + 1
  !CONTINUE
ELSE
  !CX = CellInfo%Geom%cx; CY = CellInfo%Geom%cy;
  CX = CellInfo%Geom%x(CCentIdx(1, CellInfo%Geom%CCentType))
  CY = CellInfo%Geom%y(CCentIdx(2, CellInfo%Geom%CCentType))
  x0 = X(1) - CX; y0 = X(2) - CY
  R = 0
  nDivAzi = CellInfo%nDivAzi
  Del_ang = PI / 4._8
  nCircle = CellInfo%geom%nCircle;
  R(0) = CellInfo%Geom%lx * 4
  R(1:nCircle) = CellInfo%Geom%Circle(3,:)
  radius = SQRT(x0 * x0 + y0 * y0)
  !Ring Number
  DO ir = 1, nCircle
    IF((R(ir)-Radius)*(R(ir-1)-Radius) .LT. 0._8) EXIT
  ENDDO
  SELECT CASE(CellInfo%Geom%CCentType)
    CASE(1)
      ang = 0
      cosv = x0 / radius
    CASE(2)
      ang = 0
      cosv = -y0 / radius      
    CASE(3)
      ang = 0
      cosv = -x0 / radius
    CASE(4)
      ang = 0
      cosv = y0 / radius      
  END SELECT
  ang = ang + ACOS(COSV)
  iang = INT(ang/Del_ang)
  ireg = nDivAzi*(ir-1) + iang + 1  
ENDIF
RectAnnularCellFsrIdx = iReg
END FUNCTION

FUNCTION Find_Surface(Geom, x)
USE TYPEDEF, ONLY: BasicGeom
IMPLICIT NONE

TYPE(BasicGeom) :: Geom
REAL :: x(2)
INTEGER :: Find_Surface
INTEGER :: NbdLine                       !Number of boundary line 
INTEGER :: I
REAL :: BdLine(3, 100)

NbdLine = Geom%nbd                                    !Number of Boundary Line
BdLine(1:3, 1:NbdLine) = Geom%Bdline(1:3, 1:NbdLine)
DO i = 1, NbdLine
  IF(abs(Bdline(1, i) * x(1) + Bdline(2, i)*x(2) + Bdline(3, i)) .LT. epsm6) EXIT
ENDDO
IF(i .GT. Nbdline) THEN
   CALL Terminate('Find_Surface : CAN NOT FIND SURFACE INDEX')
ENDIF
Find_Surface = i

END FUNCTION

FUNCTION Line_Line_Intersection(line1,line2)
IMPLICIT NONE
REAL :: Line_Line_Intersection(2)
REAL :: line1(3), line2(3)
REAL :: A(4), AINV(4), B(2), SOL(2) 
A(1) = line1(1); A(3) = line1(2); B(1) = -line1(3)
A(2) = line2(1); A(4) = line2(2); B(2) = -line2(3)
AINV = MatInv(A);
SOL = MatVecProduct(AINV, B)
Line_Line_Intersection= SOL
ENDFUNCTION

FUNCTION MatMatProduct(A,B)
IMPLICIT NONE
REAL :: MatMatProduct(4)
REAL :: A(4), B(4)

MatMatProduct(1) = A(1)*B(1)+A(3)*B(2)
MatMatProduct(2) = A(1)*B(3)+A(3)*B(4)
MatMatProduct(3) = A(2)*B(1)+A(4)*B(2)
MatMatProduct(4) = A(2)*B(3)+A(4)*B(4)

END FUNCTION

FUNCTION MatVecProduct(A,B)
IMPLICIT NONE
REAL :: MatVecProduct(2)
REAL :: A(4), B(2)

MatVecProduct(1) = A(1)*B(1)+A(3)*B(2)
MatVecProduct(2) = A(2)*B(1)+A(4)*B(2)

END FUNCTION

FUNCTION MatInv(A)
IMPLICIT NONE
REAL :: MatInv(4)
REAL :: A(4), DET

DET = A(1) * A(4) - A(2) * A(3)
MatInv(1) = A(4);MatInv(4) = A(1)
MatInv(2) = -A(2); MatInv(3) = -A(3) 
MatInv = MatInv / Det
END FUNCTION

!FUNCTION RayPassingCell(GeomInfo, Ray)
!!Determine a given ray is passing a cell or not
!IMPLICIT NONE
!TYPE(BasicGeom) :: GeomInfo
!REAL :: RAY(3)
!LOGICAL :: RayPassingCell
!
!REAL :: xlim(2), ylim(2)
!REAL :: Flag
!REAL :: pt(2)
!REAL :: a, b, c, m,y0
!INTEGER :: npt
!a = RAY(1); b = RAY(2); c = RAY(3)
!xlim = GeomInfo%x; ylim = GeomInfo%y
!m = -a / b; y0 = -c / b
!pts(1, 1) = xlim(1); pts(2, 1) = m*xlim(1) + y0
!Flag = pts(2,1)
!pts(1, 2) = xlim(2); pts(2, 2) = m*xlim(2) + y0
!
!END FUNCTION

FUNCTION xAxis_symPoint(x)
IMPLICIT NONE
REAL :: x(2)
REAL :: xAxis_symPoint(2)

xAxis_symPoint(1) = - X(1)
xAxis_symPoint(2) =  X(2)

END FUNCTION

FUNCTION yAxis_symPoint(x)
IMPLICIT NONE
REAL :: x(2)
REAL :: yAxis_symPoint(2)

yAxis_symPoint(1) =  X(1)
yAxis_symPoint(2) = - X(2)

END FUNCTION

FUNCTION AnnularFsrArea(r,l,n)
IMPLICIT NONE
INTEGER, intent(in) :: n
INTEGER ::  i, j
REAL, intent(in) :: r(n), l
REAL ::  AnnularFsrArea(n+1)
REAL :: theta,T(n+1),S(n+1),rr(n),ll
AnnularFsrArea =0

DO i = 1, n
  T(i)=0
  S(i)=0
  rr(i)=r(i)*r(i)
ENDDO
T(n+1)=0
S(n+1)=0
ll=l*l
IF ( r(1) <= l) THEN
  S(1)=0.25*pi*rr(1)
ELSE
  theta=0.5*pi-2*acos(l/r(1))
  S(1)=l*sqrt(rr(1)-ll)+0.5*rr(1)*theta
ENDIF
T(1)=S(1);

DO i=2,n
  IF (r(i)<=l) THEN
    S(i)=0.25*pi*rr(i)
    T(i)=S(i)
  ELSE
    theta=0.5*pi-2*acos(l/r(i))
    T(i)=l*sqrt(rr(i)-ll)+0.5*rr(i)*theta
  ENDIF
  S(i)=T(i)-T(i-1)
ENDDO
S(n+1)=ll-T(n)
DO j = 1, n + 1
  i = n + 2 - j
	AnnularFsrArea(j) = S(i) / 2
ENDDO
!AREA=S
END FUNCTION

FUNCTION AnnularRegionDivision(HalfPitch, rr, ndiv)
REAL :: AnnularRegionDivision(ndiv)
REAL :: rr(2), halfpitch
INTEGER :: ndiv

REAL :: subrr(2000), subvol(2000),rdiv(2000)
REAL :: dvol, delr, volsum
INTEGER :: i, j, ibeg, iend
INTEGER :: ridx(2000)
subrr(1:2) = (/rr(2), rr(1)/); delr = rr(1) - rr(2)
subvol = 0
subvol(1:3) = AnnularFsrArea(subrr(1:2), halfpitch, 2)
dvol = subvol(2) / ndiv
subrr(1) = rr(2)
DO i = 1, 500
  subrr(i + 1) = Subrr(i) + delr / 500
ENDDO
SubVol(1:502) = AnnularFsrArea(subrr(1:501), halfpitch, 501)
subrr(1) = rr(1)
DO i = 1, 500
  subrr(i + 1) = Subrr(i) - delr / 500
ENDDO
iend = 1
DO j = 1 , ndiv - 1
  ibeg = iend + 1
  volsum = 0
  DO i = ibeg, 501
    volsum = volsum + SubVol(i)
    IF(volsum .GT. dvol) EXIT
  ENDDO
  iend = i
  ridx(j) = iend
ENDDO
ridx(ndiv) = 501
AnnularRegionDivision(1) = subrr(1)
DO j = 2, ndiv
  AnnularRegionDivision(j) = Subrr(ridx(j-1))
ENDDO
CONTINUE
END FUNCTION
END MODULE