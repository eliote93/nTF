MODULE CadGeom_Mod
USE TypeDef, ONLY : Element_Type, CadGeom_Type
IMPLICIT NONE
!TYPE(CadGeom_Type) :: CadGeom

CONTAINS
SUBROUTINE SetUpElementInfo(CadGeom)
IMPLICIT NONE
TYPE(CadGeom_Type) :: CadGeom

TYPE(Element_Type), POINTER :: Element(:)
INTEGER :: nelement

nElement = CadGeom%nElement
Element => CadGeom%Element
CALL SetElementVertex(CadGeom)
CALL SetElementVolume(CadGeom)
CALL SetNeighElement(Element(1:nElement), nElement)
CALL SetTriangleLine(Element(1:nElement), nElement)

NULLIFY(Element)

END SUBROUTINE

SUBROUTINE CadCell_RayIntersection(CadGeom, xin, yin, line, npts, pts, ElementIdx)
TYPE(CadGeom_Type) :: CadGeom

TYPE(Element_Type), POINTER :: Element(:)
INTEGER :: nelement
INTEGER :: i
INTEGER, PARAMETER :: MaxPts = 1000
!REAL :: m, xin , yin, xout, yout
REAL :: line(3), xin, yin
INTEGER :: npts
REAL :: pts(2, MaxPts)

REAL :: xout, yout

!INTEGER :: Surf(MaxPts)
INTEGER :: ElementIdx(MaxPts)

INTEGER :: iElement, iSurf
INTEGER :: inElement, inSurf, OutSurf

nElement = CadGeom%nElement
Element => CadGeom%Element

!m=0.34; xin=-0.63; yin=0.01;
!line(1)=m; line(2)=-1; line(3)=-m*xin+yin;

inElement = FIndElement(xin, yin, Element, nElement, .TRUE.)
inSurf = FindSurface(xin, yin, Element(inElement), .TRUE.) 
CONTINUE
CALL VertexCheck(xin, yin, line, Element(inElement), insurf)

npts = 0
!Surf(npts) = inSurf; 
!Pts(1, npts) = xin; Pts(2, npts) = yin
DO
  CALL LineElementIntersection(xout, yout, OutSurf, line, &
                               Element(inElement), InSurf)

  ElementIdx(npts+1) = inElement
  IF(Element(inelement)%NeighSurf(OutSurf) .EQ. 0) EXIT
  npts = npts + 1
  Pts(:, npts) = (/xout, yout/)
  
  InSurf = Element(inelement)%NeighSurf(OutSurf)
  inElement = Element(inelement)%NeighElement(OutSurf)
ENDDO

END SUBROUTINE

SUBROUTINE TestCadGeom(CadGeom)
TYPE(CadGeom_Type) :: CadGeom

TYPE(Element_Type), POINTER :: Element(:)
INTEGER :: nelement
INTEGER :: i

REAL :: m, xin , yin, xout, yout
REAL :: line(3)

INTEGER, PARAMETER :: MaxPts = 1000
INTEGER :: npts
REAL :: pts(2, MaxPts)
INTEGER :: Surf(MaxPts)
INTEGER :: ElementIdx(MaxPts)

INTEGER :: iElement, iSurf
INTEGER :: inElement, inSurf, OutSurf

nElement = CadGeom%nElement
Element => CadGeom%Element

m=0.34; xin=-0.63; yin=0.01;
line(1)=m; line(2)=-1; line(3)=-m*xin+yin;

inElement = FIndElement(xin, yin, Element, nElement, .TRUE.)
inSurf = FindSurface(xin, yin, Element(inElement), .TRUE.) 
CONTINUE

npts = 1
Surf(npts) = inSurf; 
Pts(1, npts) = xin; Pts(2, npts) = yin
DO
  CALL LineElementIntersection(xout, yout, OutSurf, line, &
                               Element(inElement), InSurf)
  npts = npts + 1
  Surf(npts) = OutSurf
  Pts(:, npts) = (/xout, yout/)
  ElementIdx(npts-1) = inElement
  IF(Element(inelement)%NeighSurf(OutSurf) .EQ. 0) EXIT
  InSurf = Element(inelement)%NeighSurf(OutSurf)
  inElement = Element(inelement)%NeighElement(OutSurf)
ENDDO

END SUBROUTINE


FUNCTION FindSurface(x1, y1, Element, lBoundary)
IMPLICIT NONE
INTEGER :: FindSurface
REAL :: x1, y1
TYPE(Element_Type) :: Element
LOGICAL :: lBoundary

INTEGER :: nnode
REAL :: line(3)
REAL :: crit 
INTEGER :: i

nnode = Element%nnode

DO i = 1, nnode
  IF(lBoundary .and. Element%NeighSurf(i) .NE. 0) CYCLE
  line = Element%line(:, i)
  crit = line(1) * x1 + line(2) * y1 + line(3)
  IF(ABS(CRIT) < 1.0E-5_8) THEN
    FindSurface = i
    EXIT
  ENDIF
ENDDO
CONTINUE
END FUNCTION

Function FindElement(x1, y1, Element, nElement, lBoundary)
IMPLICIT NONE
REAL :: FindElement
REAL :: x1, y1
TYPE(Element_Type) :: Element(nElement)
INTEGER :: nElement
LOGICAL :: lBoundary
REAL :: flag1, flag2, flag3, flag
INTEGER :: i
DO i = 1, nElement
  IF(lBoundary .AND. .NOT. Element(i)%lBoundary) CYCLE
  flag1 = LineRfunction(Element(i)%line(1:3, 1), Element(i)%line(1:3, 2), x1, y1)
  flag2 = LineRfunction(Element(i)%line(1:3, 2), Element(i)%line(1:3, 3), x1, y1)
  flag3 = LineRfunction(Element(i)%line(1:3, 1), Element(i)%line(1:3, 3), x1, y1)
  flag = flag1 + flag2 + flag3
  IF(flag>1) THEN
    FindElement = i
    EXIT
  ENDIF
ENDDO
END FUNCTION

FUNCTION LineRFunction(line1, line2, x1, y1)
IMPLICIT NONE
REAL :: LineRFunction
REAL :: line1(3), line2(3)
REAL :: x1, y1
REAL :: m1, m2, crit

m1 = line1(1)*x1+line1(2)*y1+line1(3);
m2 = line2(1)*x1+line2(2)*y1+line2(3);
crit= m1 + m2 -sqrt(m1**2+m2**2);
IF(crit .GE. 0._8) THEN
!IF(crit > -1.0E-10_8) THEN
  LineRFunction = 1
ELSE
  LineRFunction = 0
ENDIF
END FUNCTION

SUBROUTINE SetTriangleLine(Element, nElement)
IMPLICIT NONE
TYPE(Element_Type) :: Element(nElement)
INTEGER :: nElement

INTEGER :: nnode
INTEGER :: i, j
REAL :: xc, yc, crit
REAL :: line(3)

DO i = 1, nElement
  nnode = Element(i)%nnode
  xc = SUM(Element(i)%x(1:nnode))/DBLE(nnode)
  yc = SUM(Element(i)%y(1:nnode))/DBLE(nnode)
  DO j = 1, nnode
    line = SetLine(Element(i)%x(j:j+1), Element(i)%y(j:j+1))
    crit = line(1) * xc + line(2) * yc + line(3);
    if(crit < 0._8) line = - line
    Element(i)%line(1:3, j) = line
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE SetNeighElement(Element, nElement)
IMPLICIT NONE
TYPE(Element_Type) :: Element(nElement)
INTEGER :: nElement

INTEGER :: nnode
INTEGER :: i, j ,k, l
INTEGER :: flag1, flag2
INTEGER :: mysurfnode(2)
LOGICAL :: lneigh

DO i = 1, nElement
  nnode = Element(i)%nnode
  DO j =1, nnode-1
    Element(i)%SurfNode(1:2, j) = Element(i)%NodeIdx(j:j+1)
  ENDDO  
  Element(i)%SurfNode(1:2, nnode) = (/Element(i)%NodeIdx(nnode), Element(i)%NodeIdx(1)/)
ENDDO

DO i = 1, nElement
  nnode = Element(i)%nnode
  DO j = 1, nnode
    lneigh = .FALSE.
    mysurfnode = Element(i)%surfnode(1:2,j)
    DO k = 1, nElement
      IF(k .EQ. i) CYCLE
      DO l = 1, Element(k)%nnode
        flag1 = mysurfnode(1) + mysurfnode(2) - Element(k)%SurfNode(1, l) - Element(k)%SurfNode(2, l)
        flag2 = mysurfnode(1) * mysurfnode(2) - Element(k)%SurfNode(1, l) * Element(k)%SurfNode(2, l)
        IF(flag1 == 0 .AND. flag2 == 0) THEN
          lneigh = .TRUE.
          Element(i)%NeighElement(j) = k
          Element(i)%NeighSurf(j) = l
          EXIT
        ENDIF
      ENDDO
      IF(lneigh) EXIT
    ENDDO
    IF(.NOT. lneigh) Element(i)%lBoundary = .TRUE.
  ENDDO

ENDDO

END SUBROUTINE


SUBROUTINE SetElementVertex(CadGeom)
IMPLICIT NONE
TYPE(CadGeom_Type) :: CadGeom

TYPE(Element_Type), POINTER :: Element(:)
REAL, POINTER :: x(:), y(:)
INTEGER :: nelement, nnode
INTEGER :: i, j

Element => CadGeom%Element
x => CadGeom%x; y => CadGeom%y
nelement = CadGeom%nelement

DO i = 1, nElement
  nnode= Element(i)%nnode
  DO j = 1, nnode
    Element(i)%x(j) = x(Element(i)%nodeidx(j))
    Element(i)%y(j) = y(Element(i)%nodeidx(j))
  ENDDO
  Element(i)%x(nnode+1) = x(Element(i)%nodeidx(1))
  Element(i)%y(nnode+1) = y(Element(i)%nodeidx(1))  
ENDDO

NULLIFY(Element)
NULLIFY(X)
NULLIFY(Y)

END SUBROUTINE

SUBROUTINE SetElementVolume(CadGeom)
IMPLICIT NONE
TYPE(CadGeom_Type) :: CadGeom
TYPE(Element_Type), POINTER :: Element
REAL :: A(2,2), area
INTEGER :: nelement
INTEGER :: i
nelement = Cadgeom%nelement

DO i = 1, nelement
  Element => CadGeom%Element(i)
  A(1,1) = Element%x(2)-Element%x(1)
  A(1,2) = Element%x(3)-Element%x(1)
  A(2,1) = Element%y(2)-Element%y(1)
  A(2,2) = Element%y(3)-Element%y(1)
  Area= 0.5_8 * abs(A(1, 1)*A(2,2)-A(1,2)*A(2,1))
  Element%vol = area
ENDDO
NULLIFY(Element)
END SUBROUTINE

FUNCTION SetLine(x,y)
REAL :: SetLine(3)
REAL :: x(2), y(2)
REAL :: m

m = (y(2) - y(1))
SetLine(1) = m
SetLine(2) = -(x(2) - x(1))
SetLine(3) = -m * x(1) + y(1) * (x(2) - x(1))

END FUNCTION

SUBROUTINE LineElementIntersection(x, y, OutSurf, line, Element, InSurf)
IMPLICIT NONE
TYPE(Element_Type) :: Element
REAL :: Line(3)
INTEGER :: InSurf

REAL :: x, y
INTEGER :: OutSurf

INTEGER :: i
INTEGER :: nnode, nsurf


REAL :: l1, l2
REAL, PARAMETER :: del = 1.0E-5_8
REAL :: x1, y1, vertex1(2), vertex2(2)
LOGICAL :: lexist, lLineIn

nnode = Element%nnode; nsurf = Element%nsurf
OutSurf = 0
DO i = 1, nsurf
  IF(i .EQ. InSurf) CYCLE 
  CALL LineLineIntersection(x1, y1, lexist, Line, Element%line(1:3, i))
  IF(.NOT. lexist) CYCLE
  Vertex1 = (/Element%x(i), Element%y(i)/)
  Vertex2 = (/Element%x(i+1), Element%y(i+1)/)
  !LineInOut
  lLineIn = LineInOut(x1, y1, Vertex1, Vertex2)
  IF(lLineIn) THEN
     CALL VertexCheck(x1, y1, line, Element, i)
!    l1 = (Vertex1(1) -x1)**2+(Vertex1(2) -y1)**2; l1 = sqrt(l1)
!    l2 = (Vertex2(1) -x1)**2+(Vertex2(2) -y1)**2; l2 = sqrt(l2)
!    IF(l1 < 1.0E-6 .OR. l2 < 1.0E-6) THEN
!      IF(l1 < 1.0E-6) THEN
!        x1 = (1._8 - del) * Vertex1(1) + del * Vertex2(1)
!        y1 = (1._8 - del) * Vertex1(2) + del * Vertex2(2)
!      ELSE
!        x1 = del * Vertex1(1) + (1._8 - del) * Vertex2(1)
!        y1 = del * Vertex1(2) + (1._8 - del) * Vertex2(2)      
!      ENDIF
!      line(3) = -x1 * line(1) - y1 * line(2)
!    ENDIF  
    x = x1; y = y1
    OutSurf = i
    EXIT
  ENDIF
ENDDO

END SUBROUTINE

SUBROUTINE LineLineIntersection(x1, y1, lexist, Line1, Line2)
IMPLICIT NONE
REAL :: x1, y1
LOGICAL :: lexist
REAL :: line1(3), line2(3)
REAL :: A(2, 2), AInv(2, 2), C(2)
REAL :: detA

A(1, 1:2) = (/Line1(1), Line1(2)/)
A(2, 1:2) = (/Line2(1), Line2(2)/)
C(1:2) = (/-Line1(3), -Line2(3)/)

detA = A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)
IF(ABS(DetA) .LT. 10.E-6_8) THEN
  lexist = .FALSE.
  RETURN
ENDIF

lexist = .TRUE.
AInv = MatInv2x2(A)
x1 = AInv(1, 1) * C(1) + AInv(1, 2) * C(2)
y1 = AInv(2, 1) * C(1) + AInv(2, 2) * C(2)

END SUBROUTINE

FUNCTION LineInOut(x1, y1, vertex1, vertex2)
IMPLICIT NONE
LOGICAL :: LineInOut
REAL :: x1, y1
REAL :: Vertex1(2), Vertex2(2)
REAL :: flag1, flag2
LineInOut = .FALSE.

flag1=(x1-vertex1(1))*(x1-vertex2(1))
flag2=(y1-vertex1(2))*(y1-vertex2(2))
IF(flag1 .LE. 1.0E-10_8 .AND. flag2 .LE. 1.0E-10_8) LineInOut = .TRUE.
END FUNCTION

FUNCTION MatInv2x2(A)
IMPLICIT NONE
REAL :: A(2,2), MatInv2x2(2,2)
REAL :: AInv(2,2)
REAL :: detA

detA = A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)

AInv(1, 1) = A(2, 2); AInv(2, 2) = A(1, 1)
AInv(1, 2) = -A(1, 2); AInv(2, 1) = -A(2, 1)
MatInv2x2 = AInv / detA
END FUNCTION

SUBROUTINE VertexCheck(x1, y1, line, Element, isurf)
IMPLICIT NONE
TYPE(Element_Type) :: Element
REAL :: x1, y1, line(3)
INTEGER :: isurf
REAL :: vertex1(2), vertex2(2)
REAL :: l1, l2
REAL, PARAMETER :: del = 1.0E-5_8
Vertex1 = (/Element%x(isurf), Element%y(isurf)/)
Vertex2 = (/Element%x(isurf+1), Element%y(isurf+1)/)
l1 = (Vertex1(1) -x1)**2+(Vertex1(2) -y1)**2; l1 = sqrt(l1)
l2 = (Vertex2(1) -x1)**2+(Vertex2(2) -y1)**2; l2 = sqrt(l2)
IF(l1 < 1.0E-6 .OR. l2 < 1.0E-6) THEN
  IF(l1 < 1.0E-6) THEN
    x1 = (1._8 - del) * Vertex1(1) + del * Vertex2(1)
    y1 = (1._8 - del) * Vertex1(2) + del * Vertex2(2)
  ELSE
    x1 = del * Vertex1(1) + (1._8 - del) * Vertex2(1)
    y1 = del * Vertex1(2) + (1._8 - del) * Vertex2(2)      
  ENDIF
  line(3) = -x1 * line(1) - y1 * line(2)
ENDIF  
END SUBROUTINE
END MODULE