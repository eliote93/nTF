SUBROUTINE SetAsyGeom()
!Set Up Assembly Geometry Information
USE PARAM
USE typedef, ONLY : BasicGeom
USE GEOM, ONLY : lEdge,      lGap,                                 &
                 nCellx0,    nCellx,     AsyPitch,    CellPitch,   &
                 AsyGeom
IMPLICIT NONE
REAL :: AsyHalfPitch
INTEGER :: i ,j ,k

AsyHalfPitch = HALF * AsyPitch
AsyGeom(0)%nbd = 4
AsyGeom(0)%x(1) = -AsyHalfPitch; AsyGeom(0)%x(2) = AsyHalfPitch
AsyGeom(0)%y(1) = -AsyHalfPitch; AsyGeom(0)%y(2) = AsyHalfPitch

IF(.NOT. lEdge) THEN
  !TYPE 1
  AsyGeom(1)%nbd = 4
  AsyGeom(1)%x(1) = -AsyHalfPitch; AsyGeom(1)%x(2) = AsyHalfPitch
  AsyGeom(1)%y(1) = -AsyHalfPitch; AsyGeom(1)%y(2) = 0
  !TYPE 2
  AsyGeom(2)%nbd = 4
  AsyGeom(2)%x(1) =  0;            AsyGeom(2)%x(2) = AsyHalfPitch
  AsyGeom(2)%y(1) = -AsyHalfPitch; AsyGeom(2)%y(2) = AsyHalfPitch
  !TYPE 3
  AsyGeom(3)%nbd = 4
  AsyGeom(3)%x(1) =  0;            AsyGeom(3)%x(2) = AsyHalfPitch
  AsyGeom(3)%y(1) = -AsyHalfPitch; AsyGeom(3)%y(2) = 0
ENDIF
CALL SetRectAsyGeomBoundary(AsyGeom(0))
CALL SetRectAsyInternalGeom(AsyGeom(0), AsyPitch, CellPitch, nCellx0, lGap)
 !Set the Partial Assemblies
IF(.NOT. lEdge) THEN
  DO i = 1, 3
    CALL SetRectAsyGeomBoundary(AsyGeom(i))
    CALL SetPartailRectAsyInternalGeom(AsyGeom(i), AsyGeom(0), AsyPitch, nCellx0, lGap, i)
  ENDDO
ENDIF
DO i = 0, 3
  IF(i .NE. 0 .AND. LEdge) CYCLE
  AsyGeom(i)%lx = AsyGeom(i)%x(2) - AsyGeom(i)%x(1)
  AsyGeom(i)%ly = AsyGeom(i)%y(2) - AsyGeom(i)%y(1)
ENDDO
END SUBROUTINE

SUBROUTINE SetRectAsyGeomBoundary(AsyGeom)
USE PARAM
USE typedef, ONLY : BasicGeom
IMPLICIT NONE
TYPE(BasicGeom), INTENT(INOUT) :: AsyGeom
!REAL, INTENT(IN) :: AsyPitch
!REAL :: AsyHalfPitch
REAL :: XR, XL, YL, YR

!AsyHalfPitch = HALF * AsyPitch

XL = AsyGeom%x(1); XR = AsyGeom%x(2)
YL = AsyGeom%y(1); YR = AsyGeom%y(2)

Allocate(AsyGeom%Bdline(3, AsyGeom%nbd))
AsyGeom%Bdline(1, SOUTH) = 0.; AsyGeom%Bdline(2, SOUTH) = 1.; AsyGeom%Bdline(3, SOUTH) = -YL
AsyGeom%Bdline(1, NORTH) = 0.; AsyGeom%Bdline(2, NORTH) = 1.; AsyGeom%Bdline(3, NORTH) = -YR
AsyGeom%Bdline(1, WEST) = 1.; AsyGeom%Bdline(2, WEST) = 0.; AsyGeom%Bdline(3, WEST) = -XL
AsyGeom%Bdline(1, EAST) = 1.; AsyGeom%Bdline(2, EAST) = 0.; AsyGeom%Bdline(3, EAST) = -XR

Allocate(AsyGeom%BdlineRange(4, AsyGeom%nbd))

AsyGeom%BdlineRange(1, SOUTH) =  XL; AsyGeom%BdlineRange(2, SOUTH) =  XR    !X-Range
AsyGeom%BdlineRange(3, SOUTH) =  YL; AsyGeom%BdlineRange(4, SOUTH) =  YL    !Y-Range


AsyGeom%BdlineRange(1, NORTH) =  XL; AsyGeom%BdlineRange(2, NORTH) =  XR   !X-Range
AsyGeom%BdlineRange(3, NORTH) =  YR; AsyGeom%BdlineRange(4, NORTH) =  YR   !Y-Range

AsyGeom%BdlineRange(1, WEST) = XL; AsyGeom%BdlineRange(2, WEST) = XL   !X-Range
AsyGeom%BdlineRange(3, WEST) = YL; AsyGeom%BdlineRange(4, WEST) = YR   !Y-Range

AsyGeom%BdlineRange(1, EAST) = XR; AsyGeom%BdlineRange(2, EAST) = XR   !X-Range
AsyGeom%BdlineRange(3, EAST) = YL; AsyGeom%BdlineRange(4, EAST) = YR   !Y-Range

END SUBROUTINE

SUBROUTINE SetRectAsyInternalGeom(AsyGeom, AsyPitch, CellPitch, nCellx0, lGap)
USE PARAM
USE typedef, ONLY : BasicGeom
IMPLICIT NONE
TYPE(BasicGeom) :: AsyGeom
REAL :: AsyPitch, CellPitch
INTEGER :: nCellx0
LOGICAL :: lGap

INTEGER :: i, j, k
INTEGER :: nx, ny, nx0, ny0
REAL :: hx(0: 400), hy(0: 400)
REAL :: xlim(2), ylim(2), lx, ly, GapPitch

xlim = AsyGeom%x; ylim = AsyGeom%y
nx = nCellX0; ny = nCellX0
nx0 = nx; ny0 = ny

IF(lGap) THEN
  nx = nx + 2;  ny = ny + 2
  GapPitch = Half*(AsyPitch - CellPitch * nCellx0)
ENDIF

AsyGeom%nline = nx + ny - 2 
AsyGeom%LCircle = FALSE
AsyGeom%lRect = TRUE
ALLOCATE(AsyGeom%line(3, AsyGeom%nline))

hx(0) = xlim(1); hy(0) = ylim(1)
hx(nx) = xlim(2); hy(ny) = ylim(2)
IF(.NOT. lGap) THEN
  DO i = 1, nx - 1
    hx(i) = hx(i-1) + CellPitch
  ENDDO
  DO i = 1, ny - 1
    hy(i) = hy(i-1) + CellPitch  
  ENDDO  
ELSE
  hx(1) = hx(0) + GapPitch; hy(1) = hy(0) + GapPitch
  hx(nx-1) = hx(nx) - GapPitch; hy(ny-1) = hy(ny) - GapPitch
  DO i = 2, nx - 2
    hx(i) = hx(i-1) + CellPitch
  ENDDO
  
  DO i = 2, ny - 2
    hy(i) = hy(i-1) + CellPitch
  ENDDO

ENDIF

DO i= 1, nx-1
  AsyGeom%line(1, i) = 1.;   AsyGeom%line(2, i) = 0.
  AsyGeom%line(3, i) = -hx(i)
ENDDO

DO i= 1, ny-1
  AsyGeom%line(1, nx+i-1) = 0.;   AsyGeom%line(2, nx+i-1) = 1.
  AsyGeom%line(3, nx+i-1) = -hy(i)
ENDDO
AsyGeom%nx = nx; AsyGeom%ny = ny 
END SUBROUTINE

SUBROUTINE SetPartailRectAsyInternalGeom(AsyGeom, AsyGeom_ref, AsyPitch, nCellx0, lGap, itype)
USE PARAM
USE typedef, ONLY : BasicGeom
IMPLICIT NONE

TYPE(BasicGeom) :: AsyGeom, AsyGeom_ref
REAL :: AsyPitch
LOGICAL :: lGap
INTEGER :: iType, nCellx0
INTEGER :: i, j, k
INTEGER :: nx, ny
REAL :: hx(0: 400), hy(0: 400)
REAL :: BDLINES(3,2), lineX(3, 500), lineY(3, 500)

AsyGeom%LCircle = FALSE
BDLINES(:, 1) = AsyGeom%BDline(:, WEST)
BDLINES(:, 2) = AsyGeom%BDline(:, NORTH)
!nx0 = AsyGeom%nx; ny0 = AsyGeom%ny
nx = AsyGeom_Ref%nx; ny = AsyGeom_Ref%ny

!Y Axis 
j = 0;
DO i = 1, nx - 1
  IF((BdLines(3, 1) - AsyGeom_Ref%line(3, i)) .LT. epsm5) CYCLE
  IF(abs(BdLines(3, 1) - AsyGeom_Ref%line(3, i)) .LT. epsm5) CYCLE
  j = j + 1
  lineX(:, j) = AsyGeom_Ref%line(:, i)
ENDDO
AsyGeom%nx = j + 1
!X Axis
j = 0
DO i = 1, ny - 1
  IF((BdLines(3, 2) - AsyGeom_Ref%line(3, nx + i - 1)) .GT. epsm5) CYCLE
  IF(abs(BdLines(3, 2) - AsyGeom_Ref%line(3, nx + i - 1)) .LT. epsm5) CYCLE
  j = j + 1
  lineY(:, j) = AsyGeom_Ref%line(:, nx + i - 1)
ENDDO
AsyGeom%nY = j + 1

nx = AsyGeom%nx; ny = AsyGeom%ny
AsyGeom%lRect = TRUE; AsyGeom%LCircle = FALSE
AsyGeom%nline = nx + ny -2

ALLOCATE(AsyGeom%line(3, AsyGeom%nline)) 
DO i = 1, nx - 1
  AsyGeom%line(:, i) = LineX(:, i)
ENDDO

DO j = 1, ny - 1
  AsyGeom%line(:, nx + j - 1) = LineY(:, j)
ENDDO

END SUBROUTINE

SUBROUTINE SetUsrAsyGap()
USE PARAM
USE TYPEDEF
USE ALLOCS 
USE CNTL
USE PE_MOD,             ONLY : PE
USE inputcards ,        ONLY : oneline,      probe
USE GEOM,               ONLY : AsyGap,       nAsyGapType,                           &
                               AsyInfo,      PinInfo,                               &
                               nPinType,     nPinType0,    nCellX0,       nCellX,   &
                               nAsyType,     nAsyType0

USE BasicOperation,     ONLY : CP_CA
USE ioutil,             ONLY : terminate
IMPLICIT NONE

INTEGER :: BaseGapPin(3)
INTEGER :: BaseAsyGap2D(nCellXMax, nCellXMax)
INTEGER :: BaseAsyGap1D(nCellMax)
INTEGER :: GapRange(4, 4, 3)
INTEGER :: ix, iy, iasy, itype 
INTEGER :: i

GapRange(:, :, 1) = RESHAPE((/        2, nCellX-1,        1,        1,  &
                                      2, nCellX-1,   nCellX,   nCellX,  &
                                      0,       -1,        0,       -1, &
                                      0,       -1,        0,       -1/),&
                             (/4, 4/))
!Left-Right
GapRange(:, :, 2) = RESHAPE((/        1,        1,        2, nCellX-1,  &
                                  nCellX,   nCellX,        2, nCellX-1,  &
                                       0,       -1,        0,       -1, &
                                       0,       -1,        0,       -1/),&
                             (/4, 4/))
!Corner Gap
GapRange(:, :, 3) = RESHAPE((/      1,      1,      1,     1,  &
                               nCellX, nCellX,      1,     1,  &
                                    1,      1, nCellX, nCellX, &
                               nCellX, nCellX, nCellX, nCellX/),&
                            (/4, 4/))

BaseGapPin = (/nPinType0-2, nPinType-1, nPinType/)
CALL CP_CA(BaseAsyGap2D, 0, nCellXMax, nCellXMax)
DO itype = 1, 3
  DO i = 1, 4
    DO iy = GapRange(3, i, itype), GapRange(4, i, itype)
      DO ix = GapRange(1, i, itype), GapRange(2, i, itype)
        BaseAsyGap2D(ix, iy) = BaseGapPin(itype)
      ENDDO
    ENDDO
  ENDDO
ENDDO

i = 0
DO iy = 1, nCellX
  DO ix = 1, nCellX
    i = i + 1
    BaseAsyGap1D(i) = BaseAsyGap2D(ix, iy)
  ENDDO
ENDDO

DO iasy = 1, nAsyType0
  IF(AsyInfo(iasy)%lEmpty) CYCLE
  DO i = 1, AsyInfo(iasy)%nxy
    IF(BaseAsyGap1D(i) .EQ. 0) CYCLE
    AsyInfo(iasy)%Pin(i) = BaseAsyGap1D(i)
  ENDDO
ENDDO

IF(.NOT. nTracerCntl%lUsrGap) RETURN

DO iasy = 1, nAsyType0
  IF(AsyInfo(iasy)%lEmpty) CYCLE
  itype = AsyInfo(iasy)%GapType
  IF(itype .EQ. 0) CYCLE
  IF(itype .GT. nAsyGapType) THEN
    WRITE(mesg,'(2(a,i3),a)') 'ASYGAP (',itype,') in ASSEMBLY (',iasy,') is not defined'
    CALL terminate(mesg)        
  ENDIF
  
  DO i = 1, AsyInfo(iasy)%nxy
    IF(AsyGap(itype)%GapPin(i) .EQ. 0) CYCLE
    AsyInfo(iasy)%Pin(i) = AsyGap(itype)%GapPin(i)
  ENDDO
ENDDO

END SUBROUTINE
