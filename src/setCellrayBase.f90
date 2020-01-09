SUBROUTINE SetCellRayBase(CellInfo, RayInfo, CellRayBase)
!Make sets of Cellrays which
USE PARAM
USE typedef, ONLY : Cell_type,  CellRayBase_Type,  RayInfo_Type,             &
                    CellRayInfo_type, AziAngleInfo_Type
!USE GEOM, ONLY : CellInfo, nCellType, CellPitch
!USE rays, ONLY :  AziAngle, RayInfo
USE GeomTreatment, ONLY : GeomBoundary_RayIntersection,                       &
                          Cell_RayIntersection,                               &
                          RectAnnularCellFsrIdx,                              &
                          RectSquareCellFsrIdx
USE CadGeom_Mod,   ONLY : CadCell_RayIntersection,                            &
                          VertexCheck
USE SetRay, ONLY : BaseCellRayVolCor
USE ALLOCS
IMPLICIT NONE
TYPE(Cell_type) :: CellInfo
TYPE(RayInfo_type) :: RayInfo
TYPE(CellRayBase_Type) :: CellRayBase

TYPE(AziAngleInfo_Type), POINTER :: AziAngle(:)
TYPE(CellRayInfo_type), POINTER :: CellRay
REAL :: del, del0
!REAL, PARAMETER :: del = 0.001_8
INTEGER :: i, j, jbeg, jend, k
INTEGER :: nAziAngle, nAziAngle90, nDiv
INTEGER :: nline, nlines(0:5000), npt, nseg, ireg
INTEGER :: RegIdx(5000)
INTEGER :: InOutSurf(2)
INTEGER, ALLOCATABLE :: CellInOutSurf(:, :)
REAL :: ang, l, sinv, cosv, tanv, del_x, del_y, len
REAL :: xin, yin
REAL :: xlim(2), ylim(2), CellInOut(2,2)
REAL, ALLOCATABLE :: CellInOutPts(:, :, :)
REAL :: pts(2, 5000)       !Internal intersection Points
REAL :: intercept0, intercept, RAY(3)
EQUIVALENCE(RAY(3), intercept)
REAL :: xlen, dx
REAL, POINTER :: dyazi(:)
LOGICAL :: symbcr

symbcr=.TRUE.
symbcr=.FALSE.
!REAL ::

ALLOCATE(CellInOutSurf(2, 1200000))
ALLOCATE(CellInOutPts(2, 2, 1200000))

IF(.NOT. CellInfo%luse) THEN
 !RETURN !---- SHOULD BE MODIFIED
ENDIF
nAziAngle =RayInfo%nAziAngle
AziAngle => RayInfo%AziAngle

allocate(dyazi(naziangle))
xlim(1) = CellInfo%geom%x(1); xlim(2) = CellInfo%geom%x(2)
ylim(1) = CellInfo%geom%y(1); ylim(2) = CellInfo%geom%y(2)
xlen= abs(xlim(2)-xlim(1))
CALL DMALLOC(CellRayBase%nlines, nAziAngle)
CALL DMALLOC(CellRayBase%AziAngIdx, nAziAngle)

nLines(0) = 0; nline = 0
CellRayBase%AziAngIdx(1) = 1
! IF( RayInfo%del .LT. 0.01 )THEN
!     ndiv=nint(RayInfo%del / 0.0005_8) ! 0.05 > 0.0001 : 500, 0.0005 : 100
! ELSE
    nDiv = 50
! ENDIF
!Number of line Search and
DO i = 1, nAziAngle
  del0 = AziAngle(i)%del
  del = del0 / nDiv
  !del = 0.01/200._8;
  ang = AziAngle(i)%ang
  sinv = AziAngle(i)%sinv; cosv = AziAngle(i)%cosv
  tanv = AziAngle(i)%tanv
  del_x = abs(del / sinv)
  del_y = abs(del / cosv)
  IF(tanv .GT. 0) THEN
    intercept0 = -ABS(xlim(2)*tanv)+ylim(1)
  ELSE
    intercept0 = -ABS(xlim(1)*tanv)+ylim(1)
  ENDIF

  RAY(1) = tanv; RAY(2) = -1
  if( symbcr )THEN
    IF(i .GT. nAziAngle*3/4)THEN
        j=3.0/2.0*nAziAngle+1.0-i
        dx=dyazi(j)- del_y*HALF
    ELSE
        dx=0
    ENDIF
    intercept = (ceiling((intercept0- del_y*HALF)/del_y))*del_y - del_y*HALF
  ELSE
    dx=0
    intercept = intercept0 - del_y*HALF+dx
  ENDIF
  !number of basic Cell Ray
  DO WHILE(TRUE)
    !Ray Settup
    intercept = intercept + del_y
    !RAY(3) = intercept;
    IF(ABS(intercept) .LT. epsm5) THEN
      intercept = intercept + del_y/500.
    ENDIF
    CALL GeomBoundary_RayIntersection(CellInfo%Geom, CellInOut, InOutSurf, RAY, npt)
    IF(npt .LT. 2) THEN   !No intersection Point or
        IF(i .GT. nAziAngle/2 .AND. i .LE. nAziAngle*3/4)THEN
        dyazi(i)=xlim(2)-cellinout(1,2)
        ENDIF
      EXIT !nline = nline + 1
    ELSE
      nline = nline + 1
      CellInOutSurf(:, nline) = InOutSurf
      CellInOutPts(:, :, nline) = CellInOut
    ENDIF
  ENDDO
  nlines(i) = nline - sum(nlines(0:i-1)); CellRayBase%nlines(i) = nlines(i)
  IF(i .LT. nAziAngle) CellRayBase%AziAngIdx(i+1) = nline + 1
ENDDO
!ALLOCATTION
CellRayBase%nline = nline
ALLOCATE(CellRayBase%CellRay(nline))
!Internal Intersection Points
DO i = 1, nAziAngle
  del0 = AziAngle(i)%del
  del = del0 / nDiv
  ang = AziAngle(i)%ang
  sinv = AziAngle(i)%sinv; cosv = AziAngle(i)%cosv
  tanv = AziAngle(i)%tanv
  del_x = abs(del / sinv)
  del_y = abs(del / cosv)
  IF(tanv .GT. 0) THEN
    intercept0 = -ABS(xlim(2)*tanv)+ylim(1)
  ELSE
    intercept0 = -ABS(xlim(1)*tanv)+ylim(1)
  ENDIF
  RAY(1) = tanv; RAY(2) = -1.
  if( symbcr )THEN
    IF(i .GT. nAziAngle*3/4)THEN
        j=3.0/2.0*nAziAngle+1.0-i
        dx=dyazi(j)- del_y*HALF
    ELSE
        dx=0
    ENDIF
    intercept = (ceiling((intercept0- del_y*HALF)/del_y))*del_y - del_y*HALF
  ELSE
    dx=0
    intercept = intercept0 - del_y*HALF+dx
  ENDIF
  jbeg = CellrayBase%AziAngIdx(i); jend = CellrayBase%AziAngIdx(i) + nlines(i) - 1
  DO j = jbeg, jend
    CellRay => CellRayBase%CellRay(j)
    intercept = intercept + del_y
    !RAY(3) = intercept;
    IF(ABS(intercept) .LT. epsm5) THEN
      intercept = intercept + del_y/500.
    ENDIF
    pts = 0
    !CALL GeomBoundary_RayIntersection(CellInfo%Geom, CellInOut, InOutSurf, RAY, npt)
    IF(.NOT. CellInfo%lCad) THEN
      CALL Cell_RayIntersection(CellInfo%Geom, Ray, pts, npt, epsm4)
    ELSE
      xin = CellInOutPts(1, 1, j); yin = CellInOutPts(2, 1, j)
      CALL CadCell_RayIntersection(CellInfo%CadGeom, xin, yin,        &
                                  Ray, npt, pts, RegIdx)
!      pts(:, 1:(npt-2)) = pts(:, 2:npt-1)
!      npt = npt - 2

    ENDIF
    nSeg = npt + 1
    CellRay%nseg = nSeg                      !Number of Segment
    CALL DMALLOC(CellRay%LocalFsrIdx,nSeg)   !Fsr Index
    CALL DMALLOC(CellRay%LenSeg,nSeg)        !Length of Segments
    CALL DMALLOC(CellRay%pts, 2, nSeg+1)     !Intersection Points

    CellRay%pts(:, 1) = CellInOutPts(:, 1, j)
    CellRay%pts(:, 2:nseg) = pts(:, 1:npt)
    CellRay%pts(:, nseg+1) = CellInOutPts(:, 2, j)
    IF(npt .GT. 0) THEN
      IF(abs(pts(1,1) - CellInOutPts(1, 1, j)) .LT. epsm5) THEN
        CellRay%pts(:, 2:nseg) = CellRay%pts(:, 3:nseg+1)
        nseg = nseg - 1
      ENDIF
      IF(abs(pts(1,npt) - CellInOutPts(1, 2, j)) .LT. epsm5) THEN
        CellRay%pts(:, nseg+1) = 0
        nseg = nseg - 1
      ENDIF
       CellRay%nseg = nSeg
    ENDIF
    IF(CellInfo%lCad) THEN
      CellRay%LocalFsrIdx(1:nSeg) = RegIdx(1:nSeg)
    ENDIF
  ENDDO
ENDDO

!!!Debug
!WRITE(86, '(e20.7)') (CellInfo%Geom%circle(3,i), i=1,CellInfo%Geom%ncircle)
!
!IF(CellInfo%Geom%lrect) THEN
!  DO i = 1, CellInfo%geom%nx-1
!    WRITE(81, *)  -CellInfo%Geom%Line(3,i), -10.0
!    WRITE(81, *)  -CellInfo%Geom%Line(3,i),  10.0
!  ENDDO
!  DO i = CellInfo%geom%nx, CellInfo%geom%nx+CellInfo%geom%ny-2
!     WRITE(81, *)   -10.0, -CellInfo%Geom%Line(3,i)
!     WRITE(81, *)    10.0, -CellInfo%Geom%Line(3,i)
!  ENDDO
!ENDIF
!DO i = 5, 5
!  jbeg = CellrayBase%AziAngIdx(i); jend = CellrayBase%AziAngIdx(i) + nlines(i) - 1
!  DO j = jbeg, jend
!    CellRay => CellRayBase%CellRay(j)
!    nseg = CellRay%nSeg
!    DO k = 1, nseg
!      WRITE(85,'(2e20.5)') CellRay%pts(1, k), CellRay%pts(2, k)
!      WRITE(85,'(2e20.5)') CellRay%pts(1, k+1), CellRay%pts(2, k+1)
!    ENDDO
!  ENDDO
!ENDDO
!
!STOP

!Set CellRay Searching Tree
ALLOCATE(CellRayBase%SearchTree(nAziAngle))
DO i = 1, nAziAngle
  jbeg = CellrayBase%AziAngIdx(i); jend = CellrayBase%AziAngIdx(i) + nlines(i) - 1
  CALL MakeSearchTree(CellInOutPts(:,:,Jbeg:Jend), CellInOutSurf(:,Jbeg:Jend), nlines(i), 4, CellRayBase%SearchTree(i))
  CellRayBase%SearchTree(i)%ndiv=ndiv
ENDDO
!CALL MakeSearchTree(CellInOutPts, CellInOutSurf, nline, 4, CellRayBase%SearchTree)


!Set Fsr InDex Searching and segment
DO i = 1, nAziAngle
  jbeg = CellrayBase%AziAngIdx(i); jend = CellrayBase%AziAngIdx(i) + nlines(i) - 1
  DO j = jbeg, jend
    CellRay => CellRayBase%CellRay(j)
    nseg = CellRay%nSeg
    DO k = 1, nseg
      IF(.NOT. CellInfo%lCAD) THEN
        pts(1, 1) = (CellRay%pts (1, k) + CellRay%pts(1, k+1))*Half
        pts(2, 1) = (CellRay%pts(2, k) + CellRay%pts(2, k+1))*Half

        !Local Fsr Index Searching
        IF(.NOT. CellInfo%Geom%lrect) THEN
          ireg = RectAnnularCellFsrIdx(CellInfo, pts(:,1))
        ELSE  !Sqaure Cell
          ireg = RectSquareCellFsrIdx(CellInfo, pts(:, 1))
        ENDIF

        CellRay%LocalFsrIdx(k) = ireg
      ENDIF
      !Length of segments calulating
      len = (CellRay%pts(1, k) - CellRay%pts(1, k+1))**2 + (CellRay%pts(2, k) - CellRay%pts(2, k+1))**2
      len = sqrt(len)
      CellRay%LenSeg(k) = len
    ENDDO
  ENDDO
ENDDO

CALL BaseCellRayVolCor(RayInfo, CellInfo, CellRayBase, nDiv)
CellInfo%CellRay => CellRayBase%CellRay
NULLIFY(AziAngle)
NULLIFY(CellRay)

DEALLOCATE(CellInOutSurf)
DEALLOCATE(CellInOutPts)

END SUBROUTINE


SUBROUTINE MakeSearchTree(Pts, InOutSurf, nline, nbd, SearchTree)
!Make Information Tree in sake of easy searching for CellRay
USE PARAM
USE TYPEDEF, ONLY : CellRayBase_SearchTree
USE ALLOCS
IMPLICIT NONE
REAL :: Pts(2, 2, nline)
INTEGER :: InOutSurf(2, nline)
INTEGER :: nLine, nbd                   !Number of Cell ray and number of boundary.
TYPE(CellRayBase_SearchTree) :: SearchTree

REAL :: Mult
INTEGER :: i, j, k
INTEGER :: jbeg, jend
LOGICAL :: OverLap

CALL Dmalloc(SearchTree%nPtsInSruf, nbd)
CALL Dmalloc(SearchTree%InSurfIdxSt, nbd)
CALL Dmalloc(SearchTree%InOutSurf, 2, nline)
CALL Dmalloc(SearchTree%pts, 2, nline, 2)
SearchTree%nline = nline
!Array Starting Point Index
DO i = 1, nbd
  !Starting Points Search
  DO j = 1, nline
    IF (InoutSurf(1, j) .EQ. i) EXIT
  ENDDO
  !IF a Givne in-surf exists
  jbeg = j
  IF( jbeg .LE. nline) THEN
    SearchTree%InSurfIdxSt(i) = jbeg
    k = 0
    DO j = jbeg, nline
      IF(InoutSurf(1, j) .NE. i) EXIT
      k = k+1
    ENDDO
    SearchTree%nPtsInSruf(i) = k
  ENDIF
ENDDO
MULT = SearchTree%MULT
DO WHILE(TRUE)
  DO i = 1, nline
    SearchTree%Pts(1, i, 1) =  INT(MULT * Pts(1, 1, i))
    SearchTree%Pts(2, i, 1) =  INT(MULT * Pts(2, 1, i))
    SearchTree%Pts(1, i, 2) =  INT(MULT * Pts(1, 2, i))
    SearchTree%Pts(2, i, 2) =  INT(MULT * Pts(2, 2, i))
  ENDDO
  OverLap = TRUE
  DO i = 1, nline - 1
    IF((SearchTree%Pts(1, i, 1) - SearchTree%Pts(1, i + 1, 1)) .EQ. 0  &
      .AND. (SearchTree%Pts(2, i, 1) - SearchTree%Pts(2, i + 1, 1)) .EQ. 0) THEN
      EXIT
    ENDIF
    IF((SearchTree%Pts(1, i, 2) - SearchTree%Pts(1, i + 1, 2)) .EQ. 0 &
      .AND. (SearchTree%Pts(2, i, 2) - SearchTree%Pts(2, i + 1, 2)) .EQ. 0) THEN
      EXIT
    ENDIF
  ENDDO
  IF(i .EQ. nline) OverLap = FALSE
  IF(.NOT. OverLap) EXIT
  MULT = 10 * MULT
ENDDO

END SUBROUTINE


SUBROUTINE ChkSameGeomStr(iCellinfo, jCellinfo, lsame)
use param
use allocs
USE TYPEDEF, ONLY : Cell_TYPE, basicgeom
IMPLICIT NONE
TYPE(Cell_Type) :: iCellInfo, jCellInfo
LOGICAL :: lsame
TYPE(basicgeom) :: iGeom, jGeom

INTEGER :: i, j

iGeom=iCellInfo%geom
jGeom=jCellInfo%geom

lsame=.TRUE.
IF( igeom%nbd .NE. jgeom%nbd )THEN
    lsame=.FALSE.
    RETURN
ENDIF
IF( igeom%lcircle .NEQV. jgeom%lcircle )THEN
    lsame=.FALSE.
    RETURN
ENDIF
IF( igeom%lrect .NEQV. jgeom%lrect )THEN
    lsame=.FALSE.
    RETURN
ENDIF
IF( igeom%lccent .NEQV. jgeom%lccent )THEN
    lsame=.FALSE.
    RETURN
ENDIF
IF( igeom%ncircle .NE. jgeom%ncircle )THEN
    lsame=.FALSE.
    RETURN
ENDIF
IF( igeom%nline .NE. jgeom%nline )THEN
    lsame=.FALSE.
    RETURN
ENDIF
IF( igeom%ccenttype .NE. jgeom%ccenttype )THEN
    lsame=.FALSE.
    RETURN
ENDIF
IF( igeom%nx .NE. jgeom%nx )THEN
    lsame=.FALSE.
    RETURN
ENDIF
IF( igeom%ny .NE. jgeom%ny )THEN
    lsame=.FALSE.
    RETURN
ENDIF
IF( abs(igeom%cx - jgeom%cx) .GT. epsm10 )THEN
    lsame=.FALSE.
    RETURN
ENDIF
IF( abs(igeom%cy - jgeom%cy) .GT. epsm10 )THEN
    lsame=.FALSE.
    RETURN
ENDIF
IF( abs(igeom%x(1) - jgeom%x(1)) .GT. epsm10 )THEN
    lsame=.FALSE.
    RETURN
ENDIF
IF( abs(igeom%x(2) - jgeom%x(2)) .GT. epsm10 )THEN
    lsame=.FALSE.
    RETURN
ENDIF
IF( abs(igeom%y(1) - jgeom%y(1)) .GT. epsm10 )THEN
    lsame=.FALSE.
    RETURN
ENDIF
IF( abs(igeom%y(2) - jgeom%y(2)) .GT. epsm10 )THEN
    lsame=.FALSE.
    RETURN
ENDIF

DO i=1, igeom%nbd
    DO j=1,3
        IF( abs(igeom%bdline(j,i) - jgeom%bdline(j,i)) .GT. epsm10 )THEN
            lsame=.FALSE.
            RETURN
        ENDIF
    ENDDO
ENDDO

DO i=1, igeom%nline
    DO j=1,3
        IF( abs(igeom%line(j,i) - jgeom%line(j,i)) .GT. epsm10 )THEN
            lsame=.FALSE.
            RETURN
        ENDIF
    ENDDO
ENDDO
IF( igeom%lrect )THEN
DO i=1, igeom%nx
    IF( abs(igeom%delx(i) - jgeom%delx(i)) .GT. epsm10 )THEN
        lsame=.FALSE.
        RETURN
    ENDIF
ENDDO
DO i=1, igeom%ny
    IF( abs(igeom%dely(i) - jgeom%dely(i)) .GT. epsm10 )THEN
        lsame=.FALSE.
        RETURN
    ENDIF
ENDDO
ELSE
DO i=1, igeom%ncircle
    DO j=1,3
        IF( abs(igeom%circle(j,i) - jgeom%circle(j,i)) .GT. epsm10 )THEN
            lsame=.FALSE.
            RETURN
        ENDIF
    ENDDO
ENDDO
ENDIF

ENDSUBROUTINE
