SUBROUTINE SetPinGeom()
USE Param
USE Allocs
USE CNTL
USE GEOM,  ONLY : CellInfo,    nCellType,     nCellType0
IMPLICIT NONE
INTEGER :: iCel, iType
LOGICAL :: lCent, lRect, lCad, lCCell

DO icel = 1,  nCellType
  lCent = CellInfo(icel)%lCentX .or. CellInfo(icel)%lCentY .or. CellInfo(icel)%lCentXY
  lRect = CellInfo(icel)%lRect; lCCell = CellInfo(icel)%lCCell
  lCad = CellInfo(icel)%lCad
  IF(lCENT .AND. lCCell) CYCLE  ! Not available to set up partial cells for Corner center cell
  IF(lCad) CYCLE
  IF(.NOT. lRect) THEN
    !Annular RIng Geometries
    IF(.NOT. lCent) THEN
      CALL SetAnnularCellGeom(CellInfo(icel)%Geom, CellInfo(icel)%nDivAzi)
      CALL SetAnularCellVol(CellInfo(icel))
    ELSE  !Half Cell
      IF(CellInfo(icel)%lCentXY) THEN
        iType = 3
        CALL SetCentAnnularCellGeom(CellInfo(icel)%Geom, iType, CellInfo(icel)%nDivAzi)
        CALL SetCentAnularCellVol(CellInfo(icel), iType)
      ELSEIF(CellInfo(icel)%lCentX) THEN
        iType = 1
        CALL SetCentAnnularCellGeom(CellInfo(icel)%Geom, iType, CellInfo(icel)%nDivAzi)
        CALL SetCentAnularCellVol(CellInfo(icel), iType)
      ELSE
        iType = 2
        CALL SetCentAnnularCellGeom(CellInfo(icel)%Geom, iType, CellInfo(icel)%nDivAzi)
        CALL SetCentAnularCellVol(CellInfo(icel), iType)
      ENDIF
    ENDIF
  ELSE
    CALL SetSquareCellGeom(CellInfo(icel)%Geom)
    CALL SetRectularCellVol(CellInfo(icel))
  ENDIF
ENDDO
END SUBROUTINE

!--- JSR Edit : nTIG Restart

SUBROUTINE SetBasePinGeom()
USE Param
USE Allocs
USE CNTL
USE GEOM,  ONLY : BaseCellInfo,    nBaseCell
IMPLICIT NONE
INTEGER :: iCel, iType
LOGICAL :: lCent, lRect, lCad, lCCell

DO icel = 1,  nBaseCell
  lCent = BaseCellInfo(icel)%lCentX .or. BaseCellInfo(icel)%lCentY .or. BaseCellInfo(icel)%lCentXY
  lRect = BaseCellInfo(icel)%lRect; lCCell = BaseCellInfo(icel)%lCCell
  lCad = BaseCellInfo(icel)%lCad
  IF(lCENT .AND. lCCell) CYCLE  ! Not available to set up partial cells for Corner center cell
  IF(lCad) CYCLE
  IF(.NOT. lRect) THEN
    !Annular RIng Geometries
    IF(.NOT. lCent) THEN
      CALL SetAnnularCellGeom(BaseCellInfo(icel)%Geom, BaseCellInfo(icel)%nDivAzi)
      CALL SetAnularCellVol(BaseCellInfo(icel))
    ELSE  !Half Cell
      IF(BaseCellInfo(icel)%lCentXY) THEN
        iType = 3
        CALL SetCentAnnularCellGeom(BaseCellInfo(icel)%Geom, iType, BaseCellInfo(icel)%nDivAzi)
        CALL SetCentAnularCellVol(BaseCellInfo(icel), iType)
      ELSEIF(BaseCellInfo(icel)%lCentX) THEN
        iType = 1
        CALL SetCentAnnularCellGeom(BaseCellInfo(icel)%Geom, iType, BaseCellInfo(icel)%nDivAzi)
        CALL SetCentAnularCellVol(BaseCellInfo(icel), iType)
      ELSE
        iType = 2
        CALL SetCentAnnularCellGeom(BaseCellInfo(icel)%Geom, iType, BaseCellInfo(icel)%nDivAzi)
        CALL SetCentAnularCellVol(BaseCellInfo(icel), iType)
      ENDIF
    ENDIF
  ELSE
    CALL SetSquareCellGeom(BaseCellInfo(icel)%Geom)
    CALL SetRectularCellVol(BaseCellInfo(icel))
  ENDIF
ENDDO

CALL CopyBaseCellGeom()

END SUBROUTINE

SUBROUTINE SetAnnularCellGeom(GeomInfo, nDivAzi)
USE PARAM
USE TypeDef,  ONLY : basicgeom
USE Allocs
IMPLICIT NONE
TYPE(basicGeom) :: GeomInfo
INTEGER :: nDivAzi
INTEGER :: i, j, k
INTEGER, PARAMETER :: nbd = 4
REAL :: lx, ly, delAng, m

lx = GeomInfo%lx; ly = GeomInfo%ly
!Boundary Line Setting
!ax + by + c == 0 (a, b, c)
CALL Dmalloc(GeomInfo%bdLine, 3, 4)
GeomInfo%nbd = nbd
!SOUTH
GeomInfo%bdline(1, SOUTH) = ZERO; GeomInfo%bdline(2, SOUTH) = 1._8
GeomInfo%bdline(3, SOUTH) = -GeomInfo%y(1)!ly*Half
!WEST
GeomInfo%bdline(1, WEST) = 1._8; GeomInfo%bdline(2, WEST) = ZERO
GeomInfo%bdline(3, WEST) = -GeomInfo%x(1)
!NORTH
GeomInfo%bdline(1, NORTH) = ZERO; GeomInfo%bdline(2, NORTH) =  1._8
GeomInfo%bdline(3, NORTH) = -GeomInfo%y(2)
!EAST
GeomInfo%bdline(1, EAST) = 1._8; GeomInfo%bdline(2, EAST) = ZERO
GeomInfo%bdline(3, EAST) = -GeomInfo%x(2)

!Azimuthal Line Setting  (90, 45, 0, -45)
IF(.NOT. GeomInfo%lCCent) THEN
  IF(nDivAzi .EQ. 8) THEN
    GeomInfo%nline = 4
    CALL Dmalloc(GeomInfo%Line, 3, 4)
    GeomInfo%line(1, 1) = 1.; GeomInfo%line(2, 1) =0.; GeomInfo%line(3, 1) = 0.
    GeomInfo%line(1, 2) = 1.; GeomInfo%line(2, 2) = -1.; GeomInfo%line(3, 2) = 0.
    GeomInfo%line(1, 3) = 0.; GeomInfo%line(2, 3) =1.; GeomInfo%line(3, 3) = 0.
    GeomInfo%line(1, 4) = 1.; GeomInfo%line(2, 4) = 1.; GeomInfo%line(3, 4) = 0.
  ELSE
    GeomInfo%nline = 8
    CALL Dmalloc(GeomInfo%Line, 3, 8)
    GeomInfo%line(1, 1) =1.0; GeomInfo%line(2, 1) =0.; GeomInfo%line(3, 1) = 0.
    DelAng = PI/8.
    DO i = 2, 8
      m = PI/2 - DelAng * (i - 1); m = tan(m)
      GeomInfo%line(1, i) = m; GeomInfo%line(2, i) = -1; GeomInfo%line(3, i) = 0.
    ENDDO
  ENDIF
ELSE
  GeomInfo%nline = 1
  CALL Dmalloc(GeomInfo%Line, 3, 1)
  IF(GeomInfo%CCentType .EQ. 1 .or. GeomInfo%CCentType .EQ. 3) THEN
    GeomInfo%line(1, 1) = 1.; GeomInfo%line(2, 1) = -1.; GeomInfo%line(3, 1) = 0.
  ELSE
    GeomInfo%line(1, 1) = 1.; GeomInfo%line(2, 1) = 1.; GeomInfo%line(3, 1) = 0.
  ENDIF
ENDIF

END SUBROUTINE

SUBROUTINE SetCentAnnularCellGeom(GeomInfo, iType, nDivAzi)
USE PARAM
USE TypeDef,  ONLY : basicgeom
USE Allocs
IMPLICIT NONE
TYPE(basicGeom) :: GeomInfo
INTEGER :: iType                  !Cente Cutted Cell(? Partial Pincell)
INTEGER :: nDivAzi                !Azimuthal Angle Division

INTEGER :: i, j, k
INTEGER, PARAMETER :: nbd = 4
REAL :: lx, ly

REAL :: DelAng, m, line_data(3,8)
DATA  line_data / 1.0_8,                0._8,      0.0_8,     &
                  2.41421368121387_8,  -1.0_8,     0.0_8,     &
                  1.0_8,               -1.0_8,     0.0_8,     &
                  0.414213569169713_8, -1.0_8,     0.0_8,     &
                  0.0_8,               -1.0_8,     0.0_8,     &
                 -0.414213569169713_8, -1.0_8,     0.0_8,     &
                 -1.0_8,               -1.0_8,     0.0_8,     &
                 -2.41421368121387_8,  -1.0_8,     0.0_8 /


lx = GeomInfo%lx; ly = GeomInfo%ly
!Boundary Line Setting
!ax + by + c == 0 (a, b, c)
CALL Dmalloc(GeomInfo%bdLine, 3, 4)
GeomInfo%nbd = nbd
!SOUTH
!SOUTH
GeomInfo%bdline(1, SOUTH) = ZERO; GeomInfo%bdline(2, SOUTH) = 1._8
GeomInfo%bdline(3, SOUTH) = -GeomInfo%y(1)!ly*Half
!WEST
GeomInfo%bdline(1, WEST) = 1._8; GeomInfo%bdline(2, WEST) = ZERO
GeomInfo%bdline(3, WEST) = -GeomInfo%x(1)
!NORTH
GeomInfo%bdline(1, NORTH) = ZERO; GeomInfo%bdline(2, NORTH) =  1._8
GeomInfo%bdline(3, NORTH) = -GeomInfo%y(2)
!EAST
GeomInfo%bdline(1, EAST) = 1._8; GeomInfo%bdline(2, EAST) = ZERO
GeomInfo%bdline(3, EAST) = -GeomInfo%x(2)

IF(iType .EQ. 1) THEN
!  GeomInfo%nline = 3
!  CALL Dmalloc(GeomInfo%Line, 3, GeomInfo%nline)
!  !Azimuthal Line Setting  (-45, 90, 45)
!  GeomInfo%line(1, 1) = 1.; GeomInfo%line(2, 1) = 1.;  GeomInfo%line(3, 1) = 0.
!  GeomInfo%line(1, 2) = 1.; GeomInfo%line(2, 2) =0.;   GeomInfo%line(3, 2) = 0.
!  GeomInfo%line(1, 3) = 1.; GeomInfo%line(2, 3) = -1.; GeomInfo%line(3, 3) = 0.
  IF(nDivAzi .EQ. 4) THEN
    GeomInfo%nline = 3
    CALL Dmalloc(GeomInfo%Line, 3, GeomInfo%nline)
    GeomInfo%line(:, 1) = line_data(:, 7); GeomInfo%line(:, 2) = line_data(:, 1)
    GeomInfo%line(:, 3) = line_data(:, 3); 
  ELSE  
    GeomInfo%nline = 7
    CALL Dmalloc(GeomInfo%Line, 3, GeomInfo%nline)
    GeomInfo%line(:, 1) = line_data(:, 6); GeomInfo%line(:, 2) = line_data(:, 7)
    GeomInfo%line(:, 3) = line_data(:, 8); GeomInfo%line(:, 4) = line_data(:, 1)
    GeomInfo%line(:, 5) = line_data(:, 2); GeomInfo%line(:, 6) = line_data(:, 3)
    GeomInfo%line(:, 7) = line_data(:, 4)
  ENDIF
ELSEIF(iType .EQ. 2) THEN
  IF(nDivAzi .EQ. 4) THEN
    GeomInfo%nline = 3
    CALL Dmalloc(GeomInfo%Line, 3, GeomInfo%nline)
    GeomInfo%line(:, 1) = line_data(:, 3); GeomInfo%line(:, 2) = line_data(:, 5)
    GeomInfo%line(:, 3) = line_data(:, 7); 
  ELSE   
    GeomInfo%nline = 7
    CALL Dmalloc(GeomInfo%Line, 3, GeomInfo%nline)
    GeomInfo%line(:, 1) = line_data(:, 2); GeomInfo%line(:, 2) = line_data(:, 3)
    GeomInfo%line(:, 3) = line_data(:, 4); GeomInfo%line(:, 4) = line_data(:, 5)
    GeomInfo%line(:, 5) = line_data(:, 6); GeomInfo%line(:, 6) = line_data(:, 7)
    GeomInfo%line(:, 7) = line_data(:, 8)
  ENDIF
ELSE
  IF(nDivAzi .EQ. 2) THEN
    GeomInfo%nline = 1
    CALL Dmalloc(GeomInfo%Line, 3, GeomInfo%nline)
    !GeomInfo%line(:, 1) = line_data(:, 3)
    GeomInfo%line(:, 1) = line_data(:, 7)
  ELSE
    GeomInfo%nline = 3
    CALL Dmalloc(GeomInfo%Line, 3, GeomInfo%nline)
    !Azimuthal Line Setting  (-45)
    GeomInfo%line(:, 1) = line_data(:, 6); GeomInfo%line(:, 2) = line_data(:, 7)
    GeomInfo%line(:, 3) = line_data(:, 8)
  ENDIF
ENDIF
END SUBROUTINE


SUBROUTINE SetSquareCellGeom(GeomInfo)
USE PARAM
USE TypeDef,  ONLY : basicgeom
USE Allocs
IMPLICIT NONE

TYPE(basicGeom) :: GeomInfo
INTEGER :: i, j, k
INTEGER, PARAMETER :: nbd = 4
INTEGER :: nx, ny
REAL :: lx, ly
REAL :: delx, dely, x, y

lx = GeomInfo%lx; ly = GeomInfo%ly
nx = GeomInfo%nx; ny = GeomInfo%ny
!Boundary Line Setting
!ax + by + c == 0 (a, b, c)
CALL Dmalloc(GeomInfo%bdLine, 3, 4)
GeomInfo%nbd = nbd
!SOUTH
GeomInfo%bdline(1, SOUTH) = ZERO; GeomInfo%bdline(2, SOUTH) = 1._8
GeomInfo%bdline(3, SOUTH) = -GeomInfo%y(1)!ly*Half
!WEST
GeomInfo%bdline(1, WEST) = 1._8; GeomInfo%bdline(2, WEST) = ZERO
GeomInfo%bdline(3, WEST) = -GeomInfo%x(1)
!NORTH
GeomInfo%bdline(1, NORTH) = ZERO; GeomInfo%bdline(2, NORTH) =  1._8
GeomInfo%bdline(3, NORTH) = -GeomInfo%y(2)
!EAST
GeomInfo%bdline(1, EAST) = 1._8; GeomInfo%bdline(2, EAST) = ZERO
GeomInfo%bdline(3, EAST) = -GeomInfo%x(2)

GeomInfo%nline = nx + ny - 2
CALL Dmalloc(GeomInfo%line, 3, GeomInfo%nline)
delx = lx/nx; dely = ly/ny
x = GeomInfo%x(1)
DO i = 1, nx - 1
  !x = -lx*half + delx*i
  !x = GeomInfo%x(1) + delx*i
  x = x + GeomInfo%delx(i)
  GeomInfo%line(1, i) = 1; GeomInfo%line(2, i) = 0; GeomInfo%Line(3, i) = -x
ENDDO

y = GeomInfo%y(1)
DO j = 1, ny -1
  i = nx -1 + j
!  x = GeomInfo%y(1) + dely*j
  y = y + GeomInfo%dely(ny-j+1)
  !y = y + GeomInfo%dely(j)
  GeomInfo%line(1, i) = 0; GeomInfo%line(2, i) = 1; GeomInfo%Line(3, i) = -y
ENDDO
END SUBROUTINE

SUBROUTINE SetAnularCellVol(CellInfo)
USE PARAM
USE TypeDef,       Only : cell_type
USE GeomTreatment, ONLY : AnnularFsrArea
USE Allocs
IMPLICIT NONE
TYPE(Cell_Type) :: CellInfo
INTEGER :: ir, nDivAzi
INTEGER :: nFsr, nFxr, nCircle
REAL :: FsrVol, CellVol,Vol, Vol1, Vol2
REAL :: R(50), BaseArea(50)

nFsr = CellInfo%nFsr; nCircle = CellInfo%Geom%nCircle
nDivAzi =   CellInfo%nDivAzi 
CALL Dmalloc(CellInfo%vol, nFsr)
CellVol = CellInfo%Geom%lx * CellInfo%Geom%ly

Vol = CellVol


IF(.NOT. CellInfo%lCCell .AND. nDivAzi .EQ. 16) THEN
  R = 0
  R(1:nCircle) = CellInfo%Geom%circle(3,1:nCircle)
  DO ir = 1, nCircle + 1
    FsrVol = Vol - PI * R(ir) * R(ir); FsrVol = FsrVol/nDivAzi
    CellInfo%vol(nDivAzi*(ir-1)+1:nDivAzi*ir) = FsrVol
    Vol = PI * R(ir) * R(ir)
  ENDDO
  Vol1 = CellInfo%Geom%lx*CellInfo%Geom%lx*tan(PI/8._8)/8._8
  Vol2 = CellInfo%Geom%lx*CellInfo%Geom%lx*tan(PI/4._8)/8._8 - Vol1
  Vol1 = Vol1 - R(1) * R(1) * PI / nDivAzi
  Vol2 = Vol2 - R(1) * R(1) * PI / nDivAzi
!  vol = vol1; vol1 = vol2; vol2 = vol
  CellInfo%Vol(1:8) = (/Vol1, Vol2, Vol2, Vol1, Vol1, Vol2, Vol2, Vol1/)
  CellInfo%Vol(9:16) = (/Vol1, Vol2, Vol2, Vol1, Vol1, Vol2, Vol2, Vol1/)
ELSEIF(.NOT. CellInfo%lCCell) THEN
  DO ir = nCircle, 1, -1
    R(ir) = CellInfo%Geom%circle(3, ncircle - ir + 1)
  ENDDO
  BaseArea(1:nCircle+1) = AnnularFsrArea(R(1:nCircle), CellInfo%Geom%lx/2, nCircle)
  DO ir = 1, nCircle + 1
    CellInfo%vol(nDivAzi*(ir-1)+1:nDivAzi*ir) =  BaseArea(ir)
  ENDDO
  CONTINUE
ELSEIF(CellInfo%lCCell) THEN
  DO ir = nCircle, 1, -1
    R(ir) = CellInfo%Geom%circle(3, ncircle - ir + 1)
  ENDDO
  BaseArea(1:nCircle+1) = AnnularFsrArea(R(1:nCircle), CellInfo%Geom%lx, nCircle)
  DO ir = 1, nCircle + 1
    CellInfo%vol(nDivAzi*(ir-1)+1:nDivAzi*ir) =  BaseArea(ir)
  ENDDO
  CONTINUE  
ENDIF
Vol = SUM(CellInfo%vol(:))
END SUBROUTINE

SUBROUTINE SetCentAnularCellVol(CellInfo, iType)
USE PARAM
USE TypeDef,   Only : cell_type
USE GeomTreatment, ONLY : AnnularFsrArea
USE Allocs
IMPLICIT NONE
TYPE(Cell_Type) :: CellInfo
INTEGER :: iType
INTEGER :: ir, nDivAzi
INTEGER :: nFsr, nFxr, nCircle
REAL :: FsrVol, CellVol,Vol, vol1, vol2
REAL :: R(50), BaseArea(50)

nFsr = CellInfo%nFsr; nCircle = CellInfo%Geom%nCircle

CALL Dmalloc(CellInfo%vol, nFsr)
CellVol = CellInfo%Geom%lx * CellInfo%Geom%ly

Vol = CellVol
R = 0
R(1:nCircle) = CellInfo%Geom%circle(3,1:nCircle)
IF(iType .NE. 3) THEN
  nDivAzi = CellInfo%nDivAzi
  DO ir = nCircle, 1, -1
    R(ir) = CellInfo%Geom%circle(3, ncircle - ir + 1)
  ENDDO
  IF(CellInfo%Geom%lx > CellInfo%Geom%ly) THEN
    BaseArea(1:nCircle+1) = AnnularFsrArea(R(1:nCircle), CellInfo%Geom%ly, nCircle)
  ELSE
    BaseArea(1:nCircle+1) = AnnularFsrArea(R(1:nCircle), CellInfo%Geom%lx, nCircle)
  ENDIF
  DO ir = 1, nCircle + 1
    CellInfo%vol(nDivAzi*(ir-1)+1:nDivAzi*ir) = BaseArea(ir)
    !FsrVol = Vol - PI * R(ir) * R(ir)*0.5; FsrVol = FsrVol/nDivAzi
    !CellInfo%vol(nDivAzi*(ir-1)+1:nDivAzi*ir) = FsrVol
    !Vol = PI * R(ir) * R(ir)*0.5
  ENDDO
  IF(nDivAzi .EQ. 8) THEN
    IF(itype .eq. 1) THEN
      Vol1 = CellInfo%Geom%ly*CellInfo%Geom%ly*tan(PI/8._8)/2._8
      Vol2 = CellInfo%Geom%ly*CellInfo%Geom%ly*tan(PI/4._8)/2._8 - Vol1     
      CONTINUE        
    ELSE
      Vol1 = CellInfo%Geom%lx*CellInfo%Geom%lx*tan(PI/8._8)/2._8
      Vol2 = CellInfo%Geom%lx*CellInfo%Geom%lx*tan(PI/4._8)/2._8 - Vol1  
      CONTINUE  
    ENDIF

    Vol1 = Vol1 - R(1) * R(1) * PI / nDivAzi*half
    Vol2 = Vol2 - R(1) * R(1) * PI / nDivAzi*half
    CellInfo%Vol(1:8) = (/Vol1, Vol2, Vol2, Vol1, Vol1, Vol2, Vol2, Vol1/)
  ENDIF
ELSE
  nDivAzi = CellInfo%nDivAzi
  DO ir = nCircle, 1, -1
    R(ir) = CellInfo%Geom%circle(3, ncircle - ir + 1)
  ENDDO
  BaseArea(1:nCircle+1) = AnnularFsrArea(R(1:nCircle), CellInfo%Geom%lx, nCircle)
  DO ir = 1, nCircle + 1
    CellInfo%vol(nDivAzi*(ir-1)+1:nDivAzi*ir) = BaseArea(ir)
    !FsrVol = Vol - PI * R(ir) * R(ir)*0.25; FsrVol = FsrVol/nDivAzi
    !CellInfo%vol(nDivAzi*(ir-1)+1:nDivAzi*ir) = FsrVol
    !Vol = PI * R(ir) * R(ir) * 0.25
  ENDDO  
  IF(nDivAzi .EQ. 4) THEN
    Vol1 = CellInfo%Geom%lx*CellInfo%Geom%lx*tan(PI/8._8)/2._8
    Vol2 = CellInfo%Geom%lx*CellInfo%Geom%lx*tan(PI/4._8)/2._8 - Vol1
    Vol1 = Vol1 - R(1) * R(1) * PI / (nDivAzi * 4)
    Vol2 = Vol2 - R(1) * R(1) * PI / (nDivAzi * 4)
    CellInfo%Vol(1:4) = (/Vol1, Vol2, Vol2, Vol1/)
  ENDIF  
ENDIF
Vol = SUM(CellInfo%vol(:))
continue
END SUBROUTINE

SUBROUTINE SetRectularCellVol(CellInfo)
USE PARAM
USE TypeDef,   Only : Cell_type
USE Allocs
IMPLICIT NONE

TYPE(Cell_Type) :: CellInfo
INTEGER :: i, j, k
INTEGER :: nFsr
REAL :: FsrVol, CellVol

nFsr = CellInfo%nFsr
CellVol = CellInfo%Geom%lx * CellInfo%Geom%ly
FsrVol = CellVol/nFsr

CALL Dmalloc(CellInfo%Vol, nFsr)
!DO i = 1, nFsr
!  CellInfo%Vol(i) = FsrVol
!ENDDO
k = 0
DO j = 1, CellInfo%Geom%ny
  DO i = 1, CellInfo%Geom%nx
    k = k + 1
    FsrVol = CellInfo%Geom%delx(i) * CellInfo%Geom%dely(j)
    CellInfo%vol(k) = FsrVol
  ENDDO
ENDDO
END SUBROUTINE

