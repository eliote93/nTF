SUBROUTINE ReadNInit_CadCell(Dataline0)
USE PARAM
USE TYPEDEF
USE ALLOCS
USE cntl,         ONLY : nTracerCntl
USE geom,         ONLY : CellInfo,  CellPitch,    nCellType
USE ioutil,       ONLY : toupper,   IFnumeric,    nfields,   fndchara,    nfieldto
USE BenchXs,      ONLY : MacXsBen
USE Material_Mod, ONLY : Mixture
USE files,        ONLY : io_quick
USE ioutil,       ONLY : openfile
USE CadGeom_mod,  ONLY : SetUpElementInfo
IMPLICIT NONE

character(256),intent(in) :: dataline0

TYPE(CadGeom_Type), POINTER :: CadGeom
TYPE(Element_Type), POINTER :: Element(:)
character(256) :: dataline !,chline
character(256) :: CadFilename

REAL :: lx, ly
REAL :: vol(100), areasum
INTEGER :: icel, ipos(100), iReg(100)
INTEGER :: nDataField, ndata, nspt
INTEGER :: i, j, K
INTEGER :: ifxr

LOGICAL :: lvol

!DATA vol(1:4) / 0.8912455511249344_8,    &
!           0.1604957062832235_8,    &
!           0.025485981906688063_8,  &
!           0.5255287606851541/
!INTEGER :: nFsrInFxr

dataline = dataline0
READ(dataline, *) icel
CellInfo(icel)%lCad = TRUE
CellInfo(icel)%lempty = FALSE; CellInfo(icel)%lrect = FALSE
CellInfo(icel)%geom%lCCent  = FALSE
CellInfo(icel)%lFuel = FALSE; CellInfo(icel)%lRes = FALSE
CellInfo(icel)%lgap = FALSE
CellInfo(icel)%lMox = FALSE
CellInfo(icel)%lCentX = FALSE; CellInfo(icel)%lCentY = FALSE; CellInfo(icel)%lCentXY = FALSE; 
lvol = .FALSE.

nDataField = len_trim(dataline)
CALL fndchara(dataline,ipos,nspt,SLASH)
ndata = nfieldto(dataline, SLASH)-1
READ(dataline, *) icel, ireg(1:ndata)

READ(dataline(ipos(nspt)+1:nDataField), *) CadFilename
ALLOCATE(CellInfo(icel)%CadGeom)

CALL OpenFile(io_quick, TRUE, FALSE, FALSE, CadFileName)
CALL ReadMesh(io_quick, CellInfo(icel)%CadGeom)
CALL SetUpElementInfo(CellInfo(icel)%CadGeom)
CLOSE(io_quick)

CadGeom => CellInfo(icel)%CadGeom
Element => CadGeom%Element

CellInfo(icel)%Geom%lx = CellPitch; CellInfo(icel)%geom%ly = CellPitch
CellInfo(icel)%Geom%x(1) = -CellPitch*Half; CellInfo(icel)%Geom%x(2) = CellPitch*Half
CellInfo(icel)%Geom%y(1) = -CellPitch*Half; CellInfo(icel)%Geom%y(2) = CellPitch*Half

CellInfo(icel)%nFSR = CadGeom%nElement
CellInfo(icel)%nFXR = ndata
CellInfo(icel)%nBd = 4

DO i = 1, ndata
  IF(nTracerCntl%lXsLib) THEN
    CellInfo(icel)%lfuel = CellInfo(icel)%lfuel .or. Mixture(ireg(i))%lfuel
    CellInfo(icel)%lGd = CellInfo(icel)%lGd .or. Mixture(ireg(i))%lGd
    CellInfo(icel)%lMOX = CellInfo(icel)%lMOX .or. Mixture(ireg(i))%lMOX
  ELSE
    CellInfo(icel)%lfuel = CellInfo(icel)%lfuel .or. MacXsBen(ireg(i))%lfuel
  ENDIF  
ENDDO

ndata = nfieldto(dataline(ipos(1)+1:nDataField), SLASH)
IF(ndata .GT. 0) THEN
  READ(dataline(ipos(1)+1:nDataField), *) vol(1:ndata)
  lvol = .TRUE.
ENDIF

CALL dmalloc(CellInfo(icel)%iReg, CellInfo(icel)%nFSR)
CALL dmalloc(CellInfo(icel)%FxrIdxSt,CellInfo(icel)%nFxr)
CALL dmalloc(CellInfo(icel)%nFsrInFxr,CellInfo(icel)%nFxr)
CALL dmalloc(CellInfo(icel)%MapFxr2FsrIdx, CellInfo(icel)%nFSR, CellInfo(icel)%nFxr)

DO i = 1, CellInfo(icel)%nFSR
  j = Element(i)%GeomType
  CellInfo(icel)%nFsrInFxr(j) = CellInfo(icel)%nFsrInFxr(j) + 1
  ifxr = CellInfo(icel)%nFsrInFxr(j)
  CellInfo(icel)%MapFxr2FsrIdx(ifxr, j) = i
  CellInfo(icel)%iReg(i) = ireg(j)
ENDDO
lx = CellInfo(icel)%Geom%lx; ly = CellInfo(icel)%Geom%ly

CALL Dmalloc(CellInfo(icel)%vol, CellInfo(icel)%nFSR)
DO i = 1, CellInfo(icel)%nFSR
  CellInfo(icel)%vol(i) = Element(i)%vol
ENDDO


!Boundary Line Setting
CALL Dmalloc(CellInfo(icel)%Geom%bdLine, 3, 4)
CellInfo(icel)%Geom%nbd = 4
!SOUTH
CellInfo(icel)%Geom%bdline(1, SOUTH) = ZERO;
CellInfo(icel)%Geom%bdline(2, SOUTH) = 1._8
CellInfo(icel)%Geom%bdline(3, SOUTH) = -CellInfo(icel)%Geom%y(1)
!WEST
CellInfo(icel)%Geom%bdline(1, WEST) = 1._8;
CellInfo(icel)%Geom%bdline(2, WEST) = ZERO
CellInfo(icel)%Geom%bdline(3, WEST) = -CellInfo(icel)%Geom%x(1)
!NORTH
CellInfo(icel)%Geom%bdline(1, NORTH) = ZERO; 
CellInfo(icel)%Geom%bdline(2, NORTH) =  1._8
CellInfo(icel)%Geom%bdline(3, NORTH) = -CellInfo(icel)%Geom%y(2)
!EAST
CellInfo(icel)%Geom%bdline(1, EAST) = 1._8; 
CellInfo(icel)%Geom%bdline(2, EAST) = ZERO
CellInfo(icel)%Geom%bdline(3, EAST) = -CellInfo(icel)%Geom%x(2)
IF(lvol) THEN
  DO i = 1, CellInfo(icel)%nFXR
    areasum = 0
    DO j=1, CellInfo(icel)%nFsrInFxr(i)
      k = CellInfo(icel)%MapFxr2FsrIdx(j, i)
      areasum = areasum + CellInfo(icel)%vol(k)
    ENDDO
    areasum = vol(i)/areasum
    DO j=1, CellInfo(icel)%nFsrInFxr(i)
      k = CellInfo(icel)%MapFxr2FsrIdx(j, i)
      Element(k)%vol=Element(k)%vol*areasum
      CellInfo(icel)%vol(k)=CellInfo(icel)%vol(k)*areasum
    ENDDO  
  ENDDO
ENDIF
!  j = Element(i)%GeomType
!  CellInfo(icel)%nFsrInFxr(j) = CellInfo(icel)%nFsrInFxr(j) + 1
!  ifxr = CellInfo(icel)%nFsrInFxr(j)
!  CellInfo(icel)%MapFxr2FsrIdx(ifxr, j) = i
!  CellInfo(icel)%iReg(i) = ireg(j)

NULLIFY(CadGeom)
NULLIFY(Element)
END SUBROUTINE

SUBROUTINE ReadMesh(indev, CadGeom)
USE TYPEDEF, ONLY : CadGeom_Type
IMPLICIT NONE
TYPE(CadGeom_Type) :: CadGeom
INTEGER, INTENT(in) :: indev
CHARACTER(1024) :: OneLine
CHARACTER(1) :: Probe
EQUIVALENCE(Probe, OneLine(1:1))
DO 
  READ(indev, '(a1024)', END=5) OneLine
  IF(PROBE .EQ. '$') THEN
    SELECT CASE(OneLine)
      CASE('$MeshFormat')
        CALL ReadMEshFormat(indev, CadGeom%MeshFormat)
      CASE('$Nodes')
        CALL ReadNodes(indev, CadGeom)
      CASE('$Elements')
        CALL ReadElements(indev, CadGeom)
    END SELECT
  ENDIF
ENDDO
5 CONTINUE
END SUBROUTINE

SUBROUTINE ReadMeshFormat(indev, MeshFormat)
IMPLICIT NONE
INTEGER, INTENT(in) :: indev
CHARACTER(50) :: MEshFormat
READ(indev, '(a50)') MeshFormat
READ(indev, *)
END SUBROUTINE

SUBROUTINE ReadNodes(indev, CadGeom)
USE TYPEDEF, ONLY : CadGeom_Type
USE allocs
IMPLICIT NONE
INTEGER :: indev
TYPE(CadGeom_Type) :: CadGeom
!CHARACTER(1024) :: OneLine
INTEGER :: i,inode
INTEGER :: nnode

READ(indev, *) CadGeom%nnode
nnode = CadGeom%nnode
CALL Dmalloc(CadGeom%x, nnode)
CALL Dmalloc(CadGeom%y, nnode)
DO i = 1, nnode
  READ(indev, *) inode, CadGeom%x(i), CadGeom%y(i)
ENDDO
READ(indev, *)
CONTINUE
END SUBROUTINE

SUBROUTINE ReadElements(indev, CadGeom)
USE TYPEDEF, ONLY : CadGeom_Type, Element_Type
USE allocs
IMPLICIT NONE
INTEGER :: indev
TYPE(CadGeom_Type) :: CadGeom
TYPE(Element_Type), POINTER :: Element(:)
INTEGER :: i,j,k,itype, idum
INTEGER :: nElement, ntag, nnode
CHARACTER(1024) :: OneLine

READ(indev, *) nElement
ALLOCATE(CadGeom%Element(nElement))
Element => CadGeom%Element
CadGeom%nElement=nElement
DO i = 1, nElement
  READ(indev, '(a1024)')  OneLine
  READ(OneLine, *) j, itype, ntag
  IF(itype .NE. 2) THEN
    PRINT '(A8,I7, A25)', 'Element', i, 'is not a triangle Element'
    STOP
  ENDIF
  nnode = 3
  READ(OneLine, *) j, itype, ntag, Element(i)%PhyType, Element(i)%GeomType, idum, Element(i)%NodeIdx(1:nnode)
ENDDO

READ(indev, *)
NULLIFY(Element)
END SUBROUTINE