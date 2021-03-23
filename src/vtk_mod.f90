#include <defines.h>
#define CENT_OPT
MODULE VTK_mod
USE PARAM
USE TYPEDEF,   ONLY : CoreInfo_Type, CmInfo_Type, FmInfo_Type, Pin_Type, Cell_Type,   &
                     basicgeom, PowerDist_Type, PE_TYPE, GroupInfo_Type,  &
                     ThInfo_Type
USE files,     ONLY : io9, io10, io11, io12, caseid, localfn
USE CNTL,      ONLY : nTracerCntl_TYPE
USE Depl_Mod,  ONLY : DeplCntl
USE Tran_MOD,  ONLY : TranCntl
IMPLICIT NONE
TYPE BaseCellVtk_TYPE

  REAL, POINTER :: pt(:, :)
  INTEGER, POINTER :: elements(:, :)
  INTEGER, POINTER :: types(:)
  INTEGER :: npt, nelmt
END TYPE

TYPE CellVtk_TYPE
  REAL, POINTER :: pt(:, :)
  INTEGER, POINTER :: elements(:, :)
  INTEGER, POINTER :: types(:)
  INTEGER :: npt = 0
  INTEGER :: nelmt = 0
  LOGICAL, POINTER :: lBdPt(:)
  REAL, POINTER :: BdDel(:, :), BdPtList(:, :, :)
  INTEGER :: nbdpt = 0
  INTEGER  :: ninpt = 0
  INTEGER :: nbddel(4) = 0
END TYPE

TYPE PlaneConfig_TYPE
  REAL, POINTER :: CellBdPt(:, :, :, :) !(x & y, idx , surf, ixy)
  INTEGER, POINTER :: PtIdx(:, :, :)
  INTEGER, POINTER :: NODE(:, :)
  REAL, POINTER :: CX(:), CY(:)
  INTEGER :: npt
END TYPE

TYPE VisCntl_TYPE
  INTEGER :: VisMod = 1       ! 1 : Plane-wise 2D visualizatoin 2: 3-D visualization 3:
  INTEGER :: FluxOutMod = 1   ! 1 : Entire Flux, 2 : 2 Group Condensed
  INTEGER :: ThermalIdxSt     !
  INTEGER :: ThermalIdxSt47G = 36
  INTEGER :: ThermalIdxSt190G = 159
END TYPE


TYPE(PlaneConfig_TYPE), POINTER, PRIVATE :: PlaneConfig(:)
TYPE(CellVtk_TYPE), POINTER, PRIVATE :: CellVtk(:)
TYPE(CellVtk_TYPE), POINTER, PRIVATE :: CellVtk3D(:)
INTEGER :: nGlobalPt(100), nGlobalElmt(100)
CHARACTER(256) :: fn_base
CONTAINS

FUNCTION GetId(CycleNo, ProcNo)
CHARACTER(8) :: GetId
INTEGER :: CycleNo, ProcNo
CHARACTER(4) :: CycleId, ProcID
INTEGER :: CycleNo0, ProcNo0

INTEGER :: i, j, k
CycleId = '0000'; ProcId = '0000'
k = 1000
CycleNo0 = CycleNo
DO i = 1, 4
  j = CycleNo0 / k
  WRITE(CycleId(i:i), '(I1)') j
  CycleNo0 = CycleNo0 - k * j
  k = k/10
ENDDO

k = 1000
ProcNo0 = ProcNo
DO i = 1, 4
  j = ProcNo0 / k
  WRITE(ProcID(i:i), '(I1)') j
  ProcNo0 = ProcNo0 - k * j
  k = k/10
ENDDO
GetId(1:4) = CycleId(1:4)
GetId(5:8) = ProcId(1:4)
END FUNCTION

SUBROUTINE ProcessVTK(Core, FmInfo, CmInfo, PowerDist, ThInfo, GroupInfo, nTracerCntl, PE)
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(PowerDist_Type) ::  PowerDist
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

IF(nTRACERCntl%OutpCntl%VisMod .EQ. 1) CALL ProcessVTK3D(Core, FMInfo, CmInfo, PowerDist, ThInfo, GroupInfo, nTracerCntl, PE)
IF(nTRACERCntl%OutpCntl%VisMod .EQ. 2) CALL ProcessVTK2DPln(Core, FMInfo, CMInfo, PowerDist, ThInfo, GroupInfo, nTracerCntl, PE)
IF(nTRACERCntl%OutpCntl%VisMod .EQ. 3) CALL ProcessVTK3D_HOM(Core, FmInfo, CMInfo, PowerDist, ThInfo, GroupInfo, nTracerCntl, PE)
END SUBROUTINE

SUBROUTINE ProcessVTK3D_HOM(Core, FmInfo, CmInfo, PowerDist, ThInfo, GroupInfo, nTracerCntl, PE)
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(PowerDist_Type) ::  PowerDist
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
INTEGER :: i, j, k, iz, myrank, nproc, myzb, myze
INTEGER :: icycle
CHARACTER(256) :: fn
CHARACTER(8) :: OutputId
REAL :: fluxnorm0

IF(PE%MASTER) PRINT *, 'WRITING VTK file...'
myrank = PE%myCMFDrank
nproc = PE%nCmfdProc
myzb = PE%myzb
myze =PE%myze
DO i = 1, 80
  if(caseid(i:i) .EQ. '') EXIT
ENDDO
i = i-1

!icycle = 0;
!IF(nTracerCntl%lProblem .EQ. lTransient) THEN
!  icycle = TranCntl%NowWriteStep
!ELSEIF(nTracerCntl%lProblem .EQ. ldepletion) THEN
!  icycle = DeplCntl%NowStep
!ENDIF
icycle = nTracerCntl%CalNoId
OutputId = GetId(icycle, myrank)
WRITE(fn,'(A,A,A8,A)') caseid(1:i), '_HOM_', OutputId, '.vtk'
!IF(PE%nCMFDproc .LT. 10) THEN
!  WRITE(fn,'(A,A,A,I1,A)') caseid(1:i), '_HOM_proc','00', myrank, '.vtk'
!ELSEIF(PE%nCMFDproc .LT. 100) THEN
!  WRITE(fn,'(A,A,A,I2,A)') caseid(1:i), '_HOM_proc','0', myrank, '.vtk'
!ELSE
!  WRITE(fn,'(A,A,I3,A)') caseid(1:i), '_HOM_proc', myrank, '.vtk'
!ENDIF

OPEN(unit=io9, file = fn, status = 'replace')

CALL WriteGeomVTK3D_HOM2(Core, myzb, myze)
WRITE(io9, '(A, x, 2I10)') 'CELL_DATA', Core%nxy * (myze- myzb+1)
fluxnorm0 = FluxNorm(Core, FmInfo, CmInfo, GroupInfo%ng, nTracerCntl, PE)
CALL WriteFlux_HOM(Core, CmInfo,fluxnorm0, GroupInfo%ng, myzb, myze)
CALL WritePower_HOM(Core, CmInfo,PowerDist%Pin3DNormalizer, GroupInfo%ng, myzb, myze)
IF(nTracerCntl%lFeedback) THEN
  CALL WriteTCool_HOM(Core, ThInfo, GroupInfo%ng, myzb, myze)
  CALL WriteTFuelAvg_HOM(Core, ThInfo, myzb, myze)
ENDIF
CLOSE(io9)

OutputId = GetId(0, icycle)
IF(PE%CMFDmaster) THEN
  WRITE(fn,'(A,A,A4,A)') caseid(1:i), '_HOM_', OutputId(5:8),'.visit'
  OPEN(unit = io12, file=fn, status = 'replace')
  WRITE(io12, '(A, 1x, I5)') '!NBLOCKS', nproc
  DO k = 0, nproc-1
    fn = ''
    icycle = 0;
    !IF(nTracerCntl%lProblem .EQ. lTransient) THEN
    !  icycle = TranCntl%NowWriteStep
    !ELSEIF(nTracerCntl%lProblem .EQ. ldepletion) THEN
    !  icycle = DeplCntl%NowStep
    !ENDIF
    icycle = nTracerCntl%CalNoId
    OutputId = GetId(icycle, k)
    IF(PE%nCMFDproc .LT. 10) THEN
      WRITE(fn,'(A,A,A,I1,A)') caseid(1:i), '_HOM_proc','00', k, '.vtk'
    ELSEIF(PE%nCMFDproc .LT. 100) THEN
      WRITE(fn,'(A,A,A,I2,A)') caseid(1:i), '_HOM_proc','0', k, '.vtk'
    ELSE
      WRITE(fn,'(A,A,I3,A)') caseid(1:i), '_HOM_proc', k, '.vtk'
    ENDIF
    !WRITE(fn,'(A,A,A8,A)') caseid(1:i), '_HOM_', OutputId, '.vtk'
    DO j = 1, 100
      IF(fn(j:j)=='') exit
    ENDDO
    j = j -1
    WRITE(io12, '(A)') fn(1:j)
  ENDDO
  CLOSE(io12)
ENDIF

END SUBROUTINE

SUBROUTINE ProcessVTK3DModel(Core, FMInfo, GroupInfo, nTracerCntl, PE)
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

TYPE(Cell_Type), POINTER :: Cell(:)
INTEGER :: nCell, ng
INTEGER :: i, j,  iz,n
REAL :: FluxNorm0

CHARACTER(256) :: fn, fn0, caseid0

IF(PE%MASTER) PRINT *, 'WRITING VTK file...'

ng = GroupInfo%ng
Cell => Core%CellInfo
nCell = Core%nCellType

ALLOCATE(CellVTK(nCell))
ALLOCATE(CellVTK3D(nCell))
DO i = 1, nCell
  IF(.NOT. Cell(i)%luse) CYCLE
  IF(.NOT. Cell(i)%lRect) THEN
    CALL MakeAnnualCellVTK(Cell(i), i)
    IF(Cell(i)%lCentX .OR. Cell(i)%lCentY .OR. Cell(i)%lCentXY) THEN
      CALL ProcessParialCell(Cell(i), i)
    ENDIF
    IF(Cell(i)%lCCell) THEN
      CALL ProcessParialCell(Cell(i), i)
    ENDIF
  ENDIF
  IF(Cell(i)%lRect) THEN
    CALL MakeRectCellVTK(Cell(i), i)
  ENDIF
  CALL ConvertCellVTK3D(i)
ENDDO

DO i = 1, 80
  if(caseid(i:i) .EQ. '') EXIT
ENDDO
i = i-1
caseid0(1:i) = caseid(1:i)

DO iz = PE%myzb, PE%myze
  IF(iz .LT. 10) THEN
    WRITE(fn,'(A,A,A,I1,A)') caseid0(1:i), '_3Dpln','00', iz, '.vtk'
    WRITE(fn_base,'(A,A,A,I1,A)') caseid0(1:i), '_3Dpln','00', iz
  ELSEIF(iz .LT.100) THEN
    WRITE(fn,'(A,A,A,I2,A)') caseid0(1:i), '_3Dpln','0', iz, '.vtk'
    WRITE(fn_base,'(A,A,A,I2,A)') caseid0(1:i), '_3Dpln','0', iz
  ELSE
    WRITE(fn,'(A,A,I3,A)') caseid0(1:i), '_3Dpln', iz, '.vtk'
    WRITE(fn_base,'(A,A,I3,A)') caseid0(1:i), '_3Dpln', iz
  ENDIF
  OPEN(unit=io9, file = fn, status = 'replace')
  CALL WriteGeomVTK3D(Core, iz)
  WRITE(io9, '(A, x, 2I10)') 'CELL_DATA', nGlobalelmt(iz)

  CALL Writematerial(Core, FmInfo, iz, nTracerCntl)
  CLOSE(io9)
ENDDO

WRITE(fn,'(A,A,A)') caseid0(1:i), '_3D','.visit'
OPEN(unit = io12, file=fn, status = 'replace')
WRITE(io12, '(A, 1x, I5)') '!NBLOCKS', COre%nz
DO iz = 1, CORE%nz
  fn = ''
  IF(iz .LT. 10) THEN
    WRITE(fn,'(A,A,A,I1,A)') caseid0(1:i), '_3Dpln','00', iz, '.vtk'
  ELSEIF(iz .LT.100) THEN
    WRITE(fn,'(A,A,A,I2,A)') caseid0(1:i), '_3Dpln','0', iz, '.vtk'
  ELSE
    WRITE(fn,'(A,A,I3,A)') caseid0(1:i), '_3Dpln', iz, '.vtk'
  ENDIF
  DO j = 1, 80
    if(fn(j:j) .EQ. '') EXIT
  ENDDO
  j = j-1
  WRITE(io12, '(A)') fn(1:j)
ENDDO
close(io12)
END SUBROUTINE




SUBROUTINE ProcessVTK3D(Core, FMInfo, CmInfo, PowerDist, ThInfo, GroupInfo, nTracerCntl, PE)
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(PowerDist_Type) ::  PowerDist
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

TYPE(Cell_Type), POINTER :: Cell(:)
INTEGER :: nCell, ng
INTEGER :: i, j,  iz,n
REAL :: FluxNorm0

CHARACTER(256) :: fn, fn0, caseid0


IF(PE%MASTER) PRINT *, 'WRITING VTK file...'

ng = GroupInfo%ng
Cell => Core%CellInfo
nCell = Core%nCellType

ALLOCATE(CellVTK(nCell))
ALLOCATE(CellVTK3D(nCell))
DO i = 1, nCell
  IF(.NOT. Cell(i)%luse) CYCLE
  IF(.NOT. Cell(i)%lRect) THEN
    CALL MakeAnnualCellVTK(Cell(i), i)
    IF(Cell(i)%lCentX .OR. Cell(i)%lCentY .OR. Cell(i)%lCentXY) THEN
      CALL ProcessParialCell(Cell(i), i)
    ENDIF
    IF(Cell(i)%lCCell) THEN
      CALL ProcessParialCell(Cell(i), i)
    ENDIF
  ENDIF
  IF(Cell(i)%lRect) THEN
    CALL MakeRectCellVTK(Cell(i), i)
  ENDIF
  CALL ConvertCellVTK3D(i)

  !CALL WriteCellVtk3D(i)
ENDDO

DO i = 1, 80
  if(caseid(i:i) .EQ. '') EXIT
ENDDO
i = i-1
caseid0(1:i) = caseid(1:i)
IF(nTracerCntl%CalNoId .LT. 10) THEN
  WRITE(caseid0(i+1:256),'(A, A, I1)') '_cycle', '00', nTracerCntl%CalNoId
ELSEIF(nTracerCntl%CalNoId .LT. 100) THEN
  WRITE(caseid0(i+1:256),'(A, A, I2)') '_cycle', '0', nTracerCntl%CalNoId
ELSE
  WRITE(caseid0(i+1:256),'(A, I3)') '_cycle', nTracerCntl%CalNoId
ENDIF
i=i+9

fluxnorm0 = FluxNorm(Core, FmInfo, CmInfo, ng, nTracerCntl, PE)
DO iz = PE%myzb, PE%myze
  IF(iz .LT. 10) THEN
    WRITE(fn,'(A,A,A,I1,A)') caseid0(1:i), '_3Dpln','00', iz, '.vtk'
    WRITE(fn_base,'(A,A,A,I1,A)') caseid0(1:i), '_3Dpln','00', iz
  ELSEIF(iz .LT.100) THEN
    WRITE(fn,'(A,A,A,I2,A)') caseid0(1:i), '_3Dpln','0', iz, '.vtk'
    WRITE(fn_base,'(A,A,A,I2,A)') caseid0(1:i), '_3Dpln','0', iz
  ELSE
    WRITE(fn,'(A,A,I3,A)') caseid0(1:i), '_3Dpln', iz, '.vtk'
    WRITE(fn_base,'(A,A,I3,A)') caseid0(1:i), '_3Dpln', iz
  ENDIF
  OPEN(unit=io9, file = fn, status = 'replace')
  CALL WriteGeomVTK3D(Core, iz)
  WRITE(io9, '(A, x, 2I10)') 'CELL_DATA', nGlobalelmt(iz)

  CALL WriteFlux(Core, FmInfo, fluxnorm0, ng, iz)
  CALL WritePower(Core, FmInfo, PowerDist, ThInfo, GroupInfo, iz, nTracerCntl, PE)
  IF(nTracerCntl%lFeedback) THEN
    CALL WritetEMP(Core, FmInfo, PowerDist, ThInfo, GroupInfo, iz, nTracerCntl)
    CALL WriteCoolTemp(Core, FmInfo, PowerDist, ThInfo, GroupInfo, iz, nTracerCntl)
  ENDIF

  CLOSE(io9)
ENDDO

WRITE(fn,'(A,A,A)') caseid0(1:i), '_3D','.visit'
OPEN(unit = io12, file=fn, status = 'replace')
WRITE(io12, '(A, 1x, I5)') '!NBLOCKS', COre%nz
DO iz = 1, CORE%nz
  fn = ''
  IF(iz .LT. 10) THEN
    WRITE(fn,'(A,A,A,I1,A)') caseid0(1:i), '_3Dpln','00', iz, '.vtk'
  ELSEIF(iz .LT.100) THEN
    WRITE(fn,'(A,A,A,I2,A)') caseid0(1:i), '_3Dpln','0', iz, '.vtk'
  ELSE
    WRITE(fn,'(A,A,I3,A)') caseid0(1:i), '_3Dpln', iz, '.vtk'
  ENDIF
  DO j = 1, 80
    if(fn(j:j) .EQ. '') EXIT
  ENDDO
  j = j-1
  WRITE(io12, '(A)') fn(1:j)
ENDDO
close(io12)
END SUBROUTINE

SUBROUTINE ProcessVTK2Dpln(Core, FMInfo, CmInfo, PowerDist, ThInfo, GroupInfo, nTracerCntl, PE)
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(PowerDist_Type) ::  PowerDist
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

TYPE(Cell_Type), POINTER :: Cell(:)
INTEGER :: nCell, ng
INTEGER :: i, iz,n

CHARACTER(256) :: fn, fn0,caseid0

REAL :: fluxnorm0

IF(PE%MASTER) PRINT *, 'WRITING VTK file...'

ng = GroupInfo%ng
Cell => Core%CellInfo
nCell = Core%nCellType

ALLOCATE(CellVTK(nCell))
DO i = 1, nCell
  IF(.NOT. Cell(i)%luse) CYCLE
  IF(.NOT. Cell(i)%lRect) THEN
    CALL MakeAnnualCellVTK(Cell(i), i)
    IF(Cell(i)%lCentX .OR. Cell(i)%lCentY .OR. Cell(i)%lCentXY) THEN
      CALL ProcessParialCell(Cell(i), i)
    ENDIF
    IF(Cell(i)%lCCell) THEN
      CALL ProcessParialCell(Cell(i), i)
    ENDIF

  ENDIF
  IF(Cell(i)%lRect) THEN
    CALL MakeRectCellVTK(Cell(i), i)

  ENDIF
  CALL BdCellVTK(Cell(i), i)
ENDDO

ALLOCATE(PlaneConfig(PE%myzb:PE%Myze))
DO iz = PE%myzb, PE%myze
  CALL SetPlaneConfig(Core, iz)
ENDDO

!Find case id
DO i = 1, 80
  if(caseid(i:i) .EQ. '') EXIT
ENDDO
i = i-1
caseid0(1:i) = caseid(1:i)
IF(nTracerCntl%CalNoId .LT. 10) THEN
  WRITE(caseid0(i+1:256),'(A, A, I1)') '_cycle', '00', nTracerCntl%CalNoId
ELSEIF(nTracerCntl%CalNoId .LT. 100) THEN
  WRITE(caseid0(i+1:256),'(A, A, I2)') '_cycle', '0', nTracerCntl%CalNoId
ELSE
  WRITE(caseid0(i+1:256),'(A, I3)') '_cycle', nTracerCntl%CalNoId
ENDIF
i=i+9

iz = 1
fluxnorm0 = FluxNorm(Core, FmInfo, CmInfo, ng, nTracerCntl, PE)
DO iz = PE%myzb, PE%myze
  IF(iz .LT. 10) THEN
    !WRITE(fn,'(A,A)') caseid0(1:i), '.vtk'
    !WRITE(fn_base,'(A,A)') caseid0(1:i)
    WRITE(fn,'(A,A,A,I1,A)') caseid0(1:i), '_pln','00', iz, '.vtk'
    WRITE(fn_base,'(A,A,A,I1,A)') caseid0(1:i), '_pln','00', iz
  ELSEIF(iz .LT.100) THEN
    WRITE(fn,'(A,A,A,I2,A)') caseid0(1:i), '_pln','0', iz, '.vtk'
    WRITE(fn_base,'(A,A,A,I2,A)') caseid0(1:i), '_pln','0', iz
  ELSE
    WRITE(fn,'(A,A,I3,A)') caseid0(1:i), '_pln', iz, '.vtk'
    WRITE(fn_base,'(A,A,I3,A)') caseid0(1:i), '_pln', iz
  ENDIF
  !WRITE(fn,'(A,A,A,I1,A)') caseid0(1:i), '.vtk'
  !WRITE(fn_base,'(A,A,A,I1,A)') caseid0(1:i)
  OPEN(unit=io9, file = fn, status = 'replace')

  CALL WriteGeomVTK2D(Core, iz)
  !CALL  WriteGeomVtk(Core, iz)
  !WRITE Cell Data
  WRITE(io9, '(A, x, 2I10)') 'CELL_DATA', nGlobalelmt(iz)

  CALL WriteFlux(Core, FmInfo, fluxnorm0, ng, iz)
  CALL WritePower(Core, FmInfo, PowerDist, ThInfo, GroupInfo, iz, nTracerCntl, PE)
  IF(nTracerCntl%lFeedback) THEN
    CALL WritetEMP(Core, FmInfo, PowerDist, ThInfo, GroupInfo, iz, nTracerCntl)
    CALL WriteCoolTemp(Core, FmInfo, PowerDist, ThInfo, GroupInfo, iz, nTracerCntl)
  ENDIF
  CALL Writeburnup(Core, FmInfo, PowerDist, ThInfo, GroupInfo, iz, nTracerCntl)
  CLOSE(io9)
ENDDO
END SUBROUTINE

SUBROUTINE ProcessVTK2DplnModel(Core, FMInfo,  GroupInfo, nTracerCntl, PE)
  TYPE(CoreInfo_Type) :: Core
  TYPE(FmInfo_Type) :: FmInfo
  TYPE(GroupInfo_Type) :: GroupInfo
  TYPE(nTracerCntl_Type) :: nTracerCntl
  TYPE(PE_Type) :: PE

  TYPE(Cell_Type), POINTER :: Cell(:)
  INTEGER :: nCell, ng
  INTEGER :: i, iz,n

  CHARACTER(256) :: fn, fn0,caseid0

  REAL :: fluxnorm0

  IF(PE%MASTER) PRINT *, 'WRITING VTK file...'

  ng = GroupInfo%ng
  Cell => Core%CellInfo
  nCell = Core%nCellType

  ALLOCATE(CellVTK(nCell))
  DO i = 1, nCell
    IF(.NOT. Cell(i)%luse) CYCLE
    IF(.NOT. Cell(i)%lRect) THEN
      CALL MakeAnnualCellVTK(Cell(i), i)
      IF(Cell(i)%lCentX .OR. Cell(i)%lCentY .OR. Cell(i)%lCentXY) THEN
        CALL ProcessParialCell(Cell(i), i)
      ENDIF
      IF(Cell(i)%lCCell) THEN
        CALL ProcessParialCell(Cell(i), i)
      ENDIF

    ENDIF
    IF(Cell(i)%lRect) THEN
      CALL MakeRectCellVTK(Cell(i), i)

    ENDIF
    CALL BdCellVTK(Cell(i), i)
  ENDDO

  ALLOCATE(PlaneConfig(PE%myzb:PE%Myze))
  DO iz = PE%myzb, PE%myze
    CALL SetPlaneConfig(Core, iz)
  ENDDO

  !Find case id
  DO i = 1, 80
    if(caseid(i:i) .EQ. '') EXIT
  ENDDO
  i = i-1
  caseid0(1:i) = caseid(1:i)
  IF(nTracerCntl%CalNoId .LT. 10) THEN
    WRITE(caseid0(i+1:256),'(A, A, I1)') '_cycle', '00', nTracerCntl%CalNoId
  ELSEIF(nTracerCntl%CalNoId .LT. 100) THEN
    WRITE(caseid0(i+1:256),'(A, A, I2)') '_cycle', '0', nTracerCntl%CalNoId
  ELSE
    WRITE(caseid0(i+1:256),'(A, I3)') '_cycle', nTracerCntl%CalNoId
  ENDIF
  i=i+9

  iz = 1
  DO iz = PE%myzb, PE%myze
    IF(iz .LT. 10) THEN
      !WRITE(fn,'(A,A)') caseid0(1:i), '.vtk'
      !WRITE(fn_base,'(A,A)') caseid0(1:i)
      WRITE(fn,'(A,A,A,I1,A)') caseid0(1:i), '_pln','00', iz, '.vtk'
      WRITE(fn_base,'(A,A,A,I1,A)') caseid0(1:i), '_pln','00', iz
    ELSEIF(iz .LT.100) THEN
      WRITE(fn,'(A,A,A,I2,A)') caseid0(1:i), '_pln','0', iz, '.vtk'
      WRITE(fn_base,'(A,A,A,I2,A)') caseid0(1:i), '_pln','0', iz
    ELSE
      WRITE(fn,'(A,A,I3,A)') caseid0(1:i), '_pln', iz, '.vtk'
      WRITE(fn_base,'(A,A,I3,A)') caseid0(1:i), '_pln', iz
    ENDIF
    !WRITE(fn,'(A,A,A,I1,A)') caseid0(1:i), '.vtk'
    !WRITE(fn_base,'(A,A,A,I1,A)') caseid0(1:i)
    OPEN(unit=io9, file = fn, status = 'replace')

    CALL WriteGeomVTK2D(Core, iz)
    !CALL  WriteGeomVtk(Core, iz)
    !WRITE Cell Data
    WRITE(io9, '(A, x, 2I10)') 'CELL_DATA', nGlobalelmt(iz)

    CALL Writematerial(Core ,FmInfo , iz, nTracerCntl)

    CLOSE(io9)
  ENDDO
  END SUBROUTINE



SUBROUTINE SetPlaneConfig(Core, iz)
USE PARAM
USE ALLOCS
USE BASICOPERATION, ONLY : CP_CA
USE UtilFunction,   ONLY : Array2DSORT
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
INTEGER :: iz

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
REAL, POINTER :: CX(:), CY(:)
INTEGER, POINTER ::  NODE(:, :)

INTEGER :: nx, ny, nxy, npt, m, m1, m2,  n, l
INTEGER :: i, j, k
INTEGER :: ixy, ixy0, ixy1, ixy2, ix, iy, icel, icel1, icel2, isurf1, isurf2
REAL :: pt(2, 100)
REAL :: y1, y2, x1, x2


nxy = Core%nxy; nx = Core%nx; ny = Core%ny
Pin => Core%Pin; Cell => Core%CellInfo
ALLOCATE(PlaneConfig(iz)%NODE(0:nx+1, 0:ny+1))
ALLOCATE(PlaneConfig(iz)%CX(nx+1))
ALLOCATE(PlaneConfig(iz)%CY(ny+1))

CX => PlaneConfig(iz)%CX
CY => PlaneConfig(iz)%CY
NODE => PlaneConfig(iz)%NODE
CALL CP_CA(node(0:nx+1, 0:ny+1), 0, nx+2, ny+2)
CALL CP_CA(CX, 0._8, nx+1)
CALL CP_CA(CY, 0._8, ny+1)
DO ixy = 1, Core%nxy
  ix = Pin(ixy)%ix; iy = Pin(ixy)%iy
  node(ix, iy) = ixy
ENDDO
DO iy = 1, ny
  IF(NODE(1, iy) .NE.0 .AND. NODE(nx, iy) .NE.0) EXIT
ENDDO
DO ix = 1, nx
  ixy = NODE(ix, iy)
  icel = Pin(ixy)%Cell(iz)
  IF(Cell(icel)%lCentXY) THEN
    CX(ix) = 0
    CX(ix+1) = Cell(icel)%geom%lx
  ELSE
    CX(ix)= CX(ix) + Cell(icel)%geom%lx/2._8
    CX(ix+1) = CX(ix) + Cell(icel)%geom%lx/2._8
  ENDIF
ENDDO

DO ix = 1, nx
  IF(NODE(ix, 1) .NE.0 .AND. NODE(ix, ny) .NE.0) EXIT
ENDDO
DO iy = 1, ny
  ixy = NODE(ix, iy)
  icel = Pin(ixy)%Cell(iz)
  IF(Cell(icel)%lCentXY) THEN
    CY(iy) = 0
    CY(iy+1) = - Cell(icel)%geom%ly
  ELSE
    CY(iy)= CY(iy) - Cell(icel)%geom%ly/2._8
    CY(iy+1) = CY(iy) - Cell(icel)%geom%ly/2._8
  ENDIF

ENDDO

m = 0
DO icel = 1, Core%nCellType
  DO j = 1, 4
    m = max(m, CellVTK(icel)%nbddel(j)+1)
  ENDDO
ENDDO
CALL Dmalloc0(PlaneConfig(iz)%CellBdPt, 1, 2, 1, 2*m, 1, 4, 1, nxy)
CALL Dmalloc0(PlaneConfig(iz)%PtIdx, 0, 2*m, 1, 4, 1, nxy)

npt = 0

DO iy = 0, ny
  isurf1 = 1;   isurf2 = 3
  k = 0
  DO ix = 1, nx
    ixy1 = NODE(ix, iy); ixy2 = NODE(ix, iy +1)
#ifndef CENT_OPT
    IF(ixy1 .NE. 0) THEN
      icel1 = Pin(ixy1)%Cell(iz)
      IF(Cell(icel1)%lCentX) ixy1=0
      IF(Cell(icel1)%lCentY) ixy1=0
      IF(Cell(icel1)%lCentXY) ixy1=0
    ENDIF
    IF(ixy2 .NE. 0) THEN
      icel2 = Pin(ixy2)%Cell(iz)
      IF(Cell(icel2)%lCentX) ixy2=0
      IF(Cell(icel2)%lCentY) ixy2=0
      IF(Cell(icel2)%lCentXY) ixy2=0
    ENDIF
#endif
    IF(ixy1 .EQ. 0 .AND. ixy2 .EQ. 0) CYCLE
    icel1 = 0; icel2 = 0; m = 0

    IF(ixy1 .NE. 0) THEN
      icel1 = Pin(ixy1)%Cell(iz)
      m1 = CellVTK(icel1)%nbddel(isurf1)+1
      m = m + m1
      pt(1:2, 1:m) = CellVTK(icel1)%BdPtList(1:2, 1:m, isurf1)
      y1= pt(2, 1)
    ENDIF
    IF(ixy2 .NE. 0) THEN
      icel2 = Pin(ixy2)%Cell(iz)
      m2 = CellVTK(icel2)%nbddel(isurf2)+1
      pt(1:2, m+1 : m + m2) = CellVTK(icel2)%BdPtList(1:2, 1:m2, isurf2)
      y2= pt(2, m+1)
      m = m + m2
    ENDIF
    pt(2, 1:m) = 0
    CALL Array2DSORT(pt(1, 1:m), pt(2, 1:m), m, n, .TRUE., .TRUE., 1)
    k = k + 1
    IF(ixy1 .NE. 0) THEN
      PlaneConfig(iz)%CellBdPT(1:2, 1, isurf1, ixy1) = (/Pt(1, 1), y1/)
      PlaneConfig(iz)%PtIdx(0, isurf1, ixy1) = n
    ENDIF
    IF(ixy2 .NE. 0) THEN
      PlaneConfig(iz)%CellBdPT(1:2, 1, isurf2, ixy2) = (/Pt(1, 1), y2/)
      PlaneConfig(iz)%PtIdx(0, isurf2, ixy2) = n
    ENDIF

    IF(k .EQ. 1) npt = npt + 1
    IF(ixy1 .NE. 0) PlaneConfig(iz)%PtIdx(1, isurf1, ixy1) = npt
    IF(ixy2 .NE. 0) PlaneConfig(iz)%PtIdx(1, isurf2, ixy2) = npt

    DO i = 2, n
      k = k + 1
      npt = npt + 1
      IF(ixy1 .NE. 0) THEN
        PlaneConfig(iz)%CellBdPT(1:2, i, isurf1, ixy1) = (/Pt(1, i), y1/)
        PlaneConfig(iz)%PtIdx(i, isurf1, ixy1) = npt
      ENDIF
      IF(ixy2 .NE. 0) THEN
        PlaneConfig(iz)%CellBdPT(1:2, i, isurf2, ixy2) = (/Pt(1, i), y2/)
        PlaneConfig(iz)%PtIdx(i, isurf2, ixy2) = npt
      ENDIF
    ENDDO
  ENDDO
ENDDO

DO ix = 0, nx
  isurf1 = 4;   isurf2 = 2
  k = 0
  DO iy = 1, ny
    ixy1 = NODE(ix, iy); ixy2 = NODE(ix + 1, iy)
#ifndef CENT_OPT
    IF(ixy1 .NE. 0) THEN
      icel1 = Pin(ixy1)%Cell(iz)
      IF(Cell(icel1)%lCentX) ixy1=0
      IF(Cell(icel1)%lCentY) ixy1=0
      IF(Cell(icel1)%lCentXY) ixy1=0
    ENDIF
    IF(ixy2 .NE. 0) THEN
      icel2 = Pin(ixy2)%Cell(iz)
      IF(Cell(icel2)%lCentX) ixy2=0
      IF(Cell(icel2)%lCentY) ixy2=0
      IF(Cell(icel2)%lCentXY) ixy2=0
    ENDIF
#endif
    IF(ixy1 .EQ. 0 .AND. ixy2 .EQ. 0) CYCLE

    icel1 = 0; icel2 = 0; m = 0
    IF(ixy1 .NE. 0) THEN
      icel1 = Pin(ixy1)%Cell(iz)
      m1 = CellVTK(icel1)%nbddel(isurf1)+1
      m = m + m1
      pt(1:2, 1:m) = CellVTK(icel1)%BdPtList(1:2, 1:m, isurf1)
      x1= pt(1, 1)
    ENDIF
    IF(ixy2 .NE. 0) THEN
      icel2 = Pin(ixy2)%Cell(iz)
      m2 = CellVTK(icel2)%nbddel(isurf2)+1
      pt(1:2, m+1 : m + m2) = CellVTK(icel2)%BdPtList(1:2, 1:m2, isurf2)
      x2= pt(1, m+1)
      m = m + m2
    ENDIF
    pt(1, 1:m) = 0
    CALL Array2DSORT(pt(1, 1:m), pt(2, 1:m), m, n, .FALSE., .TRUE., 2)

    IF(ixy1 .NE. 0) THEN
      l = PlaneConfig(iz)%PtIdx(0, 3, ixy1); l = PlaneConfig(iz)%PtIdx(l, 3, ixy1)
    ELSE
      l = PlaneConfig(iz)%PtIdx(1, 3, ixy2)
    ENDIF

    IF(ixy1 .NE. 0) THEN
      PlaneConfig(iz)%CellBdPT(1:2, 1, isurf1, ixy1) = (/x1, Pt(2, 1)/)
      PlaneConfig(iz)%PtIdx(0, isurf1, ixy1) = n
      PlaneConfig(iz)%PtIdx(1, isurf1, ixy1) = l
    ENDIF
    IF(ixy2 .NE. 0) THEN
      PlaneConfig(iz)%CellBdPT(1:2, 1, isurf2, ixy2) = (/x2, Pt(2, 1)/)
      PlaneConfig(iz)%PtIdx(0, isurf2, ixy2) = n
      PlaneConfig(iz)%PtIdx(1, isurf2, ixy2) = l
    ENDIF

    DO i = 2, n-1
      k = k + 1
      npt = npt + 1
      IF(ixy1 .NE. 0) THEN
        PlaneConfig(iz)%CellBdPT(1:2, i, isurf1, ixy1) = (/x1, Pt(2, i)/)
        PlaneConfig(iz)%PtIdx(i, isurf1, ixy1) = npt
      ENDIF
      IF(ixy2 .NE. 0) THEN
        PlaneConfig(iz)%CellBdPT(1:2, i, isurf2, ixy2) = (/x2, Pt(2, i)/)
        PlaneConfig(iz)%PtIdx(i, isurf2, ixy2) = npt
      ENDIF
    ENDDO

    i = n
    IF(ixy1 .NE. 0) THEN
      l = PlaneConfig(iz)%PtIdx(0, 1, ixy1);l = PlaneConfig(iz)%PtIdx(l, 1, ixy1)
    ELSE
      l = PlaneConfig(iz)%PtIdx(1, 1, ixy2);
    ENDIF
    IF(ixy1 .NE. 0) THEN
      PlaneConfig(iz)%CellBdPT(1:2, i, isurf1, ixy1) = (/x1, Pt(2, i)/)
      PlaneConfig(iz)%PtIdx(i, isurf1, ixy1) = l
    ENDIF
    IF(ixy2 .NE. 0) THEN
      PlaneConfig(iz)%CellBdPT(1:2, i, isurf2, ixy2) = (/x2, Pt(2, i)/)
      PlaneConfig(iz)%PtIdx(i, isurf2, ixy2) = l
    ENDIF
    CONTINUE
  ENDDO
ENDDO
PlaneConfig(iz)%npt = npt
NULLIFY(CX, CY, NODE)
END SUBROUTINE

SUBROUTINE WriteGeomVTK3D_HOM(Core, myzb, myze)
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
INTEGER :: myzb, myze

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)

INTEGER :: nxy, nx, ny, nz
INTEGER :: ix, iy, iz, ixy, icel
INTEGER, POINTER :: node(:, :)
REAL, POINTER :: x(:), Y(:), Z(:)

nxy = Core%nxy
nx = Core%nx; ny = Core%ny; nz = Core%nz
Pin => CORE%Pin
Cell => Core%CellInfo

ALLOCATE(NODE(nx, ny))
ALLOCATE(X(0:nx)); ALLOCATE(Y(0:ny)); ALLOCATE(Z(0:nz))
CALL CP_CA(node, 0, nx, ny)
CALL CP_CA(X(0:nx), 0._8, nx+1)
CALL CP_CA(Y(0:nx), 0._8, ny+1)
DO ixy = 1, Core%nxy
  ix = Pin(ixy)%ix; iy = Pin(ixy)%iy
  node(ix, iy) = ixy
ENDDO
DO iy = 1, ny
  IF(NODE(1, iy) .NE.0 .AND. NODE(nx, iy) .NE.0) EXIT
ENDDO
iz = myzb
X(0) = 0
DO ix = 1, nx
  ixy = NODE(ix, iy)
  icel = Pin(ixy)%Cell(iz)
  X(ix) = X(ix-1) + Cell(icel)%Geom%lx
ENDDO

DO ix = 1, nx
  IF(NODE(ix, 1) .NE.0 .AND. NODE(ix, ny) .NE.0) EXIT
ENDDO
Y(0) = 0
DO iy = 1, ny
  ixy = NODE(ix, iy)
  icel = Pin(ixy)%Cell(iz)
  Y(iy) = Y(iy-1) - Cell(icel)%Geom%ly
ENDDO
Z(0) = 0
DO iz = 1, nz
  Z(iz) = Z(iz-1) + Core%hz(iz)
ENDDO

WRITE(io9, '(A)') '# vtk DataFile Version 3.0'
WRITE(io9, '(A)') 'nTRACER Visualization'
WRITE(io9, '(A)') 'ASCII'
WRITE(io9, '(A)') 'DATASET RECTILINEAR_GRID'
WRITE(io9, '(A, 3I7)') 'DIMENSIONS', nx + 1, ny + 1, myze - myzb + 2
WRITE(io9, '(A, I7, 2x, A)') 'X_COORDINATES', nx + 1, 'FLOAT'
WRITE(io9, '(7F12.5)') (x(ix), ix = 0, nx)
WRITE(io9, '(A, I7, 2x, A)') 'Y_COORDINATES', ny + 1, 'FLOAT'
WRITE(io9, '(7F12.5)') (y(iy), iy = 0, ny)
WRITE(io9, '(A, I7, 2x, A)') 'Z_COORDINATES', myze - myzb + 2, 'FLOAT'
WRITE(io9, '(7F12.5)') (z(iz), iz = myzb-1, myze)
DEALLOCATE(x, y, node)

END SUBROUTINE


SUBROUTINE WriteGeomVTK3D_HOM2(Core, myzb, myze)
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
INTEGER :: myzb, myze

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)

INTEGER :: nxy, nx, ny, nz, npt, ndat, ncell
INTEGER :: ix, iy, iz, ixy, icel
INTEGER, POINTER :: node(:, :)
REAL, POINTER :: x(:), Y(:), Z(:)
INTEGER, POINTER :: NodeMap(:, :, :)

nxy = Core%nxy
nx = Core%nx; ny = Core%ny; nz = Core%nz
Pin => CORE%Pin
Cell => Core%CellInfo

ALLOCATE(NODE(nx, ny))
ALLOCATE(X(0:nx)); ALLOCATE(Y(0:ny)); ALLOCATE(Z(0:nz))
ALLOCATE(NodeMap(0:nx, 0:ny, myzb-1:myze))
CALL CP_CA(node, 0, nx, ny)
CALL CP_CA(X(0:nx), 0._8, nx+1)
CALL CP_CA(Y(0:nx), 0._8, ny+1)
DO ixy = 1, Core%nxy
  ix = Pin(ixy)%ix; iy = Pin(ixy)%iy
  node(ix, iy) = ixy
ENDDO
DO iy = 1, ny
  IF(NODE(1, iy) .NE.0 .AND. NODE(nx, iy) .NE.0) EXIT
ENDDO
iz = myzb
X(0) = 0
DO ix = 1, nx
  ixy = NODE(ix, iy)
  icel = Pin(ixy)%Cell(iz)
  X(ix) = X(ix-1) + Cell(icel)%Geom%lx
ENDDO

DO ix = 1, nx
  IF(NODE(ix, 1) .NE.0 .AND. NODE(ix, ny) .NE.0) EXIT
ENDDO
Y(0) = 0
DO iy = 1, ny
  ixy = NODE(ix, iy)
  icel = Pin(ixy)%Cell(iz)
  Y(iy) = Y(iy-1) - Cell(icel)%Geom%ly
ENDDO
Z(0) = 0
DO iz = 1, nz
  Z(iz) = Z(iz-1) + Core%hz(iz)
ENDDO
npt = 0
DO iz = myzb - 1, myze
  DO iy = 0, ny
    DO ix = 0, nx
      npt = npt+1
      NodeMap(ix, iy, iz) = npt-1
    ENDDO
  ENDDO
ENDDO

WRITE(io9, '(A)') '# vtk DataFile Version 3.0'
WRITE(io9, '(A)') '2D Unstructured Grid of Linear Triangles'
WRITE(io9, '(A)') 'ASCII'
WRITE(io9, '(A)') 'DATASET UNSTRUCTURED_GRID'
WRITE(io9, '(A, x, I10, 2x, A)') 'POINTS', npt, 'float'
DO iz = myzb - 1, myze
  DO iy = 0, ny
    DO ix = 0, nx
      WRITE(io9, '(3F25.7)') x(ix), y(iy), z(iz)
    ENDDO
  ENDDO
ENDDO
!WRITE(io9, *)

ncell = nxy * (myze-myzb+1)
ndat = ncell* 9
WRITE(io9, '(A, x, 2I15)') 'CELLS', ncell, ndat
DO iz = myzb, myze
  DO ixy = 1, nxy
    ix = Pin(ixy)%ix; iy = Pin(ixy)%iy
    WRITE(io9, '(50I8)') 8, NodeMap(ix-1, iy, iz-1), NodeMap(ix, iy, iz-1), NodeMap(ix-1, iy-1, iz-1), NodeMap(ix, iy-1, iz-1),   &
                            NodeMap(ix-1, iy, iz  ), NodeMap(ix, iy, iz  ), NodeMap(ix-1, iy-1, iz), NodeMap(ix, iy-1, iz  )
  ENDDO
ENDDO
 WRITE(io9, '(50I8)')
WRITE(io9, '(A, x, 2I10)') 'CELL_TYPES', ncell
DO ixy = 1, ncell

   WRITE(io9, '(50I8)') 11
ENDDO
 WRITE(io9, '(50I8)')
!WRITE(io9, '(A, x, 2I10)') 'CELLS', nGlobalElmt(iz), ndat
!
!n = 0
!DO ixy = 1, nxy
!  ix = Pin(ixy)%ix; iy = Pin(ixy)%iy
!  icel = Pin(ixy)%Cell(iz)
!  DO i = 1, CellVTK3D(icel)%nelmt
!    m = CellVTK3D(icel)%Elements(0, i)
!    WRITE(io9, '(50I8)') CellVTK3D(icel)%elements(0, i), CellVTK3D(icel)%elements(1:m, i) + n
!  ENDDO
!  n = n + CellVTK3D(icel)%npt
!ENDDO
!
!WRITE(io9, '(A, x, 2I10)') 'CELL_TYPES', nGlobalelmt(iz)
!DO ixy = 1, nxy
!  ix = Pin(ixy)%ix; iy = Pin(ixy)%iy
!  icel = Pin(ixy)%Cell(iz)
!  DO i = 1, CellVTK3D(icel)%nelmt
!    WRITE(io9, '(50I8)') CellVTK3D(icel)%types(i)
!  ENDDO
!ENDDO
!WRITE(io9, *)
!!  m = CellVTK3D(icel)%Elements(0, i)
!!  WRITE(iunt, '(50I8)') CellVTK3D(icel)%elements(0:m, i)
!
!!DO i = 1, CellVTK3D(icel)%nelmt
!!  m = CellVTK3D(icel)%Elements(0, i)
!!  WRITE(iunt, '(50I8)') CellVTK3D(icel)%elements(0:m, i)
!!ENDDO


!!!!!!


!WRITE(io9, '(A)') '# vtk DataFile Version 3.0'
!WRITE(io9, '(A)') 'nTRACER Visualization'
!WRITE(io9, '(A)') 'ASCII'
!WRITE(io9, '(A)') 'DATASET RECTILINEAR_GRID'
!WRITE(io9, '(A, 3I7)') 'DIMENSIONS', nx + 1, ny + 1, myze - myzb + 2
!WRITE(io9, '(A, I7, 2x, A)') 'X_COORDINATES', nx + 1, 'FLOAT'
!WRITE(io9, '(7F12.5)') (x(ix), ix = 0, nx)
!WRITE(io9, '(A, I7, 2x, A)') 'Y_COORDINATES', ny + 1, 'FLOAT'
!WRITE(io9, '(7F12.5)') (y(iy), iy = 0, ny)
!WRITE(io9, '(A, I7, 2x, A)') 'Z_COORDINATES', myze - myzb + 2, 'FLOAT'
!WRITE(io9, '(7F12.5)') (z(iz), iz = myzb-1, myze)
!DEALLOCATE(x, y, node)

END SUBROUTINE

SUBROUTINE WriteGeomVTK3D(Core,iz)
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
INTEGER :: iz

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
INTEGER, POINTER :: NODE(:, :)
INTEGER :: nxy, ndat, nx, ny, m, n

INTEGER :: ixy, icel, i, ix, iy
REAL :: CX0, CY0, Hz1, HZ2
REAL, POINTER :: CX(:), CY(:)

nxy = Core%nxy
nx = Core%nx; ny = Core%ny
Pin => CORE%Pin
Cell => Core%CellInfo

hz1= 0
DO i = 1, iz-1
  hz1 = hz1 + Core%hz(i)
ENDDO
hz2 = hz1 + Core%hz(iz)

nGlobalpt(iz) = 0
nGlobalElmt(iz) = 0
ndat = 0

DO ixy = 1, nxy
  icel = Pin(ixy)%cell(iz)
  nGlobalpt(iz) = nGlobalpt(iz) + CellVTK3D(icel)%npt
  nGlobalElmt(iz) =nGlobalElmt(iz) + CellVTK3D(icel)%nelmt
  DO i = 1, CellVTK3D(icel)%nelmt
    ndat = ndat + (CellVTK3D(icel)%Elements(0, i) + 1)
  ENDDO
ENDDO

ALLOCATE(NODE(nx, ny))
ALLOCATE(CX(nx+1))
ALLOCATE(CY(ny+1))
CALL CP_CA(node, 0, nx, ny)
CALL CP_CA(CX, 0._8, nx+1)
CALL CP_CA(CY, 0._8, ny+1)
DO ixy = 1, Core%nxy
  ix = Pin(ixy)%ix; iy = Pin(ixy)%iy
  node(ix, iy) = ixy
ENDDO
DO iy = 1, ny
  IF(NODE(1, iy) .NE.0 .AND. NODE(nx, iy) .NE.0) EXIT
ENDDO
DO ix = 1, nx
  ixy = NODE(ix, iy)
  icel = Pin(ixy)%Cell(iz)
  IF(Cell(icel)%lCentXY) THEN
    CX(ix) = 0
    CX(ix+1) = Cell(icel)%geom%lx
  ELSE
    CX(ix)= CX(ix) + Cell(icel)%geom%lx/2._8
    CX(ix+1) = CX(ix) + Cell(icel)%geom%lx/2._8
  ENDIF
ENDDO

DO ix = 1, nx
  IF(NODE(ix, 1) .NE.0 .AND. NODE(ix, ny) .NE.0) EXIT
ENDDO
DO iy = 1, ny
  ixy = NODE(ix, iy)
  icel = Pin(ixy)%Cell(iz)
  IF(Cell(icel)%lCentXY) THEN
    CY(iy) = 0
    CY(iy+1) = -Cell(icel)%geom%ly
  ELSE
    CY(iy)= CY(iy) - Cell(icel)%geom%ly/2._8
    CY(iy+1) = CY(iy) - Cell(icel)%geom%ly/2._8
  ENDIF
ENDDO
WRITE(io9, '(A)') '# vtk DataFile Version 3.0'
WRITE(io9, '(A)') '2D Unstructured Grid of Linear Triangles'
WRITE(io9, '(A)') 'ASCII'
WRITE(io9, '(A)') 'DATASET UNSTRUCTURED_GRID'
WRITE(io9, '(A, x, I10, 2x, A)') 'POINTS', nGlobalPt(iz), 'float'
DO ixy = 1, nxy
  ix = Pin(ixy)%ix; iy = Pin(ixy)%iy
  CX0 = CX(ix); CY0 = CY(iy)
  icel = Pin(ixy)%Cell(iz)
  n = CellVTK3D(icel)%npt/2
  DO i = 1, n
    WRITE(io9, '(3F25.7)') CellVTK3D(icel)%pt(1, i) + CX0, CellVTK3D(icel)%pt(2, i) + CY0, hz1
  ENDDO
  DO i = n+1, CellVTK3D(icel)%npt
    WRITE(io9, '(3F25.7)') CellVTK3D(icel)%pt(1, i) + CX0, CellVTK3D(icel)%pt(2, i) + CY0, hz2
  ENDDO
ENDDO

WRITE(io9, '(A, x, 2I10)') 'CELLS', nGlobalElmt(iz), ndat

n = 0
DO ixy = 1, nxy
  ix = Pin(ixy)%ix; iy = Pin(ixy)%iy
  icel = Pin(ixy)%Cell(iz)
  DO i = 1, CellVTK3D(icel)%nelmt
    m = CellVTK3D(icel)%Elements(0, i)
    WRITE(io9, '(50I8)') CellVTK3D(icel)%elements(0, i), CellVTK3D(icel)%elements(1:m, i) + n
  ENDDO
  n = n + CellVTK3D(icel)%npt
ENDDO

WRITE(io9, '(A, x, 2I10)') 'CELL_TYPES', nGlobalelmt(iz)
DO ixy = 1, nxy
  ix = Pin(ixy)%ix; iy = Pin(ixy)%iy
  icel = Pin(ixy)%Cell(iz)
  DO i = 1, CellVTK3D(icel)%nelmt
    WRITE(io9, '(50I8)') CellVTK3D(icel)%types(i)
  ENDDO
ENDDO
WRITE(io9, *)
!  m = CellVTK3D(icel)%Elements(0, i)
!  WRITE(iunt, '(50I8)') CellVTK3D(icel)%elements(0:m, i)

!DO i = 1, CellVTK3D(icel)%nelmt
!  m = CellVTK3D(icel)%Elements(0, i)
!  WRITE(iunt, '(50I8)') CellVTK3D(icel)%elements(0:m, i)
!ENDDO
DEALLOCATE(CX, CY, node)

END SUBROUTINE

SUBROUTINE WriteGeomVTK2D(Core, iz)
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
INTEGER :: iz

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)

INTEGER :: nxy, nx, ny, npt, ndat, nelmt
INTEGER :: idum, ixy, icel, iElmtSt
INTEGER :: Dat(0:100)
CHARACTER(256) :: fn1, fn2, oneline
nxy = Core%nxy
nx = Core%nx; ny = Core%ny
Pin => CORE%Pin
Cell => Core%CellInfo

fn1 = trim(fn_base)//'.vtktmp1'
OPEN(unit=io10, file = fn1, status = 'replace')
fn2 = trim(fn_base)//'.vtktmp2'
OPEN(unit=io11, file = fn2, status = 'replace')

WRITE(io9, '(A)') '# vtk DataFile Version 3.0'
WRITE(io9, '(A)') '2D Unstructured Grid of Linear Triangles'
WRITE(io9, '(A)') 'ASCII'
WRITE(io9, '(A)') 'DATASET UNSTRUCTURED_GRID'
CALL ProcessPoints2D(Core, iz, npt)
iElmtSt = PlaneConfig(iz)%npt
ndat = 0
nelmt = 0
DO ixy = 1, nxy
  icel = Pin(ixy)%Cell(iz)
#ifndef CENT_OPT
  IF(Cell(icel)%lCentX) CYCLE
  IF(Cell(icel)%lCentY) CYCLE
  IF(Cell(icel)%lCentXY) CYCLE
#endif
  Call ProcessCell2D(Core, iz, ixy, iElmtSt, ndat)
  nelmt = nelmt + CellVTK(Pin(ixy)%Cell(iz))%nelmt
ENDDO
CLOSE(io10)
CLOSE(io11)
nGlobalelmt(iz) = nelmt
!WRITE(io9, '(A, x, 2I10)') 'CELLS', nGlobalElmt(iz), ndat
!
!n = 0
!DO ixy = 1, nxy
!  ix = Pin(ixy)%ix; iy = Pin(ixy)%iy
!  icel = Pin(ixy)%Cell(iz)
!  DO i = 1, CellVTK(icel)%nelmt
!    m = CellVTK(icel)%Elements(0, i)
!    WRITE(io9, '(50I8)') CellVtk(icel)%elements(0, i), CellVtk(icel)%elements(1:m, i) + n
!  ENDDO
!  n = n + CellVTK(icel)%npt
!ENDDO
!
!WRITE(io9, '(A, x, 2I10)') 'CELL_TYPES', nGlobalelmt(iz)
!DO ixy = 1, nxy
!  ix = Pin(ixy)%ix; iy = Pin(ixy)%iy
!  icel = Pin(ixy)%Cell(iz)
!  DO i = 1, CellVTK(icel)%nelmt
!    WRITE(io9, '(50I8)') CellVtk(icel)%types(i)
!  ENDDO
!ENDDO
!WRITE(io9, *)
!  m = CellVt
OPEN(unit=io10, file = fn1, status = 'OLD')
OPEN(unit=io11, file = fn2, status = 'OLD')
WRITE(io9, '(A, x, 2I10)') 'CELLS', nelmt, ndat
DO
  READ(io10, '(A256)', END = 10) oneline
  READ(oneline, *)  dat(0)
  READ(oneline, *) idum, dat(1:dat(0))
  WRITE(io9, '(50I8)') dat(0), dat(1:dat(0))
ENDDO
10 continue
WRITE(io9, '(A, x, 2I10)') 'CELL_TYPES', nelmt
DO
  READ(io11, *, END = 11)  dat(0)
  WRITE(io9, '(50I8)') dat(0)
ENDDO
11 continue
CLOSE(io10, STATUS='DELETE')
CLOSE(io11, STATUS='DELETE')
END SUBROUTINE

SUBROUTINE ProcessPoints2D(Core, iz, npt)
TYPE(CoreInfo_Type) :: Core
INTEGER :: iz, npt

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)

INTEGER :: nx, ny, nxy
INTEGER :: ixy, ixy1, ixy2, isurf1, isurf2, ix, iy, ipt, icel, icel1, icel2
INTEGER :: i, j
INTEGER :: m, n0, n
REAL :: CX0, CY0, X, Y, hz


nxy = Core%nxy
nx = Core%nx; ny = Core%ny
Pin => CORE%Pin
Cell => Core%CellInfo

npt = PlaneConfig(iz)%npt
hz = Core%hz(1)/2._8
DO i = 2, iz
  hz = hz + Core%hz(i-1)/2._8 +  Core%hz(i)
ENDDO
hz = 0
DO ixy = 1, nxy
  icel = Pin(ixy)%Cell(iz)
#ifndef CENT_OPT
  IF(Cell(icel)%lCentX) CYCLE
  IF(Cell(icel)%lCentY) CYCLE
  IF(Cell(icel)%lCentXY) CYCLE
#endif
  npt = npt + CellVtk(icel)%ninpt
ENDDO
nGlobalPt(IZ) = NPT
WRITE(io9, '(A, x, I10, 2x, A)') 'POINTS', nGlobalPt(IZ), 'float'
n = 0
DO iy = 0, ny
  DO ix = 1, nx
    isurf1 = 1; isurf2 = 3
    ixy1 = PlaneConfig(iz)%NODE(ix, iy); ixy2 = PlaneConfig(iz)%NODE(ix, iy+1)
#ifndef CENT_OPT
    IF(ixy1 .NE. 0) THEN
      icel1 = Pin(ixy1)%Cell(iz)
      IF(Cell(icel1)%lCentX) ixy1=0
      IF(Cell(icel1)%lCentY) ixy1=0
      IF(Cell(icel1)%lCentXY) ixy1=0
    ENDIF
    IF(ixy2 .NE. 0) THEN
      icel2 = Pin(ixy2)%Cell(iz)
      IF(Cell(icel2)%lCentX) ixy2=0
      IF(Cell(icel2)%lCentY) ixy2=0
      IF(Cell(icel2)%lCentXY) ixy2=0
    ENDIF
#endif
    IF(ixy1 .EQ. 0 .AND. ixy2 .EQ. 0) CYCLE
    icel1 = 0; icel2 = 0;
    IF(ixy1 .NE. 0) THEN
      icel1 = Pin(ixy1)%Cell(iz)
      m = PlaneConfig(iz)%PtIdx(0, isurf1, ixy1)
      CX0 = PlaneConfig(iz)%CX(ix); CY0 = PlaneConfig(iz)%CY(iy)
      DO i = 1, m
        n0 = PlaneConfig(iz)%PtIdx(i, isurf1, ixy1)
        IF(n .LT. n0) THEN
          X = CX0 + PlaneConfig(iz)%CellBdPt(1, i, isurf1, ixy1)
          Y = CY0 + PlaneConfig(iz)%CellBdPt(2, i, isurf1, ixy1)
          WRITE(io9, '(3F25.7)') X, Y, hz
          n = n0
        ENDIF
      ENDDO
    ELSE
      icel2 = Pin(ixy2)%Cell(iz)
      m = PlaneConfig(iz)%PtIdx(0, isurf2, ixy2)
      CX0 = PlaneConfig(iz)%CX(ix); CY0 = PlaneConfig(iz)%CY(iy+1)
      DO i = 1, m
        n0 = PlaneConfig(iz)%PtIdx(i, isurf2, ixy2)
        IF(n .LT. n0) THEN
          X = CX0 + PlaneConfig(iz)%CellBdPt(1, i, isurf2, ixy2)
          Y = CY0 + PlaneConfig(iz)%CellBdPt(2, i, isurf2, ixy2)
          WRITE(io9, '(3F25.7)') X, Y, hz
          n = n0
        ENDIF
      ENDDO
    ENDIF
  ENDDO
ENDDO
DO ix = 0, nx
  DO iy = 1, ny
    isurf1 = 4; isurf2 = 2
    ixy1 = PlaneConfig(iz)%NODE(ix, iy); ixy2 = PlaneConfig(iz)%NODE(ix+1, iy)
#ifndef CENT_OPT
    IF(ixy1 .NE. 0) THEN
      icel1 = Pin(ixy1)%Cell(iz)
      IF(Cell(icel1)%lCentX) ixy1=0
      IF(Cell(icel1)%lCentY) ixy1=0
      IF(Cell(icel1)%lCentXY) ixy1=0
    ENDIF
    IF(ixy2 .NE. 0) THEN
      icel2 = Pin(ixy2)%Cell(iz)
      IF(Cell(icel2)%lCentX) ixy2=0
      IF(Cell(icel2)%lCentY) ixy2=0
      IF(Cell(icel2)%lCentXY) ixy2=0
    ENDIF
#endif
    IF(ixy1 .EQ. 0 .AND. ixy2 .EQ. 0) CYCLE
    icel1 = 0; icel2 = 0;
    IF(ixy1 .NE. 0) THEN
      icel1 = Pin(ixy1)%Cell(iz)
      m = PlaneConfig(iz)%PtIdx(0, isurf1, ixy1)
      CX0 = PlaneConfig(iz)%CX(ix); CY0 = PlaneConfig(iz)%CY(iy)
      DO i = 1, m
        n0 = PlaneConfig(iz)%PtIdx(i, isurf1, ixy1)
        IF(n .LT. n0) THEN
          X = CX0 + PlaneConfig(iz)%CellBdPt(1, i, isurf1, ixy1)
          Y = CY0 + PlaneConfig(iz)%CellBdPt(2, i, isurf1, ixy1)
          WRITE(io9, '(3F25.7)') X, Y, hz
          n = n0
        ENDIF
      ENDDO
    ELSE
      icel2 = Pin(ixy2)%Cell(iz)
      m = PlaneConfig(iz)%PtIdx(0, isurf2, ixy2)
      CX0 = PlaneConfig(iz)%CX(ix+1); CY0 = PlaneConfig(iz)%CY(iy)
      DO i = 1, m
        n0 = PlaneConfig(iz)%PtIdx(i, isurf2, ixy2)
        IF(n .LT. n0) THEN
          X = CX0 + PlaneConfig(iz)%CellBdPt(1, i, isurf2, ixy2)
          Y = CY0 + PlaneConfig(iz)%CellBdPt(2, i, isurf2, ixy2)
          WRITE(io9, '(3F25.7)') X, Y, hz
          n = n0
        ENDIF
      ENDDO
    ENDIF
  ENDDO
ENDDO

DO ixy = 1, nxy
  icel = Pin(ixy)%Cell(iz)
#ifndef CENT_OPT
  IF(Cell(icel)%lCentX) CYCLE
  IF(Cell(icel)%lCentY) CYCLE
  IF(Cell(icel)%lCentXY) CYCLE
#endif
  ix = Pin(ixy)%ix; iy = Pin(ixy)%iy
  CX0 = PlaneConfig(iz)%CX(ix); CY0 = PlaneConfig(iz)%CY(iy)
  DO j = 1, CellVTK(icel)%npt
    IF(CellVTK(icel)%lbdpt(j)) CYCLE
    x = CX0 + CellVTK(icel)%pt(1, j);y = CY0 + CellVTK(icel)%pt(2, j)
    WRITE(io9, '(3F25.7)') X, Y, hz
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE ProcessCell2D(Core, iz, ixy, iPtStIdx, ndat)
USE PARAM
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
INTEGER :: ixy, iz, iPtStIdx, ndat
INTEGER :: ElmtTypes(500)


TYPE(Pin_TYPE), POINTER :: Pin(:)

REAL :: BdPtList1(2, 500), BdPtList2(2, 500)
INTEGER :: BdIdxList1(0:500), BdIdxList2(0:500)
INTEGER :: BdIdxList1_chk(0:500), BdIdxList2_chk(0:500)
INTEGER :: SurfList(4), PtDir(4)
INTEGER :: PtIdxMap(500), BdPtLoc1(500), BdPtLoc2(500)
INTEGER :: BaseElmt(0:50), Elmt(0:50, 500)
INTEGER :: i, j, k, l, n, ist, iend, isurf, idir, icel
INTEGER :: idx0, idx, ist1, ist2, iend1, iend2

REAL :: x1, y1, x2, y2, z

LOGICAL :: lBdPt0, lBdPt, lElmt_in

SurfList = (/3, 4, 1, 2/)
PtDir = (/1, 1, - 1 , -1/)
Pin => Core%Pin
!Make BdList
icel = Core%Pin(ixy)%cell(iz)
DO i = 1, 4
  isurf = SurfList(i); idir = PtDir(i)
  IF(i==1) THEN
    BdPtList1(:, 1) = PlaneConfig(iz)%CellBdPt(:, 1, isurf, ixy)
    BdIdxList1(1) = PlaneConfig(iz)%PtIdx(1, isurf, ixy);
    BdIdxList1_chk(1) = -BdIdxList1(1)
    l = 1
  ENDIF
  IF(idir == 1) THEN
    DO j = 2, PlaneConfig(iz)%PTidx(0, isurf, ixy)
      l = l + 1
      BdPtList1(:, l) = PlaneConfig(iz)%CellBdPt(:, j, isurf, ixy)
      BdIdxList1(l) = PlaneConfig(iz)%PtIdx(j, isurf, ixy)
      BdIdxList1_chk(l) = BdIdxList1(l)
    ENDDO
    BdIdxList1_chk(l) = -BdIdxList1_chk(l)
  ELSE
    DO j =  PlaneConfig(iz)%PTidx(0, isurf, ixy) - 1, 1, -1
      l = l + 1
      BdPtList1(:, l) = PlaneConfig(iz)%CellBdPt(:, j, isurf, ixy)
      BdIdxList1(l) = PlaneConfig(iz)%PtIdx(j, isurf, ixy)
      BdIdxList1_chk(l) = BdIdxList1(l)
    ENDDO
    BdIdxList1_chk(l) = -BdIdxList1_chk(l)
  ENDIF
ENDDO
BdIdxList1(0) = l
DO i =1, l
  j =  l - i + 1
  BdIdxList2(i) = BdIdxList1(j)
  BdIdxList2_chk(i) = BdIdxList1_chk(j)
ENDDO

!Internal Point Index Change
DO i = 1, CellVTK(icel)%npt
  IF(CellVTK(icel)%lBdPt(i)) CYCLE
  iPtStIdx = iPtStIdx + 1
  PtIdxMap(i) = iPtStIdx
  !x1 = CellTVTK(icel)%Pt(1, i); y1 = CellVTK(icel)%Pt(2,i)
ENDDO
! BOundary Points Index Change
DO i = 1, CellVTK(icel)%npt
  IF(.NOT. CellVTK(icel)%lBdPt(i)) CYCLE
  !Seach Idx
  x1 = CellVTK(icel)%Pt(1, i); y1 = CellVTK(icel)%Pt(2,i)
  DO j = 1, BdIdxList1(0)
    x2 = BdPtList1(1, j); y2 = BdPtList1(2, j)
    z = (x1-x2)**2 + (y2-y1)**2
    IF(z .LT. epsm5) EXIT
  ENDDO
  PtIdxMap(i) = BdIdxList1(j)
  BdPtLoc1(i) = j
  BdPtLoc2(i) = BdIdxList1(0) - j + 1
ENDDO

!Element Transformation
DO i = 1, CellVTK(icel)%nelmt
  l = CellVTK(icel)%Elements(0,i)
  BaseElmt(0:l) =  CellVTK(icel)%Elements(0:l, i)
  BaseElmt(l+1) = BaseElmt(1)
  BaseElmt(1:l+1) = BaseElmt(1:l+1) +1
  Idx0 = BaseElmt(1)
  Elmt(1, i) = PtIdxMap(idx0)
  lBdpt0 = CellVTK(icel)%lbdpt(idx0)
  n = 1
  DO j = 2, l+1
    Idx = BaseElmt(j)
    lBdPt = CellVTK(icel)%lbdpt(idx)
    lElmt_In = .FALSE.


    IF(lBdPt .AND. lBdPt0) THEN  !Check Internal Elements or not
      ist1 = BdPtLoc1(idx0); iend1 = BdPtLoc1(idx)
      ist2 = BdPtLoc2(idx0); iend2 = BdPtLoc2(idx)
      CALL Chk_InElmt(lElmt_In, ist1, iend1, ist2, iend2, BdIdxList1_chk(0:500), BdIdxList2_chk(0:500))
!      DO k = ist1+1, iend1-1,1
!        IF(BdIdxList1_chk(k) .LT. 0) THEN
!          lElmt_In = .TRUE.
!          EXIT
!        ENDIF
!      ENDDO
!
!      DO k = ist2+1, iend2-1, 1
!        IF(BdIdxList2_chk(k) .LT. 0) THEN
!          lElmt_In = .TRUE.
!          EXIT
!        ENDIF
!      ENDDO
    ENDIF
    IF(lElmt_In) THEN
      IF(ist1 .EQ. 1) ist1 = BdIdxList1(0)
      IF(ist2 .EQ. BdIdxList1(0)) ist2 = 1
      IF(iend1 .EQ. 1) iend1 = BdIdxList1(0)
      IF(iend2 .EQ. BdIdxList1(0)) iend2 = 1
      lElmt_In = .FALSE.
      CALL Chk_InElmt(lElmt_In, ist1, iend1, ist2, iend2, BdIdxList1_chk(0:500), BdIdxList2_chk(0:500))
      IF(lElmt_In) THEN
        ist1 = BdPtLoc1(idx0); iend1 = BdPtLoc1(idx)
        ist2 = BdPtLoc2(idx0); iend2 = BdPtLoc2(idx)
      ENDIF
    ENDIF
    IF(.NOT. lBdPt .OR. .NOT. lBdPt0 .OR. lElmt_In) THEN
      n = n + 1
      Elmt(n, i) = PtIdxMap(idx)
      Idx0 = Idx
      lBdPt0 = lBdPt
      CYCLE
    ENDIF
    ! Boundary Line Processing
    IF(ist1 .LT. iend1) THEN
      DO k = ist1 + 1, iend1
        n = n + 1
        Elmt(n, i) = BdIdxList1(k)
      ENDDO
      Idx0 = idx
      lBdPt0 = lBdPt
    ELSE
      DO k = ist2 + 1, iend2
        n = n + 1
        Elmt(n, i) = BdIdxList2(k)
      ENDDO
      Idx0 = idx
      lBdPt0 = lBdPt
    ENDIF
  ENDDO
  Elmt(0, i) = n-1
  Elmt(1:n, i) = Elmt(1:n, i) - 1
  IF(Elmt(0, i) .EQ. 3) THEN
    ElmtTypes(i) = 5
  ELSEIF(Elmt(0, i) .EQ. 4) THEN
    ElmtTypes(i) = 9
  ELSE
    ElmtTypes(i) = 7
  ENDIF
  WRITE(io10, '(50I8)') Elmt(0:n-1, i)
  WRITE(io11, '(50I8)') ElmtTypes(i)
  ndat = ndat + n
  !TYPE DECISION
ENDDO
NULLIFY(Pin)
END SUBROUTINE

SUBROUTINE Chk_InElmt(lElmt_In, ist1, iend1, ist2, iend2, BdIdxList1_chk, BdIdxList2_chk)
USE PARAM
IMPLICIT NONE
LOGICAL :: lElmt_In
INTEGER :: ist1, iend1, ist2, iend2
INTEGER :: BdIdxList1_chk(0:500), BdIdxList2_chk(0:500)
INTEGER :: k
lElmt_In = .FALSE.
DO k = ist1+1, iend1-1,1
  IF(BdIdxList1_chk(k) .LT. 0) THEN
    lElmt_In = .TRUE.
    EXIT
  ENDIF
ENDDO

DO k = ist2+1, iend2-1, 1
  IF(BdIdxList2_chk(k) .LT. 0) THEN
    lElmt_In = .TRUE.
    EXIT
  ENDIF
ENDDO
END SUBROUTINE

SUBROUTINE WriteGeomVtk(Core, iz)
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
INTEGER :: iz

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
INTEGER, POINTER :: NODE(:, :)
INTEGER :: nxy, ndat, nx, ny, m, n

INTEGER :: ixy, icel, i, ix, iy
REAL :: CX0, CY0
REAL, POINTER :: CX(:), CY(:)

nxy = Core%nxy
nx = Core%nx; ny = Core%ny
Pin => CORE%Pin
Cell => Core%CellInfo

nGlobalpt(iz) = 0
nGlobalElmt(iz) = 0
ndat = 0
DO ixy = 1, nxy
  icel = Pin(ixy)%cell(iz)
  nGlobalpt(iz) = nGlobalpt(iz) + CellVTK(icel)%npt
  nGlobalElmt(iz) =nGlobalElmt(iz) + CellVTK(icel)%nelmt
  DO i = 1, CellVTK(icel)%nelmt
    ndat = ndat + (CellVTK(icel)%Elements(0, i) + 1)
  ENDDO
ENDDO

ALLOCATE(NODE(nx, ny))
ALLOCATE(CX(nx+1))
ALLOCATE(CY(ny+1))
CALL CP_CA(node, 0, nx, ny)
CALL CP_CA(CX, 0._8, nx+1)
CALL CP_CA(CY, 0._8, ny+1)
DO ixy = 1, Core%nxy
  ix = Pin(ixy)%ix; iy = Pin(ixy)%iy
  node(ix, iy) = ixy
ENDDO
DO iy = 1, ny
  IF(NODE(1, iy) .NE.0 .AND. NODE(nx, iy) .NE.0) EXIT
ENDDO
DO ix = 1, nx
  ixy = NODE(ix, iy)
  icel = Pin(ixy)%Cell(iz)
  IF(Cell(icel)%lCentXY) THEN
    CX(ix) = 0
    CX(ix+1) = Cell(icel)%geom%lx
  ELSE
    CX(ix)= CX(ix) + Cell(icel)%geom%lx/2._8
    CX(ix+1) = CX(ix) + Cell(icel)%geom%lx/2._8
  ENDIF
ENDDO

DO ix = 1, nx
  IF(NODE(ix, 1) .NE.0 .AND. NODE(ix, ny) .NE.0) EXIT
ENDDO
DO iy = 1, ny
  ixy = NODE(ix, iy)
  icel = Pin(ixy)%Cell(iz)
  IF(Cell(icel)%lCentXY) THEN
    CY(iy) = 0
    CY(iy+1) = -Cell(icel)%geom%ly
  ELSE
    CY(iy)= CY(iy) - Cell(icel)%geom%ly/2._8
    CY(iy+1) = CY(iy) - Cell(icel)%geom%ly/2._8
  ENDIF
ENDDO

WRITE(io9, '(A)') '# vtk DataFile Version 3.0'
WRITE(io9, '(A)') '2D Unstructured Grid of Linear Triangles'
WRITE(io9, '(A)') 'ASCII'
WRITE(io9, '(A)') 'DATASET UNSTRUCTURED_GRID'
WRITE(io9, '(A, x, I10, 2x, A)') 'POINTS', nGlobalPt(iz), 'float'
DO ixy = 1, nxy
  ix = Pin(ixy)%ix; iy = Pin(ixy)%iy
  CX0 = CX(ix); CY0 = CY(iy)
  icel = Pin(ixy)%Cell(iz)
  DO i = 1, CellVTK(icel)%npt
    WRITE(io9, '(3F25.7)') CellVtk(icel)%pt(1, i) + CX0, CellVtk(icel)%pt(2, i) + CY0, CellVtk(icel)%pt(3, i)
  ENDDO
ENDDO

WRITE(io9, '(A, x, 2I10)') 'CELLS', nGlobalElmt(iz), ndat

n = 0
DO ixy = 1, nxy
  ix = Pin(ixy)%ix; iy = Pin(ixy)%iy
  icel = Pin(ixy)%Cell(iz)
  DO i = 1, CellVTK(icel)%nelmt
    m = CellVTK(icel)%Elements(0, i)
    WRITE(io9, '(50I8)') CellVtk(icel)%elements(0, i), CellVtk(icel)%elements(1:m, i) + n
  ENDDO
  n = n + CellVTK(icel)%npt
ENDDO

WRITE(io9, '(A, x, 2I10)') 'CELL_TYPES', nGlobalelmt(iz)
DO ixy = 1, nxy
  ix = Pin(ixy)%ix; iy = Pin(ixy)%iy
  icel = Pin(ixy)%Cell(iz)
  DO i = 1, CellVTK(icel)%nelmt
    WRITE(io9, '(50I8)') CellVtk(icel)%types(i)
  ENDDO
ENDDO
WRITE(io9, *)
!  m = CellVtk(icel)%Elements(0, i)
!  WRITE(iunt, '(50I8)') CellVtk(icel)%elements(0:m, i)

!DO i = 1, CellVtk(icel)%nelmt
!  m = CellVtk(icel)%Elements(0, i)
!  WRITE(iunt, '(50I8)') CellVtk(icel)%elements(0:m, i)
!ENDDO
DEALLOCATE(CX, CY, node)
END SUBROUTINE

SUBROUTINE WriteTemp(Core, FmInfo, PowerDist, ThInfo, GroupInfo, iz, nTracerCntl)
USE PARAM
USE TYPEDEF,          ONLY : FxrInfo_Type,        XsMac_Type
USE XsUtil_mod,       ONLY : AllocXsMac
USE BasicOperation,   ONLY : CP_CA, CP_VA, MULTI_VA
USE BenchXs,          ONLY : xskfBen
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(PowerDist_Type) :: PowerDist
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
INTEGER :: iz

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :), myFxr
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(XsMac_Type), SAVE :: XsMac
REAL, POINTER :: Phis(:, :, :)

REAL, POINTER :: EQXS(:, :)
REAL, POINTER :: xsmackf(:)
REAL, POINTER :: pnum(:), SubGrpLv(:, :)
INTEGER, POINTER :: idiso(:)
REAL :: FsrTemp(500), pwsum, volsum, volsum0
REAL :: TempRef, Temp, TempSubGrp, XsMacNfold
INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nLocalFsr, nFsrInFxr
INTEGER :: ng, nFsr, nFxr, nxy, nz, myzb, myze, norg, nchi, niso
INTEGER :: iresoGrp1, iresoGrp2
INTEGER :: ixy, ixy0, ixya, icel, ifxr, ifsr, ifsrlocal, itype, ig
LOGICAL :: lXsLib, lRes
INTEGER :: i, j

lXsLib = nTracerCntl%lXsLib
iresoGrp1 = GroupInfo%nofg + 1; iresoGrp2 = GroupInfo%nofg + GroupInfo%norg
nxy = Core%nxy;ng = GroupInfo%ng
nFsr = Core%nCoreFsr;   nFxr = Core%nCoreFxr
nxy = Core%nxy
IF(lXsLib .and. .NOT. XsMac%lalloc) THEN
  XsMac%ng = ng
  CALL AllocXsMac(XsMac)
ELSE
  ALLOCATE(xsmackf(ng))
ENDIF
WRITE(io9, '(A, 2x, A, 2x, a, 2x, I1)') 'SCALARS', 'TEMP', 'double', 1
WRITE(io9, '(A)') 'LOOKUP_TABLE default'
Fxr => FmInfo%Fxr; Phis => FmInfo%Phis
Pin => Core%Pin; CellInfo => Core%CellInfo
DO ixy = 1, nxy
  FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
  icel = Pin(ixy)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
   ! IF(.NOT. Fxr(ifxr, iz)%lFuel) CYCLE
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    myFxr => Fxr(ifxr, iz)

    DO i = 1, nFsrInFxr
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(i,j)
      ifsr = FsrIdxSt + ifsrlocal - 1
      FsrTemp(ifsrlocal) = myFxr%temp - CKELVIN

    ENDDO
  ENDDO
  DO i = 1, CellInfo(icel)%nFsr
    WRITE(io9, '(e20.6)') FsrTemp(i)
  ENDDO
ENDDO
IF(.NOT. lxsLib) Deallocate(xsmackf)
END SUBROUTINE


SUBROUTINE Writeburnup(Core, FmInfo, PowerDist, ThInfo, GroupInfo, iz, nTracerCntl)
USE PARAM
USE TYPEDEF,          ONLY : FxrInfo_Type
USE BasicOperation,   ONLY : CP_CA, CP_VA, MULTI_VA
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(PowerDist_Type) :: PowerDist
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
INTEGER :: iz

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :), myFxr
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
REAL, POINTER :: Phis(:, :, :)

REAL, POINTER :: EQXS(:, :)
REAL, POINTER :: xsmackf(:)
REAL, POINTER :: pnum(:), SubGrpLv(:, :)
INTEGER, POINTER :: idiso(:)
REAL :: FsrBurnup(500), pwsum, volsum, volsum0
REAL :: TempRef, Temp, TempSubGrp, XsMacNfold
INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nLocalFsr, nFsrInFxr
INTEGER :: ng, nFsr, nFxr, nxy, nz, myzb, myze, norg, nchi, niso
INTEGER :: iresoGrp1, iresoGrp2
INTEGER :: ixy, ixy0, ixya, icel, ifxr, ifsr, ifsrlocal, itype, ig
LOGICAL :: lXsLib, lRes
INTEGER :: i, j

lXsLib = nTracerCntl%lXsLib
iresoGrp1 = GroupInfo%nofg + 1; iresoGrp2 = GroupInfo%nofg + GroupInfo%norg
nxy = Core%nxy;ng = GroupInfo%ng
nFsr = Core%nCoreFsr;   nFxr = Core%nCoreFxr
nxy = Core%nxy

WRITE(io9, '(A, 2x, A, 2x, a, 2x, I1)') 'SCALARS', 'BURNUP', 'double', 1
WRITE(io9, '(A)') 'LOOKUP_TABLE default'
Fxr => FmInfo%Fxr; Phis => FmInfo%Phis
Pin => Core%Pin; CellInfo => Core%CellInfo
DO ixy = 1, nxy
  FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
  icel = Pin(ixy)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
   ! IF(.NOT. Fxr(ifxr, iz)%lFuel) CYCLE
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    myFxr => Fxr(ifxr, iz)

    DO i = 1, nFsrInFxr
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(i,j)
      ifsr = FsrIdxSt + ifsrlocal - 1
      FsrBurnup(ifsrlocal) = myFxr%burnup
    ENDDO
  ENDDO
  DO i = 1, CellInfo(icel)%nFsr
    WRITE(io9, '(e20.6)') FsrBurnup(i)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE Writematerial(Core, FmInfo, iz, nTracerCntl)
USE PARAM
USE TYPEDEF,          ONLY : FxrInfo_Type
USE BasicOperation,   ONLY : CP_CA, CP_VA, MULTI_VA
USE TH_Mod,           ONLY : GetPinFuelTemp
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
INTEGER :: iz

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :), myFxr
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
REAL, POINTER :: Phis(:, :, :)

REAL, POINTER :: EQXS(:, :)
REAL, POINTER :: xsmackf(:)
REAL, POINTER :: pnum(:), SubGrpLv(:, :)
INTEGER, POINTER :: idiso(:)
REAL :: pwsum, volsum, volsum0
INTEGER :: FsrMaterial(500)
REAL :: TempRef, Temp, TempSubGrp, XsMacNfold
INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nLocalFsr, nFsrInFxr
INTEGER :: ng, nFsr, nFxr, nxy, nz, myzb, myze, norg, nchi, niso
INTEGER :: ixy, ixy0, ixya, icel, ifxr, ifsr, ifsrlocal, itype, ig
LOGICAL :: lXsLib, lRes
INTEGER :: i, j

lXsLib = nTracerCntl%lXsLib
nxy = Core%nxy
nFsr = Core%nCoreFsr;   nFxr = Core%nCoreFxr
nxy = Core%nxy

WRITE(io9, '(A, 2x, A, 2x, a, 2x, I1)') 'SCALARS', 'Material', 'int', 1
WRITE(io9, '(A)') 'LOOKUP_TABLE default'
Fxr => FmInfo%Fxr; Phis => FmInfo%Phis
Pin => Core%Pin; CellInfo => Core%CellInfo
DO ixy = 1, nxy
  FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
  icel = Pin(ixy)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    myFxr => Fxr(ifxr, iz)
    DO i = 1, nFsrInFxr
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(i,j)
      ifsr = FsrIdxSt + ifsrlocal - 1
      FsrMaterial(ifsrlocal) = myFxr%imix
    ENDDO
  ENDDO
  DO i = 1, CellInfo(icel)%nFsr
    WRITE(io9, '(I3)') FsrMaterial(i)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE WriteCoolTemp(Core, FmInfo, PowerDist, ThInfo, GroupInfo, iz, nTracerCntl)
USE PARAM
USE TYPEDEF,          ONLY : FxrInfo_Type,        XsMac_Type
USE XsUtil_mod,       ONLY : AllocXsMac
USE BasicOperation,   ONLY : CP_CA, CP_VA, MULTI_VA
USE BenchXs,          ONLY : xskfBen
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(PowerDist_Type) :: PowerDist
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
INTEGER :: iz

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :), myFxr
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(XsMac_Type), SAVE :: XsMac
REAL, POINTER :: Phis(:, :, :)

REAL, POINTER :: TCOOL(:, :)
REAL, POINTER :: xsmackf(:)
REAL, POINTER :: pnum(:), SubGrpLv(:, :)
INTEGER, POINTER :: idiso(:)
REAL :: FsrTemp(500), pwsum, volsum, volsum0
INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nLocalFsr, nFsrInFxr
INTEGER :: ng, nFsr, nFxr, nxy, nz, myzb, myze, norg, nchi, niso
INTEGER :: iresoGrp1, iresoGrp2
INTEGER :: ixy, ixy0, ixya, icel, ifxr, ifsr, ifsrlocal, itype, ig
LOGICAL :: lXsLib, lRes
INTEGER :: i, j

lXsLib = nTracerCntl%lXsLib
iresoGrp1 = GroupInfo%nofg + 1; iresoGrp2 = GroupInfo%nofg + GroupInfo%norg
nxy = Core%nxy;ng = GroupInfo%ng
nFsr = Core%nCoreFsr;   nFxr = Core%nCoreFxr
nxy = Core%nxy
Tcool => ThInfo%Tcool

WRITE(io9, '(A, 2x, A, 2x, a, 2x, I1)') 'SCALARS', 'COOLTEMP', 'double', 1
WRITE(io9, '(A)') 'LOOKUP_TABLE default'
Fxr => FmInfo%Fxr; Phis => FmInfo%Phis
Pin => Core%Pin; CellInfo => Core%CellInfo
DO ixy = 1, nxy
  FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
  icel = Pin(ixy)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
   ! IF(.NOT. Fxr(ifxr, iz)%lFuel) CYCLE
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    myFxr => Fxr(ifxr, iz)

    DO i = 1, nFsrInFxr
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(i,j)
      ifsr = FsrIdxSt + ifsrlocal - 1
      FsrTemp(ifsrlocal) = TCOOL(iz, ixy)
    ENDDO
  ENDDO
  DO i = 1, CellInfo(icel)%nFsr
    WRITE(io9, '(e20.6)') FsrTemp(i)
  ENDDO
ENDDO
IF(.NOT. lxsLib) Deallocate(xsmackf)
END SUBROUTINE

SUBROUTINE WritePower(Core, FmInfo, PowerDist, ThInfo, GroupInfo, iz, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,          ONLY : FxrInfo_Type,        XsMac_Type, PE_Type
USE MacXsLib_mod,     ONLY : EffMacXS,            EffRIFPin,        MacXsBase
USE XsUtil_mod,       ONLY : AllocXsMac
USE BasicOperation,   ONLY : CP_CA, CP_VA, MULTI_VA
USE BenchXs,          ONLY : xskfBen,             xskfDynBen
!USE TH_Mod,           ONLY : GetPinFuelTemp
USE TRAN_MOD,         ONLY : TranInfo
IMPLICIT NONE

TYPE(PE_Type) :: PE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(PowerDist_Type) :: PowerDist
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
INTEGER :: iz

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :), myFxr
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(XsMac_Type), SAVE :: XsMac
REAL, POINTER :: Phis(:, :, :)

REAL, POINTER :: EQXS(:, :)
REAL, POINTER :: xsmackf(:)
REAL, POINTER :: pnum(:), SubGrpLv(:, :)
INTEGER, POINTER :: idiso(:)
REAL :: PinPower(500), pwsum, volsum, volsum0
!REAL :: PinFuelTempAvgsq, Temp, XsMacFold
INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nLocalFsr, nFsrInFxr
INTEGER :: ng, nFsr, nFxr, nxy, nz, myzb, myze, norg, nchi, niso
INTEGER :: iresoGrp1, iresoGrp2
INTEGER :: ixy, ixy0, ixya, icel, ifxr, ifsr, ifsrlocal, itype, ig
LOGICAL :: lXsLib, lRes, FlagF
INTEGER :: i, j
LOGICAL :: lAIC

lXsLib = nTracerCntl%lXsLib
iresoGrp1 = GroupInfo%nofg + 1; iresoGrp2 = GroupInfo%nofg + GroupInfo%norg
nxy = Core%nxy;ng = GroupInfo%ng
nFsr = Core%nCoreFsr;   nFxr = Core%nCoreFxr
nxy = Core%nxy
IF(lXsLib .and. .NOT. XsMac%lalloc) THEN
  XsMac%ng = ng
  CALL AllocXsMac(XsMac)
ELSE
  ALLOCATE(xsmackf(ng))
ENDIF
WRITE(io9, '(A, 2x, A, 2x, a, 2x, I1)') 'SCALARS', 'Power', 'double', 1
WRITE(io9, '(A)') 'LOOKUP_TABLE default'
Fxr => FmInfo%Fxr; Phis => FmInfo%Phis
Pin => Core%Pin; CellInfo => Core%CellInfo
DO ixy = 1, nxy
  FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
  icel = Pin(ixy)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  PinPower =0
  !lAIC = CellInfo(icel)%lAIC
  !PinFuelTempAvgsq = dsqrt(GetPinFuelTemp(Core, FmInfo%Fxr, iz, ixy))
  !IF (ResVarPin(ixy,iz)%lres.and.nTracerCntl%lRIF) CALL EffRIFPin(ResVarPin(ixy,iz), PinFuelTempAvgsq, iz, lAIC, PE)
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    IF(.NOT. Fxr(ifxr, iz)%lFuel) CYCLE
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    myFxr => Fxr(ifxr, iz)
    IF(lXsLib) Then
      niso = myFxr%niso;
      CALL MacXsBase(XsMac, myFxr, 1, ng, ng, 1._8, FALSE, TRUE)
      XsMackf => XsMac%XsMackf
      IF(myFxr%lRes) THEN
        DO ig = iResoGrp1, iResoGrp2
          !XsMacFold = XsMac%XsMacF(ig)
          !CALL EffMacXs(XsMac, ResVarPin(ixy,iz), myFxr, PinFuelTempAvgsq, niso, ig, ng, .TRUE., iz, ixy, j, PE)
          !IF(XsMacFold .gt. epsm8) THEN
            XsMac%XsMacKF(ig) = XsMac%XsMacKF(ig) * myFxr%fresokf(ig)
          !ENDIF
        ENDDO
      ENDIF
    ELSE
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1,j)
      itype = myFxr%imix
      IF(nTracerCntl%lDynamicBen) THEN
        CALL xskfDynben(itype, TranInfo%fuelTemp(ixy, iz), 1, ng, xsmackf)
      ELSE
      CALL xskfben(itype, 1, ng, xsmackf)
    ENDIF
    ENDIF
    DO i = 1, nFsrInFxr
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(i,j)
      ifsr = FsrIdxSt + ifsrlocal - 1
      PinPower(ifsrlocal) = 0
      DO ig= 1, ng
        PinPower(ifsrlocal) = PinPower(ifsrlocal) + XsMackf(ig) * phis(ifsr, iz, ig)
      ENDDO
    ENDDO
  ENDDO
!  pwsum =0; volsum =0; volsum = 0
!  DO i = 1, CellInfo(icel)%nFsr
!    pwsum = pwsum + CellInfo(icel)%vol(i) * PinPower(i)
!    volsum0 = volsum0 + CellInfo(icel)%vol(i)
!    IF(PinPower(i) .GT. 0._8) THEN
!      volsum = volsum + CellInfo(icel)%vol(i)
!    ENDIF
!  ENDDO
!  pwsum = pwsum / volsum
!  ixya = Pin(ixy)%iasy
!  ixy0 = Pin(ixy)%ipin
!  pwsum = PowerDist%PinPower3D(ixy0, ixya, iz) / pwsum
  pwsum = PowerDist%Fm3DNormalizer
  DO i = 1, CellInfo(icel)%nFsr
    WRITE(io9, '(e20.6)') PinPower(i) * pwsum
  ENDDO
ENDDO
IF(.NOT. lxsLib) Deallocate(xsmackf)
END SUBROUTINE

SUBROUTINE WriteFlux(Core, FmInfo, normalizer, ng, iz)
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
INTEGER :: ng, iz
REAL :: normalizer
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)

INTEGER:: l, j, ig, ireg, FsrIdxSt, icel

Pin => Core%Pin; Cell => Core%CellInfo

DO ig = 1, ng
  IF(ig .LT. 10) THEN
    WRITE(io9, '(A, 2x, A, I1, 2x, a, 2x, I1)') 'SCALARS', 'Group00', ig, 'double', 1
    !SCALARS Group40 double 1
  ELSEiF(ig .LT. 100) THEN
    WRITE(io9, '(A, 2x, A, I2, 2x, a, 2x, I1)') 'SCALARS', 'Group0', ig, 'double', 1
  ELSE
    WRITE(io9, '(A, 2x, A, I3, 2x, a, 2x, I1)') 'SCALARS', 'Group', ig, 'double', 1
  ENDIF
  WRITE(io9, '(A)') 'LOOKUP_TABLE default'
  DO l = 1, Core%nxy
    FsrIdxSt = Pin(l)%FsrIdxSt; icel = Pin(l)%Cell(iz);
#ifndef CENT_OPT
    IF(Cell(icel)%lCentX) CYCLE
    IF(Cell(icel)%lCentY) CYCLE
    IF(Cell(icel)%lCentXY) CYCLE
#endif
    DO j = 1, Cell(icel)%nFsr
      ireg = FsrIdxSt + j - 1
      WRITE(io9, '(e20.6)') FmInfo%Phis(ireg, iz, ig)*normalizer
    ENDDO
  ENDDO
  WRITE(io9, '(e20.6)')
ENDDO
WRITE(io9, '(e20.6)')
ENDSUBROUTINE

SUBROUTINE WritePower_HOM(Core, CmInfo, normalizer, ng, myzb, myze)
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
INTEGER :: ng, myzb, myze
REAL :: normalizer
INTEGER :: nxy, nx, ny
INTEGER :: i, ig, iz, ixy, ix, iy
INTEGER, POINTER :: NODE(:, :)
REAL :: pwsum
nxy =  COre%nxy; nx = Core%nx; ny = Core%ny
ALLOCATE(NODE(nx, ny))
CALL CP_CA(NODE, 0, nx, ny)
DO ixy = 1, nxy
  ix = Core%Pin(ixy)%ix; iy = Core%Pin(ixy)%iy
  NODE(ix, iy) = ixy
ENDDO

WRITE(io9, '(A, 2x, A, 2x, a, 2x, I1)') 'SCALARS', 'Power', 'double', 1

WRITE(io9, '(A)') 'LOOKUP_TABLE default'
DO iz = myzb, myze
  DO ixy = 1, nxy
    pwsum = 0
    DO ig = 1, ng
      pwsum = pwsum + CmInfo%PhiC(ixy, iz, ig)*CmInfo%PinXs(ixy, iz)%xskf(ig)
    ENDDO
    pwsum = pwsum * normalizer
    WRITE(io9, '(e20.6)') pwsum
  ENDDO
ENDDO
DEALLOCATE(NODE)
END SUBROUTINE


SUBROUTINE WriteFlux_HOM(Core, CmInfo, normalizer, ng, myzb, myze)
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
INTEGER :: ng, myzb, myze
REAL :: normalizer, flux
INTEGER :: nxy, nx, ny
INTEGER :: i, ig, iz, ixy, ix, iy
INTEGER, POINTER :: NODE(:, :)

nxy =  COre%nxy; nx = Core%nx; ny = Core%ny
ALLOCATE(NODE(nx, ny))
CALL CP_CA(NODE, 0, nx, ny)
DO ixy = 1, nxy
  ix = Core%Pin(ixy)%ix; iy = Core%Pin(ixy)%iy
  NODE(ix, iy) = ixy
ENDDO

DO ig = 1, ng
  IF(ig .LT. 10) THEN
    WRITE(io9, '(A, 2x, A, I1, 2x, a, 2x, I1)') 'SCALARS', 'Group00', ig, 'double', 1
    !SCALARS Group40 double 1
  ELSEiF(ig .LT. 100) THEN
    WRITE(io9, '(A, 2x, A, I2, 2x, a, 2x, I1)') 'SCALARS', 'Group0', ig, 'double', 1
  ELSE
    WRITE(io9, '(A, 2x, A, I3, 2x, a, 2x, I1)') 'SCALARS', 'Group', ig, 'double', 1
  ENDIF
  WRITE(io9, '(A)') 'LOOKUP_TABLE default'
  DO iz = myzb, myze
    DO ixy = 1, nxy
      flux = CmInfo%PhiC(ixy, iz, ig)* normalizer
      WRITE(io9, '(e20.6)') flux
    ENDDO
  ENDDO
ENDDO
DEALLOCATE(NODE)
END SUBROUTINE


SUBROUTINE WriteTCool_HOM(Core, ThInfo, ng, myzb, myze)
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(ThInfo_Type) :: ThInfo
INTEGER :: ng, myzb, myze
INTEGER :: nxy, nx, ny
INTEGER :: i, ig, iz, ixy, ix, iy
INTEGER, POINTER :: NODE(:, :)
REAL :: pwsum
nxy =  COre%nxy; nx = Core%nx; ny = Core%ny

WRITE(io9, '(A, 2x, A, 2x, a, 2x, I1)') 'SCALARS', 'Temp_Cool', 'double', 1

WRITE(io9, '(A)') 'LOOKUP_TABLE default'
DO iz = myzb, myze
  DO ixy = 1, nxy
    WRITE(io9, '(e20.6)') ThInfo%TCool(iz, ixy)
    !pwsum = 0
    !DO ig = 1, ng
    !  pwsum = pwsum + CmInfo%PhiC(ixy, iz, ig)*CmInfo%PinXs(ixy, iz)%xskf(ig)
    !ENDDO
    !pwsum = pwsum * normalizer
    !WRITE(io9, '(e20.6)') pwsum
  ENDDO
ENDDO

END SUBROUTINE


SUBROUTINE WriteTFuelAvg_HOM(Core, ThInfo, myzb, myze)
USE BasicOperation, ONLY : CP_CA
USE TH_Mod,         ONLY : ThVar
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(ThInfo_Type) :: ThInfo
INTEGER :: myzb, myze
INTEGER :: nxy, nx, ny
INTEGER :: i, ig, iz, ixy, ix, iy
INTEGER, POINTER :: NODE(:, :)
INTEGER :: navg
REAL :: temp
nxy =  COre%nxy; nx = Core%nx; ny = Core%ny
navg = THVar%npr5
WRITE(io9, '(A, 2x, A, 2x, a, 2x, I1)') 'SCALARS', 'Temp_Fuel', 'double', 1
!FuelTh(ixy)%tfuel(navg, :)
WRITE(io9, '(A)') 'LOOKUP_TABLE default'
DO iz = myzb, myze
  DO ixy = 1, nxy
    temp = ThInfo%TCool(iz, ixy)
    IF(Core%Pin(ixy)%lfuel) THEN
      Temp = ThInfo%FuelTh(ixy)%tfuel(navg, iz)
    ENDIF
    WRITE(io9, '(e20.6)') Temp
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE BdCellVTK(Cell, icel)
uSE UtilFunction,   ONLY : Array2DSORT

USE ALLOCS
IMPLICIT NONE
TYPE(Cell_Type) :: Cell
INTEGER :: icel
REAL :: BdPts(2,100, 4)
INTEGER :: n, m, l
INTEGER :: i, j, nbdpt, nbdpt0(4)
REAL :: z, x(2), line(3)
REAL :: bdpt(2, 100, 4)

!Search Boundary Points
nbdpt = 0
ALLOCATE(CellVTK(icel)%lbdPT(CellVTK(icel)%npt))
CellVTK(icel)%lbdPT = .FALSE.
DO j = 1, 4
  line = Cell%Geom%bdline(:, j)
  DO i = 1, CellVTK(icel)%npt
    x = CellVtk(icel)%pt(1:2, i)
    z = x(1) * line(1) + x(2) * line(2) + line(3)
    IF(abs(z) .LT. 1.e-5) THEN   !Checking the points in the cell boundary
      IF(.NOT. CellVTK(icel)%lbdPT(i)) nbdpt = nbdpt + 1
      CellVTK(icel)%lbdPT(i) = .TRUE.

    ENDIF
  ENDDO
ENDDO
CellVTK(icel)%nbdpt = nbdpt
CellVTK(icel)%ninpt = CellVTK(icel)%npt - nbdpt
!Serach Boundary Points
CALL DMALLOC(CellVTK(icel)%BdDel,100, 4)
CALL Dmalloc(CellVTK(icel)%BdPtList, 2, 100, 4)
DO j = 1, 4
  nbdpt0(j) = 0
  line = Cell%Geom%bdline(:, j)
  DO i=1, CellVTK(icel)%npt
    !CellVTK(icel)%lbdPT(i) = .FALSE.
    x = CellVtk(icel)%pt(1:2, i)
    z = x(1) * line(1) + x(2) * line(2) + line(3)
    IF(abs(z) .LT. 1.e-5) THEN   !Checking the points in the cell boundary
      nbdpt0(j) = nbdpt0(j) + 1
      bdpt(:, nbdpt0(j), j) = x(1:2)
    ENDIF
  ENDDO
  IF(j == 1 .OR. j == 3) THEN
    l = 1
  ELSE
    l = 2
  ENDIF

  n = nbdpt0(j)
  IF(l .EQ. 1) THEN
    CALL Array2DSORT(bdpt(1, 1:n, j), bdpt(2, 1:n, j), n, m, .TRUE., .TRUE., l)
  ELSE
    CALL Array2DSORT(bdpt(1, 1:n, j), bdpt(2, 1:n, j), n, m, .FALSE., .TRUE., l)
  ENDIF
  CellVTK(icel)%nbddel(j) = m-1
  DO i = 1, m-1
    CellVTK(icel)%BdDel(i, j) = bdpt(l,i+1, j) - bdpt(l, i, j)
  ENDDO
  DO i = 1, m
    CellVTK(icel)%BdPtList(1:2, i, j) = bdpt(1:2, i, j)
  ENDDO
  CONTINUE

ENDDO

END SUBROUTINE

SUBROUTINE ConvertCellVTK3D(icel)
USE PARAM
USE Allocs
IMPLICIT NONE
INTEGER :: icel

INTEGER :: npt0, nelmt0, npt, nelmt
INTEGER :: i, j, k


npt0 = CellVTK(icel)%npt
nelmt0 = CellVTK(icel)%nelmt
npt = npt0 * 2
nelmt = nelmt0

CALL Dmalloc(CellVTK3D(icel)%pt, 3, npt)
j = 0
DO i = 1, npt0
  j = j + 1
  CellVTK3D(icel)%pt(:, j) = CellVTK(icel)%pt(:, i)
  CellVTK3D(icel)%pt(3, j) = 1
ENDDO

DO i = 1, npt0
  j = j+ 1
  CellVTK3D(icel)%pt(:, j) = CellVTK(icel)%pt(:, i)
  CellVTK3D(icel)%pt(3, j) = 2
ENDDO
CALL Dmalloc0(CellVTK3D(icel)%Elements, 0, 30, 1, nelmt)
CALL Dmalloc0(CellVTK3D(icel)%types, 1, nelmt)
DO i = 1, nelmt
  k = 0
  DO j = 1, CellVTK(icel)%elements(0, i)
    k = k + 1
    CellVTK3D(icel)%elements(k, i) = CellVTK(icel)%elements(j, i)
  ENDDO
  DO j = 1, CellVTK(icel)%elements(0, i)
    k = k + 1
    CellVTK3D(icel)%elements(k, i) = CellVTK(icel)%elements(j, i) + npt0
  ENDDO
  CellVTK3D(icel)%elements(0, i) = k
  IF(k==6) THEN
    CellVTK3D(icel)%types(i) = 13
  ELSEIF(k==8) THEN
    CellVTK3D(icel)%types(i) = 12
  ENDIF
ENDDO
CellVTK3D(icel)%npt = npt
CellVTK3D(icel)%nelmt = nelmt
END SUBROUTINE

SUBROUTINE MakeAnnualCellVTK(Cell, icel)
USE ALLOCS
IMPLICIT NONE
TYPE(Cell_Type) :: Cell
INTEGER :: icel

INTEGER :: i, j, k
INTEGER :: ncircle, nelmt
INTEGER :: Pt_Idx(12, 0:100)
REAL :: R(100), R_pt(2, 12, 0:100)
REAL :: lx, ly, x, y, vertex(2,4)
REAL :: l, cx, cy

LOGICAL :: lreg

nCircle = Cell%geom%ncircle
lx = Cell%Geom%lx; ly = Cell%Geom%ly
R(1:nCircle) = Cell%Geom%Circle(3, 1:nCircle)


lReg = .TRUE.


Vertex(1, 4) = Cell%geom%x(2); Vertex(2, 4) = Cell%geom%y(2)
Vertex(1, 3) = Cell%geom%x(2); Vertex(2, 3) = Cell%geom%y(1)
Vertex(1, 2) = Cell%geom%x(1); Vertex(2, 2) = Cell%geom%y(1)
Vertex(1, 1) = Cell%geom%x(1); Vertex(2, 1) = Cell%geom%y(2)

IF(Cell%lCentX) THEN
  Vertex(2, 1) = -Vertex(2,2); Vertex(2, 4) = -Vertex(2,3)
  lx = abs(Vertex(1,4)-Vertex(1,1));
  ly = abs(Vertex(1,4)-Vertex(1,1));
ENDIF
IF(Cell%lCentY) THEN
  Vertex(1, 1) = -Vertex(1,4); Vertex(1, 2) = -Vertex(1,3)
  lx = abs(Vertex(1,4)-Vertex(1,1));
  ly = abs(Vertex(1,4)-Vertex(1,1));
ENDIF
IF(Cell%lCentXY) THEN
  Vertex(:, 1) = -Vertex(:, 3); Vertex(2, 4) = -Vertex(2, 3); Vertex(1, 2) =  -Vertex(1,3)
  lx = abs(Vertex(1,4)-Vertex(1,1));
  ly = abs(Vertex(1,4)-Vertex(1,1));
ENDIF

IF(Cell%lCCell) THEN
  Vertex = 2 * Vertex
  lx = 2 * lx; ly = 2 * ly
ENDIF

IF(R(1) .GT. lx/2._8) lReg = .FALSE.

R_pt(1, 1, 0) = 0.5 * (Vertex(1, 1) + Vertex(1,4))
R_pt(2, 1, 0) = 0.5 * (Vertex(2, 1) + Vertex(2,4))
R_pt(1, 2, 0) = Vertex(1, 1); R_pt(2, 2, 0) = Vertex(2, 1)

R_pt(1, 3, 0) = 0.5 * (Vertex(1, 1) + Vertex(1,2))
R_pt(2, 3, 0) = 0.5 * (Vertex(2, 1) + Vertex(2,2))
R_pt(1, 4, 0) = Vertex(1, 2); R_pt(2, 4, 0) = Vertex(2, 2)

R_pt(1, 5, 0) = 0.5 * (Vertex(1, 2) + Vertex(1,3))
R_pt(2, 5, 0) = 0.5 * (Vertex(2, 2) + Vertex(2,3))
R_pt(1, 6, 0) = Vertex(1, 3); R_pt(2, 6, 0) = Vertex(2, 3)

R_pt(1, 7, 0) = 0.5 * (Vertex(1, 3) + Vertex(1,4))
R_pt(2, 7, 0) = 0.5 * (Vertex(2, 3) + Vertex(2,4))
R_pt(1, 8, 0) = Vertex(1, 4); R_pt(2, 8, 0) = Vertex(2, 4)

IF(lReg) THEN
  DO i = 1, nCircle
    DO j = 0, 7
      R_pt(1, j+1, i) = -R(i) * sin(pi/4*j)
      R_pt(2, j+1, i) = R(i) * cos(pi/4*j)
    ENDDO
  ENDDO

  k = 0
  DO i = 0, nCircle
    DO j = 1, 8
      k = k + 1
      Pt_Idx(j, i) = k-1
    ENDDO
  ENDDO

  CellVtk(icel)%npt =k+1
  CALL DMALLOC(CellVtk(icel)%pt, 3, k+1)
  CellVtk(icel)%pt(1:3, k+1) = (/Cell%Geom%cx,Cell%Geom%cy, 0._8/)
  k = 0
  DO i = 0, nCircle
    DO j = 1, 8
      k = k + 1
      CellVtk(icel)%pt(1:3, k) = (/R_pt(1, j, i), R_pt(2, j, i), 0._8/)
    ENDDO
  ENDDO

  nelmt = 8 * (nCircle+1)
  CellVtk(icel)%nelmt = nelmt
  CALL DMALLOC0(CellVtk(icel)%elements, 0, 4, 1, nelmt)
  CALL DMALLOC0(CellVtk(icel)%types, 1, nelmt)
  k = 0
  DO i = 0, nCircle-1
    DO j = 1, 7
      k = k + 1
      CellVtk(icel)%elements(1, k) = Pt_Idx(j, i)
      CellVtk(icel)%elements(2, k) = Pt_Idx(j+1, i)
      CellVtk(icel)%elements(3, k) = Pt_Idx(j+1, i+1)
      CellVtk(icel)%elements(4, k) = Pt_Idx(j, i+1)
      CellVtk(icel)%elements(0, k) = 4
      CellVtk(icel)%types(k) = 9
    ENDDO
    j=8
    k = k + 1
    CellVtk(icel)%elements(1, k) = Pt_Idx(j, i)
    CellVtk(icel)%elements(2, k) = Pt_Idx(1, i)
    CellVtk(icel)%elements(3, k) = Pt_Idx(1, i+1)
    CellVtk(icel)%elements(4, k) = Pt_Idx(j, i+1)
    CellVtk(icel)%elements(0, k) = 4
    CellVtk(icel)%types(k) = 9
  ENDDO
  !Center Region
  i = nCircle
  DO j = 1, 7
    k = k + 1
    CellVtk(icel)%elements(1, k) = Pt_Idx(j, i)
    CellVtk(icel)%elements(2, k) = Pt_Idx(j+1, i)
    CellVtk(icel)%elements(3, k) = CellVtk(icel)%npt-1
    CellVtk(icel)%elements(0, k) = 3
    CellVtk(icel)%types(k) = 5
  ENDDO
  j = 8
  k = k + 1
  CellVtk(icel)%elements(1, k) = Pt_Idx(j, i)
  CellVtk(icel)%elements(2, k) = Pt_Idx(1, i)
  CellVtk(icel)%elements(3, k) = CellVtk(icel)%npt-1
  CellVtk(icel)%elements(0, k) = 3
  CellVtk(icel)%types(k) = 5
ELSE
  DO i = 1, nCircle
    IF(R(i) .LT. lx/2._8) THEN
      DO j = 0, 7
        R_pt(1, j+1, i) = -R(i) * sin(pi/4*j)
        R_pt(2, j+1, i) = R(i) * cos(pi/4*j)
      ENDDO
    ELSE
      l = sqrt(R(i)*R(i) - lx * lx /4)
      R_pt(1, 1, i) = R_pt(1, 1, 0) - l ; R_pt(2, 1, i) = R_pt(2, 1, 0)
      R_pt(1, 2, i) =-R(i) / sqrt(2._8);  R_pt(2, 2, i) = R(i) / sqrt(2._8)
      R_pt(1, 3, i) = R_pt(1, 3, 0);      R_pt(2, 3, i) = R_pt(2, 3, 0) + l

      R_pt(1, 4, i) = R_pt(1, 3, 0);      R_pt(2, 4, i) = R_pt(2, 3, 0) - l
      R_pt(1, 5, i) =-R(i) / sqrt(2._8);  R_pt(2, 5, i) =-R(i) / sqrt(2._8)
      R_pt(1, 6, i) = R_pt(1, 5, 0) - l;  R_pt(2, 6, i) = R_pt(2, 5, 0)

      R_pt(1, 7, i) = R_pt(1, 5, 0) + l;  R_pt(2, 7, i) = R_pt(2, 5, 0)
      R_pt(1, 8, i) = R(i) / sqrt(2._8);  R_pt(2, 8, i) =-R(i) / sqrt(2._8)
      R_pt(1, 9, i) = R_pt(1, 7, 0);      R_pt(2, 9, i) = R_pt(2, 7, 0) - l

      R_pt(1, 10, i) = R_pt(1, 7, 0);     R_pt(2, 10, i) = R_pt(2, 7, 0) + l
      R_pt(1, 11, i) = R(i) / sqrt(2._8); R_pt(2, 11, i) = R(i) / sqrt(2._8)
      R_pt(1, 12, i) = R_pt(1, 1, 0) + l; R_pt(2, 12, i) = R_pt(2, 1, 0)
    ENDIF
  ENDDO

  k = 0
  DO j = 1, 12
    k = k + 1
    Pt_Idx(j, 0) = k-1
  ENDDO
  DO i = 1, nCircle
    IF(R(i) .LT. lx/2._8) THEN
      DO j = 1, 8
        k = k + 1
        Pt_Idx(j, i) = k-1
      ENDDO
    ELSE
      DO j = 1, 12
        k = k + 1
        Pt_Idx(j, i) = k-1
      ENDDO
    ENDIF
  ENDDO

  CellVtk(icel)%npt =k+1
  CALL DMALLOC(CellVtk(icel)%pt, 3, k+1)
  CellVtk(icel)%pt(1:3, k+1) = (/Cell%Geom%cx,Cell%Geom%cy, 0._8/)
  k = 0
  DO j = 1, 12
    k = k + 1
    CellVtk(icel)%pt(1:3, k) = (/R_pt(1, j, 0), R_pt(2, j, 0), 0._8/)
  ENDDO
  DO i = 1, nCircle
    IF(R(i) .LT. lx/2._8) THEN
      DO j = 1, 8
        k = k + 1
        CellVtk(icel)%pt(1:3, k) = (/R_pt(1, j, i), R_pt(2, j, i), 0._8/)
      ENDDO
    ELSE
      DO j = 1, 12
        k = k + 1
        CellVtk(icel)%pt(1:3, k) = (/R_pt(1, j, i), R_pt(2, j, i), 0._8/)
      ENDDO
    ENDIF
  ENDDO

  nelmt = 8 * (nCircle+1)
  CellVtk(icel)%nelmt = nelmt
  CALL DMALLOC0(CellVtk(icel)%elements, 0, 50, 1, nelmt)
  CALL DMALLOC0(CellVtk(icel)%types, 1, nelmt)
  k = 0
  !Outermost region
  DO j = 0, 3
    k = k + 1
    CellVtk(icel)%elements(1, k) = Pt_Idx(3*j+1, 1)
    CellVtk(icel)%elements(2, k) = Pt_Idx(3*j+2, 1)
    CellVtk(icel)%elements(3, k) = Pt_Idx(2*j+2, 0)
    CellVtk(icel)%elements(0, k) = 3
    CellVtk(icel)%types(k) = 5
    k = k + 1
    CellVtk(icel)%elements(1, k) = Pt_Idx(3*j+2, 1)
    CellVtk(icel)%elements(2, k) = Pt_Idx(3*j+3, 1)
    CellVtk(icel)%elements(3, k) = Pt_Idx(2*j+2, 0)
    CellVtk(icel)%elements(0, k) = 3
    CellVtk(icel)%types(k) = 5
  ENDDO

  DO i = 1, nCircle-1
    IF(R(i) .LT. lx/2._8) THEN
      DO j = 1, 7
        k = k + 1
        CellVtk(icel)%elements(1, k) = Pt_Idx(j, i)
        CellVtk(icel)%elements(2, k) = Pt_Idx(j+1, i)
        CellVtk(icel)%elements(3, k) = Pt_Idx(j+1, i+1)
        CellVtk(icel)%elements(4, k) = Pt_Idx(j, i+1)
        CellVtk(icel)%elements(0, k) = 4
        CellVtk(icel)%types(k) = 9
      ENDDO
      j=8
      k = k + 1
      CellVtk(icel)%elements(1, k) = Pt_Idx(j, i)
      CellVtk(icel)%elements(2, k) = Pt_Idx(1, i)
      CellVtk(icel)%elements(3, k) = Pt_Idx(1, i+1)
      CellVtk(icel)%elements(4, k) = Pt_Idx(j, i+1)
      CellVtk(icel)%elements(0, k) = 4
      CellVtk(icel)%types(k) = 9
    ELSEIF(R(i+1) .GT. lx/2._8) THEN
      DO j = 0, 3
        k = k + 1
        CellVtk(icel)%elements(1, k) = Pt_Idx(3*j + 1, i)
        CellVtk(icel)%elements(2, k) = Pt_Idx(3*j + 2, i)
        CellVtk(icel)%elements(3, k) = Pt_Idx(3*j + 2, i+1)
        CellVtk(icel)%elements(4, k) = Pt_Idx(3*j + 1, i+1)
        CellVtk(icel)%elements(0, k) = 4
        CellVtk(icel)%types(k) = 9
        k = k + 1
        CellVtk(icel)%elements(1, k) = Pt_Idx(3*j + 2, i)
        CellVtk(icel)%elements(2, k) = Pt_Idx(3*j + 3, i)
        CellVtk(icel)%elements(3, k) = Pt_Idx(3*j + 3, i+1)
        CellVtk(icel)%elements(4, k) = Pt_Idx(3*j + 2, i+1)
        CellVtk(icel)%elements(0, k) = 4
        CellVtk(icel)%types(k) = 9
      ENDDO
    ELSE
      DO j = 0, 2
        k = k + 1
        CellVtk(icel)%elements(1, k) = Pt_Idx(3*j + 1, i)
        CellVtk(icel)%elements(2, k) = Pt_Idx(3*j + 2, i)
        CellVtk(icel)%elements(3, k) = Pt_Idx(2*j + 2, i+1)
        CellVtk(icel)%elements(4, k) = Pt_Idx(2*j + 1, i+1)
        CellVtk(icel)%elements(5, k) = Pt_Idx(2*j+1, 0)
        CellVtk(icel)%elements(0, k) = 5
        CellVtk(icel)%types(k) = 7
        k = k + 1
        CellVtk(icel)%elements(1, k) = Pt_Idx(3*j + 3, i)
        CellVtk(icel)%elements(2, k) = Pt_Idx(3*j + 2, i)
        CellVtk(icel)%elements(3, k) = Pt_Idx(2*j + 2, i+1)
        CellVtk(icel)%elements(4, k) = Pt_Idx(2*j + 3, i+1)
        CellVtk(icel)%elements(5, k) = Pt_Idx(2*j+3, 0)
        CellVtk(icel)%elements(0, k) = 5
        CellVtk(icel)%types(k) = 7
      ENDDO
      k = k + 1; j =3
      CellVtk(icel)%elements(1, k) = Pt_Idx(3*j + 1, i)
      CellVtk(icel)%elements(2, k) = Pt_Idx(3*j + 2, i)
      CellVtk(icel)%elements(3, k) = Pt_Idx(2*j + 2, i+1)
      CellVtk(icel)%elements(4, k) = Pt_Idx(2*j + 1, i+1)
      CellVtk(icel)%elements(5, k) = Pt_Idx(7, 0)
      CellVtk(icel)%elements(0, k) = 5
      CellVtk(icel)%types(k) = 7
      k = k + 1
      CellVtk(icel)%elements(1, k) = Pt_Idx(3*j + 3, i)
      CellVtk(icel)%elements(2, k) = Pt_Idx(3*j + 2, i)
      CellVtk(icel)%elements(3, k) = Pt_Idx(2*j + 2, i+1)
      CellVtk(icel)%elements(4, k) = Pt_Idx(1, i+1)
      CellVtk(icel)%elements(5, k) = Pt_Idx(1, 0)
      CellVtk(icel)%elements(0, k) = 5
      CellVtk(icel)%types(k) = 7
    ENDIF
  ENDDO

  !Center Region
  i = nCircle
  DO j = 1, 7
    k = k + 1
    CellVtk(icel)%elements(1, k) = Pt_Idx(j, i)
    CellVtk(icel)%elements(2, k) = Pt_Idx(j+1, i)
    CellVtk(icel)%elements(3, k) = CellVtk(icel)%npt-1
    CellVtk(icel)%elements(0, k) = 3
    CellVtk(icel)%types(k) = 5
  ENDDO
  j = 8
  k = k + 1
  CellVtk(icel)%elements(1, k) = Pt_Idx(j, i)
  CellVtk(icel)%elements(2, k) = Pt_Idx(1, i)
  CellVtk(icel)%elements(3, k) = CellVtk(icel)%npt-1
  CellVtk(icel)%elements(0, k) = 3
  CellVtk(icel)%types(k) = 5
ENDIF

IF(Cell%lCCell) THEN
  CX = Cell%Geom%Circle(1, 1); CY = Cell%Geom%Circle(2, 1)
  DO k = 1, CellVtk(icel)%npt
    CellVtk(icel)%pt(1, k) = CellVtk(icel)%pt(1, k) + CX
    CellVtk(icel)%pt(2, k) = CellVtk(icel)%pt(2, k) + CY
  ENDDO
ENDIF

END SUBROUTINE

SUBROUTINE ProcessParialCell(Cell, icel)
USE ALLOCS
IMPLICIT NONE
TYPE(Cell_Type) :: Cell
INTEGER :: icel

REAL :: x1, x2, y1, y2
REAL :: Pt(3, 500)
INTEGER :: PtIdxMap(500), Elmt(0:50, 500), types(500)
INTEGER :: i, j, npt, nelmt, n
LOGICAL :: lElmt_In

x1 = min(Cell%Geom%x(1), Cell%Geom%x(2))
x2 = max(Cell%Geom%x(1), Cell%Geom%x(2))
y1 = min(Cell%Geom%y(1), Cell%Geom%y(2))
y2 = max(Cell%Geom%y(1), Cell%Geom%y(2))
n =0

DO i = 1,CellVTK(icel)%npt
  PtIdxMap(i) = 0
  IF((CellVTK(icel)%pt(1, i)-x1) * (CellVTK(icel)%pt(1, i)-x2) .GT. 1.e-5) CYCLE
  IF((CellVTK(icel)%pt(2, i)-y1) * (CellVTK(icel)%pt(2, i)-y2) .GT. 1.e-5) CYCLE
  n = n + 1
  PtIdxMap(i) = n
  Pt(:, n) = CellVTK(icel)%pt(:, i)
ENDDO
npt = n; n = 0
DO i = 1, CellVTK(icel)%nelmt
  lElmt_In = .FALSE.
  DO j = 1, CellVTK(icel)%Elements(0, i)
    IF(PtIdxMap(CellVTK(icel)%Elements(j,i)+1)==0) THEN
      lElmt_In = .TRUE.
      EXIT
    ENDIF
  ENDDO
  IF(lElmt_In) CYCLE
  n = n + 1
  DO j = 1, CellVTK(icel)%Elements(0, i)
    Elmt(j, n) = PtIdxMap(CellVTK(icel)%Elements(j,i)+1)-1
  ENDDO
  Elmt(0, n) = CellVTK(icel)%Elements(0, i)
  types(n) = CellVTK(icel)%types(i)
ENDDO
nElmt = n
DO i = 1, npt
  CellVTK(icel)%pt(:, i) = Pt(:, i)
ENDDO
CellVTK(icel)%npt = npt

DO i = 1, nElmt
  CellVTK(icel)%Elements(0, i) = Elmt(0, i)
  CellVTK(icel)%Elements(1:Elmt(0, i), i) = Elmt(1:Elmt(0, i), i)
  CellVTK(icel)%types(i) = types(i)
ENDDO
CellVTK(icel)%nElmt = nElmt
CONTINUE
END SUBROUTINE

SUBROUTINE MakeAnnualCellVTK2D(Cell, icel)
USE ALLOCS
IMPLICIT NONE
TYPE(Cell_Type) :: Cell
INTEGER :: icel

INTEGER :: i, j, k
INTEGER :: ncircle, nelmt
INTEGER :: Pt_Idx(12, 0:100)
REAL :: R(100), R_pt(2, 12, 0:100)
REAL :: lx, ly, x, y, vertex(2,4)
REAL :: l

LOGICAL :: lreg

nCircle = Cell%geom%ncircle
lx = Cell%Geom%lx; ly = Cell%Geom%ly
R(1:nCircle) = Cell%Geom%Circle(3, 1:nCircle)


lReg = .TRUE.


Vertex(1, 4) = Cell%geom%x(2); Vertex(2, 4) = Cell%geom%y(2)
Vertex(1, 3) = Cell%geom%x(2); Vertex(2, 3) = Cell%geom%y(1)
Vertex(1, 2) = Cell%geom%x(1); Vertex(2, 2) = Cell%geom%y(1)
Vertex(1, 1) = Cell%geom%x(1); Vertex(2, 1) = Cell%geom%y(2)

IF(Cell%lCentX) THEN
  Vertex(2, 1) = -Vertex(2,2); Vertex(2, 4) = -Vertex(2,3)
  lx = abs(Vertex(1,4)-Vertex(1,1));
  ly = abs(Vertex(1,4)-Vertex(1,1));
ENDIF
IF(Cell%lCentY) THEN
  Vertex(1, 1) = -Vertex(1,4); Vertex(1, 2) = -Vertex(1,3)
  lx = abs(Vertex(1,4)-Vertex(1,1));
  ly = abs(Vertex(1,4)-Vertex(1,1));
ENDIF
IF(Cell%lCentXY) THEN
  Vertex(:, 1) = -Vertex(:, 3); Vertex(2, 4) = -Vertex(2, 3); Vertex(1, 2) =  -Vertex(1,3)
  lx = abs(Vertex(1,4)-Vertex(1,1));
  ly = abs(Vertex(1,4)-Vertex(1,1));
ENDIF
IF(R(1) .GT. lx/2._8) lReg = .FALSE.

Vertex(1, 4) = Cell%geom%x(2); Vertex(2, 4) = Cell%geom%y(2)
Vertex(1, 3) = Cell%geom%x(2); Vertex(2, 3) = Cell%geom%y(1)
Vertex(1, 2) = Cell%geom%x(1); Vertex(2, 2) = Cell%geom%y(1)
Vertex(1, 1) = Cell%geom%x(1); Vertex(2, 1) = Cell%geom%y(2)

R_pt(1, 1, 0) = 0.5 * (Vertex(1, 1) + Vertex(1,4))
R_pt(2, 1, 0) = 0.5 * (Vertex(2, 1) + Vertex(2,4))
R_pt(1, 2, 0) = Vertex(1, 1); R_pt(2, 2, 0) = Vertex(2, 1)

R_pt(1, 3, 0) = 0.5 * (Vertex(1, 1) + Vertex(1,2))
R_pt(2, 3, 0) = 0.5 * (Vertex(2, 1) + Vertex(2,2))
R_pt(1, 4, 0) = Vertex(1, 2); R_pt(2, 4, 0) = Vertex(2, 2)

R_pt(1, 5, 0) = 0.5 * (Vertex(1, 2) + Vertex(1,3))
R_pt(2, 5, 0) = 0.5 * (Vertex(2, 2) + Vertex(2,3))
R_pt(1, 6, 0) = Vertex(1, 3); R_pt(2, 6, 0) = Vertex(2, 3)

R_pt(1, 7, 0) = 0.5 * (Vertex(1, 3) + Vertex(1,4))
R_pt(2, 7, 0) = 0.5 * (Vertex(2, 3) + Vertex(2,4))
R_pt(1, 8, 0) = Vertex(1, 4); R_pt(2, 8, 0) = Vertex(2, 4)

IF(lReg) THEN
  DO i = 1, nCircle
    DO j = 0, 7
      R_pt(1, j+1, i) = -R(i) * sin(pi/4*j)
      R_pt(2, j+1, i) = R(i) * cos(pi/4*j)
    ENDDO
  ENDDO

  k = 0
  DO i = 0, nCircle
    DO j = 1, 8
      k = k + 1
      Pt_Idx(j, i) = k-1
    ENDDO
  ENDDO

  CellVtk(icel)%npt =k+1
  CALL DMALLOC(CellVtk(icel)%pt, 3, k+1)
  CellVtk(icel)%pt(1:3, k+1) = (/Cell%Geom%cx,Cell%Geom%cy, 0._8/)
  k = 0
  DO i = 0, nCircle
    DO j = 1, 8
      k = k + 1
      CellVtk(icel)%pt(1:3, k) = (/R_pt(1, j, i), R_pt(2, j, i), 0._8/)
    ENDDO
  ENDDO

  nelmt = 8 * (nCircle+1)
  CellVtk(icel)%nelmt = nelmt
  CALL DMALLOC0(CellVtk(icel)%elements, 0, 4, 1, nelmt)
  CALL DMALLOC0(CellVtk(icel)%types, 1, nelmt)
  k = 0
  DO i = 0, nCircle-1
    DO j = 1, 7
      k = k + 1
      CellVtk(icel)%elements(1, k) = Pt_Idx(j, i)
      CellVtk(icel)%elements(2, k) = Pt_Idx(j+1, i)
      CellVtk(icel)%elements(3, k) = Pt_Idx(j+1, i+1)
      CellVtk(icel)%elements(4, k) = Pt_Idx(j, i+1)
      CellVtk(icel)%elements(0, k) = 4
      CellVtk(icel)%types(k) = 9
    ENDDO
    j=8
    k = k + 1
    CellVtk(icel)%elements(1, k) = Pt_Idx(j, i)
    CellVtk(icel)%elements(2, k) = Pt_Idx(1, i)
    CellVtk(icel)%elements(3, k) = Pt_Idx(1, i+1)
    CellVtk(icel)%elements(4, k) = Pt_Idx(j, i+1)
    CellVtk(icel)%elements(0, k) = 4
    CellVtk(icel)%types(k) = 9
  ENDDO
  !Center Region
  i = nCircle
  DO j = 1, 7
    k = k + 1
    CellVtk(icel)%elements(1, k) = Pt_Idx(j, i)
    CellVtk(icel)%elements(2, k) = Pt_Idx(j+1, i)
    CellVtk(icel)%elements(3, k) = CellVtk(icel)%npt-1
    CellVtk(icel)%elements(0, k) = 3
    CellVtk(icel)%types(k) = 5
  ENDDO
  j = 8
  k = k + 1
  CellVtk(icel)%elements(1, k) = Pt_Idx(j, i)
  CellVtk(icel)%elements(2, k) = Pt_Idx(1, i)
  CellVtk(icel)%elements(3, k) = CellVtk(icel)%npt-1
  CellVtk(icel)%elements(0, k) = 3
  CellVtk(icel)%types(k) = 5
ELSE
  DO i = 1, nCircle
    IF(R(i) .LT. lx/2._8) THEN
      R(i) = lx/2._8-1.e-5
      EXiT
    ENDIF
  ENDDO
  DO i = 1, nCircle
    IF(R(i) .LT. lx/2._8) THEN
      DO j = 0, 7
        R_pt(1, j+1, i) = -R(i) * sin(pi/4*j)
        R_pt(2, j+1, i) = R(i) * cos(pi/4*j)
      ENDDO
    ELSE
      l = sqrt(R(i)*R(i) - lx * lx /4)
      R_pt(1, 1, i) = R_pt(1, 1, 0) - l ; R_pt(2, 1, i) = R_pt(2, 1, 0)
      R_pt(1, 2, i) =-R(i) / sqrt(2._8);  R_pt(2, 2, i) = R(i) / sqrt(2._8)
      R_pt(1, 3, i) = R_pt(1, 3, 0);      R_pt(2, 3, i) = R_pt(2, 3, 0) + l

      R_pt(1, 4, i) = R_pt(1, 3, 0);      R_pt(2, 4, i) = R_pt(2, 3, 0) - l
      R_pt(1, 5, i) =-R(i) / sqrt(2._8);  R_pt(2, 5, i) =-R(i) / sqrt(2._8)
      R_pt(1, 6, i) = R_pt(1, 5, 0) - l;  R_pt(2, 6, i) = R_pt(2, 5, 0)

      R_pt(1, 7, i) = R_pt(1, 5, 0) + l;  R_pt(2, 7, i) = R_pt(2, 5, 0)
      R_pt(1, 8, i) = R(i) / sqrt(2._8);  R_pt(2, 8, i) =-R(i) / sqrt(2._8)
      R_pt(1, 9, i) = R_pt(1, 7, 0);      R_pt(2, 9, i) = R_pt(2, 7, 0) - l

      R_pt(1, 10, i) = R_pt(1, 7, 0);     R_pt(2, 10, i) = R_pt(2, 7, 0) + l
      R_pt(1, 11, i) = R(i) / sqrt(2._8); R_pt(2, 11, i) = R(i) / sqrt(2._8)
      R_pt(1, 12, i) = R_pt(1, 1, 0) + l; R_pt(2, 12, i) = R_pt(2, 1, 0)
    ENDIF
  ENDDO

  k = 0
  DO j = 1, 12
    k = k + 1
    Pt_Idx(j, 0) = k-1
  ENDDO
  DO i = 1, nCircle
    IF(R(i) .LT. lx/2._8) THEN
      DO j = 1, 8
        k = k + 1
        Pt_Idx(j, i) = k-1
      ENDDO
    ELSE
      DO j = 1, 12
        k = k + 1
        Pt_Idx(j, i) = k-1
      ENDDO
    ENDIF
  ENDDO

  CellVtk(icel)%npt =k+1
  CALL DMALLOC(CellVtk(icel)%pt, 3, k+1)
  CellVtk(icel)%pt(1:3, k+1) = (/Cell%Geom%cx,Cell%Geom%cy, 0._8/)
  k = 0
  DO j = 1, 12
    k = k + 1
    CellVtk(icel)%pt(1:3, k) = (/R_pt(1, j, 0), R_pt(2, j, 0), 0._8/)
  ENDDO
  DO i = 1, nCircle
    IF(R(i) .LT. lx/2._8) THEN
      DO j = 1, 8
        k = k + 1
        CellVtk(icel)%pt(1:3, k) = (/R_pt(1, j, i), R_pt(2, j, i), 0._8/)
      ENDDO
    ELSE
      DO j = 1, 12
        k = k + 1
        CellVtk(icel)%pt(1:3, k) = (/R_pt(1, j, i), R_pt(2, j, i), 0._8/)
      ENDDO
    ENDIF
  ENDDO

  nelmt = 8 * (nCircle+1)
  CellVtk(icel)%nelmt = nelmt
  CALL DMALLOC0(CellVtk(icel)%elements, 0, 4, 1, nelmt)
  CALL DMALLOC0(CellVtk(icel)%types, 1, nelmt)
  k = 0
  !Outermost region
  DO j = 0, 3
    k = k + 1
    CellVtk(icel)%elements(1, k) = Pt_Idx(3*j+1, 1)
    CellVtk(icel)%elements(2, k) = Pt_Idx(3*j+2, 1)
    CellVtk(icel)%elements(3, k) = Pt_Idx(2*j+2, 0)
    CellVtk(icel)%elements(0, k) = 3
    CellVtk(icel)%types(k) = 5
    k = k + 1
    CellVtk(icel)%elements(1, k) = Pt_Idx(3*j+2, 1)
    CellVtk(icel)%elements(2, k) = Pt_Idx(3*j+3, 1)
    CellVtk(icel)%elements(3, k) = Pt_Idx(2*j+2, 0)
    CellVtk(icel)%elements(0, k) = 3
    CellVtk(icel)%types(k) = 5
  ENDDO

  DO i = 1, nCircle-1
    IF(R(i) .LT. lx/2._8) THEN
      DO j = 1, 7
        k = k + 1
        CellVtk(icel)%elements(1, k) = Pt_Idx(j, i)
        CellVtk(icel)%elements(2, k) = Pt_Idx(j+1, i)
        CellVtk(icel)%elements(3, k) = Pt_Idx(j+1, i+1)
        CellVtk(icel)%elements(4, k) = Pt_Idx(j, i+1)
        CellVtk(icel)%elements(0, k) = 4
        CellVtk(icel)%types(k) = 9
      ENDDO
      j=8
      k = k + 1
      CellVtk(icel)%elements(1, k) = Pt_Idx(j, i)
      CellVtk(icel)%elements(2, k) = Pt_Idx(1, i)
      CellVtk(icel)%elements(3, k) = Pt_Idx(1, i+1)
      CellVtk(icel)%elements(4, k) = Pt_Idx(j, i+1)
      CellVtk(icel)%elements(0, k) = 4
      CellVtk(icel)%types(k) = 9
    ELSEIF(R(i+1) .GT. lx/2._8) THEN
      DO j = 0, 3
        k = k + 1
        CellVtk(icel)%elements(1, k) = Pt_Idx(3*j + 1, i)
        CellVtk(icel)%elements(2, k) = Pt_Idx(3*j + 2, i)
        CellVtk(icel)%elements(3, k) = Pt_Idx(3*j + 2, i+1)
        CellVtk(icel)%elements(4, k) = Pt_Idx(3*j + 1, i+1)
        CellVtk(icel)%elements(0, k) = 4
        CellVtk(icel)%types(k) = 9
        k = k + 1
        CellVtk(icel)%elements(1, k) = Pt_Idx(3*j + 2, i)
        CellVtk(icel)%elements(2, k) = Pt_Idx(3*j + 3, i)
        CellVtk(icel)%elements(3, k) = Pt_Idx(3*j + 3, i+1)
        CellVtk(icel)%elements(4, k) = Pt_Idx(3*j + 2, i+1)
        CellVtk(icel)%elements(0, k) = 4
        CellVtk(icel)%types(k) = 9
      ENDDO
    ELSE
      DO j = 0, 2
        k = k + 1
        CellVtk(icel)%elements(1, k) = Pt_Idx(3*j + 1, i)
        CellVtk(icel)%elements(2, k) = Pt_Idx(3*j + 2, i)
        CellVtk(icel)%elements(3, k) = Pt_Idx(2*j + 2, i+1)
        CellVtk(icel)%elements(4, k) = Pt_Idx(2*j + 1, i+1)
        CellVtk(icel)%elements(0, k) = 4
        CellVtk(icel)%types(k) = 9
        k = k + 1
        CellVtk(icel)%elements(1, k) = Pt_Idx(3*j + 2, i)
        CellVtk(icel)%elements(2, k) = Pt_Idx(3*j + 3, i)
        CellVtk(icel)%elements(3, k) = Pt_Idx(2*j + 3, i+1)
        CellVtk(icel)%elements(4, k) = Pt_Idx(2*j + 2, i+1)
        CellVtk(icel)%elements(0, k) = 4
        CellVtk(icel)%types(k) = 9
      ENDDO
      k = k + 1; j =3
      CellVtk(icel)%elements(1, k) = Pt_Idx(3*j + 1, i)
      CellVtk(icel)%elements(2, k) = Pt_Idx(3*j + 2, i)
      CellVtk(icel)%elements(3, k) = Pt_Idx(2*j + 2, i+1)
      CellVtk(icel)%elements(4, k) = Pt_Idx(2*j + 1, i+1)
      CellVtk(icel)%elements(0, k) = 4
      CellVtk(icel)%types(k) = 9
      k = k + 1
      CellVtk(icel)%elements(1, k) = Pt_Idx(3*j + 2, i)
      CellVtk(icel)%elements(2, k) = Pt_Idx(3*j + 3, i)
      CellVtk(icel)%elements(3, k) = Pt_Idx(1, i+1)
      CellVtk(icel)%elements(4, k) = Pt_Idx(2*j + 2, i+1)
      CellVtk(icel)%elements(0, k) = 4
      CellVtk(icel)%types(k) = 9
    ENDIF
  ENDDO

  !Center Region
  i = nCircle
  DO j = 1, 7
    k = k + 1
    CellVtk(icel)%elements(1, k) = Pt_Idx(j, i)
    CellVtk(icel)%elements(2, k) = Pt_Idx(j+1, i)
    CellVtk(icel)%elements(3, k) = CellVtk(icel)%npt-1
    CellVtk(icel)%elements(0, k) = 3
    CellVtk(icel)%types(k) = 5
  ENDDO
  j = 8
  k = k + 1
  CellVtk(icel)%elements(1, k) = Pt_Idx(j, i)
  CellVtk(icel)%elements(2, k) = Pt_Idx(1, i)
  CellVtk(icel)%elements(3, k) = CellVtk(icel)%npt-1
  CellVtk(icel)%elements(0, k) = 3
  CellVtk(icel)%types(k) = 5

ENDIF
END SUBROUTINE
!
!type basicgeom
!  INTEGER :: nbd
!  logical :: lcircle,lrect,lCCent
!  INTEGER :: ncircle,nline
!  INTEGER :: nx,ny
!  REAL,pointer :: bdline(:,:)  !Boundary Line
!  REAL,pointer :: bdlineRange(:,:)
!  REAL,pointer :: line(:,:)
!  REAL,pointer :: circle(:,:)
!  REAL :: cx,cy, lx,ly, x(2), y(2)
!  REAL,POINTER :: delx(:), dely(:)
!END type
SUBROUTINE MakeRectCellVTK(Cell, icel)
USE ALLOCS
IMPLICIT NONE
TYPE(Cell_Type) :: Cell
INTEGER :: icel

INTEGER :: i, j, k
INTEGER :: nx, ny, npt, nelmt
INTEGER :: Pt_Idx(0:100, 0:100)
REAL :: x_pt(0:100), y_pt(0:100)
REAL :: lx, ly, cx, cy

nx = Cell%Geom%nx; ny = Cell%Geom%ny
cx = Cell%Geom%cx; cy = Cell%Geom%cy

x_pt(0) = Cell%Geom%x(1)
y_pt(0) = Cell%Geom%y(2)

DO i = 1, nx
  x_pt(i) = x_pt(i-1) + Cell%Geom%delx(i)
ENDDO

DO i = 1, ny
  y_pt(i) = y_pt(i-1) - Cell%Geom%dely(i)
ENDDO

npt = (nx + 1) * (ny + 1)
CALL Dmalloc(CellVtk(icel)%pt, 3, npt)
CellVtk(icel)%npt = npt
k = 0
DO i = 0, ny
  DO j = 0, nx
    k = k + 1
    Pt_Idx(j, i) = k -1
    CellVTK(icel)%pt(1:3, k) = (/x_pt(j), y_pt(i), 0._8/)
  ENDDO
ENDDO

nelmt = nx * ny
CALL Dmalloc0(CellVtk(icel)%elements, 0, 4, 1, nelmt)
CALL Dmalloc0(CellVtk(icel)%types, 1, nelmt)
CellVtk(icel)%nelmt = nelmt

k =0
DO i = 0, ny-1
  DO j = 0, nx-1
    k = k + 1
    CellVTK(icel)%elements(1, k) = Pt_Idx(j, i)
    CellVTK(icel)%elements(2, k) = Pt_Idx(j+1, i)
    CellVTK(icel)%elements(3, k) = Pt_Idx(j+1, i+1)
    CellVTK(icel)%elements(4, k) = Pt_Idx(j, i+1)
    CellVTK(icel)%elements(0, k) = 4
    CellVTK(icel)%types(k) = 9
  ENDDO
ENDDO

CONTINUE

END SUBROUTINE

SUBROUTINE WriteCellVtk(icel)
IMPLICIT NONE
INTEGER :: icel
INTEGER :: i, n,m
INTEGER :: iunt

iunt = 44
OPEN(unit=iunt, file ='test0.vtk', status = 'REPLACE')
WRITE(iunt, '(A)') '# vtk DataFile Version 3.0'
WRITE(iunt, '(A)') '2D Unstructured Grid of Linear Triangles'
WRITE(iunt, '(A)') 'ASCII'
WRITE(iunt, '(A)') 'DATASET UNSTRUCTURED_GRID'
n = CellVtk(icel)%npt
WRITE(iunt, '(A, x, I10, 2x, A)') 'POINTS', CellVtk(icel)%npt, 'float'
DO i = 1, CellVtk(icel)%npt
  WRITE(iunt, '(3F25.7)') CellVtk(icel)%pt(1:3, i)
ENDDO
n = 0
DO i = 1, CellVtk(icel)%nelmt
  n = n + CellVtk(icel)%elements(0, i)+1
ENDDO

WRITE(iunt, '(A, x, 2I10)') 'CELLS', CellVtk(icel)%nelmt, n
DO i = 1, CellVtk(icel)%nelmt
  m = CellVtk(icel)%Elements(0, i)
  WRITE(iunt, '(50I8)') CellVtk(icel)%elements(0:m, i)
ENDDO
WRITE(iunt, '(A, x, 2I10)') 'CELL_TYPES', CellVtk(icel)%nelmt
DO i = 1, CellVtk(icel)%nelmt
  WRITE(iunt, '(I8)') CellVtk(icel)%types(i)
ENDDO
WRITE(iunt, *)
CLOSE(iunt)
END SUBROUTINE

SUBROUTINE WriteCellVtk3D(icel)
IMPLICIT NONE
INTEGER :: icel
INTEGER :: i, n,m
INTEGER :: iunt

iunt = 44
OPEN(unit=iunt, file ='test1.vtk', status = 'REPLACE')
WRITE(iunt, '(A)') '# vtk DataFile Version 3.0'
WRITE(iunt, '(A)') '2D Unstructured Grid of Linear Triangles'
WRITE(iunt, '(A)') 'ASCII'
WRITE(iunt, '(A)') 'DATASET UNSTRUCTURED_GRID'
n = CellVtk3D(icel)%npt
WRITE(iunt, '(A, x, I10, 2x, A)') 'POINTS', CellVtk3D(icel)%npt, 'float'
DO i = 1, CellVtk3D(icel)%npt
  WRITE(iunt, '(3F25.7)') CellVtk3D(icel)%pt(1:3, i)
ENDDO
n = 0
DO i = 1, CellVtk3D(icel)%nelmt
  n = n + CellVtk3D(icel)%elements(0, i)+1
ENDDO

WRITE(iunt, '(A, x, 2I10)') 'CELLS', CellVtk3D(icel)%nelmt, n
DO i = 1, CellVtk3D(icel)%nelmt
  m = CellVtk3D(icel)%Elements(0, i)
  WRITE(iunt, '(50I8)') CellVtk3D(icel)%elements(0:m, i)
ENDDO
WRITE(iunt, '(A, x, 2I10)') 'CELL_TYPES', CellVtk3D(icel)%nelmt
DO i = 1, CellVtk3D(icel)%nelmt
  WRITE(iunt, '(I8)') CellVtk3D(icel)%types(i)
ENDDO
WRITE(iunt, *)
CLOSE(iunt)
END SUBROUTINE

FUNCTION FluxNorm(Core, FmInfo, CmInfo, ng, nTracerCntl, PE)
! Flux level normalization
! normalize flux such that average flux in fuel region be unity
! then update fission source and moments accordingly
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,     FmInfo_Type,       CmInfo_Type,     &
                           PE_TYPE,                                               &
                           PinXs_Type
USE CNTL,           ONLY : nTracerCntl_Type
USE BasicOperation, ONLY : MULTI_CA
#ifdef MPI_ENV
USE MPIComm_Mod,    ONLY : REDUCE
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
INTEGER :: ng
REAL :: FluxNorm


TYPE(PinXs_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: PhiC(:, :, :), PsiC(:, :), PhiFm(:, :, :), PsiFm(:, :)
REAL, POINTER :: Phis(:, :, :), Psi(:, :), Jout(:, :, :, :, :)
REAL, POINTER :: PinVol(:, :)

INTEGER :: nbd, nfsr, nxy, myzb, myze, myzbf, myzef
INTEGER :: ixy, iz, ig

REAL :: vol, volsum, pwsum, phisum, pinpw
REAL :: avgpw, avgflx, norm
REAL :: UnitPowerLevel0, Plevel

INTEGER :: comm
REAL :: buf0(3), buf(3)

COMM = PE%MPI_CMFD_COMM

nbd = 4
nxy = Core%nxy; nFsr = Core%nCoreFsr
myzb = PE%myzb; myze = PE%myze
myzbf = PE%myzbf; myzef =PE%myzef

PhiC => CmInfo%PhiC; PsiC => CmInfo%PsiC
PhiFm => CmInfo%PhiFm; PsiFm => CmInfo%PsiFm

Phis => FmInfo%Phis; Psi => FmInfo%Psi
Jout => CMInfo%RadJout
PinVol => Core%PinVol; PinXs => CmInfo%PinXs

volsum = 0;
pwsum = 0;
phisum = 0;
DO iz = myzb, myze
  DO ixy = 1, nxy
    pinpw = 0
    DO ig = 1, ng
      PinPw = PinPw + PinXs(ixy, iz)%xskf(ig) * PhiC(ixy, iz, ig)
    ENDDO
    IF(PinPw .LT. 0._8) CYCLE
    vol = PinVol(ixy, iz)
    phisum = phisum + sum(PhiC(ixy, iz, 1:ng)) * vol
    PwSum = PwSum + PinPw * vol
    volsum = volsum + vol
  ENDDO
ENDDO

#ifdef MPI_ENV
buf0 = (/PwSum, PhiSum, volsum/)
CALL REDUCE(buf0, buf, 3, COMM, .TRUE.)
PwSum = Buf(1); PhiSum = Buf(2); Volsum = Buf(3)
#endif

AvgPw = PwSum / VolSum
avgflx = PhiSum / Volsum

norm = ng / avgflx
FluxNorm = norm
IF(nTracerCntl%lProblem .EQ. lTransient) THEN
  FluxNorm  = FluxNorm* nTracerCntl%PowerLevel
ENDIF

END FUNCTION

END MODULE
