#include <defines.h>
MODULE SubChGen_mod
USE PARAM
USE TYPEDEF,                ONLY : CoreInfo_Type,        PE_Type,           &
                                   AsyInfo_Type,         Cell_Type,          Pin_Type
USE CNTL,                   ONLY : nTracerCntl_Type
USE BasicOperation,         ONLY : CP_CA,                CP_VA
USE IOUTIL,                 ONLY : Terminate
IMPLICIT NONE

TYPE ChInfo_Type
  !Each Subchannel Information
  LOGICAL :: lUse = .FALSE.
  REAL :: area = 0 
  REAL :: w_peri = 0      !Wetted Perimeter
  REAL :: h_peri = 0      !Heated Perimeter
  INTEGER :: ngap = 0     !Number of Gap (Gap: interface between gap)
  INTEGER, POINTER :: GapList(:)   
  INTEGER, POINTER :: NeighCh(:)
  
  INTEGER :: nrod = 0
  INTEGER, POINTER :: RodList(:)
  INTEGER, POINTER :: RodConf(:, :)
  
  INTEGER :: ncell = 0
  INTEGER, POINTER :: CellList(:)
  
  INTEGER :: nNode = 0
  INTEGER, POINTER :: NodeList(:, :)

  LOGICAL :: lGT = .FALSE.    !Large Guide Tube
  LOGICAL :: lLGT = .FALSE.   !
  INTEGER :: ngt = 0
END TYPE

TYPE RodInfo_Type
  LOGICAL :: lUse = .False.
  LOGICAL :: lFuel = .FALSE.
  LOGICAL :: lGT = .FALSE.
  LOGICAL :: lLGT = .FALSE.
  INTEGER :: ixnc, iync
  
  REAL :: area = 0
  REAL :: h_peri = 0
  REAL :: Radius = 0

  INTEGER :: ncell = 0
  INTEGER, POINTER :: CellList(:)
  
  INTEGER :: nch = 0
  INTEGER, POINTER :: ChList(:)
  
END TYPE

TYPE NodeInfo_Type
  INTEGER :: irod
  !Node Information
  INTEGER :: ixnbeg, ixnend, iynbeg, iynend
  INTEGER :: ixnc, iync
END TYPE



INTEGER, PRIVATE :: nch = 0
INTEGER, PRIVATE :: nrod = 0 

INTEGER, PRIVATE :: nxc, nyc, nxy, nz
INTEGER, POINTER, PRIVATE :: CellMap(:, :)
INTEGER, POINTER, PRIVATE :: RodMap(:, :)

INTEGER, PRIVATE :: iRefPln

INTEGER :: nnx, nny
REAL, POINTER :: nodex(:)
REAL, POINTER :: nodey(:)
TYPE(NodeInfo_Type), POINTER :: NodeInfo(:, :)
INTEGER, POINTER :: ActiveNodeX(:, :)
INTEGER, POINTER :: ActiveNodeY(:, :)

TYPE(RodInfo_Type), POINTER, PRIVATE :: Rod(:)
TYPE(ChInfo_Type), POINTER, PRIVATE :: Channel(:)

TYPE(Cell_Type), POINTER, PRIVATE :: CellInfo(:)
TYPE(Pin_Type), POINTER, PRIVATE :: Pin(:)
TYPE(AsyInfo_Type), POINTER, PRIVATE :: AsyType(:)
CONTAINS

SUBROUTINE ConfigSubChGeom(Core, nTRACERCntl, PE)
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

Pin => Core%Pin; CellInfo => Core%CellInfo
AsyType => Core%AsyInfo
nxy = Core%nxy; nxc = Core%nx; nyc = Core%ny
nz = Core%nz

CALL CreateCellMap()

iRefPln = GetRefPln(Core%lFuelPlane)

CALL ConfigNode()

CALL ConfigRodInfo()

CALL ConfigSubChannel()


NULLIFY(Pin)
NULLIFY(CellInfo)

STOP

END SUBROUTINE

SUBROUTINE ConfigSubChannel() 
USE ALLOCS
IMPLICIT NONE

TYPE(ChInfo_Type), POINTER :: BaseCh(:)
INTEGER :: nBaseCh

ALLOCATE(BaseCh(nxy))
CALL ConfigBaseSubCh(BaseCh, nBaseCh)
CALL RmvPseudoSubCh(BaseCh, nBaseCh)
CALL ConfigBaseSubChNode(BaseCh, nBaseCh)

END SUBROUTINE 

SUBROUTINE RmvPseudoSubCh(BaseCh, nBaseCh)
IMPLICIT NONE
TYPE(ChInfo_Type) :: BaseCh(nxy) 
INTEGER :: nBaseCh
INTEGER :: isub, isub1, irod, i, j, k
INTEGER :: RodList(4), NewROdList(4), nrodch, nrodch0
INTEGER, POINTER :: RmvList(:)
INTEGER :: nrmv 

ALLOCATE(RmvList(nBaseCh))
CALL CP_CA(RmvList, 0, nBaseCh); nrmv = 0

DO isub = 1, nBaseCh
  nrodch = BaseCh(isub)%nrod
  nrodch0 = nrodch
  RodList = BaseCh(isub)%RodList
  DO i = 1, nrodch0
    irod = RodList(i)
    IF(Rod(irod)%lLGT) nrodch = nrodch - 1
  ENDDO
  IF(nrodch .EQ. nrodch0) CYCLE
  BaseCh(isub)%lLGT = .TRUE.
  IF(nrodch .EQ. 0) THEN
    BaseCh(isub)%lUse = .FALSE.
    nrmv = nrmv + 1
    RmvList(nrmv) = isub
    CYCLE
  ENDIF
  
  j = 0
  nrodch = 0; NewRodList = 0
  DO i = 1, nrodch0
    irod = RodList(i)
    IF(Rod(irod)%lLGT) j = j + 1
    IF(Rod(irod)%lLGT .AND. j .GT. 1) Cycle
    nrodch = nrodch + 1
    NewRodList(nrodch) = RodList(i)
  ENDDO 
  BaseCh(isub)%RodList = NewRodList
  BaseCh(isub)%nrod = nrodch
ENDDO

!Remove Subchannel Connectivity
DO k = 1, nrmv
  isub = RmvList(k)
  DO i = 1, 4
    isub1 = BaseCh(isub)%neighCh(i)
    IF(isub1 .EQ. 0) CYCLE
    DO j = 1, 4
      IF(BaseCh(isub1)%NeighCh(j) .NE. isub) CYCLE
      BaseCh(isub1)%NeighCh(j) = 0
      BaseCh(isub1)%Ngap = BaseCh(isub1)%Ngap -1
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE ConfigBaseSubCh(BaseCh, nbasech)
USE UtilFunction,    ONLY : CombineIntList
IMPLICIT NONE
TYPE(ChInfo_Type) :: BaseCh(nxy) 
INTEGER :: nbasech, nrod1d, nrod1d0, nrod1d1, YRangeMax

INTEGER :: ngx, ngy
INTEGER :: RodGridX(2*nxc), RodGridY(2*nyc), SubCh0(2,2)
INTEGER, POINTER :: SubChMap(:, :), RodGrid(:, :)
!INTEGER :: SubChMap(2*nxc, 2*nyc), RodGrid(0:2*nxc, 0:2*nyc)
INTEGER :: i, ix, iy, ixy, isub, insub, icel

LOGICAL :: lCent

INTEGER :: NeighDir(2,4), RodOd(2, 4)

Data NeighDir / 0, 1,  -1, 0,   0, -1,  1, 0/
DATA RodOd  /2, 2, 1, 2,  1, 1,  2, 1/ 
ALLOCATE(SubChMap(0:2*nxc, 0:2*nyc))
ALLOCATE(RodGrid(0:2*nxc, 0:2*nyc))

lCent = .FALSE.
icel = Pin(abs(CellMap(1, 1)))%Cell(irefPln)
lCent = CellInfo(icel)%lcentx .OR. CellInfo(icel)%lcentY .OR. CellInfo(icel)%lcentxy

CALL GetRodGrid(RodGrid, RodGridX, RodGridY, ngx, ngy)

nbasech = 0
CALL CP_CA(SubChMap, 0, 2*nxc+1, 2*nyc+1)
DO iy = 1, ngy+1
  IF(lCent .AND. iy .EQ. 1) CYCLE
  DO ix = 1, ngx+1
    IF(lCent .AND. ix .EQ. 1) CYCLE
    SubCh0 = RodGrid(ix-1:ix, iy-1:iy)
    IF(SUM(SubCh0) .EQ. 0) CYCLE
    nbasech = nbasech + 1
    SubChMap(ix, iy) = nbasech
    ALLOCATE(BaseCH(nbasech)%RodConf(2, 2))
    ALLOCATE(BaseCh(nBaseCh)%RodList(4))
    BaseCh(nbaseCh)%Rodconf = SubCh0
    BaseCh(nbaseCh)%nRod = 0
    BaseCh(nBaseCh)%RodList = 0
    DO i = 1, 4
      IF(SubCh0(RodOd(1, i), RodOd(2, i)) .EQ. 0) CYCLE
      BaseCh(nbaseCh)%nRod = BaseCh(nbaseCh)%nRod + 1
      BaseCh(nbasech)%RodList(BaseCh(nbaseCh)%nRod) = SubCh0(RodOd(1, i), RodOd(2, i))
    ENDDO
    
    BaseCh(nbaseCh)%lUse = .TRUE.
  ENDDO
ENDDO
DO iy = 1, ngy+1
  DO ix = 1, ngx + 1
    IF(SubChMap(ix, iy) .EQ. 0) CYCLE
    isub = SubChMap(ix, iy)
    BaseCh(isub)%ngap = 0
    ALLOCATE(BaseCh(isub)%NeighCh(4))
    BaseCh(isub)%NeighCh = 0
    DO i = 1, 4
      insub =SubChMap(ix+NeighDir(1, i), iy+NeighDir(2, i))
      IF(insub .EQ. 0) CYCLE
      BaseCh(isub)%ngap = BaseCh(isub)%ngap + 1
      BaseCh(isub)%NeighCh(i) = insub
    ENDDO
    CONTINUE
  ENDDO
ENDDO

DEALLOCATE(SubChMap, RodGrid)
END SUBROUTINE

SUBROUTINE GetRodGrid(RodGrid, RodGridX, RodGridY, ngx, ngy)
IMPLICIT NONE

INTEGER :: RodGrid(0:2*nxc, 0:2*nyc)
INTEGER :: RodGridX(2*nxc)
INTEGER :: ROdGridY(2*nyc)
INTEGER :: ngx, ngy
INTEGER :: ix, iy, igx, igy, icel

ngx = 0
DO ix = 1, nxc
  DO iy = 1, nyc
    IF(CellMap(ix, iy) .LT. 0) CYCLE
    icel  = CellMap(ix, iy)
    IF(RodMap(ix, iy) .EQ. 0) CYCLE
    ngx = ngx + 1
    RodGridX(ngx) = ix
    EXIT
  ENDDO
ENDDO

ngy = 0
DO iy = 1, nyc
  DO ix = 1, nxc
    IF(CellMap(ix, iy) .LT. 0) CYCLE
    icel  = CellMap(ix, iy)
    IF(RodMap(ix, iy) .EQ. 0) CYCLE
    ngy = ngy + 1
    ROdGridY(ngy) = iy
    EXIT
  ENDDO
ENDDO
CALL CP_CA(RodGrid, 0, 2*nxc+1, 2*nyc+1)
DO igy = 1, ngy
  DO igx = 1, ngx
    ix = RodGridx(igx); iy = RodGridy(igy)  
    IF(RodMap(ix, iy) .EQ. 0) CYCLE
    RodGrid(igx, igy) = RodMap(ix, iy)
  ENDDO
ENDDO
DO igx = 1, ngx
  RodGridx(igx) = NodeInfo(RodGridx(igx), 1)%ixnc
ENDDO
DO igy = 1, ngy
  RodGridy(igy) = NodeInfo(1, RodGridy(igy))%iync
ENDDO
END SUBROUTINE

SUBROUTINE ConfigRodInfo()
USE ALLOCS
USE Material_Mod
IMPLICIT NONE
INTEGER :: ixy, ix, iy, icel, irod
INTEGER :: ixbeg, ixend, iybeg, iyend
INTEGER :: i, j, k, ibeg

LOGICAL :: lHalf

ALLOCATE(Rod(nxy))
ALLOCATE(RodMap(nxc, nyc))
CALL CP_CA(RodMap, 0, nxc, nyc)
DO iy = 1, nyc
  DO ix = 1, nxc
    !ix = Pin(ixy)%ix; iy = Pin(ixy)%iy
    ixy = abs(CellMap(ix, iy))
    IF(CellMap(ix, iy) .LE. 0) CYCLE     ! Not in active core
    IF(RodMap(ix, iy) .GT. 0) CYCLE
    !IF(Pin(ixy)%lGap) CYCLE
    icel = Pin(ixy)%Cell(irefpln)
    IF(CellInfo(icel)%lGap) CYCLE
    !LARGE GUIDE Tube Check
    IF(Pin(ixy)%lGT .AND. CellInfo(icel)%lCCell) THEN
      CALL ProcLGT(ix, iy, ixy, icel, nrod, RodMap)
      Rod(nrod)%lGT = .TRUE.
      Rod(nrod)%lLGT = .TRUE.
      CYCLE
    ENDIF
    IF(Pin(ixy)%lGT) THEN
      nrod = nrod + 1
      RodMap(ix, iy) = nrod
      Rod(nrod)%nCell =  1
      Rod(nrod)%lGT = .TRUE.
      CYCLE
    ENDIF
  
    IF(Pin(ixy)%lFuel) THEN
      nrod = nrod + 1
      RodMap(ix, iy) = nrod
      Rod(nrod)%nCell =  1
      Rod(nrod)%lFUel = .TRUE.
    ENDIF
  ENDDO
ENDDO 
!ENDDO

DO irod = 1, nrod
  CALL Dmalloc(Rod(irod)%CellList, Rod(irod)%nCell)
   Rod(irod)%nCell = 0
ENDDO

!Set Rod Center Position
DO iy = 1, nyc
  DO ix = 1, nxc
    IF(RodMap(ix, iy) == 0) CYCLE
    irod = RodMap(ix, iy)
    ixy = abs(CellMap(ix, iy))
    Rod(irod)%nCell = Rod(irod)%nCell + 1
    Rod(irod)%CellList(Rod(irod)%nCell) = ixy
    
    IF(Rod(irod)%lUse) CYCLE
    icel = Pin(Ixy)%Cell(irefpln)
    
    lHalf = CellInfo(icel)%lCentX .OR. CellInfo(icel)%lCentY .OR. CellInfo(icel)%lCentXY
    IF(.NOT. CellInfo(icel)%lCCell .AND. .NOT. lHalf) THEN
      Rod(irod)%iync = NodeInfo(ix, iy)%iync
      Rod(irod)%ixnc = NodeInfo(ix, iy)%ixnc
    ELSEIF(lHalf) THEN
      IF(CellInfo(icel)%lCentX) THEN
        Rod(irod)%iync = NodeInfo(ix, iy)%iynbeg
        Rod(irod)%ixnc = NodeInfo(ix, iy)%ixnc        
      ELSEIF(CellInfo(icel)%lCentY) THEN
        Rod(irod)%iync = NodeInfo(ix, iy)%iync
        Rod(irod)%ixnc = NodeInfo(ix, iy)%ixnbeg        
      ELSE
        Rod(irod)%iync = NodeInfo(ix, iy)%iynbeg
        Rod(irod)%ixnc = NodeInfo(ix, iy)%ixnbeg          
      ENDIF
    ELSE
      SELECT CASE(CellInfo(icel)%Geom%CCentType)
        CASE(1)
          Rod(irod)%ixnc = NodeInfo(ix, iy)%ixnbeg
          Rod(irod)%iync = NodeInfo(ix, iy)%iynend
        CASE(2)
          Rod(irod)%ixnc = NodeInfo(ix, iy)%ixnbeg
          Rod(irod)%iync = NodeInfo(ix, iy)%iynbeg          
        CASE(3)
          Rod(irod)%ixnc = NodeInfo(ix, iy)%ixnend
          Rod(irod)%iync = NodeInfo(ix, iy)%iynbeg          
        CASE(4)
          Rod(irod)%ixnc = NodeInfo(ix, iy)%ixnend
          Rod(irod)%iync = NodeInfo(ix, iy)%iynend          
      END SELECT
    ENDIF
  ENDDO
ENDDO
!Rod Radius Search
DO irod = 1, nrod
  ixy = Rod(irod)%CellList(1)
  icel = Pin(ixy)%Cell(iRefPln)
  
  IF(CellInfo(icel)%lRect) THEN
    CALL TERMINATE('Can not Configure Sub-channel Geometry: Found Non-Annual Cell')
  ENDIF
  !Find Moderator
  DO i = 1, CellInfo(icel)%nFxr
    j = CellInfo(icel)%MapFxr2FsrIdx(1, i)
    k = CellInfo(icel)%iReg(j)
    IF(Mixture(k)%lh2o) EXIT
  ENDDO
  IF(i .GT. CellInfo(icel)%nFxr) CALL TERMINATE('Can not Configure Sub-channel Geometry : No Moderator Region in Cell')
  ibeg = i + 1
  DO i = ibeg, CellInfo(icel)%nFxr
    j = CellInfo(icel)%MapFxr2FsrIdx(1, i)
    k = CellInfo(icel)%iReg(j)
    IF(.NOT. Mixture(k)%lh2o) EXIT   
  ENDDO
  IF(i .GT. CellInfo(icel)%nFxr) CALL TERMINATE('Can not Configure Sub-channel Geometry : No Clad Region in Cell')
  i = i - 1
  Rod(irod)%Radius = CellInfo(icel)%Geom%Circle(3, i)
ENDDO

!DEALLOCATE(RodMap)
#define debug
#ifdef debug
DO irod = 1, nrod
  ix = Rod(irod)%ixnc; iy = Rod(irod)%iync
  Write(1001, '(I8, 3F12.4)') irod, nodeX(ix), nodeY(iy), Rod(irod)%Radius
ENDDO
#endif
END SUBROUTINE

SUBROUTINE CreateCellMap()
USE allocs
IMPLICIT NONE
INTEGER :: ixy, ix, iy, iasytype


CALL Dmalloc0(Cellmap, 0,nxc+1, 0, nyc+1)

DO ixy = 1, nxy
  ix = Pin(ixy)%ix; iy = Pin(ixy)%iy
  CellMap(ix, iy) = ixy
  
  iasytype = Pin(ixy)%AsyType
  IF(Pin(ixy)%lRadRef) CellMap(ix, iy) = - CellMap(ix, iy)

ENDDO


END SUBROUTINE

SUBROUTINE ConfigNode()
USE ALLOCS
IMPLICIT NONE

INTEGER :: ix, iy, ixn, iyn, i
INTEGER :: icel
REAL :: lx, ly
ALLOCATE(NodeInfo(nxc, nyc))
CALL DMALLOC0(Nodex, 0, 3*nxc)
CALL DMALLOC0(Nodey, 0, 3*nyc)
CALL DMALLOC0(ActiveNodeX, 1,2, 0, 3*nyc)
CALL DMALLOC0(ActiveNodeY, 1,2, 0, 3*nxc)
nnx = 0; nny = 0

Nodex(0) = 0; Nodey(0) = 0


DO ix = 1, nxc
  DO iy = 1, nyc
    IF(abs(CellMap(ix, iy)) .NE. 0) EXIT
  ENDDO
  IF(iy .GT. nyc) CALL TERMINATE('Error : Subchannel Configuration can not be generated ...')
  icel = abs(CellMap(ix, iy))
  icel = Pin(icel)%cell(irefpln)
  lx = CellInfo(icel)%Geom%lx
  nnx = nnx + 1
  NodeX(nnx) = NodeX(nnx-1) + lx/2._8
  nnx = nnx + 1
  NodeX(nnx) = NodeX(nnx-1) + lx/2._8
  
  DO iy = 1, nyc
    IF(abs(CellMap(ix, iy)) .EQ. 0) CYCLE
    NodeInfo(ix, iy)%ixnbeg = nnx - 2 
    NodeInfo(ix, iy)%ixnc = nnx - 1
    NodeInfo(ix, iy)%ixnend = nnx 
  ENDDO
ENDDO


DO iy = 1, nyc
  DO ix = 1, nxc
    IF(abs(CellMap(ix, iy)) .NE. 0) EXIT
  ENDDO
  IF(ix .GT. nxc) CALL TERMINATE('Error : Subchannel Configuration can not be generated ...')
  icel = abs(CellMap(ix, iy))
  icel = Pin(icel)%cell(irefpln)
  ly = CellInfo(icel)%Geom%ly
  nny = nny + 1
  NodeY(nny) = NodeY(nny-1) - ly/2._8
  nny = nny + 1
  NodeY(nny) = NodeY(nny-1) - ly/2._8
  
  DO ix = 1, nxc
    IF(abs(CellMap(ix, iy)) .EQ. 0) CYCLE
    NodeInfo(ix, iy)%iynbeg = nny - 2
    NodeInfo(ix, iy)%iync = nny - 1
    NodeInfo(ix, iy)%iynend = nny
  ENDDO
ENDDO
!Active Node Y
CALL CP_CA(ActiveNodeX, -1,  2, 3*nxc+1)
CALL CP_CA(ActiveNodeY, -1,  2, 3*nxc+1)
DO ix = 1, nxc
  i = 0
  DO iy = 1, nyc
    IF(CellMap(ix, iy) .LT. 1) CYCLE
    i = i + 1
  ENDDO
  IF(i .EQ. 0) CYCLE
  
  DO iy = 1, nyc
    IF(CellMap(ix, iy) .GT. 0) EXIT  
  ENDDO
  DO i = NodeInfo(ix, iy)%ixnbeg, NodeInfo(ix, iy)%ixnend
    ActiveNodeY(1, i) = NodeInfo(ix, iy)%iynbeg
  ENDDO
  
  DO iy = nyc, 1, -1
    IF(CellMap(ix, iy) .GT. 0) EXIT  
  ENDDO
  DO i = NodeInfo(ix, iy)%ixnbeg, NodeInfo(ix, iy)%ixnend
    ActiveNodeY(2, i) = NodeInfo(ix, iy)%iynend
  ENDDO
ENDDO

DO iy = 1, nyc
  i = 0
  DO ix = 1, nxc
    IF(CellMap(ix, iy) .LT. 1) CYCLE
    i = i + 1
  ENDDO
  IF(i .EQ. 0) CYCLE    
  
  DO ix = 1, nxc
    IF(CellMap(ix, iy) .GT. 0) EXIT  
  ENDDO  
  
  DO i = NodeInfo(ix, iy)%iynbeg, NodeInfo(ix, iy)%iynend
    ActiveNodeX(1, i) = NodeInfo(ix, iy)%ixnbeg  
  ENDDO
  
  DO ix = nxc, 1, -1
    IF(CellMap(ix, iy) .GT. 0) EXIT  
  ENDDO  
  
  DO i = NodeInfo(ix, iy)%iynbeg, NodeInfo(ix, iy)%iynend
    ActiveNodeX(2, i) = NodeInfo(ix, iy)%ixnend  
  ENDDO
ENDDO

CONTINUE
END SUBROUTINE

SUBROUTINE ProcLGT(ix, iy, ixy, icel, nrod, RodMap)
IMPLICIT NONE
INTEGER :: ix, iy, ixy, icel, nrod
INTEGER :: RodMap(nxc, nyc)
INTEGER :: CCentType
INTEGER :: GTDomX(2, 4)
INTEGER :: GTDomY(2, 4)
INTEGER :: ix0, iy0, ixbeg, ixend, iybeg, iyend

INTEGER :: nLGT, nOut

DATA GTDomX / -1,  0,   -1,  0,    0,  1,    0,  1 /
DATA GTDomY /  0,  1,   -1,  0,   -1,  0,    0,  1 /

CCentType = CellInfo(icel)%Geom%CCentType
ixbeg = ix + GTDomX(1, CCentType); ixend = ix + GTDomX(2, CCentType)
iybeg = iy + GTDomY(1, CCentType); iyend = iy + GTDomY(2, CCentType)

nOut = 0; nLGT = 0
DO iy0 = iybeg, iyend
  DO ix0 = ixbeg, ixend
    IF(CellMap(ix0, iy0) == 0) THEN
      nOut = nOut + 1
      CYCLE
    ENDIF
    icel = Pin(abs(CellMap(ix0, iy0)))%cell(iRefPln)
  
    IF(CellInfo(icel)%geom%lCCent) nLGT = nLGT + 1
  ENDDO
ENDDO
IF((nLGT + nOut) == 4) THEN
  nrod = nrod+ 1
  DO iy0 = iybeg, iyend
    DO ix0 = ixbeg, ixend
      IF(CellMap(ix0, iy0) == 0) CYCLE
      RodMap(ix0, iy0) = nrod
    ENDDO
  ENDDO  
  Rod(nrod)%ncell = nLGT
ENDIF

END SUBROUTINE

FUNCTION GetRefPln(FuelPlane)
IMPLICIT NONE
LOGICAL :: FuelPlane(nz)
INTEGER :: GetRefPln
INTEGER :: iz, iztop, izbot

DO iz = 1, nz
  IF(FuelPlane(iz)) THEN
    izbot = iz
    EXIT
  ENDIF
ENDDO

DO iz = nz, 1, -1
  IF(FuelPlane(iz)) THEN
    iztop = iz
    EXIT
  ENDIF
ENDDO

GetRefPln = iztop

END FUNCTION

SUBROUTINE ConfigBaseSubChNode(BaseCh, nBaseCh)
USE ALLOCS
IMPLICIT NONE
TYPE(ChInfo_Type) :: BaseCh(nxy) 
INTEGER :: nBaseCh
INTEGER :: isub, i
DO isub = 1, nBaseCh
  IF(.NOT. BaseCH(isub)%lUSe) CYCLE
  IF(.NOT. BaseCh(isub)%lLGT) THEN
    SELECT CASE(BaseCh(isub)%nRod)
      CASE(4)
        CALL SubChNode4(BaseCh(isub))
      CASE(3)
        CALL SubChNode3(BaseCh(isub))
      CASE(2)
        CALL SubChNode2(BaseCh(isub))
      CASE(1)
        CALL SubChNode1(BaseCh(isub))
    END SELECT
  ELSE
    SELECT CASE(BaseCh(isub)%nRod)
       CASE(4)
         CALL SubChNode4(BaseCh(isub))
       CASE(3)
         CALL SubChNode4(BaseCh(isub))
       CASE(2)
         CALL SubChNodeLGT2(BaseCh(isub))
       CASE(1)
         CALL TERMINATE('Error')
    END SELECT
  ENDIF
  
ENDDO
#define SubChDbg

#ifdef SubChDbg
DO isub = 1, nBasech
!DO isub = 74, 74
  IF(.NOT. BaseCh(isub)%lUse) CYCLE
  IF(BaseCh(isub)%nNode .EQ. 0) CYCLE
  WRITE(1004, '(I5, I10, (10000F14.5))') isub,  BaseCh(isub)%nNode, (NodeX(BaseCh(isub)%NodeList(1, i)), NodeY(BaseCh(isub)%NodeList(2, i)), i = 1, BaseCh(isub)%nNode)
ENDDO

#endif

END SUBROUTINE

SUBROUTINE SubChNode4(BaseCh)
USE Allocs
IMPLICIT NONE
TYPE(ChInfo_Type) :: BaseCh

INTEGER :: RodOd(2, 4)
INTEGER :: NodeList(2, 100)
INTEGER :: i, j, ix, iy, irod
DATA RodOd     / 2, 2,   1, 2,  1, 1,   2, 1/
j = 0
DO i = 1, 4
  irod = BaseCh%RodConf(RodOd(1, i), RodOd(2, i))
  IF(irod .EQ. 0) CYCLE
  j = j + 1
  NodeList(:, j) = (/Rod(irod)%ixnc, Rod(irod)%iync/)
ENDDO

BaseCh%nNode = j
CALL DMALLOC(BaseCh%NodeList, 2, j)
CALL CP_VA(BaseCh%NodeList, NodeList, 2, j)

END SUBROUTINE

SUBROUTINE SubChNode1(BaseCh)
USE Allocs
IMPLICIT NONE
TYPE(ChInfo_Type) :: BaseCh
INTEGER :: RodOd(2, 4)
INTEGER :: NodeList(2, 100), CaseDat(2,4)
INTEGER :: CaseID
INTEGER :: i, j, ixn, iyn, irod

DATA RodOd     / 2, 2,   1, 2,  1, 1,  2, 1/

DO i = 1, 4
  irod = BaseCh%RodConf(RodOd(1, i), RodOd(2, i))
  IF(irod .NE. 0) EXIT
ENDDO
CaseId = i
j = 1

ixn = Rod(irod)%ixnc; iyn = Rod(irod)%iync
NodeList(:, j) = (/ixn, iyn/)
CaseDat(:, 1) = (/ixn, ActiveNodeY(2, ixn)/); CaseDat(:,3) = (/ixn, ActiveNodeY(1, ixn)/)
CaseDat(:,2) = (/ActiveNodeX(1, iyn), iyn/); CaseDat(:,4) = (/ActiveNodeX(2, iyn), iyn/)


Select Case(CaseId)
  CASE(1)
    j = 4
    NodeList(:, 2) = CaseDat(:, 2)
    NodeList(:, 4) = CaseDat(:, 3)
    NodeList(:, 3) = (/CaseDat(1, 2), CaseDat(2, 3)/)
  CASE(2)
    j = 4
    NodeList(:, 2) = CaseDat(:, 3)
    NodeList(:, 4) = CaseDat(:, 4)
    NodeList(:, 3) = (/CaseDat(1, 4), CaseDat(2, 3)/)
  CASE(3)
    j = 4
    NodeList(:, 2) = CaseDat(:, 4)
    NodeList(:, 4) = CaseDat(:, 1)
    NodeList(:, 3) = (/CaseDat(1, 4), CaseDat(2, 1)/)    
  CASE(4)
    j = 4
    NodeList(:, 2) = CaseDat(:, 1)
    NodeList(:, 4) = CaseDat(:, 2)
    NodeList(:, 3) = (/CaseDat(1, 2), CaseDat(2, 1)/)      
END SELECT
BaseCh%nNode = j
CALL DMALLOC(BaseCh%NodeList, 2, j)
CALL CP_VA(BaseCh%NodeList, NodeList, 2, j)
END SUBROUTINE


SUBROUTINE SubChNode2(BaseCh)
USE Allocs
IMPLICIT NONE
TYPE(ChInfo_Type) :: BaseCh
INTEGER :: RodOd(2, 4), NeighRod(2,4), CaseType(2, 4), CaseDat(2, 4)
INTEGER :: NodeList(2, 100)
INTEGER :: i, j, k, irod, inrod, icase, ixn, iyn


DATA RodOd     / 2, 2,   1, 2,  1, 1,  2, 1/
DATA NeighRod  / 2, 4,   1, 3,  4, 2,  3, 1/
DATA CASETYPE  / 4, 1,   2, 1,  2, 3,  4, 3/
j = 0
DO i = 1, 4
  irod = BaseCh%RodConf(RodOd(1, i), RodOd(2, i))
  IF(irod .NE. 0) THEN
    j = j + 1
    NodeList(:, j) = (/ Rod(irod)%ixnc, Rod(irod)%iync /)
    CYCLE
  ENDIF
  DO k = 1, 2
    inrod = BaseCh%RodConf(RodOd(1, NeighRod(k, i)), RodOd(2, NeighRod(k, i)))
    IF(inrod .EQ. 0) CYCLE
    icase = CaseType(k, i)
    ixn = Rod(inrod)%ixnc; iyn = Rod(inrod)%iync
    CaseDat(:, 1) = (/ixn, ActiveNodeY(2, ixn)/); CaseDat(:,3) = (/ixn, ActiveNodeY(1, ixn)/)
    CaseDat(:,2) = (/ActiveNodeX(1, iyn), iyn/); CaseDat(:,4) = (/ActiveNodeX(2, iyn), iyn/)
    j = j + 1
    NodeList(:, j) = CaseDat(:, icase)
    EXIT
  ENDDO
ENDDO

BaseCh%nNode = j
CALL DMALLOC(BaseCh%NodeList, 2, j)
CALL CP_VA(BaseCh%NodeList, NodeList, 2, j)
  !INTEGER :: nNode = 0
  !INTEGER, POINTER :: NodeList(:)
END SUBROUTINE


SUBROUTINE SubChNode3(BaseCh)
USE Allocs
IMPLICIT NONE
TYPE(ChInfo_Type) :: BaseCh
INTEGER :: RodOd(2, 4), Dat(2, 3)
INTEGER :: NodeList(2, 100)
INTEGER :: CaseId
INTEGER :: i, j, k, ixn, iyn, irod
DATA RodOd     / 2, 2,   1, 2,  1, 1,  2, 1/
DO i = 1, 4
  irod = BaseCh%RodConf(RodOd(1, i), RodOd(2, i))
  IF(irod .EQ. 0) EXIT
ENDDO
CaseId = i

SELECT CASE(CaseId)
CASE(1)
    irod = BaseCh%RodConf(RodOd(1, 4), RodOd(2, 4))
    ixn = Rod(irod)%ixnc; iyn = Rod(irod)%iync
    Dat(:, 1) = (/ixn, ActiveNodeY(2, ixn)/)
    irod = BaseCh%RodConf(RodOd(1, 2), RodOd(2, 2))
    ixn = Rod(irod)%ixnc; iyn = Rod(irod)%iync
    Dat(:, 3) = (/ActiveNodeX(2, iyn), iyn/)
    Dat(:, 2) = (/Dat(1, 3), Dat(2, 1)/)
  CASE(2)
    irod = BaseCh%RodConf(RodOd(1, 1), RodOd(2, 1))
    ixn = Rod(irod)%ixnc; iyn = Rod(irod)%iync
    Dat(:, 1) = (/ActiveNodeX(1, iyn) , iyn/)
    irod = BaseCh%RodConf(RodOd(1, 3), RodOd(2, 3))
    ixn = Rod(irod)%ixnc; iyn = Rod(irod)%iync
    Dat(:, 3) = (/ixn, ActiveNodeY(2, ixn)/)
    Dat(:, 2) = (/Dat(1, 1), Dat(2, 3)/)
  CASE(3)
    irod = BaseCh%RodConf(RodOd(1, 2), RodOd(2, 2))
    ixn = Rod(irod)%ixnc; iyn = Rod(irod)%iync
    Dat(:, 1) = (/ixn, ActiveNodeY(1, ixn)/)
    irod = BaseCh%RodConf(RodOd(1, 4), RodOd(2, 4))
    ixn = Rod(irod)%ixnc; iyn = Rod(irod)%iync
    Dat(:, 3) = (/ActiveNodeX(1, iyn), iyn/)
    Dat(:, 2) = (/Dat(1, 3), Dat(2, 1)/)     
  CASE(4)
    irod = BaseCh%RodConf(RodOd(1, 3), RodOd(2, 3))
    ixn = Rod(irod)%ixnc; iyn = Rod(irod)%iync
    Dat(:, 1) = (/ActiveNodeX(2, iyn), iyn/)
    irod = BaseCh%RodConf(RodOd(1, 1), RodOd(2, 1))
    ixn = Rod(irod)%ixnc; iyn = Rod(irod)%iync
    Dat(:, 3) = (/ixn, ActiveNodeY(1, ixn)/)
    Dat(:, 2) = (/Dat(1, 1), Dat(2, 3)/)
END SELECT

j= 0
DO i = 1, 4
  irod = BaseCh%RodConf(RodOd(1, i), RodOd(2, i))
  IF(irod .NE. 0) THEN
    j = j + 1
    NodeList(:, j) = (/ Rod(irod)%ixnc, Rod(irod)%iync /)
    CYCLE
  ENDIF
  DO k = 1, 3
    j = j + 1  
    NodeList(:, j) = Dat(:, k)
  ENDDO
  
ENDDO
BaseCh%nNode = j
CALL DMALLOC(BaseCh%NodeList, 2, j)
CALL CP_VA(BaseCh%NodeList, NodeList, 2, j)
END SUBROUTINE

SUBROUTINE SubChNodeLGT2(BaseCh)
USE Allocs
IMPLICIT NONE
TYPE(ChInfo_Type) :: BaseCh
INTEGER :: RodOd(2, 4), Dat(2, 3)
INTEGER :: NodeList(2, 100), LGTNode(2), RodNode(2), CaseDat(2, 2)
INTEGER :: CaseId
INTEGER :: i, j, k, irod
REAL :: L(2), x1, y1, x2, y2
DATA RodOd     / 2, 2,   1, 2,  1, 1,  2, 1/

DO i = 1, 4
  irod = BaseCh%RodConf(RodOd(1, i), RodOd(2, i))
  IF(irod .EQ. 0) CYCLE
  IF(Rod(irod)%LLGT) THEN
    LGTNode = (/Rod(irod)%ixnc, Rod(irod)%iync/)
  ELSE
    RodNode = (/Rod(irod)%ixnc, Rod(irod)%iync/)
  ENDIF
ENDDO

CaseDat(:, 1) = (/RodNode(1), ActiveNodeY(1, RodNode(1))/)
CaseDat(:, 2) = (/ActiveNodeX(1, RodNode(2)), RodNode(2)/)

x1 = NodeX(RodNode(1)); Y1 = NodeY(RodNode(2))
DO i = 1, 2
  x2 = NodeX(CaseDat(1, i)); y2 = NodeY(CaseDat(2,i))
  l(i) = SQRT((x2-x1)**2 + (y2-y1)**2)
ENDDO

i = 2
IF(l(1) .LT. l(2)) i = 1  

j= 3
NodeList(:, 1) = LGTNode
NodeList(:, 2) = RodNode
NodeList(:, 3) = CaseDat(:, i)

BaseCh%nNode = j
CALL DMALLOC(BaseCh%NodeList, 2, j)
CALL CP_VA(BaseCh%NodeList, NodeList, 2, j)

END SUBROUTINE

END MODULE

