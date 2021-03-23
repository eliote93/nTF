SUBROUTINE GlobalCoreMap()
USE PARAM
USE ALLOCS
!USE CNTL, ONLY : master
USE GEOM, ONLY : nAsyType,   nPinType,    nCellType,    nPinType0,     nCellType0,         &  
                 Core,       AsyInfo,     PinInfo,      CellInfo,      Pin,                &
                 Asy,        hz,          nz,           nzfm,                              &
                 nx,         ny,          nxy,                                             &
                 hzInv,      HzFmInv,                                                      &
                 nCellX0,    nCellX,                                                       &
                 AsyPitch,   lEdge,       lGap,          lRot,         lCbd,                    &
                 hzfm,       SubPlaneMap, SubPlaneRange, nSubPlane
USE BasicOperation, ONLY : CP_CA,     CP_VA
USE CNTL,           ONLY : nTRACERCntl
USE PE_MOD,         ONLY : PE            
IMPLICIT NONE
INTEGER, PARAMETER :: NBd = 4
INTEGER :: nx0,ny0
INTEGER :: nxa, nya, nxya
INTEGER :: nxc, nyc, nxyc, nxc0, nyc0
INTEGER :: ix,  iy, ixy, ixbeg, ixend, iybeg, iyend
!INTEGER :: myzb, myze, myzbf,  myzef
INTEGER :: ixa, iya, ixya
INTEGER :: iasy, ipin, icel, iz
INTEGER :: PinType, AsyType

INTEGER, POINTER :: GlobalAsyIdx(:,:), GlobalPinMap(:,:), GlobalPinIdx(:,:), LocalPinIdx(:,:), &
                    GlobalAsyMap(:,:)
INTEGER :: Asy_NX(500)
REAL :: lx, ly
nxy=0; nx=0; ny=0
nxa = Core%nxa; nya = Core%nya; nxya = Core%nxya
Core%nz = nz; Core%nzfm =nzfm
Core%nPinType = nPinType; Core%nCellType = nCellType; Core%nCellType0 = nCellType0
Core%nAsyType = nAsyType

Core%nxc0 = nCellX0; Core%nyc0 = nCellX0; 
Core%nxyc0 = nCellX0 * nCellX0

Core%nxc = nCellX; Core%nyc = nCellX; 
Core%nxyc = nCellX * nCellX

Core%lGap = lGap; Core%lEdge = lEdge; Core%lRot = lRot; Core%lCbd = lCbd

!myzb = PE%myzb; myze = PE%myze
!myzb = 1; myze = nz
!Check Used Geom
CALL SetUsedGeom()

!Set Up Total Number of Pincell
DO ixya = 1, nxya
  iasy = Core%CoreMap(ixya)
  nxy = nxy + AsyInfo(iasy)%nxy
ENDDO
Asy_NX(1:nxa) = 0
DO iya = 1, nya
  nx0 = 0
  DO ixa = 1, nxa
    ixya = Core%CoreIdx(ixa, iya)
    IF(ixya .EQ. 0) CYCLE
    iasy = Core%CoreMap(ixya)
    nx0 = nx0 + AsyInfo(iasy)%nx
    Asy_Nx(ixa) = MAX(Asy_Nx(ixa),  AsyInfo(iasy)%nx)
  ENDDO
  nx = max(nx0, nx)
ENDDO

DO ixa = 1, nxa
  ny0 = 0
  DO iya = 1, nya
    ixya = Core%CoreIdx(ixa, iya)
    IF(ixya .EQ. 0) CYCLE
    iasy = Core%CoreMap(ixya)
    ny0 = ny0 + AsyInfo(iasy)%ny
  ENDDO
  ny = max(ny0, ny)
ENDDO
Core%nx = nx; Core%ny = ny; Core%nxy = nxy

!Set Up Global Assembly Map
ALLOCATE(Asy(nxya))
ALLOCATE(GlobalAsyMap(0:nxa+1, 0:nya+1))
CALL CP_CA(GlobalAsyMap(0:nxa+1, 0:nya+1), 0, nxa+2, nya+2)

Core%nFuelPin = 0
DO iya = 1, nya
  DO ixa = 1, nxa
    ixya = Core%CoreIdx(ixa, iya)
    GlobalAsyMap(ixa, iya) = ixya
    IF(ixya .EQ. 0) CYCLE
    iasy = Core%CoreMap(ixya)  !Asy Index
    Asy(ixya)%ixa = ixa; Asy(ixya)%iya = iya  !Assembly Location
    Asy(ixya)%AsyType = iasy  !Assembly Type
    Core%nFuelPin = Core%nFuelPin + AsyInfo(iasy)%nFuelPin  !Number of Fuel Pin
    !Partial Assembly Information
    Asy(ixya)%lcentX = AsyInfo(iasy)%lCentX; Asy(ixya)%lcentY = AsyInfo(iasy)%lCentY; 
    Asy(ixya)%lcentXY = AsyInfo(iasy)%lCentXY;
    Asy(ixya)%PartialAsyFlag = 0
    IF(Asy(ixya)%lcentX) Asy(ixya)%PartialAsyFlag = 1
    IF(Asy(ixya)%lcentY) Asy(ixya)%PartialAsyFlag = 2
    IF(Asy(ixya)%lcentXY) Asy(ixya)%PartialAsyFlag = 3
    Asy(ixya)%wt = 1.0_8
    SELECT CASE(Asy(ixya)%PartialAsyFlag) 
      CASE(1, 2)
        Asy(ixya)%wt = 0.5_8  
      CASE(3)
        Asy(ixya)%wt = 0.25_8
    END SELECT
    CAll Dmalloc(Asy(ixya)%GlobalPinIdx, AsyInfo(Asy(ixya)%AsyType)%nxy)
    !Assembly Center 
  ENDDO
ENDDO

CALL AsyCenterCoordinate()
!Assembly Neighborhood Information
DO iya = 1, nya
  DO ixa = 1, nxa
    ixya = Core%CoreIdx(ixa, iya)
    IF(ixya .EQ. 0) CYCLE
    CALL DMALLOC(Asy(ixya)%NeighIdx, 4)
    Asy(ixya)%NeighIdx(SOUTH) = GlobalAsyMap(ixa, iya + 1)
    Asy(ixya)%NeighIdx(WEST) = GlobalAsyMap(ixa - 1, iya)
    Asy(ixya)%NeighIdx(NORTH) = GlobalAsyMap(ixa, iya - 1)
    Asy(ixya)%NeighIdx(EAST) = GlobalAsyMap(ixa + 1, iya)  
  ENDDO
ENDDO
Deallocate(GlobalAsyMap)


ALLOCATE(Pin(nxy))
ALLOCATE(GlobalAsyIdx(0:nx+1, 0:ny+1))
ALLOCATE(GlobalPinMap(0:nx+1, 0:ny+1))
ALLOCATE(GlobalPinIdx(-nCellX:nx+1, -nCellX:ny+1))
ALLOCATE(LocalPinIdx(0:nx+1, 0:ny+1))

CALL CP_CA(GlobalAsyIdx, 0, nx + 2, ny + 2)
CALL CP_CA(GlobalPinMap, 0, nx + 2, ny + 2)
CALL CP_CA(GlobalPinIdx, 0, nx + nCellX + 2, ny + nCellX + 2)
CALL CP_CA(LocalPinIdx, 0, nx + 2, ny + 2)


!Set Up Global Core  Pin Map
iyend = 0
ipin =0
ixy = 0
DO iya = 1, nya
  ixend = 0; iybeg = iyend 
  DO ixa = 1, nxa
    ixya = Core%CoreIdx(ixa, iya)
    IF(ixya .EQ. 0) THEN
      ixend = ixend + Asy_Nx(ixa)
      CYCLE
    ENDIF
    iasy = Core%CoreMap(ixya)
    nxc = AsyInfo(iasy)%nx; nyc = AsyInfo(iasy)%ny; nxyc = AsyInfo(iasy)%nxy
    ixbeg = ixend
    DO iy = 1, nyc
      DO ix = 1, nxc
        ixy = ixy + 1
        ipin = AsyInfo(iasy)%Pin(AsyInfo(iasy)%Pin2DIdx(ix, iy))
        GlobalAsyIdx(ixbeg + ix, iybeg + iy) = ixya
        GlobalPinMap(ixbeg + ix, iybeg + iy) = ipin
        GlobalPinIdx(ixbeg + ix, iybeg + iy) = ixy
        LocalPinIdx(ixbeg + ix, iybeg + iy) = AsyInfo(iasy)%Pin2DIdx(ix, iy);
      ENDDO
    ENDDO
    ixend = ixend + nxc
  ENDDO
  iyend = iyend + nyc
ENDDO

!Set Boundary Condition for the Global Pin Map
IF(Core%RadBC(SOUTH) .EQ. RefCell) THEN
  CALL CP_CA(GlobalPinIdx(1 : nx, ny + 1), RefCell, nx)
ENDIF
IF(Core%RadBC(WEST) .EQ. RefCell) THEN
  CALL CP_CA(GlobalPinIdx(0, 1 : ny), RefCell, ny)
ENDIF
IF(Core%RadBC(NORTH) .EQ. RefCell) THEN
  CALL CP_CA(GlobalPinIdx(1 : nx, 0), RefCell, nx)
ENDIF
IF(Core%RadBC(EAST) .EQ. RefCell) THEN
  CALL CP_CA(GlobalPinIdx(nx + 1, 1 : ny), RefCell, ny)
ENDIF
!Rotational Boundary Condition

IF(lROT) THEN
  CALL CP_VA(GlobalPinIdx(1:nx, 0), GlobalPinIdx(1,1:nx), nx)
  CALL CP_VA(GlobalPinIdx(0, 1:ny), GlobalPinIdx(1:ny, 1), ny)
ELSEIF(lCBD) THEN
  CALL CP_VA(GlobalPinIdx(1:nx, 0), GlobalPinIdx(1:nx, ny), nx)
  CALL CP_VA(GlobalPinIdx(1:nx, ny+1), GlobalPinIdx(1:nx, 1), nx)
  CALL CP_VA(GlobalPinIdx(0, 1:ny), GlobalPinIdx(nx, 1:ny), ny)
  CALL CP_VA(GlobalPinIdx(nx+1, 1:ny), GlobalPinIdx(1, 1:ny), ny)
ENDIF

!Assign the Pin Variable data
ixy =0
DO iy = 1, ny
  DO ix = 1, nx
    PinType = GlobalPinMap(ix, iy)
    iasy = GlobalAsyIdx(ix, iy)
    IF(PinType .eq. 0)  CYCLE
            
    ixy = GlobalPinIdx(ix, iy)
    AsyType = Core%CoreMap(iasy)
    Pin(ixy)%Pintype = PinType; Pin(ixy)%AsyType = Core%CoreMap(iasy)
    Pin(ixy)%ipin = LocalPinIdx(ix, iy); Pin(ixy)%iasy = iasy
    Pin(ixy)%nCell = PinInfo(PinType)%nCell
    Pin(ixy)%nBd = nBd
    Pin(ixy)%nFsrMax = PinInfo(PinType)%nFsrMax
    Pin(ixy)%nFxrMax = PinInfo(PinType)%nFxrMax

    CALL Dmalloc(Pin(ixy)%Cell, Pin(ixy)%nCell)
    CAll Dmalloc(Pin(ixy)%BdLength, nBd)
    CALL Dmalloc(Pin(ixy)%Center2SurfaceL, nBd)
    CALL Dmalloc(Pin(ixy)%NeighIdx, nBd)
    CALL Dmalloc(Pin(ixy)%NeighSurfIdx, nBd)
    
    Pin(ixy)%ix = ix; Pin(ixy)%iy = iy
    
    Pin(ixy)%Cell = PinInfo(PinType)%Cell
    Pin(ixy)%NeighIdx(SOUTH) = GlobalPinIdx(ix, iy + 1)
    Pin(ixy)%NeighIdx(WEST) = GlobalPinIdx(ix - 1, iy)
    Pin(ixy)%NeighIdx(NORTH) = GlobalPinIdx(ix, iy - 1)
    Pin(ixy)%NeighIdx(EAST) = GlobalPinIdx(ix + 1, iy)
    
    Pin(ixy)%NeighSurfIdx(SOUTH) = NORTH
    Pin(ixy)%NeighSurfIdx(NORTH) = SOUTH
    Pin(ixy)%NeighSurfIdx(EAST) = WEST
    Pin(ixy)%NeighSurfIdx(WEST) = EAST
    !Pin - Dimension 
    lx = CellInfo(Pin(ixy)%Cell(1))%Geom%lx
    ly = CellInfo(Pin(ixy)%Cell(1))%Geom%ly
    Pin(ixy)%BdLength(SOUTH) = lx; Pin(ixy)%BdLength(NORTH) = lx
    Pin(ixy)%BdLength(WEST) = ly; Pin(ixy)%BdLength(EAST) = ly
    !Center to Surface Distance Information
    Pin(ixy)%Center2SurfaceL(SOUTH) = ly*HALF; Pin(ixy)%Center2SurfaceL(NORTH) = ly*HALF
    Pin(ixy)%Center2SurfaceL(WEST) = lx*HALF; Pin(ixy)%Center2SurfaceL(EAST) = lx*HALF
    !
    Asy(iasy)%GlobalPinIdx(LocalPinIdx(ix, iy)) = ixy
  ENDDO
ENDDO

!Rotational boundary neighboring surface index mapping
IF(lROT) THEN
  DO ix = 1, nx
    ixy = GlobalPinIdx(ix, 1)
    Pin(ixy)%NeighSurfIdx(NORTH) = WEST
  ENDDO
  DO iy = 1, ny
    ixy = GlobalPinIdx(1, iy)
    Pin(ixy)%NeighSurfIdx(WEST) = NORTH
  ENDDO
ELSEIF(lCBD) THEN
  DO ix = 1, nx
    ixy = GlobalPinIdx(ix, 1)
    Pin(ixy)%NeighSurfIdx(NORTH) = SOUTH
    ixy = GlobalPinIdx(ix, ny)
    Pin(ixy)%NeighSurfIdx(SOUTH) = NORTH
  ENDDO  
  DO iy = 1, ny
    ixy = GlobalPinIdx(1, iy)
    Pin(ixy)%NeighSurfIdx(WEST) = EAST
    ixy = GlobalPinIdx(nx, iy)
    Pin(ixy)%NeighSurfIdx(EAST) = WEST
  ENDDO
ENDIF

IF(.NOT. lEdge) THEN
  iybeg = 1; ixbeg = 1
  icel = Pin(1)%Cell(1)
  IF(CellInfo(icel)%lCentXY) THEN
    ixbeg = 2; iybeg =2
  ENDIF
  nxc0= nxc-AsyInfo(Asy(1)%AsyType)%nx
  nyc0= nyc-AsyInfo(Asy(1)%AsyType)%ny
  IF(lRot) THEN
    DO iy = 0, -nyc0+1, -1
      CALL CP_VA(GlobalPinIdx(1:nx, iy), GlobalPinIdx(ixbeg, 1:nx), nx)
      ixbeg = ixbeg + 1
    ENDDO
    CONTINUE
    DO ix = 0, -nxc0 + 1, -1
      CALL CP_VA(GlobalPinIdx(ix, 1:ny), GlobalPinIdx(1:ny, iybeg), ny)
      iybeg = iybeg + 1
    ENDDO
    
    DO ix = -nxc0 + 1, 0
      iybeg = iybeg -1
      CALL CP_VA(GlobalPinIdx(ix, -nyc0+1:0), GlobalPinIdx(-nyc0+1:0, iybeg), nyc0)
    ENDDO
    CONTINUE   
  ELSE
    DO iy = 0, -nyc0+1, -1
      CALL CP_VA(GlobalPinIdx(1:nx, iy), GlobalPinIdx(1:nx, iybeg), nx)
      iybeg = iybeg + 1
    ENDDO
    CONTINUE
    DO ix = 0, -nxc0 + 1, -1
      CALL CP_VA(GlobalPinIdx(ix, 1:ny), GlobalPinIdx(ixbeg, 1:ny), ny)
      ixbeg = ixbeg + 1
    ENDDO
    
    DO ix = -nxc0 + 1, 0
      ixbeg = ixbeg -1
      CALL CP_VA(GlobalPinIdx(ix, -nyc0+1:0), GlobalPinIdx(ixbeg, -nyc0+1:0), nyc0)
    ENDDO   
  ENDIF
ENDIF


Deallocate(GlobalAsyIdx, GlobalPinMap, GlobalPinIdx, LocalPinIdx)


END SUBROUTINE


SUBROUTINE SetCoreInfo()
USE PARAM
USE ALLOCS
!USE CNTL, ONLY : master
USE GEOM, ONLY : nAsyType,   nPinType,    nCellType,    nPinType0,     nCellType0,         &  
                 Core,       AsyInfo,     PinInfo,      CellInfo,      Pin,                &
                 Asy,        hz,          nz,           nzfm,          BaseCellInfo,       &
                 hzInv,      HzFmInv,                                                      &
!                 myzb,       myze,        myzbf,        myzef,                            &
                 hzfm,       SubPlaneMap, SubPlaneRange, nSubPlane,   nSubPlanePtr          
USE PE_MOD,         ONLY : PE            
IMPLICIT NONE 
INTEGER :: ixy, iz, icel, ipin                  
Core%Pin => Pin
Core%Asy => Asy
Core%AsyInfo => AsyInfo
Core%PinInfo => PinInfo
Core%CellInfo => CellInfo
Core%BaseCellInfo => BaseCellInfo ! --- 180724 JSR
Core%hz => hz
Core%hzInv => hzInv

!SubPlane Information
Core%hzfm => hzfm
Core%hzfmInv => hzfmInv
Core%nSubPlane = nSubPlane
Core%nSubPlanePtr => nSubPlanePtr
Core%SubPlaneMap => SubPlaneMap
Core%SubPlaneRange => SubPlaneRange
ALLOCATE(Core%lFuelPlane(nz)); Core%lFuelPlane = FALSE
ALLOCATE(Core%lCladPlane(nz)); Core%lCladPlane = FALSE
ALLOCATE(Core%lAICPlane(nz)); Core%lAICPlane = FALSE
DO ixy = 1, Core%nxy
  ALLOCATE(Pin(ixy)%lMOX(Core%nz))
  DO iz =1, Core%nz
    icel = Pin(ixy)%Cell(iz)
    
    Core%lFuelPlane(iz) =  Core%lFuelPlane(iz) .OR. CellInfo(icel)%lfuel
    Core%lAICPlane(iz) =  Core%lAICPlane(iz) .OR. CellInfo(icel)%lAIC
    Pin(ixy)%lfuel = Pin(ixy)%lfuel .or. CellInfo(icel)%lfuel
    Pin(ixy)%lGd = Pin(ixy)%lGd .or. CellInfo(icel)%lGd
    
    Pin(ixy)%lMox(iz) = CellInfo(icel)%lMox
  ENDDO
ENDDO

CALL CoreVolGen()
                                   
END SUBROUTINE

SUBROUTINE SetGlobalFsr_Fxr()
USE PARAM
USE CNTL, ONLY : nTracerCntl
USE GEOM, only : core,  Pin
USE PE_MOD, ONLY : PE 
IMPLICIT NONE
INTEGER :: nCoreFsr, nCoreFxr
INTEGER :: ixy, nxy
nCoreFsr = 0; nCoreFxr = 0
nxy = Core%nxy
DO ixy = 1, nxy
  Pin(ixy)%FsrIdxSt = nCoreFsr + 1
  Pin(ixy)%FxrIdxSt = nCoreFxr + 1
  nCoreFsr = nCoreFsr + Pin(ixy)%nFsrMax
  nCoreFxr = nCoreFxr + Pin(ixy)%nFxrMax
ENDDO

Core%nCoreFsr = nCoreFsr
Core%nCoreFxr = nCoreFxr
!Prepare the Global FXR information
CALL PrepFxr(nTracerCntl%lXsLib, nTracerCntl%lfxrlib)

END SUBROUTINE

SUBROUTINE SetUsedGeom()
USE PARAM
USE GEOM,    ONLY : CellInfo,      PinInfo,    Pin,    AsyInfo,    CORE,    &
                    nAsyType,      nPinType,   nCellType
USE IOUTIL,  ONLY : TERMINATE
IMPLICIT NONE                    
INTEGER :: ipin, icel, iasy, ixy, ixya, iz
INTEGER :: nUseAsy, nUsePin, nUseCell
INTEGER :: nxya, nxy, nz
LOGICAL :: lCENT, lCCELL

!SET USED Assembly
nxya = CORE%nxya
DO ixy = 1, nxya
  iasy = Core%CoreMap(ixy)
  AsyInfo(iasy)%luse = TRUE
ENDDO

!SET USED PIN
nUseAsy = 0
DO iasy = 1, nASyType
  IF(.NOT. AsyInfo(iasy)%lUse) CYCLE
  nUseAsy = nUseAsy + 1
  nxy = AsyInfo(iasy)%nxy
  DO ixy = 1, nxy
    ipin = ASyInfo(iasy)%Pin(ixy)
    PinInfo(ipin)%luse = TRUE
  ENDDO
ENDDO

nUsePin = 0
DO ipin = 1, nPinType
  IF(.NOT. PinInfo(ipin)%luse) CYCLE
  nUsePin = nUsePin + 1
  nz = PinInfo(ipin)%nCell
  DO iz = 1, nz
    icel = PinInfo(ipin)%Cell(iz)
    CellInfo(icel)%luse = TRUE
    lCENT = CellInfo(icel)%lCENTX .OR. CellInfo(icel)%lCENTY .OR. CellInfo(icel)%lCENTXY
    lCCELL = CellInfo(icel)%lCCell
    IF(lCCell .AND. lCent) CALL TERMINATE('Error!!: Corner Cell can not be located in symmetirc boundary line.')
  ENDDO
ENDDO

nUseCell = 0
DO iCel = 1, nCellType
  IF(.NOT. CellInfo(icel)%luse) CYCLE
  nUseCell = nUseCell + 1
ENDDO

END SUBROUTINE


SUBROUTINE AsyCenterCoordinate()
USE PARAM
USE TYPEDEF, ONLY : coreinfo_type,  Asy_Type
USE GEOM,    ONLY : Core,           Asy,        lEdge,        AsyPitch
USE ALLOCS
IMPLICIT NONE
REAL :: CX(500), CY(500)
INTEGER :: ix, iy, ixy, nxa, nya, nxya

nxa = Core%nxa; nya = Core%nya
nxya = Core%nxya

CX(1) = AsyPitch * Half; CY(1) = -AsyPitch * Half
IF(.NOT. lEdge) THEN
  CX(1) = 0.; CY(1)  = 0;
ENDIF

DO ix = 2, nxa
   CX(ix) = CX(ix-1) + AsyPitch
ENDDO
DO iy = 2, nya
  CY(iy) = CY(iy-1) - AsyPitch
ENDDO

CALL Dmalloc(Core%AsyCentX, nxa)
CALL Dmalloc(Core%AsyCentY, nya)

Core%AsyCentX(1:nxa) = Cx(1:nxa)
Core%AsyCenty(1:nya) = Cy(1:nya)

DO ixy= 1, nxya
  ix = Asy(ixy)%ixa; iy = Asy(ixy)%iya
  Asy(ixy)%Cx = CX(ix); Asy(ixy)%cy = CY(iy)
ENDDO

END SUBROUTINE

SUBROUTINE SetPseudoPinMap()

END SUBROUTINE

SUBROUTINE SetThPinType()
USE PARAM
USE Typedef,     ONLY : CoreInfo_Type,     AsyInfo_type,        Pin_Type
USE GEOM,        ONLY : Pin,               AsyInfo,                          &
                        nxy
IMPLICIT NONE
!TYPE(CoreInfo_Type) :: Core
!TYPE(AsyInfo_type), POINTER :: AsyInfo(:)
!TYPE(Pin_Type), POINTER :: Pin(:)
INTEGER :: ixy, iasytype
!INTEGER :: nxy 

!nxy = Core%nxy
!Pin => Core%Pin;      AsyInfo => Core%AsyInfo

DO ixy = 1, nxy
   IF(Pin(ixy)%lfuel) CYCLE
  iasytype=Pin(ixy)%AsyType
  IF(AsyInfo(iAsyType)%lfuel) THEN
    Pin(ixy)%lGT = .TRUE.
  ELSE
    Pin(ixy)%lRadRef =.TRUE.
  ENDIF
ENDDO

END SUBROUTINE