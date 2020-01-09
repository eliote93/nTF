module geom
use param
use typedef, only : AsyInfo_type,   pininfo_type,   coreinfo_type,    &
                    cell_type,      pin_type,       Asy_Type,         &
                    AsyGap_type,    MiscStruct_Type,  BasicGeom
implicit none
save
INTEGER :: ng, ngthrm !, ngfg, norg, ntiso
INTEGER :: nPreC = 6                    !Number if precursor
INTEGER :: nxy,nx,ny,nz
INTEGER :: nbd = 4, ncbd = 4
INTEGER :: nCellX0, nCellX
INTEGER :: lCoreAng 
LOGICAL :: lgap = .false.
LOGICAL :: lhgap = .false.  !--- 16/05/23 BYS edit : Homogenized Gap cell option
LOGICAL :: lEdge = .TRUE.
LOGICAL :: lROT = .FALSE.
LOGICAL :: lCbd = .FALSE.               !Checkerboard

INTEGER :: nzfm
INTEGER :: ndir
REAL :: albeDO(10)
!INTEGER :: myzb,myze,myzbf,myzef
!
REAL :: CellPitch,AsyPitch,CoreHight
REAL,pointer :: Hz(:),Hzfm(:), HzInv(:), HzfmInv(:)
!
TYPE(coreinfo_type) :: core
INTEGER :: nAsy
TYPE(Pin_type), POINTER:: Pin(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(MiscStruct_Type), TARGET :: MiscStruct
!INTEGER,pointer :: CoreIdx(:,:) 
!INTEGER,pointer :: CoreMap(:)
!SubPlane Variables
INTEGER :: nSubPlane = 1
INTEGER, POINTER :: SubPlaneMap(:)
INTEGER, POINTER :: SubPlaneRange(:, :)

!
INTEGER :: nAsyType, nAsyType0
type(AsyInfo_type),pointer :: AsyInfo(:)

INTEGER :: nPinType, nPinType0, nGappinType
INTEGER, POINTER :: PinTypeMap(:)   !--- CNJ Edit : Flexible Cell & Pin Numbering
type(PinInfo_type),pointer :: PinInfo(:)

INTEGER :: nCellType, nCellType0, nGapType, nVssTyp
INTEGER, POINTER :: CellTypeMap(:)   !--- CNJ Edit : Flexible Cell & Pin Numbering
type(cell_type), pointer :: CellInfo(:), BaseCellInfo(:)  !--- JSR Edit : nTIG Restart

!--- 180625 JSR
INTEGER :: nBaseCell, nBasecell0

TYPE(BasicGeom) :: AsyGeom(0:3) 

INTEGER :: nAsyGapType
TYPE(AsyGap_Type), POINTER :: AsyGap(:)

REAL, POINTER :: PinVol(:,:), PinVolFM(:, :)      !PinCell Volume
REAL, POINTER :: AsyVol(:,:)



END module