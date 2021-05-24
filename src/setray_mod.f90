MODULE SetRay

INTERFACE

SUBROUTINE SetModularRayAngle(AziAngleInfo, RayInfo, AsyPitch)
USE TYPEDEF, ONLY : AziAngleInfo_Type,  PolarAngle_Type,  RayInfo_Type
TYPE(AziAngleInfo_Type), POINTER :: AziANgleInfo(:)
TYPE(RayInfo_Type) :: RayInfo
REAL :: AsyPitch
END SUBROUTINE

SUBROUTINE SetPolarRayAngle(PolarAngleInfo, RayInfo)
USE TYPEDEF, ONLY : PolarAngle_Type,  RayInfo_Type
TYPE(PolarAngle_Type), POINTER :: PolarAngleInfo(:)
TYPE(RayInfo_Type) :: RayInfo
ENDSUBROUTINE

SUBROUTINE SetModRay(ModRay, RayInfo, AsyGeom, AsyPitch, lEdge)
USE TYPEDEF, ONLY : RayInfo_Type, ModRayInfo_Type, BasicGeom
IMPLICIT NONE
TYPE(RayInfo_Type) :: RayInfo
TYPE(ModRayInfo_Type), POINTER :: ModRay(:)
TYPE(BasicGeom) :: AsyGeom(0:3)
REAL :: AsyPitch
LOGICAL :: lEdge
END SUBROUTINE

SUBROUTINE SetAsyRay(AsyRay, RayInfo, AsyGeom, lEdge)
USE TYPEDEF, ONLY : RayInfo_Type, AsyRayInfo_type, BasicGeom
IMPLICIT NONE
TYPE(RayInfo_Type) :: RayInfo
TYPE(AsyRayInfo_type), POINTER :: AsyRay(:)
TYPE(BasicGeom) :: AsyGeom(0:3)
LOGICAL :: lEdge

END SUBROUTINE

SUBROUTINE SetModRayLinking( ModRay, IdxBeg, IdxEnd)
USE TYPEDEF, ONLY : ModRayInfo_Type
IMPLICIT NONE
TYPE(ModRayInfo_Type), POINTER :: ModRay(:)
INTEGER :: IdxBeg, IdxEnd

END SUBROUTINE

SUBROUTINE CellRayGen(RayInfo, CellRayBase, AsyGeom, lEdge)
USE PARAM
USE typedef, ONLY : RayInfo_Type,      &
                    CellRayBase_Type,       &
                    BasicGeom
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CellRayBase_Type), POINTER :: CellRayBase(:)
Type(BasicGeom) :: AsyGeom(0:3)
LOGICAL :: lEdge

END SUBROUTINE

SUBROUTINE Copy_CellRayInfo(CellRay0, CellRay)
USE typedef, only : CellRayInfo_type
IMPLICIT NONE
TYPE(CellRayInfo_type), intent(IN) :: CellRay0
TYPE(CellRayInfo_type), intent(OUT) :: CellRay
END SUBROUTINE

SUBROUTINE VolCorChk(AsyGeom, RayInfo, CellInfo, PinUsed, nCellType, iType)
USE PARAM
USE TYPEDEF, ONLY : BasicGeom,          Cell_type,           RayInfo_Type

IMPLICIT NONE
TYPE(BasicGeom) :: AsyGeom
TYPE(RayInfo_type) :: RayInfo
TYPE(Cell_type), POINTER :: CellInfo(:)
LOGICAL :: PinUsed(:, :)
INTEGER :: nCellType
INTEGER :: itype
END SUBROUTINE

SUBROUTINE BaseCellRayVolCor(RayInfo, CellInfo, CellRayBase, nDiv)
!Make sets of Cellrays which 
USE typedef, ONLY : Cell_type,  CellRayBase_Type,  RayInfo_Type
!USE ALLOCS
IMPLICIT NONE
TYPE(RayInfo_type) :: RayInfo
TYPE(Cell_type) :: CellInfo
TYPE(CellRayBase_Type) :: CellRayBase
INTEGER :: nDiv

END SUBROUTINE

SUBROUTINE CoreRayGen(RayInfo, CoreRay, Core, Asy, AsyGeom, lEdge)
USE TYPEDEF,     ONLY : RayInfo_type,     CoreRayInfo_type,                       &
                        CoreInfo_type,     Asy_Type,           BasicGeom
IMPLICIT NONE                        
TYPE(RayInfo_type) :: RayInfo
TYPE(CoreRayInfo_type), POINTER :: CoreRay(:)       !Core Ray variable which will be 
TYPE(CoreInfo_type) :: Core                          !Croe Information of the problem
TYPE(Asy_Type), POINTER :: Asy(:)                             !Information of Assemblies which resides inside of core
TYPE(BasicGeom)  :: AsyGeom(0:3)                     !Assembly Geometries Information
Logical :: lEdge
END SUBROUTINE

SUBROUTINE RotationRayGen(RayInfo, RotRay, Core, Asy, lEdge)
USE TYPEDEF,        ONLY : RayInfo_type,   CoreInfo_type,   RotRayInfo_Type, &
                           Asy_Type
IMPLICIT NONE

TYPE(RayInfo_type) :: RayInfo
TYPE(RotRayInfo_Type), POINTER :: RotRay(:)  
TYPE(CoreInfo_type) :: Core          
TYPE(Asy_Type), POINTER :: Asy(:)
LOGICAL :: lEdge
END SUBROUTINE

!--- CNJ Edit : Domain Decomposition
SUBROUTINE DcmpRayGen(Core, RayInfo, DcmpAsyRay)
USE TYPEDEF,        ONLY : RayInfo_type,    CoreInfo_type,    DcmpAsyRayInfo_Type
IMPLICIT NONE
TYPE(CoreInfo_type) :: Core
TYPE(RayInfo_type) :: RayInfo
TYPE(DcmpAsyRayInfo_type), POINTER  :: DcmpAsyRay(:, :)
END SUBROUTINE

SUBROUTINE HexDcmpRayGen(Core, RayInfo, DcmpAsyRay)

USE TYPEDEF, ONLY : RayInfo_type, CoreInfo_type, DcmpAsyRayInfo_Type

IMPLICIT NONE

TYPE(CoreInfo_type) :: Core
TYPE(RayInfo_type) :: RayInfo
TYPE(DcmpAsyRayInfo_type), POINTER, DIMENSION(:,:)  :: DcmpAsyRay

END SUBROUTINE

SUBROUTINE RayInfo4CmfdGen(RayInfo, RayInfo4Cmfd, Core)
USE PARAM
USE TYPEDEF,        ONLY : RayInfo_type,      RayInfo4Cmfd_Type,  CoreInfo_type
USE ALLOCS                           
IMPLICIT NONE

TYPE(CoreInfo_type) :: Core
TYPE(RayInfo_type) :: RayInfo
TYPE(RayInfo4Cmfd_type), POINTER :: RayInfo4Cmfd

END SUBROUTINE

!--- CNJ Edit
SUBROUTINE RayInfoMaxSize(Core, RayInfo, myzb, myze)
USE PARAM
USE TYPEDEF,        ONLY : RayInfo_type, CoreInfo_type
IMPLICIT NONE

TYPE(CoreInfo_type) :: Core
TYPE(RayInfo_type) :: RayInfo
INTEGER :: myzb, myze

END SUBROUTINE

END INTERFACE



END MODULE