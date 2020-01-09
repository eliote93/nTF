MODULE RAYS
USE PARAM  
USE TYPEDEF,   ONLY : AziAngleInfo_Type,    PolarAngle_Type,          &
                      ModRayInfo_type,      AsyRayInfo_type,          &
                      CoreRayInfo_type,     RayInfo_Type,             &
                      CellRayInfo_type,     BasicGeom,                &
                      CellRayBase_Type,     RotRayInfo_Type,          &
                      RayInfo4Cmfd_Type,    FastCoreRayDat_TYpe,      &
                      !--- CNJ Edit : Domain Decomposition
                      DcmpAsyRayInfo_Type

IMPLICIT NONE
SAVE
!Azimuthal Angle : Azimuthal Angle
TYPE(AziANgleInfo_Type), POINTER :: AziAngle(:)
!Polar Angle : Polar Angle
TYPE(PolarAngle_Type), POINTER :: PolarAngle(:)
!Modular Ray : Modular Ray
TYPE(ModRayInfo_Type), POINTER :: ModRay(:)
!Assembly Ray : Assembly Boundary to Assembly Boundary
TYPE(AsyRayInfo_Type), POINTER :: AsyRay(:)
!Core Ray : Core Boundary to Core Boundary 
TYPE(CoreRayInfo_Type), POINTER :: CoreRay(:)

TYPE(RotRayInfo_Type), POINTER :: RotRay(:)

!--- CNJ Edit : Domain Decomposition
TYPE(DcmpAsyRayInfo_Type), POINTER :: DcmpAsyRay(:, :)

TYPE(CellRayInfo_Type), POINTER :: CellRay1D

TYPE(FastCoreRayDat_Type), POINTER :: FastCoreRayDat(:, :)
TYPE(RayInfo4Cmfd_Type), POINTER :: RayInfo4Cmfd
TYPE(RayInfo_Type) :: RayInfo
!Cycle Ray : 
!TYPE(CycleRayInfo_Type), POINTER :: CycleRay(:)

!INTEGER, POINTER :: nAziAngle, nPolarAngle
!INTEGER :: nPolarAngleHemi
!Assembly Geometry Information
TYPE(CellRayBase_Type), POINTER :: CellRayBase(:)

END MODULE