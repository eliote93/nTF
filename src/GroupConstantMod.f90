MODULE GroupConst_Mod
USE PARAM
USE TYPEDEF,         ONLY : CoreInfo_Type,FmInfo_Type, CMInfo_Type, PinXs_Type , GroupInfo_Type,PE_TYPE
USE BasicOperation, ONLY : CP_VA, CP_CA, MULTI_VA, MULTI_CA
USE files,           ONLY : caseid
USE XSLIB_MOD,  ONLY : enbhel
USE CNTL,         ONLY : nTracerCntl_Type
LOGICAL :: ladf
REAL, POINTER :: cellphibdry(:,:,:,:)
REAL, POINTER :: cellphi1g(:,:,:)
INTEGER :: isosize=200
    

#ifdef gcmodule
INTERFACE
SUBROUTINE GroupConstGen(Core, FmInfo, THInfo, CmInfo, GroupInfo, nTracerCntl, PE, ng)
USE PARAM
USE TYPEDEF,         ONLY : CoreInfo_Type,FmInfo_Type, THInfo_Type, CMInfo_Type, PinXs_Type , GroupInfo_Type,PE_TYPE
USE BasicOperation, ONLY : CP_VA, CP_CA, MULTI_VA, MULTI_CA
USE files,           ONLY : caseid
USE XSLIB_MOD,  ONLY : enbhel
USE CNTL,         ONLY : nTracerCntl_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(CMInfo_Type) :: CMInfo
TYPE(THInfo_Type) :: THInfo
TYPE(FMInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
INTEGER, INTENT(IN) :: ng
INTEGER :: gc_spec
REAL :: Dng(ng)
!REAL :: cellphibdry(:,:,:,:)

END SUBROUTINE

END INTERFACE
#endif
END MODULE