MODULE DcplXsGen_Mod
USE PARAM
IMPLICIT NONE
INTERFACE

SUBROUTINE RefPlnXsGeneration(Core, DcplInfo, THInfo, GroupInfo, DcplCntl, DcplItrCntl, PE, DcplPE, lSubGrp0, lMOC)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type       ,DcplInfo_Type        ,GroupInfo_Type       &
                          ,ThInfo_Type         ,PE_Type               
USE CNTL,           ONLY : nTracerCntl_Type
USE itrcntl_mod,    ONLY : ItrCntl_TYPE                                  
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(DcplInfo_Type) :: DcplInfo
TYPE(ThInfo_Type), POINTER :: THInfo(:)
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: DcplCntl
TYPE(ItrCntl_Type) :: DcplItrCntl(100)
TYPE(PE_TYPE) :: PE
TYPE(PE_TYPE) :: DcplPE(:)
LOGICAL :: lSubGrp0, lMOC

END SUBROUTINE

SUBROUTINE ThXsGeneration(Core, DcplInfo, THInfo, GroupInfo, DcplCntl, DcplItrCntl, PE, DcplPE, lSubGrp0, lMOC)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type       ,DcplInfo_Type        ,GroupInfo_Type       &
                          ,ThInfo_Type         ,PE_Type               
USE CNTL,           ONLY : nTracerCntl_Type
USE itrcntl_mod,    ONLY : ItrCntl_TYPE                                  
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(DcplInfo_Type) :: DcplInfo
TYPE(ThInfo_Type), POINTER :: THInfo(:)
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: DcplCntl
TYPE(ItrCntl_Type) :: DcplItrCntl(100)
TYPE(PE_TYPE) :: PE
TYPE(PE_TYPE) :: DcplPE(:)
LOGICAL :: lSubGrp0, lMOC

END SUBROUTINE


SUBROUTINE DcplPlnMOC(Core, FmInfo, CmInfo, ThInfo, RayInfo, GroupInfo, EigV, ng, DcplCntl, DcplItrCntl, DcplPE, lMOC)
USE PARAM
USE TYPEDEF,         ONLY : DcplInfo_Type     ,CoreInfo_Type     ,RayInfo_Type      &
                           ,FmInfo_Type       ,CmInfo_Type       ,THInfo_Type       &
                           ,GroupInfo_Type    ,PE_TYPE
USE MOC_Mod,         ONLY : MOCSweep
USE CNTL,            ONLY : nTracerCntl_Type
USE ItrCntl_Mod,     ONLY : ItrCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(RayInfo_Type) :: RayInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_TYPE) :: DcplPE
TYPE(ItrCntl_Type) :: DcplItrCntl
TYPE(nTracerCntl_Type) :: DcplCntl
INTEGER :: ng
REAL :: Eigv
LOGICAL :: lMOC
END SUBROUTINE

SUBROUTINE DcplMOCSweep(Core, RayInfo, FMInfo, eigv, ng, PE, GroupInfo, nTracerCntl, ItrCntl)
USE TYPEDEF,     ONLY : CoreInfo_Type      ,RayInfo_Type       ,FmInfo_Type         &
                       ,PE_TYPE            ,GroupInfo_Type
USE CNTL,        ONLY : nTracerCntl_Type
USE itrcntl_mod, ONLY : ItrCntl_TYPE
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(PE_TYPE) :: PE
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(ItrCntl_TYPE) :: ItrCntl
REAL :: eigv
INTEGER :: ng

END SUBROUTINE

SUBROUTINE DcplPlnCmfd(Core, CMInfo, FMInfo, eigv, ng, lcmfd, lreset, GroupInfo, PE, nTracerCntl, ItrCntl)
USE PARAM
USE TYPEDEF,     ONLY : PE_TYPE           ,CoreInfo_Type      ,CMInfo_Type       &
                       ,FmInfo_Type       ,GroupInfo_Type
USE CNTL,        ONLY : nTracerCntl_Type
USE itrcntl_mod, ONLY : ItrCntl_TYPE, CMFDItrCntl_TYPE
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CMInfo_Type) :: CMInfo
TYPE(FMInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_TYPE) :: PE
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(ItrCntl_TYPE) :: ItrCntl
REAL :: eigv
INTEGER :: ng
LOGICAL :: lreset, lcmfd

END SUBROUTINE

SUBROUTINE DcplSetMocAxEff(Core, DcplPxs, phis1g, AxSrc1g, AxPxs1g, iz)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type                                           &
                          ,Pin_Type         ,Cell_Type
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
REAL, POINTER :: AxSrc1g(:), AxPxs1g(:)
REAL :: DcplPxs(:), phis1g(:)
INTEGER :: iz

END SUBROUTINE

SUBROUTINE DcplEffXsSave(Core, Fxr, itemp, GroupInfo, DcplPE)
USE PARAM
USE TYPEDEF,    ONLY : CoreInfo_Type    ,FxrInfo_Type     &
                      ,GroupInfo_Type   ,PE_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_TYpe) :: DcplPE
INTEGER :: itemp

END SUBROUTINE

FUNCTION DcplConvChk(ItrCntl, nRefPln)
USE PARAM
USE TYPEDEF,       ONLY : 
USE ItrCntl_mod,   ONLY : ItrCntl_Type
IMPLICIT NONE
TYPE(ITrCntl_Type) :: ItrCntl(nRefPln)
INTEGER :: nRefPln
LOGICAL :: DcplConvChk

END FUNCTION


END INTERFACE
END MODULE