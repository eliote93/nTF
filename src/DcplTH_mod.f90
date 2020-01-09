MODULE DcplTh_Mod
IMPLICIT NONE

INTERFACE

SUBROUTINE CopyThCondition(Core, Fxr, ThInfo, DcplFxr, DcplThInfo, iPln, ItrCntl)
USE PARAM
USE TypeDef,      ONLY : CoreInfo_Type      ,DcplInfo_Type       ,FxrInfo_Type     &
                        ,Pin_Type           ,Cell_Type           ,ThInfo_Type
USE itrcntl_mod,  ONLY : ItrCntl_TYPE
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:), DcplFxr(:)
TYPE(ThInfo_Type) :: ThInfo
TYPE(ThInfo_Type) :: DcplThInfo
INTEGER :: iPln
TYPE(ItrCntl_Type) :: ItrCntl
END SUBROUTINE

SUBROUTINE CopyThCondition_XsFtn(Core, Fxr, ThInfo, DcplFxr, DcplThInfo, iRefPln, iPln1, ipln2, ItrCntl)
USE PARAM
USE TypeDef,      ONLY : CoreInfo_Type      ,DcplInfo_Type       ,FxrInfo_Type     &
                        ,ThInfo_Type        ,Pin_Type            ,Cell_Type
USE itrcntl_mod,  ONLY : ItrCntl_TYPE
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(FxrInfo_Type) :: DcplFxr(:)
TYPE(ThInfo_Type) :: ThInfo
TYPE(ThInfo_Type) :: DcplThInfo
INTEGER :: iRefPln, iPln1, iPln2
TYPE(ItrCntl_Type) :: ItrCntl
LOGICAL :: lThConv

END SUBROUTINE 

SUBROUTINE SetXsGenThCondition(Core, FmInfo, DcplInfo, THInfo, DcplThInfo, lXsFtn, DcplItrCntl, PE)
USE PARAM
USE TypeDef,      ONLY : CoreInfo_Type     ,FmInfo_Type     ,DcplInfo_Type        &
                        ,ThInfo_Type       ,PE_TYPE
USE itrcntl_mod,  ONLY : ItrCntl_TYPE
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(DcplInfo_Type) :: DcplInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(ThInfo_Type), POINTER :: DcplThInfo(:)
TYPE(ItrCntl_Type) :: DcplItrCntl(100)
TYPE(PE_TYPE) :: PE
LOGICAL :: lXsFtn
END SUBROUTINE 

SUBROUTINE DcplThUpdate(Core, CmInfo, FmInfo, ThInfo, DcplInfo, DcplThInfo, Eigv, ng, lUpdt, GroupInfo, nTracerCntl, DcplItrCntl, PE)
USE PARAM
USE TypeDef,          ONLY : CoreInfo_Type      ,CmInfo_Type      ,FmInfo_Type            &
                            ,ThInfo_Type        ,DcplInfo_Type    ,PE_Type                &
                            ,GroupInfo_Type
USE CNTL,             ONLY : nTracerCntl_Type
USE itrcntl_mod,      ONLY : ItrCntl_TYPE  
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CMInfo_Type) :: CmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(DcplInfo_Type) :: DcplInfo
TYPE(ThInfo_Type), POINTER :: DcplThInfo(:)
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_TYPE) :: PE
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(ItrCntl_Type) :: DcplItrCntl(100)
REAL :: Eigv
INTEGER :: ng
LOGICAL :: lUpdt
END SUBROUTINE

SUBROUTINE DcplModTChg(Core, Fxr, ThInfo, iRefPln, ipln, imod, frac)
USE PARAM
USE TYPEDEF,     ONLY : CoreInfo_Type      ,FxrInfo_Type      ,ThInfo_Type       &
                       ,Pin_Type           ,Cell_Type
USE SteamTBL_mod, ONLY : steamtbl
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(ThInfo_Type) :: ThInfo


INTEGER :: iRefPln, ipln, imod
REAL :: Frac

END SUBROUTINE

SUBROUTINE DcplFuelTChg(Core, Fxr, ThInfo, iRefPln, ipln, frac)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type     ,FxrInfo_Type     ,ThInfo_Type             &
                            ,Pin_Type          ,Cell_Type
USE TH_Mod,           ONLY : ThOpt             ,ThVar
USE FuelProperty_Mod, ONLY : fhtcoef
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(ThInfo_Type) :: ThInfo
INTEGER :: ipln, iRefPln

REAL :: Frac

END SUBROUTINE


SUBROUTINE SetThXsGenCondition(Core, Fxr, ThInfo, iRefPln, ipln, imod, frac, boronppm)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type        ,FxrInfo_Type        ,THInfo_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(ThInfo_Type) :: ThInfo
INTEGER :: ipln, iRefPln, imod
REAL :: frac
REAL :: boronppm
END SUBROUTINE

SUBROUTINE SaveXsGenInfo(Core, FmInfo, CmInfo, imode, DcplPE)
USE PARAM
USE TYPEDEF,     ONLY : CoreInfo_Type         ,FmInfo_Type         ,PE_Type       &
                       ,FxrInfo_Type          ,Pin_Type            ,Cell_Type     &
                       ,CmInfo_Type
USE BasicOperation, ONLY : CP_VA
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(PE_Type) :: DcplPE
INTEGER :: imode

END SUBROUTINE
END INTERFACE

END MODULE