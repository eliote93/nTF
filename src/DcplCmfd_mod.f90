Module DcplCMFD_mod
USE PARAM
USE TYPEDEF, ONLY : FxrInfo_Type, Cell_Type
IMPLICIT NONE


INTERFACE
SUBROUTINE AllocDcplCellXs(CellXsData, nRefTemp, nFsrMax, ng)
USE PARAM
USE Typedef,    ONLY : CellXsData_Type
!USE DcplCmfd_Mod,  ONLY : CellXSData_Type
IMPLICIT NONE
TYPE(CellXSData_Type) :: CellXsData
INTEGER :: nRefTemp, nFsrMax, ng

END SUBROUTINE

SUBROUTINE DcplHomXsGen(Core, DcplInfo, Fxr, PinXs, myzb, myze, ng, lXsGen, lXsLib, lXsFtn, XsFtnMod, GroupInfo)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type    ,DcplInfo_Type   ,PinXs_Type      &
                          ,FxrInfo_Type                                       &
                          ,GroupInfo_Type
USE BasicOperation, ONLY : CP_VA 
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(DcplInfo_Type) :: DcplInfo
TYPE(PinXs_Type), POINTER :: PinXS(:, :)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
INTEGER :: myzb, myze, ng, XsFtnMod
LOGICAL :: lXsLib, lXsFtn, lXsGen
TYPE(GroupInfo_Type) :: GroupInfo

END SUBROUTINE

!SUBROUTINE DcplCellXsGen(PinXs, CellXsData, iRefPln, iRefPln0, ng, lfeedback, lxslib, lFuelPln, GroupInfo)
SUBROUTINE DcplCellXsGen(PinXs, CellXsData, ng, lfeedback, lxslib, lFuelPln, GroupInfo)
USE PARAM
USE TYPEDEF,      ONLY : PinXs_Type      ,CellXsData_Type       ,GroupInfo_Type
IMPLICIT NONE
TYPE(PinXs_Type) :: PinXs
TYPE(CellXsData_Type) :: CellXsData
TYPE(GroupInfo_Type) :: GroupInfo
INTEGER :: ng
!INTEGER :: iRefPln0, iRefPln, ng
LOGICAL :: lFeedBack, lXsLib, lFuelPln

END SUBROUTINE


SUBROUTINE DcplCellXsGenXsFtn(PinXS, Fxr, CellInfo, phis, ng, nFsr, nFxr, lxslib, lfuelpln, GroupInfo)
USE PARAM
USE TYPEDEF,        ONLY : PinXs_Type      ,FxrInfo_Type        ,GroupInfo_Type       &
                          ,Cell_Type
IMPLICIT NONE

TYPE(PinXS_TYPE) :: PinXS(3)
TYPE(FxrInfo_Type) :: Fxr(nFxr)
TYPE(Cell_Type), TARGET :: CellInfo
TYPE(GroupInfo_Type) :: GroupInfo

REAL :: phis(nfsr, ng, 3)
INTEGER :: nFsr, nFxr, ng
LOGICAL :: lXsLib, lFuelPln

END SUBROUTINE

SUBROUTINE DcplRadCouplingCoeffGen(Pin, DcplInfo, lXsFtn, lXsGen)
USE PARAM
USE TYPEDEF,      ONLY : PinXs_TYPE     ,Pin_Type        ,DcplInfo_Type
IMPLICIT NONE
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(DcplInfo_Type) :: DcplInfo
LOGICAL :: lXsFtn, lXsGen

END SUBROUTINE


SUBROUTINE DcplAxSrcUpdt(Core, CmInfo, DcplInfo, ng, PE, lFeedBackInt)
USE PARAM
USE TYPEDEF, ONLY : CoreINfo_Type   ,CmInfo_Type      ,PE_TYPE      &
                   ,DcplInfo_Type   ,PE_TYPE
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(DcplInfo_Type) :: DcplInfo
TYPE(PE_TYPE) :: PE
INTEGER :: ng
LOGICAL :: lFeedBackInt

END SUBROUTINE

END INTERFACE
END Module