MODULE TRAN_MOD
USE PARAM
USE TYPEDEF,    ONLY : TranInfo_TYPE,       XsChange_TYPE,     TranCntl_TYPE
IMPLICIT NONE


TYPE(TranCntl_Type), TARGET :: TranCntl
TYPE(XsChange_Type), TARGET :: XsChange(100)

TYPE(TranInfo_Type) :: TranInfo



INTERFACE


SUBROUTINE TransientFsp_Driver(Core, FmInfo, CmInfo, TranInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,           ONLY : CoreInfo_Type,           FmInfo_Type,           CmInfo_Type,         &
                              TranInfo_Type,           GroupInfo_Type,        TranCntl_Type,       &
                              ThInfo_Type,             PE_Type
USE CNTL,              ONLY : nTracerCntl_Type

IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
END SUBROUTINE



SUBROUTINE InitTransient(Core, RayInfo, FmInfo, CmInfo, ThInfo, eigv, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,     FmInfo_Type,      CmInfo_Type,      TranInfo_Type,     &
                           RayInfo_Type,      GroupInfo_Type,   TranCntl_Type,    PE_TYPE,           &
                           ThInfo_Type
USE CNTL,           ONLY : nTracerCntl_Type
IMPLICIT NONE


TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE

REAL :: eigv
END SUBROUTINE

SUBROUTINE FluxNormTransient(Core, RayInfo, FmInfo, CmInfo, TranInfo, ng, nTracerCntl, PE)
!Flux level normalization
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,     FmInfo_Type,       CmInfo_Type,     &
                           RayInfo_Type,      TranInfo_Type,     PE_TYPE
USE CNTL,           ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
INTEGER :: ng

END SUBROUTINE

SUBROUTINE InitPrecursor(Core, FmInfo, CmInfo, TranInfo, GroupInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,     ONLY : CoreInfo_Type,       FmInfo_Type,           CmInfo_Type,        &
                        TranInfo_Type,       PE_Type,               GroupInfo_Type
USE CNTL,        ONLY : nTracerCntl_Type
IMPLICIT NONE
  
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

END SUBROUTINE

SUBROUTINE SetTimeStep(TranCntl)
USE PARAM
USE TYPEDEF,       ONLY : TranCntl_Type
IMPLICIT NONE

TYPE(TranCntl_Type) :: TranCntl

END SUBROUTINE

SUBROUTINE UpdtPowerLevel(Core, FmInfo, CmInfo, TranInfo, ng, nTracerCntl, PE)
! Flux level normalization
! normalize flux such that average flux in fuel region be unity
! then update fission source and moments accordingly
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,     FmInfo_Type,       CmInfo_Type,     &
                           TranInfo_Type,     PE_TYPE,                            &
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
TYPE(TranInfo_Type) :: TranInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
INTEGER :: ng

END SUBROUTINE

SUBROUTINE UpdtPrec(Core, FmInfo, CmInfo, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,             ONLY : CoreInfo_Type,          FmInfo_Type,           CmInfo_Type,         &
                                TranInfo_Type,          GroupInfo_Type,        TranCntl_Type,       &
                                PE_Type,                                                            &
                                FxrInfo_Type,           Pin_Type,               Cell_Type
USE CNTL,                ONLY : nTracerCntl_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

END SUBROUTINE

SUBROUTINE CmFmInfoSync(Core, FmInfo, CmInfo, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,             ONLY : CoreInfo_Type,          FmInfo_Type,           CmInfo_Type,         &
                                TranInfo_Type,          GroupInfo_Type,        TranCntl_Type,       &
                                PE_Type,                                                            &
                                FxrInfo_Type,           Pin_Type,               Cell_Type
USE CNTL,                ONLY : nTracerCntl_Type
USE MOC_MOD,             ONLY : PsiUpdate
USE BasicOperation,      ONLY : MULTI_CA
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

END SUBROUTINE

SUBROUTINE SaveTranSol(Core, FmInfo, CmInfo, TranInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,             ONLY : CoreInfo_Type,          FmInfo_Type,        CmInfo_Type,       &
                                TranInfo_Type,          TranCntl_Type,      PE_Type,           &
                                ThInfo_Type,            GroupInfo_Type
USE CNTL,                ONLY : nTracerCntl_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(THInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

END SUBROUTINE

SUBROUTINE UpdtResSrc(Core, FmInfo, CmInfo, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,             ONLY : CoreInfo_Type,          FmInfo_Type,         CmInfo_Type,       &
                                TranInfo_Type,          GroupInfo_Type,      TranCntl_Type,     &
                                PE_Type
USE CNTL,                ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCNtl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE

END SUBROUTINE
!
!SUBROUTINE XsPerturbation(TranInfo, TranCntl, nTracerCntl)
!USE PARAM
!USE TYPEDEF,             ONLY : TranInfo_Type,          TranCntl_Type
!USE CNTL,                ONLY : nTracerCntl_Type
!IMPLICIT NONE
!TYPE(TranInfo_Type) :: TranInfo
!TYPE(TranCntl_Type) :: TranCntl
!TYPE(nTracerCntl_Type) :: nTracerCntl
!END SUBROUTINE


SUBROUTINE KinParamGen(Core, FmInfo, TranInfo, ThInfo,  GroupInfo,  lBetaUpdt, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,    FmInfo_Type,      ThInfo_Type,          & 
                             TranInfo_Type,    GroupInfo_Type,   PE_Type,              &
                             FXRInfo_Type,     PinXs_Type,       Pin_Type,             &              
                             PinInfo_Type,     Cell_Type
USE CNTL,             ONLY : nTracerCntl_Type
USE BenchXs,          ONLY : XsBaseBen,        DnpBetaBen,        NeutVeloBen
USE BasicOperation,   ONLY : CP_VA,            CP_CA,            MULTI_VA
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) ::  GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
LOGICAL :: lBetaUpdt
!
END SUBROUTINE


SUBROUTINE SolExpExtpltn(Core, FmInfo, CmInfo, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,                 ONLY : CoreInfo_Type,          FmInfo_Type,            CmInfo_Type,       &
                                    TranInfo_Type,          TranCntl_Type,          PE_Type,           &
                                    GroupInfo_Type
USE CNTL,                    ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

END SUBROUTINE

SUBROUTINE UpdtExpTrsf(Core, FmInfo, CmInfo, TranInfo, GroupInfo, lupdt, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,                 ONLY : CoreInfo_Type,          FmInfo_Type,            CmInfo_Type,       &
                                    TranInfo_Type,          TranCntl_Type,          PE_Type,           &
                                    GroupInfo_Type
USE CNTL,                    ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
LOGICAL :: lupdt

ENDSUBROUTINE

SUBROUTINE CheckCondMOC(Core, FmInfo, CmInfo, TranInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,            ONLY : CoreInfo_Type,         FmInfo_Type,         CmInfo_Type,    &
                               TranInfo_Type,         GroupInfo_Type,      TranCntl_Type,  &
                               ThInfo_Type,           PE_Type
USE CNTL,               ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(THInfo_Type) :: THInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

END SUBROUTINE

SUBROUTINE UpdtBaseXsCondiMOC(Core, FmInfo, CmInfo, TranInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,            ONLY : CoreInfo_Type,         FmInfo_Type,         CmInfo_Type,    &
                               TranInfo_Type,         GroupInfo_Type,      TranCntl_Type,  &
                               ThInfo_Type,           PE_Type
USE CNTL,               ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(THInfo_Type) :: THInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

END SUBROUTINE

FUNCTION TranReactivityUpdt(CmInfo, eigv, TranCntl, PE)
USE PARAM
USE TYPEDEF,   ONLY : CmInfo_Type,  PE_TYPE, TranCntl_Type
IMPLICIT NONE

TYPE(CmInfo_Type) :: CmInfo
REAL :: eigv
TYPE(PE_TYPE) :: PE
TYPE(TranCntl_Type) :: TranCntl
REAL :: TranReactivityUpdt
END FUNCTION

FUNCTION UpdtHomXsTr1g(PinXS, PhiC, ng)
USE PARAM
USE TYPEDEF,       ONLY : PinXS_Type
IMPLICIT NONE
TYPE(PinXS_Type) :: PinXS
REAL :: PhiC(ng)
REAL :: UpdtHomXsTr1g
INTEGER :: ng

END FUNCTION

FUNCTION UpdtSigA1g(Fxr, PhiFxr, ng, GroupInfo, nTracerCntl)
USE PARAM
USE TYPEDEF,              ONLY : FxrInfo_Type,         ThInfo_Type,        GroupInfo_Type
USE CNTL,                 ONLY : nTracerCntl_Type
IMPLICIT NONE

REAL :: UpdtSigA1g
TYPE(FxrInfo_Type) :: Fxr
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
REAL :: PhiFxr(ng)
INTEGER :: ng


END FUNCTION


FUNCTION UpdtThermalSigA(Fxr, PhiFxr, ng, GroupInfo, nTracerCntl)
USE PARAM
USE TYPEDEF,              ONLY : FxrInfo_Type,         ThInfo_Type,        GroupInfo_Type,       &
                                 XsMAc_Type
USE CNTL,                 ONLY : nTracerCntl_Type
USE MacXsLib_mod,         ONLY : MacXsBase
USE XsUtil_mod,           ONLY : GetXsMacDat,          ReturnXsMacDat
IMPLICIT NONE

REAL :: UpdtThermalSigA

TYPE(FxrInfo_Type) :: Fxr
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
REAL :: PhiFxr(ng)
INTEGER :: ng

END FUNCTION

FUNCTION UpdtSigA4g(Fxr, PhiFxr, ng, GroupInfo, nTracerCntl)
USE PARAM
USE TYPEDEF,              ONLY : FxrInfo_Type,         ThInfo_Type,        GroupInfo_Type,       &
                                 XsMAc_Type
USE CNTL,                 ONLY : nTracerCntl_Type
USE MacXsLib_mod,         ONLY : MacXsBase
USE XsUtil_mod,           ONLY : GetXsMacDat,          ReturnXsMacDat
IMPLICIT NONE

REAL :: UpdtSigA4g(4)

TYPE(FxrInfo_Type) :: Fxr
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
REAL :: PhiFxr(ng)
INTEGER :: ng

END FUNCTION

END INTERFACE

END MODULE
