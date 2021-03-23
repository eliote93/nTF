MODULE DEPL_MOD
USE PARAM
USE TYPEDEF,   ONLY :  XsMac_Type
USE DeplType
USE MatExp_Mod, ONLY : Mat_Type
IMPLICIT NONE

!REAL, POINTER :: VISO(:)

TYPE(DeplCntl_Type) :: DeplCntl
TYPE(DeplVars_Type) :: DeplVars(nThreadMax)
TYPE(XsMac_Type) :: XsMac(nThreadMax)
TYPE(MatExp_Type) :: MatExp(nThreadMax)
REAL, POINTER :: SigFOld(:, :, :)!, SigNFOld(:, :, :)
REAL :: decayXe135,decayI135
REAL, POINTER :: yieldXe135(:),yieldI135(:)
INTERFACE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FxrBurnUp(DeplVars, DeplLib, DeplCntl)
USE PARAM
USE DeplType,    ONLY : DeplVars_Type,      DeplLib_Type,        DeplCntl_Type,          &
                        MatExp_Type
IMPLICIT NONE

TYPE(DeplVars_Type) :: DeplVars
TYPE(DeplLib_Type) :: DeplLib
TYPE(DeplCntl_Type) :: DeplCntl

END SUBROUTINE

SUBROUTINE InitDeplXS(DeplLib)
USE PARAM
USE DeplType,    ONLY : DeplLib_Type
IMPLICIT NONE
TYPE(DeplLib_Type) :: DeplLib

END SUBROUTINE

SUBROUTINE TranDepXS2DeplLib(DeplLib, BurnUpXs, IsoNum, nIsoDep)
USE PARAM
USE DeplType,     ONLY : DeplLib_Type

TYPE(DeplLib_Type) :: DeplLib
REAL :: BurnUpXs(4, nIsoDep)
REAL :: IsoNum(nIsoDep)
INTEGER :: nIsoDep

END SUBROUTINE

SUBROUTINE MakeDeplMat(DMat, DeplLib, PHI, BurnUpTime)
USE PARAM
USE DEPLTYPE,      ONLY : Mat_TYPE,         DeplLib_Type
IMPLICIT NONE
TYPE(Mat_Type) :: DMat
TYPE(DeplLib_Type) :: DeplLib
REAL :: PHI, BurnUpTime

END SUBROUTINE

SUBROUTINE CopyIsoNumVector(V1, V2, n, torr)
IMPLICIT NONE
REAL :: V1(n), V2(n)
INTEGER :: n
REAL :: torr

END SUBROUTINE

SUBROUTINE UpdateDeplFxrInfo(Fxr, DeplVars, GroupInfo, lCorrectStep, lXeDyn)
USE PARAM
USE TYPEDEF,            ONLY : FxrInfo_Type,       GroupInfo_Type
USE DeplType,           ONLY : DeplVars_Type
IMPLICIT NONE
TYPE(FxrInfo_Type) :: Fxr
TYPE(DeplVars_Type) :: DeplVars
TYPE(GroupInfo_Type) :: GroupInfo
LOGICAL :: lCorrectStep, lXeDyn

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE DepletionStep(Core, FmInfo, DeplVars, DeplLib, ThInfo, GroupInfo, nTracerCntl, DeplCntl, PE)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,       FmInfo_Type,       GroupInfo_Type,    &
                             PE_Type,             ThInfo_Type
USE DeplType,         ONLY : DeplLib_Type,        DeplVars_Type,     DeplCntl_Type
USE CNTL,             ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) ::  FmInfo
TYPE(DeplVars_Type) :: DeplVars(nThreadMax)
TYPE(DeplLib_Type) :: DeplLib(nThreadMax)
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(DeplCntl_Type) :: DeplCntl
TYPE(PE_TYPE) :: PE
ENDSUBROUTINE

SUBROUTINE MakeDeplXs1g(Fxr, ifxr, ng, DeplXs, NormFactor, Core, ipin, iz, GroupInfo, PE, DeplCntl)
USE PARAM
USE TypeDef,          ONLY : CoreInfo_Type,     FxrInfo_Type,      Cell_Type,       PE_TYPE, &
                             GroupInfo_Type
USE DeplType,         ONLY : DeplXs_Type, DeplCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type),POINTER :: Fxr(:,:)
TYPE(DeplXs_Type) :: DeplXS
TYPE(PE_TYPE) :: PE
TYPE(DeplCntl_Type) :: DeplCntl
TYPE(GroupInfo_Type) :: GroupInfo
REAL :: NormFactor
INTEGER :: ng, ipin, iz, ifxr


END SUBROUTINE


!SUBROUTINE ConstDeplVas(DeplVars, DeplXs, IdIso, pnum, nIso, nIsoLib, nIsoDepl)
!USE PARAM
!USE DeplType,    ONLY : DeplXs_Type,      DeplVars_Type
!IMPLICIT NONE
!TYPE(DeplVars_Type) ::DeplVars
!TYPE(DeplXs_Type) ::DeplXs
!INTEGER :: IdIso(nIsoLib)
!REAL :: pnum(nIsoLib)
!INTEGER :: nIsoLib, nIsoDepl, nIso
!
!END SUBROUTINE

SUBROUTINE ConstDeplVas(DeplVars, DeplXs, Fxr, nIsoLib, nIsoDepl)
USE PARAM
USE TypeDef,     ONLY : FxrInfo_Type
USE DeplType,    ONLY : DeplXs_Type,      DeplVars_Type
IMPLICIT NONE
TYPE(DeplVars_Type) ::DeplVars
TYPE(DeplXs_Type) ::DeplXs
TYPE(FxrInfo_Type) :: Fxr
INTEGER :: nIsoLib, nIsoDepl
END SUBROUTINE

SUBROUTINE ConstDeplVas_QD(DeplVars, DeplXs, Fxr, nIsoLib, nIsoDepl, DeplCntl)
USE PARAM
USE TypeDef,     ONLY : FxrInfo_Type
USE DeplType,    ONLY : DeplCntl_Type, DeplXs_Type,      DeplVars_Type
IMPLICIT NONE
TYPE(DeplCntl_TYPE) :: DeplCntl
TYPE(DeplVars_Type) ::DeplVars
TYPE(DeplXs_Type) ::DeplXs
TYPE(FxrInfo_Type) :: Fxr
INTEGER :: nIsoLib, nIsoDepl

END SUBROUTINE

SUBROUTINE SaveDeplXs_GD(DeplXs, Fxr, lCorrectStep)
USE PARAM
USE TYPEDEF,      ONLY : FxrInfo_Type
USE DeplType,      ONLY : DeplXs_Type,      DeplVars_Type
USE nuclidmap_mod, ONLY : iposiso
USE BasicOperation, ONLY : CP_CA, CP_VA
IMPLICIT NONE
TYPE(DeplXs_Type) :: DeplXs
TYPE(FxrInfo_Type) :: Fxr
LOGICAL :: lCorrectStep

END SUBROUTINE

SUBROUTINE SaveFxrIsoInfo(Fxr, nFxr)
USE PARAM
USE TYPEDEF,     ONLY : FxrInfo_Type
IMPLICIT NONE
TYPE(FxrInfo_Type) :: Fxr(1:nFxr)
INTEGER :: nFxr

END SUBROUTINE


SUBROUTINE ReadDeplFile(InDev, filename, DeplLib, NISODEP)
USE DeplType,    ONLY : DeplLib_Type
IMPLICIT NONE
TYPE(DeplLib_Type) :: DeplLib
INTEGER :: InDev, NISODEP
character(80) :: filename
END SUBROUTINE

SUBROUTINE Init_Depl(DeplLib, DeplVars, GroupInfo,  PE)
USE PARAM
USE DeplType,       ONLY : DeplLib_Type,      DeplVars_Type,    DeplCntl_Type
USE TYPEDEF,        ONLY : GroupInfo_Type,    PE_TYPE
IMPLICIT NONE
TYPE(DeplLib_Type) :: DeplLib(nThreadMax)
TYPE(DeplVars_Type) :: DeplVars(nThreadMax)
TYPE(GroupInfo_Type) :: GroupInfo

TYPE(PE_TYPE) :: PE
END SUBROUTINE

FUNCTION ElmtLocFindnUpdt(DMAT, I, IMAT)
USE MatExp_mod,  ONLY :Mat_TYPE
TYPE(Mat_Type) :: DMAT
INTEGER :: I, IMAT
INTEGER :: ElmtLocFindnUpdt

END FUNCTION

SUBROUTINE AllocDeplFXRMem(Core, Fxr, GroupInfo, PE)
USE PARAM
USE TYPEDEF,    ONLY : CoreInfo_Type,     FxrInfo_Type,         PE_TYPE,          &
                       GroupInfo_Type,    Pin_Type,             Cell_Type
USE ALLOCS
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: FXR(:, :)
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_Type) :: PE

END SUBROUTINE

!FUNCTION FxrAvgPhi(Core, Fxr, Phis, ipin, iLocalfxr, iz, ng, PE)
!USE PARAM
!USE TYPEDEF,     ONLY : CoreInfo_Type,    PE_Type,     FxrInfo_Type,               &
!                        Cell_Type,    Pin_Type
!IMPLICIT NONE
!TYPE(CoreInfo_Type) :: Core
!TYPE(FxrInfo_Type) :: Fxr(:, :)
!TYPE(PE_Type) :: PE
!REAL :: phis(:, :, :)
!INTEGER :: iLocalfxr, ipin, iz, ng
!REAL :: FxrAvgPhi(ng)
!
!END FUNCTION

SUBROUTINE SetBurnUpStepInfo(PowerCore, DeplCntl)
USE PARAM
USE DeplType,   ONLY : DeplCntl_Type
IMPLICIT NONE
TYPE(DeplCntl_Type) :: DeplCntl
REAL :: PowerCore

END SUBROUTINE

SUBROUTINE SetLocalBurnup(Core, FXR, Power, normalizer, Tsec, lCorrectStep, PE)
USE PARAM
USE TYPEDEF,       ONLY : CoreInfo_Type,   FxrInfo_Type, PE_Type,      &
                         Pin_Type,        Cell_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
REAL, POINTER :: Power(:, :)
REAL :: Normalizer, Tsec
LOGICAL :: lCorrectStep
TYPE(PE_TYPE) :: PE

END SUBROUTINE


FUNCTION GetCoreHmMass(Core, FXR, PE)  !Calculate the 
USE PARAM 
USE TYPEDEF,       ONLY : CoreInfo_Type,   FxrInfo_TYPE,   PE_TYPE
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :) 
TYPE(PE_TYPE) :: PE

REAL :: GetCoreHmMass

END FUNCTION



FUNCTION FluxNormalizeFactor(Core, FmInfo, GroupInfo, PowerCore, lCritSpec, lXsLib, PE)
USE PARAM
USE TYPEDEF,    ONLY : CoreInfo_Type,        FmInfo_Type,         GroupInfo_TYPE, &
                       PE_TYPE
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_TYPE) :: PE

LOGICAL :: lCritSpec, lXsLib
REAL :: FluxNormalizeFactor, PowerCore

END FUNCTION

SUBROUTINE EqXeUpdate(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,       FmInfo_Type,       GroupInfo_Type,    &
                             THInfo_Type,         PE_Type
USE CNTL,             ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE

END SUBROUTINE

SUBROUTINE SetDeplCoreState(ibu, CoreState, ThInfo, ThVar, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,  ONLY : ThInfo_Type,      ThVar_Type, PE_TYPE
USE DeplType, ONLY : CoreState_Type
USE CNTL,     ONLY : nTracerCntl_Type
IMPLICIT NONE

TYPE(CoreState_Type) :: CoreState
TYPE(THInfo_Type) :: ThInfo
TYPE(THVar_Type) :: ThVar
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
INTEGER :: ibu
  
END SUBROUTINE 

SUBROUTINE SaveCoreState(ibu, CoreState, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,  ONLY : ThInfo_Type,      ThVar_Type, PE_TYPE
USE DeplType, ONLY : CoreState_Type
USE CNTL,     ONLY : nTracerCntl_Type
IMPLICIT NONE

TYPE(CoreState_Type) :: CoreState
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
INTEGER :: ibu
  
END SUBROUTINE 

SUBROUTINE UpdtFluxSolution(ibu, DeplCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,         ONLY : PE_Type
USE DeplType,        ONLY : DeplCntl_Type
USE CNTL,            ONLY : nTracerCntl_Type
IMPLICIT NONE
INTEGER :: ibu
TYPE(DeplCntl_Type) :: DeplCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

END SUBROUTINE



SUBROUTINE SetDeplTimeStep(ibu, Core, DeplCntl, nTracerCntl, PE)
USE TYPEDEF,        ONLY : CoreInfo_Type,     PE_Type
USE DeplType,       ONLY : DeplCntl_Type
USE CNTL,           ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(DeplCntl_Type) :: DeplCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
INTEGER :: ibu

END SUBROUTINE 
SUBROUTINE EffBurnUpTime(Tsec, PowLv0, PowLv1)
USE PARAM
IMPLICIT NONE
REAL :: Tsec, PowLv0, PowLv1

END SUBROUTINE



SUBROUTINE DeplBanner(imode)
IMPLICIT NONE
INTEGER :: imode

END SUBROUTINE

END INTERFACE


END MODULE