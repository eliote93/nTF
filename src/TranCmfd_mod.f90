MODULE TRANCMFD_MOD
USE PARAM
IMPLICIT NONE

INTEGER :: nprec 

REAL, POINTER :: TrSrc(:, :), PrecSrc(:, :)
REAL, POINTER :: Chid(:), lambda(:) !, neut_velo(:)
REAL, POINTER :: CellOmegam(:, :, :), Cellomega0(:, :, :), Cellomegap(:, :, :)
REAL, POINTER :: Expo(:, :, :), Expo_Alpha(:, :, :)
INTERFACE

SUBROUTINE TranCmfd_Driver(Core, FmInfo, CmInfo, TranInfo, THInfo, GroupInfo, TranCntl, nTracerCntl, ItrCntl, PE)
USE PARAM
USE TYPEDEF,           ONLY : CoreInfo_Type,           FmInfo_Type,           CmInfo_Type,         &
                              TranInfo_Type,           GroupInfo_Type,        TranCntl_Type,       &
                              ThInfo_Type,             PE_Type
USE CNTL,              ONLY : nTracerCntl_Type
USE ItrCNTL_mod,   ONLY : ItrCntl_type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(ItrCntl_Type) :: ItrCntl
TYPE(PE_Type) :: PE

END SUBROUTINE

SUBROUTINE SetTranCmfdLinearSystem(lCmfd, l3dim, Axsolver, TranCntl)
  USE PARAM
USE TYPEDEF,        ONLY : CMFDLS_TYPE,      TranCntl_Type
IMPLICIT NONE
TYPE(TranCntl_Type) :: TranCntl
LOGICAL :: lcmfd, l3dim
INTEGER :: AxSolver

END SUBROUTINE

SUBROUTINE HomKineticParamGen(Core, FmInfo, CmInfo, TranInfo, GroupInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,    FmInfo_Type,      CmInfo_Type,          &
                             TranInfo_Type,    GroupInfo_Type,   PE_Type,              &
                             FXRInfo_Type,     PinXs_Type,       Pin_Type,             &
                             PinInfo_Type,     Cell_Type
USE CNTL,             ONLY : nTracerCntl_Type
                             
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) ::  GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

ENDSUBROUTINE

SUBROUTINE SetCmfdPrecParam(Core, CmInfo, TranInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,       ONLY : CoreInfo_Type,          CmInfo_Type,          TranInfo_Type,     &
                          TranCntl_Type,          PE_Type,                                 &
                          PinXs_Type
USE CNTL,          ONLY : nTracerCntl_Type
USE itrcntl_mod,   ONLY : ItrCntl_TYPE
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCNtl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl

TYPE(PE_TYPE) :: PE

END SUBROUTINE



SUBROUTINE CmfdPrecSrcUpdt(PrecSrc, Prec, TranPsi, TranPsid, TranCntl)
USE PARAM
USE TYPEDEF,      ONLY : TranCntl_Type
IMPLICIT NONE
REAL, POINTER :: PrecSrc(:, :), Prec(:, :, :), TranPsi(:, :), TranPsid(:, :)
TYPE(TranCntl_Type) :: TranCntl
END SUBROUTINE

SUBROUTINE CmfdSteadySrcUpdt(SRC, psifm, phifm, eigv, ig, nThread)
USE PARAM
IMPLICIT NONE
REAL, POINTER :: SRC(:, :), psifm(:, :), phifm(:, :, :)
REAL :: Eigv
INTEGER :: ig, nThread
END SUBROUTINE

SUBROUTINE CmfdSteadyPsiUpdt(phifm, psifm, nThread)
USE PARAM
REAL, POINTER :: psifm(:, :), phifm(:, :, :)
INTEGER :: nThread

END SUBROUTINE

SUBROUTINE CmfdTranSrc(TranSrc, Phi, TranPhi, Psi, PrecSrc, ResSrc, TranCntl, ig, nThread)
USE PARAM
USE TYPEDEF,      ONLY : TranCntl_Type
IMPLICIT NONE

REAL, POINTER :: TranSrc(:, :), PrecSrc(:, :), Psi(:, :), TranPhi(:, :, :), Phi(:, :, :), ResSrc(:, :, :)
TYPE(TranCntl_Type) :: TranCntl
INTEGER :: ig, nThread

END SUBROUTINE


SUBROUTINE CmfdTranSrc1g(TranSrc, Phi, TranPhi, Psi, PrecSrc, ResSrc, TranCntl)
USE PARAM
USE TYPEDEF,      ONLY : TranCntl_Type
USE CMFD_MOD,     ONLY : ng,            nxy,                                               &
                         myzb,          myze,            myzbf,           myzef
IMPLICIT NONE

REAL, POINTER :: TranSrc(:, :), PrecSrc(:, :), Psi(:, :), TranPhi(:, :, :), phi(:, :, :), ResSrc(:, :, :)
TYPE(TranCntl_Type) :: TranCntl
 
END SUBROUTINE

SUBROUTINE CmfdPrecUpdt(Prec, Psi, TranPsi, TranPsid, TranCntl)
USE PARAM
USE TYPEDEF,      ONLY : TranCntl_Type
IMPLICIT NONE
REAL, POINTER :: Prec(:, :, :), Psi(:, :), TranPsi(:, :), TranPsid(:, :)
TYPE(TranCntl_Type) :: TranCntl

END SUBROUTINE

FUNCTION TranResidualError(phifm, psifm, TranPhi, PrecSrc, ResSrc, TranCntl, PE)
USE PARAM
USE TYPEDEF,   ONLY : PE_TYPE, TranCntl_Type
IMPLICIT NONE
REAL :: TranResidualError
REAL, POINTER :: phifm(:, :, :), psifm(:, :)
REAL, POINTER :: PrecSrc(:, :), TranPhi(:, :, :), ResSrc(:, :, :)
TYPE(PE_TYPE) :: PE
TYPE(TranCntl_Type) :: TranCntl

END FUNCTION

FUNCTION TranResidualError_rev(phifm, psifm, TranPhi, PrecSrc, ResSrc, TranCntl, PE)
USE PARAM
USE TYPEDEF,   ONLY : PE_TYPE, TranCntl_Type
IMPLICIT NONE
REAL :: TranResidualError_rev
REAL, POINTER :: phifm(:, :, :), psifm(:, :)
REAL, POINTER :: PrecSrc(:, :), TranPhi(:, :, :), ResSrc(:, :, :)
TYPE(PE_TYPE) :: PE
TYPE(TranCntl_Type) :: TranCntl

END FUNCTION



END INTERFACE  

CONTAINS

SUBROUTINE SetTranCmfdEnv(CoreInfo, TranInfo, PE)
USE PARAM
USE TYPEDEF,     ONLY : CoreInfo_Type,  TranInfo_Type,   PE_TYPE
USE allocs
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(PE_Type) :: PE

INTEGER :: nxy, myzbf, myzef

nprec = TranInfo%nprec

nxy = CoreInfo%nxy; MYZEF = PE%myzef; myzbf = PE%myzbf

Chid => TranInfo%Chid
Lambda => TranInfo%Lambda
!neut_velo => TranInfo%neut_velo

CellOmegam => TranInfo%CellOmegam
Cellomega0 => TranInfo%CellOmega0
Cellomegap => TranInfo%CellOmegap

Expo => TranInfo%Expo
Expo_Alpha => TranInfo%Expo_Alpha
!
!CALL Dmalloc0(TrSrc, 1, nxy, myzbf, myzef)
!CALL Dmalloc0(PrecSrc, 1, nxy, myzbf, myzef)
END SUBROUTINE


END MODULE
