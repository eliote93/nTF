MODULE TRANMOC_MOD
USE PARAM
IMPLICIT NONE
SAVE
REAL, POINTER :: TrSrc(:), PrecSrc(:, :)
REAL, POINTER :: Chid(:), lambda(:) !, neut_velo(:)
REAL, POINTER :: Omegam(:, :, :), Omega0(:, :, :), Omegap(:, :, :)
REAL, POINTER :: Expo(:, :, :), Expo_Alpha(:, :, :)

INTEGER :: nprec

INTERFACE


SUBROUTINE SetPrecParam(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,         ONLY : CoreInfo_Type,          FmInfo_Type,          TranInfo_Type,     &
                            TranCntl_Type,          PE_Type
USE CNTL,            ONLY : nTracerCntl_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCNtl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE

END SUBROUTINE 

SUBROUTINE TranMOC_Driver(Core, RayInfo, FmInfo, CmInfo, TranInfo, THInfo, GroupInfo, TranCntl, nTracerCntl, ItrCntl, PE)
USE PARAM
USE TYPEDEF,           ONLY : CoreInfo_Type,           FmInfo_Type,           CmInfo_Type,         &
                              TranInfo_Type,           GroupInfo_Type,        TranCntl_Type,       &
                              RayInfo_Type,            ThInfo_Type,           PE_Type
USE CNTL,              ONLY : nTracerCntl_Type
USE itrcntl_mod,       ONLY : ItrCntl_TYPE
IMPLICIT NONE  
TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(ItrCntl_TYPE) :: ItrCntl
TYPE(PE_Type) :: PE
  
END SUBROUTINE

SUBROUTINE PrecSrcUpdt(Core, Fxr, PrecSrc, Prec, TranPsi, TranPsid, GroupInfo, TranCntl, nTRACERCntl, PE)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,           FxrInfo_Type,            TranCntl_Type,        &
                           GroupInfo_Type,          PE_Type
USE CNTL,           ONLY : nTracerCntl_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

REAL, POINTER :: Prec(:, :, :), PrecSrc(:, :)
REAL, POINTER :: TranPsi(:, :), TranPsid(:, :)

END SUBROUTINE

SUBROUTINE PrecSrckUpdt(Core, Fxr, PrecSrcK, Prec, TranPsi, TranPsid, GroupInfo, TranCntl, nTRACERCntl, PE)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,           FxrInfo_Type,            TranCntl_Type,        &
                           GroupInfo_Type,          PE_Type
USE CNTL,           ONLY : nTracerCntl_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

REAL, POINTER :: Prec(:, :, :), PrecSrcK(:, :, :)
REAL, POINTER :: TranPsi(:, :), TranPsid(:, :)

END SUBROUTINE

SUBROUTINE SetTranSrcNM(Core, FmInfo, Fxr, TranSrcnm, Phinm, TranPhinm, Psi, ResSrc, xstnm, iz, &
                        gb, ge, GroupInfo, TranInfo, TranCntl, lxslib, PE, Offset)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,         FxrInfo_Type,         GroupInfo_Type,         TranCntl_Type,      &
                           PE_Type,               TranInfo_Type,        FmInfo_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
REAL, POINTER :: TranSrcnm(:, :), Phinm(:, :), TranPhinm(:, :), Psi(:, :), ResSrc(:, :), xstnm(:, :)
INTEGER :: iz, gb, ge
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
LOGICAL :: lxslib
TYPE(PE_Type) :: PE
INTEGER, OPTIONAL :: Offset   !--- CNJ Edit : GPU Acceleration
END SUBROUTINE

SUBROUTINE SetTranSrc(Core, Fxr, TranSrc, Phi, TranPhi, Psi, PrecSrc, ResSrc, xstr, iz, ig, GroupInfo, TranInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,            GroupInfo_Type,        TranCntl_Type,            &
                           FxrInfo_Type,             PE_Type,               TranInfo_Type
USE CNTL,           ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
INTEGER :: iz, ig

REAL, POINTER :: TranSrc(:)
REAL, POINTER :: Phi(:, :, :)
REAL, POINTER :: TranPhi(:, :, :)
REAL, POINTER :: Psi(:, :)
REAL, POINTER :: PrecSrc(:, :)
REAL, POINTER :: ResSrc(:, :, :)
REAL, POINTER :: xstr(:)
END SUBROUTINE

SUBROUTINE SetExpTrsfXs(Core, Fxr, xstr, iz, ig, GroupInfo, TranInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,            GroupInfo_Type,        TranCntl_Type,            &
                           FxrInfo_Type,             PE_Type,               TranInfo_Type
USE CNTL,           ONLY : nTracerCntl_Type

IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranInfo_TYpe) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
INTEGER :: iz, ig

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
REAL, POINTER :: xstr(:)

END SUBROUTINE

FUNCTION TranMocResidualError(Core, FmInfo, TranInfo, eigv0, GroupInfo, TranCntl, nTRACERCntl, PE)
USE PARAM
USE TYPEDEF,           ONLY : CoreInfo_Type,        FmInfo_Type,         GroupInfo_Type,         &
                              PE_Type,              TranInfo_Type,       TranCntl_Type
USE CNTL,              ONLY : nTracerCntl_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
REAL :: Eigv0
REAL :: TranMocResidualError
END FUNCTION

END INTERFACE

CONTAINS

SUBROUTINE SetTranMOCEnv(CoreInfo, TranInfo, PE)
USE PARAM
USE TYPEDEF,     ONLY : CoreInfo_Type,  TranInfo_Type,   PE_TYPE
USE allocs
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(PE_Type) :: PE

INTEGER :: nCorefsr, myzb, myze

nprec = TranInfo%nprec

nCoreFsr = CoreInfo%nCoreFsr; myze = PE%myze; myzb = PE%myzb

Chid => TranInfo%Chid
Lambda => TranInfo%Lambda
!neut_velo => TranInfo%neut_velo

Omegam => TranInfo%FxrOmegam
Omega0 => TranInfo%FxrOmega0
Omegap => TranInfo%FxrOmegap

Expo => TranInfo%FmExpo
Expo_Alpha => TranInfo%FmExpo_Alpha

IF(.NOT. ASSOCIATED(TrSrc)) THEN
  CALL Dmalloc0(TrSrc, 1, nCoreFsr)
ENDIF
IF(.NOT. ASSOCIATED(PrecSrc)) THEN
  CALL Dmalloc0(PrecSrc, 1, nCoreFsr, myzb, myze)
ENDIF
END SUBROUTINE
END MODULE
