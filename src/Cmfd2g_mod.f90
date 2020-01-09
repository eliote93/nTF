MODULE Cmfd2g_Mod
USE PARAM
IMPLICIT NONE

REAL, POINTER :: PhiC2g(:, :, :)
REAL, POINTER :: Src2g(:, :, :)

INTERFACE
SUBROUTINE Cmfd2GSrcUpdt(SRC, psi, reigv)
IMPLICIT NONE
REAL, POINTER :: SRC(:, :, :), psi(:, :)
REAL :: rEigv
END SUBROUTINE

SUBROUTINE ConvertArray2G(PhiC2G, GcPhiC, nxy, iz1, iz2, imod)
USE PARAM
IMPLICIT NONE
REAL :: PhiC2G(2, nxy, iz1:iz2), GcPhiC(nxy, iz1:iz2, 2)
INTEGER :: iz1, iz2, nxy, imod
INTEGER :: iz, ixy

END SUBROUTINE

SUBROUTINE SetCmfd2GLinearSystem(lcmfd, l3dim)
IMPLICIT NONE
LOGICAL :: lCmfd, l3dim
END SUBROUTINE



FUNCTION ResidualError2G(phi, psi, reigv, reigvs, PE, constsrc)
USE PARAM
USE TYPEDEF,   ONLY : PE_TYPE
IMPLICIT NONE
REAL :: ResidualError2G
REAL, POINTER :: phi(:, :, :), psi(:, :)
REAL :: reigv, reigvs
TYPE(PE_TYPE) :: PE
REAL, OPTIONAL :: ConstSrc

END FUNCTION

SUBROUTINE WielandtUpdt(CMFD2GLS, PsiC, PsiCd, Eigv, Eigvs, PsiErr, iter, PE)
USE PARAM
USE TYPEDEF,    ONLY : PE_TYPE,      CMFDLS_TYPE
IMPLICIT NONE
TYPE(CmfdLS_TYPE) :: CMFD2GLS
REAL, POINTER :: PsiC(:, :), PsiCd(:, :)
REAL :: Eigv, Eigvs, PsiErr
INTEGER :: iter
TYPE(PE_TYPE) :: PE

END SUBROUTINE

SUBROUTINE WielandtLsShift(Cmfd2GLs, reigv, reigvsdel , PE)
USE PARAM
USE TYPEDEF,    ONLY : PE_TYPE,      CMFDLS_TYPE
USE CMFD_MOD,   ONLY : nxy,         myzbf,       myzef,               &
                       SubPlaneMap, PinVol,      PinVolFm,            &
                       hzfm,        PinNeighIdx
USE GcCmfd_Mod, ONLY : ngc,       GcPinXS
IMPLICIT NONE
TYPE(CMFDLS_TYPE) :: Cmfd2Gls
TYPE(PE_TYPE) :: PE
REAL :: reigv, reigvsdel

END SUBROUTINE

SUBROUTINE NegativeFixUp2G(PhiC2G, PsiC2G, PE)
USE PARAM
USE TYPEDEF,   ONLY : PE_TYPE
IMPLICIT NONE
REAL, POINTER :: PhiC2G(:, :, :), PsiC2G(:, :)
TYPE(PE_TYPE) :: PE

END SUBROUTINE

END INTERFACE
END MODULE