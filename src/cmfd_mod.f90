MODULE CMFD_MOD
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_TYPE, PinXs_TYPE, AxFlx_Type, XsMac_Type, CMFDLS_TYPE

SAVE
!Global Variables


!CMFD module local variables

TYPE(XsMac_Type) :: XsMac(5000)
TYPE(PinXS_Type), POINTER :: CmfdPinXS(:, :)
TYPE(CMFDLS_TYPE), POINTER ::CMFDLS1g
TYPE(CMFDLS_TYPE), POINTER :: CMFDLS(:)
!TYPE(AxFlx_Type), POINTER :: AxFlx(:, :)

REAL, POINTER :: SRC(:, :), Phic1g(:, :)
INTEGER, POINTER :: PinNeighIdx(:, :)
INTEGER, POINTER :: SubPlaneMap(:), SubPlaneRange(:, :)
INTEGER :: nSubplane
INTEGER, POINTER :: nSubPlanePtr(:)
INTEGER :: ng, nxy, myzb, myze, myzbf, myzef, nz, nzfm
REAL, POINTER :: hz(:), hzfm(:)
REAL, POINTER :: PinVol(:, :), PinVolFm(:, :)

REAL, POINTER :: AxDhat(:, :, :, :), AxPDhat(:, :, :, :), AxDtil(:, :, :, :)

INTERFACE

SUBROUTINE AllocPinXS(PinXS, ng, nbd, InScatRange)
USE PARAM
USE TYPEDEF, ONLY : PinXs_Type
IMPLICIT NONE
TYPE(PinXs_Type) :: PinXs
INTEGER :: ng, nbd
INTEGER :: InScatRange(2, ng)
END SUBROUTINE

SUBROUTINE SetCMFDEnviorment(CORE, CmInfo, ng0,PE,  lDcplPlnCmfd)
USE PARAM
USE TYPEDEF, ONLY : PE_TYPE, CoreInfo_Type, CmInfo_Type
!USE PE_MOD, ONLY : PE
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
TYPE(CmInfo_Type) :: CmInfo
TYPE(CoreInfo_Type) :: Core
INTEGER :: ng0
LOGICAL, OPTIONAL :: lDcplPlnCmfd
END SUBROUTINE


SUBROUTINE HomoXsGen(Core, FXR, Phis, PinXS, myzb, myze, ng, lXsLib, lScat1, lsigT)
USE TYPEDEF,  ONLY : CoreInfo_Type,    FXRInfo_Type,      PinXs_Type
IMPLICIT NONE

TYPE(CoreInfo_Type), INTENT(IN) :: Core
TYPE(FxrInfo_Type), POINTER, INTENT(IN) :: FXR(:, :)
TYPE(PinXs_Type), POINTER, INTENT(INOUT) :: PinXs(:, :)
REAL, POINTER, INTENT(IN) :: Phis(:, :, :)
INTEGER, INTENT(IN) :: myzb, myze, ng
LOGICAL, INTENT(IN) :: lXsLib, lScat1, lSigT

END SUBROUTINE

SUBROUTINE HomoXsGen_Cusping(Core, FmInfo, FXR, Phis, PinXS, myzb, myze, ng, lXsLib, lScat1, lsigT)
USE TYPEDEF,  ONLY : CoreInfo_Type,    FXRInfo_Type,      PinXs_Type,     FmInfo_Type
IMPLICIT NONE

TYPE(CoreInfo_Type), INTENT(IN) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(FxrInfo_Type), POINTER, INTENT(IN) :: FXR(:, :)
TYPE(PinXs_Type), POINTER, INTENT(INOUT) :: PinXs(:, :)
REAL, POINTER, INTENT(IN) :: Phis(:, :, :)
INTEGER, INTENT(IN) :: myzb, myze, ng
LOGICAL, INTENT(IN) :: lXsLib, lScat1, lSigT

END SUBROUTINE
SUBROUTINE HomoCellXsGen(CellInfo, Phis, PinXs, XsMac, ng, nFxr, nFsr, OutScatRange, lXsLib, lsigT)
USE TYPEDEF, ONLY : Cell_Type, PinXs_Type, XsMac_Type
IMPLICIT NONE
TYPE(Cell_Type) :: CellInfo
REAL :: phis(:, :)
TYPE(PinXs_Type) :: PinXs
TYPE(XsMac_Type) :: XsMac(nFxr)
INTEGER :: ng, nFsr, nFXR
INTEGER, POINTER :: OutScatRange(:, :)
LOGICAL :: lXsLib, lsigT

END SUBROUTINE
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RadCouplingCoeffGen(Core, CmfdPinXS, Jout, ng, lDhatUpdt, PE)

USE TYPEDEF, ONLY : CoreInfo_Type, PinXs_Type, PE_TYPE, Pin_Type

IMPLICIT NONE

TYPE (CoreINfo_Type) :: Core
TYPE (PE_TYpe)       :: PE
TYPE (PinXs_Type), POINTER, DIMENSION(:,:) :: CmfdPinXS

REAL, POINTER, DIMENSION(:, :, :, :, :) :: Jout

INTEGER :: ng
LOGICAL :: lDhatUpdt

END SUBROUTINE RadCouplingCoeffGen
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetCmfdLinearSystem(lcmfd, l3dim, AxSolver)
IMPLICIT NONE
LOGICAL :: lcmfd, l3dim
INTEGER :: AxSolver
END SUBROUTINE

SUBROUTINE UpdatePhiC(PinXS, PhiC)
USE TYPEDEF,  ONLY : PinXs_Type
IMPLICIT NONE
TYPE(PinXs_Type), POINTER :: PinXs(:, :)
REAL, POINTER :: PhiC(:, :, :)
END SUBROUTINE


SUBROUTINE ConvertSubPlnPhi(PhiC, PhiFm, imod)
IMPLICIT NONE
REAL, POINTER :: PhiC(:, :, :), PhiFm(:, :, :)
INTEGER :: imod

END SUBROUTINE

SUBROUTINE CmfdSrcUpdt(SRC, psifm, phifm, eigv, ig)
IMPLICIT NONE
REAL, POINTER :: SRC(:, :), psifm(:, :), phifm(:, :, :)
REAL :: Eigv
INTEGER :: ig

END SUBROUTINE


SUBROUTINE CmfdPsiUpdt(phifm, psifm)
IMPLICIT NONE
REAL, POINTER :: psifm(:, :), phifm(:, :, :)

END SUBROUTINE

SUBROUTINE CmfdEigUpdate(psifm, psifmd, eigv, psierr, PE)
USE TYPEDEF,   ONLY : PE_TYPE
IMPLICIT NONE
REAL, POINTER :: psifm(:, :), psifmd(:, :)
REAL :: psierr, eigv
TYPE(PE_TYPE), OPTIONAL :: PE
END SUBROUTINE



SUBROUTINE MOCSolUpdt(Core, FmInfo, CmInfo, myzb, myze, ng)
USE PARAM
USE TYPEDEF,  ONLY : CoreInfo_Type,     FmInfo_Type,      CmInfo_Type,   &
                     FXRInfo_Type,      PinXs_Type,                      &
                     Pin_Type,          PinInfo_Type,      Cell_Type
IMPLICIT NONE

TYPE(CoreInfo_Type), INTENT(IN) :: Core
TYPE(FmInfo_Type), INTENT(INOUT) :: FmInfo
TYPE(CmInfo_Type), INTENT(INOUT) :: CmInfo
INTEGER, INTENT(IN) :: myzb, myze, ng

END SUBROUTINE
!SUBROUTINE MOCSolUpdt(Core, PinXS, PHIC, phis, myzb, myze, ng)
!USE TYPEDEF,  ONLY : CoreInfo_Type,    PinXs_Type
!IMPLICIT NONE
!
!TYPE(CoreInfo_Type), INTENT(IN) :: Core
!TYPE(PinXs_Type), POINTER, INTENT(INOUT) :: PinXs(:, :)
!REAL :: Phis(:, :, :), PHIC(:, : , :)
!INTEGER :: myzb, myze, ng
!
!END SUBROUTINE

SUBROUTINE MOCLinSrcUpdt(Core, PinXS, PHIC, LinSrcSlope, myzb, myze, ng)
USE TYPEDEF,  ONLY : CoreInfo_Type,    PinXs_Type
IMPLICIT NONE

TYPE(CoreInfo_Type), INTENT(IN) :: Core
TYPE(PinXs_Type), POINTER, INTENT(INOUT) :: PinXs(:, :)
REAL, POINTER :: LinSrcSlope(:, :, :, :), PHIC(:, : , :)
INTEGER :: myzb, myze, ng

END SUBROUTINE

SUBROUTINE MOCPhiInUpdt(Core, CmInfo, myzb, myze, ng)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type,  CmInfo_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CMInfo
INTEGER :: myzb, myze, ng

END SUBROUTINE

SUBROUTINE AddConstCmfdSrc(Src, ConstSrc)
USE PARAM
IMPLICIT NONE
REAL, POINTER :: Src(:,:)
REAL :: ConstSrc

END SUBROUTINE

SUBROUTINE AddCmfdPxs(Pxs, iz1, iz2, ig1, ig2)
USE PARAM
IMPLICIT NONE
REAL, POINTER :: PXS(:, :, :)
INTEGER :: iz1, iz2, ig1, ig2

END SUBROUTINE
SUBROUTINE HomBuckling(Core, Fxr, Phis, PinXS, myzb, myze, ng, bsq, lxslib)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,    FXRInfo_Type,      PinXs_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(Fxrinfo_type), POINTER :: Fxr(:,:)
TYPE(PinXs_Type), POINTER :: PinXs(:, :)
REAL, POINTER :: Phis(:, :, :)
INTEGER :: myzb, myze, ng
LOGICAL :: lXsLib
REAL :: Bsq

END SUBROUTINE

FUNCTION ResidualError(phifm, psifm, eigv, ig1, ig2, PE, ConstSrc)
USE TYPEDEF,    ONLY : PE_TYPE
IMPLICIT NONE
REAL :: ResidualError
REAL, POINTER :: phifm(:, :, :), psifm(:, :)
REAL :: eigv
TYPE(PE_TYPE) :: PE
INTEGER :: ig1, ig2
REAL, OPTIONAL :: ConstSrc
END FUNCTION


END INTERFACE

END MODULE