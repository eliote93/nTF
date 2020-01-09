#include <defines.h>
#ifdef __GAMMA_TRANSPORT
MODULE GamCMFD_mod
USE PARAM
USE GammaTYPEDEF,     ONLY : GammaCMFD_Type, GamMacXS_TYPE, GPINXS_TYPE
IMPLICIT NONE
! Macroscopic XS type for homogenization
INTEGER, PARAMETER :: nMacMax = 5000
TYPE(GamMacXS_TYPE) :: GXsMac(nMacMAx)
TYPE(GPINXS_TYPE), POINTER ::CMFD_GPINXS(:,:)
TYPE(GammaCMFD_Type) :: GammaCMFD

INTEGER :: ng, nxy, myzb, myze, myzbf, myzef, nz, nzfm

INTERFACE

! SUBROUTINES FOR CMFD CALCULATION

! HOMOGENIZATION ROUTINE
SUBROUTINE GamHomoXsGen(Core, FXR, Phis, GPhis, PinXS, myzb, myze, ng, ngg, lScat1, lsigT)
USE TYPEDEF,        ONLY : CoreInfo_Type,    FXRInfo_Type
USE GammaTYPEDEF,      ONLY : GPINXS_TYPE
IMPLICIT NONE
TYPE(CoreInfo_Type), INTENT(IN) :: Core
TYPE(FxrInfo_Type), POINTER , INTENT(IN):: FXR(:, :)
REAL, POINTER, INTENT(IN) :: Phis(:, :, :), GPhis(:, :, :)
TYPE(GPINXS_TYPE), POINTER, INTENT(INOUT) :: PinXs(:, :)
INTEGER, INTENT(IN) :: myzb, myze, ng, ngg
LOGICAL, INTENT(IN) :: lScat1, lsigT
END SUBROUTINE

! CELL CROSS SECTION ROUTINE
SUBROUTINE GamCellXsGen(CellInfo, Phis, GPhis, PinXs, GXSMac, ng, ngg, nFxr, nFsr, OutScatRange, lsigT)
USE TYPEDEF,          ONLY : Cell_Type
USE GammaTYPEDEF,        ONLY : GamMacXS_TYPE, GPINXS_TYPE
IMPLICIT NONE
TYPE(Cell_Type) :: CellInfo
REAL :: phis(:, :), GPhis(:, :)
TYPE(GPINXS_TYPE) :: PinXs
TYPE(GamMacXS_TYPE) :: GXSMac(nFxr)
INTEGER :: ng, ngg, nFXR, nFsr
INTEGER, POINTER :: OutScatRange(:, :)
LOGICAL :: lsigt
END SUBROUTINE

SUBROUTINE AllocGPINXS(PINXS, ng, ngg, nbd, InScatRange)
! subroutine to allocate gamma pin-wise homogenized cross section
USE GammaTYPEDEF,     ONLY : GPINXS_TYPE
IMPLICIT NONE
! INPUT VARIABLES
TYPE(GPINXS_TYPE) :: PINXS
INTEGER :: ng, ngg, nbd
INTEGER :: InScatRange(2,ngg)
END SUBROUTINE

SUBROUTINE RadCouplingCoeffGen_Gamma(Core, CmfdPinXS, Jout, ng, lDhatUpdt, lScat1, PE)
USE TYPEDEF,      ONLY : CoreInfo_Type,       PE_TYPE
USE GammaTYPEDEF, ONLY : GPinXs_Type
IMPLICIT NONE

TYPE(CoreINfo_Type) :: Core
TYPE(GPinXs_Type),POINTER :: CmfdPinXS(:, :)
REAL, POINTER :: Jout(:, :, :, :, :)
TYPE(PE_TYpe) :: PE
INTEGER :: ng   ! Photon Energy Group
LOGICAL :: lDhatUpdt, lScat1
END SUBROUTINE

!--- CNJ Edit : Operator Setup Routine for Gamma CMFD Calculation

SUBROUTINE SetGammaCMFDSystem(CoreInfo, GroupInfo, PinXS, GammaCMFD, lDhat)
USE PARAM
USE TYPEDEF,		    ONLY : CoreInfo_Type
USE GammaTYPEDEF,       ONLY : gPinXS_Type,         GammaGroupInfo_Type,        GammaCMFD_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(GammaGroupInfo_Type) :: GroupInfo
TYPE(gPinXS_Type), POINTER :: PinXS(:, :)
TYPE(GammaCMFD_Type) :: GammaCMFD
LOGICAL :: lDhat

END SUBROUTINE

!--- CNJ Edit : Solver Routine for Gamma CMFD Calculation

SUBROUTINE GammaBiCGSTAB(M, ILU, x, b)
USE CSRMATRIX,          ONLY : CSR_DOUBLE
IMPLICIT NONE

TYPE(CSR_DOUBLE) :: M, ILU
REAL :: x(:, :), b(:)

END SUBROUTINE

!--- CNJ Edit : Homogenized Flux Mapping

SUBROUTINE UpdateCellPhi(CoreInfo, PinXS, GroupInfo, GammaCMFD)
USE PARAM
USE TYPEDEF,            ONLY : CoreInfo_Type
USE GammaTYPEDEF,       ONLY : gPinXS_Type,         GammaCMFD_Type,     GammaGroupInfo_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(gPinXS_Type), POINTER :: PinXS(:, :)
TYPE(GammaGroupInfo_Type) :: GroupInfo
TYPE(GammaCMFD_Type) :: GammaCMFD

END SUBROUTINE

!--- CNJ Edit : Gamma MOC Solution Update

SUBROUTINE UpdateFsrPhi(CoreInfo, CmInfo, PinXS, GroupInfo, GammaCMFD, gPhis, gPhiAngIn)
USE PARAM
USE TYPEDEF,            ONLY : CoreInfo_Type,       CmInfo_Type,        RayInfo_Type,                               &
                               Pin_Type,            Cell_Type
USE GammaTYPEDEF,       ONLY : gPinXS_Type,         GammaCMFD_Type,     GammaGroupInfo_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(gPinXS_Type), POINTER :: PinXS(:, :)
TYPE(GammaGroupInfo_Type) :: GroupInfo
TYPE(GammaCMFD_Type) :: GammaCMFD
REAL, POINTER :: gPhis(:, :, :), gPhiAngIn(:, :, :, :)

END SUBROUTINE

END INTERFACE

END MODULE
#endif
