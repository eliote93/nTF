#include <defines.h>
#ifdef __GAMMA_TRANSPORT
MODULE GamMOC_MOD
USE GammaTYPEDEF,       ONLY : GammaTrackingDat_Type
IMPLICIT NONE

REAL, POINTER :: gphis1g(:), gphim1g(:, :), gxst1g(:), gsrc1g(:), gsrcm1g(:, :)
REAL, POINTER :: gPhiAngIn1g(:, :), gJout1g(:, :, :)
REAL, POINTER :: gAxSrc1g(:), gAxPXS1g(:)

REAL, POINTER :: gphisnm(:, :), gxstnm(:, :), gsrcnm(:, :)
REAL, POINTER :: gPhiAngInnm(:, :, :), gJoutnm(:, :, :, :)

REAL, POINTER :: Comp(:, :, :), mwt(:, :, :), wtang(:, :), wtsurf(:, :, :)

TYPE(GammaTrackingDat_Type), POINTER :: TrackingDat(:)

INTERFACE

SUBROUTINE SetGamMacXs(core, Fxr, xstr, iz, igg, ngg, lTrCorrection, PE)
! Photo-Atomic Reaction Cross Section Generation Routine
!   Only 0 K  case is treated with Photo-atomic reaction -> No resonance treatment
!   Only Outflow Correction is provided for PHOTON TRANSPORT CROSS SECTION
USE PARAM
USE TYPEDEF,  ONLY : coreinfo_type, Fxrinfo_type, PE_TYPE

IMPLICIT NONE

TYPE(coreinfo_type) :: Core
TYPE(Fxrinfo_type) :: Fxr(:)
TYPE(PE_TYPE) :: PE
REAL, POINTER :: xstr(:)

INTEGER :: iz, igg, ngg
logical :: lGamma, lTrCorrection
END SUBROUTINE

SUBROUTINE SetGamSrc(Core, Fxr, src, phis, gphis, gaxsrc1g, xstr1g, iz, igg, ng, ngg,                  &
                        l3dim, lscat1, PE, GroupInfo)
! Photon Source Generation Routine in Photon Transport Equation
!         Production by Neutron-Isotope Reaction and Photo-Atomic Reaction
USE PARAM
USE TYPEDEF,        ONLY : coreinfo_type, Fxrinfo_type, PE_TYPE,  GROUPINFO_TYPE
IMPLICIT NONE
! INPUT VARIABLES
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:)
TYPE(PE_TYPE) :: PE
TYPE(GROUPINFO_TYPE) :: GroupInfo

REAL, POINTER :: src(:), phis(:, :, :), gphis(:, :, :), gAxSrc1g(:), xstr1g(:)
INTEGER :: igg, ng, iz, ifsr, ifxr, ngg
LOGICAL :: lscat1, l3dim

END SUBROUTINE

SUBROUTINE SetGamP1Src(Core, Fxr, srcm, phim, xstr1g, iz, igg, ngg, lscat1, ScatOd, PE)
USE PARAM
USE TYPEDEF,        ONLY : coreinfo_type,          Fxrinfo_type,          PE_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:)
REAL, POINTER :: srcm(:, :)
REAL, POINTER :: phim(:, :, :, :)
REAL, POINTER :: xstr1g(:)
INTEGER :: myzb, myze, igg, ngg, iz, ifsr, ifxr
LOGICAL :: lscat1
INTEGER :: ScatOd
TYPE(PE_Type) :: PE

END SUBROUTINE

SUBROUTINE SetGamMacXsNM(core, Fxr, xstnm, iz, ngg, lTrCorrection, PE)
! Photo-Atomic Reaction Cross Section Generation Routine
!   Only 0 K  case is treated with Photo-atomic reaction => No resonance treatment
!   Only Outflow Correction is provided for PHOTON TRANSPORT CROSS SECTION
!   For Node Major
USE PARAM
USE TYPEDEF,          ONLY : coreinfo_type, Fxrinfo_type, PE_TYPE
IMPLICIT NONE

TYPE(CoreInfo_TYPE) :: Core
TYPE(FxrInfo_TYPE) :: Fxr(:)
TYPE(PE_TYPE) :: PE
REAL, POINTER :: xstnm(:, :)
INTEGER :: iz, ngg
LOGICAL :: lGamma, lTrCorrection

END SUBROUTINE

SUBROUTINE SetGamSrcNM(Core, Fxr, srcNM, phisNM, gphisNM, gAxSrc, xstnm, iz, igb, ige,               &
                       ng, ngg, l3dim, lscat1, PE)
! Photon Source Generation Routine in Photon Transport Equation
!      FOR NODE MAJOR
USE PARAM
USE TYPEDEF,          ONLY : coreinfo_type, Fxrinfo_type, PE_TYPE
IMPLICIT NONE
! INPUT VARIABLES
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:)
REAL, POINTER :: srcNM(:,:), phisNM(:, :), gphisNM(:, :), gAxSrc(:, :, :), xstnm(:, :)
REAL :: eigv
INTEGER :: iz, igb, ige, ng, ngg, ifsr, ifxr, fsridx
LOGICAL :: l3dim, lscat1
TYPE(PE_TYPE) :: PE

END SUBROUTINE

SUBROUTINE GamPseudoAbsorption(Core, Fxr, AxPXS, xstr1g, iz, l3dim)
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,          Fxrinfo_type,          Cell_Type,      &
                         pin_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:)
REAL, POINTER :: AxPXS(:), xstr1g(:)
INTEGER :: iz
LOGICAL :: l3dim

END SUBROUTINE

SUBROUTINE GamPseudoAbsorptionNM(Core, Fxr, AxPXS, xstnm, iz, ng, l3dim)
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,          Fxrinfo_type,          Cell_Type,      &
                         pin_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:)
REAL, POINTER :: AxPXS(:, :, :), xstnm(:, :)
INTEGER :: iz, ng
LOGICAL :: l3dim

END SUBROUTINE

SUBROUTINE RayTraceGamma(RayInfo, CoreInfo, phisnm, PhiAngInnm, xstnm, srcnm,                 &
                         joutnm, iz, gb, ge, ljout)
USE PARAM
USE TYPEDEF, ONLY :   RayInfo_Type,      CoreInfo_type
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phisnm(:, :), PhiAngInnm(:, :, :), xstnm(:, :), srcnm(:, :), joutnm(:, :, :, :)
INTEGER :: iz, gb, ge
LOGICAL :: ljout

END SUBROUTINE

SUBROUTINE RayTraceGamma_Pn(RayInfo, CoreInfo, phis, phim, PhiAngIn, xst, src, srcm, jout, iz, ScatOd, lJout)
USE PARAM
USE TYPEDEF, ONLY :   RayInfo_Type,      CoreInfo_type
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phis(:), phim(:, :), PhiAngIn(:, :), xst(:)
REAL, POINTER :: src(:), srcm(:, :), jout(:, :, :)
INTEGER :: iz
LOGICAL :: ljout
INTEGER :: ScatOd

END SUBROUTINE

FUNCTION GammaMocResidual(Core, FmInfo, GroupInfo, PE, nTracerCntl)
USE PARAM
USE TYPEDEF,  ONLY : CoreInfo_Type,     FmInfo_Type,    GroupInfo_Type,     PE_Type
USE CNTL,     ONLY : nTracerCntl_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_TYPE) :: PE
TYPE(nTracerCntl_Type) :: nTracerCntl
REAL :: GammaMocResidual

END FUNCTION

SUBROUTINE GammaPowerUpdate(Core, Fxr, gphis, gpower, myzb, myze, ng, PE)
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,   Fxrinfo_type,     PE_TYPE
IMPLICIT NONE

TYPE(coreinfo_type) :: CORE
TYPE(Fxrinfo_type),POINTER :: Fxr(:, :)
REAL, POINTER :: gphis(:, :, :)
REAL, POINTER :: gpower(:, :)
INTEGER :: myzb, myze, ng
TYPE(PE_Type) :: PE

END SUBROUTINE

SUBROUTINE NeutronLocalQUpdate(Core, Fxr, phis, localpower, GPOWERGEN, myzb, myze, ng, PE)
USE PARAM
USE TYPEDEF,        ONLY : coreinfo_type,       Fxrinfo_type,       PE_TYPE
IMPLICIT NONE
TYPE(coreinfo_type) :: CORE
TYPE(Fxrinfo_type),POINTER :: Fxr(:, :)
TYPE(PE_Type) :: PE
REAL, POINTER :: phis(:, :, :)
REAL, POINTER :: localpower(:, :), GPOWERGEN(:, :)
INTEGER :: myzb, myze, ng

END SUBROUTINE


FUNCTION FxrAvgGPhi(Core, Fxr, GPhis, ipin, iLocalfxr, iz, ng, PE)
USE PARAM
USE TYPEDEF,     ONLY : CoreInfo_Type,    PE_Type,     FxrInfo_Type,               &
                        Cell_Type,    Pin_Type
USE CNTL,        ONLY : nTRACERCNTL
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(PE_Type) :: PE
REAL, POINTER :: GPhis(:, :, :)
INTEGER :: iLocalfxr, ipin, iz, ng
REAL :: FxrAvgGPhi(ng)
END FUNCTION

SUBROUTINE CompensateGPower(Core,Fxr,GPowerGen,GPower, myzb, myze,PE)
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,       Fxrinfo_type,          PE_TYPE

IMPLICIT NONE

TYPE(coreinfo_type) :: CORE
TYPE(Fxrinfo_type),POINTER :: Fxr(:, :)
TYPE(PE_Type) :: PE
REAL, POINTER, DIMENSION(:, :) :: GPowerGen, gpower
INTEGER :: myzb, myze
END SUBROUTINE

END INTERFACE

END MODULE
#endif
