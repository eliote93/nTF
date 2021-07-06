MODULE MOC_MOD
  
USE PARAM
USE TYPEDEF, ONLY : TrackingDat_Type

IMPLICIT NONE

REAL, POINTER, DIMENSION(:)     :: phis1g, xst1g, tsrc, axSrc1g, axPxs1g
REAL, POINTER, DIMENSION(:,:)   :: phiAngIn1g, phim1g, LinSrc1g, srcm
REAL, POINTER, DIMENSION(:,:,:) :: phi1a1g, phi2a1g, MocJout1g, tsrcAng1, tsrcAng2, LinPsi

! Node Major
REAL, POINTER, DIMENSION(:,:)     :: phisNM, srcNM, xstNM
REAL, POINTER, DIMENSION(:,:,:)   :: phimNM, phiAngInNM, srcmNM
REAL, POINTER, DIMENSION(:,:,:,:) :: phiaNM, srcAngNM, MocJoutNM

! Domain Decomposition
REAL, POINTER, DIMENSION(:,:,:,:,:) :: DcmpPhiAngIn, DcmpPhiAngOut
REAL, POINTER, DIMENSION(:,:,:,:)   :: DcmpPhiAngIn1g, DcmpPhiAngOut1g

INTEGER, POINTER, DIMENSION(:,:) :: DcmpColorAsy

! Track Rotational Ray
REAL, POINTER, DIMENSION(:,:)     :: wtang, hwt
REAL, POINTER, DIMENSION(:,:,:)   :: wtsurf, SrcAng1, SrcAng2, comp, mwt, mwt2, phia1g
REAL, POINTER, DIMENSION(:,:,:,:) :: SrcAng, SrcAngNM1, SrcAngNM2

INTEGER :: nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay

INTEGER, TARGET :: AziMap(360, 2)

TYPE(TrackingDat_Type), SAVE :: TrackingDat(100)

! AFSS
INTEGER :: nOmpAng

INTEGER, POINTER, DIMENSION(:)   :: OmpRayBeg, OmpRayEnd
INTEGER, POINTER, DIMENSION(:,:) :: OmpRayBegBd, OmpRayEndBd, OmpMap

! Approximation of Exponetial Function
REAL, TARGET, DIMENSION(-40000:0, 1:12) :: expa, expb

INTERFACE

! ------------------------------------------------------------------------------------------------------------
SUBROUTINE MOCSweep(Core, RayInfo, FMInfo, eigv, ng, PE, nTracerCntl, ItrCntl)

USE TYPEDEF,     ONLY : CoreInfo_Type, RayInfo_Type, FmInfo_Type, PE_TYPE
USE CNTL,        ONLY : nTracerCntl_Type
USE itrcntl_mod, ONLY : ItrCntl_TYPE

IMPLICIT NONE

TYPE(CoreInfo_Type)    :: Core
TYPE(RayInfo_Type)     :: RayInfo
TYPE(FmInfo_Type)      :: FmInfo
TYPE(PE_TYPE)          :: PE
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(ItrCntl_TYPE)     :: ItrCntl

REAL    :: eigv
INTEGER :: ng

END SUBROUTINE MOCSweep
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE DcplMOCSweep(Core, RayInfo, FMInfo, eigv, ng, PE, GroupInfo, nTracerCntl, ItrCntl)

USE TYPEDEF,     ONLY : CoreInfo_Type, RayInfo_Type, FmInfo_Type, PE_TYPE, GroupInfo_Type
USE CNTL,        ONLY : nTracerCntl_Type
USE itrcntl_mod, ONLY : ItrCntl_TYPE

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(PE_TYPE) :: PE
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(ItrCntl_TYPE) :: ItrCntl

REAL :: eigv
INTEGER :: ng

END SUBROUTINE DcplMOCSweep
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTrace_Tran(RayInfo, CoreInfo, phis, phi1a, phi2a, PhiAngIn, xst, src, srcang1, srcang2, jout, iz, ljout, FastMocLv, lAFSS)

USE PARAM
USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type

IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL :: phis(:), PhiAngIn(:, :), xst(:), src(:), jout(:, :, :)
REAL :: phi1a(:,:,:), phi2a(:,:,:)
REAL :: srcang1(:,:,:), srcang2(:,:,:)

INTEGER :: iz
LOGICAL :: ljout
INTEGER, OPTIONAL :: FastMocLv
LOGICAL, OPTIONAL :: lAFSS

END SUBROUTINE RayTrace_Tran
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmp_GM(RayInfo, CoreInfo, phis, PhiAngIn, xst, src, MocJout, iz, lJout)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:)     :: phis, xst, src
REAL, POINTER, DIMENSION(:,:)   :: PhiAngIn
REAL, POINTER, DIMENSION(:,:,:) :: Mocjout

INTEGER :: iz
LOGICAL :: lJout

END SUBROUTINE RayTraceDcmp_GM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmp_NM(RayInfo, CoreInfo, phisNM, PhiAngInNM, xstNM, srcNM, MocJoutNM, iz, gb, ge, lJout)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:,:)     :: phisNM, xstNM, srcNM
REAL, POINTER, DIMENSION(:,:,:)   :: PhiAngInNM
REAL, POINTER, DIMENSION(:,:,:,:) :: MocjoutNM

INTEGER :: iz, gb, ge
LOGICAL :: lJout

END SUBROUTINE RayTraceDcmp_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceP1_Dcmp(RayInfo, CoreInfo, phisNM, phimNM, PhiAngInNM, xstNM, srcNM, srcmNM, MocJoutNM, iz, gb, ge, lJout)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:,:)     :: phisNM, xstNM, srcNM
REAL, POINTER, DIMENSION(:,:,:)   :: phimNM, PhiAngInNM, srcmNM
REAL, POINTER, DIMENSION(:,:,:,:) :: MocjoutNM

INTEGER :: iz, gb, ge
LOGICAL :: lJout

END SUBROUTINE RayTraceP1_Dcmp
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceLin_Dcmp(RayInfo, CoreInfo, iz, gb, ge, lJout, lHybrid)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

INTEGER :: iz, gb, ge
LOGICAL :: lJout, lHybrid

END SUBROUTINE RayTraceLin_Dcmp
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceGM_OMP(RayInfo, CoreInfo, phis, PhiAngIn, xst, src, jout, iz, ljout, FastMocLv)

USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:)     :: phis, xst, src
REAL, POINTER, DIMENSION(:,:)   :: PhiAngIn
REAL, POINTER, DIMENSION(:,:,:) :: jout

INTEGER :: iz
LOGICAL :: ljout

INTEGER, OPTIONAL :: FastMocLv

END SUBROUTINE RayTraceGM_OMP
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceGM_AFSS(RayInfo, CoreInfo, phis, PhiAngIn, xst, src, jout, iz, ljout, FastMocLv)

USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:)     :: phis, xst, src
REAL, POINTER, DIMENSION(:,:)   :: PhiAngIn
REAL, POINTER, DIMENSION(:,:,:) :: jout

INTEGER :: iz
LOGICAL :: ljout

INTEGER, OPTIONAL :: FastMocLv

END SUBROUTINE RayTraceGM_AFSS
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTrace_Tran_OMP_AFSS(RayInfo, CoreInfo, phis, phi1a, phi2a, PhiAngIn, xst, src, srcang1, srcang2, jout, iz, ljout, FastMocLv, lAFSS)

USE PARAM
USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type

IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo

REAL :: phis(:), PhiAngIn(:, :), xst(:), src(:), jout(:, :, :)
REAL :: phi1a(:,:,:)
REAL :: phi2a(:,:,:)
REAL :: srcang1(:,:,:), srcang2(:,:,:)
INTEGER :: iz
LOGICAL :: ljout
INTEGER, OPTIONAL :: FastMocLv
LOGICAL, OPTIONAL :: lAFSS

END SUBROUTINE RayTrace_Tran_OMP_AFSS
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceNM_OMP(RayInfo, CoreInfo, phisnm, PhiAngInnm, xstnm, srcnm, joutnm, iz, gb, ge, ljout)

USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type

IMPLICIT NONE

TYPE (RayInfo_Type) :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:,:)     :: phisnm, xstnm, srcnm
REAL, POINTER, DIMENSION(:,:,:)   :: PhiAngInnm
REAL, POINTER, DIMENSION(:,:,:,:) :: joutnm

INTEGER :: iz, gb, ge
LOGICAL :: ljout

END SUBROUTINE RayTraceNM_OMP
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceP1NM_OMP(RayInfo, CoreInfo, phisnm, phimnm, PhiAngInnm, xstnm, srcnm, srcmnm, joutnm, iz, gb, ge, ljout)

USE TYPEDEF, ONLY : RayInfo_Type, CoreInfo_type, Pin_Type, Cell_Type, MultigridInfo_Type

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

INTEGER :: iz, gb, ge
LOGICAL :: ljout

REAL, POINTER, DIMENSION(:,:)     :: phisnm, xstnm, srcnm
REAL, POINTER, DIMENSION(:,:,:)   :: phimnm, PhiAngInnm, srcmnm
REAL, POINTER, DIMENSION(:,:,:,:) :: joutnm

END SUBROUTINE RayTraceP1NM_OMP
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmp_OMP(RayInfo, CoreInfo, iz, iasy, gb, ge, ljout)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

INTEGER :: iz, iasy, gb, ge
LOGICAL :: ljout

END SUBROUTINE RayTraceDcmp_OMP
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceP1GM_OMP(RayInfo, CoreInfo, phis, phim, PhiAngIn, xst, src, srcm, jout, iz, ljout, ScatOd, FastMocLv)

USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type

IMPLICIT NONE

TYPE (RayInfo_Type) :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:)     :: phis, xst, src
REAL, POINTER, DIMENSION(:,:)   :: PhiAngIn, srcm, phim
REAL, POINTER, DIMENSION(:,:,:) :: jout

INTEGER :: iz
LOGICAL :: ljout
INTEGER :: ScatOd

INTEGER, OPTIONAL :: FastMocLv

END SUBROUTINE RayTraceP1GM_OMP
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceP1GM_AFSS(RayInfo, CoreInfo, phis, phim, PhiAngIn, xst, src, srcm, jout, iz, ljout, ScatOd, FastMocLv)

USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type

IMPLICIT NONE

TYPE (RayInfo_Type) :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:)     :: phis, xst, src
REAL, POINTER, DIMENSION(:,:)   :: PhiAngIn, srcm, phim
REAL, POINTER, DIMENSION(:,:,:) :: jout

INTEGER :: iz, ScatOd
LOGICAL :: ljout

INTEGER, OPTIONAL :: FastMocLv

END SUBROUTINE RayTraceP1GM_AFSS
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmp_Pn(RayInfo, CoreInfo, phisnm, phimnm, srcmnm, joutnm, iz, iAsy, gb, ge, ScatOd, lJout)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type

IMPLICIT NONE

TYPE (RayInfo_Type) :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:,:)     :: phisnm
REAL, POINTER, DIMENSION(:,:,:)   :: phimnm, srcmnm
REAL, POINTER, DIMENSION(:,:,:,:) :: joutnm

INTEGER :: iz, iAsy, gb, ge, ScatOd
LOGICAL :: ljout

END SUBROUTINE RayTraceDcmp_Pn
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayPn_Dcmp(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, ljout, iz, gb, ge, krot)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type, DcmpAsyRayInfo_Type

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

LOGICAL :: lJout
INTEGER :: iz, gb, ge, krot

END SUBROUTINE RecTrackRotRayPn_Dcmp
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayPn_Dcmp(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, ljout, iz, gb, ge, krot)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type, DcmpAsyRayInfo_Type

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

LOGICAL :: lJout
INTEGER :: iz, gb, ge, krot

END SUBROUTINE HexTrackRotRayPn_Dcmp
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceLS(RayInfo, CoreInfo, phis, PhiAngIn, xst, src, LinSrc, Slope1g, jout, iz, ljout)

USE PARAM
USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type

IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo

REAL, POINTER :: phis(:), PhiAngIn(:, :), xst(:), src(:), jout(:, :, :)
REAL, POINTER :: LinSrc(:, :)
REAL :: Slope1g(:, :)
INTEGER :: iz
LOGICAL :: ljout

END SUBROUTINE RayTraceLS
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceLS_CASMO(RayInfo, CoreInfo, phisnm, phisSlope, PhiAngInnm, srcnm, srcSlope, xstnm, joutnm, iz, gb, ge, ljout)

USE PARAM
USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type

IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo

REAL, POINTER :: phisnm(:, :), phisSlope(:, :, :, :), PhiAngInnm(:, :, :)
REAL, POINTER :: srcnm(:, :), srcSlope(:, :, :, :), xstnm(:, :), joutnm(:, :, :, :)
INTEGER :: iz, gb, ge
LOGICAL :: ljout

END SUBROUTINE RayTraceLS_CASMO
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmp_LSCASMO(RayInfo, CoreInfo, phisnm, phisSlope, joutnm, iz, iasy, gb, ge, ljout)

USE PARAM
USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type

IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo

REAL, POINTER :: phisnm(:, :)
REAL, POINTER :: srcnm(:, :), xstnm(:, :), joutnm(:, :, :, :)
REAL, POINTER :: phisSlope(:, :, :, :)
INTEGER :: iz, iasy, gb, ge
LOGICAL :: ljout

END SUBROUTINE RayTraceDcmp_LSCASMO
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE TrackRotRayLSDcmp_CASMO(RayInfo, CoreInfo, TrackingDat, phisC, phimx, phimy, jout, DcmpAsyRay, ljout, iz, gb, ge)

USE PARAM
USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, DcmpAsyRayInfo_Type, TrackingDat_Type

IMPLICIT NONE

TYPE(RayInfo_Type), INTENT(INOUT) :: RayInfo
TYPE(CoreInfo_Type), INTENT(INOUT) :: CoreInfo
TYPE(TrackingDat_Type), INTENT(INOUT) :: TrackingDat

TYPE(DcmpAsyRayInfo_Type) :: DcmpAsyRay
REAL, POINTER :: phisC(:, :), phimx(:, :, :), phimy(:, :, :), jout(:, :, :, :)
LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: iz, gb, ge

END SUBROUTINE TrackRotRayLSDcmp_CASMO
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtMacXsGM(core, Fxr, xstr, iz, ig, ng, lxslib, lTrCorrection, lRST, lssph, lssphreg, PE)

USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, PE_TYPE
USE BenchXs, ONLY : GetXstrBen

IMPLICIT NONE

TYPE (coreinfo_type) :: Core
TYPE (PE_type) :: PE

TYPE (Fxrinfo_type), DIMENSION(:) :: Fxr

REAL, POINTER, DIMENSION(:) :: xstr

INTEGER :: myzb, myze, ig, iz, ng
logical :: lxslib, lTrCorrection, lRST, lssph, lssphreg

END SUBROUTINE SetRtMacXsGM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtMacXsGM_Cusping(core, FmInfo, Fxr, xstr, phis, iz, ig, ng, lxslib, lTrCorrection, lRST, lssph, lssphreg, PE)

USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, PE_TYPE, FmInfo_Type
USE BenchXs, ONLY : GetXstrBen

IMPLICIT NONE

TYPE(coreinfo_type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(PE_type) :: PE

TYPE (Fxrinfo_type), DIMENSION(:) :: Fxr

REAL, POINTER, DIMENSION(:)     :: xstr
REAL, POINTER, DIMENSION(:,:,:) :: phis

INTEGER :: myzb, myze, ig, iz, ng
logical :: lxslib, lTrCorrection, lRST, lssph, lssphreg

END SUBROUTINE SetRtMacXsGM_Cusping
! ------------------------------------------------------------------------------------------------------------ 
SUBROUTINE SetRtMacXsNM(core, Fxr, xstnm, iz, ng, lxslib, lTrCorrection, lRST, lssph, lssphreg, PE)

USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, PE_TYPE

IMPLICIT NONE

TYPE(coreinfo_type) :: Core
TYPE(PE_TYPE) :: PE

TYPE (Fxrinfo_type), DIMENSION(:) :: Fxr

REAL, POINTER, DIMENSION(:,:) :: xstnm

INTEGER :: iz, ng
logical :: lxslib, lTrCorrection, lRST, lssph, lssphreg

END SUBROUTINE SetRtMacXsNM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtMacXsNM_Cusping(core, FmInfo, Fxr, xstnm, iz, ng, lxslib, lTrCorrection, lRST, lssph, lssphreg, PE)

USE PARAM
USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, PE_TYPE, FmInfo_Type

IMPLICIT NONE

TYPE(coreinfo_type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(Fxrinfo_type) :: Fxr(:)
TYPE(PE_TYPE) :: PE

REAL, POINTER :: xstnm(:, :)
INTEGER :: iz, ng
logical :: lxslib, lTrCorrection, lRST, lssph, lssphreg

END SUBROUTINE SetRtMacXsNM_Cusping
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtSrcGM(Core, Fxr, src, phis, psi, axsrc, xstr1g, eigv, iz, ig, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix, PE)

USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, GroupInfo_Type, PE_TYPE

IMPLICIT NONE

TYPE (CoreInfo_Type)  :: Core
TYPE (GroupInfo_Type) :: GroupInfo
TYPE(PE_TYPE)         :: PE

TYPE (Fxrinfo_type), DIMENSION(:) :: Fxr

REAL, POINTER, DIMENSION(:)     :: src, AxSrc, xstr1g
REAL, POINTER, DIMENSION(:,:)   :: psi
REAL, POINTER, DIMENSION(:,:,:) :: phis

REAL :: eigv
INTEGER :: myzb, myze, ig, ng, iz, ifsr, ifxr
LOGICAL :: lxslib, lscat1, l3dim, lNegFix

END SUBROUTINE SetRtSrcGM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtSrcGM_Cusping(Core, FmInfo, Fxr, src, phis, psi, axsrc, xstr1g, eigv, iz, ig, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix, PE)

USE PARAM
USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, GroupInfo_Type, PE_TYPE, FmInfo_Type

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(FxrInfo_Type) :: Fxr(:)
TYPE(GroupInfo_Type) :: GroupInfo

REAL, POINTER :: src(:), phis(:, :, :), psi(:, :), AxSrc(:), xstr1g(:)
REAL :: eigv
INTEGER :: myzb, myze, ig, ng, iz, ifsr, ifxr
LOGICAL :: lxslib, lscat1, l3dim
LOGICAL :: lNegFix
TYPE(PE_TYPE) :: PE

END SUBROUTINE SetRtSrcGM_Cusping
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtSrcNM(Core, Fxr, srcnm, phisnm, psi, AxSrc, xstnm, eigv, iz, gb, ge, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix, PE, Offset)

USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, GroupInfo_Type,  PE_TYPE

IMPLICIT NONE

TYPE (CoreInfo_Type)  :: Core
TYPE (GroupInfo_Type) :: GroupInfo
TYPE (PE_TYPE)        :: PE

TYPE (Fxrinfo_type), DIMENSION(:) :: Fxr

REAL, POINTER, DIMENSION(:,:)   :: srcnm, phisnm, psi, xstnm
REAL, POINTER, DIMENSION(:,:,:) :: AxSrc

REAL :: eigv
INTEGER :: gb, ge, ng, iz, ifsr, ifxr
LOGICAL :: lxslib, lscat1, l3dim, lNegFix

INTEGER, OPTIONAL :: Offset

END SUBROUTINE SetRtSrcNM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtSrcNM_Cusping(Core, FmInfo, Fxr, srcnm, phisnm, psi, AxSrc, xstnm, eigv, iz, gb, ge, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix, PE, Offset)

USE PARAM
USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, GroupInfo_Type, PE_TYPE, FmInfo_Type

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(FxrInfo_Type) :: Fxr(:)
TYPE(GroupInfo_Type):: GroupInfo
TYPE(PE_TYPE) :: PE

REAL, POINTER :: srcnm(:, :), phisnm(:, :), psi(:, :), AxSrc(:, :, :), xstnm(:, :)
REAL :: eigv
INTEGER :: gb, ge, ng, iz, ifsr, ifxr
LOGICAL :: lxslib, lscat1, l3dim
LOGICAL :: lNegFix
INTEGER, OPTIONAL :: Offset

END SUBROUTINE SetRtSrcNM_Cusping
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtLinSrc(Core, Fxr, RayInfo, src,LinSrc, LinPsi, LinSrcSlope, xstr1g, eigv, iz, ig, ng, GroupInfo, l3dim, lxslib, lscat1)

USE PARAM
USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, RayInfo_Type, GroupInfo_Type, XsMac_Type

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:)
TYPE(RayInfo_Type) :: RayInfo
TYPE(GroupInfo_Type):: GroupInfo

REAL, POINTER :: Src(:), Linsrc(:, :), LinPsi(:, :, :), LinSrcSlope(:, :, :, :), xstr1g(:)
REAL :: eigv
INTEGER :: myzb, myze, ig, ng, iz
LOGICAL :: lxslib, lscat1, l3dim

END SUBROUTINE SetRtLinSrc
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtLinSrc_CASMO(Core, Fxr, RayInfo, phisnm, phisSlope, srcnm, srcSlope, psi, psiSlope, axsrc, xstnm, eigv, iz, gb, ge, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix)

USE PARAM
USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, RayInfo_Type, GroupInfo_Type

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(RayInfo_Type) :: RayInfo
TYPE(GroupInfo_Type):: GroupInfo

REAL, POINTER :: phisnm(:, :), phisSlope(:, :, :, :)
REAL, POINTER :: srcnm(:, :), srcSlope(:, :, :, :), axsrc(:, :, :)
REAL, POINTER :: psi(:, :), psiSlope(:, :, :)
REAL, POINTER :: xstnm(:, :)
REAL :: eigv
INTEGER :: myzb, myze, gb, ge, ng, iz
LOGICAL :: lxslib, lscat1, l3dim, lNegFix

END SUBROUTINE SetRtLinSrc_CASMO
! ------------------------------------------------------------------------------------------------------------                             
SUBROUTINE SetRtP1SrcGM(Core, Fxr, srcm, phim, xstr1g, iz, ig, ng, GroupInfo, l3dim, lxslib, lscat1, ScatOd, PE)

USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, RayInfo_Type, GroupInfo_Type, XsMac_Type, PE_Type

IMPLICIT NONE

TYPE (CoreInfo_Type)  :: Core
TYPE (GroupInfo_Type) :: GroupInfo
TYPE (PE_Type)        :: PE

TYPE (Fxrinfo_type), DIMENSION(:) :: Fxr

REAL, POINTER, DIMENSION(:)       :: xstr1g
REAL, POINTER, DIMENSION(:,:)     :: srcm
REAL, POINTER, DIMENSION(:,:,:,:) :: phim

INTEGER :: myzb, myze, ig, ng, iz, ScatOd
LOGICAL :: lxslib, lscat1, l3dim

END SUBROUTINE SetRtP1SrcGM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtP1SrcNM(Core, Fxr, srcmnm, phimnm, xstnm, iz, gb, ge, ng, GroupInfo, lxslib, ScatOd, PE, Offset)

USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, GroupInfo_Type, PE_Type

IMPLICIT NONE

TYPE (CoreInfo_Type)  :: Core
TYPE (GroupInfo_Type) :: GroupInfo
TYPE (PE_Type)        :: PE

TYPE (Fxrinfo_type), DIMENSION(:) :: Fxr

REAL, POINTER, DIMENSION(:,:)   :: xstnm
REAL, POINTER, DIMENSION(:,:,:) :: srcmnm, phimnm

REAL :: eigv
INTEGER :: gb, ge, ng, iz, ifsr, ifxr, ScatOd
LOGICAL :: lxslib, lscat1, l3dim

INTEGER, OPTIONAL :: Offset

END SUBROUTINE SetRtP1SrcNM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE PseudoAbsorptionGM(Core, Fxr, src, phis1g, AxPXS, xstr1g, iz, ig, ng, GroupInfo, l3dim)

USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, GroupInfo_Type

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(GroupInfo_Type) :: GroupInfo

TYPE (Fxrinfo_type), DIMENSION(:) :: Fxr

REAL :: phis1g(:), AxPXS(:)
REAL, POINTER, DIMENSION(:) :: src, xstr1g
INTEGER :: ig, ng, iz
LOGICAL :: l3dim

END SUBROUTINE PseudoAbsorptionGM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE PseudoAbsorptionNM(Core, Fxr, AxPXS, xstnm, iz, ng, GroupInfo, l3dim)

USE TYPEDEF,        ONLY : coreinfo_type, Fxrinfo_type, Cell_Type, pin_Type, GroupInfo_Type, XsMac_Type
USE BenchXs,        ONLY : GetChiBen, xssben
USE MacXsLib_mod,   ONLY : MacXsScatMatrix
USE BasicOperation, ONLY : CP_CA

IMPLICIT NONE

TYPE(CoreInfo_Type)  :: Core
TYPE(GroupInfo_Type) :: GroupInfo

TYPE (Fxrinfo_type), DIMENSION(:) :: Fxr

REAL, POINTER, DIMENSION(:,:)   :: xstnm
REAL, POINTER, DIMENSION(:,:,:) :: AxPXS

REAL :: eigv
INTEGER :: iz, ng
LOGICAL :: l3dim

END SUBROUTINE PseudoAbsorptionNM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE PsiUpdate(Core, Fxr, phis, psi, myzb, myze, ng, lxslib, GroupInfo)

USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, GroupInfo_Type

IMPLICIT NONE

TYPE (coreinfo_type)  :: CORE
TYPE (GroupInfo_Type) :: GroupInfo

TYPE (Fxrinfo_type), DIMENSION(:,:) :: Fxr

REAL, POINTER, DIMENSION(:,:)   :: Psi
REAL, POINTER, DIMENSION(:,:,:) :: phis

INTEGER :: myzb, myze, ng
LOGICAL :: lXsLib

END SUBROUTINE PsiUpdate
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE PowerUpdate(Core, Fxr, phis, Power, myzb, myze, ng, lxslib, GroupInfo, PE)

USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, GroupInfo_Type, PE_TYPE

IMPLICIT NONE

TYPE (coreinfo_type)  :: CORE
TYPE (GroupInfo_Type) :: GroupInfo
TYPE (PE_TYPE)        :: PE

TYPE (Fxrinfo_type), DIMENSION(:,:) :: Fxr

REAL, POINTER, DIMENSION(:,:)   :: Power
REAL, POINTER, DIMENSION(:,:,:) :: phis

INTEGER :: myzb, myze, ng
LOGICAL :: lXsLib

END SUBROUTINE PowerUpdate
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE LinPsiUpdate(Core, Fxr, LinPsi, LinSrcSlope, myzb, myze, ng, lxslib, GroupInfo)

USE PARAM
USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, GroupInfo_Type

IMPLICIT NONE

TYPE(coreinfo_type) :: CORE
TYPE(Fxrinfo_type), POINTER :: Fxr(:, :)
TYPE(GroupInfo_Type) :: GroupInfo

REAL, POINTER :: LinPsi(:, :, :)
REAL, POINTER :: LinSrcSlope(:, :, :, :)
INTEGER :: myzb, myze, ng
LOGICAL :: lXsLib

END SUBROUTINE LinPsiUpdate
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE LinPsiUpdate_CASMO(Core, Fxr, phisSlope, psiSlope, myzb, myze, ng, lxslib, GroupInfo)

USE PARAM
USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, GroupInfo_Type

IMPLICIT NONE

TYPE(coreinfo_type) :: CORE
TYPE(Fxrinfo_type),POINTER :: Fxr(:, :)
TYPE(GroupInfo_Type) :: GroupInfo

REAL, POINTER :: phisSlope(:, :, :, :)
REAL, POINTER :: psiSlope(:, :, :)
INTEGER :: myzb, myze, ng
LOGICAL :: lXsLib

END SUBROUTINE LinPsiUpdate_CASMO
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE CellPsiUpdate(CORE, psi, psic, myzb, myze)

USE TYPEDEF, ONLY : coreinfo_type

IMPLICIT NONE

TYPE (coreinfo_type) :: CORE

REAL, POINTER, DIMENSION(:,:) :: Psi, psiC
INTEGER :: myzb, myze

END SUBROUTINE CellPsiUpdate
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE UpdateEigv(Core, psi, psid, eigv, peigv, myzb, myze, PE)

USE TYPEDEF, ONLY : coreinfo_type, PE_TYPE
USE BenchXs, ONLY : xsnfBen

USE BasicOperation, ONLY : CP_CA

IMPLICIT NONE

TYPE (coreinfo_type) :: CORE

REAL, POINTER, DIMENSION(:,:) :: Psi, psid

INTEGER :: myzb, myze, ng
REAL :: eigv, peigv

TYPE(PE_TYPE), OPTIONAL :: PE

END SUBROUTINE UpdateEigv
! ------------------------------------------------------------------------------------------------------------
Subroutine ApproxExp(PolarAng, npr)

USE TYPEDEF, ONLY : PolarAngle_Type

IMPLICIT NONE

TYPE (PolarAngle_Type) :: PolarAng(npr)

INTEGER :: npr

END SUBROUTINE ApproxExp
! ------------------------------------------------------------------------------------------------------------
Function PsiErr(Core, psi, psid, myzb, myze, PE)

USE TYPEDEF, ONLY : coreinfo_type, PE_TYPE

IMPLICIT NONE

TYPE (coreinfo_type) :: CORE
TYPE (PE_TYPE)       :: PE

REAL :: PsiErr
REAL, POINTER, DIMENSION(:,:) :: Psi, psid
INTEGER :: myzb, myze

END FUNCTION PsiErr
! ------------------------------------------------------------------------------------------------------------
FUNCTION MocResidual(Core, FmInfo,  eigv, GroupInfo, ng, PE, nTracerCntl)

USE TYPEDEF,  ONLY : CoreInfo_Type, FmInfo_Type, GroupInfo_Type
USE cntl,     ONLY : nTracerCntl_Type
USE PE_MOD,   ONLY : PE_TYPE

IMPLICIT NONE

TYPE (CoreInfo_Type)    :: Core
TYPE (FmInfo_Type)      :: FmInfo
TYPE (GroupInfo_Type)   :: GroupInfo
TYPE (PE_TYPE)          :: PE
TYPE (nTracerCntl_Type) :: nTracerCntl

REAL :: eigv, MocResidual
INTEGER :: ng

END FUNCTION MocResidual
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE AddBucklingGM(Core, Fxr, xstr1g, bsq, iz, ig, ng, lxslib, lRST)

USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type

IMPLICIT NONE

TYPE (CoreInfo_Type) :: Core

TYPE (Fxrinfo_type), POINTER, DIMENSION(:,:) :: Fxr

REAL :: xstr1g(:)
INTEGER :: iz, ig, ng
REAL :: Bsq
LOGICAL :: lxslib, lRST

END SUBROUTINE AddBucklingGM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE AddBucklingNM(Core, Fxr, xstnm, bsq, iz, ng, lxslib, lRST)

USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type

IMPLICIT NONE

TYPE (CoreInfo_Type) :: Core

TYPE (Fxrinfo_type), POINTER, DIMENSION(:,:) :: Fxr

REAL, POINTER, DIMENSION(:,:) :: xstnm

INTEGER :: iz, ng
REAL :: Bsq
LOGICAL :: lxslib, lRST

END SUBROUTINE AddBucklingNM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE AddConstSrc(Core, Fxr, Src, xstr1g, ConstSrc, iz, ig, ng)

USE TYPEDEF, ONLY : CoreInfo_Type, Fxrinfo_type, Cell_Type, pin_Type

IMPLICIT NONE

TYPE (CoreInfo_Type) :: Core

TYPE (FxrInfo_Type), DIMENSION(:) :: Fxr

REAL, POINTER, DIMENSION(:) :: Src, xstr1g

REAL :: ConstSrc
INTEGER :: iz, ig, ng

END SUBROUTINE AddConstSrc
! ------------------------------------------------------------------------------------------------------------
FUNCTION FxrAvgPhi(Core, Fxr, Phis, ipin, iLocalfxr, iz, ng, PE)

USE TYPEDEF, ONLY : CoreInfo_Type, PE_Type, FxrInfo_Type, Cell_Type, Pin_Type

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(PE_Type) :: PE

REAL, POINTER :: phis(:, :, :)
INTEGER :: iLocalfxr, ipin, iz, ng
REAL :: FxrAvgPhi(ng)

END FUNCTION FxrAvgPhi
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE DcmpLinkBoundaryFlux(CoreInfo, RayInfo, PhiAngInNM, DcmpPhiAngIn, DcmpPhiAngOut, gb, ge, color)

USE TYPEDEF, ONLY : CoreInfo_Type, RayInfo_Type, DcmpAsyRayInfo_Type

IMPLICIT NONE

TYPE (CoreInfo_Type) :: CoreInfo
TYPE (RayInfo_Type)  :: RayInfo

REAL, POINTER, DIMENSION(:,:,:)     :: PhiAngInNM
REAL, POINTER, DIMENSION(:,:,:,:,:) :: DcmpPhiAngIn, DcmpPhiAngOut

INTEGER :: gb, ge
INTEGER, OPTIONAL :: color

END SUBROUTINE DcmpLinkBoundaryFlux
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE DcmpLinkBoundaryFlux1g(CoreInfo, RayInfo, PhiAngIn, DcmpPhiAngIn1g, DcmpPhiAngOut1g, color)

USE TYPEDEF, ONLY : CoreInfo_Type, RayInfo_Type, DcmpAsyRayInfo_Type

IMPLICIT NONE

TYPE (CoreInfo_Type) :: CoreInfo
TYPE (RayInfo_Type)  :: RayInfo

REAL, POINTER, DIMENSION(:,:)     :: PhiAngIn
REAL, POINTER, DIMENSION(:,:,:,:) :: DcmpPhiAngIn1g, DcmpPhiAngOut1g

INTEGER, OPTIONAL :: color

END SUBROUTINE DcmpLinkBoundaryFlux1g
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE DcmpScatterBoundaryFlux(RayInfo, PhiAngInNM, DcmpPhiAngIn)

USE TYPEDEF, ONLY : RayInfo_Type

IMPLICIT NONE

TYPE (RayInfo_Type) :: RayInfo

REAL, POINTER, DIMENSION(:,:,:)     :: PhiAngInNM
REAL, POINTER, DIMENSION(:,:,:,:,:) :: DcmpPhiAngIn

END SUBROUTINE DcmpScatterBoundaryFlux
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE DcmpScatterBoundaryFlux1g(RayInfo, PhiAngIn, DcmpPhiAngIn1g)

USE TYPEDEF, ONLY : RayInfo_Type

IMPLICIT NONE

TYPE (RayInfo_Type) :: RayInfo

REAL, POINTER, DIMENSION(:,:)     :: PhiAngIn
REAL, POINTER, DIMENSION(:,:,:,:) :: DcmpPhiAngIn1g

END SUBROUTINE DcmpScatterBoundaryFlux1g
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE DcmpGatherCurrent(CoreInfo, joutNM)

USE TYPEDEF, ONLY : CoreInfo_Type

IMPLICIT NONE

TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:,:,:,:) :: joutNM

END SUBROUTINE DcmpGatherCurrent
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE DcmpGatherCurrent1g(CoreInfo, jout)

USE TYPEDEF, ONLY : CoreInfo_Type

IMPLICIT NONE

TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:,:,:) :: jout

END SUBROUTINE DcmpGatherCurrent1g
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE DcmpGatherBoundaryFlux(RayInfo, DcmpPhiAngOut)

USE TYPEDEF, ONLY : RayInfo_Type

IMPLICIT NONE

TYPE (RayInfo_Type) :: RayInfo

REAL, POINTER, DIMENSION(:,:,:,:,:)  :: DcmpPhiAngOut

END SUBROUTINE DcmpGatherBoundaryFlux
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE DcmpGatherBoundaryFlux1g(RayInfo, DcmpPhiAngOut1g)

USE TYPEDEF, ONLY : RayInfo_Type

IMPLICIT NONE

TYPE (RayInfo_Type) :: RayInfo

REAL, POINTER, DIMENSION(:,:,:,:)  :: DcmpPhiAngOut1g

END SUBROUTINE DcmpGatherBoundaryFlux1g
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE FluxUnderRelaxation(Core, Phis1g, Phis, w, iz, ig, PE)

USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type, PE_TYPE, Pin_Type, Cell_Type

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core

REAL :: Phis1g(:)
REAL :: Phis(:, :, :)
REAL :: w
INTEGER :: iz, ig
TYPE(PE_TYPE) :: PE

END SUBROUTINE FluxUnderRelaxation
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE FluxInUnderRelaxation(Core, PhiAngIn1g, PhiAngIn, w, n1, n2, iz, ig, PE)

USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type, PE_TYPE

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core

REAL :: PhiAngIn1g(:, :)
REAL :: PhiAngIn(:, :, :, :)
INTEGER :: n1, n2, iz, ig
REAL :: w
TYPE(PE_TYPE) :: PE

END SUBROUTINE FluxInUnderRelaxation
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE CurrentUnderRelaxation(Core, Jout1g, Jout, w, iz, ig, PE)

USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type, PE_TYPE

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core

REAL :: Jout1g(:, :, :)
REAL :: Jout(:, :, :, :, :)
REAL :: w
INTEGER :: iz, ig
TYPE(PE_TYPE) :: PE

END SUBROUTINE CurrentUnderRelaxation
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE GetNeighborMocFlux(phis, neighphis, nFsr, myzb, myze, gb, ge, nz, AxBC)

IMPLICIT NONE

REAL, POINTER :: phis(:, :, :), neighphis(:, :, :)
INTEGER :: nfsr, myzb, myze, gb, ge, nz
INTEGER :: AxBC(2)

END SUBROUTINE GetNeighborMocFlux
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE initRT(RayInfo, CoreInfo, nTracerCntl, PE)

USE TYPEDEF, ONLY : RayInfo_Type, CoreInfo_type, PE_TYPE
USE Cntl,    ONLY : nTracerCntl_Type

IMPLICIT NONE

TYPE (RayInfo_Type)     :: RayInfo
TYPE (CoreInfo_Type)    :: CoreInfo
TYPE (nTracerCntl_Type) :: nTracerCntl
TYPE (PE_TYPE)          :: PE

END SUBROUTINE initRT
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE initAFSS(RayInfo, CoreInfo, nTracerCntl, PE)

USE TYPEDEF, ONLY : RayInfo_Type, CoreInfo_type, PE_TYPE
USE Cntl,    ONLY : nTracerCntl_Type

IMPLICIT NONE

TYPE (RayInfo_Type)     :: RayInfo
TYPE (CoreInfo_Type)    :: CoreInfo
TYPE (nTracerCntl_Type) :: nTracerCntl
TYPE (PE_TYPE)          :: PE

END SUBROUTINE initAFSS
! ------------------------------------------------------------------------------------------------------------
END INTERFACE

END MODULE MOC_MOD