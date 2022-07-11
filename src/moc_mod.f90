MODULE MOC_MOD
  
USE TYPEDEF, ONLY : TrackingDat_Type

IMPLICIT NONE

! Basic
REAL, POINTER, DIMENSION(:,:)   :: wtang, hwt
REAL, POINTER, DIMENSION(:,:,:) :: wtsurf, comp, mwt, mwt2, LinPsi

INTEGER :: nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay

INTEGER, TARGET :: AziMap(360, 2)

TYPE(TrackingDat_Type), SAVE :: TrackingDat(100)

REAL, TARGET, DIMENSION(-40000:0, 1:12) :: expa, expb ! Approximation of Exponetial Function

! Group Major
REAL, POINTER, DIMENSION(:)     :: phis1g, xst1g, src1g, axSrc1g, axPxs1g
REAL, POINTER, DIMENSION(:,:)   :: phiAngIn1g, phim1g, LinSrc1g, srcm1g
REAL, POINTER, DIMENSION(:,:,:) :: MocJout1g, SrcAng1g1, SrcAng1g2

! Node Major
REAL, POINTER, DIMENSION(:,:)     :: phisNg, srcNg, xstNg
REAL, POINTER, DIMENSION(:,:,:)   :: phimNg, PhiAngInNg, srcmNg
REAL, POINTER, DIMENSION(:,:,:,:) :: MocJoutNg, SrcAngNg1, SrcAngNg2

! Domain Decomposition
INTEGER :: nClr ! # of Colors

REAL, POINTER, DIMENSION(:,:,:,:,:) :: DcmpPhiAngInNg, DcmpPhiAngOutNg
REAL, POINTER, DIMENSION(:,:,:,:)   :: DcmpPhiAngIn1g, DcmpPhiAngOut1g

INTEGER, POINTER, DIMENSION(:,:)   :: DcmpAsyClr ! (iAsy, iClr)
INTEGER, POINTER, DIMENSION(:,:,:) :: DcmpAziRay ! (imRay, iAzi, iAsy)

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
SUBROUTINE RayTrace_Tran(RayInfo, CoreInfo, phis, phi1a, phi2a, PhiAngIn, xst, src, srcang1g1, srcang1g2, jout, iz, ljout, FastMocLv, lAFSS)

USE PARAM
USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type

IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL :: phis(:), PhiAngIn(:, :), xst(:), src(:), jout(:, :, :)
REAL :: phi1a(:,:,:), phi2a(:,:,:)
REAL :: srcang1g1(:,:,:), srcang1g2(:,:,:)

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
SUBROUTINE RtDcmpThr_GM(RayInfo, CoreInfo, TrackingLoc, phis, MocJout, jAsy, iz, lJout, lHex, lAFSS)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type

IMPLICIT NONE

TYPE (RayInfo_Type)     :: RayInfo
TYPE (CoreInfo_Type)    :: CoreInfo
TYPE (TrackingDat_Type) :: TrackingLoc

REAL, POINTER, DIMENSION(:)     :: phis
REAL, POINTER, DIMENSION(:,:,:) :: MocJout

INTEGER :: jAsy, iz
LOGICAL :: lJout, lHex, lAFSS

END SUBROUTINE RtDcmpThr_GM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayDcmp_GM(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, ljout, iz, krot, lAFSS)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type, DcmpAsyRayInfo_Type

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

LOGICAL :: ljout, lAFSS
INTEGER :: iz, krot

END SUBROUTINE HexTrackRotRayDcmp_GM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmpP1_GM(RayInfo, CoreInfo, phis1g, phim1g, PhiAngIn1g, xst1g, src1g, srcm1g, MocJout1g, iz, lJout)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:)     :: phis1g, xst1g, src1g
REAL, POINTER, DIMENSION(:,:)   :: PhiAngIn1g, phim1g, srcm1g
REAL, POINTER, DIMENSION(:,:,:) :: MocJout1g

INTEGER :: iz
LOGICAL :: lJout

END SUBROUTINE RayTraceDcmpP1_GM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RtDcmpP1Thr_GM(RayInfo, CoreInfo, TrackingLoc, phis1g, phim1g, srcm1g, MocJout1g, jAsy, iz, lJout, lHex, lAFSS, ScatOd)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type

IMPLICIT NONE

TYPE (RayInfo_Type)     :: RayInfo
TYPE (CoreInfo_Type)    :: CoreInfo
TYPE (TrackingDat_Type) :: TrackingLoc

REAL, POINTER, DIMENSION(:)     :: phis1g
REAL, POINTER, DIMENSION(:,:)   :: phim1g, srcm1g
REAL, POINTER, DIMENSION(:,:,:) :: MocJout1g

INTEGER :: jAsy, iz, ScatOd
LOGICAL :: lJout, lHex, lAFSS

END SUBROUTINE RtDcmpP1Thr_GM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayDcmpP1_GM(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, ljout, iz, krot, lAFSS, ScatOd)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type, DcmpAsyRayInfo_Type

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

LOGICAL :: ljout, lAFSS
INTEGER :: iz, krot, ScatOd

END SUBROUTINE HexTrackRotRayDcmpP1_GM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmp_NM(RayInfo, CoreInfo, phisNg, PhiAngInNg, xstNg, srcNg, MocJoutNg, iz, gb, ge, lJout)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:,:)     :: phisNg, xstNg, srcNg
REAL, POINTER, DIMENSION(:,:,:)   :: PhiAngInNg
REAL, POINTER, DIMENSION(:,:,:,:) :: MocJoutNg

INTEGER :: iz, gb, ge
LOGICAL :: lJout

END SUBROUTINE RayTraceDcmp_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RtDcmpThr_NM(RayInfo, CoreInfo, TrackingLoc, phisNg, MocJoutNg, jAsy, iz, gb, ge, lJout, lHex, lAFSS)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type

IMPLICIT NONE

TYPE (RayInfo_Type)     :: RayInfo
TYPE (CoreInfo_Type)    :: CoreInfo
TYPE (TrackingDat_Type) :: TrackingLoc

REAL, POINTER, DIMENSION(:,:)     :: phisNg
REAL, POINTER, DIMENSION(:,:,:,:) :: MocJoutNg

INTEGER :: jAsy, iz, gb, ge
LOGICAL :: lJout, lHex, lAFSS

END SUBROUTINE RtDcmpThr_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayDcmp_NM(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, ljout, iz, gb, ge, krot, lAFSS)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type, DcmpAsyRayInfo_Type

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

LOGICAL :: ljout, lAFSS
INTEGER :: iz, gb, ge, krot

END SUBROUTINE RecTrackRotRayDcmp_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayDcmp_NM(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, ljout, iz, gb, ge, krot, lAFSS)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type, DcmpAsyRayInfo_Type

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

LOGICAL :: ljout, lAFSS
INTEGER :: iz, gb, ge, krot

END SUBROUTINE HexTrackRotRayDcmp_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmpP1_NM(RayInfo, CoreInfo, phisNg, phimNg, PhiAngInNg, xstNg, srcNg, srcmNg, MocJoutNg, iz, gb, ge, lJout)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:,:)     :: phisNg, xstNg, srcNg
REAL, POINTER, DIMENSION(:,:,:)   :: phimNg, PhiAngInNg, srcmNg
REAL, POINTER, DIMENSION(:,:,:,:) :: MocJoutNg

INTEGER :: iz, gb, ge
LOGICAL :: lJout

END SUBROUTINE RayTraceDcmpP1_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RtDcmpP1Thr_NM(RayInfo, CoreInfo, TrackingLoc, phisNg, phimNg, srcmNg, JoutNg, jAsy, iz, gb, ge, lJout, lHex, lAFSS, ScatOd)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type

IMPLICIT NONE

TYPE (RayInfo_Type)     :: RayInfo
TYPE (CoreInfo_Type)    :: CoreInfo
TYPE (TrackingDat_Type) :: TrackingLoc

REAL, POINTER, DIMENSION(:,:)     :: phisNg
REAL, POINTER, DIMENSION(:,:,:)   :: phimNg, srcmNg
REAL, POINTER, DIMENSION(:,:,:,:) :: JoutNg

INTEGER :: jAsy, iz, gb, ge, ScatOd
LOGICAL :: lJout, lHex, lAFSS

END SUBROUTINE RtDcmpP1Thr_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayDcmpP1_NM(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, lJout, iz, gb, ge, krot, lAFSS, ScatOd)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type, DcmpAsyRayInfo_Type

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

LOGICAL :: lJout, lAFSS
INTEGER :: iz, gb, ge, krot, ScatOd

END SUBROUTINE HexTrackRotRayDcmpP1_NM
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
SUBROUTINE RayTrace_GM(RayInfo, CoreInfo, phis, PhiAngIn, xst, src, jout, iz, ljout, FastMocLv)

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

END SUBROUTINE RayTrace_GM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRay_GM(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, krot, FastMocLv)

USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type, TrackingDat_Type

IMPLICIT NONE

TYPE(RayInfo_Type)     :: RayInfo
TYPE(CoreInfo_Type)    :: CoreInfo
TYPE(TrackingDat_Type) :: TrackingDat

LOGICAL :: ljout
INTEGER :: irotray, iz, krot, FastMocLv

END SUBROUTINE RecTrackRotRay_GM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRay_GM(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, krot)

USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type, TrackingDat_Type

IMPLICIT NONE

TYPE(RayInfo_Type)     :: RayInfo
TYPE(CoreInfo_Type)    :: CoreInfo
TYPE(TrackingDat_Type) :: TrackingDat

LOGICAL :: ljout
INTEGER :: irotray, iz, krot

END SUBROUTINE HexTrackRotRay_GM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTrace_NM(RayInfo, CoreInfo, phisNg, PhiAngInNg, xstNg, srcNg, JoutNg, iz, gb, ge, ljout)

USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type

IMPLICIT NONE

TYPE (RayInfo_Type) :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:,:)     :: phisNg, xstNg, srcNg
REAL, POINTER, DIMENSION(:,:,:)   :: PhiAngInNg
REAL, POINTER, DIMENSION(:,:,:,:) :: JoutNg

INTEGER :: iz, gb, ge
LOGICAL :: ljout

END SUBROUTINE RayTrace_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRay_NM(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, krot, gb, ge)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type

IMPLICIT NONE

TYPE(RayInfo_Type)     :: RayInfo
TYPE(CoreInfo_Type)    :: CoreInfo
TYPE(TrackingDat_Type) :: TrackingDat

LOGICAL :: ljout
INTEGER :: irotray, iz, krot, gb, ge

END SUBROUTINE RecTrackRotRay_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRay_NM(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, krot, gb, ge)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type

IMPLICIT NONE

TYPE(RayInfo_Type)     :: RayInfo
TYPE(CoreInfo_Type)    :: CoreInfo
TYPE(TrackingDat_Type) :: TrackingDat

LOGICAL :: ljout
INTEGER :: irotray, iz, krot, gb, ge

END SUBROUTINE HexTrackRotRay_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceP1_GM(RayInfo, CoreInfo, phis1g, phim1g, PhiAngIn1g, xst1g, src1g, srcm1g, jout1g, iz, ljout, FastMocLv)

USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:)     :: phis1g, xst1g, src1g
REAL, POINTER, DIMENSION(:,:)   :: PhiAngIn1g, srcm1g, phim1g
REAL, POINTER, DIMENSION(:,:,:) :: jout1g

INTEGER :: iz
LOGICAL :: ljout

INTEGER, OPTIONAL :: FastMocLv

END SUBROUTINE RayTraceP1_GM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayP1_GM(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, krot, ScatOd, FastMocLv)

USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type, TrackingDat_Type

IMPLICIT NONE

TYPE (RayInfo_Type)     :: RayInfo
TYPE (CoreInfo_Type)    :: CoreInfo
TYPE (TrackingDat_Type) :: TrackingDat

LOGICAL :: ljout
INTEGER :: irotray, iz, krot, ScatOd, FastMocLv

END SUBROUTINE RecTrackRotRayP1_GM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayP1_GM(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, krot, ScatOd, FastMocLv)

USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type, TrackingDat_Type

IMPLICIT NONE

TYPE (RayInfo_Type)     :: RayInfo
TYPE (CoreInfo_Type)    :: CoreInfo
TYPE (TrackingDat_Type) :: TrackingDat

LOGICAL :: ljout
INTEGER :: irotray, iz, krot, ScatOd, FastMocLv

END SUBROUTINE HexTrackRotRayP1_GM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceP1_NM(RayInfo, CoreInfo, phisNg, phimNg, PhiAngInNg, xstNg, srcNg, srcmNg, JoutNg, iz, gb, ge, ljout)

USE TYPEDEF, ONLY : RayInfo_Type, CoreInfo_type, Pin_Type, Cell_Type

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

INTEGER :: iz, gb, ge
LOGICAL :: ljout

REAL, POINTER, DIMENSION(:,:)     :: phisNg, xstNg, srcNg
REAL, POINTER, DIMENSION(:,:,:)   :: phimNg, PhiAngInNg, srcmNg
REAL, POINTER, DIMENSION(:,:,:,:) :: JoutNg

END SUBROUTINE RayTraceP1_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayP1_NM(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, krot, gb, ge, ScatOd)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type

IMPLICIT NONE

TYPE(RayInfo_Type)     :: RayInfo
TYPE(CoreInfo_Type)    :: CoreInfo
TYPE(TrackingDat_Type) :: TrackingDat

LOGICAL :: ljout
INTEGER :: irotray, iz, krot, gb, ge, ScatOd

END SUBROUTINE RecTrackRotRayP1_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayP1_NM(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, krot, gb, ge, ScatOd)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type

IMPLICIT NONE

TYPE(RayInfo_Type)     :: RayInfo
TYPE(CoreInfo_Type)    :: CoreInfo
TYPE(TrackingDat_Type) :: TrackingDat

LOGICAL :: ljout
INTEGER :: irotray, iz, krot, gb, ge, ScatOd

END SUBROUTINE HexTrackRotRayP1_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceLS(RayInfo, CoreInfo, phis, PhiAngIn, xst, src1g, LinSrc, Slope1g, jout, iz, ljout)

USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type

IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo

REAL, POINTER :: phis(:), PhiAngIn(:, :), xst(:), src1g(:), jout(:, :, :)
REAL, POINTER :: LinSrc(:, :)
REAL :: Slope1g(:, :)
INTEGER :: iz
LOGICAL :: ljout

END SUBROUTINE RayTraceLS
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceLS_CASMO(RayInfo, CoreInfo, phisNg, phisSlope, PhiAngInNg, srcNg, srcSlope, xstNg, JoutNg, iz, gb, ge, ljout)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type

IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo

REAL, POINTER :: phisNg(:, :), phisSlope(:, :, :, :), PhiAngInNg(:, :, :)
REAL, POINTER :: srcNg(:, :), srcSlope(:, :, :, :), xstNg(:, :), JoutNg(:, :, :, :)
INTEGER :: iz, gb, ge
LOGICAL :: ljout

END SUBROUTINE RayTraceLS_CASMO
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmp_LSCASMO(RayInfo, CoreInfo, phisNg, phisSlope, JoutNg, iz, iasy, gb, ge, ljout)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type

IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo

REAL, POINTER :: phisNg(:, :)
REAL, POINTER :: xstNg(:, :), JoutNg(:, :, :, :)
REAL, POINTER :: phisSlope(:, :, :, :)
INTEGER :: iz, iasy, gb, ge
LOGICAL :: ljout

END SUBROUTINE RayTraceDcmp_LSCASMO
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE TrackRotRayLSDcmp_CASMO(RayInfo, CoreInfo, TrackingDat, phisC, phimx, phimy, jout, DcmpAsyRay, ljout, iz, gb, ge)

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
SUBROUTINE SetRtMacXsNM(core, Fxr, xstNg, iz, ng, lxslib, lTrCorrection, lRST, lssph, lssphreg, PE)

USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, PE_TYPE

IMPLICIT NONE

TYPE(coreinfo_type) :: Core
TYPE(PE_TYPE) :: PE

TYPE (Fxrinfo_type), DIMENSION(:) :: Fxr

REAL, POINTER, DIMENSION(:,:) :: xstNg

INTEGER :: iz, ng
logical :: lxslib, lTrCorrection, lRST, lssph, lssphreg

END SUBROUTINE SetRtMacXsNM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtMacXsNM_Cusping(core, FmInfo, Fxr, xstNg, iz, ng, lxslib, lTrCorrection, lRST, lssph, lssphreg, PE)

USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, PE_TYPE, FmInfo_Type

IMPLICIT NONE

TYPE(coreinfo_type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(Fxrinfo_type) :: Fxr(:)
TYPE(PE_TYPE) :: PE

REAL, POINTER :: xstNg(:, :)
INTEGER :: iz, ng
logical :: lxslib, lTrCorrection, lRST, lssph, lssphreg

END SUBROUTINE SetRtMacXsNM_Cusping
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtSrcGM(Core, Fxr, src1g, phis, psi, axsrc, xstr1g, eigv, iz, ig, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix, PE)

USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, GroupInfo_Type, PE_TYPE

IMPLICIT NONE

TYPE (CoreInfo_Type)  :: Core
TYPE (GroupInfo_Type) :: GroupInfo
TYPE(PE_TYPE)         :: PE

TYPE (Fxrinfo_type), DIMENSION(:) :: Fxr

REAL, POINTER, DIMENSION(:)     :: src1g, AxSrc, xstr1g
REAL, POINTER, DIMENSION(:,:)   :: psi
REAL, POINTER, DIMENSION(:,:,:) :: phis

REAL :: eigv
INTEGER :: myzb, myze, ig, ng, iz, ifsr, ifxr
LOGICAL :: lxslib, lscat1, l3dim, lNegFix

END SUBROUTINE SetRtSrcGM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtSrcGM_Cusping(Core, FmInfo, Fxr, src, phis, psi, axsrc, xstr1g, eigv, iz, ig, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix, PE)

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
SUBROUTINE SetRtSrcNM(Core, Fxr, srcNg, phisNg, psi, AxSrc, xstNg, eigv, iz, gb, ge, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix, PE, Offset)

USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, GroupInfo_Type, PE_TYPE

IMPLICIT NONE

TYPE (CoreInfo_Type)  :: Core
TYPE (GroupInfo_Type) :: GroupInfo
TYPE (PE_TYPE)        :: PE

TYPE (Fxrinfo_type), DIMENSION(:) :: Fxr

REAL, POINTER, DIMENSION(:,:)   :: srcNg, psi, xstNg, phisNg
REAL, POINTER, DIMENSION(:,:,:) :: AxSrc

REAL :: eigv
INTEGER :: gb, ge, ng, iz, ifsr, ifxr
LOGICAL :: lxslib, lscat1, l3dim, lNegFix

INTEGER, OPTIONAL :: Offset

END SUBROUTINE SetRtSrcNM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtSrcNM_Cusping(Core, FmInfo, Fxr, srcNg, phisNg, psi, AxSrc, xstNg, eigv, iz, gb, ge, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix, PE, Offset)

USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, GroupInfo_Type, PE_TYPE, FmInfo_Type

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(FxrInfo_Type) :: Fxr(:)
TYPE(GroupInfo_Type):: GroupInfo
TYPE(PE_TYPE) :: PE

REAL, POINTER :: srcNg(:, :), phisNg(:, :), psi(:, :), AxSrc(:, :, :), xstNg(:, :)
REAL :: eigv
INTEGER :: gb, ge, ng, iz, ifsr, ifxr
LOGICAL :: lxslib, lscat1, l3dim
LOGICAL :: lNegFix
INTEGER, OPTIONAL :: Offset

END SUBROUTINE SetRtSrcNM_Cusping
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtLinSrc(Core, Fxr, RayInfo, src1g, LinSrc, LinPsi, LinSrcSlope, xstr1g, eigv, iz, ig, ng, GroupInfo, l3dim, lxslib, lscat1)

USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, RayInfo_Type, GroupInfo_Type

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:)
TYPE(RayInfo_Type) :: RayInfo
TYPE(GroupInfo_Type):: GroupInfo

REAL, POINTER :: Src1g(:), Linsrc(:, :), LinPsi(:, :, :), LinSrcSlope(:, :, :, :), xstr1g(:)
REAL :: eigv
INTEGER :: myzb, myze, ig, ng, iz
LOGICAL :: lxslib, lscat1, l3dim

END SUBROUTINE SetRtLinSrc
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtLinSrc_CASMO(Core, Fxr, RayInfo, phisNg, phisSlope, srcNg, srcSlope, psi, psiSlope, axsrc, xstNg, eigv, iz, gb, ge, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix)

USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, RayInfo_Type, GroupInfo_Type

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(RayInfo_Type) :: RayInfo
TYPE(GroupInfo_Type):: GroupInfo

REAL, POINTER :: phisNg(:, :), phisSlope(:, :, :, :)
REAL, POINTER :: srcNg(:, :), srcSlope(:, :, :, :), axsrc(:, :, :)
REAL, POINTER :: psi(:, :), psiSlope(:, :, :)
REAL, POINTER :: xstNg(:, :)
REAL :: eigv
INTEGER :: myzb, myze, gb, ge, ng, iz
LOGICAL :: lxslib, lscat1, l3dim, lNegFix

END SUBROUTINE SetRtLinSrc_CASMO
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtP1SrcGM(Core, Fxr, srcm1g, phimNg, xst1g, iz, ig, ng, GroupInfo, l3dim, lxslib, lscat1, ScatOd, PE)

USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, GroupInfo_Type, PE_Type

IMPLICIT NONE

TYPE (CoreInfo_Type)  :: Core
TYPE (GroupInfo_Type) :: GroupInfo
TYPE (PE_Type)        :: PE

TYPE (Fxrinfo_type), DIMENSION(:) :: Fxr

REAL, POINTER, DIMENSION(:)     :: xst1g
REAL, POINTER, DIMENSION(:,:)   :: srcm1g
REAL, POINTER, DIMENSION(:,:,:) :: phimNg

INTEGER :: myzb, myze, ig, ng, iz, ScatOd
LOGICAL :: lxslib, lscat1, l3dim

END SUBROUTINE SetRtP1SrcGM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtP1SrcNM(Core, Fxr, srcmNg, phimNg, xstNg, iz, gb, ge, ng, GroupInfo, lxslib, ScatOd, PE, Offset)

USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, GroupInfo_Type, PE_Type

IMPLICIT NONE

TYPE (CoreInfo_Type)  :: Core
TYPE (GroupInfo_Type) :: GroupInfo
TYPE (PE_Type)        :: PE

TYPE (Fxrinfo_type), DIMENSION(:) :: Fxr

REAL, POINTER, DIMENSION(:,:)   :: xstNg
REAL, POINTER, DIMENSION(:,:,:) :: srcmNg, phimNg

REAL :: eigv
INTEGER :: gb, ge, ng, iz, ifsr, ifxr, ScatOd
LOGICAL :: lxslib, lscat1, l3dim

INTEGER, OPTIONAL :: Offset

END SUBROUTINE SetRtP1SrcNM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE PseudoAbsorptionGM(Core, Fxr, phis1g, AxPXS, xstr1g, iz, ig, ng, GroupInfo, l3dim)

USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, GroupInfo_Type

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(GroupInfo_Type) :: GroupInfo

TYPE (Fxrinfo_type), DIMENSION(:) :: Fxr

REAL :: phis1g(:), AxPXS(:)
REAL, POINTER, DIMENSION(:) :: xstr1g
INTEGER :: ig, ng, iz
LOGICAL :: l3dim

END SUBROUTINE PseudoAbsorptionGM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE PseudoAbsorptionNM(Core, Fxr, AxPXS, xstNg, iz, ng, GroupInfo, l3dim)

USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, GroupInfo_Type

IMPLICIT NONE

TYPE(CoreInfo_Type)  :: Core
TYPE(GroupInfo_Type) :: GroupInfo

TYPE (Fxrinfo_type), DIMENSION(:) :: Fxr

REAL, POINTER, DIMENSION(:,:)   :: xstNg
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
SUBROUTINE AddBucklingNM(Core, Fxr, xstNg, bsq, iz, ng, lxslib, lRST)

USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type

IMPLICIT NONE

TYPE (CoreInfo_Type) :: Core

TYPE (Fxrinfo_type), POINTER, DIMENSION(:,:) :: Fxr

REAL, POINTER, DIMENSION(:,:) :: xstNg

INTEGER :: iz, ng
REAL :: Bsq
LOGICAL :: lxslib, lRST

END SUBROUTINE AddBucklingNM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE AddConstSrc(Core, Fxr, Src, xstr1g, ConstSrc, iz, ig, ng)

USE TYPEDEF, ONLY : CoreInfo_Type, Fxrinfo_type

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
SUBROUTINE DcmpLinkBndyFluxNg(CoreInfo, RayInfo, PhiAngInNg, DcmpPhiAngInNg, DcmpPhiAngOutNg, gb, ge, iClr)

USE TYPEDEF, ONLY : CoreInfo_Type, RayInfo_Type, DcmpAsyRayInfo_Type

IMPLICIT NONE

TYPE (CoreInfo_Type) :: CoreInfo
TYPE (RayInfo_Type)  :: RayInfo

REAL, POINTER, DIMENSION(:,:,:)     :: PhiAngInNg
REAL, POINTER, DIMENSION(:,:,:,:,:) :: DcmpPhiAngInNg, DcmpPhiAngOutNg

INTEGER :: gb, ge
INTEGER, OPTIONAL :: iClr

END SUBROUTINE DcmpLinkBndyFluxNg
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE DcmpLinkBndyFlux1g(CoreInfo, RayInfo, PhiAngIn1g, DcmpPhiAngIn1g, DcmpPhiAngOut1g, iClr)

USE TYPEDEF, ONLY : CoreInfo_Type, RayInfo_Type, DcmpAsyRayInfo_Type

IMPLICIT NONE

TYPE (CoreInfo_Type) :: CoreInfo
TYPE (RayInfo_Type)  :: RayInfo

REAL, POINTER, DIMENSION(:,:)     :: PhiAngIn1g
REAL, POINTER, DIMENSION(:,:,:,:) :: DcmpPhiAngIn1g, DcmpPhiAngOut1g

INTEGER, OPTIONAL :: iClr

END SUBROUTINE DcmpLinkBndyFlux1g
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE DcmpScatterBndyFluxNg(RayInfo, PhiAngInNg, DcmpPhiAngInNg)

USE TYPEDEF, ONLY : RayInfo_Type

IMPLICIT NONE

TYPE (RayInfo_Type) :: RayInfo

REAL, POINTER, DIMENSION(:,:,:)     :: PhiAngInNg
REAL, POINTER, DIMENSION(:,:,:,:,:) :: DcmpPhiAngInNg

END SUBROUTINE DcmpScatterBndyFluxNg
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE DcmpScatterBndyFlux1g(RayInfo, PhiAngIn1g, DcmpPhiAngIn1g)

USE TYPEDEF, ONLY : RayInfo_Type

IMPLICIT NONE

TYPE (RayInfo_Type) :: RayInfo

REAL, POINTER, DIMENSION(:,:)     :: PhiAngIn1g
REAL, POINTER, DIMENSION(:,:,:,:) :: DcmpPhiAngIn1g

END SUBROUTINE DcmpScatterBndyFlux1g
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE DcmpGatherCurrentNg(CoreInfo, JoutNg)

USE TYPEDEF, ONLY : CoreInfo_Type

IMPLICIT NONE

TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:,:,:,:) :: JoutNg

END SUBROUTINE DcmpGatherCurrentNg
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE DcmpGatherCurrent1g(CoreInfo, jout1g)

USE TYPEDEF, ONLY : CoreInfo_Type

IMPLICIT NONE

TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:,:,:) :: jout1g

END SUBROUTINE DcmpGatherCurrent1g
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE DcmpGatherBndyFluxNg(RayInfo, DcmpPhiAngOutNg)

USE TYPEDEF, ONLY : RayInfo_Type

IMPLICIT NONE

TYPE (RayInfo_Type) :: RayInfo

REAL, POINTER, DIMENSION(:,:,:,:,:)  :: DcmpPhiAngOutNg

END SUBROUTINE DcmpGatherBndyFluxNg
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE DcmpGatherBndyFlux1g(RayInfo, DcmpPhiAngOut1g)

USE TYPEDEF, ONLY : RayInfo_Type

IMPLICIT NONE

TYPE (RayInfo_Type) :: RayInfo

REAL, POINTER, DIMENSION(:,:,:,:)  :: DcmpPhiAngOut1g

END SUBROUTINE DcmpGatherBndyFlux1g
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE FluxUnderRelaxation(Core, Phis1g, Phis, w, iz, ig, PE)

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
SUBROUTINE PowerUpdate_WATT(Core, Fxr, phis, power, myzb, myze, ng, lxslib, GroupInfo, PE, LBCAST_INP)
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,       Fxrinfo_type,       Cell_Type,     pin_Type, &
                         GroupInfo_Type,      XsMac_Type,         PE_TYPE
USE BenchXs,       ONLY : xskfBen
USE MacXsLib_Mod, ONLY : MacXskf
USE BasicOperation, ONLY : CP_CA, MULTI_VA
USE CNTL,             ONLY: nTracerCntl
use Material_Mod,    only: Mixture
#ifdef MPI_ENV
USE MPIComm_Mod, ONLY : BCAST
USE MPIComm_mod, ONLY : REDUCE
#endif
IMPLICIT NONE
TYPE(coreinfo_type) :: CORE
TYPE(Fxrinfo_type),POINTER :: Fxr(:, :)
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_Type) :: PE
REAL, POINTER :: phis(:, :, :)
REAL, POINTER :: Power(:, :)
REAL, POINTER :: HZ(:)
INTEGER :: myzb, myze, ng
LOGICAL :: lXsLib
LOGICAL :: LBCAST_INP

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Fxrinfo_type),POINTER :: myFxr
TYPE(XsMac_Type), SAVE :: XsMac

INTEGER :: nxy, nCoreFsr, nCoreFxr, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr
INTEGER :: ipin, icel, ifsrlocal, ifsr, ifxr, iz, itype, ig
INTEGER :: iResoGrpBeg, iResoGrpEnd, norg
INTEGER :: i, j, k, IM

REAL, POINTER :: xsmackf(:)
REAL :: HZ_LOC, vol, localpow, pwsum, F, PowerCore, PowerLevel
REAL :: Buf

END SUBROUTINE PowerUpdate_WATT
! ------------------------------------------------------------------------------------------------------------
END INTERFACE

CONTAINS
! ------------------------------------------------------------------------------------------------------------
FUNCTION setDcmpClr(lHex, l060, iClr, mocit)

IMPLICIT NONE

LOGICAL :: lHex, l060
INTEGER :: iClr, iit, mocit, setDcmpClr

INTEGER, PARAMETER :: AuxRec(2, 0:1) = [2, 1,  1, 2]
INTEGER, PARAMETER :: AuxHex(3, 0:1) = [3, 1, 2,  2, 1, 3]

iit = mod(mocit, 2)

IF (lHex) THEN
  setDcmpClr = AuxHex(iClr, iit)
  IF (l060) setDcmpClr = iClr
ELSE
  setDcmpClr = AuxRec(iClr, iit)
END IF

END FUNCTION setDcmpClr
! ------------------------------------------------------------------------------------------------------------

END MODULE MOC_MOD