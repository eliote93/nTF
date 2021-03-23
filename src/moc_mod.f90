MODULE MOC_MOD
USE PARAM
USE TYPEDEF, ONLY : TrackingDat_Type
IMPLICIT NONE
REAL, POINTER :: phis1g(:), PhiAngIn1g(:, :), phim1g(:,:)
REAL, POINTER :: phi1a1g(:, :, :), phi2a1g(:, :, :)
REAL, POINTER :: MocJout1g(:, :, :)
REAL, POINTER :: xst1g(:)                           !transport XS
REAL, POINTER :: tSrc(:)
REAL, POINTER :: tSrcAng1(:,:,:), tSrcAng2(:,:,:)
REAL, POINTER :: AxSrc1g(:), AxPxs1g(:)
REAL, POINTER :: LinSrc1g(:, :), LinPsi(:, :, :)       !Linear Source Approximation
REAL, POINTER :: srcm(:, :)                            !P1 Scattering Source term

!--- CNJ Edit : Node Majors
REAL, POINTER :: phisnm(:, :), phimnm(:, :, :), phianm(:, :, :, :), PhiAngInnm(:, :, :)
REAL, POINTER :: srcnm(:, :), srcmnm(:, :, :), SrcAngnm(:, :, :, :)
REAL, POINTER :: xstnm(:, :)
REAL, POINTER :: MocJoutnm(:, :, :, :)

!--- CNJ Edit : Domain Decomposition
REAL, POINTER :: DcmpPhiAngIn(:, :, :, :, :)
REAL, POINTER :: DcmpPhiAngOut(:, :, :, :, :)

REAL, POINTER :: wtang(:, :), wtsurf(:, :, :)
REAL, POINTER :: SrcAng(:, :, :, :)
REAL, POINTER :: SrcAng1(:, :, :), SrcAng2(:, :, :)
REAL, POINTER :: comp(:, :, :), mwt(:, :, :), mwt2(:, :, :)
REAL, POINTER :: wPhim(:) !, wPhim3D(:)
REAL, POINTER :: phia1g(:, :, :) !(ifsr, Azi, Pol) 
REAL, POINTER :: phia2g(:, :, :)
INTEGER, POINTER :: nAziRotRay(:), nAziRotRay0(:), nSegRotRay(:), OmpRayBeg(:), OmpRayEnd(:)
INTEGER, POINTER :: OmpRayBegBd(:,:), OmpRayEndBd(:,:)
INTEGER :: OmpTemp, OmpAng
INTEGER :: nOmpAng
INTEGER, POINTER :: OmpMap(:,:)
INTEGER, POINTER :: OmpRayList(:)
INTEGER, POINTER :: OmpRayNum(:,:)
!Approximation of Exponetial Function
!REAL, POINTER, SAVE :: expa(:, :), expb(:, :)
REAL, TARGET :: expa(-40000:0, 1:12), expb(-40000:0, 1:12)
REAL, TARGET :: EXPA_p(1:12, -40000:0), EXPB_p(1:12, -40000:0)
INTEGER, TARGET :: AziMap(360, 2)
INTEGER, PRIVATE :: nFsrMoc, MyZbMoc, MyZeMoc, nPinsMoc, ngMoc
INTEGER :: nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay
INTEGER :: nMaxDcmpRaySeg, nMaxDcmpCellRay, nMaxDcmpAsyRay   !--- CNJ Edit : Domain Decomposition

TYPE(TrackingDat_Type), SAVE :: TrackingDat(100)

INTERFACE

SUBROUTINE MOCSweep(Core, RayInfo, FMInfo, eigv, ng, PE, nTracerCntl, ItrCntl)
USE TYPEDEF,     ONLY : CoreInfo_Type,      RayInfo_Type,       FmInfo_Type,                                        &
                        PE_TYPE
USE CNTL,        ONLY : nTracerCntl_Type
USE itrcntl_mod, ONLY : ItrCntl_TYPE
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(PE_TYPE) :: PE
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(ItrCntl_TYPE) :: ItrCntl
REAL :: eigv
INTEGER :: ng

END SUBROUTINE

SUBROUTINE DcplMOCSweep(Core, RayInfo, FMInfo, eigv, ng, PE, GroupInfo, nTracerCntl, ItrCntl)
USE TYPEDEF,     ONLY : CoreInfo_Type,      RayInfo_Type,       FmInfo_Type,                                        &
                        PE_TYPE,            GroupInfo_Type
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

END SUBROUTINE

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
END SUBROUTINE

SUBROUTINE RayTrace(RayInfo, CoreInfo, phis, PhiAngIn, xst, src, jout, iz, ljout, FastMocLv, lAFSS)
USE PARAM
USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type
IMPLICIT NONE
TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phis(:), PhiAngIn(:, :), xst(:), src(:), jout(:, :, :)
INTEGER :: iz
LOGICAL :: ljout
INTEGER, OPTIONAL :: FastMocLv
LOGICAL, OPTIONAL :: lAFSS
END SUBROUTINE

!--- CNJ Edit : Domain Decomposition
SUBROUTINE RayTrace_Dcmp(RayInfo, CoreInfo, phisnm, phisSlope, PhiAngInnm, srcnm, srcSlope,                 &
                         xstnm, joutnm, iz, gb, ge, lJout, lLinSrcCASMO, lHybrid)
USE PARAM
USE TYPEDEF,  ONLY :  RayInfo_Type,      Coreinfo_type
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phisnm(:, :), PhiAngInnm(:, :, :)
REAL, POINTER :: srcnm(:, :), xstnm(:, :), joutnm(:, :, :, :)
REAL, POINTER :: phisSlope(:, :, :, :), srcSlope(:, :, :, :)
INTEGER :: iz, gb, ge
LOGICAL :: lJout, lLinSrcCASMO, lHybrid
END SUBROUTINE
                         
SUBROUTINE RayTrace_OMP(RayInfo, CoreInfo, phis, PhiAngIn, xst, src, jout, iz, ljout, FastMocLv, lAFSS)
USE PARAM
USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type
IMPLICIT NONE
TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phis(:), PhiAngIn(:, :), xst(:), src(:), jout(:, :, :)
INTEGER :: iz
LOGICAL :: ljout
INTEGER, OPTIONAL :: FastMocLv
LOGICAL, OPTIONAL :: lAFSS
END SUBROUTINE

SUBROUTINE RayTrace_OMP_AFSS(RayInfo, CoreInfo, phis, PhiAngIn, xst, src, jout, iz, ljout, FastMocLv, lAFSS)
USE PARAM
USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type
IMPLICIT NONE
TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phis(:), PhiAngIn(:, :), xst(:), src(:), jout(:, :, :)
INTEGER :: iz
LOGICAL :: ljout
INTEGER, OPTIONAL :: FastMocLv
LOGICAL, OPTIONAL :: lAFSS
END SUBROUTINE

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
END SUBROUTINE

!--- CNJ Edit : Node Majors, Domain Decomposition
SUBROUTINE RayTraceNM_OMP(RayInfo, CoreInfo, phisnm, PhiAngInnm, xstnm, srcnm, joutnm,                              &
                          iz, gb, ge, ljout, lDomainDcmp, FastMocLv)
USE PARAM
USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type
IMPLICIT NONE
TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phisnm(:, :), PhiAngInnm(:, :, :), xstnm(:, :), srcnm(:, :), joutnm(:, :, :, :)
INTEGER :: iz, gb, ge
LOGICAL :: ljout, lDomainDcmp
INTEGER, OPTIONAL :: FastMocLv
END SUBROUTINE

SUBROUTINE RayTraceDcmp_NM(RayInfo, CoreInfo, phisnm, PhiAngInnm, xstnm, srcnm,                                     &
                           joutnm, iz, iasy, gb, ge, ljout)
USE PARAM
USE TYPEDEF, ONLY :   RayInfo_Type,      Coreinfo_type
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phisnm(:, :), PhiAngInnm(:, :, :), xstnm(:, :), srcnm(:, :), joutnm(:, :, :, :)
INTEGER :: iz, iasy, gb, ge
LOGICAL :: ljout
END SUBROUTINE
                           
SUBROUTINE TrackRotRayNM_Dcmp(RayInfo, CoreInfo, TrackingDat, phis, jout, ljout, DcmpAsyRay, iz, gb, ge)
USE PARAM
USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type, TrackingDat_Type, DcmpAsyRayInfo_Type
TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
TYPE(TrackingDat_Type) :: TrackingDat
TYPE(DcmpAsyRayInfo_Type) :: DcmpAsyRay
REAL, POINTER :: phis(:, :), jout(:, :, :, :)
LOGICAL :: ljout
INTEGER :: iz, gb, ge
END SUBROUTINE

SUBROUTINE RayTraceP1(RayInfo, CoreInfo, phis, phim, PhiAngIn, xst, src, srcm, jout, iz, ljout, ScatOd, FastMocLv, lAFSS)
USE PARAM
USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type
IMPLICIT NONE
!Input Arguments
TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phis(:), PhiAngIn(:, :), xst(:), src(:), jout(:, :, :), srcm(:, :)
REAL, POINTER :: phim(:, :) ! BYS Correct
INTEGER :: iz
LOGICAL :: ljout
INTEGER :: ScatOd
INTEGER, OPTIONAL :: FastMocLv
LOGICAL, OPTIONAL :: lAFSS
END SUBROUTINE

SUBROUTINE RayTraceP1_OMP(RayInfo, CoreInfo, phis, phim, PhiAngIn, xst, src, srcm, jout, iz, ljout, ScatOd, FastMocLv, lAFSS)
USE PARAM
USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type
IMPLICIT NONE
!Input Arguments
TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phis(:), PhiAngIn(:, :), xst(:), src(:), jout(:, :, :), srcm(:, :)
REAL, POINTER :: phim(:, :) ! BYS Correct
INTEGER :: iz
LOGICAL :: ljout
INTEGER :: ScatOd
INTEGER, OPTIONAL :: FastMocLv
LOGICAL, OPTIONAL :: lAFSS
END SUBROUTINE

SUBROUTINE RayTraceP1_OMP_AFSS(RayInfo, CoreInfo, phis, phim, PhiAngIn, xst, src, srcm, jout, iz, ljout, ScatOd, FastMocLv, lAFSS)
USE PARAM
USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type
IMPLICIT NONE
!Input Arguments
TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phis(:), PhiAngIn(:, :), xst(:), src(:), jout(:, :, :)
REAL, POINTER :: phim(:, :), srcm(:, :) ! BYS Correct >> pointer
INTEGER :: iz
LOGICAL :: ljout
INTEGER :: ScatOd
INTEGER, OPTIONAL :: FastMocLv
LOGICAL, OPTIONAL :: lAFSS
END SUBROUTINE

!--- CNJ Edit : Angular Multigrid Ray Tracing
SUBROUTINE RayTraceP1_Multigrid(RayInfo, CoreInfo, phis, phim, PhiAngIn, xst, src, srcm, jout, iz, ScatOd, lJout)
USE PARAM
USE TYPEDEF,        ONLY : RayInfo_Type,    CoreInfo_type
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phis(:), phim(:, :), PhiAngIn(:, :), xst(:)
REAL, POINTER :: src(:), srcm(:, :), jout(:, :, :)
INTEGER :: iz
LOGICAL :: ljout
INTEGER :: ScatOd
END SUBROUTINE

!--- CNJ Edit : Node Majors, Angular Flux Saving
SUBROUTINE RayTraceDcmp_Pn(RayInfo, CoreInfo, phisnm, phimnm, PhiAngInnm, xstnm, srcnm, srcmnm,                     &
                           joutnm, iz, iAsy, gb, ge, ScatOd, lJout)
USE PARAM
USE TYPEDEF, ONLY :   RayInfo_Type,      Coreinfo_type
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phisnm(:, :), phimnm(:, :, :), PhiAngInnm(:, :, :), xstnm(:, :)
REAL, POINTER :: srcnm(:, :), srcmnm(:, :, :), joutnm(:, :, :, :)
INTEGER :: iz, iAsy, gb, ge
LOGICAL :: ljout
INTEGER :: ScatOd
END SUBROUTINE

SUBROUTINE TrackRotRayPn_Dcmp(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, Jout, ljout, iz, gb, ge)
USE PARAM
USE TYPEDEF,  ONLY :  RayInfo_Type,         Coreinfo_type,       TrackingDat_Type,    DcmpAsyRayInfo_Type
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
TYPE(TrackingDat_Type) :: TrackingDat
TYPE(DcmpAsyRayInfo_Type) :: DcmpAsyRay
REAL, POINTER :: Jout(:, :, :, :)
LOGICAL :: lJout
INTEGER :: iz, gb, ge, ScatOd
END SUBROUTINE

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
END SUBROUTINE

!--- CNJ Edit : CASMO Linear Source, Domain Decomposition
SUBROUTINE RayTraceLS_CASMO(RayInfo, CoreInfo, phisnm, phisSlope, PhiAngInnm, srcnm, srcSlope,                      &
                            xstnm, joutnm, iz, gb, ge, ljout, lDomainDcmp)
USE PARAM
USE TYPEDEF,    ONLY :  RayInfo_Type,       Coreinfo_type
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phisnm(:, :), phisSlope(:, :, :, :), PhiAngInnm(:, :, :)
REAL, POINTER :: srcnm(:, :), srcSlope(:, :, :, :), xstnm(:, :), joutnm(:, :, :, :)
INTEGER :: iz, gb, ge
LOGICAL :: ljout, lDomainDcmp
END SUBROUTINE

SUBROUTINE RayTraceDcmp_LSCASMO(RayInfo, CoreInfo, phisnm, phisSlope, PhiAngInnm, srcnm, srcSlope,                  &
                                xstnm, joutnm, iz, iasy, gb, ge, ljout)
USE PARAM
USE TYPEDEF,    ONLY :  RayInfo_Type,       Coreinfo_type
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phisnm(:, :), PhiAngInnm(:, :, :)
REAL, POINTER :: srcnm(:, :), xstnm(:, :), joutnm(:, :, :, :)
REAL, POINTER :: phisSlope(:, :, :, :), srcSlope(:, :, :, :)
INTEGER :: iz, iasy, gb, ge
LOGICAL :: ljout
END SUBROUTINE
                                
SUBROUTINE TrackRotRayLSDcmp_CASMO(RayInfo, CoreInfo, TrackingDat, phisC, phimx, phimy, jout, DcmpAsyRay, ljout, iz, gb, ge)
USE PARAM
USE TYPEDEF,    ONLY :  RayInfo_Type,       Coreinfo_type,      DcmpAsyRayInfo_Type,   TrackingDat_Type
IMPLICIT NONE

TYPE(RayInfo_Type), INTENT(INOUT) :: RayInfo
TYPE(CoreInfo_Type), INTENT(INOUT) :: CoreInfo
TYPE(TrackingDat_Type), INTENT(INOUT) :: TrackingDat
TYPE(DcmpAsyRayInfo_Type) :: DcmpAsyRay
REAL, POINTER :: phisC(:, :), phimx(:, :, :), phimy(:, :, :), jout(:, :, :, :)
LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: iz, gb, ge
END SUBROUTINE

SUBROUTINE SetRtMacXs(core, Fxr, xstr, iz, ig, ng, lxslib, lTrCorrection, lRST, lssph, lssphreg, PE)
USE PARAM
USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, PE_TYPE
USE BenchXs, ONLY : GetXstrBen
IMPLICIT NONE
TYPE(coreinfo_type) :: Core
TYPE(Fxrinfo_type) :: Fxr(:)
TYPE(PE_type) :: PE
REAL, POINTER :: xstr(:)

INTEGER :: myzb, myze, ig, iz, ng
logical :: lxslib, lTrCorrection, lRST, lssph, lssphreg
END SUBROUTINE

SUBROUTINE SetRtMacXs_Cusping(core, FmInfo, Fxr, xstr, phis, iz, ig, ng, lxslib, lTrCorrection, lRST, lssph, lssphreg, PE)
USE PARAM
USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, PE_TYPE, FmInfo_Type
USE BenchXs, ONLY : GetXstrBen
IMPLICIT NONE
TYPE(coreinfo_type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(Fxrinfo_type) :: Fxr(:)
TYPE(PE_type) :: PE
REAL, POINTER :: xstr(:)
REAL, POINTER :: phis(:, :, :)

INTEGER :: myzb, myze, ig, iz, ng
logical :: lxslib, lTrCorrection, lRST, lssph, lssphreg
END SUBROUTINE

!--- CNJ Edit : Node Majors    
SUBROUTINE SetRtMacXsNM(core, Fxr, xstnm, iz, ng, lxslib, lTrCorrection, lRST, lssph, lssphreg, PE)
USE PARAM
USE TYPEDEF,  ONLY : coreinfo_type,    Fxrinfo_type,    PE_TYPE
IMPLICIT NONE
TYPE(coreinfo_type) :: Core
TYPE(Fxrinfo_type) :: Fxr(:)
TYPE(PE_TYPE) :: PE
REAL, POINTER :: xstnm(:, :)
INTEGER :: iz, ng
logical :: lxslib, lTrCorrection, lRST, lssph, lssphreg
END SUBROUTINE

SUBROUTINE SetRtMacXsNM_Cusping(core, FmInfo, Fxr, xstnm, iz, ng, lxslib, lTrCorrection, lRST, lssph, lssphreg, PE)
USE PARAM
USE TYPEDEF,  ONLY : coreinfo_type,    Fxrinfo_type,    PE_TYPE,   FmInfo_Type
IMPLICIT NONE
TYPE(coreinfo_type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(Fxrinfo_type) :: Fxr(:)
TYPE(PE_TYPE) :: PE
REAL, POINTER :: xstnm(:, :)
INTEGER :: iz, ng
logical :: lxslib, lTrCorrection, lRST, lssph, lssphreg
END SUBROUTINE

SUBROUTINE SetRtSrc(Core, Fxr, src, phis, psi, axsrc, xstr1g, eigv, iz, ig, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix, PE)
USE PARAM
USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, GroupInfo_Type, PE_TYPE
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:)
TYPE(GroupInfo_Type) :: GroupInfo
REAL, POINTER :: src(:), phis(:, :, :), psi(:, :), AxSrc(:), xstr1g(:)
REAL :: eigv
INTEGER :: myzb, myze, ig, ng, iz, ifsr, ifxr
LOGICAL :: lxslib, lscat1, l3dim
LOGICAL :: lNegFix
TYPE(PE_TYPE) :: PE

END SUBROUTINE

SUBROUTINE SetRtSrc_Cusping(Core, FmInfo, Fxr, src, phis, psi, axsrc, xstr1g, eigv, iz, ig, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix, PE)
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

END SUBROUTINE

!--- CNJ Edit : Node Majors 
SUBROUTINE SetRtSrcNM(Core, Fxr, srcnm, phisnm, psi, AxSrc, xstnm, eigv, iz,                                        &
                      gb, ge, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix, PE, Offset)
USE PARAM
USE TYPEDEF,        ONLY : coreinfo_type,   Fxrinfo_type,   GroupInfo_Type,  PE_TYPE
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:)
TYPE(GroupInfo_Type):: GroupInfo
TYPE(PE_TYPE) :: PE
REAL, POINTER :: srcnm(:, :), phisnm(:, :), psi(:, :), AxSrc(:, :, :), xstnm(:, :)
REAL :: eigv
INTEGER :: gb, ge, ng, iz, ifsr, ifxr
LOGICAL :: lxslib, lscat1, l3dim
LOGICAL :: lNegFix
INTEGER, OPTIONAL :: Offset
END SUBROUTINE

SUBROUTINE SetRtSrcNM_Cusping(Core, FmInfo, Fxr, srcnm, phisnm, psi, AxSrc, xstnm, eigv, iz,                                        &
                      gb, ge, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix, PE, Offset)
USE PARAM
USE TYPEDEF,        ONLY : coreinfo_type,   Fxrinfo_type,   GroupInfo_Type,  PE_TYPE, FmInfo_Type
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
END SUBROUTINE

SUBROUTINE SetRtLinSrc(Core, Fxr, RayInfo, src,LinSrc, LinPsi, LinSrcSlope, xstr1g,                                 &
                    eigv, iz, ig, ng, GroupInfo, l3dim, lxslib, lscat1)
USE PARAM
USE TYPEDEF, ONLY : coreinfo_type,      Fxrinfo_type,      RayInfo_Type,   &
                    GroupInfo_Type,     XsMac_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:)
TYPE(RayInfo_Type) :: RayInfo
TYPE(GroupInfo_Type):: GroupInfo
REAL, POINTER :: Src(:), Linsrc(:, :), LinPsi(:, :, :), LinSrcSlope(:, :, :, :), xstr1g(:)
REAL :: eigv
INTEGER :: myzb, myze, ig, ng, iz
LOGICAL :: lxslib, lscat1, l3dim

END SUBROUTINE

!--- CNJ Edit : CASMO Linear Source
SUBROUTINE SetRtLinSrc_CASMO(Core, Fxr, RayInfo, phisnm, phisSlope, srcnm, srcSlope, psi, psiSlope, axsrc, xstnm,   &
                             eigv, iz, gb, ge, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix)
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,   Fxrinfo_type,   RayInfo_Type,   GroupInfo_Type
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

END SUBROUTINE
                             
SUBROUTINE SetRtP1Src(Core, Fxr, srcm, phim, xstr1g, iz, ig, ng, GroupInfo, l3dim, lxslib, lscat1, ScatOd, PE)
USE PARAM
USE TYPEDEF, ONLY : coreinfo_type,      Fxrinfo_type,      RayInfo_Type,   &
                    GroupInfo_Type,     XsMac_Type,        PE_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:)
TYPE(GroupInfo_Type):: GroupInfo
REAL, POINTER :: srcm(:, :), phim(:, :, :, :), xstr1g(:)
INTEGER :: myzb, myze, ig, ng, iz
LOGICAL :: lxslib, lscat1, l3dim
INTEGER :: ScatOd
TYPE(PE_Type) :: PE
END SUBROUTINE

!--- CNJ Edit : Node Majors
SUBROUTINE SetRtP1SrcNM(Core, Fxr, srcmnm, phimnm, xstnm, iz, gb, ge, ng, GroupInfo, lxslib, ScatOd, PE, Offset)
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,     Fxrinfo_type,   GroupInfo_Type,     PE_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:)
TYPE(GroupInfo_Type):: GroupInfo
REAL, POINTER :: srcmnm(:, :, :)
REAL, POINTER :: phimnm(:, :, :)
REAL, POINTER :: xstnm(:, :)
REAL :: eigv
INTEGER :: gb, ge, ng, iz, ifsr, ifxr
LOGICAL :: lxslib, lscat1, l3dim
INTEGER :: ScatOd
TYPE(PE_Type) :: PE
INTEGER, OPTIONAL :: Offset
END SUBROUTINE

SUBROUTINE PseudoAbsorption(Core, Fxr, src, phis1g, AxPXS, xstr1g, iz, ig, ng, GroupInfo, l3dim)
USE PARAM
USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, GroupInfo_Type

IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:)
TYPE(GroupInfo_Type) :: GroupInfo
REAL :: phis1g(:), AxPXS(:)
REAL, POINTER :: src(:), xstr1g(:)
INTEGER :: ig, ng, iz
LOGICAL :: l3dim

END SUBROUTINE

!--- CNJ Edit : Node Majors
SUBROUTINE PseudoAbsorptionNM(Core, Fxr, AxPXS, xstnm, iz, ng, GroupInfo, l3dim)
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,          Fxrinfo_type,          Cell_Type,      & 
                         pin_Type,               GroupInfo_Type,        XsMac_Type
USE BenchXs,      ONLY : GetChiBen,              xssben
USE MacXsLib_mod, ONLY : MacXsScatMatrix
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:)
TYPE(GroupInfo_Type):: GroupInfo
REAL, POINTER :: AxPXS(:, :, :), xstnm(:, :)
REAL :: eigv
INTEGER :: iz, ng
LOGICAL :: l3dim
END SUBROUTINE

SUBROUTINE PsiUpdate(Core, Fxr, phis, psi, myzb, myze, ng, lxslib, GroupInfo)
USE PARAM
USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, GroupInfo_Type
IMPLICIT NONE
TYPE(coreinfo_type) :: CORE
TYPE(Fxrinfo_type), POINTER :: Fxr(:, :)
TYPE(GroupInfo_Type) :: GroupInfo
REAL, POINTER :: phis(:, :, :), Psi(:, :)
INTEGER :: myzb, myze, ng
LOGICAL :: lXsLib

END SUBROUTINE

SUBROUTINE PowerUpdate(Core, Fxr, phis, Power, myzb, myze, ng, lxslib, GroupInfo, PE)
USE PARAM
USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, GroupInfo_Type, PE_TYPE
IMPLICIT NONE
TYPE(coreinfo_type) :: CORE
TYPE(Fxrinfo_type), POINTER :: Fxr(:, :)
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_TYPE) :: PE
REAL, POINTER :: phis(:, :, :), Power(:, :)
INTEGER :: myzb, myze, ng
LOGICAL :: lXsLib

END SUBROUTINE

SUBROUTINE LinPsiUpdate(Core, Fxr, LinPsi, LinSrcSlope, myzb, myze, ng, lxslib, GroupInfo)
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,       Fxrinfo_type,    GroupInfo_Type
IMPLICIT NONE
TYPE(coreinfo_type) :: CORE
TYPE(Fxrinfo_type), POINTER :: Fxr(:, :)
TYPE(GroupInfo_Type) :: GroupInfo
REAL, POINTER :: LinPsi(:, :, :)
REAL, POINTER :: LinSrcSlope(:, :, :, :)
INTEGER :: myzb, myze, ng
LOGICAL :: lXsLib

END SUBROUTINE 

!--- CNJ Edit : CASMO Linear Source
SUBROUTINE LinPsiUpdate_CASMO(Core, Fxr, phisSlope, psiSlope, myzb, myze, ng, lxslib, GroupInfo)
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,       Fxrinfo_type,       GroupInfo_Type
IMPLICIT NONE
TYPE(coreinfo_type) :: CORE
TYPE(Fxrinfo_type),POINTER :: Fxr(:, :)
TYPE(GroupInfo_Type) :: GroupInfo
REAL, POINTER :: phisSlope(:, :, :, :)
REAL, POINTER :: psiSlope(:, :, :)
INTEGER :: myzb, myze, ng
LOGICAL :: lXsLib

END SUBROUTINE

SUBROUTINE CellPsiUpdate(CORE, psi, psic, myzb, myze)
USE PARAM
USE TYPEDEF, ONLY : coreinfo_type
IMPLICIT NONE
TYPE(coreinfo_type) :: CORE
REAL, POINTER :: Psi(:, :), psiC(:, :)
INTEGER :: myzb, myze

END SUBROUTINE


SUBROUTINE UpdateEigv(Core, psi, psid, eigv, peigv, myzb, myze, PE)
USE PARAM
USE TYPEDEF, ONLY : coreinfo_type, PE_TYPE
USE BenchXs, ONLY : xsnfBen
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
TYPE(coreinfo_type) :: CORE
REAL, POINTER :: Psi(:, :), psid(:, :)
INTEGER :: myzb, myze, ng
REAL :: eigv, peigv
TYPE(PE_TYPE), OPTIONAL :: PE

END SUBROUTINE

!Approximation of Exponential Function
Subroutine ApproxExp(PolarAng, npr)
!Subroutine ApproxExp(PolarAng, npr, Expa, Expb)
USE PARAM
USE TYPEDEF, ONLY : PolarAngle_Type
IMPLICIT NONE
TYPE(PolarAngle_Type) :: PolarAng(npr)
!REAL :: Expa(:, :), Expb(:, :)
INTEGER :: npr

END SUBROUTINE

Function PsiErr(Core, psi, psid, myzb, myze, PE)
USE PARAM
USE TYPEDEF, ONLY : coreinfo_type, PE_TYPE
IMPLICIT NONE
TYPE(coreinfo_type) :: CORE
REAL :: PsiErr
REAL, POINTER :: Psi(:, :), psid(:, :)
INTEGER :: myzb, myze
TYPE(PE_TYPE) :: PE

END FUNCTION

FUNCTION MocResidual(Core, FmInfo,  eigv, GroupInfo, ng, PE, nTracerCntl)
USE TYPEDEF,  ONLY : CoreInfo_Type,     FmInfo_Type,     GroupInfo_Type
USE cntl,     ONLY : nTracerCntl_Type
USE PE_MOD,   ONLY : PE_TYPE
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
REAL :: eigv
TYPE(GroupInfo_Type) :: GroupInfo
INTEGER :: ng
TYPE(PE_TYPE) :: PE
TYPE(nTracerCntl_Type) :: nTracerCntl
REAL :: MocResidual

END FUNCTION

SUBROUTINE AddBuckling(Core, Fxr, xstr1g, bsq, iz, ig, ng, lxslib, lRST)
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,          Fxrinfo_type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(Fxrinfo_type), POINTER :: Fxr(:,:)
REAL :: xstr1g(:)
INTEGER :: iz, ig, ng
REAL :: Bsq
LOGICAL :: lxslib, lRST

END SUBROUTINE

!--- CNJ Edit : Node Majors    
SUBROUTINE AddBucklingNM(Core, Fxr, xstnm, bsq, iz, ng, lxslib, lRST)
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,     Fxrinfo_type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(Fxrinfo_type), POINTER :: Fxr(:,:)
REAL, POINTER :: xstnm(:, :)
INTEGER :: iz, ng
REAL :: Bsq
LOGICAL :: lxslib, lRST
END SUBROUTINE

SUBROUTINE AddConstSrc(Core, Fxr, Src, xstr1g, ConstSrc, iz, ig, ng)
USE PARAM
USE TYPEDEF,       ONLY : CoreInfo_Type,        Fxrinfo_type,                          &
                          Cell_Type,            pin_Type
TYPE(CoreInfo_Type) :: COre
TYPE(FxrInfo_Type) :: Fxr(:)
REAL, POINTER :: Src(:), xstr1g(:)
REAL :: ConstSrc
INTEGER :: iz, ig, ng

END SUBROUTINE

FUNCTION FxrAvgPhi(Core, Fxr, Phis, ipin, iLocalfxr, iz, ng, PE)
!USE PARAM
USE TYPEDEF,     ONLY : CoreInfo_Type,    PE_Type,     FxrInfo_Type,               &
                        Cell_Type,    Pin_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(PE_Type) :: PE
REAL, POINTER :: phis(:, :, :)
INTEGER :: iLocalfxr, ipin, iz, ng
REAL :: FxrAvgPhi(ng)

END FUNCTION

!--- CNJ Edit : Domain Decomposition
SUBROUTINE DcmpLinkBoundaryFlux(CoreInfo, RayInfo, PhiAngIn, DcmpPhiAngIn, DcmpPhiAngOut, gb, ge, color)
USE PARAM
USE TYPEDEF,    ONLY : CoreInfo_Type,   RayInfo_Type,    DcmpAsyRayInfo_Type
USE PE_Mod,     ONLY : PE
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(RayInfo_Type) :: RayInfo
REAL, POINTER :: PhiAngIn(:, :, :)
REAL, POINTER :: DcmpPhiAngIn(:, :, :, :, :), DcmpPhiAngOut(:, :, :, :, :)
INTEGER :: gb, ge
INTEGER, OPTIONAL :: color
END SUBROUTINE

!--- CNJ Edit : Domain Decomposition + MPI
SUBROUTINE DcmpScatterXS(CoreInfo, xst)
USE PARAM
USE TYPEDEF,    ONLY : CoreInfo_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: xst(:, :)
END SUBROUTINE

SUBROUTINE DcmpScatterSource(CoreInfo, src, srcm)
USE PARAM
USE TYPEDEF,    ONLY : CoreInfo_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: src(:, :)
REAL, POINTER, OPTIONAL :: srcm(:, :, :)
END SUBROUTINE

SUBROUTINE DcmpScatterBoundaryFlux(RayInfo, PhiAngIn, DcmpPhiAngIn)
USE PARAM
USE TYPEDEF,    ONLY : RayInfo_Type
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
REAL, POINTER :: PhiAngIn(:, :, :), DcmpPhiAngIn(:, :, :, :, :)
END SUBROUTINE

SUBROUTINE DcmpGatherFlux(CoreInfo, phis, phim)
USE PARAM
USE TYPEDEF,    ONLY : CoreInfo_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phis(:, :)
REAL, POINTER, OPTIONAL :: phim(:, :, :)
END SUBROUTINE

SUBROUTINE DcmpGatherCurrent(CoreInfo, jout)
USE PARAM
USE TYPEDEF,    ONLY : CoreInfo_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: jout(:, :, :, :)
END SUBROUTINE

SUBROUTINE DcmpGatherBoundaryFlux(RayInfo, DcmpPhiAngOut)
USE PARAM
USE TYPEDEF,    ONLY : RayInfo_Type
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
REAL, POINTER :: DcmpPhiAngOut(:, :, :, :, :)
END SUBROUTINE    

SUBROUTINE FluxUnderRelaxation(Core, Phis1g, Phis, w, iz, ig, PE)
USE PARAM
USE TYPEDEF,  ONLY : CoreInfo_Type,  PE_TYPE,  Pin_Type, Cell_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
REAL :: Phis1g(:)
REAL :: Phis(:, :, :)
REAL :: w
INTEGER :: iz, ig
TYPE(PE_TYPE) :: PE

END SUBROUTINE

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

END SUBROUTINE

SUBROUTINE CurrentUnderRelaxation(Core, Jout1g, Jout, w, iz, ig, PE)
USE PARAM
USE TYPEDEF,  ONLY : CoreInfo_Type,  PE_TYPE
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
REAL :: Jout1g(:, :, :)
REAL :: Jout(:, :, :, :, :)
REAL :: w
INTEGER :: iz, ig
TYPE(PE_TYPE) :: PE

END SUBROUTINE

SUBROUTINE GetNeighborMocFlux(phis, neighphis, nFsr, myzb, myze, gb, ge, nz, AxBC)
IMPLICIT NONE
REAL, POINTER :: phis(:, :, :), neighphis(:, :, :)
INTEGER :: nfsr, myzb, myze, gb, ge, nz
INTEGER :: AxBC(2)

END SUBROUTINE

END INTERFACE

CONTAINS 

SUBROUTINE SetMocEnvironmentVariables(nFsr, nPins, myzb, myze, ng )
IMPLICIT NONE
INTEGER :: nFsr, nPins, myzb, myze, ng
nFsrMoc = nFsr
MyZbMoc = myzb
MyZeMoc = myze
nPinsMoc = npins 
ngMoc = ng
END SUBROUTINE

END MODULE
