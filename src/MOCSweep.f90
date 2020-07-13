#include <defines.h>
SUBROUTINE MOCSweep(Core, RayInfo, FMInfo, eigv, ng, PE, nTracerCntl, ItrCntl)
USE PARAM
USE TYPEDEF,     ONLY : CoreInfo_Type,      RayInfo_Type,       FmInfo_Type,        &
                        PE_TYPE,            FxrInfo_Type
USE CNTL,        ONLY : nTracerCntl_Type
USE itrcntl_mod, ONLY : ItrCntl_TYPE
USE CORE_MOD,    ONLY : GroupInfo,                                                  &
                        !--- CNJ Edit : CASMO Linear Source
                        srcSlope,           phisSlope,         psiSlope
USE MOC_MOD,     ONLY : RayTrace,           SetRtMacXs,        SetRtSrc,            &
                        RayTrace_OMP,       RayTraceP1,                             &
                        !--- CNJ Edit : Angular Multigrid Ray Tracing
                        RayTraceP1_Multigrid,                                       &
                        RayTraceLS,         LinPsiUpdate,      SetRtLinSrc,         &
                        SetRtP1Src,         AddBuckling,                            &
                        PsiUpdate,          CellPsiUpdate,     UpdateEigv,          &
                        MocResidual,        PsiErr,            PseudoAbsorption,    &
                        PowerUpdate,                                                &
                        FluxUnderRelaxation,                                        &
                        FluxInUnderRelaxation,                                      &
                        CurrentUnderRelaxation,                                     &
                        phis1g,             phim1g,            MocJout1g,           &
                        xst1g,              tSrc,              AxSrc1g,             &
                        LinSrc1g,           LinPsi,            PhiAngin1g,          &
                        srcm,                                                       &
                        !--- CNJ Edit : CASMO Linear Source
                        RayTraceLS_CASMO,   SetRTLinSrc_CASMO, LinPsiUpdate_CASMO,  &
                        !--- CNJ Edit : Node Majors
                        SetRtMacXsNM,       SetRtSrcNM,        AddBucklingNM,       &
                        SetRtP1SrcNM,       PseudoAbsorptionNM,                     &
                        phisnm,             PhiAngInnm,        MocJoutnm,           &
                        xstnm,              srcnm,             phimnm,              &
                        srcmnm,             RayTraceNM_OMP,                         &
                        !--- CNJ Edit : Domain Decomposition
                        DcmpPhiAngIn,               DcmpPhiAngOut,                  &
                        DcmpScatterXS,              DcmpScatterBoundaryFlux,        &
                        DcmpScatterSource,          DcmpGatherBoundaryFlux,         &
                        DcmpGatherFlux,             DcmpGatherCurrent,              &
                        DcmpLinkBoundaryFlux
USE SUbGrp_Mod,  ONLY : FxrChiGen
USE IOUTIL,      ONLY : message
USE Timer,       ONLY : nTracer_dclock, TimeChk
USE FILES,       ONLY : io8
#ifdef __INTEL_COMPILER
USE IFPORT,      ONLY : HOSTNM
#endif
#ifdef MPI_ENV
USE MPICOMM_MOD, ONLY : BCAST,              MPI_SYNC,           MPI_MAX_REAL
#endif
use SubChCoupling_mod,    only: last_TH
USE XSLIB_MOD,   ONLY : igresb,igrese
USE SPH_mod,     ONLY : ssphf,ssphfnm,calcPinSSPH
USE HexData,     ONLY : nInnMOCItr
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(PE_TYPE) :: PE
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(ItrCntl_TYPE) :: ItrCntl
REAL :: eigv
INTEGER :: ng


TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
REAL, POINTER :: PHIS(:,:,:)                          !Scalar Flux Vectors
REAL, POINTER :: PSI(:,:), PSID(:,:)                  !Fsr FIssion Scource
REAL, POINTER :: PSIC(:,:), PSICD(:,:)                !Cell Fission Source Term
REAL, POINTER :: AxSrc(:, :, :), AxPXS(:, :, :)
REAL, POINTER :: RadJout(:, :, :, :, :)

REAL, POINTER :: LinSrcSlope(:, :, :, :)
REAL, POINTER :: phim(:, :, :, :)

!Iteration Control
INTEGER :: nitermax, nitermin
REAL :: eigconv, psiconv, resconv

!Local Variable
INTEGER :: ig, iz, ib, ie, i, l, ITER, jsweep, InIter, nInIter
INTEGER :: myzb, myze
INTEGER :: nFsr, nxy, nxya, nbd, nModRay, nPhiAngSv, nPolarAngle, nGroupInfo
INTEGER :: GrpBeg, GrpEnd

REAL :: psipsi, psipsid, eigerr,fiserr, peigv, reserr, errdat(3)
REAL :: TimeRt1gBeg, TimeRt1gEnd, TimeRtElapsed, TimeNodeElapsed, MocTimeBeg, MocTimeEnd
REAL :: TimeRtngBeg, TimeRtngEnd, TimeRtngElapsed(2)
LOGICAL :: lJout, lxslib, l3dim, lscat1, lTrCorrection, lRST, lAFSS, lssph, lssphreg
LOGICAL :: MASTER, SLAVE, RTMaster, RTSlave, CMFDMaster, CMFDSlave
LOGICAL :: lDmesg                  !Detailed Message

#ifdef __PGI
#ifdef __linux
INTEGER :: hostnm
#endif
#endif
CHARACTER(80) :: hostname

!CMFD Time Check
MOCTimeBeg = nTracer_dclock(FALSE, FALSE)

!Pointing
FXR => FmInfo%Fxr
PHIS => FmInfo%PHIS; PSI => FmInfo%PSI
PSID => FmInfo%PSID; PSIC => FmInfo%PSIC
PSICD => FmInfo%PSICD; !PhiAngin => FmInfo%PhiAngin
RadJout => FmInfo%RadJout !ninout/ 1:in 2:out 3:surfphi 16/02/11 BYS edit;
AxSrc => FmInfo%AxSrc; AxPXS => FmInfo%AxPXS

!Parallel Environment In
Master = PE%Master; Slave = PE%Slave
RTMaster = PE%RTMaster; RTSlave = PE%RTSlave


IF(nTracerCntl%lscat1) phim => FmInfo%phim
IF(nTracerCntl%lLinSrc) LinSrcSlope => FmInfo%LinSrcSlope

IF(nTracerCntl%lsSPHreg) call calcPinSSPH(Core, Fxr, PE)

CALL omp_set_num_threads(PE%nThread)

!Local Variable
nbd = 4
myzb = PE%myzb; myze = PE%myze
nFsr = Core%nCoreFsr; nxy = Core%nxy; nxya = Core%nxya
nPolarAngle = RayInfo%nPolarAngle
nModRay = RayInfo%nModRay
nPhiAngSv = RayInfo%nPhiAngSv
iter = 0; nGroupInfo = 2
nInIter = 2

IF (nTracerCntl%lHex) nInIter = nInnMOCItr

!Control Variable
lJout = TRUE
lxslib = nTracerCntl%lXslib; l3dim = nTracerCntl%l3dim
lscat1 = nTracerCntl%lscat1
lTrCorrection = nTracerCntl%lTrCorrection
lRST = nTracerCntl%lRST
lssph = nTracerCntl%lsSPH
lssphreg = nTracerCntl%lsSPHreg

!Iteration Control Variable
nitermax = itrcntl%MOCItrCntl%nitermax
nitermin = itrcntl%MOCItrCntl%nitermin
psiconv = itrcntl%psiconv; eigconv = itrcntl%eigconv; resconv = itrcntl%resconv
IF(RTMASTER) THEN
  psid = psi; psicd = psic
  CALL PsiUpdate(Core, Fxr, phis, psi, myzb, myze, ng, lxslib, GroupInfo)
  CALL CellPsiUpdate(Core, Psi, psic, myzb, myze)
ENDIF
!Output Message Control
lDmesg = TRUE
IF(ng .gt. 10) lDmesg = FALSE
IF(.NOT. GroupInfo%lUpScat) nGroupInfo = 1

DO iter = 1, ItrCntl%MocItrCntl%nitermax

IF (.NOT. nTracerCntl%lNodeMajor) THEN

  itrcntl%mocit = itrcntl%mocit + 1; TimeRtElapsed = 0.; TimeNodeElapsed = 0.
  WRITE(mesg, '(a22,I5,a3)') 'Performing Ray Tracing', itrcntl%mocit, '...'
  IF(Master) CALL message(io8, TRUE, TRUE, mesg)
  IF(RTMaster) CALL FxrChiGen(Core, Fxr, FmInfo, GroupInfo, PE, nTracerCntl)
  IF(nTracerCntl%lLinSrc .AND. RTMaster) THEN
    CALL LinPsiUpdate(Core, Fxr, LinPsi, LinSrcSlope, myzb, myze, ng, lxslib, GroupInfo)
  ENDIF
  DO jsweep =1, nGroupInfo
    GrpBeg = 1; GrpEnd = ng
    IF(jsweep .GT. 1) THEN
      GrpBeg = GroupInfo%UpScatRange(1); GrpEnd = GroupInfo%UpScatRange(2)
    ENDIF
    DO ig = GrpBeg, GrpEnd
      TimeRt1gBeg = nTracer_dclock(FALSE, FALSE)
      !Set Transport XS
      DO iz = myzb, myze
        ljout = FALSE
        IF (.NOT. Core%lFuelPlane(iz) .AND. nTracerCntl%lAxRefFDM) CYCLE
        IF (nTracerCntl%lMOCUR) ljout = TRUE
        IF (RTMASTER) AxSrc1g = AxSrc(:, iz, ig)
        DO InIter = 1, nInIter
          IF(InIter .EQ. nInIter) ljout = TRUE
          IF(RTMASTER) THEN
            CALL SetRtMacXs(Core, Fxr(:, iz), xst1g, iz, ig, ng, lxslib, lTrCorrection, lRST, lssph, lssphreg, PE)
#ifdef LkgSplit
            CALL PseudoAbsorption(Core, Fxr(:, iz), tsrc, phis(:, iz, ig), AxPXS(:, iz, ig), xst1g,                 &
                                  iz, ig, ng, GroupInfo, l3dim)
#endif
#ifdef Buckling
            IF(nTracerCntl%lBsq) CALL AddBuckling(Core, Fxr, xst1g, nTracerCntl%Bsq, iz, ig, ng, lxslib, lRST)
#endif
          !Set Source
            CALL SetRtSrc(Core, Fxr(:, iz), tsrc, phis, psi, AxSrc1g, xst1g,                                        &
                          eigv, iz, ig, ng, GroupInfo, l3dim, lXslib, lscat1, FALSE, PE)
            PhiAngin1g = FmInfo%PhiAngin(:, : ,iz, ig)
            IF (lScat1) phim1g = phim(:, : ,iz, ig)
          ENDIF
          IF(.NOT. nTracerCntl%lLinSrc) THEN
            IF(.NOT. lscat1) THEN
              CALL RayTrace(RayInfo, Core, phis1g, PhiAngIn1g, xst1g, tsrc, MocJout1g, iz, lJout,                   &
                            nTracerCntl%FastMocLv, nTracerCntl%lAFSS)
            ELSE
              CALL SetRtP1Src(Core, Fxr(:, iz), srcm, phim, xst1g, iz, ig,                                          &
                              ng, GroupInfo, l3dim, lXsLib, lscat1, nTracerCntl%ScatOd, PE)
              IF (nTracerCntl%lMultigrid) THEN
                CALL RayTraceP1_Multigrid(RayInfo, Core, phis1g, phim1g, PhiAngIn1g, xst1g, tsrc, srcm,             &
                                          MocJout1g, iz, nTracerCntl%ScatOd, lJout)
              ELSE
                CALL RayTraceP1(RayInfo, Core, phis1g, phim1g, PhiAngIn1g,                                          &
                                xst1g, tsrc, Srcm, MocJout1g, iz, lJout, nTracerCntl%ScatOd,                        &
                                nTracerCntl%FastMocLv, nTracerCntl%lAFSS)
              ENDIF
            ENDIF
          ELSE
            CALL LinPsiUpdate(Core, Fxr, LinPsi, LinSrcSlope, myzb, myze, ng, lxslib, GroupInfo)
            CALL SetRtLinSrc(Core, Fxr(:, iz), RayInfo, tsrc, LinSrc1g, LinPsi, LinSrcSlope,                        &
                             xst1g, eigv, iz, ig, ng, GroupInfo, l3dim,                                             &
                             lXslib, TRUE)
            CALL RayTraceLS(RayInfo, Core, phis1g, PhiAngIn1g, xst1g, tsrc, LinSrc1g,                               &
                            LinSrcSlope(:, :, iz, ig), MocJout1g, iz, lJout)
          ENDIF
          IF (lssph) THEN
            IF (ig.ge.igresb.and.ig.le.igrese) phis1g=phis1g*ssphf(:,iz,ig)
          ENDIF
          IF (RTMASTER) THEN
            IF (.NOT. nTracerCntl%lMOCUR) THEN
              phis(:, iz, ig) = phis1g
              IF (lScat1) phim(:, :, iz, ig) = phim1g
              FMInfo%PhiAngIn(:, :, iz, ig) = PhiAngIn1g
            ELSE
              CALL FluxUnderRelaxation(Core, Phis1g, Phis, FmInfo%w(ig), iz, ig, PE)
              CALL FluxInUnderRelaxation(Core, PhiANgIn1g, FmInfo%PhiAngIn, FmInfo%w(ig),                           &
                                         RayInfo%nPolarAngle, RayInfo%nPhiAngSv, iz, ig, PE)
              CALL CurrentUnderRelaxation(Core, MocJout1g, RadJout, FmInfo%w(ig),iz, ig, PE)
            ENDIF
          ENDIF
        ENDDO
        IF (lJout .AND. RTMASTER .AND. .NOT. nTracerCntl%lMocUR) RadJout(:, :, :, iz, ig) = MocJout1g
      ENDDO ! plane sweep
      TimeRt1gEnd = nTracer_dclock(FALSE, FALSE)
      TimeNodeElapsed = TimeNodeElapsed + TimeRt1gEnd - TimeRt1gBeg
      TimeRtElapsed = TimeRtElapsed + TimeRt1gEnd - TimeRt1gBeg
      IF(lDmesg .or. MOD(iG, 10) .eq. 0 .OR. ig .EQ. nG) THEN
        CALL MPI_MAX_REAL(TimeRtElapsed, PE%MPI_RTMASTER_COMM, TRUE)
        write(mesg,'(10x, a, i4, 2x, a, f10.3, 2x, a)') 'Group ', ig, ' finished in ', TimeRtElapsed, 'Sec'
        IF(master) CALL message(io8, FALSE, TRUE, mesg)
        TimeRtElapsed = 0.
      ENDIF
    ENDDO !--- END OF group sweep
  ENDDO !--- END OF Upscatter sweep
#ifdef __INTEL_COMPILER
  IF (nTracerCntl%lNodeTime) THEN
    i = hostnm(hostname)
    PRINT *, trim(hostname), TimeNodeElapsed
  ENDIF
#endif
#ifdef __PGI
#ifdef __linux
  IF (nTracerCntl%lNodeTime) THEN
    i = hostnm(hostname)
    PRINT *, trim(hostname), TimeNodeElapsed
  ENDIF
#endif
#endif

ELSE   !--- CNJ Edit : Node Majors

itrcntl%mocit = itrcntl%mocit + 1; TimeRtngElapsed = 0.
WRITE(mesg, '(a22,I5,a3)') 'Performing Ray Tracing', itrcntl%mocit, '...'
IF(Master) CALL message(io8, TRUE, TRUE, mesg)
IF(RTMaster) CALL FxrChiGen(Core, Fxr, FmInfo, GroupInfo, PE, nTracerCntl)
IF(nTracerCntl%lLinSrcCASMO .AND. RTMaster) THEN
  CALL LinPsiUpdate_CASMO(Core, Fxr, phisSlope, psiSlope, myzb, myze, ng, lxslib, GroupInfo)
END IF
DO iz = myzb, myze
  IF (.NOT. Core%lFuelPlane(iz) .AND. nTracerCntl%lAxRefFDM) CYCLE
  IF (RTMASTER) THEN
    DO ig = 1, ng
      phisnm(ig, :) = phis(:, iz, ig)
      PhiAngInnm(:, ig, :) = FmInfo%PhiAngIn(:, :, iz, ig)
    ENDDO
    IF (lScat1) phimnm => phim(:, :, :, iz)
    IF (nTracerCntl%lDomainDcmp) DcmpPhiAngIn => FMInfo%AsyPhiAngIn(:, :, :, :, :, iz)
    CALL SetRtMacXsNM(Core, Fxr(:, iz), xstnm, iz, ng, lxslib, lTrCorrection, lRST, lssph, lssphreg, PE)
#ifdef LkgSplit
    CALL PseudoAbsorptionNM(Core, Fxr(:, iz), AxPXS, xstnm, iz, ng, GroupInfo, l3dim)
#endif
#ifdef Buckling
    IF (nTracerCntl%lBsq) CALL AddBucklingNM(Core, Fxr, xstnm, nTracerCntl%Bsq, iz, ng, lxslib, lRST)
#endif
  ENDIF
  DO jsweep = 1, nGroupInfo
    GrpBeg = 1; GrpEnd = ng
    IF(jsweep .GT. 1) THEN
      GrpBeg = GroupInfo%UpScatRange(1); GrpEnd = GroupInfo%UpScatRange(2)
    ENDIF
    TimeRtngBeg = nTracer_dclock(FALSE, FALSE)
    DO InIter = 1, nInIter
      ljout = FALSE
      IF(InIter .EQ. nInIter) ljout = TRUE
      IF (RTMASTER) THEN
        IF (.NOT. nTracerCntl%lLinSrcCASMO) THEN
          CALL SetRtSrcNM(Core, Fxr(:, iz), srcnm, phisnm, psi, AxSrc, xstnm, eigv, iz,                             &
                          GrpBeg, GrpEnd, ng, GroupInfo, l3dim, lXslib, lscat1, FALSE, PE)
          IF (lScat1) THEN
            CALL SetRtP1SrcNM(Core, Fxr(:, iz), srcmnm, phimnm, xstnm, iz, GrpBeg, GrpEnd,                          &
                              ng, GroupInfo, lXsLib, nTracerCntl%ScatOd, PE)
          ENDIF
        ELSE
          CALL SetRtLinSrc_CASMO(Core, Fxr, RayInfo, phisnm, phisSlope, srcnm, srcSlope, psi, psiSlope,             &
                                 AxSrc, xstnm, eigv, iz, GrpBeg, GrpEnd, ng, GroupInfo,                             &
                                 l3dim, lxslib, lscat1, FALSE)
        ENDIF
      ENDIF
      IF (.NOT. nTracerCntl%lDomainDcmp) THEN
        IF (.NOT. nTracerCntl%lLinSrcCASMO) THEN
          IF (.NOT. lscat1) THEN
            CALL RayTraceNM_OMP(RayInfo, Core, phisnm, PhiAngInnm, xstnm, srcnm, MocJoutnm,                         &
                                iz, GrpBeg, GrpEnd, ljout, nTracerCntl%lDomainDcmp, nTracerCntl%FastMocLv)
          ENDIF
        ELSE
          CALL RayTraceLS_CASMO(RayInfo, Core, phisnm, phisSlope, PhiAngInnm, srcnm, srcSlope,                      &
                                xstnm, MocJoutnm, iz, GrpBeg, GrpEnd, lJout, nTracerCntl%lDomainDcmp)
        ENDIF
      ELSE
        CALL RayTrace_Dcmp(RayInfo, Core, iz, GrpBeg, GrpEnd, lJout, FALSE, lScat1, nTracerCntl%lLinSrcCASMO,       &
                           nTracerCntl%lHybrid)
      ENDIF
      IF (nTracerCntl%lsSPH) THEN
        IF (GrpBeg .LE. igrese) THEN
          ib = max(igresb, GrpBeg); ie = min(igrese, GrpEnd)
          phisnm(ib : ie, :) = phisnm(ib : ie, :) * ssphfnm(ib : ie, :, iz)
        ENDIF
      ENDIF
    ENDDO
    TimeRtngEnd = nTracer_dclock(FALSE, FALSE)
    TimeRtngElapsed(jsweep) = TimeRtngElapsed(jsweep) + (TimeRtngEnd - TimeRtngBeg)
  ENDDO
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) THEN
    CALL DcmpGatherCurrent(Core, MocJoutnm)
  ENDIF
#endif
  IF (RTMASTER) THEN
    DO ig = 1, ng
      phis(:, iz, ig) = phisnm(ig, :)
      FmInfo%PhiAngIn(:, :, iz, ig) = PhiAngInnm(:, ig, :)
      RadJout(:, :, :, iz, ig) = MocJoutnm(:, ig, :, :)
    ENDDO
  ENDIF
ENDDO
#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif

DO jsweep = 1, nGroupInfo
  TimeRtElapsed = TimeRtngElapsed(jsweep)
  CALL MPI_MAX_REAL(TimeRtElapsed, PE%MPI_NTRACER_COMM, TRUE)
  GrpBeg = 1; GrpEnd = ng
  IF (jsweep .GT. 1) THEN
    GrpBeg = GroupInfo%UpScatRange(1); GrpEnd = GroupInfo%UpScatRange(2)
  ENDIF
  write(mesg,'(10x, a, i4, 2x, a, i4, 2x, a, f10.3, 2x, a)') 'Group ', GrpBeg, ' to ', GrpEnd, ' finished in ', TimeRtElapsed, 'Sec'
  IF(master) CALL message(io8, FALSE, TRUE, mesg)
ENDDO
  
ENDIF

ENDDO

!Update Fission Source
IF(RTMASTER) THEN
  psid = psi; psicd = psic
  CALL PsiUpdate(Core, Fxr, phis, psi, myzb, myze, ng, lxslib, GroupInfo)
  CALL CellPsiUpdate(Core, Psi, psic, myzb, myze)
ENDIF

CALL UpdateEigv(Core, psi, psid, eigv, peigv, myzb, myze, PE)
IF(RTmaster) THEN
  fiserr = PsiErr(Core, Psi, PsiD, myzb, myze, PE)
  eigerr = abs((eigv - peigv))/eigv
  reserr = MocResidual(Core, FmInfo, eigv, GroupInfo, ng, PE, nTracerCntl)
ENDIF
#ifdef MPI_ENV
errdat = (/fiserr, eigerr, reserr/)
CALL BCAST(errdat, 3, PE%MPI_RT_COMM)
fiserr = errdat(1); eigerr = errdat(2); reserr = errdat(3)
#endif
write(mesg ,'(A5,I7,F15.6, 3(1pE12.3))')  'RT',itrcntl%mocit, EIGV, eigerr, fiserr, reserr
IF(MASTER) CALL message(io8, TRUE, TRUE, mesg)
ItrCntl%MocItrCntl%ResErr = ResErr

!Power Update
IF(RTMASTER) CALL PowerUpdate(Core, Fxr, phis, FmInfo%Power, myzb, myze, ng, lxslib, GroupInfo, PE)

IF(fiserr .lt. psiconv .and. eigerr .lt. eigconv .and. reserr .lt. resconv) THEN
  itrcntl%lconv = TRUE
  IF(nTracerCntl%lFeedBack) THEN
    last_TH = .true.
	ENDIF
ENDIF

MocTimeEnd = nTracer_dclock(FALSE, FALSE)
TimeRtElapsed = MocTimeEnd - MocTimeBeg
TimeChk%MocTime = TimeChk%MocTime + TimeRtElapsed

!Free the pointing
NULLIFY(FXR);
NULLIFY(PHIS); NULLIFY(PHIM)
NULLIFY(PSI); NULLIFY(PSID); NULLIFY(PSIC)
NULLIFY(PSICD); !NULLIFY(PhiAngin)
NULLIFY(RadJout); NULLIFY(AxSrc)

END SUBROUTINE