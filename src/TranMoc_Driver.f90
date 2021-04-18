#include <defines.h>
SUBROUTINE TranMOC_Driver(Core, RayInfo, FmInfo, CmInfo, TranInfo, THInfo, GroupInfo, TranCntl, nTracerCntl, ItrCntl, PE)
USE PARAM
USE TYPEDEF,           ONLY : CoreInfo_Type,           FmInfo_Type,           CmInfo_Type,         &
                              TranInfo_Type,           GroupInfo_Type,        TranCntl_Type,       &
                              ThInfo_Type,             PE_Type,               RayInfo_Type,        &
                              FxrInfo_Type
USE CNTL,              ONLY : nTracerCntl_Type
USE itrcntl_mod,       ONLY : ItrCntl_TYPE, CMFDItrCntl_TYPE
USE SUBGRP_MOD,        ONLY : FxrChiGen
USE BasicOperation,    ONLY : CP_CA,                   CP_VA,                 AD_VA,                &
                              MULTI_CA
USE IOUTIL,            ONLY : message
USE TRANMOC_MOD,       ONLY : TrSrc,                   PrecSrc,                                     &
                              SetPrecParam,            SetTranMOCEnv,         &
                              PrecSrcUpdt,             SetTranSrc,           TranMocResidualError,  &
                              SetExpTrsfXs
USE MOC_MOD,           ONLY : RayTraceGM_OMP,                SetRtMacXsGM_Cusping,   SetRtSrcGM_Cusping,      &
                              RayTraceLS,              PsiUpdate,            CellPsiUpdate,         &
                              MocResidual,             PsiErr,               PseudoAbsorptionGM,      &
                              PowerUpdate,                                                          &
                              phis1g,                  MocJout1g,            xst1g,                 &
                              tSrc,                    AxSrc1g,              PhiAngin1g,            &
                              srcm
USE Timer,             ONLY : nTracer_dclock, TimeChk
USE FILES,             ONLY : io8

#ifdef MPI_ENV
USE MPICOMM_MOD,       ONLY : BCAST,              MPI_SYNC,           MPI_MAX_REAL
#endif
USE XSLIB_MOD,   ONLY : igresb,igrese
USE SPH_mod,     ONLY : ssphf,ssphfnm,calcPinSSPH
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

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
REAL, POINTER :: PHIS(:,:,:)                          !Scalar Flux Vectors
REAL, POINTER :: PSI(:,:), PSID(:,:)                  !Fsr FIssion Scource
REAL, POINTER :: PSIC(:,:), PSICD(:,:)                !Cell Fission Source Term
REAL, POINTER :: AxSrc(:, :, :), AxPXS(:, :, :)
REAL, POINTER :: RadJout(:, :, :, :, :)

REAL, POINTER :: Prec(:, :, :), ResSrc(:, :, :)
REAL, POINTER :: TranPsid(:, :), TranPsi(:, :), TranPhi(:, :, :)

INTEGER :: ig, iz, l, iter, jsweep

REAL :: Eigv0

INTEGER :: myzb, myze
INTEGER :: ng, nPrec
INTEGER :: nFSR, nxy, nbd, GrpBeg, GrpEnd, nGroupInfo
!Time Variables
REAL :: TimeRt1gBeg, TimeRt1gEnd, TimeElasped, MocTimeBeg, MocTimeEnd
!Error
REAL :: ResErr, FisErr, PsiConv, EigErr, errdat(3)

INTEGER :: InIter, nInIter

LOGICAL :: lJout, lxslib, l3dim, lscat1, lTrCorrection, lRST, lssph, lssphreg
LOGICAL :: lDmesg                  !Detailed Message
LOGICAL :: MASTER, SLAVE, RTMaster, RTSlave, CMFDMaster, CMFDSlave

LOGICAL, SAVE :: lFirst
DATA lfirst /TRUE/

FXR => FmInfo%Fxr
PHIS => FmInfo%PHIS; PSI => FmInfo%PSI
PSID => FmInfo%PSID; PSIC => FmInfo%PSIC
PSICD => FmInfo%PSICD; !PhiAngin => FmInfo%PhiAngin
RadJout => FmInfo%RadJout;
AxSrc => FmInfo%AxSrc; AxPXS => FmInfo%AxPXS
!Parallel Environment In
Master = PE%Master; Slave = PE%Slave
RTMaster = PE%RTMaster; RTSlave = PE%RTSlave
!Transient Variables
Prec => FmInfo%Prec; ResSrc => FmInfo%ResSrc
TranPhi => FmInfo%TranPhi
TranPsi => FmInfo%TranPsi; TranPsid => FmInfo%TranPsid
!
lXsLib = nTracerCntl%lXsLib;

!Iteration Control
psiconv = itrcntl%psiconv;

!Local Variable
Eigv0 = TranInfo%Eigv0

nbd = 4
myzb = PE%myzb; myze = PE%myze
nFsr = Core%nCoreFsr; nxy = Core%nxy
ng = GroupInfo%ng; nprec = GroupInfo%ng
iter = 0; nGroupInfo =2
nInIter = 2
!Control variables
lJout = TRUE
lxslib = nTracerCntl%lXslib; l3dim = nTracerCntl%l3dim
lscat1 = nTracerCntl%lscat1
lTrCorrection=nTracerCntl%lTrCorrection
lRST=nTracerCntl%lRST
lssph = nTracerCntl%lsSPH
lssphreg = nTracerCntl%lsSPHreg

IF(lFirst) THEN
  CALL SetTranMOCEnv(Core, TranInfo, PE)
  lFirst = .FALSE.
ENDIF


IF(RTMASTER) THEN
  CALL CP_VA(psid(1:nFsr, myzb:myze), psi(1:nFsr, myzb:myze), nFsr, myze - myzb +1)
  CALL CP_VA(psicd(1:nxy, myzb:myze), psic(1:nxy, myzb:myze), nxy,  myze - myzb +1)
  CALL PsiUpdate(Core, Fxr, phis, psi, myzb, myze, ng, lxslib, GroupInfo)
  CALL CellPsiUpdate(Core, Psi, psic, myzb, myze)
  CALL MULTI_CA(1._8 /eigv0, psi(1:nFsr, myzb:myze), nFsr, myze - myzb +1)
  CALL MULTI_CA(1._8 /eigv0, psic(1:nxy, myzb:myze), nxy, myze - myzb +1)
ENDIF

IF(RTMASTER) THEN
  !CALL KinParamGen(Core, FmInfo, ThInfo, TranInfo, GroupInfo, nTracerCntl, PE)
  CALL SetPrecParam(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, PE)
  CALL PrecSrcUpdt(Core, Fxr, PrecSrc, Prec, TranPsi, TranPsid, GroupInfo, TranCntl, nTRACERCntl, PE)
ENDIF

!Output Message Control
lDmesg = TRUE
IF(ng .gt. 10) lDmesg = FALSE
IF(.NOT. GroupInfo%lUpScat) nGroupInfo = 1

itrcntl%mocit = itrcntl%mocit + 1; TimeElasped = 0.
WRITE(mesg, '(a22,I5,a3)') 'Performing Ray Tracing for Transient', itrcntl%mocit, '...'
IF(Master) CALL message(io8, TRUE, TRUE, mesg)
IF(RTMaster) CALL FxrChiGen(Core, Fxr, FmInfo, GroupInfo, PE, nTracerCntl)
IF(nTracerCntl%lCusping_MPI) CALL GetNeighborMocFlux(phis, FmInfo%neighphis, nFsr, myzb, myze, 1, ng, Core%nz, Core%AxBC)
DO jsweep =1, nGroupInfo
  GrpBeg = 1; GrpEnd = ng
  !IF(.NOT. GroupInfo%lUpScat .AND. jsweep .GT. 1) EXIT
  IF(jsweep .GT. 1) THEN
    GrpBeg = GroupInfo%UpScatRange(1); GrpEnd = GroupInfo%UpScatRange(2)
  ENDIF

  DO ig = GrpBeg, GrpEnd
    TimeRt1gBeg = nTracer_dclock(FALSE, FALSE)
    !Set Transport XS
    DO iz = myzb, myze
      ljout = FALSE
      IF(RTMaster) CALL CP_VA(AxSrc1g(1:nxy), AxSrc(1:nxy, iz, ig), nxy)
      DO InIter = 1, nInIter
        IF(InIter .EQ. nInIter) ljout = TRUE
        IF(RTMASTER) THEN
          CALL SetRtMacXsGM_Cusping(Core, FmInfo, Fxr(:, iz), xst1g, phis, iz, ig, ng, lxslib, lTrCorrection, lRST, lssph, lssphreg, PE)
          CALL PseudoAbsorptionGM(Core, Fxr(:, iz), tsrc, phis(:, iz, ig),                             &
                                AxPXS(:, iz, ig), xst1g, iz, ig, ng, GroupInfo, l3dim)
          !CALL SetExpTrsfXs(Core, Fxr, xst1g, iz, ig, GroupInfo, TranInfo, TranCntl, nTracerCntl, PE)
          CALL SetRtSrcGM_Cusping(Core, FmInfo, Fxr(:, iz), tsrc, phis, psi, axSrc1g, xst1g,                           &
                        1._8, iz, ig, ng, GroupInfo, l3dim, lXslib, lscat1, FALSE, PE)
          CALL SetTranSrc(Core, Fxr, TrSrc, Phis, TranPhi, Psi, PrecSrc, ResSrc, xst1g,              &
                         iz, ig, GroupInfo, TranInfo, TranCntl, nTracerCntl, PE)
          CALL CP_VA(PhiAngin1g, FmInfo%PhiAngin(:,: ,iz, ig), RayInfo%nPolarAngle, RayInfo%nPhiAngSv)
          CALL AD_VA(TrSrc(1:nfsr), TrSrc(1:nfsr), tsrc(1:nfsr), nfsr)
        ENDIF
        
        CALL RayTraceGM_OMP(RayInfo, Core, phis1g, PhiAngIn1g, xst1g, trsrc, MocJout1g, iz, lJout)
        IF (lssph) THEN
          IF (ig.ge.igresb.and.ig.le.igrese) phis1g=phis1g*ssphf(:,iz,ig)
        ENDIF
        IF(RTMASTER) CALL CP_VA(phis(1:nFsr, iz, ig), phis1g(1:nFsr), nFsr)
        IF(RTMASTER) CALL CP_VA(FmInfo%PhiAngin(:,:,iz, ig), PhiAngin1g, RayInfo%nPolarAngle, RayInfo%nPhiAngSv)
        IF(lJout .AND. RTMASTER) CALL CP_VA(RadJout(1:3, 1:nbd, 1:nxy, iz, ig), MocJout1g(1:3, 1:nbd, 1:nxy), 3, nbd, nxy)   ! > 16/02/18
        !IF(lJout .AND. RTMASTER) CALL CP_VA(RadJout(1:2, 1:nbd, 1:nxy, iz, ig), MocJout1g(1:2, 1:nbd, 1:nxy), 2, nbd, nxy)
        IF(nTracerCntl%lCusping_MPI) THEN
          IF(iz .EQ. myzb .OR. iz .EQ. myze) CALL GetNeighborMocFlux(phis, FmInfo%neighphis, nFsr, myzb, myze, ig, ig, Core%nz, Core%AxBC)
        END IF
      ENDDO
    ENDDO
#ifdef MPI_ENV
    CALL MPI_SYNC(PE%MPI_RTMASTER_COMM)
#endif
    !Edit Print out Log messages
    TimeRt1gEnd = nTracer_dclock(FALSE, FALSE)
    TimeElasped = TimeElasped + TimeRt1gEnd - TimeRt1gBeg
    IF(lDmesg .or. MOD(iG, 10) .eq. 0 .OR. ig .EQ. nG) THEN
      CALL MPI_MAX_REAL(TimeElasped, PE%MPI_RTMASTER_COMM, TRUE)
      write(mesg,'(10x, a, i4, 2x, a, f10.3, 2x, a)') 'Group ', ig, ' finished in ', TimeElasped, 'Sec'
      IF(master) CALL message(io8, FALSE, TRUE, mesg)
      TimeElasped = 0.
    ENDIF
  ENDDO
ENDDO

!Update Fission Source
IF(RTMASTER) THEN
  CALL CP_VA(psid(1:nFsr, myzb:myze), psi(1:nFsr, myzb:myze), nFsr, myze - myzb +1)
  CALL CP_VA(psicd(1:nxy, myzb:myze), psic(1:nxy, myzb:myze), nxy,  myze - myzb +1)
  CALL PsiUpdate(Core, Fxr, phis, psi, myzb, myze, ng, lxslib, GroupInfo)
  CALL CellPsiUpdate(Core, Psi, psic, myzb, myze)
  CALL MULTI_CA(1._8 /eigv0, psi(1:nFsr, myzb:myze), nFsr, myze - myzb +1)
  CALL MULTI_CA(1._8 /eigv0, psic(1:nxy, myzb:myze), nxy, myze - myzb +1)
ENDIF

IF(RTMASTER) CALL PowerUpdate(Core, Fxr, phis, FmInfo%Power, myzb, myze, ng, lxslib, GroupInfo, PE)

eigerr = 0
fiserr = PsiErr(Core, Psi, PsiD, myzb, myze, PE)
reserr = TranMocResidualError(Core, FmInfo, TranInfo, eigv0, GroupInfo, TranCntl, nTRACERCntl, PE)
#ifdef MPI_ENV
errdat = (/fiserr, eigerr, reserr/)
CALL BCAST(errdat, 3, PE%MPI_RT_COMM)
fiserr = errdat(1); eigerr = errdat(2); reserr = errdat(3)
#endif
WRITE(mesg ,'(A5,I7,F15.6, 3(1pE12.3))')  'RT',itrcntl%mocit, EIGV0, eigerr, fiserr, reserr
IF(MASTER) CALL message(io8, TRUE, TRUE, mesg)
!IF(fiserr .lt. psiconv) itrcntl%lconv = TRUE
IF(fisErr .LT. TranCntl%psi_conv) ItrCntl%lConv = .TRUE.
END SUBROUTINE
