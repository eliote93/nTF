#include <defines.h>

SUBROUTINE Transient_Driver()
USE PARAM
USE GEOM,          ONLY : Core, ng
USE Core_mod,      ONLY : CMInfo,            FmInfo,               THInfo,        &
                          GroupInfo,         eigv
USE RAYS,          ONLY : RayInfo
USE PE_MOD,        ONLY : PE
USE CNTL,          ONLY : nTracerCntl
USE ItrCNTL_mod,   ONLY : ItrCntl
USE FILES,         ONLY : io8
USE IOUTIL,        ONLY : message
USE Timer,         ONLY : nTracer_dclock, TimeChk
#ifdef MPI_ENV
USE MPIComm_mod,   ONLY : MPI_SYNC
#endif

USE TRAN_MOD,      ONLY : TranInfo,          TranCntl,                            &
                          InitTransient,     InitPrecursor,       UpdtPrec,       &
                          KinParamGen,       UpdtResSrc,          SaveTranSol,    &
                          UpdtExpTrsf,       SolExpExtpltn,                       &
                          TranReactivityUpdt
USE TH_Mod,        ONLY : TransientTH,       THFeedback
USE SubGrp_mod,    ONLY : SubGrpEffXsGen, SubGrpFsp
USE XsPerturb_mod, ONLY : XsPerturbation
USE Boron_mod,      ONLY : SetBoronCoolant
IMPLICIT NONE

INTEGER :: istep
LOGICAL :: lout

CALL SSEIG()
CALL OutputEdit()

WRITE(mesg, '(A)') hbar2(1:77)
CALL message(io8, FALSE, TRUE, MESG)
!CALL NULLTransient_Driver()
!RETURN

CALL InitTransient(Core, RayInfo, FmInfo, CmInfo, THInfo, eigv, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)
CALL KinParamGen(Core, FmInfo, TranInfo,  ThInfo, GroupInfo, .TRUE., nTracerCntl, PE)
CALL InitPrecursor(Core, FmInfo, CmInfo, TranInfo, GroupInfo, nTracerCntl, PE)
TranCntl%lCondiMOC = .TRUE.
IF(TranCntl%lCondiMOC) THEN
  CALL CheckCondMOC(Core, FmInfo, CmInfo, TranInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
ENDIF


TranCntl%NowStep = 1
CALL UpdtResSrc(core, FmInfo, CmInfo, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)
TranCntl%lExpTrsf = .FALSE.
!CALL PrintHeatPower(Core, CmInfo, FmInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
DO istep = 1, TranCntl%nStep
  nTracerCntl%CalNoId = nTracerCntl%CalNoId + 1
  CALL XsPerturbation(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, PE)
  CALL SetBoronCoolant(Core, FmInfo%Fxr, nTracerCntl%BoronPPM , PE%myzb, PE%myze)
  CALL KinParamGen(Core, FmInfo, TranInfo, ThInfo, GroupInfo, .TRUE., nTracerCntl, PE)

  !Solution Exponential Extrapolation
  !IF(TranCntl%lExpTrsf) CALL SolExpExtpltn(Core, FmInfo, CmInfo, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)
  IF(nTracerCntl%lFeedback) THEN
    CALL TransientTH(Core, CmInfo, FmInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
    CALL THFeedBack(Core, CmInfo, FmInfo, ThInfo, nTRACERCntl, PE)
    CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, TranInfo%Eigv0, GroupInfo, nTracerCntl, PE)
  ENDIF

  IF(TranCntl%lCondiMOC) THEN
    CALL CheckCondMOC(Core, FmInfo, CmInfo, TranInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
  ENDIF

  IF(nTracerCntl%lFeedback .AND. TranCntl%lSgFspUpdt) THEN
    CALL SubGrpFsp(Core, FmInfo%Fxr, THInfo, RayInfo, GroupInfo, nTracerCntl, PE)
  ENDIF

  CALL TransientFsp_Driver(Core, FmInfo, CmInfo, TranInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)


  !Thetha method source term
  CALL UpdtResSrc(core, FmInfo, CmInfo, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)
  !Precursor Update
  CALL UpdtPrec(Core, FmInfo, CmInfo, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)
  !Prepare Exponental Transformation
  IF(TranCntl%lExpTrsf) CALL UpdtExpTrsf(Core, FmInfo, CmInfo, TranInfo, GroupInfo, .FALSE., TranCntl, nTracerCntl, PE)
  !Calculate
  CALL UpdtPowerLevel(Core, FmInfo, CmInfo, TranInfo, ng, nTracerCntl, PE)
  nTracerCntl%PowerLevel = TranInfo%PowerLevel; ThInfo%PowLv = nTracerCntl%PowerLevel
  TranInfo%reactivity = TranReactivityUpdt(CmInfo, eigv, TranCntl, PE)
  !Save solution of previous time step
  CALL SaveTranSol(Core, FmInfo, CmInfo, TranInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)

  CALL ChkOutputWrite(lout, TranCntl)
  IF(lout) THEN
    CALL OutputEdit()
    !CALL PrintHeatPower(Core, CmInfo, FmInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
    TranCntl%NowWriteStep = TranCntl%NowWriteStep + 1
  ENDIF
  IF(PE%MASTER) THEN
    WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, 1e10.3, F10.5)') '@', TranCntl%T(istep)," s, P=", TranInfo%PowerLevel * 100, TranInfo%reactivity
    CALL message(io8, TRUE, TRUE, MESG)
    WRITE(89, '(a)') mesg
    WRITE(mesg, '(A)') hbar2(1:77)
    CALL message(io8, FALSE, TRUE, MESG)
  ENDIF
  CALL AdaptiveTheta(TranInfo, TranCntl)
  TranCntl%NowStep = TranCntl%NowStep + 1
  TranInfo%PowerLeveld = TranInfo%PowerLevel
ENDDO
END SUBROUTINE

SUBROUTINE TransientFsp_Driver(Core, FmInfo, CmInfo, TranInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,           ONLY : CoreInfo_Type,           FmInfo_Type,           CmInfo_Type,         &
                              TranInfo_Type,           GroupInfo_Type,        TranCntl_Type,       &
                              ThInfo_Type,             PE_Type
USE CNTL,              ONLY : nTracerCntl_Type
USE SubGrp_mod,        ONLY : SubGrpEffXsGen,          SubGrpFsp
USE TRAN_MOD,          ONLY : UpdtPowerLevel,          CmFmInfoSync,           UpdtExpTrsf,        &
                              CheckCondMOC,      UpdtBaseXsCondiMOC
USE TRANMOC_MOD,       ONLY : TranMOC_Driver
USE TRANCMFD_MOD,      ONLY : TranCMFD_Driver
USE Boron_mod,         ONLY : UpdtBoronPPM,           SetBoronCoolant
USE RAYS,              ONLY : RayInfo
USE ItrCNTL_mod,       ONLY : ItrCntl
USE FILES,             ONLY : io8
USE IOUTIL,            ONLY : message
USE Timer,             ONLY : nTracer_dclock, TimeChk
USE TH_Mod,            ONLY : TransientTH,       THFeedback

USE SubGrp_mod,        ONLY : SubGrpEffXsGen, SubGrpFsp
USE BasicOperation,    ONLY : MULTI_CA
#ifdef MPI_ENV
USE MPIComm_mod,       ONLY : MPI_SYNC
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE

INTEGER :: iter
LOGICAL :: Master, Slave, RTMaster, RTSlave, CMFDMaster, CMFDslave, lCmfdGrp
LOGICAL :: lCMFD, lMOC

Master = PE%Master; Slave = PE%Slave
RTMaster = PE%RTMaster; RTSlave = PE%RTSlave
CMFDMaster = PE%CMFDMaster; CMFDSlave = PE%CMFDSlave
lCmfdGrp = PE%lCmfdGrp

#ifdef MPI_ENV
nTracerCntl%lcmfd = .TRUE.
CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif

lCMFD = nTracerCntl%lcmfd
lMOC = .TRUE.
ItrCntl%lconv = .FALSE.

IF(RTmaster) CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, TranInfo%Eigv0, GroupInfo, nTracerCntl, PE)
#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif

CALL  InitIterVar(Core, FMInfo, CmInfo, GroupInfo,  .FALSE., ItrCntl, nTracerCntl, PE)

!TranCntl%lMOCUpdt = .FALSE.
TranCntl%lAxNUpdt = .TRUE.
lMOC = TranCntl%lMOCUpdt
CALL TranCmfd_Driver(Core, FmInfo, CmInfo, TranInfo, THInfo, GroupInfo, TranCntl, nTracerCntl, ItrCntl, PE)

DO iter = 1, 8
  !TranCntl%lMOCUpdt = .TRUE.
  ItrCntl%SrcIt = ItrCntl%SrcIt + 1
  IF(Iter .GT. 3) lMOC = .FALSE.
  !IF(Iter .GT. 4) TranCntl%lAxNUpdt = .FALSE.
  IF(lMOC) CALL TranMOC_Driver(Core, RayInfo, FmInfo, CmInfo, TranInfo, THInfo, GroupInfo, TranCntl, nTracerCntl, ItrCntl, PE)
  !IF(nTracerCntl%lFeedback) THEN
  !  CALL UpdtPowerLevel(Core, FmInfo, CmInfo, TranInfo, GroupInfo%ng, nTracerCntl, PE)
  !  nTracerCntl%PowerLevel = TranInfo%PowerLevel; ThInfo%PowLv = nTracerCntl%PowerLevel
  !  CALL TransientTH(Core, CmInfo, FmInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
  !  CALL THFeedBack(Core, CmInfo, FmInfo, ThInfo, nTRACERCntl, PE)
  !  CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, TranInfo%Eigv0, GroupInfo, nTracerCntl, PE)
  !ENDIF
  CALL TranCmfd_Driver(Core, FmInfo, CmInfo, TranInfo, THInfo, GroupInfo, TranCntl, nTracerCntl, ItrCntl, PE)
  IF(iter .LT. 4) CYCLE
  IF(ItrCntl%lconv) EXIT
  IF(MOD(iter, 2) .EQ. 1 .AND. TranCntl%lExpTrsf) CALL UpdtExpTrsf(Core, FmInfo, CmInfo, TranInfo, GroupInfo, .FALSE., TranCntl, nTracerCntl, PE)
ENDDO
IF(TranCntl%lMocUpdt .AND. TranCntl%lCondiMOC) THEN
  CALL UpdtBaseXsCondiMOC(Core, FmInfo, CmInfo, TranInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
ENDIF
CALL CmFmInfoSync(Core, FmInfo, CmInfo, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)

END SUBROUTINE

