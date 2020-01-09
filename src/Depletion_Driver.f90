#include <defines.h>
  SUBROUTINE Depletion_Driver()
USE PARAM
USE GEOM,          ONLY : Core, ng
USE RAYS,          ONLY : RayInfo
USE Core_mod,      ONLY : CMInfo,            FmInfo,               THInfo,        &
                          GroupInfo,         eigv,                 peigv,         &
                          xkconv,            xkcrit
USE DeplLib_Mod,   ONLY : DeplLib 
USE Depl_Mod,      ONLY : DeplVars,          DeplCntl,                            &
                          Init_Depl,         AllocDeplFxrMem,      GetCoreHmMass, &
                          SetBurnUpStepInfo, DepletionStep,        SetLocalBurnup,&
                          SetDeplCoreState,  SaveCoreState,        EffBurnUpTime, &
                          SetDeplTimeStep,                                        &
                          UpdtFluxSolution,  DeplBanner
USE SubGrp_mod,    ONLY : SubGrpEffXsGen
USE Th_Mod,        ONLY : SteadyStateTH,     THFeedBack,          ThVar
USE Boron_mod,     ONLY : CalB10XS,          SaveBoronPPM,        PostB10Depl,   &
                          PostB10Depl0
USE CritSpec_mod,  ONLY : XsKinf
USE XeDyn_Mod,     ONLY : SetDeplXeDynMode,  UpdtXeDyn,                           &
                          InitXe,            FinalXe,              UpdtXe            
USE PE_MOD,        ONLY : PE
USE CNTL,          ONLY : nTracerCntl
USE ItrCNTL_mod,   ONLY : ItrCntl
USE FILES,         ONLY : io8
USE IOUTIL,        ONLY : message,           ShowHbar1,            ShowHbar2
USE Timer,         ONLY : nTracer_dclock, TimeChk
USE Restart_mod,   ONLY : WriteRestartFile
#ifdef MPI_ENV
USE MPIComm_mod,   ONLY : MPI_SYNC
#endif
IMPLICIT NONE

REAL :: PrevBurnUp(2)
REAL :: ppm
REAL :: efftime, delpow
REAL :: Tbeg, Tend
INTEGER :: nBurnUpStep
INTEGER :: iBurnUp, itype

LOGICAL :: Master, Slave
LOGICAL, SAVE :: lFirst


Data lFirst /.TRUE./


Master = PE%Master; Slave = PE%Slave

WRITE(mesg, '(A)') 'Performing Depletion ...'
IF(Master) CALL message(io8, TRUE, TRUE, mesg)

IF(lFirst) THEN
  IF(.NOT. DeplCntl%lInitDepl) THEN
    DeplCntl%lInitDepl = .TRUE.
    CALL Init_Depl(DeplLib, DeplVars, GroupInfo, PE)
    CALL AllocDeplFxrMem(Core, FmInfo%Fxr, GroupInfo, PE)
  ENDIF
  CALL DeplBanner(1)
  lFirst = .FALSE.
ENDIF


itype = DeplCntl%BurnUpType
nBurnUpStep = DeplCntl%nBurnUpStep
ppm = 0

!Set Xenon Dynamics Mode
CALL SetDeplXeDynMode(Core, FmInfo, DeplCntl, nTracerCntl, PE)

CALL SetBurnUpStepInfo(nTracerCntl%PowerCore, DeplVars(1), DeplCntl)
IF(DeplCntl%lCoreFollow) CALL InitCoreFollow(DeplCntl, nTracerCntl)

!
IF(DeplCntl%NowBurnUp(1) .LE. epsm8) THEN
  !Core%DelT = DeplCntl%TSec
  CALL SetDeplTimeStep(0, Core, DeplCntl, nTracerCntl ,PE)
  !Initialize Xenon Dynamics
  CALL UpdtXeDyn(Core, FmInfo, THInfo, GroupInfo, nTracerCntl, DeplCntl, PE, Mode = InitXe)
  !Solve Depletion Chain 
  CALL DepletionStep(Core, FmInfo, DeplVars, DeplLib, ThInfo, GroupInfo, nTracerCntl, DeplCntl, PE)
  !Effective Cross Section Updates 
  CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, Eigv, GroupInfo, nTracerCntl, PE)
  !Update Xe Dynamics 
  CALL UpdtXeDyn(Core, FmInfo, THInfo, GroupInfo, nTracerCntl, DeplCntl, PE, Mode = FinalXe)
ENDIF

DO iBurnUp = 1, nBurnUpStep
  nTracerCntl%CalNoId = nTracerCntl%CalNoId + 1
  IF(DeplCntl%lCoreFollow) THEN
    CALL SetDeplCoreState(iburnup, DeplCntl%CoreState, ThInfo, ThVar, nTracerCntl, PE)
    IF(DeplCntl%CoreState%lStateChg(iburnup)) CALL UpdtFluxSolution(iburnup, DeplCntl, nTracerCntl, PE)
  ENDIF
  !Burnup Time 
  CALL SetDeplTimeStep(iburnup, Core,  DeplCntl, nTracerCntl, PE)
  
  WRITE(mesg, '(a)') 'Predictor Step...'  
  IF(Master) CALL message(io8, TRUE, TRUE, mesg) 
  
  !Time Check
  Tbeg = nTracer_dclock(FALSE, FALSE)
  
  
  !--- BYS edit 15/12/30 --- depletion test
  
  !Prediction Step
  DeplCntl%lPredict = .TRUE.
  !Effective Cross Section
  !CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, Eigv, GroupInfo, nTracerCntl, PE)
  !Depletion Matrix Calculation
  CALL DepletionStep(Core, FmInfo, DeplVars, DeplLib, ThInfo, GroupInfo, nTracerCntl, DeplCntl, PE)
  !Time Check
  Tend = nTracer_dclock(FALSE, FALSE)
  TimeChk%DeplTime = TimeChk%DeplTime + Tend - Tbeg
  
  !IF(DeplCntl%B10DeplMod .EQ. 1 .AND. nTracerCntl%BoronPPM .GT. 0) CALL OnlineB10Depl(iBurnUp, .TRUE., DeplCntl, PE)
  
  !Flux Calculation
  CALL SSEIG
  IF (nTracerCntl%OutpCntl%nDeplRegXsOut.gt.0) CALL ProcessEffXsDepl(Core, FmInfo, CmInfo, ThInfo, GroupInfo, nTracerCntl, DeplCntl,PE,iburnup)
  
  IF(Master) CALL ShowHbar1(io8)
  
  !Correction Step
  WRITE(mesg, '(a)') 'Corrector Step'
  IF(Master) CALL message(io8, TRUE, TRUE, mesg)  
  DeplCntl%lPredict = .FALSE.
  
  Tbeg = nTracer_dclock(FALSE, FALSE)
  CALL DepletionStep(Core, FmInfo, DeplVars, DeplLib, ThInfo, GroupInfo, nTracerCntl, DeplCntl, PE)
  Tend = nTracer_dclock(FALSE, FALSE)
  TimeChk%DeplTime = TimeChk%DeplTime + Tend - Tbeg
  !Changhyun Test - Number Density Out for Each Fxr
  !CALL Write_ND_Debug(Core, FmInfo, DeplCntl, PE)
  !Effective Cross Section
  CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, Eigv, GroupInfo, nTracerCntl, PE)
  CALL UpdtXeDyn(Core, FmInfo, THInfo, GroupInfo, nTracerCntl, DeplCntl, PE, Mode = FinalXe)

  !--- BYS edit END 15/12/30 --- depletion test
  
  !2nd Flux Calculation
  IF(DeplCntl%PC_OPT .EQ. 2) CALL SSEIG() !Changhyun Test
  IF(DeplCntl%lCoreFollow .OR. DeplCntl%PC_OPT .EQ. 2) CALL SaveCoreState(iburnup, DeplCntl%CoreState, nTracerCntl, PE)
  
  xkcrit = XsKinf(Core, FmInfo, GroupInfo, .TRUE., nTracerCntl%lCritSpec, PE)
  CALL DeplBanner(2)
  IF(PE%lCmfdGrp) CALL OutputEdit()
  IF(nTRACERCntl%lWriteDeplRst) CALL  WriteRestartFile(Core, FmInfo, .TRUE., DeplCntl%NowStep, nTracerCntl, PE)
  IF(MASTER) CALL ShowHbar2(io8)
ENDDO

!IF(DeplCntl%B10DeplMod .EQ. 2) CALL PostB10Depl(DeplCntl, PE) 
END SUBROUTINE
  

  !--- BYS edit 15/12/30 --- depletion test / back-up
  
  !!Prediction Step
  !DeplCntl%lPredict = .TRUE.
  !!Effective Cross Section
  !CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, Eigv, GroupInfo, nTracerCntl, PE)
  !!Depletion Matrix Calculation
  !CALL DepletionStep(Core, FmInfo, DeplVars, DeplLib, ThInfo, GroupInfo, nTracerCntl, DeplCntl, PE)
  !!Time Check
  !Tend = nTracer_dclock(FALSE, FALSE)
  !TimeChk%DeplTime = TimeChk%DeplTime + Tend - Tbeg
  !
  !!IF(DeplCntl%B10DeplMod .EQ. 1 .AND. nTracerCntl%BoronPPM .GT. 0) CALL OnlineB10Depl(iBurnUp, .TRUE., DeplCntl, PE)
  !
  !!Flux Calculation
  !CALL SSEIG
  !
  !IF(Master) CALL ShowHbar1(io8)
  !
  !!Correction Step
  !WRITE(mesg, '(a)') 'Corrector Step'
  !IF(Master) CALL message(io8, TRUE, TRUE, mesg)  
  !DeplCntl%lPredict = .FALSE.
  !
  !Tbeg = nTracer_dclock(FALSE, FALSE)
  !CALL DepletionStep(Core, FmInfo, DeplVars, DeplLib, ThInfo, GroupInfo, nTracerCntl, DeplCntl, PE)
  !CALL UpdtXeDyn(Core, FmInfo, THInfo, GroupInfo, nTracerCntl, DeplCntl, PE, Mode = FinalXe)
  !Tend = nTracer_dclock(FALSE, FALSE)
  !TimeChk%DeplTime = TimeChk%DeplTime + Tend - Tbeg
  !
  !
  !--- BYS edit END 15/12/30 --- depletion test
  