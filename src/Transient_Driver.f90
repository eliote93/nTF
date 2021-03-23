#include <defines.h>
SUBROUTINE Transient_Driver()
USE PARAM
USE GEOM,          ONLY : Core
USE Core_mod,      ONLY : FmInfo
USE PE_Mod,        ONLY : PE
USE files,         ONLY : io8
USE ioutil,        ONLY : message
USE timer,         ONLY : nTracer_dclock,     TimeChk
USE TRAN_MOD,      ONLY : TranInfo, TranCntl 
USE XsPerturb_mod, ONLY : XsPerturbation_T
USE CNTL,          ONLY : nTracerCntl
IMPLICIT NONE

#ifdef __PGI
IF (PE%lCUDACMFD) THEN
  !CALL XsPerturbation(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, PE, 0., -1.)
  CALL SSEIG()
  CALL OutputEdit()
  IF(TranCntl%lfixtmprw) THEN 
    CALL XsPerturbation_T(Core, FmInfo, TranInfo, Trancntl, nTracerCntl, PE, 0.1, 0.)
    nTracerCntl%lFeedBack = .FALSE.
    CALL SSEIG()
    RETURN
  END IF

  WRITE(mesg, '(A)') hbar2(1:77)
  IF(PE%MASTER) CALL message(io8, FALSE, TRUE, MESG)
  IF(TranCntl%lAdptT) THEN
    CALL CUDATransient_Adpt_Driver()
  ELSE
    CALL CUDATransient_Driver()
  END IF
  GOTO 4
ENDIF
#endif
#ifdef __INTEL_MKL
!IF (PE%lMKL) THEN
!  CALL Transient_Driver_Parallel()
!  GOTO 4
!ENDIF
#endif
CALL Transient_Driver_old()
4 CONTINUE

END SUBROUTINE

SUBROUTINE Transient_Driver_old()
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
USE PromptFeedback_mod, ONLY : UpdtFuelTemp
USE TranMacXsLib_Mod, ONLY : InitChidLib
IMPLICIT NONE

INTEGER :: istep
LOGICAL :: lout

CALL SSEIG()
CALL OutputEdit()

WRITE(mesg, '(A)') hbar2(1:77)
CALL message(io8, FALSE, TRUE, MESG)
!CALL NULLTransient_Driver()
!RETURN

WRITE(mesg, '(a)') 'Initialize Transient Calculation...'
IF (PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)
CALL InitTransient(Core, RayInfo, FmInfo, CmInfo, THInfo, eigv, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)
CALL KinParamGen(Core, FmInfo, TranInfo,  ThInfo, GroupInfo, .TRUE., nTracerCntl, PE)
IF(nTracerCntl%lchidgen) CALL InitChidLib(Core, FmInfo, GroupInfo, PE, nTracerCntl, TranInfo%Chid, ng)
CALL InitPrecursor(Core, FmInfo, CmInfo, TranInfo, GroupInfo, nTracerCntl, PE)

IF(TranCntl%lCondiMOC) THEN
  CALL CheckCondMOC(Core, FmInfo, CmInfo, TranInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
ENDIF


TranCntl%NowStep = 1
CALL UpdtResSrc(core, FmInfo, CmInfo, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)

!CALL PrintHeatPower(Core, CmInfo, FmInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
DO istep = 1, TranCntl%nStep
  ItrCntl%innerit = 0
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
  !Feedback for C5G7-TD DynamicBenchmark
  IF(TranCntl%lDynamicBen) THEN
    CALL UpdtFuelTemp(Core, CmInfo, TranInfo, TranCntl, PE, GroupInfo%ng, .TRUE.)
  END IF

  IF(TranCntl%lCondiMOC) THEN
    CALL CheckCondMOC(Core, FmInfo, CmInfo, TranInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
  ENDIF

  IF(nTracerCntl%lFeedback .AND. TranCntl%lSgFspUpdt) THEN
    IF(nTracerCntl%lrestrmt) CALL SubGrpFsp(Core, FmInfo%Fxr, THInfo, RayInfo, GroupInfo, nTracerCntl, PE)
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
  CALL TranReactivityUpdt(CmInfo, eigv, TranCntl, PE, TranInfo)
  !Save solution of previous time step
  CALL SaveTranSol(Core, FmInfo, CmInfo, TranInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)

  CALL ChkOutputWrite(lout, TranCntl)
  IF(lout) THEN
    CALL OutputEdit()
    !CALL PrintHeatPower(Core, CmInfo, FmInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
    TranCntl%NowWriteStep = TranCntl%NowWriteStep + 1
  ENDIF
  IF(PE%MASTER) THEN
    WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, 1e11.4, x, a, 1e10.3, x, a, 1e12.5)') '@1', TranCntl%T(istep),"s, P=", TranInfo%PowerLevel * 100, ', $= ', TranInfo%reactivity, ', G=', TranInfo%lifetime
    CALL message(io8, TRUE, TRUE, MESG)        
    WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, F7.1, x, a, F8.5, x, a, 1e10.3)') '@2', TranCntl%T(istep),"s, r(pcm)= ", TranInfo%delrho/10._8, ', k= ', TranInfo%TranEig/10._8, ', B= ', TranInfo%CoreBeta
    CALL message(io8, TRUE, TRUE, MESG)
    !WRITE(mesg, '(a, 1p, e12.4, x, a, p, e12.4)') '@3', TranCntl%T(istep),"s. Theta= ", TranCntl%Theta
    !CALL message(io8, TRUE, TRUE, MESG)
    IF(nTracerCntl%lAdjoint) THEN
      WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, 1e11.4, x, a, 1e10.3, x, a, 1e12.5)') '@3', TranCntl%T(istep),"s, P=", TranInfo%PowerLevel * 100, ', $= ', TranInfo%reactivity_Dynamic, ', G=', TranInfo%lifetime_dynamic_MOC
      CALL message(io8, TRUE, TRUE, MESG)    
      WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, F7.1, x, a, F8.5, x, a, 1e10.3, x, 1e10.3)') '@4', TranCntl%T(istep),"s, r(pcm)= ", &
        TranInfo%delrho_Dynamic/10._8, ', k= ', TranInfo%TranEig_Dynamic/10._8, ', B= ', TranInfo%CoreBeta_Dynamic_moc, TranInfo%CoreBeta_Dynamic
      CALL message(io8, TRUE, TRUE, MESG)
    ENDIF
    IF(nTracerCntl%lFeedback) THEN
      WRITE(mesg, '(a, 1p, e12.4, x, a, f7.1, x, a, f7.1, x, a)') '@5', TranCntl%T(istep),"s, Max Tf= ", ThInfo%TfMax/10, "C, Max TfAvg= ", ThInfo%Tfavg_max/10, "C"
      CALL message(io8, TRUE, TRUE, MESG)
      WRITE(mesg, '(a, 1p, e12.4, x, a, f7.1, x, a, f7.1, x, a)') '@6', TranCntl%T(istep),"s, Max Tcl= ", ThInfo%Tfcl_max/10, "C, Max TMOA= ", ThInfo%TModoutAvg/10, "C"
      CALL message(io8, TRUE, TRUE, MESG)
    ENDIF
    !601 format(11x,"Max Tf=",f7.1, "/" f7.1, "C,  Avg Outlet Temp",f7.2,"C", ", Dop. Change=",1p,e9.2,l2,i3)
    !WRITE(mesg, '(a, 1p, e12.4, x, a, p, F7.1, x, a, F8.5, x, a, 1e10.3)') '#1', TranCntl%T(istep),"s, r(pcm)= ", TranInfo%delrho_Dynamic/10._8, ', k= ', TranInfo%TranEig_Dynamic/10._8, ', B= ', TranInfo%CoreBeta_Dynamic
    
    IF(TranCntl%lSCM) THEN
      WRITE(mesg, '(a, 1p, e12.4, x, a, f10.3, a, f10.3)') '@7', TranCntl%T(istep),"s, AmpFrq= ", TranCntl%AmpFrq/10., "  Avg. Val. = ", TranCntl%AvgAmpFrq/10.
      CALL message(io8, TRUE, TRUE, MESG)
      WRITE(mesg, '(A, 1p, e12.4, x, a, f10.6)') '@8', TranCntl%T(istep),"s, Estimated Dynamic k-eff= ", TranCntl%EigD/10.      
      CALL message(io8, TRUE, TRUE, MESG)
    ENDIF
    
    WRITE(mesg, '(A, 1p, e12.4, x, a, i8)') '@9', TranCntl%T(istep), "s, MOC = ", ItrCntl%Mocit - ItrCntl%Mocit0
    CALL message(io8, TRUE, TRUE, MESG)
    WRITE(mesg, '(A, 1p, e12.4, x, a, i8, i8)') '@10', TranCntl%T(istep), "s, MGCMFD = ", ItrCntl%Cmfdit - ItrCntl%Cmfdit0, ItrCntl%innerit
    CALL message(io8, TRUE, TRUE, MESG)
    WRITE(mesg, '(A, 1p, e12.4, x, a, i8)') '@11', TranCntl%T(istep), "s, GCCMFD = ", ItrCntl%GcCmfdIt - ItrCntl%GcCmfdit0
    CALL message(io8, TRUE, TRUE, MESG)
    
    !CALL message(io8, TRUE, TRUE, MESG)
    WRITE(89, '(a)') mesg
    WRITE(mesg, '(A)') hbar2(1:77)
    CALL message(io8, FALSE, TRUE, MESG) 
  ENDIF
  TranCntl%NowStep = TranCntl%NowStep + 1
  TranInfo%PowerLeveld = TranInfo%PowerLevel
  TranCntl%lStepApprox = .FALSE.
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

DO iter = 1, TranCntl%nMaxOuter
  !TranCntl%lMOCUpdt = .TRUE.
  ItrCntl%SrcIt = ItrCntl%SrcIt + 1
  
  IF(nTracerCntl%lFeedback) THEN
    CALL UpdtPowerLevel(Core, FmInfo, CmInfo, TranInfo, GroupInfo%ng, nTracerCntl, PE)
    nTracerCntl%PowerLevel = TranInfo%PowerLevel; ThInfo%PowLv = nTracerCntl%PowerLevel
    CALL TransientTH(Core, CmInfo, FmInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
    CALL THFeedBack(Core, CmInfo, FmInfo, ThInfo, nTRACERCntl, PE)
    CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, TranInfo%Eigv0, GroupInfo, nTracerCntl, PE)
  ENDIF
  IF(lMOC) THEN
    CALL MOCDriverSwitch
  END IF
  IF(TranCntl%lExpTrsf) CALL UpdtExpTrsf(Core, FmInfo, CmInfo, TranInfo, GroupInfo, .FALSE., TranCntl, nTracerCntl, PE)
  CALL TranCmfd_Driver(Core, FmInfo, CmInfo, TranInfo, THInfo, GroupInfo, TranCntl, nTracerCntl, ItrCntl, PE)
  IF(iter .LT. 2) CYCLE
  IF(ItrCntl%lconv) EXIT
END DO
IF(TranCntl%lMocUpdt .AND. TranCntl%lCondiMOC) THEN
  CALL UpdtBaseXsCondiMOC(Core, FmInfo, CmInfo, TranInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
ENDIF
CALL CmFmInfoSync(Core, FmInfo, CmInfo, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)

  CONTAINS
  
  SUBROUTINE MOCDriverSwitch
#ifdef __PGI
  USE IEEE_ARITHMETIC
#endif

#ifdef __PGI
  IF (PE%lCUDA) THEN
    CALL CUDATranMOCSweep(Core, RayInfo, FMInfo, TranInfo, TranCntl); GOTO 2
  ENDIF
#endif
  CALL TranMOC_Driver(Core, RayInfo, FmInfo, CmInfo, TranInfo, THInfo, GroupInfo, TranCntl, nTracerCntl, ItrCntl, PE)

2 CONTINUE

!#ifdef __PGI
!  IF (ieee_is_nan(eigv)) CALL terminate("The Calculation Has Diverged!")
!#else
!  IF (isnan(eigv)) CALL terminate("The Calculation Has Diverged!")
!#endif
  END SUBROUTINE

  SUBROUTINE CmfdDriverSwitch
#ifdef __PGI
  USE IEEE_ARITHMETIC
#endif

  IF (lCMFD .AND. lCmfdGrp) THEN

#ifdef __PGI
    IF (PE%lCUDACMFD) THEN
      !CALL CUDACmfdAcc(Core, CMInfo, FMInfo, eigv); GOTO 3
    ENDIF
#endif
#ifdef __INTEL_MKL
    IF (PE%lMKL) THEN
      !IF (.NOT. mklCntl%lDavidson) THEN
      !  CALL MKL_CmfdPower(Core, CmInfo, FmInfo, eigv)
      !ELSE
      !  CALL MKL_CmfdDavidson(Core, CmInfo, FmInfo, eigv)
      !ENDIF
      GOTO 3
    ENDIF
#endif
    CALL TranCmfd_Driver(Core, FmInfo, CmInfo, TranInfo, THInfo, GroupInfo, TranCntl, nTracerCntl, ItrCntl, PE)

3   CONTINUE

!#ifdef __PGI
!    IF (ieee_is_nan(eigv)) CALL terminate("The Calculation Has Diverged!")
!#else
!    IF (isnan(eigv)) CALL terminate("The Calculation Has Diverged!")
!#endif

  ENDIF
  END SUBROUTINE

END SUBROUTINE

! SUBROUTINE Transient_Driver()
! USE PARAM
! USE GEOM,          ONLY : Core
! USE Core_mod,      ONLY : FmInfo
! USE PE_Mod,        ONLY : PE
! USE files,         ONLY : io8
! USE ioutil,        ONLY : message, terminate
! USE timer,         ONLY : nTracer_dclock,     TimeChk
! USE TRAN_MOD,      ONLY : TranInfo, TranCntl 
! USE XsPerturb_mod, ONLY : XsPerturbation
! USE CNTL,          ONLY : nTracerCntl
! IMPLICIT NONE

! #ifdef __PGI
! IF (PE%lCUDACMFD) THEN
!   CALL SSEIG()
!   CALL OutputEdit()
!   WRITE(mesg, '(A)') hbar2(1:77)
!   IF(PE%MASTER) CALL message(io8, FALSE, TRUE, MESG)
!   CALL CUDANNFSP_Driver()
!   GOTO 5
! ENDIF
! #endif
! #ifdef __INTEL_MKL
! #endif
! CALL terminate('NNFSP is only supported for CUDA version')
! 5 CONTINUE

! END SUBROUTINE

