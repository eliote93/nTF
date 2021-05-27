#include <defines.h>
#include <DefDBG.h>
SUBROUTINE sseig()

USE FXRVAR_MOD
USE EFT_MOD
USE XS_COMMON

USE PARAM,             ONLY : TRUE, FALSE, ZERO, mesg, lDepletion, lsseigv, leftsearch, ltransient, lbranch
USE GEOM,              ONLY : Core, ng
USE RAYS,              ONLY : RayInfo
USE Core_mod,          ONLY : CMInfo, FmInfo, THInfo, GroupInfo, GcGroupInfo, eigv, xkconv, xkcrit, eigv_adj
USE PE_MOD,            ONLY : PE
USE CNTL,              ONLY : nTracerCntl
USE Depl_Mod,          ONLY : DeplCntl
USE SubGrp_mod,        ONLY : SubGrpEffXsGen, SubGrpFsp
USE Th_Mod,            ONLY : SteadyStateTH, THFeedBack
USE CritSpec_mod,      ONLY : XsKinf
USE ItrCNTL_mod,       ONLY : ItrCntl, ConvItrCntl
USE FILES,             ONLY : io8
USE IOUTIL,            ONLY : message, ShowHBar1, terminate
USE CritSpec_mod,      ONLY : GetCriticalSpectrum
USE Boron_mod,         ONLY : UpdtBoronPPM, SetBoronCoolant, SetBoronCoolantCTF, MixBoronPPMCal
USE XeDyn_Mod,         ONLY : UpdtXeDyn, XeDynRelaxContrl
USE CNTLROD_mod,       ONLY : SetCrBankPosition, CntlRodSearch
USE MOC_COMMON,        ONLY : SetCoreMacXs, SetCoreMacXs_Cusping
USE MOC_MOD,           ONLY : GetNeighborMocFlux
use SubChCoupling_mod, ONLY : last_TH
USE BenchXs,           ONLY : SetBoronCoolant_NEACRP, DnpLambdaBen, ChidBen
USE TRAN_MOD,          ONLY : TranInfo
USE MacXsLib_Mod,      ONLY : PSMEffXSGen, GetfresoFXR
USE TranMacXsLib_Mod,  ONLY : InitChidLib, InitLambdaLib, initchidklib
USE HexEdit,           ONLY : HexOutputEdit

#ifdef inflow
USE Inflow_mod,  ONLY : InflowTrXS
#endif
#ifdef MPI_ENV
USE MPIComm_mod, ONLY : MPI_SYNC
#endif
#ifdef __INTEL_MKL
USE MKL_3D,      ONLY : mklCntl
USE MKL_CMFD,    ONLY : MKL_CmfdPower, MKL_CmfdDavidson
#ifdef __GAMMA_TRANSPORT
USE GAMMA_CMFD
#endif
#endif
#ifdef __PGI
USE AuxilDriver, ONLY : CUDASubGrpEffXSGen
#endif

IMPLICIT NONE

INTEGER :: iout
LOGICAL :: lTHconv, lsubgrp
LOGICAL :: Master, RTMaster, lCmfdGrp
LOGICAL :: lCMFD, lMOC
LOGICAL :: lSearchConv
! ----------------------------------------------------

Master   = PE%Master
RTMaster = PE%RTMaster
lCmfdGrp = PE%lCmfdGrp
lSubGrp  = FALSE
lCMFD    = nTracerCntl%lCMFD
lMOC     = TRUE

lSearchConv = FALSE
IF (.NOT. nTracerCntl%lSearch) lSearchConv = TRUE

ItrCntl%lconv      = FALSE
ItrCntl%OuterMax   = ConvItrCntl%OuterMax
ItrCntl%eigconv    = ConvItrCntl%eigconv
ItrCntl%resconv    = ConvItrCntl%resconv
ItrCntl%psiconv    = ConvItrCntl%psiconv
ItrCntl%decuspconv = ConvItrCntl%decuspconv

IF (nTracerCntl%lSSPH .AND.(nTracerCntl%lProblem.EQ.lDepletion .OR. nTracerCntl%lXeDyn .OR. (nTracerCntl%BoronPPM.GT.0._8) .OR. nTracerCntl%lCrInfo)) nTracerCntl%lSSPHreg = TRUE

#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif

! Initialize Iteration Variables
CALL InitIterVar(Core, FMInfo, CmInfo, GroupInfo, TRUE, ItrCntl, nTracerCntl, PE)

! Set CR Position
IF (nTracerCntl%lCrInfo .AND. Core%CrInfo%lCrChg) CALL SetCrBankPosition(Core, FmInfo, CmInfo, GroupInfo, nTracerCntl, PE)

! Set Boron Concentration
IF (nTracerCntl%lXsLib) THEN
  IF (nTracerCntl%BoronPPM .GT. ZERO) THEN
    WRITE (mesg,'(A, F12.2, A)') ' =======  Applied Boron PPM : ',nTracerCntl%boronppm, ' PPM'
    IF (master) CALL message(io8, TRUE, TRUE, mesg)
    
    CALL SetBoronCoolant(Core, FmInfo%Fxr, nTracerCntl%boronppm, PE%myzb, PE%myze)
  ELSE
    CALL MixBoronPPMCal(Core,FmInfo%Fxr, nTracerCntl%boronppm)
  END IF
END IF

WRITE (mesg, '(A, L3)') 'MOC Under-relaxtion Option : ', nTracerCntl%lMOCUR
IF (master) CALL message(io8, TRUE, TRUE, mesg)
WRITE (mesg, '(A, L3)') 'Axial Reflector FDM Option : ', nTracerCntl%lAxRefFDM
IF (master) CALL message(io8, TRUE, TRUE, mesg)

! Xe Dyn
IF (nTracerCntl%lXeDyn) CALL XeDynRelaxContrl(FALSE)

! Align XS
IF (nTracerCntl%lXsAlign) CALL XsLinIntPolTemp(nTracerCntl%ScatOd)
! ----------------------------------------------------
IF (nTracerCntl%lFeedBack) THEN
  lThConv = FALSE
  
  IF (nTracerCntl%lSSPH) nTracerCntl%lSSPHreg = TRUE
  
  IF (nTracerCntl%lMacro) CALL SetCoreMacXs(Core, FmInfo)
  
  CALL CmfdDriverSwitch()
  
#ifdef __INTEL_MKL
#ifdef __GAMMA_TRANSPORT
  IF (nTracerCntl%lGamma) THEN
    WRITE (mesg, '(A12, I1, A31)') 'Performing P', nTracerCntl%GammaScatOd, ' Gamma Transport Calculation...'
    IF (Master) CALL message(io8, TRUE, TRUE, mesg)
    
    ItrCntl%lGammaConv = FALSE
    
    DO WHILE (.NOT. ItrCntl%lGammaConv)
      IF (core%nxy .NE. 1) CALL GammaCMFDAcc(Core, CmInfo, FmInfo)
      
      CALL GammaMOCSweep(Core, RayInfo, FmInfo, PE, nTracerCntl, ItrCntl)
    END DO
  END IF
#endif
#endif

  CALL SteadyStateTH(Core, CmInfo, FmInfo, ThInfo, Eigv, ng, GroupInfo, nTracerCntl, PE)
  CALL THFeedBack(Core, CmInfo, FmInfo, ThInfo, nTracerCntl, PE)
  
  IF (nTracerCntl%lXsAlign) THEN
    CALL UpdateBlockFxr(FmInfo%Fxr, GroupInfo, THFeed)
    CALL XsLinIntPolTemp(nTracerCntl%ScatOd)
    CALL UpdateResPinFuelTemp
  END IF

  IF (nTRACERCntl%lDeplVarInit) CALL InitIterVar(Core, FMInfo, CmInfo, GroupInfo, TRUE, ItrCntl, nTracerCntl, PE)

  ! Update Xenon Dynamics
  IF (nTracerCntl%lXeDyn) THEN
    CALL UpdtXeDyn(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, DeplCntl, PE)
    
    IF (nTracerCntl%lXsAlign) THEN
      CALL UpdateBlockFxr(FmInfo%Fxr, GroupInfo, XeDynFeed)
      CALL XsLinIntPolTemp(nTracerCntl%ScatOd)
    END IF
  END IF
  
  ! Generate Self-Shielded XS
  IF (nTracerCntl%lrestrmt) CALL SubgrpDriverSwitch()
  
#ifdef MPI_ENV
  CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif
  
  ! Effective XS Calculation
  IF (.NOT. nTRACERCntl%lPSM) THEN
#ifdef __PGI
    IF (nTracerCntl%lXsAlign) THEN
      CALL CUDASubGrpEffXSGen(core,FmInfo%fxr,nTracerCntl,PE)
    ELSE
      CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
    END IF
#else
    CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, Eigv, GroupInfo, nTracerCntl, PE)
#endif
  ELSE
    CALL GetfresoFXR(Core, FmInfo%Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
  END IF
  
  IF (nTracerCntl%lMacro) CALL SetCoreMacXs(Core, FmInfo)

  CALL CmfdDriverSwitch()
! ----------------------------------------------------
ELSE
#ifdef MPI_ENV
  CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif

  IF (nTRACERCntl%lDeplVarInit) CALL InitIterVar(Core, FMInfo, CmInfo, GroupInfo, TRUE, ItrCntl, nTracerCntl, PE)

  ! Update Xenon Dynamics
  IF (nTracerCntl%lXeDyn) THEN
    CALL UpdtXeDyn(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, DeplCntl, PE)

    IF (nTracerCntl%lXsAlign) CALL UpdateBlockFxr(FmInfo%Fxr, GroupInfo, XeDynFeed)
  END IF
  
  IF (nTracerCntl%lXsAlign) THEN
    CALL UpdateBlockFxr(FmInfo%Fxr, GroupInfo, THFeed)
    CALL XsLinIntPolTemp(nTracerCntl%ScatOd)
    CALL UpdateResPinFuelTemp
  END IF

  ! Generate Self-Shielded XS
  IF (nTracerCntl%lrestrmt) THEN
    CALL SubgrpDriverSwitch()
    
    lSubGrp = TRUE
  END IF

  ! Effective XS Calculation
  IF (.NOT. nTRACERCntl%lPSM) THEN
#ifdef __PGI
    IF (nTracerCntl%lXsAlign) THEN
      CALL CUDASubGrpEffXSGen(core,FmInfo%fxr,nTracerCntl,PE)
    ELSE
      CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
    END IF
#else
    CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, Eigv, GroupInfo, nTracerCntl, PE)
#endif
  ELSE
    CALL GetfresoFXR(Core, FmInfo%Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
  END IF
  
  IF (nTracerCntl%lMacro) CALL SetCoreMacXs(Core, FmInfo)

  CALL CmfdDriverSwitch()
END IF
! ----------------------------------------------------
IF (nTracerCntl%lBoronSearch) THEN
  CALL UpdtBoronPPM(nTracerCntl%target_eigv, eigv, nTracerCntl%BoronPPM, TRUE, MASTER)
  
  IF (nTracerCntl%lXslib) THEN
    CALL SetBoronCoolant(Core, FmInfo%Fxr, nTracerCntl%boronppm, PE%myzb, PE%myze)
    
    IF (nTracerCntl%lMacro) CALL SetCoreMacXs(Core, FmInfo)
    
    IF (nTracerCntl%lXsAlign) THEN
      CALL UpdateBlockFxr(FmInfo%Fxr, GroupInfo, PPMFeed)
      CALL XsLinIntPolTemp(nTracerCntl%ScatOd)
    END IF
  ELSE IF(nTracerCntl%libtyp .EQ. 11) THEN
    CALL SetBoronCoolant_NEACRP(nTracerCntl%boronppm)
  END IF
END IF

IF (nTracerCntl%lXeDyn) CALL XeDynRelaxContrl(TRUE)

#ifdef HISTORY_OUT
CALL PowerHistoryOutputEdit(TRUE)
#endif
! ----------------------------------------------------
DO iout = 1, ItrCntl%OuterMax
#ifdef MPI_ENV
  CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif
  
  ItrCntl%SrcIt = ItrCntl%SrcIt + 1
  
  ! MoC
  IF (nTracerCntl%lMOCUR .AND. RTMASTER) CALL MOCUnderRelaxationFactor(Core, FmInfo, CmInfo, GroupInfo, nTracerCntl, PE)
  
  IF (lMOC) CALL MOCDriverSwitch()
  
  ! Xe Dyn
  IF (nTracerCntl%lXeDyn) THEN
    CALL UpdtXeDyn(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, DeplCntl, PE)
    
    IF (nTracerCntl%lXsAlign) THEN
      CALL UpdateBlockFxr(FmInfo%Fxr, GroupInfo, XeDynFeed)
      CALL XsLinIntPolTemp(nTracerCntl%ScatOd)
    END IF
    
    IF (nTRACERCntl%lPSM) THEN
      CALL GetfresoFXR(Core, FmInfo%Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
    ELSE
#ifdef __PGI
      IF (nTracerCntl%lXsAlign) THEN
        CALL CUDASubGrpEffXSGen(core,FmInfo%fxr,nTracerCntl,PE)
      ELSE
        CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, Eigv, GroupInfo, nTracerCntl, PE)
      END IF
#else
      CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
#endif
    END IF
    
    IF (nTracerCntl%lMacro .AND. .NOT. nTracerCntl%lCusping_MPI) CALL SetCoreMacXs(Core, FmInfo)
  END IF
  
  ! CMFD
  CALL CmfdDriverSwitch()
  
#ifdef HISTORY_OUT
  CALL PowerHistoryOutputEdit(FALSE)
#endif
  
  ! T/H Update
  IF (iout .EQ. ItrCntl%OuterMax) last_TH = TRUE
  
  IF (nTracerCntl%lFeedBack .AND. .NOT.lThConv) THEN
    IF (nTracerCntl%lSSPH) nTracerCntl%lSSPHreg = TRUE
    
    CALL SteadyStateTH(Core, CmInfo, FmInfo, ThInfo, Eigv, ng, GroupInfo, nTracerCntl, PE)
    CALL THFeedBack(Core, CmInfo, FmInfo, ThInfo, nTracerCntl, PE)
    
    IF (nTracerCntl%lXsAlign) THEN
      CALL UpdateBlockFxr(FmInfo%Fxr, GroupInfo, THFeed)
      CALL XsLinIntPolTemp(nTracerCntl%ScatOd)
      CALL UpdateResPinFuelTemp
    END IF
    
#ifdef MPI_ENV
    CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif
    
    IF (nTracerCntl%lrestrmt) THEN
      IF (ThInfo%TdopChg .LT. 1.D-03 .AND. .NOT. lSubGrp) THEN
        CALL SubgrpDriverSwitch()
        lSubGrp = TRUE
      END IF
      
      IF (nTRACERCntl%lPSM) THEN
        CALL GetfresoFXR(Core, FmInfo%Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
      ELSE
#ifdef __PGI
        IF (nTracerCntl%lXsAlign) THEN
          CALL CUDASubGrpEffXSGen(core,FmInfo%fxr,nTracerCntl,PE)
        ELSE
          CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, Eigv, GroupInfo, nTracerCntl, PE)
        END IF
#else
        CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
#endif
      END IF
    END IF

    IF (nTracerCntl%lMacro .AND. .NOT. nTracerCntl%lCusping_MPI) CALL SetCoreMacXs(Core, FmInfo)
  END IF
  
  ! CHK : Cnv.
  IF (nTracerCntl%lSearch) THEN
    IF (abs(eigv - nTracerCntl%target_eigv) / nTracerCntl%target_eigv .LT. 1E-5) lSearchConv = TRUE
    IF (nTracerCntl%lBoronSearch .AND. nTracerCntl%BoronPPM .EQ. ZERO) lSearchConv = TRUE
  END IF
  
  IF (nTracerCntl%lFeedBack .AND. ThInfo%TdopChg .LT. 5E-5) lThConv = TRUE

  IF (MASTER) THEN
    WRITE(mesg, '(2x, A,L3, 3x, A, L3, 3x, A, L3)') 'Convergence - Flux :', ItrCntl%lconv, 'TH :',lThConv, 'B-Search:', lSearchConv
    CALL message(io8, FALSE, TRUE, mesg)
  END IF
  
  IF (ItrCntl%lconv .AND. lSearchConv) THEN
#ifdef __GAMMA_TRANSPORT
    IF (nTracerCntl%lGamma) THEN
      WRITE(mesg, '(A12, I1, A31)') 'Performing P', nTracerCntl%GammaScatOd, ' Gamma Transport Calculation...'
      IF (Master) CALL message(io8, TRUE, TRUE, mesg)
      
      ItrCntl%lGammaConv = FALSE
      
      DO WHILE (.NOT. ItrCntl%lGammaConv)
#ifdef __INTEL_MKL
        IF (Core%nxy .NE. 1) CALL GammaCMFDAcc(Core, CmInfo, FmInfo)
#endif
        CALL GammaMOCSweep(Core, RayInfo, FmInfo, PE, nTracerCntl, ItrCntl)
      END DO
    END IF
#endif
    
    IF (.NOT. nTracerCntl%lFeedBack) EXIT
    IF (nTracerCntl%lFeedBack .AND. lThConv) EXIT
  END IF
  
  ! Search
  IF (nTracerCntl%lSearch) THEN
    IF (nTracerCntl%lBoronSearch) THEN
      CALL UpdtBoronPPM(nTracerCntl%target_eigv, eigv, nTracerCntl%BoronPPM, FALSE, MASTER)
      
      IF (nTracerCntl%lXslib) THEN
        CALL SetBoronCoolant(Core, FmInfo%Fxr, nTracerCntl%boronppm, PE%myzb, PE%myze)
        
        IF (nTracerCntl%lXsAlign) THEN
          CALL UpdateBlockFxr(FmInfo%Fxr, GroupInfo, PPMFeed)
          CALL XsLinIntPolTemp(nTracerCntl%ScatOd)
        END IF
      ELSE IF(nTracerCntl%libtyp .EQ. 11) THEN
        CALL SetBoronCoolant_NEACRP(nTracerCntl%boronppm)
      END IF
      
      IF (nTracerCntl%lMacro .AND. .NOT.nTracerCntl%lCusping_MPI) CALL SetCoreMacXs(Core, FmInfo)
    END IF
    
    IF (nTracerCntl%lCrSearch .AND. iout.GT.1) THEN
      CALL CntlRodSearch(Core, FmInfo, CmInfo, eigv, GroupInfo, GcGroupInfo, nTracerCntl, PE)
      CALL SetBoronCoolant(Core, FmInfo%Fxr, nTracerCntl%boronppm, PE%myzb, PE%myze)
    END IF
  END IF
END DO
! ----------------------------------------------------
IF (nTracerCntl%lXeDyn) CALL XeDynRelaxContrl(FALSE)

xkconv = eigv

IF (nTracerCntl%lCritSpec .AND. RTMASTER) CALL GetCriticalSpectrum(Core, FmInfo, GroupInfo, nTracerCntl, PE)

IF (nTracerCntl%lAdjoint) THEN
  IF (MASTER) THEN
    IF (Master) CALL ShowHbar1(io8)
    
    WRITE(mesg, *) 'Start to Calculate Adjoint Flux at Initial State...'
    CALL message(io8, TRUE, TRUE, mesg)
  END IF
  
  ItrCntl%lConv = FALSE
  
  IF (nTracerCntl%lCMFDAdj) CALL CMFD_adj_Switch()
END IF

IF (nTracerCntl%lproblem.EQ.lsseigv .AND. RTMASTER) THEN
  xkcrit = xkconv
  IF (nTracerCntl%lXsLib) xkcrit = XsKinf(Core, FmInfo, GroupInfo, nTracerCntl%lXsLib, nTracerCntl%lCritSpec, PE)
  
  WRITE(mesg, '(A,2(F9.3,A),2(A,F9.5),A,F8.2,A)') 'Burnup = ',DeplCntl%NowBurnUp(2), ' MWD/kgHM', DeplCntl%NowBurnUp(1), ' EPFD','  keff = ', xkconv, &
                                                  '  kinf(BK) = ', xkcrit, '  ppm = ', nTracerCntl%boronppm, ' ppm'
  IF (MASTER) CALL message(io8, FALSE, TRUE, MESG)
END IF

! Print Out Eigenvalue
IF (nTracerCntl%lproblem.EQ.lsseigv .OR. nTracerCntl%lproblem.EQ.lEFTsearch ) THEN
  IF (Master) CALL ShowHbar1(io8)
  
  WRITE (mesg, '(A27, 5X, 1F8.5)')  'k-eff =', Eigv
  IF (Master) CALL message(io8, FALSE, TRUE, MESG)

  IF (nTracerCntl%lAdjoint) THEN
    WRITE (mesg, '(A27, 5X, 1F8.5)')  'Adjoint k-eff =', Eigv_adj
    IF (Master) CALL message(io8, FALSE, TRUE, MESG)
  END IF
  
  IF (PE%lCmfdGrp) THEN
    IF (nTracerCntl%lHex) THEN
      CALL HexOutputEdit()
    ELSE
      CALL OutputEdit()
    END IF
  END IF
ELSE IF (nTracerCntl%lproblem .EQ. lTransient) THEN
  IF (Master) CALL ShowHbar1(io8)
  
  WRITE (mesg, '(A27, 5X, 1F8.5)')  'k-eff =', Eigv
  IF (Master) CALL message(io8, FALSE, TRUE, MESG)

  IF (nTracerCntl%lAdjoint) THEN
    WRITE(mesg, '(A27, 5X, 1F8.5)')  'Adjoint k-eff =', Eigv_adj
    IF(Master) CALL message(io8, FALSE, TRUE, MESG)
  END IF
END IF

IF (nTracerCntl%lproblem .EQ. lBranch) THEN
  IF (Master) CALL ShowHbar1(io8)
  
  SELECT CASE(CurrentVarId)
  CASE (branch_tmod)
    WRITE (mesg, '(A27, 5X, 1F8.5,a,i2,a,i2)')  'branch k-eff =', Eigv, '  Tmod  CaseId : ',CurrentCaseId, ' /', ntmod
  CASE (branch_tfuel)
    WRITE (mesg, '(A27, 5X, 1F8.5,a,i2,a,i2)')  'branch k-eff =', Eigv, '  Tfuel CaseId : ',CurrentCaseId, ' /', ntfuel
  CASE (branch_boron)
    WRITE (mesg, '(A27, 5X, 1F8.5,a,i2,a,i2)')  'branch k-eff =', Eigv, '  Boron CaseId : ',CurrentCaseId, ' /', nboron
  CASE (branch_rho)
    WRITE (mesg, '(A27, 5X, 1F8.5,a,i2,a,i2)')  'branch k-eff =', Eigv, '  Rho   CaseId : ',CurrentCaseId, ' /', nrho
  END SELECT
  IF (Master) CALL message(io8, FALSE, TRUE, MESG)
  
  IF (PE%lCmfdGrp) THEN
    IF (nTracerCntl%lHex) THEN
      CALL HexOutputEdit()
    ELSE
      CALL OutputEdit()
    END IF
  END IF
END IF

IF (nTracerCntl%lGC .AND. (lEFT_GCGEN.OR..NOT.nTracerCntl%lEFTSearch)) THEN
    WRITE(mesg, '(A)') 'Generating Group Constants ...'
    CALL message(io8, TRUE, TRUE, mesg)
    
    IF (nTracerCntl%lTranON) THEN
      CALL AllocTransient()
      CALL KinParamGen(Core, FmInfo, TranInfo, ThInfo, GroupInfo, TRUE, nTracerCntl, PE)
      
      IF (nTRACERCntl%lXsLib) THEN
        CALL InitLambdaLib(TranInfo%Lambda, nTracerCntl%refdcy_del, nTracerCntl%llibdcy_del, nTracerCntl%lfitbeta)
        
        IF (nTracerCntl%lchidkgen) THEN
          CALL InitChidkLib(Core, FmInfo, GroupInfo, PE, nTracerCntl, TranInfo%chid, TranInfo%Chidk, ng, TranInfo%nprec)
        ELSE
          CALL InitChidLib(Core, FmInfo, GroupInfo, PE, nTracerCntl, TranInfo%Chid, ng)
        END IF
      ELSE
        CALL ChidBen(TranInfo%Chid)
        CALL DnpLambdaBen(TranInfo%Lambda)
      END IF
    END IF
    
    IF (nTracerCntl%gc_spec .GE. 3) THEN
        IF(RTmaster) CALL GroupConstGenMac(Core, FmInfo, CmInfo, GroupInfo, nTracerCntl, PE, ng)
    ELSE
        IF(RTmaster) CALL GCGen(Core, FmInfo, THInfo, CmInfo, GroupInfo, nTracerCntl, PE, ng)
    END IF
    
    WRITE(mesg, '(A)') 'Generating Group Constants ... END'
    CALL message(io8, TRUE, TRUE, mesg)
END IF
! ------------------------------------------------------------------------------------------------------------
CONTAINS
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SubgrpDriverSwitch

IF (nTRACERCntl%lPSM) THEN
  CALL PSMEffXSGen(Core, FmInfo%Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
ELSE
#ifdef __PGI
  IF (PE%lCUDA) THEN
    IF (nTracerCntl%lXsAlign) THEN
      CALL CUDASubGrpFsp_cuMLG(Core, FmInfo, THInfo, RayInfo, GroupInfo)
    ELSE
      CALL CUDASubGrpFsp(Core, FmInfo, THInfo, RayInfo, GroupInfo)
    END IF
    
    GOTO 1
  END IF
#endif
  CALL SubGrpFsp(Core, FmInfo%Fxr, THInfo, RayInfo, GroupInfo, nTracerCntl, PE)
END IF

1 CONTINUE

END SUBROUTINE SubgrpDriverSwitch
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE MOCDriverSwitch

#ifdef __PGI
USE IEEE_ARITHMETIC
#endif

#ifdef __PGI
IF (PE%lCUDA) THEN
  IF (nTracerCntl%lXsAlign) THEN
    CALL CUDAMOCSweep_cuMac(Core, RayInfo, FMInfo, eigv)
  ELSE
    CALL CUDAMOCSweep(Core, RayInfo, FMInfo, eigv)
  END IF
  
  GOTO 2
END IF
#endif

CALL MOCSweep(Core, RayInfo, FMInfo, eigv, ng, PE, nTracerCntl, ItrCntl)

2 CONTINUE

#ifdef __PGI
IF (ieee_is_nan(eigv)) CALL terminate("The Calculation Has Diverged!")
#else
IF (isnan(eigv)) CALL terminate("The Calculation Has Diverged!")
#endif

END SUBROUTINE MOCDriverSwitch
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE CmfdDriverSwitch

#ifdef __PGI
USE IEEE_ARITHMETIC
#endif

IF (.NOT.lCMFD .OR. .NOT.lCmfdGrp) RETURN

IF(nTracerCntl%lMacro .AND. nTracerCntl%lCusping_MPI .AND. (ItrCntl%srcit0 .NE. ItrCntl%srcit)) THEN
  CALL GetNeighborMocFlux(FmInfo%phis, FmInfo%neighphis, Core%nCoreFsr, PE%myzb, PE%myze, 1, GroupInfo%ng, Core%nz, Core%AxBC)
  CALL SetCoreMacXS_Cusping(Core, FmInfo)
END IF

#ifdef __PGI
IF (PE%lCUDACMFD) THEN
  IF (nTracerCntl%lXsAlign) THEN
    CALL CUDACmfdAcc_cuMAC(Core, CMInfo, FMInfo, eigv)
  ELSE
    CALL CUDACmfdAcc(core,CmInfo,FmInfo,eigv)
  END IF
  
  GOTO 3
END IF
#endif

#ifdef __INTEL_MKL
IF (PE%lMKL .AND. mklCntl%lEnter) THEN
  IF (.NOT. mklCntl%lDavidson) THEN
    CALL MKL_CmfdPower(Core, CmInfo, FmInfo, eigv)
  ELSE
    CALL MKL_CmfdDavidson(Core, CmInfo, FmInfo, eigv)
  END IF
  
  GOTO 3
END IF
#endif

IF (nTracerCntl%lHex) CALL terminate("HEX MUST USE MKL CMFD!")
CALL CmfdAcc(Core, CMInfo, FMInfo, THInfo, eigv, ng, TRUE, FALSE, PE, nTracerCntl, ItrCntl)

3 CONTINUE

#ifdef __PGI
IF (ieee_is_nan(eigv)) CALL terminate("The Calculation Has Diverged!")
#else
IF (isnan(eigv)) CALL terminate("The Calculation Has Diverged!")
#endif

END SUBROUTINE CmfdDriverSwitch
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE CMFD_adj_Switch

#ifdef __PGI
USE IEEE_ARITHMETIC
#endif

IF (.NOT.lCMFD .OR. .NOT.lCmfdGrp) RETURN

#ifdef __PGI
IF (PE%lCUDACMFD) THEN
  CALL CUDACmfd_Adj(Core, CMInfo, eigv_adj)
  GOTO 4
END IF
#endif

#ifdef __INTEL_MKL
IF (PE%lMKL) THEN
  CALL terminate("CMFD Adjoint Flux Calculation with MKL is under construction")
  GOTO 4
END IF
#endif

CALL terminate("CMFD Adjoint Flux Calculation is under construction")

4 CONTINUE

#ifdef __PGI
IF (ieee_is_nan(eigv_adj)) CALL terminate("The Calculation Has Diverged!")
#else
IF (isnan(eigv_adj)) CALL terminate("The Calculation Has Diverged!")
#endif

END SUBROUTINE CMFD_adj_Switch
! ------------------------------------------------------------------------------------------------------------

END SUBROUTINE sseig