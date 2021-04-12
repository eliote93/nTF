#include <defines.h>
#include <DefDBG.h>

SUBROUTINE SSEIG()
USE PARAM
USE GEOM,             ONLY : Core, ng
USE RAYS,             ONLY : RayInfo
USE Core_mod,         ONLY : CMInfo, FmInfo, THInfo, GroupInfo, GcGroupInfo,   &
                             eigv, peigv, xkconv, xkcrit, eigv_adj
USE PE_MOD,           ONLY : PE
USE CNTL,             ONLY : nTracerCntl
USE Depl_Mod,         ONLY : DeplCntl
USE SubGrp_mod,       ONLY : SubGrpEffXsGen, SubGrpFsp
USE Th_Mod,           ONLY : SteadyStateTH, THFeedBack
USE CritSpec_mod,     ONLY : XsKinf
USE ItrCNTL_mod,      ONLY : ItrCntl, ConvItrCntl
USE FILES,            ONLY : io8
USE IOUTIL,           ONLY : message, ShowHBar1,  terminate
USE CritSpec_mod,     ONLY : GetCriticalSpectrum
USE Boron_mod,        ONLY : UpdtBoronPPM, SetBoronCoolant, SetBoronCoolantCTF, &
                              MixBoronPPMCal
USE XeDyn_Mod,        ONLY : UpdtXeDyn, XeDynRelaxContrl
USE CNTLROD_mod,      ONLY : SetCrBankPosition,      CntlRodSearch
USE FXRVAR_MOD
#ifdef inflow
USE Inflow_mod,       ONLY : InflowTrXS
#endif
#ifdef MPI_ENV
USE MPIComm_mod,      ONLY : MPI_SYNC
#endif
#ifdef __INTEL_MKL
USE MKL_3D,           ONLY : mklCntl
USE MKL_CMFD,         ONLY : MKL_CmfdPower,          MKL_CmfdDavidson
#ifdef __GAMMA_TRANSPORT
USE GAMMA_CMFD
#endif
#endif
USE MOC_COMMON,       ONLY : SetCoreMacXs,           SetCoreMacXs_Cusping
USE MOC_MOD,          ONLY : GetNeighborMocFlux
use SubChCoupling_mod,    only: last_TH
USE EFT_MOD
USE XS_COMMON
#ifdef __PGI
USE AuxilDriver,      ONLY : CUDASubGrpEffXSGen
#endif
USE BenchXs,          ONLY : SetBoronCoolant_NEACRP
USE TRAN_MOD,         ONLY : TranCntl, TranInfo
USE XsPerturb_mod,    ONLY : XsCntlRodMix
USE PointXSRT_MOD,    ONLY : CalcPSM_ISOPIN, CalcShadowingCorrection
USE MacXsLib_Mod,     ONLY : PSMEffXSGen, GetfresoFXR

USE BenchXs,          ONLY : DnpLambdaBen,      ChidBen
USE TranMacXsLib_Mod, ONLY : InitChidLib,       InitLambdaLib,    initchidklib
IMPLICIT NONE
INTEGER :: i, j, n
LOGICAL :: lTHconv, lsubgrp
LOGICAL :: Master, Slave, RTMaster, RTSlave, CMFDMaster, CMFDslave, lCmfdGrp
LOGICAL :: lCMFD, lMOC
LOGICAL :: lSearchConv

Master = PE%Master; Slave = PE%Slave
RTMaster = PE%RTMaster; RTSlave = PE%RTSlave
CMFDMaster = PE%CMFDMaster; CMFDSlave = PE%CMFDSlave
lCmfdGrp = PE%lCmfdGrp
lSubGrp = .FALSE.
lCMFD = nTracerCntl%lCMFD
lMOC = .TRUE.
lSearchConv = .FALSE.
IF (.NOT. nTracerCntl%lSearch) lSearchConv = .TRUE.
ItrCntl%lconv = .FALSE.
ItrCntl%OuterMax = ConvItrCntl%OuterMax
ItrCntl%eigconv = ConvItrCntl%eigconv
ItrCntl%resconv = ConvItrCntl%resconv
ItrCntl%psiconv = ConvItrCntl%psiconv
ItrCntl%decuspconv = ConvItrCntl%decuspconv

IF (nTracerCntl%lProblem .EQ. lDepletion .or. nTracerCntl%lXeDyn .or. (nTracerCntl%BoronPPM.GT.0._8) .OR. nTracerCntl%lCrInfo) THEN
   IF (nTracerCntl%lSSPH) nTracerCntl%lSSPHreg=.TRUE.
ENDIF

#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif

! Initialize Iteration Variables
CALL InitIterVar(Core, FMInfo, CmInfo, GroupInfo, .TRUE., ItrCntl, nTracerCntl, PE)

IF (nTracerCntl%lCrInfo) THEN
  IF (Core%CrInfo%lCrChg) THEN
    CALL SetCrBankPosition(Core, FmInfo, CmInfo, GroupInfo, nTracerCntl, PE)
  ENDIF
ENDIF

! Set Boron Concentration
!IF (nTracerCntl%lXsLib .AND. nTracerCntl%BoronPPM .GT. 0._8) THEN
IF (nTracerCntl%lXsLib) THEN
  IF (nTracerCntl%BoronPPM .GT. 0._8) THEN
  WRITE(mesg,'(a,f12.2,a)') ' =======  Applied Boron PPM : ',nTracerCntl%boronppm, ' PPM'
  IF(master) CALL message(io8, TRUE, TRUE, mesg)
  CALL SetBoronCoolant(Core, FmInfo%Fxr, nTracerCntl%boronppm, PE%myzb, PE%myze)
  ELSE
    CALL MixBoronPPMCal(Core,FmInfo%Fxr, nTracerCntl%boronppm)
  END IF
ENDIF
WRITE(mesg, '(A, L3)') 'MOC Under-relaxtion Option : ', nTracerCntl%lMOCUR
IF(master) CALL message(io8,TRUE, TRUE, mesg)
WRITE(mesg, '(A, L3)') 'Axial Reflector FDM Option : ', nTracerCntl%lAxRefFDM
IF(master) CALL message(io8,TRUE, TRUE, mesg)

IF (nTracerCntl%lXeDyn) THEN
  CALL XeDynRelaxContrl(.FALSE.)
END IF

IF (nTracerCntl%lXsAlign) THEN
  CALL XsLinIntPolTemp(nTracerCntl%ScatOd)
END IF

IF (nTracerCntl%lFeedBack) THEN

  lThConv = FALSE
  IF (nTracerCntl%lSSPH) nTracerCntl%lSSPHreg=.TRUE.

 IF (nTracerCntl%lMacro) CALL SetCoreMacXs(Core, FmInfo)
  CALL CmfdDriverSwitch()

#if defined (__INTEL_MKL) && defined (__GAMMA_TRANSPORT)
  IF (nTracerCntl%lGamma) THEN
    WRITE(mesg, '(a12, i1, a31)') 'Performing P', nTracerCntl%GammaScatOd, ' Gamma Transport Calculation...'
    IF (Master) CALL message(io8, TRUE, TRUE, mesg)
    ItrCntl%lGammaConv = FALSE
    DO WHILE (.NOT. ItrCntl%lGammaConv)
      IF (core%nxy .NE. 1) CALL GammaCMFDAcc(Core, CmInfo, FmInfo)
      CALL GammaMOCSweep(Core, RayInfo, FmInfo, PE, nTracerCntl, ItrCntl)
    ENDDO
  ENDIF
#endif

  CALL SteadyStateTH(Core, CmInfo, FmInfo, ThInfo, Eigv, ng, GroupInfo, nTracerCntl, PE)
  CALL THFeedBack(Core, CmInfo, FmInfo, ThInfo, nTracerCntl, PE)
  IF (nTracerCntl%lXsAlign) THEN
    CALL UpdateBlockFxr(FmInfo%Fxr, GroupInfo, THFeed)
    CALL XsLinIntPolTemp(nTracerCntl%ScatOd)
    CALL UpdateResPinFuelTemp
  END IF

  IF (nTRACERCntl%lDeplVarInit) CALL InitIterVar(Core, FMInfo, CmInfo, GroupInfo, .TRUE., ItrCntl, nTracerCntl, PE)

  ! Update Xenon Dynamics
  IF (nTracerCntl%lXeDyn) THEN
    CALL UpdtXeDyn(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, DeplCntl, PE)
    IF (nTracerCntl%lXsAlign) THEN
      CALL UpdateBlockFxr(FmInfo%Fxr, GroupInfo, XeDynFeed)
      CALL XsLinIntPolTemp(nTracerCntl%ScatOd)
    END IF
  ENDIF

  ! Generate Self-Shielded XS
  IF (nTracerCntl%lrestrmt) CALL SubgrpDriverSwitch()

#ifdef MPI_ENV
  CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif

! Effective XS Calculation
IF (.NOT. nTRACERCntl%lPSM) THEN
#ifdef __PGI
    IF (nTracerCntl%lXsAlign) THEN
!      CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
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

ELSE
#ifdef MPI_ENV
  CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif

  IF (nTRACERCntl%lDeplVarInit) CALL InitIterVar(Core, FMInfo, CmInfo, GroupInfo, .TRUE., ItrCntl, nTracerCntl, PE)

  ! Update Xenon Dynamics
  IF (nTracerCntl%lXeDyn) THEN
    CALL UpdtXeDyn(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, DeplCntl, PE)
!    print*, 'Updt Xe Done!'
    IF (nTracerCntl%lXsAlign) THEN
      CALL UpdateBlockFxr(FmInfo%Fxr, GroupInfo, XeDynFeed)
    END IF
!    print*, 'Post Updt Xe Done!'
  ENDIF

  IF (nTracerCntl%lXsAlign) THEN
    CALL UpdateBlockFxr(FmInfo%Fxr, GroupInfo, THFeed)
    CALL XsLinIntPolTemp(nTracerCntl%ScatOd)
    CALL UpdateResPinFuelTemp
  ENDIF

  ! Generate Self-Shielded XS
  IF (nTracerCntl%lrestrmt) THEN
    CALL SubgrpDriverSwitch()
    lSubGrp = .TRUE.
  ENDIF

! Effective XS Calculation
IF (.NOT. nTRACERCntl%lPSM) THEN
#ifdef __PGI
  IF (nTracerCntl%lXsAlign) THEN
!    CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
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

ENDIF

IF (nTracerCntl%lBoronSearch) THEN
  CALL UpdtBoronPPM(nTracerCntl%target_eigv, eigv, nTracerCntl%BoronPPM, .TRUE., MASTER)
  IF(nTracerCntl%lXslib) THEN
    CALL SetBoronCoolant(Core, FmInfo%Fxr, nTracerCntl%boronppm, PE%myzb, PE%myze)
    IF (nTracerCntl%lMacro) CALL SetCoreMacXs(Core, FmInfo)
    IF (nTracerCntl%lXsAlign) THEN
      CALL UpdateBlockFxr(FmInfo%Fxr, GroupInfo, PPMFeed)
      CALL XsLinIntPolTemp(nTracerCntl%ScatOd)
    END IF
  ELSE IF(nTracerCntl%libtyp .EQ. 11) THEN
    CALL SetBoronCoolant_NEACRP(nTracerCntl%boronppm)
  END IF
ENDIF

! IF(nTracerCntl%lfxrlib)THEN
!   CALL ReadFSRXs(Core,THInfo, FmInfo, CmInfo, GroupInfo, nTracerCntl)
! ENDIF

! IF(nTracerCntl%lfsrXS)THEN
!   write(*,*) '                  entering CEA_XS ... BYS edit 16/06/09'
!   CALL PrintFSRXs(Core,THInfo, FmInfo, CmInfo, GroupInfo, nTracerCntl)
!   CALL PrintFXRXs(Core,THInfo, FmInfo, CmInfo, GroupInfo, nTracerCntl)
!   CALL PrintFXRXs_CEA(Core,THInfo, FmInfo, CmInfo, GroupInfo, nTracerCntl, PE)
!   write(*,*) '                   exiting CEA_XS ... BYS edit 16/06/09'
! ENDIF
IF (nTracerCntl%lXeDyn) THEN
  CALL XeDynRelaxContrl(.TRUE.)
END IF

#ifdef HISTORY_OUT
CALL PowerHistoryOutputEdit(TRUE)
#endif

DO i = 1, ItrCntl%OuterMax

#ifdef MPI_ENV
  CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif

  ItrCntl%SrcIt = ItrCntl%SrcIt + 1
  IF (nTracerCntl%lMOCUR .AND. RTMASTER) CALL MOCUnderRelaxationFactor(Core, FmInfo, CmInfo, GroupInfo, nTracerCntl, PE)
  IF (lMOC) THEN
    CALL MOCDriverSwitch()
  ENDIF

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
!        CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
        CALL CUDASubGrpEffXSGen(core,FmInfo%fxr,nTracerCntl,PE)
      ELSE
        CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, Eigv, GroupInfo, nTracerCntl, PE)
      ENDIF
#else
      CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
#endif
    END IF
    IF (nTracerCntl%lMacro .AND. .NOT. nTracerCntl%lCusping_MPI) CALL SetCoreMacXs(Core, FmInfo)
  ENDIF

  CALL CmfdDriverSwitch()
#ifdef HISTORY_OUT
  CALL PowerHistoryOutputEdit(FALSE)
#endif

  if(i == ItrCntl%OuterMax) then
    last_TH=.true.
  endif

  ! TH Update
  IF (nTracerCntl%lFeedBack .AND. .NOT. lThConv) THEN

    IF (nTracerCntl%lSSPH) nTracerCntl%lSSPHreg=.TRUE.

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
        lSubGrp = .TRUE.
      ENDIF
      IF (nTRACERCntl%lPSM) THEN
        CALL GetfresoFXR(Core, FmInfo%Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
      ELSE
#ifdef __PGI
        IF (nTracerCntl%lXsAlign) THEN
!          CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
          CALL CUDASubGrpEffXSGen(core,FmInfo%fxr,nTracerCntl,PE)
        ELSE
        CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, Eigv, GroupInfo, nTracerCntl, PE)
      ENDIF
#else
        CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
#endif
      END IF
    ENDIF

    IF (nTracerCntl%lMacro .AND. .NOT. nTracerCntl%lCusping_MPI) CALL SetCoreMacXs(Core, FmInfo)

  ENDIF

! #ifdef inflow
!   IF (RTMASTER) CALL InflowTrXS(Core, FmInfo, CmInfo, eigv, GroupInfo, nTracerCntl, PE)
! #endif

  IF (nTracerCntl%lSearch) THEN
    IF (abs(eigv - nTracerCntl%target_eigv) / nTracerCntl%target_eigv .LT. 1.0E-5) lSearchConv = .TRUE.
    IF (nTracerCntl%lBoronSearch .AND. nTracerCntl%BoronPPM .EQ. 0._8) lSearchConv = .TRUE.
  ENDIF

  IF(nTracerCntl%lFeedBack .AND. ThInfo%TdopChg .LT. 5.e-5_8) lThConv = TRUE

  IF (MASTER) THEN
    WRITE(mesg, '(2x, A,L3, 3x, A, L3, 3x, A, L3)') 'Convergence - Flux :', ItrCntl%lconv, 'TH :',lThConv, 'B-Search:', lSearchConv
    CALL message(io8, FALSE, TRUE, mesg)
  ENDIF

  IF(ItrCntl%lconv .AND. lSearchConv) THEN
#ifdef __GAMMA_TRANSPORT
    IF (nTracerCntl%lGamma) THEN
      WRITE(mesg, '(a12, i1, a31)') 'Performing P', nTracerCntl%GammaScatOd, ' Gamma Transport Calculation...'
      IF (Master) CALL message(io8, TRUE, TRUE, mesg)
      ItrCntl%lGammaConv = FALSE
      DO WHILE (.NOT. ItrCntl%lGammaConv)
#ifdef __INTEL_MKL
        IF (Core%nxy .NE. 1) CALL GammaCMFDAcc(Core, CmInfo, FmInfo)
#endif
        CALL GammaMOCSweep(Core, RayInfo, FmInfo, PE, nTracerCntl, ItrCntl)
      ENDDO
    ENDIF
#endif
    IF (.NOT. nTracerCntl%lFeedBack) EXIT
    IF (nTracerCntl%lFeedBack .AND. lThConv) EXIT
  ENDIF

  IF (nTracerCntl%lSearch) THEN
!  IF (nTracerCntl%lSearch .AND. (.NOT. lSearchConv)) THEN
    IF (nTracerCntl%lBoronSearch) THEN
      CALL UpdtBoronPPM(nTracerCntl%target_eigv, eigv, nTracerCntl%BoronPPM, .FALSE., MASTER)
      IF(nTracerCntl%lXslib) THEN
        CALL SetBoronCoolant(Core, FmInfo%Fxr, nTracerCntl%boronppm, PE%myzb, PE%myze)
        IF (nTracerCntl%lXsAlign) THEN
          CALL UpdateBlockFxr(FmInfo%Fxr, GroupInfo, PPMFeed)
          CALL XsLinIntPolTemp(nTracerCntl%ScatOd)
        END IF
      ELSE IF(nTracerCntl%libtyp .EQ. 11) THEN
        CALL SetBoronCoolant_NEACRP(nTracerCntl%boronppm)
      END IF
      IF (nTracerCntl%lMacro .AND. .NOT. nTracerCntl%lCusping_MPI) CALL SetCoreMacXs(Core, FmInfo)
    ENDIF
    IF (nTracerCntl%lCrSearch .AND. i .GT. 1) THEN
      CALL CntlRodSearch(Core, FmInfo, CmInfo, eigv, GroupInfo, GcGroupInfo, nTracerCntl, PE)
      CALL SetBoronCoolant(Core, FmInfo%Fxr, nTracerCntl%boronppm, PE%myzb, PE%myze)
    ENDIF
  ENDIF

ENDDO

IF (nTracerCntl%lXeDyn) THEN
  CALL XeDynRelaxContrl(.FALSE.)
END IF
xkconv = eigv

IF (nTracerCntl%lCritSpec .AND. RTMASTER) THEN
  CALL GetCriticalSpectrum(Core, FmInfo, GroupInfo, nTracerCntl, PE)
ENDIF

IF(nTracerCntl%lAdjoint) THEN
  IF(MASTER) THEN
    IF (Master) CALL ShowHbar1(io8)
    WRITE(mesg, *) 'Start to Calculate Adjoint Flux at Initial State...'
    CALL message(io8, TRUE, TRUE, mesg)
  END IF
  ItrCntl%lConv = FALSE
  IF(nTracerCntl%lCMFDAdj) THEN
    CALL CMFD_adj_Switch()
  END IF
END IF

IF (nTracerCntl%lproblem .EQ. lsseigv .AND. RTMASTER) THEN

  xkcrit = xkconv
  IF (nTracerCntl%lXsLib) xkcrit = XsKinf(Core, FmInfo, GroupInfo, nTracerCntl%lXsLib, nTracerCntl%lCritSpec, PE)
  WRITE(mesg, '(A,2(F9.3,A),2(A,F9.5),A,F8.2,A)') 'Burnup = ',DeplCntl%NowBurnUp(2),                                &
              ' MWD/kgHM',DeplCntl%NowBurnUp(1),' EPFD','  keff = ',xkconv,'  kinf(BK) = ',                         &
              xkcrit,'  ppm = ',nTracerCntl%boronppm,' ppm'
  IF (MASTER) CALL message(io8, FALSE, TRUE, MESG)

ENDIF

!Print Out Eigenvalue
IF (nTracerCntl%lproblem .EQ. lsseigv .OR. nTracerCntl%lproblem .EQ. lEFTsearch ) THEN
  IF (Master) CALL ShowHbar1(io8)
  write(mesg, '(A27, 5X, 1F8.5)')  'k-eff =', Eigv
  IF (Master) CALL message(io8, FALSE, TRUE, MESG)

  IF(nTracerCntl%lAdjoint) THEN
    WRITE(mesg, '(A27, 5X, 1F8.5)')  'Adjoint k-eff =', Eigv_adj
    IF(Master) CALL message(io8, FALSE, TRUE, MESG)
  ENDIF
  IF (PE%lCmfdGrp) THEN
    IF (nTracerCntl%lHex) THEN
      CALL HexOutputEdit()
    ELSE
      CALL OutputEdit()
    END IF
  END IF
ENDIF
IF (nTracerCntl%lproblem .EQ. lTransient) THEN
  IF (Master) CALL ShowHbar1(io8)
  write(mesg, '(A27, 5X, 1F8.5)')  'k-eff =', Eigv
  IF (Master) CALL message(io8, FALSE, TRUE, MESG)

  IF(nTracerCntl%lAdjoint) THEN
    WRITE(mesg, '(A27, 5X, 1F8.5)')  'Adjoint k-eff =', Eigv_adj
    IF(Master) CALL message(io8, FALSE, TRUE, MESG)
  ENDIF
ENDIF
IF (nTracerCntl%lproblem .EQ. lBranch) THEN
  IF (Master) CALL ShowHbar1(io8)
  SELECT CASE(CurrentVarId)
  CASE(branch_tmod)
    write(mesg, '(A27, 5X, 1F8.5,a,i2,a,i2)')  'branch k-eff =', Eigv, '  Tmod  CaseId : ',CurrentCaseId, ' /', ntmod
  CASE(branch_tfuel)
    write(mesg, '(A27, 5X, 1F8.5,a,i2,a,i2)')  'branch k-eff =', Eigv, '  Tfuel CaseId : ',CurrentCaseId, ' /', ntfuel
  CASE(branch_boron)
    write(mesg, '(A27, 5X, 1F8.5,a,i2,a,i2)')  'branch k-eff =', Eigv, '  Boron CaseId : ',CurrentCaseId, ' /', nboron
  CASE(branch_rho)
    write(mesg, '(A27, 5X, 1F8.5,a,i2,a,i2)')  'branch k-eff =', Eigv, '  Rho   CaseId : ',CurrentCaseId, ' /', nrho
  ENDSELECT
  IF (Master) CALL message(io8, FALSE, TRUE, MESG)
  IF (PE%lCmfdGrp) THEN
    IF (nTracerCntl%lHex) THEN
      CALL HexOutputEdit()
    ELSE
      CALL OutputEdit()
    END IF
  END IF
ENDIF

!--- BYS edit
IF(nTracerCntl%lGC .AND. (lEFT_GCGEN.OR..NOT.nTracerCntl%lEFTSearch)) THEN
    WRITE(mesg, '(A)') 'Generating Group Constants ...'
    CALL message(io8, TRUE, TRUE, mesg)
    IF (nTracerCntl%lTranON) THEN
    CALL AllocTransient()
    CALL KinParamGen(Core, FmInfo, TranInfo, ThInfo, GroupInfo, .TRUE., nTracerCntl, PE)
    IF(nTRACERCntl%lXsLib) THEN
      CALL InitLambdaLib(TranInfo%Lambda, nTracerCntl%refdcy_del, nTracerCntl%llibdcy_del, nTracerCntl%lfitbeta)
      IF(nTracerCntl%lchidkgen) THEN
        CALL InitChidkLib(Core, FmInfo, GroupInfo, PE, nTracerCntl, TranInfo%chid, TranInfo%Chidk, ng, TranInfo%nprec)
      ELSE
        CALL InitChidLib(Core, FmInfo, GroupInfo, PE, nTracerCntl, TranInfo%Chid, ng)
      END IF
    ELSE
      !Chid, Lambda
       CALL ChidBen(TranInfo%Chid)
       CALL DnpLambdaBen(TranInfo%Lambda)
    ENDIF
END IF
    IF(nTracerCntl%gc_spec .GE. 3 )THEN
        IF(RTmaster) CALL GroupConstGenMac(Core, FmInfo, CmInfo, GroupInfo, nTracerCntl, PE, ng)
    ELSE
        IF(RTmaster) CALL GCGen(Core, FmInfo, THInfo, CmInfo, GroupInfo, nTracerCntl, PE, ng)
    ENDIF
    WRITE(mesg, '(A)') 'Generating Group Constants ... END'
    CALL message(io8, TRUE, TRUE, mesg)
ENDIF
!--- BYS edit end

CONTAINS

!--- CNJ Edit : Conditionally Compiled Drivers

SUBROUTINE SubgrpDriverSwitch

IF (nTRACERCntl%lPSM) THEN
  CALL PSMEffXSGen(Core, FmInfo%Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
ELSE
#ifdef __PGI
  IF (PE%lCUDA) THEN
    IF (nTracerCntl%lXsAlign) THEN
      CALL CUDASubGrpFsp_cuMLG(Core, FmInfo, THInfo, RayInfo, GroupInfo)
!      CALL CUDASubGrpFsp(Core, FmInfo, THInfo, RayInfo, GroupInfo)
    ELSE
      CALL CUDASubGrpFsp(Core, FmInfo, THInfo, RayInfo, GroupInfo)
    END IF
    GOTO 1
  ENDIF
#endif
  CALL SubGrpFsp(Core, FmInfo%Fxr, THInfo, RayInfo, GroupInfo, nTracerCntl, PE)
END IF
1 CONTINUE

END SUBROUTINE

SUBROUTINE MOCDriverSwitch
#ifdef __PGI
USE IEEE_ARITHMETIC
#endif

!IF(nTracerCntl%lMacro .AND. nTracerCntl%lCusping_MPI) THEN
!  CALL GetNeighborMocFlux(FmInfo%phis, FmInfo%neighphis, Core%nCoreFsr, PE%myzb, PE%myze, 1, GroupInfo%ng, Core%nz, Core%AxBC)
!  CALL SetCoreMacXS_Cusping(Core, FmInfo)
!END IF

#ifdef __PGI
IF (PE%lCUDA) THEN
  IF (nTracerCntl%lXsAlign) THEN
    CALL CUDAMOCSweep_cuMac(Core, RayInfo, FMInfo, eigv)
!    CALL CUDAMOCSweep(Core, RayInfo, FMInfo, eigv)
  ELSE
    CALL CUDAMOCSweep(Core, RayInfo, FMInfo, eigv)
  END IF
  GOTO 2
ENDIF
#endif
CALL MOCSweep(Core, RayInfo, FMInfo, eigv, ng, PE, nTracerCntl, ItrCntl)

2 CONTINUE

#ifdef __PGI
IF (ieee_is_nan(eigv)) CALL terminate("The Calculation Has Diverged!")
#else
IF (isnan(eigv)) CALL terminate("The Calculation Has Diverged!")
#endif

END SUBROUTINE

SUBROUTINE CmfdDriverSwitch
#ifdef __PGI
USE IEEE_ARITHMETIC
#endif

IF (lCMFD .AND. lCmfdGrp) THEN

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
ENDIF
#endif
#ifdef __INTEL_MKL
IF (PE%lMKL) THEN
  IF (.NOT. mklCntl%lDavidson) THEN
    CALL MKL_CmfdPower(Core, CmInfo, FmInfo, eigv)
  ELSE
    CALL MKL_CmfdDavidson(Core, CmInfo, FmInfo, eigv)
  ENDIF
  GOTO 3
ENDIF
#endif
IF (nTracerCntl%lHex) CALL terminate("HEX MUST USE MKL CMFD!")
CALL CmfdAcc(Core, CMInfo, FMInfo, THInfo, eigv, ng, TRUE, FALSE, PE, nTracerCntl, ItrCntl)

3 CONTINUE

#ifdef __PGI
IF (ieee_is_nan(eigv)) CALL terminate("The Calculation Has Diverged!")
#else
IF (isnan(eigv)) CALL terminate("The Calculation Has Diverged!")
#endif

ENDIF

END SUBROUTINE

SUBROUTINE CMFD_adj_Switch
#ifdef __PGI
USE IEEE_ARITHMETIC
#endif

IF (lCMFD .AND. lCmfdGrp) THEN

#ifdef __PGI
IF (PE%lCUDACMFD) THEN
  CALL CUDACmfd_Adj(Core, CMInfo, eigv_adj); GOTO 3
ENDIF
#endif
#ifdef __INTEL_MKL
IF (PE%lMKL) THEN
  CALL terminate("CMFD Adjoint Flux Calculation with MKL is under construction")
  GOTO 3
ENDIF
#endif
CALL terminate("CMFD Adjoint Flux Calculation is under construction")

3 CONTINUE

#ifdef __PGI
IF (ieee_is_nan(eigv_adj)) CALL terminate("The Calculation Has Diverged!")
#else
IF (isnan(eigv_adj)) CALL terminate("The Calculation Has Diverged!")
#endif

ENDIF

END SUBROUTINE

END SUBROUTINE
