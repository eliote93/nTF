#include <defines.h>

SUBROUTINE SSEIG()
USE PARAM
USE GEOM,             ONLY : Core, ng
USE RAYS,             ONLY : RayInfo
USE Core_mod,         ONLY : CMInfo, FmInfo, THInfo, GroupInfo, GcGroupInfo,   &
                             eigv, peigv, xkconv, xkcrit
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
USE Boron_mod,        ONLY : UpdtBoronPPM,           SetBoronCoolant, SetBoronCoolantCTF
USE XeDyn_Mod,        ONLY : UpdtXeDyn
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
#endif
USE MOC_COMMON,       ONLY : SetCoreMacXs
use SubChCoupling_mod,    only: last_TH
USE EFT_MOD

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

IF (nTracerCntl%lProblem .EQ. lDepletion .or. nTracerCntl%lXeDyn .or. (nTracerCntl%BoronPPM.GT.0._8)) THEN
   IF (nTracerCntl%lSSPH) nTracerCntl%lSSPHreg=.TRUE.
ENDIF

#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif

! Initialize Iteration Variables
CALL InitIterVar(Core, FMInfo, CmInfo, GroupInfo, .TRUE., ItrCntl, nTracerCntl, PE)

IF (Core%CrInfo%lCrChg) THEN
  CALL SetCrBankPosition(Core, FmInfo, CmInfo, GroupInfo, nTracerCntl, PE)
ENDIF

! Set Boron Concentration
IF (nTracerCntl%lXsLib .AND. nTracerCntl%BoronPPM .GT. 0._8) THEN
  WRITE(mesg,'(a,f12.2,a)') ' =======  Applied Boron PPM : ',nTracerCntl%boronppm, ' PPM'
  IF(master) CALL message(io8, TRUE, TRUE, mesg)
  CALL SetBoronCoolant(Core, FmInfo%Fxr, nTracerCntl%boronppm, PE%myzb, PE%myze)
ENDIF
WRITE(mesg, '(A, L3)') 'MOC Under-relaxtion Option : ', nTracerCntl%lMOCUR
IF(master) CALL message(io8,TRUE, TRUE, mesg)
WRITE(mesg, '(A, L3)') 'Axial Reflector FDM Option : ', nTracerCntl%lAxRefFDM
IF(master) CALL message(io8,TRUE, TRUE, mesg)

IF (nTracerCntl%lFeedBack) THEN

  lThConv = FALSE
  IF (nTracerCntl%lSSPH) nTracerCntl%lSSPHreg=.TRUE.

  CALL CmfdDriverSwitch()
  
#ifdef __GAMMA_TRANSPORT
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

  IF (nTRACERCntl%lDeplVarInit) CALL InitIterVar(Core, FMInfo, CmInfo, GroupInfo, .TRUE., ItrCntl, nTracerCntl, PE)

  ! Update Xenon Dynamics
  IF (nTracerCntl%lXeDyn) THEN
    CALL UpdtXeDyn(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, DeplCntl, PE)
  ENDIF

  ! Generate Self-Shielded XS
  IF (nTracerCntl%lrestrmt) CALL SubgrpDriverSwitch()

#ifdef MPI_ENV
  CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif

  ! Effective XS Calculation
  CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, Eigv, GroupInfo, nTracerCntl, PE)
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
  ENDIF

  ! Generate Self-Shielded XS
  IF (nTracerCntl%lrestrmt) THEN
    CALL SubgrpDriverSwitch()
    lSubGrp = .TRUE.
  ENDIF

  ! Effective XS Calculation
  CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, Eigv, GroupInfo, nTracerCntl, PE)
  IF (nTracerCntl%lMacro) CALL SetCoreMacXs(Core, FmInfo)
  
  CALL CmfdDriverSwitch()

ENDIF

IF (nTracerCntl%lBoronSearch) THEN
  CALL UpdtBoronPPM(nTracerCntl%target_eigv, eigv, nTracerCntl%BoronPPM, .FALSE., MASTER)
  CALL SetBoronCoolant(Core, FmInfo%Fxr, nTracerCntl%boronppm, PE%myzb, PE%myze)
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

DO i = 1, ItrCntl%OuterMax

#ifdef MPI_ENV
  CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif

  ItrCntl%SrcIt = ItrCntl%SrcIt + 1
  IF (nTracerCntl%lMOCUR .AND. RTMASTER) CALL MOCUnderRelaxationFactor(Core, FmInfo, CmInfo, GroupInfo, nTracerCntl, PE)
  IF (lMOC) CALL MOCDriverSwitch()
  
  IF (nTracerCntl%lXeDyn) THEN
    CALL UpdtXeDyn(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, DeplCntl, PE)
    CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, Eigv, GroupInfo, nTracerCntl, PE)
  ENDIF
  
  CALL CmfdDriverSwitch()

  if(i == ItrCntl%OuterMax) then
    last_TH=.true.
  endif

  ! TH Update
  IF (nTracerCntl%lFeedBack .AND. .NOT. lThConv) THEN

    IF (nTracerCntl%lSSPH) nTracerCntl%lSSPHreg=.TRUE.

    CALL SteadyStateTH(Core, CmInfo, FmInfo, ThInfo, Eigv, ng, GroupInfo, nTracerCntl, PE)
    CALL THFeedBack(Core, CmInfo, FmInfo, ThInfo, nTracerCntl, PE)

#ifdef MPI_ENV
    CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif

    IF (nTracerCntl%lrestrmt) THEN
      IF (ThInfo%TdopChg .LT. 1.D-03 .AND. .NOT. lSubGrp) THEN
        CALL SubgrpDriverSwitch()
        lSubGrp = .TRUE.
      ENDIF
      CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, Eigv, GroupInfo, nTracerCntl, PE)
    ENDIF

    IF (nTracerCntl%lMacro) CALL SetCoreMacXs(Core, FmInfo)
    
  ENDIF

! #ifdef inflow
!   IF (RTMASTER) CALL InflowTrXS(Core, FmInfo, CmInfo, eigv, GroupInfo, nTracerCntl, PE)
! #endif

  IF (nTracerCntl%lSearch) THEN
    IF (abs(eigv - nTracerCntl%target_eigv) / nTracerCntl%target_eigv .LT. 5.0E-5) lSearchConv = .TRUE.
    IF (nTracerCntl%lBoronSearch .AND. nTracerCntl%BoronPPM .EQ. 0._8) lSearchConv = .TRUE.
  ENDIF

  IF(nTracerCntl%lFeedBack .AND. ThInfo%TdopChg .LT. 5.e-4_8) lThConv = TRUE

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
        IF (Core%nxy .NE. 1) CALL GammaCMFDAcc(Core, CmInfo, FmInfo)
        CALL GammaMOCSweep(Core, RayInfo, FmInfo, PE, nTracerCntl, ItrCntl)
      ENDDO
    ENDIF
#endif
    IF (.NOT. nTracerCntl%lFeedBack) EXIT
    IF (nTracerCntl%lFeedBack .AND. lThConv) EXIT
  ENDIF

  IF (nTracerCntl%lSearch) THEN
    IF (nTracerCntl%lBoronSearch) THEN
      CALL UpdtBoronPPM(nTracerCntl%target_eigv, eigv, nTracerCntl%BoronPPM, .FALSE., MASTER)
      CALL SetBoronCoolant(Core, FmInfo%Fxr, nTracerCntl%boronppm, PE%myzb, PE%myze)
    ENDIF
    IF (nTracerCntl%lCrSearch .AND. i .GT. 1) THEN
      CALL CntlRodSearch(Core, FmInfo, CmInfo, eigv, GroupInfo, GcGroupInfo, nTracerCntl, PE)
      CALL SetBoronCoolant(Core, FmInfo%Fxr, nTracerCntl%boronppm, PE%myzb, PE%myze)
    ENDIF
  ENDIF
  
ENDDO

xkconv = eigv

IF (nTracerCntl%lCritSpec .AND. RTMASTER) THEN
  CALL GetCriticalSpectrum(Core, FmInfo, GroupInfo, nTracerCntl, PE)
ENDIF

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
  IF (PE%lCmfdGrp) THEN
    IF (nTracerCntl%lHex) THEN
      CALL HexOutputEdit()
    ELSE
      CALL OutputEdit()
    END IF
  END IF
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

#ifdef __PGI
IF (PE%lCUDA) THEN
  CALL CUDASubGrpFsp(Core, FmInfo, THInfo, RayInfo, GroupInfo); GOTO 1
ENDIF
#endif
CALL SubGrpFsp(Core, FmInfo%Fxr, THInfo, RayInfo, GroupInfo, nTracerCntl, PE)

1 CONTINUE
  
END SUBROUTINE

SUBROUTINE MOCDriverSwitch
#ifdef __PGI
USE IEEE_ARITHMETIC
#endif

#ifdef __PGI
IF (PE%lCUDA) THEN
  CALL CUDAMOCSweep(Core, RayInfo, FMInfo, eigv); GOTO 2
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

#ifdef __PGI
IF (PE%lCUDACMFD) THEN
  CALL CUDACmfdAcc(Core, CMInfo, FMInfo, eigv); GOTO 3
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

END SUBROUTINE