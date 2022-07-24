#include <defines.h>
SUBROUTINE CmfdAcc_new(Core, CMInfo, FMInfo, THInfo, eigv, ng, lcmfd, lreset, PE, nTracerCntl, ItrCntl)
USE PARAM
USE TYPEDEF,        ONLY : PE_TYPE,         CoreInfo_Type,      Pin_Type,                   &
                           CMInfo_Type,     FmInfo_Type,        ThInfo_Type,                &
                           FxrInfo_Type,    CmfdLs_TYPE,        PinXs_Type,                 &
                           AxFlx_Type
USE CMFD_mod,       ONLY : CmfdPinXS,        CmfdLS,                                        &
                           PhiC1g,            SRC,                                          &
                           SetCMFDEnviorment, HomoXsGen,         RadCouplingCoeffGen,       &
                           UpdatePhiC,        ConvertSubPlnPhi,                             &
                           CmfdSrcUpdt,       CmfdPsiUpdt,                                  &
                           CmfdEigUpdate,     ResidualError,     MOCSolUpdt,                &
                           MOCPhiInUpdt,      MOCLinSrcUpdt,     HomBuckling
USE CORE_MOD,       ONLY : GroupInfo,         GcGroupInfo
!USE AxSolver_mod,   ONLY : AxialSolver,       AxSrcUpdate
USE BiLU_Mod
USE MpiAxSolver_mod,ONLY : MpiAxialSolver,    AxSrcUpdate,       VoidLkgCorrection, MpiAxialSolverNew
#ifdef MPI_ENV
USE MpiComm_mod,    ONLY : MPI_SYNC
#endif
USE Boron_Mod,      ONLY : UpdtBoronPPM,      UpdtBoronCmfdXS
USE BiCGSTAB_mod,   ONLY : BiCGSTAB
USE SUBGRP_MOD,     ONLY : FxrChiGen
USE BasicOperation, ONLY : CP_CA, CP_VA
USE IOUTIL,         ONLY : message, terminate
USE FILES,          ONLY : io8
USE CNTL,           ONLY : nTracerCntl_Type
USE itrcntl_mod,    ONLY : ItrCntl_TYPE, CMFDItrCntl_TYPE
use timer,          only : nTracer_dclock, TimeChk
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CMInfo_Type) :: CMInfo
TYPE(FMInfo_Type) :: FmInfo
TYPE(THInfo_Type) :: ThInfo
TYPE(PE_TYPE) :: PE
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(ItrCntl_TYPE) :: ItrCntl
REAL :: eigv
REAL :: CmfdTimeBeg, CmfdTimeEnd
INTEGER :: ng
LOGICAL :: lreset, lcmfd, lDtilOnly
LOGICAL, SAVE :: lfirst

!Pointing Varibables
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Fxrinfo_type), POINTER :: Fxr(:,:)
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
!TYPE(AxFlx_Type), POINTER :: AxFlx(:, :)
REAL, POINTER :: PHIS(:,:,:)
REAL, POINTER :: PhiC(:, :, :), PSIC(:,:), PSICD(:,:)
REAL, POINTER :: PhiFm(:, :, :), PsiFm(:, :), PsiFmD(:, :)
REAL, POINTER :: Mmt(:, :, :, :), MmtFm(:, :, :, :)
REAL, POINTER :: RadJout(:, :, :, :, :)
!
INTEGER :: ig, l, ITER, jsweep
INTEGER :: GrpBeg, GrpEnd, nGroupInfo
INTEGER :: myzb, myze, myzbf, myzef, nxy

REAL :: psierr, eigerr, peigv, ResErr, ResErr0, ResErr00, ResErrN, ResErrN0
LOGICAL :: lXsLib, lscat1, lsubplane, laxial
!Iteration COntrol
INTEGER :: nitermax, nitermin, AxSolver, ncmfdpnodal
INTEGER :: nAxUpdt, iterCMFDnodal
REAL :: convcrit, nodalcrit, nodaloffcrit, nodalURcrit
LOGICAL :: lconv, lAxDhatConv, lExit, lDhat, lsigt, lnodal, lur
LOGICAL :: CmfdMaster, CmfdSlave
DATA lfirst /TRUE/

CMFDMaster = PE%CMFDMaster; CMFDSlave = PE%CMFDSlave
!CMFD Time Check
CmfdTImeBeg = nTracer_dclock(FALSE, FALSE)

!Pointing
Pin => Core%Pin;
Fxr => FmInfo%Fxr; phis => FmInfo%phis
PinXS => CMInfo%PinXS
PhiC => CMInfo%PhiC; PhiFM => CMInfo%PhiFM
PsiC => CMInfo%PsiC; PsiCD => CMInfo%PsiCD
PsiFm => CMInfo%PsiFm; PsiFmD => CMInfo%PsiFmD
RadJout => CMinfo%RadJout !ninout/ 1:in 2:out 3:surfphi ! BYS edit 16/02/11

nitermax = itrcntl%CMFDItrCntl%nitermax
nitermin = itrcntl%CMFDItrCntl%nitermin
convcrit = itrcntl%CMFDItrCntl%convcrit
ncmfdpnodal=itrcntl%CMFDItrCntl%ncmfdpnodal


nxy = Core%nxy; myzb = PE%myzb; myze = PE%myze
myzbf = PE%myzbf; myzef = PE%myzef
lXsLib = nTracerCntl%lXsLib; lscat1 = nTracerCntl%lScat1
lsubplane = nTracerCntl%lSubPlane
AxSolver = nTracerCntl%AxSolver

lDhat = .TRUE.
IF(ItrCntl%Cmfdit .EQ. ItrCntl%Cmfdit0) lDhat = .FALSE.

!IF(ItrCntl%MocIt .LT. 1) AxSolver = LP1SENM
IF(lfirst .OR. lReset) THEN
  CALL SetCMFDEnviorment(Core, CmInfo, ng, PE)
  CALL SetCmfdMpiOmpEnv(Core, CMInfo%CoreCmfdLs , ng, PE)
  lfirst = FALSE
ENDIF
!Homogenization ProCedure
WRITE(mesg,'(a)') 'Cell Homogenization (H)...'
IF(CMFDMaster) CALL message(io8, TRUE, TRUE, mesg)
CALL FxrChiGen(Core, Fxr, FmInfo, GroupInfo, nTracerCntl, myzb, myze)
lsigT=FALSE
!IF( ItrCntl%Axit .LT. 1 )THEN
!    lsigT=TRUE
!    IF(CMFDMaster) WRITE(*,*) 'SigT Diffusion'
!ENDIF
CALL HomoXsGen(Core, Fxr, Phis, PinXs, myzb, myze, ng, lXsLib, lScat1, lsigT)
#ifdef Buckling
IF(nTracerCntl%lBsq) THEN
  CALL HomBuckling(Core, Fxr, Phis, PinXs, myzb, myze, ng, nTracerCntl%bsq, lXsLib)
ENDIF
#endif
CALL RadCouplingCoeffGen(Core, PinXS, RadJout, ng, lDhat, PE)

WRITE(mesg, '(a)') 'Performing CMFD Calculation...'
IF(CMFDMaster) CALL message(io8, TRUE, TRUE, mesg)

!Update Homo PHI
CALL UpdatePhiC(PinXs, PhiC)
CALL ConvertSubPlnPhi(PhiC, PhiFm, 1)

    nAxUpdt = 0;
    ResErr=10._8
    NodalCrit=1E-1;
    NodalCrit=1E-3;
    NodalURCrit=1E-4;
    nodaloffcrit=1E-4;
    nodaloffcrit=1E-100;
iterCMFDnodal=0
!Axail Dhat Update
#ifndef MPI_ENV
IF(nTracerCntl%l3dim) THEN
  WRITE(mesg, '(2x, A32, 2x,A10)') 'Performing Axial Calculation : ',AxSolverName(AxSolver)
  CALL message(io8, FALSE, TRUE, mesg)
  lDtilOnly = .FALSE.;
  IF(ItrCntl%Cmfdit .EQ. 0) lDtilOnly = .TRUE.
  CALL AxialSolver(Core, CMInfo, Eigv, ng, AxSolver, lreset, lDtilOnly, PE)
  !nAxUpdt = 1;
ENDIF
#else
IF(nTracerCntl%l3dim) THEN
  lDtilOnly = .FALSE.;
  IF(ItrCntl%Cmfdit .EQ. ItrCntl%Cmfdit0) lDtilOnly = .TRUE.
  IF( ldtilonly )THEN
  !IF( ldtilonly .OR. itrcntl%lnodal )THEN
  !IF( ldtilonly .OR. ResErr .LT. NodalCrit)THEN
    !IF(.NOT.ldtilonly) CALL MpiAxialSolver(Core, CMInfo, GroupInfo, Eigv, ng, AxSolver, .TRUE., .TRUE., PE)
    CALL MpiAxialSolver(Core, CMInfo, GroupInfo, Eigv, ng, AxSolver, .TRUE., lDtilOnly, PE, nTracerCntl%lnonFuelpinFDM)
    !IF(CMFDMaster) WRITE(mesg, '(2x, A32, 2x,A10)') 'Performing Axial Calculation : ',AxSolverName(AxSolver)
    !IF(CMFDMaster) CALL message(io8, FALSE, TRUE, mesg)
    !CALL MpiAxialSolver(Core, CMInfo, GroupInfo, Eigv, ng, AxSolver, .TRUE., FALSE, PE)
    !ItrCntl%Axit = ItrCntl%Axit + 1
    !nAxUpdt = 1;
  ENDIF
ENDIF
#endif
!LINEAR SYSTEM
CALL SetCmfdLinearSystem(TRUE,  nTracerCntl%l3dim, AxSolver)

CALL MakeBiLU(Core, CmInfo%CoreCMFDLS(1:ng), PhiFm, ng, PE)

CALL CmfdPsiUpdt(PhiFm, PsiFm)

lconv = FALSE; lAxDhatConv = FALSE
nGroupInfo = 2;
IF(.NOT. GroupInfo%lUpScat) nGroupInfo = 1
ResErr0 =  ResidualError(phifm, psifm, eigv, 1, ng, PE)
lAxial=.FALSE.
DO iter = 1, nitermax
  ItrCntl%Cmfdit = ItrCntl%Cmfdit + 1
  iterCMFDnodal=iterCMFDnodal+1
  DO jsweep =1, nGroupInfo
    GrpBeg = 1; GrpEnd = ng
    IF(jsweep .GT. 1) THEN
      GrpBeg = GroupInfo%UpScatRange(1); GrpEnd = GroupInfo%UpScatRange(2)
    ENDIF
    DO ig = GrpBeg, GrpEnd
      !Source Update
      CALL CmfdSrcUpdt(Src, Psifm, phifm, eigv, ig)
      CALL CP_VA(Phic1g(1:nxy, myzbf : myzef), Phifm(1:nxy, myzbf : myzef, ig), nxy, myzef - myzbf + 1)
      CALL BiCGSTAB(CmfdLs(ig), Phic1g, SRC, itrcntl%InSolverItrCntl)
#ifndef MPI_ENV
      CALL CP_VA(Phifm(1:nxy, myzbf : myzef, ig), Phic1g(1:nxy, myzbf : myzef), nxy, myzef - myzbf + 1)
#else
      CALL CP_VA(Phifm(1:nxy, myzbf-1 : myzef+1, ig), Phic1g(1:nxy, myzbf-1 : myzef+1), nxy, myzef - myzbf + 3)
#endif
    ENDDO
  END DO

  CALL CP_VA(PsifmD(1 : nxy, myzbf : myzef), Psifm(1 : nxy, myzbf : myzef), nxy, myzef - myzbf + 1)
  !Fission Source Update
  CALL CmfdPsiUpdt(phifm, psifm)
  !EigenValue Update
  peigv = eigv
  CALL CmfdEigUpdate(psifm, psifmd, eigv, psierr, PE)
  !Error Estimation
  eigerr = (eigv - peigv)/eigv
  ResErr = ResidualError(phifm, psifm, eigv, 1, ng, PE)
  !Log Out
  IF(iter .eq. 1) ResErr0 = ResErr
  IF(iter .eq. 1) ResErr00 = ResErr
  IF( nAxUpdt .EQ. 1) ResErrN0=ResErr
  IF(itrcntl%lnodal.AND.lAxial)THEN
       ResErr00 = ResErr
       laxial=.false.
       IF( ResErr .LT. nodaloffcrit )THEN
        IF(CMFDMaster) WRITE(mesg,'(A)') 'Axial Nodal Off'
        IF(CMFDMaster) CALL message(io8, FALSE, TRUE, mesg)
        itrcntl%lnodal=FALSE
       ENDIF
  endif
  IF(CMFDMaster) WRITE(mesg,'(A9, I9, F22.6, 3x, F10.5, 1pe15.3)') 'MGOUTER', ItrCntl%Cmfdit, eigv, ResErr/ResErr0, ResErr
  IF(CMFDMaster) CALL message(io8, FALSE, TRUE, mesg)
  !Convergence Check
  lExit = .TRUE.; lConv = .FALSE.
  IF((ResErr/ResErr00) .lt. convcrit ) lconv = TRUE
  IF((ResErr/ResErr0) .gt. convcrit ) lconv = FALSE
  IF( nAxUpdt .GT. 0 )THEN
      IF((ResErrN/ResErrN0) .gt. convcrit ) lconv = FALSE
  ENDIF
  lExit = lExit .AND. lConv
  IF(iter .LE. nitermin) lExit = .FALSE.
  IF(iter .GE. nitermin .AND. ResErr .LT. 1.0E-7_8) lExit = .TRUE.
  !IF( nAxUpdt .LT. nCmfdpNodal) lexit=.FALSE.
  IF( itrcntl%lnodal .AND. nAxUpdt .LT. 1 ) lexit=.FALSE.
  IF(lExit) EXIT
  IF(nitermax .EQ. iter) EXIT
  !IF(mod(iter, ncmfdpnodal) .eq. 0 .and. nTracerCntl%l3dim .AND. ResErr .LT. NodalCrit) THEN
  IF( itrcntl%lnodal .AND. (iterCMFDnodal .GE. ncmfdpnodal .and. nTracerCntl%l3dim .AND. ResErr .LT. NodalCrit .AND. lconv)) THEN
  !IF(lConv .and. nTracerCntl%l3dim) THEN
    nAxUpdt = nAxUpdt  + 1
    IF(lAxDhatConv) EXIT
    IF(CMFDMaster) WRITE(mesg, '(2x, A32, 2x,A10)') 'Performing Axial Calculation : ',AxSolverName(AxSolver)
    IF(CMFDMaster) CALL message(io8, FALSE, TRUE, mesg)
    !Axail Dhat Update
 !     CALL MpiAxialSolver(Core, CMInfo, GroupInfo, Eigv, ng, AxSolver, lreset, .TRUE., PE)
      !CALL MpiAxialSolver(Core, CMInfo, GroupInfo, Eigv, ng, AxSolver, lreset, .FALSE., PE)  !old
      lUR=FALSE
      !IF( ResErr .GT. NodalURCrit ) lUR = TRUE
      CALL MpiAxialSolverNew(Core, CMInfo, GroupInfo, Eigv, ng, AxSolver, lreset, .FALSE., PE, lUR)   !test
    ItrCntl%Axit = ItrCntl%Axit + 1
    CALL SetCmfdLinearSystem(TRUE, TRUE, AxSolver)
    CALL MakeBiLU(Core, CmInfo%CoreCMFDLS(1:ng), PhiFm, ng, PE)
    laxial=.TRUE.
    iterCMFDnodal=0
    !EXIT
  ENDIF
  IF(mod(iter, 5) .EQ. 0 .and. nTracerCntl%lGcCmfd) THEN
    CALL GcCmfdAcc(Core, CmInfo, Eigv, TRUE, GroupInfo, GcGroupInfo, nTracerCntl, ItrCntl, PE)
    CALL CP_VA(PsifmD(1 : nxy, myzbf : myzef), Psifm(1 : nxy, myzbf : myzef), nxy, myzef - myzbf + 1)
    CALL CmfdPsiUpdt(phifm, psifm)
  ENDIF
END DO

CALL ConvertSubPlnPhi(PhiC, PhiFm, 2)

!Update Radial MOC Solution
CALL MOCSolUpdt(Core, FmInfo, CmInfo, myzb, myze, ng)
CALL MOCPhiInUpdt(Core, CmInfo, myzb, myze, ng)
!Axial Source For MOC
IF(nTracerCntl%l3dim) THEN
  CALL AxSrcUpdate(Core, CmInfo, myzb, myze, 1, ng, PE, AxSolver)
  IF(.NOT. nTracerCntl%lBenchXs)THEN ! benchmark XS bug fixed 160314
    CALL VoidLkgCorrection(Core, FmInfo, CmInfo, myzb, myze, 1, ng, PE) !150416 BYS edit off
  ENDIF
ENDIF


IF(nTracerCntl%lLinSrc) CALL MOCLinSrcUpdt(Core, PinXS, PHIC, FmInfo%LinSrcSlope, myzb, myze, ng)


CmfdTImeEnd = nTracer_dclock(FALSE, FALSE)
TimeChk%CmfdTime = TimeChk%CmfdTime + (CmfdTimeEnd - CmfdTimeBeg)
!Free
NULLIFY(Pin);
NULLIFY(Fxr); NULLIFY(phis)
NULLIFY(PinXS);

NULLIFY(PhiC); NULLIFY(PhiFM)
NULLIFY(PsiC); NULLIFY(PsiCD)
NULLIFY(PsiFm); NULLIFY(PsiFmD)
NULLIFY(MmtFm); NULLIFY(Mmt)
NULLIFY(RadJout)

END SUBROUTINE

SUBROUTINE CmfdAcc(Core, CMInfo, FMInfo, THInfo, eigv, ng, lcmfd, lreset, PE, nTracerCntl, ItrCntl)
USE PARAM
USE TYPEDEF,        ONLY : PE_TYPE,         CoreInfo_Type,      Pin_Type,                   &
                           CMInfo_Type,     FmInfo_Type,        ThInfo_Type,                &
                           FxrInfo_Type,    CmfdLs_TYPE,        PinXs_Type,                 &
                           AxFlx_Type
USE CMFD_mod,       ONLY : CmfdPinXS,        CmfdLS,                                        &
                           PhiC1g,            SRC,                                          &
                           SetCMFDEnviorment, HomoXsGen,         RadCouplingCoeffGen,       &
                           UpdatePhiC,        ConvertSubPlnPhi,                             &
                           CmfdSrcUpdt,       CmfdPsiUpdt,                                  &
                           CmfdEigUpdate,     ResidualError,     MOCSolUpdt,                &
                           MOCPhiInUpdt,      MOCLinSrcUpdt,     HomBuckling
USE CORE_MOD,       ONLY : GroupInfo,         GcGroupInfo
!USE AxSolver_mod,   ONLY : AxialSolver,       AxSrcUpdate
USE BiLU_Mod
USE MpiAxSolver_mod,ONLY : MpiAxialSolver,    AxSrcUpdate,       VoidLkgCorrection
#ifdef MPI_ENV
USE MpiComm_mod,    ONLY : MPI_SYNC
#endif
USE Boron_Mod,      ONLY : UpdtBoronPPM,      UpdtBoronCmfdXS
USE BiCGSTAB_mod,   ONLY : BiCGSTAB
USE SUBGRP_MOD,     ONLY : FxrChiGen
USE BasicOperation, ONLY : CP_CA, CP_VA
USE IOUTIL,         ONLY : message, terminate
USE FILES,          ONLY : io8
USE CNTL,           ONLY : nTracerCntl_Type
USE itrcntl_mod,    ONLY : ItrCntl_TYPE, CMFDItrCntl_TYPE
use timer,          only : nTracer_dclock, TimeChk
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CMInfo_Type) :: CMInfo
TYPE(FMInfo_Type) :: FmInfo
TYPE(THInfo_Type) :: ThInfo
TYPE(PE_TYPE) :: PE
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(ItrCntl_TYPE) :: ItrCntl
REAL :: eigv
REAL :: CmfdTimeBeg, CmfdTimeEnd
INTEGER :: ng
LOGICAL :: lreset, lcmfd, lDtilOnly
LOGICAL, SAVE :: lfirst

!Pointing Varibables
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Fxrinfo_type), POINTER :: Fxr(:,:)
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
!TYPE(AxFlx_Type), POINTER :: AxFlx(:, :)
REAL, POINTER :: PHIS(:,:,:)
REAL, POINTER :: PhiC(:, :, :), PSIC(:,:), PSICD(:,:)
REAL, POINTER :: PhiFm(:, :, :), PsiFm(:, :), PsiFmD(:, :)
REAL, POINTER :: Mmt(:, :, :, :), MmtFm(:, :, :, :)
REAL, POINTER :: RadJout(:, :, :, :, :)
!
INTEGER :: ig, l, ITER, jsweep
INTEGER :: GrpBeg, GrpEnd, nGroupInfo
INTEGER :: myzb, myze, myzbf, myzef, nxy

REAL :: psierr, eigerr, peigv, ResErr, ResErr0
LOGICAL :: lXsLib, lscat1, lsubplane
!Iteration COntrol
INTEGER :: nitermax, nitermin, AxSolver
INTEGER :: nAxUpdt
REAL :: convcrit
LOGICAL :: lconv, lAxDhatConv, lExit, lDhat
LOGICAL :: CmfdMaster, CmfdSlave
DATA lfirst /TRUE/

IF (nTracerCntl%lResetResErr) THEN
  CALL CmfdAcc_new(Core, CMInfo, FMInfo, THInfo, eigv, ng, TRUE, FALSE, PE, nTracerCntl, ItrCntl)
  RETURN
ENDIF

CMFDMaster = PE%CMFDMaster; CMFDSlave = PE%CMFDSlave
!CMFD Time Check
CmfdTImeBeg = nTracer_dclock(FALSE, FALSE)

!Pointing
Pin => Core%Pin;
Fxr => FmInfo%Fxr; phis => FmInfo%phis
PinXS => CMInfo%PinXS
PhiC => CMInfo%PhiC; PhiFM => CMInfo%PhiFM
PsiC => CMInfo%PsiC; PsiCD => CMInfo%PsiCD
PsiFm => CMInfo%PsiFm; PsiFmD => CMInfo%PsiFmD
RadJout => CMinfo%RadJout !ninout/ 1:in 2:out 3:surfphi ! BYS edit 16/02/11

nitermax = itrcntl%CMFDItrCntl%nitermax
nitermin = itrcntl%CMFDItrCntl%nitermin
convcrit = itrcntl%CMFDItrCntl%convcrit

nxy = Core%nxy; myzb = PE%myzb; myze = PE%myze
myzbf = PE%myzbf; myzef = PE%myzef
lXsLib = nTracerCntl%lXsLib; lscat1 = nTracerCntl%lScat1
lsubplane = nTracerCntl%lSubPlane
AxSolver = nTracerCntl%AxSolver

lDhat = .TRUE.
IF(ItrCntl%Cmfdit .EQ. ItrCntl%Cmfdit0) lDhat = .FALSE.

!IF(ItrCntl%MocIt .LT. 1) AxSolver = LP1SENM
IF(lfirst .OR. lReset) THEN
  CALL SetCMFDEnviorment(Core, CmInfo, ng, PE)
  CALL SetCmfdMpiOmpEnv(Core, CMInfo%CoreCmfdLs , ng, PE)
  lfirst = FALSE
ENDIF
!Homogenization ProCedure
WRITE(mesg,'(a)') 'Cell Homogenization (H)...'
IF(CMFDMaster) CALL message(io8, TRUE, TRUE, mesg)
CALL FxrChiGen(Core, Fxr, FmInfo, GroupInfo, nTracerCntl, myzb, myze)
CALL HomoXsGen(Core, Fxr, Phis, PinXs, myzb, myze, ng, lXsLib, lScat1, FALSE)
#ifdef Buckling
IF(nTracerCntl%lBsq) THEN
  CALL HomBuckling(Core, Fxr, Phis, PinXs, myzb, myze, ng, nTracerCntl%bsq, lXsLib)
ENDIF
#endif
CALL RadCouplingCoeffGen(Core, PinXS, RadJout, ng, lDhat, PE)

WRITE(mesg, '(a)') 'Performing CMFD Calculation...'
IF(CMFDMaster) CALL message(io8, TRUE, TRUE, mesg)

!Update Homo PHI
CALL UpdatePhiC(PinXs, PhiC)
CALL ConvertSubPlnPhi(PhiC, PhiFm, 1)

!Axail Dhat Update
#ifndef MPI_ENV
IF(nTracerCntl%l3dim) THEN
  WRITE(mesg, '(2x, A32, 2x,A10)') 'Performing Axial Calculation : ',AxSolverName(AxSolver)
  CALL message(io8, FALSE, TRUE, mesg)
  lDtilOnly = .FALSE.;
  IF(ItrCntl%Cmfdit .EQ. 0) lDtilOnly = .TRUE.
  CALL AxialSolver(Core, CMInfo, Eigv, ng, AxSolver, lreset, lDtilOnly, PE)
ENDIF
#else
IF(nTracerCntl%l3dim) THEN
  IF(CMFDMaster) WRITE(mesg, '(2x, A32, 2x,A10)') 'Performing Axial Calculation : ',AxSolverName(AxSolver)
  IF(CMFDMaster) CALL message(io8, FALSE, TRUE, mesg)
  lDtilOnly = .FALSE.;
  IF(ItrCntl%Cmfdit .EQ. ItrCntl%Cmfdit0) lDtilOnly = .TRUE.
  CALL MpiAxialSolver(Core, CMInfo, GroupInfo, Eigv, ng, AxSolver, .TRUE., lDtilOnly, PE, nTracerCntl%lnonFuelpinFDM)
ENDIF
#endif
!Preconditioning
!CALL SetCmfdLinearSystem(TRUE,  nTracerCntl%l3dim, AxSolver)

!CALL SetPcondiLinearSystem(PhiFm, TRUE,  nTracerCntl%l3dim, AxSolver, PE)

!LINEAR SYSTEM
CALL SetCmfdLinearSystem(TRUE,  nTracerCntl%l3dim, AxSolver)

CALL MakeBiLU(Core, CmInfo%CoreCMFDLS(1:ng), PhiFm, ng, PE)

CALL CmfdPsiUpdt(PhiFm, PsiFm)

lconv = FALSE; lAxDhatConv = FALSE
nGroupInfo = 2;
nAxUpdt = 1;
IF(.NOT. GroupInfo%lUpScat) nGroupInfo = 1
ResErr0 =  ResidualError(phifm, psifm, eigv, 1, ng, PE)

DO iter = 1, nitermax
  ItrCntl%Cmfdit = ItrCntl%Cmfdit + 1
  DO jsweep =1, nGroupInfo
    GrpBeg = 1; GrpEnd = ng
    IF(jsweep .GT. 1) THEN
      GrpBeg = GroupInfo%UpScatRange(1); GrpEnd = GroupInfo%UpScatRange(2)
    ENDIF
    DO ig = GrpBeg, GrpEnd
      !Source Update
      CALL CmfdSrcUpdt(Src, Psifm, phifm, eigv, ig)
      CALL CP_VA(Phic1g(1:nxy, myzbf : myzef), Phifm(1:nxy, myzbf : myzef, ig), nxy, myzef - myzbf + 1)
      CALL BiCGSTAB(CmfdLs(ig), Phic1g, SRC, itrcntl%InSolverItrCntl)
#ifndef MPI_ENV
      CALL CP_VA(Phifm(1:nxy, myzbf : myzef, ig), Phic1g(1:nxy, myzbf : myzef), nxy, myzef - myzbf + 1)
#else
      CALL CP_VA(Phifm(1:nxy, myzbf-1 : myzef+1, ig), Phic1g(1:nxy, myzbf-1 : myzef+1), nxy, myzef - myzbf + 3)
#endif
    ENDDO
  END DO

  CALL CP_VA(PsifmD(1 : nxy, myzbf : myzef), Psifm(1 : nxy, myzbf : myzef), nxy, myzef - myzbf + 1)
  !Fission Source Update
  CALL CmfdPsiUpdt(phifm, psifm)
  !EigenValue Update
  peigv = eigv
  CALL CmfdEigUpdate(psifm, psifmd, eigv, psierr, PE)
  !Error Estimation
  eigerr = (eigv - peigv)/eigv
  ResErr = ResidualError(phifm, psifm, eigv, 1, ng, PE)
  !Log Out
  IF(iter .eq. 1) ResErr0 = ResErr
  IF(CMFDMaster) WRITE(mesg,'(A9, I9, F22.6, 3x, F10.5, 1pe15.3)') 'MGOUTER', ItrCntl%Cmfdit, eigv, ResErr/ResErr0, ResErr
  IF(CMFDMaster) CALL message(io8, FALSE, TRUE, mesg)
  !Convergence Check
  lExit = .TRUE.; lConv = .FALSE.
  IF((ResErr/ResErr0) .lt. convcrit ) lconv = TRUE
  lExit = lExit .AND. lConv
  IF(iter .LE. nitermin) lExit = .FALSE.
  IF(iter .GE. nitermin .AND. ResErr .LT. 1.0E-7_8) lExit = .TRUE.
  IF(lExit) EXIT
  IF(nitermax .EQ. iter) EXIT
  !IF(iter .gt. nitermin .and. lconv .AND. nAxUpdt .GE. 2) EXIT
!  IF(mod(iter, 5) .EQ. 0 .and. nTracerCntl%lGcCmfd) THEN
!    CALL GcCmfdAcc(Core, CmInfo, Eigv, GroupInfo, GcGroupInfo, nTracerCntl, ItrCntl, PE)
!    CALL CP_VA(PsifmD(1 : nxy, myzbf : myzef), Psifm(1 : nxy, myzbf : myzef), nxy, myzef - myzbf + 1)
!    CALL CmfdPsiUpdt(phifm, psifm)
!  ENDIF
  !CAlling 1D-Axial Kernel
  !IF(mod(iter,5) .eq.0 .and. (nTracerCntl%lBoronSearch .or. nTracerCntl%lBoronSearch)) THEN
  !  CALL UpdtBoronPPM(nTracerCntl%target_eigv, eigv, nTracerCntl%BoronPPM, .FALSE., PE%MASTER)
  !  CALL UpdtBoronCmfdXs(Core, Fxr, Phis, PinXS,  nTracerCntl%BoronPPM, myzb, myze, ng)
  !  IF(.not. nTracerCntl%l3dim .AND. .NOT. mod(iter, 5) .eq. 0) CALL SetCmfdLinearSystem(TRUE, TRUE, AxSolver)
  !ENDIF
  IF(mod(iter, itrcntl%CMFDItrCntl%ncmfdpnodal) .eq. 0 .and. nTracerCntl%l3dim) THEN
    nAxUpdt = nAxUpdt  + 1
    IF(lAxDhatConv) EXIT
    IF(CMFDMaster) WRITE(mesg, '(2x, A32, 2x,A10)') 'Performing Axial Calculation : ',AxSolverName(AxSolver)
    IF(CMFDMaster) CALL message(io8, FALSE, TRUE, mesg)
    !Axail Dhat Update
!    IF(PE%nProc .GT. 1) THEN
!      CALL MpiAxialSolver(Core, CMInfo, GroupInfo, Eigv, ng, AxSolver, .TRUE., .FALSE., PE)
!    ELSE
      CALL MpiAxialSolver(Core, CMInfo, GroupInfo, Eigv, ng, AxSolver, lreset, .FALSE., PE, nTracerCntl%lnonFuelpinFDM)
!    ENDIF
   ! lAxDhatConv = AxDhatConvCheck(ItrCntl%CMFDItrCntl%AxDhatConvCrit)
    !LINEAR SYSTEM
    !CALL SetCmfdLinearSystem(TRUE,  nTracerCntl%l3dim, AxSolver)
    !CALL SetPcondiLinearSystem(PhiFm, TRUE,  nTracerCntl%l3dim, AxSolver, PE)

    CALL SetCmfdLinearSystem(TRUE, TRUE, AxSolver)
    CALL MakeBiLU(Core, CmInfo%CoreCMFDLS(1:ng), PhiFm, ng, PE)
  ENDIF
  IF(mod(iter, 5) .EQ. 0 .and. nTracerCntl%lGcCmfd) THEN
    CALL GcCmfdAcc(Core, CmInfo, Eigv, TRUE, GroupInfo, GcGroupInfo, nTracerCntl, ItrCntl, PE)
    CALL CP_VA(PsifmD(1 : nxy, myzbf : myzef), Psifm(1 : nxy, myzbf : myzef), nxy, myzef - myzbf + 1)
    CALL CmfdPsiUpdt(phifm, psifm)
  ENDIF
END DO

CALL ConvertSubPlnPhi(PhiC, PhiFm, 2)

!Update Radial MOC Solution
!CALL MOCSolUpdt(Core, PinXS, PHIC, phis, myzb, myze, ng)
!CALL MOCPhiInUpdt(Core, PinXs, PhiC, CmInfo%RayInfo4Cmfd, myzb, myze, ng)
CALL MOCSolUpdt(Core, FmInfo, CmInfo, myzb, myze, ng)
CALL MOCPhiInUpdt(Core, CmInfo, myzb, myze, ng)
!Axial Source For MOC
IF(nTracerCntl%l3dim) THEN
  CALL AxSrcUpdate(Core, CmInfo, myzb, myze, 1, ng, PE, AxSolver)
  IF(.NOT. nTracerCntl%lBenchXs)THEN ! benchmark XS bug fixed 160314
    CALL VoidLkgCorrection(Core, FmInfo, CmInfo, myzb, myze, 1, ng, PE) !150416 BYS edit off
  ENDIF
ENDIF


IF(nTracerCntl%lLinSrc) CALL MOCLinSrcUpdt(Core, PinXS, PHIC, FmInfo%LinSrcSlope, myzb, myze, ng)
!CMFD Time Check
!IF(nTracerCntl%lGcCmfd) THEN
!  CALL GcCmfdAcc(Core, CmInfo, Eigv, GroupInfo, GcGroupInfo, nTracerCntl, ItrCntl, PE)
!ENDIF


CmfdTImeEnd = nTracer_dclock(FALSE, FALSE)
TimeChk%CmfdTime = TimeChk%CmfdTime + (CmfdTimeEnd - CmfdTimeBeg)
!Free
NULLIFY(Pin);
NULLIFY(Fxr); NULLIFY(phis)
NULLIFY(PinXS);

NULLIFY(PhiC); NULLIFY(PhiFM)
NULLIFY(PsiC); NULLIFY(PsiCD)
NULLIFY(PsiFm); NULLIFY(PsiFmD)
NULLIFY(MmtFm); NULLIFY(Mmt)
NULLIFY(RadJout)

END SUBROUTINE

SUBROUTINE DcplPlnCmfd(Core, CMInfo, FMInfo, eigv, ng, lcmfd, lreset, GroupInfo, PE, nTracerCntl, ItrCntl)
USE PARAM
USE TYPEDEF,     ONLY : PE_TYPE,      CoreInfo_Type,      Pin_Type,         &
                        CMInfo_Type,  FmInfo_Type,                          &
                        FxrInfo_Type, CmfdLs_TYPE,        PinXs_Type,       &
                        AxFlx_Type,   GroupInfo_Type
USE CMFD_mod,    ONLY : CmfdPinXS,    CmfdLS,                              &
                        PhiC1g,       SRC,                                 &
                        SetCMFDEnviorment, HomoXsGen,       RadCouplingCoeffGen,       &
                        UpdatePhiC,        CmfdSrcUpdt,     CmfdPsiUpdt,               &
                        CmfdEigUpdate,     ResidualError,   MOCSolUpdt,                &
                        MOCPhiInUpdt,      AddConstCmfdSrc, AddCmfdPxs
!USE CORE_MOD,    ONLY : GroupInfo
USE BiCGSTAB_mod, ONLY : BiCGSTAB
USE SUBGRP_MOD,  ONLY : FxrChiGen
USE BasicOperation, ONLY : CP_CA, CP_VA
USE IOUTIL,         ONLY : message
USE FILES,       ONLY : io8
USE CNTL,        ONLY : nTracerCntl_Type
USE itrcntl_mod, ONLY : ItrCntl_TYPE, CMFDItrCntl_TYPE
use timer, only : nTracer_dclock, TimeChk
#ifdef MPI_ENV
USE MPIDcpl_Mod, ONLY : MPI_DcplMessage, PrintDcplCmfdLog
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CMInfo_Type) :: CMInfo
TYPE(FMInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_TYPE) :: PE
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(ItrCntl_TYPE) :: ItrCntl
REAL :: eigv
REAL :: CmfdTimeBeg, CmfdTimeEnd
INTEGER :: ng
LOGICAL :: lreset, lcmfd
LOGICAL, SAVE :: lfirst

!Pointing Varibables
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Fxrinfo_type), POINTER :: Fxr(:,:)
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
!TYPE(AxFlx_Type), POINTER :: AxFlx(:, :)
REAL, POINTER :: PHIS(:,:,:)
REAL, POINTER :: PhiC(:, :, :), PSIC(:,:), PSICD(:,:)
REAL, POINTER :: PhiFm(:, :, :), PsiFm(:, :), PsiFmD(:, :)
REAL, POINTER :: RadJout(:, :, :, :, :)
!
INTEGER :: ig, l, ITER, jsweep
INTEGER :: GrpBeg, GrpEnd, nGroupInfo
INTEGER :: myzb, myze, myzbf, myzef, nxy

REAL :: psierr, eigerr, peigv, ResErr, ResErr0
LOGICAL :: lXsLib, lscat1, lsubplane
!Iteration COntrol
INTEGER :: nitermax, nitermin
REAL :: convcrit
LOGICAL :: lconv, lEigProb
LOGICAL :: LParallel

DATA lfirst /TRUE/

!CMFD Time Check
CmfdTImeBeg = nTracer_dclock(FALSE, FALSE)

!Pointing
Pin => Core%Pin;
Fxr => FmInfo%Fxr; phis => FmInfo%phis
PinXS => CMInfo%PinXS

PhiC => CMInfo%PhiC; PhiFM => CMInfo%PhiFM
PsiC => CMInfo%PsiC; PsiCD => CMInfo%PsiCD
PsiFm => CMInfo%PsiFm; PsiFmD => CMInfo%PsiFmD
RadJout => CMinfo%RadJout !ninout/ 1:in 2:out 3:surfphi ! BYS edit 16//02/11

nitermax = itrcntl%CMFDItrCntl%nitermax
nitermin = itrcntl%CMFDItrCntl%nitermin
convcrit = itrcntl%CMFDItrCntl%convcrit

nxy = Core%nxy; myzb = PE%myzb; myze = PE%myze
myzbf = PE%myzbf; myzef = PE%myzef
lXsLib = nTracerCntl%lXsLib; lscat1 = nTracerCntl%lScat1
lsubplane = nTracerCntl%lSubPlane
!Determine Eig problem or not
lEigProb = TRUE

LParallel = PE%lDcplParallel

IF(.NOT. Core%lFuelPlane(myzb)) lEigProb = FALSE

IF(lfirst .OR. lReset) THEN
  CALL SetCMFDEnviorment(Core, CmInfo, ng, PE,  TRUE)
  CALL SetCmfdMpiOmpEnv(Core, CMInfo%CoreCmfdLs , ng, PE)
  lfirst = FALSE
ENDIF

!Homogenization ProCedure\
WRITE(mesg,'(a)') 'Cell Homogenization (H)...'
IF(.NOT. lParallel) CALL message(io8, TRUE, TRUE, mesg)
#ifdef MPI_ENV
IF(lParallel) CALL MPI_DcplMessage(mesg, PE%myrank, TRUE, FALSE, FALSE)
#endif
CALL FxrChiGen(Core, Fxr, FmInfo, GroupInfo, nTracerCntl, myzb, myze)
CALL HomoXsGen(Core, Fxr, Phis, PinXs, myzb, myze, ng, lXsLib, lScat1, FALSE)
!
!CALL RadCouplingCoeffGen(Pin, RadJout)
CALL RadCouplingCoeffGen(Core, PinXs, RadJout, ng, .TRUE., PE)
CALL AddCmfdPxs(FmInfo%AxPxs, myzb, myze, 1, ng)

WRITE(mesg, '(a)') 'Performing CMFD Calculation...'
IF(.NOT. lParallel) CALL message(io8, TRUE, TRUE, mesg)
#ifdef MPI_ENV
IF(lParallel) THEN
  CALL MPI_DcplMessage(mesg, PE%myrank, TRUE, FALSE, FALSE)
  IF(PE%myrank .eq. 0) CALL message(io8, TRUE, TRUE, mesg)
ENDIF
#endif
!Update Homo PHI
CALL UpdatePhiC(PinXs, PhiC)

CALL SetCmfdLinearSystem(TRUE,  nTracerCntl%l3dim, lP1SENM)
CALL CmfdPsiUpdt(PhiFm, PsiFm)

lconv = FALSE;
nGroupInfo = 2
IF(.NOT. GroupInfo%lUpScat) nGroupInfo = 1
ResErr0 =  ResidualError(phifm, psifm, eigv, 1, ng, PE)
DO iter = 1, nitermax
  ItrCntl%Cmfdit = ItrCntl%Cmfdit + 1
  DO jsweep =1, nGroupInfo
    GrpBeg = 1; GrpEnd = ng
    IF(jsweep .GT. 1) THEN
      GrpBeg = GroupInfo%UpScatRange(1); GrpEnd = GroupInfo%UpScatRange(2)
    ENDIF
    DO ig = GrpBeg, GrpEnd
      !Source Update
      CALL CmfdSrcUpdt(Src, Psifm, phifm, eigv, ig)
      IF(.NOT. lEigProb) CALL AddConstCmfdSrc(Src, nTracerCntl%ConstSrc)
      !Solving Linear System
      CALL CP_VA(Phic1g(1:nxy, myzbf : myzef), Phifm(1:nxy, myzbf : myzef, ig), nxy, myzef - myzbf + 1)
      CALL BiCGSTAB(CmfdLs(ig), Phic1g, SRC, itrcntl%InSolverItrCntl)
      CALL CP_VA(Phifm(1:nxy, myzbf : myzef, ig), Phic1g(1:nxy, myzbf : myzef), nxy, myzef - myzbf + 1)
    ENDDO
  END DO
  CALL CP_VA(PsifmD(1 : nxy, myzbf : myzef), Psifm(1 : nxy, myzbf : myzef), nxy, myzef - myzbf + 1)
  !Fission Source Update
  CALL CmfdPsiUpdt(phifm, psifm)
  !EigenValue Update
  peigv = eigv

  IF(lEigProb) THEN
    CALL CmfdEigUpdate(psifm, psifmd, eigv, psierr)
    !Error Estimation
    eigerr = (eigv - peigv)/eigv
    ResErr =  ResidualError(phifm, psifm, eigv, 1, ng, PE)
  ELSE
    ResErr =  ResidualError(phifm, psifm, eigv, 1, ng, PE, nTracerCntl%ConstSrc)
  ENDIF


  !Log Out
  IF(iter .eq. 1) ResErr0 = ResErr
  WRITE(mesg,'(A9, I9, F20.6, 3x, F10.5, 1pe15.3)') 'MGOUTER', ItrCntl%Cmfdit, eigv, ResErr/ResErr0, ResErr
  IF(.NOT. lParallel) CALL message(io8, FALSE, TRUE, mesg)
#ifdef MPI_ENV
  IF(lParallel) CALL MPI_DcplMessage(mesg, PE%myrank, FALSE, FALSE, FALSE)
#endif
  !Convergence Check
  IF((ResErr/ResErr0) .lt. convcrit ) lconv = TRUE
  IF(iter .gt. nitermin .and. lconv) EXIT
  !CAlling 1D-Axial Kernel
ENDDO
!
#ifdef MPI_ENV
  IF(lParallel) THEN
    CALL PrintDcplCmfdLog(ItrCntl%CmfdIt, myzb, eigv,  ResErr/ResErr0, ResErr, &
                          PE%myrank, PE%nProc, PE%MPI_DCPLMASTER_COMM)
  ENDIF
#endif
!Update Radial MOC Solution
CALL MOCSolUpdt(Core, FmInfo, CmInfo, myzb, myze, ng)
!CALL MOCSolUpdt(Core, PinXS, PHIC, phis, myzb, myze, ng)

!CALL MOCPhiInUpdt(Core, PinXs, PhiC, CmInfo%RayInfo4Cmfd, myzb, myze, ng)
CALL MOCPhiInUpdt(Core, CmInfo, myzb, myze, ng)

!CMFD Time Check
CmfdTImeEnd = nTracer_dclock(FALSE, FALSE)
TimeChk%CmfdTime = TimeChk%CmfdTime + (CmfdTimeEnd - CmfdTimeBeg)

!Free
NULLIFY(Pin);
NULLIFY(Fxr); NULLIFY(phis)
NULLIFY(PinXS)

NULLIFY(PhiC); NULLIFY(PhiFM)
NULLIFY(PsiC); NULLIFY(PsiCD)
NULLIFY(PsiFm); NULLIFY(PsiFmD)
NULLIFY(RadJout)

END SUBROUTINE
