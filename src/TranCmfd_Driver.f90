#include <defines.h>
SUBROUTINE TranCmfd_Driver(Core, FmInfo, CmInfo, TranInfo, THInfo, GroupInfo, TranCntl, nTracerCntl, ItrCntl, PE)
USE PARAM
USE TYPEDEF,           ONLY : CoreInfo_Type,           FmInfo_Type,           CmInfo_Type,           &
                              TranInfo_Type,           GroupInfo_Type,        TranCntl_Type,         &
                              ThInfo_Type,             PE_Type,                                      &
                              FxrInfo_Type,            CmfdLs_TYPE,           PinXs_Type,            &
                              Pin_Type,                AxFlx_Type
USE CNTL,              ONLY : nTracerCntl_Type
USE itrcntl_mod,       ONLY : ItrCntl_TYPE,            CMFDItrCntl_TYPE
USE CMFD_mod,          ONLY : CmfdPinXS,               CmfdLS,                                       &
                              PhiC1g,                  SRC,                                          &
                              SetCMFDEnviorment,       HomoXsGen_Cusping,     RadCouplingCoeffGen,   &
                              UpdatePhiC,              ConvertSubPlnPhi,                             &
                              CmfdSrcUpdt,             CmfdPsiUpdt,                                  &
                              CmfdEigUpdate,           ResidualError,         MOCSolUpdt,            &
                              MOCPhiInUpdt,            MOCLinSrcUpdt,         HomBuckling
USE CORE_MOD,          ONLY : GcGroupInfo 
USE TRANCMFD_MOD,      ONLY : TrSrc,                   PrecSrc,                                      &
                              SetTranCmfdEnv,          HomKineticParamGen,    SetCmfdPrecParam,      &
                              SetTranCmfdLinearSystem, CmfdPrecSrcUpdt,       CmfdTranSrc,           &
                              TranResidualError,       CmfdPrecUpdt,          CmfdSteadySrcUpdt,     &
                              CmfdSteadyPsiUpdt,       TranResidualError_rev
USE TranGcCmfd_mod,    ONLY : TranGcCmfdAcc
USE BiCGSTAB_mod,      ONLY : BiCGSTAB
USE MPIAxSolver_Mod,   ONLY : AxNCntl_Type,                                                          &
                              AxSrcUpdate,             MpiTranAxialSolver
USE BiLU_MOD,          ONLY : MakeBiLU
USE SUBGRP_MOD,        ONLY : FxrChiGen  
USE BasicOperation,    ONLY : CP_CA,                   CP_VA,                 AD_VA,                 &
                              MULTI_CA,                                                              &
                              CP_CA_OMP,               CP_VA_OMP,             AD_VA_OMP,             &
                              MULTI_CA_OMP
USE IOUTIL,            ONLY : message
USE FILES,             ONLY : io8
USE timer,          only : nTracer_dclock, TimeChk
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(ItrCntl_TYPE) :: ItrCntl
TYPE(PE_Type) :: PE
  
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Fxrinfo_type), POINTER :: Fxr(:,:)
TYPE(AxNCntl_Type) :: AxNCntl
REAL, POINTER :: PHIS(:,:,:)     
REAL, POINTER :: PhiC(:, :, :), PSIC(:,:), PSICD(:, :)
REAL, POINTER :: TranPsiFm(: ,:), TranPsiFmd(:, :)
REAL, POINTER :: TranPhiFm(:, : ,:), TranPhiCm(:, :, :)
REAL, POINTER :: PrecFm(:, :, :), PrecCm(:, :, :), ResSrcFm(:, :, :)
REAL, POINTER :: PhiFm(:, :, :), PsiFm(:, :), PsiFmD(:, :)
REAL, POINTER :: RadJout(:, :, :, :, :)

REAL :: eigv0
REAL :: ResErr, ResErr0, Convcrit
REAL :: CmfdTImeBeg, CmfdTImeEnd

!TYPE(AxNCntl_Type) :: AxNCntl

INTEGER :: ng, nprec, ngroupinfo, GrpBeg, GrpEnd
INTEGER :: myzb, myze, myzbf, myzef
INTEGER :: nxy
INTEGER :: AxSolver
INTEGER :: nitermax, nitermin
INTEGER :: ig, iter, jsweep

LOGICAL :: lXsLib, lSubPlane, lScat1, lExit, lconv
LOGICAL :: lRadDhatUpdt
LOGICAL :: CmfdMaster, CmfdSlave
LOGICAL, SAVE :: lFirst
DATA lfirst /TRUE/


CmfdTImeBeg = nTracer_dclock(FALSE, FALSE)

CMFDMaster = PE%CMFDMaster; CMFDSlave = PE%CMFDSlave



IF(lfirst) THEN
  CALL SetCMFDEnviorment(Core, CmInfo, GroupInfo%ng, PE)
  CALL SetCmfdMpiOmpEnv(Core, CMInfo%CoreCmfdLs , GroupInfo%ng, PE)
  CALL SetTranCmfdEnv(Core, TranInfo, PE)
  lfirst = FALSE
ENDIF

myzb = PE%myzb; myze = PE%myze
myzbf = PE%myzbf; myzef = PE%myzef
ng = GroupInfo%ng
nxy = Core%nxy
lXsLib = nTracerCntl%lXsLib; lscat1 = nTracerCntl%lScat1
lsubplane = nTracerCntl%lSubPlane
AxSolver = nTracerCntl%AxSolver


nitermax = itrcntl%CMFDItrCntl%nitermax
nitermin = itrcntl%CMFDItrCntl%nitermin
convcrit = itrcntl%CMFDItrCntl%convcrit

eigv0 = TranInfo%eigv0

Fxr => FmInfo%Fxr
Phis => FmInfo%Phis
PhiC => CMInfo%PhiC; PhiFM => CMInfo%PhiFM
PsiC => CMInfo%PsiC; PsiCD => CMInfo%PsiCD 
PsiFm => CMInfo%PsiFm; PsiFmD => CMInfo%PsiFmD
RadJout => CMinfo%RadJout

TranPsiFm => CmInfo%TranPsiFm; TranPsiFmd => CmInfo%TranPsiFmd
TranPhiFm => CmInfo%TranPhiFm; TranPhiCm => CMInfo%TranPhiCm
PrecCm => CmInfo%PrecCm
PrecFm => CmInfo%PrecFm
ResSrcFm => CmInfo%ResSrcFm

WRITE(mesg,'(a)') 'Cell Homogenization (H)...'
IF(CMFDMaster) CALL message(io8, TRUE, TRUE, mesg)    

CALL FxrChiGen(Core, Fxr, FmInfo, GroupInfo, nTracerCntl, myzb, myze)
CALL HomoXsGen_Cusping(Core, FmInfo, Fxr, Phis, CmfdPinXS, myzb, myze, ng, lXsLib, lScat1, FALSE)

lRadDhatUpdt = .TRUE.
IF(ItrCntl%Cmfdit0 .EQ. ItrCntl%Cmfdit) lRadDhatUpdt = .FALSE.
IF(.NOT. TranCntl%lMocUpdt) lRadDhatUpdt = .FALSE.
  
CALL RadCouplingCoeffGen(Core, CmfdPinXS, RadJout, ng, lRadDhatUpdt, PE)

CALL UpdatePhiC(CmfdPinXS, PhiC)
!Calculate homgenized chid and beta
CALL HomKineticParamGen(Core, FmInfo, CmInfo, TranInfo, GroupInfo, nTracerCntl, PE)
!Calculate the Precursor Quadratic approx. Coefficient
CALL SetCmfdPrecParam(Core, CmInfo, TranInfo, TranCntl, nTracerCntl, PE)
CALL CmfdPrecSrcUpdt(PrecSrc, PrecFm, TranPsiFm, TranPsiFmd, TranCntl)

!Solve Axial 1-D
IF(nTracerCntl%l3dim .AND. TranCntl%lAxNUpdt) THEN
  IF(CMFDMaster) WRITE(mesg, '(2x, A32, 2x,A10)') 'Performing Axial Calculation : ',AxSolverName(AxSolver)
  IF(CMFDMaster) CALL message(io8, FALSE, TRUE, mesg)
  AxNCntl%lDtilUpdt = .FALSE.;
  IF(ItrCntl%Cmfdit0 .EQ. ItrCntl%Cmfdit) AxNCntl%lDtilUpdt = .TRUE.
  AxNCntl%AxSolverMod = AxSolver;AxNCntl%lReset = .FALSE.
  AxNCntl%ng0 = ng; AxNCntl%lTransient = .TRUE.
  AxNCntl%lXsComm = .TRUE.; AxNCntl%eigv = eigv0
!  
  CALL MpiTranAxialSolver(Core, CmInfo, TranInfo, GroupInfo, TranCntl, AxNCntl, PE)
ENDIF 

!LINEAR SYSTEM
CALL SetTranCmfdLinearSystem(TRUE,  nTracerCntl%l3dim, AxSolver, TranCntl)
CALL MakeBiLU(Core, CmInfo%CoreCMFDLS(1:ng), PhiFm, ng, PE)

lconv = FALSE
nGroupInfo = 2;
IF(.NOT. GroupInfo%lUpScat) nGroupInfo = 1

DO iter = 1, nitermax
  ItrCntl%Cmfdit = ItrCntl%Cmfdit + 1
  DO jsweep =1, nGroupInfo
    GrpBeg = 1; GrpEnd = ng
    IF(jsweep .GT. 1) THEN
      GrpBeg = GroupInfo%UpScatRange(1); GrpEnd = GroupInfo%UpScatRange(2)
    ENDIF
    DO ig = GrpBeg, GrpEnd
      CALL CmfdSteadySrcUpdt(Src, Psifm, phifm, 1._8, ig, PE%nCmfdThread)
      CALL CmfdTranSrc(TrSrc, PhiFm, TranPhiFm, PsiFm, PrecSrc, ResSrcFm, TranCntl, ig, PE%nCmfdThread)
      CALL AD_VA_OMP(TrSrc(1:nxy, myzbf:myzef), TrSrc(1:nxy, myzbf:myzef), Src(1:nxy, myzbf:myzef), nxy, myzef-myzbf+1, PE%nCmfdThread)
      CALL CP_VA_OMP(Phic1g(1:nxy, myzbf : myzef), Phifm(1:nxy, myzbf : myzef, ig), nxy, myzef - myzbf + 1, PE%nCmfdThread)    
      CALL BiCGSTAB(CmfdLs(ig), Phic1g, TrSRC, itrcntl%InSolverItrCntl, itrcntl%innerit)
#ifndef MPI_ENV
      CALL CP_VA_OMP(Phifm(1:nxy, myzbf : myzef, ig), Phic1g(1:nxy, myzbf : myzef), nxy, myzef - myzbf + 1, PE%nCmfdThread)    
#else
      CALL CP_VA_OMP(Phifm(1:nxy, myzbf-1 : myzef+1, ig), Phic1g(1:nxy, myzbf-1 : myzef+1), nxy, myzef - myzbf + 3, PE%nCmfdThread)    
#endif
    ENDDO
  ENDDO
  
  CALL CP_VA_OMP(PsifmD(1 : nxy, myzbf : myzef), Psifm(1 : nxy, myzbf : myzef), nxy, myzef - myzbf + 1, PE%nCmfdThread)    
  !Fission Source Update
  CALL CmfdSteadyPsiUpdt(phifm, psifm, PE%nCmfdThread)
  CALL MULTI_CA_OMP(1._8 / eigv0, PsiFm(1:nxy, myzbf:myzef), nxy, myzef - myzbf + 1, PE%nCmfdThread)
  ResErr = TranResidualError_rev(phifm, psifm, TranPhiFm, PrecSrc, ResSrcFm, TranCntl, PE)
  IF(iter .eq. 1) ResErr0 = ResErr
  IF(CMFDMaster) WRITE(mesg,'(A9, I9, F22.6, 3x, F10.5, 1pe15.3)') 'MGOUTER', ItrCntl%Cmfdit, eigv0, ResErr/ResErr0, ResErr
  IF(CMFDMaster) CALL message(io8, FALSE, TRUE, mesg)
  !Convergence Check
  lExit = .TRUE.; lConv = .FALSE.
  IF((ResErr/ResErr0) .lt. convcrit ) lconv = TRUE
  lExit = lExit .AND. lConv
  IF(iter .LE. nitermin) lExit = .FALSE.
  IF(iter .GE. nitermin .AND. ResErr .LT. 1.0E-7_8) lExit = .TRUE.
  IF(lExit) EXIT 
  IF(nitermax .EQ. iter) EXIT
   
  IF(mod(iter, 5) .eq. 0 .and. nTracerCntl%l3dim .AND. TranCntl%lAxNUpdt) THEN
    !IF(CMFDMaster) WRITE(mesg, '(2x, A32, 2x,A10)') 'Performing Axial Calculation : ',AxSolverName(AxSolver)
    !IF(CMFDMaster) CALL message(io8, FALSE, TRUE, mesg) 
    !AxNCntl%AxSolverMod = AxSolver;AxNCntl%lReset = .FALSE.
    !AxNCntl%ng0 = ng; AxNCntl%lTransient = .TRUE.
    !AxNCntl%lXsComm = .TRUE.; AxNCntl%eigv = eigv0
!  
    !CALL MpiTranAxialSolver(Core, CmInfo, TranInfo, GroupInfo, TranCntl, AxNCntl, PE)    
    !CALL SetTranCmfdLinearSystem(TRUE,  nTracerCntl%l3dim, AxSolver, TranCntl)
    !CALL MakeBiLU(Core, CmInfo%CoreCMFDLS(1:ng), PhiFm, ng, PE)
  ENDIF  
  IF(mod(iter, 5) .EQ. 0 .and. nTracerCntl%lGcCmfd) THEN
    CALL TranGcCmfdAcc(Core, CmInfo, TranInfo, Eigv0, .TRUE., GroupInfo, GcGroupInfo, TranCntl, nTracerCntl, ItrCntl, PE)
    CALL CP_VA(PsifmD(1 : nxy, myzbf : myzef), Psifm(1 : nxy, myzbf : myzef), nxy, myzef - myzbf + 1)    
    CALL CmfdSteadyPsiUpdt(phifm, psifm, PE%nCmfdThread)
    CALL MULTI_CA_OMP(1._8 / eigv0, PsiFm(1:nxy, myzbf:myzef), nxy, myzef - myzbf + 1, PE%nCmfdThread)   
  ENDIF
ENDDO


!CALL ConvertSubPlnPhi(PhiC, PhiFm, 2)
!ItrCntl%lconv  = .FALSE.
!IF(ResErr .LT. 1.E-4) ItrCntl%lconv = .TRUE.
!Axial Source For MOC
IF(nTracerCntl%l3dim) CALL AxSrcUpdate(Core, CmInfo, myzb, myze, 1, ng, PE, AxSolver)

!Update Radial MOC Solution
CALL MOCSolUpdt(Core, FmInfo, CmInfo, myzb, myze, ng)
CALL MOCPhiInUpdt(Core, CmInfo, myzb, myze, ng)

CmfdTImeEnd = nTracer_dclock(FALSE, FALSE)
TimeChk%CmfdTime = TimeChk%CmfdTime + (CmfdTimeEnd - CmfdTimeBeg)

CmfdTImeEnd = nTracer_dclock(FALSE, FALSE)
TimeChk%CmfdTime = TimeChk%CmfdTime + (CmfdTimeEnd - CmfdTimeBeg)

END SUBROUTINE
