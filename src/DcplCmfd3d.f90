#include <defines.h>
SUBROUTINE DcplCmfd3d(Core, CmInfo, FmInfo, DcplInfo, Eigv, ng, GroupInfo, GcGroupInfo, nTracerCntl, ItrCntl, PE)
USE PARAM
USE TYPEDEF,        ONLY : PE_TYPE         ,CoreInfo_Type      ,Pin_Type          &
                          ,CMInfo_Type     ,FmInfo_Type        ,DcplInfo_Type     &
                          ,PE_Type         ,GroupInfo_Type                        &
                          ,FxrInfo_Type    ,CmfdLs_TYPE        ,PinXs_Type        &
                          ,AxFlx_Type      ,Pin_Type
USE CNTL,           ONLY : nTracerCntl_Type
USE ItrCntl_Mod,    ONLY : ItrCntl_TYPE

USE CMFD_mod,       ONLY : CmfdPinXS        ,CmfdLS                                &
                          ,PhiC1g           ,SRC                                   &
                          ,SetCMFDEnviorment                                       &
                          ,ConvertSubPlnPhi ,CmfdSrcUpdt        ,CmfdPsiUpdt       &
                          ,CmfdEigUpdate    ,ResidualError      
USE DcplCmfd_mod,   ONLY : DcplHomXsGen     ,DcplRadCouplingCoeffGen               &
                          ,DcplAxSrcUpdt
!USE AxSolver_mod,   ONLY : AxialSolver      ,AxSrcUpdate
USE MpiAxSolver_mod,ONLY : MpiAxialSolver,    AxSrcUpdate
USE BiCGSTAB_mod,   ONLY : BiCGSTAB
USE TIMER,          ONLY : nTracer_dclock   ,TimeChk
USE IOUTIL,         ONLY : message
USE FILES,          ONLY : io8
USE BasicOperation, ONLY : CP_CA, CP_VA
#ifdef MPI_ENV
USE MpiAxSolver_mod,ONLY : MpiAxialSolver
USE MpiComm_mod,    ONLY : MPI_SYNC
#endif 
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(DcplInfo_Type) :: DcplInfo
TYPE(GroupInfo_Type) :: GroupInfo, GcGroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(ItrCntl_TYPE) :: ItrCntl
TYPE(PE_TYPE) :: PE
REAL :: EigV
INTEGER :: ng

TYPE(PinXs_Type), POINTER :: PinXs(:, :)
TYPE(Pin_Type), POINTER :: Pin(:)

REAL, POINTER :: PhiC(:, :, :), PSIC(:,:), PSICD(:,:)
REAL, POINTER :: PhiFm(:, :, :), PsiFm(:, :), PsiFmD(:, :)

REAL :: CmfdTimeBeg, CmfdTimeEnd

!
INTEGER :: nitermax, nitermin, AxSolver, nGroupInfo
REAL :: ConvCrit, EigConvCrit
REAL :: psierr, eigerr, peigv, ResErr, ResErr0, ResReduction
LOGICAL :: lconv, lAxDhatConv
!
INTEGER :: jsweep, iter, ig
!
INTEGER :: myzb, myze, myzef, myzbf, nxy
INTEGER :: GrpBeg, GrpEnd
LOGICAL :: lXslib, lsubplane, lFeedBack, lXsGen, lDtilOnly
LOGICAL :: CmfdMaster, CmfdSlave

CMFDMaster = PE%CMFDMaster; CMFDSlave = PE%CMFDSlave
CmfdTimeBeg = nTracer_dclock(FALSE, FALSE)

Pin => Core%Pin; 
PinXS => CMInfo%PinXS

PhiC => CMInfo%PhiC; PhiFM => CMInfo%PhiFM
PsiC => CMInfo%PsiC; PsiCD => CMInfo%PsiCD 
PsiFm => CMInfo%PsiFm; PsiFmD => CMInfo%PsiFmD

nxy = Core%nxy; myzb = PE%myzb; myze = PE%myze
myzbf = PE%myzbf; myzef = PE%myzef

lXsLib = nTracerCntl%lXsLib
lsubplane = nTracerCntl%lSubPlane
AxSolver = nTracerCntl%AxSolver
lFeedBack = nTracerCntl%lFeedBack
lXsGen = nTracerCntl%lDcplCmfdXsGen

!IF(nTracerCntl%AxSolver .EQ. 2) AxSolver = 1  
!CMFD Enviorment Variable
CALL SetCMFDEnviorment(Core, CmInfo, ng, PE)

! 
!ITERATION Control Variables
nitermax = itrcntl%DcplCMFD3dItrCntl%nitermax
nitermin = itrcntl%DcplCMFD3dItrCntl%nitermin
convcrit = itrcntl%DcplCMFD3dItrCntl%Res3dCmfdConv
ResReduction = itrcntl%DcplCMFD3dItrCntl%convcrit
EigConvCrit= itrcntl%DcplCMFD3dItrCntl%eigv3dCmfdconv
!XS Generation
WRITE(mesg,'(a)') 'Cell Homogenization (H)...'
IF(CMFDMaster)  CALL message(io8, TRUE, TRUE, mesg)      
CALL DcplHomXsGen(Core, DcplInfo, FmInfo%Fxr, PinXs, myzb, myze, ng, lXsGen,      &
                  lXslib, nTracerCntl%lXsFtn, nTracerCntl%XsFtnMod,GroupInfo)
CALL DcplRadCouplingCoeffGen(Pin, DcplInfo, .NOT. lXsGen .AND. nTracerCntl%lXsFtn, lXsGen)

WRITE(mesg, '(a)') 'Performing CMFD Calculation...'
IF(CMFDMaster)  CALL message(io8, TRUE, TRUE, mesg)      

!Update Homo PHI
!CALL UpdatePhiC(PinXs, PhiC)
CALL ConvertSubPlnPhi(PhiC,PhiFm, 1)

#ifndef MPI_ENV
IF(nTracerCntl%l3dim) THEN
  WRITE(mesg, '(2x, A32, 2x,A10)') 'Performing Axial Calculation : ',AxSolverName(AxSolver)
  IF(CmfdMaster) CALL message(io8, FALSE, TRUE, mesg)
  lDtilOnly = .FALSE.;
  IF(ItrCntl%Cmfdit .NE. 0) CALL AxialSolver(Core, CMInfo, Eigv, ng, AxSolver, TRUE, lDtilOnly,PE)
  ItrCntl%Axit = ItrCntl%Axit + 1
ENDIF
#else
IF(nTracerCntl%l3dim) THEN
  IF(CMFDMaster) WRITE(mesg, '(2x, A32, 2x,A10)') 'Performing Axial Calculation : ',AxSolverName(AxSolver)
  IF(CMFDMaster) CALL message(io8, FALSE, TRUE, mesg)
  lDtilOnly = .FALSE.;
  IF(ItrCntl%Cmfdit .EQ. 0) lDtilOnly = .TRUE.
!  IF(PE%nProc .GT. 1) THEN
!    CALL MpiAxialSolver(Core, CMInfo, GroupInfo, Eigv, ng, AxSolver, .TRUE., lDtilOnly, PE)
!  ELSE
    CALL MpiAxialSolver(Core, CMInfo, GroupInfo, Eigv, ng, AxSolver, .TRUE., lDtilOnly, PE, nTracerCntl%lnonFuelpinFDM) 
!  ENDIF
ENDIF 
#endif
CALL SetCmfdLinearSystem(TRUE,  nTracerCntl%l3dim, AxSolver)
CALL CmfdPsiUpdt(PhiFm, PsiFm)

lconv = FALSE; lAxDhatConv = FALSE
nGroupInfo = 2
IF(.NOT. GroupInfo%lUpScat) nGroupInfo = 1
ResErr0 =  ResidualError(phifm, psifm, eigv, 1, ng, PE)
CONTINUE
DO iter = 1, nIterMax
  ItrCntl%Cmfdit = ItrCntl%Cmfdit + 1
  DO jsweep =1, nGroupInfo
    GrpBeg = 1; GrpEnd = ng
    IF(jsweep .GT. 1) THEN
      GrpBeg = GroupInfo%UpScatRange(1); GrpEnd = GroupInfo%UpScatRange(2)
    ENDIF
    DO ig = GrpBeg, GrpEnd
    !Source Update
      CALL CmfdSrcUpdt(Src, Psifm, phifm, eigv, ig)
    !Solving Linear System
      !Initial Guess   
      CALL CP_VA(Phic1g(1:nxy, myzbf : myzef), Phifm(1:nxy, myzbf : myzef, ig), nxy, myzef - myzbf + 1)
      CALL BiCGSTAB(CmfdLs(ig), Phic1g, SRC, itrcntl%InSolverItrCntl)
!      CALL CP_VA(Phifm(1:nxy, myzbf : myzef, ig), Phic1g(1:nxy, myzbf : myzef), nxy, myzef - myzbf + 1)    
#ifndef MPI_ENV
      CALL CP_VA(Phifm(1:nxy, myzbf : myzef, ig), Phic1g(1:nxy, myzbf : myzef), nxy, myzef - myzbf + 1)    
#else
      CALL CP_VA(Phifm(1:nxy, myzbf-1 : myzef+1, ig), Phic1g(1:nxy, myzbf-1 : myzef+1), nxy, myzef - myzbf + 3)    
#endif
    ENDDO    
  ENDDO
  CALL CP_VA(PsifmD(1 : nxy, myzbf : myzef), Psifm(1 : nxy, myzbf : myzef), nxy, myzef - myzbf + 1) 
  !Fission Source Update
  CALL CmfdPsiUpdt(phifm, psifm)
  !EigenValue Update
  peigv = eigv
  CALL CmfdEigUpdate(psifm, psifmd, eigv, psierr, PE)
  !Error Estimati
  eigerr = abs(eigv - peigv)/eigv
  ResErr =  ResidualError(phifm, psifm, eigv, 1, ng, PE)
  IF(iter .eq. 1) ResErr0 = ResErr
  WRITE(mesg,'(A9, I9, F22.6, 3(1pE12.3))') 'MGOUTER', ItrCntl%Cmfdit, eigv, eigerr, psierr, ResErr
  IF(CMFDMaster) CALL message(io8, FALSE, TRUE, mesg)  

  IF(.NOT. lfeedback .AND. ResErr .lt. convcrit .AND. eigerr .lt. EigConvCrit) lconv = TRUE
  IF(lfeedback .and. (ResErr/ResErr0) .lt. ResReduction) lconv = TRUE
  IF(iter .gt. nitermin .and. lconv) EXIT

  IF(mod(iter, 5) .EQ. 0 .and. nTracerCntl%lGcCmfd) THEN
    CALL GcCmfdAcc(Core, CmInfo, Eigv, TRUE, GroupInfo, GcGroupInfo, nTracerCntl, ItrCntl, PE) 
    CALL CP_VA(PsifmD(1 : nxy, myzbf : myzef), Psifm(1 : nxy, myzbf : myzef), nxy, myzef - myzbf + 1)    
    CALL CmfdPsiUpdt(phifm, psifm)
  ENDIF 
  
  IF(mod(iter, 5) .eq. 0 .and. nTracerCntl%l3dim) THEN
    WRITE(mesg, '(2x, A32, 2x,A10)') 'Performing Axial Calculation : ',AxSolverName(AxSolver)
    IF(CMFDMaster) CALL message(io8, FALSE, TRUE, mesg) 
    !Axail Dhat Update
    CALL MpiAxialSolver(Core, CMInfo, GroupInfo, Eigv, ng, AxSolver, .FALSE., .FALSE., PE, nTracerCntl%lnonFuelpinFDM)
!#ifndef MPI_ENV
!    CALL MpiAxialSolver(Core, CMInfo, Eigv, ng, AxSolver, FALSE, .FALSE., PE)
!#else
!  IF(PE%nProc .GT. 1) THEN
!    CALL MpiAxialSolver(Core, CMInfo, GroupInfo, Eigv, ng, AxSolver, .FALSE., .FALSE., PE)
!  ELSE
!    CALL AxialSolver(Core, CMInfo, GroupInfo, Eigv, ng, AxSolver, .FALSE., FALSE, PE) 
!  ENDIF
!#endif
    ItrCntl%Axit = ItrCntl%Axit + 1
    !LINEAR SYSTEM
    CALL SetCmfdLinearSystem(TRUE, TRUE, AxSolver)
  ENDIF  
ENDDO

CALL ConvertSubPlnPhi(PhiC, PhiFm, 2)

IF(nTracerCntl%l3dim) THEN
  CALL AxSrcUpdate(Core, CmInfo, myzb, myze, 1, ng, PE, AxSolver)
  CALL DcplAxSrcUpdt(Core, CmInfo, DcplInfo, ng, PE, FALSE)
! CALL DcplMocSolUpdt(Core, CmInfo, DcplInfo, ng, PE, FALSE)
ENDIF
!
IF(lFeedBack) THEN
  IF(ResErr .lt. convcrit .AND. eigerr .lt. EigConvCrit) ItrCntl%lConv = TRUE
ENDIF
END SUBROUTINE

SUBROUTINE DcplMocSolUpdt(Core, CmInfo, DcplInfo, ng, PE, lXsFtn)
USE PARAM
USE TYPEDEF, ONLY : CoreINfo_Type   ,CmInfo_Type      ,PE_TYPE      &
                   ,DcplInfo_Type   ,PE_Type                        &    
                   ,Pin_Type        ,Cell_Type
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(DcplInfo_Type) :: DcplInfo
TYPE(PE_Type) :: PE
INTEGER :: ng
LOGICAL :: lXsFtn

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
REAL, POINTER :: PhiC(:, :, :), PHIS(:, :, :)

REAL :: phisum, volsum, frac

INTEGER :: nxy, nRefPln
INTEGER :: myRefPlnBeg, myRefPlnEnd
INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nLocalFsr
INTEGER :: iRefpln, ixy, ig, iz0, icel, ifsr
INTEGER :: j

nxy = Core%nxy; nRefPln = DcplInfo%nRefPln 
PhiC => CmInfo%PhiFM; Cell => Core%CellInfo
Pin => Core%pin
#ifndef MPI_ENV
myRefPlnBeg = 1; myRefPlnEnd = nRefPln
#else
myRefPlnBeg = PE%myRefPlnBeg; myRefPlnEnd = PE%myRefPlnEnd
#endif
IF(.NOT. lXsFtn) THEN
  DO iRefPln = myRefPlnBeg, myRefPlnEnd
    Phis => DcplInfo%DcplFmInfo(1, iRefPln)%phis
    iz0 = DcplInfo%RefPln(iRefpln)
    IF(.NOT. Core%lFuelPlane(iz0)) CYCLE
    DO ig = 1, ng
      DO ixy = 1, nxy
        FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
        icel = Pin(ixy)%Cell(iz0)
        nlocalFxr = Cell(icel)%nFxr; nlocalFsr = Cell(icel)%nFsr     
        volsum = 0; phisum = 0
        DO j = 1, nLocalFsr
          ifsr = FsrIdxSt + j -1
          volsum = volsum + Cell(icel)%vol(j)
          phisum = phisum + Cell(icel)%vol(j) * phis(ifsr, iz0, ig)          
        ENDDO
        phisum = phisum / volsum
        frac = PhiC(ixy, iz0, ig) / phisum
        DO j = 1, nLocalFsr
          ifsr = FsrIdxSt + j -1
          phis(ifsr, iz0, ig) = phis(ifsr, iz0, ig) * frac
        ENDDO        
      ENDDO
    ENDDO
  ENDDO
ENDIF
END SUBROUTINE