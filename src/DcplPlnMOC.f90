#include <defines.h>
SUBROUTINE RefPlnXsGeneration(Core, DcplInfo, THInfo, GroupInfo, DcplCntl, DcplItrCntl, PE, DcplPE, lSubGrp0, lMOC)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type       ,DcplInfo_Type        ,GroupInfo_Type    &
                          ,ThInfo_Type         ,PE_Type                                 &
                          ,FmInfo_Type         ,CmInfo_Type          ,FxrInfo_Type
USE CNTL,            ONLY : nTracerCntl_Type
USE itrcntl_mod,     ONLY : ItrCntl_TYPE                                                &
                           ,DcplXsGenCntl_Type
USE RAYS,            ONLY : RayInfo
USE DcplXsGen_Mod,   ONLY : DcplPlnMOC
USE SubGrp_mod,      ONLY : SubGrpFsp
USE FILES,           ONLY : io8
USE IOUTIL,          ONLY : message
USE MPIComm_Mod,     ONLY : MPI_SYNC
USE MPIDcpl_Mod,     ONLY : MPI_DcplMessage     ,MPI_FinalizeDcplMesg
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(DcplInfo_Type) :: DcplInfo
TYPE(THInfo_Type), POINTER :: THInfo(:)
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: DcplCntl
TYPE(ItrCntl_Type) :: DcplItrCntl(100)
TYPE(PE_TYPE) :: PE
TYPE(PE_TYPE) :: DcplPE(:)
LOGICAL :: lSubGrp0, lMOC
INTEGER :: myrank, myRefPlnBeg, myRefPlnEnd
LOGICAL :: lMPI, Master

TYPE(CmInfo_Type), POINTER :: DcplCmInfo(:, :)
TYPE(FmInfo_Type), POINTER :: DcplFmInfo(:, :)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)

INTEGER :: iRefPln, iRefTemp
INTEGER :: nRefpln, nRefTemp, ng

LOGICAL :: lSubGrp, lXsLib

INTEGER :: i, j, k, i0, j0, k0
lXslib = DcplCntl%lXsLib;
lSubGrp = lXsLib .AND. lSubGrp0

ng = GroupInfo%ng
nRefTemp = DcplInfo%nRefTemp; nRefPln = DcplInfo%nRefPln
DcplCmInfo => DcplInfo%DcplCmInfo; DcplFmInfo => DcplInfo%DcplFmInfo
!
myrank = PE%myCMFDrank; lMPI = PE%lDcplParallel
#ifndef MPI_ENV
myRefPlnBeg = 1; myRefPlnEnd = nRefPln
#else
myRefPlnBeg = PE%myRefPlnBeg; myRefPlnEnd = PE%myRefPlnEnd
Master = PE%Master
IF(PE%Master) THEN
  WRITE(mesg, '(A)') 'Decoupled Planar MOC'
  CALL message(io8, TRUE, TRUE, mesg)
ENDIF
#endif
DO iRefPln = myRefPlnBeg, myRefPlnEnd
  iRefTemp = 1
  WRITE(mesg, '(A, I5, A, I5)') 'Decoupled Planar MOC - Ref. Plane No.', iRefPln ,'/', nRefPln
#ifndef MPI_ENV
  CALL message(io8, TRUE, TRUE, mesg)
#else
  IF(.NOT. lMPI) CALL message(io8, TRUE, TRUE, mesg)
  IF(lMPI) CALL MPI_DcplMessage(mesg, myrank, TRUE, TRUE, FALSE)
#endif
  IF(lSubGrp) CALL SubGrpFsp(Core, DcplFmInfo(iRefTemp, iRefPln)%Fxr, THInfo(iRefPln), RayInfo, GroupInfo, DcplCntl, &
                             DcplPE(iRefPln))  
  DcplItrCntl(iRefPln)%DcplXsGenCntl%iRefPln = iRefPln
  DcplItrCntl(iRefPln)%DcplXsGenCntl%iRefTemp = iRefTemp
  CALL DcplPlnMOC(Core, DcplFmInfo(iRefTemp, iRefPln), DcplCmInfo(iRefTemp, iRefPln),                   &
                  ThInfo(iRefPln), RayInfo, GroupInfo, DcplInfo%eigv(iRefTemp, iRefPln),        &
                  ng, DcplCntl, DcplItrCntl(iRefPln), DcplPE(iRefPln), .TRUE.)
  WRITE(mesg, '(A)') hbar1(1:77)
#ifndef MPI_ENV
  CALL message(io8, FALSE, TRUE, MESG)
#else
  IF(.NOT. lMPI) CALL message(io8, FALSE, TRUE, mesg)
  IF(lMPI) THEN
    CALL MPI_DcplMessage(mesg, myrank, FALSE, FALSE, TRUE)
    IF(Master) CALL message(io8, FALSE, TRUE, mesg)
  ENDIF
  CALL MPI_SYNC(PE%MPI_DcplMaster_Comm)
  IF(lMPI .AND. Master) CALL MPI_FinalizeDcplMesg(PE%nCmfdProc)
  CALL MPI_SYNC(PE%MPI_DcplMaster_Comm)
#endif
ENDDO


END SUBROUTINE

FUNCTION DcplConvChk(ItrCntl, nRefPln)
USE PARAM
USE ItrCntl_mod,   ONLY : ItrCntl_Type
IMPLICIT NONE
TYPE(ITrCntl_Type) :: ItrCntl(nRefPln)
INTEGER :: nRefPln
INTEGER :: i

LOGICAL :: DcplConvChk
LOGICAL :: ThConv, MocConv

DcplConvChk = TRUE
ThConv = TRUE; MocConv = TRUE
DO i = 1, nRefPln
  ThConv = ThConv .AND. ItrCntl(i)%lThConv
  MocConv = MocConv .AND. ItrCntl(i)%lConv
ENDDO

DcplConvChk = MocConv .AND. ThConv

END FUNCTION


SUBROUTINE DcplSetMocAxEff(Core, DcplPxs, phis1g, AxSrc1g, AxPxs1g, iz)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type                                           &
                          ,Pin_Type         ,Cell_Type
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
REAL, POINTER :: AxSrc1g(:), AxPxs1g(:), DcplAxSrc(:)
REAL :: DcplPxs(:), phis1g(:)
INTEGER :: iz
!
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
REAL :: phisum, volsum
INTEGER :: nxy, nFsr
INTEGER :: FsrIdxSt
INTEGER :: ifsr, icel, ipin
INTEGER :: j, i
nxy = core%nxy
CellInfo => Core%CellInfo; Pin => Core%Pin

!CALL CP_CA(DcplPxs(1:nxy), -0.001_8, nxy)
CALL CP_CA(AxSrc1g(1:nxy), 0._8, nxy); CALL CP_CA(AxPxs1g(1:nxy), 0._8, nxy)
!Core, DcplPxs, DcplAxSrc, phis, iz1, iz2 ,ng)
!IF(.NOT. lEigProb) RETURN
DO ipin=1,nxy
  IF(DcplPxs(ipin) .LT. 0._8) THEN
    FsrIdxSt = Pin(ipin)%FsrIdxSt; icel = Pin(ipin)%Cell(iz);
    volsum = 0; phisum = 0
    DO j = 1, CellInfo(icel)%nFsr
      ifsr = FsrIdxSt + j - 1
      volsum = volsum + CellInfo(icel)%vol(j)
      phisum = phisum + CellInfo(icel)%vol(j) * phis1g(ifsr)
      !xstr1g(ifsr) = xstr1g(ifsr)+AxPXS(ipin)
    ENDDO
    phisum = phisum / volsum
    AxSrc1g(ipin) =  DcplPxs(ipin) * phisum
  ELSE
    AxPxs1g(ipin) = DcplPxs(ipin)
  ENDIF
  !AxPxs1g(ipin) = DcplPxs1g(ipin)
ENDDO
END SUBROUTINE



SUBROUTINE DcplPlnMOC(Core, FmInfo, CmInfo, ThInfo, RayInfo, GroupInfo, EigV,       &
                      ng, DcplCntl, ItrCntl, DcplPE, lMOC)
USE PARAM
USE TYPEDEF,         ONLY : DcplInfo_Type     ,CoreInfo_Type     ,RayInfo_Type      &
                           ,FmInfo_Type       ,CmInfo_Type       ,THInfo_Type       &
                           ,GroupInfo_Type    ,PE_TYPE
!USE MOC_Mod,         ONLY : DcplMocSweep
USE CMFD_mod,        ONLY : HomoXsGen         ,RadCouplingCoeffGen
USE CNTL,            ONLY : nTracerCntl_Type
USE ItrCntl_Mod,     ONLY : ItrCntl_Type
!Subroutine or Function Definition
USE SubGrp_mod,      ONLY : SubGrpEffXsGen    ,FxrChiGen
USE DcplXsGen_Mod,   ONLY : DcplMOCSweep      ,DcplPlnCmfd       ,DcplEffXsSave
#ifdef MPI_ENV
USE MPIDcpl_mod,     ONLY : Idle_DcplMoc      ,Idle_DcplCmfd
USE MPIComm_mod,     ONLY : MPIConvChk
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(ThInfo_Type) :: ThInfo

TYPE(RayInfo_Type) :: RayInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_TYPE) :: DcplPE
TYPE(nTracerCntl_Type) :: DcplCntl
TYPE(ItrCntl_Type) :: ItrCntl
LOGICAL :: lMOC
INTEGER :: ng
REAL :: Eigv, RefTemp

INTEGER :: iter
INTEGER :: iRefPln, nXsGenIterMaxfPln, iRefTemp
INTEGER :: nXsGenIterMax
#ifdef MPI_ENV
INTEGER :: Comm, nproc, myrank
LOGICAL :: GlobalConv, lMPI
#endif

nXsGenIterMax =  ItrCntl%DcplXsGenCntl%nXsGenIterMax
iRefPln = ItrCntl%DcplXsGenCntl%iRefPln
iRefTemp = ItrCntl%DcplXsGenCntl%iRefTemp

#ifdef MPI_ENV
Comm = DcplPE%MPI_DCPLMASTER_COMM
myrank = DcplPE%myrank; nproc = DcplPE%nproc
GlobalConv = .FALSE.
lMPI = DcplPE%lDcplParallel
#endif

CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, Eigv, GroupInfo, DcplCntl, DcplPE)
CALL DcplEffXsSave(Core, FmInfo%Fxr, iRefTemp, GroupInfo, DcplPE)
ItrCntl%lconv = FALSE
!IF(ItrCntl%MocIt .EQ. 0) CALL DcplPlnCmfd(Core, CmInfo, FmInfo, Eigv, ng, TRUE, TRUE, GroupInfo, DcplPE, DcplCntl, ItrCntl)
DO iter = 1, nXsGenIterMax
#ifndef MPI_ENV
  IF(lMOC) CALL DcplMocSweep(Core, RayInfo, FmInfo, Eigv, ng, DcplPE, GroupInfo, DcplCntl, ItrCntl)
  CALL DcplPlnCmfd(Core, CmInfo, FmInfo, Eigv, ng, TRUE, TRUE, GroupInfo, DcplPE, DcplCntl, ItrCntl)

  IF(ItrCntl%lconv) EXIT
#else
  IF(.NOT. ItrCntl%lConv) THEN
    IF(lMOC) CALL DcplMocSweep(Core, RayInfo, FmInfo, Eigv, ng, DcplPE, GroupInfo, DcplCntl, ItrCntl)
    CALL DcplPlnCmfd(Core, CmInfo, FmInfo, Eigv, ng, TRUE, TRUE, GroupInfo, DcplPE, DcplCntl, ItrCntl)
    !Global Conv Check
    GlobalConv = MPIConvChk(ItrCntl%lConv, nproc, comm)
    IF(GlobalConv) EXIT
  ELSEIF(.NOT. GlobalConv) THEN
    IF(lMOC) CALL Idle_DcplMoc(Eigv, ng, GroupInfo, ItrCntl, DcplPE)
    CALL Idle_DcplCmfd(EIgv, ItrCntl, DcplPE)
    GlobalConv = MPIConvChk(ItrCntl%lConv, nproc, comm)
    IF(GlobalConv) EXIT
  ENDIF
#endif
ENDDO
CALL HomoXsGen(Core, FmInfo%Fxr, FmInfo%Phis, CmInfo%PinXS, DcplPE%myzb,     &
               DcplPE%myzb, ng, DcplCntl%lXsLib, FALSE, FALSE)
CALL RadCouplingCoeffGen(Core, CmInfo%PinXs, CmInfo%RadJout, ng, .TRUE., DcplPE)
!CALL DcplPlnCmfd(Core, CmInfo, FmInfo, Eigv, ng, TRUE, TRUE, GroupInfo, DcplPE, DcplCntl, ItrCntl)



CALL FxrChiGen(Core, FmInfo%Fxr, FmInfo, GroupInfo, DcplPE, DcplCntl)

END SUBROUTINE

SUBROUTINE DcplMOCSweep(Core, RayInfo, FMInfo, eigv, ng, PE, GroupInfo, nTracerCntl, ItrCntl)
USE PARAM
USE TYPEDEF,     ONLY : CoreInfo_Type,      RayInfo_Type,       FmInfo_Type,        &
                        PE_TYPE,            FxrInfo_Type,       GroupInfo_Type
USE CNTL,        ONLY : nTracerCntl_Type
USE itrcntl_mod, ONLY : ItrCntl_TYPE
!USE CORE_MOD,    ONLY : GroupInfo
USE MOC_MOD,     ONLY : RayTrace,           SetRtMacXs,        SetRtSrc,          &
                        PsiUpdate,          CellPsiUpdate,     UpdateEigv,        &
                        MocResidual,        PsiErr,            PseudoAbsorption,  &
                        AddConstSrc,                                              &
                        phis1g,             MocJout1g,         xst1g,             &
                        tSrc,               AxSrc1g,           AxPxs1g,           &
                        PhiAngIn1g,         srcm
USE DcplXsGen_Mod,  ONLY : DcplSetMocAxEff
USE BasicOperation, ONLY : CP_VA
USE SUbGrp_Mod,     ONLY : FxrChiGen
USE IOUTIL,     ONLY : message
USE Timer,      ONLY : nTracer_dclock, TimeChk
USE FILES,      ONLY : io8
#ifdef MPI_ENV
USE MPIDcpl_Mod, ONLY : MPI_DcplMessage,   PrintDcplRTLog,   PrintDcplMOCLog
#endif
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(PE_TYPE) :: PE
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(ItrCntl_TYPE) :: ItrCntl
REAL :: eigv
INTEGER :: ng

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
REAL, POINTER :: PHIS(:,:,:)                          !Scalar Flux Vectors
REAL, POINTER :: PSI(:,:), PSID(:,:)                  !Fsr FIssion Scource
REAL, POINTER :: PSIC(:,:), PSICD(:,:)                !Cell Fission Source Term
REAL, POINTER :: PhiAngin(:,:,:,:)                    !Angular Flux data at the boundary
REAL, POINTER :: AxSrc(:, :, :), AxPXS(:, :, :)
REAL, POINTER :: RadJout(:, :, :, :, :)

REAL, POINTER :: LinSrcSlope(:, :, :, :)
REAL, POINTER :: phim(:, :, :, :)

!Iteration Control
INTEGER :: nitermax, nitermin
REAL :: eigconv, psiconv, ResConv

!Local Variable
INTEGER :: ig, iz, l, ITER, jsweep, InIter, nInIter
INTEGER :: myzb, myze
INTEGER :: nFSR, nxy, nbd, GrpBeg, GrpEnd, nGroupInfo
REAL :: psipsi, psipsid, eigerr,fiserr, peigv, reserr
REAL :: TimeRt1gBeg, TimeRt1gEnd, TimeElasped, MocTimeBeg, MocTimeEnd
LOGICAL :: lJout, lxslib, l3dim, lscat1, lTrCorrection, lRST, lssph, lssphreg
LOGICAL :: lDmesg                  !Detailed Message
INTEGER :: myrank, comm
LOGICAL :: lEigProb, lNegFix, lMPI
!CMFD Time Check
MOCTimeBeg = nTracer_dclock(FALSE, FALSE)

!Pointing
FXR => FmInfo%Fxr
PHIS => FmInfo%PHIS; PSI => FmInfo%PSI
PSID => FmInfo%PSID; PSIC => FmInfo%PSIC
PSICD => FmInfo%PSICD; PhiAngin => FmInfo%PhiAngin
RadJout => FmInfo%RadJout;
AxSrc => FmInfo%AxSrc; AxPXS => FmInfo%AxPXS

!AxSrc => FmInfo%AxSrc; AxPXS => FmInfo%AxPXS
!MPI Variables
lMPI = PE%lDcplParallel
myrank = PE%myrank; comm = PE%MPI_DCPLMASTER_COMM
!Local Variable
nbd = 4
myzb = PE%myzb; myze = PE%myze
nFsr = Core%nCoreFsr; nxy = Core%nxy
iter = 0; nGroupInfo =2
nInIter = 2

!Control Variable
lJout = TRUE
lxslib = nTracerCntl%lXslib; l3dim = nTracerCntl%l3dim
lscat1 = nTracerCntl%lscat1
lTrCorrection= nTracerCntl%lTrCorrection
lssph = nTracerCntl%lsSPH
lRST = nTracerCntl%lRST
lssphreg = nTracerCntl%lsSPHreg
!Eigen Value Proglem?
lEigProb = TRUE
IF(.NOT. Core%lFuelPlane(myzb)) lEigProb = FALSE
lNegFix = .NOT. lEigProb
!Iteration Control Variable
nitermax = itrcntl%MOCItrCntl%nitermax;nitermin = itrcntl%MOCItrCntl%nitermin
psiconv = itrcntl%DcplXsGenCntl%psiconv; eigconv = itrcntl%DcplXsGenCntl%eigconv
ResConv = itrcntl%DcplXsGenCntl%ResConv

CALL CP_VA(psid(1:nFsr, myzb:myze), psi(1:nFsr, myzb:myze), nFsr, myze - myzb +1)
CALL CP_VA(psicd(1:nxy, myzb:myze), psic(1:nxy, myzb:myze), nxy,  myze - myzb +1)
CALL PsiUpdate(Core, Fxr, phis, psi, myzb, myze, ng, lxslib, GroupInfo)
CALL CellPsiUpdate(Core, Psi, psic, myzb, myze)
!Output Message Control
lDmesg = TRUE
IF(ng .gt. 10) lDmesg = FALSE
IF(.NOT. GroupInfo%lUpScat) nGroupInfo = 1
#ifdef MPI_ENV
IF(myrank .EQ. 0 .AND. lMPI) THEN
  WRITE(mesg, '(a26)') 'Performing Ray Tracing ...'
  CALL message(io8, TRUE, TRUE, mesg)
ENDIF
#endif

DO iter = 1, 1
  itrcntl%mocit = itrcntl%mocit + 1; TimeElasped = 0.
  WRITE(mesg, '(a22,I5,a3)') 'Performing Ray Tracing', itrcntl%mocit, '...'
#ifndef MPI_ENV
  CALL message(io8, TRUE, TRUE, mesg)
#else
  IF(.NOT. lMPI) CALL message(io8, FALSE, TRUE, mesg)
  IF(lMPI) CALL MPI_DcplMessage(mesg, myrank, TRUE, FALSE, FALSE)
#endif
  CALL FxrChiGen(Core, Fxr, FmInfo, GroupInfo, PE, nTracerCntl)
  DO jsweep =1, nGroupInfo
    GrpBeg = 1; GrpEnd = ng
    !IF(.NOT. GroupInfo%lUpScat .AND. jsweep .GT. 1) EXIT
    IF(jsweep .GT. 1) THEN
      GrpBeg = GroupInfo%UpScatRange(1); GrpEnd = GroupInfo%UpScatRange(2)
    ENDIF

    DO ig = GrpBeg, GrpEnd

      TimeRt1gBeg = nTracer_dclock(FALSE, FALSE)
      iz = myzb
      ljout = FALSE

      DO InIter = 1, 2 !InIter
        CALL DcplSetMocAxEff(Core, AxPxs(:, iz, ig), phis(:, iz, ig), AxSrc1g, AxPxs1g, iz)
        IF(InIter .EQ. nInIter) ljout = TRUE

       !IF(InIter .EQ. 1 .OR. lNegFix) THEN
          CALL SetRtMacXs(Core, Fxr(:, iz), xst1g, iz, ig, ng, lxslib, lTrCorrection, lRST, lsSPH, lsSPHreg, PE)
          CALL PseudoAbsorption(Core, Fxr(:, iz), tsrc, phis(:, iz, ig),                 &
                                AxPXS1g(:), xst1g, iz, ig, ng, GroupInfo, TRUE)
       ! ENDIF
        !SET MOC Source

        CALL SetRtSrc(Core, Fxr(:, iz), tsrc, phis, psi, axSrc1g, xst1g,                           &
                        eigv, iz, ig, ng, GroupInfo, TRUE, lXslib, lscat1, lNegFix, PE)
        IF(.NOT. lEigProb) THEN
          !nTracerCntl%ConstSrc = 1.0
          CALL AddConstSrc(Core, Fxr(:, iz), tsrc, xst1g, nTracerCntl%ConstSrc, iz, ig, ng)
        ENDIF
        !Add Constant source terms for a Reflector Region
        CALL CP_VA(PhiAngin1g, PhiAngin(:,: ,iz, ig), RayInfo%nPolarAngle, RayInfo%nPhiAngSv)
        CALL RayTrace(RayInfo, Core, phis1g, PhiAngIn1g, xst1g, tsrc, MocJout1g, iz, lJout)
        CALL CP_VA(phis(1:nFsr, iz, ig), phis1g(1:nFsr), nFsr)
      ENDDO
      !IF(lJout) CALL CP_VA(RadJout(1:2, 1:nbd, 1:nxy, iz, ig), MocJout1g(1:2, 1:nbd, 1:nxy), 2, nbd, nxy)
      IF(lJout) CALL CP_VA(RadJout(1:3, 1:nbd, 1:nxy, iz, ig), MocJout1g(1:3, 1:nbd, 1:nxy), 3, nbd, nxy)   !---BYS edit / 150612 Surface flux
      !ENDDO
      !Edit Print out Log messages
      TimeRt1gEnd = nTracer_dclock(FALSE, FALSE)
      TimeElasped = TimeElasped + TimeRt1gEnd - TimeRt1gBeg
      IF(lDmesg .or. MOD(iG, 10) .eq. 0 .OR. ig .EQ. nG) THEN
        write(mesg,'(10x, a, i4, 2x, a, f10.3, 2x, a)') 'Group ', ig, ' finished in ', TimeElasped, 'Sec'
#ifndef MPI_ENV
        CALL message(io8, FALSE, TRUE, mesg)
#else
        IF(.NOT. lMPI) CALL message(io8, FALSE, TRUE, mesg)
        IF(lMPI) THEN
          CALL MPI_DcplMessage(mesg, myrank, FALSE, FALSE, FALSE)
          CALL PrintDcplRTLog(ig, TimeElasped, myrank, PE%nproc, Comm)
        ENDIF
#endif
        TimeElasped = 0.
      ENDIF

    ENDDO
    CONTINUE
  ENDDO
  !Update Fission Source
    CALL CP_VA(psid(1:nFsr, myzb:myze), psi(1:nFsr, myzb:myze), nFsr, myze - myzb +1)
    CALL CP_VA(psicd(1:nxy, myzb:myze), psic(1:nxy, myzb:myze), nxy,  myze - myzb +1)
    CALL PsiUpdate(Core, Fxr, phis, psi, myzb, myze, ng, lxslib, GroupInfo)
  IF(lEigProb) THEN
    CALL UpdateEigv(Core, psi, psid, eigv, peigv, myzb, myze)
    CALL CellPsiUpdate(Core, Psi, psic, myzb, myze)
    fiserr = PsiErr(Core, Psi, PsiD, myzb, myze, PE)
    eigerr = abs((eigv - peigv))/eigv
  ELSE
    fiserr = 0; eigerr = 0
  ENDIF
  reserr = MocResidual(Core, FmInfo, eigv, GroupInfo, ng, PE, nTracerCntl)
  write(mesg ,'(A5,I7,F15.6, 3(1pE12.3))')  'RT',itrcntl%mocit, EIGV, eigerr, fiserr, reserr
 ! IF(iter .gt. 5) exit
#ifndef MPI_ENV
  CALL message(io8, TRUE, TRUE, mesg)
#else
  IF(.NOT. lMPI) CALL message(io8, TRUE, TRUE, mesg)
  IF(lMPI) THEN
    CALL MPI_DcplMessage(mesg, myrank, TRUE, FALSE, FALSE)
    CALL PrintDcplMOCLog(itrcntl%mocit, myzb, EIGV, eigerr, fiserr, reserr, myrank, PE%nproc, Comm)
  ENDIF
#endif
  CONTINUE
 ! IF(eigerr .lt. epsm6 .and. fiserr .lt. epsm5) exit
ENDDO
IF(lEigProb .AND. fiserr .lt. psiconv .and. eigerr .lt. eigconv) itrcntl%lconv = TRUE
IF(.NOT. lEigProb .AND. ResErr .lt. ResConv) itrcntl%lconv = TRUE
MocTimeEnd = nTracer_dclock(FALSE, FALSE)
TimeElasped = MocTimeEnd - MocTimeBeg
TimeChk%MocTime = TimeChk%MocTime + TimeElasped
END SUBROUTINE
