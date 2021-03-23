#include <defines.h>
MODULE CrCspGen_MOD
USE PARAM
USE TYPEDEF,           ONLY : PinXs_Type,      CoreInfo_Type,      PE_TYPE,    &
                              GroupInfo_Type
IMPLICIT NONE
TYPE CrCspGenCntl_Type
  INTEGER :: npinxsdat = 0                      !
  INTEGER :: BankID                             !
  INTEGER :: CrPln(2)                           !Cr Insertion Range
  CHARACTER(256) :: PinXsFn(100)
  INTEGER, POINTER :: RodOutMap(:), RodInMap(:)
END TYPE

TYPE IntSpec_TYPE
  INTEGER :: nFxr, ng
  REAL, POINTER :: phi(:, :) !(nfxr, ig)
END TYPE

TYPE CrCspGenDat_Type
  REAL :: eigv
  INTEGER :: iCrPos(2)
  REAL :: CrPos
  REAL, POINTER :: Phi(:, :), Phi4th(:, :, :)
END TYPE

TYPE(CrCspGenCntl_Type) :: CrCspGenCntl
!TYPE(CoreInfo_Type), POINTER :: Core
TYPE(IntSpec_Type), POINTER, PRIVATE :: IntSpec(:)
TYPE(PinXs_Type), POINTER, PRIVATE :: PinXsDat(:, :)

TYPE(CrCspGenDat_Type), POINTER, PRIVATE :: CrCspGenDat(:)
INTEGER, PRIVATE :: nxy, nz, ng, myzb, myze, myzbf, myzef


CONTAINS


SUBROUTINE InitCrCspGen(Core, GroupInfo, PE)

IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_Type) :: PE


nxy = Core%nxy; nz = Core%nz
myzb = PE%myzb; myze = PE%myze
myzbf = PE%myzbf; myzef = PE%myzef
ng = GroupInfo%ng
CALL Alloc_CrCspGen(Core, GroupInfo)
CALL ReadPinXsDat(GroupInfo)
END SUBROUTINE

SUBROUTINE Alloc_CrCspGen(Core, GroupInfo)
USE CMFD_mod,        ONLY : AllocPinXS
USE ALLOCS
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(GroupInfo_Type) :: GroupInfo

REAL :: h(0:nz)

INTEGER :: i, j
INTEGER :: nbd

!ALLOCATE PinXS DATA
nbd = 4
ALLOCATE(PinXsDat(nxy, CrCspGenCntl%nPinXsDat))
ALLOCATE(IntSpec(CrCspGenCntl%nPinXsDat))
DO i = 1, CrCspGenCntl%nPinXsDat
  DO j = 1, nxy
    CALL AllocPinXs(PinXsDat(j, i), ng, nbd, GroupInfo%InScatRange)
  ENDDO
ENDDO

!ALLOCATE CrCspGenDat
h(0) = 0
DO i = 1, nz
  h(i) = h(i-1) + Core%hz(i)
ENDDO


ALLOCATE(CrCspGenDat(CrCspGenCntl%CrPln(1):CrCspGenCntl%CrPln(2)))
DO i = CrCspGenCntl%CrPln(1), CrCspGenCntl%CrPln(2)
  CALL DMALLOC(CrCspGenDat(i)%Phi, nz, ng)
  CALL DMALLOC0(CrCspGenDat(i)%Phi4th, 0, 4, 1, nz, 1, ng)
  CrCspGenDat(i)%iCrPos = (/i, CrCspGenCntl%CrPln(2)/)
  CrCspGenDat(i)%CrPos = h(i-1)
ENDDO

END SUBROUTINE

SUBROUTINE ReadPinXsDat(GroupInfo)
USE files,         ONLY : localfn,             io15
USE IOUTIL,        ONLY : OpenFile,            TERMINATE
USE ALLOCS
IMPLICIT NONE
TYPE(GroupInfo_Type) :: GroupInfo
INTEGER :: i, j, k, ixy
INTEGER :: io
INTEGER :: ig, ib, ie, ib1, ie1
INTEGER :: nxy0, ng0, iz0

io = io15
DO i = 1, CrCspGenCntl%nPinXsDat
  localfn = CrCspGenCntl%PinXsFn(i)
  CALL OpenFile(io15, .TRUE., .TRUE., .FALSE., localFn )
  READ(io) nxy0, ng0, iz0
  IF(nxy0 .NE. nxy) CALL TERMINATE('ReadPinXsDat : nxy inconsistency')
  IF(ng0 .NE. ng) CALL TERMINATE('ReadPinXsDat : ng inconsistency')
  !IF(nxy0 .NE. nxy) CALL TERMINATE('ReadPinXsDat : iz0 inconsistency')
  DO ig = 1, ng
    READ(io) ib, ie, ib1, ie1
    IF(ib .NE. GroupInfo%InScatRange(1, ig)) CALL TERMINATE('ReadPinXsDat : Scattering Range Inconsistency')
    IF(ie .NE. GroupInfo%InScatRange(2, ig)) CALL TERMINATE('ReadPinXsDat : Scattering Range Inconsistency')
    IF(ib1 .NE. GroupInfo%OutScatRange(1, ig)) CALL TERMINATE('ReadPinXsDat : Scattering Range Inconsistency')
    IF(ie1 .NE. GroupInfo%OutScatRange(2, ig)) CALL TERMINATE('ReadPinXsDat : Scattering Range Inconsistency')
    
  ENDDO
  DO ixy = 1, nxy
    READ(io) PinXSDat(ixy, i)%XSD(1:ng)
    READ(io) PinXSDat(ixy, i)%XST(1:ng)
    READ(io) PinXSDat(ixy, i)%XSTR(1:ng)
    READ(io) PinXSDat(ixy, i)%XSR(1:ng)
    READ(io) PinXSDat(ixy, i)%XSNF(1:ng)
    READ(io) PinXSDat(ixy, i)%XSKF(1:ng)
    READ(io) PinXSDat(ixy, i)%CHI(1:ng)
    DO j = 1, 4
      READ(io) PinXsDat(ixy, i)%DTIL(j, 1:ng)
    ENDDO
    
    DO j = 1, 4
      READ(io) PinXsDat(ixy, i)%DHAT(j, 1:ng)
    ENDDO
    DO j = 1, ng
        READ(io) ib, ie
        READ(io) PinXSDat(ixy, i)%xss(j)%WithInGroupScat
        READ(io) PinXSDat(ixy, i)%xss(j)%from(ib:ie)
    ENDDO
  ENDDO
  READ(io) IntSpec(i)%nFxr
  CALL Dmalloc(IntSpec(i)%phi, IntSpec(i)%nFxr, ng)
  DO j = 1, ng
    READ(io) ig, (IntSpec(i)%phi(k ,j), k = 1, IntSpec(i)%nFxr)
  ENDDO
  CLOSE(io)
ENDDO
END SUBROUTINE

SUBROUTINE SetCspPinXS(PINXS, CspXsMAP, PE)
IMPLICIT NONE
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
TYPE(PE_Type) :: PE
INTEGER :: CspXsMAP(nz)
INTEGER :: iz, ixy, idat

DO iz = myzb, myze
  idat = CspXsMap(iz)
  DO ixy = 1, nxy
    CALL CopyPinXS(PinXS(ixy, iz), PinXSDat(ixy, idat))
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE MovingCRPosition(CmInfo, iCrLoc, PE)
USE TYPEDEF,               ONLY : CMInfo_Type
USE BasicOperation,        ONLY : CP_VA
IMPLICIT NONE
TYPE(CmInfo_Type) :: CmInfo
TYPE(PE_TYPE) :: PE
INTEGER :: iCrLoc

INTEGER :: XSMAP(nz)
INTEGER :: iz
CALL CP_VA(XSMAP(1:nz), CrCspGenCntl%RodOutMap(1:nz), nz)
DO iz = iCrLoc, nz
  XSMAP(iz) = CrCspGenCntl%RodInMap(iz)
ENDDO

CALL SetCspPinXS(CmInfO%PinXs, XSMAP, PE)
END SUBROUTINE

SUBROUTINE CopyPinXS(PinXS1, PinXS2)
TYPE(PinXS_TYPE), INTENT(INOUT) :: PinXS1
TYPE(PinXS_TYPE), INTENT(in) :: PinXS2
INTEGER :: ig, ib, ie

PinXS1%XSD(1:ng) = PinXS2%XSD(1:ng); !PinXS1%XSD2(1:ng) = PinXS2%XSD(1:ng)
PinXS1%XST(1:ng) = PinXS2%XST(1:ng)
PinXS1%XSTR(1:ng) = PinXS2%XSTR(1:ng); PinXS1%XSR(1:ng) = PinXS2%XSR(1:ng)
PinXS1%XSNF(1:ng) = PinXS2%XSNF(1:ng); PinXS1%XSKF(1:ng) = PinXS2%XSKF(1:ng)
PinXS1%CHI(1:ng) = PinXS2%CHI(1:ng)

DO ig = 1, ng
  PinXS1%DHAT(1:4, ig) = PinXS2%DHAT(1:4, ig)
  PinXS1%DTIL(1:4, ig) = PinXS2%DTIL(1:4, ig)
  PinXS1%XSD2(ig) = 3._8 / (7._8 * PinXS1%XSD(ig))
ENDDO

!
DO ig = 1, ng
  PinXS1%Xss(ig)%withinGroupScat = PinXS2%Xss(ig)%withinGroupScat
  ib = PinXS2%Xss(ig)%ib; ie = PinXS2%Xss(ig)%ie
  PinXS1%xss(ig)%from(ib:ie) = PinXS2%xss(ig)%from(ib:ie)
ENDDO
END SUBROUTINE

SUBROUTINE CrCspGenCmfd(Core, CMInfo, eigv, PE, nTracerCntl, ItrCntl)
USE PARAM
USE TYPEDEF,        ONLY : PE_TYPE,         CoreInfo_Type,      Pin_Type,                   &
                           CMInfo_Type,     FmInfo_Type,                                    &
                           FxrInfo_Type,    CmfdLs_TYPE,        PinXs_Type,                 &
                           AxFlx_Type
USE CMFD_mod,       ONLY : CmfdPinXS,        CmfdLS,                                        &
                           PhiC1g,            SRC,                                          &
                           SetCMFDEnviorment, HomoXsGen,         RadCouplingCoeffGen,       &
                           UpdatePhiC,        ConvertSubPlnPhi,                             &
                           CmfdSrcUpdt,       CmfdPsiUpdt,                                  &
                           CmfdEigUpdate,     ResidualError
USE CORE_MOD,       ONLY : GroupInfo,         GcGroupInfo 
USE BiLU_Mod
USE MpiAxSolver_mod,ONLY : MpiAxialSolver
#ifdef MPI_ENV
USE MpiComm_mod,    ONLY : MPI_SYNC
#endif
!USE Boron_Mod,      ONLY : UpdtBoronPPM,      UpdtBoronCmfdXS
USE BiCGSTAB_mod,   ONLY : BiCGSTAB
USE BasicOperation, ONLY : CP_CA, CP_VA
USE IOUTIL,         ONLY : message
USE FILES,          ONLY : io8
USE CNTL,           ONLY : nTracerCntl_Type
USE itrcntl_mod,    ONLY : ItrCntl_TYPE, CMFDItrCntl_TYPE
use timer,          only : nTracer_dclock, TimeChk
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CMInfo_Type) :: CMInfo
TYPE(PE_TYPE) :: PE
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(ItrCntl_TYPE) :: ItrCntl
REAL :: eigv
REAL :: CmfdTimeBeg, CmfdTimeEnd


INTEGER :: iter, ig, jsweep, GrpBeg, GrpEnd
INTEGER :: nitermax, nitermin, AxSolver, nGroupInfo, nAxUpdt
REAL :: ResErr0, ResErr, convcrit, reigv, peigv, psierr, eigerr

LOGICAL :: lscat1, lsubplane, lDtilOnly, lconv, lexit
LOGICAL :: CmfdMaster, CmfdSlave
LOGICAL, SAVE :: lfirst

!Pointing Varibables
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
!TYPE(AxFlx_Type), POINTER :: AxFlx(:, :)
REAL, POINTER :: PHIS(:,:,:)     
REAL, POINTER :: PhiC(:, :, :), PSIC(:,:), PSICD(:,:)
REAL, POINTER :: PhiFm(:, :, :), PsiFm(:, :), PsiFmD(:, :)

DATA lfirst /TRUE/

CMFDMaster = PE%CMFDMaster; CMFDSlave = PE%CMFDSlave
!CMFD Time Check
CmfdTImeBeg = nTracer_dclock(FALSE, FALSE)

!Pointing
Pin => Core%Pin; 
PinXS => CMInfo%PinXS
PhiC => CMInfo%PhiC; PhiFM => CMInfo%PhiFM
PsiC => CMInfo%PsiC; PsiCD => CMInfo%PsiCD 
PsiFm => CMInfo%PsiFm; PsiFmD => CMInfo%PsiFmD

lscat1 = nTracerCntl%lScat1
lsubplane = nTracerCntl%lSubPlane
AxSolver = nTracerCntl%AxSolver

nitermax = 200
nitermin = 15
!convcrit = itrcntl%CMFDItrCntl%convcrit

!IF(ItrCntl%MocIt .LT. 1) AxSolver = LP1SENM
IF(lfirst) THEN
  CALL SetCMFDEnviorment(Core, CmInfo, ng, PE)
  CALL SetCmfdMpiOmpEnv(Core, CMInfo%CoreCmfdLs , ng, PE)
  lfirst = FALSE
ENDIF
WRITE(mesg, '(a)') 'Performing CMFD Calculation...'
IF(CMFDMaster) CALL message(io8, TRUE, TRUE, mesg)   
!CALL CP_CA(PhiFM(1:nxy, myzb-1:myze+1,1:ng), 1._8,  nxy, myze-myzb+3, ng)
!CALL ConvertSubPlnPhi(PhiC, PhiFm, 1)
!CALL  InitIterVar(Core, FMInfo, CmInfo, GroupInfo,   .TRUE., ItrCntl, PE)
IF(nTracerCntl%l3dim) THEN
  IF(CMFDMaster) WRITE(mesg, '(2x, A32, 2x,A10)') 'Performing Axial Calculation : ',AxSolverName(AxSolver)
  IF(CMFDMaster) CALL message(io8, FALSE, TRUE, mesg)
  lDtilOnly = .FALSE.;
  IF(ItrCntl%Cmfdit .EQ. 0) lDtilOnly = .TRUE.
  CALL MpiAxialSolver(Core, CMInfo, GroupInfo, Eigv, ng, AxSolver, .TRUE., lDtilOnly, PE, nTracerCntl%lnonFuelpinFDM)
ENDIF


!LINEAR SYSTEM
CALL SetCmfdLinearSystem(TRUE,  nTracerCntl%l3dim, AxSolver)
CALL MakeBiLU(Core, CmInfo%CoreCMFDLS(1:ng), PhiFm, ng, PE)
CALL CmfdPsiUpdt(PhiFm, PsiFm) 
lconv = FALSE;
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
  !IF((ResErr/ResErr0) .lt. convcrit ) lconv = TRUE
  IF(ResErr .LT. 1.0e-4) lConv = .TRUE.
  lExit = lExit .AND. lConv
  IF(iter .LE. nitermin) lExit = .FALSE.
  IF(iter .GE. nitermin .AND. ResErr .LT. 1.0E-7_8) lExit = .TRUE.
  IF(lExit) EXIT 
  IF(nitermax .EQ. iter) EXIT

  !CAlling 1D-Axial Kernel
  IF(mod(iter, 7) .eq. 0 .and. nTracerCntl%l3dim) THEN
    nAxUpdt = nAxUpdt  + 1
    IF(CMFDMaster) WRITE(mesg, '(2x, A32, 2x,A10)') 'Performing Axial Calculation : ',AxSolverName(AxSolver)
    IF(CMFDMaster) CALL message(io8, FALSE, TRUE, mesg) 
    !Axail Dhat Update
     CALL MpiAxialSolver(Core, CMInfo, GroupInfo, Eigv, ng, AxSolver, FALSE, .FALSE., PE, nTracerCntl%lnonFuelpinFDM) 

    
    CALL SetCmfdLinearSystem(TRUE, TRUE, AxSolver)
    CALL MakeBiLU(Core, CmInfo%CoreCMFDLS(1:ng), PhiFm, ng, PE)
  ENDIF
  IF(mod(iter, 7) .EQ. 0 .and. nTracerCntl%lGcCmfd) THEN
    CALL GcCmfdAcc(Core, CmInfo, Eigv, TRUE, GroupInfo, GcGroupInfo, nTracerCntl, ItrCntl, PE) 
    CALL CP_VA(PsifmD(1 : nxy, myzbf : myzef), Psifm(1 : nxy, myzbf : myzef), nxy, myzef - myzbf + 1)    
    CALL CmfdPsiUpdt(phifm, psifm)
  ENDIF   
END DO
CmfdTImeEnd = nTracer_dclock(FALSE, FALSE)
TimeChk%CmfdTime = TimeChk%CmfdTime + (CmfdTimeEnd - CmfdTimeBeg)
!Free
NULLIFY(Pin); 
NULLIFY(PinXS);

NULLIFY(PhiC); NULLIFY(PhiFM)
NULLIFY(PsiC); NULLIFY(PsiCD)
NULLIFY(PsiFm); NULLIFY(PsiFmD)
END SUBROUTINE

SUBROUTINE SaveCrCspFlux(Core, CmInfo, eigv, iBank, iCrPos, nTracerCntl, PE)
USE TYPEDEF,          ONLY : CmInfo_Type,          Pin_Type
USE CNTL,             ONLY : nTracerCntl_Type
USE BasicOperation,   ONLY : CP_CA,                CP_VA
USE MPIAxSolver_Mod,  ONLY : Get1dFlx4th
#ifdef MPI_ENV
USE MpiComm_mod,      ONLY : REDUCE
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
REAL :: eigv
INTEGER :: iBank, iCrPos

REAL, POINTER :: Phi4th(:, :, :)

TYPE(Pin_TYPE), POINTER :: Pin(:)
REAL, POINTER :: PinVol(:, :), PhiC(:, :, :)
REAL, POINTER :: buf(:, :)
REAL :: VolSum, fnorm
INTEGER :: ixy, iz, ig
INTEGER :: n

Pin => Core%Pin; PinVol => Core%PinVol; PhiC => CmInfo%PhiC
CrCspGenDat(iCrPos)%eigv = eigv
n = 0

fnorm = FluxNorm(Core, CmInfo, nTracerCntl, PE)
CALL CP_CA(CrCspGenDat(iCrPos)%Phi, 0._8, nz, ng)
DO ig = 1, ng
  DO iz = myzb, myze
    n = 0; volsum =0
    DO ixy = 1, nxy
      IF(.NOT. Pin(ixy)%lCrPin) CYCLE
      IF(Pin(ixy)%CrBankID .NE. iBank) CYCLE 
      n = n + 1
      volsum = volsum + PinVol(ixy, iz)
      CrCspGenDat(iCrPos)%Phi(iz, ig) = CrCspGenDat(iCrPos)%Phi(iz, ig) + PhiC(ixy, iz, ig) * PinVol(ixy, iz)
    ENDDO
    CrCspGenDat(iCrPos)%Phi(iz, ig) = fnorm * CrCspGenDat(iCrPos)%Phi(iz, ig) / volsum
  ENDDO
ENDDO
#ifdef MPI_ENV
ALLOCATE(Buf(nz, ng))
CALL CP_CA(Buf, 0._8, nz, ng)
CALL REDUCE(CrCspGenDat(iCrPos)%Phi, buf, nz, ng, PE%MPI_CMFD_COMM, .TRUE.)
CALL CP_VA(CrCspGenDat(iCrPos)%Phi, buf, nz, ng)
DEALLOCATE(Buf)
#endif

ALLOCATE(Phi4th(0:4, nz, ng))
n = 0
CALL CP_CA(CrCspGenDat(iCrPos)%Phi4th(0:4, 1:nz, 1:ng), 0._8, 5, nz, ng)
DO ixy = 1, nxy
  IF(.NOT. Pin(ixy)%lCrPin) CYCLE
  IF(Pin(ixy)%CrBankID .NE. iBank) CYCLE 
  CALL Get1dFlx4th(Phi4th, ixy, nTracerCntl, PE)
  n = n + 1
  DO ig = 1, ng
    DO iz = 1, nz
      CrCspGenDat(iCrPos)%Phi4th(0:4, iz, ig) = CrCspGenDat(iCrPos)%Phi4th(0:4, iz, ig) + PinVol(ixy, iz) * Phi4th(0:4, iz, ig)
    ENDDO
  ENDDO
ENDDO

DO ixy = 1, nxy
  IF(.NOT. Pin(ixy)%lCrPin) CYCLE
  IF(Pin(ixy)%CrBankID .NE. iBank) CYCLE 
  DO ig = 1, ng
    DO iz = 1, nz
      fnorm = CrCspGenDat(iCrPos)%Phi(iz, ig) / CrCspGenDat(iCrPos)%Phi4th(0, iz, ig)
      CrCspGenDat(iCrPos)%Phi4th(0:4, iz, ig) = fnorm * CrCspGenDat(iCrPos)%Phi4th(0:4, iz, ig) 
    ENDDO
  ENDDO
ENDDO
DEALLOCATE(Phi4th)
NULLIFY(PinVol, Pin, PhiC)
END SUBROUTINE

SUBROUTINE WriteCrCspFile_BIN(Core, PE)
!File Format
USE files,            ONLY : CaseID,              localfn,             io15
USE IOUTIL,           ONLY : OpenFile,            TERMINATE
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(PE_TYPE) :: PE

REAL :: h(0:nz)

INTEGER :: i, j, k
INTEGER :: icrin, icrout
INTEGER :: icrbeg, icrend, icase, iz, ig

IF(.NOT. PE%MASTER) RETURN
icrbeg = CrCspGenCntl%CrPln(1); icrend = CrCspGenCntl%CrPln(2)

h(0) = 0
DO i = 1, nz
  h(i) = h(i-1) + Core%hz(i)
ENDDO

localfn = trim(caseid)//'.crcsp'
CALL OpenFile(io15, .FALSE., .TRUE., .FALSE., localfn)


WRITE(io15) icrend - icrbeg + 1, icrbeg, icrend, nz, ng
WRITE(io15) (h(i), i = 0, nz)

icase = 0
DO i = icrbeg, icrend
  icase = icase+1
  WRITE(io15) icase, i, h(i-1)
  DO ig = 1, ng
    WRITE(io15) (CrCspGenDat(i)%PHI(iz, ig), iz = 1, nz)
  ENDDO
ENDDO

icase = 0
DO i = icrbeg, icrend
  icase = icase+1
  WRITE(io15) icase, i, h(i-1)
  DO ig = 1, ng
    WRITE(io15) (CrCspGenDat(i)%PHI4th(0:4, iz, ig), iz = 1, nz)
  ENDDO
ENDDO

DO i = icrbeg, icrend
  icrin = CrCspGenCntl%RodInMap(i)
  icrout = CrCspGenCntl%RodOutMap(i)
  WRITE(io15) IntSpec(icrout)%nfxr
  DO ig = 1, ng
    WRITE(io15) (IntSpec(icrin)%phi(j, ig), j = 1, IntSpec(icrin)%nfxr)
    WRITE(io15) (IntSpec(icrout)%phi(j, ig), j = 1, IntSpec(icrout)%nfxr)
  ENDDO
ENDDO
CLOSE(io15)

END SUBROUTINE

SUBROUTINE WriteCrCspFile_ascii(Core, PE)
!File Format
USE files,         ONLY : CaseID,              localfn,             io15
USE IOUTIL,        ONLY : OpenFile,            TERMINATE
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(PE_TYPE) :: PE

REAL :: h(0:nz)

INTEGER :: i, j, k
INTEGER :: icrbeg, icrend, icase, iz, ig


IF(.NOT. PE%MASTER) RETURN
icrbeg = CrCspGenCntl%CrPln(1); icrend = CrCspGenCntl%CrPln(2)

h(0) = 0
DO i = 1, nz
  h(i) = h(i-1) + Core%hz(i)
ENDDO

localfn = trim(caseid)//'.crcsp'
CALL OpenFile(io15, .FALSE., .FALSE., .FALSE., localfn)
!Writing Header!

WRITE(io15, '(A)') 'Control Rod Cusping Function Data File V1.0'
WRITE(io15, '(A)')
WRITE(io15, '(A11, x, F15.5, 3x, A2, F15.5)') 'CR Region :', h(icrbeg-1), 'to', h(icrend)
WRITE(io15,  '(A11, x , I5)') 'No Case :' , icrend - icrbeg + 1

WRITE(io15, '(A)') 'DATA'
WRITE(io15, '(100I5)') icrend - icrbeg + 1, icrbeg, icrend, nz, ng
WRITE(io15, '(10F15.5)') (h(i), i = 0, nz)
WRITE(io15, *)
icase = 0
DO i = icrbeg, icrend
  icase = icase+1
  WRITE(io15, '(A10, 2I5, F10.5)') 'CASE', icase, i, h(i-1)
  DO iz = 1, nz
    WRITE(io15, '(100(1pe20.8))') (CrCspGenDat(i)%PHI(iz, ig), ig = 1, ng)
  ENDDO
  WRITE(io15, *)
ENDDO
icase = 0
WRITE(io15, '(A)') '4TH Order Flux'
DO i = icrbeg, icrend
  icase = icase+1
  WRITE(io15, '(A10, 2I5, F10.5)') 'CASE', icase, i, h(i-1)
  DO ig = 1 ,ng
    DO iz = 1, nz
      WRITE(io15, '(100(1pe20.8))') CrCspGenDat(i)%PHI4TH(0:4, iz, ig)
    ENDDO
  ENDDO
  WRITE(io15, *)
ENDDO

CLOSE(io15)
END SUBROUTINE

FUNCTION FluxNorm(Core, CmInfo, nTracerCntl, PE)
! Flux level normalization
! normalize flux such that average flux in fuel region be unity
! then update fission source and moments accordingly

USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,     FmInfo_Type,       CmInfo_Type,     &
                           PE_TYPE,                                               &
                           PinXs_Type
USE CNTL,           ONLY : nTracerCntl_Type
USE BasicOperation, ONLY : MULTI_CA
#ifdef MPI_ENV
USE MPIComm_Mod,    ONLY : REDUCE
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
REAL :: FluxNorm


TYPE(PinXs_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: PhiC(:, :, :), PsiC(:, :)
REAL, POINTER :: PinVol(:, :)

INTEGER :: nbd
INTEGER :: ixy, iz, ig

REAL :: vol, volsum, pwsum, phisum, pinpw
REAL :: avgpw, avgflx, norm
REAL :: UnitPowerLevel0, Plevel

INTEGER :: comm
REAL :: buf0(3), buf(3)

COMM = PE%MPI_CMFD_COMM

nbd = 4
PhiC => CmInfo%PhiC; PsiC => CmInfo%PsiC

PinVol => Core%PinVol; PinXs => CmInfo%PinXs

volsum = 0;
pwsum = 0;
phisum = 0;
DO iz = myzb, myze
  DO ixy = 1, nxy 
    pinpw = 0
    DO ig = 1, ng
      PinPw = PinPw + PinXs(ixy, iz)%xskf(ig) * PhiC(ixy, iz, ig) 
    ENDDO
    IF(PinPw .LT. 0._8) CYCLE
    vol = PinVol(ixy, iz)
    phisum = phisum + sum(PhiC(ixy, iz, 1:ng)) * vol
    PwSum = PwSum + PinPw * vol
    volsum = volsum + vol
  ENDDO
ENDDO

#ifdef MPI_ENV
buf0 = (/PwSum, PhiSum, volsum/)
CALL REDUCE(buf0, buf, 3, COMM, .TRUE.)
PwSum = Buf(1); PhiSum = Buf(2); Volsum = Buf(3)
#endif

AvgPw = PwSum / VolSum
avgflx = PhiSum / Volsum

norm = ng / avgflx
FluxNorm = norm

END FUNCTION
END MODULE
