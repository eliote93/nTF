#include <defines.h>
MODULE MPIAxSolver_Mod
USE PARAM
USE TYPEDEF,        ONLY : AxFlx_TYPE,   PinXS_TYPE,     GroupInfo_TYPE,   PE_TYPE, AxGeom_Type
!USE Core_Mod,       ONLY : GroupInfo
USE TIMER,       ONLY : nTracer_dclock, TimeChk
#ifdef MPI_ENV
USE MPIComm_Mod,    ONLY : MPI_SYNC  !,     InitAxNodalComm
USE CMFDComm_Mod,   ONLY : InitAxNodalComm
#endif
!USE AxSolver_Mod,   ONLY : AxGeom_Type
IMPLICIT NONE

TYPE AxNCntl_TYPE
  INTEGER :: ng0
  INTEGER :: AxSolverMod = 1
  REAL :: eigv
  LOGICAL :: lReset = .FALSE.
  LOGICAL :: lDtilUpdt = .FALSE.
  LOGICAL :: lTransient = .FALSE.
  LOGICAL :: lXsComm = .FALSE.
  
END TYPE

INTEGER, PRIVATE :: nxy, myzb, myze, myzbf, myzef, nzfm, nz, ng
INTEGER, PRIVATE :: myNxyBeg, myNxyEnd
INTEGER, SAVE, PRIVATE :: sendrecvmax = 100
INTEGER, SAVE, PRIVATE :: nSendPin
REAL, POINTER :: hz(:), hzfm
INTEGER, POINTER :: SubPlaneMap(:)

TYPE(AxGEOM_TYPE), PRIVATE :: AxGeom
TYPE(PinXS_TYPE), POINTER, SAVE, PRIVATE :: AxPinXS(:, :)
TYPE(AXFlx_Type), POINTER, SAVE, PRIVATE :: AxFlx(:, :)
REAL, POINTER, SAVE, PRIVATE :: LKG(:, :, :)
REAL, POINTER, SAVE, PRIVATE :: MpiLkg(:, :, :), MpiPhi(:, :, :)
REAL, POINTER, SAVE, PRIVATE :: SendBufD(:, :, :, :), RecvBufD(:, :, :, :)

TYPE(AxFlx_Type), POINTER, PRIVATE :: FlxOMP(:, :)
TYPE(PinXS_Type), POINTER, PRIVATE :: XS(:, :)
REAL, POINTER, PRIVATE :: AxBeta(:, :, :), AxOmega0(: ,:, :), AxOmegap(:, :, :), AxOmegam(:, :, :)
REAL, POINTER, PRIVATE :: TranSrc2nd(:, :, :), MpiTranSrc2nd(:, :, :)
#ifdef MPI_ENV
CONTAINS

SUBROUTINE AllocMpiAxSol(AxSolver, GroupInfo, ng0, PE)
USE P1SENM_MOD,   ONLY : AllocP1Solver
USE Sp3Senm_Mod,  ONLY : AllocSP3Solver
USE CMFD_mod,     ONLY : AllocPinXs
USE CNTL,         ONLY : nTracerCntl
USE geom,         ONLY : nbd
USE ALLOCS
IMPLICIT NONE
INTEGER :: AxSolver, ng0
TYPE(GroupInfo_TYPE) :: GroupInfo
TYPE(PE_TYPE) :: PE

INTEGER :: ixybeg, ixyend
!INTEGER ::  izbeg, izend
INTEGER :: nz, nzfm, nxy, myzbf, myzef
INTEGER :: iz, ixy

nz = PE%nz
nzfm = PE%nzfm
nxy = PE%nxy
myzbf = PE%myzbf; myzef = PE%myzef
ixybeg = PE%mynxybeg; ixyend = PE%mynxyend
!izbeg = 1; izend = nz

ALLOCATE(AxFlx(1:nzfm, ixybeg:ixyend))
CALL Dmalloc0(Lkg, 1, nxy, myzbf, myzef, 1, ng0)

SELECT CASE(AxSolver)
  CASE(lP1SENM)
    CALL AllocP1Solver(AxFlx, ixybeg, ixyend, 1, nzfm, ng0, .TRUE.)
  CASE(lP3SENM)
    CALL AllocSP3Solver(AxFlx, ixybeg, ixyend, 1, nzfm, ng0, .TRUE.)
END SELECT
IF(PE%nCmfdProc .GT. 1) THEN
  ALLOCATE(AxPinXS(ixybeg:ixyend, 1:nz))
  DO iz = 1, nz
    DO ixy = ixybeg, ixyend
      CALL AllocPinXS(AxPinXs(ixy, iz), ng0, nbd, GroupInfo%InScatRange)
    ENDDO 
  ENDDO
  CALL InitAxNodalComm(GroupInfo, ng0, 4, SendRecvMax, PE, PE%NonBlockPinXsComm)
  CALL Dmalloc0(MpiLkg, ixybeg, ixyend, 1, nzfm, 1, ng0)
  CALL Dmalloc0(MpiPhi, ixybeg, ixyend, 1, nzfm, 1, ng0)
ELSE

ENDIF


END SUBROUTINE

SUBROUTINE InitMPIAxNSolver(Core, GroupInfo, ng0, AxSolver, PE)
USE PARAM
USE Typedef,      ONLY : PE_TYPE, CoreInfo_Type, AxFlx_Type, GroupInfo_Type, PE_TYPE
USE P1SENM_MOD,   ONLY : AllocP1AxFlxType
USE SP3SENM_MOD,  ONLY : AllocSP3AxFlxType
USE GEOM,         ONLY : nbd
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_TYPE) :: PE
INTEGER :: ng0, AxSolver
INTEGER :: iz, itid

!LOGICAL :: lfirst 
!DATA lfirst / TRUE /
!
!IF(.NOT. lfirst) THEN
!  DEALLOCATE(XS)
!ENDIF
!DATA lfirst  /.false./

ng = ng0; nxy = Core%nxy
nz = PE%nz; nzfm = PE%nzfm
myzb = PE%myzb; myze = PE%myze
myzbf = PE%myzbf; myzef = PE%myzef
SubPlaneMap => Core%SubPlaneMap
!
AxGeom%ng = ng0
AxGeom%ncomp = nzfm; AxGeom%nMesh =nzfm
AxGeom%myzb = 1; AxGeom%myze = nz
AxGeom%myzbf = 1; AxGeom%myzef = nzfm
AxGeom%BC = Core%AxBC
AxGeom%H => Core%Hzfm; AxGeom%HInv => Core%HzfmInv
AxGeom%Comp => Core%SubPlaneMap

CALL AllocMpiAxSol(AxSolver, GroupInfo, ng0, PE)

ALLOCATE(XS(Core%nzfm, PE%nAxThread))
ALLOCATE(FLXOMP(Core%nzfm, PE%nAxThread))
DO itid = 1, PE%nAxThread
  DO iz = 1, Core%nzfm
    IF(AxSolver .EQ. lP1SENM) THEN
      CALL AllocP1AxFlxType(FLXOMP(iz, itid), ng0, .TRUE.)
    ELSEIF(AxSolver .EQ. lP3SENM) THEN
      CALL AllocSp3AxFlxType(FLXOMP(iz, itid), ng0, .TRUE.)
    ENDIF
    CALL AllocPinXS(XS(iz, itid), ng0, nbd, GroupInfo%InScatRange)
  ENDDO
ENDDO
!IF(PE%nCmfdProc .EQ. 1) MpiLkg => Lkg


END SUBROUTINE



SUBROUTINE MpiAxialSolverNew(Core, CmInfo, GroupInfo, Eigv, ng0, AxSolverMod, lXsComm, lDtilGen, PE, lUr)
USE PARAM
USE TYPEDEF,          ONLY : CMInfo_Type,       CoreInfo_Type,          PE_TYPE,         &
                             GroupInfo_Type
!USE AxSolver_Mod,     ONLY : AxialXsMapping,    FinalizeAxialXsMapping, CalAxDtil
USE Tlkg_mod,         ONLY : RadTlkgUpdate,     SetTlkgShape
USE P1SENM_MOD,       ONLY : SetP1SenmEnv
USE Sp3Senm_Mod,      ONLY : SetSP3SenmEnv
USE CMFDComm_mod,     ONLY : CommPinXS_Type,                                             &
                             CommGlobalPinXsDat,   CommAxNodalInPut,                     &
                             CommAxNodalOutPut, CommAxNodalOutPutNew
USE TIMER,            ONLY : nTracer_dclock, TimeChk
USE BasicOperation,   ONLY : CP_VA
USE IOUTIL,            ONLY : TERMINATE
!USE itrcntl_mod,    ONLY : ItrCntl_TYPE
USE OMP_LIB
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CMInfo_Type) :: CMInfo
TYPE(GroupInfo_Type) :: GroupInfo
!TYPE(AxFlx_Type), POINTER :: Flx1D(:)
TYPE(PE_TYPE) :: PE
REAL :: eigv, peigv
INTEGER :: ng0, AxSolverMod
REAL :: AxNTimeBeg, AxNTimeEnd, ElaspedTime
LOGICAL :: lXsComm, lDtilGen, lUR

TYPE(CommPinXS_Type) :: CommPinXS

!TYPE(AxFlx_TYPE), POINTER :: FLX(:, :)
!TYPE(PinXS_TYPE), POINTER :: PinXS(:, :)

INTEGER :: iz, iz0, ixy, ig

INTEGER,SAVE :: iter

INTEGER :: lkgmod, tid

LOGICAL, SAVE :: lFirst
DATA lfirst /TRUE/
DATA iter /0/

iter = iter + 1

IF(PE%CMFDMaster) AxNTimeBeg = nTracer_dclock(FALSE, FALSE)

IF(lfirst) THEN
  !CALL SetMpiAxSolverEnvironment(CORE, GroupInfo, ng0, AxSolverMod, PE)
  IF(AxSolverMod .EQ. lP1SENM) CALL SetP1SenmEnv(AxGeom)
  IF(AxSolverMod .EQ. lP3SENM) CALL SetSP3SenmEnv(AxGeom, PE)
  !Communication module initialization
  !CALL InitAxNodalComm(GroupInfo, ng, 4, SendRecvMax, PE, .FALSE.)
  lfirst = FALSE
ENDIF

myNxyBeg = PE%myNxyBeg; myNxyEnd = PE%myNxyEnd

!Determine Transverse Leakage Information
CALL RadTlkgUpdate(Core, CmInfo, Lkg, 1, nxy, myzbf, myzef, ng)

!IF(PE%CMFDMaster)THEN
!    WRITE(175,*) myNxyBeg, myNxyEnd
!DO iz = 1, nz
!  DO ixy = myNxyBeg, myNxyEnd
!!DO iz = myzbf, myzef
!!  DO ixy = myNxyBeg, myNxyEnd
!      WRITE(175,*) iz, ixy
!      DO ig = 1, ng
!          WRITE(175,*) ig, AxFlx(iz, ixy)%Dhat(1, ig), AxFlx(iz, ixy)%Dhat(2, ig)
!      ENDDO
!  ENDDO
!ENDDO
!ELSE
!    WRITE(176,*) myNxyBeg, myNxyEnd
!DO iz = 1, nz
!  DO ixy = myNxyBeg, myNxyEnd
!!DO iz = myzbf, myzef
!!  DO ixy = myNxyBeg, myNxyEnd
!      WRITE(176,*) iz, ixy
!      DO ig = 1, ng
!          WRITE(176,*) ig, AxFlx(iz, ixy)%Dhat(1, ig), AxFlx(iz, ixy)%Dhat(2, ig)
!      ENDDO
!  ENDDO
!ENDDO
!ENDIF

IF(PE%nCmfdProc .GT. 1) THEN
!IF(PE%CMFDMaster)  WRITE(*,*) PE%nCmfdProc
  CommPinXS%SendPhi => CmInfo%PhiFm; CommPinXS%SendLkg => Lkg
  CommPinXS%SendPinXs => CmInfo%PinXS
  CommPinXS%RecvPhi => MpiPhi; CommPinXs%RecvLkg => MpiLkg
  CommPinXS%RecvPinXs => AxPinXS        !
  IF(lXsComm) CALL CommGlobalPinXsDat(CommPinXS, ng, GroupInfo, PE)
  !Destroy the CommPinXS variables
  NULLIFY(CommPinXS%SendPhi, CommPinXS%SendLkg, CommPinXS%SendPinXS)
  NULLIFY(CommPinXS%RecvPhi, CommPinXS%RecvLkg, CommPinXS%RecvPinXS)
  !Get Lkg and flux for Radial Domain
  CALL CommAxNodalInPut(CmInfo%PhiFm, Lkg, MpiPhi, MpiLkg, ng, PE)
ELSE
  MpiPhi => CmInfo%PhiFm 
  AxPinXs => CmInfo%PinXS
  MpiLkg => Lkg
ENDIF
!Update Average Flux info to the Axial Nodal structure
CALL PhiAvgUpdate(MpiPhi, AxFlx)
CALL TlkgAvgUpdate(Mpilkg, AxFlx)

lkgmod = 0
tid = 1 

!print *, tid
!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(PE%nAxThread) 
!$OMP PARALLEL DEFAULT(SHARED)  &
!$OMP PRIVATE(tid, ixy)
!$  tid = omp_get_thread_num()+1
!$OMP DO 
DO ixy = myNxyBeg, myNxyEnd
  !IF( .NOT.Core%Pin(ixy)%lfuel ) CYCLE
  !IF( iter .LT. 10 ) CYCLE
      CALL MpiNodal1D_OMP(AxFlx(:, ixy), AxPinXS(ixy, :), Eigv, AxGeom, AxSolverMod, lDtilGen, tid, PE)
  !ENDIF  
ENDDO
!$OMP END DO
!$OMP END PARALLEL

!IF(PE%CMFDMaster)THEN
!    WRITE(175,*) myNxyBeg, myNxyEnd
!DO iz = 1, nz
!  DO ixy = myNxyBeg, myNxyEnd
!!DO iz = myzbf, myzef
!!  DO ixy = myNxyBeg, myNxyEnd
!      WRITE(175,*) iz, ixy
!      DO ig = 1, ng
!          WRITE(175,*) ig, AxFlx(iz, ixy)%Dhat(1, ig), AxFlx(iz, ixy)%Dhat(2, ig)
!      ENDDO
!  ENDDO
!ENDDO
!ELSE
!    WRITE(176,*) myNxyBeg, myNxyEnd
!DO iz = 1, nz
!  DO ixy = myNxyBeg, myNxyEnd
!!DO iz = myzbf, myzef
!!  DO ixy = myNxyBeg, myNxyEnd
!      WRITE(176,*) iz, ixy
!      DO ig = 1, ng
!          WRITE(176,*) ig, AxFlx(iz, ixy)%Dhat(1, ig), AxFlx(iz, ixy)%Dhat(2, ig)
!      ENDDO
!  ENDDO
!ENDDO
!ENDIF

CALL MPI_SYNC(PE%MPI_nTRACER_COMM)
!
CALL StabGapPinLkg(Core, AxFlx,  PE) 

!Axial Nodal Calculation Output Distributions


 ! DO ixy = myNxyBeg, myNxyEnd
 !   DO iz = myzbf, myzef
 !     DO ig = 1, ng
 !      IF(ISNAN(AxFlx(iz, ixy)%Dhat(1, ig))) CALL TERMINATE('NaN in Dhat')
 !      IF(ISNAN(AxFlx(iz, ixy)%Dhat(2, ig))) CALL TERMINATE('NaN in Dhat')
 !     ENDDO
 !   ENDDO
 ! ENDDO
  
! AxDhat UR
IF( lUR )THEN
IF(PE%CMFDMaster) WRITE(*,*) 'AxDhat UR'
  DO ixy = myNxyBeg, myNxyEnd
    DO iz = myzbf, myzef
      DO ig = 1, ng
        !IF(AxFlx(iz, ixy)%PDhat(1, ig).NE.0)THEN
            AxFlx(iz, ixy)%Dhat(1:2, ig) = AxFlx(iz, ixy)%Dhat(1:2, ig)*0.5+AxFlx(iz, ixy)%PDhat(1:2, ig)*0.5
            !AxFlx(iz, ixy)%Dhat(1:2, ig) = AxFlx(iz, ixy)%PDhat(1:2, ig)
            !AxFlx(iz, ixy)%Dhat(1:2, ig) = 0._8
        !ELSE
        !    AxFlx(iz, ixy)%Dhat(1:2, ig)=AxFlx(iz, ixy)%Dhat(1:2, ig)*0.1
        !ENDIF   
        ! AxFlx(iz, ixy)%Dhat(1:2, ig) = AxFlx(iz, ixy)%Dhat(1:2, ig)*0.01+AxFlx(iz, ixy)%PDhat(1:2, ig)*0.09
      ENDDO
    ENDDO
ENDDO
ENDIF

IF(PE%nCmfdProc .GT. 1) THEN
!WRITE(*,*) 'AxDHat UR ??'
  CALL CommAxNodalOutPutNew(CmInfo, AxFlx, ng, PE, lUR)

ELSE
!IF( iter .LT. 5 )THEN
!WRITE(*,*) 'AxDHat UR On'
!ELSE
!WRITE(*,*) 'AxDHat UR OFF'
!ENDIF
  DO ixy = 1, nxy
    DO iz = myzbf, myzef
      DO ig = 1, ng
        CmInfo%AxDtil(1:2, ixy, iz, ig) = AxFlx(iz, ixy)%Dtil(1:2, ig)
        !IF( iter .GT. 5 )THEN
        !    !CmInfo%AxDhat(1:2, ixy, iz, ig) = AxFlx(iz, ixy)%Dhat(1:2, ig)*0.5+AxFlx(iz, ixy)%PDhat(1:2, ig)*0.5
        !    CmInfo%AxDhat(1:2, ixy, iz, ig) = AxFlx(iz, ixy)%Dhat(1:2, ig)*0.01+AxFlx(iz, ixy)%PDhat(1:2, ig)*0.09
        !ENDIF        
        CmInfo%AxDhat(1:2, ixy, iz, ig) = AxFlx(iz, ixy)%Dhat(1:2, ig)
        AxFlx(iz, ixy)%PDhat(1:2, ig) = AxFlx(iz, ixy)%Dhat(1:2, ig)
      ENDDO
    ENDDO
  ENDDO
ENDIF

IF(PE%CMFDMaster) THEN
  AxNTimeEnd = nTracer_dclock(FALSE, FALSE)
  ElaspedTime = AxNTimeEnd - AxNTimeBeg
  TimeChk%AxialNodalTime = TimeChk%AxialNodalTime + ElaspedTime
ENDIF
CALL MPI_SYNC(PE%MPI_CMFD_COMM)

END SUBROUTINE


SUBROUTINE MpiAxialSolver(Core, CmInfo, GroupInfo, Eigv, ng0, AxSolverMod, lXsComm, lDtilGen, PE, lnonFuelpinFDM)
USE PARAM
USE TYPEDEF,          ONLY : CMInfo_Type,       CoreInfo_Type,          PE_TYPE,         &
                             GroupInfo_Type
!USE AxSolver_Mod,     ONLY : AxialXsMapping,    FinalizeAxialXsMapping, CalAxDtil
USE Tlkg_mod,         ONLY : RadTlkgUpdate,     SetTlkgShape
USE P1SENM_MOD,       ONLY : SetP1SenmEnv
USE Sp3Senm_Mod,      ONLY : SetSP3SenmEnv
USE CMFDComm_mod,     ONLY : CommPinXS_Type,                                             &
                             CommGlobalPinXsDat,   CommAxNodalInPut,                     &
                             CommAxNodalOutPut
USE TIMER,            ONLY : nTracer_dclock, TimeChk
USE BasicOperation,   ONLY : CP_VA
USE CNTL,             ONLY : nTracerCntl
USE HexTLkg,          ONLY : HexRadTlkgUpdt, HexStabGapPinLkg
USE OMP_LIB
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CMInfo_Type) :: CMInfo
TYPE(GroupInfo_Type) :: GroupInfo
!TYPE(AxFlx_Type), POINTER :: Flx1D(:)
TYPE(PE_TYPE) :: PE
REAL :: eigv, peigv
INTEGER :: ng0, AxSolverMod
REAL :: AxNTimeBeg, AxNTimeEnd, ElaspedTime
LOGICAL :: lXsComm, lDtilGen, lnonFuelpinFDM

TYPE(CommPinXS_Type) :: CommPinXS

!TYPE(AxFlx_TYPE), POINTER :: FLX(:, :)
!TYPE(PinXS_TYPE), POINTER :: PinXS(:, :)

INTEGER :: iz, iz0, ixy, ig

INTEGER,SAVE :: iter

INTEGER :: lkgmod, tid

LOGICAL, SAVE :: lFirst
DATA lfirst /TRUE/
DATA iter /0/

iter = iter + 1

IF(PE%CMFDMaster) AxNTimeBeg = nTracer_dclock(FALSE, FALSE)

IF(lfirst) THEN
  !CALL SetMpiAxSolverEnvironment(CORE, GroupInfo, ng0, AxSolverMod, PE)
  IF(AxSolverMod .EQ. lP1SENM) CALL SetP1SenmEnv(AxGeom)
  IF(AxSolverMod .EQ. lP3SENM) CALL SetSP3SenmEnv(AxGeom, PE)
  !Communication module initialization
  !CALL InitAxNodalComm(GroupInfo, ng, 4, SendRecvMax, PE, .FALSE.)
  lfirst = FALSE
ENDIF

myNxyBeg = PE%myNxyBeg; myNxyEnd = PE%myNxyEnd

!Determine Transverse Leakage Information
IF (nTracerCntl%lHex) THEN
  CALL HexRadTlkgUpdt(Core, CmInfo, Lkg, 1, nxy, myzbf, myzef, ng)
ELSE
  CALL RadTlkgUpdate(Core, CmInfo, Lkg, 1, nxy, myzbf, myzef, ng)
END IF

IF(PE%nCmfdProc .GT. 1) THEN
  CommPinXS%SendPhi => CmInfo%PhiFm; CommPinXS%SendLkg => Lkg
  CommPinXS%SendPinXs => CmInfo%PinXS
  CommPinXS%RecvPhi => MpiPhi; CommPinXs%RecvLkg => MpiLkg
  CommPinXS%RecvPinXs => AxPinXS        !
  IF(lXsComm) CALL CommGlobalPinXsDat(CommPinXS, ng, GroupInfo, PE)
  !Destroy the CommPinXS variables
  NULLIFY(CommPinXS%SendPhi, CommPinXS%SendLkg, CommPinXS%SendPinXS)
  NULLIFY(CommPinXS%RecvPhi, CommPinXS%RecvLkg, CommPinXS%RecvPinXS)
  !Get Lkg and flux for Radial Domain
  CALL CommAxNodalInPut(CmInfo%PhiFm, Lkg, MpiPhi, MpiLkg, ng, PE)
ELSE
  MpiPhi => CmInfo%PhiFm 
  AxPinXs => CmInfo%PinXS
  MpiLkg => Lkg
ENDIF
!Update Average Flux info to the Axial Nodal structure
CALL PhiAvgUpdate(MpiPhi, AxFlx)
CALL TlkgAvgUpdate(Mpilkg, AxFlx)

lkgmod = 0
tid = 1

!print *, tid
!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(PE%nAxThread) 
!$OMP PARALLEL DEFAULT(SHARED)  &
!$OMP PRIVATE(tid, ixy)
!$  tid = omp_get_thread_num()+1
!$OMP DO 
DO ixy = myNxyBeg, myNxyEnd
  IF( .NOT.Core%Pin(ixy)%lfuel .AND. lnonFuelpinFDM ) CYCLE  ! FDM for nonfuel
  !IF( iter .LT. 10 ) CYCLE
  AxGeom%ix=ixy
      CALL MpiNodal1D_OMP(AxFlx(:, ixy), AxPinXS(ixy, :), Eigv, AxGeom, AxSolverMod, lDtilGen, tid, PE)
  !ENDIF  
ENDDO
!$OMP END DO
!$OMP END PARALLEL

CALL MPI_SYNC(PE%MPI_nTRACER_COMM)
!
IF (nTracerCntl%lHex) THEN
  CALL HexStabGapPinLkg(Core, AxFlx, PE)
ELSE
  CALL StabGapPinLkg(Core, AxFlx,  PE)
END IF

!Axial Nodal Calculation Output Distributions

! AxDhat UR
!IF( iter .LT. 100 )THEN
!!IF(PE%CMFDMaster) WRITE(*,*) 'AxDhat UR'
!  DO ixy = myNxyBeg, myNxyEnd
!    DO iz = myzbf, myzef
!      DO ig = 1, ng
!        !IF(AxFlx(iz, ixy)%PDhat(1, ig).NE.0)THEN
!            !AxFlx(iz, ixy)%Dhat(1:2, ig) = AxFlx(iz, ixy)%Dhat(1:2, ig)*0.1+AxFlx(iz, ixy)%PDhat(1:2, ig)*0.9
!        !ELSE
!        !    AxFlx(iz, ixy)%Dhat(1:2, ig)=AxFlx(iz, ixy)%Dhat(1:2, ig)*0.1
!        !ENDIF   
!        ! AxFlx(iz, ixy)%Dhat(1:2, ig) = AxFlx(iz, ixy)%Dhat(1:2, ig)*0.01+AxFlx(iz, ixy)%PDhat(1:2, ig)*0.09
!      ENDDO
!    ENDDO
!ENDDO
!ENDIF

IF(PE%nCmfdProc .GT. 1) THEN
!WRITE(*,*) 'AxDHat UR ??'
  CALL CommAxNodalOutPut(CmInfo, AxFlx, ng, PE)
ELSE
!IF( iter .LT. 5 )THEN
!WRITE(*,*) 'AxDHat UR On'
!ELSE
!WRITE(*,*) 'AxDHat UR OFF'
!ENDIF
  DO ixy = 1, nxy
    DO iz = myzbf, myzef
      DO ig = 1, ng
        CmInfo%AxDtil(1:2, ixy, iz, ig) = AxFlx(iz, ixy)%Dtil(1:2, ig)
        !IF( iter .GT. 5 )THEN
        !    !CmInfo%AxDhat(1:2, ixy, iz, ig) = AxFlx(iz, ixy)%Dhat(1:2, ig)*0.5+AxFlx(iz, ixy)%PDhat(1:2, ig)*0.5
        !    CmInfo%AxDhat(1:2, ixy, iz, ig) = AxFlx(iz, ixy)%Dhat(1:2, ig)*0.01+AxFlx(iz, ixy)%PDhat(1:2, ig)*0.09
        !ENDIF        
        CmInfo%AxDhat(1:2, ixy, iz, ig) = AxFlx(iz, ixy)%Dhat(1:2, ig)
        AxFlx(iz, ixy)%PDhat(1:2, ig) = AxFlx(iz, ixy)%Dhat(1:2, ig)
      ENDDO
    ENDDO
  ENDDO
ENDIF

IF(PE%CMFDMaster) THEN
  AxNTimeEnd = nTracer_dclock(FALSE, FALSE)
  ElaspedTime = AxNTimeEnd - AxNTimeBeg
  TimeChk%AxialNodalTime = TimeChk%AxialNodalTime + ElaspedTime
ENDIF
CALL MPI_SYNC(PE%MPI_CMFD_COMM)

END SUBROUTINE


SUBROUTINE MpiNodal1D_OMP(Flx0, PinXS, Eigv, AxGeom, AxSolverMod, lDtilUpdt, tid, PE)
USE PARAM
USE TYPEDEF,          ONLY : PE_TYPE
USE UtilFunction,     ONLY : CopyP1AxFlx,         CopySP3AxFlx,            CopyAxXs
!USE AxSolver_Mod,     ONLY : CalAxDtil
USE Tlkg_mod,         ONLY : SetTlkgShape,        SetTlkgShape0
USE P1SENM_MOD,       ONLY : P1SENM
USE Sp3Senm_Mod,      ONLY : DirectSP3SENM
USE OMP_LIB

IMPLICIT NONE
TYPE(AxFlx_TYPE) :: Flx0(nzfm)
TYPE(PinXS_TYPE) :: PinXS(nz)
TYPE(AxGeom_TYPE) :: AxGeom
TYPE(PE_TYPE) :: PE
REAL :: Eigv
INTEGER :: AxSolverMod, tid
LOGICAL :: lDtilUpdt
INTEGER :: iz, iz0

!XSMapping Between Axial Fine mesh and Coarse mesh(MOC mesh)
DO iz = 1, nzfm
  iz0 = AxGeom%COMP(iz)
  CALL CopyAxXs(Xs(iz, tid), PinXS(iz0), ng)
  IF(AxSolverMod .EQ. lP1SENM) CALL CopyP1AxFlx(FlxOMP(iz, tid), Flx0(iz), AxGeom%lTransient)
  IF(AxSolverMod .EQ. lP3SENM) CALL CopySP3AxFlx(FlxOMP(iz, tid),Flx0(iz), AxGeom%lTransient)
  !CALL AxialXsMapping(CMInfo%PinXS(ixy, iz0), Xs(iz, tid))
ENDDO

IF(lDtilUpdt) THEN
  CALL CalAxDtil(FlxOMP(1:nzfm, tid), XS(1:nzfm, tid), AxGeom, nzfm, AxSolverMod)
ELSE
  CALL SetTlkgShape(FlxOMP(1:nzfm, tid), AxGeom%bc, AxGeom%h(1:nzfm), nzfm, ng, AxSolverMod, 0)
  SELECT CASE(AxSolverMod)
    CASE(lP1SENM)
      CALL P1SENM(FlxOMP(1:nzfm, tid), XS(1:nzfm, tid), Eigv, PE)
    CASE(lP3SENM)
      CALL DirectSP3SENM(FlxOMP(1:nzfm, tid), XS(1:nzfm, tid), EigV, AxGeom, tid)
  END SELECT
ENDIF
  
DO iz = 1, nzfm
  IF(AxSolverMod .EQ. lP1SENM) CALL CopyP1AxFlx(Flx0(iz), FlxOMP(iz, tid), AxGeom%lTransient)
  IF(AxSolverMod .EQ. lP3SENM) CALL CopySP3AxFlx(Flx0(iz), FlxOMP(iz, tid), AxGeom%lTransient)
ENDDO
ENDSUBROUTINE

SUBROUTINE MpiAxialSolver_Mod1(Core, CmInfo, Eigv, ng0, AxSolverMod, lreset, lDtilGen, PE)
USE PARAM
USE TYPEDEF,     ONLY : CMInfo_Type, CoreInfo_Type, PE_TYPE
USE P1SENM_MOD,  ONLY : P1SENM
USE Sp3Senm_Mod, ONLY : DirectSP3SENM
USE TIMER,       ONLY : nTracer_dclock, TimeChk
USE CMFDComm_mod, ONLY : CommPinXS_Type
USE BasicOperation, ONLY : CP_VA
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(CMInfo_Type) :: CMInfo
!TYPE(AxFlx_Type), POINTER :: Flx1D(:)
TYPE(PE_TYPE) :: PE
REAL :: eigv, peigv
INTEGER :: ng0, AxSolverMod
REAL :: AxNTimeBeg, AxNTimeEnd, ElaspedTime
LOGICAL :: lreset, lDtilGen

TYPE(CommPinXS_Type) :: CommPinXS

INTEGER :: iz, iz0, ixy
INTEGER,SAVE :: iter

LOGICAL, SAVE :: lFirst
DATA lfirst /TRUE/
DATA iter /0/
AxNTimeBeg = nTracer_dclock(FALSE, FALSE)
IF(lfirst .OR. lReset) THEN
  !CALL SetMpiAxSolverEnvironment(CORE, ng0, PE)
  !Communication module initialization
  !CALL InitAxNodalComm(GroupInfo, ng, 4, SendRecvMax, PE, .TRUE.)
  lfirst = FALSE
ENDIF

END SUBROUTINE


SUBROUTINE PhiAvgUpdate(PhiC, Flx0)
USE PARAM
USE TYPEDEF, ONLY : CMInfo_Type, AxFlx_Type
IMPLICIT NONE


REAL, POINTER :: PhiC(:, :, :)
TYPE(AxFlx_Type), POINTER :: Flx0(:, :)

INTEGER :: ixy, iz, ig
DO iz = 1, nzfm
  DO ixy = myNxybeg, mynxyend
    !DO ig = 1, ng
      Flx0(iz, ixy)%Phi(0, 1, 1:ng) = PhiC(ixy, iz, 1:ng)
    !ENDDO
  ENDDO
ENDDO
END SUBROUTINE


SUBROUTINE TlkgAvgUpdate(Tlkg, Flx0)
USE PARAM
USE TYPEDEF, ONLY : CMInfo_Type
IMPLICIT NONE

!TYPE(CMInfo_Type) :: CmInfo
!INTEGER :: nxy, iz1, iz2, ig1, ig2, ixy1, ixy2

REAL, POINTER :: Tlkg(:, :, :)
TYPE(AxFlx_TYPE), POINTER :: Flx0(:, :)
INTEGER :: ixy, iz, ig

DO iz = 1, nzfm
  DO ixy = myNxyBeg, myNxyEnd
    Flx0(iz, ixy)%Tlkg(0, 1, 1:ng) = Tlkg(ixy, iz, 1:ng)
  ENDDO
ENDDO
END SUBROUTINE

#endif
!
SUBROUTINE AxFluxOut(Flx0, AxGeom, nz)
USE PARAM
USE TYPEDEF,    ONLY : AxFlx_Type
TYPE(AxFlx_Type) :: Flx0(nz)
TYPE(AxGeom_Type) :: AxGeom
INTEGER :: nz

INTEGER :: iz 
WRITE(102, '(I5, 100F10.5)') AxGeom%H(1:nz)
DO iz = 1, nz 
  WRITE(102, '(100F10.5)') Flx0(iz)%psi(0:4) 
ENDDO
WRITE(102, *)
END SUBROUTINE


SUBROUTINE AxSrcUpdate(Core, CmInfo, iz1, iz2, ig1, ig2, PE, AxSolver)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type,   CmInfo_Type,      PE_TYPE,     &
                    Pin_Type,        AxFlx_Type,       PinXs_Type
USE CNTL,    ONLY : nTracerCntl
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
Type(CmInfo_Type) :: CmInfo
TYPE(PE_TYPE) :: PE
INTEGER :: AxSolver

TYPE(PinXS_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: PhiFm(:, :, :), PhiC(:, :, :)
REAL, POINTER :: hz(:), hzInv(:)
REAL, POINTER :: AxSrc(:, :, :), AxPXS(:, :, :)
REAL, POINTER :: AxDtil(:, :, :, :), AxDhat(:, :, :, :)
REAL :: LkgAvg(ig1:ig2)
REAL :: neighphi, myphi, jnet, Dtil, Dhat
INTEGER :: iz1, iz2, ig1, ig2
INTEGER :: i, j, k, ig
INTEGER :: ixy, iz, npin

PinXS => CmInfo%PinXS
AxSrc => CmInfo%AxSrc
AxPXS => CmInfo%AxPXS
PhiFm => CmInfo%PhiFM 
PhiC => CmInfo%PhiC
hz => Core%hz
HzInv => Core%HzInv
nPin = Core%nxy
!
AxDhat => CmInfo%AxDhat
AxDtil => CmInfo%AxDTil
!IF(PE%myrank .eq. 1) THEN

DO ixy = 1, nPin
  DO iz = iz1, iz2
    DO ig = ig1, ig2
      AxSrc(ixy, iz, ig) = 0; AxPxs(ixy, iz, ig) = 0;
      ! leakage source term from lower plane
      i = Core%SubPlaneRange(1, iz); j = i - 1
      myPhi = PhiFm(ixy, i, ig); neighPhi = PhiFm(ixy, j, ig) 
      Dtil = AxDtil(1, ixy, i, ig); Dhat =AxDhat(1, ixy, i, ig)
      Jnet = (Dtil - Dhat) * MyPhi - (Dtil + Dhat) * NeighPhi
      LkgAvg(ig) = Jnet
      
      ! leakage source term from upper plane
      i = Core%SubPlaneRange(2, iz); j = i + 1
      myPhi = PhiFm(ixy, i, ig); neighPhi = PhiFm(ixy, j, ig)
      Dtil = AxDtil(2, ixy, i, ig); Dhat =AxDhat(2, ixy, i, ig)
      Jnet = (Dtil - Dhat) * MyPhi - (Dtil + Dhat) * NeighPhi
      LkgAvg(ig) = LkgAvg(ig) + Jnet

      LkgAvg(ig) = LkgAvg(ig) * HzInv(iz)
    ENDDO
    DO ig = ig1, ig2
      AxSrc(ixy, iz, ig) = LkgAvg(ig)
      AxPXS(ixy, iz, ig) = 0
      SELECT CASE(nTracerCntl%LkgSplitLv)
      CASE(0)
          IF(LkgAvg(ig) .LT. 0)THEN
            !WRITE(*,*) 'NegLeakage'
            CYCLE
        ENDIF      
        IF( PHIC(ixy, iz, ig) .GT. 0 )THEN
            AxPXS(ixy, iz, ig) = LkgAvg(ig) / PhiC(ixy, iz, ig)
            !WRITE(*,*) 'NegLeakage at', iz, ixy, ig
        ENDIF      
        !AxSrc(ixy, iz, ig) =  0 
      CASE(1)
        !--- Why was it turned on?
!        IF( PHIC(ixy, iz, ig) .NE. 0 )THEN
!            AxPXS(ixy, iz, ig) = LkgAvg(ig) / PhiC(ixy, iz, ig)
!        ENDIF
      ENDSELECT      
    ENDDO
  ENDDO
ENDDO
!ENDIF
NULLIFY(PinXS)
NULLIFY(AxSrc)
NULLIFY(AxPXS)
NULLIFY(PhiC)
NULLIFY(hz)
NULLIFY(HzInv)
END SUBROUTINE

SUBROUTINE VoidLkgCorrection(Core, FmInfo, CmInfo, iz1, iz2, ig1, ig2, PE)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type,   FmInfo_Type,      CmInfo_Type,      PE_TYPE,     &
                    Pin_Type,        Cell_Type,        FxrInfo_Type
USE IOUTIL,         ONLY : message
USE FILES,          ONLY : io8
USE MPIComm_Mod,    ONLY : Reduce
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(PE_TYPE) :: PE

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
REAL, POINTER :: AxSrc(:, :, :), AxPXS(:, :, :)
INTEGER :: iz1, iz2, ig1, ig2
INTEGER :: i, j, k, ig, icel, ifxr
INTEGER :: ixy, iz, npin, niso
INTEGER :: FxrIdxSt, nLocalFxr
INTEGER :: ivoidWarning, ivoidWarning0
REAL :: vol, voidvol, fvoid, nd, voidcrit

LOGICAL, SAVE :: lFirst = .TRUE.
LOGICAL :: lVoidWarning
LOGICAL :: lvoidCorrect

AxSrc => CmInfo%AxSrc; AxPXS => CmInfo%AxPXS
Fxr => FmInfo%Fxr; Pin => Core%Pin; CellInfo => Core%CellInfo
nPin = Core%nxy
lVoidWarning = .FALSE.
lvoidCorrect = .FALSE.; iVoidWarning = 0
voidcrit = 5.E-4_8
DO iz = iz1, iz2
  DO ixy = 1, nPin
    !Calculate Void Fraction
    icel = Pin(ixy)%cell(iz)
    FxrIdxSt = Pin(ixy)%FxrIdxSt; nLocalFxr = CellInfo(icel)%nFxr
    voidvol = 0; vol = 0
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      vol = vol + Fxr(ifxr, iz)%area
      niso = Fxr(ifxr, iz)%niso
      nd = sum(Fxr(ifxr,iz)%pnum(1:niso))
      IF(nd .LT. voidcrit ) THEN
        voidvol =voidvol + Fxr(ifxr, iz)%area
        Fxr(ifxr, iz)%lVoid = .TRUE.
      ELSE
        Fxr(ifxr, iz)%lVoid = .FALSE.
      ENDIF
    ENDDO
    fvoid = voidvol / vol
    IF(fvoid .GT. 1.e-5_8) lVoidCorrect = .TRUE.
    IF(fvoid .GT. 0.9) iVoidWarning = 1
    
    fvoid = 1._8 - fvoid
    IF(fvoid .LT. 0.05_8) fvoid = 1
    fvoid = 1._8/fvoid
    
    DO ig = ig1, ig2
       AxPXS(ixy, iz, ig) =  AxPXS(ixy, iz, ig) * fvoid
       AxSrc(ixy, iz, ig) =  AxSrc(ixy, iz, ig) * fvoid
    ENDDO
  ENDDO
ENDDO

IF(lFirst) THEN
#ifdef MPI_ENV  
  CALL REDUCE(iVoidWarning, ivoidWarning0, PE%MPI_CMFD_COMM, .TRUE.)
#else
  ivoidWarning0 = iVoidWarnig
#endif
  IF(iVoidWarning0 .GT. 0) lVoidWarning = .TRUE.
  IF(lVoidWarning .AND. PE%CMFDMaster) THEN
    WRITE(mesg,'(a)') 'Warning : There exist large void region in the core'
    CALL message(io8, TRUE, TRUE, mesg)     
  ENDIF
  lFirst = .FALSE.
ENDIF
NULLIFY(AxSrc); NULLIFY(AxPXS)
NULLIFY(Fxr); NULLIFY(CellInfo)
!AxSrc => CmInfo%AxSrc; AxPXS => CmInfo%AxPXS
!Fxr => FmInfo%Fxr; Pin => Core%Pin; CellInfo => Core%CellInfo
END SUBROUTINE

SUBROUTINE AxialXsMapping(XsTarget, XsDestination)
USE PARAM
USE TYPEDEF, ONLY : PInXS_TYPE
IMPLICIT NONE
TYPE(PinXs_TYPE) :: XsTarget, XsDestination

XsDestination%Dtil => XsTarget%Dtil;  XsDestination%Dhat => XsTarget%Dhat
XsDestination%XSD => XsTarget%XSD;    XsDestination%XSD2 => XsTarget%XSD2
XsDestination%XST => XsTarget%XST;    XsDestination%XSTR => XsTarget%XSTR
XsDestination%XSR => XsTarget%XSR;    XsDestination%XSNF => XsTarget%XSNF
XsDestination%XSKF => XsTarget%XSKF;  XsDestination%CHI => XsTarget%CHI
XsDestination%Xss => XsTarget%Xss
END SUBROUTINE

SUBROUTINE FinalizeAxialXsMapping(XsDestination, nz)
USE PARAM
USE TYPEDEF, ONLY : PInXS_TYPE
IMPLICIT NONE
TYPE(PinXs_TYPE) :: XsDestination(nz)
INTEGER :: nz
INTEGER :: iz

DO iz = 1, nz
  NULLIFY(XsDestination(iz)%Dtil); NULLIFY(XsDestination(iz)%Dhat)
  NULLIFY(XsDestination(iz)%XSD); NULLIFY(XsDestination(iz)%XSD2)
  NULLIFY(XsDestination(iz)%XST); NULLIFY(XsDestination(iz)%XSTR)
  NULLIFY(XsDestination(iz)%XSR); NULLIFY(XsDestination(iz)%XSNF)
  NULLIFY(XsDestination(iz)%XSKF); NULLIFY(XsDestination(iz)%CHI)
  NULLIFY(XsDestination(iz)%Xss)
ENDDO
END SUBROUTINE

SUBROUTINE CalAxDtil(Flx0, XS, Geom, nz, AxSolver)
USE PARAM
USE TYPEDEF,  ONLY : PinXs_Type, AxFlx_type
IMPLICIT NONE
TYPE(AxFlx_Type) :: Flx0(nz)
TYPE(PinXs_Type) :: Xs(nz)
TYPE(AxGEOM_TYPE) :: Geom
INTEGER :: nz
INTEGER :: AxSolver

INTEGER :: mp(1:2)
INTEGER :: i
INTEGER :: iz, ineigh,ig 
REAL :: Dtil, mybeta, neighbeta
mp = (/-1,1/)
DO iz = 1, nz
  DO i = 1, 2
    ineigh = iz + mp(i)
    neighbeta = 0
    IF(ineigh .EQ. 0 .AND. Geom%BC(i) .EQ. VoidCell) neighbeta = 0.25 
    IF(ineigh .EQ. nz+1 .AND. Geom%BC(i) .EQ. VoidCell) neighbeta = 0.25
    DO ig = 1, Geom%ng
      mybeta = XS(iz)%XSD(ig) / Geom%h(iz)
      IF(ineigh .NE. 0 .AND. ineigh .NE. nz+1) THEN
        neighbeta = XS(ineigh)%XSD(ig) / Geom%h(ineigh)
      ENDIF
      Dtil = 2._8 * neighbeta * mybeta / (mybeta + neighbeta)
      Flx0(iz)%Dtil(i, ig) = Dtil
      Flx0(iz)%Dhat(i, ig) = 0    
      Flx0(iz)%Phi(1:4, :, :) = 0
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE NormMpiAxNVar(norm, AxSolMod)
!Normalize Axail Nodal Method Variables
USE PARAM
USE TYPEDEF,     ONLY : AxFlx_Type
USE BasicOperation, ONLY : MULTI_CA
IMPLICIT NONE
REAL :: norm
INTEGER :: AxSolMod
INTEGER :: ixy, iz

INTEGER :: nod 
nod = 1
IF(AxSolMod .EQ. lP3SENM) nod = 2
DO ixy = myNxyBeg, myNxyEnd
  DO iz = 1, nzfm
    CALL MULTI_CA(norm, AxFlx(iz ,ixy)%PHI(0:4, 1:nod, 1:ng), 5, nod, ng)
    CALL MULTI_CA(norm, AxFlx(iz ,ixy)%PSI(0:4), 5)
    CALL MULTI_CA(norm, AxFlx(iz, ixy)%JOUT(1:nod, 1:2, 1:ng), 2, nod, ng)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE UpdtTranSol2AxNvar(CmInfo, PinVol, eigv, PE)
USE PARAM
USE TYPEDEF,          ONLY : CmInfo_Type,   PE_TYPE
USE BasicOperation,   ONLY : CP_VA
USE CMFDComm_mod,     ONLY : CommAxNodalInPut
IMPLICIT NONE

TYPE(CmInfo_Type) :: CMInfo
REAL, POINTER :: PinVol(:, :)
REAL :: eigv
TYPE(PE_TYPE) :: PE
REAL, POINTER :: PsiC(:, :), PhiC(: ,:, :)
INTEGER :: ixy, iz, ig

PhiC => CmInfo%PhiFm; PsiC => CmInfo%PsiFm
DO iz = myzbf, myzef
  DO ixy = 1, nxy
    LKG(ixy ,iz, 1) = eigv * PsiC(ixy, iz) / PinVol(ixy, iz)
  ENDDO
ENDDO
IF(PE%nCmfdProc .GT. 1) THEN
  CALL CommAxNodalInPut(CmInfo%PhiFm, Lkg, MpiPhi, MpiLkg, ng, PE)
ELSE
  MpiLkg => Lkg;   MpiPhi => PhiC
ENDIF

DO ixy = mynxybeg, mynxyend
  DO iz = 1, nzfm
    DO ig = 1, ng
      AxFlx(iz, ixy)%Phi(0, 1, ig) = MpiPhi(ixy, iz, ig)
    ENDDO
    AxFlx(iz, ixy)%Psi(0) = MpiLkg(ixy, iz, 1)
  ENDDO
ENDDO

NULLIFY(PsiC, PhiC)

END SUBROUTINE

SUBROUTINE SetPrecSrc(PinXS, TranInfo, TranCntl, PE)
USE PARAM
USE TYPEDEF,     ONLY : PinXS_Type,     TranInfo_Type,      TranCntl_Type,     PE_Type
USE TranAxNUtil_Mod, ONLY : UpdtAxNodePrecSrc
IMPLICIT NONE
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(PE_TYPE) :: PE

REAL, POINTER :: Chid(:)

INTEGER :: nprec, NowStep
INTEGER :: ixy, iz, iz0, ig, iprec

REAL :: lambda(maxPrec), kappa(MaxPrec), omegam, omega0
REAL :: PrecSrc(0:4)

Chid => TranInfo%Chid

nprec = TranInfo%nprec
NowStep = TranCntl%NowStep
DO iprec = 1, TranInfo%nprec
  lambda(iprec) = TranInfo%Lambda(iprec)  
  kappa(iprec) = exp(-lambda(iprec) * TranCntl%DelT(NowStep))
ENDDO

DO ixy = mynxybeg, mynxyend
  DO iz = 1, nzfm
    iz0 = SubPlaneMap(iz)
    AxFlx(iz, ixy)%TranSrc = 0
    omegam = AxOmegam(0, ixy, iz0); omega0 = AxOmega0(0, ixy, iz0)
    PrecSrc = UpdtAxNodePrecSrc(AxFlx(iz, ixy), lambda(1:nprec), kappa(1:nprec), omegam, omega0)
    DO ig = 1, ng
      AxFlx(iz, ixy)%TranSrc(0:4, ig) = Chid(ig) * PrecSrc(0:4)
    ENDDO

  ENDDO
ENDDO
END SUBROUTINE



SUBROUTINE SetAxNKinParam(Core, PinXS, TranInfo, TranCntl, PE)
USE PARAM
USE TYPEDEF,      ONLY : CoreInfo_Type,    TranInfo_Type,     TranCntl_Type
USE CMFDComm_Mod, ONLY : CommKineticParam, CommAxNBeta
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(PinXs_Type), POINTER :: PiNXS(:, :)
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(PE_Type) :: PE

INTEGER :: iz, ixy

IF(PE%nCmfdProc .EQ. 1) THEN
  AxOmega0 => TranInfo%CellOmega0
  AxOmegap => TranInfo%CellOmegap
  AxOmegam => TranInfo%CellOmegam
ELSE
  CALL CommKineticParam(TranInfo%CellOmega0, AxOmega0, 0, TranInfo%nPrec, ng, SubPlaneMap, PE)
  CALL CommKineticParam(TranInfo%CellOmegap, AxOmegap, 0, TranInfo%nPrec, ng, SubPlaneMap, PE)
  CALL CommKineticParam(TranInfo%CellOmegam, AxOmegam, 0, TranInfo%nPrec, ng, SubPlaneMap, PE)
  CALL CommAxNBeta(AxBeta, PinXS, TranInfo%nprec, ng, SubPlaneMap, PE)
  DO ixy = mynxybeg, mynxyend
    DO iz = 1, nz
      AxPinXS(ixy,iz)%Betat = sum(AxBeta(1:TranInfo%nPrec, ixy, iz))
    ENDDO
  ENDDO
  
ENDIF


END SUBROUTINE


SUBROUTINE InitTranAxNSolver(AxSolver, Eigv, TranInfo, TranCntl, PE)
USE PARAM
USE TYPEDEF,               ONLY : TranInfo_Type,      TranCntl_Type,         PE_TYPE
USE TranAxNUtil_Mod,       ONLY : InitTranAxNUtil,    AllocTranAxNodeVar 
USE ALLOCS
IMPLICIT NONE

TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(PE_Type) :: PE
INTEGER :: AxSolver
REAL :: Eigv

INTEGER :: ixy, iz, itid

CALL InitTranAxNUtil(AxSolver, ng, TranInfo%nPrec, PE)
DO ixy = PE%mynxybeg, PE%mynxyend
  DO iz = 1, nzfm
    !CALL AllocTranAxNodeVar(AxFlx(iz, ixy))
    CALL AllocTranAxNodeVar(AxFlx(iz, ixy), ng, TranInfo%nPrec)
  ENDDO
ENDDO

TranSrc2nd => Lkg
IF(PE%nCmfdProc .GT. 1) THEN
  CALL Dmalloc0(AxOmegam, 0, TranInfo%nPrec, mynxybeg, mynxyend, 1, nz)
  CALL Dmalloc0(AxOmega0, 0, TranInfo%nPrec, mynxybeg, mynxyend, 1, nz)
  CALL Dmalloc0(AxOmegap, 0, TranInfo%nPrec, mynxybeg, mynxyend, 1, nz)
  CALL Dmalloc0(AxBeta, 1, TranInfo%nPrec, mynxybeg, mynxyend, 1, nz)
  CALL Dmalloc0(MpiTranSrc2nd, mynxybeg, mynxyend, 1, nz, 1, ng)
ELSE
  MpiTranSrc2nd => MpiLkg  
ENDIF

!CALL DMALLOC0()



CALL InitTranAxNVar(eigv)

DO itid = 1, PE%nAxThread
  DO iz = 1, nzfm
    CALL AllocTranAxNodeVar(FLXOMP(iz, itid), ng, TranInfo%nPrec)
  ENDDO
ENDDO
END SUBROUTINE



!SUBROUTINE AllocTranAxNSolver(ng0, nPrec0, PE)
!USE PARAM
!USE ALLOCS
!IMPLICIT NONE
!TYPE(PE_TYPE) :: PE
!INTEGER :: ng0, nprec0
!
!INTEGER :: ixy, iz
!DO ixy = PE%mynxybeg, PE%mynxyend
!  DO iz = 1, PE%nzfm
!    !CALL AllocTranAxNodeVar(AxFlx(iz, ixy))
!    CALL DMALLOC0(AxFlx(iz, ixy)%TranPsi, 0, 4)
!    CALL DMALLOC0(AxFlx(iz, ixy)%TranPsid, 0, 4)
!    CALL DMALLOC0(AxFlx(iz, ixy)%Prec, 0, 4, 1, nprec0)
!    CALL DMALLOC0(AxFlx(iz, ixy)%TranSrc, 0, 4, 1, ng0)
!  ENDDO
!ENDDO
!
!END SUBROUTINE

SUBROUTINE InitTranAxNVar(eigv)
USE PARAM
IMPLICIT NONE
INTEGER :: ixy, iz, ig
REAL :: eigv, reigv

reigv = 1._8 / eigv
DO ixy = mynxybeg, mynxyend
  DO iz = 1, nzfm
    AxFlx(iz, ixy)%TranPsi(0:4) = reigv * AxFlx(iz, ixy)%Psi(0:4)
    AxFlx(iz, ixy)%TranPsid(0:4) = reigv * AxFlx(iz, ixy)%Psi(0:4)
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE InitAxNPrec(PinXs, TranInfo, PE)
USE PARAM
USE TYPEDEF,              ONLY : PinXs_Type, TranInfo_Type
USE TranAxNUtil_Mod,      ONLY : InitAxNodePrecursor
USE CMFDComm_Mod,         ONLY : CommAxNBeta
IMPLICIT NONE
TYPE(PinXs_Type), POINTER :: PinXS(:, :)
TYPE(TranInfo_Type) :: TranInfo
TYPE(PE_TYPE) :: PE
REAL ::eigv
INTEGER :: ixy, iz, iz0, iprec

REAL :: reigv
REAL :: Invlambda(maxprec), InvLamBeta(maxprec)


eigv = TranInfo%eigv0
reigv = 1._8 / eigv

DO iprec = 1, TranInfo%nPrec
  InvLambda(iprec) = 1._8 / TranInfo%lambda(iprec)
ENDDO
IF(PE%nCMFDProc .GT. 1) THEN
  CALL CommAxNBeta(AxBeta, PinXS, TranInfo%nprec, ng, SubPlaneMap, PE)
ENDIF
DO ixy = mynxybeg, mynxyend
  DO iz = 1, nzfm
    iz0 = AxGeom%Comp(iz)
    IF(PE%nCmfdProc .GT. 1) THEN
      DO iprec = 1, TranInfo%nprec
        InvLamBeta(iprec) = AxBeta(iprec, ixy, iz) * InvLambda(iprec)
      ENDDO      
    ELSE
      DO iprec = 1, TranInfo%nprec
        InvLamBeta(iprec) = PinXS(ixy, iz0)%beta(iprec) * InvLambda(iprec)
      ENDDO
    ENDIF
    CALL InitAxNodePrecursor(AxFlx(iz, ixy), InvLamBeta, reigv)
  
  ENDDO
ENDDO
END SUBROUTINE


SUBROUTINE SetTransientChi(PinXS, TranInfo, TranCntl, lreset, PE)
USE PARAM
USE TYPEDEF,                ONLY : PinXS_TYPE,      TranInfo_TYPE,          TranCntl_TYPE
IMPLICIT NONE
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(PE_Type) :: PE
REAL :: eigv
LOGICAL :: lreset

REAL, POINTER :: Chid(:)
REAL :: ChiEff
INTEGER :: iz, ixy, ig

Chid => TranInfo%ChiD
IF(.NOT. lreset) THEN
  IF(PE%nCmfdProc .GT. 1) THEN
    DO iz = 1, nz
      DO ixy = mynxybeg, mynxyend
        DO ig = 1, ng
          ChiEff = Chid(ig) * (AxOmegap(0,ixy,iz) - PinXs(ixy, iz)%betat) + PinXS(ixy, iz)%Chi(ig)
          PinXS(ixy, iz)%Chi(ig) = ChiEff
        ENDDO
      ENDDO
    ENDDO
  ELSE
    DO iz = 1, nz
      DO ixy = mynxybeg, mynxyend
        DO ig = 1, ng
          ChiEff = Chid(ig) * (PinXs(ixy, iz)%omega - PinXs(ixy, iz)%betat) + PinXS(ixy, iz)%Chi(ig)
          PinXS(ixy, iz)%Chi(ig) = ChiEff
        ENDDO
      ENDDO
    ENDDO
  ENDIF
ELSE
  IF(PE%nCmfdProc .GT. 1) THEN
    DO iz = 1, nz
      DO ixy = mynxybeg, mynxyend
        DO ig = 1, ng
          ChiEff = PinXS(ixy, iz)%Chi(ig) - Chid(ig) * (AxOmegap(0,ixy,iz) - PinXs(ixy, iz)%betat)
          PinXS(ixy, iz)%Chi(ig) = ChiEff
        ENDDO
      ENDDO
    ENDDO
  ELSE
    DO iz = 1, nz
      DO ixy = mynxybeg, mynxyend
        DO ig = 1, ng
          ChiEff = PinXS(ixy, iz)%Chi(ig) - Chid(ig) * (PinXs(ixy, iz)%omega - PinXs(ixy, iz)%betat)
          PinXS(ixy, iz)%Chi(ig) = ChiEff
        ENDDO
      ENDDO
    ENDDO
  ENDIF
ENDIF

NULLIFY(ChiD)
END SUBROUTINE

SUBROUTINE UpdtAxNTranSrc2d(Core, CmInfo, TranInfo, TranCntl, PE)
!Update 2nd order transient source term component
USE PARAM
USE TYPEDEF,                 ONLY : CoreInfo_Type,       CmInfo_Type,       TranInfo_Type,       &
                                    TranCntl_Type,       PE_TYPE,                                &
                                    PinXs_Type
USE CMFDComm_Mod,            ONLY : CommTranSrc2nd
USE TranAxNUtil_Mod,         ONLY : UpdtAxNodeTranSrc2d
USE Tlkg_mod,                ONLY : Convtlkg2nd
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(PE_TYPE) :: PE

TYPE(PinXs_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: PinVolFm(:, :)
REAL, POINTER :: PhiC(:, :, :), TranPhiC(:, :, :), ResSrc(:, :, :)
REAL, POINTER :: Expo(:, :, :), Expo_Alpha(:, :, :) 


INTEGER :: NowStep
REAL :: DelT, Theta, ThetaH, rvdt, AlphaRV
REAL :: PrevSrc, ExpSrc

INTEGER :: ig, ixy, iz, iz0
PinXs => CmInfo%PinXs; PinVolFm => Core%PinVolFm
PhiC => CmInfo%PhiFm; TranPhiC => CmInfo%TranPhiFm
ResSrc => CmInfo%ResSrcFm
Expo => TranInfo%Expo; Expo_Alpha => TranInfo%Expo_Alpha
NowStep = TranCntl%NowStep
DelT = TranCntl%DelT(NowStep); Theta = TranCntl%Theta
thetah = 1._8 / theta - 1
DO ig = 1, ng
  DO iz = myzbf, myzef
    iz0 = SubPlaneMap(iz)
    DO ixy = 1, nxy
      rvdt = 1._8 / (Delt * PinXs(ixy, iz0)%velo(ig) * Theta)
      PrevSrc = Expo(ixy, iz0, ig)  * ThetaH * ResSrc(ixy, iz, ig) / PinVolFm(ixy, iz0) 
      PrevSrc = PrevSrc - rvdt * (PhiC(ixy, iz, ig) - Expo(ixy, iz0, ig) * TranPhiC(ixy, iz, ig))
      AlphaRV =  Expo_Alpha(ixy, iz0, ig) / PinXs(ixy, iz0)%velo(ig)
      ExpSrc = - AlphaRv * (PhiC(ixy, iz, ig)+Thetah * Expo(ixy, iz0, ig) * TranPhiC(ixy, iz, ig))
      TranSrc2nd(ixy, iz, ig) = PrevSrc + ExpSrc
      
      !Theta Method
      !TranSrc2nd(ixy, iz, ig) = Thetah * ResSrc(ixy, iz, ig) / PinVolFm(ixy, iz0)
      !rvdt = 1._8 / (Delt * PinXs(ixy, iz0)%velo(ig) * Theta)
      !TranSrc2nd(ixy, iz, ig) = TranSrc2nd(ixy, iz, ig)  - rvdt *(PhiC(ixy, iz, ig)-TranPhiC(ixy, iz, ig))
    ENDDO
  ENDDO
ENDDO

IF(PE%nCmfdProc .GT. 1) THEN
  CALL CommTranSrc2nd(TranSrc2nd, MpiTranSrc2nd, ng, PE)
ENDIF

DO ixy = mynxybeg, mynxyend
  CALL UpdtAxNodeTranSrc2d(AxFlx(1:nzfm, ixy), MpiTranSrc2nd(ixy, 1:nzfm, 1:ng),  AxGeom%BC)
ENDDO
END SUBROUTINE

SUBROUTINE UpdtAxNPrecSrc(TranInfo, TranCntl, PE)
!Update 2nd order transient source term component
USE PARAM
USE TYPEDEF,                 ONLY : CoreInfo_Type,       CmInfo_Type,       TranInfo_Type,       &
                                    TranCntl_Type,       PE_TYPE,                                &
                                    PinXs_Type
USE TranAxNUtil_Mod,         ONLY : UpdtAxNodePrecursor
USE Tlkg_mod,                ONLY : Convtlkg2nd
IMPLICIT NONE
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(PE_TYPE) :: PE

TYPE(PinXs_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: PinVolFm(:, :)
REAL, POINTER :: Omegam(:, :, :), Omegap(:, :, :), Omega0(:, :, :)


REAL :: kappa(maxprec), lambda(maxprec)
REAL :: eigv, reigv, delT
INTEGER :: nprec, nowstep
INTEGER :: ixy, iz, ig, iprec

NowStep = TranCntl%NowStep
Delt = TranCntl%DelT(NowStep)

Omegam => TranInfo%CellOmegam
Omegap => TranInfo%CellOmegap
Omega0 => TranInfo%CellOmega0
eigv = 1._8 / eigv

DO iprec = 1, TranInfo%nPrec
  lambda(iprec) = TranInfo%lambda(iprec)
  kappa(iprec) = exp(-delt * lambda(iprec))    
ENDDO

IF(PE%nCmfdProc .EQ. 1) THEN
  AxOmegam => Omegam
  AxOmegap => Omegap
  AxOmega0 => Omega0
ENDIF

DO ixy = mynxybeg, mynxyend
  DO iz = 1, nzfm
    CALL UpdtAxNodePrecursor(AxFlx(iz, ixy), kappa(1:nprec), AxOmegam(:, ixy, iz),       &
                             AxOmega0(:, ixy, iz), AxOmegap(:, ixy, iz), reigv)
  ENDDO
ENDDO


END SUBROUTINE

SUBROUTINE SaveTranAxNSol(eigv, PE)
USE PARAM
USE TYPEDEF,            ONLY :  PE_TYPE
USE TranAxNUtil_Mod,    ONLY : SaveTranAxNodeSol
IMPLICIT NONE
REAL :: eigv
TYPE(PE_TYPE) :: PE

REAL :: reigv
INTEGER :: iz, ixy

reigv = 1._8 / eigv
DO ixy = myNxyBeg, myNxyEnd
  DO iz = 1, nzfm
    CALL SaveTranAxNodeSol(AxFlx(iz, ixy), reigv)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE MpiTranAxialSolver(Core, CmInfo, TranInfo, GroupInfo, TranCntl, AxNCntl, PE)
USE PARAM
USE TYPEDEF,          ONLY : CMInfo_Type,       CoreInfo_Type,          PE_TYPE,         &
                             TranInfo_Type,     TranCntl_Type,          GroupInfo_Type
!USE AxSolver_Mod,     ONLY : AxialXsMapping,    FinalizeAxialXsMapping, CalAxDtil
USE Tlkg_mod,         ONLY : RadTlkgUpdate,     SetTlkgShape,   SetTlkgShape0
USE P1SENM_MOD,       ONLY : SetTranP1SenmEnv
USE Sp3Senm_Mod,      ONLY : SetSP3SenmEnv
USE CMFDComm_mod,     ONLY : CommPinXS_Type,                                             &
                             CommGlobalPinXsDat,   CommAxNodalInPut,                     &
                             CommAxNodalOutPut
USE TIMER,            ONLY : nTracer_dclock, TimeChk
USE BasicOperation,   ONLY : CP_VA
USE OMP_LIB
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CMInfo_Type) :: CMInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(AxNCntl_Type) :: AxNCntl
TYPE(PE_TYPE) :: PE


REAL :: eigv, peigv
INTEGER :: ng0, AxSolverMod
REAL :: AxNTimeBeg, AxNTimeEnd, ElaspedTime
LOGICAL :: lXsComm, lDtilGen, lReset

TYPE(CommPinXS_Type) :: CommPinXS

INTEGER :: iz, iz0, ixy, ig

INTEGER,SAVE :: iter

INTEGER :: lkgmod, tid

LOGICAL, SAVE :: lFirst
DATA lfirst /TRUE/
DATA iter /0/

lreset = AxNCntl%lreset; lDtilGen = AxNCntl%lDtilUpdt
AxSolverMod = AxNCntl%AxSolverMod; lXsComm = AxNCntl%lXsComm
eigv = AxNCntl%eigv; ng0 = AxNCntl%ng0

AxGeom%lTransient = .TRUE.

iter = iter + 1

IF(PE%CMFDMaster) AxNTimeBeg = nTracer_dclock(FALSE, FALSE)

IF(lfirst) THEN
  IF(AxSolverMod .EQ. lP1SENM) CALL SetTranP1SenmEnv(TranCntl)
  lfirst = FALSE
ENDIF

myNxyBeg = PE%myNxyBeg; myNxyEnd = PE%myNxyEnd

!Determine Transverse Leakage Information
CALL RadTlkgUpdate(Core, CmInfo, Lkg, 1, nxy, myzbf, myzef, ng)


IF(PE%nCmfdProc .GT. 1) THEN
  CommPinXS%SendPhi => CmInfo%PhiFm; CommPinXS%SendLkg => Lkg
  CommPinXS%SendPinXs => CmInfo%PinXS
  CommPinXS%RecvPhi => MpiPhi; CommPinXs%RecvLkg => MpiLkg
  CommPinXS%RecvPinXs => AxPinXS        !
  IF(lXsComm) CALL CommGlobalPinXsDat(CommPinXS, ng, GroupInfo, PE)
  !Destroy the CommPinXS variables
  NULLIFY(CommPinXS%SendPhi, CommPinXS%SendLkg, CommPinXS%SendPinXS)
  NULLIFY(CommPinXS%RecvPhi, CommPinXS%RecvLkg, CommPinXS%RecvPinXS)
  !Get Lkg and flux for Radial Domain
  CALL CommAxNodalInPut(CmInfo%PhiFm, Lkg, MpiPhi, MpiLkg, ng, PE)
ELSE
  MpiLkg => Lkg
  MpiPhi => CmInfo%PhiFm 
  AxPinXs => CmInfo%PinXS
ENDIF
IF(lXsComm) CALL SetAxNKinParam(Core, CMInfo%PinXs, TranInfo, TranCntl, PE)
!Update Average Flux info to the Axial Nodal structure
CALL PhiAvgUpdate(MpiPhi, AxFlx)
CALL TlkgAvgUpdate(Mpilkg, AxFlx)


CALL SetPrecSrc(AxPinXS, TranInfo, TranCntl, PE)
CALL UpdtAxNTranSrc2d(Core, CmInfo, TranInfo, TranCntl, PE)
CALL SetTransientChi(AxPinXS, TranInfo, TranCntl, .FALSE., PE)
lkgmod = 0
tid = 1

!print *, tid
!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(PE%nAxThread) 
!$OMP PARALLEL DEFAULT(SHARED)  &
!$OMP PRIVATE(tid, ixy)
!$  tid = omp_get_thread_num()+1
!$OMP DO ORDERED SCHEDULE(DYNAMIC)
DO ixy = myNxyBeg, myNxyEnd
  CALL MpiNodal1D_OMP(AxFlx(:, ixy), AxPinXS(ixy, :), Eigv, AxGeom, AxSolverMod, lDtilGen, tid, PE)
ENDDO
!$OMP END DO
!$OMP END PARALLEL
CALL SetTransientChi(AxPinXS, TranInfo, TranCntl, .TRUE., PE)
CALL MPI_SYNC(PE%MPI_nTRACER_COMM)

!Axial Nodal Calculation Output Distributions
IF(PE%nCmfdProc .GT. 1) THEN
  CALL CommAxNodalOutPut(CmInfo, AxFlx, ng, PE)
ELSE
  DO ixy = 1, nxy
    DO iz = myzbf, myzef
      DO ig = 1, ng
        CmInfo%AxDtil(1:2, ixy, iz, ig) = AxFlx(iz, ixy)%Dtil(1:2, ig)
        CmInfo%AxDhat(1:2, ixy, iz, ig) = AxFlx(iz, ixy)%Dhat(1:2, ig)
        AxFlx(iz, ixy)%PDhat(1:2, ig) = AxFlx(iz, ixy)%Dhat(1:2, ig)
      ENDDO
    ENDDO
  ENDDO
ENDIF

IF(PE%CMFDMaster) THEN
  AxNTimeEnd = nTracer_dclock(FALSE, FALSE)
  ElaspedTime = AxNTimeEnd - AxNTimeBeg
  TimeChk%AxialNodalTime = TimeChk%AxialNodalTime + ElaspedTime
ENDIF
CALL MPI_SYNC(PE%MPI_CMFD_COMM)

END SUBROUTINE

SUBROUTINE Get1dFlx4th(Flx4th, ixy, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,          ONLY : PE_TYPE
USE CNTL,             ONLY : nTracerCntl_Type
USE BasicOperation,   ONLY : CP_CA,            CP_VA
#ifdef MPI_ENV
USE MPIComm_Mod,      ONLY : REDUCE
#endif
IMPLICIT NONE

INTEGER :: ixy
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
REAL :: FLX4TH(0:4, nzfm, ng), Buf(0:4, nzfm, ng)
INTEGER :: iz, ig
INTEGER :: COMM

COMM = PE%MPI_CMFD_COMM
CALL CP_CA(Flx4TH(0:4, 1:nzfm, 1:ng), 0._8, 5, nzfm, ng)

IF(ixy .GE. myNxyBeg .AND. ixy .LE. myNxyEnd) THEN
  DO iz = 1, nzfm
    DO ig = 1, ng
      Flx4TH(0:4, iz ,ig) = AxFlx(iz, ixy)%phi(0:4, 1, ig)
    ENDDO
  ENDDO
ENDIF
#ifdef MPI_ENV
CALL REDUCE(Flx4th(0:4, 1:nzfm, 1:ng), Buf(0:4, 1:nzfm, 1:ng), 5, nzfm, ng, COMM, .TRUE.)
CALL CP_VA(Flx4th(0:4, 1:nzfm, 1:ng), Buf(0:4, 1:nzfm, 1:ng), 5, nzfm, ng)
#endif
END SUBROUTINE 

SUBROUTINE StabGapPinLkg(Core, Flx0,  PE)
USE TYPEDEF,    ONLY :  CoreInfo_Type,   AxFlx_Type,       Pin_Type,         &
                        PE_Type,         cell_type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(AxFlx_Type), POINTER :: Flx0(:, :)
TYPE(PE_TYPE) :: PE

TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Pin_TYpe), POINTER :: Pin(:)
INTEGER :: ixy, iz, icel, ipin

Pin => Core%Pin
CellInfo => Core%CellInfo
ipin=0
DO ixy = mynxybeg, mynxyend
 ! IF(.NOT. PIN(ixy)%lGap) CYCLE
  DO iz = 1, nz  
    icel = Pin(ixy)%Cell(iz)
    IF(.NOT. CellInfo(icel)%lgap) CYCLE
    ipin=ipin+1;
    Flx0(iz, ixy)%Dhat = 0
    !Flx0(iz, ixy)%Dtil = 0
    Flx0(iz, ixy)%PDhat = 0
  ENDDO
ENDDO
IF(ipin.NE.0) WRITE(*,*) ipin, '/',mynxyend-mynxybeg+1, 'GapPins Ax Dhat = 0'


END SUBROUTINE

SUBROUTINE StabCrPinLkg(Core, Flx0,  PE)
USE TYPEDEF,    ONLY :  CoreInfo_Type,   AxFlx_Type,       Pin_Type,         &
                        PE_Type,         cell_type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(AxFlx_Type), POINTER :: Flx0(:, :)
TYPE(PE_TYPE) :: PE

TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Pin_TYpe), POINTER :: Pin(:)
INTEGER :: ixy, iz, icel, ipin, ncel

Pin => Core%Pin
CellInfo => Core%CellInfo
ipin=0;ncel=0;
DO ixy = mynxybeg, mynxyend
    ipin=ipin+1;
    DO iz = 1, nz  
        icel = Pin(ixy)%Cell(iz)
        ncel=ncel+1;
        Flx0(iz, ixy)%Dhat = 0
        Flx0(iz, ixy)%PDhat = 0
    ENDDO
ENDDO
!WRITE(*,*) ipin, '/',mynxyend-mynxybeg+1, 'CrPins Ax Dhat = 0'
IF( ncel .NE. 0 ) WRITE(*,*) ncel, '/',(mynxyend-mynxybeg+1)*nz, 'CrCells Ax Dhat = 0'
END SUBROUTINE
END MODULE
