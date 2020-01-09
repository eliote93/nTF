#include <defines.h>
MODULE MPIConfig_Mod
USE PARAM
USE TypeDef,      ONLY : PE_TYPE
USE MPICOMM_MOD,  ONLY : MPIWaitTurn, REDUCE, BCAST, MPI_SYNC
USE IOUTIL,       ONLY : Terminate
USE UtilFunction, ONLY : GetRangeDecomp
IMPLICIT NONE

#ifdef MPI_ENV
INCLUDE 'mpif.h'
#endif

CONTAINS

SUBROUTINE SetMPIEnv(PE)
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
!Communicator Setting
CALL SetMPIGlobalComm(PE)
!Set Ray Tracing Communicator
CALL SetMPIRTComm(PE)
!Set CMFD Communicator
CALL SetCMFDComm(PE)
PE%MPI_RTMASTER_COMM = PE%MPI_CMFD_COMM
END SUBROUTINE

SUBROUTINE SetMpiRTComm(PE)
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
INTEGER :: ikey, icolor,ierr
INTEGER :: nproc, nproc0, myrank, myrank0, comm

PE%MPI_RT_COMM = MPI_COMM_NULL
IF(.NOT. PE%lIdle) THEN
  myrank = PE%myrank; nproc = PE%nproc
  icolor = myrank / PE%nRtProc
  ikey = myrank
  CALL MPI_COMM_SPLIT(PE%MPI_NTRACER_COMM, icolor, ikey, comm, ierr)
  PE%MPI_RT_COMM = comm
  !CALL MPInproc0
  CALL MPI_COMM_RANK(COMM, myrank0, ierr)
  CALL MPI_COMM_SIZE(COMM, nproc0, ierr)
  PE%MPI_RT_COMM = COMM; PE%myRTrank = myrank0
  PE%RTMaster = .TRUE.
  IF(PE%myRTrank .NE. 0) PE%RTMaster = .FALSE.
  PE%RTSlave = .NOT. PE%RTMaster
ENDIF
IF(PE%MPI_RT_COMM .EQ. MPI_COMM_NULL) PE%lRTgrp = .FALSE.
END SUBROUTINE

FUNCTION GetSubComm(comm, myrank, nproc,lIn)
INTEGER :: GetSubComm
INTEGER :: COMM, myrank, nproc
LOGICAL :: lIn
INTEGER :: Group0, NewGroup, NEW_COMM
INTEGER :: i, ierr
INTEGER :: InclProc0(0:1000), InclProc(0:1000), nGrProc
CALL MPI_SYNC(comm)
GetSubComm = MPI_COMM_NULL
InclProc = 0; InclProc0 = 0
IF(lIn) InclProc(myrank) = 1
CALL REDUCE(InclProc(0:nproc-1), InclProc0(0:nproc-1), nproc, COMM, .TRUE.)
InclProc = 0; nGrProc = 0
DO i = 0, nproc-1
  IF(InclProc0(i) .NE. 0) THEN
    nGrProc = nGrProc + 1;  InclProc(nGrProc) = i
  ENDIF
ENDDO
CALL MPI_COMM_GROUP(COMM, Group0, ierr)
CALL MPI_GROUP_INCL(Group0, nGrProc, InclProc(1:nGrProc), NewGroup, ierr)
CALL MPI_COMM_CREATE(COMM, NewGroup, New_Comm, ierr)
GetSubComm = New_Comm
END FUNCTION

SUBROUTINE FreeComm(comm)
INTEGER :: comm
INTEGER :: ierr
CALL MPI_COMM_FREE(comm, ierr)
END SUBROUTINE

SUBROUTINE SetCMFDComm(PE)
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
INTEGER :: InclProc0(0:1000), InclProc(0:1000), nGrProc
INTEGER :: Group0, NewGroup, NEW_COMM
INTEGER :: i, ierr
INTEGER :: myrank, nproc, myrank0, nproc0
PE%MPI_CMFD_COMM = MPI_COMM_NULL
IF(.NOT. PE%lIdle) THEN
  myrank = PE%myrank; nproc = PE%nproc
  InclProc = 0
  IF(PE%RTMaster) InclProc(myrank) = 1
  !Get RT master Processor Rank
  InclProc0 = 0
  CALL REDUCE(InclProc(0:nproc-1), InclProc0(0:nproc-1), nproc, PE%MPI_NTRACER_COMM, .TRUE.)
  InclProc = 0; nGrProc = 0
  DO i = 0, nproc-1
    IF(InclProc0(i) .NE. 0) THEN
      nGrProc = nGrProc + 1
      InclProc(nGrProc) = i
    ENDIF
  ENDDO
  CALL MPI_COMM_GROUP(PE%MPI_NTRACER_COMM, Group0, ierr)
  CALL MPI_GROUP_INCL(Group0, nGrProc, InclProc(1:nGrProc), NewGroup, ierr)
  CALL MPI_COMM_CREATE(PE%MPI_NTRACER_COMM, NewGroup, New_Comm, ierr)
  PE%MPI_CMFD_COMM = New_Comm;
  IF(PE%RTMaster) THEN
    CALL MPI_COMM_RANK(New_Comm, myrank0, ierr)
    CALL MPI_COMM_SIZE(New_Comm, nproc0, ierr)
    IF(nproc0 .NE. PE%nCmfdProc) CALL Terminate('MPI Communicator Setting Error')
    PE%myCMFDRank = myrank0
    PE%nCMFDGrp = PE%myCMFDRank
    PE%CmfdMaster = .TRUE.
    IF(myrank0 .NE. 0) PE%CmfdMaster = .FALSE.
    PE%CmfdSlave = .NOT.PE%CmfdMaster
  ENDIF
  CALL BCAST(PE%nCMFDGrp, PE%MPI_RT_COMM)
  !CALL MPI_COMM_SIZE(PE%MPI_RT_COMM, nproc0, ierr)
  !CALL MPIWaitTurn(PE%MPI_NTRACER_COMM, PE%myrank, PE%nproc, .FALSE.)
  !PRINT *, PE%myrank, PE%myRTrank, PE%nCMFDGrp, nproc0, PE%RTMaster
  !CALL MPIWaitTurn(PE%MPI_NTRACER_COMM, PE%myrank, PE%nproc, .TRUE.)
ENDIF
IF(PE%MPI_CMFD_COMM .EQ. MPI_COMM_NULL) PE%lCMFDgrp = .FALSE.

END SUBROUTINE


SUBROUTINE SetMPIGlobalComm(PE)
USE TypeDef,    ONLY : CoreInfo_Type
USE IOUTIL,     ONLY : TERMINATE
USE GEOM,   ONLY : Core
IMPLICIT NONE
!TYPE(CoreInfo_Type) :: Core
TYPE(PE_TYPE) :: PE
INTEGER :: ierr
INTEGER :: i, j, k
INTEGER :: n, n1, n2
INTEGER :: nz, nxy, nproc0, nproc, nCmfdProc
INTEGER :: MPI_GROUP, MPI_NTRACER_GROUP, MPI_NTRACER_COMM
INTEGER :: ExclProc(100)
!PROCESSOR Number Control
!PE%nRTproc = PE%nproc
PE%nRTproc = 1
nproc0 = PE%nproc
n = nproc0 / PE%nRTproc
IF(n .EQ. 0) THEN
  CALL MPI_FINALIZE(ierr)
  CALL TERMINATE('NOT ENOUGH PROCESSOR')
ENDIF
nproc = n * PE%nRTproc; nCMFDProc = n

nz = Core%nz; nxy = Core%nxy
n1 = Core%nz / nCMFDProc; n2 = mod(Core%nz, nCMFDProc)
IF(n2 .NE. 0  .AND. .NOT. PE%lUsrAxDcp) THEN
  CALL MPI_FINALIZE(ierr)
  CALL TERMINATE('NOT PROPER # PROCESSOR')
ENDIF
PE%nCmfdProc = nCMFDProc
IF(nCMFDProc .GT. 1) THEN
  PE%AxialParallel = .TRUE.
ENDIF
!Set Global Communicator
IF(PE%myrank .GT. nproc-1) THEN
  PE%Lidle = .TRUE.
ENDIF
n = 0
DO i = nProc+1, nProc0
  n = n + 1
  ExclProc(n) = i - 1
ENDDO

CALL MPI_COMM_GROUP(MPI_COMM_WORLD, MPI_GROUP, ierr)
CALL MPI_GROUP_EXCL(MPI_Group, n, ExclProc, MPI_NTRACER_GROUP, ierr)
CALL MPI_COMM_CREATE(MPI_COMM_WORLD, MPI_NTRACER_GROUP, MPI_NTRACER_COMM, ierr)
PE%nproc0 = PE%nproc; PE%myrank0 = PE%myrank
PE%nproc = 0;  PE%myrank = 0
!Get New Communicator
IF(.NOT. PE%Lidle) THEN
  PE%MPI_NTRACER_COMM = MPI_NTRACER_COMM
  CALL MPI_COMM_SIZE(MPI_NTRACER_COMM, PE%nproc, ierr)
  CALL MPI_COMM_RANK(MPI_NTRACER_COMM, PE%myrank, ierr)
ENDIF

END SUBROUTINE


SUBROUTINE PEInitialize()
USE PARAM
USE TypeDef,     ONLY : PE_TYPE
USE PE_Mod,      ONLY : PE

IMPLICIT NONE
#ifdef MPI_ENV
include 'mpif.h'
#endif


INTEGER :: ierr

#ifdef MPI_ENV
CALL MPI_INIT(iErr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, PE%myrank, ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, PE%nproc, ierr)
PE%MPI_COMM = MPI_COMM_WORLD

PE%master = TRUE
IF(PE%myrank .NE. 0) PE%master = FALSE
PE%SLAVE = .NOT. PE%master
#else

#endif
ENDSUBROUTINE

SUBROUTINE PeFinalize()
USE PARAM
IMPLICIT NONE
INTEGER :: ierr
#ifdef MPI_ENV
CALL MPI_FINALIZE(ierr)
#endif
END SUBROUTINE



SUBROUTINE SetGeomPEVariables(PE)
USE PARAM
USE TYPEDEF, ONLY : PE_TYPE
USE GEOM,    ONLY : Core, AsyInfo, nz, nzfm, nSubPlane
USE CNTL,    ONLY : nTracerCntl
USE ALLOCS
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
INTEGER :: n, nxy, CMFDGrp, Buf(2, 0:1000)
INTEGER :: AsyBeg, AsyEnd, PinBeg, PinEnd, FsrBeg, FsrEnd
INTEGER :: ixya, iasy, ixy, iz
INTEGER :: i, j
#ifdef MPI_ENV
n = PE%nz / PE%nCMFDproc
CMFDGrp = PE%nCMFDGrp
IF(PE%lUsrAxDcp .AND. PE%UsrAxDecomp(0) .EQ. PE%nCMFDproc) THEN
  PE%myzb = 0; PE%myzbf = 0
  DO j = 1, CMFDGrp
    PE%myzb = PE%myzb + PE%UsrAxDecomp(j)
  ENDDO
  PE%myzbf = PE%myzb * nSubPlane; PE%myzef = PE%myzbf +  PE%UsrAxDecomp(CmfdGrp+1)*nSubPlane;
  PE%myzbf = PE%myzbf + 1
  PE%myze = PE%myzb + PE%UsrAxDecomp(CmfdGrp+1); PE%myzb = PE%myzb + 1
  !PRINT *, CmfdGrp, PE%myzbf, PE%myzef
ELSE
  PE%myzb = CMFDGrp * n+1; PE%myze = (CMFDGrp+1) * n
  !PE%myzbf = CMFDGrp * nSubPlane * n+1; PE%myzef = (CMFDGrp + 1) * nSubPlane * n
  PE%myzbf = (PE%myzb - 1) * nSubPlane + 1;
  PE%myzef = PE%myze * nSubPlane
ENDIF

PE%nzfm = PE%nz * nSubPlane

!GetRangeDecomp
IF (nTracerCntl%lHex) THEN ! KSC
  nxy = Core%nxy
ELSE
  nxy = 0
  DO ixya = 1, Core%nxya
    iasy = Core%CoreMap(ixya)
    nxy = nxy + AsyInfo(iasy)%nxy
  ENDDO
END IF

PE%nxy = nxy

PE%AxDomRange(1:2, 0) = (/1, PE%nz/); PE%RadDomRange(1:2, 0) = (/1, nxy/)
CALL GetRangeDecomp(1, nxy, PE%nCMFDproc, CMFDGrp, PE%myNxyBeg, PE%myNxyEnd)
IF(PE%lCmfdGrp) THEN
  PE%AxDomRange(1:2, 0:PE%nCMFDproc) = 0; PE%AxDomRange(1:2, 0:PE%nCMFDproc) = 0
  Buf(1:2, 0:PE%nCMFDproc) = 0;   Buf(1:2, PE%myCMFDrank) = (/PE%myzb, PE%myze/)
  CALL Reduce(Buf(1:2, 0:PE%nCMFDproc-1), PE%AxDomRange(1:2, 0:PE%nCMFDproc-1), 2, PE%nCMFDproc, PE%MPI_CMFD_COMM, TRUE)
  Buf(1:2, 0:PE%nCMFDproc) = 0; Buf(1:2, PE%myCMFDrank) = (/PE%MynxyBeg, PE%MyNxyEnd/)
  CALL Reduce(Buf(1:2, 0:PE%nCMFDproc-1), PE%RadDomRange(1:2, 0:PE%nCMFDproc-1), 2, PE%nCMFDproc, PE%MPI_CMFD_COMM, TRUE)
  CALL Dmalloc(PE%AxDomList,nz); CALL Dmalloc(PE%RadDomList,nxy)
  DO j = 0, PE%nCMFDProc - 1
    DO i = PE%RadDomRange(1, j), PE%RadDomRange(2, j)
      PE%RadDomList(i) = j
    ENDDO
    DO i = PE%AxDomRange(1, j), PE%AxDomRange(2, j)
      PE%AxDomList(i) = j
    ENDDO
  ENDDO

  DO j = 0, PE%nCMFDProc - 1
    PE%SubPlnDomRange(1, j) = nSubPlane * (PE%AxDomRange(1, j) - 1) + 1
    PE%SubPlnDomRange(2, j) = nSubPlane * PE%AxDomRange(2, j)
  ENDDO
ENDIF
!
IF(PE%nCMFDproc .GT. 1) PE%lAxSolParallel = .TRUE.
!PE%myzb = 1; PE%myze = nz
!PE%myzbf = 1; PE%myzef = nz*nSubPlane
#else
PE%myzb = 1; PE%myze = nz
PE%myzbf = 1; PE%myzef = nz*nSubPlane
PE%nzfm = PE%myzef
#endif

END SUBROUTINE

!--- CNJ Edit : Domain Decomposition + MPI
SUBROUTINE SetDcmpPEVariables(PE)
USE PARAM
USE TYPEDEF, ONLY : PE_TYPE
USE GEOM,    ONLY : Core
USE CNTL,    ONLY : nTracerCntl
USE ALLOCS
USE HexData, ONLY : RodPin, hAsy
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
INTEGER :: AsyBeg, AsyEnd, PinBeg, PinEnd, FsrBeg, FsrEnd
INTEGER :: i, j

CALL Dmalloc0(PE%nAsy, 0, PE%nRTProc - 1)
CALL Dmalloc0(PE%nPin, 0, PE%nRTProc - 1)
CALL Dmalloc0(PE%nFsr, 0, PE%nRTProc - 1)
CALL Dmalloc0(PE%Asy_displs, 0, PE%nRTProc - 1)
CALL Dmalloc0(PE%Pin_displs, 0, PE%nRTProc - 1)
CALL Dmalloc0(PE%Fsr_displs, 0, PE%nRTProc - 1)
DO i = 0, PE%nRTProc - 1
  CALL GetRangeDecomp(1, Core%nxya, PE%nRTProc, i, AsyBeg, AsyEnd)
  PinBeg = Core%Asy(AsyBeg)%GlobalPinIdx(1)
  
  IF (nTracerCntl%lHex) THEN
    PinEnd = hAsy(AsyEnd)%PinIdxSt + hAsy(AsyEnd)%nTotPin - 1
  ELSE
    PinEnd = Core%Asy(AsyEnd)%GlobalPinIdx(Core%AsyInfo(Core%Asy(AsyEnd)%AsyType)%nxy)
  END IF
  
  FsrBeg = Core%Pin(PinBeg)%FsrIdxSt
  FsrEnd = Core%Pin(PinEnd)%FsrIdxSt + Core%PinInfo(Core%Pin(PinEnd)%PinType)%nFsrMax - 1
  
  PE%nAsy(i) = AsyEnd - AsyBeg + 1
  PE%nPin(i) = PinEnd - PinBeg + 1
  PE%nFsr(i) = FsrEnd - FsrBeg + 1
  IF (i .EQ. PE%myRTRank) THEN
    PE%myAsyBeg = AsyBeg; PE%myAsyEnd = AsyEnd
    PE%myPinBeg = PinBeg; PE%myPinEnd = PinEnd
    PE%myFsrBeg = FsrBeg; PE%myFsrEnd = FsrEnd
  ENDIF
  IF (i .NE. 0) THEN
    PE%Asy_displs(i) = PE%Asy_displs(i - 1) + PE%nAsy(i - 1)
    PE%Pin_displs(i) = PE%Pin_displs(i - 1) + PE%nPin(i - 1)
    PE%Fsr_displs(i) = PE%Fsr_displs(i - 1) + PE%nFsr(i - 1)
  ENDIF
ENDDO

END SUBROUTINE

SUBROUTINE SetRayPEVariables(PE, Core, RayInfo)
USE PARAM

USE TYPEDEF, ONLY : PE_TYPE, RayInfo_Type, CoreInfo_Type
USE UtilFunction, ONLY : Array2DSORT
USE Allocs
USE HexData, ONLY : hRotRay
USE CNTL,    ONLY : nTracerCntl
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: Core
INTEGER :: nRotray
INTEGER :: ibeg, iend, irank
INTEGER :: iray, nseg
INTEGER :: i
REAL, POINTER :: SegList(:, :)
INTEGER, POINTER :: RayList(:, :)
INTEGER :: nRay(0:1000)

!--- CNJ Edit : Useless

! #ifdef MPI_ENV
! nRotRay = RayInfo%nRotRay
! CALL GetRangeDecomp(1, nRotRay, PE%nRTProc, PE%myRTrank, ibeg, iend)
! PE%myRayBeg = ibeg; PE%myRayEnd = iEnd
! 
! ALLOCATE(segList(2, nRotRay))
! ALLOCATE(RayList(nRotRay, 0:PE%nRTproc))
! IF (nTracerCntl%lHex) THEN
!   DO iray = 1, nRotRay
!     SegList(1, iray) = iray
!     SegList(2, iray) = hRotRay(iray)%nSegMax
!   ENDDO
! ELSE
!   DO iray = 1, nRotRay
!     SegList(1, iray) = iray
!     SegList(2, iray) = RayInfo%RotRay(iray)%nseg
!   ENDDO
! END IF
! CALL Array2DSort(SegList(1, :), SegList(2, :), nRotRay, nseg,.TRUE., False, 2)
! nRay(0:PE%nRTProc) = 0
! DO i = 0, nROtRay-1
!   irank = mod(i, PE%nRtProc)
!   iray = NINT(SegList(1, i+1))
!   nRay(irank) = nRay(irank) + 1
!   RayList(nRay(irank), irank) = iray
! ENDDO
! PE%nRay = nRay(PE%myRTrank)
! CALL Dmalloc(PE%RayList, PE%nRay)
! PE%RayList(1:PE%nRay)=RayList(1:PE%nRay, PE%myRTrank)
! DEALLOCATE(RayList)
! CALL Dmalloc0(PE%OmpRayList, 0, nRotRay, 1, PE%nThread)
! DO i = 0, nROtRay-1
!   irank = mod(i, PE%nThread)+1
!   iray = NINT(SegList(1, i+1))
!   PE%OmpRayList(0, irank) = PE%OmpRayList(0, irank) + 1
!   PE%OmpRayList(PE%OmpRayList(0, irank), irank) = iray
! ENDDO
! 
! #else
! PE%myRayBeg = 1; PE%myRayEnd = RayInfo%nRotRay
! #endif

IF(Core%nCoreFsr .LT. PE%nThread) THEN
  PE%myOmpFsrBeg(1) = 1; PE%myOmpFsrEnd(1) = Core%nCoreFsr
  PE%myOmpFsrBeg(2:100) = 1; PE%myOmpFsrEnd(2:100) = -1
ELSE
  DO i = 1, PE%nThread
    CALL GetRangeDecomp(1, Core%nCoreFsr, PE%nThread, i-1, PE%myOmpFsrBeg(i), PE%myOmpFsrEnd(i))
  ENDDO
ENDIF
IF(Core%nxy .LT. PE%nThread) THEN
  PE%myOmpNxyBeg(1) = 1; PE%myOmpNxyEnd(1) = Core%nxy
  PE%myOmpNxyBeg(2:100) = 1; PE%myOmpNxyEnd(2:100) = -1
ELSE
  DO i = 1, PE%nThread
    CALL GetRangeDecomp(1, Core%nxy, PE%nThread, i-1, PE%myOmpNxyBeg(i), PE%myOmpNxyEnd(i))
  ENDDO

ENDIF

END SUBROUTINE



END MODULE
