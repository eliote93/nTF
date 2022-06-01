#include <defines.h>
MODULE MPIConfig_Mod

USE PARAM,        ONLY : TRUE, FALSE
USE TypeDef,      ONLY : PE_TYPE
USE MPICOMM_MOD,  ONLY : REDUCE
USE IOUTIL,       ONLY : Terminate
USE UtilFunction, ONLY : GetDcmpRange

IMPLICIT NONE

#ifdef MPI_ENV
INCLUDE 'mpif.h'
#endif

CONTAINS
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE PEInitialize()

USE TypeDef, ONLY : PE_TYPE
USE PE_Mod,  ONLY : PE

IMPLICIT NONE

#ifdef MPI_ENV
include 'mpif.h'
#endif

INTEGER :: ierr
! ----------------------------------------------------

#ifdef MPI_ENV
CALL MPI_INIT(ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, PE%myrank, ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, PE%nproc,  ierr)

PE%MPI_COMM = MPI_COMM_WORLD

PE%master = TRUE
IF (PE%myrank .NE. 0) PE%master = FALSE
PE%SLAVE = .NOT. PE%master

CALL hostnm(PE%hostname)
#endif
! ----------------------------------------------------

END SUBROUTINE PEInitialize
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE PeFinalize()

IMPLICIT NONE

INTEGER :: ierr

#ifdef MPI_ENV
CALL MPI_FINALIZE(ierr)
#endif

END SUBROUTINE PeFinalize
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetMPIEnv(PE)

IMPLICIT NONE

TYPE (PE_TYPE) :: PE

CALL SetMPIGlobalComm(PE)
CALL SetMPIRTComm(PE)
CALL SetCMFDComm(PE)

PE%MPI_RTMASTER_COMM = PE%MPI_CMFD_COMM

END SUBROUTINE SetMPIEnv
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetMPIGlobalComm(PE)

USE TypeDef, ONLY : CoreInfo_Type
USE GEOM,    ONLY : Core

IMPLICIT NONE

TYPE (PE_TYPE) :: PE

INTEGER :: ierr, i, j, k, n, n1, n2, nz, nxy, nproc0, nproc, nCmfdProc
INTEGER :: MPI_GROUP, MPI_NTRACER_GROUP, MPI_NTRACER_COMM
INTEGER :: ExclProc(100)
! ----------------------------------------------------

! PROCESSOR Number Control
PE%nRTproc = 1
nproc0 = PE%nproc
n = nproc0 / PE%nRTproc

IF (n .EQ. 0) THEN
  CALL MPI_FINALIZE(ierr)
  CALL TERMINATE('NOT ENOUGH PROCESSOR')
END IF

nproc = n * PE%nRTproc
nCMFDProc = n

nz  = Core%nz
nxy = Core%nxy
n1  = Core%nz / nCMFDProc
n2  = mod(Core%nz, nCMFDProc)

IF (n2 .NE. 0  .AND. .NOT.PE%lUsrAxDcp) THEN
  CALL MPI_FINALIZE(ierr)
  CALL TERMINATE('NOT PROPER # PROCESSOR')
END IF

PE%nCmfdProc = nCMFDProc

IF (nCMFDProc .GT. 1) PE%AxialParallel = TRUE

! Set Global Communicator
IF (PE%myrank .GT. nproc-1) PE%Lidle = TRUE

n = 0
DO i = nProc+1, nProc0
  n = n + 1
  ExclProc(n) = i - 1
END DO

CALL MPI_COMM_GROUP(PE%MPI_COMM, MPI_GROUP, ierr)
CALL MPI_GROUP_EXCL(MPI_Group, n, ExclProc, MPI_NTRACER_GROUP, ierr)
CALL MPI_COMM_CREATE(PE%MPI_COMM, MPI_NTRACER_GROUP, MPI_NTRACER_COMM, ierr)

PE%nproc0  = PE%nproc
PE%myrank0 = PE%myrank
PE%nproc   = 0
PE%myrank  = 0

! Get New Communicator
IF (.NOT. PE%Lidle) THEN
  PE%MPI_NTRACER_COMM = MPI_NTRACER_COMM
  
  CALL MPI_COMM_SIZE(MPI_NTRACER_COMM, PE%nproc, ierr)
  CALL MPI_COMM_RANK(MPI_NTRACER_COMM, PE%myrank, ierr)
END IF
! ----------------------------------------------------

END SUBROUTINE SetMPIGlobalComm
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetMpiRTComm(PE)

IMPLICIT NONE

TYPE (PE_TYPE) :: PE

INTEGER :: ikey, icolor, ierr, nproc, nproc0, myrank, myrank0, comm
! ----------------------------------------------------

PE%MPI_RT_COMM = MPI_COMM_NULL

IF (.NOT. PE%lIdle) THEN
  myrank = PE%myrank
  nproc  = PE%nproc
  icolor = myrank / PE%nRtProc
  ikey   = myrank
  
  CALL MPI_COMM_SPLIT(PE%MPI_NTRACER_COMM, icolor, ikey, comm, ierr)
  
  PE%MPI_RT_COMM = comm
  
  CALL MPI_COMM_RANK(COMM, myrank0, ierr)
  CALL MPI_COMM_SIZE(COMM, nproc0, ierr)
  
  PE%MPI_RT_COMM = COMM
  PE%myRTrank    = myrank0
  
  PE%RTMaster = TRUE
  IF (PE%myRTrank .NE. 0) PE%RTMaster = FALSE
  
  PE%RTSlave = .NOT.PE%RTMaster
END IF

IF (PE%MPI_RT_COMM .EQ. MPI_COMM_NULL) PE%lRTgrp = FALSE
! ----------------------------------------------------

END SUBROUTINE SetMpiRTComm
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetCMFDComm(PE)

USE MPICOMM_MOD, ONLY : BCAST

IMPLICIT NONE

TYPE(PE_TYPE) :: PE

INTEGER :: InclProc0(0:1000), InclProc(0:1000)
INTEGER :: Group0, NewGroup, NEW_COMM, i, ierr, myrank, nproc, myrank0, nproc0, nGrProc
! ----------------------------------------------------

PE%MPI_CMFD_COMM = MPI_COMM_NULL

IF (.NOT. PE%lIdle) THEN
  myrank = PE%myrank
  nproc  = PE%nproc
  
  InclProc = 0
  IF (PE%RTMaster) InclProc(myrank) = 1
  
  ! Get RT master Processor Rank
  InclProc0 = 0
  
  CALL REDUCE(InclProc(0:nproc-1), InclProc0(0:nproc-1), nproc, PE%MPI_NTRACER_COMM, TRUE)
  
  InclProc = 0
  nGrProc  = 0
  
  DO i = 0, nproc-1
    IF (InclProc0(i) .NE. 0) THEN
      nGrProc = nGrProc + 1
      
      InclProc(nGrProc) = i
    END IF
  END DO
  
  CALL MPI_COMM_GROUP(PE%MPI_NTRACER_COMM, Group0, ierr)
  CALL MPI_GROUP_INCL(Group0, nGrProc, InclProc(1:nGrProc), NewGroup, ierr)
  CALL MPI_COMM_CREATE(PE%MPI_NTRACER_COMM, NewGroup, New_Comm, ierr)
  
  PE%MPI_CMFD_COMM = New_Comm
  
  IF (PE%RTMaster) THEN
    CALL MPI_COMM_RANK(New_Comm, myrank0, ierr)
    CALL MPI_COMM_SIZE(New_Comm, nproc0, ierr)
    
    IF (nproc0 .NE. PE%nCmfdProc) CALL Terminate('MPI Communicator Setting Error')
    
    PE%myCMFDRank = myrank0
    PE%nCMFDGrp   = PE%myCMFDRank
    PE%CmfdMaster = TRUE
    IF (myrank0 .NE. 0) PE%CmfdMaster = FALSE
    
    PE%CmfdSlave = .NOT.PE%CmfdMaster
  END IF
  
  CALL BCAST(PE%nCMFDGrp, PE%MPI_RT_COMM)
END IF

IF (PE%MPI_CMFD_COMM .EQ. MPI_COMM_NULL) PE%lCMFDgrp = FALSE
! ----------------------------------------------------

END SUBROUTINE SetCMFDComm
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetGeomPEVariables(PE)

USE ALLOCS
USE TYPEDEF, ONLY : PE_TYPE
USE GEOM,    ONLY : Core, AsyInfo, nz, nzfm, nSubPlane
USE CNTL,    ONLY : nTracerCntl

IMPLICIT NONE

TYPE (PE_TYPE) :: PE
INTEGER :: Buf(2, 0:1000)
INTEGER :: n, nxy, CMFDGrp, AsyBeg, AsyEnd, PinBeg, PinEnd, FsrBeg, FsrEnd, ixya, iasy, ixy, iz, i, j
! ----------------------------------------------------

#ifdef MPI_ENV
n = PE%nz / PE%nCMFDproc
CMFDGrp = PE%nCMFDGrp

! GET : Ax. Range
IF (PE%lUsrAxDcp .AND. PE%UsrAxDecomp(0).EQ.PE%nCMFDproc) THEN
  PE%myzb  = 0
  PE%myzbf = 0
  
  DO j = 1, CMFDGrp
    PE%myzb = PE%myzb + PE%UsrAxDecomp(j)
  END DO
  
  PE%myzbf = PE%myzb  * nSubPlane
  PE%myzef = PE%myzbf + PE%UsrAxDecomp(CmfdGrp+1)*nSubPlane
  PE%myzbf = PE%myzbf + 1
  PE%myze  = PE%myzb  + PE%UsrAxDecomp(CmfdGrp+1)
  PE%myzb  = PE%myzb  + 1
ELSE
  PE%myzb  = n * CMFDGrp + 1
  PE%myze  = n * (CMFDGrp+1)
  PE%myzbf = nSubPlane * (PE%myzb - 1) + 1
  PE%myzef = nSubPlane * PE%myze
END IF

PE%nzfm = PE%nz * nSubPlane

! Get : Dcmp Range
IF (nTracerCntl%lHex) THEN ! KSC
  nxy = Core%nxy
ELSE
  nxy = 0
  
  DO ixya = 1, Core%nxya
    iasy = Core%CoreMap(ixya)
    nxy = nxy + AsyInfo(iasy)%nxy
  END DO
END IF

PE%nxy = nxy

PE%AxDomRange (1:2, 0) = (/1, PE%nz/)
PE%RadDomRange(1:2, 0) = (/1, nxy/)

CALL GetDcmpRange(1, nxy, PE%nCMFDproc, CMFDGrp, PE%myNxyBeg, PE%myNxyEnd)

! CMFD Grp.
IF (PE%lCmfdGrp) THEN
  PE%AxDomRange(1:2, 0:PE%nCMFDproc) = 0
  PE%AxDomRange(1:2, 0:PE%nCMFDproc) = 0
  
  Buf(1:2, 0:PE%nCMFDproc) = 0
  Buf(1:2, PE%myCMFDrank) = (/PE%myzb, PE%myze/)
  
  CALL Reduce(Buf(1:2, 0:PE%nCMFDproc-1), PE%AxDomRange(1:2, 0:PE%nCMFDproc-1), 2, PE%nCMFDproc, PE%MPI_CMFD_COMM, TRUE)
  
  Buf(1:2, 0:PE%nCMFDproc) = 0
  Buf(1:2, PE%myCMFDrank) = (/PE%MynxyBeg, PE%MyNxyEnd/)
  
  CALL Reduce(Buf(1:2, 0:PE%nCMFDproc-1), PE%RadDomRange(1:2, 0:PE%nCMFDproc-1), 2, PE%nCMFDproc, PE%MPI_CMFD_COMM, TRUE)
  
  CALL dmalloc(PE%AxDomList,  nz)
  CALL dmalloc(PE%RadDomList, nxy)
  
  DO j = 0, PE%nCMFDProc - 1
    DO i = PE%RadDomRange(1, j), PE%RadDomRange(2, j)
      PE%RadDomList(i) = j
    END DO
    
    DO i = PE%AxDomRange(1, j), PE%AxDomRange(2, j)
      PE%AxDomList(i) = j
    END DO
  END DO

  DO j = 0, PE%nCMFDProc - 1
    PE%SubPlnDomRange(1, j) = nSubPlane * (PE%AxDomRange(1, j) - 1) + 1
    PE%SubPlnDomRange(2, j) = nSubPlane * PE%AxDomRange(2, j)
  END DO
END IF

IF (PE%nCMFDproc .GT. 1) PE%lAxSolParallel = TRUE
#else
PE%myzb  = 1
PE%myze  = nz
PE%myzbf = 1
PE%myzef = nz*nSubPlane
PE%nzfm  = PE%myzef
#endif
! ----------------------------------------------------

END SUBROUTINE SetGeomPEVariables
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetGeomPEVariables_subpln(PE)

USE ALLOCS
USE TYPEDEF, ONLY : PE_TYPE
USE GEOM,    ONLY : Core, AsyInfo, nz, nzfm, nSubPlane, SubPlaneRange
USE CNTL,    ONLY : nTracerCntl
USE HexData, ONLY : nHexPin

IMPLICIT NONE

TYPE (PE_TYPE) :: PE

INTEGER :: Buf(2, 0:1000)
INTEGER :: n, nxy, CMFDGrp, AsyBeg, AsyEnd, PinBeg, PinEnd, FsrBeg, FsrEnd, ixya, iasy, ixy, iz, i, j
! ----------------------------------------------------

#ifdef MPI_ENV
n = PE%nz / PE%nCMFDproc
CMFDGrp = PE%nCMFDGrp

! GET : Ax. Range
IF (PE%lUsrAxDcp .AND. PE%UsrAxDecomp(0) .EQ. PE%nCMFDproc) THEN
  PE%myzb = 0
  PE%myzbf = 0
  
  DO j = 1, CMFDGrp
    PE%myzb = PE%myzb + PE%UsrAxDecomp(j)
  END DO
  
  PE%myze  = PE%myzb + PE%UsrAxDecomp(CmfdGrp+1)
  PE%myzb  = PE%myzb + 1
  PE%myzbf = SubPlaneRange(1, PE%myzb)
  PE%myzef = SubPlaneRange(2, PE%myze)
ELSE
  PE%myzb  = CMFDGrp * n+1
  PE%myze  = (CMFDGrp+1) * n
  PE%myzbf = SubPlaneRange(1, PE%myzb)
  PE%myzef = SubPlaneRange(2, PE%myze)
END IF

PE%nzfm = nzfm

! GET : Dcmp Range
IF (nTracerCntl%lHex) THEN ! KSC
  nxy = nHexPin
ELSE
  nxy = 0
  
  DO ixya = 1, Core%nxya
    iasy = Core%CoreMap(ixya)
    nxy = nxy + AsyInfo(iasy)%nxy
  END DO
END IF

PE%nxy = nxy

PE%AxDomRange(1:2, 0) = (/1, PE%nz/)
PE%RadDomRange(1:2, 0) = (/1, nxy/)

CALL GetDcmpRange(1, nxy, PE%nCMFDproc, CMFDGrp, PE%myNxyBeg, PE%myNxyEnd)

! CMFD Grp
IF (PE%lCmfdGrp) THEN
  PE%AxDomRange(1:2, 0:PE%nCMFDproc) = 0
  PE%AxDomRange(1:2, 0:PE%nCMFDproc) = 0
  
  Buf(1:2, 0:PE%nCMFDproc) = 0
  Buf(1:2, PE%myCMFDrank) = (/PE%myzb, PE%myze/)
  
  CALL Reduce(Buf(1:2, 0:PE%nCMFDproc-1), PE%AxDomRange(1:2, 0:PE%nCMFDproc-1), 2, PE%nCMFDproc, PE%MPI_CMFD_COMM, TRUE)
  
  Buf(1:2, 0:PE%nCMFDproc) = 0
  Buf(1:2, PE%myCMFDrank) = (/PE%MynxyBeg, PE%MyNxyEnd/)
  
  CALL Reduce(Buf(1:2, 0:PE%nCMFDproc-1), PE%RadDomRange(1:2, 0:PE%nCMFDproc-1), 2, PE%nCMFDproc, PE%MPI_CMFD_COMM, TRUE)
  
  CALL dmalloc(PE%AxDomList,  nz)
  CALL dmalloc(PE%RadDomList, nxy)
  
  DO j = 0, PE%nCMFDProc - 1
    DO i = PE%RadDomRange(1, j), PE%RadDomRange(2, j)
      PE%RadDomList(i) = j
    END DO
    
    DO i = PE%AxDomRange(1, j), PE%AxDomRange(2, j)
      PE%AxDomList(i) = j
    END DO
  END DO

  DO j = 0, PE%nCMFDProc - 1
    PE%SubPlnDomRange(1, j) = SubPlaneRange(1, PE%AxDomRange(1, j))
    PE%SubPlnDomRange(2, j) = SubPlaneRange(2, PE%AxDomRange(2, j))
  END DO
END IF

IF (PE%nCMFDproc .GT. 1) PE%lAxSolParallel = TRUE
#else
PE%myzb  = 1
PE%myze  = nz
PE%myzbf = 1
PE%myzef = nzfm !nz*nSubPlane
PE%nzfm  = PE%myzef
#endif
! ----------------------------------------------------

END SUBROUTINE SetGeomPEVariables_subpln
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetDcmpPEVariables(PE)

USE ALLOCS
USE TYPEDEF, ONLY : PE_TYPE
USE GEOM,    ONLY : Core
USE CNTL,    ONLY : nTracerCntl
USE HexData, ONLY : RodPin, hAsy

IMPLICIT NONE

TYPE (PE_TYPE) :: PE

INTEGER :: AsyBeg, AsyEnd, PinBeg, PinEnd, FsrBeg, FsrEnd
INTEGER :: iprc
! ----------------------------------------------------

CALL dmalloc0(PE%nAsy, 0, PE%nRTProc - 1)
CALL dmalloc0(PE%nPin, 0, PE%nRTProc - 1)
CALL dmalloc0(PE%nFsr, 0, PE%nRTProc - 1)

CALL dmalloc0(PE%Asy_displs, 0, PE%nRTProc - 1)
CALL dmalloc0(PE%Pin_displs, 0, PE%nRTProc - 1)
CALL dmalloc0(PE%Fsr_displs, 0, PE%nRTProc - 1)

DO iprc = 0, PE%nRTProc - 1
  CALL GetDcmpRange(1, Core%nxya, PE%nRTProc, iprc, AsyBeg, AsyEnd)
  
  PinBeg = Core%Asy(AsyBeg)%GlobalPinIdx(1)

  IF (nTracerCntl%lHex) THEN
    PinEnd = hAsy(AsyEnd)%PinIdxSt + hAsy(AsyEnd)%nTotPin - 1
  ELSE
    PinEnd = Core%Asy(AsyEnd)%GlobalPinIdx(Core%AsyInfo(Core%Asy(AsyEnd)%AsyType)%nxy)
  END IF
  
  FsrBeg = Core%Pin(PinBeg)%FsrIdxSt
  FsrEnd = Core%Pin(PinEnd)%FsrIdxSt + Core%PinInfo(Core%Pin(PinEnd)%PinType)%nFsrMax - 1
  
  PE%nAsy(iprc) = AsyEnd - AsyBeg + 1
  PE%nPin(iprc) = PinEnd - PinBeg + 1
  PE%nFsr(iprc) = FsrEnd - FsrBeg + 1
  
  IF (iprc .EQ. PE%myRTRank) THEN
    PE%myAsyBeg = AsyBeg; PE%myAsyEnd = AsyEnd
    PE%myPinBeg = PinBeg; PE%myPinEnd = PinEnd
    PE%myFsrBeg = FsrBeg; PE%myFsrEnd = FsrEnd
  END IF
  
  IF (iprc .NE. 0) THEN
    PE%Asy_displs(iprc) = PE%Asy_displs(iprc - 1) + PE%nAsy(iprc - 1)
    PE%Pin_displs(iprc) = PE%Pin_displs(iprc - 1) + PE%nPin(iprc - 1)
    PE%Fsr_displs(iprc) = PE%Fsr_displs(iprc - 1) + PE%nFsr(iprc - 1)
  END IF
END DO
! ----------------------------------------------------
  
END SUBROUTINE SetDcmpPEVariables
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRayPEVariables(PE, Core, RayInfo)

USE Allocs
USE TYPEDEF, ONLY : PE_TYPE, RayInfo_Type, CoreInfo_Type

IMPLICIT NONE

TYPE (PE_TYPE) :: PE
TYPE (RayInfo_Type) :: RayInfo
TYPE (CoreInfo_Type) :: Core

INTEGER :: nRotray, ibeg, iend, irank, iray, nseg, i
INTEGER :: nRay(0:1000)

REAL, POINTER, DIMENSION(:,:) :: SegList
INTEGER, POINTER, DIMENSION(:,:) :: RayList
! ----------------------------------------------------

IF (Core%nCoreFsr .LT. PE%nThread) THEN
  PE%myOmpFsrBeg(1) = 1
  PE%myOmpFsrEnd(1) = Core%nCoreFsr
  
  PE%myOmpFsrBeg(2:100) = 1
  PE%myOmpFsrEnd(2:100) = -1
ELSE
  DO i = 1, PE%nThread
    CALL GetDcmpRange(1, Core%nCoreFsr, PE%nThread, i-1, PE%myOmpFsrBeg(i), PE%myOmpFsrEnd(i))
  END DO
END IF

IF (Core%nxy .LT. PE%nThread) THEN
  PE%myOmpNxyBeg(1) = 1
  PE%myOmpNxyEnd(1) = Core%nxy
  
  PE%myOmpNxyBeg(2:100) = 1
  PE%myOmpNxyEnd(2:100) = -1
ELSE
  DO i = 1, PE%nThread
    CALL GetDcmpRange(1, Core%nxy, PE%nThread, i-1, PE%myOmpNxyBeg(i), PE%myOmpNxyEnd(i))
  END DO
END IF
! ----------------------------------------------------

END SUBROUTINE SetRayPEVariables
! ------------------------------------------------------------------------------------------------------------

END MODULE MPIConfig_Mod