#include <defines.h>
MODULE MPIDcpl_Mod
USE PARAM
USE TYPEDEF,      ONLY : PE_TYPE, DcplInfo_TYpe
IMPLICIT NONE

#ifdef MPI_ENV
INCLUDE 'mpif.h'
#endif

CONTAINS

SUBROUTINE InitDcplParallel()
USE PARAM
USE DcplCore_Mod,   ONLY : DcplInfo
USE PE_MOD,         ONLY : PE,           DcplPE
IMPLICIT NONE
!TYPE(PE_TYPE) :: PE,      DcplPE(100)
!TYPE(DcplInfo_Type) :: DcplInfo

CALL DcplDomainDecomp(PE, DcplInfo)
CALL SetDcplMPIComm(PE, DcplPE)

END SUBROUTINE

SUBROUTINE DcplDomainDeComp(PE, DcplInfo)
IMPLICIT NONE

TYPE(PE_TYPE) :: PE
TYPE(DcplInfo_Type) :: DcplInfo

INTEGER :: myrank
INTEGER :: i

myrank = PE%myCmfdRank

PE%myRefPlnBeg = DcplInfo%Pln_Map0(PE%myzb); PE%myRefPlnEnd = DcplInfo%Pln_Map0(PE%myze)
PE%RefPlnRange = PE%AxDomRange

PE%lDcplParallel = TRUE
IF(PE%nproc .GT. 1) PE%lDcplParallel = TRUE

DO i = PE%myRefPlnBeg, PE%myRefPlnEnd
  PE%lmyRefPln(i) = .TRUE.
ENDDO
END SUBROUTINE

SUBROUTINE SetDcplMPIComm(PE, DcplPE)
TYPE(PE_TYPE) :: PE, DcplPE(100)
INTEGER :: iz1, iz2, iz
INTEGER :: Group0, NewGroup, NEW_COMM
INTEGER :: Comm, MyRank, nProc
INTEGER :: InclProc(0:1000)
INTEGER :: I, IERR

iz1 = PE%myRefPlnBeg; iz2 = PE%myRefPlnEnd
!Set nThread
DO iz = iz1, iz2
  DcplPE(iz)%nThread = PE%nThread
ENDDO
!

COMM = PE%MPI_RT_Comm
Myrank = PE%myRTRank; nProc = PE%nRTproc
InclProc(1) = 0;
!
#ifdef MPI_ENV
PE%MPI_DCPLMASTER_COMM = PE%MPI_RTMASTER_COMM
CALL MPI_COMM_GROUP(Comm, Group0, ierr)
CALL MPI_GROUP_INCL(Group0, 1, InclProc(1:1), NewGroup, ierr)
CALL MPI_COMM_CREATE(COMM, NewGroup, New_Comm, ierr)
CALL MPI_COMM_RANK(New_Comm, myrank, iz)
!Set MPI Communicator
DO iz = iz1, iz2
  DcplPE(iz)%lDcplParallel = PE%lDcplParallel
  DcplPE(iz)%myrank = PE%myCmfdRank
  DcplPE(iz)%nproc = PE%nCmfdProc
  DcplPE(iz)%nCmfdProc = 1
  DcplPE(iz)%nRtProc = PE%nRtProc
  DcplPE(iz)%MPI_NTRACER_Comm = PE%MPI_RT_Comm
  DcplPE(iz)%MPI_CMFD_COMM = New_Comm
  DcplPE(iz)%MPI_RT_COMM = PE%MPI_RT_COMM
  DcplPE(iz)%MPI_RTMASTER_COMM = New_Comm
  DcplPE(iz)%nCmfdProc = 1;
  DcplPE(iz)%nRtProc = PE%nRtProc
  DcplPE(iz)%myRTrank = PE%myRTrank
  DcplPE(iz)%myCmfdRank = 0
  DcplPE(iz)%Master = PE%Master
  IF(PE%RTMaster) THEN
    DcplPE(iz)%RTMaster = .TRUE.; DcplPE(iz)%RTSlave = .FALSE.
    DcplPE(iz)%CmfdMaster = .TRUE.; DcplPE(iz)%CmfdSlave = .FALSE.
    DcplPE(iz)%lCmfdGrp = .TRUE.; DcplPE(iz)%lRtGrp = .TRUE.
  ELSE
    DcplPE(iz)%RTMaster = .FALSE.; DcplPE(iz)%RTSlave = .FALSE.
    DcplPE(iz)%CmfdMaster = .FALSE.; DcplPE(iz)%CmfdSlave = .FALSE.
    DcplPE(iz)%lCmfdGrp = .FALSE.; DcplPE(iz)%lRtGrp = .TRUE.
  ENDIF

  DcplPE(iz)%MPI_DCPLMaster_Comm = PE%MPI_DCPLMASTER_COMM
ENDDO
#endif
END SUBROUTINE

SUBROUTINE MPI_DcplMessage(amesg, myrank, ltime, lopen, lclose)
USE PARAM
USE FILES, ONLY : io30
USE IOUTIL, ONLY : Message
USE UtilFunction, ONLY : StrLen
IMPLICIT NONE
CHARACTER  amesg*(*)
INTEGER :: myrank
LOGICAL :: lopen, lclose, ltime
CHARACTER(256) :: filename
INTEGER :: io
INTEGER :: caseidsize

io = io30 + myrank
IF(lopen) THEN
  filename = GetDcplLogFileName(myrank)
  OPEN(unit = io, FILE = Filename, status = 'REPLACE')
ENDIF
IF(lClose) THEN
  CLOSE(io)
ENDIF
CALL MESSAGE(io, ltime, FALSE, amesg)
END SUBROUTINE

SUBROUTINE MPI_FinalizeDcplMesg(nproc)
USE PARAM
USE FILES,          ONLY : io8,        io30
USE inputcards ,    ONLY : oneline
USE IOUTIL,         ONLY : Message
USE UtilFunction,   ONLY : StrLen
IMPLICIT NONE
INTEGER :: nproc
INTEGER :: io
INTEGER :: i

CHARACTER(256) :: filename

io = io30
DO i = 0, nproc-1
  filename = GetDcplLogFileName(i)
  OPEN(UNIT = io, file = filename, status = 'OLD')
    DO
      READ(io, '(a256)', END = 1000) oneline
      WRITE(io8, '(a256)') oneline
    ENDDO
1000 continue
  CLOSE(io)
ENDDO
WRITE(io8, '(A)') hbar2(1:77)

END SUBROUTINE

FUNCTION GetDcplLogFileName(myrank)
USE FILES, ONLY : CASEID
USE UtilFunction, ONLY : StrLen
IMPLICIT NONE
INTEGER :: myrank
CHARACTER(256) :: GetDcplLogFileName

CHARACTER(256) :: file,Ext
INTEGER :: CaseIdSize, ExtSize

IF(myrank .LT. 10) THEN
  WRITE(Ext, '(A4,I1,A4)') 'Proc', myrank, '.tmp'
ELSEIF(myrank .LT. 100) THEN
  WRITE(Ext, '(A4,I2,A4)') 'Proc', myrank, '.tmp'
ELSE
  WRITE(Ext, '(A4,I3,A4)') 'Proc', myrank, '.tmp'
ENDIF

file = trim(caseid)
caseidsize = StrLen(file(1:256), 256)
ExtSize = StrLen(Ext(1:255), 256)
File(caseidsize+1:caseidsize+ExtSize) = Ext(1:ExtSize)

GetDcplLogFileName = File

END FUNCTION

SUBROUTINE PrintDcplRTLog(ig, time, myrank, nproc, COMM)
USE PARAM
USE FILES,            ONLY : io8
USE MpiComm_mod,      ONLY : MPI_SYNC,    REDUCE
USE IOUTIL,           ONLY : message
IMPLICIT NONE
REAL :: Time
INTEGER :: ig, myrank, nproc, COMM
REAL :: Buf(0:100), Buf0(0:100), time0

CALL MPI_SYNC(COMM)
buf = 0; buf0 = 0;
buf0(myrank) = time
CALL Reduce(Buf0(0:nproc-1), Buf(0:nproc-1), nproc, Comm, .FALSE.)
time0 = maxval(Buf(0:nproc-1))
write(mesg,'(10x, a, i4, 2x, a, f10.3, 2x, a)') 'Group ', ig, ' finished in ', Time0, 'Sec'
IF(myrank .eq. 0) CALL message(io8, FALSE, TRUE, mesg)
END SUBROUTINE


SUBROUTINE PrintDcplMOCLog(iter, iz, eigv, eigerr, fiserr, reserr, myrank, nproc, COMM)
USE PARAM
USE FILES,            ONLY : io8
USE MpiComm_mod,      ONLY : MPI_SYNC,      REDUCE
USE IOUTIL,           ONLY : message
IMPLICIT NONE
INTEGER :: iz, iter, myrank, COMM, nproc
INTEGER :: i
REAL :: eigv, reserr, fiserr, eigerr
REAL :: Buf0(1:6,300), Buf(1:6,100)

CALL MPI_SYNC(COMM)
Buf(1:6, 1:nproc) = 0
Buf0(1:6, myrank + 1) = (/real(iter,8), real(iz, 8), eigv, eigerr, fiserr, reserr/)
CALL REDUCE(BUf0(1:5, 1:nproc), BUf(1:5, 1:nproc), 5, nproc, COMM, FALSE)
IF(myrank .eq. 0) THEN
  write(mesg ,'(A5,I4, I7, F13.6, 3(1pE12.3))')  'RT', iter, NINT(Buf(2, 1)), Buf(3,1), buf(4,1),buf(5,1) !, buf(6,1)
  CALL message(io8, TRUE, TRUE, mesg)
  DO i = 2, nproc
    write(mesg ,'(18x, I4, I7, F13.6, 3(1pE12.3))') NINT(Buf(1, i)), NINT(Buf(2, i)), Buf(3,i), buf(4,i),buf(5,i) !,buf(6,i)
    CALL message(io8, FALSE, TRUE, mesg)
  ENDDO
ENDIF
END SUBROUTINE

SUBROUTINE PrintDcplCmfdLog(iter, iz, eigv, Ratio, ResErr, myrank, nproc, Comm)
USE PARAM
USE FILES,            ONLY : io8
USE MpiComm_mod,      ONLY : MPI_SYNC,      REDUCE
USE IOUTIL,           ONLY : message
IMPLICIT NONE
INTEGER :: iter, iz, myrank, nproc, Comm
INTEGER :: i
REAL :: Eigv, Ratio, ResErr
REAL :: Buf0(5, 300), Buf(5, 300)
Buf = 0; Buf0 = 0
CALL MPI_SYNC(Comm)
Buf0(1:5, myrank + 1) = (/REAL(iter, 8), REAL(iz, 8), Eigv, Ratio, ResErr/)
CALL REDUCE(Buf0(1:5, 1:nproc), Buf(1:5, 1:nproc), 5, nproc, comm, FALSE)
IF(myrank .eq. 0) THEN
  WRITE(mesg,'(A9, I9, I5, F15.6, 3x, F10.5, 1pe15.3)') 'MGOUTER', NINT(Buf(1, 1)), NINT(Buf(2,1)), Buf(3,1), Buf(4,1), Buf(5,1)
  CALL message(io8, FALSE, TRUE, mesg)
  DO i = 2, nproc
    WRITE(mesg,'(9X, I9, I5, F15.6, 3x, F10.5, 1pe15.3)') NINT(Buf(1, i)), NINT(Buf(2,i)), Buf(3,i), Buf(4,i), Buf(5,i)
    CALL message(io8, FALSE, TRUE, mesg)
  ENDDO
ENDIF
END SUBROUTINE

SUBROUTINE Idle_DcplPlnMOC(ng, GroupInfo, ItrCntl, DcplPE, lMOC)
USE PARAM
USE TypeDef,         ONLY : GroupInfo_Type,          PE_TYPE
USE ItrCntl_Mod,     ONLY : ItrCntl_Type
#ifdef MPI_ENV
USE MPIComm_mod,     ONLY : MPIConvChk
#endif
IMPLICIT NONE
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(ItrCntl_Type) :: ItrCntl
TYPE(PE_Type) :: DcplPE
INTEGER :: ng
LOGICAL :: lMOC

INTEGER :: nXsGenIterMax
INTEGER :: iRefPln, iRefTemp, iter
INTEGER :: Comm, myRank, nproc
REAL :: Eigv
LOGICAL :: lMPI, GlobalConv

nXsGenIterMax =  ItrCntl%DcplXsGenCntl%nXsGenIterMax
iRefPln = ItrCntl%DcplXsGenCntl%iRefPln
iRefTemp = ItrCntl%DcplXsGenCntl%iRefTemp

#ifdef MPI_ENV
Comm = DcplPE%MPI_DCPLMASTER_COMM
myrank = DcplPE%myrank; nproc = DcplPE%nproc
GlobalConv = .FALSE.; lMPI = DcplPE%lDcplParallel
Eigv = 1._8
DO iter = 1, nXsGenIterMax
  IF(lMOC) CALL Idle_DcplMoc(Eigv, ng, GroupInfo, ItrCntl, DcplPE)
  CALL Idle_DcplCmfd(EIgv, ItrCntl, DcplPE)
  GlobalConv = MPIConvChk(ItrCntl%lConv, nproc, comm)
  IF(GlobalConv) EXIT
ENDDO
#endif
END SUBROUTINE

SUBROUTINE Idle_DcplMOC(Eigv, ng, GroupInfo, ItrCntl, PE)
USE PARAM
USE TYPEDEF,     ONLY : GroupInfo_Type,       PE_Type
USE itrcntl_mod, ONLY : ItrCntl_TYPE
USE Files,        ONLY : io8
USE ioutil,      ONLY : message
IMPLICIT NONE
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_Type) :: PE
TYPE(ItrCntl_Type) :: ItrCntl
INTEGER :: ng
REAL :: Eigv

INTEGER :: iter, jsweep, iz, ig
INTEGER :: nGroupInfo, GrpBeg, GrpEnd
INTEGER :: myrank, Comm

LOGICAL :: lDmesg
!
lDmesg = .TRUE.
IF(ng .gt. 10) lDmesg = FALSE
iz = PE%myzb
!MPI
myrank = PE%myrank; comm = PE%MPI_DCPLMASTER_COMM
!
nGroupInfo = 2
IF(.NOT. GroupInfo%lUpScat) nGroupInfo = 1
IF(myrank .EQ. 0) THEN
  WRITE(mesg, '(a26)') 'Performing Ray Tracing ...'
  CALL message(io8, TRUE, TRUE, mesg)
ENDIF

DO iter = 1, 1
  DO jsweep = 1, nGroupInfo
    GrpBeg = 1; GrpEnd = ng
    IF(jsweep .GT. 1) THEN
      GrpBeg = GroupInfo%UpScatRange(1); GrpEnd = GroupInfo%UpScatRange(2)
    ENDIF
    DO ig = GrpBeg, GrpEnd
      IF(lDmesg .or. MOD(iG, 10) .eq. 0 .OR. ig .EQ. nG) THEN
        CALL PrintDcplRTLog(ig, 0._8, myrank, PE%nproc, Comm)
      ENDIF
    ENDDO
  ENDDO
ENDDO

CALL PrintDcplMOCLog(ItrCntl%MocIt, iz, EigV, 0._8, 0._8, 0._8, myrank, PE%nproc, Comm)
END SUBROUTINE

SUBROUTINE Idle_DcplCmfd(Eigv, ItrCntl, PE)
USE PARAM
USE TYPEDEF,       ONLY : PE_TYPE
USE ItrCntl_mod,   ONLY : ItrCntl_Type
USE FILES,         ONLY : io8
USE IOUTIL,        ONLY : message
IMPLICIT NONE
REAL :: EigV
TYPE(ItrCntl_Type) :: ItrCntl
TYPE(PE_TYPE) :: PE
CALL PrintDcplCmfdLog(ItrCntl%CmfdIt, PE%myzb, eigv, 0._8, 0._8, &
                      PE%myrank, PE%nProc, PE%MPI_DCPLMASTER_COMM)
END SUBROUTINE

END MODULE
