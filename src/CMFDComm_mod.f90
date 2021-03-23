MODULE CMFDComm_Mod

USE TYPEDEF,     ONLY : PE_TYPE, PinXs_Type, GroupInfo_Type, AxFlx_Type
USE MPICOMM_Mod, ONLY : NonBlockSend, NonBlockRecv, ChkNonBlockComm, MPI_SYNC,        &
                        SendRecv
!USE MPI
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER, PARAMETER :: ngmax = 1000, nbdmax = 6, nCommMax = 5000
INTEGER, PRIVATE :: nAxMax, nRadMax, nRadCommMax, ng0
TYPE PinXSComm_Type
  !INTEGER, DIMENSION(nbdmax) :: DTil, DHat
  INTEGER, DIMENSION(nbdmax) :: AxDtil, AxDhat
  INTEGER :: XST, XSTR, XSNF, CHI
  INTEGER :: Xss(ngmax), XssRange(2, ngmax)
  INTEGER :: PHI, LKG
  INTEGER :: ng, nbd
  INTEGER :: nTypedat
  INTEGER :: tag
END TYPE

TYPE XsQueueInfo_Type
  TYPE(PinXs_Type), POINTER ::  SendPinXs(:, :), RecvPinXs(:, :)
  REAL, POINTER :: SendPhi(:, :, :), SendLkg(:, :, :)
  REAL, POINTER :: RecvPhi(:, :, :), RecvLkg(:, :, :)
!  REAL, POINTER :: SendDhat(:, :, :, :), SendDTil(:, :, :, :)
!  REAL, POINTER :: RecvDhat(:, :, :, :), RecvDTil(:, :, :, :)
  LOGICAL :: lSend = .FALSE.
  LOGICAL :: lRecv = .FALSE.
  LOGICAL :: n1, n2, n, Dest
  REAL, POINTER :: SendBuf(:, :), RecvBuf(:, :)
END TYPE

TYPE AxNqueueInfo_Type
  REAL, POINTER :: SendDhat(:, :, :, :), SendDTil(:, :, :, :)
  REAL, POINTER :: RecvDhat(:, :, :, :), RecvDTil(:, :, :, :)
  LOGICAL :: lSend = .FALSE.
  LOGICAL :: lRecv = .FALSE.
  LOGICAL :: n1, n2, n, Dest
  REAL, POINTER :: SendBuf(:, :), RecvBuf(:, :)
END TYPE

TYPE CommPinXS_Type   !Axial Nodal Solver Recv and Send Data
  LOGICAL :: lAlloc = .FALSE.
  REAL, POINTER :: SendPhi(:, :, :), SendLkg(:, :, :)
  REAL, POINTER :: RecvPhi(:, :, :), RecvLkg(:, :, :)
  TYPE(PinXS_Type), POINTER :: SendPinXS(:, :), RecvPinXS(:, :)
  INTEGER :: SendDatRange(2,100), RecvDatRange(2)
END TYPE

INTEGER, PARAMETER :: Send = 1
INTEGER, PARAMETER :: Recv = 2

TYPE(XsQueueInfo_Type), PRIVATE :: XsQueueInfo(nCommMax)
TYPE(AxNQueueInfo_Type), PRIVATE :: AxNqueueInfo(nCommMax)

INTEGER, PRIVATE :: XsQueue(nCommMax), AxNqueue(nCommMax)
INTEGER, PRIVATE :: nXsQueue, nXsSend, nXsRecv
INTEGER, PRIVATE :: nAxNQueue, nAxNSend, nAxNRecv
LOGICAL, PRIVATE :: lInit = .FALSE.

!
INTEGER, PRIVATE :: MpiCommOrder(200, 0:200)
INTEGER, PRIVATE :: nCommStep
REAL, POINTER, SAVE, PRIVATE :: RadVar(:, :, :), AxVar(:, :, :), RadAxCommBuf(:)

TYPE(PinXSComm_Type), PRIVATE :: PinXsComm
INTEGER, PRIVATE :: nTypedat

INTEGER, PRIVATE :: myNxyBeg, myNxyEnd, nxy
INTEGER, PRIVATE :: nz, nzfm, myzb, myze, myzbf, myzef, nSubPlane
INTERFACE ASSIGNMENT (=)
  MODULE PROCEDURE PinXS_Assignment
END INTERFACE

CONTAINS

SUBROUTINE PinXs_Assignment(PinXS1, PinXS2)
TYPE(PinXS_TYPE), INTENT(INOUT) :: PinXS1
TYPE(PinXS_TYPE), INTENT(in) :: PinXS2
INTEGER :: ig, ib, ie
PinXS1%XSD(1:ng0) = PinXS2%XSD(1:ng0); PinXS1%XSD2(1:ng0) = PinXS2%XSD2(1:ng0)
PinXS1%XST(1:ng0) = PinXS2%XST(1:ng0)
PinXS1%XSTR(1:ng0) = PinXS2%XSTR(1:ng0); PinXS1%XSR(1:ng0) = PinXS2%XSR(1:ng0)
PinXS1%XSNF(1:ng0) = PinXS2%XSNF(1:ng0); PinXS1%XSKF(1:ng0) = PinXS2%XSKF(1:ng0)
PinXS1%CHI(1:ng0) = PinXS2%CHI(1:ng0)
!
DO ig = 1, ng0
  PinXS1%Xss(ig)%withinGroupScat = PinXS2%Xss(ig)%withinGroupScat
  ib = PinXS2%Xss(ig)%ib; ie = PinXS2%Xss(ig)%ie
  PinXS1%xss(ig)%from(ib:ie) = PinXS2%xss(ig)%from(ib:ie)
ENDDO
END SUBROUTINE

FUNCTION GetnCommStep(n, nmax)
INTEGER :: n, nmax, GetnCommStep
INTEGER :: i, j
i = n / nmax
j = mod(n, nmax)
IF(J .EQ. 0) THEN
 GetnCommStep = i
ELSE
 GetnCommStep = i + 1
ENDIF
END FUNCTION


SUBROUTINE InitAxNComm(ng, n1max, n2max)
USE ALLOCS
USE IOUTIL,    ONLY : TERMINATE
IMPLICIT NONE
INTEGER :: ng, n1max, n2max
INTEGER :: n 
INTEGER :: i 
n = n2max
IF(n2max .GT. nCommMax) CALL TERMINATE('Exceed Maximum number of MPI Nonblock Communication')
DO i = 1, n
  CALL Dmalloc(AxNQueueInfo(i)%SendBuf, 5*ng, n1max)
ENDDO

END SUBROUTINE

SUBROUTINE InitAxNodalComm(GroupInfo, ng, nbd, SendRecvmax, PE, lNonBlock)
USE ALLOCS
USE IOUTIL,    ONLY : TERMINATE
IMPLICIT NONE
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_TYPE) :: PE
INTEGER :: ng, nbd, SendRecvMax
LOGICAL :: lNonBlock

INTEGER :: nproc
!INTEGER :: myNxyBeg, myNxyEnd, myzb, myze, nz
INTEGER :: nradpins(0:500), nAxpln(0:500), n1, n2
INTEGER :: n1max, n2max
INTEGER :: i


CALL MakePinXsCommType(GroupInfo, ng, nbd)

myzb =PE%myzb; myze=PE%myze; 
myzbf = PE%myzbf; myzef = PE%myzef
nz = PE%nz; nzfm = PE%nzfm
nSubPlane = nzfm / nz
myNxyBeg = PE%myNxyBeg; myNxyEnd = PE%myNxyEnd

ng0 = ng
nxy = PE%nxy; nproc = PE%nCMFDproc
nRadPins(0:nproc-1) = PE%RadDomRange(2, 0:nproc-1) - PE%RadDomRange(1, 0:nproc-1) + 1
nAxPln(0:nproc-1) = PE%AxDomRange(2, 0:nproc-1) - PE%AxDomRange(1, 0:nproc-1) + 1
!Find Maximum assigned radial mesh
n1 = maxval(nRadPins(0:nproc-1))       !Maximum of Radial mesh
n2 = maxval(nAxPln(0:nproc-1))         !Maximun of Axial Mesh


n2max = nproc * 2 
n1max = n1 * n2
n1max = Min(n1max, SendRecvMax * n2)
SendRecvMax = n1max / n2 

nRadCommMax = SendRecvMax; 
nRadMax = n1; nAxMax = n2 !Assign Global variable inside of this module 
IF(n2max .GT. nCommMax) CALL TERMINATE('Exceed Maximum number of MPI Nonblock Communication')

IF(.NOT. lNonBlock) THEN
  CALL SetMpiCommOrder(MpiCommOrder, nCommStep, PE%nCmfdProc)
  CALL Dmalloc0(AxVar, 1, 4 * ng, myNxyBeg, myNxyEnd, 1, nzfm)
  CALL Dmalloc0(RadVar, 1, 4 * ng, 1, nxy, myzbf, myzef)
  CALL Dmalloc(RadAxCommBuf, 4 * ng * nRadMax*nAxMax * nSubPlane)
  !PRINT *, 'Alloc Buf',  4 * ng * nRadMax*nAxMax
ELSE
  DO i = 1, n2max
    CALL Dmalloc(XsQueueInfo(i)%SendBuf, nTypedat, n1max)
    CALL Dmalloc(XsQueueInfo(i)%RecvBuf, nTypedat, n1max)
  ENDDO
ENDIF
nXsQueue = 0
nXsSend = 0
nXsRecv = 0
lInit = .TRUE.

END SUBROUTINE

SUBROUTINE MakePinXsCommType(GroupInfo, ng, nbd)
TYPE(GroupInfo_Type) :: GroupInfo
INTEGER :: ng, nbd
INTEGER :: ig, ibd, ib, ie
ng = GroupInfo%ng

nTypedat = 0

PinXsComm%ng = ng; PinXsComm%nbd = nbd
PinXsComm%XST = nTypeDat + 1; nTypeDat = nTypeDat + ng
PinXsComm%XSTR = nTypeDat + 1; nTypeDat = nTypeDat + ng
PinXsComm%XSNF = nTypeDat + 1; nTypeDat = nTypeDat + ng
PinXsComm%CHI = nTypeDat + 1; nTypeDat = nTypeDat + ng
!Scattering Term
DO ig = 1, ng
  ib = GroupInfo%InScatRange(1, ig); ie = GroupInfo%InScatRange(2, ig)
  PinXsComm%XssRange(1:2, ig) = (/ib, ie/)
  PinXsComm%Xss(ig) = nTypeDat + 1 ; nTypeDat = nTypeDat + (ie-ib+1)
ENDDO
PinXsComm%phi = nTypeDat + 1 ; nTypeDat = nTypeDat + (ng)
PinXsComm%LKG = nTypeDat + 1 ; nTypeDat = nTypeDat + (ng)
PinXsComm%nTypeDat = nTypeDat
END SUBROUTINE

SUBROUTINE UpdtArrayLoc(LocArray, n1, nsize, nTypeDat0)
INTEGER :: LocArray(n1)
INTEGER :: n1, n2, nsize, nTypeDat0
INTEGER :: i
DO i = 1, n1
  LocArray(i) = nTypeDat0 + 1
  nTypeDat0 = nTypeDat0 + nsize
ENDDO

END SUBROUTINE

SUBROUTINE WritePinXSBuf(SendBuf, PinXS, npin, npln, ng)
!SUBROUTINE WritePinXSBuf(SendBuf, PinXS, Phi, LKG, npin, npln, ng)
TYPE(PinXS_TYPE) :: PinXS(npin, npln)
REAL :: SendBuf(nTypeDat, npin * npln)
!REAL :: PHI(npin, npln, ng)
!REAL :: LKG(npin, npln, ng)
INTEGER :: npin, npln, ng
INTEGER :: i, j, k
k = 0
DO i = 1, npln
  DO j = 1, npin
    k = k + 1 
    !CALL PackingPinXS(SendBuf(:, k), PinXS(j, i), Phi(j, i, :), Lkg(j, i, :), ng)
    CALL PackingPinXS(SendBuf(:, k), PinXS(j, i), ng)
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE PackingPinXs(PackingDat, PinXS,  ng)
!SUBROUTINE PackingPinXs(PackingDat, PinXS, PHI, LKG, ng)
TYPE(PinXs_Type) :: PinXS
REAL :: PackingDat(nTypeDat)
!REAL :: PHI(NG), LKG(NG)
INTEGER :: ng
REAL :: Temp(ngmax)

INTEGER :: nbd
INTEGER :: ib, ie, ndat
INTEGER :: igb, ige, ig

nbd = PinXsComm%nbd
ib = PinXsComm%XST; ie = ib + ng - 1
PackingDat(ib:ie) = PinXs%XST(1:ng)

ib = PinXsComm%XSTR; ie = ib + ng - 1
PackingDat(ib:ie) = PinXs%XSTR(1:ng)

ib = PinXsComm%XSNF; ie = ib + ng - 1
PackingDat(ib:ie) = PinXs%XSNF(1:ng)

ib = PinXsComm%CHI; ie = ib + ng - 1
PackingDat(ib:ie) = PinXs%CHI(1:ng)

DO ig = 1, ng
  igb = PinXsComm%XssRange(1, ig); ige = PinXsComm%XssRange(2, ig)
  ndat = ige - igb + 1 
  ib = PinXsComm%Xss(ig); ie = ib + ndat - 1
  Temp(igb:ige) = PinXs%Xss(ig)%from(igb:ige); Temp(ig) = PinXS%Xss(ig)%WithInGroupScat
  PackingDat(ib:ie) = Temp(igb:ige)
ENDDO

ib = PinXsComm%PHI; ie = ib + ng - 1
!PackingDat(ib:ie) = PHI(1:ng)
ib = PinXsComm%LKG; ie = ib + ng - 1
!PackingDat(ib:ie) = LKG(1:ng)
END SUBROUTINE

!SUBROUTINE ReadPinXSbuf(RecvBuf, PinXS, Phi, LKG, npin, npln, ng)
SUBROUTINE ReadPinXSbuf(RecvBuf, PinXS, npin, npln, ng)
TYPE(PinXS_TYPE) :: PinXS(npin, npln)
REAL :: RecvBuf(nTypeDat, npin * npln)
!REAL :: PHI(npin, npln, ng)
!REAL :: LKG(npin, npln, ng)
INTEGER :: npin, npln, ng
INTEGER :: i, j, k
k = 0
DO i = 1, npln
  DO j = 1, npin
    k = k + 1
   ! CALL UnPackingPinXS(RecvBuf(:, k), PinXS(j, i), Phi(j, i, :), Lkg(j, i, :), ng)
    CALL UnPackingPinXS(RecvBuf(:, k), PinXS(j, i), ng)
  ENDDO
ENDDO

END SUBROUTINE

!SUBROUTINE UnPackingPinXs(PackingDat, PinXS, PHI, LKG, ng)
SUBROUTINE UnPackingPinXs(PackingDat, PinXS, ng)
USE CNTL,        only : nTracerCntl
TYPE(PinXs_Type) :: PinXS
REAL :: PackingDat(nTypeDat)
!REAL :: PHI(NG), LKG(NG)
INTEGER :: ng

REAL :: Temp(ngmax)
INTEGER :: nbd
INTEGER :: ib, ie, ndat
INTEGER :: igb, ige, ig, ig2
nbd = PinXsComm%nbd
ib = PinXsComm%XST; ie = ib + ng - 1
PinXs%XST(1:ng) = PackingDat(ib:ie)
ib = PinXsComm%XSTR; ie = ib + ng - 1
PinXs%XSTR(1:ng) = PackingDat(ib:ie)

ib = PinXsComm%XSNF; ie = ib + ng - 1
PinXs%XSNF(1:ng) = PackingDat(ib:ie)

ib = PinXsComm%CHI; ie = ib + ng - 1
PinXs%CHI(1:ng) = PackingDat(ib:ie)

DO ig = 1, ng
  igb = PinXsComm%XssRange(1, ig); ige = PinXsComm%XssRange(2, ig)
  ndat = ige - igb + 1 
  ib = PinXsComm%Xss(ig); ie = ib + ndat - 1
  Temp(igb:ige) = PackingDat(ib:ie)
  PinXS%Xss(ig)%WithInGroupScat = Temp(ig); Temp(ig) = 0
  PinXs%Xss(ig)%from(igb:ige) = Temp(igb:ige)
ENDDO

ib = PinXsComm%PHI; ie = ib + ng - 1
!PHI(1:ng) = PackingDat(ib:ie)
ib = PinXsComm%LKG; ie = ib + ng - 1
!LKG(1:ng) = PackingDat(ib:ie)

!Fill Out the XS
IF( nTracerCntl%libtyp .EQ. 2 )THEN
    DO ig = 1, ng
      PinXS%XSD(ig) = 1._8 / (3._8 * PinXS%XSTR(ig))
      PinXS%XSD2(ig) = 3._8 / (7._8 * PinXS%XST(ig))
      !PinXS%XSD2(ig) = 3._8 / (7._8 * PinXS%XSTR(ig))
      !- removal independent to transport correction      
      PinXS%XSR(ig) = PinXS%XSTR(ig) - PinXS%Xss(ig)%WithInGroupScat
    ENDDO
ELSE
    DO ig = 1, ng
      PinXS%XSD(ig) = 1._8 / (3._8 * PinXS%XSTR(ig))
      PinXS%XSD2(ig) = 3._8 / (7._8 * PinXS%XST(ig))
      !PinXS%XSD2(ig) = 3._8 / (7._8 * PinXS%XSTR(ig))
      !- removal independent to transport correction      
      PinXS%XSR(ig) = PinXS%XST(ig) - PinXS%Xss(ig)%WithInGroupScat
    ENDDO
ENDIF
    
END SUBROUTINE

SUBROUTINE SetMpiCommOrder(CommOrderList, nStep, nproc)
IMPLICIT NONE

TYPE MpiCommOrder_Type
  LOGICAL :: lComm(200)
  LOGICAL :: lPartner, lFinish
  INTEGER :: CommDestList(200)
  INTEGER :: CommDest(200)
  INTEGER :: nComm
END TYPE

INTEGER :: CommOrderList(200, 0:200)
INTEGER :: nproc, nstep
TYPE(MpiCommOrder_Type) :: MpiCommVar(0:nproc-1)

INTEGER :: i, j, k, m1, m2
INTEGER :: ProcStart, jrank, jdest
INTEGER :: ncomm1, ncomm2

LOGICAL :: lExit 
#ifdef ChkCompleteComm
LOGICAL :: ChkComm(0:nproc-1, 0:nproc-1)
#endif
!initialize the variable

DO i = 0, nproc - 1
  MpiCommVar(i)%lComm(1:nproc-1) = .FALSE.
  MpiCommVar(i)%CommDest(1:nproc-1) = -2
  k = 0
  DO j = 0, nproc - 1
    IF(j .EQ. i) CYCLE; k = k + 1
    MpiCommVar(i)%CommDestList(k) = j
  ENDDO
  MpiCommVar(i)%nComm = 0
ENDDO
DO i = 0, 10000
  ProcStart = mod(i, nproc)
  lExit = .TRUE.
  DO j = 0, nproc - 1
    MpiCommVar(j)%lPartner = .FALSE.
    MpiCommVar(j)%lFinish = .TRUE.
    DO k = 1, nproc-1
      MpiCommVar(j)%lFinish = MpiCommVar(j)%lFinish .AND. MpiCommVar(j)%lComm(k)
    ENDDO
  ENDDO
  
  DO j = ProcStart, Procstart + nproc - 1
    jrank = Mod(j, nproc)
    IF(MpiCommVar(jrank)%lFinish) CYCLE
    lExit = .FALSE.
    IF(MpiCommVar(jrank)%lPartner) CYCLE
    DO m1 = 1, nproc - 1
      IF(MpiCommVar(jrank)%lComm(m1)) CYCLE  !Cycle if this communication were done before
      jdest = MpiCommVar(jrank)%CommDestList(m1)
      IF(MpiCommVar(jdest)%lPartner) CYCLE !Pass if destination processor is communicating with other
      MpiCommVar(jdest)%lPartner = .TRUE.
      MpiCommVar(jrank)%lPartner = .TRUE.
      DO m2 = 1, nproc-1
        IF(MpiCommVar(jdest)%CommDestList(m2) .EQ. jrank) EXIT
      ENDDO
      !Set Variable
      MpiCommVar(jrank)%lComm(m1) = .TRUE.
      MpiCommVar(jdest)%lComm(m2) = .TRUE.
      
      MpiCOmmVar(jrank)%ncomm = MpiCOmmVar(jrank)%ncomm + 1
      MpiCOmmVar(jdest)%ncomm = MpiCOmmVar(jdest)%ncomm + 1
      ncomm1 = MpiCOmmVar(jrank)%ncomm; ncomm2 = MpiCOmmVar(jdest)%ncomm
      MpiCOmmVar(jrank)%CommDest(ncomm1) = jdest
      MpiCOmmVar(jdest)%CommDest(ncomm2) = jrank
      EXIT
    ENDDO
    !in case of not find the partner)
    IF( m1 .EQ. nproc) THEN
      MpiCommVar(jrank)%nComm = MpiCommVar(jrank)%nComm + 1
      ncomm1 = MpiCommVar(jrank)%nComm
      MpiCommVar(jrank)%CommDest(ncomm1) = -1
    ENDIF
  ENDDO
  CONTINUE
  IF(lExit) EXIT
ENDDO
nstep = i

DO j = 0, nproc - 1
  nComm1 = MpiCommVar(j)%nComm
  CommOrderList(1:nstep, j) = -1
  DO i = 1, nComm1
    CommOrderList(i, j) = MpiCommVar(j)%CommDest(i)
  ENDDO
ENDDO

#ifdef ChkCompleteComm
ChkComm = .FALSE.
DO i = 0, nProc - 1
  ChkComm(i, i) = .TRUE.
  ncomm1 = MpiCommVar(i)%nComm 
  DO j = 1, ncomm1
    jdest = MpiCOmmVar(i)%CommDest(j)
    IF(jDest .EQ. -1) CYCLE
    ChkCOmm(JDest, i) = .TRUE.
  ENDDO
ENDDO
#endif
CONTINUE
END SUBROUTINE

SUBROUTINE SendPinXS(PinXS, PHI, LKG, n1, n2, ng, Dest, Comm)
USE PARAM
USE TYPEDEF, ONLY : PinXS_Type
IMPLICIT NONE
TYPE(PinXS_Type), TARGET :: PinXS(n1, n2)
REAL, TARGET :: PHI(n1, n2, ng), LKG(n1, n2, ng)
INTEGER :: Dest, Comm
INTEGER :: n1, n2

INTEGER :: ng
INTEGER :: i, j, k, ierr

nXsQueue = nXsQueue + 1
nXsSend = nXsSend + 1
XsQueueInfo(nXsQueue)%lSend = TRUE
XsQueueInfo(nXsQueue)%lRecv = FALSE
k = 0
DO j = 1, n2
  DO i = 1, n1
    k = k + 1
!    CALL PackingPinXS(XsQueueInfo(nXsQueue)%SendBuf(:,k), PinXs(i, j), PHI(i,j, :),   &
!                      LKG(i, j, :), ng)
    CALL PackingPinXS(XsQueueInfo(nXsQueue)%SendBuf(:,k), PinXs(i, j), ng)                      
  ENDDO
ENDDO
XsQueueInfo(nXsQueue)%SendPinXS => PinXS
XsQueueInfo(nXsQueue)%SendPHI => PHI
XsQueueInfo(nXsQueue)%SendLKG => LKG
XsQueueInfo(nXsQueue)%Dest = Dest
XsQueueInfo(nXsQueue)%n1 = n1
XsQueueInfo(nXsQueue)%n2 = n2

CALL NonBlockSend(XsQueueInfo(nXsQueue)%SendBuf(:, 1:k), nTypeDat, k, Dest, Comm,       &
                  XsQueue(nXsQueue))

END SUBROUTINE

SUBROUTINE RecvPinXS(PinXS, PHI, LKG, n1, n2, ng, Dest, Comm)
USE PARAM
USE TYPEDEF, ONLY : PinXS_Type
IMPLICIT NONE
TYPE(PinXS_Type), TARGET :: PinXS(n1, n2)
REAL, TARGET :: PHI(n1, n2, ng), LKG(n1, n2, ng)
INTEGER :: Dest, Comm
INTEGER :: n1, n2

INTEGER :: ng
INTEGER :: i, j, k, ierr

k = n1*n2
nXsQueue = nXsQueue + 1
nXsRecv = nXsRecv + 1
XsQueueInfo(nXsQueue)%lSend = FALSE
XsQueueInfo(nXsQueue)%lRecv = TRUE
XsQueueInfo(nXsQueue)%RecvPinXS => PinXS
XsQueueInfo(nXsQueue)%RecvPHI => PHI
XsQueueInfo(nXsQueue)%RecvLkg => LKG
XsQueueInfo(nXsQueue)%Dest = Dest
XsQueueInfo(nXsQueue)%n1 = n1
XsQueueInfo(nXsQueue)%n2 = n2

CALL NonBlockRecv(XsQueueInfo(nXsQueue)%RecvBuf(:, 1:k), nTypeDat, k, Dest, Comm,       &
                  XsQueue(nXsQueue))
END SUBROUTINE


SUBROUTINE CommGlobalPinXsDat(CommDat, ng, GroupInfo, PE)
USE ALLOCS
USE BasicOperation,  ONLY : CP_VA
IMPLICIT NONE
TYPE(CommPinXS_Type) :: CommDat 
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_TYPE) :: PE
INTEGER :: ng

REAL, POINTER :: SendPhiDat(:, :, :), SendLkgDat(:, :, :)
REAL, POINTER :: RecvPhiDat(:, :, :), RecvLkgDat(:, :, :)
TYPE(PinXS_Type), POINTER :: SendPinXSDat(:, :), RecvPinXSDat(:, :)

REAL, POINTER, SAVE, DIMENSION(:, :, :) :: LKG, PHI, LKG0, PHI0

REAL, POINTER, SAVE :: SendBuf(:, :), RecvBuf(:, :)

!INTEGER :: nradpins(0:500), nAxpln(0:500), n1, n2
!INTEGER :: nxy, nz, myzb, myze, myNxyBeg, myNxyEnd
INTEGER :: nproc, myrank, nComm

INTEGER :: i, j, k
INTEGER :: ipartner, nsendrecvstep
INTEGER :: ixybeg, ixyend, izbeg, izend, istep, nsendpin, npln, npin, ndat
INTEGER :: ixy, iz, ixy1, ixy2
INTEGER :: COMM

LOGICAL, SAVE :: lFirst
DATA lFirst /.TRUE./



IF(lFirst) THEN
  CALL Dmalloc(SendBuf, nTypeDat, nRadCommMax * nAxMax)
  CALL Dmalloc(RecvBuf, nTypeDat, nRadCommMax * nAxMax)
  !IF(.NOT. lSetCommOrder) CALL SetMpiCommOrder(MpiCommOrder, nCommStep, PE%nCmfdProc)
  lFirst = .FALSE.
ENDIF

Comm = PE%MPI_CMFD_COMM
nproc = PE%nCmfdProc; myrank = PE%myrank
!nRadCommMax = SendRecvMax; nAxMax = n2 !Assign Global variable inside of this module 
!SendPhiDat => CommDat%SendPhi; RecvPhiDat => CommDat%RecvPhi
!SendLkgDat => CommDat%SendLkg; RecvLkgDat => CommDat%RecvLkg
SendPinXSDat => CommDat%SendPinXS; RecvPinXsDat => CommDat%RecvPinXS

!Copy own domain data
DO iz = myzb,myze
  DO ixy = myNxyBeg, myNxyEnd
    RecvPinXSDat(ixy, iz) = SendPinXSDat(ixy, iz)
  ENDDO
ENDDO

!CALL CP_VA(RecvPhiDat(myNxyBeg:myNxyEnd, myzb:myze, 1:ng), SendPhiDat(myNxyBeg:myNxyEnd, myzb:myze, 1:ng), &
!           myNxyEnd - myNxyBeg + 1, myze - myzb + 1, ng)
!CALL CP_VA(RecvLkgDat(myNxyBeg:myNxyEnd, myzb:myze, 1:ng), SendLkgDat(myNxyBeg:myNxyEnd, myzb:myze, 1:ng), &
!           myNxyEnd - myNxyBeg + 1, myze - myzb + 1, ng)

DO istep = 1, nCommStep
  ipartner = MpiCommOrder(istep, myrank)
  IF(ipartner .LT. 0) THEN
    CALL MPI_SYNC(comm)
    CYCLE
  ENDIF
  IF(myrank .LT. ipartner) THEN  !Send -> Recv
    ixybeg = PE%RadDomRange(1, ipartner); ixyend = PE%RadDomRange(2, ipartner)
    izbeg = PE%myzb;izend = PE%myze
    npin = ixyend - ixybeg + 1; npln = izend - izbeg + 1
    nSendRecvStep = GetnCommStep(npin, nRadCommMax) 
    DO j = 1, nSendRecvStep
      ixy1 = ixybeg + nRadCommMax * (j - 1)
      ixy2 = ixy1 + nRadCommMax - 1; ixy2 = min(ixy2, ixyend)
      npin = ixy2 - ixy1 + 1; ndat = npin * npln
      !Packing the data
      CALL WritePinXSBuf(SendBuf(:, 1:ndat), SendPinXSDat(ixy1:ixy2, izbeg:izend),                        & 
                         npin, npln, ng)
!      CALL WritePinXSBuf(SendBuf(:, 1:ndat), SendPinXSDat(ixy1:ixy2, izbeg:izend),                        & 
!                         SendPhiDat(ixy1:ixy2, izbeg:izend, :), SendLKGDat(ixy1:ixy2, izbeg:izend, :),  &
!                         npin, npln, ng)
      !Send Data
      CALL SendRecv(SendBuf(:, 1:ndat), nTypeDat, ndat, ipartner, .TRUE., COMM)
    ENDDO 
  ELSE  !Recv -> Send
    ixybeg = PE%RadDomRange(1, myrank); ixyend = PE%RadDomRange(2, myrank)
    izbeg = PE%AxDomRange(1, ipartner); izend = PE%AxDomRange(2, ipartner)
    npin = ixyend - ixybeg + 1; npln = izend - izbeg + 1
    nSendRecvStep = GetnCommStep(npin, nRadCommMax) 
    DO j = 1, nSendRecvStep
      ixy1 = ixybeg + nRadCommMax * (j - 1)
      ixy2 = ixy1 + nRadCommMax - 1; ixy2 = min(ixy2, ixyend)
      npin = ixy2 - ixy1 + 1; ndat = npin * npln
      CALL SendRecv(RecvBuf(:, 1:ndat), nTypeDat, ndat, ipartner, .FALSE., COMM)
!      CALL ReadPinXSBuf(RecvBuf(:, 1:ndat), RecvPinXSDat(ixy1:ixy2, izbeg:izend),                         & 
!                        RecvPhiDat(ixy1:ixy2, izbeg:izend, :), RecvLKGDat(ixy1:ixy2, izbeg:izend, :),  &
!                        npin, npln, ng) 
      CALL ReadPinXSBuf(RecvBuf(:, 1:ndat), RecvPinXSDat(ixy1:ixy2, izbeg:izend),                         & 
                        npin, npln, ng) 
    ENDDO
  ENDIF
  IF(myrank .LT. ipartner) THEN 
    ixybeg = PE%RadDomRange(1, myrank); ixyend = PE%RadDomRange(2, myrank)
    izbeg = PE%AxDomRange(1, ipartner); izend = PE%AxDomRange(2, ipartner)
    npin = ixyend - ixybeg + 1; npln = izend - izbeg + 1
    nSendRecvStep = GetnCommStep(npin, nRadCommMax) 
    DO j = 1, nSendRecvStep
      ixy1 = ixybeg + nRadCommMax * (j - 1)
      ixy2 = ixy1 + nRadCommMax - 1; ixy2 = min(ixy2, ixyend)
      npin = ixy2 - ixy1 + 1; ndat = npin * npln
      CALL SendRecv(RecvBuf(:, 1:ndat), nTypeDat, ndat, ipartner, .FALSE., COMM)
!      CALL ReadPinXSBuf(RecvBuf(:, 1:ndat), RecvPinXSDat(ixy1:ixy2, izbeg:izend),                         & 
!                        RecvPhiDat(ixy1:ixy2, izbeg:izend, :), RecvLKGDat(ixy1:ixy2, izbeg:izend, :),  &
!                        npin, npln, ng)
      CALL ReadPinXSBuf(RecvBuf(:, 1:ndat), RecvPinXSDat(ixy1:ixy2, izbeg:izend),                         & 
                        npin, npln, ng)
    ENDDO
  ELSE !Write
    ixybeg = PE%RadDomRange(1, ipartner); ixyend = PE%RadDomRange(2, ipartner)
    izbeg = PE%myzb;izend = PE%myze
    npin = ixyend - ixybeg + 1; npln = izend - izbeg + 1
    nSendRecvStep = GetnCommStep(npin, nRadCommMax) 
    DO j = 1, nSendRecvStep
      ixy1 = ixybeg + nRadCommMax * (j - 1)
      ixy2 = ixy1 + nRadCommMax - 1; ixy2 = min(ixy2, ixyend)
      npin = ixy2 - ixy1 + 1; ndat = npin * npln
      !Packing the data
      CALL WritePinXSBuf(SendBuf(:, 1:ndat), SendPinXSDat(ixy1:ixy2, izbeg:izend),                        & 
                         npin, npln, ng)
!      CALL WritePinXSBuf(SendBuf(:, 1:ndat), SendPinXSDat(ixy1:ixy2, izbeg:izend),                        & 
!                         SendPhiDat(ixy1:ixy2, izbeg:izend, :), SendLKGDat(ixy1:ixy2, izbeg:izend, :),  &
!                         npin, npln, ng)
      !Send Data
      CALL SendRecv(SendBuf(:, 1:ndat), nTypeDat, ndat, ipartner, .TRUE., COMM)      
    ENDDO 
  ENDIF
  CALL MPI_SYNC(comm)
ENDDO

!NULLIFY(SendPhiDat, RecvPhiDat)
!NULLIFY(SendLkgDat, RecvLkgDat)
NULLIFY(SendPinXsDat, RecvPinXsDat)

END SUBROUTINE

SUBROUTINE CommPinXsDat(CommDat, ng, GroupInfo, PE, SendRecv)

TYPE(CommPinXS_Type) :: CommDat 
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_TYPE) :: PE
INTEGER :: ng, SendRecv


REAL, POINTER :: SendPhiDat(:, :, :), SendLkgDat(:, :, :)
REAL, POINTER :: RecvPhiDat(:, :, :), RecvLkgDat(:, :, :)
TYPE(PinXS_Type), POINTER :: SendPinXSDat(:, :), RecvPinXSDat(:, :)

REAL, POINTER, SAVE, DIMENSION(:, :, :) :: LKG, PHI, LKG0, PHI0

!INTEGER :: nradpins(0:500), nAxpln(0:500), n1, n2
!INTEGER :: nxy, nz, myzb, myze, myNxyBeg, myNxyEnd
INTEGER :: nproc

INTEGER :: i, j, k
INTEGER :: ixybeg, ixyend, izbeg, izend, nsendpin, npln

INTEGER :: COMM


nproc = PE%nCmfdProc

Comm = PE%MPI_CMFD_COMM
CommDat%SendPhi => SendPhiDat; CommDat%RecvPhi => RecvPhiDat
CommDat%SendLkg => SendLkgDat; CommDat%RecvLkg => RecvLkgDat
CommDat%SendPinXS => SendPinXSDat; CommDat%RecvPinXS => RecvPinXsDat

IF(SendRecv .EQ. SEND) THEN
  DO i = 0, nproc-1
    ixybeg = CommDat%SendDatRange(1, i); ixyend =  CommDat%SendDatRange(2, i)
    npln = myze - myze + 1
    nsendpin = ixyend - ixybeg + 1
    CALL SendPinXS(SendPinXsDat(ixybeg:ixyend, myzb:myze), SendPhiDat(ixybeg:ixyend, myzb:myze, :),        &
                   SendLkgDat(ixybeg:ixyend, myzb:myze, :), nsendpin, npln, ng, i, Comm)
  ENDDO
  DO i = 0, nproc - 1
    izbeg = PE%AxDomRange(1, i); izend = PE%AxDomRange(2, i)
    npln = izend - izbeg + 1
    ixybeg = CommDat%RecvDatRange(1); ixyend =  CommDat%RecvDatRange(2)
    nsendpin = ixyend - ixybeg + 1
    CALL RecvPinXS(RecvPinXSDat(ixybeg:ixyend, myzb:myze), RecvPhiDat(ixybeg:ixyend, myzb:myze, :),         &
                   RecvLkgDat(ixybeg:ixyend, myzb:myze, :), nsendpin, npln, ng, i, Comm)
    !CALL RecvPinXS(PinXS0(1:1, iz1:iz2), PHI(1:ng, 1:1, iz1:iz2), LKG(1:ng, 1:1, iz1:iz2) , 1, iz2 - iz1 + 1, ng, i-1, COMM) 
  ENDDO
ELSE
  CALL FinalizePinXScomm()
ENDIF

NULLIFY(SendPhiDat, RecvPhiDat)
NULLIFY(SendLkgDat, RecvLkgDat)
NULLIFY(SendPinXsDat, RecvPinXsDat)
END SUBROUTINE

!Send radially distributed data to axial domain variable
SUBROUTINE CommRadDom2AxDom(ndatbase, PE)
IMPLICIT NONE
!TYPE(CommAxRad_Type) :: CommAxRad
TYPE(PE_TYPE) :: PE

INTEGER :: nproc, myrank, nComm

REAL, POINTER :: Buf(:)
!REAL, POINTER :: AxVar(:, :, :), RadVar(:, :, :)

INTEGER :: i, j, k
INTEGER :: ipartner, nsendrecvstep
INTEGER :: ixybeg, ixyend, izbeg, izend, istep, nsendpin, npln, npin
INTEGER :: ndatbase, nDatMax, ndat, ndat0
INTEGER :: ixy, iz, ixy1, ixy2
INTEGER :: COMM

!RadVar => CommAxRad%RadVar; AxVar => CommAxRad%AxVar;

nproc = PE%nCmfdProc; myrank = PE%myrank
Comm = PE%MPI_CMFD_COMM

!ndatbase = CommAxRad%npindat
nDatMax = nAxMax * nRadMax
Buf => RadAxCommBuf
DO istep = 1, nCommStep
  ipartner = MpiCommOrder(istep, myrank)
  IF(ipartner .LT. 0) CYCLE
  IF(myrank .LT. ipartner) THEN  !Send -> Recv
    ixybeg = PE%RadDomRange(1, ipartner); ixyend = PE%RadDomRange(2, ipartner)
    izbeg = PE%myzbf;izend = PE%myzef
    npin = ixyend - ixybeg + 1; npln = izend - izbeg + 1
    ndat = npin * npln; ndat0 = ndat * ndatBase
    CALL WriteAxRadCommBuf(RadVar(1:nDatBase, ixybeg:ixyend, izbeg:izend), Buf(1:ndat0), ndat0, nDatBase, npin, npln)
    CALL SendRecv(Buf(1:ndat0), ndat0, ipartner, .TRUE., COMM)
  ELSE  !Recv -> Send
    ixybeg = PE%RadDomRange(1, myrank); ixyend = PE%RadDomRange(2, myrank)
    !izbeg = PE%AxDomRange(1, ipartner); izend = PE%AxDomRange(2, ipartner)
    izbeg = PE%SubPlnDomRange(1, ipartner); izend = PE%SubPlnDomRange(2, ipartner)
    npin = ixyend - ixybeg + 1; npln = izend - izbeg + 1
    ndat = npin * npln; ndat0 = ndat * ndatBase
    CALL SendRecv(Buf(1:ndat0), ndat0, ipartner, .FALSE., COMM)
    CALL ReadAxRadCommBuf(AxVar(1:nDatBase, ixybeg:ixyend, izbeg:izend), Buf(1:ndat0), ndat0, nDatBase, npin, npln)
  ENDIF
  IF(myrank .LT. ipartner) THEN 
    ixybeg = PE%RadDomRange(1, myrank); ixyend = PE%RadDomRange(2, myrank)
    izbeg = PE%SubPlnDomRange(1, ipartner); izend = PE%SubPlnDomRange(2, ipartner)
    npin = ixyend - ixybeg + 1; npln = izend - izbeg + 1
    ndat = npin * npln; ndat0 = ndat * ndatBase
    CALL SendRecv(Buf(1:ndat0), ndat0, ipartner, .FALSE., COMM)
    CALL ReadAxRadCommBuf(AxVar(1:nDatBase, ixybeg:ixyend, izbeg:izend), Buf(1:ndat0), nDat0, nDatBase, npin, npln)
  ELSE !Write
    ixybeg = PE%RadDomRange(1, ipartner); ixyend = PE%RadDomRange(2, ipartner)
    izbeg = PE%myzbf;izend = PE%myzef
    npin = ixyend - ixybeg + 1; npln = izend - izbeg + 1
    ndat = npin * npln; ndat0 = ndat * ndatBase
    CALL WriteAxRadCommBuf(RadVar(1:nDatBase, ixybeg:ixyend, izbeg:izend), Buf(1:ndat0), ndat0, nDatBase, npin, npln)
    CALL SendRecv(Buf(1:ndat0), ndat0, ipartner, .TRUE., COMM)
  ENDIF
ENDDO

CALL MPI_SYNC(comm)

END SUBROUTINE

SUBROUTINE CommAxDom2RadDom(nDatBase, PE)
IMPLICIT NONE
TYPE(PE_TYPE) :: PE

!INTEGER :: nxy, nz, myzb, myze, myNxyBeg, myNxyEnd
INTEGER :: nproc, myrank, nComm

REAL, POINTER :: Buf(:)
!REAL, POINTER :: AxVar(:, :, :), RadVar(:, :, :)

INTEGER :: i, j, k
INTEGER :: ipartner, nsendrecvstep
INTEGER :: ixybeg, ixyend, izbeg, izend, istep, nsendpin, npln, npin
INTEGER :: ndatbase, nDatMax, ndat, ndat0
INTEGER :: ixy, iz, ixy1, ixy2
INTEGER :: COMM

!RadVar => CommAxRad%RadVar; AxVar => CommAxRad%AxVar;

nproc = PE%nCmfdProc; myrank = PE%myrank
Comm = PE%MPI_CMFD_COMM

nDatMax = nAxMax * nRadMax
Buf => RadAxCommBuf
DO istep = 1, nCommStep
  ipartner = MpiCommOrder(istep, myrank)
  IF(ipartner .LT. 0) THEN
    CYCLE
  ENDIF
  IF(myrank .LT. ipartner) THEN  !Send -> Recv
    ixybeg = PE%RadDomRange(1, ipartner); ixyend = PE%RadDomRange(2, ipartner)
    izbeg = PE%myzbf;izend = PE%myzef
    npin = ixyend - ixybeg + 1; npln = izend - izbeg + 1
    ndat = npin * npln; ndat0 = ndat * ndatBase
    CALL SendRecv(Buf(1:ndat0), ndat0, ipartner, .FALSE., COMM)
    CALL ReadAxRadCommBuf(RadVar(1:nDatBase, ixybeg:ixyend, izbeg:izend), Buf(1:ndat0), ndat0, nDatBase, npin, npln)
  ELSE  !Recv -> Send
    ixybeg = PE%RadDomRange(1, myrank); ixyend = PE%RadDomRange(2, myrank)
    izbeg = PE%SubPlnDomRange(1, ipartner); izend = PE%SubPlnDomRange(2, ipartner)
    npin = ixyend - ixybeg + 1; npln = izend - izbeg + 1
    ndat = npin * npln; ndat0 = ndat * ndatBase
    CALL WriteAxRadCommBuf(AxVar(:, ixybeg:ixyend, izbeg:izend), Buf(1:ndat0), ndat0, nDatBase, npin, npln)
    CALL SendRecv(Buf(1:ndat0), ndat0, ipartner, .TRUE., COMM)
  ENDIF
  IF(myrank .LT. ipartner) THEN 
    ixybeg = PE%RadDomRange(1, myrank); ixyend = PE%RadDomRange(2, myrank)
    izbeg = PE%SubPlnDomRange(1, ipartner); izend = PE%SubPlnDomRange(2, ipartner)
    npin = ixyend - ixybeg + 1; npln = izend - izbeg + 1
    ndat = npin * npln; ndat0 = ndat * ndatBase
    CALL WriteAxRadCommBuf(AxVar(:, ixybeg:ixyend, izbeg:izend), Buf(1:ndat0), ndat0, nDatBase, npin, npln)
    CALL SendRecv(Buf(1:ndat0), ndat0, ipartner, .TRUE., COMM)
  ELSE !Write
    ixybeg = PE%RadDomRange(1, ipartner); ixyend = PE%RadDomRange(2, ipartner)
    izbeg = PE%myzbf;izend = PE%myzef
    npin = ixyend - ixybeg + 1; npln = izend - izbeg + 1
    ndat = npin * npln; ndat0 = ndat * ndatBase
    CALL SendRecv(Buf(1:ndat0), ndat0, ipartner, .FALSE., COMM)
    CALL ReadAxRadCommBuf(RadVar(1:nDatBase, ixybeg:ixyend, izbeg:izend), Buf(1:ndat0), ndat0, nDatBase, npin, npln)
  ENDIF
ENDDO
CALL MPI_SYNC(comm)
NULLIFY(Buf)

END SUBROUTINE


SUBROUTINE WriteAxRadCommBuf(Dat,Buf, ndat0, ndat, n1, n2)
REAL :: Dat(ndat, n1, n2)
REAL :: Buf(ndat0)
INTEGER :: n1, n2, ndat, ndat0

INTEGER :: i, j, k, l 
l = 0
DO i = 1, n2
  DO j = 1, n1
    DO k = 1, ndat
      l = l + 1
      Buf(l) = Dat(k, j, i)
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE ReadAxRadCommBuf(Dat,Buf, ndat0, ndat, n1, n2)
REAL :: Dat(ndat, n1, n2)
REAL :: Buf(ndat0)
INTEGER :: n1, n2, ndat, ndat0
INTEGER :: i, j, k, l 
l = 0
DO i = 1, n2
  DO j = 1, n1
    DO k = 1, ndat
      l = l + 1
      Dat(k, j, i) = Buf(l)
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE
!
!SUBROUTINE GetPinXSData(SendRecvDat, ng, GroupInfo, PE)
!TYPE(AxNSendRecvDat_Type) :: SendRecvDat
!INTEGER :: ng
!TYPE(GroupInfo_Type) :: GroupInfo
!TYPE(PE_TYPE) :: PE
!
!REAL, POINTER :: SendPhi(:, :, :), SendLkg(:, :, :)
!REAL, POINTER :: RecvPhi(:, :, :), RecvLkg(:, :, :)
!TYPE(PinXS_Type), POINTER :: SendPinXS(:, :), RecvPinXS(:, :)
!
!INTEGER :: nxy, nz
!INTEGER :: myzb, myze, myNxyBeg, myNxyEnd
!INTEGER :: nmax
!
!LOGICAL, SAVE :: lFirst
!
!DATA lFirst /.TRUE./
!
!nxy = PE%nxy; nz = PE%nz
!myzb = PE%myzb; myze = PE%myze
!myNxyBeg = PE%myNxyBeg; myNxyEnd = PE%myNxyEnd
!
!nmax = MIN(100, nxy)
!
!IF(lFirst) THEN
!  CALL InitAxNodalComm(GroupInfo, ng, 4, nz, 100)
!  lFirst = .FALSE.
!ENDIF
!
!SendRecvDat%SendPhi => SendPhi; SendRecvDat%SendLkg => SendLkg
!SendRecvDat%RecvPhi => RecvPhi; SendRecvDat%RecvLkg => RecvLkg
!SendRecvDat%RecvPinXs => RecvPinXs; SendRecvDat%SendPinXs => SendPinXs
!
!
!NULLIFY(SendPhi, SendLkg)
!NULLIFY(RecvPhi, RecvLkg)
!NULLIFY(RecvPinXs, SendPinXs)
!END SUBROUTINE

SUBROUTINE FinalizePinXScomm()
IMPLICIT NONE
INTEGER :: iQueue
INTEGER :: i, j, k
INTEGER :: n1, n2, ng
REAL :: LKG(ngmax), PHI(ngmax)
TYPE(PinXS_Type), POINTER :: PinXS
INTEGER :: status(1000)

CALL MPI_WAITALL(nXsQueue, XsQueue(1:nXsQueue), status(1:nXsQueue), j)

ng = PinXsComm%ng
DO iQueue = 1, nXsQueue
  IF(XsQueueInfo(iQueue)%lSend) THEN
    NULLIFY(XsQueueInfo(iQueue)%SendPinXS)
    NULLIFY(XsQueueInfo(iQueue)%SendPhi)
    NULLIFY(XsQueueInfo(iQueue)%SendLkg)
  ELSE
    n1 = XsQueueInfo(iQueue)%n1
    n2 = XsQueueInfo(iQueue)%n2
    k = 0
    DO j = 1, n2
      DO i = 1, n1
        k = k + 1
        PHI(1:ng) = 0; LKG(1:ng) = 0
        PinXs => XsQueueInfo(iQueue)%RecvPinXS(i, j)
        CALL UnPackingPinXS(XsQueueInfo(iQueue)%RecvBuf(1:nTypeDat,k), PinXS, ng)
        XsQueueInfo(iQueue)%RecvPhi(1:ng, i, j) = Phi(1:ng)  
        XsQueueInfo(iQueue)%RecvLkg(1:ng, i, j) = Lkg(1:ng)  
      ENDDO
    ENDDO

    NULLIFY(XsQueueInfo(iQueue)%RecvPinXS);     NULLIFY(XsQueueInfo(iQueue)%RecvPhi)
    NULLIFY(XsQueueInfo(iQueue)%RecvLkg)
  ENDIF
ENDDO
nXsQueue = 0; nXsSend = 0; nXsRecv = 0
END SUBROUTINE


SUBROUTINE CommAxNodalInPut(Phi, Lkg, MpiPhi, MpiLkg, ng, PE)
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
INTEGER :: ng


REAL :: Phi(nxy, myzbf-1:myzef+1, ng), LKG(nxy, myzbf:myzef, ng)
REAL :: MpiPhi(myNxyBeg:myNxyEnd, nzfm, ng), MpiLkg(myNxyBeg:myNxyEnd, nzfm, ng)
INTEGER :: ixy, iz, K 

INTEGER :: ndat
INTEGER :: comm


comm = PE%MPI_CMFD_COMM
ndat = ng * 2

!Packing Data
DO iz = myzbf, myzef
  DO ixy = 1, nxy
    RadVar(1:ng, ixy, iz) = Phi(ixy, iz, 1:ng) 
    RadVar(ng + 1: 2 * ng, ixy, iz) = Lkg(ixy, iz, 1:ng) 
  ENDDO
ENDDO

CALL CommRadDom2AxDom(nDat, PE)

!UnPacking
DO iz = 1, nzfm
  DO ixy = myNxyBeg, myNxyEnd
    MpiPhi(ixy, iz, 1 : ng) = AxVar(1:ng, ixy, iz)
    MpiLkg(ixy, iz, 1 : ng) = AxVar(ng + 1 : 2 * ng, ixy, iz)
  ENDDO
ENDDO

!Copy Own Data
DO iz = myzbf, myzef
  DO ixy = myNxyBeg, myNxyEnd
    MpiPhi(ixy, iz, 1 : ng) = Phi(ixy, iz, 1 : ng)
    MpiLkg(ixy, iz, 1 : ng) = Lkg(ixy, iz, 1 : ng)
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE CommKineticParam(RadParam, AxParam, n1, n2, ng, SubPlaneMap, PE)
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
INTEGER :: n1, n2, ng


REAL, POINTER :: RadParam(:, :, :), AxParam(:, :, :)
INTEGER, POINTER :: SubPlaneMap(:)
INTEGER :: ixy, iz, iz0, k, i

INTEGER :: ndat, ndat0, nstep, ibeg, iend
INTEGER :: comm


comm = PE%MPI_CMFD_COMM
ndat = min(2*ng, n2-n1+1)
nstep = INT((n2-n1+1) / ndat) + 1
!Packing Data
DO i = 1, nstep
  ibeg = n1 + ndat * (i-1)
  iend = ibeg + ndat -1;   iend = min(iend, n2)
  ndat0 = iend - ibeg + 1
  DO iz = myzbf, myzef
    DO ixy = 1, nxy
      iz0 = SubPlaneMap(iz)
      RadVar(1:ndat0, ixy, iz) = RadParam(ibeg:iend, ixy, iz0) 
    ENDDO
  ENDDO

  CALL CommRadDom2AxDom(nDat0, PE)

  !UnPacking
  DO iz = 1, nzfm
    DO ixy = myNxyBeg, myNxyEnd!
      iz0 = SubPlaneMap(iz)
      AxParam(ibeg:iend, ixy, iz0) = AxVar(1:ndat0, ixy, iz)
    ENDDO
  ENDDO

  !Copy Own Data
  DO iz = myzbf, myzef
    DO ixy = myNxyBeg, myNxyEnd
      iz0 = SubPlaneMap(iz)
      AxParam(ibeg:iend, ixy, iz0) = RadParam(ibeg:iend, ixy, iz0) 
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE


SUBROUTINE CommAxNBeta(AxNBeta, PinXS, nprec, ng, SubPlaneMap,PE)
IMPLICIT NONE
REAL, POINTER :: AxNBeta(:, :, :)
INTEGER, POINTER :: SubPlaneMap(:)
TYPE(PinXS_TYPE), POINTER :: PinXS(:, :)
TYPE(PE_TYPE) :: PE
INTEGER :: nPrec, ng
INTEGER :: nstep, ndat, ndat0
INTEGER :: ixy, iz, iz0, i, ibeg, iend

 
ndat = min(2*ng, nprec)
nstep = INT((nprec) / ndat) + 1
!Packing Data
DO i = 1, nstep
  ibeg = 1 + ndat * (i-1)
  iend = ibeg + ndat -1;   iend = min(iend, nprec)
  ndat0 = iend - ibeg + 1
  DO iz = myzbf, myzef
    iz0 = SubPlaneMap(iz)
    DO ixy = 1, nxy
      RadVar(1:ndat0, ixy, iz) = PinXS(ixy, iz)%beta(ibeg:iend)
    ENDDO
  ENDDO
  
  CALL CommRadDom2AxDom(ndat0, PE)
  DO iz = 1, nzfm
    DO ixy = myNxyBeg, myNxyEnd!
      iz0 = SubPlaneMap(iz)
      AxNBeta(ibeg:iend, ixy, iz0) = AxVar(1:ndat0, ixy, iz)
    ENDDO
  ENDDO

  !Copy Own Data
  DO iz = myzbf, myzef
    DO ixy = myNxyBeg, myNxyEnd
      iz0 = SubPlaneMap(iz)
      AxNBeta(ibeg:iend, ixy, iz0) = PinXS(ixy, iz)%beta(ibeg:iend)
    ENDDO
  ENDDO  

ENDDO

END SUBROUTINE


SUBROUTINE CommTranSrc2nd(TranSrc, MpiTranSrc, ng, PE)
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
INTEGER :: ng

REAL, POINTER :: TranSrc(:, :, :)
REAL, POINTER :: MpiTranSrc(:, :, :)
INTEGER :: ixy, iz, K 

INTEGER :: ndat
INTEGER :: comm


comm = PE%MPI_CMFD_COMM
ndat = ng

!Packing Data
DO iz = myzbf, myzef
  DO ixy = 1, nxy
    RadVar(1:ng, ixy, iz) = TranSrc(ixy, iz, 1:ng) 
  ENDDO
ENDDO

CALL CommRadDom2AxDom(nDat, PE)

!UnPacking
DO iz = 1, nzfm
  DO ixy = myNxyBeg, myNxyEnd
    MpiTranSrc(ixy, iz, 1 : ng) = AxVar(1:ng, ixy, iz)
  ENDDO
ENDDO

!Copy Own Data
DO iz = myzbf, myzef
  DO ixy = myNxyBeg, myNxyEnd
    MpiTranSrc(ixy, iz, 1 : ng) = TranSrc(ixy, iz, 1 : ng)
  ENDDO
ENDDO

END SUBROUTINE


SUBROUTINE CommAxNodalOutPut(CmInfo, MpiFlx, ng, PE)
USE TYPEDEF,   ONLY : CmInfo_Type
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
INTEGER :: ng

TYPE(AxFlx_Type) :: MpiFlx(1:nzfm, myNxyBeg : myNxyEnd)
TYPE(CmInfo_Type) :: CmInfo
INTEGER :: ixy, iz, K 
INTEGER :: ib, ie
INTEGER :: ndat
INTEGER :: comm


comm = PE%MPI_CMFD_COMM
ndat = ng * 4

!Packing Data
DO ixy = myNxyBeg, myNxyEnd
  DO iz = 1, nzfm
    ib = 1; ie = ng
    AxVar(ib:ie, ixy, iz) = MpiFlx(iz, ixy)%Dtil(1, 1:ng)
    ib = ng + 1; ie = 2 * ng
    AxVar(ib:ie, ixy, iz) = MpiFlx(iz, ixy)%Dtil(2, 1:ng)
    ib = 2 * ng + 1; ie = 3 * ng
    AxVar(ib:ie, ixy, iz) = MpiFlx(iz, ixy)%Dhat(1, 1:ng)
    ib = 3 * ng + 1; ie = 4 * ng
    AxVar(ib:ie, ixy, iz) = MpiFlx(iz, ixy)%Dhat(2, 1:ng)
    MpiFlx(iz, ixy)%PDhat(1, 1:ng) = MpiFlx(iz, ixy)%Dhat(1, 1:ng)
    MpiFlx(iz, ixy)%PDhat(2, 1:ng) = MpiFlx(iz, ixy)%Dhat(2, 1:ng)
  ENDDO
ENDDO

CALL CommAxDom2RadDom(nDat, PE)

!UnPacking
DO ixy = 1, nxy
  DO iz = myzbf, myzef
    ib = 1; ie = ng
    CmInfo%AxDtil(1, ixy, iz, 1:ng) = RadVar(ib:ie, ixy, iz)
    !AxFlx(iz, ixy)%Dtil(1, 1:ng) = RadVar(ib:ie, ixy, iz)
    ib = ng + 1; ie = 2 * ng
    CmInfo%AxDtil(2, ixy, iz, 1:ng) = RadVar(ib:ie, ixy, iz)
    !AxFlx(iz, ixy)%Dtil(2, 1:ng) = RadVar(ib:ie, ixy, iz)
    ib = 2 * ng + 1; ie = 3 * ng
    CmInfo%AxDhat(1, ixy, iz, 1:ng) = RadVar(ib:ie, ixy, iz)
    !AxFlx(iz, ixy)%Dhat(1, 1:ng) = RadVar(ib:ie, ixy, iz)
    ib = 3 * ng + 1; ie = 4 * ng
    CmInfo%AxDhat(2, ixy, iz, 1:ng) = RadVar(ib:ie, ixy, iz)
    !AxFlx(iz, ixy)%Dhat(2, 1:ng) = RadVar(ib:ie, ixy, iz)
  ENDDO
ENDDO

!Copy Own Data
DO iz = myzbf, myzef
  DO ixy = myNxyBeg, myNxyEnd
    CmInfo%AxDtil(1, ixy, iz, 1:ng) = MpiFlx(iz, ixy)%Dtil(1, 1:ng)
    CmInfo%AxDtil(2, ixy, iz, 1:ng) = MpiFlx(iz, ixy)%Dtil(2, 1:ng)
    CmInfo%AxDhat(1, ixy, iz, 1:ng) = MpiFlx(iz, ixy)%Dhat(1, 1:ng)
    CmInfo%AxDhat(2, ixy, iz, 1:ng) = MpiFlx(iz, ixy)%Dhat(2, 1:ng)
    
!    AxFlx(iz, ixy)%Dtil(1, 1:ng) = MpiFlx(iz, ixy)%Dtil(1, 1:ng)
!    AxFlx(iz, ixy)%Dtil(2, 1:ng) = MpiFlx(iz, ixy)%Dtil(2, 1:ng)
!    AxFlx(iz, ixy)%Dhat(1, 1:ng) = MpiFlx(iz, ixy)%Dhat(1, 1:ng)
!    AxFlx(iz, ixy)%Dhat(2, 1:ng) = MpiFlx(iz, ixy)%Dhat(2, 1:ng)
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE CommAxNodalOutPutNew(CmInfo, MpiFlx, ng, PE, lUR)
USE TYPEDEF,   ONLY : CmInfo_Type
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
INTEGER :: ng

TYPE(AxFlx_Type) :: MpiFlx(1:nzfm, myNxyBeg : myNxyEnd)
TYPE(CmInfo_Type) :: CmInfo
INTEGER :: ixy, iz, K 
INTEGER :: ib, ie, ig
INTEGER :: ndat
INTEGER :: comm
LOGICAL :: lUR


comm = PE%MPI_CMFD_COMM
ndat = ng * 4

!Under-Relaxation
IF(lUR)THEN
!IF(PE%CMFDMaster) WRITE(*,*) 'AxDhat UR - inside Comm'
!DO ixy = myNxyBeg, myNxyEnd
!  DO iz = 1, nzfm
!      !WRITE(199,*) iz,ixy
!      !DO ig=1, ng
!      !WRITE(199,*) ig, MpiFlx(iz, ixy)%Dhat(1:2, ig)
!      !ENDDO
!    MpiFlx(iz, ixy)%Dhat(1, 1:ng)=MpiFlx(iz, ixy)%PDhat(1, 1:ng)*0.99+MpiFlx(iz, ixy)%Dhat(1, 1:ng)*0.01
!    MpiFlx(iz, ixy)%Dhat(2, 1:ng)=MpiFlx(iz, ixy)%PDhat(2, 1:ng)*0.99+MpiFlx(iz, ixy)%Dhat(2, 1:ng)*0.01
!  ENDDO
!ENDDO
ENDIF
!Packing Data
DO ixy = myNxyBeg, myNxyEnd
  DO iz = 1, nzfm
    ib = 1; ie = ng
    AxVar(ib:ie, ixy, iz) = MpiFlx(iz, ixy)%Dtil(1, 1:ng)
    ib = ng + 1; ie = 2 * ng
    AxVar(ib:ie, ixy, iz) = MpiFlx(iz, ixy)%Dtil(2, 1:ng)
    ib = 2 * ng + 1; ie = 3 * ng
    AxVar(ib:ie, ixy, iz) = MpiFlx(iz, ixy)%Dhat(1, 1:ng)
    ib = 3 * ng + 1; ie = 4 * ng
    AxVar(ib:ie, ixy, iz) = MpiFlx(iz, ixy)%Dhat(2, 1:ng)
    MpiFlx(iz, ixy)%PDhat(1, 1:ng) = MpiFlx(iz, ixy)%Dhat(1, 1:ng)
    MpiFlx(iz, ixy)%PDhat(2, 1:ng) = MpiFlx(iz, ixy)%Dhat(2, 1:ng)
  ENDDO
ENDDO

CALL CommAxDom2RadDom(nDat, PE)

!UnPacking
DO ixy = 1, nxy
  DO iz = myzbf, myzef
    ib = 1; ie = ng
    CmInfo%AxDtil(1, ixy, iz, 1:ng) = RadVar(ib:ie, ixy, iz)
    !AxFlx(iz, ixy)%Dtil(1, 1:ng) = RadVar(ib:ie, ixy, iz)
    ib = ng + 1; ie = 2 * ng
    CmInfo%AxDtil(2, ixy, iz, 1:ng) = RadVar(ib:ie, ixy, iz)
    !AxFlx(iz, ixy)%Dtil(2, 1:ng) = RadVar(ib:ie, ixy, iz)
    ib = 2 * ng + 1; ie = 3 * ng
    CmInfo%AxDhat(1, ixy, iz, 1:ng) = RadVar(ib:ie, ixy, iz)
    !AxFlx(iz, ixy)%Dhat(1, 1:ng) = RadVar(ib:ie, ixy, iz)
    ib = 3 * ng + 1; ie = 4 * ng
    CmInfo%AxDhat(2, ixy, iz, 1:ng) = RadVar(ib:ie, ixy, iz)
    !AxFlx(iz, ixy)%Dhat(2, 1:ng) = RadVar(ib:ie, ixy, iz)
  ENDDO
ENDDO

!Copy Own Data
DO iz = myzbf, myzef
  DO ixy = myNxyBeg, myNxyEnd
    CmInfo%AxDtil(1, ixy, iz, 1:ng) = MpiFlx(iz, ixy)%Dtil(1, 1:ng)
    CmInfo%AxDtil(2, ixy, iz, 1:ng) = MpiFlx(iz, ixy)%Dtil(2, 1:ng)
    CmInfo%AxDhat(1, ixy, iz, 1:ng) = MpiFlx(iz, ixy)%Dhat(1, 1:ng)
    CmInfo%AxDhat(2, ixy, iz, 1:ng) = MpiFlx(iz, ixy)%Dhat(2, 1:ng)
    
!    AxFlx(iz, ixy)%Dtil(1, 1:ng) = MpiFlx(iz, ixy)%Dtil(1, 1:ng)
!    AxFlx(iz, ixy)%Dtil(2, 1:ng) = MpiFlx(iz, ixy)%Dtil(2, 1:ng)
!    AxFlx(iz, ixy)%Dhat(1, 1:ng) = MpiFlx(iz, ixy)%Dhat(1, 1:ng)
!    AxFlx(iz, ixy)%Dhat(2, 1:ng) = MpiFlx(iz, ixy)%Dhat(2, 1:ng)
  ENDDO
ENDDO

END SUBROUTINE

END MODULE

