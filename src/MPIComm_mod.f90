MODULE MPIComm_Mod
!USE MPI
USE TIMER,            ONLY : nTracer_dclock, TimeChk
USE UtilFunction,     ONLY : ArrayLoc
IMPLICIT NONE
Include 'mpif.h'
INTEGER :: MPI_RT_COMM, MPI_CMFD_COMM
INTEGER, PARAMETER, PRIVATE :: nmax = 5000000
REAL,POINTER, SAVE :: SendBuf(:, :), RecvBuf(:, :), SendBuf2D(:, :, :), RecvBuf2D(:, :, :)
INTEGER :: ierr

INTERFACE BCAST
  MODULE PROCEDURE BCAST_REAL0D
  MODULE PROCEDURE BCAST_REAL1D
  MODULE PROCEDURE BCAST_REAL2D
  MODULE PROCEDURE BCAST_REAL3D
  MODULE PROCEDURE BCAST_INT0D
  MODULE PROCEDURE BCAST_INT1D
  MODULE PROCEDURE BCAST_INT2D
  MODULE PROCEDURE BCAST_INT3D
  MODULE PROCEDURE BCAST_LOGICAL
END INTERFACE

INTERFACE REDUCE
  MODULE PROCEDURE REDUCE_REAL0D
  MODULE PROCEDURE REDUCE_REAL1D
  MODULE PROCEDURE REDUCE_REAL2D
  MODULE PROCEDURE REDUCE_REAL3D
  MODULE PROCEDURE REDUCE_INT0D
  MODULE PROCEDURE REDUCE_INT1D
  MODULE PROCEDURE REDUCE_INT2D
  MODULE PROCEDURE REDUCE_INT3D
END INTERFACE

INTERFACE SENDRECV
  MODULE PROCEDURE SENDRECV_INT3D
  MODULE PROCEDURE SENDRECV_INT2D
  MODULE PROCEDURE SENDRECV_INT1D
  MODULE PROCEDURE SENDRECV_INT0D
  MODULE PROCEDURE SENDRECV_REAL4D
  MODULE PROCEDURE SENDRECV_REAL3D
  MODULE PROCEDURE SENDRECV_REAL2D
  MODULE PROCEDURE SENDRECV_REAL1D
  MODULE PROCEDURE SENDRECV_REAL0D
END INTERFACE

INTERFACE GetNeighDat
  MODULE PROCEDURE GetNeighDat_REAL4D
  MODULE PROCEDURE GetNeighDat_REAL3D
  MODULE PROCEDURE GetNeighDat_REAL2D  
  MODULE PROCEDURE GetNeighDat_REAL1D
  MODULE PROCEDURE GetNeighDat_INT4D
  MODULE PROCEDURE GetNeighDat_INT3D
  MODULE PROCEDURE GetNeighDat_INT2D
  MODULE PROCEDURE GetNeighDat_INT1D
END INTERFACE

INTERFACE GetNeighDatFast
  MODULE PROCEDURE GetNeighDatFast_REAL3D
  MODULE PROCEDURE GetNeighDatFast_REAL2D
END INTERFACE

INTERFACE NonBlockRecv
  MODULE PROCEDURE NonBlockRECV_REAL1D
  MODULE PROCEDURE NonBlockRECV_REAL2D
  MODULE PROCEDURE NonBlockRECV_REAL3D
  MODULE PROCEDURE NonBlockRECV_INT1D
  MODULE PROCEDURE NonBlockRECV_INT2D
  MODULE PROCEDURE NonBlockRECV_INT3D
END INTERFACE 

INTERFACE NonBlockSend
  MODULE PROCEDURE NonBlockSend_REAL1D
  MODULE PROCEDURE NonBlockSend_REAL2D
  MODULE PROCEDURE NonBlockSend_REAL3D
  MODULE PROCEDURE NonBlockSend_INT1D    
  MODULE PROCEDURE NonBlockSend_INT2D    
  MODULE PROCEDURE NonBlockSend_INT3D
END INTERFACE  

INTERFACE REDUCEnBCAST
  MODULE PROCEDURE REDUCEnBCAST_Real0D
  MODULE PROCEDURE REDUCEnBCAST_Int0D
END INTERFACE

CONTAINS
FUNCTION nCommStep(n)
INTEGER :: n, nCommStep
INTEGER :: i, j
i = n / nmax
j = mod(n, nmax)
IF(J .EQ. 0) THEN
 nCommStep = i
ELSE
 nCommStep = i + 1
ENDIF
END FUNCTION

SUBROUTINE REDUCE_REAL3D(Dat, DatBuf, n1, n2, n3, comm, lall)
IMPLICIT NONE
REAL :: DAT(n1, n2, n3)
REAL :: DatBuf(n1, n2, n3)
INTEGER :: n1, n2, n3, ndat, n, m

REAL:: Tbeg, Tend
INTEGER :: ib, ie, ie1, ie2, ie3, ib1, ib2, ib3
INTEGER :: comm
INTEGER :: i
INTEGER :: Dims(3), Locs(3)
LOGICAL :: lAll

tbeg = nTracer_dclock(.false., .false.)
ndat = n1 * n2 * n3
m = nCommStep(ndat)
Dims = (/ n1, n2, n3 /)
DO i = 1, m
  ib = nmax * (i - 1) + 1
  ie = nmax * i; ie = min(ie, ndat)
  n = ie - ib + 1
  CALL ArrayLoc(ib, 3, Dims, Locs)
  ib1 = Locs(1); ib2 = Locs(2); ib3 = Locs(3) 
!  CALL ArrayLoc(ie, 3, Dims, Locs)
!  ie1 = Locs(1); ie2 = Locs(2); ie3 = Locs(3)
  IF(.NOT. lAll) THEN
    CALL MPI_REDUCE(Dat(ib1, ib2, ib3), DatBuf(ib1, ib2, ib3), n, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, Comm, ierr)
  ELSE
    CALL MPI_ALLREDUCE(Dat(ib1, ib2, ib3), DatBuf(ib1, ib2, ib3), n, &
                      MPI_DOUBLE_PRECISION, MPI_SUM, Comm, ierr)
  ENDIF
ENDDO

tend = nTracer_dclock(.false., .false.)
TimeChk%CommTime = TimeChk%CommTime + tend - tbeg
END SUBROUTINE

SUBROUTINE REDUCE_INT3D(DAT, DatBuf, n1, n2, n3, comm, lAll)
IMPLICIT NONE
INTEGER :: DAT(n1, n2, n3)
INTEGER :: DatBuf(n1, n2, n3)
INTEGER :: n1, n2, n3, ndat, n, m

REAL:: Tbeg, Tend
INTEGER :: ib, ie, ie1, ie2, ie3, ib1, ib2, ib3
INTEGER :: comm
INTEGER :: i
INTEGER :: Dims(3), Locs(3)
LOGICAL :: lAll 

tbeg = nTracer_dclock(.false., .false.)
ndat = n1 * n2 * n3
m = nCommStep(ndat)
Dims = (/ n1, n2, n3 /)
DO i = 1, m
  ib = nmax * (i - 1) + 1
  ie = nmax * i; ie = min(ie, ndat)
  n = ie - ib + 1
  CALL ArrayLoc(ib, 3, Dims, Locs)
  ib1 = Locs(1); ib2 = Locs(2); ib3 = Locs(3) 
!  CALL ArrayLoc(ie, 3, Dims, Locs)
!  ie1 = Locs(1); ie2 = Locs(2); ie3 = Locs(3)
  IF(.NOT. lALl) THEN
    CALL MPI_REDUCE(Dat(ib1, ib2, ib3), DatBuf(ib1, ib2, ib3), n, &
                    MPI_INTEGER, MPI_SUM, 0, Comm, ierr)
  ELSE
    CALL MPI_ALLREDUCE(Dat(ib1, ib2, ib3), DatBuf(ib1, ib2, ib3), n, &
                    MPI_INTEGER, MPI_SUM, Comm, ierr)  
  ENDIF
ENDDO
tend = nTracer_dclock(.false., .false.)
TimeChk%CommTime = TimeChk%CommTime + tend - tbeg
END SUBROUTINE


SUBROUTINE REDUCE_REAL2D(DAT, DatBuf, n1, n2, comm, lAll)
IMPLICIT NONE
REAL :: Dat(n1, n2)
REAL :: DatBuf(n1, n2)
INTEGER :: n1, n2, ndat, n, m
INTEGER :: ib, ie, ie1, ie2, ib1, ib2
INTEGER :: comm
INTEGER :: i
INTEGER :: Dims(2), Locs(2)
LOGICAL :: lAll
REAL:: Tbeg, Tend

tbeg = nTracer_dclock(.false., .false.)

ndat = n1 * n2
m = nCommStep(ndat)
Dims = (/ n1, n2 /)

DO i = 1, m
  ib = nmax * (i - 1) + 1
  ie = nmax * i; ie = min(ie, ndat)
  n = ie - ib + 1
  CALL ArrayLoc(ib, 2, Dims, Locs)
  ib1 = Locs(1); ib2 = Locs(2)
!  CALL ArrayLoc(ie, 2, Dims, Locs)
!  ie1 = Locs(1); ie2 = Locs(2)
  IF(.NOT. lAll) THEN
    CALL MPI_REDUCE(Dat(ib1, ib2), DatBuf(ib1, ib2), n, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, 0, Comm, ierr)
  ELSE
    CALL MPI_ALLREDUCE(Dat(ib1, ib2), DatBuf(ib1, ib2), n, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, Comm, ierr)  
  ENDIF
ENDDO
tend = nTracer_dclock(.false., .false.)
TimeChk%CommTime = TimeChk%CommTime + tend - tbeg
END SUBROUTINE

SUBROUTINE REDUCE_INT2D(DAT, DatBuf, n1, n2, comm, lAll)
IMPLICIT NONE
INTEGER :: DAT(n1, n2)
INTEGER :: DatBuf(n1, n2)
INTEGER :: n1, n2, ndat, n, m
INTEGER :: ib, ie, ie1, ie2, ib1, ib2
INTEGER :: comm
INTEGER :: i
INTEGER :: Dims(2), Locs(2)
LOGICAL :: lAll
REAL:: Tbeg, Tend

tbeg = nTracer_dclock(.false., .false.)

ndat = n1 * n2
m = nCommStep(ndat)
Dims = (/ n1, n2 /)
DO i = 1, m
  ib = nmax * (i - 1) + 1
  ie = nmax * i; ie = min(ie, ndat)
  n = ie - ib + 1
  CALL ArrayLoc(ib, 2, Dims, Locs)
  ib1 = Locs(1); ib2 = Locs(2)
  IF(.NOT. lAll) THEN
    CALL MPI_REDUCE(Dat(ib1, ib2), DatBuf(ib1, ib2), n, &
                    MPI_INTEGER, MPI_SUM, 0, Comm, ierr)
  ELSE
    CALL MPI_ALLREDUCE(Dat(ib1, ib2), DatBuf(ib1, ib2), n, &
                    MPI_INTEGER, MPI_SUM, Comm, ierr)  
  ENDIF
ENDDO
tend = nTracer_dclock(.false., .false.)
TimeChk%CommTime = TimeChk%CommTime + tend - tbeg
END SUBROUTINE

SUBROUTINE REDUCE_REAL1D(Dat, DatBuf, ndat, comm, lAll)
REAL :: DAT(ndat)
REAL :: DATBuf(ndat)
INTEGER :: ndat, comm
INTEGER :: i, j, ib, ie
INTEGER :: n, m
LOGICAL :: lAll
REAL:: Tbeg, Tend

tbeg = nTracer_dclock(.false., .false.)

m = nCommStep(ndat)

DO i = 1, m
  ib = nmax * (i - 1) + 1
  ie = nmax * i; ie = min(ie, ndat)
  n = ie - ib + 1
  IF(.NOT. lAll) THEN
    CALL MPI_REDUCE(dat(ib), DatBuf(ib), n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
  ELSE
    CALL MPI_ALLREDUCE(dat(ib), DatBuf(ib), n, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
  ENDIF
  
ENDDO

tend = nTracer_dclock(.false., .false.)
TimeChk%CommTime = TimeChk%CommTime + tend - tbeg


END SUBROUTINE

SUBROUTINE REDUCE_INT1D(Dat, DatBuf, ndat, comm, lAll)
INTEGER :: DAT(ndat)
INTEGER :: DatBuf(ndat)
INTEGER :: ndat, comm
INTEGER :: i, j, ib, ie
INTEGER :: n, m
LOGICAL :: lAll
REAL:: Tbeg, Tend

tbeg = nTracer_dclock(.false., .false.)

m = nCommStep(ndat)

DO i = 1, m
  ib = nmax * (i - 1) + 1
  ie = nmax * i; ie = min(ie, ndat)
  n = ie - ib + 1
  IF(.NOT. lAll) THEN
    CALL MPI_REDUCE(dat(ib), DatBuf(ib), n, MPI_INTEGER, MPI_SUM, 0, comm, ierr)
  ELSE
    CALL MPI_ALLREDUCE(dat(ib), DatBuf(ib), n, MPI_INTEGER, MPI_SUM, comm, ierr)
  ENDIF
ENDDO

tend = nTracer_dclock(.false., .false.)
TimeChk%CommTime = TimeChk%CommTime + tend - tbeg


END SUBROUTINE

SUBROUTINE REDUCE_REAL0D(Dat, DatBuf, comm, lAll)
REAL :: DAT
REAL :: DATBuf
INTEGER :: ndat, comm
INTEGER :: i, j, ib, ie
INTEGER :: n
LOGICAL :: lAll
REAL:: Tbeg, Tend

tbeg = nTracer_dclock(.false., .false.)

n = 1
IF(.NOT. lAll) THEN
  CALL MPI_REDUCE(dat, DatBuf, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
ELSE
  CALL MPI_ALLREDUCE(dat, DatBuf, n, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
ENDIF

tend = nTracer_dclock(.false., .false.)
TimeChk%CommTime = TimeChk%CommTime + tend - tbeg


END SUBROUTINE

SUBROUTINE REDUCE_INT0D(Dat, DatBuf, comm, lAll)
INTEGER :: DAT
INTEGER :: DatBuf
INTEGER :: ndat, comm
INTEGER :: i, j, ib, ie
INTEGER :: n
LOGICAL :: lAll
REAL:: Tbeg, Tend

tbeg = nTracer_dclock(.false., .false.)
n = 1
IF(.NOT. lAll) THEN
  CALL MPI_REDUCE(dat, DatBuf, n, MPI_INTEGER, MPI_SUM, 0, comm, ierr)
ELSE
  CALL MPI_ALLREDUCE(dat, DatBuf, n, MPI_INTEGER, MPI_SUM, comm, ierr)
ENDIF


tend = nTracer_dclock(.false., .false.)
TimeChk%CommTime = TimeChk%CommTime + tend - tbeg

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BCAST_REAL3D(DAT, n1, n2, n3, comm, imaster0)
IMPLICIT NONE
REAL :: DAT(n1, n2, n3)
INTEGER :: n1, n2, n3, ndat, n, m
INTEGER :: ib, ie, ie1, ie2, ie3, ib1, ib2, ib3
INTEGER :: comm
INTEGER, OPTIONAL :: imaster0
INTEGER :: i, imaster
INTEGER :: Dims(3), Locs(3)
REAL:: Tbeg, Tend

tbeg = nTracer_dclock(.false., .false.)

imaster = 0
IF(PRESENT(imaster0)) imaster = imaster0

ndat = n1 * n2 * n3
m = nCommStep(ndat)
Dims = (/ n1, n2, n3 /)
DO i = 1, m
  ib = nmax * (i - 1) + 1
  ie = nmax * i; ie = min(ie, ndat)
  n = ie - ib + 1
  CALL ArrayLoc(ib, 3, Dims, Locs)
  ib1 = Locs(1); ib2 = Locs(2); ib3 = Locs(3) 
  CALL MPI_BCAST(Dat(ib1, ib2, ib3), n, MPI_DOUBLE_PRECISION, imaster, Comm, ierr)
ENDDO

tend = nTracer_dclock(.false., .false.)
TimeChk%CommTime = TimeChk%CommTime + tend - tbeg

END SUBROUTINE

SUBROUTINE BCAST_INT3D(DAT, n1, n2, n3, comm, imaster0)
IMPLICIT NONE
INTEGER :: DAT(n1, n2, n3)
INTEGER :: n1, n2, n3, ndat, n, m
INTEGER :: ib, ie, ie1, ie2, ie3, ib1, ib2, ib3
INTEGER :: comm
INTEGER, OPTIONAL :: imaster0

INTEGER :: i, imaster
INTEGER :: Dims(3), Locs(3)
REAL:: Tbeg, Tend

tbeg = nTracer_dclock(.false., .false.)

imaster = 0
IF(PRESENT(imaster0)) imaster = imaster0

ndat = n1 * n2 * n3
m = nCommStep(ndat)
Dims = (/ n1, n2, n3 /)
DO i = 1, m
  ib = nmax * (i - 1) + 1
  ie = nmax * i; ie = min(ie, ndat)
  n = ie - ib + 1
  CALL ArrayLoc(ib, 3, Dims, Locs)
  ib1 = Locs(1); ib2 = Locs(2); ib3 = Locs(3) 
  CALL MPI_BCAST(Dat(ib1, ib2, ib3), n, MPI_INTEGER, imaster, Comm, ierr)
ENDDO

tend = nTracer_dclock(.false., .false.)
TimeChk%CommTime = TimeChk%CommTime + tend - tbeg

END SUBROUTINE


SUBROUTINE BCAST_REAL2D(DAT, n1, n2, comm, imaster0)
IMPLICIT NONE
REAL :: DAT(n1, n2)
INTEGER :: n1, n2, ndat, n, m
INTEGER :: ib, ie, ie1, ie2, ib1, ib2
INTEGER :: comm
INTEGER, OPTIONAL :: imaster0

INTEGER :: i, imaster
INTEGER :: Dims(2), Locs(2)
REAL:: Tbeg, Tend

tbeg = nTracer_dclock(.false., .false.)
imaster = 0
IF(PRESENT(imaster0)) THEN
  imaster = imaster0
ENDIF

ndat = n1 * n2
m = nCommStep(ndat)
Dims = (/ n1, n2 /)
DO i = 1, m
  ib = nmax * (i - 1) + 1
  ie = nmax * i; ie = min(ie, ndat)
  n = ie - ib + 1
  CALL ArrayLoc(ib, 2, Dims, Locs)
  ib1 = Locs(1); ib2 = Locs(2)
  CALL MPI_BCAST(Dat(ib1, ib2), n, MPI_DOUBLE_PRECISION, imaster, Comm, ierr)
ENDDO

tend = nTracer_dclock(.false., .false.)
TimeChk%CommTime = TimeChk%CommTime + tend - tbeg


END SUBROUTINE

SUBROUTINE BCAST_INT2D(DAT, n1, n2, comm, imaster0)
IMPLICIT NONE
INTEGER :: DAT(n1, n2)
INTEGER :: n1, n2, ndat, n, m
INTEGER :: ib, ie, ie1, ie2, ib1, ib2
INTEGER :: comm
INTEGER, OPTIONAL :: imaster0
INTEGER :: i, imaster
INTEGER :: Dims(2), Locs(2)
REAL:: Tbeg, Tend

tbeg = nTracer_dclock(.false., .false.)

imaster = 0
IF(PRESENT(imaster0)) imaster = imaster0

ndat = n1 * n2
m = nCommStep(ndat)
Dims = (/ n1, n2 /)
DO i = 1, m
  ib = nmax * (i - 1) + 1
  ie = nmax * i; ie = min(ie, ndat)
  n = ie - ib + 1
  CALL ArrayLoc(ib, 2, Dims, Locs)
  ib1 = Locs(1); ib2 = Locs(2)
  CALL MPI_BCAST(Dat(ib1, ib2), n, MPI_INTEGER, imaster, Comm, ierr)
ENDDO

tend = nTracer_dclock(.false., .false.)
TimeChk%CommTime = TimeChk%CommTime + tend - tbeg

END SUBROUTINE

SUBROUTINE BCAST_REAL1D(Dat, ndat, comm, imaster0)
REAL :: DAT(ndat)
INTEGER :: ndat, comm
INTEGER, OPTIONAL :: imaster0
INTEGER :: i, imaster, j, ib, ie
INTEGER :: n, m
REAL:: Tbeg, Tend

tbeg = nTracer_dclock(.false., .false.)

imaster = 0
IF(PRESENT(imaster0)) imaster = imaster0

m = nCommStep(ndat)
DO i = 1, m
  ib = nmax * (i - 1) + 1
  ie = nmax * i; ie = min(ie, ndat)
  n = ie - ib + 1
  CALL MPI_BCAST(dat(ib), n, MPI_DOUBLE_PRECISION, imaster, comm, ierr)
ENDDO
tend = nTracer_dclock(.false., .false.)
TimeChk%CommTime = TimeChk%CommTime + tend - tbeg
END SUBROUTINE

SUBROUTINE BCAST_INT1D(Dat, ndat, comm, imaster0)
INTEGER :: DAT(ndat)
INTEGER :: ndat, comm
INTEGER, OPTIONAL :: imaster0
INTEGER :: i, j, imaster, ib, ie
INTEGER :: n, m
REAL:: Tbeg, Tend

tbeg = nTracer_dclock(.false., .false.)

imaster = 0
IF(PRESENT(imaster0)) imaster = imaster0

m = nCommStep(ndat)
DO i = 1, m
  ib = nmax * (i - 1) + 1
  ie = nmax * i; ie = min(ie, ndat)
  n = ie - ib + 1
  CALL MPI_BCAST(dat(ib), n, MPI_INTEGER, imaster, comm, ierr)
ENDDO

tend = nTracer_dclock(.false., .false.)
TimeChk%CommTime = TimeChk%CommTime + tend - tbeg

END SUBROUTINE

SUBROUTINE BCAST_REAL0D(Dat, comm, imaster0)
REAL :: DAT
INTEGER :: ndat, comm
INTEGER, OPTIONAL :: imaster0
INTEGER :: i, imaster, j, ib, ie
INTEGER :: n
REAL:: Tbeg, Tend

tbeg = nTracer_dclock(.false., .false.)

imaster = 0
IF(PRESENT(imaster0)) imaster = imaster0

ndat = 1
n=1
CALL MPI_BCAST(dat, n, MPI_DOUBLE_PRECISION, imaster, comm, ierr)

tend = nTracer_dclock(.false., .false.)
TimeChk%CommTime = TimeChk%CommTime + tend - tbeg

END SUBROUTINE

SUBROUTINE BCAST_INT0D(Dat, comm, imaster0)
INTEGER :: DAT
INTEGER :: ndat, comm
INTEGER, OPTIONAL :: imaster0
INTEGER :: i, imaster, j, ib, ie
INTEGER :: n
REAL:: Tbeg, Tend

tbeg = nTracer_dclock(.false., .false.)

imaster = 0
IF(PRESENT(imaster0)) imaster = imaster0

ndat = 1
n=1
CALL MPI_BCAST(dat, n, MPI_INTEGER, imaster, comm, ierr)

tend = nTracer_dclock(.false., .false.)
TimeChk%CommTime = TimeChk%CommTime + tend - tbeg
END SUBROUTINE

SUBROUTINE BCAST_LOGICAL(Dat, comm, imaster0)
LOGICAL :: DAT
INTEGER :: ndat, comm
INTEGER, OPTIONAL :: imaster0
INTEGER :: i, imaster, j, ib, ie
INTEGER :: n
REAL:: Tbeg, Tend

tbeg = nTracer_dclock(.false., .false.)

imaster = 0
IF(PRESENT(imaster0)) imaster = imaster0

ndat = 1
n=1
CALL MPI_BCAST(dat, n, MPI_LOGICAL, imaster, comm, ierr)

tend = nTracer_dclock(.false., .false.)
TimeChk%CommTime = TimeChk%CommTime + tend - tbeg

END SUBROUTINE



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SENDRECV_REAL4D(DAT, n1, n2, n3, n4, ipartner, lSend, comm)
IMPLICIT NONE
REAL :: DAT(n1, n2, n3, n4)
INTEGER :: n1, n2, n3, n4, ndat, n, m, ipartner
INTEGER :: comm
LOGICAL :: lSend

INTEGER :: ib, ie, ie1, ie2, ie3, ie4, ib1, ib2, ib3, ib4
INTEGER :: i, ierr, istatus(MPI_STATUS_SIZE)
INTEGER :: Dims(4), Locs(4)
REAL:: Tbeg, Tend

tbeg = nTracer_dclock(.false., .false.)

ndat = n1 * n2 * n3 * n4
m = nCommStep(ndat)
Dims = (/ n1, n2, n3, n4 /)
DO i = 1, m
  ib = nmax * (i - 1) + 1
  ie = nmax * i; ie = min(ie, ndat)
  n = ie - ib + 1
  CALL ArrayLoc(ib, 4, Dims, Locs)
  ib1 = Locs(1); ib2 = Locs(2); ib3 = Locs(3); ib4 = Locs(4)
  IF(lSend) THEN
    CALL MPI_SEND(Dat(ib1, ib2, ib3, ib4), n, MPI_DOUBLE_PRECISION, ipartner, 1, comm, ierr)
  ELSE
    CALL MPI_RECV(Dat(ib1, ib2, ib3, ib4), n, MPI_DOUBLE_PRECISION, ipartner, 1, comm, istatus, ierr)
  ENDIF  
ENDDO

tend = nTracer_dclock(.false., .false.)
!1TimeChk%CommTime = TimeChk%CommTime + tend - tbeg
END SUBROUTINE


SUBROUTINE SENDRECV_REAL3D(DAT, n1, n2, n3, ipartner, lSend, comm)
IMPLICIT NONE
REAL :: DAT(n1, n2, n3)
INTEGER :: n1, n2, n3, ndat, n, m, ipartner
INTEGER :: comm
LOGICAL :: lSend

INTEGER :: ib, ie, ie1, ie2, ie3, ib1, ib2, ib3
INTEGER :: i, ierr, istatus(MPI_STATUS_SIZE)
INTEGER :: Dims(3), Locs(3)
REAL:: Tbeg, Tend

tbeg = nTracer_dclock(.false., .false.)

ndat = n1 * n2 * n3
m = nCommStep(ndat)
Dims = (/ n1, n2, n3 /)
DO i = 1, m
  ib = nmax * (i - 1) + 1
  ie = nmax * i; ie = min(ie, ndat)
  n = ie - ib + 1
  CALL ArrayLoc(ib, 3, Dims, Locs)
  ib1 = Locs(1); ib2 = Locs(2); ib3 = Locs(3) 
  IF(lSend) THEN
    CALL MPI_SEND(Dat(ib1, ib2, ib3), n, MPI_DOUBLE_PRECISION, ipartner, 1, comm, ierr)
  ELSE
    CALL MPI_RECV(Dat(ib1, ib2, ib3), n, MPI_DOUBLE_PRECISION, ipartner, 1, comm, istatus, ierr)
  ENDIF  
ENDDO

tend = nTracer_dclock(.false., .false.)
TimeChk%CommTime = TimeChk%CommTime + tend - tbeg

END SUBROUTINE

SUBROUTINE SENDRECV_INT3D(DAT, n1, n2, n3, ipartner, lSend, comm)
IMPLICIT NONE
INTEGER :: DAT(n1, n2, n3)
INTEGER :: n1, n2, n3, ndat, n, m, ipartner
INTEGER :: comm
LOGICAL :: lSend

INTEGER :: ib, ie, ie1, ie2, ie3, ib1, ib2, ib3
INTEGER :: i, ierr, istatus(MPI_STATUS_SIZE)
INTEGER :: Dims(3), Locs(3)
REAL:: Tbeg, Tend

tbeg = nTracer_dclock(.false., .false.)

ndat = n1 * n2 * n3
m = nCommStep(ndat)
Dims = (/ n1, n2, n3 /)
DO i = 1, m
  ib = nmax * (i - 1) + 1
  ie = nmax * i; ie = min(ie, ndat)
  n = ie - ib + 1
  CALL ArrayLoc(ib, 3, Dims, Locs)
  ib1 = Locs(1); ib2 = Locs(2); ib3 = Locs(3) 
  IF(lSend) THEN
    CALL MPI_SEND(Dat(ib1, ib2, ib3), n, MPI_INTEGER, ipartner, 1, comm, ierr)
  ELSE
    CALL MPI_RECV(Dat(ib1, ib2, ib3), n, MPI_INTEGER, ipartner, 1, comm, istatus, ierr)
  ENDIF  
ENDDO

tend = nTracer_dclock(.false., .false.)
TimeChk%CommTime = TimeChk%CommTime + tend - tbeg

END SUBROUTINE

SUBROUTINE SENDRECV_REAL2D(DAT, n1, n2, ipartner, lSend, comm)
IMPLICIT NONE
REAL :: DAT(n1, n2)
INTEGER :: n1, n2, ndat, n, m, ipartner
LOGICAL :: lSend
INTEGER :: ib, ie, ie1, ie2, ib1, ib2
INTEGER :: comm
INTEGER :: i, ierr, istatus(MPI_STATUS_SIZE)
INTEGER :: Dims(2), Locs(2)
REAL:: Tbeg, Tend

tbeg = nTracer_dclock(.false., .false.)

ndat = n1 * n2
m = nCommStep(ndat)
Dims = (/ n1, n2 /)
DO i = 1, m
  ib = nmax * (i - 1) + 1
  ie = nmax * i; ie = min(ie, ndat)
  n = ie - ib + 1
  CALL ArrayLoc(ib, 2, Dims, Locs)
  ib1 = Locs(1); ib2 = Locs(2)
  IF(lSend) THEN
    CALL MPI_SEND(Dat(ib1, ib2), n, MPI_DOUBLE_PRECISION, ipartner, 1, comm, ierr)
  ELSE
    CALL MPI_RECV(Dat(ib1, ib2), n, MPI_DOUBLE_PRECISION, ipartner, 1, comm, istatus, ierr)
  ENDIF  
ENDDO
tend = nTracer_dclock(.false., .false.)
TimeChk%CommTime = TimeChk%CommTime + tend - tbeg
END SUBROUTINE

SUBROUTINE SENDRECV_INT2D(DAT, n1, n2, ipartner,lSend, comm)
IMPLICIT NONE
INTEGER :: DAT(n1, n2)
INTEGER :: n1, n2, ndat, n, m, ipartner
LOGICAL :: lSend
INTEGER :: ib, ie, ie1, ie2, ib1, ib2
INTEGER :: comm
INTEGER :: i, ierr, istatus(MPI_STATUS_SIZE)
INTEGER :: Dims(2), Locs(2)
REAL:: Tbeg, Tend

tbeg = nTracer_dclock(.false., .false.)

ndat = n1 * n2
m = nCommStep(ndat)
Dims = (/ n1, n2 /)
DO i = 1, m
  ib = nmax * (i - 1) + 1
  ie = nmax * i; ie = min(ie, ndat)
  n = ie - ib + 1
  CALL ArrayLoc(ib, 2, Dims, Locs)
  ib1 = Locs(1); ib2 = Locs(2)
  IF(lSend) THEN
    CALL MPI_SEND(Dat(ib1, ib2), n, MPI_INTEGER, ipartner, 1, comm, ierr)
  ELSE
    CALL MPI_RECV(Dat(ib1, ib2), n, MPI_INTEGER, ipartner, 1, comm, istatus, ierr)
  ENDIF  
ENDDO
tend = nTracer_dclock(.false., .false.)
TimeChk%CommTime = TimeChk%CommTime + tend - tbeg
END SUBROUTINE

SUBROUTINE SENDRECV_REAL1D(Dat, ndat, ipartner, lSend, comm)
REAL :: DAT(ndat)
INTEGER :: ndat, ipartner, comm
INTEGER :: i, j, ib, ie, ierr, istatus(MPI_STATUS_SIZE)
INTEGER :: n, m
LOGICAL :: lSend
REAL:: Tbeg, Tend


tbeg = nTracer_dclock(.false., .false.)
m = nCommStep(ndat)
DO i = 1, m
  ib = nmax * (i - 1) + 1
  ie = nmax * i; ie = min(ie, ndat)
  n = ie - ib + 1
  IF(lSend) THEN
    CALL MPI_SEND(Dat(ib), n, MPI_DOUBLE_PRECISION, ipartner, 1, comm, ierr)
  ELSE
    CALL MPI_RECV(Dat(ib), n, MPI_DOUBLE_PRECISION, ipartner, 1, comm, istatus, ierr)
  ENDIF
ENDDO
tend = nTracer_dclock(.false., .false.)
TimeChk%CommTime = TimeChk%CommTime + tend - tbeg
END SUBROUTINE

SUBROUTINE SENDRECV_INT1D(Dat, ndat, ipartner, lSend, comm)
INTEGER :: DAT(ndat)
INTEGER :: ndat, ipartner, comm
INTEGER :: i, j, ib, ie, ierr, istatus(MPI_STATUS_SIZE)
INTEGER :: n, m
LOGICAL :: lSend
REAL:: Tbeg, Tend

tbeg = nTracer_dclock(.false., .false.)

m = nCommStep(ndat)
DO i = 1, m
  ib = nmax * (i - 1) + 1
  ie = nmax * i; ie = min(ie, ndat)
  n = ie - ib + 1
  IF(lSend) THEN
    CALL MPI_SEND(Dat(ib), n, MPI_INTEGER, ipartner, 1, comm, ierr)
  ELSE
    CALL MPI_RECV(Dat(ib), n, MPI_INTEGER, ipartner, 1, comm, istatus, ierr)
  ENDIF
ENDDO
tend = nTracer_dclock(.false., .false.)
!TimeChk%CommTime = TimeChk%CommTime + tend - tbeg
END SUBROUTINE

SUBROUTINE SENDRECV_REAL0D(Dat, ipartner, lsend, comm)
REAL :: DAT
INTEGER :: ipartner, comm
LOGICAL :: lsend
INTEGER :: ierr, istatus(MPI_STATUS_SIZE)
REAL:: Tbeg, Tend

tbeg = nTracer_dclock(.false., .false.)

IF(lSend) THEN
  CALL MPI_SEND(Dat, 1, MPI_DOUBLE_PRECISION, ipartner, 1, comm, ierr)
ELSE
  CALL MPI_RECV(Dat, 1, MPI_DOUBLE_PRECISION, ipartner, 1, comm, istatus, ierr)
ENDIF

tend = nTracer_dclock(.false., .false.)
TimeChk%CommTime = TimeChk%CommTime + tend - tbeg

END SUBROUTINE

SUBROUTINE SENDRECV_INT0D(Dat, ipartner, lsend, comm)
INTEGER :: DAT
INTEGER :: ipartner, comm
LOGICAL :: lsend
INTEGER :: ierr, istatus(MPI_STATUS_SIZE)
REAL:: Tbeg, Tend

tbeg = nTracer_dclock(.false., .false.)

IF(lSend) THEN
  CALL MPI_SEND(Dat, 1, MPI_INTEGER, ipartner, 1, comm, ierr)
ELSE
  CALL MPI_RECV(Dat, 1, MPI_INTEGER, ipartner, 1, comm, istatus, ierr)
ENDIF

tend = nTracer_dclock(.false., .false.)
!TimeChk%CommTime = TimeChk%CommTime + tend - tbeg
!PRINT *, 'Send recv_int0', TimeChk%CommTime
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE NonBlockRECV_INT3D(DAT, n1, n2, n3, Dest, comm, request)
INTEGER :: DAT(n1, n2, n3)
INTEGER :: n1, n2, n3
INTEGER :: Dest, comm, request
INTEGER :: ierr, istatus
INTEGER :: n 

n = n1 * n2 * n3
CALL MPI_IRECV(Dat(1, 1, 1), n, MPI_INTEGER, Dest, 0, Comm, Request, ierr)

END SUBROUTINE

SUBROUTINE NonBlockRECV_REAL3D(DAT, n1, n2, n3, Dest, comm, request)
REAL :: DAT(n1, n2, n3)
INTEGER :: n1, n2, n3
INTEGER :: Dest, comm, request
INTEGER :: ierr, istatus
INTEGER :: n

n = n1 * n2 * n3
CALL MPI_IRECV(Dat(1, 1, 1), n, MPI_DOUBLE_PRECISION, Dest, 0, Comm, Request, ierr)

END SUBROUTINE

SUBROUTINE NonBlockRECV_INT2D(DAT, n1, n2, Dest, comm, request)
INTEGER :: DAT(n1, n2)
INTEGER :: n1, n2
INTEGER :: Dest, comm, request
INTEGER :: ierr, istatus
INTEGER :: n 

n = n1 * n2
CALL MPI_IRECV(Dat(1, 1), n, MPI_INTEGER, Dest, 0, Comm, Request, ierr)

END SUBROUTINE

SUBROUTINE NonBlockRECV_REAL2D(DAT, n1, n2, Dest, comm, request)
REAL :: DAT(n1, n2)
INTEGER :: n1, n2
INTEGER :: Dest, comm, request
INTEGER :: ierr, istatus
INTEGER :: n

n = n1 * n2
CALL MPI_IRECV(Dat(1, 1), n, MPI_DOUBLE_PRECISION, Dest, 0, Comm, Request, ierr)

END SUBROUTINE

SUBROUTINE NonBlockRECV_INT1D(DAT, n1, Dest, comm, request)
INTEGER :: DAT(n1)
INTEGER :: n1
INTEGER :: Dest, comm, request
INTEGER :: ierr, istatus

CALL MPI_IRECV(Dat(1), n1, MPI_INTEGER, Dest, 0, Comm, Request, ierr)

END SUBROUTINE

SUBROUTINE NonBlockRECV_REAL1D(DAT, n1, Dest, comm, request)
REAL :: DAT(n1)
INTEGER :: n1
INTEGER :: Dest, comm, request
INTEGER :: ierr, istatus

CALL MPI_IRECV(Dat(1), n1, MPI_DOUBLE_PRECISION, Dest, 0, Comm, Request, ierr)

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE NonBlockSend_INT3D(DAT, n1, n2, n3, Dest, comm, request)
INTEGER :: DAT(n1, n2, n3)
INTEGER :: n1, n2, n3
INTEGER :: Dest, comm, request
INTEGER :: ierr, istatus
INTEGER :: n 

n = n1 * n2 * n3
CALL MPI_ISEND(Dat(1, 1, 1), n, MPI_INTEGER, Dest, 0, Comm, Request, ierr)

END SUBROUTINE

SUBROUTINE NonBlockSend_REAL3D(DAT, n1, n2, n3, Dest, comm, request)
REAL :: DAT(n1, n2, n3)
INTEGER :: n1, n2, n3
INTEGER :: Dest, comm, request
INTEGER :: ierr, istatus
INTEGER :: n

n = n1 * n2 * n3
CALL MPI_ISEND(Dat(1, 1, 1), n, MPI_DOUBLE_PRECISION, Dest, 0, Comm, Request, ierr)

END SUBROUTINE

SUBROUTINE NonBlockSend_INT2D(DAT, n1, n2, Dest, comm, request)
INTEGER :: DAT(n1, n2)
INTEGER :: n1, n2
INTEGER :: Dest, comm, request
INTEGER :: ierr, istatus
INTEGER :: n 

n = n1 * n2
CALL MPI_ISEND(Dat(1, 1), n1, MPI_INTEGER, Dest, 0, Comm, Request, ierr)

END SUBROUTINE

SUBROUTINE NonBlockSend_REAL2D(DAT, n1, n2, Dest, comm, request)
REAL :: DAT(n1, n2)
INTEGER :: n1, n2
INTEGER :: Dest, comm, request
INTEGER :: ierr, istatus
INTEGER :: n

n = n1 * n2
CALL MPI_ISEND(Dat(1, 1), n, MPI_DOUBLE_PRECISION, Dest, 0, Comm, Request, ierr)

END SUBROUTINE

SUBROUTINE NonBlockSend_INT1D(DAT, n1, Dest, comm, request)
INTEGER :: DAT(n1)
INTEGER :: n1
INTEGER :: Dest, comm, request
INTEGER :: ierr, istatus

CALL MPI_ISEND(Dat(1), n1, MPI_INTEGER, Dest, 0, Comm, Request, ierr)

END SUBROUTINE

SUBROUTINE NonBlockSend_REAL1D(DAT, n1, Dest, comm, request)
REAL :: DAT(n1)
INTEGER :: n1
INTEGER :: Dest, comm, request
INTEGER :: ierr, istatus

CALL MPI_ISEND(Dat(1), n1, MPI_DOUBLE_PRECISION, Dest, 0, Comm, Request, ierr)

END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ChkNonBlockComm(nRequest, RequestList)
IMPLICIT NONE
INTEGER :: nRequest
INTEGER :: RequestList(nRequest)
INTEGER :: status(nRequest), ierr
!CALL MPI_Waitall(nRequest, RequestList(1:nRequest), status(1:nRequest), ierr)
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GetNeighDatFast_REAL3D(Dat, n1, n2, myzb, myze, myrank, nProc, Comm, imod)
USE BasicOperation, only :  CP_VA, CP_CA
IMPLICIT NONE
REAL, TARGET :: Dat(n1, n2, myzb-1:myze+1)
INTEGER :: n1, n2, myzb, myze, myrank, nproc, Comm, imod
INTEGER, SAVE :: Queue(4), nQueue
Integer :: ierr, status(MPI_STATUS_SIZE, 4)

INTEGER :: i, j, istep, igrp, ipartner
LOGICAL :: lSend
LOGICAL, SAVE :: lfirst
data lfirst /.TRUE./

IF(lFirst) THEN
  ALLOCATE(SEndBuf2D(n1, n2, 2))
  ALLOCATE(RecvBuf2D(n1, n2, 2))
  lfirst = .FALSE.
ENDIF
IF(imod .eq. 0) THEN
  nQueue = 0
  CALL CP_CA(RecvBuf2D(1:n1, 1:n2, 1:2), 0._8, n1, n2, 2)
  CALL CP_VA(SendBuf2D(1:n1, 1:n2, 1), Dat(1:n1, 1:n2, myzb), n1, n2)
  CALL CP_VA(SendBuf2D(1:n1, 1:n2, 2), Dat(1:n1, 1:n2, myze), n1, n2)
  DO istep = 0, 1
    i = MOD(myrank+istep, 2)
    IF(i .EQ. 0) THEN
      ipartner = myrank + 1; lSend = .TRUE.
      IF( ipartner .LT. 0 .OR. ipartner .GT. nproc-1) CYCLE    
      nQueue = nQueue + 1
      CALL NonBlockSend(SendBuf2D(:, :, 2), n1, n2, ipartner, Comm, Queue(nQueue))
      nQueue = nQueue + 1
      CALL NonBlockRecv(RecvBuf2D(:, :, 2), n1, n2, ipartner, Comm, Queue(nQueue))
    ENDIF
    IF(i .EQ. 1) THEN
      ipartner = myrank - 1; lSend = .FALSE.
      IF( ipartner .LT. 0 .OR. ipartner .GT. nproc-1) CYCLE    
      nQueue = nQueue + 1
      CALL NonBlockRecv(RecvBuf2D(:, :,1), n1, n2, ipartner, Comm, Queue(nQueue))
      nQueue = nQueue + 1
      CALL NonBlockSend(SendBuf2D(:, :, 1), n1, n2, ipartner, Comm, Queue(nQueue))
    ENDIF
  ENDDO
ELSE
  CALL MPI_WAITALL(nQueue, Queue(1:nQueue), Status, ierr)
  CALL CP_VA(Dat(1:n1, 1:n2, myzb-1), RecvBuf2D(1:n1, 1:n2, 1), n1, n2)
  CALL CP_VA(Dat(1:n1, 1:n2, myze+1), RecvBuf2D(1:n1, 1:n2, 2), n1, n2)
ENDIF

END SUBROUTINE


SUBROUTINE GetNeighDatFast_REAL2D(Dat, n1, myzb, myze, myrank, nProc, Comm, imod)
USE BasicOperation, only :  CP_VA, CP_CA
IMPLICIT NONE
REAL, TARGET :: Dat(n1, myzb-1:myze+1)
INTEGER :: n1, myzb, myze, myrank, nproc, Comm, imod
INTEGER, SAVE :: Queue(4), nQueue
Integer :: ierr, status(MPI_STATUS_SIZE, 4)

INTEGER :: i, j, istep, igrp, ipartner
LOGICAL :: lSend
LOGICAL, SAVE :: lfirst
data lfirst /.TRUE./

IF(lFirst) THEN
  ALLOCATE(SEndBuf(n1, 2))
  ALLOCATE(RecvBuf(n1, 2))
  lfirst = .FALSE.
ENDIF

if(imod .eq. 0) THEN
  nQueue = 0
  CALL CP_CA(RecvBuf(1:n1, 1:2), 0._8, n1, 2)
  CALL CP_VA(SendBuf(1:n1, 1), Dat(1:n1, myzb), n1)
  CALL CP_VA(SendBuf(1:n1, 2), Dat(1:n1, myze), n1)
  DO istep = 0, 1
    i = MOD(myrank+istep, 2)
    IF(i .EQ. 0) THEN
      ipartner = myrank + 1; lSend = .TRUE.
      IF( ipartner .LT. 0 .OR. ipartner .GT. nproc-1) CYCLE    
      nQueue = nQueue + 1
      CALL NonBlockSend(SendBuf(:, 2), n1,  ipartner, Comm, Queue(nQueue))
      nQueue = nQueue + 1
      CALL NonBlockRecv(RecvBuf(:, 2), n1, ipartner, Comm, Queue(nQueue))
    ENDIF
    IF(i .EQ. 1) THEN
      ipartner = myrank - 1; lSend = .FALSE.
      IF( ipartner .LT. 0 .OR. ipartner .GT. nproc-1) CYCLE    
      nQueue = nQueue + 1
      CALL NonBlockRecv(RecvBuf(:, 1), n1, ipartner, Comm, Queue(nQueue))
      nQueue = nQueue + 1
      CALL NonBlockSend(SendBuf(:, 1), n1,  ipartner, Comm, Queue(nQueue))
    ENDIF
  ENDDO
ELSE
  CALL MPI_WAITALL(nQueue, Queue(1:nQueue), Status, ierr)
  CALL CP_VA(Dat(1:n1, myzb-1), RecvBuf(1:n1, 1), n1)
  CALL CP_VA(Dat(1:n1, myze+1), RecvBuf(1:n1, 2), n1)
ENDIF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GetNeighDat_REAL4D(Dat, n1, n2, n3, myzb, myze, myrank, nProc, Comm)
REAL, TARGET :: Dat(n1, n2, n3, myzb-1:myze+1)
REAL, POINTER :: Dat0(:, :, :), Dat1(:, :, :)
INTEGER :: n1, n2, n3, myzb, myze, myrank, nproc, Comm

INTEGER :: i, j, istep, igrp, ipartner
LOGICAL :: lSend

DO istep = 0, 1 
  i = MOD(myrank+istep, 2)
  IF(i .EQ. 0) THEN
    ipartner = myrank + 1; lSend = .TRUE.
    Dat0 => Dat(:, :, :, myze); Dat1 => Dat(:, :, :, myze + 1)
  ENDIF
  IF(i .EQ. 1) THEN
    ipartner = myrank - 1; lSend = .FALSE.
    Dat1 => Dat(:, :, :, myzb); Dat0 => Dat(:, :, :, myzb -1)
  ENDIF
  IF( ipartner .LT. 0 .OR. ipartner .GT. nproc-1) CYCLE
  CALL SENDRECV(Dat0, n1, n2, n3, ipartner, lSend, Comm)
  CALL SENDRECV(Dat1, n1, n2, n3, ipartner, .NOT. lSend, Comm)
ENDDO

END SUBROUTINE

SUBROUTINE GetNeighDat_INT4D(Dat, n1, n2, n3, myzb, myze, myrank, nProc, Comm)
INTEGER, TARGET :: Dat(n1, n2, n3, myzb-1:myze+1)
INTEGER, POINTER :: Dat0(:, :, :), Dat1(:, :, :)
INTEGER :: n1, n2, n3, myzb, myze, myrank, nproc, Comm

INTEGER :: i, j, istep, igrp, ipartner
LOGICAL :: lSend

DO istep = 0, 1 
  i = MOD(myrank+istep, 2)
  IF(i .EQ. 0) THEN
    ipartner = myrank + 1; lSend = .TRUE.
    Dat0 => Dat(:, :, :, myze); Dat1 => Dat(:, :, :, myze + 1)
  ENDIF
  IF(i .EQ. 1) THEN
    ipartner = myrank - 1; lSend = .FALSE.
    Dat1 => Dat(:, :, :, myzb); Dat0 => Dat(:, :, :, myzb -1)
  ENDIF
  IF( ipartner .LT. 0 .OR. ipartner .GT. nproc-1) CYCLE
  CALL SENDRECV(Dat0, n1, n2, n3, ipartner, lSend, Comm)
  CALL SENDRECV(Dat1, n1, n2, n3, ipartner, .NOT. lSend, Comm)
ENDDO

END SUBROUTINE


SUBROUTINE GetNeighDat_REAL3D(Dat, n1, n2, myzb, myze, myrank, nProc, Comm)
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
REAL, TARGET :: Dat(n1, n2, myzb-1:myze+1)
REAL, POINTER :: Dat0(:, :), Dat1(:, :)
INTEGER :: n1, n2, myzb, myze, myrank, nproc, Comm

INTEGER :: i, j, istep, igrp, ipartner
LOGICAL :: lSend

DO istep = 0, 1 
  i = MOD(myrank+istep, 2)
  IF(i .EQ. 0) THEN
    ipartner = myrank + 1; lSend = .TRUE.
    IF( ipartner .LT. 0 .OR. ipartner .GT. nproc-1) CYCLE    
    CALL CP_CA(Dat(:, :, myze+1), 0._8, n1, n2)
    CALL SENDRECV(Dat(:, :, myze), n1, n2, ipartner, .true., Comm)
    CALL SENDRECV(Dat(:, :, myze+1), n1, n2, ipartner,.false., Comm)    
    !Dat0 => Dat(:, :, myze); Dat1 => Dat(:, :, myze + 1)
  ENDIF
  IF(i .EQ. 1) THEN
    ipartner = myrank - 1; lSend = .FALSE.
    IF( ipartner .LT. 0 .OR. ipartner .GT. nproc-1) CYCLE    
    CALL CP_CA(Dat(:, :, myzb-1), 0._8, n1, n2)
    CALL SENDRECV(Dat(:, :, myzb-1), n1, n2, ipartner,.false., Comm)   
    CALL SENDRECV(Dat(:, :, myzb), n1, n2, ipartner, .true., Comm)    
   !Dat1 => Dat(:, :, myzb); Dat0 => Dat(:, :, myzb -1)
  ENDIF
ENDDO

!  i = MOD(myrank+istep, 2)
!  IF(i .EQ. 0) THEN
!    ipartner = myrank + 1; lSend = .TRUE.
!    CALL CP_CA(Dat(:, myze+1), 0._8, n)
!    IF( ipartner .LT. 0 .OR. ipartner .GT. nproc-1) CYCLE    
!    !Dat0 => Dat(:, myze); Dat1 => Dat(:, myze + 1)
!    CALL SENDRECV(Dat(:, myze), n, ipartner, .true., Comm)
!    CALL SENDRECV(Dat(:, myze+1), n, ipartner,.false., Comm)
!  ELSE
!    ipartner = myrank - 1; lSend = .FALSE.
!    CALL CP_CA(Dat(:, myzb-1), 0._8, n)
!    IF( ipartner .LT. 0 .OR. ipartner .GT. nproc-1) CYCLE    
!    !Dat1 => Dat(:, myzb); Dat0 => Dat(:, myzb -1)
!    CALL SENDRECV(Dat(:, myzb-1), n, ipartner,.false., Comm)   
!    CALL SENDRECV(Dat(:, myzb), n, ipartner, .true., Comm)
!  ENDIF

END SUBROUTINE

SUBROUTINE GetNeighDat_INT3D(Dat, n1, n2, myzb, myze, myrank, nProc, Comm)
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
INTEGER, TARGET :: Dat(n1, n2, myzb-1:myze+1)
INTEGER, POINTER :: Dat0(:, :), Dat1(:, :)
INTEGER :: n1, n2, myzb, myze, myrank, nproc, Comm

INTEGER :: i, j, istep, igrp, ipartner
LOGICAL :: lSend

DO istep = 0, 1 
  i = MOD(myrank+istep, 2)
  IF(i .EQ. 0) THEN
    ipartner = myrank + 1; lSend = .TRUE.
    IF( ipartner .LT. 0 .OR. ipartner .GT. nproc-1) CYCLE    
    CALL CP_CA(Dat(:, :, myze+1), 0, n1, n2)
    CALL SENDRECV(Dat(:, :, myze), n1, n2, ipartner, .true., Comm)
    CALL SENDRECV(Dat(:, :, myze+1), n1, n2, ipartner,.false., Comm)    
    !Dat0 => Dat(:, :, myze); Dat1 => Dat(:, :, myze + 1)
  ENDIF
  IF(i .EQ. 1) THEN
    ipartner = myrank - 1; lSend = .FALSE.
    IF( ipartner .LT. 0 .OR. ipartner .GT. nproc-1) CYCLE    
    CALL CP_CA(Dat(:, :, myzb-1), 0, n1, n2)
    CALL SENDRECV(Dat(:, :, myzb-1), n1, n2, ipartner,.false., Comm)   
    CALL SENDRECV(Dat(:, :, myzb), n1, n2, ipartner, .true., Comm)    
   !Dat1 => Dat(:, :, myzb); Dat0 => Dat(:, :, myzb -1)
  ENDIF
ENDDO
END SUBROUTINE

SUBROUTINE GetNeighDat_REAL2D(Dat, n, myzb, myze, myrank, nProc, Comm)
USE BasicOperation, ONLY : CP_CA
REAL, TARGET :: Dat(n, myzb-1:myze+1)
REAL, POINTER :: Dat0(:), Dat1(:)
INTEGER :: n, myzb, myze, myrank, nproc, Comm

INTEGER :: i, j, istep, igrp, ipartner
LOGICAL :: lSend

DO istep = 0, 1 
  i = MOD(myrank+istep, 2)
  IF(i .EQ. 0) THEN
    ipartner = myrank + 1; lSend = .TRUE.
    CALL CP_CA(Dat(:, myze+1), 0._8, n)
    IF( ipartner .LT. 0 .OR. ipartner .GT. nproc-1) CYCLE    
    !Dat0 => Dat(:, myze); Dat1 => Dat(:, myze + 1)
    CALL SENDRECV(Dat(:, myze), n, ipartner, .true., Comm)
    CALL SENDRECV(Dat(:, myze+1), n, ipartner,.false., Comm)
  ELSE
    ipartner = myrank - 1; lSend = .FALSE.
    CALL CP_CA(Dat(:, myzb-1), 0._8, n)
    IF( ipartner .LT. 0 .OR. ipartner .GT. nproc-1) CYCLE    
    !Dat1 => Dat(:, myzb); Dat0 => Dat(:, myzb -1)
    CALL SENDRECV(Dat(:, myzb-1), n, ipartner,.false., Comm)   
    CALL SENDRECV(Dat(:, myzb), n, ipartner, .true., Comm)
  ENDIF
ENDDO

END SUBROUTINE

SUBROUTINE GetNeighDat_INT2D(Dat, n, myzb, myze, myrank, nProc, Comm)
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
INTEGER, TARGET :: Dat(n, myzb-1:myze+1)
INTEGER, POINTER :: Dat0(:), Dat1(:)
INTEGER :: n, myzb, myze, myrank, nproc, Comm

INTEGER :: i, j, istep, igrp, ipartner
LOGICAL :: lSend

DO istep = 0, 1 
  i = MOD(myrank+istep, 2)
  IF(i .EQ. 0) THEN
    ipartner = myrank + 1; lSend = .TRUE.
    CALL CP_CA(Dat(:, myze+1), 0, n)
    IF( ipartner .LT. 0 .OR. ipartner .GT. nproc-1) CYCLE    
    !Dat0 => Dat(:, myze); Dat1 => Dat(:, myze + 1)
    CALL SENDRECV(Dat(:, myze), n, ipartner, .true., Comm)
    CALL SENDRECV(Dat(:, myze+1), n, ipartner,.false., Comm)
  ELSE
    ipartner = myrank - 1; lSend = .FALSE.
    CALL CP_CA(Dat(:, myzb-1), 0, n)
    IF( ipartner .LT. 0 .OR. ipartner .GT. nproc-1) CYCLE    
    !Dat1 => Dat(:, myzb); Dat0 => Dat(:, myzb -1)
    CALL SENDRECV(Dat(:, myzb-1), n, ipartner,.false., Comm)   
    CALL SENDRECV(Dat(:, myzb), n, ipartner, .true., Comm)
  ENDIF
ENDDO
END SUBROUTINE

SUBROUTINE GetNeighDat_REAL1D(Dat, myzb, myze, myrank, nProc, Comm)
REAL, TARGET :: Dat(myzb-1:myze+1)
REAL, POINTER :: Dat0, Dat1
INTEGER :: myzb, myze, myrank, nproc, Comm

INTEGER :: i, j, istep, igrp, ipartner
LOGICAL :: lSend

DO istep = 0, 1 
  i = MOD(myrank+istep, 2)
  IF(i .EQ. 0) THEN
    ipartner = myrank + 1; lSend = .TRUE.
    Dat0 => Dat(myze); Dat1 => Dat(myze + 1)
  ENDIF
  IF(i .EQ. 1) THEN
    ipartner = myrank - 1; lSend = .FALSE.
    Dat1 => Dat(myzb); Dat0 => Dat(myzb -1)
  ENDIF
  IF( ipartner .LT. 0 .OR. ipartner .GT. nproc-1) CYCLE
  CALL SENDRECV(Dat0, ipartner, lSend, Comm)
  CALL SENDRECV(Dat1, ipartner, .NOT. lSend, Comm)
ENDDO
END SUBROUTINE

SUBROUTINE GetNeighDat_INT1D(Dat, myzb, myze, myrank, nProc, Comm)
INTEGER, TARGET :: Dat(myzb-1:myze+1)
INTEGER, POINTER :: Dat0, Dat1
INTEGER :: myzb, myze, myrank, nproc, Comm

INTEGER :: i, j, istep, igrp, ipartner
LOGICAL :: lSend

DO istep = 0, 1 
  i = MOD(myrank+istep, 2)
  IF(i .EQ. 0) THEN
    ipartner = myrank + 1; lSend = .TRUE.
    Dat0 => Dat(myze); Dat1 => Dat(myze + 1)
  ENDIF
  IF(i .EQ. 1) THEN
    ipartner = myrank - 1; lSend = .FALSE.
    Dat1 => Dat(myzb); Dat0 => Dat(myzb -1)
  ENDIF
  IF( ipartner .LT. 0 .OR. ipartner .GT. nproc-1) CYCLE
  CALL SENDRECV(Dat0, ipartner, lSend, Comm)
  CALL SENDRECV(Dat1, ipartner, .NOT. lSend, Comm)
ENDDO
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GetMPIFile(Dat, Master, comm)
IMPLICIT NONE
CHARACTER(80) :: Dat 
LOGICAL :: Master
INTEGER:: comm, ierr
REAL:: Tbeg, Tend

tbeg = nTracer_dclock(.false., .false.)

CALL MPI_BCAST(DAT, 80, MPI_CHARACTER, 0, comm, ierr)

tend = nTracer_dclock(.false., .false.)
TimeChk%CommTime = TimeChk%CommTime + tend - tbeg

END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE MPI_SYNC(comm)
IMPLICIT NONE
INTEGER :: comm, ierr
REAL:: Tbeg, Tend

tbeg = nTracer_dclock(.false., .false.)
CALL MPI_BARRIER(comm, ierr)
tend = nTracer_dclock(.false., .false.)
TimeChk%CommTime = TimeChk%CommTime + tend - tbeg

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MPIWaitTurn(comm, myrank, nproc, lPass)
IMPLICIT NONE
INTEGER :: COMM, nProc, Myrank
LOGICAL :: lPass
INTEGER :: i, j, ierr, ito, ifr, istat
REAL:: Tbeg, Tend, CommT
IF(.NOT. lPass) THEN
  CALL MPI_BARRIER(Comm, ierr)
  i=0
  DO j = 0, myrank
    !CALL MPI_BCAST(i, 1, MPI_INTEGER, j, comm, ierr)
    CALL MPI_BARRIER(Comm, ierr)
  ENDDO
ELSE
  DO j = myrank+1, nproc-1
    !CALL MPI_BCAST(i, 1, MPI_INTEGER, j, comm, ierr)
    CALL MPI_BARRIER(Comm, ierr)
  ENDDO
  CALL MPI_BARRIER(Comm, ierr)
ENDIF
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE MPI_MAX_REAL(dat, comm, lall)
REAL :: dat
INTEGER :: comm
LOGICAL :: lall
INTEGER :: ierr
REAL :: buf
REAL:: Tbeg, Tend

Tbeg = nTracer_dclock(.FALSE., .FALSE.)
IF(.NOT. lAll) THEN
  CALL MPI_REDUCE(dat, Buf, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)
ELSE
  CALL MPI_ALLREDUCE(dat, Buf, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)
ENDIF
dat = buf
Tend = nTracer_dclock(.FALSE., .FALSE.)
TimeChk%CommTime = TimeChk%CommTime + (tend - tbeg)
END SUBROUTINE

SUBROUTINE MPI_MAX_INT(dat, comm, lall)
INTEGER :: dat
INTEGER :: comm
LOGICAL :: lall
INTEGER :: ierr
INTEGER :: buf
REAL:: Tbeg, Tend

Tbeg = nTracer_dclock(.FALSE., .FALSE.)
IF(.NOT. lAll) THEN
  CALL MPI_REDUCE(dat, Buf, 1, MPI_INTEGER, MPI_MAX, 0, comm, ierr)
ELSE
  CALL MPI_ALLREDUCE(dat, Buf, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
ENDIF
dat = buf
Tend = nTracer_dclock(.FALSE., .FALSE.)
TimeChk%CommTime = TimeChk%CommTime + (tend - tbeg)
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE REDUCEnBCAST_Real0D(dat, Comm1, Comm2, lMember1, lMember2, lall)
IMPLICIT NONE
REAL :: dat
INTEGER :: comm1, comm2
LOGICAL :: lMember1, lMember2, lall
REAL :: temp
REAL:: Tbeg, Tend
!!!Master of comm1 and Master of Comm2 should be identical
IF(lMember1) THEN
  CALL REDUCE(dat, temp, comm1, lall)
ENDIF

IF(lMember2) THEN
  CALL BCAST(temp, comm2)
  dat = temp
ENDIF
END SUBROUTINE

SUBROUTINE REDUCEnBCAST_INT0D(dat, Comm1, Comm2, lMember1, lMember2, lall)
IMPLICIT NONE
INTEGER :: dat
INTEGER :: comm1, comm2
LOGICAL :: lMember1, lMember2, lall
INTEGER :: temp
REAL:: Tbeg, Tend
!!!Master of comm1 and Master of Comm2 should be identical
IF(lMember1) THEN
  CALL REDUCE(dat, temp, comm1, lall)
ENDIF

IF(lMember2) THEN
  CALL BCAST(temp, comm2)
  dat = temp
ENDIF

END SUBROUTINE

FUNCTION MPIConvChk(lconv, nproc, comm)
IMPLICIT NONE
LOGICAL :: lConv, MPIConvChk
INTEGER :: nproc, comm
REAL :: Buf0, Buf
Buf0 = 0
IF(.NOT. lConv) Buf0 = 1
CALL Reduce(Buf0, Buf, Comm, .TRUE.)
MPIConvChk = .TRUE.
IF(Buf .NE. 0) MPIConvChk = .FALSE.
END FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MPI_SeqMessage(mesg0, io, ltime, ldisp, myrank, NPROC, COMM)
USE PARAM
USE IOUTIL,   ONLY : message
IMPLICIT NONE
CHARACTER(132) :: mesg0
INTEGER :: IO, MYRANK, COMM, nproc
LOGICAL :: ltime, ldisp
INTEGER :: i, ierr, istatus(MPI_STATUS_SIZE)

IF(myrank .EQ. 0) THEN
  CALL message(io, ltime, ldisp, mesg0)
  DO i = 1, nproc-1
    CALL MPI_RECV(mesg(1:132), 132, MPI_CHARACTER, i, 2, comm, istatus, ierr)
    mesg(14:132) = mesg(1:119); mesg(1:13) =''
    CALL message(io, false, ldisp, mesg)
    
  ENDDO
ELSE
  CALL MPI_SEND(mesg0(1:132), 132, MPI_CHARACTER, 0, 2, comm, ierr)
ENDIF
CALL MPI_SYNC(Comm)
END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE

MODULE MPIGetNeighbor

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, PRIVATE :: Request(4), nRequest
INTEGER, PARAMETER :: bottom = 1, top = 2

CONTAINS

!--- Non-blocking Communication Routines

SUBROUTINE InitFastComm()

IMPLICIT NONE

nRequest = 0

END SUBROUTINE

SUBROUTINE GetNeighborFast(n, send, recv, dir, rank, comm, nproc)

IMPLICIT NONE

REAL :: send(*), recv(*)
INTEGER :: n, nproc, dir, rank
INTEGER :: i, ipartner
INTEGER :: tag = 1, info, comm

SELECT CASE (dir)
CASE (top)
  ipartner = rank + 1
  IF (ipartner .LT. nproc) THEN
    nRequest = nRequest + 1
    CALL MPI_ISEND(send, n, MPI_DOUBLE_PRECISION, ipartner, tag, comm, Request(nRequest), info)
    nRequest = nRequest + 1
    CALL MPI_IRECV(recv, n, MPI_DOUBLE_PRECISION, ipartner, tag, comm, Request(nRequest), info)
  ENDIF
CASE (bottom)
  ipartner = rank - 1
  IF (ipartner .GE. 0) THEN
    nRequest = nRequest + 1
    CALL MPI_ISEND(send, n, MPI_DOUBLE_PRECISION, ipartner, tag, comm, Request(nRequest), info)
    nRequest = nRequest + 1
    CALL MPI_IRECV(recv, n, MPI_DOUBLE_PRECISION, ipartner, tag, comm, Request(nRequest), info)
  ENDIF
END SELECT

END SUBROUTINE

SUBROUTINE FinalizeFastComm(comm)
IMPLICIT NONE
INTEGER :: comm

INTEGER :: status(MPI_STATUS_SIZE, 4), info

CALL MPI_WAITALL(nRequest, Request, status, info)
CALL MPI_BARRIER(comm, info)

END SUBROUTINE
   
END MODULE