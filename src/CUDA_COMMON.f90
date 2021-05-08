#include <CUDADEFINES.h>
#include <Depletion.h>

#ifdef __PGI

MODULE CUDA_CONST

IMPLICIT NONE

INTEGER, CONSTANT :: nfsr, nFxr, nxy, nxyc, nz, ng, nofg, norg, nchi
INTEGER, CONSTANT :: nRotRay, nModRay, nAziAngle, nPolarAngle, nPolar1D, nPhiAngSv, nMoment, ScatOd, nAziMap
INTEGER, CONSTANT :: AziMap(GPU_MAX_AZI, 2), dir(2, 2), inc(2)
INTEGER, CONSTANT :: InScatRange(2, GPU_MAX_GROUP), BegGrpScat(GPU_MAX_GROUP+1)
REAL(GPU_PRECISION), CONSTANT :: rsinv(GPU_MAX_POLAR), sinv(GPU_MAX_POLAR)
REAL(GPU_PRECISION), CONSTANT :: rcosv1D(GPU_MAX_POLAR_1D), cosv1D(GPU_MAX_POLAR_1D)
REAL(GPU_PRECISION), CONSTANT :: wt(GPU_MAX_POLAR, GPU_MAX_AZI), wt1D(GPU_MAX_POLAR_1D)
REAL(GPU_PRECISION), CONSTANT :: wtsurf(4, GPU_MAX_POLAR, GPU_MAX_AZI), wtsurf1D(GPU_MAX_POLAR_1D)
REAL(GPU_PRECISION), CONSTANT :: Comp(GPU_MAX_ORDER, GPU_MAX_POLAR, GPU_MAX_AZI), Comp1D(GPU_MAX_POLAR_1D)
REAL(GPU_PRECISION), CONSTANT :: mwt(GPU_MAX_ORDER, GPU_MAX_POLAR, GPU_MAX_AZI), mwt1D(GPU_MAX_POLAR_1D)

INTEGER(4), CONSTANT :: rowptr(1024), colIdx(4096), eyeIdx(1024)
INTEGER, CONSTANT :: CRAM_Order
#ifdef CRAM_14
COMPLEX(8), CONSTANT :: Pole(7)
COMPLEX(8), CONSTANT :: Res(7)
#else CRAM_16
COMPLEX(8), CONSTANT :: Pole(8)
COMPLEX(8), CONSTANT :: Res(8)
#endif
REAL(8), CONSTANT :: Res0

INTEGER, CONSTANT :: nprec
REAL(8), CONSTANT :: theta
REAL(8), CONSTANT :: chid(GPU_MAX_GROUP)
REAL(8), CONSTANT :: chidk(GPU_MAX_GROUP, 6)
INTEGER, CONSTANT :: ifxrbegD(32), ifxrbegT(32)
LOGICAL, CONSTANT :: lFuelPlane(64)
END MODULE

MODULE CUDA_MASTER

USE CUDATypeDef
USE OPENACC
USE OMP_LIB
IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, PARAMETER :: in = 1, out = 2, surf = 3
INTEGER, PARAMETER :: bottom = 1, top = 2
INTEGER, PARAMETER :: red = 1, black = 2
INTEGER, PARAMETER :: Void = 0, Reflective = -1
INTEGER, PARAMETER :: SOUTH = 1, WEST = 2, NORTH = 3, EAST = 4
INTEGER, PARAMETER :: NODAL = 1, MOC = 3

TYPE(cuGeometry_Type) :: cuGeometry

TYPE(cuMOC_Type) :: cuMOC
TYPE(cuCMFD_Type) :: cuCMFD, cuGcCMFD
TYPE(cuCMFD_TYPE) :: cuCMFD_adj, cuGcCMFD_adj
TYPE(cuAxial_Type) :: cuAxial

TYPE(cuTranCMInfo_Type) :: cuTranCMInfo
TYPE(cuPwDist_Type) :: cuPwDist

TYPE(cuDevice_Type) :: cuDevice
TYPE(cuCntl_Type) :: cuCntl

TYPE(cuRotRay1D_Type), POINTER :: cuFastRay1D(:)
TYPE(cuDepl_Type) :: cuDepl

INTEGER :: MPI_CUDA_COMM, MPI_CUDA_RANK, NUM_CUDA_PROC
INTEGER :: MPI_CUDA_SHARED_COMM, MPI_CUDA_SHARED_RANK, NUM_CUDA_SHARED_PROC

END MODULE

MODULE CUDA_UTIL

USE CUDA_MASTER
IMPLICIT NONE

!--- Global Communication Buffers

INTEGER :: nRequest, Request(4)

REAL(4), ALLOCATABLE, PINNED, PRIVATE :: reduceBuffer4(:)
REAL(8), ALLOCATABLE, PINNED, PRIVATE :: reduceBuffer8(:)

INTERFACE MPIComm
  MODULE PROCEDURE MPIComm4
  MODULE PROCEDURE MPIComm8
  MODULE PROCEDURE MPICommHost4
  MODULE PROCEDURE MPICommHost8
END INTERFACE

INTERFACE cuReduce
  MODULE PROCEDURE cuReduceFloat
  MODULE PROCEDURE cuReduceDouble
  MODULE PROCEDURE cuReduceHostFloat
  MODULE PROCEDURE cuReduceHostDouble
END INTERFACE

INTERFACE cuGetNeighbor
  MODULE PROCEDURE cuGetNeighbor4
  MODULE PROCEDURE cuGetNeighbor8
  MODULE PROCEDURE cuGetNeighborHost4
  MODULE PROCEDURE cuGetNeighborHost8
  MODULE PROCEDURE cuGetNeighborDirHost4
  MODULE PROCEDURE cuGetNeighborDirHost8
END INTERFACE

INTERFACE cuMatMul3D
  MODULE PROCEDURE cuMatMulFloat
  MODULE PROCEDURE cuMatMulFloatRB
  MODULE PROCEDURE cuMatMulDouble
  MODULE PROCEDURE cuMatMulDoubleRB
  MODULE PROCEDURE cuMatMulMixed
  MODULE PROCEDURE cuMatMulMixedRB
END INTERFACE

INTERFACE normalizeMulti
  MODULE PROCEDURE normalizeMulti4
  MODULE PROCEDURE normalizeMulti8
END INTERFACE

INTERFACE normMulti
  MODULE PROCEDURE normMulti4
  MODULE PROCEDURE normMulti8
END INTERFACE

INTERFACE dotMulti
  MODULE PROCEDURE dotMulti4
  MODULE PROCEDURE dotMulti8
END INTERFACE

INTERFACE asumMulti
  MODULE PROCEDURE asumMulti4
  MODULE PROCEDURE asumMulti8
END INTERFACE

INTERFACE cuVectorOp
  MODULE PROCEDURE cuVectorOp4
  MODULE PROCEDURE cuVectorOp8
END INTERFACE

INTERFACE cuTypecast
  MODULE PROCEDURE cuHtoDCopy8
  MODULE PROCEDURE cuDtoHCopy8
  MODULE PROCEDURE cuHtoDTypecastFrom8
  MODULE PROCEDURE cuDtoHTypecastFrom4
  MODULE PROCEDURE cuTypecastFrom4
  MODULE PROCEDURE cuTypecastFrom8
END INTERFACE

INTERFACE cuInitArray
  MODULE PROCEDURE cuInitArray4
  MODULE PROCEDURE cuInitArray8
END INTERFACE

CONTAINS

!--- Non-Blocking MPI Communication Routines

SUBROUTINE InitMPIComm()

IMPLICIT NONE

nRequest = 0

END SUBROUTINE

SUBROUTINE MPIComm4(ns, nr, send, recv, dir, comm)
USE PE_MOD,         ONLY : PE
IMPLICIT NONE

REAL(4), DEVICE :: send(*), recv(*)
INTEGER :: ns, nr, dir, comm
INTEGER :: irank, ipartner
INTEGER :: tag = 1, info

irank = MPI_CUDA_RANK

SELECT CASE (dir)
CASE (bottom)
  ipartner = irank - 1
  IF (ipartner .GE. 0) THEN
    nRequest = nRequest + 1
    CALL MPI_ISEND(send, ns, MPI_FLOAT, ipartner, tag, comm, Request(nRequest), info)
    nRequest = nRequest + 1
    CALL MPI_IRECV(recv, nr, MPI_FLOAT, ipartner, tag, comm, Request(nRequest), info)
  ENDIF
CASE (top)
  ipartner = irank + 1
  IF (ipartner .LT. NUM_CUDA_PROC) THEN
    nRequest = nRequest + 1
    CALL MPI_ISEND(send, ns, MPI_FLOAT, ipartner, tag, comm, Request(nRequest), info)
    nRequest = nRequest + 1
    CALL MPI_IRECV(recv, nr, MPI_FLOAT, ipartner, tag, comm, Request(nRequest), info)
  ENDIF
END SELECT

END SUBROUTINE

SUBROUTINE MPIComm8(ns, nr, send, recv, dir, comm)
USE PE_MOD,         ONLY : PE
IMPLICIT NONE

REAL(8), DEVICE :: send(*), recv(*)
INTEGER :: ns, nr, dir, comm
INTEGER :: irank, ipartner
INTEGER :: tag = 1, info

irank = MPI_CUDA_RANK

SELECT CASE (dir)
CASE (bottom)
  ipartner = irank - 1
  IF (ipartner .GE. 0) THEN
    nRequest = nRequest + 1
    CALL MPI_ISEND(send, ns, MPI_DOUBLE_PRECISION, ipartner, tag, comm, Request(nRequest), info)
    nRequest = nRequest + 1
    CALL MPI_IRECV(recv, nr, MPI_DOUBLE_PRECISION, ipartner, tag, comm, Request(nRequest), info)
  ENDIF
CASE (top)
  ipartner = irank + 1
  IF (ipartner .LT. NUM_CUDA_PROC) THEN
    nRequest = nRequest + 1
    CALL MPI_ISEND(send, ns, MPI_DOUBLE_PRECISION, ipartner, tag, comm, Request(nRequest), info)
    nRequest = nRequest + 1
    CALL MPI_IRECV(recv, nr, MPI_DOUBLE_PRECISION, ipartner, tag, comm, Request(nRequest), info)
  ENDIF
END SELECT

END SUBROUTINE

SUBROUTINE MPICommHost4(ns, nr, send, recv, dir, comm)
USE PE_MOD,         ONLY : PE
IMPLICIT NONE

REAL(4) :: send(*), recv(*)
INTEGER :: ns, nr, dir, comm
INTEGER :: irank, ipartner
INTEGER :: tag = 1, info

irank = MPI_CUDA_RANK

SELECT CASE (dir)
CASE (bottom)
  ipartner = irank - 1
  IF (ipartner .GE. 0) THEN
    nRequest = nRequest + 1
    CALL MPI_ISEND(send, ns, MPI_FLOAT, ipartner, tag, comm, Request(nRequest), info)
    nRequest = nRequest + 1
    CALL MPI_IRECV(recv, nr, MPI_FLOAT, ipartner, tag, comm, Request(nRequest), info)
  ENDIF
CASE (top)
  ipartner = irank + 1
  IF (ipartner .LT. NUM_CUDA_PROC) THEN
    nRequest = nRequest + 1
    CALL MPI_ISEND(send, ns, MPI_FLOAT, ipartner, tag, comm, Request(nRequest), info)
    nRequest = nRequest + 1
    CALL MPI_IRECV(recv, nr, MPI_FLOAT, ipartner, tag, comm, Request(nRequest), info)
  ENDIF
END SELECT

END SUBROUTINE

SUBROUTINE MPICommHost8(ns, nr, send, recv, dir, comm)
USE PE_MOD,         ONLY : PE
IMPLICIT NONE

REAL(8) :: send(*), recv(*)
INTEGER :: ns, nr, dir, comm
INTEGER :: i, irank, ipartner
INTEGER :: tag = 1, info

irank = MPI_CUDA_RANK

SELECT CASE (dir)
CASE (bottom)
  ipartner = irank - 1
  IF (ipartner .GE. 0) THEN
    nRequest = nRequest + 1
    CALL MPI_ISEND(send, ns, MPI_DOUBLE_PRECISION, ipartner, tag, comm, Request(nRequest), info)
    nRequest = nRequest + 1
    CALL MPI_IRECV(recv, nr, MPI_DOUBLE_PRECISION, ipartner, tag, comm, Request(nRequest), info)
  ENDIF
CASE (top)
  ipartner = irank + 1
  IF (ipartner .LT. NUM_CUDA_PROC) THEN
    nRequest = nRequest + 1
    CALL MPI_ISEND(send, ns, MPI_DOUBLE_PRECISION, ipartner, tag, comm, Request(nRequest), info)
    nRequest = nRequest + 1
    CALL MPI_IRECV(recv, nr, MPI_DOUBLE_PRECISION, ipartner, tag, comm, Request(nRequest), info)
  ENDIF
END SELECT

END SUBROUTINE

SUBROUTINE FinalizeMPIComm()

IMPLICIT NONE

INTEGER :: status(MPI_STATUS_SIZE, 4), info

CALL MPI_WAITALL(nRequest, Request, status, info)

END SUBROUTINE

!--- Multi-Device Reduction Routines

SUBROUTINE cuReduceFloat(n, arr, stream, comm, lMPI)

IMPLICIT NONE

INTEGER :: n, comm
INTEGER(KIND = cuda_stream_kind) :: stream
REAL(4), DEVICE :: arr(*)
LOGICAL :: lMPI

INTEGER :: ierr
REAL(4), ALLOCATABLE, PINNED :: myBuffer(:)

!$OMP MASTER
ALLOCATE(reduceBuffer4(n)); reduceBuffer4 = 0.0
!$OMP END MASTER

ALLOCATE(myBuffer(n))

ierr = cudaMemcpyAsync(myBuffer, arr, n, cudaMemcpyDeviceToHost, stream)
ierr = cudaStreamSynchronize(stream)

!$OMP CRITICAL
reduceBuffer4 = reduceBuffer4 + myBuffer
!$OMP END CRITICAL

IF (lMPI) THEN
  !$OMP MASTER
  myBuffer = reduceBuffer4
  CALL MPI_ALLREDUCE(myBuffer, reduceBuffer4, n, MPI_FLOAT, MPI_SUM, comm, ierr)
  !$OMP END MASTER
ENDIF

ierr = cudaMemcpyAsync(arr, reduceBuffer4, n, cudaMemcpyHostToDevice, stream)
ierr = cudaStreamSynchronize(stream)

DEALLOCATE(myBuffer)

!$OMP MASTER
DEALLOCATE(reduceBuffer4)
!$OMP END MASTER

END SUBROUTINE

SUBROUTINE cuReduceDouble(n, arr, stream, comm, lMPI)

IMPLICIT NONE

INTEGER :: n, comm
INTEGER(KIND = cuda_stream_kind) :: stream
REAL(8), DEVICE :: arr(*)
LOGICAL :: lMPI

INTEGER :: ierr
REAL(8), ALLOCATABLE, PINNED :: myBuffer(:)

!$OMP MASTER
ALLOCATE(reduceBuffer8(n)); reduceBuffer8 = 0.0
!$OMP END MASTER

ALLOCATE(myBuffer(n))

ierr = cudaMemcpyAsync(myBuffer, arr, n, cudaMemcpyDeviceToHost, stream)
ierr = cudaStreamSynchronize(stream)

!$OMP CRITICAL
reduceBuffer8 = reduceBuffer8 + myBuffer
!$OMP END CRITICAL

IF (lMPI) THEN
  !$OMP MASTER
  myBuffer = reduceBuffer8
  CALL MPI_ALLREDUCE(myBuffer, reduceBuffer8, n, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
  !$OMP END MASTER
ENDIF

ierr = cudaMemcpyAsync(arr, reduceBuffer8, n, cudaMemcpyHostToDevice, stream)
ierr = cudaStreamSynchronize(stream)

DEALLOCATE(myBuffer)

!$OMP MASTER
DEALLOCATE(reduceBuffer8)
!$OMP END MASTER

END SUBROUTINE

SUBROUTINE cuReduceHostFloat(n, arrDev, arrHost, stream, comm, lMPI)

IMPLICIT NONE

INTEGER :: n, comm
INTEGER(KIND = cuda_stream_kind) :: stream
REAL(4), DEVICE :: arrDev(*)
REAL(4) :: arrHost(*)
LOGICAL :: lMPI

INTEGER :: ierr
REAL(4), ALLOCATABLE, PINNED :: myBuffer(:)

!$OMP MASTER
ALLOCATE(reduceBuffer4(n)); reduceBuffer4 = 0.0
!$OMP END MASTER

ALLOCATE(myBuffer(n))

ierr = cudaMemcpyAsync(myBuffer, arrDev, n, cudaMemcpyDeviceToHost, stream)
ierr = cudaStreamSynchronize(stream)

!$OMP CRITICAL
reduceBuffer4 = reduceBuffer4 + myBuffer
!$OMP END CRITICAL

IF (lMPI) THEN
  !$OMP MASTER
  myBuffer = reduceBuffer4
  CALL MPI_ALLREDUCE(myBuffer, reduceBuffer4, n, MPI_FLOAT, MPI_SUM, comm, ierr)
  !$OMP END MASTER
ENDIF

arrHost(1 : n) = reduceBuffer4

DEALLOCATE(myBuffer)

!$OMP MASTER
DEALLOCATE(reduceBuffer4)
!$OMP END MASTER

END SUBROUTINE

SUBROUTINE cuReduceHostDouble(n, arrDev, arrHost, stream, comm, lMPI)

IMPLICIT NONE

INTEGER :: n, comm
INTEGER(KIND = cuda_stream_kind) :: stream
REAL(8), DEVICE :: arrDev(*)
REAL(8) :: arrHost(*)
LOGICAL :: lMPI

INTEGER :: ierr
REAL(8), ALLOCATABLE, PINNED :: myBuffer(:)

!$OMP MASTER
ALLOCATE(reduceBuffer8(n)); reduceBuffer8 = 0.0
!$OMP END MASTER

ALLOCATE(myBuffer(n))

ierr = cudaMemcpyAsync(myBuffer, arrDev, n, cudaMemcpyDeviceToHost, stream)
ierr = cudaStreamSynchronize(stream)

!$OMP CRITICAL
reduceBuffer8 = reduceBuffer8 + myBuffer
!$OMP END CRITICAL

IF (lMPI) THEN
  !$OMP MASTER
  myBuffer = reduceBuffer8
  CALL MPI_ALLREDUCE(myBuffer, reduceBuffer8, n, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
  !$OMP END MASTER
ENDIF

arrHost(1 : n) = reduceBuffer8

DEALLOCATE(myBuffer)

!$OMP MASTER
DEALLOCATE(reduceBuffer8)
!$OMP END MASTER

END SUBROUTINE

!--- Multi-Device Data Exchange Routines

SUBROUTINE cuGetNeighbor4(n, send, recv, comm, stream)

IMPLICIT NONE

INTEGER :: n, comm
INTEGER(KIND = cuda_stream_kind) :: stream
REAL(4), DEVICE :: send(*), recv(*)

IF (.NOT. cuCntl%lMulti) RETURN

CALL InitMPIComm()
CALL MPIComm(n, n, send(1 : n), recv(1 : n), bottom, comm)
CALL MPIComm(n, n, send(n + 1 : n + n), recv(n + 1 : n + n), top, comm)
CALL FinalizeMPIComm()

END SUBROUTINE

SUBROUTINE cuGetNeighbor8(n, send, recv, comm, stream)

IMPLICIT NONE

INTEGER :: n, comm
INTEGER(KIND = cuda_stream_kind) :: stream
REAL(8), DEVICE :: send(*), recv(*)

IF (.NOT. cuCntl%lMulti) RETURN

CALL InitMPIComm()
CALL MPIComm(n, n, send(1 : n), recv(1 : n), bottom, comm)
CALL MPIComm(n, n, send(n + 1 : n + n), recv(n + 1 : n + n), top, comm)
CALL FinalizeMPIComm()

END SUBROUTINE

SUBROUTINE cuGetNeighborHost4(n, send, recv, comm)

IMPLICIT NONE

INTEGER :: n, comm
REAL(4) :: send(*), recv(*)

IF (.NOT. cuCntl%lMulti) RETURN

CALL InitMPIComm()
CALL MPIComm(n, n, send(1 : n), recv(1 : n), bottom, comm)
CALL MPIComm(n, n, send(n + 1 : n + n), recv(n + 1 : n + n), top, comm)
CALL FinalizeMPIComm()

END SUBROUTINE

SUBROUTINE cuGetNeighborHost8(n, send, recv, comm)

IMPLICIT NONE

INTEGER :: n, comm
REAL(8) :: send(*), recv(*)

IF (.NOT. cuCntl%lMulti) RETURN

CALL InitMPIComm()
CALL MPIComm(n, n, send(1 : n), recv(1 : n), bottom, comm)
CALL MPIComm(n, n, send(n + 1 : n + n), recv(n + 1 : n + n), top, comm)
CALL FinalizeMPIComm()

END SUBROUTINE

SUBROUTINE cuGetNeighborDirHost4(n, send, recv, comm, dir)

IMPLICIT NONE

INTEGER :: n, comm, dir
REAL(4) :: send(*), recv(*)

IF (.NOT. cuCntl%lMulti) RETURN

CALL InitMPIComm()
CALL MPIComm(n, n, send, recv, dir, comm)
CALL FinalizeMPIComm()

END SUBROUTINE

SUBROUTINE cuGetNeighborDirHost8(n, send, recv, comm, dir)

IMPLICIT NONE

INTEGER :: n, comm, dir
REAL(8) :: send(*), recv(*)

IF (.NOT. cuCntl%lMulti) RETURN

CALL InitMPIComm()
CALL MPIComm(n, n, send, recv, dir, comm)
CALL FinalizeMPIComm()

END SUBROUTINE

!--- Multi-Device Matrix Multiplication Routines

SUBROUTINE cuMatMulFloat(Diag, offDiag, x, y, br, tr, sparseHandle, stream, comm)

IMPLICIT NONE

TYPE(CSR_MIXED) :: Diag
REAL(4), DEVICE :: offDiag(*), x(*), y(*)
INTEGER :: br(2), tr(2)
TYPE(cusparseHandle) :: sparseHandle
INTEGER(KIND = cuda_stream_kind) :: stream
INTEGER :: comm

REAL(4), ALLOCATABLE, DEVICE :: xNeigh(:), yTemp(:)
REAL(4), POINTER, DEVICE :: csrVal(:)
INTEGER, POINTER, DEVICE :: csrRowPtr(:), csrColIdx(:)
INTEGER :: n, nb, nt, nr, nc, nnz
INTEGER :: ierr

csrVal => Diag%d_csrVal
csrRowPtr => Diag%d_csrRowPtr
csrColIdx => Diag%d_csrColIdx
nr = Diag%nr
nc = Diag%nc
nnz = Diag%nnz
nb = br(2) - br(1) + 1
nt = tr(2) - tr(1) + 1
n = max(nb, nt)

IF (.NOT. cuCntl%lMulti .OR. cuCntl%lDcpl) THEN

  ierr = cusparseScsrmv(sparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, 1.0_4, Diag%descr,             &
                        csrVal, csrRowPtr, csrColIdx, x, 0.0_4, y)
  RETURN

ENDIF

ALLOCATE(xNeigh(nb + nt))
ALLOCATE(yTemp(nb + nt))

CALL cuInitArray(nb + nt, xNeigh, stream)

ierr = cudaStreamSynchronize(stream)

CALL InitMPIComm()
CALL MPIComm(nb, nb, x(br(1) : br(2)), xNeigh(1 : nb), bottom, comm)
CALL MPIComm(nt, nt, x(tr(1) : tr(2)), xNeigh(nb + 1 : nb + nt), top, comm)

ierr = cusparseScsrmv(sparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, 1.0_4, Diag%descr,               &
                      csrVal, csrRowPtr, csrColIdx, x, 0.0_4, y)

CALL FinalizeMPIComm()

CALL cuVectorOp('*', nb + nt, offDiag, xNeigh, yTemp, stream)
CALL cuVectorOp('+', nb, yTemp(1 : nb), y(br(1) : br(2)), y(br(1) : br(2)), stream)
CALL cuVectorOp('+', nt, yTemp(nb + 1 : nb + nt), y(tr(1) : tr(2)), y(tr(1) : tr(2)), stream)

ierr = cudaStreamSynchronize(stream)

DEALLOCATE(xNeigh)
DEALLOCATE(yTemp)

END SUBROUTINE

SUBROUTINE cuMatMulFloatRB(Diag, offDiag, x, y, br, tr, sparseHandle, stream, comm)

IMPLICIT NONE

TYPE(CSR_MIXED) :: Diag
REAL(4), DEVICE :: offDiag(*), x(*), y(*)
INTEGER :: br(2, 2), tr(2, 2)
TYPE(cusparseHandle) :: sparseHandle
INTEGER(KIND = cuda_stream_kind) :: stream
INTEGER :: comm

REAL(4), ALLOCATABLE, DEVICE :: xSelf(:, :), xNeigh(:, :), yTemp(:, :)
REAL(4), POINTER, DEVICE :: csrVal(:)
INTEGER, POINTER, DEVICE :: csrRowPtr(:), csrColIdx(:)
INTEGER :: n, nb(2), nt(2), nr, nc, nnz
INTEGER :: ierr

csrVal => Diag%d_csrVal
csrRowPtr => Diag%d_csrRowPtr
csrColIdx => Diag%d_csrColIdx
nr = Diag%nr
nc = Diag%nc
nnz = Diag%nnz
nb = br(2, :) - br(1, :) + 1
nt = tr(2, :) - tr(1, :) + 1
n = sum(nb)

IF (.NOT. cuCntl%lMulti .OR. cuCntl%lDcpl) THEN

  ierr = cusparseScsrmv(sparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, 1.0_4, Diag%descr,             &
                        csrVal, csrRowPtr, csrColIdx, x, 0.0_4, y)
  RETURN

ENDIF

ALLOCATE(xSelf(n, 2))
ALLOCATE(xNeigh(n, 2))
ALLOCATE(yTemp(n, 2))

CALL cuInitArray(n + n, xNeigh, stream)

ierr = cudaMemcpyAsync(xSelf(1 : nb(black), bottom), x(br(1, black) : br(2, black)), nb(black),                     &
                       cudaMemcpyDeviceToDevice, stream)
ierr = cudaMemcpyAsync(xSelf(nb(black) + 1 : n, bottom), x(br(1, red) : br(2, red)), nb(red),                       &
                       cudaMemcpyDeviceToDevice, stream)
ierr = cudaMemcpyAsync(xSelf(1 : nt(black), top), x(tr(1, black) : tr(2, black)), nt(black),                        &
                       cudaMemcpyDeviceToDevice, stream)
ierr = cudaMemcpyAsync(xSelf(nt(black) + 1 : n, top), x(tr(1, red) : tr(2, red)), nt(red),                          &
                       cudaMemcpyDeviceToDevice, stream)
ierr = cudaStreamSynchronize(stream)

CALL InitMPIComm()
CALL MPIComm(n, n, xSelf(:, bottom), xNeigh(:, bottom), bottom, comm)
CALL MPIComm(n, n, xSelf(:, top), xNeigh(:, top), top, comm)

ierr = cusparseScsrmv(sparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, 1.0_4, Diag%descr,               &
                      csrVal, csrRowPtr, csrColIdx, x, 0.0_4, y)

CALL FinalizeMPIComm()

CALL cuVectorOp('*', n + n, offDiag, xNeigh, yTemp, stream)
CALL cuVectorOp('+', nb(red), yTemp(1 : nb(red), bottom), y(br(1, red) : br(2, red)),                               &
                y(br(1, red) : br(2, red)), stream)
CALL cuVectorOp('+', nb(black), yTemp(nb(red) + 1 : n, bottom), y(br(1, black) : br(2, black)),                     &
                y(br(1, black) : br(2, black)), stream)
CALL cuVectorOp('+', nt(red), yTemp(1 : nt(red), top), y(tr(1, red) : tr(2, red)),                                  &
                y(tr(1, red) : tr(2, red)), stream)
CALL cuVectorOp('+', nt(black), yTemp(nt(red) + 1 : n, top), y(tr(1, black) : tr(2, black)),                        &
                y(tr(1, black) : tr(2, black)), stream)

ierr = cudaStreamSynchronize(stream)

DEALLOCATE(xSelf)
DEALLOCATE(xNeigh)
DEALLOCATE(yTemp)

END SUBROUTINE

SUBROUTINE cuMatMulDouble(Diag, offDiag, x, y, br, tr, sparseHandle, stream, comm)

IMPLICIT NONE

TYPE(CSR_DOUBLE) :: Diag
REAL(8), DEVICE :: offDiag(*), x(*), y(*)
INTEGER :: br(2), tr(2)
TYPE(cusparseHandle) :: sparseHandle
INTEGER(KIND = cuda_stream_kind) :: stream
INTEGER :: comm

REAL(8), ALLOCATABLE, DEVICE :: xNeigh(:), yTemp(:)
REAL(8), POINTER, DEVICE :: csrVal(:)
INTEGER, POINTER, DEVICE :: csrRowPtr(:), csrColIdx(:)
INTEGER :: n, nb, nt, nr, nc, nnz
INTEGER :: ierr

csrVal => Diag%d_csrVal
csrRowPtr => Diag%d_csrRowPtr
csrColIdx => Diag%d_csrColIdx
nr = Diag%nr
nc = Diag%nc
nnz = Diag%nnz
nb = br(2) - br(1) + 1
nt = tr(2) - tr(1) + 1
n = max(nb, nt)

IF (.NOT. cuCntl%lMulti .OR. cuCntl%lDcpl) THEN

  ierr = cusparseDcsrmv(sparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, 1.0_8, Diag%descr,             &
                        csrVal, csrRowPtr, csrColIdx, x, 0.0_8, y)
  RETURN

ENDIF

ALLOCATE(xNeigh(nb + nt))
ALLOCATE(yTemp(nb + nt))

CALL cuInitArray(nb + nt, xNeigh, stream)

ierr = cudaStreamSynchronize(stream)

CALL InitMPIComm()
CALL MPIComm(nb, nb, x(br(1) : br(2)), xNeigh(1 : nb), bottom, comm)
CALL MPIComm(nt, nt, x(tr(1) : tr(2)), xNeigh(nb + 1 : nb + nt), top, comm)

ierr = cusparseDcsrmv(sparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, 1.0_8, Diag%descr,               &
                      csrVal, csrRowPtr, csrColIdx, x, 0.0_8, y)

CALL FinalizeMPIComm()

CALL cuVectorOp('*', nb + nt, offDiag, xNeigh, yTemp, stream)
CALL cuVectorOp('+', nb, yTemp(1 : nb), y(br(1) : br(2)), y(br(1) : br(2)), stream)
CALL cuVectorOp('+', nt, yTemp(nb + 1 : nb + nt), y(tr(1) : tr(2)), y(tr(1) : tr(2)), stream)

ierr = cudaStreamSynchronize(stream)

DEALLOCATE(xNeigh)
DEALLOCATE(yTemp)

END SUBROUTINE

SUBROUTINE cuMatMulDoubleRB(Diag, offDiag, x, y, br, tr, sparseHandle, stream, comm)

IMPLICIT NONE

TYPE(CSR_DOUBLE) :: Diag
REAL(8), DEVICE :: offDiag(*), x(*), y(*)
INTEGER :: br(2, 2), tr(2, 2)
TYPE(cusparseHandle) :: sparseHandle
INTEGER(KIND = cuda_stream_kind) :: stream
INTEGER :: comm

REAL(8), ALLOCATABLE, DEVICE :: xSelf(:, :), xNeigh(:, :), yTemp(:, :)
REAL(8), POINTER, DEVICE :: csrVal(:)
INTEGER, POINTER, DEVICE :: csrRowPtr(:), csrColIdx(:)
INTEGER :: n, nb(2), nt(2), nr, nc, nnz
INTEGER :: ierr

csrVal => Diag%d_csrVal
csrRowPtr => Diag%d_csrRowPtr
csrColIdx => Diag%d_csrColIdx
nr = Diag%nr
nc = Diag%nc
nnz = Diag%nnz
nb = br(2, :) - br(1, :) + 1
nt = tr(2, :) - tr(1, :) + 1
n = sum(nb)

IF (.NOT. cuCntl%lMulti .OR. cuCntl%lDcpl) THEN

  ierr = cusparseDcsrmv(sparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, 1.0_8, Diag%descr,             &
                        csrVal, csrRowPtr, csrColIdx, x, 0.0_8, y)
  RETURN

ENDIF

ALLOCATE(xSelf(n, 2))
ALLOCATE(xNeigh(n, 2))
ALLOCATE(yTemp(n, 2))

CALL cuInitArray(n + n, xNeigh, stream)

ierr = cudaMemcpyAsync(xSelf(1 : nb(black), bottom), x(br(1, black) : br(2, black)), nb(black),                     &
                       cudaMemcpyDeviceToDevice, stream)
ierr = cudaMemcpyAsync(xSelf(nb(black) + 1 : n, bottom), x(br(1, red) : br(2, red)), nb(red),                       &
                       cudaMemcpyDeviceToDevice, stream)
ierr = cudaMemcpyAsync(xSelf(1 : nt(black), top), x(tr(1, black) : tr(2, black)), nt(black),                        &
                       cudaMemcpyDeviceToDevice, stream)
ierr = cudaMemcpyAsync(xSelf(nt(black) + 1 : n, top), x(tr(1, red) : tr(2, red)), nt(red),                          &
                       cudaMemcpyDeviceToDevice, stream)
ierr = cudaStreamSynchronize(stream)

CALL InitMPIComm()
CALL MPIComm(n, n, xSelf(:, bottom), xNeigh(:, bottom), bottom, comm)
CALL MPIComm(n, n, xSelf(:, top), xNeigh(:, top), top, comm)

ierr = cusparseDcsrmv(sparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, 1.0_8, Diag%descr,               &
                      csrVal, csrRowPtr, csrColIdx, x, 0.0_8, y)

CALL FinalizeMPIComm()

CALL cuVectorOp('*', n + n, offDiag, xNeigh, yTemp, stream)
CALL cuVectorOp('+', nb(red), yTemp(1 : nb(red), bottom), y(br(1, red) : br(2, red)),                               &
                y(br(1, red) : br(2, red)), stream)
CALL cuVectorOp('+', nb(black), yTemp(nb(red) + 1 : n, bottom), y(br(1, black) : br(2, black)),                     &
                y(br(1, black) : br(2, black)), stream)
CALL cuVectorOp('+', nt(red), yTemp(1 : nt(red), top), y(tr(1, red) : tr(2, red)),                                  &
                y(tr(1, red) : tr(2, red)), stream)
CALL cuVectorOp('+', nt(black), yTemp(nt(red) + 1 : n, top), y(tr(1, black) : tr(2, black)),                        &
                y(tr(1, black) : tr(2, black)), stream)

ierr = cudaStreamSynchronize(stream)

DEALLOCATE(xSelf)
DEALLOCATE(xNeigh)
DEALLOCATE(yTemp)

END SUBROUTINE

SUBROUTINE cuMatMulMixed(Diag, offDiag, x, y, br, tr, sparseHandle, stream, comm)

IMPLICIT NONE

TYPE(CSR_MIXED) :: Diag
REAL(8), DEVICE :: offDiag(*), x(*), y(*)
INTEGER :: br(2), tr(2)
TYPE(cusparseHandle) :: sparseHandle
INTEGER(KIND = cuda_stream_kind) :: stream
INTEGER :: comm

REAL(8), ALLOCATABLE, DEVICE :: xNeigh(:), yTemp(:)
REAL(8), POINTER, DEVICE :: csrVal(:)
INTEGER, POINTER, DEVICE :: csrRowPtr(:), csrColIdx(:)
INTEGER :: n, nb, nt, nr, nc, nnz
INTEGER :: ierr

csrVal => Diag%d_csrVal8
csrRowPtr => Diag%d_csrRowPtr
csrColIdx => Diag%d_csrColIdx
nr = Diag%nr
nc = Diag%nc
nnz = Diag%nnz
nb = br(2) - br(1) + 1
nt = tr(2) - tr(1) + 1
n = max(nb, nt)

IF (.NOT. cuCntl%lMulti .OR. cuCntl%lDcpl) THEN

  ierr = cusparseDcsrmv(sparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, 1.0_8, Diag%descr,             &
                        csrVal, csrRowPtr, csrColIdx, x, 0.0_8, y)
  RETURN

ENDIF

ALLOCATE(xNeigh(nb + nt))
ALLOCATE(yTemp(nb + nt))

CALL cuInitArray(nb + nt, xNeigh, stream)

ierr = cudaStreamSynchronize(stream)

CALL InitMPIComm()
CALL MPIComm(nb, nb, x(br(1) : br(2)), xNeigh(1 : nb), bottom, comm)
CALL MPIComm(nt, nt, x(tr(1) : tr(2)), xNeigh(nb + 1 : nb + nt), top, comm)

ierr = cusparseDcsrmv(sparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, 1.0_8, Diag%descr,               &
                      csrVal, csrRowPtr, csrColIdx, x, 0.0_8, y)

CALL FinalizeMPIComm()

CALL cuVectorOp('*', nb + nt, offDiag, xNeigh, yTemp, stream)
CALL cuVectorOp('+', nb, yTemp(1 : nb), y(br(1) : br(2)), y(br(1) : br(2)), stream)
CALL cuVectorOp('+', nt, yTemp(nb + 1 : nb + nt), y(tr(1) : tr(2)), y(tr(1) : tr(2)), stream)

ierr = cudaStreamSynchronize(stream)

DEALLOCATE(xNeigh)
DEALLOCATE(yTemp)

END SUBROUTINE

SUBROUTINE cuMatMulMixedRB(Diag, offDiag, x, y, br, tr, sparseHandle, stream, comm)

IMPLICIT NONE

TYPE(CSR_MIXED) :: Diag
REAL(8), DEVICE :: offDiag(*), x(*), y(*)
INTEGER :: br(2, 2), tr(2, 2)
TYPE(cusparseHandle) :: sparseHandle
INTEGER(KIND = cuda_stream_kind) :: stream
INTEGER :: comm

REAL(8), ALLOCATABLE, DEVICE :: xSelf(:, :), xNeigh(:, :), yTemp(:, :)
REAL(8), POINTER, DEVICE :: csrVal(:)
INTEGER, POINTER, DEVICE :: csrRowPtr(:), csrColIdx(:)
INTEGER :: n, nb(2), nt(2), nr, nc, nnz
INTEGER :: ierr

csrVal => Diag%d_csrVal8
csrRowPtr => Diag%d_csrRowPtr
csrColIdx => Diag%d_csrColIdx
nr = Diag%nr
nc = Diag%nc
nnz = Diag%nnz
nb = br(2, :) - br(1, :) + 1
nt = tr(2, :) - tr(1, :) + 1
n = sum(nb)

IF (.NOT. cuCntl%lMulti .OR. cuCntl%lDcpl) THEN

  ierr = cusparseDcsrmv(sparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, 1.0_8, Diag%descr,             &
                        csrVal, csrRowPtr, csrColIdx, x, 0.0_8, y)
  RETURN

ENDIF

ALLOCATE(xSelf(n, 2))
ALLOCATE(xNeigh(n, 2))
ALLOCATE(yTemp(n, 2))

CALL cuInitArray(n + n, xNeigh, stream)

ierr = cudaMemcpyAsync(xSelf(1 : nb(black), bottom), x(br(1, black) : br(2, black)), nb(black),                     &
                       cudaMemcpyDeviceToDevice, stream)
ierr = cudaMemcpyAsync(xSelf(nb(black) + 1 : n, bottom), x(br(1, red) : br(2, red)), nb(red),                       &
                       cudaMemcpyDeviceToDevice, stream)
ierr = cudaMemcpyAsync(xSelf(1 : nt(black), top), x(tr(1, black) : tr(2, black)), nt(black),                        &
                       cudaMemcpyDeviceToDevice, stream)
ierr = cudaMemcpyAsync(xSelf(nt(black) + 1 : n, top), x(tr(1, red) : tr(2, red)), nt(red),                          &
                       cudaMemcpyDeviceToDevice, stream)
ierr = cudaStreamSynchronize(stream)

CALL InitMPIComm()
CALL MPIComm(n, n, xSelf(:, bottom), xNeigh(:, bottom), bottom, comm)
CALL MPIComm(n, n, xSelf(:, top), xNeigh(:, top), top, comm)

ierr = cusparseDcsrmv(sparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, 1.0_8, Diag%descr,               &
                      csrVal, csrRowPtr, csrColIdx, x, 0.0_8, y)

CALL FinalizeMPIComm()

CALL cuVectorOp('*', n + n, offDiag, xNeigh, yTemp, stream)
CALL cuVectorOp('+', nb(red), yTemp(1 : nb(red), bottom), y(br(1, red) : br(2, red)),                               &
                y(br(1, red) : br(2, red)), stream)
CALL cuVectorOp('+', nb(black), yTemp(nb(red) + 1 : n, bottom), y(br(1, black) : br(2, black)),                     &
                y(br(1, black) : br(2, black)), stream)
CALL cuVectorOp('+', nt(red), yTemp(1 : nt(red), top), y(tr(1, red) : tr(2, red)),                                  &
                y(tr(1, red) : tr(2, red)), stream)
CALL cuVectorOp('+', nt(black), yTemp(nt(red) + 1 : n, top), y(tr(1, black) : tr(2, black)),                        &
                y(tr(1, black) : tr(2, black)), stream)

ierr = cudaStreamSynchronize(stream)

DEALLOCATE(xSelf)
DEALLOCATE(xNeigh)
DEALLOCATE(yTemp)

END SUBROUTINE

!--- CUDA Collective Vector Operation Routines

SUBROUTINE normalizeMulti4(x, n, comm, blasHandle, lMPI)

IMPLICIT NONE

INTEGER :: n, comm, ierr
TYPE(cublasHandle) :: blasHandle
REAL(4), DEVICE :: x(*)
REAL(4) :: norm
LOGICAL :: lMPI

norm = normMulti4(x, n, comm, blasHandle, lMPI)
ierr = cublasSscal_v2(blasHandle, n, 1.0_4 / norm, x, 1)

END SUBROUTINE

FUNCTION normMulti4(x, n, comm, blasHandle, lMPI) RESULT(norm)

IMPLICIT NONE

INTEGER :: n, comm
TYPE(cublasHandle) :: blasHandle
REAL(4), DEVICE :: x(*)
REAL(4) :: dot, norm
LOGICAL :: lMPI

dot = dotMulti4(x, x, n, comm, blasHandle, lMPI)
norm = sqrt(dot)

END FUNCTION

FUNCTION dotMulti4(x, y, n, comm, blasHandle, lMPI) RESULT(sum_dot)

IMPLICIT NONE

INTEGER :: n, comm, ierr
TYPE(cublasHandle) :: blasHandle
INTEGER(KIND = cuda_stream_kind) :: stream
REAL(4), DEVICE :: x(*), y(*)
REAL(4) :: dot, sum_dot
LOGICAL :: lMPI

ierr = cublasGetStream(blasHandle, stream)
ierr = cublasSdot_v2(blasHandle, n, x, 1, y, 1, dot)
ierr = cudaStreamSynchronize(stream)

IF (lMPI) THEN
  CALL MPI_ALLREDUCE(dot, sum_dot, 1, MPI_FLOAT, MPI_SUM, comm, ierr)
ELSE
  sum_dot = dot
ENDIF

END FUNCTION

FUNCTION asumMulti4(x, n, comm, blasHandle, lMPI) RESULT(sum_asum)
IMPLICIT NONE
INTEGER :: n, comm, ierr
TYPE(cublasHandle) :: blasHandle
INTEGER(KIND=cuda_stream_kind) :: stream
REAL(4), DEVICE :: x(*)
REAL(4) :: asum, sum_asum
LOGICAL :: lMPI

ierr = cublasGetStream(blasHandle, stream)
ierr = cublasSasum_v2(blasHandle, n, x, 1, asum)
ierr = cudaStreamSynchronize(stream)

IF (lMPI) THEN
  CALL MPI_ALLREDUCE(asum, sum_asum, 1, MPI_FLOAT, MPI_SUM, comm, ierr)
ELSE
  sum_asum = asum
END IF

END FUNCTION

SUBROUTINE normalizeMulti8(x, n, comm, blasHandle, lMPI)

IMPLICIT NONE

INTEGER :: n, comm, ierr
TYPE(cublasHandle) :: blasHandle
REAL(8), DEVICE :: x(*)
REAL(8) :: norm
LOGICAL :: lMPI

norm = normMulti8(x, n, comm, blasHandle, lMPI)
ierr = cublasDscal_v2(blasHandle, n, 1.0_8 / norm, x, 1)

END SUBROUTINE

FUNCTION normMulti8(x, n, comm, blasHandle, lMPI) RESULT(norm)

IMPLICIT NONE

INTEGER :: n, comm
TYPE(cublasHandle) :: blasHandle
REAL(8), DEVICE :: x(*)
REAL(8) :: dot, norm
LOGICAL :: lMPI

dot = dotMulti8(x, x, n, comm, blasHandle, lMPI)
norm = sqrt(dot)

END FUNCTION

FUNCTION dotMulti8(x, y, n, comm, blasHandle, lMPI) RESULT(sum_dot)

IMPLICIT NONE

INTEGER :: n, comm, ierr
TYPE(cublasHandle) :: blasHandle
INTEGER(KIND = cuda_stream_kind) :: stream
REAL(8), DEVICE :: x(*), y(*)
REAL(8) :: dot, sum_dot
LOGICAL :: lMPI

ierr = cublasGetStream(blasHandle, stream)
ierr = cublasDdot_v2(blasHandle, n, x, 1, y, 1, dot)
ierr = cudaStreamSynchronize(stream)

IF (lMPI) THEN
  CALL MPI_ALLREDUCE(dot, sum_dot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
ELSE
  sum_dot = dot
ENDIF

END FUNCTION

FUNCTION asumMulti8(x, n, comm, blasHandle, lMPI) RESULT(sum_asum)
IMPLICIT NONE
INTEGER :: n, comm, ierr
TYPE(cublasHandle) :: blasHandle
INTEGER(KIND=cuda_stream_kind) :: stream
REAL(8), DEVICE :: x(*)
REAL(8) :: asum, sum_asum
LOGICAL :: lMPI

ierr = cublasGetStream(blasHandle, stream)
ierr = cublasDasum_v2(blasHandle, n, x, 1, asum)
ierr = cudaStreamSynchronize(stream)

IF (lMPI) THEN
  CALL MPI_ALLREDUCE(asum, sum_asum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
ELSE
  sum_asum = asum
END IF

END FUNCTION

!--- CUDA Element-wise Vector Operation Routines

SUBROUTINE cuVectorOp4(op, n, x, y, z, stream)

IMPLICIT NONE

CHARACTER(*) :: op
INTEGER :: n
REAL(4), DEVICE :: x(*), y(*), z(*)
INTEGER(KIND = cuda_stream_kind) :: stream

TYPE(dim3) :: Blocks, Threads

Threads = dim3(1024, 1, 1)
Blocks = dim3(n / Threads%x + 1, 1, 1)

SELECT CASE (op)

CASE ('+')
  CALL VectorAdd4 <<< Blocks, Threads, 0, stream >>> (n, x, y, z)

CASE ('-')
  CALL VectorSub4 <<< Blocks, Threads, 0, stream >>> (n, x, y, z)

CASE ('*')
  CALL VectorMul4 <<< Blocks, Threads, 0, stream >>> (n, x, y, z)

CASE ('/')
  CALL VectorDiv4 <<< Blocks, Threads, 0, stream >>> (n, x, y, z)

CASE ('*+')
  CALL VectorMulAdd4 <<< Blocks, Threads, 0, stream >>> (n, x, y, z)

CASE ('-+')
  CALL VectorMulSub4 <<< Blocks, Threads, 0, stream >>> (n, x, y, z)

END SELECT

END SUBROUTINE

SUBROUTINE cuVectorOp8(op, n, x, y, z, stream)

IMPLICIT NONE

CHARACTER(*) :: op
INTEGER :: n
REAL(8), DEVICE :: x(*), y(*), z(*)
INTEGER(KIND = cuda_stream_kind) :: stream

TYPE(dim3) :: Blocks, Threads

Threads = dim3(1024, 1, 1)
Blocks = dim3(n / Threads%x + 1, 1, 1)

SELECT CASE (Op)

CASE ('+')
  CALL VectorAdd8 <<< Blocks, Threads, 0, stream >>> (n, x, y, z)

CASE ('-')
  CALL VectorSub8 <<< Blocks, Threads, 0, stream >>> (n, x, y, z)

CASE ('*')
  CALL VectorMul8 <<< Blocks, Threads, 0, stream >>> (n, x, y, z)

CASE ('/')
  CALL VectorDiv8 <<< Blocks, Threads, 0, stream >>> (n, x, y, z)

CASE ('*+')
  CALL VectorMulAdd8 <<< Blocks, Threads, 0, stream >>> (n, x, y, z)

CASE ('*-')
  CALL VectorMulSub8 <<< Blocks, Threads, 0, stream >>> (n, x, y, z)

END SELECT

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE VectorAdd4(n, x, y, z)
USE CUDAFOR
IMPLICIT NONE

INTEGER, VALUE :: n
REAL(4), DEVICE :: x(n), y(n), z(n)

INTEGER :: idx

idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
IF (idx .GT. n) RETURN

z(idx) = x(idx) + y(idx)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE VectorSub4(n, x, y, z)
USE CUDAFOR
IMPLICIT NONE

INTEGER, VALUE :: n
REAL(4), DEVICE :: x(n), y(n), z(n)

INTEGER :: idx

idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
IF (idx .GT. n) RETURN

z(idx) = x(idx) - y(idx)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE VectorMul4(n, x, y, z)
USE CUDAFOR
IMPLICIT NONE

INTEGER, VALUE :: n
REAL(4), DEVICE :: x(n), y(n), z(n)

INTEGER :: idx

idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
IF (idx .GT. n) RETURN

z(idx) = x(idx) * y(idx)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE VectorDiv4(n, x, y, z)
USE CUDAFOR
IMPLICIT NONE

INTEGER, VALUE :: n
REAL(4), DEVICE :: x(n), y(n), z(n)

INTEGER :: idx

idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
IF (idx .GT. n) RETURN

z(idx) = x(idx) / y(idx)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE VectorMulAdd4(n, x, y, z)
USE CUDAFOR
IMPLICIT NONE

INTEGER, VALUE :: n
REAL(4), DEVICE :: x(n), y(n), z(n)

INTEGER :: idx

idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
IF (idx .GT. n) RETURN

z(idx) = z(idx) + x(idx) * y(idx)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE VectorMulSub4(n, x, y, z)
USE CUDAFOR
IMPLICIT NONE

INTEGER, VALUE :: n
REAL(4), DEVICE :: x(n), y(n), z(n)

INTEGER :: idx

idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
IF (idx .GT. n) RETURN

z(idx) = z(idx) - x(idx) * y(idx)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE VectorAdd8(n, x, y, z)
USE CUDAFOR
IMPLICIT NONE

INTEGER, VALUE :: n
REAL(8), DEVICE :: x(n), y(n), z(n)

INTEGER :: idx

idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
IF (idx .GT. n) RETURN

z(idx) = x(idx) + y(idx)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE VectorSub8(n, x, y, z)
USE CUDAFOR
IMPLICIT NONE

INTEGER, VALUE :: n
REAL(8), DEVICE :: x(n), y(n), z(n)

INTEGER :: idx

idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
IF (idx .GT. n) RETURN

z(idx) = x(idx) - y(idx)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE VectorMul8(n, x, y, z)
USE CUDAFOR
IMPLICIT NONE

INTEGER, VALUE :: n
REAL(8), DEVICE :: x(n), y(n), z(n)

INTEGER :: idx

idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
IF (idx .GT. n) RETURN

z(idx) = x(idx) * y(idx)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE VectorDiv8(n, x, y, z)
USE CUDAFOR
IMPLICIT NONE

INTEGER, VALUE :: n
REAL(8), DEVICE :: x(n), y(n), z(n)

INTEGER :: idx

idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
IF (idx .GT. n) RETURN

z(idx) = x(idx) / y(idx)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE VectorMulAdd8(n, x, y, z)
USE CUDAFOR
IMPLICIT NONE

INTEGER, VALUE :: n
REAL(8), DEVICE :: x(n), y(n), z(n)

INTEGER :: idx

idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
IF (idx .GT. n) RETURN

z(idx) = z(idx) + x(idx) * y(idx)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE VectorMulSub8(n, x, y, z)
USE CUDAFOR
IMPLICIT NONE

INTEGER, VALUE :: n
REAL(8), DEVICE :: x(n), y(n), z(n)

INTEGER :: idx

idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
IF (idx .GT. n) RETURN

z(idx) = z(idx) - x(idx) * y(idx)

END SUBROUTINE

SUBROUTINE cuHtoDCopy8(n, x, y)

IMPLICIT NONE

INTEGER :: n
REAL(8) :: x(*)
REAL(8), DEVICE :: y(*)

INTEGER :: ierr

ierr = cudaMemcpy(y, x, n, cudaMemcpyHostToDevice)

END SUBROUTINE

SUBROUTINE cuDtoHCopy8(n, x, y)

IMPLICIT NONE

INTEGER :: n
REAL(8), DEVICE :: x(*)
REAL(8) :: y(*)

INTEGER :: ierr

ierr = cudaMemcpy(y, x, n, cudaMemcpyDeviceToHost)

END SUBROUTINE

SUBROUTINE cuHtoDTypecastFrom8(n, x, y)

IMPLICIT NONE

INTEGER :: n
REAL(8) :: x(*)
REAL(4), DEVICE :: y(*)

REAL(4), ALLOCATABLE :: x4(:)
INTEGER :: i, ierr

ALLOCATE(x4(n))

!$OMP PARALLEL DO
DO i = 1, n
  x4(i) = x(i)
ENDDO
!$OMP END PARALLEL DO

ierr = cudaMemcpy(y, x4, n, cudaMemcpyHostToDevice)

DEALLOCATE(x4)

END SUBROUTINE

SUBROUTINE cuDtoHTypecastFrom4(n, x, y)

IMPLICIT NONE

INTEGER :: n
REAL(4), DEVICE :: x(*)
REAL(8) :: y(*)

REAL(4), ALLOCATABLE :: y4(:)
INTEGER :: i, ierr

ALLOCATE(y4(n))

ierr = cudaMemcpy(y4, x, n, cudaMemcpyDeviceToHost)

!$OMP PARALLEL DO
DO i = 1, n
  y(i) = y4(i)
ENDDO
!$OMP END PARALLEL DO

DEALLOCATE(y4)

END SUBROUTINE

SUBROUTINE cuTypecastFrom4(n, x, y, stream)

IMPLICIT NONE

INTEGER :: n
REAL(4), DEVICE :: x(*)
REAL(8), DEVICE :: y(*)
INTEGER(KIND = cuda_stream_kind) :: stream

TYPE(dim3) :: Blocks, Threads

Threads = dim3(1024, 1, 1)
Blocks = dim3(n / Threads%x + 1, 1, 1)

CALL TypecastCopyFrom4 <<< Blocks, Threads, 0, stream >>> (n, x, y)

END SUBROUTINE

SUBROUTINE cuTypecastFrom8(n, x, y, stream)

IMPLICIT NONE

INTEGER :: n
REAL(8), DEVICE :: x(*)
REAL(4), DEVICE :: y(*)
INTEGER(KIND = cuda_stream_kind) :: stream

TYPE(dim3) :: Blocks, Threads

Threads = dim3(1024, 1, 1)
Blocks = dim3(n / Threads%x + 1, 1, 1)

CALL TypecastCopyFrom8 <<< Blocks, Threads, 0, stream >>> (n, x, y)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE TypecastCopyFrom4(n, x, y)
USE CUDAFOR
IMPLICIT NONE

INTEGER, VALUE :: n
REAL(4), DEVICE :: x(n)
REAL(8), DEVICE :: y(n)

INTEGER :: idx

idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
IF (idx .GT. n) RETURN

y(idx) = x(idx)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE TypecastCopyFrom8(n, x, y)
USE CUDAFOR
IMPLICIT NONE

INTEGER, VALUE :: n
REAL(8), DEVICE :: x(n)
REAL(4), DEVICE :: y(n)

INTEGER :: idx

idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
IF (idx .GT. n) RETURN

y(idx) = x(idx)

END SUBROUTINE

SUBROUTINE cuInitArray4(n, Arr, stream)

IMPLICIT NONE

REAL(4), DEVICE :: Arr(*)
INTEGER :: n
INTEGER(KIND = cuda_stream_kind) :: stream

TYPE(dim3) :: Blocks, Threads

Threads = dim3(1024, 1, 1)
Blocks = dim3(n / Threads%x + 1, 1, 1)

CALL ZeroArray4 <<< Blocks, Threads, 0, stream >>> (n, Arr)

END SUBROUTINE

SUBROUTINE cuInitArray8(n, Arr, stream)

IMPLICIT NONE

REAL(8), DEVICE :: Arr(*)
INTEGER :: n
INTEGER(KIND = cuda_stream_kind) :: stream

TYPE(dim3) :: Blocks, Threads

Threads = dim3(1024, 1, 1)
Blocks = dim3(n / Threads%x + 1, 1, 1)

CALL ZeroArray8 <<< Blocks, Threads, 0, stream >>> (n, Arr)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE ZeroArray4(n, Arr)
USE CUDAFOR
IMPLICIT NONE

INTEGER, VALUE :: n
REAL(4), DEVICE :: Arr(n)

INTEGER :: idx

idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
IF (idx .GT. n) RETURN

Arr(idx) = 0.0

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE ZeroArray8(n, Arr)
USE CUDAFOR
IMPLICIT NONE

INTEGER, VALUE :: n
REAL(8), DEVICE :: Arr(n)

INTEGER :: idx

idx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
IF (idx .GT. n) RETURN

Arr(idx) = 0.0

END SUBROUTINE


ATTRIBUTES(GLOBAL) SUBROUTINE cuComputeTranCMFDSrc(trSrc, resSrc, precSrc, &
                         TranPhi, VolInvVel, expo_alpha, expo, delt, nzCMFD)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE
REAL(8), DEVICE :: trSrc(:,:)
REAL(8), DEVICE, INTENT(IN) :: resSrc(:,:), precSrc(:), TranPhi(:,:)
REAL(8), DEVICE, INTENT(IN) :: VolInvVel(:,:), expo_alpha(:,:), expo(:,:)
REAL(8), VALUE :: delt
INTEGER, VALUE :: nzCMFD

REAL(8) :: thetah, prevSrc
INTEGER :: ncel
INTEGER :: ig, icel

ncel = nxyc * nzCMFD

ig = threadIdx%x
icel = threadIdx%y + (blockIdx%x - 1) * blockDim%y
IF(icel .GT. ncel) RETURN

thetah = 1./ theta - 1.

trSrc(ig, icel) = chid(ig) * PrecSrc(icel)
prevSrc = Thetah * (ResSrc(ig, icel) - expo_alpha(ig, icel) * VolInvVel(ig, icel) * TranPhi(ig, icel))
prevSrc = prevSrc + VolInvVel(ig, icel) * TranPhi(ig, icel) / (delt * theta)
prevSrc = prevSrc * Expo(ig, icel)
trSrc(ig, icel) = trSrc(ig, icel) + prevSrc

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuComputeTranCMFDSrc11(trSrc, resSrc, precSrc, &
                         TranPhi, VolInvVel, expo_alpha, expo, delt, nzCMFD)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE
REAL(8), DEVICE :: trSrc(:,:)
REAL(8), DEVICE, INTENT(IN) :: resSrc(:,:), precSrc(:), TranPhi(:,:)
REAL(8), DEVICE, INTENT(IN) :: VolInvVel(:,:), expo_alpha(:,:), expo(:,:)
REAL(8), VALUE :: delt
INTEGER, VALUE :: nzCMFD

REAL(8) :: thetah, prevSrc, volphid
INTEGER :: ncel
INTEGER :: ig, icel

IF(threadIdx%x .EQ. 1 .AND. threadIdx%y .EQ. 1) Print*, 'compute transrc'
ncel = nxyc * nzCMFD

ig = threadIdx%x
icel = threadIdx%y + (blockIdx%x - 1) * blockDim%y
IF(icel .GT. ncel) RETURN
thetah = 1./ theta - 1.

volphid = VolInvVel(ig,icel) * TranPhi(ig,icel)
prevSrc = expo_alpha(ig, icel) * volphid
prevSrc = ResSrc(ig,icel) - prevSrc
prevSrc = thetah* prevSrc
!prevSrc = prevSrc + volphid / (delt*theta)
!prevSrc = prevSrc * Expo(ig,icel)
trSrc(ig,icel) = prevSrc
!trSrc(ig,icel) = chid(ig) * precSrc(icel)
!trSrc(ig,icel) = trSrc(ig,icel) + prevSrc

IF(ig .EQ. 1 .AND. icel .EQ. 1) PRINT*, Trsrc(ig, icel)

END SUBROUTINE
ATTRIBUTES(GLOBAL) SUBROUTINE cuComputeTranCMFDSrc22(trSrc, resSrc, &
                         TranPhi, VolInvVel, ex_alpha, delt, nzCMFD)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE
REAL(8), DEVICE :: trSrc(:,:)
REAL(8), DEVICE, INTENT(IN) :: resSrc(:,:), TranPhi(:,:)
REAL(8), DEVICE, INTENT(IN) :: VolInvVel(:,:), ex_alpha(:,:)
REAL(8), VALUE :: delt
INTEGER, VALUE :: nzCMFD

REAL(8) :: thetah
INTEGER :: ncel
INTEGER :: ig, icel

IF(threadIdx%x .EQ. 1 .AND. threadIdx%y .EQ. 1) Print*, 'compute transrc'
ncel = nxyc * nzCMFD

ig = threadIdx%x
icel = threadIdx%y + (blockIdx%x - 1) * blockDim%y
IF(icel .GT. ncel) RETURN
thetah = 1./ theta - 1.

trSrc(ig,icel) = ResSrc(ig,icel)
TrSrc(ig,icel) = TrSrc(ig,icel) * VolInvVel(ig,icel) * TranPhi(ig,icel)

IF(ig .EQ. 1 .AND. icel .EQ. 1) PRINT*, Trsrc(ig, icel)

END SUBROUTINE

END MODULE

#endif
