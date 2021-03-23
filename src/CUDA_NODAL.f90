#include <CUDADEFINES.h>
    
#ifdef __PGI

MODULE CUDA_NODAL

USE CUDA_MASTER
USE CUDA_UTIL
IMPLICIT NONE

REAL(GPU_NODAL_PRECISION), PARAMETER, PRIVATE :: zero = 0.0, one = 1.0

PRIVATE
PUBLIC :: cuAllocNodal, cuNodalDriver


CONTAINS

!--- Public Routines ------------------------------------------------------------------------------

SUBROUTINE cuAllocNodal(cuAxial, cuDevice)

IMPLICIT NONE

TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice

INTEGER :: ng, nxy, nzCMFD

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD

ALLOCATE(cuAxial%phiCoeff(ng, nxy, nzCMFD, 0 : 4)); cuAxial%phiCoeff = 0.0
ALLOCATE(cuAxial%srcCoeff(ng, nxy, nzCMFD, 0 : 4)); cuAxial%srcCoeff = 0.0
ALLOCATE(cuAxial%lkgCoeff(ng, nxy, nzCMFD, 0 : 2)); cuAxial%lkgCoeff = 0.0
ALLOCATE(cuAxial%psiCoeff(nxy, nzCMFD, 0 : 4)); cuAxial%psiCoeff = 0.0
ALLOCATE(cuAxial%D(ng, nxy, nzCMFD))
ALLOCATE(cuAxial%sigR(ng, nxy, nzCMFD))

END SUBROUTINE

SUBROUTINE cuNodalDriver(cuCMFD, cuAxial, cuDevice, eigv)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice
REAL :: eigv

TYPE(dim3) :: Blocks, Threads
INTEGER :: ng, nxy, nzCMFD
INTEGER :: iter, order
INTEGER(KIND = cuda_stream_kind) :: myStream

ng = cuCMFD%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
myStream = cuDevice%myStream

CALL cuSetSourceOperator(cuCMFD, cuAxial, cuDevice)
CALL cuSetCrossSection(cuCMFD, cuAxial, cuDevice)
CALL cuSetAverageFlux(cuCMFD, cuAxial, cuDevice)
CALL cuSetPsiCoeff(cuAxial, cuDevice, 0)
CALL cuSetTransverseLeakage(cuCMFD, cuAxial, cuDevice)

Threads = dim3(cuDevice%cuThreadPerBlock, 1, 1)
Blocks = dim3(ng * nxy / Threads%x + 1, nzCMFD, 1)

DO iter = 1, 50
  DO order = 0, 4
    IF (order .NE. 0) CALL cuSetPsiCoeff(cuAxial, cuDevice, order)
    CALL cuSetSourceCoeff(cuAxial, cuDevice, eigv, order)
  ENDDO
  !$ACC HOST_DATA USE_DEVICE(cuGeometry, cuDevice)
  IF (cuCntl%lSENM) THEN
    CALL cuSENM <<< Blocks, Threads, 0, myStream >>>                                                                &
                (cuGeometry, cuDevice, cuAxial%phiCoeff, cuAxial%srcCoeff, cuAxial%Jout,                            &
                 cuAxial%D, cuAxial%sigR)
  ELSE
    CALL cuNEM <<< Blocks, Threads, 0, myStream >>>                                                                 &
               (cuGeometry, cuDevice, cuAxial%phiCoeff, cuAxial%srcCoeff, cuAxial%Jout,                             &
                cuAxial%D, cuAxial%sigR)
  ENDIF
  !$ACC END HOST_DATA          
  CALL cuGetIncomingCurrent(cuAxial, cuDevice)
ENDDO

END SUBROUTINE

!--- Private Routines -----------------------------------------------------------------------------

SUBROUTINE cuSetSourceOperator(cuCMFD, cuAxial, cuDevice)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice

INTEGER, POINTER :: pinMap(:), planeMap(:)
INTEGER :: idx, ir, ic, ig, igs, igb, ige, ipin, ipin_map, iz, izf
INTEGER :: ng, nxy, nzCMFD
INTEGER :: myzbf, myzef
REAL(GPU_NODAL_PRECISION) :: val

ng = cuCMFD%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
pinMap => cuGeometry%pinMap
planeMap => cuGeometry%planeMap

CALL createCsr(cuAxial%S, ng * ng * nxy * nzCMFD, ng * nxy * nzCMFD, ng * nxy * nzCMFD)
DO izf = myzbf, myzef
  iz = planeMap(izf)
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    DO ig = 1, ng
      ir = ig + (ipin - 1) * ng + (izf - myzbf) * ng * nxy
      igb = cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%ib
      ige = cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%ie
      DO igs = igb, ige
        ic = ir + (igs - ig)
        val = cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%from(igs)
        CALL pushCsr(cuAxial%S, val, ir, ic)
      ENDDO
    ENDDO
  ENDDO
ENDDO
CALL finalizeCsr(cuAxial%S, .TRUE.)

IF (.NOT. cuDevice%lFuel) RETURN

CALL createCsr(cuAxial%F, ng * nxy * nzCMFD, nxy * nzCMFD, ng * nxy * nzCMFD)
DO izf = myzbf, myzef
  iz = planeMap(izf)
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    ir = ipin + (izf - myzbf) * nxy
    DO ig = 1, ng
      ic = ig + (ipin - 1) * ng + (izf - myzbf) * ng * nxy
      val = cuCMFD%PinXS(ipin_map, iz)%XSnf(ig)
      CALL pushCsr(cuAxial%F, val, ir, ic)
    ENDDO
  ENDDO
ENDDO
CALL finalizeCsr(cuAxial%F, .TRUE.)

CALL createCsr(cuAxial%Chi, ng * nxy * nzCMFD, ng * nxy * nzCMFD, nxy * nzCMFD)
DO izf = myzbf, myzef
  iz = planeMap(izf)
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    ic = ipin + (izf - myzbf) * nxy
    DO ig = 1, ng
      ir = ig + (ipin - 1) * ng + (izf - myzbf) * ng * nxy
      val = cuCMFD%PinXS(ipin_map, iz)%Chi(ig)
      CALL pushCsr(cuAxial%Chi, val, ir, ic)
    ENDDO
  ENDDO
ENDDO
CALL finalizeCsr(cuAxial%Chi, .TRUE.)

END SUBROUTINE

SUBROUTINE cuSetCrossSection(cuCMFD, cuAxial, cuDevice)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice

INTEGER :: ig, ipin, ipin_map, iz, izf, ierr
INTEGER :: n, ng, nxy, nzCMFD
INTEGER :: myzbf, myzef
INTEGER, POINTER :: pinMap(:), planeMap(:)
REAL(GPU_NODAL_PRECISION), ALLOCATABLE, PINNED :: sigR(:, :, :), D(:, :, :)

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
pinMap => cuGeometry%pinMap
planeMap => cuGeometry%planeMap
n = ng * nxy * nzCMFD

ALLOCATE(sigR(ng, nxy, myzbf : myzef))
ALLOCATE(D(ng, nxy, myzbf : myzef))

DO izf = myzbf, myzef
  iz = planeMap(izf)
  !$OMP PARALLEL PRIVATE(ipin_map)
  !$OMP DO SCHEDULE(GUIDED)
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    DO ig = 1, ng
      sigR(ig, ipin, izf) = cuCMFD%PinXS(ipin_map, iz)%XSr(ig)
      D(ig, ipin, izf) = cuCMFD%PinXS(ipin_map, iz)%XSD(ig)
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDDO

ierr = cudaMemcpyAsync(cuAxial%sigR, sigR, n, cudaMemcpyHostToDevice, cuDevice%myStream)
ierr = cudaMemcpyAsync(cuAxial%D, D, n, cudaMemcpyHostToDevice, cuDevice%myStream) 
ierr = cudaStreamSynchronize(cuDevice%myStream)

DEALLOCATE(sigR)
DEALLOCATE(D)

END SUBROUTINE

SUBROUTINE cuSetAverageFlux(cuCMFD, cuAxial, cuDevice)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice

INTEGER :: n, ng, nxy, nzCMFD
INTEGER :: myzbf, myzef
INTEGER :: ig, ipin, iz, ierr
REAL(GPU_NODAL_PRECISION), ALLOCATABLE, PINNED :: phis(:, :, :)

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
n = ng * nxy * nzCMFD

ALLOCATE(phis(ng, nxy, myzbf : myzef))

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO iz = myzbf, myzef
  DO ipin = 1, nxy
    DO ig = 1, ng
      phis(ig, ipin, iz) = cuCMFD%h_phis8(ig, ipin, iz)
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

ierr = cudaMemcpyAsync(cuAxial%phiCoeff(:, :, :, 0), phis, n, cudaMemcpyHostToDevice, cuDevice%myStream)
ierr = cudaStreamSynchronize(cuDevice%myStream)

DEALLOCATE(phis)

END SUBROUTINE

SUBROUTINE cuSetPsiCoeff(cuAxial, cuDevice, order)

IMPLICIT NONE

TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice
INTEGER :: order

REAL(GPU_NODAL_PRECISION), POINTER, DEVICE :: csrVal(:)
INTEGER, POINTER, DEVICE :: csrRowPtr(:), csrColIdx(:)
INTEGER :: nr, nc, nnz
INTEGER :: ierr

IF (.NOT. cuDevice%lFuel) RETURN

csrVal => cuAxial%F%d_csrVal
csrRowPtr => cuAxial%F%d_csrRowPtr
csrColIdx => cuAxial%F%d_csrColIdx
nr = cuAxial%F%nr; nc = cuAxial%F%nc; nnz = cuAxial%F%nnz

ierr = cusparseCsrmv(cuDevice%mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, one,                   &
                     cuAxial%F%descr, csrVal, csrRowPtr, csrColIdx, cuAxial%phiCoeff(:, :, :, order),               &
                     zero, cuAxial%psiCoeff(:, :, order))
                
END SUBROUTINE

SUBROUTINE cuSetTransverseLeakage(cuCMFD, cuAxial, cuDevice)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice

TYPE(superPin_Type), POINTER :: Pin(:)
INTEGER :: n, ng, nxy, nzCMFD
INTEGER :: myzbf, myzef
INTEGER :: ig, ibd, ipin, ipin_map, ineighpin, iz, izf, ierr
INTEGER, POINTER :: pinMap(:), pinMapRev(:), planeMap(:)
REAL :: myphi, neighphi, lkg
REAL :: Dtil, Dhat
REAL(GPU_NODAL_PRECISION), ALLOCATABLE, PINNED :: lkgCoeff(:, :, :, :)

Pin => cuGeometry%superPin
ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
pinMap => cuGeometry%pinMap
pinMapRev => cuGeometry%pinMapRev
planeMap => cuGeometry%planeMap

ALLOCATE(lkgCoeff(ng, nxy, myzbf - 1 : myzef + 1, 0 : 2))

DO izf = myzbf, myzef
  iz = planeMap(izf)
  !$OMP PARALLEL PRIVATE(ipin_map, ineighpin, myphi, neighphi, lkg, Dtil, Dhat)
  !$OMP DO SCHEDULE(GUIDED)
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    DO ig = 1, ng
      myphi = cuCMFD%h_phis8(ig, ipin, izf)
      lkg = 0.0
      DO ibd = 1, 4
        ineighpin = Pin(ipin_map)%NeighIdx(ibd)
        ineighpin = pinMapRev(ineighpin)
        IF (ineighpin .EQ. Void) THEN
          neighphi = 0.0
        ELSEIF (ineighpin .EQ. Reflective) THEN
          neighphi = myphi
        ELSE
          neighphi = cuCMFD%h_phis8(ig, ineighpin, izf)
        ENDIF
        Dtil = cuCMFD%PinXS(ipin_map, iz)%Dtil(ibd, ig)
        Dhat = cuCMFD%PinXS(ipin_map, iz)%Dhat(ibd, ig)
        lkg = lkg - Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
      ENDDO
      lkgCoeff(ig, ipin, izf, 0) = lkg * cuGeometry%hzfm(izf) / cuGeometry%PinVolFm(ipin, izf)
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDDO

CALL cuGetNeighbor(ng * nxy, lkgCoeff(:, :, myzbf, 0), lkgCoeff(:, :, myzbf - 1, 0), MPI_CUDA_COMM, bottom)
CALL cuGetNeighbor(ng * nxy, lkgCoeff(:, :, myzef, 0), lkgCoeff(:, :, myzef + 1, 0), MPI_CUDA_COMM, top)

IF (cuDevice%lBottom) THEN
  IF (cuGeometry%AxBC(bottom) .EQ. Void) lkgCoeff(:, :, myzbf - 1, 0) = 0.0
  IF (cuGeometry%AxBC(bottom) .EQ. Reflective) lkgCoeff(:, :, myzbf - 1, 0) = lkgCoeff(:, :, myzbf, 0)
ENDIF

IF (cuDevice%lTop) THEN
  IF (cuGeometry%AxBC(top) .EQ. Void) lkgCoeff(:, :, myzef + 1, 0) = 0.0
  IF (cuGeometry%AxBC(top) .EQ. Reflective) lkgCoeff(:, :, myzef + 1, 0) = lkgCoeff(:, :, myzef, 0)
ENDIF

IF (cuCntl%lSENM) THEN
  CALL cuSetSENMLeakageCoeff(cuDevice, lkgCoeff)
ELSE
  CALL cuSetNEMLeakageCoeff(cuDevice, lkgCoeff)
ENDIF

n = ng * nxy * nzCMFD
ierr = cudaMemcpyAsync(cuAxial%lkgCoeff(:, :, :, 0), lkgCoeff(:, :, myzbf : myzef, 0), n,                           &
                       cudaMemcpyHostToDevice, cuDevice%myStream)
ierr = cudaMemcpyAsync(cuAxial%lkgCoeff(:, :, :, 1), lkgCoeff(:, :, myzbf : myzef, 1), n,                           &
                       cudaMemcpyHostToDevice, cuDevice%myStream)
ierr = cudaMemcpyAsync(cuAxial%lkgCoeff(:, :, :, 2), lkgCoeff(:, :, myzbf : myzef, 2), n,                           &
                       cudaMemcpyHostToDevice, cuDevice%myStream)
ierr = cudaStreamSynchronize(cuDevice%myStream)

DEALLOCATE(lkgCoeff)

END SUBROUTINE

SUBROUTINE cuSetNEMLeakageCoeff(cuDevice, L)

IMPLICIT NONE

TYPE(cuDevice_Type) :: cuDevice
REAL(GPU_NODAL_PRECISION), TARGET :: L(:, :, cuDevice%myzbf - 1 :, 0 :)
REAL(GPU_NODAL_PRECISION), POINTER :: L0(:, :), L1(:, :), L2(:, :)

REAL :: h0, h2, n1(2), n2(2), d1, d2
REAL, POINTER :: hzfm(:)
INTEGER :: ng, nxy
INTEGER :: myzbf, myzef
INTEGER :: ig, ipin, iz

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
hzfm => cuGeometry%hzfm

DO iz = myzbf, myzef
  h0 = hzfm(iz - 1) / hzfm(iz); h2 = hzfm(iz + 1) / hzfm(iz)
  n1(1) = 1.0 + 3.0 * h2 + 2.0 * h2 ** 2; n2(1) = 1.0 + h2
  n1(2) = 1.0 + 3.0 * h0 + 2.0 * h0 ** 2; n2(2) = 1.0 + h0
  d1 = 2.0 * (1.0 + h0) * (1.0 + h2) * (1.0 + h0 + h2)
  d2 = 2.0 * (1.0 + h2) * (1.0 + 2.0 * h0 + h0 ** 2 + h2 + h0 * h2)
  L0 => L(:, :, iz - 1, 0); L1 => L(:, :, iz, 0); L2 => L(:, :, iz + 1, 0)
  !$OMP PARALLEL
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
  DO ipin = 1, nxy
    DO ig = 1, ng
      L(ig, ipin, iz, 1) = - (L0(ig, ipin) * n1(1) + L1(ig, ipin) * (n1(2) - n1(1)) - L2(ig, ipin) * n1(2)) / d1
      L(ig, ipin, iz, 2) = - (L0(ig, ipin) * n2(1) - L1(ig, ipin) * (n2(2) + n2(1)) + L2(ig, ipin) * n2(2)) / d2
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDDO

END SUBROUTINE

SUBROUTINE cuSetSENMLeakageCoeff(cuDevice, L)

IMPLICIT NONE

TYPE(cuDevice_Type) :: cuDevice
REAL(GPU_NODAL_PRECISION), TARGET :: L(:, :, cuDevice%myzbf - 1 :, 0 :)
REAL(GPU_NODAL_PRECISION), POINTER :: L0(:, :), L1(:, :), L2(:, :)

REAL :: h0, h2, n1(2), n2(2), d1, d2
REAL, POINTER :: hzfm(:)
INTEGER :: ng, nxy
INTEGER :: myzbf, myzef
INTEGER :: ig, ipin, iz

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
hzfm => cuGeometry%hzfm

DO iz = myzbf, myzef
  h0 = hzfm(iz - 1) / hzfm(iz) * 2.0; h2 = hzfm(iz + 1) / hzfm(iz) * 2.0
  n1(1) = 2.0 + 3.0 * h2 + h2 ** 2; n2(1) = 2.0 + h2
  n1(2) = 2.0 + 3.0 * h0 + h0 ** 2; n2(2) = 2.0 + h0
  d1 = (2.0 + h0) * (2.0 + h2) * (2.0 + h0 + h2)
  d2 = (2.0 + h2) * (4.0 + 4.0 * h0 + h0 ** 2 + 2.0 * h2 + h0 * h2)
  L0 => L(:, :, iz - 1, 0); L1 => L(:, :, iz, 0); L2 => L(:, :, iz + 1, 0)
  !$OMP PARALLEL
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
  DO ipin = 1, nxy
    DO ig = 1, ng
      L(ig, ipin, iz, 1) = - 2.0 * (L0(ig, ipin) * n1(1) + L1(ig, ipin) * (n1(2) - n1(1)) - L2(ig, ipin) * n1(2)) / d1
      L(ig, ipin, iz, 2) = 2.0 * (L0(ig, ipin) * n2(1) - L1(ig, ipin) * (n2(2) + n2(1)) + L2(ig, ipin) * n2(2)) / d2
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDDO

END SUBROUTINE

SUBROUTINE cuSetSourceCoeff(cuAxial, cuDevice, eigv, order)

IMPLICIT NONE

TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice
REAL :: eigv
INTEGER :: order

REAL(GPU_NODAL_PRECISION), POINTER, DEVICE :: csrVal(:)
INTEGER, POINTER, DEVICE :: csrRowPtr(:), csrColIdx(:)
INTEGER :: nr, nc, nnz
INTEGER :: ierr
REAL(GPU_NODAL_PRECISION) :: reigv

reigv = 1.0 / eigv

csrVal => cuAxial%S%d_csrVal
csrRowPtr => cuAxial%S%d_csrRowPtr
csrColIdx => cuAxial%S%d_csrColIdx
nr = cuAxial%S%nr; nc = cuAxial%S%nc; nnz = cuAxial%S%nnz
  
ierr = cusparseCsrmv(cuDevice%mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, one,                   &
                     cuAxial%S%descr, csrVal, csrRowPtr, csrColIdx, cuAxial%phiCoeff(:, :, :, order),               &
                     zero, cuAxial%srcCoeff(:, :, :, order))

IF (order .LE. 2) THEN
  CALL cuVectorOp('-', nr, cuAxial%srcCoeff(:, :, :, order), cuAxial%lkgCoeff(:, :, :, order),                      &
                  cuAxial%srcCoeff(:, :, :, order), cuDevice%myStream)
ENDIF

IF (.NOT. cuDevice%lFuel) RETURN

csrVal => cuAxial%Chi%d_csrVal
csrRowPtr => cuAxial%Chi%d_csrRowPtr
csrColIdx => cuAxial%Chi%d_csrColIdx
nr = cuAxial%Chi%nr; nc = cuAxial%Chi%nc; nnz = cuAxial%Chi%nnz

ierr = cusparseCsrmv(cuDevice%mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, reigv,                 &
                     cuAxial%Chi%descr, csrVal, csrRowPtr, csrColIdx, cuAxial%psiCoeff(:, :, order),                &
                     one, cuAxial%srcCoeff(:, :, :, order))
                     
END SUBROUTINE

SUBROUTINE cuGetIncomingCurrent(cuAxial, cuDevice)

IMPLICIT NONE

TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice

INTEGER :: n, ng, nxy, nzCMFD
INTEGER :: ig, ipin, iz, ierr
REAL, ALLOCATABLE, DEVICE :: myJout(:, :, :), neighJin(:, :, :)

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
n = ng * nxy

ALLOCATE(myJout(ng, nxy, 2))
ALLOCATE(neighJin(ng, nxy, 2))

ierr = cudaMemcpyAsync(myJout(:, :, bottom), cuAxial%Jout(:, :, 1, bottom, out), n, cudaMemcpyDeviceToDevice,       &
                       cuDevice%myStream)
ierr = cudaMemcpyAsync(myJout(:, :, top), cuAxial%Jout(:, :, nzCMFD, top, out), n, cudaMemcpyDeviceToDevice,        &
                       cuDevice%myStream)

DO iz = 1, nzCMFD
  IF (iz .NE. 1) THEN
    ierr = cudaMemcpyAsync(cuAxial%Jout(:, :, iz - 1, top, in), cuAxial%Jout(:, :, iz, bottom, out), n,             &
                           cudaMemcpyDeviceToDevice, cuDevice%myStream)
  END IF
  IF (iz .NE. nzCMFD) THEN
    ierr = cudaMemcpyAsync(cuAxial%Jout(:, :, iz + 1, bottom, in), cuAxial%Jout(:, :, iz, top, out), n,             &
                           cudaMemcpyDeviceToDevice, cuDevice%myStream)
  END IF
END DO

CALL cuGetNeighbor(n, myJout, neighJin, MPI_CUDA_COMM, cuDevice%myStream)

IF (cuDevice%lBottom) THEN
  IF (cuGeometry%AxBC(bottom) .EQ. Void) neighJin(:, :, bottom) = 0.0
  IF (cuGeometry%AxBC(bottom) .EQ. Reflective) THEN
    ierr = cudaMemcpyAsync(neighJin(:, :, bottom), myJout(:, :, bottom), n, cudaMemcpyDeviceToDevice,               &
                           cuDevice%myStream)
  END IF
END IF

IF (cuDevice%lTop) THEN
  IF (cuGeometry%AxBC(top) .EQ. Void) neighJin(:, :, top) = 0.0
  IF (cuGeometry%AxBC(top) .EQ. Reflective) THEN
    ierr = cudaMemcpyAsync(neighJin(:, :, top), myJout(:, :, top), n, cudaMemcpyDeviceToDevice,                     &
                           cuDevice%myStream)
  END IF
END IF

ierr = cudaMemcpyAsync(cuAxial%Jout(:, :, 1, bottom, in), neighJin(:, :, bottom), n, cudaMemcpyDeviceToDevice,      &
                       cuDevice%myStream)
ierr = cudaMemcpyAsync(cuAxial%Jout(:, :, nzCMFD, top, in), neighJin(:, :, top), n, cudaMemcpyDeviceToDevice,       &
                       cuDevice%myStream)
                       
ierr = cudaStreamSynchronize(cuDevice%myStream)

DEALLOCATE(myJout)
DEALLOCATE(neighJin)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuSENM(cuGeometry, cuDevice, phiCoeff, srcCoeff, Jout, D, sigR)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuGeometry_Type) :: cuGeometry
TYPE(cuDevice_Type) :: cuDevice
REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: srcCoeff(:, :, :, 0 :)
REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: D(:, :, :)
REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: sigR(:, :, :)
REAL(GPU_NODAL_PRECISION), DEVICE :: phiCoeff(:, :, :, 0 :)
REAL, DEVICE :: Jout(:, :, :, :, :)

INTEGER :: ig, iz, izf, ipin
REAL(GPU_NODAL_PRECISION) :: k, kinv, kinv2, kinv4, h, alpha(2), beta, tau, sigD, xsD, xsR
REAL(GPU_NODAL_PRECISION) :: a(0 : 4), q(0 : 4), c(6)
REAL(GPU_NODAL_PRECISION) :: phi(2), j(2), jL(2), jR(2)
REAL(GPU_NODAL_PRECISION) :: coshk, sinhk

ig = mod((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1, ng) + 1
ipin = ((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1) / ng + 1
izf = blockIdx%y

IF (ipin .GT. nxyc) RETURN

a(0) = phiCoeff(ig, ipin, izf, 0)
q = srcCoeff(ig, ipin, izf, :)
jL(in) = Jout(ig, ipin, izf, bottom, in)
jL(out) = Jout(ig, ipin, izf, bottom, out)
jR(in) = Jout(ig, ipin, izf, top, in)
jR(out) = Jout(ig, ipin, izf, top, out)
xsD = D(ig, ipin, izf)
xsR = sigR(ig, ipin, izf)
h = cuGeometry%hzfm(izf + cuDevice%myzbf - 1)

!--- Constants
beta = xsD / h
sigD = 4.0 * xsD / h ** 2
k = h / 2.0 * sqrt(xsR / xsD)
kinv = 1.0 / k; kinv2 = kinv ** 2; kinv4 = kinv ** 4
coshk = cosh(k); sinhk = sinh(k)
tau = k * sinhk / (4.0 * beta * k * sinhk + coshk - sinhk * kinv)

!--- Legendre Expansion Coefficients
c(1) = (1.0 / xsR) * (q(1) + 15.0 * kinv2 * q(3))
c(2) = (1.0 / xsR) * (q(2) + 35.0 * kinv2 * q(4))
c(3) = q(3) / xsR
c(4) = q(4) / xsR
!--- Average Flux
a(0) = q(0) + sigD * (2.0 * tau * (jL(in) + jR(in)) + (3.0 - tau * (1.0 + 12.0 * beta)) * c(2) + (10.0 - tau * (1.0 + 40.0 * beta)) * c(4))
a(0) = a(0) / (xsR + sigD * tau)
!--- Hyperbolic Term Coefficients
c(5) = 2.0 * (jR(in) - jL(in)) - (1.0 + 4.0 * beta) * c(1) - (1.0 + 24.0 * beta) * c(3)
c(5) = c(5) / (4.0 * beta * k * coshk + sinhk)
c(6) = - a(0) + 2.0 * (jR(in) + jL(in)) - (1.0 + 12.0 * beta) * c(2) - (1.0 + 40.0 * beta) * c(4)
c(6) = c(6) / (4.0 * beta * k * sinhk + coshk - sinhk * kinv)
!--- Approximated Flux Coefficients
a(1) = c(1) + (3.0 * kinv) * (coshk - sinhk * kinv) * c(5)
a(2) = c(2) - 5.0 * (3.0 * coshk * kinv2 - (1.0 + 3.0 * kinv2) * sinhk * kinv) * c(6)
a(3) = c(3) + (7.0 * kinv) * ((1.0 + 15.0 * kinv2) * coshk - (6.0 + 15.0 * kinv2) * sinhk * kinv) * c(5)
a(4) = c(4) - 9.0 * (5.0 * kinv2 * (2.0 + 21.0 * kinv2) * coshk - (1.0 + 45.0 * kinv2 + 105.0 * kinv4) * sinhk * kinv) * c(6)
!--- Outgoing Current
alpha(1) = - beta * k * coshk + 0.25 * sinhk
alpha(2) = - beta * k * sinhk + 0.25 * (coshk - sinhk * kinv)
phi(1) = - c(1) + c(2) - c(3) + c(4)
phi(2) = c(1) + c(2) + c(3) + c(4)
j(1) = - 2.0 * beta * (c(1) - 3.0 * c(2) + 6.0 * c(3) - 10.0 * c(4))
j(2) = - 2.0 * beta * (c(1) + 3.0 * c(2) + 6.0 * c(3) + 10.0 * c(4))
jL(out) = - alpha(1) * c(5) + alpha(2) * c(6) + 0.25 * (a(0) + phi(1)) - 0.5 * j(1)
jR(out) = alpha(1) * c(5) + alpha(2) * c(6) + 0.25 * (a(0) + phi(2)) + 0.5 * j(2)

phiCoeff(ig, ipin, izf, :) = a
Jout(ig, ipin, izf, bottom, out) = jL(out)
Jout(ig, ipin, izf, top, out) = jR(out)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuNEM(cuGeometry, cuDevice, phiCoeff, srcCoeff, Jout, D, sigR)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuGeometry_Type) :: cuGeometry
TYPE(cuDevice_Type) :: cuDevice
REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: srcCoeff(:, :, :, 0 :)
REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: D(:, :, :)
REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: sigR(:, :, :)
REAL(GPU_NODAL_PRECISION), DEVICE :: phiCoeff(:, :, :, 0 :)
REAL, DEVICE :: Jout(:, :, :, :, :)

INTEGER :: ig, izf, ipin
REAL(GPU_NODAL_PRECISION) :: hinv, beta, sigD, xsD, xsR
REAL(GPU_NODAL_PRECISION) :: a(0 : 4), q(0 : 4), c(4)
REAL(GPU_NODAL_PRECISION) :: jL(2), jR(2)

ig = mod((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1, ng) + 1
ipin = ((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1) / ng + 1
izf = blockIdx%y

IF (ipin .GT. nxyc) RETURN

a(0) = phiCoeff(ig, ipin, izf, 0)
q = srcCoeff(ig, ipin, izf, :)
jL(in) = Jout(ig, ipin, izf, bottom, in)
jL(out) = Jout(ig, ipin, izf, bottom, out)
jR(in) = Jout(ig, ipin, izf, top, in)
jR(out) = Jout(ig, ipin, izf, top, out)
xsD = D(ig, ipin, izf)
xsR = sigR(ig, ipin, izf)
hinv = 1.0 / cuGeometry%hzfm(izf + cuDevice%myzbf - 1)

!--- Constants
beta = xsD * hinv
sigD = xsD * hinv ** 2
c(1) = 6.0 * beta / (1.0 + 12.0 * beta)
c(2) = - 8.0 * beta / ((1.0 + 4.0 * beta) * (1.0 + 12.0 * beta))
c(3) = (1.0 - 48.0 * beta ** 2) / ((1.0 + 4.0 * beta) * (1.0 + 12.0 * beta))
c(4) = 6.0 * beta / (1.0 + 4.0 * beta)

!--- Flux Coefficients
a(1) = jR(out) + jR(in) - jL(out) - jL(in)
a(2) = a(0) - (jR(out) + jR(in) + jL(out) + jL(in))
a(3) = (5.0 * q(1) + 3.0 * q(3) - 5.0 * a(1) * xsR) / (3.0 * (60.0 * sigD + xsR))
a(4) = (- 7.0 * q(2) + 3.0 * q(4) + 7.0 * a(2) * xsR) / (3.0 * (140.0 * sigD + xsR))
!--- Average Flux
a(0) = q(0) - (2.0 * a(4) * c(1) - (1.0 - c(2) - c(3)) * (jL(in) + jR(in))) * hinv
a(0) = a(0) / (xsR + 2.0 * c(1) * hinv)
!--- Outgoing Current
jL(out) = c(1) * (a(0) + a(4)) + c(3) * jL(in) + c(2) * jR(in) - c(4) * a(3)
jR(out) = c(1) * (a(0) + a(4)) + c(2) * jL(in) + c(3) * jR(in) + c(4) * a(3)

phiCoeff(ig, ipin, izf, :) = a
Jout(ig, ipin, izf, bottom, out) = jL(out)
Jout(ig, ipin, izf, top, out) = jR(out)

END SUBROUTINE

END MODULE

#endif