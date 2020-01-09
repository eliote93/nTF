#include <CUDADEFINES.h>

#ifdef __PGI

MODULE CUDA_AXMOC

USE TYPEDEF,        ONLY : XsMac_Type
USE CUDA_MASTER
USE CUDA_UTIL
IMPLICIT NONE

REAL(GPU_NODAL_PRECISION), PARAMETER, PRIVATE :: zero = 0.0, one = 1.0
TYPE(XsMac_Type), POINTER, PRIVATE :: XsMac(:)

PRIVATE
PUBLIC :: cuAllocFlatMOC, cuAllocLinearMOC
PUBLIC :: cuFlatMOCDriver, cuLinearMOCDriver
PUBLIC :: cuSetSourceOperator

CONTAINS

!--- Public Routines ------------------------------------------------------------------------------

SUBROUTINE cuAllocFlatMOC(cuAxial, cuDevice)

IMPLICIT NONE

TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice

INTEGER :: ng, nxy, nzSub, nPolar1D

CALL cuSetSubmesh(cuDevice)

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzSub = cuDevice%nzSub
nPolar1D = cuGeometry%nPolar1D

ALLOCATE(XsMac(cuCntl%nCMFDHelper))

ALLOCATE(cuAxial%phiShape(nxy, ng, nzSub)); cuAxial%phiShape = 1.0
ALLOCATE(cuAxial%phiCoeff(nxy, ng, nzSub, 0 : 0)); cuAxial%phiCoeff = 0.0
ALLOCATE(cuAxial%phim(nxy, ng, nzSub)); cuAxial%phim = 0.0
ALLOCATE(cuAxial%PhiAngIn(nPolar1D, nxy, ng, 2))
ALLOCATE(cuAxial%PhiAngOut(nPolar1D, nxy, ng, 2))

ALLOCATE(cuAxial%psiCoeff(nxy, nzSub, 0 : 0))
ALLOCATE(cuAxial%lkgCoeff(nxy, ng, nzSub, 0 : 0))
ALLOCATE(cuAxial%srcCoeff(nxy, ng, nzSub, 0 : 0))
ALLOCATE(cuAxial%srcm(nxy, ng, nzSub))

ALLOCATE(cuAxial%xst(nxy, ng, nzSub))

END SUBROUTINE

SUBROUTINE cuDeallocFlatMOC(cuAxial)

IMPLICIT NONE

TYPE(cuAxial_Type) :: cuAxial

DEALLOCATE(cuAxial%phiCoeff)
DEALLOCATE(cuAxial%psiCoeff)
DEALLOCATE(cuAxial%lkgCoeff)
DEALLOCATE(cuAxial%srcCoeff)
DEALLOCATE(cuAxial%srcm)
DEALLOCATE(cuAxial%xst)

END SUBROUTINE

SUBROUTINE cuAllocLinearMOC(cuAxial, cuDevice)

IMPLICIT NONE

TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice

INTEGER :: ng, nxy, nzSub, nPolar1D

! CALL cuSetSubmesh(cuDevice)
!
! ng = cuGeometry%ng
! nxy = cuGeometry%nxyc
! nzSub = cuDevice%nzSub
! nPolar1D = cuGeometry%nPolar1D
!
! ALLOCATE(cuAxial%phiShape(ng, nxy, nzSub)); cuAxial%phiShape = 1.0
! ALLOCATE(cuAxial%phiCoeff(ng, nxy, nzSub, 0 : 1)); cuAxial%phiCoeff = 0.0
! ALLOCATE(cuAxial%phim(ng, nxy, nzSub))
! ALLOCATE(cuAxial%srcCoeff(ng, nxy, nzSub, 0 : 1))
! ALLOCATE(cuAxial%lkgCoeff(ng, nxy, nzSub, 0 : 0))
! ALLOCATE(cuAxial%psiCoeff(nxy, nzSub, 0 : 1))
! ALLOCATE(cuAxial%PhiAngIn(nPolar1D, ng, nxy, 2))
! ALLOCATE(cuAxial%PhiAngOut(nPolar1D, ng, nxy, 2))
! ALLOCATE(cuAxial%xst(ng, nxy, nzSub))
! ALLOCATE(cuAxial%E1(ng, nxy, nPolar1D, nzSub))
! ALLOCATE(cuAxial%E3(ng, nxy, nPolar1D, nzSub))
! ALLOCATE(cuAxial%R1(ng, nxy, nPolar1D, nzSub))
! ALLOCATE(cuAxial%R3(ng, nxy, nPolar1D, nzSub))

END SUBROUTINE

SUBROUTINE cuFlatMOCDriver(cuCMFD, cuAxial, cuDevice, eigv)
USE CUDA_FLATMOC
USE TIMER,          ONLY : nTracer_dclock,      TimeChk
IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice
REAL :: eigv

TYPE(dim3) :: Blocks, Threads
INTEGER :: ng, nxy
INTEGER :: sharedMemSize
INTEGER :: iter, ierr
INTEGER(KIND = cuda_stream_kind) :: myStream
REAL :: AxNSolverBeg, AxNSolverEnd

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
myStream = cuDevice%myStream

sharedMemSize = GPU_NODAL_PRECISION * cuGeometry%nPolar1D * cuDevice%sharedMemoryDim

Threads = dim3(cuDevice%cuThreadPerBlock, 1, 1)
Blocks = dim3(ng * nxy / Threads%x + 1, 1, 1)

AxNSolverBeg = nTracer_dclock(.FALSE., .FALSE.)

CALL cuSetCrossSection(cuCMFD, cuAxial, cuDevice)
CALL cuSetFlux(cuCMFD, cuAxial, cuDevice)
CALL cuSetPsi(cuAxial, cuDevice, 0)
CALL cuSetLeakage(cuCMFD, cuAxial, cuDevice)

!$ACC HOST_DATA USE_DEVICE(cuGeometry, cuDevice)
CALL cuLeakageSplit <<< Blocks, Threads, 0, myStream >>>                                                            &
                    (cuGeometry, cuDevice, cuAxial%phiCoeff(:, :, :, 0), cuAxial%lkgCoeff(:, :, :, 0), cuAxial%xst)
!$ACC END HOST_DATA

ierr = cudaStreamSynchronize(cuDevice%myStream)
AxNSolverEnd = nTracer_dclock(.FALSE., .FALSE.)
TimeChk%AxNSolverTime = TimeChk%AxNSolverTime + (AxNSolverEnd - AxNSolverBeg)

DO iter = 1, 10
  AxNSolverBeg = nTracer_dclock(.FALSE., .FALSE.)
  CALL cuSetSource(cuAxial, cuDevice, eigv, 0)
  !$ACC HOST_DATA USE_DEVICE(cuDevice)
  CALL cuFlatRayTrace1st <<< Blocks, Threads, sharedMemSize, myStream >>>                                           &
                         (cuDevice, cuAxial%phiCoeff(:, :, :, 0), cuAxial%phim, cuAxial%xst,                        &
                          cuAxial%srcCoeff(:, :, :, 0), cuAxial%srcm, cuAxial%PhiAngIn, cuAxial%PhiAngOut)
  CALL cuSetFluxMoment <<< Blocks, Threads, 0, myStream >>>                                                         &
                       (cuDevice, cuAxial%phiCoeff(:, :, :, 0), cuAxial%phim, cuAxial%srcCoeff(:, :, :, 0),         &
                        cuAxial%srcm, cuAxial%xst)
  !$ACC END HOST_DATA
  ierr = cudaStreamSynchronize(cuDevice%myStream)
  AxNSolverEnd = nTracer_dclock(.FALSE., .FALSE.)
  TimeChk%AxNSolverTime = TimeChk%AxNSolverTime + (AxNSolverEnd - AxNSolverBeg)
  CALL cuCommBoundaryFlux(cuAxial, cuDevice)
ENDDO

AxNSolverBeg = nTracer_dclock(.FALSE., .FALSE.)
!$ACC HOST_DATA USE_DEVICE(cuDevice)
CALL cuFlatRayTrace2nd <<< Blocks, Threads, sharedMemSize, myStream >>>                                             &
                       (cuDevice, cuAxial%xst, cuAxial%srcCoeff(:, :, :, 0), cuAxial%srcm, cuAxial%PhiAngIn,        &
                        cuAxial%PhiAngOut, cuAxial%Jout)
CALL cuSetFluxShape <<< Blocks, Threads, 0, myStream >>>                                                            &
                    (cuDevice, cuAxial%phiCoeff(:, :, :, 0), cuAxial%phic, cuAxial%phiShape)
!$ACC END HOST_DATA
ierr = cudaStreamSynchronize(cuDevice%myStream)
AxNSolverEnd = nTracer_dclock(.FALSE., .FALSE.)
TimeChk%AxNSolverTime = TimeChk%AxNSolverTime + (AxNSolverEnd - AxNSolverBeg)
CALL cuCommBoundaryFlux(cuAxial, cuDevice)

END SUBROUTINE

SUBROUTINE cuLinearMOCDriver(cuCMFD, cuAxial, cuDevice, eigv)
USE CUDA_LINMOC
IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice
REAL :: eigv

! TYPE(dim3) :: Blocks, Threads
! INTEGER :: ng, nxy
! INTEGER :: sharedMemSize
! INTEGER :: iter
! INTEGER(KIND = cuda_stream_kind) :: myStream
!
! ng = cuGeometry%ng
! nxy = cuGeometry%nxyc
! myStream = cuDevice%myStream
!
! sharedMemSize = GPU_NODAL_PRECISION * cuGeometry%nPolar1D * cuDevice%sharedMemoryDim
!
! Threads = dim3(cuDevice%cuThreadPerBlock, 1, 1)
! Blocks = dim3(ng * nxy / Threads%x + 1, 1, 1)
!
! ! CALL cuSetSourceOperator(cuCMFD, cuAxial, cuDevice)
! CALL cuSetCrossSection(cuCMFD, cuAxial, cuDevice)
! CALL cuSetFlux(cuCMFD, cuAxial, cuDevice)
! CALL cuSetPsi(cuAxial, cuDevice, 0)
! CALL cuSetLeakage(cuCMFD, cuAxial, cuDevice)
!
! !$ACC HOST_DATA USE_DEVICE(cuDevice)
! CALL cuSetRTCoeff <<< Blocks, Threads, 0, myStream >>>                                                              &
!                   (cuDevice, cuAxial%xst, cuAxial%E1, cuAxial%E3, cuAxial%R1, cuAxial%R3)
! !$ACC END HOST_DATA
!
! DO iter = 1, 10
!   CALL cuSetPsi(cuAxial, cuDevice, 1)
!   CALL cuSetSource(cuAxial, cuDevice, eigv, 0)
!   CALL cuSetSource(cuAxial, cuDevice, eigv, 1)
!   !$ACC HOST_DATA USE_DEVICE(cuDevice)
!   CALL cuLinearRayTrace1st <<< Blocks, Threads, sharedMemSize, myStream >>>                                         &
!                            (cuDevice, cuAxial%phiCoeff, cuAxial%phim, cuAxial%srcCoeff, cuAxial%E1, cuAxial%E3,     &
!                             cuAxial%R1, cuAxial%R3, cuAxial%PhiAngIn, cuAxial%PhiAngOut)
!   CALL cuSetFluxMoment <<< Blocks, Threads, 0, myStream >>>                                                         &
!                        (cuDevice, cuAxial%phiCoeff, cuAxial%phim, cuAxial%xst, cuAxial%srcCoeff(:, :, :, 0))
!   !$ACC END HOST_DATA
!   CALL cuCommBoundaryFlux(cuAxial, cuDevice)
! ENDDO
!
! !$ACC HOST_DATA USE_DEVICE(cuDevice)
! CALL cuLinearRayTrace2nd <<< Blocks, Threads, sharedMemSize, myStream >>>                                           &
!                          (cuDevice, cuAxial%phiCoeff, cuAxial%phim, cuAxial%srcCoeff, cuAxial%E1, cuAxial%E3,       &
!                           cuAxial%R1, cuAxial%R3, cuAxial%PhiAngIn, cuAxial%PhiAngOut, cuAxial%Jout)
! CALL cuSetFluxMoment <<< Blocks, Threads, 0, myStream >>>                                                           &
!                      (cuDevice, cuAxial%phiCoeff, cuAxial%phim, cuAxial%xst, cuAxial%srcCoeff(:, :, :, 0))
! CALL cuGetShape <<< Blocks, Threads, 0, myStream >>> (cuDevice, cuAxial%phiCoeff(:, :, :, 0), cuAxial%phiShape)
! !$ACC END HOST_DATA
! CALL cuCommBoundaryFlux(cuAxial, cuDevice)

END SUBROUTINE

SUBROUTINE cuSetSourceOperator(CoreInfo, FmInfo, GroupInfo, cuCMFD, cuAxial, cuDevice)
USE TYPEDEF,        ONLY : CoreInfo_Type,       FmInfo_Type,        GroupInfo_Type,                                 &
                           FxrInfo_Type,        Pin_Type
USE MacXsLib_Mod,   ONLY : MacP1XsScatMatrix
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice

TYPE(superPin_Type), POINTER :: superPin(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
INTEGER :: ig, igs, igb, ige, ifxr, ixy, ixy_map, ipin, ipin_map, iz, ierr, tid
INTEGER :: n, ng, nxy
INTEGER :: myzb, myze
INTEGER, POINTER :: pinMap(:)
REAL(GPU_NODAL_PRECISION), ALLOCATABLE :: S0(:, :, :, :), S1(:, :, :, :), F(:, :, :), Chi(:, :, :)

Pin => CoreInfo%Pin
Fxr => FmInfo%Fxr
superPin => cuGeometry%superPin
ng = cuCMFD%ng
nxy = cuGeometry%nxyc
myzb = cuDevice%myzb
myze = cuDevice%myze
pinMap => cuGeometry%pinMap

ALLOCATE(S0(nxy, ng, myzb : myze, ng)); S0 = 0.0
ALLOCATE(S1(nxy, ng, myzb : myze, ng)); S1 = 0.0
ALLOCATE(F(nxy, myzb : myze, ng))
ALLOCATE(Chi(nxy, ng, myzb : myze))

!$OMP PARALLEL PRIVATE(ipin_map)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO iz = myzb, myze
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    DO ig = 1, ng
      igb = cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%ib
      ige = cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%ie
      DO igs = igb, ige
        IF (igs .EQ. ig) THEN
          IF (cuGeometry%lH2OCell(iz, ipin)) THEN
            S0(ipin, ig, iz, igs) = cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%self                                         &
                                    + (cuCMFD%PinXS(ipin_map, iz)%XSt(ig) - cuCMFD%PinXS(ipin_map, iz)%XStr(ig))
          ELSE
            S0(ipin, ig, iz, igs) = cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%self
            IF (cuGeometry%lRefCell(iz, ipin)) THEN
              IF (S0(ipin, ig, iz, igs) .LT. 0.0) S0(ipin, ig, iz, igs) = 0.0
            ENDIF
          ENDIF
        ELSE
          S0(ipin, ig, iz, igs) = cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%from(igs)
        ENDIF
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

IF (.NOT. ASSOCIATED(XsMac)) ALLOCATE(XsMac(cuCntl%nCMFDHelper))

!$OMP PARALLEL PRIVATE(tid, ifxr, ipin, ixy_map)
tid = omp_get_thread_num() + 1

!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO iz = myzb, myze
  DO ixy = 1, nxy
    IF (.NOT. cuGeometry%lH2OCell(iz, ixy)) CYCLE
    ixy_map = pinMap(ixy)
    ipin = superPin(ixy_map)%pin(1)
    ifxr = Pin(ipin)%FxrIdxSt
    CALL MacP1XsScatMatrix(XsMac(tid), Fxr(ifxr, iz), 1, ng, ng, GroupInfo)
    DO ig = 1, ng
      DO igs = 1, ng
        S1(ixy, ig, iz, igs) = XsMac(tid)%XsMacP1Sm(igs, ig)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(ipin_map)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO iz = myzb, myze
    DO ipin = 1, nxy
      ipin_map = pinMap(ipin)
      F(ipin, iz, ig) = cuCMFD%PinXS(ipin_map, iz)%XSnf(ig)
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(ipin_map)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO iz = myzb, myze
  DO ig = 1, ng
    DO ipin = 1, nxy
      ipin_map = pinMap(ipin)
      Chi(ipin, ig, iz) = cuCMFD%PinXS(ipin_map, iz)%Chi(ig)
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

n = ng * ng * nxy * (myze - myzb + 1)
ierr = cudaMemcpy(cuAxial%S0, S0, n, cudaMemcpyHostToDevice)
ierr = cudaMemcpy(cuAxial%S1, S1, n, cudaMemcpyHostToDevice)

n = ng * nxy * (myze - myzb + 1)
ierr = cudaMemcpy(cuAxial%Fa, F, n, cudaMemcpyHostToDevice)
ierr = cudaMemcpy(cuAxial%Chia, Chi, n, cudaMemcpyHostToDevice)

DEALLOCATE(S0, S1, F, Chi)

END SUBROUTINE

!--- Private Routines -----------------------------------------------------------------------------

SUBROUTINE cuSetSubmesh(cuDevice)

IMPLICIT NONE

TYPE(cuDevice_Type) :: cuDevice

INTEGER :: nzCMFD, nzSub
INTEGER :: iz, izc, izf
INTEGER :: myzb, myze, myzbf, myzef
REAL, POINTER :: hzfm(:)

nzCMFD = cuDevice%nzCMFD
myzb = cuDevice%myzb
myze = cuDevice%myze
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
hzfm => cuGeometry%hzfm

nzSub = 0
DO iz = myzbf, myzef
  cuDevice%nDiv(iz) = INT(cuGeometry%hzfm(iz) / cuCntl%MOCHeight) + 1
  nzSub = nzSub + cuDevice%nDiv(iz)
ENDDO

ALLOCATE(cuDevice%hzSub(nzSub))
ALLOCATE(cuDevice%cmMap(nzSub), cuDevice%fmMap(nzSub))
ALLOCATE(cuDevice%cmSubrange(myzb : myze, 2), cuDevice%fmSubrange(myzbf : myzef, 2))

izf = 0
DO izc = myzb, myze
  cuDevice%cmSubrange(izc, 1) = izf + 1
  DO iz = cuGeometry%fmRange(izc, 1), cuGeometry%fmRange(izc, 2)
    cuDevice%cmMap(izf + 1 : izf + cuDevice%nDiv(iz)) = izc
    cuDevice%fmMap(izf + 1 : izf + cuDevice%nDiv(iz)) = iz
    cuDevice%fmSubrange(iz, 1) = izf + 1
    cuDevice%fmSubrange(iz, 2) = izf + cuDevice%nDiv(iz)
    cuDevice%hzSub(izf + 1 : izf + cuDevice%nDiv(iz)) = hzfm(iz) / cuDevice%nDiv(iz)
    izf = izf + cuDevice%nDiv(iz)
  ENDDO
  cuDevice%cmSubrange(izc, 2) = izf
ENDDO

cuDevice%nzSub = nzSub

END SUBROUTINE

SUBROUTINE cuSetCrossSection(cuCMFD, cuAxial, cuDevice)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice

INTEGER :: ig, ipin, ipin_map, iz, izf, ierr
INTEGER :: n, ng, nxy, nzSub
INTEGER :: myzb, myze
INTEGER, POINTER :: pinMap(:)
REAL(GPU_NODAL_PRECISION), ALLOCATABLE :: xst(:, :, :)

ng = cuCMFD%ng
nxy = cuGeometry%nxyc
nzSub = cuDevice%nzSub
myzb = cuDevice%myzb
myze = cuDevice%myze
pinMap => cuGeometry%pinMap
n = ng * nxy * nzSub

ALLOCATE(xst(nxy, ng, nzSub))

DO iz = myzb, myze
  !$OMP PARALLEL PRIVATE(ipin_map)
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
  DO izf = cuDevice%cmSubrange(iz, 1), cuDevice%cmSubrange(iz, 2)
    DO ig = 1, ng
      DO ipin = 1, nxy
        ipin_map = pinMap(ipin)
        IF (cuGeometry%lH2OCell(iz, ipin)) THEN
          xst(ipin, ig, izf) = cuCMFD%PinXS(ipin_map, iz)%XSt(ig)
        ELSE
          xst(ipin, ig, izf) = cuCMFD%PinXS(ipin_map, iz)%XStr(ig)
          IF (cuGeometry%lRefCell(iz, ipin)) THEN
            IF (cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%self .LT. 0.0) THEN
              xst(ipin, ig, izf) = xst(ipin, ig, izf) - cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%self
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDDO

ierr = cudaMemcpy(cuAxial%xst, xst, n, cudaMemcpyHostToDevice)

DEALLOCATE(xst)

END SUBROUTINE

SUBROUTINE cuSetFlux(cuCMFD, cuAxial, cuDevice)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice

INTEGER :: ig, ipin, iz, izf, ierr
INTEGER :: n, ng, nxy, nzSub
INTEGER :: myzbf, myzef
REAL(GPU_NODAL_PRECISION), ALLOCATABLE :: phis(:, :, :)

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzSub = cuDevice%nzSub
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
n = ng * nxy * nzSub

ALLOCATE(phis(nxy, ng, nzSub))

DO iz = myzbf, myzef
  !$OMP PARALLEL
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
  DO izf = cuDevice%fmSubrange(iz, 1), cuDevice%fmSubrange(iz, 2)
    DO ipin = 1, nxy
      DO ig = 1, ng
        phis(ipin, ig, izf) = cuCMFD%h_phis8(ig, ipin, iz)
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDDO

ierr = cudaMemcpy(cuAxial%phiCoeff(:, :, :, 0), phis, n, cudaMemcpyHostToDevice)

DEALLOCATE(phis)

CALL cuVectorOp('*', n, cuAxial%phiShape, cuAxial%phiCoeff(:, :, :, 0), cuAxial%phiCoeff(:, :, :, 0), cuDevice%myStream)

END SUBROUTINE

SUBROUTINE cuSetPsi(cuAxial, cuDevice, order)

IMPLICIT NONE

TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice
INTEGER :: order

TYPE(dim3) :: Blocks, Threads
INTEGER(KIND = cuda_stream_kind) :: myStream
INTEGER :: nxy, nzSub

nxy = cuGeometry%nxyc
nzSub = cuDevice%nzSub
myStream = cuDevice%myStream

Threads = dim3(cuDevice%cuThreadPerBlock, 1, 1)
Blocks = dim3(nxy * nzSub / Threads%x + 1, 1, 1)

!$ACC HOST_DATA USE_DEVICE(cuDevice)
CALL cuSetPsiKernel <<< Blocks, Threads, 0, myStream >>>                                                            &
                    (cuDevice, cuAxial%Fa, cuAxial%phiCoeff, cuAxial%psiCoeff, order)
!$ACC END HOST_DATA

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuSetPsiKernel(cuDevice, F, phi, psi, order)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuDevice_Type), INTENT(IN) :: cuDevice
REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: F(:, cuDevice%myzb :, :)
REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: phi(:, :, :, 0 :)
REAL(GPU_NODAL_PRECISION), DEVICE :: psi(:, :, 0 :)
INTEGER, VALUE :: order

REAL(GPU_NODAL_PRECISION) :: psiSrc
INTEGER :: ig, ipin, iz, izf

ipin = mod((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1, nxyc) + 1
izf = ((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1) / nxyc + 1

IF (izf .GT. cuDevice%nzSub) RETURN

iz = cuDevice%cmMap(izf)

psiSrc = 0.0
DO ig = 1, ng
  psiSrc = psiSrc + F(ipin, iz, ig) * phi(ipin, ig, izf, order)
ENDDO

psi(ipin, izf, order) = psiSrc

END SUBROUTINE

SUBROUTINE cuSetLeakage(cuCMFD, cuAxial, cuDevice)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice

TYPE(superPin_Type), POINTER :: Pin(:)
INTEGER :: n, ng, nxy, nzSub
INTEGER :: myzbf, myzef
INTEGER :: ig, ibd, ipin, ipin_map, ineighpin, iz, izf, izc, ierr
INTEGER, POINTER :: pinMap(:), pinMapRev(:), planeMap(:)
REAL :: myphi, neighphi, lkg
REAL :: Dtil, Dhat
REAL, ALLOCATABLE :: radLkg(:, :, :)
REAL(GPU_NODAL_PRECISION), ALLOCATABLE :: lkgCoeff(:, :, :)

Pin => cuGeometry%superPin
ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzSub = cuDevice%nzSub
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
pinMap => cuGeometry%pinMap
pinMapRev => cuGeometry%pinMapRev
planeMap => cuGeometry%planeMap

ALLOCATE(radLkg(nxy, ng, myzbf - 1 : myzef + 1))
ALLOCATE(lkgCoeff(nxy, ng, nzSub))

DO iz = myzbf, myzef
  izc = planeMap(iz)
  !$OMP PARALLEL PRIVATE(ipin_map, ineighpin, myphi, neighphi, lkg, Dtil, Dhat)
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
  DO ig = 1, ng
    DO ipin = 1, nxy
      ipin_map = pinMap(ipin)
      myphi = cuCMFD%h_phis8(ig, ipin, iz)
      lkg = 0.0
      DO ibd = 1, 4
        ineighpin = Pin(ipin_map)%NeighIdx(ibd)
        ineighpin = pinMapRev(ineighpin)
        IF (ineighpin .EQ. Void) THEN
          neighphi = 0.0
        ELSEIF (ineighpin .EQ. Reflective) THEN
          neighphi = myphi
        ELSE
          neighphi = cuCMFD%h_phis8(ig, ineighpin, iz)
        ENDIF
        Dtil = cuCMFD%PinXS(ipin_map, izc)%Dtil(ibd, ig)
        Dhat = cuCMFD%PinXS(ipin_map, izc)%Dhat(ibd, ig)
        lkg = lkg - Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
      ENDDO
      radLkg(ipin, ig, iz) = lkg * cuGeometry%hzfm(iz) / cuGeometry%PinVolFm(ipin_map, iz)
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDDO

CALL cuGetNeighbor(ng * nxy, radLkg(:, :, myzbf), radLkg(:, :, myzbf - 1), MPI_CUDA_COMM, bottom)
CALL cuGetNeighbor(ng * nxy, radLkg(:, :, myzef), radLkg(:, :, myzef + 1), MPI_CUDA_COMM, top)

IF (cuDevice%lBottom) THEN
  IF (cuGeometry%AxBC(bottom) .EQ. Void) radLkg(:, :, myzbf - 1) = 0.0
  IF (cuGeometry%AxBC(bottom) .EQ. Reflective) radLkg(:, :, myzbf - 1) = radLkg(:, :, myzbf)
ENDIF

IF (cuDevice%lTop) THEN
  IF (cuGeometry%AxBC(top) .EQ. Void) radLkg(:, :, myzef + 1) = 0.0
  IF (cuGeometry%AxBC(top) .EQ. Reflective) radLkg(:, :, myzef + 1) = radLkg(:, :, myzef)
ENDIF

DO iz = myzbf, myzef
  CALL LeakageExpansion(cuDevice, radLkg(:, :, iz - 1), radLkg(:, :, iz), radLkg(:, :, iz + 1), lkgCoeff, iz)
ENDDO

n = ng * nxy * nzSub
ierr = cudaMemcpy(cuAxial%lkgCoeff, lkgCoeff, n, cudaMemcpyHostToDevice)

DEALLOCATE(radLkg)
DEALLOCATE(lkgCoeff)

END SUBROUTINE

SUBROUTINE LeakageExpansion(cuDevice, L0, L1, L2, lkgCoeff, iz)

IMPLICIT NONE

TYPE(cuDevice_Type) :: cuDevice
REAL :: L0(:, :), L1(:, :), L2(:, :)
REAL(GPU_NODAL_PRECISION) :: lkgCoeff(:, :, :)
INTEGER :: iz

REAL :: n0(3), n1(3), n2(3)
REAL :: d0, d1, d2
REAL :: h0, h1, h2, dh
REAL :: x0, x1
REAL :: a, b, c
INTEGER :: ig, ipin, izf
INTEGER :: ng, nxy

ng = cuGeometry%ng
nxy = cuGeometry%nxyc

h0 = cuGeometry%hzfm(iz - 1)
h1 = cuGeometry%hzfm(iz)
h2 = cuGeometry%hzfm(iz + 1)
dh = h1 / cuDevice%nDiv(iz)

n0(1) = h1 ** 3 + 2.0 * h1 ** 2 * h2 + h1 * h2 ** 2
n0(2) = 2.0 * h0 ** 2 * h1 + 3.0 * h0 * h1 ** 2 + h0 ** 2 * h2 + 3.0 * h0 * h1 * h2 + h0 * h2 ** 2
n0(3) = - h0 ** 2 * h1 - h0 * h1 ** 2

n1(1) = 2.0 * h1 ** 2 + 3.0 * h1 * h2 + h2 ** 2
n1(2) = h0 ** 2 - 3.0 * h1 ** 2 - 3.0 * h1 * h2 - h2 ** 2
n1(3) = - h0 ** 2 + h1 ** 2

n2(1) = h1 + h2
n2(2) = - h0 - 2.0 * h1 - h2
n2(3) = h0 + h1

d0 = (h1 + h2) * (h0 ** 2 + 2.0 * h0 * h1 + h1 ** 2 + h0 * h2 + h1 * h2)
d1 = (h0 + h1) * (h1 + h2) * (h0 + h1 + h2)
d2 = (h0 + h1) * (h1 + h2) * (h0 + h1 + h2)

!$OMP PARALLEL PRIVATE(a, b, c, x0, x1)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ig = 1, ng
  DO ipin = 1, nxy
    a = 3.0 * (L0(ipin, ig) * n2(1) + L1(ipin, ig) * n2(2) + L2(ipin, ig) * n2(3)) / d2
    b = - 2.0 * (L0(ipin, ig) * n1(1) + L1(ipin, ig) * n1(2) + L2(ipin, ig) * n1(3)) / d1
    c = (L0(ipin, ig) * n0(1) + L1(ipin, ig) * n0(2) + L2(ipin, ig) * n0(3)) / d0
    DO izf = cuDevice%fmSubrange(iz, 1), cuDevice%fmSubrange(iz, 2)
      x0 = (izf - cuDevice%fmSubrange(iz, 1)) * dh; x1 = x0 + dh
      lkgCoeff(ipin, ig, izf) = Integral(a, b, c, x0, x1) / dh
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

CONTAINS

FUNCTION Integral(a, b, c, x0, x1) RESULT(val)

IMPLICIT NONE

REAL :: a, b, c, x0, x1
REAL :: val

val = a * (x1 ** 3 - x0 ** 3) / 3.0 + b * (x1 ** 2 - x0 ** 2) / 2.0 + c * (x1 - x0)

END FUNCTION

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuLeakageSplit(cuGeometry, cuDevice, phis, lkg, xst)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuGeometry_Type) :: cuGeometry
TYPE(cuDevice_Type) :: cuDevice
REAL(GPU_NODAL_PRECISION) :: phis(:, :, :), lkg(:, :, :), xst(:, :, :)

INTEGER :: ig, ipin, iz, izf

ipin = mod((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1, nxyc) + 1
ig = ((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1) / nxyc + 1

IF (ig .GT. ng) RETURN

DO iz = cuDevice%myzb, cuDevice%myze
  !--- Skip Fuel Cells
  IF (.NOT. cuGeometry%lRefCell(iz, ipin)) CYCLE
  DO izf = cuDevice%cmSubrange(iz, 1), cuDevice%cmSubrange(iz, 2)
    IF (lkg(ipin, ig, izf) .GT. 0.0 .AND. phis(ipin, ig, izf) .GT. 0.0) THEN
      xst(ipin, ig, izf) = xst(ipin, ig, izf) + lkg(ipin, ig, izf) / phis(ipin, ig, izf)
      lkg(ipin, ig, izf) = 0.0
    ENDIF
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE cuSetSource(cuAxial, cuDevice, eigv, order)

IMPLICIT NONE

TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice
REAL :: eigv
INTEGER :: order

TYPE(dim3) :: Blocks, Threads
INTEGER(KIND = cuda_stream_kind) :: myStream
INTEGER :: ng, nxy, nzSub
REAL(GPU_NODAL_PRECISION) :: reigv

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzSub = cuDevice%nzSub
myStream = cuDevice%myStream
reigv = 1.0 / eigv

Threads = dim3(cuDevice%cuThreadPerBlock, 1, 1)
Blocks = dim3(ng * nxy / Threads%x + 1, nzSub, 1)

!$ACC HOST_DATA USE_DEVICE(cuDevice)
CALL cuSetSourceKernel <<< Blocks, Threads, 0, myStream >>>                                                         &
                       (cuDevice, cuAxial%S0, cuAxial%S1, cuAxial%Chia, cuAxial%phiCoeff, cuAxial%phim,             &
                        cuAxial%psiCoeff, cuAxial%lkgCoeff, cuAxial%xst, cuAxial%srcCoeff, cuAxial%srcm,            &
                        reigv, order)
!$ACC END HOST_DATA

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuSetSourceKernel(cuDevice, S0, S1, Chi, phis, phim, psi, lkg, xst, src, srcm, reigv, order)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuDevice_Type), INTENT(IN) :: cuDevice
REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: S0(:, :, cuDevice%myzb :, :), S1(:, :, cuDevice%myzb :, :), Chi(:, :, cuDevice%myzb :)
REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: phis(:, :, :, 0 :), phim(:, :, :), psi(:, :, 0 :), lkg(:, :, :, 0 :), xst(:, :, :)
REAL(GPU_NODAL_PRECISION), DEVICE :: src(:, :, :, 0 :), srcm(:, :, :)
REAL(GPU_NODAL_PRECISION), VALUE :: reigv
INTEGER, VALUE :: order

REAL(GPU_NODAL_PRECISION) :: scatSrc, fisSrc, momentSrc
INTEGER :: ig, igs, ipin, iz, izf
INTEGER :: gb, ge

ipin = mod((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1, nxyc) + 1
ig = ((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1) / nxyc + 1
izf = blockIdx%y

IF (ig .GT. ng) RETURN

iz = cuDevice%cmMap(izf)
gb = InScatRange(1, ig)
ge = InScatRange(2, ig)

scatSrc = 0.0
DO igs = gb, ge
  scatSrc = scatSrc + S0(ipin, ig, iz, igs) * phis(ipin, igs, izf, order)
ENDDO
fisSrc = Chi(ipin, ig, iz) * psi(ipin, izf, order)
src(ipin, ig, izf, order) = scatSrc + reigv * fisSrc

IF (order .EQ. 0) THEN
  momentSrc = 0.0
  DO igs = gb, ge
    momentSrc = momentSrc + S1(ipin, ig, iz, igs) * phim(ipin, igs, izf)
  ENDDO
  srcm(ipin, ig, izf) = momentSrc
ENDIF

IF (order .EQ. 0) THEN
  src(ipin, ig, izf, 0) = src(ipin, ig, izf, 0) - lkg(ipin, ig, izf, 0)
  src(ipin, ig, izf, 0) = src(ipin, ig, izf, 0) / xst(ipin, ig, izf)
  srcm(ipin, ig, izf) = 3.0 * srcm(ipin, ig, izf) / xst(ipin, ig, izf)
ELSEIF (order .EQ. 1) THEN
  src(ipin, ig, izf, 1) = src(ipin, ig, izf, 1) / xst(ipin, ig, izf) ** 2
ENDIF

END SUBROUTINE

SUBROUTINE cuCommBoundaryFlux(cuAxial, cuDevice)

IMPLICIT NONE

TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice

CALL cuInitCommBoundaryFlux(cuAxial, cuDevice)
CALL cuFinalizeCommBoundaryFlux(cuAxial, cuDevice)

END SUBROUTINE

SUBROUTINE cuInitCommBoundaryFlux(cuAxial, cuDevice)

IMPLICIT NONE

TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice

INTEGER :: n, ng, nxy, nPolar1D
INTEGER :: ierr

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nPolar1D = cuGeometry%nPolar1D
n = ng * nxy * nPolar1D

CALL InitMPIComm()

CALL MPIComm(n, n, cuAxial%PhiAngOut(:, :, :, bottom), cuAxial%PhiAngIn(:, :, :, bottom), bottom, MPI_CUDA_COMM)
CALL MPIComm(n, n, cuAxial%PhiAngOut(:, :, :, top), cuAxial%PhiAngIn(:, :, :, top), top, MPI_CUDA_COMM)

END SUBROUTINE

SUBROUTINE cuFinalizeCommBoundaryFlux(cuAxial, cuDevice)

IMPLICIT NONE

TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice

INTEGER :: n, ng, nxy, nPolar1D
INTEGER :: ierr

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nPolar1D = cuGeometry%nPolar1D
n = ng * nxy * nPolar1D

CALL FinalizeMPIComm()

IF (cuDevice%lBottom) THEN
  IF (cuGeometry%AxBC(bottom) .EQ. Void) cuAxial%PhiAngIn(:, :, :, bottom) = 0.0
  IF (cuGeometry%AxBC(bottom) .EQ. Reflective) THEN
    ierr = cudaMemcpy(cuAxial%PhiAngIn(:, :, :, bottom), cuAxial%PhiAngOut(:, :, :, bottom), n, cudaMemcpyDeviceToDevice)
  ENDIF
ENDIF

IF (cuDevice%lTop) THEN
  IF (cuGeometry%AxBC(top) .EQ. Void) cuAxial%PhiAngIn(:, :, :, top) = 0.0
  IF (cuGeometry%AxBC(top) .EQ. Reflective) THEN
    ierr = cudaMemcpy(cuAxial%PhiAngIn(:, :, :, top), cuAxial%PhiAngOut(:, :, :, top), n, cudaMemcpyDeviceToDevice)
  ENDIF
ENDIF

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuSetFluxShape(cuDevice, phis, phic, phiShape)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuDevice_Type) :: cuDevice
REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: phis(:, :, :)
REAL(GPU_NODAL_PRECISION), DEVICE :: phiShape(:, :, :)
REAL, DEVICE :: phic(:, :, :)

INTEGER :: ig, ipin, iz, izf
REAL(GPU_NODAL_PRECISION) :: phisum
LOGICAL :: lNegative

ipin = mod((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1, nxyc) + 1
ig = ((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1) / nxyc + 1

IF (ig .GT. ng) RETURN

DO iz = cuDevice%myzbf, cuDevice%myzef

  lNegative = .FALSE.

  phisum = 0.0
  DO izf = cuDevice%fmSubrange(iz, 1), cuDevice%fmSubrange(iz, 2)
    phisum = phisum + phis(ipin, ig, izf)
    IF (phis(ipin, ig, izf) .LT. 0.0) lNegative = .TRUE.
  ENDDO

  phisum = phisum / cuDevice%nDiv(iz)
  phic(ig, ipin, iz) = phisum

  IF (lNegative) THEN
    DO izf = cuDevice%fmSubrange(iz, 1), cuDevice%fmSubrange(iz, 2)
      phiShape(ipin, ig, izf) = 1.0
    ENDDO
  ELSE
    DO izf = cuDevice%fmSubrange(iz, 1), cuDevice%fmSubrange(iz, 2)
      phiShape(ipin, ig, izf) = phis(ipin, ig, izf) / phisum
    ENDDO
  ENDIF

ENDDO

END SUBROUTINE

END MODULE

#endif
