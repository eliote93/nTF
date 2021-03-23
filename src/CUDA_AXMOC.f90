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

PUBLIC :: cuTranFlatMOCDriver
PUBLIC :: cuSetTranSourceOperator
PUBLIC :: cuSetAxPsiNowStep

  CONTAINS

!--- Public Routines ------------------------------------------------------------------------------

SUBROUTINE cuSetAxPsiNowStep(inv_eigv0)
IMPLICIT NONE
REAL :: inv_eigv0

REAL(GPU_NODAL_PRECISION) :: ax_inv_eigv0
INTEGER :: nxy, nzSub
INTEGER :: n, ierr

nxy = cuGeometry%nxyc
nzSub = cuDevice%nzSub

CALL cuSetFlux(cuCMFD, cuAxial, cuDevice)
CALL cuSetPsi(cuAxial, cuDevice, 0)
n = nxy * nzSub
ax_inv_eigv0 = inv_eigv0
ierr = cublasScal(cuDevice%myblasHandle, n, ax_inv_eigv0, cuAxial%PsiCoeff(:,:,0), 1)
END SUBROUTINE

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

cuAxial%phic = 1.0
!ALLOCATE(cuAxial%phiShape(nxy, ng, nzSub)); cuAxial%phiShape = 1.0
ALLOCATE(cuAxial%phiCoeff(nxy, ng, nzSub, 0 : 0)); cuAxial%phiCoeff = 1.0
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

DO iter = 1, 20
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

SUBROUTINE cuTranFlatMOCDriver(TranInfo, TranCntl, eigv0)
USE TYPEDEF,        ONLY : TranInfo_Type,     TranCntl_Type
USE CUDA_FLATMOC
USE TIMER,          ONLY : nTracer_dclock,      TimeChk
IMPLICIT NONE
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
REAL :: eigv0

TYPE(dim3) :: Blocks, Threads
REAL(GPU_NODAL_PRECISION) :: inv_eigv0
INTEGER :: ng, nxy, nzSub
INTEGER :: sharedMemSize
INTEGER :: iter, ierr, n
INTEGER(KIND = cuda_stream_kind) :: myStream
REAL :: AxNSolverBeg, AxNSolverEnd

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzSub = cuDevice%nzSub
myStream = cuDevice%myStream

inv_eigv0 = 1./ eigv0

sharedMemSize = GPU_NODAL_PRECISION * cuGeometry%nPolar1D * cuDevice%sharedMemoryDim

Threads = dim3(cuDevice%cuThreadPerBlock, 1, 1)
Blocks = dim3(ng * nxy / Threads%x + 1, 1, 1)

AxNSolverBeg = nTracer_dclock(.FALSE., .FALSE.)

CALL cuSetCrossSection(cuCMFD, cuAxial, cuDevice)
CALL cuSetFlux(cuCMFD, cuAxial, cuDevice)
CALL cuSetPsi(cuAxial, cuDevice, 0)
n = nxy * nzSub
ierr = cublasScal(cuDevice%myblasHandle, n, inv_eigv0, cuAxial%PsiCoeff(:,:,0), 1)
CALL cuSetLeakage(cuCMFD, cuAxial, cuDevice)

!$ACC HOST_DATA USE_DEVICE(cuGeometry, cuDevice)
CALL cuLeakageSplit <<< Blocks, Threads, 0, myStream >>>                                                            &
                    (cuGeometry, cuDevice, cuAxial%phiCoeff(:, :, :, 0), cuAxial%lkgCoeff(:, :, :, 0), cuAxial%xst)
!$ACC END HOST_DATA

IF(TranCntl%lchidk) THEN
  CALL cuSetPrecSrcK(TranInfo, TranCntl)
ELSE
  CALL cuSetPrecSrc(TranInfo, TranCntl)
END IF
CALL cuSetTranFixedSrc(TranCntl)

ierr = cudaStreamSynchronize(cuDevice%myStream)
AxNSolverEnd = nTracer_dclock(.FALSE., .FALSE.)
TimeChk%AxNSolverTime = TimeChk%AxNSolverTime + (AxNSolverEnd - AxNSolverBeg)


DO iter = 1, 10
  AxNSolverBeg = nTracer_dclock(.FALSE., .FALSE.)
  IF(iter .NE. 1) THEN
    CALL cuSetPsi(cuAxial, cuDevice, 0)
    ierr = cublasScal(cuDevice%myblasHandle, n, inv_eigv0, cuAxial%PsiCoeff(:,:,0), 1)
  END IF
  CALL cuSetTranSource(TranCntl, 0)
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

SUBROUTINE cuSetTranSourceOperator(CoreInfo, FmInfo, GroupInfo, TranInfo, TranCntl, cuCMFD, cuAxial, cuDevice)
USE TYPEDEF,        ONLY : CoreInfo_Type,       FmInfo_Type,        GroupInfo_Type,                                 &
                           TranInfo_Type,       FxrInfo_Type,        Pin_Type,                                      &
                           TranCntl_Type
USE MacXsLib_Mod,   ONLY : MacP1XsScatMatrix
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice

TYPE(superPin_Type), POINTER :: superPin(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
INTEGER :: ig, igs, igb, ige, ifxr, ixy, ixy_map, ipin, ipin_map, iz, ierr, tid, iprec
INTEGER :: n, ng, nxy, nprec
INTEGER :: myzb, myze
INTEGER, POINTER :: pinMap(:)
REAL(GPU_NODAL_PRECISION), ALLOCATABLE :: S0(:, :, :, :), S1(:, :, :, :), F(:, :, :), Chi(:, :, :)
REAL(GPU_NODAL_PRECISION), ALLOCATABLE :: InvVel(:, :, :), Omegalm(:, :), Omegal0(:, :)
REAL(GPU_NODAL_PRECISION), ALLOCATABLE :: OmegalmK(:, :, :), Omegal0K(:, :, :)

Pin => CoreInfo%Pin
Fxr => FmInfo%Fxr
superPin => cuGeometry%superPin
ng = cuCMFD%ng
nxy = cuGeometry%nxyc
nprec = cuGeometry%nprec
myzb = cuDevice%myzb
myze = cuDevice%myze
pinMap => cuGeometry%pinMap

ALLOCATE(S0(nxy, ng, myzb : myze, ng)); S0 = 0.0
ALLOCATE(S1(nxy, ng, myzb : myze, ng)); S1 = 0.0
ALLOCATE(F(nxy, myzb : myze, ng))
ALLOCATE(Chi(nxy, ng, myzb : myze))
ALLOCATE(InvVel(nxy, ng, myzb : myze))
IF(TranCntl%lchidk) THEN
  ALLOCATE(OmegalmK(nxy, myzb : myze, nprec))
  ALLOCATE(Omegal0K(nxy, myzb : myze, nprec))
ELSE
  ALLOCATE(Omegalm(nxy, myzb : myze))
  ALLOCATE(Omegal0(nxy, myzb : myze))
END IF

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
      IF(TranCntl%lchidk) THEN
        Chi(ipin, ig, iz) = cuCMFD%PinXS(ipin_map, iz)%Chi(ig)
        DO iprec = 1, nprec
          Chi(ipin, ig, iz) = Chi(ipin, ig, iz) + TranInfo%chidk(ig, iprec) * (cuTranCMInfo%CellOmegap(iprec, ipin_map, iz) * TranInfo%lambda(iprec) - cuCMFD%PinXS(ipin_map, iz)%beta(iprec))
        END DO
      ELSE
        Chi(ipin, ig, iz) = cuCMFD%PinXS(ipin_map, iz)%Chi(ig) + &
          TranInfo%chid(ig) * (cuCMFD%PinXS(ipin_map, iz)%omega - cuCMFD%PinXS(ipin_map, iz)%betat)
      END IF
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
      InvVel(ipin, ig, iz) = 1./cuCMFD%PinXS(ipin_map, iz)%velo(ig)
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
ierr = cudaMemcpy(cuAxial%InvVel, InvVel, n, cudaMemcpyHostToDevice)

IF(TranCntl%lchidk) THEN
  DO iprec = 1, nprec
    !$OMP PARALLEL PRIVATE(ipin_map)
    !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
    DO iz = myzb, myze
      DO ipin = 1, nxy
        ipin_map = pinMap(ipin)
        OmegalmK(ipin, iz, iprec) = cuTranCMInfo%CellOmegam(iprec, ipin_map, iz) * TranInfo%lambda(iprec)
        Omegal0K(ipin, iz, iprec) = cuTranCMInfo%CellOmega0(iprec, ipin_map, iz) * TranInfo%lambda(iprec)
      ENDDO
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL
  END DO
  n = nxy * (myze - myzb + 1) * nprec
  ierr = cudaMemcpy(cuAxial%OmegalmK, OmegalmK, n, cudaMemcpyHostToDevice)
  ierr = cudaMemcpy(cuAxial%Omegal0K, Omegal0K, n, cudaMemcpyHostToDevice)
ELSE
  !$OMP PARALLEL PRIVATE(ipin_map)
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
  DO iz = myzb, myze
    DO ipin = 1, nxy
      ipin_map = pinMap(ipin)
      Omegalm(ipin, iz) = cuTranCMInfo%CellOmegam(0, ipin_map, iz)
      Omegal0(ipin, iz) = cuTranCMInfo%CellOmega0(0, ipin_map, iz)
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
  n = nxy * (myze - myzb + 1)
  ierr = cudaMemcpy(cuAxial%Omegalm, Omegalm, n, cudaMemcpyHostToDevice)
  ierr = cudaMemcpy(cuAxial%Omegal0, Omegal0, n, cudaMemcpyHostToDevice)
END IF

DEALLOCATE(S0, S1, F, Chi)
DEALLOCATE(InvVel)
IF(TranCntl%lchidk) THEN
  DEALLOCATE(OmegalmK)
  DEALLOCATE(Omegal0K)
ELSE
  DEALLOCATE(Omegalm)
  DEALLOCATE(Omegal0)
END IF

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
REAL, ALLOCATABLE :: phicd(:, :, :)

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzSub = cuDevice%nzSub
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
n = ng * nxy * nzSub

ALLOCATE(phis(nxy, ng, nzSub), phicd(ng, nxy, myzbf:myzef))

ierr = cudaMemCpy(phicd, cuAxial%phic, nxy*ng*(myzef-myzbf+1), cudaMemcpyDeviceToHost)
ierr = cudaMemCpy(phis, cuAxial%phiCoeff, n, cudaMemcpyDeviceToHost)

DO iz = myzbf, myzef
  !$OMP PARALLEL
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
  DO izf = cuDevice%fmSubrange(iz, 1), cuDevice%fmSubrange(iz, 2)
    DO ipin = 1, nxy
      DO ig = 1, ng
        phis(ipin, ig, izf) = phis(ipin, ig, izf) * cuCMFD%h_phis8(ig, ipin, iz)/phicd(ig, ipin, iz)
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDDO

ierr = cudaMemcpy(cuAxial%phiCoeff(:, :, :, 0), phis, n, cudaMemcpyHostToDevice)

DEALLOCATE(phis, phicd)

!CALL cuVectorOp('*', n, cuAxial%phiShape, cuAxial%phiCoeff(:, :, :, 0), cuAxial%phiCoeff(:, :, :, 0), cuDevice%myStream)

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

SUBROUTINE cuSetTranSource(TranCntl, order)
USE TYPEDEF,          ONLY : TranCntl_Type
IMPLICIT NONE
TYPE(TranCntl_Type) :: TranCntl

TYPE(dim3) :: Blocks, Threads
REAL(GPU_NODAL_PRECISION) :: delt
INTEGER(KIND = cuda_stream_kind) :: myStream
INTEGER :: ng, nxy, nzSub
INTEGER :: order

INTEGER :: ierr

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzSub = cuDevice%nzSub
myStream = cuDevice%myStream

delt = TranCntl%delt(TranCntl%nowstep)

Threads = dim3(cuDevice%cuThreadPerBlock, 1, 1)
Blocks = dim3(ng * nxy / Threads%x + 1, nzSub, 1)

!$ACC HOST_DATA USE_DEVICE(cuDevice)
CALL cuSetTranSourceKernel <<< Blocks, Threads, 0, myStream >>>                                                         &
                          (cuDevice, cuAxial%S0, cuAxial%S1, cuAxial%Chia, cuAxial%InvVel,  cuAxial%Expo_Alpha,         &
                           cuAxial%phiCoeff, cuAxial%phim, cuAxial%psiCoeff, cuAxial%lkgCoeff,                          &
                           cuAxial%FixedSrc, cuAxial%xst, cuAxial%srcCoeff, cuAxial%srcm, delt, order)
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

ATTRIBUTES(GLOBAL) SUBROUTINE cuSetTranSourceKernel(cuDevice, S0, S1, Chi, InvVel, Expo_Alpha, phis, phim, psi, &
                                                    lkg, FixedSrc , xst, src, srcm, delt, order)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuDevice_Type), INTENT(IN) :: cuDevice
REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: S0(:, :, cuDevice%myzb :, :), S1(:, :, cuDevice%myzb :, :), Chi(:, :, cuDevice%myzb :)
REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: InvVel(:, :, cuDevice%myzb :), Expo_Alpha(:, :, cuDevice%myzb :)
REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: phis(:, :, :, 0 :), phim(:, :, :), psi(:, :, 0 :), lkg(:, :, :, 0 :), xst(:, :, :)
REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: FixedSrc(:, :, :)
REAL(GPU_NODAL_PRECISION), DEVICE :: src(:, :, :, 0 :), srcm(:, :, :)
REAL(GPU_NODAL_PRECISION), VALUE :: delt
INTEGER, VALUE :: order

REAL(GPU_NODAL_PRECISION) :: scatSrc, fisSrc, momentSrc, rvdt, rvalpha
REAL(GPU_NODAL_PRECISION) :: invtheta
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

invtheta = 1./theta
rvdt = InvVel(ipin, ig, iz) * invtheta / delt
rvalpha = InvVel(ipin, ig, iz) * Expo_Alpha(ipin, ig, iz)

src(ipin, ig, izf, order) = scatSrc + fisSrc - (rvdt + rvalpha) * phis(ipin, ig, izf, order)

IF (order .EQ. 0) THEN
  momentSrc = 0.0
  DO igs = gb, ge
    momentSrc = momentSrc + S1(ipin, ig, iz, igs) * phim(ipin, igs, izf)
  ENDDO
  srcm(ipin, ig, izf) = momentSrc
ENDIF

IF (order .EQ. 0) THEN
  src(ipin, ig, izf, 0) = src(ipin, ig, izf, 0) + FixedSrc(ipin, ig, izf)
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
!LOGICAL :: lNegative

ipin = mod((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1, nxyc) + 1
ig = ((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1) / nxyc + 1

IF (ig .GT. ng) RETURN

DO iz = cuDevice%myzbf, cuDevice%myzef

!  lNegative = .FALSE.

  phisum = 0.0
  DO izf = cuDevice%fmSubrange(iz, 1), cuDevice%fmSubrange(iz, 2)
    phisum = phisum + phis(ipin, ig, izf)
!    IF (phis(ipin, ig, izf) .LT. 0.0) lNegative = .TRUE.
  ENDDO

  phisum = phisum / cuDevice%nDiv(iz)
  phic(ig, ipin, iz) = phisum

!  IF (lNegative) THEN
!    DO izf = cuDevice%fmSubrange(iz, 1), cuDevice%fmSubrange(iz, 2)
!      phiShape(ipin, ig, izf) = 1.0
!    ENDDO
!  ELSE
!    DO izf = cuDevice%fmSubrange(iz, 1), cuDevice%fmSubrange(iz, 2)
!      phiShape(ipin, ig, izf) = phis(ipin, ig, izf) / phisum
!    ENDDO
!  ENDIF

ENDDO

END SUBROUTINE

SUBROUTINE cuSetPrecSrc(TranInfo, TranCntl)
USE TYPEDEF,          ONLY : TranInfo_Type,     TranCntl_Type
IMPLICIT NONE
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl

TYPE(dim3) :: Blocks, Threads
REAL(GPU_NODAL_PRECISION), ALLOCATABLE :: lkappa(:)
REAL :: delt
INTEGER :: nowstep
INTEGER(KIND = cuda_stream_kind) :: myStream
INTEGER :: nxy, nzSub, nprec
INTEGER :: iprec
INTEGER :: n, ierr


nowstep = TranCntl%nowstep
delt = TranCntl%delt(nowstep)

nprec = cuGeometry%nprec
nxy = cuGeometry%nxyc
nzSub = cuDevice%nzSub
myStream = cuDevice%myStream

Threads = dim3(cuDevice%cuThreadPerBlock, 1, 1)
Blocks = dim3(nxy * nzSub / Threads%x + 1, 1, 1)

ALLOCATE(lkappa(nprec))

DO iprec = 1, nprec
  lkappa(iprec) = exp(-delt * TranInfo%lambda(iprec)) * TranInfo%lambda(iprec)
END DO

n = nxy * nzSub
CALL cuInitArray(n, cuAxial%PrecSrc, cuDevice%myStream)
DO iprec = 1, nprec
  ierr = cublasAxpy(cuDevice%myblasHandle, n, lkappa(iprec), cuAxial%Prec(:, :, iprec), 1, cuAxial%PrecSrc, 1)
END DO


!$ACC HOST_DATA USE_DEVICE(cuDevice)
CALL cuSetPrecSrcKernel <<< Blocks, Threads, 0, myStream >>>                                                            &
                        (cuDevice, cuAxial%Omegalm, cuAxial%TranPsid, cuAxial%PrecSrc)
CALL cuSetPrecSrcKernel <<< Blocks, Threads, 0, myStream >>>                                                            &
                        (cuDevice, cuAxial%Omegal0, cuAxial%TranPsi, cuAxial%PrecSrc)
!$ACC END HOST_DATA

DEALLOCATE(lkappa)

END SUBROUTINE

SUBROUTINE cuSetPrecSrcK(TranInfo, TranCntl)
USE TYPEDEF,          ONLY : TranInfo_Type,     TranCntl_Type
IMPLICIT NONE
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl

TYPE(dim3) :: Blocks, Threads
REAL(GPU_NODAL_PRECISION), ALLOCATABLE :: lkappa(:)
REAL :: delt
INTEGER :: nowstep
INTEGER(KIND = cuda_stream_kind) :: myStream
INTEGER :: nxy, nzSub, nprec
INTEGER :: iprec
INTEGER :: n, ierr


nowstep = TranCntl%nowstep
delt = TranCntl%delt(nowstep)

nprec = cuGeometry%nprec
nxy = cuGeometry%nxyc
nzSub = cuDevice%nzSub
myStream = cuDevice%myStream

Threads = dim3(cuDevice%cuThreadPerBlock, 1, 1)
Blocks = dim3(nxy * nzSub / Threads%x + 1, 1, 1)

ALLOCATE(lkappa(nprec))

DO iprec = 1, nprec
  lkappa(iprec) = exp(-delt * TranInfo%lambda(iprec)) * TranInfo%lambda(iprec)
END DO

n = nxy * nzSub * nprec
CALL cuInitArray(n, cuAxial%PrecSrcK, cuDevice%myStream)
n = nxy * nzSub
DO iprec = 1, nprec
  ierr = cublasAxpy(cuDevice%myblasHandle, n, lkappa(iprec), cuAxial%Prec(:, :, iprec), 1, cuAxial%PrecSrcK(:,:, iprec), 1)
END DO

!$ACC HOST_DATA USE_DEVICE(cuDevice)
DO iprec = 1, nprec
  CALL cuSetPrecSrcKernel <<< Blocks, Threads, 0, myStream >>>                                                            &
    (cuDevice, cuAxial%OmegalmK(:,:,iprec), cuAxial%TranPsid, cuAxial%PrecSrcK(:,:,iprec))
  CALL cuSetPrecSrcKernel <<< Blocks, Threads, 0, myStream >>>                                                            &
    (cuDevice, cuAxial%Omegal0K(:,:,iprec), cuAxial%TranPsi, cuAxial%PrecSrcK(:,:,iprec))
END DO
!$ACC END HOST_DATA

DEALLOCATE(lkappa)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuSetPrecSrcKernel(cuDevice, Omega, TranPsi, PrecSrc)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuDevice_Type), INTENT(IN) :: cuDevice
REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: Omega(:, cuDevice%myzb :)
REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: TranPsi(:, :)
REAL(GPU_NODAL_PRECISION), DEVICE :: PrecSrc(:, :)

INTEGER :: ipin, iz, izf

ipin = mod((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1, nxyc) + 1
izf = ((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1) / nxyc + 1

IF (izf .GT. cuDevice%nzSub) RETURN

iz = cuDevice%cmMap(izf)
PrecSrc(ipin, izf) = PrecSrc(ipin, izf) + Omega(ipin, iz) * TranPsi(ipin, izf)

END SUBROUTINE

SUBROUTINE cuSetTranFixedSrc(TranCntl)
USE TYPEDEF,          ONLY : TranCntl_Type
IMPLICIT NONE
TYPE(TranCntl_Type) :: TranCntl

TYPE(dim3) :: Blocks, Threads
REAL(GPU_NODAL_PRECISION) :: delt
INTEGER(KIND = cuda_stream_kind) :: myStream
INTEGER :: ng, nxy, nzSub

INTEGER :: ierr

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzSub = cuDevice%nzSub
myStream = cuDevice%myStream

delt = TranCntl%delt(TranCntl%nowstep)

Threads = dim3(cuDevice%cuThreadPerBlock, 1, 1)
Blocks = dim3(ng * nxy / Threads%x + 1, nzSub, 1)

IF(TranCntl%lchidk) THEN
!$ACC HOST_DATA USE_DEVICE(cuDevice)
CALL cuSetTranFixedSrcKernel_chidk <<< Blocks, Threads, 0, myStream >>>                                     &
                             (cuDevice, cuAxial%FixedSrc, cuAxial%PrecSrcK, cuAxial%ResSrc,            &
                              cuAxial%TranPhi, cuAxial%InvVel, cuAxial%Expo, cuAxial%Expo_Alpha, delt)
!$ACC END HOST_DATA
ELSE
!$ACC HOST_DATA USE_DEVICE(cuDevice)
CALL cuSetTranFixedSrcKernel <<< Blocks, Threads, 0, myStream >>>                                     &
                             (cuDevice, cuAxial%FixedSrc, cuAxial%PrecSrc, cuAxial%ResSrc,            &
                              cuAxial%TranPhi, cuAxial%InvVel, cuAxial%Expo, cuAxial%Expo_Alpha, delt)
!$ACC END HOST_DATA
END IF


END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuSetTranFixedSrcKernel(cuDevice, FixedSrc, PrecSrc, ResSrc,  &
                                                      TranPhi, InvVel, expo, expo_alpha, delt)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuDevice_Type), INTENT(IN) :: cuDevice
REAL(GPU_NODAL_PRECISION), DEVICE :: FixedSrc(:, :, :)
REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: PrecSrc(:, :), ResSrc(:, :, :), TranPhi(:, :, :)
REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: InvVel(:, :, cuDevice%myzb :), expo(:, :, cuDevice%myzb:), expo_alpha(:, :, cuDevice%myzb:)
REAL(GPU_NODAL_PRECISION), VALUE :: delt

REAL(GPU_NODAL_PRECISION) :: invtheta, thetah, prevSrc
INTEGER :: ipin, ig, izf, iz

ipin = mod((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1, nxyc) + 1
ig = ((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1) / nxyc + 1
izf = blockIdx%y

IF (ig .GT. ng) RETURN

iz = cuDevice%cmMap(izf)

invtheta = 1./ theta
thetah = invtheta - 1.

prevSrc = thetah * (ResSrc(ipin, ig, izf) - Expo_Alpha(ipin, ig, iz) * InvVel(ipin, ig, iz) * TranPhi(ipin, ig, izf))
prevSrc = prevSrc + InvVel(ipin, ig, iz) * TranPhi(ipin, ig, izf) * invtheta / delt
prevSrc = prevSrc * Expo(ipin, ig, iz)
FixedSrc(ipin, ig, izf) = prevSrc + chid(ig) * PrecSrc(ipin, izf)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuSetTranFixedSrcKernel_chidk(cuDevice, FixedSrc, PrecSrcK, ResSrc,  &
                                                      TranPhi, InvVel, expo, expo_alpha, delt)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuDevice_Type), INTENT(IN) :: cuDevice
REAL(GPU_NODAL_PRECISION), DEVICE :: FixedSrc(:, :, :)
REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: PrecSrcK(:, :, :), ResSrc(:, :, :), TranPhi(:, :, :)
REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: InvVel(:, :, cuDevice%myzb :), expo(:, :, cuDevice%myzb:), expo_alpha(:, :, cuDevice%myzb:)
REAL(GPU_NODAL_PRECISION), VALUE :: delt

REAL(GPU_NODAL_PRECISION) :: invtheta, thetah, prevSrc
INTEGER :: ipin, ig, izf, iz, iprec

ipin = mod((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1, nxyc) + 1
ig = ((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1) / nxyc + 1
izf = blockIdx%y

IF (ig .GT. ng) RETURN

iz = cuDevice%cmMap(izf)

invtheta = 1./ theta
thetah = invtheta - 1.

prevSrc = thetah * (ResSrc(ipin, ig, izf) - Expo_Alpha(ipin, ig, iz) * InvVel(ipin, ig, iz) * TranPhi(ipin, ig, izf))
prevSrc = prevSrc + InvVel(ipin, ig, iz) * TranPhi(ipin, ig, izf) * invtheta / delt
prevSrc = prevSrc * Expo(ipin, ig, iz)
DO iprec = 1, nprec
  FixedSrc(ipin, ig, izf) = prevSrc + chidK(ig, iprec) * PrecSrcK(ipin, izf, iprec)
END DO

END SUBROUTINE

END MODULE

#endif
