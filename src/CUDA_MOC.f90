#include <CUDADEFINES.h>

#ifdef __PGI

MODULE CUDA_MOC

USE CUDA_MASTER
USE CUDA_UTIL
IMPLICIT NONE

PRIVATE
PUBLIC :: CUDASourceUpdate, CUDARayTrace, CUDARayTraceAsync, CUDAP1RayTrace, CUDAInitPhiAngIn

CONTAINS

!--- Host Subroutine ------------------------------------------------------------------------------

SUBROUTINE CUDAInitPhiAngIn(ngtot, nsv)
USE CUDAFOR
IMPLICIT NONE
INTEGER :: ngtot, nsv

INTEGER :: ydim
TYPE(dim3) :: Blocks, Threads
INTEGER(KIND = cuda_stream_kind) :: myStream

myStream = cuDevice%myStream

ydim = cuDevice%cuMaxThreadPerBlock/P0_BLOCK_SIZE
Threads = dim3(P0_BLOCK_SIZE, ydim, 1)
Blocks = dim3(ngtot/P0_BLOCK_SIZE + 1, nsv/ydim + 1, 1)

CALL cuInitializePhiAngIn<<<Blocks,Threads,0,myStream>>>(cuMOC%PhiAngIn, ngtot)

END SUBROUTINE
SUBROUTINE CUDASourceUpdate(iz, gb, ge, eigv)
USE CUDAFOR
IMPLICIT NONE

INTEGER :: iz, gb, ge
REAL :: eigv
TYPE(dim3) :: Blocks, Threads
INTEGER(KIND = cuda_stream_kind) :: myStream
INTEGER :: ng

ng = ge - gb + 1
myStream = cuDevice%myStream

Threads = dim3(ng, cuDevice%cuMaxThreadPerBlock / ng, 1)
Blocks = dim3(cuGeometry%nfsr / Threads%y + 1, 1, 1)

!$ACC HOST_DATA USE_DEVICE(cuGeometry)
CALL cuComputeSource <<< Blocks, Threads, 0, myStream >>>                                                           &
                     (cuGeometry, cuMOC%phis, cuMOC%psi, cuMOC%src, cuMOC%xst, cuMOC%xssm, cuMOC%chi,               &
                      iz, gb, ge, eigv)
!$ACC END HOST_DATA

END SUBROUTINE

SUBROUTINE CUDARayTrace(cuMOC, cuDevice, jout, lJout, iz, gb, ge, time)
USE PARAM
USE CUDAFOR
IMPLICIT NONE

TYPE(cuMOC_Type) :: cuMOC
TYPE(cuDevice_Type) :: cuDevice
REAL, POINTER :: jout(:, :, :, :)
LOGICAL :: lJout
INTEGER :: iz, gb, ge
REAL, OPTIONAL :: time

TYPE(dim3) :: Blocks, Threads
TYPE(cudaEvent) :: startEvent, stopEvent
INTEGER(KIND = cuda_stream_kind) :: myStream
INTEGER :: copyMemSize, sharedMemSize
INTEGER :: iAx, ierr
INTEGER :: ng
REAL(4) :: elapsedTime

ierr = cudaEventCreate(startEvent)
ierr = cudaEventCreate(stopEvent)

ng = ge - gb + 1
myStream = cuDevice%myStream

cuMOC%phis(gb : ge, :) = zero
IF (lJout) cuMOC%jout(:, gb : ge, :, :) = zero

ierr = cudaEventRecord(startEvent, myStream)

sharedMemSize = GPU_PRECISION * cuDevice%sharedMemoryDim * cuGeometry%nPolarAngle

IF (cuDevice%lFullWarp) THEN
  Threads = dim3(cuDevice%cuThreadPerBlock, 1, 1)
  Blocks = dim3(cuDevice%nRotRay * ng / Threads%x + 1, 2, 1)
ELSE
  Threads = dim3(ng, 2, 1)
  Blocks = dim3(cuDevice%nRotRay, 1, 1)
ENDIF

iAx = cuGeometry%AxType(iz)

!$ACC HOST_DATA USE_DEVICE(cuFastRay1D, cuDevice)
CALL cu1DFastRayTrace <<< Blocks, Threads, sharedMemSize, myStream >>>                                              &
                      (cuFastRay1D(iAx), cuDevice, cuMOC%phis, cuMOC%xst, cuMOC%src, cuMOC%jout, cuMOC%PhiAngIn,    &
                       lJout, iz, gb, ge)
!$ACC END HOST_DATA
                       
Threads = dim3(ng, cuDevice%cuMaxThreadPerBlock / ng, 1)
Blocks = dim3(cuGeometry%nfsr / Threads%y + 1, 1, 1)

!$ACC HOST_DATA USE_DEVICE(cuGeometry)
CALL cuComputeScalarFlux <<< Blocks, Threads, 0, myStream >>>                                                       &
                         (cuGeometry, cuMOC%phis, cuMOC%src, cuMOC%xst, iz, gb, ge)
!$ACC END HOST_DATA

ierr = cudaEventRecord(stopEvent, myStream)
ierr = cudaEventSynchronize(stopEvent)

IF (lJout) jout = cuMOC%jout
  
ierr = cudaEventElapsedTime(elapsedTime, startEvent, stopEvent)
ierr = cudaEventDestroy(startEvent)
ierr = cudaEventDestroy(stopEvent)

IF (PRESENT(time)) time = time + elapsedTime / 1000.0

END SUBROUTINE

SUBROUTINE CUDARayTraceAsync(cuMOC, cuDevice, phis, src, xst, jout, PhiAngIn, lJout, iz, ng)
USE PARAM
USE CUDAFOR
IMPLICIT NONE

TYPE(cuMOC_Type) :: cuMOC
TYPE(cuDevice_Type) :: cuDevice
REAL(GPU_PRECISION), PINNED :: PhiAngIn(:, :, :)
REAL(GPU_FLUX_PRECISION), PINNED :: phis(:, :), jout(:, :, :, :)
REAL(GPU_SOURCE_PRECISION), PINNED :: src(:, :), xst(:, :)
LOGICAL :: lJout
INTEGER :: iz, ng

TYPE(dim3) :: Blocks, Threads
INTEGER(KIND = cuda_stream_kind) :: myStream
INTEGER :: copyMemSize, sharedMemSize
INTEGER :: iAx, ierr

myStream = cuDevice%myStream

CALL cuInitArray(cuGeometry%nFsr * P0_BLOCK_SIZE, cuMOC%phisMg, myStream)
IF (lJout) CALL cuInitArray(cuGeometry%nxy * P0_BLOCK_SIZE * 3 * 4, cuMOC%joutMg, myStream)

copyMemSize = P0_BLOCK_SIZE * cuGeometry%nPolarAngle * cuGeometry%nPhiAngSv
ierr = cudaMemcpyAsync(cuMOC%PhiAngInMg, PhiAngIn, copyMemSize, cudaMemcpyHostToDevice, myStream)

copyMemSize = P0_BLOCK_SIZE * cuGeometry%nFsr
ierr = cudaMemcpyAsync(cuMOC%xstMg, xst, copyMemSize, cudaMemcpyHostToDevice, myStream)

copyMemSize = P0_BLOCK_SIZE * cuGeometry%nFsr
ierr = cudaMemcpyAsync(cuMOC%srcMg, src, copyMemSize, cudaMemcpyHostToDevice, myStream)

sharedMemSize = GPU_PRECISION * cuDevice%sharedMemoryDim * cuGeometry%nPolarAngle

IF (cuDevice%lFullWarp) THEN
  Threads = dim3(cuDevice%cuThreadPerBlock, 1, 1)
  Blocks = dim3(cuDevice%nRotRay * ng / Threads%x + 1, 2, 1)
ELSE
  Threads = dim3(ng, 2, 1)
  Blocks = dim3(cuDevice%nRotRay, 1, 1)
ENDIF

iAx = cuGeometry%AxType(iz)

!$ACC HOST_DATA USE_DEVICE(cuFastRay1D, cuDevice)
CALL cu1DFastRayTraceAsync <<< Blocks, Threads, sharedMemSize, myStream >>>                                         &
                           (cuFastRay1D(iAx), cuDevice, cuMOC%phisMg, cuMOC%xstMg, cuMOC%srcMg, cuMOC%joutMg,       &
                            cuMOC%PhiAngInMg, lJout, iz, ng)
!$ACC END HOST_DATA

Threads = dim3(ng, cuDevice%cuMaxThreadPerBlock / ng, 1)
Blocks = dim3(cuGeometry%nfsr / Threads%y + 1, 1, 1)

!$ACC HOST_DATA USE_DEVICE(cuGeometry)
CALL cuComputeScalarFluxAsync <<< Blocks, Threads, 0, myStream >>>                                                  &
                              (cuGeometry, cuMOC%phisMg, cuMOC%srcMg, cuMOC%xstMg, iz)
!$ACC END HOST_DATA

copyMemSize = P0_BLOCK_SIZE * cuGeometry%nfsr
ierr = cudaMemcpyAsync(phis, cuMOC%phisMg, copyMemSize, cudaMemcpyDeviceToHost, myStream)

copyMemSize = P0_BLOCK_SIZE * cuGeometry%nPolarAngle * cuGeometry%nPhiAngSv
ierr = cudaMemcpyAsync(PhiAngIn, cuMOC%PhiAngInMg, copyMemSize, cudaMemcpyDeviceToHost, myStream)

IF (lJout) THEN
  copyMemSize = P0_BLOCK_SIZE * 3 * 4 * cuGeometry%nxy
  ierr = cudaMemcpyAsync(jout, cuMOC%joutMg, copyMemSize, cudaMemcpyDeviceToHost, myStream)
ENDIF

END SUBROUTINE

SUBROUTINE CUDAP1RayTrace(cuMOC, cuDevice, phis, phim, xst, src, srcm, jout, PhiAngIn, lJout, iz, ng)
USE PARAM
USE CNTL,			ONLY : nTracerCntl
USE CUDAFOR
USE OPENACC
USE OMP_LIB
USE CUDA_CONST, ONLY : nAziMap
IMPLICIT NONE

TYPE(cuMOC_Type) :: cuMOC
TYPE(cuDevice_Type) :: cuDevice
REAL(GPU_PRECISION), PINNED :: PhiAngIn(:, :, :)
REAL(GPU_FLUX_PRECISION), PINNED :: phis(:, :), phim(:, :, :), jout(:, :, :, :)
REAL(GPU_SOURCE_PRECISION), PINNED :: xst(:, :), src(:, :), srcm(:, :, :)
LOGICAL :: lJout
INTEGER :: iz, ng

TYPE(dim3) :: Blocks, Threads
INTEGER(KIND = cuda_stream_kind) :: myStream
INTEGER :: copyMemSize
INTEGER :: i, j, k, l, iAx, iAzi, ierr
INTEGER :: nMoment, ScatOd, nAzi
LOGICAL :: lRot

ScatOd = nTracerCntl%ScatOd
myStream = cuDevice%myStream

IF (ScatOd .EQ. 1) nMoment = 2
IF (ScatOd .EQ. 2) nMoment = 5
IF (ScatOd .EQ. 3) nMoment = 9

lRot = cuGeometry%lRot

CALL cuInitArray(cuGeometry%nFsr * PN_BLOCK_SIZE, cuMOC%phisMg, myStream)
CALL cuInitArray(cuGeometry%nFsr * PN_BLOCK_SIZE * nMoment, cuMOC%phimMg, myStream)
IF (lJout) CALL cuInitArray(cuGeometry%nxy * PN_BLOCK_SIZE * 3 * 4, cuMOC%joutMg, myStream)

!--- Data Copy to Device

copyMemSize = PN_BLOCK_SIZE * cuGeometry%nPolarAngle * cuGeometry%nPhiAngSv
ierr = cudaMemcpyAsync(cuMOC%PhiAngInMg, PhiAngIn, copyMemSize, cudaMemcpyHostToDevice, myStream)

copyMemSize = PN_BLOCK_SIZE * cuGeometry%nfsr
ierr = cudaMemcpyAsync(cuMOC%xstMg, xst, copyMemSize, cudaMemcpyHostToDevice, myStream)

copyMemSize = PN_BLOCK_SIZE * cuGeometry%nfsr
ierr = cudaMemcpyAsync(cuMOC%srcMg, src, copyMemSize, cudaMemcpyHostToDevice, myStream)

copyMemSize = PN_BLOCK_SIZE * nMoment * cuGeometry%nfsr
ierr = cudaMemcpyAsync(cuMOC%srcmMg, srcm, copyMemSize, cudaMemcpyHostToDevice, myStream)


nAzi = cuGeometry%nAziAngle / 2
IF(cuGeometry%lRot) nAzi = cuGeometry%nAziAngle / 4
DO iAzi = 1, nAzi

  !--- Prepare Angle Dependent Source

  Threads = dim3(cuGeometry%nPolarAngle, ng, cuDevice%cuThreadPerBlock / cuGeometry%nPolarAngle / ng)
  Blocks = dim3(cuGeometry%nfsr / Threads%z + 1, 1, 1)

  !$ACC HOST_DATA USE_DEVICE(cuDevice)
  CALL cuPrepareAngularSource <<< Blocks, Threads, 0, myStream >>>                                                  &
                              (cuDevice, cuMOC%srcMg, cuMOC%srcmMg, cuMOC%SrcAngMg, iAzi, lRot)
  !$ACC END HOST_DATA

  !--- Initialize Angular Flux

  CALL cuInitArray(cuGeometry%nFsr * cuGeometry%nPolarAngle * PN_BLOCK_SIZE * nAziMap, cuMOC%phiaMg, myStream)

  !--- Ray Tracing with Angular Flux Saving

  Threads = dim3(cuGeometry%nPolarAngle, ng, cuDevice%cuThreadPerBlock / cuGeometry%nPolarAngle / ng)
  Blocks = dim3(cuGeometry%RotRayCount(iAzi) / Threads%z + 1, 2, 1)
  
  iAx = cuGeometry%AxType(iz)
  
  !$ACC HOST_DATA USE_DEVICE(cuFastRay1D, cuGeometry, cuDevice)
  CALL cu1DFastRayTraceP1 <<< Blocks, Threads, 0, myStream >>>                                                      &
                          (cuFastRay1D(iAx), cuGeometry, cuDevice, cuMOC%phiaMg, cuMOC%xstMg, cuMOC%SrcAngMg,       &
                           cuMOC%joutMg, cuMOC%PhiAngInMg, lJout, iAzi)
  !$ACC END HOST_DATA

  !--- Weight Sum Angular Flux

  Threads = dim3(ng, cuDevice%cuThreadPerBlock / ng, 1)
  Blocks = dim3(cuGeometry%nfsr / Threads%y + 1, 1, 1)

  !$ACC HOST_DATA USE_DEVICE(cuDevice)
  CALL cuIncrementPseudoMoment <<< Blocks, Threads, 0, myStream >>>                                                 &
                               (cuDevice, cuMOC%phisMg, cuMOC%phimMg, cuMOC%phiaMg, iAzi, lRot)
  !$ACC END HOST_DATA

ENDDO

Threads = dim3(ng, cuDevice%cuThreadPerBlock / ng, 1)
Blocks = dim3(cuGeometry%nfsr / Threads%y + 1, 1, 1)

!$ACC HOST_DATA USE_DEVICE(cuGeometry)
CALL cuComputeFluxMoment <<< Blocks, Threads, 0, myStream >>>                                                       &
                         (cuGeometry, cuMOC%phisMg, cuMOC%phimMg, cuMOC%srcMg, cuMOC%srcmMg, cuMOC%xstMg, iz)
!$ACC END HOST_DATA

!--- Data Copy to Host

copyMemSize = PN_BLOCK_SIZE * cuGeometry%nfsr
ierr = cudaMemcpyAsync(phis, cuMOC%phisMg, copyMemSize, cudaMemcpyDeviceToHost, myStream)

copyMemSize = PN_BLOCK_SIZE * nMoment * cuGeometry%nfsr
ierr = cudaMemcpyAsync(phim, cuMOC%phimMg, copyMemSize, cudaMemcpyDeviceToHost, myStream)

copyMemSize = PN_BLOCK_SIZE * cuGeometry%nPolarAngle * cuGeometry%nPhiAngSv
ierr = cudaMemcpyAsync(PhiAngIn, cuMOC%PhiAngInMg, copyMemSize, cudaMemcpyDeviceToHost, myStream)

IF (lJout) THEN
  copyMemSize = PN_BLOCK_SIZE * 3 * 4 * cuGeometry%nxy
  ierr = cudaMemcpyAsync(jout, cuMOC%joutMg, copyMemSize, cudaMemcpyDeviceToHost, myStream)
ENDIF

END SUBROUTINE

!--- Utils ----------------------------------------------------------------------------------------

ATTRIBUTES(GLOBAL) SUBROUTINE cuInitializePhiAngIn(PhiAngIn, ngtot)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE
REAL(GPU_FLUX_PRECISION), DEVICE :: PhiAngIn(:, :, :)
INTEGER, VALUE :: ngtot

INTEGER :: ig, isv, ipol

ig = threadIdx%x+blockDim%x*(blockIdx%x-1)
isv = threadIdx%y+blockDim%y*(blockIdx%y-1)
IF (ig .GT. ngtot) RETURN
IF (isv .GT. nPhiAngSv) RETURN

IF (isv.EQ.1) THEN
  DO ipol = 1, nPolarAngle
    PhiAngIn(ipol,ig,isv) = 0.
  END DO
ELSE
  DO ipol = 1, nPolarAngle
    PhiAngIn(ipol,ig,isv) = 1.
  END DO
END IF

END SUBROUTINE
!--- P0 Domain Connected --------------------------------------------------------------------------

ATTRIBUTES(GLOBAL) SUBROUTINE cuComputeScalarFlux(cuGeometry, phis, src, xst, iz, gb, ge)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuGeometry_Type), INTENT(IN) :: cuGeometry
REAL(GPU_SOURCE_PRECISION), DEVICE, INTENT(IN) :: src(:, :)
REAL(GPU_SOURCE_PRECISION), DEVICE, INTENT(IN) :: xst(:, :)
REAL(GPU_FLUX_PRECISION), DEVICE :: phis(:, :)
INTEGER, VALUE :: iz, gb, ge

INTEGER :: ig, ifsr, itype
REAL(GPU_PRECISION) :: rsigv

ig = threadIdx%x + gb - 1
ifsr = threadIdx%y + (blockIdx%x - 1) * blockDim%y
IF (ifsr .GT. nfsr) RETURN

itype = cuGeometry%AxType(iz)
rsigv = 1.0 / (xst(ig, ifsr) * cuGeometry%FsrVol(ifsr, itype))
phis(ig, ifsr) = phis(ig, ifsr) * rsigv + src(ig, ifsr)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuComputeScalarFluxAsync(cuGeometry, phis, src, xst, iz)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuGeometry_Type), INTENT(IN) :: cuGeometry
REAL(GPU_SOURCE_PRECISION), DEVICE, INTENT(IN) :: src(P0_BLOCK_SIZE, nfsr)
REAL(GPU_SOURCE_PRECISION), DEVICE, INTENT(IN) :: xst(P0_BLOCK_SIZE, nfsr)
REAL(GPU_FLUX_PRECISION), DEVICE :: phis(P0_BLOCK_SIZE, nfsr)
INTEGER, VALUE :: iz

INTEGER :: ig, ifsr, itype
REAL(GPU_PRECISION) :: rsigv

ig = threadIdx%x
ifsr = threadIdx%y + (blockIdx%x - 1) * blockDim%y
IF (ifsr .GT. nfsr) RETURN

itype = cuGeometry%AxType(iz)
rsigv = 1.0 / (xst(ig, ifsr) * cuGeometry%FsrVol(ifsr, itype))
phis(ig, ifsr) = phis(ig, ifsr) * rsigv + src(ig, ifsr)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuComputeSource(cuGeometry, phis, psi, src, xst, xssm, chi, iz, gb, ge, eigv)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuGeometry_Type), INTENT(IN) :: cuGeometry
REAL(GPU_FLUX_PRECISION), DEVICE, INTENT(IN) :: phis(:, :)
REAL(GPU_SOURCE_PRECISION), DEVICE, INTENT(IN) :: psi(:)
REAL(GPU_SOURCE_PRECISION), DEVICE, INTENT(IN) :: xst(:, :), xssm(:, :, :), chi(:, :)
REAL(GPU_SOURCE_PRECISION), DEVICE :: src(:, :)
INTEGER, VALUE :: iz, gb, ge
REAL, VALUE :: eigv

INTEGER :: ig, igf, ifsr, ifxr

ig = threadIdx%x + gb - 1
ifsr = threadIdx%y + (blockIdx%x - 1) * blockDim%y
IF (ifsr .GT. nfsr) RETURN

ifxr = cuGeometry%Fsr2Fxr(ifsr, iz)
IF (ifxr .EQ. 0) RETURN

src(ig, ifsr) = chi(ig, ifxr) * psi(ifsr) / eigv

DO igf = 1, ng
  src(ig, ifsr) = src(ig, ifsr) + xssm(ig, igf, ifxr) * phis(igf, ifsr)
ENDDO

src(ig, ifsr) = src(ig, ifsr) / xst(ig, ifsr)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cu1DFastRayTrace(cuRotRay1D, cuDevice, phis, xst, src, jout, PhiAngIn,                &
											   lJout, iz, gb, ge)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuDevice_Type), INTENT(IN) :: cuDevice
TYPE(cuRotRay1D_Type), INTENT(IN) :: cuRotRay1D
REAL(GPU_SOURCE_PRECISION), DEVICE, INTENT(IN) :: xst(:, :), src(:, :)
REAL(GPU_PRECISION), DEVICE :: PhiAngIn(:, :, :)
REAL(GPU_FLUX_PRECISION), DEVICE :: phis(:, :), jout(:, :, :, :)
LOGICAL, VALUE :: lJout
INTEGER, VALUE :: iz, gb, ge

INTEGER :: iRotRay, ig, irot, itrack

IF (cuDevice%lFullWarp) THEN
  iRotRay = ((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1) / (ge - gb + 1) + cuDevice%RotRayBeg
  ig = mod((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1, ge - gb + 1) + gb
  irot = blockIdx%y
  itrack = threadIdx%x
ELSE
  iRotRay = blockIdx%x + cuDevice%RotRayBeg - 1
  ig = threadIdx%x + gb - 1
  irot = threadIdx%y
  itrack = ig + (irot - 1) * ng
ENDIF

IF (iRotRay .GT. cuDevice%RotRayEnd) RETURN

IF (.NOT. lJout) THEN
  CALL cuTrack1DFastRay1st(cuRotRay1D, cuDevice, phis, xst, src, PhiAngIn,                                          &
                           ig, irot, itrack, iRotRay)
ELSE
  CALL cuTrack1DFastRay2nd(cuRotRay1D, cuDevice, phis, xst, src, jout, PhiAngIn,                                    &
                           ig, irot, itrack, iRotRay)
ENDIF

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cu1DFastRayTraceAsync(cuRotRay1D, cuDevice, phis, xst, src, jout, PhiAngIn,           &
											        lJout, iz, ng)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuDevice_Type), INTENT(IN) :: cuDevice
TYPE(cuRotRay1D_Type), INTENT(IN) :: cuRotRay1D
REAL(GPU_SOURCE_PRECISION), DEVICE, INTENT(IN) :: xst(P0_BLOCK_SIZE, nFsr), src(P0_BLOCK_SIZE, nFsr)
REAL(GPU_PRECISION), DEVICE :: PhiAngIn(nPolarAngle, P0_BLOCK_SIZE, nPhiAngSv)
REAL(GPU_FLUX_PRECISION), DEVICE :: phis(P0_BLOCK_SIZE, nFsr), jout(3, P0_BLOCK_SIZE, 4, nxy)
LOGICAL, VALUE :: lJout
INTEGER, VALUE :: iz, ng

INTEGER :: iRotRay, ig, irot, itrack

IF (cuDevice%lFullWarp) THEN
  iRotRay = ((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1) / ng + cuDevice%RotRayBeg
  ig = mod((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1, ng) + 1
  irot = blockIdx%y
  itrack = threadIdx%x
ELSE
  iRotRay = blockIdx%x + cuDevice%RotRayBeg - 1
  ig = threadIdx%x
  irot = threadIdx%y
  itrack = ig + (irot - 1) * ng
ENDIF

IF (iRotRay .GT. cuDevice%RotRayEnd) RETURN

IF (.NOT. lJout) THEN
  CALL cuTrack1DFastRay1st(cuRotRay1D, cuDevice, phis, xst, src, PhiAngIn,                                          &
                           ig, irot, itrack, iRotRay)
ELSE
  CALL cuTrack1DFastRay2nd(cuRotRay1D, cuDevice, phis, xst, src, jout, PhiAngIn,                                    &
                           ig, irot, itrack, iRotRay)
ENDIF

END SUBROUTINE

ATTRIBUTES(DEVICE) SUBROUTINE cuTrack1DFastRay1st(cuRotRay1D, cuDevice, phis, xst, src, PhiAngIn,                   &
                                                  ig, irot, itrack, iRotRay)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuRotRay1D_Type), INTENT(IN) :: cuRotRay1D
TYPE(cuDevice_Type), INTENT(IN) :: cuDevice
REAL(GPU_SOURCE_PRECISION), DEVICE, INTENT(IN) :: xst(:, :), src(:, :)
REAL(GPU_PRECISION), DEVICE :: PhiAngIn(:, :, :)
REAL(GPU_FLUX_PRECISION), DEVICE :: phis(:, :)
INTEGER :: ig, irot, itrack, iRotRay

INTEGER :: RaySegBeg, RaySegEnd
INTEGER :: iazi, ireg
INTEGER :: j, ir, ipol
REAL(GPU_PRECISION) :: fsr_src, delta_phi, track_phi, fsr_phi
REAL(GPU_PRECISION) :: tau, ExpApp

REAL(GPU_PRECISION), SHARED :: s_track_phi(cuDevice%sharedMemoryDim, nPolarAngle)

s_track_phi(itrack, :) = PhiAngIn(:, ig, cuRotRay1D%PhiAngInSvIdx(irot, iRotRay))

IF (irot .EQ. 1) THEN

DO j = cuRotRay1D%RotRayIdxSt(iRotRay), cuRotRay1D%RotRayIdxSt(iRotRay + 1) - 1
  iazi = cuRotRay1D%iAzi(j)
  RaySegBeg = cuRotRay1D%PinRayIdxSt(cuRotRay1D%CoreRayIdxSt(j))
  RaySegEnd = cuRotRay1D%PinRayIdxSt(cuRotRay1D%CoreRayIdxSt(j + 1)) - 1
  DO ir = RaySegBeg, RaySegEnd
    ireg = cuRotRay1D%FsrIdx(ir)
    tau = - cuRotRay1D%LenSeg(ir) * xst(ig, ireg)
	fsr_src = src(ig, ireg); fsr_phi = 0.0
	DO ipol = 1, nPolarAngle
	  ExpApp = 1.0 - __EXPF(tau * rsinv(ipol))
	  track_phi = s_track_phi(itrack, ipol)
	  delta_phi = (track_phi - fsr_src) * ExpApp
	  s_track_phi(itrack, ipol) = track_phi - delta_phi
	  fsr_phi = fsr_phi + wt(ipol, iazi) * delta_phi
	ENDDO
	!$ACC ATOMIC
	phis(ig, ireg) = phis(ig, ireg) + fsr_phi
  ENDDO
ENDDO

ELSE

DO j = cuRotRay1D%RotRayIdxSt(iRotRay + 1) - 1, cuRotRay1D%RotRayIdxSt(iRotRay), -1
  iazi = cuRotRay1D%iAzi(j)
  RaySegBeg = cuRotRay1D%PinRayIdxSt(cuRotRay1D%CoreRayIdxSt(j + 1)) - 1
  RaySegEnd = cuRotRay1D%PinRayIdxSt(cuRotRay1D%CoreRayIdxSt(j))
  DO ir = RaySegBeg, RaySegEnd, -1
    ireg = cuRotRay1D%FsrIdx(ir)
    tau = - cuRotRay1D%LenSeg(ir) * xst(ig, ireg)
	fsr_src = src(ig, ireg); fsr_phi = 0.0
	DO ipol = 1, nPolarAngle
	  ExpApp = 1.0 - __EXPF(tau * rsinv(ipol))
	  track_phi = s_track_phi(itrack, ipol)
	  delta_phi = (track_phi - fsr_src) * ExpApp
	  s_track_phi(itrack, ipol) = track_phi - delta_phi
	  fsr_phi = fsr_phi + wt(ipol, iazi) * delta_phi
	ENDDO
	!$ACC ATOMIC
	phis(ig, ireg) = phis(ig, ireg) + fsr_phi
  ENDDO
ENDDO

ENDIF

PhiAngIn(:, ig, cuRotRay1D%PhiAngOutSvIdx(irot, iRotRay)) = s_track_phi(itrack, :)

END SUBROUTINE

ATTRIBUTES(DEVICE) SUBROUTINE cuTrack1DFastRay2nd(cuRotRay1D, cuDevice, phis, xst, src, jout, PhiAngIn,             &
												  ig, irot, itrack, iRotRay)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuRotRay1D_Type), INTENT(IN) :: cuRotRay1D
TYPE(cuDevice_Type), INTENT(IN) :: cuDevice
REAL(GPU_SOURCE_PRECISION), DEVICE, INTENT(IN) :: xst(:, :), src(:, :)
REAL(GPU_PRECISION), DEVICE :: PhiAngIn(:, :, :)
REAL(GPU_FLUX_PRECISION), DEVICE :: phis(:, :), jout(:, :, :, :)
INTEGER :: ig, irot, itrack, iRotRay

INTEGER :: PinRayBeg, PinRayEnd
INTEGER :: RaySegBeg, RaySegEnd
INTEGER :: iazi, ireg, ipin, isurf(2)
INTEGER :: j, k, ir, ipol
REAL(GPU_PRECISION) :: pin_jout
REAL(GPU_PRECISION) :: fsr_src, delta_phi, track_phi, fsr_phi
REAL(GPU_PRECISION) :: tau, ExpApp

REAL(GPU_PRECISION), SHARED :: s_track_phi(cuDevice%sharedMemoryDim, nPolarAngle)

s_track_phi(itrack, :) = PhiAngIn(:, ig, cuRotRay1D%PhiAngInSvIdx(irot, iRotRay))

IF (irot .EQ. 1) THEN

DO j = cuRotRay1D%RotRayIdxSt(iRotRay), cuRotRay1D%RotRayIdxSt(iRotRay + 1) - 1
  iazi = cuRotRay1D%iAzi(j)
  PinRayBeg = cuRotRay1D%CoreRayIdxSt(j)
  PinRayEnd = cuRotRay1D%CoreRayIdxSt(j + 1) - 1
  DO k = PinRayBeg, PinRayEnd
	ipin = cuRotRay1D%PinIdx(k)
	isurf = cuRotRay1D%SurfIdx(:, k)
	pin_jout = 0.0
	DO ipol = 1, nPolarAngle
	  track_phi = s_track_phi(itrack, ipol)
	  pin_jout = pin_jout + wt(ipol, iazi) * track_phi
	ENDDO
	!$ACC ATOMIC
    jout(1, ig, isurf(1), ipin) = jout(1, ig, isurf(1), ipin) + pin_jout

	RaySegBeg = cuRotRay1D%PinRayIdxSt(k)
	RaySegEnd = cuRotRay1D%PinRayIdxSt(k + 1) - 1
    DO ir = RaySegBeg, RaySegEnd
      ireg = cuRotRay1D%FsrIdx(ir)
      tau = - cuRotRay1D%LenSeg(ir) * xst(ig, ireg)
	  fsr_src = src(ig, ireg); fsr_phi = 0.0
	  DO ipol = 1, nPolarAngle
	    ExpApp = 1.0 - __EXPF(tau * rsinv(ipol))
	    track_phi = s_track_phi(itrack, ipol)
	    delta_phi = (track_phi - fsr_src) * ExpApp
		s_track_phi(itrack, ipol) = track_phi - delta_phi
		fsr_phi = fsr_phi + wt(ipol, iazi) * delta_phi
	  ENDDO
	  !$ACC ATOMIC
	  phis(ig, ireg) = phis(ig, ireg) + fsr_phi
	ENDDO

	pin_jout = 0.0
	DO ipol = 1, nPolarAngle
	  track_phi = s_track_phi(itrack, ipol)
	  pin_jout = pin_jout + wt(ipol, iazi) * track_phi
	ENDDO
	!$ACC ATOMIC
    jout(2, ig, isurf(2), ipin) = jout(2, ig, isurf(2), ipin) + pin_jout
  ENDDO
ENDDO

ELSE

DO j = cuRotRay1D%RotRayIdxSt(iRotRay + 1) - 1, cuRotRay1D%RotRayIdxSt(iRotRay), -1
  iazi = cuRotRay1D%iAzi(j)
  PinRayBeg = cuRotRay1D%CoreRayIdxSt(j + 1) - 1
  PinRayEnd = cuRotRay1D%CoreRayIdxSt(j)
  DO k = PinRayBeg, PinRayEnd, -1
	ipin = cuRotRay1D%PinIdx(k)
	isurf = cuRotRay1D%SurfIdx(:, k)
	pin_jout = 0.0
	DO ipol = 1, nPolarAngle
	  track_phi = s_track_phi(itrack, ipol)
	  pin_jout = pin_jout + wt(ipol, iazi) * track_phi
	ENDDO
	!$ACC ATOMIC
    jout(1, ig, isurf(2), ipin) = jout(1, ig, isurf(2), ipin) + pin_jout

	RaySegBeg = cuRotRay1D%PinRayIdxSt(k + 1) - 1
	RaySegEnd = cuRotRay1D%PinRayIdxSt(k)
	DO ir = RaySegBeg, RaySegEnd, -1
      ireg = cuRotRay1D%FsrIdx(ir)
      tau = - cuRotRay1D%LenSeg(ir) * xst(ig, ireg)
	  fsr_src = src(ig, ireg); fsr_phi = 0.0
	  DO ipol = 1, nPolarAngle
	    ExpApp = 1.0 - __EXPF(tau * rsinv(ipol))
	    track_phi = s_track_phi(itrack, ipol)
	    delta_phi = (track_phi - fsr_src) * ExpApp
		s_track_phi(itrack, ipol) = track_phi - delta_phi
		fsr_phi = fsr_phi + wt(ipol, iazi) * delta_phi
	  ENDDO
	  !$ACC ATOMIC
	  phis(ig, ireg) = phis(ig, ireg) + fsr_phi
	ENDDO

	pin_jout = 0.0
	DO ipol = 1, nPolarAngle
	  track_phi = s_track_phi(itrack, ipol)
	  pin_jout = pin_jout + wt(ipol, iazi) * track_phi
	ENDDO
	!$ACC ATOMIC
    jout(2, ig, isurf(1), ipin) = jout(2, ig, isurf(1), ipin) + pin_jout
  ENDDO
ENDDO

ENDIF

PhiAngIn(:, ig, cuRotRay1D%PhiAngOutSvIdx(irot, iRotRay)) = s_track_phi(itrack, :)

END SUBROUTINE

!--- Pn Domain Connected (Block Node Major) -------------------------------------------------------

ATTRIBUTES(GLOBAL) SUBROUTINE cuPrepareAngularSource(cuDevice, src, srcm, SrcAng, iazi, lRot)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuDevice_Type) :: cuDevice
REAL(GPU_SOURCE_PRECISION), DEVICE, INTENT(IN) :: src(PN_BLOCK_SIZE, nfsr)
REAL(GPU_SOURCE_PRECISION), DEVICE, INTENT(IN) :: srcm(nMoment, PN_BLOCK_SIZE, nfsr)
REAL(GPU_SOURCE_PRECISION), DEVICE :: SrcAng(nPolarAngle, PN_BLOCK_SIZE, nAziMap, nfsr)
INTEGER, VALUE :: iazi
LOGICAL, VALUE :: lRot

INTEGER :: ipol, ig, ireg
INTEGER :: AziIdx
REAL(GPU_PRECISION) :: temp_src, fsr_SrcAng(2)

ipol = threadIdx%x
ig = threadIdx%y
ireg = threadIdx%z + (blockIdx%x - 1) * blockDim%z
IF (ireg .GT. nfsr) RETURN

IF (ScatOd .EQ. 1) THEN

  AziIdx = iAzi
  fsr_SrcAng = src(ig, ireg)
  temp_src = sum(Comp(1:2, ipol, AziIdx) * srcm(1:2, ig, ireg))
  fsr_SrcAng(1) = fsr_SrcAng(1) + temp_src
  fsr_SrcAng(2) = fsr_SrcAng(2) - temp_src
  SrcAng(ipol, ig, AziMap(AziIdx, 1), ireg) = fsr_SrcAng(1)
  SrcAng(ipol, ig, AziMap(AziIdx, 2), ireg) = fsr_SrcAng(2)

  AziIdx = nAziAngle - iAzi + 1
  fsr_SrcAng = src(ig, ireg)
  temp_src = sum(Comp(1:2, ipol, AziIdx) * srcm(1:2, ig, ireg))
  fsr_SrcAng(1) = fsr_SrcAng(1) + temp_src
  fsr_SrcAng(2) = fsr_SrcAng(2) - temp_src
  SrcAng(ipol, ig, AziMap(AziIdx, 1), ireg) = fsr_SrcAng(1)
  SrcAng(ipol, ig, AziMap(AziIdx, 2), ireg) = fsr_SrcAng(2)

  IF(lRot) THEN 
    AziIdx = nAziAngle/2 + iAzi 
    fsr_SrcAng = src(ig, ireg)
    temp_src = sum(Comp(1:2, ipol, AziIdx) * srcm(1:2, ig, ireg))
    fsr_SrcAng(1) = fsr_SrcAng(1) + temp_src
    fsr_SrcAng(2) = fsr_SrcAng(2) - temp_src
    SrcAng(ipol, ig, AziMap(AziIdx, 1)+4, ireg) = fsr_SrcAng(1)
    SrcAng(ipol, ig, AziMap(AziIdx, 2)+4, ireg) = fsr_SrcAng(2)

    AziIdx = nAziAngle/2 - iAzi + 1 
    fsr_SrcAng = src(ig, ireg)
    temp_src = sum(Comp(1:2, ipol, AziIdx) * srcm(1:2, ig, ireg))
    fsr_SrcAng(1) = fsr_SrcAng(1) + temp_src
    fsr_SrcAng(2) = fsr_SrcAng(2) - temp_src
    SrcAng(ipol, ig, AziMap(AziIdx, 1)+4, ireg) = fsr_SrcAng(1)
    SrcAng(ipol, ig, AziMap(AziIdx, 2)+4, ireg) = fsr_SrcAng(2)
  END IF

ELSEIF (ScatOd .EQ. 2) THEN

  AziIdx = iAzi
  fsr_SrcAng = src(ig, ireg)
  temp_src = sum(Comp(1:2, ipol, AziIdx) * srcm(1:2, ig, ireg))
  fsr_SrcAng(1) = fsr_SrcAng(1) + temp_src
  fsr_SrcAng(2) = fsr_SrcAng(2) - temp_src
  temp_src = sum(Comp(3:5, ipol, AziIdx) * srcm(3:5, ig, ireg))
  fsr_SrcAng = fsr_SrcAng + temp_src
  SrcAng(ipol, ig, AziMap(AziIdx, 1), ireg) = fsr_SrcAng(1)
  SrcAng(ipol, ig, AziMap(AziIdx, 2), ireg) = fsr_SrcAng(2)

  AziIdx = nAziAngle - iAzi + 1
  fsr_SrcAng = src(ig, ireg)
  temp_src = sum(Comp(1:2, ipol, AziIdx) * srcm(1:2, ig, ireg))
  fsr_SrcAng(1) = fsr_SrcAng(1) + temp_src
  fsr_SrcAng(2) = fsr_SrcAng(2) - temp_src
  temp_src = sum(Comp(3:5, ipol, AziIdx) * srcm(3:5, ig, ireg))
  fsr_SrcAng = fsr_SrcAng + temp_src
  SrcAng(ipol, ig, AziMap(AziIdx, 1), ireg) = fsr_SrcAng(1)
  SrcAng(ipol, ig, AziMap(AziIdx, 2), ireg) = fsr_SrcAng(2)
  
  IF(lRot) THEN 
    AziIdx = nAziAngle/2 + iAzi 
    fsr_SrcAng = src(ig, ireg)
    temp_src = sum(Comp(1:2, ipol, AziIdx) * srcm(1:2, ig, ireg))
    fsr_SrcAng(1) = fsr_SrcAng(1) + temp_src
    fsr_SrcAng(2) = fsr_SrcAng(2) - temp_src
    temp_src = sum(Comp(3:5, ipol, AziIdx) * srcm(3:5, ig, ireg))
    fsr_SrcAng = fsr_SrcAng + temp_src
    SrcAng(ipol, ig, AziMap(AziIdx, 1)+4, ireg) = fsr_SrcAng(1)
    SrcAng(ipol, ig, AziMap(AziIdx, 2)+4, ireg) = fsr_SrcAng(2)

    AziIdx = nAziAngle/2 - iAzi + 1 
    fsr_SrcAng = src(ig, ireg)
    temp_src = sum(Comp(1:2, ipol, AziIdx) * srcm(1:2, ig, ireg))
    fsr_SrcAng(1) = fsr_SrcAng(1) + temp_src
    fsr_SrcAng(2) = fsr_SrcAng(2) - temp_src
    temp_src = sum(Comp(3:5, ipol, AziIdx) * srcm(3:5, ig, ireg))
    fsr_SrcAng = fsr_SrcAng + temp_src
    SrcAng(ipol, ig, AziMap(AziIdx, 1)+4, ireg) = fsr_SrcAng(1)
    SrcAng(ipol, ig, AziMap(AziIdx, 2)+4, ireg) = fsr_SrcAng(2)
  END IF

ELSEIF (ScatOd .EQ. 3) THEN

  AziIdx = iAzi
  fsr_SrcAng = src(ig, ireg)
  temp_src = sum(Comp(1:2, ipol, AziIdx) * srcm(1:2, ig, ireg))
  fsr_SrcAng(1) = fsr_SrcAng(1) + temp_src
  fsr_SrcAng(2) = fsr_SrcAng(2) - temp_src
  temp_src = sum(Comp(3:5, ipol, AziIdx) * srcm(3:5, ig, ireg))
  fsr_SrcAng = fsr_SrcAng + temp_src
  temp_src = sum(Comp(6:9, ipol, AziIdx) * srcm(6:9, ig, ireg))
  fsr_SrcAng(1) = fsr_SrcAng(1) + temp_src
  fsr_SrcAng(2) = fsr_SrcAng(2) - temp_src
  SrcAng(ipol, ig, AziMap(AziIdx, 1), ireg) = fsr_SrcAng(1)
  SrcAng(ipol, ig, AziMap(AziIdx, 2), ireg) = fsr_SrcAng(2)

  AziIdx = nAziAngle - iAzi + 1
  fsr_SrcAng = src(ig, ireg)
  temp_src = sum(Comp(1:2, ipol, AziIdx) * srcm(1:2, ig, ireg))
  fsr_SrcAng(1) = fsr_SrcAng(1) + temp_src
  fsr_SrcAng(2) = fsr_SrcAng(2) - temp_src
  temp_src = sum(Comp(3:5, ipol, AziIdx) * srcm(3:5, ig, ireg))
  fsr_SrcAng = fsr_SrcAng + temp_src
  temp_src = sum(Comp(6:9, ipol, AziIdx) * srcm(6:9, ig, ireg))
  fsr_SrcAng(1) = fsr_SrcAng(1) + temp_src
  fsr_SrcAng(2) = fsr_SrcAng(2) - temp_src
  SrcAng(ipol, ig, AziMap(AziIdx, 1), ireg) = fsr_SrcAng(1)
  SrcAng(ipol, ig, AziMap(AziIdx, 2), ireg) = fsr_SrcAng(2)

  IF(lRot) THEN 
    AziIdx = nAziAngle/2 + iAzi 
    fsr_SrcAng = src(ig, ireg)
    temp_src = sum(Comp(1:2, ipol, AziIdx) * srcm(1:2, ig, ireg))
    fsr_SrcAng(1) = fsr_SrcAng(1) + temp_src
    fsr_SrcAng(2) = fsr_SrcAng(2) - temp_src
    temp_src = sum(Comp(3:5, ipol, AziIdx) * srcm(3:5, ig, ireg))
    fsr_SrcAng = fsr_SrcAng + temp_src
    temp_src = sum(Comp(6:9, ipol, AziIdx) * srcm(6:9, ig, ireg))
    fsr_SrcAng(1) = fsr_SrcAng(1) + temp_src
    fsr_SrcAng(2) = fsr_SrcAng(2) - temp_src
    SrcAng(ipol, ig, AziMap(AziIdx, 1)+4, ireg) = fsr_SrcAng(1)
    SrcAng(ipol, ig, AziMap(AziIdx, 2)+4, ireg) = fsr_SrcAng(2)

    AziIdx = nAziAngle/2 - iAzi + 1 
    fsr_SrcAng = src(ig, ireg)
    temp_src = sum(Comp(1:2, ipol, AziIdx) * srcm(1:2, ig, ireg))
    fsr_SrcAng(1) = fsr_SrcAng(1) + temp_src
    fsr_SrcAng(2) = fsr_SrcAng(2) - temp_src
    temp_src = sum(Comp(3:5, ipol, AziIdx) * srcm(3:5, ig, ireg))
    fsr_SrcAng = fsr_SrcAng + temp_src
    temp_src = sum(Comp(6:9, ipol, AziIdx) * srcm(6:9, ig, ireg))
    fsr_SrcAng(1) = fsr_SrcAng(1) + temp_src
    fsr_SrcAng(2) = fsr_SrcAng(2) - temp_src
    SrcAng(ipol, ig, AziMap(AziIdx, 1)+4, ireg) = fsr_SrcAng(1)
    SrcAng(ipol, ig, AziMap(AziIdx, 2)+4, ireg) = fsr_SrcAng(2)
  END IF
ENDIF

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuIncrementPseudoMoment(cuDevice, phis, phim, phia, iazi, lRot)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuDevice_Type) :: cuDevice
REAL(GPU_FLUX_PRECISION), DEVICE, INTENT(IN) :: phia(nPolarAngle, PN_BLOCK_SIZE, nAziMap, nfsr)
REAL(GPU_FLUX_PRECISION), DEVICE :: phis(PN_BLOCK_SIZE, nfsr)
REAL(GPU_FLUX_PRECISION), DEVICE :: phim(nMoment, PN_BLOCK_SIZE, nfsr)
INTEGER, VALUE :: iazi
LOGICAL, VALUE :: lRot

INTEGER :: ipol, ig, ireg
INTEGER :: AziIdx
REAL(GPU_PRECISION) :: phia_sum, phia_diff, fsr_phis, fsr_phim(9)

ig = threadIdx%x
ireg = threadIdx%y + (blockIdx%x - 1) * blockDim%y
IF (ireg .GT. nfsr) RETURN

IF (ScatOd .EQ. 1) THEN

  fsr_phis = 0.0; fsr_phim = 0.0

  AziIdx = iAzi
  DO ipol = 1, nPolarAngle
    phia_sum = (phia(ipol, ig, AziMap(AziIdx, 1), ireg) + phia(ipol, ig, AziMap(AziIdx, 2), ireg))
	  phia_diff = (phia(ipol, ig, AziMap(AziIdx, 1), ireg) - phia(ipol, ig, AziMap(AziIdx, 2), ireg))
    fsr_phis = fsr_phis + wt(ipol, AziIdx) * phia_sum
    fsr_phim(1:2) = fsr_phim(1:2) + mwt(1:2, ipol, AziIdx) * phia_diff
  ENDDO

  AziIdx = nAziAngle - iAzi + 1
  DO ipol = 1, nPolarAngle
    phia_sum = (phia(ipol, ig, AziMap(AziIdx, 1), ireg) + phia(ipol, ig, AziMap(AziIdx, 2), ireg))
	  phia_diff = (phia(ipol, ig, AziMap(AziIdx, 1), ireg) - phia(ipol, ig, AziMap(AziIdx, 2), ireg))
    fsr_phis = fsr_phis + wt(ipol, AziIdx) * phia_sum
    fsr_phim(1:2) = fsr_phim(1:2) + mwt(1:2, ipol, AziIdx) * phia_diff
  ENDDO

  IF(lRot) THEN 
    AziIdx = nAziAngle/2 + iAzi
    DO ipol = 1, nPolarAngle
      phia_sum = (phia(ipol, ig, AziMap(AziIdx, 1)+4, ireg) + phia(ipol, ig, AziMap(AziIdx, 2)+4, ireg))
	    phia_diff = (phia(ipol, ig, AziMap(AziIdx, 1)+4, ireg) - phia(ipol, ig, AziMap(AziIdx, 2)+4, ireg))
      fsr_phis = fsr_phis + wt(ipol, AziIdx) * phia_sum
      fsr_phim(1:2) = fsr_phim(1:2) + mwt(1:2, ipol, AziIdx) * phia_diff
    ENDDO

    AziIdx = nAziAngle/2 - iAzi + 1
    DO ipol = 1, nPolarAngle
      phia_sum = (phia(ipol, ig, AziMap(AziIdx, 1)+4, ireg) + phia(ipol, ig, AziMap(AziIdx, 2)+4, ireg))
	    phia_diff = (phia(ipol, ig, AziMap(AziIdx, 1)+4, ireg) - phia(ipol, ig, AziMap(AziIdx, 2)+4, ireg))
      fsr_phis = fsr_phis + wt(ipol, AziIdx) * phia_sum
      fsr_phim(1:2) = fsr_phim(1:2) + mwt(1:2, ipol, AziIdx) * phia_diff
    ENDDO
  END IF

  phis(ig, ireg) = phis(ig, ireg) + fsr_phis
  phim(1:2, ig, ireg) = phim(1:2, ig, ireg) + fsr_phim(1:2)

ELSEIF (ScatOd .EQ. 2) THEN

  fsr_phis = 0.0; fsr_phim = 0.0

  AziIdx = iAzi
  DO ipol = 1, nPolarAngle
    phia_sum = (phia(ipol, ig, AziMap(AziIdx, 1), ireg) + phia(ipol, ig, AziMap(AziIdx, 2), ireg))
	  phia_diff = (phia(ipol, ig, AziMap(AziIdx, 1), ireg) - phia(ipol, ig, AziMap(AziIdx, 2), ireg))
    fsr_phis = fsr_phis + wt(ipol, AziIdx) * phia_sum
    fsr_phim(1:2) = fsr_phim(1:2) + mwt(1:2, ipol, AziIdx) * phia_diff
    fsr_phim(3:5) = fsr_phim(3:5) + mwt(3:5, ipol, AziIdx) * phia_sum
  ENDDO

  AziIdx = nAziAngle - iAzi + 1
  DO ipol = 1, nPolarAngle
    phia_sum = (phia(ipol, ig, AziMap(AziIdx, 1), ireg) + phia(ipol, ig, AziMap(AziIdx, 2), ireg))
	  phia_diff = (phia(ipol, ig, AziMap(AziIdx, 1), ireg) - phia(ipol, ig, AziMap(AziIdx, 2), ireg))
    fsr_phis = fsr_phis + wt(ipol, AziIdx) * phia_sum
    fsr_phim(1:2) = fsr_phim(1:2) + mwt(1:2, ipol, AziIdx) * phia_diff
    fsr_phim(3:5) = fsr_phim(3:5) + mwt(3:5, ipol, AziIdx) * phia_sum
  ENDDO

  IF(lRot) THEN 
    AziIdx = nAziAngle/2 + iAzi
    DO ipol = 1, nPolarAngle
      phia_sum = (phia(ipol, ig, AziMap(AziIdx, 1)+4, ireg) + phia(ipol, ig, AziMap(AziIdx, 2)+4, ireg))
	    phia_diff = (phia(ipol, ig, AziMap(AziIdx, 1)+4, ireg) - phia(ipol, ig, AziMap(AziIdx, 2)+4, ireg))
      fsr_phis = fsr_phis + wt(ipol, AziIdx) * phia_sum
      fsr_phim(1:2) = fsr_phim(1:2) + mwt(1:2, ipol, AziIdx) * phia_diff
      fsr_phim(3:5) = fsr_phim(3:5) + mwt(3:5, ipol, AziIdx) * phia_sum
    ENDDO

    AziIdx = nAziAngle/2 - iAzi + 1
    DO ipol = 1, nPolarAngle
      phia_sum = (phia(ipol, ig, AziMap(AziIdx, 1)+4, ireg) + phia(ipol, ig, AziMap(AziIdx, 2)+4, ireg))
	    phia_diff = (phia(ipol, ig, AziMap(AziIdx, 1)+4, ireg) - phia(ipol, ig, AziMap(AziIdx, 2)+4, ireg))
      fsr_phis = fsr_phis + wt(ipol, AziIdx) * phia_sum
      fsr_phim(1:2) = fsr_phim(1:2) + mwt(1:2, ipol, AziIdx) * phia_diff
      fsr_phim(3:5) = fsr_phim(3:5) + mwt(3:5, ipol, AziIdx) * phia_sum
    ENDDO
  END IF

  phis(ig, ireg) = phis(ig, ireg) + fsr_phis
  phim(1:2, ig, ireg) = phim(1:2, ig, ireg) + fsr_phim(1:2)
  phim(3:5, ig, ireg) = phim(3:5, ig, ireg) + fsr_phim(3:5)

ELSEIF (ScatOd .EQ. 3) THEN

  fsr_phis = 0.0; fsr_phim = 0.0

  AziIdx = iAzi
  DO ipol = 1, nPolarAngle
    phia_sum = (phia(ipol, ig, AziMap(AziIdx, 1), ireg) + phia(ipol, ig, AziMap(AziIdx, 2), ireg))
	  phia_diff = (phia(ipol, ig, AziMap(AziIdx, 1), ireg) - phia(ipol, ig, AziMap(AziIdx, 2), ireg))
    fsr_phis = fsr_phis + wt(ipol, AziIdx) * phia_sum
    fsr_phim(1:2) = fsr_phim(1:2) + mwt(1:2, ipol, AziIdx) * phia_diff
    fsr_phim(3:5) = fsr_phim(3:5) + mwt(3:5, ipol, AziIdx) * phia_sum
	  fsr_phim(6:9) = fsr_phim(6:9) + mwt(6:9, ipol, AziIdx) * phia_diff
  ENDDO

  AziIdx = nAziAngle - iAzi + 1
  DO ipol = 1, nPolarAngle
    phia_sum = (phia(ipol, ig, AziMap(AziIdx, 1), ireg) + phia(ipol, ig, AziMap(AziIdx, 2), ireg))
	  phia_diff = (phia(ipol, ig, AziMap(AziIdx, 1), ireg) - phia(ipol, ig, AziMap(AziIdx, 2), ireg))
    fsr_phis = fsr_phis + wt(ipol, AziIdx) * phia_sum
    fsr_phim(1:2) = fsr_phim(1:2) + mwt(1:2, ipol, AziIdx) * phia_diff
    fsr_phim(3:5) = fsr_phim(3:5) + mwt(3:5, ipol, AziIdx) * phia_sum
	  fsr_phim(6:9) = fsr_phim(6:9) + mwt(6:9, ipol, AziIdx) * phia_diff
  ENDDO

  IF(lRot) THEN 
    AziIdx = nAziAngle/2 + iAzi
    DO ipol = 1, nPolarAngle
      phia_sum = (phia(ipol, ig, AziMap(AziIdx, 1)+4, ireg) + phia(ipol, ig, AziMap(AziIdx, 2)+4, ireg))
	    phia_diff = (phia(ipol, ig, AziMap(AziIdx, 1)+4, ireg) - phia(ipol, ig, AziMap(AziIdx, 2)+4, ireg))
      fsr_phis = fsr_phis + wt(ipol, AziIdx) * phia_sum
      fsr_phim(1:2) = fsr_phim(1:2) + mwt(1:2, ipol, AziIdx) * phia_diff
      fsr_phim(3:5) = fsr_phim(3:5) + mwt(3:5, ipol, AziIdx) * phia_sum
	    fsr_phim(6:9) = fsr_phim(6:9) + mwt(6:9, ipol, AziIdx) * phia_diff
    ENDDO

    AziIdx = nAziAngle/2 - iAzi + 1
    DO ipol = 1, nPolarAngle
      phia_sum = (phia(ipol, ig, AziMap(AziIdx, 1)+4, ireg) + phia(ipol, ig, AziMap(AziIdx, 2)+4, ireg))
	    phia_diff = (phia(ipol, ig, AziMap(AziIdx, 1)+4, ireg) - phia(ipol, ig, AziMap(AziIdx, 2)+4, ireg))
      fsr_phis = fsr_phis + wt(ipol, AziIdx) * phia_sum
      fsr_phim(1:2) = fsr_phim(1:2) + mwt(1:2, ipol, AziIdx) * phia_diff
      fsr_phim(3:5) = fsr_phim(3:5) + mwt(3:5, ipol, AziIdx) * phia_sum
	    fsr_phim(6:9) = fsr_phim(6:9) + mwt(6:9, ipol, AziIdx) * phia_diff
    ENDDO
  END IF

  phis(ig, ireg) = phis(ig, ireg) + fsr_phis
  phim(1:2, ig, ireg) = phim(1:2, ig, ireg) + fsr_phim(1:2)
  phim(3:5, ig, ireg) = phim(3:5, ig, ireg) + fsr_phim(3:5)
  phim(6:9, ig, ireg) = phim(6:9, ig, ireg) + fsr_phim(6:9)

ENDIF

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuComputeFluxMoment(cuGeometry, phis, phim, src, srcm, xst, iz)
USE PARAM
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuGeometry_Type), INTENT(IN) :: cuGeometry
REAL(GPU_SOURCE_PRECISION), DEVICE, INTENT(IN) :: src(PN_BLOCK_SIZE, nfsr)
REAL(GPU_SOURCE_PRECISION), DEVICE, INTENT(IN) :: srcm(nMoment, PN_BLOCK_SIZE, nfsr)
REAL(GPU_SOURCE_PRECISION), DEVICE, INTENT(IN) :: xst(PN_BLOCK_SIZE, nfsr)
REAL(GPU_FLUX_PRECISION), DEVICE :: phis(PN_BLOCK_SIZE, nfsr)
REAL(GPU_FLUX_PRECISION), DEVICE :: phim(nMoment, PN_BLOCK_SIZE, nfsr)
INTEGER, VALUE :: iz

INTEGER :: ig, ifsr, itype
REAL(GPU_PRECISION) :: rsigv

ig = threadIdx%x
ifsr = threadIdx%y + (blockIdx%x - 1) * blockDim%y
IF (ifsr .GT. nfsr) RETURN

itype = cuGeometry%AxType(iz)
rsigv = 1.0 / (xst(ig, ifsr) * cuGeometry%fsrVol(ifsr, itype))
phis(ig, ifsr) = phis(ig, ifsr) * rsigv + src(ig, ifsr)
IF (ScatOd .GE. 1) phim(1:2, ig, ifsr) = phim(1:2, ig, ifsr) * rsigv + srcm(1:2, ig, ifsr) * rthree
IF (ScatOd .GE. 2) phim(3:5, ig, ifsr) = phim(3:5, ig, ifsr) * rsigv + srcm(3:5, ig, ifsr) * rfive
IF (ScatOd .EQ. 3) phim(6:9, ig, ifsr) = phim(6:9, ig, ifsr) * rsigv + srcm(6:9, ig, ifsr) * rseven

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cu1DFastRayTraceP1(cuRotRay1D, cuGeometry, cuDevice, phia, xst, SrcAng, jout,         &
                                                 PhiAngIn, lJout, iazi)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuDevice_Type), INTENT(IN) :: cuDevice
TYPE(cuGeometry_Type), INTENT(IN) :: cuGeometry
TYPE(cuRotRay1D_Type), INTENT(IN) :: cuRotRay1D
REAL(GPU_SOURCE_PRECISION), DEVICE, INTENT(IN) :: xst(PN_BLOCK_SIZE, nfsr)
REAL(GPU_SOURCE_PRECISION), DEVICE, INTENT(IN) :: SrcAng(nPolarAngle, PN_BLOCK_SIZE, nAziMap, nfsr)
REAL(GPU_PRECISION), DEVICE :: PhiAngIn(nPolarAngle, PN_BLOCK_SIZE, nPhiAngSv)
REAL(GPU_FLUX_PRECISION), DEVICE :: phia(nPolarAngle, PN_BLOCK_SIZE, nAziMap, nfsr)
REAL(GPU_FLUX_PRECISION), DEVICE :: jout(3, PN_BLOCK_SIZE, 4, nxy)
LOGICAL, VALUE :: lJout
INTEGER, VALUE :: iazi

INTEGER :: iRotRay, iray, ipol, ig, irot, itrack
INTEGER :: iazir

ipol = threadIdx%x
ig = threadIdx%y
iray = threadIdx%z + (blockIdx%x - 1) * blockDim%z
irot = blockIdx%y
itrack = threadIdx%x + (threadIdx%y - 1) * blockDim%x + (threadIdx%z - 1) * blockDim%x * blockDim%y
IF (iray .GT. cuGeometry%RotRayCount(iazi)) RETURN

iRotRay = cuGeometry%RotRayList(iray, iazi)

IF(cuGeometry%lRot) THEN 
  iazir = cuGeometry%nAziAngle - iazi + 1
  IF (.NOT. lJout) THEN
    CALL cuP1Track1DFastRay1stRot(cuRotRay1D, cuDevice, phia, xst, SrcAng, PhiAngIn,								        &
                             ipol, ig, irot, itrack, iRotRay, iazi, iazir)
  ELSE
    CALL cuP1Track1DFastRay2ndRot(cuRotRay1D, cuDevice, phia, xst, SrcAng, jout, PhiAngIn,                               &
                             ipol, ig, irot, itrack, iRotRay, iazi, iazir)
  ENDIF
ELSE
  IF (.NOT. lJout) THEN
    CALL cuP1Track1DFastRay1st(cuRotRay1D, cuDevice, phia, xst, SrcAng, PhiAngIn,								        &
                             ipol, ig, irot, itrack, iRotRay)
  ELSE
    CALL cuP1Track1DFastRay2nd(cuRotRay1D, cuDevice, phia, xst, SrcAng, jout, PhiAngIn,                               &
                             ipol, ig, irot, itrack, iRotRay)
  ENDIF
END IF

END SUBROUTINE

ATTRIBUTES(DEVICE) SUBROUTINE cuP1Track1DFastRay1st(cuRotRay1D, cuDevice, phia, xst, SrcAng, PhiAngIn,              &
                                                    ipol, ig, irot, itrack, iRotRay)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuRotRay1D_Type), INTENT(IN) :: cuRotRay1D
TYPE(cuDevice_Type), INTENT(IN) :: cuDevice
REAL(GPU_SOURCE_PRECISION), DEVICE, INTENT(IN) :: xst(PN_BLOCK_SIZE, nfsr)
REAL(GPU_SOURCE_PRECISION), DEVICE, INTENT(IN) :: SrcAng(nPolarAngle, PN_BLOCK_SIZE, nAziMap, nfsr)
REAL(GPU_PRECISION), DEVICE :: PhiAngIn(nPolarAngle, PN_BLOCK_SIZE, nPhiAngSv)
REAL(GPU_FLUX_PRECISION), DEVICE :: phia(nPolarAngle, PN_BLOCK_SIZE, nAziMap, nfsr)
INTEGER :: ipol, ig, irot, itrack, iRotRay

INTEGER :: RaySegBeg, RaySegEnd
INTEGER :: iazi, idir, ireg, AziSvIdx
INTEGER :: j, ir
REAL(GPU_PRECISION) :: delta_phi, track_phi
REAL(GPU_PRECISION) :: tau, ExpApp

track_phi = PhiAngIn(ipol, ig, cuRotRay1D%PhiAngInSvIdx(irot, iRotRay))

IF (irot .EQ. 1) THEN

DO j = cuRotRay1D%RotRayIdxSt(iRotRay), cuRotRay1D%RotRayIdxSt(iRotRay + 1) -1
  iazi = cuRotRay1D%iAzi(j)
  idir = cuRotRay1D%iDir(j)
  idir = dir(idir, irot)
  AziSvIdx = AziMap(iazi, idir)
  RaySegBeg = cuRotRay1D%PinRayIdxSt(cuRotRay1D%CoreRayIdxSt(j))
  RaySegEnd = cuRotRay1D%PinRayIdxSt(cuRotRay1D%CoreRayIdxSt(j + 1)) - 1
  DO ir = RaySegBeg, RaySegEnd
    ireg = cuRotRay1D%FsrIdx(ir)
    tau = - cuRotRay1D%LenSeg(ir) * xst(ig, ireg)
	ExpApp = 1.0 - __EXPF(tau * rsinv(ipol))
	delta_phi = (track_phi - SrcAng(ipol, ig, AziSvIdx, ireg)) * ExpApp
	track_phi = track_phi - delta_phi
	!$ACC ATOMIC
	phia(ipol, ig, AziSvIdx, ireg) = phia(ipol, ig, AziSvIdx, ireg) + delta_phi
  ENDDO
ENDDO

ELSE

DO j = cuRotRay1D%RotRayIdxSt(iRotRay + 1) - 1, cuRotRay1D%RotRayIdxSt(iRotRay), -1
  iazi = cuRotRay1D%iAzi(j)
  idir = cuRotRay1D%iDir(j)
  idir = dir(idir, irot)
  AziSvIdx = AziMap(iazi, idir)
  RaySegBeg = cuRotRay1D%PinRayIdxSt(cuRotRay1D%CoreRayIdxSt(j + 1)) - 1
  RaySegEnd = cuRotRay1D%PinRayIdxSt(cuRotRay1D%CoreRayIdxSt(j))
  DO ir = RaySegBeg, RaySegEnd, -1
    ireg = cuRotRay1D%FsrIdx(ir)
    tau = - cuRotRay1D%LenSeg(ir) * xst(ig, ireg)
	ExpApp = 1.0 - __EXPF(tau * rsinv(ipol))
	delta_phi = (track_phi - SrcAng(ipol, ig, AziSvIdx, ireg)) * ExpApp
	track_phi = track_phi - delta_phi
	!$ACC ATOMIC
	phia(ipol, ig, AziSvIdx, ireg) = phia(ipol, ig, AziSvIdx, ireg) + delta_phi
  ENDDO
ENDDO

ENDIF

PhiAngIn(ipol, ig, cuRotRay1D%PhiAngOutSvIdx(irot, iRotRay)) = track_phi

END SUBROUTINE

ATTRIBUTES(DEVICE) SUBROUTINE cuP1Track1DFastRay2nd(cuRotRay1D, cuDevice, phia, xst, SrcAng, jout, PhiAngIn,        &
                                                    ipol, ig, irot, itrack, iRotRay)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuRotRay1D_Type), INTENT(IN) :: cuRotRay1D
TYPE(cuDevice_Type), INTENT(IN) :: cuDevice
REAL(GPU_SOURCE_PRECISION), DEVICE, INTENT(IN) :: xst(PN_BLOCK_SIZE, nfsr)
REAL(GPU_SOURCE_PRECISION), DEVICE, INTENT(IN) :: SrcAng(nPolarAngle, PN_BLOCK_SIZE, nAziMap, nfsr)
REAL(GPU_PRECISION), DEVICE :: PhiAngIn(nPolarAngle, PN_BLOCK_SIZE, nPhiAngSv)
REAL(GPU_FLUX_PRECISION), DEVICE :: phia(nPolarAngle, PN_BLOCK_SIZE, nAziMap, nfsr)
REAL(GPU_FLUX_PRECISION), DEVICE :: jout(3, PN_BLOCK_SIZE, 4, nxy)
INTEGER :: ipol, ig, irot, itrack, iRotRay

INTEGER :: PinRayBeg, PinRayEnd
INTEGER :: RaySegBeg, RaySegEnd
INTEGER :: iazi, idir, ireg, ipin, isurf(2), AziSvIdx
INTEGER :: j, k, ir
REAL(GPU_PRECISION) :: pin_jout
REAL(GPU_PRECISION) :: delta_phi, track_phi
REAL(GPU_PRECISION) :: tau, ExpApp

track_phi = PhiAngIn(ipol, ig, cuRotRay1D%PhiAngInSvIdx(irot, iRotRay))

IF (irot .EQ. 1) THEN

DO j = cuRotRay1D%RotRayIdxSt(iRotRay), cuRotRay1D%RotRayIdxSt(iRotRay + 1) - 1
  iazi = cuRotRay1D%iAzi(j)
  idir = cuRotRay1D%iDir(j)
  idir = dir(idir, irot)
  AziSvIdx = AziMap(iazi, idir)
  PinRayBeg = cuRotRay1D%CoreRayIdxSt(j)
  PinRayEnd = cuRotRay1D%CoreRayIdxSt(j + 1) - 1
  DO k = PinRayBeg, PinRayEnd
	ipin = cuRotRay1D%PinIdx(k)
	isurf = cuRotRay1D%SurfIdx(:, k)
	pin_jout = wt(ipol, iazi) * track_phi
	!$ACC ATOMIC
    jout(1, ig, isurf(1), ipin) = jout(1, ig, isurf(1), ipin) + pin_jout

	RaySegBeg = cuRotRay1D%PinRayIdxSt(k)
	RaySegEnd = cuRotRay1D%PinRayIdxSt(k + 1) - 1
    DO ir = RaySegBeg, RaySegEnd
      ireg = cuRotRay1D%FsrIdx(ir)
      tau = - cuRotRay1D%LenSeg(ir) * xst(ig, ireg)
	  ExpApp = 1.0 - __EXPF(tau * rsinv(ipol))
	  delta_phi = (track_phi - SrcAng(ipol, ig, AziSvIdx, ireg)) * ExpApp
	  track_phi = track_phi - delta_phi
	  !$ACC ATOMIC
	  phia(ipol, ig, AziSvIdx, ireg) = phia(ipol, ig, AziSvIdx, ireg) + delta_phi
	ENDDO

	pin_jout = wt(ipol, iazi) * track_phi
	!$ACC ATOMIC
    jout(2, ig, isurf(2), ipin) = jout(2, ig, isurf(2), ipin) + pin_jout
  ENDDO
ENDDO

ELSE

DO j = cuRotRay1D%RotRayIdxSt(iRotRay + 1) - 1, cuRotRay1D%RotRayIdxSt(iRotRay), -1
  iazi = cuRotRay1D%iAzi(j)
  idir = cuRotRay1D%iDir(j)
  idir = dir(idir, irot)
  AziSvIdx = AziMap(iazi, idir)
  PinRayBeg = cuRotRay1D%CoreRayIdxSt(j + 1) - 1
  PinRayEnd = cuRotRay1D%CoreRayIdxSt(j)
  DO k = PinRayBeg, PinRayEnd, -1
	ipin = cuRotRay1D%PinIdx(k)
	isurf = cuRotRay1D%SurfIdx(:, k)
	pin_jout = wt(ipol, iazi) * track_phi
	!$ACC ATOMIC
    jout(1, ig, isurf(2), ipin) = jout(1, ig, isurf(2), ipin) + pin_jout

	RaySegBeg = cuRotRay1D%PinRayIdxSt(k + 1) - 1
	RaySegEnd = cuRotRay1D%PinRayIdxSt(k)
	DO ir = RaySegBeg, RaySegEnd, -1
      ireg = cuRotRay1D%FsrIdx(ir)
      tau = - cuRotRay1D%LenSeg(ir) * xst(ig, ireg)
	  ExpApp = 1.0 - __EXPF(tau * rsinv(ipol))
	  delta_phi = (track_phi - SrcAng(ipol, ig, AziSvIdx, ireg)) * ExpApp
	  track_phi = track_phi - delta_phi
	  !$ACC ATOMIC
	  phia(ipol, ig, AziSvIdx, ireg) = phia(ipol, ig, AziSvIdx, ireg) + delta_phi
	ENDDO

	pin_jout = wt(ipol, iazi) * track_phi
	!$ACC ATOMIC
    jout(2, ig, isurf(1), ipin) = jout(2, ig, isurf(1), ipin) + pin_jout
  ENDDO
ENDDO

ENDIF

PhiAngIn(ipol, ig, cuRotRay1D%PhiAngOutSvIdx(irot, iRotRay)) = track_phi

END SUBROUTINE

ATTRIBUTES(DEVICE) SUBROUTINE cuP1Track1DFastRay1stRot(cuRotRay1D, cuDevice, phia, xst, SrcAng, PhiAngIn,              &
                                                       ipol, ig, irot, itrack, iRotRay, iazi0, iazi0r)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuRotRay1D_Type), INTENT(IN) :: cuRotRay1D
TYPE(cuDevice_Type), INTENT(IN) :: cuDevice
REAL(GPU_SOURCE_PRECISION), DEVICE, INTENT(IN) :: xst(PN_BLOCK_SIZE, nfsr)
REAL(GPU_SOURCE_PRECISION), DEVICE, INTENT(IN) :: SrcAng(nPolarAngle, PN_BLOCK_SIZE, nAziMap, nfsr)
REAL(GPU_PRECISION), DEVICE :: PhiAngIn(nPolarAngle, PN_BLOCK_SIZE, nPhiAngSv)
REAL(GPU_FLUX_PRECISION), DEVICE :: phia(nPolarAngle, PN_BLOCK_SIZE, nAziMap, nfsr)
INTEGER :: ipol, ig, irot, itrack, iRotRay
INTEGER :: iazi0, iazi0r

INTEGER :: RaySegBeg, RaySegEnd
INTEGER :: iazi, idir, ireg, AziSvIdx
INTEGER :: j, ir
REAL(GPU_PRECISION) :: delta_phi, track_phi
REAL(GPU_PRECISION) :: tau, ExpApp

track_phi = PhiAngIn(ipol, ig, cuRotRay1D%PhiAngInSvIdx(irot, iRotRay))

IF (irot .EQ. 1) THEN

DO j = cuRotRay1D%RotRayIdxSt(iRotRay), cuRotRay1D%RotRayIdxSt(iRotRay + 1) -1
  iazi = cuRotRay1D%iAzi(j)
  idir = cuRotRay1D%iDir(j)
  idir = dir(idir, irot)
  AziSvIdx = AziMap(iazi, idir)
  IF(iazi .NE. iazi0 .AND. iazi .NE. iazi0r) THEN 
    AziSvIdx = AziSvIdx+4
  END IF
  RaySegBeg = cuRotRay1D%PinRayIdxSt(cuRotRay1D%CoreRayIdxSt(j))
  RaySegEnd = cuRotRay1D%PinRayIdxSt(cuRotRay1D%CoreRayIdxSt(j + 1)) - 1
  DO ir = RaySegBeg, RaySegEnd
    ireg = cuRotRay1D%FsrIdx(ir)
    tau = - cuRotRay1D%LenSeg(ir) * xst(ig, ireg)
	ExpApp = 1.0 - __EXPF(tau * rsinv(ipol))
	delta_phi = (track_phi - SrcAng(ipol, ig, AziSvIdx, ireg)) * ExpApp
	track_phi = track_phi - delta_phi
	!$ACC ATOMIC
	phia(ipol, ig, AziSvIdx, ireg) = phia(ipol, ig, AziSvIdx, ireg) + delta_phi
  ENDDO
ENDDO

ELSE

DO j = cuRotRay1D%RotRayIdxSt(iRotRay + 1) - 1, cuRotRay1D%RotRayIdxSt(iRotRay), -1
  iazi = cuRotRay1D%iAzi(j)
  idir = cuRotRay1D%iDir(j)
  idir = dir(idir, irot)
  AziSvIdx = AziMap(iazi, idir)
  IF(iazi .NE. iazi0 .AND. iazi .NE. iazi0r) THEN 
    AziSvIdx = AziSvIdx+4
  END IF
  RaySegBeg = cuRotRay1D%PinRayIdxSt(cuRotRay1D%CoreRayIdxSt(j + 1)) - 1
  RaySegEnd = cuRotRay1D%PinRayIdxSt(cuRotRay1D%CoreRayIdxSt(j))
  DO ir = RaySegBeg, RaySegEnd, -1
    ireg = cuRotRay1D%FsrIdx(ir)
    tau = - cuRotRay1D%LenSeg(ir) * xst(ig, ireg)
	ExpApp = 1.0 - __EXPF(tau * rsinv(ipol))
	delta_phi = (track_phi - SrcAng(ipol, ig, AziSvIdx, ireg)) * ExpApp
	track_phi = track_phi - delta_phi
	!$ACC ATOMIC
	phia(ipol, ig, AziSvIdx, ireg) = phia(ipol, ig, AziSvIdx, ireg) + delta_phi
  ENDDO
ENDDO

ENDIF

PhiAngIn(ipol, ig, cuRotRay1D%PhiAngOutSvIdx(irot, iRotRay)) = track_phi

END SUBROUTINE

ATTRIBUTES(DEVICE) SUBROUTINE cuP1Track1DFastRay2ndRot(cuRotRay1D, cuDevice, phia, xst, SrcAng, jout, PhiAngIn,        &
                                                      ipol, ig, irot, itrack, iRotRay, iazi0, iazi0r)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuRotRay1D_Type), INTENT(IN) :: cuRotRay1D
TYPE(cuDevice_Type), INTENT(IN) :: cuDevice
REAL(GPU_SOURCE_PRECISION), DEVICE, INTENT(IN) :: xst(PN_BLOCK_SIZE, nfsr)
REAL(GPU_SOURCE_PRECISION), DEVICE, INTENT(IN) :: SrcAng(nPolarAngle, PN_BLOCK_SIZE, nAziMap, nfsr)
REAL(GPU_PRECISION), DEVICE :: PhiAngIn(nPolarAngle, PN_BLOCK_SIZE, nPhiAngSv)
REAL(GPU_FLUX_PRECISION), DEVICE :: phia(nPolarAngle, PN_BLOCK_SIZE, nAziMap, nfsr)
REAL(GPU_FLUX_PRECISION), DEVICE :: jout(3, PN_BLOCK_SIZE, 4, nxy)
INTEGER :: ipol, ig, irot, itrack, iRotRay
INTEGER :: iazi0, iazi0r

INTEGER :: PinRayBeg, PinRayEnd
INTEGER :: RaySegBeg, RaySegEnd
INTEGER :: iazi, idir, ireg, ipin, isurf(2), AziSvIdx
INTEGER :: j, k, ir
REAL(GPU_PRECISION) :: pin_jout
REAL(GPU_PRECISION) :: delta_phi, track_phi
REAL(GPU_PRECISION) :: tau, ExpApp

track_phi = PhiAngIn(ipol, ig, cuRotRay1D%PhiAngInSvIdx(irot, iRotRay))

IF (irot .EQ. 1) THEN

DO j = cuRotRay1D%RotRayIdxSt(iRotRay), cuRotRay1D%RotRayIdxSt(iRotRay + 1) - 1
  iazi = cuRotRay1D%iAzi(j)
  idir = cuRotRay1D%iDir(j)
  idir = dir(idir, irot)
  AziSvIdx = AziMap(iazi, idir)
  IF(iazi .NE. iazi0 .AND. iazi .NE. iazi0r) THEN 
    AziSvIdx = AziSvIdx+4
  END IF
  PinRayBeg = cuRotRay1D%CoreRayIdxSt(j)
  PinRayEnd = cuRotRay1D%CoreRayIdxSt(j + 1) - 1
  DO k = PinRayBeg, PinRayEnd
	ipin = cuRotRay1D%PinIdx(k)
	isurf = cuRotRay1D%SurfIdx(:, k)
	pin_jout = wt(ipol, iazi) * track_phi
	!$ACC ATOMIC
    jout(1, ig, isurf(1), ipin) = jout(1, ig, isurf(1), ipin) + pin_jout

	RaySegBeg = cuRotRay1D%PinRayIdxSt(k)
	RaySegEnd = cuRotRay1D%PinRayIdxSt(k + 1) - 1
    DO ir = RaySegBeg, RaySegEnd
      ireg = cuRotRay1D%FsrIdx(ir)
      tau = - cuRotRay1D%LenSeg(ir) * xst(ig, ireg)
	  ExpApp = 1.0 - __EXPF(tau * rsinv(ipol))
	  delta_phi = (track_phi - SrcAng(ipol, ig, AziSvIdx, ireg)) * ExpApp
	  track_phi = track_phi - delta_phi
	  !$ACC ATOMIC
	  phia(ipol, ig, AziSvIdx, ireg) = phia(ipol, ig, AziSvIdx, ireg) + delta_phi
	ENDDO

	pin_jout = wt(ipol, iazi) * track_phi
	!$ACC ATOMIC
    jout(2, ig, isurf(2), ipin) = jout(2, ig, isurf(2), ipin) + pin_jout
  ENDDO
ENDDO

ELSE

DO j = cuRotRay1D%RotRayIdxSt(iRotRay + 1) - 1, cuRotRay1D%RotRayIdxSt(iRotRay), -1
  iazi = cuRotRay1D%iAzi(j)
  idir = cuRotRay1D%iDir(j)
  idir = dir(idir, irot)
  AziSvIdx = AziMap(iazi, idir)
  IF(iazi .NE. iazi0 .AND. iazi .NE. iazi0r) THEN 
    AziSvIdx = AziSvIdx+4
  END IF
  PinRayBeg = cuRotRay1D%CoreRayIdxSt(j + 1) - 1
  PinRayEnd = cuRotRay1D%CoreRayIdxSt(j)
  DO k = PinRayBeg, PinRayEnd, -1
	ipin = cuRotRay1D%PinIdx(k)
	isurf = cuRotRay1D%SurfIdx(:, k)
	pin_jout = wt(ipol, iazi) * track_phi
	!$ACC ATOMIC
    jout(1, ig, isurf(2), ipin) = jout(1, ig, isurf(2), ipin) + pin_jout

	RaySegBeg = cuRotRay1D%PinRayIdxSt(k + 1) - 1
	RaySegEnd = cuRotRay1D%PinRayIdxSt(k)
	DO ir = RaySegBeg, RaySegEnd, -1
      ireg = cuRotRay1D%FsrIdx(ir)
      tau = - cuRotRay1D%LenSeg(ir) * xst(ig, ireg)
	  ExpApp = 1.0 - __EXPF(tau * rsinv(ipol))
	  delta_phi = (track_phi - SrcAng(ipol, ig, AziSvIdx, ireg)) * ExpApp
	  track_phi = track_phi - delta_phi
	  !$ACC ATOMIC
	  phia(ipol, ig, AziSvIdx, ireg) = phia(ipol, ig, AziSvIdx, ireg) + delta_phi
	ENDDO

	pin_jout = wt(ipol, iazi) * track_phi
	!$ACC ATOMIC
    jout(2, ig, isurf(1), ipin) = jout(2, ig, isurf(1), ipin) + pin_jout
  ENDDO
ENDDO

ENDIF

PhiAngIn(ipol, ig, cuRotRay1D%PhiAngOutSvIdx(irot, iRotRay)) = track_phi

END SUBROUTINE

END MODULE

#endif
