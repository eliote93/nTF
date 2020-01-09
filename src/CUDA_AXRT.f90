#include <CUDADEFINES.h>
    
#ifdef __PGI

MODULE CUDA_FLATMOC

USE CUDA_MASTER
USE CUDA_UTIL
IMPLICIT NONE

CONTAINS

! ATTRIBUTES(GLOBAL) SUBROUTINE cuFlatRayTrace1st(cuDevice, phis, phim, xst, src, srcm, PhiAngIn, PhiAngOut)
! USE CUDA_CONST
! USE CUDAFOR
! IMPLICIT NONE
! 
! TYPE(cuDevice_Type) :: cuDevice
! REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: xst(:, :, :), src(:, :, :), srcm(:, :, :)
! REAL(GPU_NODAL_PRECISION), DEVICE :: phis(:, :, :), phim(:, :, :), PhiAngIn(:, :, :, :), PhiAngOut(:, :, :, :)
! 
! REAL(GPU_NODAL_PRECISION) :: del_phi, local_phi, local_phim, local_src, local_srcm
! REAL(4) :: tau
! INTEGER :: ig, ipin, ipol, izf, itrack
! 
! REAL(GPU_NODAL_PRECISION), SHARED :: track_phi(nPolar1D, cuDevice%sharedMemoryDim)
! 
! itrack = threadIdx%x
! ig = mod((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1, ng) + 1
! ipin = ((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1) / ng + 1
! 
! IF (ipin .GT. nxyc) RETURN
! 
! !--- Upward Sweep
! 
! track_phi(:, itrack) = PhiAngIn(:, ig, ipin, bottom)
! 
! DO izf = 1, cuDevice%nzSub
!   tau = xst(ig, ipin, izf) * cuDevice%hzSub(izf)
!   local_src = src(ig, ipin, izf); local_srcm = srcm(ig, ipin, izf)
!   local_phi = 0.0; local_phim = 0.0
!   DO ipol = 1, nPolar1D
!     del_phi = (local_src + Comp1D(ipol) * local_srcm - track_phi(ipol, itrack))                                     &
!               * (1.0 - __EXPF(- tau * rcosv1D(ipol)))
!     local_phi = local_phi - wtsurf1D(ipol) * del_phi
!     local_phim = local_phim - mwt1D(ipol) * del_phi
!     track_phi(ipol, itrack) = track_phi(ipol, itrack) + del_phi
!   ENDDO
!   phis(ig, ipin, izf) = local_phi
!   phim(ig, ipin, izf) = local_phim
! ENDDO
! 
! PhiAngOut(:, ig, ipin, top) = track_phi(:, itrack)
! 
! !--- Downward Sweep
! 
! track_phi(:, itrack) = PhiAngIn(:, ig, ipin, top)
! 
! DO izf = cuDevice%nzSub, 1, -1
!   tau = xst(ig, ipin, izf) * cuDevice%hzSub(izf)
!   local_src = src(ig, ipin, izf); local_srcm = srcm(ig, ipin, izf)
!   local_phi = 0.0; local_phim = 0.0
!   DO ipol = 1, nPolar1D
!     del_phi = (local_src - Comp1D(ipol) * local_srcm - track_phi(ipol, itrack))                                     &
!               * (1.0 - __EXPF(- tau * rcosv1D(ipol)))
!     local_phi = local_phi - wtsurf1D(ipol) * del_phi
!     local_phim = local_phim + mwt1D(ipol) * del_phi
!     track_phi(ipol, itrack) = track_phi(ipol, itrack) + del_phi
!   ENDDO
!   phis(ig, ipin, izf) = phis(ig, ipin, izf) + local_phi
!   phim(ig, ipin, izf) = phim(ig, ipin, izf) + local_phim
! ENDDO
! 
! PhiAngOut(:, ig, ipin, bottom) = track_phi(:, itrack)
! 
! END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuFlatRayTrace1st(cuDevice, phis, phim, xst, src, srcm, PhiAngIn, PhiAngOut)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuDevice_Type) :: cuDevice
REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: xst(:, :, :), src(:, :, :), srcm(:, :, :)
REAL(GPU_NODAL_PRECISION), DEVICE :: phis(:, :, :), phim(:, :, :), PhiAngIn(:, :, :, :), PhiAngOut(:, :, :, :)

REAL(GPU_NODAL_PRECISION) :: del_phi, local_phi, local_phim, local_src, local_srcm
REAL(4) :: tau
INTEGER :: ig, ipin, ipol, izf, itrack

REAL(GPU_NODAL_PRECISION), SHARED :: track_phi(nPolar1D, cuDevice%sharedMemoryDim)

itrack = threadIdx%x
ipin = mod((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1, nxyc) + 1
ig = ((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1) / nxyc + 1

IF (ig .GT. ng) RETURN

!--- Upward Sweep

track_phi(:, itrack) = PhiAngIn(:, ipin, ig, bottom)

DO izf = 1, cuDevice%nzSub
  tau = xst(ipin, ig, izf) * cuDevice%hzSub(izf)
  local_src = src(ipin, ig, izf); local_srcm = srcm(ipin, ig, izf)
  local_phi = 0.0; local_phim = 0.0
  DO ipol = 1, nPolar1D
    del_phi = (local_src + Comp1D(ipol) * local_srcm - track_phi(ipol, itrack))                                     &
              * (1.0 - __EXPF(- tau * rcosv1D(ipol)))
    local_phi = local_phi - wtsurf1D(ipol) * del_phi
    local_phim = local_phim - mwt1D(ipol) * del_phi
    track_phi(ipol, itrack) = track_phi(ipol, itrack) + del_phi
  ENDDO
  phis(ipin, ig, izf) = local_phi
  phim(ipin, ig, izf) = local_phim
ENDDO

PhiAngOut(:, ipin, ig, top) = track_phi(:, itrack)

!--- Downward Sweep

track_phi(:, itrack) = PhiAngIn(:, ipin, ig, top)

DO izf = cuDevice%nzSub, 1, -1
  tau = xst(ipin, ig, izf) * cuDevice%hzSub(izf)
  local_src = src(ipin, ig, izf); local_srcm = srcm(ipin, ig, izf)
  local_phi = 0.0; local_phim = 0.0
  DO ipol = 1, nPolar1D
    del_phi = (local_src - Comp1D(ipol) * local_srcm - track_phi(ipol, itrack))                                     &
              * (1.0 - __EXPF(- tau * rcosv1D(ipol)))
    local_phi = local_phi - wtsurf1D(ipol) * del_phi
    local_phim = local_phim + mwt1D(ipol) * del_phi
    track_phi(ipol, itrack) = track_phi(ipol, itrack) + del_phi
  ENDDO
  phis(ipin, ig, izf) = phis(ipin, ig, izf) + local_phi
  phim(ipin, ig, izf) = phim(ipin, ig, izf) + local_phim
ENDDO

PhiAngOut(:, ipin, ig, bottom) = track_phi(:, itrack)

END SUBROUTINE

! ATTRIBUTES(GLOBAL) SUBROUTINE cuFlatRayTrace2nd(cuDevice, xst, src, srcm, PhiAngIn, PhiAngOut, Jout)
! USE CUDA_CONST
! USE CUDAFOR
! IMPLICIT NONE
! 
! TYPE(cuDevice_Type) :: cuDevice
! REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: xst(:, :, :), src(:, :, :), srcm(:, :, :)
! REAL(GPU_NODAL_PRECISION), DEVICE :: PhiAngIn(:, :, :, :), PhiAngOut(:, :, :, :)
! REAL, DEVICE :: Jout(:, :, cuDevice%myzbf :, :, :)
! 
! REAL(GPU_NODAL_PRECISION) :: del_phi, local_jout(2), local_src, local_srcm
! REAL(4) :: tau
! INTEGER :: ig, ipin, ipol, iz, izf, itrack
! 
! REAL(GPU_NODAL_PRECISION), SHARED :: track_phi(nPolar1D, cuDevice%sharedMemoryDim)
! 
! itrack = threadIdx%x
! ig = mod((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1, ng) + 1
! ipin = ((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1) / ng + 1
! 
! IF (ipin .GT. nxyc) RETURN
! 
! !--- Upward Sweep
! 
! track_phi(:, itrack) = PhiAngIn(:, ig, ipin, bottom)
! 
! DO iz = cuDevice%myzbf, cuDevice%myzef
!   local_jout = 0.0
!   DO ipol = 1, nPolar1D
!     local_jout(in) = local_jout(in) + wtsurf1D(ipol) * track_phi(ipol, itrack)
!   ENDDO
!   DO izf = cuDevice%fmSubrange(iz, 1), cuDevice%fmSubrange(iz, 2)
!     tau = xst(ig, ipin, izf) * cuDevice%hzSub(izf)
!     local_src = src(ig, ipin, izf); local_srcm = srcm(ig, ipin, izf)
!     DO ipol = 1, nPolar1D
!       del_phi = (local_src + Comp1D(ipol) * local_srcm - track_phi(ipol, itrack))                                   &
!                 * (1.0 - __EXPF(- tau * rcosv1D(ipol)))
!       track_phi(ipol, itrack) = track_phi(ipol, itrack) + del_phi
!     ENDDO
!   ENDDO
!   DO ipol = 1, nPolar1D
!     local_jout(out) = local_jout(out) + wtsurf1D(ipol) * track_phi(ipol, itrack)
!   ENDDO
!   Jout(ig, ipin, iz, bottom, in) = local_jout(in)
!   Jout(ig, ipin, iz, top, out) = local_jout(out)
! ENDDO
! 
! PhiAngOut(:, ig, ipin, top) = track_phi(:, itrack)
! 
! !--- Downward Sweep
! 
! track_phi(:, itrack) = PhiAngIn(:, ig, ipin, top)
! 
! DO iz = cuDevice%myzef, cuDevice%myzbf, -1
!   local_jout = 0.0
!   DO ipol = 1, nPolar1D
!     local_jout(in) = local_jout(in) + wtsurf1D(ipol) * track_phi(ipol, itrack)
!   ENDDO
!   DO izf = cuDevice%fmSubrange(iz, 2), cuDevice%fmSubrange(iz, 1), -1
!     tau = xst(ig, ipin, izf) * cuDevice%hzSub(izf)
!     local_src = src(ig, ipin, izf); local_srcm = srcm(ig, ipin, izf)
!     DO ipol = 1, nPolar1D
!       del_phi = (local_src - Comp1D(ipol) * local_srcm - track_phi(ipol, itrack))                                   &
!                 * (1.0 - __EXPF(- tau * rcosv1D(ipol)))
!       track_phi(ipol, itrack) = track_phi(ipol, itrack) + del_phi
!     ENDDO
!   ENDDO
!   DO ipol = 1, nPolar1D
!     local_jout(out) = local_jout(out) + wtsurf1D(ipol) * track_phi(ipol, itrack)
!   ENDDO
!   Jout(ig, ipin, iz, top, in) = local_jout(in)
!   Jout(ig, ipin, iz, bottom, out) = local_jout(out)
! ENDDO
! 
! PhiAngOut(:, ig, ipin, bottom) = track_phi(:, itrack)
! 
! END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuFlatRayTrace2nd(cuDevice, xst, src, srcm, PhiAngIn, PhiAngOut, Jout)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuDevice_Type) :: cuDevice
REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: xst(:, :, :), src(:, :, :), srcm(:, :, :)
REAL(GPU_NODAL_PRECISION), DEVICE :: PhiAngIn(:, :, :, :), PhiAngOut(:, :, :, :)
REAL, DEVICE :: Jout(:, :, cuDevice%myzbf :, :, :)

REAL(GPU_NODAL_PRECISION) :: del_phi, local_jout(2), local_src, local_srcm
REAL(4) :: tau
INTEGER :: ig, ipin, ipol, iz, izf, itrack

REAL(GPU_NODAL_PRECISION), SHARED :: track_phi(nPolar1D, cuDevice%sharedMemoryDim)

itrack = threadIdx%x
ipin = mod((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1, nxyc) + 1
ig = ((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1) / nxyc + 1

IF (ig .GT. ng) RETURN

!--- Upward Sweep

track_phi(:, itrack) = PhiAngIn(:, ipin, ig, bottom)

DO iz = cuDevice%myzbf, cuDevice%myzef
  local_jout = 0.0
  DO ipol = 1, nPolar1D
    local_jout(in) = local_jout(in) + wtsurf1D(ipol) * track_phi(ipol, itrack)
  ENDDO
  DO izf = cuDevice%fmSubrange(iz, 1), cuDevice%fmSubrange(iz, 2)
    tau = xst(ipin, ig, izf) * cuDevice%hzSub(izf)
    local_src = src(ipin, ig, izf); local_srcm = srcm(ipin, ig, izf)
    DO ipol = 1, nPolar1D
      del_phi = (local_src + Comp1D(ipol) * local_srcm - track_phi(ipol, itrack))                                   &
                * (1.0 - __EXPF(- tau * rcosv1D(ipol)))
      track_phi(ipol, itrack) = track_phi(ipol, itrack) + del_phi
    ENDDO
  ENDDO
  DO ipol = 1, nPolar1D
    local_jout(out) = local_jout(out) + wtsurf1D(ipol) * track_phi(ipol, itrack)
  ENDDO
  Jout(ig, ipin, iz, bottom, in) = local_jout(in)
  Jout(ig, ipin, iz, top, out) = local_jout(out)
ENDDO

PhiAngOut(:, ipin, ig, top) = track_phi(:, itrack)

!--- Downward Sweep

track_phi(:, itrack) = PhiAngIn(:, ipin, ig, top)

DO iz = cuDevice%myzef, cuDevice%myzbf, -1
  local_jout = 0.0
  DO ipol = 1, nPolar1D
    local_jout(in) = local_jout(in) + wtsurf1D(ipol) * track_phi(ipol, itrack)
  ENDDO
  DO izf = cuDevice%fmSubrange(iz, 2), cuDevice%fmSubrange(iz, 1), -1
    tau = xst(ipin, ig, izf) * cuDevice%hzSub(izf)
    local_src = src(ipin, ig, izf); local_srcm = srcm(ipin, ig, izf)
    DO ipol = 1, nPolar1D
      del_phi = (local_src - Comp1D(ipol) * local_srcm - track_phi(ipol, itrack))                                   &
                * (1.0 - __EXPF(- tau * rcosv1D(ipol)))
      track_phi(ipol, itrack) = track_phi(ipol, itrack) + del_phi
    ENDDO
  ENDDO
  DO ipol = 1, nPolar1D
    local_jout(out) = local_jout(out) + wtsurf1D(ipol) * track_phi(ipol, itrack)
  ENDDO
  Jout(ig, ipin, iz, top, in) = local_jout(in)
  Jout(ig, ipin, iz, bottom, out) = local_jout(out)
ENDDO

PhiAngOut(:, ipin, ig, bottom) = track_phi(:, itrack)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuSetFluxMoment(cuDevice, phis, phim, src, srcm, xst)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE

TYPE(cuDevice_Type) :: cuDevice
REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: src(:, :, :), srcm(:, :, :), xst(:, :, :)
REAL(GPU_NODAL_PRECISION), DEVICE :: phis(:, :, :), phim(:, :, :)

INTEGER :: ig, ipin, izf

ipin = mod((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1, nxyc) + 1
ig = ((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1) / nxyc + 1

IF (ig .GT. ng) RETURN

DO izf = 1, cuDevice%nzSub
  phis(ipin, ig, izf) = phis(ipin, ig, izf) / cuDevice%hzSub(izf)
  phim(ipin, ig, izf) = phim(ipin, ig, izf) / cuDevice%hzSub(izf)
  phis(ipin, ig, izf) = phis(ipin, ig, izf) / xst(ipin, ig, izf) + src(ipin, ig, izf)
  phim(ipin, ig, izf) = phim(ipin, ig, izf) / xst(ipin, ig, izf) + 1.0 / 3.0 * srcm(ipin, ig, izf)
ENDDO

END SUBROUTINE

END MODULE

MODULE CUDA_LINMOC

USE CUDA_MASTER
USE CUDA_UTIL
IMPLICIT NONE

! CONTAINS
! 
! ATTRIBUTES(GLOBAL) SUBROUTINE cuSetRTCoeff(cuDevice, xst, E1, E3, R1, R3)
! USE CUDA_CONST
! USE CUDAFOR
! IMPLICIT NONE
! 
! TYPE(cuDevice_Type) :: cuDevice
! REAL(GPU_NODAL_PRECISION), DEVICE :: xst(:, :, :)
! REAL(GPU_NODAL_PRECISION), DEVICE :: E1(:, :, :, :), E3(:, :, :, :), R1(:, :, :, :), R3(:, :, :, :)
! 
! INTEGER :: ipol, ig, ipin, izf
! REAL(GPU_PRECISION) :: s, tau, ex
! 
! ig = mod((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1, ng) + 1
! ipin = ((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1) / ng + 1
! 
! IF (ipin .GT. nxyc) RETURN
! 
! DO izf = 1, cuDevice%nzSub
!   s = cuDevice%hzSub(izf) * xst(ig, ipin, izf)
!   DO ipol = 1, nPolar1D
!     tau = s * rcosv1D(ipol)
!     ex = 1.0 - __EXPF(- tau)
!     E1(ig, ipin, ipol, izf) = ex
!     R1(ig, ipin, ipol, izf) = (1.0 + tau / 2.0 - (1.0 + 1.0 / tau) * ex) / tau
!     E3(ig, ipin, ipol, izf) = cosv1D(ipol) / 2.0 * (2.0 * (tau - ex) - tau * ex)
!     R3(ig, ipin, ipol, izf) = cosv1D(ipol) / 2.0 * (tau / 6.0 - 2.0 / tau - 2.0 + (1.0 + 1.0 / tau) * (1.0 + 2.0 / tau) * ex)
!   ENDDO
! ENDDO
! 
! END SUBROUTINE
! 
! ATTRIBUTES(GLOBAL) SUBROUTINE cuLinearRayTrace1st(cuDevice, phis, phim, src, E1, E3, R1, R3, PhiAngIn, PhiAngOut)
! USE CUDA_CONST
! USE CUDAFOR
! IMPLICIT NONE
! 
! TYPE(cuDevice_Type) :: cuDevice
! REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: src(:, :, :, 0 :), PhiAngIn(:, :, :, :)
! REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: E1(:, :, :, :), E3(:, :, :, :), R1(:, :, :, :), R3(:, :, :, :)
! REAL(GPU_NODAL_PRECISION), DEVICE :: phis(:, :, :, :), phim(:, :, :), PhiAngOut(:, :, :, :)
! 
! REAL(GPU_NODAL_PRECISION) :: del_phi, del_phim, local_phi, local_phim, local_src, local_srcm
! INTEGER :: ig, ipin, ipol, izf, itrack
! 
! REAL(GPU_NODAL_PRECISION), SHARED :: track_phi(nPolar1D, cuDevice%sharedMemoryDim)
! 
! itrack = threadIdx%x
! ig = mod((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1, ng) + 1
! ipin = ((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1) / ng + 1
! 
! IF (ipin .GT. nxyc) RETURN
! 
! !--- Upward Sweep
! 
! track_phi(:, itrack) = PhiAngIn(:, ig, ipin, bottom)
! 
! DO izf = 1, cuDevice%nzSub
!   local_src = src(ig, ipin, izf, 0); local_srcm = src(ig, ipin, izf, 1)
!   local_phi = 0.0; local_phim = 0.0
!   DO ipol = 1, nPolar1D
!     del_phi = (local_src - track_phi(ipol, itrack)) * E1(ig, ipin, ipol, izf)                                       &
!               + local_srcm * E3(ig, ipin, ipol, izf)
!     del_phim = track_phi(ipol, itrack) * 0.5 + (local_src - track_phi(ipol, itrack)) * R1(ig, ipin, ipol, izf)      &
!                + local_srcm * R3(ig, ipin, ipol, izf)
!     local_phi = local_phi - wtsurf1D(ipol) * del_phi
!     local_phim = local_phim + wt1D(ipol) * del_phim
!     track_phi(ipol, itrack) = track_phi(ipol, itrack) + del_phi
!   ENDDO
!   phis(ig, ipin, izf, 1) = local_phi
!   phim(ig, ipin, izf) = local_phim
! ENDDO
! 
! PhiAngOut(:, ig, ipin, top) = track_phi(:, itrack)
! 
! !--- Downward Sweep
! 
! track_phi(:, itrack) = PhiAngIn(:, ig, ipin, top)
! 
! DO izf = cuDevice%nzSub, 1, -1
!   local_src = src(ig, ipin, izf, 0); local_srcm = src(ig, ipin, izf, 1)
!   local_phi = 0.0; local_phim = 0.0
!   DO ipol = 1, nPolar1D
!     del_phi = (local_src - track_phi(ipol, itrack)) * E1(ig, ipin, ipol, izf)                                       &
!               - local_srcm * E3(ig, ipin, ipol, izf)
!     del_phim = track_phi(ipol, itrack) * 0.5 + (local_src - track_phi(ipol, itrack)) * R1(ig, ipin, ipol, izf)      &
!                - local_srcm * R3(ig, ipin, ipol, izf)
!     local_phi = local_phi - wtsurf1D(ipol) * del_phi
!     local_phim = local_phim - wt1D(ipol) * del_phim
!     track_phi(ipol, itrack) = track_phi(ipol, itrack) + del_phi
!   ENDDO
!   phis(ig, ipin, izf, 2) = local_phi
!   phim(ig, ipin, izf) = phim(ig, ipin, izf) + local_phim
! ENDDO
! 
! PhiAngOut(:, ig, ipin, bottom) = track_phi(:, itrack)
! 
! END SUBROUTINE
! 
! ATTRIBUTES(GLOBAL) SUBROUTINE cuLinearRayTrace2nd(cuDevice, phis, phim, src, E1, E3, R1, R3, PhiAngIn, PhiAngOut, Jout)
! USE CUDA_CONST
! USE CUDAFOR
! IMPLICIT NONE
! 
! TYPE(cuDevice_Type) :: cuDevice
! REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: src(:, :, :, 0 :), PhiAngIn(:, :, :, :)
! REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: E1(:, :, :, :), E3(:, :, :, :), R1(:, :, :, :), R3(:, :, :, :)
! REAL(GPU_NODAL_PRECISION), DEVICE :: phis(:, :, :, :), phim(:, :, :), PhiAngOut(:, :, :, :)
! REAL, DEVICE :: Jout(:, :, cuDevice%myzbf :, :, :)
! 
! REAL(GPU_NODAL_PRECISION) :: del_phi, del_phim, local_phi, local_phim, local_src, local_srcm, local_jout(3)
! INTEGER :: ig, ipin, ipol, iz, izf, itrack
! 
! REAL(GPU_NODAL_PRECISION), SHARED :: track_phi(nPolar1D, cuDevice%sharedMemoryDim)
! 
! itrack = threadIdx%x
! ig = mod((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1, ng) + 1
! ipin = ((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1) / ng + 1
! 
! IF (ipin .GT. nxyc) RETURN
! 
! !--- Upward Sweep
! 
! track_phi(:, itrack) = PhiAngIn(:, ig, ipin, bottom)
! 
! DO iz = cuDevice%myzbf, cuDevice%myzef
!   local_jout = 0.0
!   DO ipol = 1, nPolar1D
!     local_jout(in) = local_jout(in) + wtsurf1D(ipol) * track_phi(ipol, itrack)
!     local_jout(surf) = local_jout(surf) + wt1D(ipol) * track_phi(ipol, itrack)
!   ENDDO
!   Jout(ig, ipin, iz, bottom, in) = local_jout(in)
!   Jout(ig, ipin, iz, bottom, surf) = local_jout(surf)
!   DO izf = cuDevice%fmSubrange(iz, 1), cuDevice%fmSubrange(iz, 2)
!     local_src = src(ig, ipin, izf, 0); local_srcm = src(ig, ipin, izf, 1)
!     local_phi = 0.0; local_phim = 0.0
!     DO ipol = 1, nPolar1D
!       del_phi = (local_src - track_phi(ipol, itrack)) * E1(ig, ipin, ipol, izf)                                     &
!                 + local_srcm * E3(ig, ipin, ipol, izf)
!       del_phim = track_phi(ipol, itrack) * 0.5 + (local_src - track_phi(ipol, itrack)) * R1(ig, ipin, ipol, izf)    &
!                  + local_srcm * R3(ig, ipin, ipol, izf)
!       local_phi = local_phi - wtsurf1D(ipol) * del_phi
!       local_phim = local_phim + wt1D(ipol) * del_phim
!       track_phi(ipol, itrack) = track_phi(ipol, itrack) + del_phi
!     ENDDO
!     phis(ig, ipin, izf, 1) = local_phi
!     phim(ig, ipin, izf) = local_phim
!   ENDDO
!   local_jout = 0.0
!   DO ipol = 1, nPolar1D
!     local_jout(out) = local_jout(out) + wtsurf1D(ipol) * track_phi(ipol, itrack)
!     local_jout(surf) = local_jout(surf) + wt1D(ipol) * track_phi(ipol, itrack)
!   ENDDO
!   Jout(ig, ipin, iz, top, out) = local_jout(out)
!   Jout(ig, ipin, iz, top, surf) = local_jout(surf)
! ENDDO
! 
! PhiAngOut(:, ig, ipin, top) = track_phi(:, itrack)
! 
! !--- Downward Sweep
! 
! track_phi(:, itrack) = PhiAngIn(:, ig, ipin, top)
! 
! DO iz = cuDevice%myzef, cuDevice%myzbf, -1
!   local_jout = 0.0
!   DO ipol = 1, nPolar1D
!     local_jout(in) = local_jout(in) + wtsurf1D(ipol) * track_phi(ipol, itrack)
!     local_jout(surf) = local_jout(surf) + wt1D(ipol) * track_phi(ipol, itrack)
!   ENDDO
!   Jout(ig, ipin, iz, top, in) = local_jout(in)
!   Jout(ig, ipin, iz, top, surf) = local_jout(surf)
!   DO izf = cuDevice%fmSubrange(iz, 2), cuDevice%fmSubrange(iz, 1), -1
!     local_src = src(ig, ipin, izf, 0); local_srcm = src(ig, ipin, izf, 1)
!     local_phi = 0.0; local_phim = 0.0
!     DO ipol = 1, nPolar1D
!       del_phi = (local_src - track_phi(ipol, itrack)) * E1(ig, ipin, ipol, izf)                                     &
!                 - local_srcm * E3(ig, ipin, ipol, izf)
!       del_phim = track_phi(ipol, itrack) * 0.5 + (local_src - track_phi(ipol, itrack)) * R1(ig, ipin, ipol, izf)    &
!                  - local_srcm * R3(ig, ipin, ipol, izf)
!       local_phi = local_phi - wtsurf1D(ipol) * del_phi
!       local_phim = local_phim - wt1D(ipol) * del_phim
!       track_phi(ipol, itrack) = track_phi(ipol, itrack) + del_phi
!     ENDDO
!     phis(ig, ipin, izf, 2) = local_phi
!     phim(ig, ipin, izf) = phim(ig, ipin, izf) + local_phim
!   ENDDO
!   local_jout = 0.0
!   DO ipol = 1, nPolar1D
!     local_jout(out) = local_jout(out) + wtsurf1D(ipol) * track_phi(ipol, itrack)
!     local_jout(surf) = local_jout(surf) + wt1D(ipol) * track_phi(ipol, itrack)
!   ENDDO
!   Jout(ig, ipin, iz, bottom, out) = local_jout(out)
!   Jout(ig, ipin, iz, bottom, surf) = local_jout(surf)
! ENDDO
! 
! PhiAngOut(:, ig, ipin, bottom) = track_phi(:, itrack)
! 
! END SUBROUTINE
! 
! ATTRIBUTES(GLOBAL) SUBROUTINE cuSetFluxMoment(cuDevice, phis, phim, xst, src)
! USE CUDA_CONST
! USE CUDAFOR
! IMPLICIT NONE
! 
! TYPE(cuDevice_Type) :: cuDevice
! REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: phim(:, :, :), xst(:, :, :), src(:, :, :)
! REAL(GPU_NODAL_PRECISION), DEVICE :: phis(:, :, :, 0 :)
! 
! INTEGER :: ig, ipin, izf
! REAL(GPU_NODAL_PRECISION) :: myphis(2), myphim
! 
! ig = mod((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1, ng) + 1
! ipin = ((blockIdx%x - 1) * blockDim%x + threadIdx%x - 1) / ng + 1
! 
! IF (ipin .GT. nxyc) RETURN
! 
! DO izf = 1, cuDevice%nzSub
!   myphim = phim(ig, ipin, izf)
!   myphis = phis(ig, ipin, izf, :) / cuDevice%hzSub(izf)
!   myphis = myphis / xst(ig, ipin, izf)
!   phis(ig, ipin, izf, 0) = myphis(1) + myphis(2) + src(ig, ipin, izf)
!   phis(ig, ipin, izf, 1) = (- 6.0 * (myphis(1) - myphis(2)) + 12.0 * myphim) / cuDevice%hzSub(izf)
! ENDDO
! 
! END SUBROUTINE

END MODULE

#endif