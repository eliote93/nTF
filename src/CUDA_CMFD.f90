#include <CUDADEFINES.h>

#ifdef __PGI

MODULE CUDA_GCCMFD

USE CUDA_MASTER
USE CUDA_UTIL
USE CUDA_SYSTEM
USE CUDA_SOLVER
IMPLICIT NONE

CONTAINS

SUBROUTINE cuHomogenizeGcXS(cuCMFD, cuGcCMFD, cuDevice)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD, cuGcCMFD
TYPE(cuDevice_Type) :: cuDevice

TYPE(PinXS_Type), POINTER :: PinXS(:, :), GcPinXS(:, :)
REAL :: localphis(cuGeometry%ng)
INTEGER :: nxy
INTEGER :: ipin, ipin_map, iz, izf
INTEGER :: myzbf, myzef
INTEGER, POINTER :: pinMap(:), planeMap(:)

nxy = cuGeometry%nxyc
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
pinMap => cuGeometry%pinMap
planeMap => cuGeometry%planeMap
PinXS => cuCMFD%PinXS
GcPinXS => cuGcCMFD%PinXS

DO izf = myzbf, myzef
  iz = planeMap(izf)
  !$OMP PARALLEL PRIVATE(localphis, ipin_map)
  !$OMP DO SCHEDULE(GUIDED)
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    localphis = cuCMFD%h_phis8(:, ipin, izf)
    CALL cuHomogenizeCellGcXS(PinXS(ipin_map, iz), GcPinXS(ipin_map, izf), localphis)
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDDO

END SUBROUTINE

SUBROUTINE cuHomogenizeCellGcXS(PinXS, GcPinXS, phis)
USE PARAM
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type) :: PinXS, GcPinXS
REAL :: phis(:)

REAL :: RR(3), RRS(cuGeometry%ngc, cuGeometry%ngc)
REAL :: phisum, chisum
INTEGER :: ng, ngc
INTEGER :: igc, ig, igb, ige
INTEGER :: igs, igsb, igse, ig0

ng = cuGeometry%ng
ngc = cuGeometry%ngc

DO igc = 1, ngc
  igb = cuGeometry%GcStruct(1, igc); ige = cuGeometry%GcStruct(2, igc)
  RR = 0; phisum = 0; chisum = 0
  DO ig = igb, ige
    RR(1) = RR(1) + PinXS%XStr(ig) * phis(ig)
    RR(2) = RR(2) + PinXS%XSnf(ig) * phis(ig)
    RR(3) = RR(3) + PinXS%XSD(ig) * phis(ig)
    phisum = phisum + phis(ig)
    chisum = chisum + PinXS%Chi(ig)
  ENDDO
  RR = RR / phisum
  GcPinXS%XStr(igc) = RR(1)
  GcPinXS%XSnf(igc) = RR(2)
  GcPinXS%XSD(igc) = RR(3)
  GcPinXS%Chi(igc) = chisum
  GcPinXS%Phi(igc) = phisum
ENDDO

RRS = 0
DO ig = 1, ng
  igc = cuGeometry%GcStructInv(ig)
  igsb = PinXS%XSs(ig)%ib; igse = PinXS%XSs(ig)%ie
  ig0 = igc
  RRS(ig0, igc) = RRS(ig0, igc) + PinXS%XSs(ig)%WithInGroupScat * phis(ig)
  DO igs = igsb, igse
    ig0 = cuGeometry%GcStructInv(igs)
    RRS(ig0, igc) = RRS(ig0, igc) + PinXS%XSs(ig)%from(igs) * phis(igs)
  ENDDO
ENDDO

DO igc = 1, ngc
  RRS(igc, :) = RRS(igc, :) / GcPinXS%Phi(igc)
ENDDO

DO ig = 1, ngc
  igsb = GcPinXS%XSs(ig)%ib; igse = GcPinXS%XSs(ig)%ie
  DO igs = igsb, igse
    GcPinXS%XSs(ig)%from(igs) = RRS(igs, ig)
  ENDDO
  GcPinXS%XSs(ig)%WithInGroupScat = GcPinXS%XSs(ig)%from(ig)
  GcPinXS%XSs(ig)%from(ig) = 0.0
ENDDO

DO ig = 1, ngc
  GcPinXS%XSr(ig) = GcPinXS%XStr(ig) - GcPinXS%XSs(ig)%WithInGroupScat
ENDDO

END SUBROUTINE

SUBROUTINE cuSetRadialGcCoupling(cuCMFD, cuGcCMFD, cuDevice)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD, cuGcCMFD
TYPE(cuDevice_Type) :: cuDevice

TYPE(PinXS_Type), POINTER :: PinXS(:, :), GcPinXS(:, :)
TYPE(superPin_Type), POINTER :: Pin(:)
INTEGER, POINTER :: pinMap(:), pinMapRev(:), planeMap(:)
INTEGER :: ng, ngc, nxy
INTEGER :: ig, igc, ipin, ipin_map, ineighpin, iz, izf, ibd, inbd
INTEGER :: myzbf, myzef
REAL :: Dtil, Dhat, atil, myphi, neighphi, mybeta, neighbeta, jfdm, surfphifdm, smy
REAL, POINTER :: Jnet(:, :, :, :)

Pin => cuGeometry%superPin
ng = cuGeometry%ng
ngc = cuGeometry%ngc
nxy = cuGeometry%nxyc
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
pinMap => cuGeometry%pinMap
pinMapRev => cuGeometry%pinMapRev
planeMap => cuGeometry%planeMap
PinXS => cuCMFD%PinXS
GcPinXS => cuGcCMFD%PinXS

ALLOCATE(Jnet(4, ngc, nxy, myzbf : myzef)); Jnet = 0.0

!--- Condense Currents

DO izf = myzbf, myzef
  iz = planeMap(izf)
  !$OMP PARALLEL PRIVATE(ipin_map, ineighpin, myphi, neighphi, Dtil, Dhat, jfdm)
  !$OMP DO SCHEDULE(GUIDED)
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    DO igc = 1, ngc
      DO ig = cuGeometry%GcStruct(1, igc), cuGeometry%GcStruct(2, igc)
        myphi = cuCMFD%h_phis8(ig, ipin, izf)
        DO ibd = 1, 4
          ineighpin = Pin(ipin_map)%Neighidx(ibd)
          ineighpin = pinMapRev(ineighpin)
          IF (ineighpin .LE. 0) THEN
            neighphi = 0.0
          ELSE
            neighphi = cuCMFD%h_phis8(ig, ineighpin, izf)
          ENDIF
          Dtil = PinXS(ipin_map, iz)%Dtil(ibd, ig)
          Dhat = PinXS(ipin_map, iz)%Dhat(ibd, ig)
          jfdm = - Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
          Jnet(ibd, igc, ipin, izf) = Jnet(ibd, igc, ipin, izf) + jfdm
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDDO

!--- Compute Coupling Coefficients

!$OMP PARALLEL PRIVATE(ipin_map, ineighpin, inbd, myphi, neighphi, mybeta, neighbeta, Dtil, Dhat, jfdm, smy)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO igc = 1, ngc
  DO izf = myzbf, myzef
    DO ipin = 1, nxy
      ipin_map = pinMap(ipin)
      DO ibd = 1, 4
        ineighpin = Pin(ipin_map)%NeighIdx(ibd)
        smy = Pin(ipin_map)%BdLength(ibd)
        myphi = GcPinXS(ipin_map, izf)%Phi(igc)
        mybeta = GcPinXS(ipin_map, izf)%XSD(igc) / Pin(ipin_map)%Center2SurfaceL(ibd)
        IF (ineighpin .GT. 0) THEN
          inbd = Pin(ineighpin)%NeighSurfIdx(ibd)
          neighphi = GcPinXS(ineighpin, izf)%Phi(igc)
          neighbeta = GcPinXS(ineighpin, izf)%XSD(igc) / Pin(ineighpin)%Center2SurfaceL(inbd)
          Dtil = mybeta * neighbeta / (mybeta + neighbeta) * smy
          jfdm = - Dtil * (neighphi - myphi)
          Dhat = - (Jnet(ibd, igc, ipin, izf) - jfdm) / (myphi + neighphi)
        ELSE
          IF (ineighpin .EQ. Void) THEN
            neighbeta = 0.5; neighphi = 0.0
          ELSEIF (ineighpin .EQ. Reflective) THEN
            neighbeta = 0.0; neighphi = myphi
          ENDIF
          Dtil = mybeta * neighbeta / (mybeta + neighbeta) * smy
          jfdm = - Dtil * (neighphi - myphi)
          Dhat = - (Jnet(ibd, igc, ipin, izf) - jfdm) / (myphi + neighphi)
        ENDIF
        GcPinXS(ipin_map, izf)%Dtil(ibd, igc) = Dtil
        GcPinXS(ipin_map, izf)%Dhat(ibd, igc) = Dhat
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

DEALLOCATE(Jnet)

END SUBROUTINE

SUBROUTINE cuGcReconstruction(cuCMFD, cuGcCMFD, cuDevice)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD, cuGcCMFD
TYPE(cuDevice_Type) :: cuDevice

TYPE(PinXS_Type), POINTER :: GcPinXS(:, :)
INTEGER :: ng, ngc, nxy, nzCMFD
INTEGER :: ig, igc, ipin, ipin_map, iz, ierr
INTEGER :: myzbf, myzef
INTEGER, POINTER :: pinMap(:)
REAL :: fmult

ng = cuGeometry%ng
ngc = cuGeometry%ngc
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
pinMap => cuGeometry%pinMap
GcPinXS => cuGcCMFD%PinXS

!$OMP PARALLEL PRIVATE(ipin_map, igc, fmult)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO iz = myzbf, myzef
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    DO ig = 1, ng
      igc = cuGeometry%GcStructInv(ig)
      fmult = cuGcCMFD%h_phis8(igc, ipin, iz) / GcPinXS(ipin_map, iz)%Phi(igc)
      cuCMFD%h_phis8(ig, ipin, iz) = cuCMFD%h_phis8(ig, ipin, iz) * fmult
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

! ierr = cublasDcopy_v2(cuDevice%myblasHandle, nxy * nzCMFD, cuGcCMFD%psi8, 1, cuCMFD%psi8, 1)

END SUBROUTINE

END MODULE

MODULE CUDA_CMFD

USE CUDA_MASTER
USE CUDA_UTIL
USE CUDA_SYSTEM
USE CUDA_SOLVER
IMPLICIT NONE

LOGICAL :: lFirstCMFD = .TRUE.

CONTAINS

SUBROUTINE cuSetCMFDPhis(cuCMFD, cuDevice, lCopy)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
LOGICAL :: lCopy

INTEGER :: ng, nxy
INTEGER :: ig, ipin, ipin_map, iz, izf
INTEGER :: myzb, myze, myzbf, myzef
INTEGER, POINTER :: pinMap(:), planeMap(:), fmRange(:, :)
REAL :: fmult(cuCMFD%ng)

ng = cuCMFD%ng
nxy = cuGeometry%nxyc
myzb = cuDevice%myzb
myze = cuDevice%myze
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
pinMap => cuGeometry%pinMap
planeMap => cuCMFD%planeMap
fmRange => cuGeometry%fmRange

IF (lCopy) THEN
  DO izf = myzbf, myzef
    iz = planeMap(izf)
    !$OMP PARALLEL PRIVATE(ipin_map)
    !$OMP DO SCHEDULE(GUIDED)
    DO ipin = 1, nxy
      ipin_map = pinMap(ipin)
      DO ig = 1, ng
        cuCMFD%h_phis8(ig, ipin, izf) = cuCMFD%PinXS(ipin_map, iz)%Phi(ig)
      ENDDO
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL
  ENDDO
ELSE
  !$OMP PARALLEL PRIVATE(ipin_map, fmult)
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
  DO iz = myzb, myze
    DO ipin = 1, nxy
      ipin_map = pinMap(ipin)
      fmult = cuCMFD%PinXS(ipin_map, iz)%Phi / cuCMFD%h_phic8(:, ipin, iz)
      DO izf = fmRange(iz, 1), fmRange(iz, 2)
        DO ig = 1, ng
          cuCMFD%h_phis8(ig, ipin, izf) = cuCMFD%h_phis8(ig, ipin, izf) * fmult(ig)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDIF

END SUBROUTINE

SUBROUTINE cuSetMOCPhis(Core, cuCMFD, cuDevice, phis, phim, phic, lScat1)
USE TYPEDEF,        ONLY : CoreInfo_Type,       Pin_Type,           Cell_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
REAL, POINTER :: phis(:, :, :), phim(:, :, :, :), phic(:, :, :)
LOGICAL :: lScat1

TYPE(superPin_Type), POINTER :: superPin(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
INTEGER :: ng, nxy, nLocalFsr, FsrIdxSt
INTEGER :: i, j, ixy, ixy_map, ifsr, icel, ig, ipin, iz, izf, ierr
INTEGER :: myzb, myze
INTEGER, POINTER :: pinMap(:), fmRange(:, :)
REAL :: phisum(cuCMFD%ng), fmult(cuCMFD%ng)
REAL, POINTER :: hz(:), hzfm(:)

Pin => Core%Pin
Cell => Core%CellInfo
superPin => cuGeometry%superPin
ng = cuCMFD%ng
nxy = cuGeometry%nxyc
myzb = cuDevice%myzb
myze = cuDevice%myze
pinMap => cuGeometry%pinMap
fmRange => cuGeometry%fmRange
hz => cuGeometry%hz
hzfm => cuGeometry%hzfm

!$OMP PARALLEL PRIVATE(ixy_map, ifsr, ipin, icel, FsrIdxSt, nLocalFsr, phisum, fmult)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO iz = myzb, myze
  DO ixy = 1, nxy
    phisum = 0.0
    ixy_map = pinMap(ixy)
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      phisum = phisum + cuCMFD%h_phis8(:, ixy, izf) * (hzfm(izf) / hz(iz))
    ENDDO
    fmult = phisum / cuCMFD%PinXS(ixy_map, iz)%Phi
    DO j = 1, superPin(ixy_map)%nxy
      ipin = superPin(ixy_map)%pin(j)
      FsrIdxSt = Pin(ipin)%FsrIdxSt
      icel = Pin(ipin)%Cell(iz)
      nLocalFsr = Cell(icel)%nFsr
      DO i = 1, nLocalFsr
        ifsr = FsrIdxSt + i - 1
        phis(ifsr, iz, :) = phis(ifsr, iz, :) * fmult
      ENDDO
      IF (lScat1) THEN
        DO i = 1, nLocalFsr
          ifsr = FsrIdxSt + i - 1
          DO ig = 1, ng
            phim(:, ig, ifsr, iz) = phim(:, ig, ifsr, iz) * fmult(ig)
          ENDDO
        ENDDO
      ENDIF
      phic(ipin, iz, :) = phisum
    ENDDO
    cuCMFD%h_phic8(:, ixy, iz) = phisum
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE cuSetAxialSrc(cuCMFD, cuDevice, AxSrc, AxPXS, phic)
USE CNTL,           ONLY : nTracerCntl
IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
REAL, POINTER :: AxSrc(:, :, :), AxPXS(:, :, :), phic(:, :, :)

TYPE(superPin_Type), POINTER :: Pin(:)
INTEGER :: i, ig, iz, izf, ipin, ixy, ixy_map
INTEGER :: ng, nxy
INTEGER :: myzb, myze
INTEGER, POINTER :: pinMap(:), fmRange(:, :)
REAL :: myphi, neighphi, Dtil, Dhat, Jnet
REAL, POINTER :: hz(:)

Pin => cuGeometry%superPin
ng = cuCMFD%ng
nxy = cuGeometry%nxyc
myzb = cuDevice%myzb
myze = cuDevice%myze
pinMap => cuGeometry%pinMap
fmRange => cuGeometry%fmRange
hz => cuGeometry%hz

!$OMP PARALLEL PRIVATE(ixy_map, ipin, myphi, neighphi, Dtil, Dhat, Jnet)

!--- Axial Source from Bottom
DO iz = myzb, myze
  izf = fmRange(iz, bottom)
  !$OMP DO SCHEDULE(GUIDED)
  DO ixy = 1, nxy
    ixy_map = pinMap(ixy)
    DO ig = 1, ng
      myphi = cuCMFD%h_phis8(ig, ixy, izf)
      IF (iz .EQ. myzb) THEN
        neighphi = cuCMFD%h_neighphis8(ig, ixy, bottom)
      ELSE
        neighphi = cuCMFD%h_phis8(ig, ixy, izf - 1)
      ENDIF
      Dtil = cuCMFD%AxDtil(bottom, ig, ixy, izf)
      Dhat = cuCMFD%AxDhat(bottom, ig, ixy, izf)
      Jnet = - Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
      DO i = 1, Pin(ixy_map)%nxy
        ipin = Pin(ixy_map)%pin(i)
        AxSrc(ipin, iz, ig) = Jnet
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
ENDDO

!--- Axial Source from Top
DO iz = myzb, myze
  izf = fmRange(iz, top)
  !$OMP DO SCHEDULE(GUIDED)
  DO ixy = 1, nxy
    ixy_map = pinMap(ixy)
    DO ig = 1, ng
      myphi = cuCMFD%h_phis8(ig, ixy, izf)
      IF (iz .EQ. myze) THEN
        neighphi = cuCMFD%h_neighphis8(ig, ixy, top)
      ELSE
        neighphi = cuCMFD%h_phis8(ig, ixy, izf + 1)
      ENDIF
      Dtil = cuCMFD%AxDtil(top, ig, ixy, izf)
      Dhat = cuCMFD%AxDhat(top, ig, ixy, izf)
      Jnet = - Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
      DO i = 1, Pin(ixy_map)%nxy
        ipin = Pin(ixy_map)%pin(i)
        AxSrc(ipin, iz, ig) = AxSrc(ipin, iz, ig) + Jnet
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
ENDDO

!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(ipin)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO iz = myzb, myze
    DO ixy = 1, nxy
      DO i = 1, Pin(ixy)%nxy
        ipin = Pin(ixy)%pin(i)
        AxSrc(ipin, iz, ig) = AxSrc(ipin, iz, ig) / hz(iz)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

AxPXS(:, myzb : myze, :) = 0.0

IF (nTracerCntl%LkgSplitLv .EQ. 0) THEN
  !$OMP PARALLEL PRIVATE(ipin)
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
  DO ig = 1, ng
    DO iz = myzb, myze
      DO ixy = 1, nxy
        DO i = 1, Pin(ixy)%nxy
          ipin = Pin(ixy)%pin(i)
          IF (AxSrc(ipin, iz, ig) .LT. 0.0) CYCLE
          IF (phic(ipin, iz, ig) .LT. 0.0) CYCLE
          AxPXS(ipin, iz, ig) = AxSrc(ipin, iz, ig) / phic(ipin, iz, ig)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDIF

END SUBROUTINE

SUBROUTINE cuCopyFlux(cuCMFD, cuDevice, dir)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
INTEGER :: dir

INTEGER :: ng, nxy, nzCMFD, memSize
INTEGER :: ierr

ng = cuCMFD%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD

memSize = ng * nxy * nzCMFD

SELECT CASE (dir)

CASE (1)   !--- Host to Device

ierr = cudaMemcpy(cuCMFD%phis8, cuCMFD%h_phis8, memSize, cudaMemcpyHostToDevice)

CASE (2)   !--- Device to Host

ierr = cudaMemcpy(cuCMFD%h_phis8, cuCMFD%phis8, memSize, cudaMemcpyDeviceToHost)

END SELECT

END SUBROUTINE

SUBROUTINE cuReorderFlux(cuCMFD, cuDevice, dir)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
INTEGER :: dir

INTEGER :: n, ng, nxy, nzCMFD
INTEGER :: ig, ixy, ipin, iz, color
INTEGER :: myzbf, myzef, rbBeg(2), rbEnd(2)
INTEGER, POINTER :: pinMap(:), pinMapRev(:), pinMapRB(:), pinMapRevRB(:, :), planeMapRB(:)
REAL, POINTER :: phisRB(:, :)

ng = cuCMFD%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
rbBeg = cuDevice%rbBeg
rbEnd = cuDevice%rbEnd
pinMap => cuGeometry%pinMap
pinMapRev => cuGeometry%pinMapRev
pinMapRB => cuDevice%pinMapRB
pinMapRevRB => cuDevice%pinMapRevRB
planeMapRB => cuDevice%planeMapRB

n = ng * nxy * nzCMFD

ALLOCATE(phisRB(ng, nxy * nzCMFD))

SELECT CASE (dir)

CASE (1)   !--- Natural to Red-Black

!$OMP PARALLEL PRIVATE(ipin)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO iz = myzbf, myzef
  DO ixy = 1, nxy
    ipin = pinMap(ixy)
    ipin = pinMapRevRB(ipin, iz)
    DO ig = 1, ng
      phisRB(ig, ipin) = cuCMFD%h_phis8(ig, ixy, iz)
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

CALL ArrayCopy(phisRB, cuCMFD%h_phis8, n)

CASE (2)   !--- Red-Black to Natural

CALL ArrayCopy(cuCMFD%h_phis8, phisRB, n)

!$OMP PARALLEL PRIVATE(iz, ipin)
DO color = 1, 2
  !$OMP DO SCHEDULE(GUIDED)
  DO ixy = rbBeg(color), rbEnd(color)
    iz = planeMapRB(ixy)
    ipin = pinMapRB(ixy)
    ipin = pinMapRev(ipin)
    DO ig = 1, ng
      cuCMFD%h_phis8(ig, ipin, iz) = phisRB(ig, ixy)
    ENDDO
  ENDDO
  !$OMP END DO
ENDDO
!$OMP END PARALLEL

END SELECT

DEALLOCATE(phisRB)

CONTAINS

SUBROUTINE ArrayCopy(A, B, n)

IMPLICIT NONE

REAL :: A(*), B(*)
INTEGER :: n

B(1 : n) = A(1 : n)

END SUBROUTINE

END SUBROUTINE

SUBROUTINE cuGetNeighborFlux(cuCMFD, cuDevice)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice

INTEGER :: n, ng, nxy
INTEGER :: myzbf, myzef
INTEGER :: ierr
REAL(8), ALLOCATABLE :: myphis8(:, :, :)

ng = cuCMFD%ng
nxy = cuGeometry%nxyc
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
n = ng * nxy

! ALLOCATE(myphis8(ng, nxy, 2))

! myphis8(:, :, bottom) = cuCMFD%h_phis8(:, :, myzbf)
! myphis8(:, :, top) = cuCMFD%h_phis8(:, :, myzef)
! 
! CALL cuGetNeighbor(n, myphis8, cuCMFD%h_neighphis8, MPI_CUDA_COMM)

CALL cuGetNeighbor(n, cuCMFD%h_phis8(:, :, myzbf), cuCMFD%h_neighphis8(:, :, bottom), MPI_CUDA_COMM, bottom)
CALL cuGetNeighbor(n, cuCMFD%h_phis8(:, :, myzef), cuCMFD%h_neighphis8(:, :, top), MPI_CUDA_COMM, top)

IF (cuDevice%lBottom) THEN
  IF (cuGeometry%AxBC(bottom) .EQ. Void) cuCMFD%h_neighphis8(:, :, bottom) = 0.0
  IF (cuGeometry%AxBC(bottom) .EQ. Reflective) cuCMFD%h_neighphis8(:, :, bottom) = cuCMFD%h_phis8(:, :, myzbf)
ENDIF

IF (cuDevice%lTop) THEN
  IF (cuGeometry%AxBC(top) .EQ. Void) cuCMFD%h_neighphis8(:, :, top) = 0.0
  IF (cuGeometry%AxBC(top) .EQ. Reflective) cuCMFD%h_neighphis8(:, :, top) = cuCMFD%h_phis8(:, :, myzef)
ENDIF

! DEALLOCATE(myphis8)

END SUBROUTINE

SUBROUTINE cuCMFDSrcUpdt(cuCMFD, cuDevice, eigv, lUpscat)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
REAL :: eigv
LOGICAL :: lUpscat

INTEGER :: n

n = cuCMFD%ng * cuGeometry%nxyc * cuDevice%nzCMFD

CALL cuInitArray(n, cuCMFD%src8, cuDevice%myStream)
CALL cuCMFDScatSrcUpdt(cuCMFD, cuDevice, lUpscat)
CALL cuCMFDFisSrcUpdt(cuCMFD, cuDevice, eigv)
CALL cuCMFDDcplSrcUpdt(cuCMFD, cuDevice)
  
END SUBROUTINE

SUBROUTINE cuCMFDScatSrcUpdt(cuCMFD, cuDevice, lUpscat)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
LOGICAL :: lUpscat

REAL(8), POINTER, DEVICE :: csrVal(:)
INTEGER, POINTER, DEVICE :: csrRowPtr(:), csrColIdx(:)
INTEGER :: nr, nc, nnz
INTEGER :: ierr

IF (cuCntl%lNatural) RETURN
IF (.NOT. lUpscat) RETURN

csrVal => cuCMFD%S%d_csrVal
csrRowPtr => cuCMFD%S%d_csrRowPtr
csrColIdx => cuCMFD%S%d_csrColIdx
nr = cuCMFD%S%nr
nc = cuCMFD%S%nc
nnz = cuCMFD%S%nnz

ierr = cusparseDcsrmv(cuDevice%mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, 1.0_8,                &
                      cuCMFD%S%descr, csrVal, csrRowPtr, csrColIdx, cuCMFD%phis8, 0.0_8, cuCMFD%src8)

END SUBROUTINE

SUBROUTINE cuCMFDFisSrcUpdt(cuCMFD, cuDevice, eigv)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
REAL :: eigv

REAL(8), POINTER, DEVICE :: csrVal(:)
INTEGER, POINTER, DEVICE :: csrRowPtr(:), csrColIdx(:)
INTEGER :: nr, nc, nnz
INTEGER :: ierr
REAL(8) :: reigv

IF (.NOT. cuDevice%lFuel) RETURN

csrVal => cuCMFD%Chi%d_csrVal
csrRowPtr => cuCMFD%Chi%d_csrRowPtr
csrColIdx => cuCMFD%Chi%d_csrColIdx
nr = cuCMFD%Chi%nr
nc = cuCMFD%Chi%nc
nnz = cuCMFD%Chi%nnz
reigv = 1.0 / eigv

ierr = cusparseDcsrmv(cuDevice%mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, reigv,                &
                      cuCMFD%Chi%descr, csrVal, csrRowPtr, csrColIdx, cuCMFD%psi8, 1.0_8, cuCMFD%src8)

END SUBROUTINE

SUBROUTINE cuCMFDDcplSrcUpdt(cuCMFD, cuDevice)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice

INTEGER :: n, nb(2), nt(2), ng, nxy
INTEGER :: myzbf, myzef
INTEGER :: bottomRange(2), bottomRangeRB(2, 2), topRange(2), topRangeRB(2, 2)
INTEGER :: ierr
REAL(8), ALLOCATABLE, DEVICE :: myphis8(:, :), neighphis8(:, :), dcplSrc8(:, :)

IF (.NOT. cuCntl%lDcpl .OR. .NOT. cuCntl%lMulti) RETURN

ng = cuCMFD%ng
nxy = cuGeometry%nxyc
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
bottomRange = cuDevice%bottomRange
bottomRangeRB = cuDevice%rbRange(:, myzbf, :)
topRange = cuDevice%topRange
topRangeRB = cuDevice%rbRange(:, myzef, :)

n = ng * nxy
nb = ng * (bottomRangeRB(2, :) - bottomRangeRB(1, :) + 1)
nt = ng * (topRangeRB(2, :) - topRangeRB(1, :) + 1)

ALLOCATE(myphis8(n, 2))
ALLOCATE(neighphis8(n, 2))
ALLOCATE(dcplSrc8(n, 2))

CALL cuInitArray(n * 2, neighphis8, cuDevice%myStream)

IF (cuCntl%lNatural) THEN
  ierr = cudaMemcpy(myphis8(:, bottom), cuCMFD%phis8(:, bottomRange(1) : bottomRange(2)), n,                        &
                    cudaMemcpyDeviceToDevice)
  ierr = cudaMemcpy(myphis8(:, top), cuCMFD%phis8(:, topRange(1) : topRange(2)), n,                                 &
                    cudaMemcpyDeviceToDevice)
ELSE
  ierr = cudaMemcpy(myphis8(1 : nb(black), bottom),                                                                 &
                    cuCMFD%phis8(:, bottomRangeRB(1, black) : bottomRangeRB(2, black)),                             &
                    nb(black), cudaMemcpyDeviceToDevice)
  ierr = cudaMemcpy(myphis8(nb(black) + 1 : n, bottom),                                                             &
                    cuCMFD%phis8(:, bottomRangeRB(1, red) : bottomRangeRB(2, red)),                                 &
                    nb(red), cudaMemcpyDeviceToDevice)
  ierr = cudaMemcpy(myphis8(1 : nt(black), top),                                                                    &
                    cuCMFD%phis8(:, topRangeRB(1, black) : topRangeRB(2, black)),                                   &
                    nt(black), cudaMemcpyDeviceToDevice)
  ierr = cudaMemcpy(myphis8(nt(black) + 1 : n, top),                                                                &
                    cuCMFD%phis8(:, topRangeRB(1, red) : topRangeRB(2, red)),                                       &
                    nt(red), cudaMemcpyDeviceToDevice)     
ENDIF
  
CALL InitMPIComm()
CALL MPIComm(n, n, myphis8(:, bottom), neighphis8(:, bottom), bottom, MPI_CUDA_COMM)
CALL MPIComm(n, n, myphis8(:, top), neighphis8(:, top), top, MPI_CUDA_COMM)
CALL FinalizeMPIComm()

CALL AddDcplSrc(bottom)
CALL AddDcplSrc(top)

ierr = cudaStreamSynchronize(cuDevice%myStream)

DEALLOCATE(myphis8)
DEALLOCATE(neighphis8)
DEALLOCATE(dcplSrc8)

CONTAINS

SUBROUTINE AddDcplSrc(dir)

IMPLICIT NONE

INTEGER :: dir

SELECT CASE (dir)

CASE (bottom)

  CALL cuVectorOp('*', n, cuCMFD%offDiag8(:, :, bottom), neighphis8(:, bottom),                                     &
                  dcplSrc8(:, bottom), cuDevice%myStream)
  IF (cuCntl%lNatural) THEN
    CALL cuVectorOp('-', n, cuCMFD%src8(:, bottomRange(1) : bottomRange(2)), dcplSrc8(:, bottom),                   &
                    cuCMFD%src8(:, bottomRange(1) : bottomRange(2)), cuDevice%myStream)
  ELSE
    CALL cuVectorOp('-', nb(red),                                                                                   &
                    cuCMFD%src8(:, bottomRangeRB(1, red) : bottomRangeRB(2, red)),                                  &
                    dcplSrc8(1 : nb(red), bottom),                                                                  &
                    cuCMFD%src8(:, bottomRangeRB(1, red) : bottomRangeRB(2, red)),                                  &
                    cuDevice%myStream)
    CALL cuVectorOp('-', nb(black),                                                                                 &
                    cuCMFD%src8(:, bottomRangeRB(1, black) : bottomRangeRB(2, black)),                              &
                    dcplSrc8(nb(red) + 1 : n, bottom),                                                              &
                    cuCMFD%src8(:, bottomRangeRB(1, black) : bottomRangeRB(2, black)),                              &
                    cuDevice%myStream)
  ENDIF

CASE (top)

  CALL cuVectorOp('*', n, cuCMFD%offDiag8(:, :, top), neighphis8(:, top),                                           &
                  dcplSrc8(:, top), cuDevice%myStream)
  IF (cuCntl%lNatural) THEN
    CALL cuVectorOp('-', n, cuCMFD%src8(:, topRange(1) : topRange(2)), dcplSrc8(:, top),                            &
                    cuCMFD%src8(:, topRange(1) : topRange(2)), cuDevice%myStream)
  ELSE
    CALL cuVectorOp('-', nt(red),                                                                                   &
                    cuCMFD%src8(:, topRangeRB(1, red) : topRangeRB(2, red)),                                        &
                    dcplSrc8(1 : nt(red), top),                                                                     &
                    cuCMFD%src8(:, topRangeRB(1, red) : topRangeRB(2, red)),                                        &
                    cuDevice%myStream)
    CALL cuVectorOp('-', nt(black),                                                                                 &
                    cuCMFD%src8(:, topRangeRB(1, black) : topRangeRB(2, black)),                                    &
                    dcplSrc8(nt(red) + 1 : n, top),                                                                 &
                    cuCMFD%src8(:, topRangeRB(1, black) : topRangeRB(2, black)),                                    &
                    cuDevice%myStream)
  ENDIF

END SELECT

END SUBROUTINE

END SUBROUTINE

SUBROUTINE cuCMFDPsiUpdt(cuCMFD, cuDevice)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice

REAL(8), POINTER, DEVICE :: csrVal(:)
INTEGER, POINTER, DEVICE :: csrRowPtr(:), csrColIdx(:)
INTEGER :: nr, nc, nnz
INTEGER :: ierr

IF (.NOT. cuDevice%lFuel) RETURN

csrVal => cuCMFD%F%d_csrVal
csrRowPtr => cuCMFD%F%d_csrRowPtr
csrColIdx => cuCMFD%F%d_csrColIdx
nr = cuCMFD%F%nr
nc = cuCMFD%F%nc
nnz = cuCMFD%F%nnz

ierr = cublasDcopy_v2(cuDevice%myblasHandle, nr, cuCMFD%psi8, 1, cuCMFD%psid8, 1)
ierr = cusparseDcsrmv(cuDevice%mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, 1.0_8,                &
                      cuCMFD%F%descr, csrVal, csrRowPtr, csrColIdx, cuCMFD%phis8, 0.0_8, cuCMFD%psi8)

END SUBROUTINE

SUBROUTINE cuCMFDEigUpdt(cuCMFD, cuDevice, eigv)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
REAL :: eigv

REAL(8) :: psipsi, psipsid
INTEGER :: n, nxy, nzCMFD

nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
n = nxy * nzCMFD

psipsi = dotMulti(cuCMFD%psi8, cuCMFD%psi8, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
psipsid = dotMulti(cuCMFD%psi8, cuCMFD%psid8, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)

eigv = eigv * psipsi / psipsid

END SUBROUTINE

FUNCTION cuCMFDResidual(cuCMFD, cuDevice) RESULT(res)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice

INTEGER :: n, ng, nxy, nzCMFD
INTEGER :: bottomRange(2), topRange(2), bottomRangeRB(2, 2), topRangeRB(2, 2)
REAL(8), ALLOCATABLE, DEVICE :: r(:)
REAL :: res

ng = cuCMFD%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
bottomRange = cuCMFD%bottomRange
bottomRangeRB = cuCMFD%bottomRangeRB
topRange = cuCMFD%topRange
topRangeRB = cuCMFD%topRangeRB
n = ng * nxy * nzCMFD

ALLOCATE(r(n))

IF (cuCntl%lNatural) THEN
  CALL cuMatMul3D(cuCMFD%M, cuCMFD%offDiag8, cuCMFD%phis8, r, bottomRange, topRange,                                &
                  cuDevice%mySparseHandle, cuDevice%myStream, MPI_CUDA_COMM)
ELSE
  CALL cuMatMul3D(cuCMFD%M, cuCMFD%offDiag8, cuCMFD%phis8, r, bottomRangeRB, topRangeRB,                            &
                  cuDevice%mySparseHandle, cuDevice%myStream, MPI_CUDA_COMM)
ENDIF
CALL cuVectorOp('-', n, cuCMFD%src8, r, r, cuDevice%myStream)
res = normMulti(r, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
res = res / normMulti(cuCMFD%src8, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)

DEALLOCATE(r)

END FUNCTION

FUNCTION cuCMFDPsiErr(cuCMFD, cuDevice) RESULT(err)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice

INTEGER :: n, nxy, nzCMFD
REAL(8), ALLOCATABLE, DEVICE :: e(:)
REAL :: err

nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
n = nxy * nzCMFD

ALLOCATE(e(n))

CALL cuVectorOp('-', n, cuCMFD%psi8, cuCMFD%psid8, e, cuDevice%myStream)
err = normMulti(e, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
err = err / normMulti(cuCMFD%psi8, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)

DEALLOCATE(e)

END FUNCTION

END MODULE

SUBROUTINE CUDAReorderPinXS(CoreInfo, CmInfo)
USE TYPEDEF,        ONLY : CoreInfo_Type,       CmInfo_Type,         PinXS_Type
USE CMFD_COMMON
USE CUDA_MASTER
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(CmInfo_Type) :: CmInfo

TYPE(PinXS_Type), POINTER :: PinXS(:, :)
TYPE(superPin_Type), POINTER :: Pin(:)
REAL, POINTER :: PinVol(:, :)
INTEGER :: ng, nxy, myzb, myze
INTEGER :: i, ipin, iFuelPin, ixy, iz

PinXS => CmInfo%PinXS
Pin => cuGeometry%superPin
PinVol => CoreInfo%PinVol
ng = cuGeometry%ng
nxy = cuGeometry%nxyc
myzb = cuDevice%myzb
myze = cuDevice%myze

DO iz = myzb, myze
  DO ixy = 1, nxy
    DO i = 1, Pin(ixy)%nxy
      ipin = Pin(ixy)%pin(i)
      CALL CopyPinXS(cuCMFD%PinXS(ixy, iz), PinXS(ipin, iz), ng)
    ENDDO
  ENDDO
ENDDO

DO iz = myzb, myze
  DO ixy = 1, nxy
    IF (.NOT. Pin(ixy)%lFuel(iz)) CYCLE
    iFuelPin = Pin(ixy)%iFuelPin
    PinXS(iFuelPin, iz)%XSnf = PinXS(iFuelPin, iz)%XSnf * PinVol(iFuelPin, iz)
    PinXS(iFuelPin, iz)%XSkf = PinXS(iFuelPin, iz)%XSkf * PinVol(iFuelPin, iz)
    DO i = 1, Pin(ixy)%nxy
      ipin = Pin(ixy)%pin(i); IF (ipin .EQ. iFuelPin) CYCLE
      PinXS(iFuelPin, iz)%XSnf = PinXS(iFuelPin, iz)%XSnf + PinXS(ipin, iz)%XSnf * PinVol(ipin, iz)
      PinXS(iFuelPin, iz)%XSkf = PinXS(iFuelPin, iz)%XSkf + PinXS(ipin, iz)%XSkf * PinVol(ipin, iz)
      PinXS(ipin, iz)%XSnf = 0.0; PinXS(ipin, iz)%XSkf = 0.0
    ENDDO
    PinXS(iFuelPin, iz)%XSnf = PinXS(iFuelPin, iz)%XSnf / PinVol(iFuelPin, iz)
    PinXS(iFuelPin, iz)%XSkf = PinXS(iFuelPin, iz)%XSkf / PinVol(iFuelPin, iz)
  ENDDO
ENDDO

END SUBROUTINE
    
#endif
