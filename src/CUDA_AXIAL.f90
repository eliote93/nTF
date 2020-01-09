#include <CUDADEFINES.h>

#ifdef __PGI

MODULE CUDA_AXIAL

USE CUDA_MASTER
USE CUDA_UTIL
IMPLICIT NONE

LOGICAL :: lFirstAxial = .TRUE.

CONTAINS

SUBROUTINE cuAxialSolver(cuCMFD, cuAxial, cuDevice, eigv)
USE PARAM
USE TIMER,          ONLY : nTracer_dclock,      TimeChk
USE PE_MOD,         ONLY : PE
USE IOUTIL,         ONLY : message
USE FILES,          ONLY : io8
USE CUDA_NODAL,     ONLY : cuNodalDriver
USE CUDA_AXMOC,     ONLY : cuFlatMOCDriver,     cuLinearMOCDriver
IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice
REAL :: eigv

TYPE(dim3) :: Blocks, Threads
INTEGER :: ng, nxy, nzCMFD
INTEGER :: ierr, iter, order
INTEGER(KIND = cuda_stream_kind) :: myStream
REAL :: AxNTimeBeg, AxNTimeEnd

AxNTimeBeg = nTracer_dclock(.FALSE., .FALSE.)

SELECT CASE (cuCntl%AxSolver)
CASE (NODAL)

  IF (cuCntl%lSENM) THEN

    WRITE(mesg, '(a)') 'Performing Axial Calculation :  P1 SENM'
    IF (PE%Master) CALL message(io8, .TRUE., .TRUE., mesg)

  ELSE

    WRITE(mesg, '(a)') 'Performing Axial Calculation :  P1 NEM'
    IF (PE%Master) CALL message(io8, .TRUE., .TRUE., mesg)

  ENDIF

  CALL cuNodalDriver(cuCMFD, cuAxial, cuDevice, eigv)

CASE (MOC)

  IF (lFirstAxial) CALL cuSetBoundaryFlux(cuCMFD, cuAxial, cuDevice)

  IF (cuCntl%lCASMO) THEN

    WRITE(mesg, '(a)') 'Performing Axial Calculation :  Linear Source MOC'
    IF (PE%Master) CALL message(io8, .TRUE., .TRUE., mesg)

    CALL cuLinearMOCDriver(cuCMFD, cuAxial, cuDevice, eigv)

  ELSE

    WRITE(mesg, '(a)') 'Performing Axial Calculation :  Flat Source MOC'
    IF (PE%Master) CALL message(io8, .TRUE., .TRUE., mesg)

    CALL cuFlatMOCDriver(cuCMFD, cuAxial, cuDevice, eigv)

  ENDIF

END SELECT

CALL cuSetAxialDhat(cuCMFD, cuAxial, cuDevice)

CALL MPI_BARRIER(MPI_CUDA_COMM, ierr)

lFirstAxial = .FALSE.

AxNTimeEnd = nTracer_dclock(.FALSE., .FALSE.)
TimeChk%AxialNodalTime = TimeChk%AxialNodalTime + (AxNTimeEnd - AxNTimeBeg)

END SUBROUTINE

SUBROUTINE cuSetAxialSourceOperator(CoreInfo, FmInfo, GroupInfo, cuCMFD, cuAxial, cuDevice)
USE TYPEDEF,        ONLY : CoreInfo_Type,       FmInfo_Type,        GroupInfo_Type
USE CUDA_AXMOC,     ONLY : cuSetSourceOperator
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice

CALL cuSetSourceOperator(CoreInfo, FmInfo, GroupInfo, cuCMFD, cuAxial, cuDevice)

END SUBROUTINE

SUBROUTINE cuSetBoundaryFlux(cuCMFD, cuAxial, cuDevice)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice

INTEGER :: n, ng, nxy, nPolar1D
INTEGER :: ig, ipin, ierr
INTEGER :: myzbf, myzef
REAL :: myphi
REAL(GPU_NODAL_PRECISION), ALLOCATABLE :: PhiAngIn(:, :, :, :)

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nPolar1D = cuGeometry%nPolar1D
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
n = ng * nxy * nPolar1D * 2

ALLOCATE(PhiAngIn(nPolar1D, nxy, ng, 2))

!$OMP PARALLEL PRIVATE(myphi)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ig = 1, ng
  DO ipin = 1, nxy
    !--- Update Bottom Boundary Flux
    myphi = cuCMFD%h_phis8(ig, ipin, myzbf)
    IF (cuDevice%lBottom .AND. cuGeometry%AxBC(bottom) .EQ. Void) THEN
      PhiAngIn(:, ipin, ig, bottom) = 0.0
    ELSE
      PhiAngIn(:, ipin, ig, bottom) = myphi
    ENDIF
    !--- Update Top Boundary Flux
    myphi = cuCMFD%h_phis8(ig, ipin, myzef)
    IF (cuDevice%lTop .AND. cuGeometry%AxBC(top) .EQ. Void) THEN
      PhiAngIn(:, ipin, ig, top) = 0.0
    ELSE
      PhiAngIn(:, ipin, ig, top) = myphi
    ENDIF
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

ierr = cudaMemcpy(cuAxial%PhiAngIn, PhiAngIn, n, cudaMemcpyHostToDevice)

DEALLOCATE(PhiAngIn)

END SUBROUTINE

SUBROUTINE cuSetAxialDtil(cuCMFD, cuAxial, cuDevice)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice

INTEGER :: ng, nxy
INTEGER :: ig, izf, ipin, ipin_map
INTEGER :: myzb, myze, myzbf, myzef
INTEGER, POINTER :: pinMap(:), planeMap(:)
REAL :: Dtil, mybeta, neighbeta
REAL, POINTER :: hzfm(:)
REAL, POINTER :: neighD(:, :, :)

ng = cuCMFD%ng
nxy = cuGeometry%nxyc
myzb = cuDevice%myzb
myze = cuDevice%myze
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
pinMap => cuGeometry%pinMap
planeMap => cuGeometry%planeMap
hzfm => cuGeometry%hzfm

ALLOCATE(neighD(ng, nxy, 2))

CALL cuGetNeighborD(cuDevice, cuCMFD%PinXS, neighD, myzb, myze, ng)

!$OMP PARALLEL PRIVATE(ipin_map, Dtil, mybeta, neighbeta)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO izf = myzbf, myzef
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    DO ig = 1, ng
      mybeta = cuCMFD%PinXS(ipin_map, planeMap(izf))%xsD(ig) / hzfm(izf)
      !--- Coupling with Bottom Plane
      IF (izf .EQ. myzbf) THEN
        IF (cuDevice%lBottom) THEN
          IF (cuGeometry%AxBC(bottom) .EQ. Reflective) neighbeta = 0.0
          IF (cuGeometry%AxBC(bottom) .EQ. Void) neighbeta = 0.25
        ELSE
          neighbeta = neighD(ig, ipin_map, bottom) / hzfm(izf - 1)
        ENDIF
      ELSE
        neighbeta = cuCMFD%PinXS(ipin_map, planeMap(izf - 1))%xsD(ig) / hzfm(izf - 1)
      ENDIF
      Dtil = 2.0 * mybeta * neighbeta / (mybeta + neighbeta)
      cuCMFD%AxDtil(bottom, ig, ipin, izf) = Dtil
      !--- Coupling with Top Plane
      IF (izf .EQ. myzef) THEN
        IF (cuDevice%lTop) THEN
          IF (cuGeometry%AxBC(top) .EQ. Reflective) neighbeta = 0.0
          IF (cuGeometry%AxBC(top) .EQ. Void) neighbeta = 0.25
        ELSE
          neighbeta = neighD(ig, ipin_map, top) / hzfm(izf + 1)
        ENDIF
      ELSE
        neighbeta = cuCMFD%PinXS(ipin_map, planeMap(izf + 1))%xsD(ig) / hzfm(izf + 1)
      ENDIF
      Dtil = 2.0 * mybeta * neighbeta / (mybeta + neighbeta)
      cuCMFD%AxDtil(top, ig, ipin, izf) = Dtil
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

DEALLOCATE(neighD)

END SUBROUTINE

SUBROUTINE cuSetAxialGcDtil(cuGcCMFD, cuDevice)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuGcCMFD
TYPE(cuDevice_Type) :: cuDevice

INTEGER :: ng, nxy
INTEGER :: ig, izf, ipin, ipin_map
INTEGER :: myzbf, myzef
INTEGER, POINTER :: pinMap(:)
REAL :: Dtil, mybeta, neighbeta
REAL, POINTER :: hzfm(:)
REAL, POINTER :: neighD(:, :, :)

ng = cuGcCMFD%ng
nxy = cuGeometry%nxyc
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
hzfm => cuGeometry%hzfm
pinMap => cuGeometry%pinMap

ALLOCATE(neighD(ng, nxy, 2))

CALL cuGetNeighborD(cuDevice, cuGcCMFD%PinXS, neighD, myzbf, myzef, ng)

!$OMP PARALLEL PRIVATE(ipin_map, Dtil, mybeta, neighbeta)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO izf = myzbf, myzef
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    DO ig = 1, ng
      mybeta = cuGcCMFD%PinXS(ipin_map, izf)%xsD(ig) / hzfm(izf)
      !--- Coupling with Bottom Plane
      IF (izf .EQ. myzbf) THEN
        IF (cuDevice%lBottom) THEN
          IF (cuGeometry%AxBC(bottom) .EQ. Reflective) neighbeta = 0.0
          IF (cuGeometry%AxBC(bottom) .EQ. Void) neighbeta = 0.25
        ELSE
          neighbeta = neighD(ig, ipin_map, bottom) / hzfm(izf - 1)
        ENDIF
      ELSE
        neighbeta = cuGcCMFD%PinXS(ipin_map, izf - 1)%xsD(ig) / hzfm(izf - 1)
      ENDIF
      Dtil = 2.0 * mybeta * neighbeta / (mybeta + neighbeta)
      cuGcCMFD%AxDtil(bottom, ig, ipin, izf) = Dtil
      !--- Coupling with Top Plane
      IF (izf .EQ. myzef) THEN
        IF (cuDevice%lTop) THEN
          IF (cuGeometry%AxBC(top) .EQ. Reflective) neighbeta = 0.0
          IF (cuGeometry%AxBC(top) .EQ. Void) neighbeta = 0.25
        ELSE
          neighbeta = neighD(ig, ipin_map, top) / hzfm(izf + 1)
        ENDIF
      ELSE
        neighbeta = cuGcCMFD%PinXS(ipin_map, izf + 1)%xsD(ig) / hzfm(izf + 1)
      ENDIF
      Dtil = 2.0 * mybeta * neighbeta / (mybeta + neighbeta)
      cuGcCMFD%AxDtil(top, ig, ipin, izf) = Dtil
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

DEALLOCATE(neighD)

END SUBROUTINE

SUBROUTINE cuSetAxialDhat(cuCMFD, cuAxial, cuDevice)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice

INTEGER :: n, ng, nxy, nzCMFD
INTEGER :: myzbf, myzef
INTEGER :: ig, izf, ipin, ierr
REAL :: Dtil, Dhat, myphi, neighphi, jfdm, jref
REAL, ALLOCATABLE :: neighphic(:, :, :), phic(:, :, :), Jout(:, :, :, :, :)

ng = cuCMFD%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
n = ng * nxy * nzCMFD

ALLOCATE(neighphic(ng, nxy, 2)); neighphic = 0.0
ALLOCATE(phic(ng, nxy, myzbf : myzef))
ALLOCATE(Jout(ng, nxy, myzbf : myzef, 2, 2))

ierr = cudaMemcpy(phic, cuAxial%phic, n, cudaMemcpyDeviceToHost)
ierr = cudaMemcpy(Jout, cuAxial%Jout, n * 4, cudaMemcpyDeviceToHost)

CALL InitMPIComm()
CALL MPIComm(ng * nxy, ng * nxy, phic(:, :, myzbf), neighphic(:, :, bottom), bottom, MPI_CUDA_COMM)
CALL MPIComm(ng * nxy, ng * nxy, phic(:, :, myzef), neighphic(:, :, top), top, MPI_CUDA_COMM)
CALL FinalizeMPIComm()

!$OMP PARALLEL PRIVATE(Dtil, Dhat, myphi, neighphi, jfdm, jref)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO izf = myzbf, myzef
  DO ipin = 1, nxy
    DO ig = 1, ng
      myphi = phic(ig, ipin, izf)
      !--- Coupling with Bottom Plane
      IF (izf .EQ. myzbf) THEN
        neighphi = neighphic(ig, ipin, bottom)
      ELSE
        neighphi = phic(ig, ipin, izf - 1)
      ENDIF
      Dtil = cuCMFD%AxDtil(bottom, ig, ipin, izf)
      jfdm = - Dtil * (neighphi - myphi)
      jref = Jout(ig, ipin, izf, bottom, out) - Jout(ig, ipin, izf, bottom, in)
      Dhat = - (jref - jfdm) / (myphi + neighphi)
      IF (myphi .LT. 0.0 .OR. neighphi .LT. 0.0) Dhat = 0.0
      cuCMFD%AxDhat(bottom, ig, ipin, izf) = Dhat
      !--- Coupling with Top Plane
      IF (izf .EQ. myzef) THEN
        neighphi = neighphic(ig, ipin, top)
      ELSE
        neighphi = phic(ig, ipin, izf + 1)
      ENDIF
      Dtil = cuCMFD%AxDtil(top, ig, ipin, izf)
      jfdm = - Dtil * (neighphi - myphi)
      jref = Jout(ig, ipin, izf, top, out) - Jout(ig, ipin, izf, top, in)
      Dhat = - (jref - jfdm) / (myphi + neighphi)
      IF (myphi .LT. 0.0 .OR. neighphi .LT. 0.0) Dhat = 0.0
      cuCMFD%AxDhat(top, ig, ipin, izf) = Dhat
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

DEALLOCATE(neighphic, phic, Jout)

END SUBROUTINE

SUBROUTINE cuSetAxialGcDhat(cuCMFD, cuGcCMFD, cuDevice)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD, cuGcCMFD
TYPE(cuDevice_Type) :: cuDevice

INTEGER :: ngc, nxy, nzCMFD
INTEGER :: myzbf, myzef
INTEGER :: ig, igc, izf, ipin
REAL :: Dtil, Dhat, myphi, neighphi, jfdm
REAL, ALLOCATABLE :: Jnet(:, :, :, :)

ngc = cuGcCMFD%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef

ALLOCATE(Jnet(2, ngc, nxy, myzbf : myzef)); Jnet = 0.0

!$OMP PARALLEL PRIVATE(Dtil, Dhat, myphi, neighphi, jfdm)

!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO izf = myzbf, myzef
  DO ipin = 1, nxy
    DO igc = 1, ngc
      DO ig = cuGeometry%GcStruct(1, igc), cuGeometry%GcStruct(2, igc)
        myphi = cuCMFD%h_phis8(ig, ipin, izf)
        !--- Condense Bottom Currents
        IF (izf .EQ. myzbf) THEN
          neighphi = cuCMFD%h_neighphis8(ig, ipin, bottom)
        ELSE
          neighphi = cuCMFD%h_phis8(ig, ipin, izf - 1)
        ENDIF
        Dtil = cuCMFD%AxDtil(bottom, ig, ipin, izf)
        Dhat = cuCMFD%AxDhat(bottom, ig, ipin, izf)
        jfdm = - Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
        Jnet(bottom, igc, ipin, izf) = Jnet(bottom, igc, ipin, izf) + jfdm
        !--- Condense Top Currents
        IF (izf .EQ. myzef) THEN
          neighphi = cuCMFD%h_neighphis8(ig, ipin, top)
        ELSE
          neighphi = cuCMFD%h_phis8(ig, ipin, izf + 1)
        ENDIF
        Dtil = cuCMFD%AxDtil(top, ig, ipin, izf)
        Dhat = cuCMFD%AxDhat(top, ig, ipin, izf)
        jfdm = - Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
        Jnet(top, igc, ipin, izf) = Jnet(top, igc, ipin, izf) + jfdm
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO

!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO izf = myzbf, myzef
  DO ipin = 1, nxy
    DO igc = 1, ngc
      myphi = cuGcCMFD%h_phis8(igc, ipin, izf)
      !--- Coupling with Bottom Plane
      IF (izf .EQ. myzbf) THEN
        neighphi = cuGcCMFD%h_neighphis8(igc, ipin, bottom)
      ELSE
        neighphi = cuGcCMFD%h_phis8(igc, ipin, izf - 1)
      ENDIF
      Dtil = cuGcCMFD%AxDtil(bottom, igc, ipin, izf)
      jfdm = - Dtil * (neighphi - myphi)
      Dhat = - (Jnet(bottom, igc, ipin, izf) - jfdm) / (myphi + neighphi)
      cuGcCMFD%AxDhat(bottom, igc, ipin, izf) = Dhat
      !--- Coupling with Top Plane
      IF (izf .EQ. myzef) THEN
        neighphi = cuGcCMFD%h_neighphis8(igc, ipin, top)
      ELSE
        neighphi = cuGcCMFD%h_phis8(igc, ipin, izf + 1)
      ENDIF
      Dtil = cuGcCMFD%AxDtil(top, igc, ipin, izf)
      jfdm = - Dtil * (neighphi - myphi)
      Dhat = - (Jnet(top, igc, ipin, izf) - jfdm) / (myphi + neighphi)
      cuGcCMFD%AxDhat(top, igc, ipin, izf) = Dhat
    ENDDO
  ENDDO
ENDDO
!$OMP END DO

!$OMP END PARALLEL

DEALLOCATE(Jnet)

END SUBROUTINE

SUBROUTINE cuGetNeighborD(cuDevice, PinXS, neighD, myzb, myze, ng)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(cuDevice_Type) :: cuDevice
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: myD(:, :, :), neighD(:, :, :)

INTEGER :: ipin
INTEGER :: n, ng, nxy
INTEGER :: myzb, myze

nxy = cuGeometry%nxyc
n = ng * nxy

ALLOCATE(myD(ng, nxy, 2))

DO ipin = 1, nxy
  myD(:, ipin, bottom) = PinXS(ipin, myzb)%xsD
ENDDO

DO ipin = 1, nxy
  myD(:, ipin, top) = PinXS(ipin, myze)%xsD
ENDDO

CALL cuGetNeighbor(n, myD, neighD, MPI_CUDA_COMM)

DEALLOCATE(myD)

END SUBROUTINE

END MODULE

#endif
