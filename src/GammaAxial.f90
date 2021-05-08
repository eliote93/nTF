#include <defines.h>
!--- CNJ Edit : Driver Routine for Gamma Axial Calculation
#ifdef __INTEL_MKL
#ifdef __GAMMA_TRANSPORT

MODULE GAMMA_AXIAL

USE MKL_3D
IMPLICIT NONE

LOGICAL :: lFirstAxial = TRUE

CONTAINS

SUBROUTINE GammaAxialSolver
USE PARAM
USE PE_MOD,             ONLY : PE
USE IOUTIL,             ONLY : message
USE FILES,              ONLY : io8
USE TIMER,              ONLY : nTracer_dclock,      TimeChk
USE GAMMA_FLATMOC
IMPLICIT NONE

REAL :: AxNTimeBeg, AxNTimeEnd
INTEGER :: ierr

AxNTImeBeg = nTracer_dclock(.FALSE., .FALSE.)

CALL omp_set_num_threads(PE%nAxThread)

WRITE(mesg, '(a)') 'Performing Axial Calculation :  Flat Source MOC'
IF (PE%Master) CALL message(io8, TRUE, TRUE, mesg)

IF (lFirstAxial) CALL UpdateBoundaryFlux(mklGammaCMFD, mklGammaAxial)
CALL FlatMOCDriver(mklGammaCMFD, mklAxial, mklGammaAxial)
CALL SetAxialDhat(mklGammaCMFD, mklGammaAxial)

IF (lFirstAxial) lFirstAxial = FALSE

CALL MPI_BARRIER(PE%MPI_CMFD_COMM, ierr)

AxNTimeEnd = nTracer_dclock(.FALSE., .FALSE.)
TimeChk%AxialNodalTime = TimeChk%AxialNodalTime + (AxNTimeEnd - AxNTimeBeg)

END SUBROUTINE

SUBROUTINE UpdateBoundaryFlux(CMFD, Axial)

IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
TYPE(mklAxial_Type) :: Axial

INTEGER :: ng, nxy, nzCMFD
INTEGER :: ig, ipin
REAL :: atil, myphi, neighphi, surfphifdm

ng = CMFD%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

!$OMP PARALLEL PRIVATE(atil, myphi, neighphi, surfphifdm)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ig = 1, ng
  DO ipin = 1, nxy
    !--- Update Bottom Boundary Flux
    atil = Axial%atil(bottom, ipin, 1, ig)
    myphi = CMFD%phis(ipin, 1, ig)
    neighphi = CMFD%neighphis(ipin, ig, bottom)
    surfphifdm = atil * myphi + (1.0 - atil) * neighphi
    IF (mklGeom%lBottom .AND. mklGeom%AxBC(bottom) .EQ. VoidCell) THEN
      Axial%PhiAngIn(:, ig, ipin, bottom) = 0.0
    ELSE
      Axial%PhiAngIn(:, ig, ipin, bottom) = surfphifdm
    ENDIF
    !--- Update Top Boundary Flux
    atil = Axial%atil(top, ipin, nzCMFD, ig)
    myphi = CMFD%phis(ipin, nzCMFD, ig)
    neighphi = CMFD%neighphis(ipin, ig, top)
    surfphifdm = atil * myphi + (1.0 - atil) * neighphi
    IF (mklGeom%lTop .AND. mklGeom%AxBC(top) .EQ. VoidCell) THEN
      Axial%PhiAngIn(:, ig, ipin, top) = 0.0
    ELSE
      Axial%PhiAngIn(:, ig, ipin, top) = surfphifdm
    ENDIF
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetAxialDtil(CMFD, Axial)

IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
TYPE(mklAxial_Type) :: Axial

INTEGER :: ng, nxy, nz, nzCMFD
INTEGER :: ig, iz, ipin
INTEGER :: myzb, myze
REAL :: Dtil, atil, mybeta, neighbeta
REAL, POINTER :: hzfm(:)
REAL, POINTER :: myD(:, :, :), neighD(:, :, :)

ng = CMFD%ng
nxy = mklGeom%nxy
nz = mklGeom%nz
nzCMFD = mklGeom%nzCMFD
myzb = mklGeom%myzb
myze = mklGeom%myze
hzfm => mklGeom%hzfm

ALLOCATE(myD(nxy, ng, nzCMFD), neighD(nxy, ng, 2))

CALL GetDiffusionCoeff(CMFD, myD, neighD)

!$OMP PARALLEL PRIVATE(Dtil, atil, mybeta, neighbeta)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO iz = 1, nzCMFD
    DO ipin = 1, nxy
      mybeta = myD(ipin, ig, iz) / hzfm(iz)
      !--- Coupling with Bottom Plane
      IF (iz .EQ. 1) THEN
        IF (mklGeom%lBottom) THEN
          IF (mklGeom%AxBC(bottom) .EQ. RefCell) neighbeta = 0.0
          IF (mklGeom%AxBC(bottom) .EQ. VoidCell) neighbeta = 0.25
        ELSE
          neighbeta = neighD(ipin, ig, bottom) / hzfm(iz - 1)
        ENDIF
      ELSE
        neighbeta = myD(ipin, ig, iz - 1) / hzfm(iz - 1)
      ENDIF
      Dtil = 2.0 * mybeta * neighbeta / (mybeta + neighbeta)
      atil = mybeta / (mybeta + neighbeta)
      CMFD%AxDtil(bottom, ipin, iz, ig) = Dtil
      Axial%atil(bottom, ipin, iz, ig) = atil
      !--- Coupling with Top Plane
      IF (iz .EQ. nzCMFD) THEN
        IF (mklGeom%lTop) THEN
          IF (mklGeom%AxBC(top) .EQ. RefCell) neighbeta = 0.0
          IF (mklGeom%AxBC(top) .EQ. VoidCell) neighbeta = 0.25
        ELSE
          neighbeta = neighD(ipin, ig, top) / hzfm(iz + 1)
        ENDIF
      ELSE
        neighbeta = myD(ipin, ig, iz + 1) / hzfm(iz + 1)
      ENDIF
      Dtil = 2.0 * mybeta * neighbeta / (mybeta + neighbeta)
      atil = mybeta / (mybeta + neighbeta)
      CMFD%AxDtil(top, ipin, iz, ig) = Dtil
      Axial%atil(top, ipin, iz, ig) = atil
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

DEALLOCATE(neighD)

END SUBROUTINE

SUBROUTINE SetAxialDhat(CMFD, Axial)

IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
TYPE(mklAxial_Type) :: Axial

INTEGER :: ng, nxy, nz, nzCMFD
INTEGER :: ig, iz, ipin
REAL :: Dtil, Dhat, myphi, neighphi, jfdm, jmoc
REAL, POINTER :: neighphic(:, :, :)

ng = CMFD%ng
nxy = mklGeom%nxy
nz = mklGeom%nz
nzCMFD = mklGeom%nzCMFD

ALLOCATE(neighphic(ng, nxy, 2)); neighphic = 0.0

CALL InitFastComm()
CALL GetNeighborFast(ng * nxy, Axial%phic(:, :, 1), neighphic(:, :, top), bottom)
CALL GetNeighborFast(ng * nxy, Axial%phic(:, :, nzCMFD), neighphic(:, :, bottom), top)
CALL FinalizeFastComm()

!$OMP PARALLEL PRIVATE(Dtil, Dhat, myphi, neighphi, jfdm, jmoc)
DO iz = 1, nzCMFD
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
  DO ig = 1, ng
    DO ipin = 1, nxy
      IF (mklCntl%lRefPinFDM) THEN
        IF (mklGeom%lRefPin(ipin)) CYCLE
      ENDIF
      myphi = Axial%phic(ig, ipin, iz)
      !--- Coupling with Bottom Plane
      IF (iz .EQ. 1) THEN
        neighphi = neighphic(ig, ipin, bottom)
      ELSE
        neighphi = Axial%phic(ig, ipin, iz - 1)
      ENDIF
      Dtil = CMFD%AxDtil(bottom, ipin, iz, ig)
      jfdm = - Dtil * (neighphi - myphi)
      jmoc = Axial%Jout(out, ig, bottom, iz, ipin) - Axial%Jout(in, ig, bottom, iz, ipin)
      Dhat = - (jmoc - jfdm) / (myphi + neighphi)
      ! IF (myphi .LT. 0.0 .OR. neighphi .LT. 0.0) Dhat = 0.0
      CMFD%AxDhat(bottom, ipin, iz, ig) = Dhat
      !--- Coupling with Top Plane
      IF (iz .EQ. nzCMFD) THEN
        neighphi = neighphic(ig, ipin, top)
      ELSE
        neighphi = Axial%phic(ig, ipin, iz + 1)
      ENDIF
      Dtil = CMFD%AxDtil(top, ipin, iz, ig)
      jfdm = - Dtil * (neighphi - myphi)
      jmoc = Axial%Jout(out, ig, top, iz, ipin) - Axial%Jout(in, ig, top, iz, ipin)
      Dhat = - (jmoc - jfdm) / (myphi + neighphi)
      ! IF (myphi .LT. 0.0 .OR. neighphi .LT. 0.0) Dhat = 0.0
      CMFD%AxDhat(top, ipin, iz, ig) = Dhat
    ENDDO
  ENDDO
  !$OMP END DO
ENDDO
!$OMP END PARALLEL

DEALLOCATE(neighphic)

END SUBROUTINE

SUBROUTINE GetDiffusionCoeff(CMFD, myD, neighD)

IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
REAL, POINTER :: myD(:, :, :), neighD(:, :, :)

TYPE(GPinXS_Type), POINTER :: GPinXS(:, :)
INTEGER :: ig, ipin, ipin_map, iz, izf
INTEGER :: ng, nxy, nzCMFD
INTEGER, POINTER :: pinMap(:), planeMap(:)

GPinXS => CMFD%GPinXS
ng = CMFD%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
pinMap => mklGeom%pinMap
planeMap => CMFD%planeMap

!$OMP PARALLEL PRIVATE(ipin_map)
DO izf = 1, nzCMFD
  iz = planeMap(izf)
  !$OMP DO SCHEDULE(GUIDED)
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    DO ig = 1, ng
      myD(ipin, ig, izf) = GPinXS(ipin_map, iz)%XSD(ig)
    ENDDO
  ENDDO
  !$OMP END DO
ENDDO
!$OMP END PARALLEL

CALL InitFastComm()
CALL GetNeighborFast(ng * nxy, myD(:, :, 1), neighD(:, :, top), bottom)
CALL GetNeighborFast(ng * nxy, myD(:, :, nzCMFD), neighD(:, :, bottom), top)
CALL FinalizeFastComm()

END SUBROUTINE

END MODULE

#endif
#endif
