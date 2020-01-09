#include <defines.h>
!--- CNJ Edit : 1D Axial MOC Modules with Intel MKL
#ifdef __INTEL_MKL

MODULE MKL_AXIAL

USE CMFD_COMMON,    ONLY : odCMFD
USE MKL_3D
IMPLICIT NONE

#define mpCMFD

LOGICAL :: lFirstAxial = TRUE

CONTAINS

SUBROUTINE MKL_AxialSolver(CoreInfo, CmInfo, PinXS, eigv)
USE PARAM
USE TYPEDEF,            ONLY : CoreInfo_Type,       CmInfo_Type,        PinXS_Type
USE PE_MOD,             ONLY : PE
USE IOUTIL,             ONLY : message
USE FILES,              ONLY : io8
USE TIMER,              ONLY : nTracer_dclock,      TimeChk
USE MKL_LINMOC
USE MKL_FLATMOC
USE MKL_NODAL
USE MKL_FDM
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
REAL :: eigv

REAL :: AxNTimeBeg, AxNTimeEnd
INTEGER :: ierr

AxNTImeBeg = nTracer_dclock(.FALSE., .FALSE.)

CALL omp_set_num_threads(PE%nAxThread)

SELECT CASE (mklCntl%AxSolver)

CASE (NODAL)

!  IF (mklCntl%lSP3) THEN
!    IF (mklCntl%lSENM) THEN
!      WRITE(mesg, '(a)') 'Performing Axial Calculation :  SP3 SENM'
!      IF (PE%Master) CALL message(io8, TRUE, TRUE, mesg)
!    ELSE
!      WRITE(mesg, '(a)') 'Performing Axial Calculation :  SP3 NEM'
!      IF (PE%Master) CALL message(io8, TRUE, TRUE, mesg)
!    ENDIF
!  ELSE
!    IF (mklCntl%lSENM) THEN
!      WRITE(mesg, '(a)') 'Performing Axial Calculation :  P1 SENM'
!      IF (PE%Master) CALL message(io8, TRUE, TRUE, mesg)
!    ELSE
!      WRITE(mesg, '(a)') 'Performing Axial Calculation :  P1 NEM'
!      IF (PE%Master) CALL message(io8, TRUE, TRUE, mesg)
!    ENDIF
!  ENDIF
!  
!  CALL NodalDriver(PinXS, eigv)
!  CALL SetAxialDhat()
  
CASE (FDM)

!  IF (mklCntl%lSP3) THEN
!    WRITE(mesg, '(a)') 'Performing Axial Calculation :  SP3 FDM'
!    IF (PE%Master) CALL message(io8, TRUE, TRUE, mesg)
!  ELSE
!    WRITE(mesg, '(a)') 'Performing Axial Calculation :  P1 FDM'
!    IF (PE%Master) CALL message(io8, TRUE, TRUE, mesg)
!  ENDIF
!  
!  CALL FDMDriver(PinXS, eigv)
!  CALL SetAxialDhat()
  
CASE (MOC)
  
  IF (mklCntl%lCASMO) THEN
  
    WRITE(mesg, '(a)') 'Performing Axial Calculation :  Linear Source MOC'
    IF (PE%Master) CALL message(io8, TRUE, TRUE, mesg)
    
    IF (lFirstAxial) CALL UpdateBoundaryFlux()
    CALL LinearMOCDriver(PinXS, eigv)
    CALL SetAxialDhat()
    
  ELSE
  
    WRITE(mesg, '(a)') 'Performing Axial Calculation :  Flat Source MOC'
    IF (PE%Master) CALL message(io8, TRUE, TRUE, mesg)
    
    IF (lFirstAxial) CALL UpdateBoundaryFlux()
    CALL FlatMOCDriver(PinXS, eigv)
    CALL SetAxialDhat()
    
  ENDIF
  
END SELECT
      
IF (lFirstAxial) lFirstAxial = FALSE

CALL MPI_BARRIER(PE%MPI_CMFD_COMM, ierr)

AxNTimeEnd = nTracer_dclock(.FALSE., .FALSE.)
TimeChk%AxialNodalTime = TimeChk%AxialNodalTime + (AxNTimeEnd - AxNTimeBeg)

END SUBROUTINE

SUBROUTINE UpdateBoundaryFlux()

IMPLICIT NONE

INTEGER :: ng, nxy, nzCMFD
INTEGER :: ig, ipin
REAL :: atil, myphi, neighphi, surfphifdm

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
  
!$OMP PARALLEL PRIVATE(atil, myphi, neighphi, surfphifdm)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ig = 1, ng
  DO ipin = 1, nxy
    !--- Update Bottom Boundary Flux
    atil = mklAxial%atil(bottom, ipin, 1, ig)
    myphi = mklCMFD%phis(ipin, 1, ig)
    neighphi = mklCMFD%neighphis(ipin, ig, bottom)
    surfphifdm = atil * myphi + (1.0 - atil) * neighphi
    IF (mklGeom%lBottom .AND. mklGeom%AxBC(bottom) .EQ. VoidCell) THEN
      mklAxial%PhiAngIn(:, ig, ipin, bottom) = 0.0
    ELSE
      mklAxial%PhiAngIn(:, ig, ipin, bottom) = surfphifdm
    ENDIF
    !--- Update Top Boundary Flux
    atil = mklAxial%atil(top, ipin, nzCMFD, ig)
    myphi = mklCMFD%phis(ipin, nzCMFD, ig)
    neighphi = mklCMFD%neighphis(ipin, ig, top)
    surfphifdm = atil * myphi + (1.0 - atil) * neighphi
    IF (mklGeom%lTop .AND. mklGeom%AxBC(top) .EQ. VoidCell) THEN
      mklAxial%PhiAngIn(:, ig, ipin, top) = 0.0
    ELSE
      mklAxial%PhiAngIn(:, ig, ipin, top) = surfphifdm
    ENDIF
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetAxialDtil(PinXS)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)

INTEGER :: ng, nxy, nz, nzCMFD
INTEGER :: ig, iz, ipin
INTEGER :: myzb, myze
REAL :: Dtil, atil, mybeta, neighbeta
REAL, POINTER :: hzfm(:)
REAL, POINTER :: myD(:, :, :), neighD(:, :, :)

ng = mklGeom%ng
nxy = mklGeom%nxy
nz = mklGeom%nz
nzCMFD = mklGeom%nzCMFD
myzb = mklGeom%myzb
myze = mklGeom%myze
hzfm => mklGeom%hzfm

ALLOCATE(myD(nxy, ng, nzCMFD), neighD(nxy, ng, 2))

IF (mklCntl%odCMFD) CALL SetTheta()
CALL GetDiffusionCoeff(mklCMFD, myD, neighD)

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
      mklCMFD%AxDtil(bottom, ipin, iz, ig) = Dtil
      mklAxial%atil(bottom, ipin, iz, ig) = atil
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
      mklCMFD%AxDtil(top, ipin, iz, ig) = Dtil
      mklAxial%atil(top, ipin, iz, ig) = atil
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

DEALLOCATE(neighD)

END SUBROUTINE

SUBROUTINE SetAxialGcDtil(GcPinXS)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: GcPinXS(:, :)

INTEGER :: ngc, nxy, nz, nzCMFD
INTEGER :: ig, iz, ipin
REAL :: Dtil, mybeta, neighbeta
REAL, POINTER :: hzfm(:)
REAL, POINTER :: myD(:, :, :), neighD(:, :, :)

ngc = mklGeom%ngc
nxy = mklGeom%nxy
nz = mklGeom%nz
nzCMFD = mklGeom%nzCMFD
hzfm => mklGeom%hzfm

ALLOCATE(myD(nxy, ngc, nzCMFD), neighD(nxy, ngc, 2))

IF (mklCntl%odCMFD) CALL SetGcTheta()
CALL GetDiffusionCoeff(mklGcCMFD, myD, neighD)

!$OMP PARALLEL PRIVATE(Dtil, mybeta, neighbeta)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ngc
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
      mklGcCMFD%AxDtil(bottom, ipin, iz, ig) = Dtil
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
      mklGcCMFD%AxDtil(top, ipin, iz, ig) = Dtil
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

DEALLOCATE(myD, neighD)

END SUBROUTINE

SUBROUTINE SetTheta()
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)
INTEGER :: ig, ipin, ipin_map, iz, izf
INTEGER :: ng, nxy, nzCMFD
INTEGER, POINTER :: pinMap(:), planeMap(:)
REAL :: h
REAL, POINTER :: hzfm(:)

PinXS => mklCMFD%PinXS
ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
pinMap => mklGeom%pinMap
planeMap => mklGeom%planeMap
hzfm => mklGeom%hzfm

!$OMP PARALLEL PRIVATE(ipin_map)
DO izf = 1, nzCMFD
  h = hzfm(izf)
  iz = planeMap(izf)
  !$OMP DO SCHEDULE(GUIDED)
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    DO ig = 1, ng
      mklCMFD%theta(ig, ipin, izf) = odCMFD(h, PinXS(ipin_map, iz)%XSt(ig))
    ENDDO
  ENDDO
  !$OMP END DO
ENDDO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetGcTheta()
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)
INTEGER :: ig, igc, ipin, iz, izf
INTEGER :: ng, ngc, nxy, nzCMFD
REAL :: theta, phisum

ng = mklGeom%ng
ngc = mklGeom%ngc
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

!$OMP PARALLEL PRIVATE(theta, phisum)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO izf = 1, nzCMFD
  DO ipin = 1, nxy
    DO igc = 1, ngc
      theta = 0.0; phisum = 0.0
      DO ig = mklGeom%GcStruct(1, igc), mklGeom%GcStruct(2, igc)
        theta = theta + mklCMFD%theta(ig, ipin, izf) * mklCMFD%phis(ipin, izf, ig)
        phisum = phisum + mklCMFD%phis(ipin, izf, ig)
      ENDDO
      mklGcCMFD%theta(igc, ipin, izf) = theta / phisum
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetAxialDhat()

IMPLICIT NONE

INTEGER :: ng, nxy, nz, nzCMFD
INTEGER :: ig, iz, ipin
REAL :: Dtil, Dhat, myphi, neighphi, jfdm, jmoc
REAL, POINTER :: neighphic(:, :, :)

ng = mklGeom%ng
nxy = mklGeom%nxy
nz = mklGeom%nz
nzCMFD = mklGeom%nzCMFD

ALLOCATE(neighphic(ng, nxy, 2)); neighphic = 0.0

CALL InitFastComm()
CALL GetNeighborFast(ng * nxy, mklAxial%phic(:, :, 1), neighphic(:, :, top), bottom)
CALL GetNeighborFast(ng * nxy, mklAxial%phic(:, :, nzCMFD), neighphic(:, :, bottom), top)
CALL FinalizeFastComm()

!$OMP PARALLEL PRIVATE(Dtil, Dhat, myphi, neighphi, jfdm, jmoc)
DO iz = 1, nzCMFD
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
  DO ig = 1, ng
    DO ipin = 1, nxy
      IF (mklCntl%lRefPinFDM) THEN
        IF (mklGeom%lRefPin(ipin)) CYCLE
      ENDIF
      myphi = mklAxial%phic(ig, ipin, iz)
      !--- Coupling with Bottom Plane
      IF (iz .EQ. 1) THEN
        neighphi = neighphic(ig, ipin, bottom)
      ELSE
        neighphi = mklAxial%phic(ig, ipin, iz - 1)
      ENDIF
      Dtil = mklCMFD%AxDtil(bottom, ipin, iz, ig)
      jfdm = - Dtil * (neighphi - myphi)
      jmoc = mklAxial%Jout(out, ig, bottom, iz, ipin) - mklAxial%Jout(in, ig, bottom, iz, ipin)
      Dhat = - (jmoc - jfdm) / (myphi + neighphi)
      IF (myphi .LT. 0.0 .OR. neighphi .LT. 0.0) Dhat = 0.0
      mklCMFD%AxDhat(bottom, ipin, iz, ig) = Dhat
      !--- Coupling with Top Plane
      IF (iz .EQ. nzCMFD) THEN
        neighphi = neighphic(ig, ipin, top)
      ELSE
        neighphi = mklAxial%phic(ig, ipin, iz + 1)
      ENDIF
      Dtil = mklCMFD%AxDtil(top, ipin, iz, ig)
      jfdm = - Dtil * (neighphi - myphi)
      jmoc = mklAxial%Jout(out, ig, top, iz, ipin) - mklAxial%Jout(in, ig, top, iz, ipin)
      Dhat = - (jmoc - jfdm) / (myphi + neighphi)
      IF (myphi .LT. 0.0 .OR. neighphi .LT. 0.0) Dhat = 0.0
      mklCMFD%AxDhat(top, ipin, iz, ig) = Dhat
    ENDDO
  ENDDO
  !$OMP END DO
ENDDO
!$OMP END PARALLEL

DEALLOCATE(neighphic)

END SUBROUTINE

SUBROUTINE SetAxialGcDhat()

IMPLICIT NONE

INTEGER :: ng, ngc, nxy, nz, nzCMFD
INTEGER :: ig, igc, iz, ipin
REAL :: Dtil, Dhat, myphi, neighphi, jfdm
REAL, POINTER :: Jnet(:, :, :, :)

ng = mklGeom%ng
ngc = mklGeom%ngc
nxy = mklGeom%nxy
nz = mklGeom%nz
nzCMFD = mklGeom%nzCMFD

ALLOCATE(Jnet(2, nxy, nzCMFD, ngc)); Jnet = 0.0

!$OMP PARALLEL PRIVATE(Dtil, Dhat, myphi, neighphi, jfdm)

DO igc = 1, ngc
  DO ig = mklGeom%GcStruct(1, igc), mklGeom%GcStruct(2, igc)
    !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
    DO iz = 1, nzCMFD
      DO ipin = 1, nxy
        myphi = mklCMFD%phis(ipin, iz, ig)
        !--- Condense Bottom Currents
        IF (iz .EQ. 1) THEN
          neighphi = mklCMFD%neighphis(ipin, ig, bottom)
        ELSE
          neighphi = mklCMFD%phis(ipin, iz - 1, ig)
        ENDIF
        Dtil = mklCMFD%AxDtil(bottom, ipin, iz, ig)
        Dhat = mklCMFD%AxDhat(bottom, ipin, iz, ig)
        jfdm = - Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
        Jnet(bottom, ipin, iz, igc) = Jnet(bottom, ipin, iz, igc) + jfdm
        !--- Condense Top Currents
        IF (iz .EQ. nzCMFD) THEN
          neighphi = mklCMFD%neighphis(ipin, ig, top)
        ELSE
          neighphi = mklCMFD%phis(ipin, iz + 1, ig)
        ENDIF
        Dtil = mklCMFD%AxDtil(top, ipin, iz, ig)
        Dhat = mklCMFD%AxDhat(top, ipin, iz, ig)
        jfdm = - Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
        Jnet(top, ipin, iz, igc) = Jnet(top, ipin, iz, igc) + jfdm
      ENDDO
    ENDDO
    !$OMP END DO
  ENDDO
ENDDO

!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO igc = 1, ngc
  DO iz = 1, nzCMFD
    DO ipin = 1, nxy
      myphi = mklGcCMFD%phis(ipin, iz, igc)
      !--- Coupling with Bottom Plane
      IF (iz .EQ. 1) THEN
        neighphi = mklGcCMFD%neighphis(ipin, igc, bottom)
      ELSE
        neighphi = mklGcCMFD%phis(ipin, iz - 1, igc)
      ENDIF
      Dtil = mklGcCMFD%AxDtil(bottom, ipin, iz, igc)
      jfdm = - Dtil * (neighphi - myphi)
      Dhat = - (Jnet(bottom, ipin, iz, igc) - jfdm) / (myphi + neighphi)
      mklGcCMFD%AxDhat(bottom, ipin, iz, igc) = Dhat
      !--- Coupling with Top Plane
      IF (iz .EQ. nzCMFD) THEN
        neighphi = mklGcCMFD%neighphis(ipin, igc, top)
      ELSE
        neighphi = mklGcCMFD%phis(ipin, iz + 1, igc)
      ENDIF
      Dtil = mklGcCMFD%AxDtil(top, ipin, iz, igc)
      jfdm = - Dtil * (neighphi - myphi)
      Dhat = - (Jnet(top, ipin, iz, igc) - jfdm) / (myphi + neighphi)
      mklGcCMFD%AxDhat(top, ipin, iz, igc) = Dhat
    ENDDO
  ENDDO
ENDDO
!$OMP END DO

!$OMP END PARALLEL

DEALLOCATE(Jnet)
  
END SUBROUTINE

SUBROUTINE GetDiffusionCoeff(CMFD, myD, neighD)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
REAL, POINTER :: myD(:, :, :), neighD(:, :, :)

TYPE(PinXS_Type), POINTER :: PinXS(:, :)
INTEGER :: ig, ipin, ipin_map, iz, izf
INTEGER :: ng, nxy, nzCMFD
INTEGER, POINTER :: pinMap(:), planeMap(:)

PinXS => CMFD%PinXS
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
      myD(ipin, ig, izf) = PinXS(ipin_map, iz)%XSD(ig) + CMFD%theta(ig, ipin, izf)
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