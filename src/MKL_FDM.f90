#include <defines.h>
!--- CNJ Edit : 1D Axial FDM Modules with Intel MKL
#ifdef __INTEL_MKL

MODULE MKL_FDM

USE MKL_3D
IMPLICIT NONE

REAL, POINTER, PRIVATE :: phiShape(:, :, :), phi(:, :, :), psi(:, :), src(:, :), lkg(:, :, :)
REAL, POINTER, PRIVATE :: l(:, :, :), d(:, :, :), u(:, :, :)
REAL, POINTER, PRIVATE :: S(:, :, :, :), F(:, :, :), Chi(:, :, :)
REAL, POINTER, PRIVATE :: Dtil(:, :, :, :), atil(:, :, :, :)

REAL, POINTER, PRIVATE :: hzFDM(:)
INTEGER, POINTER, PRIVATE :: cmRange(:, :), fmRange(:, :)
INTEGER, POINTER, PRIVATE :: cmMap(:), fmMap(:)
INTEGER, PRIVATE :: nzFDM, nDiv = 20

PRIVATE
PUBLIC :: AllocFDM, FDMDriver

CONTAINS

!--- Public Routines ------------------------------------------------------------------------------

SUBROUTINE AllocFDM

IMPLICIT NONE

INTEGER :: ng, nxy

ng = mklGeom%ng
nxy = mklGeom%nxy

CALL SetSubmesh()

ALLOCATE(phi(nzFDM, nxy, ng))
ALLOCATE(phiShape(nzFDM, nxy, ng)); phiShape = 1.0
ALLOCATE(psi(nzFDM, nxy))
ALLOCATE(src(nzFDM, nxy))
ALLOCATE(lkg(nzFDM, nxy, ng))

ALLOCATE(l(2 : nzFDM, nxy, ng))
ALLOCATE(d(nzFDM, nxy, ng))
ALLOCATE(u(1 : nzFDM - 1, nxy, ng))

ALLOCATE(S(nzFDM, nxy, ng, ng)); S = 0.0
ALLOCATE(F(nzFDM, nxy, ng)); F = 0.0
ALLOCATE(Chi(nzFDM, nxy, ng)); Chi = 0.0

ALLOCATE(Dtil(2, nzFDM, nxy, ng))
ALLOCATE(atil(2, nzFDM, nxy, ng))

END SUBROUTINE

SUBROUTINE FDMDriver(PinXS, eigv)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)
REAL :: eigv

INTEGER :: ig, ng

ng = mklGeom%ng

CALL SetSystem(PinXS)
CALL SetFlux()
CALL SetPsi()
CALL SetLeakage(PinXS)

DO ig = 1, ng
  CALL SetSource(ig, eigv)
  CALL SolveFDM(ig)
ENDDO

CALL GetShape()
CALL GetCurrent()

END SUBROUTINE

!--- Private Routines -----------------------------------------------------------------------------

SUBROUTINE SetSubmesh()

IMPLICIT NONE

INTEGER :: nzCMFD
INTEGER :: iz, izc, izf
INTEGER :: myzb, myze
REAL, POINTER :: hzfm(:)

nzCMFD = mklGeom%nzCMFD
myzb = mklGeom%myzb
myze = mklGeom%myze
hzfm => mklGeom%hzfm

nzFDM = nzCMFD * nDiv

ALLOCATE(hzFDM(0 : nzFDM + 1))
ALLOCATE(cmRange(myzb : myze, 2), fmRange(nzCMFD, 2))
ALLOCATE(cmMap(nzFDM), fmMap(nzFDM))

izf = 0
DO izc = myzb, myze
  cmRange(izc, 1) = izf + 1
  DO iz = mklGeom%fmRange(izc, 1), mklGeom%fmRange(izc, 2)
    fmRange(iz, 1) = izf + 1
    fmRange(iz, 2) = izf + nDiv
    fmMap(fmRange(iz, 1) : fmRange(iz, 2)) = iz
    hzFDM(izf + 1 : izf + nDiv) = hzfm(iz) / nDiv
    izf = izf + nDiv
  ENDDO
  cmRange(izc, 2) = izf
  cmMap(cmRange(izc, 1) : cmRange(izc, 2)) = izc
ENDDO
hzFDM(0) = hzfm(0) / nDiv
hzFDM(nzFDM + 1) = hzfm(nzCMFD + 1) / nDiv

END SUBROUTINE

SUBROUTINE SetSystem(PinXS)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)

INTEGER :: ig, igs, ipin, ipin_map, iz, izf, ierr
INTEGER :: ng, nxy
INTEGER :: myzb, myze
INTEGER, POINTER :: pinMap(:)
REAL :: mybeta, neighbeta
REAL, POINTER :: myD(:, :, :), neighD(:, :, :)

ng = mklGeom%ng
nxy = mklGeom%nxy
myzb = mklGeom%myzb
myze = mklGeom%myze
pinMap => mklGeom%pinMap

ALLOCATE(myD(nxy, ng, 2), neighD(nxy, ng, 2))

DO ig = 1, ng
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    myD(ipin, ig, bottom) = PinXS(ipin_map, myzb)%XSD(ig)
    myD(ipin, ig, top) = PinXS(ipin_map, myze)%XSD(ig)
  ENDDO
ENDDO

CALL InitFastComm()
CALL GetNeighborFast(ng * nxy, myD(:, :, bottom), neighD(:, :, top), bottom)
CALL GetNeighborFast(ng * nxy, myD(:, :, top), neighD(:, :, bottom), top)
CALL FinalizeFastComm()

IF (mklGeom%lBottom) THEN
  IF (mklGeom%AxBC(bottom) .EQ. VoidCell) neighD(:, :, bottom) = 0.5
  IF (mklGeom%AxBC(bottom) .EQ. RefCell) neighD(:, :, bottom) = 0.0
ENDIF
IF (mklGeom%lTop) THEN
  IF (mklGeom%AxBC(top) .EQ. VoidCell) neighD(:, :, top) = 0.5
  IF (mklGeom%AxBC(top) .EQ. RefCell) neighD(:, :, top) = 0.0
ENDIF

!$OMP PARALLEL PRIVATE(ipin_map, mybeta, neighbeta)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO iz = myzb, myze
    DO ipin = 1, nxy
      ipin_map = pinMap(ipin)
      DO izf = cmRange(iz, 1), cmRange(iz, 2)
        mybeta = PinXS(ipin_map, iz)%XSD(ig) / hzFDM(izf)
        d(izf, ipin, ig) = PinXS(ipin_map, iz)%XSr(ig) * hzFDM(izf)
        !--- Coupling with Bottom
        IF (izf .EQ. 1) THEN
          IF (iz .EQ. myzb .AND. mklGeom%lBottom) THEN
            neighbeta = neighD(ipin, ig, bottom) / 2.0
          ELSE
            neighbeta = neighD(ipin, ig, bottom) / hzFDM(izf - 1)
          ENDIF
        ELSE
          neighbeta = PinXS(ipin_map, cmMap(izf - 1))%XSD(ig) / hzFDM(izf - 1)
        ENDIF
        Dtil(bottom, izf, ipin, ig) = 2.0 * mybeta * neighbeta / (mybeta + neighbeta)
        atil(bottom, izf, ipin, ig) = mybeta / (mybeta + neighbeta)
        IF (izf .NE. 1) l(izf, ipin, ig) = - Dtil(bottom, izf, ipin, ig)
        d(izf, ipin, ig) = d(izf, ipin, ig) + Dtil(bottom, izf, ipin, ig)
        !--- Coupling with Top
        IF (izf .EQ. nzFDM) THEN
          IF (iz .EQ. myze .AND. mklGeom%lTop) THEN
            neighbeta = neighD(ipin, ig, top) / 2.0
          ELSE
            neighbeta = neighD(ipin, ig, top) / hzFDM(izf + 1)
          ENDIF
        ELSE
          neighbeta = PinXS(ipin_map, cmMap(izf + 1))%XSD(ig) / hzFDM(izf + 1)
        ENDIF
        Dtil(top, izf, ipin, ig) = 2.0 * mybeta * neighbeta / (mybeta + neighbeta)
        atil(top, izf, ipin, ig) = mybeta / (mybeta + neighbeta)
        IF (izf .NE. nzFDM) u(izf, ipin, ig) = - Dtil(top, izf, ipin, ig)
        d(izf, ipin, ig) = d(izf, ipin, ig) + Dtil(top, izf, ipin, ig)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ig = 1, ng
  DO ipin = 1, nxy
!    CALL ddttrf(nzFDM, l(:, ipin, ig), d(:, ipin, ig), u(:, ipin, ig), ierr)
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(ipin_map)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(4)
DO ig = 1, ng
  DO igs = 1, ng
    DO iz = myzb, myze
      DO ipin = 1, nxy
        ipin_map = pinMap(ipin)
        DO izf = cmRange(iz, 1), cmRange(iz, 2)
          IF (PinXS(ipin_map, iz)%XSs(ig)%ib .GT. igs) CYCLE
          IF (PinXS(ipin_map, iz)%XSs(ig)%ie .LT. igs) CYCLE
          S(izf, ipin, igs, ig) = PinXS(ipin_map, iz)%XSs(ig)%from(igs) * hzFDM(izf)
        ENDDO
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
      DO izf = cmRange(iz, 1), cmRange(iz, 2)
        F(izf, ipin, ig) = PinXS(ipin_map, iz)%XSnf(ig) * hzFDM(izf)
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
      DO izf = cmRange(iz, 1), cmRange(iz, 2)
        Chi(izf, ipin, ig) = PinXS(ipin_map, iz)%Chi(ig)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

DEALLOCATE(myD, neighD)

END SUBROUTINE

SUBROUTINE SetFlux()

IMPLICIT NONE

INTEGER :: ig, ipin, iz, izf
INTEGER :: ng, nxy, nzCMFD

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO ipin = 1, nxy
    DO iz = 1, nzCMFD
      DO izf = fmRange(iz, 1), fmRange(iz, 2)
        phi(izf, ipin, ig) = mklCMFD%phis(ipin, iz, ig)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

CALL vdmul(ng * nxy * nzFDM, phiShape, phi, phi)

END SUBROUTINE

SUBROUTINE SetPsi()

IMPLICIT NONE

INTEGER :: ig, igs, ipin, izf
INTEGER :: ng, nxy

ng = mklGeom%ng
nxy = mklGeom%nxy

psi = 0.0
DO ig = 1, ng
  !$OMP PARALLEL
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
  DO ipin = 1, nxy
    DO izf = 1, nzFDM
      psi(izf, ipin) = psi(izf, ipin) + F(izf, ipin, ig) * phi(izf, ipin, ig)
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDDO

END SUBROUTINE

SUBROUTINE SetLeakage(PinXS)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)

TYPE(superPin_Type), POINTER :: Pin(:)
INTEGER :: ig, ibd, ipin, ineighpin, ipin_map, iz, izc, izf
INTEGER :: ng, nxy, nzCMFD
INTEGER, POINTER :: pinMap(:), pinMapRev(:), planeMap(:)
REAL :: Dtil, Dhat, myphi, neighphi, radLkg

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
Pin => mklGeom%superPin
pinMap => mklGeom%pinMap
pinMapRev => mklGeom%pinMapRev
planeMap => mklGeom%planeMap

DO ig = 1, ng
  DO iz = 1, nzCMFD
    izc = planeMap(iz)
    !$OMP PARALLEL PRIVATE(ipin_map, ineighpin, myphi, neighphi, radLkg, Dtil, Dhat)
    !$OMP DO SCHEDULE(GUIDED)
    DO ipin = 1, nxy
      ipin_map = pinMap(ipin)
      myphi = mklCMFD%phis(ipin, iz, ig)
      radLkg = 0.0
      DO ibd = 1, 4
        ineighpin = Pin(ipin_map)%NeighIdx(ibd)
        ineighpin = pinMapRev(ineighpin)
        IF (ineighpin .EQ. VoidCell) THEN
          neighphi = 0.0
        ELSEIF (ineighpin .EQ. RefCell) THEN
          neighphi = myphi
        ELSE
          neighphi = mklCMFD%phis(ineighpin, iz, ig)
        ENDIF
        Dtil = PinXS(ipin_map, izc)%Dtil(ibd, ig)
        Dhat = PinXS(ipin_map, izc)%Dhat(ibd, ig)
        radLkg = radLkg - Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
      ENDDO
      radLkg = radLkg * mklGeom%hzfm(iz) / mklGeom%PinVolFm(ipin, iz)
      DO izf = fmRange(iz, 1), fmRange(iz, 2)
        lkg(izf, ipin, ig) = radLkg * hzFDM(izf)
      ENDDO
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE SetSource(ig, eigv)

IMPLICIT NONE

REAL :: eigv, reigv

INTEGER :: ig, igs, ipin, izf
INTEGER :: ng, nxy

ng = mklGeom%ng
nxy = mklGeom%nxy

reigv = 1.0 / eigv

src = 0.0

DO igs = 1, ng
  !$OMP PARALLEL
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
  DO ipin = 1, nxy
    DO izf = 1, nzFDM
      src(izf, ipin) = src(izf, ipin) + S(izf, ipin, igs, ig) * phi(izf, ipin, igs)
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDDO

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ipin = 1, nxy
  DO izf = 1, nzFDM
    src(izf, ipin) = src(izf, ipin) + reigv * Chi(izf, ipin, ig) * psi(izf, ipin)
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

CALL vdsub(nxy * nzFDM, src, lkg(:, :, ig), src)

END SUBROUTINE

SUBROUTINE SolveFDM(ig)

IMPLICIT NONE

INTEGER :: ig, ipin, ierr
INTEGER :: nxy

nxy = mklGeom%nxy

CALL dcopy(nxy * nzFDM, src, 1, phi(:, :, ig), 1)

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED)
DO ipin = 1, nxy
!  CALL ddttrsv('L', 'n', nzFDM, 1, l(:, ipin, ig), d(:, ipin, ig), u(:, ipin, ig), phi(:, ipin, ig), nzFDM, ierr)
!  CALL ddttrsv('U', 'n', nzFDM, 1, l(:, ipin, ig), d(:, ipin, ig), u(:, ipin, ig), phi(:, ipin, ig), nzFDM, ierr)
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE GetShape()

IMPLICIT NONE

INTEGER :: ig, ipin, iz, izf
INTEGER :: ng, nxy, nzCMFD
REAL :: phisum

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

!$OMP PARALLEL PRIVATE(phisum)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO iz = 1, nzCMFD
    DO ipin = 1, nxy
      phisum = 0.0
      DO izf = fmRange(iz, 1), fmRange(iz, 2)
        phisum = phisum + phi(izf, ipin, ig)
      ENDDO
      phisum = phisum / nDiv
      DO izf = fmRange(iz, 1), fmRange(iz, 2)
        phiShape(izf, ipin, ig) = phi(izf, ipin, ig) / phisum
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE GetCurrent()

IMPLICIT NONE

INTEGER :: ig, igs, ipin, iz, izf
INTEGER :: n, ng, nxy, nzCMFD
INTEGER :: myzb, myze
REAL, POINTER :: myphi(:, :, :), neighphi(:, :, :)
REAL :: myphis, neighphis, phisurf, jnet

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

ALLOCATE(myphi(nxy, ng, 2), neighphi(nxy, ng, 2))

DO ig = 1, ng
  CALL dcopy(nxy, phi(:, 1, ig), 1, myphi(:, ig, bottom), 1)
  CALL dcopy(nxy, phi(:, nzFDM, ig), 1, myphi(:, ig, top), 1)
ENDDO

CALL InitFastComm()
CALL GetNeighborFast(ng * nxy, myphi(:, :, bottom), neighphi(:, :, top), bottom)
CALL GetNeighborFast(ng * nxy, myphi(:, :, top), neighphi(:, :, bottom), top)
CALL FinalizeFastComm()

IF (mklGeom%lBottom) THEN
  IF (mklGeom%AxBC(bottom) .EQ. VoidCell) neighphi(:, :, bottom) = 0.0
  IF (mklGeom%AxBC(bottom) .EQ. RefCell) neighphi(:, :, bottom) = myphi(:, :, bottom)
ENDIF
IF (mklGeom%lTop) THEN
  IF (mklGeom%AxBC(top) .EQ. VoidCell) neighphi(:, :, top) = 0.0
  IF (mklGeom%AxBC(top) .EQ. RefCell) neighphi(:, :, top) = myphi(:, :, top)
ENDIF

!$OMP PARALLEL PRIVATE(izf, myphis, neighphis, phisurf, jnet)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO ipin = 1, nxy
    DO iz = 1, nzCMFD
      !--- Bottom Current
      izf = fmRange(iz, 1)
      myphis = phi(izf, ipin, ig)
      IF (izf .EQ. 1) THEN
        neighphis = neighphi(ipin, ig, bottom)
      ELSE
        neighphis = phi(izf - 1, ipin, ig)
      ENDIF
      phisurf = atil(bottom, izf, ipin, ig) * myphis + (1.0 - atil(bottom, izf, ipin, ig)) * neighphis
      jnet = - Dtil(bottom, izf, ipin, ig) * (neighphis - myphis)
      mklAxial%Jout(in, ig, bottom, iz, ipin) = 0.25 * phisurf - 0.5 * jnet
      mklAxial%Jout(out, ig, bottom, iz, ipin) = 0.25 * phisurf + 0.5 * jnet
      mklAxial%Jout(surf, ig, bottom, iz, ipin) = phisurf
      !--- Top Current
      izf = fmRange(iz, 2)
      myphis = phi(izf, ipin, ig)
      IF (izf .EQ. nzFDM) THEN
        neighphis = neighphi(ipin, ig, top)
      ELSE
        neighphis = phi(izf + 1, ipin, ig)
      ENDIF
      phisurf = atil(top, izf, ipin, ig) * myphis + (1.0 - atil(top, izf, ipin, ig)) * neighphis
      jnet = - Dtil(top, izf, ipin, ig) * (neighphis - myphis)
      mklAxial%Jout(in, ig, top, iz, ipin) = 0.25 * phisurf - 0.5 * jnet
      mklAxial%Jout(out, ig, top, iz, ipin) = 0.25 * phisurf + 0.5 * jnet
      mklAxial%Jout(surf, ig, top, iz, ipin) = phisurf
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

DEALLOCATE(myphi, neighphi)

END SUBROUTINE

END MODULE

#endif