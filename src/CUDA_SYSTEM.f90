#include <CUDADEFINES.h>

#ifdef __PGI

MODULE CUDA_SYSTEM

USE CUDA_MASTER
USE CUDA_UTIL
IMPLICIT NONE

CONTAINS

SUBROUTINE cuSetNaturalBiCGSystem(cuCMFD, cuDevice, l3dim, lShift, seigv)
USE CUDA_PCR
IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
LOGICAL :: l3dim, lShift
REAL :: seigv

TYPE(superPin_Type), POINTER :: Pin(:)
INTEGER, POINTER :: planeMap(:), pinMap(:), pinMapRev(:)
INTEGER, PARAMETER :: DOWN = 5, UP = 6, SELF = 7
INTEGER :: surf(7) = (/ UP, NORTH, WEST, SELF, EAST, SOUTH, DOWN /)
INTEGER :: ir, ic, iz, izf, ig, igs, ibd, isurf, ipin, ipin_map, ineighpin, ierr, tid
INTEGER :: dz, gb, ge
INTEGER :: ng, nxy, nzCMFD
INTEGER :: myzbf, myzef
REAL(GPU_CMFD_PRECISION), ALLOCATABLE :: offDiag(:, :, :)
REAL(8), ALLOCATABLE :: offDiag8(:, :, :)
REAL, POINTER :: PinVolFm(:, :), hzfm(:)
REAL :: diagVal(7), Dtil, Dhat, val

Pin => cuGeometry%superPin
ng = cuCMFD%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
PinVolFm => cuGeometry%PinVolFm
hzfm => cuGeometry%hzfm
planeMap => cuCMFD%planeMap
pinMap => cuGeometry%pinMap
pinMapRev => cuGeometry%pinMapRev

CALL createCsr(cuCMFD%M, (ng + 6) * ng * nxy * nzCMFD, ng * nxy * nzCMFD, ng * nxy * nzCMFD)
IF (cuCntl%lSPAI) CALL createCsr(cuCMFD%D, 7 * ng * nxy * nzCMFD, ng * nxy * nzCMFD, ng * nxy * nzCMFD)

DO izf = myzbf, myzef
  DO ipin = 1, nxy
    iz = planeMap(izf)
    ipin_map = pinMap(ipin)
    DO ig = 1, ng
      ir = ig + (ipin - 1) * ng + (izf - myzbf) * ng * nxy
	  diagVal = 0.0
      DO ibd = 1, 4
	    Dtil = cuCMFD%PinXS(ipin_map, iz)%Dtil(ibd, ig)
	    Dhat = cuCMFD%PinXS(ipin_map, iz)%Dhat(ibd, ig)
	    diagVal(ibd) = - (Dtil + Dhat) * hzfm(izf)
	    diagVal(SELF) = diagVal(SELF) + (Dtil - Dhat) * hzfm(izf)
	  ENDDO
	  IF (l3dim) THEN
	    DO ibd = 5, 6
	      Dtil = cuCMFD%AxDtil(ibd - 4, ig, ipin, izf)
          Dhat = cuCMFD%AxDhat(ibd - 4, ig, ipin, izf)
          diagVal(ibd) = - (Dtil + Dhat) * PinVolFm(ipin_map, izf) / hzfm(izf)
          diagVal(SELF) = diagVal(SELF) + (Dtil - Dhat) * PinVolFm(ipin_map, izf) / hzfm(izf)
	    ENDDO
      ENDIF
	  diagVal(SELF) = diagVal(SELF) + PinVolFm(ipin_map, izf) * cuCMFD%PinXS(ipin_map, iz)%XSr(ig)
      IF (lShift) THEN
        diagVal(SELF) = diagVal(SELF) - PinVolFm(ipin_map, izf)                                                     &
                                      * cuCMFD%PinXS(ipin_map, iz)%XSnf(ig)                                         &
                                      * cuCMFD%PinXS(ipin_map, iz)%Chi(ig) * seigv
      ENDIF
	  DO ibd = 1, 7
	    isurf = surf(ibd)
		SELECT CASE (isurf)
		CASE (SELF)
		  ineighpin = ipin
		  dz = 0
		CASE (UP)
		  ineighpin = ipin
		  dz = 1
		CASE (DOWN)
		  ineighpin = ipin
		  dz = -1
		CASE (NORTH, WEST, EAST, SOUTH)
		  ineighpin = Pin(ipin_map)%NeighIdx(isurf)
          ineighpin = pinMapRev(ineighpin)
          IF (ineighpin .LE. 0) diagVal(isurf) = 0.0
          IF (ineighpin .EQ. ipin) THEN
            Dtil = cuCMFD%PinXS(ipin_map, iz)%Dtil(isurf, ig)
            Dhat = cuCMFD%PinXS(ipin_map, iz)%Dhat(isurf, ig)
            diagVal(SELF) = diagVal(SELF) - (Dtil + Dhat) * hzfm(izf)
            diagVal(isurf) = 0.0
          ENDIF
          dz = 0
        END SELECT
        ic = ir + (ineighpin - ipin) * ng + dz * ng * nxy
        val = diagVal(isurf)
        CALL pushCsr(cuCMFD%M, val, ir, ic)
        IF (cuCntl%lSPAI) CALL pushCsr(cuCMFD%D, val, ir, ic)
      ENDDO
      gb = cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%ib
      ge = cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%ie
      IF (lShift) THEN
        DO igs = 1, ng
          IF (igs .EQ. ig) CYCLE
          ic = ir + (igs - ig)
          val = - cuCMFD%PinXS(ipin_map, iz)%XSnf(igs) * cuCMFD%PinXS(ipin_map, iz)%Chi(ig) * seigv
          IF (igs .GE. gb .AND. igs .LE. ge) THEN
            val = val - cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%from(igs)
          ENDIF
          val = val * PinVolFm(ipin_map, izf)
          CALL pushCsr(cuCMFD%M, val, ir, ic)
        ENDDO
      ELSE
        DO igs = gb, ge
          IF (igs .EQ. ig) CYCLE
          ic = ir + (igs - ig)
          val = - cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%from(igs) * PinVolFm(ipin_map, izf)
          CALL pushCsr(cuCMFD%M, val, ir, ic)
        ENDDO
      ENDIF
    ENDDO
  ENDDO
ENDDO

CALL finalizeSortCsr(cuCMFD%M, .TRUE.)
IF (cuCntl%lSPAI) CALL finalizeCsr(cuCMFD%D, .FALSE.)

IF (cuCMFD%lFirst) CALL cuSetPreconditioner(cuCMFD, cuDevice)

IF (l3dim) THEN

  ALLOCATE(offDiag(ng, nxy, 2), offDiag8(ng, nxy, 2))
  offDiag = 0.0; offDiag8 = 0.0

  IF (.NOT. cuDevice%lBottom) THEN
    !$OMP PARALLEL PRIVATE(Dtil, Dhat, ipin_map, val)
    !$OMP DO SCHEDULE(GUIDED)
    DO ipin = 1, nxy
      ipin_map = pinMap(ipin)
      DO ig = 1, ng
        Dtil = cuCMFD%AxDtil(bottom, ig, ipin, myzbf)
        Dhat = cuCMFD%AxDhat(bottom, ig, ipin, myzbf)
        val = - (Dtil + Dhat) * PinVolFm(ipin_map, myzbf) / hzfm(myzbf)
        offDiag(ig, ipin, bottom) = val
        offDiag8(ig, ipin, bottom) = val
      ENDDO
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL
  ENDIF
  IF (.NOT. cuDevice%lTop) THEN
    !$OMP PARALLEL PRIVATE(Dtil, Dhat, ipin_map, val)
    !$OMP DO SCHEDULE(GUIDED)
    DO ipin = 1, nxy
      ipin_map = pinMap(ipin)
      DO ig = 1, ng
        Dtil = cuCMFD%AxDtil(top, ig, ipin, myzef)
        Dhat = cuCMFD%AxDhat(top, ig, ipin, myzef)
        val = - (Dtil + Dhat) * PinVolFm(ipin_map, myzef) / hzfm(myzef)
        offDiag(ig, ipin, top) = val
        offDiag8(ig, ipin, top) = val
      ENDDO
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL
  ENDIF

  ierr = cudaMemcpy(cuCMFD%offDiag, offDiag, ng * nxy * 2, cudaMemcpyHostToDevice)
  ierr = cudaMemcpy(cuCMFD%offDiag8, offDiag8, ng * nxy * 2, cudaMemcpyHostToDevice)

  DEALLOCATE(offDiag, offDiag8)

ENDIF

cuCMFD%lFirst = .FALSE.

END SUBROUTINE

SUBROUTINE cuSetNaturalJacobiSystem(cuCMFD, cuDevice, l3dim)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
LOGICAL :: l3dim

TYPE(superPin_Type), POINTER :: Pin(:)
INTEGER, POINTER :: planeMap(:), pinMap(:), pinMapRev(:)
INTEGER, PARAMETER :: DOWN = 5, UP = 6, SELF = 7
INTEGER :: surf(7) = (/ UP, NORTH, WEST, SELF, EAST, SOUTH, DOWN /)
INTEGER :: ir, ic, iz, izf, ig, igs, ibd, isurf, ipin, ipin_map, ineighpin, ierr
INTEGER :: dz, gb, ge
INTEGER :: ng, nxy, nzCMFD
INTEGER :: myzbf, myzef
REAL(GPU_CMFD_PRECISION), ALLOCATABLE :: invDiag(:)
REAL(GPU_CMFD_PRECISION), ALLOCATABLE :: offDiag(:, :, :)
REAL(8), ALLOCATABLE :: offDiag8(:, :, :)
REAL, POINTER :: PinVolFm(:, :), hzfm(:)
REAL :: diagVal(7), Dtil, Dhat, val

Pin => cuGeometry%superPin
ng = cuCMFD%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
PinVolFm => cuGeometry%PinVolFm
hzfm => cuGeometry%hzfm
planeMap => cuCMFD%planeMap
pinMap => cuGeometry%pinMap
pinMapRev => cuGeometry%pinMapRev

CALL createCsr(cuCMFD%M, (ng + 6) * ng * nxy * nzCMFD, ng * nxy * nzCMFD, ng * nxy * nzCMFD)
CALL createCsr(cuCMFD%jM, (ng + 6) * ng * nxy * nzCMFD, ng * nxy * nzCMFD, ng * nxy * nzCMFD)

ALLOCATE(invDiag(ng * nxy * nzCMFD))

DO izf = myzbf, myzef
  iz = planeMap(izf)
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    DO ig = 1, ng
      !--- Diffusion
      ir = ig + (ipin - 1) * ng + (izf - myzbf) * ng * nxy
	  diagVal = 0.0
      DO ibd = 1, 4
	    Dtil = cuCMFD%PinXS(ipin_map, iz)%Dtil(ibd, ig)
	    Dhat = cuCMFD%PinXS(ipin_map, iz)%Dhat(ibd, ig)
	    diagVal(ibd) = - (Dtil + Dhat) * hzfm(izf)
	    diagVal(SELF) = diagVal(SELF) + (Dtil - Dhat) * hzfm(izf)
	  ENDDO
	  IF (l3dim) THEN
	    DO ibd = 5, 6
	      Dtil = cuCMFD%AxDtil(ibd - 4, ig, ipin, izf)
          Dhat = cuCMFD%AxDhat(ibd - 4, ig, ipin, izf)
          diagVal(ibd) = - (Dtil + Dhat) * PinVolFm(ipin_map, izf) / hzfm(izf)
          diagVal(SELF) = diagVal(SELF) + (Dtil - Dhat) * PinVolFm(ipin_map, izf) / hzfm(izf)
	    ENDDO
      ENDIF
	  diagVal(SELF) = diagVal(SELF) + PinVolFm(ipin_map, izf) * cuCMFD%PinXS(ipin_map, iz)%XSr(ig)
	  DO ibd = 1, 7
	    isurf = surf(ibd)
		SELECT CASE (isurf)
		CASE (SELF)
		  ineighpin = ipin
		  dz = 0
		CASE (UP)
		  ineighpin = ipin
		  dz = 1
		CASE (DOWN)
		  ineighpin = ipin
		  dz = -1
		CASE (NORTH, WEST, EAST, SOUTH)
		  ineighpin = Pin(ipin_map)%NeighIdx(isurf)
          ineighpin = pinMapRev(ineighpin)
          IF (ineighpin .LE. 0) diagVal(isurf) = 0.0
          dz = 0
        END SELECT
        ic = ir + (ineighpin - ipin) * ng + dz * ng * nxy
        val = diagVal(isurf)
        CALL pushCsr(cuCMFD%M, val, ir, ic)
        IF (isurf .EQ. SELF) THEN
          invDiag(ir) = 1.0 / val; CYCLE
        ENDIF
        CALL pushCsr(cuCMFD%jM, val, ir, ic)
      ENDDO
      !--- Scattering
      gb = cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%ib
      ge = cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%ie
      DO igs = gb, ge
        IF (igs .EQ. ig) CYCLE
        ic = ir + (igs - ig)
        val = - cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%from(igs) * PinVolFm(ipin_map, izf)
        CALL pushCsr(cuCMFD%M, val, ir, ic)
        CALL pushCsr(cuCMFD%jM, val, ir, ic)
      ENDDO
    ENDDO
  ENDDO
ENDDO

CALL finalizeSortCsr(cuCMFD%M, .TRUE.)
CALL finalizeSortCsr(cuCMFD%jM, .TRUE.)

ierr = cudaMemcpy(cuCMFD%invDiag, invDiag, ng * nxy * nzCMFD, cudaMemcpyHostToDevice)

DEALLOCATE(invDiag)

IF (l3dim) THEN

  ALLOCATE(offDiag(ng, nxy, 2), offDiag8(ng, nxy, 2))
  offDiag = 0.0; offDiag8 = 0.0

  IF (.NOT. cuDevice%lBottom) THEN
    !$OMP PARALLEL PRIVATE(Dtil, Dhat, ipin_map, val)
    !$OMP DO SCHEDULE(GUIDED)
    DO ipin = 1, nxy
      ipin_map = pinMap(ipin)
      DO ig = 1, ng
        Dtil = cuCMFD%AxDtil(bottom, ig, ipin, myzbf)
        Dhat = cuCMFD%AxDhat(bottom, ig, ipin, myzbf)
        val = - (Dtil + Dhat) * PinVolFm(ipin_map, myzbf) / hzfm(myzbf)
        offDiag(ig, ipin, bottom) = val
        offDiag8(ig, ipin, bottom) = val
      ENDDO
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL
  ENDIF
  IF (.NOT. cuDevice%lTop) THEN
    !$OMP PARALLEL PRIVATE(Dtil, Dhat, ipin_map, val)
    !$OMP DO SCHEDULE(GUIDED)
    DO ipin = 1, nxy
      ipin_map = pinMap(ipin)
      DO ig = 1, ng
        Dtil = cuCMFD%AxDtil(top, ig, ipin, myzef)
        Dhat = cuCMFD%AxDhat(top, ig, ipin, myzef)
        val = - (Dtil + Dhat) * PinVolFm(ipin_map, myzef) / hzfm(myzef)
        offDiag(ig, ipin, top) = val
        offDiag8(ig, ipin, top) = val
      ENDDO
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL
  ENDIF

  ierr = cudaMemcpy(cuCMFD%offDiag, offDiag, ng * nxy * 2, cudaMemcpyHostToDevice)
  ierr = cudaMemcpy(cuCMFD%offDiag8, offDiag8, ng * nxy * 2, cudaMemcpyHostToDevice)

  DEALLOCATE(offDiag, offDiag8)

ENDIF

cuCMFD%lFirst = .FALSE.

END SUBROUTINE

SUBROUTINE cuSetRedblackSORSystem(cuCMFD, cuDevice, l3dim)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
LOGICAL :: l3dim, lAxial

TYPE(superPin_Type), POINTER :: Pin(:)
INTEGER, POINTER :: pinMap(:), pinMapRev(:), planeMap(:), planeMapRB(:), pinMapRB(:), pinMapRevRB(:, :)
INTEGER, PARAMETER :: DOWN = 5, UP = 6, SELF = 7
INTEGER :: surf(6) = (/ UP, NORTH, WEST, EAST, SOUTH, DOWN /)
INTEGER :: ir, ic, irg, icg, iz, izf, ig, igs, ibd, isurf, ipin, ipin_map, ipin_global, ineighpin, ierr
INTEGER :: color, gb, ge
INTEGER :: nr, nc, nnz, ng, nxy, nxyRB(100, 2), nxyzRB(2), nzCMFD
INTEGER :: myzbf, myzef, rbBeg(2), rbEnd(2), rbRange(2, 100, 2)
INTEGER :: mp(2) = (/ 2, 1 /)
REAL(GPU_CMFD_PRECISION), ALLOCATABLE :: rbOffDiag(:, :, :)
REAL(8), ALLOCATABLE :: offDiag8(:, :, :)
REAL, POINTER :: PinVolFm(:, :), hzfm(:)
REAL :: diagVal(7), Dtil, Dhat, val

Pin => cuGeometry%superPin
ng = cuCMFD%ng
nxy = cuGeometry%nxyc
nxyRB = cuDevice%nxyRB
nxyzRB = cuDevice%nxyzRB
nzCMFD = cuDevice%nzCMFD
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
rbBeg = cuDevice%rbBeg
rbEnd = cuDevice%rbEnd
rbRange = cuDevice%rbRange
PinVolFm => cuGeometry%PinVolFm
hzfm => cuGeometry%hzfm
pinMap => cuGeometry%pinMap
pinMapRev => cuGeometry%pinMapRev
planeMap => cuCMFD%planeMap
planeMapRB => cuDevice%planeMapRB
pinMapRB => cuDevice%pinMapRB
pinMapRevRB => cuDevice%pinMapRevRB

CALL createCsr(cuCMFD%M, (ng + 6) * ng * nxy * nzCMFD, ng * nxy * nzCMFD, ng * nxy * nzCMFD)
CALL createCsr(cuCMFD%rbM(red), 6 * ng * nxyzRB(red), ng * nxyzRB(red), ng * nxyzRB(black))
CALL createCsr(cuCMFD%rbM(black), 6 * ng * nxyzRB(black), ng * nxyzRB(black), ng * nxyzRB(red))

IF (cuCMFD%lFirst) THEN
  CALL createCsr(cuCMFD%rbDiag(red), ng * ng * nxyzRB(red), ng * nxyzRB(red), ng * nxyzRB(red))
  CALL createCsr(cuCMFD%rbDiag(black), ng * ng * nxyzRB(black), ng * nxyzRB(black), ng * nxyzRB(black))
ELSE
  CALL clearCsr(cuCMFD%rbDiag(red))
  CALL clearCsr(cuCMFD%rbDiag(black))
ENDIF

! IF (cuCMFD%lFirst) THEN
!   CALL createBsr(cuCMFD%rbDiagB(red), nxyzRB(red), nxyzRB(red), nxyzRB(red), ng)
!   CALL createBsr(cuCMFD%rbDiagB(black), nxyzRB(black), nxyzRB(black), nxyzRB(black), ng)
! ELSE
!   CALL clearBsr(cuCMFD%rbDiagB(red))
!   CALL clearBsr(cuCMFD%rbDiagB(black))
! ENDIF
!
! DO color = red, black
!   ir = 1
!   DO ipin = rbBeg(color), rbEnd(color)
!     cuCMFD%rbDiagB(color)%bsrColIdx(ir) = ir
!     cuCMFD%rbDiagB(color)%bsrRowPtr(ir) = ir
!     ir = ir + 1
!   ENDDO
!   cuCMFD%rbDiagB(color)%bsrRowPtr(ir) = ir
! ENDDO

DO color = red, black
  DO ipin = rbBeg(color), rbEnd(color)
    ipin_map = pinMapRB(ipin)
    ipin_global = pinMapRev(ipin_map)
    izf = planeMapRB(ipin)
    iz = planeMap(izf)
    DO ig = 1, ng
      ir = ig + (ipin - rbBeg(color)) * ng
      irg = ig + (ipin - 1) * ng
      diagVal = 0.0
  	  DO ibd = 1, 4
  	    Dtil = cuCMFD%PinXS(ipin_map, iz)%Dtil(ibd, ig)
  	    Dhat = cuCMFD%PinXS(ipin_map, iz)%Dhat(ibd, ig)
  	    diagVal(ibd) = - (Dtil + Dhat) * hzfm(izf)
  	    diagVal(SELF) = diagVal(SELF) + (Dtil - Dhat) * hzfm(izf)
  	  ENDDO
  	  IF (l3dim) THEN
  	    DO ibd = 5, 6
  	      Dtil = cuCMFD%AxDtil(ibd - 4, ig, ipin_global, izf)
  	      Dhat = cuCMFD%AxDhat(ibd - 4, ig, ipin_global, izf)
          diagVal(ibd) = - (Dtil + Dhat) * PinVolFm(ipin_map, izf) / hzfm(izf)
          diagVal(SELF) = diagVal(SELF) + (Dtil - Dhat) * PinVolFm(ipin_map, izf) / hzfm(izf)
  	    ENDDO
      ENDIF
  	  diagVal(SELF) = diagVal(SELF) + PinVolFm(ipin_map, izf) * cuCMFD%PinXS(ipin_map, iz)%XSr(ig)
  	  DO ibd = 1, 6
  	    isurf = surf(ibd)
  		SELECT CASE (isurf)
  	  	CASE (UP)
          IF (izf .EQ. myzef) CYCLE
  		  ineighpin = pinMapRevRB(ipin_map, izf + 1)
  		CASE (DOWN)
          IF (izf .EQ. myzbf) CYCLE
  		  ineighpin = pinMapRevRB(ipin_map, izf - 1)
  		CASE (NORTH, WEST, EAST, SOUTH)
 	      ineighpin = Pin(ipin_map)%NeighIdx(isurf)
          ineighpin = pinMapRevRB(ineighpin, izf)
          IF (ineighpin .LE. 0) diagVal(isurf) = 0.0
          IF (ineighpin .EQ. ipin) THEN
            Dtil = cuCMFD%PinXS(ipin_map, iz)%Dtil(isurf, ig)
            Dhat = cuCMFD%PinXS(ipin_map, iz)%Dhat(isurf, ig)
            diagVal(SELF) = diagVal(SELF) - (Dtil + Dhat) * hzfm(izf)
            diagVal(isurf) = 0.0
          ENDIF
        END SELECT
        ic = ig + (ineighpin - rbBeg(mp(color))) * ng
        icg = ig + (ineighpin - 1) * ng
        val = diagVal(isurf)
        CALL pushCsr(cuCMFD%rbM(color), val, ir, ic)
        CALL pushCsr(cuCMFD%M, val, irg, icg)
      ENDDO
      gb = cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%ib; ge = ig - 1
      DO igs = gb, ge
        ic = ir + (igs - ig)
        icg = irg + (igs - ig)
        val = - cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%from(igs) * PinVolFm(ipin_map, izf)
        CALL pushCsr(cuCMFD%rbDiag(color), val, ir, ic)
!        cuCMFD%rbDiagB(color)%bsrVal(ig, igs, ipin - rbBeg(color) + 1) = val
        CALL pushCsr(cuCMFD%M, val, irg, icg)
      ENDDO
      ic = ir; icg = irg
      val = diagVal(SELF)
      CALL pushCsr(cuCMFD%rbDiag(color), val, ir, ic)
!      cuCMFD%rbDiagB(color)%bsrVal(ig, ig, ipin - rbBeg(color) + 1) = val
      CALL pushCsr(cuCMFD%M, val, irg, icg)
    ENDDO
  ENDDO
ENDDO

CALL finalizeSortCsr(cuCMFD%M, .TRUE.)
CALL finalizeCsr(cuCMFD%rbM(red), .TRUE.)
CALL finalizeCsr(cuCMFD%rbM(black), .TRUE.)
CALL finalizeCsr(cuCMFD%rbDiag(red), .TRUE.)
CALL finalizeCsr(cuCMFD%rbDiag(black), .TRUE.)

IF (cuCMFD%lFirst) THEN
  CALL SetBlockDiagSolveInfo(cuCMFD%rbDiag(red), cuDevice%mySparseHandle)
  CALL SetBlockDiagSolveInfo(cuCMFD%rbDiag(black), cuDevice%mySparseHandle)
ENDIF

! CALL copyDeviceBsr(cuCMFD%rbDiagB(red))
! CALL copyDeviceBsr(cuCMFD%rbDiagB(black))
!
! IF (cuCMFD%lFirst) THEN
!   CALL SetBlockDiagSolveInfoB(cuCMFD%rbDiagB(red), cuDevice%mySparseHandle)
!   CALL SetBlockDiagSolveInfoB(cuCMFD%rbDiagB(black), cuDevice%mySparseHandle)
! ENDIF

IF (l3dim) THEN

  ALLOCATE(rbOffDiag(ng, nxy * nzCMFD, 2), offDiag8(ng, nxy * nzCMFD, 2))
  rbOffDiag = 0.0; offDiag8 = 0.0

  IF (.NOT. cuDevice%lBottom) THEN
    !$OMP PARALLEL PRIVATE(Dtil, Dhat, ipin_map, ipin_global, val)
    DO color = red, black
      !$OMP DO SCHEDULE(GUIDED)
      DO ipin = rbRange(1, myzbf, color), rbRange(2, myzbf, color)
        ipin_map = pinMapRB(ipin)
        ipin_global = pinMapRev(ipin_map)
        DO ig = 1, ng
          Dtil = cuCMFD%AxDtil(bottom, ig, ipin_global, myzbf)
          Dhat = cuCMFD%AxDhat(bottom, ig, ipin_global, myzbf)
          val = - (Dtil + Dhat) * PinVolFm(ipin_map, myzbf) / hzfm(myzbf)
          rbOffDiag(ig, ipin, bottom) = val
          offDiag8(ig, ipin, bottom) = val
        ENDDO
      ENDDO
      !$OMP END DO
    ENDDO
    !$OMP END PARALLEL
  ENDIF
  IF (.NOT. cuDevice%lTop) THEN
    !$OMP PARALLEL PRIVATE(Dtil, Dhat, ipin_map, ipin_global, val)
    DO color = red, black
      !$OMP DO SCHEDULE(GUIDED)
      DO ipin = rbRange(1, myzef, color), rbRange(2, myzef, color)
        ipin_map = pinMapRB(ipin)
        ipin_global = pinMapRev(ipin_map)
        DO ig = 1, ng
          Dtil = cuCMFD%AxDtil(top, ig, ipin_global, myzef)
          Dhat = cuCMFD%AxDhat(top, ig, ipin_global, myzef)
          val = - (Dtil + Dhat) * PinVolFm(ipin_map, myzef) / hzfm(myzef)
          rbOffDiag(ig, ipin, top) = val
          offDiag8(ig, ipin, top) = val
        ENDDO
      ENDDO
      !$OMP END DO
    ENDDO
    !$OMP END PARALLEL
  ENDIF

  ierr = cudaMemcpy(cuCMFD%rOffDiag(:, 1 : nxyRB(myzbf, red)),                                                      &
                    rbOffDiag(:, rbRange(1, myzbf, red) : rbRange(2, myzbf, red), bottom),                          &
                    ng * nxyRB(myzbf, red), cudaMemcpyHostToDevice)
  ierr = cudaMemcpy(cuCMFD%rOffDiag(:, nxyRB(myzbf, red) + 1 : nxyRB(myzbf, red) + nxyRB(myzef, red)),              &
                    rbOffDiag(:, rbRange(1, myzef, red) : rbRange(2, myzef, red), top),                             &
                    ng * nxyRB(myzef, red), cudaMemcpyHostToDevice)
  ierr = cudaMemcpy(cuCMFD%bOffDiag(:, 1 : nxyRB(myzbf, black)),                                                    &
                    rbOffDiag(:, rbRange(1, myzbf, black) : rbRange(2, myzbf, black), bottom),                      &
                    ng * nxyRB(myzbf, black), cudaMemcpyHostToDevice)
  ierr = cudaMemcpy(cuCMFD%bOffDiag(:, nxyRB(myzbf, black) + 1 : nxyRB(myzbf, black) + nxyRB(myzef, black)),        &
                    rbOffDiag(:, rbRange(1, myzef, black) : rbRange(2, myzef, black), top),                         &
                    ng * nxyRB(myzef, black), cudaMemcpyHostToDevice)

  ierr = cudaMemcpy(cuCMFD%offDiag8(:, 1 : nxyRB(myzbf, red), bottom),                                              &
                    offDiag8(:, rbRange(1, myzbf, red) : rbRange(2, myzbf, red), bottom),                           &
                    ng * nxyRB(myzbf, red), cudaMemcpyHostToDevice)
  ierr = cudaMemcpy(cuCMFD%offDiag8(:, 1 : nxyRB(myzef, red), top),                                                 &
                    offDiag8(:, rbRange(1, myzef, red) : rbRange(2, myzef, red), top),                              &
                    ng * nxyRB(myzef, red), cudaMemcpyHostToDevice)
  ierr = cudaMemcpy(cuCMFD%offDiag8(:, nxyRB(myzbf, red) + 1 : nxy, bottom),                                        &
                    offDiag8(:, rbRange(1, myzbf, black) : rbRange(2, myzbf, black), bottom),                       &
                    ng * nxyRB(myzbf, black), cudaMemcpyHostToDevice)
  ierr = cudaMemcpy(cuCMFD%offDiag8(:, nxyRB(myzef, red) + 1 : nxy, top),                                           &
                    offDiag8(:, rbRange(1, myzef, black) : rbRange(2, myzef, black), top),                          &
                    ng * nxyRB(myzef, black), cudaMemcpyHostToDevice)

  DEALLOCATE(rbOffDiag, offDiag8)

ENDIF

CONTAINS

SUBROUTINE SetBlockDiagSolveInfo(Diag, sparseHandle)

IMPLICIT NONE

TYPE(CSR_CMFD_PRECISION) :: Diag
TYPE(cusparseHandle) :: sparseHandle

REAL(GPU_CMFD_PRECISION), POINTER, DEVICE :: csrVal(:)
INTEGER, POINTER, DEVICE :: csrRowPtr(:), csrColIdx(:)
INTEGER :: nr, nc, nnz
INTEGER :: ierr

csrVal => Diag%d_csrVal
csrRowPtr => Diag%d_csrRowPtr
csrColIdx => Diag%d_csrColIdx
nr = Diag%nr
nc = Diag%nc
nnz = Diag%nnz

ierr = cusparseCreateSolveAnalysisInfo(Diag%infoLU(1))
ierr = cusparsePcsrsv_analysis(sparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nnz,							    &
							   Diag%descrLU(1), csrVal, csrRowPtr, csrColIdx, Diag%infoLU(1))

END SUBROUTINE

! SUBROUTINE SetBlockDiagSolveInfoB(Diag, sparseHandle)
!
! IMPLICIT NONE
!
! TYPE(BSR_CMFD_PRECISION) :: Diag
! TYPE(cusparseHandle) :: sparseHandle
!
! REAL(GPU_CMFD_PRECISION), POINTER, DEVICE :: bsrVal(:, :, :)
! INTEGER, POINTER, DEVICE :: bsrRowPtr(:), bsrColIdx(:)
! INTEGER :: nbr, nbc, nBlock, blockSize, pBufferSize
! INTEGER :: ierr
!
! bsrVal => Diag%d_bsrVal
! bsrRowPtr => Diag%d_bsrRowPtr
! bsrColIdx => Diag%d_bsrColIdx
! nbr = Diag%nbr
! nbc = Diag%nbc
! nBlock = Diag%nBlock
! blockSize = Diag%blockSize
!
! ierr = cusparsePbsrsv2_bufferSize(sparseHandle, CUSPARSE_DIRECTION_ROW, CUSPARSE_OPERATION_NON_TRANSPOSE,           &
!                                   nbr, nBlock, Diag%descrLU(1), bsrVal, bsrRowPtr, bsrColIdx, blockSize,            &
!                                   Diag%infoLU(1), pBufferSize)
!
! print *, pBufferSize
! ALLOCATE(Diag%pBuffer(pBufferSize))
!
! ierr = cusparsePbsrsv2_analysis(sparseHandle, CUSPARSE_DIRECTION_ROW, CUSPARSE_OPERATION_NON_TRANSPOSE,             &
!                                 nbr, nBlock, Diag%descrLU(1), bsrVal, bsrRowPtr, bsrColIdx, blockSize,              &
!                                 Diag%infoLU(1), CUSPARSE_SOLVE_POLICY_NO_LEVEL, Diag%pBuffer)
!
! END SUBROUTINE

END SUBROUTINE

SUBROUTINE cuSetCMFDSourceOperator(cuCMFD, cuDevice, lUpscat)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
LOGICAL :: lUpscat

INTEGER, POINTER :: planeMap(:), planeMapRB(:), pinMap(:), pinMapRB(:)
INTEGER :: ir, ic, ig, igs, ipin, ipin_map, iz, izf
INTEGER :: color, gb, ge, rbBeg(2), rbEnd(2)
INTEGER :: ng, nxy, nzCMFD
INTEGER :: myzbf, myzef
REAL, POINTER :: PinVolFm(:, :)
REAL :: val

ng = cuCMFD%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
rbBeg = cuDevice%rbBeg
rbEnd = cuDevice%rbEnd
pinMap => cuGeometry%pinMap
pinMapRB => cuDevice%pinMapRB
planeMap => cuCMFD%planeMap
planeMapRB => cuDevice%planeMapRB
PinVolFm => cuGeometry%PinVolFm

IF (cuCntl%lNatural) THEN

  !--- Natural Ordered Operator

  IF (.NOT. cuDevice%lFuel) RETURN

  CALL createCsr(cuCMFD%F, ng * nxy * nzCMFD, nxy * nzCMFD, ng * nxy * nzCMFD)

  DO izf = myzbf, myzef
    DO ipin = 1, nxy
      ipin_map = pinMap(ipin)
      iz = planeMap(izf)
      ir = ipin + (izf - myzbf) * nxy
      DO ig = 1, ng
        ic = ig + (ipin - 1) * ng + (izf - myzbf) * ng * nxy
        val = cuCMFD%PinXS(ipin_map, iz)%XSnf(ig) * PinVolFm(ipin_map, izf)
        CALL pushCsr(cuCMFD%F, val, ir, ic)
      ENDDO
    ENDDO
  ENDDO

  CALL finalizeCsr(cuCMFD%F, .TRUE.)

  CALL createCsr(cuCMFD%Chi, ng * nxy * nzCMFD, ng * nxy * nzCMFD, nxy * nzCMFD)

  DO izf = myzbf, myzef
    DO ipin = 1, nxy
      ipin_map = pinMap(ipin)
      iz = planeMap(izf)
      ic = ipin + (izf - myzbf) * nxy
      DO ig = 1, ng
        ir = ig + (ipin - 1) * ng + (izf - myzbf) * ng * nxy
        val = cuCMFD%PinXS(ipin_map, iz)%Chi(ig)
        CALL pushCsr(cuCMFD%Chi, val, ir, ic)
      ENDDO
    ENDDO
  ENDDO

  CALL finalizeCsr(cuCMFD%Chi, .TRUE.)

ELSE

  !--- Red-Black Ordered Operator

  IF (lUpscat) THEN

    CALL createCsr(cuCMFD%S, ng * ng * nxy * nzCMFD, ng * nxy * nzCMFD, ng * nxy * nzCMFD)

    DO ipin = rbBeg(red), rbEnd(black)
      izf = planeMapRB(ipin)
      iz = planeMap(izf)
      ipin_map = pinMapRB(ipin)
      DO ig = 1, ng
        ir = ig + (ipin - 1) * ng
        gb = ig + 1; ge = cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%ie
        DO igs = gb, ge
          ic = ir + (igs - ig)
          val = cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%from(igs) * PinVolFm(ipin_map, izf)
          CALL pushCsr(cuCMFD%S, val, ir, ic)
        ENDDO
      ENDDO
    ENDDO

    CALL finalizeCsr(cuCMFD%S, .TRUE.)

  ENDIF

  IF (.NOT. cuDevice%lFuel) RETURN

  CALL createCsr(cuCMFD%F, ng * nxy * nzCMFD, nxy * nzCMFD, ng * nxy * nzCMFD)

  DO ipin = rbBeg(red), rbEnd(black)
    izf = planeMapRB(ipin)
    iz = planeMap(izf)
    ipin_map = pinMapRB(ipin)
    ir = ipin
    DO ig = 1, ng
      ic = ig + (ipin - 1) * ng
      val = cuCMFD%PinXS(ipin_map, iz)%XSnf(ig) * PinVolFm(ipin_map, izf)
      CALL pushCsr(cuCMFD%F, val, ir, ic)
    ENDDO
  ENDDO

  CALL finalizeCsr(cuCMFD%F, .TRUE.)

  CALL createCsr(cuCMFD%Chi, ng * nxy * nzCMFD, ng * nxy * nzCMFD, nxy * nzCMFD)

  DO ipin = rbBeg(red), rbEnd(black)
    izf = planeMapRB(ipin)
    iz = planeMap(izf)
    ipin_map = pinMapRB(ipin)
    ic = ipin
    DO ig = 1, ng
      ir = ig + (ipin - 1) * ng
      val = cuCMFD%PinXS(ipin_map, iz)%Chi(ig)
      CALL pushCsr(cuCMFD%Chi, val, ir, ic)
    ENDDO
  ENDDO

  CALL finalizeCsr(cuCMFD%Chi, .TRUE.)

ENDIF

END SUBROUTINE

END MODULE

#endif
