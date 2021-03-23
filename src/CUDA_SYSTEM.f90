#include <CUDADEFINES.h>

#ifdef __PGI

MODULE CUDA_SYSTEM

USE CUDA_MASTER
USE CUDA_UTIL
IMPLICIT NONE

CONTAINS

SUBROUTINE cuSetNaturalBiCGSparsity(cuCMFD, cuDevice)
USE CUDA_PCR
IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice

TYPE(superPin_Type), POINTER :: Pin(:)
INTEGER, POINTER :: planeMap(:), pinMap(:), pinMapRev(:)
INTEGER, PARAMETER :: DOWN = 5, UP = 6, SELF = 7
INTEGER :: surf(7) = (/ UP, NORTH, WEST, SELF, EAST, SOUTH, DOWN /)
INTEGER :: ir, ic, iz, izf, ig, igs, ibd, isurf, ipin, ipin_map, ineighpin
INTEGER :: dz, gb, ge
INTEGER :: ng, nxy, nzCMFD
INTEGER :: myzbf, myzef

IF (.NOT. cuCMFD%lFirst) RETURN

Pin => cuGeometry%superPin
ng = cuCMFD%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
planeMap => cuCMFD%planeMap
pinMap => cuGeometry%pinMap
pinMapRev => cuGeometry%pinMapRev

CALL createSparsity(cuCMFD%spM, (ng + 6) * ng * nxy * nzCMFD, ng * nxy * nzCMFD, ng * nxy * nzCMFD)
IF (cuCntl%lSPAI) CALL createSparsity(cuCMFD%spD, 7 * ng * nxy * nzCMFD, ng * nxy * nzCMFD, ng * nxy * nzCMFD)

DO izf = myzbf, myzef
  DO ipin = 1, nxy
    iz = planeMap(izf)
    ipin_map = pinMap(ipin)
    DO ig = 1, ng
      ir = ig + (ipin - 1) * ng + (izf - myzbf) * ng * nxy
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
          IF (ineighpin .LE. 0) CYCLE
          IF (ineighpin .EQ. ipin) CYCLE
          dz = 0
        END SELECT
        ic = ir + (ineighpin - ipin) * ng + dz * ng * nxy
        CALL setSparsity(cuCMFD%spM, ir, ic)
        IF (cuCntl%lSPAI) CALL setSparsity(cuCMFD%spD, ir, ic)
      ENDDO
      gb = cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%ib
      ge = cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%ie
      DO igs = gb, ge
        IF (igs .EQ. ig) CYCLE
        ic = ir + (igs - ig)
        CALL setSparsity(cuCMFD%spM, ir, ic)
      ENDDO
    ENDDO
  ENDDO
ENDDO

CALL finalizeSortSparsity(cuCMFD%spM)
IF (cuCntl%lSPAI) CALL finalizeSortSparsity(cuCMFD%spD)

END SUBROUTINE


SUBROUTINE cuSetNaturalBiCGSystem(cuCMFD, cuDevice, l3dim, lShift, seigv)
USE CUDA_PCR
USE OMP_LIB
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

LOGICAL, PARAMETER :: lPreSp = .TRUE. ! Preset sparsity pattern

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

IF (.NOT.lPreSp) THEN
  CALL createCsr(cuCMFD%M, (ng + 6) * ng * nxy * nzCMFD, ng * nxy * nzCMFD, ng * nxy * nzCMFD)
  IF (cuCntl%lSPAI) CALL createCsr(cuCMFD%D, 7 * ng * nxy * nzCMFD, ng * nxy * nzCMFD, ng * nxy * nzCMFD)
ELSE
  CALL createCsr(cuCMFD%M, cuCMFD%spM)
  IF (cuCntl%lSPAI) CALL createCsr(cuCMFD%D, cuCMFD%spD)
END IF

IF (.NOT. lPreSp) call omp_set_num_threads(1)
!$OMP PARALLEL PRIVATE(ir, ic, iz, isurf, ipin_map, ineighpin, dz, gb, ge, val, Dtil, Dhat, diagVal)
!$OMP DO COLLAPSE(3) SCHEDULE(GUIDED)
DO izf = myzbf, myzef
  DO ipin = 1, nxy
    DO ig = 1, ng
  iz = planeMap(izf)
    ipin_map = pinMap(ipin)
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
            diagVal(SELF) = diagVal(SELF) + diagVal(isurf) !-(Dtil + Dhat) * hzfm(izf)
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
!$OMP END DO
!$OMP END PARALLEL
IF (.NOT. lPreSp) call omp_set_num_threads(cuCntl%nCMFDHelper)

CALL finalizeSortCsr(cuCMFD%M, .TRUE.)
IF (cuCntl%lSPAI) CALL finalizeCsr(cuCMFD%D, .FALSE.)

!CALL cuSetPreconditioner(cuCMFD, cuDevice)
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
        CALL pushCsr(cuCMFD%M, val, irg, icg)
      ENDDO
      ic = ir; icg = irg
      val = diagVal(SELF)
      CALL pushCsr(cuCMFD%rbDiag(color), val, ir, ic)
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

SUBROUTINE cuSetNaturalTranBiCGSystem(cuCMFD, cuDevice, TranInfo, TranCntl, l3dim, lpcond, ldt)
USE TYPEDEF,          ONLY : TranInfo_Type,     TranCntl_Type
USE CUDA_PCR
IMPLICIT NONE
TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
LOGICAL :: l3dim, lpcond, ldt

TYPE(superPin_Type), POINTER :: Pin(:)
REAL(GPU_CMFD_PRECISION), ALLOCATABLE :: offDiag(:, :, :)
REAL(8), ALLOCATABLE :: offDiag8(:, :, :)
REAL, ALLOCATABLE, PINNED :: Expo_Alpha(:, :, :)
REAL, POINTER :: PinVolFm(:, :), hzfm(:)
INTEGER, POINTER :: planeMap(:), pinMap(:), pinMapRev(:)
REAL :: diagVal(7), Dtil, Dhat, val
REAL :: delt, theta, velo, rvdelt, rvalpha
REAL :: chieff, inveigv0
INTEGER, PARAMETER :: DOWN = 5, UP = 6, SELF = 7
INTEGER :: surf(7) = (/UP, NORTH, WEST, SELF, EAST, SOUTH, DOWN/)
INTEGER :: ng, nxy, nzCMFD, nprec
INTEGER :: myzbf, myzef, gb, ge
INTEGER :: izf, iz, ipin, ipin_map, ig, jg, iprec
INTEGER :: ir, ic, isurf, ineighpin, dz, ibd
INTEGER :: memSize, ierr

ng = cuCMFD%ng
nxy = cuGeometry%nxyc
nprec = cuGeometry%nprec
nzCMFD = cuDevice%nzCMFD
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef

delt = TranCntl%Delt(TranCntl%nowstep)
theta = TranCntl%Theta
inveigv0 = 1./TranInfo%eigv0

Pin => cuGeometry%superPin
PinVolFm => cuGeometry%PinVolFm
hzfm => cuGeometry%hzfm
planeMap => cuCMFD%planeMap
pinMap => cuGeometry%pinMap
pinMapRev => cuGeometry%pinMapRev

ALLOCATE(Expo_Alpha(ng, nxy, myzbf:myzef))
memSize = ng * nxy * nzCMFD
ierr = cudaMemcpyAsync(Expo_Alpha, cuTranCMInfo%Expo_Alpha, memSize, cudaMemcpyDeviceToHost, &
                       cuDevice%myStream)
ierr = cudaStreamSynchronize(cuDevice%myStream)

CALL createCsr(cuCMFD%M, (ng+6) * ng * nxy * nzCMFD, ng * nxy * nzCMFD, ng * nxy * nzCMFD)
DO izf = myzbf, myzef
  iz = planeMap(izf)
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    !IF(cuDevice%myzb .EQ. 1 .AND. ipin .EQ. 1) THEN
    !  DO ig = 1, ng
    !    PRINT*, ipin_map, ig
    !    PRINT'(4es14.6)', cuCMFD%PinXs(ipin_map, iz)%Dtil(:, ig)
    !    PRINT'(4es14.6)', cuCMFD%PinXs(ipin_map, iz)%Dhat(:, ig)
    !    PRINT'(2es14.6)', cuCMFD%AxDtil(:, ig, ipin, izf)
    !    PRINT'(2es14.6)', cuCMFD%AxDhat(:, ig, ipin, izf)
    !    PRINT'(4es14.6)', cuCMFD%PinXs(ipin_map, iz)%xsr(ig), cuCMFD%PinXS(ipin_map, iz)%velo(ig)
    !    PRINT'(4es14.6)', cuCMFD%PinXS(ipin_map, iz)%omega,  cuCMFD%PinXS(ipin_map, iz)%betat, cuCMFD%PinXS(ipin_map, iz)%Chi(ig), cuCMFD%PinXS(ipin_map, iz)%XSNF(ig)
    !  END DO
    !END IF
    !STOP
    DO ig = 1, ng
      ir = ig + (ipin-1) * ng + (izf - myzbf) * ng * nxy
      diagVal = 0.0
      DO ibd = 1,4
        Dtil = cuCMFD%PinXS(ipin_map, iz)%Dtil(ibd, ig)
        Dhat = cuCMFD%PinXS(ipin_map, iz)%Dhat(ibd, ig)
        diagVal(ibd) = - (Dtil + Dhat) * hzfm(izf)
        diagVal(SELF) = diagVal(SELF) + (Dtil - Dhat) * hzfm(izf)
      END DO
      IF (l3dim) THEN
        DO ibd = 5, 6
          Dtil = cuCMFD%AxDtil(ibd - 4, ig, ipin, izf)
          Dhat = cuCMFD%AxDhat(ibd - 4, ig, ipin, izf)
          diagVal(ibd) = - (Dtil + Dhat) * PinVolFm(ipin_map, izf) / hzfm(izf)
          diagVal(SELF) = diagVal(SELF) + (Dtil - Dhat) * PinVolFm(ipin_map, izf) / hzfm(izf)
        ENDDO
      ENDIF
      diagVal(SELF) = diagVal(SELF) + PinVolFm(ipin_map, izf) * cuCMFD%PinXS(ipin_map, iz)%XSr(ig)
      ! Transient Terms
      IF(ldt) THEN
        velo = cuCMFD%PinXS(ipin_map, iz)%velo(ig)
        rvdelt = 1./(velo * delt * theta)
        rvalpha = Expo_Alpha(ig, ipin, izf) / velo
        diagVal(SELF) = diagVal(SELF) + (rvdelt + rvalpha) * PinVolFm(ipin_map, izf)
      END IF

      IF(TranCntl%lchidk) THEN
        chieff = cuCMFD%PinXS(ipin_map, iz)%Chi(ig)
        DO iprec = 1, nprec
         chieff = chieff + (cuTranCMInfo%CellOmegap(iprec, ipin_map, iz) * TranInfo%lambda(iprec) - cuCMFD%PinXS(ipin_map,iz)%beta(iprec))* TranInfo%chidk(ig, iprec)
        END DO
      ELSE
        chieff = cuCMFD%PinXS(ipin_map, iz)%Chi(ig) + &
          TranInfo%chid(ig) * (cuCMFD%PinXS(ipin_map, iz)%omega - cuCMFD%PinXS(ipin_map, iz)%betat)
      END IF
       diagVal(SELF) = diagVal(SELF) &
                     - inveigv0 * chieff * cuCMFD%PinXS(ipin_map, iz)%XSNF(ig) * PinVolFm(ipin_map, izf)
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
      END DO
      gb = cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%ib
      ge = cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%ie
      DO jg = 1, ng
        IF(jg .EQ. ig) CYCLE
        ic = ir + (jg - ig)
        val = -inveigv0 * chieff * cuCMFD%PinXS(ipin_map, iz)%XSNF(jg)
        IF(jg .GE. gb .AND. jg .LE. ge) THEN
          val = val - cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%from(jg)
        END IF
        val = val * PinVolFm(ipin_map, izf)
        CALL pushCsr(cuCMFD%M, val, ir, ic)
      END DO
    END DO
  END DO
END DO

CALL finalizeSortCsr(cuCMFD%M, .TRUE.)

DEALLOCATE(Expo_Alpha)

IF (lpcond) CALL cuSetPreconditioner(cuCMFD, cuDevice)

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

END SUBROUTINE

SUBROUTINE cuSetTranCMFDSourceOperator(cuCMFD, cuDevice, TranInfo, TranCntl, lUpscat)
USE TYPEDEF,        ONLY : TranInfo_Type,       TranCntl_Type
IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
LOGICAL :: lUpscat

INTEGER, POINTER :: planeMap(:), planeMapRB(:), pinMap(:), pinMapRB(:)
INTEGER :: ir, ic, ig, igs, ipin, ipin_map, iz, izf, iprec
INTEGER :: color, gb, ge, rbBeg(2), rbEnd(2)
INTEGER :: ng, nxy, nzCMFD, nprec
INTEGER :: myzbf, myzef
REAL, POINTER :: PinVolFm(:, :)
REAL :: val

ng = cuCMFD%ng
nxy = cuGeometry%nxyc
nprec = cuGeometry%nprec
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

  CALL createCsr(cuCMFD%S, ng * ng * nxy * nzCMFD, ng * nxy * nzCMFD, ng * nxy * nzCMFD)
  DO izf = myzbf, myzef
    DO ipin = 1, nxy
      ipin_map = pinMap(ipin)
      iz = planeMap(izf)
      DO ig = 1, ng
        ir = ig + (ipin - 1) * ng + (izf - myzbf) * ng * nxy
        gb = cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%ib
        ge = cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%ie
        DO igs = gb, ge
          ic  = ir + (igs - ig)
          val = cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%From(igs) * PinVolFm(ipin_map, izf)
          CALL pushCsr(cuCMFD%S, val, ir, ic)
        END DO
      ENDDO
    ENDDO
  ENDDO

  CALL finalizeCsr(cuCMFD%S, .TRUE.)

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
        IF(TranCntl%lchidk) THEN
          val = cuCMFD%PinXS(ipin_map, iz)%Chi(ig) - TranInfo%Chid(ig) * cuCMFD%PinXS(ipin_map, iz)%betat
          DO iprec = 1, nprec
            val = val + cuTranCMInfo%CellOmegap(iprec, ipin_map, iz) * TranInfo%lambda(iprec) * TranInfo%chidk(ig, iprec)
          END DO
        ELSE
          val = cuCMFD%PinXS(ipin_map, iz)%Chi(ig) + &
            TranInfo%chid(ig) * (cuCMFD%PinXS(ipin_map, iz)%omega - cuCMFD%PinXS(ipin_map, iz)%betat)
        END IF
        CALL pushCsr(cuCMFD%Chi, val, ir, ic)
      ENDDO
    ENDDO
  ENDDO

  CALL finalizeCsr(cuCMFD%Chi, .TRUE.)

ELSE

  !!--- Red-Black Ordered Operator
  !
  !IF (lUpscat) THEN
  !
  !  CALL createCsr(cuCMFD%S, ng * ng * nxy * nzCMFD, ng * nxy * nzCMFD, ng * nxy * nzCMFD)
  !
  !  DO ipin = rbBeg(red), rbEnd(black)
  !    izf = planeMapRB(ipin)
  !    iz = planeMap(izf)
  !    ipin_map = pinMapRB(ipin)
  !    DO ig = 1, ng
  !      ir = ig + (ipin - 1) * ng
  !      gb = ig + 1; ge = cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%ie
  !      DO igs = gb, ge
  !        ic = ir + (igs - ig)
  !        val = cuCMFD%PinXS(ipin_map, iz)%XSs(ig)%from(igs) * PinVolFm(ipin_map, izf)
  !        CALL pushCsr(cuCMFD%S, val, ir, ic)
  !      ENDDO
  !    ENDDO
  !  ENDDO
  !
  !  CALL finalizeCsr(cuCMFD%S, .TRUE.)
  !
  !ENDIF
  !
  !IF (.NOT. cuDevice%lFuel) RETURN
  !
  !CALL createCsr(cuCMFD%F, ng * nxy * nzCMFD, nxy * nzCMFD, ng * nxy * nzCMFD)
  !
  !DO ipin = rbBeg(red), rbEnd(black)
  !  izf = planeMapRB(ipin)
  !  iz = planeMap(izf)
  !  ipin_map = pinMapRB(ipin)
  !  ir = ipin
  !  DO ig = 1, ng
  !    ic = ig + (ipin - 1) * ng
  !    val = cuCMFD%PinXS(ipin_map, iz)%XSnf(ig) * PinVolFm(ipin_map, izf)
  !    CALL pushCsr(cuCMFD%F, val, ir, ic)
  !  ENDDO
  !ENDDO
  !
  !CALL finalizeCsr(cuCMFD%F, .TRUE.)
  !
  !CALL createCsr(cuCMFD%Chi, ng * nxy * nzCMFD, ng * nxy * nzCMFD, nxy * nzCMFD)
  !
  !DO ipin = rbBeg(red), rbEnd(black)
  !  izf = planeMapRB(ipin)
  !  iz = planeMap(izf)
  !  ipin_map = pinMapRB(ipin)
  !  ic = ipin
  !  DO ig = 1, ng
  !    ir = ig + (ipin - 1) * ng
  !    val = cuCMFD%PinXS(ipin_map, iz)%Chi(ig) + &
  !          TranInfo%chid(ig) * (cuCMFD%PinXS(ipin_map, iz)%omega - cuCMFD%PinXS(ipin_map, iz)%betat)
  !    CALL pushCsr(cuCMFD%Chi, val, ir, ic)
  !  ENDDO
  !ENDDO
  !
  !CALL finalizeCsr(cuCMFD%Chi, .TRUE.)

ENDIF

END SUBROUTINE
SUBROUTINE PrintMatrix(lsseig, TranCntl, TranInfo)
USE TYPEDEF,      ONLY : PinXs_Type, TranCntl_TYPE, TranInfo_TYPE
USE CUDA_PCR
IMPLICIT NONE
LOGICAL :: lsseig
TYPE(TranCntl_TYPE), OPTIONAL :: TranCntl
TYPE(TranInfo_TYPE), OPTIONAL :: TranInfo

TYPE(PinXs_Type), POINTER :: PinXS(:,:)
INTEGER, POINTER :: pinMap(:), planeMap(:)
REAL, POINTER :: tmpb(:,:,:)
REAL, POINTER :: PinVolFm(:,:), hzfm(:)
TYPE(CSR_DOUBLE) :: M1g, P1g, SGM, FGM, ChiGM
CHARACTER(LEN = 20) :: filename
CHARACTER(LEN=2) :: grname
REAL :: diagVal(5), Dtil, Dhat, val
INTEGER :: ng, nxy, nzcmfd, gb, ge, nprec
INTEGER :: iz, izf, ipin_map, ixy, ibd, ig, ir, ic, jg, iprec
INTEGER :: ineighpin, isurf
INTEGER :: surfarray(5) = (/ 3, 2, 5, 4, 1 /)
INTEGER :: io
INTEGER, SAVE :: nedit = 0
INTEGER :: i, ierr

nedit = nedit + 1
ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzcmfd = cuDevice%nzCMFD
pinMap => cuGeometry%pinMap
planeMap => cuGeometry%planeMap
pinXs => cuCMFD%PinXS

hzfm => cuGeometry%hzfm
pinvolfm => cuGeometry%PinVolFm

io = 129

IF(lsseig) THEN
  write(filename, '(a5,i1)') 'NMSS_', nedit
ELSE
  write(filename, '(a5,i1)') 'NMTR_', nedit-1
END IF

CALL printCsr(cuCMFD%M, filename, io)
CALL cuPrepareSPAI(cuCMFD%M, cuCMFD%SPAI)

IF(lsseig) THEN
  write(filename, '(a6,i1)') 'NMSSP_', nedit
ELSE
  write(filename, '(a6,i1)') 'NMTRP_', nedit-1
END IF
CALL printCsr(cuCMFD%SPAI, filename, io)

IF(lsseig) THEN
  write(filename, '(a6,i1)') 'NMSSF_', nedit
  CALL printCsr(cuCMFD%F, filename, io)
  write(filename, '(a6,i1)') 'NMSSC_', nedit
  CALL printCsr(cuCMFD%Chi, filename, io)
ELSE
  ALLOCATE(tmpb(ng, nxy, nzcmfd))
  ierr = cudaMemcpy(tmpb, cuCMFD%src8, ng*nxy*nzcmfd, cudaMemcpyDeviceToHost)
  write(filename, '(a6,i1)') 'NMTRb_', nedit-1
  open(io, file = filename, status = 'replace')
  i = 0
  DO izf = 1, nzcmfd
    DO ixy = 1, nxy
      DO ig = 1, ng
        i = i + 1
        write(io,'(4i5es14.6)') i, ig, ixy, izf, tmpb(ig, ixy, izf)
      END DO
    END DO
  END DO
  close(io)
  write(filename, '(a6,i1)') 'GMTRb_', nedit-1
  open(io, file = filename, status = 'replace')
  i = 0
  DO ig = 1, ng
    DO izf = 1, nzcmfd
      DO ixy = 1, nxy
        i = i + 1
        write(io,'(4i5es14.6)') i, ig, ixy, izf, tmpb(ig, ixy, izf)
      END DO
    END DO
  END DO
  close(io)
  DEALLOCATE(tmpb)
END IF

DO ig = 1, ng
  CALL createCsr(M1g, 7*nxy*nzCMFD, nxy*nzcmfd, nxy*nzcmfd)
  DO izf = 1, nzCMFD
    iz = planeMap(izf)
    DO ixy = 1, nxy
      ir = ixy + (izf-1) * nxy
      ipin_map = pinMap(ixy)
      diagVal = 0.
      DO ibd = 1, 4
        Dtil = PinXS(ipin_map, iz)%Dtil(ibd, ig)
        Dhat = PinXS(ipin_map, iz)%Dhat(ibd, ig)
        diagVal(ibd) = - (Dtil + Dhat) * hzfm(izf)
        diagVal(5) = diagVal(5) + (Dtil - Dhat) * pinVolfm(ixy, izf) / hzfm(izf)
      END DO
      diagval(5) = diagval(5) + pinvolfm(ixy, izf) * pinxs(ipin_map, iz)%xsr(ig)
      IF(.NOT. lsseig) THEN
        diagVal(5) = diagVal(5) + Pinvolfm(ixy, izf) * (1./ (pinxs(ipin_map, iz)%velo(ig) * TranCntl%DelT(TranCntl%nowstep)))
      END IF
      DO ibd = 1,5
        isurf = surfarray(ibd)
        IF(isurf .NE. 5) THEN
          ineighpin = cuGeometry%superPin(ipin_map)%NeighIdx(isurf)
          ineighpin = cuGeometry%pinMapRev(ineighpin)
          IF(ineighpin .EQ. 0) diagVal(isurf) = 0.
          IF(ineighpin .EQ. ixy) THEN
            diagVal(isurf) = 0.
            Dtil = pinXS(ipin_map, iz)%DTIL(isurf,ig)
            Dhat = pinXS(ipin_map, iz)%DHAT(isurf,ig)
            diagVal(5) = diagVal(5) - (Dtil+Dhat) * hzfm(izf)
          END IF
        ELSE
          ineighpin = ixy
        END IF
        ic = ir + (ineighpin - ixy)
        CALL pushCsr(M1g, diagVal(isurf), ir, ic)
      END DO
    END DO
  END DO
  CALL finalizeSortCsr(M1g, .FALSE.)
  CALL cuPrepareSPAI(M1g, P1g)
  IF(ig .LT. 10) THEN
    write(grname, '(2i1)') 0, ig
  ELSE
    write(grname, '(i2)') ig
  END IF
  IF(lsseig) THEN
    write(filename, '(a5,a2,a1,i1)') 'GMSS_', grname, '_', nedit
  ELSE
    write(filename, '(a5,a2,a1,i1)') 'GMTR_', grname, '_', nedit-1
  END IF
  CALL printCsr(M1g, filename, io)
  IF(lsseig) THEN
    write(filename, '(a6,a2,a1,i1)') 'GMSSP_', grname, '_', nedit
  ELSE
    write(filename, '(a6,a2,a1,i1)') 'GMTRP_', grname, '_', nedit-1
  END IF
  CALL printCsr(P1g, filename, io)
END DO

CALL createCsr(FGM, ng*nxy*nzCMFD, nxy*nzCMFD, ng*nxy*nzCMFD)
DO ig = 1, ng
  CALL createCsr(SGM, ng*nxy*nzCMFD, nxy*nzCMFD, ng*nxy*nzCMFD)
  DO izf = 1, nzCMFD
    iz = planeMap(izf)
    DO ixy = 1, nxy
      ir = ixy + (izf - 1) * nxy  !+ (ig - 1) * nxy * nzCMFD
      ipin_map = pinMap(ixy)
      gb = cuCMFD%PinXS(ipin_map, iz)%Xss(ig)%ib
      ge = cuCMFD%PinXS(ipin_map, iz)%Xss(ig)%ie
      DO jg = gb, ge
        ic = ir + (jg - 1) * nxy * nzCMFD
        val = PinXS(ipin_map, iz)%Xss(ig)%From(jg) * PinVolFm(ipin_map, izf)
        CALL pushCsr(SGM, val, ir, ic)
      END DO
    END DO
  END DO
  CALL FinalizeCsr(SGM, .FALSE.)
  IF(ig .LT. 10) THEN
    write(grname, '(2i1)') 0, ig
  ELSE
    write(grname, '(i2)') ig
  END IF
  IF(lsseig) THEN
    write(filename, '(a6,a2,a1,i1)') 'GMSSS_', grname, '_', nedit
  ELSE
    write(filename, '(a6,a2,a1,i1)') 'GMTRS_', grname, '_', nedit-1
  END IF
  CALL printCsr(SGM, filename, io)
END DO

DO izf = 1, nzCMFD
  iz = planeMap(izf)
  DO ixy = 1, nxy
    ir = ixy + (izf - 1) * nxy
    ipin_map = pinMap(ixy)
    DO ig = 1, ng
      ic = ir + (ig - 1) * nxy * nzCMFD
      val = PinXS(ipin_map, iz)%XSNF(ig) * PinVolFm(ipin_map, izf)
      CALL pushCsr(FGM, val, ir, ic)
    END DO
  END DO
END DO
CALL FinalizeCsr(FGM, .FALSE.)
IF(lsseig) THEN
  write(filename, '(a6,a2,a1,i1)') 'GMSSF_', grname, '_', nedit
ELSE
  write(filename, '(a6,a2,a1,i1)') 'GMTRF_', grname, '_', nedit-1
END IF
CALL printCsr(FGM, filename, io)

DO ig = 1, ng
  CALL createCsr(ChiGM, nxy*nzCMFD, nxy*nzCMFD, nxy*nzCMFD)
  DO izf = 1, nzCMFD
    iz = planeMap(izf)
    DO ixy = 1, nxy
      ir = ixy + (izf - 1) * nxy !+ (ig - 1) * nxy * nzCMFD
      ic = ixy + (izf - 1) * nxy
      ipin_map = pinMap(ixy)
      IF(lsseig) THEN
        val = PinXS(ipin_map, iz)%Chi(ig) * PinVolFm(ipin_map, izf)
      ELSE
        IF(TranCntl%lchidk) THEN
          val = cuCMFD%PinXS(ipin_map, iz)%Chi(ig) - TranInfo%Chid(ig) * cuCMFD%PinXS(ipin_map, iz)%betat
          DO iprec = 1, nprec
            val = val + cuTranCMInfo%CellOmegap(iprec, ipin_map, iz) * TranInfo%lambda(iprec) * TranInfo%chidk(ig, iprec)
          END DO
        ELSE
          val = cuCMFD%PinXS(ipin_map, iz)%Chi(ig) + &
            TranInfo%chid(ig) * (cuCMFD%PinXS(ipin_map, iz)%omega - cuCMFD%PinXS(ipin_map, iz)%betat)
        END IF
      END IF
      CALL pushCsr(ChiGM, val, ir, ic)
    END DO
  END DO
  CALL FinalizeCsr(ChiGM, .FALSE.)
  IF(ig .LT. 10) THEN
    write(grname, '(2i1)') 0, ig
  ELSE
    write(grname, '(i2)') ig
  END IF
  IF(lsseig) THEN
    write(filename, '(a6,a2,a1,i1)') 'GMSSC_', grname, '_', nedit
  ELSE
    write(filename, '(a6,a2,a1,i1)') 'GMTRC_', grname, '_', nedit-1
  END IF
  CALL printCsr(ChiGM, filename, io)
END DO

CALL destroyCsr(M1g)
CALL destroyCsr(P1g)
CALL destroyCsr(SGM)
CALL destroyCsr(FGM)
CALL destroyCsr(ChiGM)


END SUBROUTINE

END MODULE

#endif
