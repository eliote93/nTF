#include <defines.h>
!--- CNJ Edit : 3D CMFD Acceleration Module with Intel MKL
#ifdef __INTEL_MKL

!--- Generalized Davidson Method Routines ---------------------------------------------------------

MODULE MKL_DAVIDSON

USE MKL_3D
IMPLICIT NONE

CONTAINS

SUBROUTINE ReorderFlux(dir)
USE PE_MOD,         ONLY : PE
IMPLICIT NONE

INTEGER :: dir

INTEGER :: n, ng, nxy, nzCMFD
INTEGER :: ig, ipin, iz, idx
REAL, POINTER :: phis(:, :, :), u(:)

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
n = ng * nxy * nzCMFD
phis => mklCMFD%phis
u => mklCMFD%Davidson%u

SELECT CASE (dir)

CASE (1)   !--- Map Group Major to Node Major

!$OMP PARALLEL PRIVATE(idx)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO iz = 1, nzCMFD
    DO ipin = 1, nxy
      idx = ig + (ipin - 1) * ng + (iz - 1) * ng * nxy
      u(idx) = phis(ipin, iz, ig)
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

CALL normalizeMPI(u, n, PE%MPI_CMFD_COMM)

CASE (2)   !--- Map Node Major to Group Major

CALL normalizeMPI(u, n, PE%MPI_CMFD_COMM)

!$OMP PARALLEL PRIVATE(idx)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO iz = 1, nzCMFD
  DO ipin = 1, nxy
    DO ig = 1, ng
      idx = ig + (ipin - 1) * ng + (iz - 1) * ng * nxy
      phis(ipin, iz, ig) = u(idx)
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SELECT

END SUBROUTINE

SUBROUTINE SetDavidsonOperator(PinXS, lDhat, l3dim)
USE PARAM
USE TYPEDEF,        ONLY : PinXS_Type
USE PE_MOD,         ONLY : PE
USE MKL_BILU,       ONLY : MKL_PrepareILU
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)
LOGICAL :: l3dim, lDhat

TYPE(mklDavidson_Type), POINTER :: Davidson
TYPE(superPin_Type), POINTER :: Pin(:)
INTEGER, POINTER :: pinMap(:), pinMapRev(:), planeMap(:)
INTEGER, PARAMETER :: DOWN = 5, UP = 6, SELF = 7
INTEGER :: surf(7) = (/ UP, NORTH, WEST, SELF, EAST, SOUTH, DOWN /)
INTEGER :: ir, ic, iz, izf, ig, igs, ibd, isurf, ipin, ipin_map, ineighpin, tid
INTEGER :: dz, gb, ge
INTEGER :: ng, nxy, nzCMFD, nbd
REAL, POINTER :: PinVolFm(:, :), hzfm(:)
REAL, POINTER :: send(:, :), recv(:, :)
REAL :: diagVal(7), Dtil, Dhat, val

Davidson => mklCMFD%Davidson
Pin => mklGeom%superPin
ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
hzfm => mklGeom%hzfm
planeMap => mklGeom%planeMap
pinMap => mklGeom%pinMap
pinMapRev => mklGeom%pinMapRev
PinVolFm => mklGeom%PinVolFm

!--- Set Node Major Fission Source Operator
CALL createCsr(Davidson%F, ng * ng * nxy * nzCMFD, ng * nxy * nzCMFD, ng * nxy * nzCMFD)
DO izf = 1, nzCMFD
  iz = planeMap(izf)
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    DO ig = 1, ng
      ir = ig + (ipin - 1) * ng + (izf - 1) * ng * nxy
      DO igs = 1, ng
        ic = ir + igs - ig
        val = PinXS(ipin_map, iz)%XSnf(igs) * PinXS(ipin_map, iz)%Chi(ig) * PinVolFm(ipin, izf)
        CALL pushCsr(Davidson%F, val, ir, ic)
      ENDDO
    ENDDO
  ENDDO
ENDDO
CALL finalizeCsr(Davidson%F, FALSE)

!--- Set Node Major Migration Operator
CALL createCsr(Davidson%M, (ng + 6) * ng * nxy * nzCMFD, ng * nxy * nzCMFD, ng * nxy * nzCMFD, PE%nCMFDThread)

!$OMP PARALLEL PRIVATE(tid, iz, ir, ic, isurf, ipin_map, ineighpin, Dtil, Dhat, dz, gb, ge, diagVal, val)
tid = omp_get_thread_num() + 1

!$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
DO izf = 1, nzCMFD
  DO ipin = 1, nxy
    DO ig = 1, ng
      ipin_map = pinMap(ipin)
      iz = planeMap(izf)
      ir = ig + (ipin - 1) * ng + (izf - 1) * ng * nxy
      diagVal = 0.0
      DO ibd = 1, 4
        Dtil = PinXS(ipin_map, iz)%Dtil(ibd, ig)
        Dhat = PinXS(ipin_map, iz)%Dhat(ibd, ig)
        diagVal(ibd) = - (Dtil + Dhat) * hzfm(izf)
        diagVal(SELF) = diagVal(SELF) + (Dtil - Dhat) * hzfm(izf)
      ENDDO
      IF (l3dim) THEN
        DO ibd = 5, 6
          Dtil = mklCMFD%AxDtil(ibd - 4, ipin, izf, ig)
          Dhat = mklCMFD%AxDhat(ibd - 4, ipin, izf, ig)
          diagVal(ibd) = - (Dtil + Dhat) * PinVolFm(ipin, izf) / hzfm(izf)
          diagVal(SELF) = diagVal(SELF) + (Dtil - Dhat) * PinVolFm(ipin, izf) / hzfm(izf)
        ENDDO
      ENDIF
      diagVal(SELF) = diagVal(SELF) + PinVolFm(ipin, izf) * PinXS(ipin_map, iz)%XSr(ig)
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
        CALL pushCsr(Davidson%M, val, ir, ic, tid)
      ENDDO
      gb = PinXS(ipin_map, iz)%XSs(ig)%ib
      ge = PinXS(ipin_map, iz)%XSs(ig)%ie
      DO igs = gb, ge
        IF (igs .EQ. ig) CYCLE
        ic = ir + (igs - ig)
        val = - PinXS(ipin_map, iz)%XSs(ig)%from(igs) * PinVolFm(ipin, izf)
        CALL pushCsr(Davidson%M, val, ir, ic, tid)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

CALL finalizeSortCsr(Davidson%M, FALSE)

!--- Set Preconditioner
CALL MKL_PrepareILU(Davidson%M, Davidson%ILU)

!--- Set Axial Off-diagonals for MPI
IF (l3dim) THEN

  CALL InitFastComm()

  IF (.NOT. mklGeom%lBottom) THEN
    !$OMP PARALLEL PRIVATE(Dtil, Dhat, val)
    !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
    DO ig = 1, ng
      DO ipin = 1, nxy
        Dtil = mklCMFD%AxDtil(bottom, ipin, 1, ig)
        Dhat = mklCMFD%AxDhat(bottom, ipin, 1, ig)
        val = - (Dtil + Dhat) * PinVolFm(ipin, 1) / hzfm(1)
        mklCMFD%AxOffDiag(ig, ipin, bottom) = val
      ENDDO
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL
  ENDIF
  !--- Transpose
  send => mklCMFD%AxOffDiag(:, :, bottom)
  recv => mklCMFD%trAxOffDiag(:, :, top)
  CALL GetNeighborFast(ng * nxy, send, recv, bottom)

  IF (.NOT. mklGeom%lTop) THEN
    !$OMP PARALLEL PRIVATE(Dtil, Dhat, val)
    !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
    DO ig = 1, ng
      DO ipin = 1, nxy
        Dtil = mklCMFD%AxDtil(top, ipin, nzCMFD, ig)
        Dhat = mklCMFD%AxDhat(top, ipin, nzCMFD, ig)
        val = - (Dtil + Dhat) * PinVolFm(ipin, nzCMFD) / hzfm(nzCMFD)
        mklCMFD%AxOffDiag(ig, ipin, top) = val
      ENDDO
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL
  ENDIF
  !--- Transpose
  send => mklCMFD%AxOffDiag(:, :, top)
  recv => mklCMFD%trAxOffDiag(:, :, bottom)
  CALL GetNeighborFast(ng * nxy, send, recv, top)

  CALL FinalizeFastComm()

ENDIF

END SUBROUTINE

! =================================================== !
!               Rayleigh-Ritz Procedure               !
!  http://en.wikipedia.org/wiki/Rayleigh-Ritz_method  !
! =================================================== !

SUBROUTINE RayleighRitz(Davidson, m, eigv)
USE PE_MOD,         ONLY : PE
IMPLICIT NONE

TYPE(mklDavidson_Type), POINTER :: Davidson
INTEGER :: n, m
REAL :: eigv

INTEGER :: i, ierr
INTEGER :: ng, nxy, nzCMFD
INTEGER :: bottomRange(2), topRange(2)
INTEGER, POINTER :: csrRowPtrM(:), csrRowPtrF(:), csrColIdxM(:), csrColIdxF(:)
REAL, POINTER :: csrValM(:), csrValF(:)
REAL, POINTER :: t(:), u(:), v(:), vM(:), vF(:), Mr(:, :), Fr(:, :), Vmat(:, :), Rmat(:, :)
REAL, POINTER :: s(:)
REAL :: proj, dot, buf_dot, norm

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
n = ng * nxy * nzCMFD
bottomRange = (/ 1, ng * nxy /)
topRange = (/ ng * nxy * (nzCMFD - 1) + 1, ng * nxy * nzCMFD /)

IF (m .GT. 1) THEN
  Mr => Davidson%Mr; ALLOCATE(Davidson%Mr(m, m))
  Fr => Davidson%Fr; ALLOCATE(Davidson%Fr(m, m))
  Vmat => Davidson%Vmat; ALLOCATE(Davidson%Vmat(n, m))
  Davidson%Mr(1 : m - 1, 1 : m - 1) = Mr; DEALLOCATE(Mr)
  Davidson%Fr(1 : m - 1, 1 : m - 1) = Fr; DEALLOCATE(Fr)
  CALL dcopy(n * (m - 1), Vmat, 1, Davidson%Vmat, 1); DEALLOCATE(Vmat)
ELSE
  IF (ASSOCIATED(Davidson%Mr)) DEALLOCATE(Davidson%Mr)
  IF (ASSOCIATED(Davidson%Fr)) DEALLOCATE(Davidson%Fr)
  IF (ASSOCIATED(Davidson%Vmat)) DEALLOCATE(Davidson%Vmat)
  ALLOCATE(Davidson%Mr(m, m))
  ALLOCATE(Davidson%Fr(m, m))
  ALLOCATE(Davidson%Vmat(n, m))
ENDIF

t => Davidson%t
u => Davidson%u
v => Davidson%Vmat(:, m)
csrRowPtrM => Davidson%M%csrRowPtr; csrRowPtrF => Davidson%F%csrRowPtr
csrColIdxM => Davidson%M%csrColIdx; csrColIdxF => Davidson%F%csrColIdx
csrValM => Davidson%M%csrVal; csrValF => Davidson%F%csrVal

IF (m .EQ. 1) CALL dcopy(n, u, 1, t, 1)

!--- Gram-Schmidt Orthogonalization

DO i = 1, m - 1
  proj = dotMPI(Davidson%Vmat(:, i), t, n, PE%MPI_CMFD_COMM)
  CALL daxpy(n, -proj, Davidson%Vmat(:, i), 1, t, 1)
ENDDO

CALL dcopy(n, t, 1, v, 1)
CALL normalizeMPI(v, n, PE%MPI_CMFD_COMM)

ALLOCATE(vM(n), vF(n))

!--- Construct Reduced Matrices

CALL MatMulMPIFast(Davidson%M, mklCMFD%AxOffDiag, v, vM, bottomRange, topRange)
CALL mkl_dcsrgemv('n', n, csrValF, csrRowPtrF, csrColIdxF, v, vF)

DO i = 1, m
  Davidson%Mr(i, m) = ddot(n, Davidson%Vmat(:, i), 1, vM, 1)
  Davidson%Fr(i, m) = ddot(n, Davidson%Vmat(:, i), 1, vF, 1)
ENDDO

CALL MatMulMPIFast(Davidson%M, mklCMFD%trAxOffDiag, v, vM, bottomRange, topRange, 't')
CALL mkl_dcsrgemv('t', n, csrValF, csrRowPtrF, csrColIdxF, v, vF)

DO i = 1, m
  Davidson%Mr(m, i) = ddot(n, Davidson%Vmat(:, i), 1, vM, 1)
  Davidson%Fr(m, i) = ddot(n, Davidson%Vmat(:, i), 1, vF, 1)
ENDDO

DEALLOCATE(vM, vF)

ALLOCATE(Mr(m, m), Fr(m, m), Rmat(m, m))
ALLOCATE(s(m))

CALL MPI_REDUCE(Davidson%Mr, Mr, m ** 2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, PE%MPI_CMFD_COMM, ierr)
CALL MPI_REDUCE(Davidson%Fr, Fr, m ** 2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, PE%MPI_CMFD_COMM, ierr)

IF (m .EQ. 1) THEN

  IF (PE%CMFDMaster) eigv = Fr(1, 1) / Mr(1, 1)

  CALL MPI_BCAST(eigv, 1, MPI_DOUBLE_PRECISION, 0, PE%MPI_CMFD_COMM, ierr)

ELSE

  IF (PE%CMFDMaster) THEN

    !--- Construct Reduced Eigenvalue Problem

    CALL MatInv(Mr, m)
    CALL dgemm('n', 'n', m, m, m, 1.0D0, Mr, m, Fr, m, 0.0D0, Rmat, m)

    !--- Find the Largest Ritz Pair

    CALL PowerIter(Rmat, s, eigv, m)

  ENDIF

  CALL MPI_BCAST(s, m, MPI_DOUBLE_PRECISION, 0, PE%MPI_CMFD_COMM, ierr)
  CALL MPI_BCAST(eigv, 1, MPI_DOUBLE_PRECISION, 0, PE%MPI_CMFD_COMM, ierr)

  CALL dgemv('n', n, m, 1.0D0, Davidson%Vmat, n, s, 1, 0.0D0, u, 1)

ENDIF

DEALLOCATE(Mr, Fr)
DEALLOCATE(s)

CONTAINS

SUBROUTINE MatInv(A, n)

IMPLICIT NONE

REAL, POINTER :: A(:, :)
INTEGER :: n

INTEGER :: ipiv(n), ierr, worksize
REAL, POINTER :: workspace(:)

CALL dgetrf(n, n, A, n, ipiv, ierr)

worksize = n * 64; ALLOCATE(workspace(worksize))
CALL dgetri(n, A, n, ipiv, workspace, worksize, ierr)
DEALLOCATE(workspace)

END SUBROUTINE

END SUBROUTINE

SUBROUTINE PowerIter(A, x, lambda, n)

IMPLICIT NONE

REAL, POINTER :: A(:, :), x(:), x_prev(:), e(:)
REAL :: lambda, norm, err, tol = 1.0D-8
INTEGER :: n

ALLOCATE(x_prev(n), e(n))

x = 1.0; err = 1.0
DO WHILE (err .GT. tol)
  x_prev = x
  CALL dgemv('n', n, n, 1.0D0, A, n, x_prev / lambda, 1, 0.0D0, x, 1)
  lambda = lambda * ddot(n, x, 1, x, 1) / ddot(n, x, 1, x_prev, 1)
  e = x - x_prev
  err = dnrm2(n, e, 1) / dnrm2(n, x, 1)
ENDDO

norm = dnrm2(n, x, 1)
x = x / norm

DEALLOCATE(x_prev, e)

END SUBROUTINE

SUBROUTINE SolveDavidsonEq(Davidson, eigv)
USE MKL_BILU,       ONLY : MKL_SolveILU
IMPLICIT NONE

TYPE(mklDavidson_Type), POINTER :: Davidson

INTEGER :: n
REAL, POINTER :: u(:), t(:), r(:)
REAL :: res, eigv

u => Davidson%u
t => Davidson%t
r => Davidson%r
n = mklGeom%ng * mklGeom%nxy * mklGeom%nzCMFD

CALL MKL_SolveILU(Davidson%ILU, t, r)

END SUBROUTINE

FUNCTION GetResidual(Davidson, eigv) RESULT(res)
USE PE_MOD,         ONLY : PE
IMPLICIT NONE

TYPE(mklDavidson_Type), POINTER :: Davidson
REAL :: eigv

INTEGER :: ng, nxy, nzCMFD, nr, nc, ierr
INTEGER :: bottomRange(2), topRange(2)
INTEGER, POINTER :: csrRowPtr(:), csrColIdx(:)
REAL, POINTER :: csrVal(:)
REAL, POINTER :: u(:), r(:), b(:)
REAL :: res

ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
bottomRange = (/ 1, ng * nxy /)
topRange = (/ ng * nxy * (nzCMFD - 1) + 1, ng * nxy * nzCMFD /)
nr = Davidson%M%nr
nc = Davidson%M%nc
u => Davidson%u
r => Davidson%r

ALLOCATE(b(nc))

csrRowPtr => Davidson%F%csrRowPtr
csrColIdx => Davidson%F%csrColIdx
csrVal => Davidson%F%csrVal
CALL mkl_dcsrgemv('n', nr, csrVal, csrRowPtr, csrColIdx, u, b)
CALL dscal(nc, 1.0D0 / eigv, b, 1)

CALL MatMulMPIFast(Davidson%M, mklCMFD%AxOffDiag, u, r, bottomRange, topRange)
CALL vdsub(nc, b, r, r)

DEALLOCATE(b)

res = normMPI(r, nc, PE%MPI_CMFD_COMM)

END FUNCTION

END MODULE

#endif
