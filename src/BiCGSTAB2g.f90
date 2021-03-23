#include <defines.h>
#define pre
#define FastComm

MODULE BICGSTAB2G_mod
USE PARAM
USE TYPEDEF, ONLY : CMFDLS_TYPE
USE itrcntl_mod, ONLY : InSolverItrCntl_TYPE
USE BasicOperation, ONLY : CP_CA, CP_VA, SUB_VA
USE OMP_LIB
USE BiLU_MOD,       ONLY : FreeBilu2gSolverEnv,     SetBilu2gSolverEnv,     SolveBilu2g
IMPLICIT NONE

PRIVATE

TYPE(CMFDLS_TYPE),POINTER, PRIVATE :: LS
REAL, POINTER, PRIVATE :: resv(:, :, :), resvhat(:, :, :)
REAL, POINTER, PRIVATE :: p(:, :, :), t(:, :, :), s(:, :, :), v(:, :, :), vy(:, :, :), vz(:, :, :)
REAL, POINTER, PRIVATE :: temp(:, :, :)
INTEGER, PRIVATE :: myzbf, myzef, nxy, nbd
#ifdef MPI_ENV
INTEGER, PRIVATE :: comm, myrank, nproc
#endif
INTEGER, PRIVATE :: nThread

PUBLIC :: AllocBiCGSTAB2G,  BiCGSTAB2G

CONTAINS

SUBROUTINE BiCGSTAB2G(A, Sol, RHS, itrcntl)
USE PARAM
#ifdef MPI_ENV
USE MPIComm_Mod, ONLY : GetNeighDat
#endif
USE TIMER,          only : nTracer_dclock, TimeChk

IMPLICIT NONE
TYPE(CMFDLS_TYPE), TARGET :: A
REAL, POINTER :: Sol(:, :, :), RHS(:, :, :)
TYPE(InSolverItrCntl_TYPE) :: itrcntl
INTEGER :: iter, itermax, itermin
INTEGER :: ixy, iz
INTEGER :: i, j, k, n
REAL :: rho,prho,alpha,omega,gamma,beta
REAL :: reserr, reserr0, temp
REAL :: Tbeg, Tend
LOGICAL :: lconv


TBeg = nTracer_dclock(FALSE, FALSE)
LS => A;
#ifdef pre
CALL SetBilu2gSolverEnv(A)
#endif
myzbf = LS%myzbf; myzef = LS%myzef; nxy = LS%nxy
n = nxy * (myzef - myzbf + 1)
nbd = LS%nbd
#ifdef MPI_ENV
comm = Ls%comm; myrank = Ls%myrank
nproc = Ls%nproc
#endif

itermax = itrcntl%ninmax
itermin = itrcntl%ninmin
!Initialized Solution related Array set
rho = 1._8; prho = 1._8; alpha = 1._8; omega = 1._8;

P(1:2, 1:nxy, myzbf-1:myzef+1) = 0.
v(1:2, 1:nxy, myzbf-1:myzef+1) = 0.

CALL Residual2g(Sol, RHS, ResV)

reserr0 = DotProduct2g(resv, resv)
reserr0 = sqrt(reserr0)

IF(reserr0 .lt. 1.0e-8_8) then
  reserr0 = 1._8
  itermax = 1
  RETURN
ENDIF

resvhat(1:2, 1:nxy, myzbf-1:myzef+1) = resv(1:2, 1:nxy, myzbf-1:myzef+1)
lconv = .FALSE.
DO iter = 1, 100
  prho = rho
  rho = DotProduct2G(resvhat, resv)
  beta = (rho*alpha)/(prho*omega)
  !P vector Update
  DO iz = myzbf, myzef
    DO ixy = 1, nxy
      P(:, ixy, iz) = resv(:, ixy, iz) + beta * (P(:, ixy, iz) - omega * v(:, ixy, iz))
    ENDDO
  ENDDO
#ifdef pre
  CALL SolveBilu2g(p, vy)
  CALL MatVecOP2g(vy, v)
#else
  CALL MatVecOP2G(p, v)
#endif
  gamma = DOtProduct2G(resvhat, v)
  alpha = rho/gamma
  !Update s = r(i-1)-alpha*v
  DO iz = myzbf, myzef
    DO ixy = 1, nxy
      s(:, ixy, iz) = resv(:, ixy, iz) - alpha * v(:, ixy, iz)
    ENDDO
  ENDDO
#ifdef pre
  CALL SolveBilu2g(s, vz)
  CALL MatVecOP2g(vz, t)
#else
  CALL MatVecOP2G(s, t)
#endif
  omega = DotProduct2G(t, s); temp = DotProduct2G(t, t)
  IF(temp .lt. 1.0E-20_8) THEN
    DO iz = myzbf, myzef
      DO ixy = 1, nxy
        sol(:, ixy, iz) = sol(:, ixy, iz) + alpha * p(:, ixy, iz)
      ENDDO
    ENDDO
    CALL Residual2G(Sol, RHS, ResV)
    reserr = DotProduct2G(resv, resv)
    reserr = sqrt(reserr)
    CONTINUE
    EXIT
  ELSE

    omega = omega/temp
 ENDIF

  !Update r = s - omega * t
  DO iz = myzbf, myzef
    DO ixy = 1, nxy
      resv(:, ixy, iz) = s(:, ixy, iz) - omega * t(:, ixy, iz)
    ENDDO
  ENDDO

  reserr = DotProduct2G(resv, resv)
  reserr = sqrt(reserr)
  !IF(myrank .eq. 0) PRINT *, alpha, omega, reserr
  !Solution Update
  DO iz = myzbf, myzef
    DO ixy = 1, nxy
#ifdef pre
     sol(:, ixy, iz) = sol(:, ixy, iz) + alpha * vy(:, ixy, iz) + omega * vz(:, ixy, iz)
#else
      sol(:, ixy, iz) = sol(:, ixy, iz) + alpha * p(:, ixy, iz) + omega * s(:, ixy, iz)
#endif
    ENDDO
  ENDDO
  !write(99, *) iter, reserr
  IF((reserr/reserr0) .lt. itrcntl%convcrit ) lconv = .TRUE.
  IF(iter .GT. itermin .AND. lCONV) EXIT
  !IF(reserr .lt. 1.0E-9_8) EXIT
  !IF(lCONV) EXIT
ENDDO

Tend = nTracer_dclock(FALSE, FALSE)
TimeChk%AxBTime = TimeChk%AxBTime + (Tend - Tbeg)

END SUBROUTINE

SUBROUTINE Residual2G(x, y, residualVec)
IMPLICIT NONE
REAL, POINTER :: x(:, :, :), y(:, :, :) , residualVec(:, :, :)
INTEGER :: ixy, iz
CALL MatVecOP2g(x, residualVec)
DO iz = myzbf, myzef
  DO ixy = 1, nxy
    ResidualVec(1:2, ixy, iz) = y(1:2, ixy, iz) - ResidualVec(1:2, ixy, iz)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE AllocBiCGSTAB2G(nxy0, myzbf0, myzef0)
USE ALLOCS
IMPLICIT NONE
INTEGER :: nxy0, myzbf0, myzef0

CALL Dmalloc0(resv, 1, 2, 1, nxy0, myzbf0 - 1, myzef0 + 1)
CALL Dmalloc0(resvhat, 1, 2, 1, nxy0, myzbf0 - 1, myzef0 + 1)
CALL Dmalloc0(p, 1, 2, 1, nxy0, myzbf0 - 1, myzef0 + 1)
CALL Dmalloc0(t, 1, 2, 1, nxy0, myzbf0 - 1, myzef0 + 1)
CALL Dmalloc0(s, 1, 2, 1, nxy0, myzbf0 - 1, myzef0 + 1)
CALL Dmalloc0(v, 1, 2, 1, nxy0, myzbf0 - 1, myzef0 + 1)
CALL Dmalloc0(vy, 1, 2, 1, nxy0, myzbf0 - 1, myzef0 + 1)
CALL Dmalloc0(vz, 1, 2, 1, nxy0, myzbf0 - 1, myzef0 + 1)
CALL Dmalloc0(temp, 1, 2, 1, nxy0, myzbf0 - 1, myzef0 + 1)
END SUBROUTINE

SUBROUTINE MatVecOP2G(x, y)
#ifdef MPI_ENV
USE MAT2x2OP
USE MPIComm_Mod, ONLY : GetNeighDat,  GetNeighDatFast
#endif
INTEGER :: ixy, ineigh, ibd, iz, iz0
REAL, POINTER :: X(:, :, :), Y(:, :, :)
REAL :: LMNT(2)
#ifndef MPI_ENV
DO iz = myzbf, myzef
  iz0 = LS%AxialPlaneMap(iz)
  DO ixy = 1, nxy
    !lmnt = LS%Diag(ixy, iz) * x(ixy, iz)
    lmnt = SUBMATVECOP(LS%Diag2g(:, ixy, iz), x(:, ixy, iz))
    DO ibd = 1, nbd
      ineigh = LS%NeighIdx(ibd, ixy)
      IF(ineigh .LE. 0) CYCLE
      !lmnt = lmnt + Ls%RadOffDiag(ibd, ixy, iz0) * X(ineigh, iz)
      lmnt = lmnt + SUBMATVECOP(LS%RadOffDiag2g(:, ibd, ixy, iz), x(:, ineigh, iz))
    ENDDO
    !Axial Contribution
    lmnt = lmnt + SUBMATVECOP(Ls%AxOffDiag(2, ixy, iz), X(:, ixy, iz + 1))  !--- BYS edit / 160223
    lmnt = lmnt + SUBMATVECOP(Ls%AxOffDiag(1, ixy, iz), X(:, ixy, iz - 1))  !--- BYS edit / 160223
    Y(:, ixy, iz) = lmnt
  ENDDO
ENDDO
#else

#ifdef FastComm
IF(nproc .GT. 1) CALL  GetNeighDatFast(X(1:2,1:nxy, myzbf-1:myzef+1), 2, nxy, myzbf, &
                                       myzef, myrank, nproc, comm, 0)
#endif
DO iz = myzbf, myzef
  iz0 = LS%AxialPlaneMap(iz)
  DO ixy = 1, nxy
    !lmnt = LS%Diag(ixy, iz) * x(ixy, iz)
    lmnt = SUBMATVECOP(LS%Diag2g(:, ixy, iz), x(:, ixy, iz))
    DO ibd = 1, nbd
      ineigh = LS%NeighIdx(ibd, ixy)
      IF(ineigh .LE. 0) CYCLE
      !lmnt = lmnt + Ls%RadOffDiag(ibd, ixy, iz0) * X(ineigh, iz)
      lmnt = lmnt + SUBMATVECOP(LS%RadOffDiag2g(:, ibd, ixy, iz0), x(:, ineigh, iz))
    ENDDO
    Y(:, ixy, iz) = lmnt
  ENDDO
ENDDO
IF(nproc .GT. 1) THEN
#ifndef FastComm
  CALL GetNeighDat(X(1:2, 1:nxy, myzbf-1:myzef+1), 2, nxy, myzbf, &
                        myzef, myrank, nproc, comm)
#else
 CALL  GetNeighDatFast(X(1:2, 1:nxy, myzbf-1:myzef+1), 2, nxy, myzbf, &
                        myzef, myrank, nproc, comm, 1)
#endif
ENDIF

DO iz = myzbf, myzef
  DO ixy = 1, nxy
    lmnt = SUBMATVECOP(Ls%AxOffDiag2g(:, 2, ixy, iz), X(:, ixy, iz + 1))
    lmnt = lmnt + SUBMATVECOP(Ls%AxOffDiag2g(:, 1, ixy, iz), X(:, ixy, iz - 1))
    Y(:, ixy, iz) = Y(:, ixy, iz) + lmnt
  ENDDO
ENDDO
#endif
END SUBROUTINE

FUNCTION DotProduct2G(x, y)
#ifdef MPI_ENV
USE MpiComm_mod, ONLY : REDUCE
#endif
IMPLICIT NONE
INTEGER :: ixy, iz
REAL, POINTER :: X(:, :, :), Y(:, :, :)
REAL :: DotProduct2G
#ifdef MPI_ENV
REAL :: PDotProduct
#endif
DotProduct2G = 0
DO iz = myzbf, myzef
  DO ixy = 1, nxy
    DotProduct2G = DotProduct2G + X(1, ixy, iz) * Y(1, ixy, iz)+ X(2, ixy, iz) * Y(2, ixy, iz)
  ENDDO
ENDDO
#ifdef MPI_ENV
PDotProduct = DotProduct2G
CALL REDUCE(PDotProduct, DotProduct2G, comm, .TRUE.)
#endif
END FUNCTION

END MODULE
