#include <defines.h>
#define pre
!#define bilu_pre
!#define ssor_pre
#define FastComm
!#define FastComm
MODULE BiCGSTAB_mod
USE TYPEDEF, ONLY : CMFDLS_TYPE
USE itrcntl_mod, ONLY : InSolverItrCntl_TYPE
USE BasicOperation, ONLY : CP_CA, CP_VA, SUB_VA
USE OMP_LIB
USE CNTL,           ONLY : nTracerCntl
USE BiLU_MOD,       ONLY : FreeBiluSolverEnv,     SetBiluSolverEnv,     SolveBilu,     &
                           SolveBiLU_OMP
                           
IMPLICIT NONE
TYPE(CMFDLS_TYPE),POINTER, PRIVATE :: LS
REAL, POINTER, PRIVATE :: resv(:, :), resvhat(:, :)
REAL, POINTER, PRIVATE :: p(:, :), t(:, :), s(:, :), v(:, :), vy(:, :), vz(:, :)
REAL, POINTER, PRIVATE :: temp(:, :)
INTEGER, PRIVATE :: myzbf, myzef, nxy, nbd
#ifdef MPI_ENV
INTEGER, PRIVATE :: comm, myrank, nproc
#endif
INTEGER, PRIVATE :: nThread
REAL, PRIVATE :: GlobalBuf(0:100)
#ifdef pre
LOGICAL, PRIVATE :: lssor = .FALSE.
LOGICAL, PRIVATE :: lbilu = .TRUE.
#endif
CONTAINS 

SUBROUTINE AllocBiCGSTAB(nxy0, myzbf0, myzef0)
USE ALLOCS
IMPLICIT NONE 
INTEGER :: nxy0, myzbf0, myzef0

CALL Dmalloc0(resv, 1, nxy0, myzbf0 - 1, myzef0 + 1)
CALL Dmalloc0(resvhat, 1, nxy0, myzbf0 - 1, myzef0 + 1)
CALL Dmalloc0(p, 1, nxy0, myzbf0 - 1, myzef0 + 1)
CALL Dmalloc0(t, 1, nxy0, myzbf0 - 1, myzef0 + 1)
CALL Dmalloc0(s, 1, nxy0, myzbf0 - 1, myzef0 + 1)
CALL Dmalloc0(v, 1, nxy0, myzbf0 - 1, myzef0 + 1)
CALL Dmalloc0(vy, 1, nxy0, myzbf0 - 1, myzef0 + 1)
CALL Dmalloc0(vz, 1, nxy0, myzbf0 - 1, myzef0 + 1)
CALL Dmalloc0(temp, 1, nxy0, myzbf0 - 1, myzef0 + 1)
END SUBROUTINE

SUBROUTINE BiCGSTAB(A, Sol, RHS, itrcntl)
USE PARAM
#ifdef MPI_ENV
USE MPIComm_Mod, ONLY : GetNeighDat
#endif
USE TIMER,          only : nTracer_dclock, TimeChk
IMPLICIT NONE
TYPE(CMFDLS_TYPE), TARGET :: A
REAL, POINTER :: Sol(:, :), RHS(:, :)
TYPE(InSolverItrCntl_TYPE) :: itrcntl
INTEGER :: iter, itermax, itermin
INTEGER :: ixy, iz
INTEGER :: i, j, k, n
REAL :: rho,prho,alpha,omega,gamma,beta
REAL :: reserr, reserr0, temp, convcrit
REAL :: Tbeg, Tend
LOGICAL :: lconv


IF (A%nThread .GT. 1) THEN
!  CALL BiCGSTAB_OMP(A, Sol, RHS, itrcntl)
!  RETURN
ENDIF

TBeg = nTracer_dclock(FALSE, FALSE)

LS => A;
#ifdef pre 
IF (.NOT. nTracerCntl%lHex) THEN ! KSC edit
  CALL SetBiluSolverEnv(A)
END IF
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
convcrit = ItrCntl%convcrit

!Initialized Solution related Array set
rho = 1._8; prho = 1._8; alpha = 1._8; omega = 1._8;
CALL CP_CA(P(1:nxy, myzbf - 1:myzef + 1), 0._8, nxy, myzef - myzbf + 3)
CALL CP_CA(v(1:nxy, myzbf - 1:myzef + 1), 0._8, nxy, myzef - myzbf + 3)

CALL Residual(Sol, RHS, ResV)

reserr0 = DotProduct(resv, resv)
reserr0 = sqrt(reserr0)
!
IF(reserr0 .lt. 1.0e-8_8) then
  reserr0 = 1._8
  itermax = 1
  RETURN
ENDIF

CALL CP_VA(resvhat(1:nxy, myzbf - 1:myzef + 1), RESV(1:nxy, myzbf - 1:myzef + 1), nxy, myzef - myzbf + 3)
!CALL CP_VA(resvhat(1:nxy, myzbf - 1:myzef + 1), RESV(1:nxy, myzbf - 1:myzef + 1), nxy, myzef - myzbf + 3)

lconv = .FALSE.
#define bicgstab2
DO iter = 1, itermax
  prho = rho
  rho = DotProduct(resvhat, resv)
  beta = (rho*alpha)/(prho*omega)
  !P vector Update
  DO iz = myzbf, myzef
    DO ixy = 1, nxy
      P(ixy, iz) = resv(ixy, iz) + beta * (P(ixy, iz) - omega * v(ixy, iz))
    ENDDO
  ENDDO
#ifdef pre

  IF (nTracerCntl%lHex) THEN ! KSC edit
    CALL MatVecOP(p, v) 
  ELSE
    CALL SolveBilu(p, vy)
    CALL MatVecOP(vy, v)
  END IF
#else
  CALL MatVecOP(p, v)  
#endif  

  !Update Gamma  = 
  gamma = DOtProduct(resvhat, v)
  alpha = rho/gamma
  !Update s = r(i-1)-alpha*v
  DO iz = myzbf, myzef
    DO ixy = 1, nxy
      s(ixy, iz) = resv(ixy, iz) - alpha * v(ixy, iz)
    ENDDO
  ENDDO
#ifdef pre  

  IF (nTracerCntl%lHex) THEN ! KSC edit
    CALL MatVecOP(s, t)
  ELSE
    CALL SolveBilu(s, vz)
    Call MatVecOP(vz, t)
  END IF
#else
  CALL MatVecOP(s, t)
#endif  
  !Update Omega = <t, s>/ <t, t>
  omega = DotProduct(t, s); temp = DotProduct(t, t)
!  IF(temp .NE. 0) omega = omega/temp
  IF(temp .lt. 1.0E-20_8) THEN
    DO iz = myzbf, myzef
      DO ixy = 1, nxy
#ifdef pre    
        IF (nTracerCntl%lHex) THEN ! KSC edit
          sol(ixy, iz) = sol(ixy, iz) + alpha * p(ixy, iz)
        ELSE
          sol(ixy, iz) = sol(ixy, iz) + alpha * vy(ixy, iz)
        END IF
#else
        sol(ixy, iz) = sol(ixy, iz) + alpha * p(ixy, iz)
#endif
      ENDDO
    ENDDO    
    CALL Residual(Sol, RHS, ResV)
    reserr = DotProduct(resv, resv)
    reserr = sqrt(reserr)
    CONTINUE    
    EXIT
  ELSE
    omega = omega/temp
  ENDIF
  
  !Update r = s - omega * t
  DO iz = myzbf, myzef
    DO ixy = 1, nxy
      resv(ixy, iz) = s(ixy, iz) - omega * t(ixy, iz)
    ENDDO
  ENDDO
  
  reserr = DotProduct(resv, resv)
  reserr = sqrt(reserr)
  !IF(myrank .eq. 0) PRINT *, alpha, omega, reserr    
  !Solution Update
  DO iz = myzbf, myzef
    DO ixy = 1, nxy
#ifdef pre    
      IF (nTracerCntl%lHex) THEN ! KSC edit
        sol(ixy, iz) = sol(ixy, iz) + alpha * p(ixy, iz) + omega * s(ixy, iz)
      ELSE
        sol(ixy, iz) = sol(ixy, iz) + alpha * vy(ixy, iz) + omega * vz(ixy, iz)
      END IF
#else
      sol(ixy, iz) = sol(ixy, iz) + alpha * p(ixy, iz) + omega * s(ixy, iz)
#endif      
    ENDDO
  ENDDO

  !write(99, *) iter, reserr
  IF((reserr/reserr0) .lt. convcrit ) lconv = .TRUE.
  IF(iter .GT. itermin .AND. lCONV) EXIT
  !IF(reserr .lt. 1.0E-9_8) EXIT
  !IF(lCONV) EXIT
ENDDO
!
#ifdef MPI_ENV
IF(nproc .GT. 1) THEN
  CALL GetNeighDat(Sol(1:nxy, myzbf-1:myzef+1), nxy, myzbf, myzef, myrank, nproc, comm)
ENDIF
#endif
CALL FreeBiluSolverEnv()
Tend = nTracer_dclock(FALSE, FALSE)
TimeChk%AxBTime = TimeChk%AxBTime + (Tend - Tbeg)
END SUBROUTINE

SUBROUTINE MatVecOP(x, y)
#ifdef MPI_ENV
USE MPIComm_Mod, ONLY : GetNeighDat,  GetNeighDatFast
#endif
INTEGER :: ixy, ineigh, ibd, iz, iz0
REAL, POINTER :: X(:, :), Y(:, :)
REAL :: LMNT

#ifndef MPI_ENV
DO iz = myzbf, myzef
  iz0 = LS%AxialPlaneMap(iz)
  DO ixy = 1, nxy
    lmnt = LS%Diag(ixy, iz) * x(ixy, iz)
    DO ibd = 1, nbd
      ineigh = LS%NeighIdx(ibd, ixy)
      IF(ineigh .LE. 0) CYCLE
      lmnt = lmnt + Ls%RadOffDiag(ibd, ixy, iz0) * X(ineigh, iz)
    ENDDO
    !Axial Contribution
    lmnt = lmnt + Ls%AxOffDiag(2, ixy, iz) * X(ixy, iz + 1)
    lmnt = lmnt + Ls%AxOffDiag(1, ixy, iz) * X(ixy, iz - 1)
    Y(ixy, iz) = lmnt
  ENDDO
ENDDO
#else

#ifdef FastComm
IF(nproc .GT. 1) CALL  GetNeighDatFast(X(1:nxy, myzbf-1:myzef+1), nxy, myzbf, &
                                       myzef, myrank, nproc, comm, 0)
#endif
DO iz = myzbf, myzef
  iz0 = LS%AxialPlaneMap(iz)
  DO ixy = 1, nxy
    lmnt = LS%Diag(ixy, iz) * x(ixy, iz)
    DO ibd = 1, nbd
      ineigh = LS%NeighIdx(ibd, ixy)
      IF(ineigh .LE. 0) CYCLE
      lmnt = lmnt + Ls%RadOffDiag(ibd, ixy, iz0) * X(ineigh, iz)
    ENDDO
    Y(ixy, iz) = lmnt
  ENDDO
ENDDO

IF(nproc .GT. 1) THEN
#ifndef FastComm
  CALL GetNeighDat(X(1:nxy, myzbf-1:myzef+1), nxy, myzbf, &
                        myzef, myrank, nproc, comm)
#else
 CALL  GetNeighDatFast(X(1:nxy, myzbf-1:myzef+1), nxy, myzbf, &
                        myzef, myrank, nproc, comm, 1)
#endif                        
ENDIF

DO iz = myzbf, myzef
  DO ixy = 1, nxy
    lmnt = Ls%AxOffDiag(2, ixy, iz) * X(ixy, iz + 1)
    lmnt = lmnt + Ls%AxOffDiag(1, ixy, iz) * X(ixy, iz - 1)    
    Y(ixy, iz) = Y(ixy, iz) + lmnt
  ENDDO
ENDDO
#endif
END SUBROUTINE
FUNCTION DotProduct(x, y)
#ifdef MPI_ENV
USE MpiComm_mod, ONLY : REDUCE 
#endif
IMPLICIT NONE
INTEGER :: ixy, iz
REAL, POINTER :: X(:, :), Y(:, :)
REAL :: DotProduct
#ifdef MPI_ENV
REAL :: PDotProduct
#endif
DotProduct = 0
DO iz = myzbf, myzef
  DO ixy = 1, nxy
    DotProduct = DotProduct + X(ixy, iz) * Y(ixy, iz)
  ENDDO
ENDDO
#ifdef MPI_ENV
PDotProduct = DotProduct
CALL REDUCE(PDotProduct, DotProduct, comm, .TRUE.)
#endif
END FUNCTION

SUBROUTINE Residual(x, y, residualVec)
IMPLICIT NONE
REAL, POINTER :: x(:, :), y(:, :) , residualVec(:, :)
INTEGER :: ixy, iz
CALL MatVecOP(x, residualVec)
DO iz = myzbf, myzef
  DO ixy = 1, nxy
    ResidualVec(ixy, iz) = y(ixy, iz) - ResidualVec(ixy, iz)
  ENDDO
ENDDO  
END SUBROUTINE

SUBROUTINE Minv(b, x)
INTEGER :: ixy, ineigh, ibd, iz, iz0
REAL, POINTER :: B(:, :), X(:, :)
REAL :: LMNT
REAL :: w

w = 1.4
DO iz = myzbf, myzef
  iz0 = LS%AxialPlaneMap(iz)
  DO ixy = 1, nxy
     temp(ixy, iz) = B(ixy, iz)
!    lmnt = LS%Diag(ixy, iz) * x(ixy, iz)
    lmnt = 0
    DO ibd = 1, nbd
      ineigh = LS%NeighIdx(ibd, ixy)
      IF(ineigh .LE. 0) CYCLE
      IF(ineigh .GT. ixy) CYCLE
      lmnt = lmnt +  Ls%RadOffDiag(ibd, ixy, iz0) * temp(ineigh, iz)
!      lmnt = lmnt + Ls%RadOffDiag(ibd, ixy, iz0) * X(ineigh, iz)
    ENDDO
!    !Axial Contribution
     IF(iz .NE. myzbf) THEN
      lmnt = lmnt + Ls%AxOffDiag(1, ixy, iz) * temp(ixy, iz - 1)
     ENDIF

     Temp(ixy, iz) = (Temp(ixy, iz) - w * lmnt) / LS%Diag(ixy, iz)
  ENDDO
ENDDO

DO iz = myzef, myzbf, -1
  iz0 = LS%AxialPlaneMap(iz)
  DO ixy = nxy, 1, -1
     temp(ixy, iz) = B(ixy, iz)
    lmnt = 0
    DO ibd = 1, nbd
      ineigh = LS%NeighIdx(ibd, ixy)
      IF(ineigh .LE. 0) CYCLE
      IF(ineigh .LT. ixy) CYCLE
      lmnt = lmnt +  Ls%RadOffDiag(ibd, ixy, iz0) * x(ineigh, iz)
    ENDDO
     IF(iz .NE. myzef) THEN
       lmnt = lmnt + Ls%AxOffDiag(2, ixy, iz) * x(ixy, iz + 1)
     ENDIF
     x(ixy, iz) = temp(ixy, iz) - lmnt / LS%Diag(ixy, iz)
  ENDDO
ENDDO

END SUBROUTINE

FUNCTION DotProduct_OMP(x, y, idom)
USE MPIComm_Mod,       ONLY : REDUCE
USE LsRadDcmp_mod,     ONLY : RadDcmp
IMPLICIT NONE
INTEGER :: ixy, iz, idom, i, nxylocal
INTEGER, POINTER :: list(:)
REAL, POINTER :: X(:, :), Y(:, :)
REAL :: DotProduct_OMP, PDotProduct


DotProduct_OMP = 0
nxylocal = RadDcmp%nxylocal(idom)
list => RadDcmp%PinIdx(:, idom)
!$OMP Barrier
DO iz = myzbf, myzef
  DO i = 1, nxylocal
    ixy = list(i)
    DotProduct_OMP = DotProduct_OMP + X(ixy, iz) * Y(ixy, iz)
  ENDDO
ENDDO
Globalbuf(idom) = DotProduct_OMP
!$OMP Barrier
!$OMP MASTER
PDotProduct = sum(GlobalBuf(1:nthread)); DotProduct_OMP = PDotProduct
CALL REDUCE(PDotProduct, DotProduct_OMP, comm, .TRUE.)
Globalbuf(0)=DotProduct_OMP
!$OMP END MASTER
!$OMP Barrier
DotProduct_OMP = GlobalBuf(0)
END FUNCTION


SUBROUTINE MatVecOP_OMP(x, y, idom)
#ifdef MPI_ENV
USE MPIComm_Mod,    ONLY : GetNeighDat,  GetNeighDatFast
#endif
USE LsRadDcmp_MOD,  ONLY : RadDcmp
INTEGER :: idom
REAL, POINTER :: X(:, :), Y(:, :)

INTEGER :: ixy, ineigh, ibd, iz, iz0, i, nxylocal
INTEGER, POINTER :: list(:)
REAL :: LMNT

nxylocal = RadDcmp%nxylocal(idom)
List => RadDcmp%PinIdx(:, idom)
#ifndef MPI_ENV
DO iz = myzbf, myzef
  iz0 = LS%AxialPlaneMap(iz)
  DO i = 1, nxylocal
    ixy = list(i)
    lmnt = LS%Diag(ixy, iz) * x(ixy, iz)
    DO ibd = 1, nbd
      ineigh = LS%NeighIdx(ibd, ixy)
      IF(ineigh .LE. 0) CYCLE
      lmnt = lmnt + Ls%RadOffDiag(ibd, ixy, iz0) * X(ineigh, iz)
    ENDDO
    !Axial Contribution
    lmnt = lmnt + Ls%AxOffDiag(2, ixy, iz) * X(ixy, iz + 1)
    lmnt = lmnt + Ls%AxOffDiag(1, ixy, iz) * X(ixy, iz - 1)
    Y(ixy, iz) = lmnt
  ENDDO
ENDDO
!$OMP BARRIER
#else


#ifdef FastComm
!$OMP BARRIER
!$OMP MASTER
IF(nproc .GT. 1) CALL  GetNeighDatFast(X(1:nxy, myzbf-1:myzef+1), nxy, myzbf, &
                                       myzef, myrank, nproc, comm, 0)
!$OMP END MASTER
!$OMP BARRIER
#endif

DO iz = myzbf, myzef
  iz0 = LS%AxialPlaneMap(iz)
  DO i = 1, nxylocal
    ixy = list(i)
    lmnt = LS%Diag(ixy, iz) * x(ixy, iz)
    DO ibd = 1, nbd
      ineigh = LS%NeighIdx(ibd, ixy)
      IF(ineigh .LE. 0) CYCLE
      lmnt = lmnt + Ls%RadOffDiag(ibd, ixy, iz0) * X(ineigh, iz)
    ENDDO
    Y(ixy, iz) = lmnt
  ENDDO
ENDDO
!$OMP BARRIER
IF(nproc .GT. 1) THEN
!$OMP MASTER
#ifndef FastComm
  CALL GetNeighDat(X(1:nxy, myzbf-1:myzef+1), nxy, myzbf, &
                        myzef, myrank, nproc, comm)
#else
 CALL  GetNeighDatFast(X(1:nxy, myzbf-1:myzef+1), nxy, myzbf, &
                        myzef, myrank, nproc, comm, 1)
#endif      
!$OMP END MASTER    
!$OMP BARRIER
ENDIF

DO iz = myzbf, myzef
  DO i = 1, nxylocal
    ixy = list(i)
    lmnt = Ls%AxOffDiag(2, ixy, iz) * X(ixy, iz + 1)
    lmnt = lmnt + Ls%AxOffDiag(1, ixy, iz) * X(ixy, iz - 1)    
    Y(ixy, iz) = Y(ixy, iz) + lmnt
  ENDDO
ENDDO
!$OMP BARRIER
#endif

END SUBROUTINE


SUBROUTINE BiCGSTAB_OMP(A, Sol, RHS, itrcntl)
USE PARAM
#ifdef MPI_ENV
USE MPIComm_Mod, ONLY : GetNeighDat, REDUCE
#endif

USE LsRadDcmp_MOD, ONLY : RadDcmp
USE TIMER,          only : nTracer_dclock, TimeChk
IMPLICIT NONE
TYPE(CMFDLS_TYPE), TARGET :: A
REAL, POINTER :: Sol(:, :), RHS(:, :)
TYPE(InSolverItrCntl_TYPE) :: itrcntl
INTEGER :: iter, itermax, itermin, tid
INTEGER :: ixy, iz
INTEGER :: i, j, k, n
REAL :: rho,prho,alpha,omega,gamma,beta
REAL :: reserr, reserr0, temp
REAL :: tmp, tmp2
REAL :: tbeg, tend
LOGICAL :: lconv

INTEGER :: nxylocal
INTEGER, POINTER :: List(:, :)

tbeg = nTracer_dclock(FALSE, FALSE)
LS => A;
#ifdef pre 
CALL SetBiluSolverEnv(A)
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
CALL CP_CA(P(1:nxy, myzbf - 1:myzef + 1), 0._8, nxy, myzef - myzbf + 3)
CALL CP_CA(v(1:nxy, myzbf - 1:myzef + 1), 0._8, nxy, myzef - myzbf + 3)

CALL Residual(Sol, RHS, ResV)

reserr0 = DotProduct(resv, resv)
reserr0 = sqrt(reserr0)

IF(reserr0 .lt. 1.0e-8_8) then
  reserr0 = 1._8
  itermax = 1
  RETURN
ENDIF

CALL CP_VA(resvhat(1:nxy, myzbf - 1:myzef + 1), RESV(1:nxy, myzbf - 1:myzef + 1), nxy, myzef - myzbf + 3)
!CALL CP_VA(resvhat(1:nxy, myzbf - 1:myzef + 1), RESV(1:nxy, myzbf - 1:myzef + 1), nxy, myzef - myzbf + 3)

lconv = .FALSE.
nThread = A%nThread
!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(nThread) 

list => RadDcmp%PinIdx
!t1 = OMP_GET_WTIME()
!$OMP PARALLEL DEFAULT(SHARED)      &
!!$OMP PRIVATE(tid, iter, tmp, tmp2, iz, ixy, i, nxylocal)
!$OMP PRIVATE(tid, iter, tmp, tmp2, iz, ixy, i, nxylocal, rho, prho, beta, alpha, gamma, omega)
!$  tid = omp_get_thread_num()+1
nxylocal = RadDcmp%nxylocal(tid)
rho = 1._8; prho = 1._8; alpha = 1._8; omega = 1._8;
DO iter = 1, itermax
  prho = rho
  
  tmp = DotProduct_OMP(resvhat, resv, tid)
  rho = tmp
  beta = (rho*alpha)/(prho*omega)
!P vector Update
  DO iz = myzbf, myzef
    DO i = 1, nxylocal
      ixy = list(i,tid)
      P(ixy, iz) = resv(ixy, iz) + beta * (P(ixy, iz) - omega * v(ixy, iz))
    ENDDO
  ENDDO
!$OMP BARRIER
#ifdef pre
  CALL SolveBilu_OMP(p, vy, tid)
!!$OMP MASTER
!  CALL SolveBilu(p, vy)
!!$OMP END MASTER
  CALL MatVecOP_OMP(vy, v, tid)
#else
  CALL MatVecOP_OMP(p, v, tid)  
#endif    

  tmp = DOtProduct_OMP(resvhat, v, tid)
  gamma = tmp
  alpha = rho/gamma
  DO iz = myzbf, myzef
    DO i = 1, nxylocal
      ixy = list(i, tid)
      s(ixy, iz) = resv(ixy, iz) - alpha * v(ixy, iz)
    ENDDO
  ENDDO
#ifdef pre  
  CALL SolveBilu_OMP(s, vz, tid)
  Call MatVecOP_OMP(vz, t, tid)
#else
  CALL MatVecOP_OMP(s, t, tid)
#endif 
!Update Omega = <t, s>/ <t, t>
  tmp  = DotProduct_OMP(t, s, tid)
  tmp2 = DotProduct_OMP(t, t, tid)

  omega = tmp; temp = tmp2
  IF(temp .NE. 0) omega = omega/temp
!  IF(temp .lt. 1.0E-15_8) THEN
!    !print *, tid
!    DO iz = myzbf, myzef
!      DO i = 1, nxylocal
!        ixy = list(i, tid)
!#ifdef pre    
!        sol(ixy, iz) = sol(ixy, iz) + alpha * vy(ixy, iz)
!#else
!        sol(ixy, iz) = sol(ixy, iz) + alpha * p(ixy, iz)
!#endif
!      ENDDO
!    ENDDO   
!!$OMP BARRIER 
!!$OMP Master
!    CALL Residual(Sol, RHS, ResV)
!    reserr = DotProduct(resv, resv)
!    reserr = sqrt(reserr)
!!$OMP END Master 
!!$OMP BARRIER
!    CONTINUE    
!    EXIT
!  ENDIF
  !Update r = s - omega * t
  DO iz = myzbf, myzef
    DO i = 1, nxylocal
      ixy = list(i, tid)
      resv(ixy, iz) = s(ixy, iz) - omega * t(ixy, iz)
    ENDDO
  ENDDO

  tmp = DotProduct_OMP(resv, resv, tid)
!$OMP MASTER
  reserr = sqrt(tmp)
!$OMP END MASTER
!$OMP BARRIER   
  DO iz = myzbf, myzef
    DO i = 1, nxylocal
      ixy = list(i, tid)
#ifdef pre    
      sol(ixy, iz) = sol(ixy, iz) + alpha * vy(ixy, iz) + omega * vz(ixy, iz)
#else
      sol(ixy, iz) = sol(ixy, iz) + alpha * p(ixy, iz) + omega * s(ixy, iz)
#endif      
    ENDDO
  ENDDO
!$OMP BARRIER
  !IF(tid .EQ. 1) WRITE(99, *) iter, reserr
  IF((reserr/reserr0) .lt. itrcntl%convcrit ) lconv = .TRUE.
  IF(iter .GT. itermin .AND. lCONV) EXIT
  !IF(reserr .lt. 1.0E-9_8) EXIT
  IF(lCONV) EXIT
ENDDO
!$OMP END PARALLEL

#ifdef MPI_ENV
IF(nproc .GT. 1) THEN
  CALL GetNeighDat(Sol(1:nxy, myzbf-1:myzef+1), nxy, myzbf, myzef, myrank, nproc, comm)
ENDIF
#endif
CALL FreeBiluSolverEnv()
Tend = nTracer_dclock(FALSE, FALSE)
TimeChk%AxBTime = TimeChk%AxBTime + (Tend - Tbeg)
!PRINT *, (Tend - Tbeg)
!STOP
END SUBROUTINE
END MODULE

