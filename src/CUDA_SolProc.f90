#include <Depletion.h>
! Solve System
!   : Predictor/Corrector, Gd QD routines
! Post Process
!   : Burnup for each fxr, Final ND afte PC, Gd Post-correction, [Xe,Sm] Eq/Tr
#ifdef __PGI
SUBROUTINE cuSolveDeplSys(DeplSysBundle, DeplFxrBundle, lCorrector, lGd, SolverTyp, Nths)
! *-------------------- Made by LHG, Aug. 2018 ---------------------* !
!   This funtion is designed to utilize iterative CRAM solver. Also,  !
!  batching the systems is included.                                  !
! *-----------------------------------------------------------------* !
USE HPDeplType
USE CSRMATRIX
USE OMP_LIB
USE CUDA_MASTER
USE CUDA_Workspace, ONLY : NNZ, NR, AllocBatchMem, CopyOutBatchMem
USE cuMatExponential, ONLY : MatExpCRAM_Iter
USE Timer,          ONLY : nTracer_dclock, TimeChk
IMPLICIT NONE

TYPE(DeplSysBundle_Type) :: DeplSysBundle
TYPE(DeplFxrBundle_Type) :: DeplFxrBundle
LOGICAL :: lCorrector, lGd
INTEGER :: Nths, SolverTyp

TYPE(DeplFxr_Type), POINTER :: Fxrs(:)

REAL(8), POINTER :: pnums(:,:)
REAL(8), POINTER :: pnums1(:)
REAL(8), POINTER :: pnum0(:), pnum1(:)

INTEGER, POINTER :: IdTrue(:)
INTEGER, ALLOCATABLE :: nnzs(:)

INTEGER :: NofIso
INTEGER :: Nsys, ifxrbeg
INTEGER :: i, j, k, l

REAL :: Tb, Te

Tb = nTracer_dclock(.FALSE., .FALSE.)
NofIso = DeplSysBundle%NofIso; Nsys = DeplSysBundle%Nsys; ifxrbeg = DeplSysBundle%ifxrbeg
pnums(1:NofIso,1:Nsys) => DeplSysBundle%pnums; pnums1 => DeplSysBundle%pnums1
Fxrs => DeplFxrBundle%FxrBundle

CALL AllocBatchMem(Nsys)
CALL CopyOutBatchMem(pnums,DeplSysBundle%Diag,DeplSysBundle%OffDiag)

Te = nTracer_dclock(.FALSE., .FALSE.)
TimeChk%DeplcuSysTime = TimeChk%DeplcuSysTime + (Te-Tb)
TimeChk%DeplSolTime = TimeChk%DeplSolTime - (Te-Tb)
CALL MatExpCRAM_Iter(pnums1,cuDevice%myStream)
!print*, pnums1(Nsys*343), pnums1(Nsys*342+1)
END SUBROUTINE

!SUBROUTINE cuSolveDeplSys(DeplSysBundle, DeplFxrBundle, lCorrector, lGd, SolverTyp, Nths)
!! *-------------------- Made by LHG, Aug. 2018 ---------------------* !
!!   This funtion is designed to utilize iterative CRAM solver. Also,  !
!!  batching the systems is included.                                  !
!! *-----------------------------------------------------------------* !
!USE HPDeplType
!USE CSRMATRIX
!USE OMP_LIB
!USE CUDA_MASTER
!USE CUDA_Workspace, ONLY : NNZ, NR, AllocBatchMem, CopyOutBatchMem
!USE cuMatExponential, ONLY : MatExpCRAM_Iter
!USE Timer,          ONLY : nTracer_dclock, TimeChk
!IMPLICIT NONE
!
!TYPE(DeplSysBundle_Type) :: DeplSysBundle
!TYPE(DeplFxrBundle_Type) :: DeplFxrBundle
!LOGICAL :: lCorrector, lGd
!INTEGER :: Nths, SolverTyp
!
!TYPE(DeplFxr_Type), POINTER :: Fxrs(:)
!
!REAL(8), POINTER :: pnums(:)
!REAL(8), POINTER :: pnums1(:)
!REAL(8), POINTER :: pnum0(:), pnum1(:)
!
!COMPLEX(8), ALLOCATABLE, SAVE :: val(:,:)
!REAL(8), ALLOCATABLE, SAVE :: vec(:,:)
!TYPE(CSR_DOUBLE_COMPLEX) :: DeplMat
!
!INTEGER, POINTER :: IdTrue(:)
!INTEGER, ALLOCATABLE :: nnzs(:)
!
!INTEGER :: NofIso
!INTEGER :: Nsys, ifxrbeg
!INTEGER :: i, j, k, l
!
!REAL :: Tb, Te
!
!Tb = nTracer_dclock(.FALSE., .FALSE.)
!NofIso = DeplSysBundle%NofIso; Nsys = DeplSysBundle%Nsys; ifxrbeg = DeplSysBundle%ifxrbeg
!pnums => DeplSysBundle%pnums; pnums1 => DeplSysBundle%pnums1
!!DeplMats => DeplSysBundle%DeplMats;
!Fxrs => DeplFxrBundle%FxrBundle
!!IF (.NOT. ALLOCATED(Val)) ALLOCATE(Val(Nsys, nnz),vec(Nsys, NofIso))
!!call OMP_SET_NUM_THREADS(Nths);
!!!$OMP PARALLEL DO PRIVATE(j) SCHEDULE(GUIDED)
!!DO i = 1, Nsys
!!  Val(i,:) = DeplSysBundle%DeplMats(i)%csrVal(:)
!!  j = (i-1)*NofIso
!!  vec(i,:) = pnums(j+1:j+NofIso)
!!END DO
!!!$OMP END PARALLEL DO
!IF (.NOT. ALLOCATED(Val)) ALLOCATE(Val(nnz,Nsys),vec(NofIso,Nsys))
!call OMP_SET_NUM_THREADS(Nths);
!!$OMP PARALLEL DO PRIVATE(j) SCHEDULE(GUIDED)
!DO i = 1, Nsys
!  Val(:,i) = DeplSysBundle%DeplMats(i)%csrVal(:)
!  j = (i-1)*NofIso
!  vec(:,i) = pnums(j+1:j+NofIso)
!END DO
!!$OMP END PARALLEL DO
!CALL AllocBatchMem(Nsys)
!CALL CopyOutBatchMem(vec,Val)
!!print*, vec(Nsys, 343), vec(1, 343)
!!DEALLOCATE(Val,vec)
!Te = nTracer_dclock(.FALSE., .FALSE.)
!TimeChk%DeplcuSysTime = TimeChk%DeplcuSysTime + (Te-Tb)
!TimeChk%DeplSolTime = TimeChk%DeplSolTime - (Te-Tb)
!CALL MatExpCRAM_Iter(pnums1,cuDevice%myStream)
!!print*, pnums1(Nsys*343), pnums1(Nsys*342+1)
!END SUBROUTINE

#endif
