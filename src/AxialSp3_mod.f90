#include <SP3DEFINES.H>
MODULE Sp3Senm_Mod
USE TYPEDEF, ONLY : PE_TYPE, PinXS_Type, AxFLX_TYPE
USE SP3DEF
IMPLICIT NONE
LOGICAL :: lFirst = .TRUE.

REAL,PARAMETER, PRIVATE :: pi=3.1415926535897932384626433832795_8
REAL,PARAMETER, PRIVATE :: tor=1.0E-6
REAL,PARAMETER, PRIVATE :: infinite=1.0e+9
!REAL, PRIVATE :: Eigv
TYPE(AxSolItrCntl_TYPE), PRIVATE :: AxSolItrCntl
INTEGER, PRIVATE :: NG_
INTEGER, PRIVATE :: NMESH_

TYPE(ModalFlx_Type), POINTER, PRIVATE, SAVE :: MFLX_OMP(:, :)
TYPE(SenmMAT_Type), POINTER,  PRIVATE, SAVE :: SenmMat_OMP(:, :)
REAL, POINTER,  PRIVATE, SAVE :: COEFF_OMP(:, :, :), Q_OMP(:, :, :, :), RHS_OMP(:, :, :), PSI_OMP(:, :)
REAL, POINTER,  PRIVATE, SAVE :: PHI0_OMP(:, :, :, :)

CONTAINS

SUBROUTINE SetSP3SenmEnv(Geom, PE)
USE SenmMatSolver_Mod,   ONLY : Alloc_SENMMAT
USE ALLOCS
IMPLICIT NONE
TYPE(AxGeom_TYPE) :: Geom
TYPE(PE_TYPE) :: PE
INTEGER :: NTHREAD_
INTEGER :: i
NMESH_ = Geom%NMESH
NG_ = Geom%NG
NTHREAD_ = PE%NAXTHREAD
ALLOCATE(MFLX_OMP(NMESH_, NTHREAD_)); !CALL ALLOC_MFLX(MFLX, NMESH, NG)
ALLOCATE(SenmMAT_OMP(NG_, NTHREAD_)); !CALL ALLOC_SENMMAT(SenmMAT, NMESH, NG)
CALL Dmalloc(Coeff_OMP, 2, 2*NMESH_, NTHREAD_);  CALL Dmalloc(RHS_OMP, 2, 2*NMESH_, NTHREAD_)
CALL Dmalloc0(Q_OMP, 0, 4, 1, 2, 1, NMESH_, 1, NTHREAD_)
CALL Dmalloc0(PHI0_OMP, 1, 2, 1, NMESH_, 1, ng_, 1, NTHREAD_)
CALL Dmalloc(PSI_OMP, NMESH_, NTHREAD_)

DO i = 1, NTHREAD_
  CALL ALLOC_MFLX(MFLX_OMP(:, i), NMESH_, NG_)
  CALL Alloc_SENMMAT(SenmMat_OMP(:, i), NMESH_, NG_)
ENDDO
END SUBROUTINE

SUBROUTINE DirectSP3SENM(FLX, XS, EIGV0, GEOM, TID)
USE SENMOP_MOD,          ONLY : INIT_SP3SENM
USE SenmMatSolver_Mod,   ONLY : Alloc_SENMMAT,  SenmMatSolve_LU, &
                                Init_SenmMAT
USE DirectSENM_MOD,      ONLY : SetSENMMAT,     UpdtAxialPsi,   &
                                UpdtAxialQ,     UpdtAxialPSOL,  &
                                UpdtCurrent,    UpdtP1Dhat,     &
                                SetRHS,         UpdtFlux,       &
                                PrintFlux
USE ALLOCS
IMPLICIT NONE
TYPE(AXFLX_TYPE) :: FLX(NMESH_)
TYPE(PINXS_TYPE) :: XS(NMESH_)
TYPE(AxGeom_Type) :: Geom
REAL :: eigv0
INTEGER :: TID

TYPE(ModalFlx_Type), POINTER :: MFLX(:)
TYPE(SenmMAT_Type), POINTER :: SenmMat(:)
REAL, POINTER :: COEFF(:, :), Q(:, :, :), RHS(:, :), PSI(:)
REAL, POINTER :: PHI0(:, :, :)
!REAL, POINTER :: HZ(:)
!INTEGER, POINTER :: COMP(:)
INTEGER :: BC(2), NMESH, NG
INTEGER :: IT, IG, I
LOGICAL :: LCONV

NMESH = GEOM%NMESH; NG = GEOM%NG
BC = GEOM%BC

MFLX => MFLX_OMP(:, TID); SenmMat => SenmMat_OMP(:, TID)
Coeff => Coeff_OMP(:, :, TID); RHS => RHS_OMP(:, :, TID)
Q => Q_OMP(:, :, :, TID); PHI0 => PHI0_OMP(:, :, :, TID)
PSI => PSI_OMP(:, TID)
!IF(lFirst) THEN
!  ALLOCATE(MFLX(NMESH)); CALL ALLOC_MFLX(MFLX, NMESH, NG)
!  ALLOCATE(SenmMAT(NG)); CALL ALLOC_SENMMAT(SenmMAT, NMESH, NG)
!  CALL Dmalloc(Coeff, 2, 2*NMESH);  CALL Dmalloc(RHS, 2, 2*NMESH)
!  CALL Dmalloc0(Q, 0, 4, 1, 2, 1, NMESH)
!  CALL Dmalloc0(PHI0, 1, 2, 1, NMESH, 1, ng)
!  CALL Dmalloc(PSI, NMESH)
!  lFirst = .FALSE.
!ENDIF

CALL INIT_SP3SENM(FLX,XS,GEOM,NMESH, NG)
CALL Init_SenmMAT(SenmMAT, NMESH, NG)

CALL SetSENMMAT(SENMMAT, FLX, XS, GEOM, NMESH, NG)

DO IT = 1, 5
  IF(IT .EQ. 1) THEN
    CALL UpdtAxialPsi(PSI, FLX, XS, GEOM, NMESH, .TRUE.)
  ELSE
!    CALL UpdtAxialPsi(PSI, FLX, XS, GEOM, NMESH, .FALSE.)
  ENDIF
  DO IG = 1, NG
    !Update 4th order Flux Shape
    CALL UpdtAxialQ(Q, Flx, MFLX, XS, GEOM, EIGV0, NMESH, IG, NG, Geom%lTransient)
    IF(IT .EQ. 1) THEN
      DO I = 1, NMESH
        PHI0(1:2, I, ig) =  FLX(I)%PHI(0, 1:2, ig)
      ENDDO
    ENDIF

    !Update Particular Solutin terms
    CALL UpdtAxialPSOL(FLX, MFLX, XS, GEOM, NMESH, IG, NG)
    !Determine RHS of Linear System using particular solution
    CALL SetRHS(RHS, FLX, MFLX, XS, GEOM, NMESH, IG, NG)
    !Solve Linear System
    CALL SenmMatSolve_LU(SenmMat(IG), Coeff(:, :), RHS, 2*NMESH)
    !Update Average Moment and 4th order flux shape fluctions
    CALL UpdtFlux(Coeff(:, :), FLX, MFLX, XS, GEOM, NMESH, IG, NG)
  ENDDO
ENDDO
!Current Update at the node interfaces
CALL UpdtCurrent(FLX, MFLX, XS, GEOM, NMESH, NG)
DO I = 1, NMESH
 ! FLX(I)%PHI(0, 1:2, 1:NG) = PHI0(1:2, I, 1:NG)
ENDDO
CALL UpdtP1Dhat(FLX, MFLX, XS, GEOM, NMESH, NG)
!CALL PrintFlux(FLX, MFLX, GEOM, NMESH, 1, NG, NG)
END SUBROUTINE

SUBROUTINE DirectSP3SENM_old(FLX, XS, EIGV0, GEOM, TID)
USE SENMOP_MOD,          ONLY : INIT_SP3SENM
USE SenmMatSolver_Mod,   ONLY : Alloc_SENMMAT,  SenmMatSolve_LU, &
                                Init_SenmMAT
USE DirectSENM_MOD,      ONLY : SetSENMMAT,     UpdtAxialPsi,   &
                                UpdtAxialQ,     UpdtAxialPSOL,  &
                                UpdtCurrent,    UpdtP1Dhat,     &
                                SetRHS,         UpdtFlux,       &
                                PrintFlux
USE ALLOCS
IMPLICIT NONE
TYPE(AXFLX_TYPE) :: FLX(NMESH_), FLX0(nmesh_)
TYPE(PINXS_TYPE) :: XS(NMESH_)
TYPE(AxGeom_Type) :: Geom
REAL :: eigv0
INTEGER :: TID

TYPE(ModalFlx_Type), POINTER :: MFLX(:)
TYPE(SenmMAT_Type), POINTER :: SenmMat(:)
REAL, POINTER :: COEFF(:, :), Q(:, :, :), RHS(:, :), PSI(:)
REAL, POINTER :: PHI0(:, :, :)
!REAL, POINTER :: HZ(:)
!INTEGER, POINTER :: COMP(:)
INTEGER :: BC(2), NMESH, NG
INTEGER :: IT, IG, I
LOGICAL :: LCONV

NMESH = GEOM%NMESH; NG = GEOM%NG
BC = GEOM%BC

MFLX => MFLX_OMP(:, TID); SenmMat => SenmMat_OMP(:, TID)
Coeff => Coeff_OMP(:, :, TID); RHS => RHS_OMP(:, :, TID)
Q => Q_OMP(:, :, :, TID); PHI0 => PHI0_OMP(:, :, :, TID)
PSI => PSI_OMP(:, TID)
!IF(lFirst) THEN
!  ALLOCATE(MFLX(NMESH)); CALL ALLOC_MFLX(MFLX, NMESH, NG)
!  ALLOCATE(SenmMAT(NG)); CALL ALLOC_SENMMAT(SenmMAT, NMESH, NG)
!  CALL Dmalloc(Coeff, 2, 2*NMESH);  CALL Dmalloc(RHS, 2, 2*NMESH)
!  CALL Dmalloc0(Q, 0, 4, 1, 2, 1, NMESH)
!  CALL Dmalloc0(PHI0, 1, 2, 1, NMESH, 1, ng)
!  CALL Dmalloc(PSI, NMESH)
!  lFirst = .FALSE.
!ENDIF

CALL INIT_SP3SENM(FLX,XS,GEOM,NMESH, NG)
CALL Init_SenmMAT(SenmMAT, NMESH, NG)

CALL SetSENMMAT(SENMMAT, FLX, XS, GEOM, NMESH, NG)

DO IT = 1, 5
  IF(IT .EQ. 1) THEN
    CALL UpdtAxialPsi(PSI, FLX, XS, GEOM, NMESH, .TRUE.)
  ELSE
!    CALL UpdtAxialPsi(PSI, FLX, XS, GEOM, NMESH, .FALSE.)
  ENDIF
  DO IG = 1, NG
    !Update 4th order Flux Shape
    !CALL UpdtAxialQ(Q, Flx, MFLX, XS, GEOM, EIGV0, NMESH, IG, NG, Geom%lTransient)
  ENDDO
  DO IG = 1, NG
    !Update 4th order Flux Shape
    CALL UpdtAxialQ(Q, Flx, MFLX, XS, GEOM, EIGV0, NMESH, IG, NG, Geom%lTransient)
    IF(IT .EQ. 1) THEN
      DO I = 1, NMESH
        PHI0(1:2, I, ig) =  FLX(I)%PHI(0, 1:2, ig)
      ENDDO
    ENDIF

    !Update Particular Solutin terms
    CALL UpdtAxialPSOL(FLX, MFLX, XS, GEOM, NMESH, IG, NG)
    !Determine RHS of Linear System using particular solution
    CALL SetRHS(RHS, FLX, MFLX, XS, GEOM, NMESH, IG, NG)
    !Solve Linear System
    CALL SenmMatSolve_LU(SenmMat(IG), Coeff(:, :), RHS, 2*NMESH)
    !Update Average Moment and 4th order flux shape fluctions
    CALL UpdtFlux(Coeff(:, :), FLX, MFLX, XS, GEOM, NMESH, IG, NG)
  ENDDO
ENDDO
!Current Update at the node interfaces
CALL UpdtCurrent(FLX, MFLX, XS, GEOM, NMESH, NG)
DO I = 1, NMESH
 ! FLX(I)%PHI(0, 1:2, 1:NG) = PHI0(1:2, I, 1:NG)
ENDDO
CALL UpdtP1Dhat(FLX, MFLX, XS, GEOM, NMESH, NG)
!CALL PrintFlux(FLX, MFLX, GEOM, NMESH, 1, NG, NG)
END SUBROUTINE

SUBROUTINE ALLOC_MFLX(MFLX, NMESH, NG)
USE ALLOCS
IMPLICIT NONE
TYPE(MODALFLX_TYPE) :: MFLX(NMESH)
INTEGER :: NMESH, NG
INTEGER :: I

DO I = 1, NMESH
  CALL DMALLOC0(MFLX(I)%PHI, 0, 4, 1, 2, 1,NG)
  CALL DMALLOC0(MFLX(I)%QH, 0, 4, 1, 2)
  CALL DMALLOC(MFLX(I)%A, 2, NG)
  CALL DMALLOC(MFLX(I)%B, 2, NG)
  CALL DMALLOC0(MFLX(I)%PSOL ,0, 4, 1, 2, 1, NG)
ENDDO
ENDSUBROUTINE


SUBROUTINE AllocSP3Solver(AxFlx, ixybeg, ixyend, izbeg, izend, ng0, lUse)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type, PE_TYPE
IMPLICIT NONE
TYPE(AxFlx_Type), POINTER :: AxFlx(:, :)
TYPE(CoreInfo_Type) :: Core
TYPE(PE_TYPE) :: PE
INTEGER :: ng0, ixybeg, ixyend, izbeg, izend
INTEGER :: nxy, myzbf, myzef
INTEGER :: ixy, iz
LOGICAL :: lUse

!nxy = Core%nxy
!myzbf = PE%myzbf; myzef = PE%myzef
DO ixy =  ixybeg, ixyend
  DO iz = izbeg, izend
    CALL AllocSP3AxFlxType(AxFlx(iz, ixy), ng0, lUse)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE AllocSP3AxFlxType(AxFlx, ng0, lUse)
USE PARAM
USE ALLOCS
IMPLICIT NONE
TYPE(AxFlx_Type) :: AxFlx
INTEGER :: ng0
LOGICAL :: luse

IF(lUse) THEN
  CALL DMALLOC0(AXFLX%PHI, 0, 4, 1, 2, 1, NG0); CALL DMALLOC0(AXFLX%PSI, 0, 4)
  CALL DMALLOC0(AXFLX%TLKG, 0, 2, 1, 2, 1, NG0); CALL DMALLOC0(AXFLX%JOUT, 1, 2, 1, 2, 1, NG0)

  CALL DMALLOC0(AXFLX%S,1,4,1,NG0);CALL DMALLOC0(AXFLX%KSQ,1,4,1,NG0);
  CALL DMALLOC0(AXFLX%QT,1,4,1,NG0);
  CALL Dmalloc0(AxFlx%LkgShape, 0, 4, 1, ng0)
ENDIF

CALL DMALLOC0(AXFLX%DHAT, 1, 2, 1, NG0); CALL DMALLOC0(AXFLX%DTIL, 1, 2, 1, NG0)
CALL DMALLOC0(AXFLX%PDHAT, 1, 2, 1, NG0);

END SUBROUTINE

END MODULE