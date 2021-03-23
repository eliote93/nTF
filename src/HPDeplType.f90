MODULE HPDeplType
  USE CSRMATRIX
  IMPLICIT NONE
  INTEGER, PARAMETER :: RctIdCAP = 1, RctIdN2N = 2, RctIdNA = 3, RctIdNP = 4, RctIdCAPm = 5, RctIdN2Nm = 6
  INTEGER, PARAMETER :: RctIdN3N = 3, RctIdFis = 4
  INTEGER, PARAMETER :: KrySolTyp = 1, CRAMSolTyp = 2
  TYPE DeplFxr_Type
    INTEGER :: iDeplFxr, NisoEig
    LOGICAL :: lDepl, lFuel, lGd ! (T, F, F) => Structural Material
    REAL(8), POINTER :: pnum_sseig(:), pnum_depl(:)
    REAL(8), POINTER :: pnum_pre(:), pnum_cor(:)
    REAL(8), POINTER :: xs1g(:,:) ! (NofRct, NofIso)
    REAL(8), POINTER :: Kappa(:) ! (NofIso)
    INTEGER, POINTER :: IdIsoEig(:)
    REAL(8) :: burnup, burnup0
    REAL(8) :: Tmpt, Vol, Hmkg0, NormFlux1g
  END TYPE
  TYPE GdFxrInfo_Type ! Child, Parent : DeplFxr
    TYPE(DeplFxr_Type), POINTER :: aFxr ! Size of just one
    REAL(8), POINTER :: GdRR(:,:) ! size : (-1:1, NofGdIso)
    REAL(8), POINTER :: f_pc(:) ! post-correction factor, (NofGdIso)
    REAL(8), POINTER :: Gdpnum(:,:,:) ! PC Gd ND (2, -1:1, NofGdIso)
    REAL(8), POINTER :: Gd155(:)
    REAL(8), POINTER :: c_qd(:,:) ! Gd Quadratic Coefficients (3, NofGdIso)
  END TYPE
  TYPE DeplFxrBundle_Type
    INTEGER :: nfxr, nTrueDepl, nTrueGd, nIsoGd, i155
    INTEGER, POINTER :: IdTrueDepl(:), IdTrueGd(:)
    INTEGER, POINTER :: IdIsoGd(:)
    TYPE(DeplFxr_Type), POINTER :: FxrBundle(:)
    TYPE(GdFxrInfo_Type), POINTER :: GdFxrBundle(:)
    REAL(8) :: delT, Hmkg0, Power ! (sec), (kg), (MW)
  END TYPE
  TYPE DeplLib_Type
    INTEGER :: NofIso, NofFis, NofYld
    INTEGER :: NofRct, NofDec
    INTEGER :: FisRctId
    INTEGER, POINTER :: Idiso(:)
    LOGICAL, POINTER :: lActinide(:)
    REAL(8), POINTER :: RctXs(:,:)       ! (NofRct, NofIso)
    REAL(8), POINTER :: FisFrac(:,:)     ! (NofFis, NofYld)
    REAL(8), POINTER :: DecFrac(:,:)     ! (NofDec, NofIso)
    REAL(8), POINTER :: Kappa(:)         ! (NofIso)
    INTEGER, ALLOCATABLE :: delIdRct(:), delIdRctAct(:), delIdDec(:)   ! Change of Index from reaction or decay, (NofRct), (NofDec)
    INTEGER, POINTER :: SubYldRct(:), SubYldRctAct(:), SubYldDec(:)
    INTEGER, POINTER :: IdAftRct(:,:), IdAftDec(:,:) ! (NofRct, NofIso), (NofDec, NofIso)
    INTEGER, POINTER :: IdFis(:), IdYld(:) ! (NofFis), (NofYld)
    INTEGER, POINTER :: IzSYR(:,:), IzSYRAct(:,:), IzSYD(:,:), IzRct(:,:), IzDec(:,:), IzFY(:,:)
    INTEGER, ALLOCATABLE :: YldMapColIdx(:), YldMapRowPtr(:) ! Diagonal filled
  END TYPE
  TYPE DeplSysBundle_Type
    INTEGER :: NofIso, Nsys, ifxrbeg
    TYPE(CSR_DOUBLE), POINTER :: DeplMats(:)
    ! ----------- For cuDEPL Solvers --------------- !
!    COMPLEX(8), POINTER :: Diag(:,:) ! (Nsys,NR)
    REAL(8), POINTER :: Diag(:,:) ! (Nsys,NR)
    REAL(8), POINTER :: OffDiag(:,:) ! (Nsys,NNZ-NR)
    ! ---------------------------------------------- !
    REAL(8), POINTER :: pnums(:) ! (NofIso*Nsys)
    REAL(8), POINTER :: pnums1(:) ! (NofIso*Nsys)
  END TYPE
END MODULE
