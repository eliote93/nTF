MODULE DeplType
!     --------------------------------------------------------------------------------
!     @ NISODEP                (I)    THE NUMBER OF ISOTOPES IN LIBRARY              +
!     @ NATOMDEP               (I)    MAXIMUM ATOMIC NUMBER                          +
!     @ NFISYLD                (I)    THE NUMBER OF ACTINIDE W/ F.P.Y.               +
!     @ STATEKIND(:)%          (T)    TYPE FOR BURNABLE ISOTOPES                     +
!     @        %IGRP           (I)    GROUP NUMBER (1:ACTIVATION,2:ACTINIDE,3:F.P.)  +
!     @        %IMAT           (I)    MATRIX ID.                                     +
!     @        %NTO1           (I)    NUMBER OF DECAY CHILDREN                       +
!     @        %NTO2           (I)    NUMBER OF REACTION CHILDREN                    +
!     @        %NFR3           (I)    NUMBER OF FISSION SOURCES                      +
!     @        %ITO1(I,J)      (I)    DECAY CHILDS                                   +
!                                      (1,:) DECAY BRANCH NUMBER                     +
!                                      (2,:) CHILD ISOTOPE ID (MATID)                +
!                                      (3,:) POSITION IN CHILD MATRIX ELEMENTS       +
!     @        %ITO2(I,J)      (I)    REACTION CHILDS                                +
!                                      (1,:) REACTION BRANCH NUMBER                  +
!                                      (2,:) CHILD ISOTOPE ID (MATID)                +
!                                      (3,:) POSITION IN CHILD MATRIX ELEMENTS       +
!     @        %IFR3(I,J)      (I)    FISSION SOURCES                                +
!                                      (1,:) FISSIONABLE ISOTOPE ID IN FISYLD        +
!                                      (2,:) POSITION IN MY MATRIX ELEMENTS          +
!     @        %RAMBDA         (R)    DECAY CONSTANT                                 +
!     @        %FRAC(I)        (R)    FRACTION FOR EACH DECAY BRANCH(I=1,8)          +
!                                      1 : BETA(-) DECAY TO METASTABLE               +
!                                      2 : POSITRON EMISSION OR ELECTRON CAPTURE     +
!                                      3 : POSITRON EMISSION OR ELECTRON CAPTURE TO  +
!                                          METASTABLE                                +
!                                      4 : ALPHA DECAY FRACTION                      +
!                                      5 : METASTABLE TO STABLE                      +
!                                      6 : SPONTANEOUS FISSION                       +
!                                      7 : BETA(-)+NEUTRON DECAY                     +
!                                      8 : BETA(-) DECAY TO STABLE                   +
!     @        %XS(I)          (R)    CROSS SECTIONS (I=0,6)                         +
!                                      1 : (N,R)  (Z,A,0)->(Z,A+1,0)                 +
!                                      2 : (N,2N) (Z,A,0)->(Z,A-1,0)                 +
!                                      3 : (N,A)  (Z,A,0)->(Z-2,A-3,0)               +
!                                      4 : (N,P)  (Z,A,0)->(Z-1,A,0)                 +
!                                      5 : (N,R') (Z,A,0)->(Z,A+1,1)                 +
!                                      6 : (N,2N') (Z,A,0)->(Z,A-1,1)                +
!                                      3 : (N,3N) (Z,A,0)->(Z,A-2,0) FOR ACTINIDE (J=3)
!                                      4 : FISSION FOR ACTINIDE (J=4)                +
!     @        %Y(I)           (R)    FISSION PRODUCT YIELDS                         +
!     @ ISOTOPEKIND(:)%        (T)    TYPE FOR BURNABLE ISOTOPES                     +
!     @        %NSTATE         (I)    NUMBER OF STATES                               +
!     @        %STAT(I)        (T)    TYPE STATEKIND                                 +
!     @ ATOMKIND(:)%           (T)    TYPE FOR BURNABLE ATOM                         +
!     @        %IB             (I)    ISOTOPE BEGINNING NUMBER                       +
!     @        %IE             (I)    ISOTOPE ENDING NUMBER                          +
!     @        %A(I)           (T)    TYPE ISOTOPEKIND                               +
!     @ MATRIXDEPL(:)%         (T)    DEPLETION MATRIX                               +
!     @        %N              (I)    NUMBER OF ELEMENT                              +
!     @        %IFR(I)         (I)    INDICATOR IN ISOTOPE VECTOR                    +
!     @        %ELE(I)         (R)    ELEMENT IN IFR(I)                              +
!     @ IFISYD(:,:)            (I)    INFORMATION FOR FISSION YIELD ISOTOPE          +
!                                      (0,:) - MATRIX ID                             +
!                                      (1,:) - ATOMIC NUMBER                         +
!                                      (2,:) - ATOMIC WEIGHT                         +
!     @ Z0(:),Z1(:)       (ATOMKIND)  ATOMS IN LIBRARY, Z0 FOR KEEP, Z1 FOR USE      +
!     @ VISO(:)                (R)    ISOTOPE NUMBER DENSITY                         +
!     @ DMAT(:)          (MATRIXDEPL) DEPLETION MATRIX                               +
!     @ IM2Z(:)                (I)    MAPPING DATA FROM MATRIX ID TO ATOMIC NUMBER   +
!     @ IM2A(:)                (I)    MAPPING DATA FROM MATRIX ID TO ISOTOPE WEIGHT  +
!     @ IM2ST(:)               (I)    MAPPING DATA FROM MATRIX ID TO ISOTOPE STATE   +
!     --------------------------------------------------------------------------------
USE MatExp_Mod, ONLY : Mat_Type,  MatExp_Type
IMPLICIT NONE

INTEGER, PARAMETER, PRIVATE :: nMAXSTATE=100

TYPE STATEKIND
  INTEGER :: IGRP, IMAT, NTO1, NTO2, NFR3, idum
  INTEGER, POINTER :: ITO1(:,:), ITO2(:,:), IFR3(:,:)
  REAL :: RAMBDA, FRAC(8)
  REAL :: XS(0:6),kappa
  REAL, POINTER :: Y(:)
END TYPE

TYPE ISOTOPEKIND
  INTEGER :: NSTATE, idum
  TYPE(STATEKIND),POINTER :: STAT(:)
END TYPE

TYPE ATOMKIND
  INTEGER :: IB,IE
  TYPE(ISOTOPEKIND), POINTER :: A(:)  ! ATOMIC WEIGHT
END TYPE

TYPE FisYield_TYPE
  INTEGER :: MatID
  INTEGER :: AtomNum
  INTEGER :: AtomWeight
  INTEGER :: Stat
END TYPE

TYPE DeplLib_Type
  TYPE(ATOMKIND), POINTER :: AtomLib0(:), AtomLib1(:)
  TYPE(FisYield_TYPE), POINTER :: FisYield(:)
  INTEGER :: NISODEP, NFISYLD, NATOMDEP
  INTEGER, POINTER :: MapMatId2ANum(:)   !IM2Z
  INTEGER, POINTER :: MapMatId2IsoWt(:)  !IM2A
  INTEGER, POINTER :: MapMatId2State(:)  !IM2ST
  CHARACTER(80) :: FileName
END TYPE

TYPE CoreState_Type
  INTEGER :: nBurnupStep=0
  !LOGICAL :: lBoronSearch = .FALSE.
  LOGICAL :: LinStateChg = .FALSE.
  LOGICAL :: StepStateChg = .TRUE.
  LOGICAL :: lStateChg(0:nMAXSTATE) = .FALSE.
  LOGICAL :: lBoronSearch(0:nMAXSTATE) = .FALSE.
  REAL :: RelPow(0:nMAXSTATE)=0
  REAL :: FlowRate(0:nMAXSTATE)=0
  REAL :: BoronPPM(0:nMAXSTATE)=0
  REAL :: Target_Keff(0:nMAXSTATE) = 1.0_8
  REAL :: T_efpd(0:nMAXSTATE)=0
  REAL :: T_mwdkghm(0:nMAXSTATE)=0
END TYPE

TYPE DeplCntl_Type
  INTEGER :: NowStep = 0
  INTEGER :: BurnUpType = 1
  INTEGER :: nBurnUpStep
  INTEGER :: PC_OPT = 1
  INTEGER :: SOLVER = 1
  LOGICAL :: lInitDepl = .FALSE.
  LOGICAL :: lDeplFile = .FALSE.
  LOGICAL :: lB10Depl = .FALSE.    !Boron-10 Depletion
  INTEGER :: B10DeplMod = 0       !Boron-10 Depletion   0 : TurnOff 1: Online  2 : Post Boron
  REAL :: vratio = 0.05
  INTEGER :: nSubStep = 8         !SubStep for Gd Depletion

  LOGICAL :: lXeDyn = .FALSE.
  LOGICAL :: lTrXe = .TRUE.
  LOGICAL :: lEqXe = .FALSE.

  REAL, POINTER :: T_efpd(:)
  REAL, POINTER :: T_mwdkghm(:)
  REAL :: Hm_Mass0_ton = 0.  !Initially Loaded Heavy Metal in Ton Unit
  REAL :: Hm_Mass0_kg = 0.
  REAL :: PowerCore = 0.
  REAL :: Tsec = 0.        !Burn Up Step Time size
  REAL :: PTsec = 0.
  REAL :: NowBurnUp(2) = 0 !1: EFPD, 2 : mwd/kghm
  !LOGICAL :: lAutoOrder = .FALSE. !  which was originally .TRUE.   !BYS test 15/12/31
  LOGICAL :: lAutoOrder = .TRUE.
  LOGICAL :: lPredict = .TRUE.

  LOGICAL :: lCoreFollow = .FALSE.
  TYPE(CoreState_Type) :: CoreState
END TYPE



TYPE DeplVars_Type
  INTEGER :: nIsoDepl, nIsoXs, tid
  TYPE(Mat_TYPE), POINTER :: DMat
  REAL, POINTER :: IsoNum(:), BurnUpXs(:, :)
  REAL :: phi1g
  INTEGER, POINTER :: MapXs2Dep(:), MapDep2XS(:)

  REAL :: GdXsFtn(0:2, 64152:64160)
ENDTYPE

!TYPE MatExp_Type
!
!  INTEGER :: nisodepl
!  TYPE(Mat_TYPE), POINTER :: Mat
!  LOGICAL :: lAllocVec = .FALSE.
!  REAL, POINTER :: Viso0(:), Viso(:)
!ENDTYPE

TYPE DeplXs_Type
  LOGICAL :: lAlloc = .FALSE.
  INTEGER :: ndim
  REAL, POINTER :: xsa(:), xsf(:)      !Absorption, Fission
  REAL, POINTER :: xsn2n(:), xsn3n(:)  !n2n, n3n reaction
  REAL, POINTER :: AvgPhi(:)
  REAL :: Phi1g
  INTEGER :: tid
END TYPE

CONTAINS

SUBROUTINE AllocDeplXS(DeplXs, niso, ng)
USE Allocs
IMPLICIT NONE
TYPE(DeplXs_Type) :: DeplXs
INTEGER :: niso, ng
DeplXs%lAlloc = .TRUE.
DeplXs%ndim = niso
CALL Dmalloc(DeplXs%xsa, niso)
CALL Dmalloc(DeplXs%xsf, niso)
CALL Dmalloc(DeplXs%xsn2n, niso)
CALL Dmalloc(DeplXs%xsn3n, niso)
CALL Dmalloc(DeplXs%AvgPhi, ng)
END SUBROUTINE

SUBROUTINE AllocMatExpVec(MatExp, ndim)
USE Allocs
IMPLICIT NONE
TYPE(MatExp_Type) :: MatExp
INTEGER :: ndim
CALL Dmalloc(MatExp%Viso0, ndim)
CALL Dmalloc(MatExp%Viso, ndim)
MatExp%lAllocVec = .TRUE.
END SUBROUTINE

END MODULE
