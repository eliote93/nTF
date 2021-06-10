#include <defines.h>
#ifdef __GAMMA_TRANSPORT
!--- JSU Edit
MODULE GammaTYPEDEF    
USE CSRMATRIX,        ONLY : CSR_DOUBLE
IMPLICIT NONE


TYPE PSMAT_TYPE
  SEQUENCE
  INTEGER :: IB, IE, ioutsb, ioutse
  REAL, POINTER :: FROM(:)
  REAL :: self
END TYPE


TYPE GamMacXS_TYPE
  INTEGER :: id = 0               ! FOR TEMPORARY DATA IN CMFD CUSPING PIN CALCULATION
  INTEGER :: ng = 0, ngg = 0, niso = 0
!  LOGICAL :: lprodalloc = .FALSE.
!  REAL, POINTER :: GMacProdMTot(:,:)                   ! total production matrix
!  REAL, POINTER :: GMacProdMFis(:,:)                   ! fission production matrix
!  REAL, POINTER :: GMacProdMRad(:,:)                   ! radioactuve capture production matrix
!  REAL, POINTER :: GMacProdMIne(:,:)                   ! inelastic scattering production matrix
  ! PHOTOATOMIC REACTION CROSS SECTION           
  LOGICAL :: lalloc = .FALSE.
  REAL, POINTER :: MacKERMA(:)                  ! Total KERMA of photon reaction
  REAL, POINTER :: XsMacT(:)                    ! Total XS
  REAL, POINTER :: XsMacTR(:)                   ! Transport
  REAL, POINTER :: XsMacA(:)                    ! Absorption
  REAL, POINTER :: XsMacS(:)                    ! Total Scattering
  REAL, POINTER :: XsMacSTR(:)                  ! Total Scattering
  REAL, POINTER :: XsMacSM(:,:)                 ! Scattering Matrix
  LOGICAL :: lallocsm = .FALSE.
  REAL, POINTER :: MacGSM1(:,:)               ! P1Scattering Matrices
  REAL, POINTER :: MacGSM2(:,:)               ! P2Scattering Matrices
  REAL, POINTER :: MacGSM3(:,:)               ! P3Scattering Matrices
  ! Element-wise scattering matrix..
  LOGICAL :: lIsoSMAlloc = .FALSE.
  REAL, POINTER :: IsoSM(:, :, :)             ! Element-wise scattering matrices..
  ! ELEMENT-WISE CROSS SECTION
  LOGICAL :: lIsoAlloc = .FALSE.
  REAL, POINTER :: IsoKERMA(:,:)                  ! Total KERMA of photon reaction
  REAL, POINTER :: IsoXsMacT(:,:)                    ! Total XS
  REAL, POINTER :: IsoXsMacTR(:,:)                   ! Transport
  REAL, POINTER :: IsoXsMacA(:,:)                    ! Absorption
  REAL, POINTER :: IsoXsMacS(:,:)                    ! Total Scattering
  REAL, POINTER :: IsoXsMacSTR(:,:)                  ! Total Scattering
! PHOTON PRODUCTION DATA
  LOGICAL :: lIsoProdAlloc = .FALSE.
  LOGICAL, POINTER :: lexist(:, :)
  REAL, POINTER :: IsoGProdFis(:,:,:)
  REAL, POINTER :: IsoGProdRad(:,:,:)
  REAL, POINTER :: IsoGProdInel(:,:,:)
  REAL, POINTER :: IsoGProdNnel(:,:,:)
  LOGICAL :: lProdAlloc = .FALSE.
  REAL, POINTER :: GProdTot(:,:)
! FOR EXPLICIT KAPPA CALCULATION                                         |-- JSU EDIT 2017.09.14. |
!  LOGICAL :: lInelProd = .FALSE.
!  LOGICAL :: lKappaAlloc = .FALSE.
!  REAL, POINTER :: FisLocal(:)         ! LOCAL DEPOSIT FISSION ENERGY
!  REAL, POINTER :: QN2N(:), QN3N(:)    ! N2N, N3N REACTION ENERGY
!  REAL, POINTER :: InelLocal(:)        ! 
!  REAL, POINTER :: LocalQ(:)      ! INELASTIC SCATTERING ENERGY LOSS MATRIX
!  REAL, POINTER :: GProdInel(:, :)     ! GAMMA ENERGY FROM INELASTIC SCATTERING MATRIX
END TYPE

TYPE GPINXS_TYPE
  SEQUENCE
  TYPE(PSMAT_TYPE), POINTER, DIMENSION(:) :: XSS  ! SCATTERING MATRIX
  REAL, POINTER, DIMENSION(:,:) :: DTIL, DHAT, PDHAT, DTIL2
  REAL, POINTER, DIMENSION(:,:) :: atil, ahat
  REAL, POINTER, DIMENSION(:,:) :: AxDtil, AxDhat
  REAL, POINTER, DIMENSION(:) :: XSD, XSD2, XST, XSTR, XSR, XSA, XSDA
  REAL, POINTER, DIMENSION(:) :: NPHI, GPHI, FAC, XSRD
  REAL, POINTER, DIMENSION(:,:,:) :: Dcp1_DHAT
  REAL :: FuelTemp, ModTemp, PinTemp
  !PRODUCTION MATRIX
  REAL, POINTER, DIMENSION(:,:) :: XSP ! PRODUCTION MATRIX
  !Kinetic Parameter
  REAL :: BETAT, BETAP, OMEGA 
  REAL :: xstr1g, xstr1g0                  !Absoprtion XS 1group
  !REAL :: BETAT,OMEGA,BETAP
  REAL,POINTER,DIMENSION(:) :: BETA, CHIP, VELO,RVDELT
END TYPE

TYPE PINHEATXS_TYPE
  SEQUENCE
  REAL, POINTER, DIMENSION(:) :: KERMA_T, KERMA_F, PhotonGen, kFis ! Neutron induced reaction
  REAL, POINTER, DIMENSION(:) :: KERMA_P  
  REAL, POINTER, DIMENSION(:) :: NPHI, GPHI
  REAL :: FuelTemp, ModTemp, PinTemp                     ! Photoatomic Reaction
END TYPE

TYPE GammaTrackingDat_Type
  REAL, POINTER :: phia1g(:, :, :), SrcAng1g(:, :, :), xst1g(:)
  REAL, POINTER :: phisnm(:, :), srcnm(:, :), xstnm(:, :)
  REAL, POINTER :: PhiAngIn1g(:, :), Jout1g(:, :, :)
  REAL, POINTER :: PhiAngInnm(:, :, :), Joutnm(:, :, :, :)
  REAL, POINTER :: Expa(:, :), Expb(:, :), wtang(:, :), wtsurf(:, :, :)
  INTEGER, POINTER :: AziMap(:, :)
END TYPE

TYPE GammaCMFD_Type
  TYPE(CSR_DOUBLE), ALLOCATABLE :: M(:), ILU(:)
  REAL, POINTER :: Scat(:, :, :), Prod(:, :, :)
  REAL, POINTER :: Phic(:, :, :), gPhic(:, :, :), gProd(:, :), gSrc(:)
END TYPE

END MODULE
  
    
MODULE GammaCore_mod
USE GammaTYPEDEF,  ONLY : GPINXS_TYPE, PINHEATXS_TYPE

! PIN XS
TYPE(GPINXS_TYPE), POINTER ::GPINXS(:,:)

REAL, POINTER :: gphis(:, :, :)         ! Gamma Scalar Flux
REAL, POINTER :: gphim(:, :, :, :)      ! Gamma Angular Flux Moment
REAL, POINTER :: gphic(:, :, :)
REAL, POINTER :: gPhiAngIn(:, :, :, :)  ! Gamma Incoming Angular Flux
REAL, POINTER :: gJout(:, :, :, :, :)   ! Gamma Pin Current
REAL, POINTER :: gAxSrc(:, :, :)
REAL, POINTER :: gAxPXS(:, :, :)
REAL, POINTER :: gPower(:, :)           ! Gamma Power
REAL, POINTER :: LocalNPower(:, :)        ! NEUTRON LOCAL DEPOSITED ENERGY   -- JSU EDIT 2017.09.15.
REAL, POINTER :: LocalNPower_KERMA(:, :)  ! NEUTRON LOCAL DEPOSITED ENERGY Calculated with HEATR KERMA  -- JSU EDIT 2019.08.12.
REAL, POINTER :: GPowerGen(:, :)          ! Photon energy generation in the region from neutron induced reaction   -- JSU EDIT 2019.08.12.

TYPE(PINHEATXS_TYPE), POINTER :: PinHeatXS(:,:)

END MODULE
#endif