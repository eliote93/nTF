#include <defines.h>
#ifdef __GAMMA_TRANSPORT
!--- JSU Edit
MODULE GammaTYPEDEF    
USE CSRMATRIX,        ONLY : CSR_DOUBLE
IMPLICIT NONE

TYPE GAMMAGROUPINFO_TYPE
  SEQUENCE 
  LOGICAL :: lUpScat, ldum
  INTEGER :: ng, ngg, nofg, norg, ntiso
  INTEGER :: nFxSub, nCatMg
  INTEGER :: UpscatRange(2)
  INTEGER, POINTER :: InScatRange(:, :), OutScatRange(:, :)
  INTEGER :: nGC                 !Number of Condensed Group
  INTEGER :: GCStruct(2,300)     !Coarse Group => Fine Group
  INTEGER :: InvGCStruct(300)    !FineGroup => Coarse Group
  INTEGER :: nprec = 6
  ! AVERAGE (GEOMATRICAL AVERAGE) OF ENERGY IN THE GROUP                  |-- JSU EDIT 2017.09.13. |
  REAL, POINTER :: GamAvgE(:), NeuAvgE(:)   ! in Joule (1.602*1.e-19eV)
END TYPE

TYPE PMAT_TYPE
  SEQUENCE
  INTEGER :: IB, IE
  REAL, POINTER :: FROM(:)
END TYPE

TYPE PSMAT_TYPE
  SEQUENCE
  INTEGER :: IB, IE, ioutsb, ioutse
  REAL, POINTER :: FROM(:)
  REAL :: WithInGroupScat
END TYPE

TYPE GAMLIBDATA
  SEQUENCE
  INTEGER :: nid
  CHARACTER*8  :: aid
  REAL :: aw
  REAL, POINTER, DIMENSION(:) :: TEMP
  INTEGER :: ityp         !  0:ELEMENT(REACTION), 1:ISOTOPE(PRODUCTION), 2:NATURAL(REACTION+PRODUCTION)
  ! PHOTON PRODUCTION DATA (NEUTRON INDUCED)      !--- JSU EDIT 20170720       FOR ITYP = 1
  !        REACTIONS(IMT) 0:TOTAL, 1:FISSION, 2: RADIOACTIVE CAPTURE, 3:INELEATIC
  INTEGER :: ifis = 0, iradcap = 0, iinel = 0     ! INTEGER VARIABLE INDICATES PHOTON SOURCE REACTION
  INTEGER :: ntemp
  TYPE(PMAT_TYPE), POINTER, DIMENSION(:, :, :) :: GPM  ! GAMMA PRODUCTION MATRIX, (GAMMA GROUP, TEMPERATURE, REACTION)
  LOGICAL :: LPHOTON(0:3) = .FALSE.               ! LOGICAL VARIABLE INDICATES PHOTON SOURCE REACTION
  INTEGER :: PPGRUP(0:3) = 0, PPGRLOW(0:3) = 0    ! GROUP UPPER AND LOWER BOUND EXISTING PHOTON SOURCE CORRESPONDING REACTION
  
  ! NEUTRON GROUP BOUNDARY CONTRIBUTING PHOTON GENERATION, (NEUTRON GROUP, TEMPERATURE, REACTION)  <-- JSU EDIT 2017.09.14.
  INTEGER, POINTER, DIMENSION(:, :, :) :: GPMOUTB, GPMOUTE
  
  ! PHOTON REACTION DATA                          !--- JSU EDIT 20170720       FOR ITYP ~= 1 (0, 2)
  !        ONLY FOR ELEMENT **000 (PHOTO-ATOMIC REACTION)
  REAL, POINTER, DIMENSION(:) :: GSIGA, GSIGTR, GSIGS, GSIGSTR, KERMA  ! TRANSPORT CORRECTED XS FOR PHOTON TRANSPORT
  TYPE(PSMAT_TYPE), POINTER, DIMENSION(:) :: GSM, GSM1, GSM2, GSM3   ! SCATTERING MATRIX UPTO 3rd ORDERS
  REAL, POINTER, DIMENSION(:) :: GSIGSP1, GSIGSP2, GSIGSP3
  ! KAPPA DATA                                                            |-- JSU EDIT 2017.09.12. |
  REAL :: CAPKAPPA = 0. ! RADIO
  REAL :: CapQ, n2nQ, n3nQ
  LOGICAL :: lfis = .FALSE.
  REAL :: EXKAPPA(6) = 0., EXKAPPA0(6) = 0.  ! 1: FISSION FRAGMENTS, 2: PROMPT NEUTRON, 3: DELAYED NEUTRON, 4: PROMPT GAMMA 5: DELAYED GAMMA, 6: DELAYED BETA

END TYPE

TYPE GamMacXS_TYPE
  SEQUENCE
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
  ! ELEMENT-WISE CROSS SECTION
  LOGICAL :: lIsoAlloc = .FALSE.
  REAL, POINTER :: IsoKERMA(:,:)                  ! Total KERMA of photon reaction
  REAL, POINTER :: IsoXsMacT(:,:)                    ! Total XS
  REAL, POINTER :: IsoXsMacTR(:,:)                   ! Transport
  REAL, POINTER :: IsoXsMacA(:,:)                    ! Absorption
  REAL, POINTER :: IsoXsMacS(:,:)                    ! Total Scattering
  REAL, POINTER :: IsoXsMacSTR(:,:)                  ! Total Scattering
! PHOTON PRODUCTION DATA
  LOGICAL :: lProdAlloc = .FALSE.
  LOGICAL, POINTER :: lfis(:), lrad(:), linel(:)
  INTEGER, POINTER :: ifisb(:), ifise(:)  ! Lower and Upper bound of gamma group data exist 
  REAL, POINTER :: IsoGProdFis(:,:,:)
  INTEGER, POINTER :: iradb(:), irade(:)
  REAL, POINTER :: IsoGProdRad(:,:,:)
  INTEGER, POINTER :: iinelb(:), iinele(:)  
  REAL, POINTER :: IsoGProdInel(:,:,:)

!  REAL, POINTER :: GProdFis(:,:)
!  REAL, POINTER :: GProdRad(:,:)
  REAL, POINTER :: GProdTot(:,:)
! FOR EXPLICIT KAPPA CALCULATION                                         |-- JSU EDIT 2017.09.14. |
  LOGICAL :: lInelProd = .FALSE.
  LOGICAL :: lKappaAlloc = .FALSE.
  REAL, POINTER :: FisLocal(:)         ! LOCAL DEPOSIT FISSION ENERGY
  REAL, POINTER :: QN2N(:), QN3N(:)    ! N2N, N3N REACTION ENERGY
  REAL, POINTER :: InelLocal(:)        ! 
  REAL, POINTER :: LocalQ(:)      ! INELASTIC SCATTERING ENERGY LOSS MATRIX
!  REAL, POINTER :: GProdInel(:, :)     ! GAMMA ENERGY FROM INELASTIC SCATTERING MATRIX
END TYPE

TYPE GPINXS_TYPE
  SEQUENCE
  TYPE(PSMAT_TYPE), POINTER, DIMENSION(:) :: XSS  ! SCATTERING MATRIX
  REAL, POINTER, DIMENSION(:,:) :: DTIL, DHAT, PDHAT, DTIL2
  REAL, POINTER, DIMENSION(:,:) :: atil, ahat
  REAL, POINTER, DIMENSION(:,:) :: AxDtil, AxDhat
  REAL, POINTER, DIMENSION(:) :: XSD, XSD2, XST, XSTR, XSR, KERMA, XSA, XSDA
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
  REAL, ALLOCATABLE :: Scat(:, :, :), Prod(:, :, :)
  REAL, ALLOCATABLE :: Phic(:, :, :), gPhic(:, :, :), gProd(:, :), gSrc(:)
END TYPE

END MODULE
  
MODULE GammaLibdata  
USE GammaTYPEDEF,    ONLY : GAMLIBDATA
! IMT INDEX
INTEGER, PARAMETER :: imtTOT = 0, imtFIS = 1, imtRAD = 2, imtINEL = 3

! PHOTON LIBRARY INFORMATION
CHARACTER*255 :: GFILE                                ! PHOTON LIBRARY FILE NAME
REAL(4), POINTER, DIMENSION(:) :: NeuENB, NeutronU    ! NEUTRON ENERGY GROUP BOUNDARY (UPPER)
REAL(4), POINTER, DIMENSION(:) :: GamENB, GammaU      ! GAMMA ENERGY GROUP BOUNDARY (UPPER)
INTEGER :: NOGGAM, NOGGGAM                            ! A # OF NEUTRON AND PHOTON ENERGY GROUP
! JUST FOR THE CASE IT IS NEEDED...
INTEGER :: nofgGAM,norgGAM,notgGAM,neltGAM,nfisGAM,nradGAM,ninelGAM,nelmGAM
INTEGER :: igresb, igrese
! PHOTON LIBRARY DATA
TYPE(GAMLIBDATA), POINTER, dimension(:) :: Gldiso

! xs library (helios) data variables
INTEGER, POINTER, DIMENSION(:) :: nuclidgam,idfisGAM,idradGAM,idinelGAM,idelmGAM,elementGAM,iso2elm
! mapping of nuclide ID to helios library order
INTEGER, PARAMETER :: nxtempmap=5000
INTEGER, POINTER :: itempmapGAM(:, :)
INTEGER :: mapnuclGAM(100000), mapnuclELM(100000)

INTEGER :: TEMPMax

END MODULE
    
MODULE GammaCore_mod
USE GammaTYPEDEF,  ONLY : GAMMAGROUPINFO_TYPE, GPINXS_TYPE

! PIN XS
TYPE(GPINXS_TYPE), POINTER ::GPINXS(:,:)

! GAMMA GROUP INFO
TYPE(GAMMAGROUPINFO_TYPE) :: GamGroupInfo

REAL, POINTER :: gphis(:, :, :)         ! Gamma Scalar Flux
REAL, POINTER :: gphim(:, :, :, :)      ! Gamma Angular Flux Moment
REAL, POINTER :: gPhiAngIn(:, :, :, :)  ! Gamma Incoming Angular Flux
REAL, POINTER :: gJout(:, :, :, :, :)   ! Gamma Pin Current
REAL, POINTER :: gPower(:, :)           ! Gamma Power
REAL, POINTER :: LocalNPower(:, :)      ! NEUTRON LOCAL DEPOSITED ENERGY   -- JSU EDIT 2017.09.15.
END MODULE
#endif