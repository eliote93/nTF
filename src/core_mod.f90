MODULE Core_mod
  
USE PARAM,   ONLY : ONE
USE TYPEDEF, ONLY : Fxrinfo_type, PinXs_Type, GroupInfo_Type, CMFDLS_TYPE, AXFLX_TYPE, CMInfo_Type, FMInfo_Type, THInfo_Type

INTEGER :: nCoreFsr, nCoreFxr, nPhiAngSv

REAL :: eigv      = ONE
REAL :: peigv     = ONE
REAL :: eigv_adj  = ONE
REAL :: peigv_adj = ONE

REAL :: xkconv ! Converged K
REAL :: xkcrit ! Critical K

! Basic 
REAL, POINTER, SAVE, DIMENSION(:)           :: wmoc
REAL, POINTER, SAVE, DIMENSION(:,:)         :: psi, psid    ! FSR FIssion Scource
REAL, POINTER, SAVE, DIMENSION(:,:)         :: power
REAL, POINTER, SAVE, DIMENSION(:,:)         :: psic, psicd  ! Cell Fission Source
REAL, POINTER, SAVE, DIMENSION(:,:,:)       :: phis         ! Scalar Flux
REAL, POINTER, SAVE, DIMENSION(:,:,:)       :: AxSrc, AxPxs ! Axial Source, Psuedo Absorption XS
REAL, POINTER, SAVE, DIMENSION(:,:,:,:)     :: PhiAngIn     ! Angular Flux at the Boundary
REAL, POINTER, SAVE, DIMENSION(:,:,:,:)     :: phim
REAL, POINTER, SAVE, DIMENSION(:,:,:,:,:)   :: RadJout      ! Moc Jout (idir, ibd, ixy, iz, ig)
REAL, POINTER, SAVE, DIMENSION(:,:,:,:,:,:) :: phia         ! Angular Flux (idir, ipol, iazi, ifsr, ig, iz)

! Linear Source
REAL, POINTER, SAVE, DIMENSION(:,:,:,:) :: LinSrcSlope ! Linear Source Approx. Slope

REAL, POINTER, DIMENSION(:,:,:)   :: psiSlope
REAL, POINTER, DIMENSION(:,:,:,:) :: srcSlope, phisSlope

! Domain Decomposition
REAL, POINTER, DIMENSION(:,:,:,:,:,:) :: AsyPhiAngIn

! CMFD
TYPE (PinXS_Type), POINTER, DIMENSION(:,:) :: PinXS

REAL, POINTER, SAVE, DIMENSION(:,:)     :: PsiFm, PsiFmD
REAL, POINTER, SAVE, DIMENSION(:,:,:)   :: PhiC, PhiFm
REAL, POINTER, SAVE, DIMENSION(:,:,:,:) :: AxDtil, AxDhat, AxPDhat

! Group Condensed CMFD
TYPE (PinXS_TYPE),  POINTER, DIMENSION(:,:) :: GcPinXS
TYPE (CMFDLS_TYPE), POINTER, DIMENSION(:)   :: GcCMFDLS

REAL, POINTER, DIMENSION(:,:)   :: GcPsiC, GcPsicD
REAL, POINTER, DIMENSION(:,:,:) :: GcPhiC

! Critical Spectrum
REAL, POINTER, SAVE, DIMENSION(:) :: PhiCrit  ! Critical Spectrum
REAL, POINTER, SAVE, DIMENSION(:) :: SpecConv ! Conversion Factor
! ------------------------------------------------------------------------------------------------------------
! Transient Variables
REAL, POINTER, SAVE, DIMENSION(:,:) :: TranPsi, TranPsid ! Transient Fission Source
REAL, POINTER, SAVE, DIMENSION(:,:) :: TranPsiCm, TranPsiCmd, TranPsiFm, TranPsiFmd

REAL, POINTER, SAVE, DIMENSION(:,:,:) :: TranPhi, TranPhiCm, TranPhiFm
REAL, POINTER, SAVE, DIMENSION(:,:,:) :: Prec, PrecCm, PrecFm ! Precursor Number Density
REAL, POINTER, SAVE, DIMENSION(:,:,:) :: ResSrcCm, ResSrcFm   ! Theta Method Source 
REAL, POINTER, SAVE, DIMENSION(:,:,:) :: GcTranSrc

! Adaptive Theta
REAL, POINTER, SAVE, DIMENSION(:,:,:) :: ThetaCM

! BDF Variables
REAL, POINTER, SAVE, DIMENSION(:,:,:) :: TranPhi2, TranPhi3, TranPhi4, TranPhi5
REAL, POINTER, SAVE, DIMENSION(:,:,:) :: TranPhiCm2, TranPhiCm3, TranPhiCm4, TranPhiCm5
REAL, POINTER, SAVE, DIMENSION(:,:,:) :: TranPhiFm2, TranPhiFm3, TranPhiFm4, TranPhiFm5

! SCM Variables
REAL, POINTER, SAVE, DIMENSION(:,:) :: PsiSCM, PsiCSCM

REAL, POINTER, SAVE, DIMENSION(:,:,:) :: TranPrec, ShpFrqCM, ShpFrqCMd, ShpFrqFM, ShpFrqFMd
REAL, POINTER, SAVE, DIMENSION(:,:,:) :: AvgShpFrqCM, AvgShpFrqFM, PrecFrqCM, PrecFrqFM
REAL, POINTER, SAVE, DIMENSION(:,:,:) :: PhiSCM, PhiCSCM, xsnfSCM

! AM3
REAL, POINTER, SAVE, DIMENSION(:,:,:) :: ResSrcCMd, ResSrcFMd
!-------------------------------------------------------------------------------------------------------
TYPE (FMINfo_Type) :: FmInfo
TYPE (CMInfo_TYpe) :: CmInfo
TYPE (THInfo_Type) :: THInfo

TYPE (GroupInfo_Type), SAVE :: GroupInfo, GcGroupInfo

TYPE (CMFDLS_TYPE), POINTER, DIMENSION(:) :: CoreCMFDLs

TYPE (Fxrinfo_type), POINTER, SAVE, DIMENSION(:,:) :: Fxr
! ------------------------------------------------------------------------------------------------------------
END MODULE core_mod