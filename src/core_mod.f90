MODULE Core_mod
USE PARAM
USE TYPEDEF, ONLY : Fxrinfo_type,   PinXs_Type,    GroupInfo_Type,       &
                    CMFDLS_TYPE,    AXFLX_TYPE,                          &
                    CMInfo_Type,    FMInfo_Type,                         &
                    THInfo_Type

INTEGER :: nCoreFsr, nCoreFxr, nPhiAngSv

REAL :: eigv = 1.0_8
REAL :: peigv = 1.0_8
REAL :: eigv_adj = 1.0_8
REAL :: peigv_adj = 1.0_8
REAL :: xkconv, xkcrit                                      !Converged K, Critical K

REAL, POINTER, SAVE :: PHIS(:,:,:)                          !Scalar Flux Vectors
REAL, POINTER, SAVE :: PSI(:,:), PSID(:,:)                  !Fsr FIssion Scource
REAL, POINTER, SAVE :: Power(:, :)                           !Power
REAL, POINTER, SAVE :: PSIC(:,:), PSICD(:,:)                !Cell Fission Source Term

REAL, POINTER, SAVE :: LinSrcSlope(:, :, :, :)                 !Linear Source Approx. Slope

!--- CNJ Edit : CASMO Linear Source
REAL, POINTER :: srcSlope(:, :, :, :)
REAL, POINTER :: phisSlope(:, :, :, :)
REAL, POINTER :: psiSlope(:, :, :)

!REAL, POINTER, SAVE :: tSrc(:,:)                           !
!REAL, POINTER, SAVE :: xst1g(:,:)                          !transport XS
REAL, POINTER, SAVE :: PhiAngin(:,:,:,:)                    !Angular Flux data at the boundary

!--- CNJ Edit : Domain Decomposition
REAL, POINTER :: AsyPhiAngIn(:, :, :, :, :, :)

REAL, POINTER, SAVE :: RadJout(:, :, :, :, :)               !Moc Jout (Surface, Pin, Plane, Group)
REAL, POINTER, SAVE :: AxSrc(:, :, :), AxPXS(:, :, :)       !Axial Source, Psuedo Absorption XS
REAL, POINTER, SAVE :: wmoc(:)
!REAL, POINTER, SAVE :: MocJout1g(:, :, :, :)
!REAL, POINTER, SAVE ::
 
!REAL, POINTER, SAVE :: PHIS1g(:, :)

TYPE(PinXS_Type), POINTER :: PinXS(:, :)
REAL, POINTER, SAVE :: PhiC(:, :, :)
REAL, POINTER, SAVE :: PhiFm(:, :, :), PsiFm(:, :), PsiFmD(:, :)
REAL, POINTER, SAVE :: phim(:, :, :, :)
REAL, POINTER, SAVE :: AxDtil(:, :, :, :), AxDhat(:, :, :, :), AxPDhat(:, :, :, :)
!Group Condensed CMFD
TYPE(PinXS_TYPE), POINTER :: GcPinXS(:, :)      !Group Condensed PinXS
TYPE(CMFDLS_TYPE), POINTER :: GcCMFDLS(:)
!TYPE(CMFDLS_TYPE), POINTER :: CMFD2GLS(:)
REAL, POINTER :: GcPsiC(:, :), GcPsicD(:, :)
REAL, POINTER :: GcPhiC(:, :, :)
!Critical Spectrum
REAL, POINTER, SAVE :: PhiCrit(:), SpecConv(:)              ! Critical Spectrum and Conversion Factor

!Transient Variables------------------------------------------------------------------------------------
REAL, POINTER, SAVE :: TranPhi(:, :, :)
REAL, POINTER, SAVE :: TranPhiCm(:, :, :), TranPhiFm(:, :, :)   

REAL, POINTER, SAVE :: Prec(:, :, :)                         ! Precursor Number Density
REAL, POINTER, SAVE :: PrecCm(:, :, :) ,PrecFm(:, :, :)      ! Precursor Number Density

REAL, POINTER, SAVE :: TranPsi(:, :), TranPsid(:, :)            ! Transient Fission Source
REAL, POINTER, SAVE :: TranPsiCm(:, :), TranPsiCmd(:, :)
REAL, POINTER, SAVE :: TranPsiFm(:, :), TranPsiFmd(:, :)

REAL, POINTER, SAVE :: ResSrcCm(:, :, :), ResSrcFm(:, :, :)   !Theta Method Source 

REAL, POINTER, SAVE :: GcTranSrc(:, :, :)

!Adaptive Theta
REAL, POINTER, SAVE :: ThetaCM(:, :, :)

!BDF Variables
REAL, POINTER, SAVE :: TranPhi2(:, :, :), TranPhi3(:, :, :), TranPhi4(:, :, :), TranPhi5(:, :, :)
REAL, POINTER, SAVE :: TranPhiCm2(:, :, :), TranPhiCm3(:, :, :), TranPhiCm4(:, :, :), TranPhiCm5(:, :, :)
REAL, POINTER, SAVE :: TranPhiFm2(:, :, :), TranPhiFm3(:, :, :), TranPhiFm4(:, :, :), TranPhiFm5(:, :, :)

!SCM Variables
REAL, POINTER, SAVE :: TranPrec(:, :, :)
REAL, POINTER, SAVE :: ShpFrqCM(:, :, :), ShpFrqCMd(:, :, :)
REAL, POINTER, SAVE :: ShpFrqFM(:, :, :), ShpFrqFMd(:, :, :)
REAL, POINTER, SAVE :: AvgShpFrqCM(:, :, :), AvgShpFrqFM(:, :, :)
REAL, POINTER, SAVE :: PrecFrqCM(:, :, :), PrecFrqFM(:, :, :)

REAL, POINTER, SAVE :: PhiSCM(:, :, :), PhiCSCM(:, :, :)
REAL, POINTER, SAVE :: PsiSCM(:, :), PsiCSCM(:, :)
REAL, POINTER, SAVE :: xsnfSCM(:, :, :)

!AfSrc
REAL, POINTER, SAVE :: TranPhi1a(:, :, :, :, :), TranPhi2a(:, :, :, :, :)

!AM3
REAL, POINTER, SAVE :: ResSrcCMd(:, :, :), ResSrcFMd(:, :, :)
!-------------------------------------------------------------------------------------------------------

TYPE(CMFDLS_TYPE), POINTER :: CoreCMFDLs(:)
TYPE(Fxrinfo_type), POINTER, SAVE :: Fxr(:,:)
TYPE(GroupInfo_Type),SAVE :: GroupInfo, GcGroupInfo
!TYPE(AXFLX_TYPE), POINTER, SAVE :: AxFlx(:, :)
TYPE(FMINfo_Type) :: FmInfo
TYPE(CMInfo_TYpe) :: CmInfo
TYPE(THInfo_Type) :: THInfo

END MODULE