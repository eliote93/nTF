#include <CUDADEFINES.h>

#ifdef __PGI

MODULE CUDATypeDef

USE TYPEDEF,	    ONLY : PinXS_Type
USE TYPEDEF_COMMON
USE CSRMATRIX
USE CUDAFOR
USE CUBLAS
USE CUSPARSE
IMPLICIT NONE

TYPE cuGeometry_Type
  !--- Core Geometry Variables
  TYPE(superPin_Type), POINTER :: superPin(:)
  TYPE(cuCell_Type), POINTER :: Cell(:)
  REAL, POINTER :: PinVolFm(:, :), hz(:), hzfm(:), FsrVol(:, :)
  REAL, POINTER :: PinVolFmF(:, :) ! Fuel Pin Volume
  INTEGER, POINTER :: pinMap(:), pinMapRev(:), planeMap(:)
  INTEGER, POINTER :: fmRange(:, :)
  INTEGER, POINTER :: GcStruct(:, :), GcStructInv(:)
  INTEGER, ALLOCATABLE :: AxType(:)
  INTEGER, ALLOCATABLE :: Fsr2Fxr(:, :)
  INTEGER :: nfxr, nfsr, ng, ngc = 2, nx, ny, nxy, nxyc, nxya, nzCMFD, nAxType, nSubplane = 1
  INTEGER :: myzb, myze
  INTEGER :: RadBC(4), AxBC(2)
  LOGICAL :: l3dim
  LOGICAL :: lRot
  LOGICAL, POINTER :: lRefPin(:), lRefCell(:, :), lRefPlane(:)
  LOGICAL, POINTER :: lH2OCell(:, :)
  !--- Ray Geometry Variables
  INTEGER, ALLOCATABLE :: RotRayCount(:), RotRayList(:, :)
  INTEGER :: nRotRay, nModRay, nAziAngle, nPolarAngle, nPolar1D = 5, nPhiAngSv
  INTEGER :: nprec
END TYPE

TYPE cuCntl_Type
  !--- Calculation Option Variables
  LOGICAL :: lAxial = .TRUE.
  LOGICAL :: lRay1D = .TRUE.
  LOGICAL :: lDcpl = .FALSE.
  LOGICAL :: lNatural = .TRUE.
  LOGICAL :: lSuperpin = .TRUE.
  LOGICAL :: lSENM = .TRUE.
  LOGICAL :: lGcCMFD = .FALSE.
  LOGICAL :: lCASMO = .TRUE.
  LOGICAL :: lAsync = .TRUE.
  LOGICAL :: lSubplane = .FALSE.
  LOGICAL :: lShift = .TRUE.
  LOGICAL :: lSPAI = .FALSE.
  LOGICAL :: lPwDist = .FALSE.
  INTEGER :: AxSolver = 1     !--- 1 : Nodal,    2 : MOC
  INTEGER :: DcplLv = 0
  INTEGER :: PCRtype = 1
  REAL :: MOCHeight, CMFDHeight
  !--- Parallel Control Variables
  INTEGER :: nDevice = 1, nMOCHelper, nCMFDHelper
  LOGICAL :: lMulti = .FALSE.
  REAL :: Shift = 0.25
  !--- cuXS related parameter
  LOGICAL :: lGrpBlock = .TRUE.
END TYPE

TYPE cuMOC_Type
  REAL(GPU_FLUX_PRECISION), ALLOCATABLE, DEVICE :: phis(:, :), phia(:, :, :, :), phim(:, :, :), jout(:, :, :, :)
  REAL(GPU_SOURCE_PRECISION), ALLOCATABLE, DEVICE :: psi(:), src(:, :), srcm(:, :, :), SrcAng(:, :, :, :)
  REAL(GPU_PRECISION), ALLOCATABLE, DEVICE :: xst(:, :), xssm(:, :, :), chi(:, :)
  REAL(GPU_PRECISION), ALLOCATABLE, DEVICE :: PhiAngIn(:, :, :)
  !--- Domain Connected Pn Ray Tracing Variables
  REAL(GPU_PRECISION), POINTER, DEVICE :: PhiAngInMg(:, :, :)
  REAL(GPU_FLUX_PRECISION), POINTER, DEVICE :: phisMg(:, :), phimMg(:, :, :), phiaMg(:, :, :, :), joutMg(:, :, :, :)
  REAL(GPU_SOURCE_PRECISION), ALLOCATABLE, DEVICE :: xstMg(:, :), srcMg(:, :), srcmMg(:, :, :), SrcAngMg(:, :, :, :)
END TYPE

TYPE cuCMFD_Type
  TYPE(CSR_SPARSITY) :: spM, spD
  TYPE(CSR_CMFD_PRECISION) :: M, D, ILU, SPAI
  TYPE(CSR_CMFD_PRECISION) :: rbM(2), rbDiag(2)   !--- Internal Matrices for Solvers
  TYPE(CSR_DOUBLE) :: S, F, Chi
  REAL(GPU_CMFD_PRECISION), POINTER, DEVICE :: invDiag(:, :), offDiag(:, :, :), rOffDiag(:, :), bOffDiag(:, :)
  REAL(8), ALLOCATABLE :: AxDtil(:, :, :, :), AxDhat(:, :, :, :)
  REAL(GPU_CMFD_PRECISION) :: wSOR
  TYPE(PinXS_Type), POINTER :: PinXS(:, :)
  INTEGER :: ng, nIter
  LOGICAL :: lFirst = .TRUE.
  !--- MPI Control Variables
  INTEGER :: bottomRange(2), bottomRangeRB(2, 2)
  INTEGER :: topRange(2), topRangeRB(2, 2)
  INTEGER, POINTER :: planeMap(:)
  !--- Double Precision Variables for Residual Correction (Iterative Refinement)
  REAL(8), POINTER, DEVICE :: offDiag8(:, :, :)
  REAL(8), POINTER, DEVICE :: phis8(:, :), src8(:, :), srcd8(:,:), psi8(:), psid8(:)
  REAL(8), ALLOCATABLE :: h_phis8(:, :, :), h_phic8(:, :, :), h_neighphis8(:, :, :)
END TYPE

TYPE cuTranCMInfo_Type
  !TYPE(CSR_DOUBLE) :: Omegalm, Omegal0
  REAL(8) :: Inv_FuelVol
  REAL(8), POINTER :: CellOmegam(:, :, :), CellOmega0(:, :, :), CellOmegap(:, :, :)
  REAL(8), POINTER :: Prec(:, :, :)
  REAL(8), POINTER :: PhiC(:, :, :), TranPhiC(:, :, :)
  REAL(8), POINTER, DEVICE :: dPrec(:, :), dOmegam(:, :), dOmega0(:, :), dOmegap(:, :)
  REAL(8), POINTER, DEVICE :: TranPhi(:, :), TranPsi(:), TranPsid(:)
  REAL(8), POINTER, DEVICE :: Expo_Alpha(:, :), Expo(:, :), VolinvVel(:, :)
  REAL(8), POINTER, DEVICE :: ResSrc(:, :), PrecSrc(:), PrecSrcK(:,:)
  REAL(8), POINTER, DEVICE :: Omegalm(:), Omegal0(:)
  !REAL(8), POINTER, DEVICE :: OmegalmK(:), Omegal0K(:)
END TYPE

TYPE cuPwDist_Type
  REAL(8), POINTER, DEVICE :: pinPw(:)   ! normalized power distribution(=RelPower)
  REAL(8), POINTER, DEVICE :: pinPwf(:)  ! fission power (not normalized)
  REAL(8), POINTER, DEVICE :: hPrec(:, :)! decay heat precursor
  REAL(8), POINTER, DEVICE :: invPinVolFm(:)
  LOGICAL :: lTran = .FALSE.
END TYPE
TYPE cuAxial_Type
  TYPE(CSR_NODAL_PRECISION) :: S, F, Chi
  REAL(GPU_NODAL_PRECISION), ALLOCATABLE, DEVICE :: S0(:, :, :, :), S1(:, :, :, :), Fa(:, :, :), Chia(:, :, :)
  REAL(GPU_NODAL_PRECISION), ALLOCATABLE, DEVICE :: phiShape(:, :, :), phiCoeff(:, :, :, :), phim(:, :, :)
  REAL(GPU_NODAL_PRECISION), ALLOCATABLE, DEVICE :: PhiAngIn(:, :, :, :), PhiAngOut(:, :, :, :)
  REAL(GPU_NODAL_PRECISION), ALLOCATABLE, DEVICE :: srcCoeff(:, :, :, :), srcm(:, :, :)
  REAL(GPU_NODAL_PRECISION), ALLOCATABLE, DEVICE :: lkgCoeff(:, :, :, :), psiCoeff(:, :, :)
  REAL(GPU_NODAL_PRECISION), ALLOCATABLE, DEVICE :: D(:, :, :), sigR(:, :, :), xst(:, :, :)
  REAL(GPU_NODAL_PRECISION), ALLOCATABLE, DEVICE :: E1(:, :, :, :), E3(:, :, :, :), R1(:, :, :, :), R3(:, :, :, :)
  REAL, ALLOCATABLE, DEVICE :: phic(:, :, :), Jout(:, :, :, :, :)
  !--- Transient Variables
  REAL(GPU_NODAL_PRECISION), ALLOCATABLE, DEVICE :: InvVel(:, :, :)
  REAL(GPU_NODAL_PRECISION), ALLOCATABLE, DEVICE :: TranPhi(:, :, :)
  REAL(GPU_NODAL_PRECISION), ALLOCATABLE, DEVICE :: TranPsi(:, :), TranPsid(:, :)
  REAL(GPU_NODAL_PRECISION), ALLOCATABLE, DEVICE :: Omegam(:, :, :), Omega0(:, :, :), Omegap(:, :, :)
  REAL(GPU_NODAL_PRECISION), ALLOCATABLE, DEVICE :: Omegalm(:, :), Omegal0(:, :)
  REAL(GPU_NODAL_PRECISION), ALLOCATABLE, DEVICE :: OmegalmK(:, :, :), Omegal0K(:, :, :)
  REAL(GPU_NODAL_PRECISION), ALLOCATABLE, DEVICE :: Prec(:, :, :)
  REAL(GPU_NODAL_PRECISION), ALLOCATABLE, DEVICE :: PrecSrc(:, :), ResSrc(:, :, :), FixedSrc(:, :, :)
  REAL(GPU_NODAL_PRECISION), ALLOCATABLE, DEVICE :: PrecSrcK(:, :, :)
  REAL(GPU_NODAL_PRECISION), ALLOCATABLE, DEVICE :: Expo(:, :, :), Expo_Alpha(:, :, :)
END TYPE

TYPE cuDevice_Type
  TYPE(cublasHandle) :: myblasHandle
  TYPE(cusparseHandle) :: mySparseHandle
  INTEGER(KIND = cuda_stream_kind) :: myStream
  !--- Device Property Variables
  INTEGER :: cuSMXCount
  INTEGER :: cuArchitecture
  INTEGER :: cuWarpSize
  INTEGER :: cuMaxThreadPerSMX
  INTEGER :: cuMaxThreadPerBlock
  INTEGER :: cuMaxBlockPerSMX
  INTEGER :: cuMaxWarpPerSMX
  INTEGER :: cuWarpPerBlock
  INTEGER :: cuThreadPerBlock
  INTEGER :: sharedMemoryDim
  !--- 1D MOC Control Variables
  INTEGER :: nzSub, nDiv(1000)
  INTEGER, ALLOCATABLE :: cmSubrange(:, :), fmSubrange(:, :)
  INTEGER, ALLOCATABLE :: cmMap(:), fmMap(:)
  REAL(GPU_PRECISION), ALLOCATABLE :: hzSub(:)
  !--- 2D MOC Control Variables
  INTEGER :: myzb, myze, myAxType, nz
  INTEGER :: nRotRay, RotRayBeg, RotRayEnd
  LOGICAL :: lRayStatic = .TRUE.
  !--- 3D CMFD Control Variables
  INTEGER, POINTER :: pinMapRB(:), pinMapRevRB(:, :), planeMapRB(:)
  INTEGER :: myzbf, myzef, nzCMFD
  INTEGER :: rbBeg(2), rbEnd(2), rbRange(2, 100, 2), nxyRB(100, 2), nxyzRB(2)
  INTEGER :: bottomRange(2), bottomRangeRB(2, 2)
  INTEGER :: topRange(2), topRangeRB(2, 2)
  LOGICAL :: lTop, lBottom, lFuel
  !--- Options
  LOGICAL :: lFullWarp = .TRUE.
END TYPE

TYPE cuRotRay1D_Type
  INTEGER, ALLOCATABLE :: PhiAngInSvIdx(:, :), PhiAngOutSvIdx(:, :)
  INTEGER, ALLOCATABLE :: RotRayIdxSt(:), CoreRayIdxSt(:), PinRayIdxSt(:)
  INTEGER, ALLOCATABLE :: iAzi(:), iDir(:)
  INTEGER, ALLOCATABLE :: PinIdx(:), SurfIdx(:, :)
  INTEGER, ALLOCATABLE :: FsrIdx(:)
  REAL(GPU_PRECISION), ALLOCATABLE :: LenSeg(:)
END TYPE

TYPE cuDepl_Type
  INTEGER :: SysByte = 500, Scale = 6
  INTEGER :: nSubStep = 8
END TYPE

TYPE cuCoreXsMac_Type
  REAL(GPU_XS_PRECISION), DEVICE, ALLOCATABLE :: chi(:, :)
  REAL(GPU_XS_PRECISION), DEVICE, ALLOCATABLE :: XSt(:, :)
  REAL(GPU_XS_PRECISION), DEVICE, ALLOCATABLE :: XSS(:, :)
  REAL(GPU_XS_PRECISION), DEVICE, ALLOCATABLE :: XSStr(:, :)
  REAL(GPU_XS_PRECISION), DEVICE, ALLOCATABLE :: XStr(:, :)
  REAL(GPU_XS_PRECISION), DEVICE, ALLOCATABLE :: XSa(:, :)
  REAL(GPU_XS_PRECISION), DEVICE, ALLOCATABLE :: XSnf(:, :)
  REAL(GPU_XS_PRECISION), DEVICE, ALLOCATABLE :: XSkf(:, :)
  REAL(GPU_SM_PRECISION), DEVICE, ALLOCATABLE :: XSsm(:, :), XSsmP1(:, :), XSsmP2(:, :), XSsmP3(:, :)
END TYPE

TYPE cuBlockFxr_Type
  INTEGER :: nFxrD, nFxrTri, NtTri, nResD, nResTri
  INTEGER, DEVICE, ALLOCATABLE :: mapglobalidD(:), mapglobalidT(:)
  INTEGER, DEVICE, ALLOCATABLE :: mapR2GD(:), mapR2GT(:), mapG2RD(:), mapG2RT(:)
  INTEGER, DEVICE, ALLOCATABLE :: mapC2GD(:), mapC2GT(:), mapG2CD(:), mapG2CT(:)
  INTEGER, DEVICE, ALLOCATABLE :: mapGR2GD(:), mapGR2GT(:), mapG2GRD(:), mapG2GRT(:)
  INTEGER, DEVICE, ALLOCATABLE :: mapA2GT(:), mapG2AT(:)
  INTEGER, DEVICE, ALLOCATABLE :: nisoD(:), nisoT(:), AccNT(:), ACCNtR(:)
  INTEGER, DEVICE, ALLOCATABLE :: MapNuclD(:,:), MapNuclT(:), MapNuclResD(:,:), MapNuclResT(:)
  INTEGER, DEVICE, ALLOCATABLE :: fsrstD(:), fsrstT(:), nfsrD(:), nfsrT(:), pinidD(:), pinidT(:), mapfsrD(:), mapfsrT(:)
  INTEGER, DEVICE, ALLOCATABLE :: itTempD(:,:), itTempPxD(:,:), itTempT(:,:), itTempPxT(:,:)
  REAL(GPU_XS_PRECISION), DEVICE, ALLOCATABLE :: wtTempD(:,:), wtTempPxD(:,:), wtTempT(:,:), wtTempPxT(:,:)
  REAL(GPU_PNUM_PRECISION), DEVICE, ALLOCATABLE :: pnumD(:,:), pnumT(:), indD(:,:)
  REAL(GPU_XS_PRECISION), DEVICE, ALLOCATABLE :: chi(:,:), areaD(:), areaT(:), fsrvolD(:), fsrvolT(:)
  REAL(GPU_XS_PRECISION), DEVICE, ALLOCATABLE :: tempD(:), tempT(:)
  REAL(GPU_RES_PRECISION), DEVICE, ALLOCATABLE :: FresoAD(:,:), FresoAT(:,:), FresoF(:,:), FresoNF(:,:), FresoKF(:,:)
  REAL(GPU_RES_PRECISION), DEVICE, ALLOCATABLE :: FresoAIsoD(:,:,:), FresoFIsoD(:,:,:), FresoAIsoT(:)
  LOGICAL, DEVICE, ALLOCATABLE :: lvoidD(:), lvoidT(:)
  LOGICAL, DEVICE, ALLOCATABLE :: lresD(:), lresT(:), lcldD(:), lcldT(:)
  REAL(GPU_RES_PRECISION), DEVICE, ALLOCATABLE :: FnAdjD(:,:), FnAdjT(:,:), FtAdjD(:,:,:), FtAdjT(:,:,:)
  REAL(GPU_RES_PRECISION), DEVICE, ALLOCATABLE :: xseq_c_1gD(:,:), xseq_c_1gT(:,:), xseq_f_1gT(:,:)
  REAL(GPU_RES_PRECISION), DEVICE, ALLOCATABLE :: xseq_f_mgD(:,:,:), xseq_f_mgT(:,:,:)
END TYPE

TYPE cuIsoData_Type
  INTEGER :: ntiso, ng
  REAL(GPU_XS_PRECISION), DEVICE, ALLOCATABLE :: ptsTemp(:), ptsTempPx(:)
  INTEGER, DEVICE, ALLOCATABLE :: ptrTemp(:), ptrTempPx(:)

  REAL(GPU_XS_PRECISION), DEVICE, DIMENSION(:,:), ALLOCATABLE :: siga, sigf, sigkf, signf, sigs, sigstr
  INTEGER, DEVICE, DIMENSION(:,:), ALLOCATABLE :: gidSmp0, gidSmp1, gidSmp2, gidSmp3
  INTEGER, DEVICE, DIMENSION(:,:), ALLOCATABLE :: ptrSmp0, ptrSmp1, ptrSmp2, ptrSmp3
  REAL(GPU_XS_PRECISION), DEVICE, ALLOCATABLE :: smp0(:), smp1(:), smp2(:), smp3(:)
END TYPE

TYPE cuResIsoData_Type
  INTEGER :: ntrtemp, ntsig0, ntlv, nresiso, ntempPot, ntempLv
  LOGICAL, DEVICE, ALLOCATABLE :: lclad(:)
  INTEGER, DEVICE, ALLOCATABLE :: ityp(:), rifid(:,:), IsoG2R(:)
  INTEGER, DEVICE, ALLOCATABLE :: ptrRTemp(:), ptrSig0(:), ptrNlv(:)
  INTEGER, DEVICE, ALLOCATABLE :: ptrTempPot(:), ptrTempLv(:), ptrRifRat(:), ptrRifRx(:)
  REAL(GPU_RES_PRECISION), DEVICE, ALLOCATABLE :: xsalog(:,:), ri_a(:,:)
  REAL(GPU_XS_PRECISION), DEVICE, ALLOCATABLE :: sigp(:), lamsigp(:,:), sig0sq(:), rtempsq(:)
  REAL(GPU_XS_PRECISION), DEVICE, ALLOCATABLE :: lvabs(:,:), lvfis(:,:), wgtabs(:,:), wgtfis(:,:), ratlog(:)
  REAL(GPU_RES_PRECISION), DEVICE, ALLOCATABLE :: rifabs(:,:), riffis(:,:)
END TYPE

TYPE cuMLGLib_Type
  INTEGER :: f_nmaclv, f_nmaclv1G, c_nmaclv1G
  REAL(GPU_XS_PRECISION), DEVICE, ALLOCATABLE :: c_maclv1G(:), c_maclv1G_log(:), f_maclv1G(:,:), f_maclv1G_log(:,:)
  REAL(GPU_XS_PRECISION), DEVICE, ALLOCATABLE :: f_maclv(:,:), f_maclv_log(:,:)
END TYPE

TYPE cuResPin_Type
  LOGICAL, DEVICE, ALLOCATABLE :: lCellRes(:), lAIC(:), lPinRes(:), lresA(:), lresC(:)
  REAL(GPU_RES_PRECISION), DEVICE, ALLOCATABLE :: avgxseq_mg(:,:,:), avgxseq_1g(:,:)
  REAL(GPU_RES_PRECISION), DEVICE, ALLOCATABLE :: rifa(:,:,:), rifs(:,:,:), riff(:,:,:)
  REAL(GPU_RES_PRECISION), DEVICE, ALLOCATABLE :: FnAdj(:,:)
  REAL(GPU_PNUM_PRECISION), DEVICE, ALLOCATABLE :: ind(:,:)
  INTEGER, DEVICE, ALLOCATABLE :: niso(:), igt(:), idiso(:,:), mapnucl(:,:)
  REAL(GPU_PNUM_PRECISION), DEVICE, ALLOCATABLE :: pnum(:,:), temp(:)
END TYPE

TYPE cuPinInfo_Type
  INTEGER, DEVICE, ALLOCATABLE :: FxrStD(:), FxrStT(:), GenStD(:), GenStT(:), AICStT(:)
END TYPE

END MODULE

#endif
