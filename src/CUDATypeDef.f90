#include <CUDADEFINES.h>

#ifdef __PGI

MODULE CUDATypeDef

USE TYPEDEF,	    ONLY : PinXS_Type
USE TYPEDEF_COMMON
USE CSRMATRIX
USE BSRMATRIX
USE CUDAFOR
USE CUBLAS
USE CUSPARSE
IMPLICIT NONE

TYPE cuGeometry_Type
  !--- Core Geometry Variables
  TYPE(superPin_Type), POINTER :: superPin(:)
  TYPE(cuCell_Type), POINTER :: Cell(:)
  REAL, POINTER :: PinVolFm(:, :), hz(:), hzfm(:), FsrVol(:, :)
  INTEGER, POINTER :: pinMap(:), pinMapRev(:), planeMap(:)
  INTEGER, POINTER :: fmRange(:, :)
  INTEGER, POINTER :: GcStruct(:, :), GcStructInv(:)
  INTEGER, ALLOCATABLE :: AxType(:)
  INTEGER, ALLOCATABLE :: Fsr2Fxr(:, :)
  INTEGER :: nfxr, nfsr, ng, ngc = 2, nx, ny, nxy, nxyc, nxya, nzCMFD, nAxType, nSubplane = 1
  INTEGER :: myzb, myze
  INTEGER :: RadBC(4), AxBC(2)
  LOGICAL :: l3dim
  LOGICAL, POINTER :: lRefPin(:), lRefCell(:, :), lRefPlane(:)
  LOGICAL, POINTER :: lH2OCell(:, :)
  !--- Ray Geometry Variables
  INTEGER, ALLOCATABLE :: RotRayCount(:), RotRayList(:, :)
  INTEGER :: nRotRay, nModRay, nAziAngle, nPolarAngle, nPolar1D = 5, nPhiAngSv
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
  INTEGER :: CMFDSolver = 1   !--- 1 : BiCGSTAB, 2 : Stationary
  INTEGER :: AxSolver = 1     !--- 1 : Nodal,    2 : MOC
  INTEGER :: DcplLv = 0
  REAL :: MOCHeight, CMFDHeight
  !--- Parallel Control Variables
  INTEGER :: nDevice = 1, nMOCHelper, nCMFDHelper
  LOGICAL :: lMulti = .FALSE.
  REAL :: Shift = 0.25
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
  TYPE(CSR_CMFD_PRECISION) :: M, D, ILU, SPAI
  TYPE(CSR_CMFD_PRECISION) :: jM, rbM(2), rbDiag(2)   !--- Internal Matrices for Solvers
  TYPE(BSR_CMFD_PRECISION) :: rbDiagB(2)
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
  REAL(8), POINTER, DEVICE :: phis8(:, :), src8(:, :), psi8(:), psid8(:)
  REAL(8), ALLOCATABLE :: h_phis8(:, :, :), h_phic8(:, :, :), h_neighphis8(:, :, :)
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

END MODULE

#endif
