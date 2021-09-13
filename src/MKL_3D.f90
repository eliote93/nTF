#include <defines.h>
#ifdef __INTEL_MKL
MODULE MKL_3D

USE PARAM
USE TYPEDEF,        ONLY : PinXS_Type
USE TYPEDEF_COMMON

#ifdef __GAMMA_TRANSPORT
USE GammaTYPEDEF,   ONLY : GPinXS_Type
#endif

USE CSRMATRIX
USE OMP_LIB

IMPLICIT NONE

INCLUDE 'mpif.h'
INCLUDE 'mkl_blas.fi'

INTEGER, PARAMETER :: in = 1, out = 2, surf = 3
INTEGER, PARAMETER :: NODAL = 1, FDM = 2, MOC = 3

TYPE triLU_Type
  INTEGER :: nx
  REAL, POINTER :: lower(:), diag(:), upper(:)
END TYPE

TYPE blockMat_Type
  INTEGER :: ny
  TYPE(triLU_Type), POINTER :: blockDiag(:), Del(:)
  REAL, POINTER :: lower(:, :), upper(:, :)
END TYPE

TYPE mklCMFD_Type
  TYPE(CSR_DOUBLE), POINTER, DIMENSION(:) :: M
  TYPE(CSR_DOUBLE), POINTER, DIMENSION(:) :: ILU
  TYPE(CSR_DOUBLE), POINTER, DIMENSION(:) :: SPAI
  
  TYPE(mklDavidson_Type), POINTER :: Davidson
  
  TYPE(mklJacobi_Type),   POINTER, DIMENSION(:)   :: Jacobi
  TYPE(PinXS_Type),       POINTER, DIMENSION(:,:) :: PinXS
  
#ifdef __GAMMA_TRANSPORT
  TYPE(GPinXS_Type), POINTER, DIMENSION(:,:) :: GPinXS
#endif

  REAL, POINTER, DIMENSION(:,:)       :: F, Chi, src, psi, psid
  REAL, POINTER, DIMENSION(:,:,:)     :: S, phis, phic, neighphis, theta, AxOffDiag, trAxOffDiag
  REAL, POINTER, DIMENSION(:,:,:,:)   :: phisd, AxDtil, AxDhat, dcplAxOffDiag
  REAL, POINTER, DIMENSION(:,:,:,:,:) :: partialAxDhat, superJout
  
  INTEGER, POINTER, DIMENSION(:)   :: planeMap
  INTEGER, POINTER, DIMENSION(:,:) :: InScatRange
  
  INTEGER, DIMENSION(2) :: bottomRange, topRange
  INTEGER :: ng
  
END TYPE mklCMFD_Type

TYPE mklDavidson_Type
  TYPE(CSR_DOUBLE) :: M, ILU, F
  !--- Generalized Davidson Method Variables
  REAL, POINTER :: u(:), t(:), r(:)
  !--- Rayleigh-Ritz Procedure Variables
  REAL, POINTER :: Vmat(:, :), Mr(:, :), Fr(:, :)
END TYPE

TYPE mklJacobi_Type
  TYPE(CSR_DOUBLE) :: M
  REAL, POINTER :: invDiag(:), offDiag(:, :)
END TYPE

TYPE mklBILU_Type
  TYPE(blockMat_Type), POINTER :: blockMat(:)
  REAL, POINTER :: lower(:, :), upper(:, :)
END TYPE

TYPE mklAxial_Type
  TYPE(mklAngle_Type), POINTER :: Angle(:)
  REAL, POINTER :: SmP1(:, :, :, :), SmP2(:, :, :, :), SmP3(:, :, :, :)
  REAL, POINTER :: phic(:, :, :), PhiAngIn(:, :, :, :), PhiAngOut(:, :, :, :), Jout(:, :, :, :, :)
  REAL, POINTER :: atil(:, :, :, :)
  INTEGER :: ng
END TYPE

TYPE mklAngle_Type
  REAL :: cosv, rcosv
  REAL :: wt, wtsurf
  REAL :: wt2D(2)
  INTEGER :: polar2D(2)
END TYPE

TYPE mklGeom_Type
  INTEGER :: myzb, myze
  INTEGER :: ng, ngc = 2, ngg, nFsr, nxmax, ny, nxy, nz, nzCMFD, nSubplane = 1, nPolar1D = 10, nPolar2D
  INTEGER :: AxBC(2)
  INTEGER, POINTER :: nx(:), nzfm(:)
  INTEGER, POINTER :: ixRange(:, :), pinRange(:, :), fmRange(:, :), fxrRange(:, :)
  INTEGER, POINTER :: planeMap(:), pinMap(:), pinMapRev(:)
  INTEGER, POINTER :: GcStruct(:, :), GcStructInv(:)
  INTEGER, POINTER :: InScatRange(:, :), OutScatRange(:, :)
  REAL, POINTER :: hz(:), hzfm(:), PinVolFm(:, :)
  TYPE(superPin_Type), POINTER :: superPin(:)
  LOGICAL :: lTop, lBottom, l3dim
  LOGICAL, POINTER :: lRefPin(:), lRefCell(:, :), lRefPlane(:)
  LOGICAL, POINTER :: lH2OCell(:, :)
END TYPE

TYPE mklDepl_Type
  INTEGER :: SysByte = 500, Scale = 6
  INTEGER :: nSubStep = 8
END TYPE
TYPE mklCntl_Type
  LOGICAL :: lAxial = TRUE                          !--- Axial Calculation
  LOGICAL :: lAxRefFDM = FALSE                      !--- Axial Reflector FDM
  LOGICAL :: lRefPinFDM = FALSE                     !--- Reflector Pin FDM
  LOGICAL :: lDcpl = FALSE                          !--- Axial Coupling Iteration between MPI Nodes
  LOGICAL :: lSPAI = FALSE                          !--- Sparse Approximate Inverse Preconditioning
  LOGICAL :: lDavidson = FALSE                      !--- Generalized Davidson Eigensolver
  INTEGER :: AxSolver = MOC                         !--- Axial Solver Type
  LOGICAL :: lMOC = FALSE                           !--- Axial MOC Solver
  LOGICAL :: lSENM = FALSE                          !--- Axial SENM Solver
  LOGICAL :: lFDM = FALSE                           !--- Axial FDM Solver
  LOGICAL :: lSP3 = FALSE                           !--- Axial SP3 Diffusion
  LOGICAL :: lCASMO = FALSE                         !--- CASMO Linear Source for Axial MOC Calculation
  LOGICAL :: lPolarXS = FALSE                       !--- Polar XS Homogenization for Axial MOC Calculation
  LOGICAL :: lShift = FALSE                         !--- Wielandt Shift
  LOGICAL :: lChebyshev = FALSE                     !--- Two Parameter Chebyshev Acceleration
  LOGICAL :: lJacobi = FALSE                        !--- Anderson-Jacobi Scheme
  LOGICAL :: lDirect = FALSE                        !--- Direct Solution Employing PARDISO
  LOGICAL :: lPardiso                               !--- Internal Flag for PARDISO
  LOGICAL :: lGcCMFD = FALSE                        !--- Group Condensed CMFD
  LOGICAL :: lSuperpin = TRUE                       !--- Super Cell Homogenization of Gap Pins
  LOGICAL :: lScat1 = FALSE                         !--- Internal Flag for Axial PN MOC
  LOGICAL :: lHyperbolic = TRUE                     !--- Axial MOC Hyperbolic Expansion
  LOGICAL :: lSubplane = FALSE                      !--- CMFD Subplane Scheme
  LOGICAL :: lGamma = FALSE                         !--- Gamma Transport with MKL
  LOGICAL :: pCMFD = FALSE, odCMFD = FALSE          !--- Variants of CMFD
  INTEGER :: DcplLv = 1                             !--- Decoupling Level
  INTEGER :: polyOrder = 0                          !--- Axial MOC Expansion Order
  INTEGER :: scatOrder = 0                          !--- Axial MOC Scattering Order
  INTEGER :: chebyOrder                             !--- Chebyshev Polynomial Order
  INTEGER :: DcplInner = 2                          !--- Inner Iteration for Decoupled CMFD
  INTEGER :: maxOuter = 100, maxInner = 100         !--- Maximum Iterations
  INTEGER :: minOuter = 7, minInner = 3             !--- Minimum Iterations
  INTEGER :: nAxIter = 10                           !--- # of 1-D Ax. MoC Inner Iter.
  REAL :: CMFDHeight, MOCHeight = 0.5               !--- Maximum Height of CMFD and MOC Subplane
  REAL :: outerConv = 0.1, innerConv = 0.01         !--- Convergence Criteria
  REAL :: Shift = 0.25                              !--- Wielandt Shift Value
  LOGICAL :: lDepl = .FALSE.
  LOGICAL :: lEnter = .FALSE.                       !--- Enter MKL CMFD Driver
END TYPE

TYPE(mklCMFD_Type) :: mklCMFD, mklGcCMFD, mklGammaCMFD
TYPE(mklAxial_Type) :: mklAxial, mklGammaAxial
TYPE(mklGeom_Type) :: mklGeom
TYPE(mklCntl_Type) :: mklCntl
TYPE(mklDepl_Type) :: mklDepl

INTEGER, PRIVATE :: Request(4), nRequest

CONTAINS

!--- Blocking Communication Routines

SUBROUTINE GetNeighbor(n, send, recv, dir)
USE PE_MOD,         ONLY : PE
IMPLICIT NONE

REAL :: send(*), recv(*)
INTEGER :: n, dir
INTEGER :: i, irank, ipartner, imod
INTEGER :: tag = 1, info, comm, status(MPI_STATUS_SIZE)

irank = PE%myCMFDrank
comm = PE%MPI_CMFD_COMM
imod = mod(irank + 1, 2)

SELECT CASE (dir)
CASE (top)
  DO i = 1, 0, -1
    IF (imod .EQ. i) THEN
      ipartner = irank + 1
      IF (ipartner .LT. PE%nCMFDproc) THEN
        CALL MPI_SEND(send, n, MPI_DOUBLE_PRECISION, ipartner, tag, comm, info)
      ENDIF
    ELSE
      ipartner = irank - 1
      IF (ipartner .GE. 0) THEN
        CALL MPI_RECV(recv, n, MPI_DOUBLE_PRECISION, ipartner, tag, comm, status, info)
      ENDIF
    ENDIF
    CALL MPI_BARRIER(PE%MPI_CMFD_COMM, info)
  ENDDO
CASE (bottom)
  DO i = 0, 1
    IF (imod .EQ. i) THEN
      ipartner = irank - 1
      IF (ipartner .GE. 0) THEN
        CALL MPI_SEND(send, n, MPI_DOUBLE_PRECISION, ipartner, tag, comm, info)
      ENDIF
    ELSE
      ipartner = irank + 1
      IF (ipartner .LT. PE%nCMFDproc) THEN
        CALL MPI_RECV(recv, n, MPI_DOUBLE_PRECISION, ipartner, tag, comm, status, info)
      ENDIF
    ENDIF
    CALL MPI_BARRIER(PE%MPI_CMFD_COMM, info)
  ENDDO
END SELECT

END SUBROUTINE

SUBROUTINE MatMulMPI(Diag, AxOffDiag, x, y, bottomRange, topRange, op)

IMPLICIT NONE

TYPE(CSR_DOUBLE) :: Diag
REAL :: AxOffDiag(*), x(*), y(*)
INTEGER :: bottomRange(2), topRange(2)
CHARACTER(1), OPTIONAL :: op

REAL, POINTER :: recv(:, :), temp(:, :)
REAL, POINTER :: csrVal(:)
INTEGER, POINTER :: csrRowPtr(:), csrColIdx(:)
INTEGER :: nSelf, nNeigh
CHARACTER(1) :: operate = 'n'

IF (PRESENT(op)) operate = op

csrVal => Diag%csrVal
csrRowPtr => Diag%csrRowPtr
csrColIdx => Diag%csrColIdx

nSelf = Diag%nr
nNeigh = bottomRange(2) - bottomRange(1) + 1

IF (mklCntl%lDcpl .OR. .NOT. mklGeom%l3dim) THEN

  CALL mkl_dcsrgemv(operate, nSelf, csrVal, csrRowPtr, csrColIdx, x, y)

ELSE

  ALLOCATE(recv(nNeigh, 2), temp(nNeigh, 2))

  CALL GetNeighbor(nNeigh, x(bottomRange(1) : bottomRange(2)), recv(:, top), bottom)
  CALL GetNeighbor(nNeigh, x(topRange(1) : topRange(2)), recv(:, bottom), top)

  CALL mkl_dcsrgemv(operate, nSelf, csrVal, csrRowPtr, csrColIdx, x, y)

  CALL vdmul(nNeigh * 2, AxOffDiag, recv, temp)
  IF (.NOT. mklGeom%lBottom) CALL daxpy(nNeigh, 1.0D0, temp(:, bottom), 1, y(bottomRange(1) : bottomRange(2)), 1)
  IF (.NOT. mklGeom%lTop) CALL daxpy(nNeigh, 1.0D0, temp(:, top), 1, y(topRange(1) : topRange(2)), 1)

  DEALLOCATE(recv, temp)

ENDIF

END SUBROUTINE

!--- Non-blocking Communication Routines

SUBROUTINE InitFastComm()

IMPLICIT NONE

nRequest = 0

END SUBROUTINE

SUBROUTINE GetNeighborFast(n, send, recv, dir)
USE PE_MOD,         ONLY : PE
IMPLICIT NONE

REAL :: send(*), recv(*)
INTEGER :: n, dir
INTEGER :: i, irank, ipartner, imod
INTEGER :: tag = 1, info, comm

irank = PE%myCMFDrank
comm = PE%MPI_CMFD_COMM
imod = mod(irank + 1, 2)

SELECT CASE (dir)
CASE (top)
  DO i = 1, 0, -1
    IF (imod .EQ. i) THEN
      ipartner = irank + 1
      IF (ipartner .LT. PE%nCMFDproc) THEN
        nRequest = nRequest + 1
        CALL MPI_ISEND(send, n, MPI_DOUBLE_PRECISION, ipartner, tag, comm, Request(nRequest), info)
      ENDIF
    ELSE
      ipartner = irank - 1
      IF (ipartner .GE. 0) THEN
        nRequest = nRequest + 1
        CALL MPI_IRECV(recv, n, MPI_DOUBLE_PRECISION, ipartner, tag, comm, Request(nRequest), info)
      ENDIF
    ENDIF
  ENDDO
CASE (bottom)
  DO i = 0, 1
    IF (imod .EQ. i) THEN
      ipartner = irank - 1
      IF (ipartner .GE. 0) THEN
        nRequest = nRequest + 1
        CALL MPI_ISEND(send, n, MPI_DOUBLE_PRECISION, ipartner, tag, comm, Request(nRequest), info)
      ENDIF
    ELSE
      ipartner = irank + 1
      IF (ipartner .LT. PE%nCMFDproc) THEN
        nRequest = nRequest + 1
        CALL MPI_IRECV(recv, n, MPI_DOUBLE_PRECISION, ipartner, tag, comm, Request(nRequest), info)
      ENDIF
    ENDIF
  ENDDO
END SELECT

END SUBROUTINE

SUBROUTINE FinalizeFastComm()
USE PE_MOD,         ONLY : PE
IMPLICIT NONE

INTEGER :: status(MPI_STATUS_SIZE, 4), info

CALL MPI_WAITALL(nRequest, Request, status, info)
CALL MPI_BARRIER(PE%MPI_CMFD_COMM, info)

END SUBROUTINE

SUBROUTINE MatMulMPIFast(Diag, AxOffDiag, x, y, bottomRange, topRange, op)

IMPLICIT NONE

TYPE(CSR_DOUBLE) :: Diag
REAL :: AxOffDiag(*), x(*), y(*)
INTEGER :: bottomRange(2), topRange(2)
CHARACTER(1), OPTIONAL :: op

REAL, POINTER :: recv(:, :), temp(:, :)
REAL, POINTER :: csrVal(:)
INTEGER, POINTER :: csrRowPtr(:), csrColIdx(:)
INTEGER :: nSelf, nNeigh
CHARACTER(1) :: operate

operate = 'n'
IF (PRESENT(op)) operate = op

csrVal => Diag%csrVal
csrRowPtr => Diag%csrRowPtr
csrColIdx => Diag%csrColIdx

nSelf = Diag%nr
nNeigh = bottomRange(2) - bottomRange(1) + 1

IF (mklCntl%lDcpl .OR. .NOT. mklGeom%l3dim) THEN

  CALL mkl_dcsrgemv(operate, nSelf, csrVal, csrRowPtr, csrColIdx, x, y)

ELSE

  ALLOCATE(recv(nNeigh, 2), temp(nNeigh, 2))

  CALL InitFastComm()

  CALL GetNeighborFast(nNeigh, x(bottomRange(1) : bottomRange(2)), recv(:, top), bottom)
  CALL GetNeighborFast(nNeigh, x(topRange(1) : topRange(2)), recv(:, bottom), top)

  CALL mkl_dcsrgemv(operate, nSelf, csrVal, csrRowPtr, csrColIdx, x, y)

  CALL FinalizeFastComm()

  CALL vdmul(nNeigh * 2, AxOffDiag, recv, temp)
  IF (.NOT. mklGeom%lBottom) CALL daxpy(nNeigh, 1.0D0, temp(:, bottom), 1, y(bottomRange(1) : bottomRange(2)), 1)
  IF (.NOT. mklGeom%lTop) CALL daxpy(nNeigh, 1.0D0, temp(:, top), 1, y(topRange(1) : topRange(2)), 1)

  DEALLOCATE(recv, temp)

ENDIF

END SUBROUTINE

!--- MPI Collective Routines

SUBROUTINE normalizeMPI(x, n, comm)

IMPLICIT NONE

INTEGER :: n, comm
REAL :: x(n)
REAL :: norm

norm = normMPI(x, n, comm)
CALL dscal(n, 1.0 / norm, x, 1)

END SUBROUTINE

FUNCTION normMPI(x, n, comm) RESULT(norm)

IMPLICIT NONE

INTEGER :: n, comm
REAL :: x(n)
REAL :: dot, norm

dot = dotMPI(x, x, n, comm)
norm = sqrt(dot)

END FUNCTION

FUNCTION dotMPI(x, y, n, comm) RESULT(buf_dot)

IMPLICIT NONE

INTEGER :: n, comm, ierr
REAL :: x(n), y(n)
REAL :: dot, buf_dot

dot = ddot(n, x, 1, y, 1)

CALL MPI_ALLREDUCE(dot, buf_dot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)

END FUNCTION

END MODULE
#endif