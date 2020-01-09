#include <defines.h>
!--- CNJ Edit : Driver Routine for Gamma CMFD Calculation
#ifdef __GAMMA_TRANSPORT
SUBROUTINE GammaCMFDAcc(CoreInfo, CmInfo, FmInfo)
USE PARAM
USE TYPEDEF,            ONLY : CoreInfo_Type,       CmInfo_Type,        FmInfo_Type,        FxrInfo_Type
USE GamCMFD_MOD,        ONLY : CMFD_gPinXS,         GammaCMFD,          GamHomoXsGen,       GammaBiCGSTAB,          &
                               UpdateCellPhi,       UpdateFsrPhi,       SetGammaCMFDSystem,                         &
                               RadCouplingCoeffGen_Gamma
USE GammaCore_MOD,      ONLY : gPhiAngIn,           gPhis,              gJout,              GamGroupInfo
USE PE_MOD,             ONLY : PE
USE CNTL,               ONLY : nTracerCntl
USE TIMER,              ONLY : nTracer_dclock,      TimeChk
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(FmInfo_Type) :: FmInfo

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
REAL, POINTER :: Phis(:, :, :)
REAL :: CmfdTimeBeg, CmfdTimeEnd
INTEGER :: ig, iter
INTEGER :: myzb, myze
INTEGER :: ng, ngg
INTEGER :: GrpBeg, GrpEnd, nGroupInfo = 2
LOGICAL :: lDhat, lScat1
LOGICAL, SAVE :: lFirst = TRUE

CmfdTImeBeg = nTracer_dclock(FALSE, FALSE)

Fxr => FmInfo%Fxr
Phis => FmInfo%Phis
ng = GamGroupInfo%ng
ngg = GamGroupInfo%ngg
myzb = PE%myzb
myze = PE%myze
lScat1 = nTracerCntl%lGammaScat1
lDhat = .NOT. lFirst
IF (.NOT. GamGroupInfo%lUpScat) nGroupInfo = 1

CALL GamHomoXsGen(CoreInfo, Fxr, Phis, gPhis, CMFD_gPinXS, myzb, myze, ng, ngg, lScat1, FALSE)
CALL RadCouplingCoeffGen_Gamma(CoreInfo, CMFD_gPinXS, gJout, ngg, lDhat, lScat1, PE)

CALL SetGammaCMFDSystem(CoreInfo, GamGroupInfo, CMFD_gPinXS, GammaCMFD, lDhat)
CALL UpdateCellPhi(CoreInfo, CMFD_gPinXS, GamGroupInfo, GammaCMFD)
CALL GammaCMFDProdUpdt(CoreInfo, GamGroupInfo, GammaCMFD)

DO iter = 1, nGroupInfo
  GrpBeg = 1; GrpEnd = ngg
  IF (iter .GT. 1) THEN
    GrpBeg = GamGroupInfo%UpScatRange(1); GrpEnd = GamGroupInfo%UpScatRange(2)
  ENDIF
  DO ig = 1, ngg
    CALL GammaCMFDScatUpdt(CoreInfo, GamGroupInfo, GammaCMFD, ig)
    CALL GammaBiCGSTAB(GammaCMFD%M(ig), GammaCMFD%ILU(ig), GammaCMFD%gPhic(:, :, ig), GammaCMFD%gSrc)
  ENDDO
ENDDO

CALL UpdateFsrPhi(CoreInfo, CmInfo, CMFD_gPinXS, GamGroupInfo, GammaCMFD, gPhis, gPhiAngIn)

IF (lFirst) lFirst = FALSE

CmfdTImeEnd = nTracer_dclock(FALSE, FALSE)
TimeChk%CmfdTime = TimeChk%CmfdTime + (CmfdTimeEnd - CmfdTimeBeg)

END SUBROUTINE

!--- CNJ Edit : Operator Setup Routine for Gamma CMFD Calculation

SUBROUTINE SetGammaCMFDSystem(CoreInfo, GroupInfo, PinXS, GammaCMFD, lDhat)
USE PARAM
USE TYPEDEF,		    ONLY : CoreInfo_Type,       Pin_Type
USE GammaTYPEDEF,       ONLY : gPinXS_Type,         GammaGroupInfo_Type,        GammaCMFD_Type
USE PE_MOD,             ONLY : PE
USE UTILFUNCTION,       ONLY : Array2DSort
USE MKL_BILU,           ONLY : MKL_PrepareILU
USE CSRMATRIX
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(GammaGroupInfo_Type) :: GroupInfo
TYPE(gPinXS_Type), POINTER :: PinXS(:, :)
TYPE(GammaCMFD_Type) :: GammaCMFD
LOGICAL :: lDhat

TYPE(Pin_Type), POINTER :: Pin(:)
INTEGER, PARAMETER :: SELF = 5
INTEGER :: surf(5) = (/ NORTH, WEST, SELF, EAST, SOUTH /)
INTEGER :: ir, ic, iz, ig, igf, igt, ibd, isurf, ipin, ineighpin
INTEGER :: gb, ge
INTEGER :: nxy, nbd, ng, ngg
REAL, POINTER :: PinVol(:, :), hz(:)
REAL, POINTER :: diagVal(:), diagCol(:)
REAL :: Dtil, Dhat
REAL :: val

Pin => CoreInfo%Pin
nxy = CoreInfo%nxy
hz => CoreInfo%hz
PinVol => CoreInfo%PinVol
ng = GroupInfo%ng
ngg = GroupInfo%ngg

iz = 1

!--- Set Group Major Diffusion Operator

!$OMP PARALLEL PRIVATE(diagVal, diagCol, Dtil, Dhat, isurf, ineighpin, ir, ic, val)

ALLOCATE(diagVal(5), diagCol(5))

!$OMP DO SCHEDULE(DYNAMIC)
DO ig = 1, ngg
  CALL createCsr(GammaCMFD%M(ig), 5 * nxy, nxy, nxy)
  DO ipin = 1, nxy
    diagVal = 0.0; diagCol = 0.0
    DO ibd = 1, 4
      Dtil = PinXS(ipin, iz)%Dtil(ibd, ig)
      Dhat = PinXS(ipin, iz)%Dhat(ibd, ig)
      diagVal(ibd) = - (Dtil + Dhat) * hz(iz)
      diagVal(SELF) = diagVal(SELF) + (Dtil - Dhat) * hz(iz)
    ENDDO
    diagVal(SELF) = diagVal(SELF) + PinVol(ipin, iz) * PinXS(ipin, iz)%XSr(ig)
    DO ibd = 1, 5
      isurf = surf(ibd)
      SELECT CASE (isurf)
      CASE (SELF)
        ineighpin = ipin
      CASE (NORTH, WEST, EAST, SOUTH)
        ineighpin = Pin(ipin)%NeighIdx(isurf)
        IF (ineighpin .EQ. VoidCell .OR. ineighpin .EQ. RefCell) diagVal(isurf) = 0.0
      END SELECT
      ir = ipin
      ic = ir + (ineighpin - ipin)
      diagCol(isurf) = ic
    ENDDO
    CALL Array2DSort(diagVal, diagCol, 5, nbd, TRUE, FALSE, 2)
    DO ibd = 1, nbd
      ir = ipin
      ic = diagCol(ibd)
      val = diagVal(ibd)
      CALL pushCsr(GammaCMFD%M(ig), val, ir, ic)
    ENDDO
  ENDDO
  CALL finalizeCsr(GammaCMFD%M(ig), TRUE)
ENDDO
!$OMP END DO

DEALLOCATE(diagVal, diagCol)

!$OMP END PARALLEL

!--- Set Preconditioner

DO ig = 1, ngg
  CALL MKL_PrepareILU(GammaCMFD%M(ig), GammaCMFD%ILU(ig))
ENDDO

!--- Set Source Operator

GammaCMFD%Scat = 0.0

!$OMP PARALLEL PRIVATE(ir, val)
!$OMP DO SCHEDULE(DYNAMIC) COLLAPSE(3)
DO igt = 1, ngg
  DO igf = 1, ngg
    DO ipin = 1, nxy
      ir = ipin
      IF (PinXS(ipin, iz)%XSs(igt)%ib .GT. igf) CYCLE
      IF (PinXS(ipin, iz)%XSs(igt)%ie .LT. igf) CYCLE
      val = PinXS(ipin, iz)%XSs(igt)%from(igf) * PinVol(ipin, iz)
      GammaCMFD%Scat(ir, igf, igt) = val
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

GammaCMFD%Prod = 0.0

!$OMP PARALLEL PRIVATE(ir, val)
!$OMP DO SCHEDULE(DYNAMIC) COLLAPSE(3)
DO igt = 1, ngg
  DO igf = 1, ng
    DO ipin = 1, nxy
      ir = ipin
      val = PinXS(ipin, iz)%XSp(igf, igt) * PinVol(ipin, iz)
      GammaCMFD%Prod(ir, igf, igt) = val
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

!--- CNJ Edit : Source Calculation Routine for Gamma CMFD Calculation

SUBROUTINE GammaCMFDScatUpdt(CoreInfo, GroupInfo, GammaCMFD, ig)
USE PARAM
USE TYPEDEF,            ONLY : CoreInfo_Type
USE GammaTYPEDEF,       ONLY : GammaCMFD_Type,      GammaGroupInfo_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(GammaGroupInfo_Type) :: GroupInfo
TYPE(GammaCMFD_Type) :: GammaCMFD

INTEGER :: ng, nxy
INTEGER :: ig, igf
REAL, ALLOCATABLE :: scatSrc(:)

ng = GroupInfo%ngg
nxy = CoreInfo%nxy

ALLOCATE(scatSrc(nxy))

GammaCMFD%gSrc = 0.0

DO igf = 1, ng
  CALL vdmul(nxy, GammaCMFD%Scat(:, igf, ig), GammaCMFD%gPhic(:, :, igf), scatSrc)
  CALL vdadd(nxy, scatSrc, GammaCMFD%gSrc, GammaCMFD%gSrc)
ENDDO

CALL vdadd(nxy, GammaCMFD%gProd(:, ig), GammaCMFD%gSrc, GammaCMFD%gSrc)

DEALLOCATE(scatSrc)

END SUBROUTINE

SUBROUTINE GammaCMFDProdUpdt(CoreInfo, GroupInfo, GammaCMFD)
USE PARAM
USE TYPEDEF,            ONLY : CoreInfo_Type
USE GammaTYPEDEF,       ONLY : GammaCMFD_Type,      GammaGroupInfo_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(GammaGroupInfo_Type) :: GroupInfo
TYPE(GammaCMFD_Type) :: GammaCMFD

INTEGER :: ng, ngg, nxy
INTEGER :: ig, igf
REAL, ALLOCATABLE :: prodSrc(:)

ng = GroupInfo%ng
ngg = GroupInfo%ngg
nxy = CoreInfo%nxy

ALLOCATE(prodSrc(nxy))

GammaCMFD%gProd = 0.0

DO ig = 1, ngg
  DO igf = 1, ng
    CALL vdmul(nxy, GammaCMFD%Prod(:, igf, ig), GammaCMFD%Phic(:, :, igf), prodSrc)
    CALL vdadd(nxy, prodSrc, GammaCMFD%gProd(:, ig), GammaCMFD%gProd(:, ig))
  ENDDO
ENDDO

DEALLOCATE(prodSrc)

END SUBROUTINE

!--- CNJ Edit : Homogenized Flux Mapping

SUBROUTINE UpdateCellPhi(CoreInfo, PinXS, GroupInfo, GammaCMFD)
USE PARAM
USE TYPEDEF,            ONLY : CoreInfo_Type
USE GammaTYPEDEF,       ONLY : gPinXS_Type,         GammaCMFD_Type,     GammaGroupInfo_Type
USE PE_MOD,             ONLY : PE
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(gPinXS_Type), POINTER :: PinXS(:, :)
TYPE(GammaGroupInfo_Type) :: GroupInfo
TYPE(GammaCMFD_Type) :: GammaCMFD

INTEGER :: iz, ig, ipin
INTEGER :: ng, ngg, nxy
INTEGER :: myzb, myze

ng = GroupInfo%ng
ngg = GroupInfo%ngg
nxy = CoreInfo%nxy
myzb = PE%myzb
myze = PE%myze

DO ig = 1, ng
  DO iz = myzb, myze
    DO ipin = 1, nxy
      GammaCMFD%Phic(ipin, iz, ig) = PinXS(ipin, iz)%nPhi(ig)
    ENDDO
  ENDDO
ENDDO

DO ig = 1, ngg
  DO iz = myzb, myze
    DO ipin = 1, nxy
      GammaCMFD%gPhic(ipin, iz, ig) = PinXS(ipin, iz)%gPhi(ig)
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE

!--- CNJ Edit : Gamma MOC Solution Update

SUBROUTINE UpdateFsrPhi(CoreInfo, CmInfo, PinXS, GroupInfo, GammaCMFD, gPhis, gPhiAngIn)
USE PARAM
USE TYPEDEF,            ONLY : CoreInfo_Type,       CmInfo_Type,        RayInfo_Type,                               &
                               Pin_Type,            Cell_Type
USE GammaTYPEDEF,       ONLY : gPinXS_Type,         GammaCMFD_Type,     GammaGroupInfo_Type
USE PE_MOD,             ONLY : PE
USE CNTL,               ONLY : nTracerCntl
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(gPinXS_Type), POINTER :: PinXS(:, :)
TYPE(GammaGroupInfo_Type) :: GroupInfo
TYPE(GammaCMFD_Type) :: GammaCMFD
REAL, POINTER :: gPhis(:, :, :), gPhiAngIn(:, :, :, :)

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
INTEGER, POINTER :: RotRayInOutCell(:, :)
INTEGER, POINTER :: PhiAngInSvIdx(:, :)
INTEGER :: ng, nxy, nLocalFsr, nRotRay
INTEGER :: myzb, myze
INTEGER :: i, ig, iz, ipin, isv, icel, ifsr, irot
INTEGER :: FsrIdxSt
REAL :: fmult, fmultMg(GroupInfo%ngg)
LOGICAL :: lScat1

Pin => CoreInfo%Pin
Cell => CoreInfo%CellInfo
ng = GroupInfo%ngg
nxy = CoreInfo%nxy
myzb = PE%myzb
myze = PE%myze
nRotRay = CmInfo%RayInfo4CMFD%nRotRay
RotRayInOutCell => CmInfo%RayInfo4CMFD%RotRayInOutCell
PhiAngInSvIdx => CmInfo%RayInfo4CMFD%PhiAngInSvIdx
lScat1 = nTracerCntl%lGammaScat1

IF (lScat1) THEN

  DO ig = 1, ng
    DO iz = myzb, myze
      DO ipin = 1, nxy
        FsrIdxSt = Pin(ipin)%FsrIdxSt
        icel = Pin(ipin)%Cell(iz)
        nLocalFsr = Cell(icel)%nFsr
        fmult = GammaCMFD%gPhic(ipin, iz, ig) / PinXS(ipin, iz)%gPhi(ig)
        DO i = 1, nLocalFsr
          ifsr = FsrIdxSt + i - 1
          gPhis(ifsr, ig, iz) = gPhis(ifsr, ig, iz) * fmult
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  DO irot = 1, 2
    DO i = 1, nRotRay
      ipin = RotRayInOutCell(i, irot)
      IF (ipin .EQ. 0) CYCLE
      isv = PhiAngInSvIdx(i, irot)
      DO iz = myzb, myze
        DO ig = 1, ng
          fmult = GammaCMFD%gPhic(ipin, iz, ig) / PinXS(ipin, iz)%gPhi(ig)
          gPhiAngIn(:, isv, ig, iz) = gPhiAngIn(:, isv, ig, iz) * fmult
        ENDDO
      ENDDO
    ENDDO
  ENDDO

ELSE

  DO iz = myzb, myze
    DO ipin = 1, nxy
      FsrIdxSt = Pin(ipin)%FsrIdxSt
      icel = Pin(ipin)%Cell(iz)
      nLocalFsr = Cell(icel)%nFsr
      fmultMg = GammaCMFD%gPhic(ipin, iz, :) / PinXS(ipin, iz)%gPhi
      DO i = 1, nLocalFsr
        ifsr = FsrIdxSt + i - 1
        gPhis(:, ifsr, iz) = gPhis(:, ifsr, iz) * fmultMg
      ENDDO
    ENDDO
  ENDDO

  DO irot = 1, 2
    DO i = 1, nRotRay
      ipin = RotRayInOutCell(i, irot)
      IF (ipin .EQ. 0) CYCLE
      isv = PhiAngInSvIdx(i, irot)
      DO iz = myzb, myze
        fmultMg = GammaCMFD%gPhic(ipin, iz, :) / PinXS(ipin, iz)%gPhi
        DO ig = 1, ng
          gPhiAngIn(:, ig, isv, iz) = gPhiAngIn(:, ig, isv, iz) * fmultMg(ig)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

ENDIF

END SUBROUTINE

!--- CNJ Edit : Solver Routine for Gamma CMFD Calculation

SUBROUTINE GammaBiCGSTAB(M, ILU, x, b)
USE PARAM
USE CSRMATRIX,      ONLY : CSR_DOUBLE
USE MKL_BILU,       ONLY : MKL_SolveILU
IMPLICIT NONE

TYPE(CSR_DOUBLE) :: M, ILU
REAL :: x(:, :), b(:)

INTEGER :: i, ierr, iter, itermin = 3, itermax = 100
INTEGER :: nr, nc, nnz
REAL :: err, tol = 1.0D-06
REAL :: rho0, rho1, w0, w1, norm0, norm1, dot0, dot1, alpha, beta
REAL, POINTER :: h(:), s(:), t(:), y(:), z(:)
REAL, POINTER :: r0(:), r1(:), rhat(:), v0(:), v1(:), p0(:), p1(:)
REAL, POINTER :: csrVal(:)
INTEGER, POINTER :: csrRowPtr(:), csrColIdx(:)
LOGICAL :: lConv

INCLUDE 'mkl_blas.fi'

nr = M%nr
nc = M%nc
nnz = M%nnz

csrVal => M%csrVal
csrRowPtr => M%csrRowPtr
csrColIdx => M%csrColIdx

! =========================================================================== !
!  Intel MKL Preconditioned Bi-Conjugate Gradient Stablized(BiCGSTAB) Solver  !
!     http://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method     !
! =========================================================================== !

ALLOCATE(h(nc), s(nc), t(nc), y(nc), z(nc))
ALLOCATE(r0(nc), r1(nc), rhat(nc))
ALLOCATE(v0(nc), v1(nc))
ALLOCATE(p0(nc), p1(nc))

!--- 1 ----------------------
CALL mkl_dcsrgemv('n', nr, csrVal, csrRowPtr, csrColIdx, x, r1)
CALL vdsub(nc, b, r1, r1)
!--- 2 ----------------------
CALL dcopy(nc, r1, 1, rhat, 1)
!--- 3 ----------------------
norm0 = dnrm2(nc, r1, 1)
rho1 = 1.0; w1 = 1.0; alpha = 1.0
!--- 4 ----------------------
v1 = 0.0; p1 = 0.0
!--- 5 ----------------------
lConv = FALSE; iter = 0
DO WHILE (.NOT. lConv)
  !--- 5.0 ------------------
  iter = iter + 1
  rho0 = rho1; w0 = w1
  CALL dcopy(nc, r1, 1, r0, 1)
  CALL dcopy(nc, p1, 1, p0, 1)
  CALL dcopy(nc, v1, 1, v0, 1)
  !--- 5.1 ------------------
  rho1 = ddot(nc, rhat, 1, r0, 1)
  !--- 5.2 ------------------
  beta = (rho1 / rho0) * (alpha / w0)
  !--- 5.3 ------------------
  CALL dcopy(nc, r0, 1, p1, 1)
  CALL daxpy(nc, -w0, v0, 1, p0, 1)
  CALL daxpy(nc, beta, p0, 1, p1, 1)
  !--- 5.4 ------------------
  CALL MKL_SolveILU(ILU, y, p1)
  !--- 5.5 ------------------
  CALL mkl_dcsrgemv('n', nr, csrVal, csrRowPtr, csrColIdx, y, v1)
  !--- 5.6 ------------------
  dot1 = ddot(nc, rhat, 1, v1, 1)
  alpha = rho1 / dot1
  !--- 5.7 ------------------
  CALL dcopy(nc, x, 1, h, 1)
  CALL daxpy(nc, alpha, y, 1, h, 1)
  !--- 5.9 ------------------
  CALL dcopy(nc, r0, 1, s, 1)
  CALL daxpy(nc, -alpha, v1, 1, s, 1)
  !--- 5.10 -----------------
  CALL MKL_SolveILU(ILU, z, s)
  !--- 5.11 -----------------
  CALL mkl_dcsrgemv('n', nr, csrVal, csrRowPtr, csrColIdx, z, t)
  !--- 5.12 -----------------
  dot0 = ddot(nc, t, 1, s, 1)
  dot1 = ddot(nc, t, 1, t, 1)
  w1 = dot0 / dot1
  !--- 5.13 -----------------
  CALL dcopy(nc, h, 1, x, 1)
  CALL daxpy(nc, w1, z, 1, x, 1)
  !--- 5.15 -----------------
  CALL dcopy(nc, s, 1, r1, 1)
  CALL daxpy(nc, -w1, t, 1, r1, 1)
  !--- Convergence ----------
  norm1 = dnrm2(nc, r1, 1)
  err = norm1 / norm0; lConv = (err .LE. tol) .AND. (iter .GE. itermin)
  lConv = lConv .OR. (iter .GE. itermax)
ENDDO

DEALLOCATE(h, s, t, y, z)
DEALLOCATE(r0, r1, rhat)
DEALLOCATE(v0, v1)
DEALLOCATE(p0, p1)

END SUBROUTINE
#endif
