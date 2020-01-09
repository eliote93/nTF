#include <defines.h>
!--- CNJ Edit : BILU Preconditioner Module for 3D CMFD Acceleration with Intel MKL
#ifdef __INTEL_MKL

MODULE MKL_BILU

USE MKL_3D
IMPLICIT NONE

CONTAINS

!--- ILU Routines ---------------------------------------------------------------------------------

SUBROUTINE MKL_PrepareILU(M, ILU)

IMPLICIT NONE

TYPE(CSR_DOUBLE) :: M, ILU

REAL, POINTER :: csrVal(:)
INTEGER, POINTER :: csrRowPtr(:), csrColIdx(:)
INTEGER :: nr, nc
INTEGER :: ierr, ipar(128)
REAL :: dpar(128)

nr = M%nr
nc = M%nc
ipar = 0; ipar(31) = -1
dpar(31 : 32) = 1.0D-08

csrVal => M%csrVal
csrRowPtr => M%csrRowPtr
csrColIdx => M%csrColIdx
ILU = M

CALL dcsrilu0(nr, csrVal, csrRowPtr, csrColIdx, ILU%csrVal, ipar, dpar, ierr)

END SUBROUTINE

SUBROUTINE MKL_SolveILU(ILU, x, b)

IMPLICIT NONE

TYPE(CSR_DOUBLE) :: ILU
REAL, POINTER :: x(:), y(:), b(:)

REAL, POINTER :: csrVal(:)
INTEGER, POINTER :: csrRowPtr(:), csrColIdx(:)
INTEGER :: nr, nc

csrVal => ILU%csrVal
csrRowPtr => ILU%csrRowPtr
csrColIdx => ILU%csrColIdx
nr = ILU%nr
nc = ILU%nc

ALLOCATE(y(nr))

CALL mkl_dcsrtrsv('L', 'n', 'u', nr, csrVal, csrRowPtr, csrColIdx, b, y)
CALL mkl_dcsrtrsv('U', 'n', 'n', nr, csrVal, csrRowPtr, csrColIdx, y, x)

DEALLOCATE(y)

END SUBROUTINE

!--- SPAI Routines --------------------------------------------------------------------------------

SUBROUTINE MKL_PrepareSPAI(M, SPAI)

IMPLICIT NONE 

TYPE(CSR_DOUBLE) :: M, SPAI

TYPE(CSR_DOUBLE) :: trM, trSPAI
TYPE(CSR_DOUBLE) :: sub_trM
REAL, POINTER :: subM(:, :), tau(:), work(:), ek(:, :), Pk(:)
INTEGER, POINTER :: sub_colIdx(:), sub_csrRowIdx(:), sub_rowIdx(:)

REAL :: tmpVal
INTEGER :: tmpInt
INTEGER :: lwork, info, lda

INTEGER :: ntau, kIdx
INTEGER :: nc, nr, nnz 
INTEGER :: sub_nc, sub_nr, max_subr, sub_nnz
INTEGER :: subr, subc
INTEGER :: icol, tr_col
INTEGER :: offset
INTEGER :: i, j, k

nnz = M%nnz
nr = M%nr
nc = M%nc

CALL transposeCsr(M, trM)

CALL createCsr(trSPAI, nnz, nr, nc)
DO i = 1, nc
  sub_nc = trM%csrRowPtr(i + 1) - trM%csrRowPtr(i)
  ALLOCATE(sub_colIdx(sub_nc))
  offset = trM%csrRowPtr(i) - 1
  DO j = trM%csrRowPtr(i), trM%csrRowPtr(i + 1) - 1
    sub_colIdx(j - offset) = trM%csrColIdx(j)
  END DO             
  
  max_subr = 0
  DO j = 1, sub_nc
    icol = sub_colIdx(j)
    max_subr = max_subr + trM%csrRowPtr(icol + 1) - trM%csrRowPtr(icol)
  END DO
    
  !--- Save Transposed Submatrix in CSR format
  CALL createCsr(sub_trM, max_subr, sub_nc, nc)
  DO j = 1, sub_nc
    icol = sub_colIdx(j)
    DO k = trM%csrRowPtr(icol), trM%csrRowPtr(icol + 1) - 1
      tr_col = trM%csrColIdx(k)
      tmpVal = trM%csrVal(k)
      CALL pushCsr(sub_trM, tmpVal, j, tr_col)
    END DO 
  END DO
  CALL finalizeCsr(sub_trM, .FALSE.)
    
  !--- Sorting
  sub_nnz = sub_trM%nnz
  ALLOCATE(sub_csrRowIdx(sub_nnz), sub_rowIdx(sub_nnz))
  DO j = 1, sub_nc
    DO k = sub_trM%csrRowPtr(j), sub_trM%csrRowPtr(j + 1) - 1
      sub_csrRowIdx(k) = j
    END DO
  END DO
    
  DO j = sub_nnz, 1, -1
    DO k = 1, j-1
      IF(sub_trM%csrColIdx(k) .GT. sub_trM%csrColIdx(k + 1)) THEN
        tmpInt = sub_trM%csrColIdx(k)
        sub_trM%csrColIdx(k) = sub_trM%csrColIdx(k + 1)
        sub_trM%csrColIdx(k+1) = tmpInt
              
        tmpInt = sub_csrRowIdx(k)
        sub_csrRowIdx(k) = sub_csrRowIdx(k + 1)
        sub_csrRowIdx(k + 1) = tmpInt
                
        tmpVal = sub_trM%csrVal(k)
        sub_trM%csrVal(k) = sub_trM%csrVal(k + 1)
        sub_trM%csrVal(k + 1) = tmpVal
      END IF
    END DO              
  END DO 
    
  sub_rowIdx(1) = 1
  IF(sub_trM%csrColIdx(1) .EQ. i) kIdx =  1
  DO j = 2, sub_nnz
    IF(sub_trM%csrColIdx(j) .EQ. sub_trM%csrColIdx(j - 1)) THEN 
      sub_rowIdx(j) = sub_rowIdx(j - 1)
    ELSE
      sub_rowIdx(j) = sub_rowIdx(j - 1) + 1
    END IF
    IF(sub_trM%csrColIdx(j) .EQ. i) kIdx =  sub_rowIdx(j)
  END DO
  sub_nr = sub_rowIdx(sub_nnz)        
  
  ntau = min(sub_nr, sub_nc)
      
  ALLOCATE(subM(sub_nr, sub_nc), ek(sub_nr, 1), tau(ntau), work(sub_nc), Pk(sub_nc))
  subM = 0.
  ek = 0.
  ek(kIdx, 1) = 1.
      
  DO j = 1, sub_nnz
    subr = sub_rowIdx(j)
    subc = sub_csrRowIdx(j)
    subM(subr, subc) = sub_trM%csrVal(j)             
  END DO 
      
  CALL dgeqrf(sub_nr, sub_nc, subM, sub_nr, tau, work, sub_nc, info)
  CALL dormqr('L', 't', sub_nr, 1, ntau, subM, sub_nr, tau, ek, sub_nr, work, sub_nc, info)
  
  DO j = sub_nc, 1, -1
    tmpVal = ek(j, 1)
    DO k = sub_nc, j + 1, -1
      tmpVal = tmpVal - subM(j, k) * Pk(k)
    END DO 
    Pk(j) = tmpVal / subM(j, j)                        
  END DO 
      
  DO j = 1, sub_nc
    icol = sub_colIdx(j)
    CALL pushCsr(trSPAI, Pk(j), i, icol)            
  END DO 
             
  CALL destroyCsr(sub_trM)
  DEALLOCATE(subM, ek, tau, work, Pk)
  DEALLOCATE(sub_colIdx, sub_csrRowIdx, sub_rowIdx)
      
END DO 

CALL finalizeCsr(trSPAI, .FALSE.)
CALL transposeCsr(trSPAI, SPAI)
CALL finalizeCsr(SPAI, .TRUE.)
CALL destroyCsr(trSPAI)
CALL destroyCsr(trM)

END SUBROUTINE 

SUBROUTINE MKL_ApplySPAI(SPAI, x, y)

IMPLICIT NONE

TYPE(CSR_DOUBLE) :: SPAI
REAL, POINTER :: x(:), y(:)

REAL, POINTER :: csrVal(:)
INTEGER, POINTER :: csrRowPtr(:), csrColIdx(:)
INTEGER :: nr

csrVal => SPAI%csrVal
csrRowPtr => SPAI%csrRowPtr
csrColIdx => SPAI%csrColIdx

nr = SPAI%nr

CALL mkl_dcsrgemv('n', nr, csrVal, csrRowPtr, csrColIdx, x, y)

END SUBROUTINE

!--- BILU Setup Routines --------------------------------------------------------------------------

! SUBROUTINE MKL_PrepareBILU(CMFD, PinXS, l3dim)
! USE TYPEDEF,        ONLY : PinXS_Type
! IMPLICIT NONE
! 
! TYPE(mklCMFD_Type) :: CMFD
! TYPE(PinXS_Type), POINTER :: PinXS(:, :)
! LOGICAL :: l3dim
! 
! INTEGER :: ig, iz
! INTEGER :: ng, nzCMFD
! 
! ng = CMFD%ng
! nzCMFD = mklGeom%nzCMFD
! 
! CALL SetPackedCMFDSystem(CMFD, PinXS, l3dim)
! 
! DO ig = 1, ng
!   DO iz = 1, nzCMFD
!     CALL LUFactorize2D(CMFD%BILU(ig)%blockMat(iz))
!   ENDDO
! ENDDO
! 
! END SUBROUTINE
! 
! SUBROUTINE SetPackedCMFDSystem(CMFD, PinXS, l3dim)
! USE TYPEDEF,        ONLY : PinXS_Type
! IMPLICIT NONE
! 
! TYPE(mklCMFD_Type) :: CMFD
! TYPE(PinXS_Type), POINTER :: PinXS(:, :)
! LOGICAL :: l3dim
! 
! TYPE(mklBILU_Type), POINTER :: BILU
! TYPE(blockMat_Type), POINTER :: blockMat
! TYPE(triLU_Type), POINTER :: blockDiag, Del
! TYPE(superPin_Type), POINTER :: Pin(:)
! REAL, POINTER :: hzfm(:), PinVolFm(:, :)
! INTEGER, POINTER :: pinMap(:), pinMapRev(:), planeMap(:)
! INTEGER :: ng, nx, ny, nxy, nzCMFD
! INTEGER :: ig, idx, ibd, ix, iy, iz, izf, ipin, ipin_map, ineighpin
! REAL :: Dtil, Dhat
! REAL :: diag, offdiag(6)
! 
! ng = CMFD%ng
! ny = mklGeom%ny
! nxy = mklGeom%nxy
! nzCMFD = mklGeom%nzCMFD
! Pin => mklGeom%superPin
! hzfm => mklGeom%hzfm
! PinVolFm => mklGeom%PinVolFm
! pinMap => mklGeom%pinMap
! pinMapRev => mklGeom%pinMapRev
! planeMap => mklGeom%planeMap
! 
! !$OMP PARALLEL PRIVATE(BILU, blockMat, blockDiag, Del, nx, idx, iz, ipin, ipin_map, ineighpin, Dtil, Dhat, diag, offdiag)
! !$OMP DO SCHEDULE(DYNAMIC)
! DO ig = 1, ng
!   BILU => CMFD%BILU(ig)
!   ipin = 0
!   DO izf = 1, nzCMFD
!     blockMat => BILU%blockMat(izf)
!     iz = planeMap(izf)
!     DO iy = 1, ny
!       blockDiag => blockMat%blockDiag(iy)
!       Del => blockMat%Del(iy)
!       nx = mklGeom%nx(iy)
!       DO ix = 1, nx
!         ipin = ipin + 1
!         ipin_map = pinMap(ipin)
!         diag = 0
!         DO ibd = 1, 4
!           Dtil = PinXS(ipin_map, iz)%Dtil(ibd, ig)
!           Dhat = PinXS(ipin_map, iz)%Dhat(ibd, ig)
!           diag = diag + (Dtil - Dhat) * hzfm(izf)
!           offdiag(ibd) = - (Dtil + Dhat) * hzfm(izf)
!         ENDDO
!         IF (l3dim) THEN
!           DO ibd = 5, 6
!             Dtil = CMFD%AxDtil(ibd - 4, ig, ipin, izf)
!             Dhat = CMFD%AxDhat(ibd - 4, ig, ipin, izf)
!             offdiag(ibd) = - (Dtil + Dhat) * PinVolFm(ipin, izf) / hzfm(izf)
! 	        diag = diag + (Dtil - Dhat) * PinVolFm(ipin, izf) / hzfm(izf)
!           ENDDO
!         ENDIF
!         diag = diag + PinVolFm(ipin, izf) * PinXS(ipin_map, iz)%XSr(ig)
!         blockDiag%diag(ix) = diag
!         IF (ix .NE. 1) blockDiag%lower(ix) = offdiag(WEST)
!         IF (ix .NE. nx) blockDiag%upper(ix) = offdiag(EAST)
!         IF (iy .NE. 1) blockMat%lower(ix, iy) = offdiag(NORTH)
!         IF (iy .NE. ny) blockMat%upper(ix, iy) = offdiag(SOUTH)
!         IF (l3dim) THEN
!           idx = (iy - 1) * nx + ix
!           IF (iz .NE. 1) BILU%lower(idx, izf) = offdiag(BOTTOM + 4)
!           IF (iz .NE. nzCMFD) BILU%upper(idx, izf) = offdiag(TOP + 4)
!         ENDIF
!       ENDDO
!       CALL dcopy(nx - 1, blockDiag%lower, 1, Del%lower, 1)
!       CALL dcopy(nx, blockDiag%diag, 1, Del%diag, 1)
!       CALL dcopy(nx - 1, blockDiag%upper, 1, Del%upper, 1)
!     ENDDO
!   ENDDO
! ENDDO
! !$OMP END DO
! !$OMP END PARALLEL
! 
! END SUBROUTINE
! 
! SUBROUTINE LUFactorize2D(blockMat)
! 
! IMPLICIT NONE
! 
! TYPE(blockMat_Type) :: blockMat
! TYPE(triLU_Type), POINTER :: Del
! INTEGER :: nx, ny, iy, info
! 
! ny = blockMat%ny
! 
! DO iy = 2, ny
!   CALL Recursion(blockMat, iy)
! ENDDO
! 
! DO iy = 1, ny
!   Del => blockMat%Del(iy)
!   nx = mklGeom%nx(iy)
! !  CALL ddttrf(nx, Del%lower, Del%diag, Del%upper, info)
! ENDDO
! 
! END SUBROUTINE
! 
! SUBROUTINE Recursion(blockMat, iy)
! 
! IMPLICIT NONE
! 
! TYPE(blockMat_Type) :: blockMat
! TYPE(triLU_Type), POINTER :: blockDiag(:), Del(:)
! INTEGER :: nx, iy
! 
! blockDiag => blockMat%blockDiag
! Del => blockMat%Del
! nx = blockDiag(iy)%nx
! 
! CALL ApproxLUInv(Del(iy - 1), Del(iy), nx)
! 
! CALL vdmul(nx - 1, blockMat%lower(2 : nx, iy), Del(iy)%lower, Del(iy)%lower)
! CALL vdmul(nx, blockMat%lower(1 : nx, iy), Del(iy)%diag, Del(iy)%diag)
! CALL vdmul(nx - 1, blockMat%lower(1 : nx - 1, iy), Del(iy)%upper, Del(iy)%upper)
! 
! CALL vdmul(nx - 1, blockMat%upper(2 : nx, iy - 1), Del(iy)%lower, Del(iy)%lower)
! CALL vdmul(nx, blockMat%upper(1 : nx, iy - 1), Del(iy)%diag, Del(iy)%diag)
! CALL vdmul(nx - 1, blockMat%upper(1 : nx - 1, iy - 1), Del(iy)%upper, Del(iy)%upper)
! 
! CALL vdsub(nx - 1, blockDiag(iy)%lower, Del(iy)%lower, Del(iy)%lower)
! CALL vdsub(nx, blockDiag(iy)%diag, Del(iy)%diag, Del(iy)%diag)
! CALL vdsub(nx - 1, blockDiag(iy)%upper, Del(iy)%upper, Del(iy)%upper)
! 
! END SUBROUTINE
! 
! SUBROUTINE ApproxLUInv(LU, LUinv, n)
! 
! IMPLICIT NONE
! 
! TYPE(triLU_Type) :: LU, LUinv
! INTEGER :: i, j, n
! REAL :: phi(1 : n + 1), theta(0 : n)
! 
! ! =========================================================== !
! !     Sparsity Preserving Inversion of Tridiagonal Matrix     !
! !  http://en.wikipedia.org/wiki/Tridiagonal_matrix#Inversion  !
! ! =========================================================== !
! 
! !--- Recurrence Relations
! 
! theta(0) = 1.0; theta(1) = LU%diag(1)
! phi(n + 1) = 1.0; phi(n) = LU%diag(n)
! 
! DO i = 2, n
!   theta(i) = LU%diag(i) * theta(i - 1) - LU%upper(i - 1) * LU%lower(i) * theta(i - 2)
! ENDDO
! 
! DO i = n - 1, 1, -1
!   phi(i) = LU%diag(i) * phi(i + 1) - LU%upper(i) * LU%lower(i + 1) * phi(i + 2)
! ENDDO
! 
! !--- Inversion
! 
! DO i = 1, n
!   j = i
!   LUinv%diag(i) = theta(i - 1) * phi(j + 1) / theta(n)
!   j = i - 1
!   IF (i .NE. 1) LUinv%lower(i) = - LU%lower(i) * theta(j - 1) * phi(i + 1) / theta(n)
!   j = i + 1
!   IF (i .NE. n) LUinv%upper(i) = - LU%upper(j - 1) * theta(i - 1) * phi(j + 1) / theta(n) 
! ENDDO
! 
! END SUBROUTINE
! 
! !--- BILU Solver Routines -------------------------------------------------------------------------
! 
! SUBROUTINE MKL_SolveBILU(BILU, p, s)
! 
! IMPLICIT NONE
! 
! TYPE(mklBILU_Type) :: BILU
! REAL, POINTER :: p(:), s(:)
! 
! INTEGER :: iz, ib, ie
! INTEGER :: nxy, nzCMFD
! REAL, ALLOCATABLE :: x(:, :), y(:, :), xsi(:, :), b(:)
! 
! nxy = mklGeom%nxy
! nzCMFD = mklGeom%nzCMFD
! 
! ALLOCATE(x(nxy, nzCMFD), y(nxy, nzCMFD), xsi(nxy, nzCMFD), b(nxy))
! 
! CALL dcopy(nxy, s, 1, b, 1)
! CALL ForwardSolve2D(BILU%blockMat(1), y(:, 1), b)
! CALL dcopy(nxy, y(:, 1), 1, b, 1)
! CALL BackwardSolve2D(BILU%blockMat(1), y(:, 1), b)
! 
! DO iz = 2, nzCMFD
!   CALL vdmul(nxy, BILU%lower(:, iz), y(:, iz - 1), b)
!   ib = (iz - 1) * nxy + 1; ie = iz * nxy
!   CALL vdsub(nxy, s(ib : ie), b, b)
!   CALL ForwardSolve2D(BILU%blockMat(iz), y(:, iz), b)
!   CALL dcopy(nxy, y(:, iz), 1, b, 1)
!   CALL BackwardSolve2D(BILU%blockMat(iz), y(:, iz), b)
! ENDDO
! 
! CALL dcopy(nxy, y(:, nzCMFD), 1, x(:, nzCMFD), 1)
! 
! DO iz = nzCMFD - 1, 1, -1
!   CALL vdmul(nxy, BILU%upper(:, iz), x(:, iz + 1), b)
!   CALL ForwardSolve2D(BILU%blockMat(iz), xsi(:, iz), b)
!   CALL dcopy(nxy, xsi(:, iz), 1, b, 1)
!   CALL BackwardSolve2D(BILU%blockMat(iz), xsi(:, iz), b)
!   CALL vdsub(nxy, y(:, iz), xsi(:, iz), x(:, iz))
! ENDDO
! 
! CALL dcopy(nxy * nzCMFD, x, 1, p, 1)
! 
! DEALLOCATE(x, y, xsi, b)
! 
! END SUBROUTINE
! 
! SUBROUTINE ForwardSolve2D(blockMat, x, b)
! 
! IMPLICIT NONE
! 
! TYPE(blockMat_Type) :: blockMat
! REAL :: x(:), b(:)
! 
! INTEGER :: iy, nx, npx, ny, xb, xe, xpb, xpe, ixb, ixe
! REAL :: bhat(mklGeom%nxmax)
! 
! ny = blockMat%ny
! 
! nx = mklGeom%nx(1)
! xb = mklGeom%pinRange(1, 1)
! xe = mklGeom%pinRange(2, 1)
! CALL dcopy(nx, b(xb : xe), 1, bhat, 1)
! CALL LUSolve1D(blockMat%Del(1), x(xb : xe), bhat)
! 
! DO iy = 2, ny
!   nx = mklGeom%nx(iy); npx = mklGeom%nx(iy - 1)
!   xb = mklGeom%pinRange(1, iy); xe = mklGeom%pinRange(2, iy)
!   xpb = mklGeom%pinRange(1, iy - 1); xpe = mklGeom%pinRange(2, iy - 1)
!   ixb = mklGeom%ixRange(1, iy - 1); ixe = mklGeom%ixRange(2, iy - 1)
!   CALL vdmul(npx, blockMat%lower(ixb : ixe, iy), x(xpb : xpe), bhat)
!   CALL vdsub(nx, b(xb : xe), bhat, bhat)
!   CALL LUSolve1D(blockMat%Del(iy), x(xb : xe), bhat)
! ENDDO
! 
! END SUBROUTINE
! 
! SUBROUTINE BackwardSolve2D(blockMat, x, y)
! 
! IMPLICIT NONE
! 
! TYPE(blockMat_Type) :: blockMat
! REAL :: x(:), y(:)
! 
! INTEGER :: iy, nx, npx, ny, xb, xe, xpb, xpe, ixb, ixe
! REAL :: xsi(mklGeom%nxmax), bhat(mklGeom%nxmax)
! 
! ny = blockMat%ny
! 
! nx = mklGeom%nx(ny)
! xb = mklGeom%pinRange(1, ny)
! xe = mklGeom%pinRange(2, ny)
! CALL dcopy(nx, y(xb : xe), 1, x(xb : xe), 1)
! 
! DO iy = ny - 1, 1, -1
!   nx = mklGeom%nx(iy); npx = mklGeom%nx(iy + 1)
!   xb = mklGeom%pinRange(1, iy); xe = mklGeom%pinRange(2, iy)
!   xpb = mklGeom%pinRange(1, iy + 1); xpe = mklGeom%pinRange(2, iy + 1)
!   ixb = mklGeom%ixRange(1, iy + 1); ixe = mklGeom%ixRange(2, iy + 1)
!   CALL vdmul(npx, blockMat%upper(ixb : ixe, iy), x(xpb : xpe), bhat)
!   CALL LUSolve1D(blockMat%Del(iy), xsi, bhat)
!   CALL vdsub(nx, y(xb : xe), xsi, x(xb : xe))
! ENDDO
! 
! END SUBROUTINE
! 
! SUBROUTINE LUSolve1D(triLU, x, b)
! 
! IMPLICIT NONE
! 
! TYPE(triLU_Type) :: triLU
! REAL :: x(:), b(:)
! INTEGER :: nx, info
! 
! nx = triLU%nx
! 
! CALL dcopy(nx, b, 1, x, 1)
! ! CALL ddttrsv('L', 'n', nx, 1, triLU%lower, triLU%diag, triLU%upper, x, nx, info)
! ! CALL ddttrsv('U', 'n', nx, 1, triLU%lower, triLU%diag, triLU%upper, x, nx, info)
! 
! END SUBROUTINE

END MODULE
    
#endif