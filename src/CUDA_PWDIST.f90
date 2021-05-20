#include <CUDADEFINES.h>
#include <defines.h>

#ifdef __PGI
MODULE CUDA_PWDIST
USE CUDA_MASTER
USE CUDA_UTIL
IMPLICIT NONE
TYPE(CSR_DOUBLE), PRIVATE :: CMConverter
REAL, PRIVATE :: inv_totvolf
LOGICAL, PRIVATE :: linitfin = .FALSE.

CONTAINS
SUBROUTINE cuUpdtPinPw(PwDist, lDcyHeat, lfpowersave, AvgPwOut)
USE TranMacXsLib_Mod,     ONLY : Alpha, Zeta
IMPLICIT NONE
TYPE(cuPwDist_Type) :: PwDist
LOGICAL :: lDcyHeat, lfpowersave
REAL, OPTIONAL :: AvgPwOut

TYPE(CSR_DOUBLE) :: kF
REAL :: locVol, PwSum, AvgPw, fnorm
REAL :: fpcoeff
INTEGER :: nxy, nzCMFD
INTEGER :: iprec
INTEGER :: ierr, n

nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
n = nxy * nzCMFD

CALL cuSetFissionHeatOperator(kF, cuCMFD)
IF(cuDevice%lFuel) THEN
  ierr = cusparseDcsrmv(cuDevice%mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, kF%nr, kF%nc, kF%nnz, 1.0_8,   &
                        kF%descr, kF%d_csrVal, kF%d_csrRowPtr, kF%d_csrColIdx, cuCMFD%phis8, 0.0_8, PwDist%pinPw)
  CALL destroyCsr(kF)
ELSE
  PwDist%pinPw = 0.
END IF
IF(lfpowersave) THEN
  ierr = cublasDcopy_v2(cuDevice%myblasHandle, n, PwDist%pinPw, 1, PwDist%pinPwf, 1)
END IF
IF(lDcyHeat) THEN
  fpcoeff = 1._8 - sum(Alpha(1:6))
  ierr = cublasDscal_v2(cuDevice%myblasHandle, n, fpcoeff, PwDist%pinPw, 1)
  DO iprec = 1, 6
    ierr = cublasDaxpy_v2(cuDevice%myblasHandle, n, Zeta(iprec), PwDist%hPrec(:,iprec), 1, PwDist%pinPw, 1)
  END DO
END IF

PwSum = asumMulti(PwDist%pinPw, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
AvgPw = PwSum * inv_totvolf
fnorm = 1._8 / AvgPw
ierr = cublasDscal_v2(cuDevice%myblasHandle, n, fnorm, PwDist%pinPw, 1)
CALL cuVectorOp('*', n, PwDist%pinPw, PwDist%invPinVolFm, PwDist%pinPw, cuDevice%myStream)
IF(PRESENT(AvgPwOut)) AvgPwOut = AvgPw

END SUBROUTINE

SUBROUTINE cuSetFissionHeatOperator(kF, CMFD)
IMPLICIT NONE
TYPE(CSR_DOUBLE) :: kF
TYPE(cuCMFD_Type) :: CMFD

REAL, POINTER :: PinVolFm(:,:)
INTEGER, POINTER :: planeMap(:), pinMap(:)
REAL :: locVol, val
INTEGER :: myzbf, myzef, ng, nxy, nzCMFD
INTEGER :: izf, ixy, ixy_map, iz, ig, ir, ic

IF(.NOT. cuDevice%lFuel) RETURN

ng = CMFD%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef

PinVolFm => cuGeometry%PinVolFmF
pinMap => cuGeometry%pinMap
planeMap => cuGeometry%planeMap

CALL createCsr(kF, ng*nxy*nzCMFD, nxy*nzCMFD, ng*nxy*nzCMFD)
DO izf = myzbf, myzef
  iz = planeMap(izf)
  DO ixy = 1, nxy
    ixy_map = pinMap(ixy)
    ir = ixy + (izf - myzbf) * nxy
    locVol = PinVolFm(ixy_map, izf)
    DO ig = 1, ng
      ic = ig + (ixy - 1) * ng + (izf - myzbf) * ng * nxy
      val = CMFD%PinXS(ixy_map, iz)%XSKF(ig) * locVol
      CALL pushCsr(kF, val, ir, ic)
    END DO
  END DO
END DO
CALL finalizeCsr(kF, .TRUE.)
NULLIFY(PlaneMap, PinMap, PinVolFm)

END SUBROUTINE

SUBROUTINE cuInithPrec(PWDist)
USE TranMacXsLib_Mod,     ONLY : Alpha, Zeta
IMPLICIT NONE
TYPE(cuPwDist_Type) :: PwDist

REAL :: alpinvzet
INTEGER :: iprec
INTEGER :: ierr, n

n = cuGeometry%nxyc * cuDevice%nzCMFD
PwDist%hPrec = 0.
DO iprec = 1, 6
  alpinvzet = Alpha(iprec) / Zeta(iprec)
  ierr = cublasDaxpy_v2(cuDevice%myblasHandle, n, alpinvzet, PwDist%pinPwf, 1, PwDist%hPrec(:,iprec), 1)
END DO

END SUBROUTINE

SUBROUTINE cuUpdthPrec(PWDist, delt)
USE TranMacXsLib_Mod, ONLY : Alpha, Zeta
IMPLICIT NONE
TYPE(cuPwDist_Type) :: PwDist
REAL :: delt

REAL :: kappa, fpcoeff
INTEGER :: iprec
INTEGER :: ierr, n

n = cuGeometry%nxyc * cuDevice%nzCMFD
DO iprec = 1, 6
  kappa = exp(-delt * Zeta(iprec))
  ierr = cublasDscal_v2(cuDevice%myblasHandle, n, kappa, PWDist%hPrec(:,iprec), 1)
  fpcoeff = (1._8-kappa)*Alpha(iprec)/Zeta(iprec)
  ierr = cublasDaxpy_v2(cuDevice%myblasHandle, n, fpcoeff, PwDist%pinPwf, 1, PWDist%hPrec(:,iprec), 1 )
END DO

END SUBROUTINE

SUBROUTINE cuSetCMRelPower(PWDist, RelPower, PE, corenz, lGather)
USE TYPEDEF,    ONLY : PE_TYPE
IMPLICIT NONE
TYPE(cuPwDist_Type) :: PwDist
REAL, POINTER :: RelPower(:, :)
TYPE(PE_TYPE) :: PE
INTEGER :: corenz
LOGICAL :: lGather

REAL, ALLOCATABLE, DEVICE :: dRelPower(:)
REAL, ALLOCATABLE :: hRelPower(:)
INTEGER :: nxy, myzb, myze
INTEGER :: nr, nc, nnz
INTEGER :: nb, ne, nproc
INTEGER :: i, ixy, iz, cnt
INTEGER :: ierr, n

nxy = cuGeometry%nxy
myzb = cuDevice%myzb
myze = cuDevice%myze

nr = CMConverter%nr
nc = CMConverter%nc
nnz = CMConverter%nnz

n = nxy * corenz
nb = (myzb-1) * nxy + 1
ne = (myze) * nxy

ALLOCATE(hRelPower(n))
IF(cuDevice%lFuel) THEN
  ALLOCATE(dRelPower(nr))
  ierr = cusparseDcsrmv(cuDevice%mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, 1.0_8,   &
                        CMConverter%descr, CMConverter%d_csrVal, CMConverter%d_csrRowPtr, CMConverter%d_csrColIdx,&
                        PwDist%pinPw, 0.0_8, dRelPower)
  ierr = cudaMemcpy(hRelPower(nb:ne), dRelPower, nr, cudaMemcpyDeviceToHost)
  DEALLOCATE(dRelPower)
ELSE
  hRelPower(nb:ne) = 0.
END IF
IF(lGather) THEN
  nproc = PE%nCmfdProc
  DO i = 0, nproc-1
    nb = (PE%AxDomRange(1,i)-1) * nxy + 1
    ne = (PE%AxDomRange(2,i)) * nxy
    CALL MPI_BCAST(hRelPower(nb:ne), ne-nb+1, MPI_DOUBLE_PRECISION, i, PE%MPI_CMFD_COMM, ierr)
  END DO
  cnt = 0
ELSE
  nproc = 1
  cnt = nb-1
END IF
DO i = 0, nproc-1
  DO ixy = 1, nxy
    DO iz = PE%AxDomRange(1, i), PE%AxDomRange(2,i)
      cnt = cnt + 1
      RelPower(iz, ixy)  = hRelPower(cnt)
    END DO
  END DO
END DO

DEALLOCATE(hRelPower)
END SUBROUTINE

SUBROUTINE cuInitPWDist()
IMPLICIT NONE
TYPE(superPin_Type), POINTER :: superPin(:)
REAL :: volsum, tmpvolsum
INTEGER :: nxy, nzCMFD
INTEGER :: ixy, iz, izf
INTEGER :: ierr

IF(linitfin) RETURN
linitfin = .TRUE.

nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
superPin => cuGeometry%superPin

volsum = 0.
!$OMP PARALLEL PRIVATE(iz)
DO izf = 1, nzCMFD
  iz = cuGeometry%planeMap(izf)
  !$OMP DO SCHEDULE(GUIDED) REDUCTION(+:volsum)
  DO ixy = 1, nxy
    IF(.NOT. superPin(ixy)%lFuel(iz)) CYCLE
    volsum = volsum + cuGeometry%PinVolFmF(ixy, izf)
  END DO
  !$OMP END DO
END DO
!$OMP END PARALLEL
#ifdef MPI_ENV
CALL MPI_ALLREDUCE(volsum, tmpvolsum, 1, MPI_DOUBLE_PRECISION,, MPI_SUM, MPI_CUDA_COMM, ierr)
volsum = tmpvolsum
#endif

inv_totvolf = 1._8/volsum

NULLIFY(superPin)
IF(cuDevice%lFuel) CALL cuSetCMConverter()

END SUBROUTINE

SUBROUTINE cuSetCMConverter()
USE geom, ONLY : core
IMPLICIT NONE

TYPE(superPin_Type), POINTER :: superPin(:)
REAL, POINTER :: hz(:), hzfm(:)
INTEGER, POINTER :: PinMap(:), fmRange(:,:)
REAL :: invhz, val
INTEGER :: ng, nxy, myzb, myze, nz
INTEGER :: nxyc, myzbf, myzef, nzf
INTEGER :: iz, ixy, ipin, izf, ixy_map
INTEGER :: ir, ic

ng = cuGeometry%ng
nxy = cuGeometry%nxy
nxyc = cuGeometry%nxyc

myzb = cuDevice%myzb
myze = cuDevice%myze
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
nz = myze - myzb + 1
nzf = cuDevice%nzCMFD

superPin => cuGeometry%superPin
hz => cuGeometry%hz
hzfm => cuGeometry%hzfm
PinMap => cuGeometry%pinMap
fmRange => cuGeometry%fmRange
CALL createCsr(CMConverter, nxy*nzf, nxy*nz, nxyc*nzf)

! DO ixy = 1, nxyc
!   ixy_map = pinMap(ixy)
!   ipin = superPin(ixy_map)%iFuelPin
!   DO iz = myzb, myze
!     invhz = 1._8 / hz(iz)
!     IF(.NOT. superPin(ixy_map)%lFuel(iz)) CYCLE
!     ir = (ipin-1) * nz + iz - myzb + 1
!     DO izf = fmRange(iz,1), fmRange(iz,2)
!       ic = (izf - myzbf) * nxyc + ixy
!       val = hzfm(izf) * invhz
!       CALL pushCsr(CMConverter, val, ir, ic)
!     END DO
!   END DO
! END DO
DO ipin = 1, nxy
  ixy_map = Core%Pin(ipin)%isuperPin
  ixy = cuGeometry%pinMapRev(ixy_map)
  DO iz = myzb, myze
    invhz = 1._8 / hz(iz)
    IF(.NOT. superPin(ixy_map)%lFuel(iz)) CYCLE
    ir = (ipin-1) * nz + iz - myzb + 1
    DO izf = fmRange(iz,1), fmRange(iz,2)
      ic = (izf - myzbf) * nxyc + ixy
      val = hzfm(izf) * invhz
      CALL pushCsr(CMConverter, val, ir, ic)
    END DO
  END DO
END DO
CALL finalizeSortCsr(CMConverter, .TRUE.)

NULLIFY(superPin, hz, hzfm, PinMap, fmRange)
END SUBROUTINE

END MODULE
#endif
