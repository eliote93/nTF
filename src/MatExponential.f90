  !#include <KrylovDefines.h>
#include <Depletion.h>

#ifdef CRAM_DIRECT
!#pragma once
!  INCLUDE 'mkl_pardiso.f90'
#endif
#ifdef __INTEL_MKL
MODULE MatExponential
  USE CSRMATRIX
  !USE mkl_service
  USE OMP_LIB
#ifdef __PGI
  USE IEEE_ARITHMETIC   !--- CNJ Edit : F2003 Standard
#endif
  IMPLICIT NONE
#ifdef CRAM_14
  INTEGER, PARAMETER::CRAM_Order = 14;
  COMPLEX(8), PARAMETER::Pole(7) = (/ (-8.897773186468888, 16.630982619902085), (-3.703275049423448, 13.656371871483268), &
    (-0.208758638250130, 10.991260561901260), (3.993369710578568, 6.004831642235037), &
    (5.089345060580624, 3.588824029027006), (5.623142572745977, 1.194069046343966), &
    (2.269783829231112, 8.461797973040221) /);
  COMPLEX(8), PARAMETER::Res(7) = (/ (-7.154288063589067e-5, 1.436104334854130e-4), (9.439025310736168e-3, -1.718479195848301e-2)  ,&
    (-3.763600387822696e-1, 3.351834702945010e-1), (-2.349823209108270e+1, -5.808359129714207)    ,&
    (4.693327448883129e+1, 4.564364976882776e+1), (-2.787516194014564e+1, -1.0214733999015645e+2) ,&
    (4.807112098832508, -1.320979383742872) /);
  REAL(8), PARAMETER::Res0 = 1.832174378254041e-14
#endif
#ifdef CRAM_16
  INTEGER, PARAMETER::CRAM_Order = 16;
  COMPLEX(8), PARAMETER::Pole(8) = (/(-1.0843917078696988026e1, 1.9277446167181652284e1), (-5.2649713434426468895, 1.6220221473167927305e1)  ,&
    (5.9481522689511774808, 3.5874573620183222829), (3.5091036084149180974, 8.4361989858843750826)          ,&
    (6.4161776990994341923, 1.1941223933701386874), (1.4193758971856659786, 1.0925363484496722585e1)        ,&
    (4.9931747377179963991, 5.9968817136039422260), (-1.4139284624888862114, 1.3497725698892745389e1)/);
  COMPLEX(8), PARAMETER::Res(8) = (/(-5.0901521865224915650e-7, -2.4220017652852287970e-5), (2.1151742182466030907e-4, 4.3892969647380673918e-3) ,&
    (1.1339775178483930527e2, 1.0194721704215856450e2), (1.5059585270023467528e1, -5.7514052776421819979)        ,&
    (-6.4500878025539646595e1, -2.2459440762652096056e2), (-1.4793007113557999718, 1.7686588323782937906)        ,&
    (-6.2518392463207918892e1, -1.1190391094283228480e1), (4.1023136835410021273e-2, -1.5743466173455468191e-1)/);
  REAL(8), PARAMETER::Res0 = 2.1248537104952237488e-16
#endif
  INTEGER, PARAMETER :: iLU0 = 1, Diag = 2

  CONTAINS

  SUBROUTINE MatExpKrylov_full(A, x0, x1, rank, numthread)
  IMPLICIT NONE
  REAL(8),INTENT(IN) :: A(:,:), x0(:)
  REAL(8),INTENT(INOUT) :: x1(:)
  INTEGER, INTENT(IN) :: rank
  INTEGER, INTENT(IN), OPTIONAL :: numthread

  INTEGER :: n
  REAL(8) :: beta
  REAL(8), ALLOCATABLE :: e1(:)
  REAL(8), ALLOCATABLE :: v(:,:), h(:,:)       ! v : normalized vectors, h : Hessenberg matrix
  REAL(8), ALLOCATABLE :: ExpH(:,:)
  REAL(8), ALLOCATABLE :: yapprox(:)
  INTEGER :: m          ! reduced order
  INTEGER :: i
  INTEGER :: live_nt

  IF (PRESENT(numthread)) THEN
    live_nt = min(omp_get_max_threads(), numthread)
  ELSE
    live_nt = omp_get_max_threads()
  END IF

  CALL omp_set_num_threads(live_nt)

  n = rank
  ALLOCATE(v(n,n), h(n,n))

  CALL ArnoldiProcess(m)

  ALLOCATE(ExpH(m,m), yapprox(m), e1(m))

  e1(:) = 0._8; e1(1) = 1._8;

  CALL MatExpSns_full(.FALSE., h(1:m,1:m), ExpH(1:m,1:m), e1, yapprox(1:m), m, live_nt)

  yapprox(1:m) = beta*yapprox(1:m)
  x1(1:n) = 0._8
  DO i = 1, m
    x1(1:n) = x1(1:n)+v(:,i)*yapprox(i)
  END DO
  DEALLOCATE(v, h, ExpH, yapprox)
  CONTAINS

  SUBROUTINE ArnoldiProcess(reduced_rank)
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: reduced_rank
  REAL(8), ALLOCATABLE :: p(:)
  INTEGER :: i, j, icount
  REAL(8) :: HmOut

  ALLOCATE(p(n))

  beta = dot_product(x0, x0); beta = sqrt(beta)
  v(:,1) = x0(1:n)/beta;

  icount = 0;
  reduced_rank = 0;
  h(:,:) = 0.;
  DO i = 1, n-1 ! n-1 for last h(i+1, i)
    p(:) = 0._8
    DO j = 1, n
      p(:) = p(:)+A(:,j)*v(j,i)
    END DO
    DO j = 1, i
      h(j,i) = dot_product(v(:,j), p(:))
      p(:) = p(:) - h(j,i)*v(:,j)
    END DO
    h(i+1, i) = dot_product(p,p); h(i+1,i) = sqrt(h(i+1,i))
    IF (i .EQ. 25) HmOut = ABS(h(i+1,i))
    IF (i .GT. 25) THEN
      HmOut = 0.5_8*(HmOut+ABS(h(i+1,i)))
      IF (HmOut .LT. 4._8) icount = icount+1
      IF (icount .GE. 5) THEN
        reduced_rank = i; exit
      END IF
    END IF
    v(:,i+1) = p(:)/h(i+1,i)
  END DO
  IF(reduced_rank .EQ. 0) reduced_rank = n-1

  DEALLOCATE(p)
  END SUBROUTINE

  END SUBROUTINE

  SUBROUTINE MatExpSns_full(lOutMat, A, ExpA, x0, x1, rank, numthreads)
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: lOutMat
  REAL(8), INTENT(IN) :: A(:,:), x0(:)
  INTEGER, INTENT(IN) :: rank
  REAL(8), INTENT(INOUT) :: ExpA(:,:), x1(:)
  INTEGER, INTENT(IN), OPTIONAL :: numthreads

  REAL(8), ALLOCATABLE :: ExpMat(:,:), Ahat(:,:)           ! Used for first Taylor Exponential and Scailing
  INTEGER :: i, j
  INTEGER :: m, SqrK, RecMulK, K                        ! K = SqrK + log_2(RecMulK) -> {exp(2**(-K)*A)}**(2**K)
  INTEGER :: O_Taylor
  REAL(8), ALLOCATABLE :: SqrdMat(:,:), SaveSqrdMat(:,:), x_temp(:)
  REAL(8) :: MaxInA, PrevCond
  INTEGER :: live_nt

  IF (present(numthreads)) THEN
    live_nt = min(omp_get_max_threads(), numthreads)
  ELSE
    live_nt = omp_get_max_threads()
  END IF

  CALL omp_set_num_threads(live_nt)

  m = rank
  ALLOCATE(ExpMat(m,m), Ahat(m,m), SqrdMat(m,m), SaveSqrdMat(m,m))
  IF (.NOT. lOutMat) ALLOCATE(x_temp(m))

  MaxInA = maxval(ABS(A))

  K = CEILING(log(DBLE(MaxInA))/log(2.))
  K = MIN(26,K)
  RecMulK = CEILING(log(DBLE(m))/log(2.))

  IF (K .LT. RecMulK) THEN
    SqrK = 0
    IF (K .LE. 0) THEN
      K = 1;
    END IF
    RecMulK = 2*K
  ELSE
    SqrK = K - RecMulK
    RecMulK = 2**RecMulK
  END IF

  Ahat(:,:) = A(:,:)/DBLE(2**K)

  !print*, 1./DBLE(2**K), SqrK, RecMulK
  !WRITE(58,*) m, SqrK, RecMulK, 1./DBLE(2**K)

  O_Taylor = 4; PrevCond = 1.;
  DO WHILE(.TRUE.)
    CALL MatExpTaylor_full(.TRUE., Ahat, ExpMat, x0, x1, O_Taylor, m, live_nt)
    SaveSqrdMat(:,:) = ExpMat
    !print*, ExpMat
    ! Squaring the matrix exponentials, {{exp(A)}**2}**2 ...
    DO i = 1, SqrK
      CALL dgemm('N','N', m, m, m, 1._8, SaveSqrdMat, m, SaveSqrdMat, m, 0._8, SqrdMat, m)
      SaveSqrdMat = SqrdMat
    END DO
    !DO i = 0, m-1
      !print*, SaveSqrdMat(:,i+1), '//'
    !END DO
    IF (lOutMat) THEN
      ! Recrusively multiplying, exp(A)*exp(A)*exp(A)*...
      ExpMat = SaveSqrdMat
      IF (SqrK .EQ. 0) SqrdMat = SaveSqrdMat
      DO i = 1, RecMulK-1
        SaveSqrdMat = SqrdMat
        CALL dgemm('N','N',m,m,m,1._8, SaveSqrdMat, m, ExpMat, m, 0._8, SqrdMat, m)
      END DO
      IF (.NOT. SqrdMat(1,1).LT.1.2) THEN
        IF (ABS(PrevCond-Sqrdmat(1,1)).LT. 1.e-8) EXIT
        O_Taylor = O_Taylor + 4
        PrevCond = SqrdMat(1,1)
        !print*, m, K
      ELSE
        ExpA(1:m,1:m) = SqrdMat(:,:)
        !print*, O_Taylor, m
        EXIT
      END IF
      IF (O_Taylor .GT. m) THEN
        print*, 'Non Convergable Exponential Matrix with SnS'
        STOP
      END IF
    ELSE
      x1(1:m) = x0(1:m)
      DO i = 1, RecMulK
        x_temp(:) = x1(1:m)
        CALL dgemv('N', m, m, 1._8, SaveSqrdMat, m, x_temp, 1, 0._8, x1, 1)
      END DO
#ifdef __PGI
      IF (ANY(ieee_is_NAN(x1)).OR.ANY(x1.LT.0.)) THEN
        O_Taylor = O_Taylor + 4
      ELSE
#else
      IF (ANY(ISNAN(x1)).OR.ANY(x1.LT.0.)) THEN
        O_Taylor = O_Taylor + 4
      ELSE
#endif
        EXIT
      END IF
      IF (O_Taylor .GT. m) THEN
        print*, 'Non Convergable Exponential Matrix with SnS'
        STOP
      END IF
    END IF
  END DO

  DEALLOCATE(ExpMat, Ahat, SqrdMat, SaveSqrdMat)
  IF (.NOT. lOutMat) DEALLOCATE(x_temp)

  END SUBROUTINE

  SUBROUTINE MatExpTaylor_full(lOutMat, A, ExpA, x0, x1, order, rank, numthreads)
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: lOutMat
  REAL(8), INTENT(IN) :: A(:,:), x0(:)
  REAL(8), INTENT(INOUT) :: ExpA(:,:), x1(:)
  INTEGER, INTENT(IN) :: order, rank
  INTEGER, INTENT(IN), OPTIONAL :: numthreads

  INTEGER :: i, j
  REAL(8) :: Coeff
  REAL(8), ALLOCATABLE :: SqrA(:,:), SaveSqrA(:,:), x_temp(:,:)
  INTEGER :: live_nt

  IF (present(numthreads)) THEN
    live_nt = min(omp_get_max_threads(), numthreads)
  ELSE
    live_nt = omp_get_max_threads()
  END IF

  CALL omp_set_num_threads(live_nt)

  IF (lOutMat) ALLOCATE(SqrA(rank,rank), SaveSqrA(rank,rank))
  IF (.NOT. lOutMat) ALLOCATE(x_temp(rank,2))

  ! Initialize Matrices, ExpA = I + A + ...
  IF (lOutMat) THEN
    ExpA(:,:) = 0._8; SqrA(:,:) = 0._8
    DO i = 1, rank
      ExpA(i,i) = 1._8
      SqrA(i,i) = 1._8
    END DO
  END IF

  ! Coeff = (Factorial of i), SqrA = A**i
  Coeff = 1._8; SaveSqrA(:,:) = SqrA(:,:)
  IF (.NOT. lOutMat) x1(1:rank) = x0(1:rank); x_temp(:,2) = x0(1:rank)
  DO i = 1, order
    Coeff = 1/dble(i)
    IF (lOutMat) THEN
      IF (mod(i,2).EQ.1) THEN
      CALL dgemm('N','N', rank, rank, rank, Coeff, SaveSqrA, rank, A, rank, 0._8, SqrA, rank)
      ExpA(:,:) = ExpA(:,:) + SqrA(:,:)
      ELSE
      CALL dgemm('N','N', rank, rank, rank, Coeff, SqrA, rank, A, rank, 0._8, SaveSqrA, rank)
      ExpA(:,:) = ExpA(:,:) + SaveSqrA(:,:)
      END IF
    ELSE
      x_temp(:,1) = Coeff*x_temp(:,2)
      x_temp(:,2) = 0._8
      DO j = 1, rank
        x_temp(:,2) = x_temp(:,2)+x_temp(j,1)*A(:,j)
      END DO
      x1(1:rank) = x1(1:rank)+x_temp(:,2)
    END IF
    !IF(i.eq.4)print*, ExpA(:,1)
  END DO

  IF (lOutMat) DEALLOCATE(SqrA, SaveSqrA)
  IF (.NOT. lOutMat) DEALLOCATE(x_temp)
  END SUBROUTINE

  SUBROUTINE MatExpCRAM_full(lOutMat, A, ExpA, x0, x1, rank)
  LOGICAL, INTENT(IN) :: lOutMat
  REAL(8), INTENT(IN) :: A(:,:), x0(:)
  REAL(8), INTENT(INOUT) :: ExpA(:,:), x1(:)
  INTEGER, INTENT(IN) :: rank

  INTEGER :: Half_Order, ierr, I_Work_LPCK, info_LPCK
  INTEGER :: i, j
  COMPLEX(8), ALLOCATABLE :: A_inv(:,:), Work_LPCK(:), x1_z(:)
  INTEGER, ALLOCATABLE :: ipiv(:)

  Half_Order = CRAM_Order/2
  ALLOCATE(A_inv(rank,rank),ipiv(rank))
  IF (lOutMat) THEN
    I_Work_LPCK = 64*rank;
    ALLOCATE(Work_LPCK(I_Work_LPCK))
    ExpA(:,:) = 0._8;
    DO i = 1, rank
      ExpA(i,i) = Res0
    END DO
    DO i = 1, Half_Order
      A_inv(:,:) = CMPLX(A(:,:))
      DO j  = 1, rank
        A_inv(j,j) = A_inv(j,j)-Pole(i)
      END DO
      CALL zgetrf(rank,rank,A_inv,rank,ipiv,ierr)
      CALL zgetri(rank,A_inv,rank,ipiv,Work_LPCK,I_Work_LPCK,info_LPCK)
      ExpA(:,:) = ExpA(:,:)+2*DBLE(Res(i)*A_inv(:,:))
    END DO
    DEALLOCATE(Work_LPCK)
  ELSE
    ALLOCATE(x1_z(rank))
    x1(1:rank) = Res0*x0(1:rank)
    DO i = 1, Half_Order
      A_inv(:,:) = CMPLX(A(:,:))
      DO j = 1, rank
        A_inv(j,j) = A_inv(j,j)-Pole(i)
      END DO
      x1_z(:) = x0(1:rank)
      CALL zgesv(rank, 1, A_inv, rank, ipiv, x1_z, rank, ierr)
      x1(1:rank) = x1(1:rank) + 2*DBLE(Res(i)*x1_z(:))
    END DO
    DEALLOCATE(x1_z)
  END IF
  DEALLOCATE(A_inv, ipiv)

  END SUBROUTINE

  SUBROUTINE MatExpKrylov_CSR(A_csr, x0, x1, numthreads)
  IMPLICIT NONE
  INCLUDE 'mkl_blas.fi'
  TYPE(CSR_DOUBLE),INTENT(IN) :: A_csr
  REAL(8), POINTER :: x0(:)
  REAL(8), POINTER :: x1(:)
  INTEGER, INTENT(IN), OPTIONAL :: numthreads

  INTEGER :: n, rank
  REAL(8) :: beta
  REAL(8), ALLOCATABLE :: e1(:)
  REAL(8), ALLOCATABLE :: v(:,:), h(:,:), h_small(:,:)       ! v : normalized vectors, h : Hessenberg matrix
  REAL(8), ALLOCATABLE :: ExpH(:,:)
  REAL(8), ALLOCATABLE :: yapprox(:)
  INTEGER :: m          ! reduced order
  INTEGER :: i
  INTEGER :: live_nt

  LOGICAL, PARAMETER :: lOutMat = .TRUE.
  IF (.NOT. A_csr%lFinalized) RETURN
  IF (present(numthreads)) THEN
    live_nt = min(omp_get_max_threads(), numthreads)
  ELSE
    live_nt = omp_get_max_threads()
  END IF

  rank = A_Csr%nr;
  n = rank
  ALLOCATE(v(n,n), h(n,n))

  CALL ArnoldiProcess_CSR(A_csr, x0, rank, v, h, m, beta)

  !print*, m

  ALLOCATE(h_small(m,m), ExpH(m,m), yapprox(m), e1(m))
  h_small(:,:) = h(1:m,1:m)
  DEALLOCATE(h)
!#define HGDBG
#ifdef HGDBG
  DO i = 1, m
    WRITE(120,'(100ES10.3)') h_small(:,i)
  END DO
  WRITE(120,*)
#endif
  e1(:) = 0._8; e1(1) = 1._8;

  CALL MatExpSns_full(lOutMat, h_small(1:m,1:m), ExpH(1:m,1:m), e1, yapprox(1:m), m, live_nt)
  !print*, 'yapprox'
  !print*, yapprox
  IF (lOutMat) THEN
    CALL dgemv('N', m, m, 1._8, ExpH, m, e1, 1, 0._8, yapprox, 1)
  END IF

  yapprox(1:m) = beta*yapprox(1:m)
  CALL dgemv('N', n, m, 1._8, v(:,1:m), n, yapprox, 1, 0._8, x1, 1)
  DEALLOCATE(v, h_small, ExpH, yapprox, e1)
  END SUBROUTINE

  SUBROUTINE ArnoldiProcess_CSR(A_Csr, x0, rank, v, h, reduced_rank, beta)
  IMPLICIT NONE
  INCLUDE 'mkl_blas.fi'
  TYPE(CSR_DOUBLE),INTENT(IN) :: A_csr
  INTEGER :: rank
  REAL(8) :: v(:,:), h(:,:)
  REAL(8), POINTER :: x0(:)
  INTEGER, INTENT(INOUT) :: reduced_rank
  REAL(8), ALLOCATABLE :: p(:)
  INTEGER :: n, i, j, icount
  REAL(8) :: HmOut, hval, beta
  REAL(8), POINTER :: val(:)
  INTEGER, POINTER :: rowptr(:), colIdx(:)
  INTEGER :: nnz, nr

  nnz = A_Csr%nnz; nr = A_Csr%nr;
  n = nr;
  ALLOCATE(p(n))

  beta = dnrm2(n, x0, 1)
  v(:,1) = x0(1:n)/beta;
  !print*, beta

  icount = 0;
  reduced_rank = 0;
  h(:,:) = 0.;
  val => A_Csr%CsrVal; rowptr => A_Csr%csrRowPtr; colIdx => A_Csr%csrColIdx
  DO i = 1, n-1 ! n-1 for last h(i+1, i)
    CALL mkl_dcsrgemv('N', rank, val, rowptr, colIdx, v(:,i), p)
    !IF(i.EQ.1) print*, p
    DO j = 1, i
      hval = ddot(n, v(:,j), 1, p, 1)
      h(j,i) = hval
      !print*, i,j,hval
      p(:) = p(:) - hval*v(:,j)
    END DO
    !print*, i, h(i,i)
    !IF(i.eq.1) print*, p, '//'
    hval = dnrm2(n, p, 1)
    !print*, i, hval
    h(i+1, i) = hval
    IF (i .GT. 25) THEN
    IF (i .EQ. 26) HmOut = hval
      HmOut = 0.5_8*(HmOut+hval)
      IF (HmOut .LT. 4._8) icount = icount+1
      IF (icount .GE. 5) THEN
        reduced_rank = i; exit
      END IF
    END IF
    v(:,i+1) = p(:)/h(i+1,i)
  END DO
  IF(reduced_rank .EQ. 0) reduced_rank = n-1

  DEALLOCATE(p); NULLIFY(val, rowptr, colIdx)
  END SUBROUTINE

  SUBROUTINE MatExpSns_CSR(lOutMat, A_csr, ExpA_csr, x0, x1, numthread)
  USE AuxilCSR
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: lOutMat
  TYPE(CSR_DOUBLE), INTENT(IN) :: A_csr
  TYPE(CSR_DOUBLE), INTENT(INOUT) :: ExpA_csr
  REAL(8), INTENT(IN) :: x0(:)
  REAL(8), INTENT(INOUT) :: x1(:)
  INTEGER, INTENT(IN), OPTIONAL :: numthread

  INTEGER :: rank, live_nt
  TYPE(CSR_DOUBLE) :: ExpMat, Ahat          ! Used for first Taylor Exponential and Scailing
  INTEGER :: i, j
  INTEGER :: SqrK, RecMulK, K                        ! K = SqrK + log_2(RecMulK) -> {exp(2**(-K)*A)}**(2**K)
  INTEGER :: O_Taylor
  REAL(8), ALLOCATABLE :: SqrdMat(:,:), SaveSqrdMat(:,:), x_temp(:)
  REAL(8), ALLOCATABLE :: val(:)
  INTEGER, ALLOCATABLE :: rowptr(:), colIdx(:)
  INTEGER :: nr, nc, nnz
  REAL(8) :: MaxInA
  IF (.NOT. A_csr%lFinalized) RETURN
  IF (.NOT. present(numthread)) THEN
    live_nt = omp_get_max_threads()
  ELSE
    live_nt = min(omp_get_max_threads(), numthread)
  END IF
  CALL omp_set_num_threads(live_nt)

  rank = A_csr%nr;

  ALLOCATE(SqrdMat(rank,rank), SaveSqrdMat(rank,rank))
  IF (.NOT. lOutMat) ALLOCATE(x_temp(rank))
  Ahat = A_csr

  MaxInA = maxval(ABS(A_csr%csrVal))

  K = CEILING(log(DBLE(MaxInA))/log(2.))
  RecMulK = CEILING(log(DBLE(rank))/log(2.))

  IF (K .LT. RecMulK) THEN
    SqrK = 0
    IF (K .LE. 0) THEN
      K = 1;
    END IF
    RecMulK = 2**K
  ELSE
    SqrK = K - RecMulK
    RecMulK = 2**RecMulK
  END IF

  Ahat%csrVal(:) = Ahat%csrVal(:)/DBLE(2**K)

  O_Taylor = 8
  DO WHILE(.TRUE.)
    CALL MatExpTaylor_CSR(.TRUE., Ahat, ExpMat, x0, x1, O_Taylor, live_nt)
    !DO i = 1, rank
      !print*, ExpMat
    !END DO
    ! Squaring the matrix exponentials, {{exp(A)}**2}**2 ...
    CALL CSR2Full(SqrdMat, ExpMat)
    DO i = 1, SqrK
      CALL dgemm('N','N', rank,rank,rank, 1._8, SqrdMat, rank, SqrdMat, rank, 0._8, SaveSqrdMat, rank)
      SqrdMat(:,:) = SaveSqrdMat(:,:)
    END DO
    CALL Full2CSR(SqrdMat, ExpMat)
    nnz = ExpMat%nnz; nr = ExpMat%nr; nc = ExpMat%nc;
    ALLOCATE(val(nnz), rowptr(nr+1), colIdx(nnz))
    val = ExpMat%csrVal(:); rowptr(:) = ExpMat%csrRowPtr(:); colIdx(:) = ExpMat%csrColIdx(:)
    IF (lOutMat) THEN
      ! Recrusively multiplying, exp(A)*exp(A)*exp(A)*...
      DO i = 1, RecMulK-1
        CALL mkl_dcsrmm('N', rank,rank,rank, 1._8, 'G', val, colIdx, rowptr(1:nr), rowptr(2:nr+1), &
          SqrdMat, rank, 0._8, SaveSqrdMat, rank)
        SqrdMat = SaveSqrdMat
      END DO
#ifdef __PGI
      IF (ANY(ieee_IS_NAN(SqrdMat))) THEN
        O_Taylor = O_Taylor + 4
      ELSE
#else
      IF (ANY(ISNAN(SqrdMat))) THEN
        O_Taylor = O_Taylor + 4
      ELSE
#endif
        CALL Full2CSR(SqrdMat, ExpA_csr)
        CALL mkl_dcsrgemv('N', rank, ExpA_csr%csrVal, ExpA_csr%csrRowPtr, ExpA_csr%csrColIdx, x0, x1)
        !print*, O_Taylor, m
        EXIT
      END IF
      IF (O_Taylor .GT. rank) THEN
        print*, 'Non Convergable Exponential Matrix with SnS'
        STOP
      END IF
    ELSE
      x_temp(:) = x0(1:rank)
      DO i = 1, RecMulK
        CALL mkl_dcsrgemv('N', rank, val, rowptr, colIdx, x_temp, x1)
        x_temp(:) = x1(1:rank)
      END DO
#ifdef __PGI
      IF (ANY(ieee_IS_NAN(SqrdMat))) THEN
        O_Taylor = O_Taylor + 4
      ELSE
#else
      IF (ANY(ISNAN(SqrdMat))) THEN
        O_Taylor = O_Taylor + 4
      ELSE
#endif
        EXIT
      END IF
      IF (O_Taylor .GT. rank) THEN
        print*, 'Non Convergable Exponential Matrix with SnS'
        STOP
      END IF
    END IF
    DEALLOCATE(val, rowptr, colIdx)
  END DO

  CALL destroyCsr(Ahat)
  CALL destroyCsr(ExpMat)
  DEALLOCATE(SqrdMat, SaveSqrdMat)
  IF (.NOT. lOutMat) DEALLOCATE(x_temp)

  END SUBROUTINE

  SUBROUTINE MatExpTaylor_CSR(lOutMat, A_Csr, ExpA_Csr, x0, x1, order, numthread)
  USE AuxilCSR
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: lOutMat
  TYPE(CSR_DOUBLE), INTENT(IN) :: A_Csr
  TYPE(CSR_DOUBLE), INTENT(INOUT) :: ExpA_Csr
  REAL(8), INTENT(IN) :: x0(:)
  REAL(8), INTENT(INOUT) :: x1(:)
  INTEGER, INTENT(IN) :: order
  INTEGER, INTENT(IN), OPTIONAL :: numthread

  INTEGER :: rank, live_nt
  INTEGER :: i, j, ierr
  INTEGER :: nnz, nr, nc
  INTEGER, ALLOCATABLE :: pntrb(:), pntre(:)
  REAL(8), ALLOCATABLE :: A_val(:), ExpFull(:,:), B(:,:), C(:,:), x_temp(:), x_temp1(:)
  INTEGER, ALLOCATABLE :: A_rowptr(:), A_colIdx(:)
  REAL(8) :: Coeff
  IF (.NOT. A_csr%lFinalized) RETURN
  IF (.NOT. present(numthread)) THEN
    live_nt = omp_get_max_threads()
  ELSE
    live_nt = min(omp_get_max_threads(), numthread)
  END IF
  CALL omp_set_num_threads(live_nt)
  nnz = A_Csr%nnz; nr = A_Csr%nr; nc = A_Csr%nc
  rank = nr
  ALLOCATE(A_val(nnz), A_rowptr(nr+1), A_colIdx(nnz))
  A_val(:) = A_Csr%CsrVal(:); A_rowptr(:) = A_Csr%csrRowPtr(:); A_colIdx(:) = A_Csr%csrColIdx(:)

  IF (lOutMat) THEN
    ALLOCATE(pntrb(rank), pntre(rank))
    pntrb(1:rank) = A_rowPtr(1:rank)
    pntre(1:rank) = A_rowptr(2:rank+1)
    ALLOCATE(ExpFull(rank,rank), B(rank, rank), C(rank, rank))
    ExpFull(:,:) = 0._8
    B(:,:) = 0._8
    DO i = 1, rank
      ExpFull(i,i) = 1._8
      B(i,i) = 1._8
    END DO
    DO i = 1, order
      Coeff = 1._8/dble(i)
      CALL mkl_dcsrmm('N', rank,rank,rank, Coeff, 'G', A_val, A_ColIdx, pntrb, pntre, &
        B, rank, 0._8, C, rank)
      ExpFull(:,:) = ExpFull(:,:) + C(:,:)
      B(:,:) = C(:,:)
    END DO
    CALL Full2CSR(ExpFull, ExpA_csr)
    DEALLOCATE(pntrb, pntre, ExpFull)
    ! -----------------------------------------------------------------------------------------------
    !CALL mkl_dcsrgemv('N', rank, ExpA_Csr%csrVal, ExpA_csr%csrrowptr, ExpA_csr%csrcolIdx, x0, x1)   !
    ! -----------------------------------------------------------------------------------------------
  ELSE
    ALLOCATE(x_temp(rank),x_temp1(rank))
    x1(1:rank) = x0(1:rank)
    x_temp(:) = x0(1:rank)
    DO i = 1, order
      Coeff = 1./dble(i)
      x_temp1(:) = x_temp(:)
      CALL mkl_dcsrgemv('N', rank, A_val, A_rowptr, A_colIdx, x_temp1(:), x_temp)
      x_temp(:) = x_temp(:)*Coeff
      x1(1:rank) = x1(1:rank) + x_temp(:)
    END DO
    DEALLOCATE(x_temp, x_temp1)
  END IF

  DEALLOCATE(A_val, A_rowptr, A_colIdx)
  END SUBROUTINE

  SUBROUTINE MatExpCRAM_CSR(lOutMat, A_Csr, ExpA_Csr, x0, x1, numthread)
  USE AuxilCSR
  USE MKL_PARDISO
  IMPLICIT NONE

  LOGICAL, INTENT(IN) :: lOutMat
  TYPE(CSR_DOUBLE), INTENT(IN) :: A_Csr
  TYPE(CSR_DOUBLE), INTENT(INOUT) :: ExpA_Csr
  REAL(8), INTENT(IN) :: x0(:)
  REAL(8), INTENT(INOUT)  :: x1(:)
  INTEGER, INTENT(IN), OPTIONAL :: numthread

  INTEGER :: i, j, k
  INTEGER :: nr, nc, nnz, rank
  INTEGER :: Half_Order
  COMPLEX(8), ALLOCATABLE :: AZ(:,:)
  TYPE(CSR_DOUBLE_COMPLEX) :: A_CsrZ
  COMPLEX(8), ALLOCATABLE :: valZ(:), valZ_batch(:)
  INTEGER, ALLOCATABLE :: rowptr(:), colIdx(:), eyeIdx(:)
  INTEGER, ALLOCATABLE :: rowptr_batch(:), colIdx_batch(:)

  INTEGER :: live_nt
  TYPE(MKL_PARDISO_HANDLE) :: pt(64)
  INTEGER :: iparm(64)
  INTEGER :: maxfct, mnum, mtype, phase, n, nrhs, error, msglvl
  INTEGER, ALLOCATABLE :: perm(:)
  COMPLEX(8), ALLOCATABLE :: x_temp(:), x_temp_sol(:)
  COMPLEX(8) :: temp_pole
  IF (.NOT. A_csr%lFinalized) RETURN
  IF (lOutMat) THEN
    print*, 'Exponential Matrix Output with CSR Format Is Not Serviced on CRAM_CSR'
    print*, '( lOutMat = .TRUE.  is not available )'
    return
  END IF

  IF (PRESENT(numthread)) THEN
    live_nt = min(omp_get_max_threads(), numthread)
  ELSE
    live_nt = omp_get_max_threads()
  END IF
  CALL omp_set_num_threads(live_nt)

  nr = A_Csr%nr;  rank = nr
  Half_Order = CRAM_Order/2

  IF (lOutMat) THEN
    nrhs = nr
  ELSE
    nrhs = 1
  END IF
  maxfct = 1; mnum = 1; mtype = 13; phase = 13; msglvl = 0; error = 0
  ALLOCATE(perm(nr))
  iparm(:) = 0; perm(:) = 0
  DO i = 1, 64
    pt(i)%DUMMY = 0
  END DO

  ! iparm configuration
  iparm(1) = 1; iparm(2) = 0; iparm(8) = 1
  iparm(10) = 8; iparm(18) = 0; iparm(19) = -1

  !ALLOCATE(AZ(rank,rank),x_temp(rank), x_temp_sol(rank))     ! one by one
  ALLOCATE(AZ(rank,rank), x_temp(rank*Half_Order), x_temp_sol(rank*Half_Order)) ! batched

  nnz = A_Csr%nnz; nr = A_Csr%nr; nc = A_Csr%nc;
  ALLOCATE(valZ(nnz), rowptr(nr+1), colIdx(nnz), eyeIdx(nr))
  rowptr(:) = A_Csr%csrRowPtr(:); colIdx(:) = A_Csr%csrColIdx(:); valZ(:) = A_Csr%csrVal(:)
  DO i = 1, nr
    DO j = rowptr(i),rowptr(i+1)-1
      IF (colIdx(j) .EQ. i) THEN
        eyeIdx(i) = j; CYCLE
      END IF
    END DO
  END DO

  ALLOCATE(valZ_batch(nnz*Half_Order), rowptr_batch(nr*Half_Order+1), colIdx_batch(nnz*Half_Order))  ! batched
  x1(1:rank) = x0(1:rank)*Res0

  ! batched
  DO i = 1, Half_Order
    rowptr_batch((i-1)*nr+1: i*nr)=rowptr(:)+(i-1)*nnz; colIdx_batch((i-1)*nnz+1:i*nnz)=colIdx(:)+(i-1)*rank
    valZ_batch((i-1)*nnz+1:i*nnz) = valZ(:)
    x_temp((i-1)*rank+1:i*rank) = 2._8*Res(i)*x0(1:rank)
    temp_pole = Pole(i)!-1._8
    DO j = 1, nr
      valZ_batch((i-1)*nnz+eyeIdx(j)) = valZ(eyeIdx(j))-temp_pole
    END DO
  END DO
  rowptr_batch(nr*Half_Order+1) = nnz*Half_Order+1

  CALL createcsr(A_csrZ, nnz*Half_Order, nr*Half_Order, nr*Half_Order)
  A_csrZ%nnz = nnz*Half_Order; A_csrZ%nr = nr*Half_Order; A_csrZ%nc = nc*Half_Order;
  A_csrZ%csrval = valZ_batch; A_csrZ%csrrowptr = rowptr_batch; A_csrZ%csrcolidx = colIdx_batch
  A_csrZ%lFinalized = .true.

  !print*, 'Direct CRAM'
  CALL pardiso(pt, maxfct, mnum, mtype, phase, rank*Half_Order, valZ_batch, rowptr_batch, colIdx_batch, &
    perm, nrhs, iparm, msglvl, x_temp, x_temp_sol, error)
  phase = -1
  CALL pardiso(pt, maxfct, mnum, mtype, phase, rank*Half_Order, valZ_batch, rowptr_batch, colIdx_batch, &
    perm, nrhs, iparm, msglvl, x_temp, x_temp_sol, error)



  DO i = 1, Half_Order
    x1(1:rank) = x1(1:rank) + DBLE(x_temp_sol((i-1)*rank+1:i*rank))
  END DO


  ! one by one
  !DO i = 1, Half_Order
  !  valZ(:) = A_CsrZ%csrVal(:); temp_pole = Pole(i) - 1.;
  !  DO j = 1, nr
  !    valZ(eyeIdx(j)) = valZ(eyeIdx(j))-temp_pole
  !  END DO
  !  x_temp(:) = 2.*x0(:)*Res(i)
  !  CALL pardiso(pt, maxfct, mnum, mtype, phase, rank, valZ, rowptr, colIdx, &
  !    perm, nrhs, iparm, msglvl, x_temp, x_temp_sol, error)
  !  x1(:) = x1(:)+DBLE(x_temp_sol)
  !END DO


  !DEALLOCATE(valZ, rowptr, colIdx)
  DEALLOCATE(perm)
  DEALLOCATE(AZ, x_temp, x_temp_sol)
  DEALLOCATE(valZ, rowptr, colIdx, eyeIdx)
  DEALLOCATE(valZ_batch, rowptr_batch, colIdx_batch)
  CALL destroycsr(A_csrZ);

  END SUBROUTINE

  SUBROUTINE MatExpCRAM_Iter(lOutMat, A_Csr, ExpA_Csr, x0, x1, PreCond, numthread)
  USE AuxilCSR
  IMPLICIT NONE
  INCLUDE 'mkl_rci.fi'

  LOGICAL, INTENT(IN) :: lOutMat
  TYPE(CSR_DOUBLE), INTENT(IN) :: A_Csr
  TYPE(CSR_DOUBLE), INTENT(INOUT) :: ExpA_Csr
  REAL(8), INTENT(IN) :: x0(:)
  REAL(8), INTENT(INOUT)  :: x1(:)
  INTEGER :: PreCond
  INTEGER, INTENT(IN), OPTIONAL :: numthread

  INTEGER :: i, j, k
  INTEGER :: nr, nc, nnz, rank
  INTEGER :: Half_Order
  COMPLEX(8), ALLOCATABLE :: AZ(:,:)
  TYPE(CSR_DOUBLE_COMPLEX) :: A_CsrZ
  COMPLEX(8), ALLOCATABLE :: valZ(:), valZ_batch(:)
  INTEGER, ALLOCATABLE :: rowptr(:), colIdx(:), eyeIdx(:)
  INTEGER, ALLOCATABLE :: rowptr_batch(:), colIdx_batch(:)

  INTEGER :: live_nt
  TYPE(CSR_DOUBLE_COMPLEX) :: LU_csrZ
  INTEGER :: nr_bch
  COMPLEX(8), ALLOCATABLE :: r(:), p(:), phat(:), s(:), shat(:), t(:), v(:), rtilde(:), y(:)
  COMPLEX(8) :: rho0, rho1, rho2, a, b, w
  REAL(8) :: IterNum

  COMPLEX(8), ALLOCATABLE :: x_temp(:), x_temp_sol(:)
  COMPLEX(8) :: temp_pole
  IF (.NOT. A_csr%lFinalized) RETURN
  IF (lOutMat) THEN
    print*, 'Exponential Matrix Output with CSR Format Is Not Serviced on CRAM_CSR'
    print*, '( lOutMat = .TRUE.  is not available )'
    return
  END IF

  IF (PRESENT(numthread)) THEN
    live_nt = min(omp_get_max_threads(), numthread)
  ELSE
    live_nt = omp_get_max_threads()
  END IF
  CALL omp_set_num_threads(live_nt)

  nr = A_Csr%nr;  rank = nr
  Half_Order = CRAM_Order/2

  !ALLOCATE(AZ(rank,rank),x_temp(rank), x_temp_sol(rank))     ! one by one
  ALLOCATE(AZ(rank,rank), x_temp(rank*Half_Order), x_temp_sol(rank*Half_Order)) ! batched

  nnz = A_Csr%nnz; nr = A_Csr%nr; nc = A_Csr%nc;
  ALLOCATE(valZ(nnz), rowptr(nr+1), colIdx(nnz), eyeIdx(nr))
  rowptr(:) = A_Csr%csrRowPtr(:); colIdx(:) = A_Csr%csrColIdx(:); valZ(:) = A_Csr%csrVal(:)
  DO i = 1, nr
    DO j = rowptr(i),rowptr(i+1)-1
      IF (colIdx(j) .EQ. i) THEN
        eyeIdx(i) = j; CYCLE
      END IF
    END DO
  END DO

  ALLOCATE(valZ_batch(nnz*Half_Order), rowptr_batch(nr*Half_Order+1), colIdx_batch(nnz*Half_Order))  ! batched

  nr_bch = nr*Half_Order
  ALLOCATE(r(nr_bch), rtilde(nr_bch), p(nr_bch), phat(nr_bch), s(nr_bch), shat(nr_bch), t(nr_bch), v(nr_bch), y(nr_bch))

  x1(1:rank) = x0(1:rank)*Res0

  ! batched
  DO i = 1, Half_Order
    rowptr_batch((i-1)*nr+1: i*nr)=rowptr(:)+(i-1)*nnz; colIdx_batch((i-1)*nnz+1:i*nnz)=colIdx(:)+(i-1)*rank
    valZ_batch((i-1)*nnz+1:i*nnz) = valZ(:)
    x_temp((i-1)*rank+1:i*rank) = 2._8*Res(i)*x0(1:rank)
    temp_pole = Pole(i)!-1._8
    DO j = 1, nr
      valZ_batch((i-1)*nnz+eyeIdx(j)) = valZ(eyeIdx(j))-temp_pole
    END DO
  END DO
  rowptr_batch(nr*Half_Order+1) = nnz*Half_Order+1

  CALL createcsr(A_csrZ, nnz*Half_Order, nr*Half_Order, nr*Half_Order)
  A_csrZ%nnz = nnz*Half_Order; A_csrZ%nr = nr*Half_Order; A_csrZ%nc = nc*Half_Order;
  A_csrZ%csrval = valZ_batch; A_csrZ%csrrowptr = rowptr_batch; A_csrZ%csrcolidx = colIdx_batch
  A_csrZ%lFinalized = .true.

  !print*, 'Iterative CRAM'
  !print*, '0', LU_csrZ%nnz, LU_csrZ%lAlloc
  IF (PreCond .EQ. iLU0) CALL DeplZcsr_Bilu0(A_csrZ, LU_csrZ, eyeidx, nr, Half_Order)
  !print*, '1', LU_csrZ%nnz, LU_csrZ%lAlloc

  x_temp_sol = 0.; r = 0.
  !CALL mkl_zcsrgemv('N', nr_bch, valZ_batch, rowptr_batch, colidx_batch, x_temp_sol, r)
  r(:) = x_temp(:)-r(:)
  rtilde = r

  rho1 = DOT_PRODUCT(rtilde,r)
  DO i = 1, 10
    IF (i.EQ.1) THEN
      p = r; rho0 = rho1
    ELSE
      b = (rho1/rho2)*(a/w)
      p(:) = r(:) + b*(p(:)-w*v(:))
    END IF
    IF (PreCond .EQ. iLU0) THEN
      CALL mkl_zcsrtrsv('L','n','u',nr_bch,LU_csrZ%csrval,LU_csrZ%csrrowptr,LU_csrZ%csrcolidx,p,y)
      CALL mkl_zcsrtrsv('U','n','n',nr_bch,LU_csrZ%csrval,LU_csrZ%csrrowptr,LU_csrZ%csrcolidx,y,phat)
    ELSE
      DO j = 1, Half_Order
        phat(nr*(j-1)+1:nr*j) = p(nr*(j-1)+1:nr*j)/valZ_batch((eyeidx(:)+nnz*(j-1)))
      END DO
    END IF
    CALL mkl_zcsrgemv('N', nr_bch, valZ_batch, rowptr_batch, colidx_batch, phat, v)
    !print*, phat
    a = rho1/(DOT_PRODUCT(rtilde,v))
    s(:) = r(:) - a*v(:)
    rho2 = DOT_PRODUCT(rtilde,s)
    IF (.NOT.(CDABS(rho2) .LT. 0. .OR. CDABS(rho2) .GE. 0)) THEN
      IterNum = i-1
      EXIT
    END IF

    x_temp_sol(:) = x_temp_sol(:)+a*phat(:)
    IF (CDABS(rho2/rho0) < 1.e-30) THEN
      IterNum = i-1+0.5
      EXIT
    END IF

    IF (PreCond .EQ. iLU0) THEN
      CALL mkl_zcsrtrsv('L','n','u',nr_bch,LU_csrZ%csrval,LU_csrZ%csrrowptr,LU_csrZ%csrcolidx,s,y)
      CALL mkl_zcsrtrsv('U','n','n',nr_bch,LU_csrZ%csrval,LU_csrZ%csrrowptr,LU_csrZ%csrcolidx,y,shat)
    ELSE
      DO j = 1, Half_Order
        shat(nr*(j-1)+1:nr*j) = s(nr*(j-1)+1:nr*j)/valZ_batch((eyeidx(:)+nnz*(j-1)))
      END DO
    END IF
    CALL mkl_zcsrgemv('N', nr_bch, valZ_batch, rowptr_batch, colidx_batch, shat, t)

    w = DOT_PRODUCT(t,s)/DOT_PRODUCT(t,t)
    r(:) = s(:)-w*t(:)
    rho1 = DOT_PRODUCT(rtilde,r)
    rho2 = rho1
    IF (.NOT.(CDABS(rho2) .LT. 0. .OR. CDABS(rho2) .GE. 0)) THEN
      IterNum = i-1+0.5
      EXIT
    END IF

    x_temp_sol(:) = x_temp_sol(:)+w*shat(:)
    IF (CDABS(rho2/rho0) < 1.e-30) THEN
      IterNum = i
      EXIT
    END IF
  END DO
  !WRITE(97,*) IterNum, CDABS(rho2), ANY(ISNAN(DBLE(x_temp_sol)))
  !WRITE(97,*) 'Iter : ', IterNum
  DO i = 1, Half_Order
    x1(1:rank) = x1(1:rank) + DBLE(x_temp_sol((i-1)*rank+1:i*rank))
  END DO

  !DEALLOCATE(valZ, rowptr, colIdx)
  DEALLOCATE(r, rtilde, p, phat, s, shat, t, v, y)
  CALL destroycsr(LU_csrZ);
  !print*, '3', LU_csrZ%nnz, LU_csrZ%lAlloc
  DEALLOCATE(AZ, x_temp, x_temp_sol)
  DEALLOCATE(valZ, rowptr, colIdx, eyeIdx)
  DEALLOCATE(valZ_batch, rowptr_batch, colIdx_batch)
  CALL destroycsr(A_csrZ);

  END SUBROUTINE
! ***************************** Reference Page ****************************.
! https://www.cfd-online.com/Wiki/Incomplete_LU_factorization_-_ILU        |
! https://www.cc.gatech.edu/~echow/pubs/parilu-sisc.pdf                    |
! https://www.cfd-online.com/Wiki/Biconjugate_gradient_stabilized_method   |
! *************************************************************************'

  SUBROUTINE DeplZcsr_Bilu0(Z_csr, ZLU_csr, eyeidx, nr4one, NofB)
  !-----------------------------------------------------------!
  !                   No Off-Diagonal Block                   !
  !-----------------------------------------------------------!
  IMPLICIT NONE
  TYPE(CSR_DOUBLE_COMPLEX) :: Z_csr, ZLU_csr
  INTEGER :: eyeidx(:), nr4one, NofB

  COMPLEX(8), POINTER :: val(:)
  INTEGER, POINTER :: colidx(:), rowptr(:)

  INTEGER, ALLOCATABLE :: colptr(:), rowidx(:)
  INTEGER :: iB, nnz, nnz1
  INTEGER :: i, ir, ic

  ZLU_csr = Z_csr;
  nnz = ZLU_csr%nnz;
  val => ZLU_csr%csrval
  colidx => ZLU_csr%csrColIdx(1:nnz1)
  rowptr => ZLU_csr%csrRowPtr(1:nr4one+1)

  nnz1 = nnz/NofB

  ALLOCATE(colptr(nr4one+1), rowidx(nnz1))
  colptr = 0.; rowidx = 0.;
  DO i = 1, nnz1
    colptr(colIdx(i)+1) = colptr(colIdx(i)+1)+1;
  END DO
  colptr(1) = 1;
  DO i = 1, nr4one
    colptr(i+1) = colptr(i+1)+colptr(i)
  END DO
  ir = 0
  DO i = 1, nnz1
    IF(rowptr(ir+1) .EQ. i) ir=ir+1
    ic = colidx(i)
    rowidx(colptr(ic)) = ir
    colptr(ic)=colptr(ic)+1
  END DO
  DO ic = nr4one,2,-1
    colptr(ic) = colptr(ic-1)
  END DO
  colptr(1) = 1;

  ! Need to sweep to csc

  DO iB = 0, NofB-1
    CALL zcsr_ilu0(val(iB*nnz1+1:(iB+1)*nnz1), colidx, rowptr, eyeidx, colptr, rowidx, nnz1, nr4one)
  END DO

  DEALLOCATE(colptr, rowidx)
  END SUBROUTINE

  SUBROUTINE zcsr_ilu0(val, colidx, rowptr, eyeidx, colptr, rowidx, nnz, nr)
  IMPLICIT NONE
  INTEGER :: nnz, nr
  COMPLEX(8) :: val(:)
  INTEGER :: colidx(:), rowptr(:), eyeidx(:)
  INTEGER :: colptr(:), rowidx(:)

  COMPLEX(8) :: invD, e
  INTEGER :: i, j, k, l
  INTEGER :: ir, ic, jc
  INTEGER :: irL, irU, icL, icU
  INTEGER :: jrL, jrU, jcL, jcU

  DO i = 1, nr-1
    invD = 1./val(eyeidx(i))
    irL = colptr(i); irU = colptr(i+1)-1
    icL = eyeidx(i)+1; icU = rowptr(i+1)-1
    DO j = irL, irU
      ir = rowidx(j)
      IF (ir .LE. i) CYCLE
      jcL = rowptr(ir); jcU = rowptr(ir+1)-1
      DO k = jcL, jcU
        ic = colidx(k)
        IF (ic .EQ. i) THEN
          e=val(k)*invD; val(k) = e
        ELSEIF (ic .GT. i) THEN
          IF (dble(e) .LE. 1.e-30) EXIT
          DO l = icL, icU
            jc = colidx(l)
            IF (jc .EQ. ic) val(l) = val(l) - e*val(k)
            IF (jc .GE. ic) EXIT
          END DO
        END IF
      END DO
      e = 0.
    END DO
  END DO

  END SUBROUTINE

END MODULE
#endif
