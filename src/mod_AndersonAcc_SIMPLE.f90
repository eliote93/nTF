MODULE Anderson_Acceleration_SIMPLETH
  USE PARAM
  USE TYPEDEF,   ONLY: THINFO_TYPE, PE_TYPE, FuelTh_Type, CoolantTH_Type
  IMPLICIT NONE

  private
  public :: init_AA, Anderson_acc, m_AA
  integer :: m_AA, ndim_aa, nz_aa, nrad_aa, nrf
  real*8, allocatable, dimension(:) :: g_new, f_old, delf, x_solaa, g_old, gam, f_new
  real*8, allocatable, dimension(:,:) :: Q,R,DG,DF,DF_dag
  real*8 :: beta_AA = 1.0_8
  LOGICAL :: lDirectInv

  interface matrixmul
    module procedure matrixmul2D2D
    module procedure matrixmul2D1D
  end interface matrixmul

  contains

  SUBROUTINE init_AA(is_AA,nrad,nz,nrfp4)
    implicit none
    logical, intent(in):: is_AA
    integer, intent(in) :: nrad,nz,nrfp4

    nz_aa = nz
    nrad_aa = nrad
    nrf = nrfp4 - 4
    if (is_AA) then
      IF (m_AA.LE.2) lDirectInv = .TRUE.
      ndim_aa = nrad*nz*2 + nrad*nz*(nrfp4)
      allocate(g_new(ndim_aa), source = 0._8)
      allocate(g_old(ndim_aa), source = 0._8)
      allocate(x_solaa(ndim_aa), source = 0._8)
      allocate(gam(m_AA), source = 0._8)
      allocate(f_old(ndim_aa), source = 0._8)
      allocate(f_new(ndim_aa), source = 0._8)
      allocate(DG(m_AA, ndim_aa), source = 0._8)
      IF (lDirectInv) THEN
        ALLOCATE(DF(m_AA, ndim_aa), source = 0.)
        ALLOCATE(DF_dag(m_AA, ndim_aa), source = 0.)
      ELSE
        allocate(delf(ndim_aa), source = 0._8)
        allocate(Q(ndim_AA,m_AA), source = 0._8)
        allocate(R(m_AA,m_AA), source = 0._8)
      END IF
    else
      return
    end if

  END SUBROUTINE init_AA

  SUBROUTINE Anderson_acc(ncall,ThInfo,PE)
    USE IOUTIL,         ONLY : message
    USE FILES,          ONLY : IO8
    integer, intent (in) :: ncall
    TYPE(ThInfo_Type) :: ThInfo
    TYPE(PE_Type) :: PE
    integer :: ixy, iz, i, j, m, k, ii, nrp4, mk_cond
    integer, save :: mk = 0
    real(8) :: tol_cond = 1000.
    real(8) :: temp, rnorm, stemp, condn
    logical :: lqrdel, Master
    real(8),save :: rho_init, Tf_init, Tc_init, rho_temp, Tc_temp, Tf_temp
    TYPE(FuelTh_Type), POINTER :: FuelTH(:)
    TYPE(CoolantTH_Type), POINTER :: CoolantTH(:)

    REAL(8) :: Finv_d(2,2)

    MASTER = PE%CMFDMASTER
    FuelTH => ThInfo%FuelTH
    CoolantTH => ThInfo%CoolantTH
    nrp4=nrf+4
    lqrdel = .false.

    if (ncall == 0) then
      mk = 0
      rho_init = 0._8
      Tc_init = 0._8
      Tf_init = 0._8
      do ixy=1,nrad_aa
        rho_temp = maxval(CoolantTH(ixy)%DenCool(:))
        rho_init = max(rho_temp,rho_init)
        Tc_temp = maxval(CoolantTH(ixy)%Tcool(:))
        Tc_init = max(Tc_temp,Tc_init)
        Tf_temp = maxval(FuelTH(ixy)%Tfuel(:,:))
        Tf_init=max(Tf_temp,Tf_init)
      end do
    endif
    !print*, rho_init, Tc_init, Tf_init

    k=0
    do ixy=1,nrad_aa
      do iz=1,nz_aa
        k=k+1
        IF(.NOT. CoolantTH(ixy)%lfuel) then
          g_new(k)=1.D-12
        ELSE
          g_new(k)=CoolantTH(ixy)%DenCool(iz)/rho_init
        ENDIF
      end do
    end do

    do ixy=1,nrad_aa
      do iz=1,nz_aa
        k=k+1
        IF(.NOT. CoolantTH(ixy)%lfuel) then
          g_new(k)=1.D-12
        ELSE
          g_new(k)=CoolantTH(ixy)%Tcool(iz)/Tc_init
        END IF
      end do
    end do

    do ixy=1,nrad_aa
      do iz=1,nz_aa
        do m=1,nrp4
          k=k+1
          IF (CoolantTH(ixy)%lfuel) then
            g_new(k)=FuelTH(ixy)%Tfuel(m,iz)/Tf_init
          else
            g_new(k)=1.D-12
          end if
        end do
      end do
    end do

    do i = 1,ndim_aa
      f_new(i)=g_new(i)-x_solaa(i)
    enddo

    if (ncall>0) then
      mk = mk +1
      IF (.NOT.lDirectInv) THEN
        do i=1,ndim_aa
          delf(i) = f_new(i)-f_old(i)
        enddo
      END IF
      if (mk > m_AA) then
        mk = mk -1
        lqrdel = .true.
        do i=1,ndim_aa
          do m=2,mk
            DG(m-1,i)=DG(m,i)
          enddo
          IF (lDirectInv) THEN
            do m=2,mk
              DF(m-1,i)=DF(m,i)
            enddo
          END IF
        enddo
      endif
      IF (lDirectInv) THEN
        do i =1, ndim_aa
          DF(mk,i) = f_new(i)-f_old(i)
          DG(mk,i) = g_new(i)-g_old(i)
        enddo
      ELSE
        DO i = 1, ndim_aa
          DG(mk,i) = g_new(i)-g_old(i)
        ENd DO
      END IF

    endif

    do i=1,ndim_aa
      g_old(i) = g_new(i)
      f_old(i) = f_new(i)
    enddo

    if (mk == 0) then
      do i=1,ndim_aa
        x_solaa(i) = g_new(i)
      enddo
      return
    endif

    IF (lDirectInv) THEN !! Edit by LHG, 2021.01.14
      mk_cond = 1
      DO WHILE ((mk-mk_cond).GE.0)
        Finv_d = 0.

        ! ---- Determine F_reduced = tr(DF)*DF
        DO i = mk_cond, mk
          DO j = mk_cond, mk
            Finv_d(i,j) = dot_product(DF(i,:),DF(j,:))
          END DO
        END DO

        ! ---- Get the inverse of F_reduced
        IF ((mk-mk_cond).EQ.1) THEN
          temp = Finv_d(1,1)*Finv_d(2,2)-Finv_d(2,1)*Finv_d(1,2)
          Finv_d = Finv_d/temp
          temp = Finv_d(1,1)
          Finv_d(1,1) = Finv_d(2,2)
          Finv_d(2,2) = temp
          Finv_d(1,2) = -Finv_d(1,2); Finv_d(2,1) = -Finv_d(2,1)
        ELSE
          Finv_d(1,1) = 1./Finv_d(1,1)
          Finv_d(2,2) = 1./Finv_d(2,2)
        END IF

        ! ---- Obtain DF^(dagger)
        DF_dag = MATMUL(Finv_d,DF)

        ! ---- Condition number calculate, cond(DF) = ||DF||*||DF_dag||
        condn = 0.; temp = 0.;
        DO i = mk_cond, mk
          condn = condn+dot_product(DF_dag(i,:),DF_dag(i,:))
          temp = temp+dot_product(DF(i,:),DF(i,:))
        END DO
        condn = condn*temp

        IF (condn.LT.tol_cond) EXIT
        mk_cond = mk_cond+1
      END DO

      ! ---- Picard iteration solution
      x_solaa(:) = g_new(:)

      IF ((mk-mk_cond+1)*mk.NE.0) THEN
        ! ---- Solving least-square for min_{gam} (||f_new - DG*gam||)
        DO i = mk_cond, mk
          gam(i) = dot_product(DF_dag(i,:),f_new(:))
        END DO
        WRITE(MESG, '(A,I2,A,20f10.5)') 'AA-',(mk-mk_cond+1), ": Gamma = ", gam(mk_cond:mk)
        if(MASTER) CALL MESSAGE(io8,TRUE,TRUE,MESG)

        ! ---- Anderson acc. solution without damping as (1-beta_AA)
        DO i = mk_cond, mk
          x_solaa(:) = x_solaa(:)-gam(i)*DG(i,:)
        END DO

        ! ---- Damping the solution
        IF (beta_AA < 1.0_8) THEN
          x_solaa(:) = beta_AA*x_solaa(:)-(1.-beta_AA)*(f_new(:)-gam(mk_cond)*DF(mk_cond,:))
          DO i = mk_cond+1, mk
            x_solaa(:) = x_solaa(:)+(1.-beta_AA)*gam(i)*DF(i,:)
          END DO
        END IF

      END IF
    END IF

    IF (.NOT. lDirectInv) THEN
      if (lqrdel) call qrdelete(Q,R,mk)

  !Single modified Gram-Schmidt for newly added column
      do i=1,mk-1
        temp = 0.
        do ii =1,ndim_aa
          temp = temp + Q(ii,i)*delf(ii)
        enddo
        R(i,mk) = temp
        do ii= 1,ndim_aa
          delf(ii) = delf(ii) - temp*Q(ii,i)
        enddo
      enddo

      R(mk,mk) = norm2(delf)
      rnorm = 1./R(mk,mk)
      do i=1,ndim_aa
        Q(i,mk) = delf(i)*rnorm
      enddo

  !Monitoring condition number of R

      condn = cond_num(R,mk)
      do while (condn>tol_cond .and. mk>1)
        call qrdelete(Q,R,mk)
        do i=1,ndim_aa
          do m=2,mk
            DG(m-1,i)=DG(m,i)
          enddo
        enddo
        mk = mk - 1
        condn = cond_num(R,mk)
      enddo

  ! Solve Rr=Q'f by backward substitution
      stemp =0.
      do i=1,ndim_aa
        stemp = stemp+Q(i,mk)*f_new(i)
      enddo

      gam(mk) = stemp/R(mk,mk)
      do i=mk-1,1,-1
        temp = 0.
        do j = i+1,mk
          temp = temp + R(i,j)*gam(j)
        enddo

        stemp =0.
        do j=1,ndim_aa
          stemp = stemp+Q(j,i)*f_new(j)
        enddo
        gam(i) = (stemp - temp)/R(i,i)
      enddo

      WRITE(MESG, '(A,I2,A,20f10.5)') 'AA-',mk, ": Gamma = ", gam(1:mk)
      if(MASTER) CALL MESSAGE(io8,TRUE,TRUE,MESG)

  ! Update solutions

      do i=1,ndim_aa
        x_solaa(i) = g_new(i)
        do ii=1,mk
          x_solaa(i) = x_solaa(i) -DG(ii,i)*gam(ii)
        enddo
      enddo

      if (beta_AA < 1.0 ) then
        call matrixmul(R(1:mk,1:mk), gam(1:mk),gam(1:mk))
        call matrixmul(Q(1:ndim_aa,1:mk), gam(1:mk),g_new)
        do i=1,ndim_aa
          x_solaa(i) = x_solaa(i) - (1.-beta_AA)*(f_new(i)-g_new(i))
        enddo
      endif
    END IF

    k=0
    do ixy=1,nrad_aa
      do iz=1,nz_aa
        k=k+1
        IF(.NOT. CoolantTH(ixy)%lfuel) cycle
        CoolantTH(ixy)%DenCool(iz)=x_solaa(k)*rho_init
      end do
    end do

    do ixy=1,nrad_aa
      do iz=1,nz_aa
        k=k+1
        IF(.NOT. CoolantTH(ixy)%lfuel) cycle
        CoolantTH(ixy)%Tcool(iz)=x_solaa(k)*Tc_init
      end do
    end do

    do ixy=1,nrad_aa
      do iz=1,nz_aa
        do m=1,nrp4
          k=k+1
          IF(.NOT. CoolantTH(ixy)%lfuel) cycle
          FuelTH(ixy)%Tfuel(m,iz)=x_solaa(k)*Tf_init
        end do
      end do
    end do

  END SUBROUTINE Anderson_acc

  subroutine qrdelete(Q,R, mk)
   real(8), intent (inout) :: Q(:,:), R(:,:)
   integer, intent (inout) :: mk
   integer :: i, j, ii
   real(8) :: temp, c, s
!Upper-Hessenberg to Upper triangle Matrix
    do i = 1,mk-1
      temp = sqrt(R(i,i+1)*R(i,i+1)+R(i+1,i+1)*R(i+1,i+1))
      c = R(i,i+1)/temp
      s = -R(i+1,i+1)/temp
      do j = i+1,mk
        temp = c*R(i,j)-s*R(i+1,j)
        R(i+1,j) = s*R(i,j)+c*R(i+1,j)
        R(i,j) = temp
      enddo

      do ii = 1,ndim_aa
        temp = c*Q(ii,i)-s*Q(ii,i+1)
        Q(ii,i+1) = s*Q(ii,i)+c*Q(ii,i+1)
        Q(ii,i) = temp
      enddo
    enddo
!Shift R
    do i=1,mk
      do j=i+1,mk
        R(i,j-1) = R(i,j)
      enddo
    enddo

  end subroutine


  real(8) function cond_num(R,ndim)
    ! Function to calculate conditioni number of matrix R
    ! Matrix R is an upper triangular matrix, and ndim is dimension of R
    real(8), intent (in) :: R(:,:)
    integer, intent (in) :: ndim
    real(8) :: norm_inf1, norm_inf2, temp(ndim), inv_R(ndim,ndim), inv_temp, div
    integer :: i, j, ii

    do i = 1,ndim
      temp(i) = 0.
      do j=i,ndim
        temp(i) = temp(i) + abs(R(i,j))
      enddo
    enddo
    norm_inf1 = maxval(temp)


    do j=1,ndim
      do i = 1,ndim
        if (i==j) then
          inv_R(i,j) = 1.
        else
          inv_R(i,j) = 0.
        endif
      enddo
    enddo

    do i=ndim,1,-1
      inv_temp = 1./(R(i,i))

      do ii = i-1,1,-1
        div = R(ii,i)*inv_temp
        do j = i,ndim
          inv_R(ii,j) = inv_R(ii,j) - div*inv_R(i,j)
        enddo
      enddo

      do j = i,ndim
        inv_R(i,j) = inv_R(i,j)*inv_temp
      enddo
    enddo

    do i = 1,ndim
      temp(i) = 0.
      do j=i,ndim
        temp(i) = temp(i) + abs(inv_R(i,j))
      enddo
    enddo
    norm_inf2 = maxval(temp)

    cond_num = norm_inf1*norm_inf2

  end function cond_num

  subroutine matrixmul2D1D(Mat1,vec,ans)
    real(8), intent (in) :: Mat1(:,:), vec(:)
    real(8), dimension(:), intent (inout) :: ans
    real(8), allocatable :: Mat1temp(:,:), vectemp(:)
    real(8) :: a
    integer :: I,rank1, rank2, dimx, k, rank3

    dimx=size(Mat1,1)
    rank1=size(Mat1,2)
    rank2=size(vec,1)
    rank3=size(ans,1)
    if (rank1/=rank2) stop "The rank of Mat1 and vector is not matched"
    if (dimx/=rank3) stop "The rank of Mat1 and solution vector is not matched"

    allocate(Mat1temp(rank1,dimx),vectemp(rank2))

    do k=1,rank1
      do i =1,dimx
        Mat1temp(k,i)=Mat1(I,k)
      enddo
      vectemp(k)=vec(k)
    enddo


    do I=1,dimx
      a=0.
      do k=1,rank1
        a=a+Mat1temp(k,i)*vectemp(k)
      enddo
      ans(I)=a
    enddo

    deallocate(Mat1temp,vectemp)
  end subroutine

  subroutine matrixmul2D2D(Mat1,Mat2,ans)
    real(8), intent (in) :: Mat1(:,:), Mat2(:,:)
    real(8), dimension(:,:), intent (inout) :: ans
    real(8), allocatable :: Mat1temp(:,:), Mat2temp(:,:)
    real(8) :: a
    integer :: I,J, rank1, rank2, dimx, dimy, k

    dimx=size(Mat1,1)
    dimy=size(Mat2,2)
    rank1=size(Mat1,2)
    rank2=size(Mat2,1)
    if (rank1/=rank2) stop "The rank of Mat1 and Mat2 is not matched"



    allocate(Mat1temp(dimx,rank1),Mat2temp(rank2,dimy))
    Mat1temp=Mat1
    Mat2temp=Mat2

    do J=1,dimy
      do I=1,dimx
        a=0.
        do k=1,rank1
          a=a+Mat1temp(I,k)*Mat2temp(k,J)
        enddo
        ans(I,J)=a
      enddo
    enddo

    deallocate(Mat1temp,Mat2temp)
  end subroutine

END MODULE
