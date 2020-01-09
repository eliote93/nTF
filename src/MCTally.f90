module MCTally
!private
 !   public :: InitMCTally, ResetAcc, AcmPinPow, WritePinPow, nthrd, nfxr, ng, Fxrphi, Fxrphiomp
    integer :: nx, ny, nz, ng, nthrd, nfxr
    real(8), pointer, dimension(:,:,:) :: sumpsi, sqsumpsi
    real(8), pointer, dimension(:,:) :: sumpsir, sqsumpsir
    real(8), pointer, dimension(:) :: sumpsiz, sqsumpsiz
    real(8), pointer, dimension(:,:,:) :: psi, stdpsi
    real(8), pointer, dimension(:,:) :: psir, stdpsir
    real(8), pointer, dimension(:) :: psiz, stdpsiz
    integer :: iacc
    logical :: lxslib
contains

subroutine InitMCTally(nxyz, nfxr, ng, nthrd, lxslib, lmicxs)
    use allocs
    use MCXsLib, ONLY : Fxrphi, FxrPhiOmp
    implicit none
    INTEGER :: nthrd, nfxr, ng
    LOGICAL :: lxslib, lmicxs
    integer, intent(in) :: nxyz(3)
    nx=nxyz(1);ny=nxyz(2);nz=nxyz(3)
    call dmalloc(sumpsi, nx,ny,nz)
    call dmalloc(sqsumpsi, nx,ny,nz)
    call dmalloc(sumpsir, nx,ny)
    call dmalloc(sqsumpsir, nx,ny)
    call dmalloc(sumpsiz, nz)
    call dmalloc(sqsumpsiz, nz)
    call dmalloc(psi, nx,ny,nz)
    call dmalloc(psir, nx,ny)
    call dmalloc(psiz, nz)
    call dmalloc(stdpsi, nx,ny,nz)
    call dmalloc(stdpsir, nx,ny)
    call dmalloc(stdpsiz, nz)
    IF( lxslib .AND. .NOT. lmicxs )THEN
        ALLOCATE(FxrPhi(nfxr,nz,ng))
        ALLOCATE(FxrPhiomp(nfxr,nz,ng,nthrd))
    !CALL Dmalloc(Fxrphi, nfxr, nz, ng)
    !CALL Dmalloc(Fxrphiomp, nfxr, nz, ng, nthrd)
        Fxrphi=1._8
        Fxrphiomp=0._8
    ENDIF
    iacc=0
end subroutine

subroutine ResetAcc
    implicit none
    integer :: k, j, tid
    iacc=0
    sumpsiz=0
    sqsumpsiz=0
    do k=1, nz
        sumpsir(:,k)=0
        sqsumpsir(:,k)=0
        do j=1, ny
            sumpsi(:,j,k)=0
            sqsumpsi(:,j,k)=0
        enddo   
    enddo
end subroutine

subroutine AcmPinPow(psi)
    implicit none
    real(8), pointer, dimension(:,:,:), intent(in) :: psi
    integer :: i,j,k
    
    iacc=iacc+1
    psiz=0
    psir=0
    do k=1, nz
        do j=1, ny
            do i=1, nx
                sumpsi(i,j,k)=sumpsi(i,j,k)+psi(i,j,k)
                sqsumpsi(i,j,k)=sqsumpsi(i,j,k)+psi(i,j,k)**2
                psiz(k)=psiz(k)+psi(i,j,k)
                psir(i,j)=psir(i,j)+psi(i,j,k)
            enddo
        enddo
    enddo

    do k=1, nz
        sumpsiz(k)=sumpsiz(k)+psiz(k)
        sqsumpsiz(k)=sqsumpsiz(k)+psiz(k)**2
    enddo
    do j=1, ny
        do i=1, nx
            sumpsir(i,j)=sumpsir(i,j)+psir(i,j)
            sqsumpsir(i,j)=sqsumpsir(i,j)+psir(i,j)**2
        enddo
    enddo
end subroutine

subroutine WritePinPow(ctrlMC)
    use MCMap
    USE files, only : caseid
    implicit none
    type(MCCTRL), intent(in) :: ctrlMC
    integer :: i, j, k, l
    real(8) :: psum, pnrm, vsum, vol
    real(8), pointer, dimension(:,:,:) :: avg, std
    real(8), pointer, dimension(:) :: avg1d, std1d
    real(8) :: racc, raccm
    integer :: favg=5001, fstd=5002
    integer :: ncnt
    character*256 :: fn
    
    allocate(avg(nx,ny,nz), std(nx,ny,nz))
    allocate(avg1d(nx), std1d(nx))
    psum=0;vsum=0
    racc=1./iacc
    raccm=1./(iacc-1)
    ncnt=0
    do k=1, nz
        do j=1, ny
            do i=1, nx
                if (ctrlMC%lwrt(i,j) .and. sumpsi(i,j,k)>0) then
                    avg(i,j,k)=sumpsi(i,j,k)*racc
                    psum=psum+avg(i,j,k)
                    ! relative standard deviation
                    std(i,j,k)=sqrt((sqsumpsi(i,j,k)-sumpsi(i,j,k)**2*racc)*racc*raccm)/avg(i,j,k)
                    vol=GetVolPin(i,j,k)
                    avg(i,j,k)=avg(i,j,k)/vol
                    vsum=vsum+vol
                    ncnt=ncnt+1
                else
                    avg(i,j,k)=0.
                    std(i,j,k)=0.
                endif
            enddo
        enddo
    enddo
    pnrm=vsum/psum
    WRITE(fn,'(A,A,A)') 'MC_',TRIM(caseid),'.tly'
    open(favg, file=fn, status='unknown')    
    WRITE(fn,'(A,A,A)') 'MC_',TRIM(caseid),'.std'
    open(fstd, file=fn, status='unknown')    
    do k=1, nz
!        do j=ny, 1, -1
        do j=ctrlMC%lwrty(2), ctrlMC%lwrty(1), -1
            avg(:,j,k)=avg(:,j,k)*pnrm
            l=0
            do i=1, nx
                if (ctrlMC%lgap(i,j)) cycle
                l=l+1
                avg1d(l)=avg(i,j,k)
                std1d(l)=std(i,j,k)
            enddo
            if (l .eq. 0 ) cycle
            write(favg, '(2000(es12.5,1x))') avg1d(1:l)
            write(fstd, '(2000(es12.5,1x))') std1d(1:l)
            !write(favg, '(2000(es12.5,1x))') avg1d(ctrlMC.lwrtx(1):ctrlMC.lwrtx(2))
            !write(fstd, '(2000(es12.5,1x))') std1d(ctrlMC.lwrtx(1):ctrlMC.lwrtx(2))
            
        !    write(favg, '(2000(es12.5,1x))') avg(:,j,k)
        !    write(fstd, '(2000(es12.5,1x))') std(:,j,k)
        !    write(favg, '(2000(es12.5,1x))') avg(ctrlMC.lwrtx(1):ctrlMC.lwrtx(2),j,k)
        !    write(fstd, '(2000(es12.5,1x))') std(ctrlMC.lwrtx(1):ctrlMC.lwrtx(2),j,k)
        enddo
    enddo
    close(favg)
    close(fstd)
    deallocate(avg, std)
    deallocate(avg1d, std1d)
    
end subroutine
!
subroutine SpecTally(lmn, g, xs, trkl, tid)
    USE MCXsLib, ONLY : fxrphiomp
    use MCDefine
    implicit none
    integer, intent(in) :: lmn(3), g, tid
    real(8), intent(in) :: trkl
    integer :: ifxr, iz
    type(XsecSet), pointer :: xs
    
    ifxr=xs%ifxr;iz=xs%iz;
    fxrphiomp(ifxr,iz,g,tid)=fxrphiomp(ifxr,iz,g,tid)+trkl    
ENDSUBROUTINE


subroutine CollectSpecTally()
    USE MCXsLib, ONLY : fxrphi, fxrphiomp
    implicit none
    integer :: ifxr, iz, ig
    DO ifxr=1, nfxr
        do iz=1, nz
            DO ig=1, ng
                fxrphi(ifxr,iz,ig)=sum(fxrphiomp(ifxr,iz,ig,:))
                fxrphiomp(ifxr,iz,ig,:)=0._8;
            ENDDO
        ENDDO
    ENDDO
ENDSUBROUTINE

    
end module

