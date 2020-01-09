module mtbl2g
    integer, parameter                          :: ng2=2, nrdir2=4, nzdir2=2
    integer, private, dimension(4) :: LTDIR_INNER ! change direction order : bottom-top to top-bottom
    data LTDIR_INNER  /1, 2, 4, 3/ 
    integer, parameter :: indm24(2,2)=reshape((/1,3,2,4/),shape=(/2,2/))

    type mb2g
        real(8), pointer, dimension(:,:,:)      :: am           ! (4, nxy, nz) 2x2 diagonal block
        real(8), pointer, dimension(:,:,:,:)    :: ccr          ! (4,2, nxy, nz) coupling for radial
        real(8), pointer, dimension(:,:,:,:)    :: ccz          ! (2,2, nxy, nz) coupling for z-direction
        integer, pointer, dimension(:,:,:)      :: iccr         ! (4,nxy,nz) coupling index for radial direction
        integer, pointer, dimension(:,:,:)      :: iccz         ! (2,nxy,nz) coupling index for axial direction
        integer, pointer, dimension(:)          :: nxs, nxe     ! starting and ending x-coordinates of the core region for each y-direction
        integer, pointer, dimension(:,:)        :: nodel        ! radial node index corresponding to the x- and y- cordinate
        integer :: nx, ny, nxy, nz
    contains 
        procedure :: settozero
        procedure :: newmb2g
        procedure :: copymb2g
        procedure :: delmb2g
        procedure :: setnxse
        procedure :: setnodel
        procedure :: setmb2g
        procedure :: accmb2g
        procedure :: axb2g
    end type
    
contains    

subroutine settozero(mt)
    implicit none
    class(mb2g) :: mt
    mt%am = 0._8
    mt%ccr = 0._8
    mt%ccz = 0._8
    mt%iccr = 0
    mt%iccz = 0
    mt%nxs = 0
    mt%nxe = 0
    mt%nodel = 0
end subroutine

subroutine newmb2g(mt, nxy, nx, ny, nz)
    use allocs
    implicit none
    class(mb2g) :: mt
    integer, intent(in)         :: nx, ny, nz
    integer                     :: i, nxy
    ! caution
    ! nxy != nx * ny
    ! nxy means the # of real mesh  nxy <= nx * ny 
    call dmalloc(mt%am,4,nxy,nz)
    call dmalloc(mt%ccr,4,2,nxy,nz)
    call dmalloc(mt%ccz,2,2,nxy,nz)
    call dmalloc(mt%iccr,4,nxy,nz)
    call dmalloc(mt%iccz,2,nxy,nz)
    call dmalloc(mt%nxs,ny)
    call dmalloc(mt%nxe,ny)
    call dmalloc0(mt%nodel,0,nx+1,0,ny+1)
    mt%nx=nx
    mt%ny=ny
    mt%nxy=nxy
    mt%nz=nz
endsubroutine 

subroutine copymb2g(mt,mt1)
    implicit none
    class(mb2g) :: mt
    type(mb2g), intent(in) :: mt1
    integer :: ij, z, j
    if (mt1%nz .ne. mt%nz .or. mt1%nxy .ne. mt%nxy) return
    do z=1, mt1%nz
        do ij=1, mt1%nxy
            mt%am(:,ij,z)=mt1%am(:,ij,z)
            mt%ccr(:,:,ij,z)=mt1%ccr(:,:,ij,z)
            mt%ccz(:,:,ij,z)=mt1%ccz(:,:,ij,z)
            mt%iccr(:,ij,z)=mt1%iccr(:,ij,z)
            mt%iccz(:,ij,z)=mt1%iccz(:,ij,z)
        enddo
    enddo
    
    mt%nodel=mt1%nodel
    mt%nxs=mt1%nxs
    mt%nxe=mt1%nxe
endsubroutine   

subroutine delmb2g(mt)
    implicit none
    class(mb2g) :: mt
    deallocate(mt%am)
    deallocate(mt%ccr)
    deallocate(mt%ccz)
endsubroutine

! set an element
subroutine setnxse(mt, j, s, e)
    implicit none
    class(mb2g) :: mt
    integer, intent(in) :: j, s, e
    mt%nxs(j)=s
    mt%nxe(j)=e
endsubroutine

subroutine setnodel(mt, i, j, l)
    implicit none
    class(mb2g) :: mt
    integer, intent(in) :: i, j, l
    mt%nodel(i,j) =l
endsubroutine

subroutine setmb2g(mt,l,z,dir, g, value, neigh)
    implicit none
    class(mb2g) :: mt
    integer, intent(in) :: l        ! row in matrix
    integer, intent(in) :: z        ! index for plane
    integer, intent(in) :: dir      ! W E N S T B
    integer, intent(in) :: g        ! g11=1 g12=2 g21=3 g22=4 g1=1 g2=2
    real(8), intent(in) :: value
    integer, intent(in), optional :: neigh   
    
    if(dir .eq. 0) then
        mt%am(g,l,z)=value
    elseif(dir<5) then
        mt%ccr(LTDIR_INNER(dir),g,l,z)=value
        if (PRESENT(neigh)) mt%iccr(LTDIR_INNER(dir),l,z)=neigh
    else
        mt%ccz(dir-4,g,l,z)=value
        if (PRESENT(neigh)) mt%iccz(dir-4,l,z)=neigh
    endif
endsubroutine 

! accumulate an element
subroutine accmb2g(mt,l,z,dir,g,value,neigh) 
    implicit none
    class(mb2g) :: mt
    integer, intent(in) :: l        ! row in matrix
    integer, intent(in) :: z        ! index for plane
    integer, intent(in) :: dir    ! N S W E T B
    integer, intent(in) :: g        ! g11=1 g12=2 g21=3 g22=4 g1=1 g2=2
    real(8), intent(in) :: value    
    integer, intent(in), optional   :: neigh    
    
    if(dir .eq. 0) then
        mt%am(g,l,z)=mt%am(g,l,z)+value
    elseif(dir<5) then
        mt%ccr(LTDIR_INNER(dir),g,l,z)=mt%ccr(LTDIR_INNER(dir),g,l,z)+value
        if (PRESENT(neigh)) mt%iccr(LTDIR_INNER(dir),l,z)=neigh
    else
        mt%ccz(dir-4,g,l,z)=mt%ccz(dir-4,g,l,z)+value
        if (PRESENT(neigh)) mt%iccz(dir-4,l,z)=neigh
    endif
endsubroutine 

subroutine axb2g(mt,phi,aphi)
!    use const
    implicit none
    class(mb2g)                 :: mt
    real(8), pointer            :: phi(:,:,:)
    real(8), pointer            :: aphi(:,:,:)
    integer                     :: k, l, m, m1, m2, idir, idirz
    real(8)                     :: aphil
    real(8), pointer, dimension(:,:,:)      :: am       
    real(8), pointer, dimension(:,:,:,:)    :: ccr      
    real(8), pointer, dimension(:,:,:,:)    :: ccz      
    integer, pointer, dimension(:,:,:)      :: iccr     
    integer, pointer, dimension(:,:,:)      :: iccz   
    integer :: nz, nxy  
    
    am=>mt%am;ccr=>mt%ccr;ccz=>mt%ccz;iccr=>mt%iccr;iccz=>mt%iccz
    nz=mt%nz;nxy=mt%nxy
    do k=1,nz
        do l=1,nxy
            do m=1,ng2  ! 1 to 2
                m1=indm24(m,1)
                m2=indm24(m,2)
                !am(:,l,k) is 2 x 2 matrix and phi(:,l,k) is 2 x 1.
                aphil=am(m1,l,k)*phi(1,l,k)+am(m2,l,k)*phi(2,l,k)
                do idir=1,nrdir2
                    aphil=aphil+ccr(idir,m,l,k)*phi(m,iccr(idir,l,k),k)
                enddo    
                do idir=1,nzdir2
                    aphil=aphil+ccz(idir,m,l,k)*phi(m,l, iccz(idir,l,k))
                enddo
                           
                aphi(m,l,k)=aphil
            enddo
        enddo
    enddo
endsubroutine

endmodule