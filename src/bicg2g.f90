module bicg2g
! lu factors obtained from incomplete lu factorization
    use allocs
!    use geom,   only : nx,ny,nz,nxy
    implicit none

    integer, parameter :: ng2=2
    INTEGER, parameter :: indm24(2,2)=reshape((/1,3,2,4/),shape=(/2,2/))
    
    real(8)                            :: calpha,cbeta,crho,comega
    real(8),pointer,dimension(:,:)     :: del, au, ainvl, ainvu, ainvd !(4,nx)
    real(8),pointer,dimension(:,:,:)   :: delinv, al, deliau           !(4,nxy,nz)
    real(8),pointer,dimension(:,:,:)   :: vr,vr0,vp,vv,vs,vt,vy,vz     !(ng,nxy,nz)
    ! matrix variables
    integer,pointer                         :: nxs(:), nxe(:), nodel(:,:)
    real(8), pointer, dimension(:,:,:)      :: am           ! (4, nxy, nz) 2x2 diagonal block
    real(8), pointer, dimension(:,:,:,:)    :: ccr          ! (4,2, nxy, nz) coupling for radial
    real(8), pointer, dimension(:,:,:,:)    :: ccz          ! (2,2, nxy, nz) coupling for z-direction
    integer :: nx, ny, nz, nxy
    
    contains
    subroutine mallocbicg2g(nxy_, nx_, ny_, nz_)
        implicit none
        integer :: nxy_, nx_, ny_, nz_
        nx=nx_
        ny=ny_
        nz=nz_
        nxy=nxy_
        if (ASSOCIATED (del)) then
            deallocate(del)
            deallocate(au)
            deallocate(ainvl)
            deallocate(ainvu)
            deallocate(ainvd)
            deallocate(delinv)
            deallocate(al)
            deallocate(deliau)
            
            deallocate(vr)
            deallocate(vr0)
            deallocate(vp)
            deallocate(vv)
            deallocate(vs)
            deallocate(vt)
            deallocate(vr)
        endif
        call dmalloc(del,4,nx)
        call dmalloc(au,4,nx)
        call dmalloc(ainvl,4,nx)
        call dmalloc(ainvu,4,nx)
        call dmalloc(ainvd,4,nx)        
        call dmalloc(delinv,4,nxy,nz)
        call dmalloc(al,4,nxy,nz)
        call dmalloc(deliau,4,nxy,nz)

        call dmalloc(vr,ng2,nxy,nz)      
        call dmalloc(vr0,ng2,nxy,nz)
        call dmalloc(vp,ng2,nxy,nz)
        call dmalloc(vv,ng2,nxy,nz)
        call dmalloc(vs,ng2,nxy,nz)
        call dmalloc(vt,ng2,nxy,nz)
        call dmalloc0(vy,1,ng2,0,nxy,0,nz)
        call dmalloc0(vz,1,ng2,0,nxy,0,nz)      
    end subroutine
    
end module
    
module bicg2g_interface

    interface
    subroutine facilu2g(mt)
    use mtbl2g
    type(mb2g) :: mt
    end subroutine

    subroutine facilu1d2g(irow,k)
    integer :: irow, k
    end subroutine

    subroutine abi1d2g(irow,k)
    integer      :: irow, k
    end subroutine
    
    subroutine initbicg2g(mt, phi, rhs, r20)
    use mtbl2g
    type(mb2g) :: mt
    real(8)    :: rhs(:,:,:)
    real(8)    :: phi(:,:,:)
    real(8)    :: r20    
    end subroutine    
    
    subroutine solbicg2g(mt, r20, r2,phi)
    use mtbl2g
    type(mb2g) :: mt
    real(8)    :: r20
    real(8)    :: r2
    real(8)    :: phi(:,:,:)
    end subroutine
        
    subroutine minv2g(b,x)
! solve Mx=b for x given b
    real(8)        :: x(:,:,:),b(:,:,:) 
    end subroutine

    subroutine minvhex(b,x)
! solve Mx=b for x given b
    real(8)        :: x(:,:,:),b(:,:,:) 
    end subroutine

    
    subroutine sol1d2g(irow,k,b,x)
!  solve 1D problem using predetermined LU factors
    integer         :: irow,k
    real(8)         :: x(:,:),b(:,:)
    end subroutine
    
    subroutine sol2d2g(k,b,x)
    integer                 :: k
    real(8),pointer         :: b(:,:)
    real(8),pointer         :: x(:,:)
    end subroutine 

    subroutine sol2dhex(df,xlu,b,x)
    real(8)        :: df(:,:)
    real(8)        :: xlu(:,:,:)
    real(8)        :: b(:,:)
    real(8)        :: x(:,:)
    end subroutine 	
	   

!    subroutine solbicg2g(r20,r2)
!    real(8)  :: r20
!    real(8) :: r2
!    end subroutine

    end interface
    
end module