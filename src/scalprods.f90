module scalprods
    implicit none
    interface scalprod
        module procedure scalprod1g
        module procedure scalprodmg
    end interface
    
    contains
    function scalprod1g(nxy,nz,x,y)
        integer                 :: nxy, nz
        real(8),pointer         :: x(:,:),y(:,:)
        integer                 :: k, l, m
        real(8)                 :: scalprod1g
        
        scalprod1g=0
        do k=1,nz
        do l=1,nxy
            scalprod1g=scalprod1g+x(l,k)*y(l,k)
        enddo
        enddo

        return
    end function

    function scalprodmg(nxy,nz,ng,x,y)
        integer                 :: nxy, nz
        real(8),pointer        :: x(:,:,:),y(:,:,:)
        integer             :: ng
        integer             :: k, l, m
        real(8)                :: scalprodmg
        
        scalprodmg=0
        do k=1,nz
        do l=1,nxy
        do m=1,ng
            scalprodmg=scalprodmg+x(m,l,k)*y(m,l,k)
        enddo
        enddo
        enddo

        return
    end function    
end module