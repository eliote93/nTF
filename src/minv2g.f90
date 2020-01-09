subroutine minv2g(b,x)
! solve Mx=b for x given b
    use allocs
    use bicg2g
    use bicg2g_interface, only : sol2d2g
!    use geom,   only : nxy,nx,ny,nz
!    use cmfd2g, only : ccz
    implicit none

    real(8),pointer         :: b(:,:,:)
    real(8),pointer         :: x(:,:,:)
    integer                 :: l,k,kp1
    real(8),pointer, save   :: s(:,:), b0(:,:)
    logical,save            :: first=.true.


    if (first) then
        first = .false.
        call dmalloc0(s, 1, ng2, -nx+1, nxy+nx)
        call dmalloc(b0,ng2,nxy)
    endif
    s = 0


    if (nz .eq. 1) then
        ! forward solve ONLY
        do l=1,nxy
            b0(1,l)=b(1,l,1)
            b0(2,l)=b(2,l,1)
        enddo
        call sol2d2g(1,b0,s)
        do l=1,nxy
            x(1,l,1)=s(1,l)
            x(2,l,1)=s(2,l)
        enddo
!            b0=>b(1:,1:,1)
!            s=>x(1:,1:,1)
!        call sol2d2g(1,b0,s)
    else
        ! forward solve
        do k=1,nz
            do l=1,nxy
                b0(1,l)=b(1,l,k)-ccz(1,1,l,k)*s(1,l)
                b0(2,l)=b(2,l,k)-ccz(1,2,l,k)*s(2,l)
            enddo

            call sol2d2g(k,b0,s)

            do l=1,nxy
                x(1,l,k)=s(1,l)
                x(2,l,k)=s(2,l)
            enddo
        enddo

        ! backward solve
        do k=nz-1,1,-1
            kp1=k+1
            do l=1,nxy
                b0(1,l)=x(1,l,kp1)*ccz(2,1,l,k)
                b0(2,l)=x(2,l,kp1)*ccz(2,2,l,k)
            enddo

            call sol2d2g(k,b0,s)

            do l=1,nxy
                x(1,l,k)=x(1,l,k)-s(1,l)
                x(2,l,k)=x(2,l,k)-s(2,l)
            enddo
        enddo
    endif
    return
end subroutine
