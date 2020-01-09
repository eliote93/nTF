subroutine sol2d2g(k,b,x)
! solve a 2d problem using precalculated LU factors
    use bicg2g
    use bicg2g_interface, only : sol1d2g
!    use geom, only : nxs,nxe,nodel
!    use cmfd2g, only : ccr
    implicit none

    integer                 :: k
    real(8),pointer         :: b(:,:)
    real(8),pointer         :: x(:,:) !dmm

    integer                 :: i,j,jp1,lout,l,ls
    real(8)                 :: s1dl(ng2,nx), b01d(ng2,nx)

    j=1

    do i=1,nx !dmm
        s1dl(1,i)=0
        s1dl(2,i)=0
    enddo

!  forward solve
    lout=0
    b01d=0
    do j=1,ny
!        l=lout
        do i=nxs(j),nxe(j)
            l=nodel(i,j)
            b01d(1,i)=b(1,l)-ccr(3,1,l,k)*s1dl(1,i)
            b01d(2,i)=b(2,l)-ccr(3,2,l,k)*s1dl(2,i)
        enddo
        call sol1d2g(j,k,b01d,s1dl)

!        l=lout
        do i=nxs(j),nxe(j)
            l=nodel(i,j)
            x(1,l)=s1dl(1,i)
            x(2,l)=s1dl(2,i)
        enddo
!        lout=lout+(nxe(j)-nxs(j)+1)
    enddo

!  backward solve
!    lout=nxy-(nxe(ny)-nxs(ny)+1)
    jp1=ny
    do j=ny-1,1,-1
!        l=lout
        do i=nxe(j),nxs(j),-1
            l=nodel(i,j)
            ls=nodel(i,jp1)
            b01d(1,i)=x(1,ls)*ccr(4,1,l,k)
            b01d(2,i)=x(2,ls)*ccr(4,2,l,k)
!            l=l-1
        enddo

        call sol1d2g(j,k,b01d,s1dl)


!        l=lout
        do i=nxe(j),nxs(j),-1
            l=nodel(i,j)
            x(1,l)=x(1,l)-s1dl(1,i)
            x(2,l)=x(2,l)-s1dl(2,i)
!            l=l-1
        enddo
!        lout=lout-(nxe(j)-nxs(j)+1)
        jp1=j
    enddo

    return
end subroutine
