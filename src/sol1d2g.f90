subroutine sol1d2g(irow,k,b,x)
!  solve 1D problem using predetermined LU factors
    use bicg2g
!    use geom, only : nxs,nxe,nodel
    implicit none

    integer                 :: irow,k
    real(8),pointer            :: b(:,:)
    real(8),pointer            :: x(:,:)
    
    integer                 :: ibeg,iend,l,i,im1,ip1
    real(8)                 :: y(ng2,nx),b1i,b2i

    ibeg=nxs(irow)
    iend=nxe(irow)
    l=nodel(ibeg,irow)

!  forward substitution
    i=ibeg
    y(1,i)=delinv(1,l,k)*b(1,i)+delinv(2,l,k)*b(2,i)
    y(2,i)=delinv(3,l,k)*b(1,i)+delinv(4,l,k)*b(2,i)
    im1=i
    do i=ibeg+1,iend
        l=l+1
        b1i=b(1,i)-al(1,l,k)*y(1,im1)-al(2,l,k)*y(2,im1)
        b2i=b(2,i)-al(3,l,k)*y(1,im1)-al(4,l,k)*y(2,im1)
        y(1,i)=delinv(1,l,k)*b1i+delinv(2,l,k)*b2i
        y(2,i)=delinv(3,l,k)*b1i+delinv(4,l,k)*b2i
        im1=i
    enddo

!  backward substitution
    x(1,iend)=y(1,iend)
    x(2,iend)=y(2,iend)
    ip1=iend
    do i=iend-1,ibeg,-1
        l=l-1
        x(1,i)=y(1,i)-deliau(1,l,k)*x(1,ip1)-deliau(2,l,k)*x(2,ip1)
        x(2,i)=y(2,i)-deliau(3,l,k)*x(1,ip1)-deliau(4,l,k)*x(2,ip1)
        ip1=i
    enddo

    return
end subroutine