subroutine facilu1d2g(irow,k)
! incomplete factorization of block-tridiagonal matices
    use bicg2g
!    use geom,   only : nxs,nxe,nodel
    implicit none
    
    integer      :: irow,k

    integer                 :: i,l,im1,lm1
    real(8)                 :: f,ald1,ald2,ald3,ald4
    
    i=nxs(irow)
    l=nodel(i,irow)

! first column, inverse of a 2x2 matrix
    f=1/(del(1,i)*del(4,i)-del(2,i)*del(3,i))
    delinv(1,l,k)=del(4,i)*f
    delinv(2,l,k)=-del(2,i)*f
    delinv(3,l,k)=-del(3,i)*f
    delinv(4,l,k)=del(1,i)*f

!   calc. inv(del)*u for later use in backsub
    deliau(1,l,k)=delinv(1,l,k)*au(1,i)+delinv(2,l,k)*au(3,i)
    deliau(2,l,k)=delinv(1,l,k)*au(2,i)+delinv(2,l,k)*au(4,i)
    deliau(3,l,k)=delinv(3,l,k)*au(1,i)+delinv(4,l,k)*au(3,i)
    deliau(4,l,k)=delinv(3,l,k)*au(2,i)+delinv(4,l,k)*au(4,i)

    im1=i
    do i=nxs(irow)+1,nxe(irow)
        lm1=l; l=l+1
        
        ald1=al(1,l,k)*delinv(1,lm1,k)+al(2,l,k)*delinv(3,lm1,k)
        ald2=al(1,l,k)*delinv(2,lm1,k)+al(2,l,k)*delinv(4,lm1,k)
        ald3=al(3,l,k)*delinv(1,lm1,k)+al(4,l,k)*delinv(3,lm1,k)
        ald4=al(3,l,k)*delinv(2,lm1,k)+al(4,l,k)*delinv(4,lm1,k)
        
        del(1,i)=del(1,i)-ald1*au(1,im1)-ald2*au(3,im1)
        del(2,i)=del(2,i)-ald1*au(2,im1)-ald2*au(4,im1)
        del(3,i)=del(3,i)-ald3*au(1,im1)-ald4*au(3,im1)
        del(4,i)=del(4,i)-ald3*au(2,im1)-ald4*au(4,im1)
        
        f=1._8/(del(1,i)*del(4,i)-del(2,i)*del(3,i))
        delinv(1,l,k)=del(4,i)*f
        delinv(2,l,k)=-del(2,i)*f
        delinv(3,l,k)=-del(3,i)*f
        delinv(4,l,k)=del(1,i)*f
        
!   calc. inv(del)*u for later use in backsub
        deliau(1,l,k)=delinv(1,l,k)*au(1,i)+delinv(2,l,k)*au(3,i)
        deliau(2,l,k)=delinv(1,l,k)*au(2,i)+delinv(2,l,k)*au(4,i)
        deliau(3,l,k)=delinv(3,l,k)*au(1,i)+delinv(4,l,k)*au(3,i)
        deliau(4,l,k)=delinv(3,l,k)*au(2,i)+delinv(4,l,k)*au(4,i)
        im1=i
    enddo

    return
end subroutine