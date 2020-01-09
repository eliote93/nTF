subroutine abi1d2g(irow,k)
! approximate block inverse from the LU factors
    use bicg2g
!    use geom,   only : nxs,nxe,nodel
    implicit none
    integer                 :: irow, k
    integer                 :: lp1,mp1,i,ix,l
    real(8)                 :: al1,al2,al3,al4,au1,au2,au3,au4


    ix=nxe(irow)
    l=nodel(ix,irow)

    do i=1,4
        ainvd(i,ix)=delinv(i,l,k)
    enddo

    do i=nxe(irow)-1,nxs(irow),-1
        lp1=l
        l=l-1
        mp1=ix
        ix=ix-1
!     lower part of the inverse
        al1=ainvd(1,mp1)*al(1,lp1,k)+ainvd(2,mp1)*al(3,lp1,k)
        al2=ainvd(1,mp1)*al(2,lp1,k)+ainvd(2,mp1)*al(4,lp1,k)
        al3=ainvd(3,mp1)*al(1,lp1,k)+ainvd(4,mp1)*al(3,lp1,k)
        al4=ainvd(3,mp1)*al(2,lp1,k)+ainvd(4,mp1)*al(4,lp1,k)
        ainvl(1,mp1)=-al1*delinv(1,l,k)-al2*delinv(3,l,k)
        ainvl(2,mp1)=-al1*delinv(2,l,k)-al2*delinv(4,l,k)
        ainvl(3,mp1)=-al3*delinv(1,l,k)-al4*delinv(3,l,k)
        ainvl(4,mp1)=-al3*delinv(2,l,k)-al4*delinv(4,l,k)
!     upper part of the inverse
        au1=delinv(1,l,k)*au(1,ix)+delinv(2,l,k)*au(3,ix)
        au2=delinv(1,l,k)*au(2,ix)+delinv(2,l,k)*au(4,ix)
        au3=delinv(3,l,k)*au(1,ix)+delinv(4,l,k)*au(3,ix)
        au4=delinv(3,l,k)*au(2,ix)+delinv(4,l,k)*au(4,ix)
        ainvu(1,ix)=-au1*ainvd(1,mp1)-au2*ainvd(3,mp1)
        ainvu(2,ix)=-au1*ainvd(2,mp1)-au2*ainvd(4,mp1)
        ainvu(3,ix)=-au3*ainvd(1,mp1)-au4*ainvd(3,mp1)
        ainvu(4,ix)=-au3*ainvd(2,mp1)-au4*ainvd(4,mp1)
!     diagonal part
        ainvd(1,ix)=delinv(1,l,k)-au1*ainvl(1,mp1)-au2*ainvl(3,mp1)
        ainvd(2,ix)=delinv(2,l,k)-au1*ainvl(2,mp1)-au2*ainvl(4,mp1)
        ainvd(3,ix)=delinv(3,l,k)-au3*ainvl(1,mp1)-au4*ainvl(3,mp1)
        ainvd(4,ix)=delinv(4,l,k)-au3*ainvl(2,mp1)-au4*ainvl(4,mp1)
    enddo
return
end subroutine
