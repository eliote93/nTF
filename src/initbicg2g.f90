subroutine initbicg2g(mt,phi,rhs,r20)
    use bicg2g
!    use cmfd2g, only : axb2g,axbhex
!MZ
    use mtbl2g, only : mb2g, axb2g
!!
!    use geom,   only : nxy,nz
!	use sfam_cntl, only : ifrect
    implicit none

    real(8),pointer            :: rhs(:,:,:)
    real(8),pointer            :: phi(:,:,:)
    real(8)                    :: r20

    integer                    :: l,k
    real(8)                    :: b2t
    real(8)                    :: aphi(ng2,nxy,nz)

!MZ
    type(mb2g)                              :: mt
!!

    calpha=1.
    crho=1.
    comega=1.


!	if(ifrect) then
	call axb2g(mt,phi,vp)
!	else
!		call axbhex(phi,aphi)
!	endif
    r20=0
    b2t=0

    do k=1,nz
    do l=1,nxy
        vr(1,l,k)=rhs(1,l,k)-vp(1,l,k)
        vr(2,l,k)=rhs(2,l,k)-vp(2,l,k)
        vr0(1,l,k)=vr(1,l,k)
        vr0(2,l,k)=vr(2,l,k)
        vp(1,l,k)=0
        vp(2,l,k)=0
        vv(1,l,k)=0
        vv(2,l,k)=0
        r20=r20+vr(1,l,k)*vr(1,l,k)+vr(2,l,k)*vr(2,l,k)
        b2t=b2t+rhs(1,l,k)*rhs(1,l,k)+rhs(2,l,k)*rhs(2,l,k)
    enddo
    enddo
    r20=sqrt(r20)

    return
end subroutine
