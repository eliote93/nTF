subroutine solbicg2g(mt, r20, r2,phi)
    use bicg2g
    use bicg2g_interface, only : minv2g
!    use cmfd2g,     only : axb2g,axbhex
!MZ
    use mtbl2g,     only : mb2g, axb2g
!!
    use scalprods,  only : scalprod
!	use sfam_cntl,  only : ifrect
    implicit none

    real(8), intent(in)         :: r20
    real(8), intent(out)        :: r2
    real(8), pointer            :: phi(:,:,:)

    real(8)                     :: crhod, r0v, pts, ptt
	integer                     :: l, k
!MZ
    type(mb2g),intent(in)       :: mt
!!

! solves the linear system by preconditioned BiCGSTAB Algorithm
    crhod=crho
    crho=scalprod(nxy,nz,ng2,vr0,vr)
    cbeta=crho*calpha/(crhod*comega)
!    print *, crho, crhod, cbeta
   do k=1,nz
       do l=1,nxy
            vp(1,l,k)=vr(1,l,k)+cbeta*(vp(1,l,k)-comega*vv(1,l,k))
            vp(2,l,k)=vr(2,l,k)+cbeta*(vp(2,l,k)-comega*vv(2,l,k))
        enddo
    enddo
!    print *, sum(vp), sum(vy)
!	if(ifrect) then
		call minv2g(vp,vy)
!    print *, sum(vp), sum(vy)
		call axb2g(mt,vy,vv)
!	else
!		call minvhex(vp,vy)
!		call axbhex(vy,vv)
!	endif

    r0v=scalprod(nxy,nz,ng2,vr0,vv)
    calpha=crho/r0v

    do k=1,nz
    do l=1,nxy
        vs(1,l,k)=vr(1,l,k)-calpha*vv(1,l,k)
        vs(2,l,k)=vr(2,l,k)-calpha*vv(2,l,k)
    enddo
    enddo

!	if(ifrect) then
		call minv2g(vs,vz)
		call axb2g(mt,vz,vt)
!	else
!		call minvhex(vs,vz)
!		call axbhex(vz,vt)
!	endif

    pts=scalprod(nxy,nz,ng2,vs,vt)
    ptt=scalprod(nxy,nz,ng2,vt,vt)

    comega=0
    if(ptt .ne. 0) comega=pts/ptt

    r2=0
    do k=1,nz
    do l=1,nxy
        phi(1,l,k)=phi(1,l,k)+calpha*vy(1,l,k)+comega*vz(1,l,k)
        phi(2,l,k)=phi(2,l,k)+calpha*vy(2,l,k)+comega*vz(2,l,k)
        vr(1,l,k)=vs(1,l,k)-comega*vt(1,l,k)
        vr(2,l,k)=vs(2,l,k)-comega*vt(2,l,k)
        r2=r2+vr(1,l,k)*vr(1,l,k)+vr(2,l,k)*vr(2,l,k)
    enddo
    enddo
    r2=sqrt(r2)/r20

    return
end subroutine
