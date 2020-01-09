subroutine facilu2g(mt)
! perform incomplete LU factorization for the 2D coefficient matrices
    use bicg2g
!    use geom,   only : nxs, nxe, nodel    
!    use cmfd2g, only : am, ccr
!MZ
    use mtbl2g,     only : mb2g
!!    
    implicit none
    real(8), parameter :: zero=0.
    
    integer                 :: m, i, j, k, l, ln, jm1, lnm1, lnp1    
!MZ
    type(mb2g)                              :: mt
    ! set geometric module variables
    nxs=>mt%nxs;nxe=>mt%nxe;nodel=>mt%nodel
    ! set matrix variables
    am=>mt%am;ccr=>mt%ccr;ccz=>mt%ccz
!!
    ainvl=0._8; ainvu=0._8
    
! loop over planes
    do k=1,nz
! first row
        j=1; l=0
        do i=nxs(j),nxe(j)
!            l=l+1
            l=nodel(i,j)
            del(1,i)=am(indm24(1,1),l,k)
            del(2,i)=am(indm24(1,2),l,k)
            del(3,i)=am(indm24(2,1),l,k)
            del(4,i)=am(indm24(2,2),l,k)

            al(1,l,k)=ccr(1,1,l,k)
            al(2,l,k)=zero
            al(3,l,k)=zero
            al(4,l,k)=ccr(1,2,l,k)

            au(1,i)=ccr(2,1,l,k)
            au(2,i)=zero
            au(3,i)=zero
            au(4,i)=ccr(2,2,l,k)
        enddo
! loop over rows
        do j=2,ny
            jm1=j-1
!           obtain incomplete lu factor for the 1d matrix of the row
            call facilu1d2g(jm1,k)
!           obtain the inverse of the 1d matrix
            call abi1d2g(jm1,k)

!           d_j+1 = a_j+1 - l_j+1*jnv(d_j)*u_j
            do i=nxs(j),nxe(j)
                l=nodel(i,j)
                
                if (i.ge.nxs(jm1) .and. i.le.nxe(jm1)) then ! coupling of previous row
                    ln=nodel(i,jm1)
                    del(1,i)=  am(indm24(1,1),l,k)-ccr(3,1,l,k)*ainvd(1,i)*ccr(4,1,ln,k)
                    del(2,i)=  am(indm24(1,2),l,k)-ccr(3,1,l,k)*ainvd(2,i)*ccr(4,2,ln,k)
                    del(3,i)=  am(indm24(2,1),l,k)-ccr(3,2,l,k)*ainvd(3,i)*ccr(4,1,ln,k)
                    del(4,i)=  am(indm24(2,2),l,k)-ccr(3,2,l,k)*ainvd(4,i)*ccr(4,2,ln,k)
                else
                    del(1,i)=  am(indm24(1,1),l,k)
                    del(2,i)=  am(indm24(1,2),l,k)
                    del(3,i)=  am(indm24(2,1),l,k)
                    del(4,i)=  am(indm24(2,2),l,k)
                endif

                if(i.ne.nxs(j) .and. (i-1).ge.nxs(jm1) .and. (i-1).le.nxe(jm1)) then
                    lnm1=nodel(i-1,jm1)
                    al(1,l,k)=ccr(1,1,l,k)-ccr(3,1,l,k)*ainvl(1,i)*ccr(4,1,lnm1,k)
                    al(2,l,k)=            -ccr(3,1,l,k)*ainvl(2,i)*ccr(4,2,lnm1,k)
                    al(3,l,k)=            -ccr(3,2,l,k)*ainvl(3,i)*ccr(4,1,lnm1,k)
                    al(4,l,k)=ccr(1,2,l,k)-ccr(3,2,l,k)*ainvl(4,i)*ccr(4,2,lnm1,k)
                else
                    al(1,l,k)=ccr(1,1,l,k)
                    al(2,l,k)=zero
                    al(3,l,k)=zero
                    al(4,l,k)=ccr(1,2,l,k)
                endif
                
                if(i.ne.nxe(j) .and. (i+1).ge.nxs(jm1) .and. (i+1).le.nxe(jm1)) then
                    lnp1=nodel(i+1,jm1)
                    au(1,i)=ccr(2,1,l,k)-ccr(3,1,l,k)*ainvu(1,i)*ccr(4,1,lnp1,k)
                    au(2,i)=            -ccr(3,1,l,k)*ainvu(2,i)*ccr(4,2,lnp1,k)
                    au(3,i)=            -ccr(3,2,l,k)*ainvu(3,i)*ccr(4,1,lnp1,k)
                    au(4,i)=ccr(2,2,l,k)-ccr(3,2,l,k)*ainvu(4,i)*ccr(4,2,lnp1,k)
                else
                    au(1,i)=ccr(2,1,l,k)
                    au(2,i)=zero
                    au(3,i)=zero
                    au(4,i)=ccr(2,2,l,k)
                endif
            enddo ! i
        enddo ! j
! obtain incomplete lu factor for the 1d matrix of the last row
        call facilu1d2g(ny,k)
    enddo ! k, plane 

    return
end subroutine