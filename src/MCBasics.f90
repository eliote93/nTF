module MCBasics
    use MCDefine
    use PARAM, only : NORTH, EAST, SOUTH, WEST, TOP, BOTTOM
    private
    public :: InitMCBasics, AdrFSD, AdrPPow, DTSRect, SMPDTC, Collision, FindRing, DTCL, SMPNeutron
    public :: GetNQueueCur, GetMCK, ShannonEntropy, GetFSNxt, GetFS, SmpGroupiso, SMPGroup, GetSurfOP
    public :: SMPDirIso, Move, TrackFSD, Tally, Track_Keff, TallyTurnOn, TallyTurnOff
    public :: BankNeutron, BankNeutronIso, ResetKeffTally, SetDirRefl, Move2SurfRect, ResetFSD, ResetQueue
    public :: ReflCent, CollectTallies
    real(8), private :: trk_gk, trk_gkomp(16)
    real(8), private, pointer, dimension(:,:,:) :: FSD
    real(8), private, pointer, dimension(:,:,:,:) :: FSDOMP
    real(8), private, pointer, dimension(:,:,:) :: ppow
    real(8), private, pointer, dimension(:,:,:,:) :: ppowomp
    type(nfs), pointer, private, dimension(:,:,:) :: nqueue
    type(nfs), pointer, private, dimension(:,:) :: nqcur, nqnxt
    integer, private :: nthrd, nqompcur(16), nqompnxt(16), iqomp(16)
    integer, pointer, dimension(:,:,:) :: nseqomp
    integer, private :: nqsum(0:16)
    integer, private :: nqidx, icur, inxt
    integer, private :: ng, nqmax, nqmaxloc, nht, ndim, ncur, nnxt, nx,ny,nz, nchi
    logical, private :: mobility(3)
    logical, private :: lic=.true.   !--???
    integer, private, parameter :: COLABSR=1, COLSCT=2
    integer, private, parameter :: LTNQ(3)=(/1,2,1/)
    real(8), private :: MCkeff
    logical, private :: ltally
    logical, private :: lomp
    integer, private, parameter, dimension(3,4,2) :: LTREFL=RESHAPE((/ 1, -1, 1, &
                                                                       -1, 1, 1, &
                                                                       1, -1, 1, &
                                                                       -1, 1, 1, &
                                                                       1, 1, -1, &
                                                                       1, 1, -1, &
                                                                       1, 1, 1,  & ! not used from here
                                                                       1, 1, 1/), SHAPE(LTREFL))


contains

subroutine InitMCBasics(nhtMC, nxyzMC, albedo, lomp_, nthrd_)
    use MCMap, only : GetNBaseFS, GetBaseFS
    use BenchXs, only : ngben
    USE Core_mod,   ONLY: GroupInfo
    use rng
    implicit none
    integer, intent(in) :: nhtMC, nxyzMC(3), nthrd_
    real(8), intent(in) :: albedo(10)
    logical, intent(in) :: lomp_
    integer :: nbase, fs, i, ir, ithrd, iloc

    nthrd=nthrd_
    nht=nhtMC
    ng=ngben
    nchi=ngben
    if( ng .EQ. 0 )THEN !lxslib
        ng=groupInfo%ng
        nchi=GroupInfo%nchi
    ENDIF
    nqmax=nht*5
    nqmaxloc=nqmax/nthrd
!    allocate(nqueue(nqmax,2))
    allocate(nqueue(nqmaxloc,nthrd,2))
    nqidx=1
    icur=LTNQ(nqidx)
    inxt=LTNQ(nqidx+1)
    nqcur=>nqueue(:,:,icur)
    nqnxt=>nqueue(:,:,inxt)
    lomp=lomp_

    nbase=GetNBaseFS()
    ncur=0
    iloc=0
    ithrd=1
    nqompcur=0
    do i=1, nht
        iloc=iloc+1
        ncur=ncur+1
        if (iloc>nqmaxloc) then
            iloc=1
            ithrd=ithrd+1
        endif
        if (i>nbase) then
            ir=nbase*getRN()+1
        else
            ir=i
        endif
        call CopyNFS(nqcur(iloc, ithrd), GetBaseFS(ir))
        nqompcur(ithrd)=nqompcur(ithrd)+1
    enddo
    nqsum(1)=nqompcur(1)
    do i=2, nthrd
        nqsum(i)=nqsum(i-1)+nqompcur(i)
    enddo


    nx=nxyzMC(1); ny=nxyzMC(2); nz=nxyzMC(3)
    allocate(FSD(nx,ny,nz))
!$  allocate(FSDOMP(nx,ny,nz,nthrd))
    allocate(ppow(nx,ny,nz))
!$  allocate(ppowomp(nx,ny,nz,nthrd))
    call ResetFSD()
    call ResetKeffTally()

    mobility(1)=.true.
    mobility(2)=.true.
    mobility(3)=.true.
    if (nz .eq. 1 .and. albedo(5) .eq. 0 .and. albedo(6) .eq. 0) mobility(3)=.false.

    MCkeff=1.0
    ltally=.false.
end subroutine

subroutine AdrFSD(pFSD)
    real(8), pointer, dimension(:,:,:) :: pFSD
    pFSD=>FSD
end subroutine

subroutine AdrPPow(pPPow)
    real(8), pointer, dimension(:,:,:) :: pPPow
     pPPow=>ppow
end subroutine

subroutine TallyTurnOn()
    ltally=.true.
end subroutine

subroutine TallyTurnOff()
    ltally=.false.
end subroutine

function GetSurfOp(surf, axis) result(ret)
    integer :: ret
    integer, intent(in) :: surf, axis
    ret=SurfOp(surf,axis)
end function

subroutine SetDirRefl(nst)
    implicit none
    type(nstat), intent(inout) :: nst
    nst%dir=nst%dir*LTREFL(:, nst%surf, nst%axis)
end subroutine

function GetMCK() result(keff)
    implicit none
    real(8) :: keff
!$  if (lomp) trk_gk=sum(trk_gkomp)
    keff=trk_gk/nht
    MCkeff=keff
end function

function GetNQueueCur()
    implicit none
    integer :: GetNQueueCur
    GetNQueueCur=ncur
end function

function GetNQueueNxt()
    implicit none
    integer :: GetNQueueNxt
    GetNQueueNxt=nnxt
end function

function GetFS(idx)
    implicit none
    type(nfs) :: GetFS
    integer, intent(in) :: idx
    GetFS=nqcur(idx,1)
end function

function GetFSNxt(idx)
    implicit none
    type(nfs) :: GetFSNxt
    integer, intent(in) :: idx
    integer :: ithrd, nsum, in
    do ithrd=1, nthrd-1
        if (nqsum(ithrd)>=idx) exit
    enddo
    in=idx-nqsum(ithrd-1)
    GetFSNxt=nqcur(in, ithrd)
end function

subroutine ResetFSD()
    implicit none
    integer :: t, k, j
!$  if (lomp) then
!$      do t=1, nthrd
!$          do k=1, nz
!$              do j=1, ny
!$                  FSDOMP(:,j,k,t)=0.
!$                  ppowomp(:,j,k,t)=0.
!$              enddo
!$          enddo
!$      enddo
!$  else
        do k=1, nz
            do j=1, ny
                FSD(:,j,k)=0.
                ppow(:,j,k)=0.
            enddo
        enddo
!$  endif
end subroutine


!subroutine GenInitialSource()
!    use RNG
!    use MCDefine
!    implicit none
!    integer :: i, j, nadd, icp
!    nadd=nht-ncur
!
!    do i=1, nadd
!        icp=GetRn()*ncur+1
!        call CopyNFS(nqcur(ncur+i), nqcur(icp))
!    end do
!!$  if (lomp) then
!!$      nqompcur=0
!!$      nqompcur(1)=nht/nthrd
!!$      nseqsum(1)=nht/nthrd
!!$      do i=2, nthrd-1
!!$          nqompcur(i)=nht/nthrd
!!$          nseqsum(i)=nseqsum(i-1)+nqompcur(i)
!!$      enddo
!!$      nseqsum(i)=nht
!!$      nqompcur(i)=nseqsum(i)-nseqsum(i-1)
!!$      do i=1, nthrd
!!$          do j=1, nqompcur(i)
!!$              nseqcur(j,i)=nseqsum(i-1)+j
!!$          enddo
!!$      enddo
!!$  endif
!
!    ncur=nht
!
!!!$  do i=1, nht
!!!$      nseqcur(i,1)=i
!!!$  enddo
!end subroutine


subroutine Move2SurfRect(pa, cent, slen, surf, axis)
    use PARAM
    use MCDefine
    implicit none
    real(8), intent(inout) :: pa(3)
    real(8), intent(in) :: cent(3), slen(3)
    integer, intent(in) :: surf, axis

    if (axis .eq. RADIAL) then
        select  case(surf)
            case(SOUTH)
                pa(2)=cent(2)-slen(2)/2
            case(WEST)
                pa(1)=cent(1)-slen(1)/2
            case(NORTH)
                pa(2)=cent(2)+slen(2)/2
            case(EAST)
                pa(1)=cent(1)+slen(1)/2
        end select
    else ! axial
        select  case(surf)
            case(TOP)
                pa(3)=cent(3)+slen(3)/2
            case(BOTTOM)
                pa(3)=cent(3)-slen(3)/2
        end select
    endif
end subroutine

subroutine CollectTallies()
    implicit none
    integer :: i,j,k
    do k=1,nz;do j=1,ny; do i=1,nx
!$      if (lomp) then
!$          FSD(i,j,k)=sum(FSDomp(i,j,k,:))
!$          ppow(i,j,k)=sum(ppowomp(i,j,k,:))
!$       endif
    enddo;enddo;enddo
end subroutine

function ShannonEntropy()
    implicit none
    integer :: i,j,k
    real(8) :: ShannonEntropy, sumFS, rsumFS, S_i

    sumFS=0
    do k=1,nz;do j=1,ny; do i=1,nx
        sumFS=sumFS+FSD(i,j,k)
    enddo;enddo;enddo
    rsumFS=1._8/sumFS
    ShannonEntropy=0
    do k=1,nz;do j=1,ny; do i=1,nx
        S_i=FSD(i,j,k)*rsumFS
        if (S_i>0) ShannonEntropy=ShannonEntropy-S_i*log(S_i)
    enddo;enddo;enddo
    ShannonEntropy=ShannonEntropy/log(2.0_8)
end function


subroutine ResetQueue()
    implicit none
    integer :: i
    nqidx=LTNQ(nqidx+1)
    icur=LTNQ(nqidx)
    inxt=LTNQ(nqidx+1)
    nqcur=>nqueue(:,:,icur)
    nqnxt=>nqueue(:,:,inxt)

    nqsum(1)=nqompnxt(1)
    do i=2, nthrd
        nqsum(i)=nqsum(i-1)+nqompnxt(i)
    enddo
    nqompcur=nqompnxt
    nqompnxt=0
    ncur=nqsum(nthrd)
    nnxt=0
end subroutine

subroutine ResetKeffTally()
    implicit none
    trk_gk=0
!$  trk_gkomp=0
end subroutine

!subroutine NormalizeQueue()
!    implicit none
!    real(8) :: wtsum, rnorm
!    integer :: i
!
!    wtsum=0.
!    do i=1, ncur
!        wtsum=wtsum+nqcur(i)%wt
!    enddo
!
!    rnorm=nht/wtsum
!    do i=1, ncur
!        nqcur(i)%wt=nqcur(i)%wt*rnorm
!    enddo
!end subroutine

subroutine track_keff(nst, xs)
    implicit none
    type(nstat), intent(inout) :: nst
    type(XsecSet), intent(in) :: xs
    real(8) :: lk
    lk=xs%xsnf(nst%g)/xs%xst(nst%g)*nst%wt
    !lk=xs%xsnf(nst%g)/xs%xstr(nst%g)*nst%wt
!$  if (lomp) then
!$      trk_gkomp(nst%tid)=trk_gkomp(nst%tid)+lk
!$  else
        trk_gk=trk_gk+lk
!$  endif
!!$omp critical
!    trk_gk=trk_gk+lk
!!$omp end critical
end subroutine

subroutine move(trk, nst)
    real(8), intent(in) :: trk
    type(nstat), intent(inout) :: nst
    real(8) :: inc(3)
    inc=trk*nst%dir
    nst%pr=nst%pr+inc
    nst%pa=nst%pa+inc
end subroutine

subroutine trackFSD(lmn, g, xs, trk, tid)
    integer, intent(in) :: lmn(3), g
    type(XsecSet), pointer :: xs
    real(8), intent(in) :: trk
    integer, intent(in), optional :: tid
    real(8) :: lf, lp
    lf=xs%xsnf(g)*trk
    lp=xs%xskf(g)*trk
!$  if (lomp) then
!$      FSDOMP(lmn(1), lmn(2), lmn(3), tid)=FSDOMP(lmn(1), lmn(2), lmn(3), tid)+lf
!$      ppowomp(lmn(1), lmn(2), lmn(3), tid)=ppowomp(lmn(1), lmn(2), lmn(3), tid)+lp
!$  else
        FSD(lmn(1), lmn(2), lmn(3))=FSD(lmn(1), lmn(2), lmn(3))+lf
        ppow(lmn(1), lmn(2), lmn(3))=ppow(lmn(1), lmn(2), lmn(3))+lp
!$  endif
end subroutine

subroutine tally(nst, xs, trk, tallyset)
    type(nstat), intent(inout) :: nst
    type(XsecSet), intent(in) :: xs
    real(8), intent(in) :: trk
    type(MCTallySet) :: tallyset
    integer :: ir, g
    if (.not. ltally) return
    ! only flux is accumulated so far
    !g=nst%g
    !tallyset.trkl(g,nst%ir)=tallyset.trkl(g,nst%ir)+trk*nst%wt
    !tallyset.pinpow=tallyset.pinpow+xs%xsnf(g)*trk*nst%wt
end subroutine

function findring(nr, ringsq, x, y)
    integer :: findring
    integer :: nr
    real(8), pointer, dimension(:) :: ringsq
    real(8) :: x, y, r
    integer :: i
    r = x**2+y**2
    do i=nr,1, -1
        if (r<ringsq(i)) exit
    enddo
    findring=i
end function

function collision(nst, xs) result(term)
    use RNG
    implicit none
    type(nstat), intent(inout) :: nst
    type(XsecSet), intent(in) :: xs
    logical :: term

    if (lic) then

        nst%wt=nst%wt*(xs%xst(nst%g)-xs%xsa(nst%g))/xs%xst(nst%g) !implicit capture
!       MJ Fixing Negative Scattering
!        nst%wt=nst%wt*(xs%xst(nst%g)-xs%xsa(nst%g))/(xs%xst(nst%g)+xs%xsaa(nst%g))
        nst%g=SmpSctGroup(nst,xs)
        if (nst%wt<0.25_8) then
!            if (GetRN()<nst%wt*2) then
            if (GetRNS(nst%seed)<nst%wt*2) then
                nst%wt=0.5
                term=.false.
            else ! Absorbed
                nst%wt=0.
                term=.true.
                return
            endif
        else
            term=.false.
        endif
    else
        ! Required to be implemented
        term=.true.
    endif
end function

function SmpSctGroup(nst, xs) result (gn)
    use RNG
    implicit none
    type(nstat), intent(in) :: nst
    type(XsecSet), intent(in) :: xs
    integer :: gn
    real(8), pointer, dimension(:) :: cdf
    real(8) :: xi

!    xi=GetRN()
    xi=GetRNS(nst%seed)
    cdf=>xs%SctCDF(:,nst%g)     ! BYSMC EDIT CDF scattering
    do gn=1, ng-1
        if (xi<cdf(gn)) exit
    enddo
end function

function SmpNeutron(nst, xs) result(nn)
    use RNG
    implicit none
    type(nstat), intent(in) :: nst
    type(XsecSet), pointer, intent(in) :: xs
    real(8) :: xtot, xi, xnuf, lk
    integer :: g, nn
    xtot=xs%xst(nst%g)
!    xtot=xs%xstr(nst%g)
    xnuf=xs%xsnf(nst%g)
    lk=xnuf/xtot

!    nn=1._8/MCkeff*lk*nst%wt+GetRN()
    nn=1._8/MCkeff*lk*nst%wt+GetRNS(nst%seed)
end function

function SmpDTC(nst, xs)
    use RNG
    real(8) :: SmpDTC
    type(nstat), intent(in) :: nst
    type(XsecSet), pointer, intent(in) :: xs
    SmpDTC= -1._8/xs%xst(nst%g)*log(getRNS(nst%seed))
!    SmpDTC= -1._8/xs%xstr(nst%g)*log(getRNS(nst%seed))
end function

function SmpGroup(nst, chi)
    use RNG
    implicit none
    integer :: SmpGroup
    type(nstat), intent(in) :: nst
    real(8), pointer, dimension(:) :: chi
    real(8) :: sum, rn
    integer :: g
!    rn=GetRN()
    rn=GetRNS(nst%seed)
    sum=0.
    !do g=1, ng
    do g=1, nchi
        sum=sum+chi(g)
        if (rn<sum) then
            exit
        endif
    enddo
    SmpGroup=g
end function
function SmpGroupIso(nst, xs) ! microscopic fission neutron selection
    use RNG
    USE MCXsLib, ONLY : MicChiCDF, nchi
    USE XSLIB_MOD, ONLY : MapNucl
    implicit none
    integer :: SmpGroupIso
    type(nstat), intent(in) :: nst
    type(XsecSet), pointer, intent(in) :: xs
    real(8) :: rn
    !REAL(8), POINTER :: chi(:)
    integer :: g, iso

    iso=nst%iso
    rn=GetRNS(nst%seed)
    do g=1, nchi-1
        if (rn<MicChiCDF(MapNucl(xs%idiso(iso)),g)) exit
    enddo
    SmpGroupiso=g
end function
subroutine SmpDirIso(nst, dir)
    use RNG
    use param
    implicit none
    type(nstat), intent(in) :: nst
    real(8), intent(out) :: dir(3)
    real(8) :: mu, ms, alpha
!    mu=1-2*GetRN()
!    alpha=2*pi*GetRN()
    mu=1-2*GetRNS(nst%seed)
    alpha=2*pi*GetRNS(nst%seed)
    ms=sqrt(1._8-mu**2)
    dir(1)=ms*cos(alpha)
    dir(2)=ms*sin(alpha)
    dir(3)=mu
end subroutine

!SUBROUTINE SmpDirAnIso(nst, dir)
!END SUBROUTINE

subroutine BankNeutronIso(nst, xs, nadd)
    use MCDefine
    implicit none
    type(nstat), intent(in) :: nst
    integer, intent(in) :: nadd
    type(XsecSet), pointer, intent(in) :: xs
    integer :: i, iso
    !return
    IF( nadd .EQ. 0 ) RETURN
    iso=SmpFisIso(nst, xs%FisCDF, xs%niso)  ! sample iso
    do i=1, nadd
        nqompnxt(nst%tid)=nqompnxt(nst%tid)+1
        nqnxt(nqompnxt(nst%tid),nst%tid)%lmn=nst%lmn
        nqnxt(nqompnxt(nst%tid),nst%tid)%ir=nst%ir
        nqnxt(nqompnxt(nst%tid),nst%tid)%pr=nst%pr
        nqnxt(nqompnxt(nst%tid),nst%tid)%wt=1.
        nqnxt(nqompnxt(nst%tid),nst%tid)%iso=iso
        IF( iso .EQ. 0 )THEN
            write(*,*) 'iso = 0'
        ENDIF

    enddo
end subroutine

function SmpFisIso(nst, FisCdf, niso)
    use RNG
    implicit none
    integer :: SmpFisIso
    type(nstat), intent(in) :: nst
    real(8), pointer, dimension(:,:) :: FisCdf
    INTEGER :: niso
    real(8) :: rn
    integer :: iso
    rn=GetRNS(nst%seed)
    IF( niso .EQ. 1 )THEN
        SmpFisIso=1
        RETURN
    ENDIF
    do iso=1, niso-1
        if (rn<FisCdf(nst%g, iso)) exit
    enddo
    SmpFisIso=iso
end function
subroutine BankNeutron(nst, nadd)
    use MCDefine
    implicit none
    type(nstat), intent(in) :: nst
    integer, intent(in) :: nadd
    integer :: i
    !return
    do i=1, nadd
        nqompnxt(nst%tid)=nqompnxt(nst%tid)+1
        nqnxt(nqompnxt(nst%tid),nst%tid)%lmn=nst%lmn
        nqnxt(nqompnxt(nst%tid),nst%tid)%ir=nst%ir
        nqnxt(nqompnxt(nst%tid),nst%tid)%pr=nst%pr
        nqnxt(nqompnxt(nst%tid),nst%tid)%wt=1.
    enddo
end subroutine

function dtsRect(nst, lxyz) result(ret)
    use PARAM
    use MCDefine
    implicit none
    type(nstat), intent(inout) :: nst
    real(8), intent(in) :: lxyz(3)
    real(8) :: ret
    integer :: id, isg(3), iout
    integer, parameter :: MCPLUS=1, MCMINUS=2
    integer, parameter :: AXISX=1, AXISY=2, AXISZ=3
    real(8) :: dts(3)
    real(8), dimension(3) :: dir, pr
    integer, parameter :: LTSURF(2,3) = RESHAPE((/EAST, WEST, NORTH, SOUTH, TOP, BOTTOM/), SHAPE(LTSURF))

    dir=nst%dir
    pr=nst%pr

    dts=1e10_8
    ! compare x-direction
    do id=1, 3
        if (mobility(id)) then
            if (dir(id)>0.) then
                isg(id)=MCPLUS
                dts(id)=(lxyz(id)/2-pr(id))/dir(id)
            else
                isg(id)=MCMINUS
                dts(id)=(-lxyz(id)/2-pr(id))/dir(id)
            endif
        endif
    enddo

    if (dts(1)>dts(2)) then
        if (dts(2)>dts(3)) then
            iout = AXISZ
            nst%axis=AXIAL
        else
            iout = AXISY
            nst%axis=RADIAL
        endif
    else
        if (dts(1)>dts(3)) then
            iout = AXISZ
            nst%axis=AXIAL
        else
            iout = AXISX
            nst%axis=RADIAL
        endif
    endif

    ret = dts(iout)
    nst%surf=LTSURF(isg(iout), iout)
    if( ret<0) then
        print *, "ERROR"
    endif
end function

function dtsRectCentX(nst, lxyz, lcentx, lcenty) result(ret)
    use PARAM
    use MCDefine
    implicit none
    type(nstat), intent(inout) :: nst
    real(8), intent(in) :: lxyz(3)
    logical, intent(in) :: lcentx, lcenty
    real(8) :: lbndry(2)
    real(8) :: ret
    integer :: id, isg(3), iout
    integer, parameter :: MCPLUS=1, MCMINUS=2
    integer, parameter :: AXISX=1, AXISY=2, AXISZ=3
    real(8) :: dts(3)
    real(8), dimension(3) :: dir, pr
    integer, parameter :: LTSURF(2,3) = RESHAPE((/EAST, WEST, NORTH, SOUTH, TOP, BOTTOM/), SHAPE(LTSURF))

    dir=nst%dir
    pr=nst%pr

    dts=1e10_8

    id=1
    if (mobility(id)) then
        if (lcentx) then
            lbndry(1)=lxyz(id)
            lbndry(2)=0
        else
            lbndry(1)=lxyz(id)/2
            lbndry(2)=-lxyz(id)/2
        endif
        if (dir(id)>0.) then
            isg(id)=MCPLUS
            dts(id)=(lxyz(id)/2-pr(id))/dir(id)
        else
            isg(id)=MCMINUS
            dts(id)=(-lxyz(id)/2-pr(id))/dir(id)
        endif
    endif
    id=2
    if (mobility(id)) then
        if (lcenty) then
            lbndry(1)=0
            lbndry(2)=-lxyz(id)
        else
            lbndry(1)=lxyz(id)/2
            lbndry(2)=-lxyz(id)/2
        endif
        if (dir(id)>0.) then
            isg(id)=MCPLUS
            dts(id)=(lbndry(1)-pr(id))/dir(id)
        else
            isg(id)=MCMINUS
            dts(id)=(lbndry(2)-pr(id))/dir(id)
        endif
    endif
    id=3
    if (mobility(id)) then
        if (dir(id)>0.) then
            isg(id)=MCPLUS
            dts(id)=(lxyz(id)/2-pr(id))/dir(id)
        else
            isg(id)=MCMINUS
            dts(id)=(-lxyz(id)/2-pr(id))/dir(id)
        endif
    endif

    if (dts(1)>dts(2)) then
        if (dts(2)>dts(3)) then
            iout = AXISZ
            nst%axis=AXIAL
        else
            iout = AXISY
            nst%axis=RADIAL
        endif
    else
        if (dts(1)>dts(3)) then
            iout = AXISZ
            nst%axis=AXIAL
        else
            iout = AXISX
            nst%axis=RADIAL
        endif
    endif

    ret = dts(iout)
    nst%surf=LTSURF(isg(iout), iout)
    if( ret<0) then
        print *, "ERROR"
        pause
    endif
end function

function dtsRectCentY(nst, lxyz) result(ret)
    use PARAM
    use MCDefine
    implicit none
    type(nstat), intent(inout) :: nst
    real(8), intent(in) :: lxyz(3)
    real(8) :: ret
    integer :: id, isg(3), iout
    integer, parameter :: MCPLUS=1, MCMINUS=2
    integer, parameter :: AXISX=1, AXISY=2, AXISZ=3
    real(8) :: dts(3)
    real(8), dimension(3) :: dir, pr
    integer, parameter :: LTSURF(2,3) = RESHAPE((/EAST, WEST, NORTH, SOUTH, TOP, BOTTOM/), SHAPE(LTSURF))

    dir=nst%dir
    pr=nst%pr

    dts=1e10_8
    ! compare x-direction
    do id=1, 3
        if (mobility(id)) then
            if (dir(id)>0.) then
                isg(id)=MCPLUS
                dts(id)=(lxyz(id)/2-pr(id))/dir(id)
            else
                isg(id)=MCMINUS
                dts(id)=(-lxyz(id)/2-pr(id))/dir(id)
            endif
        endif
    enddo

    if (dts(1)>dts(2)) then
        if (dts(2)>dts(3)) then
            iout = AXISZ
            nst%axis=AXIAL
        else
            iout = AXISY
            nst%axis=RADIAL
        endif
    else
        if (dts(1)>dts(3)) then
            iout = AXISZ
            nst%axis=AXIAL
        else
            iout = AXISX
            nst%axis=RADIAL
        endif
    endif

    ret = dts(iout)
    nst%surf=LTSURF(isg(iout), iout)
    if( ret<0) then
        print *, "ERROR"
        pause
    endif
end function

function dtsRectCentXY(nst, lxyz) result(ret)
    use PARAM
    use MCDefine
    implicit none
    type(nstat), intent(inout) :: nst
    real(8), intent(in) :: lxyz(3)
    real(8) :: ret
    integer :: id, isg(3), iout
    integer, parameter :: MCPLUS=1, MCMINUS=2
    integer, parameter :: AXISX=1, AXISY=2, AXISZ=3
    real(8) :: dts(3)
    real(8), dimension(3) :: dir, pr
    integer, parameter :: LTSURF(2,3) = RESHAPE((/EAST, WEST, NORTH, SOUTH, TOP, BOTTOM/), SHAPE(LTSURF))

    dir=nst%dir
    pr=nst%pr

    dts=1e10_8
    ! compare x-direction
    do id=1, 3
        if (mobility(id)) then
            if (dir(id)>0.) then
                isg(id)=MCPLUS
                dts(id)=(lxyz(id)/2-pr(id))/dir(id)
            else
                isg(id)=MCMINUS
                dts(id)=(-lxyz(id)/2-pr(id))/dir(id)
            endif
        endif
    enddo

    if (dts(1)>dts(2)) then
        if (dts(2)>dts(3)) then
            iout = AXISZ
            nst%axis=AXIAL
        else
            iout = AXISY
            nst%axis=RADIAL
        endif
    else
        if (dts(1)>dts(3)) then
            iout = AXISZ
            nst%axis=AXIAL
        else
            iout = AXISX
            nst%axis=RADIAL
        endif
    endif

    ret = dts(iout)
    nst%surf=LTSURF(isg(iout), iout)
    if( ret<0) then
        print *, "ERROR"
        pause
    endif
end function
function dtcl(pr, dir, r, laxial, rd, surf, leakaxis)
    use MCDefine
    real(8), intent(in) :: pr(3), dir(3), r, laxial
    integer, intent(in) :: rd ! in or out direction
    integer, intent(out) :: leakaxis, surf
    real(8) :: rx, ry, rsq
    real(8) :: sgn(2)
    data sgn / -1.0, +1.0 /
    real(8) :: a, b, c, d, t
    real(8) :: ox, oy, oz
    real(8) :: dtcl, dtcmp
    integer :: mk
    real(8) :: lz

    dtcl = 0._8
    rx=pr(1)
    ry=pr(2)
    rsq = sqrt(rx**2+ry**2)
    ox=dir(1);oy=dir(2);oz=dir(3);
    a = ox*ox+oy*oy
    b = (ox*rx+oy*ry)
    c = rx*rx+ry*ry-r**2
    d = b*b-a*c
    if (d>0._8) then
        t=(sgn(rd)*sqrt(d)-b)/a
        if (rd .eq. OUTDIR .or. b <= 0) dtcl = sqrt((ox*t)**2+(oy*t)**2+(oz*t)**2)
    endif

    leakaxis=RADIAL
    surf=0
    if (mobility(3)) then
        rz=pr(3)
        if (oz>0) then
            lz=(laxial*0.5-rz)/oz
            if (lz<dtcl)then
                dtcl=lz
                surf=TOP
                leakaxis=AXIAL
            else
            endif
        else
            lz=(-laxial*0.5-rz)/oz
            if (lz<dtcl)then
                dtcl=lz
                surf=BOTTOM
                leakaxis=AXIAL
            else
            endif
        endif
    endif

    ! if the neutron meet the plane surface before cylinder
!    if (oz>0) then; mk=0
!    else; mk=-1; endif
!    lz = (hacz(n+mk)-z)/mu
!    if (lz < dtcl_) then
!        dtcl_=lz
!        if (mu>0) then; dcyl=CYLT;
!        else; dcyl=CYLB; endif
!    endif
    if( dtcl .LE. 1E-8 )THEN
        dtcl=0._8
    ENDIF

endfunction


subroutine ReflCent(cutxy, nst)
    use MCDefine
    logical, intent(in) :: cutxy(2)
    type(nstat), intent(inout) :: nst
    real(8) :: adjpa

    ! for left boundary cells
    if (cutxy(1) .and. nst%pr(1)<0.) then
        nst%pr(1)=-nst%pr(1)
        adjpa=nst%pr(1)*2
        nst%pa(1)=nst%pa(1)+adjpa
        nst%dir(1)=-nst%dir(1)
        if(nst%axis .eq. 1 .and. nst%surf .eq. WEST) then
            nst%surf=EAST
        endif
    endif

    ! for top boundary cells
    if (cutxy(2) .and. nst%pr(2)>0.) then
        nst%pr(2)=-nst%pr(2)
        adjpa=nst%pr(2)*2
        nst%pa(2)=nst%pa(2)+adjpa
        nst%dir(2)=-nst%dir(2)
        if(nst%axis .eq. 1 .and. nst%surf  .eq. NORTH) then
            nst%surf=SOUTH
        endif
    endif
end subroutine



end module
