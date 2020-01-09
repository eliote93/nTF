module MCCMFD
    use mtbl2g
    private
    public :: SetCMFDSet, InitMCCMFD, CMFD4MC, GetAdjFactor, CalWtAdj, CMFDCTRL

    integer, parameter :: CMFDSET_NXT=1, CMFDSET_CONTINUE=2, CMFDSET_STOP=3

    type cmfd_set
        integer, private :: period=5
        integer, private :: skip=1
        integer, private :: flag=CMFDSET_CONTINUE
    contains
        procedure :: init
    end type
    integer, parameter:: max_cmfdset = 100
    type(cmfd_set), private :: cmfdset(max_cmfdset)

    integer, private :: dimcmfd(3)
    integer, private :: ng, nx, ny, nz, nxy

    real(8), private, pointer, dimension(:,:,:) :: phi
    real(8), private, pointer, dimension(:,:) :: psi, ppsi
    real(8), private, pointer, dimension(:,:,:) :: pdf, wtadj
    real(8), private, pointer, dimension(:,:,:) :: q
    real(8), private, pointer, dimension(:,:,:,:,:) :: dtilde, dhat
    type(mb2g), private :: mt, mtorg
    real(8), private :: albedo(6)

    integer, private :: iterconv, maxin
    real(8), private :: errconv, ccin
    integer, private :: skip=0, nfdb=0, ninc=0, iset=1, istep=0
    logical, private :: lcmfdfdb(2), lmultiset
    integer, private :: nact, ninact, nskip
    integer, parameter :: INACTIVE=1, ACTIVE=2
    real(8), private :: fnorm

contains
subroutine SetCMFDSet(idx, period, skip, flag)
    integer, intent(in) :: idx, period, skip, flag
    call cmfdset(idx)%Init(period, skip, flag)
end subroutine

subroutine Init(this, period, skip, flag)
    class(cmfd_set) :: this
    integer, intent(in) :: period, skip, flag
    this%period=period
    this%skip=skip
    this%flag=flag
end subroutine

subroutine ReadCMFDInput(io)
    implicit none
    integer, intent(in) :: io
end subroutine

subroutine NormalizePhi()
    implicit none
    integer :: i, j, k, l
    real(8) :: sumFS, norm
    sumFS=sum(psi)
    fnorm=1./sumFS
    do k=1, nz
        l=0
        do j=1, ny
            do i=1, nx
                l=l+1
                phi(:,l,k)=phi(:,l,k)*fnorm
            enddo
        enddo
    enddo
end subroutine

 subroutine CMFDCTRL(lact, batch, icmfd, lfdb, lreset)
    implicit none
    logical, intent(in) :: lact
    integer, intent(in) :: batch
    integer, intent(out) :: icmfd
    logical, intent(out) :: lfdb, lreset
    integer, save :: istep=0, iskip=0, iset=1
    integer :: period
    logical :: oncmfd, onfdb

    lfdb=.false.
    lreset=.false.

    istep = istep + 1
    icmfd=istep
    if (.NOT. lact) then ! for inacitve cycles
        if (lcmfdfdb(1)) then ! fdb for inactive
            if (istep <= cmfdset(iset)%skip) then
                lfdb=.false.
            else
                lfdb=.true.
            endif

            if (istep .eq. cmfdset(iset)%period) then
                if (cmfdset(iset)%flag .eq. CMFDSET_NXT) then
                    iset = iset + 1
                    istep = 0
                    lreset = .true.
                elseif (cmfdset(iset)%flag .eq. CMFDSET_CONTINUE) then
                    iset = 0
                elseif (cmfdset(iset)%flag .eq. CMFDSET_STOP) then
                    iset = -1
                    istep = 0
                    lreset = .true.
                endif
            endif
        endif
    else
        lreset=.false.
        if (lcmfdfdb(2)) then
            lfdb=.true.
        else
            lfdb=.false.
        endif
    endif
end subroutine

subroutine InitMCCMFD(dimcmfd_, albedo_, lfdb)
    use allocs
    use bicg2g, only : mallocbicg2g
    implicit none
    integer, intent(in) :: dimcmfd_(3)
    real(8), intent(in) :: albedo_(6)
    logical, intent(in) :: lfdb(2)

    dimcmfd=dimcmfd_
    albedo(1:6)=albedo_(1:6)
    lcmfdfdb=lfdb

    nx=dimcmfd(1)
    ny=dimcmfd(2)
    nz=dimcmfd(3)
    nxy=nx*ny
    ng=ng2
    call dmalloc0(phi,1,ng,0,nxy,0,nz)
    call dmalloc(psi,nxy,nz)
    call dmalloc(ppsi,nxy,nz)
    call dmalloc(q,ng,nxy,nz)
    call dmalloc(dtilde,6,ng,nx,ny,nz)
    call dmalloc(dhat,6,ng,nx,ny,nz)
    call dmalloc(pdf,nx,ny,nz)
    call dmalloc(wtadj,nx,ny,nz)
    wtadj=1.
    call mt%newmb2g(nxy, nx, ny, nz)
    call mtorg%newmb2g(nxy, nx, ny, nz)
    ! Required to be set by input sheet
    iterconv=100
    errconv=1e-5_8
    maxin=10
    ccin=1e-2_8
    call mallocbicg2g(nxy, nx, ny, nz)
end subroutine

subroutine CalDtilde_hat()
    use CMFDTally
    implicit none
    type(TCMFDCELL), pointer :: pm, po
    integer :: i,j,k,idx,ixy,g,gp,dir
    integer :: lmn(3), olmn(3)
    real(8) :: dttemp
    idx=0
    do k=1, nz
        ixy=0
        do j=1, ny
            do i=1, nx
                idx=idx+1
                ixy=ixy+1
                lmn=(/i,j,k/)
                call AdrCMFDCell(lmn, pm)
                if (pm%idx .eq. 0) cycle
                do g=1, ng
                    do dir=1, 6
                        olmn=GetNeighCMFD(lmn,dir)
                        call AdrCMFDCell(olmn, po)
                        if (po%idx .eq. 0) then ! boundary
                            dtilde(dir,g,i,j,k)=CalDtildeBoundary()
                            dhat(dir,g,i,j,k)=CalDhatBoundary()
                            ! use typical FDM for zero-scoring node
                            if (dhat(dir,g,i,j,k) .eq. -1.) dhat(dir,g,i,j,k)=dtilde(dir,g,i,j,k)
                        else ! inner cell
                            dttemp=CalDtildeInterface()
                            dtilde(dir,g,i,j,k)=dttemp
                            dhat(dir,g,i,j,k)=CalDhatInterface(dttemp)
                        endif
                    enddo
                    if (phi(g,ixy,k) .eq. 0) phi(g,ixy,k)=pm%phi(g)
                enddo
            enddo
        enddo
    enddo
contains

function CalDtildeBoundary() result(dt)
    implicit none
    real(8) :: dt
    real(8) :: df, bt
    dt=0.
    df=1./3/pm%xs%xst(g)
    bt=df/pm%h(dir)
    dt=2*albedo(dir)*bt/(albedo(dir)+2*bt)
    ! do nothing
end function

function CalDtildeInterface() result(dt)
    implicit none
    real(8) :: dt
    real(8) :: dl, dr, bl, br, th
    dl=1./3/pm%xs%xst(g)
    bl=dl/pm%h(dir)
    dr=1./3/po%xs%xst(g)
    br=dl/po%h(dir)
    dt=2*bl*br/(bl+br)
end function

function CalDhatBoundary() result(dh)
    implicit none
    real(8) :: dh
    real(8) :: jnet, phim
    phim=pm%phi(g)
    if (phim .eq. 0.) then
        dh=-1.0
!        dh=0.
    else
        jnet=pm%j(dir, g)
        dh=jnet/phim
    endif

end function

function CalDhatInterface(dt) result(dh)
    implicit none
    real(8) :: dh
    real(8), intent(in) :: dt
    real(8) :: jnet, phim, phio

    phim=pm%phi(g)
    phio=po%phi(g)
    jnet=pm%j(dir, g)

!    if(phim .eq. 0. .or. phio .eq. 0.) then
    if (jnet .eq. 0) then
        dh=0.
    else
        dh=-(jnet+dt*(phio-phim))/(phio+phim)
    endif
end function
end subroutine

subroutine MakeMT()
    use mtbl2g
    use MCDefine
    use CMFDTally
    implicit none
    type(TCMFDCELL), pointer :: pm, po
    type(XsecSet), pointer :: xs
    integer :: i,j,k,l,idx,g,gp,dir
    integer :: lmn(3), olmn(3)
    real(8) :: vol, area, dt, dh, ineigh
    integer :: nxs, nxe

    do k=1, nz
        l=0
        do j=1, ny
            nxs=0
            nxe=1
            do i=1, nx
                l=l+1
                lmn=(/i,j,k/)
                call AdrCMFDCell(lmn, pm)
                if (pm%idx .eq. 0) cycle
                if (nxs .eq. 0) nxs=i
                nxe=i

                vol=pm%vol
                xs=>pm%xs
                call mt%setnodel(i,j,l)
                do g=1, ng
                    call mt%setmb2g(l,k,0,indm24(g,g), xs%xsrmv(g)*vol)
                    do dir=1, 6
                        area=pm%area(dir)
                        olmn=GetNeighCMFD(lmn,dir)
                        call AdrCMFDCell(olmn, po)
                        dt=dtilde(dir,g,i,j,k)
                        dh=dhat(dir,g,i,j,k)
                        if (po%idx .ne. 0) then
                            call mt%accmb2g(l,k,0,indm24(g,g),  (dt-dh)*area)
                            call mt%setmb2g(l,k,dir,g,         -(dt+dh)*area, GetiNeigh())
                        else ! boundary
                            call mt%accmb2g(l,k,0,indm24(g,g),   dh*area)
                        endif
                    enddo
                enddo
                ! set group coupling
                call mt%setmb2g(l,k,0,indm24(1,2), -xs%xssm(2,1)*vol) ! scattering from thermal to fast
                call mt%setmb2g(l,k,0,indm24(2,1), -xs%xssm(1,2)*vol) ! scattering from fast to thermal
            enddo
            call mt%setnxse(j,nxs,nxe)
        enddo
    enddo
    call mtorg%copymb2g(mt)
    call facilu2g(mt)
contains

function GetiNeigh() result(ineigh)
    implicit none
    integer :: ineigh
    if (dir>4) then
        ineigh=po%lmn(3)
    else
        ineigh=po%ixy
    endif
end function
end subroutine

subroutine CMFD4MC(kcmfd, itercmfd, errcmfd, shncmfd)
    implicit none
    real(8), intent(out) :: kcmfd, errcmfd, shncmfd
    integer :: itercmfd
    logical :: lwldt
    real(8) :: keff, err
    integer :: iter
    lwldt=.true.
    keff=1
    call CalDtilde_hat()
    call MakeMT()
    call SolveCMFD(lwldt, keff, iter, err)
!    call SolveCMFD()
    kcmfd=keff
    itercmfd=iter
    errcmfd=err
    shncmfd=ShannonEntropy()
    ! normalization
    call NormalizePhi()
end subroutine

function ShannonEntropy()
    implicit none
    real(8), parameter :: rlog2=3.321928094887_8
    integer :: k,l
    real(8) :: ShannonEntropy, sumFS, rsumFS, S_i
    sumFS=0
    do k=1, nz
        do l=1, nxy
            sumFS=sumFS+psi(l,k)
        enddo
    enddo
    rsumFS=1./sumFS
    ShannonEntropy=0.
    do k=1, nz
        do l=1, nxy
            S_i=psi(l,k)*rsumFS
            if (S_i>0) ShannonEntropy=ShannonEntropy-S_i*log(S_i)*rlog2
        enddo
    enddo
end function

subroutine MakeFS()
    use CMFDTally
    use MCDefine
    implicit none
    type(TCMFDCELL), pointer :: pm
    type(XsecSet), pointer :: xs
    real(8) :: rk, vol
    integer :: k,j,i,ixy,g,lmn(3)
    do k=1, nz
        ixy=0
        do j=1, ny
            do i=1, nx
                ixy=ixy+1
                lmn=(/i,j,k/)
                call AdrCMFDCell(lmn, pm)
                if (pm%idx .eq. 0) cycle
                xs=>pm%xs
                ppsi(ixy,k)=psi(ixy,k)
                psi(ixy,k)=pm%vol*(xs%xsnf(1)*phi(1,ixy,k) &
                                  +xs%xsnf(2)*phi(2,ixy,k))
            enddo
        enddo
    enddo
end subroutine

subroutine MakeQ(keff)
    implicit none
    real(8), intent(in) :: keff
    real(8) :: rk, vol
    integer :: k,ixy
    rk=1./keff
    do k=1, nz
        ixy=0
        do ixy=1, nxy
            q(1,ixy,k)=psi(ixy,k)*rk
            q(2,ixy,k)=0.
        enddo
    enddo
end subroutine

subroutine Shiftmt(ke)
    use CMFDTally
    use MCDefine
    implicit none
    real(8), intent(in) :: ke
    type(TCMFDCELL), pointer :: pm
    type(XsecSet), pointer :: xs
    real(8) :: rke, vol
    integer :: k,j,i,ixy,g,lmn(3)
    rke=1./ke
    call mt%copymb2g(mtorg)
    do k=1, nz
        ixy=0
        do j=1, ny
            do i=1, nx
                ixy=ixy+1
                lmn=(/i,j,k/)
                call AdrCMFDCell(lmn, pm)
                if (pm%idx .eq. 0) cycle
                xs=>pm%xs
                call mt%accmb2g(ixy,k,0,indm24(1,1), -rke*xs%xsnf(1)*pm%vol)
                call mt%accmb2g(ixy,k,0,indm24(1,2), -rke*xs%xsnf(2)*pm%vol)
            enddo
        enddo
    enddo
end subroutine

subroutine SolveCMFD(lwldt, keff, iter, err)
    use bicg2g
    use bicg2g_interface
    implicit none
    logical, intent(in) :: lwldt
    real(8), intent(inout) :: keff
    integer, intent(out) :: iter
    real(8), intent(out) :: err

    real(8) :: kr, kd, ke
    real(8) :: r20, r2
    integer :: iti
    integer :: i,j,k,g

    err=1e10
    iter=0
    if (lwldt) kr=keff
    call MakeFS()
    do while (iter .lt. iterconv)
        iter=iter+1
        call MakeQ(keff)
        call initbicg2g(mt,phi,q,r20)
        if (r20<1e-10_8 .and. iter > 2) exit
        r2=r20;iti=0;
        do while (iti .lt. maxin .and. r2/r20 .gt. ccin)
            iti=iti+1
            call solbicg2g(mt,r20,r2,phi)
        enddo
        call MakeFS()
        call UpdateKeffNErr()
        if (lwldt) then
            kd=shiftk(iter, 0.01_8)
            ke=kr+kd
            keff = 1._8/(1._8/kr-1._8/ke)
            call shiftmt(ke)
        endif
        if (err<errconv) exit
    enddo
    if (lwldt) keff=kr
contains

subroutine UpdateKeffNErr()
    implicit none
    real(8) :: upper, upper2, lower
    integer :: ixy, k
    upper=0
    lower=0
    upper2=0
    do k=1, nz
        do ixy=1, nxy
            upper = upper + psi(ixy,k)*psi(ixy,k)
            lower = lower + psi(ixy,k)*ppsi(ixy,k)
            upper2 = upper2 + (psi(ixy,k)- ppsi(ixy,k))**2
        enddo
    enddo
    if(lwldt .and. iter>=2) then
        keff=keff*upper/lower
        kr=(keff*ke)/(keff+ke)
    else
        keff=keff*upper/lower
    endif
    lower = sqrt(upper)
    upper = sqrt(upper2)
    err = upper/lower
end subroutine

function shiftk(iter, delk) result (ret)
!    real(8) :: SHIFT_GUIDE(17)
!    data SHIFT_GUIDE /1.0_8, 0.75_8, 0.50_8, 0.35_8, 0.25_8, &
!                    0.20_8, 0.15_8, 0.125_8, 0.100_8, 0.075_8, 0.050_8, &
!                    0.035_8, 0.025_8, 0.020_8, 0.015_8, 0.0125_8, 0.0100_8/
    integer, parameter :: nguide = 5
    real(8) :: SHIFT_GUIDE(nguide)
    data SHIFT_GUIDE /0.25_8, 0.125_8, 0.075_8, 0.035_8, 0.01_8/
    integer, intent(in) :: iter
    real(8), intent(in) :: delk
    real(8) :: ret
    ret = delk
    if (iter .le. nguide) then
        if (ret<SHIFT_GUIDE(iter)) ret = SHIFT_GUIDE(iter)
    endif
end function
end subroutine

subroutine GetPDF()
    implicit none
    integer :: i, j, k, l
    real(8) :: sumFS, norm
    do k=1, nz
        l=0
        do j=1, ny
            do i=1, nx
                l=l+1
                pdf(i,j,k)=psi(l,k)*fnorm
            enddo
        enddo
    enddo
end subroutine

subroutine CalWtAdj(nht)
    use CMFDTally
    use MCDefine
    implicit none
    integer, intent(in) :: nht
    type(TCMFDCELL), pointer :: pm
    integer :: k,j,i,ixy,g,lmn(3)
    real(8) :: nexp
    call GetPDF()
    do k=1, nz
        ixy=0
        do j=1, ny
            do i=1, nx
                ixy=ixy+1
                lmn=(/i,j,k/)
                call AdrCMFDCell(lmn, pm)
                if (pm%idx .eq. 0) cycle
                nexp=nht*pdf(i,j,k)
                wtadj(i,j,k)=nexp/pm%nsrc
            enddo
        enddo
    enddo
end subroutine

function GetAdjFactor(lmn) result(fwt)
    implicit none
    integer, intent(in) :: lmn(3)
    real(8) :: fwt
    fwt=wtadj(lmn(1), lmn(2), lmn(3))
!    print *, fwt
end function

end module

