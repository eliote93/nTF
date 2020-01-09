module CMFDTally
    use MCDefine

    private
    public :: GetNeighCMFD, GetDimCMFD, GetNGCMFD, SetKeff4NBal
    public :: InitCMFDTally, AdrCMFDCell, CMFDRRTally, CurrentTallyOut, CurrentTallyIn, SrcDistTally
    public :: ResetCMFDTally, CollectTally, GenHomXsec, SetZeroCMFDTally
    public :: TCMFDCELL

    integer, private :: dimcmfd(3), ngcmfd, ngMC, nthrd
    integer, private :: nx, ny, nz
    equivalence (nx, dimcmfd(1))
    equivalence (ny, dimcmfd(2))
    equivalence (nz, dimcmfd(3))

    type TCMFDCELL
        integer :: idx=0, ixy, lmn(3)
        real(8) :: vol, rvol, area(6), rarea(6), h(6), rh(6)
        type(XsecSet) :: xs, xsb
        real(8) :: phi(2), j(6,2)
        integer :: nsrc
    contains
!        procedure:: InitializeNull
        procedure:: Initialize
        procedure:: SetVacuum
        procedure:: DumpXS
        procedure:: RestoreXS
    end type

    type(TCMFDCELL), private, pointer, dimension(:,:,:) :: CellCMFD

    ! To be accumulated
    real(8), private, pointer, dimension(:,:,:,:) :: mcphi
    real(8), private, pointer, dimension(:,:,:,:,:) :: mcj
    real(8), private, pointer, dimension(:,:,:,:) :: rtot, rnuf, rabsr
    real(8), private, pointer, dimension(:,:,:,:,:) :: rsct
    integer, private, pointer, dimension(:,:,:) :: nsrc
    ! For Cycle-wise tally
    real(8), private, pointer, dimension(:,:,:,:,:) :: mcphiomp
    real(8), private, pointer, dimension(:,:,:,:,:,:) :: mcjomp
    real(8), private, pointer, dimension(:,:,:,:,:) :: rtotomp, rnufomp, rabsromp
    real(8), private, pointer, dimension(:,:,:,:,:,:) :: rsctomp
    integer, private, pointer, dimension(:,:,:,:) :: nsrcomp
    integer, private, pointer, dimension(:) :: g2gc
    integer, private :: LTSURF(4,2) = RESHAPE((/1,2,3,4,5,6,0,0/), SHAPE(LTSURF))

    real(8), private, pointer, dimension(:,:,:,:) :: nbal
    real(8), private :: keff

contains

!subroutine InitializeNull(this)
!    use PARAM
!    implicit none
!    class (TCMFDCELL) :: this
!    this%idx=0
!end  subroutine

subroutine Initialize(this, asypitch, hz, lmn, idx, ixy, lcutxy)
    use PARAM
    implicit none
    class (TCMFDCELL) :: this
    real(8), intent(in) :: asypitch, hz
    integer, intent(in) :: lmn(3), idx, ixy
    logical, intent(in) :: lcutxy(2)

    this%idx=idx
    this%ixy=ixy
    this%lmn=lmn

    this%h(SOUTH)=asypitch
    this%h(WEST)=asypitch
    this%h(NORTH)=asypitch
    this%h(EAST)=asypitch
    if (lcutxy(1)) then
        this%h(WEST)=this%h(WEST)/2
        this%h(EAST)=this%h(EAST)/2
    endif
    if (lcutxy(2)) then
        this%h(NORTH)=this%h(NORTH)/2
        this%h(SOUTH)=this%h(SOUTH)/2
    endif
    this%h(BOTTOM+4)=hz
    this%h(TOP+4)=hz

    this%area(SOUTH)=this%h(EAST)*hz
    this%area(WEST)=this%h(NORTH)*hz
    this%area(NORTH)=this%h(EAST)*hz
    this%area(EAST)=this%h(NORTH)*hz
    this%area(BOTTOM+4)=this%h(EAST)*this%h(NORTH)
    this%area(TOP+4)=this%area(BOTTOM+4)
    this%vol=this%h(EAST)*this%h(NORTH)*hz
    this%rvol=1./this%vol
    this%rarea=1./this%area
    this%rh=1./this%h
    allocate(this%xs%xst(ngcmfd), this%xs%xsrmv(ngcmfd), this%xs%xsnf(ngcmfd), this%xs%xsa(ngcmfd), this%xs%chi(ngcmfd))
    allocate(this%xs%xssm(ngcmfd, ngcmfd))
    allocate(this%xsb%xst(ngcmfd), this%xsb%xsrmv(ngcmfd), this%xsb%xsnf(ngcmfd), this%xsb%xsa(ngcmfd), this%xsb%chi(ngcmfd))
    allocate(this%xsb%xssm(ngcmfd, ngcmfd))
end subroutine

subroutine DumpXS(this)
    use PARAM
    implicit none
    class (TCMFDCELL) :: this
    this%xsb%xstr=this%xs%xstr
    this%xsb%xsrmv=this%xs%xsrmv
    this%xsb%xsnf=this%xs%xsnf
    this%xsb%xsa=this%xs%xsa
    this%xsb%chi=this%xs%chi
    this%xsb%xssm=this%xs%xssm
end subroutine

subroutine RestoreXS(this)
    use PARAM
    implicit none
    class (TCMFDCELL) :: this
    this%xs%xstr=this%xsb%xstr
    this%xs%xsrmv=this%xsb%xsrmv
    this%xs%xsnf=this%xsb%xsnf
    this%xs%xsa=this%xsb%xsa
    this%xs%chi=this%xsb%chi
    this%xs%xssm=this%xsb%xssm
end subroutine

subroutine SetVacuum(this)
    implicit none
    class (TCMFDCELL) :: this
    this%idx=0
    this%vol=0
    this%area=0
    this%h=0
    this%rvol=0
    this%rarea=0
    this%rh=0
end subroutine

function GetNeighCMFD(lmn, dir) result(olmn)
    integer, intent(in) :: lmn(3), dir
    integer :: olmn(3)
    integer, dimension(3,6) :: LTNEIGH
    data LTNEIGH /0,-1,0, -1,0,0, 0,+1,0, +1,0,0, 0,0,-1, 0,0,+1/
    olmn=lmn+LTNEIGH(:,dir)
end function

function GetDimCMFD() result(cmfddim)
    integer :: cmfddim(3)
    cmfddim=dimcmfd
end function

function GetNGCMFD() result(cmfdng)
    integer :: cmfdng
    ngcmfd=cmfdng
end function

subroutine SetKeff4NBal(keff_)
    implicit none
    real(8), intent(in) :: keff_
    keff=keff_
end subroutine

subroutine InitCMFDTally(dimcmfd_, ngcmfd_, nthrd_, ngMC_, asypitch, hz, cmfdmap, cutxymap)
    use Xsec4MC
    implicit none
    type(TCMFDCELl), pointer :: mycell
    integer, intent(in) :: dimcmfd_(3), ngcmfd_, nthrd_, ngMC_
    real(8), intent(in) :: asypitch
    real(8), intent(in), pointer, dimension(:) :: hz
    integer, dimension(:,:), intent(in) :: cmfdmap
    logical, dimension(:,:,:), intent(in) :: cutxymap
    integer :: i, j, k, idx, ixy
    ! ngcmfd shoul be 2
    dimcmfd=dimcmfd_
    nx=dimcmfd(1)
    ny=dimcmfd(2)
    nz=dimcmfd(3)
    ngcmfd=ngcmfd_
    nthrd=nthrd_
    allocate(mcphi(ngcmfd, nx, ny, nz))
    allocate(mcj(6, ngcmfd, nx, ny, nz))
    allocate(rtot(ngcmfd, nx, ny, nz))
    allocate(rnuf(ngcmfd, nx, ny, nz))
    allocate(rabsr(ngcmfd, nx, ny, nz))
    allocate(rsct(ngcmfd, ngcmfd, nx, ny, nz))
    allocate(nsrc(nx,ny,nz))
    allocate(nbal(ngcmfd, nx, ny, nz))

    allocate(mcphiomp(ngcmfd, nx, ny, nz, nthrd))
    allocate(mcjomp(6, ngcmfd, nx, ny, nz, nthrd))
    allocate(rtotomp(ngcmfd, nx, ny, nz, nthrd))
    allocate(rnufomp(ngcmfd, nx, ny, nz, nthrd))
    allocate(rabsromp(ngcmfd, nx, ny, nz, nthrd))
    allocate(rsctomp(ngcmfd, ngcmfd, nx, ny, nz, nthrd))
    allocate(nsrcomp(nx,ny,nz,nthrd))
    nsrcomp=0

    ngMC=ngMC_
    !allocate(g2gc(ngMC))
    !call SetG2GC(g2gc)

    ! initialize mycell structure
    allocate(CellCMFD(0:nx+1,0:ny+1,0:nz+1))
    idx=0
    do k=1, nz
        ixy=0
        do j=1, ny
            do i=1, nx
                idx=idx+1
                ixy=ixy+1
                call CellCMFD(i,j,k)%Initialize(asypitch, hz(k), (/i,j,k/), idx, ixy, cutxymap(:,i,j))
                if (cmfdmap(i,j) .eq. 0 ) call CellCMFD(i,j,k)%SetVacuum()
            enddo
        enddo
    enddo
    idx=0
    ixy=0


    call ResetCMFDTally()
    call SetZeroCMFDTally()

end subroutine

subroutine ResetCMFDTally()
    mcphi=0
    mcj=0
    rtot=0
    rnuf=0
    rabsr=0
    rsct=0
    nbal=0
    ! dump cross sections

    do k=1, nz
        do j=1, ny
            do i=1, nx
#ifndef __GFORTRAN__
                if (CellCMFD(i,j,k)%idx .ne. 0) call CellCMFD(i,j,k)%DumpXS()
#else
                stop 'Unresolved failure in GNU Fortran'
#endif
            enddo
        enddo
    enddo

end subroutine


subroutine SetZeroCMFDTally()
    mcphiomp=0
    mcjomp=0
    rtotomp=0
    rnufomp=0
    rabsromp=0
    rsctomp=0
end subroutine

subroutine CollectTally()
    implicit none
    integer :: i, j, k, g, gp, dir
    real(8) :: lhs, rhs
    real(8) :: chi(2)
    chi(1)=1.
    chi(2)=0.
    do k=1, nz
        do j=1, ny
            do i=1, nx
                do g=1, ngcmfd
                    mcphi(g,i,j,k)=mcphi(g,i,j,k)+sum(mcphiomp(g,i,j,k,:))
                    rtot(g,i,j,k) =rtot(g,i,j,k) +sum(rtotomp(g,i,j,k,:))
                    rnuf(g,i,j,k) =rnuf(g,i,j,k) +sum(rnufomp(g,i,j,k,:))
                    rabsr(g,i,j,k)=rabsr(g,i,j,k)+sum(rabsromp(g,i,j,k,:))
                    do gp=1, ngcmfd
                        rsct(gp,g,i,j,k)=rsct(gp,g,i,j,k)+sum(rsctomp(gp,g,i,j,k,:))
                    enddo
                    do dir=1, 6
                        mcj(dir,g,i,j,k)=mcj(dir,g,i,j,k)+sum(mcjomp(dir,g,i,j,k,:))
                    enddo
                enddo
                do g=1, ngcmfd
                    ! balance check
                    if (rtot(g,i,j,k)>0) then
                        lhs=sum(mcj(:,g,i,j,k))+rtot(g,i,j,k)
                        rhs=chi(g)/keff*sum(rnuf(:,i,j,k))
                        do gp=1, ngcmfd
                            rhs=rhs+rsct(gp,g,i,j,k)
                        enddo
                        nbal(g,i,j,k)=(rhs-lhs)/rhs
                    endif
                enddo
                nsrc(i,j,k)=sum(nsrcomp(i,j,k,:))
                nsrcomp(i,j,k,:)=0
            enddo
        enddo
    enddo
end subroutine

subroutine GenHomXsec()
    implicit none
    type (TCMFDCELL), pointer :: pcmfd
    type (XsecSet), pointer :: xs
    integer :: i, j, k, g, gp, dir
    real(8) :: rphi

    do k=1, nz
        do j=1, ny
            do i=1, nx
                pcmfd=>CellCMFD(i,j,k)
                if (pcmfd%idx .eq. 0) cycle
                xs=>pcmfd%xs
                xs%chi(1)=1.
                xs%chi(2)=0.
                do g=1, ngcmfd
                    if (mcphi(g,i,j,k) .ne. 0.) then
                        rphi=1./mcphi(g,i,j,k)
                        xs%xst(g)=rtot(g,i,j,k)*rphi
                        xs%xsnf(g)=rnuf(g,i,j,k)*rphi
                        xs%xsa(g)=rabsr(g,i,j,k)*rphi
                        do gp=1, ngcmfd
                            xs%xssm(g,gp)=rsct(g,gp,i,j,k)*rphi
                        enddo
                        xs%xsrmv(g)=xs%xst(g)-xs%xssm(g,g)
                        pcmfd%phi(g)=mcphi(g,i,j,k)*pcmfd%rvol
                        do dir=1, 6
                            pcmfd%j(dir,g)=mcj(dir,g,i,j,k)*pcmfd%rarea(dir)
                        enddo
                    else
                        pcmfd%phi(g)=0.
                        pcmfd%j(:,g)=0.
#ifndef __GFORTRAN__
                        call pcmfd%RestoreXS()
#else
                        stop 'Unresolved failure in GNU Fortran'
#endif
                    endif
                enddo
                pcmfd%nsrc=nsrc(i,j,k)
            enddo
        enddo
    enddo
end subroutine

subroutine AdrCMFDCell(lmn, pcmfd)
    implicit none
    integer, intent(in) :: lmn(3)
    type(TCMFDCELL), pointer :: pcmfd
    pcmfd=>CellCMFD(lmn(1), lmn(2), lmn(3))
end subroutine

subroutine SrcDistTally(lmn, nsrc, tid)
    implicit none
    integer, intent(in) :: lmn(3), nsrc, tid
    integer :: i, j, k
    i=lmn(1);j=lmn(2);k=lmn(3)
    nsrcomp(i,j,k,tid)=nsrcomp(i,j,k,tid)+nsrc
end subroutine

subroutine CMFDRRTally(lmn, g, xs, trkl, tid)
    use MCDefine
    USE Xsec4MC, ONLY : g2gc
    implicit none
    integer, intent(in) :: lmn(3), g, tid
    real(8), intent(in) :: trkl
    integer :: i, j, k, gc, gp
    type(XsecSet), pointer :: xs

    gc=g2gc(g)
    i=lmn(1);j=lmn(2);k=lmn(3)

    mcphiomp(gc,i,j,k,tid)=mcphiomp(gc,i,j,k,tid)+trkl
!    rtotomp(gc,i,j,k,tid)=rtotomp(gc,i,j,k,tid)+trkl*xs%xst(g)
    rtotomp(gc,i,j,k,tid)=rtotomp(gc,i,j,k,tid)+trkl*xs%xstr(g)
    rnufomp(gc,i,j,k,tid)=rnufomp(gc,i,j,k,tid)+trkl*xs%xsnf(g)
    rabsromp(gc,i,j,k,tid)=rabsromp(gc,i,j,k,tid)+trkl*xs%xsa(g)
!    mcphiomp(gc,i,j,k,tid)=mcphiomp(gc,i,j,k,tid)+1
    if (ngCMFD .eq. 2) then
        do gp=1, ngCMFD
            rsctomp(gc,gp,i,j,k,tid)=rsctomp(gc,gp,i,j,k,tid)+trkl*xs%xssm2g(g,gp)
        enddo
    else
        do gp=1, ngMC
            rsctomp(gc,gp,i,j,k,tid)=rsctomp(gc,gp,i,j,k,tid)+trkl*xs%xssm(g,gp)
        enddo
    endif
end subroutine

subroutine CurrentTallyOut(nst, infoCMFD)
    use MCDefine
    USE Xsec4MC, ONLY : g2gc
    implicit none
    type(nstat), intent(in) :: nst
    type(CMFDTallyHelper), intent(in) :: infoCMFD
    integer :: i, j, k, is, g, tid

    if (infoCMFD%lsurf(nst%surf, nst%axis)) then
        i=infoCMFD%lmn(1);j=infoCMFD%lmn(2);k=infoCMFD%lmn(3);
        g=g2gc(nst%g)
        tid=nst%tid
        is=LTSURF(nst%surf, nst%axis)
        mcjomp(is,g,i,j,k,tid)=mcjomp(is,g,i,j,k,tid)+nst%wt
    endif
end subroutine

subroutine CurrentTallyIn(nst, infoCMFD)
    use MCDefine
    USE Xsec4MC, ONLY : g2gc
    implicit none
    type(nstat), intent(in) :: nst
    type(CMFDTallyHelper), intent(in) :: infoCMFD
    integer :: insurf
    integer :: i, j, k, is, g, tid

    if (nst%surf .eq. 0) return

    insurf=SurfOp(nst%surf, nst%axis)
    if (infoCMFD%lsurf(insurf, nst%axis)) then
        i=infoCMFD%lmn(1);j=infoCMFD%lmn(2);k=infoCMFD%lmn(3);
        g=g2gc(nst%g)
        tid=nst%tid
        is=LTSURF(insurf, nst%axis)
        mcjomp(is,g,i,j,k,tid)=mcjomp(is,g,i,j,k,tid)-nst%wt
    endif
end subroutine

end module
