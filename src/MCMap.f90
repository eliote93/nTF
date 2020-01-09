module MCMap
    use MCDefine
    USE MCXsLib, ONLY: MicXS, MCXS, Fxr2XsMap
    type(cellmap_elmt), pointer, private, dimension(:,:,:) :: map_mc
    type(STR_RectFuel), pointer, private, dimension(:) :: frectarray
    type(STR_RectGap), pointer, private, dimension(:) :: grectarray
    type(STR_RectVac), pointer, private, dimension(:) :: vrectarray
    type tlmn
        integer:: l, m, n
    end type

    integer, private :: nxm, nym, nzm
    integer, private, pointer, dimension(:,:,:) ::  idx3to1, valid
    type(tlmn), private, pointer, dimension(:) :: idx1to3
    integer, private :: nTotMod
    type(nfs), private, pointer, dimension(:) :: basefs
    integer, private :: nbasefs
    LOGICAL :: lbench, lmic
contains
function GetTotNMod()
    integer :: GetTotNMod
    GetTotNMod=nTotMod
end function
function maplmn(lmn)
    integer :: maplmn
    integer, intent(in) :: lmn(3)
    maplmn=idx3to1(lmn(1), lmn(2), lmn(3))
end function

function GetNBaseFS()
    integer :: GetNBaseFS
    GetNBaseFS=nbasefs
end function

function GetBaseFS(idx)
    type(nfs) :: GetBaseFS
    integer, intent(in) :: idx
    GetBaseFS=basefs(idx)
end function

subroutine mapidx(idx, lmn)
    integer, intent(in) :: idx
    integer, intent(out) :: lmn(3)
    lmn=(/ idx1to3(idx)%l, idx1to3(idx)%m, idx1to3(idx)%n /)
end subroutine


subroutine InitMCMap(core, asyinfo, pininfo, cellinfo, FmInfo, ctrlMC, asypitch, albedo, ledge, lBenchxs, lmicxs, ng)
    use typedef
    use BenchXs
    use geom, only : nz, hz
!    use geom,           only: asypitch, albedo
!    use CMFDTally, only: InitCMFDTally
    implicit none
!    use MCModular, only : AddInitialSource
    type(coreinfo_type), intent(in) :: core
    type(asyinfo_type), pointer, intent(in), dimension(:) :: asyinfo
    type(pininfo_type), pointer, intent(in), dimension(:) :: pininfo
    type(cell_type), pointer, intent(in), dimension(:) :: cellinfo
    type(FmInfo_type) :: Fminfo
    type(MCCTRL), intent(inout) :: ctrlMC
    !integer, intent(in) :: nz
    !real(8), pointer, dimension(:), intent(in) :: hz
    real(8), intent(in) :: asypitch, albedo(10)
    logical, intent(in) :: ledge

    integer :: ngMC, gThrm
    type(asyinfo_type), pointer :: myasy
    type(cell_type), pointer :: mycell
    type(pininfo_type), pointer :: mypin
    type(nfs), pointer :: nf
    real(8), pointer, dimension(:) :: ptest
    real(8) :: pax, pay, pcx, pcy
    real(8) :: hzacc(0:nz)

    integer :: ix, iy, icx, iax, ipx, icy, iay, ipy, icz, iaz, ipz
    integer :: icyt, nv, ishft
    integer :: idx, ir, ic
    integer :: iasy, ipin, icell, tasy, tpin, tcell, txs, imod
    integer :: icmfd
    integer :: nfrect, ifrect
    integer :: ngrect, igrect
    integer :: nvrect, ivrect
    integer :: l, m, n
    integer :: nxt, nyt
    real(8) :: slen
    logical :: lcentx, lcenty, lodd, lcutx, lcuty
    integer :: ladj, madj
    type(STR_RectFuel), pointer :: rfcell
    type(STR_RectGap), pointer :: rgcell
    type(STR_RectVac), pointer :: rvcell
    type(MCTallySet), pointer :: tally
    integer, pointer, dimension(:) :: l2d, m2d
    real(8), pointer, dimension(:) :: pax2d, pay2d
    LOGICAL :: lbenchxs, lmicxs
    INTEGER :: ixy, iz
    INTEGER :: ng
    TYPE(FxrInfo_Type), POINTER :: myFXR
    INTEGER :: FxrIdxSt, ifxr
    ! Analyze the core configuration
    allocate(ctrlMC%CoreMap(Core%nxa, Core%nya))
    allocate(ctrlMC%cutxy(2,Core%nxa, Core%nya))
    lbench=lbenchxs
    lmic=lmicxs
    iasy=0
    nxm=0;nym=0
    do icyt=1, Core%nya
        icy=Core%nya+1-icyt
        nv=0
        do icx=1, Core%nxa
            iasy=Core%coreidx(icx,icy)
            if (iasy>0) then
                ctrlMC%CoreMap(icx,icyt)=Core%coremap(iasy)
                nv=nv+1
            else
                ctrlMC%CoreMap(icx,icyt)=0
            endif
        enddo
        ishft=(Core%nxa-nv)/2
        if (albedo(WEST) .ne. 0 .and. ishft>0) then
            do icx=nv, 1, -1
                ctrlMC%CoreMap(ishft+icx,icyt)=ctrlMC%CoreMap(icx, icyt)
            enddo
            ctrlMC%CoreMap(1:ishft, icyt)=0
        endif

        nxt=0
        do icx=1, Core%nxa
            tasy=ctrlMC%CoreMap(icx,icyt)
            if (tasy .eq. 0) cycle
            myasy=>asyinfo(tasy)
            ctrlMC%cutxy(1,icx,icyt)=myasy%lcenty .or. myasy%lcentxy
            ctrlMC%cutxy(2,icx,icyt)=myasy%lcentx .or. myasy%lcentxy
            nxt=nxt+myasy%nx
            nyt=myasy%ny
        enddo
        if (nxt>nxm) nxm=nxt
        nym=nym+nyt
    enddo

    allocate(l2d(Core%nxa+1), m2d(Core%nya+1))
    allocate(pax2d(Core%nxa), pay2d(Core%nya))

    m2d=0
    pay2d=0
    nym=0
    do icy=1, Core%nya
        do icx=1, Core%nxa
            tasy=ctrlMC%CoreMap(icx,icy)
            if (tasy .ne. 0) then
                myasy=>asyinfo(tasy)
                nyt=myasy%ny
                slen=asypitch
                if (myasy%lcentxy .or. myasy%lcentx) slen=slen/2
                exit
            endif
        enddo
        m2d(icy+1)=m2d(icy)+nyt
        if (icy .ne. Core%nya) then
            pay2d(icy+1)=pay2d(icy)+slen
        endif
        nym=nym+nyt
    enddo

    l2d=0
    pax2d=0
    nxm=0
    do icx=1, Core%nxa
        do icy=1, Core%nya
            tasy=ctrlMC%CoreMap(icx,icy)
            if (tasy .ne. 0) then
                myasy=>asyinfo(tasy)
                nxt=myasy%nx
                slen=asypitch
                if (myasy%lcentxy .or. myasy%lcenty) slen=slen/2
                exit
            endif
        enddo
        l2d(icx+1)=l2d(icx)+nxt
        if (icx .ne. Core%nxa) then
            pax2d(icx+1)=pax2d(icx)+slen
        endif
        nxm=nxm+nxt
    enddo

    allocate(ctrlMC%lwrt(nxm, nym))
    allocate(ctrlMC%lgap(nxm, nym))
    ctrlMC%lwrt = .false.
    ctrlMC%lgap = .false.
    ctrlMC%lwrtx(1)=nxm
    ctrlMC%lwrtx(2)=0
    ctrlMC%lwrty(1)=nym
    ctrlMC%lwrty(2)=0

    lcentx=.false.
    lcenty=.false.
    lodd=.false.
    ladj=0
    madj=0
    if (.not. ledge) then
        if (albedo(WEST) .eq. 0.) lcenty=.true.
        if (albedo(NORTH) .eq. 0.) lcentx=.true.
        if (mod(asyinfo(1)%nx, 2) .eq. 1) lodd=.true.
    endif


    !pax2d=0
    !do icx=1, Core%nxa-1
    !    tasy=ctrlMC%CoreMap(icx, Core%nya)
    !    l2d(icx)=l2d(icx-1)+asyinfo(tasy)%nx
    !    pax2d(icx)=pax2d(icx-1)+asypitch
    !enddo
    !if (lcenty) pax2d(2:)=pax2d(2:)-asypitch/2
    !m2d=0
    !pay2d=0
    !do icy=2, Core%nya
    !    tasy=ctrlMC%CoreMap(1, icy)
    !    m2d(icy)=l2d(icy-1)+asyinfo(tasy)%ny
    !    pay2d(icy)=pay2d(icy-1)+asypitch
    !enddo
    !if (lcentx) pay2d(Core%nya)=pay2d(Core%nya)-asypitch/2

    ngMC=ngben

    !nxm=Core%nxa*asyinfo(1)%nx
    !nym=Core%nya*asyinfo(1)%ny
    nzm=nz



    nbasefs=0
    allocate(basefs(nxm*nym*nzm))
    allocate(map_mc(nxm,nym,nzm))

    idx = 0
    nfrect=0

    allocate(valid(0:nxm+1, 0:nym+1, 0:nzm+1))

    ! Set boundary condition according to albedo
    if (albedo(SOUTH) .eq. 0.5) then;valid(:,0,:)=VAC
    elseif (albedo(SOUTH) .eq. 0.0) then;valid(:,0,:)=REFL
    else;
    endif
    if (albedo(WEST) .eq. 0.5) then;valid(0,:,:)=VAC
    elseif (albedo(WEST) .eq. 0.0) then;valid(0,:,:)=REFL
    else;
    endif
    if (albedo(NORTH) .eq. 0.5) then;valid(:,nym+1,:)=VAC
    elseif (albedo(NORTH) .eq. 0.0) then;valid(:,nym+1,:)=REFL
    else;
    endif
    if (albedo(EAST) .eq. 0.5) then;valid(nxm+1,:,:)=VAC
    elseif (albedo(EAST) .eq. 0.0) then;valid(nxm+1,:,:)=REFL
    else;
    endif
    if (albedo(BOTTOM+4) .eq. 0.5) then;valid(:,:,0)=VAC
    elseif (albedo(BOTTOM+4) .eq. 0.0) then;valid(:,:,0)=REFL
    else;
    endif
    if (albedo(TOP+4) .eq. 0.5) then;valid(:,:,nzm+1)=VAC
    elseif (albedo(TOP+4) .eq. 0.0) then;valid(:,:,nzm+1)=REFL
    else;
    endif

!    valid=REFL ! Set boundary condition

    allocate(idx3to1(nxm, nym, nzm))
    allocate(idx1to3(nxm*nym*nzm))

    hzacc(0)=0.
    do icz=1, nz
        hzacc(icz)=hzacc(icz-1)+hz(icz)
    enddo


    ! Assigning global map
    nfrect=0
    ngrect=0
    nvrect=0
!    iasy=0
    imod=0
    do icy=1, Core%nya
        do icx=1, Core%nxa
            iasy=iasy+1
            tasy=ctrlMC%CoreMap(icx, icy)
            if (tasy .eq. 0) cycle
            myasy=>asyinfo(tasy)
            ipin=0
            m=m2d(icy)
            if (myasy%lfuel) then
                if (ctrlMC%lwrtx(1)>l2d(icx))   ctrlMC%lwrtx(1)=l2d(icx)
                if (ctrlMC%lwrtx(2)<l2d(icx+1)) ctrlMC%lwrtx(2)=l2d(icx+1)
                if (ctrlMC%lwrty(1)>m2d(icy))   ctrlMC%lwrty(1)=m2d(icy)
                if (ctrlMC%lwrty(2)<m2d(icy+1)) ctrlMC%lwrty(2)=m2d(icy+1)
            endif

            do iay=1, myasy%ny
                m=m+1
                l=l2d(icx)
                do iax=1, myasy%nx
                    l=l+1
                    ipin=(myasy%ny-iay)*myasy%nx+iax ! top2bot -> bot2top
                    !ipin=(iay-1)*myasy%nx+iax  ! BYS edit 16/12/13
                    tpin=myasy%pin(ipin)
                    mypin=>pininfo(tpin)
                    if (mypin%lfuel) ctrlMC%lwrt(l,m)=.true.
                    !do icell=1, pininfo(1)%ncell
                    do icell=1, mypin%ncell
                        n=icell
                        imod=imod+1
                        tcell=mypin%cell(icell)
                        mycell=>cellinfo(tcell)
                        if (mycell%lgap) ctrlMC%lgap(l,m)=.true.

                        idx1to3(imod)%l=l
                        idx1to3(imod)%m=m
                        idx1to3(imod)%n=n
                        idx3to1(l,m,n)=imod
                        if(tasy .eq. 0) then ! vacuum
                            nvrect=nvrect+1
                            map_mc(l,m,n)%telmt=VRECT
                            map_mc(l,m,n)%ielmt=nvrect
                            valid(l,m,n)=VRECT
                        elseif (mycell%geom%nbd .eq. 4 .and. mycell%geom%lcircle) then ! including non-fuel but have annular geometries
                            nfrect=nfrect+1
                            map_mc(l,m,n)%telmt=FRECT
                            map_mc(l,m,n)%ielmt=nfrect
                            ! assigning global map
                            valid(l,m,n)=FRECT
!                        elseif (mycell%lfuel .eq. .false. .and. mycell%geom%nbd .eq. 4  .and. mycell%geom%ncircle .eq. 0) then ! gap cell
                        elseif (mycell%geom%nbd .eq. 4  .and. .not. mycell%geom%lcircle) then ! gap cell
                            ngrect=ngrect+1
                            map_mc(l,m,n)%telmt=GRECT ! Duplicated information. Required to be modified.
                            map_mc(l,m,n)%ielmt=ngrect
                            valid(l,m,n)=GRECT
                        else
                            print *, "ERR] not assigned MC module:", l,m,n
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo
    ctrlMC%lwrtx(1)=ctrlMC%lwrtx(1)+1
    ctrlMC%lwrty(1)=ctrlMC%lwrty(1)+1
    allocate(frectarray(nfrect))
    allocate(grectarray(ngrect))
    allocate(vrectarray(nvrect))
    nTotMod=imod



    ! assuming the same asy config.
!    iasy=0
    ipin=0
    imod=0
    pay=0
    do icy=1, Core%nya
        pay=pay2d(icy)
        do icx=1, Core%nxa
!            iasy=iasy+1
!            tasy=Core%coremap(iasy)
            tasy=ctrlMC%CoreMap(icx, icy)
            if (tasy .eq. 0) cycle
            !if (tasy .eq. 0) then
            !    myasy=>asyinfo(1)
            !else
            !    myasy=>asyinfo(tasy)
            !endif
            myasy=>asyinfo(tasy)
            pax=pax2d(icx)
            pcy=0
            ipin=0
            do iay=1, myasy%ny
                pcx=0
                !if (lcenty .and. (icy .eq. Core%nya) ) then
                !    if (iasy > (myasy%ny+1)/2 ) then
                !        cycle ! check validity for y-dir
                !    elseif( lodd .and. (myasy%ny+1)/2 .eq. iasy) then
                !        lcuty=.true.
                !    else
                !    endif
                !endif
                do iax=1, myasy%nx
!                    ipin=ipin+1
                    ipin=(myasy%ny-iay)*myasy%nx+iax
                    !ipin=(iay-1)*myasy%nx+iax  ! BYS edit 16/12/13
                    !if (lcentx) then
                    !    if (l <= myasy%nx/2)  then
                    !        cycle ! check validity for x-dir
                    !    elseif( lodd .and. l .eq. (myasy%nx+1)/2) then
                    !        lcutx=.true.
                    !    else
                    !    endif
                    !endif
                    tpin=myasy%pin(ipin)
                    mypin=>pininfo(tpin)
                    !lcutx=mypin%lcenty .or. mypin%lcentxy
                    !lcuty=mypin%lcentx .or. mypin%lcentxy
                    do icell=1, mypin%ncell
                        imod=imod+1
                        tcell=mypin%cell(icell)
                        mycell=>cellinfo(tcell)

                        lcutx=mycell%lcenty .or. mycell%lcentxy   !---BYS edit 16/12/13 : Gapcell CentXY bug in pininfo
                        lcuty=mycell%lcentx .or. mycell%lcentxy
                        l=idx1to3(imod)%l
                        m=idx1to3(imod)%m
                        n=idx1to3(imod)%n
                        !--- FOR lXSlib
                        iz=n
                        !ixy=(m-1)*Core%nx+l
                        !ixy=(Core%ny-m)*Core%nx+l
                        !FxrIdxSt=Core%Pin(ixy)%FxrIdxSt
                        ixy=(myasy%ny-iay)*myasy%nx+iax
!                        iasy=(Core%nya-icy)*Core%nxa+icx   !--- CNJ Edit : Segmentation Fault Fix
                        iasy=Core%CoreIdx(icx,Core%nya-icy+1)
                        ipin=Core%Asy(iasy)%GlobalPinIdx(ixy)
                        FxrIdxSt=Core%Pin(ipin)%FxrIdxSt

                        if (map_mc(l,m,n)%telmt .eq. VRECT) then
                            call AssignVRect()
                        elseif (map_mc(l,m,n)%telmt .eq. FRECT) then
                            call AssignFRect()
                        elseif (map_mc(l,m,n)%telmt .eq. GRECT) then
                            if (tasy .eq. 0) then
                                print *, "TTT"
                            endif
                            call AssignGRect()
                        endif
                    enddo
                    pcx=pcx+mycell%geom%lx
                enddo
                pcy=pcy+mycell%geom%ly
            enddo
        enddo
    enddo

    ! assigning MCCTRL
    ctrlMC%nxyz=(/nxm, nym, nzm/)
    ctrlMC%nxyzCMFD=(/Core%nxa, Core%nya, nzm/)
    deallocate(l2d, m2d, pax2d, pay2d)
contains

subroutine AssignVRect()
    ! assigning vacuum cells
    ivrect=map_mc(l,m,n)%ielmt
    rvcell=>vrectarray(ivrect)
    rvcell%idx=imod
    rvcell%lmn=(/l,m,n/)
 !   allocate(rgcell%tally(1))
    rvcell%slen(1)=mycell%geom%lx
    rvcell%slen(2)=mycell%geom%ly
    rvcell%slen(3)=hz(icell)

    rvcell%cut=.false.
    rvcell%cutxy=.false.
    if (lcutx) then
        rvcell%cut=.true.
        rvcell%cutxy(1)=.true.
        rvcell%slen(1)=rvcell%slen(1)*2
        rvcell%cent(1)=pax+pcx
    else
        rvcell%cent(1)=pax+pcx+rvcell%slen(1)/2
    endif
    if (lcuty) then
        rvcell%cut=.true.
        rvcell%cutxy(2)=.true.
        rvcell%slen(2)=rvcell%slen(2)*2
        rvcell%cent(2)=pay+pcy+rvcell%slen(2)/2
    endif
    rvcell%cent(2)=pay+pcy+rvcell%slen(2)/2
    rvcell%cent(3)=hzacc(icell-1)+rvcell%slen(3)/2

end subroutine

subroutine AssignFRect()
    ifrect=map_mc(l,m,n)%ielmt
    rfcell=>frectarray(ifrect)
    rfcell%FxrIdxSt=FxrIdxSt
    rfcell%idx=imod
    rfcell%idxt=ifrect
    rfcell%lmn=(/l,m,n/)
    rfcell%nr=mycell%geom%ncircle
    allocate(rfcell%rr(rfcell%nr))
    allocate(rfcell%rrsq(rfcell%nr))
    allocate(rfcell%xs(0:rfcell%nr))
    allocate(rfcell%tally%trkl(ngMC, 0:rfcell%nr))
!    allocate(rfcell%tally(0:rfcell%nr))
    rfcell%slen(1)=mycell%geom%lx
    rfcell%slen(2)=mycell%geom%ly
    rfcell%slen(3)=hz(icell)
    rfcell%vol=rfcell%slen(1)*rfcell%slen(2)*rfcell%slen(3)
    ! Add initial source ---------------------------
    if(mycell%lfuel) then
        nbasefs=nbasefs+1
        nf=>basefs(nbasefs)
        nf%lmn=rfcell%lmn;nf%idx=rfcell%idx;nf%ir=rfcell%nr
        nf%pr=0;nf%wt=1
    endif
!   call AddInitialSource(nf)
    ! ----------------------------------------------

    txs=mycell%ireg(mycell%ndivazi)
    IF( lbench )THEN !--- BYSMC EDIT 16/11/09
        call AdrXsBen(rfcell%xs(0), txs)
    ELSE
        ifxr=FxrIdxSt
        iz=n
        IF( lmic )THEN
            rfcell%xs(0)=micXs(fxr2xsmap(ifxr,iz))
        ELSE
            rfcell%xs(0)=mcXs(ifxr,iz)
        ENDIF

        !myFxr=>FmInfo%FXR(ifxr,iz)
        !call AdrXsLib(rfcell%xs(0), myFxr, ng)
    ENDIF

!    call AllocTallySet(ngMC, rfcell%nr)
    rfcell%nor=0;
    do ir=1, rfcell%nr
        ! assining radius
        rfcell%rr(ir)=mycell%geom%circle(3,ir)
        rfcell%rrsq(ir)=rfcell%rr(ir)**2
        txs=mycell%ireg((ir+1)*mycell%ndivazi)
        !IF( rfcell%CUTXY(1) )THEN
        IF( lcutx )THEN
            IF( rfcell%rr(ir) .GT. rfcell%slen(1) )THEN
                rfcell%nor=rfcell%nor+1;  ! # of radius larger than half pitch
            ENDIF
        ELSE
            IF( 2*rfcell%rr(ir) .GT. rfcell%slen(1) )THEN
                rfcell%nor=rfcell%nor+1;  ! # of radius larger than half pitch
            ENDIF
        ENDIF
        IF( lbench )THEN !--- BYSMC EDIT 16/11/09
            call AdrXsBen(rfcell%xs(ir), txs)
        ELSE
            ifxr=FxrIdxSt+ir
            iz=n
            IF( lmic )THEN
                rfcell%xs(ir)=micXs(fxr2xsmap(ifxr,iz))
            ELSE
                rfcell%xs(ir)=mcXs(ifxr,iz)
            ENDIF

            !ifxr=FxrIdxSt+ir
            !myFxr=>FmInfo%FXR(ifxr,iz)
            !rfcell%xs(ir)%ifxr=ifxr
            !rfcell%xs(ir)%iz=iz
            !call AdrXsLib(rfcell%xs(ir), myFxr, ng)
        ENDIF

        !call AllocTallySet(rfcell%tally(ir), ngMC)
    enddo

    ! Assigning CMFD info
    rfcell%infoCMFD%lmn=(/icx, icy, icell/)
    call checkCMFDSurf(rfcell%infoCMFD%lsurf)

    rfcell%cut=.false.
    rfcell%cutxy=.false.
    if (lcutx) then
        rfcell%cut=.true.
        rfcell%cutxy(1)=.true.
        rfcell%slen(1)=rfcell%slen(1)*2
        rfcell%cent(1)=pax+pcx
    else
        rfcell%cent(1)=pax+pcx+rfcell%slen(1)/2
    endif
    if (lcuty) then
        rfcell%cut=.true.
        rfcell%cutxy(2)=.true.
        rfcell%slen(2)=rfcell%slen(2)*2
    endif
    rfcell%cent(2)=pay+pcy+rfcell%slen(2)/2
    rfcell%cent(3)=hzacc(icell-1)+rfcell%slen(3)/2
end subroutine

subroutine AssignGRect()

    ! assigning gap cells
    igrect=map_mc(l,m,n)%ielmt
    rgcell=>grectarray(igrect)
    rgcell%FxrIdxSt=FxrIdxSt
    rgcell%idx=imod
    rgcell%idxt=igrect
    rgcell%lmn=(/l,m,n/)
    allocate(rgcell%xs(1))
    allocate(rgcell%tally%trkl(ngMC,1))
 !   allocate(rgcell%tally(1))
    rgcell%slen(1)=mycell%geom%lx
    rgcell%slen(2)=mycell%geom%ly
    rgcell%slen(3)=hz(icell)
    rgcell%vol=rgcell%slen(1)*rgcell%slen(2)*rgcell%slen(3)

    txs=mycell%ireg(1)
    IF( lbench )THEN !--- BYSMC EDIT 16/11/09
        call AdrXsBen(rgcell%xs(1), txs)
    ELSE
        ifxr=FxrIdxSt
        iz=n
        IF( lmic )THEN
            rgcell%xs(1)=micXs(fxr2xsmap(ifxr,iz))
        ELSE
            rgcell%xs(1)=mcXs(ifxr,iz)
        ENDIF
        !ifxr=FxrIdxSt
        !myFxr=>FmInfo%FXR(ifxr,iz)
        !rfcell%xs(1)%ifxr=ifxr
        !rfcell%xs(1)%iz=iz
        !call AdrXsLib(rfcell%xs(1), myFxr, ng)
    ENDIF
!    call AllocTallySet(rgcell%tally(1), ngMC)


rgcell%infoCMFD%lmn=(/icx, icy, icell/)
    call checkCMFDSurf(rgcell%infoCMFD%lsurf)

    rgcell%cut=.false.
    rgcell%cutxy=.false.
    if (lcutx) then
        rgcell%cut=.true.
        rgcell%cutxy(1)=.true.
        rgcell%slen(1)=rgcell%slen(1)*2
        rgcell%cent(1)=pax+pcx
    else
        rgcell%cent(1)=pax+pcx+rgcell%slen(1)/2
    endif
    if (lcuty) then
        rgcell%cut=.true.
        rgcell%cutxy(2)=.true.
        rgcell%slen(2)=rgcell%slen(2)*2
    endif
    rgcell%cent(2)=pay+pcy+rgcell%slen(2)/2
    rgcell%cent(3)=hzacc(icell-1)+rgcell%slen(3)/2
end subroutine

subroutine checkCMFDSurf(lsurf)
    use PARAM
    implicit none
    logical, intent(out) :: lsurf(4,2)

    lsurf=.false.

    if(iax .eq. 1)              lsurf(WEST, RADIAL)=.true.
    if(iax .eq. myasy%nx)       lsurf(EAST, RADIAL)=.true.
    if(iay .eq. 1)              lsurf(SOUTH, RADIAL)=.true.
    if(iay .eq. myasy%ny)       lsurf(NORTH, RADIAL)=.true.
    if(icell .eq. 1)            lsurf(BOTTOM, AXIAL)=.true.
    if(icell .eq. mypin%ncell)  lsurf(TOP, AXIAL)=.true.

end subroutine



!subroutine checkCMFDSurf(ncmfdsurf, surfrd, surfax)
!    use param
!    integer, intent(out) :: ncmfdsurf(2)
!    integer, intent(out) :: surfrd(2), surfax
!    integer :: nrd, nax
!
!    nrd=0
!    nax=0
!    surfrd=0
!    surfax=0
!    if(iax .eq. 1) then
!        nrd=nrd+1
!        surfrd(nrd)=EAST
!    elseif(iax .eq. myasy%nx) then
!        nrd=nrd+1
!        surfrd(nrd)=WEST
!    endif
!    if(iay .eq. 1) then
!        nrd=nrd+1
!        surfrd(nrd)=SOUTH
!    elseif(iay .eq. myasy%ny) then
!        nrd=nrd+1
!        surfrd(nrd)=NORTH
!    endif
!    if(icell .eq. 1) then
!        nax=1
!        surfax=TOP
!    elseif(icell .eq. mypin%ncell) then
!        nax=1
!        surfax=BOTTOM
!    endif
!
!    ncmfdsurf=(/nrd, nax/)
!
!end subroutine
end subroutine

function GetVolPin(l,m,n) result(volcell)
    integer, intent(in) :: l,m,n
    real(8) :: volcell
    integer :: idx
    type(STR_RectFuel), pointer :: rfcell
    type(STR_RectGap), pointer :: rgcell
    if (map_mc(l,m,n)%telmt .eq. FRECT) then
        idx=map_mc(l,m,n)%ielmt
        rfcell=>frectarray(idx)
        volcell=rfcell%vol
    elseif (map_mc(l,m,n)%telmt .eq. GRECT) then
        idx=map_mc(l,m,n)%ielmt
        rgcell=>grectarray(idx)
        volcell=rgcell%vol
    endif
end function

! PRIVATE FUNCTION=====================================================
subroutine AdrXsBen(xs, txs)
    use BenchXs
    USE HighOrderSC
    use Xsec4MC
    implicit none
    type(XsecSet), intent(inout) :: xs
    integer, intent(in) :: txs
    integer :: ig, ng, igs
    ng=ngben
    xs%idx=txs
    xs%xst=>MacXsBen(txs)%xst
    xs%xstr=>MacXsBen(txs)%xstr
    xs%xsnf=>MacXsBen(txs)%xsnf
    xs%xskf=>MacXsBen(txs)%xskf
    xs%chi=>MacXsBen(txs)%chi
    xs%xsa=>MacXsBen(txs)%xsa

    xs%xssm=>MacXsBen(txs)%xss
    xs%xss0=>MacXsBen(txs)%xss0;
    xs%xss1=>MacXsBen(txs)%xss1;
    xs%xss2=>MacXsBen(txs)%xss2;
    xs%xss3=>MacXsBen(txs)%xss3;

    IF( Scatod .EQ. 1)THEN
        xs%rt1=>MacXsBen(txs)%rt1
    ELSEIF( ScatOd .EQ. 3)THEN
        xs%rt1=>MacXsBen(txs)%rt1
        xs%rt2=>MacXsBen(txs)%rt2
        xs%w1=>MacXsBen(txs)%w1
        xs%w2=>MacXsBen(txs)%w2
    ELSEIF( ScatOd .NE. 0 )THEN
    ENDIF
    !CALL SetAnIsoGaussian(xs, scatod)
    call AdrSctCDF(xs, txs)
    call AdrSct2G(xs, txs)
    call AdrAbsrA(xs, txs)

!    xs%SctCDF=>SctCDF(:,:,txs)
!    xs%xssm2g=>xssm2g(:,:,txs)
end subroutine
!subroutine AdrXsLib(xs, Fxr, ng)
!    USE MacXsLib_Mod,  ONLY : MacXsBase
!    USE TYPEDEF, ONLY : XsMac_Type, FXRInfo_Type
!    USE Xsec4MC
!    implicit none
!    type(XsecSet), intent(inout) :: xs
!    TYPE(XsMac_Type), Save :: XsMac
!    TYPE(FXRInfo_Type) :: FXR
!    INTEGER :: ng, g, gp
!    real(8) :: xsct
!    real(8), pointer, dimension(:,:) :: SctCDF1, xssm2g1
!    real(8), pointer, dimension(:) :: xsaa1
!    xs%idx=Fxr%imix
!    CALL MacXsBase(XSMac, Fxr, 1, ng, ng, 1._8, FALSE, FALSE, FALSE)
!    xs%xst=XsMac.XsMact
!    xs%xstr=XsMac.XsMactr
!    xs%xsnf=XsMac.XsMacnf
!    xs%xskf=XsMac.XsMackf
!    xs%chi=XsMac.chi
!    xs%xsa=XsMac.XsMaca
!    xs%xssm=XsMac.XsMacSm
!    xs%xss0=XsMac.XsMacSm
!    xs%xss1=XsMac.XsMacP1Sm
!    xs%xss2=XsMac.XsMacP2Sm
!    xs%xss3=XsMac.XsMacP3Sm
!
!
!    !call AdrSctCDF(xs, txs)
!    allocate(SctCDF1(ng,ng))
!    SctCDF1=0.
!    do g=1, ng
!        xsct=xs%xss0(g,1)
!        if (xsct<0) then
!            xsct=0.
!        endif
!
!        SctCDF1(1,g)=xsct
!        do gp=2, ng
!            xsct=xs%xss0(g,gp)
!            if (xsct<0) then
!                xsct=0.
!            endif
!
!            SctCDF1(gp,g)=SctCDF1(gp-1,g)+xsct
!        enddo
!        SctCDF1(:,g)=SctCDF1(:,g)/SctCDF1(ng,g)
!    enddo
!    xs%SctCDF=SctCDF1
!    !call AdrSct2G(xs, txs)
!    allocate(xssm2g1(ng,2))
!    xssm2g1=0.
!    do g=1, ng
!        do gp=1, ng
!            xssm2g1(g,g2gc(gp))=xssm2g1(g,g2gc(gp))+xs%xss0(g,gp)
!        enddo
!    enddo
!    xs%xssm2g=xssm2g1
!    !call AdrAbsrA(xs, txs)
!    allocate(xsaa1(ng))
!    do g=1, ng
!        do gp=1, ng
!            if (xs%xss0(g,gp)<0) then
!                if (g==gp) xsaa1(g)=-xs%xss0(g,gp)
!                xs%xss0(g,gp)=0
!            endif
!        enddo
!    enddo
!    xs%xsaa=xsaa1
!
!end subroutine
subroutine AllocTallySet(tally, ng, nr)
    implicit none
    type(MCTallySet), intent(out) :: tally
    integer :: ng, nr
    allocate(tally%trkl(ng, nr))
    tally%trkl=0.
end subroutine
! =====================================================================

function TypePin(lmn)
    implicit none
    integer, intent(in) :: lmn(3)
    integer :: TypePin
    TypePin=map_mc(lmn(1),lmn(2), lmn(3))%telmt
end function

subroutine AdrRectVac(lmn, mpVR)
    use MCDefine
    implicit none
    integer, intent(in) :: lmn(3)
    type(STR_RectVac), pointer, intent(out) :: mpVR
    integer :: idx
    idx=map_mc(lmn(1),lmn(2),lmn(3))%ielmt
    mpVR=>vrectarray(idx)
end subroutine

subroutine AdrRectFuel(lmn, mpFR)
    use MCDefine
    implicit none
    integer, intent(in) :: lmn(3)
    type(STR_RectFuel), pointer, intent(out) :: mpFR
    integer :: idx
    idx=map_mc(lmn(1),lmn(2),lmn(3))%ielmt
    mpFR=>frectarray(idx)
end subroutine

subroutine AdrRectGap(lmn, mpGR)
    use MCDefine
    implicit none
    integer, intent(in) :: lmn(3)
    type(STR_RectGap), pointer, intent(out) :: mpGR
    integer :: idx
    idx=map_mc(lmn(1),lmn(2),lmn(3))%ielmt
    mpGR=>grectarray(idx)
end subroutine

function MCNeighbor(nst)
    use PARAM
    implicit none
    integer :: MCNeighbor
    type(nstat), intent(inout) :: nst

    !integer :: AXIS(6), INC(6)
    !data AXIS /1,1,2,2,3,3/
    !data INC /-1,+1,+1,-1,+1,-1/
    integer,parameter :: LTAXIS(4,2)=RESHAPE((/2,1,2,1,3,3,0,0/), SHAPE(LTAXIS))
    integer,parameter :: LTINC(4,2)=RESHAPE((/-1,-1,+1,+1,-1,+1,0,0/), SHAPE(LTINC))
    !integer,parameter :: LTINC(4,2)=(/+1,-1,-1,+1,-1,+1,0,0/)
    integer :: axis, inc
    axis=LTAXIS(nst%surf, nst%axis)
    inc=LTINC(nst%surf, nst%axis)

    nst%lmn(axis)=nst%lmn(axis)+inc
    MCNeighbor=valid(nst%lmn(1), nst%lmn(2), nst%lmn(3))
end function

function MCNeighborRepeated(lmn, dirneigh, axisneigh)
    use PARAM
    implicit none
    integer :: MCNeighborRepeated
    integer, intent(inout) :: lmn(3)
    integer, intent(in) :: dirneigh, axisneigh
    if (axisneigh .eq. RADIAL) then
        select  case(dirneigh)
            case(SOUTH)
                lmn(2)=nym
            case(WEST)
                lmn(1)=nxm
            case(NORTH)
                lmn(2)=1
            case(EAST)
                lmn(1)=1
        end select
    else ! axial
        select  case(dirneigh)
            case(TOP)
                lmn(3)=1
            case(BOTTOM)
                lmn(3)=nzm
        end select
    endif
    MCNeighborRepeated=valid(lmn(1), lmn(2), lmn(3))
end function


end module MCmap
