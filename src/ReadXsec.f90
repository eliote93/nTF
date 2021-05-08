Module ReadXsec
!--------------------------------------------------------------------------------
!1. ReadMLX
!Function   : Read the nTRACER multi-group library (.mlb or .mla)
!2. ReadRLB
!Function   : Read the nTRACER resonance library (.rlb)
!3. ReadSLB
!Function   : Read the nTRACER spectral SPH library (.slb)
!Institute  : Seoul National University
!Authors    : Hansol Park
!Revision   : 7/07/2017  Initial version
!--------------------------------------------------------------------------------
implicit none
contains

subroutine ReadMLX(indev,mlbfile,ng,nofg,norg,ntiso,ngg)
    USE ioutil,         only : terminate, toupper, openfile
    USE XSLIB_MOD,      only : libdata, ldiso,  uhel,    enbhel,  igresb,  igrese, &
                               noghel, nofghel, norghel, nelthel, nelrhel, nchihel, &
                               enbgam, nogghel, nmlver
    USE nuclidmap_mod,  only : nuclidmap
    USE CNTL,           ONLY : nTracerCntl
    USE ALLOCS
    implicit none
    integer,intent(in) :: indev
    character*256,intent(in) :: mlbfile
    integer,intent(out) :: ng,nofg,norg,ntiso,ngg 
    integer :: iiso,ig,jg,it,igx,itx,notg,dumi,i,igg,imt, iprec
    integer :: nid,ityp,ifis,ibur,inmn,ntemp,np1temp, ipp,nx1,nx2
    real,parameter :: neutronmass = 1.008664916d0
    real :: aw,xtemp(1000)
    character*3 :: ext
    character*20 :: aid
    logical :: lbin
    TYPE(libdata),POINTER :: lib
    
    INTEGER :: nver   ! VERSION INDEX FOR multigroup library

    i = len(TRIM(mlbfile))
    ext = mlbfile(i-2:i)
    call toupper(ext)
    if (ext.eq.'MLB') then
        lbin = .true.
    elseif (ext.eq.'MLA') then
        lbin = .false.
    else
        call terminate('Extension of the multigroup library should be MLA or MLB.')
    endif

    close(indev)
    if (lbin) then
        call openfile(indev,.TRUE.,.TRUE.,.FALSE.,mlbfile)
        read(indev) dumi
        IF (dumi .GE. 100) THEN  ! Libraries controled with a integer version index (start with nver=100)
          ! ver 100 : Libaray with photon production and Neutron induced KERMA
          nver = dumi
          read(indev) ng,nofg,norg,notg,ntiso,ngg

          !!!!!!!!!
          noghel = ng
          nofghel = nofg
          norghel = norg
          nelthel = ntiso
          nelrhel = ntiso !! for future usage
          nchihel = ng    !! for future usage
          igresb = nofghel + 1
          igrese = nofghel + norghel
          nmlver = nver
          nTRACERCNTL%lExplicitKappa = nTracerCntl%lGamma
          nogghel = ngg
          !!!!!!!!!

          allocate(ldiso(ntiso))
          ! Group Boundary
          call dmalloc(enbhel,ng+1)
          call dmalloc(uhel,ng+1)
          read(indev) xtemp(1:ng+1)
          enbhel(1:ng+1) = xtemp(1:ng+1)
          ! Gamma Group boundary
          CALL dmalloc(enbgam,ngg+1)
          read(indev) xtemp(1:ngg+1)
          enbgam(1:ngg+1) = xtemp(1:ngg+1)
        
          do ig = 1, ng+1
              uhel(ig) = dlog(1.0e7_8/xtemp(ig))
          enddo

          do iiso = 1, ntiso
              lib=>ldiso(iiso)
              read(indev) dumi,nid,aw,ityp,ifis,ibur,inmn,ntemp,np1temp,ipp
              lib%nid = nid
              lib%aw = aw
              lib%mu = neutronmass * 2. / 3./ aw
              lib%ityp = ityp
              lib%ifis = ifis
              lib%ichi = ifis
              lib%ibur = ibur
              lib%inmn = inmn
              lib%ipp = ipp
              lib%ntemp = ntemp
              lib%np1temp = np1temp
              !! Memory Allocation
              call dmalloc(lib%temp,ntemp)
              call dmalloc(lib%siga,ng,ntemp)
              call dmalloc(lib%sigs,ng,ntemp)
              call dmalloc(lib%sigtr,ng,ntemp)
              call dmalloc(lib%sigstr,ng,ntemp)
              call dmalloc(lib%sigss,ng,ntemp)
              CALL dmalloc(lib%kerma_t,ng,ntemp)
              CALL dmalloc(lib%kerma_s,ng,ntemp)
              CALL dmalloc(lib%kerma_d,ng,ntemp)
              CALL dmalloc(lib%kerma_p,ng,ntemp)
              allocate(lib%sm(ng,ntemp))
              if (ifis.eq.1) then
                  call dmalloc(lib%sigf,ng,ntemp)
                  call dmalloc(lib%signf,ng,ntemp)
                  call dmalloc(lib%chi,ng)
                  call dmalloc0(lib%beta,0,6)
                  IF (nver .GE. 101) THEN
                    call dmalloc(lib%chip,ng)
                    call dmalloc(lib%chid,ng)
                    call dmalloc(lib%chidk,ng,6)
                    call dmalloc0(lib%dcy_del,1,6)
                  END IF
                  CALL dmalloc(lib%exkappa,6)     ! for explicit kappa...
                  CALL dmalloc(lib%exkappa0,6)    ! for explicit kappa...
                  CALL dmalloc(lib%kerma_f,ng,ntemp)
              endif
              if (np1temp.gt.0) then
                  call dmalloc(lib%p1temp,np1temp)
                  call dmalloc(lib%sigsp1,ng,np1temp)
                  call dmalloc(lib%sigsp2,ng,np1temp)
                  call dmalloc(lib%sigsp3,ng,np1temp)
                  allocate(lib%smp1(ng,np1temp))
                  allocate(lib%smp2(ng,np1temp))
                  allocate(lib%smp3(ng,np1temp))
              else
                  call dmalloc(lib%sigsp1,ng,ntemp)
                  lib%sigsp1 = 0._8
              endif
              call dmalloc0(lib%lamsigp,igresb,igrese)
              if (inmn.ge.2) call dmalloc(lib%sign2n,ng)
              if (inmn.eq.3) call dmalloc(lib%sign3n,ng)

              read(indev) lib%temp(1:ntemp)
              do ig = 1, ng
                  do it = 1, ntemp
                      read(indev) nx1, nx2
                      lib%sm(ig,it)%ib = nx1
                      lib%sm(ig,it)%ie = nx2
                      allocate(lib%sm(ig,it)%from(nx1:nx2))
                      if (ifis.eq.0) then
                          read(indev) lib%siga(ig,it),lib%sigtr(ig,it), lib%sigs(ig,it),&
                              lib%sm(ig,it)%from(nx1:nx2)
                      else
                          read(indev) lib%siga(ig,it),lib%sigf(ig,it), lib%signf(ig,it),&
                              lib%sigtr(ig,it),lib%sigs(ig,it),lib%sm(ig,it)%from(nx1:nx2)
                      endif
                  end do
              end do
              do it = 1, ntemp
                  lib%sigstr(:,it) = 0._8
                  do ig = 1, ng
                      do igx = 1, ng
                          if(ig.ge.lib%sm(igx,it)%ib .and. ig.le.lib%sm(igx,it)%ie) then
                              lib%sm(ig,it)%ioutsb = igx
                              exit
                          endif
                      enddo
                      do igx = ng, 1, -1
                          if(ig.ge.lib%sm(igx,it)%ib .and. ig.le.lib%sm(igx,it)%ie) then
                              lib%sm(ig,it)%ioutse = igx
                              exit
                          endif
                      enddo
                      do igx = lib%sm(ig,it)%ib, lib%sm(ig,it)%ie
                          lib%sigstr(igx,it) = lib%sigstr(igx,it) + lib%sm(ig,it)%from(igx)
                      enddo
                  enddo
              enddo
              do it = 1, ntemp
                  do ig = 1, ng
                      lib%sigss(ig,it) = lib%sm(ig,it)%from(ig) + (lib%sigs(ig,it) - lib%sigstr(ig,it))
                  enddo
              enddo
              if (np1temp.gt.0) then
                  read(indev) lib%p1temp(1:np1temp)
                  do ig = 1, ng
                      do it = 1, np1temp
                          read(indev) nx1, nx2
                          lib%smp1(ig,it)%ib = nx1
                          lib%smp1(ig,it)%ie = nx2
                          allocate(lib%smp1(ig,it)%from(nx1:nx2))
                          read(indev) lib%sigsp1(ig,it), lib%smp1(ig,it)%from(nx1:nx2)
                      end do
                  end do
                  do ig = 1, ng
                      do it = 1, np1temp
                          read(indev) nx1, nx2
                          lib%smp2(ig,it)%ib = nx1
                          lib%smp2(ig,it)%ie = nx2
                          allocate(lib%smp2(ig,it)%from(nx1:nx2))
                          read(indev) lib%sigsp2(ig,it), lib%smp2(ig,it)%from(nx1:nx2)
                      end do
                  end do
                  do ig = 1, ng
                      do it = 1, np1temp
                          read(indev) nx1, nx2
                          lib%smp3(ig,it)%ib = nx1
                          lib%smp3(ig,it)%ie = nx2
                          allocate(lib%smp3(ig,it)%from(nx1:nx2))
                          read(indev) lib%sigsp3(ig,it), lib%smp3(ig,it)%from(nx1:nx2)
                      end do
                  end do
              else
                  call dmalloc(lib%sigsp1,ng,ntemp)
                  lib%sigsp1 = 1._8
              endif
              read(indev) lib%sigp, lib%lamsigp1G
              read(indev) lib%lamsigp(igresb:igrese)
              if (ifis.eq.1) THEN
                  read(indev) lib%chi(1:ng)
                  IF (nver.GE.101) THEN
                    read(indev) lib%chip(1:ng)
                    read(indev) lib%chid(1:ng)
                    DO iprec = 1, 6
                      read(indev) lib%chidk(1:ng, iprec)
                    END DO 
                  END IF
                  read(indev) lib%kappa,lib%kappa0
                  read(indev) lib%exkappa(1:6)
                  read(indev) lib%exkappa0(1:6)
                  IF (nver.GE.101) read(indev) lib%dcy_del(1:6)
                  read(indev) lib%beta(0:6)
              endif
              IF (nver.GE.101) THEN
                ALLOCATE(lib%InvV(ig,it))
                DO ig = 1, ng
                  DO it = 1, ntemp
                    READ(indev) lib%InvV(ig,it)
                  END DO
                END DO
              END IF
              READ(indev) lib%capkappa, lib%capQ, lib%n2nQ, lib%n3nQ
              if (ibur.eq.1) read(indev) lib%dcy
              if (inmn.ge.2) read(indev) lib%sign2n(1:ng)
              if (inmn.eq.3) read(indev) lib%sign3n(1:ng)
              IF (ifis.eq.0) THEN
                DO ig = 1,ng
                  DO it = 1, ntemp
                    READ(indev) lib%kerma_t(ig,it), lib%kerma_s(ig,it), lib%kerma_d(ig,it), lib%kerma_p(ig,it)
                  END DO
                END DO
              ELSE
                DO ig = 1,ng
                  DO it = 1, ntemp
                    READ(indev) lib%kerma_t(ig,it), lib%kerma_s(ig,it), lib%kerma_d(ig,it), lib%kerma_p(ig,it), lib%kerma_f(ig,it)
                  END DO
                END DO
              END IF
              IF (ipp .eq.1) THEN
                ALLOCATE(lib%ppm(0:4))
                READ(indev) lib%lphoton(1:4)
                DO imt = 1, 4
                  IF (.NOT.lib%lphoton(imt)) CYCLE
                  READ(indev) lib%ppm(imt)%iglow, lib%ppm(imt)%igup
                  ALLOCATE(lib%ppm(imt)%mat(lib%ppm(imt)%iglow:lib%ppm(imt)%igup,ntemp))
                  DO igg = lib%ppm(imt)%iglow, lib%ppm(imt)%igup
                    DO it = 1, ntemp
                      READ(indev) nx1, nx2
                      lib%ppm(imt)%mat(igg,it)%ib = nx1
                      lib%ppm(imt)%mat(igg,it)%ie = nx2
                      ALLOCATE(lib%ppm(imt)%mat(igg,it)%from(nx1:nx2))
                      READ(indev) lib%ppm(imt)%mat(igg,it)%from(nx1:nx2)
                    END DO
                  END DO
                END DO
              END IF
          enddo
        ELSE !IF (dumi.EQ.47) THEN ! library without photon generation and KERMA
          nTracerCntl%lGamma = .FALSE.
          REWIND(indev)
        read(indev) ng,nofg,norg,notg,ntiso

        !!!!!!!!!
        noghel = ng
        nofghel = nofg
        norghel = norg
        nelthel = ntiso
        nelrhel = ntiso !! for future usage
        nchihel = ng    !! for future usage
        igresb = nofghel + 1
        igrese = nofghel + norghel
        !!!!!!!!!

        allocate(ldiso(ntiso))
        ! Group Boundary
        call dmalloc(enbhel,ng+1)
        call dmalloc(uhel,ng+1)
        read(indev) xtemp(1:ng+1)
        enbhel(1:ng+1) = xtemp(1:ng+1)
        do ig = 1, ng+1
            uhel(ig) = dlog(1.0e7_8/xtemp(ig))
        enddo

        do iiso = 1, ntiso
            lib=>ldiso(iiso)

            read(indev) dumi,nid,aw,ityp,ifis,ibur,inmn,ntemp,np1temp
            lib%nid = nid
            lib%aw = aw
            lib%mu = neutronmass * 2. / 3./ aw
            lib%ityp = ityp
            lib%ifis = ifis
            lib%ichi = ifis
            lib%ibur = ibur
            lib%inmn = inmn
            lib%ntemp = ntemp
            lib%np1temp = np1temp
            !! Memory Allocation
            call dmalloc(lib%temp,ntemp)
            call dmalloc(lib%siga,ng,ntemp)
            call dmalloc(lib%sigs,ng,ntemp)
            call dmalloc(lib%sigtr,ng,ntemp)
            call dmalloc(lib%sigstr,ng,ntemp)
            call dmalloc(lib%sigss,ng,ntemp)
            allocate(lib%sm(ng,ntemp))
            if (ifis.eq.1) then
                call dmalloc(lib%sigf,ng,ntemp)
                call dmalloc(lib%signf,ng,ntemp)
                call dmalloc(lib%chi,ng)
                call dmalloc0(lib%beta,0,6)
            endif
            if (np1temp.gt.0) then
                call dmalloc(lib%p1temp,np1temp)
                call dmalloc(lib%sigsp1,ng,np1temp)
                call dmalloc(lib%sigsp2,ng,np1temp)
                call dmalloc(lib%sigsp3,ng,np1temp)
                allocate(lib%smp1(ng,np1temp))
                allocate(lib%smp2(ng,np1temp))
                allocate(lib%smp3(ng,np1temp))
            else
                call dmalloc(lib%sigsp1,ng,ntemp)
                lib%sigsp1 = 0._8
            endif
            call dmalloc0(lib%lamsigp,igresb,igrese)
            if (inmn.ge.2) call dmalloc(lib%sign2n,ng)
            if (inmn.eq.3) call dmalloc(lib%sign3n,ng)

            read(indev) lib%temp(1:ntemp)
            do ig = 1, ng
                do it = 1, ntemp
                    read(indev) nx1, nx2
                    lib%sm(ig,it)%ib = nx1
                    lib%sm(ig,it)%ie = nx2
                    allocate(lib%sm(ig,it)%from(nx1:nx2))
                    if (ifis.eq.0) then
                        read(indev) lib%siga(ig,it),lib%sigtr(ig,it), lib%sigs(ig,it),&
                            lib%sm(ig,it)%from(nx1:nx2)
                    else
                        read(indev) lib%siga(ig,it),lib%sigf(ig,it), lib%signf(ig,it),&
                            lib%sigtr(ig,it),lib%sigs(ig,it),lib%sm(ig,it)%from(nx1:nx2)
                    endif
                end do
            end do
            do it = 1, ntemp
                lib%sigstr(:,it) = 0._8
                do ig = 1, ng
                    do igx = 1, ng
                        if(ig.ge.lib%sm(igx,it)%ib .and. ig.le.lib%sm(igx,it)%ie) then
                            lib%sm(ig,it)%ioutsb = igx
                            exit
                        endif
                    enddo
                    do igx = ng, 1, -1
                        if(ig.ge.lib%sm(igx,it)%ib .and. ig.le.lib%sm(igx,it)%ie) then
                            lib%sm(ig,it)%ioutse = igx
                            exit
                        endif
                    enddo
                    do igx = lib%sm(ig,it)%ib, lib%sm(ig,it)%ie
                        lib%sigstr(igx,it) = lib%sigstr(igx,it) + lib%sm(ig,it)%from(igx)
                    enddo
                enddo
            enddo
            do it = 1, ntemp
                do ig = 1, ng
                    lib%sigss(ig,it) = lib%sm(ig,it)%from(ig) + (lib%sigs(ig,it) - lib%sigstr(ig,it))
                enddo
            enddo
            if (np1temp.gt.0) then
                read(indev) lib%p1temp(1:np1temp)
                do ig = 1, ng
                    do it = 1, np1temp
                        read(indev) nx1, nx2
                        lib%smp1(ig,it)%ib = nx1
                        lib%smp1(ig,it)%ie = nx2
                        allocate(lib%smp1(ig,it)%from(nx1:nx2))
                        read(indev) lib%sigsp1(ig,it), lib%smp1(ig,it)%from(nx1:nx2)
                    end do
                end do
                do ig = 1, ng
                    do it = 1, np1temp
                        read(indev) nx1, nx2
                        lib%smp2(ig,it)%ib = nx1
                        lib%smp2(ig,it)%ie = nx2
                        allocate(lib%smp2(ig,it)%from(nx1:nx2))
                        read(indev) lib%sigsp2(ig,it), lib%smp2(ig,it)%from(nx1:nx2)
                    end do
                end do
                do ig = 1, ng
                    do it = 1, np1temp
                        read(indev) nx1, nx2
                        lib%smp3(ig,it)%ib = nx1
                        lib%smp3(ig,it)%ie = nx2
                        allocate(lib%smp3(ig,it)%from(nx1:nx2))
                        read(indev) lib%sigsp3(ig,it), lib%smp3(ig,it)%from(nx1:nx2)
                    end do
                end do
            else
                call dmalloc(lib%sigsp1,ng,ntemp)
                lib%sigsp1 = 1._8
            endif
            read(indev) lib%sigp, lib%lamsigp1G
            read(indev) lib%lamsigp(igresb:igrese)
            if (ifis.eq.1) THEN
                read(indev) lib%chi(1:ng)
                read(indev) lib%kappa,lib%kappa0
                read(indev) lib%beta(0:6)
            endif
            if (ibur.eq.1) read(indev) lib%dcy
            if (inmn.ge.2) read(indev) lib%sign2n(1:ng)
            if (inmn.eq.3) read(indev) lib%sign3n(1:ng)
        enddo
        !ELSE
        !  STOP 'Currently unassigned MLB -- ReadXsec.f90'
        END IF
        close(indev)
    else ! for ASCII mla
        call openfile(indev,.TRUE.,.FALSE.,.FALSE.,mlbfile)
        read(indev,*)
        !read(indev,'(i15,6i5)') ng,nofg,norg,notg,ntiso,ngg,nver
        read(indev,'(i15,6i5)') nver,ng,nofg,norg,notg,ntiso,ngg

        !!!!!!!!!
        noghel=ng
        nofghel=nofg
        norghel=norg
        nelthel=ntiso
        nelrhel=ntiso !! for future usage
        nchihel=ng    !! for future usage
        igresb=nofghel+1
        igrese=nofghel+norghel
        nogghel = ngg
        !!!!!!!!!

        allocate(ldiso(ntiso))
        ! Group Boundary
        call dmalloc(enbhel,ng+1)
        call dmalloc(uhel,ng+1)
        read(indev,*)
        read(indev,'(10ES14.7)') (xtemp(ig),ig=1,ng+1)
        enbhel(1:ng+1) = xtemp(1:ng+1)
        ! Gamma Group boundary
        CALL dmalloc(enbgam,ngg+1)
        read(indev,*)
        read(indev,'(10ES14.7)') xtemp(1:ngg+1)
        enbgam(1:ngg+1) = xtemp(1:ngg+1)
        do ig = 1, ng+1
            uhel(ig) = dlog(1.0e7_8/xtemp(ig))
        enddo
        read(indev,*)
        do iiso = 1, ntiso
            read(indev,*)
        enddo
        do iiso = 1, ntiso
            lib=>ldiso(iiso)
            read(indev,*) aid,i,nid,aw,ityp,ifis,ibur,inmn,ntemp,np1temp,ipp,aid
            lib%nid = nid
            lib%aw = aw
            lib%mu = neutronmass * 2. / 3./ aw
            lib%ityp = ityp
            lib%ifis = ifis
            lib%ibur = ibur
            lib%inmn = inmn
            lib%ntemp = ntemp
            lib%np1temp = np1temp
            lib%ipp = ipp
            lib%aid=aid

            !! Memory Allocation
            call dmalloc(lib%temp,ntemp)
            call dmalloc(lib%siga,ng,ntemp)
            call dmalloc(lib%sigs,ng,ntemp)
            call dmalloc(lib%sigtr,ng,ntemp)
            call dmalloc(lib%sigstr,ng,ntemp)
            call dmalloc(lib%sigss,ng,ntemp)
            allocate(lib%sm(ng,ntemp))
            if (ifis.eq.1) then
                call dmalloc(lib%sigf,ng,ntemp)
                call dmalloc(lib%signf,ng,ntemp)
                call dmalloc(lib%chi,ng)
                call dmalloc0(lib%dcy_del,1,6)
                call dmalloc0(lib%beta,0,6)
                IF (nver.GE.101) THEN
                  call dmalloc(lib%chip,ng)
                  call dmalloc(lib%chid,ng)
                  call dmalloc(lib%chidk,ng,6)
                  call dmalloc0(lib%dcy_del,1,6)
                END IF
                CALL dmalloc(lib%exkappa,6)     ! for explicit kappa...
                CALL dmalloc(lib%exkappa0,6)    ! for explicit kappa...
            endif
            if (np1temp.gt.0) then
                call dmalloc(lib%p1temp,np1temp)
                call dmalloc(lib%sigsp1,ng,np1temp)
                call dmalloc(lib%sigsp2,ng,np1temp)
                call dmalloc(lib%sigsp3,ng,np1temp)
                allocate(lib%smp1(ng,np1temp))
                allocate(lib%smp2(ng,np1temp))
                allocate(lib%smp3(ng,np1temp))
            endif
            call dmalloc0(lib%lamsigp,igresb,igrese)
            if (inmn.ge.2) call dmalloc(lib%sign2n,ng)
            if (inmn.eq.3) call dmalloc(lib%sign3n,ng)

            read(indev,*) ! 'TP0+'
            read(indev,*) lib%temp(1:lib%ntemp)
            read(indev,*) ! 'XSD+'
            do ig=1,ng
                do it=1,ntemp
                    if (ifis.eq.0) then
                        read(indev,*) igx,itx,lib%siga(ig,it),lib%sigtr(ig,it), &
                                      lib%sigs(ig,it),nx1,nx2,(xtemp(jg),jg=nx1,nx2)
                    else
                        read(indev,*) igx,itx,lib%siga(ig,it),lib%sigf(ig,it), &
                            lib%signf(ig,it),lib%sigtr(ig,it),lib%sigs(ig,it), &
                            nx1,nx2,(xtemp(jg),jg=nx1,nx2)
                    endif
                    lib%sm(ig,it)%ib=nx1
                    lib%sm(ig,it)%ie=nx2
                    allocate(lib%sm(ig,it)%from(nx1:nx2))
                    do jg=nx1,nx2
                       lib%sm(ig,it)%from(jg)=xtemp(jg)
                    enddo
                end do
            end do
            do it=1,ntemp
                lib%sigstr(:,it)=0._8
                do ig=1,ng
                    do igx=1,ng
                        if(ig.ge.lib%sm(igx,it)%ib .and. ig.le.lib%sm(igx,it)%ie) then
                            lib%sm(ig,it)%ioutsb=igx
                            exit
                        endif
                    enddo
                    do igx=ng,1,-1
                        if(ig.ge.lib%sm(igx,it)%ib .and. ig.le.lib%sm(igx,it)%ie) then
                            lib%sm(ig,it)%ioutse=igx
                            exit
                        endif
                    enddo
                    do igx=lib%sm(ig,it)%ib,lib%sm(ig,it)%ie
                        lib%sigstr(igx,it)=lib%sigstr(igx,it)+lib%sm(ig,it)%from(igx)
                    enddo
                enddo
            enddo
            do it = 1, ntemp
                do ig = 1, ng
                    lib%sigss(ig,it) = lib%sm(ig,it)%from(ig) + (lib%sigs(ig,it) - lib%sigstr(ig,it))
                enddo
            enddo
            if (np1temp.gt.0) then
                read(indev,*)
                read(indev,*) (lib%p1temp(it),it=1,np1temp)
                read(indev,*)
                do ig=1,ng
                    do it=1,np1temp
                        read(indev,*) igx,itx,lib%sigsp1(ig,it),nx1,nx2,(xtemp(jg),jg=nx1,nx2)
                        lib%smp1(ig,it)%ib=nx1
                        lib%smp1(ig,it)%ie=nx2
                        allocate(lib%smp1(ig,it)%from(nx1:nx2))
                        do jg=nx1,nx2
                           lib%smp1(ig,it)%from(jg)=xtemp(jg)
                        enddo
                    end do
                end do
                read(indev,*)
                do ig=1,ng
                    do it=1,np1temp
                        read(indev,*) igx,itx,lib%sigsp2(ig,it),nx1,nx2,(xtemp(jg),jg=nx1,nx2)
                        lib%smp2(ig,it)%ib=nx1
                        lib%smp2(ig,it)%ie=nx2
                        allocate(lib%smp2(ig,it)%from(nx1:nx2))
                        do jg=nx1,nx2
                           lib%smp2(ig,it)%from(jg)=xtemp(jg)
                        enddo
                    end do
                end do
                read(indev,*)
                do ig=1,ng
                    do it=1,np1temp
                        read(indev,*) igx,itx,lib%sigsp3(ig,it),nx1,nx2,(xtemp(jg),jg=nx1,nx2)
                        lib%smp3(ig,it)%ib=nx1
                        lib%smp3(ig,it)%ie=nx2
                        allocate(lib%smp3(ig,it)%from(nx1:nx2))
                        do jg=nx1,nx2
                           lib%smp3(ig,it)%from(jg)=xtemp(jg)
                        enddo
                    end do
                end do
            else
                call dmalloc(lib%sigsp1,ng,ntemp)
                lib%sigsp1 = 1._8
            endif
            read(indev,*)
            read(indev,*) lib%sigp,lib%lamsigp1g
            read(indev,*) lib%lamsigp(igresb:igrese)
            if (ifis.eq.1) then
                read(indev,*)
                read(indev,'(10ES14.6)') lib%chi(1:ng)
                IF (nver .GE. 101) THEN
                  read(indev,'(10ES14.6)') lib%chip(1:ng)
                  read(indev,'(10ES14.6)') lib%chid(1:ng)
                  DO iprec = 1,6
                    read(indev,'(10ES14.6)') lib%chidk(1:ng,iprec)
                  END DO 
                END IF
                read(indev,*)
                read(indev,'(2ES14.6)') lib%kappa,lib%kappa0
                read(indev,*)
                read(indev,'(6ES14.6)') lib%exkappa(1:6)
                read(indev,*)
                read(indev,'(6ES14.6)') lib%exkappa0(1:6)
                read(indev,*)
                IF (nver .GE. 101) read(indev,'(6ES14.6)') lib%dcy_del(0:6)
                read(indev,'(7ES14.6)') lib%beta(0:6)
            endif
            IF (nver .GE. 101) THEN
              read(indev,*)
              ALLOCATE(lib%InvV(ng,ntemp))
              DO ig = 1, ng
                READ(indev,*) (lib%InvV(ig,it),it=1,ntemp)
              END DO
            END IF
            
            read(indev,*)
            READ(indev,'(1ES14.6)') lib%capkappa
            read(indev,*)
            READ(indev,'(3ES14.6)') lib%capQ, lib%n2nQ, lib%n3nQ
            if (ibur.eq.1) then
                read(indev,*)
                read(indev,*) lib%dcy
            endif
            if (inmn.ge.2) then
                read(indev,*)
                read(indev,*) lib%sign2n(1:ng)
            endif
            if (inmn.eq.3) then
                read(indev,*)
                read(indev,*) lib%sign3n(1:ng)
            endif
            IF (ipp .eq.1) THEN  ! Photon Production Data
              READ(indev, *)
              ALLOCATE(lib%ppm(0:4))
              READ(indev, '(4L2)') (lib%lphoton(imt), imt = 1,4)
              DO imt = 1, 4
                IF (.NOT.lib%lphoton(imt)) CYCLE
                READ(indev, '(3I6)') igx, lib%ppm(imt)%iglow, lib%ppm(imt)%igup
                ALLOCATE(lib%ppm(imt)%mat(lib%ppm(imt)%iglow:lib%ppm(imt)%igup,ntemp))
                DO igg = lib%ppm(imt)%iglow, lib%ppm(imt)%igup
                  DO it = 1, ntemp
                    read(indev,*) igx,itx,nx1,nx2,(xtemp(jg),jg=nx1,nx2)
                    lib%ppm(imt)%mat(igg,it)%ib = nx1
                    lib%ppm(imt)%mat(igg,it)%ie = nx2
                    ALLOCATE(lib%ppm(imt)%mat(igg,it)%from(nx1:nx2))
                    do jg=nx1,nx2
                       lib%ppm(imt)%mat(igg,it)%from(jg)=xtemp(jg)
                    enddo
                  END DO
                END DO
              END DO
            END IF
        enddo
        close(indev)

    endif

    !- Photon PRODUCTION RANGE
    !      RANGE OF PHORON GROUP CORRESPONDING NEUTRON GROUP
    DO iiso = 1, ntiso
      lib => ldiso(iiso)
      IF (lib%ipp .EQ. 0) CYCLE
      DO imt = 1, 4
        IF (.NOT. lib%lphoton(imt)) CYCLE
        CALL dmalloc(lib%ppm(imt)%ioutpb,ng,lib%ntemp)
        CALL dmalloc(lib%ppm(imt)%ioutpe,ng,lib%ntemp)
        lib%ppm(imt)%ioutpb = lib%ppm(imt)%igup
        lib%ppm(imt)%ioutpe = lib%ppm(imt)%iglow
        DO it = 1, lib%ntemp
          DO ig = 1, ng
            DO igg = lib%ppm(imt)%iglow, lib%ppm(imt)%igup         ! ig -> igg, beginning
              IF(ig .GE. lib%ppm(imt)%mat(igg,it)%ib .AND. ig .LE. lib%ppm(imt)%mat(igg,it)%ie) THEN
                lib%ppm(imt)%ioutpb(ig,it) = igg  ! upper(energy) group boundary from igg
                EXIT
              END IF
            END DO
            DO igg=lib%ppm(imt)%igup,lib%ppm(imt)%iglow,-1         ! ig -> igg, ending
              IF(ig .GE. lib%ppm(imt)%mat(igg,it)%ib .AND. ig .LE. lib%ppm(imt)%mat(igg,it)%ie) THEN
                lib%ppm(imt)%ioutpe(ig,it) = igg  ! upper(energy) group boundary from igg
                EXIT
              END IF
            END DO
          END DO ! NEUTRON GROUP LOOP (ig)
        END DO ! TEMPERATURE LOOP (it)
      END DO ! REACTION TYPE LOOP (imt)
    END DO
    !!! mapping !!!
    call nuclidmap

end subroutine

subroutine ReadRLB(indev,rlbfile)
    USE ioutil,    ONLY : terminate, openfile
    USE XSLIB_MOD, ONLY : libdata, ldiso, norghel, nreshel, nelthel, igresb, igrese, idres_cor, nlvmax, mlgdata0, nmaxz
    USE ALLOCS
    USE cntl,      ONLY : nTracerCntl
    implicit none
    integer,intent(in) :: indev
    character*256,intent(in) :: rlbfile
    integer :: nrg, f_nmaclv, f_nmaclv1G, c_nmaclv, c_nmaclv1G, dumi
    integer :: irl,nid,njd,ireg,ifis,ntemp,nlv,nrif,nsig0,ir,nidrif,nrat,i,j,xxx,mlbnid,list(1000),is,ns,i1,isig,it,ig,nlvflx
    real,allocatable :: dum1_(:),dum2_(:,:),dum3_(:,:,:),dum4_(:,:,:,:)
    real(4),allocatable :: dum1(:),dum3(:,:,:),dum4(:,:,:,:)
    TYPE(libdata),POINTER :: lib,ljb

    close(indev)
    call openfile(indev,.TRUE.,.TRUE.,.FALSE.,rlbfile)
    read(indev) nrg,nreshel
    if (nrg.ne.norghel) then
        call terminate('# of resonance groups in resonance library differs from that in multigroup library.')
    endif
    do irl=1,nreshel
        IF (irl.eq.8) THEN
            CONTINUE
        ENDIF
        read(indev) nid,ireg,ifis
        ! mapping with mlb
        ns=0
        list=0
        do i=1,nelthel
            mlbnid=ldiso(i)%nid
            xxx=mod(mlbnid,1000)
            if (xxx.gt.500) mlbnid=mlbnid-500
            if (nid.eq.mlbnid) then
                ns=ns+1
                list(ns)=i
            endif
        enddo
        !
        if (ns.eq.0) then  ! IF "nid" in rlb doesn't exist in mlb
            read(indev) ntemp,nlv,nlvflx
            allocate(dum1_(ntemp))
            read(indev) dum1_(1:ntemp)
            deallocate(dum1_)
            allocate(dum2_(1:nlv,igresb:igrese),dum3_(1:nlv,1:ntemp,igresb:igrese))
            read(indev) dum2_(1:nlv,igresb:igrese),dum3_(1:nlv,1:ntemp,igresb:igrese)
            if (ifis.ne.0) then
                read(indev) dum2_(1:nlv,igresb:igrese),dum3_(1:nlv,1:ntemp,igresb:igrese)
            endif
            deallocate(dum2_,dum3_)
            allocate(dum2_(1:nlvflx,igresb:igrese))
            read(indev) dum2_(1:nlvflx,igresb:igrese)
            deallocate(dum2_)
            read(indev) ntemp,nrif,nsig0
            allocate(dum1(1:ntemp))
            read(indev) dum1(1:ntemp)
            deallocate(dum1)
            allocate(dum1(1:nsig0))
            read(indev) dum1(1:nsig0)
            deallocate(dum1)
            allocate(dum3(1:nsig0,igresb:igrese,1:ntemp))
            read(indev) dum3(1:nsig0,igresb:igrese,1:ntemp)
            read(indev) dum3(1:nsig0,igresb:igrese,1:ntemp)
            allocate(dum1(igresb:igrese))
            read(indev) dum1(igresb:igrese)
            deallocate(dum1)
            if (ifis.ne.0) then
                read(indev) dum3(1:nsig0,igresb:igrese,1:ntemp)
            endif
            deallocate(dum3)
            do ir=1,nrif
                read(indev) nidrif,nrat
                allocate(dum1(1:nrat))
                read(indev) dum1(1:nrat)
                deallocate(dum1)
                allocate(dum4(1:nrat,1:nsig0,igresb:igrese,1:ntemp))
                read(indev) dum4(1:nrat,1:nsig0,igresb:igrese,1:ntemp)
                read(indev) dum4(1:nrat,1:nsig0,igresb:igrese,1:ntemp)
                if (ifis.ne.0) then
                    read(indev) dum4(1:nrat,1:nsig0,igresb:igrese,1:ntemp)
                endif
                deallocate(dum4)
            enddo
            cycle
        endif
        i1=list(1)
        do is=1,ns
            i=list(is)
            lib=>ldiso(i)
            lib%lreso=.true.
            if (ireg.eq.0) lib%lfuel=.true.
            if (ireg.eq.1) then
                xxx=mod(ldiso(i)%nid,1000)
                if (xxx.le.500) then
                    lib%lclad=.true.
                else
                    lib%lreso=.false.
                endif
            endif
            if (lib%ityp.eq.3) then
                if (ifis.eq.0) then
                    call terminate('Fissile isotope in multigroup library is not fissile in resonance library.')
                endif
            else
                ifis=0
            endif
            if (is.eq.1) then
                read(indev) ntemp,nlv,nlvflx
            else
                ntemp=ldiso(i1)%nrtemp
                nlv=ldiso(i1)%nlv
                nlvflx=ldiso(i1)%nlvflx
            endif
            lib%nrtemp=ntemp
            lib%nlv=nlv
            nlvmax=max(nlvmax,nlv)
            lib%nlvflx=nlvflx
            call dmalloc(lib%rtemp,ntemp)
            call dmalloc(lib%rtempsq,ntemp)
            call dmalloc0(lib%lvabs,1,nlv,igresb,igrese)
            call dmalloc0(lib%wgtabs,1,nlv,1,ntemp,igresb,igrese)
            if (is.eq.1) then
                read(indev) lib%rtemp(1:ntemp)
                lib%rtempsq(1:ntemp)=dsqrt(lib%rtemp(1:ntemp))
                read(indev) lib%lvabs(1:nlv,igresb:igrese),lib%wgtabs(1:nlv,1:ntemp,igresb:igrese)
            else
                lib%rtemp(1:ntemp)=ldiso(i1)%rtemp(1:ntemp)
                lib%lvabs(1:nlv,igresb:igrese)=ldiso(i1)%lvabs(1:nlv,igresb:igrese)
                lib%wgtabs(1:nlv,1:ntemp,igresb:igrese)=ldiso(i1)%wgtabs(1:nlv,1:ntemp,igresb:igrese)
                lib%rtempsq(1:ntemp)=ldiso(i1)%rtempsq(1:ntemp)
            endif
            if (ifis.ne.0) then
                call dmalloc0(lib%lvfis,1,nlv,igresb,igrese)
                call dmalloc0(lib%wgtfis,1,nlv,1,ntemp,igresb,igrese)
                if (is.eq.1) then
                    read(indev) lib%lvfis(1:nlv,igresb:igrese),lib%wgtfis(1:nlv,1:ntemp,igresb:igrese)
                else
                    lib%lvfis(1:nlv,igresb:igrese)=ldiso(i1)%lvfis(1:nlv,igresb:igrese)
                    lib%wgtfis(1:nlv,1:ntemp,igresb:igrese)=ldiso(i1)%wgtfis(1:nlv,1:ntemp,igresb:igrese)
                endif
            endif
            call dmalloc0(lib%lvflx,1,nlvflx,igresb,igrese)
            if (is.eq.1) then
                read(indev) lib%lvflx(1:nlvflx,igresb:igrese)
            else
                lib%lvflx(1:nlvflx,igresb:igrese)=ldiso(i1)%lvflx(1:nlvflx,igresb:igrese)
            endif
            if (is.eq.1) then
                read(indev) ntemp,nrif,nsig0
            else
                nrif=ldiso(i1)%nrif
                nsig0=ldiso(i1)%nsig0
            endif
            lib%nrif=nrif
            allocate(lib%rif(nrif),lib%rifid(nelthel))
            lib%rifid=0
            if (is.eq.1) then
                lib%nsig0=nsig0
            else
                lib%nsig0=ldiso(i1)%nsig0
            endif
            call dmalloc(lib%sig0sq,nsig0)
            call dmalloc0(lib%ri_a,1,nsig0,igresb,igrese,1,ntemp)
            call dmalloc0(lib%xsalog,1,nsig0,igresb,igrese,1,ntemp)
            call dmalloc0(lib%xss,1,nsig0,igresb,igrese,1,ntemp)
            call dmalloc0(lib%ssf,igresb,igrese)
!            print*, lib%nid
            if (is.eq.1) then
                allocate(dum1(ntemp))
                read(indev) dum1(1:ntemp)
                lib%rtemp(1:ntemp)=dum1(1:ntemp)
                deallocate(dum1)
                allocate(dum1(nsig0))
                read(indev) dum1(1:nsig0)
                lib%sig0sq(1:nsig0)=dum1(1:nsig0)
                deallocate(dum1)
                allocate(dum3(1:nsig0,igresb:igrese,1:ntemp))
                read(indev) dum3(1:nsig0,igresb:igrese,1:ntemp)
                lib%xsalog(1:nsig0,igresb:igrese,1:ntemp)=dum3(1:nsig0,igresb:igrese,1:ntemp)
                do it=1,ntemp
                  do ig=igresb,igrese
!                    print*, it, ig
                    do isig=1,nsig0
                      lib%ri_a(isig,ig,it)=lib%xsalog(isig,ig,it)*(lib%sig0sq(isig)+lib%lamsigp(ig))/(lib%xsalog(isig,ig,it)+lib%sig0sq(isig)+lib%lamsigp(ig))
                      lib%xsalog(isig,ig,it)=log(lib%xsalog(isig,ig,it))
                    enddo
                    do isig=1,nsig0-1
                      if (lib%xsalog(isig,ig,it).GT.lib%xsalog(isig+1,ig,it)) exit
                    end do
                    if (isig.LT.nsig0) then
                      lib%xsalog(:,ig,it) = -10.
                    end if
!                    print*, lib%xsalog(:,ig,it)
                  enddo
                enddo
                lib%sig0sq(1:nsig0)=dsqrt(lib%sig0sq(1:nsig0))
                read(indev) dum3(1:nsig0,igresb:igrese,1:ntemp)
                lib%xss(1:nsig0,igresb:igrese,1:ntemp)=dum3(1:nsig0,igresb:igrese,1:ntemp)
                deallocate(dum3)
                allocate(dum1(igresb:igrese))
                read(indev) dum1(igresb:igrese)
                lib%ssf(igresb:igrese)=dum1(igresb:igrese)
                deallocate(dum1)
                if (ifis.ne.0) then
                    allocate(dum3(1:nsig0,igresb:igrese,1:ntemp))
                    read(indev) dum3(1:nsig0,igresb:igrese,1:ntemp)
                deallocate(dum3)
              endif
            else
                lib%rtemp(1:ntemp)=ldiso(i1)%rtemp(1:ntemp)
                lib%sig0sq(1:nsig0)=ldiso(i1)%sig0sq(1:nsig0)
                lib%xsalog(1:nsig0,igresb:igrese,1:ntemp)=ldiso(i1)%xsalog(1:nsig0,igresb:igrese,1:ntemp)
                lib%ri_a(1:nsig0,igresb:igrese,1:ntemp)=ldiso(i1)%ri_a(1:nsig0,igresb:igrese,1:ntemp)
                lib%xss(1:nsig0,igresb:igrese,1:ntemp)=ldiso(i1)%xss(1:nsig0,igresb:igrese,1:ntemp)
                lib%ssf(igresb:igrese)=ldiso(i1)%ssf(igresb:igrese)
            endif
            do ir=1,nrif
                if (is.eq.1) then
                    read(indev) nidrif,nrat
                else
                    nidrif=ldiso(i1)%rif(ir)%nid
                    nrat=ldiso(i1)%rif(ir)%nrat
                endif
                lib%rif(ir)%nid=nidrif
                lib%rif(ir)%nrat=nrat
                call dmalloc(lib%rif(ir)%rat,nrat)
                call dmalloc(lib%rif(ir)%ratlog,nrat)
                if (is.eq.1) then
                    allocate(dum1(nrat))
                    read(indev) dum1(1:nrat)
                    lib%rif(ir)%rat(1:nrat)=dum1(1:nrat)
                    lib%rif(ir)%ratlog(1:nrat)=log(lib%rif(ir)%rat(1:nrat))
                    deallocate(dum1)
                else
                    lib%rif(ir)%rat(1:nrat)=ldiso(i1)%rif(ir)%rat(1:nrat)
                    lib%rif(ir)%ratlog(1:nrat)=ldiso(i1)%rif(ir)%ratlog(1:nrat)
                endif
                allocate(lib%rif(ir)%abs(1:nrat,1:nsig0,igresb:igrese,1:ntemp))
                allocate(lib%rif(ir)%sct(1:nrat,1:nsig0,igresb:igrese,1:ntemp))
                if (is.eq.1) then
                    allocate(dum4(1:nrat,1:nsig0,igresb:igrese,1:ntemp))
                    read(indev) dum4(1:nrat,1:nsig0,igresb:igrese,1:ntemp)
                    lib%rif(ir)%abs(1:nrat,1:nsig0,igresb:igrese,1:ntemp)=dum4(1:nrat,1:nsig0,igresb:igrese,1:ntemp)
                    read(indev) dum4(1:nrat,1:nsig0,igresb:igrese,1:ntemp)
                    lib%rif(ir)%sct(1:nrat,1:nsig0,igresb:igrese,1:ntemp)=dum4(1:nrat,1:nsig0,igresb:igrese,1:ntemp)
                    deallocate(dum4)
                else
                    lib%rif(ir)%abs(1:nrat,1:nsig0,igresb:igrese,1:ntemp)=ldiso(i1)%rif(ir)%abs(1:nrat,1:nsig0,igresb:igrese,1:ntemp)
                    lib%rif(ir)%sct(1:nrat,1:nsig0,igresb:igrese,1:ntemp)=ldiso(i1)%rif(ir)%sct(1:nrat,1:nsig0,igresb:igrese,1:ntemp)
                endif
                if (ifis.ne.0) then
                    allocate(lib%rif(ir)%fis(1:nrat,1:nsig0,igresb:igrese,1:ntemp))
                    if (is.eq.1) then
                        allocate(dum4(1:nrat,1:nsig0,igresb:igrese,1:ntemp))
                        read(indev) dum4(1:nrat,1:nsig0,igresb:igrese,1:ntemp)
                        lib%rif(ir)%fis(1:nrat,1:nsig0,igresb:igrese,1:ntemp)=dum4(1:nrat,1:nsig0,igresb:igrese,1:ntemp)
                        deallocate(dum4)
                    else
                        lib%rif(ir)%fis(1:nrat,1:nsig0,igresb:igrese,1:ntemp)=ldiso(i1)%rif(ir)%fis(1:nrat,1:nsig0,igresb:igrese,1:ntemp)
                    endif
                endif
            enddo
        enddo
    enddo
    close(indev)

    do i = 1, nelthel
        lib => ldiso(i)
        if (.not.lib%lreso) cycle
        !IF (nTracerCntl%lMLG) THEN
        !    IF (lib%lclad) then
        !        lib%lvabs=lib%lvabs_mlg1G
        !        lib%wgtabs=lib%wgtabs_mlg1G
        !    ENDIF
        !    nid = lib%nid
        !    IF ((nid.eq.47107).or.(nid.eq.47109).or.(nid.eq.49113).or.(nid.eq.49115).or.(nid.eq.48000)) THEN
        !        lib%lvabs=lib%lvabs_mlg1G
        !        lib%wgtabs=lib%wgtabs_mlg1G
        !    ENDIF
        !ENDIF
        allocate(lib%lvabs_log(1:lib%nlv,igresb:igrese))
        allocate(lib%lvflx_log(1:lib%nlvflx,igresb:igrese))
        lib%lvabs_log=dlog(lib%lvabs)
        lib%lvflx_log=dlog(lib%lvflx)
        xxx = mod(nid,1000)
        if (xxx.ge.500) nid = nid - 500
        do j = 1, nelthel
            ljb => ldiso(j)
            if (.not.ljb%lreso) cycle
            if (i.eq.j) cycle
            njd = ljb%nid
            xxx = mod(njd,1000)
            if (xxx.ge.500) njd = njd - 500
            if (nid.eq.njd) cycle
            do ir = 1, lib%nrif
                if (lib%rif(ir)%nid.eq.njd) exit
            enddo
            if (ir.gt.lib%nrif) cycle
            lib%rifid(j)=ir
        enddo
    enddo

end subroutine

subroutine ReadSLB(indev,slbfile,scatord)
    USE ioutil,    ONLY : openfile,  terminate
    USE XSLIB_MOD, ONLY : norghel
    USE SPH_MOD
    implicit none
    integer,intent(in) :: indev,scatord
    character*256,intent(in) :: slbfile
    integer :: i,ig,im,ip,iu,it,ng,nisointer,ninter,nm,np,nu,nt,dg,dt,dp,du,dm,ispf
    integer :: nTEMP,nPLTRAD,nU238ND,nMODXSV,idx_start,dum,ntinteriso
    integer :: nSubring0,nSubring,nSPHreg_srd,SPHreg_srd(0:19)
    integer :: idiso,jj,ii,k,j,xxx,z,nid
    INTEGER :: nSPHreg,nrestfxr,nsphdata_tot,nsphreg_tot
    !real(4) :: ndrat2u238
    real :: ndrat2u238  ! EDIT by JSU ('21 MAR)
    real :: r,dumr(100)
    type SPHINTER_temp
        integer :: idiso,nndrat2u238,idx
        real :: ndrat2u238(100)
    end type
    type(SPHINTER_temp) :: tempINTER(100)
    TYPE(SPHvar_type),POINTER :: svr

    close(indev)
    call openfile(indev,.TRUE.,.TRUE.,.FALSE.,slbfile)
    read(indev) ng,ninter,nt,np,nu,nm,nSPHreg,nSubring0,nSPHreg_srd,nspf
    allocate(SPHvar(0:nspf))
    allocate(SPHdata(0:nspf))
    allocate(DI(0:nspf))
    allocate(idxmap(0:nspf))
    svr=>SPHvar(0)
    nrestfxr=nSPHreg-nSubring0
    if (nrestfxr.gt.2) call terminate('Something is wrong in the SPH library.(nrestfxr)')
    if (nSPHreg_srd.gt.0) then
        read(indev) SPHreg_srd(1:nSPHreg_srd) ! fuel subdevision...
    else
        read(indev)
    endif
    SPHreg_srd(0)=nSPHreg ! number of entire region (fuel + rest)
    SPHreg_srd(1:nSPHreg_srd)=SPHreg_srd(1:nSPHreg_srd)+nrestfxr
    svr%nrestfxr=nrestfxr
    svr%nTEMP=nt
    svr%nPLTRAD=np
    svr%nU238ND=nu
    svr%nMODXSV=nm
    svr%nSubring0=nSubring0
    svr%nSPHreg=nSPHreg
    svr%nSPHreg_srd=nSPHreg_srd
    svr%ninter=ninter
    allocate(svr%SPHreg_srd(0:nSPHreg_srd))
    svr%SPHreg_srd=SPHreg_srd(0:nSPHreg_srd)
    IF (ng.ne.norghel) call terminate('ERROR in SPH library (Wrong # of Groups)')
    allocate(svr%TEMP(nt),svr%PLTRAD(np),svr%U238ND(nu),svr%MODXSV(nm))
    nsphreg_tot=sum(SPHreg_srd(0:nSPHreg_srd))
    nsphdata_tot=(ninter+1)*nt*np*nu*nm*ng*nsphreg_tot
    allocate(idxmap(0)%idx(ng,nm,nu,np,nt,0:ninter))
    k=1
    do ii=0,ninter
        do it=1,nt
            do ip=1,np
                do iu=1,nu
                    do im=1,nm
                        do ig=1,ng
                            idxmap(0)%idx(ig,im,iu,ip,it,ii)=k
                            k=k+nsphreg_tot
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    allocate(SPHdata(0)%sphf(nsphdata_tot))
    do it=1,nt
        read(indev) dumr(1:1)
        svr%TEMP(it)=dumr(1)
    enddo
    read(indev) svr%PLTRAD(1:np)
    read(indev) svr%U238ND(1:nu)
    read(indev) svr%MODXSV(1:nm)
    read(indev) ntinteriso
    do ii=1,ntinteriso
        read(indev) dum,dum
        read(indev) dumr(1:dum)
    enddo
    if (ninter.gt.0) read(indev)
    do it=1,nt
        do ip=1,np
            do iu=1,nu
                do im=1,nm
                    do ig=1,ng
                        idx_start=idxmap(0)%idx(ig,im,iu,ip,it,0)
                        if (scatord.eq.0) then
                            read(indev)
                            read(indev) SPHdata(0)%sphf(idx_start:idx_start+nsphreg_tot-1)
                        else
                            read(indev) SPHdata(0)%sphf(idx_start:idx_start+nsphreg_tot-1)
                            read(indev)
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo
    ! Check # of isotopes
    if (ninter.eq.0) then
        svr%nISOINTER=0
        svr%MAPISOINTER=0
    else
        k=1; j=1
        do ii=1,ninter
            read(indev) idiso,ndrat2u238
            tempINTER(k)%idiso=idiso
            tempINTER(k)%idx=1
            if (k.gt.1) then
                if (idiso.eq.tempINTER(k-1)%idiso) then
                    k=k-1
                    j=j+1
                else
                    tempINTER(k-1)%nndrat2u238=j
                    tempINTER(k)%idx=tempINTER(k-1)%idx+j
                    j=1
                endif
            endif
            tempINTER(k)%ndrat2u238(j)=ndrat2u238
            do it=1,nt
                do ip=1,np
                    do iu=1,nu
                        do im=1,nm
                            do ig=1,ng
                                read(indev)
                                read(indev)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            k=k+1
        enddo
        tempINTER(k-1)%nndrat2u238=j
        nisointer=k-1
        svr%nISOINTER=nisointer
        allocate(svr%ISOINTER(nisointer))
        svr%MAPISOINTER=0
        do ii=1,nisointer
            svr%ISOINTER(ii)%idiso=tempINTER(ii)%idiso
            svr%ISOINTER(ii)%nndrat2u238=tempINTER(ii)%nndrat2u238
            allocate(svr%ISOINTER(ii)%ndrat2u238(svr%ISOINTER(ii)%nndrat2u238))
            do jj=1,svr%ISOINTER(ii)%nndrat2u238
                svr%ISOINTER(ii)%ndrat2u238(jj)=tempINTER(ii)%ndrat2u238(jj)
                svr%ISOINTER(ii)%idx=tempINTER(ii)%idx
            enddo
            svr%MAPISOINTER(tempINTER(ii)%idiso)=ii

            xxx=mod(tempINTER(ii)%idiso,1000)
            z=(tempINTER(ii)%idiso-xxx)/1000
            nid=tempINTER(ii)%idiso
            if (xxx.gt.500) then
                nid = nid - 500
                svr%MAPISOINTER(nid)=ii
            else
                nid = nid + 500
                svr%MAPISOINTER(nid)=ii
            endif
        enddo
        rewind(indev)
        do i=1,nt+ng*nt*np*nu*nm*2+6+1+ntinteriso*2
            read(indev)
        enddo

        do ii=1,ninter
            read(indev)
            do it=1,nt
                do ip=1,np
                    do iu=1,nu
                        do im=1,nm
                            do ig=1,ng
                                idx_start=idxmap(0)%idx(ig,im,iu,ip,it,ii)
                                if (scatord.eq.0) then
                                    read(indev)
                                    read(indev) SPHdata(0)%sphf(idx_start:idx_start+nsphreg_tot-1)
                                else
                                    read(indev) SPHdata(0)%sphf(idx_start:idx_start+nsphreg_tot-1)
                                    read(indev)
                                endif
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
    endif
    allocate(DI(0)%spDI(np))
    do ip=1,np
        allocate(DI(0)%spDI(ip)%DIpos(nsubring0))
        do i=1,nsubring0
            allocate(DI(0)%spDI(ip)%DIpos(i)%rc(i))
            r=svr%PLTRAD(ip)
            call calRingCenter(r,i,DI(0)%spDI(ip)%DIpos(i)%rc)
        enddo
    enddo
    nullify(svr)
    ! special SPH library like AIC
    do ispf=1,nspf
        svr=>SPHvar(ispf)
        read(indev) ng,nt,np,nm,nSPHreg,nSubring0,nSPHreg_srd
        nrestfxr=nSPHreg-nSubring0
        if (nrestfxr.gt.2) call terminate('Something is wrong in the special SPH library.(nrestfxr)')
        if (nSPHreg_srd.gt.0) then
            read(indev) SPHreg_srd(1:nSPHreg_srd)
        else
            read(indev)
        endif
        SPHreg_srd(0)=nSPHreg
        SPHreg_srd(1:nSPHreg_srd)=SPHreg_srd(1:nSPHreg_srd)+nrestfxr
        svr%nrestfxr=nrestfxr
        svr%nTEMP=nt
        svr%nPLTRAD=np
        svr%nMODXSV=nm
        svr%nSubring0=nSubring0
        svr%nSPHreg=nSPHreg
        svr%nSPHreg_srd=nSPHreg_srd
        allocate(svr%SPHreg_srd(0:nSPHreg_srd))
        svr%SPHreg_srd=SPHreg_srd(0:nSPHreg_srd)
        IF (ng.ne.norghel) call terminate('ERROR in special SPH library (Wrong # of Groups)')
        allocate(svr%TEMP(nt),svr%PLTRAD(np),svr%MODXSV(nm))
        nsphreg_tot=sum(SPHreg_srd(0:nSPHreg_srd))
        nsphdata_tot=nt*np*nm*ng*nsphreg_tot
        allocate(idxmap(ispf)%spidx(ng,nm,np,nt))
        k=1
        do it=1,nt
            do ip=1,np
                do im=1,nm
                    do ig=1,ng
                        idxmap(ispf)%spidx(ig,im,ip,it)=k
                        k=k+nsphreg_tot
                    enddo
                enddo
            enddo
        enddo
        allocate(SPHdata(ispf)%sphf(nsphdata_tot))
        do it=1,nt
            read(indev) dumr(1:1)
            svr%TEMP(it)=dumr(1)
        enddo
        read(indev) svr%PLTRAD(1:np)
        read(indev) svr%MODXSV(1:nm)
        do it=1,nt
            do ip=1,np
                do im=1,nm
                    do ig=1,ng
                        idx_start=idxmap(ispf)%spidx(ig,im,ip,it)
                        if (scatord.eq.0) then
                            read(indev)
                            read(indev) SPHdata(ispf)%sphf(idx_start:idx_start+nsphreg_tot-1)
                        else
                            read(indev) SPHdata(ispf)%sphf(idx_start:idx_start+nsphreg_tot-1)
                            read(indev)
                        endif
                    enddo
                enddo
            enddo
        enddo
        do i=1,nsphdata_tot
            if (SPHdata(ispf)%sphf(i).ge.100._4) then
                k=i+1
                do while (SPHdata(ispf)%sphf(k).ge.100._4)
                    k=k+1
                enddo
                SPHdata(ispf)%sphf(i)=SPHdata(ispf)%sphf(k)
            endif
        enddo
        allocate(DI(ispf)%spDI(np))
        do ip=1,np
            allocate(DI(ispf)%spDI(ip)%DIpos(nsubring0))
            do i=1,nsubring0
                allocate(DI(ispf)%spDI(ip)%DIpos(i)%rc(i))
                r=svr%PLTRAD(ip)
                call calRingCenter(r,i,DI(ispf)%spDI(ip)%DIpos(i)%rc)
            enddo
        enddo
        nullify(svr)
    enddo
    close(indev)
    contains

    subroutine calRingCenter(r,n,arr)
        use PARAM, only : pi
        implicit none
        integer,intent(in) :: n
        real,intent(in) :: r
        real,intent(out) :: arr(n)
        integer :: j
        real :: area,b(n)
        area=pi*r**2
        area=area/real(n,8)
        do j=1,n
           b(j)=dsqrt(real(j,8)*area/pi)
        enddo
        arr(1)=b(1)/2
        do j=2,n
            arr(j)=(b(j-1)+b(j))/2
        enddo
    end subroutine

end subroutine

SUBROUTINE ReadXSL(indev,xslfile,ng,nofg,norg,ntiso)   !Binary File, HELIOS format
  USE ioutil,         ONLY : terminate, toupper, openfile
  USE XSLIB_MOD,      ONLY : libdata, ldiso, uhel, enbhel, igresb, igrese, &
                             noghel, nofghel, norghel, notghel, nelthel, nelrhel, nchihel, &
                             nreshel, nburhel, nbahel, nfishel, np1hel, &
                             indn2nhel, indn3nhel, n2nhel, n3nhel, indchxhel, &
                             nchixhel, chixhel, indp2mhel, indp3mhel, &
                             indsghel, nsubghel, nflxhel, nxsghel, nwsghel, nwasgfhel, &
                             nuclidhel, idcodehel, awhel, ntemphel, npothel, ntabhel, &
                             rstab1hel, resdathel, rstab2hel, &
                             nxsthel, xstemphel, istarthel, xsdatahel, &
                             idburhel, idfishel, xenhel, xen0hel, betadnhel, decayhel, yieldhel, &
                             idnp1hel, np1thel, p1temphel, ip1stahel, p1datahel, ip2stahel, p2datahel, ip3stahel, p3datahel, &
                             infchxhel, idn2nhel, xsn2nhel, idn3nhel, xsn3nhel, &
                             nwaddhel, xsghel, wsghel, xasgfhel, wasgfhel, &
                             nlvmax
  USE NUCLIDMAP_MOD,  ONLY : nuclidmap
  USE ALLOCS
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: indev
  CHARACTER,INTENT(IN) :: xslfile*256
  INTEGER,INTENT(OUT) :: ng,nofg,norg,ntiso
  REAL,PARAMETER :: neutronmass = 1.008664916d0

  INTEGER :: i,j,ig,jg,it,ip,ilv
  INTEGER :: iadr,iadr0,iadr1,iadr2,iadr3,itadr,ib,ie
  INTEGER :: idum,ihel(10),mtot
  INTEGER :: ntemp,np1temp
  REAL(4) :: fdum,sigb
  REAL(8) :: nu
  CHARACTER*120 :: SSW(3)
  TYPE(libdata),POINTER :: lib

  CLOSE(indev)
  CALL openfile(indev,.TRUE.,.TRUE.,.FALSE.,xslfile)

  !Library Information
  READ(indev) SSW(1)
  READ(indev) SSW(2)
  READ(indev) SSW(3)

  !Numbers
  READ(indev) ihel(1:10)
  READ(indev) (idum, i=1,3), noghel, nofghel, norghel, notghel, nchihel, nelthel, nelrhel
  READ(indev) nreshel, nburhel, nbahel, nfishel, np1hel, indn2nhel, indn3nhel, n2nhel, n3nhel, indchxhel
  READ(indev) nchixhel, indp2mhel, indp3mhel, indsghel
  READ(indev) nsubghel, nflxhel, nxsghel, nwsghel, nwasgfhel

  ng = noghel
  nofg = nofghel
  norg = norghel
  ntiso = nelthel
  igresb = nofghel + 1
  igrese = nofghel + norghel

  ALLOCATE(ldiso(ntiso))
  !Init.
  DO i = 1, nelthel
    ldiso(i)%ityp = 0
    ldiso(i)%ichi = 0
    ldiso(i)%ifis = 0
    ldiso(i)%ibur = 0
    ldiso(i)%inmn = 0
    ldiso(i)%ntemp = 0
    ldiso(i)%np1temp = 0
    ldiso(i)%nrtemp = 0

    ldiso(i)%kappa = 0.
    ldiso(i)%kappa0 = 0.
  END DO ! of i

  !Group Boundary
  CALL dmalloc(enbhel,ng+1)
  CALL dmalloc(uhel,ng+1)
  READ(indev) uhel(1:ng)
  DO i = 1, ng
    enbhel(i) = 1.e7/exp(uhel(i))
  END DO ! of i
  enbhel(ng+1) = 1.0e-4

  !Global Chi - IGNORE
  READ(indev) (fdum, i=1,nchihel)

  !NID
  ALLOCATE(nuclidhel(nelthel))
  READ(indev) nuclidhel(1:nelthel)
  !ITYP
  ALLOCATE(idcodehel(nelthel))
  READ(indev) idcodehel(1:nelthel)
  !Atomic Weight
  ALLOCATE(awhel(nelthel))
  READ(indev) awhel(1:nelthel)
  DO i = 1, nelthel
    ldiso(i)%nid = nuclidhel(i)
    ldiso(i)%ityp = idcodehel(i)
    ldiso(i)%aw = dble(awhel(i))
  END DO ! of i
  DEALLOCATE(nuclidhel,idcodehel,awhel)

  !Resonance Data for Resonance Nuclides
  ALLOCATE(ntemphel(nreshel),npothel(nreshel),ntabhel(nreshel))
  READ(indev) ntemphel(1:nreshel) !Number of Resonance Temp.
  READ(indev) npothel(1:nreshel)  !Number of Sigma Potential
  READ(indev) ntabhel(1:nreshel)  !Number of

  !General Resonance Data
  ALLOCATE(rstab1hel(ihel(2)))
  ALLOCATE(resdathel(norghel,ihel(3)))
  ALLOCATE(rstab2hel(norghel,ihel(4)))
  !Square roots or resonance temperature and potential XSs
  READ(indev) rstab1hel(1:ihel(2))
  !Resonance Data such as Lamda*SigP, Resonance Scattering XS
  DO ig = 1, norghel
    READ(indev) resdathel(ig,1:ihel(3))
    READ(indev) rstab2hel(ig,1:ihel(4))
  END DO ! of ig
  DO i = 1, ntiso
    CALL dmalloc0(ldiso(i)%lamsigp,igresb,igrese)
    IF(i > nelrhel) EXIT
    do j = 1, norghel
      ig = j + nofghel
      ldiso(i)%lamsigp(ig) = resdathel(j,4*i-3)
    END DO ! of j
    ldiso(i)%sigp = 0.d0
  END DO ! of i

  !Resonance Data for Resonance Nuclides
  DO i = 1, nreshel
    !Resonance Temp.
    ldiso(i)%nrtemp = mod(ntemphel(i),100)
    CALL dmalloc(ldiso(i)%rtemp,ldiso(i)%nrtemp)
    CALL dmalloc(ldiso(i)%rtempsq,ldiso(i)%nrtemp)
    iadr = ntemphel(i)/100 - 1
    DO j = 1, ldiso(i)%nrtemp
      ldiso(i)%rtempsq(j) = rstab1hel(iadr+j)
      ldiso(i)%rtemp(j) = ldiso(i)%rtempsq(j)**2
    END DO ! of j
    !Background XS
    ldiso(i)%nsig0 = mod(npothel(i),100)
    CALL dmalloc(ldiso(i)%sig0sq,ldiso(i)%nsig0)
    CALL dmalloc0(ldiso(i)%ri_a,1,ldiso(i)%nsig0,igresb,igrese,1,ldiso(i)%nrtemp)
    CALL dmalloc0(ldiso(i)%xsalog,1,ldiso(i)%nsig0,igresb,igrese,1,ldiso(i)%nrtemp)
    iadr = npothel(i)/100 - 1
    DO j = 1, ldiso(i)%nsig0
      ldiso(i)%sig0sq(j) = rstab1hel(iadr+j)
    END DO ! of j
  END DO ! of i
  !Resonance Integral and Absorption XS
  DO j = 1, norghel
    iadr = 0
    ig = j + nofghel
    DO i = 1, nreshel
      DO it = 1, ldiso(i)%nrtemp
        DO ip = 1, ldiso(i)%nsig0
          iadr = iadr+1
          ldiso(i)%ri_a(ip,ig,it) = rstab2hel(j,iadr)
          sigb = ldiso(i)%sig0sq(ip)**2 + ldiso(i)%lamsigp(ig)
          ldiso(i)%xsalog(ip,ig,it) = ldiso(i)%ri_a(ip,ig,it)*sigb/(sigb-ldiso(i)%ri_a(ip,ig,it))
          ldiso(i)%xsalog(ip,ig,it) = log(ldiso(i)%xsalog(ip,ig,it))
        END DO ! of ip
      END DO ! of it
      !IGNORE - Resonance Integral for Fission
      if(ldiso(i)%ityp == 3) iadr = iadr + ldiso(i)%nrtemp*ldiso(i)%nsig0
    END DO ! of i
  END DO ! of j
  DEALLOCATE(ntemphel,npothel,ntabhel)
  DEALLOCATE(rstab1hel,resdathel,rstab2hel)

  !Temperature info. for Cross Section
  ALLOCATE(nxsthel(nelthel))
  ALLOCATE(xstemphel(ihel(5)))
  READ(indev) nxsthel(1:nelthel)
  READ(indev) xstemphel(1:ihel(5))

  !Cross Section Data
  ALLOCATE(istarthel(ng,3*ihel(5)))
  ALLOCATE(xsdatahel(ng,ihel(7)))
  DO ig = 1, ng
    READ(indev) istarthel(ig,1:3*ihel(5))
    mtot = istarthel(ig,3*ihel(5)-2) + istarthel(ig,3*ihel(5)) -1
    READ(indev) xsdatahel(ig,1:mtot)
  END DO ! of ig

  !Beta and Fission Index
  ALLOCATE(idburhel(nburhel),idfishel(nfishel))
  ALLOCATE(xenhel(nfishel),xen0hel(nfishel),betadnhel(6,nfishel))
  ALLOCATE(decayhel(nburhel),yieldhel(nfishel,nburhel))

  READ(indev) idburhel(1:nburhel)
  READ(indev) idfishel(1:nfishel)
  READ(indev) xenhel(1:nfishel)
  READ(indev) xen0hel(1:nfishel)
  READ(indev) betadnhel(1:6,1:nfishel)
  READ(indev) decayhel(1:nburhel)
  READ(indev) yieldhel(1:nfishel,1:nburhel)

  !Rearrange Fission Data
  DO j = 1, nfishel
    i = idfishel(j)
    ldiso(i)%ifis = 1
    ldiso(i)%ichi = 1
    CALL dmalloc0(ldiso(i)%beta,0,6)
    ldiso(i)%kappa = xenhel(j)
    ldiso(i)%kappa0 = xen0hel(j)
    ldiso(i)%beta(1:6) = betadnhel(:,j)
    ldiso(i)%beta(0) = sum(betadnhel(:,j))
  END DO ! of j
  DEALLOCATE(xenhel,xen0hel,betadnhel)

  !Rearrange Decay Data - IGNORE
  DEALLOCATE(idburhel,idfishel)
  DEALLOCATE(decayhel,yieldhel)

  !Rearrange Cross Section Data
  iadr = 0
  DO i = 1, nelthel
    ntemp = mod(nxsthel(i),100)
    itadr = nxsthel(i)/100
    !Temperature for Cross Section Data
    ldiso(i)%ntemp = ntemp
    CALL dmalloc0(ldiso(i)%temp,1,ntemp)
    ldiso(i)%temp(1:ntemp) = xstemphel(itadr:itadr+ntemp-1)
    CALL dmalloc0(ldiso(i)%siga,1,ng,1,ntemp)
    CALL dmalloc0(ldiso(i)%sigs,1,ng,1,ntemp)
    CALL dmalloc0(ldiso(i)%sigtr,1,ng,1,ntemp)
    CALL dmalloc0(ldiso(i)%sigstr,1,ng,1,ntemp)
    CALL dmalloc0(ldiso(i)%sigss,1,ng,1,ntemp)
    IF(ldiso(i)%ifis >= 1) THEN
      CALL dmalloc0(ldiso(i)%sigf,1,ng,1,ntemp)
      CALL dmalloc0(ldiso(i)%signf,1,ng,1,ntemp)
    END IF
    ALLOCATE(ldiso(i)%sm(ng,ntemp))

    DO it = 1, ntemp
      iadr1 = iadr + 1
      iadr2 = iadr + 2
      iadr3 = iadr + 3
      DO ig = 1, ng
        j = istarthel(ig,iadr1)
        ldiso(i)%siga(ig,it) = xsdatahel(ig,j)
        IF(i <= nelrhel) THEN
          IF(ldiso(i)%ifis == 1) THEN
            ldiso(i)%sigf(ig,it) = xsdatahel(ig,j+1)
            ldiso(i)%signf(ig,it) = xsdatahel(ig,j+2)
          END IF
          ldiso(i)%sigtr(ig,it) = xsdatahel(ig,j+3)
          ldiso(i)%sigs(ig,it) = xsdatahel(ig,j+4)
          ib = ig + 6 - istarthel(ig,iadr2)
          ie = ig + istarthel(ig,iadr3) - istarthel(ig,iadr2)
          ldiso(i)%sm(ig,it)%ib = ib
          ldiso(i)%sm(ig,it)%ie = ie
          CALL dmalloc0(ldiso(i)%sm(ig,it)%from,ib,ie)
          ldiso(i)%sm(ig,it)%from(ib:ie) = xsdatahel(ig,j+5:j+5+(ie-ib))
        ELSE
          ldiso(i)%sm(ig,it)%ib = ig
          ldiso(i)%sm(ig,it)%ie = ig
          CALL dmalloc0(ldiso(i)%sm(ig,it)%from,ig,ig)
          ldiso(i)%sm(ig,it)%from(ig) = 0.
        END IF
      END DO ! of ig
      iadr = iadr + 3
      !Scattering Matrix Post-Processing
      DO jg = 1, ng
        DO ig = 1, ng
          IF(jg >= ldiso(i)%sm(ig,it)%ib .AND. jg <= ldiso(i)%sm(ig,it)%ie) THEN
            ldiso(i)%sm(jg,it)%ioutsb = ig
            exit
          END IF
        END DO ! of ig
        DO ig = ng, 1, -1
          IF(jg >= ldiso(i)%sm(ig,it)%ib .AND. jg <= ldiso(i)%sm(ig,it)%ie) THEN
            ldiso(i)%sm(jg,it)%ioutse = ig
            exit
          END IF
        END DO ! of ig
        DO ig = ldiso(i)%sm(jg,it)%ib, ldiso(i)%sm(jg,it)%ie
          ldiso(i)%sigstr(ig,it) = ldiso(i)%sigstr(ig,it) + ldiso(i)%sm(jg,it)%from(ig)
        END DO ! of ig
      END DO ! of jg
      DO ig = 1, ng
        ldiso(i)%sigss(ig,it) = ldiso(i)%sm(ig,it)%from(ig) + (ldiso(i)%sigs(ig,it) - ldiso(i)%sigstr(ig,it))
      END DO ! of ig
    END DO ! of it
  END DO ! of i
  DEALLOCATE(nxsthel,xstemphel)
  DEALLOCATE(istarthel,xsdatahel)

  !Pn Scattering Data
  ALLOCATE(idnp1hel(np1hel),np1thel(np1hel),p1temphel(ihel(8)))
  READ(indev) idnp1hel(1:np1hel)
  READ(indev) np1thel(1:np1hel)
  READ(indev) p1temphel(1:ihel(8))

  !P1 Scattering Cross Section
  ALLOCATE(ip1stahel(ng,3*ihel(8)))
  ALLOCATE(p1datahel(ng,ihel(9)))
  DO ig = 1, ng
    READ(indev) ip1stahel(ig,1:3*ihel(8))
    mtot = ip1stahel(ig,3*ihel(8)-2) + ip1stahel(ig,3*ihel(8)) -1
    READ(indev) p1datahel(ig,1:mtot)
  END DO ! of ig

  !P2 Scattering Cross Section
  ALLOCATE(ip2stahel(ng,3*ihel(8)))
  ALLOCATE(p2datahel(ng,ihel(9)))
  DO ig = 1, ng
    READ(indev) ip2stahel(ig,1:3*ihel(8))
    mtot = ip2stahel(ig,3*ihel(8)-2) + ip2stahel(ig,3*ihel(8)) -1
    READ(indev) p2datahel(ig,1:mtot)
  END DO ! of ig

  !P3 Scattering Cross Section
  ALLOCATE(ip3stahel(ng,3*ihel(8)))
  ALLOCATE(p3datahel(ng,ihel(9)))
  DO ig = 1, ng
    READ(indev) ip3stahel(ig,1:3*ihel(8))
    mtot = ip3stahel(ig,3*ihel(8)-2) + ip3stahel(ig,3*ihel(8)) -1
    READ(indev) p3datahel(ig,1:mtot)
  END DO ! of ig

  !Rearrange Pn Scattering Data
  iadr = 0
  DO j = 1, np1hel
    i = idnp1hel(j)
    np1temp = mod(np1thel(j),100)
    itadr = np1thel(j)/100
    ldiso(i)%np1temp = np1temp
    CALL dmalloc0(ldiso(i)%p1temp,1,np1temp)
    ldiso(i)%p1temp = p1temphel(itadr:itadr+np1temp-1)
    ALLOCATE(ldiso(i)%smp1(ng,np1temp))
    ALLOCATE(ldiso(i)%smp2(ng,np1temp))
    ALLOCATE(ldiso(i)%smp3(ng,np1temp))
    CALL dmalloc0(ldiso(i)%sigsp1,1,ng,1,np1temp)
    CALL dmalloc0(ldiso(i)%sigsp2,1,ng,1,np1temp)
    CALL dmalloc0(ldiso(i)%sigsp3,1,ng,1,np1temp)

    DO it = 1, np1temp
      iadr1 = iadr + 1
      iadr2 = iadr + 2
      iadr3 = iadr + 3
      DO ig = 1, ng
        !P1
        ib = ig + 2 - ip1stahel(ig,iadr2)
        ie = ig + ip1stahel(ig,iadr3) - ip1stahel(ig,iadr2)
        ldiso(i)%smp1(ig,it)%ib = ib
        ldiso(i)%smp1(ig,it)%ie = ie
        ldiso(i)%sigsp1(ig,it) = p1datahel(ig,ip1stahel(ig,iadr1))
        CALL dmalloc0(ldiso(i)%smp1(ig,it)%from,ib,ie)
        ldiso(i)%smp1(ig,it)%from(ib:ie) = p1datahel(ig,ip1stahel(ig,iadr1)+1:ip1stahel(ig,iadr1)+1+ie-ib)
        !P2
        ib = ig + 2 - ip2stahel(ig,iadr2)
        ie = ig + ip2stahel(ig,iadr3) - ip2stahel(ig,iadr2)
        ldiso(i)%smp2(ig,it)%ib = ib
        ldiso(i)%smp2(ig,it)%ie = ie
        ldiso(i)%sigsp2(ig,it) = p2datahel(ig,ip2stahel(ig,iadr1))
        CALL dmalloc0(ldiso(i)%smp2(ig,it)%from,ib,ie)
        ldiso(i)%smp2(ig,it)%from(ib:ie) = p2datahel(ig,ip2stahel(ig,iadr1)+1:ip2stahel(ig,iadr1)+1+ie-ib)
        !P3
        ib = ig + 2 - ip3stahel(ig,iadr2)
        ie = ig + ip3stahel(ig,iadr3) - ip3stahel(ig,iadr2)
        ldiso(i)%smp3(ig,it)%ib = ib
        ldiso(i)%smp3(ig,it)%ie = ie
        ldiso(i)%sigsp3(ig,it) = p3datahel(ig,ip3stahel(ig,iadr1))
        CALL dmalloc0(ldiso(i)%smp3(ig,it)%from,ib,ie)
        ldiso(i)%smp3(ig,it)%from(ib:ie) = p3datahel(ig,ip3stahel(ig,iadr1)+1:ip3stahel(ig,iadr1)+1+ie-ib)
      END DO ! of ig
      iadr = iadr3
    END DO ! of it
  END DO ! of j
  DO i = 1, nelthel
    IF(ldiso(i)%np1temp > 0) CYCLE
    CALL dmalloc0(ldiso(i)%sigsp1,1,ng,1,ntemp)
  END DO ! of i
  DEALLOCATE(idnp1hel,np1thel,p1temphel)
  DEALLOCATE(ip1stahel,p1datahel)
  DEALLOCATE(ip2stahel,p2datahel)
  DEALLOCATE(ip3stahel,p3datahel)

  !Isotopewise Fission Spectrum
  ALLOCATE(infchxhel(3,nchixhel))
  READ(indev) infchxhel(1:3,1:nchixhel)
  j = infchxhel(2,1) - 1
  DO i = 1, nchixhel
    infchxhel(2,i) = infchxhel(2,i) - j
    infchxhel(3,i) = infchxhel(3,i) - j
  END DO ! of i
  ALLOCATE(chixhel(nchixhel,nchihel))
  READ(indev) ((chixhel(i,j), j=1,infchxhel(3,i)-infchxhel(2,i)+1),i=1,nchixhel)
  DO j = 1, nchixhel
    i = infchxhel(1,j)
    CALL dmalloc0(ldiso(i)%chi,1,ng)
    fdum = 1./sum(chixhel(j,1:nchihel))
    DO ig = 1, nchihel
      ldiso(i)%chi(ig) = chixhel(j,ig)*fdum
    END DO ! of ig
  END DO ! of j
  DEALLOCATE(infchxhel)
  DEALLOCATE(chixhel)

  !n2n, n3n Cross Section
  ALLOCATE(idn2nhel(n2nhel),xsn2nhel(n2nhel,ng))
  ALLOCATE(idn3nhel(n3nhel),xsn3nhel(n3nhel,ng))
  READ(indev) idn2nhel(1:n2nhel)
  READ(indev) xsn2nhel(1:n2nhel,1:ng)
  READ(indev) idn3nhel(1:n3nhel)
  READ(indev) xsn3nhel(1:n3nhel,1:ng)

  !Rearrange n2n, n3n Data
  DO j = 1, n2nhel
    i = idn2nhel(j)
    ldiso(i)%inmn = ldiso(i)%inmn + 1
    CALL dmalloc0(ldiso(i)%sign2n,1,ng)
    ldiso(i)%sign2n = xsn2nhel(j,:)
  END DO ! of j
  DO j = 1, n3nhel
    i = idn3nhel(j)
    ldiso(i)%inmn = ldiso(i)%inmn + 2
    CALL dmalloc0(ldiso(i)%sign3n,1,ng)
    ldiso(i)%sign3n = xsn3nhel(j,:)
  END DO ! of j
  DEALLOCATE(idn2nhel,xsn2nhel)
  DEALLOCATE(idn3nhel,xsn3nhel)

  !Subgroup Data
  ALLOCATE(nwaddhel(1:3*nreshel))
  ALLOCATE(xsghel(norg,nxsghel))
  ALLOCATE(wsghel(norg,nwsghel))
  ALLOCATE(xasgfhel(ng,4*nreshel))
  ALLOCATE(wasgfhel(ng,nwasgfhel))
  READ(indev) nwaddhel(1:3*nreshel)
  DO ig = 1, norg
    READ(indev) xsghel(ig,1:nxsghel)
  END DO ! of ig
  DO ig = 1, norg
    READ(indev) wsghel(ig,1:nwsghel)
  END DO ! of ig
  DO ig = 1, norg+1
    READ(indev) xasgfhel(ig,1:4*nreshel)
  END DO ! of ig
  DO ig = 1, norg+1
    READ(indev) wasgfhel(ig,1:nwasgfhel)
  END DO ! of ig

  !Rearrange Subgroup Data
  DO i = 1, nreshel
    ntemp = ldiso(i)%nrtemp
    ldiso(i)%lreso = .true.
    ldiso(i)%lfuel = .true.
    IF(indsghel /= 0) THEN
      CALL dmalloc0(ldiso(i)%lvabs,1,nsubghel,igresb,igrese)
      CALL dmalloc0(ldiso(i)%lvabs_log,1,nsubghel,igresb,igrese)
      CALL dmalloc0(ldiso(i)%wgtabs,1,nsubghel,1,ntemp,igresb,igrese)
      CALL dmalloc0(ldiso(i)%lvflx,1,nflxhel,igresb,igrese)
      IF(ldiso(i)%ityp > 0) THEN
        CALL dmalloc0(ldiso(i)%lvfis,1,nsubghel,igresb,igrese)
        CALL dmalloc0(ldiso(i)%wgtfis,1,nsubghel,1,ntemp,igresb,igrese)
      END IF
      !Subgroup Level for FSP Calc.
      ldiso(i)%nlvflx = nflxhel
      DO jg = 1, norg
        ig = jg + igresb - 1
        ldiso(i)%lvflx(:,ig) = xasgfhel(jg,(i-1)*nflxhel+1:i*nflxhel)
      END DO ! of jg
      !Subgroup Level and Weight
      DO jg = 1, norg
        iadr = nwaddhel(i)
        iadr1 = nwaddhel(i+nreshel)

        ig = jg + igresb - 1
        ldiso(i)%lvabs(:,ig) = xsghel(jg,iadr:iadr+nsubghel-1)
        ldiso(i)%lvabs_log(:,ig) = dlog(ldiso(i)%lvabs(:,ig))
        iadr = iadr + nsubghel
        IF(ldiso(i)%ityp == 3) THEN
          ldiso(i)%lvfis(:,ig) = xsghel(jg,iadr:iadr+nsubghel-1)
          iadr = iadr + nsubghel
          !divide by nu
          nu = 0.d0
          do it=1, ldiso(i)%ntemp
            nu = nu + ldiso(i)%signf(ig,it)/ldiso(i)%sigf(ig,it)
          end do ! of it
          nu = nu	/ ldiso(i)%ntemp
          ldiso(i)%lvfis(:,ig) = ldiso(i)%lvfis(:,ig) / nu
        END IF

        !Subgroup Weight
        DO it = 1, ntemp
            ldiso(i)%wgtabs(:,it,ig) = wsghel(jg,iadr1:iadr1+nsubghel-1)
            iadr1 = iadr1 + nsubghel
            IF(ldiso(i)%ityp == 3) THEN
              ldiso(i)%wgtfis(:,it,ig) = wsghel(jg,iadr1:iadr1+nsubghel-1)
              iadr1 = iadr1 + nsubghel
            END IF
        END DO ! of it
      END DO ! of jg
      ldiso(i)%nlv = nsubghel

      CALL dmalloc0(ldiso(i)%lvabs_log,1,nsubghel,igresb,igrese)
      CALL dmalloc0(ldiso(i)%lvflx_log,1,nflxhel,igresb,igrese)
      ldiso(i)%lvabs_log = dlog(ldiso(i)%lvabs)
      ldiso(i)%lvflx_log = dlog(ldiso(i)%lvflx)
    END IF
  END DO ! of i
  nlvmax = nsubghel

  DEALLOCATE(nwaddhel,xsghel)
  DEALLOCATE(wsghel,xasgfhel,wasgfhel)

  CLOSE(indev)

  CALL nuclidmap
  RETURN
END SUBROUTINE ReadXSL  !Binary File, HELIOS format

SUBROUTINE ReadPLC(indev,xslfile,ng,nofg,norg,ntiso)
  USE ioutil,         ONLY : terminate, toupper, openfile
  USE XSLIB_MOD,      ONLY : libdata,ldiso,uhel,enbhel,igresb,igrese,&
														 noghel, nofghel, norghel, notghel, nchihel, nelthel, nelrhel,&
                             nreshel, nburhel, nfishel, nsubghel, nflxhel, &
                             nlvmax
  USE NUCLIDMAP_MOD,  ONLY : nuclidmap
  USE ALLOCS
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: indev
  CHARACTER*(*),INTENT(in) :: xslfile
  INTEGER,INTENT(OUT) :: ng,nofg,norg,ntiso
  REAL,PARAMETER :: neutronmass = 1.008664916d0

  INTEGER :: i,j,ig,jg,it,ip,id,igx,itx
  INTEGER :: ix,nid,aw,ityp,ifis,ibur,ichi,inmn,ntemp,nrtemp,np1temp,npot
  INTEGER :: ib,ie
  REAL(4) :: fdum,farr1(5),farr2(500),sigb,nu
  CHARACTER :: aid*20, blockname*20
  CHARACTER :: oneline*300,onelinel*1000,probe,probel
  TYPE(libdata),POINTER :: lib
  EQUIVALENCE(PROBE,ONELINE)
  EQUIVALENCE(PROBEL,ONELINEL)

  CLOSE(indev)
  CALL openfile(indev,.TRUE.,.FALSE.,.FALSE.,xslfile)

  !Read Basic Information
  DO WHILE(.TRUE.)
    READ(indev,'(a1000)', END = 100) onelinel
    onelinel = trim(adjustl(onelinel))
    IF(probel == '/') exit
    IF(probel == ' ' .or. probel == '' .or. onelinel == ' ' .or. onelinel == '\t') CYCLE
    READ(onelinel,*) blockname
    CALL toupper(blockname)

    SELECT CASE(blockname)
      CASE('!DIMENSION')
        READ(indev,*) noghel, nofghel, norghel, notghel, nchihel, nelthel, nelrhel,&
                      nreshel, nburhel, nfishel, nsubghel, nflxhel
        igresb = nofghel + 1
        igrese = nofghel + norghel

        ng = noghel
        nofg = nofghel
        norg = norghel
        ntiso = nelthel
      CASE('!GROUP')
        CALL dmalloc(enbhel,noghel+1)
        CALL dmalloc(uhel,noghel)
        READ(indev,'(10E14.7)') uhel(1:noghel)
        DO ig	= 1, ng
          enbhel(ig) = 1.E7/exp(uhel(ig))
        END DO ! of ig
        enbhel(ng+1) = 1.0E-4
      CASE('!CHI')
        READ(indev,*) (fdum, i=1,nchihel)
      CASE('!DIR')
        ALLOCATE(ldiso(nelthel))
        DO i = 1, nelthel
          READ(indev,*) ix,nid,aw,ityp,ifis,ibur,ichi,inmn,ntemp,nrtemp,np1temp,npot,aid
          IF(ix /= i) CALL terminate('Check index of !DIR')
          ldiso(i)%nid = nid
          ldiso(i)%aw = aw
          ldiso(i)%ityp = ityp
          ldiso(i)%ifis = ifis
          ldiso(i)%ibur = ibur
          ldiso(i)%ichi = ichi
          ldiso(i)%inmn = inmn
          ldiso(i)%ntemp = ntemp
          ldiso(i)%nrtemp = nrtemp
          ldiso(i)%np1temp = np1temp
          ldiso(i)%nsig0 = npot
          ldiso(i)%aid = aid
        END DO ! of i
      CASE('!NUC:')
        backspace(indev)
        exit
      CASE DEFAULT
        CALL terminate(' Block Name '//trim(blockname)//' NOT Allowed...')
    END SELECT
  END DO
100 CONTINUE

  !Allocate
  DO i = 1, ntiso
    ntemp = ldiso(i)%ntemp
    nrtemp = ldiso(i)%nrtemp
    np1temp = ldiso(i)%np1temp
    npot = ldiso(i)%nsig0
    !Temperature
    CALL dmalloc0(ldiso(i)%temp,1,ntemp)
    !Cross Section Data
    CALL dmalloc0(ldiso(i)%siga,1,ng,1,ntemp)
    CALL dmalloc0(ldiso(i)%sigs,1,ng,1,ntemp)
    CALL dmalloc0(ldiso(i)%sigtr,1,ng,1,ntemp)
    CALL dmalloc0(ldiso(i)%sigstr,1,ng,1,ntemp)
    CALL dmalloc0(ldiso(i)%sigss,1,ng,1,ntemp)
    IF(ldiso(i)%ifis >= 1) THEN
      CALL dmalloc0(ldiso(i)%sigf,1,ng,1,ntemp)
      CALL dmalloc0(ldiso(i)%signf,1,ng,1,ntemp)
    END IF
    !Scattering Matrix
    ALLOCATE(ldiso(i)%sm(ng,ntemp))
    !Resonance Data
    CALL dmalloc0(ldiso(i)%lamsigp,igresb,igrese)
    IF(i <= nreshel) THEN
      CALL dmalloc(ldiso(i)%rtemp,nrtemp)
      CALL dmalloc(ldiso(i)%rtempsq,nrtemp)
      CALL dmalloc(ldiso(i)%sig0sq,npot)
      CALL dmalloc0(ldiso(i)%ri_a,1,npot,igresb,igrese,1,nrtemp)
      CALL dmalloc0(ldiso(i)%xsalog,1,npot,igresb,igrese,1,nrtemp)
      CALL dmalloc0(ldiso(i)%lvabs,1,nsubghel,igresb,igrese)
      CALL dmalloc0(ldiso(i)%lvabs_log,1,nsubghel,igresb,igrese)
      CALL dmalloc0(ldiso(i)%wgtabs,1,nsubghel,1,nrtemp,igresb,igrese)
      CALL dmalloc0(ldiso(i)%lvflx,1,nflxhel,igresb,igrese)
      CALL dmalloc0(ldiso(i)%lvflx_log,1,nflxhel,igresb,igrese)
      IF(ldiso(i)%ityp > 0) THEN
        CALL dmalloc0(ldiso(i)%lvfis,1,nsubghel,igresb,igrese)
        CALL dmalloc0(ldiso(i)%wgtfis,1,nsubghel,1,nrtemp,igresb,igrese)
      END IF
    END IF
    !Fission Spectrum
    IF(ldiso(i)%ichi > 0) CALL dmalloc(ldiso(i)%chi,ng)
    !Beta
    IF(ldiso(i)%ifis > 0) CALL dmalloc0(ldiso(i)%beta,0,6)
    !(n,2n), (n,3n) Scattering
    IF(ldiso(i)%inmn == 1 .or. ldiso(i)%inmn == 3) CALL dmalloc0(ldiso(i)%sign2n,1,ng)
    IF(ldiso(i)%inmn == 2 .or. ldiso(i)%inmn == 3) CALL dmalloc0(ldiso(i)%sign3n,1,ng)
    !Pn Scattering
    IF(np1temp > 0) THEN
      CALL dmalloc0(ldiso(i)%p1temp,1,np1temp)
      CALL dmalloc0(ldiso(i)%sigsp1,1,noghel,1,np1temp)
      CALL dmalloc0(ldiso(i)%sigsp2,1,noghel,1,np1temp)
      CALL dmalloc0(ldiso(i)%sigsp3,1,noghel,1,np1temp)
      ALLOCATE(ldiso(i)%smp1(noghel,np1temp))
      ALLOCATE(ldiso(i)%smp2(noghel,np1temp))
      ALLOCATE(ldiso(i)%smp3(noghel,np1temp))
    END IF
  END DO ! of i

  DO WHILE(.TRUE.)
    READ(indev,'(a1000)', END = 200) onelinel
    onelinel = trim(adjustl(onelinel))
    IF(probel == '/') EXIT
    IF(probel == ' ' .or. probel == '' .or. onelinel == ' ' .or. onelinel == '\t') CYCLE
    READ(onelinel,*) blockname
    CALL toupper(blockname)

    READ(indev,'(a256)') oneline
    oneline = trim(adjustl(oneline))
    backspace(indev)
    IF((oneline == ' ' .or. probe == '!') .and. blockname /= '!NUC:') CYCLE

    SELECT CASE(blockname)
      CASE('!NUC:')
        READ(onelinel,*) aid,id,nid,aw,ityp,ifis,ibur,ichi,inmn,ntemp,nrtemp,np1temp,npot,aid
      CASE('!TP1+')
        READ(indev,*) ldiso(id)%temp(1:ntemp)
      CASE('!XSD+')
        DO ig = 1, ng
          DO it = 1, ntemp
            READ(indev,*) igx,itx,farr1(1:5),ib,ie,farr2(ib:ie)
            ldiso(id)%siga(ig,it)  = farr1(1)
            IF(ldiso(id)%ifis >= 1) THEN
              ldiso(id)%sigf(ig,it)  = farr1(2)
              ldiso(id)%signf(ig,it) = farr1(3)
            END IF
            ldiso(id)%sigtr(ig,it) = farr1(4)
            ldiso(id)%sigs(ig,it)  = farr1(5)
            ldiso(id)%sm(ig,it)%ib = ib
            ldiso(id)%sm(ig,it)%ie = ie
            ALLOCATE(ldiso(id)%sm(ig,it)%from(ib:ie))
            ldiso(id)%sm(ig,it)%from(ib:ie) = farr2(ib:ie)
          END DO ! of it
        END DO ! of ig
      CASE('!ABS+')
        DO ig = 1, ng
          DO it = 1, ntemp
            READ(indev,*) igx,itx,ldiso(id)%siga(ig,it)
          END DO ! of it
        END DO ! of ig
      CASE('!POT+')
        READ(indev,*) ldiso(id)%lamsigp(igresb:igrese)
        READ(indev,*) farr2(igresb:igrese)
        READ(indev,*) farr2(igresb:igrese)
      CASE('!TP2+')
        READ(indev,*) ldiso(id)%rtempsq(1:nrtemp)
        DO it = 1, nrtemp
          ldiso(id)%rtemp(it) = ldiso(id)%rtemp(it)**2
        END DO ! of it
      CASE('!RGP+')
        stop '!RPG+ CARD in ReadPLC'
      CASE('!XS0+')
        READ(indev,*) ldiso(id)%sig0sq(1:npot)
      CASE('!RIA+')
        DO ig = igresb, igrese
          DO it = 1, nrtemp
            READ(indev,*) igx,itx,ldiso(id)%ri_a(1:npot,ig,it)
            DO ip = 1, npot
              sigb = ldiso(id)%sig0sq(ip)**2 + ldiso(id)%lamsigp(ig)
              ldiso(id)%xsalog(ip,ig,it) = ldiso(id)%ri_a(ip,ig,it)*sigb/(sigb-ldiso(id)%ri_a(ip,ig,it))
              ldiso(id)%xsalog(ip,ig,it) = log(ldiso(id)%xsalog(ip,ig,it))
            END DO ! of ip
          END DO ! it
        END DO ! of ig
      CASE('!RNF+')
        DO ig = igresb, igrese
          DO it = 1, nrtemp
            READ(indev,*) igx,itx,farr2(1:npot)
          END DO ! it
        END DO ! of ig
      CASE('!SA7S+')
        ldiso(id)%nlv = nsubghel
        ldiso(id)%lreso = .true.
        ldiso(id)%lfuel = .true.
        DO ig = igresb, igrese
          READ(indev,*) igx,ldiso(id)%lvabs(1:nsubghel,ig)
        END DO ! of ig
        ldiso(id)%lvabs_log = dlog(ldiso(id)%lvabs)
      CASE('!SA7W+')
        DO ig = igresb, igrese
          DO it = 1, nrtemp
            READ(indev,*) igx,itx,ldiso(id)%wgtabs(1:nsubghel,it,ig)
          END DO ! of it
        END DO ! of ig
      CASE('!SF7S+')
        DO ig = igresb, igrese
          READ(indev,*) igx,ldiso(id)%lvfis(1:nsubghel,ig)
          !divide by nu
          nu = 0.d0
          do it=1, ldiso(id)%ntemp
            nu = nu + ldiso(id)%signf(ig,it)/ldiso(id)%sigf(ig,it)
          end do ! of it
          nu = nu	/ ldiso(id)%ntemp
          ldiso(id)%lvfis(:,ig) = ldiso(id)%lvfis(:,ig) / nu
        END DO ! of ig
      CASE('!SF7W+')
        DO ig = igresb, igrese
          DO it = 1, nrtemp
            READ(indev,*) igx,itx,ldiso(id)%wgtfis(1:nsubghel,it,ig)
          END DO ! of it
        END DO ! of ig
      CASE('!SA4S+')
        ldiso(id)%nlvflx = nflxhel
        DO ig = igresb, igrese
          READ(indev,*) igx,ldiso(id)%lvflx(1:nflxhel,ig)
        END DO ! of ig
        ldiso(id)%lvflx_log = dlog(ldiso(id)%lvflx)
      CASE('!SA4W+')
        DO ig = igresb, igrese
          DO it = 1, nrtemp
            READ(indev,*) igx,itx,farr2(1:nflxhel)
          END DO ! of it
        END DO ! of ig
      CASE('!CHI+')
        READ(indev,*) ldiso(id)%chi(1:nchihel)
        fdum = 1./sum(ldiso(id)%chi(1:nchihel))
        ldiso(id)%chi(1:nchihel) = ldiso(id)%chi(1:nchihel)*fdum
      CASE('!FIS+')
        READ(indev,*) ldiso(id)%kappa,ldiso(id)%kappa0
        READ(indev,*) ldiso(id)%beta(0:6)
        READ(indev,*) farr2(1:nburhel)
      CASE('!DCY+')
        READ(indev,*) ldiso(id)%dcy
      CASE('!TP3+')
        READ(indev,*) ldiso(id)%p1temp(1:np1temp)
      CASE('!SP1+')
        DO ig = 1, ng
          DO it = 1, np1temp
            READ(indev,*) igx,itx,ldiso(id)%sigsp1(ig,it),ib,ie,farr2(ib:ie)
            ldiso(id)%smp1(ig,it)%ib = ib
            ldiso(id)%smp1(ig,it)%ie = ie
            ALLOCATE(ldiso(id)%smp1(ig,it)%from(ib:ie))
            ldiso(id)%smp1(ig,it)%from(ib:ie) = farr2(ib:ie)
          END DO ! of it
        END DO ! of ig
      CASE('!SP2+')
        DO ig = 1, ng
          DO it = 1, np1temp
            READ(indev,*) igx,itx,ldiso(id)%sigsp2(ig,it),ib,ie,farr2(ib:ie)
            ldiso(id)%smp2(ig,it)%ib = ib
            ldiso(id)%smp2(ig,it)%ie = ie
            ALLOCATE(ldiso(id)%smp2(ig,it)%from(ib:ie))
            ldiso(id)%smp2(ig,it)%from(ib:ie) = farr2(ib:ie)
          END DO ! of it
        END DO ! of ig
      CASE('!SP3+')
        DO ig = 1, ng
          DO it = 1, np1temp
            READ(indev,*) igx,itx,ldiso(id)%sigsp3(ig,it),ib,ie,farr2(ib:ie)
            ldiso(id)%smp3(ig,it)%ib = ib
            ldiso(id)%smp3(ig,it)%ie = ie
            ALLOCATE(ldiso(id)%smp3(ig,it)%from(ib:ie))
            ldiso(id)%smp3(ig,it)%from(ib:ie) = farr2(ib:ie)
          END DO ! of it
        END DO ! of ig
      CASE('!N2N+')
        READ(indev,*) ldiso(id)%sign2n(1:ng)
      CASE('!N3N+')
        READ(indev,*) ldiso(id)%sign3n(1:ng)
      CASE DEFAULT
        CALL terminate('BLOCK NAME '//trim(blockname)//' NOT Allowed...')
    END SELECT
  END DO
200 CONTINUE
  nlvmax = nsubghel
  CLOSE(indev)

  !Scattering Matrix Post-Processing
  DO i = 1, ntiso
    IF(ldiso(i)%np1temp == 0) CALL dmalloc0(ldiso(i)%sigsp1,1,ng,1,ntemp)
    IF(i > nelrhel) THEN
      DO it = 1, ldiso(i)%ntemp
        DO ig = 1, ng
          ldiso(i)%sm(ig,it)%ib = ig
          ldiso(i)%sm(ig,it)%ie = ig
          CALL dmalloc0(ldiso(i)%sm(ig,it)%from,ig,ig)
          ldiso(i)%sm(ig,it)%from(ig) = 0.
        END DO ! of ig
      END DO ! of it
    END IF
    DO it = 1, ldiso(i)%ntemp
      DO jg = 1, ng
        DO ig = 1, ng
          IF(jg >= ldiso(i)%sm(ig,it)%ib .AND. jg <= ldiso(i)%sm(ig,it)%ie) THEN
            ldiso(i)%sm(jg,it)%ioutsb = ig
            exit
          END IF
        END DO ! of ig
        DO ig = ng, 1, -1
          IF(jg >= ldiso(i)%sm(ig,it)%ib .AND. jg <= ldiso(i)%sm(ig,it)%ie) THEN
            ldiso(i)%sm(jg,it)%ioutse = ig
            exit
          END IF
        END DO ! of ig
        DO ig = ldiso(i)%sm(jg,it)%ib, ldiso(i)%sm(jg,it)%ie
          ldiso(i)%sigstr(ig,it) = ldiso(i)%sigstr(ig,it) + ldiso(i)%sm(jg,it)%from(ig)
        END DO ! of ig
      END DO ! of jg
      DO ig = 1, ng
        ldiso(i)%sigss(ig,it) = ldiso(i)%sm(ig,it)%from(ig) + (ldiso(i)%sigs(ig,it) - ldiso(i)%sigstr(ig,it))
      END DO ! of ig
    END DO ! of it
  END DO ! of i

  CALL nuclidmap
  RETURN
END SUBROUTINE ReadPLC

SUBROUTINE ReadPHL(indev,phlbfile,nele,ntiso,ngg)
!--- JSU EDIT 20170721

USE PARAM
USE ALLOCS
USE IOUTIL,            ONLY : OpenFile, ToUpper, Terminate
USE XSLIB_MOD,         ONLY : enbgam, phatom, GAMLIBDATA, nelmGAM, noggphl, nogghel
USE GamXSUtil,         ONLY : mappingphatom
IMPLICIT NONE

INTEGER :: indev ! input device
CHARACTER*256 :: phlbfile
INTEGER :: ngg, nele,nver, ntiso

CHARACTER*4 :: ext

CHARACTER*20 :: blockname

INTEGER :: nid,ntemp
INTEGER :: ityp,ifis,iradcap,iinel
CHARACTER*20  :: aid
REAL :: aw
LOGICAL :: lbin

INTEGER :: i,ig,ix,it,id,ic,iel
INTEGER :: igx,itx,nx1,nx2,jg,ip,ib,iy

REAL :: xtemp(1000)
REAL, POINTER :: PhoEnergybdry(:)

TYPE(GAMLIBDATA), POINTER :: elem

! VARIABLES TO HANDLE GAMMA TRANSPORT VARIABLES
INTEGER :: igg, imt, nmt = 4

CHARACTER*1 :: probe, probel
CHARACTER(1000) :: onelinel
CHARACTER(256) :: oneline

EQUIVALENCE(probe, oneline)
EQUIVALENCE(probel, onelinel)

PRINT*, 'PHL READ'
!!    1. open file
i = len(TRIM(phlbfile))
ext = phlbfile(i-3:i)
CALL toupper(ext)
IF (ext.eq.'PHLB') THEN
    lbin = .true.
ELSEIF (ext.eq.'PHLA') THEN
    lbin = .false.
ELSE
    CALL terminate('Extension of the multigroup library should be PHLA or PHLB.')
END IF

CLOSE(indev)

IF (lbin) THEN
  CALL openfile(indev,TRUE,TRUE,FALSE,phlbfile)
  READ(indev) noggphl, nele, nver
  nelmGAM = nele
  IF (ngg.NE.noggphl) CALL terminate('The Photon Group Number is Different in PHLA from MLX')
  ALLOCATE(PhoEnergybdry(ngg+1))
  ALLOCATE(phatom(nele))
  READ(indev) (PhoEnergybdry(ig), ig=1,ngg+1)
  DO ig = 1, ngg+1
    IF ((PhoEnergybdry(ig)-enbgam(ig))/enbgam(ig).GT.1E-5) CALL terminate('The Photon Group Boundary is Different in PHLA from MLX')
  END DO
  DO iel = 1,nele
    elem => phatom(iel)
    READ(indev) i, nid, aw
    elem%nid = nid
    elem%aw = aw
    ! Memory Allocation
    CALL dmalloc(elem%siga, ngg)
    CALL dmalloc(elem%sigtr, ngg)
    CALL dmalloc(elem%sigs, ngg)
    CALL dmalloc(elem%sigss, ngg)
    CALL dmalloc(elem%sigstr, ngg)
    CALL dmalloc(elem%kerma, ngg)
    ALLOCATE(elem%sm(ngg),elem%smp1(ngg),elem%smp2(ngg),elem%smp3(ngg))
    CALL dmalloc(elem%sigsp1,ngg)
    CALL dmalloc(elem%sigsp2,ngg)
    CALL dmalloc(elem%sigsp3,ngg)
    ! XSP+
    DO ig = 1, ngg
      READ(indev) nx1, nx2
      CALL dmalloc0(elem%sm(ig)%from,nx1,nx2)
      READ(indev) elem%kerma(ig), elem%siga(ig), elem%sigtr(ig), elem%sigs(ig), elem%sm(ig)%from(nx1:nx2)
      elem%sm(ig)%ib = nx1
      elem%sm(ig)%ie = nx2
    END DO
    ! PS1+
    DO ig=1,ngg
      READ(indev) nx1,nx2
      CALL dmalloc0(elem%smp1(ig)%from,nx1,nx2)
      READ(indev) elem%sigsp1(ig), elem%smp1(ig)%from(nx1:nx2)
      elem%smp1(ig)%ib=nx1;      elem%smp1(ig)%ie=nx2
    END DO
    ! PS2+
    DO ig=1,ngg
      READ(indev) nx1,nx2
      CALL dmalloc0(elem%smp2(ig)%from,nx1,nx2)
      READ(indev) elem%sigsp2(ig), elem%smp2(ig)%from(nx1:nx2)
      elem%smp2(ig)%ib=nx1;      elem%smp2(ig)%ie=nx2
    END DO
    ! PS3+
    DO ig=1,ngg
      READ(indev) nx1,nx2
      CALL dmalloc0(elem%smp3(ig)%from,nx1,nx2)
      READ(indev) elem%sigsp3(ig),elem%smp3(ig)%from(nx1:nx2)
      elem%smp3(ig)%ib=nx1;      elem%smp3(ig)%ie=nx2
    END DO
  END DO
ELSE
  CALL openfile(indev,TRUE,FALSE,FALSE,phlbfile)
  READ(indev, *) ! !Dimension [ngg, nele, nver]
  READ(indev, '(I15, 2I6)') noggphl, nele, nver
  IF (ngg.NE.noggphl) CALL terminate('The Photon Group Number is Different in PHLA from MLX')
  nelmGAM = nele
  ALLOCATE(PhoEnergybdry(ngg+1))
  ALLOCATE(phatom(nele))
  READ(indev, *) ! Photon Group Boundary
  READ(indev, '(10ES14.7)') (PhoEnergybdry(ig), ig=1,ngg+1)
  DO ig = 1, ngg+1
    IF ((PhoEnergybdry(ig)-enbgam(ig))/enbgam(ig).GT.1E-5) CALL terminate('The Photon Group Boundary is Different in PHLA from MLX')
  END DO
  READ(indev, *) ! DIR  [ i,    nid,        amass,            aid]
  DO iel = 1, nele
    READ(indev, *) ! Elements List
  END DO
  DO iel = 1,nele
    elem => phatom(iel)
    READ(indev, *) aid, i, nid, aw, aid
    elem%nid = nid
    elem%aw = aw
    elem%aid = aid
    ! Memory Allocation
    CALL dmalloc(elem%siga, ngg)
    CALL dmalloc(elem%sigtr, ngg)
    CALL dmalloc(elem%sigs, ngg)
    CALL dmalloc(elem%sigss, ngg)
    CALL dmalloc(elem%sigstr, ngg)
    CALL dmalloc(elem%kerma, ngg)
    ALLOCATE(elem%sm(ngg),elem%smp1(ngg),elem%smp2(ngg),elem%smp3(ngg))
    CALL dmalloc(elem%sigsp1,ngg)
    CALL dmalloc(elem%sigsp2,ngg)
    CALL dmalloc(elem%sigsp3,ngg)
    READ(indev, *) ! XSP+
    DO ig = 1, ngg
      READ(indev, *) igx, elem%kerma(ig), elem%siga(ig), elem%sigtr(ig), elem%sigs(ig), nx1, nx2, (xtemp(jg), jg = nx1, nx2)
      CALL dmalloc0(elem%sm(ig)%from,nx1,nx2)
      elem%sm(ig)%ib = nx1
      elem%sm(ig)%ie = nx2
      elem%sm(ig)%from(nx1:nx2)=xtemp(nx1:nx2)
    END DO
    READ(indev, *) ! PS1+
    DO ig=1,ngg
      READ(indev,*) igx,elem%sigsp1(ig),nx1,nx2,(xtemp(jg),jg=nx1,nx2)
      elem%smp1(ig)%ib=nx1;      elem%smp1(ig)%ie=nx2
      CALL dmalloc0(elem%smp1(ig)%from,nx1,nx2)
      DO jg=nx1,nx2
        elem%smp1(ig)%from(jg)=xtemp(jg)
      END DO
    END DO
    READ(indev, *) ! PS2+
    DO ig=1,ngg
      READ(indev,*) igx,elem%sigsp2(ig),nx1,nx2,(xtemp(jg),jg=nx1,nx2)
      elem%smp2(ig)%ib=nx1;      elem%smp2(ig)%ie=nx2
      CALL dmalloc0(elem%smp2(ig)%from,nx1,nx2)
      DO jg=nx1,nx2
        elem%smp2(ig)%from(jg)=xtemp(jg)
      END DO
    END DO
    READ(indev, *) ! PS3+
    DO ig=1,ngg
      READ(indev,*) igx,elem%sigsp3(ig),nx1,nx2,(xtemp(jg),jg=nx1,nx2)
      elem%smp3(ig)%ib=nx1;      elem%smp3(ig)%ie=nx2
      CALL dmalloc0(elem%smp3(ig)%from,nx1,nx2)
      DO jg=nx1,nx2
        elem%smp3(ig)%from(jg)=xtemp(jg)
      END DO
    END DO
  END DO
END IF
CLOSE(indev)
! End of Reading Data... ************

! scattering correction data calculation.
DO iel = 1, nele
  elem => phatom(iel)
  elem%sigstr=0._8
  DO ig=1,ngg
    DO igx=1,ngg
      IF(ig.GE.elem%sm(igx)%ib .AND. ig.LE.elem%sm(igx)%ie) THEN
        elem%sm(ig)%ioutsb=igx
        EXIT
      END IF
    END DO
    DO igx=ngg,1,-1
      IF(ig.GE.elem%sm(igx)%ib .AND. ig.LE.elem%sm(igx)%ie) THEN
        elem%sm(ig)%ioutse=igx
        EXIT
      END IF
    END DO
    DO igx=elem%sm(ig)%ib,elem%sm(ig)%ie
      elem%sigstr(igx)=elem%sigstr(igx)+elem%sm(ig)%from(igx)
    END DO
  END DO
  DO ig = 1, ngg
      elem%sigss(ig) = elem%sm(ig)%from(ig) + (elem%sigs(ig) - elem%sigstr(ig)) ! Un-corrected Scattering XS
  END DO
END DO

!!!    7. mapping
CALL mappingphatom(nele,ntiso)
!
RETURN

END SUBROUTINE ReadPHL

end Module
