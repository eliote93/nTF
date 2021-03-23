SUBROUTINE ReadISOTXS(indev,filenm,ng_out,ntiso_out,scatord_out)
    USE ioutil,         ONLY : OpenFile
    USE CNTL,           ONLY : nTracerCntl
    USE XSLIB_MOD,      ONLY : ldiso,nelthel,nelrhel,nreshel,nburhel,mapnucl,mapnuclp13,mapfis,&
                               noghel,nofghel,norghel,notghel,nchihel,enbhel, SCATMAT
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: indev
    CHARACTER*(*),INTENT(IN) :: filenm
    INTEGER,INTENT(OUT) :: ng_out,ntiso_out,scatord_out

    INTEGER :: ng,ntiso,maxup,maxdn,maxord,ichist,nscmax,nsblok,sctord
    INTEGER :: ig,jg,iso,ib,ie,il,i,j,k,scttyp,iord,it,igx
    INTEGER :: kbr,ichi,ifis,ialf,inp,in2n,ind,int,ltot,ltrn,istrpd,nscr
    INTEGER,POINTER :: isspec(:),loca(:),idsct(:),lord(:),jband(:,:),ijj(:,:),isopec(:)
    REAL(4) :: amass,efiss,ecapt,temp,sigpot,numden,fissum,tmp
    REAL(4),POINTER :: chi(:,:),erg(:),vel(:)
    REAL(4),POINTER :: strpl(:,:),stotpl(:,:),sngam(:),sfis(:),snutot(:),schiso(:),&
                       snalf(:),snp(:),sn2n(:),snd(:),snt(:),strpd(:,:),scr(:),chiiso(:,:)
    REAL(4),POINTER :: sctmat(:,:,:)
    REAL(8),POINTER :: sctxs(:)
    CHARACTER*6 :: hid(12),absid,ident,mat
    CHARACTER*6,POINTER :: aid(:)
    LOGICAL,POINTER :: ltotmat(:)
    TYPE(SCATMAT),POINTER :: SM(:)

    ng_out = 0
    ntiso_out = 0
    scatord_out = 0

    DO i = 1, 100000
        mapnucl(i) = 0
        mapnuclp13(i) = 0
        mapfis(i) = 0
    ENDDO

    !1. Open a XS Library File
    CLOSE(InDev)
    CALL OpenFile(indev,.TRUE.,.FALSE.,.FALSE.,filenm)

    !Card Type 0, Skip this Line
    READ(indev,*)

    !Card Type 1, File Control (1D Record)
100 FORMAT (' 1D ',8i6)
    READ(indev,100) ng,ntiso,maxup,maxdn,maxord,ichist,nscmax,nsblok

    !Card Type 2, File Data (2D Record)
200 FORMAT (' 2D ','*',11a6,'*'/'*',a6,'*',9(1x,a6)/(10(1x,a6)))
210 FORMAT(1P,6E12.5)
220 FORMAT(12I6)
    ALLOCATE(aid(ntiso))
    ALLOCATE(erg(ng+1),vel(ng),loca(1:ntiso))
    READ(indev,200) hid(1:12),aid(1:ntiso)
    IF(ichist == 1) THEN
      ALLOCATE(chi(ng,ichist))
      READ(indev,210) chi(1:ng,ichist)
    ENDIF
    READ(indev,210) vel(1:ng),erg(1:ng+1)
    READ(indev,220) loca(1:ntiso)

    ALLOCATE(enbhel(ng+1))
    enbhel(1:ng+1) = erg(1:ng+1)

    !Card Type 3, File-wide Chi Data (3D Record)
300 FORMAT (' 3D ',1p5e12.5/1p,(6e12.5))
    IF(ichist > 1) THEN
      ALLOCATE(chi(ng,ichist),isspec(ng))
      READ(indev,300) ((chi(ichi,ig),ichi=1,ichist),ig=1,ng),(isspec(ig),ig=1,ng)
    ENDIF

    ALLOCATE(idsct(nscmax),lord(nscmax),jband(ng,nscmax),ijj(ng,nscmax))
    ALLOCATE(strpl(ng,10),stotpl(ng,10),sngam(ng),sfis(ng),snutot(ng))
    ALLOCATE(schiso(ng),snalf(ng),snp(ng),sn2n(ng),snd(ng),snt(ng))
    ALLOCATE(ldiso(ntiso))
    ALLOCATE(sctmat(ng,ng,0:maxord))
    ALLOCATE(ltotmat(0:maxord))
    DO iso = 1, ntiso
        !Card Type 4, Isotope Control and Group Independent Datg (4D Record)
400     FORMAT (' 4D ',3(1x,a6)/1p,(6e12.5)/0p,(12i6))
        READ(indev,400) absid,ident,mat,amass,efiss,ecapt,temp,sigpot,numden,&
                        kbr,ichi,ifis,ialf,inp,in2n,ind,int,ltot,ltrn,istrpd,&
                        idsct(1:nscmax),lord(1:nscmax),jband(1:ng,1:nscmax),ijj(1:ng,1:nscmax)

        ldiso(iso)%aid = aid(iso)
        ldiso(iso)%aw  = amass
        ldiso(iso)%ntemp = 1
        ldiso(iso)%nrtemp = 0
        ldiso(iso)%ifis = ifis
        ALLOCATE(ldiso(iso)%temp(1))
        ldiso(iso)%temp(1) = temp
        IF(ifis == 1) ldiso(iso)%ityp = 3
        ldiso(iso)%np1temp = 0
        ldiso(iso)%kappa = efiss


        !Card Type 5, Principal Cross Section
        nscr = ltrn*ng + ltot*ng + ng + 2*ifis*ng + ialf*ng + inp*ng + in2n*ng + ind*ng + int*ng + istrpd*ng
        IF(ichi == 1) nscr = nscr + ng
        ALLOCATE(scr(nscr))
500     FORMAT (' 5d ',1p5e12.5/1p,(6e12.5))
        READ(indev,500) scr(1:nscr)
        ie = 0
        !PL weighted transport cross section
        DO il = 1, ltrn
            ib = ie + 1
            ie = ib + ng - 1
            strpl(1:ng,il) = scr(ib:ie)
        ENDDO
        !PL weighted total cross section
        DO il = 1, ltot
            ib = ie + 1
            ie = ib + ng - 1
            stotpl(1:ng,il) = scr(ib:ie)
        ENDDO
        !(n,gamma) cross section, capture
        ib = ie + 1
        ie = ib + ng - 1
        sngam(1:ng) = scr(ib:ie)
        IF(ifis > 0) THEN
            !fission cross section
            ib = ie + 1
            ie = ib + ng - 1
            sfis(1:ng) = scr(ib:ie)
            !total neutron yield per fission
            ib = ie + 1
            ie = ib + ng - 1
            snutot(1:ng) = scr(ib:ie)
        ENDIF
        !isotope chi
        IF(ichi == 1) THEN
            ib = ie + 1
            ie = ib + ng - 1
            schiso(1:ng) = scr(ib:ie)
        ENDIF
        !(n,alpha) cross section
        IF(ialf == 1) THEN
            ib = ie + 1
            ie = ib + ng - 1
            snalf(1:ng) = scr(ib:ie)
        ENDIF
        !(n,proton) cross section
        IF(inp == 1) THEN
            ib = ie + 1
            ie = ib + ng - 1
            snp(1:ng) = scr(ib:ie)
        ENDIF
        !(n,2n) cross section
        IF(in2n == 1) THEN
            ib = ie + 1
            ie = ib + ng - 1
            sn2n(1:ng) = scr(ib:ie)
        ENDIF
        !(n,d) cross section
        IF(ind == 1) THEN
            ib = ie + 1
            ie = ib + ng - 1
            snd(1:ng) = scr(ib:ie)
        ENDIF
        !(n,t) cross section
        IF(int == 1) THEN
            ib = ie + 1
            ie = ib + ng - 1
            snt(1:ng) = scr(ib:ie)
        ENDIF
        !coordinate direction i trnasport cross section
        IF(istrpd > 0) ALLOCATE(strpd(ng,istrpd))
        DO il = 1, istrpd
            ib = ie + 1
            ie = ib + ng - 1
            strpd(1:ng,il) = scr(ib:ie)
        ENDDO
        DEALLOCATE(scr)

        !Save cross section to ldiso(iso)
        ALLOCATE(ldiso(iso)%sigtr(ng,1),ldiso(iso)%siga(ng,1),ldiso(iso)%sigs(ng,1),&
                 ldiso(iso)%sigf(ng,1),ldiso(iso)%signf(ng,1),ldiso(iso)%chi(ng),&
                 ldiso(iso)%sigstr(ng,1), ldiso(iso)%sigss(ng,1))
        FORALL(ig=1:ng)
            ldiso(iso)%sigtr(ig,1) = 0.
            ldiso(iso)%siga (ig,1) = 0.
            ldiso(iso)%sigs (ig,1) = 0.
            ldiso(iso)%sigf (ig,1) = 0.
            ldiso(iso)%signf(ig,1) = 0.
            ldiso(iso)%chi  (ig)   = 0.
            ldiso(iso)%sigstr(ig,1) = 0.
            ldiso(iso)%sigss(ig,1) = 0.
        END FORALL
        ldiso(iso)%siga (1:ng,1) = stotpl(1:ng,1)
        ldiso(iso)%sigtr(1:ng,1) = strpl(1:ng,1)
        IF(ifis > 0) THEN
            ldiso(iso)%sigf(1:ng,1) = sfis(1:ng)
            ldiso(iso)%signf(1:ng,1) = snutot(1:ng)*sfis(1:ng)
        ENDIF
        IF(ichi == 1) ldiso(iso)%chi(1:ng) = schiso(1:ng)

        !Card Type 6, Isotope Chi Data
        IF(ichi > 1) THEN
600         FORMAT(' 6d ',1p5e12.5/1p,(6e12.5))
610         FORMAT(12i6)
            ALLOCATE(CHIISO(ichi,ng),ISOPEC(ng))
            READ(indev,600) CHIISO(1:ichi,1:ng)
            READ(indev,610) ISOPEC(1:ng)

            fissum = 0.
            DO ig = 1, ng
                tmp = snutot(ig)*sfis(ig)
                fissum = fissum + tmp
                DO jg = 1, ng
                    ldiso(iso)%chi(jg) = ldiso(iso)%chi(jg) + tmp*chiiso(isopec(ig),jg)
                ENDDO
            ENDDO
            tmp = 0.
            DO ig = 1, ng
                ldiso(iso)%chi(ig) = ldiso(iso)%chi(ig) / fissum
                tmp = tmp + ldiso(iso)%chi(ig)
            ENDDO
            ldiso(iso)%chi = ldiso(iso)%chi/tmp
        ENDIF

        !Card Type 7, Scattering Sub-block
        forall(ig=1:ng,jg=1:ng,il=0:maxord) sctmat(ig,jg,il) = 0.
        ltotmat = .FALSE.
        sctord = 0
        DO i = 1, nscmax
            IF(lord(i) == 0) CYCLE
700         FORMAT(' 7d ',1p5e12.5/1p,(6e12.5))
            nscr = 0
            DO ig = 1, ng
                nscr = nscr + jband(ig,i)
            ENDDO
            nscr = nscr*lord(i)
            ALLOCATE(scr(nscr))
            READ(indev,700) scr(1:nscr)

            !save scattering matrix
            scttyp = idsct(i)/100
            il = idsct(i) - scttyp*100 - 1
            DO iord = 1, lord(i)
                il = il + 1
                sctord = max(sctord,il)
                IF(scttyp == 0) THEN
                    ltotmat(il) = .TRUE.
                    k = 0
                    DO ig = 1, ng
                        jg = ig - ijj(ig,i) + 1
                        DO j = 1, jband(ig,i)
                            k = k + 1
                            sctmat(jg,ig,il) = scr(k)
                            jg = jg -1
                        ENDDO
                    ENDDO
                ELSE
                    IF(ltotmat(il)) THEN
                        DEALLOCATE(scr)
                        CYCLE
                    ENDIF
                    k = 0
                    DO ig = 1, ng
                        jg = ig - ijj(ig,i) + 1
                        DO j = 1, jband(ig,i)
                            k = k + 1
                            sctmat(jg,ig,il) = sctmat(jg,ig,il) + scr(k)
                            jg = jg - 1
                        ENDDO
                    ENDDO
                ENDIF
            ENDDO
            DEALLOCATE(scr)
        ENDDO

        IF(sctord > 0) THEN
            ldiso(iso)%np1temp = 1
            ALLOCATE(ldiso(iso)%p1temp(1))
            ldiso(iso)%p1temp(1) = temp
        ENDIF

        DO il = 0, sctord
            IF(il == 0) THEN
                ALLOCATE(ldiso(iso)%sm(ng,1))
                SM => ldiso(iso)%sm(:,1)
                sctxs => ldiso(iso)%sigs(:,1)
            ELSEIF(il == 1) THEN
                ALLOCATE(ldiso(iso)%smp1(ng,1))
                SM => ldiso(iso)%smp1(:,1)
                ALLOCATE(ldiso(iso)%sigsp1(ng,1))
                sctxs => ldiso(iso)%sigsp1(:,1)
            ELSEIF(il == 2) THEN
                ALLOCATE(ldiso(iso)%smp2(ng,1))
                SM => ldiso(iso)%smp2(:,1)
                ALLOCATE(ldiso(iso)%sigsp2(ng,1))
                sctxs => ldiso(iso)%sigsp2(:,1)
            ELSEIF(il == 3) THEN
                ALLOCATE(ldiso(iso)%smp3(ng,1))
                SM => ldiso(iso)%smp3(:,1)
                ALLOCATE(ldiso(iso)%sigsp3(ng,1))
                sctxs => ldiso(iso)%sigsp3(:,1)
            ELSEIF(il == 4) THEN
                ALLOCATE(ldiso(iso)%smp4(ng,1))
                SM => ldiso(iso)%smp4(:,1)
                ALLOCATE(ldiso(iso)%sigsp4(ng,1))
                sctxs => ldiso(iso)%sigsp4(:,1)
            ELSEIF(il == 5) THEN
                ALLOCATE(ldiso(iso)%smp5(ng,1))
                SM => ldiso(iso)%smp5(:,1)
                ALLOCATE(ldiso(iso)%sigsp5(ng,1))
                sctxs => ldiso(iso)%sigsp5(:,1)
            ELSE
                stop 'readISOTXS.f90, legendre expansion order is larger than 5 !!'
            ENDIF

            DO ig = 1, ng
                DO jg = 1, ng
                    IF(sctmat(jg,ig,il) == 0.) CYCLE
                    ib = jg
                    EXIT
                END DO
                DO jg = ng, 1, -1
                    IF(sctmat(jg,ig,il) == 0.) CYCLE
                    ie = jg
                    EXIT
                END DO
                SM(ig)%ib = ib
                SM(ig)%ie = ie
                ALLOCATE(SM(ig)%from(ib:ie))
                SM(ig)%from(ib:ie) = sctmat(ib:ie,ig,il)

                sctxs(ig) = sum(sctmat(ig,1:ng,il))
                IF(il == 0) ldiso(iso)%siga(ig,1) = ldiso(iso)%siga(ig,1) - ldiso(iso)%sigs(ig,1)
            ENDDO
        ENDDO

        !!Transport Correction for P0 Scattering Matrix
        DO ig = 1, ng
            ldiso(iso)%sm(ig,1)%from(ig) = ldiso(iso)%sm(ig,1)%from(ig) - (stotpl(ig,1) - strpl(ig,1))
        ENDDO

        CONTINUE
!        DO ig = 1, ng
!            WRITE(50,'(100ES16.6)') sctmat(ig,1:ng,0)
!        ENDDO
        IF(istrpd > 0) DEALLOCATE(strpd)

        do it = 1, 1
            ldiso(iso)%sigstr(:,it) = 0._8
            do ig = 1, ng
                do igx = 1, ng
                    if(ig.ge.ldiso(iso)%sm(igx,it)%ib .and. ig.le.ldiso(iso)%sm(igx,it)%ie) then
                        ldiso(iso)%sm(ig,it)%ioutsb = igx
                        exit
                    endif
                enddo
                do igx = ng, 1, -1
                    if(ig.ge.ldiso(iso)%sm(igx,it)%ib .and. ig.le.ldiso(iso)%sm(igx,it)%ie) then
                        ldiso(iso)%sm(ig,it)%ioutse = igx
                        exit
                    endif
                enddo
                do igx = ldiso(iso)%sm(ig,it)%ib, ldiso(iso)%sm(ig,it)%ie
                    ldiso(iso)%sigstr(igx,it) = ldiso(iso)%sigstr(igx,it) + ldiso(iso)%sm(ig,it)%from(igx)
                enddo
            enddo
        enddo
        do it = 1, 1
            do ig = 1, ng
                ldiso(iso)%sigss(ig,it) = ldiso(iso)%sm(ig,it)%from(ig) + (ldiso(iso)%sigs(ig,it) - ldiso(iso)%sigstr(ig,it))
            enddo
        enddo

    ENDDO !of iso
    DEALLOCATE(strpl,stotpl,sngam,sfis,snutot,schiso,snalf,snp,sn2n,snd,snt)
    DEALLOCATE(sctmat,ltotmat)

    ng_out = ng
    ntiso_out = ntiso
    scatord_out = maxord

    nelthel = ntiso
    nelrhel = ntiso
    nreshel = 0
    nburhel = 0

    noghel = ng
    DO ig = 1, ng
        IF(erg(ig) <= 5.0) EXIT
    ENDDO
    nofghel = ig-1
    notghel = noghel - nofghel
    norghel = 0
    nchihel = noghel

    IF(ichist == 1) DEALLOCATE(chi)
    IF(ichist > 1)  DEALLOCATE(chi,isspec)
    DEALLOCATE(erg,vel,loca)
    DEALLOCATE(idsct,lord,jband,ijj)
    RETURN
END SUBROUTINE

SUBROUTINE ReadGAMISO(indev,filenm,ng_out,ntiso_out)
    USE ioutil,         ONLY : OpenFile
    USE CNTL,           ONLY : nTracerCntl
    USE XSLIB_MOD,      ONLY : enbgam, phatom, GAMLIBDATA, nelmGAm, noggphl, scatmat
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: indev
    CHARACTER*(*),INTENT(IN) :: filenm
    INTEGER,INTENT(OUT) :: ng_out,ntiso_out

    INTEGER :: ng,ntiso,maxup,maxdn,maxord,ichist,nscmax,nsblok,sctord
    INTEGER :: ig,jg,iso,ib,ie,il,i,j,k,scttyp,iord,it,igx
    INTEGER :: kbr,ichi,ifis,ialf,inp,in2n,ind,int,ltot,ltrn,istrpd,nscr
    INTEGER,POINTER :: isspec(:),loca(:),idsct(:),lord(:),jband(:,:),ijj(:,:),isopec(:)
    REAL(4) :: amass,efiss,ecapt,temp,sigpot,numden,fissum,tmp
    REAL(4),POINTER :: chi(:,:),erg(:),vel(:)
    REAL(4),POINTER :: strpl(:,:),stotpl(:,:),sngam(:),sfis(:),snutot(:),schiso(:),&
                       snalf(:),snp(:),sn2n(:),snd(:),snt(:),strpd(:,:),scr(:),chiiso(:,:)
    REAL(4),POINTER :: sctmat(:,:,:)
    REAL(8),POINTER :: sctxs(:)
    CHARACTER*6 :: hid(12),absid,ident,mat
    CHARACTER*6,POINTER :: aid(:)
    LOGICAL,POINTER :: ltotmat(:)
    TYPE(SCATMAT),POINTER :: SM(:)

          nTRACERCNTL%lExplicitKappa = nTracerCntl%lGamma
    ng_out = 0
    ntiso_out = 0

    !1. Open a XS Library File
    CLOSE(InDev)
    CALL OpenFile(indev,.TRUE.,.FALSE.,.FALSE.,filenm)

    !Card Type 0, Skip this Line
    READ(indev,*)

    !Card Type 1, File Control (1D Record)
100 FORMAT (' 1D ',8i6)
    READ(indev,100) ng,ntiso,maxup,maxdn,maxord,ichist,nscmax,nsblok

    !Card Type 2, File Data (2D Record)
200 FORMAT (' 2D ','*',11a6,'*'/'*',a6,'*',9(1x,a6)/(10(1x,a6)))
210 FORMAT(1P,6E12.5)
220 FORMAT(12I6)
    ALLOCATE(aid(ntiso))
    ALLOCATE(erg(ng+1),vel(ng),loca(1:ntiso))
    READ(indev,200) hid(1:12),aid(1:ntiso)
    IF(ichist == 1) THEN
      ALLOCATE(chi(ng,ichist))
      READ(indev,210) chi(1:ng,ichist)
    ENDIF
    READ(indev,210) vel(1:ng),erg(1:ng+1)
    READ(indev,220) loca(1:ntiso)

    ALLOCATE(enbgam(ng+1))
    enbgam(1:ng+1) = erg(1:ng+1)

    ALLOCATE(idsct(nscmax),lord(nscmax),jband(ng,nscmax),ijj(ng,nscmax))
    ALLOCATE(strpl(ng,10),stotpl(ng,10),sngam(ng),sfis(ng),snutot(ng))
    ALLOCATE(schiso(ng),snalf(ng),snp(ng),sn2n(ng),snd(ng),snt(ng))
    ALLOCATE(phatom(ntiso))
    ALLOCATE(sctmat(ng,ng,0:maxord))
    ALLOCATE(ltotmat(0:maxord))
    DO iso = 1, ntiso
        !Card Type 4, Isotope Control and Group Independent Datg (4D Record)
400     FORMAT (' 4D ',3(1x,a6)/1p,(6e12.5)/0p,(12i6))
        READ(indev,400) absid,ident,mat,amass,efiss,ecapt,temp,sigpot,numden,&
                        kbr,ichi,ifis,ialf,inp,in2n,ind,int,ltot,ltrn,istrpd,&
                        idsct(1:nscmax),lord(1:nscmax),jband(1:ng,1:nscmax),ijj(1:ng,1:nscmax)

        phatom(iso)%aid = aid(iso)
        phatom(iso)%aw  = amass

        !Card Type 5, Principal Cross Section
        nscr = ltrn*ng + ltot*ng + ng + 2*ifis*ng + ialf*ng + inp*ng + in2n*ng + ind*ng + int*ng + istrpd*ng
        ALLOCATE(scr(nscr))
500     FORMAT (' 5d ',1p5e12.5/1p,(6e12.5))
        READ(indev,500) scr(1:nscr)
        ie = 0
        !PL weighted transport cross section
        DO il = 1, ltrn
            ib = ie + 1
            ie = ib + ng - 1
            strpl(1:ng,il) = scr(ib:ie)
        ENDDO
        !PL weighted total cross section
        DO il = 1, ltot
            ib = ie + 1
            ie = ib + ng - 1
            stotpl(1:ng,il) = scr(ib:ie)
        ENDDO
        !(n,gamma) cross section, capture
        ib = ie + 1
        ie = ib + ng - 1
        sngam(1:ng) = scr(ib:ie)
        IF(ifis > 0) THEN
            !fission cross section
            ib = ie + 1
            ie = ib + ng - 1
            sfis(1:ng) = scr(ib:ie)
            !total neutron yield per fission
            ib = ie + 1
            ie = ib + ng - 1
            snutot(1:ng) = scr(ib:ie)
        ENDIF
        !isotope chi
        IF(ichi == 1) THEN
            ib = ie + 1
            ie = ib + ng - 1
            schiso(1:ng) = scr(ib:ie)
        ENDIF
        !(n,alpha) cross section
        IF(ialf == 1) THEN
            ib = ie + 1
            ie = ib + ng - 1
            snalf(1:ng) = scr(ib:ie)
        ENDIF
        !(n,proton) cross section
        IF(inp == 1) THEN
            ib = ie + 1
            ie = ib + ng - 1
            snp(1:ng) = scr(ib:ie)
        ENDIF
        !(n,2n) cross section
        IF(in2n == 1) THEN
            ib = ie + 1
            ie = ib + ng - 1
            sn2n(1:ng) = scr(ib:ie)
        ENDIF
        !(n,d) cross section
        IF(ind == 1) THEN
            ib = ie + 1
            ie = ib + ng - 1
            snd(1:ng) = scr(ib:ie)
        ENDIF
        !(n,t) cross section
        IF(int == 1) THEN
            ib = ie + 1
            ie = ib + ng - 1
            snt(1:ng) = scr(ib:ie)
        ENDIF
        !coordinate direction i trnasport cross section
        IF(istrpd > 0) ALLOCATE(strpd(ng,istrpd))
        DO il = 1, istrpd
            ib = ie + 1
            ie = ib + ng - 1
            strpd(1:ng,il) = scr(ib:ie)
        ENDDO
        DEALLOCATE(scr)

        !Save cross section to phatom(iso)
        ALLOCATE(phatom(iso)%sigtr(ng),phatom(iso)%siga(ng),phatom(iso)%sigs(ng),&
                 phatom(iso)%sigstr(ng),phatom(iso)%sigss(ng))
        DO ig=1,ng
            phatom(iso)%sigtr(ig) = 0.
            phatom(iso)%siga (ig) = 0.
            phatom(iso)%sigs (ig) = 0.
            phatom(iso)%sigstr(ig) = 0.
            phatom(iso)%sigss(ig) = 0.
        END DO 
        phatom(iso)%siga (1:ng) = stotpl(1:ng,1)
        phatom(iso)%sigtr(1:ng) = strpl(1:ng,1)

        !Card Type 7, Scattering Sub-block
        forall(ig=1:ng,jg=1:ng,il=0:maxord) sctmat(ig,jg,il) = 0.
        ltotmat = .FALSE.
        sctord = 0
        DO i = 1, nscmax
            IF(lord(i) == 0) CYCLE
700         FORMAT(' 7d ',1p5e12.5/1p,(6e12.5))
            nscr = 0
            DO ig = 1, ng
                nscr = nscr + jband(ig,i)
            ENDDO
            nscr = nscr*lord(i)
            ALLOCATE(scr(nscr))
            READ(indev,700) scr(1:nscr)

            !save scattering matrix
            scttyp = idsct(i)/100
            il = idsct(i) - scttyp*100 - 1
            DO iord = 1, lord(i)
                il = il + 1
                sctord = max(sctord,il)
                IF(scttyp == 0) THEN
                    ltotmat(il) = .TRUE.
                    k = 0
                    DO ig = 1, ng
                        jg = ig - ijj(ig,i) + 1
                        DO j = 1, jband(ig,i)
                            k = k + 1
                            sctmat(jg,ig,il) = scr(k)
                            jg = jg -1
                        ENDDO
                    ENDDO
                ELSE
                    IF(ltotmat(il)) THEN
                        DEALLOCATE(scr)
                        CYCLE
                    ENDIF
                    k = 0
                    DO ig = 1, ng
                        jg = ig - ijj(ig,i) + 1
                        DO j = 1, jband(ig,i)
                            k = k + 1
                            sctmat(jg,ig,il) = sctmat(jg,ig,il) + scr(k)
                            jg = jg - 1
                        ENDDO
                    ENDDO
                ENDIF
            ENDDO
            DEALLOCATE(scr)
        ENDDO

        DO il = 0, sctord
            IF(il == 0) THEN
                ALLOCATE(phatom(iso)%sm(ng))
                SM => phatom(iso)%sm(:)
                sctxs => phatom(iso)%sigs(:)
            ELSEIF(il == 1) THEN
                ALLOCATE(phatom(iso)%smp1(ng))
                SM => phatom(iso)%smp1(:)
                ALLOCATE(phatom(iso)%sigsp1(ng))
                sctxs => phatom(iso)%sigsp1(:)
            ELSEIF(il == 2) THEN
                ALLOCATE(phatom(iso)%smp2(ng))
                SM => phatom(iso)%smp2(:)
                ALLOCATE(phatom(iso)%sigsp2(ng))
                sctxs => phatom(iso)%sigsp2(:)
            ELSEIF(il == 3) THEN
                ALLOCATE(phatom(iso)%smp3(ng))
                SM => phatom(iso)%smp3(:)
                ALLOCATE(phatom(iso)%sigsp3(ng))
                sctxs => phatom(iso)%sigsp3(:)
            ELSE
                stop 'readGAMISO, legendre expansion order is larger than 3 !!'
            ENDIF

            DO ig = 1, ng
                DO jg = 1, ng
                    IF(sctmat(jg,ig,il) == 0.) CYCLE
                    ib = jg
                    EXIT
                END DO
                DO jg = ng, 1, -1
                    IF(sctmat(jg,ig,il) == 0.) CYCLE
                    ie = jg
                    EXIT
                END DO
                SM(ig)%ib = ib
                SM(ig)%ie = ie
                ALLOCATE(SM(ig)%from(ib:ie))
                SM(ig)%from(ib:ie) = sctmat(ib:ie,ig,il)

                sctxs(ig) = sum(sctmat(ig,1:ng,il))
                IF(il == 0) phatom(iso)%siga(ig) = phatom(iso)%siga(ig) - phatom(iso)%sigs(ig)
            ENDDO
        ENDDO

        !!Transport Correction for P0 Scattering Matrix
        DO ig = 1, ng
            phatom(iso)%sm(ig)%from(ig) = phatom(iso)%sm(ig)%from(ig) - (stotpl(ig,1) - strpl(ig,1))
        ENDDO

        CONTINUE
!        DO ig = 1, ng
!            WRITE(50,'(100ES16.6)') sctmat(ig,1:ng,0)
!        ENDDO
        IF(istrpd > 0) DEALLOCATE(strpd)

        phatom(iso)%sigstr(:) = 0._8
        do ig = 1, ng
            do igx = 1, ng
                if(ig.ge.phatom(iso)%sm(igx)%ib .and. ig.le. phatom(iso)%sm(igx)%ie) then
                    phatom(iso)%sm(ig)%ioutsb = igx
                    exit
                endif
            enddo
            do igx = ng, 1, -1
                if(ig.ge.phatom(iso)%sm(igx)%ib .and. ig.le.phatom(iso)%sm(igx)%ie) then
                    phatom(iso)%sm(ig)%ioutse = igx
                    exit
                endif
            enddo
            do igx = phatom(iso)%sm(ig)%ib, phatom(iso)%sm(ig)%ie
                phatom(iso)%sigstr(igx) = phatom(iso)%sigstr(igx) + phatom(iso)%sm(ig)%from(igx)
            enddo
        enddo
        do it = 1, 1
            do ig = 1, ng
                phatom(iso)%sigss(ig) = phatom(iso)%sm(ig)%from(ig) + (phatom(iso)%sigs(ig) - phatom(iso)%sigstr(ig))
            enddo
        enddo

    ENDDO !of iso
    DEALLOCATE(strpl,stotpl,sngam,sfis,snutot,schiso,snalf,snp,sn2n,snd,snt)
    DEALLOCATE(sctmat,ltotmat)

    ng_out = ng
    ntiso_out = ntiso

    noggphl = ng
    nelmGAM = ntiso

    DEALLOCATE(erg,vel,loca)
    DEALLOCATE(idsct,lord,jband,ijj)
    RETURN
END SUBROUTINE

SUBROUTINE ReadPMATRX(indev, filenm)
USE ioutil,         ONLY : OpenFile
USE CNTL,           ONLY : nTracerCntl
USE XSLIB_MOD,      ONLY : libdata, GAMLIBDATA, ldiso, nelthel, phatom
IMPLICIT NONE 
INTEGER, INTENT(IN) :: indev
CHARACTER*(*), INTENT(IN) :: filenm

TYPE(libdata), POINTER :: lib
TYPE(GAMLIBDATA), POINTER :: glib
REAL, POINTER ::  prod(:,:), gh(:), heat(:), damage(:), en(:), eg(:)
INTEGER, POINTER :: mat(:)
CHARACTER(8), POINTER :: cname(:)
INTEGER :: nclpse, ggroup, ngroup, ipltfg, niso, idose, maxlrd, ncomp, nism, nregm, nregim, idum(2)
INTEGER :: iheat, lrd, igheat, iact, icp
INTEGER :: ingb, inge
INTEGER :: i, iiso, ird, k, jiso, igg, ing
LOGICAL :: liglow

! Open a XS Library File
CLOSE(indev)
CALL openfile(indev, .TRUE., .FALSE., .FALSE., filenm)

! Card Type 1
100 FORMAT(' 1D ',10i7/(4x,10i7))
READ(indev, 100) nclpse, ggroup, ngroup, ipltfg, niso, idose, maxlrd, ncomp, nism, nregm, nregim, (idum(i), i=1,2)

! Card Type 2: group structure
ALLOCATE(en(ngroup+1), eg(ggroup+1))
200 FORMAT(' 2D ',1p5e15.8/(4x,5e15.8))
READ(indev, 200) (en(i), i=1,ngroup+1), (eg(i), i=1,ggroup+1)

! Card Type 4: isotopic names for cell average production matrix
ALLOCATE(cname(niso), mat(niso))
400 FORMAT(' 4D ',8(1x,a8)/(4x,8(1x,a8)))
410 FORMAT(' 4D ',10i7/(4x,10i7))
READ(indev, 400) (cname(i), i=1,niso)
READ(indev, 410) (mat(i), i=1,niso)

DO iiso = 1, niso
  DO jiso = 1, nelthel
   IF(cname(iiso) .EQ. ldiso(jiso)%aid) EXIT  
  END DO 
  IF(jiso .EQ. nelthel + 1)  STOP 'ReadPMATRX, Inconsistency with ISOTXS'
  lib => ldiso(jiso)
  lib%ipp = 1
  DO jiso = 1, nelthel
   IF(cname(iiso) .EQ. phatom(jiso)%aid) EXIT  
  END DO 
  IF(jiso .EQ. nelthel + 1)  STOP 'ReadPMATRX, Inconsistency with GAMISO'
  glib => phatom(jiso)
  ! Card Type 5: isotope heating
  !      IHEAT   NEUTRON HEATING AND DAMAGE FLAG          
  !              IHEAT =0, NO HEATING OR DAMAGE PRESENT FOR ISOTOPE                         
  !                    =1, HEATING AND DAMAGE PRESENT      
  !      LRD     MAXIMUM ORDER OF SCATTERING PRESENT      
  !              LRD   =0, NO SCATTERING PRESENT             
  !                    =1, ISOTROPIC SCATTERING              
  !                    =2, LINEAR ANISOTROPIC SCATTERING     
  !                    =N, ORDER N-1 ANISOTROPIC SCATTERING  
  !      IGHEAT  GAMMA HEATING FLAG                       
  !              IGHEAT=0, NO GAMMA HEATING PRESENT       
  !                    =1, GAMMA HEATING PRESENT          
  !      IACT    NUMBER OF NEUTRON REACTION CROSS SECTIONS
  !      ICP     COLLAPSING REGION NUMBER OF CURRENT ISOTOPE
  READ(indev, '(" 5D ",5i7)') iheat, lrd, igheat, iact, icp

  ! Card Type 6: neutron heating and damage
  ALLOCATE(heat(ngroup), damage(ngroup))
  600 FORMAT(' 6D ',1p5e15.8/(4x,5e15.8))
  READ(indev, 600) (heat(i), i=1,ngroup), (damage(i), i=1,ngroup)

  ALLOCATE(lib%exkappa(1:6), lib%exkappa0(1:6))
  ALLOCATE(lib%kerma_t(ngroup,1))
  ALLOCATE(lib%kerma_s(ngroup,1))
  ALLOCATE(lib%kerma_d(ngroup,1))
  ALLOCATE(lib%kerma_p(ngroup,1))
  IF(lib%ifis .EQ. 1) ALLOCATE(lib%kerma_f(ngroup,1))
  lib%kerma_t(:,1) = heat(:)
  lib%kerma_s = 0.
  IF(lib%ifis .EQ. 1) THEN 
    lib%kerma_f(:,1) = heat(:)
    lib%kerma_d = 0.
  ELSE
    lib%kerma_d(:,1) = heat(:)
    lib%kerma_f = 0.
  END IF

  ! Card Type 8: gamma heating
  ALLOCATE(gh(ggroup))
  800 FORMAT(' 8D ',1p5e15.8/(5x,5e15.8))
  READ(indev, 800) (gh(i), i=1,ggroup)

  ALLOCATE(glib%KERMA(ggroup))
  glib%KERMA(:) = gh(:)

  ! Card Type 9: cell average production matrix
  ALLOCATE(prod(ggroup,ngroup))
  900 FORMAT(' 9D ',5e15.8/(4x,5e15.8))
  DO ird = 1, lrd
    READ(indev, 900) ((prod(i,k), i=1,ggroup), k=1,ngroup)
  END DO 

  ALLOCATE(lib%ppm(0:4))
  lib%lphoton = .FALSE.
  lib%lphoton(3) = .TRUE.

  liglow = .TRUE.
  DO igg = 1, ggroup
    ingb = 0
    DO ing = 1, ngroup
      IF(prod(igg,ing) .NE. 0.) THEN
        ingb = ing
        EXIT
      END IF
    END DO 
    IF(ingb .NE. 0) THEN 
      IF(liglow) THEN 
        lib%ppm(3)%iglow = igg
        liglow = .FALSE. 
      ELSE
        lib%ppm(3)%igup = igg
      END IF
    END IF
  END DO 

  ALLOCATE(lib%ppm(3)%mat(lib%ppm(3)%iglow:lib%ppm(3)%igup,1))
  DO igg = lib%ppm(3)%iglow, lib%ppm(3)%igup
    ingb = 0
    DO ing = 1, ngroup
      IF(prod(igg,ing) .NE. 0.) THEN 
        ingb = ing
        EXIT 
      END IF
    END DO 
    DO ing = ngroup, 1, -1
      IF(prod(igg,ing) .NE. 0.) THEN
        inge = ing
        EXIT
      END IF
    END DO 
    lib%ppm(3)%mat(igg,1)%ib = ingb
    lib%ppm(3)%mat(igg,1)%ie = inge
    ALLOCATE(lib%ppm(3)%mat(igg,1)%from(ingb:inge))
    DO ing = ingb, inge
      lib%ppm(3)%mat(igg,1)%from(ing) = prod(igg,ing)
    END DO 
  END DO 

END DO 

DEALLOCATE(prod, damage, heat, gh, en, eg, mat, cname)
NULLIFY(lib, glib)

END SUBROUTINE
