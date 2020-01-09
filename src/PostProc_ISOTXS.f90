SUBROUTINE PostProc_ISOTXS
    USE XSLIB_MOD,      ONLY : ldiso, nelthel, noghel, np1hel, idnp1hel, mapnuclp13, nfishel, mapfis,&
                               itempmap, itempmapp1, nxtempmap
!    USE MultiCat_Mod,   ONLY : resoset, nresoset
    USE NuclidMap_mod,  ONLY : lExistIsotope
    USE ALLOCS
    IMPLICIT NONE
    INTEGER :: iso,jso,ig,jg,idiso,it,it1,it2,i
    INTEGER :: ng,ntemp
    REAL(4) :: numden
    
    ng = noghel
    
    DO iso = 1, nelthel
        IF(ldiso(iso)%nid /= 0) CYCLE
        idiso = nint(ldiso(iso)%aw)
        idiso = idiso*10
        DO WHILE(.TRUE.)
            IF(.NOT. lExistIsotope(idiso)) EXIT
            idiso = idiso + 1
        ENDDO
        ldiso(iso)%nid = idiso
    ENDDO
    
    !List of isotopes contained P1 scattering XS
    np1hel = 0
    DO iso = 1, nelthel
        IF(ldiso(iso)%nid == 0) CYCLE
        IF(ldiso(iso)%np1temp > 0) np1hel = np1hel + 1
    ENDDO
    ALLOCATE(idnp1hel(np1hel))
    CALL dmalloc(idnp1hel,np1hel)
    jso = 0
    DO iso = 1, nelthel
        IF(ldiso(iso)%np1temp == 0) CYCLE
        IF(ldiso(iso)%nid == 0) CYCLE
        jso = jso + 1
        idnp1hel(jso) = iso
        idiso = ldiso(iso)%nid
        mapnuclp13(idiso) = jso
    ENDDO
    
    !List of isotopes contained Fission XS
    nfishel = 0
    DO iso = 1, nelthel
        IF(ldiso(iso)%nid == 0) CYCLE
        IF(ldiso(iso)%ityp /= 3) CYCLE
        nfishel = nfishel + 1
        mapfis(ldiso(iso)%nid) = nfishel
        !ldiso(iso)%kappa = 1.e-6
    ENDDO
    
    !Temperatrue Map for P0 Scattering
    CALL dmalloc(itempmap,nxtempmap,nelthel)
    CALL dmalloc(itempmapp1,nxtempmap,np1hel)
    DO iso = 1, nelthel
        ntemp = ldiso(iso)%ntemp
        it2 = 0
        IF(it1 > 1) THEN
            it2 = NINT(ldiso(iso)%temp(1))
            DO it = 1, it2
                itempmap(it,iso) = 1
            ENDDO
            DO i = 2, ntemp
                it1 = it2 + 1
                it2 = NINT(ldiso(iso)%temp(i))
                DO it = it1, it2
                    itempmap(it,iso) = i - 1
                ENDDO
            ENDDO
            it1 = it2 + 1
            it2 = nxtempmap
            DO it = it1, it2
                itempmap(it,iso) = ntemp - 1
            ENDDO
        ELSE
            DO it = 1, nxtempmap
                itempmap(it,iso) = 1
            ENDDO
        ENDIF
    ENDDO
    
    !Temperatrue Map for P1 Scattering
    DO iso = 1, np1hel
        jso = idnp1hel(iso)
        ntemp = ldiso(jso)%np1temp
        it2 = 0
        IF(it1 > 1) THEN
            it2 = NINT(ldiso(iso)%p1temp(1))
            DO it = 1, it2
                itempmapp1(it,iso) = 1
            ENDDO
            DO i = 2, ntemp
                it1 = it2 + 1
                it2 = NINT(ldiso(iso)%p1temp(i))
                DO it = it1, it2
                    itempmapp1(it,iso) = i - 1
                ENDDO
            ENDDO
            it1 = it2 + 1
            it2 = nxtempmap
            DO it = it1, it2
                itempmapp1(it,iso) = ntemp - 1
            ENDDO
        ELSE
            DO it = 1, nxtempmap
                itempmapp1(it,iso) = 1
            ENDDO
        ENDIF
    ENDDO
    
    !Out Scattering Range
    DO iso = 1, nelthel
        DO it = 1, ldiso(iso)%ntemp
            DO ig = 1, noghel
                DO jg = 1, noghel
                    IF(ig >= ldiso(iso)%sm(jg,it)%ib .and. ig <= ldiso(iso)%sm(jg,it)%ie) THEN
                        ldiso(iso)%sm(ig,it)%ioutsb = jg
                        EXIT
                    END IF
                ENDDO
                DO jg = noghel, 1, -1
                    IF(ig >= ldiso(iso)%sm(jg,it)%ib .and. ig <= ldiso(iso)%sm(jg,it)%ie) THEN
                        ldiso(iso)%sm(ig,it)%ioutse = jg
                        EXIT
                    END IF
                ENDDO
            ENDDO
        ENDDO
        
        ALLOCATE(ldiso(iso)%sigstr(noghel,ldiso(iso)%ntemp))
        DO it = 1, ldiso(iso)%ntemp
            DO ig = 1, noghel
                ldiso(iso)%sigstr(ig,it) = 0.
            ENDDO
            DO ig = 1, noghel
                DO jg = ldiso(iso)%sm(ig,it)%ib,ldiso(iso)%sm(ig,it)%ie
                    ldiso(iso)%sigstr(jg,it) = ldiso(iso)%sigstr(jg,it) + ldiso(iso)%sm(ig,it)%from(jg)
                ENDDO
            ENDDO
        ENDDO
    ENDDO
    
!   allocate(resoset(nresoset))
    
    RETURN
END SUBROUTINE PostProc_ISOTXS