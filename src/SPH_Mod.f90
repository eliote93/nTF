Module SPH_Mod
    
    implicit none
    
    TYPE SPHINTER
        INTEGER :: idiso,nndrat2u238,idx
        REAL,POINTER :: ndrat2u238(:)
    END TYPE
        
    TYPE SPHvar_type
        INTEGER :: nISOINTER,nTEMP,nU238ND,nPLTRAD,nMODXSV
        REAL(8),POINTER :: TEMP(:),U238ND(:),PLTRAD(:),MODXSV(:)
        TYPE(SPHINTER),POINTER :: ISOINTER(:)
        INTEGER :: MAPISOINTER(110000) ! EX> MAPISOINTER(94239)=1
        INTEGER :: nSubring0,nSPHreg_srd,nrestfxr,nSPHreg
        INTEGER,ALLOCATABLE :: SPHreg_srd(:)
    END TYPE
    
    ! Direct Interpolation
    TYPE RC_type
        REAL,POINTER :: rc(:)
    END TYPE
    TYPE DIpos_type
        TYPE(RC_type),POINTER :: DIpos(:)
    END TYPE
    TYPE spDIpos_type
        TYPE(DIpos_type),POINTER :: spDI(:)
    END TYPE
    
    TYPE(SPHvar_type),POINTER :: SPHvar(:)
    TYPE SPHdata_type
        REAL(4),POINTER :: sphf(:)
    END TYPE
    TYPE(SPHdata_type),POINTER :: SPHdata(:)
    TYPE idxmaptype
        INTEGER,POINTER :: idx(:,:,:,:,:,:)
        INTEGER,POINTER :: spidx(:,:,:,:)
    END TYPE
    TYPE (idxmaptype),POINTER :: idxmap(:)
    
    CHARACTER :: FILE_SPH*80
    LOGICAL :: lssphndavg
    REAL(4),POINTER :: ssphf(:,:,:), ssphfnm(:,:,:)   !--- CNJ Edit : Node Majors
        
     ! Direct Interpolation
    TYPE (spDIpos_type),POINTER :: DI(:)

    INTEGER :: nspf
    
    CONTAINS
    
    SUBROUTINE findSRDidx(svr,nSPHreg_interp,srdidx,geomchflag)
        IMPLICIT NONE
        
        INTEGER,INTENT(IN) :: nSPHreg_interp
        INTEGER,INTENT(OUT) :: srdidx
        LOGICAL,INTENT(OUT) :: geomchflag
        TYPE(SPHvar_type),POINTER :: svr
        INTEGER :: isrd,nr
        LOGICAL :: flag
        geomchflag=.FALSE.
        IF (nSPHreg_interp.ge.svr%SPHreg_srd(0)) THEN
            srdidx=0
            IF (nSPHreg_interp.gt.svr%SPHreg_srd(0)) geomchflag=.TRUE.
            RETURN
        ENDIF
        IF (svr%nSPHreg_srd.eq.0) THEN
            srdidx=0
            RETURN
        ENDIF
        flag=.false.
        DO isrd=1,svr%nSPHreg_srd
            nr=svr%SPHreg_srd(isrd)
            IF (nr.le.nSPHreg_interp) THEN
                IF (nr.eq.nSPHreg_interp) THEN
                    srdidx=isrd
                ELSEIF (nr.lt.nSPHreg_interp) THEN
                    srdidx=isrd-1
                ENDIF
                flag=.true.
                EXIT
            ENDIF
        ENDDO    
        IF (.not.flag) THEN
            srdidx=svr%nSPHreg_srd
        ENDIF
        
    END SUBROUTINE

    SUBROUTINE calcSPH(ispf,pltrad,modxsv,u238nd,temp,idiso,pnum,niso,nFxr,nSPHreg_interp,srdidx,ig,sphf_localfxr)
        USE TYPEDEF,  ONLY : Cell_Type
        IMPLICIT NONE
        
        INTEGER,INTENT(IN) :: niso,nFxr,srdidx,ig,nSPHreg_interp,ispf
        INTEGER,INTENT(IN) :: idiso(niso)
        REAL,INTENT(IN) :: pltrad,modxsv,u238nd,temp,pnum(niso)
        REAL(4),INTENT(OUT) :: sphf_localfxr(nFxr)
        TYPE(SPHvar_type),POINTER :: svr
        TYPE(Cell_Type),POINTER :: CellInfo(:)
        
        INTEGER :: idxm(2),idxu(2),idxp(2),idxt(2),idxi(2)
        REAL(4) :: sphfm_interp(25,8),sphfu_interp(25,4),sphfp_interp(25,2)
        REAL(4) :: sphcorfm_interp(25,16),sphcorfu_interp(25,8),sphcorfp_interp(25,4)
        REAL(4) :: sphcorft_interp(25,2),sphcorfi_interp(25)
        REAL(4) :: basesphf(25),corsum(25)
        INTEGER :: k,m,it,ip,iu,im,j,ii,iiso
        INTEGER :: ninterfiso,jdiso,iin,nSPHreg_in
        REAL(4) :: ndratiso
        
        svr=>SPHvar(ispf)
        basesphf=1._4
        CALL findidxSPHP(ispf,pltrad,idxp)
        CALL findidxSPHM(ispf,modxsv,idxm)
        if (ispf.eq.0) CALL findidxSPHU(ispf,u238nd,idxu)
        CALL findidxSPHT(ispf,temp,idxt)
        nSPHreg_in=svr%SPHreg_srd(srdidx)
        if (ispf.eq.0) then
            sphfm_interp=1._4; sphfu_interp=1._4; sphfp_interp=1._4
            k=1
            DO it=1,2; DO ip=1,2; DO iu=1,2
                CALL interpBASESPHM(ispf,nSPHreg_interp,nSPHreg_in,srdidx,ig,idxt(it),idxp(ip),idxu(iu),idxm,modxsv,sphfm_interp(1:nSPHreg_interp,k))
                k=k+1
            ENDDO; ENDDO; ENDDO
            k=1
            DO it=1,2; DO ip=1,2
                m=2*k-1
                CALL interpSPHU(ispf,nSPHreg_interp,ig,idxu,u238nd,sphfu_interp(1:nSPHreg_interp,k),sphfm_interp(1:nSPHreg_interp,m:m+1)) 
                k=k+1
            ENDDO; ENDDO
            k=1
            DO it=1,2
                m=2*k-1
                CALL interpSPHP(ispf,nSPHreg_interp,ig,idxp,pltrad,sphfp_interp(1:nSPHreg_interp,k),sphfu_interp(1:nSPHreg_interp,m:m+1)) 
                k=k+1
            ENDDO
        else
            sphfm_interp=1._4; sphfp_interp=1._4
            k=1
            DO it=1,2; DO ip=1,2
                CALL interpBASEspSPHM(ispf,nSPHreg_interp,nSPHreg_in,srdidx,ig,idxt(it),idxp(ip),idxm,modxsv,sphfm_interp(1:nSPHreg_interp,k))
                k=k+1
            ENDDO; ENDDO
            k=1
            DO it=1,2
                m=2*k-1
                CALL interpSPHP(ispf,nSPHreg_interp,ig,idxp,pltrad,sphfp_interp(1:nSPHreg_interp,k),sphfm_interp(1:nSPHreg_interp,m:m+1)) 
                k=k+1
            ENDDO
        endif
        CALL interpSPHT(ispf,nSPHreg_interp,ig,idxt,temp,basesphf(1:nSPHreg_interp),sphfp_interp(1:nSPHreg_interp,1:2))

        if (ispf.eq.0) then
            corsum=0._4
            DO iin=1,niso
                jdiso=idiso(iin)
                if (jdiso.eq.92238) cycle
                if (jdiso.eq.0) cycle
                ndratiso=pnum(iin)/U238ND
                iiso=svr%MAPISOINTER(jdiso)
                IF (iiso.eq.0) CYCLE
                CALL findSPHInterfIdx(ispf,ndratiso,iiso,idxi)
                sphcorfm_interp=1._4; sphcorfu_interp=1._4; sphcorfp_interp=1._4; sphcorft_interp=1._4; sphcorfi_interp=1._4
                k=1
                DO ii=1,2; DO it=1,2; DO ip=1,2; DO iu=1,2;
                    CALL interpCORSPHM(ispf,nSPHreg_interp,nSPHreg_in,srdidx,iiso,ig,idxi(ii),idxt(it),idxp(ip),idxu(iu),idxm,modxsv,sphcorfm_interp(1:nSPHreg_interp,k)) 
                    k=k+1
                ENDDO; ENDDO; ENDDO; ENDDO
                k=1
                DO ii=1,2; DO it=1,2; DO ip=1,2
                    m=2*k-1
                    CALL interpSPHU(ispf,nSPHreg_interp,ig,idxu,u238nd,sphcorfu_interp(1:nSPHreg_interp,k),sphcorfm_interp(1:nSPHreg_interp,m:m+1)) 
                    k=k+1
                ENDDO; ENDDO; ENDDO
                k=1
                DO ii=1,2; DO it=1,2
                    m=2*k-1
                    CALL interpSPHP(ispf,nSPHreg_interp,ig,idxp,pltrad,sphcorfp_interp(1:nSPHreg_interp,k),sphcorfu_interp(1:nSPHreg_interp,m:m+1)) 
                    k=k+1
                ENDDO; ENDDO
                k=1
                DO ii=1,2
                    m=2*k-1
                    CALL interpSPHT(ispf,nSPHreg_interp,ig,idxt,temp,sphcorft_interp(1:nSPHreg_interp,k),sphcorfp_interp(1:nSPHreg_interp,m:m+1))
                    k=k+1
                ENDDO
                CALL interpSPHI(ispf,nSPHreg_interp,iiso,ig,idxi,ndratiso,sphcorfi_interp(1:nSPHreg_interp),sphcorft_interp(1:nSPHreg_interp,1:2))
                corsum=corsum+(sphcorfi_interp-1._4)
            ENDDO
            corsum=corsum+1._4
        
            sphf_localfxr=1._4
            DO j=1,nSPHreg_interp
                sphf_localfxr(nFxr-j+1)=corsum(j)*basesphf(j)
            ENDDO
        else
            sphf_localfxr=1._4
            DO j=1,nSPHreg_interp
                sphf_localfxr(nFxr-j+1)=basesphf(j)
            ENDDO
        endif
        
    END SUBROUTINE
    
    SUBROUTINE findidxSPHM(ispf,modxsv_c,idxm) ! moderator xsv interpolation index
        IMPLICIT NONE        
        INTEGER,INTENT(IN) :: ispf
        REAL,INTENT(IN) :: modxsv_c
        INTEGER,INTENT(OUT) :: idxm(2)
        TYPE(SPHvar_type),POINTER :: svr
        INTEGER :: nm,im
        svr=>SPHvar(ispf)
        nm=svr%nMODXSV
        DO im=1,nm
            IF (modxsv_c.lt.svr%MODXSV(im)) EXIT
        ENDDO
        IF (im.eq.1) THEN
            idxm=1
        ELSEIF (im.gt.nm) THEN
            idxm=nm
        ELSE
            idxm(1)=im-1
            idxm(2)=im
        ENDIF
        nullify(svr)
    END SUBROUTINE
    SUBROUTINE findidxSPHU(ispf,u238nd_c,idxu) ! U238 ND Index
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: ispf
        REAL,INTENT(IN) :: u238nd_c
        INTEGER,INTENT(OUT) :: idxu(2)
        TYPE(SPHvar_type),POINTER :: svr
        INTEGER :: nu,iu
        svr=>SPHvar(ispf)
        nu=svr%nU238ND
        DO iu=1,nu
            IF (u238nd_c.lt.svr%U238ND(iu)) EXIT
        ENDDO
        IF (iu.eq.1) THEN
            idxu=1
        ELSEIF (iu.gt.nu) THEN
            idxu=nu
        ELSE
            idxu(1)=iu-1
            idxu(2)=iu
        ENDIF
        nullify(svr)
    END SUBROUTINE   
    SUBROUTINE findidxSPHP(ispf,pltrad_c,idxp) ! Fuel Pellet Radius
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: ispf
        REAL,INTENT(IN) :: pltrad_c
        INTEGER,INTENT(OUT) :: idxp(2)
        TYPE(SPHvar_type),POINTER :: svr
        INTEGER :: np,ip
        svr=>SPHvar(ispf)
        np=svr%nPLTRAD
        DO ip=1,np
            IF (pltrad_c.lt.svr%PLTRAD(ip)) EXIT
        ENDDO
        IF (ip.eq.1) THEN
            idxp=1
        ELSEIF (ip.gt.np) THEN
            idxp=np            
        ELSE
            idxp(1)=ip-1
            idxp(2)=ip
        ENDIF
        nullify(svr)
    END SUBROUTINE       
    SUBROUTINE findidxSPHT(ispf,temp_c,idxt) ! Temperature interpolation...
        IMPLICIT NONE 
        INTEGER,INTENT(IN) :: ispf
        REAL,INTENT(IN) :: temp_c
        INTEGER,INTENT(OUT) :: idxt(2)
        TYPE(SPHvar_type),POINTER :: svr
        INTEGER :: nt,it
        svr=>SPHvar(ispf)
        nt=svr%nTEMP
        DO it=1,nt
            IF (temp_c.lt.svr%TEMP(it)) EXIT
        ENDDO
        IF (it.eq.1) THEN
            idxt=1
        ELSEIF (it.gt.nt) THEN
            idxt=nt
        ELSE
            idxt(1)=it-1
            idxt(2)=it
        ENDIF
        nullify(svr)
    END SUBROUTINE 
    
    SUBROUTINE interpBASESPHM(ispf,nSPHreg_interp,nSPHreg_in,srdidx,ig,idxt,idxp,idxu,idxm,modxsv_c,sphfm_interp)
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: ispf,nSPHreg_interp,nSPHreg_in,srdidx,ig,idxt,idxp,idxu,idxm(2)
        REAL,INTENT(IN) :: modxsv_c
        REAL(4),INTENT(OUT) :: sphfm_interp(nSPHreg_interp)
        TYPE(SPHvar_type),POINTER :: svr
        REAL(4),POINTER :: sphf(:)
        INTEGER,POINTER :: idmp(:,:,:,:,:,:)
        REAL(4) :: sph1(nSPHreg_in),sph2(nSPHreg_in),sphfm1_interp(nSPHreg_interp),sphfm2_interp(nSPHreg_interp)
        REAL :: m1,m2
        INTEGER :: idx_start,srd_start
        svr=>SPHvar(ispf)
        sphf=>SPHdata(ispf)%sphf
        idmp=>idxmap(ispf)%idx
        srd_start=sum(svr%SPHreg_srd(0:srdidx-1))
        IF (idxm(1).ne.idxm(2)) THEN
            m1=svr%MODXSV(idxm(1))
            m2=svr%MODXSV(idxm(2))
            idx_start=idmp(ig,idxm(1),idxu,idxp,idxt,0)
            sph1=sphf(idx_start+srd_start:idx_start+srd_start+svr%SPHreg_srd(srdidx)-1)
            CALL InterpSPH(ispf,sph1,sphfm1_interp,nSPHreg_in,nSPHreg_interp,idxp)
            idx_start=idmp(ig,idxm(2),idxu,idxp,idxt,0)
            sph2=sphf(idx_start+srd_start:idx_start+srd_start+svr%SPHreg_srd(srdidx)-1)
            CALL InterpSPH(ispf,sph2,sphfm2_interp,nSPHreg_in,nSPHreg_interp,idxp)
            sphfm_interp=(sphfm2_interp-sphfm1_interp)/(m2-m1)*(modxsv_c-m1)+sphfm1_interp
        ELSE
            idx_start=idmp(ig,idxm(1),idxu,idxp,idxt,0)
            sph1=sphf(idx_start+srd_start:idx_start+srd_start+svr%SPHreg_srd(srdidx)-1)
            CALL InterpSPH(ispf,sph1,sphfm_interp,nSPHreg_in,nSPHreg_interp,idxp)
        ENDIF
        nullify(svr,sphf,idmp)
    END SUBROUTINE
    SUBROUTINE interpBASEspSPHM(ispf,nSPHreg_interp,nSPHreg_in,srdidx,ig,idxt,idxp,idxm,modxsv_c,sphfm_interp)
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: ispf,nSPHreg_interp,nSPHreg_in,srdidx,ig,idxt,idxp,idxm(2)
        REAL,INTENT(IN) :: modxsv_c
        REAL(4),INTENT(OUT) :: sphfm_interp(nSPHreg_interp)
        TYPE(SPHvar_type),POINTER :: svr
        REAL(4),POINTER :: sphf(:)
        INTEGER,POINTER :: idmp(:,:,:,:)
        REAL(4) :: sph1(nSPHreg_in),sph2(nSPHreg_in),sphfm1_interp(nSPHreg_interp),sphfm2_interp(nSPHreg_interp)
        REAL :: m1,m2
        INTEGER :: idx_start,srd_start
        svr=>SPHvar(ispf)
        sphf=>SPHdata(ispf)%sphf
        idmp=>idxmap(ispf)%spidx
        srd_start=sum(svr%SPHreg_srd(0:srdidx-1))
        IF (idxm(1).ne.idxm(2)) THEN
            m1=svr%MODXSV(idxm(1))
            m2=svr%MODXSV(idxm(2))
            idx_start=idmp(ig,idxm(1),idxp,idxt)
            sph1=sphf(idx_start+srd_start:idx_start+srd_start+svr%SPHreg_srd(srdidx)-1)
            CALL InterpSPH(ispf,sph1,sphfm1_interp,nSPHreg_in,nSPHreg_interp,idxp)
            idx_start=idmp(ig,idxm(2),idxp,idxt)
            sph2=sphf(idx_start+srd_start:idx_start+srd_start+svr%SPHreg_srd(srdidx)-1)
            CALL InterpSPH(ispf,sph2,sphfm2_interp,nSPHreg_in,nSPHreg_interp,idxp)
            sphfm_interp=(sphfm2_interp-sphfm1_interp)/(m2-m1)*(modxsv_c-m1)+sphfm1_interp
        ELSE
            idx_start=idmp(ig,idxm(1),idxp,idxt)
            sph1=sphf(idx_start+srd_start:idx_start+srd_start+svr%SPHreg_srd(srdidx)-1)
            CALL InterpSPH(ispf,sph1,sphfm_interp,nSPHreg_in,nSPHreg_interp,idxp)
        ENDIF
        nullify(svr,sphf,idmp)
    END SUBROUTINE
    SUBROUTINE interpSPHU(ispf,nSPHreg_interp,ig,idxu,u238nd_c,sphfu_interp,sphfm_interp)
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: ispf,nSPHreg_interp,ig,idxu(2) 
        REAL(4),INTENT(IN) :: sphfm_interp(nSPHreg_interp,2)
        REAL,INTENT(IN) :: u238nd_c
        REAL(4),INTENT(OUT) :: sphfu_interp(nSPHreg_interp)
        TYPE(SPHvar_type),POINTER :: svr
        REAL(4) :: sph1(nSPHreg_interp),sph2(nSPHreg_interp)
        REAL :: u1,u2
        svr=>SPHvar(ispf)
        IF (idxu(1).ne.idxu(2)) THEN
            u1=svr%U238ND(idxu(1))
            u2=svr%U238ND(idxu(2))
            sph1=sphfm_interp(:,1)
            sph2=sphfm_interp(:,2)
            sphfu_interp=(sph2-sph1)/(u2-u1)*(u238nd_c-u1)+sph1
        ELSE
            sphfu_interp=sphfm_interp(:,1)
        ENDIF
        nullify(svr)
    END SUBROUTINE
    SUBROUTINE interpSPHP(ispf,nSPHreg_interp,ig,idxp,pltrad_c,sphfp_interp,sphfu_interp) 
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: ispf,nSPHreg_interp,ig,idxp(2)
        REAL(4),INTENT(IN) :: sphfu_interp(nSPHreg_interp,2)
        REAL,INTENT(IN) :: pltrad_c
        REAL(4),INTENT(OUT) :: sphfp_interp(nSPHreg_interp)
        TYPE(SPHvar_type),POINTER :: svr
        REAL(4) :: sph1(nSPHreg_interp),sph2(nSPHreg_interp)
        REAL :: p1,p2
        svr=>SPHvar(ispf)
        IF (idxp(1).ne.idxp(2)) THEN
            p1=svr%PLTRAD(idxp(1))
            p2=svr%PLTRAD(idxp(2))
            sph1=sphfu_interp(:,1)
            sph2=sphfu_interp(:,2)
            sphfp_interp=(sph2-sph1)/(p2-p1)*(pltrad_c-p1)+sph1
        ELSE
            sphfp_interp=sphfu_interp(:,1)
        ENDIF
        nullify(svr)
    END SUBROUTINE      
    SUBROUTINE interpSPHT(ispf,nSPHreg_interp,ig,idxt,temp_c,sphft_interp,sphfp_interp)
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: ispf,nSPHreg_interp,ig,idxt(2)
        REAL,INTENT(IN) :: temp_c
        REAL(4),INTENT(IN) :: sphfp_interp(nSPHreg_interp,2)
        REAL(4),INTENT(OUT) :: sphft_interp(nSPHreg_interp)
        TYPE(SPHvar_type),POINTER :: svr
        REAL(4) :: sph1(nSPHreg_interp),sph2(nSPHreg_interp)
        REAL :: t1,t2
        svr=>SPHvar(ispf)
        IF (idxt(1).ne.idxt(2)) THEN
            t1=sqrt(svr%TEMP(idxt(1)))
            t2=sqrt(svr%TEMP(idxt(2)))
            sph1=sphfp_interp(:,1)
            sph2=sphfp_interp(:,2)
            sphft_interp=(sph2-sph1)/(t2-t1)*(dsqrt(temp_c)-t1)+sph1
        ELSE
            sphft_interp=sphfp_interp(:,1)
        ENDIF
        nullify(svr)
    END SUBROUTINE
    
    SUBROUTINE InterpSPH(ispf,sphin,sphout,nSPHregIN,nSPHregOUT,idxp)
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: ispf,nSPHregIN,nSPHregOUT,idxp
        REAL(4),INTENT(IN) :: sphin(nSPHregIN)
        REAL(4),INTENT(OUT) :: sphout(nSPHregOUT)
        TYPE(SPHvar_type),POINTER :: svr
        TYPE(RC_type),POINTER :: DIpos(:)
        INTEGER :: n1
        svr=>SPHvar(ispf)
        DIpos=>DI(ispf)%spDI(idxp)%DIpos
        IF (nSPHregIN.eq.nSPHregOUT) THEN
            sphout=sphin
            RETURN
        ENDIF
        n1=svr%nrestfxr
        sphout(nSPHregOUT-n1+1:nSPHregOUT)=sphin(nSPHregIN-n1+1:nSPHregIN)
        CALL DI_SPH(DIpos(nSPHregIN-n1)%rc(1:nSPHregIN-n1),sphin,DIpos(nSPHregOUT-n1)%rc(1:nSPHregOUT-n1),sphout(1:nSPHregOUT-n1),nSPHregIN-n1,nSPHregOUT-n1)
        nullify(svr,DIpos)
    END SUBROUTINE
    
    SUBROUTINE findSPHInterfIdx(ispf,ndratiso_c,iiso,idxi)
        IMPLICIT NONE
        REAL(4),INTENT(IN) :: ndratiso_c
        INTEGER,INTENT(IN) :: ispf,iiso
        INTEGER,INTENT(OUT) :: idxi(2)
        TYPE(SPHvar_type),POINTER :: svr
        TYPE(SPHINTER),POINTER :: sphisointer
        INTEGER :: nndrat,ii
        svr=>SPHvar(ispf)
        sphisointer=>svr%ISOINTER(iiso)
        nndrat=sphisointer%nndrat2u238
        DO ii=1,nndrat
            IF (ndratiso_c.lt.sphisointer%ndrat2u238(ii)) EXIT
        ENDDO
        IF (ii.gt.nndrat) THEN
            idxi=sphisointer%idx+nndrat-1
        ELSE
            idxi(1)=sphisointer%idx+ii-2
            idxi(2)=sphisointer%idx+ii-1
        ENDIF
        nullify(svr)
    END SUBROUTINE
    
    SUBROUTINE interpCORSPHM(ispf,nSPHreg_interp,nSPHreg_in,srdidx,iiso,ig,idxi,idxt,idxp,idxu,idxm,modxsv_c,corsphfm_interp)
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: ispf,nSPHreg_interp,nSPHreg_in,srdidx,iiso,ig,idxi,idxt,idxp,idxu,idxm(2)
        REAL,INTENT(IN) :: modxsv_c
        REAL(4),INTENT(OUT) :: corsphfm_interp(nSPHreg_interp)
        TYPE(SPHvar_type),POINTER :: svr
        REAL(4),POINTER :: sphf(:)
        INTEGER,POINTER :: idmp(:,:,:,:,:,:)
        REAL(4) :: sph1(nSPHreg_in),sph2(nSPHreg_in),cor1(nSPHreg_in),cor2(nSPHreg_in)
        REAL :: m1,m2
        REAL(4) :: sphfm1_interp(nSPHreg_interp),sphfm2_interp(nSPHreg_interp),corsphfm1_interp(nSPHreg_interp),corsphfm2_interp(nSPHreg_interp)
        INTEGER :: idx_start,srd_start
        svr=>SPHvar(ispf)
        idmp=>idxmap(ispf)%idx
        sphf=>SPHdata(ispf)%sphf
        m1=svr%MODXSV(idxm(1))
        m2=svr%MODXSV(idxm(2))
        srd_start=sum(svr%SPHreg_srd(0:srdidx-1))
        IF (idxi.ge.svr%ISOINTER(iiso)%idx) THEN
            IF (idxm(1).ne.idxm(2)) THEN    
                idx_start=idmp(ig,idxm(1),idxu,idxp,idxt,idxi)
                cor1=sphf(idx_start+srd_start:idx_start+srd_start+svr%SPHreg_srd(srdidx)-1)
                idx_start=idmp(ig,idxm(1),idxu,idxp,idxt,0)
                sph1=sphf(idx_start+srd_start:idx_start+srd_start+svr%SPHreg_srd(srdidx)-1)
                CALL InterpSPH(ispf,cor1*sph1,corsphfm1_interp,nSPHreg_in,nSPHreg_interp,idxp)
                CALL InterpSPH(ispf,sph1,sphfm1_interp,nSPHreg_in,nSPHreg_interp,idxp)
                corsphfm1_interp=corsphfm1_interp/sphfm1_interp
                idx_start=idmp(ig,idxm(2),idxu,idxp,idxt,idxi)
                cor2=sphf(idx_start+srd_start:idx_start+srd_start+svr%SPHreg_srd(srdidx)-1)
                idx_start=idmp(ig,idxm(2),idxu,idxp,idxt,0)
                sph2=sphf(idx_start+srd_start:idx_start+srd_start+svr%SPHreg_srd(srdidx)-1)
                CALL InterpSPH(ispf,cor2*sph2,corsphfm2_interp,nSPHreg_in,nSPHreg_interp,idxp)
                CALL InterpSPH(ispf,sph2,sphfm2_interp,nSPHreg_in,nSPHreg_interp,idxp)
                corsphfm2_interp=corsphfm2_interp/sphfm2_interp
                corsphfm_interp=(corsphfm2_interp-corsphfm1_interp)/(m2-m1)*(modxsv_c-m1)+corsphfm1_interp
            ELSE
                idx_start=idmp(ig,idxm(1),idxu,idxp,idxt,idxi)
                cor1=sphf(idx_start+srd_start:idx_start+srd_start+svr%SPHreg_srd(srdidx)-1)
                idx_start=idmp(ig,idxm(1),idxu,idxp,idxt,0)
                sph1=sphf(idx_start+srd_start:idx_start+srd_start+svr%SPHreg_srd(srdidx)-1)
                CALL InterpSPH(ispf,cor1*sph1,corsphfm1_interp,nSPHreg_in,nSPHreg_interp,idxp)
                CALL InterpSPH(ispf,sph1,sphfm1_interp,nSPHreg_in,nSPHreg_interp,idxp)
                corsphfm_interp=corsphfm1_interp/sphfm1_interp
            ENDIF
        ELSEIF (idxi.lt.svr%ISOINTER(iiso)%idx) THEN
            corsphfm_interp=1._4
        ENDIF
        nullify(svr,sphf,idmp)
    END SUBROUTINE
    SUBROUTINE interpSPHI(ispf,nSPHreg_interp,iiso,ig,idxi,ndratiso_c,sphfi_interp,sphft_interp)
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: ispf,nSPHreg_interp,iiso,ig,idxi(2)
        REAL(4),INTENT(IN) :: sphft_interp(nSPHreg_interp,2)
        REAL(4),INTENT(IN) :: ndratiso_c
        REAL(4),INTENT(OUT) :: sphfi_interp(nSPHreg_interp)
        TYPE(SPHvar_type),POINTER :: svr
        REAL(4) :: sph1(nSPHreg_interp),sph2(nSPHreg_interp)
        REAL :: i1,i2
        TYPE(SPHINTER),POINTER :: sphisointer
        svr=>SPHvar(ispf)
        sphisointer=>svr%ISOINTER(iiso)
        IF (idxi(1).ne.idxi(2)) THEN
            IF (idxi(1).lt.sphisointer%idx) THEN
                i1=0._4
                sph1=1._4
            ELSE
                i1=sphisointer%ndrat2u238(idxi(1)-sphisointer%idx+1)
                sph1=sphft_interp(:,1)
            ENDIF
            i2=sphisointer%ndrat2u238(idxi(2)-sphisointer%idx+1)
            sph2=sphft_interp(:,2)
            sphfi_interp=(sph2-sph1)/(i2-i1)*(ndratiso_c-i1)+sph1
        ELSE
            sphfi_interp=sphft_interp(:,1)
        ENDIF
        nullify(svr,sphisointer)
    END SUBROUTINE
    
    SUBROUTINE DI_SPH(xin,yin,xout,yout,nin,nout)
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: nin,nout
        REAL,INTENT(IN) :: xin(nin),xout(nout)
        REAL(4),INTENT(IN) :: yin(nin)
        REAL(4),INTENT(OUT) :: yout(nout)
        INTEGER :: i,j
        REAL :: x1,x2,y1,y2
        DO i=1,nout
            DO j=1,nin
                IF (xin(j).gt.xout(i)) EXIT
            ENDDO
            IF (j.eq.nin+1) THEN
                yout(i)=yin(nin)
            ELSE
                x1=xin(j-1); x2=xin(j);
                y1=yin(j-1); y2=yin(j);
                yout(i)=(y2-y1)/(x2-x1)*(xout(i)-x1)+y1          
            ENDIF
        ENDDO
    END SUBROUTINE

! Subroutine to make SPH factors in all the FXR regions in fuel cell with hole in the center of pellet, which is
! specialized for air(the region doesn't contain U238)
    SUBROUTINE calcCellSSPH(CellInfo,icel,ndata,iReg,rr,ndiv,ibfuel,iefuel,nfueldiv,lhole) 
    ! EDIT JSU EDIT 20190711 for Hole problem
! CAUTION-------------------------------------------------------------------------------------------------------*
!   With this routine, the spectral SPH method can be applied to air-hole containing problems, however,         *
!   the ring-tickness of each subregions for SPH library are generated with the assumption that the subrings    *
!   are the results of making equi-volumic rings of one fuel region. Therefore, if one try to use more than one *
!   air hole regions in input's cell information or use more than one fuel regions, it would deteriorate the    *
!   result due to unaccurate SPH factor.                                                                        *
!---------------------------------------------------------------------------------------------------------------*
    !** Cell information
    ! iReg       : (INT(:)) material index of regions in cell
    ! rr         : (DBL(:)) radii of regions
    ! ndiv       : (INT(:)) equivolumetric division number of each cell
    ! iefuel     : (INT)    region index of outtermost fuel region
    ! ibfuel     : (INT)    region index of innermost fuel region
    ! nfueldiv   : (INT)    number of fuel division...
    ! lhole      : (logic)  for existence of air hole
        USE param ,ONLY : sigpH, sigpO, sigpB10, sigpB11
        USE TYPEDEF ,ONLY : Cell_Type
        USE Material_Mod, ONLY : Mixture
        USE XSLIB_MOD ,ONLY : igresb,igrese,nofghel
        IMPLICIT NONE
        INTEGER,intent(in) :: icel,ndata,iReg(0:300),ibfuel,iefuel
        INTEGER,intent(inout) :: ndiv(0:300)
        REAL, intent(in) :: rr(300)
        LOGICAL, INTENT(IN) :: lhole
        TYPE(Cell_Type), POINTER :: CellInfo(:)
        INTEGER,PARAMETER :: ispf=0
        INTEGER :: i,j,k,idiso_c_avg(300), niso, l, idiso, imat
        INTEGER :: nSPHreg_interp,nSPHreg0,ig,nfueldiv,ngapdiv,ncladdiv,srdidx,srdmapidx,nfxrdiff,ifxr
        REAL :: vol,temptemp,Temp,sphpnum_avg,pltrad_c,modxsv_c,U238ND,u238nd_c,temp_c,pnum_c_avg(300)  ! BYS EDIT 100-> 300
        REAL(4) :: sphf_interp(40)
        LOGICAL :: lgeomchflag
        TYPE(SPHvar_type),POINTER :: svr
        
        IF (.not.CellInfo(icel)%lfuel) return ! Only Fuel Pin (Do not consider surrounding effect)
        svr=>SPHvar(ispf)
!  1. nuclides(iso) IDs and average number density in fuel region...
        idiso_c_avg=0; niso=1
        idiso_c_avg(1)=Mixture(ireg(ibfuel))%idiso(1)
        DO j=ibfuel,iefuel,-1
            do l=1,Mixture(ireg(j))%niso
                do k=1,niso
                    if (idiso_c_avg(k).eq.Mixture(ireg(j))%idiso(l)) exit
                enddo
                if (k.gt.niso) then
                    niso=niso+1
                    idiso_c_avg(niso)=Mixture(ireg(j))%idiso(l)
                endif
            enddo
        ENDDO
        pnum_c_avg=0._8
        do k=1,niso
            vol=0._8
            DO j=ibfuel,iefuel,-1
                vol=vol+CellInfo(icel)%vol(j)
                DO l=1,Mixture(ireg(j))%niso
                    IF (idiso_c_avg(k).eq.Mixture(ireg(j))%idiso(l)) exit
                ENDDO
                if (l.gt.Mixture(ireg(j))%niso) cycle
                pnum_c_avg(k)=pnum_c_avg(k)+Mixture(ireg(j))%pnum(l)*CellInfo(icel)%vol(j)
            ENDDO
            pnum_c_avg(k)=pnum_c_avg(k)/vol
        enddo
!  2. average number density of spectralSPH factor nuclides(U238 and nuclides indicated with #)
        sphpnum_avg=0._8; temptemp=0._8; vol=0._8
        DO j=ibfuel,iefuel,-1
            vol=vol+CellInfo(icel)%vol(j)
            DO l=1,Mixture(ireg(j))%niso
                IF (Mixture(ireg(j))%lSSPH(l)) exit ! Material contains U238 or isotope starts with #
            ENDDO
            if (l.gt.Mixture(ireg(j))%niso) then
                write(*,'(a10,i4,a21)') 'Fuel Cell ',icel,' does not have U238.'
                cycle
            endif
            sphpnum_avg=sphpnum_avg+Mixture(ireg(j))%pnum(l)*CellInfo(icel)%vol(j)
            temptemp=temptemp+Mixture(ireg(j))%temp*CellInfo(icel)%vol(j)
        ENDDO
        sphpnum_avg=sphpnum_avg/vol
        temptemp=temptemp/vol      ! average temperature [K]
        CellInfo(icel)%U238ND0 = sphpnum_avg
        CellInfo(icel)%FuelRefTemp0 = temptemp
!  3. MODXSV (Macroscopic XS of potential XS and volume
        DO j=ndata,0,-1
            IF (.not.Mixture(ireg(j))%lh2o) CYCLE
            DO k=1,Mixture(ireg(j))%niso
                idiso=Mixture(ireg(j))%idiso(K)
                IF (idiso.eq.1001) THEN 
                    CellInfo(icel)%MODXSV0=CellInfo(icel)%MODXSV0+Mixture(ireg(j))%pnum(k)*sigpH*CellInfo(icel)%vol(j)
                ELSEIF (idiso.eq.8016) THEN
                    CellInfo(icel)%MODXSV0=CellInfo(icel)%MODXSV0+Mixture(ireg(j))%pnum(k)*sigpO*CellInfo(icel)%vol(j) 
                ELSEIF (idiso.eq.5010) THEN
                    CellInfo(icel)%MODXSV0=CellInfo(icel)%MODXSV0+Mixture(ireg(j))%pnum(k)*sigpB10*CellInfo(icel)%vol(j) 
                ELSEIF (idiso.eq.5011) THEN
                    CellInfo(icel)%MODXSV0=CellInfo(icel)%MODXSV0+Mixture(ireg(j))%pnum(k)*sigpB11*CellInfo(icel)%vol(j) 
                ENDIF
            ENDDO
        ENDDO 
!  4. Region info for non-fuel materials surrounding pellet (gap and cladding)
        ngapdiv=0
        DO j=iefuel-1,0,-1
            IF (Mixture(ireg(j))%lcld) EXIT
            ngapdiv=ngapdiv+ndiv(j)
        ENDDO
        if (j.eq.-1) then
            ngapdiv=0
            ncladdiv=0
        else
            ncladdiv=ndiv(j)
        endif
        IF (lhole) nfueldiv = nfueldiv  + SUM(ndiv(ibfuel+1:ndata)) ! for air gap...
        nSPHreg_interp=nfueldiv+svr%nrestfxr
!  5. Fine sSPH fuel radius index whose #of radii is closest larger or equal data to given problem(nTRACER input).
        CALL findSRDidx(svr,nSPHreg_interp,srdidx,lgeomchflag)
        IF (lgeomchflag) THEN
            !IF (fuelregidx.ne.ndata) STOP '# of fuel sub-divisions is greater than the maximum division of SPH library in procgeom.f90!'
          IF (.NOT. lhole) THEN
            ndiv(ndata)=svr%SPHreg_srd(srdidx)-svr%nrestfxr
            nfxrdiff=ndiv(ndata)-nfueldiv
            nfueldiv=ndiv(ndata)
            nSPHreg_interp=nSPHreg_interp+nfxrdiff
            CellInfo(icel)%geom%nCircle=CellInfo(icel)%geom%nCircle+nfxrdiff
            CellInfo(icel)%nFXR=CellInfo(icel)%nFXR+nfxrdiff
            CellInfo(icel)%nFSR=CellInfo(icel)%nFSR+nfxrdiff*CellInfo(icel)%nDivAzi  
          ELSE
            !STOP '# of fuel sub-divisions is greater than the maximum division of SPH library in SPH_Mod.f90!' ! remaining task for air hole 20190710
            ndiv(ndata-1)=svr%SPHreg_srd(srdidx)-svr%nrestfxr - SUM(ndiv(ibfuel+1:ndata))
            nfxrdiff=ndiv(ndata-1) - nfueldiv + SUM(ndiv(ibfuel+1:ndata))
            nfueldiv=ndiv(ndata-1) + SUM(ndiv(ibfuel+1:ndata))
            nSPHreg_interp=nSPHreg_interp+nfxrdiff
            CellInfo(icel)%geom%nCircle=CellInfo(icel)%geom%nCircle + nfxrdiff
            CellInfo(icel)%nFXR=CellInfo(icel)%nFXR + nfxrdiff
            CellInfo(icel)%nFSR=CellInfo(icel)%nFSR + nfxrdiff*CellInfo(icel)%nDivAzi  
          END IF
        ENDIF
        CellInfo(icel)%nfueldiv=nfueldiv
        CellInfo(icel)%ngapdiv=ngapdiv
        CellInfo(icel)%ncladdiv=ncladdiv
        CellInfo(icel)%srdidx=srdidx 
        ! calculate default Cell SPH factor
        ALLOCATE(CellInfo(icel)%SPHfactor(1:CellInfo(icel)%nFXR,igresb:igrese))
        pltrad_c=CellInfo(icel)%FuelRad0
        modxsv_c=CellInfo(icel)%MODXSV0
        u238nd_c=CellInfo(icel)%U238ND0
        temp_c=CellInfo(icel)%FuelRefTemp0
        DO ig=igresb,igrese
            CellInfo(icel)%SPHfactor(:,ig)=1._4
            ifxr = 0
            DO j=ndata, ibfuel+1, -1 ! HOLE REGION meaningless if ndata.EQ.ibfuel
                imat = ireg(j)
                DO k = 1, ndiv(j)
                    ifxr = ifxr + 1
                ENDDO
            END DO
            DO j=ibfuel,iefuel,-1
                imat = ireg(j)
                DO l=1,Mixture(imat)%niso
                    IF (Mixture(imat)%lSSPH(l)) exit
                ENDDO
                U238ND = Mixture(imat)%pnum(l)
                Temp = Mixture(imat)%temp
                DO k = 1, ndiv(j)
                    ifxr = ifxr + 1
                    i = CellInfo(icel)%nFXR-ifxr+1 ! from inside to outside
                    CALL calcSPH(ispf,pltrad_c,modxsv_c,U238ND,Temp,Mixture(imat)%idiso(1:niso),Mixture(imat)%pnum(1:niso),&
                         Mixture(imat)%niso,CellInfo(icel)%nFXR,nSPHreg_interp,srdidx,ig-nofghel,sphf_interp(1:CellInfo(icel)%nFXR))
                    CellInfo(icel)%SPHfactor(i,ig)=sphf_interp(i)
                ENDDO
            ENDDO
            CALL calcSPH(ispf,pltrad_c,modxsv_c,u238nd_c,temp_c,idiso_c_avg,pnum_c_avg,niso,CellInfo(icel)%nFXR,nSPHreg_interp,srdidx,ig-nofghel,sphf_interp(1:CellInfo(icel)%nFXR))
            IF (svr%nrestfxr.ge.1) THEN
                DO ifxr=nfueldiv+1,nfueldiv+ngapdiv
                    CellInfo(icel)%SPHfactor(CellInfo(icel)%nFXR-ifxr+1,ig)=sphf_interp(CellInfo(icel)%nFXR-nfueldiv)
                ENDDO
                IF (svr%nrestfxr.eq.2) THEN
                    DO ifxr=nfueldiv+ngapdiv+1,nfueldiv+ngapdiv+ncladdiv
                        CellInfo(icel)%SPHfactor(CellInfo(icel)%nFXR-ifxr+1,ig)=sphf_interp(CellInfo(icel)%nFXR-nfueldiv-1)
                    ENDDO
                ENDIF
            ENDIF
        ENDDO
        nullify(svr)
    END SUBROUTINE 
    SUBROUTINE calcAICCellSSPH(CellInfo,icel,ndata,iReg,rr,ndiv,fuelregidx,nfueldiv)
        USE param ,ONLY : sigpH, sigpO, sigpB10, sigpB11
        USE TYPEDEF ,ONLY : Cell_Type
        USE Material_Mod, ONLY : Mixture
        USE XSLIB_MOD ,ONLY : igresb,igrese,nofghel
        IMPLICIT NONE
        INTEGER,intent(in) :: icel,ndata,iReg(0:300),fuelregidx
        INTEGER,intent(inout) :: ndiv(0:300)
        REAL, intent(in) :: rr(300)
        TYPE(Cell_Type), POINTER :: CellInfo(:)
        INTEGER,PARAMETER :: ispf=1
        INTEGER :: i,j,k, niso, l, idiso, imat, idiso_c_avg(300)
        INTEGER :: nSPHreg_interp,nSPHreg0,ig,nfueldiv,ngapdiv,ncladdiv,srdidx,srdmapidx,nfxrdiff,ifxr
        REAL :: vol,Temp,pltrad_c,modxsv_c,temp_c,pnum_c_avg(300)
        REAL(4) :: sphf_interp(40)
        LOGICAL :: geomchflag
        TYPE(SPHvar_type),POINTER :: svr
        
        svr=>SPHvar(ispf)
        CellInfo(icel)%U238ND0 = 0._8
        CellInfo(icel)%FuelRefTemp0 = Mixture(ireg(fuelregidx))%temp
        niso = Mixture(ireg(fuelregidx))%niso
        DO j=ndata,0,-1
            IF (.not.Mixture(ireg(j))%lh2o) CYCLE
            DO k=1,Mixture(ireg(j))%niso
                idiso=Mixture(ireg(j))%idiso(K)
                IF (idiso.eq.1001) THEN 
                    CellInfo(icel)%MODXSV0=CellInfo(icel)%MODXSV0+Mixture(ireg(j))%pnum(k)*sigpH*CellInfo(icel)%vol(j)
                ELSEIF (idiso.eq.8016) THEN
                    CellInfo(icel)%MODXSV0=CellInfo(icel)%MODXSV0+Mixture(ireg(j))%pnum(k)*sigpO*CellInfo(icel)%vol(j) 
                ELSEIF (idiso.eq.5010) THEN
                    CellInfo(icel)%MODXSV0=CellInfo(icel)%MODXSV0+Mixture(ireg(j))%pnum(k)*sigpB10*CellInfo(icel)%vol(j) 
                ELSEIF (idiso.eq.5011) THEN
                    CellInfo(icel)%MODXSV0=CellInfo(icel)%MODXSV0+Mixture(ireg(j))%pnum(k)*sigpB11*CellInfo(icel)%vol(j) 
                ENDIF
            ENDDO
        ENDDO 
        !nfueldiv=0
        !DO j=ndata,fuelregidx,-1
        !    nfueldiv=nfueldiv+ndiv(j)
        !ENDDO
        DO j=fuelregidx-1,fuelregidx-2,-1
            IF (Mixture(ireg(j))%lh2o) EXIT
        ENDDO
        if (j.eq.fuelregidx-1) then
            ngapdiv=0
            ncladdiv=0
        elseif (j.eq.fuelregidx-2) then
            ngapdiv=0
            ncladdiv=1
        else
            ngapdiv=1
            ncladdiv=1
        endif
        nSPHreg_interp=nfueldiv+svr%nrestfxr
        CALL findSRDidx(svr,nSPHreg_interp,srdidx,geomchflag)
        IF (geomchflag) THEN
            !IF (fuelregidx.ne.ndata) STOP '# of AIC sub-divisions is greater than the maximum division of SPH library in procgeom.f90!'
            ndiv(ndata)=svr%SPHreg_srd(srdidx)-svr%nrestfxr
            nfxrdiff=ndiv(ndata)-nfueldiv
            nfueldiv=ndiv(ndata)
            nSPHreg_interp=nSPHreg_interp+nfxrdiff
            CellInfo(icel)%geom%nCircle=CellInfo(icel)%geom%nCircle+nfxrdiff
            CellInfo(icel)%nFXR=CellInfo(icel)%nFXR+nfxrdiff
            CellInfo(icel)%nFSR=CellInfo(icel)%nFSR+nfxrdiff*CellInfo(icel)%nDivAzi      
        ENDIF
        CellInfo(icel)%nfueldiv=nfueldiv
        CellInfo(icel)%ngapdiv=ngapdiv
        CellInfo(icel)%ncladdiv=ncladdiv
        CellInfo(icel)%srdidx=srdidx 
        ! calculate default Cell SPH factor
        ALLOCATE(CellInfo(icel)%SPHfactor(1:CellInfo(icel)%nFXR,igresb:igrese))
        pltrad_c=CellInfo(icel)%FuelRad0
        modxsv_c=CellInfo(icel)%MODXSV0
        temp_c=CellInfo(icel)%FuelRefTemp0
        DO ig=igresb,igrese
            CellInfo(icel)%SPHfactor(:,ig)=1._4
            CALL calcSPH(ispf,pltrad_c,modxsv_c,1._8,temp_c,idiso_c_avg,pnum_c_avg,niso,CellInfo(icel)%nFXR,nSPHreg_interp,srdidx,ig-nofghel,sphf_interp(1:CellInfo(icel)%nFXR))
            CellInfo(icel)%SPHfactor(CellInfo(icel)%nFXR-nfueldiv+1:CellInfo(icel)%nFXR,ig)=sphf_interp(CellInfo(icel)%nFXR-nfueldiv+1:CellInfo(icel)%nFXR)
            IF (svr%nrestfxr.ge.1) THEN
                DO ifxr=nfueldiv+1,nfueldiv+ngapdiv
                    CellInfo(icel)%SPHfactor(CellInfo(icel)%nFXR-ifxr+1,ig)=sphf_interp(CellInfo(icel)%nFXR-nfueldiv)
                ENDDO
                IF (svr%nrestfxr.eq.2) THEN
                    DO ifxr=nfueldiv+ngapdiv+1,nfueldiv+ngapdiv+ncladdiv
                        CellInfo(icel)%SPHfactor(CellInfo(icel)%nFXR-ifxr+1,ig)=sphf_interp(CellInfo(icel)%nFXR-nfueldiv-1)
                    ENDDO
                ENDIF
            ENDIF
        ENDDO
        nullify(svr)
    END SUBROUTINE 
    
    ! sph calculation w.r.t. the ND change in each FXR
    SUBROUTINE calcPinSSPH(Core, Fxr, PE)
        USE param ,ONLY : sigpH, sigpO, sigpB10, sigpB11
        USE TYPEDEF ,ONLY : CoreInfo_Type, FxrInfo_Type, PE_Type, Cell_Type, Pin_Type
        USE XSLIB_MOD ,ONLY : igresb,igrese,nofghel
        IMPLICIT NONE
        
        TYPE(CoreInfo_Type) :: Core
        TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
        TYPE(PE_Type) :: PE
        !Pointing Variables
        TYPE(FxrInfo_Type), POINTER :: myFxr
        TYPE(Cell_Type), POINTER :: CellInfo(:)
        TYPE(Pin_Type), POINTER :: Pin(:)
        INTEGER :: myzb, myze, nxy,ispf
        INTEGER :: nlocalFxr, FxrIdxSt
        INTEGER :: ipin ,icel, iz, ifsr, ifxr, ig, niso, icel0
        INTEGER :: j, k, l, idiso, idiso_c_avg(1000), nSPHreg_interp, nFsrinfxr
        REAL :: MODXSV,PLTRAD,TEMP_avg,U238ND_avg,volsum,U238ND,pnum_c_avg(1000),vol
        REAL(4),ALLOCATABLE :: sphf_localfxr(:)
        TYPE(SPHvar_type),POINTER :: svr
        
        myzb = PE%myzb; myze = PE%myze
        Pin => Core%Pin
        CellInfo => Core%CellInfo
        nxy = Core%nxy
        DO iz=myzb, myze
            IF(.NOT. (Core%lFuelPlane(iz) .OR. Core%lAICPlane(iz))) CYCLE
            ispf=0; svr=>SPHvar(ispf)
            DO ipin = 1, nxy
                FxrIdxSt = Pin(ipin)%FxrIdxSt
                icel = Pin(ipin)%Cell(iz)
                icel0 = CellInfo(icel)%icel0
                IF (.not.CellInfo(icel)%lfuel) CYCLE
                nlocalFxr = CellInfo(icel)%nFxr 
                ALLOCATE(sphf_localfxr(nlocalFxr))
                nSPHreg_interp=CellInfo(icel)%nfueldiv+svr%nrestfxr
                ! calculate MODXSV
                MODXSV=0._8
                DO j= 1, nLocalFxr
                    ifxr = FxrIdxSt + j -1;  myFxr => Fxr(ifxr, iz)
                    ifsr=CellInfo(icel0)%MapFxr2FsrIdx(1,j)
                    vol=CellInfo(icel0)%vol(ifsr)
                    nfsrinfxr=CellInfo(icel0)%nFsrInFxr(j)
                    vol=vol*nfsrinfxr
                    IF (myFxr%lh2o) THEN
                        DO k=1,myFxr%niso
                            idiso=myFxr%idiso(k)
                            IF (idiso.eq.1001) THEN 
                                MODXSV=MODXSV+myFxr%pnum(k)*sigpH*vol
                            ELSEIF (idiso.eq.8016) THEN
                                MODXSV=MODXSV+myFxr%pnum(k)*sigpO*vol
                            ELSEIF (idiso.eq.5010) THEN
                                MODXSV=MODXSV+myFxr%pnum(k)*sigpB10*vol
                            ELSEIF (idiso.eq.5011) THEN
                                MODXSV=MODXSV+myFxr%pnum(k)*sigpB11*vol 
                            ENDIF
                        ENDDO
                    ENDIF
                ENDDO !DO j= 1, nLocalFxr
                PLTRAD=CellInfo(icel)%FuelRad0
                ! calculate U238 average number density in a fuel pellet
                ! collect all isotopes in a fuel pellet
                ! calculate average temperature in a fuel pellet
                ! calculate region-wise spectral sph factors in a fuel pellet using region-wise parameters
                U238ND_avg=0._8
                TEMP_avg=0._8
                volsum=0._8
                idiso_c_avg=0; niso=1
                idiso_c_avg(1)=Fxr(FxrIdxSt + nLocalFxr - 1, iz)%idiso(1)
                DO j = 1, nLocalFxr
                    ifxr = FxrIdxSt + j - 1;  myFxr => Fxr(ifxr, iz)
                    IF (.not.myFxr%lSSPHalloc) THEN
                        ALLOCATE(myFxr%SPHfactor(igresb:igrese))
                        myFxr%lSSPHalloc=.true.
                    ENDIF
                    myFxr%SPHfactor=1._4
                    IF (.not.myFxr%lfuel) CYCLE
                    U238ND=0._8
                    ifsr=CellInfo(icel)%MapFxr2FsrIdx(1,j)
                    volsum=volsum+CellInfo(icel)%vol(ifsr)
                    DO k=1,myFxr%niso
                        IF (myFxr%idiso(k).eq.92238) U238ND=myFxr%pnum(k) 
                        do l=1,niso
                            if (idiso_c_avg(l).eq.myFxr%idiso(k)) exit
                        enddo
                        if (l.gt.niso) then
                            niso=niso+1
                            idiso_c_avg(niso)=myFxr%idiso(k)
                        endif
                    ENDDO    
                    IF (U238ND.eq.0) CYCLE
                    U238ND_avg=U238ND_avg+U238ND*CellInfo(icel)%vol(ifsr)
                    TEMP_avg=TEMP_avg+myFxr%Temp*CellInfo(icel)%vol(ifsr)
                    DO ig=igresb,igrese
                        CALL calcSPH(ispf,PLTRAD,MODXSV,U238ND,myFxr%Temp,myFxr%idiso(1:niso),myFxr%pnum(1:niso),niso,nlocalFxr,nSPHreg_interp,CellInfo(icel)%srdidx,ig-nofghel,sphf_localfxr)
                        myFxr%SPHfactor(ig)=sphf_localfxr(j)
                    ENDDO
                ENDDO !DO j= 1, nLocalFxr
                U238ND_avg=U238ND_avg/volsum
                TEMP_avg=TEMP_avg/volsum
                ! calculate isotope-wise average number density in a fuel pellet
                pnum_c_avg=0._8
                do l=1,niso
                    volsum=0._8
                    DO j = 1, nLocalFxr
                        ifxr = FxrIdxSt + j - 1;  myFxr => Fxr(ifxr, iz)
                        IF (.not.myFxr%lfuel) CYCLE
                        ifsr=CellInfo(icel)%MapFxr2FsrIdx(1,j)
                        DO k=1,myFxr%niso
                            IF (idiso_c_avg(l).eq.myFxr%idiso(k)) exit
                        ENDDO
                        if (k.gt.myFxr%niso) cycle
                        pnum_c_avg(l)=pnum_c_avg(l)+myFxr%pnum(k)*CellInfo(icel)%vol(ifsr)
                        volsum=volsum+CellInfo(icel)%vol(ifsr)
                    ENDDO
                    pnum_c_avg(l)=pnum_c_avg(l)/volsum
                enddo
                IF (svr%nrestfxr.ge.1) THEN
                    DO k=1,CellInfo(icel)%ngapdiv
                        j=CellInfo(icel)%nfueldiv+k
                        ifxr = FxrIdxSt + nLocalFxr - j;  myFxr => Fxr(ifxr, iz)
                        DO ig=igresb,igrese
                            CALL calcSPH(ispf,PLTRAD,MODXSV,U238ND_avg,TEMP_avg,idiso_c_avg(1:niso),pnum_c_avg(1:niso),niso,nlocalFxr,nSPHreg_interp,CellInfo(icel)%srdidx,ig-nofghel,sphf_localfxr)
                            myFxr%SPHfactor(ig)=sphf_localfxr(nLocalFxr - CellInfo(icel)%nfueldiv)
                        ENDDO
                    ENDDO
                ENDIF
                IF (svr%nrestfxr.eq.2) THEN
                    DO k=1,CellInfo(icel)%ncladdiv
                        j=CellInfo(icel)%nfueldiv+CellInfo(icel)%ngapdiv+k
                        ifxr = FxrIdxSt + nLocalFxr - j;  myFxr => Fxr(ifxr, iz)
                        DO ig=igresb,igrese
                            CALL calcSPH(ispf,PLTRAD,MODXSV,U238ND_avg,TEMP_avg,idiso_c_avg(1:niso),pnum_c_avg(1:niso),niso,nlocalFxr,nSPHreg_interp,CellInfo(icel)%srdidx,ig-nofghel,sphf_localfxr)
                            myFxr%SPHfactor(ig)=sphf_localfxr(nLocalFxr - CellInfo(icel)%nfueldiv - 1)
                        ENDDO
                    ENDDO
                ENDIF
                DEALLOCATE(sphf_localfxr)
            ENDDO ! DO ipin = 1, nxy
            ! AIC
            IF(.NOT. Core%lAICPlane(iz)) CYCLE
            ispf=1; svr=>SPHvar(ispf)
            DO ipin = 1, nxy
                FxrIdxSt = Pin(ipin)%FxrIdxSt
                icel = Pin(ipin)%Cell(iz)
                icel0 = CellInfo(icel)%icel0
                IF (.not.CellInfo(icel)%lAIC) CYCLE
                nlocalFxr = CellInfo(icel)%nFxr 
                ALLOCATE(sphf_localfxr(nlocalFxr))
                nSPHreg_interp=CellInfo(icel)%nfueldiv+svr%nrestfxr
                ! calculate MODXSV
                MODXSV=0._8
                DO j= 1, nLocalFxr
                    ifxr = FxrIdxSt + j -1;  myFxr => Fxr(ifxr, iz)
                    ifsr=CellInfo(icel0)%MapFxr2FsrIdx(1,j)
                    vol=CellInfo(icel0)%vol(ifsr)
                    nfsrinfxr=CellInfo(icel0)%nFsrInFxr(j)
                    vol=vol*nfsrinfxr
                    IF (myFxr%lh2o) THEN
                        DO k=1,myFxr%niso
                            idiso=myFxr%idiso(k)
                            IF (idiso.eq.1001) THEN 
                                MODXSV=MODXSV+myFxr%pnum(k)*sigpH*vol
                            ELSEIF (idiso.eq.8016) THEN
                                MODXSV=MODXSV+myFxr%pnum(k)*sigpO*vol
                            ELSEIF (idiso.eq.5010) THEN
                                MODXSV=MODXSV+myFxr%pnum(k)*sigpB10*vol
                            ELSEIF (idiso.eq.5011) THEN
                                MODXSV=MODXSV+myFxr%pnum(k)*sigpB11*vol 
                            ENDIF
                        ENDDO
                    ENDIF
                ENDDO !DO j= 1, nLocalFxr
                PLTRAD=CellInfo(icel)%FuelRad0
                ! calculate average temperature in a fuel pellet
                ! calculate region-wise spectral sph factors in a fuel pellet using region-wise parameters
                TEMP_avg=0._8; volsum=0._8
                DO j = 1, nLocalFxr
                    ifxr = FxrIdxSt + j - 1;  myFxr => Fxr(ifxr, iz)
                    IF (.not.myFxr%lSSPHalloc) THEN
                        ALLOCATE(myFxr%SPHfactor(igresb:igrese))
                        myFxr%lSSPHalloc=.true.
                    ENDIF
                    myFxr%SPHfactor=1._4
                    IF (.not.myFxr%lAIC) CYCLE
                    ifsr=CellInfo(icel)%MapFxr2FsrIdx(1,j)
                    volsum=volsum+CellInfo(icel)%vol(ifsr)
                    TEMP_avg=TEMP_avg+myFxr%Temp*CellInfo(icel)%vol(ifsr)
                ENDDO !DO j= 1, nLocalFxr
                TEMP_avg=TEMP_avg/volsum
                DO ig=igresb,igrese
                    CALL calcSPH(ispf,PLTRAD,MODXSV,1._8,TEMP_avg,idiso_c_avg(1:niso),pnum_c_avg(1:niso),niso,nlocalFxr,nSPHreg_interp,CellInfo(icel)%srdidx,ig-nofghel,sphf_localfxr)
                    DO j = 1, nLocalFxr
                        ifxr = FxrIdxSt + j - 1;  myFxr => Fxr(ifxr, iz)
                        if (myFxr%lh2o) cycle
                        myFxr%SPHfactor(ig)=sphf_localfxr(j)
                    ENDDO
                ENDDO
                DEALLOCATE(sphf_localfxr)
            ENDDO ! DO ipin = 1, nxy
        ENDDO ! DO iz=myzb, myze
        NULLIFY(Pin,CellInfo,myFxr,svr)
    END SUBROUTINE
    
END Module