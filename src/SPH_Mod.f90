  Module SPH_Mod
    use param,        ONLY :  FALSE
    USE Timer,        ONLY : nTracer_dclock, TimeChk
    IMPLICIT NONE
    TYPE SPHINTER
        INTEGER :: idiso,nndrat2u238,idx
        REAL,POINTER :: ndrat2u238(:)
    END TYPE
    TYPE SPHvar_type
        INTEGER :: nISOINTER,nTEMP,nU238ND,nPLTRAD,nMODXSV, nINTER
        REAL(8),POINTER :: TEMP(:),U238ND(:),PLTRAD(:),MODXSV(:)
        TYPE(SPHINTER),POINTER :: ISOINTER(:)
        INTEGER :: MAPISOINTER(110000) ! EX> MAPISOINTER(94239)=1
        INTEGER :: nSubring0,nSPHreg_srd,nrestfxr,nSPHreg
        INTEGER,ALLOCATABLE :: SPHreg_srd(:)
    END TYPE
    ! Geometrical SPH table
    INTEGER :: nfuelSPH, nAICSPH
    INTEGER, POINTER :: idxSPH2BaseCell(:), idxAICSPH2BaseCell(:)
    TYPE FUEL_SPH_type
        INTEGER,POINTER :: idx(:,:,:,:,:) ! for fuel   | ig,MODXSV,U238ND,TEMP,ISO)
        INTEGER,POINTER :: idx_NM(:,:,:,:) ! for fuel, fxr major | MODXSV,U238ND,TEMP,ISO | ordered by ifxr-ig
        REAL(4),POINTER :: sphf(:)
        REAL(4),POINTER :: sphf_NM(:)
    END TYPE
    TYPE (FUEL_SPH_type),POINTER :: FUEL_SPH(:)
    TYPE Special_SPH_TYPE
        INTEGER,POINTER :: idx(:,:) ! for special SPH | ig,MODXSV
        INTEGER,POINTER :: idx_NM(:)! for special SPH fxr major | MODXSV
        REAL(4),POINTER :: sphf(:)
        REAL(4),POINTER :: sphf_NM(:)
    END TYPE
    TYPE (Special_SPH_TYPE),POINTER :: SPECIAL_SPH(:)
    ! Direct Interpolation
    TYPE RC_type
        REAL,POINTER :: rc(:) ! radial center position
    END TYPE
    TYPE DIpos_type
        TYPE(RC_type),POINTER :: DIpos(:) ! by total subring number
    END TYPE
    TYPE spDIpos_type
        TYPE(DIpos_type),POINTER :: spDI(:) ! pellet radius array
    END TYPE
    TYPE (spDIpos_type),POINTER :: DI(:) ! DATA for each SPH library
    ! DI(0)%spDI(1)%DIpos(10)%rc(7) : center position of 7th ring when fuel is divided for 10 subrings(from innermost)  &
    !                                 with the radius of the 1st pellet radius data
    TYPE(SPHvar_type),POINTER :: SPHvar(:)
    TYPE SPHdata_type
        REAL(4),POINTER :: sphf(:)
    END TYPE
    TYPE(SPHdata_type),POINTER :: SPHdata(:)
    TYPE idxmaptype
        INTEGER,POINTER :: idx(:,:,:,:,:,:) ! ig,MODXSV,U238ND,PLTRAD,TEMP,ISO)
        INTEGER,POINTER :: spidx(:,:,:,:)
    END TYPE
    TYPE (idxmaptype),POINTER :: idxmap(:)
    
    CHARACTER :: FILE_SPH*80
    LOGICAL :: lssphndavg
    REAL(4),POINTER :: ssphf(:,:,:), ssphfnm(:,:,:)   !--- CNJ Edit : Node Majors

    INTEGER :: nspf

  CONTAINS
!
SUBROUTINE SetSSPHBasicCellInfo(CellInfo,ndata,iReg,ndiv,ibfuel,iefuel,nfueldiv,lhole,lAIC)
  USE param ,ONLY : sigpH, sigpO, sigpB10, sigpB11
  USE TYPEDEF ,ONLY : Cell_Type
  USE Material_Mod, ONLY : Mixture
  USE XSLIB_MOD ,ONLY : igresb,igrese,nofghel
  IMPLICIT NONE
  INTEGER,intent(in) :: ndata,iReg(0:300),ibfuel,iefuel
  INTEGER,intent(inout) :: ndiv(0:300)
  LOGICAL, INTENT(IN) :: lhole, lAIC
  TYPE(Cell_Type) :: CellInfo
  INTEGER :: i,j,k
  INTEGER :: nFXR_interpol,nfueldiv,ngapdiv,ncladdiv,srdidx,nfxrdiff
  LOGICAL :: lgeomchflag
  TYPE(SPHvar_type),POINTER :: svr
  REAL :: Tbeg2_SSPH, Tend2_ssph
  !
  IF (lAIC) THEN
      CellInfo%FuelRefTEMP0 = Mixture(ireg(iefuel))%temp
      DO j=iefuel-1,iefuel-2,-1
          IF (Mixture(ireg(j))%lh2o) EXIT
      ENDDO
      if (j.eq.iefuel-1) then
          ngapdiv=0;  ncladdiv=0
      elseif (j.eq.iefuel-2) then
          ngapdiv=0;  ncladdiv=1
      else
          ngapdiv=1;  ncladdiv=1
      endif
      svr => sphvar(1)
  ELSE !fuel
      IF (CellInfo%lhole) nfueldiv = nfueldiv + SUM(ndiv(ibfuel+1:ndata))
      ngapdiv=0
      DO j=iefuel-1,0,-1
        IF (Mixture(ireg(j))%lcld) EXIT
        ngapdiv = ngapdiv + ndiv(j)
      END DO
      if (j.eq.-1) then
        ngapdiv=0; ncladdiv=0
      else
        ncladdiv=ndiv(j)
      endif
      svr => sphvar(0)
  ENDIF
  nFXR_interpol = nfueldiv+svr%nrestfxr
  CALL findSRDidx(svr, nFXR_interpol, srdidx, lgeomchflag)
  IF (lgeomchflag) THEN
      IF (lAIC) THEN
          ndiv(ndata)=svr%SPHreg_srd(srdidx)-svr%nrestfxr
          nfxrdiff=ndiv(ndata)-nfueldiv
          nfueldiv=ndiv(ndata)
          CellInfo%geom%nCircle=CellInfo%geom%nCircle+nfxrdiff
          CellInfo%nFXR=CellInfo%nFXR+nfxrdiff
          CellInfo%nFSR=CellInfo%nFSR+nfxrdiff*CellInfo%nDivAzi      
      ELSE !fuel
          IF (.NOT. CellInfo%lhole) THEN
              ndiv(ndata)=svr%SPHreg_srd(srdidx)-svr%nrestfxr
              nfxrdiff=ndiv(ndata)-nfueldiv !# of fuel FXRs to be reduced to use PSSL (sign will be negative)
              nfueldiv=ndiv(ndata)          !# of fuel FXRs to be used
          ELSE
              ndiv(ndata-1)=svr%SPHreg_srd(srdidx)-svr%nrestfxr - SUM(ndiv(ibfuel+1:ndata))
              nfxrdiff=ndiv(ndata-1) - nfueldiv + SUM(ndiv(ibfuel+1:ndata))
              nfueldiv=ndiv(ndata-1) + SUM(ndiv(ibfuel+1:ndata))
          END IF
          CellInfo%geom%nCircle=CellInfo%geom%nCircle + nfxrdiff
          CellInfo%nFXR=CellInfo%nFXR + nfxrdiff
          CellInfo%nFSR=CellInfo%nFSR + nfxrdiff*CellInfo%nDivAzi  
        ENDIF
  ENDIF
  CellInfo%nfueldiv = nfueldiv
  CellInfo%ngapdiv  = ngapdiv
  CellInfo%ncladdiv = ncladdiv
  CellInfo%srdidx   = srdidx 
END SUBROUTINE SetSSPHBasicCellInfo
!
SUBROUTINE findSRDidx(svr,nSPHreg_interp,srdidx,geomchflag)
  ! Subroutine to find SRD index
  ! INPUT  -----------------------------------------------------------------------------
  ! svr            : Spectral SPH variable data
  ! nSPHreg_interp :# of FXRs used for the interpolation of SPH factor, nfueldiv+nrest
  !                  nrest is determined by spectral SPH library (gap, cladding)
  ! OUTPUT -----------------------------------------------------------------------------
  ! srdidx         : Radial division index 
  ! geomchflag     : Flag if (# of FXRs > Max.# of FXRs in svr)geometry change is required for PSSL
  ! REMARK -----------------------------------------------------------------------------
  !  Only called at the CalcCellSSPH & CalcAICCellSSPH (initialization step)
  !
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
        RETURN !# of FXRs is too large to use PSSL
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
        srdidx=svr%nSPHreg_srd !# of FXRs used for PSSL is less than Min.# of FXRs in svr
    ENDIF
END SUBROUTINE
!
SUBROUTINE findidxSPHM(ispf,modxsv_c,idxm, wgtm) ! moderator xsv interpolation index
  IMPLICIT NONE        
  INTEGER,INTENT(IN) :: ispf
  REAL,INTENT(IN) :: modxsv_c
  INTEGER,INTENT(OUT) :: idxm(2)
  REAL, INTENT(OUT) :: wgtm(2)
  TYPE(SPHvar_type),POINTER :: svr
  INTEGER :: nm,im
  svr=>SPHvar(ispf)
  nm=svr%nMODXSV
  DO im=1,nm
      IF (modxsv_c.lt.svr%MODXSV(im)) EXIT
  ENDDO
  IF (im.eq.1) THEN
      idxm=1
      wgtm(1) = 1._8
      wgtm(2) = 0._8
  ELSEIF (im.gt.nm) THEN
      idxm=nm
      wgtm(1) = 1._8
      wgtm(2) = 0._8
  ELSE
      idxm(1)=im-1
      idxm(2)=im
      wgtm(1)=(svr%MODXSV(im)-modxsv_c)/(svr%MODXSV(im)-svr%MODXSV(im-1))
      wgtm(2)=(modxsv_c-svr%MODXSV(im-1))/(svr%MODXSV(im)-svr%MODXSV(im-1))
  ENDIF
  nullify(svr)
END SUBROUTINE
!
SUBROUTINE findidxSPHU(ispf,u238nd_c,idxu, wgtu) ! U238 ND Index
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ispf
    REAL,INTENT(IN) :: u238nd_c
    INTEGER,INTENT(OUT) :: idxu(2)
    REAL, INTENT(OUT) :: wgtu(2)
    TYPE(SPHvar_type),POINTER :: svr
    INTEGER :: nu,iu
    svr=>SPHvar(ispf)
    nu=svr%nU238ND
    DO iu=1,nu
        IF (u238nd_c.lt.svr%U238ND(iu)) EXIT
    ENDDO
    IF (iu.eq.1) THEN
        idxu=1
        wgtu(1) = 1._8
        wgtu(2) = 0._8
    ELSEIF (iu.gt.nu) THEN
        idxu=nu
        wgtu(1) = 1._8
        wgtu(2) = 0._8
    ELSE
        idxu(1)=iu-1
        idxu(2)=iu
        wgtu(1)=(svr%U238ND(iu)-u238nd_c)/(svr%U238ND(iu)-svr%U238ND(iu-1))
        wgtu(2)=(u238nd_c-svr%U238ND(iu-1))/(svr%U238ND(iu)-svr%U238ND(iu-1))
    ENDIF
    nullify(svr)
END SUBROUTINE
!
SUBROUTINE findidxSPHP(ispf,pltrad_c,idxp,wgtp) ! Fuel Pellet Radius
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ispf
    REAL,INTENT(IN) :: pltrad_c
    INTEGER,INTENT(OUT) :: idxp(2)
    REAL, INTENT(OUT) :: wgtp(2)
    TYPE(SPHvar_type),POINTER :: svr
    INTEGER :: np,ip
    svr=>SPHvar(ispf)
    np=svr%nPLTRAD
    DO ip=1,np
        IF (pltrad_c.lt.svr%PLTRAD(ip)) EXIT
    ENDDO
    IF (ip.eq.1) THEN
        idxp=1
        wgtp(1) = 1._8
        wgtp(2) = 0._8
    ELSEIF (ip.gt.np) THEN
        idxp=np     
        wgtp(1) = 1._8
        wgtp(2) = 0._8       
    ELSE
        idxp(1)=ip-1
        idxp(2)=ip
        wgtp(1)=(svr%PLTRAD(ip)-pltrad_c)/(svr%PLTRAD(ip)-svr%PLTRAD(ip-1))
        wgtp(2)=(pltrad_c-svr%PLTRAD(ip-1))/(svr%PLTRAD(ip)-svr%PLTRAD(ip-1))
    ENDIF
    nullify(svr)
END SUBROUTINE
!
SUBROUTINE findidxSPHT(ispf,temp_c,idxt,wgtt) ! Temperature interpolation...
    IMPLICIT NONE 
    INTEGER,INTENT(IN) :: ispf
    REAL,INTENT(IN) :: temp_c
    INTEGER,INTENT(OUT) :: idxt(2)
    REAL, INTENT(OUT) :: wgtt(2)
    TYPE(SPHvar_type),POINTER :: svr
    INTEGER :: nt,it
    REAL :: t1, t2, tin
    svr=>SPHvar(ispf)
    nt=svr%nTEMP
    DO it=1,nt
        IF (temp_c.lt.svr%TEMP(it)) EXIT
    ENDDO
    IF (it.eq.1) THEN
        idxt=1
        wgtt(1) = 1._8
        wgtt(2) = 0._8
    ELSEIF (it.gt.nt) THEN
        idxt=nt
        wgtt(1) = 1._8
        wgtt(2) = 0._8   
    ELSE
        tin = SQRT(temp_c)
        t1 = SQRT(svr%TEMP(it-1))
        t2 = SQRT(svr%TEMP(it))
        idxt(1)=it-1
        idxt(2)=it
        wgtt(1)=(t2-tin)/(t2-t1)
        wgtt(2)=(tin-t1)/(t2-t1)
    ENDIF
    nullify(svr)
END SUBROUTINE 
!
SUBROUTINE findSPHInterfIdx(ispf,ndratiso_c,iiso,idxi,wgti)
    IMPLICIT NONE
    REAL(4),INTENT(IN) :: ndratiso_c
    INTEGER,INTENT(IN) :: ispf,iiso
    INTEGER,INTENT(OUT) :: idxi(2)
    REAL,INTENT(OUT) :: wgti(2)
    TYPE(SPHvar_type),POINTER :: svr
    TYPE(SPHINTER),POINTER :: sphisointer
    INTEGER :: nndrat,ii
    svr=>SPHvar(ispf)
    sphisointer=>svr%ISOINTER(iiso)
    nndrat=sphisointer%nndrat2u238
    DO ii=1,nndrat
        IF (ndratiso_c.lt.sphisointer%ndrat2u238(ii)) EXIT
    ENDDO
    IF (ii.eq.1) THEN
        idxi=sphisointer%idx
        wgti(1) = 1._8
        wgti(2) = 0._8
    ELSE IF (ii.gt.nndrat) THEN
        idxi=sphisointer%idx+nndrat-1
        wgti(1) = 1._8
        wgti(2) = 0._8
    ELSE
        idxi(1)=sphisointer%idx+ii-2
        idxi(2)=sphisointer%idx+ii-1
        wgti(1)=(sphisointer%ndrat2u238(ii)-ndratiso_c)/(sphisointer%ndrat2u238(ii)-sphisointer%ndrat2u238(ii-1))
        wgti(2)=(ndratiso_c-sphisointer%ndrat2u238(ii-1))/(sphisointer%ndrat2u238(ii)-sphisointer%ndrat2u238(ii-1))
    ENDIF
    nullify(svr)
END SUBROUTINE
! Calculating spectral SPH factor based on the geometrical table...
SUBROUTINE calcSPH_AllG(ispf,isphset,modxsv,u238nd,temp,idiso,pnum,niso,nFxr,nSPHreg_interp,srdidx,igresb,igrese, &
  sphf_localfxr)
    USE TYPEDEF,  ONLY : Cell_Type
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: niso,nFxr,srdidx,nSPHreg_interp,ispf, isphset, igresb, igrese
    INTEGER,INTENT(IN) :: idiso(niso)
    REAL,INTENT(IN) :: modxsv,u238nd,temp,pnum(niso)
    REAL(4),INTENT(OUT) :: sphf_localfxr(igresb:igrese,nFxr)
    TYPE(SPHvar_type),POINTER :: svr
    TYPE(Cell_Type),POINTER :: CellInfo(:)
    
    INTEGER :: idxm(2),idxu(2),idxt(2),idxi(2)
    REAL    :: wgtm(2),wgtu(2),wgtt(2),wgti(2)
    REAL(4) :: sphfm_interp(300,4),sphfp_interp(300,2)
    REAL(4) :: sphcorfm_interp(300,8),sphcorfu_interp(300,4)
    REAL(4) :: sphcorft_interp(300,2),sphcorfi_interp(300)
    REAL(4) :: basesphf(300),corsum(300)
    REAL(4), POINTER :: sphf(:)
    INTEGER :: k,m,it,ip,iu,im,ii,iiso,J,srtidx, ig, ng,srd_start
    INTEGER :: jdiso,iin,nSPHreg_in, nSPHFactor
    REAL(4) :: ndratiso
    REAL    :: Tbeg_ssph, Tend_ssph  ! time measurement for Spectral SPH
    Tbeg_SSPH = nTracer_dclock(FALSE, FALSE)
    svr=>SPHvar(ispf)
    IF (ispf .EQ. 0) THEN
      sphf=>FUEL_SPH(isphset)%sphf_NM
    ELSE
      sphf=>Special_SPH(isphset)%sphf_NM
    END IF
    basesphf=1._4
    CALL findidxSPHM(ispf,modxsv,idxm,wgtm)
    srd_start=sum(svr%SPHreg_srd(0:srdidx-1))
    !
    nSPHreg_in=svr%SPHreg_srd(srdidx)
    ng =igrese - igresb + 1
    nSPHFactor = nSPHreg_interp * ng
    basesphf(1:nSPHFactor) = 0._4
    if (ispf.eq.0) then
      if (ispf.eq.0) CALL findidxSPHU(ispf,u238nd,idxu,wgtu)
      CALL findidxSPHT(ispf,temp,idxt,wgtt)
      DO it=1,2; DO iu =1,2; DO im = 1,2
        srtidx = FUEL_SPH(isphset)%idx_NM(idxm(im),idxu(iu),idxt(it),0)
        basesphf(1:nSPHFactor) = basesphf(1:nSPHFactor)+sphf(srtidx:srtidx+nSPHFactor-1)*wgtm(im)*wgtu(iu)*wgtt(it)
      END DO; END DO; END DO;
    ELSE
      DO im = 1,2
        srtidx = Special_SPH(isphset)%idx_NM(idxm(im))
        basesphf(1:nSPHFactor) = basesphf(1:nSPHFactor)+sphf(srtidx:srtidx+nSPHFactor-1)*wgtm(im)
      END DO
    endif
    ! Calculating correction factor and final spectral SPH factor
    if (ispf.eq.0) then
        corsum=0._4
        DO iin=1,niso
            jdiso=idiso(iin)
            if (jdiso.eq.92238) cycle
            if (jdiso.eq.0) cycle
            ndratiso=pnum(iin)/U238ND
            iiso=svr%MAPISOINTER(jdiso)
            IF (iiso.eq.0) CYCLE
            CALL findSPHInterfIdx(ispf,ndratiso,iiso,idxi,wgti)
            sphcorfm_interp=1._4; sphcorfu_interp=1._4; sphcorft_interp=1._4; sphcorfi_interp=1._4
            sphcorfi_interp(1:nSPHFactor) = 0._4
            DO ii=1,2; DO it=1,2; DO iu =1,2; DO im = 1,2
              srtidx = FUEL_SPH(isphset)%idx_NM(idxm(im),idxu(iu),idxt(it),idxi(ii))
              sphcorfi_interp(1:nSPHFactor) =sphcorfi_interp(1:nSPHFactor) + &
                                        sphf(srtidx:srtidx+nSPHFactor-1)*wgtm(im)*wgtu(iu)*wgtt(it)*wgti(ii)
            END DO; END DO; END DO; END DO;
            corsum=corsum+(sphcorfi_interp-1._4)
        ENDDO ! iin=1,niso
        corsum=corsum+1._4
        sphf_localfxr=1._4
        DO j=1,nSPHreg_interp
          DO ig=igresb,igrese
            sphf_localfxr(ig,nFxr-j+1)=corsum(ig-igresb+1+(j-1)*ng)*basesphf(ig-igresb+1+(j-1)*ng)
          END DO
        ENDDO
    ELSE
        sphf_localfxr=1._4
        DO j=1,nSPHreg_interp
          DO ig=igresb,igrese
            sphf_localfxr(ig,nFxr-j+1)=basesphf(ig-igresb+1+(j-1)*ng)
          END DO
        ENDDO
    endif
    Tend_SSPH = nTracer_dclock(FALSE, FALSE)
    TimeChk%SSPHTime(4) = TimeChk%SSPHTime(4) + Tend_SSPH - Tbeg_SSPH !        calcSPH =
    TimeChk%SSPHCounter(4) = TimeChk%SSPHCounter(4) + 1 ! CalcSPH
END SUBROUTINE
! sph calculation w.r.t. the XS change in each FXR (temperature, burnup)
SUBROUTINE calcPinSSPH(Core, Fxr, PE)
  USE param ,ONLY : sigpH, sigpO, sigpB10, sigpB11
  USE TYPEDEF ,ONLY : CoreInfo_Type, FxrInfo_Type, PE_Type, Cell_Type, Pin_Type
  USE XSLIB_MOD ,ONLY : igresb,igrese,nofghel, nelthel
  USE Timer,       ONLY : nTracer_dclock, TimeChk
  !USE OMP_LIB
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
  INTEGER :: ipin ,icel, iz, ifsr, ifxr, ig, niso, icel0, niso_PIN, isphcell
  INTEGER :: j, k, l, idiso, nSPHreg_interp, nFsrinfxr
  REAL :: MODXSV,TEMP_avg,U238ND_avg,volsum,U238ND,vol
  REAL(4),ALLOCATABLE :: sphf_localfxr(:),sphf_allgfxr(:,:)
  REAL, ALLOCATABLE :: U238ND_PIN(:),pnum_c_avg(:)
  REAL ::  Tbeg2, Tend2
  INTEGER, ALLOCATABLE :: idiso_c_avg(:)
  TYPE(SPHvar_type),POINTER :: svr
  Tbeg2 = nTracer_dclock(FALSE, FALSE)
  ALLOCATE(idiso_c_avg(nelthel))
  ALLOCATE(pnum_c_avg(nelthel))
  myzb = PE%myzb; myze = PE%myze
  Pin => Core%Pin
  CellInfo => Core%CellInfo
  nxy = Core%nxy
  !!$  call omp_set_dynamic(.FALSE.)
  !!$  call omp_set_num_threads(PE%nThread) 
  DO iz=myzb, myze
      IF(.NOT. (Core%lFuelPlane(iz) .OR. Core%lAICPlane(iz))) CYCLE
      ispf=0; svr=>SPHvar(ispf)
      DO ipin = 1, nxy
          !time_beg = nTracer_dclock(FALSE, FALSE)
          FxrIdxSt = Pin(ipin)%FxrIdxSt
          icel = Pin(ipin)%Cell(iz)
          icel0 = CellInfo(icel)%icel0
          isphcell = CellInfo(icel0)%icelSSPH
          IF (.not.CellInfo(icel)%lfuel) CYCLE
          nlocalFxr = CellInfo(icel)%nFxr 
          ALLOCATE(sphf_localfxr(nlocalFxr))
          ALLOCATE(sphf_allgfxr(igresb:igrese,nlocalFXR))
          ALLOCATE(U238ND_PIN(nlocalFXR))
          nSPHreg_interp=CellInfo(icel)%nfueldiv+svr%nrestfxr
          ! calculate U238 average number density in a fuel pellet
          ! collect all isotopes in a fuel pellet
          U238ND_PIN = 0._8
          idiso_c_avg=0; niso_PIN=1
          idiso_c_avg(1)=Fxr(FxrIdxSt + nLocalFxr - 1, iz)%idiso(1)
          DO j = 1, nLocalFXR
            ifxr = FXRIdxSt + j - 1; myFXR => FXR(ifxr, iz)
            DO k = 1, myFXR%niso
              IF (myFxr%idiso(k).eq.92238) U238ND_PIN(j)=myFxr%pnum(k) 
              DO l=1,niso_PIN
                  if (idiso_c_avg(l).eq.myFxr%idiso(k)) EXIT
              END DO
              IF (l.gt.niso_PIN) THEN
                  niso_PIN=niso_PIN+1
                  idiso_c_avg(niso_PIN)=myFxr%idiso(k)
              END IF
            END DO
          END DO
          ! calculate isotope-wise average number density in a fuel pellet
          pnum_c_avg=0._8
          do l=1,niso_PIN
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
          enddo ! l = 1, niso_PIN
          ! calculate MODXSV
          ! calculate average temperature and U238 ND in a fuel pellet
          MODXSV=0._8
          volsum=0._8
          U238ND_avg=0._8
          TEMP_avg=0._8
          !!$OMP PARALLEL DEFAULT(SHARED) &
          !!$OMP PRIVATE(j, ifxr, myFXR, ifsr, vol, nfsrinfxr, k, idiso) REDUCTION(+: MODXSV, volsum, U238ND_avg)
          !!$OMP DO SCHEDULE(GUIDED)
          DO j= 1, nLocalFxr
              ifxr = FxrIdxSt + j -1;  myFxr => Fxr(ifxr, iz)
              ifsr=CellInfo(icel0)%MapFxr2FsrIdx(1,j)
              vol=CellInfo(icel0)%vol(ifsr)
              nfsrinfxr=CellInfo(icel0)%nFsrInFxr(j)
              vol=vol*nfsrinfxr
              IF (.not.myFxr%lSSPHalloc) THEN
                  ALLOCATE(myFxr%SPHfactor(igresb:igrese))
                  myFxr%lSSPHalloc=.true.
              ENDIF
              myFxr%SPHfactor = 1._4
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
            IF (.not.myFxr%lfuel) CYCLE
            IF (U238ND_PIN(j).eq.0) CYCLE
            volsum=volsum+CellInfo(icel)%vol(ifsr)
            U238ND_avg=U238ND_avg+U238ND_PIN(j)*CellInfo(icel)%vol(ifsr)
            TEMP_avg=TEMP_avg+myFxr%Temp*CellInfo(icel)%vol(ifsr) 
          ENDDO !DO j= 1, nLocalFxr
          !!$OMP END DO
          !!$OMP END PARALLEL
          U238ND_avg=U238ND_avg/volsum
          TEMP_avg=TEMP_avg/volsum
          
          DO j = 1, nLocalFxr
            ifxr = FxrIdxSt + j - 1;  myFxr => Fxr(ifxr, iz)
            IF (.not.myFxr%lfuel) CYCLE
            IF (U238ND_PIN(j).eq.0) CYCLE  
            CALL calcSPH_AllG(ispf,isphcell,MODXSV,U238ND_PIN(j),myFxr%Temp,myFxr%idiso(1:myFxr%niso),&
            myFxr%pnum(1:myFxr%niso),myFxr%niso,nlocalFxr,nSPHreg_interp,CellInfo(icel)%srdidx,igresb,igrese,sphf_allgfxr)
            myFxr%SPHfactor(igresb:igrese)=sphf_allgfxr(igresb:igrese,j)
          END DO
          IF (svr%nrestfxr.GE.1) THEN
            CALL calcSPH_AllG(ispf,isphcell,MODXSV,U238ND_avg,TEMP_avg,idiso_c_avg(1:niso_PIN),pnum_c_avg(1:niso_PIN), &
            niso_PIN,nlocalFxr,nSPHreg_interp,CellInfo(icel)%srdidx,igresb,igrese,sphf_allgfxr)
            DO k=1,CellInfo(icel)%ngapdiv
              j=CellInfo(icel)%nfueldiv+k
              ifxr = FxrIdxSt + nLocalFxr - j;  myFxr => Fxr(ifxr, iz)
              myFxr%SPHfactor(igresb:igrese)=sphf_allgfxr(igresb:igrese,nLocalFxr - CellInfo(icel)%nfueldiv)
            ENDDO
            DO k=1,CellInfo(icel)%ncladdiv
              j=CellInfo(icel)%nfueldiv+CellInfo(icel)%ngapdiv+k
              ifxr = FxrIdxSt + nLocalFxr - j;  myFxr => Fxr(ifxr, iz)
              myFxr%SPHfactor(igresb:igrese)=sphf_allgfxr(igresb:igrese,nLocalFxr - CellInfo(icel)%nfueldiv - 1)
            ENDDO
          END IF
          DEALLOCATE(sphf_localfxr, sphf_allgfxr)
          DEALLOCATE(U238ND_PIN)
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
          CALL calcSPH_AllG(ispf,isphcell,MODXSV,U238ND_avg,TEMP_avg,idiso_c_avg(1:niso_PIN),&
          pnum_c_avg(1:niso_PIN),niso_PIN,nlocalFxr,nSPHreg_interp,CellInfo(icel)%srdidx,igresb,igrese,sphf_allgfxr)
          DO j = 1, nLocalFxr
              ifxr = FxrIdxSt + j - 1;  myFxr => Fxr(ifxr, iz)
              if (myFxr%lh2o) cycle
              myFxr%SPHfactor(igresb:igrese)=sphf_allgfxr(igresb:igrese,j)
          ENDDO
          DEALLOCATE(sphf_localfxr)
      ENDDO ! DO ipin = 1, nxy
  ENDDO ! DO iz=myzb, myze
  NULLIFY(Pin,CellInfo,myFxr,svr)
  DEALLOCATE(idiso_c_avg,pnum_c_avg)
  TimeChk%SSPHCounter(3) = TimeChk%SSPHCounter(3) + 1 ! CalcPinSSPH
  Tend2 = nTracer_dclock(FALSE, FALSE)
  TimeChk%SSPHTime(1) = TimeChk%SSPHTime(1) + Tend2 - Tbeg2  ! calcPinSSPH
END SUBROUTINE
!
SUBROUTINE CalcCellSSPH_GEOM(CellInfo)
    !(CellInfo,icel,ndata,iReg,rr,ndiv,ibfuel,iefuel,nfueldiv,lhole) 
    ! EDIT JSU EDIT 20210420 for SSPH optimization
    ! CAUTION-------------------------------------------------------------------------------------------------------*
    !   With this routine, the spectral SPH method can be applied to air-hole containing problems, however,         *
    !   the ring-tickness of each subregions for SPH library are generated with the assumption that the subrings    *
    !   are the results of making equi-volumic rings of one fuel region. Therefore, if one try to use more than one *
    !   air hole regions in input's cell information or use more than one fuel regions, it would deteriorate the    *
    !   result due to unaccurate SPH factor.                                                                        *
    !---------------------------------------------------------------------------------------------------------------*
      USE param ,ONLY : sigpH, sigpO, sigpB10, sigpB11
      USE TYPEDEF ,ONLY : Cell_Type
      USE Material_Mod, ONLY : Mixture
      USE XSLIB_MOD ,ONLY : igresb,igrese,nofghel
      !USE timer,  only : TimeChk
      IMPLICIT NONE
      TYPE(Cell_Type)  ::  CellInfo
      INTEGER,PARAMETER :: ispf=0
      INTEGER :: i,j,k, niso, l, idiso, imat, icel, nfxr, nmodfxr, idiso_c_avg(300),ibfuelFXR,iefuelFXR
      INTEGER :: nSPHreg_interp,nSPHreg0,ig,nfueldiv,ngapdiv,ncladdiv,srdidx,srdmapidx,nfxrdiff,ifxr,idxm(2)
      REAL :: vol,temptemp,Temp,sphpnum_avg,modxsv_c,U238ND,u238nd_c,temp_c, pnum_c_avg(300),wgtm(2)
      REAL(4) :: sphf_interp(40)
      REAL(4), POINTER :: sphf_allgfxr(:,:)
      TYPE(SPHvar_type),POINTER :: svr
      INTEGER :: ndata,ibfuel,iefuel
      INTEGER, POINTER :: matidx(:)
      IF (.not.CellInfo%lfuel) return ! Only Fuel Pin (Do not consider surrounding effect)
      svr=>SPHvar(ispf)
      
      matidx => CellInfo%matfxridx ! from our to inside
      nfxr   = CellInfo%nFXR;        nmodfxr= CellInfo%nmodfxr;
      nfueldiv  = CellInfo%nfueldiv;  ngapdiv   = CellInfo%ngapdiv;
      ncladdiv  = CellInfo%ncladdiv;  srdidx    = CellInfo%srdidx;
      nSPHreg_interp = nfueldiv + svr%nrestfxr
      !
      ALLOCATE(CellInfo%SPHfactor(1:CellInfo%nFXR,igresb:igrese))
      !  1. nuclides(iso) IDs and average number density in fuel region...
      ! finding beginning and ending fuel material...
      DO ibfuel = 1, CellInfo%nmat
        IF (Mixture(CellInfo%matidx(ibfuel))%lfuel) EXIT
      END DO
      !
      DO iefuel = ibfuel, CellInfo%nmat-1
        IF (.NOT.Mixture(CellInfo%matidx(iefuel+1))%lfuel) EXIT
      END DO
      !
      DO ibfuelFXR = 1, CellInfo%nFXR
        IF (Mixture(CellInfo%matfxridx(ibfuelFXR))%lfuel) EXIT
      END DO
      DO iefuelFXR = ibfuelFXR, CellInfo%nFXR-1
        IF (.NOT.Mixture(CellInfo%matfxridx(iefuelFXR+1))%lfuel) EXIT
      END DO
      idiso_c_avg=0; niso=1
      idiso_c_avg(1)=Mixture(CellInfo%matidx(ibfuel))%idiso(1)
      pnum_c_avg=0._8
      DO imat=ibfuel,iefuel
          do l=1,Mixture(CellInfo%matidx(imat))%niso
              do k=1,niso
                  if (idiso_c_avg(k).eq.Mixture(CellInfo%matidx(imat))%idiso(l)) exit
              enddo
              if (k.gt.niso) then
                  niso=niso+1
                  idiso_c_avg(niso)=Mixture(CellInfo%matidx(imat))%idiso(l)
              endif
          enddo
      ENDDO
      do k=1,niso
          vol=0._8
          DO imat=ibfuel,iefuel
              vol=vol+CellInfo%matvol(imat)
              DO l=1,Mixture(CellInfo%matidx(imat))%niso
                  IF (idiso_c_avg(k).eq.Mixture(CellInfo%matidx(imat))%idiso(l)) exit
              ENDDO ! l (iso in Mixture)
              if (l.gt.Mixture(CellInfo%matidx(imat))%niso) cycle
              pnum_c_avg(k)=pnum_c_avg(k)+Mixture(CellInfo%matidx(imat))%pnum(l)*CellInfo%matvol(imat)
          END DO !vol
          pnum_c_avg(k)=pnum_c_avg(k)/vol
      enddo ! k (iso in fuel cell)
      !  2. average number density of spectralSPH factor nuclides(U238 ND and nuclides indicated with#)
      sphpnum_avg=0._8; temptemp=0._8; vol=0._8
      DO imat=ibfuel,iefuel
          vol=vol+CellInfo%matvol(imat)
          DO l=1,Mixture(CellInfo%matidx(imat))%niso
              IF (Mixture(CellInfo%matidx(imat))%lSSPH(l)) exit ! Material contains U238 or isotope starts with#
          ENDDO
          if (l.gt.Mixture(CellInfo%matidx(imat))%niso) then
              write(*,'(a10,i4,a21)') 'Fuel Cell ',icel,' does not have U238.'
              cycle
          endif
          sphpnum_avg=sphpnum_avg+Mixture(CellInfo%matidx(imat))%pnum(l)*CellInfo%matvol(imat)
          temptemp=temptemp+Mixture(CellInfo%matidx(imat))%temp*CellInfo%matvol(imat)
      ENDDO ! fuel ring
      sphpnum_avg=sphpnum_avg/vol
      temptemp=temptemp/vol      ! average temperature [K]
      CellInfo%U238ND0 = sphpnum_avg
      CellInfo%FuelRefTemp0 = temptemp
      !  3. MODXSV (Macroscopic XS of potential XS and volume
      CellInfo%MODXSV0 = 0._8
      DO imat=iefuel+1,CellInfo%nmat
          IF (.not.Mixture(CellInfo%matidx(imat))%lh2o) CYCLE
          DO k=1,Mixture(CellInfo%matidx(imat))%niso
              idiso=Mixture(CellInfo%matidx(imat))%idiso(K)
              IF (idiso.eq.1001) THEN 
                  CellInfo%MODXSV0=CellInfo%MODXSV0+Mixture(CellInfo%matidx(imat))%pnum(k)*sigpH*CellInfo%matvol(imat)
              ELSEIF (idiso.eq.8016) THEN
                  CellInfo%MODXSV0=CellInfo%MODXSV0+Mixture(CellInfo%matidx(imat))%pnum(k)*sigpO*CellInfo%matvol(imat)
              ELSEIF (idiso.eq.5010) THEN
                  CellInfo%MODXSV0=CellInfo%MODXSV0+Mixture(CellInfo%matidx(imat))%pnum(k)*sigpB10*CellInfo%matvol(imat)
              ELSEIF (idiso.eq.5011) THEN
                  CellInfo%MODXSV0=CellInfo%MODXSV0+Mixture(CellInfo%matidx(imat))%pnum(k)*sigpB11*CellInfo%matvol(imat)
              ENDIF
          ENDDO
      ENDDO 
      modxsv_c = CellInfo%MODXSV0;
      u238nd_c = CellInfo%U238ND0;  temp_c   = CellInfo%FuelRefTemp0;
      CellInfo%SPHfactor=1._4
      ALLOCATE(sphf_allgfxr(igresb:igrese,CellInfo%nFXR))
      DO ifxr = ibfuelfxr, iefuelfxr
          imat = CellInfo%matfxridx(ifxr)
          DO l=1,Mixture(imat)%niso
              IF (Mixture(imat)%lSSPH(l)) exit ! finding nuclide using SPH factor (U238 for fuel & nuclides tagged with#)
          ENDDO
          U238ND = Mixture(imat)%pnum(l)
          Temp = Mixture(imat)%temp
          CALL calcSPH_AllG(ispf,CellInfo%icelSSPH,modxsv_c,U238ND,Temp,Mixture(imat)%idiso(1:niso),Mixture(imat)%pnum(1:niso),&
                      Mixture(imat)%niso,CellInfo%nFXR,nSPHreg_interp,srdidx,igresb,igrese,sphf_allgfxr)
          DO ig = igresb, igrese
            CellInfo%SPHfactor(ifxr,ig)=sphf_allgfxr(ig, ifxr)
          END DO
      ENDDO ! j
      CALL calcSPH_AllG(ispf,CellInfo%icelSSPH,modxsv_c,u238nd_c,temp_c,idiso_c_avg,pnum_c_avg,niso,CellInfo%nFXR,&
                        nSPHreg_interp,srdidx,igresb,igrese,sphf_allgfxr)
      IF (svr%nrestfxr.ge.1) THEN
          ! GAP
          DO ig = igresb, igrese
            DO ifxr=nfueldiv+1,nfueldiv+ngapdiv
                CellInfo%SPHfactor(CellInfo%nFXR-ifxr+1,ig)=sphf_allgfxr(ig,CellInfo%nFXR-nfueldiv)
            ENDDO
            ! CLADDING
            IF (svr%nrestfxr.eq.2) THEN
                DO ifxr=nfueldiv+ngapdiv+1,nfueldiv+ngapdiv+ncladdiv
                    CellInfo%SPHfactor(CellInfo%nFXR-ifxr+1,ig)=sphf_allgfxr(ig,CellInfo%nFXR-nfueldiv-1)
                ENDDO
            ENDIF
          END DO
      ENDIF
      DEALLOCATE(sphf_allgfxr)
END SUBROUTINE
!
SUBROUTINE CalcAICCellSSPH_GEOM(CellInfo)
  USE param ,       ONLY : sigpH, sigpO, sigpB10, sigpB11
  USE TYPEDEF ,     ONLY : Cell_Type
  USE Material_Mod, ONLY : Mixture
  USE XSLIB_MOD ,   ONLY : igresb,igrese,nofghel
  IMPLICIT NONE
  TYPE(Cell_Type) :: CellInfo
  INTEGER,PARAMETER :: ispf=1
  INTEGER :: i,j,k, niso, l, idiso, imat, idiso_c_avg(300), nfxr, fuelregidx, ibfuelFXR, iefuelFXR
  INTEGER :: nSPHreg_interp,ig,nfueldiv,ngapdiv,ncladdiv,srdidx,ifxr
  REAL :: vol,Temp,modxsv_c,temp_c,pnum_c_avg(300)
  REAL(4), POINTER :: sphf_allgfxr(:,:)
  TYPE(SPHvar_type),POINTER :: svr
  svr=>SPHvar(ispf)
  !
  DO fuelregidx = 1, CellInfo%nmat
    IF (Mixture(CellInfo%matidx(fuelregidx))%lAIC) EXIT
  END DO
  !
  DO ibfuelFXR = 1, CellInfo%nFXR
    IF (Mixture(CellInfo%matfxridx(ibfuelFXR))%lAIC) EXIT
  END DO
  DO iefuelFXR = ibfuelFXR, CellInfo%nFXR-1
    IF (.NOT.Mixture(CellInfo%matfxridx(iefuelFXR+1))%lAIC) EXIT
  END DO
  !
  niso = Mixture(CellInfo%matidx(fuelregidx))%niso
  CellInfo%MODXSV0 = 0._8
  DO imat=fuelregidx+1,CellInfo%nmat
      IF (.not.Mixture(CellInfo%matidx(imat))%lh2o) CYCLE
      DO k=1,Mixture(CellInfo%matidx(imat))%niso
          idiso=Mixture(CellInfo%matidx(imat))%idiso(K)
          IF (idiso.eq.1001) THEN 
              CellInfo%MODXSV0=CellInfo%MODXSV0+Mixture(CellInfo%matidx(imat))%pnum(k)*sigpH*CellInfo%matvol(imat)
          ELSEIF (idiso.eq.8016) THEN
              CellInfo%MODXSV0=CellInfo%MODXSV0+Mixture(CellInfo%matidx(imat))%pnum(k)*sigpO*CellInfo%matvol(imat)
          ELSEIF (idiso.eq.5010) THEN
              CellInfo%MODXSV0=CellInfo%MODXSV0+Mixture(CellInfo%matidx(imat))%pnum(k)*sigpB10*CellInfo%matvol(imat)
          ELSEIF (idiso.eq.5011) THEN
              CellInfo%MODXSV0=CellInfo%MODXSV0+Mixture(CellInfo%matidx(imat))%pnum(k)*sigpB11*CellInfo%matvol(imat)
          ENDIF
      ENDDO
  ENDDO 
  nfueldiv  =  CellInfo%nfueldiv
  ngapdiv   =  CellInfo%ngapdiv 
  ncladdiv  =  CellInfo%ncladdiv
  srdidx    =  CellInfo%srdidx
  nfxr      =  CellInfo%nFXR  
  nSPHreg_interp = nfueldiv + svr%nrestfxr
  ! calculate default Cell SPH factor
  ALLOCATE(CellInfo%SPHfactor(1:nfxr,igresb:igrese))
  ALLOCATE(sphf_allgfxr(igresb:igrese,CellInfo%nFXR))
  modxsv_c=CellInfo%MODXSV0
  temp_c=CellInfo%FuelRefTemp0
  CellInfo%SPHfactor=1._4
  CALL calcSPH_AllG(ispf,CellInfo%icelSSPH,modxsv_c,0.,temp_c,idiso_c_avg,pnum_c_avg,niso,nfxr,nSPHreg_interp, &
                    srdidx,igresb,igrese,sphf_allgfxr)
  DO ifxr = ibfuelfxr, iefuelfxr
    DO ig = igresb, igrese
      CellInfo%SPHfactor(ifxr,ig)=sphf_allgfxr(ig, ifxr)
    END DO
  END DO
  IF (svr%nrestfxr.ge.1) THEN
    DO ig = igresb, igrese
      DO ifxr=nfueldiv+1,nfueldiv+ngapdiv
          CellInfo%SPHfactor(CellInfo%nFXR-ifxr+1,ig)=sphf_allgfxr(ig,CellInfo%nFXR-nfueldiv)
      ENDDO
      ! CLADDING
      IF (svr%nrestfxr.eq.2) THEN
          DO ifxr=nfueldiv+ngapdiv+1,nfueldiv+ngapdiv+ncladdiv
              CellInfo%SPHfactor(CellInfo%nFXR-ifxr+1,ig)=sphf_allgfxr(ig,CellInfo%nFXR-nfueldiv-1)
          ENDDO
      ENDIF
    END DO
  ENDIF
  DEALLOCATE(sphf_allgfxr)
  nullify(svr)
END SUBROUTINE
!
SUBROUTINE ChkSameFuelStr(iCellinfo, jCellinfo, lsame)
  ! SUBROUTINE TO VERIGY ONLY IF BOTH THE CELL HAVE THE SAME FUEL RADIUS WITH FXR...
  USE PARAM,   ONLY : epsm10
  USE TYPEDEF, ONLY : Cell_TYPE, basicgeom
  IMPLICIT NONE
  TYPE(Cell_Type) :: iCellInfo, jCellInfo
  LOGICAL :: lsame
  TYPE(basicgeom) :: iGeom, jGeom
  INTEGER :: i, j
  iGeom=iCellInfo%geom
  jGeom=jCellInfo%geom
  lsame=.TRUE.
  IF( igeom%ncircle .NE. jgeom%ncircle )THEN
      lsame=.FALSE.
      RETURN
  ENDIF
  IF ( iCellInfo%ibfuel .NE. jCellInfo%ibfuel) THEN
  lsame=.FALSE.
  RETURN
  END IF
  IF ( iCellInfo%iefuel .NE. jCellInfo%iefuel) THEN
  lsame=.FALSE.
  RETURN
  END IF
  IF ( iCellInfo%FuelRad0 .NE. jCellInfo%FuelRad0) THEN
  lsame=.FALSE.
  RETURN
  END IF
  IF ( iCellInfo%nfueldiv .NE. jCellInfo%nfueldiv) THEN
  lsame=.FALSE.
  RETURN
  END IF
  DO i=1, iCellInfo%nfueldiv
  IF( abs(iCellInfo%Rad_CP(i) - jCellInfo%Rad_CP(i)) .GT. epsm10 )THEN
      lsame=.FALSE.
      RETURN
  ENDIF
  ENDDO
END SUBROUTINE
!
SUBROUTINE PreTabGeomSPH_FUEL(Core, nCellTYPE0, PE)
  !USE GEOM,     ONLY : core,       nCellTYPE0
  USE TYPEDEF,      ONLY : cell_type, coreinfo_type, PE_Type
  USE XSLIB_MOD,    ONLY : igresb,igrese
  USE Material_Mod, ONLY : Mixture
  IMPLICIT NONE
  !
  TYPE(coreinfo_type) :: Core
  TYPE(PE_Type) :: PE
  INTEGER :: nCellType0
  !
  TYPE(cell_type), POINTER :: CellInfo(:)
  TYPE(sphvar_type), pointer :: svr
  TYPE(RC_type),POINTER :: DIpos(:)
  TYPE(FUEL_SPH_type), POINTER :: stbl
  INTEGER, POINTER :: idxlib(:,:,:,:,:,:)
  REAL(4),POINTER :: sphf(:), sphf_temp(:,:)
  REAL(4) :: sphin(30),sphout(30),basesphin(30),basesphout(30)
  REAL    :: x1, y1, x2, y2
  REAL    :: wgtp(2), wgtt(2)
  INTEGER :: idxp(2), idxt(2)
  INTEGER :: icel, celstr, isph, isphl
  INTEGER :: ngapdiv, ncladdiv, i, j, nfueldiv, srdidx
  INTEGER :: ntemp, nUnd, nmod, nFXR_cell, nFXR_lib, ng, ninter, nrestfxr
  INTEGER :: itemp, iUND, iMOD, iso, irat, idx_start, iinter, nf, irad, ig, jcel, idx_start2, nf2
  INTEGER :: srd_start, ifxr
  LOGICAL :: lsame
  REAL    :: Tbeg_ssph, Tend_ssph ! time measurement for Spectral SPH
  Tbeg_SSPH = nTracer_dclock(FALSE, FALSE)
  !
  CellInfo => Core%CellInfo
  ALLOCATE(idxSPH2BaseCell(nCellTYPE0))
  ! Fuel ***
  isphl = 0
  svr => sphvar(isphl)
  idxlib=> idxmap(isphl)%idx
  sphf=>SPHdata(isphl)%sphf
  nrestfxr = svr%nrestfxr
  ntemp = svr%nTEMP
  ! nrad  = svr%nPLTRAD
  nund  = svr%nU238ND
  nmod  = svr%nMODXSV
  ninter = svr%ninter
  ng    = igrese - igresb + 1
  DO icel = 1, nCellType0
      IF (CellInfo(icel)%lfuel) EXIT
  END DO
  IF (icel .GT. nCellTYPE0) RETURN ! for the case there is no fuel in all the cell 
  celstr = icel
  nfuelsph = 1
  CellInfo(icel)%icelSSPH = 1
  idxSPH2BaseCell(nfuelsph) = icel
  !
  DO icel = celstr+1, nCellType0
      IF (.NOT.CellInfo(icel)%lfuel) CYCLE
      !--start searching
      DO jcel = 1, nfuelsph
          CALL ChkSameFuelStr(CellInfo(icel), CellInfo(idxSPH2BaseCell(jcel)), lsame)
          IF (lSame) EXIT
      ENDDO
      IF (lSame) THEN
          CellInfo(iCel)%icelSSPH = jcel
      ELSE
          nfuelsph = nfuelsph + 1
          CellInfo(iCel)%icelSSPH = nfuelsph
          idxSPH2BaseCell(nfuelsph) = icel
      ENDIF
  ENDDO
  ! generate geometrical SPH table for Fuel cells
  ALLOCATE(FUEL_SPH(nfuelsph))
  DO isph = 1, nfuelsph
      icel = idxSPH2BaseCell(isph)
      nfueldiv = CellInfo(icel)%nfueldiv;        ngapdiv  =  CellInfo(icel)%ngapdiv 
      ncladdiv = CellInfo(icel)%ncladdiv;        srdidx   = CellInfo(icel)%srdidx   
      nFXR_lib = svr%SPHreg_srd(srdidx)
      nFXR_cell = nfueldiv + ngapdiv + ncladdiv
      srd_start = SUM(svr%SPHreg_srd(0:srdidx-1))
      ALLOCATE(sphf_temp(nFXR_cell,2))
      ALLOCATE(FUEL_SPH(isph)%idx(igresb:igrese,nmod,nUnd,ntemp,0:ninter))
      ALLOCATE(FUEL_SPH(isph)%sphf(ng*nmod*ntemp*nund*(ninter+1)*nFXR_cell))
      ALLOCATE(FUEL_SPH(isph)%idx_NM(nmod,nUnd,ntemp,0:ninter))
      ALLOCATE(FUEL_SPH(isph)%sphf_NM(ng*nmod*ntemp*nund*(ninter+1)*nFXR_cell))
      ! find Pellet radius index
      CALL findidxSPHP(isphl,CellInfo(icel)%FuelRad0,idxp,wgtp)
      nf = 0; nf2 = 0
      DO itemp = 1, ntemp
        DO iUND = 1, nUnd
          DO iMOD = 1, nmod
            FUEL_SPH(isph)%idx_NM(imod,iund,itemp,0)=nf2 + 1
            DO ig=igresb,igrese
                FUEL_SPH(isph)%idx(ig,iMOD,iUND,itemp,0) = nf + 1
                DO irad = 1, 2
                    DIpos => DI(isphl)%spDI(idxp(irad))%DIpos
                    idx_start = idxlib(ig-igresb+1,imod,iund,idxp(irad),itemp,0)
                    sphin(1:nFXR_lib)=sphf(idx_start+srd_start:idx_start+srd_start+nFXR_lib-1)
                    IF (nfxr_cell .EQ. nFXR_lib) THEN
                        sphf_temp(1:nFXR_cell,irad)= sphin(1:nFXR_lib)
                        CYCLE
                    END IF
                    sphout(nFXR_cell-nrestfxr+1:nFXR_cell) = sphin(nFXR_lib-nrestfxr+1:nFXR_lib)
                    DO i=1,nFXR_cell-nrestfxr
                        DO j=1,nFXR_lib-nrestfxr
                            IF (DIpos(nFXR_lib-nrestfxr)%rc(j).gt.DIpos(nFXR_cell-nrestfxr)%rc(i)) EXIT
                        ENDDO
                        IF (j.eq.nFXR_lib-nrestfxr+1) THEN
                            sphout(i)=sphin(nFXR_lib-nrestfxr)
                        ELSE
                            x1=DIpos(nFXR_lib-nrestfxr)%rc(j-1); x2=DIpos(nFXR_lib-nrestfxr)%rc(j);
                            y1=sphin(j-1); y2=sphin(j);
                            sphout(i)=(y2-y1)/(x2-x1)*(DIpos(nFXR_cell-nrestfxr)%rc(i)-x1)+y1          
                        ENDIF ! j
                    ENDDO ! i
                    sphf_temp(1:nFXR_cell,irad)= sphout(1:nFXR_cell)
                END DO ! irad
                FUEL_SPH(isph)%sphf(nf+1:nf+nFXR_cell) = sphf_temp(1:nFXR_cell,1)*wgtp(1)+sphf_temp(1:nFXR_cell,2)*wgtp(2)
                nf = nf + nFXR_Cell
            END DO ! ig
            ! FXR major table...
            DO ifxr = 1, nFXR_Cell
              DO ig = igresb, igrese
                idx_start = FUEL_SPH(isph)%idx(ig,iMOD,iUND,itemp,0)
                nf2 = nf2 + 1
                FUEL_SPH(isph)%sphf_NM(nf2) = FUEL_SPH(isph)%sphf(idx_start+ifxr-1)
              END DO ! ig
            END DO ! ifxr 
          END DO ! imod
        END DO ! iund
      END DO ! itemp
      DO iinter = 1, ninter
        DO itemp = 1, ntemp
          DO iUND = 1, nUnd
            DO iMOD = 1, nmod
              FUEL_SPH(isph)%idx_NM(imod,iund,itemp,iinter)=nf2 + 1
              DO ig=igresb,igrese
                  FUEL_SPH(isph)%idx(ig,iMOD,iUND,itemp,iinter) = nf + 1
                  DO irad = 1, 2
                      DIpos => DI(isphl)%spDI(idxp(irad))%DIpos
                      idx_start = idxlib(ig-igresb+1,imod,iund,idxp(irad),itemp,iinter)
                      idx_start2 = idxlib(ig-igresb+1,imod,iund,idxp(irad),itemp,0)
                      basesphin(1:nFXR_lib)=sphf(idx_start2+srd_start:idx_start2+srd_start+nFXR_lib-1)
                      sphin(1:nFXR_lib)=sphf(idx_start+srd_start:idx_start+srd_start+nFXR_lib-1)
                      IF (nfxr_cell .EQ. nFXR_lib) THEN
                          sphf_temp(1:nFXR_cell,irad)= sphin(1:nFXR_lib)
                          CYCLE
                      END IF
                      sphin(1:nFXR_lib)=sphin(1:nFXR_lib)*basesphin(1:nFXR_lib)
                      sphout(nFXR_cell-nrestfxr+1:nFXR_cell) = sphin(nFXR_lib-nrestfxr+1:nFXR_lib)
                      basesphout(nFXR_cell-nrestfxr+1:nFXR_cell) = basesphin(nFXR_lib-nrestfxr+1:nFXR_lib)
                      DO i=1,nFXR_cell-nrestfxr
                          DO j=1,nFXR_lib-nrestfxr
                              IF (DIpos(nFXR_lib-nrestfxr)%rc(j).gt.DIpos(nFXR_cell-nrestfxr)%rc(i)) EXIT
                          ENDDO
                          IF (j.eq.nFXR_lib-nrestfxr+1) THEN
                              sphout(i)=sphin(nFXR_lib-nrestfxr)
                              basesphout(i)=basesphin(nFXR_lib-nrestfxr)
                          ELSE
                              x1=DIpos(nFXR_lib-nrestfxr)%rc(j-1); x2=DIpos(nFXR_lib-nrestfxr)%rc(j);
                              y1=sphin(j-1); y2=sphin(j);
                              sphout(i)=(y2-y1)/(x2-x1)*(DIpos(nFXR_cell-nrestfxr)%rc(i)-x1)+y1   
                              y1=basesphin(j-1); y2=basesphin(j);       
                              basesphout(i)=(y2-y1)/(x2-x1)*(DIpos(nFXR_cell-nrestfxr)%rc(i)-x1)+y1          
                          ENDIF ! j
                      ENDDO ! i
                      sphf_temp(1:nFXR_cell,irad)= sphout(1:nFXR_cell)/basesphout(1:nFXR_cell)
                  END DO ! irad
                  FUEL_SPH(isph)%sphf(nf+1:nf+nFXR_cell) = sphf_temp(1:nFXR_cell,1)*wgtp(1)+sphf_temp(1:nFXR_cell,2)*wgtp(2)
                  nf = nf + nFXR_Cell
              END DO ! ig
              ! FXR major table...
              DO ifxr = 1, nFXR_Cell
                DO ig = igresb, igrese
                  idx_start = FUEL_SPH(isph)%idx(ig,iMOD,iUND,itemp,iinter)
                  nf2 = nf2 + 1
                  FUEL_SPH(isph)%sphf_NM(nf2) = FUEL_SPH(isph)%sphf(idx_start+ifxr-1)
                END DO ! ig
              END DO ! ifxr 
            END DO ! imod
          END DO ! iund
        END DO ! itemp
      END DO ! iinter
      DEALLOCATE(sphf_temp)
  END DO ! isph
  ! Calculate SPH factor for fuel cell
  DO icel = 1, nCellType0
    IF (.NOT.CellInfo(icel)%lfuel) CYCLE
    CALL CalcCellSSPH_GEOM(Cellinfo(icel))
  END DO ! icel
  Tend_SSPH = nTracer_dclock(FALSE, FALSE)
  TimeChk%SSPHCounter(1) = TimeChk%SSPHCounter(1) + 1
  TimeChk%SSPHTime(2) = TimeChk%SSPHTime(2) + Tend_SSPH - Tbeg_SSPH  ! PreTabFuelSSPH (moderator xsv)
END SUBROUTINE
!
SUBROUTINE PreTabGeomSPH_AIC(Core, nCellTYPE0, PE)
  USE TYPEDEF,      ONLY : cell_type, coreinfo_type, PE_Type
  USE XSLIB_MOD,    ONLY : igresb,igrese
  USE Material_Mod, ONLY : Mixture
  IMPLICIT NONE
  !
  TYPE(coreinfo_type) :: Core
  TYPE(PE_Type) :: PE
  INTEGER :: nCellType0
  !
  TYPE(cell_type), POINTER :: CellInfo(:)
  TYPE(sphvar_type), pointer :: svr
  TYPE(RC_type),POINTER :: DIpos(:)
  INTEGER,POINTER :: idxlib_sp(:,:,:,:)
  REAL(4),POINTER :: sphf(:), sphf_temp(:,:)
  REAL(4) :: sphin(30),sphout(30)
  REAL    :: x1, y1, x2, y2
  REAL    :: wgtp(2), wgtt(2)
  INTEGER :: idxp(2), idxt(2)
  INTEGER :: icel, celstr, isph, isphl
  INTEGER :: ngapdiv, ncladdiv, i, j, k, nfueldiv, srdidx
  INTEGER :: ntemp, nmod, nFXR_cell, nFXR_lib, ng, nrestfxr
  INTEGER :: itemp, iMOD, idx_start, nf, irad, ig, jcel
  INTEGER :: srd_start, ifxr, nf2
  LOGICAL :: lsame
  REAL    :: Tbeg_ssph, Tend_ssph ! time measurement for Spectral SPH
  Tbeg_SSPH = nTracer_dclock(FALSE, FALSE)
  !
  CellInfo => Core%CellInfo
  ALLOCATE(idxAICSPH2BaseCell(nCellTYPE0))
  ! AIC ***
  isphl = 1
  svr => sphvar(isphl)
  idxlib_sp => idxmap(isphl)%spidx
  sphf=>SPHdata(isphl)%sphf
  nrestfxr = svr%nrestfxr
  ntemp = svr%nTEMP
  nmod  = svr%nMODXSV
  ng    = igrese - igresb + 1
  DO icel = 1, nCellType0
      IF (CellInfo(icel)%lAIC) EXIT
  END DO
  IF (icel .GT. nCellTYPE0) RETURN ! for the case there is no AIC pin in all the cell 
  celstr = icel
  nAICSPH = 1
  CellInfo(icel)%icelSSPH = 1
  idxAICSPH2BaseCell(nAICSPH) = icel
  !
  DO icel = celstr+1, nCellType0
      IF (.NOT.CellInfo(icel)%lAIC) CYCLE
      !--start searching
      DO jcel = 1, nAICSPH
          CALL ChkSameFuelStr(CellInfo(icel), CellInfo(idxAICSPH2BaseCell(jcel)), lsame)
          EXIT
      ENDDO
      IF (lSame) THEN
          CellInfo(iCel)%icelSSPH = jcel
      ELSE
          nAICSPH = nAICSPH + 1
          CellInfo(iCel)%icelSSPH = nAICSPH
          idxAICSPH2BaseCell(nAICSPH) = icel
      ENDIF
  ENDDO
  ALLOCATE(SPECIAL_SPH(nAICSPH))
  DO isph = 1, nAICSPH
      icel = idxAICSPH2BaseCell(isph)
      nfueldiv = CellInfo(icel)%nfueldiv;        ngapdiv  =  CellInfo(icel)%ngapdiv 
      ncladdiv = CellInfo(icel)%ncladdiv;        srdidx   = CellInfo(icel)%srdidx   
      nFXR_lib = svr%SPHreg_srd(srdidx)
      nFXR_cell = nfueldiv + ngapdiv + ncladdiv
      srd_start = SUM(svr%SPHreg_srd(0:srdidx-1))
      ALLOCATE(sphf_temp(nFXR_cell,4))
      ALLOCATE(SPECIAL_SPH(isph)%idx(igresb:igrese,nmod))
      ALLOCATE(SPECIAL_SPH(isph)%sphf(ng*nmod*nFXR_cell))
      ALLOCATE(SPECIAL_SPH(isph)%idx_nm(nmod))
      ALLOCATE(SPECIAL_SPH(isph)%sphf_NM(ng*nmod*nFXR_cell))
      ! find Pellet radius index
      CALL findidxSPHP(isphl,CellInfo(icel)%FuelRad0,idxp,wgtp)
      CALL findidxSPHT(isphl,CellInfo(icel)%FuelRefTemp0,idxT,wgtT)
      nf = 0; nf2 = 0;
      DO iMOD = 1, nmod
          SPECIAL_SPH(isph)%idx_NM(imod)=nf2 + 1
          DO ig=igresb,igrese
              SPECIAL_SPH(isph)%idx(ig,iMOD) = nf + 1
              SPECIAL_SPH(isph)%sphf(nf+1:nf+nFXR_cell) = 0._4
              DO irad = 1, 2
                  DO itemp = 1, 2
                      DIpos => DI(isphl)%spDI(idxp(irad))%DIpos
                      idx_start = idxlib_sp(ig-igresb+1,imod,idxp(irad),idxt(itemp))
                      sphin(1:nFXR_lib)=sphf(idx_start+srd_start:idx_start+srd_start+nFXR_lib-1)
                      IF (nfxr_cell .EQ. nFXR_lib) THEN
                          sphf_temp(1:nFXR_cell,irad+(itemp-1)*2) = sphin(1:nFXR_cell)
                      ELSE
                        sphout(nFXR_cell-nrestfxr+1:nFXR_cell) = sphin(nFXR_lib-nrestfxr+1:nFXR_lib)
                        DO i=1,nFXR_cell-nrestfxr
                            DO j=1,nFXR_lib-nrestfxr
                                IF (DIpos(nFXR_lib-nrestfxr)%rc(j).gt.DIpos(nFXR_cell-nrestfxr)%rc(i)) EXIT
                            ENDDO
                            IF (j.eq.nFXR_lib-nrestfxr+1) THEN
                                sphout(i)=sphin(nFXR_lib-nrestfxr)
                            ELSE
                                x1=DIpos(nFXR_lib-nrestfxr)%rc(j-1); x2=DIpos(nFXR_lib-nrestfxr)%rc(j);
                                y1=sphin(j-1); y2=sphin(j);
                                sphout(i)=(y2-y1)/(x2-x1)*(DIpos(nFXR_cell-nrestfxr)%rc(i)-x1)+y1          
                            ENDIF ! j
                        ENDDO ! i
                        sphf_temp(1:nFXR_cell,irad+(itemp-1)*2)= sphout(1:nFXR_cell)
                        SPECIAL_SPH(isph)%sphf(nf+1:nf+nFXR_cell) = SPECIAL_SPH(isph)%sphf(nf+1:nf+nFXR_cell)+sphout(1:nFXR_cell)*wgtp(irad)*wgtt(itemp)!+sphf_temp(1:nFXR_cell,2)*wgtp(2)
                      END IF
                  END DO ! itemp
              END DO ! irad
              SPECIAL_SPH(isph)%sphf(nf+1:nf+nFXR_cell) = sphf_temp(1:nFXR_cell,1)*wgtp(1)*wgtt(1) + sphf_temp(1:nFXR_cell,2)*wgtp(2)*wgtt(1) + &
                                                          sphf_temp(1:nFXR_cell,3)*wgtp(1)*wgtt(2) + sphf_temp(1:nFXR_cell,4)*wgtp(2)*wgtt(2)
              nf = nf + nFXR_Cell
          END DO ! ig
          ! FXR major table...
          DO ifxr = 1, nFXR_Cell
              DO ig = igresb, igrese
                idx_start = SPECIAL_SPH(isph)%idx(ig,iMOD)
                nf2 = nf2 + 1
                SPECIAL_SPH(isph)%sphf_NM(nf2) = SPECIAL_SPH(isph)%sphf(idx_start+ifxr-1)
              END DO ! ig
          END DO ! ifxr 
      END DO ! imod
      DEALLOCATE(sphf_temp)
  END DO ! isph
  ! Calculate SPH factor for AIC cell
  DO icel = 1, nCellType0
    IF (.NOT.CellInfo(icel)%lAIC) CYCLE
    CALL CalcAICCellSSPH_GEOM(Cellinfo(icel))
  END DO ! icel
  Tend_SSPH = nTracer_dclock(FALSE, FALSE)
  TimeChk%SSPHCounter(2) = TimeChk%SSPHCounter(2) + 1
  TimeChk%SSPHTime(3) = TimeChk%SSPHTime(3) + Tend_SSPH - Tbeg_SSPH  ! PreTabFuelSSPH (moderator xsv)
END SUBROUTINE
END Module
