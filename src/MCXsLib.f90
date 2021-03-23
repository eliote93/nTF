MODULE MCXsLib
    USE MCDefine
    USE TYPEDEF,    ONLY : FXRInfo_Type
    IMPLICIT NONE
    INTEGER :: nfxr, nz, nchi
    TYPE(XsecSet), POINTER :: mcxs(:, :)
    TYPE(XsecSet), POINTER :: micxs(:)
    TYPE(FXRInfo_Type), POINTER :: CoreFXR(:,:)
    REAL(8), POINTER, Dimension(:,:,:) :: Fxrphi        !spectrum for chi generation (ifxr,iz,ig)
    REAL(8), POINTER, Dimension(:,:,:,:) :: Fxrphiomp   ! (ifxr, iz, ig, tid)
    INTEGER, POINTER :: fxr2xsmap(:,:), xs2fxrmap(:,:)    !
    REAL(8), POINTER :: MicChiCDF(:,:) ! iso, g
    !ALLOCATE(Fxr2XsMap(nCoreFxr, nz)) ! FXR idx -> xs idx
    !ALLOCATE(Xs2FxrMap(nCoreFxr*nz,2))  ! xs idx -> FXR idx

CONTAINS
SUBROUTINE SetMCXsLib(Core, Fxr, PE, CtrlMC)
    USE TYPEDEF,    ONLY : CoreInfo_Type, FXRInfo_Type, Pin_Type, Cell_Type, XsMac_Type, PE_TYPE
    USE Xsec4MC, ONLY : g2gc
    USE MacXsLib_Mod,  ONLY : MacXsBase, MacXsScatMatrix, MacP1XsScatMatrix, MacP2XsScatMatrix, MacP3XsScatMatrix
    USE Core_mod,   ONLY: GroupInfo
    USE BasicOperation, ONLY : MULTI_VA
    USE CNTL, ONLY : nTracerCntl
    USE HighOrderSC, ONLY : SetAnIsoGaussian
    USE XsUtil_mod, ONLY : FreeXsMac, FreeXsIsoMac
    IMPLICIT NONE
    type(coreinfo_type), intent(in) :: core
    TYPE(FXRInfo_Type), POINTER :: FXR(:,:)
    TYPE(Pin_Type), POINTER :: Pin(:)
    TYPE(Cell_Type), POINTER :: CellInfo(:)
    TYPE(PE_Type) :: PE
    TYPE(FXRInfo_Type), POINTER :: myfxr, tfxr
    TYPE(XsMac_Type) :: XsMac
    TYPE(XsecSet), POINTER :: xs
    TYPE(MCCTRL) :: ctrlMC
    real(8), pointer, dimension(:,:) :: SctCDF1, xssm2g1, FisCdf
    real(8), pointer, dimension(:) :: xsaa1
    REAL(8) :: xsct, xsnf, sphfac
    INTEGER :: g, gp, ng, niso, ipin, icel, FxrIdxSt, j
    INTEGER :: nCoreFXR, scatod
    INTEGER :: ifxr, iz, iResoGrp1, iResoGrp2, ixs, jfxr, jz, iso
    LOGICAL :: lres, lfuel, lfuelcel
    INTEGER :: nxstype
    LOGICAL :: lexist
    CoreFxr=>Fxr
    Pin => Core%Pin; CellInfo => Core%CellInfo
    scatod=nTracerCntl%scatod
    ng = GroupInfo%ng
    iresoGrp1 = GroupInfo%nofg + 1; iresoGrp2 = GroupInfo%nofg + GroupInfo%norg
    allocate(g2gc(ctrlMC%ngMC))
    g2gc(1:ctrlMC%gThrm)=1
    g2gc(ctrlMC%gThrm+1:ctrlMC%ngMC)=2
    nCoreFXR=Core%nCoreFXR
    nfxr=nCoreFxr
    nz=Core%nz
    ALLOCATE(Fxr2XsMap(nCoreFxr, nz)) ! FXR idx -> xs idx
    ALLOCATE(Xs2FxrMap(nCoreFxr*nz,2))  ! xs idx -> FXR idx
    fxr2xsmap=0; xs2fxrmap=0;
    nxstype=0;
    DO iz = 1, nz
        DO ifxr = 1, nCoreFxr
            myFXR=> FXR(ifxr, iz)
            IF( myFXR%niso .EQ. 0 ) CYCLE !luse=.FALSE.
            lexist=.FALSE.  !exist in the list
            DO ixs = 1, nxstype
                jfxr=Xs2Fxrmap(ixs,1)
                jz=Xs2Fxrmap(ixs,2)
                tfxr=>FXR(jfxr,jz)
                IF( CheckFXR(myFxr, tfxr) )THEN
                    lexist=.TRUE.
                    fxr2xsmap(ifxr,iz)=ixs
                    EXIT
                ENDIF
            ENDDO
            IF( .NOT. lexist )THEN
                nxstype=nxstype+1
                fxr2xsmap(ifxr,iz)=nxstype
                xs2fxrmap(ixs,1)=ifxr
                xs2fxrmap(ixs,2)=iz
            ENDIF
        ENDDO
    ENDDO
    ALLOCATE(micxs(nxstype))
    WRITE(*,'(a11, i3)') ' nXStype = ', nxstype
    DO ixs = 1, nxstype
        ifxr=xs2fxrmap(ixs,1)
        iz=xs2fxrmap(ixs,2)
        myfxr=>FXR(ifxr,iz)
        ipin=myFxr%ipin; icel = Pin(ipin)%Cell(iz); FxrIdxSt = Pin(ipin)%FxrIdxSt; j = ifxr - FxrIdxSt + 1
        lres=myFxr%lres
        lfuel=myFxr%lfuel; lfuelcel=CellInfo(icel)%lfuel
        xs=>micxs(ixs)
        xs%lfuel=lfuel
        xs%niso=myFXR%niso
        xs%idiso=>myFXR%idiso
        xs%pnum=>myFXR%pnum
        niso=xs%niso
        CALL ALLOCMCXS(xs, ng)
        xs%ifxr=ifxr; xs%iz=iz
        xs%ixs=ixs
        CALL MacXsBase(XSMac, myFxr, 1, ng, ng, 1._8, FALSE, TRUE)
        IF( nTracerCntl%lSSPH .and. lfuelcel) THEN
            DO g=iresoGrp1,iresoGrp2
                sphfac = CellInfo(icel)%SPHfactor(j,g)
                XsMac%xsmact(g)  =XsMac%xsmact(g)  *sphfac
                XsMac%xsmaca(g)  =XsMac%xsmaca(g)  *sphfac
                XsMac%xsmacs(g)  =XsMac%xsmacs(g)  *sphfac
                XsMac%xsmactr(g) =XsMac%xsmactr(g) *sphfac
                XsMac%xsmacstr(g)=XsMac%xsmacstr(g)*sphfac
                XsMac%xsmacf(g)  =XsMac%xsmacf(g)  *sphfac
                XsMac%xsmacnf(g) =XsMac%xsmacnf(g) *sphfac
                XsMac%xsmackf(g) =XsMac%xsmackf(g) *sphfac
                DO iso = 1, niso
                    XsMac%IsoXsMacT(iso, g) =XsMac%IsoXsMacT(iso, g) *sphfac
                    XsMac%IsoXsMacA(iso, g) =XsMac%IsoXsMacA(iso, g) *sphfac
                    XsMac%IsoXsMacS0(iso, g)=XsMac%IsoXsMacS0(iso, g)*sphfac
                    XsMac%IsoXsMacS1(iso, g)=XsMac%IsoXsMacS1(iso, g)*sphfac
                    XsMac%IsoXsMacSS(iso, g)=XsMac%IsoXsMacSS(iso, g)*sphfac
                    XsMac%IsoXsMactr(iso, g)=XsMac%IsoXsMactr(iso, g)*sphfac
                    XsMac%IsoXsMacf(iso, g) =XsMac%IsoXsMacf(iso, g) *sphfac
                    XsMac%IsoXsMacNf(iso, g)=XsMac%IsoXsMacNf(iso, g)*sphfac
                    XsMac%IsoXsMacKf(iso, g)=XsMac%IsoXsMacKf(iso, g)*sphfac
                ENDDO
            ENDDO
        ENDIF
        xs%xsa=XsMac%XsMaca
        xs%xss=XsMac%XsMacs
        IF( lres .AND. CtrlMc%lsgfsp )THEN
            do g = iresoGrp1,iresoGrp2
              xs%xsa(g) = xs%xsa(g) * myFxr%fresoa(g)
            enddo
            IF (nTracerCntl%lRST) THEN
              do g = iresoGrp1,iresoGrp2
                xs%xss(g) = xs%xss(g) * myFxr%fresos(g)
              enddo
            ENDIF
            DO g= 1, ng
                xs%xst(g)=xs%xsa(g)+xs%xss(g)
            ENDDO
        ELSE
            xs%xst=XsMac%XsMact
        ENDIF
        xs%xstr=XsMac%XsMactr

        xs%xsnf=XsMac%XsMacnf
        xs%xskf=XsMac%XsMackf

        IF( lfuel )THEN
            IF( lres .AND. CtrlMc%lsgfsp )THEN
              do g = iresoGrp1,iresoGrp2
                xs%xsnf(g) = xs%xsnf(g) * myFxr%fresoNF(g)  
                xs%xskf(g) = xs%xskf(g) * myFxr%fresokF(g)  
              enddo
            ENDIF
            !call AdrFisCDF(xs, txs)-------------------------------------
            ALLOCATE(FisCDF(ng,niso))
            FisCDF=0.
            IF( lres .AND. CtrlMc%lsgfsp )THEN
                DO g=1, ng
                    iso=1
                    xsnf=XsMac%IsoXsMacnf(iso, g)
                    IF( g .GE. iresoGrp1 .AND. g .LE. iresoGrp2 )THEN
                        xsnf = xsnf * myFxr%fresoFiso(iso, g)
                    ENDIF
                    FisCDF(g, iso)=xsnf
                    DO iso = 2, niso
                        xsnf=XsMac%IsoXsMacnf(iso, g)
                        IF( g .GE. iresoGrp1 .AND. g .LE. iresoGrp2 )THEN
                            xsnf=xsnf*myFxr%fresoFiso(iso, g)
                        ENDIF
                        FisCDF(g, iso)=FisCDF(g, iso-1)+xsnf
                    ENDDO
                    FisCDF(g,:)=FisCDF(g,:)/FisCDF(g,niso)
                ENDDO
            ELSE
                DO g=1, ng
                    iso=1
                    xsnf=XsMac%IsoXsMacnf(iso, g)
                    FisCDF(g, iso)=xsnf
                    DO iso = 2, niso
                        xsnf=XsMac%IsoXsMacnf(iso, g)
                        FisCDF(g, iso)=FisCDF(g, iso-1)+xsnf
                    ENDDO
                    FisCDF(g,:)=FisCDF(g,:)/FisCDF(g,niso)
                ENDDO
            ENDIF
            xs%FisCDF=FisCDF
            CALL InitialEffChi(xs, myFxr) ! for first cycle
            DEALLOCATE(FisCdf)
        ENDIF

        CALL MacXsScatMatrix(XsMac, myFxr, 1, ng, ng, GroupInfo, TRUE)

        xs%xssm=XsMac%XsMacSm
        xs%xss0=XsMac%XsMacSm
        IF(ScatOd .GE. 1) THEN
            CALL MacP1XsScatMatrix(XsMac, myFxr, 1, ng, ng, GroupInfo)
            xs%xss1=XsMac%XsMacP1Sm
        ENDIF
        IF(ScatOd .GE. 2)THEN
            CALL MacP2XsScatMatrix(XsMac, myFxr, 1, ng, ng, GroupInfo)
            xs%xss2=XsMac%XsMacP2Sm
        ENDIF
        IF(ScatOd .EQ. 3)THEN
            CALL MacP3XsScatMatrix(XsMac, myFxr, 1, ng, ng, GroupInfo)
            xs%xss3=XsMac%XsMacP3Sm
        ENDIF

        !call AdrSctCDF(xs, txs)-------------------------------------
        allocate(SctCDF1(ng,ng))
        SctCDF1=0.
        do g=1, ng
            xsct=xs%xss0(g,1)
            if (xsct<0) then
                xsct=0.
            endif

            SctCDF1(1,g)=xsct
            do gp=2, ng
                xsct=xs%xss0(g,gp)
                if (xsct<0) then
                    xsct=0.
                endif

                SctCDF1(gp,g)=SctCDF1(gp-1,g)+xsct
            enddo
            SctCDF1(:,g)=SctCDF1(:,g)/SctCDF1(ng,g)
        enddo
        xs%SctCDF=SctCDF1

        !call AdrSct2G(xs, txs)--------------------------------------
        allocate(xssm2g1(ng,2))
        xssm2g1=0.
        do g=1, ng
            do gp=1, ng
                xssm2g1(g,g2gc(gp))=xssm2g1(g,g2gc(gp))+xs%xss0(g,gp)
            enddo
        enddo
        xs%xssm2g=xssm2g1

        CALL SetAnIsoGaussian(xs, scatod)

        !call AdrAbsrA(xs, txs)--------------------------------------
        allocate(xsaa1(ng))
        do g=1, ng
            do gp=1, ng
                if (xs%xss0(g,gp)<0) then
                    if (g==gp) xsaa1(g)=-xs%xss0(g,gp)
                    xs%xss0(g,gp)=0
                endif
            enddo
        enddo
        xs%xsaa=xsaa1
        DEALLOCATE(xsaa1, SctCDF1, xssm2g1)
        !DEALLOCATE(xs%xss0)
        CALL FreeXsMac(XsMac)
        CALL FreeXsIsoMac(XsMac)
    ENDDO
    CALL SetMicChiCDF()

    !CHI update
    !CALL InitEffChi()
ENDSUBROUTINE
FUNCTION CheckFXR(fxr1, fxr2) ! fxr1 : subject to check its uniqueness, fxr2 : object of comparison
    USE TYPEDEF, ONLY : FXRInfo_Type
    LOGICAL :: CheckFXR
    LOGICAL :: lsg
    TYPE(FxrInfo_TYPE) :: fxr1, fxr2
    INTEGER :: i
    checkfxr=.TRUE. ! fxr1 == fxr2
    IF( fxr1%lres )THEN
        checkfxr=.FALSE.
        RETURN
    ELSE
        IF( fxr1%niso .NE. fxr2%niso )THEN
            checkfxr=.FALSE.
            RETURN
        ELSE
            DO i=1, fxr1%niso
                IF( fxr1%idiso(i) .NE. fxr2%idiso(i) )THEN
                    checkfxr=.FALSE.
                    RETURN
                ELSE
                    IF( fxr1%pnum(i) .NE. fxr2%pnum(i) )THEN
                        checkfxr=.FALSE.
                        RETURN
                    ENDIF
                ENDIF
            ENDDO
        ENDIF
    ENDIF
ENDFUNCTION
SUBROUTINE SetMicChiCDF()
    USE XSLIB_MOD, ONLY : nelthel, ldiso
    USE Core_mod, ONLY : GroupInfo
    IMPLICIT NONE
    INTEGER :: ng, niso
    INTEGER :: g, iso
    REAL :: sum
    niso=nelthel
    ng=GroupInfo%ng
    nchi=GroupInfo%nchi
    ALLOCATE(MicChiCdf(niso, nchi))
    DO iso = 1, niso
        IF( ldiso(iso)%ichi .EQ. 0 ) CYCLE
        g=1
        MicChiCdf(iso,g)=ldiso(iso)%chi(g)
        DO g = 2, nchi
            MicChiCdf(iso,g)=MicChiCdf(iso,g-1)+ldiso(iso)%chi(g)
        ENDDO
        MicChiCDF(iso,:)=MicChiCDF(iso,:)/MicChiCDF(iso,nchi)
    ENDDO
ENDSUBROUTINE

SUBROUTINE SetMCXsLib_old(Core, Fxr, PE, CtrlMC) ! save all macroscopic for all region
    USE TYPEDEF,    ONLY : CoreInfo_Type, FXRInfo_Type, XsMac_Type, PE_TYPE
    USE Xsec4MC, ONLY : g2gc
    USE MacXsLib_Mod,  ONLY : MacXsBase, MacXsScatMatrix, MacP1XsScatMatrix, MacP2XsScatMatrix, MacP3XsScatMatrix
    USE Core_mod,   ONLY: GroupInfo
    USE BasicOperation, ONLY : MULTI_VA
    USE CNTL, ONLY : nTracerCntl
    IMPLICIT NONE
    type(coreinfo_type), intent(in) :: core
    TYPE(FXRInfo_Type), POINTER :: FXR(:,:)
    TYPE(PE_Type) :: PE
    TYPE(FXRInfo_Type), POINTER :: myfxr
    TYPE(XsMac_Type) :: XsMac
    TYPE(XsecSet), POINTER :: xs
    TYPE(MCCTRL) :: ctrlMC
    real(8), pointer, dimension(:,:) :: SctCDF1, xssm2g1
    real(8), pointer, dimension(:) :: xsaa1
    REAL(8) :: xsct
    INTEGER :: g, gp, ng
    INTEGER :: nCoreFXR, scatod
    INTEGER :: ifxr, iz, iResoGrp1, iResoGrp2
    LOGICAL :: lres, lfuel
    CoreFxr=>Fxr
    scatod=nTracerCntl%scatod
    ng = GroupInfo%ng
    iresoGrp1 = GroupInfo%nofg + 1; iresoGrp2 = GroupInfo%nofg + GroupInfo%norg
    allocate(g2gc(ctrlMC%ngMC))
    g2gc(1:ctrlMC%gThrm)=1
    g2gc(ctrlMC%gThrm+1:ctrlMC%ngMC)=2
    nCoreFXR=Core%nCoreFXR
    nfxr=nCoreFxr
    nz=Core%nz
    ALLOCATE(mcxs(nCoreFxr, nz))

    DO iz= 1, nz
    DO ifxr=1, nCoreFXR
        myFxr=>FXR(ifxr,iz)
        xs=>mcxs(ifxr, iz)
        CALL ALLOCMCXS(xs, ng)
        lres=myFxr%lres
        lfuel=myFxr%lfuel
        xs%ifxr=ifxr; xs%iz=iz
        CALL MacXsBase(XSMac, myFxr, 1, ng, ng, 1._8, FALSE, TRUE, FALSE)
        xs%xst=XsMac%XsMact
        xs%xstr=XsMac%XsMactr
        xs%xsnf=XsMac%XsMacnf
        xs%xskf=XsMac%XsMackf
        IF( lfuel )THEN
          do g = iresoGrp1,iresoGrp2  
            xs%xsnf(g) = xs%xsnf(g) * myFxr%fresoNF(g)
            xs%xskf(g) = xs%xskf(g) * myFxr%fresoKF(g)
          enddo
        ENDIF
        !CALL UpdtEffChi()
        xs%xsa=XsMac%XsMaca
        IF( lres )THEN
            do g = iresoGrp1,iresoGrp2
              xs%xsa(g) = xs%xsa(g) * myFxr%fresoa(g)
            enddo
            IF (nTracerCntl%lRST) THEN
              do g = iresoGrp1,iresoGrp2
                xs%xss(g) = xs%xss(g) * myFxr%fresos(g)
              enddo
            ENDIF
        ENDIF
        CALL MacXsScatMatrix(XsMac, myFxr, 1, ng, ng, GroupInfo, TRUE, TRUE)
        xs%xssm=XsMac%XsMacSm
        xs%xss0=XsMac%XsMacSm
        IF(ScatOd .GE. 1) CALL MacP1XsScatMatrix(XsMac, myFxr, 1, ng, ng, GroupInfo)
        IF(ScatOd .GE. 2) CALL MacP2XsScatMatrix(XsMac, myFxr, 1, ng, ng, GroupInfo)
        IF(ScatOd .EQ. 3) CALL MacP3XsScatMatrix(XsMac, myFxr, 1, ng, ng, GroupInfo)
        xs%xss1=XsMac%XsMacP1Sm
        xs%xss2=XsMac%XsMacP2Sm
        xs%xss3=XsMac%XsMacP3Sm

        !call AdrSctCDF(xs, txs)
        allocate(SctCDF1(ng,ng))
        SctCDF1=0.
        do g=1, ng
            xsct=xs%xss0(g,1)
            if (xsct<0) then
                xsct=0.
            endif

            SctCDF1(1,g)=xsct
            do gp=2, ng
                xsct=xs%xss0(g,gp)
                if (xsct<0) then
                    xsct=0.
                endif

                SctCDF1(gp,g)=SctCDF1(gp-1,g)+xsct
            enddo
            SctCDF1(:,g)=SctCDF1(:,g)/SctCDF1(ng,g)
        enddo
        xs%SctCDF=SctCDF1

        !call AdrSct2G(xs, txs)
        allocate(xssm2g1(ng,2))
        xssm2g1=0.
        do g=1, ng
            do gp=1, ng
                xssm2g1(g,g2gc(gp))=xssm2g1(g,g2gc(gp))+xs%xss0(g,gp)
            enddo
        enddo
        xs%xssm2g=xssm2g1

        !call AdrAbsrA(xs, txs)
        allocate(xsaa1(ng))
        do g=1, ng
            do gp=1, ng
                if (xs%xss0(g,gp)<0) then
                    if (g==gp) xsaa1(g)=-xs%xss0(g,gp)
                    xs%xss0(g,gp)=0
                endif
            enddo
        enddo
        xs%xsaa=xsaa1
    ENDDO
    ENDDO
    !CHI update
    CALL InitEffChi()
ENDSUBROUTINE
SUBROUTINE ALLOCMCXS(xs, ng)
    USE MCDefine
    USE Core_mod,   ONLY: GroupInfo
    USE allocs
    IMPLICIT NONE
    TYPE(XsecSet) :: xs
    INTEGER :: ng, nchi, niso
    xs%idx=1
    niso=xs%niso
    CALL Dmalloc0(xs%xst, 1, ng)
    CALL Dmalloc0(xs%xstr, 1, ng)
    CALL Dmalloc0(xs%xsnf, 1, ng)
    CALL Dmalloc0(xs%xskf, 1, ng)
    CALL Dmalloc0(xs%xsa, 1, ng)
    CALL Dmalloc0(xs%xss, 1, ng)
    CALL Dmalloc0(xs%xssm, 1, ng, 1, ng)
    CALL Dmalloc0(xs%xss0, 1, ng, 1, ng)
    CALL Dmalloc0(xs%xss1, 1, ng, 1, ng)
    CALL Dmalloc0(xs%xss2, 1, ng, 1, ng)
    CALL Dmalloc0(xs%xss3, 1, ng, 1, ng)
    CALL Dmalloc0(xs%Sctcdf, 1, ng, 1, ng)
    CALL Dmalloc0(xs%xssm2g, 1, ng, 1, 2)
    CALL Dmalloc0(xs%xsaa, 1, ng)
    IF( xs%lfuel )THEN
        CALL Dmalloc0(xs%fiscdf, 1, ng, 1, niso)
        CALL Dmalloc0(xs%FisCDF, 1, ng, 1, niso)
        nchi=GroupInfo%nchi
        CALL Dmalloc0(xs%chi, 1, nchi)
    ENDIF

ENDSUBROUTINE

SUBROUTINE UpdtEffChi()
    USE MCDefine
    USE Core_mod,   ONLY: GroupInfo
    USE TYPEDEF,    ONLY : FXRInfo_Type
    USE MacXsLib_mod,   ONLY : GetMacChi
    IMPLICIT NONE
    INTEGER :: iz, ifxr, nchi, ig, ng
    REAL(8), POINTER :: Spectrum(:)
    TYPE(XsecSet), Pointer :: xs
    TYPE(FXRInfo_Type), POINTER :: myfxr
    nchi = GroupInfo%nchi
    ng=GroupInfo%ng
    ALLOCATE(spectrum(ng))
    DO iz= 1, nz
        DO ifxr=1, nFXR
            myFxr=>CoreFXR(ifxr,iz)
            IF( myFxr%lfuel )THEN
                xs=>mcxs(ifxr, iz)
                DO ig = 1, ng
                    Spectrum(ig)=Fxrphi(ifxr,iz,ig)
                ENDDO
                CALL GetMacChi(myFXR, Spectrum, 1, nchi, nchi, ng, .FALSE.)
                xs%chi=>myFxr%chi
            ENDIF
        ENDDO
    ENDDO
ENDSUBROUTINE
!
SUBROUTINE InitEffChi()
    USE MCDefine
    USE Core_mod,   ONLY: GroupInfo
    USE TYPEDEF,    ONLY : FXRInfo_Type
    USE MacXsLib_mod,   ONLY : GetMacChi
    IMPLICIT NONE
    INTEGER :: iz, ifxr, nchi, ig, ng
    REAL(8), POINTER :: Spectrum(:)
    TYPE(XsecSet), Pointer :: xs
    TYPE(FXRInfo_Type), POINTER :: myfxr
    nchi = GroupInfo%nchi
    ng=GroupInfo%ng
    ALLOCATE(spectrum(ng))
    DO iz= 1, nz
        DO ifxr=1, nFXR
            myFxr=>CoreFXR(ifxr,iz)
            IF( myFxr%lfuel)THEN
                xs=>mcxs(ifxr, iz)
                DO ig = 1, ng
                    Spectrum(ig)=1._8
                ENDDO
                ALLOCATE(myFxr%chi(nchi))
                CALL GetMacChi(myFXR, Spectrum, 1, nchi, nchi, ng, .FALSE.)
                xs%chi=>myFxr%chi
            ENDIF
        ENDDO
    ENDDO
ENDSUBROUTINE
SUBROUTINE InitialEffChi(xs, myfxr)
    USE MCDefine
    USE Core_mod,   ONLY: GroupInfo
    USE TYPEDEF,    ONLY : FXRInfo_Type
    USE MacXsLib_mod,   ONLY : GetMacChi
    IMPLICIT NONE
    INTEGER :: iz, ifxr, nchi, ig, ng
    REAL(8), POINTER :: Spectrum(:)
    TYPE(XsecSet), Pointer :: xs
    TYPE(FXRInfo_Type), POINTER :: myfxr
    nchi = GroupInfo%nchi
    ng=GroupInfo%ng
    ALLOCATE(Spectrum(ng))
    DO ig = 1, ng
        Spectrum(ig)=1._8
    ENDDO
    ALLOCATE(myFxr%chi(nchi))
    CALL GetMacChi(myFXR, Spectrum, 1, nchi, nchi, ng, .FALSE.)
    xs%chi=myFxr%chi
    DEALLOCATE(spectrum)
ENDSUBROUTINE

ENDMODULE
