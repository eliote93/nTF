SUBROUTINE GCHom_FAST(Core, THInfo, FmInfo, GroupInfo, nTracerCntl, PE, ng)
    USE TYPEDEF,      ONLY : CoreInfo_Type, FmInfo_Type, CMInfo_Type, THInfo_Type, GroupInfo_Type,PE_TYPE, XsMac_Type
    USE CNTL,         ONLY : nTracerCntl_Type
    USE GC_mod
    USE XsUtil_mod, ONLY : FreeXsIsoMac
    USE BenchXs,      ONLY : MacXsBen
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FMInfo_Type) :: FmInfo
    TYPE(THInfo_Type) :: THInfo
    TYPE(GroupInfo_Type) :: GroupInfo
    TYPE(nTracerCntl_Type) :: nTracerCntl
    TYPE(PE_TYPE) :: PE
    TYPE(XsMac_Type) :: XsMac
    INTEGER :: ng
    
    INTEGER :: i, j, k 
    INTEGER :: ig, ig2
    
    INTEGER :: iso, isoidx, id, idiso, nisoinFxr, imix
    INTEGER :: ixy, ixyl
    INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, iCel, iFxr, ifsr, nFsrInFxr
    REAL :: myPhiVol, Vol, myfxrvol(100)
    REAL :: FSV, isoFSV ! Fission Source*Volume
    
INTERFACE
    SUBROUTINE GenFxrIsoMacXS(XsMac, IsoMacXsSm, IsoMacXsSmP1, Fxr, ifxr, Tempref, ig1, ig2, Core, ipin, iz, GroupInfo, nTRACERCntl, PE)
    USE PARAM
    USE TYPEDEF,      ONLY : CoreInfo_Type, XsMac_Type,  FxrInfo_Type,  GroupInfo_Type,  PE_Type
    USE CNTL,          ONLY : nTracerCntl_Type
    USE MacXsLib_mod
    USE BasicOperation, ONLY : CP_VA
    USE GC_mod, ONLY : isosize
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(XsMac_Type) :: XsMac
    TYPE(FxrInfo_Type),pointer :: Fxr(:,:)
    TYPE(GroupInfo_Type) :: GroupInfo
    TYPE(nTracerCntl_Type) :: nTracerCntl
    TYPE(PE_TYPE) :: PE
    INTEGER :: ig1, ig2, ipin, iz, ifxr
    REAL :: TempRef
    REAL :: IsoMacXsSm(isosize,GroupInfo%ng,GroupInfo%ng)
    REAL :: IsoMacXsSmP1(isosize,GroupInfo%ng,GroupInfo%ng)
    END SUBROUTINE
END INTERFACE

    FSV=0; !fission source *Volume
    IF( nTracerCntl%lXsLib )tempref = THInfo%RefFuelTemp(iz) + CKELVIN
    Vol=0
    DO ixyl = 1, nxy 
        ixy = Asy(iasy)%GlobalPinIdx(ixyl) 
        icel = Pin(ixy)%Cell(iz)
        Vol=Vol+Core%PinVol(ixy,iz)
        FxrIdxSt = Pin(ixy)%FxrIdxSt
        FsrIdxSt = Pin(ixy)%FsrIdxSt
        nlocalFxr = Cell(icel)%nFxr
        DO j = 1, nLocalFxr
            ifxr = FxrIdxSt + j -1
            myFxr => FmInfo%FXR(ifxr,iz)
            nFsrInFxr = myFxr%nFsrInFxr
            !FsrIdxSt = myFxr%FsrIdxSt
            nisoInFxr = myFxr%niso
            IF( nTracerCntl%lXsLib ) CALL GenFxrIsoMacXS(XsMac, IsoMacXsSm, IsoMacXsSmP1, FmInfo%FXR, ifxr, Tempref, 1, ng, Core, ixy, iz, GroupInfo, nTRACERCntl, PE)
            ! 13/11/01
            DO k = 1, nFsrInFxr
                !ifsr = FsrIdxSt + k -1
                !--- EDIT 14/12/23 Sq. cell problem
                ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                DO ig = 1, ng
                    PinPhiVol(ig,ixy)=PinPhiVol(ig,ixy)+FsrPhiVol(ifsr,ig)
                ENDDO
            ENDDO        
            ! 13/11/01 end
            myFXRvol=0.0
            DO k = 1, nFsrInFxr
                ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                DO ig = 1, ng
                    U238phi=U238phi+      FsrPhiVol(ifsr,ig)
                    myFXRvol(ig)=myFXRvol(ig)+FsrPhiVol(ifsr,ig)
                ENDDO !---END of G sweep
            ENDDO !---END of Fsr sweep
            IF( nTracerCntl%lXsLib )THEN
                DO iso = 1, nisoInFxr
                    idiso=myFxr%idiso(iso)
                    IF( idiso .EQ. 8001 )THEN
                        idiso=8016
                    ENDIF
                    id = MapNucl(idiso)  !id = Map(92235) = 37
                    !isodata => ldiso(id)
                    isoidx=isoList(id,1)
                    IF( idiso .EQ. 8016 )THEN
                        IF( myFxr%lFuel )THEN
                            isoidx=isoidx+1
                        ENDIF
                    ENDIF
                    isoFSV=0._8
                    DO ig2 = 1, ng
                        isoFSV=isoFSV+XsMac%isoXsMacKF(iso,ig2) *myFXRvol(ig2)
                    ENDDO
                    FSV=FSV+isoFSV
                    DO ig = 1, ng   
                        myphivol=myFXRvol(ig)
                        !--- MacroScopic XS ----------
                        !MAC : (total), D, a, r, nf, kf, kf                        
                        isoMacXs(0,0,ig)=isoMacXs(0,0,ig)+XsMac%isoXsMacT(iso,ig) *myPhiVol   !Total
                        isoMacXs(0,1,ig)=isoMacXs(0,1,ig)+XsMac%isoXsMacTR(iso,ig) *myPhiVol  !D
                        isoMacXs(0,2,ig)=isoMacXs(0,2,ig)+XsMac%isoXsMacA(iso,ig) *myPhiVol   !absorption                        
                        isoMacXs(0,3,ig)=isoMacXs(0,3,ig)+XsMac%isoXsMacA(iso,ig) *myPhiVol   !removal base
                        
                        DO ig2 = 1, ng
                            IF( ig2 .NE. ig )THEN
                                 isoMacXs(0,3,ig)=isoMacXs(0,3,ig)+IsoMacXsSm(iso,ig,ig2) *myPhiVol   !removal
                            ENDIF 
                        ENDDO
                        isoMacXs(0,4,ig)=isoMacXs(0,4,ig)+XsMac%isoXsMacNF(iso,ig) *myPhiVol  !nu-fission
                        isoMacXs(0,5,ig)=isoMacXs(0,5,ig)+XsMac%isoXsMacKF(iso,ig) *myPhiVol  !kappa-fission
                        IF( ig .LE. groupinfo%nchi .AND. isoFSV .NE. 0)THEN
                            isoMacXs(0,6,ig)=isoMacXs(0,6,ig)+myFXR%CHI(ig) *isoFSV           !CHI
                        ENDIF
                        DO ig2 = 1, ng
                            isoMacSm(0,ig,ig2)  =isoMacSm(0,ig,ig2)  +IsoMacXsSm(iso,ig,ig2)  *myPhiVol !MacroScattering
                            isoMacSmP1(0,ig,ig2)=isoMacSmP1(0,ig,ig2)+IsoMacXsSmP1(iso,ig,ig2)*myPhiVol
                        ENDDO
                        
                        !--- MicroScopic XS ----------
                        !MIC : tr, a, r, f, nu, k
                        isoMacXs(isoidx,0,ig)=isoMacXs(isoidx,0,ig)+XsMac%isoXsMacT(iso,ig) *myPhiVol  !Total
                        isoMacXs(isoidx,1,ig)=isoMacXs(isoidx,1,ig)+XsMac%isoXsMacTR(iso,ig) *myPhiVol  !D                        
                        isoMacXs(isoidx,2,ig)=isoMacXs(isoidx,2,ig)+XsMac%isoXsMacA(iso,ig)  *myPhiVol  !absorption                        
                        !--- transport removal
                        !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+(IsoMacXsSm(iso,ig,ig))*myPhiVol  !ss
                        !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+(XsMac%isoXsMacTR(iso,ig)-IsoMacXsSm(iso,ig,ig))*myPhiVol  !removal
                        !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+(XsMac%isoXsMacTR(iso,ig))*myPhiVol  !tr
                        !--- total removal
                        !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+XsMac%isoXsMacA(iso,ig)  *myPhiVol  ! abs
                        DO ig2 = 1, ng
                            IF( ig2 .NE. ig )THEN
                                 !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+IsoMacXsSm(iso,ig,ig2) *myPhiVol  !removal
                            ENDIF 
                        ENDDO
                        isoMacXs(isoidx,4,ig)=isoMacXs(isoidx,4,ig)+XsMac%isoXsMacF(iso,ig)  *myPhiVol  !fission
                        isoMacXs(isoidx,5,ig)=isoMacXs(isoidx,5,ig)+XsMac%isoXsMacNF(iso,ig) *myPhiVol  !nu-fission
                        isoMacXs(isoidx,6,ig)=isoMacXs(isoidx,6,ig)+XsMac%isoXsMacKF(iso,ig) *myPhiVol  !kappa-fission
                        DO ig2 = 1, ng
                            isoMacSm(isoidx,ig,ig2)  =isoMacSm(isoidx,ig,ig2)  +IsoMacXsSm(iso,ig,ig2)  *myPhiVol
                            isoMacSmP1(isoidx,ig,ig2)=isoMacSmP1(isoidx,ig,ig2)+IsoMacXsSmP1(iso,ig,ig2)*myPhiVol
                        ENDDO
                    ENDDO !---END of G sweep
                    IF( idiso.EQ.92238 )THEN !--- U238 N2N
                        DO k = 1, nFsrInFxr
                            ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                            DO ig = 1, ng
                                CALL GenU238N2N(myn2n, TempRef, ig)
                                U238n2n=U238n2n+myn2n*myfxr%pnum(iso)*FsrPhiVol(ifsr,ig)
                            ENDDO !---END of G sweep
                        ENDDO !---END of Fsr sweep
                    ELSE
                    ENDIF ! END OF U238 N2N
                ENDDO !---END of Iso sweep
            ELSE
                iso=1
                imix=myFxr%imix
                myfxrvol=0._8
                DO k = 1, nFsrInFxr
                    !ifsr = FsrIdxSt + k -1
                    !--- EDIT 14/12/23 Sq. cell problem
                    ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                    isoFSV=0._8
                    DO ig2 = 1, ng
                        isoFSV=isoFSV+MacXsBen(imix)%XsKF(ig2) *FsrPhiVol(ifsr,ig2)
                    ENDDO
                    FSV=FSV+isoFSV
                    DO ig = 1, ng
                        myFXRvol(ig)=myFXRvol(ig)+FsrPhiVol(ifsr,ig)
                    ENDDO
                ENDDO   
                DO ig = 1, ng
                    myphivol=myFXRvol(ig)
                    !--- MacroScopic XS ----------
                    !MAC : (total), D, a, r, nf, kf, kf                        
                    isoMacXs(0,0,ig)=isoMacXs(0,0,ig)+MacXsBen(imix)%XsT(ig) *myPhiVol   !Total
                    isoMacXs(0,1,ig)=isoMacXs(0,1,ig)+MacXsBen(imix)%XsTR(ig) *myPhiVol  !D
                    isoMacXs(0,2,ig)=isoMacXs(0,2,ig)+MacXsBen(imix)%Xsa(ig) *myPhiVol   !absorption                        
                    isoMacXs(0,3,ig)=isoMacXs(0,3,ig)+MacXsBen(imix)%Xsa(ig)  *myPhiVol   !removal base
                    
                    DO ig2 = 1, ng
                        IF( ig2 .NE. ig )THEN
                             isoMacXs(0,3,ig)=isoMacXs(0,3,ig)+MacXsBen(imix)%XSS0(ig,ig2)  *myPhiVol   !removal
                        ENDIF 
                    ENDDO
                    isoMacXs(0,4,ig)=isoMacXs(0,4,ig)+MacXsBen(imix)%Xsnf(ig) *myPhiVol  !nu-fission
                    isoMacXs(0,5,ig)=isoMacXs(0,5,ig)+MacXsBen(imix)%Xskf(ig) *myPhiVol  !kappa-fission
                    IF( isoFSV .NE. 0)THEN
                        isoMacXs(0,6,ig)=isoMacXs(0,6,ig)+MacXsBen(imix)%chi(ig) *isoFSV           !CHI
                    ENDIF
                    DO ig2 = 1, ng
                        isoMacSm(0,ig,ig2)  =isoMacSm(0,ig,ig2)  +MacXsBen(imix)%XSS0(ig,ig2)  *myPhiVol !MacroScattering
                        !isoMacSmP1(0,ig,ig2)=isoMacSmP1(0,ig,ig2)+MacXsBen(imix)%XSS1(ig,ig2)  *myPhiVol !MacroScattering
                    ENDDO
                    isoidx=1
                    !--- MicroScopic XS ----------
                    !MIC : tr, a, r, f, nu, k
                    isoMacXs(isoidx,0,ig)=isoMacXs(isoidx,0,ig)+MacXsBen(imix)%XsT(ig) *myPhiVol  !Total
                    isoMacXs(isoidx,1,ig)=isoMacXs(isoidx,1,ig)+MacXsBen(imix)%XsTr(ig) *myPhiVol  !D                        
                    isoMacXs(isoidx,2,ig)=isoMacXs(isoidx,2,ig)+MacXsBen(imix)%Xsa(ig) *myPhiVol  !absorption                        
                    !--- transport removal
                    !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+(IsoMacXsSm(iso,ig,ig))*myPhiVol  !ss
                    !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+(XsMac%isoXsMacTR(iso,ig)-IsoMacXsSm(iso,ig,ig))*myPhiVol  !removal
                    !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+(XsMac%isoXsMacTR(iso,ig))*myPhiVol  !tr
                    !--- total removal
                    !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+XsMac%isoXsMacA(iso,ig)  *myPhiVol  ! abs
                    DO ig2 = 1, ng
                        IF( ig2 .NE. ig )THEN
                             !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+IsoMacXsSm(iso,ig,ig2) *myPhiVol  !removal
                        ENDIF 
                    ENDDO
                    isoMacXs(isoidx,4,ig)=isoMacXs(isoidx,4,ig)+MacXsBen(imix)%Xsnf(ig)  *myPhiVol  !fission
                    isoMacXs(isoidx,5,ig)=isoMacXs(isoidx,5,ig)+MacXsBen(imix)%Xsnf(ig) *myPhiVol  !nu-fission
                    isoMacXs(isoidx,6,ig)=isoMacXs(isoidx,6,ig)+MacXsBen(imix)%Xskf(ig) *myPhiVol  !kappa-fission
                    DO ig2 = 1, ng
                        isoMacSm(isoidx,ig,ig2)  =isoMacSm(isoidx,ig,ig2)  +MacXsBen(imix)%XSS0(ig,ig2)  *myPhiVol
                        !isoMacSmP1(isoidx,ig,ig2)=isoMacSmP1(isoidx,ig,ig2)+MacXsBen(imix)%XSS1(ig,ig2)*myPhiVol
                    ENDDO
                ENDDO !---END of G sweep
            ENDIF            
            !DEALLOCATE(IsoMacXsSm)
            CALL FreeXsIsoMac(XsMac)
        ENDDO !---END of Fxr sweep
    ENDDO !---END of nPin sweep
    DO ig = 1, ng
        DO ixyl = 1, nxy 
            ixy = Asy(iasy)%GlobalPinIdx(ixyl) 
            Phi(ig)       =Phi(ig)+PinPhiVol(ig,ixy)  !PhiVol sum
        ENDDO
        Phi(ig)=Phi(ig)/Vol
    ENDDO
    
    iso=0; !homogenized Macro XS 47G iso=0
    DO ig = 1, ng
        DO i = 0, nxs-1
            isoMacXs(iso,i,ig)=isoMacXs(iso,i,ig)/PhiVol(ig)
        ENDDO
        XSt(ig)=isoMacXs(iso,0,ig)
        IF( FSV .NE. 0 )THEN
            isoMacXs(iso,6,ig)=isoMacXs(iso,6,ig)/FSV !chi
        ELSE
            isoMacXs(iso,6,ig)=0.0 !chi for non-fissile material
        ENDIF            
        DO ig2 = 1, ng
            isoMacSm(iso,ig,ig2)  =isoMacSm(iso,ig,ig2)  /PhiVol(ig)
            isoMacSmP1(iso,ig,ig2)=isoMacSmP1(iso,ig,ig2)/PhiVol(ig)
        ENDDO
    ENDDO
    
    
    !--- MicroScopic XS ----------
    DO iso = 1, nisotot
        DO ig = 1, ng
            DO i = 0, nxs
                isoMacXs(iso,i,ig)=isoMacXs(iso,i,ig)/PhiVol(ig)
            ENDDO
            DO ig2 = 1, ng
                isoMacSm(iso,ig,ig2)  =isoMacSm(iso,ig,ig2)  /PhiVol(ig)
                isoMacSmP1(iso,ig,ig2)=isoMacSmP1(iso,ig,ig2)/PhiVol(ig)
            ENDDO
        ENDDO
    ENDDO
    
ENDSUBROUTINE

SUBROUTINE GCHom_SA(Core, THInfo, FmInfo, GroupInfo, nTracerCntl, PE, ng)
    USE TYPEDEF,      ONLY : CoreInfo_Type, FmInfo_Type, CMInfo_Type, THInfo_Type, GroupInfo_Type,PE_TYPE, XsMac_Type
    USE CNTL,         ONLY : nTracerCntl_Type
    USE GC_mod
    USE XsUtil_mod, ONLY : FreeXsIsoMac
    USE BenchXs,      ONLY : MacXsBen
    USE TH_Mod,           ONLY : GetPinFuelTemp
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FMInfo_Type) :: FmInfo
    TYPE(THInfo_Type) :: THInfo
    TYPE(GroupInfo_Type) :: GroupInfo
    TYPE(nTracerCntl_Type) :: nTracerCntl
    TYPE(PE_TYPE) :: PE
    TYPE(XsMac_Type) :: XsMac
    INTEGER :: ng
    
    INTEGER :: i, j, k 
    INTEGER :: ig, ig2
    
    INTEGER :: iso, isoidx, id, idiso, nisoinFxr, imix
    INTEGER :: ixy, ixyl
    INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, iCel, iFxr, ifsr, nFsrInFxr
    REAL :: myPhiVol, Vol
    REAL :: FSV, isoFSV ! Fission Source*Volume
    
INTERFACE
    SUBROUTINE GenFxrIsoMacXS(XsMac, IsoMacXsSm, IsoMacXsSmP1, Fxr, ifxr, Tempref, ig1, ig2, Core, ipin, iz, GroupInfo, nTRACERCntl, PE)
    USE PARAM
    USE TYPEDEF,      ONLY : CoreInfo_Type, XsMac_Type,  FxrInfo_Type,  GroupInfo_Type,  PE_Type
    USE CNTL,          ONLY : nTracerCntl_Type
    USE MacXsLib_mod
    USE BasicOperation, ONLY : CP_VA
    USE GC_mod, ONLY : isosize
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(XsMac_Type) :: XsMac
    TYPE(FxrInfo_Type),pointer :: Fxr(:,:)
    TYPE(GroupInfo_Type) :: GroupInfo
    TYPE(nTracerCntl_Type) :: nTracerCntl
    TYPE(PE_TYPE) :: PE
    INTEGER :: ig1, ig2, ipin, iz, ifxr
    REAL :: TempRef
    REAL :: IsoMacXsSm(isosize,GroupInfo%ng,GroupInfo%ng)
    REAL :: IsoMacXsSmP1(isosize,GroupInfo%ng,GroupInfo%ng)
    END SUBROUTINE
END INTERFACE
    
    FSV=0; !fission source *Volume
    Vol=0
    DO ixyl = 1, nxy 
        ixy = Asy(iasy)%GlobalPinIdx(ixyl) 
        icel = Pin(ixy)%Cell(iz)
        Vol=Vol+Core%PinVol(ixy,iz)
        FxrIdxSt = Pin(ixy)%FxrIdxSt
        FsrIdxSt = Pin(ixy)%FsrIdxSt
        nlocalFxr = Cell(icel)%nFxr
        IF( nTracerCntl%lXsLib ) Tempref = GetPinFuelTemp(Core, FmInfo%FXR, iz, ixy)
        DO j = 1, nLocalFxr
            ifxr = FxrIdxSt + j -1
            myFxr => FmInfo%FXR(ifxr,iz)
            nFsrInFxr = myFxr%nFsrInFxr
            !FsrIdxSt = myFxr%FsrIdxSt
            nisoInFxr = myFxr%niso
            IF( nTracerCntl%lXsLib ) CALL GenFxrIsoMacXS(XsMac, IsoMacXsSm, IsoMacXsSmP1, FmInfo%FXR, ifxr, Tempref, 1, ng, Core, ixy, iz, GroupInfo, nTRACERCntl, PE)
            ! 13/11/01
            DO k = 1, nFsrInFxr
                !ifsr = FsrIdxSt + k -1
                !--- EDIT 14/12/23 Sq. cell problem
                ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                DO ig = 1, ng
                    PinPhiVol(ig,ixy)=PinPhiVol(ig,ixy)+FsrPhiVol(ifsr,ig)
                ENDDO
            ENDDO        
            ! 13/11/01 end
            
            IF( nTracerCntl%lXsLib )THEN
                DO iso = 1, nisoInFxr
                    idiso=myFxr%idiso(iso)
                    IF( idiso .EQ. 8001 )THEN
                        idiso=8016
                    ENDIF
                    id = MapNucl(idiso)  !id = Map(92235) = 37
                    !isodata => ldiso(id)
                    isoidx=isoList(id,1)
                    IF( idiso .EQ. 8016 )THEN
                        IF( myFxr%lFuel )THEN
                            isoidx=isoidx+1
                        ENDIF
                    ENDIF
                    DO k = 1, nFsrInFxr
                        !ifsr = FsrIdxSt + k -1
                        !--- EDIT 14/12/23 Sq. cell problem
                        ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                        isoFSV=0
                        DO ig2 = 1, ng
                            isoFSV=isoFSV+XsMac%isoXsMacKF(iso,ig2) *FsrPhiVol(ifsr,ig2)
                        ENDDO
                        FSV=FSV+isoFSV
                        DO ig = 1, ng
                            myPhiVol=FsrPhiVol(ifsr,ig)
                            
                            !--- MacroScopic XS ----------
                            !MAC : (total), D, a, r, nf, kf, kf                        
                            isoMacXs(0,0,ig)=isoMacXs(0,0,ig)+XsMac%isoXsMacT(iso,ig) *myPhiVol   !Total
                            isoMacXs(0,1,ig)=isoMacXs(0,1,ig)+XsMac%isoXsMacTR(iso,ig) *myPhiVol  !D
                            isoMacXs(0,2,ig)=isoMacXs(0,2,ig)+XsMac%isoXsMacA(iso,ig) *myPhiVol   !absorption                        
                            isoMacXs(0,3,ig)=isoMacXs(0,3,ig)+XsMac%isoXsMacA(iso,ig) *myPhiVol   !removal base
                            
                            DO ig2 = 1, ng
                                IF( ig2 .NE. ig )THEN
                                     isoMacXs(0,3,ig)=isoMacXs(0,3,ig)+IsoMacXsSm(iso,ig,ig2) *myPhiVol   !removal
                                ENDIF 
                            ENDDO
                            isoMacXs(0,4,ig)=isoMacXs(0,4,ig)+XsMac%isoXsMacNF(iso,ig) *myPhiVol  !nu-fission
                            isoMacXs(0,5,ig)=isoMacXs(0,5,ig)+XsMac%isoXsMacKF(iso,ig) *myPhiVol  !kappa-fission
                            IF( ig .LE. groupinfo%nchi .AND. isoFSV .NE. 0)THEN
                                isoMacXs(0,6,ig)=isoMacXs(0,6,ig)+myFXR%CHI(ig) *isoFSV           !CHI
                            ENDIF
                            DO ig2 = 1, ng
                                isoMacSm(0,ig,ig2)  =isoMacSm(0,ig,ig2)  +IsoMacXsSm(iso,ig,ig2)  *myPhiVol !MacroScattering
                                isoMacSmP1(0,ig,ig2)=isoMacSmP1(0,ig,ig2)+IsoMacXsSmP1(iso,ig,ig2)*myPhiVol
                            ENDDO
                            
                            !--- MicroScopic XS ----------
                            !MIC : tr, a, r, f, nu, k
                            isoMacXs(isoidx,0,ig)=isoMacXs(isoidx,0,ig)+XsMac%isoXsMacT(iso,ig) *myPhiVol  !Total
                            isoMacXs(isoidx,1,ig)=isoMacXs(isoidx,1,ig)+XsMac%isoXsMacTR(iso,ig) *myPhiVol  !D                        
                            isoMacXs(isoidx,2,ig)=isoMacXs(isoidx,2,ig)+XsMac%isoXsMacA(iso,ig)  *myPhiVol  !absorption                        
                            !--- transport removal
                            !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+(IsoMacXsSm(iso,ig,ig))*myPhiVol  !ss
                            !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+(XsMac%isoXsMacTR(iso,ig)-IsoMacXsSm(iso,ig,ig))*myPhiVol  !removal
                            !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+(XsMac%isoXsMacTR(iso,ig))*myPhiVol  !tr
                            !--- total removal
                            !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+XsMac%isoXsMacA(iso,ig)  *myPhiVol  ! abs
                            DO ig2 = 1, ng
                                IF( ig2 .NE. ig )THEN
                                     !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+IsoMacXsSm(iso,ig,ig2) *myPhiVol  !removal
                                ENDIF 
                            ENDDO
                            isoMacXs(isoidx,4,ig)=isoMacXs(isoidx,4,ig)+XsMac%isoXsMacF(iso,ig)  *myPhiVol  !fission
                            isoMacXs(isoidx,5,ig)=isoMacXs(isoidx,5,ig)+XsMac%isoXsMacNF(iso,ig) *myPhiVol  !nu-fission
                            isoMacXs(isoidx,6,ig)=isoMacXs(isoidx,6,ig)+XsMac%isoXsMacKF(iso,ig) *myPhiVol  !kappa-fission
                            DO ig2 = 1, ng
                                isoMacSm(isoidx,ig,ig2)  =isoMacSm(isoidx,ig,ig2)  +IsoMacXsSm(iso,ig,ig2)  *myPhiVol
                                isoMacSmP1(isoidx,ig,ig2)=isoMacSmP1(isoidx,ig,ig2)+IsoMacXsSmP1(iso,ig,ig2)*myPhiVol
                            ENDDO
                            
                        ENDDO !---END of G sweep
                    ENDDO !---END of Fsr sweep
                ENDDO !---END of Iso sweep
            ELSE
                iso=1
                imix=myFxr%imix
                DO k = 1, nFsrInFxr
                    !ifsr = FsrIdxSt + k -1
                    !--- EDIT 14/12/23 Sq. cell problem
                    ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                    isoFSV=0._8
                    DO ig2 = 1, ng
                        isoFSV=isoFSV+MacXsBen(imix)%XsKF(ig2) *FsrPhiVol(ifsr,ig2)
                    ENDDO
                    FSV=FSV+isoFSV
                    DO ig = 1, ng
                        myPhiVol=FsrPhiVol(ifsr,ig)
                        
                        !--- MacroScopic XS ----------
                        !MAC : (total), D, a, r, nf, kf, kf                        
                        isoMacXs(0,0,ig)=isoMacXs(0,0,ig)+MacXsBen(imix)%XsT(ig) *myPhiVol   !Total
                        isoMacXs(0,1,ig)=isoMacXs(0,1,ig)+MacXsBen(imix)%XsTR(ig) *myPhiVol  !D
                        isoMacXs(0,2,ig)=isoMacXs(0,2,ig)+MacXsBen(imix)%Xsa(ig) *myPhiVol   !absorption                        
                        isoMacXs(0,3,ig)=isoMacXs(0,3,ig)+MacXsBen(imix)%Xsa(ig)  *myPhiVol   !removal base
                        
                        DO ig2 = 1, ng
                            IF( ig2 .NE. ig )THEN
                                 isoMacXs(0,3,ig)=isoMacXs(0,3,ig)+MacXsBen(imix)%XSS0(ig,ig2)  *myPhiVol   !removal
                            ENDIF 
                        ENDDO
                        isoMacXs(0,4,ig)=isoMacXs(0,4,ig)+MacXsBen(imix)%Xsnf(ig) *myPhiVol  !nu-fission
                        isoMacXs(0,5,ig)=isoMacXs(0,5,ig)+MacXsBen(imix)%Xskf(ig) *myPhiVol  !kappa-fission
                        IF( isoFSV .NE. 0)THEN
                            isoMacXs(0,6,ig)=isoMacXs(0,6,ig)+MacXsBen(imix)%chi(ig) *isoFSV           !CHI
                        ENDIF
                        DO ig2 = 1, ng
                            isoMacSm(0,ig,ig2)  =isoMacSm(0,ig,ig2)  +MacXsBen(imix)%XSS0(ig,ig2)  *myPhiVol !MacroScattering
                            !isoMacSmP1(0,ig,ig2)=isoMacSmP1(0,ig,ig2)+MacXsBen(imix)%XSS1(ig,ig2)  *myPhiVol !MacroScattering
                        ENDDO
                        isoidx=1
                        !--- MicroScopic XS ----------
                        !MIC : tr, a, r, f, nu, k
                        isoMacXs(isoidx,0,ig)=isoMacXs(isoidx,0,ig)+MacXsBen(imix)%XsT(ig) *myPhiVol  !Total
                        isoMacXs(isoidx,1,ig)=isoMacXs(isoidx,1,ig)+MacXsBen(imix)%XsTr(ig) *myPhiVol  !D                        
                        isoMacXs(isoidx,2,ig)=isoMacXs(isoidx,2,ig)+MacXsBen(imix)%Xsa(ig) *myPhiVol  !absorption                        
                        !--- transport removal
                        !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+(IsoMacXsSm(iso,ig,ig))*myPhiVol  !ss
                        !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+(XsMac%isoXsMacTR(iso,ig)-IsoMacXsSm(iso,ig,ig))*myPhiVol  !removal
                        !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+(XsMac%isoXsMacTR(iso,ig))*myPhiVol  !tr
                        !--- total removal
                        !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+XsMac%isoXsMacA(iso,ig)  *myPhiVol  ! abs
                        DO ig2 = 1, ng
                            IF( ig2 .NE. ig )THEN
                                 !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+IsoMacXsSm(iso,ig,ig2) *myPhiVol  !removal
                            ENDIF 
                        ENDDO
                        isoMacXs(isoidx,4,ig)=isoMacXs(isoidx,4,ig)+MacXsBen(imix)%Xsnf(ig)  *myPhiVol  !fission
                        isoMacXs(isoidx,5,ig)=isoMacXs(isoidx,5,ig)+MacXsBen(imix)%Xsnf(ig) *myPhiVol  !nu-fission
                        isoMacXs(isoidx,6,ig)=isoMacXs(isoidx,6,ig)+MacXsBen(imix)%Xskf(ig) *myPhiVol  !kappa-fission
                        DO ig2 = 1, ng
                            isoMacSm(isoidx,ig,ig2)  =isoMacSm(isoidx,ig,ig2)  +MacXsBen(imix)%XSS0(ig,ig2)  *myPhiVol
                            !isoMacSmP1(isoidx,ig,ig2)=isoMacSmP1(isoidx,ig,ig2)+MacXsBen(imix)%XSS1(ig,ig2)*myPhiVol
                        ENDDO
                        
                    ENDDO !---END of G sweep
                ENDDO !---END of Fsr sweep
            ENDIF            
            !DEALLOCATE(IsoMacXsSm)
            CALL FreeXsIsoMac(XsMac)
        ENDDO !---END of Fxr sweep
    ENDDO !---END of nPin sweep
    DO ig = 1, ng
        DO ixyl = 1, nxy 
            ixy = Asy(iasy)%GlobalPinIdx(ixyl) 
            Phi(ig)       =Phi(ig)+PinPhiVol(ig,ixy)  !PhiVol sum
        ENDDO
        Phi(ig)=Phi(ig)/Vol
    ENDDO
    
    iso=0; !homogenized Macro XS 47G iso=0
    DO ig = 1, ng
        DO i = 0, nxs-1
            isoMacXs(iso,i,ig)=isoMacXs(iso,i,ig)/PhiVol(ig)
        ENDDO
        XSt(ig)=isoMacXs(iso,0,ig)
        IF( FSV .NE. 0 )THEN
            isoMacXs(iso,6,ig)=isoMacXs(iso,6,ig)/FSV !chi
        ELSE
            isoMacXs(iso,6,ig)=0.0 !chi for non-fissile material
        ENDIF            
        DO ig2 = 1, ng
            isoMacSm(iso,ig,ig2)  =isoMacSm(iso,ig,ig2)  /PhiVol(ig)
            isoMacSmP1(iso,ig,ig2)=isoMacSmP1(iso,ig,ig2)/PhiVol(ig)
        ENDDO
    ENDDO
    
    
    !--- MicroScopic XS ----------
    DO iso = 1, nisotot
        DO ig = 1, ng
            DO i = 0, nxs
                isoMacXs(iso,i,ig)=isoMacXs(iso,i,ig)/PhiVol(ig)
            ENDDO
            DO ig2 = 1, ng
                isoMacSm(iso,ig,ig2)  =isoMacSm(iso,ig,ig2)  /PhiVol(ig)
                isoMacSmP1(iso,ig,ig2)=isoMacSmP1(iso,ig,ig2)/PhiVol(ig)
            ENDDO
        ENDDO
    ENDDO
    
ENDSUBROUTINE

SUBROUTINE GCHom_pin(Core, THInfo, FmInfo, GroupInfo, nTracerCntl, PE, ng)
    USE TYPEDEF,      ONLY : CoreInfo_Type, FmInfo_Type, CMInfo_Type, THInfo_Type, GroupInfo_Type,PE_TYPE, XsMac_Type
    USE CNTL,         ONLY : nTracerCntl_Type
    USE GC_mod
    USE XsUtil_mod, ONLY : FreeXsIsoMac
    USE BenchXs,      ONLY : MacXsBen
    USE TH_Mod,           ONLY : GetPinFuelTemp
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FMInfo_Type) :: FmInfo
    TYPE(THInfo_Type) :: THInfo
    TYPE(GroupInfo_Type) :: GroupInfo
    TYPE(nTracerCntl_Type) :: nTracerCntl
    TYPE(PE_TYPE) :: PE
    TYPE(XsMac_Type) :: XsMac
    INTEGER :: ng
    
    INTEGER :: i, j, k
    INTEGER :: ig, ig2
    
    INTEGER :: iso, isoidx, id, idiso, nisoinFxr, imix
    INTEGER :: ix, iy, ixy, ixyl
    INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, iCel, iFxr, ifsr, nFsrInFxr
    REAL :: myPhiVol, Vol
    REAL :: FSV, isoFSV ! Fission Source*Volume
    
INTERFACE
    SUBROUTINE GenFxrIsoMacXS(XsMac, IsoMacXsSm, IsoMacXsSmP1, Fxr, ifxr, Tempref, ig1, ig2, Core, ipin, iz, GroupInfo, nTRACERCntl, PE)
    USE PARAM
    USE TYPEDEF,      ONLY : CoreInfo_Type, XsMac_Type,  FxrInfo_Type,  GroupInfo_Type,  PE_Type
    USE CNTL,          ONLY : nTracerCntl_Type
    USE MacXsLib_mod
    USE BasicOperation, ONLY : CP_VA
    USE GC_mod, ONLY : isosize
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(XsMac_Type) :: XsMac
    TYPE(FxrInfo_Type),pointer :: Fxr(:,:)
    TYPE(GroupInfo_Type) :: GroupInfo
    TYPE(nTracerCntl_Type) :: nTracerCntl
    TYPE(PE_TYPE) :: PE
    INTEGER :: ig1, ig2, ipin, iz, ifxr
    REAL :: TempRef
    REAL :: IsoMacXsSm(isosize,GroupInfo%ng,GroupInfo%ng)
    REAL :: IsoMacXsSmP1(isosize,GroupInfo%ng,GroupInfo%ng)
    END SUBROUTINE
END INTERFACE
    
    FSV=0; !fission source *Volume
    Vol=0
    DO iy = yst, yed
    DO ix = xbg, xed
        ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy)
        ixy = Asy(iasy)%GlobalPinIdx(ixyl) 
        icel = Pin(ixy)%Cell(iz)
        Vol=Vol+Core%PinVol(ixy,iz)
        FxrIdxSt = Pin(ixy)%FxrIdxSt
        FsrIdxSt = Pin(ixy)%FsrIdxSt
        nlocalFxr = Cell(icel)%nFxr
        DO j = 1, nLocalFxr
            ifxr = FxrIdxSt + j -1
            myFxr => FmInfo%FXR(ifxr,iz)
            nFsrInFxr = myFxr%nFsrInFxr
            !FsrIdxSt = myFxr%FsrIdxSt
            nisoInFxr = myFxr%niso
            IF( nTracerCntl%lXsLib ) CALL GenFxrIsoMacXS(XsMac, IsoMacXsSm, IsoMacXsSmP1, FmInfo%FXR, ifxr, Tempref, 1, ng, Core, ixy, iz, GroupInfo, nTRACERCntl, PE)
            ! 13/11/01
            DO k = 1, nFsrInFxr
                !ifsr = FsrIdxSt + k -1
                !--- EDIT 14/12/23 Sq. cell problem
                ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                DO ig = 1, ng
                    PinPhiVol(ig,ixy)=PinPhiVol(ig,ixy)+FsrPhiVol(ifsr,ig)
                ENDDO
            ENDDO        
            ! 13/11/01 end
            IF( nTracerCntl%lXsLib )THEN
                DO iso = 1, nisoInFxr
                    idiso=myFxr%idiso(iso)
                    IF( idiso .EQ. 8001 )THEN
                        idiso=8016
                    ENDIF
                    id = MapNucl(idiso)  !id = Map(92235) = 37
                    !isodata => ldiso(id)
                    isoidx=isoList(id,1)
                    IF( idiso .EQ. 8016 )THEN
                        IF( myFxr%lFuel )THEN
                            isoidx=isoidx+1
                        ENDIF
                    ENDIF
                    DO k = 1, nFsrInFxr
                        !ifsr = FsrIdxSt + k -1
                        !--- EDIT 14/12/23 Sq. cell problem
                        ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                        isoFSV=0
                        DO ig2 = 1, ng
                            isoFSV=isoFSV+XsMac%isoXsMacKF(iso,ig2) *FsrPhiVol(ifsr,ig2)
                        ENDDO
                        FSV=FSV+isoFSV
                        DO ig = 1, ng
                            myPhiVol=FsrPhiVol(ifsr,ig)
                            
                            !--- MacroScopic XS ----------
                            !MAC : (total), D, a, r, nf, kf, kf                        
                            isoMacXs(0,0,ig)=isoMacXs(0,0,ig)+XsMac%isoXsMacT(iso,ig) *myPhiVol   !Total
                            isoMacXs(0,1,ig)=isoMacXs(0,1,ig)+XsMac%isoXsMacTR(iso,ig) *myPhiVol  !D
                            isoMacXs(0,2,ig)=isoMacXs(0,2,ig)+XsMac%isoXsMacA(iso,ig) *myPhiVol   !absorption                        
                            isoMacXs(0,3,ig)=isoMacXs(0,3,ig)+XsMac%isoXsMacA(iso,ig) *myPhiVol   !removal base
                            
                            DO ig2 = 1, ng
                                IF( ig2 .NE. ig )THEN
                                     isoMacXs(0,3,ig)=isoMacXs(0,3,ig)+IsoMacXsSm(iso,ig,ig2) *myPhiVol   !removal
                                ENDIF 
                            ENDDO
                            isoMacXs(0,4,ig)=isoMacXs(0,4,ig)+XsMac%isoXsMacNF(iso,ig) *myPhiVol  !nu-fission
                            isoMacXs(0,5,ig)=isoMacXs(0,5,ig)+XsMac%isoXsMacKF(iso,ig) *myPhiVol  !kappa-fission
                            IF( ig .LE. groupinfo%nchi .AND. isoFSV .NE. 0)THEN
                                isoMacXs(0,6,ig)=isoMacXs(0,6,ig)+myFXR%CHI(ig) *isoFSV           !CHI
                            ENDIF
                            DO ig2 = 1, ng
                                isoMacSm(0,ig,ig2)  =isoMacSm(0,ig,ig2)  +IsoMacXsSm(iso,ig,ig2)  *myPhiVol !MacroScattering
                                isoMacSmP1(0,ig,ig2)=isoMacSmP1(0,ig,ig2)+IsoMacXsSmP1(iso,ig,ig2)*myPhiVol
                            ENDDO
                            
                            !--- MicroScopic XS ----------
                            !MIC : tr, a, r, f, nu, k
                            isoMacXs(isoidx,0,ig)=isoMacXs(isoidx,0,ig)+XsMac%isoXsMacT(iso,ig) *myPhiVol  !Total
                            isoMacXs(isoidx,1,ig)=isoMacXs(isoidx,1,ig)+XsMac%isoXsMacTR(iso,ig) *myPhiVol  !D                        
                            isoMacXs(isoidx,2,ig)=isoMacXs(isoidx,2,ig)+XsMac%isoXsMacA(iso,ig)  *myPhiVol  !absorption                        
                            !--- transport removal
                            !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+(IsoMacXsSm(iso,ig,ig))*myPhiVol  !ss
                            !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+(XsMac%isoXsMacTR(iso,ig)-IsoMacXsSm(iso,ig,ig))*myPhiVol  !removal
                            !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+(XsMac%isoXsMacTR(iso,ig))*myPhiVol  !tr
                            !--- total removal
                            !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+XsMac%isoXsMacA(iso,ig)  *myPhiVol  ! abs
                            DO ig2 = 1, ng
                                IF( ig2 .NE. ig )THEN
                                     !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+IsoMacXsSm(iso,ig,ig2) *myPhiVol  !removal
                                ENDIF 
                            ENDDO
                            isoMacXs(isoidx,4,ig)=isoMacXs(isoidx,4,ig)+XsMac%isoXsMacF(iso,ig)  *myPhiVol  !fission
                            isoMacXs(isoidx,5,ig)=isoMacXs(isoidx,5,ig)+XsMac%isoXsMacNF(iso,ig) *myPhiVol  !nu-fission
                            isoMacXs(isoidx,6,ig)=isoMacXs(isoidx,6,ig)+XsMac%isoXsMacKF(iso,ig) *myPhiVol  !kappa-fission
                            DO ig2 = 1, ng
                                isoMacSm(isoidx,ig,ig2)  =isoMacSm(isoidx,ig,ig2)  +IsoMacXsSm(iso,ig,ig2)  *myPhiVol
                                isoMacSmP1(isoidx,ig,ig2)=isoMacSmP1(isoidx,ig,ig2)+IsoMacXsSmP1(iso,ig,ig2)*myPhiVol
                            ENDDO
                            
                        ENDDO !---END of G sweep
                    ENDDO !---END of Fsr sweep
                ENDDO !---END of Iso sweep
            ELSE
                iso=1
                imix=myFxr%imix
                DO k = 1, nFsrInFxr
                    !ifsr = FsrIdxSt + k -1
                    !--- EDIT 14/12/23 Sq. cell problem
                    ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                    isoFSV=0._8
                    DO ig2 = 1, ng
                        isoFSV=isoFSV+MacXsBen(imix)%XsKF(ig2) *FsrPhiVol(ifsr,ig2)
                    ENDDO
                    FSV=FSV+isoFSV
                    DO ig = 1, ng
                        myPhiVol=FsrPhiVol(ifsr,ig)
                        
                        !--- MacroScopic XS ----------
                        !MAC : (total), D, a, r, nf, kf, kf                        
                        isoMacXs(0,0,ig)=isoMacXs(0,0,ig)+MacXsBen(imix)%XsT(ig) *myPhiVol   !Total
                        isoMacXs(0,1,ig)=isoMacXs(0,1,ig)+MacXsBen(imix)%XsTR(ig) *myPhiVol  !D
                        isoMacXs(0,2,ig)=isoMacXs(0,2,ig)+MacXsBen(imix)%Xsa(ig) *myPhiVol   !absorption                        
                        isoMacXs(0,3,ig)=isoMacXs(0,3,ig)+MacXsBen(imix)%Xsa(ig)  *myPhiVol   !removal base
                        
                        DO ig2 = 1, ng
                            IF( ig2 .NE. ig )THEN
                                 isoMacXs(0,3,ig)=isoMacXs(0,3,ig)+MacXsBen(imix)%XSS0(ig,ig2)  *myPhiVol   !removal
                            ENDIF 
                        ENDDO
                        isoMacXs(0,4,ig)=isoMacXs(0,4,ig)+MacXsBen(imix)%Xsnf(ig) *myPhiVol  !nu-fission
                        isoMacXs(0,5,ig)=isoMacXs(0,5,ig)+MacXsBen(imix)%Xskf(ig) *myPhiVol  !kappa-fission
                        IF( isoFSV .NE. 0)THEN
                            isoMacXs(0,6,ig)=isoMacXs(0,6,ig)+MacXsBen(imix)%chi(ig) *isoFSV           !CHI
                        ENDIF
                        DO ig2 = 1, ng
                            isoMacSm(0,ig,ig2)  =isoMacSm(0,ig,ig2)  +MacXsBen(imix)%XSS0(ig,ig2)  *myPhiVol !MacroScattering
                            !isoMacSmP1(0,ig,ig2)=isoMacSmP1(0,ig,ig2)+MacXsBen(imix)%XSS1(ig,ig2)  *myPhiVol !MacroScattering
                        ENDDO
                        isoidx=1
                        !--- MicroScopic XS ----------
                        !MIC : tr, a, r, f, nu, k
                        isoMacXs(isoidx,0,ig)=isoMacXs(isoidx,0,ig)+MacXsBen(imix)%XsT(ig) *myPhiVol  !Total
                        isoMacXs(isoidx,1,ig)=isoMacXs(isoidx,1,ig)+MacXsBen(imix)%XsTr(ig) *myPhiVol  !D                        
                        isoMacXs(isoidx,2,ig)=isoMacXs(isoidx,2,ig)+MacXsBen(imix)%Xsa(ig) *myPhiVol  !absorption                        
                        !--- transport removal
                        !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+(IsoMacXsSm(iso,ig,ig))*myPhiVol  !ss
                        !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+(XsMac%isoXsMacTR(iso,ig)-IsoMacXsSm(iso,ig,ig))*myPhiVol  !removal
                        !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+(XsMac%isoXsMacTR(iso,ig))*myPhiVol  !tr
                        !--- total removal
                        !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+XsMac%isoXsMacA(iso,ig)  *myPhiVol  ! abs
                        DO ig2 = 1, ng
                            IF( ig2 .NE. ig )THEN
                                 !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+IsoMacXsSm(iso,ig,ig2) *myPhiVol  !removal
                            ENDIF 
                        ENDDO
                        isoMacXs(isoidx,4,ig)=isoMacXs(isoidx,4,ig)+MacXsBen(imix)%Xsnf(ig)  *myPhiVol  !fission
                        isoMacXs(isoidx,5,ig)=isoMacXs(isoidx,5,ig)+MacXsBen(imix)%Xsnf(ig) *myPhiVol  !nu-fission
                        isoMacXs(isoidx,6,ig)=isoMacXs(isoidx,6,ig)+MacXsBen(imix)%Xskf(ig) *myPhiVol  !kappa-fission
                        DO ig2 = 1, ng
                            isoMacSm(isoidx,ig,ig2)  =isoMacSm(isoidx,ig,ig2)  +MacXsBen(imix)%XSS0(ig,ig2)  *myPhiVol
                            !isoMacSmP1(isoidx,ig,ig2)=isoMacSmP1(isoidx,ig,ig2)+MacXsBen(imix)%XSS1(ig,ig2)*myPhiVol
                        ENDDO
                        
                    ENDDO !---END of G sweep
                ENDDO !---END of Fsr sweep
            ENDIF            
            !DEALLOCATE(IsoMacXsSm)
            CALL FreeXsIsoMac(XsMac)
        ENDDO !---END of Fxr sweep
    ENDDO    
    ENDDO !---END of nPin sweep
    
    DO ig = 1, ng
        DO iy = yst, yed
        DO ix = xbg, xed
            ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy)
            ixy = Asy(iasy)%GlobalPinIdx(ixyl) 
            Phi(ig)       =Phi(ig)+PinPhiVol(ig,ixy)  !PhiVol sum
        ENDDO
        ENDDO
        Phi(ig)=Phi(ig)/Vol
    ENDDO
    
    iso=0; !homogenized Macro XS 47G iso=0
    DO ig = 1, ng
        DO i = 0, nxs-1
            isoMacXs(iso,i,ig)=isoMacXs(iso,i,ig)/PhiVol(ig)
        ENDDO
        XSt(ig)=isoMacXs(iso,0,ig)
        IF( FSV .NE. 0 )THEN
            isoMacXs(iso,6,ig)=isoMacXs(iso,6,ig)/FSV !chi
        ELSE
            isoMacXs(iso,6,ig)=0.0 !chi for non-fissile material
        ENDIF            
        DO ig2 = 1, ng
            isoMacSm(iso,ig,ig2)  =isoMacSm(iso,ig,ig2)  /PhiVol(ig)
            isoMacSmP1(iso,ig,ig2)=isoMacSmP1(iso,ig,ig2)/PhiVol(ig)
        ENDDO
    ENDDO
    
    
    !--- MicroScopic XS ----------
    DO iso = 1, nisotot
        DO ig = 1, ng
            DO i = 0, nxs
                isoMacXs(iso,i,ig)=isoMacXs(iso,i,ig)/PhiVol(ig)
            ENDDO
            DO ig2 = 1, ng
                isoMacSm(iso,ig,ig2)  =isoMacSm(iso,ig,ig2)  /PhiVol(ig)
                isoMacSmP1(iso,ig,ig2)=isoMacSmP1(iso,ig,ig2)/PhiVol(ig)
            ENDDO
        ENDDO
    ENDDO
    
ENDSUBROUTINE

SUBROUTINE GenU238N2N(xsn2n, temp, ig)
USE XsUtil_mod,      ONLY : XsTempInterpolation
USE XSLIB_MOD
IMPLICIT NONE
INTEGER :: niso, ng
REAL :: xsn2n
REAL :: temp
INTEGER :: idiso = 92238
TYPE(LIBDATA), POINTER :: isodata
INTEGER :: ig, iso, id, idn2n, it1, it2
REAL :: wt1, wt2, phisum                       !Temperature interpolation

id = mapnucl(idiso);   isodata => ldiso(id)
idn2n=mapn2n(id);
!Temperature Interpolation
CALL XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2) 
IF(libtyp == 0) THEN
    IF(idn2n .NE. 0) xsn2n = xsn2nhel(idn2n,ig)
ELSE
    IF(idn2n .NE. 0) xsn2n = isodata%sign2n(ig)
ENDIF
NULLIFY(isodata) !Free the pointing variable
END SUBROUTINE

SUBROUTINE GCHom_residual(Core, THInfo, FmInfo, GroupInfo, nTracerCntl, PE, ng)
    USE TYPEDEF,      ONLY : CoreInfo_Type, FmInfo_Type, CMInfo_Type, THInfo_Type, GroupInfo_Type,PE_TYPE, XsMac_Type
    USE CNTL,         ONLY : nTracerCntl_Type
    USE GC_mod
    USE GCpin_mod
    USE TH_Mod,           ONLY : GetPinFuelTemp
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FMInfo_Type) :: FmInfo
    TYPE(THInfo_Type) :: THInfo
    TYPE(GroupInfo_Type) :: GroupInfo
    TYPE(nTracerCntl_Type) :: nTracerCntl
    TYPE(PE_TYPE) :: PE
    TYPE(XsMac_Type) :: XsMac
    INTEGER :: ng
    
    INTEGER :: i, j, k 
    INTEGER :: ig, ig2
    
    INTEGER :: iso, isoidx, id, idiso, nisoinFxr
    INTEGER :: ixy, ixyl
    INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, iCel, iFxr, ifsr, nFsrInFxr
    REAL :: myPhiVol, Vol
    REAL :: FSV, isoFSV ! Fission Source*Volume
    
INTERFACE
    SUBROUTINE GenFxrIsoMacXS(XsMac, IsoMacXsSm, IsoMacXsSmP1, Fxr, ifxr, Tempref, ig1, ig2, Core, ipin, iz, GroupInfo, nTRACERCntl, PE)
    USE PARAM
    USE TYPEDEF,      ONLY : CoreInfo_Type, XsMac_Type,  FxrInfo_Type,  GroupInfo_Type,  PE_Type
    USE CNTL,          ONLY : nTracerCntl_Type
    USE MacXsLib_mod
    USE BasicOperation, ONLY : CP_VA
    USE GC_mod, ONLY : isosize
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(XsMac_Type) :: XsMac
    TYPE(FxrInfo_Type),pointer :: Fxr(:,:)
    TYPE(GroupInfo_Type) :: GroupInfo
    TYPE(nTracerCntl_Type) :: nTracerCntl
    TYPE(PE_TYPE) :: PE
    INTEGER :: ig1, ig2, ipin, iz, ifxr
    REAL :: TempRef
    REAL :: IsoMacXsSm(isosize,GroupInfo%ng,GroupInfo%ng)
    REAL :: IsoMacXsSmP1(isosize,GroupInfo%ng,GroupInfo%ng)
    END SUBROUTINE
END INTERFACE
    !ALLOCATE(pinXsr(nxy,ng),pinXSnf(nxy,ng),pinXSss(nxy,ng,ng),pinChi(nxy,ng))
    !pinXsr=0;pinXsnf=0;pinXsss=0;pinchi=0;
    
    
    FSV=0; !fission source *Volume
    Vol=0
    DO ixyl = 1, nxy 
        ixy = Asy(iasy)%GlobalPinIdx(ixyl) 
        icel = Pin(ixy)%Cell(iz)
        Vol=Vol+Core%PinVol(ixy,iz)
        FxrIdxSt = Pin(ixy)%FxrIdxSt
        FsrIdxSt = Pin(ixy)%FsrIdxSt
        nlocalFxr = Cell(icel)%nFxr
        Tempref = GetPinFuelTemp(Core, FmInfo%FXR, iz, ixy)
        DO j = 1, nLocalFxr
            ifxr = FxrIdxSt + j -1
            myFxr => FmInfo%FXR(ifxr,iz)
            nFsrInFxr = myFxr%nFsrInFxr
            !FsrIdxSt = myFxr%FsrIdxSt
            nisoInFxr = myFxr%niso
            CALL GenFxrIsoMacXS(XsMac, IsoMacXsSm, IsoMacXsSmP1, FmInfo%FXR, ifxr, Tempref, 1, ng, Core, ixy, iz, GroupInfo, nTRACERCntl, PE)
            ! 13/11/01
            DO k = 1, nFsrInFxr
                !ifsr = FsrIdxSt + k -1
                !--- EDIT 14/12/23 Sq. cell problem
                ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                DO ig = 1, ng
                    PinPhiVol(ig,ixy)=PinPhiVol(ig,ixy)+FsrPhiVol(ifsr,ig)
                    pinFlux(ixy,ig)=pinFlux(ixy,ig)+FsrPhiVol(ifsr,ig)
                ENDDO
            ENDDO        
            ! 13/11/01 end
            DO iso = 1, nisoInFxr
                idiso=myFxr%idiso(iso)
                IF( idiso .EQ. 8001 )THEN
                    idiso=8016
                ENDIF
                id = MapNucl(idiso)  !id = Map(92235) = 37
                !isodata => ldiso(id)
                isoidx=isoList(id,1)
                IF( idiso .EQ. 8016 )THEN
                    IF( myFxr%lFuel )THEN
                        isoidx=isoidx+1
                    ENDIF
                ENDIF
                DO k = 1, nFsrInFxr
                    !ifsr = FsrIdxSt + k -1
                    !--- EDIT 14/12/23 Sq. cell problem
                    ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                    isoFSV=0
                    DO ig2 = 1, ng
                        !isoFSV=isoFSV+XsMac%isoXsMacKF(iso,ig2) *FsrPhiVol(ifsr,ig2)
                        isoFSV=isoFSV+XsMac%isoXsMacNF(iso,ig2) *FsrPhiVol(ifsr,ig2)
                    ENDDO
                    FSV=FSV+isoFSV
                    DO ig = 1, ng
                        myPhiVol=FsrPhiVol(ifsr,ig)
                        
                        !--- MacroScopic XS ----------
                        !MAC : (total), D, a, r, nf, kf, kf                        
                        isoMacXs(0,0,ig)=isoMacXs(0,0,ig)+XsMac%isoXsMacT(iso,ig) *myPhiVol   !Total
                        isoMacXs(0,1,ig)=isoMacXs(0,1,ig)+XsMac%isoXsMacTR(iso,ig) *myPhiVol  !D
                        isoMacXs(0,2,ig)=isoMacXs(0,2,ig)+XsMac%isoXsMacA(iso,ig) *myPhiVol   !absorption                        
                        isoMacXs(0,3,ig)=isoMacXs(0,3,ig)+XsMac%isoXsMacA(iso,ig) *myPhiVol   !removal base
                        pinRmv(ixy,ig)=pinRmv(ixy,ig)+XsMac%isoXsMacA(iso,ig) *myPhiVol   !removal base
                        pinAbs(ixy,ig)=pinAbs(ixy,ig)+XsMac%isoXsMacA(iso,ig) *myPhiVol   !removal base
                        DO ig2 = 1, ng
                            IF( ig2 .NE. ig )THEN
                                 isoMacXs(0,3,ig)=isoMacXs(0,3,ig)+IsoMacXsSm(iso,ig,ig2) *myPhiVol   !removal
                            ENDIF 
                            pinRmv(ixy,ig)=pinRmv(ixy,ig)+IsoMacXsSm(iso,ig,ig2) *myPhiVol   !total
                        ENDDO
                        isoMacXs(0,4,ig)=isoMacXs(0,4,ig)+XsMac%isoXsMacNF(iso,ig) *myPhiVol  !nu-fission
                        !pinFis(ixyl,ig)=pinFis(ixyl,ig)+XsMac%isoXsMacNF(iso,ig) *myPhiVol  !nu-fission
                        
                        isoMacXs(0,5,ig)=isoMacXs(0,5,ig)+XsMac%isoXsMacKF(iso,ig) *myPhiVol  !kappa-fission
                        IF( ig .LE. groupinfo%nchi .AND. isoFSV .NE. 0)THEN
                            isoMacXs(0,6,ig)=isoMacXs(0,6,ig)+myFXR%CHI(ig) *isoFSV           !CHI
                            pinFis(ixy,ig)=pinFis(ixy,ig)+myFXR%CHI(ig) *isoFSV           !CHI
                        ENDIF
                        DO ig2 = 1, ng
                            isoMacSm(0,ig,ig2)=isoMacSm(0,ig,ig2)+IsoMacXsSm(iso,ig,ig2)*myPhiVol !MacroScattering
                            pinSS(ixy,ig,ig2)=pinSS(ixy,ig,ig2)+IsoMacXsSm(iso,ig,ig2)*myPhiVol !MacroScattering
                            !Sm1(ig,ig2)=Sm1(ig,ig2)+IsoMacXsSm1(iso,ig,ig2)*myPhiVol
                        ENDDO
                        
                        !--- MicroScopic XS ----------
                        !MIC : tr, a, r, f, nu, k
                        isoMacXs(isoidx,0,ig)=isoMacXs(isoidx,0,ig)+XsMac%isoXsMacT(iso,ig) *myPhiVol  !Total
                        isoMacXs(isoidx,1,ig)=isoMacXs(isoidx,1,ig)+XsMac%isoXsMacTR(iso,ig) *myPhiVol  !D                        
                        isoMacXs(isoidx,2,ig)=isoMacXs(isoidx,2,ig)+XsMac%isoXsMacA(iso,ig)  *myPhiVol  !absorption                        
                        !--- transport removal
                        !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+(IsoMacXsSm(iso,ig,ig))*myPhiVol  !ss
                        !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+(XsMac%isoXsMacTR(iso,ig)-IsoMacXsSm(iso,ig,ig))*myPhiVol  !removal
                        !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+(XsMac%isoXsMacTR(iso,ig))*myPhiVol  !tr
                        !--- total removal
                        !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+XsMac%isoXsMacA(iso,ig)  *myPhiVol  ! abs
                        DO ig2 = 1, ng
                            IF( ig2 .NE. ig )THEN
                                 !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+IsoMacXsSm(iso,ig,ig2) *myPhiVol  !removal
                            ENDIF 
                        ENDDO
                        isoMacXs(isoidx,4,ig)=isoMacXs(isoidx,4,ig)+XsMac%isoXsMacF(iso,ig)  *myPhiVol  !fission
                        isoMacXs(isoidx,5,ig)=isoMacXs(isoidx,5,ig)+XsMac%isoXsMacNF(iso,ig) *myPhiVol  !nu-fission
                        isoMacXs(isoidx,6,ig)=isoMacXs(isoidx,6,ig)+XsMac%isoXsMacKF(iso,ig) *myPhiVol  !kappa-fission
                        DO ig2 = 1, ng
                            isoMacSm(isoidx,ig,ig2)=isoMacSm(isoidx,ig,ig2)+IsoMacXsSm(iso,ig,ig2)*myPhiVol
                        ENDDO
                        
                    ENDDO !---END of G sweep
                ENDDO !---END of Fsr sweep
            ENDDO !---END of Iso sweep
            !DEALLOCATE(IsoMacXsSm)
        ENDDO !---END of Fxr sweep
        
    ENDDO !---END of nPin sweep
    
    DO ig = 1, ng
        DO ixyl = 1, nxy 
            ixy = Asy(iasy)%GlobalPinIdx(ixyl) 
            Phi(ig)       =Phi(ig)+PinPhiVol(ig,ixy)  !PhiVol sum
        ENDDO
        Phi(ig)=Phi(ig)/Vol
    ENDDO
    
    iso=0; !homogenized Macro XS 47G iso=0
    DO ig = 1, ng
        DO i = 0, nxs-1
            isoMacXs(iso,i,ig)=isoMacXs(iso,i,ig)/PhiVol(ig)
        ENDDO
        XSt(ig)=isoMacXs(iso,0,ig)
        IF( FSV .NE. 0 )THEN
            isoMacXs(iso,6,ig)=isoMacXs(iso,6,ig)/FSV !chi
        ELSE
            isoMacXs(iso,6,ig)=0.0 !chi for non-fissile material
        ENDIF            
        DO ig2 = 1, ng
            isoMacSm(iso,ig,ig2)=isoMacSm(iso,ig,ig2)/PhiVol(ig)
            !Sm1(ig,ig2)=Sm1(ig,ig2)/PhiVol(ig)
        ENDDO
    ENDDO
    
    
    !--- MicroScopic XS ----------
    DO iso = 1, nisotot
        DO ig = 1, ng
            DO i = 0, nxs
                isoMacXs(iso,i,ig)=isoMacXs(iso,i,ig)/PhiVol(ig)
            ENDDO
            DO ig2 = 1, ng
                isoMacSm(iso,ig,ig2)=isoMacSm(iso,ig,ig2)/PhiVol(ig)
            ENDDO
        ENDDO
    ENDDO
    
ENDSUBROUTINE


SUBROUTINE GCB1D(ng)
    USE GC_mod
    IMPLICIT NONE
    INTEGER :: ng
    
    INTEGER :: iso
    INTEGER :: ig
    
    Xstr=0
    DO ig = 1, ng        
        DO iso = 1, nisotot
           Xstr(ig)= Xstr(ig)+isoMacXs(iso,1,ig)
        ENDDO
        IF( Dng(ig) .NE. 0 )THEN !---if infinite spectrum, no correction
            XstrB1(ig)=1.0/(3*Dng(ig))
        ELSE
            Dng(ig)=1.0/(3*Xstr(ig))
            XstrB1(ig)=Xstr(ig)
        ENDIF
        DO iso = 0, nisotot
            isoMacXs(iso,1,ig) = XstrB1(ig) * isoMacXs(iso,1,ig)/Xstr(ig)
        ENDDO
    ENDDO
ENDSUBROUTINE


SUBROUTINE GCMac2Mic(ng)
    USE GC_mod
    IMPLICIT NONE
    INTEGER :: ng
    
    INTEGER :: iso
    INTEGER :: ig, ig2
    
    DO ig = 1, ng
        DO iso = 1, nisotot
            !MAC : D, a, r, f, nf, kf
            !MIC : tr, a, r, f, nu, k
            IF( isoNumden(iso) .EQ. 0 )THEN 
                isoMicXs(iso,1,ig)=0._8      ! tr
                isoMicXs(iso,2,ig)=0._8      ! absorption
                isoMicXs(iso,4,ig)=0._8      ! fission
                isoMicXs(iso,5,ig)=0  ! nu
                isoMicXs(iso,6,ig)=0  ! kappa
                DO ig2 = 1, ng
                    isoMicSm(iso,ig,ig2)=0._8
                ENDDO
            ELSE
                isoMicXs(iso,1,ig)=isoMacXs(iso,1,ig)/isoNumden(iso)      ! tr
                isoMicXs(iso,2,ig)=isoMacXs(iso,2,ig)/isoNumden(iso)      ! absorption
                !isoMicXs(iso,3,ig)=isoMacXs(iso,3,ig)/isoNumden(iso)      ! removal
                isoMicXs(iso,4,ig)=isoMacXs(iso,4,ig)/isoNumden(iso)      ! fission
                IF( isoMacXs(iso,4,ig) .NE. 0 )THEN
                    isoMicXs(iso,5,ig)=isoMacXs(iso,5,ig)/isoMacXs(iso,4,ig)  ! nu
                    isoMicXs(iso,6,ig)=isoMacXs(iso,6,ig)/isoMacXs(iso,4,ig)  ! kappa
                ELSE
                    isoMicXs(iso,5,ig)=0  ! nu
                    isoMicXs(iso,6,ig)=0  ! kappa
                ENDIF
                DO ig2 = 1, ng
                    isoMicSm(iso,ig,ig2)=isoMacSm(iso,ig,ig2)/isoNumden(iso)
                ENDDO
            ENDIF
        ENDDO
    ENDDO    
ENDSUBROUTINE
