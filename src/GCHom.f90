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
  
#define H2H
SUBROUTINE GCHom_pin(Core, THInfo, FmInfo, GroupInfo, nTracerCntl, PE, ng)
    USE TYPEDEF,      ONLY : CoreInfo_Type, FmInfo_Type, CMInfo_Type, THInfo_Type, GroupInfo_Type,PE_TYPE, XsMac_Type
    USE CNTL,         ONLY : nTracerCntl_Type
    USE GC_mod
    USE GCpin_mod
    USE XsUtil_mod, ONLY : FreeXsIsoMac
    USE BenchXs,      ONLY : MacXsBen
    USE TH_Mod,           ONLY : GetPinFuelTemp
    USE TRAN_MOD, ONLY : TranInfo
#ifdef H2H    
    USE DEPL_MOD,         ONLY : FxrBurnUp, ConstDeplVas, MakeDeplXs1g, fluxnormalizefactor,DeplCntl
    USE MOC_Mod,          ONLY : FxrAvgPhi
#endif    
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
    INTEGER :: ix, iy, ixy, ixyl, iprec
    INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, iCel, iFxr, ifsr, nFsrInFxr
    REAL :: myPhiVol, Vol, areasum
    REAL :: FSV, isoFSV, FSV_nu, isoFSV_nu ! Fission Source*Volume
    REAL :: locFSV_nu(0:6), FSV_nu_Adj, FSV_chi(ng)
#ifdef H2H
    REAL :: BUTime
    REAL :: normFactor
    INTEGER :: ngDep, nIsoLib, nIsoDepl, nFxrPin, ihetfxr, isogd
    REAL(8), ALLOCATABLE :: BurnupXS(:,:), avgphi(:), pinND(:), refND(:)
    REAL(8) :: fxrA, pinA, numinfxr
    LOGICAL :: lDeplPin
    INTEGER :: GdIsoIdx(7) = (/64152, 64154, 64155, 64156, 64157, 64158, 64160/)
#endif
    
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
    
#ifdef H2H
    nFxrPin = 0; ngDep = GroupInfo%ng; nIsoDepl = DeplLibPin%nIsoDep
    DO iy = yst, yed
      DO  ix = xbg, xed
        ixyl = AsyInfo(IasyType)%Pin2DIdx(ix,iy)
        ixy = Asy(iasy)%GlobalPinIdx(ixyl)
        icel = Pin(ixy)%Cell(iz)
        nlocalFxr = Cell(icel)%nFxr
        nFxrPin = nFxrPin+nlocalFxr
      END DO
    END DO
    hetNreg = nFxrPin
    ALLOCATE(hetRxFrac(nFxrPin,7), hetNumFrac(nFxrPin,7))
    hetRxFrac = 0.; hetNumFrac = 0.;
    IF (DeplCntl%lInitDepl) THEN
      ngDep = GroupInfo%ng; BUTime = DeplCntl%Tsec
      nIsoLib = GroupInfo%ntiso; nIsoDepl = DeplLibPin%nIsoDep
      NormFactor = FluxNormalizeFactor(Core, FmInfo, GroupInfo, DeplCntl%PowerCore, nTracerCntl%lCritSpec, TRUE, PE)
      NormFactor = NormFactor*nTracerCntl%PowerLevel
      ALLOCATE(BurnupXS(4,nIsoDepl),avgphi(ngDep),pinND(nIsoDepl),refND(nIsoDepl))
      DO iy = yst, yed
        DO ix = xbg, xed
          ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy)
          ixy = Asy(iasy)%GlobalPinIdx(ixyl) 
          icel = Pin(ixy)%Cell(iz)
          Vol=Vol+Core%PinVol(ixy,iz)
          FxrIdxSt = Pin(ixy)%FxrIdxSt
          FsrIdxSt = Pin(ixy)%FsrIdxSt
          nlocalFxr = Cell(icel)%nFxr
          BurnupXS = 0.; avgphi = 0.; pinND = 0.; refND = 0.; pinA = 0.; lDeplPin = .FALSE.
          DO j = 1, nLocalFxr
            ifxr = FxrIdxSt + j -1
            myFxr => FmInfo%FXR(ifxr,iz)
            fxrA = myFxr%area; pinA = pinA+fxrA
            DeplXSPin%AvgPhi = FxrAvgPhi(Core, FmInfo%Fxr, FmInfo%Phis, ixy, j, iz, ngDep, PE)
            avgphi(:) = avgphi(:)+fxrA*DeplXsPin%AvgPhi(:)
            IF (.NOT. myFxr%lDepl) CYCLE
            lDeplPin = .TRUE.
            nFsrInFxr = myFxr%nFsrInFxr
            CALL MakeDeplXs1g(FmInfo%Fxr, ifxr, ngDep, DeplXsPin, NormFactor, Core, ixy, iz, GroupInfo, PE, DeplCntl)
            CALL ConstDeplVas(DeplVarPin, DeplXsPin, myFxr, nIsoLib, nIsoDepl)
            DO  k = 1, 4
            BurnupXS(k,:) = BurnupXS(k,:) + fxrA*DeplVarPin%BurnupXs(k,:)*DeplXSPin%phi1g*DeplVarPin%IsoNum(:)
            END DO
            pinND(:) = pinND(:)+DeplVarPin%IsoNUM(:)*fxrA
            CALL FxrBurnUp(DeplVarPin, DeplLibPin, DeplCntl)
            refND(:) = refND(:)+DeplVarPin%IsoNum(:)*fxrA
          END DO
          IF (.NOT. lDeplPin) CYCLE
          avgphi(:) = avgphi(:)/pinA; DeplXSPin%phi1g = SUM(avgphi)*NormFactor
          DO k = 1, nIsoDepl
            pinND(k) = pinND(k)/pinA; refND(k) = refND(k)/pinA; DeplVarPin%IsoNum(k) = pinND(k)
            IF (pinND(k) .LT. 1.e-30) CYCLE
            BurnupXS(:,k) = BurnupXS(:,k)/DeplXSPin%phi1g/pinA/pinND(k)
            DeplVarPin%BurnupXS(:,k) = BurnupXS(:,k)
          END DO
          CALL FxrBurnUp(DeplVarPin, DeplLibPin, DeplCntl)
          pinND(:) = DeplVarPin%IsoNum(:)
          DO k = 1, nIsoDepl
            IF (pinND(k) .LT. 1.e-30) CYCLE
            h2hfactorpin(k) = LOG(refND(k)/pinND(k))/BUTime
          END DO
        END DO
      END DO
      DEALLOCATE(BurnupXS, avgphi, pinND, refND)
    END IF
#endif   
#ifdef H2H    
    ihetfxr = 0
#endif    
    FSV=0; !fission source *Volume
    Vol=0
    areasum = 0
#ifdef H2H    
    ihetfxr = 0
#endif

    KinParMac_velo = 0._8
    TranInfo%ChiD = (/ 4*0._8, 0.005_8, 0.021_8, 0.269_8, 0.247_8, 0.429_8, 0.029_8, 37*0._8/) 
    FSV=0; FSV_nu=0._8!fission source *Volume
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
#ifdef H2H          
            ihetfxr = ihetfxr+1
#endif            
            ifxr = FxrIdxSt + j -1
            myFxr => FmInfo%FXR(ifxr,iz)
            nFsrInFxr = myFxr%nFsrInFxr
            IF( Fminfo%FXR(ifxr,iz)%lfuel )THEN
                areasum=areasum+Fminfo%FXR(ifxr,iz)%area
                bupin=bupin+Fminfo%FXR(ifxr,iz)%burnup*Fminfo%FXR(ifxr,iz)%area
            ENDIF
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
            DO k = 1, nFsrInFxr
                ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                DO ig = 1, ng
                    U238phi=U238phi+      FsrPhiVol(ifsr,ig)
                    u238phiMG(ig) = u238phiMG(ig) + FsrPhiVol(ifsr,ig)
                ENDDO !---END of G sweep
            ENDDO !---END of Fsr sweep
#ifdef H2H
            IF (nTracerCntl%lXsLib) THEN
              DO iso = 1, nisoInFxr
                idiso = myFxr%idiso(iso)
                isogd = 0;
                DO  k = 1, 7
                  IF (idiso.EQ. gdisoidx(k)) isogd = k
                END DO
                IF (isogd.GT.0) THEN
                  hetNumFrac(ihetfxr,isogd) = hz*myFxr%area*myFxr%pnum(iso)
                END IF
              END DO
            END IF
#endif
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
                        isoFSV=0; isoFSV_nu=0._8
                        DO ig2 = 1, ng
                            isoFSV=isoFSV+XsMac%isoXsMacKF(iso,ig2) *FsrPhiVol(ifsr,ig2)
                            isoFSV_nu=isoFSV_nu+XsMac%isoXsMacNF(iso,ig2)*FsrPhiVol(ifsr,ig2)
                        ENDDO
                        FSV=FSV+isoFSV; FSV_nu=FSV_nu+isoFSV_nu
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
                            IF (iso .EQ. 1 .AND. nTracerCntl%lTranON) KinParMac_velo(0,ig) = KinParMac_velo(0,ig) + 1._8/myFXR%velo(ig) * myPhiVol
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
                    IF( idiso.EQ.92238 )THEN !--- U238 N2N
                        DO k = 1, nFsrInFxr
                            ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                            DO ig = 1, ng
                                CALL GenU238N2N(myn2n, TempRef, ig)
                                U238n2n=U238n2n+myn2n*myfxr%pnum(iso)*FsrPhiVol(ifsr,ig)
                                u238n2nMG(ig) = u238n2nMG(ig) + myn2n*myfxr%pnum(iso)*FsrPhiVol(ifsr,ig)
                            ENDDO !---END of G sweep
                        ENDDO !---END of Fsr sweep
                    ENDIF ! END OF U238 N2N
#ifdef H2H                    
                    isogd = 0;
                    DO k = 1, 7
                      IF (GdIsoIdx(k) .EQ. idiso) isogd = k;
                    END DO
                    IF (isogd.GT.0) THEN
                      DO ig = 1, ng
                        hetRxFrac(ihetfxr, isogd) = hetRxFrac(ihetfxr,isogd)+isoMacXS(isoidx,2,ig)
                      END DO
                      !hetRxFrac(ihetfxr,isogd) = hetRxFrac(ihetfxr,isogd)/hetNumFrac(ihetfxr,isogd)
                    END IF
#endif                    
                ENDDO !---END of Iso sweep
            ELSE
                iso=1
                imix=myFxr%imix
                DO k = 1, nFsrInFxr
                    !ifsr = FsrIdxSt + k -1
                    !--- EDIT 14/12/23 Sq. cell problem
                    ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                    isoFSV=0._8; isoFSV_nu=0._8
                    DO ig2 = 1, ng
                        isoFSV=isoFSV+MacXsBen(imix)%XsKF(ig2) *FsrPhiVol(ifsr,ig2)
                        isoFSV_nu=isoFSV_nu+MacXsBen(imix)%Xsnf(ig2)*FsrPhiVol(ifsr,ig2)
                    ENDDO
                    FSV=FSV+isoFSV; FSV_nu=FSV_nu+isoFSV_nu
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

    IF (nTracerCntl%lTranON) THEN
      FSV_nu=0._8
      FSV_nu_Adj = 0._8
      KinParMac_beta = 0._8
      KinParMac_betaeff = 0._8
      DO iy=yst, yed
        DO ix=xbg, xed
          ixyl = AsyInfo(iasytype)%Pin2DIdx(ix,iy)
          ixy = Asy(iasy)%GlobalPinIdx(ixyl)
          icel = Pin(ixy)%Cell(iz)
          Vol = Vol + Core%PinVol(ixy,iz)
          FxrIdxSt = Pin(ixy)%FxrIdxSt
          nlocalFxr = Cell(icel)%nFxr
          locFSV_nu = 0._8; FSV_chi = 0._8
          DO j=1, nLocalFxr
            ifxr = FxrIdxSt + j -1
            myFxr => FmInfo%FXR(ifxr,iz)
            nFsrInFxr = myFxr%nFsrInFxr
            nisoInFxr = myFxr%niso
            CALL GenFxrIsoMacXS(XsMac, IsoMacXsSM, IsoMacXsSMP1, FMInfo%FXR, ifxr, Tempref, 1, ng, Core, ixy, iz, GroupInfo, nTracerCntl, PE)
            isoFSV_nu = 0._8
            DO iso=1, nisoInFxr
              DO k=1, nFsrInFxr
                DO ig2=1, ng
                  isoFSV_nu = isoFSV_nu + XsMac%isoXsMacNF(iso,ig2)*FsrPhiVol(ifsr,ig2)
                END DO
              END DO
            END DO
            FSV_nu = FSV_nu + isoFSV_nu
            IF (myFxr%lFuel) THEN
              DO ig=1, ng
                FSV_chi(ig) = FSV_chi(ig) + myFxr%Chi(ig) * isoFSV_nu
              END DO
            END IF
            DO iprec=1, 6
              KinParMac_beta(0,iprec) = KinParMac_beta(0,iprec) + myFxr%beta(iprec) * isoFSV_nu
              locFSV_nu(iprec) = locFSV_nu(iprec) + myFxr%beta(iprec) * isoFSV_nu
            END DO
          END DO
          DO iprec=1, 6
            DO ig=1, ng
              !KinParMac_betaeff(0,iprec) = KinParMac_betaeff(0,iprec) + TranInfo%ChiD(ig) * CMInfo%Phic_Adj(ixy,iz,ig) * locFSV_nu(iprec)
            END DO
          END DO
          DO ig=1, ng
            !FSV_nu_Adj = FSV_nu_Adj + CMInfo%Phic_Adj(ixy,iz,ig) * FSV_chi(ig)
          END DO
        END DO
      END DO
      iso=0;
      DO iprec = 1, 6
        DO ig = 1, ng 
          KinParMac_ChiDg(iso, ig, iprec) = TranInfo%chidk(ig, iprec)
        END DO 
      END DO 
      KinParMac_fisrate(iso) = FSV_nu
      IF (FSV_nu .NE. 0._8) KinParMac_beta(iso,:)=KinParMac_beta(iso,:)/FSV_nu
      IF (FSV_nu_Adj .NE. 0._8) KinParMac_betaeff(iso,:)=KinParMac_betaeff(iso,:)/FSV_nu_Adj
    END IF

#ifdef H2H
    DO j = 1, HetNreg-1
      hetRxFrac(HetNreg-j+1,:) = hetRxFrac(HetNreg-j+1,:)-hetRxFrac(HetNreg-j,:)
      hetRxFrac(HetNreg-j+1,:) = hetRxFrac(HetNreg-j+1,:)/hetNumFrac(HetNreg-j+1,:)
      DO i = 1, 7
        IF (hetNumFrac(HetNreg-j+1,i) .LE. 1.e-60) hetRxFrac(HetNreg-j+1,i) = 0.
      END DO
    END DO
    
    DO i = 1, 7
      hetRxFrac(:,i) = hetRxFrac(:,i)/SUM(hetRxFrac(:,i))
      hetNumFrac(:,i) = hetNumFrac(:,i)/SUM(hetNumFrac(:,i))
    END DO
#endif    
    bupin = bupin/areasum ! Edit by LHG 19/03/14
    IF (areasum .LT. 1.e-40) bupin = 0.
    
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
        IF (nTracerCntl%lTranON) THEN
          KinParMac_velo(iso,ig)=KinParMac_velo(iso,ig)/PhiVol(ig)
          KinParMac_velo(iso,ig)=1._8/KinParMac_velo(iso,ig)
        END IF
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
