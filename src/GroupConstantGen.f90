#include <defines.h>
#define gcgen
#ifdef gcgen
SUBROUTINE GroupConstGen(Core, FmInfo, THInfo, CmInfo, GroupInfo, nTracerCntl, PE, ng)
USE PARAM
USE TYPEDEF,         ONLY : CoreInfo_Type,FmInfo_Type, THInfo_Type, CMInfo_Type,Cell_Type,Pin_Type, &
                                FxrInfo_Type, PinXs_Type, XsMac_Type, GroupInfo_Type,PE_TYPE, Powerdist_TYPE, &
                              Asy_Type, Asyinfo_type
USE BasicOperation, ONLY : CP_VA, CP_CA, MULTI_VA, MULTI_CA
USE files,           ONLY : caseid
USE XSLIB_MOD
USE XsUtil_mod, ONLY : FreeXsIsoMac
USE CNTL,         ONLY : nTracerCntl_Type
USE CritSpec_mod,     ONLY : GetDiffusionCoeff
USE GroupConst_Mod
USE Core_mod,         ONLY : eigv
USE TH_Mod,           ONLY : GetPinFuelTemp
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(LIBDATA), POINTER :: isodata
REAL :: TempRef   !THinfo%RefFuelTemp(:)
TYPE(XsMac_Type) :: XsMac

TYPE(FxrInfo_Type), POINTER :: myFXR
TYPE(CMInfo_Type) :: CMInfo
TYPE(THInfo_Type) :: THInfo
TYPE(FMInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
INTEGER, INTENT(IN) :: ng

REAL :: cellbphi(4,core%nxy,ng), cellbJ(4,core%nxy,ng), cellavgphi(core%nxy,ng)
REAL :: cellbphi2g(4,core%nxy,2), cellbJ2g(4,core%nxy,2), cellavgphi2g(core%nxy,2)
LOGICAL :: noString

INTEGER :: isoidx, iso, idiso, id, nisoinFxr, nisotot , isoName(isosize)
!INTEGER : isoList(500,2), isoName(500) ! 500 > 100 14/01/29 > 500 14/02/27
!LOGICAL :: checkiso(500)! 500 > 100 14/01/29  > 500 14/02/27
INTEGER, POINTER :: isoList(:,:)
LOGICAL, POINTER :: checkiso(:)
!REAL :: IsoMacXsSm(200,ng,ng) !(10,ng,ng)
REAL, POINTER :: IsoMacXsSm(:,:,:)
REAL, POINTER :: isoMacXsSmP1(:,:,:)  !for inflow correction 14/06/01
REAL :: XSt(ng) !, Sm1(ng,ng)

REAL,POINTER :: isoMicXs(:,:,:),isoMicSm(:,:,:), isoMicXs2g(:,:,:), isoMicSm2g(:,:,:) !tr, a, r, f, nu, k
REAL,POINTER :: isoMacXs(:,:,:),isoMacSm(:,:,:), isoMacXs2g(:,:,:), isoMacSm2g(:,:,:) !D, a, r, f, nf, kf
REAL :: isoNumden(isosize), FsrVolsum, RFsrVolsum, h2oNumden, uo2Numden ! O+H, O-H number density  ! 500 > 100  14/02/27
REAL :: FSV, isoFSV
REAL,POINTER :: asyFSV(:)

!REAL :: FsrPhi(Core%nCoreFsr,ng), FsrPhiVol(Core%nCoreFsr,ng), FsrVol(Core%nCoreFsr), PhiVol(ng), Vol
!REAL :: FsrPhi(1,ng), FsrPhiVol(1,ng), FsrVol(1), PhiVol(ng), Vol
REAL :: PhiVol(ng), Vol
REAL, POINTER :: FsrPhi(:,:), FsrPhiVol(:,:), FsrVol(:)
REAL :: myVol, myPhiVol, myNumden, hz
INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nLocalFsr, localFsr, iCel, iFxr, ifsr, nFsrInFxr, nCoreFsr, nCoreFxr

TYPE(PinXS_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: PhiC(:, :, :)
REAL :: Phi(ng,Core%nxy), Phi2g(2,Core%nxy)
REAL :: Bsq, Msq, kinf, keff
REAL :: f_ratio(ng)
REAL :: phicrit(ng), Dng(ng), XstrB1(ng), Xstr(ng)
REAL :: critphisum, infphisum, Rcritphisum, Rinfphisum

INTEGER :: npins, pidx(2,100), gidx, gidx2
REAL :: ADF(4,ng), ADFden,ADF2g(4,2), phiavg(2), Rnxy
REAL :: avgADF2g(2), avgbphi2g(2)
REAL :: conphi2g(3,2) !corner, center, avg/ ig
LOGICAL :: lhomADF = .TRUE.

INTEGER :: nxs = 6 !0 for total
INTEGER :: i, j, k
INTEGER :: ig, ig2, iz, ixy, igb, ige
INTEGER :: myzb, myze, nxy, nx, ny
INTEGER :: io, io2, io3

REAL :: phisum(ng), volsum
CHARACTER(256) :: fn,fn2,fn3

REAL :: fphi,tphi, rphi, p_ratio, mod_ratio

REAL :: EStruct(ng+1),Eavg(ng)
REAL :: ELB, temp
INTEGER :: EIdx

REAL, POINTER :: buexp(:), pwr(:)
REAL :: areasum
TYPE(PowerDist_Type) :: PowerDist
!--- definitions for TA(two-asy) problem---
LOGICAL :: lSA, lCB !lSingle Assembly, lCheckerBoard
INTEGER :: nxya, iasy, neiso
REAL :: AsyPhiVol(Core%nxya,ng), AsyPhiVolt(Core%nxya,ng)
REAL :: Asyphi2g(Core%nxya,2), Asyphi2gt(Core%nxya,2)
REAL :: surflux2g(2), surfJ2g(2), surflux2gt(2), surfJ2gt(2)
REAL,POINTER :: asyisoMicXs(:,:,:,:),asyisoMicSm(:,:,:,:), asyisoMicXs2g(:,:,:,:), asyisoMicSm2g(:,:,:,:) !tr, a, r, f, nu, k
REAL,POINTER :: asyisoMacXs(:,:,:,:),asyisoMacSm(:,:,:,:), asyisoMacXs2g(:,:,:,:), asyisoMacSm2g(:,:,:,:) !D, a, r, f, nf, kf
REAL,POINTER :: asyisoNumden(:,:), asyh2onumden(:), asyuo2numden(:)
REAL :: asyfphi, asytphi
REAL,POINTER :: asyDng(:,:), asyXstr(:,:)
REAL :: w1, w2, f1, f2, u, sigr, D1, D2

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
!----------INITIALIZE-----------
!14/02/27
CONTINUE
ALLOCATE(FsrPhi(Core%nCoreFsr,ng), FsrPhiVol(Core%nCoreFsr,ng), FsrVol(Core%nCoreFsr))
ALLOCATE(isoList(isosize,2), checkiso(isosize))
ALLOCATE(IsoMacXsSm(isosize,ng,ng))
ALLOCATE(IsoMacXsSmP1(isosize,ng,ng)) ! S1 scattering Matrix 14/06/01
ALLOCATE(cellphibdry(4,Core%nxy,3,ng)) !---bdry direction/cell#/1:numerator 2:denominator/group idx.

checkiso=.FALSE. ! 0=not Used. 1=used
nisotot=0; isoNumden=0; FsrVolsum=0; h2oNumden=0; uo2Numden=0;
XSt=0.0; !Sm1=0.0;
iz=1;
PinXS => CMInfo%PinXS; PHIC => CmInfo%PhiC
AsyInfo => Core%AsyInfo; Asy => Core%Asy
Pin => Core%Pin; Cell => Core%CellInfo
nxy = Core%nxy; nx = Core%nx; ny = Core%ny;
nCoreFsr=Core%nCoreFsr; nCoreFxr = 0
PhiVol=0; Phi=0; Vol=0; Asyphivol=0;
hz=Core%HZ(iz)
nxya=Core%nxya ! n-assemblies
lSA=.TRUE.
lCB=.FALSE.
IF( Core%nxya .NE. 1 )THEN
    lSA=.FALSE.
    IF( Core%nxya .EQ. 4 )THEN
        lCB=.TRUE.
        nxya=4
    ENDIF
    ALLOCATE(asyisonumden(nxya,isosize),asyh2onumden(nxya),asyuo2numden(nxya))
    asyisonumden=0; asyh2onumden=0; asyuo2numden=0;
ENDIF
!---------------------POWER and Burn-up Exposure calculation -----------
ALLOCATE(buexp(nxy),pwr(nxy))
iz = 1
DO ixy = 1, nxy
    areasum=0
    buexp(ixy)=0
    fxridxst=Pin(ixy)%fxridxst
    nLocalFxr=Pin(ixy)%nFXRmax
    DO ifxr = fxridxst, fxridxst+nlocalfxr-1
        IF( Fminfo%FXR(ifxr,iz)%lfuel )THEN
            areasum=areasum+Fminfo%FXR(ifxr,iz)%area
            buexp(ixy)=buexp(ixy)+Fminfo%FXR(ifxr,iz)%burnup*Fminfo%FXR(ifxr,iz)%area
        ENDIF
    ENDDO
    IF( areasum .NE. 0 )THEN
        buexp(ixy)=buexp(ixy)/areasum
    ELSE
        buexp(ixy)=0
    ENDIF
ENDDO

pwr=0
IF( lSA )THEN
CALL CorePowerCal(Core, CmInfo, PowerDist, ng, nTracerCntl, PE)
DO ixy = 1, nxy
    pwr(ixy)=powerdist%pinpower2d(ixy,iz)
ENDDO
ENDIF
CONTINUE



!--- Check Isotope LIST
DO ixy = 1, nxy
    icel = Pin(ixy)%Cell(iz)
    iasy = Pin(ixy)%iasy
    FxrIdxSt = Pin(ixy)%FxrIdxSt
    nlocalFxr = Cell(icel)%nFxr
    localFsr = 0
    DO j = 1, nLocalFxr
        nCoreFxr=nCoreFxr+1
        ifxr = FxrIdxSt + j -1
        myFxr => FmInfo%FXR(ifxr,iz)
        nFsrInFxr = myFxr%nFsrInFxr
        FsrIdxSt = myFxr%FsrIdxSt
        nisoInFxr = myFxr%niso
        DO k = 1, nFsrInFxr
            ifsr = FsrIdxSt + k -1
            localFsr=localFsr+1
            FsrVol(ifsr)=Cell(icel)%vol(localFsr)*hz
            myVol=FsrVol(ifsr)
            Vol=Vol+myVol
            DO ig = 1, ng
                FsrPhi(ifsr,ig)=FmInfo%phis(ifsr,iz,ig)
                FsrPhiVol(ifsr,ig)=FmInfo%phis(ifsr,iz,ig)*myVol
                PhiVol(ig)=PhiVol(ig)+FsrPhiVol(ifsr,ig)
            ENDDO
            FsrVolsum=FsrVolsum+myVol
            DO iso = 1, nisoInFxr
                idiso=myFxr%idiso(iso)
                IF( idiso .EQ. 8001 )THEN
                    idiso=8016
                ENDIF
                id = MapNucl(idiso)  !id = Map(92235) = 37
                isodata => ldiso(id)
                IF(.NOT. checkiso(id))THEN
                    nisotot=nisotot+1
                    checkiso(id) = .TRUE.
                    isoList(id,1)=nisotot
                    isoList(id,2)=idiso
                    isoidx=nisotot
                    isoName(isoidx)=idiso
                    !--- one more space for H2O-O and UO2-O
                    IF( idiso .EQ. 8016 )THEN
                        nisotot=nisotot+1
                        isoName(nisotot)=idiso
                    ENDIF
                ELSE
                    isoidx=isoList(id,1)
                ENDIF
                IF( idiso .EQ. 8016 )THEN
                    !--- originally O+H
                    IF( myFxr%lFuel )THEN !--- O-H
                        isoidx=isoidx+1
                    ENDIF
                ENDIF
                continue
                isoNumden(isoidx)=isoNumden(isoidx)+myVol*myFxr%pnum(iso)
                IF( idiso .EQ. 8016 )THEN
                    IF( myFxr%lH2O )THEN
                        h2oNumden=h2oNumden+myVol*myFxr%pnum(iso)
                    ELSEIF( myFxr%lFuel) THEN
                        uo2Numden=uo2Numden+myVol*myFxr%pnum(iso)
                    ENDIF
                !ELSEIF( myFxr%idiso(iso) .EQ. 8001 )THEN
                !    uo2Numden=uo2Numden+myVol*myFxr%pnum(iso)
                ENDIF
                IF( .NOT. lSA )THEN
                    asyisoNumden(iasy,isoidx)=asyisoNumden(iasy,isoidx)+myVol*myFxr%pnum(iso)
                    IF( idiso .EQ. 8016 )THEN
                        IF( myFxr%lH2O )THEN
                            asyh2oNumden(iasy)=asyh2oNumden(iasy)+myVol*myFxr%pnum(iso)
                        ELSEIF( myFxr%lFuel) THEN
                            asyuo2Numden(iasy)=asyuo2Numden(iasy)+myVol*myFxr%pnum(iso)
                        ENDIF
                    ENDIF
                ENDIF
            ENDDO
        ENDDO
    ENDDO
ENDDO
IF( lSA )THEN
    ALLOCATE(isoMicXs(nisotot,0:nxs,ng), isoMicSm(nisotot,ng,ng),isoMicXs2g(nisotot,0:nxs,2),isoMicSm2g(nisotot,2,2)) !tr, a, r, f, nu, k
    ALLOCATE(isoMacXs(0:nisotot,0:nxs,ng), isoMacSm(0:nisotot,ng,ng),isoMacXs2g(0:nisotot,0:nxs,2),isoMacSm2g(nisotot,2,2)) !D, a, r, f, nf, kf
    isoMacXs=0;isoMacSm=0;isoMacXs2g=0;isoMicXs=0;isoMicSm=0;isoMicXs2g=0; isoMicSm2g=0; isoMacSm2g=0;
ELSE
    ALLOCATE(isoMicXs(nisotot,0:nxs,ng), isoMicSm(nisotot,ng,ng),isoMicXs2g(nisotot,0:nxs,2),isoMicSm2g(nisotot,2,2)) !tr, a, r, f, nu, k
    ALLOCATE(isoMacXs(0:nisotot,0:nxs,ng), isoMacSm(0:nisotot,ng,ng),isoMacXs2g(0:nisotot,0:nxs,2),isoMacSm2g(nisotot,2,2)) !D, a, r, f, nf, kf
    isoMacXs=0;isoMacSm=0;isoMacXs2g=0;isoMicXs=0;isoMicSm=0;isoMicXs2g=0;isoMicSm2g=0; isoMacSm2g=0;
    ALLOCATE(asyisoMicXs(nxya,nisotot,0:nxs,ng), asyisoMicSm(nxya,nisotot,ng,ng),asyisoMicXs2g(nxya,nisotot,0:nxs,2),asyisoMicSm2g(nxya,nisotot,2,2)) !tr, a, r, f, nu, k
    ALLOCATE(asyisoMacXs(nxya,0:nisotot,0:nxs,ng), asyisoMacSm(nxya,0:nisotot,ng,ng),asyisoMacXs2g(nxya,0:nisotot,0:nxs,2),asyisoMacSm2g(nxya,nisotot,2,2)) !D, a, r, f, nf, kf
    asyisoMicXs=0;asyisoMicSm=0;asyisoMicXs2g=0; asyisoMicSm2g=0;
    asyisoMacXs=0;asyisoMacSm=0;asyisoMacXs2g=0; asyisoMacSm2g=0;
    ALLOCATE(asyDng(nxya,ng),asyXstr(nxya,ng))
    ALLOCATE(asyFSV(nxya))
    asyFSV=0;
ENDIF
IF( nisotot .GT. isosize )THEN
    write(*,*) '                                EXCEED. BYS edit 14/02/27'
ENDIF
RFsrVolSum=1.0/FsrVolsum
CALL MULTI_CA(RFsrVolsum, isoNumden(:), nisotot)!isoNumden=isoNumden/FsrVolsum ! isoNumden = avg(n)_iso
h2onumden=rFsrVolsum*h2onumden
uo2numden=rFsrVolsum*uo2numden
IF( .NOT. lSA )THEN
    DO iasy = 1, nxya
        CALL MULTI_CA(rFsrVolsum*nxya, asyisoNumden(iasy,:), nisotot)!isoNumden=isoNumden/FsrVolsum ! isoNumden = avg(n)_iso
        asyh2onumden(iasy)=rFsrVolsum*asyh2onumden(iasy)*nxya
        asyuo2numden(iasy)=rFsrVolsum*asyuo2numden(iasy)*nxya
    ENDDO
ENDIF
!------------------------------ GET flux spectrum -----------------------------------------------
SELECT CASE(nTracerCntl%gc_spec) !--- flux spectrum selection
CASE(0) !---critical spectrum
    CALL GetDiffusionCoeff(Core, FmInfo, GroupInfo, nTracerCntl, PE, Dng, phicrit, Bsq, Msq, kinf, keff)
    critphisum=0; infphisum=0;
    DO ig = 1, ng
        f_ratio(ig)=0
        DO ifsr = 1, nCoreFsr
            f_ratio(ig)=f_ratio(ig)+FsrPhiVol(ifsr,ig)
        ENDDO
        infphisum=infphisum+f_ratio(ig)
        critphisum=critphisum+phicrit(ig)
    ENDDO
    Rinfphisum=one/infphisum
    Rcritphisum=one/critphisum
    CALL MULTI_CA(Rinfphisum, f_ratio(:), ng)
    CALL MULTI_CA(Rcritphisum, phicrit(:), ng)
    DO ig = 1, ng
        f_ratio(ig)=phicrit(ig)/f_ratio(ig)
        PhiVol(ig)=PhiVol(ig)*f_ratio(ig) ! BYS_edit 13/11/01
        DO ifsr = 1, nCoreFsr
            FsrPhi(ifsr,ig)=FsrPhi(ifsr,ig)*f_ratio(ig)
            FsrPhivol(ifsr,ig)=FsrPhi(ifsr,ig)*FsrVol(ifsr)
        ENDDO
    ENDDO
CASE(1) !---Infinit medium spectrum _ inflow corection
    CALL GetDiffusionCoeff(Core, FmInfo, GroupInfo, nTracerCntl, PE, Dng, phicrit, Bsq, Msq, kinf, keff)
    kinf=Eigv
    DO ig = 1, ng
        f_ratio(ig)=1
    ENDDO
CASE(2) !---Infinit medium spectrum (no critical search calculation)
    Dng=0;    Bsq=0;  Msq=0
    kinf=Eigv
    DO ig = 1, ng
        f_ratio(ig)=1
    ENDDO

ENDSELECT

CONTINUE

!----------HOMOGENIZATION-------------------------------------------------------
FSV=0; !fission source *Volume
DO ixy = 1, nxy
    icel = Pin(ixy)%Cell(iz)
    iasy = Pin(ixy)%iasy
    FxrIdxSt = Pin(ixy)%FxrIdxSt
    nlocalFxr = Cell(icel)%nFxr
    Tempref = GetPinFuelTemp(Core, FmInfo%FXR, iz, ixy)
    DO j = 1, nLocalFxr
        ifxr = FxrIdxSt + j -1
        myFxr => FmInfo%FXR(ifxr,iz)
        nFsrInFxr = myFxr%nFsrInFxr
        FsrIdxSt = myFxr%FsrIdxSt
        nisoInFxr = myFxr%niso
        CALL GenFxrIsoMacXS(XsMac, IsoMacXsSm, IsoMacXsSmP1, FmInfo%FXR, ifxr, Tempref, 1, ng, Core, ixy, iz, GroupInfo, nTRACERCntl, PE)
        ! 13/11/01
        DO k = 1, nFsrInFxr
            ifsr = FsrIdxSt + k -1
            DO ig = 1, ng
                Phi(ig,ixy)=Phi(ig,ixy)+FsrPhiVol(ifsr,ig)
            ENDDO
            IF( .NOT. lSA )THEN
                DO ig = 1, ng
                    AsyPhiVol(iasy,ig)=AsyPhiVol(iasy,ig)+FsrPhiVol(ifsr,ig)
                ENDDO
            ENDIF
        ENDDO
        ! 13/11/01 end
        DO iso = 1, nisoInFxr
            idiso=myFxr%idiso(iso)
            IF( idiso .EQ. 8001 )THEN
                idiso=8016
            ENDIF
            id = MapNucl(idiso)  !id = Map(92235) = 37
            isodata => ldiso(id)
            isoidx=isoList(id,1)
            IF( idiso .EQ. 8016 )THEN
                IF( myFxr%lFuel )THEN
                    isoidx=isoidx+1
                ENDIF
            ENDIF
            myNumden=myFxr%pnum(iso)
            DO k = 1, nFsrInFxr
                ifsr = FsrIdxSt + k -1
                isoFSV=0
                DO ig2 = 1, ng
                    isoFSV=isoFSV+XsMac%isoXsMacKF(iso,ig2) *FsrPhiVol(ifsr,ig2)
                ENDDO
                FSV=FSV+isoFSV
                IF( .NOT. lSA )THEN !--- two assembly calculation
                    asyFSV(iasy)= asyFSV(iasy)+isoFSV
                ENDIF
                DO ig = 1, ng
                    myPhiVol=FsrPhiVol(ifsr,ig)
                    !Phi(ig,ixy)=Phi(ig,ixy)+myPhiVol
                    !MAC : D, a, r, f, nf, kf
                    !MIC : tr, a, r, f, nu, k
                    isoMacXs(isoidx,0,ig)=isoMacXs(isoidx,0,ig)+XsMac%isoXsMacT(iso,ig) *myPhiVol  !Total
                    isoMacXs(0,0,ig)=isoMacXs(0,0,ig)+XsMac%isoXsMacT(iso,ig) *myPhiVol
                    isoMacXs(isoidx,1,ig)=isoMacXs(isoidx,1,ig)+XsMac%isoXsMacTR(iso,ig) *myPhiVol  !D
                    isoMacXs(0,1,ig)=isoMacXs(0,1,ig)+XsMac%isoXsMacTR(iso,ig) *myPhiVol
                    isoMacXs(isoidx,2,ig)=isoMacXs(isoidx,2,ig)+XsMac%isoXsMacA(iso,ig)  *myPhiVol  !absorption
                    isoMacXs(0,2,ig)=isoMacXs(0,2,ig)+XsMac%isoXsMacA(iso,ig) *myPhiVol
                    !--- transport removal
                    !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+(IsoMacXsSm(iso,ig,ig))*myPhiVol  !ss
                    !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+(XsMac%isoXsMacTR(iso,ig)-IsoMacXsSm(iso,ig,ig))*myPhiVol  !removal
                    !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+(XsMac%isoXsMacTR(iso,ig))*myPhiVol  !tr
                    !--- total removal
                    !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+XsMac%isoXsMacA(iso,ig)  *myPhiVol  ! abs
                    isoMacXs(0,3,ig)=isoMacXs(0,3,ig)+XsMac%isoXsMacA(iso,ig) *myPhiVol
                    DO ig2 = 1, ng
                        IF( ig2 .NE. ig )THEN
                             !isoMacXs(isoidx,3,ig)=isoMacXs(isoidx,3,ig)+IsoMacXsSm(iso,ig,ig2) *myPhiVol  !removal
                             isoMacXs(0,3,ig)=isoMacXs(0,3,ig)+IsoMacXsSm(iso,ig,ig2) *myPhiVol  !removal
                        ENDIF
                    ENDDO

                    isoMacXs(isoidx,4,ig)=isoMacXs(isoidx,4,ig)+XsMac%isoXsMacF(iso,ig)  *myPhiVol  !fission
                    isoMacXs(0,4,ig)=isoMacXs(0,4,ig)+XsMac%isoXsMacF(iso,ig) *myPhiVol
                    !--- no fission -> no nu or not?
                    !IF( isoMacXs(isoidx,4,ig) .NE. 0 )THEN
                        isoMacXs(isoidx,5,ig)=isoMacXs(isoidx,5,ig)+XsMac%isoXsMacNF(iso,ig) *myPhiVol  !nu-fission
                        isoMacXs(0,5,ig)=isoMacXs(0,5,ig)+XsMac%isoXsMacNF(iso,ig) *myPhiVol
                        isoMacXs(isoidx,6,ig)=isoMacXs(isoidx,6,ig)+XsMac%isoXsMacKF(iso,ig) *myPhiVol  !kappa-fission
                        IF( ig .LE. groupinfo%nchi .AND. isoFSV .NE. 0)THEN
                        isoMacXs(0,6,ig)=isoMacXs(0,6,ig)+myFXR%CHI(ig) *isoFSV   ! CHI
                        ENDIF
                    !ENDIF
                    DO ig2 = 1, ng
                        isoMacSm(isoidx,ig,ig2)=isoMacSm(isoidx,ig,ig2)+IsoMacXsSm(iso,ig,ig2)*myPhiVol
                        isoMacSm(0,ig,ig2)=isoMacSm(0,ig,ig2)+IsoMacXsSm(iso,ig,ig2)*myPhiVol
                        !Sm1(ig,ig2)=Sm1(ig,ig2)+IsoMacXsSm1(iso,ig,ig2)*myPhiVol
                    ENDDO
                ENDDO !---END of G sweep
                IF( .NOT. lSA )THEN !--- two assembly calculation
                    DO ig = 1, ng
                        myPhiVol=FsrPhiVol(ifsr,ig)
                        !MAC : D, a, r, f, nf, kf
                        !MIC : tr, a, r, f, nu, k
                        asyisoMacXs(iasy,0,0,ig)=asyisoMacXs(iasy,0,0,ig)+XsMac%isoXsMacT(iso,ig) *myPhiVol  !total
                        asyisoMacXs(iasy,isoidx,1,ig)=asyisoMacXs(iasy,isoidx,1,ig)+XsMac%isoXsMacTR(iso,ig) *myPhiVol  !D
                        asyisoMacXs(iasy,0,1,ig)=asyisoMacXs(iasy,0,1,ig)+XsMac%isoXsMacTR(iso,ig) *myPhiVol  !D
                        asyisoMacXs(iasy,isoidx,2,ig)=asyisoMacXs(iasy,isoidx,2,ig)+XsMac%isoXsMacA(iso,ig)  *myPhiVol  !absorption
                        asyisoMacXs(iasy,0,2,ig)=asyisoMacXs(iasy,0,2,ig)+XsMac%isoXsMacA(iso,ig)  *myPhiVol  !absorption
                        asyisoMacXs(iasy,0,3,ig)=asyisoMacXs(iasy,0,3,ig)+XsMac%isoXsMacA(iso,ig)  *myPhiVol
                        DO ig2 = 1, ng
                            IF( ig2 .NE. ig )THEN
                                asyisoMacXs(iasy,0,3,ig)=asyisoMacXs(iasy,0,3,ig)+IsoMacXsSm(iso,ig,ig2) *myPhiVol  !removal
                            ENDIF
                        ENDDO
                        asyisoMacXs(iasy,isoidx,4,ig)=asyisoMacXs(iasy,isoidx,4,ig)+XsMac%isoXsMacF(iso,ig)  *myPhiVol  !fission
                        asyisoMacXs(iasy,0,4,ig)=asyisoMacXs(iasy,0,4,ig)+XsMac%isoXsMacF(iso,ig)  *myPhiVol  !fission
                        asyisoMacXs(iasy,isoidx,5,ig)=asyisoMacXs(iasy,isoidx,5,ig)+XsMac%isoXsMacNF(iso,ig) *myPhiVol  !nu-fission
                        asyisoMacXs(iasy,0,5,ig)=asyisoMacXs(iasy,0,5,ig)+XsMac%isoXsMacNF(iso,ig) *myPhiVol  !nu-fission !14/05/01
                        asyisoMacXs(iasy,isoidx,6,ig)=asyisoMacXs(iasy,isoidx,6,ig)+XsMac%isoXsMacKF(iso,ig) *myPhiVol  !kappa-fission
                        IF( ig .LE. groupinfo%nchi .AND. isoFSV .NE. 0)THEN
                        asyisoMacXs(iasy,0,6,ig)=asyisoMacXs(iasy,0,6,ig)+myFXR%CHI(ig) *isoFSV   ! CHI
                        ENDIF
                        DO ig2 = 1, ng
                            asyisoMacSm(iasy,isoidx,ig,ig2)=asyisoMacSm(iasy,isoidx,ig,ig2)+IsoMacXsSm(iso,ig,ig2)*myPhiVol
                            asyisoMacSm(iasy,0,ig,ig2)=asyisoMacSm(iasy,0,ig,ig2)+IsoMacXsSm(iso,ig,ig2)*myPhiVol
                        ENDDO
                    ENDDO
                ENDIF
            ENDDO !---END of Fsr sweep
        ENDDO !---END of Iso sweep
        !DEALLOCATE(IsoMacXsSm)
        CALL FreeXsIsoMac(XsMac)
    ENDDO !---END of Fxr sweep
ENDDO !---END of nPin sweep
DO ig = 1, ng
    DO ixy = 1, nxy
        Phi(ig,ixy)=Phi(ig,ixy)/Core%PinVol(ixy,iz)
    ENDDO
ENDDO

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
iso=0; !homogenized Macro XS 47G iso=0
DO ig = 1, ng
    DO i = 0, nxs-1
        isoMacXs(iso,i,ig)=isoMacXs(iso,i,ig)/PhiVol(ig)
    ENDDO
    XSt(ig)=isoMacXs(iso,0,ig)
    isoMacXs(iso,6,ig)=isoMacXs(iso,6,ig)/FSV !chi
    DO ig2 = 1, ng
        isoMacSm(iso,ig,ig2)=isoMacSm(iso,ig,ig2)/PhiVol(ig)
        !Sm1(ig,ig2)=Sm1(ig,ig2)/PhiVol(ig)
    ENDDO
ENDDO
IF( .NOT. lSA )THEN !--- two assembly calculation
    DO iasy = 1, nxya
        DO iso = 0, nisotot !0 for homogenized Macro XS iso=0 14/05/01
            DO ig = 1, ng
                DO i = 0, nxs
                    IF(iso .EQ. 0 .AND. i .EQ. 6)THEN
                        asyisoMacXs(iasy,iso,6,ig)=asyisoMacXs(iasy,iso,6,ig)/asyFSV(iasy) ! 14/09/16 corrected : kappa fission error
                    ELSE
                        asyisoMacXs(iasy,iso,i,ig)=asyisoMacXs(iasy,iso,i,ig)/AsyPhiVol(iasy,ig)
                    ENDIF
                ENDDO
                DO ig2 = 1, ng
                    asyisoMacSm(iasy,iso,ig,ig2)=asyisoMacSm(iasy,iso,ig,ig2)/AsyPhiVol(iasy,ig)
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDIF
!--- inflow corrected Diffusion Coefficient
!CALL GetDiffusionCoeff_IFC(Dng, phivol, XSt, Sm1, ng)
!--- B1 diffusion correction
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
    DO iso = 1, nisotot
        isoMacXs(iso,1,ig) = XstrB1(ig) * isoMacXs(iso,1,ig)/Xstr(ig)
    ENDDO
ENDDO
IF( .NOT. lSA )THEN
    DO iasy = 1, nxya
        asyXstr(iasy,:)=0
        DO ig = 1, ng
            DO iso = 1, nisotot
                asyXstr(iasy,ig)= asyXstr(iasy,ig)+asyisoMacXs(iasy,iso,1,ig)
            ENDDO
            asyDng(iasy,ig)=1.0/(3*asyXstr(iasy,ig))
        ENDDO
    ENDDO
ENDIF

DO ig = 1, ng
    DO iso = 1, nisotot
        !MAC : D, a, r, f, nf, kf
        !MIC : tr, a, r, f, nu, k
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
    ENDDO
ENDDO
IF( .NOT. lSA )THEN !--- two assembly calculation
    DO iasy = 1, nxya
    DO ig = 1, ng
        DO iso = 1, nisotot
            asyisoMicXs(iasy,iso,1,ig)=asyisoMacXs(iasy,iso,1,ig)/asyisoNumden(iasy,iso)      ! tr
            asyisoMicXs(iasy,iso,2,ig)=asyisoMacXs(iasy,iso,2,ig)/asyisoNumden(iasy,iso)      ! absorption
            asyisoMicXs(iasy,iso,4,ig)=asyisoMacXs(iasy,iso,4,ig)/asyisoNumden(iasy,iso)      ! fission
            IF( asyisoMacXs(iasy,iso,4,ig) .NE. 0 )THEN
                asyisoMicXs(iasy,iso,5,ig)=asyisoMacXs(iasy,iso,5,ig)/asyisoMacXs(iasy,iso,4,ig)  ! nu
                asyisoMicXs(iasy,iso,6,ig)=asyisoMacXs(iasy,iso,6,ig)/asyisoMacXs(iasy,iso,4,ig)  ! kappa
            ELSE
                asyisoMicXs(iasy,iso,5,ig)=0  ! nu
                asyisoMicXs(iasy,iso,6,ig)=0  ! kappa
            ENDIF
            DO ig2 = 1, ng
                asyisoMicSm(iasy,iso,ig,ig2)=asyisoMacSm(iasy,iso,ig,ig2)/asyisoNumden(iasy,iso)
            ENDDO
        ENDDO
    ENDDO
    ENDDO
ENDIF


!-------------- 2group condensing-----------------------------------------------
!--------Energy spectrum selection
IF (ng .eq. 47) THEN
    EStruct(1:ng) = enbhel(1:ng)
    EStruct(ng+1) = 1.0E-4_8
    DO ig = 1, ng
      Eavg(ig) = (EStruct(ig) + EStruct(ig+1))/2._8
    ENDDO
    ELB = 6.2506E-01  !---- thermal neutron criteria
    DO ig = 1, ng-1
        temp = (ELB-Eavg(ig))*(ELB-Eavg(ig+1))
        IF(temp .LE. 0._8) THEN
            EIdx=ig !---1~EIdx-1 : fast , EIdx~ng : thermal / EIdx=35 / 1~34, 35~47
            EXIT
        ENDIF
    ENDDO
    EIdx=EIdx+1
ELSEIF (ng .eq. 7) THEN
    EIdx=4
ELSEIF (ng .eq. 2) THEN
    EIdx=2
ENDIF


fphi=0; tphi=0;
DO ig = 1, ng
    IF (ig .LT. EIdx) THEN
        gidx=1
        fphi=fphi+PhiVol(ig)
    ELSE
        gidx=2
        tphi=tphi+PhiVol(ig)
    ENDIF
ENDDO
!D, a, r, f, nf, kf
!ALLOCATE(isoMacXs(nisotot,6,ng), isoMacSm(nisotot,ng,ng),isoMacXs2g(nisotot,6,2))
!tr, a, r, f, nu, k
!ALLOCATE(isoMicXs(nisotot,6,ng), isoMicSm(nisotot,ng,ng),isoMicXs2g(nisotot,6,2))

!--- Diffusion coefficient & Transport XS
!--- Diffusion coeff.
DO ig = 1, ng
    IF (ig .LT. EIdx) THEN
        gidx=1
    ELSE
        gidx=2
    ENDIF
    isoMacXs2g(0,1,gidx)=isoMacXs2g(0,1,gidx)+Dng(ig)*PhiVol(ig)
ENDDO
DO gidx = 1, 2
    IF( gidx .EQ. 1 )THEN
        rphi = 1.0/fphi
    ELSE
        rphi = 1.0/tphi
    ENDIF
    isoMacXs2g(0,1,gidx)=isoMacXs2g(0,1,gidx)*rphi
ENDDO

!--- Microscopic 2G XS
DO iso = 1, nisotot
    DO ig = 1, ng
        IF (ig .LT. EIdx) THEN
            gidx=1
        ELSE
            gidx=2
        ENDIF
        isoMacXs2g(iso,1,gidx)=isoMacXs2g(iso,1,gidx)+ 1/isoMacXs(iso,1,ig)*PhiVol(ig) ! 1/tr
        isoMacXs2g(iso,2,gidx)=isoMacXs2g(iso,2,gidx)+ isoMacXs(iso,2,ig)*PhiVol(ig) ! a
        isoMacXs2g(iso,3,gidx)=isoMacXs2g(iso,3,gidx)+ isoMacXs(iso,2,ig)*PhiVol(ig) ! r= !a! + s_out
        isoMacXs2g(iso,4,gidx)=isoMacXs2g(iso,4,gidx)+ isoMacXs(iso,4,ig)*PhiVol(ig) ! f
        isoMacXs2g(iso,5,gidx)=isoMacXs2g(iso,5,gidx)+ isoMacXs(iso,5,ig)*PhiVol(ig) ! nf
        isoMacXs2g(iso,6,gidx)=isoMacXs2g(iso,6,gidx)+ isoMacXs(iso,6,ig)*PhiVol(ig) ! kf
        DO ig2 = 1, ng
            IF (ig2 .LT. EIdx) THEN
                gidx2=1
            ELSE
                gidx2=2
            ENDIF
            isoMacSm2g(iso,gidx,gidx2)=isoMacSm2g(iso,gidx,gidx2)+isoMacSm(iso,ig,ig2)*PhiVol(ig)
            IF( gidx .NE. gidx2 )THEN
                isoMacXs2g(iso,3,gidx)=isoMacXs2g(iso,3,gidx)+ isoMacSm(iso,ig,ig2)*PhiVol(ig) ! r= !a! + s_out
            ENDIF
        ENDDO
    ENDDO
ENDDO
DO iso = 1, nisotot
    DO gidx = 1, 2
        IF( gidx .EQ. 1 )THEN
            rphi = 1.0/fphi
        ELSE
            rphi = 1.0/tphi
        ENDIF
        isoMacXs2g(iso,1,gidx)=isoMacXs2g(iso,1,gidx)*rphi ! 1/tr_iso_G
        isoMacXs2g(iso,1,gidx)=1/(isoMacXs2g(iso,1,gidx)) ! 1/tr -> tr_iso_G
        isoMicXs2g(iso,1,gidx)=isoMacXs2g(iso,1,gidx)/isoNumden(iso)  ! sig_tr

        isoMacXs2g(iso,2,gidx)=isoMacXs2g(iso,2,gidx)*rphi  ! a
        isoMicXs2g(iso,2,gidx)=isoMacXs2g(iso,2,gidx)/isoNumden(iso)
        isoMacXs2g(0,2,gidx) = isoMacXs2g(0,2,gidx) + isoMacXs2g(iso,2,gidx)

        isoMacXs2g(iso,3,gidx)=isoMacXs2g(iso,3,gidx)*rphi  ! r
        isoMicXs2g(iso,3,gidx)=isoMacXs2g(iso,3,gidx)/isoNumden(iso)
        isoMacXs2g(0,3,gidx) = isoMacXs2g(0,3,gidx) + isoMacXs2g(iso,3,gidx)

        IF( isoMacXs2g(iso,4,gidx) .NE. 0 )THEN
            isoMacXs2g(iso,4,gidx)=isoMacXs2g(iso,4,gidx)*rphi  ! f
            isoMicXs2g(iso,4,gidx)=isoMacXs2g(iso,4,gidx)/isoNumden(iso) !f
            isoMacXs2g(0,4,gidx) = isoMacXs2g(0,4,gidx) + isoMacXs2g(iso,4,gidx)

            isoMacXs2g(iso,5,gidx)=isoMacXs2g(iso,5,gidx)*rphi  ! isoMac nuf
            isoMicXs2g(iso,5,gidx)=isoMacXs2g(iso,5,gidx)/isoMacXs2g(iso,4,gidx) ! micro nu
            isoMacXs2g(0,5,gidx) = isoMacXs2g(0,5,gidx) + isoMacXs2g(iso,5,gidx) ! Macro nuf

            isoMacXs2g(iso,6,gidx)=isoMacXs2g(iso,6,gidx)*rphi  ! kf
            isoMicXs2g(iso,6,gidx)=isoMacXs2g(iso,6,gidx)/isoMacXs2g(iso,4,gidx) ! kappa
            isoMacXs2g(0,6,gidx) = isoMacXs2g(0,6,gidx) + isoMacXs2g(iso,6,gidx)
        ELSE
            isoMacXs2g(iso,4,gidx)=0 ! f
            isoMicXs2g(iso,4,gidx)=0 ! f
            isoMacXs2g(iso,5,gidx)=0 ! nuf
            isoMicXs2g(iso,5,gidx)=0 ! nu
            isoMacXs2g(iso,6,gidx)=0 ! kf
            isoMicXs2g(iso,6,gidx)=0 ! kappa
        ENDIF

    ENDDO
ENDDO
IF( .NOT. lSA )THEN
DO iasy = 1, nxya
    asyfphi=0; asytphi=0;
    DO ig = 1, ng
        IF (ig .LT. EIdx) THEN
            gidx=1
            asyfphi=asyfphi+asyPhiVol(iasy,ig)
        ELSE
            gidx=2
            asytphi=asytphi+asyPhiVol(iasy,ig)
        ENDIF
    ENDDO
    !--- Diffusion coefficient & Transport XS
    !--- Diffusion coeff.
    DO ig = 1, ng
        IF (ig .LT. EIdx) THEN
            gidx=1
        ELSE
            gidx=2
        ENDIF
        asyisoMacXs2g(iasy,0,1,gidx)=asyisoMacXs2g(iasy,0,1,gidx)+asyDng(iasy,ig)*asyPhiVol(iasy,ig)
    ENDDO
    DO gidx = 1, 2
        IF( gidx .EQ. 1 )THEN
            rphi = 1.0/asyfphi
        ELSE
            rphi = 1.0/asytphi
        ENDIF
        asyisoMacXs2g(iasy,0,1,gidx)=asyisoMacXs2g(iasy,0,1,gidx)*rphi
    ENDDO

    !--- Microscopic 2G XS
    DO iso = 1, nisotot
        DO ig = 1, ng
            IF (ig .LT. EIdx) THEN
                gidx=1
            ELSE
                gidx=2
            ENDIF
            IF( ig .EQ. 34 )THEN
                CONTINUE
            ENDIF

            asyisoMacXs2g(iasy,iso,1,gidx)=asyisoMacXs2g(iasy,iso,1,gidx)+ 1/asyisoMacXs(iasy,iso,1,ig)*asyPhiVol(iasy,ig) ! 1/tr
            asyisoMacXs2g(iasy,iso,2,gidx)=asyisoMacXs2g(iasy,iso,2,gidx)+ asyisoMacXs(iasy,iso,2,ig)*asyPhiVol(iasy,ig) ! a
            asyisoMacXs2g(iasy,iso,3,gidx)=asyisoMacXs2g(iasy,iso,3,gidx)+ asyisoMacXs(iasy,iso,2,ig)*asyPhiVol(iasy,ig) ! r= !a! + s_out
            asyisoMacXs2g(iasy,iso,4,gidx)=asyisoMacXs2g(iasy,iso,4,gidx)+ asyisoMacXs(iasy,iso,4,ig)*asyPhiVol(iasy,ig) ! f
            asyisoMacXs2g(iasy,iso,5,gidx)=asyisoMacXs2g(iasy,iso,5,gidx)+ asyisoMacXs(iasy,iso,5,ig)*asyPhiVol(iasy,ig) ! nf
            asyisoMacXs2g(iasy,iso,6,gidx)=asyisoMacXs2g(iasy,iso,6,gidx)+ asyisoMacXs(iasy,iso,6,ig)*asyPhiVol(iasy,ig) ! kf
            DO ig2 = 1, ng
                IF (ig2 .LT. EIdx) THEN
                    gidx2=1
                ELSE
                    gidx2=2
                ENDIF
                asyisoMacSm2g(iasy,iso,gidx,gidx2)=asyisoMacSm2g(iasy,iso,gidx,gidx2)+asyisoMacSm(iasy,iso,ig,ig2)*asyPhiVol(iasy,ig)
                IF( gidx .NE. gidx2 )THEN
                    asyisoMacXs2g(iasy,iso,3,gidx)=asyisoMacXs2g(iasy,iso,3,gidx)+ asyisoMacSm(iasy,iso,ig,ig2)*asyPhiVol(iasy,ig) ! r= !a! + s_out
                ENDIF
            ENDDO
        ENDDO
    ENDDO
    DO iso = 1, nisotot
        DO gidx = 1, 2
            IF( gidx .EQ. 1 )THEN
                rphi = 1.0/asyfphi
            ELSE
                rphi = 1.0/asytphi
            ENDIF
            asyisoMacXs2g(iasy,iso,1,gidx)=asyisoMacXs2g(iasy,iso,1,gidx)*rphi ! 1/tr_iso_G
            asyisoMacXs2g(iasy,iso,1,gidx)=1/(asyisoMacXs2g(iasy,iso,1,gidx)) ! 1/tr -> tr_iso_G
            asyisoMicXs2g(iasy,iso,1,gidx)=asyisoMacXs2g(iasy,iso,1,gidx)/asyisoNumden(iasy,iso)  ! sig_tr

            asyisoMacXs2g(iasy,iso,2,gidx)=asyisoMacXs2g(iasy,iso,2,gidx)*rphi  ! a
            asyisoMicXs2g(iasy,iso,2,gidx)=asyisoMacXs2g(iasy,iso,2,gidx)/asyisoNumden(iasy,iso)
            asyisoMacXs2g(iasy,0,2,gidx) = asyisoMacXs2g(iasy,0,2,gidx) + asyisoMacXs2g(iasy,iso,2,gidx)

            asyisoMacXs2g(iasy,iso,3,gidx)=asyisoMacXs2g(iasy,iso,3,gidx)*rphi  ! r
            asyisoMicXs2g(iasy,iso,3,gidx)=asyisoMacXs2g(iasy,iso,3,gidx)/asyisoNumden(iasy,iso)
            asyisoMacXs2g(iasy,0,3,gidx) = asyisoMacXs2g(iasy,0,3,gidx) + asyisoMacXs2g(iasy,iso,3,gidx)

            IF( asyisoMacXs2g(iasy,iso,4,gidx) .NE. 0 )THEN
                asyisoMacXs2g(iasy,iso,4,gidx)=asyisoMacXs2g(iasy,iso,4,gidx)*rphi  ! f
                asyisoMicXs2g(iasy,iso,4,gidx)=asyisoMacXs2g(iasy,iso,4,gidx)/asyisoNumden(iasy,iso) !f
                asyisoMacXs2g(iasy,0,4,gidx) = asyisoMacXs2g(iasy,0,4,gidx) + asyisoMacXs2g(iasy,iso,4,gidx)

                asyisoMacXs2g(iasy,iso,5,gidx)=asyisoMacXs2g(iasy,iso,5,gidx)*rphi  ! isoMac nuf
                asyisoMicXs2g(iasy,iso,5,gidx)=asyisoMacXs2g(iasy,iso,5,gidx)/asyisoMacXs2g(iasy,iso,4,gidx) ! micro nu
                asyisoMacXs2g(iasy,0,5,gidx) = asyisoMacXs2g(iasy,0,5,gidx) + asyisoMacXs2g(iasy,iso,5,gidx) ! Macro nuf

                asyisoMacXs2g(iasy,iso,6,gidx)=asyisoMacXs2g(iasy,iso,6,gidx)*rphi  ! kf
                asyisoMicXs2g(iasy,iso,6,gidx)=asyisoMacXs2g(iasy,iso,6,gidx)/asyisoMacXs2g(iasy,iso,4,gidx) ! kappa
                asyisoMacXs2g(iasy,0,6,gidx) = asyisoMacXs2g(iasy,0,6,gidx) + asyisoMacXs2g(iasy,iso,6,gidx) ! Macro kf
            ELSE
                asyisoMacXs2g(iasy,iso,4,gidx)=0 ! f
                asyisoMicXs2g(iasy,iso,4,gidx)=0 ! f
                asyisoMacXs2g(iasy,iso,5,gidx)=0 ! nuf
                asyisoMicXs2g(iasy,iso,5,gidx)=0 ! nu
                asyisoMacXs2g(iasy,iso,6,gidx)=0 ! kf
                asyisoMicXs2g(iasy,iso,6,gidx)=0 ! kappa
            ENDIF
        ENDDO
    ENDDO
ENDDO
ENDIF



continue

!----------Assembly Discontinuity Factor-------------------------
!-critical spectrum correction
DO ig = 1, ng
    DO ixy = 1, nxy
        DO i = 1, 4
            DO j = 1, 3
                CellPhiBdry(i,ixy,j,ig)=CellPhiBdry(i,ixy,j,ig)*f_ratio(ig)
            ENDDO
        ENDDO
    ENDDO
ENDDO !correction ends

!---ADF for single assembly-----------------------------------------------------
IF( lSA )THEN
DO ig = 1, ng
    ADFden=0
    DO ixy = 1, nxy
        ADFden = ADFden + cellphibdry(1,ixy,2,ig)
        Cellavgphi(ixy,ig) = cellphibdry(1,ixy,2,ig)
        DO i = 1, 4
            CellbPhi(i,ixy,ig) = cellphibdry(i,ixy,1,ig)
            CellbJ(i,ixy,ig) = cellphibdry(i,ixy,3,ig)
        ENDDO
    ENDDO
    ADFden=ADFden/nxy
    DO i = 1, 4
        ADF(i, ig)=0
        SELECT CASE(i)
        CASE(1) !---SOUTH
            npins=nx
            DO j = 1, npins
                pidx(1,j)=nxy-nx+j
                pidx(2,j)=pidx(1,j)-nx
            ENDDO
        CASE(2) !---WEST
            npins=ny
            DO j = 1, npins
                pidx(1,j)=(j-1)*nx+1
                pidx(2,j)=pidx(1,j)+1
            ENDDO
        CASE(3) !---NORTH
            npins=nx
            DO j = 1, npins
                pidx(1,j)=j
                pidx(2,j)=pidx(1,j)+nx
            ENDDO
        CASE(4) !---EAST
            npins=ny
            DO j = 1, npins
                pidx(1,j)=j*nx
                pidx(2,j)=pidx(1,j)-1
            ENDDO
        ENDSELECT
        IF(.FALSE.)THEN ! T:het ADF / F:hom ADF
            DO j = 1, npins
                ADF(i,ig)=ADF(i,ig) + cellphibdry(i,j,1,ig) !surf/pin/1/g
            ENDDO
            ADF(i,ig)=ADF(i,ig)/(ADFden*npins)
        ELSE
            DO j = 1, npins
                ADF(i, ig) = ADF(i, ig) + ( 3*Phi(ig,pidx(1,j)) - Phi(ig,pidx(2,j)) )*0.5
            ENDDO
            ADF(i, ig)=ADF(i, ig)*volsum/(npins*phisum(ig))
        ENDIF
    ENDDO
ENDDO
!----------Assembly Discontinuity Factor in 2-Groups
phiavg=0
phi2g=0
Cellavgphi2g=0
CellbPhi2g=0
CellbJ2g=0
DO ig = 1, ng
    IF (ig .LT. EIdx) THEN
        gidx=1
    ELSE
        gidx=2
    ENDIF
    DO ixy = 1, nxy
        IF( .NOT.lhomADF )THEN
            phiavg(gidx) = phiavg(gidx) + cellphibdry(1,ixy,2,ig)
        ELSE
            phiavg(gidx)=phiavg(gidx) + Phi(ig,ixy)
        ENDIF
        phi2g(gidx,ixy) = phi2g(gidx,ixy) + Phi(ig,ixy)
        !--- ADF verfication work
        Cellavgphi2g(ixy,gidx) = Cellavgphi2g(ixy,gidx) + Cellavgphi(ixy,ig)
        DO i = 1, 4
            CellbPhi2g(i,ixy,gidx) = CellbPhi2g(i,ixy,gidx) + CellbPhi(i,ixy,ig)
            CellbJ2g(i,ixy,gidx) = CellbJ2g(i,ixy,gidx) + CellbJ(i,ixy,ig)
        ENDDO
    ENDDO
ENDDO

Rnxy=one/nxy
CALL MULTI_CA(Rnxy, phiavg(1:2), 2)
avgADF2g=0
DO i = 1, 4
    SELECT CASE(i)
    CASE(1) !---SOUTH
        npins=nx
        DO j = 1, npins
            pidx(1,j)=nxy-nx+j
            pidx(2,j)=pidx(1,j)-nx
        ENDDO
    CASE(2) !---WEST
        npins=ny
        DO j = 1, npins
            pidx(1,j)=(j-1)*nx+1
            pidx(2,j)=pidx(1,j)+1
        ENDDO
    CASE(3) !---NORTH
        npins=nx
        DO j = 1, npins
            pidx(1,j)=j
            pidx(2,j)=pidx(1,j)+nx
        ENDDO
    CASE(4) !---EAST
        npins=ny
        DO j = 1, npins
            pidx(1,j)=j*nx
            pidx(2,j)=pidx(1,j)-1
        ENDDO
    ENDSELECT
    DO ig = 1, ng
        IF (ig .LT. EIdx) THEN
            gidx=1
        ELSE
            gidx=2
        ENDIF
        DO j = 1, npins
            IF( .NOT.lhomADF )THEN
                ADF2g(i, gidx) = ADF2g(i, gidx) + cellphibdry(i,pidx(1,j),1,ig)
            ELSE
                ADF2g(i, gidx) = ADF2g(i, gidx) + ( 3*Phi(ig,pidx(1,j)) - Phi(ig,pidx(2,j)) )*0.5
            ENDIF
        ENDDO
    ENDDO
    DO ig = 1, 2
        avgbphi2g(ig)=avgbphi2g(ig)+ADF2g(i, ig)/(4*npins)
        ADF2g(i, ig)=ADF2g(i, ig)/(npins*phiavg(ig))
        avgADF2g(ig)=avgADF2g(ig)+ADF2g(i,ig)*0.25
    ENDDO
ENDDO
!Corner flux
conphi2g=0
DO ig = 1, ng
    IF (ig .LT. EIdx) THEN
        gidx=1
    ELSE
        gidx=2
    ENDIF
    IF( nx .GT. 1 )THEN ! expection for single pin problem
        conphi2g(1,gidx)=conphi2g(1,gidx)+cellphibdry(3,1,1,ig)
        conphi2g(2,gidx)=conphi2g(2,gidx)+cellphibdry(3,nx/2,1,ig)
    ELSE
        conphi2g(1,gidx)=0
        conphi2g(2,gidx)=0
    ENDIF
ENDDO
DO gidx=1,2
    conphi2g(3,gidx)=(conphi2g(1,gidx)+conphi2g(2,gidx))/2
ENDDO

CONTINUE
mod_ratio=(fphi+tphi)/(phiavg(1)+phiavg(2))
p_ratio=nTracerCntl%PowerCore/(isoMacXs2g(0,6,1)*fphi+isoMacXs2g(0,6,2)*tphi)!*1.0E6
fphi=fphi*p_ratio; tphi=tphi*p_ratio
p_ratio=p_ratio*mod_ratio
DO ig = 1, 2
    !rphi=1
    avgbPhi2g(ig)=avgbPhi2g(ig)*p_ratio
    DO i = 1, 3
        conphi2g(i,ig)=conphi2g(i,ig)*p_ratio
    ENDDO
ENDDO

ELSE !-TWO-ASSEMBLIES-------------------------------------------------------------------------------------
    Asyphivolt=0;
    Asyphi2g=0; Asyphi2gt=0;
    surflux2g=0; surfJ2g=0;
    surflux2gt=0; surfJ2gt=0;
    DO ig = 1, ng
        IF (ig .LT. EIdx) THEN
            gidx=1
        ELSE
            gidx=2
        ENDIF
        DO iasy = 1, nxya
            !Asyphi2g(iasy,gidx)=AsyPhiVol(iasy,ig)
        ENDDO
        DO ixy = 1, nxy
            iasy=Pin(ixy)%iasy
            AsyPhivolt(iasy,ig)=AsyPhivolt(iasy,ig)+cellphibdry(1,ixy,2,ig)
            !AsyPhi2gt(iasy,gidx)=AsyPhi2gt(iasy,gidx)+cellphibdry(1,ixy,2,ig)
            !AsyPhi2g(iasy,gidx)=AsyPhi2g(iasy,gidx)+cellphibdry(1,ixy,2,ig) << instead of .. / 14/05/15
            AsyPhi2g(iasy,gidx)=AsyPhi2g(iasy,gidx)+Phi(ig,ixy)
        ENDDO
        ! surflux % surfJ
        i=4;
        npins=Asyinfo(Asy(1)%asytype)%ny !-east
        DO j = 1, npins
            pidx(1,j)=j*Asyinfo(Asy(1)%asytype)%nx
            pidx(2,j)=pidx(1,j)-1
        ENDDO
        DO j = 1, npins
            surflux2g(gidx)=surflux2g(gidx)+cellphibdry(i,pidx(1,j),1,ig)  !surf#/pinidx/1:surflx/ig
            surfJ2g(gidx)=surfJ2g(gidx)+cellphibdry(i,pidx(1,j),3,ig)      !surf#/pinidx/3:surfJ/ig
        ENDDO

        i=2;
        npins=Asyinfo(Asy(2)%asytype)%ny !-west
        DO j = 1, npins
            pidx(1,j)=Asyinfo(Asy(1)%asytype)%nxy+(j-1)*Asyinfo(Asy(1)%asytype)%nx+1
            pidx(2,j)=pidx(1,j)+1
        ENDDO
        DO j = 1, npins
            surflux2gt(gidx)=surflux2gt(gidx)+cellphibdry(i,pidx(1,j),1,ig)  !surf#/pinidx/1:surflx/ig
            surfJ2gt(gidx)=surfJ2gt(gidx)+cellphibdry(i,pidx(1,j),3,ig)      !surf#/pinidx/3:surfJ/ig
        ENDDO
    ENDDO
    DO iasy = 1, nxya
        DO gidx = 1,2
            !AsyPhi2g(iasy,gidx)=AsyPhi2g(iasy,gidx)/AsyInfo(iasy)%nxy
            AsyPhi2g(iasy,gidx)=AsyPhi2g(iasy,gidx)/AsyInfo(Core%Asy(iasy)%asytype)%nxy !14/05/15 << index difference
        ENDDO
    ENDDO
    DO gidx = 1,2
        surflux2g(gidx)=surflux2g(gidx)/Asyinfo(Asy(1)%asytype)%ny
        surfJ2g(gidx)=surfJ2g(gidx)/Asyinfo(Asy(1)%asytype)%ny
    ENDDO


    mod_ratio=0
    DO iasy = 1, nxya
        DO gidx = 1, 2
            mod_ratio=mod_ratio+asyphi2g(iasy,gidx)
        ENDDO
    ENDDO
    mod_ratio=(fphi+tphi)/mod_ratio
    !p_ratio=nxya*nTracerCntl%PowerCore/(isoMacXs2g(0,6,1)*fphi+isoMacXs2g(0,6,2)*tphi)!*1.0E6
    p_ratio=nTracerCntl%PowerCore/(isoMacXs2g(0,6,1)*fphi+isoMacXs2g(0,6,2)*tphi)!*1.0E6 ! 14/05/15 nxya already multiplied in PowerCore
    fphi=fphi*p_ratio; tphi=tphi*p_ratio
    p_ratio=p_ratio*mod_ratio
    DO gidx = 1,2
        DO iasy = 1, nxya
            Asyphi2g(iasy,gidx)=Asyphi2g(iasy,gidx)*p_ratio
        ENDDO
        Surflux2g(gidx)=Surflux2g(gidx)*p_ratio
        SurfJ2g(gidx)=SurfJ2g(gidx)*p_ratio
    ENDDO

    !--------pin power & flux
    phi2g=0
    nx=nx/2
    DO ig = 1, ng
        IF (ig .LT. EIdx) THEN
            gidx=1
        ELSE
            gidx=2
        ENDIF
        DO iasy = 1, 2
            DO j = 1, ny
                DO i = 1, nx
                    ixy = i+(j-1)*nx+(iasy-1)*nx*ny
                    k = i+(j-1)*nx*2+(iasy-1)*nx
                    !phi2g(gidx,k) = phi2g(gidx,k) + cellphibdry(1,ixy,2,ig)*p_ratio
                    phi2g(gidx,k) = phi2g(gidx,k) + Phi(ig,ixy) ! << 14/05/15 : to print 2G_PinFlux in more than 2 ASY
                ENDDO
            ENDDO
        ENDDO
    ENDDO
    nx=nx*2
    !--- Reflector XS modification
!D, a, r, f, nf, kf
!ALLOCATE(isoMacXs(nisotot,6,ng), isoMacSm(nisotot,ng,ng),isoMacXs2g(nisotot,6,2))
    iasy = 2
    sigr=(asyisoMacXs2g(iasy,0,3,1)-asyisoMacXs2g(iasy,0,2,1))-(asyisoMacXs2g(iasy,0,3,2)-asyisoMacXs2g(iasy,0,2,2))*asyPhi2g(iasy,2)/asyPhi2g(iasy,1)
    D1=asyisoMacXs2g(iasy,0,1,1)
    D2=asyisoMacXs2g(iasy,0,1,2)
    w1=SQRT((asyisoMacXs2g(iasy,0,2,1)+sigr)/D1)
    w2=SQRT((asyisoMacXs2g(iasy,0,2,2))/D2)
    u=Core%AsyCentX(2)-Core%AsyCentX(1)
    f1=D1*asyphi2g(iasy,1)*w1*TANH(w1*u)/surfJ2g(1)
    f2=D2*asyphi2g(iasy,2)/( sigr*surfJ2g(1)/(D1*(w2*w2-w1*w1))*(1/(w1*TANH(w1*u))-1/(w2*TANH(w2*u))) + surfJ2g(2)/(w2*TANH(w2*u)) )
ENDIF
CONTINUE




!---------------PDQ file edit------------------------------------------

io2=41
noString=.TRUE.
noString=.FALSE.

WRITE(fn2,'(A,A)') TRIM(caseid),'_PDQ.xslib'
OPEN(unit=io2, file = fn2, status = 'replace')

WRITE(io2,'(a)') '! 2G_PinFlux'
WRITE(io2,'(i4,i4)') nx, 360
DO ig = 1, 2
    CALL FluxNormalize(Phi2g(ig,:),nxy,rphi)  ! this is the replacement when there is no ADF calc.
    !CALL FluxNormalize(CellAvgPhi2g(:,ig),nxy,rphi) !  14/01/27 speed up test

    WRITE(io2,'(a4,i1)') '! g=', ig
    DO i = 1, ny
        j=(i-1)*nx+1
        WRITE(io2,'(t2,1p,100e12.5)') Phi2g(ig,j:j+nx-1)
        !WRITE(io2,'(t2,1p,100e12.5)') cellavgphi2g(j:j+nx-1,ig)
    ENDDO
ENDDO

WRITE(io2,'(a12)') '! Power_Dist'
DO i = 1, ny
    j=(i-1)*nx+1
    WRITE(io2,'(t2,1p,100e12.5)') pwr(j:j+nx-1)
ENDDO
WRITE(io2,'(a13)') '! BU-Exposure'
DO i = 1, ny
    j=(i-1)*nx+1
    WRITE(io2,'(t2,1p,100e12.5)') buexp(j:j+nx-1)
ENDDO


IF( lSA )THEN
    WRITE(io2,'(a)') '! 2G_ADF'
    !avgADF2g=1 !adf off
    !WRITE(io2,'(t2,100f12.5)') avgADF2g(1), avgADF2g(2)
    DO ig=1,2
        WRITE(io2,'(t2,i3,100f12.5)') ig, ADF2g(1:4,ig) !--- 14/10/06 edit for non-symmetric assembly : SWNE
    ENDDO
    WRITE(io2,'(a)') '! 2G_SurFlx'
    avgbphi2g(1)=fphi
    avgbphi2g(2)=tphi
    WRITE(io2,'(t2,1p,100e12.5)') avgbphi2g(1), avgbphi2g(2)
    WRITE(io2,'(a)') '! 2G_ConFlx'
    conphi2g(:,1)=fphi
    conphi2g(:,2)=tphi
    WRITE(io2,'(t2,1p,100e12.5)') conphi2g(1,1), conphi2g(2,1), conphi2g(3,1)
    WRITE(io2,'(t2,1p,100e12.5)') conphi2g(1,2), conphi2g(2,2), conphi2g(3,2)
ELSE
    WRITE(io2,'(a)') '! 2G_asy_flux'
    DO gidx = 1, 2
        WRITE(io2,'(t2,1p,100e12.5)') Asyphi2g(:,gidx)
    ENDDO
    WRITE(io2,'(a)') '! 2G_asy_surflux'
    WRITE(io2,'(t2,1p,100e12.5)') surflux2g(1:2)
    WRITE(io2,'(a)') '! 2G_asy_surfJ'
    WRITE(io2,'(t2,1p,100e12.5)') surfJ2g(1:2)
    IF( lCB )THEN
        WRITE(io2,'(a)') '! 47G_asy_flux'
        DO ig = 1, ng
            WRITE(io2,'(t2,1p,100e12.5)') AsyphiVol(:,ig)
        ENDDO
        WRITE(io2,'(a)') '! 47G_asy_nusig_f'
        DO ig = 1, ng
            WRITE(io2,'(t2,1p,100e12.5)') asyisoMacXs(:,0,5,ig)
        ENDDO
        WRITE(io2,'(a)') '! 47G_asy_sig_a'
        DO ig = 1, ng
            WRITE(io2,'(t2,1p,100e12.5)') asyisoMacXs(:,0,2,ig)
        ENDDO

    ENDIF
ENDIF

WRITE(io2,'(a)') '! 2G_flux'
WRITE(io2,'(t2,1p,100e12.5)') fphi, tphi

!D, a, r, f, nf, kf
!ALLOCATE(isoMacXs(nisotot,6,ng), isoMacSm(nisotot,ng,ng),isoMacXs2g(nisotot,6,2))
!tr, a, r, f, nu, k
!ALLOCATE(isoMicXs(nisotot,6,ng), isoMicSm(nisotot,ng,ng),isoMicXs2g(nisotot,6,2))
WRITE(io2,'(a)') '! 2G_MAC'
#ifdef totalxs
WRITE(io2,'(t2,1p,100e12.5)') isoMacXs2g(0,:,1)
WRITE(io2,'(t2,1p,100e12.5)') isoMacXs2g(0,:,2)
#else
WRITE(io2,'(t2,1p,100e12.5)') isoMacXs2g(0,1:6,1)
WRITE(io2,'(t2,1p,100e12.5)') isoMacXs2g(0,1:6,2)
#endif
IF( .NOT. lSA )THEN
    WRITE(io2,'(a)') '! TA_ADF'
    WRITE(io2,'(t2,1p,100e12.5)') sigr, w1, w2, f1, f2
    DO iasy = 1, nxya
        WRITE(io2,'(a14,i2)') '! 2G_MAC_asy# ', iasy
        WRITE(io2,'(t2,1p,100e12.5)') asyisoMacXs2g(iasy,0,:,1)
        WRITE(io2,'(t2,1p,100e12.5)') asyisoMacXs2g(iasy,0,:,2)
    ENDDO
    DO iso = 1, nisotot
        !WRITE(io2,'(i5,1p,e12.5)') isoname(iso),isonumden(iso)
        !WRITE(io2,'(a)') '! isoMacSm'
        DO ig = 1, ng
            !WRITE(io2,'(t2,1p,100e12.5)') isoMacSm(iso,ig,:)
        ENDDO
        !WRITE(io2,'(a)') '! isoMacXsSm' !--- ??? 14/01/29
        !DO ig = 1, ng
        !    WRITE(io2,'(t2,1p,100e12.5)') isoMacXsSm(iso,ig,:)
        !ENDDO
    ENDDO
ENDIF
WRITE(io2,'(a)') '! k_inf'
WRITE(io2,'(t2,100f12.9)') kinf
WRITE(io2,'(a)') '! Msquare'
WRITE(io2,'(t2,1p,100e12.5)') Msq
WRITE(io2,'(a)') '! Bsquare'
WRITE(io2,'(t2,1p,100e12.5)') Bsq

WRITE(io2,'(a)') '! 2G_Mic'
WRITE(io2,'(i4)') nisotot
DO iso = 1, nisotot
    WRITE(io2,'(i5,1p,e12.5)') isoname(iso),isonumden(iso)
#ifdef totalxs
    WRITE(io2,'(t2,1p,100e12.5)') isoMicXs2g(iso,:,1)
    WRITE(io2,'(t2,1p,100e12.5)') isoMicXs2g(iso,:,2)
#else
    WRITE(io2,'(t2,1p,100e12.5)') isoMicXs2g(iso,1:6,1)
    WRITE(io2,'(t2,1p,100e12.5)') isoMicXs2g(iso,1:6,2)
#endif
ENDDO
IF( .NOT. lSA )THEN
    DO iasy = 1, nxya
        WRITE(io2,'(a14,i2)') '! 2G_Mic_asy# ', iasy
        !---non-exist isotope deletion
        neiso=0
        DO iso = 1, nisotot
            IF( asyisonumden(iasy,iso) .EQ. 0 )THEN
                neiso=neiso+1
            ENDIF
        ENDDO
        WRITE(io2,'(i4)') nisotot-neiso
        DO iso = 1, nisotot
            IF( asyisonumden(iasy,iso) .NE. 0 )THEN
            WRITE(io2,'(i5,1p,e12.5)') isoname(iso),asyisonumden(iasy,iso)
            WRITE(io2,'(t2,1p,100e12.5)') asyisoMicXs2g(iasy,iso,:,1)
            WRITE(io2,'(t2,1p,100e12.5)') asyisoMicXs2g(iasy,iso,:,2)
            ENDIF
        ENDDO
    ENDDO
ENDIF

WRITE(io2,'(a)') '.'  !---END of official PDQ EDIT LINE -----------------------------------------------
#define editxs
#ifdef editxs
io2=42 !additionally added files to FBXS function

WRITE(fn2,'(A,A)') TRIM(caseid),'_FBX.xsl'
OPEN(unit=io2, file = fn2, status = 'replace')

IF( .NOT. lSA )THEN
WRITE(io2,'(t2,1p,100e12.5)') kinf
DO iasy = 1, nxya
    WRITE(io2,'(a14)') ' '
    WRITE(io2,'(t2,1p,100e12.5)') asyisoMacXs2g(iasy,0,:,1)
    WRITE(io2,'(t2,1p,100e12.5)') asyisoMacXs2g(iasy,0,:,2)
ENDDO
!WRITE(io2,'(a)') '! 2G_PinFlux'
!WRITE(io2,'(a)') '! 47G_flux'
DO ig = 1, ng
    WRITE(io2,'(t2,1p,100e12.5)') AsyphiVol(:,ig) !, isoMacXs(0,5,ig), isoMacXs(0,2,ig)
ENDDO

!WRITE(io2,'(a)') '! 47G_Mac'
!DO ig = 1, ng
!    WRITE(io2,'(t2,1p,100e13.5)') isoMacXs(0,:,ig)
!ENDDO
!DO ig = 1, ng
!    WRITE(io2,'(t2,1p,100e13.5)') isoMacSm(0,ig,:)
!ENDDO

    DO iasy = 1, nxya
!        WRITE(io2,'(a14,i2)') '! 47G_Mac_asy# ', iasy
        DO ig = 1, ng
            WRITE(io2,'(t2,1p,100e13.5)') asyisoMacXs(iasy,0,:,ig)
        ENDDO
        DO ig = 1, ng
            WRITE(io2,'(t2,1p,100e13.5)') asyisoMacSm(iasy,0,ig,:)
        ENDDO
    ENDDO
ENDIF


#endif


#ifdef printf2
DO iso = 1, nisotot
    WRITE(io2,'(i5,1p,e12.5)') isoname(iso),isonumden(iso)
    DO ig = 1, ng
        WRITE(io2,'(t2,1p,100e12.5)') isoMacXs(iso,:,ig)
    ENDDO
ENDDO

WRITE(io2,'(a)') '! 47G_Mic'
DO iso = 1, nisotot
    WRITE(io2,'(i5,1p,e12.5)') isoname(iso),isonumden(iso)
    DO ig = 1, ng
        WRITE(io2,'(t2,1p,100e12.5)') isoMicXs(iso,:,ig)
    ENDDO
ENDDO


!#ifdef printf2
WRITE(io2,'(a)') ' 2group flux distribution (normalized)'
DO ig = 1, 2
    WRITE(io2,'(a4,i1)') 'g = ', ig
    DO i = 1, ny
        j=(i-1)*nx+1
        WRITE(io2,'(t2,1p,100e12.5)') phi2g(ig,j:j+nx-1)
    ENDDO
ENDDO
WRITE(io2,'(a)') ' 2group boundary flux distribution '
DO k = 1, 4 !direction
    WRITE(io2,'(a4,i1)') 'd = ', k
    DO ig = 1, 2
        WRITE(io2,'(a6,i1)') '  g = ', ig
        DO i = 1, ny
            j=(i-1)*nx+1
            WRITE(io2,'(t2,1p,100e12.5)') cellbphi2g(k,j:j+nx-1,ig)
        ENDDO
    ENDDO
ENDDO


WRITE(io2,'(a)') ' '
WRITE(io2,'(a)') ' '
WRITE(io2,'(a)') ' '
WRITE(io2,'(a)') ' 2group boundary current distribution '
DO k = 1, 4 !direction
    WRITE(io2,'(a4,i1)') 'd = ', k
    DO ig = 1, 2
        WRITE(io2,'(a6,i1)') '  g = ', ig
        DO i = 1, ny
            j=(i-1)*nx+1
            WRITE(io2,'(t2,1p,100e12.5)') cellbJ2g(k,j:j+nx-1,ig)
        ENDDO
    ENDDO
ENDDO
WRITE(io2,'(a)') ' '

io2=42
noString=.TRUE.
WRITE(fn2,'(A,A,A)') 'H_',TRIM(caseid),'_2G_noStr.xslib'
OPEN(unit=io2, file = fn2, status = 'replace')
IF(.NOT. noString) WRITE(io2,'(a)') ' base_micro 1'
WRITE(io2,'(t2,1p,6e12.5)') FRR(1),FRR(2),FRR(3),FRR(3),FRR(7),FRR(4)
WRITE(io2,'(t2,1p,6e12.5)') TRR(1),TRR(2),TRR(3),TRR(3),TRR(7),TRR(4)
WRITE(io2,'(t2,1p,100e12.5)') FRR(5), FRR(6)
WRITE(io2,'(t2,1p,100e12.5)') TRR(5), TRR(6)
IF(.NOT. noString) WRITE(io2,'(a)') ' 2G_ADF : SF,T  WF,T  NF,T  EF,T '
WRITE(io2,'(t2,1p,100e12.5)') ADF2g(1,1), ADF2g(1,2), ADF2g(2,1), ADF2g(2,2), ADF2g(3,1), ADF2g(3,2), ADF2g(4,1), ADF2g(4,2)
!WRITE(io2,'(a)') ' Bsquare '
!WRITE(io2,'(t2,1p,100e12.5)') Bsq
WRITE(io2,'(a)') ' '
IF(.NOT. noString) WRITE(io2,'(a)') ' 2group flux distribution'
DO ig = 1, 2
    IF(.NOT. noString) WRITE(io2,'(a4,i1)') 'g = ', ig
    DO i = 1, ny
        j=(i-1)*nx+1
        WRITE(io2,'(t2,1p,100e12.5)') cellavgphi2g(j:j+nx-1,ig)
    ENDDO
ENDDO
WRITE(io2,'(a)') ' '
IF(.NOT. noString) WRITE(io2,'(a)') ' 2group boundary flux distribution '
DO k = 1, 4 !direction
    IF(.NOT. noString) WRITE(io2,'(a4,i1)') 'd = ', k
    DO ig = 1, 2
        IF(.NOT. noString) WRITE(io2,'(a6,i1)') '  g = ', ig
        DO i = 1, ny
            j=(i-1)*nx+1
            WRITE(io2,'(t2,1p,100e12.5)') cellbphi2g(k,j:j+nx-1,ig)
        ENDDO
    ENDDO
ENDDO
WRITE(io2,'(a)') ' '
IF(.NOT. noString) WRITE(io2,'(a)') ' 2group boundary current distribution '
DO k = 1, 4 !direction
    IF(.NOT. noString) WRITE(io2,'(a4,i1)') 'd = ', k
    DO ig = 1, 2
        IF(.NOT. noString) WRITE(io2,'(a6,i1)') '  g = ', ig
        DO i = 1, ny
            j=(i-1)*nx+1
            WRITE(io2,'(t2,1p,100e12.5)') cellbJ2g(k,j:j+nx-1,ig)
        ENDDO
    ENDDO
ENDDO
WRITE(io2,'(a)') ' '
IF(.NOT. noString) WRITE(io2,'(a)') ' 2group flux distribution (normalized)'
DO ig = 1, 2
    IF(.NOT. noString) WRITE(io2,'(a4,i1)') 'g = ', ig
    DO i = 1, ny
        j=(i-1)*nx+1
        WRITE(io2,'(t2,1p,100e12.5)') phi2g(ig,j:j+nx-1)
    ENDDO
ENDDO
#endif




END SUBROUTINE
#endif


SUBROUTINE FluxNormalize(phi,nxy,rphi)
    IMPLICIT NONE
    REAL :: sumphi, rphi
    INTEGER :: nxy, ixy
    REAL :: phi(nxy)
    sumphi=0
    DO ixy = 1, nxy
        sumphi=sumphi+phi(ixy)
    ENDDO
    rphi=nxy/sumphi
    DO ixy = 1, nxy
        phi(ixy)=phi(ixy)*rphi
    ENDDO
    sumphi=0
    DO ixy = 1, nxy
        sumphi=sumphi+phi(ixy)
    ENDDO
    IF( abs(sumphi-nxy) .GT. 1E-3 )THEN
        continue
    ENDIF
ENDSUBROUTINE

SUBROUTINE FluxNormalizePower(phi,nxy,power,rphi)
    IMPLICIT NONE
    REAL :: sumphi, power, rphi
    INTEGER :: nxy, ixy
    REAL :: phi(nxy)
    sumphi=0
    DO ixy = 1, nxy
        sumphi=sumphi+phi(ixy)
    ENDDO
    rphi=nxy/sumphi
    DO ixy = 1, nxy
        phi(ixy)=phi(ixy)*rphi
    ENDDO
    sumphi=0
    DO ixy = 1, nxy
        sumphi=sumphi+phi(ixy)
    ENDDO
    IF( abs(sumphi-nxy) .GT. 1E-3 )THEN
        continue
    ENDIF
ENDSUBROUTINE
