SUBROUTINE GCInit(Core, CmInfo, ng)
    USE GC_mod
    USE GCpin_mod
    USE TYPEDEF,         ONLY : CoreInfo_Type, CMInfo_Type
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(CMInfo_Type) :: CMInfo
    INTEGER :: ng

    INTEGER :: ig
    !--- ALLOCATE ---
    ALLOCATE(Phi(ng), PhiVol(ng), PinPhiVol(ng,Core%nxy))
    ALLOCATE(GrpIdx(ng), PhiG(ngrp))
    ALLOCATE(FsrPhi(Core%nCoreFsr,ng), FsrPhiVol(Core%nCoreFsr,ng), FsrVol(Core%nCoreFsr))

    ALLOCATE(IsoMacXsSm(isosize,ng,ng))
    ALLOCATE(IsoMacXsSmP1(isosize,ng,ng)) ! S1 scattering Matrix 14/06/01

    !For B1
    ALLOCATE(f_ratio(ng), phicrit(ng), Dng(ng), XstrB1(ng), Xstr(ng), Xst(ng))
    ALLOCATE(ADF(4,ng),ADF2g(4,ngrp), avgADF2g(ngrp), ADF2GSurf(4,ngrp))
    ALLOCATE(SurfJ(0:2, 4, ng),SurfJG(0:2, 4, ngrp), SurfPhiG(4, ngrp), ChiG(ngrp))
    ALLOCATE(avgBPhiG(ngrp), conPhiG(3,ngrp))

    !--- from Core
    AsyInfo => Core%AsyInfo; Asy => Core%Asy
    Pin => Core%Pin; Cell => Core%CellInfo
    nxya=Core%nxya ! n-assemblies
    nCorexy=Core%nxy !nxy pin
    nCoreFsr=Core%nCoreFsr;
    nz=Core%nz

    !--- from CmfdInfo
    PinXS => CMInfo%PinXS; PHIC => CmInfo%PhiC

    !--- initialize
    checkiso=.FALSE. ! 0=not Used. 1=used
    nisotot=0; isoNumden=0; h2oNumden=0; uo2Numden=0;
    XSt=0.0; !Sm1=0.0;
    nCoreFxr = 0
    Phi=0; PhiVol=0; PinPhiVol=0;


    Dng=0;    Bsq=0;  Msq=0
    kinf=Eigv
    DO ig = 1, ng
        f_ratio(ig)=1
    ENDDO

    !--- Check whether it is SA, TA, CB or else
    lSA=.TRUE.
    lZdir=.FALSE.
    lCB=.FALSE.
    IF( nz .NE. 1 )THEN
        lZdir=.TRUE.
    ENDIF
    IF( nxya .NE. 1 )THEN
        lSA=.FALSE.
        IF( nxya .EQ. 4 )THEN
            lCB=.TRUE.
        ENDIF
    ENDIF

    lPinwise=.FALSE.
    lPin=.FALSE.


    ALLOCATE(pinRmv(nCorexy,ng),pinFis(nCorexy,ng),pinSS(nCorexy,ng,ng),pinR(nCorexy,ng),pinFlux(nCorexy,ng),pinJ(nCorexy,ng))
    ALLOCATE(pinSin(nCorexy,ng),pinSout(nCorexy,ng),pinAbs(nCorexy,ng))
    ALLOCATE(pinSinG(nCorexy,ng),pinSoutG(nCorexy,ng))
    ALLOCATE(pinJdir(nCorexy,ng,4))
    !ALLOCATE(SrcR(nCorexy,ng),lossr(nCorexy,ng))
    pinrmv=0;pinfis=0;pinss=0;pinr=0;pinFlux=0;pinJ=0;
    pinSin=0;pinSout=0;pinAbs=0;
    pinSinG=0;pinSoutG=0;
    pinJdir=0;
    
    u238n2n=0.;u238phi=0.0

ENDSUBROUTINE

SUBROUTINE GCBU(FmInfo)
    USE TYPEDEF,         ONLY : FmInfo_Type
    USE GC_mod
    IMPLICIT NONE
    TYPE(FMInfo_Type) :: FmInfo

    INTEGER :: ix, iy, ixy
    INTEGER :: ifxr, fxridxst, nLocalFxr
    REAL :: areasum

    ALLOCATE(buexp(nx, ny))
    DO ix = 1, nx
        DO iy = 1, ny
            ixy = AsyInfo(iasytype)%Pin2DIdx(ix,iy)
            ixy = Asy(iasy)%GlobalPinIdx(ixy)
            areasum=0
            buexp(ix,iy)=0
            fxridxst=Pin(ixy)%fxridxst
            nLocalFxr=Pin(ixy)%nFXRmax
            DO ifxr = fxridxst, fxridxst+nlocalfxr-1
                IF( Fminfo%FXR(ifxr,iz)%lfuel )THEN
                    areasum=areasum+Fminfo%FXR(ifxr,iz)%area
                    buexp(ix,iy)=buexp(ix,iy)+Fminfo%FXR(ifxr,iz)%burnup*Fminfo%FXR(ifxr,iz)%area
                ENDIF
            ENDDO
            IF( areasum .NE. 0 )THEN
                buexp(ix,iy)=buexp(ix,iy)/areasum
            ELSE
                buexp(ix,iy)=0
            ENDIF
        ENDDO
    ENDDO
ENDSUBROUTINE

SUBROUTINE GCPWR(Core, CmInfo, nTracerCntl, PE, ng)
    USE TYPEDEF,         ONLY : CoreInfo_Type, CMInfo_Type, PE_TYPE
    USE CNTL,         ONLY : nTracerCntl_Type
    USE GC_mod
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(CMInfo_Type) :: CMInfo
    TYPE(nTracerCntl_Type) :: nTracerCntl
    TYPE(PE_TYPE) :: PE
    INTEGER :: ng

    INTEGER :: ix, iy, ixy

    ALLOCATE(pwr(nx, ny))

    pwr=0
    IF( lzDir )THEN
        DO ix = 1, nx
            DO iy = 1, ny
                ixy=AsyInfo(iasytype)%Pin2DIdx(ix, iy)
                pwr(ix,iy)=powerdist%pinpower3d(ixy,iasy,iz)  !!!CHECK IX, IY ORDER
            ENDDO
        ENDDO
    ELSE
        DO ix = 1, nx
            DO iy = 1, ny
                ixy=AsyInfo(iasytype)%Pin2DIdx(ix, iy)
                pwr(ix,iy)=powerdist%pinpower2d(ixy,iasy)
            ENDDO
        ENDDO
    ENDIF


ENDSUBROUTINE

SUBROUTINE GCISO_fast(Core, FmInfo, nTracerCntl, ng)
    !GCIso : Search the list of isotope in FSR
    USE TYPEDEF,         ONLY : CoreInfo_Type, FmInfo_Type
    USE GC_mod
    USE GCpin_mod
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FMInfo_Type) :: FmInfo
    TYPE(nTracerCntl_Type) :: nTracerCntl
    INTEGER :: ng

    INTEGER :: ig
    INTEGER :: j, k

    INTEGER :: iso, isoidx, id, idiso, nisoinFxr
    INTEGER :: ixy, ixyl, ix, iy ! local ixy
    INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nLocalFsr, localFsr, iCel, iFxr, ifsr, nFsrInFxr
    REAL :: myVol, myfxrvol
    REAL :: FsrVolsum, RFsrVolsum
    LFP_ID=38
    isolist(:,1)=LFP_ID;
    isoidx=0;
    isoidx=isoidx+1;  idiso=92234; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 1
    isoidx=isoidx+1;  idiso=92235; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 2
    isoidx=isoidx+1;  idiso=92236; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 3
    isoidx=isoidx+1;  idiso=92237; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 4
    isoidx=isoidx+1;  idiso=92238; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 5
    isoidx=isoidx+1;  idiso=93237; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 6
    isoidx=isoidx+1;  idiso=93238; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 7
    isoidx=isoidx+1;  idiso=93239; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 8
    isoidx=isoidx+1;  idiso=94238; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 9
    isoidx=isoidx+1;  idiso=94239; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 10
    isoidx=isoidx+1;  idiso=94240; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 11
    isoidx=isoidx+1;  idiso=94241; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 12
    isoidx=isoidx+1;  idiso=94242; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 13
    isoidx=isoidx+1;  idiso=95241; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 14
    isoidx=isoidx+1;  idiso=95342; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 15
    isoidx=isoidx+1;  idiso=95243; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 16
    isoidx=isoidx+1;  idiso=96242; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 17
    isoidx=isoidx+1;  idiso=96243; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 18
    isoidx=isoidx+1;  idiso=96244; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 19
    isoidx=isoidx+1;  idiso=60647; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 20
    isoidx=isoidx+1;  idiso=61647; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 21
    isoidx=isoidx+1;  idiso=61648; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 22
    isoidx=isoidx+1;  idiso=61748; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 23
    isoidx=isoidx+1;  idiso=61649; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 24
    isoidx=isoidx+1;  idiso=62649; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 25
    isoidx=isoidx+1;  idiso=53635; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 26
    isoidx=isoidx+1;  idiso=54635; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 27
    isoidx=isoidx+1;  idiso=64152; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 28
    isoidx=isoidx+1;  idiso=64154; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 29
    isoidx=isoidx+1;  idiso=64155; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 30
    isoidx=isoidx+1;  idiso=64156; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 31
    isoidx=isoidx+1;  idiso=64157; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 32
    isoidx=isoidx+1;  idiso=64158; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 33
    isoidx=isoidx+1;  idiso=64160; isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 34
    isoidx=isoidx+1;  idiso=5000;  isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! 35
    isoidx=isoidx+1;  idiso=8016;  isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! H2O 36
    isoidx=isoidx+1;  idiso=8016;  isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx; isolist(id,2)=id; ! UO2 37
    isoidx=isoidx+1;  idiso=0;     isoname(isoidx)=idiso;                                                            ! LFP 38
    nisotot=isoidx;                
    isoidx=isoidx+1;  idiso=8001;  isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx-2; isolist(id,2)=id; ! H2O 37
    isoidx=isoidx+1;  idiso=8001;  isoname(isoidx)=idiso; id=MapNucl(idiso); isolist(id,1)=isoidx-2; isolist(id,2)=id; ! UO2 38
    
    FsrVolsum=0;
    !--- Check Isotope LIST
    IF( .NOT. lPin )THEN
        IF( nTracerCntl%lXsLib )THEN
        DO ixyl = 1, nxy
            ixy = Asy(iasy)%GlobalPinIdx(ixyl) !global Pin Index
            icel = Pin(ixy)%Cell(iz)
            FxrIdxSt = Pin(ixy)%FxrIdxSt
            FsrIdxSt = Pin(ixy)%FsrIdxSt
            nlocalFxr = Cell(icel)%nFxr
            DO j = 1, nLocalFxr
                nCoreFxr=nCoreFxr+1
                ifxr = FxrIdxSt + j -1
                myFxr => FmInfo%FXR(ifxr,iz)
                nFsrInFxr = myFxr%nFsrInFxr
                nisoInFxr = myFxr%niso
                myfxrvol = myFxr%area*hz
                DO k = 1, nFsrInFxr
                    ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                    localFsr=core%cellinfo(icel)%mapfxr2fsridx(k,j)
                    FsrVol(ifsr)=Cell(icel)%vol(localFsr)*hz
                    myVol=FsrVol(ifsr)
                    DO ig = 1, ng
                        FsrPhi(ifsr,ig)=FmInfo%phis(ifsr,iz,ig)
                        FsrPhiVol(ifsr,ig)=FmInfo%phis(ifsr,iz,ig)*myVol
                        PhiVol(ig)=PhiVol(ig)+FsrPhiVol(ifsr,ig)
                    ENDDO
                    FsrVolsum=FsrVolsum+myVol
                ENDDO
                IF( myFxr%lFuel )THEN !--- O-H
                    DO iso = 1, nisoInFxr
                        idiso=myFxr%idiso(iso)
                        id = MapNucl(idiso)  !id = Map(92235) = 37
                        !isoidx=isoList(id,1)
                        !IF( idiso .EQ. 8016 )THEN
                        !    isoidx=isoidx+1
                        !    isoNumden(isoidx)=isoNumden(isoidx)+myfxrvol*myFxr%pnum(iso)
                        !    uo2Numden=uo2Numden+myfxrvol*myFxr%pnum(iso)
                        !ELSE
                        !    isoNumden(isoidx)=isoNumden(isoidx)+myfxrvol*myFxr%pnum(iso)
                        !ENDIF
                    ENDDO
                ELSE
                    DO iso = 1, nisoInFxr
                        idiso=myFxr%idiso(iso)
                        id = MapNucl(idiso)  !id = Map(92235) = 37
                        !isoidx=isoList(id,1)
                        !isoNumden(isoidx)=isoNumden(isoidx)+myfxrvol*myFxr%pnum(iso)
                        !IF( idiso .EQ. 8016 )THEN
                        !    h2oNumden=h2oNumden+myfxrvol*myFxr%pnum(iso)
                        !ENDIF
                    ENDDO
                ENDIF
            ENDDO
        ENDDO
        ELSE ! macro
            DO ixyl = 1, nxy
                ixy = Asy(iasy)%GlobalPinIdx(ixyl) !global Pin Index
                icel = Pin(ixy)%Cell(iz)
                FxrIdxSt = Pin(ixy)%FxrIdxSt
                FsrIdxSt = Pin(ixy)%FsrIdxSt
                nlocalFxr = Cell(icel)%nFxr
                DO j = 1, nLocalFxr
                    nCoreFxr=nCoreFxr+1
                    ifxr = FxrIdxSt + j -1
                    myFxr => FmInfo%FXR(ifxr,iz)
                    nFsrInFxr = myFxr%nFsrInFxr
                    nisoInFxr = myFxr%niso
                    myfxrvol = myFxr%area*hz
                    DO k = 1, nFsrInFxr
                        ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                        localFsr=core%cellinfo(icel)%mapfxr2fsridx(k,j)
                        FsrVol(ifsr)=Cell(icel)%vol(localFsr)*hz
                        myVol=FsrVol(ifsr)
                        DO ig = 1, ng
                            FsrPhi(ifsr,ig)=FmInfo%phis(ifsr,iz,ig)
                            FsrPhiVol(ifsr,ig)=FmInfo%phis(ifsr,iz,ig)*myVol
                            PhiVol(ig)=PhiVol(ig)+FsrPhiVol(ifsr,ig)
                        ENDDO
                        FsrVolsum=FsrVolsum+myVol
                    ENDDO
                ENDDO
            ENDDO
        ENDIF
    ELSE
IF( nTracerCntl%lXsLib )THEN
      DO iy = yst, yed
        DO ix = xbg, xed
            ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy)
            ixy = Asy(iasy)%GlobalPinIdx(ixyl) !global Pin Index
            icel = Pin(ixy)%Cell(iz)
            FxrIdxSt = Pin(ixy)%FxrIdxSt
            FsrIdxSt = Pin(ixy)%FsrIdxSt
            nlocalFxr = Cell(icel)%nFxr
            DO j = 1, nLocalFxr
                nCoreFxr=nCoreFxr+1
                ifxr = FxrIdxSt + j -1
                myFxr => FmInfo%FXR(ifxr,iz)
                nFsrInFxr = myFxr%nFsrInFxr
                nisoInFxr = myFxr%niso
                myfxrvol = myFxr%area*hz
                DO k = 1, nFsrInFxr
                    ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                    localFsr=core%cellinfo(icel)%mapfxr2fsridx(k,j)
                    FsrVol(ifsr)=Cell(icel)%vol(localFsr)*hz
                    myVol=FsrVol(ifsr)
                    DO ig = 1, ng
                        FsrPhi(ifsr,ig)=FmInfo%phis(ifsr,iz,ig)
                        FsrPhiVol(ifsr,ig)=FmInfo%phis(ifsr,iz,ig)*myVol
                        PhiVol(ig)=PhiVol(ig)+FsrPhiVol(ifsr,ig)
                    ENDDO
                    FsrVolsum=FsrVolsum+myVol
                ENDDO
                IF( myFxr%lFuel )THEN !--- O-H
                    DO iso = 1, nisoInFxr
                        idiso=myFxr%idiso(iso)
                        id = MapNucl(idiso)  !id = Map(92235) = 37
                        isoidx=isoList(id,1)
                        IF( idiso .EQ. 8016 )THEN
                            isoidx=isoidx+1
                            isoNumden(isoidx)=isoNumden(isoidx)+myfxrvol*myFxr%pnum(iso)
                            uo2Numden=uo2Numden+myfxrvol*myFxr%pnum(iso)
                        ELSE
                            isoNumden(isoidx)=isoNumden(isoidx)+myfxrvol*myFxr%pnum(iso)
                        ENDIF
                    ENDDO
                ELSE
                    DO iso = 1, nisoInFxr
                        idiso=myFxr%idiso(iso)
                        id = MapNucl(idiso)  !id = Map(92235) = 37
                        isoidx=isoList(id,1)
                        isoNumden(isoidx)=isoNumden(isoidx)+myfxrvol*myFxr%pnum(iso)
                        IF( idiso .EQ. 8016 )THEN
                            h2oNumden=h2oNumden+myfxrvol*myFxr%pnum(iso)
                        ENDIF
                    ENDDO
                ENDIF
            ENDDO  ! FXR
        ENDDO ! IX
    ENDDO ! IY
ELSE ! macro
    DO iy = yst, yed
        DO ix = xbg, xed
            ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy)
            ixy = Asy(iasy)%GlobalPinIdx(ixyl) !global Pin Index
            icel = Pin(ixy)%Cell(iz)
            FxrIdxSt = Pin(ixy)%FxrIdxSt
            FsrIdxSt = Pin(ixy)%FsrIdxSt
            nlocalFxr = Cell(icel)%nFxr
            DO j = 1, nLocalFxr
                nCoreFxr=nCoreFxr+1
                ifxr = FxrIdxSt + j -1
                myFxr => FmInfo%FXR(ifxr,iz)
                nFsrInFxr = myFxr%nFsrInFxr
                nisoInFxr = myFxr%niso
                myfxrvol = myFxr%area*hz
                DO k = 1, nFsrInFxr
                    ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                    localFsr=core%cellinfo(icel)%mapfxr2fsridx(k,j)
                    FsrVol(ifsr)=Cell(icel)%vol(localFsr)*hz
                    myVol=FsrVol(ifsr)
                    DO ig = 1, ng
                        FsrPhi(ifsr,ig)=FmInfo%phis(ifsr,iz,ig)
                        FsrPhiVol(ifsr,ig)=FmInfo%phis(ifsr,iz,ig)*myVol
                        PhiVol(ig)=PhiVol(ig)+FsrPhiVol(ifsr,ig)
                    ENDDO
                    FsrVolsum=FsrVolsum+myVol
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDIF
ENDIF
    IF( .NOT. nTracerCntl%lXsLib ) nisotot=1
    ALLOCATE(isoMicXs(nisotot,0:nxs,ng), isoMicSm(nisotot,ng,ng),isoMicXs2g(nisotot,0:nxs,ngrp),isoMicSm2g(nisotot,ngrp,ngrp)) !tr, a, r, f, nu, k
    ALLOCATE(isoMacXs(0:nisotot,0:nxs,ng), isoMacSm(0:nisotot,ng,ng),isoMacXs2g(0:nisotot,0:nxs,ngrp),isoMacSm2g(0:nisotot,ngrp,ngrp)) !D, a, r, f, nf, kf
    isoMacXs=0;isoMacSm=0;isoMacXs2g=0;isoMicXs=0;isoMicSm=0;isoMicXs2g=0; isoMicSm2g=0; isoMacSm2g=0;
    !--- P1 variable
    ALLOCATE(isoMacSmP1(0:nisotot,ng,ng),isoMicSmP1(nisotot,ng,ng))
    ALLOCATE(isoMacSm2gP1(0:nisotot,ngrp,ngrp),isoMicSm2gP1(nisotot,ngrp,ngrp))
    isoMacSmP1=0;isoMicSmP1=0;isoMacSm2gP1=0;isoMicSm2gP1=0

    IF( nisotot .GT. isosize )THEN
        write(*,*) '                                EXCEED. BYS edit 14/02/27'
    ENDIF
    RFsrVolSum=1.0/FsrVolsum
    CALL MULTI_CA(RFsrVolsum, isoNumden(:), nisotot) !isoNumden=isoNumden/FsrVolsum ! isoNumden = avg(n)_iso
    IF( .NOT. nTracerCntl%lXsLib ) isonumden=1.0_8
    h2onumden=rFsrVolsum*h2onumden
    uo2numden=rFsrVolsum*uo2numden


ENDSUBROUTINE
SUBROUTINE GCISO_new(Core, FmInfo, nTracerCntl, ng)
    !GCIso : Search the list of isotope in FSR
    USE TYPEDEF,         ONLY : CoreInfo_Type, FmInfo_Type
    USE GC_mod
    USE GCpin_mod
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FMInfo_Type) :: FmInfo
    TYPE(nTracerCntl_Type) :: nTracerCntl
    INTEGER :: ng

    INTEGER :: ig
    INTEGER :: j, k

    INTEGER :: iso, isoidx, id, idiso, nisoinFxr
    INTEGER :: ixy, ixyl, ix, iy ! local ixy
    INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nLocalFsr, localFsr, iCel, iFxr, ifsr, nFsrInFxr
    REAL :: myVol, myfxrvol
    REAL :: FsrVolsum, RFsrVolsum

    FsrVolsum=0;
    !--- Check Isotope LIST
    IF( .NOT. lPin )THEN
        DO ixyl = 1, nxy
            ixy = Asy(iasy)%GlobalPinIdx(ixyl) !global Pin Index
            icel = Pin(ixy)%Cell(iz)
            FxrIdxSt = Pin(ixy)%FxrIdxSt
            FsrIdxSt = Pin(ixy)%FsrIdxSt
            nlocalFxr = Cell(icel)%nFxr
            !localFsr = 0
            DO j = 1, nLocalFxr
                nCoreFxr=nCoreFxr+1
                ifxr = FxrIdxSt + j -1
                myFxr => FmInfo%FXR(ifxr,iz)
                nFsrInFxr = myFxr%nFsrInFxr
                !FsrIdxSt = myFxr%FsrIdxSt
                nisoInFxr = myFxr%niso
                myfxrvol = myFxr%area*hz
                DO k = 1, nFsrInFxr
                    !ifsr = FsrIdxSt + k -1
                    !localFsr=localFsr+1
                    !--- EDIT 14/12/23 Sq. cell problem
                    ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                    localFsr=core%cellinfo(icel)%mapfxr2fsridx(k,j)
                    FsrVol(ifsr)=Cell(icel)%vol(localFsr)*hz
                    myVol=FsrVol(ifsr)
                    DO ig = 1, ng
                        FsrPhi(ifsr,ig)=FmInfo%phis(ifsr,iz,ig)
                        FsrPhiVol(ifsr,ig)=FmInfo%phis(ifsr,iz,ig)*myVol
                        !---BYS edit 15/02/23 Ryumin feedback
                        !FsrPhiVol(ifsr,ig)=fsrphir(ifsr,ig)*myVol
                        !FsrPhi(ifsr,ig)=fsrphir(ifsr,ig)
                        PhiVol(ig)=PhiVol(ig)+FsrPhiVol(ifsr,ig)
                    ENDDO
                    FsrVolsum=FsrVolsum+myVol
                ENDDO
                IF( nTracerCntl%lXsLib )THEN
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
                        isoNumden(isoidx)=isoNumden(isoidx)+myfxrvol*myFxr%pnum(iso)
                        IF( idiso .EQ. 8016 )THEN
                            IF( myFxr%lH2O )THEN
                                h2oNumden=h2oNumden+myfxrvol*myFxr%pnum(iso)
                            ELSEIF( myFxr%lFuel) THEN
                                uo2Numden=uo2Numden+myfxrvol*myFxr%pnum(iso)
                            ENDIF
                        !ELSEIF( myFxr%idiso(iso) .EQ. 8001 )THEN
                        !    uo2Numden=uo2Numden+myVol*myFxr%pnum(iso)
                        ENDIF
                    ENDDO
                ELSE

                ENDIF
            ENDDO
        ENDDO
    ELSE
      DO iy = yst, yed
        DO ix = xbg, xed
            ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy)
            ixy = Asy(iasy)%GlobalPinIdx(ixyl) !global Pin Index
            icel = Pin(ixy)%Cell(iz)
            FxrIdxSt = Pin(ixy)%FxrIdxSt
            FsrIdxSt = Pin(ixy)%FsrIdxSt
            nlocalFxr = Cell(icel)%nFxr
            !localFsr = 0
            DO j = 1, nLocalFxr
                nCoreFxr=nCoreFxr+1
                ifxr = FxrIdxSt + j -1
                myFxr => FmInfo%FXR(ifxr,iz)
                nFsrInFxr = myFxr%nFsrInFxr
                !FsrIdxSt = myFxr%FsrIdxSt
                nisoInFxr = myFxr%niso
                myfxrvol = myFxr%area*hz
                DO k = 1, nFsrInFxr
                    !ifsr = FsrIdxSt + k -1
                    !localFsr=localFsr+1
                    !--- EDIT 14/12/23 Sq. cell problem
                    ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                    localFsr=core%cellinfo(icel)%mapfxr2fsridx(k,j)
                    FsrVol(ifsr)=Cell(icel)%vol(localFsr)*hz
                    myVol=FsrVol(ifsr)
                    DO ig = 1, ng
                        FsrPhi(ifsr,ig)=FmInfo%phis(ifsr,iz,ig)
                        FsrPhiVol(ifsr,ig)=FmInfo%phis(ifsr,iz,ig)*myVol
                        !---BYS edit 15/02/23 Ryumin feedback
                        !FsrPhiVol(ifsr,ig)=fsrphir(ifsr,ig)*myVol
                        !FsrPhi(ifsr,ig)=fsrphir(ifsr,ig)
                        PhiVol(ig)=PhiVol(ig)+FsrPhiVol(ifsr,ig)
                    ENDDO
                    FsrVolsum=FsrVolsum+myVol
                ENDDO
                
                IF( nTracerCntl%lXsLib )THEN
                    DO iso = 1, nisoInFxr
                        idiso=myFxr%idiso(iso)
                        IF( idiso .EQ. 8001 )THEN
                            idiso=8016
                        ENDIF
                        id = MapNucl(idiso)  !id = Map(92235) = 37
                        !isodata => ldiso(id)
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
                        isoNumden(isoidx)=isoNumden(isoidx)+myfxrvol*myFxr%pnum(iso)
                        IF( idiso .EQ. 8016 )THEN
                            IF( myFxr%lH2O )THEN
                                h2oNumden=h2oNumden+myfxrvol*myFxr%pnum(iso)
                            ELSEIF( myFxr%lFuel) THEN
                                uo2Numden=uo2Numden+myfxrvol*myFxr%pnum(iso)
                            ENDIF
                        !ELSEIF( myFxr%idiso(iso) .EQ. 8001 )THEN
                        !    uo2Numden=uo2Numden+myVol*myFxr%pnum(iso)
                        ENDIF
                    ENDDO
                ENDIF
            ENDDO
        ENDDO
        ENDDO
    ENDIF
    IF( .NOT. nTracerCntl%lXsLib ) nisotot=1
    ALLOCATE(isoMicXs(nisotot,0:nxs,ng), isoMicSm(nisotot,ng,ng),isoMicXs2g(nisotot,0:nxs,ngrp),isoMicSm2g(nisotot,ngrp,ngrp)) !tr, a, r, f, nu, k
    ALLOCATE(isoMacXs(0:nisotot,0:nxs,ng), isoMacSm(0:nisotot,ng,ng),isoMacXs2g(0:nisotot,0:nxs,ngrp),isoMacSm2g(0:nisotot,ngrp,ngrp)) !D, a, r, f, nf, kf
    isoMacXs=0;isoMacSm=0;isoMacXs2g=0;isoMicXs=0;isoMicSm=0;isoMicXs2g=0; isoMicSm2g=0; isoMacSm2g=0;
    !--- P1 variable
    ALLOCATE(isoMacSmP1(0:nisotot,ng,ng),isoMicSmP1(nisotot,ng,ng))
    ALLOCATE(isoMacSm2gP1(0:nisotot,ngrp,ngrp),isoMicSm2gP1(nisotot,ngrp,ngrp))
    isoMacSmP1=0;isoMicSmP1=0;isoMacSm2gP1=0;isoMicSm2gP1=0

    IF( nisotot .GT. isosize )THEN
        write(*,*) '                                EXCEED. BYS edit 14/02/27'
    ENDIF
    RFsrVolSum=1.0/FsrVolsum
    CALL MULTI_CA(RFsrVolsum, isoNumden(:), nisotot) !isoNumden=isoNumden/FsrVolsum ! isoNumden = avg(n)_iso
    IF( .NOT. nTracerCntl%lXsLib ) isonumden=1.0_8
    h2onumden=rFsrVolsum*h2onumden
    uo2numden=rFsrVolsum*uo2numden


ENDSUBROUTINE

SUBROUTINE GCISO(Core, FmInfo, nTracerCntl, ng) !slow
    !GCIso : Search the list of isotope in FSR
    USE TYPEDEF,         ONLY : CoreInfo_Type, FmInfo_Type
    USE GC_mod
    USE GCpin_mod
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FMInfo_Type) :: FmInfo
    TYPE(nTracerCntl_Type) :: nTracerCntl
    INTEGER :: ng

    INTEGER :: ig
    INTEGER :: j, k

    INTEGER :: iso, isoidx, id, idiso, nisoinFxr
    INTEGER :: ixy, ixyl, ix, iy ! local ixy
    INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nLocalFsr, localFsr, iCel, iFxr, ifsr, nFsrInFxr
    REAL :: myVol
    REAL :: FsrVolsum, RFsrVolsum

    FsrVolsum=0;
    !--- Check Isotope LIST
    IF( .NOT. lPin )THEN
        DO ixyl = 1, nxy
            ixy = Asy(iasy)%GlobalPinIdx(ixyl) !global Pin Index
            icel = Pin(ixy)%Cell(iz)
            FxrIdxSt = Pin(ixy)%FxrIdxSt
            FsrIdxSt = Pin(ixy)%FsrIdxSt
            nlocalFxr = Cell(icel)%nFxr
            !localFsr = 0
            DO j = 1, nLocalFxr
                nCoreFxr=nCoreFxr+1
                ifxr = FxrIdxSt + j -1
                myFxr => FmInfo%FXR(ifxr,iz)
                nFsrInFxr = myFxr%nFsrInFxr
                !FsrIdxSt = myFxr%FsrIdxSt
                nisoInFxr = myFxr%niso
                DO k = 1, nFsrInFxr
                    !ifsr = FsrIdxSt + k -1
                    !localFsr=localFsr+1
                    !--- EDIT 14/12/23 Sq. cell problem
                    ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                    localFsr=core%cellinfo(icel)%mapfxr2fsridx(k,j)
                    FsrVol(ifsr)=Cell(icel)%vol(localFsr)*hz
                    myVol=FsrVol(ifsr)
                    DO ig = 1, ng
                        FsrPhi(ifsr,ig)=FmInfo%phis(ifsr,iz,ig)
                        FsrPhiVol(ifsr,ig)=FmInfo%phis(ifsr,iz,ig)*myVol
                        !---BYS edit 15/02/23 Ryumin feedback
                        !FsrPhiVol(ifsr,ig)=fsrphir(ifsr,ig)*myVol
                        !FsrPhi(ifsr,ig)=fsrphir(ifsr,ig)
                        PhiVol(ig)=PhiVol(ig)+FsrPhiVol(ifsr,ig)
                    ENDDO
                    FsrVolsum=FsrVolsum+myVol
                    IF( nTracerCntl%lXsLib )THEN
                        DO iso = 1, nisoInFxr
                            idiso=myFxr%idiso(iso)
                            IF( idiso .EQ. 8001 )THEN
                                idiso=8016
                            ENDIF
                            id = MapNucl(idiso)  !id = Map(92235) = 37
                            !isodata => ldiso(id)
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
                        ENDDO
                    ELSE

                    ENDIF
                ENDDO
            ENDDO
        ENDDO
    ELSE
      DO iy = yst, yed
        DO ix = xbg, xed
            ixyl= AsyInfo(iasytype)%Pin2DIdx(ix, iy)
            ixy = Asy(iasy)%GlobalPinIdx(ixyl) !global Pin Index
            icel = Pin(ixy)%Cell(iz)
            FxrIdxSt = Pin(ixy)%FxrIdxSt
            FsrIdxSt = Pin(ixy)%FsrIdxSt
            nlocalFxr = Cell(icel)%nFxr
            !localFsr = 0
            DO j = 1, nLocalFxr
                nCoreFxr=nCoreFxr+1
                ifxr = FxrIdxSt + j -1
                myFxr => FmInfo%FXR(ifxr,iz)
                nFsrInFxr = myFxr%nFsrInFxr
                !FsrIdxSt = myFxr%FsrIdxSt
                nisoInFxr = myFxr%niso
                DO k = 1, nFsrInFxr
                    !ifsr = FsrIdxSt + k -1
                    !localFsr=localFsr+1
                    !--- EDIT 14/12/23 Sq. cell problem
                    ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                    localFsr=core%cellinfo(icel)%mapfxr2fsridx(k,j)
                    FsrVol(ifsr)=Cell(icel)%vol(localFsr)*hz
                    myVol=FsrVol(ifsr)
                    DO ig = 1, ng
                        FsrPhi(ifsr,ig)=FmInfo%phis(ifsr,iz,ig)
                        FsrPhiVol(ifsr,ig)=FmInfo%phis(ifsr,iz,ig)*myVol
                        !---BYS edit 15/02/23 Ryumin feedback
                        !FsrPhiVol(ifsr,ig)=fsrphir(ifsr,ig)*myVol
                        !FsrPhi(ifsr,ig)=fsrphir(ifsr,ig)
                        PhiVol(ig)=PhiVol(ig)+FsrPhiVol(ifsr,ig)
                    ENDDO
                    FsrVolsum=FsrVolsum+myVol
                    IF( nTracerCntl%lXsLib )THEN
                        DO iso = 1, nisoInFxr
                            idiso=myFxr%idiso(iso)
                            IF( idiso .EQ. 8001 )THEN
                                idiso=8016
                            ENDIF
                            id = MapNucl(idiso)  !id = Map(92235) = 37
                            !isodata => ldiso(id)
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
                        ENDDO
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
        ENDDO
    ENDIF
    IF( .NOT. nTracerCntl%lXsLib ) nisotot=1
    ALLOCATE(isoMicXs(nisotot,0:nxs,ng), isoMicSm(nisotot,ng,ng),isoMicXs2g(nisotot,0:nxs,ngrp),isoMicSm2g(nisotot,ngrp,ngrp)) !tr, a, r, f, nu, k
    ALLOCATE(isoMacXs(0:nisotot,0:nxs,ng), isoMacSm(0:nisotot,ng,ng),isoMacXs2g(0:nisotot,0:nxs,ngrp),isoMacSm2g(0:nisotot,ngrp,ngrp)) !D, a, r, f, nf, kf
    isoMacXs=0;isoMacSm=0;isoMacXs2g=0;isoMicXs=0;isoMicSm=0;isoMicXs2g=0; isoMicSm2g=0; isoMacSm2g=0;
    !--- P1 variable
    ALLOCATE(isoMacSmP1(0:nisotot,ng,ng),isoMicSmP1(nisotot,ng,ng))
    ALLOCATE(isoMacSm2gP1(0:nisotot,ngrp,ngrp),isoMicSm2gP1(nisotot,ngrp,ngrp))
    isoMacSmP1=0;isoMicSmP1=0;isoMacSm2gP1=0;isoMicSm2gP1=0

    IF( nisotot .GT. isosize )THEN
        write(*,*) '                                EXCEED. BYS edit 14/02/27'
    ENDIF
    RFsrVolSum=1.0/FsrVolsum
    CALL MULTI_CA(RFsrVolsum, isoNumden(:), nisotot) !isoNumden=isoNumden/FsrVolsum ! isoNumden = avg(n)_iso
    IF( .NOT. nTracerCntl%lXsLib ) isonumden=1.0_8
    h2onumden=rFsrVolsum*h2onumden
    uo2numden=rFsrVolsum*uo2numden


ENDSUBROUTINE


SUBROUTINE GCRst(Core, FmInfo, nTracerCntl, ng)
    !GCIso : Search the list of isotope in FSR
    USE TYPEDEF,         ONLY : CoreInfo_Type, FmInfo_Type,  Mixture_Type
    USE PARAM, ONLY : ckelvin
    USE GC_mod
    USE GCpin_mod
    USE Depl_Mod,  ONLY : DeplCntl  !15/03/10 r534 : added for GC gen while depl_calc
    USE Material_Mod,  ONLY : Mixture
    use ioutil, ONLY : fndchara512
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FMInfo_Type) :: FmInfo
    TYPE(nTracerCntl_Type) :: nTracerCntl
    INTEGER :: ng

    INTEGER :: ig
    INTEGER :: j, k, i

    INTEGER :: iso, isoidx, id, idiso, nisoinFxr, FxrIdxSt
    INTEGER :: ixy, ixyl, ix, iy ! local ixy
    INTEGER :: nLocalFxr, iCel, iFxr, icell
    REAL :: pnum
    CHARACTER(256) :: fn, fn2, fn3
    CHARACTER(512) :: dataline
    CHARACTER(1024) :: gstr
    
    INTEGER :: io, io2
    INTEGER :: matidx, matidxst, imix
    INTEGER :: ncircle
    INTEGER :: imixlist(10000)
    
    INTEGER :: iRegFr,iRegTo, ipos(100), nspt, ireg(0:300)=0
    
    
    
    io=44
    io2=io+1
    
    IF( .NOT. nTracerCntl%lXsLib ) RETURN
    
    WRITE(fn,'(A)') TRIM(caseid)
    IF( .NOT. lSA )THEN
        WRITE(fn,'(A,A)') TRIM(caseid),'_'
        IF( iasy .LT. 10 )THEN
            WRITE(fn,'(A,A,i1)') TRIM(fn),'0', iasy
        ELSE
            WRITE(fn,'(A,i2)') TRIM(fn), iasy
        ENDIF
    ENDIF
    IF( lZdir )THEN
        WRITE(fn,'(A,A)') TRIM(fn),'_z'
        IF( iz .LT. 10 )THEN
            WRITE(fn,'(A,A,i1)') TRIM(fn),'0', iz
        ELSE
            WRITE(fn,'(A,i2)') TRIM(fn), iz
        ENDIF
    ENDIF
    IF( DeplCntl%nburnupstep .GT. 1 )THEN
        WRITE(fn,'(A,A)') TRIM(fn),'_d'
        IF( DeplCntl%nowstep .LT. 10 )THEN
            WRITE(fn,'(A,A,i1)') TRIM(fn),'00',DeplCntl%nowstep
        ELSEIF( DeplCntl%nowstep .LT. 100 )THEN
            WRITE(fn,'(A,A,i2)') TRIM(fn),'0',DeplCntl%nowstep
        ELSE
            WRITE(fn,'(A,i3)') TRIM(fn),DeplCntl%nowstep
        ENDIF
    ENDIF        
    WRITE(fn2,'(A,A)') TRIM(fn),'.mat'
    OPEN(unit=io, file = fn2, status = 'replace')
    IF(  DeplCntl%nburnupstep .LE. 1 .AND. DeplCntl%nowstep.EQ.0 )THEN
        WRITE(fn3,'(A,A)') TRIM(fn),'.geo'
        OPEN(unit=io2, file = fn3, status = 'replace')
    ENDIF
    
    
    matidx=0;
    icell=0;
    DO ixyl = 1, nxy
        ixy = Asy(iasy)%GlobalPinIdx(ixyl) !global Pin Index
        icel = Pin(ixy)%Cell(iz)
        FxrIdxSt = Pin(ixy)%FxrIdxSt
        nlocalFxr = Cell(icel)%nFxr
        IF( cell(icel)%lgap ) CYCLE
        icell=icell+1
        
        !localFsr = 0
        matidxst=matidx+1
        DO j = nLocalFxr, 1, -1
            nCoreFxr=nCoreFxr+1
            ifxr = FxrIdxSt + j -1
            myFxr => FmInfo%FXR(ifxr,iz)
            nisoInFxr = myFxr%niso
            nisoInFxr = max(myFxr%niso, myFxr%niso_depl)
                
            imix=myfxr%imix
            matidx=matidx+1
            imixlist(imix)=matidx
            
            iso=1;
            idiso=myFxr%idiso(iso)
            pnum=myfxr%pnum(iso)
            !id = MapNucl(idiso)  !id = MapNucl(92235) = 37
            WRITE(io,'(a,i7,a1,a3,i3,f7.3,f8.2,a3,i6,1p,e15.5)')  '  mixture ' , matidx, ' ', mixture(imix)%name, mixture(imix)%deplopt, mixture(imix)%dens, mixture(imix)%temp-ckelvin, ' / ', idiso, pnum
            !WRITE(io,'()')
            DO iso = 2, nisoInFxr
                idiso=myFxr%idiso(iso)
                pnum=myfxr%pnum(iso)
                WRITE(io,'(a,a7,a1,a6,a7,a,a3,i6,1p,e15.5)')  '          ' ,'       ', ' ', '      ', '       ', '        ', '   ', idiso, pnum
            ENDDO
        ENDDO
        
        IF( DeplCntl%nowstep .EQ. 0 )THEN
          WRITE(gstr, '(a,i5,i5)') ' cell ', icell, Cell(icel)%basecellstr
          IF( Cell(icel)%geom%lcircle )THEN
              ncircle=Cell(icel)%geom%ncircle
              DO j = 1, ncircle
                  WRITE(gstr, '(a,f8.5)'), TRIM(gstr), Cell(icel)%geom%circle(3,ncircle+1-j)
              ENDDO
              WRITE(gstr, '(a,a)'), TRIM(gstr), ' / '
              DO j = matidxst, matidx
                  WRITE(gstr, '(a,i5)'), TRIM(gstr), j
              ENDDO
              WRITE(gstr, '(a,a)'), TRIM(gstr), ' / '
              DO j = matidxst, matidx-1
                  WRITE(gstr, '(a,i2)'), TRIM(gstr), 1
              ENDDO
              IF( Cell(icel)%geom%lccent )THEN
                  SELECT CASE(Cell(icel)%geom%CCentType)
                  CASE(1)
                      WRITE(gstr, '(a,a)'), TRIM(gstr), ' / SW'
                  CASE(2)
                      WRITE(gstr, '(a,a)'), TRIM(gstr), ' / NW'
                  CASE(3)
                      WRITE(gstr, '(a,a)'), TRIM(gstr), ' / NE'
                  CASE(4)
                      WRITE(gstr, '(a,a)'), TRIM(gstr), ' / SE'
                  ENDSELECT
              ENDIF
          ELSEIF( Cell(icel)%geom%lrect )THEN
              IF( matidxst .EQ. matidx )THEN
                  WRITE(gstr, '(a,a,i5,a,i3,a,i3)'), TRIM(gstr), ' / / ', matidxst, ' / ', Cell(icel)%geom%nx, ' / ', Cell(icel)%geom%ny
              ELSE
                  IF( Cell(icel)%geom%inpnx .EQ. 1 )THEN
                      WRITE(gstr, '(a,a)'), TRIM(gstr), '         / '
                  ELSE
                      DO j = 1, Cell(icel)%geom%inpnx-1
                          WRITE(gstr, '(a,f8.5)'), TRIM(gstr), Cell(icel)%geom%inpdelx(j)
                      ENDDO
                      WRITE(gstr, '(a,a)'), TRIM(gstr), ' / '
                  ENDIF
                  IF( Cell(icel)%geom%inpny .EQ. 1 )THEN
                      WRITE(gstr, '(a,a)'), TRIM(gstr), '         / '
                  ELSE
                      DO j = 1, Cell(icel)%geom%inpny-1
                          WRITE(gstr, '(a,f8.5)'), TRIM(gstr), Cell(icel)%geom%inpdely(j)
                      ENDDO
                      WRITE(gstr, '(a,a)'), TRIM(gstr), ' / '
                  ENDIF
                  
                  dataline=Cell(icel)%CellInput
                  CALL fndchara512(dataline,ipos,nspt,SLASH)
                  iRegFr = 1
                  iRegTo = Cell(icel)%geom%inpnx*Cell(icel)%geom%inpny
                  read(dataline(ipos(2)+1:ipos(3)-1),*) (ireg(i),i=iRegFr,iRegTo)
                  !DO j = matidxst, matidx
                  !    WRITE(gstr, '(a,i5)'), TRIM(gstr), j
                  !ENDDO          
                  DO j = iRegFr, iRegTo
                      WRITE(gstr, '(a,i5)'), TRIM(gstr), imixlist(ireg(j))
                  ENDDO                
                  WRITE(gstr, '(a,a)'), TRIM(gstr), ' / '
                  
                  DO j = 1, Cell(icel)%geom%inpnx
                      WRITE(gstr, '(a,i3)'), TRIM(gstr), Cell(icel)%geom%inpdivx(j)
                  ENDDO
                  WRITE(gstr, '(a,a)'), TRIM(gstr), ' / '
                  DO j = 1, Cell(icel)%geom%inpny
                      WRITE(gstr, '(a,i3)'), TRIM(gstr), Cell(icel)%geom%inpdivy(j)
                  ENDDO
                  !WRITE(*,*) 'GEOM RESTART UNEXPECTED GEOM at # ', ixyl           
              ENDIF                
          ENDIF
          WRITE(io2,'(a)'), TRIM(gstr)
        ENDIF
    ENDDO

    CLOSE(io)
    IF(  DeplCntl%nburnupstep .LE. 1 .AND. DeplCntl%nowstep.EQ.0 )THEN
        CLOSE(io2)
    ENDIF
ENDSUBROUTINE



SUBROUTINE GCSpec(Core, FmInfo, GroupInfo, nTracerCntl, PE, ng, ispec)
    USE TYPEDEF,        ONLY : CoreInfo_Type, FmInfo_Type, GroupInfo_Type, PE_Type
    USE CNTL,           ONLY : nTracerCntl_Type
    USE BasicOperation, ONLY : MULTI_CA
    USE GC_mod
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FMInfo_Type) :: FmInfo
    TYPE(GroupInfo_Type) :: GroupInfo
    TYPE(nTracerCntl_Type) :: nTracerCntl
    TYPE(PE_TYPE) :: PE
    INTEGER :: ng
    INTEGER :: ispec

    INTEGER :: ifsr, ig
    REAL :: critphisum, infphisum, Rcritphisum, Rinfphisum


    SELECT CASE(ispec)
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

ENDSUBROUTINE

SUBROUTINE GCPwrNorm(Core, FmInfo, THInfo, CmInfo, nTracerCntl, GroupInfo, PE, ng)
    USE TYPEDEF,      ONLY : CoreInfo_Type, FmInfo_Type, THInfo_Type, CMInfo_Type, PE_Type, GroupInfo_Type
    USE CNTL,         ONLY : nTracerCntl_Type
    USE GC_mod
    USE GCpin_mod
    USE XsUtil_mod, ONLY : FreeXsIsoMac, FreeXsMac
    USE MacXsLib_Mod, ONLY : MacXsNf_FXR
    USE BenchXs,        ONLY : MacXsBen
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FMInfo_Type) :: FmInfo
    TYPE(THInfo_Type) :: THInfo
    TYPE(CMInfo_Type) :: CMInfo
    TYPE(nTracerCntl_Type) :: nTracerCntl
    TYPE(GroupInfo_Type) :: GroupInfo
    TYPE(PE_TYPE) :: PE
    TYPE(XsMac_Type) :: XsMac
    INTEGER :: ng

    INTEGER :: i, j, k
    INTEGER :: ig, ig2

    INTEGER :: ixy, ixyl
    INTEGER :: iCel, iso, iFxr, iFsr, LocalFSR, FxrIdxSt, FsrIdxSt, nLocalFxr, nFsrInFxr, nIsoinFXR
    INTEGER :: imix
    REAL :: myVol
    REAL :: isoFSV, FSV, CFSV ! Fission Source*Volume
    REAL :: MacXsSm(GroupInfo%ng,GroupInfo%ng)
    REAL :: MacXsSmP1(GroupInfo%ng,GroupInfo%ng)

    CALL CorePowerCal(Core, CmInfo, PowerDist, ng, nTracerCntl, PE)

    CFSV=0;NormFac=0;
    DO iz = PE%myzb, PE%myze
    IF( nTracerCntl%lXsLib ) tempref = THInfo%RefFuelTemp(iz) + CKELVIN
    DO iasy = 1, nxya
        FSV=0; !fission source *Volume
        iasytype=Asy(iasy)%AsyType
        nxy=AsyInfo(iasytype)%nxy
        DO ixyl = 1, nxy
            ixy = Asy(iasy)%GlobalPinIdx(ixyl)
            icel = Pin(ixy)%Cell(iz)
            FxrIdxSt = Pin(ixy)%FxrIdxSt
            FsrIdxSt = Pin(ixy)%FsrIdxSt
            nlocalFxr = Cell(icel)%nFxr
            !localFsr = 0
            DO j = 1, nLocalFxr
                ifxr = FxrIdxSt + j -1
                myFxr => FmInfo%FXR(ifxr,iz)
                nFsrInFxr = myFxr%nFsrInFxr
                !IF( 0 )THEN
                !    !FsrIdxSt = myFxr%FsrIdxSt
                !    nisoInFxr = myFxr%niso
                !    IF( nTracerCntl%lXsLib ) CALL GenFxrIsoMacXS(XsMac, IsoMacXsSm, IsoMacXsSmP1, myFxr, Tempref, 1, ng, GroupInfo, nTRACERCntl, PE)
                !    DO k = 1, nFsrInFxr
                !        !ifsr = FsrIdxSt + k -1
                !        !localFsr=localFsr+1
                !        !--- EDIT 14/12/23 Sq. cell problem
                !        ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                !        localFsr=core%cellinfo(icel)%mapfxr2fsridx(k,j)
                !        isoFSV=0
                !        myVol=Cell(icel)%vol(localFsr)*Core%HZ(iz)
                !        DO ig = 1, ng
                !            FsrPhiVol(ifsr,ig)=FmInfo%phis(ifsr,iz,ig)*myVol
                !        ENDDO
                !        DO iso = 1, nisoInFxr
                !            DO ig2 = 1, ng
                !                isoFSV=isoFSV+XsMac%isoXsMacKF(iso,ig2) *FsrPhiVol(ifsr,ig2)
                !            ENDDO
                !        ENDDO !---END of Iso sweep
                !        FSV=FSV+isoFSV
                !    ENDDO !---END of Fsr sweep
                !    CALL FreeXsIsoMac(XsMac)
                !    !DEALLOCATE(IsoMacXsSm)
                !ELSE
                    IF( nTracerCntl%lXsLib )THEN
                        CALL MacXsNF_FXR(XsMac, myFxr, 1, ng, ng, 1.0, .FALSE.)
                        DO k = 1, nFsrInFxr
                            !--- EDIT 14/12/23 Sq. cell problem
                            ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                            localFsr=core%cellinfo(icel)%mapfxr2fsridx(k,j)
                            isoFSV=0
                            myVol=Cell(icel)%vol(localFsr)*Core%HZ(iz)
                            DO ig = 1, ng
                                FsrPhiVol(ifsr,ig)=FmInfo%phis(ifsr,iz,ig)*myVol
                                isoFSV=isoFSV+XsMac%XsMacNf(ig)*FsrPhiVol(ifsr,ig)
                            ENDDO
                            FSV=FSV+isoFSV
                        ENDDO !---END of Fsr sweep
                        CALL FreeXsMac(XsMac)
                    ELSE
                        imix=myFXR%imix
                        DO k = 1, nFsrInFxr
                            !--- EDIT 14/12/23 Sq. cell problem
                            ifsr = FsrIdxSt-1+core%cellinfo(icel)%mapfxr2fsridx(k,j)
                            localFsr=core%cellinfo(icel)%mapfxr2fsridx(k,j)
                            isoFSV=0
                            myVol=Cell(icel)%vol(localFsr)*Core%HZ(iz)
                            DO ig = 1, ng
                                FsrPhiVol(ifsr,ig)=FmInfo%phis(ifsr,iz,ig)*myVol
                                isoFSV=isoFSV+MacXsBen(imix)%XsNf(ig)*FsrPhiVol(ifsr,ig)
                            ENDDO
                            FSV=FSV+isoFSV
                        ENDDO !---END of Fsr sweep
                    ENDIF
                !ENDIF
            ENDDO !---END of Fxr sweep
        ENDDO !---END of nPin sweep
        CFSV=CFSV+FSV
    ENDDO !---END of Assembly sweep
    ENDDO !---END of Plan
    NormFac=Core%PowerCore
    NormFac=NormFac/CFSV
    IF( Core%PowerCore .LE. 0._8 ) NormFac=1.0_8


ENDSUBROUTINE

SUBROUTINE SetPinList(lCentX,lCentY,lCentXY)
    USE GC_mod
    IMPLICIT NONE
    INTEGER :: i, ix, iy
    LOGICAL :: lCentX, lCentY, lCentXY
    i=ipin
    IF( .NOT. lGapHom )THEN
        ix=mod(i-1,nx)+1;
        iy=(i-ix)/nx+1;
        xbg=ix;        xed=xbg;
        yst=iy;        yed=yst;
    ELSE
        IF( lCentY .OR. lCentXY )THEN
            ix=mod(i-1,nx-1)+1;
            iy=(i-ix)/(nx-1)+1;
        ELSE
            ix=mod(i-1,nx-2)+1;
            iy=(i-ix)/(nx-2)+1;
        ENDIF
        
        !IF( lCentY .OR. lCentXY )THEN
        !    ix=mod(i-1,nx-1)+1;
        !ELSE
        !    ix=mod(i-1,nx-2)+1;
        !ENDIF
        !IF( lCentX .OR. lCentXY )THEN
        !    iy=(i-ix)/(nx-1)+1;
        !ELSE
        !    iy=(i-ix)/(nx-2)+1;
        !ENDIF
        
        IF( lCentY .OR. lCentXY )THEN
            IF( ix .EQ. nx-1 )THEN
                xbg=nx-1;
                xed=nx;
            ELSE
                xbg=ix;
                xed=ix;
            ENDIF       
        ELSE
            IF( ix .EQ. 1 )THEN
                xbg=1;
                xed=2; 
            ELSEIF( ix .EQ. nx-2 )THEN
                xbg=nx-1;
                xed=nx;
            ELSE
                xbg=ix+1;
                xed=xbg;
            ENDIF            
        ENDIF
        
        IF( lCentX .OR. lCentXY )THEN
            IF( iy .EQ. ny-1 )THEN
                yst=ny-1;
                yed=ny;
            ELSE
                yst=iy;
                yed=iy;
            ENDIF       
        ELSE
            IF( iy .EQ. 1 )THEN
                yst=1;
                yed=2; 
            ELSEIF( iy .EQ. ny-2 )THEN
                yst=ny-1;
                yed=ny;
            ELSE
                yst=iy+1;
                yed=yst;
            ENDIF            
        ENDIF
                
    ENDIF
ENDSUBROUTINE

SUBROUTINE GCFin
    USE GC_mod
    USE GCpin_mod
    IMPLICIT NONE
    !--- ALLOCATE ---
    DEALLOCATE(Phi, PhiVol, PinPhiVol)
    DEALLOCATE(GrpIdx, PhiG)
    DEALLOCATE(FsrPhi, FsrPhiVol, FsrVol)


    !For B1
    DEALLOCATE(f_ratio, phicrit, Dng, XstrB1, Xstr, Xst)
    DEALLOCATE(ADF,ADF2g, avgADF2g, ADF2GSurf)
    DEALLOCATE(SurfJ,SurfJG, SurfPhiG, ChiG)
    DEALLOCATE(avgBPhiG, conPhiG)

    DEALLOCATE(pinRmv,pinFis,pinSS,pinR,pinFlux,pinJ)
    DEALLOCATE(pinSin,pinSout,pinAbs)
    DEALLOCATE(pinSinG,pinSoutG)
    DEALLOCATE(pinJdir)

ENDSUBROUTINE



