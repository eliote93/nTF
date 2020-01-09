#include <defines.h>

SUBROUTINE ProcessEffMGXs(Core, FmInfo, CmInfo, ThInfo, GroupInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF
USE CNTL,             ONLY : nTracerCntl_Type,    OutpCntl_Type
USE FILES,            ONLY : io16,                caseid
USE XsUtil_mod,       ONLY : GetXsMacDat,         ReturnXsMacDat
USE MacXsLib_mod,     ONLY : EffMacXS,   EffRIFPin,   MacXsBase,   IsoMacXsScatMatrix_Gen, IsoMacXsScatMatrix_GenR
USE MOC_Mod,          ONLY : FxrAvgPhi
USE ioutil
USE BasicOperation
USE TH_Mod,           ONLY : GetPinFuelTemp
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

TYPE(Pin_Type), POINTER :: PIN(:)
TYPE(ResVarPin_Type), POINTER :: myPin
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :), myFxr
TYPE(XsMac_Type), POINTER :: XsMac

INTEGER :: iz, ixa, iya, ixc, iyc, ir
INTEGER :: ifxr, ixy, ixya, ixyc, iasytype, icel
INTEGER :: i, j, k, id, ig, iresoGrp1, iresoGrp2, ig2, icg, icg2
!Fxr Data
INTEGER :: niso, ng, iso
LOGICAL :: lres, lAIC
INTEGER, POINTER :: idiso(:), maxniso(:)
REAL, POINTER :: pnum(:), NUMDEN(:,:)
REAL :: eigv
REAL :: Tempref, Temp

INTEGER,PARAMETER :: MAXNUM = 500

!INTEGER,parameter :: GB2G=26

!Output Data
REAL, POINTER :: outMGxstot(:),outMGxsnftot(:),outMGscatmattot(:,:),fluxMGtot(:),chiMGtot(:)
REAL :: areatot,psitot
REAL, POINTER :: chiMG(:,:),psi(:),IsoXsMacAold(:,:),IsoXsMacNFold(:,:),IsoXsMacFold(:,:),IsoXsMacSold(:,:),IsoXsMacSSold(:,:),IsoXsMacS1old(:,:)
REAL, POINTER :: OutMGXS(:, :, :),OutMGXSNF(:, :, :),OUTMGSCATMAT(:,:,:,:),OUTMGCHI(:,:,:)
REAL, POINTER :: OutMGFlux(:, :), fxrisopsi(:), fxrisochi(:,:), LocalPsi(:,:)
REAL, POINTER :: Area(:)
INTEGER, POINTER :: IDMAP(:,:)  ! idmap(1,1) : 1st isotope id(c.g.92238) in the 1st tally region
INTEGER :: GBdry(0:ngmax), ncmg
REAL :: Flux(ngmax), PinFuelTempAvgsq
CHARACTER(256) :: fn
#ifdef __GFORTRAN__
INTEGER :: ncmgp3
CHARACTER(2) :: str2num
#endif

IF(.NOT. nTracerCntl%lXsLib) RETURN
IF(PE%nCmfdProc .GT. 1) RETURN

Pin => Core%Pin; CellInfo => Core%CellInfo; Fxr => FmInfo%Fxr
ng = GroupInfo%ng
iresoGrp1 = GroupInfo%nofg + 1; iresoGrp2 = GroupInfo%nofg + GroupInfo%norg
eigv = 1.0
ALLOCATE(OutMGXS(MAXNUM, ng, nTracerCntl%OutpCntl%nRegMGXsOut))
ALLOCATE(OutMGXSNF(MAXNUM, ng, nTracerCntl%OutpCntl%nRegMGXsOut))
ALLOCATE(OutMGSCATMAT(MAXNUM, ng,ng, nTracerCntl%OutpCntl%nRegMGXsOut))
ALLOCATE(OutMGCHI(MAXNUM, ng, nTracerCntl%OutpCntl%nRegMGXsOut))
ALLOCATE(OutMGFlux(ng, nTracerCntl%OutpCntl%nRegMGXsOut))
ALLOCATE(Area(nTracerCntl%OutpCntl%nRegMGXsOut))
ALLOCATE(IDMAP(MAXNUM,nTracerCntl%OutpCntl%nRegMGXsOut))
ALLOCATE(maxniso(nTracerCntl%OutpCntl%nRegMGXsOut))
ALLOCATE(fxrisopsi(MAXNUM))
ALLOCATE(fxrisochi(MAXNUM, ng))
ALLOCATE(LocalPsi(MAXNUM,nTracerCntl%OutpCntl%nRegMGXsOut))
ALLOCATE(chiMG(ng,nTracerCntl%OutpCntl%nRegMGXsOut))
ALLOCATE(psi(nTracerCntl%OutpCntl%nRegMGXsOut))
ALLOCATE(outMGxstot(ng),outMGxsnftot(ng),outMGscatmattot(ng,ng),fluxMGtot(ng),chiMGtot(ng))
ALLOCATE(NUMDEN(MAXNUM,nTracerCntl%OutpCntl%nRegMGXsOut))

CALL GetXsMacDat(XsMac, ng, TRUE)
DO i = 1, nTracerCntl%OutpCntl%nRegMGXsOut

  ncmg=nTracerCntl%OutpCntl%nMG(i)

  iz = nTracerCntl%OutpCntl%RegMGXsOutList(1, i)
  ixa = nTracerCntl%OutpCntl%RegMGXsOutList(2, i); iya = nTracerCntl%OutpCntl%RegMGXsOutList(3, i)
  ixc = nTracerCntl%OutpCntl%RegMGXsOutList(4, i); iyc = nTracerCntl%OutpCntl%RegMGXsOutList(5, i)

  ixya = Core%CoreIdx(ixa, iya); iasytype = Core%CoreMap(ixya)
  ixyc = Core%AsyInfo(iasyType)%pin2DIdx(ixc, iyc)
  ixy = Core%Asy(ixya)%GlobalPinIdx(ixyc); icel = Pin(ixy)%Cell(iz)
  !mypin => ResVarPin(ixy,iz); lAIC = CellInfo(icel)%lAIC
  !PinFuelTempAvgsq = dsqrt(GetPinFuelTemp(Core, Fxr, iz, ixy))
  !IF (mypin%lres.and.nTracerCntl%lRIF) CALL EffRIFPin(mypin, PinFuelTempAvgsq, iz, lAIC, PE)  
  CALL CP_CA(OutMGXs(1:MAXNUM, 1:ncmg, i), 0._8, MAXNUM, ncmg)
  CALL CP_CA(OutMGXsNF(1:MAXNUM, 1:ncmg, i), 0._8, MAXNUM, ncmg)
  CALL CP_CA(OutMGCHI(1:MAXNUM, 1:ncmg, i), 0._8, MAXNUM, ncmg)
  CALL CP_CA(OutMGSCATMAT(1:MAXNUM, 1:ncmg, 1:ncmg, i), 0._8, MAXNUM, ncmg, ncmg)
  CALL CP_CA(OutMGFLux(1:ncmg, i), 0._8, ncmg)
  CALL CP_CA(fxrisopsi(1:MAXNUM), 0._8, MAXNUM)
  CALL CP_CA(fxrisochi(1:MAXNUM, 1:ncmg), 0._8, MAXNUM, ncmg)
  CALL CP_CA(LocalPsi(1:MAXNUM, i), 0._8, MAXNUM)
  CALL CP_CA(NUMDEN(1:MAXNUM, i), 0._8, MAXNUM)


  GBdry=0
  GBdry(0:ncmg)=nTracerCntl%OutpCntl%MGBdry(0:ncmg, i)

  maxniso(i)=0
  DO ir= nTracerCntl%OutpCntl%RegMGXsOutList(6, i), nTracerCntl%OutpCntl%RegMGXsOutList(7, i)
    ifxr = Pin(ixy)%FxrIdxst + ir - 1;myFxr => Fxr(ifxr, iz)
    tempref = THInfo%RefFuelTemp(iz) + CKELVIN
    temp = myFxr%temp
    niso = myFxr%niso; idiso => myFxr%idiso
    pnum => myFxr%pnum
    DO j=1,niso
        DO k=1,MAXNUM
            if (IDMAP(k,i).eq.0) EXIT
            if (IDMAP(k,i).eq.idiso(j)) goto 10
        ENDDO
        if (k.gt.maxniso(i)) then
            maxniso(i)=maxniso(i)+1
            IDMAP(maxniso(i),i)=idiso(j)
        endif
     10 continue
    ENDDO
  ENDDO
  AREA(i) = 0
  DO ir= nTracerCntl%OutpCntl%RegMGXsOutList(6, i), nTracerCntl%OutpCntl%RegMGXsOutList(7, i)
    ifxr = Pin(ixy)%FxrIdxst + ir - 1;myFxr => Fxr(ifxr, iz)
    temp = myFxr%temp
    niso = myFxr%niso; idiso => myFxr%idiso
    pnum => myFxr%pnum; lres = myFxr%lres
    CALL MacXsBase(XsMac, myFxr, 1, ng, ng, eigv, FALSE, TRUE)
    ALLOCATE(XsMac%IsoXsMacSm(niso,ng,ng))
    ALLOCATE(IsoXsMacAold(niso,ng),IsoXsMacNFold(niso,ng),IsoXsMacFold(niso,ng),IsoXsMacSold(niso,ng),IsoXsMacSSold(niso,ng),IsoXsMacS1old(niso,ng))
    IsoXsMacAold=XsMac%IsoXsMacA
    IsoXsMacNFold=XsMac%IsoXsMacNF
    IsoXsMacFold=XsMac%IsoXsMacF
    IsoXsMacSold=XsMac%IsoXsMacS0
    IsoXsMacS1old=XsMac%IsoXsMacS1
    IsoXsMacSSold=XsMac%IsoXsMacSS
    AREA(i) = Area(i) + myFxr%area

    IF (lres) THEN

      DO ig = iResoGrp1, iResoGrp2

        !CALL EffMacXs(XsMac, myPin, myFxr, PinFuelTempAvgsq, niso, ig, ng, .TRUE., iz, ixy, ir, PE)

        DO k = 1, niso
          XsMac%IsoXsMacA(k, ig) = IsoXsMacAold(k,ig) * myFxr%fresoAIso(k,ig)
          IF (IsoXsMacFold(k,ig).ne.0) XsMac%IsoXsMacNF(k, ig) = IsoXsMacNFold(k,ig) * myFxr%fresoFIso(k,ig)  
        ENDDO
      ENDDO

      if (nTracerCntl%lRST) CALL IsoMacXsScatMatrix_GenR(XsMac%IsoXsMacSm, temp, niso, idiso, pnum, 1, ng, ng, GroupInfo, &
                                nTracerCntl%lscat1, lres, myFxr%fresoSIso, myFxr%fresoSSIso, myFxr%fresoS1Iso)
    else
      CALL IsoMacXsScatMatrix_Gen(XsMac%IsoXsMacSm, temp, niso, idiso, pnum, 1, ng, ng, GroupInfo, nTracerCntl%lscat1)
    ENDIF
    Flux(1:ng) = FxrAvgPhi(Core, Fxr, FMInfo%Phis, ixy, ir, iz, ng, PE) !*nTracerCNTL%PowerCoreNormalizer
    CALL IsoFxrPsiMGChi(XsMac%IsoXsMacNF, idiso, niso, Flux, ng, fxrisopsi(1:niso), fxrisochi(1:niso,1:ncmg), GBdry(0:ncmg),ncmg)

    DO icg=1,ncmg
        OutMGflux(icg, i) = OutMGFlux(icg, i) + myFxr%area * sum(Flux(GBdry(icg-1)+1:GBdry(icg)))
    ENDDO

    DO j = 1, niso
        DO k = 1, maxniso(i)
            if (IDMAP(k,i).eq.idiso(j)) EXIT
        ENDDO
        DO icg=1,ncmg
            DO ig=GBdry(icg-1)+1,GBdry(icg)
                OutMGXs(k, icg, i) = OutMGXs(k, icg, i)+ XsMac%IsoXsMacA(j, ig) * Flux(ig) * myFxr%area
                OutMGXsNF(k, icg, i) = OutMGXsNF(k, icg, i)+ XsMac%IsoXsMacNF(j, ig) * Flux(ig) * myFxr%area
                DO icg2=1,ncmg
                    DO ig2=GBdry(icg2-1)+1,GBdry(icg2)
                        OUTMGSCATMAT(k,icg,icg2,i) = OUTMGSCATMAT(k,icg,icg2,i) + XsMac%IsoXsMacSm(j,ig,ig2) * Flux(ig) * myFxr%area
                    ENDDO
                ENDDO
            ENDDO
            OutMGCHI(k,icg,i)=OutMGCHI(k,icg,i)+fxrisopsi(j)*fxrisochi(j,icg)* myFxr%area
        ENDDO
        LocalPsi(k,i) = LocalPsi(k,i) + fxrisopsi(j) * myFxr%area
        NUMDEN(k,i) = NUMDEN(k,i) + pnum(j) * myFxr%area
    ENDDO

  ENDDO !ir sweep

  DO icg=1,ncmg
      OutMGXs(1:maxniso(i), icg, i) = OutMGXs(1:maxniso(i), icg, i) / OutMGflux(icg, i)
      OutMGXsNF(1:maxniso(i), icg, i) = OutMGXsNF(1:maxniso(i), icg, i) / OutMGflux(icg, i)
      DO icg2=1,ncmg
          OUTMGSCATMAT(1:maxniso(i),icg,icg2,i) = OUTMGSCATMAT(1:maxniso(i),icg,icg2,i) / OutMGflux(icg, i)
      ENDDO
  ENDDO

  DO iso=1,maxniso(i)
    if (LocalPsi(iso,i).ne.0._8) then
        OutMGCHI(iso,1:ncmg,i)=OutMGCHI(iso,1:ncmg,i) / LocalPsi(iso,i)
        LocalPsi(iso,i)=LocalPsi(iso,i)/Area(i)
    endif
    NUMDEN(iso,i)=NUMDEN(iso,i)/Area(i)
  ENDDO
  OutMGflux(1:ncmg, i) = OutMGflux(1:ncmg, i) / Area(i)


  NULLIFY(IDISO, pnum, myFxr)
ENDDO
CALL ReturnXsMacDat(XsMac)
fn=trim(caseid)//'.fxrmg'

CALL openfile(io16,FALSE,FALSE,FALSE, fn)
areatot=0; fluxMgtot=0; outMgxstot=0; outMgxsnftot=0
outMgscatmattot=0
chiMg=0; psi=0
chiMgtot=0; psitot=0
DO i = 1, nTracerCntl%OutpCntl%nRegMGXsOut
  ncmg=nTracerCntl%OutpCntl%nMG(i)
  WRITE(io16, '(A9, I5, A3, 5(I5, I5, A3))') 'Region :',  nTracerCntl%OutpCntl%RegMGXsOutList(1, i), '/', &
                                             nTracerCntl%OutpCntl%RegMGXsOutList(2:3, i), '/', &
                                             nTracerCntl%OutpCntl%RegMGXsOutList(4:5, i), '/',&
                                             nTracerCntl%OutpCntl%RegMGXsOutList(6:7, i), '/'
  WRITE(io16, '(1X,A7,ES16.6)') 'AREA : ',Area(i)
  WRITE(io16, '(1X,A4)') 'Flux'
#ifndef __GFORTRAN__
  WRITE(io16, '(<ncmg>ES16.6)')  OutMGflux(1:ncmg, i)
#else
  READ(str2num, *) ncmg
  WRITE(io16, '(' // TRIM(str2num) // 'ES16.6)')  OutMGflux(1:ncmg, i)
#endif
  DO j=1,maxniso(i)
    WRITE(io16, '(I6,ES16.6)') IDMAP(j,i), NUMDEN(j,i)
    DO icg=1,ncmg
#ifndef __GFORTRAN__
        WRITE(io16, '(<3+ncmg>ES16.6)') OutMGXS(j,icg,i),OutMGXsNF(j,icg,i),OutMGCHI(j,icg,i),OUTMGSCATMAT(j,icg,1:ncmg,i)
#else
        ncmgp3 = ncmg + 3
        READ(str2num, *) ncmgp3
        WRITE(io16, '(' // TRIM(str2num) // 'ES16.6)') OutMGXS(j,icg,i),OutMGXsNF(j,icg,i),OutMGCHI(j,icg,i),OUTMGSCATMAT(j,icg,1:ncmg,i)
#endif
    ENDDO
    DO icg=1,ncmg
        chiMg(icg,i)=chiMg(icg,i)+OutMGCHI(j,icg,i)*LocalPsi(j,i)
    ENDDO
    psi(i)=psi(i)+LocalPsi(j,i)
  ENDDO
  if (psi(i).ne.0._8) then
      chiMg(1:ncmg,i)=chiMg(1:ncmg,i)/psi(i)
  endif
  WRITE(io16, '(1X,A5)') 'total'
  DO icg=1,ncmg
#ifndef __GFORTRAN__
    WRITE(io16, '(<3+ncmg>ES16.6)') sum(OutMGXS(1:maxniso(i),icg,i)),sum(OutMGXsNF(1:maxniso(i),icg,i)),chiMg(icg,i),(sum(OUTMGSCATMAT(1:maxniso(i),icg,icg2,i)),icg2=1,ncmg)
#else
    ncmgp3 = ncmg + 3
    READ(str2num, *) ncmgp3
    WRITE(io16, '(' // TRIM(str2num) // 'ES16.6)') sum(OutMGXS(1:maxniso(i),icg,i)),sum(OutMGXsNF(1:maxniso(i),icg,i)),chiMg(icg,i),(sum(OUTMGSCATMAT(1:maxniso(i),icg,icg2,i)),icg2=1,ncmg)
#endif
  ENDDO
  WRITE(io16, '(A)') '--------------------------------------------------------------------------------------------'

  areatot=areatot+area(i)
  DO icg=1,ncmg
      fluxMgtot(icg)=fluxMgtot(icg)+area(i)*OutMGflux(icg, i)
      outMgxstot(icg)=outMgxstot(icg)+sum(OutMGXS(1:maxniso(i),icg,i))*area(i)*OutMGflux(icg, i)
      outMgxsnftot(icg)=outMgxsnftot(icg)+sum(OutMGXSnf(1:maxniso(i),icg,i))*area(i)*OutMGflux(icg, i)
      DO icg2=1,ncmg
          outMgscatmattot(icg,icg2)=outMgscatmattot(icg,icg2)+sum(OUTMGSCATMAT(1:maxniso(i),icg,icg2,i))*area(i)*OutMGflux(icg, i)
      ENDDO
      chiMgtot(icg)=chiMgtot(icg)+chiMg(icg,i)*psi(i)*area(i)
  ENDDO
  psitot=psitot+psi(i)*area(i)
ENDDO
DO icg=1,ncmg
    outMgxstot(icg)=outMgxstot(icg)/fluxMgtot(icg)
    outMgxsnftot(icg)=outMgxsnftot(icg)/fluxMgtot(icg)
    DO icg2=1,ncmg
        outMgscatmattot(icg,icg2)=outMgscatmattot(icg,icg2)/fluxMgtot(icg)
    ENDDO
    fluxMgtot(icg)=fluxMgtot(icg)/areatot
    if (psitot.ne.0) chiMgtot(icg)=chiMgtot(icg)/psitot
ENDDO
WRITE(io16, '(A)') '--------------------------------------------------------------------------------------------'
WRITE(io16,'(1X,A12)') 'REGION TOTAL'
WRITE(io16, '(1X,A7,ES16.6)') 'AREA : ',areatot
WRITE(io16, '(1X,A4)') 'Flux'
#ifndef __GFORTRAN__
WRITE(io16, '(<ncmg>ES16.6)')  (fluxMgtot(icg),icg=1,ncmg)
#else
READ(str2num, *) ncmg
WRITE(io16, '(' // TRIM(str2num) // 'ES16.6)')  (fluxMgtot(icg),icg=1,ncmg)
#endif
WRITE(io16, '(1X,A5)') 'total'
DO icg=1,ncmg
#ifndef __GFORTRAN__
    WRITE(io16,'(<3+ncmg>ES16.6)') outMgxstot(icg),outMgxsnftot(icg),chiMgtot(icg),(outMgscatmattot(icg,icg2),icg2=1,ncmg)
#else
    ncmgp3 = ncmg + 3
    READ(str2num, *) ncmgp3
    WRITE(io16,'(' // TRIM(str2num) // 'ES16.6)') outMgxstot(icg),outMgxsnftot(icg),chiMgtot(icg),(outMgscatmattot(icg,icg2),icg2=1,ncmg)
#endif
ENDDO
close(io16)

DEALLOCATE(OutMGXS,OutMGXSNF,OutMGSCATMAT,OutMGCHI,OutMGFlux,Area,IDMAP)
DEALLOCATE(maxniso,fxrisopsi,fxrisochi,LocalPsi,chiMg,chiMgtot,psi)
DEALLOCATE(outMgxstot,outMgxsnftot,outMgscatmattot,fluxMgtot,NUMDEN)

END SUBROUTINE


SUBROUTINE IsoFxrPsiMGChi(IsoXsMacNF, idiso, niso, phi, ng, localpsi, regisochi, GBdry, ncmg)
USE PARAM, ONLY : epsm10
USE XSLIB_MOD, ONLY : ldiso,mapnucl,nchihel,nelrhel,libdata
IMPLICIT NONE
TYPE(LIBDATA), POINTER :: isodata
INTEGER :: idiso(niso)
INTEGER :: niso, ng, ncmg, GBdry(0:ncmg)
REAL :: phi(ng)

INTEGER :: iso, ig, id, icg, igi,ige
REAL :: localPsi(niso)
REAL :: regisochi(niso,ncmg)

REAL :: IsoXsMacNf(niso, ng)


!Initialize the chi
regisochi=0._8
DO iso = 1, niso
  id = mapnucl(idiso(iso))
  isodata => ldiso(id)
  IF(isodata%ityp .le. 1) CYCLE
  IF(isodata%chi(1) .le. epsm10) CYCLE
  !Obtain a Fission Source
  localpsi(iso) = 0._8
  DO ig = 1, ng
    localpsi(iso) = localpsi(iso) + IsoXsMacNf(iso, ig) * phi(ig)
  ENDDO
  !Calculate a Fission Neutron Dist.
  DO icg=1,ncmg
      igi=GBdry(icg-1)+1
      ige=min(GBdry(icg),nchihel)
      DO ig = igi, ige
        regisochi(iso,icg) = regisochi(iso,icg) + isodata%chi(ig)
      ENDDO
  ENDDO
ENDDO

IF(ASSOCIATED(isodata)) NULLIFY(isodata)
END SUBROUTINE
