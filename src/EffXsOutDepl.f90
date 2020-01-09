#include <defines.h>

SUBROUTINE ProcessEffXsDepl(Core, FmInfo, CmInfo, ThInfo, GroupInfo, nTracerCntl,DeplCntl,PE,iburnup)
USE PARAM
USE TYPEDEF
USE DeplType,         ONLY : DeplCntl_Type
USE CNTL,             ONLY : nTracerCntl_Type,    OutpCntl_Type
USE FILES,            ONLY : io16,                caseid
USE XsUtil_mod,       ONLY : GetXsMacDat,         ReturnXsMacDat
USE MacXsLib_mod,     ONLY : EffMacXS,   EffRIFPin,         MacXsBase
USE TH_Mod,           ONLY : GetPinFuelTemp
USE MOC_Mod,          ONLY : FxrAvgPhi
USE XSLIB_MOD
USE ioutil
USE BasicOperation
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
TYPE(DeplCntl_Type) :: DeplCntl

TYPE(Pin_Type), POINTER :: PIN(:)
TYPE(ResVarPin_Type), POINTER :: myPin
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(FxrInfo_Type), POINTER :: myFxr
TYPE(XsMac_Type), POINTER :: XsMac
TYPE(LIBDATA), POINTER :: isodata

INTEGER :: iz, ixa, iya, ixc, iyc, ir
INTEGER :: ifxr, ixy, ixya, ixyc, iasytype, icel
INTEGER :: i, j, k, id, ig, iresoGrp1, iresoGrp2
INTEGER :: iburnup
CHARACTER*2 :: sburnup
!Fxr Data
INTEGER :: niso, ng, iso
LOGICAL :: lres
INTEGER, POINTER :: idiso(:)
REAL, POINTER :: pnum(:)
REAL :: eigv
REAL :: Temp

!Output Data
REAL, POINTER :: OutXS(:, :, :),OutXSNF(:, :, :),OutXSF(:, :, :)
REAL, POINTER :: OutFlux(:, :),SigFOld(:, :)
REAL, POINTER :: Area(:), Num(:, :)
REAL :: Flux(ngmax), PinFuelTempAvgsq
CHARACTER(256) :: fn
#ifdef __GFORTRAN__
INTEGER :: nisop2
CHARACTER(3) :: c_niso, c_nisop2
#endif
LOGICAL :: lAIC

IF(.NOT. nTracerCntl%lXsLib) RETURN
IF(PE%nCmfdProc .GT. 1) RETURN

Pin => Core%Pin; CellInfo => Core%CellInfo; Fxr => FmInfo%Fxr
ng = GroupInfo%ng
iresoGrp1 = GroupInfo%nofg + 1; iresoGrp2 = GroupInfo%nofg + GroupInfo%norg
eigv = 1.0

ALLOCATE(OutXS(200, ng, nTracerCntl%OutpCntl%nDeplRegXsOut))
ALLOCATE(OutXSNF(200, ng, nTracerCntl%OutpCntl%nDeplRegXsOut))
ALLOCATE(OutXSF(200, ng, nTracerCntl%OutpCntl%nDeplRegXsOut))
ALLOCATE(SigFOld(200, ng))

SigFOld=0._8

write(sburnup,'(i2)') iburnup
fn=trim(caseid)//'_BU'//adjustl(sburnup)//'.exs'

CALL openfile(io16,FALSE,FALSE,FALSE, fn)

CALL GetXsMacDat(XsMac, ng, TRUE)

DO i = 1, nTracerCntl%OutpCntl%nDeplRegXsOut
  iz = nTracerCntl%OutpCntl%DeplRegXsOutList(1, i)
  ixa = nTracerCntl%OutpCntl%DeplRegXsOutList(2, i); iya = nTracerCntl%OutpCntl%DeplRegXsOutList(3, i)
  ixc = nTracerCntl%OutpCntl%DeplRegXsOutList(4, i); iyc = nTracerCntl%OutpCntl%DeplRegXsOutList(5, i)

  !ir sweep
  ixya = Core%CoreIdx(ixa, iya); iasytype = Core%CoreMap(ixya)
  ixyc = Core%AsyInfo(iasyType)%pin2DIdx(ixc, iyc)
  ixy = Core%Asy(ixya)%GlobalPinIdx(ixyc); icel = Pin(ixy)%Cell(iz)
  !myPin => ResVarPin(ixy,iz); lAIC = CellInfo(icel)%lAIC
  
  CALL CP_CA(OutXs(1:nTracerCntl%OutpCntl%DeplIsoXsOutList(0, i), 1:ng, i), 0._8, nTracerCntl%OutpCntl%DeplIsoXsOutList(0, i), ng)
  CALL CP_CA(OutXsF(1:nTracerCntl%OutpCntl%DeplIsoXsOutList(0, i), 1:ng, i), 0._8, nTracerCntl%OutpCntl%DeplIsoXsOutList(0, i), ng)
  CALL CP_CA(OutXsNF(1:nTracerCntl%OutpCntl%DeplIsoXsOutList(0, i), 1:ng, i), 0._8, nTracerCntl%OutpCntl%DeplIsoXsOutList(0, i), ng)

  !PinFuelTempAvgsq = dsqrt(GetPinFuelTemp(Core, Fxr, iz, ixy))
  !IF (myPin%lres.and.nTracerCntl%lRIF) CALL EffRIFPin(myPin, PinFuelTempAvgsq, iz, lAIC, PE)
  
  DO ir= nTracerCntl%OutpCntl%DeplRegXsOutList(6, i), nTracerCntl%OutpCntl%DeplRegXsOutList(7, i)
    ifxr = Pin(ixy)%FxrIdxst + ir - 1; myFxr => Fxr(ifxr, iz)
    temp = myFxr%temp
    niso = myFxr%niso; idiso => myFxr%idiso
    pnum => myFxr%pnum; lres = myFxr%lres

    CALL MacXsBase(XsMac, myFxr, 1, ng, ng, eigv, FALSE, TRUE)

    DO ig = iResoGrp1, iResoGrp2

        !CALL EffMacXs(XsMac, myPin, myFxr, PinFuelTempAvgsq, niso, ig, ng, .TRUE., iz, ixy, ir, PE)
        DO iso=1,niso
          XsMac%IsoXsMacA(iso, ig) =XsMac%IsoXsMacA(iso, ig) *myFxr%fresoAIso(iso,ig)
          XsMac%IsoXsMacNF(iso, ig)=XsMac%IsoXsMacNF(iso, ig)*myFxr%fresoFIso(iso,ig)
          XsMac%IsoXsMacF(iso, ig) =XsMac%IsoXsMacF(iso, ig) *myFxr%fresoFIso(iso,ig)
        ENDDO

    ENDDO

    Flux(1:ng) = FxrAvgPhi(Core, Fxr, FMInfo%Phis, ixy, ir, iz, ng, PE)

    DO j = 1, niso
      id = mapnucl(idiso(j));   isodata => ldiso(id)
      DO ig =  1, ng
        OutXs(j, ig, i) = XsMac%IsoXsMacA(j, ig) / pnum(j)
        OutXsNF(j, ig, i) = XsMac%IsoXsMacNF(j, ig) / pnum(j)
        OutXsF(j, ig, i) = XsMac%IsoXsMacF(j, ig) / pnum(j)
      ENDDO
    ENDDO


    WRITE(io16, '(A10, I5, A3, 5(I5, I5, A3),I5)') 'Region :',  iz, '/', &
                                             ixa,iya, '/', &
                                             ixc,iyc, '/',&
                                             ir,ir, '/', niso
    WRITE(io16, '(A)') 'ABSORPTION'
#ifndef __GFORTRAN__
    WRITE(io16, '(37X, <niso>ES16.6)') (pnum(j),j=1, niso)
    WRITE(io16, '(A5, 2A16, <niso>I16)') 'GRP', 'Flux', 'Area', (idiso(j),j = 1, niso)
    DO ig = 1, ng, 1
        WRITE(io16, '(I5, <niso+2>(1pe16.6))') ig, Flux(ig), myFxr%area, (OutXs(j, ig, i), j = 1, niso)
    ENDDO

    WRITE(io16, '(A)') 'NU-FISSION'
    WRITE(io16, '(A5, 2A16, <niso>I16)') 'GRP', 'Flux', 'Area', (idiso(j),j = 1, niso)
    DO ig = 1, ng, 1
        WRITE(io16, '(I5, <niso+2>(1pe16.6))') ig, Flux(ig), myFxr%area, (OutXsNF(j, ig, i), j = 1, niso)
    ENDDO

    WRITE(io16, '(A)') 'FISSION'
    WRITE(io16, '(A5, 2A16, <niso>I16)') 'GRP', 'Flux', 'Area', (idiso(j),j = 1, niso)
    DO ig = 1, ng, 1
        WRITE(io16, '(I5, <niso+2>(1pe16.6))') ig, Flux(ig), myFxr%area, (OutXsF(j, ig, i), j = 1, niso)
    ENDDO
    WRITE(io16, '(A)') '------------------------------------------------------------------------------------------'
#else
    nisop2 = niso + 2
    READ(c_niso, *) niso; READ(c_nisop2, *) nisop2

    WRITE(io16, '(37X, ' // TRIM(c_niso) // 'ES16.6)') (pnum(j),j=1, niso)
    WRITE(io16, '(A5, 2A16, ' // TRIM(c_niso) // 'I16)') 'GRP', 'Flux', 'Area', (idiso(j),j = 1, niso)
    DO ig = 1, ng, 1
        WRITE(io16, '(I5, ' // TRIM(c_nisop2) // '(1pe16.6))') ig, Flux(ig), myFxr%area, (OutXs(j, ig, i), j = 1, niso)
    ENDDO

    WRITE(io16, '(A)') 'NU-FISSION'
    WRITE(io16, '(A5, 2A16, ' // TRIM(c_niso) // 'I16)') 'GRP', 'Flux', 'Area', (idiso(j),j = 1, niso)
    DO ig = 1, ng, 1
        WRITE(io16, '(I5, ' // TRIM(c_nisop2) // '(1pe16.6))') ig, Flux(ig), myFxr%area, (OutXsNF(j, ig, i), j = 1, niso)
    ENDDO

    WRITE(io16, '(A)') 'FISSION'
    WRITE(io16, '(A5, 2A16, ' // TRIM(c_niso) // 'I16)') 'GRP', 'Flux', 'Area', (idiso(j),j = 1, niso)
    DO ig = 1, ng, 1
        WRITE(io16, '(I5, ' // TRIM(c_nisop2) // '(1pe16.6))') ig, Flux(ig), myFxr%area, (OutXsF(j, ig, i), j = 1, niso)
    ENDDO
    WRITE(io16, '(A)') '------------------------------------------------------------------------------------------'
#endif
  ENDDO !ir sweep


  NULLIFY(IDISO, pnum, myFxr)
ENDDO
CALL ReturnXsMacDat(XsMac)

close(io16)
END SUBROUTINE
