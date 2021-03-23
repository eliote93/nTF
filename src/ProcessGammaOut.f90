#include <defines.h>
SUBROUTINE ProcessGEffXs(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF
USE CNTL,             ONLY : nTracerCntl_Type,    OutpCntl_Type
USE FILES,            ONLY : io16,                caseid
USE ioutil
USE BasicOperation
!USE MacXsLib_mod,   ONLY : EffRIFPin
!USE TH_Mod,         ONLY : GetPinFuelTemp
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Pin_Type), POINTER :: PIN(:)
!TYPE(ResVarPin_Type), POINTER :: mypin

INTEGER :: iz, ixa, iya, ixc, iyc, nxc, nyc
INTEGER :: ixy, ixya, ixyc, iasytype, icel
INTEGER :: i, j, ig, ng, ir1, ir2, igg
REAL :: pinArea!, PinFuelTempAvgsq
REAL, POINTER :: pinOutXSA(:, :), pinOutXSS(:, :), pinOutXSStr(:, :), PinOutXSPP(:, :, :), PinOutPPM(:, :, :), PinOutPSM(:, :, :)
REAL, POINTER :: pinOutKERMA(:, :)
REAL, POINTER :: pinOutPFlux(:), pinOutNFlux(:)
REAL, POINTER :: pinNum(:)
!Output Data
REAL, POINTER :: OutXSA(:, :, :),OutXSS(:, :, :),OutXSStr(:, :, :),OutKERMA(:, :, :), OutXSPP(:, :, :, :), OutPPM(:, :, :, :), OutPSM(:, :, :, :)
REAL, POINTER :: OutPFlux(:, :), OutNFlux(:, :)
REAL, POINTER :: Area(:), Num(:, :)
LOGICAL :: lAIC
CHARACTER(256) :: fn
INTEGER :: ngg, imt

IF(.NOT. nTracerCntl%lXsLib) RETURN

Pin => Core%Pin; CellInfo => Core%CellInfo;
!ResVarPin => Core%ResVarPin
ng = GroupInfo%ng
ngg = GroupInfo%ngg

ALLOCATE(OutXSA(200, ngg, nTracerCntl%OutpCntl%nRegPhXsOut))
ALLOCATE(OutXSS(200, ngg, nTracerCntl%OutpCntl%nRegPhXsOut))
ALLOCATE(OutXSStr(200, ngg, nTracerCntl%OutpCntl%nRegPhXsOut))
ALLOCATE(OutKERMA(200, ngg, nTracerCntl%OutpCntl%nRegPhXsOut))
ALLOCATE(OutXSPP(200, ng, 0:4, nTracerCntl%OutpCntl%nRegPhXsOut))
ALLOCATE(OutPPM(200, ngg, ng, nTracerCntl%OutpCntl%nRegPhXsOut))
ALLOCATE(OutPSM(200, ngg, ngg, nTracerCntl%OutpCntl%nRegPhXsOut))
ALLOCATE(OutPFlux(ngg, nTracerCntl%OutpCntl%nRegPhXsOut))
ALLOCATE(OutNFlux(ng, nTracerCntl%OutpCntl%nRegPhXsOut))
ALLOCATE(Area(nTracerCntl%OutpCntl%nRegPhXsOut))
ALLOCATE(NUM(200, nTracerCntl%OutpCntl%nRegPhXsOut))

DO i = 1, nTracerCntl%OutpCntl%nRegPhXsOut
  CALL CP_CA(OutXsA(1:nTracerCntl%OutpCntl%IsoPhXsOutList(0, i), 1:ngg, i), 0._8, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i), ngg)
  CALL CP_CA(OutXsS(1:nTracerCntl%OutpCntl%IsoPhXsOutList(0, i), 1:ngg, i), 0._8, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i), ngg)
  CALL CP_CA(OutXsStr(1:nTracerCntl%OutpCntl%IsoPhXsOutList(0, i), 1:ngg, i), 0._8, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i), ngg)
  CALL CP_CA(OutKERMA(1:nTracerCntl%OutpCntl%IsoPhXsOutList(0, i), 1:ngg, i), 0._8, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i), ngg)
  OutXSPP(1:nTracerCntl%OutpCntl%IsoPhXsOutList(0, i), 1:ng, 0:4, i) = 0._8
  OutPPM(1:nTracerCntl%OutpCntl%IsoPhXsOutList(0, i), 1:ngg, 1:ng, i) = 0._8
  OutPSM(1:nTracerCntl%OutpCntl%IsoPhXsOutList(0, i), 1:ngg, 1:ngg, i) = 0._8
  CALL CP_CA(OutPFlux(1:ngg, i), 0._8, ngg)
  CALL CP_CA(OutNFlux(1:ng, i), 0._8, ng)
  CALL CP_CA(NUM(1:200, i), 0._8, 200)
  iz = nTracerCntl%OutpCntl%RegPhXsOutList(1, i)
  ixa = nTracerCntl%OutpCntl%RegPhXsOutList(2, i); iya = nTracerCntl%OutpCntl%RegPhXsOutList(3, i)
  ixya = Core%CoreIdx(ixa, iya); iasytype = Core%CoreMap(ixya)
  IF (nTracerCntl%OutpCntl%RegPhXsOutASM(i)) THEN
      ALLOCATE(pinOutXSA(200, ngg),pinOutXSS(200, ngg),pinOutXSStr(200, ngg),pinOutKERMA(200, ngg),pinOutPFlux(ngg),pinNum(200))
      ALLOCATE(pinOutXSPP(200, ng, 0:4), pinOutNFlux(ng))
      ALLOCATE(PinOutPPM(200, ngg, ng))
      ALLOCATE(PinOutPSM(200, ngg, ngg))
      nxc = Core%AsyInfo(iasytype)%nx; nyc = Core%AsyInfo(iasytype)%ny;
      AREA(i)=0
      DO iyc=1,nyc
         DO ixc=1,nxc
             ixyc = Core%AsyInfo(iasyType)%pin2DIdx(ixc, iyc)
             ixy = Core%Asy(ixya)%GlobalPinIdx(ixyc); icel = Core%Pin(ixy)%Cell(iz)
             !mypin => ResVarPin(ixy,iz); lAIC = CellInfo(icel)%lAIC
             !PinFuelTempAvgsq = dsqrt(GetPinFuelTemp(Core, FmInfo%Fxr, iz, ixy))
             !if (mypin%lres.and.nTracerCntl%lRIF) CALL EffRIFPin(mypin, PinFuelTempAvgsq, iz, lAIC, PE)
             ir1=1; ir2=CellInfo(icel)%nFXR
             CALL FxrPhXSCnds_Numerator(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE, pinOutXSA, &
             pinOutXSS,pinOutXSStr,pinOutKERMA, PinOutPSM,pinOutPFlux,pinNum,pinArea,i,ir1,ir2,ixy,iz,ngg)
             OutXsA(:,:,i)=OutXsA(:,:,i)+pinOutXSA
             OutXsS(:,:,i)=OutXsS(:,:,i)+pinOutXSS
             OutXsStr(:,:,i)=OutXsStr(:,:,i)+pinOutXSStr
             OutKERMA(:,:,i)=OutKERMA(:,:,i)+pinOutKERMA
             OutPFlux(:,i)=OutPFlux(:,i)+pinOutPFlux
             Num(:,i)=Num(:,i)+pinNum
             AREA(i)=AREA(i)+pinArea
             CALL FxrPPXSCnds_Numerator(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE, pinOutXSPP, PinOutPPM,pinOutNFlux,pinNum,pinArea,i,ir1,ir2,ixy,iz,ng,ngg)
             OutXsPP(:, :, :, i) = OutXsPP(:, :, :, i) + pinOutXSPP
             OutPPM(:, :, :, i) = OutPPM(:, :, :, i) + PinOutPPM
             OutNFlux(:,i)=OutNFlux(:,i)+pinOutNFlux
         ENDDO
      ENDDO
      DEALLOCATE(pinOutXSA,pinOutXSS,pinOutXSStr,pinOutKERMA,pinOutPFlux,pinNUM)
  ELSEIF (nTracerCntl%OutpCntl%RegPhXsOutPIN(i)) THEN
      ixc = nTracerCntl%OutpCntl%RegPhXsOutList(4, i); iyc = nTracerCntl%OutpCntl%RegPhXsOutList(5, i)
      if (ixc.ne.0.and.iyc.ne.0) then
          ixyc = Core%AsyInfo(iasyType)%pin2DIdx(ixc, iyc)
          ixy = Core%Asy(ixya)%GlobalPinIdx(ixyc); icel = Pin(ixy)%Cell(iz)
          ir1=1; ir2=CellInfo(icel)%nFXR
          CALL FxrPhXSCnds_Numerator(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE, OutXSA(:,:,i), &
          OutXSS(:,:,i), OutXSStr(:,:,i),OutKERMA(:,:,i),OutPSM(:,:,:,i),OutPFlux(:,i),Num(:,i),AREA(i),i,ir1,ir2,ixy,iz,ngg)
      elseif (ixc.eq.0) then
          ALLOCATE(pinOutXSA(200, ngg),pinOutXSS(200, ngg),pinOutXSStr(200, ngg),&
              pinOutKERMA(200, ngg),pinOutPFlux(ngg),pinNum(200))
          nxc=Core%AsyInfo(iasyType)%nx
          AREA(i)=0
          do ixc=1,nxc
              ixyc = Core%AsyInfo(iasyType)%pin2DIdx(ixc, iyc)
              ixy = Core%Asy(ixya)%GlobalPinIdx(ixyc); icel = Pin(ixy)%Cell(iz)
              ir1=1; ir2=CellInfo(icel)%nFXR
              CALL FxrPhXSCnds_Numerator(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE, pinOutXSA,pinOutXSS,&
                  pinOutXSStr,pinOutKERMA,PinOutPSM,pinOutPFlux,pinNum,pinArea,i,ir1,ir2,ixy,iz,ngg)
              OutXsA(:,:,i)=OutXsA(:,:,i)+pinOutXSA
              OutXsS(:,:,i)=OutXsS(:,:,i)+pinOutXSS
              OutXsStr(:,:,i)=OutXsStr(:,:,i)+pinOutXSStr
              OutKERMA(:,:,i)=OutKERMA(:,:,i)+pinOutKERMA
              OutPFlux(:,i)=OutPFlux(:,i)+pinOutPFlux
              Num(:,i)=Num(:,i)+pinNum
              AREA(i)=AREA(i)+pinArea
             CALL FxrPPXSCnds_Numerator(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE, pinOutXSPP, PinOutPPM,pinOutNFlux,pinNum,pinArea,i,ir1,ir2,ixy,iz,ng,ngg)
             OutXsPP(:, :, :, i) = OutXsPP(:, :, :, i) + pinOutXSPP
             OutPPM(:, :, :, i) = OutPPM(:, :, :, i) + PinOutPPM
             OutNFlux(:,i)=OutNFlux(:,i)+pinOutNFlux
          enddo
          DEALLOCATE(pinOutXSA,pinOutXSS,pinOutXSStr,pinOutKERMA,pinOutPFlux,pinNUM)
      else
          ALLOCATE(pinOutXSA(200, ngg),pinOutXSS(200, ngg),pinOutXSStr(200, ngg),pinOutKERMA(200, ngg),pinOutPFlux(ngg),pinNum(200))
          nyc=Core%AsyInfo(iasyType)%ny
          AREA(i)=0
          do iyc=1,nyc
              ixyc = Core%AsyInfo(iasyType)%pin2DIdx(ixc, iyc)
              ixy = Core%Asy(ixya)%GlobalPinIdx(ixyc); icel = Pin(ixy)%Cell(iz)
              ir1=1; ir2=CellInfo(icel)%nFXR
              CALL FxrPhXSCnds_Numerator(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE, pinOutXSA, pinOutXSS, &
                  pinOutXSStr,pinOutKERMA,PinOutPSM,pinOutPFlux,pinNum,pinArea,i,ir1,ir2,ixy,iz,ngg)
              OutXsA(:,:,i)=OutXsA(:,:,i)+pinOutXSA
              OutXsS(:,:,i)=OutXsS(:,:,i)+pinOutXSS
              OutXsStr(:,:,i)=OutXsStr(:,:,i)+pinOutXSStr
              OutKERMA(:,:,i)=OutKERMA(:,:,i)+pinOutKERMA
              OutPFlux(:,i)=OutPFlux(:,i)+pinOutPFlux
              Num(:,i)=Num(:,i)+pinNum
              AREA(i)=AREA(i)+pinArea
             CALL FxrPPXSCnds_Numerator(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE, pinOutXSPP, PinOutPPM,pinOutNFlux,pinNum,pinArea,i,ir1,ir2,ixy,iz,ng,ngg)
             OutXsPP(:, :, :, i) = OutXsPP(:, :, :, i) + pinOutXSPP
             OutPPM(:, :, :, i) = OutPPM(:, :, :, i) + PinOutPPM
             OutNFlux(:,i)=OutNFlux(:,i)+pinOutNFlux
          enddo
          DEALLOCATE(pinOutXSA,pinOutXSS,pinOutXSStr,pinOutKERMA,pinOutPFlux,pinNUM)
      endif
  ELSEIF (nTracerCntl%OutpCntl%RegPhXsOutFXR(i)) THEN
      ixc = nTracerCntl%OutpCntl%RegPhXsOutList(4, i); iyc = nTracerCntl%OutpCntl%RegPhXsOutList(5, i)
      ixyc = Core%AsyInfo(iasyType)%pin2DIdx(ixc, iyc)
      ixy = Core%Asy(ixya)%GlobalPinIdx(ixyc); icel = Core%Pin(ixy)%Cell(iz)
      !mypin => ResVarPin(ixy,iz); lAIC = CellInfo(icel)%lAIC
      !PinFuelTempAvgsq = dsqrt(GetPinFuelTemp(Core, FmInfo%Fxr, iz, ixy))
      !IF (mypin%lres.and.nTracerCntl%lRIF) CALL EffRIFPin(mypin, PinFuelTempAvgsq, iz, lAIC, PE)
      ir1=nTracerCntl%OutpCntl%RegPhXsOutList(6, i)
      ir2=nTracerCntl%OutpCntl%RegPhXsOutList(7, i)
      CALL FxrPhXSCnds_Numerator(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE, OutXSA(:,:,i), OutXSS(:,:,i),&
          OutXSStr(:,:,i),OutKERMA(:,:,i),OutPSM(:,:,:,i),OutPFlux(:,i),Num(:,i),AREA(i),i,ir1,ir2,ixy,iz,ngg)
      CALL FxrPPXSCnds_Numerator(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE, OutXSPP(:,:,:,i), OutPPM(:, :, :, i),OutNFlux(:,i),Num(:,i),AREA(i),i,ir1,ir2,ixy,iz,ng,ngg)
  ENDIF
  DO j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i)
    NUM(j, i) = NUM(j, i) / Area(i)
    DO ig =  1, ngg
      OutXsA(j, ig, i) = OutXsA(j, ig, i) / OutPFlux(ig, i)
      OutXsA(j, ig, i) = OutXsA(j, ig, i) / NUM(j, i)
      OutXsS(j, ig, i) = OutXsS(j, ig, i) / OutPFlux(ig, i)
      OutXsS(j, ig, i) = OutXsS(j, ig, i) / NUM(j, i)
      OutXsStr(j, ig, i) = OutXsStr(j, ig, i) / OutPFlux(ig, i)
      OutXsStr(j, ig, i) = OutXsStr(j, ig, i) / NUM(j, i)
      OutKERMA(j, ig, i) = OutKERMA(j, ig, i) / OutPFlux(ig, i)
      OutKERMA(j, ig, i) = OutKERMA(j, ig, i) / NUM(j, i)
    ENDDO
    DO imt = 0, 4
      DO ig = 1, ng
        OutXsPP(j, ig, imt, i) = OutXsPP(j, ig, imt, i) / OutNFlux(ig, i) / NUM(j, i)
      END DO
    END DO
    DO ig = 1, ngg
      DO igg = 1, ngg
        OutPSM(j, igg, ig, i) = OutPSM(j, igg, ig, i) / OutPFlux(ig, i) / NUM(j, i)
      END DO
    END DO
    DO ig = 1, ng
      DO igg = 1, ngg
        OutPPM(j, igg, ig, i) = OutPPM(j, igg, ig, i) / OutNFlux(ig, i) / NUM(j, i)
      END DO
    END DO
  ENDDO
  DO ig= 1, ngg
    OutPFlux(ig, i) = OutPFlux(ig, i) / Area(i)
  ENDDO
  DO ig = 1, ng
    OutNFlux(ig, i) = OutNFlux(ig, i) / Area(i)
  END DO
ENDDO
! Photoatomic XS Tally
fn=trim(caseid)//'.gexs'

CALL openfile(io16,FALSE,FALSE,FALSE, fn)
DO i = 1, nTracerCntl%OutpCntl%nRegPhXsOut
  IF (nTracerCntl%OutpCntl%RegPhXsOutASM(i)) THEN
      WRITE(io16, '(A10, I5, A3, 1(I5, I5, A3),200I7)') 'Region :',  nTracerCntl%OutpCntl%RegPhXsOutList(1, i), '/', &
                                             nTracerCntl%OutpCntl%RegPhXsOutList(2:3, i), '/', (nTracerCntl%OutpCntl%IsoPhXsOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i))
  ELSEIF (nTracerCntl%OutpCntl%RegPhXsOutPIN(i)) THEN
      WRITE(io16, '(A10, I5, A3, 2(I5, I5, A3),200I7)') 'Region :',  nTracerCntl%OutpCntl%RegPhXsOutList(1, i), '/', &
                                             nTracerCntl%OutpCntl%RegPhXsOutList(2:3, i), '/', &
                                             nTracerCntl%OutpCntl%RegPhXsOutList(4:5, i), '/', (nTracerCntl%OutpCntl%IsoPhXsOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i))
  ELSEIF (nTracerCntl%OutpCntl%RegPhXsOutFXR(i)) THEN
      WRITE(io16, '(A10, I5, A3, 3(I5, I5, A3),200I7)') 'Region :',  nTracerCntl%OutpCntl%RegPhXsOutList(1, i), '/', &
                                             nTracerCntl%OutpCntl%RegPhXsOutList(2:3, i), '/', &
                                             nTracerCntl%OutpCntl%RegPhXsOutList(4:5, i), '/', &
                                             nTracerCntl%OutpCntl%RegPhXsOutList(6:7, i), '/', (nTracerCntl%OutpCntl%IsoPhXsOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i))
  ENDIF
  WRITE(io16, '(A5, 2A16, 600I16)') 'GRP', 'PFlux', 'Area', &
      (nTracerCntl%OutpCntl%IsoPhXsOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i)),&
      (nTracerCntl%OutpCntl%IsoPhXsOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i)),&
      (nTracerCntl%OutpCntl%IsoPhXsOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i)),&
      (nTracerCntl%OutpCntl%IsoPhXsOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i))

  DO ig = 1, ngg, 1
    WRITE(io16, '(I5, 600(1pe16.6))') ig, OutPFlux(ig, i), Area(i),  &
                                (OutKERMA(j, ig, i), j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i)), &
                                (OutXsA(j, ig, i), j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i)), &
                                (OutXsS(j, ig, i), j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i)), &
                                (OutXsStr(j, ig, i), j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i))

  ENDDO
  DO j = 1, (37+16*nTracerCntl%OutpCntl%IsoPhXsOutList(0, i))
    WRITE(io16, '(A)', advance='NO') '-'
  ENDDO
  WRITE(io16, *)
  WRITE(io16, '(A35, 50(1pe15.5))') 'N.D :', (NUM(j, i), j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i))
  WRITE(io16, *)
ENDDO
close(io16)
! Scattering Matrices of Photon
fn=trim(caseid)//'.sm'
CALL openfile(io16,FALSE,FALSE,FALSE, fn)
DO i = 1, nTracerCntl%OutpCntl%nRegPhXsOut
  DO j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i)
    IF (nTracerCntl%OutpCntl%RegPhXsOutASM(i)) THEN
        WRITE(io16, '(A10, I5, A3, 1(I5, I5, A3),I7)') 'SMAT  :',  nTracerCntl%OutpCntl%RegPhXsOutList(1, i), '/', &
                                               nTracerCntl%OutpCntl%RegPhXsOutList(2:3, i), '/', nTracerCntl%OutpCntl%IsoPhXsOutList(j, i)
    ELSEIF (nTracerCntl%OutpCntl%RegPhXsOutPIN(i)) THEN
        WRITE(io16, '(A10, I5, A3, 2(I5, I5, A3),I7)') 'SMAT :',  nTracerCntl%OutpCntl%RegPhXsOutList(1, i), '/', &
                                               nTracerCntl%OutpCntl%RegPhXsOutList(2:3, i), '/', &
                                               nTracerCntl%OutpCntl%RegPhXsOutList(4:5, i), '/', nTracerCntl%OutpCntl%IsoPhXsOutList(j, i)
    ELSEIF (nTracerCntl%OutpCntl%RegPhXsOutFXR(i)) THEN
        WRITE(io16, '(A10, I5, A3, 3(I5, I5, A3),I7)') 'SMAT :',  nTracerCntl%OutpCntl%RegPhXsOutList(1, i), '/', &
                                               nTracerCntl%OutpCntl%RegPhXsOutList(2:3, i), '/', &
                                               nTracerCntl%OutpCntl%RegPhXsOutList(4:5, i), '/', &
                                               nTracerCntl%OutpCntl%RegPhXsOutList(6:7, i), '/', nTracerCntl%OutpCntl%IsoPhXsOutList(j, i)
    ENDIF
    
    WRITE(io16, '(A5, 4A16)') 'GRP', 'NFlux', 'Area', 'SIGS', 'SIGSTR'

    DO ig = 1, ngg
      WRITE(io16, '(I5, 600(1pe16.6))') ig, OutPFlux(ig, i), Area(i), OutXsS(j, ig, i), OutXsStr(j, ig, i), (OutPSM(j, ig, igg, i), igg = 1, ngg)
    ENDDO
  END DO
  DO j = 1, (37+16*ngg)
    WRITE(io16, '(A)', advance='NO') '-'
  ENDDO
  WRITE(io16, *)
  WRITE(io16, '(A35, 50(1pe15.5))') 'N.D :', (NUM(j, i), j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i))
  WRITE(io16, *)
END DO
CLOSE(io16)
! Photon Production XS of each reactions
fn=trim(caseid)//'.ppxs'

CALL openfile(io16,FALSE,FALSE,FALSE, fn)
DO i = 1, nTracerCntl%OutpCntl%nRegPhXsOut
  IF (nTracerCntl%OutpCntl%RegPhXsOutASM(i)) THEN
      WRITE(io16, '(A10, I5, A3, 1(I5, I5, A3),200I7)') 'Region :',  nTracerCntl%OutpCntl%RegPhXsOutList(1, i), '/', &
                                             nTracerCntl%OutpCntl%RegPhXsOutList(2:3, i), '/', (nTracerCntl%OutpCntl%IsoPhXsOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i))
  ELSEIF (nTracerCntl%OutpCntl%RegPhXsOutPIN(i)) THEN
      WRITE(io16, '(A10, I5, A3, 2(I5, I5, A3),200I7)') 'Region :',  nTracerCntl%OutpCntl%RegPhXsOutList(1, i), '/', &
                                             nTracerCntl%OutpCntl%RegPhXsOutList(2:3, i), '/', &
                                             nTracerCntl%OutpCntl%RegPhXsOutList(4:5, i), '/', (nTracerCntl%OutpCntl%IsoPhXsOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i))
  ELSEIF (nTracerCntl%OutpCntl%RegPhXsOutFXR(i)) THEN
      WRITE(io16, '(A10, I5, A3, 3(I5, I5, A3),200I7)') 'Region :',  nTracerCntl%OutpCntl%RegPhXsOutList(1, i), '/', &
                                             nTracerCntl%OutpCntl%RegPhXsOutList(2:3, i), '/', &
                                             nTracerCntl%OutpCntl%RegPhXsOutList(4:5, i), '/', &
                                             nTracerCntl%OutpCntl%RegPhXsOutList(6:7, i), '/', (nTracerCntl%OutpCntl%IsoPhXsOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i))
  ENDIF
  WRITE(io16, '(A5, 2A16, 600I16)') 'GRP', 'NFlux', 'Area', &
      (nTracerCntl%OutpCntl%IsoPhXsOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i)),&
      (nTracerCntl%OutpCntl%IsoPhXsOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i)),&
      (nTracerCntl%OutpCntl%IsoPhXsOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i)),&
      (nTracerCntl%OutpCntl%IsoPhXsOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i)),&
      (nTracerCntl%OutpCntl%IsoPhXsOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i))

  DO ig = 1, ng
    WRITE(io16, '(I5, 600(1pe16.6))') ig, OutNFlux(ig, i), Area(i),  &
                                (OutXSPP(j, ig, 0, i), j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i)), &
                                (OutXSPP(j, ig, 1, i), j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i)), &
                                (OutXSPP(j, ig, 2, i), j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i)), &
                                (OutXSPP(j, ig, 3, i), j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i)), &
                                (OutXSPP(j, ig, 4, i), j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i))

  ENDDO
  DO j = 1, (37+16*nTracerCntl%OutpCntl%IsoPhXsOutList(0, i))
    WRITE(io16, '(A)', advance='NO') '-'
  ENDDO
  WRITE(io16, *)
  WRITE(io16, '(A35, 50(1pe15.5))') 'N.D :', (NUM(j, i), j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i))
  WRITE(io16, *)
ENDDO
close(io16)
! Photon Production Matrices
fn=trim(caseid)//'.ppm'
CALL openfile(io16,FALSE,FALSE,FALSE, fn)
DO i = 1, nTracerCntl%OutpCntl%nRegPhXsOut
  DO j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i)
    IF (nTracerCntl%OutpCntl%RegPhXsOutASM(i)) THEN
        WRITE(io16, '(A10, I5, A3, 1(I5, I5, A3),I7)') 'PPMAT  :',  nTracerCntl%OutpCntl%RegPhXsOutList(1, i), '/', &
                                               nTracerCntl%OutpCntl%RegPhXsOutList(2:3, i), '/', nTracerCntl%OutpCntl%IsoPhXsOutList(j, i)
    ELSEIF (nTracerCntl%OutpCntl%RegPhXsOutPIN(i)) THEN
        WRITE(io16, '(A10, I5, A3, 2(I5, I5, A3),I7)') 'PPMAT :',  nTracerCntl%OutpCntl%RegPhXsOutList(1, i), '/', &
                                               nTracerCntl%OutpCntl%RegPhXsOutList(2:3, i), '/', &
                                               nTracerCntl%OutpCntl%RegPhXsOutList(4:5, i), '/', nTracerCntl%OutpCntl%IsoPhXsOutList(j, i)
    ELSEIF (nTracerCntl%OutpCntl%RegPhXsOutFXR(i)) THEN
        WRITE(io16, '(A10, I5, A3, 3(I5, I5, A3),I7)') 'PPMAT :',  nTracerCntl%OutpCntl%RegPhXsOutList(1, i), '/', &
                                               nTracerCntl%OutpCntl%RegPhXsOutList(2:3, i), '/', &
                                               nTracerCntl%OutpCntl%RegPhXsOutList(4:5, i), '/', &
                                               nTracerCntl%OutpCntl%RegPhXsOutList(6:7, i), '/', nTracerCntl%OutpCntl%IsoPhXsOutList(j, i)
    ENDIF
    
    WRITE(io16, '(A5, 3A16)') 'GRP', 'NFlux', 'Area', 'TOTPPXS'

    DO ig = 1, ng
      WRITE(io16, '(I5, 600(1pe16.6))') ig, OutNFlux(ig, i), Area(i), OutXSPP(j, ig, 0, i), (OutPPM(j, igg, ig, i), igg = 1, ngg)
    ENDDO
  END DO
  DO j = 1, (37+16*ngg)
    WRITE(io16, '(A)', advance='NO') '-'
  ENDDO
  WRITE(io16, *)
  WRITE(io16, '(A35, 50(1pe15.5))') 'N.D :', (NUM(j, i), j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i))
  WRITE(io16, *)
END DO

CLOSE(io16)

DEALLOCATE(OutXSA)
DEALLOCATE(OutXSS)
DEALLOCATE(OutXSStr)
DEALLOCATE(OutKERMA)
DEALLOCATE(OutPFlux)
DEALLOCATE(OutPPM, OutXSPP)
DEALLOCATE(Area)
DEALLOCATE(NUM)
END SUBROUTINE

SUBROUTINE FxrPhXSCnds_Numerator(Core, FmInfo, ThInfo, GroupInfo,nTracerCntl, PE, OutXSA, OutXSS, OutXSStr,OutKERMA,OutPSM,OutFlux,Num,AREA,i,ir1,ir2,ixy,iz,ngg)
    USE PARAM, only : CKELVIN
    USE TYPEDEF
    USE GammaTYPEDEF,     ONLY : GamMacXS_TYPE
    USE CNTL,             ONLY : nTracerCntl_Type
    USE GamXsLib_Mod,     ONLY : GamXsBase, GamScatMatrix
    USE GamMOC_MOD,       ONLY : FxrAvgGPhi
    USE GamXSUtil,        ONLY : GetGamXsMacDat,         ReturngamXsMacDat
    USE GammaCore_mod,    ONLY : Gphis
    USE ioutil
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FmInfo_Type) :: FmInfo
    TYPE(ThInfo_Type) :: ThInfo
    TYPE(GroupInfo_Type) :: GroupInfo
    TYPE(FxrInfo_Type), POINTER :: myFxr
    TYPE(GamMacXS_TYPE), POINTER :: XsMac
    TYPE(nTracerCntl_Type) :: nTracerCntl
    TYPE(PE_Type) :: PE

    INTEGER, INTENT(IN) :: i,ir1,ir2,ixy,iz,ngg
    !TYPE(ResVarPin_Type), INTENT(IN) :: mypin
    !REAL,INTENT(IN) :: PinFuelTempAvgsq
    REAL,INTENT(OUT) :: AREA
    REAL,INTENT(OUT) :: OutXSA(200, ngg),OutXSS(200, ngg),OutXSStr(200, ngg),OutKERMA(200, ngg)
    REAL,INTENT(OUT) :: OutFlux(ngg), OutPSM(200, ngg, ngg)
    REAL,INTENT(OUT) :: Num(200)
    REAL :: Flux(ngg)

    LOGICAL :: FlagF

    INTEGER :: ir, j, k, id, ig, iresoGrp1, iresoGrp2, ifxr, jG
    !Fxr Data
    INTEGER :: niso
    LOGICAL :: lres
    INTEGER, POINTER :: idiso(:)
    REAL, POINTER :: pnum(:)
    REAL :: eigv
    REAL :: Temp

    eigv = 1.0
    iresoGrp1 = GroupInfo%nofg + 1; iresoGrp2 = GroupInfo%nofg + GroupInfo%norg
    CALL GetGamXsMacDat(XsMac, ngg, .TRUE.)
    AREA=0; OutFlux=0; OutXSA=0; OutXSS=0; OutXSStr=0; OutKERMA=0; NUM=0; OutPSM = 0
    DO ir= ir1, ir2
        ifxr = Core%Pin(ixy)%FxrIdxst + ir - 1;myFxr => FmInfo%Fxr(ifxr, iz)
        temp = myFxr%temp
        niso = myFxr%niso; idiso => myFxr%idiso
        pnum => myFxr%pnum; lres = myFxr%lres
        CALL GamXsBase(XsMac, myFxr, 1, ngg, ngg, .TRUE.)
        CALL GamScatMatrix(XsMac, myFxr, 1, ngg, ngg, .TRUE.,.TRUE.)
        AREA = AREA + myFxr%area

        Flux = FxrAvgGPhi(Core, FmInfo%Fxr, Gphis, ixy, ir, iz, ngg, PE)
        Outflux = OutFlux + myFxr%area * Flux

        DO j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i)
          id = nTracerCntl%OutpCntl%IsoPhXsOutList(j, i)
          DO k = 1, niso
            IF(IDISO(k) .NE. id) CYCLE
            NUM(j) = NUM(j) +  myFxr%area * pnum(k)
            DO ig =  1, ngg
              OutKERMA(j, ig) = OutKERMA(j, ig)+  myFxr%area * XsMac%IsoKERMA(k, ig) * Flux(ig)
              OutXsA(j, ig) = OutXsA(j, ig)+  myFxr%area * XsMac%IsoXsMacA(k, ig) * Flux(ig)
              OutXsS(j, ig) = OutXsS(j, ig)+  myFxr%area * XsMac%IsoXsMacS(k, ig) * Flux(ig)
              OutXsStr(j, ig) = OutXsStr(j, ig)+  myFxr%area * XsMac%IsoXsMacSTR(k, ig) * Flux(ig)
              DO jg = 1, ngg
                OutPSM(j, jg, ig) = OutPSM(j, jg, ig)+  myFxr%area * XsMac%IsoSM(ig, jg, k) * Flux(ig)
              END DO
            ENDDO
          ENDDO
        ENDDO
    ENDDO
    NULLIFY(IDISO, pnum, myFxr)
    CALL ReturngamXsMacDat(XsMac)

  END SUBROUTINE

  SUBROUTINE FxrPPXSCnds_Numerator(Core, FmInfo, ThInfo, GroupInfo,nTracerCntl, PE, OutXSPP, OutPPM,OutFlux,Num,AREA,i,ir1,ir2,ixy,iz,ng,ngg)
    USE PARAM, only : CKELVIN
    USE TYPEDEF
    USE CNTL,             ONLY : nTracerCntl_Type
    USE MacXsLib_mod,     ONLY : EffMacXS,            MacXsBase
    USE MOC_Mod,          ONLY : FxrAvgPhi
    USE GamXSUtil,        ONLY : GetGamXsMacDat,         ReturnGamXsMacDat
    USE ioutil
    USE GammaTYPEDEF,     ONLY : GamMacXS_TYPE
    USE GamXsLib_Mod,     ONLY : GamProdMatrix, GamScatMatrix
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FmInfo_Type) :: FmInfo
    TYPE(ThInfo_Type) :: ThInfo
    TYPE(GroupInfo_Type) :: GroupInfo
    TYPE(FxrInfo_Type), POINTER :: myFxr
    TYPE(GamMacXS_TYPE), POINTER :: XsMac
    TYPE(nTracerCntl_Type) :: nTracerCntl
    TYPE(PE_Type) :: PE

    INTEGER, INTENT(IN) :: i,ir1,ir2,ixy,iz,ng, ngg
    REAL,INTENT(OUT) :: AREA
    REAL,INTENT(OUT) :: OutXSPP(200, ng, 0:4), OutPPM(200, ngg, ng)
    REAL,INTENT(OUT) :: OutFlux(ng)
    REAL,INTENT(OUT) :: Num(200)
    REAL :: Flux(ng)

    INTEGER :: ir, j, k, id, ig, iresoGrp1, iresoGrp2, ifxr
    !Fxr Data
    INTEGER :: niso
    INTEGER, POINTER :: idiso(:)
    REAL, POINTER :: pnum(:), IsoXsMacfold(:,:), IsoXsMacNfold(:,:)
    
    REAL, POINTER :: IsoProd(:, :, :)
    INTEGER :: imt, igg

    CALL GetGamXsMacDat(XsMac, ng, .TRUE.)
    AREA=0; OutFlux=0; OutXsPP=0; NUM=0; OutPPM = 0
    DO ir= ir1, ir2
        ifxr = Core%Pin(ixy)%FxrIdxst + ir - 1;myFxr => FmInfo%Fxr(ifxr, iz)
        niso = myFxr%niso; idiso => myFxr%idiso
        pnum => myFxr%pnum
        CALL GamProdMatrix(XsMac, myFxr, 1, ngg, ng, ngg, GroupInfo, .TRUE.)
        AREA = AREA + myFxr%area

        Flux = FxrAvgPhi(Core, FmInfo%Fxr, FMInfo%Phis, ixy, ir, iz, ng, PE)
        Outflux = OutFlux + myFxr%area * Flux

        DO j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i)
          id = nTracerCntl%OutpCntl%IsoPhXsOutList(j, i)
          DO k = 1, niso
            IF(IDISO(k) .NE. id) CYCLE
            NUM(j) = NUM(j) +  myFxr%area * pnum(k)
            DO imt = 1, 4
              IF (.NOT.XsMac%lexist(k, imt)) CYCLE
              SELECT CASE(imt)
              CASE(1)
                IsoProd => XSMac%IsoGProdFis
              CASE(2)
                IsoProd => XSMac%IsoGProdRad
              CASE(3)
                IsoProd => XSMac%IsoGProdInel
              CASE(4)
                IsoProd => XSMac%IsoGProdNnel
              END SELECT
              DO ig =  1, ng
                DO igg = 1, ngg
                  OutXsPP(j, ig, imt) = OutXsPP(j, ig, imt) +  myFxr%area * IsoProd(k, ig, igg) * Flux(ig)
                  OutPPM(j, igg, ig) = OutPPM(j, igg, ig) +  myFxr%area * IsoProd(k, ig, igg) * Flux(ig)
                END DO ! igg
              ENDDO ! ig
            END DO ! imt
          ENDDO ! k (niso)
          !PRINT '(3I, 5ES14.6)', j, OutXSPP(j, ig, 0), OutXSPP(j, ig, 1), OutXSPP(j, ig, 2), OutXSPP(j, ig, 3), OutXSPP(j, ig, 4)
        ENDDO ! j (outlist)
    ENDDO ! ir
    
    DO j = 1, nTracerCntl%OutpCntl%IsoPhXsOutList(0, i)
      id = nTracerCntl%OutpCntl%IsoPhXsOutList(j, i)
      DO k = 1, niso
        IF(IDISO(k) .NE. id) CYCLE
        DO imt = 1, 4
          DO ig =  1, ng
            OutXsPP(j, ig, 0) = OutXsPP(j, ig, 0) + OutXsPP(j, ig, imt)
          ENDDO ! ig
        END DO ! imt
      ENDDO ! k (niso)
    ENDDO ! j (outlist)
      
    CALL ReturnGamXsMacDat(XsMac)

  END SUBROUTINE
!SUBROUTINE ProcessGFlux(Core, FmInfo, ThInfo, GroupInfo, nTRACERCntl, PE)
!  END SUBROUTINE
!  
!  SUBROUTINE PrintGammaPower
!  END SUBROUTINE
!  
!  SUBROUTINE 
!  END SUBROUTINE

!SUBROUTINE PrintGammaFlux(io, core, CmInfo, PowerDist, ngg, nTRACERCntl, PE)
!USE PARAM
!USE TYPEDEF,          ONLY : CoreInfo_type,   PowerDist_Type, &
!                             PE_Type,         Asy_type,          AsyInfo_type
!USE CNTL,             ONLY : nTRACERCntl_TYPE
!USE GammaTYPEDEF,  ONLY : GammaCMFD_TYPE
!USE BasicOperation,   ONLY : CP_CA,    CP_VA
!#ifdef MPI_ENV
!USE MPIComm_Mod,      ONLY : REDUCE
!#endif
!IMPLICIT NONE
!TYPE(CoreInfo_Type) :: Core
!TYPE(GammaCMFD_TYPE) :: CmInfo
!TYPE(PowerDist_Type) :: PowerDist
!TYPE(nTracerCntl_Type) :: nTracerCntl
!TYPE(PE_Type) :: PE
!INTEGER :: io, ngg
!
!TYPE(Asy_Type), POINTER :: Asy(:)
!TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
!
!
!CALL AxAvg2DGammaFlux(io, core, CmInfo, ngg, PE)
!CALL Asy3DGammaFlux(io, core, CmInfo, ngg, PE)
!CALL RadAvg1DGammaFlux(io, core, CmInfo, ngg, PE)
!CALL CoreGammaSpectrum(io, core, CmInfo, ngg, PE)
!CALL Pin3DGammaFlux(io, core, CmInfo, ngg, nTracerCntl, PE)
!
!  CONTAINS
  
SUBROUTINE AxAvg2DGammaFlux(io, core, CmInfo, ng, PE)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type,     PE_TYPE,   &
                    AsyInfo_Type,    Asy_Type
USE IOUTIL,  ONLY : PrintReal1DarrayTo2Darray
USE GammaTYPEDEF,  ONLY : GammaCMFD_TYPE
#ifdef MPI_ENV
USE BasicOperation, ONLY : CP_CA, CP_VA
USE MpiComm_mod, ONLY : Reduce
#endif
IMPLICIT NONE
INTEGER :: io, ng
TYPE(CoreInfo_Type) :: Core
TYPE(GammaCMFD_TYPE) :: CMInfo
TYPE(PE_TYPE) :: PE

TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)
REAL, POINTER :: PhiC(:, :, :)
REAL, POINTER :: PinVol(:, :), AsyVol(:, :)
REAL, POINTER :: hz(:)
INTEGER :: nz, nxy, nxya, nxa, nya, myzb, myze
INTEGER :: nAsyType, AsyType
INTEGER :: iz, ig, iasy, ixy, ixa, iya, i, j, k
LOGICAL :: master, slave
REAL,POINTER :: Avg2DFlx(:, :,  :), row(:)
#ifdef MPI_ENV
REAL, POINTER :: Buf(:, :)
#endif

Master = PE%Cmfdmaster; Slave = .NOT. Master

PhiC => CmInfo%GPhiC
AsyInfo => Core%AsyInfo; Asy => Core%Asy
hz => Core%hz;
PinVol => Core%PinVol; AsyVol => Core%AsyVol
nz = Core%nz; myzb = PE%myzb; myze =PE%myze
nAsyType = Core%nAsyType; nxya = Core%nxya;
nxy = Core%nxy; nxa = Core%nxa; nya = Core%nya

ALLOCATE(Avg2DFlx(nxa,nya, ng))
#ifdef MPI_ENV
ALLOCATE(Buf(nxa, nya))
#endif
DO ig = 1, ng
  Avg2DFlx(:, :, ig) = 0
  DO iasy = 1, nxya
    AsyType = Asy(iasy)%AsyType
    ixa = Asy(iasy)%ixa; iya = Asy(iasy)%iya
    nxy = AsyInfo(AsyType)%nxy
    Do iz = myzb, myze  !Axial Sweep
      DO i = 1, nxy  !Cell Sweep
        ixy = Asy(iasy)%GlobalPinIdx(i)
        Avg2DFlx(ixa, iya, ig) = Avg2DFlx(ixa, iya, ig) + Pinvol(ixy, iz) * Phic(ixy, iz, ig)
      ENDDO
    ENDDO
  ENDDO !
#ifdef MPI_ENV
  CALL CP_CA(Buf, zero, nxa, nya)
  CALL REDUCE(Avg2DFlx(:, :, ig), Buf(:, :), nxa, nya, PE%MPI_CMFD_COMM, .FALSE.)
  IF(Master) CALL CP_VA(Avg2DFlx(:, :, ig), Buf, nxa, nya)
#endif
ENDDO

#ifdef MPI_ENV
DEALLOCATE(Buf)
#endif

IF(Master) THEN
  WRITE(io, '(A)') ' -  Axially Averaged 2D Gamma Flux -'
  WRITE(io, *)
  DO iya = 1, nya
    DO ig = 1, ng
      WRITE(io, '(8x, 200(1pe15.4))') (Avg2DFlx(ixa, iya, ig), ixa = 1, nxa)
    ENDDO
    WRITE(io, *)
  ENDDO
ENDIF
DEALLOCATE(Avg2DFlx)

NULLIFY(PhiC)
NULLIFY(AsyInfo); NULLIFY(Asy)
NULLIFY(hz); NULLIFY(PinVol)
NULLIFY(AsyVol)
END SUBROUTINE

SUBROUTINE RadAvg1DGammaFlux(io, core, CmInfo, ng, PE)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type,     PE_TYPE,   &
                    AsyInfo_Type,    Asy_Type
USE IOUTIL,  ONLY : PrintReal1DarrayTo2Darray
USE GammaTYPEDEF,  ONLY : GammaCMFD_TYPE
#ifdef MPI_ENV
USE MpiComm_mod, ONLY : SENDRECV
#endif
IMPLICIT NONE
INTEGER :: io, ng
TYPE(CoreInfo_Type) :: Core
TYPE(GammaCMFD_TYPE) :: CMInfo
TYPE(PE_TYPE) :: PE

TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)
REAL, POINTER :: PhiC(:, :, :)
REAL, POINTER :: PinVol(:, :), AsyVol(:, :)
REAL, POINTER :: hz(:)
INTEGER :: nz, nxy, nxya, nxa, nya, myzb, myze
INTEGER :: nAsyType, AsyType
INTEGER :: iz, ig, iasy, ixy, ixa, iya, i, j, k
LOGICAL :: master, slave
REAL,POINTER :: Avg1DFlx(:, :)
#ifdef MPI_ENV
INTEGER :: nproc, comm
INTEGER :: iz1, iz2
#endif

master = PE%CmfdMaster; slave = .NOT. master

PhiC => CmInfo%GPhiC
AsyInfo => Core%AsyInfo; Asy => Core%Asy
hz => Core%hz;
PinVol => Core%PinVol; AsyVol => Core%AsyVol
nz = Core%nz
myzb = PE%myzb; myze = PE%myze
nAsyType = Core%nAsyType; nxya = Core%nxya;
nxy = Core%nxy; nxa = Core%nxa; nya = Core%nya

ALLOCATE(Avg1DFlx(ng, nz))

Avg1DFlx = 0
DO ig = 1, ng
  DO iasy = 1, nxya
    AsyType = Asy(iasy)%AsyType
    nxy = AsyInfo(AsyType)%nxy
    Do iz = myzb, myze   !Axial Sweep
      DO i = 1, nxy  !Cell Sweep
        ixy = Asy(iasy)%GlobalPinIdx(i)
        Avg1DFlx(ig, iz) = Avg1DFlx(ig, iz) + Pinvol(ixy, iz) * Phic(ixy, iz, ig)
      ENDDO
    ENDDO
  ENDDO !
ENDDO

#ifdef MPI_ENV
nproc = PE%nCmfdProc; Comm = PE%MPI_CMFD_COMM
IF(Master) THEN
  DO i = 1, nproc - 1
    iz1 = PE%AxDomRange(1, i); iz2 = PE%AxDomRange(2, i)
    CALL SendRecv(Avg1DFlx(:, iz1:iz2), ng, iz2 - iz1 + 1, i, FALSE, COMM)
  ENDDO
ELSE
  CALL SendRecv(Avg1DFlx(:, myzb:myze), ng, myze - myzb + 1, 0, TRUE, COMM)
ENDIF
#endif

IF(master) THEN
  WRITE(io, '(A)') ' -  Radially Averaged 1D Gamma Flux -'
  WRITE(io, *)
  WRITE(io, '(A)') '                --- > Plane #'
  WRITE(io, '(2x,A,200(I10,5X))') 'Energy Group',(i, i=1,nz)
  DO ig = 1, ng
    WRITE(io, '(I10,4x, 200(1pe15.4))') ig, (Avg1DFlx(ig, iz), iz = 1, nz)
  ENDDO
  WRITE(io, *)
ENDIF

NULLIFY(PhiC)
NULLIFY(AsyInfo); NULLIFY(Asy)
NULLIFY(hz); NULLIFY(PinVol)
NULLIFY(AsyVol)

END SUBROUTINE

SUBROUTINE Asy3DGammaFlux(io, core, CmInfo, ng, PE)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type,     PE_TYPE,   &
                    AsyInfo_Type,    Asy_Type
USE IOUTIL,  ONLY : PrintReal1DarrayTo2Darray
USE GammaTYPEDEF,  ONLY : GammaCMFD_TYPE
#ifdef MPI_ENV
USE MpiComm_mod, ONLY : SENDRECV
#endif
IMPLICIT NONE
INTEGER :: io, ng
TYPE(CoreInfo_Type) :: Core
TYPE(GammaCMFD_TYPE) :: CMInfo
TYPE(PE_TYPE) :: PE

TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)
REAL, POINTER :: PhiC(:, :, :)
REAL, POINTER :: PinVol(:, :), AsyVol(:, :)
REAL, POINTER :: hz(:)
INTEGER :: nz, nxy, nxya, nxa, nya, myzb, myze
INTEGER :: nAsyType, AsyType
INTEGER :: iz, ig, iasy, ixy, ixa, iya, i, j, k
LOGICAL :: master, slave
REAL,POINTER :: Asy3DFlx(:,:,:,:), row(:)
#ifdef MPI_ENV
INTEGER :: COMM, nproc
INTEGER :: iz1, iz2
#endif

Master = PE%CmfdMaster; Slave = .NOT. Master
PhiC => CmInfo%GPhiC
AsyInfo => Core%AsyInfo; Asy => Core%Asy
hz => Core%hz;
PinVol => Core%PinVol; AsyVol => Core%AsyVol
nz = Core%nz
nAsyType = Core%nAsyType; nxya = Core%nxya;
nxy = Core%nxy; nxa = Core%nxa; nya = Core%nya
myzb = PE%myzb; myze = PE%myze

#ifndef MPI_ENV
ALLOCATE(Asy3DFlx(nxa,nya, ng, nz))
#else
IF(Master) THEN
  ALLOCATE(Asy3DFlx(nxa, nya, ng, nz))
ELSE
  ALLOCATE(Asy3DFlx(nxa, nya, ng, myzb:myze))
ENDIF
#endif
Do iz = myzb, myze   !Axial Sweep
  DO ig = 1, ng
    Asy3DFlx(:, :, ig, iz) = 0
    DO iasy = 1, nxya
      AsyType = Asy(iasy)%AsyType
      ixa = Asy(iasy)%ixa; iya = Asy(iasy)%iya
      nxy = AsyInfo(AsyType)%nxy
      DO i = 1, nxy  !Cell Sweep
        ixy = Asy(iasy)%GlobalPinIdx(i)
        Asy3DFlx(ixa, iya, ig, iz) = Asy3DFlx(ixa, iya, ig, iz) + Pinvol(ixy, iz) * Phic(ixy, iz, ig)
      ENDDO
    ENDDO
  ENDDO
ENDDO

#ifdef MPI_ENV
nproc = PE%nCmfdProc
COMM = PE%MPI_CMFD_COMM
IF(Master) THEN
  DO i = 1, nproc - 1
    iz1 = PE%AxDomRange(1, i); iz2 = PE%AxDomRange(2, i)
    CALL SendRecv(Asy3DFlx(:, :, :, iz1:iz2), nxa, nya, ng, iz2 - iz1 + 1,       &
                  i, FALSE, COMM)
  ENDDO
ELSE
  CALL SendRecv(Asy3DFlx(:, :, :, myzb:myze), nxa, nya, ng, myze - myzb + 1,       &
                0, TRUE, COMM)
  !ipartner, lSend, comm
ENDIF
#endif

IF(Master) THEN
  WRITE(io, '(A)') ' -  Planewise Radial Gamma Flux  -'
  WRITE(io, *)
  DO iz = 1, nz
    WRITE(io, '(5x, A5, I5)'), 'Plane', iz
    DO iya = 1, nya
      DO ig = 1, ng
        WRITE(io, '(8x, 200(1pe15.4))') (Asy3DFlx(ixa, iya, ig, iz), ixa = 1, nxa)
      ENDDO
      WRITE(io, *)
    ENDDO
  ENDDO
ENDIF
DEALLOCATE(Asy3DFlx)

NULLIFY(PhiC)
NULLIFY(AsyInfo); NULLIFY(Asy)
NULLIFY(hz); NULLIFY(PinVol)
NULLIFY(AsyVol)

END SUBROUTINE

SUBROUTINE CoreGammaSpectrum(io, core, CmInfo, ng, PE)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type,   &
                    PE_TYPE
USE IOUTIL,  ONLY : PrintReal1DarrayTo2Darray
USE GammaTYPEDEF,  ONLY : GammaCMFD_TYPE
#ifdef MPI_ENV
USE MpiComm_mod, ONLY : REDUCE
#endif
IMPLICIT NONE
INTEGER :: io, ng
TYPE(CoreInfo_Type) :: Core
TYPE(GammaCMFD_TYPE) :: CMInfo
TYPE(PE_TYPE) :: PE

REAL, POINTER :: PhiC(:, :, :)
REAL, POINTER :: PinVol(:, :)
REAL, POINTER :: hz(:)
INTEGER :: nz, nxy, nxya, nxa, nya, myzb, myze
INTEGER :: nAsyType, AsyType
INTEGER :: iz, ig, iasy, ixy, ixa, iya, i, j, k
LOGICAL :: master, slave
REAL :: Spectrum(ng)
#ifdef MPI_ENV
REAL :: buf(ng)
#endif

Master = PE%Master; Slave = .NOT. Master
myzb = PE%myzb; myze = PE%myze
PhiC => CmInfo%GPhiC;
hz => Core%hz;
PinVol => Core%PinVol
nxy = Core%nxy; nz = Core%nz

DO ig = 1, ng
  Spectrum(ig) = 0
  DO iz = myzb, myze
    DO ixy = 1, nxy
      Spectrum(ig) = Spectrum(ig) + Phic(ixy, iz, ig) * PinVol(ixy, iz)
    ENDDO
  ENDDO
ENDDO

#ifdef MPI_ENV
CALL REDUCE(Spectrum, buf, ng, PE%MPI_CMFD_COMM, FALSE)
Spectrum=buf
#endif

IF(Master) THEN
  WRITE(io, '(A)') ' - Core Spectrum -'
  WRITE(io, *)
  DO ig = 1, ng
    WRITE(io, '(5x,I5,1pe15.4)') ig, Spectrum(ig)
  ENDDO
  WRITE(io, *)
ENDIF
NULLIFY(PhiC)
NULLIFY(hz)
NULLIFY(PinVol)

END SUBROUTINE

SUBROUTINE Pin3DGammaFlux(io, core, CmInfo, ng, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,             ONLY : CoreInfo_Type,   PE_Type,      &
                                Asy_Type,        AsyInfo_Type
USE CNTL,                ONLY : nTracerCntl_Type

USE BasicOperation,      ONLY : CP_CA,           CP_VA
USE GammaTYPEDEF,       ONLY : GammaCMFD_TYPE
#ifdef MPI_ENV
USE MPIComm_Mod,         ONLY : REDUCE
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(GammaCMFD_TYPE) :: CmInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
INTEGER :: io, ng

TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
REAL, POINTER :: PhiC(:, :, :)

INTEGER, PARAMETER :: noutline = 10

INTEGER :: i, j, k, ibeg, iend
INTEGER :: ixy, ixya, ixa, iya, ix, iy, iz, ig
INTEGER :: iasytype

INTEGER :: nz, myzb, myze
INTEGER :: npinout, nstep, nout
INTEGER :: idx(4, noutline)

CHARACTER(5) :: GrpIdx
REAL, ALLOCATABLE :: OutDat(:, :, :), buf(:, :, :)


Asy => Core%Asy; AsyInfo => Core%AsyInfo
PhiC => CmInfo%GPhiC

nz = Core%nz; myzb = PE%myzb; myze = PE%myze
npinout = nTracerCntl%OutpCntl%FluxOutList(1, 0)

ALLOCATE(OutDat(ng, nz, noutline))
ALLOCATE(Buf(ng, nz, noutline))

nstep = INT(npinout / noutline)
IF(mod(npinout, noutline) .NE. 0) nstep = nstep + 1
IF(PE%MASTER) THEN
  WRITE(io, '(A)') ' - Gamma Flux Distribution for selected Pins -  '
  WRITE(io, '(7x, A)') 'Pin Identifier'
  DO i = 1, npinout
    WRITE(io, '(7x, I5, 3x, A2, 3x, (A1, 4I4, A3))') i, '=>', '[', nTracerCntl%OutpCntl%FluxOutList(1:4, i), ']'
  ENDDO
  WRITE(io, *)
ENDIF
DO j = 1, nstep
  ibeg = (j - 1) * noutline + 1; iend = j * noutline
  iend = min(iend, npinout)
  nout = iend - ibeg +1
  CALL CP_CA(OutDat, 0._8, ng, nz, noutline)
  DO i = ibeg, iend
    k  = i - ibeg + 1
    ixa = nTracerCntl%OutpCntl%FluxOutList(1, i); iya = nTracerCntl%OutpCntl%FluxOutList(2, i)
    ix = nTracerCntl%OutpCntl%FluxOutList(3, i);  iy = nTracerCntl%OutpCntl%FluxOutList(4, i)
    ixya = Core%CoreIdx(ixa, iya); iasytype = Core%CoreMap(ixya)
    ixy = AsyInfo(iasytype)%Pin2DIdx(ix, iy)
    ixy = Asy(ixya)%GlobalPinIdx(ixy)
    idx(:, k) = nTracerCntl%OutpCntl%FluxOutList(2, i)
    DO ig = 1, ng
      DO iz = myzb, myze
        OutDat(ig, iz, k)= PhiC(ixy, iz, ig)
      ENDDO
    ENDDO
  ENDDO
#ifdef MPI_ENV
  CALL REDUCE(OutDat, Buf, ng, nz, noutline, PE%MPI_CMFD_COMM, .FALSE.)
  CALL CP_VA(OutDat, buf, ng, nz, noutline)
#endif
  IF(.NOT. PE%MASTER) CYCLE
  WRITE(io, '(7x, A)') 'Pin Index ->'
  WRITE(io, '(4x, A5, 2x, A5, 2x, 300(I10,2x))') 'GROUP', 'PLANE', (i, i = ibeg, iend)
  DO ig = 1, ng
    WRITE(GrpIdx, '(I5)') ig
    DO iz = 1, nz
      IF(iz .NE. 1) GrpIdx=''
      WRITE(io, '(4x, A5, 2x, I5, 2x, 300(es10.3,2x))') GrpIdx, iz, (OutDat(ig, iz, k), k = 1, nout)
    ENDDO
  ENDDO
  WRITE(io, *)
ENDDO


NULLIFY(Asy, ASyInfo)

END SUBROUTINE


!END SUBROUTINE