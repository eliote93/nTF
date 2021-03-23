SUBROUTINE ProcessKERMA(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE)
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
REAL, POINTER :: pinOutXSkF(:, :), pinOutPGen(:, :)
REAL, POINTER :: pinOutKERMA_T(:, :),pinOutKERMA_F(:, :),pinOutKERMA_D(:, :)
REAL, POINTER :: pinOutNFlux(:), PinOutDelXSkf(:, :)
REAL, POINTER :: pinNum(:)
!Output Data
REAL, POINTER :: OutXSkF(:, :, :),OutPGen(:, :, :),OutKERMA_T(:, :, :),OutKERMA_F(:, :, :),OutKERMA_D(:, :, :)
REAL, POINTER :: OutDelXSkf(:, :, :)
REAL, POINTER :: OutNFlux(:, :)
REAL, POINTER :: Area(:), Num(:, :)
LOGICAL :: lAIC
CHARACTER(256) :: fn

IF(.NOT. nTracerCntl%lXsLib) RETURN

Pin => Core%Pin; CellInfo => Core%CellInfo;
!ResVarPin => Core%ResVarPin
ng = GroupInfo%ng

ALLOCATE(OutXSkF(200, ng, nTracerCntl%OutpCntl%nRegKERMAOut))
ALLOCATE(OutPGen(200, ng, nTracerCntl%OutpCntl%nRegKERMAOut))
ALLOCATE(OutKERMA_T(200, ng, nTracerCntl%OutpCntl%nRegKERMAOut))
ALLOCATE(OutKERMA_F(200, ng, nTracerCntl%OutpCntl%nRegKERMAOut))
ALLOCATE(OutKERMA_D(200, ng, nTracerCntl%OutpCntl%nRegKERMAOut))
ALLOCATE(OutDelXSkf(200, ng, nTracerCntl%OutpCntl%nRegKERMAOut))
ALLOCATE(OutNFlux(ng, nTracerCntl%OutpCntl%nRegKERMAOut))
ALLOCATE(Area(nTracerCntl%OutpCntl%nRegKERMAOut))
ALLOCATE(NUM(200, nTracerCntl%OutpCntl%nRegKERMAOut))

DO i = 1, nTracerCntl%OutpCntl%nRegKERMAOut
  CALL CP_CA(OutXSkF(1:nTracerCntl%OutpCntl%IsoKERMAOutList(0, i), 1:ng, i), 0._8, nTracerCntl%OutpCntl%IsoKERMAOutList(0, i), ng)
  CALL CP_CA(OutPGen(1:nTracerCntl%OutpCntl%IsoKERMAOutList(0, i), 1:ng, i), 0._8, nTracerCntl%OutpCntl%IsoKERMAOutList(0, i), ng)
  CALL CP_CA(OutKERMA_T(1:nTracerCntl%OutpCntl%IsoKERMAOutList(0, i), 1:ng, i), 0._8, nTracerCntl%OutpCntl%IsoKERMAOutList(0, i), ng)
  CALL CP_CA(OutKERMA_F(1:nTracerCntl%OutpCntl%IsoKERMAOutList(0, i), 1:ng, i), 0._8, nTracerCntl%OutpCntl%IsoKERMAOutList(0, i), ng)
  CALL CP_CA(OutKERMA_D(1:nTracerCntl%OutpCntl%IsoKERMAOutList(0, i), 1:ng, i), 0._8, nTracerCntl%OutpCntl%IsoKERMAOutList(0, i), ng)
  CALL CP_CA(OutDelXSkf(1:nTracerCntl%OutpCntl%IsoKERMAOutList(0, i), 1:ng, i), 0._8, nTracerCntl%OutpCntl%IsoKERMAOutList(0, i), ng)
  CALL CP_CA(OutNFlux(1:ng, i), 0._8, ng)
  CALL CP_CA(NUM(1:200, i), 0._8, 200)
  iz = nTracerCntl%OutpCntl%RegKERMAOutList(1, i)
  ixa = nTracerCntl%OutpCntl%RegKERMAOutList(2, i); iya = nTracerCntl%OutpCntl%RegKERMAOutList(3, i)
  ixya = Core%CoreIdx(ixa, iya); iasytype = Core%CoreMap(ixya)
  IF (nTracerCntl%OutpCntl%RegKERMAOutASM(i)) THEN
      ALLOCATE(pinOutXSkF(200, ng),pinOutPGen(200, ng),pinOutKERMA_T(200, ng),pinOutKERMA_F(200, ng),pinOutKERMA_D(200, ng),pinOutDelXSkf(200, ng),pinNum(200))
      nxc = Core%AsyInfo(iasytype)%nx; nyc = Core%AsyInfo(iasytype)%ny;
      AREA(i)=0
      DO iyc=1,nyc
         DO ixc=1,nxc
             ixyc = Core%AsyInfo(iasyType)%pin2DIdx(ixc, iyc)
             ixy = Core%Asy(ixya)%GlobalPinIdx(ixyc); icel = Core%Pin(ixy)%Cell(iz)
             ir1=1; ir2=CellInfo(icel)%nFXR
             CALL FxrKERMACnds_Numerator(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE, pinOutXSkF, &
             pinOutPGen,pinOutKERMA_T,pinOutKERMA_F,pinOutKERMA_D,pinOutDelXSkf,pinOutNFlux,pinNum,pinArea,i,ir1,ir2,ixy,iz,ng)
             OutXSkF(:,:,i)=OutXSkF(:,:,i)+pinOutXSkF
             OutPGen(:,:,i)=OutPGen(:,:,i)+pinOutPGen
             OutKERMA_T(:,:,i)=OutKERMA_T(:,:,i)+pinOutKERMA_T
             OutKERMA_F(:,:,i)=OutKERMA_F(:,:,i)+pinOutKERMA_F
             OutKERMA_D(:,:,i)=OutKERMA_D(:,:,i)+pinOutKERMA_D
             OutDelXSkf(:,:,i)=OutDelXSkf(:,:,i)+pinOutDelXSkf
             Num(:,i)=Num(:,i)+pinNum
             AREA(i)=AREA(i)+pinArea
             OutNFlux(:,i)=OutNFlux(:,i)+pinOutNFlux
         ENDDO
      ENDDO
      DEALLOCATE(pinOutXSkF,pinOutPGen,pinOutKERMA_T,pinOutKERMA_F,pinOutKERMA_D,pinOutDelXSkf,pinOutNFlux,pinNUM)
  ELSEIF (nTracerCntl%OutpCntl%RegKERMAOutPIN(i)) THEN
      ixc = nTracerCntl%OutpCntl%RegKERMAOutList(4, i); iyc = nTracerCntl%OutpCntl%RegKERMAOutList(5, i)
      if (ixc.ne.0.and.iyc.ne.0) then
          ixyc = Core%AsyInfo(iasyType)%pin2DIdx(ixc, iyc)
          ixy = Core%Asy(ixya)%GlobalPinIdx(ixyc); icel = Pin(ixy)%Cell(iz)
          ir1=1; ir2=CellInfo(icel)%nFXR
          CALL FxrKERMACnds_Numerator(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE, OutXSkF(:,:,i), &
          OutPGen(:,:,i),OutKERMA_T(:,:,i),OutKERMA_F(:,:,i),OutKERMA_D(:,:,i),OutDelXSkf(:,:,i),OutNFlux(:,i),Num(:,i),AREA(i),i,ir1,ir2,ixy,iz,ng)
      elseif (ixc.eq.0) then
          ALLOCATE(pinOutXSkF(200, ng),pinOutPGen(200, ng), pinOutKERMA_T(200, ng), pinOutKERMA_F(200, ng), pinOutKERMA_D(200, ng), pinOutDelXSkf(200, ng),pinOutNFlux(ng),pinNum(200))
          nxc=Core%AsyInfo(iasyType)%nx
          AREA(i)=0
          do ixc=1,nxc
              ixyc = Core%AsyInfo(iasyType)%pin2DIdx(ixc, iyc)
              ixy = Core%Asy(ixya)%GlobalPinIdx(ixyc); icel = Pin(ixy)%Cell(iz)
              ir1=1; ir2=CellInfo(icel)%nFXR
              CALL FxrKERMACnds_Numerator(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE, pinOutXSkF,pinOutPGen,&
                  pinOutKERMA_T,pinOutKERMA_F,pinOutKERMA_D,pinOutDelXSkf,pinOutNFlux,pinNum,pinArea,i,ir1,ir2,ixy,iz,ng)
              OutXSkF(:,:,i)=OutXSkF(:,:,i)+pinOutXSkF
              OutPGen(:,:,i)=OutPGen(:,:,i)+pinOutPGen
              OutKERMA_T(:,:,i)=OutKERMA_T(:,:,i)+pinOutKERMA_T
              OutKERMA_F(:,:,i)=OutKERMA_F(:,:,i)+pinOutKERMA_F
              OutKERMA_D(:,:,i)=OutKERMA_D(:,:,i)+pinOutKERMA_D
              OutKERMA_D(:,:,i)=OutKERMA_D(:,:,i)+pinOutDelXSkf
              OutNFlux(:,i)=OutNFlux(:,i)+pinOutNFlux
              Num(:,i)=Num(:,i)+pinNum
              AREA(i)=AREA(i)+pinArea
          enddo
          DEALLOCATE(pinOutXSkF,pinOutPGen,pinOutKERMA_T,pinOutKERMA_F,pinOutKERMA_D,pinOutDelXSkf,pinOutNFlux,pinNUM)
      else
          ALLOCATE(pinOutXSkF(200, ng),pinOutPGen(200, ng),pinOutKERMA_T(200, ng),pinOutKERMA_F(200, ng),pinOutKERMA_D(200, ng),pinOutDelXSkf(200, ng),pinOutNFlux(ng),pinNum(200))
          nyc=Core%AsyInfo(iasyType)%ny
          AREA(i)=0
          do iyc=1,nyc
              ixyc = Core%AsyInfo(iasyType)%pin2DIdx(ixc, iyc)
              ixy = Core%Asy(ixya)%GlobalPinIdx(ixyc); icel = Pin(ixy)%Cell(iz)
              ir1=1; ir2=CellInfo(icel)%nFXR
              CALL FxrKERMACnds_Numerator(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE, pinOutXSkF, pinOutPGen, &
                  pinOutKERMA_T,pinOutKERMA_F,pinOutKERMA_D,pinOutDelXSkf,pinOutNFlux,pinNum,pinArea,i,ir1,ir2,ixy,iz,ng)
              OutXSkF(:,:,i)=OutXSkF(:,:,i)+pinOutXSkF
              OutPGen(:,:,i)=OutPGen(:,:,i)+pinOutPGen
              OutKERMA_T(:,:,i)=OutKERMA_T(:,:,i)+pinOutKERMA_T
              OutKERMA_F(:,:,i)=OutKERMA_F(:,:,i)+pinOutKERMA_F
              OutKERMA_D(:,:,i)=OutKERMA_D(:,:,i)+pinOutKERMA_D
              OutDelXSkf(:,:,i)=OutDelXSkf(:,:,i)+pinOutDelXSkf
              OutNFlux(:,i)=OutNFlux(:,i)+pinOutNFlux
              Num(:,i)=Num(:,i)+pinNum
              AREA(i)=AREA(i)+pinArea
          enddo
          DEALLOCATE(pinOutXSkF,pinOutPGen,pinOutKERMA_T,pinOutKERMA_F,pinOutKERMA_D,pinOutDelXSkf,pinOutNFlux,pinNUM)
      endif
  ELSEIF (nTracerCntl%OutpCntl%RegKERMAOutFXR(i)) THEN
      ixc = nTracerCntl%OutpCntl%RegKERMAOutList(4, i); iyc = nTracerCntl%OutpCntl%RegKERMAOutList(5, i)
      ixyc = Core%AsyInfo(iasyType)%pin2DIdx(ixc, iyc)
      ixy = Core%Asy(ixya)%GlobalPinIdx(ixyc); icel = Core%Pin(ixy)%Cell(iz)
      !mypin => ResVarPin(ixy,iz); lAIC = CellInfo(icel)%lAIC
      !PinFuelTempAvgsq = dsqrt(GetPinFuelTemp(Core, FmInfo%Fxr, iz, ixy))
      !IF (mypin%lres.and.nTracerCntl%lRIF) CALL EffRIFPin(mypin, PinFuelTempAvgsq, iz, lAIC, PE)
      ir1=nTracerCntl%OutpCntl%RegKERMAOutList(6, i)
      ir2=nTracerCntl%OutpCntl%RegKERMAOutList(7, i)
      CALL FxrKERMACnds_Numerator(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE, OutXSkF(:,:,i), OutPGen(:,:,i),&
          OutKERMA_T(:,:,i),OutKERMA_F(:,:,i),OutKERMA_D(:,:,i),OutDelXSkf(:,:,i),OutNFlux(:,i),Num(:,i),AREA(i),i,ir1,ir2,ixy,iz,ng)
  ENDIF
  DO j = 1, nTracerCntl%OutpCntl%IsoKERMAOutList(0, i)
    NUM(j, i) = NUM(j, i) / Area(i)
    DO ig =  1, ng
      OutXSkF(j, ig, i) = OutXSkF(j, ig, i) / OutNFlux(ig, i)
      OutXSkF(j, ig, i) = OutXSkF(j, ig, i) / NUM(j, i)
      OutPGen(j, ig, i) = OutPGen(j, ig, i) / OutNFlux(ig, i)
      OutPGen(j, ig, i) = OutPGen(j, ig, i) / NUM(j, i)
      OutKERMA_T(j, ig, i) = OutKERMA_T(j, ig, i) / OutNFlux(ig, i)
      OutKERMA_T(j, ig, i) = OutKERMA_T(j, ig, i) / NUM(j, i)
      OutKERMA_F(j, ig, i) = OutKERMA_F(j, ig, i) / OutNFlux(ig, i)
      OutKERMA_F(j, ig, i) = OutKERMA_F(j, ig, i) / NUM(j, i)
      OutKERMA_D(j, ig, i) = OutKERMA_D(j, ig, i) / OutNFlux(ig, i)
      OutKERMA_D(j, ig, i) = OutKERMA_D(j, ig, i) / NUM(j, i)
      OutDelXSkf(j, ig, i) = OutDelXSkf(j, ig, i) / OutNFlux(ig, i)
      OutDelXSkf(j, ig, i) = OutDelXSkf(j, ig, i) / NUM(j, i)
    ENDDO
  ENDDO
  DO ig = 1, ng
    OutNFlux(ig, i) = OutNFlux(ig, i) / Area(i)
  END DO
ENDDO
! Photoatomic XS Tally
fn=trim(caseid)//'.kerma'

CALL openfile(io16,FALSE,FALSE,FALSE, fn)
DO i = 1, nTracerCntl%OutpCntl%nRegKERMAOut
  IF (nTracerCntl%OutpCntl%RegKERMAOutASM(i)) THEN
      WRITE(io16, '(A10, I5, A3, 1(I5, I5, A3),200I7)') 'Region :',  nTracerCntl%OutpCntl%RegKERMAOutList(1, i), '/', &
                                             nTracerCntl%OutpCntl%RegKERMAOutList(2:3, i), '/', (nTracerCntl%OutpCntl%IsoKERMAOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoKERMAOutList(0, i))
  ELSEIF (nTracerCntl%OutpCntl%RegKERMAOutPIN(i)) THEN
      WRITE(io16, '(A10, I5, A3, 2(I5, I5, A3),200I7)') 'Region :',  nTracerCntl%OutpCntl%RegKERMAOutList(1, i), '/', &
                                             nTracerCntl%OutpCntl%RegKERMAOutList(2:3, i), '/', &
                                             nTracerCntl%OutpCntl%RegKERMAOutList(4:5, i), '/', (nTracerCntl%OutpCntl%IsoKERMAOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoKERMAOutList(0, i))
  ELSEIF (nTracerCntl%OutpCntl%RegKERMAOutFXR(i)) THEN
      WRITE(io16, '(A10, I5, A3, 3(I5, I5, A3),200I7)') 'Region :',  nTracerCntl%OutpCntl%RegKERMAOutList(1, i), '/', &
                                             nTracerCntl%OutpCntl%RegKERMAOutList(2:3, i), '/', &
                                             nTracerCntl%OutpCntl%RegKERMAOutList(4:5, i), '/', &
                                             nTracerCntl%OutpCntl%RegKERMAOutList(6:7, i), '/', (nTracerCntl%OutpCntl%IsoKERMAOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoKERMAOutList(0, i))
  ENDIF
  WRITE(io16, '(A5, 2A16, 600I16)') 'GRP', 'PFlux', 'Area', &
      (nTracerCntl%OutpCntl%IsoKERMAOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoKERMAOutList(0, i)),&
      (nTracerCntl%OutpCntl%IsoKERMAOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoKERMAOutList(0, i)),&
      (nTracerCntl%OutpCntl%IsoKERMAOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoKERMAOutList(0, i))

  DO ig = 1, ng, 1
    WRITE(io16, '(I5, 600(1pe16.6))') ig, OutNFlux(ig, i), Area(i),  &
                                (OutKERMA_T(j, ig, i), j = 1, nTracerCntl%OutpCntl%IsoKERMAOutList(0, i)), &
                                (OutXSkF(j, ig, i), j = 1, nTracerCntl%OutpCntl%IsoKERMAOutList(0, i)), &
                                (OutPGen(j, ig, i), j = 1, nTracerCntl%OutpCntl%IsoKERMAOutList(0, i)), &
                                (OutKERMA_F(j, ig, i), j = 1, nTracerCntl%OutpCntl%IsoKERMAOutList(0, i)), &
                                (OutKERMA_D(j, ig, i), j = 1, nTracerCntl%OutpCntl%IsoKERMAOutList(0, i)), &
                                (OutDelXSkf(j, ig, i), j = 1, nTracerCntl%OutpCntl%IsoKERMAOutList(0, i))

  ENDDO
  DO j = 1, (37+16*nTracerCntl%OutpCntl%IsoKERMAOutList(0, i))
    WRITE(io16, '(A)', advance='NO') '-'
  ENDDO
  WRITE(io16, *)
  WRITE(io16, '(A35, 50(1pe15.5))') 'N.D :', (NUM(j, i), j = 1, nTracerCntl%OutpCntl%IsoKERMAOutList(0, i))
  WRITE(io16, *)
ENDDO
close(io16)

CLOSE(io16)

DEALLOCATE(OutXSkF)
DEALLOCATE(OutPGen)
DEALLOCATE(OutKERMA_T)
DEALLOCATE(OutNFlux)
DEALLOCATE(Area)
DEALLOCATE(NUM)
END SUBROUTINE
  
  
SUBROUTINE FxrKERMACnds_Numerator(Core, FmInfo, ThInfo, GroupInfo,nTracerCntl, PE, OutXSkF,OutPGen, OutKERMA_T,OutKERMA_F,OutKERMA_D,OutDelXSkf,OutFlux,Num,AREA,i,ir1,ir2,ixy,iz,ng)
    USE PARAM, only : CKELVIN, TRUE, FALSE
    USE TYPEDEF,          ONLY : XsMac_TYPE, CoreInfo_Type,FmInfo_Type, ThInfo_Type,GroupInfo_Type,FxrInfo_Type,PE_Type
    USE GammaTYPEDEF,     ONLY : GamMacXS_TYPE
    USE CNTL,             ONLY : nTracerCntl_Type
    USE MacXsLib_Mod,     ONLY : MacXsBase
    USE BenchXs,          ONLY : XsBaseBen
    USE GamXsLib_Mod,     ONLY : GetLocalQMat
    USE MOC_Mod,          ONLY : FxrAvgPhi
    USE XsUtil_mod,       ONLY : GetXsMacDat,         ReturnXsMacDat
    USE GammaCore_mod,    ONLY : Gphis
    USE ioutil
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FmInfo_Type) :: FmInfo
    TYPE(ThInfo_Type) :: ThInfo
    TYPE(GroupInfo_Type) :: GroupInfo
    TYPE(FxrInfo_Type), POINTER :: myFxr
    TYPE(XsMac_TYPE), POINTER :: XsMac
    TYPE(nTracerCntl_Type) :: nTracerCntl
    TYPE(PE_Type) :: PE

    INTEGER, INTENT(IN) :: i,ir1,ir2,ixy,iz,ng
    !TYPE(ResVarPin_Type), INTENT(IN) :: mypin
    !REAL,INTENT(IN) :: PinFuelTempAvgsq
    REAL,INTENT(OUT) :: AREA
    REAL,INTENT(OUT) :: OutXSkF(200, ng),OutPGen(200, ng),OutKERMA_T(200, ng),OutKERMA_F(200, ng),OutKERMA_D(200, ng), OutDelXSkf(200,ng)
    REAL,INTENT(OUT) :: OutFlux(ng)
    REAL,INTENT(OUT) :: Num(200)
    REAL :: Flux(ng)

    LOGICAL :: lXsLib

    INTEGER :: ir, j, k, id, ig, iResoGrpBeg, iResoGrpEnd, ifxr, jG, iso
    !Fxr Data
    INTEGER :: niso
    LOGICAL :: lres
    INTEGER, POINTER :: idiso(:)
    REAL, POINTER :: pnum(:)
    REAL :: eigv
    REAL :: Temp
    
    INTEGER :: ifsrlocal
    
    lXsLib = nTracerCntl%lXsLib
    iResoGrpBeg = GroupInfo%nofg + 1; iResoGrpEnd = GroupInfo%nofg + GroupInfo%norg
    CALL GetXsMacDat(XsMac, ng, .TRUE.)
    AREA=0; OutFlux=0; OutXSkF=0; OutPGen=0; OutKERMA_T=0; NUM=0;
    DO ir= ir1, ir2
        ifxr = Core%Pin(ixy)%FxrIdxst + ir - 1;myFxr => FmInfo%Fxr(ifxr, iz)
        temp = myFxr%temp
        niso = myFxr%niso; idiso => myFxr%idiso
        pnum => myFxr%pnum; lres = myFxr%lres
          CALL MacXsBase(XsMac, myFxr, 1, ng, ng, 1._8, FALSE, TRUE, TRUE)  ! CrCSPFtn On        
          CALL GetLocalQMat(XSMac, myFxr, 1, ng, ng, TRUE)
          !Self-Sheilding Effect
          IF(myFxr%lres) THEN
             DO ig = iResoGrpBeg, iResoGrpEnd
               DO iso = 1, niso
                  XsMac%IsoXsMackf(iso,ig) = XsMac%IsoXsMackf(iso,ig) * myFxr%fresoFIso(iso,ig)  
               END DO
             END DO
          ENDIF
          !Obtaining
        AREA = AREA + myFxr%area
        Flux = FxrAvgPhi(Core, FmInfo%Fxr, FMInfo%Phis, ixy, ir, iz, ng, PE)
        Outflux = OutFlux + myFxr%area * Flux

        DO j = 1, nTracerCntl%OutpCntl%IsoKERMAOutList(0, i)
          id = nTracerCntl%OutpCntl%IsoKERMAOutList(j, i)
          DO k = 1, niso
            IF(IDISO(k) .NE. id) CYCLE
            NUM(j) = NUM(j) +  myFxr%area * pnum(k)
            DO ig =  1, ng
              OutKERMA_T(j, ig) = OutKERMA_T(j, ig) +  myFxr%area * XsMac%IsoMacKERMA_t(k, ig) * Flux(ig)
              OutKERMA_F(j, ig) = OutKERMA_F(j, ig) +  myFxr%area * XsMac%IsoMacKERMA_F(k, ig) * Flux(ig)
              OutKERMA_D(j, ig) = OutKERMA_D(j, ig) +  myFxr%area * XsMac%IsoMacKERMA_D(k, ig) * Flux(ig)
              OutDelXSkf(j, ig) = OutDelXSkf(j, ig) +  myFxr%area * XsMac%IsoMacDelkf(k, ig) * Flux(ig)
              OutXskF(j, ig)    = OutXskF(j, ig)    +  myFxr%area *    XsMac%IsoXsMackf(k, ig) * Flux(ig)
              OutPGen(j, ig)    = OutPGen(j, ig)    +  myFxr%area * XsMac%IsoMacKERMA_p(k, ig) * Flux(ig)
            ENDDO
          ENDDO
        ENDDO
    ENDDO
    NULLIFY(IDISO, pnum, myFxr)
    CALL ReturnXsMacDat(XsMac)

  END SUBROUTINE