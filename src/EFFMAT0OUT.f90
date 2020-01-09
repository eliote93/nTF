#include <defines.h>

SUBROUTINE ProcessEffMAT0(Core, FmInfo, GroupInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF
USE CNTL,             ONLY : nTracerCntl_Type
USE FILES,            ONLY : io16,                caseid
USE XsUtil_mod,       ONLY : GetXsMacDat,         ReturnXsMacDat
USE MacXsLib_mod,     ONLY : IsoMacXsScatMatrix_Gen
USE ioutil
USE BasicOperation
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl

TYPE(Pin_Type), POINTER :: PIN(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(PE_Type) :: PE

INTEGER :: iz, ixa, iya, ixc, iyc,nxc, nyc
INTEGER :: ixy, ixya, ixyc, iasytype, icel
INTEGER :: i, j, k, ig, ig2, ir1, ir2, igb, ige
!Fxr Data
INTEGER :: ng

REAL :: pinArea
REAL, POINTER :: pinOutMAT(:, :, :)
REAL, POINTER :: pinOutFlux(:)
REAL, POINTER :: pinNum(:)
!Output Data
REAL, POINTER :: OutMAT0(:, :, :, :)
REAL, POINTER :: OutFlux(:, :)
REAL, POINTER :: Area(:), Num(:, :)
CHARACTER(256) :: fn
#ifdef __GFORTRAN__
CHARACTER(2) :: str2num
#endif

IF(.NOT. nTracerCntl%lXsLib) RETURN

Pin => Core%Pin; CellInfo => Core%CellInfo; Fxr => FmInfo%Fxr
ng = GroupInfo%ng
ALLOCATE(OutMAT0(200, ng, ng, nTracerCntl%OutpCntl%nRegMAT0Out))
ALLOCATE(OutFlux(ng, nTracerCntl%OutpCntl%nRegMAT0Out))
ALLOCATE(Area(nTracerCntl%OutpCntl%nRegMAT0Out))
ALLOCATE(NUM(200, nTracerCntl%OutpCntl%nRegMAT0Out))
DO i = 1, nTracerCntl%OutpCntl%nRegMAT0Out
  OutMAT0(:,:,:,i)=0._8
  CALL CP_CA(OutFLux(1:ng, i), 0._8, ng)
  CALL CP_CA(NUM(1:200, i), 0._8, 200)
  iz = nTracerCntl%OutpCntl%RegMAT0OutList(1, i)
  ixa = nTracerCntl%OutpCntl%RegMAT0OutList(2, i); iya = nTracerCntl%OutpCntl%RegMAT0OutList(3, i)
  ixya = Core%CoreIdx(ixa, iya); iasytype = Core%CoreMap(ixya)
  IF (nTracerCntl%OutpCntl%RegMAT0OutASM(i)) THEN
      ALLOCATE(pinOutMAT(200,ng,ng),pinOutFlux(ng),pinNum(200))
      nxc = Core%AsyInfo(iasytype)%nx; nyc = Core%AsyInfo(iasytype)%ny;
      AREA(i)=0
      DO iyc=1,nyc
         DO ixc=1,nxc
             ixyc = Core%AsyInfo(iasyType)%pin2DIdx(ixc, iyc)
             ixy = Core%Asy(ixya)%GlobalPinIdx(ixyc); icel = Core%Pin(ixy)%Cell(iz)
             ir1=1; ir2=CellInfo(icel)%nFXR
             CALL FxrMAT0Cnds_Numerator(Core, FmInfo, GroupInfo, nTracerCntl, PE, pinOutMAT,pinOutFlux,pinNum,pinArea,i,ir1,ir2,ixy,iz,ng)
             OutMAT0(:,:,:,i)=OutMAT0(:,:,:,i)+pinOutMAT
             OutFlux(:,i)=OutFlux(:,i)+pinOutFlux
             Num(:,i)=Num(:,i)+pinNum
             AREA(i)=AREA(i)+pinArea
         ENDDO
      ENDDO
      DEALLOCATE(pinOutMAT,pinOutFlux,pinNUM)
  ELSEIF (nTracerCntl%OutpCntl%RegMAT0OutPIN(i)) THEN
      ixc = nTracerCntl%OutpCntl%RegMAT0OutList(4, i); iyc = nTracerCntl%OutpCntl%RegMAT0OutList(5, i)
      if (ixc.ne.0.and.iyc.ne.0) then
          ixyc = Core%AsyInfo(iasyType)%pin2DIdx(ixc, iyc)
          ixy = Core%Asy(ixya)%GlobalPinIdx(ixyc); icel = Pin(ixy)%Cell(iz)
          ir1=1; ir2=CellInfo(icel)%nFXR
          CALL FxrMAT0Cnds_Numerator(Core, FmInfo, GroupInfo, nTracerCntl, PE, OutMAT0(:,:,:,i),OutFlux(:,i),Num(:,i),AREA(i),i,ir1,ir2,ixy,iz,ng)
      elseif (ixc.eq.0) then
          ALLOCATE(pinOutMAT(200, ng, ng),pinOutFlux(ng),pinNum(200))
          nxc=Core%AsyInfo(iasyType)%nx
          AREA(i)=0
          do ixc=1,nxc
              ixyc = Core%AsyInfo(iasyType)%pin2DIdx(ixc, iyc)
              ixy = Core%Asy(ixya)%GlobalPinIdx(ixyc); icel = Pin(ixy)%Cell(iz)
              ir1=1; ir2=CellInfo(icel)%nFXR
              CALL FxrMAT0Cnds_Numerator(Core, FmInfo, GroupInfo, nTracerCntl, PE, pinOutMAT,pinOutFlux,pinNum,pinArea,i,ir1,ir2,ixy,iz,ng)
              OutMAT0(:,:,:,i)=OutMAT0(:,:,:,i)+pinOutMAT
              OutFlux(:,i)=OutFlux(:,i)+pinOutFlux
              Num(:,i)=Num(:,i)+pinNum
              AREA(i)=AREA(i)+pinArea
          enddo
          DEALLOCATE(pinOutMAT,pinOutFlux,pinNUM)
      else
          ALLOCATE(pinOutMAT(200, ng, ng),pinOutFlux(ng),pinNum(200))
          nyc=Core%AsyInfo(iasyType)%ny
          AREA(i)=0
          do iyc=1,nyc
              ixyc = Core%AsyInfo(iasyType)%pin2DIdx(ixc, iyc)
              ixy = Core%Asy(ixya)%GlobalPinIdx(ixyc); icel = Pin(ixy)%Cell(iz)
              ir1=1; ir2=CellInfo(icel)%nFXR
              CALL FxrMAT0Cnds_Numerator(Core, FmInfo, GroupInfo, nTracerCntl, PE, pinOutMAT,pinOutFlux,pinNum,pinArea,i,ir1,ir2,ixy,iz,ng)
              OutMAT0(:,:,:,i)=OutMAT0(:,:,:,i)+pinOutMAT
              OutFlux(:,i)=OutFlux(:,i)+pinOutFlux
              Num(:,i)=Num(:,i)+pinNum
              AREA(i)=AREA(i)+pinArea
          enddo
          DEALLOCATE(pinOutMAT,pinOutFlux,pinNUM)
      endif
  ELSEIF (nTracerCntl%OutpCntl%RegMAT0OutFXR(i)) THEN
      ixc = nTracerCntl%OutpCntl%RegMAT0OutList(4, i); iyc = nTracerCntl%OutpCntl%RegMAT0OutList(5, i)
      ixyc = Core%AsyInfo(iasyType)%pin2DIdx(ixc, iyc)
      ixy = Core%Asy(ixya)%GlobalPinIdx(ixyc); icel = Core%Pin(ixy)%Cell(iz)
      ir1=nTracerCntl%OutpCntl%RegMAT0OutList(6, i)
      ir2=nTracerCntl%OutpCntl%RegMAT0OutList(7, i)
      CALL FxrMAT0Cnds_Numerator(Core, FmInfo, GroupInfo, nTracerCntl, PE, OutMAT0(:,:,:,i),OutFlux(:,i),Num(:,i),AREA(i),i,ir1,ir2,ixy,iz,ng)
  ENDIF
  DO j = 1, nTracerCntl%OutpCntl%IsoMAT0OutList(0, i)
    NUM(j, i) = NUM(j, i) / Area(i)
    DO ig =  1, ng
        DO ig2 = 1, ng
            OutMAT0(j, ig, ig2, i) = OutMAT0(j, ig, ig2, i) / Outflux(ig, i) / NUM(j, i)
        ENDDO
    ENDDO
  ENDDO
  DO ig= 1, ng
    Outflux(ig, i) = Outflux(ig, i) / Area(i)
  ENDDO
ENDDO
fn=trim(caseid)//'.emat0'

CALL openfile(io16,FALSE,FALSE,FALSE, fn)
DO i = 1, nTracerCntl%OutpCntl%nRegMAT0Out
  IF (nTracerCntl%OutpCntl%RegMAT0OutASM(i)) THEN
      WRITE(io16, '(A10, I5, A3, 1(I5, I5, A3),200I7)') 'Region :',  nTracerCntl%OutpCntl%RegMAT0OutList(1, i), '/', &
                                             nTracerCntl%OutpCntl%RegMAT0OutList(2:3, i), '/', (nTracerCntl%OutpCntl%IsoMAT0OutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoMAT0OutList(0, i))
  ELSEIF (nTracerCntl%OutpCntl%RegMAT0OutPIN(i)) THEN
      WRITE(io16, '(A10, I5, A3, 2(I5, I5, A3),200I7)') 'Region :',  nTracerCntl%OutpCntl%RegMAT0OutList(1, i), '/', &
                                             nTracerCntl%OutpCntl%RegMAT0OutList(2:3, i), '/', &
                                             nTracerCntl%OutpCntl%RegMAT0OutList(4:5, i), '/', (nTracerCntl%OutpCntl%IsoMAT0OutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoMAT0OutList(0, i))
  ELSEIF (nTracerCntl%OutpCntl%RegMAT0OutFXR(i)) THEN
      WRITE(io16, '(A10, I5, A3, 3(I5, I5, A3),200I7)') 'Region :',  nTracerCntl%OutpCntl%RegMAT0OutList(1, i), '/', &
                                             nTracerCntl%OutpCntl%RegMAT0OutList(2:3, i), '/', &
                                             nTracerCntl%OutpCntl%RegMAT0OutList(4:5, i), '/', &
                                             nTracerCntl%OutpCntl%RegMAT0OutList(6:7, i), '/', (nTracerCntl%OutpCntl%IsoMAT0OutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoMAT0OutList(0, i))
  ENDIF
  DO j = 1, nTracerCntl%OutpCntl%IsoMAT0OutList(0, i)
    WRITE(io16, '(3A5,I16)') 'GRP', 'G1','G2',nTracerCntl%OutpCntl%IsoMAT0OutList(j, i)
    DO ig = 1, ng
      DO ig2 = ng, 1, -1
        IF (OutMAT0(j, ig2, ig, i).NE.0.) EXIT
      ENDDO
      ige = ig2
      DO ig2 = 1, ng
        IF (OutMAT0(j, ig2, ig, i).NE.0.) EXIT
      ENDDO
      igb = ig2
#ifndef __GFORTRAN__
      WRITE(io16, '(3I5, <ng>(1pe16.6))') ig, igb, ige, OutMAT0(j, igb:ige, ig, i)
#else
      READ(str2num, *) ng
      WRITE(io16, '(3I5, ' // TRIM(str2num) // '(1pe16.6))') ig, igb, ige, OutMAT0(j, igb:ige, ig, i)
#endif
    ENDDO
  ENDDO
  DO j = 1, (35+15*nTracerCntl%OutpCntl%IsoMAT0OutList(0, i))
    WRITE(io16, '(A)', advance='NO') '-'
  ENDDO
  WRITE(io16, *)
  WRITE(io16, '(A35, 20(1pe15.5))') 'N.D :', (NUM(j, i), j = 1, nTracerCntl%OutpCntl%IsoMAT0OutList(0, i))
  WRITE(io16, *)
ENDDO

close(io16)

DEALLOCATE(OutMAT0)
DEALLOCATE(OutFlux)
DEALLOCATE(Area)
DEALLOCATE(NUM)

END SUBROUTINE

SUBROUTINE FxrMAT0Cnds_Numerator(Core, FmInfo, GroupInfo,nTracerCntl, PE, OutMAT,OutFlux,Num,AREA,i,ir1,ir2,ixy,iz,ng)
    USE PARAM, only : CKELVIN
    USE TYPEDEF
    USE CNTL,             ONLY : nTracerCntl_Type
    USE MacXsLib_mod,     ONLY : IsoMacXsScatMatrix_Gen,IsoMacXsScatMatrix_GenR
    USE MOC_Mod,          ONLY : FxrAvgPhi
    USE ioutil
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FmInfo_Type) :: FmInfo
    TYPE(GroupInfo_Type) :: GroupInfo
    TYPE(FxrInfo_Type), POINTER :: myFxr
    TYPE(nTracerCntl_Type) :: nTracerCntl
    TYPE(PE_Type) :: PE

    INTEGER, INTENT(IN) :: i,ir1,ir2,ixy,iz,ng
    REAL,INTENT(OUT) :: AREA
    REAL,INTENT(OUT) :: OutMAT(200, ng, ng)
    REAL,INTENT(OUT) :: OutFlux(ng)
    REAL,INTENT(OUT) :: Num(200)
    REAL :: Flux(ng)
    REAL,POINTER :: IsoXsmacsm(:,:,:)

    INTEGER :: ir, j, k, id, ig, ig2, ifxr
    !Fxr Data
    INTEGER :: niso
    INTEGER, POINTER :: idiso(:)
    REAL, POINTER :: pnum(:)
    REAL :: Temp

    AREA=0; OutFlux=0; OutMAT=0; NUM=0
    DO ir= ir1, ir2
        ifxr = Core%Pin(ixy)%FxrIdxst + ir - 1; myFxr => FmInfo%Fxr(ifxr, iz)
        temp = myFxr%temp;    niso = myFxr%niso
        ALLOCATE(IsoXsmacsm(niso,ng,ng))
        idiso => myFxr%idiso; pnum => myFxr%pnum
        AREA = AREA + myFxr%area
        if (nTracerCntl%lRST) then
          CALL IsoMacXsScatMatrix_GenR(IsoXsMacSm, temp, niso, idiso, pnum, 1, ng, ng, GroupInfo, &
                                .FALSE., myFxr%lres, myFxr%fresoSIso, myFxr%fresoSSIso, myFxr%fresoS1Iso)
        else
          CALL IsoMacXsScatMatrix_Gen(IsoXsMacSm, temp, niso, idiso, pnum, 1, ng, ng, GroupInfo, .FALSE.)  
        endif
        Flux = FxrAvgPhi(Core, FmInfo%Fxr, FmInfo%Phis, ixy, ir, iz, ng, PE)
        Outflux = OutFlux + myFxr%area * Flux
        DO j = 1, nTracerCntl%OutpCntl%IsoMAT0OutList(0, i)
          id = nTracerCntl%OutpCntl%IsoMAT0OutList(j, i)
          DO k = 1, niso
            IF(IDISO(k) .NE. id) CYCLE
            NUM(j) = NUM(j) +  myFxr%area * pnum(k)
            DO ig =  1, ng
                DO ig2 = 1, ng
                    OutMAT(j, ig, ig2) = OutMAT(j, ig, ig2)+  myFxr%area * IsoXsMacSm(k, ig, ig2) * Flux(ig)
                ENDDO
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(IsoXsmacsm)
    ENDDO
    NULLIFY(IDISO, pnum, myFxr)

END SUBROUTINE
