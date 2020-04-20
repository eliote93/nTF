#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
!                                     01. HEX PROCESS : Eff Xs
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexProcessEffXs(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE)

USE PARAM
USE TYPEDEF
USE ioutil
USE BasicOperation

USE CNTL,    ONLY : nTracerCntl_Type, OutpCntl_Type
USE FILES,   ONLY : io16, caseid
USE HexUtil, ONLY : HexChkInc_REAL
USE HexData, ONLY : hCore, hAsyTypInfo, hAsy

IMPLICIT NONE

TYPE(CoreInfo_Type)    :: Core
TYPE(FmInfo_Type)      :: FmInfo
TYPE(ThInfo_Type)      :: ThInfo
TYPE(GroupInfo_Type)   :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type)          :: PE

TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Pin_Type),  POINTER :: PIN(:)

INTEGER :: iz, ixa, iya, ixc, iyc, nxc, nyc
INTEGER :: ixy, ixya, ixyc, iaTyp, icel
INTEGER :: i, j, ig, ng, ir1, ir2
INTEGER :: iPin, iaGeo, nPin

REAL :: pinArea
REAL, POINTER :: pinOutXSA(:, :), pinOutXSS(:, :), pinOutXSStr(:, :)
REAL, POINTER :: pinOutXSF(:, :), pinOutXSNF(:, :)
REAL, POINTER :: pinOutFlux(:)
REAL, POINTER :: pinNum(:)

!Output Data
REAL, POINTER :: OutXSA(:, :, :),OutXSS(:, :, :),OutXSStr(:, :, :),OutXSF(:, :, :),OutXSNF(:, :, :)
REAL, POINTER :: OutFlux(:, :)
REAL, POINTER :: Area(:), Num(:, :)
LOGICAL :: lAIC
CHARACTER(256) :: fn, c
! ----------------------------------------------------

IF(.NOT. nTracerCntl%lXsLib) RETURN

Pin      => Core%Pin
CellInfo => Core%CellInfo
ng        = GroupInfo%ng

ALLOCATE(OutXSA   (200, ng, nTracerCntl%OutpCntl%nRegXsOut))
ALLOCATE(OutXSS   (200, ng, nTracerCntl%OutpCntl%nRegXsOut))
ALLOCATE(OutXSStr (200, ng, nTracerCntl%OutpCntl%nRegXsOut))
ALLOCATE(OutXSF   (200, ng, nTracerCntl%OutpCntl%nRegXsOut))
ALLOCATE(OutXSNF  (200, ng, nTracerCntl%OutpCntl%nRegXsOut))
ALLOCATE(OutFlux       (ng, nTracerCntl%OutpCntl%nRegXsOut))
ALLOCATE(Area              (nTracerCntl%OutpCntl%nRegXsOut))
ALLOCATE(NUM      (200,     nTracerCntl%OutpCntl%nRegXsOut))
! ----------------------------------------------------
!               01. SET : Eff XS
! ----------------------------------------------------
DO i = 1, nTracerCntl%OutpCntl%nRegXsOut
  CALL CP_CA(OutXsA  (1:nTracerCntl%OutpCntl%IsoXsOutList(0, i), 1:ng, i), 0._8, nTracerCntl%OutpCntl%IsoXsOutList(0, i), ng)
  CALL CP_CA(OutXsS  (1:nTracerCntl%OutpCntl%IsoXsOutList(0, i), 1:ng, i), 0._8, nTracerCntl%OutpCntl%IsoXsOutList(0, i), ng)
  CALL CP_CA(OutXsStr(1:nTracerCntl%OutpCntl%IsoXsOutList(0, i), 1:ng, i), 0._8, nTracerCntl%OutpCntl%IsoXsOutList(0, i), ng)
  CALL CP_CA(OutXsF  (1:nTracerCntl%OutpCntl%IsoXsOutList(0, i), 1:ng, i), 0._8, nTracerCntl%OutpCntl%IsoXsOutList(0, i), ng)
  CALL CP_CA(OutXsNF (1:nTracerCntl%OutpCntl%IsoXsOutList(0, i), 1:ng, i), 0._8, nTracerCntl%OutpCntl%IsoXsOutList(0, i), ng)
  
  CALL CP_CA(OutFlux(1:ng,  i), 0._8, ng)
  CALL CP_CA(NUM    (1:200, i), 0._8, 200)
  
  iz   = nTracerCntl%OutpCntl%RegXsOutList(1, i)
  ixa  = nTracerCntl%OutpCntl%RegXsOutList(2, i)
  iya  = nTracerCntl%OutpCntl%RegXsOutList(3, i)
  ixya = hCore(ixa, iya)
  
  iaTyp = hAsy(ixya)%AsyTyp
  iaGeo = hAsy(ixya)%GeoTyp
  ! ----------------------------
  !      1. CASE : Asy
  ! ----------------------------
  IF (nTracerCntl%OutpCntl%RegXsOutASM(i)) THEN
    ALLOCATE (pinOutXSA   (200, ng))
    ALLOCATE (pinOutXSS   (200, ng))
    ALLOCATE (pinOutXSStr (200, ng))
    ALLOCATE (pinOutXSF   (200, ng))
    ALLOCATE (pinOutXSNF  (200, ng))
    ALLOCATE (pinOutFlux       (ng))
    ALLOCATE (pinNum      (200))
    
    nPin    = hAsy(ixya)%nTotPin
    AREA(i) = 0
    
    DO iPin = 1, nPin
      ixy  = hAsy(ixya)%PinIdxSt + iPin - 1
      iCel = Core%Pin(ixy)%Cell(iz)
      ir1  = 1
      ir2  = CellInfo(icel)%nFXR
      
      CALL HexFxrXSCnds_Numerator(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE, pinOutXSA, &
           pinOutXSS,pinOutXSStr,pinOutXSF,pinOutXSNF,pinOutFlux,pinNum,pinArea,i,ir1,ir2,ixy,iz,ng)
      
      OutXsA  (:,:,i) = OutXsA  (:,:,i) + pinOutXSA
      OutXsS  (:,:,i) = OutXsS  (:,:,i) + pinOutXSS
      OutXsStr(:,:,i) = OutXsStr(:,:,i) + pinOutXSStr
      OutXsF  (:,:,i) = OutXsF  (:,:,i) + pinOutXSF
      OutXsNF (:,:,i) = OutXsNF (:,:,i) + pinOutXSNF
      OutFlux   (:,i) = OutFlux   (:,i) + pinOutFlux
      
      Num(:,i) = Num(:,i) + pinNum
      AREA (i) = AREA (i) + pinArea
    ENDDO
    
    DEALLOCATE(pinOutXSA,pinOutXSS,pinOutXSStr,pinOutXSF,pinOutXSNF,pinOutFlux,pinNUM)
  ! ----------------------------
  !      1. CASE : Pin
  ! ----------------------------
  ELSE IF (nTracerCntl%OutpCntl%RegXsOutPIN(i)) THEN
    ixc = nTracerCntl%OutpCntl%RegXsOutList(4, i)
    iyc = nTracerCntl%OutpCntl%RegXsOutList(5, i)
    
    ixyc = hAsyTypInfo(iaTyp)%Pin2Dto1Dmap(ixc, iyc)
    ixyc = hAsyTypInfo(iaTyp)%PinLocIdx(iaGeo, ixyc)
    ixy  = hAsy(ixya)%PinIdxSt + ixyc - 1
    
    iCel = Pin(ixy)%Cell(iz)
    ir1  = 1
    ir2  = CellInfo(icel)%nFXR
    
    CALL HexFxrXSCnds_Numerator(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE, OutXSA(:,:,i), &
        OutXSS(:,:,i), OutXSStr(:,:,i),OutXSF(:,:,i),OutXSNF(:,:,i),OutFlux(:,i),Num(:,i),AREA(i),i,ir1,ir2,ixy,iz,ng)
  ! ----------------------------
  !      1. CASE : FXR
  ! ----------------------------
  ELSE IF (nTracerCntl%OutpCntl%RegXsOutFXR(i)) THEN
    ixc = nTracerCntl%OutpCntl%RegXsOutList(4, i)
    iyc = nTracerCntl%OutpCntl%RegXsOutList(5, i)
    
    ixyc = hAsyTypInfo(iaTyp)%Pin2Dto1Dmap(ixc, iyc)
    ixyc = hAsyTypInfo(iaTyp)%PinLocIdx(iaGeo, ixyc)
    ixy  = hAsy(ixya)%PinIdxSt + ixyc - 1
    
    iCel = Core%Pin(ixy)%Cell(iz)
    ir1  = nTracerCntl%OutpCntl%RegXsOutList(6, i)
    ir2  = nTracerCntl%OutpCntl%RegXsOutList(7, i)
    
    CALL HexFxrXSCnds_Numerator(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE, OutXSA(:,:,i), OutXSS(:,:,i),&
        OutXSStr(:,:,i),OutXSF(:,:,i),OutXSNF(:,:,i),OutFlux(:,i),Num(:,i),AREA(i),i,ir1,ir2,ixy,iz,ng)
  ENDIF
  ! ----------------------------
  !      2. XS / Flx / #
  ! ----------------------------
  DO j = 1, nTracerCntl%OutpCntl%IsoXsOutList(0, i)
    CALL HexChkInc_REAL(0._8, Num(j, i), "EFF XS : ISOTOPE DOES NOT EXIST IN FXR")
    
    NUM(j, i) = NUM(j, i) / Area(i)
    
    DO ig = 1, ng
      OutXsA  (j, ig, i) = OutXsA  (j, ig, i) / Outflux(ig, i)
      OutXsA  (j, ig, i) = OutXsA  (j, ig, i) / NUM(j, i)
      OutXsS  (j, ig, i) = OutXsS  (j, ig, i) / Outflux(ig, i)
      OutXsS  (j, ig, i) = OutXsS  (j, ig, i) / NUM(j, i)
      OutXsStr(j, ig, i) = OutXsStr(j, ig, i) / Outflux(ig, i)
      OutXsStr(j, ig, i) = OutXsStr(j, ig, i) / NUM(j, i)
      OutXsF  (j, ig, i) = OutXsF  (j, ig, i) / Outflux(ig, i)
      OutXsF  (j, ig, i) = OutXsF  (j, ig, i) / NUM(j, i)
      OutXsNF (j, ig, i) = OutXsNF (j, ig, i) / Outflux(ig, i)
      OutXsNF (j, ig, i) = OutXsNF (j, ig, i) / NUM(j, i)
    END DO
  END DO
  ! ----------------------------
  !      3. Flx / Area
  ! ----------------------------
  DO ig = 1, ng
    Outflux(ig, i) = Outflux(ig, i) / Area(i)
  END DO
END DO
! ----------------------------------------------------
!               02. PRINT
! ----------------------------------------------------
fn = trim(caseid) // '.exs'

CALL openfile(io16,FALSE,FALSE,FALSE, fn)

DO i = 1, nTracerCntl%OutpCntl%nRegXsOut
  IF (nTracerCntl%OutpCntl%RegXsOutASM(i)) THEN
      WRITE(io16, '(A10, I5, A3, 1(I5, I5, A3),200I7)') 'Region :',  nTracerCntl%OutpCntl%RegXsOutList(1, i), '/', &
                                             nTracerCntl%OutpCntl%RegXsOutList(2:3, i), '/', (nTracerCntl%OutpCntl%IsoXsOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoXsOutList(0, i))
  ELSE IF (nTracerCntl%OutpCntl%RegXsOutPIN(i)) THEN
      WRITE(io16, '(A10, I5, A3, 2(I5, I5, A3),200I7)') 'Region :',  nTracerCntl%OutpCntl%RegXsOutList(1, i), '/', &
                                             nTracerCntl%OutpCntl%RegXsOutList(2:3, i), '/', &
                                             nTracerCntl%OutpCntl%RegXsOutList(4:5, i), '/', (nTracerCntl%OutpCntl%IsoXsOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoXsOutList(0, i))
  ELSE IF (nTracerCntl%OutpCntl%RegXsOutFXR(i)) THEN
      WRITE(io16, '(A10, I5, A3, 3(I5, I5, A3),200I7)') 'Region :',  nTracerCntl%OutpCntl%RegXsOutList(1, i), '/', &
                                             nTracerCntl%OutpCntl%RegXsOutList(2:3, i), '/', &
                                             nTracerCntl%OutpCntl%RegXsOutList(4:5, i), '/', &
                                             nTracerCntl%OutpCntl%RegXsOutList(6:7, i), '/', (nTracerCntl%OutpCntl%IsoXsOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoXsOutList(0, i))
  ENDIF
  
  WRITE(io16, '(A5, 2A15, 600I16)') 'GRP', 'Flux', 'Area', (nTracerCntl%OutpCntl%IsoXsOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoXsOutList(0, i)),&
      (nTracerCntl%OutpCntl%IsoXsOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoXsOutList(0, i)),&
      (nTracerCntl%OutpCntl%IsoXsOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoXsOutList(0, i)),&
      (nTracerCntl%OutpCntl%IsoXsOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoXsOutList(0, i)),&
      (nTracerCntl%OutpCntl%IsoXsOutList(j, i),j = 1, nTracerCntl%OutpCntl%IsoXsOutList(0, i))

  DO ig = 1, ng, 1
    WRITE(io16, '(I5, 600(1pe16.6))') ig, Outflux(ig, i), Area(i),  &
                                (OutXsA  (j, ig, i), j = 1, nTracerCntl%OutpCntl%IsoXsOutList(0, i)), &
                                (OutXsS  (j, ig, i), j = 1, nTracerCntl%OutpCntl%IsoXsOutList(0, i)), &
                                (OutXsStr(j, ig, i), j = 1, nTracerCntl%OutpCntl%IsoXsOutList(0, i)), &
                                (OutXsF  (j, ig, i), j = 1, nTracerCntl%OutpCntl%IsoXsOutList(0, i)), &
                                (OutXsNF (j, ig, i), j = 1, nTracerCntl%OutpCntl%IsoXsOutList(0, i))

  ENDDO
  
  DO j = 1, (35+15*nTracerCntl%OutpCntl%IsoXsOutList(0, i))
    WRITE(io16, '(A)', advance='NO') '-'
  END DO
  
  WRITE(io16, *)
  WRITE(io16, '(A35, 50(1pe15.5))') 'N.D :', (NUM(j, i), j = 1, nTracerCntl%OutpCntl%IsoXsOutList(0, i))
  WRITE(io16, *)
ENDDO

close(io16)

DEALLOCATE(OutXSA)
DEALLOCATE(OutXSS)
DEALLOCATE(OutXSStr)
DEALLOCATE(OutXSF)
DEALLOCATE(OutXSNF)
DEALLOCATE(OutFlux)
DEALLOCATE(Area)
DEALLOCATE(NUM)
! ----------------------------------------------------

END SUBROUTINE HexProcessEffXs
! ------------------------------------------------------------------------------------------------------------
!                                     02. FXR XS Cnds Numberator
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexFxrXSCnds_Numerator(Core, FmInfo, ThInfo, GroupInfo,nTracerCntl, PE, OutXSA, OutXSS, OutXSStr,OutXsF,OutXSNF,OutFlux,Num,AREA,i,ir1,ir2,ixy,iz,ng)

USE PARAM, only : CKELVIN
USE TYPEDEF
USE CNTL, ONLY : nTracerCntl_Type
USE MacXsLib_mod, ONLY : EffMacXS, MacXsBase
USE MOC_Mod, ONLY : FxrAvgPhi
USE XsUtil_mod, ONLY : GetXsMacDat, ReturnXsMacDat
USE ioutil

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(FxrInfo_Type), POINTER :: myFxr
TYPE(XsMac_Type), POINTER :: XsMac
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

INTEGER, INTENT(IN) :: i,ir1,ir2,ixy,iz,ng

REAL,INTENT(OUT) :: AREA
REAL,INTENT(OUT) :: OutXSA(200, ng),OutXSS(200, ng),OutXSStr(200, ng),OutXSF(200, ng),OutXSNF(200, ng)
REAL,INTENT(OUT) :: OutFlux(ng)
REAL,INTENT(OUT) :: Num(200)
REAL :: Flux(ng)

LOGICAL :: FlagF

INTEGER :: ir, j, k, id, ig, iresoGrp1, iresoGrp2, ifxr

!Fxr Data
INTEGER :: niso
LOGICAL :: lres
INTEGER, POINTER :: idiso(:)
REAL, POINTER :: pnum(:), IsoXsMacfold(:,:), IsoXsMacNfold(:,:)
REAL :: eigv
REAL :: Temp
! ----------------------------------------------------

eigv      = 1.0
iresoGrp1 = GroupInfo%nofg + 1
iresoGrp2 = GroupInfo%nofg + GroupInfo%norg

CALL GetXsMacDat(XsMac, ng, .TRUE.)

AREA     = 0
OutFlux  = 0
OutXSA   = 0
OutXSS   = 0
OutXSStr = 0
OutXsF   = 0
OutXSNF  = 0
NUM      = 0

DO ir = ir1, ir2
  ifxr   = Core%Pin(ixy)%FxrIdxst + ir - 1
  myFxr => FmInfo%Fxr(ifxr, iz)
  temp   = myFxr%temp
  niso   = myFxr%niso
  idiso => myFxr%idiso
  
  ALLOCATE (IsoXsMacfold  (niso,ng))
  ALLOCATE (IsoXsMacNfold (niso,ng))
  
  IsoXsMacfold  = 0
  IsoXsMacNfold = 0
  
  pnum => myFxr%pnum
  lres  = myFxr%lres
  
  CALL MacXsBase(XsMac, myFxr, 1, ng, ng, eigv, .FALSE., .TRUE.)
  
  IsoXsMacfold  = XsMac%IsoXsMacf
  IsoXsMacnfold = XsMac%IsoXsMacNf
  AREA = AREA + myFxr%area
  
  IF (lres) THEN
    DO ig = iresoGrp1, iresoGrp2
      DO k = 1, niso
        XsMac%IsoXsMacA(k, ig) = XsMac%IsoXsMacA(k, ig) * myFxr%fresoAIso(k,ig)
        
        IF (IsoXsMacFold(k,ig).ne.0) XsMac%IsoXsMacF(k, ig) = XsMac%IsoXsMacF(k, ig) * myFxr%fresoFIso(k,ig)
      ENDDO
    ENDDO
    
    IF (nTracerCntl%lRST) THEN
      DO ig = iresoGrp1, iresoGrp2
        DO k = 1, niso
          XsMac%IsoXsMacS0(k, ig) = XsMac%IsoXsMacS0(k, ig) * myFxr%fresoSIso(k,ig)
          XsMac%IsoXsMacS1(k, ig) = XsMac%IsoXsMacS1(k, ig) * myFxr%fresoS1Iso(k,ig)
        ENDDO
      ENDDO
    ENDIF
  ENDIF
  
  Flux    = FxrAvgPhi(Core, FmInfo%Fxr, FMInfo%Phis, ixy, ir, iz, ng, PE)
  Outflux = OutFlux + myFxr%area * Flux
  
  DO j = 1, nTracerCntl%OutpCntl%IsoXsOutList(0, i)
    id = nTracerCntl%OutpCntl%IsoXsOutList(j, i)
    
    DO k = 1, niso
      IF(IDISO(k) .NE. id) CYCLE
      
      NUM(j) = NUM(j) + myFxr%area * pnum(k)
      
      DO ig = 1, ng
        OutXsA  (j, ig) = OutXsA  (j, ig) + myFxr%area * Flux(ig) * XsMac%IsoXsMacA (k, ig)
        OutXsS  (j, ig) = OutXsS  (j, ig) + myFxr%area * Flux(ig) * XsMac%IsoXsMacS0(k, ig)
        OutXsStr(j, ig) = OutXsStr(j, ig) + myFxr%area * Flux(ig) * (XsMac%IsoXsMacS0(k, ig) - XsMac%IsoXsMacS1(k, ig))
        
        IF (IsoXsMacFold(k,ig) .NE. 0) THEN
          OutXsNF(j, ig) = OutXsNF(j, ig) + myFxr%area * XsMac%IsoXsMacF(k, ig) * Flux(ig)
          OutXsF (j, ig) = OutXsF (j, ig) + myFxr%area * XsMac%IsoXsMacF(k, ig) * Flux(ig)
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  
  DEALLOCATE (IsoXsMacNFold,IsoXsMacFold)
ENDDO

NULLIFY(IDISO, pnum, myFxr)

CALL ReturnXsMacDat(XsMac)
! ----------------------------------------------------

END SUBROUTINE HexFxrXSCnds_Numerator