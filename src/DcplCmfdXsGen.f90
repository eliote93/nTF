#include <defines.h>
SUBROUTINE DcplHomXsGen(Core, DcplInfo, Fxr, PinXs, myzb, myze, ng, lXsGen, lXsLib, lXsFtn, XsFtnMod, GroupInfo)   
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type    ,DcplInfo_Type   ,PinXs_Type      &
                          ,FmInfo_Type      ,CmInfo_Type                      &
                          ,CellXsData_TYPE  ,GroupInfo_Type  ,FxrInfo_Type    &
                          ,Pin_Type         ,PinInfo_Type    ,Cell_Type       &
                          ,FmInfo_Type      ,CmInfo_Type     ,DhatData_Type
USE DcplCmfd_Mod,   ONLY : DcplCellXsGen    ,DcplCellXsGenXsFtn
USE DcplMicXsFtn_mod, ONLY : UpdtMacXs_MicXsFtn
USE CMFD_MOD,       ONLY : AllocPinXS
USE Th_Mod,         ONLY : GetPinFuelTemp,   GetPinModTemp,     GetPinTemp    &
                          ,CalPlnTemp
USE BasicOperation, ONLY : CP_VA 
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(DcplInfo_Type) :: DcplInfo
TYPE(PinXs_Type), POINTER :: PinXS(:, :)

TYPE(GroupInfo_Type) :: GroupInfo
INTEGER :: myzb, myze, ng,XsFtnMod
LOGICAL :: lXsLib, lXsFtn, lXsGen

TYPE(CellXsData_Type), SAVE :: CellXsData
TYPE(FmInfo_Type), POINTER :: FmInfo(:, :)
TYPE(CmInfo_Type), POINTER :: CmInfo(:, :)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(PinInfo_Type), POINTER :: PinInfo(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(PinXs_Type), POINTER :: RefPinXs



INTEGER, PARAMETER :: nbd = 4
#ifdef DcplMemSave
TYPE(PinXs_Type), SAVE :: PinXsFtnData(3)
#else
TYPE(PinXs_Type) :: PinXsFtnData(3)
#endif
INTEGER :: nCoreFxr, nCoreFsr, nxy
INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nlocalFsr
INTEGER :: ixy, iz, icel, iRefPln, iRefPln0
INTEGER :: i, j, ibeg, iend
LOGICAL :: lFuelPln
REAL :: PlnTemp(100), RefPlnTemp(5, 100)
REAL :: CellPhis(400, ng, 3)

IF(.NOT. CellXsData%lalloc) THEN
  CALL AllocDcplCellXs(CellXsData, 3, 400, ng)
#ifdef DcplMemSave  
  CALL AllocPinXS(PinXsFtnData(1), ng, 4, GroupInfo%InScatRange)
  CALL AllocPinXS(PinXsFtnData(2), ng, 4, GroupInfo%InScatRange)
  CALL AllocPinXS(PinXsFtnData(3), ng, 4, GroupInfo%InScatRange)
#endif
ENDIF


Pin => Core%Pin
PinInfo => Core%Pininfo;   Cell => Core%CellInfo

FmInfo => DcplInfo%DcplFmInfo;  CmInfo => DcplInfo%DcplCmInfo

nxy = Core%nxy
nCoreFxr = Core%nCoreFxr; nCoreFsr = Core%nCoreFsr
!myzb = PE%myzb; myze = PE%myze

!Get Pin Temperature Information
IF(lXSLib) THEN
  DO iz = myzb, myze
    DO ixy = 1, nxy
      PinXS(ixy,iz)%PinTemp = GetPinTemp(Core, Fxr, iz, ixy)
      IF(Core%lFuelPlane(iz) .AND. Pin(ixy)%lFuel) THEN
        PinXS(ixy,iz)%FuelTemp = GetPinFuelTemp(Core, Fxr, iz, ixy)
        PinXS(ixy,iz)%ModTemp = GetPinModTemp(Core, Fxr, iz, ixy)
      ENDIF
    ENDDO
  ENDDO  
ENDIF

DO iz = myzb, myze
  iRefPln0 = DcplInfo%Pln_Map0(iz)
  iRefPln = DcplInfo%Pln_Map(iz)
  lFuelPln = Core%lFuelPlane(iz)
  DO ixy = 1, nxy
    FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
    icel = Pin(ixy)%Cell(iz)
    nlocalFxr = Cell(icel)%nFxr; nlocalFsr = Cell(icel)%nFsr
   
    IF(.NOT. lXSGen .AND. lXsFTN) THEN
       
       
#ifdef DcplMemSave
      ibeg = FsrIdxSt; iend = FsrIdxst + nlocalFsr - 1
      DO j = 1, 3
        CALL CP_VA(CellPhis(1:nLocalFsr,1:ng, j), FmInfo(j, iRefPln0)%phis(ibeg:iend, iRefPln, 1:ng), nlocalFsr, ng) 
      ENDDO
      ibeg = FxrIdxSt; iend = FxrIdxst + nlocalFxr - 1
      !Reference Cell XS Set
      CALL DcplCellXsGenXsFtn(PinXsFtnData(1:3), FmInfo(1, iRefPln0)%FXR(ibeg:iend,iRefPln), Cell(icel), &
                              CellPhis(1:nLocalFsr, 1:ng, 1:3), ng, nLocalFsr, nLocalFxr, lxslib, lfuelpln, GroupInfo)
      !Copy Base XS
      CALL CopyDcplPinXS(PinXsFtnData(1), PinXS(ixy, iz), ng, 4)                              
#else  
      DO j = 1, 3
        PinXsFtnData(j)%XSTR => CMInfo(j, iRefPln0)%PinXS(ixy, iz)%XSTR
        PinXsFtnData(j)%XST => CMInfo(j, iRefPln0)%PinXS(ixy, iz)%XST
        PinXsFtnData(j)%XSNF => CMInfo(j, iRefPln0)%PinXS(ixy, iz)%XSNF
        PinXsFtnData(j)%XSKF => CMInfo(j, iRefPln0)%PinXS(ixy, iz)%XSKF
        PinXsFtnData(j)%CHI => CMInfo(j, iRefPln0)%PinXS(ixy, iz)%CHI
        PinXsFtnData(j)%XSS => CMInfo(j, iRefPln0)%PinXS(ixy, iz)%XSS
        
        PinXsFtnData(j)%PinTemp = CMInfo(j, iRefPln0)%PinXS(ixy, iz)%PinTemp
        PinXsFtnData(j)%FuelTemp = CMInfo(j, iRefPln0)%PinXS(ixy, iz)%FuelTemp
        PinXsFtnData(j)%ModTemp = CMInfo(j, iRefPln0)%PinXS(ixy, iz)%ModTemp
      ENDDO
      CALL CopyDcplPinXS(CmInfo(1, iRefPln0)%PinXS(ixy, iz), PinXS(ixy, iz), ng, 4)
#endif     

      !updt temperature change

      IF(XsFtnMod .EQ. 2) THEN
        ibeg = FxrIdxSt; iend = FxrIdxst + nlocalFxr - 1
        !CALL UpdtMacXs_MicXsFtn(PinXS(ixy, iz), DcplInfo%MicXsFtn(ixy, iRefPln0),         &
        !                        FmInfo(1, iRefPln0)%FXR(ibeg:iend,iRefPln0), Cell(icel),  &
        !                        nLocalFxr, GroupInfo)    
        CALL UpdtMacXs_MicXsFtn(PinXS(ixy, iz), DcplInfo%MicXsFtn(ixy, iRefPln0),         &
                                FXR(ibeg:iend,iz), Cell(icel), nLocalFxr, GroupInfo)              
      ELSE
        CALL PinXsFtn(Cell(icel), PinXs(ixy, iz), PinXsFtnData, ng, ixy, iz)
      ENDIF
     
    ELSEIF(.NOT. lXsGen) THEN
      CALL CopyDcplPinXS(CmInfo(1, iRefPln0)%PinXS(ixy, iz), PinXS(ixy, iz), ng, 4)
    ELSE
      CellXsData%CellInfo => Cell(icel)
      !Cell Phi distribution
      CellXSData%ntemp = 1
      ibeg = FsrIdxSt; iend = FsrIdxst + nlocalFsr - 1
      CALL CP_VA(CellXsData%phis(1, 1:nlocalFsr, 1:ng), FmInfo(1, iRefPln0)%phis(ibeg:iend, iRefPln, 1:ng), nlocalFsr, ng )
      ibeg = FxrIdxSt; iend = FxrIdxst + nlocalFxr - 1
      CellXsData%FXR => Fxr(ibeg:iend, iz)
      CALL DcplCellXsGen(PinXs(ixy, iz), CellXsData, ng, FALSE,                          &
                         lXsLib, Core%lFuelPlane(iRefPln), GroupInfo)      
!      CALL DcplCellXsGen(PinXs(ixy, iz), CellXsData, iRefPln, iRefPln0, ng, FALSE,   &
!                         lXsLib, Core%lFuelPlane(iRefPln), GroupInfo)
    ENDIF
  ENDDO
ENDDO
continue
END SUBROUTINE

SUBROUTINE PinXsFtn(CellInfo, PinXS, PinXSFtnData, ng, ixy, iz)
USE PARAM
USE TYPEDEF,        ONLY : PinXs_Type     ,Cell_TYPE
USE BasicOperation, ONLY : AD_VA         ,CP_CA          ,Multi_CA             &
                          ,SUB_VA
IMPLICIT NONE
TYPE(Cell_Type) :: CellInfo
TYPE(PinXs_Type) :: PinXS, PinXSFtnData(3)
INTEGER :: ng, ixy ,iz

REAL :: RefT0, RefFuelT0, RefModT0
REAL :: RefT1, RefFuelT1, RefModT1
REAL :: T, FuelT, ModT
REAL :: RDFT, RDMT, DFT, DMT, FRAC
REAL :: DelXs(ng, 5),DelXss(ng), DelXs0(ng, 5), DelXss0(ng,ng)

INTEGER :: ig, ib, ie
RefT0 = PinXSFtnData(1)%PinTemp
RefFuelT0 = PinXSFtnData(1)%FuelTemp; RefModT0 = PinXSFtnData(1)%ModTemp
RefT1 = PinXSFtnData(3)%PinTemp
RefFuelT1 = PinXSFtnData(2)%FuelTemp; RefModT1 = PinXSFtnData(3)%ModTemp
T = PinXs%PinTemp; FuelT = PinXs%FuelTemp; ModT =PinXs%ModTemp

IF(.NOT. CellInfo%lFuel) THEN
  RefModT0 = RefT0
  RefModT1 = RefT1
  ModT = T
ENDIF

DO ig = 1, ng
  PinXs%Xss(ig)%from(ig) = PinXs%Xss(ig)%WithInGroupScat
ENDDO
CALL CP_CA(DelXs0, 0._8, ng, 5)
CALL CP_CA(DelXss0, 0._8, ng, ng)

DMT = (ModT - RefModT0)
RDMT = (RefModT1 - RefModT0)

DelXs(1:ng, 1) = PinXsFtnData(3)%xst - PinXsFtnData(1)%xst
DelXs(1:ng, 2) = PinXsFtnData(3)%xstr - PinXsFtnData(1)%xstr
DelXs(1:ng, 3) = PinXsFtnData(3)%xsnf - PinXsFtnData(1)%xsnf
DelXs(1:ng, 4) = PinXsFtnData(3)%xskf - PinXsFtnData(1)%xskf
FRAC = 0
IF(abs(RDMT) .GT. epsm5) THEN
  FRAC = DMT / RDMT
ENDIF

CALL MULTI_CA(FRAC, DelXS(1:ng, 1:4), ng, 4)

!WRITE(64, '(2I8, 3F15.7, 200e20.5)') iz, ixy, ModT,RefModT0, RefModT1, FRAC, DelXS(1:ng, 1)
!DelXS(1:ng, 1:2) = 0; FRAC = 0;
CALL AD_VA(DelXs0(1:ng, 1:4), DelXs0(1:ng, 1:4), DelXS(1:ng, 1:4), ng, 4)

DO ig = 1, ng
  ib = PinXs%Xss(ig)%ib; ie = PinXs%Xss(ig)%ie
  DelXss(ib:ie) = PinXsFtnData(3)%Xss(ig)%from(ib:ie) - PinXsFtnData(1)%Xss(ig)%from(ib:ie)
  DelXss(ig) = PinXsFtnData(3)%Xss(ig)%WithInGroupScat - PinXsFtnData(1)%Xss(ig)%WithInGroupScat
  CALL MULTI_CA(Frac, DelXss(ib:ie), ie-ib+1)
  CALL AD_VA(DelXss0(ib:ie, ig), DelXss0(ib:ie, ig), DelXSs(ib:ie), ie-ib+1)
ENDDO

IF(CellInfo%lFuel) THEN
  DFT = (SQRT(FuelT) - SQRT(RefFuelT0))
  RDFT = (SQRT(RefFuelT1) - SQRT(RefFuelT0))
  !FRAC = 0
  IF(abs(RDFT) .GT. epsm5) THEN
    FRAC = DFT / RDFT
  ELSE
    PRINT *, RefFuelT1, RefFuelT0
    PAUSE
  ENDIF
  DelXs(1:ng, 1) = PinXsFtnData(2)%xst - PinXsFtnData(1)%xst
  DelXs(1:ng, 2) = PinXsFtnData(2)%xstr - PinXsFtnData(1)%xstr
  DelXs(1:ng, 3) = PinXsFtnData(2)%xsnf - PinXsFtnData(1)%xsnf
  DelXs(1:ng, 4) = PinXsFtnData(2)%xskf - PinXsFtnData(1)%xskf
  CALL MULTI_CA(FRAC, DelXS(1:ng, 1:4), ng, 4)
  CALL AD_VA(DelXs0(1:ng, 1:4), DelXs0(1:ng, 1:4), DelXS(1:ng, 1:4), ng, 4)
  DO ig = 1, ng
    ib = PinXs%Xss(ig)%ib; ie = PinXs%Xss(ig)%ie
    DelXss(ib:ie) = PinXsFtnData(2)%Xss(ig)%from(ib:ie) - PinXsFtnData(1)%Xss(ig)%from(ib:ie)
    DelXss(ig) = PinXsFtnData(2)%Xss(ig)%WithInGroupScat - PinXsFtnData(1)%Xss(ig)%WithInGroupScat
    CALL MULTI_CA(Frac, DelXss(ib:ie), ie-ib+1)
    CALL AD_VA(DelXss0(ib:ie, ig), DelXss0(ib:ie, ig), DelXSs(ib:ie), ie-ib+1)
  ENDDO
ENDIF

!UPDATE XS
CALL AD_VA(PinXS%xst(1:ng), PinXS%xst(1:ng), DelXs0(1:ng, 1), ng)
CALL AD_VA(PinXS%xstr(1:ng), PinXS%xstr(1:ng), DelXs0(1:ng, 2), ng)
IF(CellInfo%lFuel) THEN
  CALL AD_VA(PinXS%xsnf(1:ng), PinXS%xsnf(1:ng), DelXs0(1:ng, 3), ng)
  CALL AD_VA(PinXS%xskf(1:ng), PinXS%xskf(1:ng), DelXs0(1:ng, 4), ng)
ENDIF

DO ig = 1, ng
  ib = PinXs%Xss(ig)%ib; ie = PinXs%Xss(ig)%ie
  CALL AD_VA(PinXs%Xss(ig)%from(ib:ie), PinXs%Xss(ig)%from(ib:ie), DelXss0(ib:ie, ig), ie - ib +1)
  PinXS%Xss(ig)%WithInGroupScat = PinXs%Xss(ig)%from(ig)
  PinXs%Xss(ig)%from(ig) = 0
ENDDO
!SET RMV XS
DO ig = 1, ng
  PinXS%xsr(ig) = PinXS%xstr(ig) - PinXs%Xss(ig)%WithInGroupScat
  PinXS%XSD(ig) = 1._8/(3._8 * PinXS%xstr(ig))
ENDDO

END SUBROUTINE

SUBROUTINE CopyDcplPinXS(DcplPinXs, PinXs, ng, nbd)
USE PARAM
USE TYPEDEF,        ONLY : PinXS_Type
USE BASICOPERATION, ONLY : CP_VA
IMPLICIT NONE
TYPE(PinXs_Type) :: DcplPinXS, PinXs
INTEGER :: ng, nbd
INTEGER :: ib, ie, ig

CALL CP_VA(PinXS%XSD(1:ng), DcplPinXS%XSD(1:ng), ng)
CALL CP_VA(PinXS%XST(1:ng), DcplPinXS%XST(1:ng), ng)
CALL CP_VA(PinXS%XSTR(1:ng), DcplPinXS%XSTR(1:ng), ng)
CALL CP_VA(PinXS%XSNF(1:ng), DcplPinXS%XSNF(1:ng), ng)
CALL CP_VA(PinXS%XSKF(1:ng), DcplPinXS%XSKF(1:ng), ng)
CALL CP_VA(PinXS%CHI(1:ng), DcplPinXS%CHI(1:ng), ng)

CALL CP_VA(PinXS%DHAT(1:nbd, 1:ng), DcplPinXS%DHAT(1:nbd, 1:ng), nbd, ng)
!CALL CP_VA(PinXS%PDHAT(1:nbd, 1:ng), DcplPinXS%PDHAT(1:nbd, 1:ng), nbd, ng)
!CALL CP_VA(PinXS%DTIL(1:nbd, 1:ng), DcplPinXS%DTIL(1:nbd, 1:ng), nbd, ng)
!Scattering Mat
DO ig = 1, ng
  ib =DcplPinXs%XSS(ig)%ib; ie = DcplPinXs%XSS(ig)%ie
  PinXs%XSS(ig)%ib = ib; PinXs%XSS(ig)%ie = ie
  PinXs%XSS(ig)%WithInGroupScat = DcplPinXs%XSS(ig)%WithInGroupScat
  PinXS%XSS(ig)%from(ib : ie) = DcplPinXs%XSS(ig)%from(ib : ie)
  PinXS%XSR(ig) = PinXS%XSTR(ig) - PinXs%XSS(ig)%WithInGroupScat
  PinXS%XSS(ig)%from(ig) = 0._8
!  DcplPinXs%XSS(ig)%from(Ig) = 0
!  DcplPinXS%xsr(ig) =DcplPinXS%xstr(ig) - DcplPinXS%XSS(ig)%WithInGroupScat
ENDDO
END SUBROUTINE


SUBROUTINE DcplCellXsGenXsFtn(PinXS, Fxr, CellInfo, phis, ng, nFsr, nFxr, lxslib, lfuelpln, GroupInfo)
USE PARAM
USE TYPEDEF,        ONLY : PinXs_Type      ,FxrInfo_Type        ,GroupInfo_Type       &
                          ,Cell_Type       ,CellXsData_Type
USE CMFD_mod,       ONLY : XsMac           ,HomoCellXsGen
USE BenchXs,        ONLY : XsBaseBen
USE MacXsLib_Mod,   ONLY : MacXsBase       ,MacXsScatMatrix
USE DcplCmfd_Mod,   ONLY : DcplCellXsGen
USE BasicOperation, ONLY : CP_VA            ,CP_CA             ,MULTI_VA
IMPLICIT NONE

TYPE(PinXS_TYPE) :: PinXS(3)
TYPE(FxrInfo_Type), TARGET :: Fxr(nFxr)
TYPE(Cell_Type), TARGET :: CellInfo
TYPE(GroupInfo_Type) :: GroupInfo

REAL :: phis(nfsr, ng, 3)
INTEGER :: nFsr, nFxr, ng
LOGICAL :: lXsLib, lFuelPln

TYPE(CellXsData_Type), SAVE :: CellXsData
REAL :: areasum, area, AreaFsum, AreaMsum
REAL :: Tsum, Tfsum, Tmsum
INTEGER :: imode, ifxr
INTEGER :: niso
INTEGER :: norg, nchi, iResoGrpBeg, iResoGrpEnd

IF(.NOT. CellXsData%lalloc)  CALL AllocDcplCellXs(CellXsData, 3, 400, ng)

IF(lXsLib) THEN
  norg = GroupInfo%norg; nCHI = GroupInfo%nCHi
  iResoGrpBeg = GroupInfo%nofg + 1
  iResoGrpEnd = GroupInfo%nofg + GroupInfo%norg 
ENDIF

DO imode = 1, 3
  IF(.NOT. lFuelPln .AND. imode .EQ. 2) CYCLE
  CellXsData%CellInfo => CellInfo
  CellXSData%ntemp = 1
  CALL CP_VA(CellXsData%phis(1, 1:nFsr, 1:ng), phis(1:nFsr, 1:ng, imode), nFsr, ng )
  DO ifxr = 1, nFxr
    niso = Fxr(ifxr)%niso
    Fxr(ifxr)%pnum(1:niso) = Fxr(ifxr)%Dcpl_pnum(1:niso, imode)
    Fxr(ifxr)%Temp = Fxr(ifxr)%Dcpl_Temp(imode)
    Fxr(ifxr)%FresoA(iResoGrpBeg:iResoGrpEnd) = Fxr(ifxr)%Dcpl_FresoA(iResoGrpBeg:iResoGrpEnd, imode)
    Fxr(ifxr)%FresoS(iResoGrpBeg:iResoGrpEnd) = Fxr(ifxr)%Dcpl_FresoS(iResoGrpBeg:iResoGrpEnd, imode)
    Fxr(ifxr)%FresoF(iResoGrpBeg:iResoGrpEnd) = Fxr(ifxr)%Dcpl_FresoF(iResoGrpBeg:iResoGrpEnd, imode)
  ENDDO
  CellXsData%Fxr => Fxr(1:nFxr)
  CALL DcplCellXsGen(PinXS(imode), CellXSData, ng, .FALSE., lXsLib, lFuelPln, GroupInfo)
  !Determine Temperature
  IF(lXsLib) THEN
    areasum = 0; Tsum = 0
    AreaMsum = 0; AreaFsum = 0
    Tfsum = 0; Tmsum = 0
    DO ifxr = 1, nFxr
      Area = Fxr(ifxr)%area; AreaSum = AreaSum + area
      Tsum = Tsum + Fxr(ifxr)%temp * area
      IF(Fxr(ifxr)%lFuel) THEN
        Tfsum = Tfsum + Fxr(ifxr)%temp * area
        AreaFsum = AreaFsum + area
      ENDIF
      IF(Fxr(ifxr)%lh2o) THEN
        Tmsum = Tmsum + Fxr(ifxr)%temp * area 
        AreaMsum = AreaMsum + area
      ENDIF
    ENDDO
    PinXS(imode)%PinTemp = Tsum / AreaSum
    PinXS(imode)%FuelTemp = 0; PinXS(imode)%ModTemp = 0
    IF(CellInfo%lfuel) THEN
      PinXS(imode)%FuelTemp = Tfsum / AreaFSum
      PinXS(imode)%ModTemp = Tmsum / AreaMSum
    ENDIF
  ENDIF
ENDDO

DO ifxr = 1, nFxr
  niso = Fxr(ifxr)%niso
  Fxr(ifxr)%pnum(1:niso) = Fxr(ifxr)%Dcpl_pnum(1:niso, 1)
  Fxr(ifxr)%Temp = Fxr(ifxr)%Dcpl_Temp(1)
  Fxr(ifxr)%FresoA(iResoGrpBeg:iResoGrpEnd) = Fxr(ifxr)%Dcpl_FresoA(iResoGrpBeg:iResoGrpEnd, 1)
  Fxr(ifxr)%FresoS(iResoGrpBeg:iResoGrpEnd) = Fxr(ifxr)%Dcpl_FresoS(iResoGrpBeg:iResoGrpEnd, 1)
  Fxr(ifxr)%FresoF(iResoGrpBeg:iResoGrpEnd) = Fxr(ifxr)%Dcpl_FresoF(iResoGrpBeg:iResoGrpEnd, 1)
ENDDO
END SUBROUTINE

SUBROUTINE DcplCellXsGen(PinXs, CellXsData, ng, lXsFtn, lxslib, lfuelpln, GroupInfo)
USE PARAM
USE TYPEDEF,        ONLY : PinXs_Type      ,CellXsData_Type     ,GroupInfo_Type        &
                          ,Cell_Type       ,FxrInfo_Type        
USE CMFD_mod,       ONLY : XsMac           ,HomoCellXsGen
USE BenchXs,        ONLY : XsBaseBen
USE MacXsLib_Mod,   ONLY : MacXsBase        ,MacXsScatMatrix  
USE BasicOperation, ONLY : CP_VA            ,CP_CA             ,MULTI_VA
IMPLICIT NONE
TYPE(PinXs_Type) :: PinXs
TYPE(CellXsData_Type) :: CellXsData
TYPE(GroupInfo_Type) :: GroupInfo
INTEGER :: ng, ig
!INTEGER :: iRefPln0, iRefPln, ng
LOGICAL :: lXsFtn, lXsLib, lFuelPln

TYPE(FxrInfo_Type), POINTER :: myFxr
TYPE(Cell_Type), POINTER :: CellInfo
INTEGER :: nFxr, nFsr, norg
INTEGER :: i, j, k
INTEGER :: itype, iResoGrpBeg, iResoGrpEnd, nCHI
REAL :: phis(1000, 400)
CellInfo => CellXsData%CellInfo
nFxr = CellInfo%nFxr; nFsr = CellInfo%nFSR
IF(lXsLib) THEN
  norg = GroupInfo%norg; nCHI = GroupInfo%nCHi
  iResoGrpBeg = GroupInfo%nofg + 1
  iResoGrpEnd = GroupInfo%nofg + GroupInfo%norg 
ENDIF
CALL CP_CA(phis(1:nfsr, 1:ng), 0._8, nFsr, ng)
DO j = 1, nFxr
  XsMac(j)%lFuel = FALSE
  myFxr => CellXsData%Fxr(j)
  IF(lXsLib) THEN
    CALL MacXsBase(XsMac(j), myFxr, 1, ng, ng, 1._8, FALSE, TRUE)
    CALL MacXsScatMatrix(XsMac(j), myFxr, 1, ng, ng, GroupInfo, FALSE)
    !Self-Sheilding Effect
    IF(myFxr%lres) THEN
       IF(.NOT. lXsFtn) THEN
         do ig = iResoGrpBeg, iResoGrpEnd
           XsMac(j)%XsMacA(ig) = XsMac(j)%XsMacA(ig) * myFxr%FresoA(ig)    
         enddo
         IF(lFuelPln) THEN
           do ig = iResoGrpBeg, iResoGrpEnd
             XsMac(j)%XsMacNf(ig) = XsMac(j)%XsMacNf(ig) * myFxr%FresoF(ig)    
             XsMac(j)%XsMacKf(ig) = XsMac(j)%XsMacKf(ig) * myFxr%FresoF(ig)    
           enddo
         ENDIF
       ELSE
         !effective XS genenration routine should be called.
       ENDIF
       XsMac(j)%XsMacTr = XsMac(j)%XsMacA + XsMac(j)%XsMacStr
       XsMac(j)%XsMacT = XsMac(j)%XsMacA + XsMac(j)%XsMacS
       
       CALL CP_CA(XsMac(j)%CHI, 0._8, ng)
       IF(myFxr%lDepl) THEN
         CALL CP_VA(XsMac(j)%CHI(1:nCHI), myFxr%CHI(1:nCHI), nCHI)
         XsMac(j)%lFuel = TRUE
       ENDIF     
    ENDIF    
  ELSE
    i = CellInfo%MapFxr2FsrIdx(1,j); !itype = CellInfo%iReg(i)
    itype=myFxr%imix
    CALL xsbaseBen(itype, 1, ng, 1, ng, FALSE, XsMac(j))
  ENDIF
ENDDO
!
IF(.NOT. lXsFtn) THEN
  CALL CP_VA(phis(1:nFsr, 1:ng), CellXsData%phis(1, 1:nfsr, 1:ng), nfsr, ng)
ELSE
! linear interpolation for given temperature
ENDIF
CALL HomoCellXsGen(CellInfo, Phis(1:nFsr, 1:ng), PinXs, XsMac(1:nFxr), ng, nFxr, nFsr,   &
                   GroupInfo%OutScatRange, lXsLib, FALSE)
!    CALL HomoCellXsGen(Cell(icel), Phis(FsrIdxSt:FsrIdxSt+nlocalFsr-1, iz, 1:ng), PinXS(ixy, iz),  &
!                       XsMac(1:nLocalFxr), ng, nLocalFxr, nLocalFsr, GroupInfo%OutScatRange, lXsLib)
END SUBROUTINE


SUBROUTINE AllocDcplCellXs(CellXsData, nRefTemp, nFsrMax, ng)
USE PARAM
USE Typedef,  ONLY : CellXSData_Type
IMPLICIT NONE
TYPE(CellXSData_Type) :: CellXsData
INTEGER :: nRefTemp, nFsrMax, ng
CellXsData%nFsrMax = nFsrMax
CellXsData%nTempMax = nRefTemp
IF(CellXsData%lAlloc) DEALLOCATE(CellXsData%phis)
ALLOCATE(CellXsData%phis(nRefTemp, nFsrMax, ng))
CellXsData%lalloc = .TRUE.
END SUBROUTINE

SUBROUTINE DcplAxSrcUpdt(Core, CmInfo, DcplInfo, ng, PE, lXsFtn)
USE PARAM
USE TYPEDEF, ONLY : CoreINfo_Type   ,CmInfo_Type      ,PE_TYPE      &
                   ,DcplInfo_Type   ,PE_Type                        &    
                   ,Pin_Type        ,AxFlx_Type       ,PinXs_Type  
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(DcplInfo_Type) :: DcplInfo
TYPE(PE_TYPE) :: PE
INTEGER :: ng
LOGICAL :: lXsFtn

TYPE(PinXS_Type), POINTER :: PinXS(:, :)
TYPE(AxFlx_Type), POINTER :: AxFlx(:)
REAL, POINTER :: PhiC(:, :, :), hz(:), hzInv(:)
REAL, POINTER :: AxSrc(:, :, :), AxPXS(:, :, :)
REAL, POINTER :: DcplAxSrc(:, :, :), DcplPXS(:, :, :)
REAL :: Lkg(400)
REAL :: neighphi, myphi, jnet, Dtil, Dhat
REAL :: hsum, srcsum, phisum
INTEGER :: i, j, k, ig
INTEGER :: iz, iz0, iz1, iz2
INTEGER :: ixy, npin, nxy

INTEGER :: iRefPln, iTempPln
INTEGER :: nRefPln, nTempPln
INTEGER :: myRefPlnBeg, myRefPlnEnd
nxy = Core%nxy
PinXS => CmInfo%PinXS; PhiC => CmInfo%PhiFM
hz => Core%hz; HzInv => Core%HzInv

AxSrc => CmInfo%AxSrc
nRefPln = DcplInfo%nRefPln; nTempPln = DcplInfo%nRefTemp
#ifndef MPI_ENV
myRefPlnBeg = 1; myRefPlnEnd = nRefPln
#else
myRefPlnBeg = PE%myRefPlnBeg; myRefPlnEnd = PE%myRefPlnEnd
#endif
IF(.NOT. lXsFtn) THEN
  DO iRefPln = myRefPlnBeg, myRefPlnEnd
    DcplPxs => DcplInfo%DcplFmInfo(1, iRefPln)%AxPXS
    !PhiC => DcplInfo%DcplCmInfo(1, iRefPln)%PhiC
    iz0 = DcplInfo%RefPln(iRefpln)
    CALL CP_CA(DcplPxs(1:nxy, iz0:iz0, 1:ng), 0._8, nxy, 1, ng)
    IF(.NOT. Core%lFuelPlane(iz0)) CYCLE
    iz1 = DcplInfo%RefPln_map(1, iRefPln); iz2 = DcplInfo%RefPln_map(2, iRefPln)
    hsum = sum(hz(iz1:iz2))
    DO ig = 1, ng
      DO ixy = 1, nxy
        srcsum = 0; phisum = 0
        DO iz = iz1, iz2
          !srcsum = srcsum + hz(iz) * AxSrc(ixy, iz , ig)
          phisum = phisum + hz(iz) * PhiC(ixy, iz, ig)
        ENDDO
        myPhi = PhiC(ixy, iz1, ig); neighPhi = PhiC(ixy, iz1-1, ig) 
      !Dtil = AxFlx(i)%Dtil(1, ig); Dhat = AxFlx(i)%Dhat(1, ig)
        Dtil = AxFlx(iz1)%Dtil(1, ig); Dhat = AxFlx(iz1)%Dhat(1, ig)
        Jnet = (Dtil - Dhat) * MyPhi - (Dtil + Dhat) * NeighPhi
        srcsum = Jnet
        myPhi = PhiC(ixy, iz2, ig); neighPhi = PhiC(ixy, iz2+1, ig)
        Dtil = AxFlx(iz2)%Dtil(2, ig); Dhat = AxFlx(iz2)%Dhat(2, ig)
        !Dtil = AxFlx(i, ixy)%Dtil(2, ig); Dhat = AxFlx(i, ixy)%Dhat(2, ig)
        Jnet = (Dtil - Dhat) * MyPhi - (Dtil + Dhat) * NeighPhi
        srcsum = srcsum + Jnet        
          
        !srcsum = srcsum / hsum; phisum = phisum / hsum
        DcplPxs(ixy, iz0, ig) = srcsum / phisum 
      ENDDO
    ENDDO
  ENDDO
  CONTINUE
ELSE

ENDIF
END SUBROUTINE

SUBROUTINE DcplRadCouplingCoeffGen(Pin, DcplInfo, lXsFtn, lXsGen)
USE PARAM
USE TYPEDEF,      ONLY : PinXs_TYPE     ,Pin_Type        ,DcplInfo_Type  &
                        ,CmInfo_Type
USE CMFD_MOD,     ONLY : ng             ,nxy             ,myzb           &
                        ,myze           ,hzfm            ,hz             &
                        ,CmfdPinXS      ,PinNeighIdx
IMPLICIT NONE
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(DcplInfo_Type) :: DcplInfo
LOGICAL :: lXsFtn, lXsGen

TYPE(PinXS_Type), POINTER :: PinXS
TYPE(PinXS_Type), POINTER :: RefPinXS
TYPE(CmInfo_Type), POINTER :: DcplCmInfo(:, :)
REAL, POINTER :: Jout(:, :, :, :, :)
INTEGER :: irefpln, irefpln0
INTEGER :: ixy, ineigh, ig, iz, ibd, inbd
INTEGER, PARAMETER :: nbd = 4
REAL :: myphi(ng), neighphi, phisum, jfdm, jnet
REAL :: PDHAT, Del, Alpha
REAL :: Dhat, Dtil, mybeta, neighbeta, smy
REAL :: DhatDat(5, 2,  ng)


DcplCmInfo => DcplInfo%DcplCmInfo

DO iz = myzb, myze
  iRefPln0 = DcplInfo%Pln_Map0(iz)
  iRefPln = DcplInfo%Pln_Map(iz)
  DO ixy = 1, nxy
    PinXS => CmfdPinXS(ixy, iz)
    myphi(1:ng) = PinXS%phi(1:ng)
    !
    DO ibd = 1, nbd
      ineigh = PinNeighIdx(ibd, ixy)
      inbd = Pin(ixy)%NeighSurfIdx(ibd)    !The Corresponding Surface Index of 
      smy = Pin(ixy)%BdLength(ibd) !* hz(iz)
      IF(ineigh .GT. 0) THEN
        DO ig = 1, ng
         
          mybeta = CmfdPinXS(ixy, iz)%XSD(ig) / Pin(ixy)%Center2SurfaceL(ibd)
          neighbeta = CmfdPinXS(ineigh, iz)%XSD(ig) / Pin(ineigh)%Center2SurfaceL(inbd)
          Dtil = mybeta * neighbeta/(mybeta + neighbeta) * smy
          PinXs%dtil(ibd, ig) = Dtil
!          IF(lXsGen) THEN
!            neighphi = CmfdPinXs(ineigh, iz)%phi(ig); phisum = neighphi + myphi(ig)
!            jfdm = Dtil * (myphi(ig) - neighphi); jnet = 0
!            dhat = (jfdm - jnet) / phisum
!            PinXs%dhat(ibd, ig) = Dhat
!          ENDIF
        ENDDO  !ENd of Group Sweep 
      ELSE     !Boundary 
        IF(ineigh .EQ. VoidCell) THEN
          neighbeta = 0.5_8
        ELSE
          neighbeta = 0
        ENDIF
        DO ig = 1, ng
          neighphi = 0; phisum = neighphi + myphi(ig)
          mybeta = CmfdPinXS(ixy, iz)%XSD(ig) / Pin(ixy)%Center2SurfaceL(ibd)
          Dtil = mybeta * neighbeta/(mybeta + neighbeta) * smy
          PinXs%dtil(ibd, ig) = Dtil
        ENDDO
      ENDIF
    ENDDO  !End of Neighborhood Sweep
    IF(lXsGen) CYCLE
    DO ibd= 1, nbd
      IF(.NOT. lxsFtn) THEN
        RefPinXS => DcplCmInfo(1, iRefPln0)%PinXs(ixy, iRefPln)
        PinXS%dhat(ibd, 1:ng) = RefPinXs%Dhat(ibd, 1:ng)
        PinXS%pdhat(ibd, 1:ng) = PinXS%dhat(ibd, 1:ng)
        !PinXS%dhat(ibd, 1:ng) = 0
        !PinXS%pdhat(ibd, 1:ng) = 0
      ELSE
        
        DhatDat = 0._8

        RefPinXS => DcplCmInfo(1, iRefPln0)%PinXs(ixy, iRefPln)
#ifdef DcplMemSave
        PinXS%dhat(ibd, 1:ng) = RefPinXs%Dcpl_Dhat(ibd, 1:ng, 1)
        PinXS%pdhat(ibd, 1:ng) = PinXS%dhat(ibd, 1:ng)          
#else
        PinXS%dhat(ibd, 1:ng) = RefPinXs%Dhat(ibd, 1:ng)
        PinXS%pdhat(ibd, 1:ng) = PinXS%dhat(ibd, 1:ng)        
#endif
      ENDIF
    ENDDO
    !CALL DcplDhatGen(PinXS, DhatData(1:nbd), nbd, ng, lXsFtn)
    NULLIFY(PinXs)
  ENDDO  !End of Radial Pin Sweep
ENDDO
CONTINUE
END SUBROUTINE

SUBROUTINE DcplDhatInt(Dtil, DHat, DhatDat, n, ng)
USE PARAM
USE UtilFunction        ,ONLY : Array2DSORT
USE XsUtil_mod          ,ONLY : LineIntPol
IMPLICIT NONE
REAL :: Dtil(ng)
REAL :: DHat(ng)
REAL :: DhatDat(5, 2, ng)
INTEGER :: ng, n
INTEGER :: i, j, ig, m, m0

DO ig = 1, ng
  m0 = n
  IF(Dtil(ig) .EQ. 0._8) CYCLE
  IF(DhatDat(3, 1, IG) .eq. 0._8) m0 = m0 -1
  CALL Array2DSORT(DhatDat(1:m0, 1, ig), DhatDat(1:3, 2, ig), m0, m, FALSE, FALSE, 1)
  Dhat(ig) = LineIntPol(Dtil(ig), m, DhatDat(1:m, 1, ig), DhatDat(1:m, 2, ig)) 
ENDDO
END SUBROUTINE

SUBROUTINE DcplEffXsSave(Core, Fxr, itemp, GroupInfo, DcplPE)
USE PARAM
USE TYPEDEF,    ONLY : CoreInfo_Type    ,FxrInfo_Type     &
                      ,GroupInfo_Type   ,PE_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_TYpe) :: DcplPE
INTEGER :: itemp

INTEGER :: nFxr, iz
INTEGER :: ifxr

nfxr = Core%nCoreFxr
iz = DcplPE%myzb

DO ifxr = 1, nFxr
  IF(ASSOCIATED(Fxr(ifxr, iz)%fresoa)) Fxr(ifxr, iz)%Dcpl_fresoa(:, itemp) = Fxr(ifxr, iz)%fresoa(:)
  IF(ASSOCIATED(Fxr(ifxr, iz)%fresoS)) Fxr(ifxr, iz)%Dcpl_fresoS(:, itemp) = Fxr(ifxr, iz)%fresoS(:)
  IF(ASSOCIATED(Fxr(ifxr, iz)%fresoF)) Fxr(ifxr, iz)%Dcpl_fresoF(:, itemp) = Fxr(ifxr, iz)%fresoF(:)
ENDDO
END SUBROUTINE
