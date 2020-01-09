#include <defines.h>

SUBROUTINE CopyThCondition_XsFTN(Core, Fxr, ThInfo, DcplFxr, DcplThInfo, iRefPln, iPln1, ipln2, ItrCntl)
USE PARAM
USE TypeDef,      ONLY : CoreInfo_Type      ,DcplInfo_Type       ,FxrInfo_Type     &
                        ,ThInfo_Type        ,Pin_Type            ,Cell_Type
USE Th_Mod,       ONLY : ThVar                        
USE itrcntl_mod,  ONLY : ItrCntl_TYPE
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(FxrInfo_Type) :: DcplFxr(:)
TYPE(ThInfo_Type) :: ThInfo
TYPE(ThInfo_Type) :: DcplThInfo
INTEGER :: iRefPln, iPln1, iPln2
TYPE(ItrCntl_Type) :: ItrCntl
LOGICAL :: lThConv

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
INTEGER :: nxy, FxrIdxSt, nLocalFxr, niso
INTEGER :: ixy, icel, ifxr, ipln
INTEGER :: i, j, k

REAL :: DelTsum, Tsum, DelT, Tavg1, Tavg2, Tavg0, ncell
REAL :: sum1, sum2, ptemp
nxy = Core%nxy
Pin => Core%Pin
CellInfo => Core%CellInfo
lThConv = .FALSE.
DelTsum = 0._8; Tsum =0
!Copy FXR Temperature.
DcplThInfo%Tin    = ThInfo%Tin; DcplThInfo%PowFa  = ThInfo%PowFa 
DcplThInfo%PExit  = ThInfo%PExit; DcplThInfo%PowLin = ThInfo%PowLin
DcplThInfo%MdotFa = ThInfo%MdotFa

DO ixy = 1, nxy
  FxrIdxSt = Pin(ixy)%FxrIdxSt; icel = Pin(ixy)%Cell(iRefPln)
  nLocalFxr = CellInfo(icel)%nFxr

  DO j = 1, nLocalFxr

    ifxr = FxrIdxSt + j - 1
    niso = Fxr(ifxr, ipln1)%niso
    
    DO i=1, niso
      sum1 = 0; sum2 = 0
      DO ipln = ipln1, ipln2
        sum1 = sum1 + Fxr(ifxr, ipln)%pnum(i) * Core%hz(ipln) * Fxr(ifxr, ipln)%area
        sum2 = sum2 + Core%hz(ipln) * Fxr(ifxr, ipln)%area
      ENDDO
      DcplFxr(ifxr)%idiso(i) = Fxr(ifxr, ipln1)%idiso(i)
      DcplFxr(ifxr)%pnum(i) = sum1/sum2
    ENDDO
    
    sum1 = 0; sum2 = 0
    DO ipln = ipln1, ipln2
      sum1 = sum1 + Fxr(ifxr, ipln)%temp * Core%hz(ipln) * Fxr(ifxr, ipln)%area
      sum2 = sum2 + Core%hz(ipln) * Fxr(ifxr, ipln)%area
    ENDDO
    ptemp = DcplFxr(ifxr)%temp
    DcplFxr(ifxr)%temp = sum1 / sum2
    
    DcplFxr(ifxr)%Dcpl_temp(1) = DcplFxr(ifxr)%Temp
    DcplFxr(ifxr)%Dcpl_pnum(1:niso, 1) = DcplFxr(ifxr)%pnum(1:niso)
    
    DelT = (DcplFxr(ifxr)%temp - ptemp) * Fxr(ifxr, iRefpln)%area
    DelTsum = DelTsum + DelT * DelT
    Tsum = Tsum + (DcplFxr(ifxr)%temp * DcplFxr(ifxr)%area)**2
      !
  ENDDO
ENDDO

DO ixy = 1, nxy
  !Coolant Temperature
  sum1 = 0; sum2 = 0
  DO ipln = ipln1, ipln2
    sum1 = sum1 + ThInfo%Tcool(ipln, ixy) * Core%hz(ipln)
    sum2 = sum2 + Core%hz(ipln)
  ENDDO
  DcplThInfo%Tcool(iRefPln, ixy) = sum1 / sum2
  
  !
  DcplThInfo%TcoolInOut(2, iRefPln, ixy) = ThInfo%TcoolInOut(2, ipln2, ixy)
  DcplThInfo%TcoolInOut(1, iRefPln, ixy) = ThInfo%TcoolInOut(1, ipln1, ixy)

  sum1 = 0; sum2 = 0
  DO ipln = ipln1, ipln2
    sum1 = sum1 +  THInfo%RelPower(ipln, ixy) * Core%hz(ipln)
    sum2 = sum2 + Core%hz(ipln)
  ENDDO
  DcplTHInfo%RelPower(iRefPln, ixy) = sum1 / sum2
!  DcplTHInfo%Tcool(ipln, ixy) = ThInfo%Tcool(ipln, ixy)
!  DcplTHInfo%TcoolInOut(1:2, ipln, ixy) = ThInfo%TcoolInOut(1:2, ipln, ixy)
!  DcplTHInfo%RelPower(ipln, ixy) = THInfo%RelPower(ipln, ixy)
ENDDO
#ifdef Tmodchg1
!DO ixy = 1, nxy
!  IF(.NOT. Pin(ixy)%lGT) CYCLE
!  Tavg1 = 0; Tavg2 = 0; ncell = 0
!  DO j = 1, Pin(ixy)%nBD
!    k = Pin(ixy)%NeighIdx(j)
!    IF(.NOT. Pin(k)%lfuel) CYCLE
!    ncell = ncell + 1
!    Tavg1 = Tavg1 + DcplTHInfo%TcoolInOut(1, iRefpln, k)
!    Tavg2 = Tavg2 + DcplTHInfo%TcoolInOut(2, iRefpln, k)
!  ENDDO
!  Tavg1 = Tavg1 / nCell; Tavg2 = Tavg2 / nCell
!  DcplThInfo%TcoolInOut(1:2, iRefPln, ixy) = (/Tavg1, Tavg2/)
!ENDDO
#endif
DO ixy = 1, nxy
  IF(Pin(ixy)%lFuel) THEN
    DcplThInfo%RhoU = ThInfo%CoolantTH(ixy)%rhou(0)
  ENDIF
ENDDO

DelTsum = SQRT(DelTsum / Tsum)
IF(DelTSum .LT. ItrCntl%ThConv) ItrCntl%lThConv = TRUE
END SUBROUTINE


SUBROUTINE CopyThCondition(Core, Fxr, ThInfo, DcplFxr, DcplThInfo, iPln, ItrCntl)
USE PARAM
USE TypeDef,      ONLY : CoreInfo_Type      ,DcplInfo_Type       ,FxrInfo_Type     &
                        ,ThInfo_Type        ,Pin_Type            ,Cell_Type

USE itrcntl_mod,  ONLY : ItrCntl_TYPE
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:), DcplFxr(:)
TYPE(ThInfo_Type) :: ThInfo
TYPE(ThInfo_Type) :: DcplThInfo
INTEGER :: iPln
TYPE(ItrCntl_Type) :: ItrCntl
LOGICAL :: lThConv

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
INTEGER :: nxy, FxrIdxSt, nLocalFxr, niso
INTEGER :: ixy, icel, ifxr
INTEGER :: i, j, k

REAL :: DelTsum, Tsum, DelT, Tavg1, Tavg2, Tavg0, ncell

nxy = Core%nxy
Pin => Core%Pin
CellInfo => Core%CellInfo
lThConv = .FALSE.
DelTsum = 0._8; Tsum =0
!Copy FXR Temperature.
DcplThInfo%Tin    = ThInfo%Tin; DcplThInfo%PowFa  = ThInfo%PowFa 
DcplThInfo%PExit  = ThInfo%PExit; DcplThInfo%PowLin = ThInfo%PowLin
DcplThInfo%MdotFa = ThInfo%MdotFa

DO ixy = 1, nxy
  FxrIdxSt = Pin(ixy)%FxrIdxSt; icel = Pin(ixy)%Cell(ipln)
  nLocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1
    niso = Fxr(ifxr)%niso
    !
    DelT = (DcplFxr(ifxr)%temp - Fxr(ifxr)%temp) * Fxr(ifxr)%area
    DelTsum = DelTsum + DelT * DelT
    Tsum = Tsum + (Fxr(ifxr)%temp * Fxr(ifxr)%area)**2
    !
    DcplFxr(ifxr)%temp = Fxr(ifxr)%temp
    DcplFxr(ifxr)%pnum(1:niso) = Fxr(ifxr)%pnum(1:niso)
    
    DcplFxr(ifxr)%Dcpl_temp(1) = Fxr(ifxr)%temp
    DcplFxr(ifxr)%Dcpl_pnum(1:niso, 1) = Fxr(ifxr)%pnum(1:niso)
  ENDDO
ENDDO

DO ixy = 1, nxy
  DcplTHInfo%Tcool(ipln, ixy) = ThInfo%Tcool(ipln, ixy)
  DcplTHInfo%TcoolInOut(1:2, ipln, ixy) = ThInfo%TcoolInOut(1:2, ipln, ixy)
  DcplTHInfo%RelPower(ipln, ixy) = THInfo%RelPower(ipln, ixy)
ENDDO
!
!DO ixy = 1, nxy
!  IF(.NOT. Pin(ixy)%lGT) CYCLE
!  Tavg1 = 0; Tavg2 = 0; ncell = 0
!  DO j = 1, Pin(ixy)%nBD
!    k = Pin(ixy)%NeighIdx(j)
!    IF(.NOT. Pin(k)%lfuel) CYCLE
!    ncell = ncell + 1
!    Tavg1 = Tavg1 + ThInfo%TcoolInOut(1, ipln, k)
!    Tavg2 = Tavg2 + ThInfo%TcoolInOut(1, ipln, k)
!  ENDDO
!  Tavg1 = Tavg1 / nCell; Tavg2 = Tavg2 / nCell
!  DcplThInfo%TcoolInOut(1:2, ipln, ixy) = (/Tavg1, Tavg2/)
!ENDDO

DO ixy = 1, nxy
  IF(Pin(ixy)%lFuel) THEN
    DcplThInfo%RhoU = ThInfo%CoolantTH(ixy)%rhou(0)
  ENDIF
ENDDO

DelTsum = SQRT(DelTsum / Tsum)
IF(DelTSum .LT. ItrCntl%ThConv) ItrCntl%lThConv = TRUE
END SUBROUTINE


SUBROUTINE SetXsGenThCondition(Core, FmInfo, DcplInfo, THInfo, DcplThInfo, lXsFtn, DcplItrCntl, PE)
USE PARAM
USE TypeDef,      ONLY : CoreInfo_Type     ,FmInfo_Type     ,DcplInfo_Type    &
                        ,FxrInfo_Type      ,ThInfo_Type     ,PE_TYPE
USE DcplTh_Mod,   ONLY : CopyThCondition   ,DcplModTChg     ,DcplFuelTChg     &
                        ,CopyThCondition_XsFtn
USE itrcntl_mod,  ONLY : ItrCntl_TYPE
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(DcplInfo_Type) :: DcplInfo
TYPE(ItrCntl_Type) :: DcplItrCntl(100)
TYPE(PE_TYPE) :: PE
TYPE(ThInfo_Type) :: ThInfo
LOGICAL :: lXsFtn

TYPE(ThInfo_Type), POINTER :: DcplThInfo(:)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :), DcplFxr(:, :)
INTEGER :: iRefPln, ipln, ipln1, ipln2
INTEGER :: myRefPlnBeg, myRefPlnEnd
INTEGER :: nRefPln
LOGICAL :: lThConv
nRefPln = DcplInfo%nRefPln
#ifndef MPI_ENV
myRefPlnBeg = 1; myRefPlnEnd = nRefPln
#else
myRefPlnBeg = PE%myRefPlnBeg; myRefPlnEnd = PE%myRefPlnEnd
#endif
DO iRefPln = myRefPlnBeg, myRefPlnEnd
  ipln = DcplInfo%RefPln(iRefPln)
  ipln1 = DcplInfo%RefPln_Map(1, iRefPln)
  ipln2 = DcplInfo%RefPln_Map(2, iRefPln)
  Fxr => FmInfo%Fxr
  DcplFxr => DcplInfo%DcplFmInfo(1, iRefPln)%Fxr
  IF(lXsFtn) THEN
    CALL CopyThCondition_XsFtn(Core, Fxr, ThInfo, DcplFxr(:, ipln), DcplTHInfo(iRefpln), ipln, ipln1, ipln2, DcplItrCntl(iRefPln))
  ELSE
    CALL CopyThCondition(Core, Fxr(:, ipln), ThInfo, DcplFxr(:, ipln), DcplTHInfo(iRefpln), ipln, DcplItrCntl(iRefPln))
  ENDIF
ENDDO

END SUBROUTINE

SUBROUTINE DcplThUpdate(Core, CmInfo, FmInfo, ThInfo, DcplInfo, DcplThInfo, Eigv, ng, lUpdt, GroupInfo,nTracerCntl, DcplItrCntl, PE)
USE PARAM
USE TypeDef,          ONLY : CoreInfo_Type      ,CmInfo_Type      ,FmInfo_Type            &
                            ,ThInfo_Type        ,DcplInfo_Type    ,PE_Type                &
                            ,GroupInfo_Type
USE CNTL,             ONLY : nTracerCntl_Type
USE Th_Mod,           ONLY : SteadyStateTH      ,THFeedBack
USE DcplTH_Mod,       ONLY : SetXsGenThCondition
USE itrcntl_mod,      ONLY : ItrCntl_TYPE  
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CMInfo_Type) :: CmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(ThInfo_Type), POINTER :: DcplThInfo(:)
TYPE(DcplInfo_Type) :: DcplInfo
TYPE(PE_TYPE) :: PE
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(ItrCntl_Type) :: DcplItrCntl(100)

REAL :: Eigv
INTEGER :: ng
LOGICAL :: lUpdt

CALL SteadyStateTH(Core, CmInfo, FmInfo, ThInfo, Eigv, ng, GroupInfo, nTracerCntl, PE)
CALL THFeedBack(Core, CmInfo, FmInfo, ThInfo, nTracerCntl, PE)    
IF(lUpdt) CALL SetXsGenThCondition(Core, FmInfo, DcplInfo, ThInfo, DcplThInfo, nTracerCntl%lXsFtn, DcplItrCntl, PE)
END SUBROUTINE

SUBROUTINE DcplModTChg(Core, Fxr, ThInfo, iRefPln, ipln, imod, frac)
USE PARAM
USE TYPEDEF,     ONLY : CoreInfo_Type      ,FxrInfo_Type      ,ThInfo_Type       &
                       ,Pin_Type           ,Cell_Type
USE SteamTBL_mod, ONLY : steamtbl
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(ThInfo_Type) :: ThInfo


INTEGER :: iRefPln, ipln, imod
REAL :: Frac

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
REAL, POINTER :: Tcool(:, :)

REAL :: Tc, TcInOut(2), TcNew, DelT, PExit
REAL :: wt, wh, hin, wrho, wvin, wxin, wbetain, wkapain, wcpin
REAL ::  nh2o, nh
REAL :: MAX_DELT
REAL, PARAMETER :: DelTmin = 5._8
INTEGER :: ixy, ifxr, icel
INTEGER :: nxy, nfxr, niso
INTEGER :: FxrIdxSt, nLocalFxr
INTEGER :: i, j, k
INTEGER :: CoolReg(2)

nxy = Core%nxy;  nFxr = Core%nCoreFxr
CellInfo => Core%CellInfo
Pin => Core%Pin; Tcool => THInfo%Tcool

PEXIT = ThInfo%PExit 

!ENDDO
MAX_DELT = 0
DO ixy = 1, nxy
  icel = Pin(ixy)%Cell(ipln)
  FxrIdxSt = Pin(ixy)%FxrIdxSt; nLocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1  
    niso = Fxr(ifxr, ipln)%niso
    Fxr(ifxr, ipln)%pnum(1:niso) = Fxr(ifxr, ipln)%Dcpl_Pnum(1:niso, 1)
    Fxr(ifxr, ipln)%temp = Fxr(ifxr, ipln)%dcpl_temp(1)
  ENDDO
ENDDO

DO ixy = 1, nxy
  IF(Pin(ixy)%lRadRef) CYCLE
  icel = Pin(ixy)%Cell(ipln)
  FxrIdxSt = Pin(ixy)%FxrIdxSt; nLocalFxr = CellInfo(icel)%nFxr
  IF(Core%lFuelPlane(ipln)) THEN
    Tc = ThInfo%Tcool(ipln, ixy)
    TcInOut = ThInfo%TcoolInOut(1:2, ipln, ixy)
    DelT = frac * Abs(TcInout(imod) - Tc)
    DelT = MAX(DelT,DelTmin) 
    Tc = Tc + DelT
    Max_DelT = MAX(Max_DelT, DelT)
    !Tc = Tc + 5
  ELSE
    IF(imod .EQ. 2) Tc = ThInfo%Tcool(ipln, ixy) + 10._8
    IF(imod .EQ. 1) Tc = ThInfo%Tcool(ipln, ixy) - 10._8
  ENDIF
  !Tc = ThInfo%Tcool(ipln, ixy)
  wh = 0; wrho = 0; wvin = 0; wxin = 0
  wbetain = 0; wkapain = 0; wcpin = 0
  TC = TC + CKELVIN
  CALL SteamTbl(TRUE, pexit, Tc, wh, wrho, wvin, wxin, wbetain, wkapain, wcpin)
  wrho = wrho * epsm3;
  nh2o = wrho / awh2o * avogadro; nh = 2. * nh2o
  
  IF(CellInfo(icel)%lFuel) THEN
    CoolReg = CellInfo(icel)%ThCell%CoolReg
    DO j = CoolReg(1), CoolReg(2)
      ifxr = FxrIdxSt + j - 1; Fxr(ifxr, ipln)%Temp = TC
      CALL H2ODensityUpdate(nH2O, nH, Fxr(ifxr, ipln)%niso, Fxr(ifxr, ipln)%IdIso, Fxr(ifxr, ipln)%pnum)
    ENDDO    
  ELSE
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j - 1
      Fxr(ifxr, ipln)%temp = TC
      IF(Fxr(ifxr, ipln)%lh2o) THEN
        CALL H2ODensityUpdate(nH2O, nH, Fxr(ifxr, ipln)%niso, Fxr(ifxr, ipln)%IdIso, Fxr(ifxr, ipln)%pnum)
      ENDIF
    ENDDO  
  ENDIF
ENDDO

END SUBROUTINE

SUBROUTINE DcplFuelTChg(Core, Fxr, ThInfo, iRefPln, ipln, frac)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type     ,FxrInfo_Type     ,ThInfo_Type             &
                            ,Pin_Type          ,Cell_Type        ,THCell_Type
USE TH_Mod,           ONLY : ThOpt             ,ThVar
USE FuelProperty_Mod, ONLY : fhtcoef
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(ThInfo_Type) :: ThInfo
INTEGER :: ipln, iRefPln

REAL :: Frac

TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(THCell_Type), POINTER :: ThCell
REAL, POINTER :: RelPower(: ,:)
REAL :: tfvol(100), tfuel(100)
REAL :: tcool, TempK, htcoef, qf
REAL :: PowLin, PowLv, RelPinPw, frac0, qprime, rhou
REAL :: deq, afp, fracdf

INTEGER :: CLDReg(2), FuelReg(2)
INTEGER :: ixy, icel, ifxr
INTEGER :: nxy, FxrIdxSt, nLocalFxr, nFxr, niso, nr, npr2
INTEGER :: i, j, k

powlin = ThInfo%PowLin; PowLv = ThInfo%PowLv
RelPower => ThInfo%RelPower
frac0 = (1.+frac)

nxy = Core%nxy; nFxr = Core%nCoreFxr
CellInfo => Core%CellInfo
Pin => Core%Pin;
IF(.NOT. Core%lFuelPlane(ipln)) RETURN

DO ixy = 1, nxy
  icel = Pin(ixy)%Cell(ipln)
  FxrIdxSt = Pin(ixy)%FxrIdxSt; nLocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1  
    niso = Fxr(ifxr, ipln)%niso
    Fxr(ifxr, ipln)%pnum(1:niso) = Fxr(ifxr, ipln)%Dcpl_Pnum(1:niso, 1)
    Fxr(ifxr, ipln)%temp = Fxr(ifxr, ipln)%dcpl_temp(1)
  ENDDO
ENDDO

rhou = ThInfo%Rhou
nr = ThVar%npr5; Deq = ThVar%Deq
npr2 = ThVar%npr2 
afp = ThVar%afp; fracdf = ThVar%fracdf
DO ixy = 1, nxy
  IF(.NOT. Pin(ixy)%lFuel) CYCLE
  !Obtain the Fuel Temp Distribution
  RelPinPw = RelPower(ipln, ixy) * frac0
  tcool = ThInfo%Tcool(ipln, ixy)
  qprime = powlv * powlin * RelPinPw
  qf = fracdf * qprime / afp
  htcoef = fhtcoef(tcool, deq, Rhou)
  CALL tfcalss(tfvol(1:nr-3), tfuel(1:nr), tcool, htcoef, qf, nr, ThVar, ThOpt, .FALSE.)
  
  
  icel = Pin(ixy)%Cell(ipln)
  ThCell => CellInfo(icel)%ThCell 
  FxrIdxSt = Pin(ixy)%FxrIdxSt; nLocalFxr = CellInfo(icel)%nFxr
  CLDReg = CellInfo(icel)%ThCell%CldReg; FuelReg = CellInfo(icel)%ThCell%FuelReg
  DO j = CLDReg(1), CLDReg(2)
    ifxr = FxrIdxSt + j - 1
    Fxr(ifxr, ipln)%Temp = Tfvol(npr2)
  ENDDO
  
  DO j = FuelReg(1), FuelReg(2)
    ifxr = FxrIdxSt + j - 1
    TempK = 0
    DO k =  ThCell%FuelMapping(1, j), ThCell%FuelMapping(2, j)
      TempK = TempK + ThCell%Frac(k, j) * tfvol(k)
    ENDDO
    Fxr(ifxr, ipln)%Temp = TempK
  ENDDO
ENDDO
END SUBROUTINE