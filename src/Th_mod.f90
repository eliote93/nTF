MODULE TH_Mod
USE PARAM
USE TypeDef,    ONLY : ThOpt_Type,  ThVar_Type
IMPLICIT NONE

TYPE(THOpt_Type) :: ThOpt
TYPE(THVar_Type) :: ThVar
#ifndef __GFORTRAN__
! conductivity, w/m-C, t in K
DATA ThOPT%kFUelCorrelation  /1.05_8, 0._8, 0._8, 0._8, 2150._8, 73.15_8/
DATA ThOPT%KCladCorrelation  /7.51_8, 2.09e-2_8, -1.45e-5_8, 7.67e-9_8/
!  volumetric heat capacity, in J/m^3-C, t in K
DATA ThOpt%CpFuelCorrelation /1668768.6_8, 3123.6716_8, -2.4584262_8, 0.000658459_8/ !rhof=10.282
DATA ThOpt%CpCladCorrelation /1666764.0_8, 757.284_8, 0., 0./ !rhoclad=6.6 g/cc
#endif

INTERFACE
  SUBROUTINE AllocTH()
  END SUBROUTINE

  SUBROUTINE Cal_RefFuelTemp(RefFuelTemp, Core, Fxr, nTracerCntl, PE)
  USE PARAM
  USE TYPEDEF,       ONLY : coreinfo_type,     Fxrinfo_type,       PE_TYPE
  USE CNTL,          ONLY : nTracerCntl_Type
  IMPLICIT NONE
  TYPE(CoreInfo_Type) :: Core
  TYPE(Fxrinfo_type), POINTER :: Fxr(:, :)
  TYPE(PE_TYPE) :: PE
  TYPE(nTracerCntl_Type) :: nTracerCntl
  REAL, POINTER :: RefFuelTemp(:)
  END SUBROUTINE

  SUBROUTINE THFeedBack(Core, CmInfo, FmInfo, ThInfo, nTRACERCntl, PE)
  USE PARAM
  USE TYPEDEF,        ONLY : CoreInfo_Type,        CMInfo_Type,      FmInfo_Type,     &
                             ThInfo_Type,          ThVar_TYpe,       ThOpt_Type,      &
                             PE_Type
  USE CNTL,           ONLY : nTracerCntl_Type
  IMPLICIT NONE

  TYPE(CoreInfo_Type) :: Core
  TYPE(CMInfo_Type) :: CmInfo
  TYPE(FMInfo_Type) :: FmInfo
  TYPE(ThInfo_Type) :: ThInfo
  TYPE(nTracerCntl_Type) :: nTRACERCntl
  TYPE(PE_Type) :: PE

  END SUBROUTINE

  SUBROUTINE TransientCoolantTH(qeff, qeffd, Pexit, CoolantTH, iz, ThVar, ThOpt)
  USE PARAM
  USE TYPEDEF,          ONLY : CoolantTH_Type,       ThVar_Type,       ThOpt_Type
  IMPLICIT NONE
  TYPE(CoolantTH_Type) :: CoolantTH
  TYPE(ThVar_Type) :: ThVar
  TYPE(ThOPT_Type) :: ThOpt
  REAL :: qeff, qeffd, PEXIT
  INTEGER :: iz
  END SUBROUTINE

  SUBROUTINE SteadyCoolantTH(powlin, plevel, Pexit, Tout, RelPW, CoolantTH, ThVar,ThOpt, PE)
  USE PARAM
  USE TYPEDEF,      ONLY : CoolantTH_Type,       ThVar_Type,     &
                           ThOpt_Type,           PE_Type
  IMPLICIT NONE
  TYPE(CoolantTH_Type) :: CoolantTH        ! One Channel Coolant TH information
  TYPE(ThVar_Type) :: ThVar                !
  TYPE(ThOPT_Type) :: ThOpt                !
  TYPE(PE_Type) :: PE                      !
  REAL :: PowLin, Plevel, PEXIT, Tout      !
  REAL :: RelPW(:)

  END SUBROUTINE

  SUBROUTINE Tfcaltr(tfvol, Tfuel, Tfueld, Tcool, htCoef, Qf, Qfd, nr, ThVar, ThOpt)
  USE PARAM
  USE TYPEDEF,           ONLY : ThVar_Type,     ThOpt_Type
  IMPLICIT NONE
  REAL :: tfvol(nr-3)
  REAL :: Tfuel(nr), Tfueld(nr), tcool, htcoef, qf, Qfd
  INTEGER :: nr
  TYPE(ThVar_Type) :: ThVar
  TYPE(ThOPT_Type) :: ThOPT
  END SUBROUTINE

  SUBROUTINE SteadyFuelConduction(powlin, plevel, Tfmax, RelPw, PwShape, BurnUp, RhoU,FuelTH, ThVar, ThOpt, nTracerCntl, PE)
  USE PARAM
  USE TYPEDEF,          ONLY :  FuelTH_Type,     ThVar_Type,     ThOpt_Type,          &
                                PE_TYPE
  USE CNTL,             ONLY : nTracerCntl_Type
  IMPLICIT NONE
  TYPE(FuelTH_Type):: FuelTH        ! One Channel Coolant TH information
  TYPE(ThVar_Type) :: ThVar                !
  TYPE(ThOPT_Type) :: ThOpt                !
  TYPE(nTracerCntl_Type) :: nTracerCntl
  TYPE(PE_Type) :: PE                      !
  REAL :: PowLin,  Plevel, Tfmax
  REAL :: RelPW(:)
  REAL, POINTER :: RhoU(:), PwShape(:, :), BurnUp(:, :)

  END SUBROUTINE

  SUBROUTINE SteadyStateTH(Core, CmInfo, FmInfo, ThInfo, Eigv, ng, GroupInfo, nTRACERCntl, PE)
  USE PARAM
  USE TYPEDEF,        ONLY : CoreInfo_Type,        CMInfo_Type,      FmInfo_Type,     &
                             ThInfo_Type,          GroupInfo_Type,   PE_Type
  USE CNTL,           ONLY : nTracerCntl_Type
  IMPLICIT NONE
  TYPE(CoreInfo_Type) :: Core
  TYPE(CMInfo_Type) :: CmInfo
  TYPE(FMInfo_Type) :: FmInfo
  TYPE(ThInfo_Type) :: ThInfo
  TYPE(GroupInfo_Type) :: GroupInfo
  TYPE(nTracerCntl_Type) :: nTRACERCntl
  TYPE(PE_Type) :: PE
  REAL :: Eigv
  INTEGER :: ng

  END SUBROUTINE

  SUBROUTINE SimpleTH(powlin, plevel, Tout, Tfmax, RelPW, Profile, FuelTh, CoolantTH, Thvar, ThOpt, PE)
  USE PARAM
  USE TYPEDEF,           ONLY : CoolantTh_Type,    FuelTh_Type,   ThVar_Type,       &
                                CoreInfo_Type,     ThOpt_Type,        PE_TYPE
  IMPLICIT NONE
  REAL :: PowLin, Plevel, Tout, Tfmax
  REAL, POINTER :: Profile(:, :)
  REAL :: RelPW(:)
  TYPE(FuelTh_Type) :: FuelTh
  TYPE(CoolantTh_Type) :: CoolantTH
  TYPE(Thvar_Type) :: ThVar
  TYPE(THOpt_Type) :: ThOPT
  TYPE(PE_TYPE) :: PE

  END SUBROUTINE

  SUBROUTINE TransientTH(Core, CmInfo, FmInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
  USE PARAM
  USE TYPEDEF,                    ONLY : CoreInfo_Type,              FmInfo_Type,               CmInfo_Type,          &
                                         ThInfo_Type,                GroupInfo_Type,            TranCntl_Type,        &
                                         PE_TYPE,                                                                     &
                                         FuelTH_Type,                CoolantTH_Type,            Pin_Type,             &
                                         FxrInfo_Type,               PinXS_Type
  USE CNTL,                       ONLY : nTracerCntl_Type
  IMPLICIT NONE
  TYPE(CoreInfo_Type) :: Core
  TYPE(CmInfo_Type) :: CmInfo
  TYPE(FmInfo_Type) :: FmInfo
  TYPE(ThInfo_Type) :: ThInfo
  TYPE(GroupInfo_Type) :: GroupInfo
  TYPE(TranCntl_Type) :: TranCntl
  TYPE(nTracerCntl_Type) :: nTracerCntl
  TYPE(PE_Type) :: PE

END SUBROUTINE

  SUBROUTINE SaveTranThSol(Core, ThInfo, TranCntl, nTracerCntl, PE)
  USE PARAM
  USE TYPEDEF,                     ONLY : CoreInfo_Type,           ThInfo_Type,             TranCntl_Type,      &
                                          PE_Type
  USE CNTL,                        ONLY : nTracerCntl_Type
  IMPLICIT NONE
  TYPE(CoreInfo_Type) :: Core
  TYPE(ThInfo_Type) :: ThInfo
  TYPE(TranCntl_Type) :: TranCntl
  TYPE(nTracerCntl_Type) :: nTracerCntl
  TYPE(PE_Type) :: PE
  END SUBROUTINE

  SUBROUTINE CalRelPower(Core, CmInfo, RelPower, ng, nTracerCntl, PE, lGather)

  !Update Normalized Power for T/H
  USE PARAM
  USE TypeDef,        ONLY : CoreInfo_Type,            CmInfo_Type,         PE_Type,          &
                             PinXS_Type,               Pin_Type
  USE CNTL,           ONLY : nTracerCntl_Type
  IMPLICIT NONE

  TYPE(CoreInfo_Type) :: Core
  TYPE(CmInfo_Type) :: CmInfo
  TYPE(nTracerCntl_Type) :: nTracerCntl
  TYPE(PE_TYPE) :: PE
  INTEGER :: ng
  REAL, POINTER :: RelPower(:, :)
   LOGICAL :: lGather
  END SUBROUTINE

  SUBROUTINE SubCellPowerProfile(Core, FmInfo, Profile, ixy, npr, nz, ng, GroupInfo, PE)
  USE PARAM
  USE TYPEDEF,      ONLY : CoreInfo_Type,        FmInfo_Type,    PE_TYPE,          &
                           GroupInfo_Type
  IMPLICIT NONE
  TYPE(CoreInfo_Type) :: Core
  TYPE(FmInfo_Type) :: FmInfo
  TYPE(GroupInfo_Type) :: GroupInfo
  TYPE(PE_Type) :: PE
  REAL, POINTER :: Profile(:, :)
  INTEGER :: ixy, ng, nz, npr

  END SUBROUTINE

  SUBROUTINE Grp_RelPower(Core, CmInfo, RelPower, ng, nTracerCntl, PE, lGather)

  !Update Normalized Power for T/H
  USE PARAM
  USE TypeDef,        ONLY : CoreInfo_Type,            CmInfo_Type,         PE_Type,          &
                             PinXS_Type,               Pin_Type
  USE CNTL,           ONLY : nTracerCntl_Type
  IMPLICIT NONE

  TYPE(CoreInfo_Type) :: Core
  TYPE(CmInfo_Type) :: CmInfo
  TYPE(nTracerCntl_Type) :: nTracerCntl
  TYPE(PE_TYPE) :: PE
  INTEGER :: ng
  REAL, POINTER :: RelPower(:, :)
   LOGICAL :: lGather
  END SUBROUTINE


  SUBROUTINE GatherRelPower(Core, RelPow, PE)
  USE PARAM
  USE TYPEDEF,  ONLY : CoreInfo_Type, PE_TYPE
  USE MPIComm_mod, ONLY : BCAST
  IMPLICIT NONE
  TYPE(CoreInfo_Type) :: Core
  TYPE(PE_TYPE) :: PE
  REAL, POINTER :: RelPow(:, :)

  END SUBROUTINE

  SUBROUTINE SetInletTHInfo(ThInfo, THVarIn, Pin, nTracerCntl, nxy, nz)
  USE TYPEDEF,      ONLY : THInfo_Type,           THVar_Type,     Pin_Type
  USE CNTL,         ONLY : nTracerCntl_Type

  IMPLICIT NONE

  TYPE(THInfo_Type) :: ThInfo
  TYPE(THVar_Type) :: ThVarIn
  TYPE(Pin_Type) :: Pin(nxy)
  TYPE(nTracerCntl_Type) :: nTracerCntl
  INTEGER :: nz, nxy

  END SUBROUTINE

  SUBROUTINE InitCoolantVar(Tink, pexit, MdotFa, ACF, CoolantTh, nrpallet, nxy, nz)
  USE PARAM
  USE TYPEDEF,      ONLY : THInfo_TYPE,    CoolantTH_Type

  IMPLICIT NONE
  REAL :: TinK, Pexit, MdotFa, ACF
  TYPE(CoolantTh_Type) :: CoolantTH(nxy)
  INTEGER :: nrpallet, nxy, nz

  END SUBROUTINE

  SUBROUTINE AllocCoolantTH(Pin, ThInfo, nrpallet, nxy, nz)
  USE PARAM
  USE TYPEDEF,      ONLY : THInfo_Type, Pin_Type
  USE ALLOCS
  IMPLICIT NONE

  TYPE(Pin_Type) :: Pin(nxy)
  TYPE(THInfo_Type) :: ThInfo
  INTEGER :: nrpallet, nxy, nz

  ENDSUBROUTINE

  SUBROUTINE InitFuelThVar(Tin, FuelTh, nrpallet, nxy, nz)
  USE PARAM
  USE TYPEDEF,    ONLY : FuelTH_Type
  IMPLICIT NONE
  INTEGER :: nrpallet, nxy, nz
  REAL :: Tin
  TYPE(FuelTH_Type) :: FuelTH(nxy)
  END SUBROUTINE

  SUBROUTINE AllocFuelTh(Pin, ThInfo, nrpallet, nxy, nz)
  USE PARAM
  USE TYPEDEF,      ONLY : ThInfo_Type, Pin_Type
  USE ALLOCS
  IMPLICIT NONE
  TYPE(Pin_Type) :: Pin(nxy)
  TYPE(ThInfo_Type) :: ThInfo
  INTEGER :: nrpallet, nxy, nz

  END SUBROUTINE


  FUNCTION CalPlnFuelTemp(Core, Fxr, ipln)
  USE PARAM
  USE TYPEDEF,          ONLY : CoreInfo_Type       ,FxrInfo_Type
  IMPLICIT NONE
  TYPE(CoreInfo_Type) :: Core
  TYPE(FxrInfo_Type), POINTER :: Fxr(:)
  INTEGER :: ipln

  REAL :: CalPlnFuelTemp

  END FUNCTION

  FUNCTION CalPlnModTemp(Core, Fxr, ipln)
  USE PARAM
  USE TYPEDEF,          ONLY : CoreInfo_Type       ,FxrInfo_Type
  IMPLICIT NONE
  TYPE(CoreInfo_Type) :: Core
  TYPE(FxrInfo_Type), POINTER :: Fxr(:)
  INTEGER :: ipln

  REAL :: CalPlnModTemp

  END FUNCTION

  FUNCTION CalPlnTemp(Core, Fxr, ipln)
  USE PARAM
  USE TYPEDEF,          ONLY : CoreInfo_Type       ,FxrInfo_Type
  IMPLICIT NONE
  TYPE(CoreInfo_Type) :: Core
  TYPE(FxrInfo_Type), POINTER :: Fxr(:)
  INTEGER :: ipln

  REAL :: CalPlnTemp
  END FUNCTION

  FUNCTION GetPinFuelTemp(Core, Fxr, ipln, ipin)
  USE PARAM
  USE TYPEDEF,       ONLY : CoreInfo_Type       ,FxrInfo_Type            &
                           ,Pin_Type            ,Cell_Type
  TYPE(CoreInfo_Type) :: Core
  TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
  INTEGER :: ipln, ipin
  REAL :: GetPinFuelTemp

  END FUNCTION

  FUNCTION GetPinModTemp(Core, Fxr, ipln, ipin)
  USE PARAM
  USE TYPEDEF,       ONLY : CoreInfo_Type       ,FxrInfo_Type            &
                           ,Pin_Type            ,Cell_Type
  TYPE(CoreInfo_Type) :: Core
  TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
  INTEGER :: ipln, ipin
  REAL :: GetPinModTemp

  END FUNCTION

  FUNCTION GetPinTemp(Core, Fxr, ipln, ipin)
  USE PARAM
  USE TYPEDEF,       ONLY : CoreInfo_Type       ,FxrInfo_Type            &
                           ,Pin_Type            ,Cell_Type
  TYPE(CoreInfo_Type) :: Core
  TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
  INTEGER :: ipln, ipin
  REAL :: GetPinTemp

  END FUNCTION

  SUBROUTINE SetPwShape(Pwshape, Core, Fxr, Power, ixy, nz, nr, Thvar, ThOpt, PE)
  USE PARAM
  USE TYPEDEF, ONLY : CoreInfo_Type,  FxrInfo_Type,     ThVar_Type,     ThOpt_Type,   &
                      PE_TYPE
  IMPLICIT NONE
  INTEGER :: ixy, nz, nr
  REAL :: Pwshape(nr, nz)        !  Heat source Shape
  TYPE(CoreInfo_Type) :: Core   !  Core Geometry Infomation
  TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
  REAL, POINTER :: Power(:, :)
  TYPE(ThVar_Type) :: ThVar     !  Th related parameters
  TYPE(ThOpt_TYPE) :: ThOpt     !  Th Option
  TYPE(PE_TYPE) :: PE


  END SUBROUTINE


  SUBROUTINE SetBurnupDist(BurnupDist, Core, FmInfo, ixy, nz, nr, Thvar, ThOpt, PE)
  USE PARAM
  USE TYPEDEF, ONLY : CoreInfo_Type,  FmInfo_Type,     ThVar_Type,     ThOpt_Type,     &
                      PE_TYPE
  IMPLICIT NONE
  INTEGER :: ixy, nz, nr
  REAL :: BurnupDist(0:nr, nz)        !  Heat source Shape
  TYPE(CoreInfo_Type) :: Core   !  Core Geometry Infomation
  TYPE(FmInfo_Type) :: FmInfo   !  Fine Mesh Information
  TYPE(ThVar_Type) :: ThVar     !  Th related parameters
  TYPE(ThOpt_TYPE) :: ThOpt     !  Th Option
  TYPE(PE_TYPE) :: PE
  END SUBROUTINE

  SUBROUTINE TranRadFuelCondFDM(tfvol, TFuel, TFuelD, Tcool, Tcoold, htcoef, htcoefd, qf, qfd, QShape, Qshaped, nr, ThVar, ThOpt)
  USE PARAM
  USE TYPEDEF,           ONLY : ThVar_Type,     ThOpt_Type
  IMPLICIT NONE

  REAL :: tfvol(nr-3)
  REAL :: Tfuel(nr), Tfueld(nr), QShape(nr), QShaped(nr)
  REAL :: tcoold, tcool, htcoef, htcoefd, qf, qfd
  INTEGER :: nr
  TYPE(ThVar_Type) :: ThVar
  TYPE(ThOPT_Type) :: ThOPT

  END SUBROUTINE

  SUBROUTINE RadFuelCondFDM(tfvol, Tfuel, tcool, htcoef, qf, qshape, nr, ThVar, ThOpt, lMox)
  USE PARAM
  USE TYPEDEF,           ONLY : ThVar_Type,     ThOpt_Type
  IMPLICIT NONE
  REAL :: tfvol(nr-3), Qshape(nr)
  REAL :: Tfuel(nr), tcool, htcoef, qf
  INTEGER :: nr
  TYPE(ThVar_Type) :: ThVar
  TYPE(ThOPT_Type) :: ThOPT
  LOGICAL :: lMox

  END SUBROUTINE

  SUBROUTINE SetFuelCondLs(Diag, L, U, b, x, tcool, htcoef, qf, qshape, nr, ThVar, ThOpt, lMox)
  USE PARAM
  USE TYPEDEF,           ONLY : ThVar_Type,     ThOpt_Type
  IMPLICIT NONE
  REAL :: diag(100), L(100), U(100), x(100), b(100), qshape(nr)
  REAL :: tcool, htcoef, qf
  INTEGER :: nr
  TYPE(ThVar_Type) :: ThVar
  TYPE(ThOPT_Type) :: ThOPT
  LOGICAL :: lMox

  END SUBROUTINE

  FUNCTION CalWallTemp(X, Tcool, htcoef, ThVar, ThOpt)
  USE PARAM
  USE TYPEDEF,           ONLY : ThVar_Type, ThOpt_Type
  IMPLICIT NONE
  REAL :: CalWallTemp
  REAL :: x(100)
  REAL :: tcool, htcoef
  TYPE(ThVar_Type) :: ThVar
  TYPE(ThOpt_Type) :: ThOpt

  END FUNCTION



  FUNCTION CalGapTemp(X, ThVar, ThOpt)
  USE PARAM
  USE TYPEDEF,           ONLY : ThVar_Type, ThOpt_Type
  IMPLICIT NONE
  REAL :: x(100)
  TYPE(ThVar_Type) :: ThVar
  TYPE(ThOpt_Type) :: ThOpt
  REAL :: CalGapTemp(2)

  END FUNCTION

  FUNCTION CalAvgFuelTemp(X, ThVar, ThOpt)
  USE PARAM
  USE TYPEDEF,            ONLY : ThVar_Type,      ThOpt_Type
  IMPLICIT NONE
  REAL :: X(100)
  TYPE(ThVar_Type) :: ThVar
  TYPE(ThOpt_Type) :: ThOpt
  REAL :: CalAvgFuelTemp

  END FUNCTION


  SUBROUTINE ReadUserDefTemp(Core, ThInfo, PE)
  USE TYPEDEF,        ONLY : CoreInfo_Type,       ThInfo_Type,          PE_Type
  TYPE(CoreInfo_Type) :: Core
  TYPE(ThInfo_Type) :: ThInfo
  TYPE(PE_Type) :: PE

  END SUBROUTINE

  SUBROUTINE SetUserDefTH(Core, CmInfo, FmInfo, ThInfo, nTRACERCntl, PE)
  USE PARAM
  USE TYPEDEF,        ONLY : CoreInfo_Type,        CMInfo_Type,      FmInfo_Type,     &
                             ThInfo_Type,          GroupInfo_Type,   PE_Type
  USE CNTL,           ONLY : nTracerCntl_Type
  IMPLICIT NONE
  TYPE(CoreInfo_Type) :: Core
  TYPE(CmInfo_Type) :: CmInfo
  TYPE(FmInfo_Type) :: FmInfo
  TYPE(ThInfo_Type) :: ThInfo
  TYPE(nTracerCntl_Type) :: nTracerCntl
  TYPE(PE_Type) :: PE

  END SUBROUTINE

  SUBROUTINE ChangeFlowRate(Core, ThInfo, Thvar, FlowRateLevel)
  USE PARAM
  USE Typedef,          ONLY : CoreInfo_Type,     ThInfo_Type,      Thvar_Type

  END SUBROUTINE

END INTERFACE
END MODULE

MODULE FuelProperty_Mod
USE PARAM
USE TH_MOD, ONLY : ThOpt
IMPLICIT NONE

CONTAINS

FUNCTION CondFuel_FRAPCON(T, Burnup)
IMPLICIT NONE
REAL :: t , burnup, CondFuel_FRAPCON
REAL :: a, b, burn, fbu, f, gbu, h, l, k
!Thermal conductivity corelation taken from FRAPCON-3
a = 4.52e-2_8 ;
b = 2.46e-4_8;
burn = burnup;
fbu = 1.87e-3_8*burn;
f=(1-0.9_8*exp(-0.04_8*burn))
gbu = 0.038_8*burn**0.28_8;
h = 1.0_8/(1.0_8+396.0_8*exp(-6380.0_8/t));
l=3.5E+9_8/(t*t)*exp(-16361_8/t);
CondFuel_FRAPCON = (1.0_8/(a+b*t+fbu+f*gbu*h)+l);
END FUNCTION

FUNCTION CondFuel(t)
!Fuel thermal conductivity Clad
REAL :: t, CondFuel
REAL :: akfuel(0:5)
akfuel = ThOpt%kFUelCorrelation
CondFuel = akfuel(0) + t*(akfuel(1)+t*(akfuel(2)+t*akfuel(3)))+akfuel(4)/(t-akfuel(5))
END FUNCTION

FUNCTION CondClad(t)
!Fuel thermal conductivity Clad
REAL :: t, CondClad
REAL :: akclad(0:3)
akclad = ThOpt%kCladCorrelation
CondClad = akclad(0) + t*(akclad(1) + t*(akclad(2) + t*akclad(3)))
END FUNCTION

FUNCTION GapCon(burnup)
USE UtilFunction,  ONLY : LineIntPol
IMPLICIT NONE
REAL :: GapCon
REAL :: burnup
REAL :: Xdat(19), Ydat(19)
INTEGER :: ndat
DATA XDat / 0.00010_8,  0.04999_8,  0.50000_8,  0.99998_8,  1.99997_8,  2.99997_8,  3.99997_8,  4.99994_8,         &
            7.49995_8,  9.99993_8, 14.99966_8, 19.99966_8, 24.99994_8, 30.00041_8, 35.00067_8, 40.00068_8,   &
           45.00016_8, 49.99991_8, 55.00042_8 /
DATA YDat / 0.77617_8,  0.81695_8,  0.85461_8,  0.90683_8,  0.97879_8,  1.04153_8,  1.11385_8,  1.19990_8,         &
            1.52150_8,  2.21212_8,  3.77267_8,  4.36162_8,  5.11622_8,  5.92855_8,  6.46608_8,  7.07055_8,         &
            7.76742_8,  8.45296_8,  8.75555_8/
ndat = 19
GapCon = LineIntPol(burnup, ndat, Xdat, Ydat)
GapCon = GapCon * 10000
END FUNCTION

FUNCTION GapCon_ROPER(pw, burnup)
USE UtilFunction,  ONLY : LineIntPol
IMPLICIT NONE
REAL :: pw, burnup
REAL :: GapCon_ROPER
REAL(8) :: burnuplv(12)
REAL(8) :: plv(15)
REAL(8) :: hcond_dat(12,15)
REAL(8) :: wt1, wt2, dat1, dat2
INTEGER, PARAMETER :: nburnup = 12 , nplv = 15
INTEGER :: i, ix1, ix2, iy1, iy2

data burnuplv    /   0.00009_8,   0.10008_8,   0.15008_8,   0.50008_8,   1.00008_8,   2.00008_8,  10.00008_8,  20.00008_8,  30.00008_8,  40.00008_8,  50.00008_8,  60.00008_8  /
data plv         /  4921.25984_8,     9842.51968_8,    14763.77953_8,    19685.03937_8,    24606.29921_8,    29527.55906_8,    34448.81890_8,    39370.07874_8,    44291.33858_8,  &
                  49212.59843_8,    54133.85827_8,    59055.11811_8,    65498.68766_8,    69724.40945_8,    74793.30709_8 /
data hcond_dat /    0.5386_8,    0.6291_8,    0.6530_8,    0.7451_8,    0.8196_8,    0.9248_8,    1.8042_8,    3.8379_8,    2.6063_8,    2.4169_8,    2.0812_8,    1.9920_8, &
                    0.6110_8,    0.7195_8,    0.7486_8,    0.8624_8,    0.9565_8,    1.0935_8,    2.5508_8,    3.9339_8,    2.6773_8,    2.4851_8,    2.1430_8,    2.0535_8, &
                    0.6976_8,    0.8310_8,    0.8674_8,    1.0126_8,    1.1363_8,    1.3237_8,    4.2402_8,    4.0269_8,    2.7448_8,    2.5495_8,    2.2006_8,    2.2216_8, &
                    0.8056_8,    0.9753_8,    1.0226_8,    1.2160_8,    1.3874_8,    1.6598_8,    5.2377_8,    4.1160_8,    2.8090_8,    2.6108_8,    2.2645_8,    3.3180_8, &
                    0.9479_8,    1.1739_8,    1.2388_8,    1.5139_8,    1.7720_8,    2.2137_8,    5.3451_8,    4.2031_8,    2.8716_8,    4.0363_8,    3.4951_8,    3.5248_8, &
                    1.1476_8,    1.4697_8,    1.5665_8,    1.9995_8,    2.4445_8,    3.3129_8,    5.4505_8,    4.2881_8,    4.4323_8,    4.2791_8,    3.6379_8,    3.7515_8, &
                    1.4728_8,    2.0042_8,    2.1772_8,    3.0452_8,    4.1396_8,    6.1274_8,    5.4639_8,    4.4435_8,    4.6260_8,    4.5378_8,    3.8813_8,    3.9950_8, &
                    2.1527_8,    3.4081_8,    3.9127_8,    6.2005_8,    6.2004_8,    6.1970_8,    5.4167_8,    6.3518_8,    4.8907_8,    4.7950_8,    4.0696_8,    4.1171_8, &
                    3.3977_8,    6.2602_8,    6.2615_8,    6.2645_8,    6.2649_8,    6.2618_8,    5.4291_8,    6.5541_8,    5.1801_8,    5.0128_8,    4.1708_8,    4.1739_8, &
                    6.3138_8,    6.3216_8,    6.3232_8,    6.3272_8,    6.3277_8,    6.3247_8,    5.4701_8,    6.8986_8,    5.4074_8,    5.0821_8,    4.2240_8,    4.2237_8, &
                    6.3748_8,    6.3841_8,    6.3857_8,    6.3900_8,    6.3907_8,    6.3878_8,    7.8617_8,    7.1301_8,    5.4778_8,    5.1385_8,    4.2703_8,    4.2707_8, &
                    6.4371_8,    6.4466_8,    6.4482_8,    6.4526_8,    6.4533_8,    6.4506_8,    8.5742_8,    7.1031_8,    5.5348_8,    5.1903_8,    4.3139_8,    4.3163_8, &
                    6.5117_8,    8.3233_8,    9.1981_8,    9.8088_8,    9.9022_8,    9.9531_8,    9.1824_8,    7.1764_8,    5.6025_8,    5.2554_8,    4.3702_8,    4.3761_8, &
                    9.4854_8,   10.0087_8,   10.0805_8,   10.2836_8,   10.3719_8,   10.4141_8,    9.5537_8,    7.2258_8,    5.6450_8,    5.2975_8,    4.4014_8,    4.4124_8, &
                   10.2496_8,   10.4339_8,   10.4684_8,   10.6029_8,   10.7157_8,   10.8386_8,    9.6640_8,    7.2802_8,    5.6941_8,    5.3469_8,    4.4347_8,    4.7556_8  /


DO i = 2, nplv
  if(pw .le. plv(i)) EXIT
ENDDO
IF(i .GT. nplv) i = nplv
iy1 = i-1; iy2 = i

DO i = 2, nburnup
  if(burnup .le. burnuplv(i)) EXIT
ENDDO
IF(i .GT. nburnup) i = nburnup
ix1 = i-1; ix2 = i

wt2 = (pw - plv(iy1)) / (plv(iy2) - plv(iy1))
wt1 = 1._8 - wt2
dat1 = wt2 * hcond_dat(ix1, iy2) + wt1 * hcond_dat(ix1, iy1)
dat2 = wt2 * hcond_dat(ix2, iy2) + wt1 * hcond_dat(ix2, iy1)

wt2 = (burnup - burnuplv(ix1)) / (burnuplv(ix2) - burnuplv(ix1))
wt1 = 1._8 - wt2
GapCon_ROPER = wt2 * dat2 + wt1 * dat1
GapCon_ROPER = GapCon_ROPER * 10000._8
END FUNCTION

FUNCTION FuelCp(T)
REAL :: FuelCp, T
REAL :: arcpfuel(0:3)
arcpfuel = ThOpt%CpFuelCorrelation
FuelCp = arcpfuel(0) + t*(arcpfuel(1) + t*(arcpfuel(2)+t*arcpfuel(3)))
END FUNCTION

FUNCTION CladCp(T)
REAL :: CladCp, T
REAL :: arcpclad(0:3)
arcpclad = ThOpt%CpCladCorrelation
CladCp = arcpclad(0) + t*(arcpclad(1) + t*(arcpclad(2)+t*arcpclad(3)))
END FUNCTION

SUBROUTINE ftfavg(x, xvol, r, nrpellet)
!
!REAL :: ftfavg
REAL :: x(nrpellet + 4),xvol(nrpellet + 5), r(nrpellet + 4) !dmm
REAL :: xvol2(20)
INTEGER :: nrpellet

REAL :: xavg, area, areaR, areaL, delr
INTEGER :: nrp1, nrp2, nrp3, nrp4, nrm1
INTEGER :: i

nrp1 = nrpellet + 1; nrp2 = nrpellet + 2
nrp3 = nrpellet + 3; nrp4 = nrpellet + 4
delR = R(2) - R(1)
!
i=1
xavg=(7*x(i)+x(i+1))*0.25
xvol(1)=xavg
do i=2,nrpellet
   area=2/3.*i*x(i+1)+44/3.*(i-1)*x(i)+2/3.*(i-2)*x(i-1)
   xavg=xavg+area
   arear=(-7/8.+4/3.*i)*x(i+1)+(21/12.+22/3.*(i-1))*x(i)+(-7/8.-2/3.*(i-2))*x(i-1)
   areal=area-arear
   xvol(i-1)=xvol(i-1)+areal
   xvol(i)=arear
enddo
nrm1=nrpellet-1
area=(16/3.*nrm1+4+5/24.)*x(nrpellet+1)+(10/3.*nrm1+2.25_8)*x(nrpellet)-(2/3.*nrm1+11/24.)*x(nrpellet-1)
xvol(nrpellet + 5)=(xavg+area)/(8.*nrpellet*nrpellet)
xvol(nrpellet)=xvol(nrpellet)+area
!
xvol(1)=xvol(1)/8._8
do i=2,nrpellet
  xvol(i)=xvol(i)/((2.*r(i)/delr+1)*8.)
enddo
xvol(nrp1)=half*(x(nrp1)+x(nrp2))  !gap
xvol(nrp2)=half*(x(nrp2)+x(nrp4))  !clad
!xvol(nrp3)=half*(x(nrp3)+x(nrp4))
RETURN
END SUBROUTINE

FUNCTION fenthalf(tcel)
!
REAL ::  fenthalf, t, tcel
t=tcel+273.15
fenthalf=(162.3_8+t*(0.1519_8+t*(-7.970e-5_8+t*1.601e-8_8)))*t
return
end function
!!!!!!!!!!!
! FUNCTIONs for water properties at 15.5 MPa
! cubic polynomial for thermal conductivity, max err=0.0204%
!  15.5 Mpa,  280 < T < 340, k in w/m-C, T in C
FUNCTION fcond(t)
REAL ::  fcond, t
fcond=8.9182016e-01+t*(-2.1996892e-03+t*(9.9347652e-06+t*(-2.0862471e-08)))
return
end FUNCTION

! cubic polynomial for density as a FUNCTION of temperature, max err=0.0692%
!  15.5 Mpa,  280 < T < 340, rho in Kg/M^3, T in C
FUNCTION fdens(t)
!
REAL :: fdens, t
fdens=5.9901166e+03+t*(-5.1618182e+01+t*(1.7541848e-01+t*(-2.0613054e-04)))
return
end FUNCTION

! cubic polynomial for density as a FUNCTION of enthalpy, max err=0.0112%
!  15.5 Mpa,  280 < T(h) < 340, rho in Kg/M^3, input h in J/Kg
FUNCTION fdensh(h)
REAL ::  fdensh, h, hk
hk=epsm3*h
fdensh=1.4532039e+03+hk*(-1.1243975e+00+hk*(7.4502004e-04+hk*(-2.3216531e-07)))
return
end FUNCTION

! cubic polynomial for enthalpy as a FUNCTION of temperature, max err=0.0306%
!  15.5 Mpa,  280 < T < 340, output h in J/Kg, T in C
FUNCTION fenthal(t)
REAL ::  fenthal, t, y
y=-5.9301427e+03+t*(6.5488800e+01+t*(-2.1237562e-01+t*(2.4941725e-04)))
fenthal=y*1000
return
end FUNCTION

! quartic polynomial for heat cappacity, max error=0.3053%
!  15.5 Mpa,  280 < T < 340, output h in J/Kg-C, T in C
FUNCTION fhcap(t)
REAL :: fhcap, t, y
y=3.0455749e+03+t*(-4.0684599e+01+t*(2.0411250e-01+t*(-4.5526705e-04+t*(3.8115453e-07))))
fhcap=y*1000
return
end FUNCTION

! wall-to-coolant heat transfer coeffcient in w/m^2-C
FUNCTION fhtcoef(t,deq,rhou)
REAL :: fhtcoef, t, deq, rhou
REAL :: k, mu, cp, pr, re
k=fcond(t)
mu=fvisco(t)
cp=fhcap(t)
pr=cp*mu/k
re=deq*rhou/mu
fhtcoef=0.023_8*k/deq*(pr**0.4)*(re**0.8)
return
end FUNCTION

! cubic polynomial for temperature as a FUNCTION of enthalpy, max err=0.0055%
!  15.5 Mpa,  280 < T < 340, T in C, input h in J/Kg
FUNCTION ftemp(h)
REAL(8) ftemp, h, hk
hk=epsm3*h
ftemp=1.4851739e+02+hk*(-1.2764991e-01+hk*(3.0781294e-04+hk*(-9.5429959e-08)))
return
end FUNCTION
! cubic polynomial for viscosity, max err=0.0641%
!  15.5 Mpa,  280 < T < 340, mu in Pa-sec, T in C
FUNCTION fvisco(t)
REAL(8) fvisco, t
fvisco=9.0836878e-04+t*(-7.4542195e-06+t*(2.3658072e-08+t*(-2.6398601e-11)))
return
end FUNCTION
! find enthalpy for given temperature, input h is the initial guess
! use fenthal instead of this unless it is REALly necessary
FUNCTION findh(t,h)
REAL :: findh, t, h
REAL :: x1, x2, y1, y2, slope, xn
INTEGER :: i, j
x1=1.01_8*h
x2=h
y1=ftemp(x1)
y2=ftemp(x2)
DO i=1,20
   IF(ABS(x2-x1).lt.epsm2) go to 100
   slope=(y2-y1)/(x2-x1)
   xn=x2+(t-y2)/slope
   x1=x2
   y1=y2
   x2=xn
   y2=ftemp(x2)
ENDDO
100 findh=xn
IF(i.eq.20) PRINT '("Fail to find enthalpy for temperature",e12.5)',t
RETURN
END FUNCTION

END MODULE

MODULE Water155_mod
USE PARAM
IMPLICIT NONE

CONTAINS
! cubic polynomial for viscosity, max err=0.0641%
!  15.5 Mpa,  280 < T < 340, mu in Pa-sec, T in C
FUNCTION FVISCO(T)
REAL FVISCO, T
FVISCO=9.0836878E-04_8+T*(-7.4542195E-06_8+T*(2.3658072E-08_8+T*(-2.6398601E-11_8)))
RETURN
END FUNCTION

! CUBIC POLYNOMIAL FOR THERMAL CONDUCTIVITY, MAX ERR=0.0204%
! 15.5 MPA,  280 < T < 340, K IN W/M-C, T IN C
FUNCTION FCOND(T)
REAL T, FCOND
FCOND=8.9182016E-01_8+T*(-2.1996892E-03_8+T*(9.9347652E-06_8+T*(-2.0862471E-08_8)))
END FUNCTION

! CUBIC POLYNOMIAL FOR DENSITY AS A FUNCTION OF TEMPERATURE, MAX ERR=0.0692%
!  15.5 MPA,  280 < T < 340, RHO IN KG/M^3, T IN C
FUNCTION FDENS(T)
REAL FDENS,T
FDENS=5.9901166E+03_8+T*(-5.1618182E+01_8+T*(1.7541848E-01_8+T*(-2.0613054E-04_8)))
!      FDENS=FDENS*0.9987546 !FOR MCCARD
RETURN
END FUNCTION

! CUBIC POLYNOMIAL FOR DENSITY AS A FUNCTION OF ENTHALPY, MAX ERR=0.0112%
!  15.5 MPA,  280 < T(H) < 340, RHO IN KG/M^3, INPUT H IN J/KG
FUNCTION FDENSH(H)
REAL FDENSH, H
REAL :: HK
HK=EPSM3*H
FDENSH=1.4532039E+03_8+HK*(-1.1243975E+00_8+HK*(7.4502004E-04_8+HK*(-2.3216531E-07_8)))
RETURN
END FUNCTION

! CUBIC POLYNOMIAL FOR ENTHALPY AS A FUNCTION OF TEMPERATURE, MAX ERR=0.0306%
!  15.5 MPA,  280 < T < 340, OUTPUT H IN J/KG, T IN C
FUNCTION FENTHAL(T)

REAL FENTHAL, T
REAL Y
Y=-5.9301427E+03_8+T*(6.5488800E+01_8+T*(-2.1237562E-01_8+T*(2.4941725E-04_8)))
FENTHAL=Y*1000_8
RETURN
END FUNCTION

! QUARTIC POLYNOMIAL FOR HEAT CAPPACITY, MAX ERROR=0.3053%
!  15.5 MPA,  280 < T < 340, OUTPUT H IN J/KG-C, T IN C
FUNCTION FHCAP(T)
REAL FHCAP, T
REAL Y
Y=3.0455749E+03_8+T*(-4.0684599E+01_8+T*(2.0411250E-01_8+T*(-4.5526705E-04_8+T*(3.8115453E-07_8))))
FHCAP=Y*1000_8
RETURN
END FUNCTION

! WALL-TO-COOLANT HEAT TRANSFER COEFFCIENT IN W/M^2-C
FUNCTION FHTCOEF(T,DEQ,RHOU)
REAL FHTCOEF, T, DEQ, RHOU
REAL K,MU, CP, PR, RE
K=FCOND(T)
MU=FVISCO(T)
CP=FHCAP(T)
PR=CP*MU/K
RE=DEQ*RHOU/MU
FHTCOEF=0.023_8*K/DEQ*(PR**0.4_8)*(RE**0.8_8)
RETURN
END FUNCTION

! CUBIC POLYNOMIAL FOR TEMPERATURE AS A FUNCTION OF ENTHALPY, MAX ERR=0.0055%
!  15.5 MPA,  280 < T < 340, T IN C, INPUT H IN J/KG
FUNCTION FTEMP(H)
REAL FTEMP, H
REAL HK
HK=EPSM3*H
FTEMP=1.4851739E+02_8+HK*(-1.2764991E-01_8+HK*(3.0781294E-04_8+HK*(-9.5429959E-08_8)))
RETURN
END FUNCTION

! FIND ENTHALPY FOR GIVEN TEMPERATURE, INPUT H IS THE INITIAL GUESS
! USE FENTHAL INSTEAD OF THIS UNLESS IT IS REALLY NECESSARY
FUNCTION FINDH(T,H)
REAL FINDH, T, H
REAL X1, X2, XN, Y1, Y2, SLOPE
INTEGER :: I
X1=1.01_8*H
X2=H
Y1=FTEMP(X1)
Y2=FTEMP(X2)
DO I=1,20
   IF(ABS(X2-X1).LT.EPSM2) GO TO 100
   SLOPE=(Y2-Y1)/(X2-X1)
   XN=X2+(T-Y2)/SLOPE
   X1=X2
   Y1=Y2
   X2=XN
   Y2=FTEMP(X2)
ENDDO
100 FINDH=XN
IF(I.EQ.20) PRINT '("FAIL TO FIND ENTHALPY FOR TEMPERATURE",E12.5)',T
RETURN
END FUNCTION

END MODULE