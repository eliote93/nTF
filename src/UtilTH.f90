#include <defines.h>
SUBROUTINE SetGTnGapCoolTemp(Core, ThInfo, PE)
USE PARAM
USE TYPEDEF,       ONLY : CoreInfo_Type,         ThInfo_Type,       PE_TYPE,      &
                          Pin_Type,              Asy_Type,          AsyInfo_Type
USE TH_Mod,        ONLY : ThOpt,                 ThVar
USE SteamTBL_mod,  ONLY : steamtbl

IMPLICIT NONE
TYPE (CoreInfo_Type) :: Core
TYPE (ThInfo_Type) :: ThInfo
TYPE (PE_Type) :: PE

TYPE (Pin_Type), POINTER :: Pin(:)
TYPE (Asy_Type), POINTER :: Asy(:)
TYPE (AsyInfo_Type), POINTER :: AsyInfo(:)

REAL, POINTER :: TCool(:, :), DenCool(:, :)
REAL :: TSum(1000), nFuelCell(1000), TsumOut(1000), TsumIn(1000)
REAL :: wt, wh, hin, wrho, wvin, wxin, wbetain, wkapain, wcpin
REAL :: PEXIT

INTEGER :: nxy, nz, nxya, nxyc
INTEGER :: iz, ixy, ixya, itype
INTEGER :: i, j

Pin => Core%Pin
Asy => Core%Asy;  AsyInfo => Core%Asyinfo
TCool => ThInfo%TCool; DenCool => THInfo%DenCool

nxy = Core%nxy; nz = ThVar%nzth
PEXIT = ThInfo%PEXIT
nxya = Core%nxya
DO iz = 1, nz
  Tsum(1:nxya) = 0; nFuelCell(1:nxya) = 0
  TsumIn(1:nxya) = 0;TsumOut(1:nxya) = 0;
  DO ixy = 1, nxy
    ixya = Pin(ixy)%iasy; itype = Asy(ixya)%AsyType
    IF (.NOT. Pin(ixy)%lFuel) CYCLE
    IF (.NOT. AsyInfo(iType)%lFuel) CYCLE
    nFuelCell(ixya) = nFuelCell(ixya) + 1
    TSum(ixya) = Tsum(ixya) + Tcool(iz, ixy)
    TsumOut(ixya) = TsumOut(ixya) + ThInfo%TCoolInOut(2, iz, ixy)
    TsumIn(ixya) = TsumIn(ixya) + ThInfo%TCoolInOut(1, iz, ixy)
  END DO
  DO ixya = 1, nxya
    IF (Tsum(ixya) .GT. epsm6) Tsum(ixya) = Tsum(ixya) / nFuelCell(ixya)
    IF (Tsum(ixya) .GT. epsm6) TsumIn(ixya) = TsumIn(ixya) / nFuelCell(ixya)
    IF (Tsum(ixya) .GT. epsm6) TsumOut(ixya) = TsumOut(ixya) / nFuelCell(ixya)
  END DO
  DO ixy = 1, nxy
    ixya = Pin(ixy)%iasy; itype = Asy(ixya)%AsyType
    !IF (Pin(ixy)%lFuel .AND. .NOT. Pin(ixy)%lGd) CYCLE
    IF (Pin(ixy)%lFuel) CYCLE
    IF (.NOT. AsyInfo(iType)%lFuel) CYCLE
    TCool(iz, ixy) = Tsum(ixya)
    wt = Tsum(ixya) + CKELVIN
    CALL steamtbl(TRUE,pEXIT,wt,wh,wrho,wvin,wxin,wbetain,wkapain,wcpin)
    DenCool(iz, ixy) = wrho
    ThInfo%TCoolInOut(2, iz, ixya) = TSumOut(ixya)
    ThInfo%TCoolInOut(1, iz, ixya) = TSumIn(ixya)
  END DO
END DO

NULLIFY(Pin, Asy, AsyInfo)
END SUBROUTINE
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE BurnupUpdate(is_reset)
  USE DEPL_MOD,     ONLY: DeplCntl
  implicit NONE
  LOGICAL, intent(inout) :: is_reset
  INTEGER, save :: bu_step_old = -1
  INTEGER :: bu_step_new, PC_OPT
  LOGICAL :: LCORRECT
  LOGICAL, save :: LCORRECT_SAVE = .false.

  bu_step_new = DeplCntl%NowStep
  PC_OPT = DeplCntl%PC_OPT

  IF (bu_step_new > bu_step_old) THEN
    IF (PC_OPT == 1) THEN
      is_reset = .true.
    else
      is_reset = .true.
      LCORRECT_SAVE = .false.
    end if
  else
    IF (PC_OPT == 1) THEN
      is_reset = .false.
    else
      LCORRECT = .not. DeplCntl%LPREDICT
      IF (LCORRECT_SAVE .NEQV. LCORRECT) THEN
        is_reset = .true.
      else
        is_reset = .false.
      END IF
      LCORRECT_SAVE = LCORRECT
    END IF
  END IF

  bu_step_old = bu_step_new

END SUBROUTINE
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SteadyCoolantTH(powlin, plevel, PEXIT, Tout, RelPW, CoolantTH, ThVar)

USE PARAM
USE TYPEDEF,      ONLY : CoolantTH_Type, ThVar_Type
USE SteamTBL_mod, ONLY : steamtbl

IMPLICIT NONE

TYPE (CoolantTH_Type) :: CoolantTH ! One Channel Coolant TH in
TYPE (ThVar_Type) :: ThVar
REAL :: PowLin, Plevel, PEXIT, Tout
REAL :: RelPW(:)

REAL, POINTER, DIMENSION(:) :: hcool, rhou, rhohu, qvol, qeff, Tcool, DenCool, hz
REAL, POINTER, DIMENSION(:,:) :: TcoolInOut
REAL :: qprime, qf, qc, qflux, qeffnew, acf, afp, xi, zeta, zetap, FracDC, Fracdf, rhouin, rhouout, rhohuin, rhohuout, hout, wt, wh, hin, wrho, wvin, wxin, wbetain, wkapain, wcpin
INTEGER :: nzth, iz
! ----------------------------------------------------

acf    = ThVar%acf
afp    = ThVar%afp
xi     = ThVar%Xi
zeta   = ThVar%zeta
zetap  = ThVar%zetap
FracDc = ThVar%FracDC
FracDf = ThVar%FracDf
nzth   = ThVar%nzth
hz    => ThVar%hz

hcool      => CoolantTH%hcool
rhou       => CoolantTH%rhou
rhohu      => CoolantTH%rhohu
qvol       => CoolantTH%qvol
Tcool      => CoolantTH%Tcool
qeff       => CoolantTH%qeff
TCoolInOut => CoolantTh%TcoolInOut
DenCool    => CoolantTH%DenCool

rhouin  = rhou (0)
rhohuin = rhohu(0)

! INIT : Temperature Calculation
wh = rhohuin / rhouin
CALL steamtbl(FALSE,pEXIT,wt,wh,wrho,wvin,wxin,wbetain,wkapain,wcpin)
TCoolInOut(1, 1) = wt - CKELVIN

DO iz = 1, nzth
  qprime = plevel * powlin * RelPW(iz)
  qf     = FracDf * qprime / afp;  qc = FracDc * qprime / acf
  
  qvol(iz) = qf ! Heat source of Fuel Pallet
  qflux    = qf * afp / zeta
  qeffnew  = qflux * zetap + qc
  qeff(iz) = qeffnew
  rhohuout = rhohuin + qeffnew * hz(iz)
  hcool(iz) = 0.5 * (rhohuout + rhohuin) / rhouin
  
  wh = hcool(iz)
  CALL steamtbl(FALSE,pEXIT,wt,wh,wrho,wvin,wxin,wbetain,wkapain,wcpin)
  
  Tcool  (iz) = wt - CKELVIN
  Dencool(iz) = wrho
  rhohuin     = rhohuout
  rhohu(iz)   = rhohuout
  
  ! Out Temperature
  wh = rhohuout / rhouin
  CALL steamtbl(FALSE,pEXIT,wt,wh,wrho,wvin,wxin,wbetain,wkapain,wcpin)
  
  TCoolInOut(2, iz)   = wt - CKELVIN
  TCoolInOut(1, iz+1) = TCoolInOut(2, iz)
END DO

hout = RhoHu(nzth)/rhou(nzth)

wh = hout
CALL steamtbl(FALSE,pEXIT,wt,wh,wrho,wvin,wxin,wbetain,wkapain,wcpin)

tout = wt - CKELVIN

NULLIFY (hcool)
NULLIFY (rhou)
NULLIFY (rhohu)
NULLIFY (qvol)
NULLIFY (Tcool)
NULLIFY (qeff)
! ----------------------------------------------------

END SUBROUTINE SteadyCoolantTH
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SteadyCoolantTH_ThCh(Core, powlin, plevel, PEXIT, Tout, RelPW, CoolantTH, ThVar,ThOpt, PE, ixy)
USE PARAM
USE TYPEDEF,      ONLY : CoreInfo_Type,        CoolantTH_Type,       ThVar_Type,     &
                         ThOpt_Type,           PE_Type
USE SteamTBL_mod, ONLY : steamtbl
IMPLICIT NONE
TYPE (CoreInfo_Type) :: Core
TYPE (CoolantTH_Type) :: CoolantTH        ! One Channel Coolant TH in
TYPE (ThVar_Type) :: ThVar
TYPE (ThOPT_Type) :: ThOpt                !
TYPE (PE_Type) :: PE                      !
REAL :: PowLin, Plevel, PEXIT, Tout           !
REAL :: RelPW(:)
INTEGER :: ixy

REAL, POINTER :: hcool(:), rhou(:), rhohu(:), u(:), ud(:), qvol(:), qeff(:), Tcool(:), DenCool(:)
REAL, POINTER :: TcoolInOut(:, :)
REAL, POINTER :: hz(:)
REAL :: qprime, qf, qc, qflux, qeffnew              !
REAL :: acf, afp, xi, zeta, zetap, FracDC, Fracdf
REAL :: rhouin, rhouout, rhohuin, rhohuout, hout
REAL :: Btemp
REAL :: wt, wh, hin, wrho, wvin, wxin, wbetain, wkapain, wcpin
INTEGER :: nzth
INTEGER :: iz, ixya, ichtyp

ixya = Core%Pin(ixy)%iasy
ichtyp = Core%ThChMap(ixya)

IF (ichtyp .EQ. 0) THEN
  acf = ThVar%acf
  xi = ThVar%xi
  zetap = ThVar%zetap
ELSE
  acf = ThVar%ThCh(ichtyp)%acf
  xi = ThVar%ThCh(ichtyp)%xi
  zetap = ThVar%ThCh(ichtyp)%zetap
END IF

afp = ThVar%afp
zeta = ThVar%zeta
FracDc = ThVar%FracDC; FracDf = ThVar%FracDf;
nzth = ThVar%nzth
Btemp = ThVar%BoilingTemp

hcool => CoolantTH%hcool; rhou => CoolantTH%rhou
rhohu => CoolantTH%rhohu; u => CoolantTH%u
ud => CoolantTH%ud; qvol => CoolantTH%qvol
Tcool => CoolantTH%Tcool; qeff => CoolantTH%qeff
TCoolInOut => CoolantTh%TcoolInOut
DenCool => CoolantTH%DenCool
hz => ThVar%hz
rhouin = rhou(0); rhohuin = rhohu(0)
!Init Temperature Calculation
wh = rhohuin / rhouin
CALL steamtbl(FALSE,pEXIT,wt,wh,wrho,wvin,wxin,wbetain,wkapain,wcpin)
TCoolInOut(1, 1) = wt -CKELVIN

DO iz = 1, nzth
  qprime = plevel * powlin * RelPW(iz)
  qf = FracDf * qprime / afp;  qc = FracDc * qprime / acf
  qvol(iz) = qf                         !Heat source of Fuel Pallet
  qflux = qf * afp / zeta
  qeffnew = qflux * zetap + qc
  qeff(iz) = qeffnew
  rhohuout = rhohuin + qeffnew * hz(iz)
  hcool(iz) = 0.5 * (rhohuout+rhohuin) / rhouin
  wh = hcool(iz)
  CALL steamtbl(FALSE,pEXIT,wt,wh,wrho,wvin,wxin,wbetain,wkapain,wcpin)
  Tcool(iz) = wt - CKELVIN
  Dencool(iz) = wrho
  RhoHUin = RhoHuOut
  RhoHU(iz) = RhoHuout
  !Out Temperature
  wh = rhohuout / rhouin
  CALL steamtbl(FALSE,pEXIT,wt,wh,wrho,wvin,wxin,wbetain,wkapain,wcpin)
  TCoolInOut(2, iz) = wt - CKELVIN
  TCoolInOut(1, iz+1) = TCoolInOut(2, iz)
END DO
hout = RhoHu(nzth)/rhou(nzth); wh = hout
CALL steamtbl(FALSE,pEXIT,wt,wh,wrho,wvin,wxin,wbetain,wkapain,wcpin)
tout = wt - CKELVIN

NULLIFY(hcool); NULLIFY(rhou)
NULLIFY(rhohu); NULLIFY(u)
NULLIFY(ud); NULLIFY(qvol)
NULLIFY(Tcool); NULLIFY(qeff)

  END SUBROUTINE
! ------------------------------------------------------------------------------------------------------------
!SUBROUTINE SteadyFuelConduction(powlin, plevel, Tfmax, RelPw, PwShape, BurnUp, RhoU, FuelTH, ThVar, ThOpt, nTracerCntl, PE, hGapArray)
SUBROUTINE SteadyFuelConduction(powlin, plevel, Tfmax, RelPw, PwShape, BurnUp, RhoU, FuelTH, ThVar, ThOpt, nTracerCntl)

USE PARAM
USE TYPEDEF,          ONLY : FuelTH_Type, ThVar_Type, ThOpt_Type
USE CNTL,             ONLY : nTracerCntl_Type
USE FuelProperty_Mod, ONLY : fhtcoef, INTGAPCOND
USE TH_Mod,           ONLY : RadFuelCondFDM

IMPLICIT NONE

TYPE (FuelTH_Type):: FuelTH  ! One Channel Coolant TH information
TYPE (ThVar_Type) :: ThVar
TYPE (ThOPT_Type) :: ThOpt
TYPE (nTracerCntl_Type) :: nTracerCntl

REAL :: PowLin,  Plevel, Tfmax
REAL :: RelPW(:)
REAL, POINTER :: RhoU(:), PwShape(:, :), BurnUp(:, :)
!REAL :: hGapArray(:)
! ----------------------------------------------------
REAL, POINTER, DIMENSION(:)   :: qvol, tcool, htcoef
REAL, POINTER, DIMENSION(:,:) :: tfuel, tfvol
REAL :: DEQ, QF, fracdf, qprime, afp
INTEGER :: iz, nzth, nr
! ----------------------------------------------------

qvol   => FuelTh%qvol
tcool  => FuelTH%tcool
htcoef => FuelTh%htcoef
tfuel  => FuelTh%tfuel
tfvol  => FuelTH%tfvol

fracdf = THVar%fracdf
nzth   = ThVar%nzth
DEQ    = ThVar%Deq
afp    = ThVar%afp
nr     = Thvar%npr5

Tfmax = ZERO

DO iz = 1, nzth
  qprime = plevel * powlin * RelPw(iz)
  !IF (ThOpt%hGapModel .EQ. 2) THEN
  !  ThOpt%hgap = INTGAPCOND(qprime)
  !END IF
  qf = fracdf * qprime / afp
  qvol(iz) = qf
  htcoef(iz) = fhtcoef(tcool(iz), deq, RhoU(iz))
  
  IF (.NOT. nTracerCntl%LIHP) Pwshape = 1.
  
  IF (nTracerCntl%lFuelCondFDM) THEN
    CALL RadFuelCondFDM(tfvol(1:nr-3, iz), tfuel(1:nr, iz), tcool(iz), htcoef(iz), qf, PwShape(:, iz), nr, ThVar, ThOpt, FuelTH%lMox(iz))
  ELSE
    CALL tfcalss(tfvol(1:nr-3, iz), tfuel(1:nr, iz), tcool(iz), htcoef(iz), qf, nr, ThVar, ThOpt, FuelTH%lMox(iz))
    !CALL tfcalss_pwshape(tfvol(1:nr-3, iz), tfuel(1:nr, iz), tcool(iz), htcoef(iz), qf, PwShape(:,iz), nr, ThVar, ThOpt, FuelTH%lMox(iz), hGapArray(iz))
    !CALL tfcalss_qshape(tfvol(1:nr-3, i0z), tfuel(1:nr, iz), tcool(iz), htcoef(iz), qf, qprime, PwShape(:,iz), Burnup(:, iz), nr, ThVar, ThOpt, FuelTH%lMox(iz))
  END IF
  
  IF (nTracerCntl%lEffTemp) CALL SetEffectiveTemp(tfvol(1:nr-3, iz),  tfuel(1:nr, iz), BurnUp(0, iz), nr, ThVar, ThOpt)
  CALL CalEffectiveTemp(FuelTh%Teff(iz), tfvol(1:nr-3, iz), tfuel(1:nr, iz), BurnUp(0, iz), nr, ThVar, ThOpt)
  
  tfmax = max(tfmax, tfuel(1, iz))
END DO

NULLIFY (qvol, tcool, htcoef, tfuel, tfvol)
! ----------------------------------------------------

END SUBROUTINE SteadyFuelConduction
! ------------------------------------------------------------------------------------------------------------
!SUBROUTINE SteadyFuelConduction_ThCh(Core, powlin, plevel, Tfmax, RelPw, PwShape, BurnUp, RhoU,FuelTH, ThVar, ThOpt, nTracerCntl, PE, hGapArray, ixy)
!USE PARAM
!USE TYPEDEF,          ONLY :  CoreInfo_Type,    FuelTH_Type,     ThVar_Type,     ThOpt_Type,          &
!                              PE_TYPE
!USE CNTL,             ONLY : nTracerCntl_Type
!USE FuelProperty_Mod, ONLY : fhtcoef,    INTGAPCOND
!USE TH_Mod,           ONLY : RadFuelCondFDM
!IMPLICIT NONE
!TYPE (CoreInfo_Type) :: Core
!TYPE (FuelTH_Type):: FuelTH        ! One Channel Coolant TH information
!TYPE (ThVar_Type) :: ThVar                !
!TYPE (ThOPT_Type) :: ThOpt                !
!TYPE (nTracerCntl_Type) :: nTracerCntl
!TYPE (PE_Type) :: PE                      !
!REAL :: PowLin,  Plevel, Tfmax
!REAL :: RelPW(:)
!REAL, POINTER :: RhoU(:), PwShape(:, :), BurnUp(:, :)
!REAL :: hGapArray(:)
!INTEGER :: ixy
!
!REAL, POINTER :: hflux(:), qvol(:), tcool(:), htcoef(:), tfuel(:, :), tfvol(:, :)
!REAL :: DEQ, QF, fracdf, qprime, afp
!INTEGER :: iz, nzth, nr, ixya, ichtyp
!
!ixya = Core%Pin(ixy)%iasy
!ichtyp = Core%ThChMap(ixya)
!
!IF (ichtyp .EQ. 0) THEN
!  Deq = ThVar%Deq
!ELSE
!  Deq = ThVar%ThCh(ichtyp)%Deq
!END IF
!
!hflux => FuelTH%hflux; qvol => FuelTh%qvol
!tcool => FuelTH%tcool; htcoef => FuelTh%htcoef
!tfuel => FuelTh%tfuel; tfvol => FuelTH%tfvol
!fracdf = THVar%fracdf
!nzth = ThVar%nzth; DEQ =ThVar%Deq; afp = ThVar%afp
!nr = Thvar%npr5
!Tfmax = 0
!DO iz = 1, nzth
!  qprime = plevel * powlin * RelPw(iz)
!  IF (ThOpt%hGapModel .EQ. 2) THEN
!    ThOpt%hgap = INTGAPCOND(qprime)
!  END IF
!  qf = fracdf * qprime / afp; qvol(iz) = qf
!  htcoef(iz) = fhtcoef(tcool(iz), deq, RhoU(iz))
!  IF (.NOT. nTracerCntl%LIHP) Pwshape=1
!  IF (nTracerCntl%lFuelCondFDM) THEN
!    CALL RadFuelCondFDM(tfvol(1:nr-3, iz), tfuel(1:nr, iz), tcool(iz), htcoef(iz), qf, PwShape(:, iz), nr, ThVar, ThOpt, FuelTH%lMox(iz))
!  ELSE
!    !CALL tfcalss(tfvol(1:nr-3, iz), tfuel(1:nr, iz), tcool(iz), htcoef(iz), qf, nr, ThVar, ThOpt, FuelTH%lMox(iz), hGapArray(iz))
!    CALL tfcalss_pwshape(tfvol(1:nr-3, iz), tfuel(1:nr, iz), tcool(iz), htcoef(iz), qf, PwShape(:,iz), nr, ThVar, ThOpt, FuelTH%lMox(iz), hGapArray(iz))
!    !CALL tfcalss_qshape(tfvol(1:nr-3, iz), tfuel(1:nr, iz), tcool(iz), htcoef(iz), qf, qprime, PwShape(:,iz), Burnup(:, iz), nr, ThVar, ThOpt, FuelTH%lMox(iz))
!  END IF
!  IF (nTracerCntl%lEffTemp) CALL SetEffectiveTemp(tfvol(1:nr-3, iz),  tfuel(1:nr, iz), BurnUp(0, iz), nr, ThVar, ThOpt)
!  CALL CalEffectiveTemp(FuelTh%Teff(iz), tfvol(1:nr-3, iz), tfuel(1:nr, iz), BurnUp(0, iz), nr, ThVar, ThOpt)
!
!  Tfmax = max(tfmax, tfuel(1, iz))
!END DO
!NULLIFY(hflux, qvol, tcool, htcoef, tfuel, tfvol)
!END SUBROUTINE
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE recalculate_vars(FuelTH, BurnUp, Tfmax, ThVar, ThOpt, nTracerCntl, PE)
  USE PARAM
  USE TYPEDEF,          ONLY :  FuelTH_Type,     ThVar_Type,     ThOpt_Type,          &
                                PE_TYPE
  USE CNTL,             ONLY : nTracerCntl_Type
  IMPLICIT NONE
  TYPE (FuelTH_Type):: FuelTH        ! One Channel Coolant TH information
  TYPE (ThVar_Type) :: ThVar                !
  TYPE (ThOPT_Type) :: ThOpt                !
  TYPE (nTracerCntl_Type) :: nTracerCntl
  TYPE (PE_Type) :: PE
  REAL, intent(inout) :: Tfmax
  INTEGER :: nr, nzth, iz
  REAL, POINTER :: BurnUp(:, :)
  REAL, POINTER :: tfuel(:, :), tfvol(:, :)
  tfuel => FuelTh%tfuel; tfvol => FuelTH%tfvol
  nzth = ThVar%nzth
  nr = Thvar%npr5
  Tfmax = 0._8

  DO iz = 1, nzth
    call Adjust_after_AA(tfvol(1:nr-3, iz), tfuel(1:nr, iz), nr, ThVar, ThOpt)
    IF (nTracerCntl%lEffTemp) CALL SetEffectiveTemp(tfvol(1:nr-3, iz),  tfuel(1:nr, iz), BurnUp(0, iz), nr, ThVar, ThOpt)
    FuelTh%Teff(iz)=0._8
    CALL CalEffectiveTemp(FuelTh%Teff(iz), tfvol(1:nr-3, iz), tfuel(1:nr, iz), BurnUp(0, iz), nr, ThVar, ThOpt)
    Tfmax = max(tfmax, tfuel(1, iz))
  END DO

  NULLIFY(tfuel, tfvol)
  END SUBROUTINE
! ------------------------------------------------------------------------------------------------------------
!SUBROUTINE tfcalss(tfvol, Tfuel, tcool, htcoef, qf, nr, ThVar, ThOpt, lMox, hgap)
!USE PARAM
!USE TYPEDEF,           ONLY : ThVar_Type,     ThOpt_Type
!USE SteamTBL_mod,      ONLY : steamtbl
!USE FuelProperty_Mod,  ONLY : FTHCON,       CTHCON,       ftfavg,  RSGAPCOND
!USE LU_Mod,            ONLY : DirLu1Dsolver
!IMPLICIT NONE
!REAL :: tfvol(nr-3)
!REAL :: Tfuel(nr), tcool, htcoef, qf
!INTEGER :: nr
!TYPE (ThVar_Type) :: ThVar
!TYPE (ThOPT_Type) :: ThOPT
!LOGICAL :: lMox
!REAL, OPTIONAL :: hgap
!
!INTEGER :: hGapModel
!REAL :: diag(100), L(100), U(100), x(100), xd(100), b(100), xvol(100)
!REAL :: kf(100), kfb(100), kfm(100), kfmb(100)
!REAL :: kmr, kml, kgap, kgap2, kgap4, kconv, kconv1, tworm
!REAL, POINTER :: R(:)
!REAL :: delR, delR2, ri,  q, alpha
!REAL :: rs, rw, tw, rgap
!REAL :: err
!INTEGER :: npr, npr1, npr2, npr3, npr4, npr5
!INTEGER :: i, m
!
!hgapmodel = ThOpt%hGapModel
!
!npr = ThVar%npr; npr1 =ThVar%npr1
!npr2 = ThVar%npr2; npr3 = Thvar%npr3
!npr4= ThVar%npr4; npr5 = ThVar%npr5
!rs = ThVar%rs; rw = ThVar%rw
!tw = ThVar%tw; rgap = ThVar%rgap;
!R => ThVar%R
!DelR = R(2) - R(1); DelR2 = DelR * DelR
!hgap = ThOpt%hgap
!kgap = ThOpt%hgap * delr;
!kgap2 = ThOpt%hgap * tw * rs / rgap
!kgap4 = ThOpt%hgap * tw * (4._8 - tw / rgap) * rs / rgap
!tworm = tw / (rgap + 0.5_8 * tw)
!q = qf * DelR2
!
!DO
!  DO i = 1, npr4
!    x(i) = tfuel(i) + CKELVIN
!  END DO
!
!  DO i = 1, npr1
!    kf(i) = FTHCON(x(i), 0.)
!  END DO
!
!  IF (lMOX) kf(1:npr1) = 0.9 * kf(1:npr1)
!
!  DO i = 1, npr
!    kfm(i) = 0.5_8 * (kf(i) + kf(i+1))
!  END DO
!
!  DO i = npr2, npr4
!    kf(i) = CTHCON(x(i))
!  END DO
!  kmr = 0.5_8 * (kf(npr3) + kf(npr4))
!  kml = 0.5_8 * (kf(npr3) + kf(npr2))
!
!  IF (PRESENT(hGap)) THEN
!    IF (hGapModel .EQ. 3) THEN
!      hgap = RSGAPCOND(kf(npr1), kf(npr2), x(npr1), x(npr2))
!      kgap = hgap * delr;
!      kgap2 = hgap * tw * rs / rgap
!      kgap4 = hgap * tw * (4._8 - tw / rgap) * rs / rgap
!    END IF
!  END IF
!
!  x(1:npr4) = tfuel(1:npr4)
!  xd(1:npr4) = x(1:npr4)
!  !Matrix Set Up
!  m = 1
!  Diag(m) = 4._8 * kf(m)
!  U(m) = -4*kf(m); !L(m) = 0;
!  b(m) = q
!  i = m
!  DO m = 2, npr
!    ri = 1/dble(i)
!    Diag(m) = kfm(i) + kfm(m) + 0.5_8 * (kfm(m) - kfm(i)) * ri
!    L(m) = - kfm(i) * (1._8 - 0.5_8 * ri); U(m) = - kfm(m) * (1._8 + 0.5_8 * ri)
!    b(m) = q
!    i = m
!  END DO
!  m = npr1
!  alpha = kgap * (1._8 -kf(m-1)/kf(m))
!  Diag(m) = 2.*(kf(m)+kgap*(one+half/npr))+alpha
!  L(m) = -2.*kf(m); U(m) = -2.*kgap*(one+half/npr)-alpha
!  b(m)=q
!
!  m = npr2
!  alpha = 2.*kgap2*(kf(m+1)/kf(m)-1)
!  Diag(m) = 8.*kf(m)+kgap4-alpha
!  L(m) = -kgap4+alpha
!  U(m) =-8.*kf(m)
!  b(m) = 0.
!
!  m=npr3
!  Diag(m) = 4.*(kmr+kml)+tworm*(kmr-kml)
!  L(m) = -kml*(4-tworm)
!  U(m) = -kmr*(4+tworm)
!  B(m) = zero
!  m=npr4
!  kconv1 = htcoef * ThVar%tw
!  alpha = 2.*kconv1*(one-kf(m-1)/kf(m))
!  kconv = htcoef * tw * (4.+tw/rw)
!  Diag(m) = 8.*kf(m)+kconv+alpha
!  L(m) = -8.*kf(m); !U(m) = 0
!  b(m) = (kconv+alpha) * tcool
!
!  CALL DirLu1Dsolver(Diag(1:npr4), L(1:npr4), U(1:npr4), x(1:npr4), b(1:npr4), npr4)
!  !Obtained temperature is C unit not k(becaUSE the bulk temp is C unit)
!  !Check Error
!  err = 0
!  DO m = 1, npr4
!    tfuel(m) = x(m)
!    err = max(err, abs(x(m)-xd(m)) )
!  END DO
!  IF (err .lt. epsm2) EXIT
!END DO
!xvol = 0
!call ftfavg(x, xvol(1:npr5), r, npr)
!tfvol(1:npr2) = xvol(1:npr2) + CKELVIN
!tfuel(npr5) = xvol(npr5)
!NULLIFY(R)
!END SUBROUTINE

SUBROUTINE tfcalss(tfvol, Tfuel, tcool, htcoef, qf, nr, ThVar, ThOpt, lMox)
USE PARAM
USE TYPEDEF,           ONLY : ThVar_Type,     ThOpt_Type
USE SteamTBL_mod,      ONLY : steamtbl
USE FuelProperty_Mod,  ONLY : CondFuel,       CondClad,       ftfavg 
USE LU_Mod,            ONLY : DirLu1Dsolver
IMPLICIT NONE
REAL :: tfvol(nr-3)
REAL :: Tfuel(nr), tcool, htcoef, qf
INTEGER :: nr
TYPE(ThVar_Type) :: ThVar
TYPE(ThOPT_Type) :: ThOPT
LOGICAL :: lMox

REAL :: diag(100), L(100), U(100), x(100), xd(100), b(100), xvol(100)
REAL :: kf(100), kfb(100), kfm(100), kfmb(100)
REAL :: kmr, kml, kgap, kgap2, kgap4, kconv, kconv1, tworm
REAL, POINTER :: R(:)
REAL :: delR, delR2, ri,  q, alpha
REAL :: rs, rw, tw, rgap
REAL :: err
INTEGER :: npr, npr1, npr2, npr3, npr4, npr5
INTEGER :: i, m

npr = ThVar%npr; npr1 =ThVar%npr1
npr2 = ThVar%npr2; npr3 = Thvar%npr3
npr4= ThVar%npr4; npr5 = ThVar%npr5
rs = ThVar%rs; rw = ThVar%rw
tw = ThVar%tw; rgap = ThVar%rgap; 
R => ThVar%R
DelR = R(2) - R(1)
IF (ThVar%rsi /= 0) DelR = R(3) - R(2) 
DelR2 = DelR * DelR
kgap = ThOpt%hgap * delr; 
kgap2 = ThOpt%hgap * tw * rs / rgap
kgap4 = ThOpt%hgap * tw * (4._8 - tw / rgap) * rs / rgap
tworm = tw / (rgap + 0.5_8 * tw) 
q = qf * DelR2

DO
  DO i = 1, npr4
    x(i) = tfuel(i) + CKELVIN
  ENDDO
  
  DO i = 1, npr1
    kf(i) = CondFuel(x(i))
  ENDDO

  IF(lMOX) kf(1:npr1) = 0.9 * kf(1:npr1) 
      
  DO i = 1, npr
    kfm(i) = 0.5_8 * (kf(i) + kf(i+1))
  ENDDO
  
  DO i = npr2, npr4
    kf(i) = CondClad(x(i))
  ENDDO
  kmr = 0.5_8 * (kf(npr3) + kf(npr4))
  kml = 0.5_8 * (kf(npr3) + kf(npr2))
  
  x(1:npr4) = tfuel(1:npr4)
  xd(1:npr4) = x(1:npr4)
  !Matrix Set Up
  m = 1
  Diag(m) = 4._8 * kf(m)
  U(m) = -4*kf(m); !L(m) = 0; 
  b(m) = q
  IF (ThVar%rsi /= 0) b(m) = 0
  i = m
  DO m = 2, npr
    ri = 1/dble(i)
    Diag(m) = kfm(i) + kfm(m) + 0.5_8 * (kfm(m) - kfm(i)) * ri
    L(m) = - kfm(i) * (1._8 - 0.5_8 * ri); U(m) = - kfm(m) * (1._8 + 0.5_8 * ri)
    b(m) = q
    i = m
  ENDDO
  m = npr1
  alpha = kgap * (1._8 -kf(m-1)/kf(m))
  Diag(m) = 2.*(kf(m)+kgap*(one+half/npr))+alpha
  L(m) = -2.*kf(m); U(m) = -2.*kgap*(one+half/npr)-alpha
  b(m)=q
  
  m = npr2
  alpha = 2.*kgap2*(kf(m+1)/kf(m)-1)
  Diag(m) = 8.*kf(m)+kgap4-alpha
  L(m) = -kgap4+alpha
  U(m) =-8.*kf(m)
  b(m) = 0.
  
  m=npr3
  Diag(m) = 4.*(kmr+kml)+tworm*(kmr-kml)
  L(m) = -kml*(4-tworm)
  U(m) = -kmr*(4+tworm)
  B(m) = zero
  m=npr4
  kconv1 = htcoef * ThVar%tw
  alpha = 2.*kconv1*(one-kf(m-1)/kf(m))
  kconv = htcoef * tw * (4.+tw/rw)
  Diag(m) = 8.*kf(m)+kconv+alpha
  L(m) = -8.*kf(m); !U(m) = 0
  b(m) = (kconv+alpha) * tcool
  
  CALL DirLu1Dsolver(Diag(1:npr4), L(1:npr4), U(1:npr4), x(1:npr4), b(1:npr4), npr4)
  !Obtained temperature is C unit not k(because the bulk temp is C unit)
  !Check Error
  err = 0
  DO m = 1, npr4
    tfuel(m) = x(m)
    err = max(err, abs(x(m)-xd(m)) )
  ENDDO
  IF(err .lt. epsm2) EXIT
ENDDO
xvol = 0
call ftfavg(x, xvol(1:npr5), r, npr)
tfvol(1:npr2) = xvol(1:npr2) + CKELVIN
tfuel(npr5) = xvol(npr5)
!print*,'tfuel',tfuel
!print*,'tfvol',tfvol
NULLIFY(R)
END SUBROUTINE
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE Adjust_after_AA(tfvol, Tfuel, nr, ThVar, ThOpt)
  USE PARAM
  USE TYPEDEF,           ONLY : ThVar_Type,     ThOpt_Type
  USE FuelProperty_Mod,  ONLY : ftfavg
  IMPLICIT NONE
  REAL :: tfvol(nr-3)
  REAL :: Tfuel(nr)
  INTEGER :: nr, npr2, npr5, npr
  TYPE (ThVar_Type) :: ThVar
  TYPE (ThOPT_Type) :: ThOPT
  REAL :: xvol(100), x(100)
  REAL, POINTER :: R(:)

  npr = ThVar%npr
  npr2 = ThVar%npr2
  npr5 = ThVar%npr5
  R => ThVar%R
  xvol = 0._8; x = 0._8
  x(1:npr5)=Tfuel(1:npr5)
  call ftfavg(x, xvol(1:npr5), r, npr)
  tfvol(1:npr2) = xvol(1:npr2) + CKELVIN
  tfuel(npr5) = xvol(npr5)
  NULLIFY(R)
END SUBROUTINE
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE tfcalss_pwshape(tfvol, Tfuel, tcool, htcoef, qf, pwshape, nr, ThVar, ThOpt, lMox, hgap)
  USE PARAM
  USE TYPEDEF,           ONLY : ThVar_Type,     ThOpt_Type
  USE SteamTBL_mod,      ONLY : steamtbl
  USE FuelProperty_Mod,  ONLY : FTHCON,       CTHCON,       ftfavg,  RSGAPCOND
  USE LU_Mod,            ONLY : DirLu1Dsolver
  IMPLICIT NONE
  REAL :: tfvol(nr-3)
  REAL :: Tfuel(nr), tcool, htcoef, qf, pwshape(nr)
  INTEGER :: nr
  TYPE (ThVar_Type) :: ThVar
  TYPE (ThOPT_Type) :: ThOPT
  LOGICAL :: lMox
  !REAL, OPTIONAL :: hgap
  REAL :: hgap
  
  INTEGER :: hGapModel
  REAL :: diag(100), L(100), U(100), x(100), xd(100), b(100), xvol(100)
  REAL :: kf(100), kfb(100), kfm(100), kfmb(100)
  REAL :: kmr, kml, kgap, kgap2, kgap4, kconv, kconv1, tworm
  REAL, POINTER :: R(:)
  REAL :: delR, delR2, ri,  q, alpha
  REAL :: rs, rw, tw, rgap
  REAL :: err
  INTEGER :: npr, npr1, npr2, npr3, npr4, npr5
  INTEGER :: i, m
  
  hgapmodel = ThOpt%hGapModel
  
  npr = ThVar%npr; npr1 =ThVar%npr1
  npr2 = ThVar%npr2; npr3 = Thvar%npr3
  npr4= ThVar%npr4; npr5 = ThVar%npr5
  rs = ThVar%rs; rw = ThVar%rw
  tw = ThVar%tw; rgap = ThVar%rgap;
  R => ThVar%R
  DelR = R(2) - R(1); DelR2 = DelR * DelR
  hgap = ThOpt%hgap
  kgap = ThOpt%hgap * delr;
  kgap2 = ThOpt%hgap * tw * rs / rgap
  kgap4 = ThOpt%hgap * tw * (4._8 - tw / rgap) * rs / rgap
  tworm = tw / (rgap + 0.5_8 * tw)
  q = qf * DelR2
  
  DO
    DO i = 1, npr4
      x(i) = tfuel(i) + CKELVIN
    END DO
  
    DO i = 1, npr1
      kf(i) = FTHCON(x(i), 0.)
    END DO
  
    IF (lMOX) kf(1:npr1) = 0.9 * kf(1:npr1)
  
    DO i = 1, npr
      kfm(i) = 0.5_8 * (kf(i) + kf(i+1))
    END DO
  
    DO i = npr2, npr4
      kf(i) = CTHCON(x(i))
    END DO
    kmr = 0.5_8 * (kf(npr3) + kf(npr4))
    kml = 0.5_8 * (kf(npr3) + kf(npr2))
  
    !IF (PRESENT(hGap)) THEN
      IF (hGapModel .EQ. 3) THEN
        hgap = RSGAPCOND(kf(npr1), kf(npr2), x(npr1), x(npr2))
        kgap = hgap * delr;
        kgap2 = hgap * tw * rs / rgap
        kgap4 = hgap * tw * (4._8 - tw / rgap) * rs / rgap
      END IF
    !END IF
  
    x(1:npr4) = tfuel(1:npr4)
    xd(1:npr4) = x(1:npr4)
    !Matrix Set Up
    m = 1
    Diag(m) = 4._8 * kf(m)
    U(m) = -4*kf(m); !L(m) = 0;
    b(m) = q * pwshape(m)
    i = m
    DO m = 2, npr
      ri = 1/dble(i)
      Diag(m) = kfm(i) + kfm(m) + 0.5_8 * (kfm(m) - kfm(i)) * ri
      L(m) = - kfm(i) * (1._8 - 0.5_8 * ri); U(m) = - kfm(m) * (1._8 + 0.5_8 * ri)
      b(m) = q * pwshape(m)
      i = m
    END DO
    m = npr1
    alpha = kgap * (1._8 -kf(m-1)/kf(m))
    Diag(m) = 2.*(kf(m)+kgap*(one+half/npr))+alpha
    L(m) = -2.*kf(m); U(m) = -2.*kgap*(one+half/npr)-alpha
    b(m)=q * pwshape(m-1)
  
    m = npr2
    alpha = 2.*kgap2*(kf(m+1)/kf(m)-1)
    Diag(m) = 8.*kf(m)+kgap4-alpha
    L(m) = -kgap4+alpha
    U(m) =-8.*kf(m)
    b(m) = 0.
  
    m=npr3
    Diag(m) = 4.*(kmr+kml)+tworm*(kmr-kml)
    L(m) = -kml*(4-tworm)
    U(m) = -kmr*(4+tworm)
    B(m) = zero
    m=npr4
    kconv1 = htcoef * ThVar%tw
    alpha = 2.*kconv1*(one-kf(m-1)/kf(m))
    kconv = htcoef * tw * (4.+tw/rw)
    Diag(m) = 8.*kf(m)+kconv+alpha
    L(m) = -8.*kf(m); !U(m) = 0
    b(m) = (kconv+alpha) * tcool
  
    CALL DirLu1Dsolver(Diag(1:npr4), L(1:npr4), U(1:npr4), x(1:npr4), b(1:npr4), npr4)
    !Obtained temperature is C unit not k(becaUSE the bulk temp is C unit)
    !Check Error
    err = 0
    DO m = 1, npr4
      tfuel(m) = x(m)
      err = max(err, abs(x(m)-xd(m)) )
    END DO
    IF (err .lt. epsm2) EXIT
  END DO
  xvol = 0
  call ftfavg(x, xvol(1:npr5), r, npr)
  tfvol(1:npr2) = xvol(1:npr2) + CKELVIN
  tfuel(npr5) = xvol(npr5)
  NULLIFY(R)
  END SUBROUTINE
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE tfcalss_qshape(tfvol, Tfuel, tcool, htcoef, qf, qprime, qshape, BurnUp, nr, ThVar, ThOpt, lMOX)
USE PARAM
USE TYPEDEF,           ONLY : ThVar_Type,     ThOpt_Type
USE SteamTBL_mod,      ONLY : steamtbl
USE FuelProperty_Mod,  ONLY : CondFuel,       CondFuel_FRAPCON,       CondClad,       ftfavg,       &
                              GapCon,         GapCon_ROPER
USE LU_Mod,            ONLY : DirLu1Dsolver
IMPLICIT NONE
REAL :: tfvol(nr-3), qshape(nr), BurnUp(0:nr)
REAL :: Tfuel(nr), tcool, htcoef, qf, qprime
INTEGER :: nr
TYPE (ThVar_Type) :: ThVar
TYPE (ThOPT_Type) :: ThOPT
LOGICAL :: lMox

REAL :: diag(100), L(100), U(100), x(100), xd(100), b(100), xvol(100)
REAL :: kf(100), kfb(100), kfm(100), kfmb(100)
REAL :: kmr, kml, kgap, kgap2, kgap4, kconv, kconv1, tworm
REAL, POINTER :: R(:)
REAL :: delR, delR2, ri,  q, alpha
REAL :: rs, rw, tw, rgap
REAL :: err
REAL :: hgap0, burnupavg
INTEGER :: npr, npr1, npr2, npr3, npr4, npr5
INTEGER :: i, m

npr = ThVar%npr; npr1 =ThVar%npr1
npr2 = ThVar%npr2; npr3 = Thvar%npr3
npr4= ThVar%npr4; npr5 = ThVar%npr5
rs = ThVar%rs; rw = ThVar%rw
tw = ThVar%tw; rgap = ThVar%rgap;
R => ThVar%R

hgap0 = ThOpt%hgap
ThOpt%hgap = GapCon_ROPER(qprime, BurnUp(0))

DelR = R(2) - R(1); DelR2 = DelR * DelR
kgap = ThOpt%hgap * delr;
kgap2 = ThOpt%hgap * tw * rs / rgap
kgap4 = ThOpt%hgap * tw * (4._8 - tw / rgap) * rs / rgap
tworm = tw / (rgap + 0.5_8 * tw)
q = qf * DelR2

DO
  DO i = 1, npr4
    x(i) = tfuel(i) + CKELVIN
  END DO

  DO i = 1, npr1
    kf(i) = CondFuel(x(i))
    kf(i) = CondFuel_FRAPCON(x(i), Burnup(i))
    continue
  END DO

  IF (lMOX) kf(1:npr1) = 0.9 * kf(1:npr1)

  DO i = 1, npr
    kfm(i) = 0.5_8 * (kf(i) + kf(i+1))
  END DO

  DO i = npr2, npr4
    kf(i) = CondClad(x(i))
  END DO
  kmr = 0.5_8 * (kf(npr3) + kf(npr4))
  kml = 0.5_8 * (kf(npr3) + kf(npr2))

  x(1:npr4) = tfuel(1:npr4)
  xd(1:npr4) = x(1:npr4)
  !Matrix Set Up
  m = 1
  Diag(m) = 4._8 * kf(m)
  U(m) = -4*kf(m); !L(m) = 0;
  b(m) = qf * qshape(1) * DelR2
  i = m
  DO m = 2, npr
    ri = 1/dble(i)
    Diag(m) = kfm(i) + kfm(m) + 0.5_8 * (kfm(m) - kfm(i)) * ri
    L(m) = - kfm(i) * (1._8 - 0.5_8 * ri); U(m) = - kfm(m) * (1._8 + 0.5_8 * ri)
    b(m) = qf * qshape(m) * DelR2
    i = m
  END DO
  m = npr1
  alpha = kgap * (1._8 -kf(m-1)/kf(m))
  Diag(m) = 2.*(kf(m)+kgap*(one+half/npr))+alpha
  L(m) = -2.*kf(m); U(m) = -2.*kgap*(one+half/npr)-alpha
  b(m)= qf * qshape(m) * DelR2

  m = npr2
  alpha = 2.*kgap2*(kf(m+1)/kf(m)-1)
  Diag(m) = 8.*kf(m)+kgap4-alpha
  L(m) = -kgap4+alpha
  U(m) =-8.*kf(m)
  b(m) = 0.

  m=npr3
  Diag(m) = 4.*(kmr+kml)+tworm*(kmr-kml)
  L(m) = -kml*(4-tworm)
  U(m) = -kmr*(4+tworm)
  B(m) = zero
  m=npr4
  kconv1 = htcoef * ThVar%tw
  alpha = 2.*kconv1*(one-kf(m-1)/kf(m))
  kconv = htcoef * tw * (4.+tw/rw)
  Diag(m) = 8.*kf(m)+kconv+alpha
  L(m) = -8.*kf(m); !U(m) = 0
  b(m) = (kconv+alpha) * tcool

  CALL DirLu1Dsolver(Diag(1:npr4), L(1:npr4), U(1:npr4), x(1:npr4), b(1:npr4), npr4)
  !Obtained temperature is C unit not k(becaUSE the bulk temp is C unit)
  !Check Error
  err = 0
  DO m = 1, npr4
    tfuel(m) = x(m)
    err = max(err, abs(x(m)-xd(m)) )
  END DO
  IF (err .lt. epsm2) EXIT
END DO

!ThOpt%hgap = hgap0
xvol = 0
call ftfavg(x, xvol(1:npr5), r, npr)
tfvol(1:npr2) = xvol(1:npr2) + CKELVIN
tfuel(npr5) = xvol(npr5)
NULLIFY(R)

END SUBROUTINE
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetBurnupDist(BurnupDist, Core, FmInfo, ixy, nz, nr, Thvar, ThOpt, PE)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type,  FmInfo_Type,     ThVar_Type,     ThOpt_Type,           &
                    PE_TYPE,                                                               &
                    FxrInfo_Type, Pin_Type,       Cell_Type
USE UtilFunction, ONLY : LineIntPol
#ifdef MPI_ENV
USE MPIComm_Mod,  ONLY : REDUCE
#endif
IMPLICIT NONE
INTEGER :: ixy, nz, nr
REAL :: BurnupDist(0:nr, nz)        !  Heat source Shape
TYPE (CoreInfo_Type) :: Core   !  Core Geometry Infomation
TYPE (FmInfo_Type) :: FmInfo   !  Fine Mesh Information
TYPE (ThVar_Type) :: ThVar     !  Th related parameters
TYPE (ThOpt_TYPE) :: ThOpt     !  Th Option
TYPE (PE_TYPE) :: PE

TYPE (Cell_Type), POINTER :: CellInfo(:)
TYPE (Pin_Type), POINTER :: Pin(:)

INTEGER :: iz, ir,  icel, i, j, k, ifxr, ifsr, ifsrlocal
INTEGER :: nCircle, FsrIdxSt, St, FxrIdxSt, nLocalFxr, nFsrInFxr
REAL :: ThR(100), R(100), BurnUpDist0(100), xdat(100), ydat(100)
REAL :: burnupavg, area
REAL :: Buf(0:nr, nz)
LOGICAL :: lIntPol

lIntPol = .TRUE.

CellInfo => Core%CellInfo
Pin => Core%Pin
FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt

IF (ThVar%rsi /= 0) THEN
  ThR(1) = 0.
  ThR(2) = ThVar%rsi*100
  DO ir = 3, ThVar%npr1
    ThR(ir) = ThVar%rsi*100 + (ir -2) * ThVar%FuelDelR*100
  ENDDO
ELSE
  DO ir = 1, ThVar%npr1
    ThR(ir) = (ir -1) * ThVar%FuelDelR*100
  ENDDO
ENDIF

DO iz = 1, nz
  BurnupDist(:, iz) = 0
  IF (iz .LT. PE%myzb .or. iz .GT. PE%myze) CYCLE
  icel = Pin(ixy)%Cell(iz)
  IF (.NOT. CellInfo(icel)%lfuel) CYCLE
  IF (CellInfo(icel)%lRect) CYCLE
  IF (CellInfo(icel)%lCad) CYCLE

  nLocalFxr = CellInfo(icel)%nFxr
  area = 0;burnupavg = 0
  DO j = CellInfo(icel)%ThCell%FuelReg(1), CellInfo(icel)%ThCell%FuelReg(2)
    ifxr = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    BurnUpDist0(j) = FmInfo%Fxr(ifxr, iz)%burnup
    burnupavg = burnupavg + FmInfo%Fxr(ifxr, iz)%burnup * FmInfo%Fxr(ifxr, iz)%area
    area = area + FmInfo%Fxr(ifxr, iz)%area
  END DO
  burnupavg = burnupavg / area
  BurnupDist(0, iz) = burnupavg
  !pwsum= sum(power0(1:nLocalFxr))
  R = 0
  nCircle = CellInfo(icel)%geom%nCircle
  DO j = 1, nCircle
    R(j) = CellInfo(icel)%Geom%Circle(3, j)
  END DO
  k = 0
  DO j = CellInfo(icel)%ThCell%FuelReg(2), CellInfo(icel)%ThCell%FuelReg(1), -1
    i = CellInfo(icel)%ThCell%FuelReg(2)-j + 1
    k = k + 1
    xdat(k) = 0.5 * (R(j) + R(j-1))
    ydat(k) = BurnUpDist0(j)

  END DO
  IF (ThVar%rsi /= 0) THEN
    DO i=k+1,2,-1
      xdat(i) = xdat(i-1) 
      ydat(i) = ydat(i-1)
    ENDDO
    xdat(1) = CellInfo(icel)%Geom%Circle(3, nCircle)
    ydat(1) = 0.
    k = k+1
  ENDIF
  DO i = 1, ThVar%npr1
    IF (.NOT. lIntPol) THEN  !Step Power Shape
      DO j = CellInfo(icel)%ThCell%FuelReg(2), CellInfo(icel)%ThCell%FuelReg(1), -1
        IF (ThR(i) .GE. R(j) .AND. ThR(i) .LE. R(j-1)) EXIT
      END DO
      j=MAX(j, CellInfo(icel)%ThCell%FuelReg(1))
      BurnUpDist(i, iz) = BurnUpDist0(j)
    ELSE  !Linear Interpolation Power Shape
      BurnUpDist(i, iz) = LineIntPol(ThR(i), k, xdat(1:k), ydat(1:k))
    END IF
    CONTINUE

  END DO

  CONTINUE
END DO
#ifdef MPI_ENV
CALL REDUCE(BurnupDist(0:nr, 1:nz), Buf(0:nr, 1:nz), nr+1, nz, PE%MPI_CMFD_COMM, TRUE)
BurnupDist(0:nr, 1:nz) = Buf(0:nr, 1:nz)
#endif
NULLIFY(CellInfo); NULLIFY(Pin)


END SUBROUTINE
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetPwShape(Pwshape, Core, Fxr, Power, ixy, nz, nr, Thvar, ThOpt, PE)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type,  FmInfo_Type,     ThVar_Type,     ThOpt_Type,           &
                    PE_Type,        &
                    FxrInfo_Type, Pin_Type,       Cell_Type
USE UtilFunction, ONLY : LineIntPol
USE MPIComm_Mod,  ONLY : REDUCE
IMPLICIT NONE
INTEGER :: ixy, nz, nr
REAL :: Pwshape(nr, nz)        !  Heat source Shape
TYPE (CoreInfo_Type) :: Core   !  Core Geometry Infomation
!TYPE (FmInfo_Type) :: FmInfo   !  Fine Mesh Information
TYPE (FxrInfo_Type), POINTER :: Fxr(:, :)
REAL, POINTER :: Power(:, :)
TYPE (ThVar_Type) :: ThVar     !  Th related parameters
TYPE (ThOpt_TYPE) :: ThOpt     !  Th Option
TYPE (PE_TYPE) :: PE



TYPE (Cell_Type), POINTER :: CellInfo(:)
TYPE (Pin_Type), POINTER :: Pin(:)

INTEGER :: iz, ir,  icel, i, j, k, ifxr, ifsr, ifsrlocal, ir1, ir2
INTEGER :: nCircle, FsrIdxSt, St, FxrIdxSt, nLocalFxr, nFsrInFxr
REAL :: ThR(100), R(100), Power0(100), xdat(100), ydat(100)
REAL :: Buf(nr, nz)
REAL :: vol, pwsum, f

LOGICAL :: lIntPol

lIntPol = .FALSE.

CellInfo => Core%CellInfo
Pin => Core%Pin
FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt

IF (ThVar%rsi /= 0) THEN
  ThR(1) = 0.
  ThR(2) = ThVar%rsi*100
  DO ir = 3, ThVar%npr1
    ThR(ir) = ThVar%rsi*100 + (ir -2) * ThVar%FuelDelR*100
  ENDDO
ELSE
  DO ir = 1, ThVar%npr1
    ThR(ir) = (ir -1) * ThVar%FuelDelR*100
  ENDDO
ENDIF
Pwshape = 0
DO iz = PE%myzb, PE%myze
  !Pwshape(:, iz) = 0
  !Pwshape(1:ThVar%npr1, iz) = 1._8
  IF (iz .LT. PE%myzb .OR. iz .GT. PE%myze) CYCLE
  icel = Pin(ixy)%Cell(iz)
  IF (.NOT. CellInfo(icel)%lfuel) CYCLE
  IF (CellInfo(icel)%lRect) CYCLE
  IF (CellInfo(icel)%lCad) CYCLE

  !Power Shape Normalization => Make Pin Cell Average power to unity
  nLocalFxr = CellInfo(icel)%nFxr
  Power0 = 0; pwsum =0; vol = 0
  DO j = CellInfo(icel)%ThCell%FuelReg(1), CellInfo(icel)%ThCell%FuelReg(2)
    ifxr = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    DO i = 1, nFsrInFxr
      ifsrlocal = Cellinfo(icel)%MapFxr2FsrIdx(i, j)
      ifsr = FsrIdxSt + ifsrlocal - 1
      Power0(j) = Power0(j) + CellInfo(icel)%vol(ifsrlocal)*Power(ifsr, iz)
    END DO
    IF (Power0(j) .GT. 0._8) vol = vol + Fxr(ifxr, iz)%area
    Pwsum = Pwsum + Power0(j)
    Power0(j) = Power0(j) / Fxr(ifxr, iz)%area
  END DO
  !pwsum= sum(power0(1:nLocalFxr))
  f = vol / pwsum
  power0 = power0 * f
  R = 0
  nCircle = CellInfo(icel)%geom%nCircle
  DO j = 1, nCircle
    R(j) = CellInfo(icel)%Geom%Circle(3, j)
  END DO
  k = 0
  DO j = CellInfo(icel)%ThCell%FuelReg(2), CellInfo(icel)%ThCell%FuelReg(1), -1
    i = CellInfo(icel)%ThCell%FuelReg(2)-j + 1
    k = k + 1
    xdat(k) = (R(j-1))
    ydat(k) = Power0(j)
  END DO
  IF (ThVar%rsi /= 0) THEN
    DO i=k+1,2,-1
      xdat(i) = xdat(i-1) 
      ydat(i) = ydat(i-1)
    ENDDO
    xdat(1) = CellInfo(icel)%Geom%Circle(3, nCircle)
    ydat(1) = 0.
    k = k+1
  ENDIF
  CONTINUE
  ir1 = 1
  DO i = 1, ThVar%npr
    DO j= 1, k
      IF (ThR(i+1) .LE. xdat(j)) EXIT
    END DO
    ir2 = min(j, k)
    IF (ir1 .EQ. ir2) THEN
      PwShape(i, iz) = ydat(ir2)
    ELSE
      PwSum = (xdat(ir1)**2 - ThR(i)**2) * ydat(ir1)
      PwSum = Pwsum + (ThR(i+1)**2 - xdat(ir2-1)**2) * ydat(ir2)
      DO j = ir1 + 1, ir2-1
        PwSum =PwSum + (xdat(j)**2 - xdat(j-1)**2) * ydat(j)
      END DO
      PwSum = PwSum / (ThR(i+1)**2 - ThR(i)**2)
      PwShape(i, iz) = PwSum
    END IF
    ir1 = ir2
  END DO

END DO

#ifdef MPI_ENV
CALL REDUCE(PwShape, Buf, nr, nz, PE%MPI_CMFD_COMM, .TRUE.)
PwShape = Buf
#endif

NULLIFY(CellInfo); NULLIFY(Pin)
END SUBROUTINE
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetEffectiveTemp(tfvol, Tfuel,  Burnupavg, nr, ThVar, ThOpt)
USE PARAM
USE TYPEDEF,           ONLY : ThVar_Type,     ThOpt_Type
IMPLICIT NONE
REAL :: tfvol(nr-3)
REAL :: Tfuel(nr)
REAL :: BurnupAvg
INTEGER :: nr
TYPE (ThVar_Type) :: ThVar
TYPE (ThOPT_Type) :: ThOPT
CHARACTER(10) :: warg

REAL, SAVE :: w
LOGICAL, SAVE :: lFirst
INTEGER :: npr, npr1, npr2, npr3, npr4, npr5
DATA lFirst /.TRUE./
!IF (lFirst) THEN
!  call getarg(2,warg)
!  READ(warg, *) w
!END IF
w =0.57;
npr = ThVar%npr; npr1 =ThVar%npr1
npr2 = ThVar%npr2; npr3 = Thvar%npr3
npr4= ThVar%npr4; npr5 = ThVar%npr5

!tfvol(1:npr2) = tfuel(npr5) + CKELVIN
CALL GetOptEffTempW(w, Burnupavg)
!tfvol(1:npr2) = (1._8- w)*tfuel(1)+w*tfuel(npr1)+CKELVIN
tfvol(1:npr2) = (w)*tfuel(npr5)+(1-w)*tfuel(npr1)+CKELVIN
CONTINUE
END SUBROUTINE
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE CalEffectiveTemp(EffTemp, tfvol, Tfuel, Burnupavg, nr,  ThVar, ThOpt)
USE PARAM
USE TYPEDEF,           ONLY : ThVar_Type,     ThOpt_Type
IMPLICIT NONE
REAL::EffTemp
REAL :: tfvol(nr-3)
REAL :: Tfuel(nr)
REAL :: Burnupavg
INTEGER :: nr
TYPE (ThVar_Type) :: ThVar
TYPE (ThOPT_Type) :: ThOPT
CHARACTER(10) :: warg

REAL, SAVE :: w
LOGICAL, SAVE :: lFirst
INTEGER :: npr, npr1, npr2, npr3, npr4, npr5
DATA lFirst /.TRUE./
!IF (lFirst) THEN
!  call getarg(2,warg)
!  READ(warg, *) w
!END IF
w =0.57;
npr = ThVar%npr; npr1 =ThVar%npr1
npr2 = ThVar%npr2; npr3 = Thvar%npr3
npr4= ThVar%npr4; npr5 = ThVar%npr5

!tfvol(1:npr2) = tfuel(npr5) + CKELVIN

!tfvol(1:npr2) = (1._8- w)*tfuel(1)+w*tfuel(npr1)+CKELVIN
!Get Optimum Weighting Factor
CALL GetOptEffTempW(w, Burnupavg)
EffTemp = (w)*tfuel(npr5)+(1-w)*tfuel(npr1)+CKELVIN
CONTINUE
END SUBROUTINE
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE GetOptEffTempW(W, Burnup)
USE PARAM
USE UtilFunction,    ONLY : LineIntPol
IMPLICIT NONE
REAL :: W, Burnup
REAL :: BurnupPt(10), Wopt(10)
DATA BurnupPt / 0.15_8, 5.0_8, 10.0_8, 15.0_8, 20.0_8, 25.0_8, 30.0_8, 35.0_8, 40.0_8, 45.0_8 /
DATA Wopt /0.865_8, 0.848_8, 0.834_8, 0.828_8, 0.818_8, 0.812_8, 0.810_8, 0.804_8, 0.801_8, 0.799_8/

W=LineIntPol(Burnup, 10, BurnupPt, Wopt)
END SUBROUTINE

SUBROUTINE WriteThDist(Pwshape, Core, FmInfo, FuelTH, ixy, nz, nr, Thvar, ThOpt, PE)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type,  FmInfo_Type,     ThVar_Type,     ThOpt_Type,           &
                    PE_Type,        FuelTH_Type,                                &
                    FxrInfo_Type, Pin_Type,       Cell_Type
USE UtilFunction, ONLY : LineIntPol
IMPLICIT NONE
INTEGER :: ixy, nz, nr
REAL :: Pwshape(nr, nz)        !  Heat source Shape
TYPE (CoreInfo_Type) :: Core   !  Core Geometry Infomation
TYPE (FmInfo_Type) :: FmInfo   !  Fine Mesh Information
TYPE (FuelTh_Type) :: FuelTH
TYPE (ThVar_Type) :: ThVar     !  Th related parameters
TYPE (ThOpt_TYPE) :: ThOpt     !  Th Option
TYPE (PE_TYPE) :: PE

TYPE (Cell_Type), POINTER :: CellInfo(:)
TYPE (Pin_Type), POINTER :: Pin(:)

INTEGER :: iz, ir,  icel, i, j, k, ifxr, ifsr, ifsrlocal
INTEGER :: nCircle, FsrIdxSt, St, FxrIdxSt, nLocalFxr, nFsrInFxr
REAL :: ThR(100), R(100), Power0(100), xdat(100), ydat(100)
REAL :: vol, pwsum, f

LOGICAL :: lIntPol

lIntPol = .TRUE.

CellInfo => Core%CellInfo
Pin => Core%Pin
FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt

DO ir = 1, ThVar%npr1
  ThR(ir) = (ir -1) * ThVar%FuelDelR*100
END DO
Pwshape = 1
  iz = 22
  IF (ixy .NE. 20) RETURN
  !Pwshape(:, iz) = 0
  !Pwshape(1:ThVar%npr1, iz) = 1._8
  IF (iz .LT. PE%myzb .OR. iz .GT. PE%myze) RETURN
  icel = Pin(ixy)%Cell(iz)
  IF (.NOT. CellInfo(icel)%lfuel) RETURN
  IF (CellInfo(icel)%lRect) RETURN
  IF (CellInfo(icel)%lCad) RETURN

  !Power Shape Normalization => Make Pin Cell Average power to unity
  nLocalFxr = CellInfo(icel)%nFxr
  Power0 = 0; pwsum =0; vol = 0
  DO j = CellInfo(icel)%ThCell%FuelReg(1), CellInfo(icel)%ThCell%FuelReg(2)
    ifxr = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    DO i = 1, nFsrInFxr
      ifsrlocal = Cellinfo(icel)%MapFxr2FsrIdx(i, j)
      ifsr = FsrIdxSt + ifsrlocal - 1
      Power0(j) = Power0(j) + CellInfo(icel)%vol(ifsrlocal)*FmInfo%Power(ifsr, iz)
    END DO
    IF (Power0(j) .GT. 0._8) vol = vol + FmInfo%Fxr(ifxr, iz)%area
    Pwsum = Pwsum + Power0(j)
    Power0(j) = Power0(j) / FmInfo%Fxr(ifxr, iz)%area
  END DO
  !DO j = CellInfo(icel)%ThCell%FuelReg(1), CellInfo(icel)%ThCell%FuelReg(2)
    WRITE(101, '(1000F15.4)') (Power0(j), j = CellInfo(icel)%ThCell%FuelReg(1), CellInfo(icel)%ThCell%FuelReg(2))
    !WRITE(101, '(1000F15.4)') FuelTH%Tfuel(:, iz)
    !WRITE(101, *)
  !END DO
NULLIFY(CellInfo); NULLIFY(Pin)
END SUBROUTINE
! ------------------------------------------------------------------------------------------------------------