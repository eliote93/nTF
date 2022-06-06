SUBROUTINE SimpleTH(powlin, plevel, Tout, Tfmax, RelPW, Profile, FuelTh, CoolantTH, Thvar, ThOpt, PE)
  
USE PARAM
USE TYPEDEF,           ONLY : CoolantTh_Type, FuelTh_Type, ThVar_Type, ThOpt_Type, PE_TYPE, Pin_Type, Cell_Type
USE Water155_mod,      ONLY : fhcap, fdens, fhtcoef
USE FuelProperty_Mod,  ONLY : CondClad, CondFuel

IMPLICIT NONE

REAL :: PowLin, Plevel, Tout, Tfmax
REAL, POINTER :: Profile(:, :)
REAL :: RelPW(:)
TYPE (FuelTh_Type) :: FuelTh
TYPE (CoolantTh_Type) :: CoolantTH
TYPE (Thvar_Type) :: ThVar
TYPE (THOpt_Type) :: ThOPT
TYPE (PE_TYPE) :: PE
! ----------------------------------------------------
TYPE (Cell_Type), POINTER :: Cell
REAL, POINTER, DIMENSION(:) :: hz, tcool, DenCool, rhou
REAL, POINTER, DIMENSION(:,:) :: tfvol, tfuel

REAL :: qprime, qvolk, cp, htcoefk, hgap, fkck, fkfk, tink, toutk, tavg, tfbar, pFracSum, vFracSum, pFracRing, vFracRing, fmdotfa, deq, twopi, rs, rw, tw, rgap, alogval
INTEGER :: npr, npr1, npr2, npr3, npr4, npr5, iz, ir, icel
! ----------------------------------------------------

tink = ThVar%Tin

fmdotfa = ThVar%MdotFA
deq     = ThVar%deq
hgap    = ThOpt%hgap
npr     = ThVar%npr
npr1    = ThVar%npr1
npr2    = ThVar%npr2
npr3    = Thvar%npr3
npr4    = ThVar%npr4
npr5    = ThVar%npr5
tink    = ThVar%Tin
rs      = ThVar%rs
rw      = ThVar%rw
tw      = ThVar%tw
rgap    = ThVar%rgap
twopi   = 2._8*pi

hz      => ThVar%hz
Tcool   => CoolantTH%Tcool
DenCool => CoolantTH%DenCool
rhou    => CoolantTH%rhou
tfuel   => FuelTh%tfuel
Tfvol   => FuelTh%Tfvol

tfmax = 0.

DO iz = 1, ThVar%nzth
  alogval = log(rw/rs)
  
  ! Cool.
  qprime  = plevel * powlin * RelPw(iz)
  qvolk   = qprime * hz(iz) ! Vol. Heat E.
  cp      = fhcap(tcool(iz))
  toutk   = tink + qvolk/(fmdotfa*cp) ! Out Temp.
  
  tcool  (iz) = 0.5_8 * (toutk + tink) ! Avg. Coo. Temp.
  DenCool(iz) = fdens(tcool(iz))       ! Avg. Coo. Dens.
  tink        = toutk
  
  ! Conduction : Clad.
  htcoefk = fhtcoef(tcool(iz), deq, rhou(0))
  tfuel(npr4, iz) = tcool(iz) + qprime/(twopi * rw * htcoefk)     ! Clad. Outer Wall
  fkck = CondClad(tfuel(npr3, iz)+CKELVIN)
  tfuel(npr2, iz) = tfuel(npr4, iz) + alogval*qprime/(twopi*fkck) ! Clad. Inner Wall
  tfuel(npr3, iz) = 0.5*(tfuel(npr2, iz)+tfuel(npr4, iz))         ! Clad Cnt.
  
  ! Conduction : Fuel
  tfuel(npr1, iz) = tfuel(npr2, iz) + qprime / (twopi*rs*hgap) ! Fuel Pellet Sfc.
  
  tavg     = 0.
  pFracSum = 0.
  vFracSum = 0.
  
  DO ir = npr, 1, -1
    tfbar = (tfuel(ir, iz) + tfuel(ir + 1, iz)) * 0.5 ! Use : Prv. tfuel(ir, iz)
    fkfk  = CondFuel(tfbar + CKELVIN)
    tfuel(ir, iz) = tfuel(ir+1, iz) + qprime * profile(ir, iz) / (2*twopi * fkfk) ! Upd. : tfuel (ir, iz)
    tavg = tavg + Profile(npr+ir, iz) * 0.5 * (tfuel(ir, iz)+tfuel(ir+1, iz))
  END DO
  tfuel(npr5, iz) = tavg ! Avg. Fuel Temp.
  
  DO ir = 1, npr
    Tfvol(ir, iz) = 0.5 * (tfuel(ir, iz) + tfuel(ir+1, iz)) + CKELVIN
  END DO
  
  Tfvol(npr1, iz) = 0.5 * (tfuel(npr1, iz) + tfuel(npr2, iz)) + CKELVIN ! Gap
  Tfvol(npr2, iz) = 0.5 * (tfuel(npr2, iz) + tfuel(npr4, iz)) + CKELVIN ! Clad Cnt.
  
  tfmax = max(tfmax, tfuel(1, iz)) ! Fuel Inner Most Temp.
  
  IF (tfmax .NE. tfmax) THEN
    STOP
  END IF
END DO

Tout = toutk

NULLIFY (hz)
! ----------------------------------------------------

END SUBROUTINE SimpleTH