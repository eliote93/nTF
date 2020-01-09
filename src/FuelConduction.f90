SUBROUTINE TranRadFuelCondFDM(tfvol, TFuel, TFuelD, Tcool, Tcoold, htcoef, htcoefd, qf, qfd, QShape, Qshaped, nr, ThVar, ThOpt)
USE PARAM
USE TYPEDEF,           ONLY : ThVar_Type,     ThOpt_Type
USE SteamTBL_mod,      ONLY : steamtbl
USE LU_Mod,            ONLY : DirLu1Dsolver
USE FuelProperty_Mod,  ONLY : CladCp, FuelCp
USE TH_MOD,            ONLY : CalWallTemp,    CalGapTemp,     CalAvgFuelTemp
IMPLICIT NONE

REAL :: tfvol(nr-3)
REAL :: Tfuel(nr), Tfueld(nr), QShape(nr), QShaped(nr)
REAL :: tcoold, tcool, htcoef, htcoefd, qf, qfd
INTEGER :: nr
TYPE(ThVar_Type) :: ThVar
TYPE(ThOPT_Type) :: ThOPT

REAL :: cetaf, cetafb, cetah
REAL :: diag(100), L(100), U(100), x(100), xd(100), b(100), btran(100)
REAL :: Cp(100)
REAL, POINTER :: R(:)
REAL :: rs, rw, tw, rgap, delt
REAL :: twall, Tgap(2)
REAL :: err
INTEGER :: npr, npr1, npr2, npr3, npr4, npr5
INTEGER :: i, m

cetaf = 0.5_8; cetafb = 1._8 - cetaf
cetah = cetafb / cetaf
DelT = ThVar%DelT
npr = ThVar%npr; npr1 =ThVar%npr1
npr2 = ThVar%npr2; npr3 = Thvar%npr3
npr4= ThVar%npr4; npr5 = ThVar%npr5
rs = ThVar%rs; rw = ThVar%rw
tw = ThVar%tw; rgap = ThVar%rgap; 
R => ThVar%R

DO i=1, npr
  xd(i) = Tfueld(i)
ENDDO
i=npr+1; xd(i) = TFueld(i+1)
i=npr+2; xd(i) = TFueld(i+1)
xd(1:npr2) = xd(1:npr2) + CKELVIN
CALL SetFuelCondLs(Diag, L, U, b, xd, tcoold, htcoefd, qfd, qshaped, nr, ThVar, ThOpt, FALSE)
xd(1:npr2) = xd(1:npr2) - CKELVIN

i = 1
btran(i) = b(i) - (Diag(i) * xd(i) + U(i) * xd(i+1))
DO i = 2, npr2 -1
  btran(i) = b(i) - (L(i) * xd(i-1) + Diag(i) * xd(i)  + U(i) * xd(i+1)) 
ENDDO
i = npr2
btran(i) = b(i) - (L(i) * xd(i-1) + Diag(i) * xd(i))


DO i=1, npr
  x(i) = Tfuel(i)
ENDDO
i=npr+1; x(i) = TFuel(i+1); i=npr+2; x(i) = TFuel(i+1)
x(1:npr2) = x(1:npr2) + CKELVIN

CALL SetFuelCondLs(Diag, L, U, b, x, tcool, htcoef, qf, qshape, nr, ThVar, ThOpt, FALSE)

DO i = 1, npr
  CP(i) =   FuelCp(x(i)) / (DelT*cetaf) * (R(i+1)**2-R(i)**2) / 2
ENDDO
DO i = npr1, npr2
  CP(i) =   CladCp(x(i)) / (DelT*cetaf) * (R(i+2)**2-R(i+1)**2) / 2
ENDDO
x(1:npr2) = x(1:npr2) - CKELVIN

DO i = 1, npr2
  btran(i) = Cetah * btran(i) + b(i) + CP(i) * xd(i)
  Diag(i) = Diag(i) + CP(i)
ENDDO

CALL DirLu1Dsolver(Diag(1:npr2), L(1:npr2), U(1:npr2), x(1:npr2), btran(1:npr2), npr2)

Twall = CalWallTemp(X, Tcool, htcoef, ThVar, ThOpt)
Tgap = CalGapTemp(X, ThVar, ThOpt)
Tfuel(1:npr) = X(1:npr)
Tfuel(npr1) = Tgap(1)
Tfuel(npr2:npr3) = x(npr1:npr2)
Tfuel(npr4) = Twall
Tfuel(npr5) = CalAvgFuelTemp(X, ThVar, ThOpt)
tfvol(1:npr) = X(1:npr) + CKELVIN
tfvol(npr1) = 0.5 * (Tgap(1) + Tgap(2)) + CKELVIN
tfvol(npr2) = 0.5 * (x(npr1) + x(npr2)) + CKELVIN


END SUBROUTINE

  
SUBROUTINE RadFuelCondFDM(tfvol, Tfuel, tcool, htcoef, qf, qshape, nr, ThVar, ThOpt, lMox)
USE PARAM
USE TYPEDEF,           ONLY : ThVar_Type,     ThOpt_Type
USE SteamTBL_mod,      ONLY : steamtbl
USE LU_Mod,            ONLY : DirLu1Dsolver
USE TH_MOD,            ONLY : CalWallTemp,    CalGapTemp,     CalAvgFuelTemp
IMPLICIT NONE
REAL :: tfvol(nr-3), Qshape(nr)
REAL :: Tfuel(nr), tcool, htcoef, qf
INTEGER :: nr
TYPE(ThVar_Type) :: ThVar
TYPE(ThOPT_Type) :: ThOPT
LOGICAL :: lMox

REAL :: diag(100), L(100), U(100), x(100), xd(100), b(100)
REAL, POINTER :: R(:)
REAL :: hg, hgt
REAL :: rs, rw, tw, rgap
REAL :: twall, Tgap(2)
REAL :: err
INTEGER :: npr, npr1, npr2, npr3, npr4, npr5
INTEGER :: i, m

npr = ThVar%npr; npr1 =ThVar%npr1
npr2 = ThVar%npr2; npr3 = Thvar%npr3
npr4= ThVar%npr4; npr5 = ThVar%npr5
rs = ThVar%rs; rw = ThVar%rw
tw = ThVar%tw; rgap = ThVar%rgap; 
R => ThVar%R

hg = ThOpt%hgap; hgt = hg * rs / rgap 

DO i=1, npr
  x(i) = Tfuel(i)
ENDDO
i=npr+1; x(i) = TFuel(i+1)
i=npr+2; x(i) = TFuel(i+1)
DO
  x(1:npr2) = x(1:npr) + CKELVIN
  CALL SetFuelCondLs(Diag, L, U, b, x, tcool, htcoef, qf, qshape, nr, ThVar, ThOpt, lMox)
  x(1:npr2) = x(1:npr) - CKELVIN
  xd(1:npr2) = x(1:npr)
  CALL DirLu1Dsolver(Diag(1:npr2), L(1:npr2), U(1:npr2), x(1:npr2), b(1:npr2), npr2)
  err = 0
  DO m = 1, npr4
    err = max(err, abs(x(m)-xd(m)) )
  ENDDO
  IF(err .lt. epsm2) EXIT
ENDDO
!Wall Temperature
Twall = CalWallTemp(X, Tcool, htcoef, ThVar, ThOpt)
Tgap = CalGapTemp(X, ThVar, ThOpt)
Tfuel(1:npr) = X(1:npr)
Tfuel(npr1) = Tgap(1)
Tfuel(npr2:npr3) = x(npr1:npr2)
Tfuel(npr4) = Twall
Tfuel(npr5) = CalAvgFuelTemp(X, ThVar, ThOpt)
tfvol(1:npr) = X(1:npr) + CKELVIN
tfvol(npr1) = 0.5 * (Tgap(1) + Tgap(2)) + CKELVIN
tfvol(npr2) = 0.5 * (x(npr1) + x(npr2)) + CKELVIN

CONTINUE
NULLIFY(R)
END SUBROUTINE 

FUNCTION CalAvgFuelTemp(X, ThVar, ThOpt)
USE PARAM
USE TYPEDEF,            ONLY : ThVar_Type,      ThOpt_Type
IMPLICIT NONE
REAL :: X(100)
TYPE(ThVar_Type) :: ThVar
TYPE(ThOpt_Type) :: ThOpt
REAL :: CalAvgFuelTemp

INTEGER :: i, npr
REAL, POINTER :: R(:)

R => ThVar%R; npr = ThVar%npr
CalAvgFuelTemp = 0
DO i = 1, npr
  CalAvgFuelTemp = CalAvgFuelTemp + (R(i+1)**2 - R(i)**2) * X(i)
ENDDO
CalAvgFuelTemp = CalAvgFuelTemp / (R(npr+1)**2)


NULLIFY(R)
END FUNCTION

FUNCTION CalGapTemp(X, ThVar, ThOpt)
USE PARAM
USE TYPEDEF,           ONLY : ThVar_Type, ThOpt_Type
USE FuelProperty_Mod,  ONLY : CondFuel
IMPLICIT NONE
REAL :: x(100)
TYPE(ThVar_Type) :: ThVar
TYPE(ThOpt_Type) :: ThOpt
REAL :: CalGapTemp(2)

REAL :: rs, rw, tw, rgap, hg, hgt
REAL :: DelL, DelR, kfl, kfr
REAL :: betal, betar
REAL :: TL, TR
INTEGER :: npr, npr1, npr2
INTEGER :: i

npr = ThVar%npr; npr1 = ThVar%npr1; npr2 = ThVar%npr2

rs = ThVar%rs; rw = ThVar%rw
tw = ThVar%tw; rgap = ThVar%rgap; 
hg = ThOpt%hgap; hgt = hg * rs / rgap 

i=npr
DelL = ThVar%R(i+1) - ThVar%R(i); DelR = ThVar%R(i+3) - ThVar%R(i+2)
kfl = Condfuel(X(npr)+CKELVIN); kfr = Condfuel(x(npr1)+CKELVIN)
betaL = kfl / DelL; betaR = kfR / DelR

TL = X(npr); TR = X(npr1)
CalGapTemp(1) = hgt * TL * betaL + hg * TR * betaR + 2 * TL * betaL * BetaR
CalGapTemp(1) = CalGapTemp(1) / (hgt * betaL + hg * betaR + 2 * betaL *BetaR)

CalGapTemp(2) =  hgt * TL * BetaL +  TR * hg * BetaR + 2 * TR * BetaL * BetaR
CalGapTemp(2) = CalGapTemp(2) / (hgt * betaL + hg * betaR + 2 * betaL *BetaR)
END FUNCTION

FUNCTION CalWallTemp(X, Tcool, htcoef, ThVar, ThOpt)
USE PARAM
USE TYPEDEF,           ONLY : ThVar_Type, ThOpt_Type
USE FuelProperty_Mod,  ONLY : CondClad
IMPLICIT NONE
REAL :: CalWallTemp
REAL :: x(100)
REAL :: tcool, htcoef
TYPE(ThVar_Type) :: ThVar
TYPE(ThOpt_Type) :: ThOpt

REAL :: kf, Del, Beta
INTEGER :: npr, npr2, npr3, npr4

npr = ThVar%npr; npr2 = ThVar%npr2
npr3 = ThVar%npr3; npr4 = ThVar%npr4


!kf = CondClad(x(npr)+CKELVIN)
kf = CondClad(x(npr+2)+CKELVIN)
Del = ThVar%R(npr4) - ThVar%R(npr3)
beta = kf / Del
CalWallTemp = (htcoef * Tcool + 2 * beta * x(npr2)) / (htcoef + 2 * beta)
END FUNCTION

SUBROUTINE SetFuelCondLs(Diag, L, U, b, x, tcool, htcoef, qf, qshape, nr, ThVar, ThOpt, lMox)
USE PARAM
USE TYPEDEF,           ONLY : ThVar_Type,     ThOpt_Type
USE SteamTBL_mod,      ONLY : steamtbl
USE FuelProperty_Mod,  ONLY : CondFuel,       CondClad,       ftfavg 
USE LU_Mod,            ONLY : DirLu1Dsolver
IMPLICIT NONE
REAL :: diag(100), L(100), U(100), x(100), b(100), qshape(nr)
REAL :: tcool, htcoef, qf
INTEGER :: nr
TYPE(ThVar_Type) :: ThVar
TYPE(ThOPT_Type) :: ThOPT
LOGICAL :: lMox

REAL :: kf(100)
REAL, POINTER :: R(:)
REAL :: delL, delR, betar, betal, ktil
REAL :: hg, hgt
REAL :: rs, rw, tw, rgap
INTEGER :: npr, npr1, npr2, npr3, npr4, npr5
INTEGER :: i, m

npr = ThVar%npr; npr1 =ThVar%npr1
npr2 = ThVar%npr2; npr3 = Thvar%npr3
npr4= ThVar%npr4; npr5 = ThVar%npr5
rs = ThVar%rs; rw = ThVar%rw
tw = ThVar%tw; rgap = ThVar%rgap; 
R => ThVar%R

hg = ThOpt%hgap; hgt = hg * rs / rgap 

DO i = 1, npr
  kf(i) =   CondFuel(x(i))
ENDDO

DO i = npr1, npr2
  kf(i) =   CondClad(x(i))
ENDDO

i = 1
DelL = R(2) - R(1); DelR = R(3) - R(2)
betaL = kf(1) / DelL; betaR = kf(2) / DelR
ktil = 2._8 * betaL * betaR / (betaL + BetaR)
Diag(1) = ktil * R(2); U(1) = - ktil * R(2); L(1) = 0
b(1) = qf *(R(2)**2-R(1)**2) * qshape(1);

DO i = 2, npr-1
  DelL = R(i) - R(i-1); DelR = R(i+1) - R(i)
  betaL = kf(i-1) / DelL; betaR = kf(i) / DelR
  ktil = 2._8 * betaL * betaR / (betaL + BetaR)   
  DelL = R(i+1) - R(i); DelR = R(i+2) - R(i+1)
  Diag(i) = ktil * R(i); L(i) = -ktil * R(i)
  
  betaL = kf(i) / DelL; betaR = kf(i+1) / DelR
  ktil = 2._8 * betaL * betaR / (betaL + BetaR)      
  Diag(i) = Diag(i) + ktil * R(i+1); U(i) = - ktil * R(i+1)
  b(i) = qf *(R(i+1)**2-R(i)**2) * qshape(i);
ENDDO

i = npr
DelL = R(i) - R(i-1); DelR = R(i+1) - R(i)
betaL = kf(i-1) / DelL; betaR = kf(i) / DelR
ktil = 2._8 * betaL * betaR / (betaL + BetaR)   
Diag(i) = ktil * R(i); L(i) = -ktil * R(i)

DelL = R(i+1) - R(i); DelR = R(i+3) - R(i+2)
betaL = kf(i) / DelL; betaR = kf(i+1) / DelR
ktil = 2 * hg * betaL * BetaR / (hgt * betaL + hg * BetaR + 2 * betaL * BetaR)
Diag(i) = Diag(i) + ktil * R(i+1); U(i) = - ktil * R(i+1)
B(i) = qf *(R(i+1)**2-R(i)**2) * qshape(i);

B(1:i) = 0.5 * B(1:i) 

ktil = 2 * hgt * betaL * BetaR / (hgt * betaL + hg * BetaR + 2 * betaL * BetaR)
Diag(i+1) = ktil * R(i+2); L(i+1) = -ktil * R(i+2)

DelL = R(i+3) - R(i+2); betaL = kf(i+1) / DelL
DelR = R(i+4) - R(i+3); betaR = kf(i+2) / DelR
ktil = 2._8 * betaL * betaR / (betaL + BetaR) 
Diag(i+1) = Diag(i+1) + ktil * R(i+3); U(i+1) = - ktil * R(i+3)
B(i+1) = 0

Diag(i+2) = ktil * R(i+3); L(i+2) = -ktil * R(i+3)
BetaL = BetaR; ktil = 2*betaL * htcoef / (2*betaL + htcoef)
Diag(i+2) = Diag(i+2) + ktil * R(i+4) ; U(i+2) = 0
B(i+2) = ktil * R(i+4) * Tcool

NULLIFY(R)
END SUBROUTINE 

