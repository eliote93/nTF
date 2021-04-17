#include <defines.h>
SUBROUTINE InitTH()
USE PARAM
USE TYPEDEF
USE Geom,         ONLY : Core,           CellInfo,        CellPitch,     &
                         AsyPitch
USE PE_MOD,       ONLY : PE                         
USE Core_MOD,     ONLY : THInfo
USE TH_MOD,       ONLY : ThOpt,          ThVar,           hGapArray
USE CNTL,         ONLY : nTracerCntl
USE IOUTIL,       ONLY : OpenFile
USE FILES,        ONLY : io4,            filename,        SteamTblIdx
USE Th_MOD,       ONLY : SetInletTHInfo, InitCoolantVar,  AllocCoolantTH,   &
                         InitFuelThVar,  AllocFuelTh
USE SteamTBL_mod, ONLY : ReadSteamTbl
USE ALLOCS
use SubChCoupling_mod,      only: is_coupled
#ifdef MPI_ENV
USE MPIComm_Mod,  ONLY : MPIWaitTurn
#endif

USE HexData, ONLY : ncTyp

IMPLICIT NONE

REAL :: rs, rw, tw, rgt, TinK
INTEGER :: icel, i, ich, nEd
INTEGER :: nxy, nz, nrpallet

IF(nTracerCntl%lHex) THEN
  nEd = ncTyp
ELSE
  nEd = Core%nCellType
ENDIF

DO icel = 1, nEd
 CALL SetThCell(CellInfo(icel), ThOpt, nTracerCntl)
ENDDO

nxy = Core%nxy; nz = Core%nz
nrpallet = ThOpt%nrpellet

CALL Dmalloc0(hGapArray, 1, nz, 1, nxy)

!Allocation
CALL Dmalloc0(ThInfo%Tdop, 0, nz, 0, nxy + 1)
CALL Dmalloc0(ThInfo%Tcool, 1, nz, 0, nxy + 1)
CALL Dmalloc0(ThInfo%TcoolInOut, 1, 2, 1, nz+1, 0, nxy + 1)
CALL Dmalloc0(ThInfo%DenCool, 1, nz, 0, nxy + 1)
CALL Dmalloc0(ThInfo%CbmCool, 1, nz, 0, nxy + 1)
CALL Dmalloc(ThInfo%Tfvol, nrpallet + 2, nz, nxy)
CALL Dmalloc0(ThInfo%RelPower, 0, nz, 0, nxy)
CALL Dmalloc(ThInfo%qvol, nz, nxy)
!Get Steam Tables
#ifdef MPI_ENV
CALL MPIWaitTurn(PE%MPI_CMFD_COMM, PE%myCmfdRank, PE%nCmfdProc, .FALSE.)
#endif
CALL OPENFILE(io4, TRUE, TRUE, FALSE, filename(SteamTblIdx))
CALL ReadSteamTbl(io4)
CLOSE(io4)
#ifdef MPI_ENV
CALL MPIWaitTurn(PE%MPI_CMFD_COMM, PE%myCmfdRank, PE%nCmfdProc, .TRUE.)
#endif

!Assign ThVar variables

!Radial Nodalization
ThVar%npr = nrpallet; ThVar%npr1 = nrpallet + 1
ThVar%npr2 = nrpallet + 2; ThVar%npr3 = nrpallet + 3
ThVar%npr4 = nrpallet + 4; ThVar%npr5 = nrpallet + 5
ThVar%rs =ThOpt%PinDim(1) * epsm3
ThVar%rw =ThOpt%PinDim(2) * epsm3
ThVar%tw =ThOpt%PinDim(3) * epsm3
ThVar%rgap = ThVar%rw - ThVar%tw
ThVar%rgto =ThOpt%PinDim(4) * epsm3
if(is_coupled) ThVar%rgti =ThOpt%PinDim(5) * epsm3
ThVar%FuelDelR = ThVar%rs/nrpallet
ThVar%CLDDelR = ThVar%tw/2._8
rw = ThVar%rw; rgt = ThVar%rgto
CALL Dmalloc(THvar%r, nrpallet+5)
IF(.NOT. nTracerCntl%lSimpleTh)  THEN
  DO i = 1, nrpallet + 1
    ThVar%r(i) = ThVar%FuelDelR * (i-1)
  ENDDO
ELSE
  DO i = 1, nrpallet + 1
    ThVar%r(i) = ThVar%rs * SQRT(1._8*(i-1)/nrpallet)
  ENDDO
ENDIF
ThVar%r(nrpallet+2:nrpallet+4) = (/ThVar%rgap, ThVar%rgap+ThVar%CLDDelR, ThVar%rw/)

ThVar%afp = PI * ThVar%rs * ThVar%rs                !Fuel Pallet Area
Thvar%zeta = 2. * PI * rw

IF (nTracerCntl%lHex) THEN
  CALL HexSetThGeo(ThVar%nAsyCh, rw, rgt, ThVar%acf, ThVar%xi)
ELSE
  ThVar%acf = (ThVar%AsyPitch**2 - PI*(ThVar%nAsyCh * rw**2 + ThVar%nAsyGT * rgt**2))/ThVar%nAsyCh  !coolant Flow Area
  ThVar%xi = 2 * PI * (ThVar%nAsyCh * rw + ThVar%nAsyGT * rgt) / ThVar%nAsyCh
END IF

Thvar%zetap = Thvar%zeta / ThVar%acf
Thvar%Deq = 4. * ThVar%acf / ThVar%xi                
!Heat Deposition Rate
!ThVar%FracDc = 0.
ThVar%Fracdf = 1._8 - ThVar%FracDc 

IF(nTracerCntl%lthchconf) THEN 
  DO ich = 1, ThVar%nChType
    ThVar%ThCh(ich)%acf = (ThVar%ThCh(ich)%AsyPitch**2 - PI*(ThVar%ThCh(ich)%nAsyCh * rw**2 + ThVar%ThCh(ich)%nAsyGT * rgt**2))/ThVar%ThCh(ich)%nAsyCh  
    ThVar%ThCh(ich)%xi = 2 * PI * (ThVar%ThCh(ich)%nAsyCh * rw + ThVar%ThCh(ich)%nAsyGT * rgt) / ThVar%ThCh(ich)%nAsyCh
    Thvar%ThCh(ich)%zetap = Thvar%zeta / ThVar%ThCh(ich)%acf
    Thvar%ThCh(ich)%Deq = 4. * ThVar%ThCh(ich)%acf / ThVar%ThCh(ich)%xi
  END DO 
END IF

!Active Height
IF(Thvar%lhact) THEN 
  Thvar%hact = Thvar%hact * epsm2
ELSE
  THvar%Hact = 0
  DO i = 1, Core%nz
    IF(Core%lFuelPlane(i)) THvar%Hact = THvar%Hact + Core%hz(i) * epsm2
  ENDDO
END IF
CALL Dmalloc0(THvar%hz, 0,THvar%nzTH+1)
THvar%hz(1:Thvar%nzTH) = Core%Hz(1:Thvar%nzTH) * epsm2

!InLet Information
CALL SetInLetTHInfo(ThInfo, THVar, Core%Pin, nTracerCntl, nxy, nz)

!Allocation & Initiation
TinK = ThInfo%Tin + CKELVIN
ALLOCATE(ThInfo%CoolantTh(nxy))
CALL AllocCoolantTH(Core%Pin, ThInfo, nrpallet, nxy ,nz)
IF(nTracerCntl%lthchconf) THEN
  CALL InitCoolantVar_thch(Core, ThVar, Tink, ThInfo%PExit, ThVar%MdotFa, ThInfo%CoolantTh, nrpallet, nxy ,nz)
ELSE
  CALL InitCoolantVar(Tink, ThInfo%PExit, ThVar%MdotFa, ThVar%acf, ThInfo%CoolantTh, nrpallet, nxy ,nz)
END IF

ALLOCATE(ThInfo%FuelTh(nxy)); 
CALL AllocFuelTh(Core%Pin, ThInfo, nrpallet, nxy ,nz)
CALL InitFuelThVar(ThInfo%Tin, ThInfo%FuelTH, nrpallet, nxy, nz)
!
!CAll SetThPinType(Core)
!
IF(nTracerCntl%ThCh_mod .GT. 0) CALL InitThChGrp(Core, nTracerCntl)

! Initial Geuss Power Shape
CALL InitTHPwShape(THInfo, THvar, Core%lfuelPlane, nxy, nz)
CALL SetBoilingTemp(Thvar%BoilingTemp, ThInfo%PEXIT, ThInfo%Tin)
IF(nTracerCntl%lProblem .EQ. lTransient) CALL AllocTransientTH(THInfo, ThVar, nxy, nz)
END SUBROUTINE

SUBROUTINE SetBoilingTemp(BoilingTemp, PEXIT, Tin)
USE PARAM
USE SteamTBL_mod, ONLY : steamtbl
IMPLICIT NONE
REAL :: BoilingTemp, PEXIT, Tin
REAL :: wt, wh, hin, wrho, wvin, wxin, wbetain, wkapain, wcpin

!wt = 270 + CKELVIN
wt = Tin + CKELVIN
CALL SteamTBL(TRUE, pexit, wt, wh, wrho, wvin, wxin, wbetain, wkapain, wcpin)

BoilingTemp = 0
wt = 10;
DO
 IF(abs(wt - BoilingTemp) .LT. epsm5) EXIT
 BoilingTemp = wt
 wh = wh * 1.05
 CALL SteamTbl(FALSE, pexit, wt, wh, wrho, wvin, wxin, wbetain, wkapain, wcpin)
 continue
ENDDO
BoilingTemp = BoilingTemp - CKELVIN
END SUBROUTINE


SUBROUTINE SetInLetTHInfo(ThInfo, THVar, Pin, nTracerCntl, nxy, nz)
USE PARAM
USE TYPEDEF,      ONLY : THInfo_Type,           Pin_Type,     ThVar_Type
USE CNTL,         ONLY : nTracerCntl_Type
USE SteamTBL_mod, ONLY : steamtbl
USE geom,         ONLY : Core
IMPLICIT NONE

TYPE(THInfo_Type) :: ThInfo
TYPE(THVar_Type) :: ThVar
TYPE(Pin_Type) :: Pin(nxy)
TYPE(nTracerCntl_Type) :: nTracerCntl
INTEGER :: nz, nxy

REAL :: Tink, pexit, hin, din, wvin, wxin, wbetain, wkapain, wcpin
REAL :: avgnch, avghnch
INTEGER :: i
INTEGER :: ixya, ncha, ichtyp, ityp

IF(nTracerCntl%lThChConf) THEN 
  ncha = 0
  avgnch = 0. ; avghnch = 0.
  DO ixya = 1, Core%nxya
    ityp = Core%CoreMap(ixya)
    IF(.NOT. Core%AsyInfo(ityp)%lfuel) CYCLE
    ichtyp = Core%ThChMap(ixya)
    ncha = ncha + 1
    IF(ichtyp .EQ. 0) THEN
      avgnch = avgnch + ThVar%nAsyCh
      avghnch = avghnch + ThVar%nAsyCh * ThVar%hact
    ELSE
      avgnch = avgnch + ThVar%ThCh(ichtyp)%nAsyCh
      avghnch = avghnch + ThVar%ThCh(ichtyp)%nAsyCh * ThVar%ThCh(ichtyp)%hact
    END IF
  END DO
  avgnch = avgnch / ncha
  avghnch = avghnch / ncha
ELSE
  avgnch = ThVar%nAsyCh
  avghnch = THvar%Hact * ThVar%nAsyCh
END IF

ThInfo%Tin = nTracerCntl%TempInlet
ThInfo%PowFa = nTracerCntl%PowerFA !* 1.E6_8
ThInfo%PExit = nTracerCntl%PExit !* 1.E6_8
ThInfo%PowLin = ThInfo%PowFa/avghnch
ThInfo%PowLv = nTracerCntl%PowerLevel
ThInfo%MdotFa = nTracerCntl%fMdotFA/avgnch

ThVar%Tin = ThInfo%Tin
ThVar%MdotFa = ThInfo%MdotFa
THVar%Pexit = ThInfo%PExit
Tink = ThInfo%Tin + CKELVIN; pexit = ThInfo%PExit
CALL SteamTBL(TRUE, pexit, tink, hin, din, wvin, wxin, wbetain, wkapain, wcpin)

ThInfo%DenIn = Din
ThInfo%TDopIn = SQRT(TinK)
DO i = 1, nxy
  ThInfo%Tcool(1:nz, i) = ThInfo%Tin
  ThInfo%DenCool(1:nz, i) = Din
  ThInfo%CbmCool(1:nz,i) =  nTracerCntl%BoronPPM
  IF(.not. Pin(i)%lfuel) CYCLE
  ThInfo%TfVol(:, 1:nz, i) = ThInfo%Tin
  ThInfo%DenCool(1:nz, i) = Din
ENDDO

END SUBROUTINE

SUBROUTINE InitCoolantVar(Tink, pexit, MdotFa, ACF, CoolantTh, nrpallet, nxy, nz)
USE PARAM
USE TYPEDEF,      ONLY : THInfo_TYPE,    CoolantTH_Type
USE SteamTBL_mod, ONLY : steamtbl

IMPLICIT NONE
REAL :: TinK, Pexit, MdotFa, ACF
TYPE(CoolantTh_Type) :: CoolantTH(nxy)

INTEGER :: nrpallet, nxy, nz, I
REAL :: hin, din, wvin, wxin, wbetain, wkapain, wcpin
REAL :: RhoIn, RhoUin, RhoHUin, uin

CALL SteamTBL(TRUE, pexit, tink, hin, din, wvin, wxin, wbetain, wkapain, wcpin)

Rhoin = Din; RhoUin = MdotFA / Acf
RhoHUin = RhoUin * hin
Uin = RhoUin/Rhoin

DO i = 1, nxy
  IF(.NOT. CoolantTH(i)%lFuel) CYCLE
  CoolantTH(i)%hcool(1:nz) = hin
  CoolantTH(i)%rhou(0:nz) = RhoUin
  CoolantTH(i)%RhoHU(0:nz) = RhoHUin
  CoolantTH(i)%u(0:nz) = Uin
  CoolantTH(i)%ud(0:nz) = Uin
  CoolantTH(i)%DenCool(1:nz) = Din
  !CoolantTH(i)%tfuel(:, :) = Tink
ENDDO
END SUBROUTINE

SUBROUTINE InitCoolantVar_ThCh(Core, ThVar, Tink, pexit, MdotFa, CoolantTh, nrpallet, nxy, nz)
USE PARAM
USE TYPEDEF,      ONLY : THInfo_TYPE,    ThVar_Type, CoolantTH_Type,    CoreInfo_Type,      Pin_Type
USE SteamTBL_mod, ONLY : steamtbl
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(ThVar_Type) :: ThVar
REAL :: TinK, Pexit, MdotFa
TYPE(CoolantTh_Type) :: CoolantTH(nxy)

TYPE(Pin_Type), POINTER :: Pin(:)
INTEGER :: nrpallet, nxy, nz, I
INTEGER :: ixya, ichtyp
REAL :: hin, din, wvin, wxin, wbetain, wkapain, wcpin, acf_i
REAL :: RhoIn, RhoUin, RhoHUin, uin

CALL SteamTBL(TRUE, pexit, tink, hin, din, wvin, wxin, wbetain, wkapain, wcpin)
Rhoin = Din; 

Pin => Core%Pin

DO i = 1, nxy
  IF(.NOT. CoolantTH(i)%lFuel) CYCLE
  ixya = Pin(i)%iasy
  ichtyp = Core%THChMap(ixya)
  IF(ichtyp .EQ. 0) THEN 
    acf_i = ThVar%acf
  ELSE
    acf_i = ThVar%ThCH(ichtyp)%acf
  END IF

  RhoUin = MdotFA / Acf_i
  RhoHUin = RhoUin * hin
  Uin = RhoUin/Rhoin

  CoolantTH(i)%hcool(1:nz) = hin
  CoolantTH(i)%rhou(0:nz) = RhoUin
  CoolantTH(i)%RhoHU(0:nz) = RhoHUin
  CoolantTH(i)%u(0:nz) = Uin
  CoolantTH(i)%ud(0:nz) = Uin
  CoolantTH(i)%DenCool(1:nz) = Din
  !CoolantTH(i)%tfuel(:, :) = Tink
ENDDO
END SUBROUTINE

SUBROUTINE AllocCoolantTH(Pin, ThInfo, nrpallet, nxy, nz)
USE PARAM
USE TYPEDEF,      ONLY : THInfo_Type, CoolantTh_Type, Pin_Type
USE ALLOCS
IMPLICIT NONE
TYPE(Pin_Type) :: Pin(nxy)
TYPE(THInfo_Type) :: THInfo
TYPE(CoolantTh_Type), POINTER :: CoolantTH(:)
INTEGER :: nrpallet, nxy, nz, I

CoolantTH => THinfo%CoolantTH
DO I = 1, nxy
  IF(.NOT. Pin(i)%lFuel) CYCLE
  CoolantTH(i)%lfuel = TRUE
  CALL Dmalloc(CoolantTH(i)%hcool, nz)
  CALL Dmalloc0(CoolantTH(i)%rhou, 0, nz)
  CALL Dmalloc0(CoolantTH(i)%rhohu, 0, nz)
  CALL Dmalloc0(CoolantTH(i)%u, 0, nz)
  CALL Dmalloc0(CoolantTH(i)%ud, 0, nz)
  CALL Dmalloc(CoolantTH(i)%qeff, nz)
  CoolantTH(i)%qvol => ThInfo%qvol(:, i) 
  CoolantTH(i)%tcool => ThInfo%tcool(:, i) 
  CoolantTH(i)%TCoolInOut => ThInfo%TcoolInOut(:, :, i)
  CoolantTH(i)%DenCool => ThInfo%DenCool(:, i)
ENDDO
NULLIFY(CoolantTH)
END SUBROUTINE

SUBROUTINE InitFuelThVar(Tin, FuelTh, nrpallet, nxy, nz)
USE PARAM
USE TYPEDEF,    ONLY : FuelTH_Type
IMPLICIT NONE
TYPE(FuelTH_Type) :: FuelTH(nxy)
INTEGER :: nrpallet, nxy, nz
REAL :: Tin
INTEGER :: i, iz

DO i = 1, nxy
  IF(.NOT. FuelTH(i)%lFuel) CYCLE
  DO iz = 1, nz
    FuelTH(i)%Tfuel(:, iz) = Tin
  ENDDO
ENDDO 
END SUBROUTINE

SUBROUTINE AllocFuelTh(Pin, THInfo, nrpallet, nxy, nz)
USE PARAM
USE TYPEDEF,      ONLY : ThInfo_Type, FuelTh_Type, Pin_Type
USE ALLOCS
IMPLICIT NONE
TYPE(Pin_Type) :: Pin(nxy)
TYPE(THInfo_Type) :: ThInfo
INTEGER :: nrpallet, nxy, nz, i

TYPE(FuelTh_Type), POINTER :: FuelTh(:)

FuelTH => ThInfo%FuelTH
DO i = 1, nxy
  IF(.NOT. Pin(i)%lFuel) CYCLE
  FuelTH(i)%lFuel = TRUE
  CALL Dmalloc(FuelTh(i)%tfuel, nrpallet + 5, nz)
  CALL Dmalloc(FuelTh(i)%hflux, nz)
  CALL Dmalloc(FuelTh(i)%teff, nz)
  CALL Dmalloc(FuelTh(i)%htcoef, nz)
  !CALL Dmalloc(FuelTh(i)%qvol, nz)
  FuelTH(i)%qvol => ThInfo%qvol(:, i)
  FuelTH(i)%tcool => ThInfo%tcool(:, i)
  FuelTH(i)%Tfvol => ThInfo%tfvol(:, :, i)
  
  CALL Dmalloc(FuelTH(i)%lMox, nz)
  FuelTH(i)%lMox(1:nz) = Pin(i)%lMox(1:nz)
ENDDO
NULLIFY(FuelTH)
END SUBROUTINE

SUBROUTINE AllocTransientTH(THInfo, ThVar, nxy, nz)
USE PARAM
USE TYPEDEF,             ONLY : ThInfo_Type,         ThVar_Type,       CoolantTH_Type,        &
                                FuelTh_Type
USE TH_mod,              ONLY : THOpt
USE ALLOCS
IMPLICIT NONE
TYPE(ThInfo_Type) :: ThInfo
TYPE(ThVar_Type) :: ThVar
INTEGER :: nxy, nz

TYPE(CoolantTH_Type), POINTER :: CoolantTH(:)
TYPE(FuelTH_Type), POINTER :: FuelTH(:)

INTEGER :: ixy, iz

CoolantTH => ThInfo%CoolantTH
FuelTH => ThInfo%FuelTH

DO ixy = 1, nxy
  CALL DMALLOC0(CoolantTH(ixy)%Tcoold, 1, nz)
  CALL DMALLOC0(CoolantTH(ixy)%DenCoold, 1, nz)
  CALL DMALLOC0(CoolantTH(ixy)%hCoold, 1, nz)
  CALL DMALLOC0(CoolantTH(ixy)%rhoud, 0, nz)
  CALL DMALLOC0(CoolantTH(ixy)%rhohud, 0, nz)
  CALL DMALLOC0(CoolantTH(ixy)%qeffd, 1, nz)
  CALL DMALLOC0(FuelTH(ixy)%qvold, 1, nz)
  CALL DMALLOC0(FuelTH(ixy)%Tfueld, 1, ThOpt%nrpellet + 5, 1, nz)
ENDDO

END SUBROUTINE



SUBROUTINE InitTHPwShape(THInfo, THvar, lfuelPlane, nxy, nz)
USE PARAM
USE TYPEDEF,      ONLY : THInfo_Type,   THVar_Type
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
TYPE(THInfo_Type) :: ThInfo 
TYPE(THvar_Type) :: ThVar
LOGICAL :: lFuelPlane(nz)
INTEGER :: nxy, nz

INTEGER :: iz
REAL :: hact, refsave, buckling, p0
REAL :: z, cosi, cosip1, PW
REAL, POINTER :: RelPower(:, :), HzTH(:)

RelPower => ThInfo%RelPower
HzTh => ThVar%Hz
hact = ThVar%hact*100.0_8
refSave = 10.0_8
buckling = PI / (2. * refsave + hact)
p0 = hact /(COS(buckling * refSave) - COS(buckling * (hact + refSave))) * 0.01_8

z = refSave; cosi = cos(buckling * refsave)
DO iz = 1, nz
  CALL CP_CA(RelPower(iz, 0:nxy), 0._8, nxy + 1)
  IF(.NOT. lFuelPlane(iz)) CYCLE
  z = z + hzTH(iz) * 100._8
  Cosip1 = COS(buckling * z)
  PW = p0 * (Cosi - Cosip1)/hzth(iz)
  IF(PW .LT. 0._8) PW = 0._8
  PW = 1._8
  CALL CP_CA(RelPower(iz, 1:nxy), PW, nxy)
  cosi = cosip1
ENDDO
!Buckling = PI / ()
NULLIFY(RelPower)
END SUBROUTINE
!
SUBROUTINE InitThChGrp(Core, nTracerCntl)
USE PARAM
USE TYPEDEF,   ONLY : CoreInfo_Type,          Asy_Type,        AsyInfo_Type
USE Cntl,      ONLY : nTracerCntl_Type
USE BasicOperation, ONLY : CP_CA, CP_VA
USE ALLOCS
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(nTracerCntl_Type) :: nTracerCntl

TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)

INTEGER :: ThCh_mod
INTEGER :: nAsyType
INTEGER :: nsubdiv, ncx, ncell, ngrp
INTEGER :: nx0, ny0, nx, ny, divtyp
INTEGER :: iasy, icase
INTEGER :: i, j, k, m, ix, iy, ix1, ix2, iy1, iy2
INTEGER :: BaseGrp(1000, 1000)
INTEGER :: Rg(0:1000,0:1)
ThCh_mod = nTracerCntl%ThCh_mod
nAsyType = Core%nAsyType
ncx = Core%nxc; nCell = ncx * ncx
Asy => Core%Asy
AsyInfo => Core%AsyInfo

CALL CP_CA(BaseGrp, 0, 1000, 1000)

divtyp = mod(ncx, ThCh_mod)
nsubdiv = ncx / ThCh_mod

Rg(1:ThCh_mod, 0:1) = nsubdiv
IF(divtyp .NE. 0) THEN
!  j = ThCh_mod/2
!  Rg(j, 0) = Rg(j, 0) + 1
!  j = ThCh_mod/2+1
!  Rg(j, 1) = Rg(j, 1) + 1
  DO j = 1, divtyp
    Rg(j, 0) = Rg(j, 0) + 1
    Rg(ThCh_mod + 1 - j, 1) = Rg(ThCh_mod + 1 - j, 1) + 1
  ENDDO
ENDIF
Rg(0, :) = 0
DO j = 2, ThCh_mod
  Rg(j, 0) = Rg(j, 0) + Rg(j-1, 0)
  Rg(j, 1) = Rg(j, 1) + Rg(j-1, 1)
ENDDO

nx = ncx; ny = ncx
m = 0
DO j = 1, ThCh_Mod
  iy1 = Rg(j-1, 1) + 1; iy2 = Rg(j, 1)
  DO i = 1, ThCh_Mod
    k = mod(j-1, 2); m = m + 1
    ix1 = Rg(i-1, k) + 1; ix2 = Rg(i, k)
    BaseGrp(ix1:ix2, iy1:iy2) = m
  ENDDO
ENDDO
ngrp = m

DO iasy = 1, nasytype
  AsyInfo(iasy)%nThChGrp = ngrp
  CALL Dmalloc(AsyInfo(iasy)%ThChGrp, AsyInfo(iasy)%nxy)
  ix1 = 1; ix2 = nx
  iy1 = 1; iy2 = ny
  IF(AsyInfo(iasy)%lCentX) THEN
    iy1 = iy2 - AsyInfo(iasy)%ny + 1
  ELSEIF(AsyInfo(iasy)%lCentY) THEN
    ix1 = ix2 - AsyInfo(iasy)%nx + 1
  ELSEIF(AsyInfo(iasy)%lCentXY) THEN
    ix1 = ix2 - AsyInfo(iasy)%nx + 1
    iy1 = iy2 - AsyInfo(iasy)%ny + 1
  ENDIF
  k = 0
  DO iy = iy1, iy2
    DO ix = ix1, ix2
      k = k + 1
      AsyInfo(iasy)%ThChGrp(k) = BaseGrp(ix, iy)
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE ChangeFlowRate(Core, ThInfo, Thvar, FlowRateLevel)
USE PARAM
USE Typedef,          ONLY : CoreInfo_Type,     ThInfo_Type,      Thvar_Type,     &
                             CoolantTH_Type
USE SteamTBL_mod, ONLY : steamtbl
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(ThInfo_Type) :: ThInfo
TYPE(ThVar_Type) :: ThVar
REAL :: FlowRateLevel

REAL :: ACF, RhoUin, RhoHUin, RhoIn, Uin, MdotFa
REAL :: PEXIT, TinK, Hin ,Din, wvin, wxin, wbetain, wkapin, wkapain, wcpin

INTEGER :: i

TinK = ThInfo%Tin + CKELVIN
PEXIT = ThInfo%PExit 

CALL SteamTBL(TRUE, pexit, tink, hin, din, wvin, wxin, wbetain, wkapain, wcpin)

ThVar%MdotFa = ThInfo%MdotFa * FlowRateLevel

Acf = ThVar%ACF
MdotFa = ThVar%MdotFa

RhoUin = MdotFA / Acf
RhoHUin = RhoUin * hin
Uin = RhoUin/Rhoin

DO i = 1, Core%nxy
  IF(.NOT. ThInfo%CoolantTH(i)%lFuel) CYCLE
  ThInfo%CoolantTH(i)%rhou(0) = RhoUin
  ThInfo%CoolantTH(i)%RhoHU(0) = RhoHUin
  ThInfo%CoolantTH(i)%u(0) = Uin
  ThInfo%CoolantTH(i)%ud(0) = Uin
ENDDO

END SUBROUTINE

SUBROUTINE ChangeFlowRate_ThCh(Core, ThInfo, Thvar, FlowRateLevel)
USE PARAM
USE Typedef,          ONLY : CoreInfo_Type,     ThInfo_Type,      Thvar_Type,     &
                             CoolantTH_Type,    Pin_Type
USE SteamTBL_mod, ONLY : steamtbl
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(ThInfo_Type) :: ThInfo
TYPE(ThVar_Type) :: ThVar
REAL :: FlowRateLevel

TYPE(Pin_Type), POINTER :: Pin(:)
REAL :: ACF, RhoUin, RhoHUin, RhoIn, Uin, MdotFa
REAL :: PEXIT, TinK, Hin ,Din, wvin, wxin, wbetain, wkapin, wkapain, wcpin

INTEGER :: i, ixya, ichtyp

TinK = ThInfo%Tin + CKELVIN
PEXIT = ThInfo%PExit 

CALL SteamTBL(TRUE, pexit, tink, hin, din, wvin, wxin, wbetain, wkapain, wcpin)

ThVar%MdotFa = ThInfo%MdotFa * FlowRateLevel
MdotFa = ThVar%MdotFa

Pin => Core%Pin
DO i = 1, Core%nxy
  IF(.NOT. ThInfo%CoolantTH(i)%lFuel) CYCLE
  ixya = Pin(i)%iasy
  ichtyp = Core%ThChMap(ixya)
  IF(ichtyp .EQ. 0) THEN 
    Acf = ThVar%ACF
  ELSE
    Acf = ThVar%ThCh(ichtyp)%acf
  END IF

  RhoUin = MdotFA / Acf
  RhoHUin = RhoUin * hin
  Uin = RhoUin/Rhoin

  ThInfo%CoolantTH(i)%rhou(0) = RhoUin
  ThInfo%CoolantTH(i)%RhoHU(0) = RhoHUin
  ThInfo%CoolantTH(i)%u(0) = Uin
  ThInfo%CoolantTH(i)%ud(0) = Uin
ENDDO

END SUBROUTINE