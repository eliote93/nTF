#include <defines.h>

SUBROUTINE SetLocalBurnup(Core, FXR, Power, normalizer, Tsec, lCorrectStep, PE)
USE PARAM
USE TYPEDEF,       ONLY : CoreInfo_Type,   FxrInfo_Type, PE_Type,      &
                         Pin_Type,        Cell_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
REAL, POINTER :: Power(:, :)
REAL :: Normalizer, Tsec
LOGICAL :: lCorrectStep
TYPE(PE_TYPE) :: PE

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(FxrInfo_Type), POINTER :: myFxr

REAL :: Area,norm_Mw, powersum, burnup, Tday
INTEGER :: FsrIdxSt, FxrIdxSt
INTEGER :: ifxr, ifsr, icel, iz, ipin, iLocalFsr, nFsrInFxr, nLocalFxr
INTEGER :: i, j

CellInfo => Core%CellInfo
Pin => Core%Pin
powersum = 0; Tday = Tsec / 86400._8
norm_Mw = normalizer / 1.e+6_8
DO iz = PE%myzb, PE%myze
  DO ipin = 1, Core%nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
    icel = Pin(ipin)%Cell(iz)
    nLocalFxr = CellInfo(icel)%nFxr
    DO j =1, nLocalFxr
      ifxr = FxrIdxSt + j - 1; nFsrInFxr = Fxr(ifxr, iz)%nFsrInFxr
      IF(.NOT. Fxr(ifxr, iz)%lfuel) CYCLE
      burnup =0
      DO i = 1, nFsrInFxr
        iLocalFsr = CellInfo(icel)%MapFxr2FsrIdx(i, j)
        iFsr = FsrIdxSt + iLocalFsr - 1
        Area = CellInfo(icel)%vol(iLocalFsr);
        !FxrAvgPhi = FxrAvgPhi + Area * Phis(ifsr, iz, 1:ng)
        burnup = burnup + area * core%hz(iz) * power(ifsr, iz)
        powersum  = powersum + area * core%hz(iz) * power(ifsr, iz) * normalizer
      ENDDO
      burnup = burnup * norm_Mw * Tday / Fxr(ifxr, iz)%Hmkg0
      IF(.NOT. lCorrectStep) THEN
        Fxr(ifxr, iz)%burnup = Fxr(ifxr, iz)%burnup_past + burnup
      ELSE
        burnup = 0.5_8*(burnup + Fxr(ifxr, iz)%burnup - Fxr(ifxr, iz)%burnup_past)
        Fxr(ifxr, iz)%burnup = Fxr(ifxr, iz)%burnup_past + burnup
        Fxr(ifxr, iz)%burnup_past = Fxr(ifxr, iz)%burnup
      ENDIF
    ENDDO
  ENDDO
ENDDO
NULLIFY(CellInfo, Pin)
END SUBROUTINE

FUNCTION GetCoreHmMass(Core, FXR, PE)
USE PARAM
USE TYPEDEF,       ONLY : CoreInfo_Type,   FxrInfo_TYPE,   PE_TYPE,    &
                          Pin_Type,        Cell_Type
USE nuclidmap_mod, ONLY : iposnucl,        AtomicWeight
#ifdef MPI_ENV
USE MPIComm_mod,   ONLY : REDUCE
#endif
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(PE_TYPE) :: PE

REAL :: GetCoreHmMass

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(FxrInfo_Type), POINTER :: myFxr

REAL :: area, HmReg, HmCell, HmPln, HmCore, buf

INTEGER :: iz, ixy, ifxr, icell, iso, id
INTEGER :: i, j, k
INTEGER :: nFSR, nFXR, nxy
INTEGER :: FxrIdxSt, nLocalFxr, niso
INTEGER :: myzb, myze

Pin => Core%Pin; Cell => Core%CellInfo
nFSR = Core%nCoreFsr; nFxr = Core%nCoreFxr
nxy = Core%nxy
myzb = PE%myzb; myze = PE%myze
GetCoreHmMass = 0
HmCore = 0
DO iz = myzb, myze
  HmPln = 0
  DO ixy = 1, nxy
    icell = Pin(ixy)%cell(iz)
    FxrIdxSt = Pin(ixy)%FxrIdxSt; nLocalFxr = Cell(icell)%nFxr
    HmCell = 0
    DO i = 1, nLocalFxr
      ifxr = FxrIdxSt + i - 1
      IF(.NOT. Fxr(ifxr, iz)%lfuel) CYCLE
      myFxr => Fxr(ifxr, iz); niso = myFxr%niso
      area = myFxr%area; HmReg = 0
      DO iso = 1, niso
        id = myFxr%idiso(iso)
        IF(id .LT. 90000) CYCLE
        HmReg = HmReg + myFxr%pnum(iso) * AtomicWeight(id) / AVOGADRO
      ENDDO  !Isotope Sweep
      HmCell = HmCell + HmReg * Area
      !myFxr%Hmkg0 = HmReg * Area * Core%Hz(iz) / 1000._8
    ENDDO  !Local Fxr Sweep
    HmPln = HmPln + HmCell
  ENDDO
  HmCore = HmCore + HmPln * Core%Hz(iz)
ENDDO

#ifdef MPI_ENV
CALL REDUCE(HmCore, buf, PE%MPI_RTMASTER_COMM, .TRUE.)
HmCore = buf
#endif

GetCoreHmMass = HmCore
GetCoreHmMass = GetCoreHmMass * epsm6    !Ton Unit Conversion

END FUNCTION

SUBROUTINE SetBurnUpStepInfo(PowerCore,DeplCntl)
USE PARAM
USE DeplType,   ONLY : DeplCntl_Type
IMPLICIT NONE
TYPE(DeplCntl_Type) :: DeplCntl
REAL :: PowerCore


REAL :: HmMass_kg0
REAL, POINTER :: T_efpd(:), T_mwdkghm(:)

INTEGER :: nBurnUpStep, BurnUpType
INTEGER :: i


HmMass_kg0 = DeplCntl%Hm_Mass0_kg
nBurnUpStep = DeplCntl%nBurnUpStep; BurnUpType = DeplCntl%BurnUpType

T_efpd => DeplCntl%T_efpd; T_mwdkghm => DeplCntl%T_mwdkghm
T_mwdkghm(0) = 0
T_efpd(0) = 0
IF(BurnUpType .EQ. 1) THEN
  DO i = 1, nBurnUpStep
    T_mwdkghm(i) = T_efpd(i) * PowerCore / HmMass_kg0
  ENDDO
ELSEIF(BurnUpType .EQ. 2) THEN
  DO i = 1, nBurnUpStep
    T_efpd(i) = T_mwdkghm(i)/PowerCore * hmMass_kg0
  ENDDO
ENDIF

DeplCntl%PowerCore = PowerCore

NULLIFY(T_efpd); NULLIFY(T_mwdkghm)

END SUBROUTINE

SUBROUTINE MakeDeplXs1g(Fxr, ifxr, ng, DeplXs, NormFactor, Core, ipin, iz, GroupInfo, PE, DeplCntl)
USE PARAM
USE TypeDef,          ONLY : CoreInfo_Type,     FxrInfo_Type,      Cell_Type,       PE_TYPE,    &
                             GroupInfo_Type,    XsMac_Type
USE DeplType,         ONLY : DeplXs_Type, DeplCntl_Type
USE Depl_Mod,         ONLY : XsMac
USE MacXsLib_Mod,     ONLY : MacXsBase,  EffMacXs,  MacXS1gBase
USE BasicOperation,   ONLY : CP_VA
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type),POINTER :: Fxr(:,:)
TYPE(FxrInfo_Type),POINTER :: myFxr
TYPE(DeplXs_Type) :: DeplXS
TYPE(PE_TYPE) :: PE
TYPE(DeplCntl_Type) :: DeplCntl
TYPE(GroupInfo_Type) :: GroupInfo
REAL :: NormFactor
INTEGER :: ng, ipin, iz, ifxr

REAL, POINTER :: phis(:)
REAL, POINTER :: xsa(:), xsf(:), xsn2n(:), xsn3n(:)
REAL, POINTER :: IsoXsMacf(:, :), IsoXsMacA(:, :)
REAL :: Area, Temp
REAL :: Phi1g
INTEGER, POINTER :: idiso(:)
REAL, POINTER :: pnum(:)
INTEGER :: nFsrInFxr, niso, iResoGrp1, iResoGrp2
INTEGER :: ifsr, ig
INTEGER:: i
LOGICAL :: FlagF

myFxr => Fxr(ifxr,iz)

!Resonance Range
iresoGrp1 = GroupInfo%nofg + 1; iresoGrp2 = GroupInfo%nofg + GroupInfo%norg
#define Depl_mod1
#ifdef Depl_mod0
Temp = myFxr%Temp; niso = myFxr%niso_past
idiso => myFxr%idiso_past; pnum => myFxr%pnum_past
#endif
#ifdef Depl_mod1
Temp = myFxr%Temp; niso = myFxr%niso
idiso => myFxr%idiso; pnum => myFxr%pnum
#endif

phis => DeplXs%AvgPhi
!TempSubgrp = Fxr%SubGrpTemp
!1G XS
phi1g = SUM(Phis(1:ng))

xsa => DeplXS%xsa; xsf => DeplXs%xsf
xsn2n => DeplXs%xsn2n; xsn3n => DeplXs%xsn3n

xsa = 0; xsf = 0;
DO i = 1, niso
  CALL MacXS1gBase(xsa(i:i), xsf(i:i), xsn2n(i:i), xsn3n(i:i), Temp, 1, idiso(i:i), pnum(i:i), phis(1:ng), ng)
ENDDO

CALL MacXsBase(XsMac(DeplXs%tid), temp, niso, idiso, pnum, 1, ng, ng, 1._8, FALSE, TRUE)
IsoXsMacf => XsMac(DeplXs%tid)%IsoXsMacf
IsoXsMacA => XsMac(DeplXs%tid)%IsoXsMacA
IF(myFxr%lRes) THEN

  DO ig = iResoGrp1, iResoGrp2

        DO i = 1, niso
            IsoXsMacF(i,ig) = IsoXsMacF(i,ig) * myFXR%fresoFIso(i,ig)
            IsoXsMacA(i,ig) = IsoXsMacA(i,ig) * myFXR%fresoAIso(i,ig)
        ENDDO
  ENDDO
ENDIF
DO i = 1, niso
  IF(pnum(i) .GT. epsm20) THEN
    xsa(i) = 0; xsf(i) = 0
    DO ig = 1, ng
      xsa(i) = xsa(i) + Phis(ig) * IsoXsMacA(i, ig)
      xsf(i) = xsf(i) + Phis(ig) * IsoXsMacF(i, ig)
    ENDDO
    xsa(i) = xsa(i) / phi1g; xsf(i) = xsf(i) / phi1g
  ENDIF
  xsa(i) = xsa(i) / pnum(i); xsf(i) = xsf(i) / pnum(i)
  xsn2n(i) = xsn2n(i) / pnum(i); xsn3n(i) = xsn3n(i) / pnum(i)
  xsa(i)  = xsa(i)  + xsn3n(i)
ENDDO
DeplXS%phi1g = phi1g * NormFactor
NULLIFY(idiso, pnum)
NULLIFY(xsa, xsf, xsn2n, xsn3n)
END SUBROUTINE

SUBROUTINE SaveFxrIsoInfo(Fxr, nFxr)
USE PARAM
USE TYPEDEF,     ONLY : FxrInfo_Type
IMPLICIT NONE
TYPE(FxrInfo_Type) :: Fxr(1:nFxr)
INTEGER :: nFxr
INTEGER :: ifxr, niso

DO ifxr = 1, nFxr
  IF(.NOT. Fxr(ifxr)%ldepl) CYCLE
  niso = Fxr(ifxr)%niso_depl
  Fxr(ifxr)%niso_past = niso
  Fxr(ifxr)%idiso_past(1:niso) = Fxr(ifxr)%idiso(1:niso)
  Fxr(ifxr)%pnum_past(1:niso) = Fxr(ifxr)%pnum(1:niso)
  IF(Fxr(ifxr)%l_pnum_all) Fxr(ifxr)%pnum_past_all = Fxr(ifxr)%pnum_all  ! 16/02/11 Depletion timestep bug fixed
ENDDO

END SUBROUTINE



FUNCTION FluxNormalizeFactor(Core, FmInfo, GroupInfo, PowerCore, lCritSpec, lXsLib, PE)
USE PARAM
USE TYPEDEF,    ONLY : CoreInfo_Type,        FmInfo_Type,         GroupInfo_TYPE, &
                       PE_TYPE,                                                   &
                       FxrInfo_Type,         Cell_Type,           Pin_Type,       &
                       XsMac_Type
USE BenchXs,       ONLY : xskfBen,           xskfDynBen
USE MacXsLib_Mod, ONLY : MacXsKf
USE BasicOperation, ONLY : CP_CA, MULTI_VA
#ifdef MPI_ENV
USE MPIComm_mod, ONLY : REDUCE
#endif
USE TRAN_MOD,     ONLY : TranInfo,            TranCntl
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_TYPE) :: PE
LOGICAL :: lXsLib, lCritSpec
REAL :: FluxNormalizeFactor, PowerCore

TYPE(Fxrinfo_type),POINTER :: Fxr(:, :)

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Fxrinfo_type),POINTER :: myFxr
TYPE(XsMac_Type), SAVE :: XsMac

REAL, POINTER :: phis(:, :, :), hz(:), SpecConv(:)
REAL, POINTER :: xsmackf(:)

REAL :: pwsum, localpwsum, area, FxrPhis(1000)
REAL :: Buf

INTEGER :: nxy, nCoreFsr, nCoreFxr, myzb, myze
INTEGER :: FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, ng
INTEGER :: ipin, icel, ifsrlocal, ifsr, ifxr, iz, itype, ig
INTEGER :: iResoGrpBeg, iResoGrpEnd, norg
INTEGER :: i, j, k



!lXsLib = nTracerCntl%lXsLib !; lCritSpec = nTracerCntl%lCritSpec

Pin => Core%Pin; hz => Core%hz
CellInfo => Core%CellInfo; nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr; nxy = Core%nxy; ng = GroupInfo%ng
phis => FmInfo%Phis; Fxr => FmInfo%Fxr
IF(lCritSpec) SpecConv => FmInfo%SpecConv

myzb = PE%myzb; myze = PE%myze


IF(lxslib) THEN
  iResoGrpBeg = GroupInfo%nofg + 1
  iResoGrpEnd = GroupInfo%nofg + GroupInfo%norg
  norg = GroupInfo%norg
ENDIF
IF(.NOT. lxsLib) ALLOCATE(xsmackf(ng))


pwsum = 0
DO iz = myzb, myze
  DO ipin = 1, nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
    icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
      myFxr => Fxr(ifxr, iz)
      IF(lXsLib) Then
        CALL MacXsKf(XsMac, myFxr, 1, ng, ng, 1._8, FALSE)
        xsmackf => XsMac%XsMackf
        IF(myFxr%lres) THEN
          do ig = iResoGrpBeg,iResoGrpEnd
            XsMackf(ig) = XsMackf(ig) * myFxr%fresokf(ig)
          enddo
        ENDIF
      ELSE
        ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1,j)
        !itype = CellInfo(icel)%iReg(ifsrlocal)
        itype = myFxr%imix
        IF(TranCntl%lDynamicBen) THEN
          CALL xskfDynben(itype, TranInfo%fuelTemp(ipin, iz), 1, ng, xsmackf)
        ELSE
          CALL xskfben(itype, 1, ng, xsmackf)
        END IF
      ENDIF

      !Get Fxr Flux
      FxrPhis(1:ng) = 0
      DO i = 1, nFsrInFxr
        ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1  !Global FSR Index
        ifsrlocal = Cellinfo(icel)%MapFxr2FsrIdx(i, j)
        area = CellInfo(icel)%vol(ifsrlocal)
        DO ig = 1, ng
          FxrPhis(ig) = FxrPhis(ig) + area * phis(ifsr, iz, ig)
        ENDDO
      ENDDO
      IF(lCritSpec) CALL MULTI_VA(SpecConv(1:ng), FxrPhis(1:ng), ng)
      localpwsum = 0
      DO ig = 1, ng
        localpwsum = localpwsum + FxrPhis(ig) * xsmackf(ig)
      ENDDO
      pwsum = pwsum + localPwSum * hz(iz)
    ENDDO
  ENDDO
ENDDO

#ifdef MPI_ENV
CALL REDUCE(pwsum, buf, PE%MPI_RTMASTER_COMM, .TRUE.)
pwsum = buf
#endif

FluxNormalizeFactor = PowerCore / PwSum * 1.E6_8

NULLIFY(Pin); NULLIFY(Hz);
NULLIFY(CellInfo); NULLIFY(phis)
NULLIFY(Fxr)

IF(.NOT. lxsLib) Deallocate(xsmackf)
IF(lXsLib) NULLIFY(XsMackf)
IF(lCritSpec) NULLIFY(SpecConv)

END FUNCTION

!SUBROUTINE ConstDeplVas(DeplVars, DeplXs, IdIso, pnum, niso, nIsoLib, nIsoDepl)
!USE PARAM
!USE DeplType,    ONLY : DeplXs_Type,      DeplVars_Type
!USE nuclidmap_mod,ONLY : iposiso
!USE BasicOperation, ONLY : CP_CA
!IMPLICIT NONE
!TYPE(DeplVars_Type) ::DeplVars
!TYPE(DeplXs_Type) ::DeplXs
!INTEGER :: IdIso(nIsoLib)
!REAL :: pnum(nIsoLib)
!INTEGER :: niso, nIsoLib, nIsoDepl
!INTEGER :: niso, nIsoLib, nIsoDepl
!
!INTEGER, POINTER :: MapXs2Dep(:)
!REAL, POINTER :: BurnUpXs(:, :), IsoNum(:)
!INTEGER :: i, j, id_lib, id_depl
!
!MapXs2Dep => DeplVars%MapXs2Dep
!BurnUpXs => DeplVars%BurnUpXs
!IsoNum => DeplVars%IsoNum
!CALL CP_CA(IsoNum, 0._8, nIsoDepl)
!CALL CP_CA(BurnUpXS, 0._8, 4, nIsoDepl)
!DO i= 1, nIso
!  j = IdIso(i); id_lib = iPosIso(j)
!  id_depl = MapXs2Dep(id_lib)
!  IsoNum(id_depl) = pnum(i)
!  BurnUpXs(1, id_depl) = DeplXs%xsa(i)
!  BurnUpXs(2, id_depl) = DeplXs%xsf(i)
!  BurnUpXs(3, id_depl) = DeplXs%xsn2n(i)
!  BurnUpXs(4, id_depl) = DeplXs%xsn3n(i)
!ENDDO
!DeplVars%phi1g = DeplXs%phi1g
!END SUBROUTINE

SUBROUTINE ConstDeplVas_QD(DeplVars, DeplXs, Fxr, nIsoLib, nIsoDepl, DeplCntl)
USE PARAM
USE TypeDef,     ONLY : FxrInfo_Type
USE DeplType,    ONLY : DeplCntl_Type, DeplXs_Type,      DeplVars_Type
USE nuclidmap_mod,ONLY : iposiso
USE BasicOperation, ONLY : CP_CA, CP_VA
IMPLICIT NONE
TYPE(DeplCntl_TYPE) :: DeplCntl
TYPE(DeplVars_Type) ::DeplVars
TYPE(DeplXs_Type) ::DeplXs
TYPE(FxrInfo_Type) :: Fxr
!INTEGER :: IdIso(nIsoLib)
!REAL :: pnum(nIsoLib)
INTEGER :: nIsoLib, nIsoDepl

real :: w, tp, t, w0(3), w_QI(3), w_LI(3), w_LE(3)

INTEGER, POINTER :: IdIso(:)
REAL, POINTER :: pnum(:)
INTEGER, POINTER :: MapXs2Dep(:)
REAL, POINTER :: BurnUpXs(:, :), IsoNum(:)
INTEGER :: i, j, id_lib, id_depl, niso
LOGICAL :: lCorrectStep

lCorrectStep = .NOT. DeplCntl%lPredict

MapXs2Dep => DeplVars%MapXs2Dep
BurnUpXs => DeplVars%BurnUpXs
IsoNum => DeplVars%IsoNum
CALL CP_CA(IsoNum, 0._8, nIsoDepl)
CALL CP_CA(BurnUpXS, 0._8, 4, nIsoDepl)
IdIso => Fxr%IdIso_past; pnum => Fxr%pnum_past
nIso = Fxr%nIso_past
DO i= 1, nIso
  j = IdIso(i); id_lib = iPosIso(j)
  id_depl = MapXs2Dep(id_lib)
  IsoNum(id_depl) = pnum(i)
ENDDO

t = DeplCntl%tsec
tp = DeplCntl%ptsec
w_LI = (/0._8, 0.5_8, 0.5_8/)
w_LE = (/-DeplCntl%tsec/DeplCntl%ptsec/2._8, DeplCntl%tsec/DeplCntl%ptsec/2._8 + 1, 0._8/)
W_QI = (/-t*t/(6*tp*(t+tp)), 0.5 + t / (6 * tp),  0.5-t/(6*(t+tp))/)

IF(DeplCntl%NowStep == 1) THEN
  IF(.NOT. lCorrectStep) THEN
    IdIso => Fxr%IdIso; nIso = Fxr%nIso_depl
    DO i= 1, nIso
      j = IdIso(i); id_lib = iPosIso(j)
      id_depl = MapXs2Dep(id_lib)
      BurnUpXs(1, id_depl) = DeplXs%xsa(i)
      BurnUpXs(2, id_depl) = DeplXs%xsf(i)
      BurnUpXs(3, id_depl) = DeplXs%xsn2n(i)
      BurnUpXs(4, id_depl) = DeplXs%xsn3n(i)
    ENDDO
    DeplVars%phi1g = DeplXs%phi1g
  ELSE
    IdIso => Fxr%DeplXs1g(0)%IdIso; nIso = Fxr%DeplXs1g(0)%niso
    DO i = 1, niso
      j = IdIso(i); id_lib = iPosIso(j)
      id_depl = MapXs2Dep(id_lib)
      BurnUpXs(1, id_depl) = 0.5_8 * Fxr%DeplXs1g(0)%xsa(i)
      BurnUpXs(2, id_depl) = 0.5_8 * Fxr%DeplXs1g(0)%xsf(i)
      BurnUpXs(3, id_depl) = 0.5_8 * Fxr%DeplXs1g(0)%xsn2n(i)
      BurnUpXs(4, id_depl) = 0.5_8 * Fxr%DeplXs1g(0)%xsn3n(i)
    ENDDO
    DeplVars%Phi1g = 0.5 * Fxr%DeplXs1g(0)%phi1g
    IdIso => Fxr%IdIso; nIso = Fxr%nIso_depl
    DO i = 1, niso
      j = IdIso(i); id_lib = iPosIso(j)
      id_depl = MapXs2Dep(id_lib)
      BurnUpXs(1, id_depl) = BurnUpXs(1, id_depl) + 0.5_8 * DeplXs%xsa(i)
      BurnUpXs(2, id_depl) = BurnUpXs(2, id_depl) + 0.5_8 * DeplXs%xsf(i)
      BurnUpXs(3, id_depl) = BurnUpXs(3, id_depl) + 0.5_8 * DeplXs%xsn2n(i)
      BurnUpXs(4, id_depl) = BurnUpXs(4, id_depl) + 0.5_8 * DeplXs%xsn3n(i)
    ENDDO
    DeplVars%phi1g = DeplVars%phi1g + 0.5*DeplXs%phi1g
  ENDIF
ELSE
  IF(.NOT. lCorrectStep) THEN
    IdIso => Fxr%DeplXs1g(-1)%IdIso; nIso = Fxr%DeplXs1g(-1)%niso
    w = - DeplCntl%tsec/DeplCntl%ptsec/2._8
    DO i = 1, niso
      j = IdIso(i); id_lib = iPosIso(j)
      id_depl = MapXs2Dep(id_lib)
      BurnUpXs(1, id_depl) = w * Fxr%DeplXs1g(-1)%xsa(i)
      BurnUpXs(2, id_depl) = w * Fxr%DeplXs1g(-1)%xsf(i)
      BurnUpXs(3, id_depl) = W * Fxr%DeplXs1g(-1)%xsn2n(i)
      BurnUpXs(4, id_depl) = w * Fxr%DeplXs1g(-1)%xsn3n(i)
    ENDDO
    DeplVars%Phi1g = w * Fxr%DeplXs1g(-1)%phi1g
    IdIso => Fxr%IdIso; nIso = Fxr%nIso_depl
    w =  DeplCntl%tsec/DeplCntl%ptsec/2._8 + 1
    DO i = 1, niso
      j = IdIso(i); id_lib = iPosIso(j)
      id_depl = MapXs2Dep(id_lib)
      BurnUpXs(1, id_depl) = BurnUpXs(1, id_depl) + w * DeplXs%xsa(i)
      BurnUpXs(2, id_depl) = BurnUpXs(2, id_depl) + w * DeplXs%xsf(i)
      BurnUpXs(3, id_depl) = BurnUpXs(3, id_depl) + w * DeplXs%xsn2n(i)
      BurnUpXs(4, id_depl) = BurnUpXs(4, id_depl) + w * DeplXs%xsn3n(i)
    ENDDO
    DeplVars%phi1g = DeplVars%phi1g + w * DeplXs%phi1g
  ELSE
    IdIso => Fxr%DeplXs1g(-1)%IdIso; nIso = Fxr%DeplXs1g(-1)%niso
#define QuadInt
#ifdef QuadInt
    w0 = w_QI
#else
    w0 = w_LI
#endif
    w  = w0(1)
    DO i = 1, niso
      j = IdIso(i); id_lib = iPosIso(j)
      id_depl = MapXs2Dep(id_lib)
      BurnUpXs(1, id_depl) = w * Fxr%DeplXs1g(-1)%xsa(i)
      BurnUpXs(2, id_depl) = w * Fxr%DeplXs1g(-1)%xsf(i)
      BurnUpXs(3, id_depl) = w * Fxr%DeplXs1g(-1)%xsn2n(i)
      BurnUpXs(4, id_depl) = w * Fxr%DeplXs1g(-1)%xsn3n(i)
    ENDDO
    DeplVars%Phi1g = w * Fxr%DeplXs1g(-1)%phi1g

    w  = w0(2)
    IdIso => Fxr%DeplXs1g(0)%IdIso; nIso = Fxr%DeplXs1g(0)%niso
    DO i = 1, niso
      j = IdIso(i); id_lib = iPosIso(j)
      id_depl = MapXs2Dep(id_lib)
      BurnUpXs(1, id_depl) = BurnUpXs(1, id_depl) + w * Fxr%DeplXs1g(0)%xsa(i)
      BurnUpXs(2, id_depl) = BurnUpXs(2, id_depl) + w * Fxr%DeplXs1g(0)%xsf(i)
      BurnUpXs(3, id_depl) = BurnUpXs(3, id_depl) + w * Fxr%DeplXs1g(0)%xsn2n(i)
      BurnUpXs(4, id_depl) = BurnUpXs(4, id_depl) + w * Fxr%DeplXs1g(0)%xsn3n(i)
    ENDDO
    DeplVars%Phi1g = DeplVars%Phi1g + w * Fxr%DeplXs1g(0)%phi1g

    w  = w0(3)
    IdIso => Fxr%IdIso; nIso = Fxr%nIso_depl
    DO i = 1, niso
      j = IdIso(i); id_lib = iPosIso(j)
      id_depl = MapXs2Dep(id_lib)
      BurnUpXs(1, id_depl) = BurnUpXs(1, id_depl) + w * DeplXs%xsa(i)
      BurnUpXs(2, id_depl) = BurnUpXs(2, id_depl) + w * DeplXs%xsf(i)
      BurnUpXs(3, id_depl) = BurnUpXs(3, id_depl) + w * DeplXs%xsn2n(i)
      BurnUpXs(4, id_depl) = BurnUpXs(4, id_depl) + w * DeplXs%xsn3n(i)
    ENDDO
    DeplVars%phi1g = DeplVars%phi1g + w * DeplXs%phi1g
  ENDIF
ENDIF

END SUBROUTINE




SUBROUTINE SaveDeplXs_GD(DeplXs, Fxr, lCorrectStep)
USE PARAM
USE TYPEDEF,      ONLY : FxrInfo_Type
USE DeplType,      ONLY : DeplXs_Type,      DeplVars_Type
USE nuclidmap_mod, ONLY : iposiso
USE BasicOperation, ONLY : CP_CA, CP_VA
IMPLICIT NONE
TYPE(DeplXs_Type) :: DeplXs
TYPE(FxrInfo_Type) :: Fxr
LOGICAL :: lCorrectStep
INTEGER :: niso, i

IF(.NOT. lCorrectStep) THEN
  Fxr%DeplXs1g(0)%idiso=0
  Fxr%DeplXs1g(0)%xsa=0; Fxr%DeplXs1g(0)%xsf=0
  Fxr%DeplXs1g(0)%xsn2n=0; Fxr%DeplXs1g(0)%xsn3n=0

  Fxr%DeplXs1g(0)%phi1g = DeplXs%phi1g
  Fxr%DeplXs1g(0)%niso = Fxr%nIso_depl
  niso=Fxr%nIso_depl
  CALL CP_VA(Fxr%DeplXs1g(0)%idiso(1:niso),Fxr%IdIso(1:niso), niso)
  CALL CP_VA(Fxr%DeplXs1g(0)%xsa(1:niso), DeplXs%xsa(1:niso), niso)
  CALL CP_VA(Fxr%DeplXs1g(0)%xsf(1:niso), DeplXs%xsf(1:niso), niso)
  CALL CP_VA(Fxr%DeplXs1g(0)%xsn2n(1:niso), DeplXs%xsn2n(1:niso), niso)
  CALL CP_VA(Fxr%DeplXs1g(0)%xsn3n(1:niso), DeplXs%xsn3n(1:niso), niso)
  DO i = 1, niso
    IF(Fxr%idiso(i) .NE. 64155) CYCLE
    Fxr%DeplXs1g(0)%n155 = Fxr%pnum(i)
  ENDDO
ELSE
  Fxr%DeplXs1g(-1)%idiso=0
  Fxr%DeplXs1g(-1)%xsa=0; Fxr%DeplXs1g(-1)%xsf=0
  Fxr%DeplXs1g(-1)%xsn2n=0; Fxr%DeplXs1g(-1)%xsn3n=0

  Fxr%DeplXs1g(-1)%phi1g = Fxr%DeplXs1g(0)%phi1g
  Fxr%DeplXs1g(-1)%niso = Fxr%DeplXs1g(0)%niso
  Fxr%DeplXs1g(-1)%n155 = Fxr%DeplXs1g(0)%n155
  niso = Fxr%DeplXs1g(-1)%niso

  CALL CP_VA(Fxr%DeplXs1g(-1)%idiso(1:niso), Fxr%DeplXs1g(0)%IdIso(1:niso), niso)
  CALL CP_VA(Fxr%DeplXs1g(-1)%xsa(1:niso), Fxr%DeplXs1g(0)%xsa(1:niso), niso)
  CALL CP_VA(Fxr%DeplXs1g(-1)%xsf(1:niso), Fxr%DeplXs1g(0)%xsf(1:niso), niso)
  CALL CP_VA(Fxr%DeplXs1g(-1)%xsn2n(1:niso), Fxr%DeplXs1g(0)%xsn2n(1:niso), niso)
  CALL CP_VA(Fxr%DeplXs1g(-1)%xsn3n(1:niso), Fxr%DeplXs1g(0)%xsn3n(1:niso), niso)

ENDIF
END SUBROUTINE


SUBROUTINE ConstDeplVas(DeplVars, DeplXs, Fxr, nIsoLib, nIsoDepl)
USE PARAM
USE TypeDef,     ONLY : FxrInfo_Type
USE DeplType,    ONLY : DeplXs_Type,      DeplVars_Type
USE nuclidmap_mod,ONLY : iposiso
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
TYPE(DeplVars_Type) ::DeplVars
TYPE(DeplXs_Type) ::DeplXs
TYPE(FxrInfo_Type) :: Fxr
!INTEGER :: IdIso(nIsoLib)
!REAL :: pnum(nIsoLib)
INTEGER :: nIsoLib, nIsoDepl

INTEGER, POINTER :: IdIso(:)
REAL, POINTER :: pnum(:)
INTEGER, POINTER :: MapXs2Dep(:)
REAL, POINTER :: BurnUpXs(:, :), IsoNum(:)
INTEGER :: i, j, id_lib, id_depl, niso

MapXs2Dep => DeplVars%MapXs2Dep
BurnUpXs => DeplVars%BurnUpXs
IsoNum => DeplVars%IsoNum
CALL CP_CA(IsoNum, 0._8, nIsoDepl)
CALL CP_CA(BurnUpXS, 0._8, 4, nIsoDepl)
IdIso => Fxr%IdIso_past; pnum => Fxr%pnum_past
nIso = Fxr%nIso_past
!--beg 16/02/11 Depletion timestep bug fixed
IF(.NOT. Fxr%l_pnum_all) THEN
  DO i= 1, nIso
    j = IdIso(i); id_lib = iPosIso(j)
    id_depl = MapXs2Dep(id_lib)
    IF(id_depl .EQ. 0) CYCLE
    Fxr%pnum_past_all(id_depl) = pnum(i)
  ENDDO
  Fxr%l_pnum_all = .TRUE.
ENDIF

DO i = 1, nIsoDepl
  IsoNum(i) = Fxr%pnum_past_all(i)
ENDDO
! DO i= 1, nIso
!   j = IdIso(i); id_lib = iPosIso(j)
!   id_depl = MapXs2Dep(id_lib)
!   IF(id_depl .EQ. 0) CYCLE
!   IsoNum(id_depl) = pnum(i)
! ENDDO
!--end 16/02/11 Depletion timestep bug fixed
#ifdef Depl_mod0

#endif
#ifdef Depl_mod1
NULLIFY(Idiso, Pnum)
IdIso => Fxr%IdIso; nIso = Fxr%nIso_depl
#endif

DO i= 1, nIso
  j = IdIso(i); id_lib = iPosIso(j)
  id_depl = MapXs2Dep(id_lib)
  IF(id_depl .EQ. 0) CYCLE
  !BurnUpXs(1, id_depl) = DeplXs%xsa(i)
  BurnUpXs(1, id_depl) = DeplXs%xsa(i) + DeplXs%xsn2n(i) + DeplXs%xsn3n(i)
  BurnUpXs(2, id_depl) = DeplXs%xsf(i)
  BurnUpXs(3, id_depl) = DeplXs%xsn2n(i)
  BurnUpXs(4, id_depl) = DeplXs%xsn3n(i)
  !BurnUpXs(1, id_depl) = BurnUpXs(1, id_depl) + DeplXs%xsn2n(i) + 2*DeplXs%xsn3n(i)
ENDDO
DeplVars%phi1g = DeplXs%phi1g
END SUBROUTINE

SUBROUTINE GdPostCorrection(Fxr, DeplCntl)
USE PARAM
USE TYPEDEF,       ONLY : FxrInfo_Type,        DeplPCP_Type
USE DeplType,      ONLY : DeplCntl_Type
IMPLICIT NONE
TYPE(FxrInfo_Type) :: Fxr
TYPE(DeplCntl_Type) :: DeplCntl

INTEGER :: niso
INTEGER :: id, iso
REAL :: f, f1, f2
REAL :: CorrectPnum(64152:64160)

IF(DeplCntl%NowStep ==1 .AND. DeplCntl%lPredict) THEN
  Fxr%DeplPCP%pnum(2, 0, :) = 1.e-40_8
  niso = Fxr%niso_past
  DO iso = 1, niso
    id = Fxr%idiso_past(iso)
    IF(id .LT. 64152 .OR. id .GT. 64160) CYCLE
    Fxr%DeplPCP%pnum(2, 0, id) = Fxr%pnum_past(iso)
  ENDDO
ENDIF

IF(DeplCntl%lPredict) THEN
  niso = Fxr%niso_depl
  Fxr%DeplPCP%pnum(1, 1, :) = 1.e-40_8
  DO iso = 1, niso
    id = Fxr%idiso(iso)
    IF(id .LT. 64152 .OR. id .GT. 64160) CYCLE
    Fxr%DeplPCP%pnum(1, 1, id) = Fxr%pnum(iso)
  ENDDO
ELSE
  niso = Fxr%niso_depl
  Fxr%DeplPCP%pnum(2, 1, :) = 1.e-40_8
  DO iso = 1, niso
    id = Fxr%idiso(iso)
    IF(id .LT. 64152 .OR. id .GT. 64160) CYCLE
    Fxr%DeplPCP%pnum(2, 1, id) = Fxr%pnum(iso)
  ENDDO
ENDIF

IF(DeplCntl%NowStep .GT. 1 .AND. DeplCntl%lPredict) THEN
  CorrectPnum = 0
  DO id = 64152, 64160
    IF(id .EQ. 64153) CYCLE
    IF(id .EQ. 64159) CYCLE
    f = Fxr%DeplPCP%pnum(1, 0, id) * Fxr%DeplPCP%pnum(2, 0, id)
    f = f * Fxr%DeplPCP%pnum(1, -1, id) * Fxr%DeplPCP%pnum(2, -1, id)
    IF(abs(f) .LT. epsm20) THEN
      Fxr%DeplPCP%f(id) = 1
      CYCLE
    ENDIF
    f = Fxr%DeplPCP%pnum(1, 0, id) / Fxr%DeplPCP%pnum(2, -1, id)
    IF(abs(f-1._8) .LT. epsm20) THEN
      Fxr%DeplPCP%f(id) = 1
      CYCLE
    ENDIF
    f1 = log(f)

    f = Fxr%DeplPCP%pnum(2, 0, id) / Fxr%DeplPCP%pnum(2, -1, id)
    IF(abs(f-1._8) .LT. epsm20) THEN
      Fxr%DeplPCP%f(id) = 1
      CYCLE
    ENDIF
    f2 = log(f)
    Fxr%DeplPCP%f(id) = f2 / f1
  ENDDO

  DO id = 64152, 64160
    IF(id .EQ. 64153) CYCLE
    IF(id .EQ. 64159) CYCLE
    !IF(Fxr%DeplPCP%pnum(2, 0, id) .LT. epsm10) CYCLE
    f = Fxr%DeplPCP%pnum(1, 1, id) * Fxr%DeplPCP%pnum(2, 0, id)
    IF(abs(f) .LT. epsm20) CYCLE
    f = Fxr%DeplPCP%pnum(1, 1, id) / Fxr%DeplPCP%pnum(2, 0, id)
    IF(abs(f-1._8) .LT. epsm20) CYCLE
    f1 = log(f) * Fxr%DeplPCP%f(id)
    CorrectPnum(id) = exp(f1) * Fxr%DeplPCP%pnum(2, 0, id)
  ENDDO

  niso = Fxr%niso_depl
  DO iso = 1, niso
    id = Fxr%idiso(iso)
    IF(id .LT. 64152 .OR. id .GT. 64160) CYCLE

    IF(id .EQ. 64153) CYCLE
    IF(id .EQ. 64159) CYCLE
    IF(CorrectPnum(id) .LT. 1.e-30) CYCLE
    f = CorrectPnum(id) / Fxr%pnum(iso)
    IF(f>1.5) THEN
      CorrectPnum(id) = 1.5 * Fxr%pnum(iso)
    ELSEif(f<0.5) THEN
      CorrectPnum(id) =  0.5 * Fxr%pnum(iso)
    ENDIF
    Fxr%pnum(iso) = CorrectPnum(id)
    !Fxr%DeplPCP%pnum(1, 1, id) = CorrectPnum(id)
  ENDDO
ENDIF
IF(.NOT. DeplCntl%lPredict) THEN
  DO id = 64152, 64160
    Fxr%DeplPCP%pnum(:, -1, id) = Fxr%DeplPCP%pnum(:, 0, id)
    Fxr%DeplPCP%pnum(:, 0, id) = Fxr%DeplPCP%pnum(:, 1, id)
  ENDDO
ENDIF
END SUBROUTINE


SUBROUTINE SetDeplCoreState(ibu, CoreState, ThInfo, ThVar, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,  ONLY : ThInfo_Type,      ThVar_Type, PE_TYPE
USE DeplType, ONLY : CoreState_Type
USE CNTL,     ONLY : nTracerCntl_Type
USE FILES,         ONLY : io8
USE IOUTIL,        ONLY : message
IMPLICIT NONE
INTEGER :: ibu
TYPE(CoreState_Type) :: CoreState
TYPE(THInfo_Type) :: ThInfo
TYPE(THVar_Type) :: ThVar
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
CHARACTER(100) :: boronppm

REAL :: powlv, flowlv
IF(CoreState%lBoronSearch(ibu)) THEN
  nTracerCntl%lBoronSearch = .TRUE.
  nTracerCntl%Target_eigv = CoreState%Target_Keff(ibu)
ELSE
  IF(CoreState%BoronPPM(ibu) .LT. -50.0) CoreState%BoronPPM(ibu) = CoreState%BoronPPM(ibu-1)
  nTracerCntl%BoronPpm = CoreState%BoronPPM(ibu)
ENDIF

nTracerCntl%PowerLevel = CoreState%RelPow(ibu)
ThInfo%PowLv = CoreState%RelPow(ibu)

IF(nTracerCntl%lBoronSearch) THEN
  WRITE(boronppm, '(x, A)')  'SEARCH'
ELSE
  WRITE(boronppm, '(F8.2 x, A)') nTracerCntl%BoronPpm, 'ppm'
ENDIF

PowLv = CoreState%RelPow(ibu) * 100
flowlv = 100.0

WRITE(mesg, '(A)') 'Set Core State ...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)
WRITE(mesg, '(5x, A, F7.2, A2, x, A, F7.2, A2, x, A, A12)') 'Power Level :', nTracerCntl%PowerLevel*100,' %',        &
      'Flow Rate Level :', flowlv , '%', 'Boron:', boronppm
IF(PE%Master) CALL message(io8, FALSE, TRUE, mesg)

END SUBROUTINE

SUBROUTINE UpdtFluxSolution(ibu, DeplCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,         ONLY : PE_Type
USE DeplType,        ONLY : DeplCntl_Type
USE CNTL,            ONLY : nTracerCntl_Type
USE geom,            ONLY : Core
USE Core_mod,        ONLY : FmInfo,             GroupInfo,            THInfo,           &
                            xkconv,             xkcrit
USE XeDyn_Mod,       ONLY : UpdtXeDyn,                                                  &
                            InitXe,             UpdtXe,               FinalXe
USE CritSpec_mod,    ONLY : XsKinf
USE FILES,           ONLY : io8
USE IOUTIL,          ONLY : message

IMPLICIT NONE
INTEGER :: ibu
TYPE(DeplCntl_Type) :: DeplCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

WRITE(mesg, '(A)') 'Update Steady-State Solution for Core State Change ...'
IF(PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)
IF(DeplCntl%lTrXe .AND. DeplCntl%lXeDyn) nTracerCntl%lXeDyn = .FALSE.
CALL SSEIG()
xkcrit = XsKinf(Core, FmInfo, GroupInfo, .TRUE., nTracerCntl%lCritSpec, PE)
IF(DeplCntl%lTrXe .AND. DeplCntl%lXeDyn) nTracerCntl%lXeDyn = .TRUE.
IF(nTracerCntl%lTrXe) CALL UpdtXeDyn(Core, FmInfo, THInfo, GroupInfo, nTracerCntl, DeplCntl, PE, Mode = InitXe)
WRITE(mesg,'(A,2(F9.3,A),2(A,F9.5),A,F8.2,A)') 'Burnup = ',DeplCntl%T_mwdkghm(IBU-1),  &
      ' MWD/kgHM',DeplCntl%T_efpd(IBU-1),' EPFD','  keff = ',xkconv,'  kinf(BK) = ', &
        xkcrit,'  ppm = ',nTracerCntl%boronppm,' ppm'
IF(PE%Master) CALL message(io8, FALSE, TRUE, mesg)

WRITE(mesg, '(A)') hbar1(1:77)
IF(PE%Master) CALL message(io8, FALSE, TRUE, MESG)
END SUBROUTINE

SUBROUTINE SaveCoreState(ibu, CoreState, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,  ONLY : ThInfo_Type,      ThVar_Type, PE_TYPE
USE DeplType, ONLY : CoreState_Type
USE CNTL,     ONLY : nTracerCntl_Type
USE FILES,         ONLY : io8
USE IOUTIL,        ONLY : message
IMPLICIT NONE
INTEGER :: ibu
TYPE(CoreState_Type) :: CoreState
TYPE(THInfo_Type) :: ThInfo
TYPE(THVar_Type) :: ThVar
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
CoreState%BoronPpm(ibu) = nTracerCntl%boronPpm
CoreState%RelPow(ibu) = nTracerCntl%PowerLevel
END SUBROUTINE

SUBROUTINE EffBurnUpTime(Tsec, PowLv0, PowLv1)
USE PARAM
IMPLICIT NONE
REAL :: Tsec, PowLv0, PowLv1
REAL :: DelPow, EffFactor
delpow = PowLv1 - PowLv0
EffFactor = 1./PowLv0
IF(DelPow .GT. epsm6) THEN
  EffFactor = (log(PowLv1/ PowLv0)) / delpow
ENDIF
Tsec = Tsec * EffFactor
PRINT *, EffFactor
END SUBROUTINE

SUBROUTINE DeplBanner(imode)
USE PARAM
USE TYPEDEF,         ONLY : CoreInfo_Type,      FmInfo_Type,        PE_Type
USE DeplType,        ONLY : DeplCntl_Type
USE geom,            ONLY : Core
USE Core_mod,        ONLY : FmInfo,                                            &
                            xkconv,             xkcrit
USE CNTL,            ONLY : nTracerCntl
USE PE_Mod,          ONLY : PE
USE Depl_Mod,        ONLY : DeplCntl,            GetCoreHmMass
USE FILES,           ONLY : io8
USE IOUTIL,          ONLY : message
IMPLICIT NONE

INTEGER :: imode
INTEGER :: iburnup

IF(imode .EQ. 1) THEN
  !WRITE(mesg, '(A)') 'Depletion Calculation ...'
  !IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)

  WRITE(mesg, '(A)') '=== Initial Isotope Summary ==='
  IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)
  DeplCntl%Hm_Mass0_kg = Core%Hm_Mass0
  IF (.NOT. nTracerCntl%lnTIGRst) THEN                                             ! --- 180827 JSR
    IF(Core%Hm_Mass0 .GT. 0._8) THEN                                               ! --- 180827 JSR
      DeplCntl%Hm_Mass0_ton = Core%Hm_Mass0                                        ! --- 180827 JSR
    ELSE                                                                           ! --- 180827 JSR
      DeplCntl%Hm_Mass0_ton = GetCoreHmMass(Core, FmInfo%Fxr, PE)                  ! --- 180827 JSR
    ENDIF                                                                          ! --- 180827 JSR
    DeplCntl%Hm_Mass0_kg = DeplCntl%Hm_Mass0_ton * 1000._8
  END IF
  write(mesg,'(a,1pe12.4)') 'Initial HM Loading (kg) = ',DeplCntl%Hm_Mass0_kg
  IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)
ELSEIF(imode .EQ. 2)THEN
  iburnup = DeplCntl%NowStep
  WRITE(mesg,'(A,2(F9.3,A),2(A,F9.5),A,F8.2,A)') 'Burnup = ',DeplCntl%T_mwdkghm(iBurnUp),  &
       ' MWD/kgHM',DeplCntl%T_efpd(iBurnUp),' EPFD','  keff = ',xkconv,'  kinf(BK) = ', &
         xkcrit,'  ppm = ',nTracerCntl%boronppm,' ppm'
  IF(PE%Master) CALL message(io8, FALSE, TRUE, mesg)
  !WRITE(mesg, '(A)') '=== Final Isotope Summary ==='
  !IF(PE%Master) CALL message(io8, FALSE, TRUE, mesg)
  !WRITE(mesg, '(A)') hbar2(1:77)
  !IF(PE%Master) CALL message(io8, FALSE, TRUE, MESG)
ENDIF
END SUBROUTINE

SUBROUTINE SetDeplTimeStep(ibu,Core, DeplCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type, PE_Type
USE DeplType,       ONLY : DeplCntl_Type
USE CNTL,           ONLY : nTracerCntl_Type
USE FILES,          ONLY : io8
USE IOUTIL,         ONLY : message
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(DeplCntl_Type) :: DeplCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
INTEGER :: ibu

INTEGER :: itype
REAL :: PrevBurnup(2)

IF(ibu .EQ. 0) THEN
  DeplCntl%NowBurnUp(2) = 0.00001_8
  DeplCntl%NowBurnUp(1) = DeplCntl%NowBurnUp(2) / DeplCntl%PowerCore * DeplCntl%Hm_Mass0_kg
  DeplCntl%TSec = DeplCntl%NowBurnUp(1) * 86400._8  !Day Unit => Sec Unit
  DeplCntl%NowStep = 0
  DeplCntl%lPredict = .TRUE.
  RETURN
ENDIF

itype = DeplCntl%BurnUpType
PrevBurnUp = DeplCntl%NowBurnUp
DeplCntl%NowBurnUp(1) = DeplCntl%T_efpd(ibu)
DeplCntl%NowBurnUp(2) = DeplCntl%T_mwdkghm(ibu)
DeplCntl%TSec =  (DeplCntl%NowBurnUp(1) - PrevBurnUp(1)) * 86400._8
DeplCntl%TSec = DeplCntl%TSec / nTracerCntl%PowerLevel
DeplCntl%NowStep = ibu
Core%DelT = DeplCntl%TSec

if(itype .EQ. 1) WRITE(mesg,'(a,f12.2,a, 5x, a, i5, a, i5, a)') 'Burnup = ',DeplCntl%NowBurnUp(1),' EPFD', '(',iBu, ' /',DeplCntl%nBurnUpStep,' )'
if(itype .EQ. 2) WRITE(mesg,'(a,f12.4,a, 5x, a, i5, a, i5, a)') 'Burnup = ',DeplCntl%NowBurnUp(2),' MWD/kgHM', '(',iBu, ' /',DeplCntl%nBurnUpStep,' )'
IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)
END SUBROUTINE
