#include <defines.h>
#include <DefDBG.h>

SUBROUTINE SetDeplXeDynMode(Core, FmInfo, DeplCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,   FmInfo_Type,  PE_Type
USE DeplType,         ONLY : DeplCntl_Type
USE CNTL,             ONLY : nTracerCntl_Type
USE XeDyn_Mod,        ONLY : InitXeDynInfo
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(DeplCntl_Type) :: DeplCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
!Xenon Dynamics
IF(nTracerCntl%lXeDyn) THEN
  DeplCntl%lXeDyn = nTracerCntl%lXeDyn
  DeplCntl%lTrXe = nTracerCntl%lTrXe
  DeplCntl%lEqXe = nTracerCntl%lEqXe
ELSE
  DeplCntl%lXeDyn = .TRUE.
  DeplCntl%lTrXe = .FALSE.
  DeplCntl%lEqXe = .TRUE.
ENDIF
nTracerCntl%lXeDyn = DeplCntl%lXeDyn
nTracerCntl%lTrXe = DeplCntl%lTrXe
nTracerCntl%lEqXe = DeplCntl%lEqXe

CALL InitXeDynInfo(FmInfo%Fxr, Core%nCoreFxr, PE%myzb, PE%myze)

END SUBROUTINE


SUBROUTINE InitXeDynInfo(Fxr, nfxr, myzb, myze)
USE PARAM
USE TYPEDEF,         ONLY : FxrInfo_Type,     XeDynFxr_Type
!USE CNTL,            ONLY : nTracerCntl_Type

TYPE(FxrInfo_Type) :: Fxr(nfxr, myzb:myze)
INTEGER :: nfxr, myzb, myze

TYPE(XeDynFxr_Type), POINTER :: XeDyn
INTEGER :: iz, ifxr, i
INTEGER :: niso, niso_depl, niso0
REAL, POINTER :: pnum(:)
INTEGER, POINTER :: idiso(:)

LOGICAL :: lXe, lI

DO iz = myzb, myze
  DO ifxr = 1, nfxr
    IF(.NOT. Fxr(ifxr, iz)%lfuel) CYCLE
    IF(.NOT. ASSOCIATED(Fxr(ifxr, iz)%XeDyn)) THEN
       ALLOCATE(Fxr(ifxr, iz)%XeDyn)
    ENDIF
    XeDyn => Fxr(ifxr, iz)%XeDyn
    pnum => Fxr(ifxr, iz)%pnum
    idiso => Fxr(ifxr, iz)%idiso
    niso0 = Fxr(ifxr, iz)%niso
    niso = niso0; niso_depl = Fxr(ifxr, iz)%niso_depl
    lXe = .FALSE.; lI = .FALSE.
    DO i = 1, niso0
      IF(idiso(i) .EQ. 53635) lI = .TRUE.
      IF(idiso(i) .EQ. 54635) lXe = .TRUE.
    ENDDO
    IF(.NOT. lI) THEN
      niso = niso + 1; niso_depl = niso_depl + 1
      idiso(niso) = 53635
      pnum(niso) = 1.E-30_8
    ENDIF
    IF(.NOT. lXe) THEN
      niso = niso + 1; niso_depl = niso_depl + 1
      idiso(niso) = 54635
      pnum(niso) = 1.E-30_8
    ENDIF
    Fxr(ifxr, iz)%niso = niso; Fxr(ifxr, iz)%niso_depl = niso_depl
    NULLIFY(XeDyn, pnum, idiso)
  ENDDO
ENDDO

END SUBROUTINE


SUBROUTINE UpdtXeDyn(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, DeplCntl, PE, Mode)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,       FmInfo_Type,       GroupInfo_Type,    &
                             THInfo_Type,         PE_Type,                              &
                             FxrInfo_Type,        Pin_Type,          Cell_Type,         &
                             XsMac_Type
USE DeplType,         ONLY : DeplCntl_Type
USE CNTL,             ONLY : nTracerCntl_Type
USE FILES,           ONLY : io8
USE IOUTIL,          ONLY : message
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(DeplCntl_Type) :: DeplCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
INTEGER, OPTIONAL :: Mode
INTEGER, PARAMETER :: DefaultMode = 1
INTEGER :: ActMode
IF(.NOT. Present(Mode)) THEN
  ActMode = DefaultMode
ELSE
  ActMode = Mode
ENDIF



WRITE(mesg, '(A)') 'Update Xenon Dynamics ...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)
SELECT CASE(nTracerCntl%lProblem)
CASE(lsseigv)
  IF(nTracerCntl%lEqXe) THEN
    CALL  EqXeUpdate(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE)
  ENDIF
CASE(lDepletion, lXenonDynamics)
  IF(.NOT. nTracerCntl%lXeDyn) RETURN
  IF(nTracerCntl%lTrXe) THEN
    CALL TrXeUpdate(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, DeplCntl, PE, ActMode)
  ELSEIF(nTracerCntl%lEqXe) THEN
    CALL EqXeUpdate(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE)
  ENDIF
END SELECT

END SUBROUTINE

SUBROUTINE TrXeUpdate(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, DeplCntl, PE, Mode)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,       FmInfo_Type,       GroupInfo_Type,    &
                             THInfo_Type,         PE_Type,                              &
                             FxrInfo_Type,        Pin_Type,          Cell_Type,         &
                             XsMac_Type
USE DeplType,         ONLY : DeplCntl_Type
USE CNTL,             ONLY : nTracerCntl_Type
USE Depl_mod,         ONLY : FluxNormalizeFactor,decayXe135,decayI135,yieldXe135,yieldI135
USE MacXsLib_mod,     ONLY : EffMacXS,            MacXsBase
USE NuclidMap_mod,    ONLY : lFissile!,            yieldxesm,                            &
                             !decayxe135,          decayI135
USE XsUtil_mod,       ONLY : AllocXsMac,          GetXsMacDat,       ReturnXsMacDat,    &
                             FreeXsIsoMac,        FreeXsMac
USE XSLIB_MOD,        ONLY : mapnucl
USE XeDyn_Mod,        ONLY : UpdtXe,              InitXe,            FinalXe
USE BasicOperation,   ONLY : CP_CA, CP_VA, MULTI_VA
USE TH_Mod,           ONLY : GetPinFuelTemp
USE OMP_LIB
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(DeplCntl_Type) :: DeplCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
INTEGER, INTENT(IN) :: Mode

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :), myFxr
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
!TYPE(XsMac_Type), SAVE :: XsMac
REAL, POINTER :: Phis(:, :, :)

REAL :: NormFactor, dt
REAL :: lamXe, lamI
INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr
INTEGER :: ng, nFsr, nFxr, nxy, nz, myzb, myze, xyb, xye, norg   !--- CNJ Edit : Domain Decomposition + MPI
INTEGER :: iresoGrp1, iresoGrp2
INTEGER :: iz, ixy, icel, ifxr
INTEGER :: j
ng = GroupInfo%ng
iresoGrp1 = GroupInfo%nofg + 1; iresoGrp2 = GroupInfo%nofg + GroupInfo%norg

myzb = PE%myzb; myze = PE%myze
nFsr = Core%nCoreFsr;   nFxr = Core%nCoreFxr
nxy = Core%nxy

!--- CNJ Edit : Domain Decomposition + MPI
xyb = PE%myPinBeg; xye = PE%myPinEnd
IF (PE%RTMASTER) THEN
  xyb = 1; xye = nxy
ENDIF

dt = Core%DelT
Fxr => FmInfo%Fxr; Phis => FmInfo%Phis
Pin => Core%Pin; CellInfo => Core%CellInfo
NormFactor = FluxNormalizeFactor(Core, FmInfo, GroupInfo, nTracerCntl%PowerCore, FALSE, TRUE, PE)
NormFactor = NormFactor * nTracerCntl%PowerLevel  !Adjusting Power Level

!lamI = DecayI135(); lamXe = DecayXe135()
lamI = decayI135; lamXe = decayXe135

!IF(.NOT. XsMac%lAlloc) THEN
!  XsMac%ng = ng
!  CALL AllocXsMac(XsMac)
!ENDIF


!CALL GetXsMacDat(XsMac,ng, .TRUE.)
CALL OMP_SET_NUM_THREADS(PE%nDeplThread)
DO iz = myzb, myze
!$OMP PARALLEL DEFAULT(SHARED)&
!$OMP PRIVATE(FsrIdxSt, FxrIdxSt, icel, nLocalFxr, j, ifxr, myFxr)
!$OMP DO SCHEDULE(GUIDED)
  DO ixy = xyb, xye
    FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
    icel = Pin(ixy)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      IF(.NOT. Fxr(ifxr, iz)%lFuel) CYCLE
      myFxr => Fxr(ifxr, iz)
      CALL FxrWiseTrXe(myFxr, iz, icel, j, FsrIdxSt)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL
ENDDO

!CALL ReturnXsMacDat(XsMac)
NULLIFY(PIN); NULLIFY(CellINfo); NULLIFY(myFxr)

CONTAINS

  SUBROUTINE FxrWiseTrXe(myFxr, iz, icel, ifxrlocal, FsrIdxSt)
  IMPLICIT NONE
  TYPE(Fxrinfo_type), POINTER :: myFxr
  INTEGER :: iz, icel, ifxrlocal, FsrIdxSt

  INTEGER :: niso, niso_past, nFsrInFxr
  INTEGER :: ig, i, ism, ixe, iI, ifsr, ifsrlocal
  TYPE(XsMac_Type) :: XsMac
  INTEGER, PARAMETER :: MaxGrp = 1000

  REAL :: PhiFxr(MaxGrp)
  REAL, POINTER :: IsoXsMacf(:, :), IsoXsMacA(:, :)
  REAL, POINTER :: pnum(:), pnum_past(:)
  INTEGER, POINTER :: idiso(:), idiso_past(:)

  REAL :: fisrate, AbsRate, Area
  REAL :: ydi, ydxe, ydpm, Prate_Xe135, Prate_I135
  REAL :: ExpI, ExpXe
  REAL :: PnI, PnXe, PnI0, PnXe0

  LOGICAL :: lres

  nFsrInFxr = CellInfo(icel)%nFsrInFxr(ifxrlocal)
  niso = myFxr%niso; niso_past = myFxr%niso_past
  CALL MacXsBase(XsMac, myFxr, 1, ng, ng, 1._8, FALSE, TRUE, TRUE)
  idiso => myFxr%idiso; pnum => myFxr%pnum; lres = myFxr%lres
  idiso_past => myFxr%idiso_past; pnum_past => myFxr%pnum_past
  IsoXsMacf => XsMac%IsoXsMacf
  IsoXsMacA => XsMac%IsoXsMacA
  IF(myFxr%lRes) THEN
      DO ig = iResoGrp1, iResoGrp2
          DO i = 1, niso
              IsoXsMacF(i,ig) = IsoXsMacF(i,ig) * myFXR%fresoFIso(i,ig)
              IsoXsMacA(i,ig) = IsoXsMacA(i,ig) * myFXR%fresoAIso(i,ig)
          ENDDO
      ENDDO
  ENDIF
  PhiFxr(1:ng) = 0
  DO i = 1, nFsrInFxr
    ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(i, ifxrlocal)
    iFsr = FsrIdxSt + ifsrlocal - 1
    Area = CellInfo(icel)%vol(ifsrlocal)
    PhiFxr(1:ng) = PhiFxr(1:ng) + Area * Phis(iFsr, iz, 1:ng)
  ENDDO
  PhiFxr(1:ng) = PhiFxr(1:ng) / myFxr%Area

  Prate_Xe135 = 0; Prate_I135 = 0
  DO i = 1, niso
    IF(idiso(i) .EQ. 54635) ixe = i
    IF(idiso(i) .EQ. 53635) ii = i
    IF(idiso(i) .EQ. 62649) ism = i
    IF(.NOT. lfissile(idiso(i))) CYCLE
    fisrate = 0
    DO ig = 1, ng
      FisRate = FisRate + IsoXsMacf(i, ig) * PhiFxr(ig)
    ENDDO
    !CALL yieldxesm(IdIso(i), ydi, ydxe, ydpm)
    ydxe = yieldXe135(mapnucl(IdIso(i)))
    ydi = yieldI135(mapnucl(IdIso(i)))
    Prate_Xe135 = Prate_Xe135 + ydxe * FisRate
    Prate_I135 = Prate_I135 + ydi * FisRate
  ENDDO
  Prate_I135 = NormFactor * Prate_I135
  Prate_Xe135 = NormFactor * Prate_Xe135

  AbsRate = 0
  DO ig = 1, ng
    AbsRate = AbsRate + IsoXsMacA(ixe, ig) * PhiFxr(ig)
  ENDDO
  AbsRate = NormFactor * AbsRate /pnum(ixe)

  !Save Production Rate Information
  myFxr%XeDyn%ARate_Xe135(2) = AbsRate            !Xe 135 Absorption Rate
  myFxr%XeDyn%Prate_I135(2) = Prate_I135
  myFxr%XeDyn%Prate_Xe135(2) = Prate_Xe135

  IF(Mode .EQ. UpdtXe .OR. Mode .EQ. FinalXe) THEN
    PnI0 = epsm30; PnXe0 = epsm30
    IF(nTracerCntl%lProblem .EQ. lDepletion) THEN
    DO i = 1, niso_past
      IF(idiso_past(i) .EQ. 53635) THEN
        PnI0 = pnum_past(i)
      ENDIF
      IF(idiso_past(i) .EQ. 54635) THEN
        PnXe0 = pnum_past(i)
      ENDIF
    ENDDO
    ELSEIF(nTracerCntl%lProblem .EQ. lXenonDynamics) THEN
      PnI0 = myFxr%XeDyn%Pn_I135(1)
      PnXe0 = myFxr%XeDyn%Pn_Xe135(1)
    ENDIF
#define AvgPrate
#ifdef AvgPrate
    AbsRate = 0.5_8 * (AbsRate + myFxr%XeDyn%ARate_Xe135(1))
    Prate_Xe135 = 0.5_8 * (Prate_Xe135 + myFxr%XeDyn%PRate_Xe135(1))
    Prate_I135 = 0.5_8 * (Prate_I135 + myFxr%XeDyn%PRate_I135(1))
#endif
    ExpI = exp(-dt*lamI)
    ExpXe = exp(-dt*(lamXe+1.e-24_8*AbsRate))
    PnI = PnI0 * ExpI + 1.e-24_8 * Prate_I135 * (1._8-ExpI) / LamI
    PnXe = PnXe0*ExpXe
    PnXe = PnXe + (Prate_I135 + Prate_Xe135) / (LamXe*1.E+24_8+AbsRate)*(1._8-ExpXe)
    PnXe = PnXe + (1.e+24_8*LamI*PnI0-Prate_I135) / (1.e+24_8*(LamXe-LamI)+AbsRate) * (ExpI-ExpXe)

    pnum(ixe) = PnXe
    pnum(ii) = PnI
  ENDIF

  IF(Mode .EQ. InitXe .OR. Mode .EQ. FinalXe) THEN
    myFxr%XeDyn%ARate_Xe135(1) = myFxr%XeDyn%ARate_Xe135(2)           !Xe 135 Absorption Rate
    myFxr%XeDyn%Prate_I135(1) = myFxr%XeDyn%Prate_I135(2)
    myFxr%XeDyn%Prate_Xe135(1) = myFxr%XeDyn%Prate_Xe135(2)
    myFxr%XeDyn%Pn_I135 = pnum(iI)
    myFxr%XeDyn%Pn_Xe135 = pnum(ixe)
  ELSEIF(Mode .EQ. UpdtXe) THEN
    myFxr%XeDyn%Pn_I135(2) = pnum(iI)
    myFxr%XeDyn%Pn_Xe135(2) = pnum(ixe)
  ENDIF

  CALL FreeXsIsoMac(XsMac)
  CALL FreeXsMac(XsMac)

  END SUBROUTINE
END SUBROUTINE

!--- CNJ Edit : Unused
!--- PHS : Also Invalid for the resonance treatment
SUBROUTINE TrXeUpdate2(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, DeplCntl, PE, Mode)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,       FmInfo_Type,       GroupInfo_Type,    &
                             THInfo_Type,         PE_Type,                              &
                             FxrInfo_Type,        Pin_Type,          Cell_Type,         &
                             XsMac_Type
USE DeplType,         ONLY : DeplCntl_Type
USE CNTL,             ONLY : nTracerCntl_Type
USE Depl_mod,         ONLY : FluxNormalizeFactor,decayXe135,decayI135,yieldXe135,yieldI135
USE MacXsLib_mod,     ONLY : EffMacXS,            MacXsBase
USE NuclidMap_mod,    ONLY : lFissile!,            yieldxesm,                            &
                             !decayxe135,          decayI135
USE XsUtil_mod,       ONLY : AllocXsMac,          GetXsMacDat,       ReturnXsMacDat
USE XSLIB_MOD,        ONLY : mapnucl
USE XeDyn_Mod,        ONLY : UpdtXe,              InitXe,            FinalXe
USE BasicOperation,   ONLY : CP_CA, CP_VA, MULTI_VA
USE TH_Mod,           ONLY : GetPinFuelTemp
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(DeplCntl_Type) :: DeplCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
INTEGER, INTENT(IN) :: Mode

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :), myFxr
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(XsMac_Type), SAVE :: XsMac
REAL, POINTER :: Phis(:, :, :)

INTEGER, PARAMETER :: MaxGrp = 1000

REAL :: PhiFxr(MaxGrp)
REAL, POINTER :: SigNfOld(:, :)
REAL, POINTER :: EQXS(:, :), IsoXsMacf(:, :), IsoXsMacA(:, :)
REAL, POINTER :: pnum(:), pnum_past(:), SubGrpLv(:, :)
INTEGER, POINTER :: idiso(:), idiso_past(:)
REAL :: NormFactor, fisrate, FisRate0, AbsRate, AbsRate0, Area, Temp
REAL :: yd, ydi, ydxe, ydpm, lamXe, lamI, Prate_Xe135, Prate_I135, Prate0_Xe135, Prate0_I135
REAL :: a0, a1, c0, c1, c2
REAL :: R0, R1, Q0, Q1, G0, dIdt, dQdT, F1, F2
REAL :: ExpI, ExpXe, dt
REAL :: PnI, PnXe, PnI0, PnXe0
INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nLocalFsr, nFsrInFxr
INTEGER :: ng, nFsr, nFxr, nxy, nz, myzb, myze, norg, nchi, niso, niso_past
INTEGER :: iresoGrp1, iresoGrp2
INTEGER :: iz, ixy, icel, ifxr, ifsr, ifsrlocal, ism, ixe, iI
INTEGER :: i, j, ig, ig2
LOGICAL :: lXsLib, lRes
ng = GroupInfo%ng
iresoGrp1 = GroupInfo%nofg + 1; iresoGrp2 = GroupInfo%nofg + GroupInfo%norg

myzb = PE%myzb; myze = PE%myze
nFsr = Core%nCoreFsr;   nFxr = Core%nCoreFxr
nxy = Core%nxy
dt = Core%DelT
Fxr => FmInfo%Fxr; Phis => FmInfo%Phis
Pin => Core%Pin; CellInfo => Core%CellInfo
NormFactor = FluxNormalizeFactor(Core, FmInfo, GroupInfo, nTracerCntl%PowerCore, FALSE, TRUE, PE)
NormFactor = NormFactor * nTracerCntl%PowerLevel  !Adjusting Power Level

!lamI = DecayI135(); lamXe = DecayXe135()
lamI = decayI135; lamXe = decayXe135

IF(.NOT. XsMac%lAlloc) THEN
  XsMac%ng = ng
  CALL AllocXsMac(XsMac)
ENDIF


!CALL GetXsMacDat(XsMac,ng, .TRUE.)

DO iz = myzb, myze
  DO ixy = 1, nxy
    FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
    icel = Pin(ixy)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      IF(.NOT. Fxr(ifxr, iz)%lFuel) CYCLE
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
      myFxr => Fxr(ifxr, iz)
      niso = myFxr%niso; niso_past = myFxr%niso_past
      CALL MacXsBase(XsMac, myFxr, 1, ng, ng, 1._8, FALSE, TRUE, TRUE)
      idiso => myFxr%idiso; pnum => myFxr%pnum; lres = myFxr%lres
      idiso_past => myFxr%idiso_past; pnum_past => myFxr%pnum_past
      IF(myFxr%lRes) THEN
        temp = myFxr%temp
        ALLOCATE(SigNfOld(niso, ng))

        SigNfOld(1:niso, 1:ng) = XsMac%IsoXsMacnf(1:niso, 1:ng)

        DO ig = iResoGrp1, iResoGrp2
          DO i = 1, niso
            IF(SigNfOld(i, ig) * XsMac%IsoXsMacNf(i, ig) .NE. 0._8) THEN
              XsMac%IsoXsMacf(i, ig) = XsMac%IsoXsMacf(i, ig) * XsMac%IsoXsMacNf(i, ig) / SigNfOld(i, ig)
            ENDIF
          ENDDO
        ENDDO
        DEALLOCATE(SigNfOld)
      ENDIF
      IsoXsMacf => XsMac%IsoXsMacf
      IsoXsMacA => XsMac%IsoXsMacA
      PhiFxr(1:ng) = 0
      DO i = 1, nFsrInFxr
        ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(i, j)
        iFsr = FsrIdxSt + ifsrlocal - 1
        Area = CellInfo(icel)%vol(ifsrlocal)
        PhiFxr(1:ng) = PhiFxr(1:ng) + Area * Phis(iFsr, iz, 1:ng)
      ENDDO
      PhiFxr(1:ng) = PhiFxr(1:ng) / myFxr%Area

      Prate_Xe135 = 0; Prate_I135 = 0
      DO i = 1, niso
        IF(idiso(i) .EQ. 54635) ixe = i
        IF(idiso(i) .EQ. 53635) ii = i
        IF(idiso(i) .EQ. 62649) ism = i
        IF(.NOT. lfissile(idiso(i))) CYCLE
        fisrate = 0
        DO ig = 1, ng
          FisRate = FisRate + IsoXsMacf(i, ig) * PhiFxr(ig)
        ENDDO
        !CALL yieldxesm(IdIso(i), ydi, ydxe, ydpm)
        ydxe = yieldXe135(mapnucl(IdIso(i)))
        ydi = yieldI135(mapnucl(IdIso(i)))
        Prate_Xe135 = Prate_Xe135 + ydxe * FisRate
        Prate_I135 = Prate_I135 + ydi * FisRate
      ENDDO
      Prate_I135 = NormFactor * Prate_I135
      Prate_Xe135 = NormFactor * Prate_Xe135

      AbsRate = 0
      DO ig = 1, ng
        AbsRate = AbsRate + IsoXsMacA(ixe, ig) * PhiFxr(ig)
      ENDDO
      AbsRate = NormFactor * AbsRate /pnum(ixe)

      !Save Production Rate Information
      Fxr(ifxr, iz)%XeDyn%ARate_Xe135(2) = AbsRate            !Xe 135 Absorption Rate
      Fxr(ifxr, iz)%XeDyn%Prate_I135(2) = Prate_I135
      Fxr(ifxr, iz)%XeDyn%Prate_Xe135(2) = Prate_Xe135

      IF(Mode .EQ. UpdtXe .OR. Mode .EQ. FinalXe) THEN
        PnI0 = epsm30; PnXe0 = epsm30
        IF(nTracerCntl%lProblem .EQ. lDepletion) THEN
        DO i = 1, niso_past
          IF(idiso_past(i) .EQ. 53635) THEN
            PnI0 = pnum_past(i)
          ENDIF
          IF(idiso_past(i) .EQ. 54635) THEN
            PnXe0 = pnum_past(i)
          ENDIF
        ENDDO
        ELSEIF(nTracerCntl%lProblem .EQ. lXenonDynamics) THEN
          PnI0 = Fxr(ifxr, iz)%XeDyn%Pn_I135(1)
          PnXe0 = Fxr(ifxr, iz)%XeDyn%Pn_Xe135(1)
        ENDIF

        PnI0 = PnI0 * 1.0E+24_8; PnXe0 = PnXe0 *1.0E+24_8
        AbsRate = AbsRate *1.E-24_8; AbsRate0 = Fxr(ifxr, iz)%XeDyn%ARate_Xe135(1) * 1.0E-24_8
        Prate0_Xe135 = Fxr(ifxr, iz)%XeDyn%PRate_Xe135(1); Prate0_I135 = Fxr(ifxr, iz)%XeDyn%PRate_I135(1)


        a0 = Prate0_I135
        a1 = (Prate_I135 - Prate0_I135) / dt
        ExpI = exp(-dt*lamI)

        c0 = (LamI * a0 - a1) / (LamI * LamI)
        c1 = a1 / LamI; c2 = PnI0 - c0
        PnI = c0 + c1 * dt + c2 * ExpI

        R0 = Prate0_Xe135 + LamI * PnI0; R1 = Prate_Xe135 + LamI * PnI
        Q0 = -(LamXe + AbsRate0); Q1 = -(LamXe + AbsRate)
        G0 = R0 + Q0 * PnXe0

        F1 = Dt * 0.5_8
        F2 = Dt * Dt / 12._8

        dIdT = Prate_I135 - Prate0_I135 - LamI * (PnI - PnI0)
        dQdT = (Q1 - Q0) / dt

        PnXe = PnXe0 + F1 * (R1 + G0)
        PnXe = PnXe + F2 * (LamI * dIdT - Q0 * G0 + Q1 * R1 -dQdT * PnXe0)
        PnXe = PnXe / (1 - F1 * Q1 - F2 *(Q1 * Q1 + dQdT))

        PnXe = PnXe * 1.0E-24_8;
        PnI = PnI * 1.0E-24_8;

        pnum(ixe) = PnXe
        pnum(Ii) = PnI
!#define AvgPrate
!#ifdef AvgPrate
!        AbsRate = 0.5_8 * (AbsRate + Fxr(ifxr, iz)%XeDyn%ARate_Xe135(1))
!        Prate_Xe135 = 0.5_8 * (Prate_Xe135 + Fxr(ifxr, iz)%XeDyn%PRate_Xe135(1))
!        Prate_I135 = 0.5_8 * (Prate_I135 + Fxr(ifxr, iz)%XeDyn%PRate_I135(1))
!#endif

!        ExpXe = exp(-dt*(lamXe+1.e-24_8*AbsRate))
!        PnI = PnI0 * ExpI + 1.e-24_8 * Prate_I135 * (1._8-ExpI) / LamI
!        PnXe = PnXe0*ExpXe
!        PnXe = PnXe + (Prate_I135 + Prate_Xe135) / (LamXe*1.E+24_8+AbsRate)*(1._8-ExpXe)
!        PnXe = PnXe + (1.e+24_8*LamI*PnI0-Prate_I135) / (1.e+24_8*(LamXe-LamI)+AbsRate) * (ExpI-ExpXe)
!
!        pnum(ixe) = PnXe
!        pnum(ii) = PnI
      ENDIF

      IF(Mode .EQ. InitXe .OR. Mode .EQ. FinalXe) THEN
        Fxr(ifxr, iz)%XeDyn%ARate_Xe135(1) = Fxr(ifxr, iz)%XeDyn%ARate_Xe135(2)           !Xe 135 Absorption Rate
        Fxr(ifxr, iz)%XeDyn%Prate_I135(1) = Fxr(ifxr, iz)%XeDyn%Prate_I135(2)
        Fxr(ifxr, iz)%XeDyn%Prate_Xe135(1) = Fxr(ifxr, iz)%XeDyn%Prate_Xe135(2)
        Fxr(ifxr, iz)%XeDyn%Pn_I135 = pnum(iI)
        Fxr(ifxr, iz)%XeDyn%Pn_Xe135 = pnum(ixe)
      ELSEIF(Mode .EQ. UpdtXe) THEN
        Fxr(ifxr, iz)%XeDyn%Pn_I135(2) = pnum(iI)
        Fxr(ifxr, iz)%XeDyn%Pn_Xe135(2) = pnum(ixe)
      ENDIF
    ENDDO
  ENDDO
ENDDO

!CALL ReturnXsMacDat(XsMac)
NULLIFY(PIN); NULLIFY(CellINfo); NULLIFY(myFxr)
NULLIFY(idiso); NULLIFY(pnum); NULLIFY(EQXS)
NULLIFY(IsoXsMacA);NULLIFY(IsoXsMacF)
END SUBROUTINE

SUBROUTINE EqXeUpdate(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,       FmInfo_Type,       GroupInfo_Type,    &
                             THInfo_Type,         PE_Type,                              &
                             FxrInfo_Type,        Pin_Type,          Cell_Type,         &
                             XsMac_Type
USE CNTL,             ONLY : nTracerCntl_Type
USE Depl_mod,         ONLY : FluxNormalizeFactor,decayXe135,decayI135,yieldXe135,yieldI135
USE MacXsLib_mod,     ONLY : EffMacXS,            MacXsBase, MacXsAF_CrCsp
USE NuclidMap_mod,    ONLY : lFissile!,            yieldxesm,         decayxe135,        &
                             !decayI135
USE XsUtil_mod,       ONLY : AllocXsMac,          GetXsMacDat,       ReturnXsMacDat,     &
                             FreeXsMac,           FreeXsIsoMac
USE XSLIB_MOD,        ONLY : mapnucl
USE BasicOperation,   ONLY : CP_CA, CP_VA, MULTI_VA
USE TH_Mod,           ONLY : GetPinFuelTemp
USE XeDyn_Mod,        ONLY : lRelaxation
USE OMP_LIB
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE


TYPE(FxrInfo_Type), POINTER :: Fxr(:, :), myFxr
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(XsMac_Type), SAVE :: XsMac
REAL, POINTER :: Phis(:, :, :)

!REAL, POINTER :: SigFOld(:, :)!SigNfOld(:, :)
REAL :: NormFactor
INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nLocalFsr
INTEGER :: ng, nFsr, nFxr, nxy, nz, myzb, myze, xyb, xye, norg, nchi  !--- CNJ Edit : Domain Decomposition + MPI
INTEGER :: iresoGrp1, iresoGrp2
INTEGER :: iz, ixy, icel, ifxr
INTEGER :: j
LOGICAL :: lXsLib, FlagF

REAL, PARAMETER :: f_relax = 0.3 ! Edit by LHG, 2021-0107

IF(nTracerCntl%lReadXeEQ) RETURN
lXsLib = nTracerCntl%lXsLib
IF(.NOT. lXsLib) RETURN

ng = GroupInfo%ng
iresoGrp1 = GroupInfo%nofg + 1; iresoGrp2 = GroupInfo%nofg + GroupInfo%norg

myzb = PE%myzb; myze = PE%myze
nFsr = Core%nCoreFsr;   nFxr = Core%nCoreFxr
nxy = Core%nxy

!--- CNJ Edit : Domain Decomposition + MPI
xyb = PE%myPinBeg; xye = PE%myPinEnd
IF (PE%RTMASTER) THEN
  xyb = 1; xye = nxy
ENDIF

Fxr => FmInfo%Fxr; Phis => FmInfo%Phis
Pin => Core%Pin; CellInfo => Core%CellInfo
NormFactor = FluxNormalizeFactor(Core, FmInfo, GroupInfo, nTracerCntl%PowerCore, FALSE, TRUE, PE)
NormFactor = NormFactor * nTracerCntl%PowerLevel  !Adjusting Power Level

!CALL GetXsMacDat(XsMac,ng, .TRUE.)

!IF(.NOT. XsMac%lAlloc) THEN
!  XsMac%ng = ng
!  CALL AllocXsMac(XsMac)
!ENDIF

CALL OMP_SET_NUM_THREADS(PE%nDeplThread)
DO iz = myzb, myze
!$OMP PARALLEL DEFAULT(SHARED)&
!$OMP PRIVATE(FsrIdxSt, FxrIdxSt, icel, nLocalFxr, j, ifxr, myFxr)
!$OMP DO SCHEDULE(GUIDED)
  DO ixy = xyb, xye
    FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
    icel = Pin(ixy)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      IF(.NOT. Fxr(ifxr, iz)%lFuel) CYCLE
      myFxr => Fxr(ifxr, iz)
      CALL FxrWiseEqXe(myFxr, iz, icel, j, FsrIdxSt)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL
ENDDO
!CALL ReturnXsMacDat(XsMac)
CONTAINS
  SUBROUTINE FxrWiseEqXe(myFxr, iz, icel, ifxrlocal, FsrIdxSt)
  IMPLICIT NONE
  TYPE(FxrInfo_Type), POINTER :: myFxr
  INTEGER :: iz, icel, ifxrlocal, FsrIdxSt

  INTEGER :: i, ig, ism, ixe, iI, ifsr, ifsrlocal
  INTEGER :: niso, nFsrInFxr
  REAL, POINTER :: IsoXsMacf(:, :), IsoXsMacA(:, :)!EQXS(:, :)
  REAL, POINTER :: pnum(:)!, SubGrpLv(:, :)
  INTEGER, POINTER :: idiso(:)
  TYPE(XsMac_Type) :: XsMac
  INTEGER, PARAMETER :: MaxGrp = 1000

  REAL :: PhiFxr(MaxGrp)
  REAL :: yd, ydi, ydxe, ydpm, pn, pnI, pneq, lambda, lambdaI, fisrate, AbsRate, Area
  REAL :: pnD, pnID, absD

  REAL :: locf_relax

  nFsrInFxr = CellInfo(icel)%nFsrInFxr(ifxrlocal)
  niso = myFxr%niso;
  ALLOCATE(IsoXsMacA(niso,ng), IsoXsMacF(niso,ng))
  !CALL MacXsBase(XsMac, myFxr, 1, ng, ng, 1._8, FALSE, TRUE, TRUE)
  CALL MacXsAF_CrCSP(IsoXsMacA,IsoXsMacf, myFxr, 1, ng, ng, TRUE)
  idiso => myFxr%idiso; pnum => myFxr%pnum

  IF (lRelaxation) THEN
    locf_relax = f_relax
  ELSE
    locf_relax = 0.0
  END IF

  !IsoXsMacf => XsMac%IsoXsMacf
  !IsoXsMacA => XsMac%IsoXsMacA
  IF(myFxr%lRes) THEN
      DO ig = iResoGrp1, iResoGrp2
          DO i = 1, niso
              IsoXsMacF(i,ig) = IsoXsMacF(i,ig) * myFXR%fresoFIso(i,ig)
              IsoXsMacA(i,ig) = IsoXsMacA(i,ig) * myFXR%fresoAIso(i,ig)
          ENDDO
      ENDDO
  ENDIF
  PhiFxr(1:ng) = 0
  DO i = 1, nFsrInFxr
    ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(i, ifxrlocal)
    iFsr = FsrIdxSt + ifsrlocal - 1
    Area = CellInfo(icel)%vol(ifsrlocal)
    PhiFxr(1:ng) = PhiFxr(1:ng) + Area * Phis(iFsr, iz, 1:ng)
  ENDDO
  PhiFxr(1:ng) = PhiFxr(1:ng) / myFxr%Area
  pn = 0; pnI = 0
  DO i = 1, niso
    IF(idiso(i) .EQ. 54635) ixe = i
    IF(idiso(i) .EQ. 62649) ism = i
    IF(idiso(i) .EQ. 53635) iI = i
    IF(.NOT. lfissile(idiso(i))) CYCLE
    fisrate = 0
    DO ig = 1, ng
      FisRate = FisRate + IsoXsMacf(i, ig) * PhiFxr(ig)
    ENDDO
    !CALL yieldxesm(IdIso(i), ydi, ydxe, ydpm)
    ydxe = yieldXe135(mapnucl(IdIso(i)))
    ydi = yieldI135(mapnucl(IdIso(i)))
    yd = ydi + ydxe  !Sum yeild of I-135 + Xe-134
    pn = pn + yd * FisRate
    pnI = pnI + ydi * FisRate
  ENDDO
  pn = Pn * NormFactor
  pnI = pnI * NormFactor
  AbsRate = 0
  DO ig = 1, ng
    AbsRate = AbsRate + IsoXsMacA(ixe, ig) * PhiFxr(ig)
  ENDDO
  AbsRate = NormFactor * AbsRate
  AbsRate = AbsRate / pnum(ixe)  ! Macro -> Micro
  lambda = decayxe135!()
  lambdaI = decayI135!()

  pnD = myFxr%pnXe; pnID = myFxr%pnI; absD = myFxr%absXe
  myFxr%pnXe = pn; myFxr%pnI = pnI; myFxr%absXe = AbsRate
  pn = pnD*locf_relax+(1.-locf_relax)*pn
  AbsRate = absD*locf_relax+(1.-locf_relax)*AbsRate
  !pnum(ixe) = locf_relax*pnum(ixe) + (1.-locf_relax) * pn / (AbsRate + 1.e24_8 * lambda)
  pnum(ixe) = pn / (AbsRate + 1.e24_8 * lambda)

  pnI = pnID*locf_relax+(1.-locf_relax)*pnI
  !pn = ydi * FisRate * NormFactor
  !pnum(iI) = locf_relax*pnum(iI) + (1.-locf_relax) * pnI / (1.e24_8 * lambdaI)
  pnum(iI) = pnI / (1.e24_8 * lambdaI)
  !Fxr(ifxr, iz)%XeDyn%Prate_I135(1) = PnI

  !CALL FreeXsIsoMac(XsMac)
  !CALL FreeXsMac(XsMac)
  DEALLOCATE(IsoXsMacA,IsoXsMacf)
  END SUBROUTINE
END SUBROUTINE

SUBROUTINE ReadXeND(Core, FmInfo, PE)
USE PARAM
USE TYPEDEF,    ONLY : CoreInfo_TYpe,    FmInfo_Type,   PE_TYPE,           &
                       FxrInfo_Type
USE IOUTIL,     ONLY : openfile,         terminate
USE Files,      ONLY : FileName,         XeEqFileIdx,   io29
#ifdef MPI_ENV
USE MPIComm_mod, ONLY : MPIWaitTurn
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(PE_Type) :: PE

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)

CHARACTER(80) :: localfn
LOGICAL :: Master
INTEGER :: i, iz, ifxr, ndum, nXeFxr0, nXeFxr
INTEGER :: idum(1:20)
REAL :: rdum(1:20)

master = PE%CmfdMaster
Fxr => FmInfo%FXR

nXeFxr0= 0
DO iz = PE%myzb, PE%myze
  DO ifxr = 1, Core%nCoreFxr
    IF(.NOT. Fxr(ifxr, iz)%lfuel) CYCLE
    DO i = 1, Fxr(ifxr, iz)%niso
      IF(Fxr(ifxr, iz)%idiso(i) .NE. 54635) CYCLE
      nXeFxr0 = nXeFxr0 + 1
    ENDDO
  ENDDO
ENDDO

localfn = FileName(XeEqFileIdx)
#ifdef MPI_ENV
CALL MPIWaitTurn(PE%MPI_Cmfd_Comm, PE%myCmfdRank, PE%nCmfdproc, .FALSE.)
#endif
OPEN(io29, file=localfn, status='OLD', form = 'unformatted')

READ(io29) idum(1), idum(2), idum(3)
IF(idum(1) .NE. Core%nz) CALL TERMINATE('Read XE EQ : Inconsistent NZ')
IF(idum(2) .NE. Core%nCoreFsr) CALL TERMINATE('Read XE EQ : Inconsistent nCoreFsr')
IF(idum(3) .NE. Core%nCoreFxr) CALL TERMINATE('Read XE EQ : Inconsistent nCoreFxr')

nXeFxr = 0
DO
  READ(io29, end=1000) idum(1), idum(2), rdum(1)
  iz = idum(1)
  IF(iz .LT. PE%myzb) CYCLE
  IF(iz .GT. PE%myze) EXIT
  ifxr = idum(2)
  IF(.NOT. Fxr(ifxr, iz)%lfuel) CALL TERMINATE('Read XE EQ : Inconsistent Geometry')
  DO i = 1, Fxr(ifxr, iz)%niso
    IF(Fxr(ifxr, iz)%idiso(i) .NE. 54635) CYCLE
    Fxr(ifxr,iz)%pnum(i) = rdum(1)
    nXeFxr = nXeFxr + 1
    EXIT
  ENDDO
ENDDO
1000 continue

IF(nXeFxr .NE. nXeFxr0) CALL TERMINATE('Read XE EQ: Inconsistent Geometry')

#ifdef MPI_ENV
CALL MPIWaitTurn(PE%MPI_Cmfd_Comm, PE%myCmfdRank, PE%nCmfdproc, .TRUE.)
#endif
END SUBROUTINE

SUBROUTINE WriteXeND(Core, FmInfo, PE)
USE PARAM
USE TYPEDEF,    ONLY : CoreInfo_TYpe,    FmInfo_Type,   PE_TYPE,           &
                       FxrInfo_Type
USE IOUTIL,     ONLY : openfile
USE Files,      ONLY : caseid,           io29
#ifdef MPI_ENV
USE MPIComm_mod, ONLY : MPIWaitTurn
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(PE_Type) :: PE

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)

CHARACTER(128) :: filename
LOGICAL :: Master
INTEGER :: i, iz, ifxr, ndum
INTEGER :: idum(1:20)
REAL :: rdum(1:20)

master = PE%CmfdMaster
Fxr => FmInfo%FXR
filename = trim(caseid)//'.xeeq'
#ifdef MPI_ENV
CALL MPIWaitTurn(PE%MPI_Cmfd_Comm, PE%myCmfdRank, PE%nCmfdproc, .FALSE.)
#endif
IF(MASTER) THEN
  OPEN(io29, file=filename, status='REPLACE', form = 'unformatted')
ELSE
  OPEN(io29, file=filename, status='OLD', form = 'unformatted')
ENDIF
IF(MASTER) WRITE(io29) Core%nz, Core%nCoreFsr, Core%nCoreFxr
IF(.NOT. MASTER) THEN
  READ(io29) idum(1), idum(2), idum(3)
  DO
    READ(io29, end=1000) idum(1), idum(2), rdum(1)
  ENDDO
1000 continue
ENDIF
DO iz = PE%myzb, PE%myze
  DO ifxr = 1, Core%nCoreFxr
    IF(.NOT. Fxr(ifxr, iz)%lfuel) CYCLE
    DO i = 1, Fxr(ifxr, iz)%niso
      IF(Fxr(ifxr, iz)%idiso(i) .NE. 54635) CYCLE
      WRITE(io29) iz, ifxr, Fxr(ifxr,iz)%pnum(i)
    ENDDO
  ENDDO
ENDDO
CLOSE(io29)
#ifdef MPI_ENV
CALL MPIWaitTurn(PE%MPI_Cmfd_Comm, PE%myCmfdRank, PE%nCmfdproc, .TRUE.)
#endif
END SUBROUTINE
