#include <defines.h>
SUBROUTINE AllocTransient()
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,     FmInfo_Type,       CmInfo_Type,       PE_TYPE,     &
                             TranCntl_Type,     TranInfo_Type
USE GEOM,             ONLY : Core
USE CORE_MOD,         ONLY : FmInfo,            CmInfo,            GroupInfo,          ThInfo,     & 
                             TranPhi,           TranPhiCm,         TranPhiFm,                      &
                             Prec,              PrecCm,            PrecFm,                         &
                             TranPsi,           TranPsiCm,         TranPsiFm,                      &
                             TranPsid,          TranPsiCmd,        TranPsiFmd,                     &
                             ResSrcCm,          ResSrcFm,                                          &
                             GCTranSrc,                                                            &
                             !Adaptive Theta                                                       &
                             ThetaCM,                                                              &
                             !BDF                                                                  &
                             TranPhi2,          TranPhi3,          TranPhi4,           TranPhi5,   &
                             TranPhiCm2,        TranPhiCm3,        TranPhiCm4,         TranPhiCm5, &
                             TranPhiFm2,        TranPhiFm3,        TranPhiFm4,         TranPhiFm5, &
                             !SCM                                                                  &
                             TranPrec,                                                             &
                             ShpFrqCm,          ShpFrqFm,          ShpFrqCmd,          ShpFrqFmd,  &
                             AvgShpFrqCm,       AvgShpFrqFm,       PrecFrqCm,          PrecFrqFm,  &
                             PhiSCM,            PhiCSCM,           PsiSCM,             PsiCSCM,    &
                             xsnfSCM,                                                              &
                             !AfSrc                                                                &
                             TranPhi1a,         TranPhi2a,                                         &
                             !AM3
                             ResSrcCmd,         ResSrcFmd
USE TRAN_MOD,         ONLY : TranCntl,          TranInfo,          XsChange,                       &
                             SetTimeStep,       XsNoise,           XsCntlRod
USE TH_MOD,           ONLY : ThVar
USE TRANCMFD_MOD,     ONLY : TrSrc,             PrecSrc
USE RAYS,             ONLY : RayInfo

!USE MPIAxSolver_Mod,  ONLY : AllocTranAxNSolver

USE CNTL,             ONLY : nTRACERCntl
USE PE_Mod,           ONLY : PE
USE ALLOCS 
IMPLICIT NONE

INTEGER :: myzbf, myzef, myzb, myze, nxy, nfsr, nfxr, nz
INTEGER :: nprec, ng, ngc
INTEGER :: iz, ixy
INTEGER :: nPolarAng, nAziAng


myzb = PE%myzb; myze = PE%myze; nz = PE%nz
myzbf = PE%myzbf; myzef = PE%myzef
nPrec = GroupInfo%nPrec; ng = GroupInfo%ng; ngc = GroupInfo%ngc
nxy = Core%nxy; nFsr = Core%nCoreFsr; nFxr = Core%nCoreFxr
nPolarAng = RayInfo%nPolarAngle; nAziAng = RayInfo%nAziAngle

!Fm Data
CALL DMALLOC0(TranPhi, 1, nFsr, myzb, myze, 1, ng)
CALL DMALLOC0(Prec, 1, nPrec, 1, nFsr, myzb, myze)
CALL DMALLOC0(TranPsi, 1, nFsr, myzb, myze)
CALL DMALLOC0(TranPsid, 1, nFsr, myzb, myze)

FmInfo%TranPhi => TranPhi
FmInfo%Prec => Prec
FmInfo%TranPsi => TranPsi
FmInfo%TranPsid => TranPsid

IF(TranCntl%MOC_BDF) THEN
  CALL DMALLOC0(TranPhi2, 1, nFsr, myzb, myze, 1, ng)
  CALL DMALLOC0(TranPhi3, 1, nFsr, myzb, myze, 1, ng)
  CALL DMALLOC0(TranPhi4, 1, nFsr, myzb, myze, 1, ng)
  CALL DMALLOC0(TranPhi5, 1, nFsr, myzb, myze, 1, ng)

  FmInfo%TranPhi2 => TranPhi2
  FmInfo%TranPhi3 => TranPhi3
  FmInfo%TranPhi4 => TranPhi4
  FmInfo%TranPhi5 => TranPhi5
END IF

CALL Dmalloc0(FmInfo%TranPower, 1, nFsr, myzb, myze)

!IF(TranCntl%lCusping) THEN
!  ALLOCATE(FmInfo%neighPhis(nFsr, ng, 2))
!END IF

!Cm Data
CALL DMALLOC0(TranPhiCm, 1, nxy, myzb, myze, 1, ng)
CALL DMALLOC0(PrecCm, 1, nPrec, 1, nxy, myzb, myze)
CALL DMALLOC0(TranPsiCm, 1, nxy, myzb, myze)
CALL DMALLOC0(TranPsiCmd, 1, nxy, myzb, myze)

IF(nTracerCntl%lSubPlane) THEN
  CALL DMALLOC0(TranPhiFm, 1, nxy, myzbf, myzef, 1, ng)
  CALL DMALLOC0(PrecFm, 1, nPrec, 1, nxy, myzbf, myzef)
  CALL DMALLOC0(TranPsiFm, 1, nxy, myzbf, myzef)
  CALL DMALLOC0(TranPsiFmd, 1, nxy, myzbf, myzef)
ELSE
  TranPhiFm => TranPhiCm
  PrecFm => PrecCm
  TranPsiFm => TranPsiCm
  TranPsiFmd => TranPsicmd
ENDIF

CmInfo%TranPhiCm => TranPhiCm; CmInfo%TranPhiFm => TranPhiFm
CmInfo%TranPsiCm => TranPsiCm; CmInfo%TranPsiFm => TranPsiFm
CmInfo%TranPsiCmd => TranPsiCmd; CmInfo%TranPsiFmd => TranPsiFmd
CmInfo%PrecFm => PrecFm; CmInfo%PrecCm => PrecCm

IF(nTracerCntl%lGcCmfd) THEN
  CALL Dmalloc0(GcTranSrc, 1, nxy, myzb, myze, 1, ngc)
  CmInfo%GcTranSrc => GcTranSrc
ENDIF

!
CALL Dmalloc0(ResSrcCm, 1, nxy, myzb, myze, 1, ng)
IF(nTracerCntl%lSubPlane) THEN
  CALL Dmalloc0(ResSrcFm, 1, nxy, myzbf, myzef, 1, ng)
ELSE
  ResSrcFm => ResSrcCm
ENDIF
CmInfo%ResSrcCm => ResSrcCm; CmInfo%ResSrcFm => ResSrcFm
FmInfo%ResSrc => ResSrcCm
!CALL Dmalloc0(CmInfo%ResSrcCm, 1, nxy, myzb, myze, 1, ng)


TranInfo%nPrec = nPrec
CALL DMALLOC(TranInfo%Chid, ng)
CALL DMALLOC(TranInfo%Chidk, ng, nprec)
CALL DMALLOC(TranInfo%Lambda, nPrec)
CALL DMALLOC(TranInfo%InvLambda, nPrec)
!CALL DMALLOC(TranInfo%Neut_velo, nPrec)
!CALL DMALLOC(TranInfo%Kappa, nPrec)

CALL DMALLOC0(TranInfo%CellOmegam, 0, nPrec, 1, nxy, myzb, myze)
CALL DMALLOC0(TranInfo%Cellomega0, 0, nPrec, 1, nxy, myzb, myze)
CALL DMALLOC0(TranInfo%CellOmegap, 0, nPrec, 1, nxy, myzb, myze)

CALL DMALLOC0(TranInfo%FxrOmegam, 0, nPrec, 1, nFxr, myzb, myze)
CALL DMALLOC0(TranInfo%FxrOmega0, 0, nPrec, 1, nFxr, myzb, myze)
CALL DMALLOC0(TranInfo%FxrOmegap, 0, nPrec, 1, nFxr, myzb, myze)

CALL AllocTranFxr(Core, FmInfo, GroupInfo, nTracerCntl, PE)

CALL DMALLOC0(TranInfo%Expo, 1, nxy, myzb, myze, 1, ng)
CALL DMALLOC0(TranInfo%Expo_alpha, 1, nxy, myzb, myze, 1, ng)
CALL DMALLOC0(TranInfo%FmExpo, 1, nFsr, myzb, myze, 1, ng)
CALL DMALLOC0(TranInfo%FmExpo_alpha, 1, NFsr, myzb, myze, 1, ng)
CALL DMALLOC0(TranInfo%RefTemp0, 1, nz); CALL DMALLOC0(TranInfo%RefTemp, 1, nz)
DO iz = myzb, myze
  DO ixy = 1, nxy
    CALL DMALLOC0(CmInfo%PinXs(ixy, iz)%beta, 1, nPrec)
    CALL DMALLOC0(CmInfo%PinXs(ixy, iz)%chip, 1, ng)
    CALL DMALLOC0(CmInfo%PinXs(ixy, iz)%velo, 1, ng)
    CALL DMALLOC0(CmInfo%PinXs(ixy, iz)%rvdelt, 1, ng)
    IF(nTracerCntl%lGcCmfd) CALL Dmalloc0(CmInfo%GcPinXs(ixy, iz)%velo, 1, ngc)
  ENDDO
ENDDO

!IF(nTracerCntl%l3dim) THEN
!  CALL AllocTranAxNSolver(ng, nPrec, PE)
!ENDIF
TranCntl%XsChange => XsChange
TranCntl%XsNoise => XsNoise
TranCntl%XsCntlRod => XsCntlRod



CALL Dmalloc0(TrSrc, 1, nxy, myzbf, myzef)
CALL Dmalloc0(PrecSrc, 1, nxy, myzbf, myzef)

IF(nTracerCntl%lchidkgen .AND. nTracerCntl%lxslib ) THEN 
  CALL Dmalloc0(FmInfo%PrecSrcK, 1, nprec, 1, nfsr, myzb, myze)
ELSE
  CALL Dmalloc0(FmInfo%PrecSrc, 1, nfsr, myzb, myze)
END IF


!Conditional MOC
CALL DMALLOC0(TranInfo%PhiShape, 1, nxy, myzb, myze)
CALL DMALLOC0(TranInfo%PhiShape0, 1, nxy, myzb, myze)

ALLOCATE(TranInfo%coreBetak(nprec))
IF(TranCntl%lCorrector .OR. TranCntl%lGuess) THEN
  ALLOCATE(TranInfo%Prev_corePrec(nprec))
  ALLOCATE(TranInfo%Prev_coreBeta(nprec))
END IF

IF(TranCntl%lDynamicBen) THEN
  ALLOCATE(TranInfo%fuelTemp(nxy, myzb:myze))
  TranInfo%fuelTemp = TranInfo%InitTemp
END IF

END SUBROUTINE

SUBROUTINE AllocTranFxr(Core, FmInfo, GroupInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,    ONLY : CoreInfo_Type,       FmInfo_Type,      GroupInfo_Type,     &
                       PE_TYPE,                                                   &
                       FxrInfo_Type,        Pin_Type,         Cell_Type
USE CNTL,       ONLY : nTracerCntl_Type
USE ALLOCS
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)

INTEGER :: myzb, myze, nxy, nLocalFxr, nprec, ng
INTEGER :: FxrIdxSt
INTEGER :: ifxr, iz, ixy, icel
INTEGER :: i, j

LOGICAL :: lXsLib


Fxr => FmInfo%Fxr
Pin => Core%Pin; CellInfo => Core%CellInfo
myzb = PE%myzb; myze = PE%myze
nxy = Core%nxy
nprec = GroupInfo%nprec; ng = GroupInfo%ng
lXsLib = nTracerCntl%lXsLib
IF(.NOT. lXsLib) RETURN

DO iz = myzb, myze
  DO ixy = 1, nxy
    FxrIdxSt = Pin(ixy)%FxrIdxSt; icel = Pin(ixy)%Cell(iz)
    nlocalFxr = CellInfo(icel)%nFxr
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j - 1
      !IF(.NOT. Fxr(ifxr, iz)%lfuel) CYCLE
      CALL DMALLOC(Fxr(ifxr, iz)%beta, nprec)
      CALL DMALLOC(Fxr(ifxr, iz)%velo, ng)
      CALL DMALLOC(Fxr(ifxr, iz)%veloh, ng)
      CALL DMALLOC(Fxr(ifxr, iz)%chid, ng)
      CALL DMALLOC(Fxr(ifxr, iz)%chip, ng)
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE


SUBROUTINE InitTransient(Core, RayInfo, FmInfo, CmInfo, ThInfo, eigv, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,     FmInfo_Type,      CmInfo_Type,      TranInfo_Type,     &
                             ThInfo_Type,                                                              &
                             RayInfo_Type,      GroupInfo_Type,    TranCntl_Type,    PE_TYPE
USE CNTL,             ONLY : nTracerCntl_Type
USE TH_MOD,           ONLY : ThVar
USE BasicOperation,   ONLY : CP_CA,             CP_VA,            MULTI_CA
USE BenchXs,          ONLY : DnpLambdaBen,      NeutVeloBen,      ChidBen
USE TRAN_MOD,         ONLY : FluxNormTransient
USE CMFD_MOD,         ONLY : CMFDPsiUpdt, ConvertSubPlnPhi
USE MPIAxSolver_Mod,  ONLY : InitTranAxNSolver
USE TranAxNUtil_Mod,  ONLY : InitTranAxNUtil
USE TranMacXsLib_Mod, ONLY : InitChidLib,       InitLambdaLib
USE PromptFeedback_mod, ONLY : InitPromptFeedback
IMPLICIT NONE


TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE

REAL :: eigv

INTEGER :: nxy, nfsr, myzb, myze, myzbf, myzef
INTEGER :: nprec, ng, nfuelcell
INTEGER :: i, iz, ifxr, ixy, icel

REAL, POINTER :: Phis(:, :, :), Psi(:, :)
REAL, POINTER :: PhiC(:, :, :), PhiFM(:, :, :)
REAL, POINTER :: PsiC(:, :), PsiFm(:, :) 
REAL, POINTER :: TranPhi(:, :, :), TranPsi(:, :), TranPsid(:, :)
REAL, POINTER :: TranPhiCm(:, :, :), TranPsiCm(:, :), TranPsiCmd(:, :)
REAL, POINTER :: TranPhiFm(:, :, :), TranPsiFm(:, :), TranPsiFmd(:, :)

nxy = Core%nxy; nfsr = Core%nCoreFsr
myzb = PE%myzb; myze = PE%myze
myzbf = PE%myzbf; myzef = PE%myzef
ng = GroupInfo%ng; nprec = GroupInfo%nprec


CALL FluxNormTransient(Core, RayInfo, FmInfo, CmInfo, TranInfo, ng, nTracerCntl, PE)
IF(TranCntl%lDynamicBen) CALL InitPromptFeedback(Core, FmInfo, CmInfo, TranInfo, PE, ng)
CALL SetTimeStep(TranCntl)
IF(TranCntl%lNNSampling) CALL SetSamplingTimeStep(TranCntl)

Phis => FmInfo%Phis; Psi => FmInfo%Psi
TranPhi => FmInfo%TranPhi; 
TranPsi => FMInfo%TranPsi; TranPsid => FmInfo%TranPsid

PhiC => CmInfo%PhiC; PsiC => CmInfo%PsiC
TranPhiCm => CmInfo%TranPhiCm
TranPsiCm => CmInfo%TranPsiCM; TranPsiCmd => CmInfo%TranPsiCmd

PhiFm => CmInfo%PhiFM; PsiFm => CmInfo%PsiFm
TranPhiFm => CmInfo%TranPhiFm
TranPsiFm => CmInfo%TranPsiFm; TranPsiFmd => CmInfo%TranPsiFmd

!IF(PE%lCUDACMFD) THEN
!  CALL SetCMFDEnviorment(Core, CmInfo, ng, PE)
!  CALL ConvertSubPlnPhi(PhiC, PhiFm, 1)
!  CALL CMFDPsiUpdt(Phifm, Psifm)
!END IF

TranInfo%eigv0 = eigv
TranInfo%PowerLevel0 = nTracerCntl%PowerLevel
TranInfo%PowerLevel = nTracerCntl%PowerLevel
TranInfo%PowerLeveld = nTracerCntl%PowerLevel
!Eigen value Normalization 
CALL MULTI_CA(1._8/eigv, PsiC(1:nxy, myzb:myze), nxy, myze-myzb+1)
CALL MULTI_CA(1._8/eigv, Psi(1:nfsr, myzb:myze), nfsr, myze-myzb+1)
IF(nTracerCntl%lSubplane) THEN
  CALL MULTI_CA(1._8/eigv, PsiFm(1:nxy, myzbf:myzef), nxy, myzef-myzbf+1)
ENDIF

!CALL MULTI_VA()
!CALL CmFmInfoSync(Core, FmInfo, CmInfo, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)

CALL CP_VA(TranPhi(1:nFsr, myzb:myze, 1:ng), Phis(1:nFsr, myzb:myze, 1:ng), nFsr, myze-myzb+1, ng)
CALL CP_VA(TranPsi(1:nFsr, myzb:myze), Psi(1:nFsr, myzb:myze), nFsr, myze-myzb+1)
CALL CP_VA(TranPsid(1:nFsr, myzb:myze), Psi(1:nFsr, myzb:myze), nFsr, myze-myzb+1)

CALL CP_VA(FmInfo%TranPower(1:nFsr, myzb:myze), FmInfo%Power(1:nFsr, myzb:myze), nFsr, myze-myzb+1)

CALL CP_VA(TranPhiCm(1:nxy, myzb:myze, 1:ng), PhiC(1:nxy, myzb:myze, 1:ng), nxy, myze-myzb+1, ng)
CALL CP_VA(TranPsiCm(1:nxy, myzb:myze), PsiC(1:nxy, myzb:myze), nxy, myze-myzb+1)
CALL CP_VA(TranPsiCmd(1:nxy, myzb:myze), PsiC(1:nxy, myzb:myze), nxy, myze-myzb+1)

CALL CP_VA(TranPhiFm(1:nxy, myzb:myze, 1:ng), PhiFm(1:nxy, myzb:myze, 1:ng), nxy, myze-myzb+1, ng)
CALL CP_VA(TranPsiFm(1:nxy, myzb:myze), PsiFm(1:nxy, myzb:myze), nxy, myze-myzb+1)
CALL CP_VA(TranPsiFmd(1:nxy, myzb:myze), PsiFm(1:nxy, myzb:myze), nxy, myze-myzb+1)

!Exponential Transformation
CALL CP_CA(TranInfo%Expo(1:nxy, myzb:myze, 1:ng), 1._8, nxy, myze - myzb + 1, ng)
CALL CP_CA(TranInfo%Expo_Alpha(1:nxy, myzb:myze, 1:ng), 0._8, nxy, myze - myzb + 1, ng)
CALL CP_CA(TranInfo%FmExpo(1:nFsr, myzb:myze, 1:ng), 1._8, nFsr, myze - myzb + 1, ng)
CALL CP_CA(TranInfo%FmExpo_Alpha(1:nFsr, myzb:myze, 1:ng), 0._8, nFsr, myze - myzb + 1, ng)

nfuelcell = 0
!$OMP PARALLEL PRIVATE(icel)
DO iz = myzb, myze
  !$OMP DO SCHEDULE(GUIDED) REDUCTION(+:nfuelcell)
  DO ixy = 1, Core%nxy
    icel = Core%Pin(ixy)%Cell(iz)
    IF(Core%CellInfo(icel)%lFuel) THEN
      nfuelcell = nfuelcell + 1
    END IF
  END DO
  !$OMP END DO
END DO
!$OMP END PARALLEL
TranInfo%nfuelcell = nfuelcell

IF(nTracerCntl%l3dim) THEN
  CALL InitTranAxNSolver(nTracerCntl%AxSolver, Eigv, TranInfo, TranCntl, PE)
ENDIF

!Init Transient TH
IF(nTracerCntl%lFeedback) THEN
  CALL InitTransientTH(THInfo, ThVar, nxy, Core%nz)
ENDIF

IF(nTRACERCntl%lXsLib) THEN
  IF(.NOT. nTracerCntl%lchidgen) CALL InitChidLib(Core, FmInfo, GroupInfo, PE, nTracerCntl, TranInfo%Chid, ng)
  CALL InitLambdaLib(TranInfo%Lambda, nTracerCntl%refdcy_del, nTracerCntl%llibdcy_del, nTracerCntl%lfitbeta)
ELSE
  !Chid, Lambda, neut_velo
   CALL ChidBen(TranInfo%Chid)
   CALL DnpLambdaBen(TranInfo%Lambda)
   !CALL NeutVeloBen(TranInfo%Neut_Velo)
ENDIF

DO i = 1, nprec
  TranInfo%InvLambda(i) = 1._8 / TranInfo%Lambda(i)
ENDDO

DO iz = myzb, myze
  DO ifxr = 1, Core%nCoreFxr
    FmInfo%Fxr(ifxr, iz)%imix0 = FmInfo%Fxr(ifxr, iz)%imix
  ENDDO
ENDDO


NULLIFY(TranPhi, TranPsi, TranPsid)
NULLIFY(TranPhiCm, TranPsiCm, TranPsiCmd)
NULLIFY(TranPhiFm, TranPsiFm, TranPsiCmd)
NULLIFY(Phis, Psi)
NULLIFY(PhiC, PsiC); NULLIFY(PhiFm, PsiFm); 

END SUBROUTINE

SUBROUTINE SetUnitPowerLevel(Core, CmInfo, TranInfo, nTracerCntl, PE, ng)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,     CmInfo_Type,      TranInfo_Type,      &
                             PE_Type,           PinXS_Type
USE CNTL,             ONLY : nTracerCntl_Type
#ifdef MPI_ENV        
USE MPIComm_Mod,      ONLY : REDUCE
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
INTEGER :: ng

TYPE(PinXs_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: PhiC(:, :, :)
REAL, POINTER :: PinVol(:, :)
REAL :: volsum, pwsum, phisum
REAL :: PinPw, vol, AvgPw
INTEGER :: myzb, myze, nxy
INTEGER :: iz, ixy, ig
INTEGER :: COMM
REAL :: buf0(3), buf(3)

COMM = PE%MPI_CMFD_COMM

nxy = Core%nxy
myzb = PE%myzb; myze = PE%myze

PinXs => CmInfo%PinXs
PhiC => CmInfo%PhiC
PinVol => Core%PinVol

volsum = 0;
pwsum = 0;
phisum = 0;
DO iz = myzb, myze
  DO ixy = 1, nxy 
    pinpw = 0
    DO ig = 1, ng
      PinPw = PinPw + PinXs(ixy, iz)%xskf(ig) * PhiC(ixy, iz, ig) 
    ENDDO
    IF(PinPw .LE. 0._8) CYCLE
    vol = PinVol(ixy, iz)
    phisum = phisum + sum(PhiC(ixy, iz, 1:ng)) * vol
    PwSum = PwSum + PinPw * vol
    volsum = volsum + vol
  ENDDO
ENDDO

#ifdef MPI_ENV
buf0 = (/PwSum, PhiSum, volsum/)
CALL REDUCE(buf0, buf, 3, COMM, .TRUE.)
PwSum = Buf(1); PhiSum = Buf(2); Volsum = Buf(3)
#endif

AvgPw = PwSum / VolSum
TranInfo%UnitPowerLevel0 = nTRACERCntl%PowerLevel / AvgPw
TranInfo%PwSum0 = PwSum

NULLIFY(PinXS, PhiC, PinVol)

END SUBROUTINE

SUBROUTINE FluxNormTransient(Core, RayInfo, FmInfo, CmInfo, TranInfo, ng, nTracerCntl, PE)
! Flux level normalization
! normalize flux such that average flux in fuel region be unity
! then update fission source and moments accordingly
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,     FmInfo_Type,       CmInfo_Type,     &
                             RayInfo_Type,      TranInfo_Type,     PE_TYPE,         &
                             AxFlx_Type,        PinXs_Type
USE CNTL,             ONLY : nTracerCntl_Type
USE MPIAxSolver_Mod,  ONLY : NormMpiAxNVar,     UpdtTranSol2AxNvar
USE BasicOperation,   ONLY : MULTI_CA
#ifdef MPI_ENV        
USE MPIComm_Mod,      ONLY : REDUCE
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
INTEGER :: ng

TYPE(PinXs_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: PhiC(:, :, :), PsiC(:, :), PhiFm(:, :, :), PsiFm(:, :)
REAL, POINTER :: Phis(:, :, :), Psi(:, :), Jout(:, :, :, :, :)
REAL, POINTER :: PinVol(:, :)

INTEGER :: nbd, nfsr, nxy, myzb, myze, myzbf, myzef
INTEGER :: npol, nAzi, nSv
INTEGER :: ixy, iz, ig

REAL :: vol, volsum, pwsum, phisum, pinpw
REAL :: avgpw, avgflx, norm
REAL :: UnitPowerLevel0, Plevel

INTEGER :: comm
REAL :: buf0(3), buf(3)

COMM = PE%MPI_CMFD_COMM

nbd = 4
nxy = Core%nxy; nFsr = Core%nCoreFsr
myzb = PE%myzb; myze = PE%myze
myzbf = PE%myzbf; myzef =PE%myzef

npol = RayInfo%nPolarAngle; nSv = RayInfo%nPhiAngSv

PhiC => CmInfo%PhiC; PsiC => CmInfo%PsiC
PhiFm => CmInfo%PhiFm; PsiFm => CmInfo%PsiFm

Phis => FmInfo%Phis; Psi => FmInfo%Psi
Jout => CMInfo%RadJout
PinVol => Core%PinVol; PinXs => CmInfo%PinXs

volsum = 0;
pwsum = 0;
phisum = 0;
DO iz = myzb, myze
  DO ixy = 1, nxy 
    pinpw = 0
    DO ig = 1, ng
      PinPw = PinPw + PinXs(ixy, iz)%xskf(ig) * PhiC(ixy, iz, ig) 
    ENDDO
    IF(PinPw .LE. 0._8) CYCLE
    vol = PinVol(ixy, iz)
    phisum = phisum + sum(PhiC(ixy, iz, 1:ng)) * vol
    PwSum = PwSum + PinPw * vol
    volsum = volsum + vol
  ENDDO
ENDDO

#ifdef MPI_ENV
buf0 = (/PwSum, PhiSum, volsum/)
CALL REDUCE(buf0, buf, 3, COMM, .TRUE.)
PwSum = Buf(1); PhiSum = Buf(2); Volsum = Buf(3)
#endif

AvgPw = PwSum / VolSum
avgflx = PhiSum / Volsum

norm = ng / avgflx
norm = 1._8

IF(nTracerCntl%l3dim) THEN
  CALL UpdtTranSol2AxNvar(CmInfo, Core%PinVolFM, 1._8, PE)
ENDIF

CALL MULTI_CA(norm, Phis(1:nfsr, myzb:myze, 1:ng), nfsr, myze - myzb + 1, ng)
CALL MULTI_CA(norm, Psi(1:nfsr, myzb:myze), nfsr, myze - myzb + 1)


CALL MULTI_CA(norm, PhiC(1:nxy, myzb:myze, 1:ng), nxy, myze - myzb + 1, ng)
CALL MULTI_CA(norm, PsiC(1:nxy, myzb:myze), nxy, myze - myzb + 1)

IF(nTRACERCntl%lSubPlane) THEN
  CALL MULTI_CA(norm, PhiFm(1:nxy, myzbf:myzef, 1:ng), nxy, myzef - myzbf + 1, ng)
  CALL MULTI_CA(norm, PsiFm(1:nxy, myzbf:myzef), nxy, myzef - myzbf + 1)
ENDIF

DO ig = 1, ng
  CALL MULTI_CA(norm, FmInfo%PhiAngIn(1:npol, 1:nSv, myzb:myze, ig), nPol, nSv, myze-myzb+1)
  CALL MULTI_CA(norm, Jout(1:2, 1:nbd, 1:nxy, myzb:myze, ig), 2, nbd, nxy,  myze-myzb + 1)
ENDDO

IF(nTracerCntl%l3dim) THEN
  CALL NormMpiAxNVar(norm, nTracerCntl%AxSolver)
ENDIF

UnitPowerLevel0 = nTRACERCntl%PowerLevel / (AvgPw * norm)
TranInfo%UnitPowerLevel0 = UnitPowerLevel0
TranInfo%PwSum0 = PwSum
volsum = 0;
pwsum = 0;
phisum = 0;
DO iz = myzb, myze
  DO ixy = 1, nxy 
    pinpw = 0
    DO ig = 1, ng
      PinPw = PinPw + PinXs(ixy, iz)%xskf(ig) * PhiC(ixy, iz, ig) 
    ENDDO
    IF(PinPw .LE. 0._8) CYCLE
    vol = PinVol(ixy, iz)
    phisum = phisum + sum(PhiC(ixy, iz, 1:ng)) * vol
    PwSum = PwSum + PinPw * vol
    volsum = volsum + vol
  ENDDO
ENDDO

#ifdef MPI_ENV
buf0 = (/PwSum, PhiSum, volsum/)
CALL REDUCE(buf0, buf, 3, COMM, .TRUE.)
PwSum = Buf(1); PhiSum = Buf(2); Volsum = Buf(3)
#endif

AvgPw = PwSum / VolSum
Plevel = AvgPw * UnitPowerLevel0

NULLIFY(PhiC, PsiC)
NULLIFY(PhiFm, PsiFm)

END SUBROUTINE

SUBROUTINE SetTimeStep(TranCntl)
USE PARAM
USE TYPEDEF,       ONLY : TranCntl_Type
IMPLICIT NONE

TYPE(TranCntl_Type) :: TranCntl

REAL :: Delt0, Tbeg, Tend
INTEGER :: i, j

Tend = TranCntl%Tend

i = 0

TBeg = 0
DO j = 1, TranCntl%TStep_inp(0)
  Tend = TranCntl%Tstep_inp(j); Delt0 = TranCntl%Tdiv_inp(j)
  DO 
    i = i + 1
    Tbeg = Tbeg + Delt0
    TranCntl%T(i) = Tbeg
    TranCntl%DelT(i) = Delt0
    IF(abs(Tend - Tbeg) .LT. + 1.E-6_8) THEN
      TranCntl%T(i) = Tend
      EXIT
    END IF
  END DO
END DO

TranCntl%nstep = i

END SUBROUTINE

SUBROUTINE SetSamplingTimeStep(TranCntl)
  USE PARAM
  USE TYPEDEF,       ONLY : TranCntl_Type
  IMPLICIT NONE
  
  TYPE(TranCntl_Type) :: TranCntl
  
  REAL :: Delt0, Tbeg, Tend
  INTEGER :: i, j
  
  TBeg = TranCntl%Sbeg
  Tend = TranCntl%Send
  Delt0 = TranCntl%Speriod
  
  TranCntl%nSstep = CEILING((Tend-Tbeg)/Delt0)+1
  ALLOCATE(TranCntl%Ssteps(TranCntl%nSstep))
  TranCntl%Ssteps(1) = Tbeg
  DO i = 2, TranCntl%nSstep
    Tbeg = Tbeg + Delt0
    TranCntl%Ssteps(i) = Tbeg
    IF(abs(Tend - Tbeg) .LT. + 1.E-6_8) THEN
      TranCntl%Ssteps(i) = Tend
      EXIT
    END IF
  END DO
  
  END SUBROUTINE

