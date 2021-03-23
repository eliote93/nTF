#include <defines.h>
  
SUBROUTINE UpdtPowerLevel(Core, FmInfo, CmInfo, TranInfo, ng, nTracerCntl, PE)
! Flux level normalization
! normalize flux such that average flux in fuel region be unity
! then update fission source and moments accordingly
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,     FmInfo_Type,       CmInfo_Type,     &
                           TranInfo_Type,     PE_TYPE,                            &
                           PinXs_Type
USE CNTL,           ONLY : nTracerCntl_Type
USE BasicOperation, ONLY : MULTI_CA
#ifdef MPI_ENV
USE MPIComm_Mod,    ONLY : REDUCE
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
INTEGER :: ng

TYPE(PinXs_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: PhiC(:, :, :), PsiC(:, :)
REAL, POINTER :: PinVol(:, :)

INTEGER :: nbd, nfsr, nxy, myzb, myze, myzbf, myzef
INTEGER :: ixy, iz, ig

REAL :: AvgPw, PwSum, VolSum, PinPw, vol
REAL :: UnitPowerLevel0, plevel

INTEGER :: comm
REAL :: buf0(2), buf(2)

COMM = PE%MPI_CMFD_COMM

nbd = 4
nxy = Core%nxy; nFsr = Core%nCoreFsr
myzb = PE%myzb; myze = PE%myze
myzbf = PE%myzbf; myzef =PE%myzef

PhiC => CmInfo%PhiC; PsiC => CmInfo%PsiC
PinVol => Core%PinVol; PinXs => CmInfo%PinXs

UnitPowerLevel0 = TranInfo%UnitPowerLevel0
volsum = 0;
pwsum = 0;
!phisum = 0;
DO iz = myzb, myze
  DO ixy = 1, nxy 
    pinpw = 0
    DO ig = 1, ng
      PinPw = PinPw + PinXs(ixy, iz)%xskf(ig) * PhiC(ixy, iz, ig) 
    ENDDO
    !IF(PinPw .LT. 0._8) CYCLE
    IF(PinPw .LE. 0._8) CYCLE
    vol = PinVol(ixy, iz)
    !phisum = phisum + sum(PhiC(ixy, iz, 1:ng)) * vol
    PwSum = PwSum + PinPw * vol
    volsum = volsum + vol
  ENDDO
ENDDO
#ifdef MPI_ENV
buf0 = (/PwSum, Volsum/)
CALL REDUCE(buf0, buf, 2, COMM, TRUE)
PWSUM = buf(1); Volsum = buf(2)
AvgPw = PwSum / VolSum
#endif

Plevel = AvgPw * UnitPowerLevel0
TranInfo%PowerLevel = Plevel
NULLIFY(PhiC, PsiC, PinVol, PinXS)
END SUBROUTINE
  
SUBROUTINE SaveTranSol(Core, FmInfo, CmInfo, TranInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,             ONLY : CoreInfo_Type,          FmInfo_Type,        CmInfo_Type,       &
                                TranInfo_Type,          TranCntl_Type,      PE_Type,           &
                                ThInfo_Type,            GroupInfo_Type
USE CNTL,                ONLY : nTracerCntl_Type
USE MPIAxSolver_Mod,     ONLY : SaveTranAxNSol,            UpdtTranSol2AxNvar
USE TH_MOD,              ONLY : SaveTranThSol
USE BasicOperation,      ONLY : CP_VA
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(THInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

INTEGER :: nxy, nfxr, nfsr,ng
INTEGER :: myzb, myze, myzbf, myzef, npln, nplnf

ng = GroupInfo%ng
nxy = Core%nxy; nfxr = Core%nCoreFxr; nfsr = Core%nCoreFsr
myzb = PE%myzb; myze = PE%myze; myzbf = PE%myzbf; myzef = PE%myzef
npln = myze-myzb+1; nplnf = myzef - myzbf + 1
!Moc Variables
CALL CP_VA(FmInfo%TranPhi(1:nfsr, myzb:myze, 1:ng), FmInfo%Phis(1:nfsr, myzb:myze, 1:ng), nfsr, npln, ng)
CALL CP_VA(FmInfo%TranPsid(1:nfsr, myzb:myze), FmInfo%TranPsi(1:nfsr, myzb:myze), nfsr, npln)
CALL CP_VA(FmInfo%TranPsi(1:nfsr, myzb:myze), FmInfo%Psi(1:nfsr, myzb:myze), nfsr, npln)
CALL CP_VA(FmInfo%TranPower(1:nfsr, myzb:myze), FmInfo%Power(1:nfsr, myzb:myze), nfsr, npln)

!Cmfd Variables
CALL CP_VA(CmInfo%TranPhiCm(1:nxy, myzb:myze, 1:ng), CmInfo%PhiC(1:nxy, myzb:myze, 1:ng), nxy, npln, ng)
CALL CP_VA(CmInfo%TranPsiCmd(1:nxy, myzb:myze), CmInfo%TranPsiCm(1:nxy, myzb:myze), nxy, npln)
CALL CP_VA(CmInfo%TranPsiCm(1:nxy, myzb:myze), CmInfo%PsiC(1:nxy, myzb:myze), nxy, npln)

IF(nTracerCntl%lSubPlane) THEN
  CALL CP_VA(CmInfo%TranPhiFm(1:nxy, myzb:myze, 1:ng), CmInfo%PhiFm(1:nxy, myzb:myze, 1:ng), nxy, npln, ng)
  CALL CP_VA(CmInfo%TranPsiFmd(1:nxy, myzb:myze), CmInfo%TranPsiFm(1:nxy, myzb:myze), nxy, npln)
  CALL CP_VA(CmInfo%TranPsiFm(1:nxy, myzb:myze), CmInfo%PsiFm(1:nxy, myzb:myze), nxy, npln)
ENDIF

IF(nTracerCntl%l3dim) THEN
  CALL UpdtTranSol2AxNvar(CmInfo, Core%PinVolFM, TranInfo%eigv0, PE)  
  CALL SaveTranAxNSol(TranInfo%eigv0, PE)
ENDIF

IF(nTracerCntl%lFeedback) THEN
  CALL SaveTranTHsol(Core, ThInfo, TranCntl, nTracerCntl, PE)
ENDIF

END SUBROUTINE

!SUBROUTINE XsPerturbation(TranInfo, TranCntl, nTracerCntl)
!USE PARAM
!USE TYPEDEF,             ONLY : TranInfo_Type,          TranCntl_Type,           &
!                                XsChange_Type
!USE CNTL,                ONLY : nTracerCntl_Type
!USE BenchXs,             ONLY : XsBenChange,            MacXsBen
!IMPLICIT NONE
!TYPE(TranInfo_Type) :: TranInfo
!TYPE(TranCntl_Type) :: TranCntl
!TYPE(nTracerCntl_Type) :: nTracerCntl
!
!TYPE(XsChange_Type), POINTER :: XsChange(:)
!
!REAL :: T, Tprev0, Tprev, DelT, Tbeg, Tend, w
!INTEGER :: NowStep, nChange
!INTEGER :: i, iso1, iso0
!XsChange => TranCntl%XsChange
!
!nChange = TranCntl%nChange
!NowStep = TranCntl%NowStep;
!T= TranCntl%T(NowStep)
!Tprev = 0
!IF(NowStep .GT. 1) Tprev = TranCntl%T(NowStep - 1)
!
!DO i = 1, nChange
!  Tbeg = XsChange(i)%Tbeg; Tend = XsChange(i)%Tend
!  IF(Tend .GT. T) CYCLE
!  IF(XsChange(i)%lComplete) CYCLE
!  XsChange(i)%lComplete = .TRUE.
!  iso0 = XsChange(i)%Iso0; iso1 = XsChange(i)%iso1
!  CALL XsBenChange(iso0, iso1, 1._8)
!ENDDO
!
!DO i = 1, nChange
!  Tbeg = XsChange(i)%Tbeg; Tend = XsChange(i)%Tend
!  IF(T .LT. Tbeg .OR. T .GT. Tend) CYCLE
!  Tprev0 = Max(Tprev, Tbeg)
!  w = (T-Tprev0) / (Tend - Tprev0)
!  iso0 = XsChange(i)%Iso0; iso1 = XsChange(i)%iso1
!  CALL XsBenChange(iso0, iso1, w)
!ENDDO
!
!END SUBROUTINE

SUBROUTINE UpdtResSrc(Core, FmInfo, CmInfo, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,             ONLY : CoreInfo_Type,          FmInfo_Type,         CmInfo_Type,       &
                                TranInfo_Type,          GroupInfo_Type,      TranCntl_Type,     &
                                PE_Type,                                                        &
                                AxFlx_Type,             Pin_Type,            PinXS_Type
USE CNTL,                ONLY : nTracerCntl_Type
USE BasicOperation,      ONLY : CP_CA
USE CMFD_MOD,            ONLY : Src
USE TRANCMFD_MOD,        ONLY : PrecSrc,                TrSrc,                                  &
                                HomKineticParamGen,     SetCmfdPrecParam
#ifdef MPI_ENV
USE MPIComm_MOD, ONLY : REDUCE, GetNeighDat
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCNtl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE

TYPE(PinXS_Type), POINTER :: PinXs(:, :)
TYPE(Pin_Type), POINTER :: Pin(:)
REAL, POINTER :: Psi(:, :), TranPsi(:, :), TranPsiD(:, :), Chid(:)
REAL, POINTER :: Prec(:, :, :), CellOmegam(:, :, :), CellOmega0(:, :, :), CellOmegap(:, :, :)
REAL, POINTER :: PinVolFm(:, :)

REAL, POINTER :: ResSrc(:, :, :)
REAL, POINTER :: TranPhi(:, :, :), Phi(:, :, :)
REAL, POINTER :: AxDtil(:, :, :, :), AxDhat(:, :, :, :)
INTEGER :: myzb, myze, myzbf, myzef
INTEGER :: NowStep, nxy, ng, nprec, nbd
INTEGER :: ig, ig2, ixy, iz, iz0, igb, ige, iprec, ineigh
INTEGER :: i, j, k 

REAL :: rvdt, chieff0, vol, area, SS
REAL :: Dtil, Dhat, Jnet, lkg
REAL :: DelT, kappa(100), Lambda(100)

nbd = 4
nxy = PE%nxy; ng = GroupInfo%ng
myzb = PE%myzb; myze = PE%myze
myzbf = PE%myzbf; myzef = PE%myzef
nPrec = GroupInfo%nPrec

NowStep = TranCntl%NowStep
Delt = TranCntl%DelT(NowStep)

Prec => CmInfo%PrecFm; Chid => TranInfo%Chid
Psi => CmInfo%PsiFm; TranPsi => CmInfo%TranPsiFm; TranPsiD => CmInfo%TranPsiFmd
CellOmegam => TranInfo%CellOmegam; CellOmegap => TranInfo%CellOmegap;
CellOmega0 => TranInfo%CellOmega0
PinVolFm => Core%PinVolFm
AxDtil => CmInfo%AxDtil;
AxDhat => CmInfo%AxDhat
Phi => CmInfo%PhiFm; TranPhi => CmInfo%TranPhiFm
ResSrc => CmInfo%ResSrcFm

Pin => Core%Pin; PinXs => CmInfo%PinXs
!Calculate homgenized chid and beta
CALL HomKineticParamGen(Core, FmInfo, CmInfo, TranInfo, GroupInfo, nTracerCntl, PE)
!Calculate the Precursor Quadratic approx. Coefficient
CALL SetCmfdPrecParam(Core, CmInfo, TranInfo, TranCntl, nTracerCntl, PE)

DO iprec = 1, nprec
  Lambda(iprec) = TranInfo%Lambda(iprec)
  kappa(iprec) = exp(-delt * lambda(iprec))
ENDDO

DO iz = myzbf, myzef
  iz0 = Core%SubPlaneMap(iz)
  DO ixy = 1, nxy
    PrecSrc(ixy, iz) = 0
    DO iprec = 1, nprec
      PrecSrc(ixy, iz) = PrecSrc(ixy, iz) + Lambda(iprec) * kappa(iprec) * Prec(Iprec, ixy, iz)
    ENDDO
    PrecSrc(ixy, iz) = PrecSrc(ixy, iz) + CellOmegam(0, ixy, iz0) * TranPsid(ixy, iz) +  CellOmega0(0, ixy, iz0) * TranPsi(ixy, iz)
  ENDDO
ENDDO

#ifdef MPI_ENV
DO ig = 1, ng
  IF(PE%nCMFDproc .GT. 1) THEN 
    CALL GetNeighDat(Phi(1:nxy, myzbf-1:myzef+1, ig), nxy, myzbf, myzef,   &
                           PE%myCmfdRank, PE%nproc, PE%MPI_CMFD_COMM)
  ENDIF
ENDDO
#endif


DO ig = 1, ng
  !Set Source
  DO iz = myzbf, myzef
    iz0 = Core%SubPlaneMap(iz)
    DO ixy = 1, nxy
      vol = PinVolFm(ixy, iz0)
      rvdt = 1._8 / (Delt * PinXs(ixy, iz)%velo(ig))
      TrSrc(ixy, iz) = chid(ig) * PrecSrc(ixy, iz)
      chieff0 = chid(ig) * (PinXs(ixy, iz0)%omega - PinXs(ixy, iz0)%betat)
      TrSrc(ixy, iz) = TrSrc(ixy, iz) + chieff0 *Psi(ixy, iz)
      TrSrc(ixy, iz) = TrSrc(ixy, iz) !+ rvdt * (TranPhi(ixy, iz, ig) - Phi(ixy, iz, ig))*Vol    
    ENDDO    
   
    DO ixy = 1, nxy
      vol = PinVolFm(ixy, iz0)
      SRC(ixy, iz) = Psi(ixy, iz) * PinXs(ixy, iz0)%chi(ig)
      igb = PinXs(ixy, iz0)%Xss(ig)%ib; ige = PinXs(ixy, iz0)%Xss(ig)%ie;
      SS = 0
      DO ig2 = igb, ige
        SS = SS + Phi(ixy, iz, ig2) * PinXs(ixy, iz0)%Xss(ig)%From(ig2)
      ENDDO
      SRC(ixy, iz) = SRC(ixy, iz) + SS * PinVolFm(ixy, iz0)
      SRC(ixy, iz)  = SRC(ixy, iz) + TrSrc(ixy, iz)
    ENDDO

    !LEAKAGE Part
    DO ixy = 1, nxy
      vol = PinVolFm(ixy, iz0)
      lkg = Vol * PinXs(ixy, iz0)%xsr(ig) * Phi(ixy, iz, ig)
      DO i = 1, nbd
        dtil = PinXs(ixy, iz0)%Dtil(i, ig)
        dhat = PinXs(ixy, iz0)%Dhat(i, ig)
        Jnet = (Dtil - Dhat) * Phi(ixy, iz, ig)
        ineigh = Pin(ixy)%NeighIdx(i)
        IF(ineigh .GT. 0) Jnet = Jnet -(dtil + dhat)  * Phi(ineigh, iz, ig)
        Jnet = Jnet * Core%hzfm(iz)
        lkg = lkg + Jnet
      ENDDO
      !
      IF(nTracerCntl%l3dim) THEN
        Area = PinVolFm(ixy, iz0)/Core%hzfm(iz)
        dtil = AxDtil(1, ixy, iz, ig)
        dhat = AxDhat(1, ixy, iz, ig)
        Jnet = (Dtil - Dhat) * Phi(ixy, iz, ig) - (Dtil + Dhat) * Phi(ixy, iz-1, ig)
        Jnet = Jnet * Area
        lkg = lkg + Jnet

        dtil = AxDtil(2, ixy, iz, ig)
        dhat = AxDhat(2, ixy, iz, ig)
        Jnet = (Dtil - Dhat) * Phi(ixy, iz, ig) - (Dtil + Dhat) * Phi(ixy, iz+1, ig)
        Jnet = Jnet * Area
        lkg = lkg + Jnet  
      ENDIF
      ResSrc(ixy, iz, ig) = Src(ixy, iz) - lkg
    ENDDO
  
  ENDDO
  

ENDDO

NULLIFY(Prec, Chid)
NULLIFY(Psi, TranPsi, TranPsiD)
NULLIFY(CellOmegam, CellOmegaP, CellOmega0)
NULLIFY(Phi, TranPhi, ResSrc)
NULLIFY(Pin, PinXs, AxDhat, AxDtil)

END SUBROUTINE

SUBROUTINE CmFmInfoSync(Core, FmInfo, CmInfo, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,             ONLY : CoreInfo_Type,          FmInfo_Type,           CmInfo_Type,         &
                                TranInfo_Type,          GroupInfo_Type,        TranCntl_Type,       &
                                PE_Type,                                                            &
                                FxrInfo_Type,           Pin_Type,               Cell_Type
USE CNTL,                ONLY : nTracerCntl_Type
USE MOC_MOD,             ONLY : PsiUpdate
USE BasicOperation,      ONLY : MULTI_CA
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)


REAL, POINTER :: Psi(:, :), Phis(:, :, :), PhiC(:, :, :), PsiC(:, :) 

INTEGER :: myzb, myze, myzbf, myzef, nxy, nfxr, nfsr, nprec, ng
INTEGER :: FxrIdxSt, FsrIdxst, nLocalFxr, nFsrInFxr, NowStep
INTEGER :: ixy, iz, iz0, icel, ifxr, ifsr, iprec, ig
INTEGER :: i, j

REAL :: DelT, reigv, vol
REAL, POINTER :: PinVol(:, :)

myzb = PE%myzb; myze = PE%myze
myzbf = PE%myzbf; myzef = PE%myzef
nxy = Core%nxy; nfxr = Core%nCoreFxr; nfsr = Core%nCoreFsr
ng = GroupInfo%ng

Fxr => FmInfo%Fxr
Pin => Core%Pin; CellInfo => Core%CellInfo
nprec = TranInfo%nprec
Phis => FmInfo%Phis; Psi => FmInfo%Psi
PhiC => CMInfo%PhiC; PsiC => CmInfo%PsiC
PinVol => COre%PinVOl
CALL PsiUpdate(Core, Fxr, phis, psi, myzb, myze, ng, nTracerCNtl%lxslib, GroupInfo)      
CALL MULTI_CA(1._8 /TranInfo%eigv0, FmInfo%psi(1:nfsr, myzb:myze), nFsr, myze - myzb +1)

DO iz = myzb, myze
  DO ixy = 1, nxy
    FxrIdxSt = Pin(ixy)%FxrIdxSt; FsrIdxSt = Pin(ixy)%FsrIdxSt
    icel = Pin(ixy)%Cell(iz); nLocalFxr = CellInfo(icel)%nFxr
    PhiC(ixy, iz, 1:ng) = 0
    PsiC(ixy, iz) =0
    DO j = 1, CellInfo(icel)%nFsr
      ifsr = FsrIdxSt + j - 1
      vol = CellInfo(icel)%vol(j)*core%hz(iz)
      DO ig = 1, ng
        PhiC(ixy, iz, ig) = PhiC(ixy, iz, ig) + PhiS(ifsr, iz, ig) * vol
      ENDDO
      PsiC(ixy, iz) = PsiC(ixy, iz) + Psi(ifsr, iz) * vol
    ENDDO   
    DO ig =1, ng
      PhiC(ixy, iz, ig) = PhiC(ixy, iz, ig) / PinVOl(ixy,iz)
    ENDDO
    
  ENDDO
ENDDO
!
!Fxr => FmInfo%Fxr
!Pin => Core%Pin; CellInfo => Core%CellInfo
!nprec = TranInfo%nprec
!Phis => FmInfo%Phis; Psi => FmInfo%Psi
!PhiC => CMInfo%PhiC; PsiC => CmInfo%PsiC
!PinVol => COre%PinVOl
NULLIFY(Fxr, Pin, CellInfo, PinVol)
NULLIFY(Phis, Psi, PhiC, PsiC)
END SUBROUTINE

SUBROUTINE ChkOutputWrite(Flag, TranCntl)
USE PARAM
USE TYPEDEF,  ONLY : TranCntl_Type
TYPE(TranCntl_Type) :: TranCntl
LOGICAL :: Flag

INTEGER :: NowWriteStep
Flag = .FALSE.

IF(TranCntl%NowWriteStep .EQ. 0) TranCntl%NowWriteStep = 1
NowWriteStep =  TranCntl%NowWriteStep
IF(NowWriteStep .GT. TranCntl%NTWriteOut) RETURN
T = TranCntl%T(TranCntl%NowStep)
IF((T -TranCntl%TWriteOut(NowWriteStep)) .GT. -1.e-5) THEN
  Flag = .TRUE.
ENDIF
END SUBROUTINE

SUBROUTINE UpdtExpTrsf(Core, FmInfo, CmInfo, TranInfo, GroupInfo, lupdt, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,                 ONLY : CoreInfo_Type,          FmInfo_Type,            CmInfo_Type,       &
                                    TranInfo_Type,          TranCntl_Type,          PE_Type,           &
                                    GroupInfo_Type,                                                    &
                                    Cell_Type,              Pin_Type
USE CNTL,                    ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
LOGICAL :: lUpdt

REAL, POINTER :: Phi(:, :, :), PHIC(:, :, :), TranPhiCm(:, :, :)
REAL, POINTER :: Expo_Alpha(:, :, :), Expo(:, :, :)
REAL, POINTER :: FmExpo_Alpha(:, :, :), FmExpo(:, :, :)
REAL, POINTER :: TranPhi(:, :, :)



TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)

INTEGER :: ng, nxy, nFsr,  myzb, myze
INTEGER :: FsrIdxSt, nLocalFsr
INTEGER :: iz, ixy, ifsr, icel, ig, j
INTEGER :: nowstep

REAL :: Delt, Deltn
REAL :: ratio
REAL :: maxratio
REAL :: phisum, vol

PhiC=> CmInfo%PhiC; TranPhiCm => CmInfo%TranPhiCm
Phi => FmInfo%Phis; TranPhi => FmInfo%TranPhi
Pin => Core%Pin; CellInfo => Core%CellInfo

Expo_Alpha => TranInfo%Expo_Alpha; Expo => TranInfO%Expo
FmExpo_Alpha => TranInfo%FmExpo_Alpha; FmExpo => TranInfO%FmExpo
ng = GroupInfo%ng
nxy = Core%nxy; nFsr = Core%nCoreFsr
myzb = PE%myzb; myze = PE%myze

nowstep = TranCntl%NowStep
Delt = TranCntl%DelT(nowstep); Deltn = Delt
IF(.NOT. lupdt .AND. NowStep .NE. TranCntl%nStep) Deltn = TranCntl%DelT(nowstep+1)

DO ig = 1, ng
  DO iz = myzb, myze
    DO ixy = 1, nxy
      Expo_Alpha(ixy, iz, ig) = 0
      ratio = PhiC(ixy, iz, ig)/TranPhiCm(ixy, iz, ig)
      IF(PHIC(ixy, iz, ig) .GT. 0 .AND. TranPhiCm(ixy, iz, ig) .GT. 0) THEN
        !ratio = min(ratio, 10._8)
        Expo_Alpha(ixy, iz, ig) = log(abs(ratio)) / delt
        Expo(ixy, iz, ig) = Exp(Expo_Alpha(ixy, iz, ig) * Deltn)
      ENDIF
      FsrIdxSt = Pin(ixy)%FsrIdxSt; icel = Pin(ixy)%Cell(iz)
      nLocalFsr = CellInfo(icel)%nFsr
      DO j = 1, nLocalFsr
        ifsr = FsrIdxSt + j - 1
        FmExpo_Alpha(ifsr, iz, ig) = Expo_Alpha(ixy, iz, ig)
        FmExpo(ifsr, iz, ig) =  Expo(ixy, iz, ig)
      ENDDO
    ENDDO  
  ENDDO
ENDDO


NULLIFY(Pin, CellInfo)
NULLIFY(Expo_Alpha, Expo)
NULLIFY(PhiC, TranPhiCm, Phi, TranPhi)

END SUBROUTINE

SUBROUTINE AdaptiveTheta(TranInfo, TranCntl)
USE PARAM
USE TYPEDEF,         ONLY : TranInfo_Type,          TranCntl_Type
IMPLICIT NONE
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCNtl_Type) :: TranCntl
REAL :: rtdblr, rtdblrd

rtdblrd = TranInfo%rtdblr
rtdblr = log(TranInfo%PowerLevel /TranInfo%PowerLeveld)/ log(2._8)
IF(rtdblr .GT. 0._8) THEN
  TranCntl%Theta = 0.5
ELSE
  IF(rtdblr .GT. rtdblrd) THEN
    IF(rtdblr .GT. -0.1) THEN
      TranCntl%Theta = 1._8
    ELSE
      TranCntl%theta = 0.5 + 1._8/3._8  
    ENDIF
  ELSE
    IF(rtdblr .GT. -0.2) THEN
      TranCntl%Theta = 0.5 + 1._8/3._8
    ELSE
      TranCntl%theta = 0.5
    ENDIF    
  ENDIF
ENDIF

TranInfo%rtdblr = rtdblr
END SUBROUTiNE


SUBROUTINE SolExpExtpltn(Core, FmInfo, CmInfo, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,                 ONLY : CoreInfo_Type,          FmInfo_Type,            CmInfo_Type,       &
                                    TranInfo_Type,          TranCntl_Type,          PE_Type,           &
                                    GroupInfo_Type,                                                    &
                                    Cell_Type,              Pin_Type
USE CNTL,                    ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

REAL, POINTER :: Phis(:, :, :), PHIC(:, :, :), TranPhiCm(:, :, :)
REAL, POINTER :: Expo_Alpha(:, :, :), Expo(:, :, :)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)

INTEGER :: ng, nxy, nFsr,  myzb, myze
INTEGER :: FsrIdxSt, nLocalFsr
INTEGER :: iz, ixy, ifsr, icel, ig, j
INTEGER :: nowstep

REAL :: Delt

PhiC=> CmInfo%PhiC; TranPhiCm => CmInfo%TranPhiCm; Phis => FmInfo%Phis
Pin => Core%Pin; CellInfo => Core%CellInfo

Expo_Alpha => TranInfo%Expo_Alpha; Expo => TranInfO%Expo

ng = GroupInfo%ng
nxy = Core%nxy; nFsr = Core%nCoreFsr
myzb = PE%myzb; myze = PE%myze

nowstep = TranCntl%NowStep
Delt = TranCntl%DelT(nowstep)
DO ig = 1, ng
  DO iz = myzb, myze
    DO ixy = 1, nxy
      FsrIdxSt = Pin(ixy)%FsrIdxSt; icel = Pin(ixy)%Cell(iz)
      nLocalFsr = CellInfo(icel)%nFsr
      !PhiC(ixy, iz, ig) = Expo(ixy, iz, ig) * PhiC(ixy, iz, ig)
      DO j = 1, nLocalFsr
        ifsr = FsrIdxSt + j - 1
        !Phis(ifsr, iz, ig) = Expo(ixy, iz, ig) * Phis(ifsr, iz, ig)
      ENDDO
    ENDDO
  ENDDO
ENDDO

NULLIFY(Pin, CellInfo)
NULLIFY(Expo_Alpha, Expo)
NULLIFY(PhiC, TranPhiCm, Phis)
!
END SUBROUTINE

SUBROUTINE CalcFisRate(FisRate, CmInfo, nxy, myzb, myze, ng)
USE TYPEDEF,          ONLY : CmInfo_Type,       PinXs_Type
IMPLICIT NONE
REAL :: FisRate(nxy, myzb:myze)
TYPE(CmInfo_Type) :: CmInfo
INTEGER :: nxy, myzb, myze, ng

TYPE(PinXs_Type), POINTER :: PinXs(:,:)
REAL, POINTER :: PhiC(:,:,:)
REAL :: pinfr
INTEGER :: ixy, iz, ig

PinXS => CmInfo%PinXS
PhiC => CmInfo%PhiC

!$OMP PARALLEL PRIVATE(pinfr)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO iz = myzb, myze
  DO ixy = 1, nxy
    pinfr = 0.
    DO ig = 1, ng
      pinfr = pinfr + PinXS(ixy, iz)%xskf(ig) * PhiC(ixy, iz, ig) 
    END DO 
    FisRate(ixy, iz) = pinfr
  END DO
END DO 
!$OMP END DO
!$OMP END PARALLEL

NULLIFY(PinXS, PhiC)

END SUBROUTINE
  

