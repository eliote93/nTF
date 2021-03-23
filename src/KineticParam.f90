#include <defines.h>
#define oneoverv
SUBROUTINE HomKineticParamGen(Core, FmInfo, CmInfo, TranInfo, GroupInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,    FmInfo_Type,      CmInfo_Type,          & 
                             TranInfo_Type,    GroupInfo_Type,   PE_Type,              &
                             FXRInfo_Type,     PinXs_Type,       Pin_Type,             &              
                             PinInfo_Type,     Cell_Type
USE CNTL,             ONLY : nTracerCntl_Type
                             

USE CMFD_mod,       ONLY : XsMac
USE BenchXs,        ONLY : XsBaseBen,        DnpBetaBen,        NeutVeloBen,           &
                           DnpBetaDynBen,    NeutVeloDynBen
USE MacXsLib_Mod,   ONLY : MacXsBase

USE BasicOperation,   ONLY : CP_VA,            CP_CA,            MULTI_VA
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) ::  GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

!
TYPE(FxrInfo_Type), POINTER :: FXR(:, :)
TYPE(PinXs_Type), POINTER :: PinXs(:, :)
REAL, POINTER:: Psi(:, :), phis(:, :, :)

!
TYPE(FxrInfo_Type), POINTER :: myFXR
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
!
INTEGER :: myzb, myze, ng, nprec
INTEGER :: nFxr, nFsr, nxy
INTEGER :: nCellType, nPinType, nlocalFxr, nlocalFsr, nFsrInFxr
INTEGER :: icel, ipin, iz, ixy ,ifxr, ifsr, ig, itype, ifsrlocal
INTEGER :: FsrIdxSt, FxrIdxSt
INTEGER :: i, j, k, l, m
INTEGER :: iResoGrpBeg, iResoGrpEnd, norg, nChi

REAL :: Beta0(1:50), Beta(1:50), chid(500)
REAL :: phisum(500), velo(500), velo0(500)
REAL :: psisum, vol, volsum, betat

LOGICAL :: lXsLib


lXsLib = nTracerCntl%lXsLib

IF(lxslib) THEN
  norg = GroupInfo%norg; nChi = GroupInfo%nChi
  iResoGrpBeg = GroupInfo%nofg + 1
  iResoGrpEnd = GroupInfo%nofg + GroupInfo%norg 
ENDIF

nxy = Core%nxy; nFsr = Core%nCoreFsr; nFxr = Core%nCoreFxr
myzb = PE%myzb; myze = PE%myze

ng = GroupInfo%ng; nprec = GroupInfo%nprec

Chid(1:ng) = TranInfo%Chid(1:ng)

Pin => Core%Pin; CellInfo => Core%CellInfo
Fxr => FmInfo%Fxr; Psi => FmInfo%Psi; phis => FmInfo%Phis
PinXs => CmInfo%PinXs
DO iz = myzb, myze
  DO ixy = 1, nxy
    beta = 0; psisum = 0
    velo(1:ng) = 0; phisum(1:ng) = 0
    FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
    icel = Pin(ixy)%Cell(iz)
    nlocalFxr = CellInfo(icel)%nFxr; nlocalFsr = CellInfo(icel)%nFsr
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1; nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
      myFxr => Fxr(ifxr, iz)
      IF(lXsLib) THEN
        beta0(1:nprec) = myFxr%Beta(1:nprec)
        velo0(1:ng) = myFxr%veloh(1:ng)
      ELSE
        IF(nTracerCntl%lDynamicBen) THEN
          CALL DnpBetaDynBen(myFxr%imix, TranInfo%fuelTemp(ixy, iz), beta0(1:nprec))
          CALL NeutVeloDynBen(myFxr%imix, TranInfo%fuelTemp(ixy, iz), velo0(1:ng))
        ELSE
          CALL DnpBetaBen(myFxr%imix, beta0(1:nprec))
          CALL NeutVeloBen(myFxr%imix, velo0(1:ng))
        END IF
      ENDIF
      DO i = 1, nFsrInFxr
        ifsrlocal = Cellinfo(icel)%MapFxr2FsrIdx(i, j)
        vol = CellInfo(icel)%vol(ifsrlocal)
        ifsr = FsrIdxSt + ifsrlocal - 1
        beta(1:nprec) = beta(1:nprec) + beta0(1:nprec) * Psi(ifsr, iz) * vol
        psisum = psisum + Psi(ifsr, iz) * vol
        DO ig = 1, ng
#ifdef oneoverv
          velo(ig) = velo(ig) + 1._8/velo0(ig) * phis(ifsr, iz, ig) * vol
#else
          velo(ig) = velo(ig) + velo0(ig) * phis(ifsr, iz, ig) * vol
#endif
          phisum(ig) = phisum(ig) + phis(ifsr, iz, ig) * vol
        ENDDO
      ENDDO
    ENDDO
    IF(psisum .GT. 0) THEN
      beta(1:nprec) = beta(1:nprec) / psisum
    ELSE 
      beta(1:nprec) = 0
    ENDIF
    PinXs(ixy, iz)%Beta(1:nprec) = beta(1:nprec)
    DO ig = 1, ng
      velo(ig) = velo(ig) / phisum(ig)
#ifdef oneoverv
      velo(ig) = 1._8 / velo(ig) 
#endif
      PinXs(ixy, iz)%velo(ig) = velo(ig)
    ENDDO
    
    betat = sum(beta(1:nprec))
    PinXs(ixy, iz)%betat = betat
    PinXs(ixy, iz)%Chip(1:ng) = (PinXs(ixy, iz)%Chi(1:ng) - betat * Chid(1:ng)) / ( 1._8 - betat)
  ENDDO
ENDDO

NULLIFY(Pin, CellInfo, Fxr, Psi, phis)
END SUBROUTINE

SUBROUTINE KinParamGen(Core, FmInfo, TranInfo, ThInfo, GroupInfo, lBetaUpdt, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,    FmInfo_Type,      ThInfo_Type,          & 
                             TranInfo_Type,    GroupInfo_Type,   PE_Type,              &
                             FXRInfo_Type,     PinXs_Type,       Pin_Type,             &              
                             ResVarPin_Type, PinInfo_Type,     Cell_Type,        XsMac_Type
USE CNTL,             ONLY : nTracerCntl_Type
USE BenchXs,          ONLY : XsBaseBen,        DnpBetaBen,        NeutVeloBen,         &
                             DnpBetaDynBen,    NeutVeloDynBen
USE MacXsLib_mod,     ONLY : EffMacXS,         EffRIFPin,      MacXsBase
USE MOC_MOD,          ONLY : FxrAvgPhi
USE TranMacXsLib_Mod, ONLY : FxrBeta,          FxrVelo
USE BasicOperation,   ONLY : CP_VA,            CP_CA,             MULTI_VA
USE TH_Mod,           ONLY : GetPinFuelTemp
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) ::  GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
LOGICAL :: lBetaUpdt
!
TYPE(FxrInfo_Type), POINTER :: FXR(:, :)
TYPE(FxrInfo_Type), POINTER :: myFXR
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(ResVarPin_Type), POINTER :: ResVarPin(:,:)
TYPE(Cell_Type), POINTER :: CellInfo(:)

TYPE(XsMac_Type), SAVE :: XsMac(nThreadMax) 

INTEGER :: myzb, myze, nxy, nLocalFxr, nprec, ng, nchi
INTEGER :: FxrIdxSt
INTEGER :: iresoGrp1, iresoGrp2
INTEGER :: ifxr, iz, ixy, icel, ig, tid
INTEGER :: i, j

REAL :: betat
REAL :: Beta0(1:50), velo0(500), chid(1:500), veloh(500)
REAL :: FxrPhi(500), Phisum(500)
REAL :: Temp!, PinFuelTempAvgsq

INTEGER :: niso
INTEGER, POINTER :: idiso(:)
REAL, POINTER :: pnum(:)

LOGICAL :: lXsLib, lRes, lAIC

lXsLib = nTracerCntl%lXsLib
IF(.NOT. lXsLib) RETURN
Fxr => FmInfo%Fxr
Pin => Core%Pin; CellInfo => Core%CellInfo

ResVarPin => Core%ResVarPin

myzb = PE%myzb; myze = PE%myze
nxy = Core%nxy
nprec = GroupInfo%nprec; ng = GroupInfo%ng
nchi = GroupInfo%nchi


iresoGrp1 = GroupInfo%nofg + 1; iresoGrp2 = GroupInfo%nofg + GroupInfo%norg
Chid(1:ng) = TranInfo%Chid(1:ng)

tid = 1
DO iz = myzb, myze
  !! The routine 'EffMacXs' is not parallelized.
  !! Instead, pin-wise parallelization is adopted.
  !! Refer to the routine, 'SubGrpEffXsGen', for parallelization. (2018-03-24 by PHS)
  DO ixy = 1, nxy
    FxrIdxSt = Pin(ixy)%FxrIdxSt; icel = Pin(ixy)%Cell(iz)
    lAIC = CellInfo(icel)%lAIC
    nlocalFxr = CellInfo(icel)%nFxr
    Phisum(1:ng) = 0; Veloh(1:ng) = 0
    !PinFuelTempAvgsq = dsqrt(GetPinFuelTemp(Core, Fxr, iz, ixy))
    !IF (ResVarPin(ixy,iz)%lres.and.nTracerCntl%lRIF) CALL EffRIFPin(ResVarPin(ixy,iz), PinFuelTempAvgsq, iz, lAIC, PE)
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j - 1
      myFxr => Fxr(ifxr, iz)

      !IF(.NOT. Fxr(ifxr, iz)%lfuel) CYCLE
      IF(lXsLib) THEN
        !Calculate Beta 
        niso = myFxr%niso; idiso => myFxr%idiso
        pnum => myFxr%pnum; lres = myFxr%lres
        temp = myFxr%temp
        FxrPhi(1:ng) = FxrAvgPhi(Core, Fxr, FmInfo%Phis,ixy, j, iz, ng, PE)
              
        IF(lBetaUpdt) THEN
          CALL MacXsBase(XsMac(tid), myFxr, 1, ng, ng, 1._8, TRUE, TRUE)
          !IF (lres) THEN            
          !  DO ig = iresoGrp1, iresoGrp2
          !    CALL EffMacXs(XsMac(tid), ResVarPin(ixy,iz), myFxr, PinFuelTempAvgsq, niso, ig, ng, .TRUE., iz, ixy, j, PE)                
          !  ENDDO
          !ENDIF
          CALL FxrBeta(XsMac(tid), myFxr, FxrPhi(1:ng), ng, iresoGrp1, iresoGrp2)  
        ENDIF
        !Calculate neutron Velocity
        CALL FxrVelo(myFxr, Temp, ng, nTracerCntl%lfixvel)
        DO ig = 1, ng
          Veloh(ig) = Veloh(ig) + 1._8 / myFxr%velo(ig) * FxrPhi(ig) * myFxr%area
          PhiSum(ig) = PhiSum(ig) + Fxr(ifxr, iz)%area * FxrPhi(ig)    
        ENDDO
        Betat = sum(myFxr%beta(1:nprec))
        NULLIFY(pnum, idiso)
      ELSE
        IF(nTracerCntl%lDynamicBen) THEN
          CALL DnpBetaDynBen(myFxr%imix, TranInfo%fuelTemp(ixy, iz), beta0(1:nprec))
          CALL NeutVeloDynBen(myFxr%imix, TranInfo%fuelTemp(ixy, iz), velo0(1:ng))
        ELSE
          CALL DnpBetaBen(myFxr%imix, beta0(1:nprec))
          CALL NeutVeloBen(myFxr%imix, velo0(1:ng))
        END IF
        MyFxr%beta(1:nprec) = beta0(1:nprec)
        MyFxr%velo(1:nprec) = velo0(1:nprec)
        betat = sum(beta0(1:nprec))
      ENDIF
      
      
      MyFxr%betat = betat
      !IF(myFxr%lfuel) myFxr%Chip(1:nchi) = (myFxr%Chi(1:nchi) - betat * Chid(1:nchi)) / ( 1._8 - betat)
    ENDDO
    IF(lXsLib) THEN
      DO ig = 1, ng
        Veloh(ig) = Veloh(ig) / phisum(ig)
        Veloh(ig) = 1._8 / Veloh(ig)
      ENDDO
      DO j = 1, nLocalFxr
        ifxr = FxrIdxSt + j - 1
        Fxr(ifxr, iz)%veloh(1:ng) = veloh(1:ng)
      ENDDO
      CONTINUE
    ENDIF
  ENDDO
ENDDO

NULLIFY(CellInfo, Pin, Fxr)

END SUBROUTINE

SUBROUTINE SetCmfdPrecParam(Core, CmInfo, TranInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,       ONLY : CoreInfo_Type,          CmInfo_Type,          TranInfo_Type,     &
                          TranCntl_Type,          PE_Type,                                 &
                          PinXs_Type
USE CNTL,          ONLY : nTracerCntl_Type
USE BasicOperation,ONLY : CP_CA
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCNtl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE

TYPE(PinXS_TYPE), POINTER :: PinXs(:, :)
REAL, POINTER :: CellOmegam(:, :, :), CellOmega0(:, :, :), CellOmegap(:, :, :)

INTEGER :: nprec, nowstep, norder
INTEGER :: nxy, myzb, myze
INTEGER :: i, ixy, iz

REAL :: Delt, Deltp, gamma, Invgamma, rgammap1
REAL :: lambda(50), InvLambda(50)
REAL :: rldt(50), rldtgp1(50), kapbrldt(50), kapbrldt2(50), Kappa(50), kappap1(50)
REAL :: omegam(50), omega0(50), omegap(50)
REAL :: omegalm, omegal0, omegalp

nprec = TranInfo%nprec; nxy = Core%nxy
myzb = PE%myzb; myze = PE%myze
norder = 2
NowStep = TranCntl%NowStep
IF(NowStep .EQ. 1 .OR.  abs(TranCntl%theta - 0.5) .GT. epsm6) norder =1 !2
Delt = TranCntl%Delt(nowstep)
Deltp = TranCntl%Delt(nowstep)
IF(NowStep .GT. 1) Deltp = TranCntl%Delt(nowstep - 1)

DO i = 1, nprec
  lambda(i) = TranInfo%lambda(i)    
  InvLambda(i) = 1._8 / TranInfo%lambda(i)
  kappa(i) = exp(-lambda(i) * DelT)
  kappap1(i) = kappa(i) + 1._8
  rldt(i) = 1._8 / (deltp * lambda(i))
  kapbrldt(i) = (1._8-kappa(i)) * rldt(i)
  IF(norder .EQ. 2) THEN
    gamma = Delt/ Deltp
    InvGamma = 1._8 / Gamma
    rgammap1 = 1._8 / (gamma + 1)
    rldtgp1(i) = rldt(i) * rgammap1
    kapbrldt2(i) = (1._8 - kappa(i)) * (1._8 - 2._8 * rldt(i))
  ENDIF
ENDDO
IF(norder .EQ. 1) THEN
  DO i = 1, nprec
    omegam(i) = 0
    omega0(i) = InvLambda(i) * (Kapbrldt(i) - Kappa(i))
    omegap(i) = InvLambda(i) * (1- Kapbrldt(i))
  ENDDO
ELSE
  DO i = 1, nprec    
    omegam(i) = InvLambda(i) * rldtgp1(i) * (2._8 *kapbrldt(i)- gamma * kappap1(i))
    omega0(i) = InvLambda(i) * (rldt(i) * (kappap1(i) + kapbrldt2(i) * InvGamma)-Kappa(i))
    omegap(i) = InvLambda(i) * (1 -rldtgp1(i) * (2._8 + kapbrldt2(i)* InvGamma))
  ENDDO
ENDIF

CellOmegam => TranInfo%CellOmegam;   CellOmega0 => TranInfo%CellOmega0
CellOmegap => TranInfo%CellOmegap
PinXs => CmInfo%PinXS

CALL CP_CA(CellOmegam(0:nprec, 1:nxy, myzb:myze), 0._8, nprec+1, nxy, myze - myzb + 1)
CALL CP_CA(CellOmega0(0:nprec, 1:nxy, myzb:myze), 0._8, nprec+1, nxy, myze - myzb + 1)
CALL CP_CA(CellOmegap(0:nprec, 1:nxy, myzb:myze), 0._8, nprec+1, nxy, myze - myzb + 1)

DO iz = myzb, myze
  DO ixy = 1, nxy
    omegalm = 0; omegal0 = 0; omegalp = 0
    DO i = 1, nprec
      CellOmegam(i, ixy, iz) = PinXs(ixy, iz)%beta(i) * omegam(i)
      CellOmega0(i, ixy, iz) = PinXs(ixy, iz)%beta(i) * omega0(i)
      CellOmegap(i, ixy, iz) = PinXs(ixy, iz)%beta(i) * omegap(i)
      
      omegalm = omegalm + CellOmegam(i, ixy, iz) * lambda(i)
      omegal0 = omegal0 + CellOmega0(i, ixy, iz) * lambda(i)
      omegalp = omegalp + CellOmegap(i, ixy, iz) * lambda(i)
    ENDDO
    CellOmegam(0, ixy, iz) = omegalm
    CellOmega0(0, ixy, iz) = omegal0
    CellOmegap(0, ixy, iz) = omegalp
    
    PinXs(ixy, iz)%omega = omegalp
  ENDDO
ENDDO
!PRINT '(5es12.3)', (PinXS(i, 1)%omega , i = 51, 55)

NULLIFY(CellOmegam, CellOmega0, CellOmegap)
NULLIFY(PinXS)
END SUBROUTINE

SUBROUTINE PrepareExpTrsf(Core, FmInfo, CmInfo, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)
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

nowstep = TranCntl%NowStep + 1
Delt = TranCntl%DelT(nowstep)

DO iz = myzb, myze
  DO ixy = 1, nxy
    DO ig = 1, ng
      Expo_Alpha(ixy, iz, ig) = log(abs(PhiC(ixy, iz, ig)/TranPhiCm(ixy, iz, ig))) / delt
      Expo(ixy, iz, ig) = Exp(Expo_Alpha(ixy, iz, ig) * Delt)
    ENDDO
    !pinv(m,l,k)=log(abs(phi(m,l,k)/phid(m,l,k)))/deltmd
  ENDDO
ENDDO

!DO iz = myzb, myze
!  DO ixy = 1, nxy
!    FsrIdxSt = Pin(ixy)%FsrIdxSt; icel = Pin(ixy)%Cell(iz)
!    nLocalFsr = CellInfo(icel)%nFsr
!    DO j = 1, nLocalFsr
!      ifsr = FsrIdxSt + j - 1
!      Phis(ifsr, iz, ig) = Expo(ixy, iz, ig) * Phis(ifsr, iz, ig)
!    ENDDO
!  ENDDO
!ENDDO
!
!
END SUBROUTINE

SUBROUTINE SetPrecParam(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,         ONLY : CoreInfo_Type,          FmInfo_Type,          TranInfo_Type,     &
                            TranCntl_Type,          PE_Type,                                 &
                            FxrInfo_Type,           Pin_Type,             Cell_Type
USE CNTL,            ONLY : nTracerCntl_Type
USE BenchXs,         ONLY : DnpBetaBen,             DnpBetaDynBen
USE BasicOperation,  ONLY : CP_CA
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCNtl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)

REAL, POINTER :: FxrOmegam(:, :, :), FxrOmega0(:, :, :), FxrOmegap(:, :, :)

INTEGER :: nprec, nowstep, norder
INTEGER :: nxy, nfxr, nFsr, myzb, myze
INTEGER :: i, ifxr, ixy, iz

REAL :: Delt, Deltp, gamma, Invgamma, rgammap1
REAL :: lambda(50), InvLambda(50)
REAL :: rldt(50), rldtgp1(50), kapbrldt(50), kapbrldt2(50), Kappa(50), kappap1(50)
REAL :: omegam(50), omega0(50), omegap(50), beta(50)
REAL :: omegalm, omegal0, omegalp

LOGICAL :: lXsLib

nprec = TranInfo%nprec; nxy = Core%nxy
nFxr = Core%nCoreFxr; nFsr = Core%nCoreFsr
myzb = PE%myzb; myze = PE%myze
lXsLib = nTracerCntl%lXsLib


norder = 2
NowStep = TranCntl%NowStep
!IF(NowStep .EQ. 1 .OR. abs(TranCntl%theta - 0.5) .GT. epsm6) norder = 1
IF(NowStep .EQ. 1) norder = 1
Delt = TranCntl%Delt(nowstep); Deltp =  TranCntl%Delt(nowstep)
IF(norder .EQ. 2) Deltp = TranCntl%Delt(nowstep - 1)

DO i = 1, nprec
  lambda(i) = TranInfo%lambda(i)    
  InvLambda(i) = 1._8 / TranInfo%lambda(i)
  kappa(i) = exp(-lambda(i) * DelT)
  kappap1(i) = kappa(i) + 1._8
  rldt(i) = 1._8 / (delTp * lambda(i))
  kapbrldt(i) = (1._8-kappa(i)) * rldt(i)
  IF(norder .EQ. 2) THEN
    gamma = Delt/ Deltp
    InvGamma = 1._8 / Gamma
    rgammap1 = 1._8 / (gamma + 1)
    rldtgp1(i) = rldt(i) * rgammap1
    kapbrldt2(i) = (1._8 - kappa(i)) * (1._8 - 2._8 * rldt(i))
  ENDIF
ENDDO
IF(norder .EQ. 1) THEN
  DO i = 1, nprec
    omegam(i) = 0
    omega0(i) = InvLambda(i) * (Kapbrldt(i) - Kappa(i))
    omegap(i) = InvLambda(i) * (1- Kapbrldt(i))
  ENDDO
ELSE
  DO i = 1, nprec
    omegam(i) = InvLambda(i) * rldtgp1(i) * (2._8 *kapbrldt(i)- gamma * kappap1(i))
    omega0(i) = InvLambda(i) * (rldt(i) * (kappap1(i) + kapbrldt2(i) * InvGamma)-Kappa(i))
    omegap(i) = InvLambda(i) * (1 -rldtgp1(i) * (2._8 + kapbrldt2(i)* InvGamma))
  ENDDO
ENDIF

Fxr => FmInfo%Fxr
Pin => Core%Pin; CellInfo => Core%CellInfo
FxrOmegam => TranInfo%FxrOmegam;   FxrOmega0 => TranInfo%FxrOmega0
FxrOmegap => TranInfo%FxrOmegap

CALL CP_CA(FxrOmegam(0:nprec, 1:nFxr, myzb:myze), 0._8, nprec+1, nfxr, myze - myzb + 1)
CALL CP_CA(FxrOmega0(0:nprec, 1:nFxr, myzb:myze), 0._8, nprec+1, nfxr, myze - myzb + 1)
CALL CP_CA(FxrOmegap(0:nprec, 1:nFxr, myzb:myze), 0._8, nprec+1, nfxr, myze - myzb + 1)

DO iz = myzb, myze
  DO ifxr = 1, nfxr
    FxrOmegam(:, ifxr, iz) = 0; FxrOmega0(:, ifxr, iz) = 0; 
    FxrOmegap(:, ifxr, iz) = 0
    IF(.NOT. Fxr(ifxr, iz)%lFuel) CYCLE

    IF(lXsLib) THEN
      Beta(1:nprec) = Fxr(ifxr, iz)%beta(1:nprec)
    ELSE
      IF(TranCntl%lDynamicBen) THEN
        ixy = Fxr(ifxr, iz)%ipin
        CALL DnpBetaDynBen(Fxr(ifxr, iz)%imix, TranInfo%fuelTemp(ixy, iz), beta(1:nprec))
      ELSE
        CALL DnpBetaBen(Fxr(ifxr, iz)%imix, beta(1:nprec))
      END IF
    ENDIF
    
    omegalm = 0; omegal0 = 0; omegalp = 0
    DO i = 1, nprec
      FxrOmegam(i, ifxr, iz) = beta(i) * omegam(i)
      FxrOmega0(i, ifxr, iz) = beta(i) * omega0(i)
      FxrOmegap(i, ifxr, iz) = beta(i) * omegap(i)
      
      omegalm = omegalm + FxrOmegam(i, ifxr, iz) * lambda(i)
      omegal0 = omegal0 + FxrOmega0(i, ifxr, iz) * lambda(i)
      omegalp = omegalp + FxrOmegap(i, ifxr, iz) * lambda(i)
    ENDDO
    FxrOmegam(0, ifxr, iz) = omegalm
    FxrOmega0(0, ifxr, iz) = omegal0
    FxrOmegap(0, ifxr, iz) = omegalp
  ENDDO
ENDDO

NULLIFY(FxrOmegam, FxrOmega0, FxrOmegap)
NULLIFY(Fxr, Pin, CellInfo)
END SUBROUTINE

SUBROUTINE SetCorrectorPrecParam(TranInfo, omegam, omega0, omegap, beta, factor_F, delt, deltp, nprec, norder)
USE TYPEDEF,        ONLY : TranInfo_Type
IMPLICIT NONE
TYPE(TranInfo_Type) :: TranInfo
REAL :: omegam(nprec), omega0(nprec), omegap(nprec), beta(nprec)
REAL :: factor_F
REAL :: delt, deltp
INTEGER :: nprec, norder

REAL :: lambda(nprec), invlambda(nprec)
REAL :: invldt(nprec), invldtgp1(nprec), kapbinvldt(nprec), kapbinvldt2(nprec), kappa(nprec), kappap1(nprec)
REAL :: gamma, invgamma, invgammap1
REAL :: coeff
INTEGER :: i

IF(norder .EQ. 2) THEN
  gamma = delt / deltp
  invgamma = 1._8 / gamma
  invgammap1 = 1._8 / (gamma + 1._8)
END IF
DO i = 1, nprec
  lambda(i) = TranInfo%lambda(i)
  invlambda(i) = 1._8 / TranInfo%lambda(i)
  kappa(i) = exp(-lambda(i) * delt)
  kappap1(i) = kappa(i) + 1._8
  invldt(i) = 1._8 / (deltp * lambda(i))
  kapbinvldt(i) = (1._8 - kappa(i)) * invldt(i)
  IF(norder .EQ. 2) THEN
    invldtgp1(i) = invldt(i) * invgammap1
    kapbinvldt2(i) = (1._8 - kappa(i)) * (1._8 - 2._8 * invldt(i))
  END IF
END DO 

IF(norder .EQ. 1) THEN
  DO i = 1, nprec
    omegam(i) = 0
    omega0(i) = invlambda(i) * (kapbinvldt(i) - kappa(i))
    omegap(i) = invlambda(i) * (1 - kapbinvldt(i))
  END DO 
ELSE
  DO i = 1, nprec
    omegam(i) = invlambda(i) * invldtgp1(i) * (2._8 * kapbinvldt(i) - gamma * kappap1(i))
    omega0(i) = invlambda(i) * (invldt(i) *(kappap1(i) + kapbinvldt2(i) * invgamma) - kappa(i)) 
    omegap(i) = invlambda(i) * (1 - invldtgp1(i) * (2._8 + kapbinvldt2(i) * invgamma))
  END DO 
END IF

coeff = factor_F * TranInfo%Inv_factor_F0
DO i = 1, nprec
  omegam(i) = coeff * beta(i) * omegam(i)
  omega0(i) = coeff * beta(i) * omega0(i)
  omegap(i) = coeff * beta(i) * omegap(i)
END DO 



END SUBROUTINE

  