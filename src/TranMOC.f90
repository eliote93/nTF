#include <defines.h>
SUBROUTINE PrecSrcUpdt(Core, Fxr, PrecSrc, Prec, TranPsi, TranPsid, GroupInfo, TranCntl, nTRACERCntl, PE)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,           FxrInfo_Type,            TranCntl_Type,        &
                           GroupInfo_Type,          PE_Type,                                       &
                           Pin_Type,                Cell_Type
USE CNTL,           ONLY : nTracerCntl_Type
USE TRANMOC_MOD,    ONLY : Omegam,                  Omegap,               &
                           Omega0,                  Chid,                    Lambda
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

REAL, POINTER :: Prec(:, :, :), PrecSrc(:, :)
REAL, POINTER :: TranPsi(:, :), TranPsid(:, :)

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)

INTEGER :: iz, ifsr, ifxr, icel, iz0, ixy, iprec, nowstep
INTEGER:: i, j
INTEGER :: myzb, myze, nfsr, nfxr, nxy, nprec
INTEGER :: nLocalFxr, nFsrInFxr, FxrIdxSt, FsrIdxSt
REAL :: DelT, reigv, kappa(100)

Pin => Core%Pin; CellInfo => Core%CellInfo
nfsr = Core%nCoreFsr; nfxr = Core%nCoreFxr
myzb = PE%myzb; myze = PE%myze; nxy = Core%nxy
nprec = GroupInfo%nprec

NowStep = TranCntl%NowStep; Delt = TranCntl%DelT(NowStep)
DO iprec = 1, nprec
  kappa(iprec) = exp(-delt * lambda(iprec))
ENDDO

CALL CP_CA(PrecSrc(1:nfsr, myzb:myze), 0._8, nfsr, myze - myzb + 1)

DO iz = myzb, myze
  DO ixy = 1, nxy
    FxrIdxSt = Pin(ixy)%FxrIdxSt; FsrIdxst = Pin(ixy)%FsrIdxSt
    icel = Pin(ixy)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr 
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j - 1   
      IF(.NOT. Fxr(ifxr, iz)%lfuel) CYCLE
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)   
      DO i = 1, nFsrInFxr
        ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
        !Prvious time step Contribution
        DO iprec = 1, nPrec
          PrecSrc(ifsr, iz) = PrecSrc(ifsr, iz)  +  lambda(iprec) * kappa(iprec) * Prec(iprec, ifsr, iz)
        ENDDO
        !Fission Source Contributio term
        PrecSrc(ifsr, iz) = PrecSrc(ifsr, iz) + Omegam(0, ifxr, iz) * TranPsid(ifsr, iz)
        PrecSrc(ifsr, iz) = PrecSrc(ifsr, iz) + Omega0(0, ifxr, iz) * TranPsi(ifsr, iz)
      ENDDO
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE 

SUBROUTINE PrecSrcKUpdt(Core, Fxr, PrecSrcK, Prec, TranPsi, TranPsid, GroupInfo, TranCntl, nTRACERCntl, PE)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,           FxrInfo_Type,            TranCntl_Type,        &
                           GroupInfo_Type,          PE_Type,                                       &
                           Pin_Type,                Cell_Type
USE CNTL,           ONLY : nTracerCntl_Type
USE TRANMOC_MOD,    ONLY : Omegam,                  Omegap,               &
                           Omega0,                  Chid,                    Lambda
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

REAL, POINTER :: Prec(:, :, :), PrecSrcK(:, :, :)
REAL, POINTER :: TranPsi(:, :), TranPsid(:, :)

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)

INTEGER :: iz, ifsr, ifxr, icel, iz0, ixy, iprec, nowstep
INTEGER:: i, j
INTEGER :: myzb, myze, nfsr, nfxr, nxy, nprec
INTEGER :: nLocalFxr, nFsrInFxr, FxrIdxSt, FsrIdxSt
REAL :: DelT, reigv, kappa(100)

Pin => Core%Pin; CellInfo => Core%CellInfo
nfsr = Core%nCoreFsr; nfxr = Core%nCoreFxr
myzb = PE%myzb; myze = PE%myze; nxy = Core%nxy
nprec = GroupInfo%nprec

NowStep = TranCntl%NowStep; Delt = TranCntl%DelT(NowStep)
DO iprec = 1, nprec
  kappa(iprec) = exp(-delt * lambda(iprec))
ENDDO

DO iz = myzb, myze
  DO ixy = 1, nxy
    FxrIdxSt = Pin(ixy)%FxrIdxSt; FsrIdxst = Pin(ixy)%FsrIdxSt
    icel = Pin(ixy)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr 
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j - 1   
      IF(.NOT. Fxr(ifxr, iz)%lfuel) CYCLE
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)   
      DO i = 1, nFsrInFxr
        ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
        DO iprec = 1, nprec
          PrecSrcK(iprec, ifsr, iz) = lambda(iprec) * kappa(iprec) * Prec(iprec, ifsr, iz)
          PrecSrcK(iprec, ifsr, iz) = PrecSrcK(iprec, ifsr, iz) + Omegam(iprec, ifxr, iz) * lambda(iprec) * TranPsid(ifsr, iz)
          PrecSrcK(iprec, ifsr, iz) = PrecSrcK(iprec, ifsr, iz) + Omega0(iprec, ifxr, iz) * lambda(iprec) * TranPsi(ifsr, iz)
        END DO 
      ENDDO
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE 

SUBROUTINE SetTranSrcNM(Core, FmInfo, Fxr, TranSrcnm, Phinm, TranPhinm, Psi, ResSrc, xstnm, iz, &
                        gb, ge, GroupInfo, TranInfo, TranCntl, lxslib, PE, Offset)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,         FxrInfo_Type,         GroupInfo_Type,         TranCntl_Type,      &
                           PE_Type,               Pin_Type,             Cell_Type,              TranInfo_Type,      &
                           FmInfo_Type
USE BenchXs,        ONLY : XsBaseBen,             DnpBetaBen,           NeutVeloBen,                                &
                           DnpBetaDynBen,         NeutVeloDynBen
USE TRANMOC_MOD,    ONLY : Omegam,                Omegap,               Omega0,                 chid,               &
                           Lambda,                Expo,                 Expo_Alpha
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
REAL, POINTER :: TranSrcnm(:, :), Phinm(:, :), TranPhinm(:, :), Psi(:, :), ResSrc(:, :, :), xstnm(:, :)
INTEGER :: iz, gb, ge
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
LOGICAL :: lxslib
TYPE(PE_Type) :: PE
INTEGER, OPTIONAL :: Offset   !--- CNJ Edit : GPU Acceleration

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
REAL, POINTER :: PinVol(:, :)
REAL :: delt, chieff, betat, omega, rvdt, theta, thetah, rvol
REAL :: Beta0(1:50), velo0(1:50)
REAL :: PrevSrc, ExpSrc, rvalpha
INTEGER :: myzb, myze, nfsr, nfxr, nxy, nprec, ng
INTEGER :: nLocalFxr, nFsrInFxr, FxrIdxSt, FsrIdxSt
INTEGER :: ifsr, ifxr, icel, ixy, iprec, ig, nowstep
INTEGER :: i, j
INTEGER :: xyb, xye, off
LOGICAL :: lmethod

lMethod = TranCntl%lMethod

Pin => Core%Pin
CellInfo => Core%CellInfo
PinVol => Core%PinVol
nprec = GroupInfo%nprec
ng = GroupInfo%ng
nfsr = Core%nCoreFsr
nfxr = Core%nCoreFxr
nxy = Core%nxy
xyb = PE%myPinBeg; xye = PE%myPinEnd   !--- CNJ Edit : Domain Decomposition + MPI
off = 0
IF (PRESENT(Offset)) off = Offset

nowstep = TranCntl%nowstep
delt = TranCntl%delt(nowstep)
theta = TranCntl%theta
thetah = 1._8 / theta - 1._8

!$OMP PARALLEL DEFAULT(SHARED)                                                               & 
!$OMP PRIVATE(ixy, FsrIdxSt, FxrIdxSt, icel, nlocalFxr, rvol, ifxr, nFsrinFxr, Beta0, Velo0, &
!$OMP         betat, i, j, ig, ifsr, chieff, rvdt, PrevSrc, rvalpha, ExpSrc)
!$OMP DO SCHEDULE(GUIDED)
DO ixy = xyb, xye
  FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
  icel = Pin(ixy)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  rvol = 1._8 / PinVol(ixy, iz)
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1
    nFsrinFxr = CellInfo(icel)%nFsrInFxr(j)
    IF(lxslib) THEN
      Beta0(1:nprec) = Fxr(ifxr, iz)%Beta(1:nprec)
      Velo0(1:ng) = Fxr(ifxr, iz)%veloh(1:ng)
    ELSE
      IF(TranCntl%lDynamicBen) THEN
        CALL DnpBetaDynBen(Fxr(ifxr, iz)%imix, TranInfo%fuelTemp(ixy, iz), beta0(1:nprec))
        CALL NeutVeloDynBen(Fxr(ifxr, iz)%imix, TranInfo%fuelTemp(ixy, iz), velo0(1:ng))
      ELSE
        CALL DnpBetaBen(Fxr(ifxr, iz)%imix, beta0(1:nprec))
        CALL NeutVeloBen(Fxr(ifxr, iz)%imix, velo0(1:ng))
      END IF
    END IF
    betat = SUM(beta0(1:nprec))
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + CellInfo(icel)%MapFxr2FsrIdx(i, j) - 1
      DO ig = gb, ge
        !chieff = (Omegap(0, ifxr, iz) - betat) * chid(ig)
        IF(TranCntl%lchidk) THEN
          chieff = 0.
          DO iprec = 1, nprec
            chieff = chieff + (Omegap(iprec, ifxr, iz) * TranInfo%lambda(iprec) - beta0(iprec))* TranInfo%chidk(ig, iprec)
          END DO 
        ELSE
          chieff = - betat * chid(ig)
          chieff = chieff + Omegap(0, ifxr, iz) * chid(ig)
        END IF
        rvdt = 1._8 / (delt * velo0(ig) * theta)
        IF(TranCntl%lchidk) THEN 
          TranSrcnm(ig-off, ifsr) = 0.
          DO iprec = 1, nprec
            TranSrcnm(ig-off, ifsr) = TranSrcnm(ig-off, ifsr) + FmInfo%PrecSrcK(iprec, ifsr, iz) * TranInfo%chidk(ig, iprec)
          END DO
        ELSE
          TranSrcnm(ig-off, ifsr) = FmInfo%PrecSrc(ifsr, iz) * chid(ig)
        END IF
        TranSrcnm(ig-off, ifsr) = TranSrcnm(ig-off, ifsr) + chieff * Psi(ifsr, iz)

        PrevSrc = Thetah * ResSrc(ixy, iz, ig) * rvol * Expo(ifsr, iz, ig)
        PrevSrc = PrevSrc - rvdt * (Phinm(ig, ifsr) - Expo(ifsr, iz, ig) * TranPhinm(ig, ifsr))

        rvalpha = expo_alpha(ifsr, iz, ig) / velo0(ig)
        ExpSrc = - rvalpha * (Phinm(ig, ifsr) + Thetah * Expo(ifsr, iz, ig) * TranPhinm(ig, ifsr))

        TranSrcnm(ig-off, ifsr) = TranSrcnm(ig-off, ifsr) + PrevSrc + ExpSrc
        TranSrcnm(ig-off, ifsr) = TranSrcnm(ig-off, ifsr) / xstnm(ig, ifsr)
      END DO 
    END DO 
  END DO 
END DO 
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE


SUBROUTINE SetTranSrc(Core, Fxr, TranSrc, Phi, TranPhi, Psi, PrecSrc, ResSrc, xstr, iz, ig, GroupInfo, TranInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,            GroupInfo_Type,        TranCntl_Type,            &
                           FxrInfo_Type,             PE_Type,               TranInfo_Type,            &
                           Pin_Type,                 Cell_Type
USE CNTL,           ONLY : nTracerCntl_Type
USE BenchXs,        ONLY : XsBaseBen,                DnpBetaBen,             NeutVeloBen,            &
                           DnpBetaDynBen,            NeutVeloDynBen
USE TRANMOC_MOD,    ONLY : Omegam,                  Omegap,                                          &
                           Omega0,                  Chid,                    Lambda,                 &
                           Expo,                    Expo_Alpha
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
INTEGER :: iz, ig

REAL, POINTER :: TranSrc(:)
REAL, POINTER :: Phi(:, :, :)
REAL, POINTER :: TranPhi(:, :, :)
REAL, POINTER :: Psi(:, :)
REAL, POINTER :: PrecSrc(:, :)
REAL, POINTER :: ResSrc(:, :, :)
REAL, POINTER :: xstr(:)

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
REAL, POINTER :: PinVol(:, :)

INTEGER :: ifsr, ifxr, icel, iz0, ixy, iprec, nowstep
INTEGER:: i, j
INTEGER :: myzb, myze, nfsr, nfxr, nxy, nprec, ng
INTEGER :: nLocalFxr, nFsrInFxr, FxrIdxSt, FsrIdxSt
REAL :: DelT, chieff, betat, omega, rvdt, theta, thetah, rvol
REAL :: Beta0(1:50), velo0(500)
REAL :: PrevSrc, ExpSrc, RVALPHA
REAL :: PhiC, TranPhiC, Area

LOGICAL :: lXsLib

lXsLib = nTracerCntl%lXsLib

Pin => Core%Pin; CellInfo => Core%CellInfo
PinVol => Core%PinVol
nfsr = Core%nCoreFsr; nfxr = Core%nCoreFxr
ng = GroupInfo%ng; nprec = GroupInfo%nprec
nxy = Core%nxy; nfsr = Core%nCoreFsr; nfxr = Core%nCoreFxr

NowStep = TranCntl%NowStep; Delt = TranCntl%DelT(NowStep)
Theta = TranCntl%Theta
ThetaH = 1._8 / Theta - 1._8
DO ixy = 1, nxy
  FxrIdxSt = Pin(ixy)%FxrIdxSt; FsrIdxst = Pin(ixy)%FsrIdxSt
  icel = Pin(ixy)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  rvol = 1._8 / PinVol(ixy, iz)
  Area = 0; PhiC = 0; TranPhiC = 0
  DO j = 1, CellInfo(icel)%nFxr
    ifsr = FsrIdxSt + j - 1
    Area = Area + CellInfo(icel)%vol(j)
    PhiC = PhiC + CellInfo(icel)%vol(j) * Phi(ifsr, iz, ig)
    TranPhiC = TranPhiC + CellInfo(icel)%vol(j) * TranPhi(ifsr, iz, ig)
  ENDDO
  PhiC = PhiC / Area; TranPhiC = TranPhiC / Area 
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1   
    !IF(.NOT. Fxr(ifxr, iz)%lfuel) CYCLE
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    IF(lXsLib) THEN
      Beta0(1:nprec) = Fxr(ifxr, iz)%Beta(1:nprec)
      Velo0(1:ng) = Fxr(ifxr, iz)%veloh(1:ng)
    ELSE
      IF(TranCntl%lDynamicBen) THEN
        CALL DnpBetaDynBen(Fxr(ifxr, iz)%imix, TranInfo%fuelTemp(ixy, iz), beta0(1:nprec))
        CALL NeutVeloDynBen(Fxr(ifxr, iz)%imix, TranInfo%fuelTemp(ixy, iz), velo0(1:ng))
      ELSE
        CALL DnpBetaBen(Fxr(ifxr, iz)%imix, beta0(1:nprec))
        CALL NeutVeloBen(Fxr(ifxr, iz)%imix, velo0(1:ng))
      END IF
    ENDIF
    betat = SUM(beta0(1:nprec))
    chieff = (Omegap(0, ifxr, iz) - betat) *chid(ig)
    rvdt = 1._8 / (DelT * velo0(ig) * theta)
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      TranSrc(ifsr) = PrecSrc(ifsr, iz)*chid(ig)
      TranSrc(ifsr) = TranSrc(ifsr) + chieff * Psi(ifsr, iz)

      PrevSrc = ThetaH * ResSrc(ixy, iz, ig) * rvol * Expo(ifsr, iz, ig)
      PrevSrc = PrevSrc - rvdt * (Phi(ifsr, iz, ig) - Expo(ifsr, iz, ig) * TranPhi(ifsr, iz, ig))
      RVALPHA  = Expo_Alpha(ifsr, iz, ig) / velo0(ig)
      ExpSrc = - RVALPHA * (Phi(ifsr, iz, ig) + Thetah * Expo(ifsr, iz, ig) * TranPhi(ifsr, iz, ig))
      !ExpSrc =  RVALPHA * (Thetah * Phi(ifsr, iz, ig) - Thetah * Expo(ifsr, iz, ig) * TranPhi(ifsr, iz, ig))

      TranSrc(ifsr) = TranSrc(ifsr) + PrevSrc + ExpSrc
      TranSrc(ifsr) = TranSrc(ifsr) / xstr(ifsr)

    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE SetExpTrsfXs(Core, Fxr, xstr, iz, ig, GroupInfo, TranInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,            GroupInfo_Type,        TranCntl_Type,            &
                           FxrInfo_Type,             PE_Type,               TranInfo_Type,            &
                           Pin_Type,                 Cell_Type
USE CNTL,           ONLY : nTracerCntl_Type
USE BenchXs,        ONLY : XsBaseBen,                DnpBetaBen,             NeutVeloBen,            &
                           DnpBetaDynBen,            NeutVeloDynBen
USE TRANMOC_MOD,    ONLY : Omegam,                  Omegap,                                          &
                           Omega0,                  Chid,                    Lambda,                 &
                           Expo,                    Expo_Alpha
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
INTEGER :: iz, ig

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
REAL, POINTER :: xstr(:)

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)

INTEGER :: ifsr, ifxr, icel, iz0, ixy, iprec, nowstep
INTEGER:: i, j
INTEGER :: myzb, myze, nfsr, nfxr, nxy, nprec, ng
INTEGER :: nLocalFxr, nFsrInFxr, FxrIdxSt, FsrIdxSt
REAL :: DelT, theta, thetah
REAL :: velo0(500)
REAL :: PrevSrc, ExpSrc, RVALPHA

LOGICAL :: lXsLib

lXsLib = nTracerCntl%lXsLib

Pin => Core%Pin; CellInfo => Core%CellInfo
nfsr = Core%nCoreFsr; nfxr = Core%nCoreFxr
ng = GroupInfo%ng; nprec = GroupInfo%nprec
nxy = Core%nxy; nfsr = Core%nCoreFsr; nfxr = Core%nCoreFxr

NowStep = TranCntl%NowStep; Delt = TranCntl%DelT(NowStep)
Theta = TranCntl%Theta
ThetaH = 1._8 / Theta - 1._8
DO ixy = 1, nxy
  FxrIdxSt = Pin(ixy)%FxrIdxSt; FsrIdxst = Pin(ixy)%FsrIdxSt
  icel = Pin(ixy)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1   
    !IF(.NOT. Fxr(ifxr, iz)%lfuel) CYCLE
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    IF(lXsLib) THEN
      Velo0(1:ng) = Fxr(ifxr, iz)%veloh(1:ng)
    ELSE
      IF(TranCntl%lDynamicBen) THEN
        CALL NeutVeloDynBen(Fxr(ifxr, iz)%imix, TranInfo%fuelTemp(ixy, iz), velo0(1:ng))
      ELSE
        CALL NeutVeloBen(Fxr(ifxr, iz)%imix, velo0(1:ng))
      END IF
    ENDIF

    
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      RVALPHA = Expo_Alpha(ifsr, iz, ig) / velo0(ig)
      xstr(ifsr) = xstr(ifsr) + rvalpha + thetah * rvalpha
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE

FUNCTION TranMocResidualError(Core, FmInfo, TranInfo, eigv0, GroupInfo, TranCntl, nTRACERCntl, PE)
USE PARAM
USE TYPEDEF,           ONLY : CoreInfo_Type,        FmInfo_Type,         GroupInfo_Type,         &
                              PE_Type,              TranInfo_Type,       TranCntl_Type,          &
                              FxrInfo_Type,         Pin_Type,            Cell_Type
USE CNTL,              ONLY : nTracerCntl_Type
USE MOC_MOD,           ONLY : tSrc,                 xst1g,               AxSrc1g,                &
                              AxPxs1g,                                                           & 
                              SetRtSrcGM,             SetRtMacXsGM,        PseudoAbsorptionGM
USE TRANMOC_MOD,       ONLY : TrSrc,                PrecSrc,                                    &
                              PrecSrcUpdt,          SetTranSrc,          SetExpTrsfXs
USE BasicOperation,    ONLY : CP_CA,                CP_VA,               AD_VA
USE MPIComm_Mod,       ONLY : REDUCE

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
REAL :: Eigv0

REAL :: TranMocResidualError

REAL, POINTER :: Phis(:, :, :), Psi(:, :)
REAL, POINTER :: MocJout(:, :, :, :, :)
REAL, POINTER :: AxSrc(:, :, :), AxPxs(:, :, :)
REAL, POINTER :: Prec(:, :, :)
REAL, POINTER :: ResSrc(:, :, :)
REAL, POINTER :: TranPhi(:, :, :)

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)

REAL, POINTER :: Res(:)
REAL :: SrcSum, LocalResidual, LocalSrc, vol, Temp

INTEGER :: ng, myzb, myze, nxy, nFsr, nFxr, nLocalFxr, nLocalFsr, FsrIdxst, FxrIdxSt, nFsrInFxr
INTEGER :: ig, ifxr, ireg, ipin, icel, ixy, iz, ifsrlocal, itype, nbd
INTEGER :: i, j, k, l, m

LOGICAL :: lXsLib, l3dim, lNegFix, lRST


l3dim = nTracerCntl%l3Dim
lXsLib = nTracerCntl%lXsLib
lRST = nTracerCntl%lRST
lNegFix = .FALSE.

ng = GroupInfo%ng
Phis => FmInfo%Phis; Psi => FmInfo%Psi
MocJout => FmInfo%RadJout
AxSrc => FmInfo%AxSrc; AxPXS => FmInfo%AxPXS

Prec => FmInfo%Prec; ResSrc => FmInfo%ResSrc
TranPhi => FmInfo%TranPhi

Pin => Core%Pin
CellInfo => Core%CellInfo
Fxr => FmInfo%Fxr

myzb = PE%myzb; myze = PE%myze
nxy = Core%nxy; nbd = 4
nFxr = Core%nCoreFxr; nFsr = Core%nCoreFsr

ALLOCATE(Res(nxy))

CALL CP_CA(AxPxs1g(1:nxy), 0._8, nxy)
CALL CP_CA(AxSrc1g(1:nxy), 0._8, nxy)

TranMocResidualError = 0; SrcSum = 0

DO iz = myzb, myze
  CALL CP_CA(RES, 0._8, nxy)
  DO ig = 1, ng
    IF(nTracerCntl%l3dim) THEN
      CALL CP_VA(AxPxs1g(1:nxy), AxPxs(1:nxy, iz, ig), nxy)
      CALL CP_VA(AxSrc1g(1:nxy), AxSrc(1:nxy, iz, ig), nxy)      
    ENDIF
    CALL CP_CA(xst1g(1:nfsr), 1._8, nfsr) 
    CALL SetTranSrc(Core, Fxr, TrSrc, Phis, TranPhi, Psi, PrecSrc, ResSrc, xst1g,       &
                    iz, ig, GroupInfo, TranInfo, TranCntl, nTracerCntl, PE) 
    CALL SetRtSrcGM(Core, Fxr(:, iz), tsrc, phis, psi, axSrc1g, xst1g,                    &
                  1._8, iz, ig, ng, GroupInfo, l3dim, lXslib, FALSE, lNegFix, PE)    
  
     
    CALL AD_VA(TrSrc(1:nfsr), TrSrc(1:nfsr), tsrc(1:nfsr), nfsr)
    DO ipin = 1, nxy
      FsrIdxSt = Pin(ipin)%FsrIdxSt
      icel = Pin(ipin)%Cell(iz); nlocalFsr = CellInfo(icel)%nFsr
      DO j = 1, nlocalFsr
        ireg = FsrIdxSt + j - 1
        vol = CellInfo(icel)%vol(j)
        Res(ipin) = Res(ipin)  + vol * Trsrc(ireg)
      ENDDO
    ENDDO !Pin Cell Sweep      
  ENDDO
  
  !Get total source term
  DO ipin = 1, nxy
    SrcSum = SrcSum + Res(ipin) * Res(ipin)
  ENDDO  
  
  DO ig = 1, ng
    IF(nTracerCntl%l3dim) THEN
      CALL CP_VA(AxPxs1g(1:nxy), AxPxs(1:nxy, iz, ig), nxy)
      CALL CP_VA(AxSrc1g(1:nxy), AxSrc(1:nxy, iz, ig), nxy)      
    ENDIF
    CALL SetRtMacXsGM(Core, Fxr(:, iz), xst1g, iz, ig, ng, lxslib, TRUE, lRST, FALSE, FALSE, PE)
    IF(l3dim)  CALL PseudoAbsorptionGM(Core, Fxr(:, iz), tsrc, phis(:, iz, ig),             &
                                     AxPXS1g(:), xst1g, iz, ig, ng, GroupInfo, true)  
    !IF(TranCntl%lExptrsf) CALL SetExpTrsfXs(Core, Fxr, xst1g, iz, ig, GroupInfo, TranInfo, TranCntl, nTracerCntl, PE)  
     
    DO ipin = 1, nxy
      FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
      icel = Pin(ipin)%Cell(iz)
      nlocalFsr = CellInfo(icel)%nFsr
      LocalResidual = 0; localsrc = 0
  !    !Current
      DO j = 1, nbd
        LocalResidual = LocalResidual + (MocJout(2, j, ipin, iz, ig) - MocJout(1, j, ipin, iz, ig))
      ENDDO
      
      DO j = 1, nlocalFsr
        ireg = FsrIdxSt + j - 1
        vol = CellInfo(icel)%vol(j)
        LocalResidual = LocalResidual + vol * phis(ireg, iz, ig) * xst1g(ireg)
      ENDDO
      Res(ipin) = Res(ipin) - LocalResidual      
    ENDDO
  ENDDO  !Group Sweep
  
  DO ipin = 1, nxy
    TranMocResidualError = TranMocResidualError + Res(ipin) * Res(ipin)
  ENDDO
ENDDO

#ifdef MPI_ENV
CALL REDUCE(TranMocResidualError, temp, PE%MPI_RTMASTER_COMM, .TRUE.)
TranMocResidualError = temp
CALL REDUCE(SrcSum, temp, PE%MPI_RTMASTER_COMM, .TRUE.)
SrcSum = temp
#endif
TranMocResidualError = TranMocResidualError / SrcSum
TranMocResidualError = SQRT(TranMocResidualError)

DEALLOCATE(Res)


NULLIFY(Phis, Psi, MocJout, AxSrc, AxPxs)
NULLIFY(Prec, TranPhi, Pin, CellInfo, Fxr)

END FUNCTION



