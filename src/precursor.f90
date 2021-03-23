#include <defines.h>
SUBROUTINE InitPrecursor(Core, FmInfo, CmInfo, TranInfo, GroupInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,           ONLY : CoreInfo_Type,       FmInfo_Type,           CmInfo_Type,        &
                              TranInfo_Type,       PE_Type,               GroupInfo_Type,     &
                              FxrInfo_Type,        Pin_Type,              Cell_Type
USE CNTL,              ONLY : nTracerCntl_Type
USE TRANCMFD_MOD,      ONLY : HomKineticParamGen
USE MPIAxSolver_Mod,   ONLY : InitAxNPrec
USE BenchXs,           ONLY : DnpBetaBen,          DnpBetaDynBen
IMPLICIT NONE
  
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
  

TYPE(FxrInfo_Type), POINTER :: FXR(:, :)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(FxrInfo_Type), POINTER :: myFxr

REAL, POINTER :: psi(:, :), Prec(:, :, :)
REAL, POINTER :: PrecFm(:, :, :), PrecCm(:, :, :)
REAL, POINTER :: PsiC(:, :), PsiFm(:, :)

LOGICAL :: lXsLib

INTEGER :: myzb, myze, myzbf, myzef
INTEGER :: ng, nprec
INTEGER :: nxy, nfsr, nfxr
INTEGER :: iprec, ixy, ifsr, ifxr, icel, iz, izf
INTEGER :: i, j

INTEGER :: FsrIdxSt, FxrIdxSt
INTEGER :: nlocalFxr, nLocalFsr, nFsrInFxr

REAL :: vol, frac
REAL :: beta(1:50), Lambda(1:50), InvLambda(1:50)

CALL HomKineticParamGen(Core, FmInfo, CmInfo, TranInfo, GroupInfo, nTracerCntl, PE)

nfsr = Core%nCoreFsr; nfxr = Core%nCoreFxr
nxy = Core%nxy
myzb = PE%myzb; myze = PE%myze
myzbf = PE%myzbf; myzef = PE%myze

ng = GroupInfo%ng; nprec = GroupInfo%nprec
lXsLib = nTracerCntl%lXsLib

Fxr => FmInfo%Fxr; Psi => FmInfo%psi
Pin => Core%Pin; CellInfo => Core%CellInfo

PsiC => CmInfo%PsiC; PsiFm => CmInfo%PsiFm
Prec => FmInfo%Prec; PrecFm => CmInfo%PrecFm; PrecCm => CmInfo%PrecCm

Lambda(1:nprec) = TranInfo%Lambda(1:nprec)
InvLambda(1:nprec) = TranInfo%InvLambda(1:nprec)
DO iz = myzb, myze
  DO ixy = 1, nxy
    FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
    icel = Pin(ixy)%Cell(iz)
    nlocalFxr = CellInfo(icel)%nFxr; nlocalFsr = CellInfo(icel)%nFsr 
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1; nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
      myFxr => Fxr(ifxr, iz)
      IF(lXsLib) THEN
        beta(1:nprec) = myFxr%Beta(1:nprec)
      ELSE
        IF(nTracerCntl%lDynamicBen) THEN
          CALL DnpBetaDynBen(myFxr%imix, TranInfo%fuelTemp(ixy, iz), beta(1:nprec))
        ELSE
          CALL DnpBetaBen(myFxr%imix, beta(1:nprec))
        END IF
      ENDIF
      DO i = 1, nFsrInFxr
        ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
        DO iprec = 1, nPrec
          Prec(iprec, ifsr, iz) = beta(iprec) * InvLambda(iprec) * psi(ifsr, iz)
        ENDDO
      ENDDO
    ENDDO
    
    PrecCm(:, ixy, iz) = 0
    DO j = 1, CellInfo(icel)%nFsr
      ifsr = FsrIdxSt + j - 1
      vol = CellInfo(icel)%vol(j) * Core%hz(iz)
      DO iprec = 1, nPrec
        PrecCm(iprec, ixy, iz) = PrecCm(iprec, ixy, iz) + Prec(iprec, ifsr, iz) * vol
      ENDDO
    ENDDO
    
    IF(.NOT. nTracerCntl%lSubPlane) CYCLE
    DO izf = Core%SubPlaneRange(1, iz), Core%SubPlaneRange(2, iz)
      frac = PsiFm(ixy, izf) / PsiC(ixy, iz)
      DO iprec = 1, nprec
        PrecFm(iprec, ixy, izf) = PrecCm(iprec, ixy, iz) * frac
      ENDDO
    ENDDO
  ENDDO
ENDDO

IF(nTracerCntl%l3dim) THEN
  CALL InitAxNPrec(CmInfo%PinXs, TranInfo, PE)
ENDIF
!
!Fxr => FmInfo%Fxr; Psi => FmInfo%psi
!Pin => Core%Pin; CellInfo => Core%CellInfo
!
!PsiC => CmInfo%PsiC; PsiFm => CmInfo%PsiFm
!Prec => FmInfo%Prec; PrecFm => CmInfo%PrecFm; PrecCm => CmInfo%PrecCm
NULLIFY(Fxr, Psi, Pin, CellInfo)
NULLIFY(PsiC, PsiFm, Prec, PrecFm, PrecCm)

END SUBROUTINE

SUBROUTINE InitFmPrecursor(Core, FmInfo, TranInfo, nTracerCntl, PE, ng, nprec)
USE PARAM
USE TYPEDEF,           ONLY : CoreInfo_Type,       FmInfo_Type,         TranInfo_Type,       PE_Type,   &
                              FxrInfo_Type,        Pin_Type,            Cell_Type
USE CNTL,              ONLY : nTracerCntl_Type
USE TRANCMFD_MOD,      ONLY : HomKineticParamGen
USE MPIAxSolver_Mod,   ONLY : InitAxNPrec
USE BenchXs,           ONLY : DnpBetaBen,          DnpBetaDynBen
IMPLICIT NONE
  
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
INTEGER :: ng, nprec
  

TYPE(FxrInfo_Type), POINTER :: FXR(:, :)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(FxrInfo_Type), POINTER :: myFxr

REAL, POINTER :: psi(:, :), Prec(:, :, :)

LOGICAL :: lXsLib

INTEGER :: myzb, myze, myzbf, myzef
INTEGER :: nxy, nfsr, nfxr
INTEGER :: iprec, ixy, ifsr, ifxr, icel, iz, izf
INTEGER :: i, j

INTEGER :: FsrIdxSt, FxrIdxSt
INTEGER :: nlocalFxr, nLocalFsr, nFsrInFxr

REAL :: vol, frac
REAL :: beta(nprec), Lambda(nprec), InvLambda(nprec)

nfsr = Core%nCoreFsr; nfxr = Core%nCoreFxr
nxy = Core%nxy
myzb = PE%myzb; myze = PE%myze

lXsLib = nTracerCntl%lXsLib

Fxr => FmInfo%Fxr; Psi => FmInfo%psi
Pin => Core%Pin; CellInfo => Core%CellInfo

Prec => FmInfo%Prec

Lambda(1:nprec) = TranInfo%Lambda(1:nprec)
InvLambda(1:nprec) = TranInfo%InvLambda(1:nprec)
DO iz = myzb, myze
  DO ixy = 1, nxy
    FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
    icel = Pin(ixy)%Cell(iz)
    nlocalFxr = CellInfo(icel)%nFxr; nlocalFsr = CellInfo(icel)%nFsr 
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1; nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
      myFxr => Fxr(ifxr, iz)
      IF(lXsLib) THEN
        beta(1:nprec) = myFxr%Beta(1:nprec)
      ELSE
        IF(nTracerCntl%lDynamicBen) THEN
          CALL DnpBetaDynBen(myFxr%imix, TranInfo%fuelTemp(ixy, iz), beta(1:nprec))
        ELSE
          CALL DnpBetaBen(myFxr%imix, beta(1:nprec))
        END IF
      ENDIF
      DO i = 1, nFsrInFxr
        ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
        DO iprec = 1, nPrec
          Prec(iprec, ifsr, iz) = beta(iprec) * InvLambda(iprec) * psi(ifsr, iz)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO

NULLIFY(Fxr, Psi, Pin, CellInfo)
NULLIFY(Prec)

END SUBROUTINE

SUBROUTINE UpdtPrec(Core, FmInfo, CmInfo, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,             ONLY : CoreInfo_Type,          FmInfo_Type,           CmInfo_Type,         &
                                TranInfo_Type,          GroupInfo_Type,        TranCntl_Type,       &
                                PE_Type,                                                            &
                                FxrInfo_Type,           Pin_Type,               Cell_Type
USE CNTL,                ONLY : nTracerCntl_Type
USE MPIAxSolver_Mod,     ONLY : UpdtAxNPrecSrc
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

REAL, POINTER :: Prec(:, :, :)
REAL, POINTER :: Psi(:, :), TranPsi(:, :), TranPsid(:, :)
REAL, POINTER :: Omegam(:, :, :), Omega0(:, :, :), Omegap(:, :, :)

INTEGER :: myzb, myze, myzbf, myzef, nxy, nfxr, nfsr, nprec
INTEGER :: FxrIdxSt, FsrIdxst, nLocalFxr, nFsrInFxr, NowStep
INTEGER :: ixy, iz, iz0, icel, ifxr, ifsr, iprec
INTEGER :: i, j

REAL :: DelT, reigv, vol, kappa(500), lambda(500)
REAL, POINTER :: Temp(:, :, :)

NowStep = TranCntl%NowStep
Delt = TranCntl%DelT(NowStep)

myzb = PE%myzb; myze = PE%myze
myzbf = PE%myzbf; myzef = PE%myzef
nxy = Core%nxy; nfxr = Core%nCoreFxr; nfsr = Core%nCoreFsr

Fxr => FmInfo%Fxr
Pin => Core%Pin; CellInfo => Core%CellInfo
nprec = TranInfo%nprec
DO iprec = 1, nprec
  lambda(iprec) = TranInfo%lambda(iprec)
  kappa(iprec) = exp(-delt * lambda(iprec))
ENDDO

Prec => FmInfo%Prec
Psi => FmInfo%Psi; TranPsi => FmInfo%TranPsi; TranPsid => FmInfo%TranPsid
Omegam => TranInfo%FxrOmegam; Omega0 => TranInfo%FxrOmega0; Omegap => TranInfo%FxrOmegap

DO iz = myzb, myze
  DO ixy = 1, nxy
    FxrIdxSt = Pin(ixy)%FxrIdxSt; FsrIdxSt = Pin(ixy)%FsrIdxSt
    icel = Pin(ixy)%Cell(iz); nLocalFxr = CellInfo(icel)%nFxr
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      IF(.NOT. Fxr(ifxr, iz)%lfuel) CYCLE
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)   
      DO i = 1, nFsrInFxr
        ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
        !Prvious time step Contribution
        DO iprec = 1, nPrec
          Prec(iprec, ifsr, iz) = kappa(iprec) * Prec(iprec, ifsr, iz)
          Prec(iprec, ifsr, iz) = Prec(iprec, ifsr, iz) + Omegam(iprec, ifxr, iz) * TranPsid(ifsr, iz)
          Prec(iprec, ifsr, iz) = Prec(iprec, ifsr, iz) + Omega0(iprec, ifxr, iz) * TranPsi(ifsr, iz)
          Prec(iprec, ifsr, iz) = Prec(iprec, ifsr, iz) + Omegap(iprec, ifxr, iz) * Psi(ifsr, iz)
          
        ENDDO
      ENDDO
    ENDDO
    !!CmInfo%PrecCm(:, ixy, iz) = 0
    !temp(:, ixy, iz) = 0
    !DO j = 1, CellInfo(icel)%nFsr
    !  ifsr = FsrIdxSt + j - 1
    !  vol = CellInfo(icel)%vol(j)
    !  DO iprec = 1, nPrec
    !    temp(iprec, ixy, iz) = temp(iprec, ixy, iz) + Prec(iprec, ifsr, iz) * vol
    !  ENDDO
    !ENDDO   
    ! 
    
  ENDDO
ENDDO
!
!
Prec => CMInfo%PrecCm
Psi => CmInfo%PsiC; TranPsi => CmInfo%TranPsiCm; TranPsid => CmInfo%TranPsiCmd
Omegam => TranInfo%CellOmegam; Omega0 => TranInfo%CellOmega0; Omegap => TranInfo%CellOmegap
DO iz = myzb, myze
  DO ixy = 1, nxy
    DO iprec = 1, nprec
      Prec(iprec, ixy, iz) = kappa(iprec) * Prec(iprec, ixy, iz)
      Prec(iprec, ixy, iz) = Prec(iprec, ixy, iz) + Omegam(iprec, ixy, iz) * TranPsid(ixy, iz)
      Prec(iprec, ixy, iz) = Prec(iprec, ixy, iz) + Omega0(iprec, ixy, iz) * TranPsi(ixy, iz)
      Prec(iprec, ixy, iz) = Prec(iprec, ixy, iz) + Omegap(iprec, ixy, iz) * Psi(ixy, iz)
      !temp(iprec, ixy, iz) = Temp(iprec, ixy, iz) - Prec(iprec, ixy, iz)
    ENDDO
  ENDDO  
ENDDO

IF(nTracerCntl%lSubPlane) THEN
  Prec => CMInfo%PrecFm
  Psi => CmInfo%PsiFm; TranPsi => CmInfo%TranPsiFm; TranPsid => CmInfo%TranPsiFmd  
  DO iz = myzbf, myzef
    iz0 = Core%SubPlaneMap(iz)
    DO ixy = 1, nxy
      DO iprec = 1, nprec
        Prec(iprec, ixy, iz) = kappa(iprec) * Prec(iprec, ixy, iz)
        Prec(iprec, ixy, iz) = Prec(iprec, ixy, iz) + Omegam(iprec, ixy, iz0) * TranPsid(ixy, iz)
        Prec(iprec, ixy, iz) = Prec(iprec, ixy, iz) + Omega0(iprec, ixy, iz0) * TranPsi(ixy, iz)
        Prec(iprec, ixy, iz) = Prec(iprec, ixy, iz) + Omegap(iprec, ixy, iz0) * Psi(ixy, iz)
      ENDDO
    ENDDO  
  ENDDO
ENDIF

IF(nTracerCntl%l3dim) CALL UpdtAxNPrecSrc(TranInfo, TranCntl, PE)

NULLIFY(Fxr, Pin, CellInfo)
NULLIFY(Prec, Psi, TranPsi, TranPsid)
NULLIFY(Omegam, Omega0, Omegap)

END SUBROUTINE

