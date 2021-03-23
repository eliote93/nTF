#include <defines.h>
SUBROUTINE DepletionStep(Core, FmInfo, DeplVars, DeplLib, ThInfo, GroupInfo, nTracerCntl, DeplCntl, PE)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,       FmInfo_Type,       GroupInfo_Type,    &
                             PE_Type,             ThInfo_Type,                          &
                             FxrInfo_Type,        Pin_Type,          Cell_Type
USE DeplType,         ONLY : DeplLib_Type,        DeplVars_Type,     DeplCntl_Type,     &
                             DeplXs_Type,                                               &
                             AllocDeplXs
USE CNTL,             ONLY : nTracerCntl_Type
USE Depl_mod,         ONLY : FluxNormalizeFactor, SaveFxrIsoInfo,                       &
                             MakeDeplXs1g       , ConstDeplVas,      FxrBurnUp,         &
                             UpdateDeplFxrInfo,   SetLocalBurnup,                       &
                             ConstDeplVas_QD,     SaveDeplXs_GD
USE GdDepl_mod,       ONLY : GdFxrBurnUp,         GdXsFunction,      UpdateDeplGdFxrInfo 
USE BasicOperation,   ONLY : MULTI_VA
USE MOC_Mod,          ONLY : FxrAvgPhi
USE TH_Mod,           ONLY : GetPinFuelTemp
USE OMP_LIB
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) ::  FmInfo
TYPE(DeplVars_Type) :: DeplVars(nThreadMax)
TYPE(DeplLib_Type) :: DeplLib(nThreadMax)
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(DeplCntl_Type) :: DeplCntl
TYPE(PE_TYPE) :: PE


TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
!TYPE(FxrInfo_Type), POINTER :: myFXR
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(DeplXs_Type), SAVE :: DeplXS(nThreadMax)

REAL, POINTER :: phis(:, :, :), SpecConv(:)

REAL, POINTER :: AvgPhi(:)
REAL :: NormFactor, BurnUpTime

INTEGER :: nxy, nFxr, nFsr, ntiso, myzb, myze, ng
INTEGER :: nIsoLib, nIsoDepl
INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nLocalFsr, nFsrInFxr
INTEGER :: iz, ipin, icel, ifxr, tid
INTEGER :: i, j

LOGICAL :: lPredictStep
LOGICAL :: lCorrectStep
LOGICAL :: lCritSpec
LOGICAL :: lHighOd, lGdQuadDepl


lPredictStep = DeplCntl%lPredict
lCorrectStep = .NOT. lPredictStep
lCritSpec = nTracerCntl%lCritSpec
lHighOd = .TRUE.
Pin => Core%Pin; CellInfo => Core%CellInfo;
Fxr => FmInfo%Fxr; Phis => FmInfo%Phis

ng = GroupInfo%ng
nxy = Core%nxy;  ntiso = GroupInfo%ntiso
nFsr = Core%nCoreFsr; nFxr = Core%nCoreFxr

myzb = PE%myzb; myze = PE%myze

BurnUpTime = DeplCntl%Tsec
NormFactor = FluxNormalizeFactor(Core, FmInfo, GroupInfo, DeplCntl%PowerCore, lCritSpec, TRUE, PE)
NormFactor = NormFactor * nTracerCntl%PowerLevel
CALL SetLocalBurnup(Core, FmInfo%FXR, FmInfo%Power, NormFactor, DeplCntl%Tsec, lCorrectStep, PE)
nIsoLib = GroupInfo%ntiso; nIsoDepl = DeplLib(1)%nIsoDep

IF(DeplCntl%NowStep .LE. 2) lHighOd = .FALSE.
DO j = 1, PE%nDeplThread
  IF(.NOT. DeplXS(j)%lAlloc) CALL AllocDeplXs(DeplXs(j), nIsoDepl, ng)
  DeplXs(j)%tid = j
ENDDO
!ALLOCATE(AvgPhi(ng))

IF(lCritSpec) SpecConv => FmInfo%SpecConv

!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(PE%nDeplThread) 
tid = 1

DO iz = myzb, myze
  
!$OMP  PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(FsrIdxSt, FxrIdxSt, nLocalFxr, nLocalFsr, nFsrInFxr, ipin, icel, ifxr, i, j, tid, lHighOd, lGdQuadDepl)
!$  tid = omp_get_thread_num()+1
!$OMP DO ORDERED SCHEDULE(DYNAMIC)
  DO ipin = 1, nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
    icel = Pin(ipin)%Cell(iz); 
    nlocalFxr = CellInfo(icel)%nFxr; nlocalFsr = CellInfo(icel)%nFsr
    !Save Isotope list and number density prior to run depletion calculation
    IF(lPredictStep) CALL SaveFxrIsoInfo(Fxr(FxrIdxSt : FxrIdxst + nLocalFxr - 1, iz), nLocalFxr) 
     
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      IF(.NOT. Fxr(ifxr, iz)%lDepl) CYCLE
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)    

      DeplXs(tid)%AvgPhi = FxrAvgPhi(Core, Fxr, Phis, ipin, j, iz, ng, PE)
      !Conversion Critical Spectrum
      IF(lCritSpec) CALL MULTI_VA(SpecConv(1:ng),DeplXs(tid)%AvgPhi, ng)
      !Calculate 1G XS for Depletion Calculation     
      CALL MakeDeplXs1g(Fxr, ifxr, ng, DeplXs(tid), NormFactor, Core, ipin, iz, GroupInfo, PE, DeplCntl)
  
      lHighOd = .TRUE.
      lGdQuadDepl = .FALSE.
      IF(DeplCntl%NowStep .LE. 2) lHighOd = .FALSE.
#ifdef GdHO   
      IF(.NOT. Fxr(ifxr, iz)%lGd) lHighOd = .FALSE.
#endif
#define lQuadFtn
#ifdef lQuadFtn
      lHighOd = .FALSE.
      IF(Fxr(ifxr, iz)%lGd) THEN
        lGdQuadDepl = .TRUE.
        IF(DeplCntl%NowStep .LE. 2) lGdQuadDepl = .FALSE.
        !IF(lGdQuadDepl .AND. .NOT. lCorrectStep) lHighOd = .TRUE.
      ENDIF
#endif
!      lHighOd = .FALSE.
!      lGdQuadDepl = .FALSE.
      IF(lHighOd) THEN
        CALL ConstDeplVas_QD(DeplVars(tid),  DeplXs(tid), Fxr(ifxr, iz), nIsoLib, nIsoDepl, DeplCntl)
        CALL SaveDeplXs_GD(DeplXs(tid), Fxr(ifxr, iz), lCorrectStep)
        CALL FxrBurnUp(DeplVars(tid), DeplLib(tid), DeplCntl)
        CALL UpdateDeplFxrInfo(Fxr(ifxr, iz), DeplVars(tid), GroupInfo, FALSE, DeplCntl%lXeDyn)
      ELSEIF(lGdQuadDepl) THEN
        CALL ConstDeplVas(DeplVars(tid),  DeplXs(tid), Fxr(ifxr, iz), nIsoLib, nIsoDepl)
        IF(DeplCntl%NowStep .GT. 2 .AND. lCorrectStep) THEN
          CALL GdXsFunction(DeplVars(tid), DeplXs(tid), Fxr(ifxr, iz), lCorrectStep)
          CALL GdFxrBurnUp(DeplVars(tid), DeplLib(tid), DeplCntl)
          CALL SaveDeplXs_GD(DeplXs(tid), Fxr(ifxr, iz), lCorrectStep)
          CALL UpdateDeplGdFxrInfo(Fxr(ifxr, iz), DeplVars(tid), GroupInfo, lCorrectStep, DeplCntl%lXeDyn)     
        ELSE
          CALL FxrBurnUp(DeplVars(tid), DeplLib(tid), DeplCntl)
          CALL SaveDeplXs_GD(DeplXs(tid), Fxr(ifxr, iz), lCorrectStep)
          CALL UpdateDeplFxrInfo(Fxr(ifxr, iz), DeplVars(tid), GroupInfo, lCorrectStep, DeplCntl%lXeDyn)
        ENDIF        
      ELSE
        CALL ConstDeplVas(DeplVars(tid),  DeplXs(tid), Fxr(ifxr, iz), nIsoLib, nIsoDepl)
#ifdef GdHO      
        IF(Fxr(ifxr, iz)%lGd) CALL SaveDeplXs_GD(DeplXs(tid), Fxr(ifxr, iz), lCorrectStep)
#else
        CALL SaveDeplXs_GD(DeplXs(tid), Fxr(ifxr, iz), lCorrectStep)
#endif
        CALL FxrBurnUp(DeplVars(tid), DeplLib(tid), DeplCntl)
        CALL UpdateDeplFxrInfo(Fxr(ifxr, iz), DeplVars(tid), GroupInfo, lCorrectStep, DeplCntl%lXeDyn)  
      ENDIF
      
      IF(Fxr(ifxr, iz)%lGd .AND. DeplCntl%PC_OPT .EQ. 1) CALL GdPostCorrection(Fxr(ifxr, iz), DeplCntl, ifxr)

    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL
ENDDO

IF(lCritSpec) NULLIFY(SpecConv)

ENDSUBROUTINE

SUBROUTINE SaveDeplXs1G(myFxr, DeplXs, lCorrectStep)
USE PARAM
USE TYPEDEF,          ONLY : FxrInfo_Type,        DeplGd_TYPE
USE DeplType,         ONLY : DeplXs_Type
TYPE(FxrInfo_Type) :: myFxr
TYPE(DeplXs_Type) :: DeplXS
LOGICAL :: lCorrectStep
END SUBROUTINE

SUBROUTINE FxrBurnUp(DeplVars, DeplLib, DeplCntl)
USE PARAM
USE DeplType,    ONLY : DeplVars_Type,      DeplLib_Type,        DeplCntl_Type,          &
                        MatExp_Type
USE Depl_mod,    ONLY : MatExp,             AllocMatExpVec,      InitDeplXS,             &
                        TranDepXS2DeplLib,  CopyIsoNumVector,    MakeDeplMat
USE MatExp_mod,  ONLY : MatExpSolver
IMPLICIT NONE

TYPE(DeplVars_Type) :: DeplVars
TYPE(DeplLib_Type) :: DeplLib
TYPE(DeplCntl_Type) :: DeplCntl

REAL, POINTER :: IsoNum(:), BurnUpXs(:, :)
REAL :: BurnUpTime, Phi1g
INTEGER :: nIsoDep

INTEGER :: Tid

Tid = DeplVars%tid
nIsoDep = DeplLib%nIsoDep
MatExp(Tid)%Mat => DeplVars%Dmat; MatExp(Tid)%nIsoDepl = nIsoDep
MatExp(TID)%Solver = DeplCntl%Solver
IsoNum => DeplVars%IsoNum; BurnUpXs => DeplVars%BurnUpXS

IF(.NOT. MatExp(Tid)%lAllocVec) THEN
  CALL AllocMatExpVec(MatExp(tid), DeplLib%nIsoDep)
ENDIF

CALL InitDeplXS(DeplLib)

CALL TranDepXS2DeplLib(DeplLib, BurnUpXs, IsoNum, nIsoDep)
CALL CopyIsoNumVector(MatExp(tid)%Viso0, IsoNum, nIsoDep, epsm30)
BurnUpTime = DeplCntl%Tsec; phi1g = DeplVars%phi1g * 1.0E-24_8
CALL MakeDeplMat(MatExp(tid)%Mat, DeplLib, Phi1g, BurnUpTime) 
!--- BYS edit for PKB 14/09/17
!CALL DeplMatOut(MatExp(tid))
!--- BYS edit for PKB 14/09/22
CALL MatExpSolver(MatExp(tid))
CALL CopyIsoNumVector(IsoNum, MatExp(tid)%VIso, nIsoDep, epsm30)
NULLIFY(MatExp(tid)%Mat)
NULLIFY(IsoNum, BurnUpXs)

END SUBROUTINE

SUBROUTINE InitDeplXS(DeplLib)
USE PARAM
USE DeplType,    ONLY : DeplLib_Type,                                  &
                        AtomKind,          ISOTOPEKIND,        STATEKIND  
IMPLICIT NONE
TYPE(DeplLib_Type) :: DeplLib

TYPE(AtomKind), POINTER :: Lib0(:), Lib1(:)
TYPE(STATEKIND), POINTER :: myStat0, myStat1

INTEGER :: nIsoDep, nAtomDep
INTEGER :: i, j, k, jb, je

Lib0 => DeplLib%AtomLib0; Lib1 => DeplLib%AtomLib1

nIsoDep = DeplLib%nIsoDep
nAtomDep = DeplLib%nAtomDep
DO i = 1, nAtomDep
  jb = Lib0(i)%ib; je =Lib0(i)%ie
  DO j = jb, je
    DO k = 1, Lib0(i)%A(j)%nSTATE
      IF(Lib0(i)%A(j)%Stat(k)%IMAT .EQ. 0) CYCLE
      myStat0 => Lib0(i)%A(j)%Stat(k); myStat1 => Lib1(i)%A(j)%Stat(k)
      myStat1%XS(0:6) = myStat0%XS(0:6)
    ENDDO
  ENDDO
ENDDO

IF(ASSOCIATED(myStat0)) NULLIFY(myStat0)
IF(ASSOCIATED(myStat1)) NULLIFY(myStat1)
END SUBROUTINE

SUBROUTINE TranDepXS2DeplLib(DeplLib, BurnUpXs, IsoNum, nIsoDep)
USE PARAM
USE DeplType,     ONLY : DeplLib_Type,  AtomKind,                &
                         STATEKIND,     ISOTOPEKIND           

TYPE(DeplLib_Type) :: DeplLib
REAL :: BurnUpXs(4, nIsoDep)
REAL :: IsoNum(nIsoDep)
INTEGER :: nIsoDep

TYPE(AtomKind), POINTER :: Lib1(:)
TYPE(StateKind), POINTER :: MyStat
TYPE(IsotopeKind), POINTER :: MyIso

INTEGER, POINTER :: MapMatId2ANum(:), MapMatId2IsoWt(:), MapMatId2State(:)

REAL :: XS1, XS2, XS3, XS4, XS5, XS6
REAL :: FR

INTEGER :: nAtomDep
INTEGER :: I, J, JB, JE
INTEGER :: IZ, IA, IST, IGRP


nAtomDep = DeplLib%nAtomDep

Lib1 => DeplLib%AtomLib1
MapMatId2ANum => DeplLib%MapMatId2ANum
MapMatId2IsoWt => DeplLib%MapMatId2IsoWt
MapMatId2State => DeplLib%MapMatId2State


DO I = 1, nIsoDep
  IF(IsoNum(I) .LT. EPSM30) CYCLE
  IF((BurnUpXs(1, I) + BurnUpXs(3, I)) .LT. EPSM30) CYCLE
  IZ = MapMatId2Anum(I)
  IA = MapMatId2IsoWt(I)
  IST = MapMatId2State(I)
  MyIso => Lib1(IZ)%A(IA); MyStat => MyIso%STAT(IST)
  IGRP = MyStat%IGRP
  XS1 = MyStat%XS(1);  XS2 = MyStat%XS(2);  XS3 = MyStat%XS(3)
  XS4 = MyStat%XS(4);  XS5 = MyStat%XS(5);  XS6 = MyStat%XS(6)
  IF((XS1+XS5) .GT. epsm30 .AND. BurnUpXs(1,I) .GT. EPSM30) THEN
    
    IF(IGRP .EQ. 2) THEN
      FR = XS1 / (XS1 + XS5)
      MyStat%XS(1) = FR * (BurnUpXs(1, I) - BurnUpXs(2,I))
      MyStat%XS(3) = BurnUpXS(4, I)     !(N,3N) (Z,A,0)->(Z,A-2,0) FOR ACTINIDE (J=3)
      MyStat%XS(4) = BurnUpXS(2,I)
      MyStat%XS(5) = (ONE-FR)*(BurnUpXS(1,I) - BurnUpXS(2,I))      
    ELSE
      !MyStat%XS(1) = FR * BurnUpXS(1,I)
      !MyStat%XS(5) = (ONE-FR) * BurnUpXS(1,I)   
      FR = XS1 / (XS1 + XS3 + XS4 + XS5); MyStat%XS(1) = FR * BurnUpXS(1,I) 
      FR = XS3 / (XS1 + XS3 + XS4 + XS5); MyStat%XS(3) = FR * BurnUpXS(1,I) 
      FR = XS4 / (XS1 + XS3 + XS4 + XS5); MyStat%XS(4) = FR * BurnUpXS(1,I) 
      FR = XS5 / (XS1 + XS3 + XS4 + XS5); MyStat%XS(5) = FR * BurnUpXS(1,I)   
    ENDIF
  ENDIF
  XS2 = MyStat%XS(2);  XS6 = MyStat%XS(6)
  IF((XS2 + XS6) .GT. EPSM30 .AND. BurnUpXs(3, I) .GT. EPSM30) THEN
    FR = XS2 / (XS2 + myStat%XS(6))
    myStat%XS(2) = FR * BurnUpXS(3,I)
    myStat%XS(6) = (ONE - FR) * BurnUpXS(3,I)  
  ENDIF
  MyStat%XS(0)=SUM(MyStat%XS(1:6))
  NULLIFY(MyStat); NULLIFY(MyIso)
ENDDO

NULLIFY(Lib1)
NULLIFY(MapMatId2Anum)
NULLIFY(MapMAtId2Isowt)
NULLIFY(MapMatId2State)
END SUBROUTINE

SUBROUTINE CopyIsoNumVector(V1, V2, n, torr)
IMPLICIT NONE
REAL :: V1(n), V2(n)
INTEGER :: n
REAL :: torr

INTEGER :: i

DO i = 1, n
  IF(V2(i) .LT. Torr) THEN
    V1(i) = 0._8
    CYCLE
  ENDIF
  V1(i) = V2(i)
ENDDO
END SUBROUTINE

SUBROUTINE MakeDeplMat(DMat, DeplLib, PHI, BurnUpTime)
USE PARAM
USE DEPLTYPE,      ONLY : Mat_TYPE,         DeplLib_Type,                    &
                          AtomKind,         STATEKIND,        ISOTOPEKIND,   &
                          FisYield_Type      
USE BasicOperation, ONLY : CP_CA,           MULTI_CA
IMPLICIT NONE
TYPE(Mat_Type) :: DMat
TYPE(DeplLib_Type) :: DeplLib
REAL :: PHI, BurnUpTime

TYPE(AtomKind), POINTER :: Lib1(:)
TYPE(StateKind), POINTER :: MyStat, MyStatFr
TYPE(IsotopeKind), POINTER :: MyIso
TYPE(FisYield_Type), POINTER :: FisYield(:)


REAL, POINTER :: DIAG(:), OffDiag(:,:)              !Diag and Off Diag Elements for Sparse Matrix
INTEGER, POINTER :: nlmnt(:), lmntIdx(:,:)          !#of element for row, non-zero element index
INTEGER :: nmaxoffdiag 

REAL :: Y

INTEGER :: nIsoDep, nAtomDep
INTEGER :: I, J, K, M, JB, JE
INTEGER :: IM, IP, IFR, ITO, IZFR, IAFR

Lib1 => DeplLib%AtomLib1
FisYield => DeplLib%FisYield

nIsoDep = DeplLib%nIsoDep; nAtomDep = DeplLib%nAtomDep

Diag => DMAT%DIAG; OffDiag => DMAT%OffDiag
nLmnt => DMAT%nLmnt; lmntIdx => DMAT%LmntIdx
nMaxOffDiag = DMAT%nMaxOffDiag

CALL CP_CA(Diag, 0._8, nIsoDep)
CALL CP_CA(OffDiag, 0._8, nMaxOffDiag, nIsoDep)

DO I = 1, nAtomDep
  DO J = Lib1(I)%IB, Lib1(I)%IE
    MyIso => Lib1(I)%A(J)
    DO K = 0, MyIso%NSTATE
      IF(MyIso%STAT(K)%IMAT .EQ. 0) CYCLE
      MyStat => MyIso%STAT(K)
      IM = MyStat%IMAT
      Diag(IM) = -PHI * MyStat%XS(0) - MyStat%Rambda   ! LOSS
      DO M = 1, MyStat%NTO1    !  SOURCE FROM DECAY CHAIN
        IFR = MyStat%ITO1(1, M); ITO = MyStat%ITO1(2, M); IP = MyStat%ITO1(3, M)
        OffDiag(IP, ITO) = OffDiag(IP, ITO) + MyStat%FRAC(IFR)
      ENDDO

      DO M = 1, MyStat%NTO2   !  SOURCE FROM NUCLEAR REACTION
        IFR = MyStat%ITO2(1, M); ITO = MyStat%ITO2(2, M); IP = MyStat%ITO2(3, M)
        OffDiag(IP, ITO) = OffDiag(IP, ITO) + PHI * MyStat%XS(IFR)
      ENDDO
      
      IF(MySTAT%IGRP .EQ. 3) THEN
        DO M = 1, MyStat%NFR3  !   SOURCE FROM FISSION
          IFR = MyStat%IFR3(1, M); IP = MyStat%IFR3(2, M); Y = MyStat%Y(IFR)
          IZFR = FisYield(IFR)%AtomNum; IAFR = FisYield(IFR)%AtomWeight
          MyStatFr => Lib1(IZFR)%A(IAFR)%STAT(0)
          OffDiag(IP, IM) = OffDiag(IP, IM) + Y * PHI * MyStatFr%XS(4)
          !IZFRIAFR
        ENDDO
      ENDIF
    ENDDO
  ENDDO
ENDDO

IF(ASSOCIATED(MyStat)) NULLIFY(MyStat)
IF(ASSOCIATED(MyStatFr)) NULLIFY(MyStatFr)

!MULTIFY Time Step Size to Depleition Matrix

CALL MULTI_CA(BurnUpTime, Diag(1:NISODEP), NISODEP)
DO I = 1, NISODEP
  M = nLmnt(I)
  CALL MULTI_CA(BurnUpTime, OffDiag(1:M, I), M)
ENDDO

NULLIFY(FisYield);  NULLIFY(Lib1)
NULLIFY(Diag); NULLIFY(OffDiag)
NULLIFY(nLmnt); NULLIFY(LmntIdx)
END SUBROUTINE

SUBROUTINE UpdateDeplFxrInfo(Fxr, DeplVars, GroupInfo, lCorrectStep, lXeDyn)
USE PARAM
USE TYPEDEF,            ONLY : FxrInfo_Type,       GroupInfo_Type
USE DeplType,           ONLY : DeplVars_Type
USE BasicOperation,     ONLY : MULTI_CA,  CP_CA
USE nuclidmap_mod,      ONLY : iposiso,   PndCrit, nuclide
USE XsUtil_mod,         ONLY : SetXeDynEnv
USE ieee_arithmetic
IMPLICIT NONE
TYPE(FxrInfo_Type) :: Fxr
TYPE(DeplVars_Type) :: DeplVars
TYPE(GroupInfo_Type) :: GroupInfo
LOGICAL :: lCorrectStep, lXeDyn

INTEGER, POINTER :: MapXs2Dep(:), MapDep2Xs(:), idiso(:)
REAL, POINTER :: IsoNum(:), pnum(:)

INTEGER :: i, j, IdXs, IdDep
INTEGER :: niso, nisoDepl, ntiso 
INTEGER :: niso_wolib, idiso_wolib(100)
REAL ::  pnum_wolib(100)

MapXs2Dep => DeplVars%MapXs2Dep
MapDep2Xs => DeplVars%MapDep2Xs
IsoNum => DeplVars%IsoNum
nisodepl = DeplVars%nIsoDepl
IdIso => Fxr%IdIso
pnum => Fxr%pnum

IF(lCorrectStep) THEN
  CALL MULTI_CA(0.5_8, IsoNum(1:nIsoDepl), nIsoDepl)
  IF(Fxr%l_pnum_all) THEN ! 16/02/11 Depletion timestep bug fixed
    niso = Fxr%nIso_depl
    DO i = 1, nIsoDepl
      IsoNum(i) = IsoNum(i) + 0.5_8 * Fxr%pnum_all(i)
    ENDDO
  ELSE
    niso = Fxr%nIso_depl
    DO i = 1, niso
      IdXs = iposiso(IdIso(i))
      IdDep = MapXs2Dep(IdXs)
      IF(IdDep .EQ. 0) CYCLE
      IsoNum(IdDep) = IsoNum(IdDep) + 0.5_8 * pnum(i)
    ENDDO
  ENDIF
ENDIF


niso_wolib=0
DO i = 1,Fxr%nIso
  IdXs = iposiso(IdIso(i))
  IdDep = MapXs2Dep(IdXs)
  IF(IdDep .NE. 0) CYCLE
  niso_wolib=niso_wolib+1
  idiso_wolib(niso_wolib) = IdIso(i)
  pnum_wolib(niso_wolib) = pnum(i)
ENDDO

ntiso = GroupInfo%ntiso

CALL CP_CA(IdIso, 0, ntiso)
CALL CP_CA(pnum, 0._8, ntiso)
niso = 0
DO i = 1, nIsoDepl
  j = MapDep2Xs(i)
  IF(j .EQ. 0) CYCLE
  IF(abs(IsoNum(i)) .LT. epsm20) CYCLE
  IF(IsoNum(i) .LT. pndcrit(j)) CYCLE
  IF(ieee_is_nan(IsoNum(i))) STOP 'DEPL NaN Wow!'
  niso = niso + 1
  IdXs = nuclide(j)
  pnum(niso) = IsoNum(i)
  idiso(niso) = IdXs
ENDDO

DO i = 1, nIsoDepl ! 16/02/11 Depletion timestep bug fixed
  Fxr%pnum_all(i) = IsoNum(i)
ENDDO


IF(niso_wolib .GT. 0) THEN
  DO i = 1, niso_wolib
    niso = niso + 1
    pnum(niso) = pnum_wolib(i)
    idiso(niso) = idiso_wolib(i)
  ENDDO
ENDIF

Fxr%niso = niso
Fxr%niso_depl = niso

IF(lXeDyn) CALL SetXeDynEnv(IdIso, pnum, Fxr%niso, Fxr%niso_depl, ntiso)

NULLIFY(MapXs2Dep, IsoNum)
NULLIFY(pnum, idiso)
END SUBROUTINE

SUBROUTINE DeplMatOut(MatExp)
USE DEPLTYPE,      ONLY : Mat_TYPE, MatExp_TYPE
IMPLICIT NONE
TYPE(MatExp_Type) :: MatExp
TYPE(Mat_Type), POINTER :: DMat
INTEGER :: nisodepl
INTEGER :: iso, lmnt
INTEGER :: ounit =99

Dmat => MatExp%Mat
nisodepl = MatExp%nisodepl

OPEN(OUNIT, file='depl.mat', STATUS='REPLACE')
WRITE(OUNIT,'(i12)') nisodepl
DO iso = 1, nisodepl
    WRITE(OUNIT,'(1p, e12.5)') Dmat%diag(iso)
ENDDO
DO iso = 1, nisodepl
    WRITE(OUNIT,'(i12)') Dmat%nlmnt(iso)
    DO lmnt = 1, Dmat%nlmnt(iso)
        WRITE(OUNIT,'(1p, e12.5)') Dmat%offdiag(lmnt,iso)
    ENDDO
ENDDO
CLOSE(OUNIT)
END SUBROUTINE