MODULE XsPerturb_mod
USE PARAM
USE TYPEDEF,            ONLY : CoreInfo_Type,           FmInfo_Type,          TranInfo_Type,      &
                               TranCntl_Type,           PE_TYPE,                                  &
                               FxrInfo_Type,            AsyInfo_Type,         XsChange_Type,      &
                               Asy_Type,                Pin_Type,             Cell_Type,          &
                               XsNoise_Type,            XsCntlRod_Type
USE CNTL,               ONLY : nTracerCntl_Type
USE BenchXs,            ONLY : XsBenChange_new,         xsDynBenChange,       XsBenNoise
IMPLICIT NONE

TYPE(FxrInfo_Type), POINTER, PRIVATE :: Fxr(:, :)
TYPE(AsyInfo_Type), PRIVATE, POINTER :: AsyType(:)
TYPE(Asy_Type), PRIVATE, POINTER :: Asy(:)
TYPE(Pin_Type), PRIVATE, POINTER :: Pin(:)
TYPE(Cell_Type), PRIVATE, POINTER :: CellInfo(:)
INTEGER, PRIVATE, POINTER :: CoreMap(:)

TYPE Fxrsave_type
  LOGICAL :: lcusping, lMixtureMixing, lh2o, lRes
  INTEGER :: iso0, iso1, imix, ndim, niso
  INTEGER, ALLOCATABLE :: IdIso(:)
  REAL :: wt, h2ofrac
  REAL, ALLOCATABLE :: pnum(:)
END TYPE

TYPE(XsChange_Type), PRIVATE, POINTER :: XsChangeSave(:)
TYPE(Fxrsave_Type), PRIVATE, POINTER :: FxrSave(:,:)

INTEGER, PRIVATE :: myzb, myze, nxy

INTERFACE XsPerturbation
  MODULE PROCEDURE XsPerturbation_step
  MODULE PROCEDURE XsPerturbation_T
END INTERFACE

public
CONTAINS
SUBROUTINE XsCntlRodMix(Core, FmInfo, TranCntl, nTracerCntl, PE)
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(TranCntl_Type)  :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

TYPE(XsCntlRod_Type), POINTER :: XsCntlRod(:)

REAL :: wt, wt0
INTEGER :: nrodded
INTEGER :: i, j, iasy, izbeg, izend

LOGICAL :: lXsLib
LOGICAL :: lComplete

INTEGER, PARAMETER :: xsp_crit = 0.05_8

XsCntlRod => TranCntl%XsCntlRod
lXsLib = nTracerCntl%lXslib
Fxr => FmInfo%Fxr; AsyType => Core%AsyInfo; CoreMap => Core%CoreMap
Asy => Core%Asy; Pin => Core%Pin; CellInfo => Core%CellInfo
myzb = PE%myzb; myze = PE%myze; nxy = Core%nxy

nrodded = nTracerCntl%nCntlRod
nTracerCntl%lCusping_MPI = .FALSE.

DO i = 1, nrodded
  IF(wt .EQ. 1.) THEN 
    lComplete = .TRUE. 
  ELSE
    lComplete = .FALSE.
  END IF
  wt = XsCntlRod(i)%wt
  wt0 = XsCntlRod(i)%wt
  nTracerCntl%lCusping_MPI = nTracerCntl%lCusping_MPI .OR. XsCntlRod(i)%lCusping
  IF(.NOT. lXsLib) THEN
    IF(TranCntl%lDynamicBen) THEN
      CALL XsDynBenChange(XsCntlRod(i)%Iso0, XsCntlRod(i)%Iso1, XsCntlRod(i)%Isonew, wt0, XsCntlRod(i)%lCusping, XsCntlRod(i)%lCuspingDirection)
    ELSE
      CALL XsBenChange_New(XsCntlRod(i)%Iso0, XsCntlRod(i)%Iso1, XsCntlRod(i)%Isonew, wt0, XsCntlRod(i)%lCusping, XsCntlRod(i)%lCuspingDirection)
    END IF
  END IF
  izbeg = XsCntlRod(i)%izbeg; izend = XsCntlRod(i)%izend
  DO j = 1, XsCntlRod(i)%nasy
    iasy = XsCntlRod(i)%AsyList(j)
    IF(lXsLib) THEN
      CALL AsyXsChangeLib(XsCntlRod(i)%Iso0, XsCntlRod(i)%Iso1, izbeg, izend, iasy, wt, wt0, lComplete, XsCntlRod(i)%lCusping, XsCntlRod(i)%lCuspingDirection)
    ELSE
      CALL AsyXsChangeBen(XsCntlRod(i)%Iso0, XsCntlRod(i)%Iso1, XsCntlRod(i)%Isonew, izbeg, izend, iasy, lComplete)
    ENDIF
  ENDDO
ENDDO

NULLIFY(XsCntlRod, Fxr, AsyType, Asy, Pin, CellInfo, CoreMap)
END SUBROUTINE

SUBROUTINE XsNoisePerturbation(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, PE, T)
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type)  :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
REAL :: T

TYPE(XsNoise_Type), POINTER :: XsNoise(:)

INTEGER :: nNoise
INTEGER :: nLocalFxr, FxrIdxSt
INTEGER :: iasytype, loc_ixy, ixy, icel, iz, ifxr
INTEGER :: inn, j
LOGICAL :: lXsLib

nNoise = TranCntl%nNoise
XsNoise => TranCntl%XsNoise
lXsLib = nTracerCntl%lXsLib

Pin => Core%Pin
CellInfo => Core%CellInfo
Fxr => FmInfo%Fxr


DO inn = 1, nNoise
  IF(XsNoise(inn)%lFirst) THEN
    iasytype = Core%Asy(XsNoise(inn)%ixya)%AsyType
    loc_ixy = Core%AsyInfo(iasytype)%Pin2DIdx(XsNoise(inn)%ix, XsNoise(inn)%iy)
    XsNoise(inn)%ixy = Core%Asy(XsNoise(inn)%ixya)%GlobalPinIdx(loc_ixy)
  END IF
  ixy = XsNoise(inn)%ixy
  FxrIdxSt = Pin(ixy)%FxrIdxSt
  IF(.NOT. lXsLib) THEN
    IF(TranCntl%lDynamicBen) THEN
    ELSE
      DO iz = XsNoise(inn)%izbeg, XsNoise(inn)%izend
        IF(iz .LT. PE%myzb) CYCLE 
        IF(iz .GT. PE%myze) CYCLE
        icel = Pin(ixy)%Cell(iz)
        nLocalFxr = CellInfo(icel)%nFxr
        DO j = 1, nLocalFxr
          ifxr = FxrIdxSt + j - 1
          IF(XsNoise(inn)%iso0 .NE. Fxr(ifxr, iz)%imix0) CYCLE
          CALL XsBenNoise(XsNoise(inn), T)
          Fxr(ifxr, iz)%imix = XsNoise(inn)%isonew
        END DO
      END DO
    END IF
  END IF
END DO

END SUBROUTINE
  
SUBROUTINE XsPerturbation_Step(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, PE)
#ifdef __PGI
USE CUDA_CONST,       ONLY : theta
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type)  :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

TYPE(XsChange_Type), POINTER :: XsChange(:)

REAL :: T, Tprev, Tprev0, Tbeg, Tend
INTEGER :: nChange, NowStep

NowStep = TranCntl%NowStep
T = TranCntl%T(NowStep); Tprev = 0
IF(NowStep .GT. 1) Tprev = TranCntl%T(NowStep - 1)

CALL XsPerturbation_T(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, PE, T, Tprev)

END SUBROUTINE

SUBROUTINE XsPerturbation_T(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, PE, T, Tprev)
#ifdef __PGI
USE CUDA_CONST,       ONLY : theta
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
REAL :: T, Tprev

TYPE(XsChange_Type), POINTER :: XsChange(:)

REAL :: Tbeg, Tend, Tprev0
REAL :: wt, wt0
INTEGER :: nChange
INTEGER :: i, j, iasy, izbeg, izend
LOGICAL :: lXsLib

REAL, PARAMETER :: xsp_crit = 0.0 !5_8

XsChange => TranCntl%XsChange
lXsLib = nTracerCntl%lXslib
Fxr => FmInfo%Fxr
AsyType => Core%AsyInfo
CoreMap => Core%CoreMap
Asy => Core%Asy
Pin => Core%Pin
CellInfo => Core%CellInfo
myzb = PE%myzb
myze = PE%myze
nxy = Core%nxy

nChange = TranCntl%nChange

IF(TranCntl%ImplicitSwitch) THEN
  TranCntl%Theta = TranCntl%Theta0
#ifdef __PGI
  theta = TranCntl%Theta
#endif
  TranCntl%ImplicitSwitch = .FALSE.
END IF

TranCntl%lXsPerturb = .FALSE.
nTracerCntl%lCusping_MPI = .FALSE.
DO i = 1, nChange
  Tbeg = XsChange(i)%Tbeg
  Tend = XsChange(i)%Tend
  IF(Tend .GE. T) CYCLE
  IF(XsChange(i)%lComplete) CYCLE
  XsChange(i)%lComplete = .TRUE.
  XsChange(i)%lStart = .TRUE.
  wt = 1._8; wt0 = 1._8
  nTracerCntl%lCusping_MPI = nTracerCntl%lCusping_MPI .OR. XsChange(i)%lCusping
  IF(.NOT. lXsLib) THEN
    IF(TranCntl%lDynamicBen) THEN
      CALL XsDynBenChange(XsChange(i)%Iso0, XsChange(i)%Iso1, XsChange(i)%Isonew, wt0, XsChange(i)%lCusping, XsChange(i)%lCuspingDirection)
    ELSE
      CALL XsBenChange_New(XsChange(i)%Iso0, XsChange(i)%Iso1, XsChange(i)%Isonew, wt0, XsChange(i)%lCusping, XsChange(i)%lCuspingDirection)
    END IF
  END IF
  izbeg = XsChange(i)%izbeg; izend = XsChange(i)%izend
  DO j = 1, XsChange(i)%nasy
    iasy = XsChange(i)%AsyList(j)
    IF(lXsLib) THEN
      CALL AsyXsChangeLib(XsChange(i)%Iso0, XsChange(i)%Iso1, izbeg, izend, iasy, wt, wt0, .TRUE., XsChange(i)%lCusping, XsChange(i)%lCuspingDirection)
    ELSE
      CALL AsyXsChangeBen(XsChange(i)%Iso0, XsChange(i)%Iso1, XsChange(i)%Isonew, izbeg, izend, iasy, .TRUE.)
    ENDIF
  ENDDO
  IF(abs(wt0 - XsChange(i)%wt) > xsp_crit) TranCntl%lXsPerturb = .TRUE.
  XsChange(i)%wt = wt0
  IF(XsChange(i)%lStepFunc) THEN
    IF(TranCntl%lStepImplicit .AND. TranCntl%TD .LE. 2) THEN
      TranCntl%Theta0 = TranCntl%Theta
      TranCntl%Theta = 1.0
#ifdef __PGI
      theta = TranCntl%Theta
#endif
      TranCntl%ImplicitSwitch = .TRUE.
    END IF
  END IF 
END DO
DO i = 1, nChange
  Tbeg = XsChange(i)%TBeg; Tend =XsChange(i)%Tend
  IF(T .LT. Tbeg .OR. T .GT. Tend) CYCLE
  IF(XsChange(i)%lStepFunc) CYCLE
  XsChange(i)%lStart = .TRUE.
  Tprev0 = Max(Tprev, Tbeg)
  wt = (T - Tprev0) / (Tend - Tprev0)
  wt0 = (T - Tbeg) / (Tend - Tbeg)
  nTracerCntl%lCusping_MPI = nTracerCntl%lCusping_MPI .OR. XsChange(i)%lCusping
  IF(.NOT. lXsLib) THEN
    IF(TranCntl%lDynamicBen) THEN
      CALL XsDynBenChange(XsChange(i)%Iso0, XsChange(i)%Iso1, XsChange(i)%Isonew, wt0, XsChange(i)%lCusping, XsChange(i)%lCuspingDirection)
    ELSE
      CALL XsBenChange_New(XsChange(i)%Iso0, XsChange(i)%Iso1, XsChange(i)%Isonew, wt0, XsChange(i)%lCusping, XsChange(i)%lCuspingDirection)
    END IF
  END IF
  izbeg = XsChange(i)%izbeg; izend = XsChange(i)%izend
  DO j = 1, XsChange(i)%nasy
    iasy = XsChange(i)%AsyList(j)
    IF(lXsLib) THEN
      CALL AsyXsChangeLib(XsChange(i)%Iso0, XsChange(i)%Iso1, izbeg, izend, iasy, wt, wt0, .FALSE., XsChange(i)%lCusping, XsChange(i)%lCuspingDirection)
    ELSE
      CALL AsyXsChangeBen(XsChange(i)%Iso0, XsChange(i)%Iso1, XsChange(i)%Isonew, izbeg, izend, iasy, .FALSE.)
    ENDIF
  ENDDO
  IF(abs(wt0 - XsChange(i)%wt) > xsp_crit) TranCntl%lXsPerturb = .TRUE.
  XsChange(i)%wt = wt0
ENDDO

NULLIFY(XsChange, Fxr, AsyType, Asy, Pin, CellInfo, CoreMap)
CALL XsNoisePerturbation(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, PE, T)

END SUBROUTINE

SUBROUTINE InitXsPerturbation(Core, TranCntl, PE)
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(TranCntl_Type) :: TranCntl
TYPE(PE_Type) :: PE
ALLOCATE(FxrSave(Core%nCoreFxr, PE%myzb:PE%myze))
ALLOCATE(XsChangeSave(TranCntl%nChange))
END SUBROUTINE

SUBROUTINE SaveXsPerturbation(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, PE, T, Tsave)
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
REAL :: T, Tsave

TYPE(XsChange_Type), POINTER :: XsChange(:)
REAL :: tbeg, tend
INTEGER :: nChange, FxrIdxSt, nLocalFxr, ndim
INTEGER :: izbeg, izend, iasy, ixy, ixy0, ifxr, iz, icel, iasytype
INTEGER :: i, j, k

IF(.NOT. nTracerCntl%lXslib) STOP
XsChange => TranCntl%XsChange
nChange = TranCntl%nChange
Fxr => FmInfo%Fxr
AsyType => Core%AsyInfo
CoreMap => Core%CoreMap
Asy => Core%Asy
Pin => Core%Pin
CellInfo => Core%CellInfo
myzb = PE%myzb
myze = PE%myze
nxy = Core%nxy

DO i = 1, nChange
  tbeg = XsChange(i)%tbeg
  tend = XsChange(i)%tend
  IF(T .LT. tbeg .OR. Tsave .GT. tend) CYCLE
  izbeg = XsChange(i)%izbeg
  izend = XsChange(i)%izend

  XsChangeSave(i)%lUse = XsChange(i)%lUse
  XsChangeSave(i)%lstart = XsChange(i)%lstart
  XsChangeSave(i)%lcomplete = XsChange(i)%lcomplete
  XsChangeSave(i)%wt = XsChange(i)%wt
  
  DO j = 1, XsChange(i)%nasy
    iasy = XsChange(i)%AsyList(j)
    iasytype = Core%CoreMap(iasy)
    DO iz = PE%myzb, PE%myze
      IF(iz .GT. izend .OR. iz .LT. izbeg) CYCLE
      DO ixy0 = 1, AsyType(iasytype)%nxy
        ixy = Asy(iasy)%GlobalPinIdx(ixy0)
        icel = Pin(ixy)%Cell(iz)
        FxrIdxSt = Pin(ixy)%FxrIdxSt
        nLocalFxr = CellInfo(icel)%nFxr
        DO k = 1, nLocalFxr
          ifxr = FxrIdxSt + k - 1
          Fxrsave(ifxr, iz)%wt = FmInfo%Fxr(ifxr, iz)%wt
          Fxrsave(ifxr, iz)%h2ofrac = FmInfo%Fxr(ifxr, iz)%h2ofrac

          Fxrsave(ifxr, iz)%iso0 = FmInfo%Fxr(ifxr, iz)%iso0
          Fxrsave(ifxr, iz)%iso1 = FmInfo%Fxr(ifxr, iz)%iso1
          Fxrsave(ifxr, iz)%imix = FmInfo%Fxr(ifxr, iz)%imix
          Fxrsave(ifxr, iz)%ndim = FmInfo%Fxr(ifxr, iz)%ndim
          Fxrsave(ifxr, iz)%niso = FmInfo%Fxr(ifxr, iz)%niso
          
          Fxrsave(ifxr, iz)%lh2o = FmInfo%Fxr(ifxr, iz)%lh2o
          Fxrsave(ifxr, iz)%lres = FmInfo%Fxr(ifxr, iz)%lres
          Fxrsave(ifxr, iz)%lCusping = FmInfo%Fxr(ifxr, iz)%lCusping
          Fxrsave(ifxr, iz)%lMixtureMixing = FmInfo%Fxr(ifxr, iz)%lMixtureMixing

          ndim = FmInfo%Fxr(ifxr, iz)%ndim
          IF(allocated(Fxrsave(ifxr,iz)%pnum)) THEN
            DEALLOCATE(Fxrsave(ifxr, iz)%pnum, Fxrsave(ifxr, iz)%IdIso)
          END IF
          ALLOCATE(Fxrsave(ifxr, iz)%pnum(ndim), Fxrsave(ifxr, iz)%IdIso(ndim))
          Fxrsave(ifxr, iz)%pnum(1:ndim) = FmInfo%Fxr(ifxr, iz)%pnum(1:ndim)
          Fxrsave(ifxr, iz)%IdIso(1:ndim) = FmInfo%Fxr(ifxr, iz)%IdIso(1:ndim)
        END DO 
      END DO 
    END DO 
  END DO 
END DO

NULLIFY(XsChange, Fxr, AsyType, Asy, Pin, CellInfo, CoreMap)
END SUBROUTINE

SUBROUTINE RecoverXsPerturbation(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, PE, T, Tsave)
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
REAL :: T, Tsave

TYPE(XsChange_Type), POINTER :: XsChange(:)
REAL :: tbeg, tend
INTEGER :: nChange, FxrIdxSt, nLocalFxr, ndim, ndim0
INTEGER :: izbeg, izend, iasy, ixy, ixy0, ifxr, iz, icel, iasytype
INTEGER :: i, j, k

IF(.NOT. nTracerCntl%lXslib) STOP
XsChange => TranCntl%XsChange
nChange = TranCntl%nChange
Fxr => FmInfo%Fxr
AsyType => Core%AsyInfo
CoreMap => Core%CoreMap
Asy => Core%Asy
Pin => Core%Pin
CellInfo => Core%CellInfo
myzb = PE%myzb
myze = PE%myze
nxy = Core%nxy

TranCntl%lXsPerturb = .FALSE.
DO i = 1, nChange
  tbeg = XsChange(i)%tbeg
  tend = XsChange(i)%tend
  IF(T .LT. tbeg .OR. Tsave .GT. tend) CYCLE
  TranCntl%lXsPerturb = .TRUE.
  izbeg = XsChange(i)%izbeg
  izend = XsChange(i)%izend

  XsChange(i)%lUse = XsChangeSave(i)%lUse
  XsChange(i)%lstart = XsChangeSave(i)%lstart
  XsChange(i)%lcomplete = XsChangeSave(i)%lcomplete
  XsChange(i)%wt = XsChangeSave(i)%wt
  
  DO j = 1, XsChange(i)%nasy
    iasy = XsChange(i)%AsyList(j)
    iasytype = Core%CoreMap(iasy)
    DO iz = PE%myzb, PE%myze
      IF(iz .GT. izend .OR. iz .LT. izbeg) CYCLE
      DO ixy0 = 1, AsyType(iasytype)%nxy
        ixy = Asy(iasy)%GlobalPinIdx(ixy0)
        icel = Pin(ixy)%Cell(iz)
        FxrIdxSt = Pin(ixy)%FxrIdxSt
        nLocalFxr = CellInfo(icel)%nFxr
        DO k = 1, nLocalFxr
          ifxr = FxrIdxSt + k - 1
          FmInfo%Fxr(ifxr, iz)%wt = Fxrsave(ifxr, iz)%wt
          FmInfo%Fxr(ifxr, iz)%h2ofrac = Fxrsave(ifxr, iz)%h2ofrac

          ndim0 = Fxr(ifxr, iz)%ndim
          ndim = Fxrsave(ifxr, iz)%ndim
          FmInfo%Fxr(ifxr, iz)%iso0 = Fxrsave(ifxr, iz)%iso0
          FmInfo%Fxr(ifxr, iz)%iso1 = Fxrsave(ifxr, iz)%iso1
          FmInfo%Fxr(ifxr, iz)%imix = Fxrsave(ifxr, iz)%imix
          FmInfo%Fxr(ifxr, iz)%ndim = Fxrsave(ifxr, iz)%ndim
          FmInfo%Fxr(ifxr, iz)%niso = Fxrsave(ifxr, iz)%niso
          
          FmInfo%Fxr(ifxr, iz)%lh2o = Fxrsave(ifxr, iz)%lh2o
          FmInfo%Fxr(ifxr, iz)%lres = Fxrsave(ifxr, iz)%lres
          FmInfo%Fxr(ifxr, iz)%lCusping = Fxrsave(ifxr, iz)%lCusping
          FmInfo%Fxr(ifxr, iz)%lMixtureMixing = Fxrsave(ifxr, iz)%lMixtureMixing

          IF(ndim0 .NE. ndim) THEN
            DEALLOCATE(FmInfo%Fxr(ifxr, iz)%pnum, FmInfo%Fxr(ifxr, iz)%IdIso)
            ALLOCATE(FmInfo%Fxr(ifxr, iz)%pnum(ndim), FmInfo%Fxr(ifxr, iz)%IdIso(ndim))
          END IF
          FmInfo%Fxr(ifxr, iz)%pnum(1:ndim) = Fxrsave(ifxr, iz)%pnum(1:ndim)
          FmInfo%Fxr(ifxr, iz)%IdIso(1:ndim) = Fxrsave(ifxr, iz)%IdIso(1:ndim)
        END DO 
      END DO 
    END DO 
  END DO 
END DO
NULLIFY(XsChange, Fxr, AsyType, Asy, Pin, CellInfo, CoreMap)
END SUBROUTINE

SUBROUTINE ChangeMixture(myFxr, iso1, iso2, wt, wt0, lComplete)
USE PARAM
USE Material_Mod,     ONLY : Mixture
USE BasicOperation,   ONLY : CP_VA
USE IOUTIL,           ONLY : TERMINATE
USE ALLOCS
IMPLICIT NONE
TYPE(FxrInfo_Type) :: myFxr
INTEGER :: iso1, iso2
REAL :: wt, wt0
LOGICAL :: lComplete

INTEGER :: Idiso_new(1024)
REAL :: pnum_new(1024)
INTEGER :: niso_new

INTEGER :: i, j
INTEGER :: id
REAL :: wtbar

LOGICAL :: lh2o_1, lh2o_2

wtbar = 1._8 - wt
niso_new = 0
IF(.NOT. lComplete .OR. (wt .NE. 1._8)) THEN
  niso_new = myFxr%niso
  CALL CP_VA(IdIso_new(1:niso_new), myFxr%IdIso(1:niso_new), niso_new)
  !CALL CP_VA(pnum_new(1:niso_new), myFxr%pnum(1:niso_new), niso_new)
  DO i = 1, niso_new
    pnum_new(i) = wtbar * myFxr%pnum(i)
  ENDDO
ENDIF
DO i = 1, Mixture(iso2)%niso
  id = Mixture(iso2)%idiso(i)
  DO j = 1, niso_new
    IF(id .EQ. IdIso_new(j)) EXIT
  ENDDO
  IF(j .EQ. niso_new + 1) THEN
    niso_new = niso_new + 1;     IdIso_new(j) = id
    pnum_new(j) = 0
  ENDIF
  pnum_new(j) = pnum_new(j) + wt * Mixture(iso2)%pnum(i)
ENDDO

IF(niso_new .GT. myFxr%ndim) THEN
  DEALLOCATE(myFxr%pnum, myFxr%Idiso)
  CALL Dmalloc(myFxr%pnum, niso_new + 4); CALL Dmalloc(myFxr%IdIso, niso_new + 4)
  myFxr%ndim = niso_new + 4
  IF(myFxr%lres .OR. myFxr%lCrRes) THEN 
    DEALLOCATE(myFxr%idiso_pastpsm, myFxr%idx_Res, myFxr%pnum_Res)
    CALL dmalloc(myFxr%idiso_pastpsm, myFxr%ndim)
    CALL dmalloc(myFxr%idx_Res, myFxr%ndim)
    CALL dmalloc(myFxr%pnum_Res, myFxr%ndim)
  END IF
ENDIF

myFxr%niso = niso_new
DO i = 1, niso_new
  myFxr%idiso(i) = idiso_new(i); myFxr%pnum(i) = pnum_new(i)
ENDDO

!IF(myFxr%lres .NEQV. mixture(iso2)%lres) THEN
!  !Change non-res isotope to res-isotope
!  CALL TERMINATE('XS Change between non-res. mixture and res. mixture is not allowed.')
!ENDIF

lh2o_1 = Mixture(iso1)%lh2o .AND. Mixture(iso2)%lh2o; lh2o_2 = Mixture(iso1)%lh2o .OR. Mixture(iso2)%lh2o

myFxr%lMixtureMixing = .TRUE.

IF(lh2o_1 .AND. lh2o_2) THEN
  myFxr%h2ofrac = Mixture(iso1)%h2ofrac0 * (1-wt0) + Mixture(iso2)%h2ofrac0 * wt0
ELSEIF(.NOT. lh2o_1 .AND. lh2o_2) THEN
  IF(Mixture(iso1)%lh2o) myFxr%h2ofrac = Mixture(iso1)%h2ofrac0 * (1-wt0)
  IF(Mixture(iso2)%lh2o) myFxr%h2ofrac = Mixture(iso2)%h2ofrac0 * wt0
  myFxr%lh2o = .TRUE.
ENDIF

myFxr%lRes = Mixture(iso1)%lRes .OR. Mixture(iso1)%lRes
myFxr%lAIC = Mixture(iso1)%lAIC .OR. Mixture(iso1)%lAIC

IF(lComplete) THEN
  myFxr%lh2o = Mixture(iso2)%lh2o
  myFxr%imix = iso2
  myFxr%lMixtureMixing = .FALSE.
  myFxr%lRes = Mixture(iso2)%lRes
  myFxr%lAIC = Mixture(iso2)%lAIC
ENDIF
END SUBROUTINE

SUBROUTINE AsyXsChangeLib(iso1, iso2, iz1, iz2, iasy, wt, wt0, lComplete, lCusping, lCuspingDirection)
INTEGER :: iso1, iso2, iz1, iz2, iasy
REAL :: wt, wt0
INTEGER :: j
INTEGER :: ixy, ixy0, ifxr, iz, icel, iasytype
INTEGER :: FxrIdxSt, nLocalFxr
LOGICAL :: lComplete
LOGICAL :: lCusping, lCuspingDirection

iasytype = CoreMap(iasy)
DO iz = myzb, myze
  IF(iz .GT. iz2 .OR. iz .LT. iz1) CYCLE
  DO ixy0 = 1, AsyType(iasytype)%nxy
    ixy = Asy(iasy)%GlobalPinIdx(ixy0); icel = Pin(ixy)%Cell(iz)
    FxrIdxSt = Pin(ixy)%FxrIdxst; nLocalFxr = CellInfo(icel)%nFxr
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j - 1
      IF(FXR(ifxr, iz)%imix .NE. iso1) CYCLE
      CALL ChangeMixture(Fxr(ifxr, iz), iso1, iso2, wt, wt0, lComplete)
      IF(lCusping) THEN
        Fxr(ifxr, iz)%lCusping = .TRUE.
        IF(lCuspingDirection) THEN
          FXR(ifxr, iz)%wt = wt0
          FXR(ifxr, iz)%iso0 = iso1
          FXR(ifxr, iz)%iso1 = iso2    
        ELSE
          FXR(ifxr, iz)%wt = 1._8 - wt0
          FXR(ifxr, iz)%iso0 = iso2
          FXR(ifxr, iz)%iso1 = iso1
        ENDIF
      ENDIF
      IF(lComplete) THEN
        FXR(ifxr, iz)%imix = iso2;
      ENDIF
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE AsyXsChangeBen(iso1, iso2, isonew, iz1, iz2, iasy, lComplete)
INTEGER :: iso1, iso2, isonew, iz1, iz2, iasy
REAL :: wt
INTEGER :: j
INTEGER :: ixy, ixy0, ifxr, iz, icel, iasytype
INTEGER :: FxrIdxSt, nLocalFxr
LOGICAL :: lComplete

iasytype = CoreMap(iasy)
DO iz = myzb, myze
  IF(iz .GT. iz2 .OR. iz .LT. iz1) CYCLE
  DO ixy0 = 1, AsyType(iasytype)%nxy
    ixy = Asy(iasy)%GlobalPinIdx(ixy0); icel = Pin(ixy)%Cell(iz)
    FxrIdxSt = Pin(ixy)%FxrIdxst; nLocalFxr = CellInfo(icel)%nFxr
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j - 1
      IF(FXR(ifxr, iz)%imix0 .NE. iso1) CYCLE
      FXR(ifxr, iz)%imix = isonew
      IF(lComplete) THEN
        FXR(ifxr, iz)%imix = iso2; FXR(ifxr, iz)%imix0 = iso2
      ENDIF
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE

END MODULE
