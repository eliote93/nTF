MODULE XsPerturb_mod
USE PARAM
USE TYPEDEF,            ONLY : CoreInfo_Type,           FmInfo_Type,          TranInfo_Type,      &
                               TranCntl_Type,           PE_TYPE,                                  &
                               FxrInfo_Type,            AsyInfo_Type,         XsChange_Type,      &
                               Asy_Type,                Pin_Type,             Cell_Type
USE CNTL,               ONLY : nTracerCntl_Type
USE BenchXs,            ONLY : XsBenChange_new
IMPLICIT NONE

TYPE(FxrInfo_Type), POINTER, PRIVATE :: Fxr(:, :)
TYPE(AsyInfo_Type), PRIVATE, POINTER :: AsyType(:)
TYPE(Asy_Type), PRIVATE, POINTER :: Asy(:)
TYPE(Pin_Type), PRIVATE, POINTER :: Pin(:)
TYPE(Cell_Type), PRIVATE, POINTER :: CellInfo(:)
INTEGER, PRIVATE, POINTER :: CoreMap(:)


INTEGER, PRIVATE :: myzb, myze, nxy

CONTAINS
SUBROUTINE XsPerturbation(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, PE)

IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type)  :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

TYPE(XsChange_Type), POINTER :: XsChange(:)

REAL :: T, Tprev, Tprev0, Tbeg, Tend
REAL :: wt, wt0
INTEGER :: nChange, NowStep
INTEGER :: i, j, iasy, izbeg, izend

LOGICAL :: lXsLib

INTEGER, PARAMETER :: xsp_crit = 0.05_8

XsChange => TranCntl%XsChange
lXsLib = nTracerCntl%lXslib
Fxr => FmInfo%Fxr; AsyType => Core%AsyInfo; CoreMap => Core%CoreMap
Asy => Core%Asy; Pin => Core%Pin; CellInfo => Core%CellInfo
myzb = PE%myzb; myze = PE%myze; nxy = Core%nxy

TranCntl%lXsPerturb = .FALSE.

nChange = TranCntl%nChange; NowStep = TranCntl%NowStep
T = TranCntl%T(NowStep); Tprev = 0


IF(NowStep .GT. 1) Tprev = TranCntl%T(NowStep - 1)

DO i = 1, nChange
  Tbeg = XsChange(i)%Tbeg; Tend = XsChange(i)%Tend
  IF(Tend .GE. T) CYCLE
  IF(XsChange(i)%lComplete) CYCLE
  XsChange(i)%lComplete = .TRUE.; XsChange(i)%lStart = .TRUE.
  wt = 1._8; wt0 = 1._8
  IF(.NOT. lXsLib) CALL XsBenChange_New(XsChange(i)%Iso0, XsChange(i)%Iso1, XsChange(i)%Isonew, wt)
  izbeg = XsChange(i)%izbeg; izend = XsChange(i)%izend
  DO j = 1, XsChange(i)%nasy
    iasy = XsChange(i)%AsyList(j)
    IF(lXsLib) THEN
      CALL AsyXsChangeLib(XsChange(i)%Iso0, XsChange(i)%Iso1, izbeg, izend, iasy, wt, wt0, .TRUE.)
    ELSE
      CALL AsyXsChangeBen(XsChange(i)%Iso0, XsChange(i)%Iso1, XsChange(i)%Isonew, izbeg, izend, iasy, .TRUE.)
    ENDIF
  ENDDO
  IF(abs(wt0 - XsChange(i)%wt) > xsp_crit) TranCntl%lXsPerturb = .TRUE.
  XsChange(i)%wt = wt0
ENDDO

DO i = 1, nChange
  Tbeg = XsChange(i)%TBeg; Tend =XsChange(i)%Tend
  IF(T .LT. Tbeg .OR. T .GT. Tend) CYCLE
  XsChange(i)%lStart = .TRUE.
  Tprev0 = Max(Tprev, Tbeg)
  wt = (T - Tprev0) / (Tend - Tprev0)
  wt0 = (T - Tbeg) / (Tend - Tbeg)
  IF(.NOT. lXsLib) CALL XsBenChange_New(XsChange(i)%Iso0, XsChange(i)%Iso1, XsChange(i)%Isonew, wt)
  izbeg = XsChange(i)%izbeg; izend = XsChange(i)%izend
  DO j = 1, XsChange(i)%nasy
    iasy = XsChange(i)%AsyList(j)
    IF(lXsLib) THEN
      CALL AsyXsChangeLib(XsChange(i)%Iso0, XsChange(i)%Iso1, izbeg, izend, iasy, wt, wt0, .FALSE.)
    ELSE
      CALL AsyXsChangeBen(XsChange(i)%Iso0, XsChange(i)%Iso1, XsChange(i)%Isonew, izbeg, izend, iasy, .FALSE.)
    ENDIF
  ENDDO
  IF(abs(wt0 - XsChange(i)%wt) > xsp_crit) TranCntl%lXsPerturb = .TRUE.
  XsChange(i)%wt = wt0
ENDDO

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
ENDIF

myFxr%niso = niso_new
DO i = 1, niso_new
  myFxr%idiso(i) = idiso_new(i); myFxr%pnum(i) = pnum_new(i)
ENDDO

IF(myFxr%lres .NEQV. mixture(iso2)%lres) THEN
  !Change non-res isotope to res-isotope
  CALL TERMINATE('XS Change between non-res. mixture and res. mixture is not allowed.')
ENDIF

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

IF(lComplete) THEN
  myFxr%lh2o = Mixture(iso2)%lh2o
  myFxr%imix = iso2
  myFxr%lMixtureMixing = .FALSE.
  myFxr%lRes = Mixture(iso2)%lRes
ENDIF
END SUBROUTINE

SUBROUTINE AsyXsChangeLib(iso1, iso2, iz1, iz2, iasy, wt, wt0, lComplete)
INTEGER :: iso1, iso2, iz1, iz2, iasy
REAL :: wt, wt0
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
      IF(FXR(ifxr, iz)%imix .NE. iso1) CYCLE
      CALL ChangeMixture(Fxr(ifxr, iz), iso1, iso2, wt, wt0, lComplete)
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
