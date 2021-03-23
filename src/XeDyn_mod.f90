#include <defines.h>
MODULE XeDyn_Mod
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,       FmInfo_Type,       GroupInfo_Type,    &
                             THInfo_Type,         PE_Type,                              &
                             XeDynInfo_Type,      XeDynState_Type
USE DeplType,         ONLY : DeplCntl_Type
USE CNTL,             ONLY : nTracerCntl_Type
IMPLICIT NONE

!INTEGER, PARAMETER :: FirXeND = 0
INTEGER, PARAMETER :: UpdtXe = 1
INTEGER, PARAMETER :: InitXe = 2
INTEGER, PARAMETER :: FinalXe = 3

LOGICAL :: lRelaxation = .FALSE.

TYPE(XeDynInfo_Type), SAVE :: XeDynInfo

INTERFACE

SUBROUTINE SetXeDynCoreState(istep, XeDynInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,         ONLY : XeDynInfo_Type,         PE_Type
USE CNTL,            ONLY : nTracerCntl_Type

IMPLICIT NONE
TYPE(XeDynInfo_Type) :: XeDynInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
INTEGER :: istep

END SUBROUTINE

SUBROUTINE InitXeDynInfo(Fxr, nfxr, myzb, myze)
USE PARAM
USE TYPEDEF,         ONLY : FxrInfo_Type,     XeDynFxr_Type
!USE CNTL,            ONLY : nTracerCntl_Type

TYPE(FxrInfo_Type) :: Fxr(nfxr, myzb:myze)
INTEGER :: nfxr, myzb, myze

END SUBROUTINE

SUBROUTINE SetDeplXeDynMode(Core, FmInfo, DeplCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,   FmInfo_Type,  PE_Type
USE DeplType,         ONLY : DeplCntl_Type
USE CNTL,             ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(DeplCntl_Type) :: DeplCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
END SUBROUTINE

SUBROUTINE UpdtXeDyn(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, DeplCntl, PE, Mode)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,       FmInfo_Type,       GroupInfo_Type,    &
                             THInfo_Type,         PE_Type
USE DeplType,         ONLY : DeplCntl_Type
USE CNTL,             ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(DeplCntl_Type) :: DeplCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
INTEGER, OPTIONAL :: Mode

END SUBROUTINE

SUBROUTINE TrXeUpdate(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, DeplCntl, PE, Mode)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,       FmInfo_Type,       GroupInfo_Type,    &
                             THInfo_Type,         PE_Type
USE DeplType,         ONLY : DeplCntl_Type
USE CNTL,             ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(DeplCntl_Type) :: DeplCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
INTEGER :: Mode
END SUBROUTINE

SUBROUTINE EqXeUpdate(Core, FmInfo, ThInfo, GroupInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,       FmInfo_Type,       GroupInfo_Type,    &
                             THInfo_Type,         PE_Type
USE CNTL,             ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE

END SUBROUTINE

SUBROUTINE WriteXeND(Core, FmInfo, PE)
USE PARAM
USE TYPEDEF,    ONLY : CoreInfo_TYpe,    FmInfo_Type,   PE_TYPE,           &
                       FxrInfo_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(PE_Type) :: PE

END SUBROUTINE

SUBROUTINE ReadXeND(Core, FmInfo, PE)
USE PARAM
USE TYPEDEF,    ONLY : CoreInfo_TYpe,    FmInfo_Type,   PE_TYPE,           &
                       FxrInfo_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(PE_Type) :: PE

END SUBROUTINE

END INTERFACE

CONTAINS

SUBROUTINE XeDynRelaxContrl(activity)
IMPLICIT NONE
LOGICAL :: activity
lRelaxation = activity
END SUBROUTINE

SUBROUTINE GetInp_TimeStep(Oneline)
USE IOUTIL,            ONLY : FndChara,        nfields,         Terminate
IMPLICIT NONE
CHARACTER(256) :: Oneline
INTEGER :: nspt, ipos(100)
REAL :: tend, delt, Tbeg, T
INTEGER :: i, j, k
INTEGER :: nsub, nField
CHARACTER(80) :: astring
CHARACTER(512) :: aline

CALL FndChara(Oneline, ipos, nspt, SLASH)

READ(oneline, *) astring, tend, DelT
nsub = nint(tend/DelT)

i = 0; XeDynInfo%Tsec(0) = 0
Tbeg = 0
DO j = 1, nsub
  i = i + 1
  T = T + delt
  XeDynInfo%Tsec(i) = T
ENDDO
Tbeg = T

DO j = 1, nspt
  aline = oneline(ipos(j)+1:256)
  nField = nFields(aline)
  IF(nField .LT. 2) CYCLE
  READ(oneline(ipos(j)+1:256), *) tend, DelT

  nsub = NINT((Tend - Tbeg)/DelT)
  DO k = 1, nsub
    i = i + 1
    T = T + DelT
    XeDynInfo%Tsec(i) = T
  ENDDO
  Tbeg = T
ENDDO

XeDynInfo%nTimeStep = i

END SUBROUTINE

SUBROUTINE GetInp_TimeUnit(oneline)
USE IOUTIL,              ONLY : TOUPPER
IMPLICIT NONE
CHARACTER(256) :: Oneline
CHARACTER(10) :: astring, Unit

READ(oneline, *) astring, Unit

CALL TOUPPER(Unit)
XeDynInfo%Unit = UNIT

END SUBROUTINE

SUBROUTINE GetInp_CoreState(oneline)
USE ioutil,              ONLY : FndChara,      nFields,        Terminate,            &
                                IFNumeric
IMPLICIT NONE
CHARACTER(256) :: Oneline
CHARACTER(256) :: ErrMesg
CHARACTER(512) :: aline

INTEGER :: id
INTEGER :: nField, nspt, ipos(100)
REAL :: t, PowLv, Val
LOGICAL :: LSearch, lCntl

CHARACTER(50) :: astring

READ(oneline, *, Err = 1000) Astring, id

CALL FndChara(Oneline, ipos, nspt, SLASH)

IF(nspt .NE. 4) THEN
  GOTO 2000
ENDIF

READ(oneline, *, Err = 2000) astring, id, T
READ(oneline(ipos(1)+1:ipos(2)-1), *, Err = 2000) PowLv
READ(oneline(ipos(2)+1:ipos(3)-1), *, Err = 2000) LSearch, astring

XeDynInfo%CoreState(id)%T = T
XeDynInfo%CoreState(id)%PowLv  = PowLv / 100

IF(IfNumeric(astring)) THEN
  READ(astring, *) val
ELSEIF(.NOT. Lsearch) THEN
  XeDynInfo%CoreState(id)%lPrevBoron = .TRUE.
ELSE
  GOTO 2000
ENDIF

IF(.NOT. XeDynInfo%CoreState(id)%lPrevBoron) THEN
  XeDynInfo%CoreState(id)%lBoronSearch = LSearch
  IF(Lsearch) THEN
    XeDynInfo%CoreState(id)%Target_keff = val
  ELSE
    XeDynInfo%CoreState(id)%BoronPPM = val
  ENDIF
ENDIF

lCntl = .TRUE.
aline = oneline(ipos(3)+1:ipos(4)-1)
nField = nFields(aline)
IF(nfield .EQ. 0) lCntl = .FALSE.
XeDynInfo%CoreState(id)%lCntlRod = lCntl
IF(lCntl) THEN
  READ(oneline(ipos(3)+1:ipos(4)-1), *) XeDynInfo%CoreState(id)%CntlRodPos
ENDIF

XeDynInfO%CoreState%Inpline = oneline
XeDynInfo%nState = max(XeDynInfo%nState, id)

RETURN
1000 CALL TERMINATE('XEDYN Input Error - CORE_STATE')
2000 WRITE(ErrMesg, '(A, I5)') 'XEDYN Input Error - CORE_STATE :', id
CALL TERMINATE(ErrMesg)
END SUBROUTINE

SUBROUTINE ProcessXeDynInp()
USE ioutil,       ONLY : TERMINATE
IMPLICIT NONE

INTEGER :: i, j
REAL :: d2s, H2s, M2s
REAL :: s2d, s2h, s2m
REAL :: f, T
!Time Step Process
d2s = 24._8 * 3600._8; H2s = 3600._8; M2s = 60
s2d = 1._8 / d2s; s2h = 1._8 / h2s; s2m = 1._8 / m2s
f = 1
SELECT CASE(XeDynInfo%Unit)
  CASE('D')
    f = d2s
  CASE('H')
    f = h2s
  CASE('M')
    f = M2s
END SELECT


DO i = 1, XeDynInfo%nTimeStep
  XeDynInfo%Tsec(i) = XeDynInfo%Tsec(i) * f
  XeDynInfo%Tmin(i) = XeDynInfo%Tsec(i) * s2m
  XeDynInfo%Th(i) = XeDynInfo%Tsec(i) * s2h
  XeDynInfo%Tday(i) = XeDynInfo%Tsec(i) * s2d
ENDDO

DO i = 1, XeDynInfo%nState
  XeDynInfo%CoreState(i)%T = XeDynInfo%CoreState(i)%T * f
ENDDO

XeDynInfo%Tsec(0) = 1.0E-6_8
!Core State Mapping
DO i = 0, XeDynInfo%nTimeStep
  XeDynInfo%lPerturb(i) = .FALSE.
  XeDynInfo%StateMapping(i) = -1
ENDDO

DO i = 1, XeDynInfo%nstate
  T = XeDynInfo%CoreState(i)%T
  DO j = 0, XeDynInfo%nTimeStep
    IF(abs(T-XeDynInfo%Tsec(j)) .LT. epsm4) THEN
      XeDynInfo%StateMapping(j) = i
      XeDynInfo%lPerturb(j) = i
    ENDIF
  ENDDO
ENDDO

END SUBROUTINE

END MODULE
