MODULE MOC_COMMON

USE TYPEDEF,        ONLY : XsMac_Type
USE TYPEDEF_COMMON
USE OMP_LIB
IMPLICIT NONE

TYPE(XsMac_Type), POINTER, PRIVATE :: XsMac(:)
TYPE(XsMac_Type), POINTER, PRIVATE :: XsMac0(:), XsMac1(:)
TYPE(CoreXsMac_Type), ALLOCATABLE :: CoreXsMac(:)
INTEGER, ALLOCATABLE :: InScatRange(:, :), InScatIdx(:, :)

CONTAINS

SUBROUTINE AllocCoreMacXs(CoreInfo)
USE TYPEDEF,        ONLY : CoreInfo_Type
USE PE_MOD,         ONLY : PE
USE CNTL,           ONLY : nTracerCntl
USE CORE_MOD,       ONLY : GroupInfo
USE XsUtil_mod,     ONLY : AllocXsMac
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
INTEGER :: ng, ngs, nFxr, nThread
INTEGER :: iz, ig, igs, tid
INTEGER :: myzb, myze

ng = GroupInfo%ng
nFxr = CoreInfo%nCoreFxr
nThread = PE%nThread
myzb = PE%myzb
myze = PE%myze

ALLOCATE(InScatRange(2, ng), InScatIdx(ng, ng))

InScatRange = GroupInfo%InScatRange
InScatIdx = 0

ngs = 0;
!DO ig = 2, ng
DO ig = 1, ng
  DO igs = InScatRange(1, ig), InScatRange(2, ig)
    ngs = ngs + 1
    InScatIdx(igs, ig) = ngs
  ENDDO
ENDDO

CALL omp_set_num_threads(nThread)

ALLOCATE(XsMac(nThread))
ALLOCATE(XsMac0(nThread))
ALLOCATE(XsMac1(nThread))

!$OMP PARALLEL PRIVATE(tid)
tid = omp_get_thread_num() + 1
  XsMac(tid)%ng = ng
XsMac0(tid)%ng = ng
XsMac1(tid)%ng = ng
  CALL AllocXsMac(XsMac(tid))
CALL AllocXsMac(XsMac0(tid))
CALL AllocXsMac(XsMac1(tid))
!$OMP END PARALLEL

ALLOCATE(CoreXsMac(myzb : myze))

DO iz = myzb, myze
  ALLOCATE(CoreXsMac(iz)%XSt(ng, nFxr))
  ALLOCATE(CoreXsMac(iz)%XStr(ng, nFxr))
  ALLOCATE(CoreXsMac(iz)%XSa(ng, nFxr))
  ALLOCATE(CoreXsMac(iz)%XSnf(ng, nFxr))
  ALLOCATE(CoreXsMac(iz)%XSkf(ng, nFxr))
  ALLOCATE(CoreXsMac(iz)%XSsm(ngs, nFxr))
  IF (nTracerCntl%ScatOd .GE. 1) ALLOCATE(CoreXsMac(iz)%XSsmP1(ngs, nFxr))
  IF (nTracerCntl%ScatOd .GE. 2) ALLOCATE(CoreXsMac(iz)%XSsmP2(ngs, nFxr))
  IF (nTracerCntl%ScatOd .EQ. 3) ALLOCATE(CoreXsMac(iz)%XSsmP3(ngs, nFxr))
ENDDO

END SUBROUTINE

SUBROUTINE SetCoreMacXs(CoreInfo, FmInfo)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,       FmInfo_Type,        FxrInfo_Type
USE PE_MOD,         ONLY : PE
USE CNTL,           ONLY : nTracerCntl
USE CORE_MOD,       ONLY : GroupInfo
USE IOUTIL,         ONLY : message
USE FILES,          ONLY : io8
USE BenchXs,        ONLY : xsbaseBen,           xsbaseDynBen,       xsbaseben_NEACRP
USE MacXsLib_Mod,   ONLY : MacXsBase,           MacXsScatMatrix,    MacP1XsScatMatrix,                              &
                           MacP2XsScatMatrix,   MacP3XsScatMatrix
USE MPIComm_Mod,    ONLY : MPI_SYNC
USE TRAN_MOD,       ONLY : TranInfo
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(FmInfo_Type) :: FmInfo

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :), myFxr
REAL, POINTER :: XSt(:, :), XStr(:, :), XSa(:, :), XSnf(:, :), XSkf(:, :), XSsm(:, :)
INTEGER :: ng, nFxr
INTEGER :: ig, igf, igs, itype, ifxr, iz, irgb, irge, tid
INTEGER :: myzb, myze
LOGICAL :: lxslib, lscat1, lTrCorrection

Fxr => FmInfo%Fxr
ng = GroupInfo%ng
nFxr = CoreInfo%nCoreFxr
myzb = PE%myzb
myze = PE%myze
lxslib = nTracerCntl%lxslib
lscat1 = nTracerCntl%lscat1
lTrCorrection = nTracerCntl%lTrCorrection

IF (lxslib) THEN
  irgb = GroupInfo%nofg + 1
  irge = GroupInfo%nofg + GroupInfo%norg
ENDIF

WRITE(mesg, '(a)') 'Calculating Macro Cross Sections...'
IF (PE%MASTER) CALL message(io8, .TRUE., .TRUE., mesg)

!$OMP PARALLEL PRIVATE(tid, igs, itype, myFxr)
tid = omp_get_thread_num() + 1

DO iz = myzb, myze
  !$OMP DO SCHEDULE(GUIDED)
  DO ifxr = 1, nFxr
    myFxr => Fxr(ifxr, iz)
    IF (lxslib) THEN
      IF (myFxr%niso .EQ. 0) CYCLE
      CALL MacXsBase(XsMac(tid), myFxr, 1, ng, ng, 1.0, .FALSE., .FALSE., .TRUE.)
      CALL MacXsScatMatrix(XsMac(tid), myFxr, 1, ng, ng, GroupInfo, .FALSE., .TRUE.)
      IF (nTracerCntl%ScatOd .GE. 1) CALL MacP1XsScatMatrix(XsMac(tid), myFxr, 1, ng, ng, GroupInfo)
      IF (nTracerCntl%ScatOd .GE. 2) CALL MacP2XsScatMatrix(XsMac(tid), myFxr, 1, ng, ng, GroupInfo)
      IF (nTracerCntl%ScatOd .EQ. 3) CALL MacP3XsScatMatrix(XsMac(tid), myFxr, 1, ng, ng, GroupInfo)
      IF (myFxr%lres) THEN
        XsMac(tid)%XsMacA(irgb : irge) = myFxr%FresoA(irgb : irge) * XsMac(tid)%XsMacA(irgb : irge)
        IF (CoreInfo%lFuelPlane(iz)) THEN
          XsMac(tid)%XsMacNf(irgb : irge) = myFxr%fresonf(irgb : irge) * XsMac(tid)%XsMacNf(irgb : irge)
          XsMac(tid)%XsMacKf(irgb : irge) = myFxr%fresokf(irgb : irge) * XsMac(tid)%XsMacKf(irgb : irge)
        ENDIF
      ENDIF
      XsMac(tid)%XsMacTr = XsMac(tid)%XsMacA + XsMac(tid)%XsMacStr
      XsMac(tid)%XsMacT = XsMac(tid)%XsMacA + XsMac(tid)%XsMacS
    ELSE
      itype = myFxr%imix
      IF(itype .GT. 0) THEN
        IF(nTracerCntl%libtyp .EQ. 11) THEN
          CALL xsbaseben_NEACRP(itype, myFxr%rho, myFxr%temp, myFxr%DopTemp, XsMac(tid))
        ELSE
          IF(nTracerCntl%lDynamicBen) THEN
            CALL xsbaseDynBen(itype, TranInfo%fuelTemp(myFxr%ipin, iz), 1, ng, 1, ng, lscat1, XsMac(tid))
          ELSE
            CALL xsbaseBen(itype, 1, ng, 1, ng, lscat1, XsMac(tid))
          ENDIF
        END IF
      END IF
    ENDIF
    CoreXsMac(iz)%XSt(:, ifxr) = XsMac(tid)%XsMacT
    CoreXsMac(iz)%XStr(:, ifxr) = XsMac(tid)%XsMacTr
    CoreXsMac(iz)%XSa(:, ifxr) = XsMac(tid)%XsMacA
    CoreXsMac(iz)%XSnf(:, ifxr) = XsMac(tid)%XsMacNf
    CoreXsMac(iz)%XSkf(:, ifxr) = XsMac(tid)%XsMacKf
    DO ig = 1, ng
      DO igf = InScatRange(1, ig), InScatRange(2, ig)
        igs = InScatIdx(igf, ig)
        CoreXsMac(iz)%XSsm(igs, ifxr) = XsMac(tid)%XsMacSm(igf, ig)
        IF (nTracerCntl%ScatOd .GE. 1) CoreXsMac(iz)%XSsmP1(igs, ifxr) = XsMac(tid)%XsMacP1Sm(igf, ig)
        IF (nTracerCntl%ScatOd .GE. 2) CoreXsMac(iz)%XSsmP2(igs, ifxr) = XsMac(tid)%XsMacP2Sm(igf, ig)
        IF (nTracerCntl%ScatOd .EQ. 3) CoreXsMac(iz)%XSsmP3(igs, ifxr) = XsMac(tid)%XsMacP3Sm(igf, ig)
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
ENDDO

!$OMP END PARALLEL

CALL MPI_SYNC(PE%MPI_CMFD_COMM)

END SUBROUTINE

SUBROUTINE SetCoreMacXs_Cusping(CoreInfo, FmInfo)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,       FmInfo_Type,        FxrInfo_Type,           &
                           Pin_Type,            Cell_Type
USE PE_MOD,         ONLY : PE
USE CNTL,           ONLY : nTracerCntl
USE CORE_MOD,       ONLY : GroupInfo
USE IOUTIL,         ONLY : message
USE FILES,          ONLY : io8
USE BenchXs,        ONLY : xsbaseBen,           MacXsBen,           xsbaseBen_Cusping,      &
                           xsbaseDynBen,        xsbaseDynBen_Cusping, DynMacXsBen,          &
                           xsbaseBen_NEACRP
USE MacXsLib_Mod,   ONLY : MacXsBase,           MacXsScatMatrix,    BaseMacStr,             &
                           MacXsBase_Cusping,   MacXsScatMatrix_Cusping, BaseMacStr_Cusping, &
                           MacP1XsScatMatrix,   MacP1XsScatMatrix_Cusping, &
                           MacP2XsScatMatrix,   MacP2XsScatMatrix_Cusping, &
                           MacP3XsScatMatrix,   MacP3XsScatMatrix_Cusping
USE MPIComm_Mod,    ONLY : MPI_SYNC
USE TRAN_MOD,       ONLY : TranInfo
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(FmInfo_Type) :: FmInfo

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :), myFxr
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
REAL, POINTER :: XSt(:, :), XStr(:, :), XSa(:, :), XSnf(:, :), XSkf(:, :), XSsm(:, :)
REAL, ALLOCATABLE :: phiz(:, :), philz(:, :), phiuz(:, :)
REAL :: vol, volsum
REAL :: wt, wtbar, wtg, wtgbar, uflux, lflux, hz, hzl, hzu
INTEGER :: ifsr, ifsrlocal, nFsrinFxr, FsrIdxSt, FxrIdxSt
INTEGER :: ng, nFxr, nxy
INTEGER :: i, j, ipin, icel
INTEGER :: ig, igf, igs, itype, ifxr, iz, irgb, irge, tid, jg
INTEGER :: myzb, myze
LOGICAL :: lxslib, lscat1, lTrCorrection
LOGICAL :: lCusping

Fxr => FmInfo%Fxr
Pin => CoreInfo%Pin
Cell => CoreInfo%CellInfo
ng = GroupInfo%ng
nFxr = CoreInfo%nCoreFxr
nxy = CoreInfo%nxy
myzb = PE%myzb
myze = PE%myze
lxslib = nTracerCntl%lxslib
lscat1 = nTracerCntl%lscat1
lTrCorrection = nTracerCntl%lTrCorrection

ALLOCATE(phiz(ng, PE%nThread), philz(ng, PE%nThread), phiuz(ng, PE%nThread))

IF (lxslib) THEN
  irgb = GroupInfo%nofg + 1
  irge = GroupInfo%nofg + GroupInfo%norg
ENDIF

WRITE(mesg, '(a)') 'Calculating Macro Cross Sections with AFW...'
IF (PE%MASTER) CALL message(io8, .TRUE., .TRUE., mesg)

!$OMP PARALLEL PRIVATE(tid, igs, itype, myFxr, FsrIdxSt, FxrIdxSt, icel, ifxr, nFsrinFxr, ifsr, ifsrlocal, vol, volsum, wt, wtbar, hz, hzl, hzu, uflux, lflux, wtg, wtgbar)
tid = omp_get_thread_num() + 1

DO iz = myzb, myze
  !$OMP DO SCHEDULE(GUIDED)
  DO ipin = 1, nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt
    FxrIdxSt = Pin(ipin)%FxrIdxSt
    icel = Pin(ipin)%Cell(iz)
    DO j = 1, Cell(icel)%nFxr
      ifxr = FxrIdxSt + j - 1
      nFsrInFxr = Cell(icel)%nFsrInFxr(j)
      myFxr => Fxr(ifxr, iz)
      IF (lxslib) THEN
        IF (myFxr%niso .EQ. 0) CYCLE
        IF(Fxr(ifxr,iz)%lcusping) THEN
          phiz(:, tid) = 0.
          philz(:, tid) = 0.
          phiuz(:, tid) = 0.
          volsum = 0.
          DO ig = 1, ng
            DO i = 1, nFsrInFxr
              ifsrlocal = Cell(icel)%MapFxr2FsrIdx(i, j)
              ifsr = FsrIdxSt + ifsrlocal - 1
              vol = Cell(icel)%Vol(ifsrlocal)
              IF(ig .EQ. 1) volsum = volsum + vol
              phiz(ig, tid) = phiz(ig, tid) + FmInfo%phis(ifsr, iz, ig) * vol
              IF(iz .EQ. myzb) THEN
                philz(ig, tid) = philz(ig, tid) + FmInfo%neighphis(ifsr, ig, BOTTOM) * vol
              ELSE
                philz(ig, tid) = philz(ig, tid) + FmInfo%phis(ifsr, iz-1, ig) * vol
              END IF
              IF(iz .EQ. myze) THEN
                phiuz(ig, tid) = phiuz(ig, tid) + FmInfo%neighphis(ifsr, ig, TOP) * vol
              ELSE
                phiuz(ig, tid) = phiuz(ig, tid) + FmInfo%phis(ifsr, iz+1, ig) * vol
              END IF
            END DO
            phiz(ig, tid) = phiz(ig, tid) / volsum
            philz(ig, tid) = philz(ig, tid) / volsum
            phiuz(ig, tid) = phiuz(ig, tid) / volsum
          END DO
          CALL MacXsBase_Cusping(XSMac0(tid), myFxr%iso0, myFxr, 1, ng, ng, 1., FALSE, FALSE)
          CALL MacXsBase_Cusping(XSMac1(tid), myFxr%iso1, myFxr, 1, ng, ng, 1., FALSE, FALSE)
          CALL MacXsBase(XsMac(tid), myFxr, 1, ng, ng, 1.0, .FALSE., .FALSE., .TRUE.)

          CALL MacXsScatMatrix_Cusping(XsMac0(tid), myFxr%iso0, myFxr, 1, ng, ng ,GroupInfo, .FALSE.)
          CALL MacXsScatMatrix_Cusping(XsMac1(tid), myFxr%iso1, myFxr, 1, ng, ng ,GroupInfo, .FALSE.)
          CALL MacXsScatMatrix(XsMac(tid), myFxr, 1, ng, ng, GroupInfo, .FALSE., .TRUE.)
          IF (lTrCorrection) THEN
            DO ig = 1, ng
              CALL BaseMacStr_Cusping(XsMac0(tid), myFxr%iso0, myFxr, ig, ng)
              CALL BaseMacStr_Cusping(XsMac1(tid), myFxr%iso1, myFxr, ig, ng)
              CALL BaseMacStr(XsMac(tid), myFxr, ig, ng, .TRUE.)
            ENDDO
          ENDIF
          IF (nTracerCntl%ScatOd .GE. 1) THEN
             CALL MacP1XsScatMatrix_Cusping(XsMac0(tid), myFxr%iso0, myFxr, 1, ng, ng, GroupInfo)
             CALL MacP1XsScatMatrix_Cusping(XsMac1(tid), myFxr%iso1, myFxr, 1, ng, ng, GroupInfo)
             CALL MacP1XsScatMatrix(XsMac(tid), myFxr, 1, ng, ng, GroupInfo)
          END IF
          IF (nTracerCntl%ScatOd .GE. 2) THEN
             CALL MacP2XsScatMatrix_Cusping(XsMac0(tid), myFxr%iso0, myFxr, 1, ng, ng, GroupInfo)
             CALL MacP2XsScatMatrix_Cusping(XsMac1(tid), myFxr%iso1, myFxr, 1, ng, ng, GroupInfo)
             CALL MacP2XsScatMatrix(XsMac(tid), myFxr, 1, ng, ng, GroupInfo)
          END IF
          IF (nTracerCntl%ScatOd .GE. 3) THEN
             CALL MacP3XsScatMatrix_Cusping(XsMac0(tid), myFxr%iso0, myFxr, 1, ng, ng, GroupInfo)
             CALL MacP3XsScatMatrix_Cusping(XsMac1(tid), myFxr%iso1, myFxr, 1, ng, ng, GroupInfo)
             CALL MacP3XsScatMatrix(XsMac(tid), myFxr, 1, ng, ng, GroupInfo)
          END IF

          XsMac0(tid)%XsMacTr = XsMac0(tid)%XsMacA + XsMac0(tid)%XsMacStr
          XsMac0(tid)%XsMacT = XsMac0(tid)%XsMacA + XsMac0(tid)%XsMacS

          XsMac1(tid)%XsMacTr = XsMac1(tid)%XsMacA + XsMac1(tid)%XsMacStr
          XsMac1(tid)%XsMacT = XsMac1(tid)%XsMacA + XsMac1(tid)%XsMacS

          wt = myFxr%wt; wtbar = 1. - wt
          hz = CoreInfo%hzfm(iz); hzl = CoreInfo%hzfm(iz-1); hzu = CoreInfo%hzfm(iz+1)
          DO ig = 1, ng
            uflux = (hzu * phiuz(ig, tid) + wt * hz * phiz(ig, tid)) / (hzu + wt * hz)
            lflux = (hzl * philz(ig, tid) + wtbar * hz * phiz(ig, tid)) / (hzl + wtbar * hz)
            wtg = wt * uflux / (wt * uflux + wtbar * lflux)
            wtgbar = 1. - wtg

            XsMac(tid)%XsMacT(ig) = XsMac0(tid)%XsMacT(ig) * wtgbar + XsMac1(tid)%XsMacT(ig) * wtg
            XsMac(tid)%XsMacTr(ig) = XsMac0(tid)%XsMacTr(ig) * wtgbar + XsMac1(tid)%XsMacTr(ig) * wtg

            XsMac(tid)%XsMacA(ig) = XsMac0(tid)%XsMacA(ig) * wtgbar + XsMac1(tid)%XsMacA(ig) * wtg
            XsMac(tid)%XsMacS(ig) = XsMac0(tid)%XsMacS(ig) * wtgbar + XsMac1(tid)%XsMacS(ig) * wtg

            XsMac(tid)%XsMacnf(ig) = XsMac0(tid)%XsMacnf(ig) * wtgbar + XsMac1(tid)%XsMacnf(ig) * wtg
            XsMac(tid)%XsMackf(ig) = XsMac0(tid)%XsMackf(ig) * wtgbar + XsMac1(tid)%XsMackf(ig) * wtg

            DO jg = 1, ng
              XsMac(tid)%XsMacSm(ig, jg) = XsMac0(tid)%XsMacSm(ig, jg) * wtgbar + XsMac1(tid)%XsMacSm(ig, jg) * wtg
            END DO
            IF (nTracerCntl%ScatOd .GE. 1) THEN
              DO jg = 1, ng
                XsMac(tid)%xsmacp1sm(ig, jg) = XsMac0(tid)%xsmacp1sm(ig, jg) * wtgbar + XsMac1(tid)%xsmacp1sm(ig, jg) * wtg
              END DO
            END IF
            IF (nTracerCntl%ScatOd .GE. 2) THEN
              DO jg = 1, ng
                XsMac(tid)%xsmacp2sm(ig, jg) = XsMac0(tid)%xsmacp2sm(ig, jg) * wtgbar + XsMac1(tid)%xsmacp2sm(ig, jg) * wtg
              END DO
            END IF
            IF (nTracerCntl%ScatOd .GE. 3) THEN
              DO jg = 1, ng
                XsMac(tid)%xsmacp3sm(ig, jg) = XsMac0(tid)%xsmacp3sm(ig, jg) * wtgbar + XsMac1(tid)%xsmacp3sm(ig, jg) * wtg
              END DO
            END IF
          END DO
          IF (myFxr%lres) THEN
            XsMac(tid)%XsMacA(irgb : irge) = myFxr%FresoA(irgb : irge) * XsMac(tid)%XsMacA(irgb : irge)
            IF (CoreInfo%lFuelPlane(iz)) THEN
              XsMac(tid)%XsMacNf(irgb : irge) = myFxr%fresonf(irgb : irge) * XsMac(tid)%XsMacNf(irgb : irge)
              XsMac(tid)%XsMacKf(irgb : irge) = myFxr%fresokf(irgb : irge) * XsMac(tid)%XsMacKf(irgb : irge)
            ENDIF
          ENDIF
        ELSE
          CALL MacXsBase(XsMac(tid), myFxr, 1, ng, ng, 1.0, .FALSE., .FALSE., .TRUE.)
          CALL MacXsScatMatrix(XsMac(tid), myFxr, 1, ng, ng, GroupInfo, .FALSE., .TRUE.)
          IF (lTrCorrection) THEN
            DO ig = 1, ng
              CALL BaseMacStr(XsMac(tid), myFxr, ig, ng, .TRUE.)
            ENDDO
          ENDIF
          IF (nTracerCntl%ScatOd .GE. 1) CALL MacP1XsScatMatrix(XsMac(tid), myFxr, 1, ng, ng, GroupInfo)
          IF (nTracerCntl%ScatOd .GE. 2) CALL MacP2XsScatMatrix(XsMac(tid), myFxr, 1, ng, ng, GroupInfo)
          IF (nTracerCntl%ScatOd .EQ. 3) CALL MacP3XsScatMatrix(XsMac(tid), myFxr, 1, ng, ng, GroupInfo)
          IF (myFxr%lres) THEN
            XsMac(tid)%XsMacA(irgb : irge) = myFxr%FresoA(irgb : irge) * XsMac(tid)%XsMacA(irgb : irge)
            IF (CoreInfo%lFuelPlane(iz)) THEN
              XsMac(tid)%XsMacNf(irgb : irge) = myFxr%fresonf(irgb : irge) * XsMac(tid)%XsMacNf(irgb : irge)
              XsMac(tid)%XsMacKf(irgb : irge) = myFxr%fresokf(irgb : irge) * XsMac(tid)%XsMacKf(irgb : irge)
            ENDIF
          ENDIF
          XsMac(tid)%XsMacTr = XsMac(tid)%XsMacA + XsMac(tid)%XsMacStr
          XsMac(tid)%XsMacT = XsMac(tid)%XsMacA + XsMac(tid)%XsMacS
        ENDIF
      ELSE
        itype = myFxr%imix
        IF(nTracerCntl%lDynamicBen) THEN
          lCusping = DynMacXsBen(itype)%lCusping
        ELSE
          lCusping = MacXsBen(itype)%lCusping
        END IF
        IF(lCusping) THEN
          phiz(:, tid) = 0.
          philz(:, tid) = 0.
          phiuz(:, tid) = 0.
          volsum = 0.
          DO ig = 1, ng
            DO i = 1, nFsrInFxr
              ifsrlocal = Cell(icel)%MapFxr2FsrIdx(i, j)
              ifsr = FsrIdxSt + ifsrlocal - 1
              vol = Cell(icel)%Vol(ifsrlocal)
              IF(ig .EQ. 1) volsum = volsum + vol
              phiz(ig, tid) = phiz(ig, tid) + FmInfo%phis(ifsr, iz, ig) * vol
              IF(iz .EQ. myzb) THEN
                philz(ig, tid) = philz(ig, tid) + FmInfo%neighphis(ifsr, ig, BOTTOM) * vol
              ELSE
                philz(ig, tid) = philz(ig, tid) + FmInfo%phis(ifsr, iz-1, ig) * vol
              END IF
              IF(iz .EQ. myze) THEN
                phiuz(ig, tid) = phiuz(ig, tid) + FmInfo%neighphis(ifsr, ig, TOP) * vol
              ELSE
                phiuz(ig, tid) = phiuz(ig, tid) + FmInfo%phis(ifsr, iz+1, ig) * vol
              END IF
            END DO
            phiz(ig, tid) = phiz(ig, tid) / volsum
            philz(ig, tid) = philz(ig, tid) / volsum
            phiuz(ig, tid) = phiuz(ig, tid) / volsum
          END DO
          IF(nTracerCntl%lDynamicBen) THEN
            CALL xsbaseDynBen_Cusping(itype, TranInfo%fuelTemp(ipin, iz), 1, ng, 1, ng, lscat1, XsMac(tid),&
              phiz(:,tid), philz(:,tid), phiuz(:,tid), CoreInfo%hzfm(iz), CoreInfo%hzfm(iz-1), CoreInfo%hzfm(iz+1))
          ELSE
            CALL xsbaseBen_Cusping(itype, 1, ng, 1, ng, lscat1, XsMac(tid),&
              phiz(:,tid), philz(:,tid), phiuz(:,tid), CoreInfo%hzfm(iz), CoreInfo%hzfm(iz-1), CoreInfo%hzfm(iz+1))
          END IF
        ELSE
          IF(nTracerCntl%lDynamicBen) THEN
            CALL xsbaseDynBen(itype, TranInfo%fuelTemp(ipin, iz), 1, ng, 1, ng, lscat1, XsMac(tid))
          ELSE
            CALL xsbaseBen(itype, 1, ng, 1, ng, lscat1, XsMac(tid))
          END IF
          IF(nTracerCntl%libtyp .EQ. 11) THEN
            CALL xsbaseben_NEACRP(itype, myFxr%rho, myFxr%temp, myFxr%DopTemp, XsMac(tid))
          ELSE
            IF(nTracerCntl%lDynamicBen) THEN
              CALL xsbaseDynBen(itype, TranInfo%fuelTemp(ipin, iz), 1, ng, 1, ng, lscat1, XsMac(tid))
            ELSE
              CALL xsbaseBen(itype, 1, ng, 1, ng, lscat1, XsMac(tid))
            END IF
          END IF
        END IF
      ENDIF
      CoreXsMac(iz)%XSt(:, ifxr) = XsMac(tid)%XsMacT
      CoreXsMac(iz)%XStr(:, ifxr) = XsMac(tid)%XsMacTr
      CoreXsMac(iz)%XSa(:, ifxr) = XsMac(tid)%XsMacA
      CoreXsMac(iz)%XSnf(:, ifxr) = XsMac(tid)%XsMacNf
      CoreXsMac(iz)%XSkf(:, ifxr) = XsMac(tid)%XsMacKf
      DO ig = 1, ng
        DO igf = InScatRange(1, ig), InScatRange(2, ig)
          igs = InScatIdx(igf, ig)
          CoreXsMac(iz)%XSsm(igs, ifxr) = XsMac(tid)%XsMacSm(igf, ig)
        ENDDO
      ENDDO
    END DO
  END DO
  !$OMP END DO
END DO
!$OMP END PARALLEL

CALL MPI_SYNC(PE%MPI_CMFD_COMM)

DEALLOCATE(phiz, philz, phiuz)

END SUBROUTINE

SUBROUTINE SetMOCtrXS(CoreInfo, Fxr, xst, iz, lsSPH, lsSPHreg, lTrCorrection)
USE TYPEDEF,        ONLY : CoreInfo_Type,       FxrInfo_Type,       Pin_Type,           Cell_Type
USE XSLIB_MOD,      ONLY : igresb,              igrese
USE CORE_MOD,       ONLY : GroupInfo
USE SPH_MOD,        ONLY : ssphfnm
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
REAL, POINTER :: xst(:, :)
INTEGER :: iz
LOGICAL :: lsSPH, lsSPHreg, lTrCorrection

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
REAL, POINTER :: tr(:, :)
REAL :: SPH(100)
INTEGER :: ng, nxy, nLocalFxr, nFsrInFxr
INTEGER :: FsrIdxSt, FxrIdxSt
INTEGER :: i, j, ipin, ifxr, icel, ifsr

Pin => CoreInfo%Pin
Cell => CoreInfo%CellInfo
ng = GroupInfo%ng
nxy = CoreInfo%nxy

IF (lTrCorrection) THEN
  tr => CoreXsMac(iz)%XStr
ELSE
  tr => CoreXsMac(iz)%XSt
ENDIF

!$OMP PARALLEL PRIVATE(ifsr, ifxr, icel, FsrIdxSt, FxrIdxSt, nLocalFxr, nFsrInFxr, SPH)
!$OMP DO SCHEDULE(GUIDED)
DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz); nLocalFxr = Cell(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1
    nFsrInFxr = Cell(icel)%nFsrInFxr(j)
    SPH = 1.0
    IF (lsSPH) THEN
      IF (Cell(icel)%lsSPH) THEN
        IF (lsSPHreg) THEN
          SPH(igresb : igrese) = Fxr(ifxr, iz)%SPHfactor(igresb : igrese)
        ELSE
          SPH(igresb : igrese) = Cell(icel)%SPHfactor(j, igresb : igrese)
        ENDIF
      ENDIF
    ENDIF
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cell(icel)%MapFxr2FsrIdx(i, j) - 1
      xst(1 : ng, ifsr) = tr(1 : ng, ifxr) * SPH(1 : ng)
      IF (lsSPH) ssphfnm(igresb : igrese, ifsr, iz) = SPH(igresb : igrese)
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetMOCPsi(CoreInfo, phis, psi)
USE TYPEDEF,        ONLY : CoreInfo_Type,       Pin_Type,           Cell_Type
USE CORE_MOD,       ONLY : GroupInfo
USE PE_MOD,         ONLY : PE
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phis(:, :, :), psi(:, :)
INTEGER :: iz

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
REAL, POINTER :: nf(:, :)
INTEGER :: ng, nxy, nLocalFxr, nFsrInFxr
INTEGER :: FsrIdxSt, FxrIdxSt
INTEGER :: i, j, ig, ipin, ifxr, icel, ifsr
INTEGER :: myzb, myze

Pin => CoreInfo%Pin
Cell => CoreInfo%CellInfo
ng = GroupInfo%ng
nxy = CoreInfo%nxy
myzb = PE%myzb
myze = PE%myze

psi = 0.0

DO ig = 1, ng
  DO iz = myzb, myze
    nf => CoreXsMac(iz)%XSnf
    !$OMP PARALLEL PRIVATE(ifsr, ifxr, icel, FsrIdxSt, FxrIdxSt, nLocalFxr, nFsrInFxr)
    !$OMP DO SCHEDULE(GUIDED)
    DO ipin = 1, nxy
      FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
      icel = Pin(ipin)%Cell(iz); nLocalFxr = Cell(icel)%nFxr
      DO j = 1, nLocalFxr
        ifxr = FxrIdxSt + j - 1
        nFsrInFxr = Cell(icel)%nFsrInFxr(j)
        DO i = 1, nFsrInFxr
          ifsr = FsrIdxSt + Cell(icel)%MapFxr2FsrIdx(i, j) - 1
          psi(ifsr, iz) = psi(ifsr, iz) + nf(ig, ifxr) * phis(ifsr, iz, ig)
        ENDDO
      ENDDO
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE SetMOCSource(CoreInfo, Fxr, src, phis, psi, AxSrc, xst, eigv, iz, gb, ge, Offset)
USE TYPEDEF,        ONLY : CoreInfo_Type,       FxrInfo_Type,       Pin_Type,           Cell_Type
USE CORE_MOD,       ONLY : GroupInfo
USE CNTL,           ONLY : nTracerCntl
USE BenchXs,        ONLY : GetChiBen,           GetChiDynBen
USE TRAN_MOD,       ONLY : TranInfo
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(FxrInfo_Type) :: Fxr(:)
REAL, POINTER :: src(:, :), phis(:, :), psi(:, :), AxSrc(:, :, :), xst(:, :)
REAL :: eigv
INTEGER :: iz, gb, ge
INTEGER, OPTIONAL :: Offset

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
REAL, POINTER :: tot(:, :), tr(:, :), sm(:, :)
REAL :: Chi(gb : ge)
REAL :: sigs, reigv
INTEGER :: nxy, nLocalFxr, nFsrInFxr, nChi
INTEGER :: FsrIdxSt, FxrIdxSt
INTEGER :: i, j, ipin, ifxr, icel, ifsr, itype, ig, igs, igf
LOGICAL :: l3dim, lxslib, lscat1

INTEGER :: off

Pin => CoreInfo%Pin
Cell => CoreInfo%CellInfo
tot => CoreXsMac(iz)%XSt
tr => CoreXsMac(iz)%XStr
sm => CoreXsMac(iz)%XSsm
nxy = CoreInfo%nxy
nChi = GroupInfo%nChi
l3dim = nTracerCntl%l3dim
lxslib = nTracerCntl%lxslib
lscat1 = nTracerCntl%lscat1
reigv = 1.0 / eigv

off = 0
IF (PRESENT(Offset)) off = Offset

!$OMP PARALLEL PRIVATE(ifsr, ifxr, icel, itype, igs, FsrIdxSt, FxrIdxSt, nLocalFxr, nFsrInFxr, Chi, sigs)
!$OMP DO SCHEDULE(GUIDED)
DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz); nLocalFxr = Cell(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1
    nFsrInFxr = Cell(icel)%nFsrInFxr(j)
    IF (lxslib) THEN
      DO ig = gb, ge
        IF (ig .GT. nChi) THEN
          Chi(ig) = 0
        ELSE
          Chi(ig) = 0
          IF (Fxr(ifxr)%ldepl) Chi(ig) = Fxr(ifxr)%Chi(ig)
        ENDIF
      ENDDO
    ELSE
      itype = Fxr(ifxr)%imix
      IF(nTracerCntl%lDynamicBen) THEN
        Chi(gb : ge) = GetChiDynBen(itype, TranInfo%fuelTemp(ipin, iz), gb, ge)
      ELSE
        Chi(gb : ge) = GetChiBen(itype, gb, ge)
      ENDIF
    ENDIF
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cell(icel)%MapFxr2FsrIdx(i, j) - 1
      DO ig = gb, ge
        src(ig - off, ifsr) = reigv * Chi(ig) * psi(ifsr, iz)
        DO igf = InScatRange(1, ig), InScatRange(2, ig)
          igs = InScatIdx(igf, ig)
          sigs = sm(igs, ifxr)
          IF (igf .EQ. ig .AND. lScat1) sigs = sigs + (tot(ig, ifxr) - tr(ig, ifxr))
          src(ig - off, ifsr) = src(ig - off, ifsr) + sigs * phis(igf, ifsr)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO

IF (l3dim) THEN
  !$OMP DO SCHEDULE(GUIDED)
  DO ipin = 1, nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
    icel = Pin(ipin)%Cell(iz); nLocalFxr = Cell(icel)%nFxr
    IF (nTracerCntl%LkgSplitLv .EQ. 0) THEN
      DO j = 1, nLocalFxr
        ifxr = FxrIdxSt + j - 1
        nFsrInFxr = Cell(icel)%nFsrInFxr(j)
        DO i = 1, nFsrInFxr
          ifsr = FsrIdxSt + Cell(icel)%MapFxr2FsrIdx(i, j) - 1
          DO ig = gb, ge
            IF (AxSrc(ipin, iz, ig) .LT. 0.0 .AND. .NOT. Fxr(ifxr)%lvoid) THEN
              src(ig - off, ifsr) = src(ig - off, ifsr) - AxSrc(ipin, iz, ig)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ELSE
      DO j = 1, nLocalFxr
        ifxr = FxrIdxSt + j - 1
        nFsrInFxr = Cell(icel)%nFsrInFxr(j)
        DO i = 1, nFsrInFxr
          ifsr = FsrIdxSt + Cell(icel)%MapFxr2FsrIdx(i, j) - 1
          DO ig = gb, ge
            src(ig - off, ifsr) = src(ig - off, ifsr) - AxSrc(ipin, iz, ig)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDDO
  !$OMP END DO
ENDIF
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(ifsr, icel, FsrIdxSt)
!$OMP DO SCHEDULE(GUIDED)
DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt
  icel = Pin(ipin)%Cell(iz)
  DO j = 1, Cell(icel)%nFsr
    ifsr = FsrIdxSt + j - 1
    src(gb - off : ge - off, ifsr) = src(gb - off : ge - off, ifsr) / xst(gb : ge, ifsr)
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetMOCSourcePN(CoreInfo, Fxr, srcm, phim, xst, iz, gb, ge, ScatOd, Offset)
USE TYPEDEF,        ONLY : CoreInfo_Type,       FxrInfo_Type,       Pin_Type,           Cell_Type
USE CORE_MOD,       ONLY : GroupInfo
USE CNTL,           ONLY : nTracerCntl
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(FxrInfo_Type) :: Fxr(:)
REAL, POINTER :: srcm(:, :, :), phim(:, :, :), xst(:, :)
INTEGER :: iz, gb, ge
INTEGER :: ScatOd
INTEGER, OPTIONAL :: Offset

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
REAL, POINTER :: smP1(:, :), smP2(:, :), smP3(:, :)
INTEGER :: nxy, nLocalFxr, nFsrInFxr
INTEGER :: FsrIdxSt, FxrIdxSt
INTEGER :: i, j, ipin, ifxr, icel, ifsr, itype, ig, igs, igf
LOGICAL :: lxslib
REAL :: scatSrc(9), xstinv

INTEGER :: off

Pin => CoreInfo%Pin
Cell => CoreInfo%CellInfo
smP1 => CoreXsMac(iz)%XSsmP1
smP2 => CoreXsMac(iz)%XSsmP2
smP3 => CoreXsMac(iz)%XSsmP3
nxy = CoreInfo%nxy

off = 0
IF (PRESENT(Offset)) off = Offset

!$OMP PARALLEL PRIVATE(ifsr, ifxr, icel, itype, igs, FsrIdxSt, FxrIdxSt, nLocalFxr, nFsrInFxr, scatSrc)
!$OMP DO SCHEDULE(GUIDED)
DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz); nLocalFxr = Cell(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1
    nFsrInFxr = Cell(icel)%nFsrInFxr(j)
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cell(icel)%MapFxr2FsrIdx(i, j) - 1
      IF (ScatOd .EQ. 1) THEN
        DO ig = gb, ge
          scatSrc = 0.0
          DO igf = InScatRange(1, ig), InScatRange(2, ig)
            igs = InScatIdx(igf, ig)
            scatSrc(1:2) = scatSrc(1:2) + smP1(igs, ifxr) * phim(1:2, igf, ifsr)
            srcm(1:2, ig - off, ifsr) = scatSrc(1:2)
          ENDDO
        ENDDO
      ELSEIF (ScatOd .EQ. 2) THEN
        DO ig = gb, ge
          scatSrc = 0.0
          DO igf = InScatRange(1, ig), InScatRange(2, ig)
            igs = InScatIdx(igf, ig)
            scatSrc(1:2) = scatSrc(1:2) + smP1(igs, ifxr) * phim(1:2, igf, ifsr)
            scatSrc(3:5) = scatSrc(3:5) + smP2(igs, ifxr) * phim(3:5, igf, ifsr)
            srcm(1:5, ig - off, ifsr) = scatSrc(1:5)
          ENDDO
        ENDDO
      ELSEIF (ScatOd .EQ. 3) THEN
        DO ig = gb, ge
          scatSrc = 0.0
          DO igf = InScatRange(1, ig), InScatRange(2, ig)
            igs = InScatIdx(igf, ig)
            scatSrc(1:2) = scatSrc(1:2) + smP1(igs, ifxr) * phim(1:2, igf, ifsr)
            scatSrc(3:5) = scatSrc(3:5) + smP2(igs, ifxr) * phim(3:5, igf, ifsr)
            scatSrc(6:9) = scatSrc(6:9) + smP3(igs, ifxr) * phim(6:9, igf, ifsr)
            srcm(1:9, ig - off, ifsr) = scatSrc(1:9)
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(ifsr, icel, FsrIdxSt, xstinv)
!$OMP DO SCHEDULE(GUIDED)
DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt
  icel = Pin(ipin)%Cell(iz)
  DO j = 1, Cell(icel)%nFsr
    ifsr = FsrIdxSt + j - 1
    IF (ScatOd .EQ. 1) THEN
      DO ig = gb, ge
        xstinv = 1.0 / xst(ig, ifsr)
        srcm(1:2, ig - off, ifsr) = 3.0 * srcm(1:2, ig - off, ifsr) * xstinv
      ENDDO
    ELSEIF (ScatOd .EQ. 2) THEN
      DO ig = gb, ge
        xstinv = 1.0 / xst(ig, ifsr)
        srcm(1:2, ig - off, ifsr) = 3.0 * srcm(1:2, ig - off, ifsr) * xstinv
        srcm(3:5, ig - off, ifsr) = 5.0 * srcm(3:5, ig - off, ifsr) * xstinv
      ENDDO
    ELSEIF (ScatOd .EQ. 3) THEN
      DO ig = gb, ge
        xstinv = 1.0 / xst(ig, ifsr)
        srcm(1:2, ig - off, ifsr) = 3.0 * srcm(1:2, ig - off, ifsr) * xstinv
        srcm(3:5, ig - off, ifsr) = 5.0 * srcm(3:5, ig - off, ifsr) * xstinv
        srcm(6:9, ig - off, ifsr) = 7.0 * srcm(6:9, ig - off, ifsr) * xstinv
      ENDDO
    ENDIF
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetMOCPower(CoreInfo, phis, Power)
USE TYPEDEF,        ONLY : CoreInfo_Type,       Pin_Type,           Cell_Type
USE CORE_MOD,       ONLY : GroupInfo
USE PE_MOD,         ONLY : PE
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phis(:, :, :), Power(:, :)

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
INTEGER :: ng, nxy, nFsr, nLocalFxr, nFsrInFxr
INTEGER :: FsrIdxSt, FxrIdxSt
INTEGER :: i, j, ig, icel, ifsr, ifxr, ipin, iz
INTEGER :: myzb, myze

Pin => CoreInfo%Pin
Cell => CoreInfo%CellInfo
ng = GroupInfo%ng
nxy = CoreInfo%nxy
nFsr = CoreInfo%nCoreFsr
myzb = PE%myzb
myze = PE%myze

Power = 0.0

!$OMP PARALLEL PRIVATE(icel, ifsr, ifxr, FsrIdxSt, FxrIdxSt, nLocalFxr, nFsrInFxr)
DO ig = 1, ng
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
  DO iz = myzb, myze
    DO ipin = 1, nxy
      FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
      icel = Pin(ipin)%Cell(iz); nLocalFxr = Cell(icel)%nFxr
      DO j = 1, nLocalFxr
        ifxr = FxrIdxSt + j - 1
        nFsrInFxr = Cell(icel)%nFsrInFxr(j)
        DO i = 1, nFsrInFxr
          ifsr = FsrIdxSt + Cell(icel)%MapFxr2FsrIdx(i, j) - 1
          Power(ifsr, iz) = Power(ifsr, iz) + CoreXsMac(iz)%XSkf(ig, ifxr) * phis(ifsr, iz, ig)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
ENDDO
!$OMP END PARALLEL

END SUBROUTINE

#ifdef __PGI

SUBROUTINE CopyXS(CoreInfo, FmInfo, iz)
USE TYPEDEF,        ONLY : CoreInfo_Type,       FmInfo_Type,         Pin_Type,           Cell_Type,                 &
                           FxrInfo_Type
USE CORE_MOD,       ONLY : GroupInfo
USE CNTL,           ONLY : nTracerCntl
USE BenchXs,        ONLY : GetChiBen,           getChiDynBen
USE CUDA_MASTER
USE TRAN_MOD,       ONLY : TranInfo
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(FxrInfo_Type), POINTER :: Fxr(:)
INTEGER :: ng, nxy, nfxr, nChi, nLocalFxr, FxrIdxSt
INTEGER :: j, ig, igf, igs, ifxr, icel, ipin, iz
REAL, ALLOCATABLE :: xssm(:, :, :), chi(:, :)
LOGICAL :: lxslib

Pin => CoreInfo%Pin
Cell => CoreInfo%CellInfo
Fxr => FmInfo%Fxr(:, iz)
ng = cuGeometry%ng
nxy = cuGeometry%nxy
nfxr = cuGeometry%nfxr
nChi = GroupInfo%nChi
lxslib = nTracerCntl%lxslib

ALLOCATE(xssm(ng, ng, nfxr)); xssm = 0.0
ALLOCATE(chi(ng, nfxr)); chi = 0.0

!$OMP PARALLEL PRIVATE(icel, ifxr, igs, FxrIdxSt, nLocalFxr)
!$OMP DO SCHEDULE(GUIDED)
DO ipin = 1, nxy
  FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz); nLocalFxr = Cell(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1
    IF (lxslib) THEN
      IF (Fxr(ifxr)%ldepl) chi(1 : nChi, ifxr) = Fxr(ifxr)%chi
    ELSE
      IF(nTracerCntl%lDynamicBen) THEN
        chi(:, ifxr) = GetChiDynBen(Fxr(ifxr)%imix, TranInfo%fuelTemp(ipin, iz), 1, ng)
      ELSE
        chi(:, ifxr) = GetChiBen(Fxr(ifxr)%imix, 1, ng)
      ENDIF
    ENDIF
    DO ig = 1, ng
      DO igf = InScatRange(1, ig), InScatRange(2, ig)
        igs = InScatIdx(igf, ig)
        xssm(ig, igf, ifxr) = CoreXsMac(iz)%XSsm(igs, ifxr)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

cuMOC%xssm = xssm
cuMOC%chi = chi

DEALLOCATE(xssm, chi)

END SUBROUTINE

SUBROUTINE SetMOCPsi_iter(CoreInfo, phis, psi, iz)
USE TYPEDEF,        ONLY : CoreInfo_Type,       Pin_Type,           Cell_Type
USE CORE_MOD,       ONLY : GroupInfo
USE PE_MOD,         ONLY : PE
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phis(:, :), psi(:, :)
INTEGER :: iz

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
REAL, POINTER :: XSnf(:, :)
INTEGER :: ng, nxy, nLocalFxr, nFsrInFxr
INTEGER :: FsrIdxSt, FxrIdxSt
INTEGER :: i, j, ig, ipin, ifxr, icel, ifsr
INTEGER :: myzb, myze

Pin => CoreInfo%Pin
Cell => CoreInfo%CellInfo
ng = GroupInfo%ng
nxy = CoreInfo%nxy
myzb = PE%myzb
myze = PE%myze

psi(:,iz) = 0.0

XSnf => CoreXsMac(iz)%XSnf
!$OMP PARALLEL PRIVATE(ifsr, ifxr, icel, FsrIdxSt, FxrIdxSt, nLocalFxr, nFsrInFxr)
!$OMP DO SCHEDULE(GUIDED)
DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz); nLocalFxr = Cell(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1
    nFsrInFxr = Cell(icel)%nFsrInFxr(j)
    DO i = 1, nFsrInFxr
      DO ig = 1, ng
        ifsr = FsrIdxSt + Cell(icel)%MapFxr2FsrIdx(i, j) - 1
        psi(ifsr, iz) = psi(ifsr, iz) + XSnf(ig, ifxr) * phis(ig, ifsr)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE
#endif

END MODULE
