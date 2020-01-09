MODULE MOC_COMMON

USE TYPEDEF,        ONLY : XsMac_Type
USE TYPEDEF_COMMON
USE OMP_LIB
IMPLICIT NONE

TYPE(XsMac_Type), POINTER, PRIVATE :: XsMac(:)
TYPE(CoreXsMac_Type), ALLOCATABLE :: CoreXsMac(:)
INTEGER, ALLOCATABLE :: InScatRange(:, :), InScatIdx(:, :)

CONTAINS

SUBROUTINE AllocCoreMacXs(CoreInfo)
USE TYPEDEF,        ONLY : CoreInfo_Type
USE PE_MOD,         ONLY : PE
USE CORE_MOD,       ONLY : GroupInfo
USE XsUtil_mod,     ONLY : AllocXsMac
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
INTEGER :: ng, ngs, nFxr
INTEGER :: iz, ig, igs, tid
INTEGER :: myzb, myze

ng = GroupInfo%ng
nFxr = CoreInfo%nCoreFxr
myzb = PE%myzb
myze = PE%myze

ALLOCATE(InScatRange(2, ng), InScatIdx(ng, ng))

InScatRange = GroupInfo%InScatRange
InScatIdx = 0

ngs = 1; InScatIdx(1, 1) = 1
DO ig = 2, ng
  DO igs = InScatRange(1, ig), InScatRange(2, ig)
    ngs = ngs + 1
    InScatIdx(igs, ig) = ngs
  ENDDO
ENDDO

ALLOCATE(XsMac(100))

DO tid = 1, 100
  XsMac(tid)%ng = ng
  CALL AllocXsMac(XsMac(tid))
ENDDO

ALLOCATE(CoreXsMac(myzb : myze))

DO iz = myzb, myze
  ALLOCATE(CoreXsMac(iz)%XSt(ng, nFxr))
  ALLOCATE(CoreXsMac(iz)%XStr(ng, nFxr))
  ALLOCATE(CoreXsMac(iz)%XSa(ng, nFxr))
  ALLOCATE(CoreXsMac(iz)%XSnf(ng, nFxr))
  ALLOCATE(CoreXsMac(iz)%XSkf(ng, nFxr))
  ALLOCATE(CoreXsMac(iz)%XSsm(ngs, nFxr))
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
USE BenchXs,        ONLY : xsbaseBen
USE MacXsLib_Mod,   ONLY : MacXsBase,           MacXsScatMatrix,    BaseMacStr
USE MPIComm_Mod,    ONLY : MPI_SYNC
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
      IF (lTrCorrection) THEN
        DO ig = 1, ng
          CALL BaseMacStr(XsMac(tid), myFxr, ig, ng, .TRUE.)
        ENDDO
      ENDIF
      IF (myFxr%lres) THEN
        XsMac(tid)%XsMacA(irgb : irge) = myFxr%FresoA(irgb : irge) * XsMac(tid)%XsMacA(irgb : irge)
        IF (CoreInfo%lFuelPlane(iz)) THEN
          XsMac(tid)%XsMacNf(irgb : irge) = myFxr%FresoF(irgb : irge) * XsMac(tid)%XsMacNf(irgb : irge)
          XsMac(tid)%XsMacKf(irgb : irge) = myFxr%FresoF(irgb : irge) * XsMac(tid)%XsMacKf(irgb : irge)
        ENDIF
      ENDIF
      XsMac(tid)%XsMacTr = XsMac(tid)%XsMacA + XsMac(tid)%XsMacStr
      XsMac(tid)%XsMacT = XsMac(tid)%XsMacA + XsMac(tid)%XsMacS
    ELSE
      itype = myFxr%imix
      IF (itype .EQ. 0) CYCLE
      CALL xsbaseBen(itype, 1, ng, 1, ng, lscat1, XsMac(tid))
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
  ENDDO
  !$OMP END DO
ENDDO

!$OMP END PARALLEL

CALL MPI_SYNC(PE%MPI_CMFD_COMM)

END SUBROUTINE

SUBROUTINE SetMOCtrXS(CoreInfo, xst, iz)
USE TYPEDEF,        ONLY : CoreInfo_Type,       Pin_Type,           Cell_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: xst(:, :)
INTEGER :: iz

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
REAL, POINTER :: XStr(:, :)
INTEGER :: nxy, nLocalFxr, nFsrInFxr
INTEGER :: FsrIdxSt, FxrIdxSt
INTEGER :: i, j, ipin, ifxr, icel, ifsr

Pin => CoreInfo%Pin
Cell => CoreInfo%CellInfo
XStr => CoreXsMac(iz)%XStr
nxy = CoreInfo%nxy

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
      xst(:, ifsr) = XStr(:, ifxr)
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

psi = 0.0

DO ig = 1, ng
  DO iz = myzb, myze
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
          ifsr = FsrIdxSt + Cell(icel)%MapFxr2FsrIdx(i, j) - 1
          psi(ifsr, iz) = psi(ifsr, iz) + XSnf(ig, ifxr) * phis(ifsr, iz, ig)
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
USE BenchXs,        ONLY : GetChiBen
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(FxrInfo_Type) :: Fxr(:)
REAL, POINTER :: src(:, :), phis(:, :), psi(:, :), AxSrc(:, :, :), xst(:, :)
REAL :: eigv
INTEGER :: iz, gb, ge
INTEGER, OPTIONAL :: Offset

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
REAL, POINTER :: XSsm(:, :)
REAL :: Chi(gb : ge)
REAL :: reigv
INTEGER :: nxy, nLocalFxr, nFsrInFxr, nChi
INTEGER :: FsrIdxSt, FxrIdxSt
INTEGER :: i, j, ipin, ifxr, icel, ifsr, itype, ig, igs, igf
LOGICAL :: l3dim, lxslib, lscat1

INTEGER :: off

Pin => CoreInfo%Pin
Cell => CoreInfo%CellInfo
XSsm => CoreXsMac(iz)%XSsm
nxy = CoreInfo%nxy
nChi = GroupInfo%nChi
l3dim = nTracerCntl%l3dim
lxslib = nTracerCntl%lxslib
lscat1 = nTracerCntl%lscat1
reigv = 1.0 / eigv

off = 0
IF (PRESENT(Offset)) off = Offset

!$OMP PARALLEL PRIVATE(ifsr, ifxr, icel, itype, igs, FsrIdxSt, FxrIdxSt, nLocalFxr, nFsrInFxr, Chi)
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
      Chi(gb : ge) = GetChiBen(itype, gb, ge)
    ENDIF
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cell(icel)%MapFxr2FsrIdx(i, j) - 1
      DO ig = gb, ge
        src(ig - off, ifsr) = reigv * Chi(ig) * psi(ifsr, iz)
        DO igf = InScatRange(1, ig), InScatRange(2, ig)
          igs = InScatIdx(igf, ig)
          src(ig - off, ifsr) = src(ig - off, ifsr) + XSsm(igs, ifxr) * phis(igf, ifsr)
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

SUBROUTINE SetMOCPower(CoreInfo, phis, Power)
USE TYPEDEF,        ONLY : CoreInfo_Type,       Pin_Type,           Cell_Type
USE CORE_MOD,       ONLY : GroupInfo
USE PE_MOD,         ONLY : PE
USE MPIComm_Mod,    ONLY : BCAST
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phis(:, :, :), Power(:, :)

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
INTEGER :: ng, nxy, nLocalFxr, nFsrInFxr
INTEGER :: FsrIdxSt, FxrIdxSt
INTEGER :: i, j, ig, icel, ifsr, ifxr, ipin, iz
INTEGER :: myzb, myze

Pin => CoreInfo%Pin
Cell => CoreInfo%CellInfo
ng = GroupInfo%ng
nxy = CoreInfo%nxy
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

DO iz = 1, CoreInfo%nz
  CALL BCAST(Power(:, iz), CoreInfo%nCoreFsr, PE%MPI_RTMASTER_COMM, PE%AxDomList(iz))
ENDDO

END SUBROUTINE

#ifdef __PGI

SUBROUTINE CopyXS(CoreInfo, FmInfo, iz)
USE TYPEDEF,        ONLY : CoreInfo_Type,       FmInfo_Type,         Pin_Type,           Cell_Type,                 &
                           FxrInfo_Type
USE CORE_MOD,       ONLY : GroupInfo
USE CNTL,           ONLY : nTracerCntl
USE BenchXs,        ONLY : GetChiBen
USE CUDA_MASTER
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
      chi(:, ifxr) = GetChiBen(Fxr(ifxr)%imix, 1, ng)
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

#endif

END MODULE