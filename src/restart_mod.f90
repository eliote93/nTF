#include <defines.h>
MODULE Restart_mod
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,         FmInfo_Type,           PE_TYPE,             &
                           FxrInfo_Type,          AsyInfo_Type,          Asy_Type,            &
                           Pin_Type,              PinInfo_Type,          Cell_Type,           &
                           LpShf_Type
USE files,          ONLY : io8,                   io13,                  io14,                &
                           caseid,                localfn,               RstFileIdx,          &
                           DeplRstFileIdx,        Filename
USE ioutil,         ONLY : openfile,              message,               terminate
USE CNTL,           ONLY : nTracerCntl_Type
USE BasicOperation, ONLY : CP_CA
#ifdef MPI_ENV
USE MPIComm_Mod,    ONLY : MPIWaitTurn
#endif
IMPLICIT NONE
SAVE
TYPE(FxrInfo_Type), PRIVATE, POINTER :: Fxr(:, :)
TYPE(AsyInfo_Type), PRIVATE, POINTER :: AsyType(:)
TYPE(Asy_Type), PRIVATE, POINTER :: Asy(:)
TYPE(Pin_Type), PRIVATE, POINTER :: Pin(:)
TYPE(Cell_Type), PRIVATE, POINTER :: CellInfo(:)

INTEGER, PRIVATE, POINTER :: CoreIdx(:, :), CoreMap(:)
INTEGER, PRIVATE :: nxc, nyc

CONTAINS

SUBROUTINE WriteRestartFile(Core, FmInfo, LDeplRst, istep, nTRACERCntl, PE)
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(nTracerCntl_Type) :: nTRACERCntl
TYPE(PE_TYPE) :: PE
INTEGER :: istep
LOGICAL :: LDeplRst

INTEGER :: ixya, ixa, iya, iasy, iz
INTEGER :: iasytype
LOGICAL :: LPartialAsy
!File IO
INTEGER :: nchar, nproc, nasywrite
CHARACTER(256) :: Fn, Fn0, AsyFn
!Assembly ID
!
!MPI Variables
INTEGER :: irank
LOGICAL :: Master, Slave

Master =  PE%CmfdMaster; Slave = PE%CmfdSlave
!
Fxr => FmInfo%Fxr; Pin => Core%Pin
AsyType => Core%AsyInfo; Asy => Core%Asy
CoreIdx => Core%CoreIdx; CoreMap => Core%CoreMap
CellInfo => Core%CellInfo

nxc = Core%nxc; nyc = Core%nyc
!DO nchar = 1, 80
!  if(caseid(nchar:nchar) .EQ. '') EXIT
!ENDDO
!nchar = nchar - 1
!nproc = PE%nCmfdProc
!WRITE(Fn, '(A,A,A)') caseid(1:nchar),'.RST'   !Restart file
!


IF(.NOT. LDeplRst) THEN
  Fn = filename(RstFileIdx)
  DO nchar = 1, 80
    if(Fn(nchar:nchar) .EQ. '') EXIT
  ENDDO
  nchar = nchar - 1
ELSE
  DO nchar = 1, 80
    if(caseid(nchar:nchar) .EQ. '') EXIT
  ENDDO
  nchar = nchar - 1
  Fn0(1:nchar) = caseid(1:nchar)
  IF(istep .LT. 10) THEN
    WRITE(fn0(nchar+1:256),'(A, A, I1)') '_cycle', '00', istep
  ELSEIF(istep .LT. 100) THEN
    WRITE(fn0(nchar+1:256),'(A, A, I2)') '_cycle', '0', istep
  ELSE
    WRITE(fn0(nchar+1:256),'(A, I3)') '_cycle', istep
  ENDIF
  nchar = nchar + 9
  !nproc = PE%nCmfdProc
  WRITE(Fn, '(A,A,A)') Fn0(1:nchar),'.DEPLRST'   !Restart file
!
  !Fn = filename(DeplRstFileIdx)
ENDIF
WRITE(mesg, '(A, A)')  'Writing Restart File : ', Fn(1:nchar+10)
IF(master) CALL message(io8, TRUE, TRUE, mesg)
#ifdef MPI_ENV
CALL MPIWaitTurn(PE%MPI_Cmfd_Comm, PE%myCmfdRank, PE%nCmfdproc, .FALSE.)
#endif
IF(Master) THEN
  CALL Openfile(io14, FALSE, TRUE, FALSE, Fn)
  WRITE(io14) Core%nxa, Core%nya, Core%nxya, Core%nz, Core%Hm_Mass0
  CLOSE(io14)
ENDIF
CALL Openfile(io14, FALSE, TRUE, TRUE, Fn)
DO iz = PE%myzb, PE%myze
  IF(.NOT. Core%lFuelPlane(iz)) CYCLE
  nasywrite  = 0
  DO iya = 1, Core%nya
    DO ixa = 1, Core%nxa
      iasy = CoreIdx(ixa, iya) !Assembly Idx
      IF(iasy .EQ. 0) CYCLE
      iasytype = CoreMap(iasy)
      IF(.NOT. AsyType(iasytype)%lFuel) CYCLE
      nasywrite  =  nasywrite  + 1
    ENDDO
  ENDDO
  WRITE(io14) iz, nasywrite
  DO iya = 1, Core%nya
    DO ixa = 1, Core%nxa
      iasy = CoreIdx(ixa, iya) !Assembly Idx
      IF(iasy .EQ. 0) CYCLE
      iasytype = CoreMap(iasy) !Assembly Type
      IF(.NOT. AsyType(iasytype)%lFuel) CYCLE
        LPartialAsy = .FALSE.
        LPartialAsy = Asy(iasy)%lCentX .OR. Asy(iasy)%lCentY .OR. Asy(iasy)%lCentXY
      IF(.NOT. lPartialAsy) THEN
        CALL WriteAsyRstData(iasytype, iasy, ixa, iya, iz)
      ELSE
        CALL WritePartialAsyRstData(iasytype, iasy, ixa, iya, iz)
      ENDIF
    ENDDO
  ENDDO
ENDDO
CLOSE(io14)
#ifdef MPI_ENV
CALL MPIWaitTurn(PE%MPI_Cmfd_Comm, PE%myCmfdRank, PE%nCmfdproc, .TRUE.)
#endif
!Fxr => FmInfo%Fxr(:, :); Pin => Core%Pin
!AsyType => Core%AsyInfo; Asy => Core%Asy
!CoreIdx => Core%CoreIdx; CoreMap => Core%CoreMap
!CellInfo => Core%CellInfo
NULLIFY(Fxr, Pin, AsyType, Asy)
NULLIFY(CoreIdx, CoreMap, CellInfo)

END SUBROUTINE
!
SUBROUTINE  WritePartialAsyRstData(iasytype, iasy, ixa, iya, iz)
INTEGER :: iasytype, iasy, ixa, iya, iz
INTEGER :: ix, ix0, iy, iy0, ixy, ipin, icel, ifxr
INTEGER :: nxy, nwrite, nwrite0
INTEGER :: FxrIdxSt, nLocalFxr, nDeplFxr
INTEGER :: i, j
INTEGER :: PinMapExt(100, 100)
LOGICAL :: lCentX, lCentY, lCentXY

lCentX = AsyType(iAsyType)%lCentX; lCentY = AsyType(iAsyType)%lCentY
lCentXY = AsyType(iAsyType)%lCentXY
iy0 = 0
CALL CP_CA(PinMapExt, 0, 100, 100)
DO iy = nyc-AsyType(iAsyType)%ny + 1, nyc
  iy0 = iy0 + 1; ix0 = 0
  DO ix = nxc-AsyType(iAsyType)%nx + 1, nxc
    ix0 = ix0 + 1
    ixy = AsyType(iAsyType)%Pin2DIdx(ix0, iy0)
    PinMapExt(ix, iy) = ixy
  ENDDO
ENDDO

DO iy = 1, nyc
  DO ix = 1, nxc
    IF(PinMapExt(ix, iy) .NE. 0) CYCLE
    ix0 = ix; iy0 = iy
    IF(lCentX .OR. lCentXY) THEN
      iy0 = nyc - iy0 + 1; ix0 = ix0
      PinMapExt(ix, iy) = PinMapExt(ix0, iy0)
    ENDIF
    IF(PinMapExt(ix, iy) .NE. 0) CYCLE
    IF(lCentY .OR. lCentXY) THEN
      iy0 = iy0; ix0 = nxc - ix0 + 1
      PinMapExt(ix, iy) = PinMapExt(ix0, iy0)
    ENDIF
  ENDDO
  !WRITE(88, '(100I5)') PinMapExt(1:nxc, iy)
ENDDO
 !WRITE(88, '(100I5)')
nxy = AsyType(iasytype)%nxy
!Calculate Number of fuel fxr in assembly
nwrite = 0
DO iy = 1, nyc
  DO ix = 1, nxc
    ixy = PinMapExt(ix, iy)
    ipin = Asy(iasy)%GlobalPinIdx(ixy)
    icel = Pin(ipin)%Cell(iz)
    FxrIdxSt = Pin(ixy)%FxrIdxSt
    nLocalFxr = CellInfo(icel)%nFxr
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      IF(.NOT. Fxr(ifxr, iz)%lDepl) CYCLE
      nwrite = nwrite + 1
    ENDDO
  ENDDO
ENDDO

WRITE(io14) ixa, iya, iz, nxc * nyc, nwrite

DO iy = 1, nyc
  DO ix = 1, nxc
    ixy = PinMapExt(ix, iy)
    ipin = Asy(iasy)%GlobalPinIdx(ixy)
    icel = Pin(ipin)%Cell(iz)
    FxrIdxSt = Pin(ipin)%FxrIdxSt
    nLocalFxr = CellInfo(icel)%nFxr
    nDeplFxr = 0
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      IF(Fxr(ifxr, iz)%lDepl) nDeplFxr = nDeplFxr + 1
    ENDDO
    WRITE(io14) ix, iy, nLocalFxr, nDeplFxr
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      WRITE(io14) ifxr, j, Fxr(ifxr, iz)%lfuel, Fxr(ifxr, iz)%ldepl, Fxr(ifxr, iz)%lres, Fxr(ifxr, iz)%lGd, Fxr(ifxr, iz)%lh2o
      IF(.NOT. Fxr(ifxr, iz)%lDepl) CYCLE
      IF(Fxr(ifxr,iz)%niso .GT. Fxr(ifxr, iz)%niso_depl) Fxr(ifxr,iz)%niso_depl = Fxr(ifxr, iz)%niso
      WRITE(io14) Fxr(ifxr, iz)%niso, Fxr(ifxr, iz)%niso_depl, Fxr(ifxr, iz)%temp, Fxr(ifxr, iz)%burnup
      DO i = 1, Fxr(ifxr, iz)%niso_depl
        WRITE(io14) Fxr(ifxr, iz)%idiso(i), Fxr(ifxr, iz)%pnum(i)
      ENDDO
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE


SUBROUTINE WriteAsyRstData(iasytype, iasy, ixa, iya, iz)
INTEGER :: iasytype, iasy, ixa, iya, iz
INTEGER :: ixy, ix, iy, ipin, icel, ifxr
INTEGER :: nxy, nwrite
INTEGER :: FxrIdxSt, nLocalFxr, nLocalFxrWrite, nDeplFxr
INTEGER :: i, j


nxy = AsyType(iasytype)%nxy
!Calculate Number of fuel fxr in assembly
nwrite = 0
DO ixy = 1, nxy
  ipin = Asy(iasy)%GlobalPinIdx(ixy)
  icel = Pin(ipin)%Cell(iz)
  FxrIdxSt = Pin(ipin)%FxrIdxSt
  nLocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    IF(.NOT. Fxr(ifxr, iz)%lDepl) CYCLE
    nwrite = nwrite + 1
  ENDDO
ENDDO
WRITE(io14) ixa, iya, iz, nxy, nwrite

DO iy = 1, nyc
  DO ix = 1, nxc
    ixy = AsyType(iAsyType)%Pin2DIdx(ix, iy)
    ipin = Asy(iasy)%GlobalPinIdx(ixy)
    icel = Pin(ipin)%Cell(iz)
    FxrIdxSt = Pin(ipin)%FxrIdxSt
    nLocalFxr = CellInfo(icel)%nFxr
    nDeplFxr = 0
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      IF(Fxr(ifxr, iz)%lDepl) nDeplFxr = nDeplFxr + 1
    ENDDO
    WRITE(io14) ix, iy, nLocalFxr, nDeplFxr
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      WRITE(io14) ifxr, j, Fxr(ifxr, iz)%lfuel, Fxr(ifxr, iz)%ldepl, Fxr(ifxr, iz)%lres,                 &
                  Fxr(ifxr, iz)%lGd, Fxr(ifxr, iz)%lh2o
      IF(.NOT. Fxr(ifxr, iz)%lDepl) CYCLE
      IF(Fxr(ifxr,iz)%niso .GT. Fxr(ifxr, iz)%niso_depl) Fxr(ifxr,iz)%niso_depl = Fxr(ifxr, iz)%niso
      WRITE(io14) Fxr(ifxr, iz)%niso, Fxr(ifxr, iz)%niso_depl, Fxr(ifxr, iz)%temp, Fxr(ifxr, iz)%burnup
      DO i = 1, Fxr(ifxr, iz)%niso_depl
        WRITE(io14) Fxr(ifxr, iz)%idiso(i), Fxr(ifxr, iz)%pnum(i)
      ENDDO
    ENDDO
  ENDDO
ENDDO
continue
END SUBROUTINE

SUBROUTINE ReadShfRstFile(Core, FmInfo, cycleid, LpShf, nTRACERCntl, PE)
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(LpShf_Type) :: LpShf
TYPE(nTracerCntl_Type) :: nTRACERCntl
TYPE(PE_TYPE) :: PE

INTEGER :: cycleid
INTEGER :: iz0, ixa0, iya0, nxy0, nuse, ixyalist(10)
INTEGER :: nxa0, nya0, nxya0, nz0, nasywrite, nfxrwrite, nreadline
INTEGER :: iz, i, j, k

REAL :: Hm_mass0

LOGICAL :: luse, lRewind
!
Fxr => FmInfo%Fxr; Pin => Core%Pin
AsyType => Core%AsyInfo; Asy => Core%Asy
CoreIdx => Core%CoreIdx; CoreMap => Core%CoreMap
CellInfo => Core%CellInfo
nxc = Core%nxc; nyc = COre%nyc
IF(PE%MASTER) THEN
  WRITE(mesg, '(A, A40)') 'Reading Restart File : ', LpSHf%RstFiles(cycleid)
  CALL message(io8, TRUE, TRUE, mesg)
ENDIF

#ifdef MPI_ENV
CALL MPIWaitTurn(PE%MPI_Cmfd_Comm, PE%myCmfdRank, PE%nCmfdproc, .FALSE.)
#endif
CALL OPENFILE(io14, TRUE, TRUE, FALSE, LpSHf%RstFiles(cycleid))
IF(PE%nCMFDPROC .GT. 1) THEN
  PRINT '(5x, A, I5)', 'Reading Proc. No. ', PE%myCmfdRank
ENDIF
READ(io14)  nxa0, nya0, nxya0, nz0, Hm_mass0
LpShf%Hm_mass0(cycleid) = Hm_mass0
lRewind = .FALSE.
DO iz = PE%myzb, PE%myze
  IF(.NOT. Core%lFuelPlane(iz)) CYCLE
  iz0=iz
  IF(LpShf%lPln_MAP) THEN
    iz0 = LpShf%Pln_MAP(iz)
    IF(iz0 .LT. iz) lRewind = .TRUE.
    IF(lRewind) THEN
      REWIND(io14);
      READ(io14)  nxa0, nya0, nxya0, nz0, Hm_mass0
    ENDIF
  ENDIF
  !PRINT *, iz, iz0
  CALL LocPln(iz0)  !Locate File Handler at the iz-th plane
  READ(io14) iz0, nasywrite
  DO i = 1, nasywrite
    !Read Assembly Data Head
    READ(io14) ixa0, iya0, iz0, nxy0, nfxrwrite
    BACKSPACE(io14)
    luse = FALSE;nuse=0
    DO j = 1, Core%nxya
      IF(LpShf%Shf(j)%lFreshfuel) CYCLE
      IF(LpShf%Shf(j)%lNoneFuel) CYCLE
      IF(Lpshf%Shf(j)%cycleid .NE. cycleid) CYCLE
      IF(LpShf%Shf(j)%ix .EQ. ixa0 .AND. LpShf%Shf(j)%iy .EQ. iya0 ) THEN
        lUse = TRUE
        nuse = nuse + 1
        ixyalist(nuse) = j
      ENDIF
    ENDDO
    IF(.NOT. lUSE) THEN
      CALL SkipAsyData()
      CYCLE
    ENDIF
    DO j = 1, nuse
      IF(LpShf%Shf(ixyalist(j))%irot .EQ. 0) THEN
        CALL ReadAsyRstfile(ixyalist(j), iz, nreadline)
        LpShf%Shf(ixyalist(j))%lRead =.TRUE.
      ELSE
        CALL ReadAsyRstfile_ROT(ixyalist(j), iz, LpShf%Shf(ixyalist(j))%irot, nreadline)
        LpShf%Shf(ixyalist(j))%lRead =.TRUE.
      ENDIF
      IF(j + 1 .LE. nuse) THEN
        DO k = 1, nreadline
          BACKSPACE(io14)
        ENDDO
      ENDIF
    ENDDO
  ENDDO
  CALL RstReadingCheck(LpShf, Core%nxya, cycleid)
ENDDO
CLOSE(io14)
#ifdef MPI_ENV
CALL MPIWaitTurn(PE%MPI_Cmfd_Comm, PE%myCmfdRank, PE%nCmfdproc, .TRUE.)
#endif
NULLIFY(Fxr, Pin, AsyType, Asy)
NULLIFY(CoreIdx, CoreMap, CellInfo)
END SUBROUTINE

SUBROUTINE ReadAsyRstfile(ixya, iz, nreadline)
INTEGER :: ixya, iz, nreadline

INTEGER :: PinMapExt(100, 100)

INTEGER :: ixa0, iya0, iz0, nxy0, ifxr0, ilocalfxr0, nwrite
INTEGER :: ix, iy, nLocalFxr0, nDeplFxr, niso, niso_depl, idiso
INTEGER :: ixy, ipin, icel, iAsyType, i, j, k
INTEGER :: iy0, ix0

INTEGER :: FxrIdxSt, nLocalFxr, ifxr
REAL :: temp, burnup, pnum

LOGICAL :: lfuel, ldepl, lres, lgd, lh2o, ldump
LOGICAL :: lCentX, lCentY, lCentXY

iAsyType =  CoreMap(ixya)
lCentX = AsyType(iAsyType)%lCentX; lCentY = AsyType(iAsyType)%lCentY
lCentXY = AsyType(iAsyType)%lCentXY

iy0 = 0
CALL CP_CA(PinMapExt, 0, 100, 100)
DO iy = nyc-AsyType(iAsyType)%ny + 1, nyc
  iy0 = iy0 + 1; ix0 = 0
  DO ix = nxc-AsyType(iAsyType)%nx + 1, nxc
    ix0 = ix0 + 1
    ixy = AsyType(iAsyType)%Pin2DIdx(ix0, iy0)
    PinMapExt(ix, iy) = ixy
  ENDDO
ENDDO

nreadline = 0
READ(io14) ixa0, iya0, iz0, nxy0, nwrite
nreadline = nreadline + 1
DO ixy = 1, nxy0
  READ(io14) ix, iy, nLocalFxr0, nDeplFxr
  nreadline = nreadline + 1
  lDump = .FALSE.
  IF(PinMapExt(ix, iy) .EQ. 0) lDump = TRUE
  !Skip the file
  IF(lDump) THEN
    DO i = 1, nLocalFxr0
      READ(io14) ifxr0, ilocalfxr0, lfuel, ldepl, lres, lgd, lh2o
      nreadline = nreadline + 1
      IF(.NOT. lDepl) CYCLE
      READ(io14) niso, niso_depl, temp, burnup
      nreadline = nreadline + 1
      DO j = 1, niso_depl
        READ(io14) idiso, pnum
        nreadline = nreadline + 1
      ENDDO
    ENDDO
    CYCLE
  ENDIF

  ipin = Asy(ixya)%GlobalPinIdx(PinMapExt(ix, iy))
  icel = Pin(ipin)%Cell(iz)
  FxrIdxSt = Pin(ipin)%FxrIdxSt
  nLocalFxr= CellInfo(icel)%nFxr
  IF(nLocalFxr .NE. nLocalFxr0) THEN
    WRITE(mesg, '(A, I10)') "LP_SHF : Mismatch Cell Type :", ipin
    CALL TERMINATE(mesg)
  ENDIF

  DO i = 1, nLocalFxr0
    READ(io14) ifxr0, ilocalfxr0, lfuel, ldepl, lres, lgd, lh2o
    nreadline = nreadline + 1
    ifxr = FxrIdxSt + i - 1
    IF(Fxr(ifxr, iz)%ldepl .NEQV. lDepl) THEN
      WRITE(mesg, '(A, I10)') "LP_SHF : Mismatch Cell Type :", ipin
      CALL TERMINATE(mesg)
    ENDIF
    Fxr(ifxr, iz)%lfuel = lfuel; Fxr(ifxr, iz)%ldepl = ldepl
    Fxr(ifxr, iz)%lres = lres; Fxr(ifxr, iz)%lgd = lgd
    Fxr(ifxr, iz)%lh2o = lh2o
    IF(.NOT. lDepl) CYCLE
    READ(io14) niso, niso_depl, temp, burnup
    nreadline = nreadline + 1
    Fxr(ifxr, iz)%niso = niso; Fxr(ifxr, iz)%niso_depl = niso_depl
    Fxr(ifxr, iz)%burnup = burnup; !Fxr(ifxr, iz0)%temp = temp
    Fxr(ifxr, iz)%burnup_past = burnup;

    DO j = 1, niso_depl
      READ(io14) idiso, pnum
      nreadline = nreadline + 1
      Fxr(ifxr, iz)%idiso(j) = idiso; Fxr(ifxr, iz)%pnum(j) = pnum
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE ReadAsyRstfile_ROT(ixya, iz, irot, nreadline)
IMPLICIT NONE
INTEGER :: ixya, iz, irot, nreadline

INTEGER :: PinMapExt(100, 100), PinMapExt_ROT(100, 100)

INTEGER :: ixa0, iya0, iz0, nxy0, ifxr0, ilocalfxr0, nwrite
INTEGER :: ix, iy, ix0, iy0, nLocalFxr0, nDeplFxr, niso, niso_depl, idiso
INTEGER :: ixy, ipin, icel, iAsyType, i, j, k

REAL :: cx, cy, x, y, x0, y0, rotmat(2, 2)

INTEGER :: FxrIdxSt, nLocalFxr, ifxr
REAL :: temp, burnup, pnum

LOGICAL :: lfuel, ldepl, lres, lgd, lh2o, ldump
LOGICAL :: lCentX, lCentY, lCentXY

iAsyType =  CoreMap(ixya)
lCentX = AsyType(iAsyType)%lCentX; lCentY = AsyType(iAsyType)%lCentY
lCentXY = AsyType(iAsyType)%lCentXY

iy0 = 0
CALL CP_CA(PinMapExt, 0, 100, 100)
DO iy = nyc-AsyType(iAsyType)%ny + 1, nyc
  iy0 = iy0 + 1; ix0 = 0
  DO ix = nxc-AsyType(iAsyType)%nx + 1, nxc
    ix0 = ix0 + 1
    ixy = AsyType(iAsyType)%Pin2DIdx(ix0, iy0)
    PinMapExt(ix, iy) = ixy
  ENDDO
ENDDO

CALL CP_CA(PinMapExt_ROT, 0, 100, 100)

cx = 0.5_8 * real(nxc+1, 8); cy = 0.5_8 * real(nyc+1, 8)
DO iy = 1, nyc
  DO ix = 1, nxc
    x0 = real(ix, 8); y0 = real(nyc + 1) - real(iy, 8)
    x0 = x0 - cx; y0 = y0 - cy
    rotmat(1, 1) = dcos( -0.5 * pi * irot); rotmat(1, 2) = -dsin( -0.5 * pi * irot)
    rotmat(2, 1) = dsin( -0.5 * pi * irot); rotmat(2, 2) = dcos( -0.5 * pi * irot)
    x = rotmat(1, 1) * x0 + rotmat(1, 2) * y0
    y = rotmat(2, 1) * x0 + rotmat(2, 2) * y0
    x = x + cx; y = y + cy
    x = x; y = real(nyc+1) - y
    ix0 = nint(x); iy0 = nint(y)
    PinMapExt_ROT(ix0, iy0) = PinMapExt(ix, iy)
    CONTINUE
  ENDDO
ENDDO

nreadline = 0
READ(io14) ixa0, iya0, iz0, nxy0, nwrite
nreadline = nreadline + 1
DO ixy = 1, nxy0
  READ(io14) ix, iy, nLocalFxr0, nDeplFxr
  nreadline = nreadline + 1
  lDump = .FALSE.
  IF(PinMapExt_ROT(ix, iy) .EQ. 0) lDump = TRUE
  !Skip the file
  IF(lDump) THEN
    DO i = 1, nLocalFxr0
      READ(io14) ifxr0, ilocalfxr0, lfuel, ldepl, lres, lgd, lh2o
      nreadline = nreadline + 1
      IF(.NOT. lDepl) CYCLE
      READ(io14) niso, niso_depl, temp, burnup
      nreadline = nreadline + 1
      DO j = 1, niso_depl
        READ(io14) idiso, pnum
        nreadline = nreadline + 1
      ENDDO
    ENDDO
    CYCLE
  ENDIF

  ipin = Asy(ixya)%GlobalPinIdx(PinMapExt_ROT(ix, iy))
  icel = Pin(ipin)%Cell(iz)
  FxrIdxSt = Pin(ipin)%FxrIdxSt
  nLocalFxr= CellInfo(icel)%nFxr
  IF(nLocalFxr .NE. nLocalFxr0) THEN
    WRITE(mesg, '(A, I10)') "LP_SHF : Mismatch Cell Type :", ipin
    CALL TERMINATE(mesg)
  ENDIF

  DO i = 1, nLocalFxr0
    READ(io14) ifxr0, ilocalfxr0, lfuel, ldepl, lres, lgd, lh2o
    nreadline = nreadline + 1
    ifxr = FxrIdxSt + i - 1
    IF(Fxr(ifxr, iz)%ldepl .NEQV. lDepl) THEN
      WRITE(mesg, '(A, I10)') "LP_SHF : Mismatch Cell Type :", ipin
      CALL TERMINATE(mesg)
    ENDIF
    Fxr(ifxr, iz)%lfuel = lfuel; Fxr(ifxr, iz)%ldepl = ldepl
    Fxr(ifxr, iz)%lres = lres; Fxr(ifxr, iz)%lgd = lgd
    Fxr(ifxr, iz)%lh2o = lh2o
    IF(.NOT. lDepl) CYCLE
    READ(io14) niso, niso_depl, temp, burnup
    nreadline = nreadline + 1
    Fxr(ifxr, iz)%niso = niso; Fxr(ifxr, iz)%niso_depl = niso_depl
    Fxr(ifxr, iz)%burnup = burnup; !Fxr(ifxr, iz0)%temp = temp
    Fxr(ifxr, iz)%burnup_past = burnup;
    DO j = 1, niso_depl
      READ(io14) idiso, pnum
      nreadline = nreadline + 1
      Fxr(ifxr, iz)%idiso(j) = idiso; Fxr(ifxr, iz)%pnum(j) = pnum
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE LocPln(iz0)
INTEGER :: iz0
INTEGER :: iz, nasywrite
INTEGER :: i
DO
  READ(io14) iz, nasywrite
  IF(iz0 .EQ. iz) EXIT
  DO i = 1, nasywrite
    CALL SkipAsyData()
  ENDDO
ENDDO
BACKSPACE(io14)
END SUBROUTINE

SUBROUTINE SkipAsyData()
INTEGER :: ixa0, iya0, iz0, nxy0, ifxr0, ilocalfxr, nwrite
INTEGER :: ix, iy, nLocalFxr, nDeplFxr, niso, niso_depl, idiso
INTEGER :: ixy, i, j
REAL :: temp, burnup, pnum
LOGICAL :: lfuel, ldepl, lres, lgd, lh2o

READ(io14) ixa0, iya0, iz0, nxy0, nwrite
DO ixy = 1, nxy0
  READ(io14) ix, iy, nLocalFxr, nDeplFxr
  DO i = 1, nLocalFxr
    READ(io14) ifxr0, ilocalfxr, lfuel, ldepl, lres, lgd, lh2o
    IF(.NOT. lDepl) CYCLE
    READ(io14) niso, niso_depl, temp, burnup
    DO j = 1, niso_depl
      READ(io14) idiso, pnum
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE
SUBROUTINE RstReadingCheck(LpShf, nxya, cycleid)
TYPE(LpShf_Type) :: LpShf
INTEGER :: nxya
INTEGER :: ixya
INTEGER :: cycleid
DO ixya = 1, nxya
  IF(LpShf%Shf(ixya)%lFreshFuel) CYCLE
  IF(LpShf%Shf(ixya)%LNoneFuel) CYCLE
  IF(LpShf%Shf(ixya)%cycleid .NE. cycleid) CYCLE
  IF(.NOT. LpShf%Shf(ixya)%LRead) THEN
    WRITE(mesg, '(A, I5)')  'Can not find reloaing assembly :', ixya
    PRINT *, LpShf%Shf(ixya)%fields(1:15)
    CALL terminate(mesg)
  ENDIF
  LpShf%Shf(ixya)%LRead = .FALSE.
ENDDO
END SUBROUTINE
END MODULE
