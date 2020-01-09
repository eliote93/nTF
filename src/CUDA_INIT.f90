#include <CUDADEFINES.h>

#ifdef __PGI

MODULE CUDA_INIT

USE CUDA_MASTER
IMPLICIT NONE

CONTAINS

!--- Ray Processing Routines ----------------------------------------------------------------------

SUBROUTINE CUDAFastRay1DGen(Core, RayInfo)
USE TYPEDEF,		ONLY : CoreInfo_Type,		RayInfo_Type,		RotRayInfo_Type,		CoreRayInfo_Type,		&
						   AsyRayInfo_Type,		CellRayInfo_Type,	Asy_Type,				Pin_Type,				&
						   Cell_Type
USE PE_MOD,			ONLY : PE
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo

TYPE(RotRayInfo_Type), POINTER :: RotRay(:)
TYPE(CoreRayInfo_Type), POINTER :: CoreRay(:)
TYPE(AsyRayInfo_Type), POINTER :: AsyRay(:)
TYPE(CellRayInfo_Type), POINTER :: CellRay
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
INTEGER :: i, j, k, l, iz, iAx
INTEGER :: iRotRay, iCoreRay, iAsyRay, iCellRay, iAsy, iPin, iCell, iBaseCell, iRaySeg, iAzi, iDir
INTEGER :: nRotRay, nCoreRay, nAsyRay, nPinRay, nRaySeg
INTEGER :: nTotCoreRay, nTotPinRay, nTotRaySeg
INTEGER :: FsrIdxSt
INTEGER, ALLOCATABLE :: CoreRayBeg(:), PinRayBeg(:), RaySegBeg(:)

ALLOCATE(cuFastRay1D(cuGeometry%nAxType))

RotRay => RayInfo%RotRay
CoreRay => RayInfo%CoreRay
AsyRay => RayInfo%AsyRay

Asy => Core%Asy
Pin => Core%Pin
Cell => Core%CellInfo

nRotRay = RayInfo%nRotRay

ALLOCATE(CoreRayBeg(nRotRay), PinRayBeg(nRotRay), RaySegBeg(nRotRay))

DO iAx = 1, cuGeometry%nAxType

  DO iz = PE%myzb, PE%myze
    IF (cuGeometry%AxType(iz) .EQ. iAx) EXIT
  ENDDO

  ALLOCATE(cuFastRay1D(iAx)%RotRayIdxSt(nRotRay + 1))
  ALLOCATE(cuFastRay1D(iAx)%PhiAngInSvIdx(2, nRotRay))
  ALLOCATE(cuFastRay1D(iAx)%PhiAngOutSvIdx(2, nRotRay))

  cuFastRay1D(iAx)%RotRayIdxSt(1) = 1

  nTotCoreRay = 0; nTotPinRay = 0; nTotRaySeg = 0
  DO i = 1, nRotRay
    iRotRay = i
    CoreRayBeg(iRotRay) = nTotCoreRay + 1
    PinRayBeg(iRotRay) = nTotPinRay + 1
    RaySegBeg(iRotRay) = nTotRaySeg + 1
    nCoreRay = RotRay(iRotRay)%nRay
    nTotCoreRay = nTotCoreRay + nCoreRay
    cuFastRay1D(iAx)%RotRayIdxSt(iRotRay + 1) = nTotCoreRay + 1
    DO j = 1, nCoreRay
      iCoreRay = RotRay(iRotRay)%RayIdx(j)
      nAsyRay = CoreRay(iCoreRay)%nRay
      DO k = 1, nAsyRay
        iAsyRay = CoreRay(iCoreRay)%AsyRayIdx(k)
        iAsy = CoreRay(iCoreRay)%AsyIdx(k)
        IF (iAsy .EQ. 0) CYCLE
        nPinRay = AsyRay(iAsyRay)%nCellRay
	    nTotPinRay = nTotPinRay + nPinRay
        DO l = 1, nPinRay
          iPin = AsyRay(iAsyRay)%PinIdx(l)
          iCellRay = AsyRay(iAsyRay)%PinRayIdx(l)
		  iPin = Asy(iAsy)%GlobalPinIdx(ipin)
          iCell = Pin(ipin)%Cell(iz)
          iBaseCell = Cell(iCell)%BaseCellStr
		  nRaySeg = Cell(iBaseCell)%CellRay(iCellRay)%nSeg
		  nTotRaySeg = nTotRaySeg + nRaySeg
        ENDDO
      ENDDO
    ENDDO
    cuFastRay1D(iAx)%PhiAngInSvIdx(:, iRotRay) = RayInfo%PhiAngInSvIdx(iRotRay, :)
    cuFastRay1D(iAx)%PhiAngOutSvIdx(:, iRotRay) = RayInfo%PhiAngOutSvIdx(iRotRay, :)
  ENDDO

  ALLOCATE(cuFastRay1D(iAx)%CoreRayIdxSt(nTotCoreRay + 1))
  ALLOCATE(cuFastRay1D(iAx)%PinRayIdxSt(nTotPinRay + 1))
  ALLOCATE(cuFastRay1D(iAx)%iAzi(nTotCoreRay))
  ALLOCATE(cuFastRay1D(iAx)%iDir(nTotCoreRay))
  ALLOCATE(cuFastRay1D(iAx)%PinIdx(nTotPinRay))
  ALLOCATE(cuFastRay1D(iAx)%SurfIdx(2, nTotPinRay))
  ALLOCATE(cuFastRay1D(iAx)%FsrIdx(nTotRaySeg))
  ALLOCATE(cuFastRay1D(iAx)%LenSeg(nTotRaySeg))

  !$OMP PARALLEL PRIVATE(iRotRay, iCoreRay, iAsyRay, iCellRay, iAsy, iPin, iCell, iBaseCell, iRaySeg, iAzi, iDir,   &
  !$OMP                  CellRay, nCoreRay, nAsyRay, nPinRay, nRaySeg, nTotCoreRay, nTotPinRay, nTotRaySeg, FsrIdxSt)
  !$OMP DO SCHEDULE(GUIDED)
  DO i = 1, nRotRay
    iRotRay = i
    nTotCoreRay = CoreRayBeg(iRotRay)
    nTotPinRay = PinRayBeg(iRotRay)
    nTotRaySeg = RaySegBeg(iRotRay)
    nCoreRay = RotRay(iRotRay)%nRay
    DO j = 1, nCoreRay
      iCoreRay = RotRay(iRotRay)%RayIdx(j)
      nAsyRay = CoreRay(iCoreRay)%nRay
	  iAzi = CoreRay(iCoreRay)%iAng
  	  iDir = RotRay(iRotRay)%Dir(j)
	  cuFastRay1D(iAx)%CoreRayIdxSt(nTotCoreRay) = nTotPinRay
	  cuFastRay1D(iAx)%iAzi(nTotCoreRay) = iAzi
	  cuFastRay1D(iAx)%iDir(nTotCoreRay) = iDir
	  IF (iDir .EQ. 1) THEN
        DO k = 1, nAsyRay
          iAsyRay = CoreRay(iCoreRay)%AsyRayIdx(k)
          iAsy = CoreRay(iCoreRay)%AsyIdx(k)
          IF (iAsy .EQ. 0) CYCLE
          nPinRay = AsyRay(iAsyRay)%nCellRay
          DO l = 1, nPinRay
            iPin = AsyRay(iAsyRay)%PinIdx(l)
  		    iPin = Asy(iAsy)%GlobalPinIdx(ipin)
		    FsrIdxSt = Pin(iPin)%FsrIdxSt
            iCell = Pin(ipin)%Cell(iz)
		    iCellRay = AsyRay(iAsyRay)%PinRayIdx(l)
            iBaseCell = Cell(iCell)%BaseCellStr
		    CellRay => Cell(iBaseCell)%CellRay(iCellRay)
		    cuFastRay1D(iAx)%PinIdx(nTotPinRay) = iPin
		    cuFastRay1D(iAx)%SurfIdx(1, nTotPinRay) = AsyRay(iAsyRay)%PinRaySurf(1, l)   !--- Incoming
		    cuFastRay1D(iAx)%SurfIdx(2, nTotPinRay) = AsyRay(iAsyRay)%PinRaySurf(2, l)   !--- Outgoing
		    cuFastRay1D(iAx)%PinRayIdxSt(nTotPinRay) = nTotRaySeg
		    nRaySeg = CellRay%nSeg
		    DO iRaySeg = 1, nRaySeg
		      cuFastRay1D(iAx)%FsrIdx(nTotRaySeg) = FsrIdxSt + CellRay%LocalFsrIdx(iRaySeg) - 1
		      cuFastRay1D(iAx)%LenSeg(nTotRaySeg) = CellRay%LenSeg(iRaySeg) * 0.001
		      nTotRaySeg = nTotRaySeg + 1
		    ENDDO
		    nTotPinRay = nTotPinRay + 1
          ENDDO
        ENDDO
	  ELSE
	    DO k = nAsyRay, 1, -1
          iAsyRay = CoreRay(iCoreRay)%AsyRayIdx(k)
          iAsy = CoreRay(iCoreRay)%AsyIdx(k)
          IF (iAsy .EQ. 0) CYCLE
          nPinRay = AsyRay(iAsyRay)%nCellRay
          DO l = nPinRay, 1, -1
            iPin = AsyRay(iAsyRay)%PinIdx(l)
  		    iPin = Asy(iAsy)%GlobalPinIdx(ipin)
		    FsrIdxSt = Pin(iPin)%FsrIdxSt
            iCell = Pin(ipin)%Cell(iz)
		    iCellRay = AsyRay(iAsyRay)%PinRayIdx(l)
            iBaseCell = Cell(iCell)%BaseCellStr
		    CellRay => Cell(iBaseCell)%CellRay(iCellRay)
		    cuFastRay1D(iAx)%PinIdx(nTotPinRay) = iPin
		    cuFastRay1D(iAx)%SurfIdx(1, nTotPinRay) = AsyRay(iAsyRay)%PinRaySurf(2, l)   !--- Incoming
		    cuFastRay1D(iAx)%SurfIdx(2, nTotPinRay) = AsyRay(iAsyRay)%PinRaySurf(1, l)   !--- Outgoing
		    cuFastRay1D(iAx)%PinRayIdxSt(nTotPinRay) = nTotRaySeg
		    nRaySeg = CellRay%nSeg
		    DO iRaySeg = nRaySeg, 1, -1
		      cuFastRay1D(iAx)%FsrIdx(nTotRaySeg) = FsrIdxSt + CellRay%LocalFsrIdx(iRaySeg) - 1
		      cuFastRay1D(iAx)%LenSeg(nTotRaySeg) = CellRay%LenSeg(iRaySeg) * 0.001
		      nTotRaySeg = nTotRaySeg + 1
		    ENDDO
		    nTotPinRay = nTotPinRay + 1
          ENDDO
        ENDDO
	  ENDIF
      nTotCoreRay = nTotCoreRay + 1
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL

  cuFastRay1D(iAx)%CoreRayIdxSt(nTotCoreRay + 1) = nTotPinRay + 1
  cuFastRay1D(iAx)%PinRayIdxSt(nTotPinRay + 1) = nTotRaySeg + 1

ENDDO

DEALLOCATE(CoreRayBeg, PinRayBeg, RaySegBeg)

END SUBROUTINE

!--- Geometry Processing Routines -----------------------------------------------------------------

SUBROUTINE AxialCheck(Core)
USE TYPEDEF,		ONLY : CoreInfo_Type,		Pin_Type,			Cell_Type
USE PE_MOD,			ONLY : PE
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
LOGICAL, ALLOCATABLE :: Visited(:)
LOGICAL, ALLOCATABLE :: RelationTable(:, :)
INTEGER, ALLOCATABLE :: CellArrangement(:, :)
INTEGER :: myzb, myze, nxy, nz, nAxType
INTEGER :: i, iAx, iz, iz1, iz2, icel, ibcel

Pin => Core%Pin
Cell => Core%CellInfo
myzb = PE%myzb
myze = PE%myze
nxy = Core%nxy
nz = PE%myze - PE%myzb + 1

ALLOCATE(Visited(myzb : myze))
ALLOCATE(RelationTable(myzb : myze, myzb : myze))
ALLOCATE(CellArrangement(nxy, myzb : myze))

DO iz = myzb, myze
  DO i = 1, nxy
	icel = Pin(i)%Cell(iz)
    ibcel = Cell(icel)%BaseCellStr
	CellArrangement(i, iz) = ibcel
  ENDDO
ENDDO

DO iz1 = myzb, myze
  DO iz2 = iz1 + 1, myze
    IF (ANY(CellArrangement(:, iz1) .NE. CellArrangement(:, iz2))) THEN
      RelationTable(iz1, iz2) = .FALSE.
    ELSE
	  RelationTable(iz1, iz2) = .TRUE.
	ENDIF
  ENDDO
ENDDO

ALLOCATE(cuGeometry%AxType(myzb : myze))

Visited = .FALSE.
nAxType = 1
DO iz1 = myzb, myze
  IF (Visited(iz1)) CYCLE
  Visited(iz1) = .TRUE.
  cuGeometry%AxType(iz1) = nAxType
  DO iz2 = iz1 + 1, myze
    IF (RelationTable(iz1, iz2)) THEN
	  Visited(iz2) = .TRUE.
	  cuGeometry%AxType(iz2) = nAxType
	ENDIF
  ENDDO
  IF (ALL(Visited)) EXIT
  nAxType = nAxType + 1
ENDDO

cuGeometry%nAxType = nAxType

DEALLOCATE(Visited)
DEALLOCATE(RelationTable)
DEALLOCATE(CellArrangement)

END SUBROUTINE

SUBROUTINE SetGeometry(Core, RayInfo, FmInfo, GroupInfo, cuDevice)
USE TYPEDEF,        ONLY : CoreInfo_Type,       RayInfo_Type,       FmInfo_Type,        GroupInfo_Type,             &
                           Pin_Type,            Cell_Type
USE CNTL,           ONLY : nTracerCntl
USE PE_MOD,         ONLY : PE
USE CMFD_COMMON,    ONLY : SetSuperPin
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(cuDevice_Type) :: cuDevice

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
INTEGER :: i, j, iz, izf, iAx, ifxr, ifsr, ireg, icel, ibcel, ipin, ipin_map, ipintype
INTEGER :: nfxr, nfsr, nxy, nRotRay, nAziAngle, nLocalFxr, nLocalFsr, nFsrInFxr
INTEGER :: myzb, myze
INTEGER :: FsrIdxSt, FxrIdxSt

Pin => Core%Pin
Cell => Core%CellInfo
nfxr = Core%nCoreFxr
nfsr = Core%nCoreFsr
nxy = Core%nxy
nRotRay = RayInfo%nRotRay
nAziAngle = RayInfo%nAziAngle
myzb = PE%myzb
myze = PE%myze

cuGeometry%nfxr = nfxr
cuGeometry%nfsr = nfsr
cuGeometry%nxy = nxy
cuGeometry%nxyc = nxy
cuGeometry%nxya = Core%nxya
cuGeometry%ng = GroupInfo%ng
cuGeometry%nx = Core%nx
cuGeometry%ny = Core%ny
cuGeometry%nRotRay = nRotRay
cuGeometry%nModRay = RayInfo%nModRay
cuGeometry%nAziAngle = nAziAngle
cuGeometry%nPolarAngle = RayInfo%nPolarAngle
cuGeometry%nPhiAngSv = RayInfo%nPhiAngSv
cuGeometry%myzb = myzb
cuGeometry%myze = myze
cuGeometry%AxBC = Core%AxBC
cuGeometry%RadBC = Core%RadBC(1 : 4)
cuGeometry%l3dim = nTracerCntl%l3dim

cuCntl%lSuperpin = cuCntl%lSuperpin .AND. Core%lGap

CALL SetPlaneMap(Core, cuDevice)
CALL SetSuperPin(Core, cuGeometry%superPin, cuGeometry%nxyc, cuGeometry%myzb, cuGeometry%myze, cuCntl%lSuperpin)
CALL SetPinMap(Core, FmInfo, cuDevice)

ALLOCATE(cuGeometry%FsrVol(nfsr, cuGeometry%nAxType))

DO iz = myzb, myze
  iAx = cuGeometry%AxType(iz)
  DO ipin = 1, nxy
    ipintype = Core%Pin(ipin)%PinType
    icel = Core%PinInfo(ipintype)%Cell(iz)
    ibcel = Core%CellInfo(icel)%basecellstr
    nLocalFsr = Core%CellInfo(ibcel)%nFsr
    FsrIdxSt = Core%Pin(ipin)%FsrIdxSt
    DO ireg = 1, nLocalFsr
      ifsr = ireg + FsrIdxSt - 1
      cuGeometry%FsrVol(ifsr, iAx) = Core%CellInfo(ibcel)%vol(ireg)
    ENDDO
  ENDDO
ENDDO

ALLOCATE(cuGeometry%Fsr2Fxr(nfsr, myzb : myze)); cuGeometry%Fsr2Fxr = 0

DO iz = myzb, myze
  DO ipin = 1, nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
    icel = Pin(ipin)%Cell(iz); nLocalFxr = Cell(icel)%nFxr
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j - 1
      nFsrInFxr = Cell(icel)%nFsrInFxr(j)
      DO i = 1, nFsrInFxr
        ifsr = FsrIdxSt + Cell(icel)%MapFxr2FsrIdx(i, j) - 1
        cuGeometry%Fsr2Fxr(ifsr, iz) = ifxr
      ENDDO
    ENDDO
  ENDDO
ENDDO

cuDevice%lFuel = .FALSE.
DO iz = cuDevice%myzb, cuDevice%myze
  IF (Core%lFuelPlane(iz)) cuDevice%lFuel = .TRUE.
ENDDO

ALLOCATE(cuGeometry%RotRayCount(nAziAngle / 2))
ALLOCATE(cuGeometry%RotRayList(nRotRay, nAziAngle / 2))
cuGeometry%RotRayCount = RayInfo%RotRayAziList(0, :)
cuGeometry%RotRayList = RayInfo%RotRayAziList(1 : nRotRay, :)

END SUBROUTINE

SUBROUTINE SetPlaneMap(Core, cuDevice)
USE TYPEDEF,        ONLY : CoreInfo_Type
USE PE_MOD,         ONLY : PE
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(cuDevice_Type) :: cuDevice

REAL, POINTER :: hz(:)
INTEGER :: nz, nzfm, nzCMFD
INTEGER :: iz, izf
INTEGER :: myzb, myze
INTEGER, POINTER :: nSubplane(:)

hz => Core%hz
nz = Core%nz
myzb = cuGeometry%myzb
myze = cuGeometry%myze

ALLOCATE(nSubplane(nz))

DO iz = 1, nz
  nSubplane(iz) = INT(hz(iz) / cuCntl%CMFDHeight) + 1
ENDDO

IF (.NOT. cuCntl%lSubplane) nSubplane = 1

nzCMFD = 0
DO iz = myzb, myze
  nzCMFD = nzCMFD + nSubplane(iz)
ENDDO

cuGeometry%nzCMFD = nzCMFD

ALLOCATE(cuGeometry%planeMap(nzCMFD))
ALLOCATE(cuGeometry%fmRange(myzb : myze, 2))

izf = 1
DO iz = myzb, myze
  nzfm = nSubplane(iz)
  cuGeometry%planeMap(izf : izf + nzfm - 1) = iz
  izf = izf + nzfm
  cuGeometry%fmRange(iz, 1) = izf - nzfm
  cuGeometry%fmRange(iz, 2) = izf - 1
ENDDO

ALLOCATE(cuGeometry%hz(myzb : myze))
ALLOCATE(cuGeometry%hzfm(0 : nzCMFD + 1))

cuGeometry%hz(myzb : myze) = hz(myzb : myze)
DO izf = 1, nzCMFD
  iz = cuGeometry%planeMap(izf)
  cuGeometry%hzfm(izf) = cuGeometry%hz(iz) / nSubplane(iz)
ENDDO

IF (myzb .EQ. 1) THEN
  cuGeometry%hzfm(0) = cuGeometry%hzfm(1)
ELSE
  cuGeometry%hzfm(0) = hz(myzb - 1) / nSubplane(myzb - 1)
ENDIF
IF (myze .EQ. nz) THEN
  cUGeometry%hzfm(nzCMFD + 1) = cuGeometry%hzfm(nzCMFD)
ELSE
  cuGeometry%hzfm(nzCMFD + 1) = hz(myze + 1) / nSubplane(myze + 1)
ENDIF

cuDevice%myzbf = cuGeometry%fmRange(cuDevice%myzb, 1)
cuDevice%myzef = cuGeometry%fmRange(cuDevice%myze, 2)
cuDevice%nzCMFD = cuDevice%myzef - cuDevice%myzbf + 1

DEALLOCATE(nSubplane)

END SUBROUTINE

SUBROUTINE SetPinMap(Core, FmInfo, cuDevice)
USE TYPEDEF,        ONLY : CoreInfo_Type,       FmInfo_Type,        Pin_Type,           Cell_Type,                  &
                           FxrInfo_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(cuDevice_Type) :: cuDevice

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(superPin_Type), POINTER :: superPin(:)
INTEGER :: node(cuGeometry%nx, cuGeometry%ny)
INTEGER :: ng, nx, ny, nxy, nzCMFD
INTEGER :: myzb, myze, myzbf, myzef
INTEGER :: i, j, ix, iy, iz, ixy, ixy_map, izf, icel, ipin, ipin_red, ipin_black, ipin_map, ifxr, ierr
INTEGER :: color(2), mp(2) = (/ 2, 1 /)
INTEGER, POINTER :: colorMap(:, :, :), nz(:), ixRange(:, :)

Fxr => FmInfo%Fxr
Pin => Core%Pin
Cell => Core%CellInfo
superPin => cuGeometry%superPin
ng = cuGeometry%ng
nx = cuGeometry%nx
ny = cuGeometry%ny
nxy = cuGeometry%nxyc
nzCMFD = cuGeometry%nzCMFD
myzb = cuGeometry%myzb
myze = cuGeometry%myze

ALLOCATE(cuGeometry%pinMap(nxy))
ALLOCATE(cuGeometry%pinMapRev(-1 : nxy))

!--- Global Pin Map

ALLOCATE(ixRange(nx, ny))

node = 0
DO ipin = 1, nxy
  ix = superPin(ipin)%ix
  iy = superPin(ipin)%iy
  node(ix, iy) = ipin
ENDDO

DO iy = 1, ny
  DO ix = 1, nx
    IF (node(ix, iy) .GT. 0) EXIT
  ENDDO
  ixRange(1, iy) = ix
  DO ix = nx, 1, -1
    IF (node(ix, iy) .GT. 0) EXIT
  ENDDO
  ixRange(2, iy) = ix
ENDDO

ipin = 0
DO iy = 1, ny
  DO ix = ixRange(1, iy), ixRange(2, iy)
    ipin = ipin + 1
    cuGeometry%pinMap(ipin) = node(ix, iy)
    cuGeometry%pinMapRev(node(ix, iy)) = ipin
  ENDDO
ENDDO
cuGeometry%pinMapRev(-1) = -1
cuGeometry%pinMapRev(0) = 0

!--- Red-Black Pin Map

ALLOCATE(colorMap(nx, ny, nzCMFD))
ALLOCATE(nz(NUM_CUDA_PROC))

CALL MPI_ALLGATHER(nzCMFD, 1, MPI_INT, nz, 1, MPI_INT, MPI_CUDA_COMM, ierr)

DO izf = 1, nzCMFD
  IF (MPI_CUDA_RANK .EQ. 0) THEN
    iz = izf
  ELSE
    iz = izf + sum(nz(1 : MPI_CUDA_RANK))
  ENDIF
  IF (mod(iz, 2) .EQ. 0) color = (/ red, black /)
  IF (mod(iz, 2) .EQ. 1) color = (/ black, red /)
  DO iy = 1, ny
    DO ix = ixRange(1, iy), ixRange(2, iy)
      IF (mod(ix, 2) .EQ. 0) CYCLE
      colorMap(ix, iy, izf) = color(1)
    ENDDO
    DO ix = ixRange(1, iy), ixRange(2, iy)
      IF (mod(ix, 2) .EQ. 1) CYCLE
      colorMap(ix, iy, izf) = color(2)
    ENDDO
    color(1) = mp(color(1))
    color(2) = mp(color(2))
  ENDDO
ENDDO

myzbf = cuDevice%myzbf; myzef = cuDevice%myzef
ALLOCATE(cuDevice%planeMapRB(nxy * cuDevice%nzCMFD))
ALLOCATE(cuDevice%pinMapRB(nxy * cuDevice%nzCMFD))
ALLOCATE(cuDevice%pinMapRevRB(-1 : nxy, myzbf : myzef))
cuDevice%pinMapRevRB(-1, :) = -1
cuDevice%pinMapRevRB(0, :) = 0
ipin = 0; cuDevice%rbBeg(red) = ipin + 1
DO izf = myzbf, myzef
  cuDevice%rbRange(1, izf, red) = ipin + 1
  cuDevice%nxyRB(izf, red) = ipin
  DO iy = 1, ny
    DO ix = ixRange(1, iy), ixRange(2, iy)
      IF (colorMap(ix, iy, izf) .NE. red) CYCLE
      ipin = ipin + 1
      cuDevice%pinMapRB(ipin) = node(ix, iy)
      cuDevice%pinMapRevRB(node(ix, iy), izf) = ipin
      cuDevice%planeMapRB(ipin) = izf
    ENDDO
  ENDDO
  cuDevice%rbRange(2, izf, red) = ipin
  cuDevice%nxyRB(izf, red) = ipin - cuDevice%nxyRB(izf, red)
ENDDO
cuDevice%rbEnd(red) = ipin
cuDevice%rbBeg(black) = ipin + 1
DO izf = myzbf, myzef
  cuDevice%rbRange(1, izf, black) = ipin + 1
  cuDevice%nxyRB(izf, black) = ipin
  DO iy = 1, ny
    DO ix = ixRange(1, iy), ixRange(2, iy)
      IF (colorMap(ix, iy, izf) .NE. black) CYCLE
      ipin = ipin + 1
      cuDevice%pinMapRB(ipin) = node(ix, iy)
      cuDevice%pinMapRevRB(node(ix, iy), izf) = ipin
      cuDevice%planeMapRB(ipin) = izf
    ENDDO
  ENDDO
  cuDevice%rbRange(2, izf, black) = ipin
  cuDevice%nxyRB(izf, black) = ipin - cuDevice%nxyRB(izf, black)
ENDDO
cuDevice%rbEnd(black) = ipin
cuDevice%nxyzRB(red) = cuDevice%rbEnd(red) - cuDevice%rbBeg(red) + 1
cuDevice%nxyzRB(black) = cuDevice%rbEnd(black) - cuDevice%rbBeg(black) + 1

!--- Natural MPI Map
cuDevice%bottomRange = (/ 1, nxy /)
cuDevice%topRange = (/ nxy * (cuDevice%nzCMFD - 1) + 1, nxy * cuDevice%nzCMFD /)
!--- Red-black MPI Map
cuDevice%bottomRangeRB = cuDevice%rbRange(:, cuDevice%myzbf, :)
cuDevice%topRangeRB = cuDevice%rbRange(:, cuDevice%myzef, :)

ALLOCATE(cuGeometry%PinVolFm(nxy, nzCMFD))

DO izf = 1, nzCMFD
  iz = cuGeometry%planeMap(izf)
  DO ipin = 1, nxy
    cuGeometry%PinVolFm(ipin, izf) = superPin(ipin)%Area * cuGeometry%hzfm(izf)
  ENDDO
ENDDO

ALLOCATE(cuGeometry%lRefPin(nxy))
ALLOCATE(cuGeometry%lRefCell(myzb : myze, nxy))

DO ipin = 1, nxy
  ipin_map = cuGeometry%pinMap(ipin)
  cuGeometry%lRefPin(ipin) = .NOT. ANY(superPin(ipin_map)%lFuel)
  DO iz = myzb, myze
    cuGeometry%lRefCell(iz, ipin) = .NOT. superPin(ipin_map)%lFuel(iz)
  ENDDO
ENDDO

ALLOCATE(cuGeometry%lH2OCell(myzb : myze, nxy))

DO ixy = 1, nxy
  ixy_map = cuGeometry%pinMap(ixy)
  DO iz = myzb, myze
    cuGeometry%lH2OCell(iz, ixy) = .TRUE.
    DO i = 1, superPin(ixy_map)%nxy
      ipin = superPin(ixy_map)%pin(i)
      icel = Pin(ipin)%Cell(iz)
      DO j = 1, Cell(icel)%nFxr
        ifxr = Pin(ipin)%FxrIdxSt + j - 1
        IF (.NOT. Fxr(ifxr, iz)%lH2O) cuGeometry%lH2OCell(iz, ixy) = .FALSE.
      ENDDO
    ENDDO
  ENDDO
ENDDO

DEALLOCATE(ixRange)
DEALLOCATE(colorMap)
DEALLOCATE(nz)

END SUBROUTINE

SUBROUTINE AssignMOCPlane(Core, cuDevice)
USE TYPEDEF,		ONLY : CoreInfo_Type
USE PE_MOD,			ONLY : PE
USE CNTL,			ONLY : nTracerCntl
USE OPENACC
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(cuDevice_Type) :: cuDevice

INTEGER :: iz, myzb, myze

myzb = PE%myzb
myze = PE%myze

cuDevice%myzb = myzb
cuDevice%myze = myze
cuDevice%nz = myze - myzb + 1
IF (myzb .EQ. 1) cuDevice%lBottom = .TRUE.
IF (myze .EQ. Core%nz) cuDevice%lTop = .TRUE.

DO iz = myzb + 1, myze
  IF (cuGeometry%AxType(myzb) .NE. cuGeometry%AxType(iz)) THEN
    cuDevice%lRayStatic = .FALSE.; EXIT
  ENDIF
ENDDO
cuDevice%myAxType = cuGeometry%AxType(myzb)

END SUBROUTINE

SUBROUTINE SetMPIEnv()
USE PE_MOD,         ONLY : PE
IMPLICIT NONE

INTEGER :: color, status(MPI_STATUS_SIZE)

color = 0; IF (PE%lCMFDGrp) color = 1

CALL MPI_COMM_SPLIT(PE%MPI_NTRACER_COMM, color, PE%myCMFDRank, MPI_CUDA_COMM, status)
CALL MPI_COMM_RANK(MPI_CUDA_COMM, MPI_CUDA_RANK, status)
CALL MPI_COMM_SIZE(MPI_CUDA_COMM, NUM_CUDA_PROC, status)

CALL MPI_COMM_SPLIT_TYPE(MPI_CUDA_COMM, MPI_COMM_TYPE_SHARED, MPI_CUDA_RANK, MPI_INFO_NULL,                         &
                         MPI_CUDA_SHARED_COMM, status)
CALL MPI_COMM_RANK(MPI_CUDA_SHARED_COMM, MPI_CUDA_SHARED_RANK, status)
CALL MPI_COMM_SIZE(MPI_CUDA_SHARED_COMM, NUM_CUDA_SHARED_PROC, status)

END SUBROUTINE

!--- Data Allocation Routines ---------------------------------------------------------------------

SUBROUTINE SetControlVar(Core, RayInfo, cuDevice)
USE TYPEDEF,		ONLY : CoreInfo_Type,		RayInfo_Type
USE CNTL,			ONLY : nTracerCntl
USE PE_MOD,			ONLY : PE
USE CUDAFOR
USE OPENACC
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo
TYPE(cuDevice_Type) :: cuDevice

TYPE(cudaDeviceProp) :: cuProperty
INTEGER :: PinBeg, PinEnd, PinCount
INTEGER :: i, ierr

ierr = cudaGetDeviceProperties(cuProperty, MPI_CUDA_SHARED_RANK)
cuDevice%cuSMXCount = cuProperty%multiProcessorCount
cuDevice%cuArchitecture = cuProperty%major
cuDevice%cuWarpSize = cuProperty%warpSize
cuDevice%cuMaxThreadPerSMX = cuProperty%maxThreadsPerMultiprocessor
cuDevice%cuMaxThreadPerBlock = cuProperty%maxThreadsPerBlock
cuDevice%cuMaxWarpPerSMX = cuProperty%maxThreadsPerMultiprocessor / cuProperty%warpSize

SELECT CASE (cuDevice%cuArchitecture)
CASE (2)   !--- Fermi
  cuDevice%cuMaxBlockPerSMX = 8
CASE (3)   !--- Kepler
  cuDevice%cuMaxBlockPerSMX = 16
CASE (5)   !--- Maxwell
  cuDevice%cuMaxBlockPerSMX = 32
CASE (6)   !--- Pascal
  cuDevice%cuMaxBlockPerSMX = 32
CASE (7)   !--- Volta, Turing
  cuDevice%cuMaxBlockPerSMX = 32
END SELECT

cuDevice%cuWarpPerBlock = cuDevice%cuMaxWarpPerSMX / cuDevice%cuMaxBlockPerSMX
cuDevice%cuThreadPerBlock = cuDevice%cuWarpPerBlock * cuDevice%cuWarpSize

IF (cuDevice%lFullWarp) THEN
  cuDevice%sharedMemoryDim = cuDevice%cuThreadPerBlock
ELSE
  cuDevice%sharedMemoryDim = 2 * cuGeometry%ng
ENDIF

cuDevice%RotRayBeg = 1
cuDevice%RotRayEnd = RayInfo%nRotRay
cuDevice%nRotRay = RayInfo%nRotRay

IF (NUM_CUDA_PROC .NE. 1) cuCntl%lMulti = .TRUE.

END SUBROUTINE

SUBROUTINE SetConstVar(Core, RayInfo, GroupInfo)
USE TYPEDEF,		ONLY : CoreInfo_Type,		RayInfo_Type,		GroupInfo_Type,									&
						   AziAngleInfo_Type,	PolarAngle_Type
USE CUDA_CONST
USE CNTL,			ONLY : nTracerCntl
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo
TYPE(GroupInfo_Type) :: GroupInfo

TYPE(AziAngleInfo_Type), POINTER :: AziAngle(:)
TYPE(PolarAngle_Type), POINTER :: PolarAngle(:)
INTEGER :: ig, iazi, ipol
REAL :: wttemp, wtsin2, wtcos, wtpolar
REAL :: wt_host(GPU_MAX_POLAR, GPU_MAX_AZI), Comp_host(GPU_MAX_ORDER, GPU_MAX_POLAR, GPU_MAX_AZI)

AziAngle => RayInfo%AziAngle
PolarAngle => RayInfo%PolarAngle
wt_host = 0; Comp_host = 0

nofg = GroupInfo%nofg
norg = GroupInfo%norg
nchi = GroupInfo%nchi
nRotRay = RayInfo%nRotRay
nModRay = RayInfo%nModRay
nAziAngle = RayInfo%nAziAngle
nPolarAngle = RayInfo%nPolarAngle
nPhiAngSv = RayInfo%nPhiAngSv
nFsr = Core%nCoreFsr
nxy = cuGeometry%nxy
nxyc = cuGeometry%nxyc
ng = GroupInfo%ng
dir = reshape((/ 1, 2, 2, 1 /), shape(dir))
inc = (/ 1, -1 /)

ScatOd = nTracerCntl%ScatOd
IF (nTracerCntl%ScatOd .EQ. 1) nMoment = 2
IF (nTracerCntl%ScatOd .EQ. 2) nMoment = 5
IF (nTracerCntl%ScatOd .EQ. 3) nMoment = 9

DO ig = 1, GroupInfo%ng
  InScatRange(:, ig) = GroupInfo%InScatRange(:, ig)
ENDDO

DO iazi = 1, RayInfo%nAziAngle / 2
  AziMap(iazi, 1) = 1
  AziMap(iazi, 2) = 2
  AziMap(RayInfo%nAziAngle - iazi + 1, 1) = 3
  AziMap(RayInfo%nAziAngle - iazi + 1, 2) = 4
ENDDO

DO ipol = 1, RayInfo%nPolarAngle
  rsinv(ipol) = 1.0 / PolarAngle(ipol)%sinv
  sinv(ipol) = PolarAngle(ipol)%sinv
ENDDO

DO iazi = 1, RayInfo%nAziAngle
  wttemp = AziAngle(iazi)%weight * AziAngle(iazi)%del
  DO ipol = 1, RayInfo%nPolarAngle
	wt_host(ipol, iazi) = wttemp * PolarAngle(ipol)%weight * PolarAngle(ipol)%sinv
  ENDDO
ENDDO

DO iazi = 1, RayInfo%nAziAngle
  wttemp = AziAngle(iazi)%weight * AziAngle(iazi)%del
  DO ipol = 1, RayInfo%nPolarAngle
    wtsurf(1, ipol, iazi) = wttemp * PolarAngle(ipol)%weight / abs(AziAngle(iazi)%sinv)
    wtsurf(3, ipol, iazi) = wttemp * PolarAngle(ipol)%weight / abs(AziAngle(iazi)%sinv)
    wtsurf(2, ipol, iazi) = wttemp * PolarAngle(ipol)%weight / abs(AziAngle(iazi)%cosv)
    wtsurf(4, ipol, iazi) = wttemp * PolarAngle(ipol)%weight / abs(AziAngle(iazi)%cosv)
  ENDDO
ENDDO

IF (nTracerCntl%lScat1) THEN
  DO ipol = 1, RayInfo%nPolarAngle
    wttemp = PolarAngle(ipol)%sinv
    DO iazi = 1, RayInfo%nAziAngle
	  Comp_host(1, ipol, iazi) = wttemp * AziAngle(iazi)%cosv
      Comp_host(2, ipol, iazi) = wttemp * AziAngle(iazi)%sinv
      mwt(1, ipol, iazi) = Comp_host(1, ipol, iazi) * wt_host(ipol, iazi)
      mwt(2, ipol, iazi) = Comp_host(2, ipol, iazi) * wt_host(ipol, iazi)
    ENDDO
  ENDDO
  IF (nTracerCntl%ScatOd .GE. 2) THEN
    DO ipol = 1, RayInfo%nPolarAngle
      wttemp = PolarAngle(ipol)%sinv
      wtsin2 = PolarAngle(ipol)%sinv * PolarAngle(ipol)%sinv
      wtcos = PolarAngle(ipol)%cosv
      wtpolar =  1.5 * PolarAngle(ipol)%cosv * PolarAngle(ipol)%cosv - 0.5
      DO iazi = 1, RayInfo%nAziAngle
        Comp_host(3, ipol, iazi) = wtpolar
        Comp_host(4, ipol, iazi) = wtsin2 * (1.0 - 2.0 * AziAngle(iazi)%sinv * AziAngle(iazi)%sinv)
        Comp_host(5, ipol, iazi) = wtsin2 * (2.0 * AziAngle(iazi)%sinv * AziAngle(iazi)%cosv)
        mwt(3, ipol, iazi) = Comp_host(3, ipol, iazi) * wt_host(ipol, iazi)
        mwt(4, ipol, iazi) = 0.75 * Comp_host(4, ipol, iazi) * wt_host(ipol, iazi)
        mwt(5, ipol, iazi) = 0.75 * Comp_host(5, ipol, iazi) * wt_host(ipol, iazi)
      ENDDO
    ENDDO
  ENDIF
  IF (nTracerCntl%ScatOd .EQ. 3) THEN
    DO ipol = 1, RayInfo%nPolarAngle
      wttemp = PolarAngle(ipol)%sinv
      DO iazi = 1, RayInfo%nAziAngle
        Comp_host(6, ipol, iazi) = (5.0 * PolarAngle(ipol)%cosv * PolarAngle(ipol)%cosv - 1.0) * wttemp * AziAngle(iazi)%cosv
        Comp_host(7, ipol, iazi) = (5.0 * PolarAngle(ipol)%cosv * PolarAngle(ipol)%cosv - 1.0) * wttemp * AziAngle(iazi)%sinv
        Comp_host(8, ipol, iazi) = (wttemp ** 3) * (4.0 * (AziAngle(iazi)%cosv ** 3) - 3.0 * AziAngle(iazi)%cosv)
        Comp_host(9, ipol, iazi) = (wttemp ** 3) * (-4.0 * (AziAngle(iazi)%sinv ** 3) + 3.0 * AziAngle(iazi)%sinv)
        mwt(6, ipol, iazi) = 0.375 * Comp_host(6, ipol, iazi) * wt_host(ipol, iazi)
        mwt(7, ipol, iazi) = 0.375 * Comp_host(7, ipol, iazi) * wt_host(ipol, iazi)
        mwt(8, ipol, iazi) = 0.625 * Comp_host(8, ipol, iazi) * wt_host(ipol, iazi)
        mwt(9, ipol, iazi) = 0.625 * Comp_host(9, ipol, iazi) * wt_host(ipol, iazi)
      ENDDO
    ENDDO
  ENDIF
ENDIF

wt = wt_host
Comp = Comp_host

END SUBROUTINE

SUBROUTINE SetGcInfo(GroupInfo, GcGroupInfo)
USE TYPEDEF,        ONLY : GroupInfo_Type
IMPLICIT NONE

TYPE(GroupInfo_Type) :: GroupInfo, GcGroupInfo

INTEGER :: ng, ngc
INTEGER :: ig, igc, ib, ie, ibc, iec

ng = cuGeometry%ng
ngc = cuGeometry%ngc

IF (.NOT. ASSOCIATED(cuGeometry%GcStruct)) THEN
  ALLOCATE(cuGeometry%GcStruct(2, 2))
  cuGeometry%GcStruct(1, 1) = 1
  cuGeometry%GcStruct(2, 1) = GroupInfo%UpScatRange(1) - 1
  cuGeometry%GcStruct(1, 2) = GroupInfo%UpScatRange(1)
  cuGeometry%GcStruct(2, 2) = GroupInfo%UpScatRange(2)
ENDIF
ALLOCATE(cuGeometry%GcStructInv(ng))
DO igc = 1, ngc
  cuGeometry%GcStructInv(cuGeometry%GcStruct(1, igc) : cuGeometry%GcStruct(2, igc)) = igc
ENDDO
ALLOCATE(GcGroupInfo%InScatRange(2, ngc))
GcGroupInfo%InScatRange(1, :) = ngc
GcGroupInfo%InScatRange(2, :) = 1
GcGroupInfo%ng = ngc
DO ig = 1, ng
  igc = cuGeometry%GcStructInv(ig)
  ib = GroupInfo%InScatRange(1, ig); ie = GroupInfo%InScatRange(2, ig)
  ibc = cuGeometry%GcStructInv(ib); iec = cuGeometry%GcStructInv(ie)
  GcGroupInfo%InScatRange(1, igc) = MIN(GcGroupInfo%InScatRange(1, igc), ibc)
  GcGroupInfo%InScatRange(2, igc) = MAX(GcGroupInfo%InScatRange(2, igc), iec)
ENDDO
GcGroupInfo%lUpScat = .FALSE.
DO igc = 1, ngc
  IF (GcGroupInfo%InScatRange(2, igc) .GT. igc) THEN
    GcGroupInfo%lUpScat = .TRUE.
    EXIT
  ENDIF
ENDDO
IF (GcGroupInfo%lUpScat) THEN
  GcGroupInfo%UpScatRange(1) = ngc; GcGroupInfo%UpScatRange(2) = 1
  DO igc = 1, ngc
    IF (GcGroupInfo%InScatRange(2, igc) .GT. igc) THEN
      GcGroupInfo%UpScatRange(1) = MIN(GcGroupInfo%UpScatRange(1), igc)
      GcGroupInfo%UpScatRange(2) = MAX(GcGroupInfo%UpScatRange(2), igc)
    ENDIF
  ENDDO
  GcGroupInfo%UpScatRange(2) = ngc
ENDIF

END SUBROUTINE

SUBROUTINE AllocMOCVar(Core, RayInfo, cuMOC, cuDevice, ng, lScat1, lAsync, lXS)
USE TYPEDEF,		ONLY : CoreInfo_Type,		RayInfo_Type
USE CNTL,			ONLY : nTracerCntl
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo
TYPE(cuMOC_Type) :: cuMOC
TYPE(cuDevice_Type) :: cuDevice
INTEGER :: ng
LOGICAL :: lScat1, lAsync, lXS

INTEGER :: nfxr, nfsr, nxy, nxya, nAziAngle, nPolarAngle, nModRay, nPhiAngSv, ScatOd
INTEGER :: Moment(3) = (/ 2, 5, 9 /)
INTEGER :: i, j

nfxr = cuGeometry%nfxr
nfsr = cuGeometry%nfsr
nxy = cuGeometry%nxy
nxya = cuGeometry%nxya
nAziAngle = RayInfo%nAziAngle
nPolarAngle = RayInfo%nPolarAngle
nModRay = RayInfo%nModRay
nPhiAngSv = RayInfo%nPhiAngSv
ScatOd = nTracerCntl%ScatOd

IF (.NOT. lScat1) THEN
  IF (lAsync) THEN
    ALLOCATE(cuMOC%phisMg(P0_BLOCK_SIZE, nfsr))
    ALLOCATE(cuMOC%PhiAngInMg(nPolarAngle, P0_BLOCK_SIZE, nPhiAngSv))
    ALLOCATE(cuMOC%joutMg(3, P0_BLOCK_SIZE, 4, nxy))
    ALLOCATE(cuMOC%xstMg(P0_BLOCK_SIZE, nfsr))
    ALLOCATE(cuMOC%srcMg(P0_BLOCK_SIZE, nfsr))
  ELSE
    ALLOCATE(cuMOC%phis(ng, nfsr))
    ALLOCATE(cuMOC%psi(nfsr))
    ALLOCATE(cuMOC%src(ng, nfsr))
    ALLOCATE(cuMOC%PhiAngIn(nPolarAngle, ng, nPhiAngSv))
    ALLOCATE(cuMOC%jout(3, ng, 4, nxy))
    ALLOCATE(cuMOC%xst(ng, nfsr))
    IF (lXS) THEN
      ALLOCATE(cuMOC%xssm(ng, ng, nfxr))
      ALLOCATE(cuMOC%chi(ng, nfxr))
    ENDIF
  ENDIF
ELSE
  ALLOCATE(cuMOC%phisMg(PN_BLOCK_SIZE, nfsr))
  ALLOCATE(cuMOC%phimMg(Moment(ScatOd), PN_BLOCK_SIZE, nfsr))
  ALLOCATE(cuMOC%phiaMg(nPolarAngle, PN_BLOCK_SIZE, 4, nfsr))
  ALLOCATE(cuMOC%PhiAngInMg(nPolarAngle, PN_BLOCK_SIZE, nPhiAngSv))
  ALLOCATE(cuMOC%joutMg(3, PN_BLOCK_SIZE, 4, nxy))
  ALLOCATE(cuMOC%xstMg(PN_BLOCK_SIZE, nfsr))
  ALLOCATE(cuMOC%srcMg(PN_BLOCK_SIZE, nfsr))
  ALLOCATE(cuMOC%srcmMg(Moment(ScatOd), PN_BLOCK_SIZE, nfsr))
  ALLOCATE(cuMOC%SrcAngMg(nPolarAngle, PN_BLOCK_SIZE, 4, nfsr))
ENDIF

END SUBROUTINE

SUBROUTINE DeallocMOCVar(cuMOC, lScat1, lAsync, lXS)

IMPLICIT NONE

TYPE(cuMOC_Type) :: cuMOC
LOGICAL :: lScat1, lAsync, lXS

DEALLOCATE(cuMOC%phis)
DEALLOCATE(cuMOC%xst)
DEALLOCATE(cuMOC%src)
DEALLOCATE(cuMOC%PhiAngIn)

IF (.NOT. lScat1) THEN
  IF (lAsync) THEN
    DEALLOCATE(cuMOC%phisMg)
    DEALLOCATE(cuMOC%PhiAngInMg)
    DEALLOCATE(cuMOC%joutMg)
    DEALLOCATE(cuMOC%xstMg)
    DEALLOCATE(cuMOC%srcMg)
  ELSE
    DEALLOCATE(cuMOC%phis)
    DEALLOCATE(cuMOC%psi)
    DEALLOCATE(cuMOC%src)
    DEALLOCATE(cuMOC%PhiAngIn)
    DEALLOCATE(cuMOC%jout)
    DEALLOCATE(cuMOC%xst)
    IF (lXS) THEN
      DEALLOCATE(cuMOC%xssm)
      DEALLOCATE(cuMOC%chi)
    ENDIF
  ENDIF
ELSE
  DEALLOCATE(cuMOC%phisMg)
  DEALLOCATE(cuMOC%phimMg)
  DEALLOCATE(cuMOC%phiaMg)
  DEALLOCATE(cuMOC%PhiAngInMg)
  DEALLOCATE(cuMOC%joutMg)
  DEALLOCATE(cuMOC%xstMg)
  DEALLOCATE(cuMOC%srcMg)
  DEALLOCATE(cuMOC%srcmMg)
  DEALLOCATE(cuMOC%SrcAngMg)
ENDIF

END SUBROUTINE

SUBROUTINE AllocCMFDVar(cuCMFD, cuDevice, GroupInfo, lGcCMFD)
USE TYPEDEF,        ONLY : GroupInfo_Type
USE CMFD_COMMON,    ONLY : AllocHomoXsVar,      AllocPinXS
IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
TYPE(GroupInfo_Type) :: GroupInfo
LOGICAL :: lGcCMFD

INTEGER :: ng, nxy, nxyRB(100, 2), nzCMFD
INTEGER :: iz, myzb, myze, myzbf, myzef

ng = cuCMFD%ng
nxy = cuGeometry%nxyc
nxyRB = cuDevice%nxyRB
nzCMFD = cuDevice%nzCMFD
myzb = cuDevice%myzb
myze = cuDevice%myze
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef

ALLOCATE(cuCMFD%phis8(ng, nxy * nzCMFD))
ALLOCATE(cuCMFD%src8(ng, nxy * nzCMFD))
ALLOCATE(cuCMFD%psi8(nxy * nzCMFD)); cuCMFD%psi8 = 0.0
ALLOCATE(cuCMFD%psid8(nxy * nzCMFD)); cuCMFD%psid8 = 0.0

ALLOCATE(cuCMFD%h_phis8(ng, nxy, myzbf : myzef))
ALLOCATE(cuCMFD%h_phic8(ng, nxy, myzb : myze))

IF (cuCntl%lNatural) THEN
  IF (cuCntl%CMFDSolver .EQ. 2) THEN
    ALLOCATE(cuCMFD%invDiag(ng, nxy * nzCMFD)); cuCMFD%invDiag = 0.0
  ENDIF
ENDIF

IF (cuGeometry%l3dim) THEN
  ALLOCATE(cuCMFD%h_neighphis8(ng, nxy, 2))
  ALLOCATE(cuCMFD%AxDtil(2, ng, nxy, myzbf : myzef))
  ALLOCATE(cuCMFD%AxDhat(2, ng, nxy, myzbf : myzef)); cuCMFD%AxDhat = 0.0
  ALLOCATE(cuCMFD%offDiag(ng, nxy, 2)); cuCMFD%offDiag = 0.0
  ALLOCATE(cuCMFD%offDiag8(ng, nxy, 2)); cuCMFD%offDiag8 = 0.0
  IF (.NOT. cuCntl%lNatural) THEN
    ALLOCATE(cuCMFD%rOffDiag(ng, nxyRB(myzbf, red) + nxyRB(myzef, red))); cuCMFD%rOffDiag = 0.0
    ALLOCATE(cuCMFD%bOffDiag(ng, nxyRB(myzbf, black) + nxyRB(myzef, black))); cuCMFD%bOffDiag = 0.0
  ENDIF
ENDIF

cuCMFD%bottomRange = (/ (cuDevice%bottomRange(1) - 1) * ng + 1, cuDevice%bottomRange(2) * ng /)
cuCMFD%bottomRangeRB(:, 1) = (/ (cuDevice%bottomRangeRB(1, 1) - 1) * ng + 1, cuDevice%bottomRangeRB(2, 1) * ng /)
cuCMFD%bottomRangeRB(:, 2) = (/ (cuDevice%bottomRangeRB(1, 2) - 1) * ng + 1, cuDevice%bottomRangeRB(2, 2) * ng /)
cuCMFD%topRange = (/ (cuDevice%topRange(1) - 1) * ng + 1, cuDevice%topRange(2) * ng /)
cuCMFD%topRangeRB(:, 1) = (/ (cuDevice%topRangeRB(1, 1) - 1) * ng + 1, cuDevice%topRangeRB(2, 1) * ng /)
cuCMFD%topRangeRB(:, 2) = (/ (cuDevice%topRangeRB(1, 2) - 1) * ng + 1, cuDevice%topRangeRB(2, 2) * ng /)

IF (.NOT. lGcCMFD) THEN
  cuCMFD%planeMap => cuGeometry%planeMap
  CALL AllocPinXS(cuCMFD%PinXS, GroupInfo, nxy, myzb, myze)
ELSE
  ALLOCATE(cuCMFD%planeMap(myzbf : myzef))
  DO iz = myzbf, myzef
    cuCMFD%planeMap(iz) = iz
  ENDDO
  CALL AllocPinXS(cuCMFD%PinXS, GroupInfo, nxy, myzbf, myzef)
ENDIF

END SUBROUTINE

SUBROUTINE DeallocCMFDVar(cuCMFD)

IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD

CALL destroyCsr(cuCMFD%M)
CALL destroyCsr(cuCMFD%jM)
CALL destroyCsr(cuCMFD%rbM(red))
CALL destroyCsr(cuCMFD%rbM(black))
CALL destroyCsr(cuCMFD%rbDiag(red))
CALL destroyCsr(cuCMFD%rbDiag(black))
CALL destroyCsr(cuCMFD%D)
CALL destroyCsr(cuCMFD%S)
CALL destroyCsr(cuCMFD%F)
CALL destroyCsr(cuCMFD%Chi)

END SUBROUTINE

SUBROUTINE AllocAxialVar(cuAxial, cuDevice)
USE CUDA_NODAL,     ONLY : cuAllocNodal
USE CUDA_AXMOC,     ONLY : cuAllocFlatMOC,      cuAllocLinearMOC
IMPLICIT NONE

TYPE(cuAxial_Type) :: cuAxial
TYPE(cuDevice_Type) :: cuDevice

INTEGER :: ng, nxy, nzCMFD, nPolar1D
INTEGER :: myzb, myze, myzbf, myzef

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
myzb = cuDevice%myzb
myze = cuDevice%myze
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef

ALLOCATE(cuAxial%phic(ng, nxy, myzbf : myzef)); cuAxial%phic = 0.0
ALLOCATE(cuAxial%Jout(ng, nxy, myzbf : myzef, 2, 2)); cuAxial%Jout = 0.0

ALLOCATE(cuAxial%S0(nxy, ng, myzb : myze, ng)); cuAxial%S0 = 0.0
ALLOCATE(cuAxial%S1(nxy, ng, myzb : myze, ng)); cuAxial%S1 = 0.0
ALLOCATE(cuAxial%Fa(nxy, myzb : myze, ng)); cuAxial%Fa = 0.0
ALLOCATE(cuAxial%Chia(nxy, ng, myzb : myze)); cuAxial%Chia = 0.0

SELECT CASE (cuCntl%AxSolver)
CASE (NODAL)
  CALL cuAllocNodal(cuAxial, cuDevice)
CASE (MOC)
  CALL SetAngle()
  IF (cuCntl%lCASMO) THEN
    CALL cuAllocLinearMOC(cuAxial, cuDevice)
  ELSE
    CALL cuAllocFlatMOC(cuAxial, cuDevice)
  ENDIF
END SELECT

CONTAINS

SUBROUTINE SetAngle()
USE CUDA_CONST
IMPLICIT NONE

INTEGER :: ipol
REAL :: mu
REAL, ALLOCATABLE :: abscissa(:), weight(:)

ALLOCATE(abscissa(cuGeometry%nPolar1D * 2), weight(cuGeometry%nPolar1D * 2))

CALL gauleg(cuGeometry%nPolar1D * 2, abscissa, weight)

nPolar1D = cuGeometry%nPolar1D

DO ipol = 1, cuGeometry%nPolar1D
  mu = - abscissa(ipol)
  cosv1D(ipol) = mu
  rcosv1D(ipol) = 1.0 / mu
  wt1D(ipol) = weight(ipol) / 2.0
  wtsurf1D(ipol) = weight(ipol) * mu / 2.0
  Comp1D(ipol) = mu
  mwt1D(ipol) = mu * (weight(ipol) * mu / 2.0)
ENDDO

DEALLOCATE(abscissa, weight)

END SUBROUTINE

!********************************************************************************
!* Calculation of GAUSS-LEGENDRE abscissas and weights for Gaussian Quadrature
!* integration of polynomial functions.
!*      For normalized lower and upper limits of integration -1.0 & 1.0, and
!* given n, this routine calculates, arrays xabsc(1:n) and  weig(1:n) of length n,
!* containing the abscissas and weights of the Gauss-Legendre n-point quadrature
!* formula.  For detailed explanations finding weights & abscissas, see
!* "Numerical Recipes in Fortran */
!********************************************************************************
	SUBROUTINE  gauleg(ngp, xabsc, weig)

      implicit none
      INTEGER, parameter :: dbp = SELECTED_REAL_KIND (15,307)
      INTEGER  i, j, m
      REAL(dbp)  p1, p2, p3, pp, z, z1
      INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
      REAL(dbp), INTENT(OUT) :: xabsc(ngp), weig(ngp)
      REAL(dbp)  :: EPS, M_PI
      PARAMETER (EPS=3.0d-15)       	!EPS is the relative precision
      PARAMETER (M_PI=3.141592654d0)      ! Pi value

	   m = (ngp + 1) / 2
!* Roots are symmetric in the interval - so only need to find half of them  */

	   do i = 1, m				! Loop over the desired roots */

     		z = cos( M_PI * (i-0.25d0) / (ngp+0.5d0) )
!*   Starting with the above approximation to the ith root,
!*          we enter the main loop of refinement by NEWTON'S method   */
100     	p1 = 1.0d0
        	p2 = 0.0d0
!*  Loop up the recurrence relation to get the Legendre
!*  polynomial evaluated at z                 */

        	do j = 1, ngp
           	p3 = p2
           	p2 = p1
           	p1 = ((2.0d0*j-1.0d0) * z * p2 - (j-1.0d0)*p3) / j
        	enddo

!* p1 is now the desired Legendre polynomial. We next compute pp,
!* its derivative, by a standard relation involving also p2, the
!* polynomial of one lower order.      */
        	pp = ngp*(z*p1-p2)/(z*z-1.0d0)
        	z1 = z
        	z = z1 - p1/pp             ! Newton's Method  */

        	if (dabs(z-z1) .gt. EPS) GOTO  100

      	xabsc(i) =  - z                    	! Roots will be bewteen -1.0 & 1.0 */
      	xabsc(ngp+1-i) =  + z                	! and symmetric about the origin  */
      	weig(i) = 2.0d0/((1.0d0-z*z)*pp*pp) ! Compute the weight and its       */
      	weig(ngp+1-i) = weig(i)               ! symmetric counterpart         */

      end do     ! i loop

   End subroutine gauleg

END SUBROUTINE

SUBROUTINE CopyCoreVar(Core, cuDevice)
USE TYPEDEF,		ONLY : CoreInfo_Type
USE PE_MOD,			ONLY : PE
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(cuDevice_Type) :: cuDevice

!$ACC ENTER DATA COPYIN(cuGeometry) ASYNC(0)
!$ACC ENTER DATA COPYIN(cuGeometry%fsrVol) ASYNC(0)
!$ACC ENTER DATA COPYIN(cuGeometry%hzfm) ASYNC(0)
!$ACC ENTER DATA COPYIN(cuGeometry%AxType) ASYNC(0)
!$ACC ENTER DATA COPYIN(cuGeometry%Fsr2Fxr) ASYNC(0)
!$ACC ENTER DATA COPYIN(cuGeometry%lRefCell) ASYNC(0)
!$ACC ENTER DATA COPYIN(cuGeometry%lRefPin) ASYNC(0)
!$ACC ENTER DATA COPYIN(cuGeometry%lH2OCell) ASYNC(0)
!$ACC ENTER DATA COPYIN(cuGeometry%RotRayCount) ASYNC(0)
!$ACC ENTER DATA COPYIN(cuGeometry%RotRayList) ASYNC(0)

END SUBROUTINE

SUBROUTINE CopyControlVar()

IMPLICIT NONE

!$ACC ENTER DATA COPYIN(cuDevice) ASYNC(0)
!$ACC ENTER DATA COPYIN(cuDevice%cmMap) ASYNC(0)
!$ACC ENTER DATA COPYIN(cuDevice%fmMap) ASYNC(0)
!$ACC ENTER DATA COPYIN(cuDevice%cmSubrange) ASYNC(0)
!$ACC ENTER DATA COPYIN(cuDevice%fmSubrange) ASYNC(0)
!$ACC ENTER DATA COPYIN(cuDevice%hzSub) ASYNC(0)

END SUBROUTINE

SUBROUTINE CopyRayVar(iAx)
USE TYPEDEF,		ONLY : CoreInfo_Type,		RayInfo_Type
USE CNTL,			ONLY : nTracerCntl
USE PE_MOD,			ONLY : PE
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo

INTEGER :: iAx

!$ACC ENTER DATA PRESENT_OR_COPYIN(cuFastRay1D) ASYNC(0)
!$ACC ENTER DATA COPYIN(cuFastRay1D(iAx)%PhiAngInSvIdx) ASYNC(0)
!$ACC ENTER DATA COPYIN(cuFastRay1D(iAx)%PhiAngOutSvIdx) ASYNC(0)
!$ACC ENTER DATA COPYIN(cuFastRay1D(iAx)%RotRayIdxSt) ASYNC(0)
!$ACC ENTER DATA COPYIN(cuFastRay1D(iAx)%CoreRayIdxSt) ASYNC(0)
!$ACC ENTER DATA COPYIN(cuFastRay1D(iAx)%PinRayIdxSt) ASYNC(0)
!$ACC ENTER DATA COPYIN(cuFastRay1D(iAx)%iAzi) ASYNC(0)
!$ACC ENTER DATA COPYIN(cuFastRay1D(iAx)%iDir) ASYNC(0)
!$ACC ENTER DATA COPYIN(cuFastRay1D(iAx)%PinIdx) ASYNC(0)
!$ACC ENTER DATA COPYIN(cuFastRay1D(iAx)%SurfIdx) ASYNC(0)
!$ACC ENTER DATA COPYIN(cuFastRay1D(iAx)%FsrIdx) ASYNC(0)
!$ACC ENTER DATA COPYIN(cuFastRay1D(iAx)%LenSeg) ASYNC(0)

END SUBROUTINE

SUBROUTINE DeleteRayVar(iAx)
USE TYPEDEF,		ONLY : CoreInfo_Type,		RayInfo_Type
USE CNTL,			ONLY : nTracerCntl
USE PE_MOD,			ONLY : PE
IMPLICIT NONE

INTEGER :: iAx

!$ACC EXIT DATA DELETE(cuFastRay1D(iAx)%PhiAngInSvIdx) ASYNC(0)
!$ACC EXIT DATA DELETE(cuFastRay1D(iAx)%PhiAngOutSvIdx) ASYNC(0)
!$ACC EXIT DATA DELETE(cuFastRay1D(iAx)%RotRayIdxSt) ASYNC(0)
!$ACC EXIT DATA DELETE(cuFastRay1D(iAx)%CoreRayIdxSt) ASYNC(0)
!$ACC EXIT DATA DELETE(cuFastRay1D(iAx)%PinRayIdxSt) ASYNC(0)
!$ACC EXIT DATA DELETE(cuFastRay1D(iAx)%iAzi) ASYNC(0)
!$ACC EXIT DATA DELETE(cuFastRay1D(iAx)%iDir) ASYNC(0)
!$ACC EXIT DATA DELETE(cuFastRay1D(iAx)%PinIdx) ASYNC(0)
!$ACC EXIT DATA DELETE(cuFastRay1D(iAx)%SurfIdx) ASYNC(0)
!$ACC EXIT DATA DELETE(cuFastRay1D(iAx)%FsrIdx) ASYNC(0)
!$ACC EXIT DATA DELETE(cuFastRay1D(iAx)%LenSeg) ASYNC(0)

END SUBROUTINE

END MODULE

#endif
