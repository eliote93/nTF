module SetCentGeom
USE param
USE allocs
USE cntl
implicit none
contains

SUBROUTINE SetCentCell(icel0, icel, iType)
use geom,   only :  CellInfo,  CellPitch,    nCellType
use BenchXs,only :  MacXsBen
USE XSLIB_MOD ,ONLY : igresb,igrese
use SPH_Mod, only :  calcCellSSPH
IMPLICIT NONE

INTEGER,INTENT(IN) :: icel0, icel, iType     !Cell Index
INTEGER :: nFXR, nFSR, nx, ny, nFsrInFxr     !
INTEGER :: nCircle, nmat, nreg
INTEGER :: i, j, k, m, n, ir0, ir, ireg
LOGICAL :: lRect                             !Rectangular FSR Flag
REAL :: temp
REAL :: HalfPitch                            !Half Pitch Cell
REAL :: lx(3), ly(3)                         !Cell Pitch Information

REAL :: DelX(100), DelY(100)
INTEGER :: FxrList(100), FxrMap(1000)

CellInfo(icel)%lempty = CellInfo(icel0)%lempty; CellInfo(icel0)%EdgeCellIdx(itype) = icel
IF(CellInfo(icel)%lempty) RETURN

CellInfo(icel)%icel0 = icel0
CellInfo(icel)%lRect = CellInfo(icel0)%lRect
CellInfo(icel)%lGap = CellInfo(icel0)%lGap
CellInfo(icel)%lFuel = CellInfo(icel0)%lFuel
CellInfo(icel)%lCCell = CellInfo(icel0)%lCCell
CellInfo(icel)%Geom%lCCent = CellInfo(icel0)%Geom%lCCent
CellInfo(icel)%lRes = CellInfo(icel0)%lRes
CellInfo(icel)%lMox = CellInfo(icel0)%lMox
CellInfo(icel)%lAIC = CellInfo(icel0)%lAIC
CellInfo(icel)%lCentX = FALSE; CellInfo(icel)%lCentY = FALSE; CellInfo(icel)%lCentXY = FALSE;
CellInfo(icel)%nbd = 4
lRect = CellInfo(icel0)%lRect
nFXR = CellInfo(icel0)%nFXR

SELECT CASE(iType)
  CASE(1)
    CellInfo(icel)%lCentX = TRUE
  CASE(2)
    CellInfo(icel)%lCentY = TRUE
  CASE(3)
    CellInfo(icel)%lCentXY = TRUE
  END SELECT

IF (CellInfo(icel)%Geom%lCCent) THEN
  CellInfo(icel)%Geom%nCircle = 0
  CellInfo(icel)%nFxr = 0
  RETURN
ENDIF

IF(lRect) THEN
  nx = CellInfo(icel0)%Geom%nx
  ny = CellInfo(icel0)%Geom%ny
  IF(itype .EQ. 1 .OR. itype .EQ. 3) THEN
    temp = CellPitch * 0.5_8; j = 0
    DO i = ny, 1, -1
      temp = temp - CellInfo(icel0)%Geom%Dely(i)
      j = j + 1
      IF(temp .LT. 1.0E-4_8) EXIT
    ENDDO
    ny = j
  ENDIF
  IF(itype .EQ. 2 .OR. itype .EQ. 3) THEN
    temp = CellPitch * 0.5_8; j = 0
    DO i = nx, 1, -1
      temp = temp - CellInfo(icel0)%Geom%Delx(i)
      j = j + 1
      IF(temp .LT. 1.0E-4_8) EXIT
    ENDDO
    nx = j
  ENDIF

  DO i = 1, CellInfo(icel0)%nFxr
    DO j = 1, CellInfo(iCel0)%nFsrInFxr(i)
      FxrMap(CellInfo(icel0)%MapFxr2FsrIdx(j, i)) = -i
    ENDDO
  ENDDO

  FxrList = 0
  DO j = CellInfo(icel0)%Geom%ny - ny + 1 , CellInfo(icel0)%Geom%ny
    DO i = CellInfo(icel0)%Geom%nx - nx + 1 , CellInfo(icel0)%Geom%nx
      k = CellInfo(icel0)%Geom%nx * (j - 1) + i
      FxrMap(k) = abs(FxrMap(k))
      FxrList(abs(FxrMap(k))) = 1
    ENDDO
  ENDDO

  i = 0
  DO j = 1, CellInfo(icel0)%nFxr
    IF(FxrList(j) .EQ. 1) THEN
      i = i+ 1; FxrList(j) = i
    ENDIF
  ENDDO
  nFxr = i

  DO j = CellInfo(icel0)%Geom%ny - ny + 1 , CellInfo(icel0)%Geom%ny
    DO i = CellInfo(icel0)%Geom%nx - nx + 1 , CellInfo(icel0)%Geom%nx
      k = CellInfo(icel0)%Geom%nx * (j - 1) + i
      FxrMap(k) = FxrList(FxrMap(k))
    ENDDO
  ENDDO

  DelX = 0; DelY = 0
  i = 0
  DO j = CellInfo(icel0)%Geom%ny - ny + 1 , CellInfo(icel0)%Geom%ny
    i = i + 1; DelY(i) = CellInfo(icel0)%Geom%DelY(j)
  ENDDO
  IF(itype .EQ. 1 .OR. itype .EQ. 3) THEN
    temp = sum(DelY(1:ny)) - CellPitch * 0.5
    IF(temp .GT. 1.E-4) THEN
      DelY(1) = DelY(1) - temp
    ENDIF
  ENDIF
  i = 0
  DO j = CellInfo(icel0)%Geom%nx - nx + 1 , CellInfo(icel0)%Geom%nx
    i = i + 1; DelX(i) = CellInfo(icel0)%Geom%DelX(j)
  ENDDO
  IF(itype .EQ. 2 .OR. itype .EQ. 3) THEN
    temp = sum(DelX(1:nx)) - CellPitch * 0.5
    IF(temp .GT. 1.E-4) THEN
      DelX(1) = DelX(1) - temp
    ENDIF
  ENDIF
  nFSR = nx*ny
ENDIF
IF(.NOT. lRECT) THEN
  nFSR = CellInfo(icel0)%nFSR/2
  IF(iType .EQ. 3) nFsr = nFsr/2
ENDIF

!ALLOCATION PART
CALL Dmalloc(CellInfo(icel)%iReg, nFSR)
CALL dmalloc(CellInfo(icel)%FxrIdxSt, nFxr)
CALL dmalloc(CellInfo(icel)%nFsrInFxr, nFxr)
CALL dmalloc(CellInfo(icel)%MapFxr2FsrIdx, nFSR, nFxr)

IF(lRect) THEN
  CALL Dmalloc(CellInfo(icel)%Geom%delx, nx)
  CALL Dmalloc(CellInfo(icel)%Geom%dely, ny)
ENDIF
HalfPitch = Half*CellPitch
lx = (/CellPitch, HalfPitch, HalfPitch/)
ly = (/HalfPitch, CellPitch, HalfPitch/)
lx = (/ONE,  HALF, HALF/)
ly = (/HALF, ONE,  HALF/)


!Center Line Pincell Setting
CellInfo(icel)%nFxr = nFxr;  Cellinfo(icel)%nFSR = nFSR
CellInfo(icel)%Geom%nbd = 4
CellInfo(icel)%Geom%lx = lx(itype)*CellInfo(icel0)%Geom%lx
CellInfo(icel)%Geom%ly = ly(itype)*CellInfo(icel0)%Geom%ly

CellInfo(icel)%Geom%cx = CellInfo(icel0)%Geom%cx
CellInfo(icel)%Geom%cy = CellInfo(icel0)%Geom%cy

SELECT CASE(iType)
  CASE(1)
    CellInfo(icel)%nDivAzi = CellInfo(icel0)%nDivAzi/2
    CellInfo(icel)%Geom%x(1) = -HalfPitch; CellInfo(icel)%Geom%x(2) = HalfPitch
    CellInfo(icel)%Geom%y(1) = -HalfPitch; CellInfo(icel)%Geom%y(2) = 0
!     ny = nint(0.5*ny)
  CASE(2)
    CellInfo(icel)%nDivAzi = CellInfo(icel0)%nDivAzi/2
    CellInfo(icel)%Geom%x(1) = 0; CellInfo(icel)%Geom%x(2) = HalfPitch
    CellInfo(icel)%Geom%y(1) = -HalfPitch; CellInfo(icel)%Geom%y(2) = HalfPitch
  CASE(3)
    CellInfo(icel)%nDivAzi = CellInfo(icel0)%nDivAzi/4
    CellInfo(icel)%Geom%x(1) = 0; CellInfo(icel)%Geom%x(2) = HalfPitch
    CellInfo(icel)%Geom%y(1) = -HalfPitch; CellInfo(icel)%Geom%y(2) = 0
    !CellInfo(icel)%Geom%x(1) = 0; CellInfo(icel)%Geom%x(2) = HalfPitch
    !CellInfo(icel)%Geom%y(1) = 0; CellInfo(icel)%Geom%y(2) = HalfPitch
END SELECT

IF(CellInfo(icel)%lGap) THEN
  lx = (/ONE,  HALF, HALF/)
  ly = (/HALF, ONE,  HALF/)
  CellInfo(icel)%Geom%lx = lx(itype)*CellInfo(icel0)%Geom%lx
  CellInfo(icel)%Geom%ly = ly(itype)*CellInfo(icel0)%Geom%ly
  CellInfo(icel)%Geom%x = CellInfo(icel0)%Geom%x
  CellInfo(icel)%Geom%y = CellInfo(icel0)%Geom%y

  SELECT CASE(iType)
    CASE(1)
      !CellInfo(icel)%Geom%x(1) = -HalfPitch; CellInfo(icel)%Geom%x(2) = HalfPitch
      CellInfo(icel)%Geom%y(1) = -CellInfo(icel)%Geom%ly; CellInfo(icel)%Geom%y(2) = 0
    CASE(2)
      CellInfo(icel)%Geom%x(1) = 0; CellInfo(icel)%Geom%x(2) = CellInfo(icel)%Geom%lx
    CASE(3)
      CellInfo(icel)%Geom%x(1) = 0; CellInfo(icel)%Geom%x(2) = CellInfo(icel)%Geom%lx
      CellInfo(icel)%Geom%y(1) = -CellInfo(icel)%Geom%ly; CellInfo(icel)%Geom%y(2) = 0
  END SELECT
ENDIF

CellInfo(icel)%lRect = lRECT;
CellInfo(icel)%Geom%lRect = lRect; CellInfo(icel)%Geom%lCircle = .NOT. lRect
CellInfo(icel)%Geom%nCircle = 0
IF(.NOT. lRect) THEN
  nFsrInFxr = CellInfo(icel)%nDivAzi
  !Geometries Setting
  CellInfo(iCel)%Geom%ncircle = CellInfo(icel0)%Geom%nCircle
  nCircle = CellInfo(iCel)%Geom%ncircle
  CellInfo(iCel)%Geom%cx = CellInfo(icel0)%Geom%cx
  CellInfo(iCel)%Geom%cy = CellInfo(icel0)%Geom%cy
  Call Dmalloc(CellInfo(icel)%Geom%circle, 3, nCircle)
  CellInfo(icel)%Geom%circle(3, :) = CellInfo(icel0)%Geom%circle(3, :)
  !Index Setting - FSR, FXR index setting
  ir = 0
  DO i = 1, nFXR
    ir0 = CellInfo(icel0)%MapFxr2FsrIdx(1,i)
    ireg = CellInfo(icel0)%ireg(ir0)
    CellInfo(icel)%ireg((nFsrinFxr*(i-1)+1):(nFsrInFxr*i)) = ireg
    CellInfo(icel)%nFsrInFxr(i) = nFsrInFxr
    CellInfo(icel)%FxrIdxSt(i) = nFsrInFxr*(i-1)+1
    DO j = 1,nFsrInFxr
      CellInfo(icel)%MapFxr2FsrIdx(j,i) = CellInfo(icel)%FxrIdxSt(i) + j - 1
    ENDDO
  ENDDO
  
  nmat = CellInfo(icel)%nmat
  allocate(CellInfo(icel)%matidx(nmat))
  CellInfo(icel)%matidx(1:nmat) = CellInfo(icel0)%matidx(1:nmat)
  allocate(CellInfo(icel)%matrad(nmat))
  CellInfo(icel)%matrad(1:nmat) = CellInfo(icel0)%matrad(1:nmat)
  
  nreg = CellInfo(icel0)%nreg_cp
  CellInfo(icel)%nreg_cp = nreg
  
  allocate(CellInfo(icel)%rad_cp(nreg),CellInfo(icel)%fxrvol(nreg))
  CellInfo(icel)%rad_cp(1:nreg) = CellInfo(icel0)%rad_cp(1:nreg)
  CellInfo(icel)%fxrvol(1:nreg) = CellInfo(icel0)%fxrvol(1:nreg)
  CellInfo(icel)%nfueldiv = CellInfo(icel0)%nfueldiv
  CellInfo(icel)%fuelgapcldvol = CellInfo(icel0)%fuelgapcldvol
  CellInfo(icel)%cldfxridx = CellInfo(icel0)%cldfxridx
  CellInfo(icel)%nmodfxr = CellInfo(icel0)%nmodfxr
  CellInfo(icel)%invnmodfxr = CellInfo(icel0)%invnmodfxr
  CellInfo(icel)%FuelRad0 = CellInfo(icel0)%FuelRad0
  CellInfo(icel)%lsSPH = CellInfo(icel0)%lsSPH
  
  IF (nTracerCntl%lSSPH) then
      IF (CellInfo(icel)%lfuel.or.CellInfo(icel)%lAIC) then
          ALLOCATE(CellInfo(icel)%SPHfactor(1:CellInfo(icel0)%nFXR,igresb:igrese))
          CellInfo(icel)%SPHfactor = CellInfo(icel0)%SPHfactor
          CellInfo(icel)%ngapdiv = CellInfo(icel0)%ngapdiv
          CellInfo(icel)%ncladdiv = CellInfo(icel0)%ncladdiv
          CellInfo(icel)%FuelRad0 = CellInfo(icel0)%FuelRad0
          CellInfo(icel)%srdidx = CellInfo(icel0)%srdidx
      ENDIF
  ENDIF
ELSE
  nFsrInFxr = nFSR
  CellInfo(icel)%Geom%nx = nx
  CellInfo(icel)%Geom%ny = ny
  ir = 0
  DO j = CellInfo(icel0)%Geom%ny - ny + 1 , CellInfo(icel0)%Geom%ny
    DO i = CellInfo(icel0)%Geom%nx - nx + 1 , CellInfo(icel0)%Geom%nx
      k = CellInfo(icel0)%Geom%nx * (j - 1) + i; ir = ir + 1
      ireg = CellInfo(icel0)%iReg(k)
      CellInfo(icel)%iReg(ir) = ireg
    ENDDO
  ENDDO

  DO m = 1, nFxr
    ir = 0; j = 0; k =0
    DO i = 1, CellInfo(icel0)%nFsr
      IF(FxrMap(i) .LT. 0) CYCLE
      j = j + 1
      IF(FxrMap(i) .NE. m) CYCLE
      !CellInfo(icel)%nFsrInFxr(m) = CellInfo(icel)%nFsrInFxr(m) + 1
      k = k + 1; CellInfo(icel)%nFsrInFxr(m) = k
      CellInfo(icel)%MapFxr2FsrIdx(k, m) = j
    ENDDO
    CellInfo(icel)%FxrIdxSt(m) = CellInfo(icel)%MapFxr2FsrIdx(1, m)
  ENDDO

  CellInfo(icel)%Geom%Delx(1:nx) = DelX(1:nx)
  CellInfo(icel)%Geom%Dely(1:ny) = DelY(1:ny)
  CONTINUE
!  iReg = CellInfo(icel0)%Ireg(1);  CellInfo(icel)%iReg(1:nFSR) = iReg   !Mixture or Composition Number
!  CellInfo(icel)%nFsrInFxr(:) = nFsrInFxr;  !CellInfo(icel)%FxrIdxst(1) = 1
!  DO j = 1, nFSR
!    CellInfo(icel)%MapFxr2FsrIdx(j, 1) = j
!  ENDDO
ENDIF

END SUBROUTINE

!--- JSR Edit : nTIG Restart

SUBROUTINE SetCentCellBase(icel0, icel, iType)
use geom,   only :  BaseCellInfo,  CellPitch,    nCellType
use BenchXs,only :  MacXsBen
USE XSLIB_MOD ,ONLY : igresb,igrese
use SPH_Mod, only :  calcCellSSPH
IMPLICIT NONE

INTEGER,INTENT(IN) :: icel0, icel, iType     !Cell Index
INTEGER :: nFXR, nFSR, nx, ny, nFsrInFxr     !
INTEGER :: nCircle
INTEGER :: i, j, k, m, n, ir0, ir
LOGICAL :: lRect                             !Rectangular FSR Flag
REAL :: temp
REAL :: HalfPitch                            !Half Pitch Cell
REAL :: lx(3), ly(3)                         !Cell Pitch Information

REAL :: DelX(100), DelY(100)
INTEGER :: FxrList(100), FxrMap(1000)

BaseCellInfo(icel)%lempty = BaseCellInfo(icel0)%lempty; BaseCellInfo(icel0)%EdgeCellIdx(itype) = icel
IF(BaseCellInfo(icel)%lempty) RETURN

BaseCellInfo(icel)%lRect = BaseCellInfo(icel0)%lRect
BaseCellInfo(icel)%lGap = BaseCellInfo(icel0)%lGap
BaseCellInfo(icel)%lFuel = BaseCellInfo(icel0)%lFuel
BaseCellInfo(icel)%lCCell = BaseCellInfo(icel0)%lCCell
BaseCellInfo(icel)%Geom%lCCent = BaseCellInfo(icel0)%Geom%lCCent
BaseCellInfo(icel)%lRes = BaseCellInfo(icel0)%lRes
BaseCellInfo(icel)%lMox = BaseCellInfo(icel0)%lMox
BaseCellInfo(icel)%lCentX = FALSE; BaseCellInfo(icel)%lCentY = FALSE; BaseCellInfo(icel)%lCentXY = FALSE;
BaseCellInfo(icel)%nbd = 4
lRect = BaseCellInfo(icel0)%lRect
nFXR = BaseCellInfo(icel0)%nFXR

SELECT CASE(iType)
  CASE(1)
    BaseCellInfo(icel)%lCentX = TRUE
  CASE(2)
    BaseCellInfo(icel)%lCentY = TRUE
  CASE(3)
    BaseCellInfo(icel)%lCentXY = TRUE
  END SELECT

IF (BaseCellInfo(icel)%Geom%lCCent) THEN
  BaseCellInfo(icel)%Geom%nCircle = 0
  BaseCellInfo(icel)%nFxr = 0
  RETURN
ENDIF

IF(lRect) THEN
  nx = BaseCellInfo(icel0)%Geom%nx
  ny = BaseCellInfo(icel0)%Geom%ny
  IF(itype .EQ. 1 .OR. itype .EQ. 3) THEN
    temp = CellPitch * 0.5_8; j = 0
    DO i = ny, 1, -1
      temp = temp - BaseCellInfo(icel0)%Geom%Dely(i)
      j = j + 1
      IF(temp .LT. 1.0E-4_8) EXIT
    ENDDO
    ny = j
  ENDIF
  IF(itype .EQ. 2 .OR. itype .EQ. 3) THEN
    temp = CellPitch * 0.5_8; j = 0
    DO i = nx, 1, -1
      temp = temp - BaseCellInfo(icel0)%Geom%Delx(i)
      j = j + 1
      IF(temp .LT. 1.0E-4_8) EXIT
    ENDDO
    nx = j
  ENDIF

  DO i = 1, BaseCellInfo(icel0)%nFxr
    DO j = 1, BaseCellInfo(iCel0)%nFsrInFxr(i)
      FxrMap(BaseCellInfo(icel0)%MapFxr2FsrIdx(j, i)) = -i
    ENDDO
  ENDDO

  FxrList = 0
  DO j = BaseCellInfo(icel0)%Geom%ny - ny + 1 , BaseCellInfo(icel0)%Geom%ny
    DO i = BaseCellInfo(icel0)%Geom%nx - nx + 1 , BaseCellInfo(icel0)%Geom%nx
      k = BaseCellInfo(icel0)%Geom%nx * (j - 1) + i
      FxrMap(k) = abs(FxrMap(k))
      FxrList(abs(FxrMap(k))) = 1
    ENDDO
  ENDDO

  i = 0
  DO j = 1, BaseCellInfo(icel0)%nFxr
    IF(FxrList(j) .EQ. 1) THEN
      i = i+ 1; FxrList(j) = i
    ENDIF
  ENDDO
  nFxr = i

  DO j = BaseCellInfo(icel0)%Geom%ny - ny + 1 , BaseCellInfo(icel0)%Geom%ny
    DO i = BaseCellInfo(icel0)%Geom%nx - nx + 1 , BaseCellInfo(icel0)%Geom%nx
      k = BaseCellInfo(icel0)%Geom%nx * (j - 1) + i
      FxrMap(k) = FxrList(FxrMap(k))
    ENDDO
  ENDDO

  DelX = 0; DelY = 0
  i = 0
  DO j = BaseCellInfo(icel0)%Geom%ny - ny + 1 , BaseCellInfo(icel0)%Geom%ny
    i = i + 1; DelY(i) = BaseCellInfo(icel0)%Geom%DelY(j)
  ENDDO
  IF(itype .EQ. 1 .OR. itype .EQ. 3) THEN
    temp = sum(DelY(1:ny)) - CellPitch * 0.5
    IF(temp .GT. 1.E-4) THEN
      DelY(1) = DelY(1) - temp
    ENDIF
  ENDIF
  i = 0
  DO j = BaseCellInfo(icel0)%Geom%nx - nx + 1 , BaseCellInfo(icel0)%Geom%nx
    i = i + 1; DelX(i) = BaseCellInfo(icel0)%Geom%DelX(j)
  ENDDO
  IF(itype .EQ. 2 .OR. itype .EQ. 3) THEN
    temp = sum(DelX(1:nx)) - CellPitch * 0.5
    IF(temp .GT. 1.E-4) THEN
      DelX(1) = DelX(1) - temp
    ENDIF
  ENDIF
  nFSR = nx*ny
ENDIF
IF(.NOT. lRECT) THEN
  nFSR = BaseCellInfo(icel0)%nFSR/2
  IF(iType .EQ. 3) nFsr = nFsr/2
ENDIF

!ALLOCATION PART
CALL dmalloc(BaseCellInfo(icel)%FxrIdxSt, nFxr)
CALL dmalloc(BaseCellInfo(icel)%nFsrInFxr, nFxr)
CALL dmalloc(BaseCellInfo(icel)%MapFxr2FsrIdx, nFSR, nFxr)

IF(lRect) THEN
  CALL Dmalloc(BaseCellInfo(icel)%Geom%delx, nx)
  CALL Dmalloc(BaseCellInfo(icel)%Geom%dely, ny)
ENDIF
HalfPitch = Half*CellPitch
lx = (/CellPitch, HalfPitch, HalfPitch/)
ly = (/HalfPitch, CellPitch, HalfPitch/)
lx = (/ONE,  HALF, HALF/)
ly = (/HALF, ONE,  HALF/)


!Center Line Pincell Setting
BaseCellInfo(icel)%nFxr = nFxr;  BaseCellInfo(icel)%nFSR = nFSR
BaseCellInfo(icel)%Geom%nbd = 4
BaseCellInfo(icel)%Geom%lx = lx(itype)*BaseCellInfo(icel0)%Geom%lx
BaseCellInfo(icel)%Geom%ly = ly(itype)*BaseCellInfo(icel0)%Geom%ly

BaseCellInfo(icel)%Geom%cx = BaseCellInfo(icel0)%Geom%cx
BaseCellInfo(icel)%Geom%cy = BaseCellInfo(icel0)%Geom%cy

SELECT CASE(iType)
  CASE(1)
    BaseCellInfo(icel)%nDivAzi = BaseCellInfo(icel0)%nDivAzi/2
    BaseCellInfo(icel)%Geom%x(1) = -HalfPitch; BaseCellInfo(icel)%Geom%x(2) = HalfPitch
    BaseCellInfo(icel)%Geom%y(1) = -HalfPitch; BaseCellInfo(icel)%Geom%y(2) = 0
!     ny = nint(0.5*ny)
  CASE(2)
    BaseCellInfo(icel)%nDivAzi = BaseCellInfo(icel0)%nDivAzi/2
    BaseCellInfo(icel)%Geom%x(1) = 0; BaseCellInfo(icel)%Geom%x(2) = HalfPitch
    BaseCellInfo(icel)%Geom%y(1) = -HalfPitch; BaseCellInfo(icel)%Geom%y(2) = HalfPitch
  CASE(3)
    BaseCellInfo(icel)%nDivAzi = BaseCellInfo(icel0)%nDivAzi/4
    BaseCellInfo(icel)%Geom%x(1) = 0; BaseCellInfo(icel)%Geom%x(2) = HalfPitch
    BaseCellInfo(icel)%Geom%y(1) = -HalfPitch; BaseCellInfo(icel)%Geom%y(2) = 0
    !BaseCellInfo(icel)%Geom%x(1) = 0; BaseCellInfo(icel)%Geom%x(2) = HalfPitch
    !BaseCellInfo(icel)%Geom%y(1) = 0; BaseCellInfo(icel)%Geom%y(2) = HalfPitch
END SELECT

IF(BaseCellInfo(icel)%lGap) THEN
  lx = (/ONE,  HALF, HALF/)
  ly = (/HALF, ONE,  HALF/)
  BaseCellInfo(icel)%Geom%lx = lx(itype)*BaseCellInfo(icel0)%Geom%lx
  BaseCellInfo(icel)%Geom%ly = ly(itype)*BaseCellInfo(icel0)%Geom%ly
  BaseCellInfo(icel)%Geom%x = BaseCellInfo(icel0)%Geom%x
  BaseCellInfo(icel)%Geom%y = BaseCellInfo(icel0)%Geom%y

  SELECT CASE(iType)
    CASE(1)
      !BaseCellInfo(icel)%Geom%x(1) = -HalfPitch; BaseCellInfo(icel)%Geom%x(2) = HalfPitch
      BaseCellInfo(icel)%Geom%y(1) = -BaseCellInfo(icel)%Geom%ly; BaseCellInfo(icel)%Geom%y(2) = 0
    CASE(2)
      BaseCellInfo(icel)%Geom%x(1) = 0; BaseCellInfo(icel)%Geom%x(2) = BaseCellInfo(icel)%Geom%lx
    CASE(3)
      BaseCellInfo(icel)%Geom%x(1) = 0; BaseCellInfo(icel)%Geom%x(2) = BaseCellInfo(icel)%Geom%lx
      BaseCellInfo(icel)%Geom%y(1) = -BaseCellInfo(icel)%Geom%ly; BaseCellInfo(icel)%Geom%y(2) = 0
  END SELECT
ENDIF

BaseCellInfo(icel)%lRect = lRECT;
BaseCellInfo(icel)%Geom%lRect = lRect; BaseCellInfo(icel)%Geom%lCircle = .NOT. lRect
BaseCellInfo(icel)%Geom%nCircle = 0
IF(.NOT. lRect) THEN
  nFsrInFxr = BaseCellInfo(icel)%nDivAzi
  !Geometries Setting
  BaseCellInfo(iCel)%Geom%ncircle = BaseCellInfo(icel0)%Geom%nCircle
  nCircle = BaseCellInfo(iCel)%Geom%ncircle
  BaseCellInfo(iCel)%Geom%cx = BaseCellInfo(icel0)%Geom%cx
  BaseCellInfo(iCel)%Geom%cy = BaseCellInfo(icel0)%Geom%cy
  Call Dmalloc(BaseCellInfo(icel)%Geom%circle, 3, nCircle)
  BaseCellInfo(icel)%Geom%circle(3, :) = BaseCellInfo(icel0)%Geom%circle(3, :)
  !Index Setting - FSR, FXR index setting
  ir = 0
  DO i = 1, nFXR
    ir0 = BaseCellInfo(icel0)%MapFxr2FsrIdx(1,i)
    BaseCellInfo(icel)%nFsrInFxr(i) = nFsrInFxr
    BaseCellInfo(icel)%FxrIdxSt(i) = nFsrInFxr*(i-1)+1
    DO j = 1,nFsrInFxr
      BaseCellInfo(icel)%MapFxr2FsrIdx(j,i) = BaseCellInfo(icel)%FxrIdxSt(i) + j - 1
    ENDDO
  ENDDO
  IF (nTracerCntl%lSSPH) then
      IF (BaseCellInfo(icel)%lfuel) then
          ALLOCATE(BaseCellInfo(icel)%SPHfactor(1:BaseCellInfo(icel0)%nFXR,igresb:igrese))
          BaseCellInfo(icel)%SPHfactor = BaseCellInfo(icel0)%SPHfactor
          BaseCellInfo(icel)%nfueldiv = BaseCellInfo(icel0)%nfueldiv
          BaseCellInfo(icel)%ngapdiv = BaseCellInfo(icel0)%ngapdiv
          BaseCellInfo(icel)%ncladdiv = BaseCellInfo(icel0)%ncladdiv
          BaseCellInfo(icel)%FuelRad0 = BaseCellInfo(icel0)%FuelRad0
          BaseCellInfo(icel)%srdidx = BaseCellInfo(icel0)%srdidx
      ENDIF
  ENDIF
ELSE
  nFsrInFxr = nFSR
  BaseCellInfo(icel)%Geom%nx = nx
  BaseCellInfo(icel)%Geom%ny = ny
  ir = 0
  DO j = BaseCellInfo(icel0)%Geom%ny - ny + 1 , BaseCellInfo(icel0)%Geom%ny
    DO i = BaseCellInfo(icel0)%Geom%nx - nx + 1 , BaseCellInfo(icel0)%Geom%nx
      k = BaseCellInfo(icel0)%Geom%nx * (j - 1) + i; ir = ir + 1
    ENDDO
  ENDDO

  DO m = 1, nFxr
    ir = 0; j = 0; k =0
    DO i = 1, BaseCellInfo(icel0)%nFsr
      IF(FxrMap(i) .LT. 0) CYCLE
      j = j + 1
      IF(FxrMap(i) .NE. m) CYCLE
      !BaseCellInfo(icel)%nFsrInFxr(m) = BaseCellInfo(icel)%nFsrInFxr(m) + 1
      k = k + 1; BaseCellInfo(icel)%nFsrInFxr(m) = k
      BaseCellInfo(icel)%MapFxr2FsrIdx(k, m) = j
    ENDDO
    BaseCellInfo(icel)%FxrIdxSt(m) = BaseCellInfo(icel)%MapFxr2FsrIdx(1, m)
  ENDDO

  BaseCellInfo(icel)%Geom%Delx(1:nx) = DelX(1:nx)
  BaseCellInfo(icel)%Geom%Dely(1:ny) = DelY(1:ny)
  CONTINUE
ENDIF

END SUBROUTINE

SUBROUTINE SetCentPin(CentPinInfo, PinInfo, CellInfo, iType, ipin)
!USE geom,     ONLY :  PinInfo,      CellInfo
USE Typedef, only : pininfo_type, Cell_type
USE Geom, only : nCellType
IMPLICIT NONE
TYPE(PinInfo_type) :: CentPinInfo, PinInfo
TYPE(Cell_type), pointer :: CellInfo(:)
INTEGER :: iPin, iType
INTEGER :: nCell, iCel0, icel
INTEGER :: nFsrMax, nFxrMax
INTEGER :: i, j, k
CentPinInfo%lEmpty = PinInfo%lEmpty
IF(CentPinInfo%lEmpty) RETURN

PinInfo%EdgePinIdx(iType) = ipin
CentPinInfo%luse = FALSE;   CentPinInfo%lFuel = PinInfo%lFuel

CentPinInfo%nCell = PinInfo%nCell;  nCell = CentPinInfo%nCell
CALL Dmalloc(CentPinInfo%Cell, nCell)
DO i = 1, nCell
  icel0 = PinInfo%cell(i);    icel = CellInfo(icel0)%EdgeCellIdx(iType)
  CentPinInfo%Cell(i) = icel
ENDDO

nFsrMax = 0; nFxrMax = 0
DO i = 1, nCell
  icel = CentPinInfo%Cell(i)
  nFsrMax = max(nFsrMax, CellInfo(icel)%nFsr)
  nFxrMax = max(nFxrMax, CellInfo(icel)%nFxr)
ENDDO

CentPinInfo%lCentX = CellInfo(icel)%lCentX;
CentPinInfo%lCentY = CellInfo(icel)%lCentY;
CentPinInfo%lCentXY = CellInfo(icel)%lCentXY;

CentPinInfo%nFsrMax = nFsrMax
CentPinInfo%nFxrMax = nFxrMax

ENDSUBROUTINE

SUBROUTINE SetCentAsy(CentAsyInfo, AsyInfo, PinInfo, iType, iasy)
USE Typedef,      only : AsyInfo_type,    pininfo_type
USE Geom,     only : nPinType
IMPLICIT NONE

TYPE(AsyInfo_type) :: CentAsyInfo, AsyInfo
TYPE(pininfo_type), pointer :: PinInfo(:)
INTEGER, INTENT(IN) :: iType, iasy

INTEGER :: i, j, k
INTEGER :: ixBeg(3), iyBeg(3), ixbeg0, iybeg0
INTEGER :: ipin, ipin0, ixy, ix, iy, nx, ny, nxy
INTEGER :: nHalfCellX, nCellX
INTEGER :: AsyConfig(200,200)
LOGICAL :: CentXYINFO(3,3)
DATA CentXYINFO /TRUE,  FALSE, FALSE, &
                 FALSE, TRUE,  FALSE, &
                 FALSE, FALSE, TRUE/
CentAsyInfo%lEmpty = AsyInfo%lEmpty
IF(CentAsyInfo%lEmpty) RETURN

AsyInfo%EdgeAsyIdx(iType) = iasy
CentAsyInfo%lFuel = AsyInfo%lFuel; CentAsyInfo%luse = FALSE
!CentAsyInfo%lCentX = AsyInfo%lCentX; CentAsyInfo%lCentY = AsyInfo%lCentY;
!CentAsyInfo%lCentXY = AsyInfo%lCentXY;
CentAsyInfo%lCentX = CentXYInfo(1, itype)
CentAsyInfo%lCentY = CentXYInfo(2, itype)
CentAsyInfo%lCentXY = CentXYInfo(3, itype)
!type AsyInfo_type
!  logical :: lempty,lfuel,lgeom,lEdgeX,lEdgeY,lEdgeXY
!  INTEGER :: nx,ny,nxy
!  INTEGER,pointer :: pin(:),pin2DIdx(:,:)
!  INTEGER :: EdgeAsyIdx(3)  !1: X-dir edge, 2:Y-dir edge, 3:both
!END type



nCellX = AsyInfo%nx
nHalfCellX = nCellX/2 + mod(NcellX, 2)
ixBeg = (/1, nHalfCellx, nHalfCellx/); iyBeg = (/nHalfCellx, 1, nHalfCellx/)

nx = nCellx - ixBeg(iType) + 1;
ny = nCellx - iyBeg(iType) + 1;
IF(nx .lt. nCellx .and. mod(nCellx, 2) .EQ. 0) nx = nx - 1
IF(ny .lt. nCellx .and. mod(nCellx, 2) .EQ. 0) ny = ny - 1
nxy = nx*ny
CentAsyInfo%nx = nx; CentAsyInfo%ny =ny; CentAsyInfo%nxy = nxy
CALL Dmalloc(CentAsyInfo%pin, nxy); CALL Dmalloc(CentAsyInfo%pin2DIdx, nx, ny)

iybeg0 = iybeg(itype); ixbeg0 = ixbeg(itype)
IF(mod(nCellX, 2) .EQ. 0) THEN
  IF(iType .EQ. 1) THEN
    iybeg0 = iybeg0 + 1
  ELSEIF(iType .EQ. 2) THEN
    ixbeg0 = ixbeg0 + 1
  ELSEIF(iType .EQ. 3) THEN
    iybeg0 = iybeg0 + 1; ixbeg0 = ixbeg0 + 1
  ENDIF
ENDIF
i=0
DO iy = iyBeg0, nCellX
  DO ix = ixBeg0, nCellX
    j = AsyInfo%pin2DIdx(ix,iy); ipin = AsyInfo%pin(j)
    AsyConfig(ix - ixBeg0 + 1, iy - iyBeg0 + 1) = ipin
  ENDDO
ENDDO

IF(mod(nCellX,2) .eq. 1) then
  IF(iType .EQ. 1) THEN
    CentAsyInfo%lCentX = TRUE; iy = 1
    DO ix = 1, nCellX
      ipin0 = AsyConfig(ix, iy); ipin = PinInfo(ipin0)%EdgePinIdx(iType)
      AsyConfig(ix, iy) = ipin
    ENDDO
  ENDIF

  IF(iType .EQ. 2) THEN
    CentAsyInfo%lCentY = TRUE; ix = 1
    DO iy = 1, nCellX
      ipin0 = AsyConfig(ix, iy); ipin = PinInfo(ipin0)%EdgePinIdx(iType)
      AsyConfig(ix, iy) = ipin
    ENDDO
  ENDIF

  IF(iType .EQ. 3) THEN
    CentAsyInfo%lCentXY = TRUE; ix = 1
    DO iy = 2, ny
      ipin0 = AsyConfig(ix, iy); ipin = PinInfo(ipin0)%EdgePinIdx(2)
      AsyConfig(ix, iy) = ipin
    ENDDO

    iy = 1
    DO ix = 2, nx
      ipin0 = AsyConfig(ix, iy); ipin = PinInfo(ipin0)%EdgePinIdx(1)
      AsyConfig(ix, iy) = ipin
    ENDDO

    ix = 1; iy = 1
    ipin0 = AsyConfig(ix, iy); ipin = PinInfo(ipin0)%EdgePinIdx(3)
    AsyConfig(ix, iy) = ipin
  ENDIF
ENDIF
i = 0
DO iy = 1, ny
  DO ix = 1, nx
    i = i + 1
    CentAsyInfo%Pin(i) = AsyCOnfig(ix, iy)
    CentAsyInfo%Pin2DIdx(ix, iy) = i
  ENDDO
ENDDO
CentAsyInfo%nFuelPin = 0
DO i = 1, nxy
  ipin = CentASyInfo%Pin(i)
  if(ipin.NE.0)then
    IF(PinInfo(ipin)%lFuel) CentAsyInfo%nFuelPin = CentAsyInfo%nFuelPin + 1
  endif
ENDDO
END SUBROUTINE

SUBROUTINE SetCentCore(CoreInfo, AsyInfo)
USE Typedef, only : CoreInfo_Type, AsyInfo_type
USE Geom, only : nAsyType, nAsyType0
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(AsyInfo_type),Pointer :: AsyInfo(:)

INTEGER :: iasy0, iasy, itype
INTEGER :: ixa, iya, ixya
INTEGER :: iasy1, iasy2, iasy3
INTEGER :: nxy, nxa, nya

nxy = CoreInfo%nxy; nxa = CoreInfo%nxa; nya = CoreInfo%nya

IF(CoreInfo%RadSym(NORTH)) THEN
  iya = 1; iType = 1
  DO ixa = 1, nxa
    ixya = CoreInfo%CoreIdx(ixa, iya)
    IF(ixya .eq. 0) CYCLE
    iasy0 = CoreInfo%CoreMap(ixya)
    iasy = AsyInfo(iasy0)%EdgeAsyIdx(iType)
    CoreInfo%CoreMap(ixya) = iasy
  ENDDO
ENDIF

IF(CoreInfo%RadSym(WEST)) THEN
  ixa = 1; iType = 2
  DO iya = 1, nya
    ixya = CoreInfo%CoreIdx(ixa, iya)
    IF(ixya .eq. 0) CYCLE
    iasy0 = CoreInfo%CoreMap(ixya)
    iasy = AsyInfo(iasy0)%EdgeAsyIdx(iType)
    CoreInfo%CoreMap(ixya) = iasy
  ENDDO
ENDIF

IF(CoreInfo%RadSym(NORTH) .AND. CoreInfo%RadSym(WEST)) THEN
  ixa = 1; iya = 1; iType = 3
  ixya = CoreInfo%CoreIdx(ixa, iya); iasy0 = CoreInfo%CoreMap(ixya)
  iasy = AsyInfo(iasy0)%EdgeAsyIdx(iType)
  CoreInfo%CoreMap(ixya) = iasy
ENDIF

END SUBROUTINE
end module
