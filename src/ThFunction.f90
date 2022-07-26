#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE FuelTempChange(Core, Fxr,  deltemp, nTracerCntl, PE)
!Calculate the reference fuel temperature and save the reference temperature
USE PARAM
USE TYPEDEF,       ONLY : coreinfo_type,     Fxrinfo_type,       PE_TYPE,   &
                          Pin_Type,          Cell_Type    
USE CNTL,          ONLY : nTracerCntl_Type
#ifdef MPI_ENV
USE MPICOMM_mod,   ONLY : REDUCE,      MPI_SYNC
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(Fxrinfo_type), POINTER :: Fxr(:, :)
TYPE(PE_TYPE) :: PE
TYPE(nTracerCntl_Type) :: nTracerCntl
REAL :: DelTemp

TYPE(Pin_Type), POINTER :: PIN(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
REAL, POINTER :: hz(:)
REAL :: sumT, sumVol, sumPlnT, sumPlnVol
REAL :: h, vol, temp

INTEGER :: j
INTEGER :: nxy, nCoreFsr, nCoreFxr, nFsrInFxr, FsrIdxSt, FxrIdxSt, nlocalFxr
INTEGER :: myzb, myze
INTEGER :: iz, ixy, icel, ifxr,ifsr
#ifdef MPI_ENV
REAL :: MPIdat(2), MPIdatbuf(2)
#endif
IF(nTracerCntl%lBenchXs) RETURN

hz => Core%hz; Pin => Core%Pin
CellInfo => Core%CellInfo;

nxy = Core%nxy; nCoreFsr = Core%nCoreFsr; nCoreFxr = Core%nCoreFxr
myzb = PE%myzb; myze = PE%myze

SumT=0._8; SumVol = 0._8;
DO iz = myzb, myze
  h = hz(iz)
  sumPlnT = 0.; sumPlnVol = 0
  DO ixy = 1, nxy
    FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
    icel = Pin(ixy)%Cell(iz)
    nlocalFxr = CellInfo(icel)%nFxr
    DO j = 1, nlocalFxr
      ifxr = FxrIdxSt + j -1
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
      IF(.NOT. Fxr(ifxr, iz)%lfuel) CYCLE
      Fxr(ifxr, iz)%temp = Fxr(ifxr, iz)%temp + deltemp
    ENDDO
  ENDDO
ENDDO


NULLIFY(hz); NULLIFY(Pin); NULLIFY(CellInfo)

END SUBROUTINE
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE Cal_RefFuelTemp(RefFuelTemp, Core, Fxr, nTracerCntl, PE)

USE PARAM
USE TYPEDEF,       ONLY : coreinfo_type, Fxrinfo_type, PE_TYPE, Pin_Type, Cell_Type
USE CNTL,          ONLY : nTracerCntl_Type

#ifdef MPI_ENV
USE MPICOMM_mod,   ONLY : REDUCE, MPI_SYNC
#endif

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(Fxrinfo_type), POINTER :: Fxr(:, :)
TYPE(PE_TYPE) :: PE
TYPE(nTracerCntl_Type) :: nTracerCntl
REAL, POINTER :: RefFuelTemp(:)
! ----------------------------------------------------
TYPE(Pin_Type), POINTER :: PIN(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
REAL, POINTER :: hz(:)
REAL :: sumT, sumVol, sumPlnT, sumPlnVol, h, vol, temp
INTEGER :: nxy, nCoreFsr, nCoreFxr, nFsrInFxr, FsrIdxSt, FxrIdxSt, nlocalFxr, myzb, myze
INTEGER :: j, iz, ixy, icel, ifxr,ifsr

#ifdef MPI_ENV
REAL :: MPIdat(2), MPIdatbuf(2)
#endif
! ----------------------------------------------------

IF (nTracerCntl%lBenchXs) RETURN

hz       => Core%hz
Pin      => Core%Pin
CellInfo => Core%CellInfo;
nxy       = Core%nxy
nCoreFsr  = Core%nCoreFsr
nCoreFxr  = Core%nCoreFxr

myzb = PE%myzb
myze = PE%myze

SumT   = ZERO
SumVol = ZERO

DO iz = myzb, myze
  h         = hz(iz)
  sumPlnT   = ZERO
  sumPlnVol = ZERO
  
  DO ixy = 1, nxy
    FsrIdxSt = Pin(ixy)%FsrIdxSt
    FxrIdxSt = Pin(ixy)%FxrIdxSt
    icel     = Pin(ixy)%Cell(iz)
    
    nlocalFxr = CellInfo(icel)%nFxr
    
    DO j = 1, nlocalFxr
      ifxr      = FxrIdxSt + j -1
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
      
      IF (.NOT. Fxr(ifxr, iz)%lfuel) CYCLE
      
      vol  = FXR(ifxr, iz)%area * h
      temp = Fxr(ifxr, iz)%temp
      
      sumPlnT   = sumPlnT   + temp*vol
      SumPlnVol = SumPlnVol + vol
    END DO
  END DO
  
  sumT   = sumT   + sumPlnT
  SumVol = SumVol + SumPlnVol
  
  IF (SumPlnVOl .GT. ZERO) THEN
    RefFuelTemp(iz) = sumPlnT / SumPlnVol - CKELVIN
  ELSE
    RefFuelTemp(iz) = ZERO
  END IF
END DO

#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_RTMASTER_COMM)
MPIdat = (/sumT, SumVol/)
CALL REDUCE(MPIdat, MPIdatbuf, 2, PE%MPI_RTMASTER_COMM, TRUE)
sumT = MPIdatBuf(1); SumVol = MPIdatBuf(2)
#endif

RefFuelTemp(0) = SumT / SumVol - CKELVIN

DO iz = myzb, myze
  IF (RefFuelTemp(iz) .LT. 1E-5) RefFuelTemp(iz) = RefFuelTemp(0)
END DO

NULLIFY (hz)
NULLIFY (Pin)
NULLIFY (CellInfo)
! ----------------------------------------------------

END SUBROUTINE Cal_RefFuelTemp
! ------------------------------------------------------------------------------------------------------------
FUNCTION CalPlnFuelTemp(Core, Fxr, ipln)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type       ,FxrInfo_Type                &
                            ,Pin_Type            ,Cell_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:)
INTEGER :: ipln

REAL :: CalPlnFuelTemp

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
INTEGER :: nLocalFxr, FxrIdxSt, nxy
INTEGER :: ixy, icel, ifxr
INTEGER :: i, j, k

REAL :: area, temp, areasum, tempsum

nxy = Core%nxy;
Pin => Core%Pin; CellInfo => Core%CellInfo

CalPlnFuelTemp = 0
areasum = 0; tempsum = 0
DO ixy = 1, nxy
  IF(.NOT. Pin(ixy)%lFuel) CYCLE
  FxrIdxSt = Pin(ixy)%FxrIdxSt
  icel = Pin(ixy)%Cell(ipln)
  nlocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nlocalFxr
    ifxr = FxrIdxSt + j - 1
    IF(.NOT. Fxr(ifxr)%lfuel) CYCLE
    area = Fxr(ifxr)%area; temp = Fxr(ifxr)%temp
    areasum = areasum + area
    tempsum = tempsum + temp * area
  ENDDO
ENDDO
tempsum = tempsum / areasum
CalPlnFuelTemp = tempsum

NULLIFY(Pin); NULLIFY(CellInfo)
END FUNCTION
! ------------------------------------------------------------------------------------------------------------
FUNCTION CalPlnModTemp(Core, Fxr, ipln)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type       ,FxrInfo_Type                &
                            ,Pin_Type            ,Cell_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:)
INTEGER :: ipln

REAL :: CalPlnModTemp

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
INTEGER :: nLocalFxr, FxrIdxSt, nxy
INTEGER :: ixy, icel, ifxr
INTEGER :: i, j, k

REAL :: area, temp, areasum, tempsum

nxy = Core%nxy;
Pin => Core%Pin; CellInfo => Core%CellInfo

CalPlnModTemp = 0
areasum = 0; tempsum = 0
DO ixy = 1, nxy
  IF(.NOT. Pin(ixy)%lFuel) CYCLE
  FxrIdxSt = Pin(ixy)%FxrIdxSt
  icel = Pin(ixy)%Cell(ipln)
  nlocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nlocalFxr
    ifxr = FxrIdxSt + j - 1
    IF(.NOT. Fxr(ifxr)%lh2o) CYCLE
    area = Fxr(ifxr)%area; temp = Fxr(ifxr)%temp
    areasum = areasum + area
    tempsum = tempsum + temp * area
  ENDDO
ENDDO
tempsum = tempsum / areasum
CalPlnModTemp = tempsum

NULLIFY(Pin); NULLIFY(CellInfo)
END FUNCTION
! ------------------------------------------------------------------------------------------------------------
FUNCTION CalPlnTemp(Core, Fxr, ipln)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type       ,FxrInfo_Type                &
                            ,Pin_Type            ,Cell_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:)
INTEGER :: ipln

REAL :: CalPlnTemp

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
INTEGER :: nLocalFxr, FxrIdxSt, nxy
INTEGER :: ixy, icel, ifxr
INTEGER :: i, j, k

REAL :: area, temp, areasum, tempsum

nxy = Core%nxy;
Pin => Core%Pin; CellInfo => Core%CellInfo

CalPlnTemp = 0
areasum = 0; tempsum = 0
DO ixy = 1, nxy
  !IF(.NOT. Pin(ixy)%lFuel) CYCLE
  FxrIdxSt = Pin(ixy)%FxrIdxSt
  icel = Pin(ixy)%Cell(ipln)
  nlocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nlocalFxr
    ifxr = FxrIdxSt + j - 1
    area = Fxr(ifxr)%area; temp = Fxr(ifxr)%temp
    areasum = areasum + area
    tempsum = tempsum + temp * area
  ENDDO
ENDDO
tempsum = tempsum / areasum
CalPlnTemp = tempsum

NULLIFY(Pin); NULLIFY(CellInfo)
END FUNCTION
! ------------------------------------------------------------------------------------------------------------
FUNCTION GetPinFuelTemp(Core, Fxr, ipln, ipin)
USE PARAM
USE TYPEDEF,       ONLY : CoreInfo_Type       ,FxrInfo_Type            &
                         ,Pin_Type            ,Cell_Type
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
INTEGER :: ipln, ipin
REAL :: GetPinFuelTemp

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
INTEGER :: nLocalFxr, FxrIdxSt, nxy
INTEGER :: ixy, icel, ifxr
INTEGER :: i, j, k

nxy = Core%nxy;
Pin => Core%Pin; CellInfo => Core%CellInfo

  !IF(.NOT. Pin(ixy)%lFuel) CYCLE
FxrIdxSt = Pin(ipin)%FxrIdxSt
icel = Pin(ipin)%Cell(ipln)
nlocalFxr = CellInfo(icel)%nFxr
TempSum = 0; areasum = 0;
DO j = 1, nlocalFxr
  ifxr = FxrIdxSt + j - 1
  !IF(.NOT. Fxr(ifxr, ipln)%lFuel) CYCLE
  IF(CellInfo(icel)%lfuel .AND. .NOT. Fxr(ifxr, ipln)%lFuel) CYCLE
  IF(.NOT. CellInfo(icel)%lfuel .AND. .NOT. Fxr(ifxr, ipln)%lres) CYCLE
  area = Fxr(ifxr, ipln)%area; temp = Fxr(ifxr, ipln)%temp
  areasum = areasum + area
  tempsum = tempsum + temp * area
ENDDO
IF( areasum .GT. 0 )THEN
GetPinFuelTemp = tempsum / areasum
ELSE
    GetPinFuelTemp = 0._8
ENDIF

NULLIFY(Pin); NULLIFY(CellInfo)
END FUNCTION
! ------------------------------------------------------------------------------------------------------------
FUNCTION GetPinModTemp(Core, Fxr, ipln, ipin)
USE PARAM
USE TYPEDEF,       ONLY : CoreInfo_Type       ,FxrInfo_Type            &
                         ,Pin_Type            ,Cell_Type
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
INTEGER :: ipln, ipin
REAL :: GetPinModTemp

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
INTEGER :: nLocalFxr, FxrIdxSt, nxy
INTEGER :: ixy, icel, ifxr
INTEGER :: i, j, k

nxy = Core%nxy;
Pin => Core%Pin; CellInfo => Core%CellInfo
FxrIdxSt = Pin(ipin)%FxrIdxSt
icel = Pin(ipin)%Cell(ipln)
nlocalFxr = CellInfo(icel)%nFxr
TempSum = 0; areasum = 0;
DO j = 1, nlocalFxr
  ifxr = FxrIdxSt + j - 1
  IF(.NOT. Fxr(ifxr, ipln)%lh2o) CYCLE
  area = Fxr(ifxr, ipln)%area; temp = Fxr(ifxr, ipln)%temp
  areasum = areasum + area
  tempsum = tempsum + temp * area
ENDDO
IF( areasum .GT. 0 )THEN
GetPinModTemp = tempsum / areasum
ELSE
    GetPinModTemp = 0._8
ENDIF
NULLIFY(Pin); NULLIFY(CellInfo)    
END FUNCTION
! ------------------------------------------------------------------------------------------------------------
FUNCTION GetPinTemp(Core, Fxr, ipln, ipin)
USE PARAM
USE TYPEDEF,       ONLY : CoreInfo_Type       ,FxrInfo_Type            &
                         ,Pin_Type            ,Cell_Type
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
INTEGER :: ipln, ipin
REAL :: GetPinTemp

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
INTEGER :: nLocalFxr, FxrIdxSt, nxy
INTEGER :: ixy, icel, ifxr
INTEGER :: i, j, k

nxy = Core%nxy;
Pin => Core%Pin; CellInfo => Core%CellInfo
FxrIdxSt = Pin(ipin)%FxrIdxSt
icel = Pin(ipin)%Cell(ipln)
nlocalFxr = CellInfo(icel)%nFxr
TempSum = 0; areasum = 0;
DO j = 1, nlocalFxr
  ifxr = FxrIdxSt + j - 1
  area = Fxr(ifxr, ipln)%area; temp = Fxr(ifxr, ipln)%temp
  areasum = areasum + area
  tempsum = tempsum + temp * area
ENDDO
IF( areasum .GT. 0 )THEN
GetPinTemp = tempsum / areasum
ELSE
    GetPinTemp = 0._8
ENDIF

NULLIFY(Pin); NULLIFY(CellInfo)
END FUNCTION
! ------------------------------------------------------------------------------------------------------------