#include <defines.h>
SUBROUTINE PrintDepletion(io, Core, FmInfo, CmInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,        FmInfo_Type,       CmInfo_Type,     &
                             PE_Type
USE CNTL,             ONLY : nTRACERCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
INTEGER :: io

INTEGER :: id, idiso

IF(.NOT. nTracerCntl%lXsLib) RETURN

CALL PrintAsyBurnUp(io, Core, FmInfo, CmInfo, nTracerCntl, PE)
IF(nTracerCntl%OutpCntl%IsoOutList(0) .EQ. 0) RETURN
IF(PE%MASTER) WRITE(io, '(A)') ' - Assembly-wise Isotopic Inventory -'
DO id = 1, nTracerCntl%OutpCntl%IsoOutList(0)
  idiso = nTracerCntl%outpCntl%IsoOutList(id)
  CALL PrintAsyND(io, IdIso, Core, FmInfo, CmInfo, nTracerCntl, PE)
ENDDO
CALL PrintAsyND2(Core, FmInfo, CmInfo, nTracerCntl, PE)
END SUBROUTINE

SUBROUTINE PrintAsyBurnUp(io, Core, FmInfo, CmInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,        FmInfo_Type,       CmInfo_Type,     &
                             PE_Type,                                                  &
                             FxrInfo_Type,         Pin_Type,          Cell_Type,       &
                             AsyInfo_Type,         Asy_Type
USE CNTL,             ONLY : nTRACERCntl_Type
USE BasicOperation,   ONLY : CP_CA
USE IOUTIL,           ONLY : PrintReal1DarrayTo2Darray
USE MPIComm_Mod,      ONLY : REDUCE
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
INTEGER :: io

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(AsyInfo_Type), POINTER :: AsyType(:)
TYPE(Asy_Type), POINTER :: Asy(:)

INTEGER,  POINTER :: CoreIdx(:, :), CoreMap(:)
REAL, POINTER :: hz(:)
REAL, POINTER :: AsyBurnup(:)


REAL :: Hmsum, burnupsum, maxv
REAL :: buf0(2), buf(2)

INTEGER :: myzb, myze, nfxr, nfsr, nxy, nxya, nxya0, nxa, nya
INTEGER :: nyc, nxc, nxy0, nLocalFxr, FxrIdxSt
INTEGER :: iz, ifxr, ifsr, icel, ipin, ixya, ixa, iya, iasy, iasytype
INTEGER :: iy, ix, ixy, i, j, k
CHARACTER(120) :: formatdata

Fxr => FmInfo%Fxr; Pin => Core%Pin; CellInfo => Core%CellInfo
AsyType => Core%AsyInfo; Asy => Core%Asy
CoreIdx => Core%CoreIdx; CoreMap => Core%CoreMap
hz => Core%hz
nxy = Core%nxy
nxa = Core%nxa; nya = Core%nya; nxya = Core%nxya
nxya0 = nxa * nya
nfsr = Core%nCoreFsr; nFxr = Core%nCoreFxr
myzb = PE%myzb; myze = PE%myze

ALLOCATE(AsyBurnUp(nxya0))
CALL CP_CA(AsyBurnUp, 0._8, nxya0)

i = 0
DO iya = 1, Core%nya
  DO ixa = 1, Core%nxa
    i = i + 1
    iasy = CoreIdx(ixa, iya) !Assembly Idx
    IF(iasy .EQ. 0) CYCLE
    iasytype = CoreMap(iasy)
    IF(.NOT. AsyType(iasytype)%lFuel) CYCLE
    !Assembly Sweeping
    nxy0 = AsyType(iasytype)%nxy
    Hmsum = 0; burnupsum = 0
    DO iz = myzb, myze
      DO ixy = 1, nxy0
        ipin = Asy(iasy)%GlobalPinIdx(ixy)
        icel = Pin(ipin)%Cell(iz)
        FxrIdxSt = Pin(ipin)%FxrIdxSt
        nLocalFxr = CellInfo(icel)%nFxr
        DO j = 1, nLocalFxr
          ifxr = FxrIdxSt + j -1
          IF(.NOT. Fxr(ifxr, iz)%lDepl) CYCLE
!          volsum = volsum + Fxr(ifxr, iz)%area * hz(iz)
!          burnupsum = burnupsum + Fxr(ifxr, iz)%burnup * Fxr(ifxr, iz)%area * hz(iz)
          Hmsum = Hmsum + Fxr(ifxr, iz)%Hmkg0
          burnupsum = burnupsum + Fxr(ifxr, iz)%burnup * Fxr(ifxr, iz)%Hmkg0
        ENDDO
      ENDDO
    ENDDO
#ifdef MPI_ENV
    buf0 = (/Hmsum, burnupsum/)
    CALL REDUCE(buf0, buf, 2, PE%MPI_CMFD_COMM, .TRUE.)
    Hmsum = buf(1); burnupsum = buf(2)
#endif
    AsyBurnUp(i) = burnupsum/Hmsum
  ENDDO
ENDDO
IF(PE%MASTER) THEN
WRITE(io, '(A)') ' - Assembly-wise Burnup(MWD/kgU) -'
WRITE(formatdata,'(a)') '(8x, 200F8.4)'
WRITE(formatdata,'(a)') '(8x, 200F8.4)'
CALL PrintReal1DarrayTo2Darray(io,  AsyBurnUp(1:nxya0), nxya0, nxa, nya, formatdata)
maxv = maxval(AsyBurnUp(1:nxya0))
WRITE(io, *)
WRITE(io, '(5x, A7, 2x, F7.3)') 'Max. = ', maxv
WRITE(io, *)
ENDIF
DEALLOCATE(AsyBurnUp)
END SUBROUTINE

SUBROUTINE PrintAsyND(io, IdIso, Core, FmInfo, CmInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,        FmInfo_Type,       CmInfo_Type,     &
                             PE_Type,                                                  &
                             FxrInfo_Type,         Pin_Type,          Cell_Type,       &
                             AsyInfo_Type,         Asy_Type
USE CNTL,             ONLY : nTRACERCntl_Type
USE BasicOperation,   ONLY : CP_CA
USE IOUTIL,           ONLY : PrintReal1DarrayTo2Darray
USE MPIComm_Mod,      ONLY : REDUCE
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
INTEGER :: io, IdIso

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(AsyInfo_Type), POINTER :: AsyType(:)
TYPE(Asy_Type), POINTER :: Asy(:)

INTEGER,  POINTER :: CoreIdx(:, :), CoreMap(:)
REAL, POINTER :: hz(:)
REAL, POINTER :: AsyND(:), AsyND_Fuel(:)


REAL :: volsum, NDsum, FuelVolSum, maxv
REAL :: buf0(3), buf(3)

INTEGER :: myzb, myze, nfxr, nfsr, nxy, nxya, nxya0, nxa, nya
INTEGER :: nyc, nxc, nxy0, nLocalFxr, FxrIdxSt
INTEGER :: iz, ifxr, ifsr, icel, ipin, ixya, ixa, iya, iasy, iasytype
INTEGER :: iy, ix, ixy, i, j, k, iso, id
CHARACTER(120) :: formatdata
LOGICAL :: lExist
IF(.NOT. nTracerCntl%lXsLib) RETURN

Fxr => FmInfo%Fxr; Pin => Core%Pin; CellInfo => Core%CellInfo
AsyType => Core%AsyInfo; Asy => Core%Asy
CoreIdx => Core%CoreIdx; CoreMap => Core%CoreMap
hz => Core%hz
nxy = Core%nxy
nxa = Core%nxa; nya = Core%nya; nxya = Core%nxya
nxya0 = nxa * nya
nfsr = Core%nCoreFsr; nFxr = Core%nCoreFxr
myzb = PE%myzb; myze = PE%myze

ALLOCATE(AsyND(nxya0))
ALLOCATE(AsyND_Fuel(nxya0))
CALL CP_CA(AsyND, 0._8, nxya0)
CALL CP_CA(AsyND_Fuel, 0._8, nxya0)

i = 0
DO iya = 1, Core%nya
  DO ixa = 1, Core%nxa
    i = i + 1
    iasy = CoreIdx(ixa, iya) !Assembly Idx
    IF(iasy .EQ. 0) CYCLE
    iasytype = CoreMap(iasy)
    IF(.NOT. AsyType(iasytype)%lFuel) CYCLE
    !Assembly Sweeping
    nxy0 = AsyType(iasytype)%nxy
    volsum = 0; NDsum = 0;FuelVolSum = 0
    DO iz = myzb, myze
      IF(.NOT. Core%lFuelPlane(iz)) CYCLE
      DO ixy = 1, nxy0
        ipin = Asy(iasy)%GlobalPinIdx(ixy)
        icel = Pin(ipin)%Cell(iz)
        FxrIdxSt = Pin(ipin)%FxrIdxSt
        nLocalFxr = CellInfo(icel)%nFxr
        volsum = volsum + Core%PinVol(ixy, iz)
        DO j = 1, nLocalFxr
          ifxr = FxrIdxSt + j -1
          !IF(.NOT. Fxr(ifxr, iz)%lDepl) CYCLE
          lExist = .FALSE.
          DO iso = 1, Fxr(ifxr, iz)%niso_depl
            id = Fxr(ifxr, iz)%idiso(iso)
            IF(id .NE. idiso) CYCLE
            lExist = .TRUE.
            NDsum = Ndsum + Fxr(ifxr, iz)%pnum(iso) * Fxr(ifxr, iz)%area * hz(iz)
          ENDDO
          IF(lExist) FuelVolSum = FuelVolSum + Fxr(ifxr, iz)%area * hz(iz)
          !IF(.NOT. Fxr(ifxr, iz)%lDepl) CYCLE
          !volsum = volsum + Fxr(ifxr, iz)%area * hz(iz)
          !NDsum = NDsum + Fxr(ifxr, iz)%burnup * Fxr(ifxr, iz)%area * hz(iz)
        ENDDO
      ENDDO
    ENDDO
#ifdef MPI_ENV
    buf0 = (/volsum, FuelVolSum, NDsum/)
    CALL REDUCE(buf0, buf, 3, PE%MPI_CMFD_COMM, .TRUE.)
    volsum = buf(1); FuelVolSum = Buf(2); NDsum = buf(3)
#endif
    AsyND(i) = NDsum/volsum
    AsyND_Fuel(i) = NDsum/FuelVolSum
  ENDDO
ENDDO
IF(PE%MASTER) THEN
WRITE(io, '(A, I7, 2x, A)') '   Isotope :', IdIso, '(#/barn/cm)'
WRITE(io, '(A)') '     Assembly Averaged Number Density'
WRITE(formatdata,'(a)') '(8x, 1p200e12.3)'
CALL PrintReal1DarrayTo2Darray(io,  AsyND(1:nxya0), nxya0, nxa, nya, formatdata)
maxv = maxval(AsyND(1:nxya0))
WRITE(io, *)
WRITE(io, '(A)') '     Loaded Region Averaged Number Density'
CALL PrintReal1DarrayTo2Darray(io,  AsyND_Fuel(1:nxya0), nxya0, nxa, nya, formatdata)
maxv = maxval(AsyND_Fuel(1:nxya0))
WRITE(io, *)
ENDIF
DEALLOCATE(AsyND)
END SUBROUTINE

SUBROUTINE PrintAsyND2(Core, FmInfo, CmInfo, nTracerCntl, PE)  !Changhyun
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,        FmInfo_Type,       CmInfo_Type,     &
                             PE_Type,                                                  &
                             FxrInfo_Type,         Pin_Type,          Cell_Type,       &
                             AsyInfo_Type,         Asy_Type
USE CNTL,             ONLY : nTRACERCntl_Type
USE files,            ONLY : caseid
USE BasicOperation,   ONLY : CP_CA
USE IOUTIL,           ONLY : PrintReal1DarrayTo2Darray
USE MPIComm_Mod,      ONLY : REDUCE
USE Depl_Mod,         ONLY : DeplCntl
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(AsyInfo_Type), POINTER :: AsyType(:)
TYPE(Asy_Type), POINTER :: Asy(:)

INTEGER,  POINTER :: CoreIdx(:, :), CoreMap(:)
REAL, POINTER :: hz(:)
REAL, POINTER :: AsyND(:), AsyND_Fuel(:)


REAL :: volsum, NDsum, FuelVolSum, maxv
REAL :: buf0(3), buf(3)

INTEGER :: idiso, nisoout, ii
INTEGER :: myzb, myze, nfxr, nfsr, nxy, nxya, nxya0, nxa, nya
INTEGER :: nyc, nxc, nxy0, nLocalFxr, FxrIdxSt
INTEGER :: iz, ifxr, ifsr, icel, ipin, ixya, ixa, iya, iasy, iasytype
INTEGER :: iy, ix, ixy, i, j, k, iso, id
CHARACTER :: localfn*100, formatdata*120
LOGICAL :: lExist
LOGICAL :: lOpened
#if defined __GFORTRAN__ || __PGI
CHARACTER(5) :: str2num
INTEGER :: n
#endif

IF(.NOT. nTracerCntl%lXsLib) RETURN


Fxr => FmInfo%Fxr; Pin => Core%Pin; CellInfo => Core%CellInfo
AsyType => Core%AsyInfo; Asy => Core%Asy
CoreIdx => Core%CoreIdx; CoreMap => Core%CoreMap
hz => Core%hz
nxy = Core%nxy
nxa = Core%nxa; nya = Core%nya; nxya = Core%nxya
nxya0 = nxa * nya
nfsr = Core%nCoreFsr; nFxr = Core%nCoreFxr
myzb = PE%myzb; myze = PE%myze

nIsoOut = nTracerCntl%OutpCntl%IsoOutList(0)
INQUIRE(UNIT=300, OPENED=lOpened)
#if defined __GFORTRAN__ || __PGI
IF(.NOT.lOpened) THEN
  n = nIsoout * nxya0
  WRITE(str2num, '(I5)') n
  localfn = trim(caseid)//'_AsyND.out'
  OPEN(UNIT=300,FILE=localfn,STATUS='REPLACE',FORM='FORMATTED')
  WRITE(300, '(a10,' // TRIM(str2num) // 'I16)') 'Isotope',((nTracerCntl%OutpCntl%IsoOutList(iso), ixy=1,nxya0), iso=1,nIsoOut)
  WRITE(300, '(a10,' // TRIM(str2num) // '(6X,"(",I3,", ",I3,")"))') 'Asy Coord.',(((iy, ix, ix=1,nxa), iy=1,nya),iso=1,nIsoOut)
  localfn = trim(caseid)//'_LRND.out'
  OPEN(UNIT=301,FILE=localfn,STATUS='REPLACE',FORM='FORMATTED')
  WRITE(301, '(a10,' // TRIM(str2num) // 'I16)') 'Isotope',((nTracerCntl%OutpCntl%IsoOutList(iso), ixy=1,nxya0), iso=1,nIsoOut)
  WRITE(301, '(a10,' // TRIM(str2num) // '(6X,"(",I3,", ",I3,")"))') 'Asy Coord.',(((iy, ix, ix=1,nxa), iy=1,nya),iso=1,nIsoOut)
ENDIF
#else
IF(.NOT.lOpened) THEN
  localfn = trim(caseid)//'_AsyND.out'
  OPEN(UNIT=300,FILE=localfn,STATUS='REPLACE',FORM='FORMATTED')
  WRITE(300, '(a10,<nIsoOut*nxya0>I16)') 'Isotope',((nTracerCntl%OutpCntl%IsoOutList(iso), ixy=1,nxya0), iso=1,nIsoOut)
  WRITE(300, '(a10,<nIsoOut*nxya0>(6X,"(",I3,", ",I3,")"))') 'Asy Coord.',(((iy, ix, ix=1,nxa), iy=1,nya),iso=1,nIsoOut)
  localfn = trim(caseid)//'_LRND.out'
  OPEN(UNIT=301,FILE=localfn,STATUS='REPLACE',FORM='FORMATTED')
  WRITE(301, '(a10,<nIsoOut*nxya0>I16)') 'Isotope',((nTracerCntl%OutpCntl%IsoOutList(iso), ixy=1,nxya0), iso=1,nIsoOut)
  WRITE(301, '(a10,<nIsoOut*nxya0>(6X,"(",I3,", ",I3,")"))') 'Asy Coord.',(((iy, ix, ix=1,nxa), iy=1,nya),iso=1,nIsoOut)
ENDIF
#endif
ALLOCATE(AsyND(nxya0))
ALLOCATE(AsyND_Fuel(nxya0))
CALL CP_CA(AsyND, 0._8, nxya0)
CALL CP_CA(AsyND_Fuel, 0._8, nxya0)

WRITE(300,'(F10.3)',ADVANCE='NO') DeplCntl%NowBurnup(2)
WRITE(301,'(F10.3)',ADVANCE='NO') DeplCntl%NowBurnup(2)
DO ii = 1, nIsoOut
  idiso = nTracerCntl%outpCntl%IsoOutList(ii)
  i = 0
  DO iya = 1, Core%nya
    DO ixa = 1, Core%nxa
      i = i + 1
      iasy = CoreIdx(ixa, iya) !Assembly Idx
      IF(iasy .EQ. 0) CYCLE
      iasytype = CoreMap(iasy)
      IF(.NOT. AsyType(iasytype)%lFuel) CYCLE
      !Assembly Sweeping
      nxy0 = AsyType(iasytype)%nxy
      volsum = 0; NDsum = 0;FuelVolSum = 0
      DO iz = myzb, myze
        IF(.NOT. Core%lFuelPlane(iz)) CYCLE
        DO ixy = 1, nxy0
          ipin = Asy(iasy)%GlobalPinIdx(ixy)
          icel = Pin(ipin)%Cell(iz)
          FxrIdxSt = Pin(ipin)%FxrIdxSt
          nLocalFxr = CellInfo(icel)%nFxr
          volsum = volsum + Core%PinVol(ixy, iz)
          DO j = 1, nLocalFxr
            ifxr = FxrIdxSt + j -1
            !IF(.NOT. Fxr(ifxr, iz)%lDepl) CYCLE
            lExist = .FALSE.
            DO iso = 1, Fxr(ifxr, iz)%niso_depl
              id = Fxr(ifxr, iz)%idiso(iso)
              IF(id .NE. idiso) CYCLE
              lExist = .TRUE.
              NDsum = Ndsum + Fxr(ifxr, iz)%pnum(iso) * Fxr(ifxr, iz)%area * hz(iz)
            ENDDO
            IF(lExist) FuelVolSum = FuelVolSum + Fxr(ifxr, iz)%area * hz(iz)
            !IF(.NOT. Fxr(ifxr, iz)%lDepl) CYCLE
            !volsum = volsum + Fxr(ifxr, iz)%area * hz(iz)
            !NDsum = NDsum + Fxr(ifxr, iz)%burnup * Fxr(ifxr, iz)%area * hz(iz)
          ENDDO
        ENDDO
      ENDDO
#ifdef MPI_ENV
      buf0 = (/volsum, FuelVolSum, NDsum/)
      CALL REDUCE(buf0, buf, 3, PE%MPI_CMFD_COMM, .TRUE.)
      volsum = buf(1); FuelVolSum = Buf(2); NDsum = buf(3)
#endif
      AsyND(i) = NDsum/volsum
      AsyND_Fuel(i) = NDsum/FuelVolSum
    ENDDO
  ENDDO
#if defined __GFORTRAN__ || __PGI
  IF(PE%MASTER) THEN
    WRITE(str2num, '(I5)') nxya0
    WRITE(300,'(' // TRIM(str2num) // 'ES16.6)',ADVANCE='NO') AsyND(1:nxya0)
    WRITE(301,'(' // TRIM(str2num) // 'ES16.6)',ADVANCE='NO') AsyND_Fuel(1:nxya0)
  ENDIF
#else
  IF(PE%MASTER) THEN
    WRITE(300,'(<nxya0>ES16.6)',ADVANCE='NO') AsyND(1:nxya0)
    WRITE(301,'(<nxya0>ES16.6)',ADVANCE='NO') AsyND_Fuel(1:nxya0)
  ENDIF
#endif
ENDDO
WRITE(300,*)
WRITE(301,*)

DEALLOCATE(AsyND,AsyND_Fuel)
END SUBROUTINE
