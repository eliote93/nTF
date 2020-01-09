#include <defines.h>
MODULE SubChannel_mod
USE PARAM
#ifdef HAVE_MATRA
USE Matra_Mod,     ONLY : MatNtData_TYPE,                                     &
                          Matra_INIT,         DriveMatra
#endif                          
USE typedef, ONLY : CoreInfo_Type,       FmInfo_Type,        Pin_Type,        &
                    CmInfo_Type,         ThInfo_Type,        GroupInfo_Type,  &
                    Asy_Type,            AsyInfo_Type,       CoolantTH_Type,  &
                    PE_TYPE
USE CNTL,    ONLY : nTRACERCntl_Type
use SubChCoupling_mod,      only: last_TH
IMPLICIT NONE
#ifdef HAVE_MATRA
TYPE(MatNtData_TYPE), PRIVATE :: MatNtData
#endif
TYPE(Pin_TYPE), PRIVATE, POINTER :: Pin(:)
TYPE(Asy_TYPE), PRIVATE, POINTER :: Asy(:)
TYPE(ASYINFO_TYPE), PRIVATE, POINTER :: AsyInfo(:)
TYPE(CoolantTh_TYPE), PRIVATE, POINTER :: CoolantTH(:)

INTEGER, PRIVATE, POINTER :: CoreMap(:)

LOGICAL, PRIVATE, POINTER :: lFuelPlane(:)
LOGICAL, PRIVATE :: lGap = .FALSE.
REAL, PRIVATE, POINTER :: RelPow(:, :)
REAL, PRIVATE, POINTER :: DenCool(:, :), Tcool(:, :)

INTEGER, PRIVATE :: myzb, myze, nz, nz_act
INTEGER, PRIVATE :: nxya, nxy, nasy_act
LOGICAL, SAVE :: lFirst = .TRUE.


CONTAINS

SUBROUTINE SubChannelTH(Core, FmInfo, CmInfo, ThInfo, GroupInfo, nTracerCntl, PE)
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

INTEGER :: iz, iasy, ixa, iya, iasytype

CHARACTER(50) :: FN 
#ifdef HAVE_MATRA
myzb = PE%myzb; myze = PE%myze
nxy = Core%nxy; nz = COre%nz;
nxya = Core%nxya
Pin => Core%Pin; Asy => Core%asy; AsyInfo => Core%AsyInfo
CoreMap => Core%Coremap; lFuelPlane => Core%lFuelPlane
CoolantTH => ThInfo%CoolantTH
TCool => ThInfo%TCool; DenCool => ThInfo%DenCool
lGap = Core%lgAp

RelPow => ThInfo%RelPower

nz_act =0; nasy_act = 0
DO iz = 1, Core%nz
  IF(.NOT. Core%lFuelPlane(iz)) CYCLE
  nz_act = nz_act + 1
ENDDO

DO iasy = 1, Core%nxya
    ixa = Asy(iasy)%ixa; iya = Asy(iasy)%iya
    iasytype = Core%CoreMap(iasy)
    IF(.NOT. AsyInfo(iasyType)%lFuel) CYCLE
    nasy_act = nasy_act + 1
ENDDO

IF(.NOT. MatNtData%lAlloc) CALL AllocMatNtData()
IF(PE%CmfdMaster) THEN
  IF(lFirst) THEN
    FN = 'ntracer.matra'
    IF(PE%CmfdMASTER) CALL Matra_INIT(FN, MatNtData)  
    lFirst = .FALSE.
  ENDIF
  CALL SetMatraData()
  CALL DriveMatra(MatNtData,last_TH)
ENDIF
CALL BcastMatraOut(PE)
CALL GetMatraData()
#endif
CONTINUE
END SUBROUTINE

SUBROUTINE BcastMatraOut(PE)
#ifdef MPI_ENV
USE MPIComm_Mod, ONLY : BCAST, MPI_SYNC
#endif
IMPLICIT NONE
TYPE(PE_TYPE) :: PE

INTEGER :: comm
INTEGER :: iasy, n1, n2, n3
REAL, POINTER :: Dat(:, :, :)
#ifdef HAVE_MATRA
#ifdef MPI_ENV
COMM = PE%MPI_CMFD_COMM
CALL MPI_SYNC(COMM)
DO iasy = 1, nasy_act
  Dat => MatntData%AsyDat(iasy)%Dat
  n1 = MatNtData%AsyDat(iasy)%nx;
  n2 = MatNtData%AsyDat(iasy)%ny;
  n3 = nz_act+2
  CALL BCAST(MatntData%AsyDat(iasy)%Dat, n1, n2, n3, COMM)
  Dat => MatntData%AsyDat(iasy)%Dat2
  CALL BCAST(MatntData%AsyDat(iasy)%Dat2, n1, n2, n3, COMM)
ENDDO
#endif
#endif
END SUBROUTINE

SUBROUTINE GetMatraData()
IMPLICIT NONE

INTEGER :: iz, iz_act, iasy_act, iasytype
INTEGER :: ix, iy, ixy0, ixy, ixa, iya, iasy
INTEGER :: nxc, nyc
#ifdef HAVE_MATRA
iz_act = 0
DO iz = 1, nz
  IF(.NOT. lFuelPlane(iz)) CYCLE
  iz_act = iz_act + 1
  iasy_act = 0
  DO iasy = 1, nxya
    ixa = Asy(iasy)%ixa; iya = Asy(iasy)%iya
    iasytype = CoreMap(iasy)
    IF(.NOT. AsyInfo(iasyType)%lFuel) CYCLE
    iasy_act = iasy_act + 1
    nxc = AsyInfo(iasyType)%nx; nyc = AsyInfo(iasyType)%ny    
    DO iy = 1, nyc
      DO ix = 1, nxc
        ixy0 = AsyInfo(iasyType)%Pin2dIdx(ix, iy)
        ixy = Asy(iasy)%GlobalPinIdx(ixy0)
        !MatNtData%AsyDat(iasy_act)%Dat(ix, iy, iz_act) = RelPow(iz, ixy)
        Tcool(iz, ixy) = MatNtData%AsyDat(iasy_act)%Dat(ix, iy, iz_act)
        DenCool(iz, ixy) = MatNtData%AsyDat(iasy_act)%Dat2(ix, iy, iz_act)
      ENDDO
    ENDDO    
  ENDDO
ENDDO

iz_act = 0
DO iz = 1, nz
  IF(lFuelPlane(iz)) EXIT
  iasy_act = 0
  DO iasy = 1, nxya
    ixa = Asy(iasy)%ixa; iya = Asy(iasy)%iya
    iasytype = CoreMap(iasy)
    IF(.NOT. AsyInfo(iasyType)%lFuel) CYCLE
    iasy_act = iasy_act + 1
    nxc = AsyInfo(iasyType)%nx; nyc = AsyInfo(iasyType)%ny    
    DO iy = 1, nyc
      DO ix = 1, nxc
        ixy0 = AsyInfo(iasyType)%Pin2dIdx(ix, iy)
        ixy = Asy(iasy)%GlobalPinIdx(ixy0)
        !MatNtData%AsyDat(iasy_act)%Dat(ix, iy, iz_act) = RelPow(iz, ixy)
        Tcool(iz, ixy) = MatNtData%AsyDat(iasy_act)%Dat(ix, iy, iz_act)
        DenCool(iz, ixy) = MatNtData%AsyDat(iasy_act)%Dat2(ix, iy, iz_act)
      ENDDO
    ENDDO    
  ENDDO  
ENDDO

iz_act = nz_act + 1
DO iz = nz, 1, -1
  IF(lFuelPlane(iz)) EXIT
  iasy_act = 0
  DO iasy = 1, nxya
    ixa = Asy(iasy)%ixa; iya = Asy(iasy)%iya
    iasytype = CoreMap(iasy)
    IF(.NOT. AsyInfo(iasyType)%lFuel) CYCLE
    iasy_act = iasy_act + 1
    nxc = AsyInfo(iasyType)%nx; nyc = AsyInfo(iasyType)%ny    
    DO iy = 1, nyc
      DO ix = 1, nxc
        ixy0 = AsyInfo(iasyType)%Pin2dIdx(ix, iy)
        ixy = Asy(iasy)%GlobalPinIdx(ixy0)
        !MatNtData%AsyDat(iasy_act)%Dat(ix, iy, iz_act) = RelPow(iz, ixy)
        Tcool(iz, ixy) = MatNtData%AsyDat(iasy_act)%Dat(ix, iy, iz_act)
        DenCool(iz, ixy) = MatNtData%AsyDat(iasy_act)%Dat2(ix, iy, iz_act)
      ENDDO
    ENDDO    
  ENDDO  
ENDDO
#endif
END SUBROUTINE


SUBROUTINE SetMatraData()
IMPLICIT NONE

INTEGER :: iz, iz_act, iasy_act, iasytype
INTEGER :: ix, iy, ixy0, ixy, ixa, iya, iasy
INTEGER :: nxc, nyc
#ifdef HAVE_MATRA
iz_act = 0
DO iz = 1, nz
  IF(.NOT. lFuelPlane(iz)) CYCLE
  iz_act = iz_act + 1
  iasy_act = 0
  DO iasy = 1, nxya
    ixa = Asy(iasy)%ixa; iya = Asy(iasy)%iya
    iasytype = CoreMap(iasy)
    IF(.NOT. AsyInfo(iasyType)%lFuel) CYCLE
    iasy_act = iasy_act + 1
    nxc = AsyInfo(iasyType)%nx; nyc = AsyInfo(iasyType)%ny    
    DO iy = 1, nyc
      DO ix = 1, nxc
        ixy0 = AsyInfo(iasyType)%Pin2dIdx(ix, iy)
        ixy = Asy(iasy)%GlobalPinIdx(ixy0)
        MatNtData%AsyDat(iasy_act)%Dat(ix, iy, iz_act) = RelPow(iz, ixy)
      ENDDO
    ENDDO    
  ENDDO
ENDDO
#endif
END SUBROUTINE

SUBROUTINE AllocMatNtData()
USE ALLOCS
IMPLICIT NONE


INTEGER :: nxc, nyc
INTEGER :: ixy, ixy0, ix, iy, iz, iasy, ixa, iya, iasytype
INTEGER :: iz_act, iasy_act

#ifdef HAVE_MATRA
MatNtData%nz =nz_act; MatNtData%nasy = nasy_act
MatNtData%lAlloc = .TRUE.
MatNtData%lGap =lGap

ALLOCATE(MatNtData%AsyDat(nasy_act))
iz = 1
iasy_act = 0
DO iasy = 1, nxya
  ixa = Asy(iasy)%ixa; iya = Asy(iasy)%iya
  iasytype = CoreMap(iasy)
  IF(.NOT. AsyInfo(iasyType)%lFuel) CYCLE
  nxc = AsyInfo(iasyType)%nx; nyc = AsyInfo(iasyType)%ny
  iasy_act = iasy_act + 1
  CALL Dmalloc0(MatNtData%AsyDat(iasy_act)%Dat, 1, nxc, 1, nyc, 0, nz_act+1)
  CALL Dmalloc0(MatNtData%AsyDat(iasy_act)%Dat2, 1, nxc, 1, nyc, 0, nz_act+1)
  MatNtData%AsyDat(iasy_act)%nx = nxc
  MatNtData%AsyDat(iasy_act)%ny = nyc
  MatNtData%AsyDat(iasy_act)%iasyx = ixa
  MatNtData%AsyDat(iasy_act)%iasyy = iya
ENDDO
#endif
END SUBROUTINE
!
END MODULE
