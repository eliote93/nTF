#include <defines.h>
SUBROUTINE SetUserDefTH(Core, CmInfo, FmInfo, ThInfo, nTRACERCntl, PE)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,        CMInfo_Type,      FmInfo_Type,     &
                           ThInfo_Type,          GroupInfo_Type,   PE_Type,         &
                           FxrInfo_Type,         PinXS_Type,       Pin_Type,        &
                           ThVar_Type,           THOpt_Type,       FuelTh_Type,     &
                           CoolantTH_Type,       AsyInfo_Type
USE TH_Mod,         ONLY : ThOpt,                ThVar,                             &
                           ReadUserDefTemp
USE BasicOperation, ONLY : CP_CA,                CP_VA
USE CNTL,           ONLY : nTracerCntl_Type
USE itrcntl_mod,    ONLY : ItrCntl
USE FILES,          ONLY : IO8
USE IOUTIL,         ONLY : message
USE TIMER,          ONLY : nTracer_dclock,       TimeChk
USE SteamTBL_mod,  ONLY : steamtbl
#ifdef MPI_ENV
USE MPIComm_Mod,    ONLY : MPI_SYNC
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(FuelTH_Type), POINTER :: FuelTH(:)
TYPE(CoolantTH_Type), POINTER :: CoolantTH(:)

INTEGER :: nxy, nz
INTEGER :: iz, ixy, iasy, iasytype


REAL :: PExit, wt, wh, hin, wrho, wvin, wxin, wbetain, wkapain, wcpin

LOGICAL, SAVE :: lFirst
DATA lFirst /.TRUE./

IF(lFirst) THEN
  CALL ReadUserDefTemp(Core, ThInfo, PE)
  lFirst = .FALSE.
ENDIF

WRITE(mesg, '(A)') 'Overwriting Temperature Distribution to User Defined Values...'
IF(PE%MASTER) CALL MESSAGE(io8, TRUE, TRUE, mesg)

nxy = Core%nxy; nz = Core%nz
Pin => Core%Pin; AsyInfo => Core%AsyInfo

CoolantTH => THInfo%CoolantTH
FuelTH => THInfo%FuelTH
PEXIT = THVar%PEXIT
DO ixy = 1, nxy
  iasy = Pin(ixy)%iasy
  iasytype = Pin(ixy)%AsyType
  IF(.NOT. AsyInfo(iasytype)%lfuel) CYCLE
  DO iz = 1, nz
    wt = ThInfo%UserDefModT(iasy, iz) + CKELVIN
    CALL SteamTBL(TRUE, pexit, wt, hin, wrho, wvin, wxin, wbetain, wkapain, wcpin)
    CoolantTH(ixy)%DenCool(iz) = wrho
    CoolantTH(ixy)%Tcool(iz) = wt - CKELVIN
    FuelTH(ixy)%tfvol(:, iz) = ThInfo%UserDefFuelT(iasy, iz)
    FuelTH(ixy)%tfuel(:, iz) = ThInfo%UserDefFuelT(iasy, iz)
  ENDDO
ENDDO


END SUBROUTINE

SUBROUTINE ReadUserDefTemp(Core, ThInfo, PE)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,       ThInfo_Type,          PE_Type,          &
                           Asy_Type,            AsyInfo_Type
USE FILES,          ONLY : Filename,            ModTFileIdx,          FuelTFileIdx,    &
                           io4,                 io8
USE IOUTIL,         ONLY : OPENFILE
USE ALLOCS
#ifdef MPI_ENV
USE MPIComm_Mod,    ONLY : MPIWaitTurn
#endif
TYPE(CoreInfo_Type) :: Core
TYPE(ThInfo_Type) :: ThInfo
TYPE(PE_Type) :: PE


TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)

REAL, POINTER :: dat(:)
INTEGER :: nxya, nz
INTEGER :: nz0, nxyaf

INTEGER :: ixyaf, iz

Asy => Core%Asy; AsyInfo => Core%AsyInfo
nxya = Core%nxya; nz = Core%nz
!
CALL DMALLOC(ThInfo%UserDefFuelT, nxya, nz)
CALL DMALLOC(ThInfo%UserDefModT, nxya, nz)
#ifdef MPI_ENV
CALL MPIWaitTurn(PE%MPI_Cmfd_Comm, PE%myCmfdRank, PE%nCmfdproc, .FALSE.)
#endif
!Read Fuel Temperature
CALL OPENFILE(io4, TRUE, FALSE, FALSE, filename(FuelTFileIdx))
READ(io4, *) nz0, nxyaf
ALLOCATE(dat(nxyaf))
DO iz = 1, nz
  READ(io4, *) (dat(ixyaf), ixyaf = 1, nxyaf)  
  IF(iz .NE. nz) READ(io4, *)
  ixyaf=0
  DO ixya = 1, nxya
    iasytype = Asy(ixya)%AsyType
    ThInfo%UserDefFuelT(ixya, iz) = 0
    IF(.NOT. AsyInfo(iasytype)%lFuel) CYCLE
    ixyaf = ixyaf + 1
    ThInfo%UserDefFuelT(ixya, iz) = Dat(ixyaf)
  ENDDO
ENDDO
DEALLOCATE(dat)
CLOSE(io4)
!Read Mod Temperature
CALL OPENFILE(io4, TRUE, TRUE, FALSE, filename(ModTFileIdx))
READ(io4, *) nz0, nxyaf
ALLOCATE(dat(nxyaf))
DO iz = 1, nz
  READ(io4, *) (dat(ixyaf), ixyaf = 1, nxyaf)  
  IF(iz .NE. nz) READ(io4, *)
  DO ixya = 1, nxya
    iasytype = Asy(ixya)%AsyType
    IF(.NOT. AsyInfo(iasytype)%lFuel) CYCLE
    ixyaf = ixyaf + 1
    ThInfo%UserDefModT(ixya, iz) = Dat(ixyaf)
  ENDDO
ENDDO
DEALLOCATE(dat)
CLOSE(io4)
#ifdef MPI_ENV
CALL MPIWaitTurn(PE%MPI_Cmfd_Comm, PE%myCmfdRank, PE%nCmfdproc, .TRUE.)

#endif
NULLIFY(ASY, AsyInfo)

END SUBROUTINE