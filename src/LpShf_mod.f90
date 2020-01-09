#include <defines.h>
MODULE LpShf_mod
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,         FmInfo_Type,           PE_TYPE,             &
                           FxrInfo_Type,          AsyInfo_Type,          Asy_Type,            &
                           Pin_Type,              PinInfo_Type,          Cell_Type,           &
                           GroupInfo_Type,        LpShf_Type,            Shf_Type
USE files,          ONLY : io8, io13, io14, caseid, localfn
USE ioutil,         ONLY : openfile,              message,               nFields
USE CNTL,           ONLY : nTracerCntl_Type
USE BasicOperation, ONLY : CP_CA,                 CP_VA
USE Restart_mod,    ONLY : ReadShfRstFile
#ifdef MPI_ENV
USE MPIComm_Mod,    ONLY : MPIWaitTurn
#endif
IMPLICIT NONE
SAVE

TYPE(LPShf_Type) :: LpSHF
INTEGER :: RmvBP(2, 200)
INTEGER :: nRmvBP = 0
TYPE(FxrInfo_Type), PRIVATE, POINTER :: Fxr(:, :)
TYPE(AsyInfo_Type), PRIVATE, POINTER :: AsyType(:)
TYPE(Asy_Type), PRIVATE, POINTER :: Asy(:)
TYPE(Pin_Type), PRIVATE, POINTER :: Pin(:)
TYPE(Cell_Type), PRIVATE, POINTER :: CellInfo(:)
TYPE(PinInfo_Type), PRIVATE, POINTER :: PinInfo(:)

INTEGER, PRIVATE, POINTER :: CoreIdx(:, :), CoreMap(:)
INTEGER, PRIVATE :: nxc, nyc
INTEGER, PRIVATE :: nxa, nya, nxya

CHARACTER(1), PRIVATE :: IdxX(23)
DATA IdxX /'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'J', 'K', 'L', 'M', 'N', &
           'P', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'/

CONTAINS

SUBROUTINE LP_Shuffling(Core, FmInfo, GroupInfo, nTRACERCntl, PE)
USE Cooling_mod, ONLY : Init_Cooling, Fin_Cooling
USE MPIComm_mod, ONLY : MPI_SYNC
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTRACERCntl
TYPE(PE_TYPE) :: PE

INTEGER :: ixya, ixa, iya, iasy, iz
INTEGER :: iasytype
LOGICAL :: LPartialAsy
!File IO
INTEGER :: nchar, nproc
CHARACTER(256) :: Fn, AsyFn
!Assembly ID
!!MPI Variables
INTEGER :: irank
LOGICAL :: Master, Slave
!
INTEGER :: icyc

Master =  PE%CmfdMaster; Slave = PE%CmfdSlave
!
Fxr => FmInfo%Fxr; Pin => Core%Pin
AsyType => Core%AsyInfo; Asy => Core%Asy
CoreIdx => Core%CoreIdx; CoreMap => Core%CoreMap
CellInfo => Core%CellInfo; PinInfo => Core%PinInfo
nxc = Core%nxc; nyc = Core%nyc
nxya = Core%nxya; nxa = Core%nxa; nya =Core%nya

IF(nTracerCntl%lCooling) THEN
  CALL Init_Cooling(Core, FmInfo, GroupInfo, PE)
  CALL ReloadFuelCooling(nTRACERCntl, PE)
ENDIF
CALL MPI_SYNC(PE%MPI_CMFD_COMM)
NULLIFY(Fxr, Pin, AsyType, Asy)
NULLIFY(CoreIdx, CoreMap, CellInfo)

IF (nTracerCntl%LPShf) THEN
  mesg = 'Loading Pattern Shuffling ...'
  IF(nTRACERCntl%lRstCal) CALL SetRstMod()
  IF(master) CALL message(io8, TRUE, TRUE, mesg)
  DO icyc = 1, LpShf%nrstfile
    CALL ReadShfRstFile(Core, FmInfo, icyc, LpShf, nTRACERCntl, PE)
  ENDDO
  IF(nTRACERCntl%lRstCal) THEN
    icyc = LpShf%RstCycleId
    Core%Hm_Mass0 = LpShf%Hm_mass0(icyc)
  ENDIF
  !CALL FreshFuelLoading(nTRACERCntl, PE)
  IF(nRmvBP .GT. 0) THEN
    mesg = 'Remove BP Cell from Previous Cycles ...'
    IF(master) CALL message(io8, TRUE, TRUE, mesg)
    CALL CoreRmvBP(nTracerCntl, PE)
  ENDIF
  CALL MPI_SYNC(PE%MPI_CMFD_COMM)
  IF(LpShf%lPUL) THEN
    CALL Init_Cooling(Core, FmInfo, GroupInfo, PE)
    CALL ReloadFuelCooling(nTRACERCntl, PE)
    CALL Fin_Cooling
  ENDIF
  CALL MPI_SYNC(PE%MPI_CMFD_COMM)
ENDIF

NULLIFY(Fxr, Pin, AsyType, Asy)
NULLIFY(CoreIdx, CoreMap, CellInfo)

END SUBROUTINE


SUBROUTINE SetRstMod()
INTEGER :: ix, iy, iasy, iasytype

DO iy = 1, nya
  DO ix = 1, nxa
    iasy = CoreIdx(ix, iy)
    IF(iasy .EQ. 0) CYCLE
    iasytype = CoreMap(iasy)
    IF(.NOT. AsyType(iasytype)%lFuel) THEN
      LpShf%Shf(iasy)%Lfreshfuel = .FALSE.
      LpShf%Shf(iasy)%LNoneFuel = .TRUE.
      CYCLE
    ENDIF
    LpShf%Shf(iasy)%lFreshFuel = .FALSE.
    LpShf%Shf(iasy)%lNoneFuel = .FALSE.
    LpShf%Shf(iasy)%ix = ix;LpShf%Shf(iasy)%iy = iy
    LpShf%Shf(iasy)%irot = 0
    LpShf%Shf(iasy)%cycleid = LpShf%RstCycleId 
  ENDDO
ENDDO


END SUBROUTINE

SUBROUTINE FreshFuelLoading(nTRACERCntl, PE)
USE Material_mod,    ONLY : Mixture
IMPLICIT NONE
TYPE(nTracerCntl_Type) :: nTRACERCntl
TYPE(PE_TYPE) :: PE
INTEGER :: iasy, iasytype, ixy, ixy0, ipin, iz
INTEGER :: ipintype, iceltype, ifxr, imix
INTEGER :: FxrIdxSt, ifsrbeg, nLocalFxr, niso
INTEGER :: j

DO iasy = 1, nxya
  IF(.NOT. LpShf%shf(iasy)%lFreshFuel) CYCLE
  iasytype = CoreMap(iasy)
  IF(.NOT. AsyType(iasytype)%lfuel) CYCLE
  DO iz = PE%myzb, PE%myze
    DO ixy = 1, AsyType(iasytype)%nxy
      ipintype = AsyType(iasytype)%Pin(ixy)
      iceltype = PinInfo(ipintype)%cell(iz)
      ipin = Asy(iasy)%GlobalPinIdx(ixy)
      FxrIdxSt = Pin(ipin)%FxrIdxSt
      nLocalFxr = CellInfo(iceltype)%nFxr
      DO j = 1, nLocalFxr
        ifxr = FxrIdxSt + j - 1
        IF(.NOT. Fxr(ifxr, iz)%ldepl) CYCLE
        ifsrbeg = CellInfo(iceltype)%MapFxr2FsrIdx(1, j)
        imix = CellInfo(iceltype)%ireg(ifsrbeg)
        niso = Mixture(imix)%niso
        Fxr(ifxr, iz)%niso = niso; Fxr(ifxr, iz)%niso_depl = niso
        CALL CP_VA(Fxr(ifxr,iz)%idiso(1:niso),Mixture(imix)%idiso(1:niso),niso)
        CALL CP_VA(Fxr(ifxr,iz)%pnum(1:niso),Mixture(imix)%pnum(1:niso),niso)        
        IF(nTRACERCntl%lXeDyn .AND. Fxr(ifxr, iz)%lfuel) THEN
          Fxr(ifxr, iz)%niso = Fxr(ifxr, iz)%niso + 1
          Fxr(ifxr, iz)%niso_depl = Fxr(ifxr, iz)%niso_depl + 1
          Fxr(ifxr, iz)%idiso(Fxr(ifxr, iz)%niso) = 54635
          Fxr(ifxr, iz)%pnum(Fxr(ifxr, iz)%niso) = epsm10
         ENDIF         
      ENDDO
    ENDDO
  ENDDO
  !
ENDDO

END SUBROUTINE
!
SUBROUTINE ReloadFuelCooling(nTRACERCntl, PE)
!Process the depletion of the Reloaded fuel in the Pool
USE GEOM,    ONLY : Core
TYPE(nTracerCntl_TYPE) :: nTRACERCNTL
TYPE(PE_TYPE) :: PE
INTEGER :: iasy, iasytype, icyc
REAL :: CoolingTime
mesg = 'Reloaded Assembly Cooling...'
IF(PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)  
DO iasy = 1, Core%nxya
  IF(LpShf%shf(iasy)%lFreshFuel) CYCLE
  iasytype = CoreMap(iasy); 
  IF(.NOT. AsyType(iasytype)%lfuel) CYCLE
  IF (nTracerCntl%lnTIGRst) THEN
    icyc = 1
  ELSE
    icyc = LpShf%Shf(iasy)%cycleid
  ENDIF
  IF(LpShf%PUL(icyc) .EQ. 0._8) CYCLE
  CoolingTime = LpShf%PUL(icyc)
  CoolingTime = CoolingTime * 24._8* 3600._8;
  CALL AssemblyCooling(iasy, iasytype, CoolingTime, nTracerCntl%lXeDyn, PE)
ENDDO

END SUBROUTINE

SUBROUTINE AssemblyCooling(iasy, iasytype, CoolingTime, lXeDyn, PE)
USE Cooling_mod,      ONLY : FxrCooling
USE OMP_LIB
IMPLICIT NONE
INTEGER :: iasy, iasytype
REAL :: CoolingTIme
LOGICAL :: lXeDyn
TYPE(PE_TYPE) :: PE

INTEGER :: i, j, k
INTEGER :: iz, ixy0, ixy, icel, ifxr
INTEGER :: FxrIdxst, nLocalFxr
INTEGER :: tid


!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(PE%nDeplThread) 

DO iz = PE%myzb, PE%myze
!$OMP  PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(j, tid, ixy0, ixy, icel, ifxr, FxrIdxst, nLocalFxr)
!$  tid = omp_get_thread_num()+1
!$OMP DO ORDERED SCHEDULE(DYNAMIC)  
  DO ixy0 = 1, AsyType(iasytype)%nxy
    ixy = Asy(iasy)%GlobalPinIdx(ixy0);    icel = Pin(ixy)%Cell(iz)
    FxrIdxst = Pin(ixy)%FxrIdxst
    nLocalFxr = CellInfo(icel)%nFxr
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j - 1
      IF(.NOT. Fxr(ifxr, iz)%ldepl) CYCLE
      CALL FxrCooling(Fxr(ifxr, iz), CoolingTime, tid, lXeDyn)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL
ENDDO


END SUBROUTINE

SUBROUTINE CoreRmvBP(nTracerCntl, PE)
IMPLICIT NONE
TYPE(nTracerCntl_TYPE) :: nTRACERCNTL
TYPE(PE_TYPE) :: PE

INTEGER :: iasy, IASYTYPE

CALL SetBPCellInfo()

DO iasy = 1, nxya
  IF(LpShf%shf(iasy)%lFreshFuel) CYCLE
  iasytype = CoreMap(iasy); 
  IF(.NOT. AsyType(iasytype)%lfuel) CYCLE
  CALL AsyRmvBP(iasy, iasytype, nTracerCntl, PE)
ENDDO

END SUBROUTINE

SUBROUTINE AsyRmvBP(iasy, iasytype, nTracerCntl, PE)

USE Material_Mod
!USE Boron_mod, ONLY : FxrBoronUpdate
USE ioutil,    ONLY : Terminate
USE allocs
USE BasicOperation, ONLY : CP_VA, CP_CA
IMPLICIT NONE
INTEGER :: iasy, iasytype
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE

INTEGER :: i, j, k
INTEGER :: iz, ixy0, ixy, icel, ifxr, ireg 
INTEGER :: FxrIdxst, nLocalFxr, tid
INTEGER :: icel1
INTEGER :: imix, imix1, niso, niso1
LOGICAL :: lRmvBP
DO iz = PE%myzb, PE%myze
  DO ixy0 = 1, AsyType(iasytype)%nxy
    ixy = Asy(iasy)%GlobalPinIdx(ixy0);    icel = Pin(ixy)%Cell(iz)
    lRmvBP = .FALSE.
    DO k = 1, nRmvBP
      IF(icel .EQ. RmvBP(1, k)) THEN
        lRmvBP = .TRUE.
        icel1 = RmvBP(2, k)
        EXIT
      ENDIF
    ENDDO
    IF(.NOT. lRmvBP) CYCLE
    FxrIdxst = Pin(ixy)%FxrIdxst
    nLocalFxr = CellInfo(icel)%nFxr
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j - 1
      ireg = CellInfo(icel)%MapFxr2FsrIdx(1, j)
      IF(CellInfo(icel)%iReg(ireg) .eq. CellInfo(icel1)%iReg(ireg)) CYCLE
      imix = CellInfo(icel)%iReg(ireg)
      imix1 = CellInfo(icel1)%iReg(ireg)
      
      IF(mixture(imix)%lfuel .or. mixture(imix1)%lfuel) THEN
         CALL TERMINATE('RMV_BP : Fuel Mixture Replace is not allowed.')  
      ENDIF
      IF(.not. mixture(imix)%lres .AND. mixture(imix1)%lres) THEN
         CALL TERMINATE('RMV_BP : Fuel Mixture Replace from non-Res. Mix. to Res. Mix. is not allowed')  
      ENDIF
      IF(.not. mixture(imix)%ldepl .AND. mixture(imix1)%ldepl) THEN
         CALL TERMINATE('RMV_BP : Fuel Mixture Replace from non-Depl. Mix. to Depl. Mix. is not allowed')  
      ENDIF
      Fxr(ifxr, iz)%lRes = mixture(imix1)%lres;       Fxr(ifxr, iz)%ldepl = mixture(imix1)%lDepl
      Fxr(ifxr, iz)%temp = mixture(imix1)%temp;       Fxr(ifxr, iz)%lh2o = mixture(imix1)%lh2o
      Fxr(ifxr, iz)%imix = imix1
      niso1 = mixture(imix1)%niso
      IF(Fxr(ifxr,iz)%ndim .LT. niso1+4) THEN
        DEALLOCATE(Fxr(ifxr,iz)%Pnum, Fxr(ifxr,iz)%idiso)
        CALL Dmalloc(Fxr(ifxr,iz)%pnum, niso1 + 4); CALL Dmalloc(Fxr(ifxr,iz)%IdIso, niso1 + 4)
        Fxr(ifxr,iz)%ndim = niso1 + 4
      ELSE
        CALL CP_CA(Fxr(ifxr, iz)%pnum(1:Fxr(ifxr,iz)%ndim), 0._8, Fxr(ifxr,iz)%ndim)
        CALL CP_CA(Fxr(ifxr, iz)%idiso(1:Fxr(ifxr,iz)%ndim), 0, Fxr(ifxr,iz)%ndim)
      ENDIF
      Fxr(ifxr,iz)%niso = Mixture(imix1)%niso
      
      CALL CP_VA(Fxr(ifxr,iz)%idiso(1:niso1), Mixture(imix1)%idiso(1:niso1), niso1)
      CALL CP_VA(Fxr(ifxr,iz)%pnum(1:niso1), Mixture(imix1)%pnum(1:niso1), niso1)
      IF(Fxr(ifxr, iz)%lH2O .AND. nTracerCntl%lInitBoron) THEN
        CALL   FxrBoronUpdate(Fxr(ifxr,iz), nTracerCntl%BoronPPM)
      ENDIF
      
      CONTINUE
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE SetBPCellInfo()
USE ioutil,    ONLY : Terminate
IMPLICIT NONE
INTEGER :: nRmvBP0
INTEGER :: i, j
INTEGER :: icel0, icel1, jcel0, jcel1
nRmvBP0 = nRmvBP
DO i = 1, nRMvBP0
  DO j = 1, 3
    icel0 = RmvBP(1, i); icel1 = RmvBP(2, i)
    IF(CellInfo(icel0)%EdgeCellIdx(j) .EQ. 0) CYCLE
    jcel0 = CellInfo(icel0)%EdgeCellIdx(j)
    jcel1 = CellInfo(icel1)%EdgeCellIdx(j)
    nRmvBP = nRmvBP + 1
    RmvBP(1, nRmvBP) = jcel0
    RmvBP(2, nRmvBP) = jcel1
    CONTINUE
  ENDDO
ENDDO

DO i = 1, nRmvBP
  icel0 = RmvBP(1, i); icel1 = RmvBP(2, i)
  IF(CellInfo(icel0)%nFxr .NE. CellInfo(icel1)%nFxr) THEN
    CALL TERMINATE('RMV_BP Error : Cell Type Mismatch')
  ENDIF
  IF(CellInfo(icel0)%nFsr .NE. CellInfo(icel1)%nFsr) THEN
    CALL TERMINATE('RMV_BP Error : Cell Type Mismatch')
  ENDIF

ENDDO

CONTINUE
END SUBROUTINE

SUBROUTINE ShfInfo(field, Shf, IdxSt)
USE IOUTIL,       ONLY : TOUPPER, TERMINATE
IMPLICIT NONE
CHARACTER(80) :: field
TYPE(Shf_Type) :: Shf
INTEGER :: IdxSt(2)

CHARACTER(10) :: dat1, dat2, dat3, dat4
INTEGER :: icyc, ix, iy, irot
INTEGER :: i

read(field, *) icyc, dat2, iy, irot
CALL TOUPPER(dat2) 
DO ix = 1, 26
   IF(dat2(1:1) .EQ. IDXX(ix)) EXIT
ENDDO
IF(ix .GT. 26) CALL TERMINATE('Wrong Shf Assembly Location')
Shf%ix = ix -Idxst(1) + 1 ; Shf%iy = iy -IdxSt(2) + 1
Shf%cycleid = icyc
Shf%irot = irot
Shf%lFreshFuel = .FALSE.
Shf%fields = Field
END SUBROUTINE

FUNCTION Xidx2No(Idx)
IMPLICIT NONE
INTEGER :: Xidx2No
CHARACTER(1) :: Idx
INTEGER :: ix 
DO ix = 1, 26
   IF(Idx(1:1) .EQ. IDXX(ix)) EXIT
ENDDO
Xidx2No = ix
END FUNCTION
END MODULE
