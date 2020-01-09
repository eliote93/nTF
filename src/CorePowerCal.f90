#include <defines.h>
SUBROUTINE CorePowerCal(Core, CmInfo, PowerDist, ng, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,   CMInfo_Type,   PowerDist_Type,     PE_TYPE,   &
                           AsyInfo_Type,    Asy_Type,      PinXs_Type
USE Files,          ONLY : io8
USE BasicOperation, ONLY : DotProduct
USE CNTL,           ONLY : nTracerCntl_Type
#ifdef MPI_ENV
USE MPIComm_mod,    ONLY : SENDRECV
#endif
USE ALLOCS
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CMInfo
TYPE(PowerDist_Type) :: PowerDist
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
INTEGER :: ng

TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(PinXS_Type), POINTER :: PinXs(:, :)
REAL, POINTER :: PhiC(:, :, :)
REAL, POINTER :: PinVol(:, :), AsyVol(:, :)
REAL, POINTER :: hz(:), hzfm(:)
INTEGER, POINTER :: SubPlaneMap(:), SubPlaneRange(:, :)
INTEGER :: ig, ixy, iz, iz1, iz2, iasy, ipin, icell
INTEGER :: i, j, k, ix, iy, jbeg, jend
INTEGER :: nx, ny, nxy, nxymax,nz, nxya, myzb, myze
INTEGER :: nAsyType, AsyType
REAL :: localpinpw, LocalAsyPw, pwsum, hactive
REAL :: fuelvol, fuelarea, vol, FuelAsyVol, FuelAsyArea
REAL :: Pin3DNormalizer, Pin2DNormalizer, Asy3DNormalizer, Asy2DNolmalizer, Axial1Dnormalizer
REAL :: Fm3DNormalizer, Fm2DNormalizer
LOGICAL :: master, slave

master = PE%Cmfdmaster; slave = PE%Cmfdslave
PhiC => CmInfo%PhiC; PinXs => CmInfo%PinXs
AsyInfo => Core%AsyInfo; Asy => Core%Asy
hz => Core%hz; hzfm => Core%hzfm
PinVol => Core%PinVol; AsyVol => Core%AsyVol
SubPlaneMap => Core%SubPlaneMap; SubPlaneRange => Core%SubPlaneRange

nAsyType = Core%nAsyType
nxya = Core%nxya; nz = Core%nz
myzb = PE%myzb; myze = PE%myze

nxy = 0
DO i = 1, nAsyType
  nxy = max(AsyInfo(i)%nxy, nxy)
ENDDO
nxymax = nxy
IF(Master) THEN
  ALLOCATE(PowerDist%PinPower2D(nxy, nxya))
  ALLOCATE(PowerDist%PinPower3D(nxy, nxya, nz))
  ALLOCATE(PowerDist%AsyPower2D(nxya))
  ALLOCATE(PowerDist%AsyPower3D(nz, nxya))
  ALLOCATE(PowerDist%Axial1DPower(nz))
ELSE
  ALLOCATE(PowerDist%PinPower3D(nxy, nxya, myzb:myze))
ENDIF

!Communicate

!--- CNJ Edit : MKL, CUDA Reorder Pin Cross Section

#ifdef __INTEL_MKL
IF (PE%lMKL) CALL MKL_ReorderPinXS(Core, CmInfo)
#endif

#ifdef __PGI
IF (PE%lCUDACMFD) CALL CUDAReorderPinXS(Core, CmInfo)
#endif

pwsum = 0
DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  nxy = AsyInfo(AsyType)%nxy
  Do iz = myzb, myze   !Axial Sweep
    IF(.NOT. Core%lFuelPlane(iz)) CYCLE
    DO i = 1, nxy  !Cell Sweep
      ixy = Asy(iasy)%GlobalPinIdx(i)
      LocalPinPw = DotProduct(PinXs(ixy, iz)%XSKF(1:ng), PhiC(ixy, iz, 1:ng), ng)
      LocalPinPw = LocalPinPw * PinVol(ixy, iz) !hz(iz)
      PowerDist%PinPower3D(i, iasy, iz) = LocalPinPw
      pwsum = pwsum + LocalPinPw 
      vol = vol + PinVol(ixy, iz)
    ENDDO
  ENDDO
ENDDO

#ifdef MPI_ENV 
IF(Master) THEN
  DO i = 1, PE%nCmfdProc - 1
    iz1 = PE%AxDomRange(1, i); iz2 = PE%AxDomRange(2, i)
    CALL SendRecv(PowerDist%PinPower3D(1:nxymax, 1:nxya, iz1:iz2), nxymax, nxya,    &
                  iz2 - iz1 + 1, i, .FALSE., PE%MPI_CMFD_COMM )
  ENDDO
ELSE
    iz1 = myzb; iz2 = myze
    CALL SendRecv(PowerDist%PinPower3D(1:nxymax, 1:nxya, iz1:iz2), nxymax, nxya,    &
                  iz2 - iz1 + 1, 0, .TRUE., PE%MPI_CMFD_COMM )
ENDIF

IF(slave) THEN
  NULLIFY(PhiC); NULLIFY(PinXs)
  NULLIFY(AsyInfo); NULLIFY(Asy)
  NULLIFY(hz); NULLIFY(PinVol)
  NULLIFY(AsyVol)
  RETURN
ENDIF
IF(master) THEN
  pwsum = 0; vol = 0
  DO iasy = 1, nxya
    AsyType = Asy(iasy)%AsyType
    IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
    nxy = AsyInfo(AsyType)%nxy
    Do iz = 1, nz   !Axial Sweep
      IF(.NOT. Core%lFuelPlane(iz)) CYCLE
      DO i = 1, nxy  !Cell Sweep
        ixy = Asy(iasy)%GlobalPinIdx(i)
        LocalPinPw = PowerDist%PinPower3D(i, iasy, iz)
        pwsum = pwsum + LocalPinPw 
        vol = vol + PinVol(ixy, iz)
      ENDDO
    ENDDO
  ENDDO !
ENDIF
#endif
!2D Pin Power
DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  nxy = AsyInfo(AsyType)%nxy
  PowerDist%PinPower2D(1:nxy, iasy) = 0
  Do iz = 1, nz   !Axial Sweep
    IF(.NOT. Core%lFuelPlane(iz)) CYCLE
    DO i = 1, nxy  !Cell Sweep
      PowerDist%PinPower2D(i, iasy) = PowerDist%PinPower2D(i, iasy) + PowerDist%PinPower3D(i, iasy, iz)
    ENDDO
  ENDDO
ENDDO !

!3D & 2D Assemby Power
DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  PowerDist%AsyPower3D(:, iasy) = 0._8; PowerDist%AsyPower2D(iasy) = 0._8
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  nxy = AsyInfo(AsyType)%nxy
  Do iz = 1, nz   !Axial Sweep
    IF(.NOT. Core%lFuelPlane(iz)) CYCLE
    LocalAsyPw = 0
    DO i = 1, nxy  !Cell Sweep
      LocalAsyPw = LocalAsyPw + PowerDist%PinPower3D(i, iasy, iz)
    ENDDO
    PowerDist%AsyPower3D(iz, iasy) = LocalAsyPw
  ENDDO
  PowerDist%AsyPower2D(iasy) = SUM(PowerDist%AsyPower3D(1:nz, iasy))
ENDDO 

!Axial Power Distribution
PowerDist%Axial1DPower = 0._8
DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  nxy = AsyInfo(AsyType)%nxy
  Do iz = 1, nz   !Axial Sweep
    IF(.NOT. Core%lFuelPlane(iz)) CYCLE
    PowerDist%Axial1DPower(iz) = PowerDist%Axial1DPower(iz) + PowerDist%AsyPower3D(iz, iasy)
  ENDDO
ENDDO

!Fuel Pin Volume
fuelvol = 0
DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  nxy = AsyInfo(AsyType)%nxy
  Do iz = 1, nz   !Axial Sweep
    IF(.NOT. Core%lFuelPlane(iz)) CYCLE
    DO i = 1, nxy  !Cell Sweep
      ixy = Asy(iasy)%GlobalPinIdx(i)
      LocalPinPw =PowerDist%PinPower3D(i, iasy, iz)
      IF(LocalPinPw .GT. 0) fuelvol = fuelvol + PinVol(ixy, iz)
    ENDDO
  ENDDO
ENDDO

!Fuel Pin area
fuelarea = 0
DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  nxy = AsyInfo(AsyType)%nxy
  DO i = 1, nxy  !Cell Sweep
    ixy = Asy(iasy)%GlobalPinIdx(i)
    LocalPinPw =PowerDist%PinPower2D(i, iasy)
    IF(LocalPinPw .GT. epsm6) fuelarea = fuelarea + PinVol(ixy, 1)/hz(1)
  ENDDO
ENDDO

!Fuel Asy Volume
FuelAsyVol = 0
DO iz =1, nz
  DO iasy = 1, nxya
    IF(PowerDist%AsyPower3D(iz, iasy) .gt. 0._8) FuelAsyVol = FuelAsyVol + AsyVol(iasy, iz)
  ENDDO
ENDDO

!Fuel Asy Area
FuelAsyArea = 0
DO iasy = 1, nxya
  IF(PowerDist%AsyPower2D(iasy) .gt. 0._8) FuelAsyArea = FuelAsyArea + AsyVol(iasy, 1) / hz(1)
ENDDO

!Active Core Height
hactive = 0
DO iz = 1, nz
  IF(Core%lFuelPlane(iz)) hactive = hactive + hz(iz)
ENDDO

Pin3DNormalizer = fuelvol / pwsum
Fm3DNormalizer = Core%FuelVolFM / pwsum
Pin2DNormalizer = Pin3DNormalizer/hactive
Fm2DNormalizer = Fm3DNormalizer / hactive
Asy3DNormalizer = FuelAsyVol / pwsum
Asy2DNolmalizer = Asy3DNormalizer / hactive
Axial1Dnormalizer = hactive / pwsum

IF(nTracerCntl%lProblem .EQ. lTransient) THEN
  Pin3DNormalizer = Pin3DNormalizer * nTracerCntl%PowerLevel
  Fm3DNormalizer = Fm3DNormalizer * nTracerCntl%PowerLevel
  Pin2DNormalizer = Pin2DNormalizer * nTracerCntl%PowerLevel
  Fm2DNormalizer = Fm2DNormalizer * nTracerCntl%PowerLevel
  Asy3DNormalizer = Asy3DNormalizer * nTracerCntl%PowerLevel
  Asy2DNolmalizer = Asy2DNolmalizer * nTracerCntl%PowerLevel
  Axial1Dnormalizer = Axial1Dnormalizer * nTracerCntl%PowerLevel
ENDIF

!Pin3DNormalizer = 1
!Pin2DNormalizer = 1
!Asy3DNormalizer = 1
!Asy2DNolmalizer = 1
!Axial1Dnormalizer = 1

nxy = Core%nxy

!3D Power
DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  nxy = AsyInfo(AsyType)%nxy
  Do iz = 1, nz   !Axial Sweep
    IF(.NOT. Core%lFuelPlane(iz)) CYCLE
    DO i = 1, nxy  !Cell Sweep
      ixy = Asy(iasy)%GlobalPinIdx(i)
      PowerDist%PinPower3D(i, iasy, iz) = PowerDist%PinPower3D(i, iasy, iz) * Pin3DNormalizer / PinVol(ixy, iz)
    ENDDO
  ENDDO
ENDDO 

!2D Power Normalization
DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  nxy = AsyInfo(AsyType)%nxy
  DO i = 1, nxy  !Cell Sweep
    ixy = Asy(iasy)%GlobalPinIdx(i)
    PowerDist%PinPower2D(i, iasy) = PowerDist%PinPower2D(i, iasy) * Pin2DNormalizer/(PinVol(ixy, 1)/hz(1))
  ENDDO
ENDDO 
CONTINUE

!3D Assembly Power Normalization
DO iasy = 1, nxya
  AsyType= Asy(iasy)%AsyType
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  DO iz = 1, nz
    IF(.NOT. Core%lFuelPlane(iz)) CYCLE
    PowerDist%AsyPower3D(iz, iasy) = PowerDist%AsyPower3D(iz, iasy) * Asy3DNormalizer / (AsyVol(iasy, iz))
  ENDDO
ENDDO

DO iasy = 1, nxya
  AsyType= Asy(iasy)%AsyType
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  PowerDist%AsyPower2D(iasy) = PowerDist%AsyPower2D(iasy) * Asy2DNolmalizer * (hz(1)/AsyVol(iasy, 1))
ENDDO

DO iz = 1, nz
  PowerDist%Axial1DPower(iz) = PowerDist%Axial1DPower(iz) * Axial1Dnormalizer/ hz(iz)
ENDDO

PowerDist%Pin3DNormalizer = Pin3DNormalizer
PowerDist%Pin2DNormalizer = Pin2DNormalizer
PowerDist%Fm3DNormalizer = Fm3DNormalizer
PowerDist%Fm2DNormalizer = Fm2DNormalizer
PowerDist%Asy3DNormalizer = Asy3DNormalizer
PowerDist%Asy2DNormalizer = Asy2DNolmalizer
PowerDist%Axial1Dnormalizer = Axial1Dnormalizer

NULLIFY(PhiC); NULLIFY(PinXs)
NULLIFY(AsyInfo); NULLIFY(Asy)
NULLIFY(hz); NULLIFY(PinVol)
NULLIFY(AsyVol)

END SUBROUTINE