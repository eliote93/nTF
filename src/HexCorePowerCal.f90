#include <defines.h>
SUBROUTINE HexCorePowerCal(Core, CmInfo, PowerDist, ng, nTracerCntl, PE)

USE PARAM,          ONLY : ZERO, lTransient, FALSE, TRUE
USE ioutil,         ONLY : terminate
USE TYPEDEF,        ONLY : CoreInfo_Type,   CMInfo_Type,   PowerDist_Type,     PE_TYPE,   &
                           AsyInfo_Type,    Asy_Type,      PinXs_Type,         FxrInfo_Type, &
                           Cell_Type,       Pin_Type
USE Files,          ONLY : io8
USE BasicOperation, ONLY : DotProduct
USE CNTL,           ONLY : nTracerCntl_Type
USE ALLOCS

USE Core_mod,       ONLY : FmInfo
USE HexData,        ONLY : hAsy, hPinInfo
!USE HexTst,         ONLY : HexPrtPhiC, HexPrtXsKF

#ifdef MPI_ENV
USE MPIComm_mod,    ONLY : SENDRECV
#endif

#ifdef __INTEL_MKL
USE MKL_3D
USE CMFD_COMMON,    ONLY : HomogenizeXS, AllocHomoXSVar
#endif

IMPLICIT NONE

TYPE(CoreInfo_Type)    :: Core
TYPE(CmInfo_Type)      :: CMInfo
TYPE(PowerDist_Type)   :: PowerDist
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type)          :: PE

INTEGER :: ng

TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type),     POINTER :: Asy(:)
TYPE(PinXS_Type),   POINTER :: PinXs(:, :)
TYPE(cell_type),    POINTER :: Cell(:)
TYPE(Fxrinfo_type), POINTER :: Fxr(:,:)
TYPE(Pin_Type),     POINTER :: Pin(:)

REAL, POINTER :: PhiC(:, :, :), PhiS(:,:,:)
REAL, POINTER :: AsyVol(:, :)
REAL, POINTER :: hz(:)

INTEGER :: j, ig, ixy, jxy, kxy, iz, iz1, iz2, iasy, ipin, iCel, ifsr, jfsr
INTEGER :: nxymax, nz, nxya, myzb, myze
INTEGER :: nAsyType, AsyType, nFxrMax, FsrIdxSt, nLocalFsr

REAL :: localpinpw, LocalAsyPw, pwsum, hactive
REAL :: fuelvol, vol, FuelAsyVol, phisum, volsum
REAL :: Pin3DNormalizer, Pin2DNormalizer, Asy3DNormalizer, Asy2DNolmalizer, Axial1Dnormalizer
REAL :: Fm3DNormalizer, Fm2DNormalizer

LOGICAL :: master, slave, lCMFD, lxslib, lscat1
! ----------------------------------------------------

master   = PE%Cmfdmaster
slave    = PE%Cmfdslave
AsyInfo => Core%AsyInfo
Cell    => Core%CellInfo
Pin     => Core%Pin
Asy     => Core%Asy
hz      => Core%hz
AsyVol  => Core%AsyVol
phis    => FmInfo%phis
Fxr     => FmInfo%Fxr

nAsyType = Core%nAsyType
nxya     = Core%nxya
nz       = Core%nz
myzb     = PE%myzb
myze     = PE%myze
lCMFD    = nTracerCntl%lCMFD
lxslib   = nTracerCntl%lxslib
lscat1   = nTracerCntl%lscat1

nxymax = 0 ! Max Number of Fuel Pin

DO iAsy = 1, nAsyType
  nxymax = max(AsyInfo(iAsy)%nxy, nxymax)
END DO

IF(Master) THEN
  ALLOCATE (PowerDist%PinPower2D   (nxymax, nxya))
  ALLOCATE (PowerDist%PinPower3D   (nxymax, nxya, nz))
  ALLOCATE (PowerDist%AsyPower2D   (nxya))
  ALLOCATE (PowerDist%AsyPower3D   (nz, nxya))
  ALLOCATE (PowerDist%Axial1DPower (nz))
  
  PowerDist%PinPower2D   = ZERO
  PowerDist%PinPower3D   = ZERO
  PowerDist%AsyPower2D   = ZERO
  PowerDist%AsyPower3D   = ZERO
  PowerDist%Axial1DPower = ZERO
ELSE
  ALLOCATE(PowerDist%PinPower3D (nxymax, nxya, myzb:myze))
  PowerDist%PinPower3D = ZERO
END IF
! ----------------------------------------------------
!               01. SET : Cm Info % Pin XS
! ----------------------------------------------------
#ifdef __INTEL_MKL
IF (PE%lMKL .AND. .NOT. lCMFD) THEN
  CALL AllocCMFD(TRUE)
  CALL AllocHomoXSVar(Core, mklGeom%ng)
  CALL HomogenizeXS(Core, mklGeom%superPin, Fxr, mklCMFD%PinXS, phis, ng, &
                    mklGeom%nxy, mklGeom%myzb, mklGeom%myze, lxslib, lscat1, FALSE)
END IF

IF (PE%lMKL) CALL MKL_ReorderPinXS(Core, CmInfo)

IF (PE%lMKL .AND. .NOT. lCMFD) THEN
  DO ig = 1, ng
    DO iz = myzb, myze
      DO ixy = 1, mklGeom%nxy
        phisum = ZERO
        volsum = ZERO
        
        DO jxy = 1, mklGeom%superPin(ixy)%nxy
          kxy       = mklGeom%superPin(ixy)%pin(jxy)
          FsrIdxSt  = Pin(kxy)%FsrIdxSt
          icel      = Pin(kxy)%Cell(iz)
          nLocalFsr = Cell(icel)%nFsr
          
          DO ifsr = 1, nLocalFsr
            jfsr = FsrIdxSt + ifsr - 1
            vol  = Cell(icel)%vol(ifsr)
            
            volsum = volsum + vol
            phisum = phisum + vol * phis(jfsr, iz, ig)
          END DO
        END DO
        
        DO jxy = 1, mklGeom%superPin(ixy)%nxy
          kxy = mklGeom%superPin(ixy)%pin(jxy)
          
          CmInfo%PhiC(kxy, iz, ig) = phisum / volsum
        END DO
      END DO
    END DO
  END DO
END IF
#endif

#ifdef __PGI
IF (PE%lCUDACMFD) CALL CUDAReorderPinXS(Core, CmInfo)
#endif

PhiC  => CmInfo%PhiC
PinXs => CmInfo%PinXs
! ----------------------------------------------------
!               02. SET : 3D Pin Power
! ----------------------------------------------------
pwsum = ZERO

DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  
  DO iz = myzb, myze   !Axial Sweep
    IF(.NOT. Core%lFuelPlane(iz)) CYCLE
    
    DO ixy = 1, hAsy(iAsy)%nTotPin  !Pin Sweep
      jxy  = hAsy(iAsy)%PinIdxSt - 1 + ixy
      iPin = hPinInfo(jxy)%OrdInAsy01
      
      LocalPinPw = DotProduct(PinXs(jxy, iz)%XsKF(1:ng), PhiC(jxy, iz, 1:ng), ng)
      LocalPinPw = LocalPinPw * hPinInfo(jxy)%Vol(iz)
      
      IF (LocalPinPw .NE. LocalPinPw) CALL terminate("NAN POWER")
      
      PowerDist%PinPower3D(iPin, iasy, iz) = LocalPinPw
      
      pwsum = pwsum + LocalPinPw
      vol   = vol   + hPinInfo(jxy)%Vol(iz)
    END DO
  END DO
END DO
! ----------------------------------------------------
!               03. MERGE in MPI
! ----------------------------------------------------
#ifdef MPI_ENV 
IF(Master) THEN
  DO j = 1, PE%nCmfdProc - 1
    iz1 = PE%AxDomRange(1, j)
    iz2 = PE%AxDomRange(2, j)
    
    CALL SendRecv(PowerDist%PinPower3D(1:nxymax, 1:nxya, iz1:iz2), nxymax, nxya, iz2 - iz1 + 1, j, FALSE, PE%MPI_CMFD_COMM )
  END DO
ELSE
    iz1 = myzb
    iz2 = myze
    CALL SendRecv(PowerDist%PinPower3D(1:nxymax, 1:nxya, iz1:iz2), nxymax, nxya, iz2 - iz1 + 1, 0, TRUE, PE%MPI_CMFD_COMM )
END IF

IF(slave) THEN
  NULLIFY (PhiC, PinXs, AsyInfo, Asy, hz, AsyVol)
  
  RETURN
END IF

IF(Master) THEN
  pwsum = ZERO
  vol   = ZERO
  
  DO iasy = 1, nxya
    AsyType = Asy(iasy)%AsyType
    
    IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
    
    DO iz = 1, nz   !Axial Sweep
      IF(.NOT. Core%lFuelPlane(iz)) CYCLE
      
      DO ixy = 1, hAsy(iAsy)%nRodPin  !Cell Sweep
        jxy  = hAsy(iAsy)%PinIdxSt - 1 + ixy
        iPin = hPinInfo(jxy)%OrdInAsy01
        
        LocalPinPw = PowerDist%PinPower3D(iPin, iAsy, iz)
        
        pwsum = pwsum + LocalPinPw
        vol   = vol + hPinInfo(jxy)%Vol(iz)
      END DO
    END DO
  END DO
END IF
#endif
! ----------------------------------------------------
!               04. SET : 2D Pin Power
! ----------------------------------------------------
DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  
  DO iz = 1, nz   !Axial Sweep
    IF(.NOT. Core%lFuelPlane(iz)) CYCLE
    
    DO iPin = 1, AsyInfo(AsyType)%nxy  !Pin Sweep
      PowerDist%PinPower2D(iPin, iasy) = PowerDist%PinPower2D(iPin, iasy) + PowerDist%PinPower3D(iPin, iasy, iz)
    END DO
  END DO
END DO
! ----------------------------------------------------
!               05. SET : 3D & 2D Asy Power
! ----------------------------------------------------
DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  
  DO iz = 1, nz   !Axial Sweep
    IF(.NOT. Core%lFuelPlane(iz)) CYCLE
    
    LocalAsyPw = ZERO
    
    DO iPin = 1, AsyInfo(AsyType)%nxy  !Pin Sweep
      LocalAsyPw = LocalAsyPw + PowerDist%PinPower3D(iPin, iasy, iz)
    END DO
    
    PowerDist%AsyPower3D(iz, iasy) = LocalAsyPw
  END DO
  
  PowerDist%AsyPower2D(iasy) = SUM(PowerDist%AsyPower3D(1:nz, iasy))
END DO 
! ----------------------------------------------------
!               06. Axial Power Distribution
! ----------------------------------------------------
DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  
  DO iz = 1, nz   !Axial Sweep
    IF(.NOT. Core%lFuelPlane(iz)) CYCLE
    PowerDist%Axial1DPower(iz) = PowerDist%Axial1DPower(iz) + PowerDist%AsyPower3D(iz, iasy)
  END DO
END DO
! ----------------------------------------------------
!               07. SUM : Fuel Pin Volume
! ----------------------------------------------------
FuelVol = ZERO

DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  
  DO iz = 1, nz   !Axial Sweep
    IF(.NOT. Core%lFuelPlane(iz)) CYCLE
    
    DO ixy = 1, hAsy(iAsy)%nRodPin  !Cell Sweep
      jxy  = hAsy(iAsy)%PinIdxSt - 1 + ixy
      iPin = hPinInfo(jxy)%OrdInAsy01
      
      LocalPinPw = PowerDist%PinPower3D(iPin, iasy, iz)
      
      IF(LocalPinPw .GT. 1E-15) THEN
        FuelVol = FuelVol + hPinInfo(jxy)%Vol(iz)
      END IF
    END DO
  END DO
END DO
! ----------------------------------------------------
!               08. SUM : Fuel Asy Volume
! ----------------------------------------------------
FuelAsyVol = ZERO

DO iz =1, nz
  DO iasy = 1, nxya
    IF(PowerDist%AsyPower3D(iz, iasy) .gt. 1E-15) FuelAsyVol = FuelAsyVol + AsyVol(iasy, iz)
  END DO
END DO
! ----------------------------------------------------
!               09. Active Core Height
! ----------------------------------------------------
hActive = ZERO

DO iz = 1, nz
  IF(Core%lFuelPlane(iz)) hactive = hactive + hz(iz)
END DO

Pin3DNormalizer   = FuelVol         / pwsum
Fm3DNormalizer    = Core%FuelVolFM  / pwsum
Pin2DNormalizer   = Pin3DNormalizer / hactive
Fm2DNormalizer    = Fm3DNormalizer  / hactive
Asy3DNormalizer   = FuelAsyVol      / pwsum
Asy2DNolmalizer   = Asy3DNormalizer / hactive
Axial1Dnormalizer = hactive         / pwsum

IF(nTracerCntl%lProblem .EQ. lTransient) THEN
  Pin3DNormalizer   = Pin3DNormalizer   * nTracerCntl%PowerLevel
  Fm3DNormalizer    = Fm3DNormalizer    * nTracerCntl%PowerLevel
  Pin2DNormalizer   = Pin2DNormalizer   * nTracerCntl%PowerLevel
  Fm2DNormalizer    = Fm2DNormalizer    * nTracerCntl%PowerLevel
  Asy3DNormalizer   = Asy3DNormalizer   * nTracerCntl%PowerLevel
  Asy2DNolmalizer   = Asy2DNolmalizer   * nTracerCntl%PowerLevel
  Axial1Dnormalizer = Axial1Dnormalizer * nTracerCntl%PowerLevel
END IF
! ----------------------------------------------------
!               10. NORM : 3D Pin Power
! ----------------------------------------------------
DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  
  DO iz = 1, nz   !Axial Sweep
    IF(.NOT. Core%lFuelPlane(iz)) CYCLE
    
    DO ixy = 1, hAsy(iAsy)%nRodPin  !Pin Sweep
      jxy  = hAsy(iAsy)%PinIdxSt - 1 + ixy
      iPin = hPinInfo(jxy)%OrdInAsy01
      
      PowerDist%PinPower3D(iPin, iasy, iz) = PowerDist%PinPower3D(iPin, iasy, iz) * Pin3DNormalizer &
                                           / hPinInfo(jxy)%Vol(iz)
    END DO
  END DO
END DO
! ----------------------------------------------------
!               11. NORM : 2D Pin Power
! ----------------------------------------------------
DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  
  DO ixy = 1, hAsy(iAsy)%nRodPin  !Cell Sweep
    jxy  = hAsy(iAsy)%PinIdxSt - 1 + ixy
    iPin = hPinInfo(jxy)%OrdInAsy01

    PowerDist%PinPower2D(iPin, iasy) = PowerDist%PinPower2D(iPin, iasy) * Pin2DNormalizer &
                                     / (hPinInfo(jxy)%Vol(1) / hz(1))
  END DO
END DO 

CONTINUE
! ----------------------------------------------------
!               12. NORM : 3D Assy Power
! ----------------------------------------------------
DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  
  DO iz = 1, nz
    IF(.NOT. Core%lFuelPlane(iz)) CYCLE
    
    PowerDist%AsyPower3D(iz, iasy) = PowerDist%AsyPower3D(iz, iasy) * Asy3DNormalizer / (AsyVol(iasy, 1))
  END DO
END DO
! ----------------------------------------------------
!               13. NORM : 2D Assy Power 
! ----------------------------------------------------
DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  
  PowerDist%AsyPower2D(iasy) = PowerDist%AsyPower2D(iasy) * Asy2DNolmalizer * (hz(1)/AsyVol(iasy, 1))
END DO

DO iz = 1, nz
  PowerDist%Axial1DPower(iz) = PowerDist%Axial1DPower(iz) * Axial1Dnormalizer/ hz(iz)
END DO
! ----------------------------------------------------
!               14. Move
! ----------------------------------------------------
PowerDist%Pin3DNormalizer   = Pin3DNormalizer
PowerDist%Pin2DNormalizer   = Pin2DNormalizer
PowerDist%Fm3DNormalizer    = Fm3DNormalizer
PowerDist%Fm2DNormalizer    = Fm2DNormalizer
PowerDist%Asy3DNormalizer   = Asy3DNormalizer
PowerDist%Asy2DNormalizer   = Asy2DNolmalizer
PowerDist%Axial1Dnormalizer = Axial1Dnormalizer

NULLIFY (PhiC, PinXs, AsyInfo, Asy, hz, AsyVol)
! ----------------------------------------------------

END SUBROUTINE HexCorePowerCal
! ------------------------------------------------------------------------------------------------------------