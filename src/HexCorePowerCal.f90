#include <defines.h>
SUBROUTINE HexCorePowerCal(Core, CmInfo, PowerDist, ng, nTracerCntl, PE)

USE allocs
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

IF (Master) THEN
  CALL dmalloc(PowerDist%PinPower2D,   nxymax, nxya)
  CALL dmalloc(PowerDist%PinPower3D,   nxymax, nxya, nz)
  CALL dmalloc(PowerDist%AsyPower2D,   nxya)
  CALL dmalloc(PowerDist%AsyPower3D,   nz, nxya)
  CALL dmalloc(PowerDist%Axial1DPower, nz)
ELSE
  CALL dmalloc0(PowerDist%PinPower3D, 1, nxymax, 1, nxya, myzb, myze)
END IF
! ----------------------------------------------------
!               01. SET : Cm Info % Pin XS
! ----------------------------------------------------
#ifdef __INTEL_MKL
IF (PE%lMKL .AND. .NOT. lCMFD) THEN
  CALL AllocCMFD(TRUE)
  CALL AllocHomoXSVar(Core, mklGeom%ng)
  CALL HomogenizeXS(Core, mklGeom%superPin, Fxr, mklCMFD%PinXS, phis, ng, mklGeom%nxy, mklGeom%myzb, mklGeom%myze, lxslib, lscat1, FALSE)
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

SUBROUTINE HexCalExplicitPower(Core, CmInfo, FmInfo, PhotonPower, NeutronPower, TotalExpPower, ng, ngg, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,        ONLY : coreinfo_type, PowerDist_Type, PE_TYPE, ASYINFO_TYPE,  &
                           Asy_Type,      FMInfo_TYPE,    Fxrinfo_type
USE files,          ONLY : io8
USE BasicOperation, ONLY : DotProduct 
USE CNTL,           ONLY : nTracerCntl_Type
#ifdef MPI_ENV
USE MPIComm_Mod,    ONLY : SENDRECV
#endif
USE allocs
USE GammaTYPEDEF,   ONLY : PINHEATXS_TYPE,  GammaCMFD_Type
USE GammaCore_mod,  ONLY : PinHeatXS, gphis
USE HexData,        ONLY : hAsy,      hPinInfo
IMPLICIT NONE
TYPE(coreinfo_type) :: Core 
TYPE(GammaCMFD_Type) :: CmInfo
TYPE(FMInfo_TYPE) :: FmInfo
TYPE(PowerDist_Type) :: PhotonPower, NeutronPower, TotalExpPower
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
INTEGER :: ng, ngg

TYPE(ASYINFO_TYPE), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(Fxrinfo_type), POINTER :: FXR(:, :)
REAL, POINTER :: NPhiC(:,:,:), GPhiC(:,:,:)
REAL, POINTER :: phis(:,:,:)
REAL :: localpinpw, localasypw, pwsum, hactive
REAL :: fuelvol, fuelarea, vol, fuelasyvol, fuelasyarea
REAL :: PGensum, Gpwsum, Npwsum
INTEGER :: nAsyType
INTEGER :: nxy, nxymax, myzb, myze, nxya, nz, i, iz1, iz2
INTEGER :: iasy, AsyType, iz, ixy, jxy, ipin
LOGICAL :: master, slave
LOGICAL :: lXsLib, lScat1

AsyInfo => Core%AsyInfo
Asy => Core%Asy
FXR => FmInfo%Fxr
phis => FmInfo%Phis
NPhiC => CmInfo%Phic
GPhiC => CmInfo%gPhic

nAsyType = Core%nAsyType
nxya = Core%nxya
nz = Core%nz
myzb = PE%myzb
myze = PE%myze
master = PE%CMFDMaster
slave = PE%CMFDSlave
lXsLib = nTracerCntl%lXsLib
lScat1 = nTracerCntl%lGammaScat1

CALL GenHomoKERMA(Core, FXR, Phis, GPhis, PinHeatXS, myzb, myze, ng, ngg, lxslib, lScat1)

nxymax = 0
DO iasy = 1, nAsyType
  nxymax = max(AsyInfo(iAsy)%nxy, nxymax)
END DO 
nxy = nxymax

CALL UpdateNGPHIC(PinHeatXS, NPhiC, GPhiC, myzb, myze, Core%nxy, ng, ngg)

IF(master) THEN 
  ALLOCATE(PhotonPower%PinPower2D(nxy, nxya))
  ALLOCATE(PhotonPower%PinPower3D(nxy, nxya, nz))
  ALLOCATE(PhotonPower%AsyPower2D(nxya))
  ALLOCATE(PhotonPower%AsyPower3D(nz, nxya))
  ALLOCATE(PhotonPower%Axial1DPower(nz))
  
  ALLOCATE(NeutronPower%PinPower2D(nxy, nxya))
  ALLOCATE(NeutronPower%PinPower3D(nxy, nxya, nz))
  ALLOCATE(NeutronPower%AsyPower2D(nxya))
  ALLOCATE(NeutronPower%AsyPower3D(nz, nxya))
  ALLOCATE(NeutronPower%Axial1DPower(nz))
  
  ALLOCATE(TotalExpPower%PinPower2D(nxy, nxya))
  ALLOCATE(TotalExpPower%PinPower3D(nxy, nxya, nz))
  ALLOCATE(TotalExpPower%AsyPower2D(nxya))
  ALLOCATE(TotalExpPower%AsyPower3D(nz, nxya))
  ALLOCATE(TotalExpPower%Axial1DPower(nz))
ELSE 
  ALLOCATE(PhotonPower%PinPower3D(nxy, nxya, myzb:myze))
  ALLOCATE(NeutronPower%PinPower3D(nxy, nxya, myzb:myze))
  ALLOCATE(TotalExpPower%PinPower3D(nxy, nxya, myzb:myze))
END IF

pwsum = 0
PGensum = 0
Gpwsum = 0
Npwsum = 0
DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  DO iz = myzb, myze
    DO ixy = 1, hAsy(iasy)%nTotPin
      jxy = hAsy(iasy)%PinIdxSt - 1 + ixy
      ipin = hPinInfo(jxy)%OrdInAsy01

      localpinpw = DotProduct(PinHeatxs(jxy,iz)%KERMA_T(1:ng), NPhiC(jxy, iz, 1:ng), ng) * hPinInfo(jxy)%Vol(iz)
      NeutronPower%PinPower3D(ipin, iasy, iz) = localpinpw
      Npwsum = Npwsum + localpinpw

      localpinpw = DotProduct(PinHeatXs(jxy,iz)%KERMA_P(1:ngg), GPhiC(jxy, iz, 1:ngg), ngg) * hPinInfo(jxy)%Vol(iz)
      PhotonPower%PinPower3D(ipin, iasy, iz) = localpinpw
      Gpwsum = Gpwsum + localpinpw

      localpinpw = DotProduct(PinHeatXs(jxy,iz)%PhotonGen(1:ng), NPhiC(jxy,iz,1:ng), ng) * hPinInfo(jxy)%Vol(iz)
      PGensum = PGensum + localpinpw

      TotalExpPower%PinPower3D(ipin, iasy, iz) = NeutronPower%PinPower3D(ipin, iasy, iz) + PhotonPower%PinPower3D(ipin, iasy, iz)
      pwsum  = pwsum + TotalExpPower%PinPower3D(ipin, iasy, iz)
    END DO 
  END DO 
END DO 

#ifdef MPI_ENV 
IF(Master) THEN
  DO i = 1, PE%nCmfdProc - 1
    iz1 = PE%AxDomRange(1, i); iz2 = PE%AxDomRange(2, i)
    CALL SendRecv(NeutronPower%PinPower3D(1:nxymax, 1:nxya, iz1:iz2), nxymax, nxya,    &
                  iz2 - iz1 + 1, i, .FALSE., PE%MPI_CMFD_COMM )
    CALL SendRecv(PhotonPower%PinPower3D(1:nxymax, 1:nxya, iz1:iz2), nxymax, nxya,    &
                  iz2 - iz1 + 1, i, .FALSE., PE%MPI_CMFD_COMM )
    CALL SendRecv(TotalExpPower%PinPower3D(1:nxymax, 1:nxya, iz1:iz2), nxymax, nxya,    &
                  iz2 - iz1 + 1, i, .FALSE., PE%MPI_CMFD_COMM )
  ENDDO
ELSE
    iz1 = myzb; iz2 = myze
    CALL SendRecv(NeutronPower%PinPower3D(1:nxymax, 1:nxya, iz1:iz2), nxymax, nxya,    &
                  iz2 - iz1 + 1, 0, .TRUE., PE%MPI_CMFD_COMM )
    CALL SendRecv(PhotonPower%PinPower3D(1:nxymax, 1:nxya, iz1:iz2), nxymax, nxya,    &
                  iz2 - iz1 + 1, 0, .TRUE., PE%MPI_CMFD_COMM )
    CALL SendRecv(TotalExpPower%PinPower3D(1:nxymax, 1:nxya, iz1:iz2), nxymax, nxya,    &
                  iz2 - iz1 + 1, 0, .TRUE., PE%MPI_CMFD_COMM )
ENDIF

IF(slave) RETURN

IF(master) THEN
  pwsum = 0; vol = 0
  Gpwsum = 0.; Npwsum = 0.
  DO iasy = 1, nxya
    AsyType = Asy(iasy)%AsyType
    Do iz = 1, nz   !Axial Sweep
      DO ixy = 1, hasy(iasy)%nTotPin  !Cell Sweep
        jxy = hAsy(iasy)%PinIdxSt -1 + ixy
        iPin = hPinInfo(jxy)%OrdInAsy01

        Npwsum = Npwsum +  NeutronPower%PinPower3D(ipin, iasy, iz)
        Gpwsum = Gpwsum + PhotonPower%PinPower3D(ipin, iasy, iz) 
        pwsum = pwsum + TotalExpPower%PinPower3D(ipin, iasy, iz)
        vol = vol + hPinInfo(jxy)%vol(iz)
      ENDDO
    ENDDO
  ENDDO !
ENDIF
#endif

NeutronPower%pwsum = Npwsum
PhotonPower%pwsum = Gpwsum
TotalExpPower%pwsum = pwsum

CALL NormPowerDist(Core, NeutronPower, CmInfo, ng, ngg, nxy, nTracerCntl, PE)
CALL NormPowerDist(Core, PhotonPower, CmInfo, ng, ngg, nxy, nTracerCntl, PE)
CALL NormPowerDist(Core, TotalExpPower, CmInfo, ng, ngg, nxy, nTracerCntl, PE)



contains

SUBROUTINE NormPowerDist(Core, PowerDist, CmInfo, ng, ngg, nxy, nTracerCntl, PE)!, lNormalization
  USE TYPEDEF,        ONLY : CoreInfo_type, PowerDist_Type, CMINFO_TYPE, PE_TYPE,  AsyInfo_Type,    Asy_Type
  USE CNTL,           ONLY : nTRACERCntl_TYPE
  USE GammaTYPEDEF,   ONLY : GammaCMFD_Type
  IMPLICIT NONE
  ! INPUT/OUTPUT VARIABLES
  TYPE(CoreInfo_TYPE) :: Core
  TYPE(PowerDist_TYPE) :: PowerDist
  TYPE(GammaCMFD_Type) :: CmInfo
  INTEGER :: ng, ngg
  REAL :: pwsum
  TYPE(nTRACERCntl_TYPE) :: nTRACERCntl
  TYPE(PE_TYPE) :: PE
  ! LOCAL VARIABLES..
  TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
  TYPE(Asy_Type), POINTER :: Asy(:)
  REAL, POINTER :: PinVol(:, :), AsyVol(:, :)
  REAL, POINTER :: hz(:)
  INTEGER :: ixy, iz, iasy, jxy, ipin
  INTEGER :: i, j, k
  INTEGER :: nx, ny, nxy,nz, nxya
  INTEGER :: AsyType
  REAL :: localpinpw, LocalAsyPw, hactive
  REAL :: fuelvol, fuelarea, vol, FuelAsyVol, FuelAsyArea
  REAL :: Pin3DNormalizer, Pin2DNormalizer, Asy3DNormalizer, Asy2DNolmalizer, Axial1Dnormalizer
  REAL :: Fm3DNormalizer, Fm2DNormalizer
  
  AsyInfo => Core%AsyInfo; Asy => Core%Asy
  hz => Core%hz
  PinVol => Core%PinVol; AsyVol => Core%AsyVol
  nxya = Core%nxya; nz = Core%nz
  pwsum = PowerDist%pwsum
  !2D Pin Power
  DO iasy = 1, nxya
    AsyType = Asy(iasy)%AsyType
    PowerDist%PinPower2D(1:hAsy(iasy)%nTotPin, iasy) = 0
    Do iz = 1, nz   !Axial Sweep
      DO ixy = 1, hAsy(iasy)%nTotPin  !Cell Sweep
        jxy = hAsy(iasy)%PinIdxSt - 1 + ixy
        ipin = hPinInfo(jxy)%OrdInAsy01
        PowerDist%PinPower2D(ipin, iasy) = PowerDist%PinPower2D(ipin, iasy) + PowerDist%PinPower3D(ipin, iasy, iz)
      ENDDO
    ENDDO
  ENDDO !
  
  !3D & 2D Assemby Power
  DO iasy = 1, nxya
    AsyType = Asy(iasy)%AsyType
    PowerDist%AsyPower3D(:, iasy) = 0._8; PowerDist%AsyPower2D(iasy) = 0._8
    Do iz = 1, nz   !Axial Sweep
      LocalAsyPw = 0
      DO ixy = 1, hAsy(iasy)%nTotPin  !Cell Sweep
        jxy = hAsy(iasy)%PinIdxSt - 1 + ixy
        ipin = hPinInfo(jxy)%OrdInAsy01
        LocalAsyPw = LocalAsyPw + PowerDist%PinPower3D(ipin, iasy, iz)
      ENDDO
      PowerDist%AsyPower3D(iz, iasy) = LocalAsyPw
    ENDDO
    PowerDist%AsyPower2D(iasy) = SUM(PowerDist%AsyPower3D(1:nz, iasy))
  ENDDO 
  
  !Axial Power Distribution
  PowerDist%Axial1DPower = 0._8
  DO iasy = 1, nxya
    AsyType = Asy(iasy)%AsyType
    Do iz = 1, nz   !Axial Sweep
      PowerDist%Axial1DPower(iz) = PowerDist%Axial1DPower(iz) + PowerDist%AsyPower3D(iz, iasy)
    ENDDO
  ENDDO
  
  !Fuel Pin Volume
  fuelvol = 0
  DO iasy = 1, nxya
    AsyType = Asy(iasy)%AsyType
    Do iz = 1, nz   !Axial Sweep
      DO ixy = 1, hasy(iasy)%nTotPin  !Cell Sweep
        jxy = hAsy(iasy)%PinIdxSt - 1 + ixy
        ipin = hPinInfo(jxy)%OrdInAsy01
        LocalPinPw = PowerDist%PinPower3D(ipin, iasy, iz)
        fuelvol = fuelvol + hPinInfo(jxy)%vol(iz)
      ENDDO
    ENDDO
  ENDDO
  
  !Fuel Pin area
  fuelarea = 0
  DO iasy = 1, nxya
    AsyType = Asy(iasy)%AsyType
    IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
    DO ixy = 1, hasy(iasy)%nTotPin  !Cell Sweep
      jxy = hAsy(iasy)%PinIdxSt - 1 + ixy
      ipin = hPinInfo(jxy)%OrdInAsy01
      LocalPinPw =PowerDist%PinPower2D(ipin, iasy)
      fuelarea = fuelarea + hPinInfo(jxy)%vol(1)/hz(1)
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
    !IF(Core%lFuelPlane(iz)) hactive = hactive + hz(iz)
     hactive = hactive + hz(iz)
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
  PowerDist%Pin3DNormalizer = Pin3DNormalizer
  PowerDist%Pin2DNormalizer = Pin2DNormalizer
  PowerDist%Fm3DNormalizer = Fm3DNormalizer
  PowerDist%Fm2DNormalizer = Fm2DNormalizer
  PowerDist%Asy3DNormalizer = Asy3DNormalizer
  PowerDist%Asy2DNormalizer = Asy2DNolmalizer
  PowerDist%Axial1Dnormalizer = Axial1Dnormalizer
  
  !! Normalization
  !IF (.NOT.lNormalization) RETURN
  nxy = Core%nxy
  !3D Power
  DO iasy = 1, nxya
    AsyType = Asy(iasy)%AsyType
    Do iz = 1, nz   !Axial Sweep
      DO ixy = 1, hAsy(iasy)%nTotPin  !Cell Sweep
        jxy = hAsy(iasy)%PinIdxSt - 1 + ixy
        ipin = hPinInfo(jxy)%OrdInAsy01
        PowerDist%PinPower3D(ipin, iasy, iz) = PowerDist%PinPower3D(ipin, iasy, iz) * Pin3DNormalizer / hPinInfo(jxy)%vol(iz)
      ENDDO
    ENDDO
  ENDDO 
  
  !2D Power Normalization
  DO iasy = 1, nxya
    AsyType = Asy(iasy)%AsyType
    DO ixy = 1, hAsy(iasy)%nTotPin  !Cell Sweep
      jxy = hAsy(iasy)%PinIdxSt - 1 + ixy
      ipin = hPinInfo(jxy)%OrdInAsy01
      PowerDist%PinPower2D(ipin, iasy) = PowerDist%PinPower2D(ipin, iasy) * Pin2DNormalizer/(hPinInfo(jxy)%vol(1)/hz(1))
    ENDDO
  ENDDO 
  CONTINUE
  
  !3D Assembly Power Normalization
  DO iasy = 1, nxya
    AsyType= Asy(iasy)%AsyType
    DO iz = 1, nz
      PowerDist%AsyPower3D(iz, iasy) = PowerDist%AsyPower3D(iz, iasy) * Asy3DNormalizer / (AsyVol(iasy, iz))
    ENDDO
  ENDDO
  
  DO iasy = 1, nxya
    AsyType= Asy(iasy)%AsyType
    PowerDist%AsyPower2D(iasy) = PowerDist%AsyPower2D(iasy) * Asy2DNolmalizer * (hz(1)/AsyVol(iasy, 1))
  ENDDO
  
  DO iz = 1, nz
    PowerDist%Axial1DPower(iz) = PowerDist%Axial1DPower(iz) * Axial1Dnormalizer/ hz(iz)
  ENDDO
  
END SUBROUTINE 

! EDIT JSU 2019/08/14 
SUBROUTINE GenHomoKERMA(Core, FXR, Phis, GPhis, PinXS, myzb, myze, ng, ngg, lXsLib, lscat1)
  ! GENERATE HOMOGENIZED CROSS SECTION
  USE PARAM
  USE TYPEDEF,        ONLY : CoreInfo_Type,    FXRInfo_Type,      Pin_Type,        Cell_Type
  USE BasicOperation, ONLY : CP_VA,            CP_CA,             MULTI_VA
  USE Th_Mod,         ONLY : GetPinFuelTemp,   GetPinModTemp,     GetPinTemp
  USE CMFD_MOD,       ONLY : XsMac
  USE GAMCMFD_mod,    ONLY : GXSMac
  USE GamXsLib_Mod,   ONLY : GamXsBase, GetLocalQMat
  USE MacXsLib_Mod,   ONLY : MacXsBase
  USE BenchXs,        ONLY : XsBaseBen
  USE GammaTYPEDEF,   ONLY : PINHEATXS_TYPE
  USE Core_mod,       ONLY : GroupInfo
  
  IMPLICIT NONE
  
  TYPE(CoreInfo_Type), INTENT(IN) :: Core
  TYPE(FxrInfo_Type), POINTER , INTENT(IN):: FXR(:, :)
  TYPE(PINHEATXS_TYPE), POINTER, INTENT(INOUT) :: PinXs(:, :)
  REAL, POINTER, INTENT(IN) :: Phis(:, :, :), GPhis(:, :, :)
  INTEGER, INTENT(IN) :: myzb, myze, ng, ngg
  LOGICAL, INTENT(IN) :: lXsLib, lscat1
  
  ! Pointing Variables
  TYPE(FxrInfo_Type), POINTER :: myFXR
  TYPE(Pin_Type), POINTER :: Pin(:)
  !TYPE(AsyInfo_Type), POINTER :: AsyInfo
  TYPE(Cell_Type), POINTER :: Cell(:)
  REAL, POINTER :: gphiIn(:,:)
  
  INTEGER :: nCoreFxr, nCoreFsr, nxy
  INTEGER :: nCellType, nPinType, nlocalFxr, nlocalFsr, nFsrInFxr
  INTEGER :: icel, ipin, iz, ixy ,ifxr, ifsr, itype, ifsrlocal, ig
  INTEGER :: FsrIdxSt, FxrIdxSt
  INTEGER :: i, j, k, l, m, igg
  INTEGER :: iResoGrpBeg, iResoGrpEnd, norg, nChi
  
  Pin => Core%Pin
  Cell => Core%CellInfo
  
  nxy = Core%nxy
  nCoreFxr = Core%nCoreFxr; nCoreFsr = Core%nCoreFsr
  
  IF(lxslib) THEN
    norg = GroupInfo%norg; nChi = GroupInfo%nChi
    iResoGrpBeg = GroupInfo%nofg + 1
    iResoGrpEnd = GroupInfo%nofg + GroupInfo%norg
  ENDIF
  
  DO iz = myzb, myze
    DO ixy = 1, nxy
      ! Index Preparation
      FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
      icel = Pin(ixy)%Cell(iz)
      nlocalFxr = Cell(icel)%nFxr; nlocalFsr = Cell(icel)%nFsr
      ALLOCATE(gphiIn(nLocalFSR, ngg))
      DO j = 1, nLocalFxr
        ! Index Preparation
        ifxr = FxrIdxSt + j -1; nFsrInFxr = Cell(icel)%nFsrInFxr(j)
        myFxr => FXR(ifxr, iz)
        ! Macroscopic XS and matrices making to generate Photoatomic KERMA
        CALL GamXsBase(GXSMac(j), myFxr, 1, ngg, ngg, TRUE)
        IF (lXsLib) THEN
          myFxr => FXR(ifxr, iz)
          CALL MacXsBase(XsMac(j), myFxr, 1, ng, ng, 1._8, FALSE, TRUE, TRUE)  ! CrCSPFtn On        
          CALL GetLocalQMat(XSMac(j), myFxr, 1, ng, ng, TRUE)
          !Self-Sheilding Effect
          IF(myFxr%lres) THEN
             DO ig = iResoGrpBeg, iResoGrpEnd
               XsMac(j)%XsMacKf(ig) = XsMac(j)%XsMacKf(ig) * myFxr%FresokF(ig)  
             END DO
          ENDIF
          !Obtaining
        ELSE
          ifsrlocal = Cell(icel)%MapFxr2FsrIdx(1,j)
          itype=Fxr(ifxr,iz)%imix
          CALL xsbaseBen(itype, 1, ng, 1, ng, lscat1, XsMac(j))    !--r554 old >> r562c 17/02/05
        END IF
      ENDDO !Fxr sweep in the Cell
      IF (lscat1) THEN  ! Manipulate Index of gphis
        gphiIn(:, :) = GPhis(FsrIdxSt:FsrIdxSt+nlocalFsr-1, 1:ngg, iz)
      ELSE
        DO igg = 1, ngg
          gphiIn(:,igg) = GPhis(igg, FsrIdxSt:FsrIdxSt+nlocalFsr-1, iz)
        END DO
      END IF
      CALL GamCellHeatXsGen(Cell(icel), Phis(FsrIdxSt:FsrIdxSt+nlocalFsr-1, iz, 1:ng),                    &
            gphiIn, PinXS(ixy, iz),XsMac(1:nLocalFXR), GXSMac(1:nLocalFxr), ng, ngg, nLocalFxr, nLocalFsr)
      DEALLOCATE(gphiIn)
    ENDDO
  ENDDO
  
  DO iz = myzb, myze
    DO ixy = 1, nxy
      PinXS(ixy,iz)%PinTemp = GetPinTemp(Core, Fxr, iz, ixy)
      IF(Core%lFuelPlane(iz) .AND. Pin(ixy)%lFuel) THEN
        PinXS(ixy,iz)%FuelTemp = GetPinFuelTemp(Core, Fxr, iz, ixy)
        PinXS(ixy,iz)%ModTemp = GetPinModTemp(Core, Fxr, iz, ixy)
      ENDIF
    ENDDO
  ENDDO
  
  
  END SUBROUTINE
  SUBROUTINE GamCellHeatXsGen(CellInfo, Phis, GPhis, PinXs, NXSMac, GXSMac, ng, ngg, nFxr, nFsr)
    ! USING MACROSCOPIC XS IN FXR'S, GENERATE CELL HOMOGENIZED CROSS SECTION
    USE PARAM
    USE BasicOperation,   ONLY : CP_CA, CP_VA, MULTI_CA
    USE CNTL,             ONLY : nTracerCntl
    USE TYPEDEF,          ONLY : Cell_Type, XsMac_Type
    USE GammaTYPEDEF,     ONLY : GamMacXS_TYPE, PINHEATXS_TYPE
    IMPLICIT NONE
    ! INPUT VARIABLES
    TYPE(Cell_Type) :: CellInfo
    REAL :: phis(:, :), GPhis(:, :)
    TYPE(PINHEATXS_TYPE) :: PinXs
    TYPE(XsMac_Type) :: NXSMac(nFxr)
    TYPE(GamMacXS_TYPE) :: GXSMac(nFxr)
    INTEGER :: ng, ngg
    INTEGER :: nFXR, nFsr
    
    ! LOCAL VARIABLES
    REAL :: phisum, vol, volsum, localphi, Rphisum, RR(4), scatsum !BYSedit
    REAL :: gphisum, localgphi, RGphisum   ! GAMMA VARIABLES
    REAL :: onethree, threeseven
    INTEGER :: nFsrInFxr
    INTEGER :: i, j, k, ig, ig2, ireg, igg, igg2
    LOGICAL :: lfuel
    ONETHREE = one/3._8
    THREESEVEN = 3._8/7._8
    
    ! for photo-atomic reaction data (Homogenization Using Gamma Flux)
    ! PHOTON REACTION XS'S, WEIGHTING FUNCTION : PHOTON FLUX  "gphis"
    DO igg = 1, ngg  ! PHOTON ENERGY GROUP
      volsum = 0.
      gphisum = 0.
      RR = 0.   ! RR : 1-TOTAL, 2-TRANSPORT, 3-KERMA, 3+1~3+ngg-SCAT, ngg+4-ABSORPTION, ngg+5-inverse
      DO i = 1, nFxr
        nFsrInFxr = CellInfo%nFsrInFxr(i)
        DO j = 1, nFsrInFxr
          ireg = CellInfo%MapFxr2FsrIdx(j, i)
          vol = CellInfo%vol(ireg)
          volsum = volsum + vol
          localgphi = gphis(ireg, igg) * vol
          gphisum = gphisum + localgphi         ! volume-weighted summation of multigroup photon flux
          ! Reaction rate and KERMA photo-atomic reactions
          RR(1) = RR(1) + localgphi * GXSMac(i)%MacKERMA(igg)    ! RR(1) : Macro KERMA
        ENDDO ! FSR loop
      ENDDO ! FXR loop
      PinXS%Gphi(igg) = gphisum / volsum
      RGphisum = one / gphisum
      RR = RR * RGphisum
      PinXS%KERMA_P(igg) = RR(1)
    ENDDO  ! PHOTON GROUP LOOP (igg) end
    
    ! for photon generation reaction data (Homogenization Using Neutron Flux)
    ! PHOTON PRODUCTION MATRIX, WEIGHTING FUNCTION : NEUTRON FLUX "phis"
    DO ig = 1, ng
      volsum = 0.
      phisum = 0.
      RR = 0.  ! Production,  RR(igg) : photon generated within E group igg
      DO i = 1, nFxr
        nFsrInFxr = CellInfo%nFsrInFxr(i)
        DO j = 1, nFsrInFxr
          ireg = CellInfo%MapFxr2FsrIdx(j, i)
          vol = CellInfo%vol(ireg)
          volsum = volsum + vol
          localphi = phis(ireg, ig) * vol
          phisum = phisum + localphi
          ! Production from (ig) to (igg), contained at RR(igg)
          RR(1) = RR(1) + localphi * NXSMac(i)%MacKERMA_t(ig)    ! RR(1) : Macro KERMA
          RR(2) = RR(2) + localphi * NXSMac(i)%MacKERMA_p(ig)    ! RR(1) : Macro KERMA
          RR(3) = RR(3) + localphi * NXSMac(i)%MacKERMA_f(ig)    ! RR(1) : Macro KERMA
          RR(4) = RR(4) + localphi * NXSMac(i)%Xsmackf(ig)       ! RR(1) : Macro KERMA
          !RR(4) = RR(4) + localphi * NXSMac(i)%Xsmackf(igg)     ! RR(1) : Macro KERMA
        END DO ! FSR loop
      END DO ! FXR loop
      PinXS%Nphi(ig) = phisum / volsum
      Rphisum = one / phisum
      RR = RR * Rphisum
      PinXS%KERMA_T(ig)  = RR(1)
      PinXS%PhotonGen(ig) = RR(2)
      PinXS%KERMA_F(ig)  = RR(3)
      PinXS%kFis(ig)     = RR(4)
    END DO ! Neutron Energy Group loop
    
    END SUBROUTINE



  SUBROUTINE UpdateNGPhiC(PinXS, NPhiC, GPhiC, myzb, myze, nxy, ng, ngg)
    USE PARAM
    USE GammaTYPEDEF,   ONLY : PINHEATXS_TYPE
    IMPLICIT NONE
    TYPE(PINHEATXS_TYPE), POINTER  :: PinXS(:,:)
    REAL, POINTER :: NPhiC(:,:,:), GPhiC(:,:,:)
    INTEGER :: myzb, myze, nxy, ng, ngg
    
    INTEGER :: iz, ig, ixy
    
    DO iz = myzb, myze
      DO ixy = 1, nxy
        NPhiC(ixy, iz, 1:ng)  = PinXS(ixy, iz)%NPhi(1:ng)
        GPhiC(ixy, iz, 1:ngg) = PinXS(ixy, iz)%GPhi(1:ngg)
      end do
    END DO
    
    END SUBROUTINE 


END SUBROUTINE