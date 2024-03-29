#include <defines.h>
MODULE HexEdit

IMPLICIT NONE

CONTAINS
! ------------------------------------------------------------------------------------------------------------
!                                     01. HEX : Output Edit
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexOutputEdit()

USE PARAM
USE TYPEDEF,     ONLY : CoreInfo_Type, PowerDist_Type, FmInfo_Type, FXRinfo_Type, PinXs_Type, DancoffDist_Type
USE GEOM,        ONLY : Core, ng, ngg
USE Core_mod,    ONLY : FmInfo, CmInfo, eigv, ThInfo, GroupInfo
USE TH_mod,      ONLY : THVar
USE CNTL,        ONLY : nTracerCntl
USE FILES,       ONLY : IO8, io15
USE PE_mod,      ONLY : PE
USE ioutil,      ONLY : message, PrintReal1DarrayTo2Darray
USE VTK_Mod,     ONLY : ProcessVTK
USE BinOutp_Mod, ONLY : BinaryOutput
USE GamCMFD_mod, ONLY : GammaCMFD
USE CNTL,        ONLY : nTracerCntl

USE HexTst,   ONLY : HexChkEffXs
USE HexData,  ONLY : hLgc
USE HexThOut, ONLY : HexPrintFuelAvgTemp, HexPrintFuelCentTemp, HexPrintFuelSurfTemp, HexPrintCoolantTemp

#ifdef MPI_ENV
USE MpiComm_mod, ONLY : BCAST, MPI_SYNC
#endif

IMPLICIT NONE

TYPE(PowerDist_Type)   :: PowerDist
TYPE(PowerDist_Type)   :: PhotonPower, NeutronPower, TotalExpPower
TYPE(DancoffDist_Type) :: DancoffDist

TYPE(Fxrinfo_type), POINTER, DIMENSION(:,:) :: Fxr
TYPE(PinXS_Type),   POINTER, DIMENSION(:,:) :: PinXS

REAL, POINTER, DIMENSION(:,:,:) :: PhiS, PhiC

INTEGER :: io
LOGICAL :: Master
! ----------------------------------------------------

Master = TRUE

#ifdef MPI_ENV
Master = PE%Cmfdmaster
#endif
! ----------------------------------------------------
!                1. PRINT : Banner
! ----------------------------------------------------
IF(Master) THEN
  WRITE(mesg, *) 'Writing Calculation Results in Output File...'
  CALL message(io8, TRUE, TRUE, mesg)
  
  WRITE(mesg, '(14X, A, I5, 4X, A, I5)') 'Case No :', nTracerCntl%CaseNo, 'Calculation No :', nTracerCntl%CalNoId 
  CALL message(io8, FALSE, TRUE, mesg)
ENDIF

IF(Master) THEN
  io = io8
  
  WRITE(io, *)
  WRITE(io, '(A)') '========================================================================'
  WRITE(io, '(5X, A)') 'Results'
  WRITE(io, '(5X, A, I5)') 'Case No        :', nTracerCntl%CaseNo
  WRITE(io, '(5X, A, I5)') 'Calculation No :', nTracerCntl%CalNoId
  WRITE(io, '(A)') '========================================================================'

  WRITE(io, *)
  WRITE(io, '(7X, A7, 3X, F8.5)') 'k-eff =', eigv
  WRITE(io, *)
END IF
! ----------------------------------------------------
!                2. SET : Power
! ----------------------------------------------------
IF (hLgc%lSngCel) RETURN

CALL HexCorePowerCal(Core, CmInfo, PowerDIst, ng, nTracerCntl, PE)
IF (nTracerCntl%lED .AND. nTracerCntl%lrestrmt) CALL TreatDancoffData(Core,DancoffDist, PE)
#ifdef MPI_ENV
  CALL BCAST(PowerDist%Fm3DNormalizer,  PE%MPI_CMFD_COMM)
  CALL BCAST(PowerDist%Pin3DNormalizer, PE%MPI_CMFD_COMM)
#endif
IF(nTracerCntl%lExplicitKappa) CALL HexCalExplicitPower(Core, GammaCMFD, FmInfo, PhotonPower, NeutronPower, TotalExpPower, ng, ngg, nTracerCntl, PE)
! ----------------------------------------------------
!                3. PRINT : Power
! ----------------------------------------------------
IF(Master) THEN
  WRITE(io, '(A)') '-------------------------Power Distribution------------------------------'
  SELECT CASE(nTracerCntl%pout_mode)
  CASE(0)
  !CALL PrintRadialPower(io, Core, PowerDist)
    WRITE(io, '(A)') '--- Fission Power' 
  CALL HexPrint2DPinPower   (io, Core, PowerDist)
  CALL HexPrintLocalPinPower(io, Core, PowerDist)
  
  !CALL Asy3DPower(io, Core, PowerDist)
  
  CALL HexAxialAvgAsy2DPower(io, Core, PowerDist)
  CALL Axial1DPower         (io, Core, PowerDist)
  
  !IF (nTracerCntl%lED) CALL PrintDancoffFactor(io, Core, DancoffDist)
  CASE(1)
  CASE(2)
    WRITE(io, '(A, ES14.6)') '--- Fission Power', PowerDist%pwsum
  CALL HexPrint2DPinPower   (io, Core, PowerDist)
  CALL HexPrintLocalPinPower(io, Core, PowerDist)
  
  !CALL Asy3DPower(io, Core, PowerDist)
  
  CALL HexAxialAvgAsy2DPower(io, Core, PowerDist)
  CALL Axial1DPower         (io, Core, PowerDist)

    WRITE(io, '(A, ES14.6)') '--- Neutron Power', NeutronPower%pwsum
  CALL HexPrint2DPinPower   (io, Core, NeutronPower)
  CALL HexPrintLocalPinPower(io, Core, NeutronPower)
  
  !CALL Asy3DPower(io, Core, PowerDist)
  
  CALL HexAxialAvgAsy2DPower(io, Core, NeutronPower)
  CALL Axial1DPower         (io, Core, NeutronPower)

    WRITE(io, '(A, ES14.6)') '--- Photon Power' , PhotonPower%pwsum
  CALL HexPrint2DPinPower   (io, Core, PhotonPower)
  CALL HexPrintLocalPinPower(io, Core, PhotonPower)
  
  !CALL Asy3DPower(io, Core, PowerDist)
  
  CALL HexAxialAvgAsy2DPower(io, Core, PhotonPower)
  CALL Axial1DPower         (io, Core, PhotonPower)

    WRITE(io, '(A, ES14.6)') '--- Explicit Power', TotalExpPower%pwsum
  CALL HexPrint2DPinPower   (io, Core, TotalExpPower)
  CALL HexPrintLocalPinPower(io, Core, TotalExpPower)
  
  !CALL Asy3DPower(io, Core, PowerDist)
  
  CALL HexAxialAvgAsy2DPower(io, Core, TotalExpPower)
  CALL Axial1DPower         (io, Core, TotalExpPower)
  END SELECT
ENDIF
! ----------------------------------------------------
!                4. PRINT : Flux
! ----------------------------------------------------
CALL HexAxAvg2DFlux    (io, core, CmInfo, PowerDist, ng, PE)
!CALL Asy3DFlux(io, core, CmInfo, PowerDist, ng, PE)
CALL HexRadAvg1DFlux   (io, core, CmInfo, PowerDist, ng, PE)
CALL HexCoreSpectrum   (io, core, CmInfo, PowerDist, ng, PE)
CALL HexPin3DFlux      (io, core, CmInfo, PowerDist, ng, PE, nTracerCntl)
!CALL PrintDepletion(io, Core, FmInfo, CmInfo, nTracerCntl, PE)
! ----------------------------------------------------
!                5. PRINT : Temp
! ----------------------------------------------------
IF(nTracerCntl%lfeedback .AND. MASTER) THEN
  CALL HexPrintFuelAvgTemp (io, Core, THInfo, THVar)
  CALL HexPrintFuelCentTemp(io, Core, THInfo)
  CALL HexPrintFuelSurfTemp(io, Core, THInfo, THVar)
  CALL HexPrintCoolantTemp (io, Core, THInfo)
END IF  
! ----------------------------------------------------
!                6. PRINT : Etc, Normally FALSE
! ----------------------------------------------------
IF(Master) WRITE(io, '(a)') '========================================================================'

IF(nTRACERCntl%OutpCntl%lBoutp)           CALL BinaryOutput()
IF(nTRACERCntl%OutpCntl%VisMod .GT. 0)    CALL ProcessVTK      (Core, FMInfo, CmInfo, PowerDist, ThInfo, GroupInfo, nTracerCntl, PE)
IF(nTracerCntl%OutpCntl%lPinXsOut)        CALL PrintPinXs(io15, Core,         CmInfo,                    GroupInfo, nTracerCntl, PE)
IF(nTracerCntl%OutpCntl%lCspGenOut)       CALL CspGenOut   (io, Core, FmInfo, CmInfo,                    GroupInfo, nTracerCntl, PE)
IF(nTracerCntl%OutpCntl%nRegXsOut .GT. 0) CALL HexChkEffXs
IF(nTracerCntl%OutpCntl%nRegXsOut .GT. 0) CALL HexProcessEffXs (Core, FmInfo,                    ThInfo, GroupInfo, nTracerCntl, PE)

IF(nTracerCntl%OutpCntl%nRegMGXsOut .GT. 0) CALL ProcessEffMGXs(Core, FmInfo, CmInfo, ThInfo, GroupInfo, nTracerCntl, PE)
IF(nTracerCntl%OutpCntl%nRegMAT0Out .GT. 0) CALL ProcessEffMAT0(Core, FmInfo, GroupInfo, nTracerCntl, PE)
IF(nTracerCntl%OutpCntl%nRegPhimOut .GT. 0) CALL ProcessPhims(Core, FmInfo, GroupInfo, nTracerCntl, PE)
IF(nTracerCntl%OutpCntl%lSSPHout)           CALL ProcessSSPHfactor(Core,nTracerCntl)

CALL HexFreePowerDist(PowerDist, PE)
! ----------------------------------------------------

END SUBROUTINE HexOutputEdit
! ------------------------------------------------------------------------------------------------------------
!                                     02. PRINT : 2D Pin Power
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexPrint2DPinPower(io, Core, PowerDist)

USE PARAM,   ONLY : ZERO, BLANK
USE TYPEDEF, ONLY : CoreInfo_Type, PowerDist_Type, Pin_Type
USE CNTL,    ONLY : nTracerCntl
USE geom,    ONLY : nAsyType0

USE HexData, ONLY : hAsy, hAsyTypInfo, hPinInfo, Asy1Dto2Dmap, nHexPin

IMPLICIT NONE

INTEGER :: io
TYPE(CoreInfo_Type) :: Core
TYPE(PowerDist_Type) :: PowerDist
TYPE(Pin_Type), POINTER, DIMENSION(:) :: Pin
! ----------------------------------------------------
INTEGER :: nxya, nz, nCx, iAsyTyp, iGeoTyp
INTEGER :: iasy, ixy, ix, iy, iPin, jPin, ivTyp

REAL :: EffPinNum
REAL :: VtxTypVolInv(7) = [6., 2., 1., 1., 1., 2., 2.]

CHARACTER(12) :: cTmp

CHARACTER(12), POINTER, DIMENSION(:) :: sTmp

301 FORMAT('Assembly', I5, X, 'at', X, '(', I7, I7, ')')
302 FORMAT(2X, 100A12)
! ----------------------------------------------------

nxya = Core%nxya
nz   = Core%nz
Pin => Core%Pin

EffPinNum = ZERO

DO iPin = 1, nHexPin
  IF (Pin(iPin)%lfuel) THEN
    EffPinNum = EffPinNum + hPinInfo(iPin)%Wt
  END IF
END DO

nCx = 0

DO iAsyTyp = 1, nAsyType0
  nCx = max(nCx, hAsyTypInfo(iAsyTyp)%nPin)
END DO

ALLOCATE (sTmp (2*nCx-1))

WRITE(io, '(A)') '- 2D Pin Power -'
WRITE(io, '(A)') '- Coordinate : 1st = from NE to SW, 2nd = from NNW to SSE'
WRITE(io, *)
WRITE(io, '(A, F12.5)') 'Eff Pin Num', EffPinNum
WRITE(io, *)
! ----------------------------------------------------
DO iAsy = 1, nxya
  WRITE(io, 301), iAsy, Asy1Dto2Dmap(1, iAsy), Asy1Dto2Dmap(2, iAsy)

  iAsyTyp = hAsy(iAsy)%AsyTyp
  iGeoTyp = hAsy(iAsy)%GeoTyp
  nCx     = hAsyTypInfo(iAsyTyp)%nPin
  
  DO iy = 1, 2*nCx - 1   ! NE  to SW
    sTmp = BLANK
    
    DO ix = 1, 2*nCx - 1 ! NNW to SSE
      iPin = hAsyTypInfo(iAsyTyp)%Pin2Dto1Dmap(ix, iy)
      
      IF (iPin < 1) THEN
        cTmp = BLANK
      ELSE
        jPin  = hAsyTypInfo(iAsyTyp)%PinLocIdx(iGeoTyp, iPin)
        ivTyp = hAsyTypInfo(iAsyTyp)%PinVtxTyp(iGeoTyp, iPin)
        
        IF (jPin < 1) THEN
          WRITE (cTmp, '(F12.5)') 0._8
        ELSE
          ixy = hAsy(iAsy)%PinIdxSt + jPin - 1 ! Global Pin Idx
          
          WRITE (cTmp, '(F12.5)') PowerDist%PinPower2D(iPin, iAsy) * hPinInfo(ixy)%Vol(1) * VtxTypVolInv(ivTyp)
        END IF
      END IF
      
      WRITE (sTmp(ix), '(A12)') cTmp
    END DO
    
    WRITE(io, 302), sTmp(1:2*nCx - 1)
  END DO
END DO

WRITE(io, *)

NULLIFY (Pin)

DEALLOCATE (sTmp)
! ----------------------------------------------------

END SUBROUTINE HexPrint2DPinPower
! ------------------------------------------------------------------------------------------------------------
!                                     04. PRINT : Local Pin Power
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexPrintLocalPinPower(io, Core, PowerDist)

USE PARAM,   ONLY : ZERO, BLANK
USE TYPEDEF, ONLY : CoreInfo_Type, PowerDist_Type, Pin_Type
USE CNTL,    ONLY : nTracerCntl
USE geom,    ONLY : nAsyType0

USE HexData, ONLY : hAsy, hAsyTypInfo, hPinInfo, Asy1Dto2Dmap, nHexPin

IMPLICIT NONE

INTEGER :: io
TYPE(CoreInfo_Type) :: Core
TYPE(PowerDist_Type) :: PowerDist
TYPE(Pin_Type), POINTER, DIMENSION(:) :: Pin
! ----------------------------------------------------
INTEGER :: nxya, nz, nCx, iAsyTyp, iGeoTyp
INTEGER :: iasy, iz, ixy, ix, iy, iPin, jPin, ivTyp

REAL :: EffPinNum
REAL :: VtxTypVolInv(7) = [6., 2., 1., 1., 1., 2., 2.]

CHARACTER(12)  :: cTmp

CHARACTER(12), POINTER, DIMENSION(:) :: sTmp

301 FORMAT('Assembly', i5, x,'at' ,x,  '(', i7, i7,' )')
302 FORMAT(2X, 100A12)
! ----------------------------------------------------

nxya = Core%nxya
nz   = Core%nz
Pin => Core%Pin

EffPinNum = ZERO

DO iPin = 1, nHexPin
  IF (Pin(iPin)%lfuel) THEN
    EffPinNum = EffPinNum + hPinInfo(iPin)%Wt
  END IF
END DO

nCx = 0

DO iAsyTyp = 1, nAsyType0
  nCx = max(nCx, hAsyTypInfo(iAsyTyp)%nPin)
END DO

ALLOCATE (sTmp (2*nCx-1))

WRITE(io, '(A)') '- Local Pin Power -'
WRITE(io, '(A)') '- Coordinate : 1st = from NE to SW, 2nd = from NNW to SSE'
WRITE(io, *)
WRITE(io, '(A, F12.5)') 'Eff Pin Num', EffPinNum
WRITE(io, *)
! ----------------------------------------------------
DO iz = 1, nz
  IF(.NOT. nTracerCntl%lExplicitKappa .AND. .NOT. Core%lFuelPlane(iz)) CYCLE
  WRITE(io, '(X, A5, I5)'), 'Plane', iz
  
  DO iAsy = 1, nxya
    WRITE(io, 301), iAsy, Asy1Dto2Dmap(1, iAsy), Asy1Dto2Dmap(2, iAsy)
    
    iAsyTyp = hAsy(iAsy)%AsyTyp
    iGeoTyp = hAsy(iAsy)%GeoTyp
    nCx     = hAsyTypInfo(iAsyTyp)%nPin
    
    DO iy = 1, 2*nCx - 1   ! NE  to SW
      sTmp = BLANK
      
      DO ix = 1, 2*nCx - 1 ! NNW to SSE
        iPin = hAsyTypInfo(iAsyTyp)%Pin2Dto1Dmap(ix, iy)
        
        IF (iPin < 1) THEN
          cTmp = BLANK
        ELSE
          jPin  = hAsyTypInfo(iAsyTyp)%PinLocIdx(iGeoTyp, iPin)
          ivTyp = hAsyTypInfo(iAsyTyp)%PinVtxTyp(iGeoTyp, iPin)
          
          IF (jPin < 1) THEN
            WRITE (cTmp, '(F12.5)') 0._8
          ELSE
            ixy = hAsy(iAsy)%PinIdxSt + jPin - 1 ! Global Pin Idx
            
            WRITE (cTmp, '(F12.5)') PowerDist%PinPower3D(iPin, iAsy, iz) * hPinInfo(ixy)%Vol(iz) * VtxTypVolInv(ivTyp)
          END IF
        END IF
        
        WRITE (sTmp(ix), '(A12)') cTmp
      END DO
      
      WRITE(io, 302), sTmp(1:2*nCx - 1)
    END DO
  END DO
ENDDO

WRITE(io, *)

NULLIFY (Pin)

DEALLOCATE (sTmp)
! ----------------------------------------------------

END SUBROUTINE HexPrintLocalPinPower
! ------------------------------------------------------------------------------------------------------------
!                                     06. Axial Avg Asy 2D Power
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexAxialAvgAsy2DPower(io, Core, PowerDist)

USE allocs
USE PARAM,   ONLY : ZERO, BLANK
USE TYPEDEF, ONLY : CoreInfo_Type, PowerDist_Type, Asy_Type
USE CNTL,    ONLY : nTracerCntl
USE HexData, ONLY : hLgc, nAsyCore, Asy2Dto1Dmap

IMPLICIT NONE

INTEGER :: io
TYPE(CoreInfo_Type) :: Core
TYPE(PowerDist_Type) :: PowerDist
TYPE(Asy_Type), POINTER, DIMENSION(:) :: Asy
! ----------------------------------------------------
REAL :: maxv

INTEGER :: nxy, nLgh, i, j, ix, iy, iAsy

CHARACTER(12)  :: cTmp

CHARACTER(12), POINTER, DIMENSION(:) :: sTmp

302 FORMAT(2X, 100A12)
! ----------------------------------------------------

nxy = Core%nxya

WRITE(io, '(A)') '- Axially Averaged 2D Assembly Power -'
WRITE(io, *)

IF (hLgc%l360) THEN
  nLgh = 2*nAsyCore - 1
ELSE
  nLgh = nAsyCore
END IF

ALLOCATE (sTmp (nLgh))
! ----------------------------------------------------
DO iy = 1, nLgh
  sTmp = BLANK
  
  DO ix = 1, nLgh
    iAsy = Asy2Dto1Dmap(ix, iy)
    
    IF (iAsy < 1) THEN
      cTmp = BLANK
    ELSE
      WRITE (cTmp, '(F12.5)') PowerDist%AsyPower2D(iAsy)
    END IF
    
    WRITE (sTmp(ix), '(A12)') cTmp
  END DO
  
  WRITE(io, 302), sTmp(1:nLgh)
END DO
! ----------------------------------------------------
maxv = maxval(PowerDist%AsyPower2D(1:nxy))

WRITE(io, *)
WRITE(io, '(5X, A7, 2X, F7.3)') 'Max. = ', maxv
WRITE(io, *)

DEALLOCATE (sTmp)
! ----------------------------------------------------

END SUBROUTINE HexAxialAvgAsy2DPower
! ------------------------------------------------------------------------------------------------------------
!                                     07. Ax Avg 2D Flux
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexAxAvg2DFlux(io, core, CmInfo, PowerDist, ng, PE)

USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type, PowerDist_Type, CMInfo_Type, PE_TYPE, AsyInfo_Type, Asy_Type, PinXs_Type
USE CNTL,    ONLY : nTracerCntl
USE HexData, ONLY : hPinInfo

#ifdef MPI_ENV
USE BasicOperation, ONLY : CP_CA, CP_VA
USE MpiComm_mod,    ONLY : Reduce
#endif

IMPLICIT NONE

INTEGER :: io, ng
TYPE(CoreInfo_Type)  :: Core
TYPE(CMInfo_Type)    :: CMInfo
TYPE(PowerDist_Type) :: PowerDist
TYPE(PE_TYPE)        :: PE
! ----------------------------------------------------
TYPE(AsyInfo_Type), POINTER, DIMENSION(:)   :: AsyInfo
TYPE(Asy_Type),     POINTER, DIMENSION(:)   :: Asy
TYPE(PinXS_Type),   POINTER, DIMENSION(:,:) :: PinXs

REAL, POINTER, DIMENSION(:)     :: hz, row
REAL, POINTER, DIMENSION(:,:)   :: PinVol, AsyVol, Avg2DFlx
REAL, POINTER, DIMENSION(:,:,:) :: PhiC

INTEGER :: nz, nxy, nxya, myzb, myze, nAsyType, AsyType
INTEGER :: iz, ig, iasy, ixy, i, j, k
LOGICAL :: master

#ifdef MPI_ENV
REAL, POINTER, DIMENSION(:) :: Buf
#endif
! ----------------------------------------------------

Master = PE%Cmfdmaster

PhiC    => CmInfo%PhiC
PinXs   => CmInfo%PinXs
AsyInfo => Core%AsyInfo
Asy     => Core%Asy
hz      => Core%hz 
PinVol  => Core%PinVol
AsyVol  => Core%AsyVol
nz       = Core%nz
myzb     = PE%myzb
myze     = PE%myze
nAsyType = Core%nAsyType
nxya     = Core%nxya
nxy      = Core%nxy

ALLOCATE(Avg2DFlx(nxya, ng))

#ifdef MPI_ENV
ALLOCATE(Buf(nxya)) 
#endif
! ----------------------------------------------------
DO ig = 1, ng
  Avg2DFlx(:, ig) = 0

  DO iasy = 1, nxya
    AsyType = Asy(iasy)%AsyType
    nxy     = AsyInfo(AsyType)%nxy

    Do iz = myzb, myze  ! Axial Sweep
      DO i = 1, nxy     ! Local Cell Sweep
        ixy = Asy(iasy)%GlobalPinIdx(i)
        
        IF (ixy < 0) CYCLE

        Avg2DFlx(iAsy, ig) = Avg2DFlx(iAsy, ig) + hPinInfo(ixy)%Vol(iz) * Phic(ixy, iz, ig)
      END DO
    END DO
  END DO

#ifdef MPI_ENV
  CALL CP_CA(Buf, zero, nxya)
  CALL REDUCE(Avg2DFlx(:, ig), Buf(:), nxya, PE%MPI_CMFD_COMM, .FALSE.)
  IF(Master) CALL CP_VA(Avg2DFlx(:,ig), Buf, nxya)
#endif  
ENDDO
! ----------------------------------------------------
#ifdef MPI_ENV
DEALLOCATE(Buf)
#endif

IF(Master) THEN
  WRITE(io, '(A)') ' -  Axially Averaged 2D Flux -'
  WRITE(io, '(A)') ' -  Groups per Assembly'
  WRITE(io, '(A)') '- Order of assembly : Start from left-top, proceed first x-dir, later y-dir -'
  WRITE(io, *)

  DO iasy = 1, nxya
    WRITE(io, '(8X, A5, 3X, A8, 3X, I5)') 'Group', 'Assembly', iAsy
    
    DO ig = 1, ng
      WRITE(io, '(8X, I5, 3X, 200(ES15.4))') ig, Avg2DFlx(iAsy, ig)
    END DO
    
    WRITE(io, *)
  END DO
END IF

NULLIFY (PhiC)
NULLIFY (PinXs)
NULLIFY (AsyInfo)
NULLIFY (Asy)
NULLIFY (hz)
NULLIFY (PinVol)
NULLIFY (AsyVol)

DEALLOCATE(Avg2DFlx)
! ----------------------------------------------------

END SUBROUTINE HexAxAvg2DFlux
! ------------------------------------------------------------------------------------------------------------
!                                     08. Rad Avg 1D Flux
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexRadAvg1DFlux(io, core, CmInfo, PowerDist, ng, PE)

USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type, PowerDist_Type, CMInfo_Type, PE_TYPE, AsyInfo_Type, Asy_Type, PinXs_Type

#ifdef MPI_ENV
USE MpiComm_mod, ONLY : SENDRECV
#endif

IMPLICIT NONE

INTEGER :: io, ng
TYPE(CoreInfo_Type)  :: Core
TYPE(CMInfo_Type)    :: CMInfo
TYPE(PowerDist_Type) :: PowerDist
TYPE(PE_TYPE)        :: PE
! ----------------------------------------------------
TYPE(AsyInfo_Type), POINTER, DIMENSION(:)   :: AsyInfo
TYPE(Asy_Type),     POINTER, DIMENSION(:)   :: Asy
TYPE(PinXS_Type),   POINTER, DIMENSION(:,:) :: PinXs

REAL, POINTER, DIMENSION(:)     :: hz
REAL, POINTER, DIMENSION(:,:)   :: PinVol, AsyVol, Avg1DFlx
REAL, POINTER, DIMENSION(:,:,:) :: PhiC

INTEGER :: nz, nxy, nxya, myzb, myze, nAsyType, AsyType
INTEGER :: iz, ig, iasy, ixy, i, j, k
LOGICAL :: master

#ifdef MPI_ENV
INTEGER :: nproc, comm, iz1, iz2
#endif
! ----------------------------------------------------

master = PE%CmfdMaster

PhiC    => CmInfo%PhiC
PinXs   => CmInfo%PinXs
AsyInfo => Core%AsyInfo
Asy     => Core%Asy
hz      => Core%hz; 
PinVol  => Core%PinVol
AsyVol  => Core%AsyVol
nz       = Core%nz
myzb     = PE%myzb
myze     = PE%myze
nAsyType = Core%nAsyType
nxya     = Core%nxya; 
nxy      = Core%nxy

ALLOCATE(Avg1DFlx(ng, nz))
! ----------------------------------------------------
Avg1DFlx = 0

DO ig = 1, ng  
  DO iasy = 1, nxya
    AsyType = Asy(iasy)%AsyType
    nxy = AsyInfo(AsyType)%nxy

    Do iz = myzb, myze   !Axial Sweep
      DO i = 1, nxy  !Cell Sweep
        ixy = Asy(iasy)%GlobalPinIdx(i) ! Global Pin Idx
        IF (ixy < 0) CYCLE

        Avg1DFlx(ig, iz) = Avg1DFlx(ig, iz) + Pinvol(ixy, iz) * Phic(ixy, iz, ig)
      ENDDO
    ENDDO
  ENDDO !
ENDDO

#ifdef MPI_ENV
  nproc = PE%nCmfdProc
  Comm  = PE%MPI_CMFD_COMM

  IF(Master) THEN
    DO i = 1, nproc - 1
      iz1 = PE%AxDomRange(1, i)
      iz2 = PE%AxDomRange(2, i)    

      CALL SendRecv(Avg1DFlx(:, iz1:iz2), ng, iz2 - iz1 + 1, i, FALSE, COMM)
    ENDDO
  ELSE
    CALL SendRecv(Avg1DFlx(:, myzb:myze), ng, myze - myzb + 1, 0, TRUE, COMM)
  ENDIF
#endif
! ----------------------------------------------------
IF(master) THEN
  WRITE(io, '(A)') ' -  Radially Averaged 1D Flux -'
  WRITE(io, *)
  WRITE(io, '(A)') '                --- > Plane #'
  WRITE(io, '(2X, A, 200(I10, 5X))') 'Energy Group',(i, i=1,nz)

  DO ig = 1, ng
    WRITE(io, '(I10, 4X, 200(ES15.4))') ig, (Avg1DFlx(ig, iz), iz = 1, nz)
  ENDDO
  WRITE(io, *)
ENDIF

NULLIFY (PhiC, PinXs, AsyInfo, Asy, hz, PinVol, AsyVol)
! ----------------------------------------------------

END SUBROUTINE HexRadAvg1DFlux
! ------------------------------------------------------------------------------------------------------------
!                                     09. Core Spectrum
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexCoreSpectrum(io, core, CmInfo, PowerDist, ng, PE)

USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type,   PowerDist_Type,  CMInfo_Type,   &
                    PE_TYPE
USE HexData, ONLY : hPinInfo

#ifdef MPI_ENV
USE MpiComm_mod, ONLY : REDUCE
#endif

IMPLICIT NONE

INTEGER :: io, ng
TYPE(CoreInfo_Type) :: Core
TYPE(CMInfo_Type) :: CMInfo
TYPE(PowerDist_Type) :: PowerDist
TYPE(PE_TYPE) :: PE

REAL, POINTER :: PhiC(:, :, :)
REAL, POINTER :: PinVol(:, :)
REAL, POINTER :: hz(:)
INTEGER :: nz, nxy, nxya, myzb, myze
INTEGER :: iz, ig, iasy, ixy, i, j, k
LOGICAL :: master
REAL :: Spectrum(ng)
#ifdef MPI_ENV
REAL :: buf(ng)
#endif
! ----------------------------------------------------

Master  = PE%Master
myzb    = PE%myzb
myze    = PE%myze
PhiC   => CmInfo%PhiC;
hz     => Core%hz
PinVol => Core%PinVol
nxy     = Core%nxy
nz      = Core%nz
! ----------------------------------------------------
DO ig = 1, ng
  Spectrum(ig) = 0

  DO iz = myzb, myze
    DO ixy = 1, nxy
      Spectrum(ig) = Spectrum(ig) + Phic(ixy, iz, ig) * hPinInfo(ixy)%Vol(iz)
    ENDDO
  ENDDO
ENDDO
! ----------------------------------------------------
#ifdef MPI_ENV
CALL REDUCE(Spectrum, buf, ng, PE%MPI_CMFD_COMM, FALSE)
Spectrum=buf
#endif

IF(Master) THEN
  WRITE(io, '(A)') ' - Core Spectrum -'
  WRITE(io, *)
  DO ig = 1, ng
    WRITE(io, '(5X, I5, ES15.4)') ig, Spectrum(ig)
  ENDDO
  WRITE(io, *)
ENDIF

NULLIFY (PhiC, hz, PinVol)
! ----------------------------------------------------

END SUBROUTINE HexCoreSpectrum
! ------------------------------------------------------------------------------------------------------------
!                                     10. Pin 3D Flux
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexPin3DFlux(io, core, CmInfo, PowerDist, ng, PE, nTracerCntl)
USE PARAM
USE TYPEDEF,             ONLY : CoreInfo_Type,   PowerDist_Type,    CmInfo_Type,  &
                                PE_Type,                                          &
                                Asy_Type,        AsyInfo_Type
USE CNTL,                ONLY : nTracerCntl_Type
USE BasicOperation,      ONLY : CP_CA,           CP_VA

#ifdef MPI_ENV
USE MPIComm_Mod,         ONLY : REDUCE
#endif

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(PowerDist_Type) :: PowerDist
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
INTEGER :: io, ng

TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
REAL, POINTER :: PhiC(:, :, :)

INTEGER, PARAMETER :: noutline = 10

INTEGER :: i, j, k, ibeg, iend
INTEGER :: ixy, ixya, ix, iy, iz, ig

INTEGER :: nz, myzb, myze
INTEGER :: npinout, nstep, nout

CHARACTER(5) :: GrpIdx
REAL, POINTER :: OutDat(:, :, :), buf(:, :, :)
! ----------------------------------------------------

Asy     => Core%Asy
AsyInfo => Core%AsyInfo
PhiC    => CmInfo%PhiC

nz      = Core%nz
myzb    = PE%myzb
myze    = PE%myze
npinout = nTracerCntl%OutpCntl%FluxOutList(1, 0)
! ----------------------------------------------------
ALLOCATE(OutDat (ng, nz, noutline))
ALLOCATE(Buf    (ng, nz, noutline))

nstep = INT(npinout / noutline)
IF(mod(npinout, noutline) .NE. 0) nstep = nstep + 1

IF(PE%MASTER) THEN
  WRITE(io, '(A)') ' -  Flux Distribution for selected Pins -  '
  WRITE(io, '(7X, A)') 'Pin Identifier'
  DO i = 1, npinout
    WRITE(io, '(7X, I5, 3X, A2, 3X, (A1, 4I4, A3))') i, '=>', '[', nTracerCntl%OutpCntl%FluxOutList(1:4, i), ']'
  ENDDO
  WRITE(io, *)
ENDIF
! ----------------------------------------------------
DO j = 1, nstep 
  ibeg = (j - 1) * noutline + 1
  iend = j * noutline
  iend = min(iend, npinout)
  nout = iend - ibeg +1
  CALL CP_CA(OutDat, 0._8, ng, nz, noutline)

  DO i = ibeg, iend
    k  = i - ibeg + 1

    ixya = nTracerCntl%OutpCntl%FluxOutList(1, i) ! Start from left-top, proceed first x-dir, later y-dir
    ixy  = nTracerCntl%OutpCntl%FluxOutList(2, i) ! Start from center of assembly, proceed first anti-clockwisely from top-most, later outer ring
    ixy  = Asy(ixya)%GlobalPinIdx(ixy)

    DO ig = 1, ng
      DO iz = myzb, myze
        OutDat(ig, iz, k) = PhiC(ixy, iz, ig)
      ENDDO
    ENDDO
  ENDDO

#ifdef MPI_ENV
  CALL REDUCE(OutDat, Buf, ng, nz, noutline, PE%MPI_CMFD_COMM, .FALSE.)
  CALL CP_VA(OutDat, buf, ng, nz, noutline)
#endif

  IF(.NOT. PE%MASTER) CYCLE
  WRITE(io, '(7X, A)') 'Pin Index ->'
  WRITE(io, '(4X, A5, 2X, A5, 2X, 100(I10, 2X))') 'GROUP', 'PLANE', (i, i = ibeg, iend)

  DO ig = 1, ng
    WRITE(GrpIdx, '(I5)') ig
    DO iz = 1, nz 
      IF(iz .NE. 1) GrpIdx=''
      WRITE(io, '(4X, A5, 2X, I5, 2X, 100(E10.3, 2X))') GrpIdx, iz, (OutDat(ig, iz, k), k = 1, nout)
    ENDDO
  ENDDO
  WRITE(io, *)
ENDDO

NULLIFY(Asy, AsyInfo)
! ----------------------------------------------------

END SUBROUTINE HexPin3DFlux
! ------------------------------------------------------------------------------------------------------------
!                                     11. FREE : Power Dist
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexFreePowerDist(PowerDist, PE)

USE PARAM
USE TypeDef, ONLY : PowerDist_Type, PE_TYPE

IMPLICIT NONE

TYPE(PowerDist_Type) :: PowerDist
TYPE(PE_TYPE)        :: PE
! ----------------------------------------------------

DEALLOCATE (PowerDist%PinPower3D)

IF (PE%CMFDmaster) THEN
  DEALLOCATE(PowerDist%PinPower2D)
  DEALLOCATE(PowerDist%AsyPower2D)
  DEALLOCATE(PowerDist%AsyPower3D)
  DEALLOCATE(PowerDist%Axial1DPower)
END IF
! ----------------------------------------------------

END SUBROUTINE HexFreePowerDist
! ------------------------------------------------------------------------------------------------------------

END MODULE HexEdit