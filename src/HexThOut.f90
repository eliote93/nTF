MODULE HexThOut
  
IMPLICIT NONE

CONTAINS
! ------------------------------------------------------------------------------------------------------------
!                                     01. PRINT : TH Results
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexPrintTHResults(io, Output, nxyp, nxp, nyp, nz)

USE allocs
USE PARAM, ONLY : BLANK

IMPLICIT NONE

INTEGER :: io, nxyp, nxp, nyp, nz
REAL :: OutPut(nz, nxyp)
INTEGER :: iyp, iz, ixp, kbeg, xbeg, xend, x1, x2

CHARACTER(512) :: sTmp
! ----------------------------------------------------

WRITE(io, '(A7, A7, x, A10)') 'Y-Idx', 'Z-Idx','--> X-Idx'
WRITE(io, '(14x, 200I8)') (ixp, ixp = 1, nyp)

kbeg = 1
xbeg = 1
xend = nxp

DO iyp = 1, nyp
  DO iz = 1, nz
    sTmp = BLANK
    
    IF (iz .EQ. 1) WRITE(sTmp(1:5), '(I5)') iyp
    
    WRITE(sTmp(6:10), '(I5)') iz
    
    DO ixp = xbeg, xend
      x1 = 14 + ixp*8 - 7
      x2 = 14 + ixp*8
      
      WRITE(sTmp(x1:x2), '(F8.2)') OutPut(iz, kbeg + ixp - xbeg)
    END DO
    
    WRITE(io, '(A)') sTmp(1:x2)
  END DO
  
  WRITE(io, *)
  
  kbeg = kbeg + xend - xbeg + 1
  
  IF (iyp < nxp) THEN
    xend = xend + 1
  ELSE
    xbeg = xbeg + 1
  END IF
END DO

WRITE(io, *)
! ----------------------------------------------------

END SUBROUTINE HexPrintTHResults
! ------------------------------------------------------------------------------------------------------------
!                                     02. PRINT : Fuel Avg Temp
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexPrintFuelAvgTemp(io, Core, THInfo, THVar)

USE allocs
USE PARAM,   ONLY : ZERO
USE TYPEDEF, ONLY : CoreInfo_Type, ThInfo_Type, FuelTH_Type, THVar_Type, AsyInfo_Type, Asy_Type
USE HexData, ONLY : hAsyTypInfo

IMPLICIT NONE

INTEGER :: io

TYPE(CoreInfo_Type) :: Core
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)

TYPE(ThInfo_Type) :: ThInfo
TYPE(FuelTH_Type), POINTER :: FuelTH(:)
TYPE(THVar_Type) :: THVar

REAL, POINTER :: OutPut(:, :)

INTEGER :: iz, ixy, ixyp, iasy, AsyType
INTEGER :: nz, nxya, nxp, nyp, nxyp, navg
! ----------------------------------------------------

nz   = Core%nz
nxya = Core%nxya
navg = THVar%npr5

AsyInfo => Core%AsyInfo
Asy     => Core%Asy
FuelTH  => THInfo%FuelTH

301 FORMAT('Assembly', I5, x, 'at', x, '(', I7, I7,' )')

WRITE(io, '(A)') '- Fuel Pin Average Temperature -'
WRITE(io, *)

DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  
  nxyp = hAsyTypInfo(AsyType)%nTotPin(1)
  nxp  = hAsyTypInfo(AsyType)%nPin
  nyp  = 2*nxp - 1
  
  CALL dmalloc(OutPut, nz, nxyp)
  
  DO ixyp = 1, nxyp
    ixy = Asy(iasy)%GlobalPinIdx(ixyp)
    
    IF(.NOT. FuelTh(ixy)%lfuel) CYCLE
    
    OutPut(1:nz, ixyp) = FuelTh(ixy)%tfuel(navg, 1:nz)
  END DO
  
  WRITE(io, 301) iasy, Asy(iasy)%ixa, Asy(iasy)%iya
  WRITE(io, *)
  
  CAll HexPrintTHResults(io, OutPut(1:nz, 1:nxyp), nxyp, nxp, nyp, nz)
  
  DEALLOCATE(OutPut)
END DO

NULLIFY (AsyInfo, Asy, FuelTH)
! ----------------------------------------------------

END SUBROUTINE HexPrintFuelAvgTemp
! ------------------------------------------------------------------------------------------------------------
!                                     03. PRINT : Fuel Cent Temp
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexPrintFuelCentTemp(io, Core, THInfo)

USE allocs
USE PARAM,   ONLY : ZERO
USE TYPEDEF, ONLY : CoreInfo_Type, ThInfo_Type, FuelTH_Type, AsyInfo_Type, Asy_Type
USE HexData, ONLY : hAsyTypInfo

IMPLICIT NONE

INTEGER :: io

TYPE(CoreInfo_Type) :: Core
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)

TYPE(ThInfo_Type) :: ThInfo
TYPE(FuelTH_Type), POINTER :: FuelTH(:)

REAL, POINTER :: OutPut(:, :)

INTEGER :: iz, ixy, ixyp, iasy, AsyType
INTEGER :: nz, nxya, nxp, nyp, nxyp
! ----------------------------------------------------

nz   = Core%nz
nxya = Core%nxya

AsyInfo => Core%AsyInfo
Asy     => Core%Asy
FuelTH  => THInfo%FuelTH

301 FORMAT('Assembly', I5, x, 'at', x, '(', I7, I7,' )')

WRITE(io, '(A)') '- Fuel Center-line Temperature -'
WRITE(io, *)

DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  
  nxyp = hAsyTypInfo(AsyType)%nTotPin(1)
  nxp  = hAsyTypInfo(AsyType)%nPin
  nyp  = 2*nxp - 1
  
  CALL dmalloc(OutPut, nz, nxyp)
  
  DO ixyp = 1, nxyp
    ixy = Asy(iasy)%GlobalPinIdx(ixyp)
    
    IF(.NOT. FuelTh(ixy)%lfuel) CYCLE
    
    OutPut(1:nz, ixyp) = FuelTh(ixy)%Tfuel(1, 1:nz)
  END DO
  
  WRITE(io, 301) iasy, Asy(iasy)%ixa, Asy(iasy)%iya
  WRITE(io, *)
  
  CAll HexPrintTHResults(io, OutPut, nxyp, nxp, nyp, nz)
  
  DEALLOCATE(OutPut)
END DO

NULLIFY (AsyInfo, Asy, FuelTH)
! ----------------------------------------------------

END SUBROUTINE HexPrintFuelCentTemp
! ------------------------------------------------------------------------------------------------------------
!                                     04. PRINT : Fuel Surf Temp
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexPrintFuelSurfTemp(io, Core, THInfo, THVar)

USE allocs
USE PARAM,   ONLY : ZERO
USE TYPEDEF, ONLY : CoreInfo_Type, ThInfo_Type, FuelTH_Type, THVar_Type, AsyInfo_Type, Asy_Type
USE HexData, ONLY : hAsyTypInfo

IMPLICIT NONE

INTEGER :: io

TYPE(CoreInfo_Type) :: Core
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)

TYPE(ThInfo_Type) :: ThInfo
TYPE(FuelTH_Type), POINTER :: FuelTH(:)
TYPE(THVar_Type) :: THVar

REAL, POINTER :: OutPut(:, :)

INTEGER :: iz, ixy, ixyp, iasy, AsyType
INTEGER :: nz, nxya, nxp, nyp, nxyp
! ----------------------------------------------------

nz   = Core%nz
nxya = Core%nxya

AsyInfo => Core%AsyInfo
Asy     => Core%Asy
FuelTH  => THInfo%FuelTH

301 FORMAT('Assembly', I5, x,'at' ,x,'(', I7, I7,' )')

WRITE(io, '(A)') '- Fuel Surface Temperature -'
WRITE(io, *) 

DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  
  nxyp = hAsyTypInfo(AsyType)%nTotPin(1)
  nxp  = hAsyTypInfo(AsyType)%nPin
  nyp  = 2*nxp - 1
  
  CALL dmalloc(OutPut, nz, nxyp)
  
  DO ixyp = 1, nxyp
    ixy = Asy(iasy)%GlobalPinIdx(ixyp)
    
    IF(.NOT. FuelTh(ixy)%lfuel) CYCLE
    
    OutPut(1:nz, ixyp) = FuelTh(ixy)%Tfuel(ThVar%npr1, 1:nz)
  END DO
  
  WRITE(io, 301) iasy, Asy(iasy)%ixa, Asy(iasy)%iya
  WRITE(io, *) 
  
  CAll HexPrintTHResults(io, OutPut, nxyp, nxp, nyp, nz)
  
  DEALLOCATE(OutPut)
END DO

NULLIFY (AsyInfo, Asy, FuelTH)
! ----------------------------------------------------

END SUBROUTINE HexPrintFuelSurfTemp
! ------------------------------------------------------------------------------------------------------------
!                                     05. PRINT : Coolant Temp
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexPrintCoolantTemp(io, Core, THInfo)

USE allocs
USE PARAM,   ONLY : ZERO
USE TYPEDEF, ONLY : CoreInfo_Type, ThInfo_Type, CoolantTh_Type, FuelTH_Type, AsyInfo_Type, Asy_Type
USE HexData, ONLY : hAsyTypInfo

IMPLICIT NONE

INTEGER :: io

TYPE(CoreInfo_Type) :: Core
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)

TYPE(ThInfo_Type) :: ThInfo
TYPE(CoolantTH_Type), POINTER :: CoolantTH(:)

REAL, POINTER :: OutPut(:, :)
INTEGER :: iz, ixy, ixyp, iasy, AsyType
INTEGER :: nz, nxy, nxya, nxp, nyp, nxyp
! ----------------------------------------------------

nz   = Core%nz
nxy  = Core%nxy
nxya = Core%nxya

AsyInfo   => Core%AsyInfo
Asy       => Core%Asy
CoolantTH => THInfo%CoolantTH

301 FORMAT('Assembly', I5, x,'at' ,x,'(', I7, I7,' )')

WRITE(io, '(A)') '- Coolant Temperature -'
WRITE(io, *) 

DO iasy = 1, nxya
  AsyType = Asy(iasy)%AsyType
  
  IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
  
  nxyp = hAsyTypInfo(AsyType)%nTotPin(1)
  nxp  = hAsyTypInfo(AsyType)%nPin
  nyp  = 2*nxp - 1
  
  CALL dmalloc(OutPut, nz, nxyp)
  
  DO ixyp = 1, nxyp
    ixy = Asy(iasy)%GlobalPinIdx(ixyp)
    
    OutPut(1:nz, ixyp) = THInfo%Tcool(1:nz, ixy)
  END DO
  
  WRITE(io, 301) iasy, Asy(iasy)%ixa, Asy(iasy)%iya
  WRITE(io, *) 
  
  CAll HexPrintTHResults(io, OutPut, nxyp, nxp, nyp, nz)
  
  DEALLOCATE(OutPut)
END DO

NULLIFY(AsyInfo, Asy, CoolantTH)
! ----------------------------------------------------

END SUBROUTINE HexPrintCoolantTemp
! ------------------------------------------------------------------------------------------------------------

END MODULE HexThOut