MODULE HexTLkg
  
USE TYPEDEF,   ONLY : AxFlx_TYPE, PinXS_TYPE, CoreInfo_TYPE, CmInfo_Type

IMPLICIT NONE

INTEGER, PARAMETER :: ConvLkg2ndMod  = 0
INTEGER, PARAMETER :: IntraLkg2ndMod = 1

CONTAINS
! ------------------------------------------------------------------------------------------------------------
!                                     01. Rad tLkg Update Type 0
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexRadTlkgUpdt(Core, CmInfo, TLKG, ixy1, ixy2, iz1, iz2, ng)

USE PARAM,   ONLY : ZERO
USE CNTL,    ONLY : nTracerCntl
USE HexType, ONLY : Type_HexCmfdPin
USE HexData, ONLY : hcPin, hPinInfo

IMPLICIT NONE

TYPE(CoreInfo_TYPE) :: Core
TYPE(CmInfo_TYPE) :: CmInfo
REAL, POINTER :: Tlkg(:, :, :)
INTEGER :: ixy1, ixy2, ng, iz1, iz2
! ----------------------------------------------------

INTEGER :: iz, iz0, ixy, jxy, ig, iNgh
INTEGER :: nxy, nNgh

REAL    :: RadLkg, Jnet, Dtil, Dhat, neighphi, MyPhi

REAL, POINTER :: PhiC(:, :, :)
REAL, POINTER :: Hz(:)

TYPE(PinXS_TYPE), POINTER :: PinXS(:, :)
TYPE(Type_HexCmfdPin), POINTER :: hcPin_Loc
! ----------------------------------------------------

PinXS  => CmInfo%PinXS
PhiC   => CmInfo%PhiFm
Hz     => Core%Hzfm
nxy     = Core%nxy

DO ig = 1, ng
  DO iz = iz1, iz2
    iz0 = Core%SubPlaneMap(iz)
    
    DO ixy = 1, nxy
      hcPin_Loc => hcPin(ixy)
      
      RadLkg = ZERO
      MyPhi  = PhiC(ixy, iz, ig)
      
      nNgh = hcPin_Loc%nNgh
      
      DO iNgh = 1, nNgh
        DHat = PinXS(ixy, iz0)%DHat(iNgh, ig)
        DTil = PinXS(ixy, iz0)%Dtil(iNgh, ig)
        
        !jxy      = hcPin_Loc%NghIdx(iNgh)
        NeighPhi = ZERO
        
        IF(jxy .GT. 0) THEN
          NeighPhi = PhiC(jxy, iz, ig)
        ENDIF
        
        Jnet   = (Dtil - Dhat) * MyPhi - (Dtil + Dhat) * NeighPhi
        RadLkg = RadLkg + Jnet
      END DO
      
      RadLkg = RadLkg * Hz(iz) / hPinInfo(ixy)%VolFm(iz0)
      
      TLKG(ixy, iz, ig) = RadLkg
    END DO
  END DO
END DO

NULLIFY(PhiC, PinXS, hcPin_Loc, hz)
! ----------------------------------------------------

END SUBROUTINE HexRadTlkgUpdt
! ------------------------------------------------------------------------------------------------------------
!                                     01. Stab Gap Pin Lkg Hex
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexStabGapPinLkg(Core, Flx0,  PE)

USE TYPEDEF, ONLY : CoreInfo_Type, AxFlx_Type, PE_Type

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(AxFlx_Type), POINTER :: Flx0(:, :)
TYPE(PE_TYPE) :: PE
! ----------------------------------------------------

RETURN

END SUBROUTINE HexStabGapPinLkg
! ------------------------------------------------------------------------------------------------------------

END MODULE HexTLkg