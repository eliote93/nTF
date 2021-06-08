! ------------------------------------------------------------------------------------------------------------
!                                     HEX SET : Asy Typ MP Ngh
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetAsyTypMPngh(iaTyp)

USE allocs
USE HexType, ONLY : Type_HexAsyTypInfo
USE HexData, ONLY : hAsyTypInfo, nGeoTyp, mpTypNumNgh

IMPLICIT NONE

INTEGER :: iaTyp
INTEGER :: iGeo, iPin, ivTyp, iBndy, nTot, nBndy

TYPE(Type_HexAsyTypInfo), POINTER :: aInf_Loc
! ----------------------------------------------------

aInf_Loc => hAsyTypInfo(iaTyp)

nTot = aInf_Loc%nTotPin(1)

CALL dmalloc(aInf_Loc%cnpSlfMPnum,       nGeoTyp, nTot)
CALL dmalloc(aInf_Loc%cnpSlfMPidx,    3, nGeoTyp, nTot)
CALL dmalloc(aInf_Loc%cnpSufMPidx, 2, 6, nGeoTyp, nTot)
CALL dmalloc(aInf_Loc%cnpSufMPsuf, 2, 6, nGeoTyp, nTot)
CALL dmalloc(aInf_Loc%cnpSufMPnum,    6, nGeoTyp, nTot)

aInf_Loc%cnpSlfMPnum = 1 ! Self
! ----------------------------------------------------
!               01. SET : Self Data
! ----------------------------------------------------
DO iGeo = 1, nGeoTyp
  DO iPin = 1, nTot
    IF (.NOT. aInf_Loc%lGeoPin(iGeo, iPin)) CYCLE
    
    ivTyp = aInf_Loc%PinVtxTyp(iGeo, iPin)
    nBndy = mpTypNumNgh(ivTyp)
    
    aInf_Loc%cnpSlfMPidx(1, iGeo, iPin) = iPin
    
    DO iBndy = 1, nBndy
      aInf_Loc%cnpSufMPnum(   iBndy, iGeo, iPin) = 1
      aInf_Loc%cnpSufMPidx(1, iBndy, iGeo, iPin) = iPin
      aInf_Loc%cnpSufMPsuf(1, iBndy, iGeo, iPin) = iBndy
    END DO
    
  END DO
END DO

NULLIFY (aInf_Loc)
! ----------------------------------------------------

END SUBROUTINE HexSetAsyTypMPngh