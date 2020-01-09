! ------------------------------------------------------------------------------------------------------------
!                                     HEX SET : Asy Typ MP Ngh
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetAsyTypMPngh(iaTyp)

USE HexType, ONLY : Type_HexAsyTypInfo
USE HexData, ONLY : hAsyTypInfo, nGeoTyp, mpTypNumNgh

IMPLICIT NONE

INTEGER :: iaTyp
INTEGER :: iGeo, iPin, ivTyp, iBndy, nTot, nBndy

TYPE(Type_HexAsyTypInfo), POINTER :: aInf_Loc
! ----------------------------------------------------

aInf_Loc => hAsyTypInfo(iaTyp)

nTot = aInf_Loc%nTotPin(1)

ALLOCATE (aInf_Loc%cpSlfMPnum       (nGeoTyp, nTot)); aInf_Loc%cpSlfMPnum = 1
ALLOCATE (aInf_Loc%cpSlfMPidx    (3, nGeoTyp, nTot)); aInf_Loc%cpSlfMPidx = 0
ALLOCATE (aInf_Loc%cpSufMPidx (2, 6, nGeoTyp, nTot)); aInf_Loc%cpSufMPidx = 0
ALLOCATE (aInf_Loc%cpSufMPsuf (2, 6, nGeoTyp, nTot)); aInf_Loc%cpSufMPsuf = 0
ALLOCATE (aInf_Loc%cpSufMPnum    (6, nGeoTyp, nTot)); aInf_Loc%cpSufMPnum = 0
! ----------------------------------------------------
!               01. SET : Self Data
! ----------------------------------------------------
DO iGeo = 1, nGeoTyp
  DO iPin = 1, nTot
    IF (.NOT. aInf_Loc%lGeoPin(iGeo, iPin)) CYCLE
    
    ivTyp = aInf_Loc%PinVtxTyp(iGeo, iPin)
    nBndy = mpTypNumNgh(ivTyp)
    
    aInf_Loc%cpSlfMPidx(1, iGeo, iPin) = iPin
    
    DO iBndy = 1, nBndy
      aInf_Loc%cpSufMPnum(   iBndy, iGeo, iPin) = 1
      aInf_Loc%cpSufMPidx(1, iBndy, iGeo, iPin) = iPin
      aInf_Loc%cpSufMPsuf(1, iBndy, iGeo, iPin) = iBndy
    END DO
    
  END DO
END DO

NULLIFY (aInf_Loc)
! ----------------------------------------------------

END SUBROUTINE HexSetAsyTypMPngh