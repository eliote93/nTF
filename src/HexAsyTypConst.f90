MODULE HexAsyTypConst
  
IMPLICIT NONE

CONTAINS

! ------------------------------------------------------------------------------------------------------------
!                                     01. HEX SET : Asy Typ Pin Map
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetAsyTypPinMap(iAsyTyp)

USE allocs
USE PARAM,   ONLY : HALF, PI, ZERO
USE HexType, ONLY : Type_HexAsyTypInfo
USE HexData, ONLY : hAsyTypInfo, aoF2F, hEps, hLgc, PI_3, PI_2
USE HexUtil, ONLY : SetEqn, ChkArrayZero

IMPLICIT NONE

INTEGER :: iAsyTyp

INTEGER :: nPin, nRng, nRod, nTot, nTmp
INTEGER :: ix, iy, jx, jy, iPin, jPin, kPin, iBndy, iDir, iCor
INTEGER :: mp(1:2) = (/-1, 1/)

REAL :: dx, dy, rr, Ang
REAL :: Cnt(2), Tmp(2)

TYPE(Type_HexAsyTypInfo), POINTER :: aInf_Loc
! ----------------------------------------------------

aInf_Loc => hAsyTypInfo(iAsyTyp)

nPin = aInf_Loc%nPin
nRng = 2 * nPin - 1
nRod = 3 * nPin * (nPin - 1) + 1
nTot = 3 * nPin * (nPin - 1) + 1 + 12 * (nPin - 1)

aInf_Loc%nRodPin(1) = nRod
aInf_Loc%nTotPin(1) = nTot

CALL dmalloc(aInf_Loc%Pin1Dto2Dmap, 2, nRod)

CALL dmalloc0(aInf_Loc%Pin2Dto1Dmap, 0, nRng+1, 0, nRng+1)
! ----------------------------------------------------
!               01. Upper Part
! ----------------------------------------------------
nTmp = nPin
iPin = 0

DO iy = 1, nPin
  DO ix = 1, nTmp
    iPin = iPin + 1
    
    aInf_Loc%Pin1Dto2Dmap(1, iPin) = ix
    aInf_Loc%Pin1Dto2Dmap(2, iPin) = iy
    
    aInf_Loc%Pin2Dto1Dmap(ix, iy) = iPin
  END DO
  
  nTmp = nTmp + 1
END DO
! ----------------------------------------------------
!               02. Lower Part
! ----------------------------------------------------
nTmp = 2

DO iy = nPin + 1, nRng
  DO ix = nTmp, nRng
    iPin = iPin + 1
    
    aInf_Loc%Pin1Dto2Dmap(1, iPin) = ix
    aInf_Loc%Pin1Dto2Dmap(2, iPin) = iy
    
    aInf_Loc%Pin2Dto1Dmap(ix, iy) = iPin
  END DO
  
  nTmp = nTmp + 1
END DO
! ----------------------------------------------------
!               03. SET : Rod Pin Cnt
! ----------------------------------------------------
CALL dmalloc(aInf_Loc%PinCnt, 2, nTot)

dx = aInf_Loc%pPch * 1.5_8
dy = aInf_Loc%pF2F * 0.5_8

DO iPin = 1, nRod
  jx = aInf_Loc%Pin1Dto2Dmap(1, iPin) - nPin
  jy = aInf_Loc%Pin1Dto2Dmap(2, iPin) - nPin
  
  aInf_Loc%PinCnt(1, iPin) =  dx * (jx - jy)
  aInf_Loc%PinCnt(2, iPin) = -dy * (jx + jy)
END DO
! ----------------------------------------------------
!               04. SET : Gap Pin Cnt
! ----------------------------------------------------
IF (hLgc%lSngCel) RETURN

aInf_Loc%gHgt = (aoF2F - aInf_Loc%aiF2F) * HALF

dx  = aInf_Loc%pF2F * HALF
rr  = (aInf_Loc%aiF2F + aInf_Loc%gHgt) * HALF
Ang = PI_3

DO iBndy = 1, 6
  Cnt(1) = rr * cos(Ang)
  Cnt(2) = rr * sin(Ang)
  
  jPin = nRod + 2 * (nPin - 1) * (iBndy - 1) + nPin
  
  DO iDir = 1, 2
    Tmp(1) = Cnt(1) + mp(iDir) * dx * HALF * cos(Ang - PI_2)
    Tmp(2) = Cnt(2) + mp(iDir) * dx * HALF * sin(Ang - PI_2)
    
    DO iPin = 1, nPin - 1
      kPin = jPin + mp(iDir) * iPin
      
      aInf_Loc%PinCnt(1, kPin) = Tmp(1) + mp(iDir) * dx * (iPin - 1) * cos(Ang - PI_2)
      aInf_Loc%PinCnt(2, kPin) = Tmp(2) + mp(iDir) * dx * (iPin - 1) * sin(Ang - PI_2)
    END DO
    
    jPin = jPin - 1
  END DO
  
  Ang = Ang - PI_3
END DO

CALL ChkArrayZero(aInf_Loc%PinCnt, 1, 2, 1, nTot)

NULLIFY (aInf_Loc)
! ----------------------------------------------------

END SUBROUTINE HexSetAsyTypPinMap
! ------------------------------------------------------------------------------------------------------------
!                                     02. HEX SET : Asy Typ Pin Loc Idx
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetAsyTypPinLocIdx(iAsyTyp)

USE allocs
USE PARAM,   ONLY : TRUE, FALSE, ZERO
USE HexType, ONLY : Type_HexAsyTypInfo, Type_HexGeoTypInfo
USE HexData, ONLY : hAsyTypInfo, hGeoTypInfo, nGeoTyp
USE HexUtil, ONLY : ChkPtEqn, CalPtLineSgn

IMPLICIT NONE

INTEGER :: iAsyTyp, iGeo, iPin, jPin, nRod, iBndy

REAL :: gCnt(2), Cnt(2), gEqn(3)

TYPE(Type_HexAsyTypInfo), POINTER :: aInf_Loc
TYPE(Type_HexGeoTypInfo), POINTER :: gInf_Loc
! ----------------------------------------------------

aInf_Loc => hAsyTypInfo(iAsyTyp)

CALL dmalloc(aInf_Loc%PinLocIdx, nGeoTyp, aInf_Loc%nTotPin(1))
CALL dmalloc(aInf_Loc%lGeoPin,   nGeoTyp, aInf_Loc%nTotPin(1))
! ----------------------------------------------------
!               01. SET : l Geo Pin
! ----------------------------------------------------
aInf_Loc%lGeoPin = TRUE

DO iGeo = 2, nGeoTyp
  gInf_Loc => hGeoTypInfo(iGeo)
  
  gCnt = gInf_Loc%Cnt
  
  DO iPin = 1, aInf_Loc%nTotPin(1)
    Cnt = aInf_Loc%PinCnt(1:2, iPin)
    
    DO iBndy = 1, gInf_Loc%nBndy
      gEqn = gInf_Loc%Eqn(1:3, iBndy)
      
      IF (ChkPtEqn    (Cnt, gEqn)) CYCLE
      IF (CalPtLineSgn(Cnt, gEqn, gCnt) > ZERO) CYCLE
      
      aInf_Loc%lGeoPin(iGeo, iPin) = FALSE
    END DO
  END DO
END DO
! ----------------------------------------------------
!               02. SET : Pin Loc Idx
! ----------------------------------------------------
aInf_Loc%PinLocIdx = -1

DO iGeo = 1, nGeoTyp
  jPin = 0
  nRod = 0
  
  DO iPin = 1, aInf_Loc%nTotPin(1)
    IF (aInf_Loc%lGeoPin(iGeo, iPin)) THEN
      jPin = jPin + 1
      
      aInf_Loc%PinLocIdx(iGeo, iPin) = jPin
      
      IF (iPin <= aInf_Loc%nRodPin(1)) nRod = nRod + 1
    END IF
  END DO
  
  aInf_Loc%nRodPin(iGeo) = nRod
  aInf_Loc%nTotPin(iGeo) = jPin
END DO

NULLIFY (aInf_Loc, gInf_Loc)
! ----------------------------------------------------

END SUBROUTINE HexSetAsyTypPinLocIdx
! ------------------------------------------------------------------------------------------------------------

END MODULE HexAsyTypConst