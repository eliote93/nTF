MODULE HexAsyConst
  
IMPLICIT NONE

CONTAINS

! ------------------------------------------------------------------------------------------------------------
!                                     01. HEX SET : Asy Bndy
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetAsyBndy()

USE PARAM,   ONLY : Pi, HALF, ZERO
USE HexUtil, ONLY : SetEqn, SetNewPt
USE HexData, ONLY : aoF2F, AsyEqn, AsyVtx, aoPch, hEps, hLgc, PI_3, PI_2

IMPLICIT NONE

INTEGER :: iBndy, iCor
REAL    :: Theta, Cnt(2)
! ----------------------------------------------------

DO iBndy = 1, 7
  Theta = PI_2 - PI_3 * (iBndy - 1)
  
  AsyVtx(1:2, iBndy) = SetNewPt(aoPch, Theta)
END DO

Cnt = ZERO

DO iBndy = 1, 6
  AsyEqn(1:3, iBndy) = SetEqn(AsyVtx(1:2, iBndy), AsyVtx(1:2, iBndy + 1), Cnt)
END DO
! ----------------------------------------------------

END SUBROUTINE HexSetAsyBndy
! ------------------------------------------------------------------------------------------------------------
!                                     02. SET : Asy Geo Typ
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetGeoTyp()

USE PARAM,   ONLY : HALF, PI, ZERO
USE HexType, ONLY : Type_HexGeoTypInfo
USE HexData, ONLY : hGeoTypInfo, aoF2F, aoPch, Sq3, AsyVtx, AsyEqn, hEps, PI_3, PI_6, PI_2
USE HexUtil, ONLY : SetEqn, SetNewPt, ChkSameVal

IMPLICIT NONE

INTEGER :: iGeo, iBndy

REAL :: Area, Ang
REAL :: Vtx(2, 6), Cnt(2)

TYPE(Type_HexGeoTypInfo), POINTER :: gInf_Loc
! ----------------------------------------------------

ALLOCATE (hGeoTypInfo (7)) ! Can be extended to 13 for 360 & REF

Area = aoF2F * aoF2F * Sq3 * HALF
! ----------------------------------------------------
!               01. Geo Typ = 1
! ----------------------------------------------------
gInf_Loc => hGeoTypInfo(1)

gInf_Loc%nBndy = 6
gInf_Loc%Area  = Area

gInf_Loc%Vtx(1:2, 1:7) = AsyVtx(1:2, 1:7)
! ----------------------------------------------------
!               02. Geo Typ = 2 ~ 4
! ----------------------------------------------------
Ang = PI * (5._8 / 6._8) ! Normal angle

DO iGeo = 2, 4
  gInf_Loc => hGeoTypInfo(iGeo)
  
  gInf_Loc%nBndy = 5
  gInf_Loc%Area  = Area * HALF
  
  CALL HexSetGeoBndy_180(Ang, gInf_Loc%Vtx)
  
  gInf_Loc%Cnt = SetNewPt(aoPch * HALF, Ang)
  
  Ang = Ang + 2._8 * PI_3
END DO
! ----------------------------------------------------
!               03. Geo Typ = 5 ~ 7
! ----------------------------------------------------
Ang = PI_2 ! Normal angle

DO iGeo = 5, 7
  gInf_Loc => hGeoTypInfo(iGeo)
  
  gInf_Loc%nBndy = 4
  gInf_Loc%Area  = Area / 6._8
  
  CALL HexSetGeoBndy_060(Ang, gInf_Loc%Vtx)
  
  gInf_Loc%Cnt = SetNewPt(aoPch * HALF, Ang)
  
  Ang = Ang + 2._8 * PI_3
END DO
! ----------------------------------------------------
!               04. SET : Bndy Eqn & Cor
! ----------------------------------------------------
Cnt = ZERO

DO iGeo = 1, 7
  gInf_Loc => hGeoTypInfo(iGeo)
  
  DO iBndy = 1, gInf_Loc%nBndy
    gInf_Loc%Eqn(1:3, iBndy) = SetEqn(gInf_Loc%Vtx(1:2, iBndy), gInf_Loc%Vtx(1:2, iBndy + 1), Cnt)
    
    IF (ChkSameVal(gInf_Loc%Eqn(1, iBndy), ZERO)) gInf_Loc%Cor(iBndy) = 1
  END DO
END DO

NULLIFY (gInf_Loc)

CONTAINS
! ----------------------------------------------------
!               SUB. HEX SET : Geo Bndy 180
! ----------------------------------------------------
SUBROUTINE HexSetGeoBndy_180(Ang, Vtx)

IMPLICIT NONE

REAL :: Ang, Vtx(2, 7)

INTEGER :: iBndy
REAL :: Tmp
! ----------------------------

Tmp = Ang - PI_2

Vtx(1:2, 1) = SetNewPt(aoF2F * HALF, Tmp)

Tmp = Tmp + PI_6

DO iBndy = 2, 4
  Vtx(1:2, iBndy) = SetNewPt(aoPch, Tmp)
  
  Tmp = Tmp + PI_3
END DO

Tmp = Tmp - PI_6

Vtx(1:2, 5) = SetNewPt(aoF2F * HALF, Tmp)
Vtx(1:2, 6) = Vtx(1:2, 1)

END SUBROUTINE HexSetGeoBndy_180
! ----------------------------------------------------
!               SUB. HEX SET : Geo Bndy 060
! ----------------------------------------------------
SUBROUTINE HexSetGeoBndy_060(Ang, Vtx)

IMPLICIT NONE

REAL :: Ang, Vtx(2, 7)

INTEGER :: iBndy
REAL :: Tmp
! ----------------------------

Vtx(1:2, 1) = SetNewPt(aoF2F * HALF, Ang - PI_6)
Vtx(1:2, 2) = SetNewPt(aoPch,        Ang)
Vtx(1:2, 3) = SetNewPt(aoF2F * HALF, Ang + PI_6)
Vtx(1:2, 4) = ZERO
Vtx(1:2, 5) = Vtx(1:2, 1)

END SUBROUTINE HexSetGeoBndy_060
! ----------------------------------------------------

END SUBROUTINE HexSetGeoTyp
! ------------------------------------------------------------------------------------------------------------
!                                     03. SET : Asy Loc
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetAsyLoc()

USE PARAM,   ONLY : ZERO, HALF
USE ioutil,  ONLY : terminate
USE HexData, ONLY : nAsyCore, nhAsy, aoF2F, aoPch, hAsy, Asy1Dto2DMap, hCore, hLgc, &
                    ncBndy, cBndyPt, cBndyEq, AsyVtx, AsyEqn, Sq3, nGeoTyp
USE HexUtil, ONLY : SetEqn, ChkPtEqn

IMPLICIT NONE

INTEGER :: iaX, iaY, jaX, jaY, caX, caY, iAsy, iaTyp, iBndy, iTmp
REAL    :: xDel, yDel, Tmp, Cnt(2)
! ----------------------------------------------------

ALLOCATE (hAsy (nhAsy))

xDel = aoF2F
yDel = aoPch * 1.5_8

SELECT CASE (hLgc%iSym)
CASE (3) ! 360
  caX = nAsyCore
  caY = nAsyCore
CASE DEFAULT
  caX = 1
  caY = 1
END SELECT
! ----------------------------------------------------
!               01. SET : Asy Loc
! ----------------------------------------------------
DO iAsy = 1, nhAsy
  iaX   = Asy1Dto2DMap(1, iAsy)
  iaY   = Asy1Dto2DMap(2, iAsy)
  jaX   = iaX - caX
  jaY   = iaY - caY
  iaTyp = hCore(iaX, iaY)
  
  hAsy(iAsy)%iaX    = iaX
  hAsy(iAsy)%iaY    = iaY
  hAsy(iAsy)%AsyTyp = iaTyp
  
  hAsy(iAsy)%Cnt(1) =  xDel * jaX - xDel * 0.5 * jaY
  hAsy(iAsy)%Cnt(2) = -yDel * jaY
END DO
! ----------------------------------------------------
!               02. SET : Core Bndy
! ----------------------------------------------------
nGeoTyp = 1

SELECT CASE (hLgc%iSym)
! ----------------------------
!      CASE : Sng Asy / Cel
! ----------------------------
CASE (4:5)
  ncBndy = 6
  
  cBndyPt(1:2, 1:7) = AsyVtx(1:2, 1:7)
  cBndyEq(1:3, 1:6) = AsyEqn(1:3, 1:6)
! ----------------------------
!      CASE : 360
! ----------------------------
CASE (3)
  IF (hLgc%lRadRef) CALL terminate("UNDER-CONSTRUCTION")

  ncBndy = 0
! ----------------------------
!      CASE : 120
! ----------------------------
CASE (2)
  CALL terminate("UNDER-CONSTRUCTION")
! ----------------------------
!      CASE : 060
! ----------------------------
CASE (1)
  Tmp = aoF2F * (nAsyCore - 1)
  
  cBndyPt(1, 1) =  Tmp
  cBndyPt(2, 1) =  ZERO
  cBndyPt(1, 2) =  ZERO
  cBndyPt(2, 2) =  ZERO
  cBndyPt(1, 3) =  Tmp * HALF
  cBndyPt(2, 3) = -Tmp * HALF * Sq3
  cBndyPt(1, 4) =  cBndyPt(1, 1)
  cBndyPt(2, 4) =  cBndyPt(2, 1)
  
  nGeoTyp = 7
  
  IF (hLgc%lRadRef) THEN
    ncBndy = 3
  ELSE
    ncBndy = 2
  END IF
  
  Cnt = ZERO
  
  DO iBndy = 1, ncBndy
    cBndyEq(1:3, iBndy) = SetEqn(cBndyPt(1:2, iBndy), cBndyPt(1:2, iBndy + 1), Cnt)
  END DO
END SELECT
! ----------------------------------------------------
!               03. SET : Geo Typ
! ----------------------------------------------------
IF (hLgc%iSym > 3) RETURN
! ----------------------------
!      1. 360
! ----------------------------
IF (hLgc%l360) THEN
  IF (hLgc%lRadRef) CALL terminate("UNDER-CONSTRUCTION")
  
  RETURN
! ----------------------------
!      2. 120
! ----------------------------
ELSE IF (hLgc%l120) THEN
  CALL terminate("UNDER-CONSTRUCTION")
! ----------------------------
!      3. 060
! ----------------------------
ELSE
  DO iAsy = 1, nhAsy
    iTmp = 0
    
    DO iBndy = 1, ncBndy
      IF (ChkPtEqn(hAsy(iAsy)%Cnt(1:2), cBndyEq(1:3, iBndy))) THEN
        iTmp = iTmp + 10 ** (iBndy - 1)
      END IF
    END DO
    
    SELECT CASE (iTmp)
    CASE (  0)
    CASE (001); hAsy(iAsy)%GeoTyp = 3 !  011 --- 001 --- 101
    CASE (010); hAsy(iAsy)%GeoTyp = 4 !     \           /
    CASE (011); hAsy(iAsy)%GeoTyp = 7 !     010  000  100
    CASE (100); hAsy(iAsy)%GeoTyp = 2 !        \     /
    CASE (101); hAsy(iAsy)%GeoTyp = 6 !          110
    CASE (110); hAsy(iAsy)%GeoTyp = 5
    CASE DEFAULT
      CALL terminate("SET ASY GEO TYP")
    END SELECT
    
    SELECT CASE (hAsy(iAsy)%GeoTyp)
    CASE (2:4); hAsy(iAsy)%wt = HALF
    CASE (5:7); hAsy(iAsy)%wt = 1._8/6._8
    END SELECT
  END DO
END IF
! ----------------------------------------------------

END SUBROUTINE HexSetAsyLoc
! ------------------------------------------------------------------------------------------------------------
!                                     04. HEX SET : Asy Pin Num
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetAsyPinNum(iAsy)

USE HexType, ONLY : Type_HexAsy, Type_HexAsyTypInfo
USE HexData, ONLY : hAsy, hAsyTypInfo, nHexPin

IMPLICIT NONE

INTEGER :: iPin, iAsy, iGeoTyp

TYPE(Type_HexAsy), POINTER :: hAsy_Loc
TYPE(Type_HexAsyTypInfo), POINTER :: aInf_Loc
! ----------------------------------------------------

hAsy_Loc => hAsy(iAsy)
aInf_Loc => hAsyTypInfo(hAsy_Loc%AsyTyp)

iGeoTyp = hAsy_Loc%GeoTyp

hAsy_Loc%PinIdxSt = nHexPin + 1
hAsy_Loc%nRodPin  = aInf_Loc%nRodPin(iGeoTyp)
hAsy_Loc%nTotPin  = aInf_Loc%nTotPin(iGeoTyp)

nHexPin = nHexPin + hAsy_Loc%nTotPin

NULLIFY (hAsy_Loc, aInf_Loc)
! ----------------------------------------------------

END SUBROUTINE HexSetAsyPinNum
! ------------------------------------------------------------------------------------------------------------

END MODULE HexAsyConst