MODULE HexVtx

IMPLICIT NONE

CONTAINS

! ------------------------------------------------------------------------------------------------------------
!                                     01. HEX SET : Asy Typ Pin Vtx
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetAsyTypVtxTyp(iAsyTyp)

USE PARAM,   ONLY : ZERO, PI, HALF
USE HexType, ONLY : Type_HexAsyTypInfo
USE HexData, ONLY : hLgc, hAsyTypInfo, Sq3Inv, nGeoTyp, PI_3, PI_6, mpTypNumNgh, spTypNumNgh
USE HexUtil, ONLY : ChkArrayZero, SetNewPt, RotPt, FindArea, FindCnt, FindPtLgh, FindLinePtLgh, SetEqn

IMPLICIT NONE

INTEGER :: iAsyTyp, iBndy, ivTyp, nBndy

REAL :: Ang, Hgt, Lgh, gHgt, pF2F, pPch, Del(2), Eqn(3), Cnt(2)

REAL :: mpTypPts(2, 7, 10) = ZERO
REAL :: spTypPts(2, 7,  7) = ZERO

TYPE(Type_HexAsyTypInfo), POINTER :: aInf_Loc
! ----------------------------------------------------

IF (hLgc%lSngCel) RETURN

aInf_Loc => hAsyTypInfo(iAsyTyp)

gHgt = aInf_Loc%gHgt
pF2F = aInf_Loc%pF2F
pPch = aInf_Loc%pPch
! ----------------------------------------------------
!               01. SET : Inn Pts
! ----------------------------------------------------
DO iBndy = 1, 7
  Ang = PI_3 * (2 - iBndy)
  
  mpTypPts(1:2, iBndy, 3) = SetNewPt(pPch, Ang) ! 360
END DO

! 060 & Geo = 5 (NN)
mpTypPts(1:2, 1, 1) = mpTypPts(1:2, 6, 3)
mpTypPts(1:2, 2, 1) = mpTypPts(1:2, 1, 3)
mpTypPts(1:2, 3, 1) = ZERO
mpTypPts(1:2, 4, 1) = mpTypPts(1:2, 1, 1)

! 180 & Geo = 2 (NW)
mpTypPts(1:2, 1, 2) = mpTypPts(1:2, 4, 3)
mpTypPts(1:2, 2, 2) = mpTypPts(1:2, 5, 3)
mpTypPts(1:2, 3, 2) = mpTypPts(1:2, 6, 3)
mpTypPts(1:2, 4, 2) = mpTypPts(1:2, 1, 3)
mpTypPts(1:2, 5, 2) = mpTypPts(1:2, 4, 3)
! ----------------------------------------------------
!               02. SET : Bndy Pts
! ----------------------------------------------------
Hgt  = HALF * (aInf_Loc%aiF2F - 3 * (aInf_Loc%nPin - 1) * pPch)
Lgh  = HALF * (aInf_Loc%aiPch -     (aInf_Loc%nPin - 2) * pF2F)
Ang  = -PI_6

! Bndy Typ = 1 & Dir = 1 (NN)
mpTypPts(1, 5, 4) = ZERO
mpTypPts(2, 5, 4) = Hgt * 2._8 * Sq3Inv

mpTypPts(1:2, 1, 4) = mpTypPts(1:2, 5, 4) + SetNewPt(Lgh, Ang)
mpTypPts(1:2, 4, 4) = mpTypPts(1:2, 5, 4) + SetNewPt(Lgh, PI - Ang)
mpTypPts(1:2, 2, 4) = mpTypPts(1:2, 3, 3) 
mpTypPts(1:2, 3, 4) = mpTypPts(1:2, 4, 3)
mpTypPts(1:2, 6, 4) = mpTypPts(1:2, 1, 4)

! Bndy Typ = 2 & Dir = 7 (NE)
mpTypPts(1:2, 2, 5) = mpTypPts(1:2, 3, 3)
mpTypPts(1:2, 3, 5) = mpTypPts(1:2, 4, 3)
mpTypPts(1:2, 4, 5) = mpTypPts(1:2, 5, 3)
mpTypPts(1:2, 1, 5) = mpTypPts(1:2, 1, 4)
mpTypPts(1:2, 5, 5) = mpTypPts(1:2, 1, 4)- SetNewPt(pF2F, Ang)
mpTypPts(1:2, 6, 5) = mpTypPts(1:2, 1, 4)- SetNewPt(pF2F, Ang) * HALF
mpTypPts(1:2, 7, 5) = mpTypPts(1:2, 1, 5)

! Bndy Typ = 2 & Dir = 7 (NE) & Spt
mpTypPts(1:2,   1, 6) = mpTypPts(1:2,   6, 5)
mpTypPts(1:2, 2:5, 6) = mpTypPts(1:2, 3:6, 5)

! Bndy Typ = 2 & Dir = 10 (SW) & Spt
mpTypPts(1:2, 2, 7) = mpTypPts(1:2, 1, 3)
mpTypPts(1:2, 3, 7) = mpTypPts(1:2, 6, 3)
mpTypPts(1:2, 4, 7) = RotPt(mpTypPts(1:2, 1, 5), PI)
mpTypPts(1:2, 1, 7) = RotPt(mpTypPts(1:2, 6, 5), PI)
mpTypPts(1:2, 5, 7) = mpTypPts(1:2, 1, 7)
! ----------------------------------------------------
!               03. SET : Gap Pts
! ----------------------------------------------------
! Gap 02 & Dir = 2 (Clock-wise)
mpTypPts(1, 1:2, 10) =  pF2F / 4._8
mpTypPts(1, 3:4, 10) = -pF2F / 4._8
mpTypPts(2, 1:4, 10) =  gHgt * HALF
mpTypPts(2, 2:3, 10) = -gHgt * HALF

DO iBndy = 1, 4
  mpTypPts(1:2, iBndy, 10) = RotPt(mpTypPts(1:2, iBndy, 10), Ang)
END DO

! Gap 01 & Dir = 1 (NE & Right)
mpTypPts(1:2, 3, 8) = mpTypPts(1:2, 3, 10)
mpTypPts(1:2, 4, 8) = mpTypPts(1:2, 4, 10)
mpTypPts(1:2, 1, 8) = mpTypPts(1:2, 4, 8) + SetNewPt(gHgt * Sq3Inv + Lgh, Ang)
mpTypPts(1:2, 2, 8) = mpTypPts(1:2, 3, 8) + SetNewPt(                Lgh, Ang)

! Gap 01 & Dir = 1 (NE & Left)
mpTypPts(1:2, 3, 9) = mpTypPts(1:2, 2, 10)
mpTypPts(1:2, 4, 9) = mpTypPts(1:2, 1, 10)
mpTypPts(1:2, 1, 9) = mpTypPts(1:2, 1, 10) - SetNewPt(gHgt * Sq3Inv + Lgh, Ang)
mpTypPts(1:2, 2, 9) = mpTypPts(1:2, 2, 10) - SetNewPt(                Lgh, Ang)

mpTypPts(1:2, 5, 8:10) = mpTypPts(1:2, 1, 8:10)

CALL ChkArrayZero(mpTypPts, 1, 2, 1, 7, 1, 10)

aInf_Loc%mpTypPts = mpTypPts
! ----------------------------------------------------
!               04. SET : Super-Pin Typ Pts
! ----------------------------------------------------
Ang = PI_3
Del = SetNewPt(gHgt, Ang)

spTypPts(:, 1:7, 1:7) = mpTypPts(:, 1:7, 1:7)

spTypPts(1:2, 4, 4) = mpTypPts(1:2, 4, 4) + SetNewPt(gHgt, Ang * 2._8)
spTypPts(1:2, 5, 4) = mpTypPts(1:2, 5, 4) + [ZERO, 2._8 * gHgt * Sq3Inv]
spTypPts(1:2, 6, 4) = mpTypPts(1:2, 1, 4) + Del(1:2)
spTypPts(1:2, 1, 4) = spTypPts(1:2, 6, 4)

spTypPts(1:2, 5, 5) = mpTypPts(1:2, 5, 5) + Del(1:2)
spTypPts(1:2, 6, 5) = mpTypPts(1:2, 1, 5) + Del(1:2)
spTypPts(1:2, 1, 5) = spTypPts(1:2, 6, 5)

spTypPts(1:2, 4, 6) = mpTypPts(1:2, 4, 6) + Del(1:2)
spTypPts(1:2, 5, 6) = mpTypPts(1:2, 1, 6) + Del(1:2)
spTypPts(1:2, 1, 6) = spTypPts(1:2, 5, 6)

spTypPts(1:2, 4, 7) = mpTypPts(1:2, 4, 7) - Del(1:2)
spTypPts(1:2, 5, 7) = mpTypPts(1:2, 1, 7) - Del(1:2)
spTypPts(1:2, 1, 7) = spTypPts(1:2, 5, 7)

CALL ChkArrayZero(spTypPts, 1, 2, 1, 7, 1, 7)

aInf_Loc%spTypPts = spTypPts
! ----------------------------------------------------
!               05. SET : Area / Lgh / C2B
! ----------------------------------------------------
DO ivTyp = 1, 10
  nBndy = mpTypNumNgh(ivTyp)
  Cnt   = FindCnt(nBndy, mpTypPts(1:2, 1:nBndy, ivTyp))
  
  aInf_Loc%mpTypAre(ivTyp) = FindArea(nBndy, mpTypPts(1:2, 1:nBndy, ivTyp))
  
  DO iBndy = 1, nBndy
    Eqn = SetEqn(mpTypPts(1:2, iBndy, ivTyp), mpTypPts(1:2, iBndy+1, ivTyp), Cnt)
    
    aInf_Loc%mpBndyLgh(iBndy, ivTyp) = FindPtLgh(mpTypPts(1:2, iBndy, ivTyp), mpTypPts(1:2, iBndy+1, ivTyp))
    aInf_Loc%mpBndyC2B(iBndy, ivTyp) = FindLinePtLgh(Eqn, Cnt)
  END DO
END DO

DO ivTyp = 1, 7
  nBndy = spTypNumNgh(ivTyp)
  Cnt   = FindCnt(nBndy, spTypPts(1:2, 1:nBndy, ivTyp))
  
  aInf_Loc%spTypAre(ivTyp) = FindArea(nBndy, spTypPts(1:2, 1:nBndy, ivTyp))
  
  DO iBndy = 1, nBndy
    Eqn = SetEqn(spTypPts(1:2, iBndy, ivTyp), spTypPts(1:2, iBndy+1, ivTyp), Cnt)
    
    aInf_Loc%spBndyLgh(iBndy, ivTyp) = FindPtLgh(spTypPts(1:2, iBndy, ivTyp), spTypPts(1:2, iBndy+1, ivTyp))
    aInf_Loc%spBndyC2B(iBndy, ivTyp) = FindLinePtLgh(Eqn, Cnt)
  END DO
END DO

NULLIFY (aInf_Loc)
! ----------------------------------------------------

END SUBROUTINE HexSetAsyTypVtxTyp
! ------------------------------------------------------------------------------------------------------------
!                                     02. HEX SET : Asy Typ Pin Vtx
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetAsyTypPinVtx(iAsyTyp)

USE allocs
USE PARAM,   ONLY : ZERO, PI, HALF
USE HexType, ONLY : Type_HexAsyTypInfo, Type_HexGeoTypInfo
USE HexData, ONLY : hAsyTypInfo, hGeoTypInfo, nGeoTyp, Sq3, PI_2, PI_3, PI_6, mpTypNumNgh, spTypNumNgh, hLgc
USE HexUtil, ONLY : ChkArrayZero, RotPt, SetSgn_REAL, ChkSameVal, ChkPtEqn, FindPtAng
USE HexTst,  ONLY : HexTstPinVtx, HexTstSpVtx

IMPLICIT NONE

INTEGER :: iAsyTyp, iPin, jMod, iGeo, iBndy, ivTyp, iDir, ix, iy, nBndy, nRod, nTot, nPin

REAL :: Ang, Tmp
REAL :: Cnt(2), Eqn(3), Pts(2, 2), cBndyEq(3, 5)

REAL :: mpTypPts(2, 7, 10) = ZERO
REAL :: spTypPts(2, 7,  7) = ZERO

REAL, PARAMETER :: TwoPI_3 = 2._8 * PI_3

INTEGER, POINTER :: pTyp
REAL,    POINTER :: pAng

TYPE(Type_HexAsyTypInfo), POINTER :: aInf_Loc
TYPE(Type_HexGeoTypInfo), POINTER :: gInf_Loc
! ----------------------------------------------------

IF (hLgc%lSngCel) THEN
  CALL HexSetSngCelVtx
  
  RETURN
END IF

aInf_Loc => hAsyTypInfo(iAsyTyp)

nRod = aInf_Loc%nRodPin(1)
nTot = aInf_Loc%nTotPin(1)
nPin = aInf_Loc%nPin

CALL dmalloc(aInf_Loc%PinVtxTyp, nGeoTyp, nTot)
CALL dmalloc(aInf_Loc%PinVtxAng, nGeoTyp, nTot)

cBndyEq = ZERO

cBndyEq(1, 1:3) = 1._8; cBndyEq(2, 2) = Sq3; cBndyEq(2, 3) = -Sq3;
cBndyEq(2, 4:5) = 1._8; cBndyEq(1, 5) = Sq3
! ----------------------------------------------------
!               01. SET : Inn Vtx Typ
! ----------------------------------------------------
DO iPin = 1, nRod
  ix = aInf_Loc%Pin1Dto2Dmap(1, iPin)
  iy = aInf_Loc%Pin1Dto2Dmap(2, iPin)
  
  IF (aInf_Loc%Pin2Dto1Dmap(ix+1, iy)   .EQ. 0) CYCLE
  IF (aInf_Loc%Pin2Dto1Dmap(ix-1, iy)   .EQ. 0) CYCLE
  IF (aInf_Loc%Pin2Dto1Dmap(ix,   iy+1) .EQ. 0) CYCLE
  IF (aInf_Loc%Pin2Dto1Dmap(ix,   iy-1) .EQ. 0) CYCLE
  
  Cnt = aInf_Loc%PinCnt(:, iPin)
  
  DO iGeo = 1, nGeoTyp
    IF (.NOT. aInf_Loc%lGeoPin(iGeo, iPin)) CYCLE
    
    gInf_Loc => hGeoTypInfo(iGeo)
    
    pTyp => aInf_Loc%PinVtxTyp(iGeo, iPin)
    pAng => aInf_Loc%PinVtxAng(iGeo, iPin)
    
    DO iBndy = 1, gInf_Loc%nBndy
      IF (ChkPtEqn(Cnt, gInf_Loc%Eqn(:, iBndy))) EXIT
    END DO
    
    ! 360
    pTyp = 3
    pAng = ZERO
    
    IF (iBndy > gInf_Loc%nBndy) CYCLE
    
    ! 180
    pTyp = 2
    
    IF (ChkPtEqn(Cnt, cBndyEq(:, 4))) pAng =  TwoPI_3
    IF (ChkPtEqn(Cnt, cBndyEq(:, 5))) pAng = -TwoPI_3
  END DO
END DO
! ----------------------------------------------------
!               02. SET : Bndy Vtx Typ
! ----------------------------------------------------
DO iPin = 1, nRod
  IF (aInf_Loc%PinVtxTyp(1, iPin) .NE. 0) CYCLE
  
  Cnt = aInf_Loc%PinCnt(:, iPin)
  Ang = FindPtAng(Cnt)
  
  DO iGeo = 1, nGeoTyp
    IF (.NOT. aInf_Loc%lGeoPin(iGeo, iPin)) CYCLE
    
    pTyp => aInf_Loc%PinVtxTyp(iGeo, iPin)
    pAng => aInf_Loc%PinVtxAng(iGeo, iPin)
    
    DO iBndy = 1, 3
      IF (ChkPtEqn(Cnt, cBndyEq(:, iBndy))) EXIT
    END DO
    
    ! Bndy 01
    aInf_Loc%PinVtxTyp(iGeo, iPin) = 4
    aInf_Loc%PinVtxAng(iGeo, iPin) = Ang - PI_2
    
    IF (iBndy < 4) CYCLE
    
    gInf_Loc => hGeoTypInfo(iGeo)
    
    DO iBndy = 1, gInf_Loc%nBndy
      IF (ChkPtEqn(Cnt, gInf_Loc%Eqn(:, iBndy))) EXIT
    END DO
    
    ! Bndy 02
    Tmp  = Ang - PI_6
    pTyp = 5
    pAng = (int(Tmp/PI_3) + HALF * (SetSgn_REAL(Tmp) - 1))* PI_3
    
    IF (iBndy > gInf_Loc%nBndy) CYCLE
    
    ! Bndy 02 - Spt
    IF (Cnt(1) > ZERO) THEN
      pTyp = 6
      
      IF (Cnt(2) > ZERO) pAng =  ZERO
      IF (Cnt(2) < ZERO) pAng = -TwoPI_3
      
      IF (.NOT. ChkSameVal(ZERO, Cnt(2))) CYCLE
      
      pTyp = 7
      pAng = TwoPI_3
    ELSE
      pTyp = 7
      
      IF (Cnt(2) > ZERO) pAng = -TwoPI_3
      IF (Cnt(2) < ZERO) pAng = ZERO
      
      IF (.NOT. ChkSameVal(ZERO, Cnt(2))) CYCLE
      
      pTyp = 6
      pAng = TwoPI_3
    END IF
  END DO
END DO
! ----------------------------------------------------
!               03. SET : Cnt Vtx Typ
! ----------------------------------------------------
iPin = aInf_Loc%Pin2Dto1Dmap(nPin, nPin)

aInf_Loc%PinVtxTyp(1, iPin) = 3
aInf_Loc%PinVtxAng(1, iPin) = ZERO

IF (nGeoTyp > 1) THEN
  aInf_Loc%PinVtxTyp(2, iPin) = 2; aInf_Loc%PinVtxAng(2, iPin) =  ZERO
  aInf_Loc%PinVtxTyp(3, iPin) = 2; aInf_Loc%PinVtxAng(3, iPin) =  TwoPI_3
  aInf_Loc%PinVtxTyp(4, iPin) = 2; aInf_Loc%PinVtxAng(4, iPin) = -TwoPi_3
  aInf_Loc%PinVtxTyp(5, iPin) = 1; aInf_Loc%PinVtxAng(5, iPin) =  ZERO
  aInf_Loc%PinVtxTyp(6, iPin) = 1; aInf_Loc%PinVtxAng(6, iPin) =  TwoPI_3
  aInf_Loc%PinVtxTyp(7, iPin) = 1; aInf_Loc%PinVtxAng(7, iPin) = -TwoPI_3
END IF
! ----------------------------------------------------
!               04. SET : Gap Vtx Typ
! ----------------------------------------------------
DO iPin = nRod + 1, nTot
  Cnt = aInf_Loc%PinCnt(:, iPin)
  Ang = FindPtAng(Cnt)
  
  pTyp => aInf_Loc%PinVtxTyp(1, iPin) ! NOTICE : iGeo = 1
  pAng => aInf_Loc%PinVtxAng(1, iPin)
  
  Tmp  = Ang - PI_6
  pAng = (int(Tmp/PI_3) + HALF * (SetSgn_REAL(Tmp) - 1))* PI_3
  
  pTyp = 10
  jMod = mod(iPin - nRod, 2*nPin - 2)
  
  IF (jMod .EQ. 0) pTyp = 8
  IF (jMod .EQ. 1) pTyp = 9
END DO

DO iPin = nRod + 1, nTot
  DO iGeo = 2, nGeoTyp
    IF (.NOT. aInf_Loc%lGeoPin(iGeo, iPin)) CYCLE
    
    aInf_Loc%PinVtxTyp(iGeo, iPin) = aInf_Loc%PinVtxTyp(1, iPin)
    aInf_Loc%PinVtxAng(iGeo, iPin) = aInf_Loc%PinVtxAng(1, iPin)
  END DO
END DO
! ----------------------------------------------------
!               05. SET : Pin Vtx
! ----------------------------------------------------
mpTypPts = aInf_Loc%mpTypPts
spTypPts = aInf_Loc%spTypPts
 
CALL dmalloc(aInf_Loc%mpVtx, 2, 7, nGeoTyp, nTot)
CALL dmalloc(aInf_Loc%spVtx, 2, 7, nGeoTyp, nRod)

DO iPin = 1, nTot
  DO iGeo = 1, nGeoTyp
    IF (.NOT. aInf_Loc%lGeoPin(iGeo, iPin)) CYCLE
    
    ivTyp = aInf_Loc%PinVtxTyp(iGeo, iPin)
    Ang   = aInf_Loc%PinVtxAng(iGeo, iPin)
    Cnt   = aInf_Loc%PinCnt    (1:2, iPin)
    
    DO iBndy = 1, mpTypNumNgh(ivTyp) + 1
      aInf_Loc%mpVtx(1:2, iBndy, iGeo, iPin) = RotPt(mpTypPts(1:2, iBndy, ivTyp), Ang) + Cnt(1:2)
    END DO
    
    IF (iPin .GT. nRod) CYCLE
    
    DO iBndy = 1, spTypNumNgh(ivTyp) + 1
      aInf_Loc%spVtx(1:2, iBndy, iGeo, iPin) = RotPt(spTypPts(1:2, iBndy, ivTyp), Ang) + Cnt(1:2)
    END DO
  END DO
END DO

CALL ChkArrayZero(aInf_Loc%mpVtx, 1, 2, 1, 7, 1, nGeoTyp, 1, nTot)
CALL ChkArrayZero(aInf_Loc%spVtx, 1, 2, 1, 7, 1, nGeoTyp, 1, nRod)

!CALL HexTstPinVtx(iAsyTyp)
!CALL HexTstSpVtx (iAsyTyp)

NULLIFY (pTyp, pAng, gInf_Loc, aInf_Loc)
! ----------------------------------------------------

END SUBROUTINE HexSetAsyTypPinVtx
! ------------------------------------------------------------------------------------------------------------
!                                     03. HEX SET : Sng Cel Vtx
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetSngCelVtx()

USE allocs
USE PARAM,   ONLY : ZERO, HALF
USE geom,    ONLY : nAsyType0
USE HexType, ONLY : Type_HexAsyTypInfo
USE HexData, ONLY : hAsyTypInfo, AsyVtx

IMPLICIT NONE

INTEGER :: iaTyp

TYPE(Type_HexAsyTypInfo), POINTER :: aInf_Loc
! ----------------------------------------------------

DO iaTyp = 1, nAsyType0
  IF (hAsyTypInfo(iaTyp)%luse) EXIT
END DO

aInf_Loc => hAsyTypInfo(iaTyp) ! # of used haTyp must be 1 for Sng Cel

CALL dmalloc(aInf_Loc%PinVtxTyp, 1, 1)
CALL dmalloc(aInf_Loc%PinVtxAng, 1, 1)

aInf_Loc%PinVtxTyp(1, 1) = 3
aInf_Loc%PinVtxAng(1, 1) = ZERO

CALL dmalloc(aInf_Loc%mpVtx, 2, 7, 1, 1)

aInf_Loc%mpVtx(1:2, 1:7, 1, 1) = AsyVtx(1:2, 1:7)

NULLIFY (aInf_Loc)
! ----------------------------------------------------

END SUBROUTINE HexSetSngCelVtx
! ------------------------------------------------------------------------------------------------------------

END MODULE HexVtx