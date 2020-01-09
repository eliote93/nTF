! ------------------------------------------------------------------------------------------------------------
!                                     HEX SET : Asy Typ Pin Ngh
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetAsyTypSPngh(iaTyp)

USE PARAM,   ONLY : ZERO
USE HexType, ONLY : Type_HexAsyTypInfo
USE HexData, ONLY : hLgc, hAsyTypInfo, nGeoTyp, hEps, spTypNumNgh
USE HexUtil, ONLY : SetEqn, ChkPtEqn, CalPtLineSgn

IMPLICIT NONE

INTEGER :: iaTyp
INTEGER :: iGeo, iPin, jPin, ivTyp, iBndy, jBndy, iNum, nRod, nTot, nBndy
LOGICAL :: lChk01, lChk02, lChk03, lChk04, lChk05
REAL    :: Cnt(2), Pts(2, 4), Eqn(3), xRng(2), yRng(2), Val

INTEGER :: BndyIdx(1:6, 4:7) = 0

TYPE(Type_HexAsyTypInfo), POINTER :: aInf_Loc
! ----------------------------------------------------

aInf_Loc => hAsyTypInfo(iaTyp)

nRod = aInf_Loc%nRodPin(1)
nTot = aInf_Loc%nTotPin(1)

ALLOCATE (aInf_Loc%cpSlfMPnum       (nGeoTyp, nRod)); aInf_Loc%cpSlfMPnum = 1
ALLOCATE (aInf_Loc%cpSlfMPidx    (3, nGeoTyp, nRod)); aInf_Loc%cpSlfMPidx = 0
ALLOCATE (aInf_Loc%cpSufMPidx (2, 6, nGeoTyp, nRod)); aInf_Loc%cpSufMPidx = 0
ALLOCATE (aInf_Loc%cpSufMPsuf (2, 6, nGeoTyp, nRod)); aInf_Loc%cpSufMPsuf = 0
ALLOCATE (aInf_Loc%cpSufMPnum    (6, nGeoTyp, nRod)); aInf_Loc%cpSufMPnum = 0
! ----------------------------------------------------
!               01. SET : Self Data
! ----------------------------------------------------
BndyIdx(1, 4) = 3; BndyIdx(1, 5) = 1; BndyIdx(1, 6) = 1; BndyIdx(1, 7) = 3
BndyIdx(3, 4) = 3; BndyIdx(4, 5) = 3; BndyIdx(3, 6) = 3; BndyIdx(3, 7) = 1
BndyIdx(4, 4) = 4; BndyIdx(5, 5) = 4; BndyIdx(4, 6) = 4; BndyIdx(4, 7) = 4
BndyIdx(5, 4) = 4

Cnt = ZERO

DO iGeo = 1, nGeoTyp
  DO iPin = 1, nRod
    IF (.NOT. aInf_Loc%lGeoPin(iGeo, iPin)) CYCLE
    
    ivTyp = aInf_Loc%PinVtxTyp(iGeo, iPin)
    nBndy = spTypNumNgh(ivTyp)
    
    aInf_Loc%cpSlfMPidx(1, iGeo, iPin) = iPin
    
    ! Inn Pin
    IF (ivTyp < 4) THEN
      DO iBndy = 1, nBndy
        aInf_Loc%cpSufMPnum(   iBndy, iGeo, iPin) = 1
        aInf_Loc%cpSufMPidx(1, iBndy, iGeo, iPin) = iPin
        aInf_Loc%cpSufMPsuf(1, iBndy, iGeo, iPin) = iBndy
      END DO
      
      CYCLE
    END IF
    
    ! Bndy Pin - Inn Bndy
    DO iBndy = 1, 3 + ivTyp / 5 - ivTyp / 6
      aInf_Loc%cpSufMPnum(   iBndy, iGeo, iPin) = 1
      aInf_Loc%cpSufMPidx(1, iBndy, iGeo, iPin) = iPin
      aInf_Loc%cpSufMPsuf(1, iBndy, iGeo, iPin) = iBndy
    END DO
    
    ! Bndy Pin - FIND : Ngh Gap Pin
    DO iBndy = 1, nBndy
      IF (iBndy .EQ. 2) CYCLE
      
      jBndy = BndyIdx(iBndy, ivTyp)
      
      IF (jBndy .EQ. 0) CYCLE
      
      Pts(1:2, 1) = aInf_Loc%spVtx(1:2, iBndy,   iGeo, iPin)
      Pts(1:2, 2) = aInf_Loc%spVtx(1:2, iBndy+1, iGeo, iPin)
      
      Eqn  = SetEqn(Pts(1:2, 1), Pts(1:2, 2), Cnt)
      iNum = aInf_Loc%cpSufMPnum(iBndy, iGeo, iPin)
      xRng = [min(Pts(1, 1), Pts(1, 2)) - hEps, max(Pts(1, 1), Pts(1, 2)) + hEps]
      yRng = [min(Pts(2, 1), Pts(2, 2)) - hEps, max(Pts(2, 1), Pts(2, 2)) + hEps]
      
      DO jPin = nRod+1, nTot
        IF (.NOT. aInf_Loc%lGeoPin(iGeo, jPin)) CYCLE
        
        Val = CalPtLineSgn(aInf_Loc%PinCnt(:, jPin), Eqn, aInf_Loc%PinCnt(:, iPin))
        
        IF (Val < ZERO) CYCLE
        
        Pts(1:2, 3) = aInf_Loc%mpVtx(1:2, jBndy,   iGeo, jPin)
        Pts(1:2, 4) = aInf_Loc%mpVtx(1:2, jBndy+1, iGeo, jPin)
        
        lChk01 = ChkPtEqn(Pts(1:2, 3), Eqn)
        lChk02 = ChkPtEqn(Pts(1:2, 4), Eqn)
        lChk05 = lChk01 .AND. lChk02
        
        IF (.NOT. lChk05) CYCLE
        
        lChk01 = min(Pts(1, 3), Pts(1, 4)) > xRng(1)
        lChk02 = max(Pts(1, 3), Pts(1, 4)) < xRng(2)
        lChk03 = min(Pts(2, 3), Pts(2, 4)) > yRng(1)
        lChk04 = max(Pts(2, 3), Pts(2, 4)) < yRng(2)
          
        lChk05 = lChk01 .AND. lChk02 .AND. lChk03 .AND. lChk04
        
        IF (.NOT. lChk05) CYCLE
        
        iNum = iNum + 1
        
        aInf_Loc%cpSufMPnum      (iBndy, iGeo, iPin) = iNum
        aInf_Loc%cpSufMPidx(iNum, iBndy, iGeo, iPin) = jPin
        aInf_Loc%cpSufMPsuf(iNum, iBndy, iGeo, iPin) = jBndy
      END DO
    END DO
    
    CALL SetSlfIdx(iGeo, iPin)
  END DO
END DO

NULLIFY (aInf_Loc)

CONTAINS
! ----------------------------------------------------
!               SUB. SET : Slf Idx
! ----------------------------------------------------
SUBROUTINE SetSlfIdx(iGeo, iPin)

USE PARAM, ONLY : FALSE, TRUE

INTEGER :: iGeo, iPin, kBndy, kPin, tPin, xPin
LOGICAL :: lFst
! ----------------------------

lFst = TRUE

DO kBndy = 1, 6
  DO kPin = 1, 2
    tPin = aInf_Loc%cpSufMPidx(kPin, kBndy, iGeo, iPin)
    
    IF ((tPin .EQ. 0).OR.(tPin .EQ. iPin)) CYCLE
    
    IF (lFst) THEN
      aInf_Loc%cpSlfMPnum   (iGeo, iPin) = 2
      aInf_Loc%cpSlfMPidx(2, iGeo, iPin) = tPin
      
      xPin = tPin
      lFst = FALSE
    ELSE IF (tPin .NE. xPin) THEN
      aInf_Loc%cpSlfMPnum   (iGeo, iPin) = 3
      aInf_Loc%cpSlfMPidx(3, iGeo, iPin) = tPin
      
      RETURN
    END IF
  END DO
END DO

END SUBROUTINE SetSlfIdx
! ----------------------------------------------------

END SUBROUTINE HexSetAsyTypSPngh