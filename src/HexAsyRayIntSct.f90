MODULE HexAsyRayIntSct

USE ioutil, ONLY : terminate
  
IMPLICIT NONE

CONTAINS
! ------------------------------------------------------------------------------------------------------------
!                                     01. SET : Ray Int Sct - Pin
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetRayIntSct_Pin(icBss, RayEqn, imRay, iGeo, nhpRay, hpRay)

USE PARAM,   ONLY : HALF, ZERO, TRUE, FALSE
USE HexType, ONLY : Type_HexRayPinInfo, Type_HexPinRay
USE HexUtil, ONLY : SolveLineEqn, ChkEqnTwoVtx, ChkSamePts, ChkPtEqn
USE HexData, ONLY : RayPinInfo, hAsyTypInfo, hCelBss, AsyEqn
USE HexTst,  ONLY : HexTsthpRay

IMPLICIT NONE

REAL    :: RayEqn(3)
INTEGER :: icBss, imRay, iGeo, nhpRay, iTmp

INTEGER, PARAMETER :: nMaxPinRay = 100 ! ARBITRARY

TYPE(Type_HexPinRay) :: hpRay(:)
TYPE(Type_HexPinRay) :: tpRay(nMaxPinRay)

INTEGER :: iBndy, iPin
INTEGER :: nTotPin, nPtTmp
REAL    :: Sol(2), EdPts(2, 2)
LOGICAL :: lSol, lChk

TYPE(Type_HexRayPinInfo), POINTER :: pInf_Loc
! ----------------------------------------------------

nTotPin = hAsyTypInfo(hCelBss(icBss)%iaTyp)%nTotPin(1)
nhpRay  = 1

DO iPin = 1, nTotPin
  pInf_Loc => RayPinInfo(iPin)
  
  nPtTmp = 0
  
  DO iBndy = 1, pInf_Loc%nBndy(iGeo)
    lSol = ChkEqnTwoVtx(RayEqn, pInf_Loc%Vtx(1:2, iBndy, iGeo), pInf_Loc%Vtx(1:2, iBndy + 1, iGeo))
    lSol = lSol .OR. ChkPtEqn(pInf_Loc%Vtx(1:2, iBndy,   iGeo), RayEqn)
    lSol = lSol .OR. ChkPtEqn(pInf_Loc%Vtx(1:2, iBndy+1, iGeo), RayEqn)
    
    IF (.NOT. lSol) CYCLE
    
    CALL SolveLineEqn(RayEqn, pInf_Loc%Eqn(1:3, iBndy, iGeo), Sol, lSol)
    
    DO iTmp = 1, nPtTmp
      IF (ChkSamePts(Sol, tpRay(nhpRay)%PinPt(1:2, iTmp))) EXIT
    END DO
    
    IF (iTmp .LE. nPtTmp) CYCLE
    
    nPtTmp = nPtTmp + 1
    
    tpRay(nhpRay)%SufIdx    (nPtTmp) = iBndy
    tpRay(nhpRay)%PinPt(1:2, nPtTmp) = Sol(1:2)
  END DO
  
  IF (nPtTmp .NE. 2) CYCLE
  
  lChk = ChkSamePts(tpRay(nhpRay)%PinPt(1:2, 1), tpRay(nhpRay)%PinPt(1:2, 2))
  
  IF (lChk) CYCLE
  
  tpRay(nhpRay)%PinIdx = iPin
  
  nhpRay = nhpRay + 1
  
  IF (nhpRay > nMaxPinRay) CALL terminate("SET : RAY INT SCT - PIN")
END DO

nhpRay = nhpRay - 1

IF (nhpRay < 1) RETURN

CALL HexSortPinRay(tpRay, hpRay, nhpRay)

EdPts(1:2, 1) = hpRay(1)%PinPt(1:2, 1)
EdPts(1:2, 2) = hpRay(nhpRay)%PinPt(1:2, 2)

DO iBndy = 1, 6
  IF (ChkPtEqn(EdPts(1:2, 1), AsyEqn(1:3, iBndy))) hpRay(1)%SufIdx(1) = 4
  IF (ChkPtEqn(EdPts(1:2, 2), AsyEqn(1:3, iBndy))) hpRay(nhpRay)%SufIdx(2) = 4
END DO

!IF ((iGeo .EQ. 7).AND.(imRay .EQ. 1369)) CALL HexTsthpRay(RayEqn, nhpRay, hpRay)
! ----------------------------------------------------

END SUBROUTINE HexSetRayIntSct_Pin
! ------------------------------------------------------------------------------------------------------------
!                                     02. SORT : Pin Ray
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSortPinRay(tpRay, hpRay, nhpRay)

USE PARAM,   ONLY : TRUE, FALSE, BIG
USE HexType, ONLY : Type_HexPinRay

IMPLICIT NONE

TYPE(Type_HexPinRay) :: tpRay(100)
TYPE(Type_HexPinRay) :: hpRay(:)
INTEGER :: nhpRay

INTEGER :: iPin, iSuf, jPin, iMn
REAL    :: Pt(2), Min
LOGICAL :: luse(nhpRay)
! ----------------------------------------------------

DO iPin = 1, nhpRay
  IF (tpRay(iPin)%PinPt(2, 1) < tpRay(iPin)%PinPt(2, 2)) CYCLE
  
  Pt(1:2) = tpRay(iPin)%PinPt(1:2, 1)
  
  tpRay(iPin)%PinPt(1:2, 1) = tpRay(iPin)%PinPt(1:2, 2)
  tpRay(iPin)%PinPt(1:2, 2) = Pt(1:2)
  
  iSuf = tpRay(iPin)%SufIdx(1)
  
  tpRay(iPin)%SufIdx(1) = tpRay(iPin)%SufIdx(2)
  tpRay(iPin)%SufIdx(2) = iSuf
END DO

luse = FALSE

DO iPin = 1, nhpRay
  Min = BIG
  
  DO jPin = 1, nhpRay
    IF (luse(jPin)) CYCLE
    IF (tpRay(jPin)%PinPt(2, 1) > Min) CYCLE
    
    Min = tpRay(jPin)%PinPt(2, 1)
    iMn = jPin
  END DO
  
  hpRay(iPin) = tpRay(iMn)
  luse (iMn)  = TRUE
END DO

END SUBROUTINE HexSortPinRay
! ------------------------------------------------------------------------------------------------------------
!                                     03. SET : Int Sct - Msh
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetRayIntSct_Msh(CelRay, RayEqn, iGeo, iPin, BndyPt)

USE PARAM,   ONLY : TRUE, HALF, ZERO
USE HexUtil, ONLY : RotPt, RotLineEqn, Array2DSORT, ChkPtVtx, SolveLineREQN, FindPtLgh
USE HexData, ONLY : RayPinInfo, RayCel, hEps
USE HexType, ONLY : Type_HexRayCel, Type_HexCelRay, Type_HexRayPinInfo

IMPLICIT NONE

REAL    :: RayEqn(3), BndyPt(2, 2)
INTEGER :: iGeo, iPin
TYPE(Type_HexCelRay) :: CelRay

INTEGER, PARAMETER :: nMaxPt = 200 ! ARBITRARY

INTEGER :: ivTyp, iEqn, iPt, iCor, iMsh, jMsh, pMsh, iRev
INTEGER :: nSolPt, nPt0, nPt
REAL    :: Ang, Cnt(2), Sol(2, 2), Eqn(5), CelRayEqn(3), Pt(2), CelBndyPt(2, 2), PtTmp(2, nMaxPt)
REAL    :: EndPts(2, 2)
LOGICAL :: lChk

TYPE(Type_HexRayCel),     POINTER :: rCel_Loc
TYPE(Type_HexRayPinInfo), POINTER :: rPin_Loc
! ----------------------------------------------------

rPin_Loc => RayPinInfo(iPin)
rCel_Loc => RayCel(rPin_Loc%iCel)
! ----------------------------------------------------
!               01. SET : Cel Ray Eqn
! ----------------------------------------------------
Cnt   = rPin_Loc%Cnt
Ang   = rPin_Loc%VtxAng(iGeo)
ivTyp = rPin_Loc%VtxTyp(iGeo)

Eqn(1:2) = RayEqn(1:2)
Eqn(3)   = RayEqn(3) - RayEqn(1) * Cnt(1) - RayEqn(2) * Cnt(2)

CelRayEqn(1:3)    = RotLineEqn(Eqn(1:3), -Ang)
CelBndyPt(1:2, 1) = RotPt(BndyPt(1:2, 1) - Cnt(1:2), -Ang)
CelBndyPt(1:2, 2) = RotPt(BndyPt(1:2, 2) - Cnt(1:2), -Ang)
! ----------------------------------------------------
!               02. CALC : Msh Int-Sct
! ----------------------------------------------------
nPt0  = 2
PtTmp = ZERO

PtTmp(1:2, 1) = CelBndyPt(1:2, 1)
PtTmp(1:2, 2) = CelBndyPt(1:2, 2)

DO iEqn = 1, rCel_Loc%nEqn
  Eqn(1:5) = rCel_Loc%Eqn(1:5, iEqn)
  iCor     = rCel_Loc%iEqnCor (iEqn)
  
  CALL SolveLineREqn(CelRayEqn, Eqn, Sol, nSolPt)
  
  IF(nSolPt .EQ. 0) CYCLE
  IF(Eqn(2) > hEps .AND. nSolPt .EQ. 1) CYCLE
  
  DO iPt = 1, nSolPt
    lChk = ChkPtVtx(Sol(1:2, iPt), CelBndyPt(1:2, 1), CelBndyPt(1:2, 2), iCor)
    
    IF(.NOT. lChk) CYCLE
    
    nPt0 = nPt0 + 1
    
    IF (nPt0 > nMaxPt) CALL terminate("SET : RAY INT SCT - MAX PT")
    
    PtTmp(1:2, nPt0) = Sol(1:2, iPt)
  END DO
END DO

CALL Array2DSORT(PtTmp(1, 1:nPt0), PtTmp(2, 1:nPt0), nPt0, nPt, TRUE, TRUE, 2)

CelRay%nSegRay = nPt - 1

ALLOCATE (CelRay%SegLgh (CelRay%nSegRay))
ALLOCATE (CelRay%MshIdx (CelRay%nSegRay))
! ----------------------------------------------------
!               03. FIND : Msh Idx
! ----------------------------------------------------
iRev = 0
pMsh = 0

EndPts(1:2, 1) = RotPt(PtTmp(1:2, 1),   Ang)
EndPts(1:2, 2) = RotPt(PtTmp(1:2, nPt), Ang)

IF (EndPts(2, 1) > EndPts(2, 2)) iRev = 1

DO iPt = 1, CelRay%nSegRay
  Pt(1:2) = HALF * (PtTmp(1:2, iPt) + PtTmp(1:2, iPt+1))
  
  iMsh = FindMshIdx(Pt(1:2))
  jMsh = rCel_Loc%MshIdx(ivTyp, iMsh)
  
  IF (iMsh .EQ. 0)    CALL terminate("RAY MSH INT SCT - NO MSH")
  IF (iMsh .EQ. pMsh) CALL terminate("RAY MSH INT SCT - SAME MSH")
  
  CelRay%MshIdx(iPt + iRev*(nPt - 2*iPt)) = jMsh
  CelRay%SegLgh(iPt + iRev*(nPt - 2*iPt)) = FindPtLgh(PtTmp(1:2, iPt), PtTmp(1:2, iPt+1)) * 1000._8 ! NOTICE
  
  pMsh = iMsh
END DO

! For DEBUG
!ALLOCATE (CelRay%SegPts (2, nPt))
!
!DO iPt = 1, nPt
!  CelRay%SegPts(1:2, iPt + iRev*(nPt + 1 - 2*iPt)) = RotPt(PtTmp(1:2, iPt), Ang) + rPin_Loc%Cnt(1:2)
!END DO

NULLIFY (rPin_Loc, rCel_Loc)

CONTAINS
! ----------------------------------------------------
!               SUB. FIND : Msh Idx
! ----------------------------------------------------
FUNCTION FindMshIdx(Pt)

IMPLICIT NONE

INTEGER :: FindMshIdx
INTEGER :: iMsh, iEqn, jEqn
REAL    :: Pt(2), Eqn(5), Val, Tmp
! ----------------------------

FindMshIdx = 0

DO iMsh = 1, rCel_Loc%nMsh
  DO iEqn = 1, rCel_Loc%nMshEqn(iMsh)
    jEqn = rCel_Loc%MshEqnLst(iEqn, iMsh)
    Eqn  = rCel_Loc%Eqn(:, jEqn)
    
    Val = Eqn(5) - Eqn(1) * Pt(1) - Eqn(2) * Pt(1) * Pt(1) &
                 - Eqn(3) * Pt(2) - Eqn(4) * Pt(2) * Pt(2)
    Tmp = Val * rCel_Loc%MshEqnVal(iEqn, iMsh)
    
    IF (Tmp < 0) EXIT
  END DO
  
  IF(Tmp < 0) CYCLE
  
  FindMshIdx = iMsh
  
  RETURN
END DO

END FUNCTION FindMshIdx
! ----------------------------------------------------

END SUBROUTINE HexSetRayIntSct_Msh
! ------------------------------------------------------------------------------------------------------------

END MODULE HexAsyRayIntSct