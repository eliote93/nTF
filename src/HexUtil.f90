MODULE HexUtil

USE PARAM,   ONLY : ZERO, HALF, TRUE, FALSE, PI
USE HexData, ONLY : hEps, PI_2, PI_6
USE ioutil,  ONLY : terminate

IMPLICIT NONE

INTERFACE ChkArrayZERO
    MODULE PROCEDURE ChkArrayZERO_1D
    MODULE PROCEDURE ChkArrayZERO_2D
    MODULE PROCEDURE ChkArrayZERO_3D
    MODULE PROCEDURE ChkArrayZERO_4D
END INTERFACE

CONTAINS

! ------------------------------------------------------------------------------------------------------------
!                                     01. SET : Pts Priority
! ------------------------------------------------------------------------------------------------------------
FUNCTION SetPtsPri(Pts)

IMPLICIT NONE

REAL    :: Pts(2, 2)
INTEGER :: SetPtsPri
! ----------------------------------------------------

SetPtsPri = 2

IF (ChkSameVal(Pts(2, 1), Pts(2, 2))) SetPtsPri = 1

END FUNCTION SetPtsPri
! ------------------------------------------------------------------------------------------------------------
!                                     01. SET : Sign - REAL
! ------------------------------------------------------------------------------------------------------------
FUNCTION SetSgn_REAL(Val)
! ASSUME : Val is not ZERO

IMPLICIT NONE

REAL :: Val, SetSgn_REAL
! ----------------------------------------------------

SetSgn_REAL = Val / abs(Val)

END FUNCTION SetSgn_REAL
! ------------------------------------------------------------------------------------------------------------
!                                     02. NORM : Vector
! ------------------------------------------------------------------------------------------------------------
FUNCTION NormVec(Pt)

IMPLICIT NONE

REAL :: Pt(2), NormVec(2), Lgh
! ----------------------------------------------------

Lgh = Pt(1)*Pt(1) + Pt(2)*Pt(2)
Lgh = sqrt(Lgh)

NormVec(1) = Pt(1) / Lgh
NormVec(2) = Pt(2) / Lgh
! ----------------------------------------------------

END FUNCTION NormVec
! ------------------------------------------------------------------------------------------------------------
!                                     17. CHK : Pt Eqn
! ------------------------------------------------------------------------------------------------------------
FUNCTION SetSgn_INT(Val)
! ASSUME : Val is not ZERO

IMPLICIT NONE

INTEGER :: Val, SetSgn_INT
! ----------------------------------------------------

SetSgn_INT = Val / abs(Val)

END FUNCTION SetSgn_INT
! ------------------------------------------------------------------------------------------------------------
!                                     03. CHK : 1d Array Zero
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE ChkArrayZero_1D(Val, n1, n2)

IMPLICIT NONE

INTEGER :: n1, n2, ii

REAL :: Val(n1:n2)
! ----------------------------------------------------

DO ii = n1, n2
  IF (abs(Val(ii)) < hEps) Val(ii) = ZERO
END DO
! ----------------------------------------------------

END SUBROUTINE ChkArrayZero_1D
! ------------------------------------------------------------------------------------------------------------
!                                     04. CHK : 2d Array Zero
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE ChkArrayZero_2D(Val, n1, n2, n3, n4)

IMPLICIT NONE

INTEGER :: n1, n2, n3, n4, ii, jj

REAL :: Val(n1:n2, n3:n4)
! ----------------------------------------------------

DO jj = n3, n4
  DO ii = n1, n2
    IF (abs(Val(ii, jj)) < hEps) Val(ii, jj) = ZERO
  END DO
END DO
! ----------------------------------------------------

END SUBROUTINE ChkArrayZero_2D
! ------------------------------------------------------------------------------------------------------------
!                                     05. CHK : 3d Array Zero
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE ChkArrayZero_3D(Val, n1, n2, n3, n4, n5, n6)

IMPLICIT NONE

INTEGER :: n1, n2, n3, n4, n5, n6, ii, jj, kk

REAL :: Val(n1:n2, n3:n4, n5:n6)
! ----------------------------------------------------

DO kk = n5, n6
  DO jj = n3, n4
    DO ii = n1, n2
      IF (abs(Val(ii, jj, kk)) < hEps) Val(ii, jj, kk) = ZERO
    END DO
  END DO
END DO
! ----------------------------------------------------

END SUBROUTINE ChkArrayZero_3D
! ------------------------------------------------------------------------------------------------------------
!                                     06. CHK : 4d Array Zero
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE ChkArrayZero_4D(Val, n1, n2, n3, n4, n5, n6, n7, n8)

IMPLICIT NONE

INTEGER :: n1, n2, n3, n4, n5, n6, n7, n8, ii, jj, kk, tt

REAL :: Val(n1:n2, n3:n4, n5:n6, n7:n8)
! ----------------------------------------------------

DO tt = n7, n8
  DO kk = n5, n6
    DO jj = n3, n4
      DO ii = n1, n2
        IF (abs(Val(ii, jj, kk, tt)) < hEps) Val(ii, jj, kk, tt) = ZERO
      END DO
    END DO
  END DO
END DO
! ----------------------------------------------------

END SUBROUTINE ChkArrayZero_4D
! ------------------------------------------------------------------------------------------------------------
!                                     07. SET : New Pt
! ------------------------------------------------------------------------------------------------------------
FUNCTION SetNewPt(Lgh, Ang)

IMPLICIT NONE

REAL :: SetNewPt(2)
REAL :: Lgh, Ang
! ----------------------------------------------------

SetNewPt(1) = Lgh * cos(Ang)
SetNewPt(2) = Lgh * sin(Ang)

CALL ChkArrayZero(SetNewPt, 1, 2)
! ----------------------------------------------------

END FUNCTION SetNewPt
! ------------------------------------------------------------------------------------------------------------
!                                     08. FIND : Line Pt Length
! ------------------------------------------------------------------------------------------------------------
FUNCTION FindLinePtLgh(Line, Pt)

IMPLICIT NONE

REAL :: FindLinePtLgh
REAL :: Line(3), Pt(2)
! ----------------------------------------------------

FindLinePtLgh = Line(1)*Pt(1) + Line(2)*Pt(2) - Line(3)
FindLinePtLgh = abs(FindLinePtLgh)
FindLinePtLgh = FindLinePtLgh / sqrt(Line(1)*Line(1) + Line(2)*Line(2))
! ----------------------------------------------------

END FUNCTION FindLinePtLgh
! ------------------------------------------------------------------------------------------------------------
!                                     09. FIND : Ref Ang
! ------------------------------------------------------------------------------------------------------------
FUNCTION FindRefAng(Ang01, Ang02)
! Ang01 is reflected on Ang02

IMPLICIT NONE

REAL :: FindRefAng(2)
REAL :: Ang01(2), Ang02(2) ! (cosv, sinv)

REAL :: Lgh, Nrm01(2), Nrm02(2)
! ----------------------------------------------------

Nrm01 = NormVec(Ang01)
Nrm02 = NormVec(Ang02)

Lgh = Nrm01(1)*Nrm02(1) + Nrm01(2)*Nrm02(2)

FindRefAng(1) = 2._8 * Lgh * Nrm02(1) - Nrm01(1) ! cosv
FindRefAng(2) = 2._8 * Lgh * Nrm02(2) - Nrm01(2) ! sinv

CALL ChkArrayZero(FindRefAng, 1, 2)
! ----------------------------------------------------

END FUNCTION FindRefAng
! ------------------------------------------------------------------------------------------------------------
!                                     10. FIND : Rot Ang
! ------------------------------------------------------------------------------------------------------------
FUNCTION FindRotAng(Ang01, Ang02, rAng)
! Ang01 is rotated on Ang02 by rAng

IMPLICIT NONE

REAL :: FindRotAng(2)
REAL :: rAng, Ang01(2), Ang02(2) ! (cosv, sinv)

REAL :: rTmp, Nrm01(2), Nrm02(2)
! ----------------------------------------------------

Nrm01 = NormVec(Ang01)
Nrm02 = NormVec(Ang02)

IF (ChkSameVal(Nrm02(2), ZERO)) THEN
  rTmp = -rAng
ELSE IF (ChkSameVal(Nrm02(1), cos(rAng))) THEN
  rTmp = rAng
ELSE
  CALL terminate("WRONG INPUT FOR FIND ROT ANG")
END IF

FindRotAng(1) = Nrm01(1)*cos(rTmp) - Nrm01(2)*sin(rTmp) ! cosv
FindRotAng(2) = Nrm01(2)*cos(rTmp) + Nrm01(1)*sin(rTmp) ! sinv

CALL ChkArrayZero(FindRotAng, 1, 2)
! ----------------------------------------------------

END FUNCTION FindRotAng
! ------------------------------------------------------------------------------------------------------------
!                                     11. FIND : Rot Pt
! ------------------------------------------------------------------------------------------------------------
FUNCTION FindRotPt(Pt, Ang)
! Pt is rotated by rAng

IMPLICIT NONE

REAL :: FindRotPt(2)
REAL :: Ang, Pt(2) ! (cosv, sinv)

REAL :: Theta, Nrm(2)
! ----------------------------------------------------

Nrm = NormVec(Pt)

IF (ChkSameVal(Nrm(2), ZERO)) THEN
  Theta = -Ang
ELSE IF (ChkSameVal(Nrm(2)/Nrm(1), -tan(Ang))) THEN
  Theta = Ang
ELSE
  CALL terminate("WRONG INPUT FOR FIND ROT PT")
END IF

FindRotPt = RotPt(Pt, Theta)
! ----------------------------------------------------

END FUNCTION FindRotPt
! ------------------------------------------------------------------------------------------------------------
!                                     12. CAL : Pt Line Sign
! ------------------------------------------------------------------------------------------------------------
FUNCTION CalPtLineSgn(Pt, Eqn, Cnt)
! Cnt : Positive
! Result = 0 : Cnt or Pt belong to line

IMPLICIT NONE

REAL :: CalPtLineSgn
REAL :: Pt(2), Eqn(3), Cnt(2)

REAL :: Tmp01, Tmp02, Tmp03
! ----------------------------------------------------

Tmp01 = Eqn(3) -  Pt(1)*Eqn(1) -  Pt(2)*Eqn(2)
Tmp02 = Eqn(3) - Cnt(1)*Eqn(1) - Cnt(2)*Eqn(2)
Tmp03 = Tmp01 * Tmp02

IF (abs(Tmp03) < hEps) THEN
  CalPtLineSgn = ZERO
ELSE
  CalPtLineSgn = SetSgn_REAL(Tmp03)
END IF
! ----------------------------------------------------

END FUNCTION CalPtLineSgn
! ------------------------------------------------------------------------------------------------------------
!                                     13. SET : Eqn
! ------------------------------------------------------------------------------------------------------------
FUNCTION SetEqn(Pt01, Pt02, Cnt)
! Eqn(1) * x + Eqn(2) * y = Eqn(3)
! Always Eqn(3) - Eqn(1) * x - Eqn(2) * y >= 0

IMPLICIT NONE

REAL :: SetEqn(3)

INTEGER :: iCor

REAL :: Pt01(2), Pt02(2), Cnt(2)
REAL :: Lgh, sinv, cosv, Val
! ----------------------------------------------------

! ----------------------------------------------------
!               01. SET : Eqn
! ----------------------------------------------------
Lgh = (Pt02(1) - Pt01(1))**2 + &
      (Pt02(2) - Pt01(2))**2
Lgh = sqrt(Lgh)

IF (Lgh < hEps) CALL terminate("SAME PTS")

sinv = (Pt02(2) - Pt01(2)) / Lgh
cosv = (Pt02(1) - Pt01(1)) / Lgh

SetEqn(1) =  sinv
SetEqn(2) = -cosv
SetEqn(3) = sinv * Pt01(1) - cosv * Pt01(2)
! ----------------------------------------------------
!               02. PROCESS
! ----------------------------------------------------
CALL ChkArrayZero(SetEqn, 1, 3)

Val = SetEqn(3) - SetEqn(1) * Cnt(1) &
                - SetEqn(2) * Cnt(2)

IF(Val < 0) THEN
  SetEqn(1:3) = - SetEqn(1:3)
END IF
! ----------------------------------------------------

END FUNCTION SetEqn
! ------------------------------------------------------------------------------------------------------------
!                                     14. SET : Rot Change
! ------------------------------------------------------------------------------------------------------------
FUNCTION RotPt(Pt, Theta)
! Rotate Pt with theta (clockwise)

IMPLICIT NONE

REAL :: Pt(2)
REAL :: Theta
REAL :: RotPt(2)
! ----------------------------------------------------

RotPt(1) = Pt(1) * cos(Theta) - Pt(2) * sin(Theta)
RotPt(2) = Pt(1) * sin(Theta) + Pt(2) * cos(Theta)

CALL ChkArrayZero(RotPt, 1, 2)
! ----------------------------------------------------

END FUNCTION RotPt
! ------------------------------------------------------------------------------------------------------------
!                                     15. ROT : Line Eqn
! ------------------------------------------------------------------------------------------------------------
FUNCTION RotLineEqn(Eqn, Theta)
! Rotate Eqn with theta (clockwise)
! ASSUME : Center equals to (0, 0)

IMPLICIT NONE

REAL :: Eqn(3), Pts(2, 2), Cnt(2)
REAL :: Theta
REAL :: RotLineEqn(3)

INTEGER :: iCor, jCor, iPt
! ----------------------------------------------------

IF (abs(Eqn(2)) > hEps) THEN
  iCor = 1; jCor = 2 ! x = -1, 1
ELSE
  iCor = 2; jCor = 1 ! y = -1, 1
END IF

DO iPt = 1, 2
  Pts(iCor, iPt) = real(-1 + 2 * mod(iPt, 2))
  Pts(jCor, iPt) = (Eqn(3) - Pts(iCor, iPt) * Eqn(iCor)) / Eqn(jCor)

  Pts(1:2, iPt) = RotPt(Pts(1:2, iPt), Theta)
END DO

Cnt = ZERO

RotLineEqn = SetEqn(Pts(1:2, 1), Pts(1:2, 2), Cnt)
! ----------------------------------------------------

END FUNCTION RotLineEqn
! ------------------------------------------------------------------------------------------------------------
!                                     16. SOLVE : Line Eqn
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SolveLineEqn(Eqn01, Eqn02, Sol, lSol)

IMPLICIT NONE

REAL    :: Eqn01(3), Eqn02(3)
REAL    :: Sol(2)
LOGICAL :: lSol

REAL :: Tmp
! ----------------------------------------------------

IF (abs(Eqn01(1)) > hEps) THEN
  Tmp = Eqn02(2) - Eqn02(1) * Eqn01(2) / Eqn01(1)

  IF (abs(Tmp) > hEps) THEN
    Sol(2) = (Eqn02(3) - Eqn02(1) * Eqn01(3) / Eqn01(1)) / Tmp
    Sol(1) = (Eqn01(3) - Eqn01(2) * Sol(2)) / Eqn01(1)
  ELSE
    !CALL terminate("PARALLEL EQN")
    lSol = .FALSE.; RETURN
  END IF
ELSE
  IF (abs(Eqn02(1)) > hEps) THEN
    Sol(2) = Eqn01(3) / Eqn01(2)
    Sol(1) = (Eqn02(3) - Eqn02(2) * Eqn01(3) / Eqn01(2)) / Eqn02(1)
  ELSE
    !CALL terminate("PARALLEL EQN")
    lSol = .FALSE.; RETURN
  END IF
END IF

CALL ChkArrayZero(Sol, 1, 2)

lSol = TRUE
! ----------------------------------------------------

END SUBROUTINE SolveLineEqn
! ------------------------------------------------------------------------------------------------------------
!                                     17. FIND : Area
! ------------------------------------------------------------------------------------------------------------
FUNCTION FindArea(nPt, Pts)

IMPLICIT NONE

INTEGER :: nPt
REAL :: Pts(2, nPt)
REAL :: FindArea

INTEGER :: iPt

REAL, POINTER :: PtsTmp(:, :)
! ----------------------------------------------------

FindArea = ZERO

ALLOCATE (PtsTmp (2, nPt+1))

PtsTmp(1:2, 1:nPt)   = Pts   (1:2, 1:nPt)
PtsTmp(1:2, nPt + 1) = PtsTmp(1:2, 1)
! ----------------------------------------------------
!               01. SET : Area
! ----------------------------------------------------
DO iPt = 1, nPt
  FindArea = FindArea + PtsTmp(1, iPt) * PtsTmp(2, iPt + 1) &
                      - PtsTmp(2, iPt) * PtsTmp(1, iPt + 1)
END DO

FindArea = abs(FindArea) * HALF

NULLIFY (PtsTmp)
! ----------------------------------------------------

END FUNCTION FindArea
! ------------------------------------------------------------------------------------------------------------
!                                     18. FIND : Cnt
! ------------------------------------------------------------------------------------------------------------
FUNCTION FindCnt(nPt, Pts)

IMPLICIT NONE

INTEGER :: nPt
REAL :: Pts(2, nPt)
REAL :: FindCnt(2)

INTEGER :: iPt, iCor
REAL    :: Tmp, Are, Tmp01, Tmp02

REAL, POINTER :: PtsTmp(:, :)
! ----------------------------------------------------

FindCnt = ZERO

ALLOCATE (PtsTmp (2, nPt+1))

PtsTmp(1:2, 1:nPt)   = Pts   (1:2, 1:nPt)
PtsTmp(1:2, nPt + 1) = PtsTmp(1:2, 1)
! ----------------------------------------------------
!               01. SET : Area
! ----------------------------------------------------
Are = ZERO

DO iPt = 1, nPt
  Are = Are + PtsTmp(1, iPt) * PtsTmp(2, iPt + 1) &
            - PtsTmp(2, iPt) * PtsTmp(1, iPt + 1)
END DO

Are = Are * HALF
! ----------------------------------------------------
!               02. CALC : Cnt
! ----------------------------------------------------
DO iCor = 1, 2
  Tmp = ZERO

  DO iPt = 1, nPt
    Tmp01 = PtsTmp(iCor, iPt) + PtsTmp(iCor, iPt + 1)
    Tmp02 = PtsTmp(1, iPt) * PtsTmp(2, iPt + 1) &
          - PtsTmp(2, iPt) * PtsTmp(1, iPt + 1)
    Tmp   = Tmp + Tmp01 * Tmp02
  END DO

  FindCnt(iCor) = Tmp / Are / 6._8
END DO

CALL ChkArrayZero(FindCnt, 1, 2)

NULLIFY (PtsTmp)
! ----------------------------------------------------

END FUNCTION FindCnt
! ------------------------------------------------------------------------------------------------------------
!                                     19. FIND : Pt Lgh
! ------------------------------------------------------------------------------------------------------------
FUNCTION FindPtLgh(Pt01, Pt02)

IMPLICIT NONE

REAL :: Pt01(2), Pt02(2)
REAL :: FindPtLgh

REAL :: Tmp01, Tmp02
! ----------------------------------------------------

Tmp01 = Pt01(1) - Pt02(1)
Tmp02 = Pt01(2) - Pt02(2)

FindPtLgh = sqrt(Tmp01*Tmp01 + Tmp02*Tmp02)
! ----------------------------------------------------

END FUNCTION FindPtLgh
! ------------------------------------------------------------------------------------------------------------
!                                     20. SOLVE : Line Circle Zero
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SolveLineCircle_Zero(Line, RadSq, Sol, nSol)
! ASSUME : Cnt of Circle = (0, 0)
! ASSUME : .NOT. ((a .EQ. 0) .OR. (b .EQ. 0))

IMPLICIT NONE

REAL    :: Line(3), RadSq, LghSq
REAL    :: Sol(2, 2)
INTEGER :: nSol

REAL    :: a, b, c, R, k, Rev
LOGICAL :: lChk
! ----------------------------------------------------

nSol = 0
Sol  = ZERO

a = Line(1)
b = Line(2)
c = Line(3)
R = RadSq
k = a*a + b*b

LghSq = c*c / k

IF (LghSq > RadSq) RETURN

Rev  = SetSgn_REAL(a*b)
nSol = 2

Sol(1, 1) = (a*c + sqrt(R*b*b*k-b*b*c*c))/k
Sol(1, 2) = (a*c - sqrt(R*b*b*k-b*b*c*c))/k
Sol(2, 1) = (b*c - sqrt(R*a*a*k-a*a*c*c) * Rev)/k
Sol(2, 2) = (b*c + sqrt(R*a*a*k-a*a*c*c) * Rev)/k

IF (R*b*b*k-b*b*c*c > hEps) RETURN

nSol = 1
! ----------------------------------------------------

END SUBROUTINE SolveLineCircle_Zero
! ------------------------------------------------------------------------------------------------------------
!                                     21. SOLVE : Line REqn
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SolveLineREqn(Line, REqn, Sol, nSol)
! ASSUME : REqn(2) is always positive

IMPLICIT NONE

REAL, INTENT(IN) :: Line(3)
REAL, INTENT(IN) :: REqn(5)

REAL,    INTENT(OUT) :: Sol(2, 2)
INTEGER, INTENT(OUT) :: nSol

REAL :: LghSq
REAL :: LineNew(3), Pt(2)

LOGICAL :: lSol
! ----------------------------------------------------

nSol = 0
Sol  = ZERO

IF (REqn(2) > 0) THEN
  CALL SolveLineCircle_Zero(Line, REqn(5), Sol, nSol)
ELSE
  LineNew(1) = REqn(1)
  LineNew(2) = REqn(3)
  LineNew(3) = REqn(5)

  CALL SolveLineEqn(Line, LineNew, Pt, lSol)

  IF (.NOT. lSol) RETURN

  nSol = 1
  Sol(1:2, 1) = Pt(1:2)
END IF
! ----------------------------------------------------

END SUBROUTINE SolveLineREqn
! ------------------------------------------------------------------------------------------------------------
!                                     22. SORT : Pts
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SortPts(nPt, Pt, iPri)

IMPLICIT NONE

INTEGER :: nPt, iPri
REAL    :: Pt(2, nPt)

INTEGER :: iPt, iM(1)
REAL    :: Tmp(nPt), PtTmp(2), TT
! ----------------------------------------------------

IF (iPri .EQ. 1) THEN
  Tmp(1:nPt) = Pt(1, 1:nPt)
ELSE
  Tmp(1:nPt) = Pt(2, 1:nPt)
END IF

DO iPt = 1, nPt - 1
  iM = minloc(Tmp(iPt:nPt)) + iPt - 1

  PtTmp(1:2)     = Pt(1:2, iPt)
  Pt(1:2, iPt)   = Pt(1:2, iM(1))
  Pt(1:2, iM(1)) = PtTmp(1:2)

  TT         = Tmp(iPt)
  Tmp(iPt)   = Tmp(iM(1))
  Tmp(iM(1)) = TT
END DO
! ----------------------------------------------------

END SUBROUTINE SortPts
! ------------------------------------------------------------------------------------------------------------
!                                     23. CHK : Pt Eqn
! ------------------------------------------------------------------------------------------------------------
FUNCTION ChkPtEqn(Pt, Eqn)

IMPLICIT NONE

REAL    :: Pt(2)
REAL    :: Eqn(3)
LOGICAL :: ChkPtEqn

REAL :: Tmp
! ----------------------------------------------------

ChkPtEqn = FALSE

Tmp = Eqn(3) - Eqn(1) * Pt(1) &
             - Eqn(2) * Pt(2)

IF (abs(Tmp) < hEps) ChkPtEqn = TRUE
! ----------------------------------------------------

END FUNCTION ChkPtEqn
! ------------------------------------------------------------------------------------------------------------
!                                     24. CHK : Same Pts
! ------------------------------------------------------------------------------------------------------------
FUNCTION ChkSameVal(Val01, Val02)

IMPLICIT NONE

REAL :: Val01, Val02, Tmp

LOGICAL :: ChkSameVal
! ----------------------------------------------------

ChkSameVal = FALSE

Tmp = Val01 - Val02

IF (abs(Tmp) < hEps) ChkSameVal = TRUE
! ----------------------------------------------------

END FUNCTION ChkSameVal
! ------------------------------------------------------------------------------------------------------------
!                                     25. CHK : Same Pts
! ------------------------------------------------------------------------------------------------------------
FUNCTION ChkSamePts(Pt01, Pt02)

IMPLICIT NONE

REAL :: Pt01(2), Pt02(2), Tmp01, Tmp02, Del

LOGICAL ::ChkSamePts
! ----------------------------------------------------

ChkSamePts = FALSE

Tmp01 = Pt01(1) - Pt02(1)
Tmp02 = Pt01(2) - Pt02(2)

Del = Tmp01 * Tmp01 + Tmp02 * Tmp02
Del = sqrt(Del)

IF (Del < hEps) ChkSamePts = TRUE
! ----------------------------------------------------

END FUNCTION ChkSamePts
! ------------------------------------------------------------------------------------------------------------
!                                     26. CHK : Pt Vtx
! ------------------------------------------------------------------------------------------------------------
FUNCTION ChkPtVtx(Pt, Vtx01, Vtx02, iPri)

IMPLICIT NONE

REAL    :: Pt(2)
REAL    :: Vtx01(2), Vtx02(2)
INTEGER :: iPri
LOGICAL :: ChkPtVtx

REAL :: Tmp01, Tmp02
! ----------------------------------------------------

ChkPtVtx = FALSE

Tmp01 = Pt(iPri) - Vtx01(iPri)
Tmp02 = Pt(iPri) - Vtx02(iPri)

IF (Tmp01 * Tmp02 < ZERO) ChkPtVtx = TRUE
! ----------------------------------------------------

END FUNCTION ChkPtVtx
! ------------------------------------------------------------------------------------------------------------
!                                     27. CHK : Pt Eqn Vtx
! ------------------------------------------------------------------------------------------------------------
FUNCTION ChkPtEqnVtx(Pt, Eqn, Vtx01, Vtx02, iPri)

IMPLICIT NONE

REAL    :: Pt(2)
REAL    :: Eqn(3)
REAL    :: Vtx01(2)
REAL    :: Vtx02(2)
INTEGER :: iPri
LOGICAL :: ChkPtEqnVtx

REAL :: Tmp, Tmp01, Tmp02
! ----------------------------------------------------

ChkPtEqnVtx = FALSE

Tmp = Eqn(3) - Eqn(1) * Pt(1) &
             - Eqn(2) * Pt(2)

IF (abs(Tmp) > hEps) RETURN

Tmp01 = Pt(iPri) - Vtx01(iPri)
Tmp02 = Pt(iPri) - Vtx02(iPri)

IF (Tmp01 * Tmp02 > ZERO) RETURN

ChkPtEqnVtx = TRUE
! ----------------------------------------------------

END FUNCTION ChkPtEqnVtx
! ------------------------------------------------------------------------------------------------------------
!                                     28. CHK : Eqn Two Vtx
! ------------------------------------------------------------------------------------------------------------
FUNCTION ChkEqnTwoVtx(Eqn, Vtx01, Vtx02)

IMPLICIT NONE

REAL :: Eqn(3)
REAL :: Vtx01(2)
REAL :: Vtx02(2)
REAL :: Tmp01, Tmp02, Tmp

LOGICAL :: ChkEqnTwoVtx
! ----------------------------------------------------

ChkEqnTwoVtx = FALSE

Tmp01 = Eqn(3) - Eqn(1) * Vtx01(1) - Eqn(2) * Vtx01(2)
Tmp02 = Eqn(3) - Eqn(1) * Vtx02(1) - Eqn(2) * Vtx02(2)
Tmp   = Tmp01 * Tmp02

!IF (abs(Tmp) < hEps) RETURN
IF (Tmp      < ZERO) ChkEqnTwoVtx = TRUE
! ----------------------------------------------------

END FUNCTION ChkEqnTwoVtx
! ------------------------------------------------------------------------------------------------------------
!                                     29. Find : Rod Rad to Vol
! ------------------------------------------------------------------------------------------------------------
FUNCTION FindRodRad2Vol(F2F, Rad)

IMPLICIT NONE

REAL :: FindRodRad2Vol
REAL :: F2F, Rad
REAL :: Ang1, Ang2

REAL, PARAMETER :: Eps = 1E-5
REAL, PARAMETER :: Sq3 = 1.73205080756888_8
REAL, PARAMETER :: rSq = 0.577350269189626_8
REAL, PARAMETER :: Pi6 = 0.52359877559830_8
! ----------------------------------------------------

IF (Rad > F2F * rSq - Eps) THEN
  FindRodRad2Vol = F2F * F2F * 0.5_8 * Sq3

ELSE IF (Rad < ZERO + Eps) THEN
  FindRodRad2Vol = ZERO

ELSE IF (Rad < F2F * 0.5_8 + Eps) THEN
  FindRodRad2Vol = Rad * Rad * Pi

ELSE
  Ang1 = acos(F2F * 0.5_8 / Rad)
  Ang2 = 2.0_8 * (Pi6 - Ang1)

  FindRodRad2Vol = 1.5_8 * F2F * F2F * tan(Ang1) & ! Triangle
                 + 3.0_8 * Rad * Rad * Ang2        ! Circular Sector
END IF
! ----------------------------------------------------

END FUNCTION FindRodRad2Vol
! ------------------------------------------------------------------------------------------------------------
!                                     30. Find : Rod Sub Rad
! ------------------------------------------------------------------------------------------------------------
FUNCTION FindRodSubRad(F2F, Rad2, Rad1, nDiv)
! Rad 2 > Rad 1

IMPLICIT NONE

INTEGER, PARAMETER :: nMaxFXR = 30

REAL :: FindRodSubRad(nMaxFXR)
REAL :: F2F, Rad2, Rad1
INTEGER :: nDiv

INTEGER, PARAMETER :: nnDiv = 10000

REAL :: dRad, tVol, dVol, pVol
REAL :: TmpRad(nnDiv)
INTEGER :: iTmp, iDiv, jDiv

REAL, PARAMETER :: rSq3 = 0.577350269189626_8
! ----------------------------------------------------

IF (Rad1 < ZERO) CALL terminate("FIND ROD SUB RAD")
IF (nDiv > nMaxFXR) CALL terminate("FIND ROD SUB RAD")

FindRodSubRad(1)      = Rad2
FindRodSubRad(nDiv+1) = Rad1
! ----------------------------------------------------
!               CASE : Rad < Half of F2F
! ----------------------------------------------------
IF (Rad2 < F2F * 0.5_8) THEN
  dVol = (Rad2 * Rad2 - Rad1 * Rad1) / real(nDiv)
  tVol = Rad2 * Rad2

  DO iDiv = 2, nDiv
    tVol = tVol - dVol

    FindRodSubRad(iDiv) = sqrt(tVol)
  END DO
  
  RETURN
END IF
! ----------------------------------------------------
!               CASE : Rad > PCH
! ----------------------------------------------------
IF (Rad2 > F2F * rSq3) THEN
  dRad = (Rad2 - Rad1) / real(nDiv)
  
  DO iDiv = 2, nDiv
    FindRodSubRad(iDiv) = FindRodSubRad(iDiv-1) - dRad
  END DO
  
  RETURN
END IF
! ----------------------------------------------------
!               CASE : Rad > Half of F2F
! ----------------------------------------------------
dRad = (Rad2 - Rad1) / real(nnDiv)

DO iTmp = 1, nnDiv
  TmpRad(iTmp) = Rad2 - iTmp * dRad
END DO

tVol = FindRodRad2Vol(F2F, Rad2)
dVol = (tVol - FindRodRad2Vol(F2F, Rad1)) / real(nDiv)
pVol = tVol - dVol

iDiv = 1
jDiv = nDiv - 1

DO iTmp = 1, nnDiv
  tVol = FindRodRad2Vol(F2F, TmpRad(iTmp))

  IF (tVol < pVol) THEN
    iDiv = iDiv + 1
    pVol = pVol - dVol

    FindRodSubRad(iDiv) = TmpRad(iTmp)
  END IF

  IF (iDiv > jDiv) EXIT
END DO

CALL ChkArrayZero(FindRodSubRad, 1, nMaxFXR)
! ----------------------------------------------------

END FUNCTION FindRodSubRad
! ------------------------------------------------------------------------------------------------------------
!                                     31. CHK ERR : Bndy - INTEGER
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexChkRange_INT(x, b1, b2, errmesg)

IMPLICIT NONE

INTEGER :: x, b1, b2
CHARACTER :: errmesg*(*)

INTEGER :: Tmp
! ----------------------------------------------------

Tmp = (x - b1) * (x - b2)

IF (Tmp > 0) CALL terminate(errmesg)

END SUBROUTINE HexChkRange_INT
! ------------------------------------------------------------------------------------------------------------
!                                     32. CHK ERR : Bndy - REAL
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexChkRange_REAL(x, b1, b2, errmesg)

IMPLICIT NONE

REAL :: x, b1, b2
CHARACTER :: errmesg*(*)

REAL :: Tmp
REAL, PARAMETER :: Eps = 1E-7
! ----------------------------------------------------

Tmp = (x - b1) * (x - b2)

IF (Tmp > Eps) CALL terminate(errmesg)

END SUBROUTINE HexChkRange_REAL
! ------------------------------------------------------------------------------------------------------------
!                                     33. CHK ERR - Equal - INTEGER
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexChkEqual_INT(x, y, i, errmesg)

IMPLICIT NONE

INTEGER :: x, y, i
CHARACTER :: errmesg*(*)
! ----------------------------------------------------

SELECT CASE (i)
CASE(1)
  IF (x .NE. y) CALL terminate(errmesg)
CASE(2)
  IF (x .EQ. y) CALL terminate(errmesg)
END SELECT
! ----------------------------------------------------

END SUBROUTINE HexChkEqual_INT
! ------------------------------------------------------------------------------------------------------------
!                                     34. CHK ERR - Equal - REAL
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexChkEqual_REAL(x, y, i, errmesg)

IMPLICIT NONE

REAL :: x, y
INTEGER :: i
CHARACTER :: errmesg*(*)

REAL :: Tmp
REAL, PARAMETER :: Eps = 1E-7
! ----------------------------------------------------

Tmp = abs(x - y)

SELECT CASE(i)
CASE(1)
  IF (Tmp > Eps) CALL terminate(errmesg)
CASE(2)
  IF (Tmp < Eps) CALL terminate(errmesg)
END SELECT
! ----------------------------------------------------

END SUBROUTINE HexChkEqual_REAL
! ------------------------------------------------------------------------------------------------------------
!                                     35. CHK ERR - INC - REAL
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexChkInc_REAL(x, y, errmesg)

IMPLICIT NONE

REAL :: x, y
CHARACTER :: errmesg*(*)

REAL :: Tmp
REAL, PARAMETER :: Eps = 1E-7
! ----------------------------------------------------

IF (x > y) CALL terminate(errmesg)

END SUBROUTINE HexChkInc_REAL
! ------------------------------------------------------------------------------------------------------------
!                                     36. FIND : Pt Ang
! ------------------------------------------------------------------------------------------------------------
FUNCTION FindPtAng(Pt)
! Center = ZERO

IMPLICIT NONE

REAL :: FindPtAng, Pt(2)
! ----------------------------------------------------

IF (ChkSameVal(Pt(1), ZERO)) THEN
  FindPtAng = PI_2 * (2._8 - SetSgn_REAL(Pt(2)))

  RETURN
END IF

! Results of atan = [-pi/2, pi/2]
! tan(pi + t) = -tan(t)
FindPtAng = atan(Pt(2)/Pt(1)) + PI_2 * (1._8 - SetSgn_REAL(Pt(1)))
! ----------------------------------------------------

END FUNCTION FindPtAng
! ------------------------------------------------------------------------------------------------------------
!                                     37. Find : Bndy Rad to Vol
! ------------------------------------------------------------------------------------------------------------
FUNCTION FindBndyRad2Vol(Rad, Vtx2, Vtx3, Vtx4, Eqn1, Eqn2)

IMPLICIT NONE

REAL :: Rad, FindBndyRad2Vol
REAL :: Vtx2(2), Vtx3(2), Vtx4(2), Eqn1(3), Eqn2(3)

REAL :: Tmp2, Tmp3, Tmp4, Ang, Bng, Vol1, Vol2, Pt(2, 2)
REAL :: Cnt(2), Vtx(2, 3)
! ----------------------------------------------------

Cnt  = ZERO
Tmp2 = FindPtLgh(Vtx2, Cnt)
Tmp3 = FindPtLgh(Vtx3, Cnt)
Tmp4 = FindPtLgh(Vtx4, Cnt)

IF (Rad > Tmp3) CALL terminate("FIND : BNDY RAD 2 VOL")

! sRad < Bndy
IF ((Rad < Tmp2).AND.(Rad < Tmp3)) THEN
  FindBndyRad2Vol = HALF * Rad * Rad * PI_6
  
  RETURN
END IF

! Vol 1
Ang = atan(Vtx3(2) / Vtx3(1))

IF (Rad < Tmp2) THEN
  Vol1 = HALF * Rad * Rad * (PI_6 - Ang)
ELSE
  Pt  = FindLineCircleIntSct(Eqn1, Rad)
  Bng = atan(Pt(2, 2) / Pt(1, 2)) - Ang
  Vtx = reshape((/ZERO, ZERO, Vtx2(1), Vtx2(2), Pt(1, 2), Pt(2, 2)/), (/2, 3/))
  
  Vol1 = HALF * Rad * Rad * Bng &
       + FindPolygonVol(3, Vtx)
END IF

! Vol2
IF (Rad < Tmp4) THEN
  Vol1 = HALF * Rad * Rad * Ang
ELSE
  Pt  = FindLineCircleIntSct(Eqn2, Rad)
  Bng = Ang - atan(Pt(2, 2) / Pt(1, 2))
  Vtx = reshape((/ZERO, ZERO, Vtx4(1), Vtx4(2), Pt(1, 2), Pt(2, 2)/), (/2, 3/))
  
  Vol1 = HALF * Rad * Rad * Bng &
       + FindPolygonVol(3, Vtx)
END IF

FindBndyRad2Vol = Vol1 + Vol2
! ----------------------------------------------------

END FUNCTION FindBndyRad2Vol
! ------------------------------------------------------------------------------------------------------------
!                                     38. Find : Bndy Rad to Vol
! ------------------------------------------------------------------------------------------------------------
FUNCTION FindPolygonVol(nVtx, Vtx)

IMPLICIT NONE

INTEGER :: nVtx
REAL    :: FindPolygonVol, Vtx(2, nVtx)

INTEGER :: iVtx
REAL    :: Tmp
! ----------------------------------------------------

Tmp = ZERO

DO iVtx = 1, nVtx - 1
  Tmp = Tmp + Vtx(1, iVtx)   * Vtx(2, iVtx+1)
  Tmp = Tmp - Vtx(1, iVtx+1) * Vtx(2, iVtx)
END DO

Tmp = Tmp + Vtx(1, nVtx) * Vtx(2, 1)
Tmp = Tmp - Vtx(1, 1)    * Vtx(2, nVtx)

FindPolygonVol = HALF * abs(Tmp)

END FUNCTION FindPolygonVol
! ------------------------------------------------------------------------------------------------------------
!                                     38. Find : Bndy Rad to Vol
! ------------------------------------------------------------------------------------------------------------
FUNCTION FindLineCircleIntSct(Eqn, Rad)
! x1 < x2

IMPLICIT NONE

REAL :: FindLineCircleIntSct(2, 2), Eqn(3), Rad, Pt(2)

REAL :: Tmp, a, b, c
! ----------------------------------------------------

IF (abs(Eqn(1)) < hEps) THEN
  IF (abs(Eqn(2)) < hEps) CALL terminate("FIND : LINE CIRCLE INT-SCT")
  
  Tmp = Eqn(3) / Eqn(2)
  
  IF (Rad < Tmp) CALL terminate("FIND : LINE CIRCLE INT-SCT")
  
  FindLineCircleIntSct(1, 1) = -sqrt(Rad*Rad - Tmp*Tmp)
  FindLineCircleIntSct(1, 2) =  sqrt(Rad*Rad - Tmp*Tmp)
  FindLineCircleIntSct(2, :) = Tmp
ELSE
  a   = 1 + Eqn(2)*Eqn(2) / Eqn(1)/Eqn(1)
  b   = - 2._8 * Eqn(2) * Eqn(3) / Eqn(1)/Eqn(1)
  c   = Eqn(3)*Eqn(3) / Eqn(1)/Eqn(1) - Rad*Rad
  Tmp = b*b - 4*a*c
  
  IF (Tmp < ZERO) CALL terminate("FIND : LINE CIRCLE INT-SCT")
  
  FindLineCircleIntSct(2, 1) = -(b + sqrt(Tmp)) / 2/a
  FindLineCircleIntSct(2, 2) = -(b - sqrt(Tmp)) / 2/a
  
  FindLineCircleIntSct(1, 1) = (Eqn(3) - Eqn(2) * FindLineCircleIntSct(2, 1)) / Eqn(1)
  FindLineCircleIntSct(1, 2) = (Eqn(3) - Eqn(2) * FindLineCircleIntSct(2, 2)) / Eqn(1)
  
  IF (FindLineCircleIntSct(1, 1) > FindLineCircleIntSct(1, 2)) THEN
    Pt = FindLineCircleIntSct(:, 1)
    
    FindLineCircleIntSct(:, 1) = FindLineCircleIntSct(:, 2)
    FindLineCircleIntSct(:, 2) = Pt
  END IF
END IF
! ----------------------------------------------------

END FUNCTION FindLineCircleIntSct
! ------------------------------------------------------------------------------------------------------------

! HexRT Routines

! ------------------------------------------------------------------------------------------------------------
!                                     01. ARRAY : 2D Sort
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE Array2DSORT(X, Y, n, npt, lAcending, lRmvRpt, iPriority)
REAL :: X(n), Y(n)
INTEGER :: n                   !array Size
INTEGER :: npt
INTEGER :: iPriority           !1 : X sorting, 2: Y sorting
LOGICAL :: lAcending           !Ascending Order Flag
LOGICAL :: lRmvRpt             !Remove Repetition Points
REAL, POINTER :: A(:), B(:)
!REAL :: A(n), B(n)

REAL :: temp
INTEGER :: i, j
ALLOCATE(A(n), B(n))
IF(iPriority .EQ. 1) THEN
  A = X; B = Y
ELSE
  A = Y; B = X;
ENDIF

IF(lAcending) THEN
  DO j = n, 1, -1
    DO i = 1, j - 1
      IF(A(i) > A(i+1)) THEN
        temp = A(i)
        A(i) = A(i+1); A(i+1) = temp
        temp = B(i)
        B(i) = B(i+1); B(i+1) = temp
      ENDIF
    ENDDO
  ENDDO
ELSE
  DO j = n, 1, -1
    DO i = 1, j - 1
      IF(A(i) < A(i+1)) THEN
        temp = A(i)
        A(i) = A(i+1); A(i+1) = temp
        temp = B(i)
        B(i) = B(i+1); B(i+1) = temp
      ENDIF
    ENDDO
  ENDDO
ENDIF

IF(lRmvRpt) THEN
  CALL Array2DRmvRpt(A, B, n, npt, 1)
  !CALL Array2DRmvRpt(A, B, n, npt, 2)
ELSE
  npt = n
ENDIF

IF(iPriority .EQ. 1) THEN
  X = A; Y = B
ELSE
  X = B; Y = A
ENDIF
DEALLOCATE(A, B)
END SUBROUTINE
! ------------------------------------------------------------------------------------------------------------
!                                     02. ARRAY : 2D Rmv Rpt
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE Array2DRmvRpt(X, Y, n, npt,  iPriority)
IMPLICIT NONE
REAL :: X(n), Y(n)
INTEGER :: n, npt, iPriority
!REAL :: tA(n), tB(n), A(n), B(n)
REAL, POINTER :: tA(:), tB(:), A(:), B(:)
INTEGER :: i

ALLOCATE(tA(n), tB(n), A(n), B(n))
tA = X; tB = Y
IF(iPriority .EQ. 2) THEN
  tA = Y; tB = X
ENDIF

A = 0; B = 0
npt = 1
A(1) = tA(1); B(1) = tB(1)
DO i = 2, n
  IF(abs(A(npt) - tA(i)) > 1.E-9) THEN
    npt = npt + 1
    A(npt) = tA(i); B(npt) = tB(i)
  ENDIF
ENDDO

X = A; Y = B
IF(iPriority .EQ. 2) THEN
  Y = A; X = B
ENDIF

DEALLOCATE(tA, tB, A, B)
END SUBROUTINE Array2DRmvRpt
! ------------------------------------------------------------------------------------------------------------

END MODULE HexUtil
