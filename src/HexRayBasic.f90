MODULE HexRayBasic

USE ioutil, ONLY : terminate

IMPLICIT NONE

INTEGER, PARAMETER :: aNxt(7, 2) = reshape((/0,  1,  1,  0, -1, -1,  0, & ! Nxt Asy xy Idx (iSurf, ix:iy)
                                            -1,  0,  1,  1,  0, -1,  0/), shape(aNxt))   ! iSurf = (NE, EE, SE, SW, WW, NW, SELF)
INTEGER, PARAMETER :: RotAux(2:7) = [3, 2, 0, 6, 5, 7]

CONTAINS

! ------------------------------------------------------------------------------------------------------------
!                                     01. HEX SET : Ray Param
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetRayParam()

USE PARAM,   ONLY : HALF, ZERO, PI
USE HexData, ONLY : nAzmAng, nPolAng, Del_Inp, Sq3Inv, aoPch, AzmAng, AzmDel, AzmTan, AzmSin, AzmCos, NumMray, &
                    AzmDel_X, AzmDel_Y, AzmWgt, hEps, AngMray, PI_3, PI_6
USE HexUtil, ONLY : ChkSameVal

IMPLICIT NONE

INTEGER :: nAzmAng_60

REAL :: AngleGuess, DelAngle, AngleAtan, PreviousAngle
REAL :: Angle, Delta, xDel, yDel

INTEGER :: NumRay, nside, ny, nside2
INTEGER :: iAng, iAng1, iAng2
! ----------------------------------------------------

IF (mod(nAzmAng, 6) .NE. 0) CALL terminate("# OF AZIMUTHAL ANGLES MUST BE MULTIPLE OF 6")

nAzmAng_60 = nAzmAng / 3
DelAngle   = PI_3 / nAzmAng_60

ALLOCATE (AzmDel    (nAzmAng)); AzmDel   = ZERO
ALLOCATE (AzmDel_X  (nAzmAng)); AzmDel_X = ZERO
ALLOCATE (AzmDel_Y  (nAzmAng)); AzmDel_Y = ZERO
ALLOCATE (AzmAng    (nAzmAng)); AzmAng   = ZERO
ALLOCATE (AzmWgt    (nAzmAng)); AzmWgt   = ZERO
ALLOCATE (AzmTan    (nAzmAng)); AzmTan   = ZERO
ALLOCATE (AzmSin    (nAzmAng)); AzmSin   = ZERO
ALLOCATE (AzmCos    (nAzmAng)); AzmCos   = ZERO
ALLOCATE (NumMray (0:nAzmAng)); NumMray  = 0
! ----------------------------------------------------
!               01. SET : Azm Ang & Del in 60
! ----------------------------------------------------
DO iAng = 1, nAzmAng_60
  AngleGuess = 0.5_8 * DelAngle + DelAngle * (iAng - 1)
  
  nside = nint(aoPch * sin(AngleGuess + PI_6) / Del_Inp)
  ny    = nint(aoPch * cos(AngleGuess)        / Del_Inp)
  
  AngleATAN = (1._8 * nside / (1._8 * ny) - sin(PI_6)) / cos(PI_6)
  
  Angle  = atan(AngleATAN)
  Delta  = aoPch * cos(angle) / (1._8 * ny)
  nside2 = nint(aoPch * sin(abs(PI_6 - Angle)) / Delta)
  
  AzmAng (iAng) = Angle
  AzmDel (iAng) = Delta
  NumMray(iAng) = ny + nside + nside2
END DO

DO iAng = 2, nAzmAng_60
  IF (ChkSameVal(AzmAng(iAng), AzmAng(iAng+1))) CALL terminate("AZIMUTHAL ANGLE OVERLAP")
END DO
! ----------------------------------------------------
!               02. CP : Azm Ang & Del
! ----------------------------------------------------
DO iAng = 1, nAzmAng_60
  iAng1 = nAzmAng_60 + iAng
  Angle = AzmAng(iAng) + PI_3
  Delta = AzmDel(iAng)
  
  NumRay = NumMray(iAng)
  
  AzmAng (iAng1) = Angle
  AzmDel (iAng1) = Delta
  NumMray(iAng1) = NumRay
  
  iAng2 = nAzmAng_60 + iang1
  Angle = Angle + PI_3
  
  AzmAng (iang2) = Angle
  AzmDel (iang2) = Delta  
  NumMray(iang2) = NumRay
END DO

NumMray(0) = SUM(NumMray(1:nAzmAng))
! ----------------------------------------------------
!               03. SET : Azm Val
! ----------------------------------------------------
DO iAng = 1, nAzmAng
  Angle = AzmAng(iAng)
  Delta = AzmDel(iAng)
  
  AzmTan  (iAng) = tan(Angle)
  AzmSin  (iAng) = sin(Angle)
  AzmCos  (iAng) = cos(Angle)
  AzmDel_Y(iAng) = Delta / abs(AzmCos(iAng))
  AzmDel_X(iAng) = Delta / abs(AzmSin(iAng))
END DO
! ----------------------------------------------------
!               04. SET : Azm Wgt
! ----------------------------------------------------
PreviousAngle = 0

DO iang = 1, nAzmAng_60
  Angle = 0.5_8 * (AzmAng(iang) + AzmAng(iAng+1))
  DelAngle = (Angle - PreviousAngle) / (2._8 * PI)
  
  iang1 = nAzmAng_60 + iang
  iang2 = nAzmAng_60 + iang1
  
  AzmWgt(iang)  = DelAngle
  AzmWgt(iang1) = DelAngle
  AzmWgt(iang2) = DelAngle
  
  PreviousAngle = Angle
END DO
! ----------------------------------------------------
!               05. Ang mRay Idx
! ----------------------------------------------------
ALLOCATE (AngMray (2, nAzmAng))

AngMray(1, 1) = 1
AngMray(2, 1) = NummRay(1)

DO iAng = 2, nAzmAng
  AngMray(1, iAng) = AngMray(2, iAng-1) + 1
  AngMray(2, iAng) = AngMray(2, iAng-1) + NumMray(iAng)
END DO
! ----------------------------------------------------

END SUBROUTINE HexSetRayParam
! ------------------------------------------------------------------------------------------------------------
!                                     02. HEX SET : Mod Ray
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetModRay()

USE PARAM,   ONLY : HALF
USE HexUtil, ONLY : SolveLineEqn, ChkPtVtx, ChkPtEqn, ChkSamePts
USE HexData, ONLY : aoF2F, aoPch, hmRay, NumMray, nAzmAng, AzmDel_Y, AzmTan, AsyEqn, AsyVtx, AngMray, Sq3

IMPLICIT NONE

INTEGER :: iAng, jAng, iRay, imRay, jmRay, iSurf, iDir, tDir, iTmp, iaX, iaY
INTEGER :: NumRay, NumRay2, nPt

REAL :: Delta_Y, Y_MAX, val1, Tmp, Tmp02, xDel, yDel
REAL :: RayEqn(3), C(2), pt(2, 2), Sol(2), OutPt(2), PtTmp(2)
REAL :: NxtAsyCnt(2) ! NE, EE, SE, SW, WW, NW

LOGICAL :: lSol, lChk
! ----------------------------------------------------

ALLOCATE(hmRay (NumMRay(0)))

imRay = 0
! ----------------------------------------------------
!               01. SET : mRay = 0 ~ 30 degree
! ----------------------------------------------------
DO iang = 1, nAzmAng / 6
  NumRay  = NumMray(iang)
  NumRay2 = NumRay / 2
  Delta_Y = AzmDel_Y(iAng)
  Y_MAX   = HALF * Delta_Y + Delta_Y * (NumRay2 - 1)
  
  CALL SetRayIntSct()
END DO
! ----------------------------------------------------
!               02. SET : mRay = 30 ~ 60 degree
! ----------------------------------------------------
DO iang = nAzmAng / 6 + 1, nAzmAng / 3
  NumRay  = NumMray(iang)
  Delta_Y = AzmDel_Y(iAng)
  
  jAng  = nAzmAng / 3 - iAng + 1
  jmRay = AngMray(2, jAng)
  
  PtTmp(1:2) = hmRay(jmRay)%Pt(1:2, 2)
  
  Tmp = AzmTan(iAng) * PtTmp(1) - PtTmp(2)
  
  DO iTmp = 1, 2*NumRay
    Tmp   = Tmp + Delta_Y
    Tmp02 = Tmp - AzmTan(iAng) * aoF2F * HALF
    
    IF (Tmp02 > aoPch * HALF) THEN
      Y_MAX = Tmp - Delta_Y
      EXIT
    END IF
  END DO
  
  IF (iTmp > 2 * NumRay) CALL terminate("FIND MRAY POINT")
  
  CALL SetRayIntSct()
END DO
! ----------------------------------------------------
!               03. SET : mRay = 60 ~ 90 degree
! ----------------------------------------------------
DO iang = nAzmAng / 3 + 1, nAzmAng / 2
  NumRay  = NumMray(iang)
  Delta_Y = AzmDel_Y(iAng)
  
  jAng  = nAzmAng / 3 - (iAng - nAzmAng / 3) + 1
  jmRay = AngMray(1, jAng)
  
  PtTmp(2) = -hmRay(jmRay)%Eq(3) / (Sq3 * hmRay(jmRay)%Eq(1) + 1)
  PtTmp(1) = -Sq3 * PtTmp(2)
  
  Tmp = AzmTan(iAng) * PtTmp(1) - PtTmp(2)
  
  DO iTmp = 1, 2*NumRay
    Tmp   = Tmp + Delta_Y
    Tmp02 = Tmp - AzmTan(iAng) * aoF2F * HALF
    
    IF (Tmp02 > aoPch * HALF) THEN
      Y_MAX = Tmp - Delta_Y
      EXIT
    END IF
  END DO
  
  IF (iTmp > 2 * NumRay) CALL terminate("FIND MRAY POINT")
  
  CALL SetRayIntSct()
END DO
! ----------------------------------------------------
!               04. SET : mRay = 90 ~ 180 degree
! ----------------------------------------------------
DO iang = nAzmAng / 2 + 1, nAzmAng
  ! ----------------------------
  !      1. Ray Position
  ! ----------------------------
  NumRay = NumMray(iang)
  jAng   = nAzmAng + 1 - iAng
  jmRay  = AngMray(1, jAng) - 1
  ! ----------------------------
  !      2. Ray Int Sct
  ! ----------------------------
  DO iRay = 1, NumRay
    imRay = imRay + 1
    jmRay = jmRay + 1
    
    hmRay(imRay)%AzmIdx = iAng
    
    hmRay(imRay)%Eq(1) = -hmRay(jmRay)%Eq(1)
    hmRay(imRay)%Eq(2) = -1._8
    hmRay(imRay)%Eq(3) =  hmRay(jmRay)%Eq(3)
    
    hmRay(imRay)%pt(2,1:2) =  hmRay(jmRay)%Pt(2, 1:2)
    hmRay(imRay)%pt(1,1:2) = -hmRay(jmRay)%Pt(1, 1:2)
  ENDDO
END DO
! ----------------------------------------------------
!               05. SET : Nxt Asy
! ----------------------------------------------------
DO imRay = 1, NumMray(0)
  Pt = hmRay(imRay)%Pt
  
  DO tDir = 1, 2
    iDir = 2*tDir - 3 ! Negative : y¢Ù, Positive : y¢Ö
    
    DO iSurf = 1, 6
      lChk = ChkPtEqn(pt(1:2, tDir), AsyEqn(1:3, iSurf))
      
      IF(lChk) EXIT
    END DO
    
    IF (.NOT. lChk) CALL terminate("SET MOD RAY SURF")
    
    hmRay(imray)%NxtAsy_Mov(1:2, iDir) = aNxt(iSurf, 1:2)
  END DO
END DO
! ----------------------------------------------------
!               06. SET : Nxt mRay
! ----------------------------------------------------
xDel = aoF2F
yDel = aoPch * 1.5_8

DO imRay = 1, NumMRay(0)
  iAng = hmRay(imRay)%AzmIdx
  
  iaX = hmRay(imray)%NxtAsy_Mov(1, 1)
  iaY = hmRay(imray)%NxtAsy_Mov(2, 1)
  
  NxtAsyCnt(1) =  xDel * iaX - xDel * 0.5 * iaY
  NxtAsyCnt(2) = -yDel * iaY
  
  OutPt = hmRay(imRay)%Pt(:, 2)
  OutPt = OutPt - NxtAsyCnt
  
  DO jmRay = AngMray(1, iAng), AngMray(2, iAng)
    lChk = ChkSamePts(OutPt(1:2), hmRay(jmRay)%pt(1:2, 1))
    
    IF (lChk) EXIT
  END DO
  
  IF(.NOT. lChk) CALL terminate("FIND NXT MRAY")
  
  hmRay(imRay)%NxtMray_Mov(1)  = jmRay
  hmRay(jmRay)%NxtMray_Mov(-1) = imRay
END DO

CONTAINS
! ----------------------------------------------------
!               SUB. SET : Ray Int Sct
! ----------------------------------------------------
SUBROUTINE SetRayIntSct()

DO iRay = 1, NumRay
  C(1) = 0._8
  C(2) = Y_MAX - Delta_Y * (iRay - 1)
  
  RayEqn(1) = AzmTan(iAng)
  RayEqn(2) = -1._8
  RayEqn(3) = C(1) * AzmTan(iAng) - C(2)
  
  nPt = 0
  
  DO iDir = 1, 6
    CALL SolveLineEqn(RayEqn(1:3), AsyEqn(1:3, iDir), Sol, lSol)
    
    IF (.NOT. lSol) CALL terminate("SET MRAY INT SCT")
    
    lSol = ChkPtVtx(Sol, AsyVtx(1:2, iDir), AsyVtx(1:2, iDir + 1), 2)
    
    IF (lSol) THEN
      nPt = nPt + 1
      
      Pt(1:2, nPt) = Sol(1:2)
    END IF
  END DO
  
  IF (nPt .NE. 2) CALL terminate("SET MRAY INT SCT")
  
  IF (Pt(2, 1) > Pt(2, 2)) THEN
    Sol(1:2)   = Pt(1:2, 1)
    Pt(1:2, 1) = Pt(1:2, 2)
    Pt(1:2, 2) = Sol(1:2)
  END IF
  
  imRay = imRay + 1
  
  hmRay(imRay)%AzmIdx = iAng
  hmRay(imRay)%Eq     = RayEqn
  hmRay(imRay)%pt     = pt
END DO

END SUBROUTINE SetRayIntSct
! ----------------------------------------------------

END SUBROUTINE HexSetModRay
! ------------------------------------------------------------------------------------------------------------
!                                     03. HEX SET : Mod Ray Nxt
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetModRayNxt()

USE HexData, ONLY : hLgc

IMPLICIT NONE

IF (hLgc%lAzmRef) THEN
  CALL HexSetModRayNxt_REF()
ELSE
  CALL HexSetModRayNxt_ROT()
END IF

END SUBROUTINE HexSetModRayNxt
! ------------------------------------------------------------------------------------------------------------
!                                     04. HEX SET : Mod Ray Nxt - REF
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetModRayNxt_REF()
! Sng Cel/Asy OR 060 Core

USE PARAM,   ONLY : TRUE, FALSE, ZERO
USE HexUtil, ONLY : SolveLineEqn, ChkPtVtx, FindRefAng, ChkPtEqn, &
                    ChkSamePts, CalPtLineSgn, ChkSameVal, SetSgn_REAL
USE HexType, ONLY : Type_HexModRay, Type_HexGeoTypInfo
USE HexData, ONLY : hmRay, hGeoTypInfo, hLgc, NumMray, AzmCos, AzmSin, AzmTan, nAzmAng, AngMray, &
                    nGeoTyp, AsyEqn, hEps

IMPLICIT NONE

INTEGER :: imRay, jmRay, iGeo, jGeo, iGeoSt, iBndy, jBndy, iDir, tDir, jDir, uDir, iAng, jAng, ixy(7)
LOGICAL :: lSol, lChk, lRef(7)
REAL    :: mp, Eqn(3), Sol(2), Ang(2), Slp(3), tMov(2), Pt01(2), Pt02(2)

INTEGER, POINTER :: RefAng(:, :, :)
LOGICAL, POINTER :: lNxtRef(:, :, :)
LOGICAL, POINTER :: lBndyPt(:, :, :)
REAL,    POINTER :: Pt(:, :, :, :)

TYPE(Type_HexModRay),     POINTER :: hmRay_Loc
TYPE(Type_HexGeoTypInfo), POINTER :: gInf_Loc
! ----------------------------------------------------

IF (hLgc%l360 .AND. hLgc%lRadVac) RETURN

ALLOCATE (Pt   (2, -1:1, 7, NumMray(0))); Pt      = ZERO
ALLOCATE (RefAng  (-1:1, 7, NumMray(0))); RefAng  = 0
ALLOCATE (lNxtRef (-1:1, 7, NumMray(0))); lNxtRef = FALSE
ALLOCATE (lBndyPt (-1:1, 7, NumMray(0))); lBndyPt = FALSE

ixy    = 1; lRef    = TRUE
ixy(1) = 0; lRef(1) = FALSE
! ----------------------------------------------------
!               01. SET : Int-Sct
! ----------------------------------------------------
DO imRay = 1, NumMray(0)
  hmRay_Loc => hmRay(imRay)
  
  Eqn    =  hmRay_Loc%Eq
  iAng   =  hmRay_Loc%AzmIdx
  Slp(1) =  AzmCos(iAng)
  Slp(2) =  AzmSin(iAng)
  Slp(3) =  AzmTan(iAng)
  mp     = -SetSgn_REAL(AzmTan(iAng))
  
  DO iGeo = 1, nGeoTyp
    gInf_Loc => hGeoTypInfo(iGeo)
    
    DO iBndy = 1, gInf_Loc%nBndy
      IF (lRef(iGeo) .EQV. (abs(gInf_Loc%Eqn(3, iBndy)) > hEps)) CYCLE
      
      ! FIND : Intersection
      CALL SolveLineEqn(Eqn(1:3), gInf_Loc%Eqn(1:3, iBndy), Sol(1:2), lSol)
      
      lSol = ChkPtVtx(Sol(1:2), gInf_Loc%Vtx(1:2, iBndy), gInf_Loc%Vtx(1:2, iBndy+1), gInf_Loc%Cor(iBndy))
      lSol = lSol .OR. ChkSamePts(Sol(1:2), gInf_Loc%Vtx(1:2, iBndy))
      lSol = lSol .OR. ChkSamePts(Sol(1:2), gInf_Loc%Vtx(1:2, iBndy+1))
      
      IF (.NOT. lSol) CYCLE
      
      ! FIND : Nxt
      tMov(1) = Sol(1) + 1._8
      tMov(2) = Sol(2) + Slp(3)
      
      iDir = mp * CalPtLineSgn(tMov, gInf_Loc%Eqn(1:3, iBndy), gInf_Loc%Cnt) ! Negative : y¢Ù, Positive : y¢Ö
      Ang  = FindRefAng(Slp(1:2), [-gInf_Loc%Eqn(2, iBndy), gInf_Loc%Eqn(1, iBndy)])
      
      IF(Ang(2) < ZERO) Ang = - Ang ! [0, Pi]
      
      DO jAng = 1, nAzmAng
        IF(jAng .EQ. iAng) CYCLE
        
        lChk = ChkSameVal(Ang(1), AzmCos(jAng))
        
        IF(lChk) EXIT
      END DO
      
      IF (jAng > nAzmAng) CALL terminate("FIND MRAY REF ANG")
      ! ----------------------------------------------------
      DO jBndy = 1, 6
        lChk = ChkPtEqn(Sol(1:2), AsyEqn(1:3, jBndy))
        
        IF (lChk) EXIT
      END DO
      
      hmRay_Loc%NxtAsy_Ref(1:2, iDir, iGeo) = aNxt(jBndy, 1:2) * ixy(iGeo)
      
      lBndyPt(iDir, iGeo, imRay) = (jBndy < 7) .AND. lRef(iGeo)
      lNxtRef(iDir, iGeo, imRay) = TRUE
      RefAng (iDir, iGeo, imRay) = jAng
      Pt(1:2, iDir, iGeo, imRay) = Sol(1:2)
    END DO
  END DO
END DO
! ----------------------------------------------------
!               02. SET : Nxt mRay
! ----------------------------------------------------
iGeoSt = 1 + nGeoTyp/7 ! 060 Core skips "iGeo = 1"

DO imRay = 1, NumMray(0)
  hmRay_Loc => hmRay(imRay)
  
  DO iGeo = iGeoSt, nGeoTyp
    DO tDir = 1, 2      ! 1 : y¢Ù, 2 : y¢Ö
      iDir = 2*tDir - 3 ! Negative : y¢Ù, Positive : y¢Ö
      
      IF (.NOT. lNxtRef(iDir, iGeo, imRay)) CYCLE
      
      jAng = RefAng(iDir, iGeo, imRay)
      
      IF (lBndyPt(iDir, iGeo, imRay)) THEN
        mp   = -1._8 ! Reflected to Next Asy
        jGeo = 1
      ELSE
        mp   = 1._8  ! Reflected to Self Asy
        jGeo = iGeo
      END IF
      
      Pt01 = mp * Pt(1:2, iDir, iGeo, imRay)
      
      Nested: DO jmRay = AngMray(1, jAng), AngMray(2, jAng)
        DO uDir = 1, 2      ! 1 : Backward, 2 : Forward
          jDir = 2*uDir - 3 ! Negative : y¢Ù, Positive : y¢Ö
          Pt02 = Pt(1:2, jDir, jGeo, jmRay)
          lChk = ChkSamePts(Pt01, Pt02)
          
          IF (lChk) EXIT Nested
        END DO
      END DO Nested
      
      IF (.NOT. lChk) CALL terminate("CONNECT MRAY")
      
      jDir = 3 - 2*uDir ! NOTICE : Sign Changes
      ! ----------------------------------------------------
      hmRay_Loc%NxtMray_Ref(iDir, iGeo) = jmRay
      hmRay_Loc%NxtDir_Ref (iDir, iGeo) = jDir
    END DO
  END DO
END DO

NULLIFY (hmRay_Loc, gInf_Loc, Pt, RefAng, lNxtRef)
! ----------------------------------------------------

END SUBROUTINE HexSetModRayNxt_REF
! ------------------------------------------------------------------------------------------------------------
!                                     04. HEX SET : Mod Ray Nxt - ROT
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetModRayNxt_ROT()
! Oonly for 60 VAC Core

USE PARAM,   ONLY : TRUE, FALSE, ZERO
USE HexUtil, ONLY : SolveLineEqn, ChkPtVtx, FindRotAng, ChkPtEqn, &
                    ChkSamePts, CalPtLineSgn, ChkSameVal, SetSgn_REAL, FindRotPt
USE HexType, ONLY : Type_HexModRay, Type_HexGeoTypInfo
USE HexData, ONLY : hmRay, hGeoTypInfo, hLgc, NumMray, AzmCos, AzmSin, AzmTan, nAzmAng, AngMray, &
                    nGeoTyp, PI_3, AsyEqn, hEps

IMPLICIT NONE

INTEGER :: imRay, jmRay, iGeo, jGeo, iBndy, jBndy, iDir, tDir, jDir, uDir, iAng, jAng
LOGICAL :: lSol, lChk, GeoCYC(3:7)
REAL    :: mp, Eqn(3), Sol(2), Ang(2), Slp(3), tMov(2), Pt01(2), Pt02(2)

INTEGER, POINTER :: RotAng(:, :, :)
LOGICAL, POINTER :: lNxtRot(:, :, :)
LOGICAL, POINTER :: lBndyPt(:, :, :)
REAL,    POINTER :: Pt(:, :, :, :)

TYPE(Type_HexModRay),     POINTER :: hmRay_Loc
TYPE(Type_HexGeoTypInfo), POINTER :: gInf_Loc
! ----------------------------------------------------

IF (hLgc%l360 .OR. hLgc%lRadRef) CALL terminate("WRONG B.C. WITH REF AZM BNDY CONDITION")

ALLOCATE (Pt   (2, -1:1, 7, NumMray(0))); Pt      = ZERO
ALLOCATE (RotAng  (-1:1, 7, NumMray(0))); RotAng  = 0
ALLOCATE (lNxtRot (-1:1, 7, NumMray(0))); lNxtRot = FALSE
ALLOCATE (lBndyPt (-1:1, 7, NumMray(0))); lBndyPt = FALSE

GeoCYC      = FALSE
GeoCYC(5:6) = TRUE
! ----------------------------------------------------
!               01. SET : Int-Sct
! ----------------------------------------------------
DO imRay = 1, NumMray(0)
  hmRay_Loc => hmRay(imRay)
  
  Eqn    =  hmRay_Loc%Eq
  iAng   =  hmRay_Loc%AzmIdx
  Slp(1) =  AzmCos(iAng)
  Slp(2) =  AzmSin(iAng)
  Slp(3) =  AzmTan(iAng)
  mp     = -SetSgn_REAL(AzmTan(iAng))
  
  DO iGeo = 3, nGeoTyp
    IF (GeoCYC(iGeo)) CYCLE
    
    gInf_Loc => hGeoTypInfo(iGeo)
    
    DO iBndy = 1, gInf_Loc%nBndy
      IF (abs(gInf_Loc%Eqn(3, iBndy)) > hEps) CYCLE
      
      ! FIND : Intersection
      CALL SolveLineEqn(Eqn(1:3), gInf_Loc%Eqn(1:3, iBndy), Sol(1:2), lSol)
      
      lSol = ChkPtVtx(Sol(1:2), gInf_Loc%Vtx(1:2, iBndy), gInf_Loc%Vtx(1:2, iBndy+1), gInf_Loc%Cor(iBndy))
      lSol = lSol .OR. ChkSamePts(Sol(1:2), gInf_Loc%Vtx(1:2, iBndy))
      lSol = lSol .OR. ChkSamePts(Sol(1:2), gInf_Loc%Vtx(1:2, iBndy+1))
      
      IF (.NOT. lSol) CYCLE
      
      ! FIND : Nxt
      tMov(1) = Sol(1) + 1._8
      tMov(2) = Sol(2) + Slp(3)
      
      iDir = mp * CalPtLineSgn(tMov, gInf_Loc%Eqn(1:3, iBndy), gInf_Loc%Cnt) ! Negative : y¢Ù, Positive : y¢Ö
      Ang  = FindRotAng(Slp(1:2), [-gInf_Loc%Eqn(2, iBndy), gInf_Loc%Eqn(1, iBndy)], PI_3)
      
      IF(Ang(2) < ZERO) Ang = - Ang ! [0, Pi]
      
      DO jAng = 1, nAzmAng
        IF(jAng .EQ. iAng) CYCLE
        
        lChk = ChkSameVal(Ang(1), AzmCos(jAng))
        
        IF(lChk) EXIT
      END DO
      
      IF (jAng > nAzmAng) CALL terminate("FIND MRAY ROT ANG")
      ! ----------------------------------------------------
      DO jBndy = 2, 6
        lChk = ChkPtEqn(Sol(1:2), AsyEqn(1:3, jBndy))
        
        IF (lChk) EXIT
      END DO
      
      hmRay_Loc%NxtAsy_Ref(1:2, iDir, iGeo) = aNxt(RotAux(jBndy), 1:2)
      
      lBndyPt(iDir, iGeo, imRay) = jBndy < 7
      lNxtRot(iDir, iGeo, imRay) = TRUE
      RotAng (iDir, iGeo, imRay) = jAng
      Pt(1:2, iDir, iGeo, imRay) = Sol(1:2)
    END DO
  END DO
END DO
! ----------------------------------------------------
!               02. SET : Nxt mRay
! ----------------------------------------------------
DO imRay = 1, NumMray(0)
  hmRay_Loc => hmRay(imRay)
  
  DO iGeo = 3, nGeoTyp
    IF (GeoCYC(iGeo)) CYCLE
    
    DO tDir = 1, 2      ! 1 : y¢Ù, 2 : y¢Ö
      iDir = 2*tDir - 3 ! Negative : y¢Ù, Positive : y¢Ö
      
      IF (.NOT. lNxtRot(iDir, iGeo, imRay)) CYCLE
      
      jAng = RotAng(iDir, iGeo, imRay)
      
      IF (lBndyPt(iDir, iGeo, imRay)) THEN
        mp = -1._8 ! Reflected to Next Asy
      ELSE
        mp = 1._8  ! Reflected to Self Asy
      END IF
      
      Pt01 = FindRotPt(Pt(1:2, iDir, iGeo, imRay), PI_3) * mp
      
      Nested: DO jmRay = AngMray(1, jAng), AngMray(2, jAng)
        DO jGeo = 3, 4
          DO uDir = 1, 2      ! 1 : Backward, 2 : Forward
            jDir = 2*uDir - 3 ! Negative : y¢Ù, Positive : y¢Ö
            Pt02 = Pt(1:2, jDir, jGeo, jmRay)
            lChk = ChkSamePts(Pt01, Pt02)
            
            IF (lChk) EXIT Nested
          END DO
        END DO
      END DO Nested
      
      IF (.NOT. lChk) CALL terminate("CONNECT MRAY")
      
      jDir = 3 - 2*uDir ! NOTICE : Sign Changes
      ! ----------------------------------------------------
      hmRay_Loc%NxtMray_Ref(iDir, iGeo) = jmRay
      hmRay_Loc%NxtDir_Ref (iDir, iGeo) = jDir
    END DO
  END DO
END DO

NULLIFY (hmRay_Loc, gInf_Loc, Pt, RotAng, lNxtRot, lBndyPt)
! ----------------------------------------------------

END SUBROUTINE HexSetModRayNxt_ROT
! ------------------------------------------------------------------------------------------------------------

END MODULE HexRayBasic