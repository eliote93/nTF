MODULE HexAsyRayConst

USE HexData, ONLY : hAsyTypInfo, hCelBss, gCelBss
  
IMPLICIT NONE

CONTAINS

! ------------------------------------------------------------------------------------------------------------
!                                     01. SET : Pin Typ Dat
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetRayPinDat(RayPinInfo, icBss)

USE HexType, ONLY : Type_HexAsyTypInfo, Type_HexRayPinInfo
USE HexData, ONLY : nGeoTyp, mpTypNumNgh
USE HexUtil, ONLY : SetEqn, ChkArrayZero, RotPt

IMPLICIT NONE

TYPE(Type_HexRayPinInfo), POINTER :: RayPinInfo(:)
TYPE(Type_HexAsyTypInfo), POINTER :: aInf_Loc
TYPE(Type_HexRayPinInfo), POINTER :: rPin_Loc

INTEGER :: icBss, iPin, iBndy, iGeo, ivTyp, nBndy
REAL    :: Ang
! ----------------------------------------------------

aInf_Loc => hAsyTypInfo(hCelBss(icBss)%iaTyp)

DO iPin = 1, aInf_Loc%nTotPin(1)
  rPin_Loc => RayPinInfo(iPin)
  
  rPin_Loc%Cnt(1:2) = aInf_Loc%PinCnt(1:2, iPin)
  rPin_Loc%Vtx(1:2, 1:7, 1:nGeoTyp) = aInf_Loc%mpVtx(1:2, 1:7, 1:nGeoTyp, iPin)
  
  IF (iPin > aInf_Loc%nRodPin(1)) rPin_Loc%iCel = 2
  
  DO iGeo = 1, nGeoTyp
    ivTyp = aInf_Loc%PinVtxTyp(iGeo, iPin)
    Ang   = aInf_Loc%PinVtxAng(iGeo, iPin)
    
    IF (ivTyp < 1) CYCLE
    
    nBndy = mpTypNumNgh(ivTyp)
    
    DO iBndy = 1, nBndy
      rPin_Loc%Eqn(1:3, iBndy, iGeo) = SetEqn(rPin_Loc%Vtx(1:2, iBndy, iGeo), rPin_Loc%Vtx(1:2, iBndy+1, iGeo), &
                                              rPin_Loc%Cnt(1:2))
    END DO
    
    rPin_Loc%nBndy (iGeo) = nBndy
    rPin_Loc%VtxTyp(iGeo) = ivTyp
    rPin_Loc%VtxAng(iGeo) = Ang
  END DO
  
  CALL ChkArrayZero(rPin_Loc%Cnt, 1, 2)
  CALL ChkArrayZero(rPin_Loc%Eqn, 1, 3, 1, 6, 1, 7)
  CALL ChkArrayZero(rPin_Loc%Vtx, 1, 2, 1, 7, 1, 7)
END DO

NULLIFY (aInf_Loc, rPin_Loc)
! ----------------------------------------------------

END SUBROUTINE HexSetRayPinDat
! ------------------------------------------------------------------------------------------------------------
!                                     05. SET : Rod Cel 12
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetRayRodCel_12(RayCel, icBss)

USE allocs
USE PARAM,   ONLY : PI, HALF, ONE, ZERO
USE ioutil,  ONLY : terminate
USE HexType, ONLY : Type_HexRayCel, Type_HexRodCelBss
USE HexData, ONLY : PI_6, PI_2
USE HexUtil, ONLY : SetNewPt, ChkArrayZero, SetSgn_REAL

IMPLICIT NONE

TYPE(Type_HexRayCel),    POINTER :: RayCel(:)
TYPE(Type_HexRodCelBss), POINTER :: hBss_Loc

INTEGER :: icBss
INTEGER :: iSub, iDir, iMsh, iEqn, jEqn
INTEGER :: nSub, nMsh, nEqn

REAL :: Ang, Tmp, Pts(2)

TYPE(Type_HexRayCel), POINTER :: RodCel
! ----------------------------------------------------

hBss_Loc => hCelBss(icBss)

IF (hBss_Loc%nSct .NE. 12) CALL terminate("# OF SECTORS MUST BE 12")

nSub = hBss_Loc%nSub
nMsh = hBss_Loc%nMsh

RodCel => RayCel(1)

CALL dmalloc(RodCel%Eqn,  5, 5 + nSub)
CALL dmalloc(RodCel%iEqnCor, 5 + nSub)

CALL dmalloc(RodCel%nMshEqn,      nMsh)
CALL dmalloc(RodCel%MshEqnLst, 4, nMsh)
CALL dmalloc(RodCel%MshEqnVal, 4, nMsh)

RodCel%iEqnCor = 1 ! NOTICE : x-cor
RodCel%nMsh    = hBss_Loc%nMsh
! ----------------------------------------------------
!               01. SET : Rod Cel Eqn
! ----------------------------------------------------
! Sct Line
DO iDir = 1, 6
  Ang = PI_2 - (iDir - 1) * PI_6
  
  RodCel%Eqn(1, iDir) = -sin(Ang)
  RodCel%Eqn(3, iDir) =  cos(Ang)
END DO

RodCel%iEqnCor(1) = 2 ! NOTICE : y-cor

! Cyl
DO iSub = 1, nSub - 1
  RodCel%Eqn(2, 6 + iSub) = 1._8
  RodCel%Eqn(4, 6 + iSub) = 1._8
  RodCel%Eqn(5, 6 + iSub) = hBss_Loc%sRad(iSub+1) * hBss_Loc%sRad(iSub+1)
END DO

CALL ChkArrayZero(RodCel%Eqn, 1, 5, 1, 5 + nSub)

RodCel%nEqn = 5 + nSub
! ----------------------------------------------------
!               02. SET : Msh Eqn Dat
! ----------------------------------------------------
iMsh = 0

DO iSub = 1, nSub
  DO iDir = 1, 12
    iMsh = iMsh + 1
    nEqn = 2
    
    ! Sector Line
    Ang  = PI_2 - (iDir - 1) * PI_6 - PI / 12._8 ! Middle Ang
    Pts  = SetNewPt(ONE, Ang)
    
    RodCel%MshEqnLst(1, iMsh) = iDir - 6 * (iDir / 7)     ! Left Sct Line
    RodCel%MshEqnLst(2, iMsh) = iDir + 1 - 6 * (iDir / 6) ! Right Sct Line
    
    DO iEqn = 1, 2
      jEqn = RodCel%MshEqnLst(iEqn, iMsh)
      Tmp  = RodCel%Eqn(1, jEqn) * Pts(1) + RodCel%Eqn(3, jEqn) * Pts(2)
      
      RodCel%MshEqnVal(iEqn, iMsh) = -SetSgn_REAL(Tmp) ! NOTICE : Minus
    END DO
    
    ! Circle
    DO iEqn = 0, 1
      jEqn = 5 + iSub + iEqn
      
      IF ((jEqn < 7).OR.(jEqn > 5 + nSub)) CYCLE
      
      nEqn = nEqn + 1
      
      RodCel%MshEqnLst(nEqn, iMsh) = jEqn
      RodCel%MshEqnVal(nEqn, iMsh) = 1._8 - 2._8 * iEqn ! 1 : Msh is under Cyl Eqn
    END DO
    
    RodCel%nMshEqn(iMsh) = nEqn
  END DO
END DO
! ----------------------------------------------------
!               03. SET : Fsr Idx
! ----------------------------------------------------
CALL dmalloc(RodCel%MshIdx, 7, nMsh)

iMsh = 0

! Inn & Bndy 01
DO iSub = 1, nSub
  RodCel%MshIdx(1, 12 * iSub - 11) = 2 * iSub - 1 ! iVtx = 1 : Inn 060
  RodCel%MshIdx(1, 12 * iSub)      = 2 * iSub
  
  DO iDir = 1, 12
    iMsh = iMsh + 1
    
    RodCel%MshIdx(3:4, iMsh) = iMsh ! iVtx = 3 & 4 : Inn 360 & Bndy 01
    
    IF ((iDir > 1).AND.(iDir < 8)) CYCLE
    
    ! iVtx = 2 : Inn 180
    RodCel%MshIdx(2, iMsh) = 6 * (iSub - 1) +  iDir - 6 * (iDir / 8)
  END DO
END DO

iMsh = 0

! Bndy 02
DO iSub = 1, nSub
  DO iDir = 1, 12
    iMsh = iMsh + 1
    
    RodCel%MshIdx(5, iMsh) = iMsh + 5 - 12 * (iDir / 8)
    
    IF ((iDir > 1).AND.(iDir < 8)) CYCLE
    
    RodCel%MshIdx(6, iMsh) = 6 * (iSub - 1) + iDir + 5 - 12 * (iDir / 8)
    RodCel%MshIdx(7, iMsh) = 6 * (iSub - 1) - iDir + 2 + 12 * (iDir / 8)
  END DO
END DO

NULLIFY (RodCel, hBss_Loc)
! ----------------------------------------------------

END SUBROUTINE HexSetRayRodCel_12
! ------------------------------------------------------------------------------------------------------------
!                                     06. SET : Gap Cel
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetRayGapCel(RayCel, icBss)

USE allocs
USE HexType, ONLY : Type_HexRayCel, Type_HexGapCelBss
USE PARAM,   ONLY : ZERO, ONE, PI, HALF
USE HexData, ONLY : Sq3, PI_6

IMPLICIT NONE

TYPE(Type_HexRayCel),    POINTER :: RayCel(:)
TYPE(Type_HexGapCelBss), POINTER :: gBss_Loc

INTEGER :: icBss
INTEGER :: iEqn, jEqn, iSub, iDir, iMsh, jMsh, kMsh, lMsh
INTEGER :: gHor, nSub, nMsh, nEqn, Bndy08, Bndy09
REAL    :: Del

REAL, PARAMETER :: Sq3Hlf = Sq3 * HALF

TYPE(Type_HexRayCel), POINTER :: GapCel
! ----------------------------------------------------

gBss_Loc => gCelBss(hCelBss(icBss)%igBss)
GapCel   => RayCel(2)

gHor = gBss_Loc%nVtxHor * 2 - 4
nSub = gBss_Loc%nSub
nMsh = gHor * nSub

CALL dmalloc(GapCel%Eqn,  5, nSub + gHor - 2)
CALL dmalloc(GapCel%iEqnCor, nSub + gHor - 2)

CALL dmalloc(GapCel%nMshEqn,      nMsh)
CALL dmalloc(GapCel%MshEqnLst, 4, nMsh)
CALL dmalloc(GapCel%MshEqnVal, 4, nMsh)

GapCel%iEqnCor = 1 ! NOTICE : x-cor
GapCel%nMsh    = nMsh
! ----------------------------------------------------
!               01. SET : Gap Cel Eqn - NE
! ----------------------------------------------------
Del = gBss_Loc%pF2F / gBss_Loc%nHor

! Vertical Line : L ¡æ R
DO iDir = -(gBss_Loc%nVtxHor - 3), gBss_Loc%nVtxHor - 3
  jEqn = iDir + gBss_Loc%nVtxHor - 2
  
  GapCel%Eqn(1, jEqn) = -Sq3Hlf
  GapCel%Eqn(3, jEqn) = HALF
  GapCel%Eqn(5, jEqn) = -Del * iDir
END DO

! Horizontal Line : T ¡æ B
DO iSub = 2, nSub
  jEqn = iSub + gHor - 2
  
  GapCel%Eqn(1, jEqn) = HALF
  GapCel%Eqn(3, jEqn) = Sq3Hlf
  GapCel%Eqn(5, jEqn) = gBss_Loc%sHgt(iSub) - gBss_Loc%sHgt(1) * HALF
END DO

GapCel%nEqn = nSub + gHor - 2
! ----------------------------------------------------
!               02. SET : Msh Eqn Dat
! ----------------------------------------------------
iMsh = 0

DO iSub = 1, nSub
  DO iDir = 1, gHor
    iMsh = iMsh + 1
    nEqn = 0
    
    ! Vertical line
    DO iEqn = 0, 1
      jEqn = iDir - 1 + iEqn
      
      IF ((jEqn < 1).OR.(jEqn > gHor - 1)) CYCLE
      
      nEqn = nEqn + 1
      
      GapCel%MshEqnLst(nEqn, iMsh) = jEqn
      GapCel%MshEqnVal(nEqn, iMsh) = 1._8 - 2._8 * iEqn ! 1 : Msh is on the Right of Line Eqn
    END DO
    
    ! Horizontal line
    DO iEqn = 0, 1
      jEqn = gHor - 2 + iEqn + iSub
      
      IF ((jEqn < gHor).OR.(jEqn > nSub + gHor - 2)) CYCLE
      
      nEqn = nEqn + 1
      
      GapCel%MshEqnLst(nEqn, iMsh) = jEqn
      GapCel%MshEqnVal(nEqn, iMsh) = 1._8 - 2._8 * iEqn ! 1 : Msh is under Line Eqn
    END DO
    
    GapCel%nMshEqn(iMsh) = nEqn
  END DO
END DO
! ----------------------------------------------------
!               04. SET : Msh Idx
! ----------------------------------------------------
CALL dmalloc0(GapCel%MshIdx, 8, 10, 1, nMsh)

iMsh = 0 ! Msh Idx
jMsh = 0 ! Msh Idx in iVtx = 8  : Gap 01 (NE & L)
kMsh = 0 ! Msh Idx in iVtx = 9  : Gap 01 (NE & R)
lMsh = 0 ! Msh Idx in iVtx = 10 : Gap 02 (Clock-wise)

Bndy08 = gHor/2 - gBss_Loc%nHor/4
Bndy09 = gHor/2 + gBss_Loc%nHor/4 + 1

DO iSub = 1, nSub
  DO iDir = 1, gHor ! L ¡æ R
    iMsh = iMsh + 1
    
    ! iVtx = 9
    IF (iDir < Bndy09) THEN
      kMsh = kMsh + 1
      
      GapCel%MshIdx(9, iMsh) = kMsh
    END IF
    
    ! iVtx = 10
    IF ((Bndy08 < iDir).AND.(iDir < Bndy09)) THEN
      lMsh = lMsh + 1
      
      GapCel%MshIdx(10, iMsh) = lMsh
    END IF
  END DO
  
  ! iVtx = 8
  DO iDir = gHor, 1, -1 ! R ¡æ L
    IF (iDir > Bndy08) THEN
      jMsh = jMsh + 1
      
      GapCel%MshIdx(8, iMsh + iDir - gHor) = jMsh
    END IF
  END DO
END DO

NULLIFY (GapCel, gBss_Loc)
! ----------------------------------------------------

END SUBROUTINE HexSetRayGapCel
! ------------------------------------------------------------------------------------------------------------

END MODULE HexAsyRayConst