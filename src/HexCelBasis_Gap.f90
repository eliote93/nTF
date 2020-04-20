MODULE HexCelBasis_Gap
  
IMPLICIT NONE
  
CONTAINS
! ------------------------------------------------------------------------------------------------------------
!                                     01. HEX SET : Gap Cel Bss
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetGapCelBss()

USE PARAM,   ONLY : TRUE
USE geom,    ONLY : nGapType, nGapPinType
USE ioutil,  ONLY : terminate

USE HexType, ONLY : Type_HexGapCel, Type_HexGapCelBss, nMaxFXR
USE HexData, ONLY : hLgc, gCel, gCelBss, ngBss, hLgc, vFxr, vAsyTyp, vMat, GapPin, hAsyTypInfo, vzSt, vzEd
USE HexUtil, ONLY : HexChkRange_INT

IMPLICIT NONE

INTEGER :: igCel, jgCel, kgCel, igBss, iMsh
LOGICAL :: lNewBss
! ----------------------------------------------------

IF (hLgc%lSngCel) RETURN

ALLOCATE (gCelBss (nGapType))
! ----------------------------------------------------
!               01. SET : Vyg Gap Cel
! ----------------------------------------------------
IF (hLgc%lvyg) THEN
  gCel(nGapType) = gCel(GapPin(hAsyTypInfo(vAsyTyp)%gTyp)%iCel(1)) ! iz is fixed as 1
  
  CALL HexChkRange_INT(vFXR, 1, gCel(nGapType)%nFXR, "WRONG VFXR")
  
  gCel(nGapType)%xMix(gCel(nGapType)%nFXR - vFXR + 1) = vMat
  gCel(nGapType)%luse = TRUE
  
  GapPin(nGapPinType)%iCel(vzSt:vzEd) = nGapType
  GapPin(nGapPinType)%luse = TRUE
END IF
! ----------------------------------------------------
!               02. SET : 1st Cel Geo Basis
! ----------------------------------------------------
DO igCel = 1, nGapType
  IF (.NOT. gCel(igCel)%luse) CYCLE
  
  CALL HexPushGapCelBss(gCelBss(1), gCel(igCel), 1, igCel)
  
  EXIT
END DO

IF (gCelBss(1)%nCel .EQ. 0) CALL terminate("NONE OF GAP CEL TYPE IS USED")
! ----------------------------------------------------
!               03. FIND : Nxt Cel Geo Basis
! ----------------------------------------------------
ngBss = 1

DO igCel = 1, nGapType
  IF (.NOT. gCel(igCel)%luse) CYCLE
    
  DO jgCel = 1, ngBss
    lNewBss = HexCompNewGapBss(gCelBss(jgCel), gCel(igCel))
    
    IF (lNewBss) CYCLE
    
    DO kgCel = 1, gCelBss(jgCel)%nCel
      IF (igCel .EQ. gCelBss(jgCel)%iCel(kgCel)) EXIT
    END DO
    
    IF (igCel .LE. gCelBss(jgCel)%nCel) CYCLE
    
    gCel(igCel)%igBss = jgCel
    
    gCelBss(jgCel)%nCel = gCelBss(jgCel)%nCel + 1
    
    gCelBss(jgCel)%iCel(gCelBss(jgCel)%nCel) = igCel
    
    EXIT
  END DO
  
  IF (.NOT. lNewBss) CYCLE
  
  ngBss = ngBss + 1
  
  CALL HexPushGapCelBss(gCelBss(ngBss), gCel(igCel), ngBss, igCel)
END DO
! ----------------------------------------------------
!               04. SET : Sub Data
! ----------------------------------------------------
DO igBss = 1, ngBss
  CALL HexSetGapBssSubRng(igBss)
END DO
! ----------------------------------------------------

END SUBROUTINE HexSetGapCelBss
! ------------------------------------------------------------------------------------------------------------
!                                     02. PUSH : Gap Cel Bss
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexPushGapCelBss(ggCelBss, ggCel, iB, iG)

USE PARAM,   ONLY : ZERO
USE GEOM,    ONLY : nGapType
USE HexType, ONLY : Type_HexGapCel, Type_HexGapCelBss

IMPLICIT NONE

TYPE(Type_HexGapCelBss) :: ggCelBss
TYPE(Type_HexGapCel)    :: ggCel

INTEGER :: iB, iG
! ----------------------------------------------------

ggCel%igBss = iB

ggCelBss%nPin = ggCel%nPin
ggCelBss%nSub = sum(ggCel%xDiv(1:ggCel%nFXR))
ggCelBss%nFXR = ggCel%nFXR

ggCelBss% pF2F = ggCel% pF2F
ggCelBss%aiF2F = ggCel%aiF2F

ggCelBss%xHgt(1:ggCel%nFXR) = ggCel%xHgt(1:ggCel%nFXR)
ggCelBss%xDiv(1:ggCel%nFXR) = ggCel%xDiv(1:ggCel%nFXR)

ALLOCATE (ggCelBss%iCel (nGapType))

ggCelBss%nCel    = 1
ggCelBss%iCel    = 0
ggCelBss%iCel(1) = iG

ALLOCATE (ggCelBss%sHgt (ggCelBss%nSub + 1))

ggCelBss%sHgt = ZERO
! ----------------------------------------------------

END SUBROUTINE HexPushGapCelBss
! ------------------------------------------------------------------------------------------------------------
!                                     03. COMP : New Gap Bss
! ------------------------------------------------------------------------------------------------------------
FUNCTION HexCompNewGapBss(ggCelBss, ggCel)

USE PARAM,   ONLY : epsm5, TRUE, FALSE
USE HexType, ONLY : Type_HexGapCel, Type_HexGapCelBss

IMPLICIT NONE

TYPE(Type_HexGapCelBss) :: ggCelBss
TYPE(Type_HexGapCel)    :: ggCel

INTEGER :: iFXR
LOGICAL :: l01, l02
LOGICAL :: HexCompNewGapBss
! ----------------------------------------------------

HexCompNewGapBss = TRUE

l01 =    (ggCelBss%nFXR  - ggCel%nFXR) .NE. 0
l02 = abs(ggCelBss%aiF2F - ggCel%aiF2F) > epsm5

IF (l01 .OR. l02) RETURN

l01 =    (ggCelBss%nPin - ggCel%nPin) .NE. 0
l02 = abs(ggCelBss%pF2F - ggCel%pF2F) > epsm5

IF (l01 .OR. l02) RETURN

DO iFXR = 1, ggCel%nFXR
  l01 =    (ggCelBss%xDiv(iFXR) - ggCel%xDiv(iFXR)) .NE. 0
  l02 = abs(ggCelBss%xHgt(iFXR) - ggCel%xHgt(iFXR)) > epsm5
  
  IF (l01.OR.l02) RETURN
END DO

HexCompNewGapBss = FALSE
! ----------------------------------------------------

END FUNCTION HexCompNewGapBss
! ------------------------------------------------------------------------------------------------------------
!                                     04. CAL : Num Vtx Hor
! ------------------------------------------------------------------------------------------------------------
FUNCTION HexCalNumVtxHor(aiF2F, pF2F, nPin, nHor)

USE PARAM,   ONLY : HALF
USE HexData, ONLY : Sq3Inv

IMPLICIT NONE

INTEGER :: HexCalNumVtxHor

REAL :: aiF2F, pF2F
INTEGER :: nPin, nHor

REAL :: Lgh, Tmp, Del
INTEGER :: iMsh
! ----------------------------------------------------

Lgh = aiF2F * Sq3Inv - (nPin - 2) * pF2F
Lgh = Lgh * HALF
Tmp = 0._8
Del = pF2F / nHor

DO iMsh = 1, nHor
  Tmp = Tmp + Del
  
  IF (.NOT.(Tmp > Lgh)) CYCLE
  
  HexCalNumVtxHor = iMsh
  
  EXIT
END DO
! ----------------------------------------------------

END FUNCTION HexCalNumVtxHor
! ------------------------------------------------------------------------------------------------------------
!                                     03. HEX SET : Gap Bss 
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetGapBssSubRng(igBss)

USE PARAM,   ONLY : HALF, ZERO
USE HexType, ONLY : Type_HexGapCelBss
USE HexData, ONLY : Sq3Inv, gCelBss

IMPLICIT NONE

INTEGER :: igBss
INTEGER :: iFXR, jFXR, iSub, jSub, iMsh, jMsh, nMsh

REAL :: Vol, Tmp, Hgt, aiPch

TYPE(Type_HexGapCelBss), POINTER :: gBss_Loc
! ----------------------------------------------------

gBss_Loc => gCelBss(igBss)
! ----------------------------------------------------
!               01. nMsh
! ----------------------------------------------------
gBss_Loc%nVtxHor = HexCalNumVtxHor(gBss_Loc%aiF2F, gBss_Loc%pF2F, gBss_Loc%nPin, gBss_Loc%nHor)
gBss_Loc%nMsh01  = gBss_Loc%nSub * gBss_Loc%nVtxHor
gBss_Loc%nMsh02  = gBss_Loc%nSub * gBss_Loc%nHor / 2

ALLOCATE (gBss_Loc%sVol (2, max(gBss_Loc%nMsh01, gBss_Loc%nMsh02))); gBss_Loc%sVol = ZERO
! ----------------------------------------------------
!               02. sHgt
! ----------------------------------------------------
jSub = gBss_Loc%nSub + 1

DO iFXR = 1, gBss_Loc%nFXR
  jFXR = gBss_Loc%nFXR - iFXR + 1
  
  DO iSub = 1, gBss_Loc%xDiv(jFXR)
    jSub = jSub - 1
    
    gBss_Loc%sHgt(jSub) = gBss_Loc%sHgt(jSub + 1) + (gBss_Loc%xHgt(jFXR) - gBss_Loc%xHgt(jFXR + 1)) / gBss_Loc%xDiv(jFXR) ! Out ¡æ Inn
  END DO
END DO
! ----------------------------------------------------
!               03. sVol 01
! ----------------------------------------------------
aiPch = gBss_Loc%aiF2F * Sq3Inv
nMsh  = gBss_Loc%nHor  * (gBss_Loc%nPin - 2) + 2 * (gBss_Loc%nVtxHor - 1) ! # of Rectangular Mshs per sRng
jMsh  = 0

DO iSub = 1, gBss_Loc%nSub
  Hgt = gBss_Loc%sHgt(iSub) - gBss_Loc%sHgt(iSub + 1)
  Tmp = Hgt * (aiPch + Sq3Inv * (gBss_Loc%sHgt(iSub) + gBss_Loc%sHgt(iSub + 1)))
  Vol = Hgt * gBss_Loc%pF2F / real(gBss_Loc%nHor)
  
  jMsh = jMsh + 1
  
  gBss_Loc%sVol(1, jMsh) = HALF * (Tmp - Vol * nMsh) ! Trapezoid Msh
  
  DO iMsh = 2, gBss_Loc%nVtxHor
    jMsh = jMsh + 1
    
    gBss_Loc%sVol(1, jMsh) = Vol
  END DO
END DO
! ----------------------------------------------------
!               04. sVol 02
! ----------------------------------------------------
jMsh = 0

DO iSub = 1, gBss_Loc%nSub
  Vol = (gBss_Loc%pF2F / gBss_Loc%nHor) * (gBss_Loc%sHgt(iSub) - gBss_Loc%sHgt(iSub + 1))
  
  DO iMsh = 1, gBss_Loc%nHor / 2
    jMsh = jMsh + 1
    
    gBss_Loc%sVol(2, jMsh) = Vol
  END DO
END DO

NULLIFY (gBss_Loc)
! ----------------------------------------------------

END SUBROUTINE HexSetGapBssSubRng
! ------------------------------------------------------------------------------------------------------------

END MODULE HexCelBasis_Gap