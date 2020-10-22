! ------------------------------------------------------------------------------------------------------------
!                                     01. HEX SET : Gap Cel Bss
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetGapCelBss()

USE PARAM,   ONLY : TRUE
USE geom,    ONLY : nGapType, nGapPinType
USE ioutil,  ONLY : terminate

USE HexType, ONLY : Type_HexGapCel, Type_HexGapCelBss, nMaxFXR
USE HexData, ONLY : hLgc, gCel, gCelBss, ngBss, hLgc, vFxr, vAsyTyp, vMat1, vMat2, GapPin, hAsyTypInfo, vzSt, vzEd
USE HexUtil, ONLY : HexChkRange_INT

IMPLICIT NONE

INTEGER :: igCel, jgCel, kgCel, igBss, iMsh, jgPin
LOGICAL :: lNewBss, HexCompNewGapBss
! ----------------------------------------------------

IF (hLgc%lSngCel) RETURN

ALLOCATE (gCelBss (nGapType))
! ----------------------------------------------------
!               01. SET : Vyg Gap Cel
! ----------------------------------------------------
IF (hLgc%lvyg) THEN
  jgPin = hAsyTypInfo(vAsyTyp)%gTyp
  
  CALL HexPushGapCel(nGapType-1, GapPin(jgPin)%iCel(1)) ! VYG Cel : Edge, iz is fixed as 1
  CALL HexPushGapCel(nGapType,   GapPin(jgPin)%iCel(1)) ! VYG Cel : Corn
  
  CALL HexChkRange_INT(vFXR, 1, gCel(nGapType-1)%nFXR, "WRONG VFXR")
  
  gCel(nGapType-1)%xMix(gCel(nGapType-1)%nFXR - vFXR + 1) = vMat1
  gCel(nGapType)  %xMix(:)                                = vMat2
  
  CALL HexPushGapPin(nGapPinType-1, jgPin)
  CALL HexPushGapPin(nGapPinType,   jgPin)
  
  GapPin(nGapPinType-1)%iCel(vzSt:vzEd) = nGapType-1 ! VYG Pin : Edge
  GapPin(nGapPinType)  %iCel(vzSt:vzEd) = nGapType   ! VYG Pin : Corn
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

USE allocs
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

CALL dmalloc(ggCelBss%iCel, nGapType)

ggCelBss%nCel    = 1
ggCelBss%iCel(1) = iG

CALL dmalloc(ggCelBss%sHgt, ggCelBss%nSub + 1)

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

USE allocs
USE PARAM,   ONLY : HALF, ZERO
USE HexType, ONLY : Type_HexGapCelBss
USE HexData, ONLY : Sq3Inv, gCelBss

IMPLICIT NONE

INTEGER :: igBss
INTEGER :: iFXR, jFXR, iSub, jSub, iMsh, jMsh, nMsh, HexCalNumVtxHor

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

CALL dmalloc(gBss_Loc%sVol, 2, max(gBss_Loc%nMsh01, gBss_Loc%nMsh02))
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
!                                     04. HEX PUSH : Gap Cel
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexPushGapCel(igCel, jgCel)

USE PARAM,   ONLY : TRUE
USE HexType, ONLY : nMaxFXR
USE HexData, ONLY : gCel

IMPLICIT NONE

INTEGER :: igCel, jgCel
! ----------------------------------------------------

gCel(igCel)%luse = TRUE

gCel(igCel)%nPin  = gCel(jgCel)%nPin
gCel(igCel)%nFXR  = gCel(jgCel)%nFXR
gCel(igCel)%pF2F  = gCel(jgCel)%pF2F
gCel(igCel)%aiF2F = gCel(jgCel)%aiF2F

gCel(igCel)%xHgt(1:nMaxFXR) = gCel(jgCel)%xHgt(1:nMaxFXR)
gCel(igCel)%xDiv(1:nMaxFXR) = gCel(jgCel)%xDiv(1:nMaxFXR)
gCel(igCel)%xMix(1:nMaxFXR) = gCel(jgCel)%xMix(1:nMaxFXR)
! ----------------------------------------------------

END SUBROUTINE HexPushGapCel
! ------------------------------------------------------------------------------------------------------------
!                                     04. HEX PUSH : Gap Pin
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexPushGapPin(igPin, jgPin)

USE allocs
USE PARAM,   ONLY : TRUE
USE geom,    ONLY : nZ
USE HexData, ONLY : GapPin

IMPLICIT NONE

INTEGER :: igPin, jgPin
! ----------------------------------------------------

GapPin(igPin)%luse = TRUE
GapPin(igPin)%lGap = TRUE

CALL dmalloc(GapPin(igPin)%iCel, nZ)

GapPin(igPin)%iCel(1:nZ) = GapPin(jgPin)%iCel(1:nZ)
! ----------------------------------------------------

END SUBROUTINE HexPushGapPin
! ------------------------------------------------------------------------------------------------------------