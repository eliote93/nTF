MODULE HexCelBasis_Rod
  
CONTAINS
! ------------------------------------------------------------------------------------------------------------
!                                     01. HEX SET : Rod Cel Bss
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetRodCelBss()

USE PARAM,   ONLY : TRUE
USE GEOM,    ONLY : nCellType, nPinType, nZ
USE ioutil,  ONLY : terminate
USE HexData, ONLY : hCel, hCelBss, ncBss, RodPin, Sq3, Sq3Inv

IMPLICIT NONE

INTEGER :: iCel, jCel, iPin, iz, nUse
LOGICAL :: lBss
! ----------------------------------------------------

! ----------------------------------------------------
!               01. REMOVE : Un-used 
! ----------------------------------------------------
nUse = 0

DO iCel = 1, nCellType
  IF (.NOT. hCel(iCel)%luse) THEN
    DO iPin = 1, nPinType
      DO iz = 1, nZ
        IF (RodPin(iPin)%iCel(iz) .EQ. iCel) RodPin(iPin)%iCel(iz) = -1
      END DO
    END DO
    
    CYCLE
  END IF
  
  nUse = nUse + 1
  
  hCel(nUse) = hCel(iCel)
  
  DO iPin = 1, nPinType
    DO iz = 1, nZ
      IF (RodPin(iPin)%iCel(iz) .EQ. iCel) RodPin(iPin)%iCel(iz) = nUse
    END DO
  END DO
END DO

nCellType = nUse

ALLOCATE (hCelBss (nCellType))
! ----------------------------------------------------
!               02. SET : 1st Cel Geo Basis
! ----------------------------------------------------
IF (.NOT. hCel(1)%luse) CALL terminate("1ST ROD CEL MUST BE USED")

CALL HexPushRodCelBss(hCelBss(1), hCel(1), 1, 1)
! ----------------------------------------------------
!               03. FIND : Nxt Cel Geo Basis
! ----------------------------------------------------
ncBss = 1

DO iCel = 2, nCellType
  lBss = TRUE
  
  DO jCel = 1, ncBss
    lBss = HexCompNewRodBss(hCelBss(jCel), hCel(iCel))
    
    IF (.NOT. lBss) THEN
      hCel(iCel)%icBss = jCel
      
      hCelBss(jCel)%nCel = hCelBss(jCel)%nCel + 1
      
      hCelBss(jCel)%iCel(hCelBss(jCel)%nCel) = iCel
      
      EXIT
    END IF
  END DO
  
  IF (.NOT. lBss) CYCLE
  
  ncBss = ncBss + 1
  
  CALL HexPushRodCelBss(hCelBss(ncBss), hCel(iCel), ncBss, iCel)
END DO
! ----------------------------------------------------
!               04. SET : Inn Sub Ring DATA
! ----------------------------------------------------
DO iCel = 1, ncBss
  CALL HexSetRodBssSubRng_Inn (iCel)
  CALL HexSetRodBssSubRng_Bndy(iCel)
END DO
! ----------------------------------------------------

END SUBROUTINE HexSetRodCelBss
! ------------------------------------------------------------------------------------------------------------
!                                     02. PUSH : Tmp Bss
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexPushRodCelBss(hhCelBss, hhCel, iB, iC)

USE PARAM,   ONLY : ZERO
USE GEOM,    ONLY : nCellType
USE HexType, ONLY : Type_HexRodCel, Type_HexRodCelBss
USE HexData, ONLY : Sq3

IMPLICIT NONE

TYPE(Type_HexRodCelBss) :: hhCelBss
TYPE(Type_HexRodCel)    :: hhCel

INTEGER :: iB, iC
INTEGER :: iFXR
! ------------------------

hhCelBss%nFXR = hhCel%nFXR
hhCelBss%nPin = hhCel%nPin
hhCelBss%nSct = hhCel%nSct
hhCelBss%nCel = 1

hhCelBss%pF2F  = hhCel%pF2F
hhCelBss%aiF2F = hhCel%aiF2F
hhCelBss%pPch  = hhCel%pF2F / Sq3
hhCelBss%cVol  = hhCelBss%pPch * hhCelBss%pPch * 1.5 * Sq3

hhCelBss%xRad(1:hhCel%nFXR) = hhCel%xRad(1:hhCel%nFXR)
hhCelBss%xDiv(1:hhCel%nFXR) = hhCel%xDiv(1:hhCel%nFXR)

hhCel%icBss = iB

ALLOCATE (hhCelBss%iCel (nCellType))

hhCelBss%iCel    = 0
hhCelBss%iCel(1) = iC

hhCelBss%nSub = sum(hhCel%xDiv(1:hhCel%nFXR))
hhCelBss%nMsh = sum(hhCel%xDiv(1:hhCel%nFXR)) * hhCel%nSct

ALLOCATE (hhCelBss%sRad (hhCelBss%nSub+1))
ALLOCATE (hhCelBss%sVol (0:2, hhCelBss%nSct, hhCelBss%nSub))

hhCelBss%sRad = ZERO
hhCelBss%sVol = ZERO

hhCelBss%luse = hhCelBss%luse .OR. hhCel%luse

END SUBROUTINE HexPushRodCelBss
! ------------------------------------------------------------------------------------------------------------
!                                     03. COMP : New Rod Bss
! ------------------------------------------------------------------------------------------------------------
FUNCTION HexCompNewRodBss(hhCelBss, hhCel)

USE PARAM,   ONLY : epsm5, TRUE, FALSE
USE HexType, ONLY : Type_HexRodCel, Type_HexRodCelBss

IMPLICIT NONE

TYPE(Type_HexRodCelBss) :: hhCelBss
TYPE(Type_HexRodCel)    :: hhCel

INTEGER :: iFXR
LOGICAL :: l01, l02
LOGICAL :: HexCompNewRodBss
! ------------------------

HexCompNewRodBss = TRUE

l01 = hhCelBss%nFXR .NE. hhCel%nFXR
l02 = hhCelBss%nSct .NE. hhCel%nSct

IF (l01.OR.l02) RETURN

l01 = abs(hhCelBss%pF2F - hhCel%pF2F) > epsm5

IF (l01) RETURN

DO iFXR = 1, hhCel%nFXR
  l01 =    (hhCelBss%xDiv(iFXR) - hhCel%xDiv(iFXR)) .NE. 0
  l02 = abs(hhCelBss%xRad(iFXR) - hhCel%xRad(iFXR)) > epsm5
  
  IF (l01.OR.l02) RETURN
END DO

l01 =    (hhCelBss%nPin  - hhCel%nPin) .NE. 0
l02 = abs(hhCelBss%aiF2F - hhCel%aiF2F) > epsm5

IF (l01.OR.l02) THEN
  IF (hhCel%nPin .EQ. 0) THEN
    hhCel%nPin  = hhCelBss%nPin
    hhCel%aiF2F = hhCelBss%aiF2F
  ELSE
    RETURN
  END IF
END IF

HexCompNewRodBss = FALSE

END FUNCTION HexCompNewRodBss
! ------------------------------------------------------------------------------------------------------------
!                                     04. HEX SET : Rod Bss Sub Rng - Inn
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetRodBssSubRng_Inn(iCel)

USE PARAM,   ONLY : HALF, PI
USE HexType, ONLY : Type_HexRodCelBss, nMaxFXR
USE HexData, ONLY : hCelBss, Sq3
USE HexUtil, ONLY : FindRodRad2Vol, FindRodSubRad

IMPLICIT NONE

INTEGER :: iCel

INTEGER :: iFXR, iSub, jSub
INTEGER :: nFXR

REAL :: pF2F, Tmp, Vol
REAL :: rSub(nMaxFXR)

TYPE(Type_HexRodCelBss), POINTER :: hBss_Loc
! ----------------------------------------------------

hBss_Loc => hCelBss(iCel)
! ----------------------------------------------------
!               01. xVol
! ----------------------------------------------------
nFXR = hBss_Loc%nFXR
pF2F = hBss_Loc%pF2F
Vol  = pF2F * pF2F * Sq3 * HALF

DO iFXR = 1, nFXR
  Tmp = FindRodRad2Vol(pF2F, hBss_Loc%xRad(iFXR+1))
  
  hBss_Loc%xVol(:, iFXR) = Vol - Tmp
  
  Vol = Tmp
END DO
! ----------------------------------------------------
!               02. sRad
! ----------------------------------------------------
jSub = 0

DO iFXR = 1, nFXR
  rSub = FindRodSubRad(pF2F, hBss_Loc%xRad(iFXR), hBss_Loc%xRad(iFXR + 1), hBss_Loc%xDiv(iFXR))
  
  DO iSub = 1, hBss_Loc%xDiv(iFXR)
    jSub = jSub + 1
    
    hBss_Loc%sRad(jSub) = rSub(iSub)
  END DO
END DO

hBss_Loc%xRad(1) = sqrt(pF2F * pF2F * Sq3 * HALF / Pi)
hBss_Loc%sRad(1) = sqrt(pF2F * pF2F * Sq3 * HALF / Pi)
! ----------------------------------------------------
!               03. sVol
! ----------------------------------------------------
Vol = pF2F * pF2F * Sq3 * HALF

DO iSub = 1, hBss_Loc%nSub
  Tmp = FindRodRad2Vol(pF2F, hBss_Loc%sRad(iSub+1))
  
  hBss_Loc%sVol(:, :, iSub) = (Vol - Tmp) / hBss_Loc%nSct
  
  Vol = Tmp
END DO

NULLIFY (hBss_Loc)
! ----------------------------------------------------

END SUBROUTINE HexSetRodBssSubRng_Inn
! ------------------------------------------------------------------------------------------------------------
!                                     05. HEX SET : Rod Bss Sub Rng - Bndy
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetRodBssSubRng_Bndy(iCel)
! # of Sct is fixed as 12

USE PARAM,   ONLY : HALF, PI, ZERO, ONE
USE ioutil,  ONLY : terminate
USE HexType, ONLY : Type_HexRodCelBss
USE HexData, ONLY : hCelBss, Sq3, Sq3Inv, hLgc
USE HexUtil, ONLY : SetEqn, SolveLineEqn, FindRodRad2Vol, FindBndyRad2Vol

IMPLICIT NONE

INTEGER :: iCel, iSub, iFXR, nSub, mSub
LOGICAL :: lSol

REAL :: pF2F, pPch, Hgt, Are1, Are2, Vol, Tmp1, Tmp2
REAL :: Cnt(2), Vtx1(2), Vtx2(2), Vtx3(2), Vtx4(2), Vtx5(2)
REAL :: Eqn1(3), Eqn2(3)

TYPE(Type_HexRodCelBss), POINTER :: hBss_Loc
! ----------------------------------------------------

IF (hLgc%lSngCel) RETURN

hBss_Loc => hCelBss(iCel)

pF2F = hBss_Loc%pF2F
pPch = hBss_Loc%pPch
! ----------------------------------------------------
!               01. SET : Math Value
! ----------------------------------------------------
Hgt  = (hBss_Loc%aiF2F - 3 * (hBss_Loc%nPin - 1) * pPch) * HALF
Are1 = Hgt * Hgt * Sq3Inv * HALF
Are2 = pF2F * HALF * (2._8 * Hgt - pPch * HALF) * HALF - Are1

Cnt     =  ZERO
Vtx1(1) =  ZERO
Vtx1(2) =  2.0_8 * Hgt * Sq3Inv
Vtx2(1) =  Hgt
Vtx2(2) =  Hgt * Sq3Inv
Vtx4(1) =  pPch
Vtx4(2) =  ZERO
Vtx5(1) =  pPch * HALF
Vtx5(2) = -pF2F * HALF

Eqn1 = SetEqn(Vtx1, Vtx2, Cnt)
Eqn2 = SetEqn(Vtx4, Vtx5, Cnt)

CALL SolveLineEqn(Eqn1, Eqn2, Vtx3, lSol)
! ----------------------------------------------------
!               02. SET : sVol
! ----------------------------------------------------
hBss_Loc%sVol(1, 4:9,  :) = hBss_Loc%sVol(0, 4:9,  :)
hBss_Loc%sVol(2, 1:4,  :) = hBss_Loc%sVol(0, 1:4,  :)
hBss_Loc%sVol(2, 9:12, :) = hBss_Loc%sVol(0, 9:12, :)

! Are1
Vol = Are1 * 12._8

DO iSub = 1, hBss_Loc%nSub
  Tmp1 = FindRodRad2Vol(2.0_8 * Hgt, hBss_Loc%sRad(iSub+1))
  Tmp2 = (Vol - Tmp1) / hBss_Loc%nSct
  
  hBss_Loc%sVol(1,  1:2,  iSub) = Tmp2
  hBss_Loc%sVol(1, 11:12, iSub) = Tmp2
  hBss_Loc%sVol(2,  6:7,  iSub) = Tmp2
  
  Vol = Tmp1
END DO

! Are2
Vol = Are2

DO iSub = 1, hBss_Loc%nSub
  Tmp1 = FindBndyRad2Vol(hBss_Loc%sRad(iSub+1), Vtx2, Vtx3, Vtx4, Eqn1, Eqn2)
  Tmp2 = Vol - Tmp1
  
  hBss_Loc%sVol(1,  3, iSub) = Tmp2
  hBss_Loc%sVol(1, 10, iSub) = Tmp2
  hBss_Loc%sVol(2,  5, iSub) = Tmp2
  hBss_Loc%sVol(2,  8, iSub) = Tmp2
  
  Vol = Tmp1
END DO
! ----------------------------------------------------
!               02. SET : xVol
! ----------------------------------------------------
nSub = 1

DO iFXR = 1, hBss_Loc%nFXR
  mSub = nSub + hBss_Loc%xDiv(iFXR) - 1
  
  hBss_Loc%xVol(1, iFXR) = sum(hBss_Loc%sVol(1, :, nSub:mSub))
  hBss_Loc%xVol(2, iFXR) = sum(hBss_Loc%sVol(2, :, nSub:mSub))
  
  nSub = mSub + 1
END DO

NULLIFY (hBss_Loc)
! ----------------------------------------------------

END SUBROUTINE HexSetRodBssSubRng_Bndy
! ------------------------------------------------------------------------------------------------------------

END MODULE HexCelBasis_Rod