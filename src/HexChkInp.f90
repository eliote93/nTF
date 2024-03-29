MODULE HexChk
  
IMPLICIT NONE

CONTAINS
! ------------------------------------------------------------------------------------------------------------  
SUBROUTINE HexChkInp()

USE CNTL,    ONLY : nTracerCntl
USE PARAM,   ONLY : TRUE, FALSE, HALF
USE geom,    ONLY : nCellType, nPinType, nGapType, nGapPinType, nAsyType0, nVssTyp, nZ
USE HexUtil, ONLY : ChkSameVal, ChkRange_INT, ChkEqual_INT, ChkEqual_REAL, ChkInc_REAL
USE HexType, ONLY : Type_HexRodCel, Type_HexGapCel, Type_HexAsyTypInfo
USE IOUTIL,  ONLY : terminate
USE HexData

USE Material_Mod, ONLY : nMixType
USE BenchXs,      ONLY : nXslType

IMPLICIT NONE

INTEGER :: iaTyp, iCel, iPin, ix, iy, iz, iBss, jBss, iFXR, nTmp, nPin, nMix
REAL    :: pF2F, aiF2F

REAL, PARAMETER :: rSq3 = 0.577350269189626_8

TYPE(Type_HexRodCel),     POINTER :: hCel_Loc
TYPE(Type_HexGapCel),     POINTER :: gCel_Loc
TYPE(Type_HexAsyTypInfo), POINTER :: aInf_Loc
! ----------------------------------------------------

hLgc%lNoRef = hLgc%l360 .AND. hLgc%lRadVac 

IF (hLgc%lSngCel .AND. nTracerCntl%lCMFD) CALL terminate("SINGLE CELL CMFD IS NOT AVAILABLE")
IF (hLgc%lSngCel .AND. .NOT.(nZ.EQ.1)) CALL terminate("3-D SINGLE CELL PROBLEM IS NOT AVAILABLE")
IF (hLgc%lSngCel) RETURN

! Core
IF (nVA.LT.0 .AND. hLgc%lRadRef) CALL terminate("VOID ASY SHOULD NOT EXIST IN RADIAL REF")

IF (hLgc%l360) THEN
  nTmp = nhAsy - (3 * nAsyCore * (nAsyCore - 1) + 1 - nVA)
ELSE IF (hLgc%l060) THEN
  nTmp = nhAsy - (nAsyCore * (nAsyCore + 1) / 2 - nVA)
END IF

IF (nTmp .NE. 0) CALL terminate("RAD CONF ASY INP")

! Asy
DO iaTyp = 1, nAsyType0
  aInf_Loc => hAsyTypInfo(iaTyp)
  
  IF (.NOT. aInf_Loc%luse) CYCLE
  
  nPin  = aInf_Loc%nPin
  pF2F  = aInf_Loc%pF2F
  aiF2F = aInf_Loc%aiF2F
  
  IF (nPin .LT. 2) CALL terminate("# of Pin")
  IF (abs(aiF2F - aoF2F) .LT. hEps) CALL terminate("UNDERC-ONSTRUCTION - NO GAP")
  
  CALL ChkRange_INT(aInf_Loc%gTyp, 1, nGapPinType, "GAP PIN 1")
  
  ! Rod Pin
  DO ix = 1, 2 * nPin - 1
    DO iy = 1, 2 * nPin - 1
      iPin = aInf_Loc%PinIdx(iy, ix)
      
      IF (iPin < 1) CYCLE
      
      CALL ChkRange_INT(iPin, 1, nPinType, "ASY PIN")
      
      DO iz = 1, nZ
        iCel = RodPin(iPin)%iCel(iz)
        
        CALL ChkEqual_REAL(pF2F,  hCel(iCel)%pF2F,  1, "LENGTH OF PIN IN ASY & CEL ARE DIFFERENT")
        CALL ChkEqual_INT (nPin,  hCel(iCel)%nPin,  1, "NUMBER OF PIN IN ASY & CEL ARE DIFFERENT")
        CALL ChkEqual_REAL(aiF2F, hCel(iCel)%aiF2F, 1, "LENGTH OF GAP IN ASY & CEL ARE DIFFERENT")
      END DO
    END DO
  END DO
  
  ! Gap Pin
  DO iz = 1, nZ
    iCel = GapPin(aInf_Loc%gTyp)%iCel(iz)
    
    CALL ChkEqual_REAL(aiF2F, gCel(iCel)%aiF2F, 1, "LENGTH OF GAP IN ASY & CEL ARE DIFFERENT")
  END DO
  
  ! Cstt.
  IF (aInf_Loc%cstTyp .NE. 0) THEN
    CALL ChkRange_INT(aInf_Loc%cstTyp, 1, nGapPinType, "CORNER STIFFENER Type")
    CALL ChkRange_INT(aInf_Loc%cstNum, 1, nPin / 2,    "CORNER STIFFENER Number")
  END IF
END DO

! Pin
DO iPin = 1, nPinType
  IF (.NOT. RodPin(iPin)%luse) CYCLE
  
  DO iz = 1, nZ
    CALL ChkRange_INT(RodPin(iPin)%iCel(iz), 1, nCellType, "ROD PIN")
  END DO
END DO

DO iPin = 1, nGapPinType
  IF (.NOT. GapPin(iPin)%luse) CYCLE
  
  DO iz = 1, nZ
    CALL ChkRange_INT(GapPin(iPin)%iCel(iz), 1, nGapType, "GAP PIN 2")
  END DO
END DO

! Rod Cel
IF (nCellType .EQ. 0) CALL terminate("CEL DOES NOT EXIST")

IF (nTracerCntl%lXsLib) THEN
  nMix = nMixType
ELSE
  nMix = nXslType
END IF

DO iCel = 1, nCellType
  hCel_Loc => hCel(iCel)
  
  IF (.NOT. hCel_Loc%luse) CYCLE
  
  IF (hCel_Loc%nPin .LT. 2) CALL terminate("# of Pin")
    
  DO iFXR = 1, hCel_Loc%nFXR
    CALL ChkRange_INT(hCel_Loc%xMix(iFXR), 1, nMix, "ROD CEL MIXTURE")
  END DO
  
  DO iFXR = 2, hCel_Loc%nFXR
    CALL ChkInc_REAL(hCel_Loc%xRad(iFXR + 1), hCel_Loc%xRad(iFXR), "ROD CEL RADIUS")
  END DO
  
  IF (hCel_Loc%xRad(2) > hCel_Loc%pF2F * rSq3) hCel_Loc%xDiv(1) = 1
  
  IF (hCel_Loc%nSct .NE. 12) CALL terminate("CEL # OF SECTORS")
END DO

!IF (nTracerCntl%lDomainDcmp) THEN
!  DO iCel = 2, nCellType
!    IF (ChkSameVal(hCel(1)%pF2F, hCel(iCel)%pF2F)) CYCLE
!    
!    CALL terminate("INVALID GEOMETRY FOR DCMP.")
!  END DO
!END IF

! Gap Cel
IF (nGapType .EQ. 0) CALL terminate("GAP DOES NOT EXIST")

DO iCel = 1, nGapType
  gCel_Loc => gCel(iCel)
  
  IF (.NOT. gCel_Loc%luse) CYCLE
  
  DO iFXR = 1, gCel_Loc%nFXR
    CALL ChkRange_INT(gCel_Loc%xMix(iFXR), 1, nMix, "GAP CEL MIXTURE")
  END DO
  
  DO iFXR = 1, gCel_Loc%nFXR
    CALL ChkInc_REAL(gCel_Loc%xHgt(iFXR + 1), gCel_Loc%xHgt(iFXR), "GAP CEL HIGHT")
  END DO
END DO

! Etc.
DO iBss = 1, nVssTyp
  CALL ChkInc_REAL(hVss(iBss)%Rad(1), hVss(iBss)%Rad(2), "VESSEL RAD")
  
  CALL ChkRange_INT(hVss(iBss)%zSt,  1, hVss(iBss)%zEd, "WRONG VSS Z")
  CALL ChkRange_INT(hVss(iBss)%zEd, hVss(iBss)%zSt, nZ,  "WRONG VSS Z")
END DO

IF (hLgc%lVyg) THEN
  CALL ChkRange_INT(vAsyTyp, 1, nAsyType0, "WRONG VASYTYPE")
  CALL ChkRange_INT(vRefTyp, 1, nAsyType0, "WRONG VREFTYPE")
  CALL ChkRange_INT(vMat1,   1, nMixType,  "WRONG VMAT")
  CALL ChkRange_INT(vMat2,   1, nMixType,  "WRONG VMAT")
  CALL ChkRange_INT(vzSt,    1, vzEd,      "WRONG VZ")
  CALL ChkRange_INT(vzEd, vzSt, nZ,        "WRONG VZ")
END IF  

IF (nZ.LT.1 .AND. .NOT.nTracerCntl%lCMFD) CALL terminate("3D CALCULATION MUST USE CMFD")

NULLIFY (hCel_Loc)
NULLIFY (gCel_Loc)
NULLIFY (aInf_Loc)
! ----------------------------------------------------

END SUBROUTINE HexChkInp
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexChkCelBss()

USE geom,    ONLY : nAsyType0, nPinType, nGapPinType, nZ
USE CNTL,    ONLY : nTracerCntl
USE ioutil,  ONLY : terminate
USE HexData, ONLY : ncBss, hCelBss, ngBss, gCelBss, hCel, gCel, RodPin, GapPin
USE HexUtil, ONLY : ChkRange_INT, ChkEqual_INT, ChkEqual_REAL, ChkInc_REAL

IMPLICIT NONE

INTEGER :: iCel, iBss, jBss, iCB, jCB, iPin, iZ
REAL    :: Tmp
! ----------------------------------------------------

! CHK : Rod Cel Bss
DO iCel = 1, ncBss
  Tmp = 3 * (hCelBss(iCel)%nPin - 1) * hCelBss(iCel)%pPch &
      + 2 * hCelBss(iCel)%sRad(hCelBss(iCel)%nSub-1)
  
  CALL ChkInc_REAL(Tmp, hCelBss(iCel)%aiF2F, "RAD EXCEEDS ASY F2F INN")
END DO

! CHK : EXTRUDED
IF (nTracerCntl%AxSolver .NE. 4) RETURN

DO iPin = 1, nPinType
  iCB = hCel(RodPin(iPin)%iCel(1))%icBss
  
  DO iZ = 2, nZ
    jCB = hCel(RodPin(iPin)%iCel(iz))%icBss
    
    CALL ChkEqual_INT(iCB, jCB, 1, "NOT EXTRUDED GEOM")
  END DO
END DO

DO iPin = 1, nGapPinType
  iCB = gCel(GapPin(iPin)%iCel(1))%igBss
  
  DO iZ = 2, nZ
    jCB = gCel(GapPin(iPin)%iCel(iz))%igBss
    
    CALL ChkEqual_INT(iCB, jCB, 1, "NOT EXTRUDED GEOM")
  END DO
END DO
! ----------------------------------------------------

END SUBROUTINE HexChkCelBss
! ------------------------------------------------------------------------------------------------------------

END MODULE HexChk