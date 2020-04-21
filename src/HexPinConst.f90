MODULE HexPinConst
  
IMPLICIT NONE

CONTAINS
! ------------------------------------------------------------------------------------------------------------
!                                     01. HEX SET : Pin Typ
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetPinTyp(iAsy)

USE PARAM,   ONLY : TRUE, FALSE
USE HexType, ONLY : Type_HexAsy, Type_HexAsyTypInfo, Type_HexPinInfo
USE HexData, ONLY : hAsy, hAsyTypInfo, hPinInfo, RodPin, hCel

IMPLICIT NONE

INTEGER :: iAsy, iPin, jPin, kPin, iGeo, nPin, nRod

LOGICAL :: lChk

TYPE(Type_HexAsy),        POINTER :: hAsy_Loc
TYPE(Type_HexPinInfo),    POINTER :: pInf_Loc
TYPE(Type_HexAsyTypInfo), POINTER :: aInf_Loc
! ----------------------------------------------------

hAsy_Loc => hAsy(iAsy)
aInf_Loc => hAsyTypInfo(hAsy_Loc%AsyTyp)

iGeo = hAsy_Loc%GeoTyp
! ----------------------------------------------------
!               01. SET : Rod Pin Data
! ----------------------------------------------------
DO iPin = 1, hAsy_Loc%nRodPin
  jPin = hAsy_Loc%PinIdxSt - 1 + iPin ! Numeric # in Core
  kPin = hPinInfo(jPin)%OrdInAsy01    ! Numeric # in Asy Typ
  
  pInf_Loc => hPinInfo(jPin)
  
  pInf_Loc%ix     = aInf_Loc%Pin1Dto2Dmap(1, kPin)
  pInf_Loc%iy     = aInf_Loc%Pin1Dto2Dmap(2, kPin)
  pInf_Loc%AsyTyp = hAsy_Loc%AsyTyp
  pInf_Loc%PinTyp = aInf_Loc%PinIdx(pInf_Loc%ix, pInf_Loc%iy)
  pInf_Loc%VtxTyp = aInf_Loc%PinVtxTyp(iGeo, kPin)
  pInf_Loc%nSct   = hCel(RodPin(pInf_Loc%PinTyp)%iCel(1))%nSct
  
  IF (pInf_Loc%VtxTyp < 4) CYCLE
  
  pInf_Loc%lInn  = FALSE
  pInf_Loc%lBndy = TRUE
END DO
! ----------------------------------------------------
!               02. SET : Gap Pin Data
! ----------------------------------------------------
DO iPin = hAsy_Loc%nRodPin + 1, hAsy_Loc%nTotPin
  jPin = hAsy_Loc%PinIdxSt - 1 + iPin ! Numeric # in Core
  kPin = hPinInfo(jPin)%OrdInAsy01    ! Numeric # in Asy Typ
  
  pInf_Loc => hPinInfo(jPin)
  
  pInf_Loc%AsyTyp = hAsy_Loc%AsyTyp
  pInf_Loc%PinTyp = aInf_Loc%gTyp
  pInf_Loc%VtxTyp = aInf_Loc%PinVtxTyp(iGeo, kPin)
  
  pInf_Loc%lInn = FALSE
  pInf_Loc%lRod = FALSE
  pInf_Loc%lGap = TRUE
END DO

NULLIFY (hAsy_Loc, pInf_Loc, aInf_Loc)
! ----------------------------------------------------

END SUBROUTINE HexSetPinTyp
! ------------------------------------------------------------------------------------------------------------
!                                     01. HEX SET : Pin Volume
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetPinVol(iAsy)

USE geom,    ONLY : Core, nZ, nSubPlane
USE HexType, ONLY : Type_HexPinInfo, Type_HexAsy, Type_HexAsyTypInfo
USE HexData, ONLY : hAsy, hAsyTypInfo, hPinInfo

IMPLICIT NONE

INTEGER :: iAsy, iAsyTyp, iPin, jPin, kPin, ivTyp, iz

TYPE(Type_HexPinInfo),    POINTER :: pInf_Loc
TYPE(Type_HexAsy),        POINTER :: hAsy_Loc
TYPE(Type_HexAsyTypInfo), POINTER :: aInf_Loc
! ----------------------------------------------------

hAsy_Loc => hAsy(iAsy)
aInf_Loc => hAsyTypInfo(hAsy_Loc%AsyTyp)

DO iPin = 1, hAsy_Loc%nTotPin
  jPin = hAsy_Loc%PinIdxSt - 1 + iPin ! Global Pin Idx
  
  pInf_Loc => hPinInfo(jPin)
  
  ivTyp = pInf_Loc%VtxTyp
  kPin  = pInf_Loc%OrdInAsy01
  
  DO iz = 1, nZ
    Core%PinVol  (jPin, iz) = aInf_Loc%mpTypAre(ivTyp) * Core%hz(iz)
    Core%PinVolFm(jPin, iz) = aInf_Loc%mpTypAre(ivTyp) * Core%hz(iz) / nSubPlane
    
    pInf_Loc%Vol  (iz) = Core%PinVol  (jPin, iz)
    pInf_Loc%VolFm(iz) = Core%PinVolFm(jPin, iz)
  END DO
  
  ! SET : Pin Cnt
  pInf_Loc%Cnt = aInf_Loc%PinCnt(:, kPin) + hAsy_Loc%Cnt
END DO

NULLIFY (pInf_Loc, hAsy_Loc, aInf_Loc)
! ----------------------------------------------------

END SUBROUTINE HexSetPinVol
! ------------------------------------------------------------------------------------------------------------
!                                     03. HEX SET : Pin FSR
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetPinFSR(iPin, FsrIdxSt, FxrIdxSt)

USE PARAM,   ONLY : HALF
USE geom,    ONLY : nZ, Core
USE HexType, ONLY : Type_HexPinInfo
USE HexData, ONLY : RodPin, GapPin, hPinInfo, hAsyTypInfo, hCel, hCelBss, gCel, gCelBss, hAsy

IMPLICIT NONE

INTEGER :: iPin, FsrIdxSt, FxrIdxSt
INTEGER :: iz, iFXR, igBss, ivTyp, nSub, nFSRinSub, nFXR, nFSRMax, nSubMax

TYPE(Type_HexPinInfo), POINTER :: pInf_Loc
! ----------------------------------------------------

pInf_Loc => hPinInfo(iPin)

ivTyp = pInf_Loc%VtxTyp
! ----------------------------------------------------
!               01. SET : nFXR & nFSR
! ----------------------------------------------------
nSubMax  = 0
nFSRMax  = 0
! ----------------------------
!      1. Rod Pin
! ----------------------------
IF (pInf_Loc%lRod) THEN
  SELECT CASE (ivTyp)
  CASE (1);   nFSRinSub = pInf_Loc%nSct / 6 ! ASSUME : nSct is same along z-axis
  CASE (2);   nFSRinSub = pInf_Loc%nSct / 2
  CASE (3:5); nFSRinSub = pInf_Loc%nSct
  CASE (6:7); nFSRinSub = pInf_Loc%nSct / 2
  END SELECT
  
  DO iz = 1, nZ
    nSub = hCelBss(hCel(RodPin(pInf_Loc%PinTyp)%iCel(iz))%icBss)%nSub
    nFXR = hCelBss(hCel(RodPin(pInf_Loc%PinTyp)%iCel(iz))%icBss)%nFXR
    
    nSubMax = max(nSubMax, nSub)
    nFSRMax = max(nFSRMax, nSub * nFSRinSub)
  END DO
! ----------------------------
!      2. Gap Pin
! ----------------------------
ELSE
  DO iz = 1, nZ
    igBss = gCel(GapPin(pInf_Loc%PinTyp)%iCel(iz))%igBss
    
    SELECT CASE (ivTyp)
    CASE (8:9); nFSRinSub = gCelBss(igBss)%nVtxHor
    CASE (10);  nFSRinSub = gCelBss(igBss)%nHor / 2
    END SELECT
    
    nSub = gCelBss(igBss)%nSub
    nFXR = gCelBss(igBss)%nFXR
    
    nSubMax = max(nSubMax, nSub)
    nFSRMax = max(nFSRMax, nSub * nFSRinSub)
  END DO
END IF
! ----------------------------------------------------
!               02. CP : nFXR / nFSR Max
! ----------------------------------------------------
Core%nCoreFXR = Core%nCoreFXR + nSubMax
Core%nCoreFSR = Core%nCoreFSR + nFSRMax

pInf_Loc%FSRIdxSt = FsrIdxSt
pInf_Loc%FXRIdxSt = FxrIdxSt

FsrIdxSt = FsrIdxSt + nFSRMax
FxrIdxSt = FxrIdxSt + nSubMax

RodPin(pInf_Loc%PinTyp)%nFsrMax = max(RodPin(pInf_Loc%PinTyp)%nFsrMax, nFsrMax)
! ----------------------------------------------------
!               03. SET : Wt
! ----------------------------------------------------
SELECT CASE (ivTyp)
CASE (1);   pInf_Loc%Wt = 1._8 / 6._8
CASE (2);   pInf_Loc%Wt = HALF
CASE (6:7); pInf_Loc%Wt = HALF
END SELECT

NULLIFY (pInf_Loc)
! ----------------------------------------------------

END SUBROUTINE HexSetPinFSR
! ------------------------------------------------------------------------------------------------------------
!                                     04. HEX SET : VSS
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetVss()

USE PARAM,   ONLY : FALSE, TRUE
USE geom,    ONLY : nVssTyp, nCellType, nGapType, nPinType, nGapPinType, nZ
USE HexType, ONLY : Type_HexRodCel, Type_HexGapCel, Type_HexPin
USE HexData, ONLY : nHexPin, hPinInfo, hVss, hLgc, RodPin, GapPin, hCel, gCel
USE HexUtil, ONLY : FindPtLgh

IMPLICIT NONE

INTEGER :: iPin, iVss, iTyp, iz
INTEGER :: nhc0, ngc0, nhp0, ngp0

REAL :: Lgh, Cnt(2)

LOGICAL, POINTER :: lvssCel(:, :, :)
LOGICAL, POINTER :: lvssPin(:, :, :)
INTEGER, POINTER :: aux01(:, :, :)
INTEGER, POINTER :: aux02(:, :, :)
! ----------------------------------------------------

IF (.NOT. hLgc%lVss) RETURN

ALLOCATE (lvssCel (nVssTyp, max(nCellType, nGapType),    2)); lvssCel = FALSE
ALLOCATE (lvssPin (nVssTyp, max(nPinType,  nGapPinType), 2)); lvssPin = FAlSE

ALLOCATE (aux01 (nVssTyp, max(nCellType, nGapType),    2)); aux01 = 0
ALLOCATE (aux02 (nVssTyp, max(nPinType,  nGapPinType), 2)); aux02 = 0

nhc0 = nCellType
ngc0 = nGapType
nhp0 = nPinType
ngp0 = nGapPinType
! ----------------------------------------------------
!               01. SEARCH : Vss Cel
! ----------------------------------------------------
DO iPin = 1, nHexPin
  Cnt = hPinInfo(iPin)%Cnt
  
  DO iVss = 1, nVssTyp
    Lgh = FindPtLgh(Cnt, hVss(iVss)%Cnt)
    
    IF (Lgh < hVss(iVss)%Rad(1)) CYCLE
    IF (Lgh > hVss(iVss)%Rad(2)) CYCLE
    
    iTyp = hPinInfo(iPin)%PinTyp
    
    IF (hPinInfo(iPin)%lRod) THEN
      lvssPin(iVss, iTyp, 1) = TRUE
      
      DO iz = 1, hVss(iVss)%zSt, hVss(iVss)%zEd
        lvssCel(iVss, RodPin(iTyp)%iCel(iz), 1) = TRUE
      END DO
    ELSE
      lvssPin(iVss, iTyp, 2) = TRUE
      
      DO iz = 1, hVss(iVss)%zSt, hVss(iVss)%zEd
        lvssCel(iVss, GapPin(iTyp)%iCel(iz), 2) = TRUE
      END DO
    END IF
  END DO
END DO
! ----------------------------------------------------
!               02. SET : Vss Cel
! ----------------------------------------------------
DO iTyp = 1, nhc0
  DO iVss = 1, nVssTyp
    IF (.NOT. lvssCel(iVss, iTyp, 1)) CYCLE
    
    nCellType = nCellType + 1
    
    hCel(nCellType)%luse  = TRUE
    hCel(nCellType)%icBss = hCel(iTyp)%icBss
    hCel(nCellType)%nPin  = hCel(iTyp)%nPin
    hCel(nCellType)%nSct  = hCel(iTyp)%nSct
    hCel(nCellType)%pF2F  = hCel(iTyp)%pF2F
    hCel(nCellType)%aiF2F = hCel(iTyp)%aiF2F
    hCel(nCellType)%nFXR  = hCel(iTyp)%nFXR
    
    hCel(nCellType)%xRad(1:hCel(iTyp)%nFXR) = hCel(iTyp)%xRad(1:hCel(iTyp)%nFXR)
    hCel(nCellType)%xDiv(1:hCel(iTyp)%nFXR) = hCel(iTyp)%xDiv(1:hCel(iTyp)%nFXR)
    
    hCel(nCellType)%xMix(1:hCel(iTyp)%nFXR) = hVss(iVss)%vMat
    
    aux01(iVss, iTyp, 1) = nCellType
  END DO
END DO

DO iTyp = 1, ngc0
  DO iVss = 1, nVssTyp
    IF (.NOT. lvssCel(iVss, iTyp, 2)) CYCLE
    
    nGapType = nGapType + 1
    
    gCel(nGapType)%luse  = TRUE
    gCel(nGapType)%igBss = gCel(iTyp)%igBss
    gCel(nGapType)%nPin  = gCel(iTyp)%nPin
    gCel(nGapType)%pF2F  = gCel(iTyp)%pF2F
    gCel(nGapType)%aiF2F = gCel(iTyp)%aiF2F
    gCel(nGapType)%nFXR  = gCel(iTyp)%nFXR
    
    gCel(nGapType)%xHgt(1:gCel(iTyp)%nFXR) = gCel(iTyp)%xHgt(1:gCel(iTyp)%nFXR)
    gCel(nGapType)%xDiv(1:gCel(iTyp)%nFXR) = gCel(iTyp)%xDiv(1:gCel(iTyp)%nFXR)
        
    gCel(nGapType)%xMix(1:gCel(iTyp)%nFXR) = hVss(iVss)%vMat
    
    aux01(iVss, iTyp, 2) = nGapType
  END DO
END DO
! ----------------------------------------------------
!               03. SET : Vss Pin
! ----------------------------------------------------
DO iTyp = 1, nhp0
  DO iVss = 1, nVssTyp
    IF (.NOT. lvssPin(iVss, iTyp, 1)) CYCLE
    
    nPinType = nPinType + 1
    
    RodPin(nPinType)%luse    = TRUE
    RodPin(nPinType)%lRod    = TRUE
    RodPin(nPinType)%lGap    = FALSE
    RodPin(nPinType)%nFsrMax = RodPin(iTyp)%nFsrMax
    
    ALLOCATE (RodPin(nPinType)%iCel (nZ))
    
    RodPin(nPinType)%iCel(1:nZ) = RodPin(iTyp)%iCel(1:nZ)
    
    DO iz = hVss(iVss)%zSt, hVss(iVss)%zEd
      RodPin(nPinType)%iCel(iz) = aux01(iVss, RodPin(nPinType)%iCel(iz), 1)
    END DO
    
    aux02(iVss, iTyp, 1) = nPinType
  END DO
END DO

DO iTyp = 1, ngp0
  DO iVss = 1, nVssTyp
    IF (.NOT. lvssPin(iVss, iTyp, 2)) CYCLE
    
    nGapPinType = nGapPinType + 1
    
    GapPin(nGapPinType)%luse    = TRUE
    GapPin(nGapPinType)%lRod    = FALSE
    GapPin(nGapPinType)%lGap    = TRUE
    GapPin(nGapPinType)%nFsrMax = GapPin(iTyp)%nFsrMax
    
    ALLOCATE (GapPin(nGapPinType)%iCel (nZ))
    
    GapPin(nGapPinType)%iCel(1:nZ) = GapPin(iTyp)%iCel(1:nZ)
    
    DO iz = hVss(iVss)%zSt, hVss(iVss)%zEd
      GapPin(nGapPinType)%iCel(iz) = aux01(iVss, GapPin(nGapPinType)%iCel(iz), 2)
    END DO
    
    aux02(iVss, iTyp, 2) = nGapPinType
  END DO
END DO
! ----------------------------------------------------
!               04. CP : Vss Cel
! ----------------------------------------------------
DO iPin = 1, nHexPin
  Cnt = hPinInfo(iPin)%Cnt
  
  DO iVss = 1, nVssTyp
    Lgh = FindPtLgh(Cnt, hVss(iVss)%Cnt)
    
    IF (Lgh < hVss(iVss)%Rad(1)) CYCLE
    IF (Lgh > hVss(iVss)%Rad(2)) CYCLE
    
    iTyp = hPinInfo(iPin)%PinTyp
    
    IF (hPinInfo(iPin)%lRod) THEN
      hPinInfo(iPin)%PinTyp = aux02(iVss, iTyp, 1)
    ELSE
      hPinInfo(iPin)%PinTyp = aux02(iVss, iTyp, 2)
    END IF
  END DO
END DO

DEALLOCATE (lvssCel, lvssPin, aux01, aux02)
! ----------------------------------------------------

END SUBROUTINE HexSetVss
! ------------------------------------------------------------------------------------------------------------

END MODULE HexPinConst