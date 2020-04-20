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

USE HexData, ONLY : nHexPin, hPinInfo, hVss, nVss, hLgc
USE HexUtil, ONLY : FindPtLgh

IMPLICIT NONE

INTEGER :: iPin, iVss

REAL :: Lgh, Cnt(2)
! ----------------------------------------------------

IF (.NOT. hLgc%lVss) RETURN

DO iPin = 1, nHexPin
  IF (hPinInfo(iPin)%lGap) CYCLE
  
  Cnt = hPinInfo(iPin)%Cnt
  
  DO iVss = 1, nVss
    Lgh = FindPtLgh(Cnt, hVss(iVss)%Cnt)
    
    IF (Lgh < hVss(iVss)%Rad(1)) CYCLE
    IF (Lgh > hVss(iVss)%Rad(2)) CYCLE
    
    hPinInfo(iPin)%PinTyp = hVss(iVss)%vPin
  END DO
END DO
! ----------------------------------------------------

END SUBROUTINE HexSetVss
! ------------------------------------------------------------------------------------------------------------

END MODULE HexPInConst