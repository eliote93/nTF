MODULE HexData

USE HexType
  
IMPLICIT NONE

REAL, PARAMETER :: hEps   = 1E-7
REAL, PARAMETER :: Sq3    = 1.73205080756888_8
REAL, PARAMETER :: Sq3Inv = 0.577350269189626_8
REAL, PARAMETER :: PI_2   = 1.5707963267948966_8
REAL, PARAMETER :: PI_3   = 1.04719755119660_8
REAL, PARAMETER :: PI_6   = 0.523598775598300_8

INTEGER, PARAMETER :: mpTypNumNgh(10) = [3, 4, 6,    & ! Inn
                                         5, 6, 4, 4, & ! Bndy
                                         4, 4, 4]      ! Gap
INTEGER, PARAMETER :: spTypNumNgh(7)  = [3, 4, 6,  & ! Inn
                                         5, 5, 4, 4] ! Bndy

INTEGER :: nAsyCore   = 0
INTEGER :: ncBndy     = 0
INTEGER :: nhAsy      = 0
INTEGER :: ncBss      = 0
INTEGER :: ngBss      = 0
INTEGER :: nInf       = 0
INTEGER :: ncRay      = 0
INTEGER :: nRotRay    = 0
INTEGER :: nAzmAng    = 0
INTEGER :: nPolAng    = 0
INTEGER :: nVss       = 0
INTEGER :: VygFxrSt   = 0
INTEGER :: nVA        = 0
INTEGER :: VygAsyTyp  = 0
INTEGER :: VygMat     = 0
INTEGER :: nHexPin    = 0
INTEGER :: nhcPin     = 0
INTEGER :: nGeoTyp    = 1
INTEGER :: nInnMOCItr = 2

REAL :: aoF2F   = ZERO
REAL :: aoPch   = ZERO
REAL :: Del_Inp = ZERO

REAL :: AsyVtx(2, 7)  = ZERO
REAL :: AsyEqn(3, 6)  = ZERO
REAL :: cBndyPt(2, 7) = ZERO ! Origin : Cnt Asy
REAL :: cBndyEq(3, 6) = ZERO ! Origin : Cnt Asy

REAL, POINTER :: AzmAng(:), AzmWgt(:)
REAL, POINTER :: AzmDel(:), AzmDel_X(:), AzmDel_Y(:)
REAL, POINTER :: AzmTan(:), AzmSin(:), AzmCos(:)

REAL, POINTER :: AxJout (:, :, :, :)

INTEGER, POINTER :: NumMray(:)
INTEGER, POINTER :: AngMray(:, :)

INTEGER, POINTER :: Asy2Dto1DMap(:, :) ! (ix, iy)
INTEGER, POINTER :: Asy1Dto2DMap(:, :) ! (ix/iy, iAsy)

TYPE(Type_HexPin), POINTER :: RodPin(:)
TYPE(Type_HexPin), POINTER :: GapPin(:)

INTEGER, POINTER :: hCore (:, :)

REAL,    POINTER :: VssRad(:, :)
INTEGER, POINTER :: va2D(:, :)
! ----------------------------------------------------
TYPE(Type_HexLogical) :: hLgc

TYPE(Type_HexRodCel), POINTER :: hCel(:)
TYPE(Type_HexGapCel), POINTER :: gCel(:)

TYPE(Type_HexRodCelBss), POINTER :: hCelBss(:)
TYPE(Type_HexGapCelBss), POINTER :: gCelBss(:)

TYPE(Type_HexPinInfo), POINTER :: hPinInfo(:)

TYPE(Type_HexGeoTypInfo), POINTER :: hGeoTypInfo(:)
TYPE(Type_HexAsyTypInfo), POINTER :: hAsyTypInfo(:)
TYPE(Type_HexAsy),        POINTER :: hAsy(:)

TYPE(Type_HexAsyRay),  POINTER :: haRay(:,:,:) ! (iGeo, icBss, imRay)
TYPE(Type_HexModRay),  POINTER :: hmRay(:)     ! (imRay)
TYPE(Type_HexCoreRay), POINTER :: hcRay(:)     ! (icRay)
TYPE(Type_HexRotRay),  POINTER :: hRotRay(:)   ! (irRay)

TYPE(Type_HexCmfdPin), POINTER :: hcPin(:)
TYPE(Type_HexVss),     POINTER :: hVss(:)

TYPE(Type_HexRayCel),     POINTER :: RayCel(:)
TYPE(Type_HexRayPinInfo), POINTER :: RayPinInfo(:)

END MODULE HexData
! ------------------------------------------------------------------------------------------------------------
!                                     01. 
! ------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------
!               01. 
! ----------------------------------------------------
! ----------------------------
!      1. 
! ----------------------------