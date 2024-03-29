MODULE HexData

USE HexType
  
IMPLICIT NONE

! Parameter
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

! Global
INTEGER :: nAsyCore   = 0 ! # of Asssemblies along Core Boundary
INTEGER :: ncBndy     = 0 ! # of Core Boundary
INTEGER :: nhAsy      = 0 ! # of Assemblies in Core
INTEGER :: ncBss      = 0 ! # of Rod Cell Basis
INTEGER :: ngBss      = 0 ! # of Gap Cell Basis
INTEGER :: ncTyp      = 0 ! # of Cell Types
INTEGER :: ncRay      = 0 ! # of Core Rays
INTEGER :: nRotRay    = 0
INTEGER :: nAzmAng    = 0
INTEGER :: nPolAng    = 0
INTEGER :: nVA        = 0 ! # of Void Asys
INTEGER :: nHexPin    = 0 ! # of MoC Pins in Core
INTEGER :: nhcPin     = 0 ! # of CMFD Pins in Core
INTEGER :: nGeoTyp    = 1
INTEGER :: nInnMOCItr = 2

! Geometry
REAL :: aoF2F = ZERO
REAL :: aoPch = ZERO

REAL :: AsyVtx (2,7) = ZERO
REAL :: AsyEqn (3,6) = ZERO
REAL :: cBndyPt(2,7) = ZERO ! Origin : Cnt Asy
REAL :: cBndyEq(3,6) = ZERO ! Origin : Cnt Asy

INTEGER, POINTER, DIMENSION(:,:) :: Asy2Dto1DMap ! (iax, iay)
INTEGER, POINTER, DIMENSION(:,:) :: Asy1Dto2DMap ! (iax/aiy, iAsy)
INTEGER, POINTER, DIMENSION(:,:) :: hCore        ! (iax, iay), Global Index of Assembly

! Vygorodka
INTEGER :: vFxr    = 0
INTEGER :: vAsyTyp = 0
INTEGER :: vRefTyp = 0
INTEGER :: vMat1   = 0
INTEGER :: vMat2   = 0
INTEGER :: vzSt    = 0
INTEGER :: vzEd    = 0

! Corner Stiffener
INTEGER :: csMat = 0

REAL :: csWdt = ZERO
REAL :: csLgh = ZERO

! Modular Ray
REAL :: Del_Inp = ZERO

INTEGER, POINTER, DIMENSION(:)   :: NumMray ! (0:iazm), 0 : # of All mRays
INTEGER, POINTER, DIMENSION(:,:) :: AngMray ! (1:2, iazm), Index of Starting & End mRay

REAL, POINTER, DIMENSION(:) :: AzmAng, AzmWgt, AzmDel, AzmDel_X, AzmDel_Y, AzmTan, AzmSin, AzmCos ! (iazm)
! ----------------------------------------------------
TYPE(Type_HexLogical) :: hLgc

TYPE(Type_HexRodCel), POINTER, DIMENSION(:) :: hCel
TYPE(Type_HexGapCel), POINTER, DIMENSION(:) :: gCel
TYPE(Type_HexGapCel), POINTER, DIMENSION(:) :: csCel

TYPE(Type_HexRodCelBss), POINTER, DIMENSION(:) :: hCelBss
TYPE(Type_HexGapCelBss), POINTER, DIMENSION(:) :: gCelBss

TYPE(Type_HexPin), POINTER, DIMENSION(:) :: RodPin
TYPE(Type_HexPin), POINTER, DIMENSION(:) :: GapPin

TYPE(Type_HexPinInfo), POINTER, DIMENSION(:) :: hPinInfo

TYPE(Type_HexGeoTypInfo), POINTER, DIMENSION(:) :: hGeoTypInfo
TYPE(Type_HexAsyTypInfo), POINTER, DIMENSION(:) :: hAsyTypInfo
TYPE(Type_HexAsy),        POINTER, DIMENSION(:) :: hAsy

TYPE(Type_HexAsyRay),  POINTER, DIMENSION(:,:,:) :: haRay   ! (iGeo, icBss, imRay)
TYPE(Type_HexModRay),  POINTER, DIMENSION(:)     :: hmRay   ! (imRay)
TYPE(Type_HexCoreRay), POINTER, DIMENSION(:)     :: hcRay   ! (icRay)
TYPE(Type_HexRotRay),  POINTER, DIMENSION(:)     :: hRotRay ! (irRay)

TYPE(Type_HexCmfdPin), POINTER, DIMENSION(:) :: hcPin
TYPE(Type_HexVss),     POINTER, DIMENSION(:) :: hVss

TYPE(Type_HexRayCel),     POINTER, DIMENSION(:) :: RayCel
TYPE(Type_HexRayPinInfo), POINTER, DIMENSION(:) :: RayPinInfo

END MODULE HexData