MODULE HexType

USE PARAM, ONLY : TRUE, FALSE, ZERO

IMPLICIT NONE

INTEGER, PARAMETER :: nMaxFXR = 30 ! ARBITRARY

! ------------------------------------------------------------------------------------------------------------
!                                     01. HEX : Cel
! ------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------
!               01. TYPE : Hex Cel
! ----------------------------------------------------
TYPE Type_HexRodCel
  SEQUENCE
  
  LOGICAL :: luse = FALSE
  
  INTEGER :: icBss = 0  ! Numeric # of Rod Cel Basis
  INTEGER :: nPin  = 0  ! # of Pins along Asy Side
  INTEGER :: nSct  = 12 ! # of Sector Lines
  
  REAL :: pF2F  = ZERO ! Pin Flat-to-Flat Length
  REAL :: aiF2F = 0    ! Asy Inner Flat-to-Flat Length
  
  INTEGER :: nFXR = 0 ! # of FXRs
  
  REAL    :: xRad(nMaxFXR) = ZERO ! (iFXR), FXR Outer Radius
  INTEGER :: xDiv(nMaxFXR) = 0    ! (iFXR), # of Sub-rings in each FXR
  INTEGER :: xMix(nMaxFXR) = 0    ! (iFXR), Numeric # of Mixture in each FXR
  
END TYPE Type_HexRodCel
! ----------------------------------------------------
!               02. TYPE : Hex Cel Basis
! ----------------------------------------------------
TYPE Type_HexRodCelBss
  SEQUENCE
  
  !LOGICAL :: luse = TRUE ! Set only for used cel type
  
  INTEGER :: igBss = 0 ! Numeric # of Gap Cel Basis
  INTEGER :: iaTyp = 0 ! Numeric # of Gap Cel Basis
  INTEGER :: nPin  = 0 ! # of Pins along Asy Side
  INTEGER :: nSct  = 6 ! # of Sector Lines
  
  REAL :: pF2F  = ZERO ! Pin Flat-to-Flat Length
  REAL :: pPch  = ZERO ! Pin Pitch Length
  REAL :: aiF2F = ZERO ! Asy Inner Flat-to-Flat Length
  REAL :: cVol  = ZERO ! Volume of Hexagonal Cell
  
  INTEGER :: nFXR = 0 ! # of FXRs
  INTEGER :: nCel = 0 ! # of Cels
  INTEGER :: nMsh = 0 ! # of Source Meshes
  INTEGER :: nSub = 0 ! # of Sub-rings
  
  REAL    :: xRad(nMaxFXR)      = ZERO ! (iFXR), FXR Outer Radius
  INTEGER :: xDiv(nMaxFXR)      = 0    ! (iFXR), # of Sub-rings in each FXR
  REAL    :: xVol(0:2, nMaxFXR) = ZERO ! (iBndy, iFXR), FXR Volume
  
  REAL, POINTER :: sRad(:)       ! (iSub), Sub-ring Outer Radius
  REAL, POINTER :: sVol(:, :, :) ! (iTyp, iFSR, iSub), FSR Volume
  
  INTEGER, POINTER :: iCel(:) ! Numeric # of Cels
  
END TYPE Type_HexRodCelBss
! ----------------------------------------------------
!               03. TYPE : Hex Gap Cel
! ----------------------------------------------------
TYPE Type_HexGapCel
  SEQUENCE
  
  LOGICAL :: luse = FALSE
  
  INTEGER :: igBss = 0 ! Numeric # of Gap Cel Basis
  INTEGER :: nPin  = 0 ! # of Pins along Asy Side
  INTEGER :: nFXR  = 0 ! # of FXRs
  
  REAL :: pF2F  = ZERO ! Pin Flat-to-Flat Length
  REAL :: aiF2F = ZERO ! Asy Inner Flat-to-Flat Length
  
  REAL    :: xHgt(nMaxFXR) = ZERO ! (iFXR), FXR Outer Height
  INTEGER :: xDiv(nMaxFXR) = 0    ! (iFXR), # of Sub-rings in each FXR
  INTEGER :: xMix(nMaxFXR) = 0    ! (iFXR), Numeric # of Mixture in each FXR
  
END TYPE Type_HexGapCel
! ----------------------------------------------------
!               04. TYPE : Hex Gap Cel Basis
! ----------------------------------------------------
TYPE Type_HexGapCelBss
  SEQUENCE
  
  !LOGICAL :: luse = TRUE ! Set only for used cel type
  
  INTEGER :: nPin = 0 ! # of Pins along Asy Side
  INTEGER :: nSub = 0 ! # of Sub-rings
  INTEGER :: nHor = 8 ! # of Horizontal regions per Pin
  INTEGER :: nFXR = 0 ! # of FXRs
  INTEGER :: nCel = 0
  
  INTEGER :: icBss   = 0 ! Numeric # of Rod Cel Basis
  INTEGER :: nVtxHor = 0 ! # of Horizontal regions in Bndy 02
  INTEGER :: nMsh01  = 0 ! # of Source Meshes in Bndy 01
  INTEGER :: nMsh02  = 0 ! # of Source Meshes in Bndy 02
  
  REAL :: pF2F  = ZERO ! Pin Flat-to-Flat Length
  REAL :: aiF2F = ZERO ! Asy Inner Flat-to-Flat Length
  
  REAL    :: xHgt(nMaxFXR) = ZERO ! (iFXR), FXR Outer Height
  INTEGER :: xDiv(nMaxFXR) = 0    ! (iFXR), # of Sub-rings in each FXR
  
  REAL, POINTER :: sHgt(:)    ! (iSub), Sub-ring Outer Height
  REAL, POINTER :: sVol(:, :) ! (iTyp, iFSR), Source Mesh Voluem
  
  INTEGER, POINTER :: iCel(:) ! (iz), Numeric # of Cel
END TYPE Type_HexGapCelBss
! ------------------------------------------------------------------------------------------------------------
!                                     02. HEX : Pin
! ------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------
!               01. TYPE : Hex Pin
! ----------------------------------------------------
TYPE Type_HexPin
  SEQUENCE
  
  LOGICAL :: luse = FALSE
  LOGICAL :: lRod = FALSE
  LOGICAL :: lGap = FALSE
  
  INTEGER :: nFsrMax = 0
  
  INTEGER, POINTER :: iCel(:) ! (iz)
END TYPE
! ----------------------------------------------------
!               02. TYPE : Hex Pin Info
! ----------------------------------------------------
TYPE Type_HexPinInfo
  SEQUENCE
  
  INTEGER :: PinTyp = 0 ! Numeric # of "HexPin" or "GapPin"
  INTEGER :: AsyIdx = 0 ! Numeric # of "hAsy"
  INTEGER :: AsyTyp = 0 ! Numeric # of "hAsyTypInfo"
  INTEGER :: VtxTyp = 0 ! Numeric # of "hAsyTypInfo"
  
  INTEGER :: OrdInAsy01 = 0 ! Numeric # of Pin in "hAsyTypInfo"
  
  INTEGER :: ix = 0 ! x-coordinate of Pin
  INTEGER :: iy = 0 ! y-coordinate of Pin
  
  INTEGER :: nSct = 12 ! # of Sector Lines
  
  INTEGER :: FsrIdxSt = 0 ! Numeric # of 1st FSR in Core
  INTEGER :: FxrIdxSt = 0 ! Numeric # of 1st FXR in Core
  
  REAL :: Wt     = 1._8
  REAL :: Cnt(2) = ZERO ! (x/y), Origin = Asy Cnt
  
  REAL, POINTER :: Vol(:)   ! Volume of Pin in Plane
  REAL, POINTER :: VolFm(:) ! Volume of Pin in Sub-plane
  
  LOGICAL :: lInn  = TRUE
  LOGICAL :: lBndy = FALSE
  LOGICAL :: lRod  = TRUE
  LOGICAL :: lGap  = FALSE
  
END TYPE Type_HexPinInfo
! ----------------------------------------------------
!               03. TYPE : Hex CMFD Pin
! ----------------------------------------------------
TYPE Type_HexCmfdPin
  SEQUENCE
  
  REAL :: Area = ZERO
  
  ! Self
  INTEGER :: nBndy = 0
  INTEGER :: aIdx  = 0
  
  REAL :: BdPts(2, 7) = ZERO ! Origin : Core Cnt
  REAL :: BdEqn(3, 6) = ZERO ! Origin : Core Cnt
  REAL :: BdLgh   (6) = ZERO
  REAL :: BdC2B   (6) = ZERO
  
  ! MOC Pin
  INTEGER :: nmPin
  INTEGER :: mpIdx(3)
  INTEGER :: nBdmPin(6)
  INTEGER :: BdMPidx(2, 6) = 0
  INTEGER :: BdMPsuf(2, 6) = 0
  
  ! CMFD Pin
  INTEGER :: nNgh       = 0
  INTEGER :: NghPin(15) = 0 ! USE RefCell & VoidCell
  INTEGER :: NghBd (15) = 0 ! Self Bndy Idx
  INTEGER :: NghSuf(15) = 0
  REAL    :: NghLgh(15) = ZERO
  
END TYPE Type_HexCmfdPin
! ------------------------------------------------------------------------------------------------------------
!                                     03. HEX : Asy
! ------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------
!               01. TYPE : Hex Asy Type Info
! ----------------------------------------------------
TYPE Type_HexAsyTypInfo
  SEQUENCE
  
  ! ivTyp =  1 : Inn 060
  ! ivTyp =  2 : Inn 180
  ! ivTyp =  3 : Inn 360
  ! ivTyp =  4 : Bndy 01
  ! ivTyp =  5 : Bndy 02
  ! ivTyp =  6 : Bndy 02 Spt (NE)
  ! ivTyp =  7 : Bndy 02 Spt (SW), End of spTyp
  ! ivTyp =  8 : Gap 01 (NE & R)
  ! ivTyp =  9 : Gap 01 (NE & L)
  ! ivTyp = 10 : Gap 02 (Clock-wise)
  
  LOGICAL :: luse = FALSE
  
  INTEGER :: nPin = 0
  INTEGER :: iBss = 0 ! Rod Cel Bss
  INTEGER :: gTyp = 0
  
  INTEGER :: nRodPin(7) = 0 ! (iGeo)
  INTEGER :: nTotPin(7) = 0 ! (iGeo)
  
  REAL :: aiF2F = ZERO
  REAL :: aiPch = ZERO
  REAL :: pF2F  = ZERO
  REAL :: pPch  = ZERO
  REAL :: gHgt  = ZERO
  
  REAL :: mpTypPts(2, 7, 10) = ZERO ! (x/y, iBndy, ivTyp)
  REAL :: spTypPts(2, 7,  7) = ZERO ! (x/y, iBndy, ivTyp)
  
  REAL :: mpTypAre(10) = ZERO ! (ivTyp)
  REAL :: spTypAre (7) = ZERO ! (ivTyp)
  
  REAL :: mpBndyLgh(6, 10) = ZERO ! (iSuf, ivTyp)
  REAL :: mpBndyC2B(6, 10) = ZERO ! (iSuf, ivTyp)
  REAL :: spBndyLgh(6,  7) = ZERO ! (iSuf, ivTyp)
  REAL :: spBndyC2B(6,  7) = ZERO ! (iSuf, ivTyp)
  
  REAL, POINTER :: mpVtx(:, :, :, :) ! (x/y, iBndy, iGeo, iPin)
  REAL, POINTER :: spVtx(:, :, :, :) ! (x/y, iBndy, iGeo, iPin)
  
  ! MOC Pin Data
  INTEGER, POINTER :: PinLocIdx(:, :)    ! (iGeo, iPin), Numeric # in each Geo
  INTEGER, POINTER :: PinIdx(:, :)       ! (ix,   iy),   Input of Pin
  
  REAL, POINTER :: PinCnt(:, :) ! (x/y,  iPin)
  
  INTEGER, POINTER :: Pin1Dto2Dmap(:, :) ! (ix/iy, iPin), Only for rod pins
  INTEGER, POINTER :: Pin2Dto1Dmap(:, :) ! (ix, iy),      Only for rod pins
  
  LOGICAL, POINTER :: lGeoPin(:, :) ! (iGeo, iPin)
  
  REAL,    POINTER :: PinVtxAng(:, :) ! (iGeo, iPin), Rotated by (Ang)
  INTEGER, POINTER :: PinVtxTyp(:, :) ! (iGeo, iPin)
  
  ! CMFD
  INTEGER, POINTER :: cpSlfMPnum(:, :)    !       (iGeo, iPin), MOC Pin
  INTEGER, POINTER :: cpSlfMPidx(:, :, :) ! (jPin, iGeo, iPin), MOC Pin
  
  INTEGER, POINTER :: cpSufMPnum(:, :, :)    !       (iSuf, iGeo, iPin), MOC Pin
  INTEGER, POINTER :: cpSufMPidx(:, :, :, :) ! (jPin, iSuf, iGeo, iPin), MOC Pin
  INTEGER, POINTER :: cpSufMPsuf(:, :, :, :) ! (jPin, iSuf, iGeo, iPin), MOC Suf
  
END TYPE Type_HexAsyTypInfo
! ----------------------------------------------------
!               02. TYPE : Hex Asy Geo Typ Info
! ----------------------------------------------------
TYPE Type_HexGeoTypInfo
  SEQUENCE
  
  INTEGER :: nBndy = 0
  
  REAL :: Area   = ZERO
  REAL :: Cnt(2) = ZERO ! (x/y), Not Exact
  
  INTEGER :: Cor(6) = 2 ! (iBndy), Not Flat Coordinate
  
  REAL :: Vtx(2, 7) = ZERO
  REAL :: Eqn(3, 6) = ZERO
  
END TYPE Type_HexGeoTypInfo
! ----------------------------------------------------
!               03. TYPE : Hex Asy
! ----------------------------------------------------
TYPE Type_HexAsy
  SEQUENCE
  
  INTEGER :: AsyTyp    = 0
  INTEGER :: GeoTyp    = 1
  INTEGER :: iaX       = 0
  INTEGER :: iaY       = 0
  INTEGER :: nTotPin   = 0
  INTEGER :: nInnPin   = 0
  INTEGER :: nRodPin   = 0
  INTEGER :: PinIdxSt  = 0
  INTEGER :: ncBss     = 0
  INTEGER :: NghIdx(6) = 0 ! 1 = NE, 2 = EE, ... (Clock-Wise)
  INTEGER :: RotNgh    = 0 ! Ngh Asy after Rotation
  
  LOGICAL :: lHom = FALSE
  
  REAL :: Cnt(2) = ZERO ! (x/y)
  REAL :: wt     = 1._8
  
  INTEGER, POINTER :: cBss(:)
  
END TYPE Type_HexAsy
! ------------------------------------------------------------------------------------------------------------
!                                     04. HEX : Ray Geo
! ------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------
!               01. TYPE : Hex Ray Cel
! ----------------------------------------------------
TYPE Type_HexRayCel
  SEQUENCE
  
  INTEGER :: nEqn = 0
  INTEGER :: nMsh = 0
  INTEGER :: nSct = 12
  
  REAL, POINTER :: Eqn(:, :) ! (1:5, iEqn)
  
  INTEGER, POINTER :: iEqnCor(:) ! (iEqn), Not flat coordinate
  INTEGER, POINTER :: nMshEqn(:) ! (iMsh), # of Eqn
  
  INTEGER, POINTER :: MshEqnLst(:, :) ! (4, nMsh), Eqn Idx
  REAL,    POINTER :: MshEqnVal(:, :) ! (4, nMsh), (c- ax - by)
  
  INTEGER, POINTER :: MshIdx(:, :) ! (10, nMsh), for Vtx Typ
  
END TYPE Type_HexRayCel
! ----------------------------------------------------
!               02. TYPE : Hex Ray Pin Info
! ----------------------------------------------------
TYPE Type_HexRayPinInfo
  SEQUENCE
  
  INTEGER :: iCel   = 1
  REAL    :: Cnt(2) = ZERO ! (x/y)
  
  INTEGER :: nBndy(7)  = 0    ! (iGeo)
  INTEGER :: VtxTyp(7) = 0    ! (iGeo), 1 ~ 10
  REAL    :: VtxAng(7) = ZERO ! (iGeo)
  
  REAL :: Vtx(2, 7, 7) = ZERO ! (x/y, iBndy, iGeo)
  REAL :: Eqn(3, 6, 7) = ZERO ! (Val, iBndy, iGeo)
                              ! (c - ax - by) > 0 for Cnt
  
END TYPE Type_HexRayPinInfo
! ------------------------------------------------------------------------------------------------------------
!                                     05. HEX : Asy Ray Base
! ------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------
!               01. TYPE : Hex Pin Ray
! ----------------------------------------------------
TYPE Type_HexPinRay
  SEQUENCE
  
  INTEGER :: PinIdx      = 0
  INTEGER :: SufIdx(2)   = 0    ! (y¢Ö)
  REAL    :: PinPt(2, 2) = ZERO ! (x/y, y¢Ö)
  
END TYPE Type_HexPinRay
! ----------------------------------------------------
!               02. TYPE : hCel Ray
! ----------------------------------------------------
TYPE Type_HexCelRay
  SEQUENCE
  
  INTEGER :: hPinIdx    = 0
  INTEGER :: SurfIdx(2) = 0 ! y¢Ö
  INTEGER :: nSegRay    = 0
  
  REAL,    POINTER :: SegLgh(:) ! (ihsRay)
  INTEGER, POINTER :: MshIdx(:) ! (ihsRay)
  
  !REAL, POINTER :: SegPts(:, :) ! (x/y, ihsRay), for DEBUG
  
END TYPE Type_HexCelRay
! ----------------------------------------------------
!               03. TYPE : Hex Asy Ray
! ----------------------------------------------------
TYPE Type_HexAsyRay
  SEQUENCE
  
  INTEGER :: nhpRay = 0
  
  TYPE(Type_HexCelRay), POINTER :: CelRay(:) ! (ihpRay)
  
END TYPE Type_HexAsyRay
! ------------------------------------------------------------------------------------------------------------
!                                     06. HEX : Ray
! ------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------
!               01. TYPE : Hex Mod Ray
! ----------------------------------------------------
TYPE Type_HexModRay
  SEQUENCE
  
  INTEGER :: AzmIdx = 0
  
  REAL :: Eq(3)    = ZERO
  REAL :: Pt(2, 2) = ZERO ! (x/y, y¢Ö)
  
  ! iDir = Negative : y¢Ù, Positive : y¢Ö
  
  ! After Move
  INTEGER :: NxtAsy_Mov (2, -1:1) = 0 ! (x/y, iDir), Nxt Asy xy Idx
  INTEGER :: NxtMray_Mov   (-1:1) = 0 ! (iDir), Nxt mRay
  
  ! After Reflection
  INTEGER :: NxtAsy_Ref (2, -1:1, 7) = 0 ! (x/y, iDir, iGeo), Nxt Asy xy Idx
  INTEGER :: NxtMray_Ref(   -1:1, 7) = 0 ! (iDir, iGeo), Nxt mRay
  INTEGER :: NxtDir_Ref (   -1:1, 7) = 0 ! (iDir, iGeo), Nxt Dir
  
END TYPE Type_HexModRay
! ----------------------------------------------------
!               02. TYPE : Hex Core Ray
! ----------------------------------------------------
TYPE Type_HexCoreRay
  SEQUENCE
  
  INTEGER :: nmRay  = 0
  INTEGER :: AzmIdx = 0
  
  INTEGER, POINTER :: mRayIdx(:) ! (imRay)
  INTEGER, POINTER ::  AsyIdx(:) ! (imRay)
  
END TYPE Type_HexCoreRay
! ----------------------------------------------------
!               03. TYPE : Hex Rot Ray
! ----------------------------------------------------
TYPE Type_HexRotRay
  SEQUENCE
  
  INTEGER :: ncRay = 0
  
  INTEGER, POINTER :: cRayIdx(:) ! (icRay)
  
END TYPE Type_HexRotRay
! ------------------------------------------------------------------------------------------------------------
!                                     08. HEX : ETC
! ------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------
!               01. TYPE : Hex Vss
! ----------------------------------------------------
TYPE Type_HexVss
  SEQUENCE
  
  INTEGER :: Mat = 0
  INTEGER :: zSt = 0
  INTEGER :: zEd = 0
  
  REAL    :: Cnt(2) = ZERO ! (x/y)
  REAL    :: Rad(2) = ZERO ! (Inn/Out)
  
END TYPE Type_HexVss
! ----------------------------------------------------
!               02. TYPE : Hex Logical
! ----------------------------------------------------
TYPE Type_HexLogical
  SEQUENCE
  
  LOGICAL :: l060      = FALSE
  LOGICAL :: l120      = FALSE
  LOGICAL :: l360      = FALSE
  LOGICAL :: lAzmRef   = FALSE
  LOGICAL :: lAzmRot   = FALSE
  LOGICAL :: lRadRef   = FALSE
  LOGICAL :: lRadVac   = FALSE
  LOGICAL :: lAxRef(2) = FALSE ! (Bottom/Top)
  LOGICAL :: lAxVac(2) = FALSE ! (Bottom/Top)
  
  LOGICAL :: lVss     = FALSE
  LOGICAL :: lVyg     = FALSE
  LOGICAL :: lSngAsy  = FALSE
  LOGICAL :: lSngCel  = FALSE
  LOGICAL :: lspCMFD  = TRUE  ! Super-pin based CMFD
  LOGICAL :: lcrnstff = FALSE ! Corner Stiffener
  
  INTEGER :: iSym = 0 ! 1 = 60 / 2 = 120 / 3 = 360
                      ! 4 = Sng Asy / 5 = Sng Cel
  
END TYPE Type_HexLogical
! ------------------------------------------------------------------------------------------------------------

END MODULE HexType