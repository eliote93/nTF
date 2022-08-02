MODULE HexType

USE PARAM, ONLY : TRUE, FALSE, ZERO

IMPLICIT NONE

INTEGER, PARAMETER :: nMaxFXR = 30 ! ARBITRARY

! ------------------------------------------------------------------------------------------------------------
!                                     01. Cel
! ------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------
!               01. Rod Cel
! ----------------------------------------------------
TYPE Type_HexRodCel
  
  LOGICAL :: luse = FALSE
  
  INTEGER :: icBss = 0  ! Index of Rod Cel Basis
  INTEGER :: nPin  = 0  ! # of Pins along Asy Side
  INTEGER :: nSct  = 12 ! # of Sector Lines
  
  REAL :: pF2F  = ZERO ! Pin Flat-to-Flat Length
  REAL :: aiF2F = 0    ! Asy Inner Flat-to-Flat Length
  
  INTEGER :: nFXR = 0 ! # of FXRs
  
  REAL    :: xRad(nMaxFXR) = ZERO ! (iFXR), FXR Outer Radius
  INTEGER :: xDiv(nMaxFXR) = 0    ! (iFXR), # of Sub-rings in each FXR
  INTEGER :: xMix(nMaxFXR) = 0    ! (iFXR), Index of Mixture in each FXR
  
END TYPE Type_HexRodCel
! ----------------------------------------------------
!               02. Rod Cel Basis
! ----------------------------------------------------
TYPE Type_HexRodCelBss
  
  !LOGICAL :: luse = TRUE ! Set only for used cel type
  
  INTEGER :: igBss = 0 ! Index of Gap Cel Basis
  INTEGER :: iaTyp = 0 ! Index of Gap Cel Basis
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
  
  REAL, POINTER, DIMENSION(:)     :: sRad ! (iSub), Sub-ring Outer Radius
  REAL, POINTER, DIMENSION(:,:,:) :: sVol ! (iTyp, iFSR, iSub), FSR Volume
  
  INTEGER, POINTER, DIMENSION(:) :: iCel ! Index of Cell
  
END TYPE Type_HexRodCelBss
! ----------------------------------------------------
!               03. Gap Cel
! ----------------------------------------------------
TYPE Type_HexGapCel
  
  LOGICAL :: luse = FALSE
  
  INTEGER :: igBss = 0 ! Index of Gap Cel Basis
  INTEGER :: nPin  = 0 ! # of Pins along Asy Side
  INTEGER :: nFXR  = 0 ! # of FXRs
  
  REAL :: pF2F  = ZERO ! Pin Flat-to-Flat Length
  REAL :: aiF2F = ZERO ! Asy Inner Flat-to-Flat Length
  
  REAL    :: xHgt(nMaxFXR) = ZERO ! (iFXR), FXR Outer Height
  INTEGER :: xDiv(nMaxFXR) = 0    ! (iFXR), # of Sub-rings in each FXR
  INTEGER :: xMix(nMaxFXR) = 0    ! (iFXR), Index of Mixture in each FXR
  
END TYPE Type_HexGapCel
! ----------------------------------------------------
!               04. Gap Cel Basis
! ----------------------------------------------------
TYPE Type_HexGapCelBss
  
  !LOGICAL :: luse = TRUE ! Set only for used cel type
  
  INTEGER :: nPin = 0 ! # of Pins along Asy Side
  INTEGER :: nSub = 0 ! # of Sub-rings
  INTEGER :: nHor = 8 ! # of Horizontal regions per Pin
  INTEGER :: nFXR = 0 ! # of FXRs
  INTEGER :: nCel = 0
  
  INTEGER :: icBss   = 0 ! Index of Rod Cel Basis
  INTEGER :: nVtxHor = 0 ! # of Horizontal regions in Bndy 02
  INTEGER :: nMsh01  = 0 ! # of Source Meshes in Bndy 01
  INTEGER :: nMsh02  = 0 ! # of Source Meshes in Bndy 02
  
  REAL :: pF2F  = ZERO ! Pin Flat-to-Flat Length
  REAL :: aiF2F = ZERO ! Asy Inner Flat-to-Flat Length
  
  REAL    :: xHgt(nMaxFXR) = ZERO ! (iFXR), FXR Outer Height
  INTEGER :: xDiv(nMaxFXR) = 0    ! (iFXR), # of Sub-rings in each FXR
  
  REAL, POINTER, DIMENSION(:)   :: sHgt ! (iSub), Sub-ring Outer Height
  REAL, POINTER, DIMENSION(:,:) :: sVol ! (iTyp, iFSR), Source Mesh Voluem
  
  INTEGER, POINTER, DIMENSION(:) :: iCel ! (iz), Index of Cell
  
END TYPE Type_HexGapCelBss
! ------------------------------------------------------------------------------------------------------------
!                                     02. Pin
! ------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------
!               01. Pin
! ----------------------------------------------------
TYPE Type_HexPin
  
  LOGICAL :: luse = FALSE
  LOGICAL :: lRod = FALSE
  LOGICAL :: lGap = FALSE
  
  INTEGER :: nFsrMax = 0
  
  INTEGER, POINTER, DIMENSION(:) :: iCel ! (iz)
  
END TYPE
! ----------------------------------------------------
!               02. Pin Info
! ----------------------------------------------------
TYPE Type_HexPinInfo
  
  INTEGER :: PinTyp = 0 ! Index of "HexPin" or "GapPin"
  INTEGER :: AsyIdx = 0 ! Global Index of "hAsy"
  INTEGER :: AsyTyp = 0 ! Index of "hAsyTypInfo"
  INTEGER :: VtxTyp = 0 ! Geometric Typ (1 ~ 10)
  INTEGER :: ihcPin = 0 ! Global Index of "hcPin"
  
  INTEGER :: OrdInAsy01 = 0 ! Local Index of Pin in "hAsyTypInfo"
  
  INTEGER :: ix = 0 ! x-coordinate of Pin
  INTEGER :: iy = 0 ! y-coordinate of Pin
  
  INTEGER :: nSct = 12 ! # of Sector Lines
  
  INTEGER :: FsrIdxSt = 0 ! Global Index of 1st FSR
  INTEGER :: FxrIdxSt = 0 ! Global Index of 1st FXR
  
  INTEGER :: DcmpMP2slfSPngh(6) = 0 ! Neighboring Index of Self CMFD Pin which MoC Pin belongs to
  INTEGER :: DcmpMP2nghSPidx(6) = 0 ! Global Index of Neighboring CMFD Pin with CMFD Pin which MoC Pin belongs to
  
  REAL :: Wt     = 1._8
  REAL :: Cnt(2) = ZERO ! (x/y), Origin = Asy Cnt
  
  REAL, POINTER, DIMENSION(:) :: Vol   ! Volume of Pin in Plane
  REAL, POINTER, DIMENSION(:) :: VolFm ! Volume of Pin in Sub-plane
  
  LOGICAL :: lInn  = TRUE
  LOGICAL :: lBndy = FALSE
  LOGICAL :: lRod  = TRUE
  LOGICAL :: lGap  = FALSE
  
END TYPE Type_HexPinInfo
! ----------------------------------------------------
!               03. CMFD Pin
! ----------------------------------------------------
TYPE Type_HexCmfdPin
  
  REAL :: Area = ZERO
  
  ! Self
  INTEGER :: nBndy = 0
  INTEGER :: aIdx  = 0
  
  REAL :: BdPts(2, 7) = ZERO ! Origin : Core Cnt
  REAL :: BdEqn(3, 6) = ZERO ! Origin : Core Cnt
  REAL :: BdLgh   (6) = ZERO
  REAL :: BdC2B   (6) = ZERO
  
  ! MoC Pin
  INTEGER :: nmPin             ! # of MoC Pins
  INTEGER :: mpIdx(3)          ! Global Index of MoC Pin
  INTEGER :: nBdmPin(6)        ! # of MoC Pins along Boundary
  INTEGER :: BdMPidx(2, 6) = 0 ! Global Index of MoC Pin at Boundary
  INTEGER :: BdMPsuf(2, 6) = 0 ! Boundary Index of MoC Pin at Boundary
  
  ! CMFD Pin
  INTEGER :: nNgh       = 0 ! # of Neighboring CMFD Pins
  INTEGER :: NghPin(15) = 0 ! Global Index of CMFD Pin using RefCell & VoidCell
  INTEGER :: NghBd (15) = 0 ! Self Boundary Index
  INTEGER :: NghSuf(15) = 0 ! Boundary Index of Neighboring CMFD Pin
  REAL    :: NghLgh(15) = ZERO
  
END TYPE Type_HexCmfdPin
! ------------------------------------------------------------------------------------------------------------
!                                     03. Asy
! ------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------
!               01. Asy Type Info
! ----------------------------------------------------
TYPE Type_HexAsyTypInfo
  
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
  
  INTEGER :: nPin   = 0
  INTEGER :: iBss   = 0 ! Rod Cel Bss
  INTEGER :: gTyp   = 0
  INTEGER :: cstTyp = 0
  INTEGER :: cstNum = 0
  
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
  
  REAL, POINTER, DIMENSION(:,:,:,:) :: mpVtx ! (x/y, iBndy, iGeo, iPin)
  REAL, POINTER, DIMENSION(:,:,:,:) :: spVtx ! (x/y, iBndy, iGeo, iPin)
  
  ! MOC Pin Data
  INTEGER, POINTER, DIMENSION(:,:) :: PinLocIdx ! (iGeo, iPin), Index in each Geo
  INTEGER, POINTER, DIMENSION(:,:) :: PinIdx    ! (ix,   iy),   Input of Pin
  
  REAL, POINTER, DIMENSION(:,:) :: PinCnt ! (x/y,  iPin)
  
  INTEGER, POINTER, DIMENSION(:,:) :: Pin1Dto2Dmap ! (ix/iy, iPin), Only for rod pins
  INTEGER, POINTER, DIMENSION(:,:) :: Pin2Dto1Dmap ! (ix, iy),      Only for rod pins
  
  LOGICAL, POINTER, DIMENSION(:,:) :: lGeoPin ! (iGeo, iPin)
  
  REAL,    POINTER, DIMENSION(:,:) :: PinVtxAng ! (iGeo, iPin), Rotated by (Ang)
  INTEGER, POINTER, DIMENSION(:,:) :: PinVtxTyp ! (iGeo, iPin)
  
  INTEGER, POINTER, DIMENSION(:) :: CstMap ! (iPin), Corner Stiffener Map
  
  ! CMFD
  INTEGER, POINTER, DIMENSION(:,:)   :: cnpSlfMPnum !       (iGeo, iPin), MOC Pin
  INTEGER, POINTER, DIMENSION(:,:,:) :: cnpSlfMPidx ! (jPin, iGeo, iPin), MOC Pin
  
  INTEGER, POINTER, DIMENSION(:,:,:)   :: cnpSufMPnum !       (iSuf, iGeo, iPin), MOC Pin
  INTEGER, POINTER, DIMENSION(:,:,:,:) :: cnpSufMPidx ! (jPin, iSuf, iGeo, iPin), MOC Pin
  INTEGER, POINTER, DIMENSION(:,:,:,:) :: cnpSufMPsuf ! (jPin, iSuf, iGeo, iPin), MOC Suf
  
END TYPE Type_HexAsyTypInfo
! ----------------------------------------------------
!               02. Asy Geo Typ Info
! ----------------------------------------------------
TYPE Type_HexGeoTypInfo
  
  INTEGER :: nBndy = 0
  
  REAL :: Area   = ZERO
  REAL :: Cnt(2) = ZERO ! (x/y), Not Exact
  
  INTEGER :: Cor(6) = 2 ! (iBndy), Not Flat Coordinate
  
  REAL :: Vtx(2, 7) = ZERO
  REAL :: Eqn(3, 6) = ZERO
  
END TYPE Type_HexGeoTypInfo
! ----------------------------------------------------
!               03. Hex Asy
! ----------------------------------------------------
TYPE Type_HexAsy
  
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
  
  INTEGER, POINTER, DIMENSION(:) :: cBss
  
END TYPE Type_HexAsy
! ------------------------------------------------------------------------------------------------------------
!                                     04. Ray Geo
! ------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------
!               01. Ray Cel
! ----------------------------------------------------
TYPE Type_HexRayCel
  
  INTEGER :: nEqn = 0
  INTEGER :: nMsh = 0
  INTEGER :: nSct = 12
  
  REAL, POINTER, DIMENSION(:,:) :: Eqn ! (1:5, iEqn)
  
  INTEGER, POINTER, DIMENSION(:) :: iEqnCor ! (iEqn), Not flat coordinate
  INTEGER, POINTER, DIMENSION(:) :: nMshEqn ! (iMsh), # of Eqn
  
  INTEGER, POINTER, DIMENSION(:,:) :: MshEqnLst ! (4, nMsh), Eqn Idx
  REAL,    POINTER, DIMENSION(:,:) :: MshEqnVal ! (4, nMsh), (c- ax - by)
  
  INTEGER, POINTER, DIMENSION(:,:) :: MshIdx ! (10, nMsh), for Vtx Typ
  
END TYPE Type_HexRayCel
! ----------------------------------------------------
!               02. Ray Pin Info
! ----------------------------------------------------
TYPE Type_HexRayPinInfo
  
  INTEGER :: iCel   = 1
  REAL    :: Cnt(2) = ZERO ! (x/y)
  
  INTEGER :: nBndy(7)  = 0    ! (iGeo)
  INTEGER :: VtxTyp(7) = 0    ! (iGeo), 1 ~ 10
  REAL    :: VtxAng(7) = ZERO ! (iGeo)
  
  REAL :: Vtx(2, 7, 7) = ZERO ! (x/y, iBndy, iGeo)
  REAL :: Eqn(3, 6, 7) = ZERO ! (Val, iBndy, iGeo), (c - ax - by) > 0 for Cnt
  
END TYPE Type_HexRayPinInfo
! ------------------------------------------------------------------------------------------------------------
!                                     05. Asy Ray Base
! ------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------
!               01. Pin Ray
! ----------------------------------------------------
TYPE Type_HexPinRay ! Temporary
  
  INTEGER :: PinIdx      = 0
  INTEGER :: SurfIdx(2)  = 0    ! (y¢Ö)
  REAL    :: PinPt(2, 2) = ZERO ! (x/y, y¢Ö)
  
  REAL :: hsn(2) = ZERO ! sinv, Memory Can be Reduced by Using Angle Index
  REAL :: hcs(2) = ZERO ! cosv
  
END TYPE Type_HexPinRay
! ----------------------------------------------------
!               02. Cel Ray
! ----------------------------------------------------
TYPE Type_HexCelRay
  
  INTEGER :: hPinIdx    = 0
  INTEGER :: hSufIdx(2) = 0 ! y¢Ö
  INTEGER :: nSegRay    = 0
  
  REAL,    POINTER, DIMENSION(:) :: SegLgh ! (ihsRay)
  INTEGER, POINTER, DIMENSION(:) :: MshIdx ! (ihsRay)
  
  REAL :: hsn(2) = ZERO ! sinv
  REAL :: hcs(2) = ZERO ! cosv
  
  !REAL, POINTER, DIMENSION(:,:) :: SegPts ! (x/y, ihsRay), for DEBUG
  
END TYPE Type_HexCelRay
! ----------------------------------------------------
!               03. Asy Ray
! ----------------------------------------------------
TYPE Type_HexAsyRay
  
  INTEGER :: nhpRay = 0
  
  TYPE(Type_HexCelRay), POINTER, DIMENSION(:) :: CelRay ! (ihpRay)
  
END TYPE Type_HexAsyRay
! ------------------------------------------------------------------------------------------------------------
!                                     06. Ray
! ------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------
!               01. Mod Ray
! ----------------------------------------------------
TYPE Type_HexModRay
  
  INTEGER :: AzmIdx = 0
  
  REAL :: Eq(3)    = ZERO ! ax * by = c
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
!               02. Core Ray
! ----------------------------------------------------
TYPE Type_HexCoreRay
  
  INTEGER :: nmRay  = 0
  INTEGER :: AzmIdx = 0
  
  INTEGER, POINTER, DIMENSION(:) :: mRayIdx ! (imRay)
  INTEGER, POINTER, DIMENSION(:) ::  AsyIdx ! (imRay)
  
END TYPE Type_HexCoreRay
! ----------------------------------------------------
!               03. Rot Ray
! ----------------------------------------------------
TYPE Type_HexRotRay
  
  INTEGER :: ncRay = 0
  
  INTEGER, POINTER, DIMENSION(:) :: cRayIdx ! (icRay)
  
END TYPE Type_HexRotRay
! ------------------------------------------------------------------------------------------------------------
!                                     08. ETC
! ------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------
!               01. Vss.
! ----------------------------------------------------
TYPE Type_HexVss
  
  INTEGER :: Mat = 0
  INTEGER :: zSt = 0
  INTEGER :: zEd = 0
  
  REAL :: Cnt(2) = ZERO ! (x/y)
  REAL :: Rad(2) = ZERO ! (Inn/Out)
  
END TYPE Type_HexVss
! ----------------------------------------------------
!               02. Logical
! ----------------------------------------------------
TYPE Type_HexLogical
  
  LOGICAL :: l060      = FALSE
  LOGICAL :: l120      = FALSE
  LOGICAL :: l360      = FALSE
  LOGICAL :: lAzmRef   = FALSE
  LOGICAL :: lAzmRot   = FALSE
  LOGICAL :: lRadRef   = FALSE
  LOGICAL :: lRadVac   = FALSE
  LOGICAL :: lAxRef(2) = FALSE ! (Bottom/Top)
  LOGICAL :: lAxVac(2) = FALSE ! (Bottom/Top)
  LOGICAL :: lNoRef    = FALSE ! 360 & Vac & .NOT.SngAsy
  
  LOGICAL :: lVss     = FALSE
  LOGICAL :: lVyg     = FALSE
  LOGICAL :: lSngAsy  = FALSE
  LOGICAL :: lSngCel  = FALSE
  LOGICAL :: lspCMFD  = TRUE  ! Super-pin based CMFD
  
  INTEGER :: iSym = 0 ! 1 = 60 / 2 = 120 / 3 = 360 / 4 = Sng Asy / 5 = Sng Cel
  INTEGER :: idcmpclr
  LOGICAL :: ldcmpad
  
END TYPE Type_HexLogical
! ------------------------------------------------------------------------------------------------------------

END MODULE HexType