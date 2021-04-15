MODULE HexCnP ! Copy & Paste

IMPLICIT NONE

INTEGER, PARAMETER :: pnRod = 6
INTEGER, PARAMETER :: pnGap = 2

! 1 :  Rod Cel in  60-degree
! 2 :  Rod Cel in 120-degree
! 3 :  Rod Cel in 360-degree
! 4 : Bndy Cel of Type 1
! 5 : Bndy Cel of Type 2
! 6 : Bndy Cel of Type 2, Spt
! 7 :  Gap Cel of Type 1
! 8 :  Gap Cel of Type 2

CONTAINS

! ------------------------------------------------------------------------------------------------------------
!                                     01. HEX CnP : Geo
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexCnPgeo()

USE geom,    ONLY : Core, Asy, AsyInfo, AsyVol, PinInfo, CellInfo, Pin, nCellType, nPinType
USE HexData, ONLY : nHexPin, RodPin

IMPLICIT NONE

INTEGER :: iPin, ipTyp, iAsy
! ----------------------------------------------------

CALL HexCnPCelinfo

DO iPin = 1, nHexPin
  CALL HexCnPPin(iPin)
END DO

ALLOCATE (PinInfo (nPinType))

DO ipTyp = 1, nPinType
  IF (.NOT. RodPin(ipTyp)%luse) CYCLE
  
  PinInfo(ipTyp)%nFsrMax = RodPin(ipTyp)%nFsrMax
END DO

Core%nxyc = 0

DO iAsy = 1, Core%nAsyType
  Core%nxyc = max(Core%nxyc, AsyInfo(iAsy)%nxy)
END DO

Core%Asy      => Asy
Core%AsyInfo  => AsyInfo
Core%AsyVol   => AsyVol
Core%Pin      => Pin
Core%PinInfo  => PinInfo
Core%CellInfo => CellInfo
Core%nCellType = nCellType
! ----------------------------------------------------

END SUBROUTINE HexCnPgeo
! ------------------------------------------------------------------------------------------------------------
!                                     02. HEX CnP : Cel Info
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexCnPCelinfo()

USE allocs
USE TYPEDEF,      ONLY : cell_type
USE param,        ONLY : TRUE, FALSE, HALF, ZERO
USE GEOM,         ONLY : CellInfo, nCellType, nGapType
USE cntl,         ONLY : nTracerCntl
USE Material_Mod, ONLY : Mixture
USE BenchXs,      ONLY : MacXsBen

USE HexType, ONLY : Type_HexRodCelBss, Type_HexGapCelBss
USE HexData, ONLY : ncTyp, hCel, gCel, hCelBss, gCelBss, hLgc

IMPLICIT NONE

INTEGER :: iCel, jCel, iDir, iSub, jSub, iFXR, iFSR, jFSR, iReg
INTEGER :: nSct, nSub, nFXR, nFSR, nFSRinSub

INTEGER :: rnFSRinSub(1:pnRod) = [6, 2, 1, 1, 1, 2]
INTEGER ::    BndyTyp(1:pnRod) = [0, 0, 0, 1, 2, 2]

LOGICAL :: lSkip(6) = FALSE

TYPE(Cell_type),         POINTER :: Cel_Loc
TYPE(Type_HexRodCelBss), POINTER :: cBs_Loc
TYPE(Type_HexGapCelBss), POINTER :: gBs_Loc
! ----------------------------------------------------

ncTyp = pnRod * nCellType + pnGap * nGapType

IF (hLgc%lSngCel) THEN
  ncTyp    = pnRod * nCellType
  lSkip    = TRUE
  lSkip(3) = FALSE
END IF

ALLOCATE (CellInfo (ncTyp))

DO iCel = 1, ncTyp
  CellInfo(iCel)%luse = FALSE
END DO
! ----------------------------------------------------
!               01. Rod Cel
! ----------------------------------------------------
DO iCel = 1, nCellType
  IF (.NOT. hCel(iCel)%luse) CYCLE
  
  cBs_Loc => hCelBss(hCel(iCel)%icBss)
  
  nSct = cBs_Loc%nSct
  nSub = cBs_Loc%nSub
  nFXR = cBs_Loc%nFXR
  ! --------------------------------------------------
  DO iDir = 1, pnRod
    IF (lSkip(iDir)) CYCLE
    
    jCel = (iCel - 1) * pnRod + iDir
    
    Cel_Loc => CellInfo(jCel)
    
    nFSRinSub = nSct / rnFSRinSub(iDir)
    
    CALL dmalloc(Cel_Loc%nFSRinFXR, nSub)
    CALL dmalloc0(Cel_Loc%Vol, 0, nSub- 1) ! Need to be changed into FSR unit
    CALL dmalloc(Cel_Loc%iReg, nSub * nFSRinSub)
    CALL dmalloc(Cel_Loc%MapFxr2FsrIdx, nFSRinSub, nSub)
    
    Cel_Loc%luse   = TRUE
    Cel_Loc%lfuel  = FALSE
    Cel_Loc%lres   = FALSE
    Cel_Loc%lGd    = FALSE
    Cel_Loc%lMOX   = FALSE
    Cel_Loc%lAIC   = FALSE
    Cel_Loc%lrect  = FALSE
    Cel_Loc%nFSR   = nSub * nFSRinSub
    Cel_Loc%nFXR   = nSub
    Cel_Loc%icel0  = jcel
    
    Cel_Loc%nFSRinFXR = nFSRinSub ! ASSUME : # of FSR in Sub is same along FXRs
    
    Cel_Loc%Geom%lCircle = TRUE
    Cel_Loc%Geom%nCircle = nSub - 1 ! # of Sub-Rings
    
    CALL dmalloc(Cel_Loc%geom%Circle, 3, Cel_Loc%geom%nCircle)
    
    DO iSub = 1, Cel_Loc%Geom%nCircle
      Cel_Loc%geom%Circle(3, iSub) = cBs_Loc%sRad(iSub + 1)
    END DO
    
    jFSR = 0
    jSub = 0
    
    DO iFXR = 1, nFXR
      Cel_Loc%Vol(iFXR-1) = cBs_Loc%xVol(BndyTyp(iDir), iFXR) / rnFSRinSub(iDir)
      
      DO iSub = 1, cBs_Loc%xDiv(iFXR)
        jSub = jSub + 1
        
        DO iFSR = 1, nFSRinSub
          jFSR = jFSR + 1
          
          Cel_Loc%iReg(jFSR) = hCel(iCel)%xMix(iFXR)
          
          Cel_Loc%MapFxr2FsrIdx(iFSR, jSub) = jFSR
        END DO
      END DO
    END DO
  END DO
  ! --------------------------------------------------
END DO
! ----------------------------------------------------
!               02. Gap Cel
! ----------------------------------------------------
DO iCel = 1, nGapType
  IF (.NOT. gCel(iCel)%luse) CYCLE
  
  gBs_Loc => gCelBss(gCel(iCel)%igBss)
  
  nSub = gBs_Loc%nSub
  nFXR = gBs_Loc%nFXR
  ! --------------------------------------------------
  DO iDir = 1, pnGap
    jCel = nCellType * pnRod + (iCel - 1) * pnGap + iDir
    
    Cel_Loc => CellInfo(jCel)
    
    nFSRinSub = (2 - iDir) * gBs_Loc%nVtxHor & ! Bndy Typ 1
              + (iDir - 1) * gBs_Loc%nHor / 2  ! Bndy Typ 2
    
    CALL dmalloc(Cel_Loc%Vol,           nFSRinSub * nSub) ! FSR unit
    CALL dmalloc(Cel_Loc%iReg,          nFSRinSub * nSub)
    CALL dmalloc(Cel_Loc%nFSRinFXR,     nSub)
    CALL dmalloc(Cel_Loc%MapFxr2FsrIdx, nFSRinSub, nSub)
    
    Cel_Loc%luse   = TRUE
    Cel_Loc%lfuel  = FALSE
    Cel_Loc%lres   = FALSE
    Cel_Loc%lGd    = FALSE
    Cel_Loc%lMOX   = FALSE
    Cel_Loc%lAIC   = FALSE
    Cel_Loc%lrect  = FALSE
    Cel_Loc%nFXR   = nSub
    Cel_Loc%nFSR   = nFSRinSub * nSub
    
    Cel_Loc%nFSRinFXR    = nFSRinSub
    Cel_Loc%Geom%lCircle = FALSE
    
    jFSR = 0
    jSub = 0
    
    DO iFXR = 1, nFXR
      DO iSub = 1, gBs_Loc%xDiv(iFXR)
        jSub = jSub + 1
        
        DO iFSR = 1, nFSRinSub
          jFSR = jFSR + 1
          
          Cel_Loc%iReg(jFSR) = gCel(iCel)%xMix(iFXR)
          Cel_Loc%Vol(jFSR)  = (2 - iDir) * gBs_Loc%sVol(1, jFSR) & ! Bndy Typ 1
                             + (iDir - 1) * gBs_Loc%sVol(2, jFSR)   ! Bndy Typ 2
          
          Cel_Loc%MapFxr2FsrIdx(iFSR, jSub) = jFSR
        END DO
      END DO
    END DO
  END DO
  ! --------------------------------------------------
END DO
! ----------------------------------------------------
!               03. SET : LOGICAL
! ----------------------------------------------------
DO iCel = 1, ncTyp
  Cel_Loc => CellInfo(iCel)
  
  IF (.NOT. Cel_Loc%luse) CYCLE
  
  DO iFSR = 1, Cel_Loc%nFSR
    iReg = Cel_Loc%iReg(iFSR)
    
    IF(nTracerCntl%lXsLib) THEN
      Cel_Loc%lfuel = Cel_Loc%lfuel .OR. Mixture(iReg)%lfuel
      Cel_Loc%lres  = Cel_Loc%lres  .OR. Mixture(iReg)%lres
      Cel_Loc%lGd   = Cel_Loc%lGd   .OR. Mixture(iReg)%lGd
      Cel_Loc%lMOX  = Cel_Loc%lMOX  .OR. Mixture(iReg)%lMOX
      Cel_Loc%lAIC  = Cel_Loc%lAIC  .OR. Mixture(iReg)%lAIC
      Cel_Loc%lsSPH = Cel_Loc%lfuel .OR. Cel_Loc%lAIC
    ELSE
      Cel_Loc%lfuel = Cel_Loc%lfuel .OR. MacXsBen(iReg)%lfuel
      Cel_Loc%lCR   = Cel_Loc%lCR   .OR. MacXsBen(iReg)%lCR
    END IF
  END DO
END DO
! ----------------------------------------------------
!               04. SET : SSPH Data
! ----------------------------------------------------
DO iCel = 1, nCellType
  IF (nTracerCntl%lRestrmt) CALL HexCalcCellSSPH(iCel)
END DO
! ----------------------------------------------------
!               05. REMEDY : Vol
! ----------------------------------------------------
DO iCel = 1, nCellType
  IF (.NOT. hCel(iCel)%luse) CYCLE
  
  cBs_Loc => hCelBss(hCel(iCel)%icBss)
  
  nSct = cBs_Loc%nSct
  nSub = cBs_Loc%nSub
  nFXR = cBs_Loc%nFXR
  ! --------------------------------------------------
  DO iDir = 1, pnRod
    IF (lSkip(iDir)) CYCLE
    
    jCel = (iCel - 1) * pnRod + iDir
    
    Cel_Loc => CellInfo(jCel)
    
    nFSRinSub = nSct / rnFSRinSub(iDir)
    
    CALL dmalloc(Cel_Loc%Vol, Cel_Loc%nFSR)
    
    jFSR = 0
    
    DO iSub = 1, nSub
      DO iFSR = 1, nFSRinSub
        jFSR = jFSR + 1
        
        Cel_Loc%Vol(jFSR) = cBs_Loc%sVol(BndyTyp(iDir), iFSR, iSub)
      END DO
    END DO
  END DO
  ! --------------------------------------------------
END DO

NULLIFY (Cel_Loc)
NULLIFY (cBs_Loc)
NULLIFY (gBs_Loc)
! ----------------------------------------------------

END SUBROUTINE HexCnPCelinfo
! ------------------------------------------------------------------------------------------------------------
!                                     03. HEX CALC : Cell SSPH
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexCalcCellSSPH(icTyp)

USE allocs
USE TYPEDEF,      ONLY : cell_type
USE param,        ONLY : ZERO, PI, INVPI, FALSE, TRUE
USE GEOM,         ONLY : CellInfo
USE cntl,         ONLY : nTracerCntl
USE Material_Mod, ONLY : Mixture
USE BenchXs,      ONLY : MacXsBen

USE HexType,      ONLY : Type_HexRodCelBss, Type_HexRodCel
USE HexData,      ONLY : hCel, hCelBss, hLgc
USE SPH_mod,      ONLY : calcCellSSPH,calcAICCellSSPH

IMPLICIT NONE

INTEGER :: icTyp, jCel, iFXR, jFXR, iDir, iReg, i, j, k
INTEGER :: nFXR, nmat, nreg, nFuelDiv, nfDiv, ibFuel, ieFuel
INTEGER :: ndiv(0:300), iRg(0:300)
INTEGER :: imat(100), idiv(100)

REAL :: vol, vol_
REAL :: rr(300), irad(0:100)

TYPE(cell_type),         POINTER :: Cel_Loc, Cel_CP
TYPE(Type_HexRodCel),    POINTER :: hCl_Loc
TYPE(Type_HexRodCelBss), POINTER :: cBs_Loc
! ----------------------------------------------------

IF (.NOT. hCel(icTyp)%luse) RETURN

hCl_Loc => hCel(icTyp)
cBs_Loc => hCelBss(hCel(icTyp)%icBss)

jCel = pnRod * (icTyp - 1) + 1
nFXR = cBs_Loc%nFXR

IF (hLgc%lSngCel) jCel = pnRod * (icTyp - 1) + 3
! ----------------------------------------------------
!               01. SET : INP Data
! ----------------------------------------------------
rr   = ZERO
nDiv = 0
iRg  = 0

DO iFXR = 1, nFXR - 1
  rr  (iFXR) = hCl_Loc%xRad(iFXR + 1)
  nDiv(iFXR) = hCl_Loc%xDiv(iFXR + 1)
END DO

DO iFXR = 0, nFXR
  iRg (iFXR) = hCl_Loc%xMix(iFXR + 1)
END DO

imat = 0
irad = 0._8
idiv = 0

imat(1) = iRg (nFXR - 1)
irad(1) = rr  (nFXR - 1)
idiv(1) = ndiv(nFXR - 1)
nmat    = 1

DO iFXR = nFXR - 2, 0, -1
  IF (iRg(iFXR) .eq. iRg(iFXR + 1)) THEN
    IF (iFXR .eq. 0) THEN
      idiv(nmat) = idiv(nmat) + 1
      irad(nmat) = dsqrt(cBs_Loc%cVol / PI)
    ELSE
      idiv(nmat) = idiv(nmat) + ndiv(iFXR)
    END IF
    
    CYCLE
  END IF
  
  nmat = nmat + 1
  imat(nmat) = iRg(iFXR)
  
  IF (iFXR .gt. 0) THEN
    irad(nmat) = rr(iFXR)
    idiv(nmat) = ndiv(iFXR)
  ELSE
    irad(nmat) = dsqrt(cBs_Loc%cVol / PI)
    idiv(nmat) = 1
  END IF
END DO
! ----------------------------------------------------
!               02. CnP : INP Data
! ----------------------------------------------------
Cel_Loc => CellInfo(jCel)

Cel_Loc%nmat = nmat

CALL dmalloc(Cel_Loc%matidx, nmat)
CALL dmalloc(Cel_Loc%matrad, nmat)

Cel_Loc%matidx(1:nmat) = imat(1:nmat)
Cel_Loc%matrad(1:nmat) = irad(1:nmat)

nreg = sum(idiv(1:nmat))

! Sub-ring Data
Cel_Loc%nreg_cp = nreg

CALL dmalloc(Cel_Loc%rad_cp, nreg)
CALL dmalloc(Cel_Loc%fxrvol, nreg)

vol = 0._8
k   = 1

DO i = 1, nmat
  vol_ = PI * (irad(i) * irad(i) - irad(i-1) * irad(i-1)) / idiv(i)
  
  DO j = 1, idiv(i)
    vol = vol + vol_
    
    Cel_Loc%rad_cp(k) = dsqrt(vol * INVPI)
    Cel_Loc%fxrvol(k) = vol_
    
    k = k + 1
  END DO
END DO

! Others
IF (Cel_Loc%lfuel) THEN
  DO i = nmat, 1, -1
    IF (nTracerCntl%lXsLib) THEN
      IF (Mixture(imat(i))%lfuel) EXIT
    END IF
  END DO
  
  DO k = i + 1, nmat
    IF (nTracerCntl%lXsLib) THEN
      IF (Mixture(imat(k))%lh2o) EXIT
    END IF
  END DO
  
  k = sum(idiv(1:k-1));
  
  Cel_Loc%fuelgapcldvol = sum(Cel_Loc%fxrvol(1:k))
  Cel_Loc%cldfxridx     = k
  Cel_Loc%nmodfxr       = nreg-k
  Cel_Loc%invnmodfxr    = 1._8 / real(Cel_Loc%nmodfxr, 8)
ELSE IF (Cel_Loc%lAIC) THEN
  DO i = nmat, 1, -1
    IF (nTracerCntl%lXsLib) THEN
      IF (Mixture(imat(i))%lAIC) EXIT 
    END IF
  END DO
  
  DO k = i + 1, nmat
    IF (nTracerCntl%lXsLib) THEN
      IF (Mixture(imat(k))%lh2o) EXIT
    END IF
  END DO
  
  k = sum(idiv(1:k-1));
  
  Cel_Loc%fuelgapcldvol = sum(Cel_Loc%fxrvol(1:k))
  Cel_Loc%cldfxridx     = k
  Cel_Loc%nmodfxr       = nreg - k
  Cel_Loc%invnmodfxr    = 1._8 / real(Cel_Loc%nmodfxr, 8)
END IF
! ----------------------------------------------------
!               03. FIND : ibFuel / ieFuel
! ----------------------------------------------------
IF (Cel_Loc%lfuel) THEN
  Cel_Loc%lhole = .NOT. Mixture(iRg(nFXR-1))%lfuel
  
  DO j = nFXR-1, 0, -1
    IF (Mixture(iRg(j))%lfuel) EXIT
  END DO
  
  Cel_Loc%ibFuel = j
  
  DO j = nFXR-2, 0, -1
    IF(nTracerCntl%lXsLib) THEN
      IF (Mixture (iRg(j+1))%lfuel .AND. .NOT.Mixture (iRg(j))%lfuel) EXIT
    ELSE
      IF (MacXsBen(iRg(j+1))%lfuel .AND. .NOT.MacXsBen(iRg(j))%lfuel) EXIT
    ENDIF
  END DO
  
  j = j + 1
  
  Cel_Loc%ieFuel   = j
  Cel_Loc%FuelRad0 = rr(Cel_Loc%ieFuel)
ELSE IF (Cel_Loc%lAIC) THEN
  DO j = nFXR-1, 0, -1
    IF (.not.Mixture(iRg(j))%lAIC) EXIT ! Based on the assumption that there's no hole in an AIC pin.
  ENDDO
  
  j = j + 1
  
  Cel_Loc%ieFuel   = j
  Cel_Loc%ibFuel   = nFXR-1
  Cel_Loc%FuelRad0 = rr(Cel_Loc%ieFuel)
END IF

IF (.NOT. nTracerCntl%lXsLib) THEN
  Cel_Loc%ibFuel = nFXR-1
  Cel_Loc%ieFuel = 1
END IF
! ----------------------------------------------------
!               04. CP : Cell Info
! ----------------------------------------------------
DO iDir = 2, pnRod
  Cel_CP => CellInfo(pnRod * (icTyp - 1) + iDir)
  
  IF (.NOT. CeL_CP%luse) CYCLE
  
  CALL dmalloc(Cel_CP%matidx, nmat)
  CALL dmalloc(Cel_CP%matrad, nmat)
  CALL dmalloc(Cel_CP%rad_cp, nreg)
  CALL dmalloc(Cel_CP%fxrvol, nreg)
  
  Cel_CP%nmat     = Cel_Loc%nmat
  Cel_CP%matidx   = Cel_Loc%matidx
  Cel_CP%matrad   = Cel_Loc%matrad
  Cel_CP%nreg_cp  = Cel_Loc%nreg_cp
  Cel_CP%rad_cp   = Cel_Loc%rad_cp
  Cel_CP%fxrvol   = Cel_Loc%fxrvol
  
  Cel_CP%fuelgapcldvol = Cel_Loc%fuelgapcldvol
  Cel_CP%cldfxridx     = Cel_Loc%cldfxridx
  Cel_CP%nmodfxr       = Cel_Loc%nmodfxr
  Cel_CP%invnmodfxr    = Cel_Loc%invnmodfxr
  
  Cel_CP%FuelRad0 = Cel_Loc%FuelRad0
  
  Cel_CP%lHole  = Cel_Loc%lHole
  Cel_CP%ibFuel = Cel_Loc%ibFuel
  Cel_CP%ieFuel = Cel_Loc%ieFuel
END DO
! ----------------------------------------------------
!               05. CAL : SSPH
! ----------------------------------------------------
IF (.NOT.(nTracerCntl%lXsLib .AND. nTracerCntl%lsSPH)) RETURN
IF (.NOT. Cel_Loc%lFuel) RETURN

ieFuel = Cel_Loc%ieFuel
ibFuel = Cel_Loc%ibFuel

nfueldiv = 0

DO k = ibFuel, ieFuel, -1
    nfueldiv = nfueldiv + ndiv(k)
END DO

DO iDir = 1, pnRod
  jCel = pnRod * (icTyp - 1) + iDir
  
  IF (.NOT. CellInfo(jCel)%luse) CYCLE
  
  nfDiv = nfueldiv
  
  IF (CellInfo(jCel)%lAIC) THEN
    CALL calcAICCellSSPH(CellInfo, jcel, nFXR - 1, iRg, rr, ndiv, iFXR - 1, nfDiv)
  ELSE
    CALL calcCellSSPH   (CellInfo, jcel, nFXR - 1, iRg, rr, ndiv, ibFuel, ieFuel, nfDiv, Cel_Loc%lhole)
  END IF
END DO

NULLIFY (Cel_Loc)
NULLIFY (hCl_Loc)
NULLIFY (cBs_Loc)
! ----------------------------------------------------

END SUBROUTINE HexCalcCellSSPH
! ------------------------------------------------------------------------------------------------------------
!                                     04. HEX CnP : Pin
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexCnPPin(iPin)

USE allocs
USE TYPEDEF, ONLY : pin_type
USE GEOM,    ONLY : Core, Pin, nZ, nCellType, nGapType, CellInfo, AsyInfo
USE HexType, ONLY : Type_HexPinInfo
USE HexData, ONLY : hPinInfo, RodPin, GapPin, hCel, hAsyTypInfo, hLgc

IMPLICIT NONE

INTEGER :: iz, iPin, ipTyp, ivTyp, iCel, iaTyp
INTEGER :: ncRod

TYPE(pin_type),        POINTER :: Pin_Loc
TYPE(Type_HexPinInfo), POINTER :: Inf_Loc
! ----------------------------------------------------

Pin_Loc =>      Pin(iPin)
Inf_Loc => hPinInfo(iPin)

ncRod = pnRod * nCellType
ipTyp = Inf_Loc%PinTyp
ivTyp = Inf_Loc%VtxTyp

! Basic
Pin_Loc%iAsy     = Inf_Loc%AsyIdx
Pin_Loc%PinType  = Inf_Loc%PinTyp ! for DEBUG
Pin_Loc%AsyType  = Inf_Loc%AsyTyp
Pin_Loc%FsrIdxSt = Inf_Loc%FsrIdxSt
Pin_Loc%FxrIdxSt = Inf_Loc%FxrIdxSt

! Cel
CALL dmalloc(Pin_Loc%Cell,    nZ)
CALL dmalloc(Pin_Loc%hCelGeo, nZ)

DO iz = 1, nZ
  ! Case : Rod
  IF (Inf_Loc%lRod) THEN
    iCel = RodPin(ipTyp)%iCel(iz)
    
    Pin_Loc%hCelGeo(iz) = hCel(iCel)%icBss
    Pin_Loc%Cell   (iz) = pnRod * (iCel - 1) + ivTyp - ivTyp / 7
  ! Case : Gap
  ELSE
    iCel = GapPin(ipTyp)%iCel(iz)
    
    ! Bndy Typ 1 = ivTyp 8, 9
    ! Bndy Typ 2 = ivTyp 10
    Pin_Loc%hCelGeo(iz) = hAsyTypInfo(Inf_Loc%AsyTyp)%iBss
    Pin_Loc%Cell   (iz) = ncRod + pnGap * (iCel - 1) + 1 + ivTyp / 10
  END IF
END DO

IF (hLgc%lSngCel) THEN
  DO iCel = 1, nCellType
    IF (hCel(iCel)%luse) EXIT
  END DO
  
  Pin_Loc%Cell = (iCel - 1) * pnRod + 3 ! # of used hCel must be 1 for Sng Cel
END IF

! LOGICAL
CALL dmalloc(Pin(iPin)%lMOX, Core%nz)

DO iz = 1, Core%nz
  iCel = Pin(iPin)%Cell(iz)
  
  Core%lFuelPlane(iz) =  Core%lFuelPlane(iz) .OR. CellInfo(icel)%lfuel
  Core% lAICPlane(iz) =  Core% lAICPlane(iz) .OR. CellInfo(icel)%lAIC
  
  Pin(iPin)%lfuel = Pin(iPin)%lfuel .OR. CellInfo(icel)%lfuel
  Pin(iPin)%lGd   = Pin(iPin)%lGd   .OR. CellInfo(icel)%lGd
  
  Pin(iPin)%lMox(iz) = CellInfo(icel)%lMox
END DO

iaTyp = hPinInfo(iPin)%AsyTyp

AsyInfo(iaTyp)%lFuel = AsyInfo(iaTyp)%lFuel .OR. Pin(iPin)%lfuel

! FIN
NULLIFY (Pin_Loc)
NULLIFY (Inf_Loc)
! ----------------------------------------------------

END SUBROUTINE HexCnPPin
! ------------------------------------------------------------------------------------------------------------
!                                     05. CONVERT : Ray
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE ConvertRay()

USE allocs
USE PARAM,    ONLY : PI, VoidCell
USE Geom,     ONLY : Core, Albedo
USE Rays,     ONLY : RayInfo, PolarAngle, AziAngle, RotRay, CoreRay, AsyRay
USE Moc_Mod,  ONLY : nMaxCoreRay, nMaxCellRay, nMaxRaySeg
USE SetRay,   ONLY : SetPolarRayAngle
USE HexData,  only : ncRay, nRotRay, hcRay, hRotRay, NumMray, nAzmAng, AzmAng, AzmDel, AzmSin, AzmCos, AzmTan, AzmWgt, &
                     nPolAng, nRotRay

IMPLICIT NONE

INTEGER :: iAzm, jAzm, iRotRay, iCoreRay, imRay, iPin
INTEGER :: iAssy, iAssyType, iAssyDatID
INTEGER :: nCoreRay, nMRay, nSeg, iSegSt, nPin, nhcPin, nhcSeg, nBndy
! ----------------------------------------------------

! Ray Info
RayInfo%nModRay  = NumMRay(0)
RayInfo%nCoreRay = ncRay
RayInfo%nRotRay  = nRotRay

RayInfo%nPhiAngSv = 2*RayInfo%nRotRay + 1

CALL dmalloc(RayInfo%PhiAngInSvIdx,  RayInfo%nRotRay, 2)
CALL dmalloc(RayInfo%PhiAngOutSvIdx, RayInfo%nRotRay, 2)

IF (Albedo(2) .EQ. VoidCell) THEN
  RayInfo%PhiAngInSvIdx = 1
  
  DO iRotRay = 1, RayInfo%nRotRay
    RayInfo%PhiAngOutSvIdx(iRotRay, 1) = 2*iRotRay + 1
    RayInfo%PhiAngOutSvIdx(iRotRay, 2) = 2*iRotRay
  END DO
ELSE
  DO iRotRay = 1, RayInfo%nRotRay
    RayInfo%PhiAngInSvIdx (iRotRay, 1) = 2*iRotRay
    RayInfo%PhiAngOutSvIdx(iRotRay, 1) = 2*iRotRay
    RayInfo%PhiAngInSvIdx (iRotRay, 2) = 2*iRotRay + 1
    RayInfo%PhiAngOutSvIdx(iRotRay, 2) = 2*iRotRay + 1
  END DO
END IF

! Azm Angle
RayInfo%nAziAngle = nAzmAng

ALLOCATE (AziANgle (RayInfo%nAziAngle * 2))

DO iAzm = 1, RayInfo%nAziAngle
  AziAngle(iAzm)%Ang    = AzmAng(iAzm)
  AziAngle(iAzm)%Del    = AzmDel(iAzm)
  AziAngle(iAzm)%sinv   = AzmSin(iAzm)
  AziAngle(iAzm)%cosv   = AzmCos(iAzm)
  AziAngle(iAzm)%tanv   = AzmTan(iAzm)
  AziAngle(iAzm)%Weight = AzmWgt(iAzm)
END DO

DO iAzm = 1, RayInfo%nAziAngle
  jAzm = 2*RayInfo%nAziAngle - iAzm + 1

  AziAngle(jAzm)%Ang    = Pi - AziAngle(iAzm)%Ang
  AziAngle(jAzm)%Del    = AziAngle(iAzm)%Del
  AziAngle(jAzm)%sinv   = SIN(AziAngle(jAzm)%ang)
  AziAngle(jAzm)%cosv   = COS(AziAngle(jAzm)%ang)
  AziAngle(jAzm)%tanv   = TAN(AziAngle(jAzm)%ang)
  AziAngle(jAzm)%Weight = AziAngle(iAzm)%Weight
END DO

! Polar Angle
RayInfo%nPolarAngle     = nPolAng
RayInfo%nPolarAngleHemi = nPolAng / 2

CALL SetPolarRayAngle(PolarAngle, RayInfo)

! RayInfo 4 CMFD - Rough
ALLOCATE (RayInfo%RayInfo4CMFD)

RayInfo%RayInfo4CMFD%nRotRay   = nRotRay
RayInfo%RayInfo4CMFD%nPolAngle = RayInfo%nPolarAngle

! Pointing
RayInfo%AziAngle => AziAngle
RayInfo%RotRay   => RotRay
RayInfo%CoreRay  => CoreRay
RayInfo%AsyRay   => AsyRay
! ----------------------------------------------------

END SUBROUTINE ConvertRay
! ------------------------------------------------------------------------------------------------------------
!                                     06. CONVERT : Xs
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE ConvertXs

USE Param,        ONLY : TRUE, FALSE
USE Material_Mod, ONLY : Mixture, nMixType
USE BenchXs,      ONLY : nxsltype
USE CNTL,         ONLY : nTracerCntl

IMPLICIT NONE
! ----------------------------------------------------

IF (nTracerCntl%libtyp .EQ. 1) THEN
  ALLOCATE (Mixture (nxsltype))

  Mixture(1:nxsltype)%lfuel = TRUE
  Mixture(1:nxsltype)%lGd   = FALSE
  !Mixture(2:7)%lMOX  = TRUE
  !Mixture(6:8)%lfuel = FALSE
END IF
! ----------------------------------------------------

END SUBROUTINE ConvertXs
! ------------------------------------------------------------------------------------------------------------

END MODULE HexCnP