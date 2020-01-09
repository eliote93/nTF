MODULE HexGeoConst

USE ioutil, ONLY : terminate

IMPLICIT NONE

CONTAINS

! ------------------------------------------------------------------------------------------------------------
!                                     01. HEX SET : Albedo
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetAlbedo()

USE PARAM,   ONLY : FALSE, TRUE, VOIDCELL, REFCELL, ROTCELL, CBDCELL, ZERO
USE geom,    ONLY : Core, lcbd, Albedo
USE HexData, ONLY : hLgc

IMPLICIT NONE

! ----------------------------------------------------
Core%RadSym = FALSE
Core%RadBC  = VoidCell
Core%AxBC   = RefCell

IF (hLgc%lAxRef(1)) THEN
  Albedo(3)    = RefCell
  Core%AxBC(1) = RefCell
ELSE
  Albedo(3)    = VoidCell
  Core%AxBC(1) = VoidCell
END IF

IF (hLgc%lAxRef(2)) THEN
  Albedo(4)    = RefCell
  Core%AxBC(2) = RefCell
ELSE
  Albedo(4)    = VoidCell
  Core%AxBC(2) = VoidCell
END IF

IF (hLgc%lAzmRot) THEN
  Albedo(1) = RotCell
ELSE
  Albedo(1) = RefCell
END IF

IF (hLgc%lRadRef) THEN
  Albedo(2)       = RefCell
  Core%RadSym     = TRUE
  Core%RadBC(1:6) = TRUE
ELSE
  Albedo(2) = VoidCell
END IF

IF(lCBD) THEN
  Core%RadBC(1:6) = CbdCell

  IF(Albedo(1) .NE. ZERO) THEN
    CALL TERMINATE('All Radial B.C should be reflective in the Checkerboard Configuration')
  ENDIF
ENDIF

IF(Albedo(7) .NE. ZERO) Core%AxBc(1) = VoidCell
IF(Albedo(8) .NE. ZERO) Core%AxBc(2) = VoidCell

END SUBROUTINE HexSetAlbedo
! ------------------------------------------------------------------------------------------------------------
!                                     02. HEX SET : Core
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetCore()

USE geom,    ONLY : nxy, nZ, nAsyType0, nPinType, Core, Asy, AsyInfo, AsyVol, PinInfo, &
                    hz, hzInv, nzFm, hzFm, hzFmInv, nSubPlane, SubPlaneMap, SubPlaneRange
USE PE_MOD,  ONLY : PE
USE HexData, ONLY : nhAsy, nHexPin, hAsy, hAsyTypInfo, Asy1Dto2DMap, hPinInfo, hGeoTypInfo

IMPLICIT NONE

INTEGER :: iAsy, iPin, jPin, iAsyTyp, iGeoTyp
INTEGER :: nPin
! ----------------------------------------------------

nxy   = nHexPin
PE%nz = nZ

Core%nxy      = nxy
Core%nz       = nZ
Core%nxya     = nhAsy
Core%nAsyType = nAsyType0
Core%nCoreFSR = 0
Core%nCoreFXR = 0
! ----------------------------------------------------
!                1. ALLOC
! ----------------------------------------------------
ALLOCATE (Asy      (nhAsy))
ALLOCATE (AsyInfo  (nAsyType0))
ALLOCATE (AsyVol   (nhAsy, nZ))
ALLOCATE (PinInfo  (nPinType))
ALLOCATE (hPinInfo (nHexPin))

DO iPin = 1, nHexPin
  ALLOCATE (hPinInfo(iPin)%Vol   (nZ))
  ALLOCATE (hPinInfo(iPin)%VolFm (nZ))
END DO
! ----------------------------------------------------
!                2. SET : Ax
! ----------------------------------------------------
CALL InitSubPlane()

Core%hz      => hz
Core%hzInv   => hzInv
Core%nzfm    =  nzfm
Core%hzfm    => hzfm
Core%hzfmInv => hzfmInv

Core%nSubPlane      = nSubPlane
Core%SubPlaneMap   => SubPlaneMap
Core%SubPlaneRange => SubPlaneRange
! ----------------------------------------------------
!                3. SET : Global Numeric #
! ----------------------------------------------------
nPin = 0

DO iAsy = 1, nhAsy
  iAsyTyp = hAsy(iAsy)%AsyTyp
  iGeoTyp = hAsy(iAsy)%GeoTyp
  ! ------------------------
  !      Asy
  ! ------------------------
  Asy(iAsy)%AsyType = iAsyTyp
  Asy(iAsy)%ixa     = Asy1Dto2DMap(1, iAsy)
  Asy(iAsy)%iya     = Asy1Dto2DMap(2, iAsy)
  Asy(iAsy)%wt      = hAsy(iAsy)%wt
  
  AsyVol(iAsy, :) = hGeoTypInfo(iGeoTyp)%Area
  
  AsyInfo(iAsyTyp)%nxy = hAsyTypInfo(iAsyTyp)%nTotPin(1)
  
  ALLOCATE (Asy(iAsy)%GlobalPinIdx (AsyInfo(iAsyTyp)%nxy))
  
  Asy(iAsy)%GlobalPinIdx(:) = -1
  ! ------------------------
  !      Pin
  ! ------------------------
  jPin = 0
  
  DO iPin = 1, hAsyTypInfo(iAsyTyp)%nTotPin(1)
    IF (hAsyTypInfo(iAsyTyp)%lGeoPin(iGeoTyp, iPin)) THEN
      nPin = nPin + 1
      jPin = jPin + 1
      
      hPinInfo(nPin)%AsyIdx     = iAsy
      hPinInfo(nPin)%OrdInAsy01 = iPin
      
      Asy(iAsy)%GlobalPinIdx(jPin) = nPin
    END IF
  END DO
END DO

END SUBROUTINE HexSetCore
! ------------------------------------------------------------------------------------------------------------
!                                     03. SET : Vyg
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetVyg(iAsy)

USE GEOM,    ONLY : nz
USE HexData, ONLY : hLgc, VygAsyTyp, VygMat, VygFxrSt, hAsy, hPinInfo, hAsyTypInfo, GapPin

IMPLICIT NONE

INTEGER :: iAsy

INTEGER :: iAsyTyp, jAsyTyp, iDir, iPin, iDirTyp, iTmp, iz
INTEGER :: PinIdxSt, nRodGlb, nTotGlb
! ----------------------------------------------------

!IF (.NOT. hLgc%lVyg) RETURN
!
!iAsyTyp = hAsy(iAsy)%AsyTyp
!
!IF (iAsyTyp .NE. VygAsyTyp) RETURN
!
!PinIdxSt = hAsy(iAsy)%PinIdxSt - 1
!nRodGlb  = hAsy(iAsy)%nRodPin + PinIdxSt
!nTotGlb  = hAsy(iAsy)%nTotPin + PinIdxSt
!
!DO iDir = 1, hAsy(iAsy)%nNgh
!  jAsyTyp = hAsy(hAsy(iAsy)%NghIdx(iDir))%AsyTyp
!  
!  IF (jAsyTyp < VygAsyTyp) THEN
!    DO iPin = nRodGlb + 1, nTotGlb
!      iDirTyp = hPinInfo(iPin)%DirTyp
!      
!      IF (iDirTyp < 19) THEN
!        iTmp = iDirTyp - 12
!      ELSE
!        iTmp = iDirTyp - 18
!      END IF
!      
!      IF (iTmp .EQ. hAsy(iAsy)%NghDir(iDir) - 6) THEN
!        DO iz = 1, nZ
!          !hPinInfo(iPin)%FxrRegIdx(VygFxrSt:nFXR, iz) = VygMat
!          hPinInfo(iPin)%lVyg = .TRUE.
!        END DO
!      END IF
!    END DO
!  END IF
!END DO

END SUBROUTINE HexSetVyg
! ------------------------------------------------------------------------------------------------------------
!                                     04. SET : TH geo
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetThGeo(nach, rw, rgt, acf, xi)

USE PARAM,   ONLY : PI, epsm4
USE HexData, ONLY : hLgc, hGeoTypInfo, hAsyTypInfo

IMPLICIT NONE

INTEGER :: nach

REAL :: rw, rgt, acf, xi
REAL :: Area
! ----------------------------------------------------

IF (.NOT. hLgc%lspCMFD) CALL terminate("SET TH GEO")

! ASSUME : Area of Fuel Pin and GT Pin are same

Area = hGeoTypInfo(1)%Area * epsm4
nach = hAsyTypInfo(1)%nRodPin(1)

! IGNORE : rgt
acf = (Area - PI * rw * rw * nach) / nach
xi  = 2._8 * PI * rw

END SUBROUTINE HexSetThGeo
! ------------------------------------------------------------------------------------------------------------

END MODULE HexGeoConst