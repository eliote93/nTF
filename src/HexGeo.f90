SUBROUTINE HexSetGeo()

USE PARAM,   ONLY : ZERO, FALSE
USE ioutil,  ONLY : terminate
USE geom,    ONLY : nAsyType0, Core, Asy, AsyInfo, AsyVol, PinInfo, CellInfo, Pin, nCellType, nZ, nPinType
USE HexData, ONLY : nHexPin, nhAsy, nHexPin, RodPin, hAsyTypInfo
USE HexTst,  ONLY : HexPrintPinTyp, HexTsthPinInfo

USE HexCP
USE HexGeoConst
USE HexAsyConst
USE HexAsyTypConst
USE HexPinConst
USE HexVtx
USE HexNgh

IMPLICIT NONE

INTEGER :: iAsy, iaTyp, iPin, ipTyp
INTEGER :: FsrIdxSt, FxrIdxSt
! ----------------------------------------------------
!                1. SET : Asy
! ----------------------------------------------------
CALL HexSetAsyBndy()
CALL HexSetGeoTyp()
CALL HexSetAsyLoc()

DO iaTyp = 1, nAsyType0
  IF (.NOT. hAsyTypInfo(iaTyp)%luse) CYCLE
  
  CALL HexSetAsyTypPinMap(iaTyp)
  CALL HexSetAsyTypPinLocIdx(iaTyp)
  CALL HexSetAsyTypVtxTyp(iaTyp)
  CALL HexSetAsyTypPinVtx(iaTyp)
  CALL HexSetAsyTypPinNgh(iaTyp)
END DO

DO iAsy = 1, nhAsy
  CALL HexSetAsyPinNum(iAsy)
  CALL HexSetAsyNgh(iAsy)
  CALL HexSetAsyRotNgh(iAsy)
END DO
! ----------------------------------------------------
!                2. SET : Core
! ----------------------------------------------------
CALL HexSetAlbedo()
CALL HexSetCore()
! ----------------------------------------------------
!                3. SET : Pin & Msh
! ----------------------------------------------------
ALLOCATE (Core%lFuelPlane        (nz)); Core%lFuelPlane = FALSE
ALLOCATE (Core%lCladPlane        (nz)); Core%lCladPlane = FALSE
ALLOCATE (Core% lAICPlane        (nz)); Core% lAICPlane = FALSE
ALLOCATE (Core%PinVol   (nHexPin, nZ)); Core%PinVol     = ZERO
ALLOCATE (Core%PinVolFm (nHexPin, nZ)); Core%PinVolFm   = ZERO

ALLOCATE (Pin (nHexPin))

DO iAsy = 1, nhAsy
  CALL HexSetPinTyp(iAsy)
  CALL HexSetPinVol(iAsy)
END DO

FsrIdxSt = 1
FxrIdxSt = 1

DO iPin = 1, nHexPin
  CALL HexSetPinFSR(iPin, FsrIdxSt, FxrIdxSt)
END DO

CALL HexSetVyg()
CALL HexSetVss() ! Vss must follow Vyg

!CALL HexTsthPinInfo()
!CALL HexPrintPinTyp()
! ----------------------------------------------------
!                4. CP
! ----------------------------------------------------
CALL HexCPCelInfo()

DO iPin = 1, nHexPin
  CALL HexCPPin(iPin)
END DO
! ----------------------------------------------------
!                5. FIN : Core
! ----------------------------------------------------
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

END SUBROUTINE HexSetGeo