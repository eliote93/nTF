SUBROUTINE HexSetGeo()

USE allocs
USE PARAM,   ONLY : ZERO, FALSE
USE geom,    ONLY : Core, Pin, nZ, nAsyType0
USE HexData, ONLY : nHexPin, nhAsy, nHexPin, hAsyTypInfo

IMPLICIT NONE

INTEGER :: iAsy, iaTyp, iPin, ipTyp
INTEGER :: FsrIdxSt, FxrIdxSt
! ----------------------------------------------------

! Asy
CALL HexSetAsyBndy
CALL HexSetGeoTyp
CALL HexSetAsyLoc

DO iaTyp = 1, nAsyType0
  IF (.NOT. hAsyTypInfo(iaTyp)%luse) CYCLE
  
  CALL HexSetAsyTypPinMap   (iaTyp)
  CALL HexSetAsyTypPinLocIdx(iaTyp)
  CALL HexSetAsyTypVtxTyp   (iaTyp)
  CALL HexSetAsyTypPinVtx   (iaTyp)
  CALL HexSetAsyTypPinNgh   (iaTyp)
  CALL HexSetAsyTypCstMap   (iaTyp)
END DO

DO iAsy = 1, nhAsy
  CALL HexSetAsyPinNum(iAsy)
  CALL HexSetAsyNgh   (iAsy)
  CALL HexSetAsyRotNgh(iAsy)
END DO

! Core
CALL HexSetAlbedo
CALL HexSetCore

! Pin & Mesh
CALL dmalloc(Core%lFuelPlane, nz)
CALL dmalloc(Core%lCladPlane, nz)
CALL dmalloc(Core% lAICPlane, nz)

CALL dmalloc(Core%PinVol,   nHexPin, nZ)
CALL dmalloc(Core%PinVolFm, nHexPin, nZ)

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
! ----------------------------------------------------

END SUBROUTINE HexSetGeo