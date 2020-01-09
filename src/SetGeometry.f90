SUBROUTINE SetGeometries()
!Set up Problem Geometries - Map(Cell, Assemblies)
USE PARAM
USE CNTL,           ONLY : nTracerCntl
USE GEOM,           ONLY : lEdge, lhgap
USE Files,          ONLY : io8
USE IOUTIL,         ONLY : message
USE PE_MOD,         ONLY : PE
USE MPIConfig_Mod,  ONLY : SetGeomPEVariables,  SetDcmpPEVariables
IMPLICIT NONE
IF(PE%Master) THEN
  mesg = 'Initialize the Problem Geometries'
  CALL message(io8,TRUE, TRUE, mesg)
ENDIF

!Initialize User-Defined Assembly Gap Structure
IF(nTracerCntl%lUsrGap) CALL SetUsrAsyGap()
!Set homogenized Gap geometries
IF(lhgap) CALL SethGapGeometries()

!Set Partial Cell & Assembly Geometry from 'CENT' OPTION
IF(.NOT. lEdge) CALL SetCentGeometries()
!Set Geometry - Intra-Cell Structure & Pin
!--- JSR Edit : nTIG Restart
IF (nTracerCntl%lnTIGRst) THEN
  CALL SetBasePinGeom()
ELSE
  CALL SetPinGeom()
ENDIF
!Set Assembly Geometry
CALL SetAsyGeom()

CALL SetGeomPEVariables(PE)
IF(nTracerCntl%lXsLib) THEN
  IF(nTracerCntl%lrestrmt) THEN
    IF(PE%RTmaster) CALL AllocResIsoInfo()
  ENDIF
ENDIF
!SubPlane
CALL InitSubPlane()
!Global Pin Cell Map
CALL GlobalCoreMap()
!Pin TH Type

!Core Information
CALL SetCoreInfo()

CALL SetThPinType()
!Set Up Baffle Information
CALL SetGlobalFsr_Fxr()

CALL SetCoreMiscStruct()

!--- CNJ Edit : Domain Decomposition + MPI
CALL SetDcmpPEVariables(PE)

END SUBROUTINE

SUBROUTINE CoreVolGen()
!Generates(or Calculate) the Pin and Assembly Volumes
USE PARAM
USE GEOM,  ONLY : Core,       CellInfo,                                          &
                  Asy,        AsyGeom,       PinInfo,          Pin,              &
                  nCellType,  nPinType,      nAsyType,         PinVol,           &
                  AsyVol,     PinVolFM
USE PE_MOD, ONLY :  PE
USE ALLOCS
IMPLICIT NONE
INTEGER :: nxy, nxya, nz, nzfm,nSUbPlane
INTEGER :: ixy, icel, ipin, ixya, iasy,iz

REAL :: Vol, H
REAL :: CelVol(500)
REAL, POINTER :: hz(:)

nz = Core%nz; nzfm = Core%nzfm
nxy = Core%nxy; nxya = Core%nxya
hz => Core%hz
nSUbPlane = Core%nSubPlane

CALL DMALLOC(PinVol, nxy, nz); CALL DMALLOC(AsyVol, nxya, nz)
IF(nz .ne. nzfm) THEN
  CALL DMALLOC(PinVolFm, nxy, nz)
ELSE
  PinVolFm => PinVOl
  !CALL DMALLOC(PinVolFm, nxy, nz)
ENDIF
!Cell Volume(area)
DO iCel = 1, nCellType
  !CelVol(icel) = CellInfo(icel)%Geom%lx*CellInfo(icel)%Geom%ly
ENDDO

  !H =Core%hz(iz)
DO ixy = 1, nxy
  DO iz = 1, nz
    icel = Pin(ixy)%Cell(iz)
    ixya = Pin(ixy)%iasy
    Vol = CellInfo(icel)%Geom%lx*CellInfo(icel)%Geom%ly * Core%hz(iz)
    PinVol(ixy, iz) = Vol
    PinVolFM(ixy, iz) = Vol/nSubPlane
    AsyVol(ixya, iz) = AsyVol(ixya, iz) + Vol
  ENDDO
ENDDO


CONTINUE
Core%PinVol => PinVol; Core%PinVolFM => PinVolFm
Core%AsyVol => AsyVol
END SUBROUTINE

SUBROUTINE SetCentGeometries()
!Set up Half pincell and half assemblies
USE PARAM
USE CNTL,  only : nTracerCntl
USE GEOM,  only : lEdge,       nz,          nCellX,     lgap,        nAsyType0,   &
                  nAsyType,    nCellType0,  nCellType,  nPinType0,   nPinType,    &
                  CellInfo,    AsyInfo,     PinInfo,    Core,                     &
                  CellPitch,   BaseCellInfo, nBaseCell, nBaseCell0
USE SetCentGeom, only : SetCentCell,  SetCentCellBase,   SetCentPin,    SetCentAsy, SetCentCore
IMPLICIT NONE

INTEGER :: icel, icel1, icel2, icel3
INTEGER :: iPin, ipin1, ipin2, ipin3
INTEGER :: iasy, iasy1, iasy2, iasy3
INTEGER :: iType

REAL :: HalfPitch

!Cent Symmetry Cell
DO icel = 1, nCellType0
  icel1 = nCellType0+icel; icel2 = icel1 + nCellType0; icel3 = icel2 + nCellType0
  iType = 1; CALL SetCentCell(icel, icel1, iType)
  iType = 2; CALL SetCentCell(icel, icel2, iType)
  iType = 3; CALL SetCentCell(icel, icel3, iType)
ENDDO

IF (nTracerCntl%lnTIGRst) THEN
  DO icel = 1, nBaseCell0
    icel1 = nBaseCell0+icel; icel2 = icel1 + nBaseCell0; icel3 = icel2 + nBaseCell0
    iType = 1; CALL SetCentCellBase(icel, icel1, iType)
    iType = 2; CALL SetCentCellBase(icel, icel2, iType)
    iType = 3; CALL SetCentCellBase(icel, icel3, iType)
  ENDDO
ENDIF

!Cent Symmetry Pin
DO ipin = 1, nPinType0
  ipin1= nPinType0+ipin; ipin2 = ipin1 + nPinType0; ipin3 = ipin2 + nPinType0
  iType = 1; CALL SetCentPin(PinInfo(ipin1), PinInfo(ipin), CellInfo, iType, ipin1)
  iType = 2; CALL SetCentPin(PinInfo(ipin2), PinInfo(ipin), CellInfo, iType, ipin2)
  iType = 3; CALL SetCentPin(PinInfo(ipin3), PinInfo(ipin), CellInfo, iType, ipin3)

  PinInfo(ipin1)%EdgePinIdx(0) = ipin;   PinInfo(ipin1)%EdgePinIdx(1) = ipin1
  PinInfo(ipin1)%EdgePinIdx(2) = ipin2;  PinInfo(ipin1)%EdgePinIdx(3) = ipin3
  PinInfo(ipin2)%EdgePinIdx(0) = ipin;   PinInfo(ipin2)%EdgePinIdx(1) = ipin1
  PinInfo(ipin2)%EdgePinIdx(2) = ipin2;  PinInfo(ipin2)%EdgePinIdx(3) = ipin3
  PinInfo(ipin3)%EdgePinIdx(0) = ipin;   PinInfo(ipin3)%EdgePinIdx(1) = ipin1
  PinInfo(ipin3)%EdgePinIdx(2) = ipin2;  PinInfo(ipin3)%EdgePinIdx(3) = ipin3
ENDDO

!Cent Symmetry Assemblies
DO iasy = 1, nAsyType0
  iasy1 = nAsyType0 + iasy; iasy2 = iasy1 + nAsyType0;  iasy3 = iasy2 + nAsyType0
  iType = 1; Call SetCentAsy(AsyInfo(iasy1), AsyInfo(iasy), PinInfo, iType, iasy1)
  iType = 2; Call SetCentAsy(AsyInfo(iasy2), AsyInfo(iasy), PinInfo, iType, iasy2)
  iType = 3; Call SetCentAsy(AsyInfo(iasy3), AsyInfo(iasy), PinInfo, iType, iasy3)

  AsyInfo(iasy1)%EdgeAsyIdx(0) = iasy;   AsyInfo(iasy1)%EdgeAsyIdx(1) = iasy1
  AsyInfo(iasy1)%EdgeAsyIdx(2) = iasy2;  AsyInfo(iasy1)%EdgeAsyIdx(3) = iasy3
  AsyInfo(iasy2)%EdgeAsyIdx(0) = iasy;   AsyInfo(iasy2)%EdgeAsyIdx(1) = iasy1;
  AsyInfo(iasy2)%EdgeAsyIdx(2) = iasy2;  AsyInfo(iasy2)%EdgeAsyIdx(3) = iasy3
  AsyInfo(iasy3)%EdgeAsyIdx(0) = iasy;   AsyInfo(iasy3)%EdgeAsyIdx(1) = iasy1;
  AsyInfo(iasy3)%EdgeAsyIdx(2) = iasy2;  AsyInfo(iasy3)%EdgeAsyIdx(3) = iasy3
ENDDO

CALL SetCentCore(Core, AsyInfo)

END SUBROUTINE


SUBROUTINE SethGapGeometries()
!Set up Half pincell and half assemblies
USE PARAM
USE GEOM,  only : lEdge,       nz,          nCellX,     lgap,        nAsyType0,   &
                  nAsyType,    nCellType0,  nCellType,  nPinType0,   nPinType,    &
                  CellInfo,    AsyInfo,     PinInfo,    Core,                     &
                  CellPitch
USE SetCentGeom, only : SetCentCell,     SetCentPin,    SetCentAsy, SetCentCore
IMPLICIT NONE

INTEGER :: icel, icel1, icel2, icel3
INTEGER :: iPin, ipin1, ipin2, ipin3
INTEGER :: iasy, iasy1, iasy2, iasy3
INTEGER :: iType
INTEGER :: nCellTypeQ

REAL :: HalfPitch

nCellTypeQ=nCellType0/4
!Cent Symmetry Cell
DO icel = 1, nCellTypeQ
  icel1 = nCellType0+icel; icel2 = icel1 + nCellType0; icel3 = icel2 + nCellType0
  iType = 1; CALL SetCentCell(icel, icel1, iType)
  iType = 2; CALL SetCentCell(icel, icel2, iType)
  iType = 3; CALL SetCentCell(icel, icel3, iType)
ENDDO

!Cent Symmetry Pin
DO ipin = 1, nPinType0
  ipin1= nPinType0+ipin; ipin2 = ipin1 + nPinType0; ipin3 = ipin2 + nPinType0
  iType = 1; CALL SetCentPin(PinInfo(ipin1), PinInfo(ipin), CellInfo, iType, ipin1)
  iType = 2; CALL SetCentPin(PinInfo(ipin2), PinInfo(ipin), CellInfo, iType, ipin2)
  iType = 3; CALL SetCentPin(PinInfo(ipin3), PinInfo(ipin), CellInfo, iType, ipin3)

  PinInfo(ipin1)%EdgePinIdx(0) = ipin;   PinInfo(ipin1)%EdgePinIdx(1) = ipin1
  PinInfo(ipin1)%EdgePinIdx(2) = ipin2;  PinInfo(ipin1)%EdgePinIdx(3) = ipin3
  PinInfo(ipin2)%EdgePinIdx(0) = ipin;   PinInfo(ipin2)%EdgePinIdx(1) = ipin1
  PinInfo(ipin2)%EdgePinIdx(2) = ipin2;  PinInfo(ipin2)%EdgePinIdx(3) = ipin3
  PinInfo(ipin3)%EdgePinIdx(0) = ipin;   PinInfo(ipin3)%EdgePinIdx(1) = ipin1
  PinInfo(ipin3)%EdgePinIdx(2) = ipin2;  PinInfo(ipin3)%EdgePinIdx(3) = ipin3
ENDDO

!Cent Symmetry Assemblies
DO iasy = 1, nAsyType0
  iasy1 = nAsyType0 + iasy; iasy2 = iasy1 + nAsyType0;  iasy3 = iasy2 + nAsyType0
  iType = 1; Call SetCentAsy(AsyInfo(iasy1), AsyInfo(iasy), PinInfo, iType, iasy1)
  iType = 2; Call SetCentAsy(AsyInfo(iasy2), AsyInfo(iasy), PinInfo, iType, iasy2)
  iType = 3; Call SetCentAsy(AsyInfo(iasy3), AsyInfo(iasy), PinInfo, iType, iasy3)

  AsyInfo(iasy1)%EdgeAsyIdx(0) = iasy;   AsyInfo(iasy1)%EdgeAsyIdx(1) = iasy1
  AsyInfo(iasy1)%EdgeAsyIdx(2) = iasy2;  AsyInfo(iasy1)%EdgeAsyIdx(3) = iasy3
  AsyInfo(iasy2)%EdgeAsyIdx(0) = iasy;   AsyInfo(iasy2)%EdgeAsyIdx(1) = iasy1;
  AsyInfo(iasy2)%EdgeAsyIdx(2) = iasy2;  AsyInfo(iasy2)%EdgeAsyIdx(3) = iasy3
  AsyInfo(iasy3)%EdgeAsyIdx(0) = iasy;   AsyInfo(iasy3)%EdgeAsyIdx(1) = iasy1;
  AsyInfo(iasy3)%EdgeAsyIdx(2) = iasy2;  AsyInfo(iasy3)%EdgeAsyIdx(3) = iasy3
ENDDO

CALL SetCentCore(Core, AsyInfo)

END SUBROUTINE
