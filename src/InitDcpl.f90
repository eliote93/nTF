#include <defines.h>
SUBROUTINE InitDcpl()
USE PARAM
USE MPIDcpl_Mod,    ONLY : InitDcplParallel
IMPLICIT NONE

CALL InitDcplPlnInfo()
CALL InitDcplCntl() 
!Set MPI Environment for Decoupled calculation
CALL InitDcplParallel()
!Set Decoupled Fxr Information
CALL DcplPrepFxr()
END SUBROUTINE

SUBROUTINE InitDcplPlnInfo()
USE PARAM
USE GEOM,           ONLY : Core
USE DcplCore_Mod,   ONLY : DcplInfo
USE Allocs
IMPLICIT NONE
INTEGER :: nz, nRefPln
INTEGER :: iz, iRefPln
INTEGER :: i, j
nz = Core%nz; nRefPln = DcplInfo%nRefPln
CALL Dmalloc(DcplInfo%Pln_map0, nz)
DO iz = 1, nz
  iRefPln = DcplInfo%Pln_map(iz)
  DO j = 1, nRefPln
    IF(iRefPln .EQ. DcplInfo%RefPln(j)) EXIT
  ENDDO
  IF(j .GT. nRefPln)  THEN
    PRINT *, 'Input Error for Ref_Map'
    STOP
  ENDIF
  DcplInfo%Pln_map0(iz) = j
ENDDO
CALL Dmalloc(DcplINfo%RefPln_map, 2, nRefPln)
DcplINfo%RefPln_map(1, :) = nz;DcplINfo%RefPln_map(2, :) = 1
DO iz = 1, nz
  j = DcplInfo%Pln_map0(iz)
  DcplInfo%RefPln_map(1, j) = MIN(DcplInfo%RefPln_map(1, j), iz)
  DcplInfo%RefPln_map(2, j) = MAX(DcplInfo%RefPln_map(2, j), iz)
ENDDO

END SUBROUTINE


SUBROUTINE InitDcplCntl()
USE PARAM
USE DcplCore_mod,  ONLY : DcplInfo 
USE Cntl,          ONLY : DcplControl    ,nTracerCntl          ,CopyCntlOption
USE ItrCntl_mod,   ONLY : ItrCntl        ,DcplItrCntl                              &
                         ,ItrCntl_TYPE   ,DcplXsGenCntl_Type   ,CMFDItrCntl_TYPE
USE PE_MOD,        ONLY : PE             ,DcplPE

IMPLICIT NONE
INTEGER :: irefpln, nrefpln

CALL CopyCntlOption(nTracerCntl, DcplControl)
DcplControl%l3dim = FALSE
DcplControl%lFeedBack = FALSE
DcplControl%lDcplCal = TRUE
nRefPln = DcplInfo%nRefPln
!nTracerCntl%lXsFtn = .FALSE.
DO iRefPln = 1, nRefPln
  DcplPE(iRefPln)%myzb = DcplInfo%RefPln(iRefPln)
  DcplPE(iRefPln)%myze = DcplPE(iRefPln)%myzb
  DcplPE(iRefPln)%myzbf = DcplPE(iRefPln)%myzb
  DcplPE(iRefPln)%myzef = DcplPE(iRefPln)%myze
ENDDO
!Iteration Control
IF(.NOT. nTracerCntl%lfeedback) THEN
  ItrCntl%DcplCMFD3dItrCntl%nIterMax = ItrCntl%DcplItrData(1, 2)
  ItrCntl%DcplCMFD3dItrCntl%nIterMin = ItrCntl%DcplItrData(2, 2)
ELSE
  ItrCntl%DcplCMFD3dItrCntl%nIterMax = ItrCntl%DcplItrData(1, 3)
  ItrCntl%DcplCMFD3dItrCntl%nIterMin = ItrCntl%DcplItrData(2, 3)
  IF(nTracerCntl%lXsFtn) THEN
    ItrCntl%nThCondiGen = ItrCntl%DcplItrData(3, 5)
    DO iRefPln = 1, nRefPln
      DcplItrCntl(iRefPln)%nThCondiGen = ItrCntl%DcplItrData(3, 5)
      DcplItrCntl(iRefPln)%DcplXsGenCntl%nXsGenFbIterMax = ItrCntl%DcplItrData(0:3, 1)
    ENDDO  
  ELSE
    DO iRefPln = 1, nRefPln
      DcplItrCntl(iRefPln)%DcplXsGenCntl%nXsGenIterMax = ItrCntl%DcplItrData(0, 1)
    ENDDO
  ENDIF
ENDIF

END SUBROUTINE

SUBROUTINE SetDcplGeom()
USE PARAM
USE TypeDef,        ONLY : CoreInfo_Type     ,Pin_Type      ,Asy_Type          &
                          ,AsyInfo_type      ,pininfo_type  ,cell_type         &
                          ,PE_Type
USE GEOM,           ONLY : Core  
USE PE_MOD,         ONLY : PE
USE DcplCore_Mod,   ONLY : DcplInfo          ,DcplCore
IMPLICIT NONE
INTEGER :: nRefPln
INTEGER :: myRefPlnBeg, myRefPlnEnd
INTEGER :: iRefPln
INTEGER :: iz

nRefPln = DcplInfo%nRefPln
#ifdef MPI_ENV
myRefPlnBeg = PE%myRefPlnBeg; myRefPlnEnd = PE%myRefPlnEnd
#else
myRefPlnBeg = 1; myRefPlnEnd = nRefPln
#endif
ALLOCATE(DcplCore(nRefPln))

DO iRefPln = myRefPlnBeg, myRefPlnEnd
  iz = DcplInfo%RefPln(iRefPln)
  Core%lDcpl = TRUE
  CALL CopyCommonCoreInfo(DcplCore(iRefPln), Core)
  CALL CopyPlnwiseCoreInfo(DcplCore(iRefPln), Core, iz)
  CALL SetDcplPin(DcplCore(iRefPln), Core, iz)
  CALL SetDcplPinInfo(DcplCore(iRefPln), Core, iz)
ENDDO
END SUBROUTINE

SUBROUTINE CopyCommonCoreInfo(DcplCore, Core)
USE PARAM
USE TypeDef,        ONLY : CoreInfo_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: DcplCore, Core

!Pin & Assembly Numer Info
DcplCore%nxya = Core%nxya; DcplCore%nxa = Core%nxa; DcplCore%nya = Core%nya
DcplCore%nxy = Core%nxy; DcplCore%nx = Core%nx; DcplCore%ny = Core%ny
DcplCore%nxy0 = Core%nxy0; DcplCore%nx0 = Core%nx0; DcplCore%ny0 = Core%ny0

!FxrInfo
DcplCore%nCoreFxr = Core%nCoreFxr; DcplCore%nCoreFsr = Core%nCoreFsr

!Radial Boundary Info
DcplCore%RadSym = Core%RadSym; DcplCore%RadBC = Core%RadBC
DcplCore%AxBC =  RefCell

DcplCore%nz = 1; DcplCore%nzfm= 1; DcplCore%nSubPlane= 1
DcplCore%nAsyCell = Core%nAsyCell; DcplCore%nAsyGT = Core%nAsyGT

DcplCore%nPinType = Core%nPinType; DcplCore%nCellType = Core%nCellType
DcplCore%nAsyType = Core%nAsyType


DcplCore%CoreIdx => Core%CoreIdx
DcplCore%CoreMap => Core%CoreMap

DcplCore%Asy => Core%Asy
!DcplCore%Pin => Core%Pin

!CellInfo Pointing
DcplCore%CellInfo => Core%CellInfo
DcplCore%AsyInfo => Core%AsyInfo
DcplCore%AsyCentX => Core%AsyCentX
DcplCore%AsyCentY => Core%AsyCentY
END SUBROUTINE

SUBROUTINE CopyPlnwiseCoreInfo(DcplCore, Core, iz)
USE PARAM
USE TypeDef,        ONLY : CoreInfo_Type
USE ALLOCS
IMPLICIT NONE
TYPE(CoreInfo_Type) :: DcplCore, Core

INTEGER :: iz

DcplCore%iRefPln = iz 
!DcplCore%
DcplCore%lFuelPlane => Core%lFuelPlane(iz:iz)
CALL Dmalloc(DcplCore%hz, 1); CALL Dmalloc(DcplCore%hzInv, 1)
DcplCore%hz(1) = Core%hz(iz); DcplCore%HzInv(1) = 1._8 / Core%hz(iz)
DcplCore%HzFm => DcplCore%hz; DcplCore%HzFmInv => DcplCore%hzInv

DcplCore%PinVol => Core%PinVol(:, iz:iz)
DcplCore%PinVolFm => Core%PinVol

DcplCore%AsyVol => Core%AsyVol(:, iz:iz)

END SUBROUTINE

SUBROUTINE SetDcplPinInfo(DcplCore, Core, iz)
USE PARAM
USE TypeDef,        ONLY : CoreInfo_Type,    pininfo_type
USE ALLOCS
IMPLICIT NONE
TYPE(CoreInfo_Type) :: DcplCore, Core
TYPE(PinInfo_Type), POINTER :: RefPinInfo, myPinInfo
INTEGER :: iz, itype, ipin
INTEGER :: nPinType
!PinInfo => Core%Pin
nPinType = Core%nPinType
ALLOCATE(DcplCore%PinInfo(nPinType))
DO ipin = 1, nPinType
  RefPinInfo => Core%PinInfo(ipin)
  myPinInfo => DcplCore%PinInfo(ipin)
  myPinInfo%lempty = RefPinInfo%lempty; myPinInfo%luse = RefPinInfo%luse
  myPinInfo%lfuel = RefPinInfo%lfuel
  myPinInfo%lCentX = RefPinInfo%lCentX; myPinInfo%lCentY = RefPinInfo%lCentY
  myPinInfo%lCentXY = RefPinInfo%lCentXY
  myPinInfo%nCell = 1
  myPinInfo%nFsrMax = RefPinInfo%nFsrMax
  myPinInfo%nFxrMax = RefPinInfo%nFxrMax
  myPinInfo%PartialAsyFlag = RefPinInfo%PartialAsyFlag
  myPinInfo%EdgePinIdx = RefPinInfo%EdgePinIdx
  CALL Dmalloc0(myPinInfo%Cell, iz, iz)
  myPinInfo%Cell(iz) = RefPinInfo%Cell(iz)
ENDDO
END SUBROUTINE

SUBROUTINE SetDcplPin(DcplCore, Core, iz)
USE PARAM
USE TypeDef,        ONLY : CoreInfo_Type,    pin_type
USE ALLOCS
IMPLICIT NONE
TYPE(CoreInfo_Type) :: DcplCore, Core
TYPE(Pin_Type), POINTER :: RefPin, myPin
INTEGER :: iz, itype, ixy
INTEGER :: nxy
!PinInfo => Core%Pin
nxy = Core%nxy
ALLOCATE(DcplCore%Pin(nxy))
DO ixy = 1, nxy
  mypin => DcplCore%Pin(ixy); refpin => Core%Pin(ixy)
  MyPin%lfuel = RefPin%lfuel
  MyPin%PinType = RefPin%PinType; MyPin%AsyType = RefPin%AsyType
  MyPin%nCell = 1
  MyPin%nFsrMax = RefPin%nFsrMax; MyPin%nFxrMax = RefPin%nFxrMax
  MyPin%iasy = RefPin%iasy; MyPin%ipin = RefPin%ipin
  MyPin%nBd = RefPin%nBd;
  MyPin%FsrIdxSt = RefPin%FsrIdxSt
  MyPin%FxrIdxSt = RefPin%FxrIdxSt
  CALL Dmalloc0(MyPin%Cell, iz,iz); MyPin%Cell(iz) = RefPin%Cell(iz)
  CALL Dmalloc(MyPin%BDLength, MyPin%nBd); CALL Dmalloc(MyPin%Center2SurfaceL, MyPin%nBd)
  CALL Dmalloc(MyPin%NeighIdx, MyPin%nBd); CALL Dmalloc(MyPin%NeighSurfIdx, MyPin%nBd)
  MyPin%BDLength = RefPin%BdLength; MyPin%Center2SurfaceL = RefPin%Center2SurfaceL
  MyPin%NeighIdx = RefPin%NeighIdx; MyPin%NeighSurfIdx = RefPin%NeighSurfIdx
ENDDO

END SUBROUTINE

SUBROUTINE AllocDcplMicXsFtn()
USE PARAM
USE TYPEDEF,       ONLY : CoreInfo_Type     ,MicXsFtn_Type,     PE_TYPE
USE GEOM,          ONLY : Core              ,ng
USE Core_mod,      ONLY : GroupInfo
USE DcplCore_mod,  ONLY : DcplInfo          ,MicXsFtn
USE PE_MOD,        ONLY : PE
IMPLICIT NONE
INTEGER :: myRefPlnBeg, myRefPlnEnd, myzb, myze, nxy
INTEGER :: nreftemp
INTEGER :: irefpln,ixy
myzb = PE%myzb; myze = PE%myze
nxy = Core%nxy
nRefTemp = DcplInfo%nRefTemp
#ifdef MPI_ENV
myRefPlnBeg = PE%myRefPlnBeg; myRefPlnEnd = PE%myRefPlnEnd
#else
myRefPlnBeg = myzb; myRefPlnEnd = PE%myze
#endif
 !Alloc_MicXsFtn(MicXsFtn, ng, GroupInfo)
ALLOCATE(MicXsFtn(nxy, myRefPlnBeg : myRefPlnEnd))
DO irefpln = myRefPlnBeg, myRefPlnEnd
  DO ixy = 1, nxy
    CALL Alloc_MicXsFtn(MicXsFtn(ixy, irefpln), ng, GroupInfo)
  ENDDO
  !CALL Alloc_MicXsFtn(MicXsFtn(
ENDDO
DcplInfo%MicXsFtn => MicXsFtn
END SUBROUTINE

SUBROUTINE AllocDcplFmVariables()
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type     ,FmInfo_Type      ,DcplInfo_Type
USE GEOM,           ONLY : Core              ,ng
USE DcplCore_Mod,   ONLY : DcplInfo          ,DcplFmInfo       ,DcplFxr
USE RAYS,           ONLY : RayInfo
USE ALLOCS
USE PE_MOD,         ONLY : PE
USE CNTL,           ONLY : nTracerCntl
IMPLICIT NONE

INTEGER :: nRefPln, nRefTemp
INTEGER :: nFsr, nz, nzfm, nCoreRay, nxy, nbd, myzb, myze
INTEGER :: nPolar, nPhiAngSv
INTEGER :: irefpln, irefTemp
INTEGER :: myRefPlnBeg, myRefPlnEnd
TYPE(FmInfo_Type), POINTER :: FmInfo
nRefPln = DcplInfo%nRefPln
nRefTemp = DcplInfo%nRefTemp
nPolar = RayInfo%nPolarAngle
nPhiAngSv = RayInfo%nPhiAngSv
#ifdef MPI_ENV
myRefPlnBeg = PE%myRefPlnBeg; myRefPlnEnd = PE%myRefPlnEnd
ALLOCATE(DcplFmInfo(nRefTemp, myRefPlnBeg:myRefPlnEnd))
#else
myRefPlnBeg = 1; myRefPlnEnd = nRefPln
ALLOCATE(DcplFmInfo(nRefTemp, nRefPln))
#endif

DO irefpln = myRefPlnBeg, myRefPlnEnd
  DO irefTemp = 1, nRefTemp
    FmInfo => DcplFmInfo(iRefTemp, iRefPln)
    FmInfo%Fxr => DcplFxr(iRefPln)%Fxr
    nFsr = Core%nCoreFsr; nz = 1
    nxy = Core%nxy
    myzb = DcplInfo%RefPln(irefpln); myze = DcplInfo%RefPln(irefpln)
    nbd = 4 ! For rectangular geometry
    CALL Dmalloc0(FmInfo%Phis, 1, nFsr, myzb, myze, 1, ng)
    CALL Dmalloc0(FmInfo%PhiAngIn, 1, nPolar, 1, nPhiAngSv, myzb, myze, 1, ng)
    !CALL Dmalloc0(FmInfo%RadJout, 1, 2, 1, nbd, 1, nxy, myzb, myze, 1, ng)
    CALL Dmalloc0(FmInfo%RadJout, 1, 3, 1, nbd, 1, nxy, myzb, myze, 1, ng)  !---BYS edit / 150612 Surface flux
    CALL Dmalloc0(FmInfo%Psi, 1, nFsr, myzb, myze); CALL Dmalloc0(FmInfo%PsiD, 1, nFsr, myzb, myze)
    CALL Dmalloc0(FmInfo%PsiC, 1, nxy, myzb, myze); CALL Dmalloc0(FmInfo%PsiCD, 1, nxy, myzb, myze)
    CALL Dmalloc0(FmInfo%AxSrc, 1, nxy, myzb, myze, 1, ng); CALL Dmalloc0(FmInfo%AxPXS, 1, nxy, myzb, myze, 1, ng)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE AllocDcplCMFD(lCMFD)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type      ,PinXs_Type       ,CMFDLs_Type            &
                          ,DcplCMFDLS_Type    ,CmInfo_TYPE      ,RayInfo4CMFD_Type
USE CNTL,           ONLY : nTracerCntl
USE GEOM,           ONLY : Core              ,ng
USE DcplCore_Mod,   ONLY : DcplInfo          ,DcplCore         ,DcplFmInfo              &
                          ,DcplCmInfo        ,DcplCmfdLs
USE CORE_Mod,       ONLY : GroupInfo
USE CMFD_Mod,       ONLY : PinNeighIdx       ,AllocPinXs
USE RAYS,           ONLY : RayInfo4Cmfd
USE PE_MOD,         ONLY : PE
USE ALLOCS 
LOGICAL :: lCMFD
INTEGER :: nxy, nbd
INTEGER :: i, ig
INTEGER :: iRefPln,iRefTemp
INTEGER :: ixy, iz
INTEGER :: myzb, myze, myzbf, myzef, nFxrMax
INTEGER :: myRefPlnBeg, myRefPlnEnd
INTEGER :: nRefTemp, nRefPln
INTEGER :: nCellType

TYPE(CmInfo_Type), POINTER :: CmInfo
TYPE(CMFDLs_Type), POINTER :: CMFDLs(:)
nRefTemp = DcplInfo%nRefTemp; nRefPln = DcplInfo%nRefPln
#ifndef MPI_ENV
myRefPlnBeg = 1; myRefPlnEnd = nRefPln
#else
myRefPlnBeg = PE%myRefPlnBeg; myRefPlnEnd = PE%myRefPlnEnd
#endif
ALLOCATE(DcplCmInfo(nRefTemp, myRefPlnBeg:myRefPlnEnd))
ALLOCATE(DcplCmfdLs(nRefTemp, myRefPlnBeg:myRefPlnEnd))
nxy  = Core%nxy; nbd = 4
DO iRefPln =  myRefPlnBeg, myRefPlnEnd
  myzb = DcplInfo%RefPln(iRefPln); myze = myzb
  DO iRefTemp = 1, nRefTemp
    ALLOCATE(DcplCmfdLs(iRefTemp, iRefPln)%CMFDLS(ng))
    CMFDLs => DcplCmfdLs(iRefTemp, iRefPln)%CMFDLS
    DO ig = 1, ng
      CALL Dmalloc0(CMFDLs(ig)%diag, 1, nxy, myzb, myze)
      CALL Dmalloc0(CMFDLs(ig)%RadOffDiag, 1, nbd, 1, nxy, myzb, myze)
      CALL Dmalloc0(CMFDLs(ig)%AxOffDiag, 1, 2, 1, nxy, myzb, myze)
      CMFDLS(ig)%NeighIdx => PinNeighIdx
      CMFDLS(ig)%myzbf = myzb; CMFDLS(ig)%myzef = myze
      CMFDLS(ig)%nxy = nxy; CMFDLS(ig)%nbd = nbd      
    ENDDO
    
    CmInfo => DcplCmInfo(iRefTemp, iRefPln)
    CMInfo%CoreCMfdLs => CmfdLS
    
    CALL Dmalloc0(CmInfo%PhiC, 1, nxy, myzb-1, myze+1, 1, ng)
    CmInfo%PhiFm => CmInfo%PhiC
    CmInfo%PsiC => DcplFmInfo(iRefTemp, iRefPln)%PsiC
    CmInfo%PsiCD => DcplFmInfo(iRefTemp, iRefPln)%PsiCD
    CmInfo%PsiFM => DcplFmInfo(iRefTemp, iRefPln)%PsiC
    CmInfo%PsiFMD => DcplFmInfo(iRefTemp, iRefPln)%PsiCD
    
    CmInfo%RadJout => DcplFmInfo(iRefTemp, iRefPln)%RadJout
    ALLOCATE(CmInfo%PinXs(nxy, myzb:myze))

#ifdef DcplMemSave
    IF(iRefTemp .EQ. 1) THEN
      DO ixy = 1, nxy
        CALL AllocPinXs(CmInfo%PinXs(ixy, myzb), ng, nbd, GroupInfo%InScatRange)
        CALL Dmalloc(CmInfo%PinXs(ixy, myzb)%Dcpl_Dhat, nbd, ng, nRefPln)
      ENDDO
    ELSE
      CmInfo%PinXs => DcplCmInfo(1, iRefPln)%PinXS
    ENDIF
#else
    DO ixy = 1, nxy
      CALL AllocPinXs(CmInfo%PinXs(ixy, myzb), ng, nbd, GroupInfo%InScatRange)
    ENDDO
#endif    

    CMInfo%RayInfo4Cmfd  => RayInfo4Cmfd

  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE AllocDcplThInfo()
USE PARAM
USE TypeDef,         ONLY : ThInfo_Type
USE Geom,            ONLY : Core
USE DcplCore_mod,    ONLY : DcplInfo        ,DcplThInfo
USE CNTL,            ONLY : nTracerCntl
USE PE_MOD,          ONLY : PE
USE ALLOCS
IMPLICIT NONE
INTEGER :: nRefPln, nRefTemp, myzb, myze
INTEGER :: iz, nz, nFxr, nxy
INTEGER :: myRefPlnBeg, myRefPlnEnd

nRefPln = DcplInfo%nRefPln; nRefTemp = DcplInfo%nRefTemp
#ifndef MPI_ENV
myRefPlnBeg = 1; myRefPlnEnd = nRefPln
#else
myRefPlnBeg = PE%myRefPlnBeg; myRefPlnEnd = PE%myRefPlnEnd
#endif
nz = Core%nz; nfxr = Core%nCoreFxr
nxy = Core%nxy
ALLOCATE(DcplThInfo(myRefPlnBeg : myRefPlnEnd))
DO iz = myRefPlnBeg, myRefPlnEnd
  !DcplThInfo(iz)
  myzb = DcplInfo%RefPln(iz); myze = myzb
  CALL dmalloc0(DcplThInfo(iz)%RefFuelTemp, 0, nz)
  !CALL dmalloc0(DcplThInfo(iz)%FxrTemp, 1, nfxr, myzb, myze)
  IF(nTracerCntl%lFeedBack) THEN
    CALL Dmalloc0(DcplThInfo(iz)%RelPower, myzb, myze, 1, nxy)
    CALL Dmalloc0(DcplThInfo(iz)%Tcool, myzb, myze, 1, nxy)
    CALL Dmalloc0(DcplThInfo(iz)%TcoolInOut, 1, 2, myzb, myze, 1, nxy)
  ENDIF 
ENDDO
IF(.NOT. nTracerCntl%lFeedBack) RETURN
END SUBROUTINE

SUBROUTINE Alloc_MicXsFtn(MicXsFtn, ng, GroupInfo)
USE PARAM
USE TYPEDEF,  ONLY : MicXsFtn_TYPE           ,GroupInfo_Type
USE ALLOCS
IMPLICIT NONE
TYPE(MicXsFtn_Type) :: MicXsFtn
TYPE(GroupInfo_Type) :: GroupInfo
INTEGER :: ng
INTEGER :: igb,ige, ig, i, J
CALL Dmalloc0(MicXsFtn%sig_T, 1, ng, 1, 4, 0, 1)
CALL Dmalloc0(MicXsFtn%sig_Nf, 1, ng, 1, 4, 0, 1)
CALL Dmalloc0(MicXsFtn%sig_kf, 1, ng, 1, 4, 0, 1)

ALLOCATE(MicXsFtn%XSS(ng, 1:4, 0:1))
DO j = 0, 1
  DO i = 1, 4
    DO ig = 1, ng
      igb = GroupInfo%InScatRange(1, ig); ige = GroupInfo%InScatRange(2, ig)
      CALL Dmalloc0(MicXsFtn%Xss(ig, i, j)%From, igb, ige)
      MicXsFtn%Xss(ig, i, j)%ib = igb; MicXsFtn%Xss(ig, i, j)%ie = ige
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE
