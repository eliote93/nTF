#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetCMFDEnviorment(CORE, CmInfo, ng0, PE, lDcplPlnCmfd)

USE geom,     ONLY : nbd
USE TYPEDEF,  ONLY : PE_TYPE, CoreInfo_Type, CmInfo_Type
USE CMFD_MOD, ONLY : nSubPlanePtr, nSubPlane, SubPlaneMap, SubPlaneRange, hz, hzfm, ng, nxy, myzb, myze, myzbf, myzef, &
                     PinNeighIdx, PinVol, PinVolFm, CmfdPinXs, CmfdLs, AxDhat, AxDtil, AxPdhat

IMPLICIT NONE

TYPE (CoreInfo_Type) :: Core
TYPE (CmInfo_Type)   :: CmInfo
TYPE (PE_TYPE)       :: PE

INTEGER :: ng0
LOGICAL, OPTIONAL :: lDcplPlnCmfd

INTEGER :: ixy, ig, nz, nzfm, i, k
! ----------------------------------------------------

nz   = Core%nz
nzfm = Core%nzfm
nxy  = Core%nxy
nz   = CORE%nz
nzfm = CORE%nzfm

ng = ng0

myzb  = PE%myzb
myze  = PE%myze
myzbf = PE%myzbf
myzef = PE%myzef

! Pointing
PinVol   => Core%PinVol
PinVolFm => Core%PinVolFm
hz       => Core%hz
hzfm     => Core%hzfm

CmfdPinXS => CmInfo%PinXS
CmfdLS    => CMInfo%CoreCmfdLs
AxDhat    => CmInfo%AxDhat
AxDtil    => CmInfo%AxDtil
AxPdhat   => CmInfo%AxpDhat

DO ixy = 1, nxy
  PinNeighIdx(1:nbd, ixy) = Core%Pin(ixy)%NeighIdx(1:nbd)
END DO

! Sub-Pln.
nSUbPlane     = Core%nSubPlane
nSubPlanePtr => Core%nSubPlanePtr

IF (.NOT. ASSOCIATED(SubPlaneMap)) THEN
  ALLOCATE (SubPlaneMap(nzfm))
  ALLOCATE (SubPlaneRange(2,nzfm))
END IF

IF(.NOT. present(lDcplPlnCmfd)) THEN
  SubPlaneMap  (1:nzfm)    = Core%SubPlaneMap  (1:nzfm)
  SubPlaneRange(1:2, 1:nz) = Core%SubPlaneRange(1:2, 1:nz)
ELSE IF(lDcplPlnCmfd) THEN
  DO k = myzb, myze
    SubPlaneMap     (k) = k
    SubPlaneRange(:, k) = k
  END DO
ELSE
  SubPlaneMap  (1:nz)        = Core%SubPlaneMap  (1:nz)
  SubPlaneRange(1:2, 1:nzfm) = Core%SubPlaneRange(1:2, 1:nzfm)
END IF
! ----------------------------------------------------

END SUBROUTINE SetCMFDEnviorment
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE AllocCMFD(lCMFD)

USE Allocs
USE PARAM
USE TYPEDEF,      ONLY : CoreInfo_Type, PinXs_Type, CMFDLs_Type, CmInfo_TYPE, RayInfo4CMFD_Type
USE CNTL,         ONLY : nTracerCntl
USE PE_MOD,       ONLY : PE
USE GEOM,         ONLY : Core, ng, nbd
USE Core_Mod,     ONLY : PhiC, PhiFm, CoreCMFDLs, PsiC, PsicD, CMInfo, RadJout, PsiFm, PsiFmD, PinXs, AxDtil, AxDhat, AxPDhat, GroupInfo
USE CMFD_Mod,     ONLY : XsMac, PinNeighIdx, SRC, Phic1g, AllocPinXs
USE HexData,      ONLY : ncTyp
USE RAYS,         ONLY : RayInfo4Cmfd
USE BiCGSTAB_mod, ONLY : AllocBiCGSTAB
USE XsUtil_mod,   ONLY : AllocXsMac

IMPLICIT NONE

LOGICAL :: lCMFD
INTEGER :: nxy, nEd, i, ig, ixy, iz, myzb, myze, myzbf, myzef, nz, nFxrMax, nCellType
! ----------------------------------------------------

nxy       = Core%nxy
nCellType = Core%nCellType
nz        = Core%nz

myzb  = PE%myzb
myze  = PE%myze
myzbf = PE%myzbf
myzef = PE%myzef

! Flux
CALL dmalloc0(Phic,   1, nxy, myzb-1,  myze+1, 1, ng)
CALL dmalloc0(SRC,    1, nxy, myzbf,   myzef)
CALL dmalloc0(Phic1g, 1, nxy, myzbf-1, myzef+1)

IF (nTracerCntl%lSubPlane) THEN
  CALL dmalloc0(PhiFm,  1, nxy, myzbf-1, myzef+1, 1, ng)
  CALL dmalloc0(PsiFm,  1, nxy, myzbf,   myzef)
  CALL dmalloc0(PsiFmD, 1, nxy, myzbf,   myzef)
ELSE
  PhiFm  => PhiC
  PsiFm  => PsiC
  PsiFmd => PsiCD
ENDIF

! Pin XS
ALLOCATE (PinXs(nxy, myzb:myze))

DO iz = myzb, myze
  DO ixy = 1, nxy
    CALL AllocPinXs(PinXs(ixy, iz), ng, nbd, GroupInfo%InScatRange)
  END DO
END DO

CMInfo%PhiC         => PhiC
CMInfo%PhiFm        => PhiFm
CMInfo%PsiC         => PsiC
CMInfo%PsiFm        => PsiFm
CMInfo%PsiCD        => PsiCD
CMInfo%PsiFmD       => PsiFmD
CMInfo%RadJout      => RadJout ! 1 : In-coming / 2 : Out-going / 3 : Surf. phi
CMInfo%PinXs        => PinXS
CMInfo%RayInfo4Cmfd => RayInfo4Cmfd

IF (PE%lMKL) RETURN ! Does not used in MKL CMFD

! LS
CALL dmalloc(PinNeighIdx, nbd, nxy)

ALLOCATE (CoreCMFDLs (ng))

DO ig = 1, ng
  CALL Dmalloc0(CoreCMFDLs(ig)%diag,               1, nxy, myzbf, myzef)
  CALL Dmalloc0(CoreCMFDLs(ig)%RadOffDiag, 1, nbd, 1, nxy, myzb,  myze)
  CALL Dmalloc0(CoreCMFDLs(ig)%AxOffDiag,  1, 2,   1, nxy, myzbf, myzef)
  
  CoreCMFDLS(ig)%NeighIdx => PinNeighIdx
  CoreCMFDLS(ig)%myzbf     = myzbf
  CoreCMFDLS(ig)%myzef     = myzef
  CoreCMFDLS(ig)%nxy       = nxy
  CoreCMFDLS(ig)%nbd       = nbd
END DO

! Env. of MPI & OMP
CALL SetCmfdMpiOmpEnv(Core, CoreCMFDLS, ng, PE)

! Temporary XS
IF (nTracerCntl%lHex) THEN
  nEd = ncTyp
ELSE
  nEd = nCellType
ENDIF

nFxrMax = 0
DO i = 1, nEd
  IF (.NOT. Core%CellInfo(i)%luse) CYCLE
  
  nFxrMax = MAX(nFxrMax, Core%CellInfo(i)%nFxr)
END DO

DO i = 1, nFxrMax
  XsMac(i)%ng = ng
  
  CALL AllocXsMac(XsMac(i))
END DO

! Axial Dhat
CALL dmalloc0(AxDtil,  1, 2, 1, nxy, myzbf, myzef, 1, ng)
CALL dmalloc0(AxDhat,  1, 2, 1, nxy, myzbf, myzef, 1, ng)
CALL dmalloc0(AxPDhat, 1, 2, 1, nxy, myzbf, myzef, 1, ng)

! Pointing
CMInfo%CoreCMfdLs   => CoreCMFDLs
CmInfo%AxDtil       => AxDtil
CmInfo%AxDhat       => AxDhat
CmInfo%AxPDhat      => AxPDhat

! BiCGSTAB
CALL AllocBiCGSTAB(nxy, myzbf, myzef)
! ----------------------------------------------------

END SUBROUTINE AllocCMFD
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE AllocPinXS(PinXS, ng, nbd, InScatRange)

USE ALLOCS
USE TYPEDEF, ONLY : PinXs_Type
USE CNTL,    ONLY : nTracerCntl

IMPLICIT NONE

TYPE (PinXs_Type) :: PinXs

INTEGER :: ng, nbd, ig, igb, ige, InScatRange(2, ng)
! ----------------------------------------------------

CALL dmalloc(PinXS%Dtil,  nbd, ng)
CALL dmalloc(PinXs%Dhat,  nbd, ng)
CALL dmalloc(PinXs%PDhat, nbd, ng)
CALL dmalloc(PinXS%AxDtil,  2, ng)
CALL dmalloc(PinXs%AxDhat,  2, ng)
CALL dmalloc(PinXS%Dtil2, nbd, ng)

CALL dmalloc(PinXs%XSD,  ng)
CALL dmalloc(PinXs%XSD2, ng)
CALL dmalloc(PinXs%XSA,  ng)
CALL dmalloc(PinXs%XST,  ng)
CALL dmalloc(PinXs%XSTR, ng)
CALL dmalloc(PinXs%XSR,  ng)
CALL dmalloc(PinXs%XSNF, ng)
CALL dmalloc(PinXs%XSKF, ng)
CALL dmalloc(PinXs%CHI,  ng)
CALL dmalloc(PInXs%phi,  ng)
CALL dmalloc(PInXs%FAC,  ng)

IF (nTracerCntl%lDomainDcmp) THEN
  CALL dmalloc(PinXS%atil, nbd, ng)
  CALL dmalloc(PinXS%ahat, nbd, ng)
END IF

CALL dmalloc(PinXS%partialDhat, 2, nbd, ng)

ALLOCATE (PinXS%XSS (ng))

DO ig = 1, ng
  igb = InScatRange(1, ig)
  ige = InScatRange(2, ig)
  
  CALL dmalloc0(PinXs%Xss(ig)%From, igb, ige)
  
  PinXs%Xss(ig)%ib = igb
  PinXs%Xss(ig)%ie = ige
END DO
! ----------------------------------------------------

END SUBROUTINE AllocPinXS
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetCmfdMpiOMPEnv(Core, CMFDLS, ng, PE)

USE TYPEDEF,       ONLY : CoreInfo_Type, CMFDLS_TYPE, PE_TYPE
USE LsRadDcmp_MOD, ONLY : SetRadDcmp
USE CNTL,          ONLY : nTracerCntl

IMPLICIT NONE

TYPE (CoreInfo_Type) :: Core
TYPE (CmfdLs_Type) :: CMFDLS(ng)
TYPE (PE_TYPE) :: PE

INTEGER :: ng, ig
! ----------------------------------------------------

DO ig = 1, ng
#ifdef MPI_ENV
  CmfdLs(ig)%nProc  = PE%nCmfdProc
  CmfdLs(ig)%comm   = PE%MPI_CMFD_COMM
  CmfdLs(ig)%myrank = PE%myCmfdrank
#endif

  CmfdLs(ig)%nThread = PE%nCmfdThread
END DO

IF ((PE%nCmfdThread .GT. 1).AND.(.NOT.nTracerCntl%lHex)) CALL SetRadDcmp(Core, PE)
! ----------------------------------------------------

END SUBROUTINE SetCmfdMpiOMPEnv
! ------------------------------------------------------------------------------------------------------------