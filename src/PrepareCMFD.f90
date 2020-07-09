#include <defines.h>
SUBROUTINE SetCMFDEnviorment(CORE, CmInfo, ng0,PE,lDcplPlnCmfd)
USE PARAM
USE geom,    ONLY : nbd
USE TYPEDEF, ONLY : PE_TYPE, CoreInfo_Type,  CmInfo_Type
USE CMFD_MOD, ONLY : nSubPlane,      SubPlaneMap,     SubPlaneRange,              &
                     hz,             hzfm,            ng,             nxy,        &
                     myzb,           myze,            myzbf,          myzef,      &
                     PinNeighIdx,    PinVol,          PinVolFm,                    &
                     CmfdPinXs,      CmfdLs,          AxDhat,         AxDtil,     &
                     AxPdhat

USE BASICOPERATION, ONLY : CP_CA
!USE PE_MOD, ONLY : PE
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(PE_TYPE) :: PE
INTEGER :: ng0
LOGICAL, OPTIONAL :: lDcplPlnCmfd
INTEGER :: ixy, ig
INTEGER :: nz, nzfm
INTEGER :: i, k

nz  = Core%nz; nzfm = Core%nzfm
ng  = ng0
nxy = Core%nxy

myzb = PE%myzb; myze = PE%myze
myzbf = PE%myzbf; myzef = PE%myzef
nz = CORE%nz; nzfm = CORE%nzfm
nSUbPlane = Core%nSubPlane
IF(.NOT. ASSOCIATED(SubPlaneMap)) THEN
  ALLOCATE(SubPlaneMap(nzfm))
  ALLOCATE(SubPlaneRange(2,nzfm))
ENDIF
IF(.NOT. present(lDcplPlnCmfd)) THEN
  SubPlaneMap(1:nzfm) = Core%SubPlaneMap(1:nzfm)
  SubPlaneRange(1:2,1:nz) = Core%SubPlaneRange(1:2,1:nz)
ELSEIF(lDcplPlnCmfd) THEN
  DO k = myzb, myze
    SubPlaneMap(k) = k;
    SubPlaneRange(:, k) = k
  ENDDO
ELSE
  SubPlaneMap(1:nz) = Core%SubPlaneMap(1:nz)
  SubPlaneRange(1:2,1:nzfm) = Core%SubPlaneRange(1:2,1:nzfm)
ENDIF
PinVol => Core%PinVol
PinVolFm => Core%PinVolFm
hz => Core%hz
hzfm => Core%hzfm
DO ixy = 1, nxy
  PinNeighIdx(1:nbd, ixy) = Core%Pin(ixy)%NeighIdx(1:nbd)
ENDDO

CmfdPinXS => CmInfo%PinXS
CmfdLS => CMInfo%CoreCmfdLs
AxDhat => CmInfo%AxDhat
AxDtil => CmInfo%AxDtil
AxPdhat => CmInfo%AxpDhat
CONTINUE
END SUBROUTINE

SUBROUTINE AllocCMFD(lCMFD)
USE PARAM
USE TYPEDEF,      ONLY : CoreInfo_Type, PinXs_Type, CMFDLs_Type, CmInfo_TYPE, &
                         RayInfo4CMFD_Type
!USE CNTL,   ONLY : lSubPlane
USE CNTL,         ONLY : nTracerCntl
USE PE_MOD,       ONLY : PE
USE GEOM,         ONLY : Core,   ng,       nbd
USE Core_Mod,     ONLY : PhiC,   PhiFm,    CoreCMFDLs,                       &
                         PsiC,   PsicD,    CMInfo,     RadJout,              &
                         PsiFm,  PsiFmD,   PinXs,      AxDtil,               &
                         AxDhat, AxPDhat,  GroupInfo
USE CMFD_Mod,     ONLY : XsMac,    PinNeighIdx,                              &
                         SRC,    Phic1g,                                     &
                         AllocPinXs
USE HexData,      ONLY : nInf
USE RAYS,         ONLY : RayInfo4Cmfd
USE BiCGSTAB_mod, ONLY : AllocBiCGSTAB
USE XsUtil_mod,   ONLY : AllocXsMac
USE Allocs
IMPLICIT NONE

LOGICAL :: lCMFD
INTEGER :: nxy, nEd
INTEGER :: i, ig
INTEGER :: ixy, iz
INTEGER :: myzb, myze, myzbf, myzef, nz, nFxrMax
INTEGER :: nCellType


nxy = Core%nxy; nCellType = Core%nCellType
nz = Core%nz
myzb = PE%myzb; myze = PE%myze
myzbf = PE%myzbf; myzef = PE%myzef

CALL Dmalloc0(Phic, 1, nxy, myzb-1, myze+1, 1, ng)
Allocate(PinXs(nxy, myzb : myze))

IF(nTracerCntl%lHex) THEN
  nEd = nInf
ELSE
  nEd = nCellType
ENDIF

CALL Dmalloc(PinNeighIdx, nbd, nxy)
!ALLOCATE CMFD Linear Systems
Allocate(CoreCMFDLs(ng))
DO ig = 1, ng
  CALL Dmalloc0(CoreCMFDLs(ig)%diag, 1, nxy, myzbf, myzef)
  CALL Dmalloc0(CoreCMFDLs(ig)%RadOffDiag, 1, nbd, 1, nxy, myzb, myze)
  CALL Dmalloc0(CoreCMFDLs(ig)%AxOffDiag, 1, 2, 1, nxy, myzbf, myzef)
!  IF(nTracerCntl%AxSolver .EQ. lP3SENM) THEN
!    CALL Dmalloc0(CoreCMFDLs(ig)%SP3AxDiag, 1, 4, 1, nxy, myzbf, myzef)
!    CALL Dmalloc0(CoreCMFDLs(ig)%SP3AxOffDiag, 1, 4, 1, 2, 1, nxy, myzbf, myzef)
!  ENDIF
  CoreCMFDLS(ig)%NeighIdx => PinNeighIdx
  CoreCMFDLS(ig)%myzbf = myzbf
  CoreCMFDLS(ig)%myzef = myzef
  CoreCMFDLS(ig)%nxy = nxy
  CoreCMFDLS(ig)%nbd = nbd
ENDDO

CALL SetCmfdMpiOmpEnv(Core, CoreCMFDLS, ng, PE)

!ALLOCATE HOMOGENIZED PinCell XS
DO iz = myzb, myze
  DO ixy = 1, nxy
    CALL AllocPinXs(PinXs(ixy, iz), ng, nbd, GroupInfo%InScatRange)
  ENDDO
ENDDO

CALL Dmalloc0(SRC, 1, nxy, myzbf, myzef)
CALL Dmalloc0(Phic1g, 1, nxy, myzbf - 1, myzef + 1)

IF(nTracerCntl%lSubPlane) THEN
  CALL Dmalloc0(PhiFm, 1, nxy, myzbf - 1, myzef + 1, 1, ng)
  CALL Dmalloc0(PsiFm, 1, nxy, myzbf, myzef)
  CALL Dmalloc0(PsiFmD, 1, nxy, myzbf, myzef)
ELSE
  PhiFm => PhiC
  PsiFm => PsiC
  PsiFmd => PsiCD
ENDIF


nFxrMax = 0
DO i = 1, nEd
  IF (.NOT. Core%CellInfo(i)%luse) CYCLE
  
  nFxrMax = MAX(nFxrMax, Core%CellInfo(i)%nFxr)
ENDDO
DO i = 1, nFxrMax
  XsMac(i)%ng = ng
  CALL AllocXsMac(XsMac(i))
ENDDO

!Allocate Axial Dhat
CALL Dmalloc0(AxDtil, 1, 2, 1, nxy, myzbf, myzef, 1, ng)
CALL Dmalloc0(AxDhat, 1, 2, 1, nxy, myzbf, myzef, 1, ng)
CALL Dmalloc0(AxPDhat, 1, 2, 1, nxy, myzbf, myzef, 1, ng)

CMInfo%PhiC => PhiC; CMInfo%PhiFm => PhiFm
CMInfo%PsiC => PsiC; CMInfo%PsiFm => PsiFm
CMInfo%PsiCD => PsiCD; CMInfo%PsiFmD => PsiFmD
CMInfo%RadJout => RadJout !ninout/ 1:in 2:out 3:surfphi ! BYS edit 16/02/11  > 16/02/18
CMInfo%CoreCMfdLs => CoreCMFDLs
CMInfo%PinXs => PinXS
CMInfo%RayInfo4Cmfd => RayInfo4Cmfd
CmInfo%AxDtil => AxDtil
CmInfo%AxDhat => AxDhat
CmInfo%AxPDhat => AxPDhat
CALL AllocBiCGSTAB(nxy, myzbf, myzef)
END SUBROUTINE

SUBROUTINE AllocPinXS(PinXS, ng, nbd, InScatRange)
USE PARAM
USE TYPEDEF, ONLY : PinXs_Type
USE ALLOCS
IMPLICIT NONE
TYPE(PinXs_Type) :: PinXs
INTEGER :: ng, nbd                    !Number of Groups and Number of Boundaries
INTEGER :: InScatRange(2, ng)
INTEGER :: ig, igb, ige

CALL Dmalloc(PinXS%Dtil, nbd, ng); CALL Dmalloc(PinXs%Dhat, nbd, ng)
CALL Dmalloc(PinXS%partialDhat, 2, nbd, ng)                            !--- CNJ Edit : p-CMFD Acceleration
CALL Dmalloc(PinXS%atil, nbd, ng); CALL Dmalloc(PinXS%ahat, nbd, ng)   !--- CNJ Edit : Domain Decomposition
CALL Dmalloc(PinXs%PDhat, nbd, ng)
CALL Dmalloc(PinXS%AxDtil, 2, ng); CALL Dmalloc(PinXs%AxDhat, 2, ng)
CALL Dmalloc(PinXS%Dtil2, nbd, ng)
CALL Dmalloc(PinXs%XSD, ng); CALL Dmalloc(PinXs%XSD2, ng); CALL Dmalloc(PinXs%XSA, ng);  !---BYS edit
CALL Dmalloc(PinXs%XST, ng); CALL Dmalloc(PinXs%XSTR, ng)
CALL Dmalloc(PinXs%XSR, ng); CALL Dmalloc(PinXs%XSNF, ng)
CALL Dmalloc(PinXs%XSKF, ng); CALL Dmalloc(PinXs%CHI, ng)
CALL Dmalloc(PInXs%phi, ng); CALL Dmalloc(PInXs%FAC, ng)

ALLOCATE(PinXS%XSS(ng))

DO ig = 1, ng
  igb = InScatRange(1, ig); ige = InScatRange(2, ig)
  CALL Dmalloc0(PinXs%Xss(ig)%From, igb, ige)
  PinXs%Xss(ig)%ib = igb; PinXs%Xss(ig)%ie = ige
ENDDO

END SUBROUTINE


SUBROUTINE SetCmfdMpiOMPEnv(Core, CMFDLS, ng, PE)
USE PARAM
USE TYPEDEF,       ONLY : CoreInfo_Type,   CMFDLS_TYPE,       PE_TYPE
USE LsRadDcmp_MOD, ONLY : SetRadDcmp
USE CNTL,          ONLY : nTracerCntl
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmfdLs_Type) :: CMFDLS(ng)
TYPE(PE_TYPE) :: PE
INTEGER :: ng

INTEGER :: ig

DO ig = 1, ng
#ifdef MPI_ENV
  CmfdLs(ig)%nProc = PE%nCmfdProc
  CmfdLs(ig)%comm = PE%MPI_CMFD_COMM
  CmfdLs(ig)%myrank = PE%myCmfdrank
#endif
  CmfdLs(ig)%nThread = PE%nCmfdThread
  !CmfdLs(ig)%nThread = 4
ENDDO

IF((PE%nCmfdThread .GT. 1).AND.(.NOT.nTracerCntl%lHex)) CALL SetRadDcmp(Core, PE)
!PRINT *, 'Environment Setting'
END SUBROUTINE
