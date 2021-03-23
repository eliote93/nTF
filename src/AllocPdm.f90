#include <defines.h>
SUBROUTINE AllocPDM()
!Allocate Problem Dependent Memomy
USE PARAM
USE CNTL,                 ONLY : nTracerCntl
USE Core_mod,             ONLY : GroupInfo
USE Geom,                 ONLY : ng
USE PE_MOD,               ONLY : PE
#ifdef MPI_ENV
USE MpiAxSolver_Mod,      ONLY : AllocMpiAxSol
!USE MPICOMM_Mod,          ONLY : MPI_SYNC
#endif
!USE AxSolver_mod, ONLY : AllocAxSol

IF(PE%RtMaster) CALL AllocFluxVariables()
CALL AllocXsVariables(nTracerCntl%lXsLib)
CALL AllocMocVariables()
IF(nTracerCntl%lCMFD .AND. PE%lCmfdGrp) THEN
  CALL AllocCMFD(TRUE)
  !IF(nTracerCntl%l3dim) CALL AllocAxSol(nTracerCntl%AxSolver)
  CALL AllocAxSol(nTracerCntl%AxSolver)
#ifdef MPI_ENV
  IF(PE%nCmfdProc .GT. 1) THEN
    CALL AllocMpiAxSol(nTracerCntl%AxSolver, GroupInfo, ng, PE)
  ENDIF
#endif
ENDIF

IF(nTracerCntl%lCritSpec) CALL Alloc_CritSpec()

IF(.NOT. nTracerCntl%lBenchXs) THEN
  CALL AllocTH()
ENDIF

IF(nTracerCntl%lDcpl) THEN
  CALL AllocDcplFmVariables()
  CALL AllocDcplCMFD(TRUE)
  CALL AllocDcplThInfo()
  IF(nTracerCntl%XsFtnMod .EQ. 2) CALL AllocDcplMicXsFtn()
ENDIF

IF(nTracerCntl%lGCCMFDGrp) THEN
  CALL Init_GCCMFD()
ENDIF
!Transient Init
IF(nTracerCntl%lProblem .EQ. lTransient) CALL AllocTransient()

END SUBROUTINE

SUBROUTINE AllocFluxVariables()

USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type
USE GEOM,         ONLY : Core,        ng,          nbd
USE Core_mod,     ONLY : Phis,        psi,         psid,     &
                         Psic,        Psicd,       Power,    &
                         PhiAngin,    RadJout,               &
                         LinSrcSlope, Phim,                  &
                         nCoreFsr,    nCoreFxr,   nPhiAngSv, &
                         CMInfo,      FmInfo,     wmoc,      &
                         !--- CNJ Edit : CASMO Linear Source
                         srcSlope,   phisSlope,   psiSlope,  &
                         !--- CNJ Edit : Domain Decomposition
                         AsyPhiAngIn
USE PE_MOD,       ONLY : PE
USE RAYS,         ONLY : RayInfo
USE ALLOCS
USE CNTL,         ONLY : nTracerCntl
USE SPH_mod,      ONLY : ssphf, ssphfnm
USE XSLIB_MOD,    ONLY : igresb,igrese
IMPLICIT NONE
INTEGER :: nFsr, nz, nzfm, nModRay, nxy, nxya, myzb, myze
INTEGER :: nPolar

nFsr = Core%nCoreFsr; nz = Core%nz
nxy = Core%nxy; nxya = Core%nxya
myzb = PE%myzb; myze = PE%myze
nPolar = RayInfo%nPolarAngle
nModRay = RayInfo%nModRay
nPhiAngSv = RayInfo%nPhiAngSv
nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr

CALL Dmalloc0(PHIS, 1, nFsr, myzb, myze, 1, ng)

CALL Dmalloc0(ssphf, 1, nFsr, myzb, myze, igresb,igrese)
ssphf=1._8

IF (nTracerCntl%lNodeMajor) THEN
  CALL Dmalloc0(ssphfnm, igresb, igrese, 1, nFsr, myzb, myze)
  ssphfnm = 1.0D0
ENDIF

!CALL Dmalloc0(Phis1g, 1, nFsr, myzb, myze)
CALL Dmalloc0(PhiAngin, 1, nPolar, 1, nPhiAngSv, myzb, myze, 1, ng)

!--- CNJ Edit : Domain Decomposition
IF (nTracerCntl%lDomainDcmp) THEN
  ALLOCATE(AsyPhiAngIn(nPolar, ng, 2, nModRay, nxya, myzb : myze))
  AsyPhiAngIn(:, :, :, :, :, :) = 1._8
ENDIF

!Current
!CALL Dmalloc0(RadJout, 1, 2, 1, nbd, 1, nxy, myzb, myze, 1, ng)
CALL Dmalloc0(RadJout, 1, 3, 1, nbd, 1, nxy, myzb, myze, 1, ng) !---BYS edit / 150612 Surface flux
CALL Dmalloc0(PSI, 1, nFsr, myzb, myze)
CALL Dmalloc0(PSID, 1, nFsr, myzb, myze)

CALL Dmalloc0(PSIC, 1, nxy, myzb, myze)
CALL Dmalloc0(PSICD, 1, nxy, myzb, myze)

CALL Dmalloc0(Power, 1, nFsr, 1, nz)
!Under-Relaxation
CALL Dmalloc(wmoc, ng)

!--- CNJ Edit : CASMO Linear Source
IF (nTracerCntl%lLinSrcCASMO) THEN
  CALL Dmalloc0(srcSlope, 1, 4, 1, ng, 1, nFsr, myzb, myze) ! 1:2 - 1/xst, 3:4 - 1/xst^2
  CALL Dmalloc0(phisSlope, 1, 2, 1, ng, 1, nFsr, myzb, myze)
  CALL Dmalloc0(psiSlope, 1, 2, 1, nFsr, myzb, myze)
ENDIF

!Linear Source Approximation
IF(nTracerCntl%lLinSrc) CALL Dmalloc0(LinSrcSlope, 1, 2, 1, nfsr, myzb, myze, 1, ng)

!--- CNJ Edit : Node Majors
IF(nTracerCntl%lScat1) THEN
  IF(nTracerCntl%SCatOd .EQ. 1) THEN
    IF (.NOT. nTracerCntl%lNodeMajor) THEN
      CALL Dmalloc0(phim, 1, 2, 1, nfsr, myzb, myze, 1, ng)
    ELSE
      CALL Dmalloc0(phim, 1, 2, 1, ng, 1, nfsr, myzb, myze)
    ENDIF
  ELSEIF (nTracerCntl%SCatOd .EQ. 2) THEN
    IF (.NOT. nTracerCntl%lNodeMajor) THEN
      CALL Dmalloc0(phim, 1, 5, 1, nfsr, myzb, myze, 1, ng)
    ELSE
      CALL Dmalloc0(phim, 1, 5, 1, ng, 1, nfsr, myzb, myze)
    ENDIF
  ELSEIF (nTracerCntl%SCatOd .EQ. 3) THEN
    IF (.NOT. nTracerCntl%lNodeMajor) THEN
      CALL Dmalloc0(phim, 1, 9, 1, nfsr, myzb, myze, 1, ng)
    ELSE
      CALL Dmalloc0(phim, 1, 9, 1, ng, 1, nfsr, myzb, myze)
    ENDIF
  ENDIF
ENDIF

IF(nTracerCntl%lDecusp) THEN
  ALLOCATE(FmInfo%neighphis(nFsr, ng, 2))
END IF

FmInfo%phis => Phis; FmInfo%PhiAngIn => PhiAngin
FMInfo%Psi => Psi;   FmInfo%PsiD => PsiD
FMInfo%PsiC => PsiC;   FmInfo%PsiCD => PsiCD
FMInfo%Power => Power
FMInfo%RadJout => RadJout !ninout/ 1:in 2:out 3:surfphi ! BYS edit 16//02/11
FmInfo%w => wmoc
RayInfo%RayInfo4Cmfd%PhiAngIn => PhiAngIn

!--- CNJ Edit : Domain Decomposition
FMInfo%AsyPhiAngIn => AsyPhiAngIn
RayInfo%RayInfo4Cmfd%AsyPhiAngIn => AsyPhiAngIn

IF(nTracerCntl%lLinSrc) FmInfo%LinSrcSlope => LinSrcSlope

IF(nTracerCntl%lScat1) FmInfo%phim => phim
END SUBROUTINE

SUBROUTINE AllocLinSrcVariables()
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type
USE GEOM,         ONLY : Core,          ng
USE Core_mod,     ONLY : FmInfo,        LinSrcSlope
USE PE_MOD,       ONLY : PE
USE ALLOCS
IMPLICIT NONE
INTEGER :: nFsr, nz, myzb, myze

nFsr = Core%nCoreFsr; nz = Core%nz
myzb = PE%myzb; myze = PE%myze


END SUBROUTINE

SUBROUTINE AllocXsVariables(lXsLib)
USE PARAM
USE TYPEDEF,    ONLY : XsMac_Type
USE BenchXs,    ONLY : MacXsBen
USE PE_MOD,     ONLY : PE
!USE LibXs,      ONLY : Fxr
USE GEOM,       ONLY : Core, ng !, myzb, myze
IMPLICIT NONE
Logical :: lxsLib
INTEGER :: nFxr, nz, myzb, myze

nFxr = Core%nCoreFxr; nz = Core%nz
myzb = PE%myzb; myze = PE%myze

!ALLOCATE(Fxr(nFxr, myzb:myze))

IF(lXsLib) THEN

ELSE

ENDIF
END SUBROUTINE

SUBROUTINE AllocMocVariables()
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type
USE GEOM,         ONLY : Core,     ng, nbd
USE PE_MOD,       ONLY : PE
USE Moc_mod,      ONLY : phis1g,     MocJout1g,   Xst1g,       tSrc,       &
                         AxSrc1g,    LinSrc1g,    LinPsi,      srcm,       &
                         AxPxs1g,    PhiAngIn1g,  phim1g,                  &
                         nMaxRaySeg, nMaxCellRay, nMaxCoreRay, nMaxAsyRay, &
                         SetMocEnvironmentVariables,                       &
                         !--- CNJ Edit : Node Majors
                         phisnm,     PhiAngInnm,  MocJoutnm,   xstnm,      &
                         srcnm,      srcmnm,      phimnm,                  &
                         !--- CNJ Edit : Domain Decomposition
                         DcmpPhiAngIn,    DcmpPhiAngOut
USE rays,         ONLY : RayInfo
USE Setray,       ONLY : RayInfoMaxSize
USE CNTL,         ONLY : nTracerCntl
USE ALLOCS
IMPLICIT NONE
INTEGER  :: nxy, nFsr, myzb, myze
INTEGER :: nPolar, nPhiAngSv

nxy = Core%nxy
nFsr = Core%nCoreFsr
nPolar = RayInfo%nPolarAngle
nPhiAngSv = RayInfo%nPhiAngSv
myzb = PE%myzb; myze = PE%myze

!CALL RayInfoMaxSize(Core, RayInfo, myzb, myze, nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay)
CALL SetMocEnvironmentVariables(nFsr, nxy, myzb, myze, ng)

CALL Dmalloc0(phis1g, 1, nFsr)
CALL Dmalloc0(PhiAngin1g, 1, nPolar, 1, nPhiAngSv)
!CALL Dmalloc0(MocJout1g, 1, 2, 1, nBd, 1, nxy)
CALL Dmalloc0(MocJout1g, 1, 3, 1, nBd, 1, nxy)  !---BYS edit / 150612 Surface flux
CALL Dmalloc0(xst1g, 1, nFsr)
CALL Dmalloc0(tsrc, 1, nFsr)
CALL Dmalloc0(AxSrc1g, 1, nxy)
CALL Dmalloc0(AxPxs1g, 1, nxy)
!Allocate Linear Source
IF(nTracerCntl%lLinSrc) Then
  CALL Dmalloc(LinSrc1g, nFsr, RayInfo%nAziAngle)
  CALL Dmalloc0(LinPsi, 1, 2, 1, nFsr, myzb, myze)
ENDIF

!--- CNJ Edit : Node Majors
IF(nTracerCntl%lscat1) THEN
  IF(nTracerCntl%lNodeMajor) THEN
    IF(nTracerCntl%ScatOd .EQ. 1) THEN
      CALL Dmalloc(srcmnm, 2, ng, nFsr)
    ELSEIF(nTracerCntl%ScatOd .EQ. 2) THEN
      CALL Dmalloc(srcmnm, 5, ng, nFsr)
    ELSEIF(nTracerCntl%ScatOd .EQ. 3) THEN
      CALL Dmalloc(srcmnm, 9, ng, nFsr)
    ENDIF
  ENDIF
  IF(nTracerCntl%ScatOd .EQ. 1) THEN
    CALL Dmalloc(srcm, 2, nFsr)
    CALL Dmalloc(phim1g, 2, nFsr)
  ELSEIF(nTracerCntl%ScatOd .EQ. 2) THEN
    CALL Dmalloc(srcm, 5, nFsr)
    CALL Dmalloc(phim1g, 5, nFsr)
  ELSEIF(nTracerCntl%ScatOd .EQ. 3) THEN
    CALL Dmalloc(srcm, 9, nFsr)
    CALL Dmalloc(phim1g, 9, nFsr)
  ENDIF
ENDIF

!--- CNJ Edit : Node Majors
IF (nTracerCntl%lNodeMajor) THEN
  CALL Dmalloc(phisnm, ng, nFsr)
  CALL Dmalloc(PhiAngInnm, nPolar, ng, nPhiAngSv)
  CALL Dmalloc(MocJoutnm, 3, ng, nbd, nxy)
  CALL Dmalloc(xstnm, ng, nFsr)
  CALL Dmalloc(srcnm, ng, nFsr)
ENDIF

!--- CNJ Edit : Domain Decomposition
IF (nTracerCntl%lDomainDcmp) THEN
  CALL Dmalloc(DcmpPhiAngOut, nPolar, ng, 2, RayInfo%nModRay, Core%nxya)
ENDIF

END SUBROUTINE

SUBROUTINE Alloc_CritSpec()
USE PARAM
USE TYPEDEF,      ONLY : CoreInfo_Type,   FmInfo_Type
USE GEOM,         ONLY : Core,            ng
USE Core_mod,     ONLY : FmInfo,          PhiCrit,  SpecConv
USE BasicOperation, ONLY : CP_CA
USE Allocs
IMPLICIT NONE
CALL Dmalloc(PhiCrit, ng)
CALL Dmalloc(SpecConv, ng)

CALL CP_CA(SpecConv(1:ng), 1._8, ng)

FmInfo%PhiCrit => PhiCrit
FmInfo%SpecConv => SpecConv

END SUBROUTINE
