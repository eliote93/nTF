#include <defines.h>
SUBROUTINE AllocPDM()
!Allocate Problem Dependent Memomy

USE PARAM,    ONLY : TRUE, mesg
USE CNTL,     ONLY : nTracerCntl
USE Core_mod, ONLY : GroupInfo
USE Geom,     ONLY : ng
USE PE_MOD,   ONLY : PE
USE files,    ONLY : io8
USE ioutil,   ONLY : message

#ifdef MPI_ENV
USE MpiAxSolver_Mod, ONLY : AllocMpiAxSol
#endif
! ----------------------------------------------------

mesg = 'Allocating Problem Dependent Memory...'
IF (PE%master) CALL message(io8, TRUE, TRUE, mesg)

IF (PE%RtMaster) CALL AllocFluxVariables

CALL AllocMocVariables

IF (nTracerCntl%lCMFD .AND. PE%lCmfdGrp) THEN
  CALL AllocCMFD(TRUE)
  CALL AllocAxSol(nTracerCntl%AxSolver)
  
#ifdef MPI_ENV
  IF (PE%nCmfdProc .GT. 1) CALL AllocMpiAxSol(nTracerCntl%AxSolver, GroupInfo, ng, PE)
#endif
END IF

IF (nTracerCntl%lCritSpec) CALL Alloc_CritSpec

IF (.NOT. nTracerCntl%lBenchXs) CALL AllocTH

IF (nTracerCntl%lDcpl) THEN
  CALL AllocDcplFmVariables
  CALL AllocDcplCMFD(TRUE)
  CALL AllocDcplThInfo
  
  IF (nTracerCntl%XsFtnMod .EQ. 2) CALL AllocDcplMicXsFtn
END IF

IF (nTracerCntl%lGCCMFDGrp) CALL Init_GCCMFD

IF (nTracerCntl%lProblem .EQ. lTransient) CALL AllocTransient
! ----------------------------------------------------

END SUBROUTINE AllocPDM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE AllocFluxVariables()

USE ALLOCS
USE PARAM,     ONLY : ONE
USE TYPEDEF,   ONLY : coreinfo_type
USE GEOM,      ONLY : Core, ng, nbd
USE Core_mod,  ONLY : Phis, psi, psid, Psic, Psicd, Power, PhiAngin, RadJout, LinSrcSlope, Phim, nCoreFsr, nCoreFxr, nPhiAngSv, CMInfo, FmInfo, wmoc, srcSlope, phisSlope, psiSlope, AsyPhiAngIn
USE PE_MOD,    ONLY : PE
USE RAYS,      ONLY : RayInfo, RayInfo4CMfd
USE CNTL,      ONLY : nTracerCntl
USE SPH_mod,   ONLY : ssphf, ssphfnm
USE XSLIB_MOD, ONLY : igresb, igrese

IMPLICIT NONE

INTEGER :: nFsr, nz, nModRay, nxy, nxya, myzb, myze, nPolar
! ----------------------------------------------------

nFsr     = Core%nCoreFsr
nz       = Core%nz
nxy      = Core%nxy
nxya     = Core%nxya
nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr

myzb = PE%myzb
myze = PE%myze

nPolar    = RayInfo%nPolarAngle
nModRay   = RayInfo%nModRay
nPhiAngSv = RayInfo%nPhiAngSv

! Spectral SPH
IF (nTracerCntl%lNodeMajor) THEN
  CALL dmalloc0(ssphfnm, igresb, igrese, 1, nFsr, myzb, myze); ssphfnm= ONE
ELSE
  CALL dmalloc0(ssphf, 1, nFsr, myzb, myze, igresb,igrese);    ssphf = ONE
END IF

! Current, Flux
CALL dmalloc0(RadJout, 1, 3, 1, nbd, 1, nxy, myzb, myze, 1, ng)

CALL dmalloc0(phis, 1, nFsr, myzb, myze, 1, ng)
CALL dmalloc0(psi,  1, nFsr, myzb, myze)
CALL dmalloc0(psid, 1, nFsr, myzb, myze)

CALL dmalloc0(psic,  1, nxy, myzb, myze)
CALL dmalloc0(psicd, 1, nxy, myzb, myze)

CALL dmalloc0(PhiAngin, 1, nPolar, 1, nPhiAngSv, myzb, myze, 1, ng)

IF (nTracerCntl%lDomainDcmp) THEN
  IF (nTracerCntl%lNodeMajor) THEN
    ALLOCATE (AsyPhiAngIn (nPolar, ng, 2, nModRay, nxya, myzb:myze))
  ELSE
    ALLOCATE (AsyPhiAngIn (nPolar, 2, nModRay, nxya, ng, myzb:myze))
  END IF
  
  AsyPhiAngIn = ONE
END IF

! Flux Moments
IF (nTracerCntl%lScat1) THEN
  SELECT CASE (nTracerCntl%SCatOd)
  CASE (1)
    IF (.NOT. nTracerCntl%lNodeMajor) THEN
      CALL dmalloc0(phim, 1, 2, 1, nfsr, myzb, myze, 1, ng)
    ELSE
      CALL dmalloc0(phim, 1, 2, 1, ng, 1, nfsr, myzb, myze)
    END IF
  CASE (2)
    IF (.NOT. nTracerCntl%lNodeMajor) THEN
      CALL dmalloc0(phim, 1, 5, 1, nfsr, myzb, myze, 1, ng)
    ELSE
      CALL dmalloc0(phim, 1, 5, 1, ng, 1, nfsr, myzb, myze)
    END IF
  CASE (3)
    IF (.NOT. nTracerCntl%lNodeMajor) THEN
      CALL dmalloc0(phim, 1, 9, 1, nfsr, myzb, myze, 1, ng)
    ELSE
      CALL dmalloc0(phim, 1, 9, 1, ng, 1, nfsr, myzb, myze)
    END IF
  END SELECT
END IF

! Power
CALL dmalloc0(Power, 1, nFsr, 1, nz)

! Under-Relaxation
CALL dmalloc(wmoc, ng)

! CASMO Linear Source
IF (nTracerCntl%lLinSrcCASMO) THEN
  CALL dmalloc0(srcSlope,  1, 4, 1, ng, 1, nFsr, myzb, myze) ! 1:2 - 1/xst, 3:4 - 1/xst^2
  CALL dmalloc0(phisSlope, 1, 2, 1, ng, 1, nFsr, myzb, myze)
  CALL dmalloc0(psiSlope,  1, 2,        1, nFsr, myzb, myze)
ENDIF

! Linear Source Approximation
IF (nTracerCntl%lLinSrc) CALL dmalloc0(LinSrcSlope, 1, 2, 1, nfsr, myzb, myze, 1, ng)

! Decusping
IF (nTracerCntl%lDecusp) ALLOCATE (FmInfo%neighphis (nFsr, ng, 2))

! Pointing
FmInfo%phis     => Phis
FmInfo%PhiAngIn => PhiAngin
FMInfo%Psi      => Psi
FmInfo%PsiD     => PsiD
FMInfo%PsiC     => PsiC
FmInfo%PsiCD    => PsiCD
FMInfo%Power    => Power
FMInfo%RadJout  => RadJout ! 1 : In-coming / 2 : Out-going / 3 : Surf. phi
FmInfo%w        => wmoc

IF (nTracerCntl%lLinSrc) FmInfo%LinSrcSlope => LinSrcSlope
IF (nTracerCntl%lScat1)  FmInfo%phim        => phim

RayInfo4CMfd%PhiAngIn => PhiAngIn

IF (nTracerCntl%lDomainDcmp) THEN
  FMInfo      %AsyPhiAngIn => AsyPhiAngIn
  RayInfo4Cmfd%AsyPhiAngIn => AsyPhiAngIn
END IF
! ----------------------------------------------------

END SUBROUTINE AllocFluxVariables
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE AllocMocVariables()

USE ALLOCS
USE GEOM,    ONLY : Core, ng, nbd
USE PE_MOD,  ONLY : PE
USE Moc_mod, ONLY : phis1g, MocJout1g, Xst1g, tSrc, AxSrc1g, LinSrc1g, LinPsi, srcm, AxPxs1g, PhiAngIn1g, &
                    phisnm, PhiAngInnm, MocJoutnm, xstnm, srcnm, srcmnm, DcmpPhiAngOutNg, DcmpPhiAngInNg, DcmpPhiAngOut1g, DcmpPhiAngIn1g
USE rays,    ONLY : RayInfo
USE CNTL,    ONLY : nTracerCntl

IMPLICIT NONE

INTEGER :: nxy, nFsr, myzb, myze, nPolar, nPhiAngSv
! ----------------------------------------------------

nxy  = Core%nxy
nFsr = Core%nCoreFsr

nPolar    = RayInfo%nPolarAngle
nPhiAngSv = RayInfo%nPhiAngSv

myzb = PE%myzb
myze = PE%myze

! Ray Tracing
CALL dmalloc(AxSrc1g, nxy)
CALL dmalloc(AxPxs1g, nxy)
  
IF (nTracerCntl%lNodeMajor) THEN
  CALL dmalloc(phisnm,     ng, nFsr)
  CALL dmalloc(PhiAngInnm, nPolar, ng, nPhiAngSv)
  CALL dmalloc(MocJoutnm,  3, ng, nbd, nxy)
  
  CALL dmalloc(xstnm, ng, nFsr)
  CALL dmalloc(srcnm, ng, nFsr)
ELSE
  CALL dmalloc(phis1g,     nFsr)
  CALL dmalloc(PhiAngin1g, nPolar, nPhiAngSv)
  CALL dmalloc(MocJout1g,  3, nBd, nxy)
  
  CALL dmalloc(xst1g, nFsr)
  CALL dmalloc(tsrc,  nFsr)
END IF

! Src. Moment
IF (nTracerCntl%lscat1) THEN
  IF (nTracerCntl%lNodeMajor) THEN
    SELECT CASE (nTracerCntl%ScatOd)
    CASE (1); CALL dmalloc(srcmnm, 2, ng, nFsr)
    CASE (2); CALL dmalloc(srcmnm, 5, ng, nFsr)
    CASE (3); CALL dmalloc(srcmnm, 9, ng, nFsr)
    END SELECT
  ELSE
    SELECT CASE (nTracerCntl%ScatOd)
    CASE (1); CALL dmalloc(srcm, 2, nFsr)
    CASE (2); CALL dmalloc(srcm, 5, nFsr)
    CASE (3); CALL dmalloc(srcm, 9, nFsr)
    END SELECT
  END IF
END IF

! Linear Source
IF (nTracerCntl%lLinSrc) Then
  CALL dmalloc (LinSrc1g, nFsr, RayInfo%nAziAngle)
  CALL dmalloc0(LinPsi, 1, 2, 1, nFsr, myzb, myze)
END IF

! Dcmp.
IF (nTracerCntl%lDomainDcmp) THEN
  IF (nTracerCntl%lNodeMajor) THEN
    CALL dmalloc(DcmpPhiAngOutNg, nPolar, ng, 2, RayInfo%nModRay, Core%nxya)
  ELSE
    CALL dmalloc(DcmpPhiAngOut1g, nPolar,     2, RayInfo%nModRay, Core%nxya)
  END IF
END IF
! ----------------------------------------------------

END SUBROUTINE AllocMocVariables
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE Alloc_CritSpec()

USE Allocs
USE GEOM,     ONLY : Core, ng
USE Core_mod, ONLY : FmInfo, PhiCrit, SpecConv


IMPLICIT NONE
! ----------------------------------------------------

CALL dmalloc(PhiCrit,  ng)
CALL dmalloc(SpecConv, ng)

SpecConv(1:ng) = 1.

FmInfo%PhiCrit  => PhiCrit
FmInfo%SpecConv => SpecConv
! ----------------------------------------------------

END SUBROUTINE Alloc_CritSpec
! ------------------------------------------------------------------------------------------------------------