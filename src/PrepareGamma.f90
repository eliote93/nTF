#include <defines.h>
#ifdef __GAMMA_TRANSPORT    
!-- JSU EDIT 20170804
! SUBROUTINES FOR GENERAL GAMMA DATA
SUBROUTINE PrepareGamma(Core, RayInfo)
USE TYPEDEF,        ONLY : CoreInfo_Type,       RayInfo_Type
! subroutine to prepare problem dependent photon data
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo

CALL InitGammaMOC(Core, RayInfo)
CALL AllocGammaCMFD(Core)

END SUBROUTINE

!***************************************************************************************************
! SUBROUTINES FOR MOC OF GAMMA CALCULATION
SUBROUTINE InitGammaMOC(Core, RayInfo)
USE TYPEDEF,        ONLY : CoreInfo_Type,       RayInfo_Type,       PolarAngle_Type,    AziAngleInfo_Type
USE GammaCore_mod
USE GamMOC_MOD
USE MOC_MOD,        ONLY : AziMap
USE ALLOCS
USE CNTL,           ONLY : nTracerCntl
USE PE_MOD,         ONLY : PE
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo

TYPE(AziAngleInfo_Type), POINTER :: AziAng(:)
TYPE(PolarAngle_Type), POINTER :: PolarAng(:)
INTEGER :: nPolarAngle, nAziAngle, nThread
INTEGER :: ipol, iazi
INTEGER :: tid, od, ScatOd
REAL :: wttemp, wtcos, wtpolar, wtsin2

AziAng => RayInfo%AziAngle
PolarAng => RayInfo%PolarAngle
nAziAngle = RayInfo%nAziAngle
nPolarAngle = RayInfo%nPolarAngle
nThread = PE%nThread
ScatOd = nTracerCntl%GammaScatOd

IF (ScatOd .EQ. 1) od = 2
IF (ScatOd .EQ. 2) od = 5
IF (ScatOd .EQ. 3) od = 9

CALL AllocGammaPDM(Core, RayInfo)

DO ipol = 1, nPolarAngle
  wttemp = PolarAng(ipol)%weight * PolarAng(ipol)%sinv
  DO iazi = 1, nAziAngle
    wtang(ipol, iazi) = wttemp * AziAng(iazi)%weight * AziAng(iazi)%del
    wtsurf(ipol, iazi, 1) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / abs(AziAng(iazi)%sinv)
    wtsurf(ipol, iazi, 3) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / abs(AziAng(iazi)%sinv)
    wtsurf(ipol, iazi, 2) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / abs(AziAng(iazi)%cosv)
    wtsurf(ipol, iazi, 4) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / abs(AziAng(iazi)%cosv)
  ENDDO
ENDDO

IF (nTracerCntl%lGammaScat1) THEN
  DO iAzi = 1, nAziAngle / 2
    AziMap(iAzi, 1) = 1
    AziMap(iAzi, 2) = 2
    AziMap(nAziAngle - iAzi + 1, 1) = 3
    AziMap(nAziAngle - iAzi + 1, 2) = 4
  ENDDO
  DO ipol = 1, nPolarAngle
    wttemp = PolarAng(ipol)%sinv
    DO iazi = 1, nAziAngle
      comp(1, ipol, iazi) = wttemp * AziAng(iazi)%cosv
      comp(2, ipol, iazi) = wttemp * AziAng(iazi)%sinv
      mwt(1:2, ipol, iazi) = comp(1:2, ipol, iazi) *  wtang(ipol, iazi)
    ENDDO
  ENDDO
  IF (ScatOd .GE. 2) THEN
    DO ipol = 1, nPolarAngle
      wttemp = PolarAng(ipol)%sinv
      wtsin2 = PolarAng(ipol)%sinv * PolarAng(ipol)%sinv
      wtcos = PolarAng(ipol)%cosv
      wtpolar = 1.5 * PolarAng(ipol)%cosv * PolarAng(ipol)%cosv - 0.5
      DO iazi = 1, nAziAngle
        Comp(3, ipol, iazi) = wtpolar
        Comp(4, ipol, iazi) = wtsin2 * (1.0 - 2.0 * AziAng(iazi)%sinv * AziAng(iazi)%sinv)
        Comp(5, ipol, iazi) = wtsin2 * (2.0 * AziAng(iazi)%sinv * AziAng(iazi)%cosv)
        mwt(3, ipol, iazi) = comp(3, ipol, iazi) *  wtang(ipol, iazi)
        mwt(4:5, ipol, iazi) = 0.75 * comp(4:5, ipol, iazi) *  wtang(ipol, iazi)
      ENDDO
    ENDDO
  ENDIF
  IF (ScatOd .EQ. 3) THEN
    DO ipol = 1, nPolarAngle
      wttemp = PolarAng(ipol)%sinv
      DO iazi = 1, nAziAngle
        Comp(6, ipol, iazi) = (5.0 * PolarAng(ipol)%cosv * PolarAng(ipol)%cosv - 1.0) * wttemp * AziAng(iazi)%cosv
        Comp(7, ipol, iazi) = (5.0 * PolarAng(ipol)%cosv * PolarAng(ipol)%cosv - 1.0) * wttemp * AziAng(iazi)%sinv
        Comp(8, ipol, iazi) = (wttemp ** 3.0) * (4.0 * (AziAng(iazi)%cosv ** 3.0) - 3.0 * AziAng(iazi)%cosv)
        Comp(9, ipol, iazi) = (wttemp ** 3.0) * (- 4.0 * (AziAng(iazi)%sinv ** 3.0) + 3.0 * AziAng(iazi)%sinv)
        mwt(6:7, ipol, iazi) = 0.375 * comp(6:7, ipol, iazi) * wtang(ipol, iazi)
        mwt(8:9, ipol, iazi) = 0.625 * comp(8:9, ipol, iazi) * wtang(ipol, iazi)
      ENDDO
    ENDDO
  ENDIF
ENDIF

END SUBROUTINE

SUBROUTINE AllocGammaPDM(Core, RayInfo)
USE TYPEDEF,        ONLY : CoreInfo_Type,       RayInfo_Type
USE GammaCore_mod
USE GamMOC_MOD
USE MOC_MOD,        ONLY : AziMap
USE ALLOCS
USE CNTL,           ONLY : nTracerCntl
USE PE_MOD,         ONLY : PE
USE Core_mod,       ONLY : GroupInfo
USE GEOM,           ONLY : nbd
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo

INTEGER :: ng, nFsr, nxy, nz, nPolarAngle, nAziAngle, nPhiAngSv, nThread
INTEGER :: myzb, myze
INTEGER :: tid, od

ng = GroupInfo%ngg
nFsr = Core%nCoreFsr
nxy = Core%nxy
nz = Core%nz
myzb = PE%myzb
myze = PE%myze
nAziAngle = RayInfo%nAziAngle
nPolarAngle = RayInfo%nPolarAngle
nPhiAngSv = RayInfo%nPhiAngSv
nThread = PE%nThread

IF (nTracerCntl%GammaScatOd .EQ. 1) od = 2
IF (nTracerCntl%GammaScatOd .EQ. 2) od = 5
IF (nTracerCntl%GammaScatOd .EQ. 3) od = 9

!--- Allocate Core Varaibles

IF (nTracerCntl%lGammaScat1) THEN
  CALL Dmalloc0(gphis, 1, nFsr, 1, ng, myzb, myze)
  CALL Dmalloc0(gphim, 1, od, 1, nFsr, 1, ng, myzb, myze)
  CALL Dmalloc0(gPhiAngIn, 1, nPolarAngle, 1, nPhiAngSv, 1, ng, myzb, myze)
  CALL Dmalloc0(gJout, 1, 3, 1, nbd, 1, nxy, 1, ng, myzb, myze)
ELSE
  CALL Dmalloc0(gphis, 1, ng, 1, nFsr, myzb, myze)
  CALL Dmalloc0(gPhiAngIn, 1, nPolarAngle, 1, ng, 1, nPhiAngSv, myzb, myze)
  CALL Dmalloc0(gJout, 1, 3, 1, ng, 1, nbd, 1, nxy, myzb, myze)
ENDIF
CALL Dmalloc0(gphic, 1, nxy, 1, ng, myzb, myze)
CALL Dmalloc0(gAxSrc, 1, nxy, 1, ng, myzb, myze)
CALL Dmalloc0(gAxPXS, 1, nxy, 1, ng, myzb, myze)

CALL Dmalloc0(gPower, 1, nFSR, 1, nz)
CALL Dmalloc0(LocalNPower, 1, nFSR, 1, nz)  !-- JSU EDIT 2017.09.15.
CALL Dmalloc0(LocalNPower_KERMA, 1, nFSR, 1, nz)  !-- JSU EDIT 2017.09.15.
CALL Dmalloc0(GPowerGen, 1, nFSR, 1, nz)  !-- JSU EDIT 2017.09.15.

gphis = 1E-10

!--- Allocate MOC Variables

IF (nTracerCntl%lGammaScat1) THEN
  CALL Dmalloc(gxst1g, nFsr)
  CALL Dmalloc(gsrc1g, nFsr)
  CALL Dmalloc(gsrcm1g, od, nFsr)
  CALL Dmalloc(Comp, od, nPolarAngle, nAziAngle)
  CALL Dmalloc(mwt, od, nPolarAngle, nAziAngle)
ELSE
  CALL Dmalloc(gxstnm, ng, nFsr)
  CALL Dmalloc(gsrcnm, ng, nFsr)
ENDIF

CALL Dmalloc(wtang, nPolarAngle, nAziAngle)
CALL Dmalloc(wtsurf, nPolarAngle, nAziAngle, 4)

ALLOCATE(TrackingDat(nThread))

IF (nTracerCntl%lGammaScat1) THEN
  DO tid = 1, nThread
    CALL Dmalloc(TrackingDat(tid)%phia1g, nPolarAngle, nFsr, 4)
    CALL Dmalloc(TrackingDat(tid)%SrcAng1g, nPolarAngle, nFsr, 4)
    CALL Dmalloc(TrackingDat(tid)%Jout1g, 3, nbd, nxy)
    TrackingDat(tid)%wtang => wtang; TrackingDat(tid)%wtsurf => wtsurf
    TrackingDat(tid)%AziMap => AziMap
  ENDDO
ELSE
  DO tid = 1, nThread
    CALL Dmalloc(TrackingDat(tid)%phisnm, ng, nFsr)
    CALL Dmalloc(TrackingDat(tid)%Joutnm, 3, ng, nbd, nxy)
    TrackingDat(tid)%wtang => wtang; TrackingDat(tid)%wtsurf => wtsurf
  ENDDO
ENDIF

END SUBROUTINE

!***************************************************************************************************
! SUBROUTINES FOR CMFD OF GAMMA CALCULATION
SUBROUTINE AllocGammaCMFD(Core)
! subroutine to allocate problem dependent CMFD data
USE PARAM
USE TYPEDEF,            ONLY : CoreInfo_Type
USE GamCMFD_MOD,        ONLY : CMFD_gPinXS,         GammaCMFD,          AllocgPinXS
USE GammaCore_mod,      ONLY : PinHeatXS
USE Core_MOD,           ONLY : GroupInfo
USE PE_MOD,             ONLY : PE
USE CNTL,               ONLY : nTRACERCntl
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
INTEGER :: ng, ngg, nxy, nz
INTEGER :: myzb, myze
INTEGER :: ipin, iz

ng = GroupInfo%ng
ngg = GroupInfo%ngg
nxy = Core%nxy
myzb = PE%myzb
myze = PE%myze
nz = myze - myzb + 1

ALLOCATE(CMFD_gPinXS(nxy, myzb : myze))

DO iz = myzb, myze
  DO ipin = 1, nxy
    CALL AllocgPinXS(CMFD_gPinXS(ipin, iz), ng, ngg, 4, GroupInfo%InScatRange_PH)
  ENDDO
ENDDO
IF (nTRACERCntl%lExplicitKappa) THEN
  ALLOCATE(PinHeatXS(nxy, myzb : myze))
  DO iz = myzb, myze
    DO ipin = 1, nxy
      CALL AllocPin_HeatXS(PinHeatXS(ipin, iz), ng, ngg)
    END DO
  END DO
END IF

ALLOCATE(GammaCMFD%M(ngg))
ALLOCATE(GammaCMFD%ILU(ngg))
ALLOCATE(GammaCMFD%Scat(nxy, ngg, ngg))
ALLOCATE(GammaCMFD%Prod(nxy, ng, ngg))
ALLOCATE(GammaCMFD%Phic(nxy, myzb : myze, ng))
ALLOCATE(GammaCMFD%gPhic(nxy, myzb : myze, ngg))
ALLOCATE(GammaCMFD%gProd(nxy * nz, ngg))
ALLOCATE(GammaCMFD%gSrc(nxy * nz))

END SUBROUTINE

SUBROUTINE AllocGPINXS(PINXS, ng, ngg, nbd, InScatRange)
! subroutine to allocate gamma pin-wise homogenized cross section
USE GammaTYPEDEF,     ONLY : GPINXS_TYPE
IMPLICIT NONE
! INPUT VARIABLES
TYPE(GPINXS_TYPE) :: PINXS
INTEGER :: ng, ngg, nbd
INTEGER :: InScatRange(2,ngg)

!LOCAL VARIABLES
INTEGER :: ig, igb, ige, igg

ALLOCATE(PINXS%NPHI(ng), PINXS%GPHI(ngg))  ! FLUX
ALLOCATE(PINXS%XST(ngg), PINXS%XSTR(ngg))  ! TOTAL & TRANSPORT
!ALLOCATE(PINXS%KERMA(ngg))  ! KERMA
ALLOCATE(PINXS%XSA(ngg), PINXS%XSD(ngg), PINXS%XSR(ngg)) ! CROSS SECTIONS
ALLOCATE(PINXS%XSP(ng,ngg)) ! PRODUCTION MATRIX
ALLOCATE(PINXS%DTIL(nbd,ngg), PINXS%DHAT(nbd,ngg), PINXS%PDHAT(nbd,ngg))

PINXS%NPHI = 0.;  PINXS%GPHI = 0.
!PINXS%KERMA = 0.
PINXS%XST = 0.;  PINXS%XSTR = 0.;
PINXS%XSA = 0.;  PINXS%XSD  = 0.;  PINXS%XSR = 0.
PINXS%XSP  = 0.
PINXS%DTIL = 0.; PINXS%DHAT = 0.; PINXS%PDHAT = 0.

ALLOCATE(PINXS%XSS(ngg))

DO igg = 1, ngg
  igb = InScatRange(1, igg); ige = InScatRange(2, igg)
  ALLOCATE(PINXS%XSS(igg)%from(igb:ige))
  PINXS%XSS(igg)%from = 0.
  PINXS%XSS(igg)%ib = igb
  PINXS%XSS(igg)%ie = ige
END DO

END SUBROUTINE

SUBROUTINE AllocPIN_HEATXS(PINXS, ng, ngg)
USE GammaTYPEDEF,       ONLY : PINHEATXS_TYPE
IMPLICIT NONE
TYPE(PINHEATXS_TYPE) :: PINXS
INTEGER :: ng, ngg

INTEGER :: ig, igg

ALLOCATE(PINXS%KERMA_T(ng))
ALLOCATE(PINXS%KERMA_F(ng))
ALLOCATE(PINXS%PhotonGen(ng)) ! Total photon energy generated by neutron reaction
ALLOCATE(PINXS%kFis(ng))      ! Fission Energy Release (kappa sigma fission)
ALLOCATE(PINXS%NPHI(ng), PINXS%GPHI(ngg))  ! FLUX

ALLOCATE(PINXS%KERMA_P(ngg)) ! photoatomic KERMA

PINXS%KERMA_T = 0.;PINXS%KERMA_F = 0.;
PINXS%PhotonGen = 0.;PINXS%kFis = 0.;
PINXS%NPHI = 0.;PINXS%GPHI = 0.;
PINXS%KERMA_P = 0.;

END SUBROUTINE
  
#endif
