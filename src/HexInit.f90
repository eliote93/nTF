#include <defines.h>
SUBROUTINE HexInit

USE ALLOCS
USE PARAM,         ONLY : mesg, TRUE
USE Geom,          ONLY : Core
USE RAYS,          ONLY : RayInfo, RayInfo4Cmfd, DcmpAsyRay
USE CNTL,          ONLY : nTracerCntl
USE PE_MOD,        ONLY : PE
USE files,         ONLY : io8
USE ioutil,        ONLY : message
USE MPIConfig_Mod, ONLY : SetGeomPEVariables, SetRayPEVariables, SetDcmpPEVariables
USE HexCmfd,       ONLY : HexRayInfo4CmfdGen
USE HexData,       ONLY : hRotRay, hcRay

IMPLICIT NONE

INTEGER :: iAng, nAziAng, nRotRay, irRay, icRay, jcRay, iAzm, nRay
! ----------------------------------------------------

! Basic
CALL SetGeomPEVariables(PE)

IF (nTracerCntl%lXsLib .AND. nTracerCntl%lrestrmt .AND. PE%RTmaster) CALL AllocResIsoInfo()

CALL PrepFxr(nTracerCntl%lXsLib, nTracerCntl%lfxrlib)

CALL SetDcmpPEVariables(PE)

IF (nTracerCntl%lDcpl) CALL initDcpl()

#ifdef MPI_ENV
CALL SetRayPEVariables(PE, Core, RayInfo)
#endif

IF (nTracerCntl%lDomainDcmp) THEN
#ifdef DetailedIO
  mesg = '  Generating Decomposed Assembly Rays...'
  IF (PE%master) CALL message(io8, TRUE, TRUE, mesg)
#endif
  CALL HexDcmpRayGen(Core, RayInfo, DcmpAsyRay)
ENDIF

CALL HexRayInfo4CmfdGen(RayInfo, RayInfo4Cmfd)

! INIT : NM
nAziAng = RayInfo%nAziAngle
nRotRay = RayInfo%nRotRay

CALL dmalloc0(RayInfo%RotRayAziList, 0, nRotRay, 1, nAziAng / 2)

nRay = nRotRay / (nAziAng / 2)

DO iAng = 1, nAziAng / 2 - 1
  RayInfo%RotRayAziList(0, iAng) = nRay
  
  DO irRay = 1, nRay
    RayInfo%RotRayAziList(irRay, iAng) = nRay*(iAng-1) + irRay
  END DO
END DO

RayInfo%RotRayAziList(0, nAziAng / 2) = nRotRay - RayInfo%RotRayAziList(nRay, nAziAng / 2 - 1)

DO irRay = 1, RayInfo%RotRayAziList(0, nAziAng / 2)
  RayInfo%RotRayAziList(irRay, iAng) = RayInfo%RotRayAziList(nRay, nAziAng / 2 - 1) + irRay
END DO

IF (.NOT. nTracerCntl%lMultigrid) THEN
  nTracerCntl%MultigridLV = 1
  nTracerCntl%gridNum     = 0
END IF

CALL MultigridRayGen(RayInfo)
! ----------------------------------------------------

END SUBROUTINE HexInit