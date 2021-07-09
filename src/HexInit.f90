#include <defines.h>
SUBROUTINE HexInit

USE PARAM,         ONLY : mesg, TRUE
USE Geom,          ONLY : Core
USE RAYS,          ONLY : RayInfo, RayInfo4Cmfd, DcmpAsyRay
USE CNTL,          ONLY : nTracerCntl
USE PE_MOD,        ONLY : PE
USE files,         ONLY : io8
USE ioutil,        ONLY : message
USE MPIConfig_Mod, ONLY : SetGeomPEVariables, SetRayPEVariables, SetDcmpPEVariables
USE setray,        ONLY : HexDcmpRayGen
USE HexCmfd,       ONLY : HexRayInfo4CmfdGen

IMPLICIT NONE
! ----------------------------------------------------

! PE
CALL SetGeomPEVariables(PE)
CALL SetDcmpPEVariables(PE)

#ifdef MPI_ENV
CALL SetRayPEVariables(PE, Core, RayInfo)
#endif

! XS
mesg = 'Allocating Resonant Isotope Information...'
IF (PE%master) CALL message(io8, TRUE, TRUE, mesg)
  
IF (nTracerCntl%lXsLib .AND. nTracerCntl%lrestrmt .AND. PE%RTmaster) CALL AllocResIsoInfo

mesg = 'Preparing FXR...'
IF (PE%master) CALL message(io8, TRUE, TRUE, mesg)

CALL PrepFxr(nTracerCntl%lXsLib, nTracerCntl%lfxrlib)

! Dcpl
IF (nTracerCntl%lDcpl) CALL initDcpl

! Dcmp
IF (nTracerCntl%lDomainDcmp) THEN
  mesg = 'Generating Decomposed Assembly Rays...'
  IF (PE%master) CALL message(io8, TRUE, TRUE, mesg)
  
  CALL HexDcmpRayGen(Core, RayInfo, DcmpAsyRay)
END IF

! Ray 4 CMFD
CALL HexRayInfo4CmfdGen(RayInfo, RayInfo4Cmfd)
! ----------------------------------------------------

END SUBROUTINE HexInit