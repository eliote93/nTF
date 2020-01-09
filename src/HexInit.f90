#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
!                                     01. HEX INIT
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexInit

USE ALLOCS
USE Geom,          ONLY : Core
USE RAYS,          ONLY : RayInfo, RayInfo4Cmfd
USE CNTL,          ONLY : nTracerCntl
USE PE_MOD,        ONLY : PE
USE MPIConfig_Mod, ONLY : SetGeomPEVariables, SetRayPEVariables, SetDcmpPEVariables
USE HexCmfd,       ONLY : HexRayInfo4CmfdGen
USE HexData,       ONLY : hRotRay, hcRay

IMPLICIT NONE

INTEGER :: iAng, nAziAng, nRotRay, irRay, icRay, jcRay, iAzm, nRay
! ----------------------------------------------------

! ----------------------------------------------------
!               01. INIT : Basic
! ----------------------------------------------------
CALL SetGeomPEVariables(PE)

IF(nTracerCntl%lXsLib .AND. nTracerCntl%lrestrmt .AND. PE%RTmaster) CALL AllocResIsoInfo()

CALL PrepFxr(nTracerCntl%lXsLib, nTracerCntl%lfxrlib)

CALL SetDcmpPEVariables(PE)

IF(nTracerCntl%lDcpl) CALL initDcpl()

#ifdef MPI_ENV
CALL SetRayPEVariables(PE, Core, RayInfo)
#endif

CALL HexRayInfo4CmfdGen(RayInfo, RayInfo4Cmfd)
! ----------------------------------------------------
!               02. INIT : NM
! ----------------------------------------------------
nAziAng = RayInfo%nAziAngle
nRotRay = RayInfo%nRotRay

CALL Dmalloc0(RayInfo%RotRayAziList, 0, nRotRay, 1, nAziAng / 2)

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