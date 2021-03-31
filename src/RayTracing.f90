#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceGM_One(RayInfo, CoreInfo, phis, PhiAngIn, xst, src, jout, iz, ljout, FastMocLv, lAFSS)

USE PARAM
USE TYPEDEF,  ONLY :  RayInfo_Type,      coreinfo_type,                                       &
                      Pin_Type,          Asy_Type,        AsyInfo_Type,     PinInfo_Type,     &
                      Cell_Type,                                                              &
                      AziAngleInfo_Type, PolarAngle_Type, ModRayInfo_type,  AsyRayInfo_type,  &
                      CoreRayInfo_Type,  RotRayInfo_Type, CellRayInfo_type
USE Moc_Mod, ONLY :   nMaxRaySeg,        nMaxCellRay,     nMaxAsyRay,       nMaxCoreRay,      &
                      Expa,              Expb,                                                &
                      ApproxExp,         RayTraceGM,    RayTraceGM_AFSS
USE PE_MOD,  ONLY :   PE

#ifdef MPI_ENV
USE MPICOMM_MOD, ONLY : REDUCE,          BCAST
#endif
USE BasicOperation, ONLY : CP_CA, CP_VA
USE ALLOCS

IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phis(:), PhiAngIn(:, :), xst(:), src(:), jout(:, :, :)

INTEGER :: iz
LOGICAL :: ljout
INTEGER, OPTIONAL :: FastMocLv
LOGICAL, OPTIONAL :: lAFSS

INTEGER :: FastMOCLv0
LOGICAL :: lAFSS0

FastMocLv0 = 0
lAFSS0 = .FALSE.
IF(Present(FastMocLv)) FastMocLv0 = FastMocLv
IF(Present(lAFSS)) lAFSS0 = lAFSS

IF (.NOT. lAFSS0) THEN
  CALL RayTraceGM     (RayInfo, CoreInfo, phis, PhiAngIn, xst, src, jout, iz, ljout, FastMocLv0, lAFSS0)
ELSE
  CALL RayTraceGM_AFSS(RayInfo, CoreInfo, phis, PhiAngIn, xst, src, jout, iz, ljout, FastMocLv0, lAFSS0)
ENDIF

END SUBROUTINE RayTraceGM_One
! ------------------------------------------------------------------------------------------------------------