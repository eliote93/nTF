SUBROUTINE RayTraceP1(RayInfo, CoreInfo, phis, phim, PhiAngIn, xst, src, srcm, jout, iz, ljout, ScatOd, FastMocLv, lAFSS)
USE PARAM
USE TYPEDEF,  ONLY :  RayInfo_Type,      coreinfo_type,                                       &
                      Pin_Type,          Asy_Type,        AsyInfo_Type,     PinInfo_Type,     &
                      Cell_Type,                                                              &
                      AziAngleInfo_Type, PolarAngle_Type, ModRayInfo_type,  AsyRayInfo_type,  &
                      CoreRayInfo_Type,  RotRayInfo_Type, CellRayInfo_type
USE Moc_Mod, ONLY :   nMaxRaySeg,        nMaxCellRay,     nMaxAsyRay,       nMaxCoreRay,      &
                      Expa,              Expb,                                                &
                      ApproxExp,                                                              &
                      RayTraceP1_OMP,    RayTraceP1_OMP_AFSS
USE BasicOperation, ONLY : CP_CA, CP_VA                    
USE ALLOCS
IMPLICIT NONE
!Input Arguments
TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phis(:), PhiAngIn(:, :), xst(:), src(:), jout(:, :, :), srcm(:, :)
REAL, POINTER :: phim(:, :)
INTEGER :: iz
LOGICAL :: ljout
INTEGER :: SCatOd
INTEGER, OPTIONAL :: FastMocLv
LOGICAL, OPTIONAL :: lAFSS

INTEGER :: FastMocLv0
LOGICAL :: lAFSS0

FastMocLv0 = 0
lAFSS0 = .FALSE.
IF(Present(FastMocLv)) FastMocLv0 = FastMocLv
IF(Present(lAFSS)) lAFSS0 = lAFSS

IF (.NOT. lAFSS0) THEN
  CALL RayTraceP1_OMP(RayInfo, CoreInfo, phis, phim, PhiAngIn, xst, src, srcm, jout, iz, ljout, ScatOd, FastMocLv0, lAFSS0)
ELSE
  CALL RayTraceP1_OMP_AFSS(RayInfo, CoreInfo, phis, phim, PhiAngIn, xst, src, srcm, jout, iz, ljout, ScatOd, FastMocLv0, lAFSS0)
ENDIF

END SUBROUTINE