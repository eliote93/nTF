#include <defines.h>
!--- CNJ Edit : Node Majors
MODULE SubGrpFspMLG_NM

USE Util_SGFSP_MLG_NM

IMPLICIT NONE

CONTAINS
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SubGrpFsp_MLG_NM(Core, Fxr, THInfo, RayInfo, GroupInfo)

USE allocs
USE OMP_LIB
USE PARAM,     ONLY : TRUE, FALSE, ZERO, ONE, EPSM3, VoidCell, mesg
USE TYPEDEF,   ONLY : CoreInfo_Type, RayInfo_Type, Fxrinfo_type, GroupInfo_Type, THInfo_Type
USE CNTL,      ONLY : nTracerCntl
USE PE_MOD,    ONLY : PE
USE FILES,     ONLY : io8
USE IOUTIL,    ONLY : message
USE TIMER,     ONLY : nTracer_dclock, TimeChk
USE xslib_mod, ONLY : mlgdata, mlgdata0
USE MOC_MOD,   ONLY : TrackingDat, RayTrace_NM
USE GEOM,      ONLY : ng, nbd

#ifdef MPI_ENV
USE MPICOMM_MOD, ONLY : REDUCE, MPI_SYNC
#endif

IMPLICIT NONE

TYPE (CoreInfo_Type)  :: Core
TYPE (THInfo_Type)    :: THInfo
TYPE (RayInfo_Type)   :: RayInfo
TYPE (GroupInfo_Type) :: GroupInfo

TYPE (FxrInfo_Type), POINTER, DIMENSION(:,:) :: Fxr
! ----------------------------------------------------
INTEGER :: ig, iz, ilv, iter, itersum, itermax, ithr, mg, nlv, nFsr, nFxr, nPhiAngSv, gb, ge, myzb, myze, nPol, nAzi, nThr, nModRay, nxy
REAL :: errmax, Tbeg, Tend, rtTbeg, rtTend
LOGICAL :: lCLD, lAIC, ldcmp

REAL, POINTER, DIMENSION(:,:)     :: phisNg, phisdNg, SiglpNg, xstNg, srcNg
REAL, POINTER, DIMENSION(:,:,:)   :: PhiAngInNg
REAL, POINTER, DIMENSION(:,:,:,:) :: JoutNg
! ----------------------------------------------------

Tbeg = nTracer_dclock(FALSE, FALSE)

nFxr = Core%nCoreFxr
nFsr = Core%nCoreFsr
nxy  = Core%nxy

myzb = PE%myzb
myze = PE%myze

nAzi      = RayInfo%nAziAngle
nPol      = RayInfo%nPolarAngle
nPhiAngSv = RayInfo%nPhiAngSv
nModRay   = RayInfo%nModRay

ldcmp = nTracerCntl%lDomainDcmp

itermax = 100
IF (any(Core%RadBC(1:4) .EQ. VoidCell)) itermax = 1

nThr = PE%nThread
CALL omp_set_num_threads(nThr)

IF (ldcmp) THEN
  DO ithr = 1, nThr
    CALL dmalloc(TrackingDat(ithr)%phisNg,    ng, nFsr)
    CALL dmalloc(TrackingDat(ithr)%JoutNg, 3, ng, nbd, nxy)
    
    TrackingDat(ithr)%lAllocNM = TRUE
  END DO
END IF

WRITE (mesg, '(A, F10.2, A)') "Reference Fuel Temperature", THInfo%RefFuelTemp(0), "C"
IF (PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)

WRITE (mesg, '(A)') 'Solving Subgroup FSP...'
IF (PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)
! ----------------------------------------------------
itersum = 0
errmax  = ZERO

DO iz = myzb, myze
  IF (.NOT. Core%lFuelPlane(iz)) CYCLE

  gb = GroupInfo%nofg + 1
  ge = GroupInfo%nofg + GroupInfo%norg
  mg = mlgdata(iz)%f_nmaclv * GroupInfo%norg
  
  CALL dmalloc(phisNg,  mg, nFsr)
  CALL dmalloc(phisdNg, mg, nFsr)
  CALL dmalloc(SiglpNg, mg, nFxr)
  CALL dmalloc(xstNg,   mg, nFsr)
  CALL dmalloc(srcNg,   mg, nFsr)
  CALL dmalloc(PhiAngInNg, nPol, nPhiAngSv, mg)
  
  phisNg     = ONE
  PhiAngInNg = ONE
  PhiAngInNg(:, 1, :) = ZERO
  
  DO ithr = 1, nThr
    DEALLOCATE (TrackingDat(ithr)%phisNg)
    CALL dmalloc(TrackingDat(ithr)%phisNg, mg, nFsr)
  END DO
    
  CALL UpdtFnAdj_NM      (Core, Fxr, iz, gb, ge)
  CALL SetPlnLsigP_MLG_NM(Core, Fxr, SiglpNg, xstNg, iz, gb, ge)
  CALL SetSubGrpSrc_NM   (Core, Fxr, SiglpNg, xstNg, srcNg, iz, 1, mg)
  
  DO iter = 1, itermax
    CALL CopyFlux(phisNg, phisdNg, nFsr, mg)
    
    rtTbeg = nTracer_dclock(FALSE, FALSE)
    CALL RayTrace_NM(RayInfo, Core, phisNg, PhiAngInNg, xstNg, srcNg, JoutNg, iz, 1, mg, FALSE)
    rtTend = nTracer_dclock(FALSE, FALSE)
    TimeChk%NetRTSubGrpTime = TimeChk%NetRTSubGrpTime + (rtTend - rtTbeg)
    
    CALL EquipXSGen_MLG_NM(Core, Fxr, SiglpNg, xstNg, phisNg, iz, gb, ge)
    CALL UpdtFtAdj_NM     (Core, Fxr, iz, gb, ge)
    
    errmax = SubGrpFspErr_NM(phisNg, phisdNg, nFsr, mg)
    
    IF (errmax .LT. EPSM3) EXIT
    IF (iter .EQ. itermax) EXIT
    
    CALL SetPlnLsigP_MLG_NM(Core, Fxr, SiglpNg, xstNg, iz, gb, ge)
    CALL SetSubGrpSrc_NM   (Core, Fxr, SiglpNg, xstNg, srcNg, iz, 1, mg)
  END DO
  
  WRITE (mesg, '(2X, A, I3, A, I5,1P, E13.3)') '[Fuel] Pln', iz, '  In_itr', iter, errmax
  IF (PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)
  
  itersum = itersum + iter

  DEALLOCATE(phisNg, phisdNg, PhiAngInNg, SiglpNg, xstNg, srcNg)
END DO

#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
CALL REDUCE(itersum, iter, PE%MPI_NTRACER_COMM, FALSE)
#endif
! ----------------------------------------------------
itersum = 0
errmax  = ZERO

DO iz = myzb, myze
  IF (.NOT. Core%lFuelPlane(iz)) CYCLE
  IF (.NOT. Core%lCladPlane(iz)) CYCLE

  lCLD = TRUE
  lAIC = FALSE
  mg   = mlgdata0%c_nmaclv1G

  CALL dmalloc(phisNg,   mg, nFsr)
  CALL dmalloc(phisdNg,  mg, nFsr)
  CALL dmalloc(SiglpNg,  mg, nFxr)
  CALL dmalloc(xstNg,    mg, nFsr)
  CALL dmalloc(srcNg,    mg, nFsr)
  CALL dmalloc(PhiAngInNg, nPol, nPhiAngSv, mg)
  
  phisNg     = ONE
  PhiAngInNg = ONE
  PhiAngInNg(:, 1, :) = ZERO
  
  DO ithr = 1, nThr
    DEALLOCATE (TrackingDat(ithr)%phisNg)
    CALL dmalloc(TrackingDat(ithr)%phisNg, mg, nFsr)
  END DO
    
  CALL SetPlnLsigP_1gMLG_NM(Core, Fxr, SiglpNg, xstNg, iz, lCLD, lAIC)
  CALL SetSubGrpSrc_NM     (Core, Fxr, SiglpNg, xstNg, srcNg, iz, 1, mg)
  
  DO iter = 1, itermax
    CALL CopyFlux(phisNg, phisdNg, nFsr, mg)
    
    rtTbeg = nTracer_dclock(FALSE, FALSE)
    CALL RayTrace_NM(RayInfo, Core, phisNg, PhiAngInNg, xstNg, srcNg, JoutNg, iz, 1, mg, FALSE)
    rtTend = nTracer_dclock(FALSE, FALSE)
    TimeChk%NetRTSubGrpTime = TimeChk%NetRTSubGrpTime + (rtTend - rtTbeg)
    
    errmax = SubGrpFspErr_NM(phisNg, phisdNg, nFsr, mg)
    
    IF (errmax .LT. EPSM3) EXIT
    IF (iter .EQ. itermax) EXIT
  END DO
  
  WRITE (mesg, '(2X, A, I3, A, I5,1P, E13.3)') '[Clad] Pln', iz, '  In_itr', iter, errmax
  IF (PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)
  
  itersum = itersum + iter
  
  CALL EquipXSGen_1gMLG_NM(Core, Fxr, SiglpNg, phisNg, xstNg, iz, mg, lCLD, lAIC)

  DEALLOCATE(phisNg, phisdNg, PhiAngInNg, SiglpNg, xstNg, srcNg)
END DO

#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
CALL REDUCE(itersum, iter, PE%MPI_NTRACER_COMM, FALSE)
#endif
! ----------------------------------------------------
itersum = 0
errmax  = ZERO

DO iz = myzb, myze
  IF (.NOT. Core%lFuelPlane(iz)) CYCLE
  IF (.NOT. Core%lAICPlane(iz)) CYCLE

  lCLD = FALSE
  lAIC = TRUE
  mg   = mlgdata(iz)%f_nmaclv1G

  CALL dmalloc(phisNg,   mg, nFsr)
  CALL dmalloc(phisdNg,  mg, nFsr)
  CALL dmalloc(SiglpNg,  mg, nFxr)
  CALL dmalloc(xstNg,    mg, nFsr)
  CALL dmalloc(srcNg,    mg, nFsr)
  CALL dmalloc(PhiAngInNg, nPol, nPhiAngSv, mg)
  
  phisNg     = ONE
  PhiAngInNg = ONE
  PhiAngInNg(:, :, 1) = ZERO
  
  DO ithr = 1, nThr
    DEALLOCATE (TrackingDat(ithr)%phisNg)
    CALL dmalloc(TrackingDat(ithr)%phisNg, mg, nFsr)
  END DO
    
  CALL SetPlnLsigP_1gMLG_NM(Core, Fxr, SiglpNg, xstNg, iz, lCLD, lAIC)
  CALL SetSubGrpSrc_NM     (Core, Fxr, SiglpNg, xstNg, srcNg, iz, 1, mg)
  
  DO iter = 1, itermax
    CALL CopyFlux(phisNg, phisdNg, nFsr, mg)
    
    rtTbeg = nTracer_dclock(FALSE, FALSE)
    CALL RayTrace_NM(RayInfo, Core, phisNg, PhiAngInNg, xstNg, srcNg, JoutNg, iz, 1, mg, FALSE)
    rtTend = nTracer_dclock(FALSE, FALSE)
    TimeChk%NetRTSubGrpTime = TimeChk%NetRTSubGrpTime + (rtTend - rtTbeg)
    
    errmax = SubGrpFspErr_NM(phisNg, phisdNg, nFsr, mg)
    
    IF (errmax .LT. EPSM3) EXIT
    IF (iter .EQ. itermax) EXIT
  END DO
  
  WRITE (mesg, '(2X, A, I3, A, I5,1P, E13.3)') '[AIC] Pln', iz, '  In_itr', iter, errmax
  IF (PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)
  
  itersum = itersum + iter
  
  CALL EquipXSGen_1gMLG_NM(Core, Fxr, SiglpNg, phisNg, xstNg, iz, mg, lCLD, lAIC)

  DEALLOCATE (phisNg, phisdNg, PhiAngInNg, SiglpNg, xstNg, srcNg)
END DO
! ----------------------------------------------------
#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
CALL REDUCE(itersum, iter, PE%MPI_NTRACER_COMM, FALSE)
#endif

DO ithr = 1, nThr
  DEALLOCATE (TrackingDat(ithr)%phisNg)
  CALL dmalloc(TrackingDat(ithr)%phisNg, ng, nFsr)
END DO

IF (ldcmp) THEN
  DO ithr = 1, nThr
    DEALLOCATE (TrackingDat(ithr)%phisNg)
    DEALLOCATE (TrackingDat(ithr)%JoutNg)
  END DO
END IF

nTracerCntl%lSubGrpSweep = TRUE

#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif

Tend = nTracer_dclock(FALSE, FALSE)
TimeChk%SubGrpTime = TimeChk%SubGrpTime + (Tend - Tbeg)
! ----------------------------------------------------

END SUBROUTINE SubGrpFsp_MLG_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE CopyFlux(phis, phisd, nFsr, ng)

IMPLICIT NONE

REAL, POINTER :: phis(:, :), phisd(:, :)
INTEGER :: nFsr, ng
INTEGER :: ifsr, ig

#ifdef __INTEL_MKL
CALL dcopy(ng * nFsr, phis, 1, phisd, 1)
#else
!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ifsr = 1, nFsr
  DO ig = 1, ng
    phisd(ig, ifsr) = phis(ig, ifsr)
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL
#endif

END SUBROUTINE
! ------------------------------------------------------------------------------------------------------------
FUNCTION SubGrpFspErr_NM(phis, phisd, nFsr, ng) RESULT(errmax)

IMPLICIT NONE

REAL, POINTER :: phis(:, :), phisd(:, :)
INTEGER :: nFsr, ng
INTEGER :: ifsr, ig
REAL :: errmax, err

errmax = 0.0
!$OMP PARALLEL DO REDUCTION(MAX : errmax) PRIVATE(err)
DO ifsr = 1, nFsr
  DO ig = 1, ng
    IF (phis(ig, ifsr) .LT. 0.0) CYCLE
    err = abs((phis(ig, ifsr) - phisd(ig, ifsr)) / phisd(ig, ifsr))
    errmax = max(err, errmax)
  ENDDO
ENDDO
!$OMP END PARALLEL DO

END FUNCTION
! ------------------------------------------------------------------------------------------------------------

END MODULE SubGrpFspMLG_NM