#include <defines.h>
#ifdef __GAMMA_TRANSPORT
!--- CNJ Edit : Driver Routine for Gamma MOC Calculation
SUBROUTINE GammaMOCSweep(Core, RayInfo, FmInfo, PE, nTracerCntl, ItrCntl)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,      RayInfo_Type,       FmInfo_Type,        &
                           PE_TYPE,            FxrInfo_Type
USE CNTL,           ONLY : nTracerCntl_Type
USE itrcntl_mod,    ONLY : ItrCntl_TYPE
USE CORE_MOD,       ONLY : GroupInfo
USE GammaCore_mod,  ONLY : gphis,              gphim,              gPhiAngIn,          &
                           gJout,              gPower,             GamGroupInfo,       &
                           LocalNPower
USE GamMOC_MOD,     ONLY : SetGamSrc,          SetGamSrcNM,        SetGamP1Src,        &
                           SetGamMacXs,        SetGamMacXsNM,                          &
                           RayTraceGamma,      RayTraceGamma_Pn,                       &
                           GammaMocResidual,   GammaPowerUpdate,                       &
                           gphis1g,            gphim1g,            gxst1g,             &
                           gsrc1g,             gsrcm1g,            gPhiAngIn1g,        &
                           gJout1g,                                                    &
                           gphisnm,            gxstnm,             gsrcnm,             &
                           gPhiAngInnm,        gJoutnm,            NeutronLocalQUpdate
USE IOUTIL,         ONLY : message
USE Timer,          ONLY : nTracer_dclock,     TimeChk
USE FILES,          ONLY : io8
#ifdef MPI_ENV
USE MPICOMM_MOD,    ONLY : BCAST,              MPI_SYNC
#endif
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(PE_TYPE) :: PE
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(ItrCntl_TYPE) :: ItrCntl

REAL, POINTER :: phis(:, :, :), phisnm(:, :)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
REAL, POINTER :: AxSrc(:, :, :), AxSrc1g(:)
INTEGER :: ig, iz, i, l, jsweep, InIter, nInIter
INTEGER :: ng, ngg, nGroupInfo, nFsr
INTEGER :: myzb, myze
INTEGER :: GrpBeg, GrpEnd
INTEGER, SAVE :: iter = 0
REAL :: resconv, reserr
REAL :: TimeRt1gBeg, TimeRt1gEnd, TimeRtngBeg, TimeRtngEnd, TimeRtElapsed
REAL :: MocTimeBeg, MocTimeEnd
LOGICAL :: lJout, lxslib, l3dim, lscat1, lTrCorrection
LOGICAL :: MASTER, RTMASTER

MOCTimeBeg = nTracer_dclock(FALSE, FALSE)

Fxr => FmInfo%Fxr
phis => FmInfo%phis
Master = PE%Master
RTMaster = PE%RTMaster
ng = GroupInfo%ng
ngg = GamGroupInfo%ngg
nFsr = Core%nCoreFsr
myzb = PE%myzb; myze = PE%myze
nGroupInfo = 2; nInIter = 1
l3dim = nTracerCntl%l3dim
lscat1 = nTracerCntl%lGammaScat1
lTrCorrection = .NOT. lscat1
resconv = itrcntl%resconv

IF (.NOT. GamGroupInfo%lUpScat) nGroupInfo = 1

IF (.NOT. lscat1) THEN
  ALLOCATE(phisnm(ng, nFsr))
ENDIF

iter = iter + 1

IF (lscat1) THEN

  TimeRtElapsed = 0.
  DO jsweep = 1, nGroupInfo
    GrpBeg = 1; GrpEnd = ngg
    IF(jsweep .GT. 1) THEN
      GrpBeg = GamGroupInfo%UpScatRange(1); GrpEnd = GamGroupInfo%UpScatRange(2)
    ENDIF
    DO ig = GrpBeg, GrpEnd
      TimeRt1gBeg = nTracer_dclock(FALSE, FALSE)
      DO iz = myzb, myze
        ljout = FALSE
        IF (.NOT. Core%lFuelPlane(iz)) CYCLE
        gphis1g => gphis(:, ig, iz)
        gphim1g => gphim(:, :, ig, iz)
        gPhiAngIn1g => gPhiAngIn(:, :, ig, iz)
        gJout1g => gJout(:, :, :, ig, iz)
        DO InIter = 1, nInIter
          IF (InIter .EQ. nInIter) ljout = TRUE
          IF (RTMASTER) THEN
            CALL SetGamMacXs(Core, Fxr(:, iz), gxst1g, iz, ig, ngg, lTrCorrection, PE)
            CALL SetGamSrc(Core, Fxr(:, iz), gsrc1g, phis, gphis, AxSrc1g, gxst1g, iz, ig, ng, ngg,                   &
                           l3dim, lscat1, PE, GamGroupInfo)
            CALL SetGamP1Src(Core, Fxr(:, iz), gsrcm1g, gphim, gxst1g, iz, ig, ngg, lscat1, nTracerCntl%GammaScatOd, PE)
          ENDIF
          CALL RayTraceGamma_Pn(RayInfo, Core, gphis1g, gphim1g, gPhiAngIn1g, gxst1g, gsrc1g, gsrcm1g, gJout1g,       &
                                iz, nTracerCntl%GammaScatOd, lJout)
        ENDDO
      ENDDO
      TimeRt1gEnd = nTracer_dclock(FALSE, FALSE)
#ifdef MPI_ENV
        CALL MPI_SYNC(PE%MPI_RTMASTER_COMM)
#endif
      TimeRt1gEnd = nTracer_dclock(FALSE, FALSE)
      TimeRtElapsed = TimeRtElapsed + TimeRt1gEnd - TimeRt1gBeg
    ENDDO
  ENDDO

ELSE

  TimeRtElapsed = 0.
  DO iz = myzb, myze
    IF (.NOT. Core%lFuelPlane(iz)) CYCLE
    IF (RTMASTER) THEN
      DO ig = 1, ng
        phisnm(ig, :) = phis(:, iz, ig)
      ENDDO
      gphisnm => gphis(:, :, iz)
      gPhiAngInnm => gPhiAngIn(:, :, :, iz)
      gJoutnm => gJout(:, :, :, :, iz)
      CALL SetGamMacXsNM(Core, Fxr(:, iz), gxstnm, iz, ngg, lTrCorrection, PE)
    ENDIF
    DO jsweep = 1, nGroupInfo
      GrpBeg = 1; GrpEnd = ngg
      IF(jsweep .GT. 1) THEN
        GrpBeg = GamGroupInfo%UpScatRange(1); GrpEnd = GamGroupInfo%UpScatRange(2)
      ENDIF
      TimeRtngBeg = nTracer_dclock(FALSE, FALSE)
      DO InIter = 1, nInIter
        ljout = FALSE
        IF (InIter .EQ. nInIter) ljout = TRUE
        IF (RTMASTER) THEN
          CALL SetGamSrcNM(Core, Fxr(:, iz), gsrcnm, phisnm, gphisnm, AxSrc, gxstnm, iz,                              &
                           GrpBeg, GrpEnd, ng, ngg, l3dim, lscat1, PE)
        ENDIF
        CALL RayTraceGamma(RayInfo, Core, gphisnm, gPhiAngInnm, gxstnm, gsrcnm, gJoutnm,                              &
                           iz, GrpBeg, GrpEnd, ljout)
      ENDDO
      TimeRtngEnd = nTracer_dclock(FALSE, FALSE)
      TimeRtElapsed = TimeRtElapsed + (TimeRtngEnd - TimeRtngBeg)
    ENDDO
#ifdef MPI_ENV
    CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif
  ENDDO

ENDIF

IF (RTMASTER) reserr = GammaMocResidual(Core, FmInfo, GroupInfo, PE, nTracerCntl)

WRITE(mesg, '(2x, a10, i2, a14, 1p, e10.4)') 'Iteration ', iter, ' / Residual = ', reserr
IF(Master) CALL message(io8, TRUE, TRUE, mesg)

ItrCntl%lGammaConv = reserr .LT. resconv

#ifdef MPI_ENV
CALL BCAST(reserr, PE%MPI_RT_COMM, 0)
#endif

IF (RTMASTER) CALL GammaPowerUpdate(Core, Fxr, gphis, gPower, myzb, myze, ngg, PE)
IF (RTMASTER) CALL NeutronLocalQUpdate(Core, Fxr, phis, LocalNPower, myzb, myze, ng, PE) !-- JSU EDIT 2017.09.15.

MocTimeEnd = nTracer_dclock(FALSE, FALSE)
TimeRtElapsed = MocTimeEnd - MocTimeBeg
TimeChk%MocTime = TimeChk%MocTime + TimeRtElapsed

IF (.NOT. lscat1) THEN
  DEALLOCATE(phisnm)
ENDIF

END SUBROUTINE
#endif
