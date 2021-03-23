#include <defines.h>
#include <CUDADEFINES.h>

#ifdef __PGI

SUBROUTINE CUDATranMOCSweep(Core, RayInfo, FmInfo, TranInfo, TranCntl)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,       RayInfo_Type,       FmInfo_Type,        TranInfo_Type,        &
                           TranCntl_Type,       FxrInfo_Type
USE PE_Mod,         ONLY : PE
USE CNTL,           ONLY : nTracerCntl
USE itrcntl_mod,    ONLY : ItrCntl
USE Core_mod,       ONLY : GroupInfo
USE TRANMOC_MOD,    ONLY : SetPrecParam,        PrecSrcUpdt,        SetTranSrcNM,       TranMocResidualError, &
                           SetTranMocEnv,       PrecSrcKUpdt
USE MOC_MOD,        ONLY : SetRtMacXsNM_Cusping,SetRtSrcNM_Cusping, SetRtP1SrcNM,       PsiUpdate,            &
                           PsiErr,              PseudoAbsorptionNM, PowerUpdate,        GetNeighborMocFlux
USE CUDA_MASTER
USE CUDA_MOC
USE CUDA_INIT,      ONLY : AllocMOCVar,         DeallocMOCVar,      CopyRayVar,         DeleteRayVar
USE SUBGRP_MOD,     ONLY : FxrChiGen
USE ioutil,         ONLY : message
USE timer,          ONLY : nTracer_dclock,      TimeChk
USE files,          ONLY : io8
#ifdef MPI_ENV
USE MPIComm_Mod,    ONLY : BCAST,               MPI_SYNC,           MPI_MAX_REAL
#endif
USE MOC_COMMON,     ONLY : SetMOCPsi,           SetMOCtrXS,         SetMOCSource,       SetMOCPower,          &
                           CopyXS,              SetCoreMacXS_Cusping, SetMOCPsi_iter
USE SPH_mod,        ONLY : ssphfnm,             calcPinSSPH
USE XSLIB_MOD,      ONLY : igresb,              igrese
USE OMP_LIB
USE CUDAFOR
USE OPENACC
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
REAL, POINTER :: psi(:, :), psid(:, :), AxSrc(:, :, :), AxPXS(:, :, :)
REAL, POINTER :: Prec(:, :, :), ResSrc(:, :, :), ResSrcD(:, :, :)
REAL, POINTER :: TranPsi(:, :), TranPsid(:, :)
REAL, POINTER :: PrecSrc(:, :)
REAL :: eigv0, inv_eigv0
REAL :: eigerr, fiserr, peigv, reserr, errdat(3)
REAL :: eigconv, psiconv
REAL :: rtTime, TimeRtBeg, TimeRtEnd, TimeRtElapsed, MocTimeBeg, MocTimeEnd
INTEGER :: ig, igb, ifsr, ipin, isurf, iz, i, ierr, InIter, nInIter, iray
INTEGER :: myzb, myze
INTEGER :: nFsr, nxy, nPhiAngSv, nPolarAngle, nGroupInfo, nMoment
INTEGER :: GrpBeg, GrpEnd, gb, ge, ng, ngBlock, ngLocal
LOGICAL :: lJout, lxslib, l3dim, lscat1, lTrCorrection, lRST, lsSPH, lsSPHreg
LOGICAL :: lPsiUpdt

!--- CUDA Host Variables --------------------------------------------------------------------------
TYPE hostMOC_Type
  REAL(GPU_PRECISION), ALLOCATABLE, PINNED :: PhiAngIn(:, :, :)
  REAL(GPU_FLUX_PRECISION), ALLOCATABLE, PINNED :: phis(:, :), phim(:, :, :), jout(:, :, :, :)
  REAL(GPU_SOURCE_PRECISION), ALLOCATABLE, PINNED :: xst(:, :), src(:, :), srcm(:, :, :)
END TYPE

TYPE(hostMOC_Type), POINTER :: hostMOC(:)
REAL, POINTER :: phis(:, :), phim(:, :, :), PhiAngIn(:, :, :)
REAL, POINTER :: xst(:, :), src(:, :), srcm(:, :, :)
REAL, POINTER :: Jout(:, :, :, :)
REAL, POINTER :: trsrc(:,:), TranPhi(:, :)
!--------------------------------------------------------------------------------------------------
LOGICAL, SAVE :: lFirst
DATA lfirst /TRUE/

MOCTimeBeg = nTracer_dclock(FALSE, FALSE)

Fxr => FmInfo%Fxr
psi => FmInfo%psi; psid => FmInfo%psid
AxSrc => FmInfo%AxSrc; AxPXS => FmInfo%AxPXS

Prec => FmInfo%Prec; ResSrc => FmInfo%ResSrc; ResSrcD => FmInfo%ResSrcD
PrecSrc => FmInfo%PrecSrc
TranPsi => FmInfo%TranPsi; TranPsid => FmInfo%TranPsid

eigv0 = TranInfo%eigv0; inv_eigv0 = 1./eigv0

myzb = PE%myzb; myze = PE%myze
nFsr = cuGeometry%nFsr; nxy = cuGeometry%nxy
nPolarAngle = RayInfo%nPolarAngle; nPhiAngSv = RayInfo%nPhiAngSv; ng = cuGeometry%ng
!nInIter = 2; nGroupInfo = 2; IF (.NOT. GroupInfo%lUpScat) nGroupInfo = 1
nInIter = 3; nGroupInfo = 1

IF (nTracerCntl%ScatOd .EQ. 1) nMoment = 2
IF (nTracerCntl%ScatOd .EQ. 2) nMoment = 5
IF (nTracerCntl%ScatOd .EQ. 3) nMoment = 9

lxslib = nTracerCntl%lXslib; l3dim = nTracerCntl%l3dim; lscat1 = nTracerCntl%lscat1
lRST = nTracerCntl%lRST; lsSPH = nTracerCntl%lsSPH; lsSPHreg = nTracerCntl%lsSPHreg
lTrCorrection = nTracerCntl%lTrCorrection

psiconv = itrcntl%psiconv; eigconv = itrcntl%eigconv

IF (lsSPHreg) THEN
  WRITE(mesg, '(a24)') 'Calculating Pin SSPH...'
  IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)
  CALL calcPinSSPH(Core, Fxr, PE)
END IF

IF(lFirst) THEN
  CALL SetTranMOCEnv(Core, TranInfo, PE)
  lFirst = .FALSE.
ENDIF

IF (nTracerCntl%lMacro) THEN
  CALL SetMOCPsi(Core, FmInfo%phis, psi)
ELSE
  CALL PsiUpdate(Core, Fxr, FmInfo%phis, psi, myzb, myze, ng, lxslib, GroupInfo)
ENDIF
psi = inv_eigv0 * psi
psid(:, myzb : myze) = psi(:, myzb : myze)
lPsiUpdt = .FALSE.

CALL SetPrecParam(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, PE)
IF(TranCntl%lchidk) THEN
  CALL PrecSrcKUpdt(Core, Fxr, FmInfo%PrecSrcK, Prec, TranPsi, TranPsid, GroupInfo, TranCntl, nTRACERCntl, PE)
ELSE
  CALL PrecSrcUpdt(Core, Fxr, PrecSrc, Prec, TranPsi, TranPsid, GroupInfo, TranCntl, nTRACERCntl, PE)
END IF

itrcntl%mocit = itrcntl%mocit + 1
WRITE(mesg, '(a36,I5,a3)') 'Performing Ray Tracing for Transient', itrcntl%mocit, '...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)
CALL FxrChiGen(Core, Fxr, FmInfo, GroupInfo, PE, nTracerCntl)
IF(nTracerCntl%lCusping_MPI) THEN
  CALL GetNeighborMocFlux(FmInfo%phis, FmInfo%neighphis, nFsr, myzb, myze, 1, ng, Core%nz, Core%AxBC)
  IF(nTracerCntl%lMacro) CALL SetCoreMacXS_Cusping(Core, FmInfo)
END  IF

CALL omp_set_nested(.TRUE.)
CALL omp_set_max_active_levels(2)

TimeRtBeg = nTracer_dclock(FALSE, FALSE)

CALL acc_set_device_num(MPI_CUDA_SHARED_RANK, acc_device_nvidia)
CALL omp_set_num_threads(cuCntl%nMOCHelper)

CALL AllocMOCVar(Core, RayInfo, cuMOC, cuDevice, ng, nTracerCntl%lScat1, cuCntl%lAsync, .NOT. cuCntl%lAsync)

IF (.NOT. nTracerCntl%lScat1) THEN

  ngBlock = (ng - 1) / P0_BLOCK_SIZE + 1
  ALLOCATE(hostMOC(ngBlock))
  DO igb = 1, ngBlock
    ALLOCATE(hostMOC(igb)%phis(P0_BLOCK_SIZE, nFsr))
    ALLOCATE(hostMOC(igb)%PhiAngIn(nPolarAngle, P0_BLOCK_SIZE, nPhiAngSv))
    ALLOCATE(hostMOC(igb)%jout(3, P0_BLOCK_SIZE, 4, nxy))
    ALLOCATE(hostMOC(igb)%xst(P0_BLOCK_SIZE, nFsr))
    ALLOCATE(hostMOC(igb)%src(P0_BLOCK_SIZE, nFsr))
  ENDDO

  ALLOCATE(phis(ng, nFsr))
  ALLOCATE(TranPhi(ng, nFsr))
  ALLOCATE(Jout(3, ng, 4, nxy))
  ALLOCATE(xst(ng, nFsr))
  ALLOCATE(src(P0_BLOCK_SIZE, nFsr))
  ALLOCATE(trsrc(P0_BLOCK_SIZE, nFsr))

  rtTime = 0.0

  DO iz = cuDevice%myzb, cuDevice%myze
    IF (.NOT. cuDevice%lRayStatic) CALL CopyRayVar(cuGeometry%AxType(iz))
    !$OMP PARALLEL DO COLLAPSE(2)
    DO ig = 1, ng
      DO ifsr = 1, nFsr
        phis(ig, ifsr) = FmInfo%phis(ifsr, iz, ig)
        TranPhi(ig, ifsr) = FmInfo%TranPhi(ifsr, iz, ig)
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    IF (nTracerCntl%lMacro) THEN
      CALL SetMOCtrXS(Core, Fxr, xst, iz, lsSPH, lsSPHreg, lTrCorrection)
    ELSE
      CALL SetRtMacXsNM_Cusping(Core, FmInfo, Fxr(:, iz), xst, iz, ng, lxslib, lTrCorrection, lRST, lsSPH, lsSPHreg, PE)
    ENDIF
    CALL PseudoAbsorptionNM(Core, Fxr(:, iz), AxPXS, xst, iz, ng, GroupInfo, l3dim)
    DO i = 1, nGroupInfo
      GrpBeg = 1; GrpEnd = ng
      IF (i .GT. 1) THEN
        GrpBeg = GroupInfo%UpScatRange(1); GrpEnd = GroupInfo%UpScatRange(2)
      ENDIF
      ngBlock = (GrpEnd - GrpBeg) / P0_BLOCK_SIZE + 1
      DO InIter = 1, nInIter
        IF(lPsiUpdt) THEN
          IF (nTracerCntl%lMacro) THEN
            CALL SetMOCPsi_iter(Core, phis, psi, iz)
          ELSE
            STOP
            !CALL PsiUpdate(Core, Fxr, FmInfo%phis, psi, myzb, myze, ng, lxslib, GroupInfo)
          ENDIF
          psi = inv_eigv0 * psi
        ELSE
         lPsiUpdt = .TRUE.
        END IF
        IF (InIter .EQ. nInIter) lJout = .TRUE.
        DO igb = 1, ngBlock
          gb = (igb - 1) * P0_BLOCK_SIZE + GrpBeg
          ge = min(igb * P0_BLOCK_SIZE + GrpBeg - 1, ng)
          ngLocal = ge - gb + 1

          CALL SetTranSrcNM(Core, FmInfo, Fxr, trsrc, phis, TranPhi, psi, ResSrc, xst, iz, gb, ge, &
            GroupInfo, TranInfo, TranCntl, lxslib, PE, gb-1)

          IF (nTracerCntl%lMacro) THEN
            CALL SetMOCSource(Core, Fxr(:, iz), src, phis, psi, AxSrc, xst, 1._8, iz, gb, ge, gb - 1)
          ELSE
            CALL SetRtSrcNM_Cusping(Core, FmInfo, Fxr(:, iz), src, phis, psi, AxSrc, xst, 1._8, iz,                               &
              gb, ge, ng, GroupInfo, l3dim, lxslib, lscat1, .FALSE., PE, gb - 1)
          ENDIF
          !$OMP PARALLEL DO COLLAPSE(2)
          DO ifsr = 1, nFsr
            DO ig = 1, ngLocal
              hostMOC(igb)%xst(ig, ifsr) = xst(ig + gb - 1, ifsr)
              hostMOC(igb)%src(ig, ifsr) = src(ig, ifsr) + trsrc(ig, ifsr)
            ENDDO
          ENDDO
          !$OMP END PARALLEL DO
          DO ig = 1, ngLocal
            hostMOC(igb)%PhiAngIn(:, ig, :) = FmInfo%PhiAngIn(:, :, iz, ig + gb - 1)
          ENDDO
          CALL CUDARayTraceAsync(cuMOC, cuDevice, hostMOC(igb)%phis, hostMOC(igb)%src,                            &
            hostMOC(igb)%xst, hostMOC(igb)%jout, hostMOC(igb)%PhiAngIn, lJout, iz, ngLocal)
        ENDDO
        ierr = cudaStreamSynchronize(cuDevice%myStream)
        DO igb = 1, ngBlock
          gb = (igb - 1) * P0_BLOCK_SIZE + GrpBeg
          ge = min(igb * P0_BLOCK_SIZE + GrpBeg - 1, ng)
          ngLocal = ge - gb + 1
          !$OMP PARALLEL DO COLLAPSE(2)
          DO ifsr = 1, nFsr
            DO ig = 1, ngLocal
              phis(ig + gb - 1, ifsr) = hostMOC(igb)%phis(ig, ifsr)
            ENDDO
          ENDDO
          !$OMP END PARALLEL DO
          !$OMP PARALLEL DO COLLAPSE(2)
          DO iray = 1, nPhiAngSv
            DO ig = 1, ngLocal
              FmInfo%PhiAngIn(:, iray, iz, ig + gb - 1) = hostMOC(igb)%PhiAngIn(:, ig, iray)
            ENDDO
          ENDDO
          !$OMP END PARALLEL DO
          IF (lJout) THEN
            !$OMP PARALLEL DO COLLAPSE(3)
            DO ipin = 1, nxy
              DO isurf = 1, 4
                DO ig = 1, ngLocal
                  Jout(:, ig + gb - 1, isurf, ipin) = hostMOC(igb)%jout(:, ig, isurf, ipin)
                ENDDO
              ENDDO
            ENDDO
            !$OMP END PARALLEL DO
          ENDIF
        ENDDO
        IF (lsSPH) THEN
          IF (GrpBeg .LE. igrese) THEN
            gb = max(igresb, GrpBeg); ge = min(igrese, GrpEnd)
            !$OMP PARALLEL DO COLLAPSE(2)
            DO ifsr = 1, nFsr
              DO ig = gb, ge
                phis(ig, ifsr) = phis(ig, ifsr) * ssphfnm(ig, ifsr, iz)
              ENDDO
            ENDDO
            !$OMP END PARALLEL DO
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    !$OMP PARALLEL DO COLLAPSE(2)
    DO ifsr = 1, nFsr
      DO ig = 1, ng
        FmInfo%phis(ifsr, iz, ig) = phis(ig, ifsr)
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    DO ig = 1, ng
      FmInfo%RadJout(:, :, :, iz, ig) = Jout(:, ig, :, :)
    ENDDO
    IF (.NOT. cuDevice%lRayStatic) CALL DeleteRayVar(cuGeometry%AxType(iz))
  ENDDO

  DEALLOCATE(phis)
  DEALLOCATE(TranPhi)
  DEALLOCATE(Jout)
  DEALLOCATE(xst)
  DEALLOCATE(src)
  DEALLOCATE(trsrc)

  ngBlock = (ng - 1) / P0_BLOCK_SIZE + 1
  DO igb = 1, ngBlock
    DEALLOCATE(hostMOC(igb)%phis)
    DEALLOCATE(hostMOC(igb)%PhiAngIn)
    DEALLOCATE(hostMOC(igb)%jout)
    DEALLOCATE(hostMOC(igb)%xst)
    DEALLOCATE(hostMOC(igb)%src)
  ENDDO
  DEALLOCATE(hostMOC)

ELSE ! Pn Scattering MOC

  ngBlock = (ng - 1) / PN_BLOCK_SIZE + 1
  ALLOCATE(hostMOC(ngBlock))
  DO igb = 1, ngBlock
    ALLOCATE(hostMOC(igb)%phis(PN_BLOCK_SIZE, nFsr))
    ALLOCATE(hostMOC(igb)%phim(nMoment, PN_BLOCK_SIZE, nFsr))
    ALLOCATE(hostMOC(igb)%PhiAngIn(nPolarAngle, PN_BLOCK_SIZE, nPhiAngSv))
    ALLOCATE(hostMOC(igb)%jout(3, PN_BLOCK_SIZE, 4, nxy))
    ALLOCATE(hostMOC(igb)%xst(PN_BLOCK_SIZE, nFsr))
    ALLOCATE(hostMOC(igb)%src(PN_BLOCK_SIZE, nFsr))
    ALLOCATE(hostMOC(igb)%srcm(nMoment, PN_BLOCK_SIZE, nFsr))
  ENDDO

  ALLOCATE(phis(ng, nFsr))
  ALLOCATE(TranPhi(ng, nFsr))
  ALLOCATE(Jout(3, ng, 4, nxy))
  ALLOCATE(xst(ng, nFsr))
  ALLOCATE(src(PN_BLOCK_SIZE, nFsr))
  ALLOCATE(srcm(nMoment, PN_BLOCK_SIZE, nFsr))
  ALLOCATE(trsrc(PN_BLOCK_SIZE, nFsr))

  DO iz = cuDevice%myzb, cuDevice%myze
    IF (.NOT. cuDevice%lRayStatic) CALL CopyRayVar(cuGeometry%AxType(iz))
    !$OMP PARALLEL DO COLLAPSE(2)
    DO ig = 1, ng
      DO ifsr = 1, nFsr
        phis(ig, ifsr) = FmInfo%phis(ifsr, iz, ig)
        TranPhi(ig, ifsr) = FmInfo%TranPhi(ifsr, iz, ig)
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    phim => FmInfo%phim(:, :, :, iz)
    CALL SetRtMacXsNM_Cusping(Core, FmInfo, Fxr(:, iz), xst, iz, ng, lxslib, lTrCorrection, lRST, lsSPH, lsSPHreg, PE)
    CALL PseudoAbsorptionNM(Core, Fxr(:, iz), AxPXS, xst, iz, ng, GroupInfo, l3dim)
    DO i = 1, nGroupInfo
      GrpBeg = 1; GrpEnd = ng
      IF (i .GT. 1) THEN
        GrpBeg = GroupInfo%UpScatRange(1); GrpEnd = GroupInfo%UpScatRange(2)
      ENDIF
      ngBlock = (GrpEnd - GrpBeg) / PN_BLOCK_SIZE + 1
      DO InIter = 1, nInIter
        IF(lPsiUpdt) THEN
          IF (nTracerCntl%lMacro) THEN
            CALL SetMOCPsi_iter(Core, phis, psi, iz)
          ELSE
            STOP
            !CALL PsiUpdate(Core, Fxr, FmInfo%phis, psi, myzb, myze, ng, lxslib, GroupInfo)
          ENDIF
          psi = inv_eigv0 * psi
        ELSE
          lPsiUpdt = .TRUE.
        END IF
        IF (InIter .EQ. nInIter) lJout = .TRUE.
        DO igb = 1, ngBlock
          gb = (igb - 1) * PN_BLOCK_SIZE + GrpBeg
          ge = min(igb * PN_BLOCK_SIZE + GrpBeg - 1, ng)
          ngLocal = ge - gb + 1
          CALL SetTranSrcNM(Core, FmInfo, Fxr, trsrc, phis, TranPhi, psi, ResSrc, xst, iz, gb, ge, &
            GroupInfo, TranInfo, TranCntl, lxslib, PE, gb-1)
          CALL SetRtSrcNM_Cusping(Core, FmInfo, Fxr(:, iz), src, phis, psi, AxSrc, xst, 1._8, iz,                                   &
                          gb, ge, ng, GroupInfo, l3dim, lxslib, lscat1, .FALSE., PE, gb - 1)
          CALL SetRtP1SrcNM(Core, Fxr(:, iz), srcm, phim, xst, iz, gb, ge, ng, GroupInfo,                           &
                            lxslib, nTracerCntl%ScatOd, PE, gb - 1)
          !$OMP PARALLEL DO COLLAPSE(2)
          DO ifsr = 1, nFsr
            DO ig = 1, ngLocal
              hostMOC(igb)%xst(ig, ifsr) = xst(ig + gb - 1, ifsr)
              hostMOC(igb)%src(ig, ifsr) = src(ig, ifsr) + trsrc(ig, ifsr)
              hostMOC(igb)%srcm(:, ig, ifsr) = srcm(:, ig, ifsr)
            ENDDO
          ENDDO
          !$OMP END PARALLEL DO
          !hostMOC(igb)%PhiAngIn = 0.0
          DO ig = 1, ngLocal
            hostMOC(igb)%PhiAngIn(:, ig, :) = FmInfo%PhiAngIn(:, :, iz, ig + gb - 1)
          ENDDO
          CALL CUDAP1RayTrace(cuMOC, cuDevice, hostMOC(igb)%phis, hostMOC(igb)%phim, hostMOC(igb)%xst,              &
                              hostMOC(igb)%src, hostMOC(igb)%srcm, hostMOC(igb)%jout, hostMOC(igb)%PhiAngIn,        &
                              lJout, iz, ngLocal)
        ENDDO
        ierr = cudaStreamSynchronize(cuDevice%myStream)
        DO igb = 1, ngBlock
          gb = (igb - 1) * PN_BLOCK_SIZE + GrpBeg
          ge = min(igb * PN_BLOCK_SIZE + GrpBeg - 1, ng)
          ngLocal = ge - gb + 1
          !$OMP PARALLEL DO COLLAPSE(2)
          DO ifsr = 1, nFsr
            DO ig = 1, ngLocal
              phis(ig + gb - 1, ifsr) = hostMOC(igb)%phis(ig, ifsr)
              phim(:, ig + gb - 1, ifsr) = hostMOC(igb)%phim(:, ig, ifsr)
            ENDDO
          ENDDO
          !$OMP END PARALLEL DO
          !$OMP PARALLEL DO COLLAPSE(2)
          DO iray = 1, nPhiAngSv
            DO ig = 1, ngLocal
              FmInfo%PhiAngIn(:, iray, iz, ig + gb - 1) = hostMOC(igb)%PhiAngIn(:, ig, iray)
            ENDDO
          ENDDO
          !$OMP END PARALLEL DO
          IF (lJout) THEN
            !$OMP PARALLEL DO COLLAPSE(3)
            DO ipin = 1, nxy
              DO isurf = 1, 4
                DO ig = 1, ngLocal
                  Jout(:, ig + gb - 1, isurf, ipin) = hostMOC(igb)%jout(:, ig, isurf, ipin)
                ENDDO
              ENDDO
            ENDDO
            !$OMP END PARALLEL DO
          ENDIF
        ENDDO
        IF (lsSPH) THEN
          IF (GrpBeg .LE. igrese) THEN
            gb = max(igresb, GrpBeg); ge = min(igrese, GrpEnd)
            !$OMP PARALLEL DO COLLAPSE(2)
            DO ifsr = 1, nFsr
              DO ig = gb, ge
                phis(ig, ifsr) = phis(ig, ifsr) * ssphfnm(ig, ifsr, iz)
              ENDDO
            ENDDO
            !$OMP END PARALLEL DO
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    !$OMP PARALLEL DO COLLAPSE(2)
    DO ifsr = 1, nFsr
      DO ig = 1, ng
        FmInfo%phis(ifsr, iz, ig) = phis(ig, ifsr)
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    DO ig = 1, ng
      FmInfo%RadJout(:, :, :, iz, ig) = Jout(:, ig, :, :)
    ENDDO
    IF (.NOT. cuDevice%lRayStatic) CALL DeleteRayVar(cuGeometry%AxType(iz))
  ENDDO

  DEALLOCATE(phis)
  DEALLOCATE(TranPhi)
  DEALLOCATE(Jout)
  DEALLOCATE(xst)
  DEALLOCATE(src)
  DEALLOCATE(srcm)
  DEALLOCATE(trsrc)

  ngBlock = (ng - 1) / PN_BLOCK_SIZE + 1
  DO igb = 1, ngBlock
    DEALLOCATE(hostMOC(igb)%phis)
    DEALLOCATE(hostMOC(igb)%phim)
    DEALLOCATE(hostMOC(igb)%PhiAngIn)
    DEALLOCATE(hostMOC(igb)%jout)
    DEALLOCATE(hostMOC(igb)%xst)
    DEALLOCATE(hostMOC(igb)%src)
    DEALLOCATE(hostMOC(igb)%srcm)
  ENDDO
  DEALLOCATE(hostMOC)

ENDIF

CALL DeallocMOCVar(cuMOC, nTracerCntl%lScat1, cuCntl%lAsync, .NOT. cuCntl%lAsync)


TimeRtEnd = nTracer_dclock(FALSE, FALSE)
TimeRtElapsed = TimeRtEnd - TimeRtBeg
#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_RTMASTER_COMM)
#endif

CALL MPI_MAX_REAL(TimeRtElapsed, PE%MPI_RTMASTER_COMM, .TRUE.)
write(mesg, '(3x, a, f10.3, 2x, a)') 'CUDA MOC Sweep Finished in ', TimeRtElapsed, 'Sec'
IF (PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)

!psid(:, myzb : myze) = psi(:, myzb : myze)
IF (nTracerCntl%lMacro) THEN
  CALL SetMOCPsi(Core, FmInfo%phis, psi)
ELSE
  CALL PsiUpdate(Core, Fxr, FmInfo%phis, psi, myzb, myze, ng, lxslib, GroupInfo)
ENDIF
psi = inv_eigv0 * psi

IF (nTracerCntl%lMacro) THEN
  CALL SetMOCPower(Core, FmInfo%phis, FmInfo%Power)
ELSE
  CALL PowerUpdate(Core, Fxr, FmInfo%phis, FmInfo%Power, myzb, myze, ng, lxslib, GroupInfo, PE)
ENDIF

fiserr = PsiErr(Core, Psi, PsiD, myzb, myze, PE)
eigerr = 0._8
reserr = 0._8
#ifdef MPI_ENV
errdat = (/fiserr, eigerr, reserr/)
CALL BCAST(errdat, 3, PE%MPI_RT_COMM)
fiserr = errdat(1); eigerr = errdat(2); reserr = errdat(3)
#endif

IF(fisErr .LT. TranCntl%psi_conv) ItrCntl%lConv = .TRUE.
WRITE(mesg ,'(A5,I7,F15.6, 2(1pE12.3), a12)')  'RT',itrcntl%mocit, EIGV0, eigerr, fiserr, '   *********'
IF(PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)

MOCTimeEnd = nTracer_dclock(FALSE, FALSE)
TimeChk%MocTime = TimeChk%MocTime + MocTimeEnd - MocTimeBeg
!TimeChk%MocRtTime = TimeChk%MocRtTime + rtTime

TranCntl%delpsifm = fiserr

END SUBROUTINE

SUBROUTINE CUDATranCmfdAcc(Core, CmInfo, FmInfo, TranInfo, TranCntl, lfirst)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,     CmInfo_Type,       FmInfo_Type,      TranCntl_Type,   &
                             TranInfo_Type,     FxrInfo_Type
USE Core_mod,         ONLY : GroupInfo
USE PE_Mod,           ONLY : PE
USE CNTL,             ONLY : nTracerCntl
USE itrcntl_mod,      ONLY : ItrCntl
USE SUBGRP_MOD,       ONLY : FxrChiGen
USE ioutil,           ONLY : message,           terminate
USE files,            ONLY : io8
USE timer,            ONLY : nTracer_dclock,    TimeChk
USE MOC_MOD,          ONLY : GetNeighborMocFlux
USE MOC_COMMON,       ONLY : SetCoreMacXS_Cusping
USE CMFD_COMMON,      ONLY : HomogenizeXS_Cusping,  SetRadialCoupling, HomogenizeKinParam,  SetCMFDPrecCoeff
USE CUDA_INIT,        ONLY : DeallocCMFDVar
USE CUDA_AXIAL,       ONLY : cuTranAxialSolver,     cuSetAxialDtil,    cuSetTranAxialSourceOperator
USE CUDA_CMFD
USE CUDA_Transient
USE CUDA_INIT
USE ieee_arithmetic
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(TranInfo_TYpe) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
LOGICAL :: lfirst

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
REAL, POINTER :: Jout(:, :, :, :, :), AxSrc(:, :, :), AxPXS(:, :, :)
REAL, POINTER :: phis(:, :, :), phim(:, :, :, :), phic(:, :, :)
REAL :: eigv0
REAL :: CmfdTimeBeg, CmfdTimeEnd, CmfdInitBeg, CmfdInitEnd
REAL :: solverRes, solverRes0, err_ratio
REAL :: srcnorm, resnorm, res0, res0in, solres0in
INTEGER :: ng, nprec, nxy
INTEGER :: myzb, myze
INTEGER :: iter, initer, axiter, iprec
LOGICAL :: l3dim, lxslib, lscat1, lDhat, lAxialUpdt, lXSUpdt
REAL :: Kerr

INTEGER :: i

CmfdTimeBeg = nTracer_dclock(FALSE, FALSE)

l3dim = nTracerCntl%l3dim
lxslib = nTracerCntl%lxslib
lscat1 = nTracerCntl%lscat1
Fxr => FmInfo%Fxr
Jout => CmInfo%RadJout
phis => FmInfo%phis
phim => FmInfo%phim
phic => CmInfo%phic
AxSrc => FmInfo%AxSrc
AxPXS => FmInfo%AxPXS

lDhat = ItrCntl%CMFDIt .NE. ItrCntl%CMFDIt0
!IF(.NOT. TranCntl%lMOCUpdt) lDhat = .TRUE.
lDhat = lDhat .AND. TranCntl%lMOCUpdt
lXSUpdt = TranCntl%lMOCUpdt .OR. lfirst
lXSUpdt = .TRUE.
lAxialUpdt = cuCntl%lAxial .AND. l3dim .AND. .NOT. lfirst
ng = cuGeometry%ng
nxy = cuGeometry%nxyc
myzb = cuDevice%myzb
myze = cuDevice%myze
eigv0 = TranInfo%eigv0
nprec = TranInfo%nprec

IF(nTracerCntl%lCusping_MPI .AND. lXSUpdt) THEN
  CALL GetNeighborMocFlux(phis, FmInfo%neighphis, cuGeometry%nFsr, myzb, myze, 1, ng, Core%nz, Core%AxBC)
  IF(nTracerCntl%lMacro .AND. lxslib) CALL SetCoreMacXS_Cusping(Core, FmInfo)
END IF

CALL omp_set_num_threads(PE%nCMFDThread)
CALL omp_set_nested(.TRUE.)
CALL omp_set_max_active_levels(2)

IF(lXSUpdt) THEN
  CmfdInitBeg = nTracer_dclock(FALSE, FALSE)
  WRITE(mesg, '(a)') 'Cell Homogenization...'
  IF (PE%Master) CALL message(io8, .TRUE., .TRUE., mesg)
  CALL FxrChiGen(Core, Fxr, FmInfo, GroupInfo, PE, nTracerCntl)
  CALL HomogenizeXS_Cusping(Core, FmInfo, cuGeometry%superPin, Fxr, cuCMFD%PinXS, phis, ng, nxy, myzb, myze, &
                            PE%nCMFDThread, lxslib, lscat1, .FALSE.)
  IF(lDhat) THEN
    WRITE(mesg, '(a)') 'Set Radial Coupling...'
    IF (PE%Master) CALL message(io8, .TRUE., .TRUE., mesg)
    CALL SetRadialCoupling(cuGeometry%superPin, cuCMFD%PinXS, Jout, ng, nxy, myzb, myze, lDhat)
  END IF

  IF(lfirst .AND. TranCntl%lGuess) THEN
    CALL cuCorrectorGuess(GroupInfo, TranInfo, TranCntl, nTracerCntl, l3dim)
  END IF

  CALL acc_set_device_num(MPI_CUDA_SHARED_RANK, acc_device_nvidia)
  CALL omp_set_num_threads(cuCntl%nCMFDHelper)

  !IF(lDhat) CALL cuSetCMFDPhis(cuCMFD, cuDevice, .FALSE.)
  CALL cuSetCMFDPhis(cuCMFD, cuDevice, .FALSE.)

  CALL HomogenizeKinParam(Core, FmInfo, TranInfo, TranCntl, GroupInfo, cuGeometry%superPin, cuCMFD%PinXS,&
    nxy, myzb, myze, ng, nprec, lxslib)

  CALL SetCMFDPrecCoeff(TranInfo, TranCntl, cuCMFD%PinXS, cuTranCMInfo%CellOmegam, cuTranCMInfo%CellOmega0, cuTranCMInfo%CellOmegap, &
    nxy, myzb, myze, nprec)

  IF (l3dim) CALL cuSetAxialDtil(cuCMFD, cuAxial, cuDevice)
  CmfdInitEnd = nTracer_dclock(FALSE, FALSE)
  TimeChk%CmfdInitTime = TimeChk%CmfdInitTime + (CmfdInitEnd - CmfdInitBeg)
END IF


IF(lAxialUpdt) THEN
  !CALL cuGetNeighborFlux(cuCMFD, cuDevice)
  CALL cuSetTranAxialSourceOperator(Core, FmInfo, GroupInfo, TranInfo, TranCntl, cuCMFD, cuAxial, cuDevice)
  CALL cuTranAxialSolver(TranInfo, TranCntl, eigv0)
ENDIF

CALL cuCopyFlux(cuCMFD, cuDevice, 1)

CmfdInitBeg = nTracer_dclock(FALSE, FALSE)
WRITE(mesg, '(a)') 'Linear System Construction...'
IF (PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)
CALL cuSetTranCMFDSourceOperator(cuCMFD, cuDevice, TranInfo, TranCntl, GroupInfo%lUpscat)
!CALL cuSetNaturalTranBiCGSystem(cuCMFD, cuDevice, TranInfo, TranCntl, l3dim, TranCntl%lMOCUpdt .AND. lfirst, .TRUE.)
CALL cuSetNaturalTranBiCGSystem(cuCMFD, cuDevice, TranInfo, TranCntl, l3dim, TranCntl%lMOCUpdt, .TRUE.)
!CALL cuSetNaturalTranBiCGSystem(cuCMFD, cuDevice, TranInfo, TranCntl, l3dim, .FALSE., .TRUE.)
CmfdInitEnd = nTracer_dclock(FALSE, FALSE)
TimeChk%CmfdInitTime = TimeChk%CmfdInitTime + (CmfdInitEnd - CmfdInitBeg)

WRITE(mesg, '(a)') 'Performing CMFD Acceleration...'
IF (PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)
ItrCntl%Cmfdit = ItrCntl%Cmfdit + 1
IF(TranCntl%lchidk) THEN
  DO iprec = 1, nprec
    IF(lXSUpdt) THEN
      CALL cuSetPrecSrcOperator(iprec, TranInfo%lambda(iprec))
    END IF
    CALL cuCMFDPrecSrcKUpdt(TranInfo, TranCntl, iprec)
  END DO
ELSE
  IF(lXSUpdt) THEN
    CALL cuSetPrecSrcOperator()
  END IF
  CALL cuCMFDPrecSrcUpdt(TranInfo, TranCntl)
END IF

CALL cuTranCMFDSrcUpdt(TranCntl)
DO iter = 1, 20
  IF(iter .EQ. 1) res0 = cuTranSrcNorm(eigv0)
  CALL cuInnerSolve_outInfo(cuCMFD, cuDevice, 0.1, solverRes, solverRes0, initer)
  ItrCntl%innerit = ItrCntl%innerit + initer
  srcnorm = cuTranSrcNorm(eigv0)
  resnorm = solverRes / srcnorm
      IF(ieee_is_nan(resnorm)) STOP 'TRAN CMFD DIverge!!'
  IF(iter .EQ. 1) THEN
    res0 = solverRes0 / res0
    solres0in = solverRes0
  END IF
  err_ratio = resnorm / res0
  IF (PE%MASTER) WRITE(mesg, '(a9, i9, i11, 3x, f10.5, 1p, 2e13.3)') 'MGOUTER', ItrCntl%CMFDIt, initer, err_ratio, resnorm, res0
  IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)

  IF(resnorm .LT. 1.E-10) EXIT
  IF(solverRes/solres0in .LT. 0.1) EXIT
END DO

IF (lAxialUpdt) THEN
  DO axiter = 1, 2
    CALL cuCopyFlux(cuCMFD, cuDevice, 2)
    CALL cuTranAxialSolver(TranInfo, TranCntl, eigv0)
    CmfdInitBeg = nTracer_dclock(FALSE, FALSE)
    WRITE(mesg, '(a)') 'Linear System Construction...'
    IF (PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)
    ItrCntl%Cmfdit = ItrCntl%Cmfdit + 1
    CALL cuSetNaturalTranBiCGSystem(cuCMFD, cuDevice, TranInfo, TranCntl, l3dim, .FALSE., .TRUE.)
    CmfdInitEnd = nTracer_dclock(FALSE, FALSE)
    TimeChk%CmfdInitTime = TimeChk%CmfdInitTime + (CmfdInitEnd - CmfdInitBeg)
    DO iter = 1, 20
      IF(iter .EQ. 1) res0in = cuTranSrcNorm(eigv0)
      CALL cuInnerSolve_outInfo(cuCMFD, cuDevice, 0.1, solverRes, solverRes0, initer)
      ItrCntl%innerit = ItrCntl%innerit + initer
      srcnorm = cuTranSrcNorm(eigv0)
      resnorm = solverRes / srcnorm
      IF(ieee_is_nan(resnorm)) STOP 'TRAN CMFD DIverge!!'
      IF(iter .EQ. 1) THEN
        res0in = solverRes0 / res0in
        solres0in = solverRes0
      END IF
      err_ratio = resnorm / res0
      IF (PE%MASTER) WRITE(mesg, '(a9, i9, i11, 3x, f10.5, 1p, 2e13.3)') 'MGOUTER', ItrCntl%CMFDIt, initer, err_ratio, resnorm, res0in
      IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)

      IF(resnorm .LT. 1.E-10) EXIT
      IF(solverRes/solres0in .LT. 0.1) EXIT
    END DO
  END DO
ENDIF

IF(.NOT. TranCntl%lMOCUpdt .AND. resnorm .LT. TranCntl%cmfd_res_conv) ItrCntl%lConv = .TRUE.
IF(TranCntl%lPCQSIter .AND. .NOT. ItrCntl%lConv) THEN
  CALL cuCorrectorIter(GroupInfo, TranInfo, TranCntl, nTracerCntl, l3dim)
END IF

CALL cuCopyFlux(cuCMFD, cuDevice, 2)
!IF(TranCntl%lMOCUpdt .OR. ItrCntl%lconv) THEN
IF(lXSUpdt .OR. ItrCntl%lconv) THEN
  CALL cuSetMOCPhis(Core, cuCMFD, cuDevice, phis, phim, phic, lScat1)
  CALL cuSetMOCPhiIn(Core, CmInfo, myzb, myze, ng)
ELSE
  CALL cuSetCoarsePhis(Core, cuCMFD, cuDevice, phic)
END IF

IF (l3dim) THEN
  CALL cuGetNeighborFlux(cuCMFD, cuDevice)
  CALL cuSetAxialSrc(cuCMFD, cuDevice, AxSrc, AxPXS, phic)
ENDIF

!CALL DeallocTranCMFD(cuCMFD, cuTranCMInfo)

CALL destroyCsr(cuCMFD%F)
CALL destroyCsr(cuCMFD%Chi)

CmfdTimeEnd = nTracer_dclock(FALSE, FALSE)
TimeChk%CmfdTime = TimeChk%CmfdTime + (CmfdTimeEnd - CmfdTimeBeg)

TranCntl%cmfdres = resnorm

END SUBROUTINE

SUBROUTINE CUDATransient_Driver()
USE PARAM
USE geom,               ONLY : Core
USE RAYS,               ONLY : RayInfo
USE Core_mod,           ONLY : FmInfo,            CmInfo,             THInfo,           GroupInfo,   &
                               eigv
USE PE_Mod,             ONLY : PE
USE CNTL,               ONLY : nTracerCntl
USE itrcntl_mod,        ONLY : ItrCntl
USE SUBGRP_MOD,         ONLY : SubGrpEffXsGen
USE TRAN_MOD,           ONLY : TranInfo,          TranCntl,                                         &
                               InitTransient,     KinParamGen,        InitPrecursor,                &
                               EDIT_NNSAMPLING
USE Boron_mod,          ONLY : SetBoronCoolant
USE MOC_COMMON,         ONLY : SetCoreMacXs
USE CMFD_COMMON,        ONLY : SetRadialCoupling, HomogenizeKinParam,  SetCMFDPrecCoeff
USE XsPerturb_mod,      ONLY : XsPerturbation
USE files,              ONLY : io8
USE ioutil,             ONLY : message,           terminate
USE CUDA_MASTER
USE CUDA_INIT,          ONLY : AllocTransientVar, deallocTranCMFD
USE CUDA_Transient
USE PromptFeedback_mod, ONLY : UpdtFuelTemp, InitPromptFeedback
USE TH_Mod,             ONLY : TransientTH,       THFeedback,         SaveTranThSol
USE CUDA_PWDIST,        ONLY : cuUpdthPrec
IMPLICIT NONE
INTEGER :: nxy, myzb, myze, nprec, ng
INTEGER :: istep, i
LOGICAL :: l3dim, lxslib, lout

INTEGER :: ierr

nxy = cuGeometry%nxyc
myzb = cuDevice%myzb
myze = cuDevice%myze
ng = cuGeometry%ng
nprec = cuGeometry%nprec
l3dim = nTracerCntl%l3dim
lxslib = nTracerCntl%lxslib

!nTracerCntl%lfeedback = .FALSE.
!nTracerCntl%lMacro = .FALSE.
IF(nTracerCntl%lDynamicBen .OR. nTracerCntl%lKineticBen .OR. nTracerCntl%libtyp .EQ. 11) THEN
  nTracerCntl%lchidgen = .FALSE.
  nTracerCntl%lchidkgen = .FALSE.
END IF
TranCntl%lchidk = nTracerCntl%lchidkgen

WRITE(mesg, '(a)') 'Initialize Transient Calculation...'
IF (PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)
TranCntl%NowStep = 1
CALL AllocTransientVar(cuDevice, PE%lCUDACMFD, l3dim, TranCntl%lchidk)
CALL cuInitTransient(Core, FmInfo, CmInfo, RayInfo, ThInfo, GroupInfo, nTracerCntl, TranInfo, TranCntl, eigv, l3dim)
IF(TranCntl%lDynamicBen) CALL InitPromptFeedback(Core, FmInfo, CmInfo, TranInfo, PE, ng)

CALL HomogenizeKinParam(Core, FmInfo, TranInfo, TranCntl, GroupInfo, cuGeometry%superPin, cuCMFD%PinXS,&
                        nxy, myzb, myze, ng, nprec, lxslib)
CALL SetCMFDPrecCoeff(TranInfo, TranCntl, cuCMFD%PinXS, cuTranCMInfo%CellOmegam, cuTranCMInfo%CellOmega0, cuTranCMInfo%CellOmegap, &
                      nxy, myzb, myze, nprec)
CALL cuUpdtResSrc(Core, FmInfo, cuCMFD, cuDevice, TranInfo, TranCntl, l3dim)

CALL cuTranReactivityUpdt_Adj(GroupInfo, TranInfo, TranCntl, l3dim, .TRUE.)
!CALL cuTranReactivityUpdt(GroupInfo, TranInfo, l3dim, .TRUE.)

IF(TranCntl%lCorrector .OR. TranCntl%lGuess) CALL cuSavePKEParameters(TranInfo, TranCntl)
IF(PE%MASTER) THEN
  IF(.NOT. nTracerCntl%lAdjoint) THEN
    !WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, 1e11.4, x, a, 1e10.3, x, a, 1e12.5)') '@1', 0.,"s, P=", TranInfo%PowerLevel * 100, ', $= ', TranInfo%reactivity, ', G=', TranInfo%lifetime
    !CALL message(io8, TRUE, TRUE, MESG)
    !WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, F7.1, x, a, F8.5, x, a, 1e10.3)') '@2', 0.,"s, r(pcm)= ", TranInfo%delrho/10._8, ', k= ', TranInfo%TranEig/10._8, ', B= ', TranInfo%CoreBeta
    !CALL message(io8, TRUE, TRUE, MESG)
    WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, 1e11.4, x, a, 1e10.3, x, a, 1e12.5)') '@1', 0.,"s, P=", TranInfo%PowerLevel * 100, ', $= ', TranInfo%reactivity_Dynamic, ', G=', TranInfo%lifetime_Dynamic
    CALL message(io8, TRUE, TRUE, MESG)
    WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, F7.1, x, a, F8.5, x, a, 1e10.3)') '@2', 0.,"s, r(pcm)= ", TranInfo%delrho_Dynamic/10._8, ', k= ', TranInfo%TranEig_Dynamic/10._8, ', B= ', TranInfo%CoreBeta_Dynamic
    CALL message(io8, TRUE, TRUE, MESG)
    IF(nTracerCntl%lFeedback) THEN
      WRITE(mesg, '(a, 1p, e12.4, x, a, f10.3, x, a, f10.3, x, a)') '@5', TranCntl%T(istep),"s, Max Tf= ", ThInfo%TfMax/10, "C, Max TfAvg= ", ThInfo%Tfavg_max/10, "C"
      CALL message(io8, TRUE, TRUE, MESG)
      WRITE(mesg, '(a, 1p, e12.4, x, a, f10.3, x, a, f10.3, x, a)') '@6', TranCntl%T(istep),"s, Max Tcl= ", ThInfo%Tfcl_max/10, "C, Max TMOA= ", ThInfo%TModoutAvg/10, "C"
      CALL message(io8, TRUE, TRUE, MESG)
      WRITE(mesg, '(a, 1p, e12.4, x, a, f10.3)') '@7', TranCntl%T(istep),"s, Max Pow= ", ThInfo%pwpeakf/10
      CALL message(io8, TRUE, TRUE, MESG)
    ENDIF
  ELSE
    WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, 1e11.4, x, a, 1e10.3, x, a, 1e12.5)') '@3', 0.,"s, P=", TranInfo%PowerLevel * 100, ', $= ', TranInfo%reactivity_Dynamic, ', G=', TranInfo%lifetime_dynamic
    CALL message(io8, TRUE, TRUE, MESG)
    WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, F7.1, x, a, F8.5, x, a, 1e10.3)') '@4', 0.,"s, r(pcm)= ", &
      TranInfo%delrho_Dynamic/10._8, ', k= ', TranInfo%TranEig_Dynamic/10._8, ', B= ', TranInfo%CoreBeta_Dynamic
    CALL message(io8, TRUE, TRUE, MESG)
    IF(nTracerCntl%lFeedback) THEN
      WRITE(mesg, '(a, 1p, e12.4, x, a, f10.3, x, a, f10.3, x, a)') '@5', TranCntl%T(istep),"s, Max Tf= ", ThInfo%TfMax/10, "C, Max TfAvg= ", ThInfo%Tfavg_max/10, "C"
      CALL message(io8, TRUE, TRUE, MESG)
      WRITE(mesg, '(a, 1p, e12.4, x, a, f10.3, x, a, f10.3, x, a)') '@6', TranCntl%T(istep),"s, Max Tcl= ", ThInfo%Tfcl_max/10, "C, Max TMOA= ", ThInfo%TModoutAvg/10, "C"
      CALL message(io8, TRUE, TRUE, MESG)
      WRITE(mesg, '(a, 1p, e12.4, x, a, f10.3)') '@7', TranCntl%T(istep),"s, Max Pow= ", ThInfo%pwpeakf/10
      CALL message(io8, TRUE, TRUE, MESG)
    ENDIF
  ENDIF
  WRITE(mesg, '(A)') hbar2(1:77)
  CALL message(io8, FALSE, TRUE, MESG)
ENDIF

DO istep = 1, TranCntl%nStep
  nTracerCntl%CalNoId = nTracerCntl%CalNoId + 1

  CALL XsPerturbation(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, PE)
  IF(nTracerCntl%BoronPPM .GT. 0.) CALL SetBoronCoolant(Core, FmInfo%Fxr, nTracerCntl%BoronPPM , PE%myzb, PE%myze)
  CALL KinParamGen(Core, FmInfo, TranInfo, ThInfo, GroupInfo, .TRUE., nTracerCntl, PE)

  IF(nTracerCntl%lDcyHeat) CALL cuUpdthPrec(cuPwDist, TranCntl%DelT(TranCntl%nowstep))

  IF(nTracerCntl%lFeedback) THEN
    CALL TransientTH(Core, CmInfo, FmInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
    CALL THFeedBack(Core, CmInfo, FmInfo, ThInfo, nTRACERCntl, PE)
  END IF
  !Feedback for C5G7-TD DynamicBenchmark
  IF(TranCntl%lDynamicBen) THEN
    CALL UpdtFuelTemp(Core, CmInfo, TranInfo, TranCntl, PE, ng, .TRUE.)
  END IF

  IF(TranCntl%lCondiMOC) THEN
    CALL CheckCondMOC(Core, FmInfo, CmInfo, TranInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
  ENDIF

  IF(nTracerCntl%lFeedback .AND. nTracerCntl%lrestrmt) THEN
    CALL CUDASubGrpFsp(Core, FmInfo,THInfo, RayInfo, GroupInfo)
  ENDIF

  ! Transient Fixed Source Problem
  CALL CUDATransientFsp_Driver(Core, FmInfo, CmInfo, THInfo, TranInfo, GroupInfo, nTracerCntl, TranCntl, PE)

  ! Homogenize Kinetic parameter
  CALL HomogenizeKinParam(Core, FmInfo, TranInfo, TranCntl, GroupInfo, cuGeometry%superPin, cuCMFD%PinXS,&
                          nxy, myzb, myze, ng, nprec, lxslib)
  CALL SetCMFDPrecCoeff(TranInfo, TranCntl, cuCMFD%PinXS, cuTranCMInfo%CellOmegam, cuTranCMInfo%CellOmega0, cuTranCMInfo%CellOmegap, &
                        nxy, myzb, myze, nprec)
  ! Predictor-Corrector method
  IF(TranCntl%lCorrector) THEN
    CALL cuTranReactivityUpdt_Adj(GroupInfo, TranInfo, TranCntl, l3dim, .FALSE.)
    CALL Transient_Corrector(TranInfo, TranCntl, nTracerCntl)
    !CALL Transient_Corrector(TranInfo, TranCntl, nTracerCntl, TranCntl%T(TranCntl%nowstep-1), TranCntl%T(TranCntl%nowstep), TranCntl%Tdiv_corrector)
    CALL cuCorrectNowStep(Core, CmInfo, FmInfo, GroupInfo, TranInfo, TranCntl, nTracerCntl)
  ELSE
    CALL cuTranReactivityUpdt_Adj(GroupInfo, TranInfo, TranCntl, l3dim, .FALSE.)
    WRITE(mesg, '(a15,ES14.7)') 'Amplitude',TranInfo%lifetime_Dynamic * TranInfo%Factor_F * TranInfo%Inv_Factor_K0
    IF (PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)
  END IF

  !Theta method source term
  CALL cuUpdtResSrc(Core, FmInfo, cuCMFD, cuDevice, TranInfo, TranCntl, l3dim)
  !Precursor Update
  CALL cuUpdtPrec(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, nprec)
  !Calculate
  CALL cuUpdtPowerLevel(TranInfo, .FALSE.,cuCntl%lPwDist,nTracerCntl%lDcyHeat)
  nTracerCntl%PowerLevel = TranInfo%PowerLevel; ThInfo%PowLv = nTracerCntl%PowerLevel

  CALL cuSaveTranSol(FmInfo, nTracerCntl)
  IF(nTracerCntl%lFeedback) THEN
    CALL SaveTranTHsol(Core, THInfo, TranCntl, nTracerCntl, PE)
  END IF
  IF(TranCntl%lCorrector .OR. TranCntl%lGuess) CALL cuSavePKEParameters(TranInfo, TranCntl)

  CALL ChkOutputWrite(lout, TranCntl)
  IF(lout) THEN
    CALL OutputEdit()
    !CALL PrintHeatPower(Core, CmInfo, FmInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
    TranCntl%NowWriteStep = TranCntl%NowWriteStep + 1
  ENDIF
  IF(TranCntl%lNNSampling) CALL EDIT_NNSAMPLING(Core, TranCntl, CmInfo, PE, ng)
  IF(PE%MASTER) THEN
    IF(.NOT. nTracerCntl%lAdjoint) THEN
      WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, 1e11.4, x, a, 1e10.3, x, a, 1e12.5)') '@1', TranCntl%T(istep),"s, P=", TranInfo%PowerLevel * 100, ', $= ', TranInfo%reactivity_Dynamic, ', G=', TranInfo%lifetime_Dynamic
      CALL message(io8, TRUE, TRUE, MESG)
      WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, F7.1, x, a, F8.5, x, a, 1e10.3)') '@2', TranCntl%T(istep),"s, r(pcm)= ", TranInfo%delrho_Dynamic/10._8, ', k= ', TranInfo%TranEig_Dynamic/10._8, ', B= ', TranInfo%CoreBeta_Dynamic
      CALL message(io8, TRUE, TRUE, MESG)
    ELSE
      WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, 1e14.7, x, a, 1e14.7, x, a, 1e14.7)') '@3', TranCntl%T(istep),"s, P=", TranInfo%PowerLevel * 100, ', $= ', TranInfo%reactivity_Dynamic, ', G=', TranInfo%lifetime_dynamic
      CALL message(io8, TRUE, TRUE, MESG)
      WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, F7.1, x, a, F8.5, x, a, 1e14.7)') '@4', TranCntl%T(istep),"s, r(pcm)= ", &
        TranInfo%delrho_Dynamic/10._8, ', k= ', TranInfo%TranEig_Dynamic/10._8, ', B= ', TranInfo%CoreBeta_Dynamic
      CALL message(io8, TRUE, TRUE, MESG)
    ENDIF
    IF(nTracerCntl%lFeedback) THEN
      WRITE(mesg, '(a, 1p, e12.4, x, a, f10.3, x, a, f10.3, x, a)') '@5', TranCntl%T(istep),"s, Max Tf= ", ThInfo%TfMax/10, "C, Max TfAvg= ", ThInfo%Tfavg_max/10, "C"
      CALL message(io8, TRUE, TRUE, MESG)
      WRITE(mesg, '(a, 1p, e12.4, x, a, f10.3, x, a, f10.3, x, a)') '@6', TranCntl%T(istep),"s, Max Tcl= ", ThInfo%Tfcl_max/10, "C, Max TMOD= ", ThInfo%TModoutAvg/10, "C"
      CALL message(io8, TRUE, TRUE, MESG)
      WRITE(mesg, '(a, 1p, e12.4, x, a, f10.3)') '@7', TranCntl%T(istep),"s, Max Pow= ", ThInfo%pwpeakf/10
      CALL message(io8, TRUE, TRUE, MESG)
    ENDIF
    IF(TranCntl%lSCM) THEN
      WRITE(mesg, '(a, 1p, e12.4, x, a, f10.3, a, f10.3)') '@7', TranCntl%T(istep),"s, AmpFrq= ", TranCntl%AmpFrq/10., "  Avg. Val. = ", TranCntl%AvgAmpFrq/10.
      CALL message(io8, TRUE, TRUE, MESG)
      WRITE(mesg, '(A, 1p, e12.4, x, a, f10.6)') '@8', TranCntl%T(istep),"s, Estimated Dynamic k-eff= ", TranCntl%EigD/10.
      CALL message(io8, TRUE, TRUE, MESG)
    ENDIF
    WRITE(mesg, '(A, 1p, e12.4, x, a, i8)') '@9', TranCntl%T(istep), "s, MOC = ", ItrCntl%Mocit - ItrCntl%Mocit0
    CALL message(io8, TRUE, TRUE, MESG)
    WRITE(mesg, '(A, 1p, e12.4, x, a, i8, x, a, i8)') '@10', TranCntl%T(istep), "s, MGCMFD = ", ItrCntl%Cmfdit - ItrCntl%Cmfdit0, "Inner = ", ItrCntl%InnerIt
    CALL message(io8, TRUE, TRUE, MESG)
    WRITE(mesg, '(A, 1p, e12.4, x, a, i8)') '@11', TranCntl%T(istep), "s, GCCMFD = ", ItrCntl%GcCmfdIt - ItrCntl%GcCmfdit0
    CALL message(io8, TRUE, TRUE, MESG)
    WRITE(mesg, '(A)') hbar2(1:77)
    CALL message(io8, FALSE, TRUE, MESG)
  ENDIF
  TranCntl%NowStep = TranCntl%NowStep + 1
  TranInfo%PowerLeveld = TranInfo%PowerLevel
  TranCntl%lStepApprox = .FALSE.
END DO

END SUBROUTINE

SUBROUTINE CUDATransientFsp_Driver(Core, FmInfo, CmInfo, THInfo, TranInfo, GroupInfo, nTracerCntl, TranCntl, PE)
USE PARAM
USE TYPEDEF,            ONLY : CoreInfo_Type,     FmInfo_Type,      CmInfo_Type,      THInfo_Type,    &
                               TranInfo_Type,     GroupInfo_Type,   TranCntl_Type,    PE_Type
USE CNTL,               ONLY : nTracerCntl_Type
USE RAYS,               ONLY : RayInfo
USE itrcntl_mod,        ONLY : ItrCntl
USE SUBGRP_MOD,         ONLY : SubGrpEffXsGen
USE MOC_COMMON,         ONLY : SetCoreMacXS,      SetCoreMacXS_Cusping
USE CUDA_Transient,     ONLY : cuExpTrsfUpdt,     cuUpdtPsiNowStep, cuUpdtPowerLevel
USE files,              ONLY : io8
USE ioutil,             ONLY : message,           terminate
USE TH_Mod,             ONLY : TransientTH,       THFeedback

USE CUDA_MASTER
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(THInfo_Type) :: THInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(TranCntl_Type) :: TranCntl
TYPE(PE_Type) :: PE

INTEGER :: iter
LOGICAL :: lMOC, lmacxs

INTEGER :: ierr, i
TYPE(cudaDeviceProp) :: cudaProperty
INTEGER(8) :: totalByte, allocByte, freeByte(2)
CHARACTER(10) :: totalMBChar, allocMBChar, freeMBChar
REAL :: totalMB, allocMB, freeMB

ItrCntl%lconv = .FALSE.
CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, TranInfo%eigv0, GroupInfo, nTracerCntl, PE)
IF (nTracerCntl%lMacro) THEN
  IF(.NOT. nTracerCntl%lCusping_MPI) CALL SetCoreMacXs(Core, FmInfo)
  lmacxs = .TRUE.
END IF

CALL InitIterVar(Core, FMInfo, CmInfo, GroupInfo,  .FALSE., ItrCntl, nTracerCntl, PE)
ItrCntl%innerit = 0
lMOC = TranCntl%lMOCUpdt

CALL CUDATranCMFDAcc(Core, CmInfo, FmInfo, TranInfo, TranCntl, .TRUE.)
IF(TranCntl%lExpTrsf) THEN
  WRITE(mesg, '(a)') 'Exponential Transformation Factor Update...'
  IF (PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)
  CALL cuExpTrsfUpdt(Core, TranInfo, TranCntl, nTracerCntl%l3dim)
END IF

ItrCntl%srcit = 0
DO iter = 1, TranCntl%nMaxOuter
  ItrCntl%SrcIt = ItrCntl%SrcIt + 1
  !IF(nTracerCntl%lFeedback .AND. TranCntl%T(TranCntl%nowstep) .LT. 0.2) THEN
  IF(nTracerCntl%lFeedback) THEN
    CALL cuUpdtPowerLevel(TranInfo, .FALSE., cuCntl%lPwDist, .FALSE.)
    nTracerCntl%PowerLevel = TranInfo%PowerLevel; ThInfo%PowLv = nTracerCntl%PowerLevel
    CALL TransientTH(Core, CmInfo, FmInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
    CALL THFeedBack(Core, CmInfo, FmInfo, ThInfo, nTRACERCntl, PE)
    IF(.NOT. TranCntl%lCondiMOC) THEN
      IF(THInfo%TDopChg .GT. 1.E-3_8) THEN
        IF(nTracerCntl%lrestrmt) CALL CUDASubGrpFsp(Core, FmInfo,THInfo, RayInfo, GroupInfo)
        CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, TranInfo%Eigv0, GroupInfo, nTracerCntl, PE)
      ENDIF
      !IF (lmacxs) CALL SetCoreMacXs(Core, FmInfo)
    ELSE
      IF(TranCntl%lSgFspUpdt .AND. THInfo%TDopChg .GT. 1.E-3_8) THEN
        IF(nTracerCntl%lrestrmt) CALL CUDASubGrpFsp(Core, FmInfo,THInfo, RayInfo, GroupInfo)
        CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, TranInfo%Eigv0, GroupInfo, nTracerCntl, PE)
      ENDIF
      !IF (lmacxs) CALL SetCoreMacXs(Core, FmInfo)
    ENDIF
    IF (lmacxs .AND. .NOT. nTracerCntl%lCusping_MPI) CALL SetCoreMacXs(Core, FmInfo)
  END IF
  IF(lMOC) THEN
    CALL CUDATranMOCSweep(Core, RayInfo, FmInfo, TranInfo, TranCntl)
  END IF

  !ierr = cudaMemGetInfo(freeByte(2), totalByte)
  !allocByte = totalByte - freeByte(2)
  !DO i = 0, NUM_CUDA_PROC - 1
  !  IF (i .EQ. MPI_CUDA_RANK) THEN
  !    totalMB = DBLE(totalByte) / 1024.0 ** 2
  !    allocMB = DBLE(allocByte) / 1024.0 ** 2
  !    freeMB = DBLE(freeByte(2)) / 1024.0 ** 2
  !    ierr = cudaGetDeviceProperties(cudaProperty, MPI_CUDA_SHARED_RANK)
  !    WRITE(totalMBChar, '(f8.2)') totalMB; totalMBChar = ADJUSTL(totalMBChar)
  !    WRITE(allocMBChar, '(f8.2)') allocMB; allocMBChar = ADJUSTL(allocMBChar)
  !    WRITE(freeMBChar, '(f8.2)') freeMB; freeMBChar = ADJUSTL(freeMBChar)
  !    WRITE(*, '(15x,i4,2x, a)'),i, TRIM(cudaProperty%name) // ' : ' // TRIM(totalMBChar) // 'MB Total, ' //                  &
  !      TRIM(allocMBChar) // 'MB Allocated, ' // TRIM(freeMBChar) // 'MB Free'
  !  ENDIF
  !  CALL MPI_BARRIER(MPI_CUDA_COMM, ierr)
  !ENDDO

  CALL CUDATranCmfdAcc(Core, CmInfo, FmInfo, TranInfo, TranCntl, .FALSE.)

  !ierr = cudaMemGetInfo(freeByte(2), totalByte)
  !allocByte = totalByte - freeByte(2)
  !DO i = 0, NUM_CUDA_PROC - 1
  !  IF (i .EQ. MPI_CUDA_RANK) THEN
  !    totalMB = DBLE(totalByte) / 1024.0 ** 2
  !    allocMB = DBLE(allocByte) / 1024.0 ** 2
  !    freeMB = DBLE(freeByte(2)) / 1024.0 ** 2
  !    ierr = cudaGetDeviceProperties(cudaProperty, MPI_CUDA_SHARED_RANK)
  !    WRITE(totalMBChar, '(f8.2)') totalMB; totalMBChar = ADJUSTL(totalMBChar)
  !    WRITE(allocMBChar, '(f8.2)') allocMB; allocMBChar = ADJUSTL(allocMBChar)
  !    WRITE(freeMBChar, '(f8.2)') freeMB; freeMBChar = ADJUSTL(freeMBChar)
  !    WRITE(*, '(15x,i4,2x, a)'),i, TRIM(cudaProperty%name) // ' : ' // TRIM(totalMBChar) // 'MB Total, ' //                  &
  !      TRIM(allocMBChar) // 'MB Allocated, ' // TRIM(freeMBChar) // 'MB Free'
  !  ENDIF
  !  CALL MPI_BARRIER(MPI_CUDA_COMM, ierr)
  !ENDDO

  IF(TranCntl%lExpTrsf) THEN
    WRITE(mesg, '(a)') 'Exponential Transformation Factor Update...'
    IF (PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)
    CALL cuExpTrsfUpdt(Core, TranInfo, TranCntl, nTracerCntl%l3dim)
  END IF

  IF(iter .LT. 2) CYCLE
  IF(ItrCntl%lConv) EXIT
END DO
IF(TranCntl%lMocUpdt .AND. TranCntl%lCondiMOC) THEN
  CALL UpdtBaseXsCondiMOC(Core, FmInfo, CmInfo, TranInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
ENDIF
CALL cuUpdtPsiNowStep(Core, FmInfo, GroupInfo, nTracerCntl, TranInfo%eigv0)

  !ierr = cudaMemGetInfo(freeByte(2), totalByte)
  !allocByte = totalByte - freeByte(2)
  !DO i = 0, NUM_CUDA_PROC - 1
  !  IF (i .EQ. MPI_CUDA_RANK) THEN
  !    totalMB = DBLE(totalByte) / 1024.0 ** 2
  !    allocMB = DBLE(allocByte) / 1024.0 ** 2
  !    freeMB = DBLE(freeByte(2)) / 1024.0 ** 2
  !    ierr = cudaGetDeviceProperties(cudaProperty, MPI_CUDA_SHARED_RANK)
  !    WRITE(totalMBChar, '(f8.2)') totalMB; totalMBChar = ADJUSTL(totalMBChar)
  !    WRITE(allocMBChar, '(f8.2)') allocMB; allocMBChar = ADJUSTL(allocMBChar)
  !    WRITE(freeMBChar, '(f8.2)') freeMB; freeMBChar = ADJUSTL(freeMBChar)
  !    WRITE(*, '(15x,i4,2x, a)'),i, TRIM(cudaProperty%name) // ' : ' // TRIM(totalMBChar) // 'MB Total, ' //                  &
  !      TRIM(allocMBChar) // 'MB Allocated, ' // TRIM(freeMBChar) // 'MB Free'
  !  ENDIF
  !  CALL MPI_BARRIER(MPI_CUDA_COMM, ierr)
  !ENDDO
  !IF(TranCntl%lExpTrsf) THEN
  !  WRITE(mesg, '(a)') 'Exponential Transformation Factor Update...'
  !  IF (PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)
  !  CALL cuExpTrsfUpdt(Core, TranInfo, TranCntl, nTracerCntl%l3dim)
  !END IF
END SUBROUTINE

SUBROUTINE CUDACmfd_Adj(Core, CmInfo, eigv)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,       Pin_Type,           Cell_Type, CmInfo_Type
USE PE_Mod,         ONLY : PE
USE CNTL,           ONLY : nTracerCntl
USE TRAN_MOD,       ONLY : TranCntl
USE IOUTIL,         ONLY : message
USE FILES,          ONLY : io8
USE CUDA_INIT,      ONLY : DeallocCMFDVar
USE CUDA_CMFD
USE TYPEDEF_COMMON, ONLY : superPin_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
REAL :: eigv

TYPE(superPin_Type), POINTER :: superPin(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
REAL, POINTER :: hz(:), hzfm(:)
REAL :: phisum(cuCMFD%ng)
REAL :: outRes
INTEGER, POINTER :: pinMap(:), fmRange(:,:)
INTEGER :: outIter
INTEGER :: ng, nxy, myzb, myze
INTEGER :: i, j, ixy, ixy_map, ipin, ig, iz, izf

INTEGER :: ierr

Pin => Core%Pin
Cell => Core%CellInfo
superPin => cuGeometry%superPin
ng = cuCMFD%ng
nxy = cuGeometry%nxyc
myzb = cuDevice%myzb
myze = cuDevice%myze
pinMap => cuGeometry%pinMap
fmRange => cuGeometry%fmRange
hz => cuGeometry%hz
hzfm => cuGeometry%hzfm

CALL omp_set_num_threads(PE%nCMFDThread)
CALL omp_set_nested(.TRUE.)
CALL omp_set_max_active_levels(2)

CALL acc_set_device_num(MPI_CUDA_SHARED_RANK, acc_device_nvidia)
CALL omp_set_num_threads(cuCntl%nCMFDHelper)

CALL cuInitAdjPhis(cuCMFD_adj, cuCMFD, .FALSE.)
CALL cuCopyFlux(cuCMFD_adj, cuDevice, 1)

CALL cuCMFDAdjSrcUpdt(cuCMFD_adj, eigv)

outIter = 0
DO WHILE(.TRUE.)
  outIter = outIter + 1
  CALL cuInnerSolve(cuCMFD_adj, cuDevice, 0.01)
  CALL cuCMFDAdjSrcUpdt(cuCMFD_adj, eigv)
  CALL cuCMFDAdjEigUpdt(cuCMFD_adj, eigv)
  outRes = cuCMFDResidual(cuCMFD_adj, cuDevice)
  IF (PE%MASTER) WRITE(mesg, '(a17, i9, f24.6, 3x, 1p, e15.3)') 'MGOUTER_ADJOINT', outIter, eigv, outRes
  IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
  !IF(outRes .LT. TranCntl%cmfd_res_conv) EXIT
  IF(outRes .LT. 1.E-6) EXIT

  !IF(mod(outIter, 5) .NE. 0) CYCLE
  !IF(cuCntl%lGcCMFD) THEN
  !  CALL cuCopyFlux(cuCMFD_adj, cuDevice, 2)
  !END IF
END DO
CALL cuCopyFlux(cuCMFD_adj, cuDevice, 2)
ALLOCATE(CmInfo%phic_adj(cuGeometry%nxy, myzb:myze, ng))
! !$OMP PARALLEL PRIVATE(ixy_map, ipin, phisum)
! !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
 DO iz = myzb, myze
   DO ixy = 1, nxy
     phisum = 0.0
     ixy_map = pinMap(ixy)
     DO izf = fmRange(iz, 1), fmRange(iz, 2)
       DO ig =1 ,ng
       phisum(ig) = phisum(ig) + cuCMFD_adj%h_phis8(ig, ixy, izf) * (hzfm(izf) / hz(iz))
       END DO
     ENDDO
     DO j = 1, superPin(ixy_map)%nxy
       ipin = superPin(ixy_map)%pin(j)
       DO ig = 1, ng
        CmInfo%phic_adj(ipin, iz, ig) = phisum(ig)
       END DO
     ENDDO
   ENDDO
 ENDDO
! !$OMP END DO
! !$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE CUDATransient_Adpt_Driver1()
USE PARAM
USE geom,               ONLY : Core
USE RAYS,               ONLY : RayInfo
USE Core_mod,           ONLY : FmInfo,            CmInfo,             THInfo,           GroupInfo,   &
                               eigv
USE PE_Mod,             ONLY : PE
USE files,              ONLY : io8
USE ioutil,             ONLY : message,           terminate,          ShowHbar2
USE CNTL,               ONLY : nTracerCntl
USE itrcntl_mod,        ONLY : ItrCntl
USE TRAN_MOD,           ONLY : TranInfo,          TranCntl,           SavePrevPKEParameters,  RecoverPrevPKEParameters, &
                               showtransientflag, CheckCondMOC,       CheckSGFSP,             CheckMocUpdt
USE BasicOperation,     ONLY : CP_CA,             CP_VA,            MULTI_CA
USE XsPerturb_mod,      ONLY : XsPerturbation,    InitXsPerturbation, SaveXsPerturbation, RecoverXsPerturbation
USE CUDA_INIT,          ONLY : AllocTransientVar
USE CUDA_Transient
USE MOC_COMMON,         ONLY : SetCoreMacXS,      SetCoreMacXS_Cusping
USE CMFD_COMMON,        ONLY : HomogenizeXS_Cusping,  SetRadialCoupling, HomogenizeKinParam,  SetCMFDPrecCoeff
USE CUDA_CMFD,          ONLY : cuCMFDPsiUpdt,     cuSetCoarsePhis,      cuCopyFlux
USE SUBGRP_MOD,         ONLY : FxrChiGen, SubGrpEffXsGen
USE TH_Mod,             ONLY : TransientTH,       THFeedback,         SaveTranTHSol,    RecoverTranTHSol ,              &
                               SaveBaseTH,        RecoverBaseTH,      UpdtTdopErr,      SaveFxrTemp,  RecoverFxrTemp
#ifdef MPI_ENV
USE MPIComm_mod,      ONLY : MPI_SYNC
#endif
IMPLICIT NONE

REAL, POINTER :: PKEparamSave(:)
REAL, POINTER :: TDopSave(:,:)
REAL :: ampsave(100)
REAL :: tcmfd, tth, tcmfd_prev, tth_prev
REAL :: dtpke, dtcmfd, dtth, dtmoc
REAL :: powd, rhod, betad, gammad, alpha, alphad
REAL :: pketherr
REAL :: inv_eigv0, ampd, ampd_cmfd, fmult, powd_cmfd
REAL :: amperr, TdopErr, shapechg
INTEGER :: nxy, myzb, myze, ng, nprec, nparam, nzCMFD, nthstep, nfsr
INTEGER :: istep, ithstep, n, ierr
CHARACTER*20 :: card
LOGICAL :: l3dim, lfinish, lxslib, lthfinish
LOGICAL :: pkethconv, lthconv, lcmfdconv, lampconv, lstepconv

nfsr = cuGeometry%nfsr
nxy = cuGeometry%nxyc
myzb = cuDevice%myzb
myze = cuDevice%myze
nzCMFD = cuDevice%nzCMFD
ng = cuGeometry%ng
nprec = cuGeometry%nprec
l3dim = nTracerCntl%l3dim
lxslib = nTracerCntl%lxslib
IF(nTracerCntl%lKineticBen) THEN
  nTracerCntl%lchidgen = .FALSE.
  nTracerCntl%lchidkgen = .FALSE.
END IF
TranCntl%lchidk = nTracerCntl%lchidkgen

nparam = 4 + 2 * nprec
ALLOCATE(PKEparamSave(nparam))
ALLOCATE(TDopSave(Core%nz, Core%nxy))
ALLOCATE(ThInfo%TDopBase(Core%nz, Core%nxy))

WRITE(mesg, '(a)') 'Initialize Transient Calculation...'
IF (PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)

tcmfd = 0
tth = 0
istep = 0
lfinish = .FALSE.

dtpke = TranCntl%Tdiv_Corrector
dtcmfd = TranCntl%Tdiv_inp(1)
dtth = dtcmfd * 0.5
dtmoc = dtcmfd

!-Initialize Transient Calculation
l3dim = nTracerCntl%l3dim
CALL AllocTransientVar(cuDevice, PE%lCUDACMFD, l3dim, TranCntl%lchidk)
CALL cuInitTransient(Core, FmInfo, CmInfo, RayInfo, ThInfo, GroupInfo, nTracerCntl, TranInfo, TranCntl, eigv, l3dim)
inv_eigv0 = 1./ TranInfo%eigv0
CALL HomogenizeKinParam(Core, FmInfo, TranInfo, TranCntl, GroupInfo, cuGeometry%superPin, cuCMFD%PinXS,&
                        nxy, myzb, myze, ng, nprec, lxslib)
CALL cuTranReactivityUpdt_Adj(GroupInfo, TranInfo, TranCntl, l3dim, .TRUE.)
IF(TranCntl%lCorrector) CALL cuSavePKEParameters(TranInfo, TranCntl)
CALL GetShape(TranInfo)
CALL SaveShape()
TDopSave(1:Core%nz, 1:Core%nxy) = THInfo%Tdop(1:Core%nz, 1:Core%nxy)
ThInfo%TDopBase(1:Core%nz, 1:Core%nxy) = ThInfo%Tdop(1:Core%nz, 1:Core%nxy)
CALL CheckMOCUpdt(Core, FmInfo, CmInfo, GroupInfo, nTracerCntl, TranCntl, PE, .TRUE.)
IF(PE%MASTER) THEN
  IF(.NOT. nTracerCntl%lAdjoint) THEN
    WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, 1e11.4, x, a, 1e10.3, x, a, 1e12.5)') '@1', 0.,"s, P=", TranInfo%PowerLevel * 100, ', $= ', TranInfo%reactivity_Dynamic, ', G=', TranInfo%lifetime_Dynamic
    CALL message(io8, TRUE, TRUE, MESG)
    WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, F7.1, x, a, F8.5, x, a, 1e10.3)') '@2', 0.,"s, r(pcm)= ", TranInfo%delrho_Dynamic/10._8, ', k= ', TranInfo%TranEig_Dynamic/10._8, ', B= ', TranInfo%CoreBeta_Dynamic
    CALL message(io8, TRUE, TRUE, MESG)
  ELSE
    WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, 1e11.4, x, a, 1e10.3, x, a, 1e12.5)') '@3', 0.,"s, P=", TranInfo%PowerLevel * 100, ', $= ', TranInfo%reactivity_Dynamic, ', G=', TranInfo%lifetime_dynamic
    CALL message(io8, TRUE, TRUE, MESG)
    WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, F7.1, x, a, F8.5, x, a, 1e10.3)') '@4', 0.,"s, r(pcm)= ", &
      TranInfo%delrho_Dynamic/10._8, ', k= ', TranInfo%TranEig_Dynamic/10._8, ', B= ', TranInfo%CoreBeta_Dynamic
    CALL message(io8, TRUE, TRUE, MESG)
  ENDIF
  WRITE(mesg, '(A)') hbar2(1:77)
  CALL message(io8, FALSE, TRUE, MESG)
ENDIF
CALL cuUpdtResSrc(Core, FmInfo, cuCMFD, cuDevice, TranInfo, TranCntl, l3dim)
CALL InitXsPerturbation(Core, TranCntl, PE)

DO WHILE(.TRUE.)
  CALL InitIterVar(Core, FMInfo, CmInfo, GroupInfo,  .FALSE., ItrCntl, nTracerCntl, PE)
  ItrCntl%innerit = 0
  istep = istep + 1
  tcmfd_prev = tcmfd
  tcmfd = tcmfd + dtcmfd
  IF(tcmfd .GT. TranCntl%Tend .OR. abs(TranCntl%Tend - tcmfd) .LT. 1.E-6) THEN
    tcmfd = TranCntl%Tend
    dtcmfd = tcmfd - tcmfd_prev
    lfinish = .TRUE.
  END IF
  TranCntl%Nowstep = istep
  TranCntl%T(istep) = tcmfd
  TranCntl%DelT(istep) = dtcmfd
  ItrCntl%lConv = .FALSE.

  lcmfdconv = .FALSE.
  lthconv = .FALSE.
  lampconv = .FALSE.
  lstepconv = .FALSE.
  !-Save Information at base time
  CALL SavePrevPKEParameters(TranInfo, PKEparamSave, nparam)
  ampd_cmfd = TranInfo%amp
  powd_cmfd = TranInfo%PowerLevel

  tth = tcmfd_prev
  lthfinish = .FALSE.
  ithstep = 0
  DO WHILE(.TRUE.) !-PKE & TH
    ithstep = ithstep + 1
    tth_prev = tth
    tth = tth + dtth
    IF(tth .GT. tcmfd .OR. abs(tcmfd - tth) .LT. 1.E-6) THEN
      tth = tcmfd
      dtth = tth - tth_prev
      lthfinish = .TRUE.
    END IF
    CALL XsPerturbation(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, PE, tth, tth_prev)
    IF(TranCntl%lXSPerturb) THEN
      IF(PE%Master) CALL ShowTransientFlag(io8, 0, tth_prev, tth)
      CALL THFeedBack(Core, CmInfo, FmInfo, ThInfo, nTRACERCntl, PE)
#ifdef MPI_ENV
      CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif
      CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, Eigv, GroupInfo, nTracerCntl, PE)
      CALL KinParamGen(Core, FmInfo, TranInfo, ThInfo, GroupInfo, .TRUE., nTracerCntl, PE)
      IF (nTracerCntl%lMacro) CALL SetCoreMacXs(Core, FmInfo)
      WRITE(mesg, '(a)') 'Cell Homogenization...'
      IF (PE%Master) CALL message(io8, .TRUE., .TRUE., mesg)
      CALL FxrChiGen(Core, FmInfo%Fxr, FmInfo, GroupInfo, PE, nTracerCntl)
      CALL HomogenizeXS_Cusping(Core, FmInfo, cuGeometry%superPin, FmInfo%Fxr, cuCMFD%PinXS, FmInfo%phis, cuGeometry%ng, cuGeometry%nxyc, &
        cuDevice%myzb, cuDevice%myze, PE%nCMFDThread, nTracerCntl%lxslib, nTracerCntl%lscat1, .FALSE.)
      CALL HomogenizeKinParam(Core, FmInfo, TranInfo, TranCntl, GroupInfo, cuGeometry%superPin, cuCMFD%PinXS,&
        nxy, myzb, myze, ng, nprec, lxslib)
      CALL cuSetCMFDSourceOperator(cuCMFD, cuDevice, GroupInfo%lUpscat)
      CALL cuCMFDPsiUpdt(cuCMFD, cuDevice)
      n = cuGeometry%nxyc * cuDevice%nzCMFD
      ierr = cublasDscal_v2(cuDevice%myblasHandle, n, inv_eigv0, cuCMFD%psi8, 1)
      CALL cuTranReactivityUpdt_Adj(GroupInfo, TranInfo, TranCntl, l3dim, .FALSE.)
    END IF
    !-1.Power Level Prediction with PKE solver
    ampd = TranInfo%amp
    powd = TranInfo%PowerLevel
    DO WHILE(.TRUE.)
      IF(PE%Master) CALL ShowTransientFlag(io8, 1, tth_prev, tth)
      CALL Transient_Corrector(TranInfo, TranCntl, nTracerCntl, tth_prev, tth, dtpke)
      fmult = TranInfo%amp / ampd_cmfd
      WRITE(mesg, '(6x,3(a11,es11.4),a)'), 'Begin Amp: ', ampd,  'End Amp: ', TranInfo%amp, 'Ratio: ', fmult*100., ' [%]'
      IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
      IF(TranCntl%lXSPerturb) THEN
        CALL cuUpdtPowerLevel(TranInfo, .FALSE.)
      ELSE
        TranInfo%PowerLevel = powd
      END IF
      TranInfo%PowerLevel = fmult * TranInfo%PowerLevel
      nTracerCntl%PowerLevel = TranInfo%PowerLevel; ThInfo%PowLv = nTracerCntl%PowerLevel
      WRITE(mesg, '(6x,3(a11,es11.4),a)'), 'Begin Pow: ', powd,  'End Pow: ', TranInfo%PowerLevel, 'Ratio: ', TranInfo%PowerLevel/powd*100., ' [%]'
      IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
      !-2.TH Update
      pkethconv = .FALSE.
      IF(PE%Master) CALL ShowTransientFlag(io8, 2, tth_prev, tth)
      CALL TransientTH(Core, CmInfo, FmInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE, dtth)
      IF(ThInfo%TDopChg .LT. 5.E-4) THEN
        pkethconv = .TRUE.
      ELSE
          CALL THFeedBack(Core, CmInfo, FmInfo, ThInfo, nTRACERCntl, PE)
#ifdef MPI_ENV
          CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif
        IF (nTracerCntl%lrestrmt) THEN
          CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, Eigv, GroupInfo, nTracerCntl, PE)
          CALL KinParamGen(Core, FmInfo, TranInfo, ThInfo, GroupInfo, .TRUE., nTracerCntl, PE)
          IF (nTracerCntl%lMacro) CALL SetCoreMacXs(Core, FmInfo)
        ENDIF
        WRITE(mesg, '(a)') 'Cell Homogenization...'
        IF (PE%Master) CALL message(io8, .TRUE., .TRUE., mesg)
        CALL FxrChiGen(Core, FmInfo%Fxr, FmInfo, GroupInfo, PE, nTracerCntl)
        CALL HomogenizeXS_Cusping(Core, FmInfo, cuGeometry%superPin, FmInfo%Fxr, cuCMFD%PinXS, FmInfo%phis, cuGeometry%ng, cuGeometry%nxyc, &
          cuDevice%myzb, cuDevice%myze, PE%nCMFDThread, nTracerCntl%lxslib, nTracerCntl%lscat1, .FALSE.)
        CALL HomogenizeKinParam(Core, FmInfo, TranInfo, TranCntl, GroupInfo, cuGeometry%superPin, cuCMFD%PinXS,&
          nxy, myzb, myze, ng, nprec, lxslib)
        CALL cuSetCMFDSourceOperator(cuCMFD, cuDevice, GroupInfo%lUpscat)
        CALL cuCMFDPsiUpdt(cuCMFD, cuDevice)
        n = cuGeometry%nxyc * cuDevice%nzCMFD
        ierr = cublasDscal_v2(cuDevice%myblasHandle, n, inv_eigv0, cuCMFD%psi8, 1)

        !--Save PKE Parameters
        rhod = TranInfo%delrho_dynamic;  betad = TranInfo%corebeta_dynamic; gammad = TranInfo%lifetime_dynamic
        CALL cuTranReactivityUpdt_Adj(GroupInfo, TranInfo, TranCntl, l3dim, .FALSE.)
        alpha = (TranInfo%delrho_dynamic * 1.E-5 - TranInfo%corebeta_dynamic) / TranInfo%lifetime_dynamic
        alphad = (rhod * 1.E-5 - betad) / gammad
        pketherr = (alpha - alphad) * dtpke
        IF(abs(pketherr) .LT. 1.E-5) pkethconv = .TRUE.
        WRITE(mesg, '(4x,a27,x,es12.4,a15,L3)'), 'ESTIMATED PKE LTE from TH: ', pketherr, ' Convergence - ', pkethconv
        IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
      END IF
      IF (pkethconv) THEN
        IF(ithstep .EQ. 1 .AND. tcmfd .NE. tth) THEN
          CALL SaveXsPerturbation(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, PE, tcmfd, tth)
        END IF
        IF(ithstep .NE. 1 .OR. tcmfd .NE. tth) THEN
          CALL SaveFxrTemp(FmInfo%Fxr, ithstep, Core%nCoreFxr, PE%myzb, PE%myze)
        END IF
        IF(ithstep .EQ. 1) THEN
          CALL SaveBaseTH(Core, THInfo)
        END IF
        ampsave(ithstep) = TranInfo%amp
        CALL cuSavePKEParameters(TranInfo, TranCntl, .TRUE.)
        CALL SaveTranTHsol(Core, ThInfo, TranCntl, nTracerCntl, PE)
        EXIT
      ELSE
        !--Recover amplitude
        TranInfo%amp = ampd
      END IF
    END DO
    IF(lthfinish) THEN
      CALL UpdtTdopErr(ThInfo, TDopErr, TDopSave, Core%nxy, Core%nz)
      IF(TDopErr .LT. 5.E-4) lthconv = .TRUE.
      nthstep = ithstep
      EXIT
    END IF
  END DO
  !-Predictor CMFD calculation
  CALL CheckSGFSP(Core, ThInfo, TranCntl, PE)
  IF(TranCntl%lSGFSPUpdt) THEN
    CALL CUDASubGrpFsp(Core, FmInfo, THInfo, RayInfo, GroupInfo)
#ifdef MPI_ENV
    CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif
    CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, Eigv, GroupInfo, nTracerCntl, PE)
    CALL KinParamGen(Core, FmInfo, TranInfo, ThInfo, GroupInfo, .TRUE., nTracerCntl, PE)
    IF (nTracerCntl%lMacro) CALL SetCoreMacXs(Core, FmInfo)
  END IF

  IF(PE%Master) CALL ShowTransientFlag(io8, 3, tcmfd_prev, tcmfd)
  !--Reflect PKE amplitude
  CALL ScaleFmFlux(FmInfo, TranInfo)
  IF(nthstep .GT. 1) TranInfo%AmpRatio = TranInfo%Amp / ampd_cmfd
  cuTranCMInfo%Expo = TranInfo%AmpRatio
  cuTranCmInfo%Expo_Alpha = TranInfo%AmpTilt
  IF(cuCntl%lAxial .AND. l3dim) THEN
    cuAxial%Expo = TranInfo%AmpRatio
    cuAxial%Expo_Alpha = TranInfo%AmpTilt
  END IF
  TranInfo%fmExpo = TranInfo%AmpRatio
  TranInfo%fmExpo_Alpha = TranInfo%AmpTilt
  CALL CUDATranCMFDAcc(Core, CmInfo, FmInfo, TranInfo, TranCntl, .TRUE.)
  !-Determine RayTracing Options
  CALL CheckMOCUpdt(Core, FmInfo, CmInfo, GroupInfo, nTracerCntl, TranCntl, PE, .FALSE.)
  !--Save PKE Parameters
  IF(TranCntl%lMOCUpdt) THEN
    CALL CUDATranMOCSweep(Core, RayInfo, FmInfo, TranInfo, TranCntl)
    CALL CUDATranCMFDAcc(Core, CmInfo, FmInfo, TranInfo, TranCntl, .FALSE.)
  END IF
  rhod = TranInfo%delrho_dynamic;  betad = TranInfo%corebeta_dynamic; gammad = TranInfo%lifetime_dynamic
  CALL cuTranReactivityUpdt_Adj(GroupInfo, TranInfo, TranCntl, l3dim, .FALSE.)
  CALL GetShape(TranInfo)
  shapechg = ShapeChange()
  ampsave(1:nthstep) = ampsave(1:nthstep) * amp_now / ampsave(nthstep)
  WRITE(mesg, '(2(a15,x, es12.4))'), 'SHAPE CHANGE: ', shapechg, 'TDOP ERROR: ', TDopErr

  IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
  WRITE(mesg, '(a15,x, 3a12)'), 'PKE Parameters', 'RHO [pcm]',  'BETA [pcm]', 'GAMMA [pcm]'
  IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
  WRITE(mesg, '(a15,x, 4es12.4)'), 'CMFD UPDATED: ', TranInfo%delrho_dynamic,  TranInfo%corebeta_dynamic*1.E5, TranInfo%lifetime_dynamic*1.E5
  IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
  WRITE(mesg, '(a15,x, 4es12.4)'), 'DIFFERENCE: ', (TranInfo%delrho_dynamic-rhod),  (TranInfo%corebeta_dynamic-betad)*1.E5, (TranInfo%lifetime_dynamic-gammad)*1.E5
  IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)


  DO WHILE(.TRUE.) !- CMFD Iteration
    IF(lampconv .AND. lthconv) GO TO 300
    !-Recover information at base time
    IF(.NOT. lampconv) THEN
      CALL RecoverPrevPKEParameters(TranInfo, PKEparamSave, nparam)
      TranInfo%amp = ampd_cmfd
    END IF
    CALL RecoverBaseTH(Core, THInfo)
    TranInfo%PowerLevel = powd_cmfd
    nTracerCntl%PowerLevel = TranInfo%PowerLevel; ThInfo%PowLv = nTracerCntl%PowerLevel

    tth = tcmfd_prev
    lthfinish = .FALSE.
    ithstep = 0
    DO WHILE(.TRUE.) !-PKE & TH Proceed
      ithstep = ithstep + 1
      tth_prev = tth
      tth = tth + dtth
      IF(tth .GT. tcmfd .OR. abs(tcmfd - tth) .LT. 1.E-6) THEN
        tth = tcmfd
        dtth = tth - tth_prev
        lthfinish = .TRUE.
      END IF
      IF(ithstep .EQ. 1) THEN
        IF(tcmfd .NE. tth) THEN
          CALL RecoverXsPerturbation(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, PE, tcmfd, tth)
          IF(TranCntl%lXSPerturb) THEN
            IF(PE%Master) CALL ShowTransientFlag(io8, 0, tth_prev, tth)
          END IF
        END IF
      ELSE
        CALL XsPerturbation(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, PE, tth, tth_prev)
        IF(TranCntl%lXSPerturb) THEN
          IF(PE%Master) CALL ShowTransientFlag(io8, 0, tth_prev, tth)
        END IF
      END IF

      ampd = TranInfo%amp
      powd = TranInfo%PowerLevel
      !- Generate Flux for TH time-step
      CALL InterpolateShape(tth, tcmfd_prev, tcmfd, ampsave(ithstep))
      CALL cuCopyFlux(cuCMFD, cuDevice, 2)
      CALL cuSetCoarsePhis(Core, cuCMFD, cuDevice, CmInfo%Phic)
      !!-TH solution with new shape
      !CALL cuUpdtPowerLevel(TranInfo, .FALSE.)
      !nTracerCntl%PowerLevel = TranInfo%PowerLevel; ThInfo%PowLv = nTracerCntl%PowerLevel
      !
      !!IF(lthconv .AND. .NOT. TranCntl%lXsPerturb) THEN
      !!  IF(nthstep .GT. 1) THEN
      !!    CALL RecoverFxrTemp(FmInfo%Fxr, ithstep, Core%nCoreFxr, PE%myzb, PE%myze)
      !!  END IF
      !!ELSE
      !!  IF(PE%Master) CALL ShowTransientFlag(io8, 2, tth_prev, tth)
      !!  CALL TransientTH(Core, CmInfo, FmInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE, dtth)
      !!  CALL THFeedBack(Core, CmInfo, FmInfo, ThInfo, nTRACERCntl, PE)
      !!END IF
      !IF(PE%Master) CALL ShowTransientFlag(io8, 2, tth_prev, tth)
      !CALL TransientTH(Core, CmInfo, FmInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE, dtth)
      !CALL THFeedBack(Core, CmInfo, FmInfo, ThInfo, nTRACERCntl, PE)
      !
      !IF(.NOT. lthconv .OR. nthstep .GT. 1) THEN
      !  CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, Eigv, GroupInfo, nTracerCntl, PE)
      !  CALL KinParamGen(Core, FmInfo, TranInfo, ThInfo, GroupInfo, .TRUE., nTracerCntl, PE)
      !END IF
      !IF (nTracerCntl%lMacro) CALL SetCoreMacXs(Core, FmInfo)
      !WRITE(mesg, '(a)') 'Cell Homogenization...'
      !IF (PE%Master) CALL message(io8, .TRUE., .TRUE., mesg)
      !!CALL FxrChiGen(Core, FmInfo%Fxr, FmInfo, GroupInfo, PE, nTracerCntl)
      !CALL HomogenizeXS_Cusping(Core, FmInfo, cuGeometry%superPin, FmInfo%Fxr, cuCMFD%PinXS, FmInfo%phis, cuGeometry%ng, cuGeometry%nxyc, &
      !  cuDevice%myzb, cuDevice%myze, PE%nCMFDThread, nTracerCntl%lxslib, nTracerCntl%lscat1, .FALSE.)
      !CALL HomogenizeKinParam(Core, FmInfo, TranInfo, TranCntl, GroupInfo, cuGeometry%superPin, cuCMFD%PinXS,&
      !  nxy, myzb, myze, ng, nprec, lxslib)
      !CALL cuSetCMFDSourceOperator(cuCMFD, cuDevice, GroupInfo%lUpscat)
      !CALL cuCMFDPsiUpdt(cuCMFD, cuDevice)
      !n = cuGeometry%nxyc * cuDevice%nzCMFD
      !ierr = cublasDscal_v2(cuDevice%myblasHandle, n, inv_eigv0, cuCMFD%psi8, 1)
      !!CALL cuTranReactivityUpdt_Adj(GroupInfo, TranInfo, TranCntl, l3dim, .FALSE.)
      !
      !!- �ӽ�
      !rhod = TranInfo%delrho_dynamic;  betad = TranInfo%corebeta_dynamic; gammad = TranInfo%lifetime_dynamic
      !CALL cuTranReactivityUpdt_Adj(GroupInfo, TranInfo, TranCntl, l3dim, .FALSE.)
      !alpha = (TranInfo%delrho_dynamic * 1.E-5 - TranInfo%corebeta_dynamic) / TranInfo%lifetime_dynamic
      !alphad = (rhod * 1.E-5 - betad) / gammad
      !pketherr = (alpha - alphad) * dtpke
      !WRITE(mesg, '(4x,a27,x,es12.4,a15,L3)'), 'ESTIMATED PKE LTE from TH: ', pketherr, ' Convergence - ', .FALSE.
      !IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)

      IF(TranCntl%lXSPerturb .OR. nthstep .GT. 1) THEN
        CALL cuUpdtPowerLevel(TranInfo, .FALSE.)
        nTracerCntl%PowerLevel = TranInfo%PowerLevel; ThInfo%PowLv = nTracerCntl%PowerLevel
        CALL TransientTH(Core, CmInfo, FmInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE, dtth)
        IF(PE%Master) CALL ShowTransientFlag(io8, 0, tth_prev, tth)
        CALL THFeedBack(Core, CmInfo, FmInfo, ThInfo, nTRACERCntl, PE)
#ifdef MPI_ENV
        CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif
        CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, Eigv, GroupInfo, nTracerCntl, PE)
        CALL KinParamGen(Core, FmInfo, TranInfo, ThInfo, GroupInfo, .TRUE., nTracerCntl, PE)
        IF (nTracerCntl%lMacro) CALL SetCoreMacXs(Core, FmInfo)
        WRITE(mesg, '(a)') 'Cell Homogenization...'
        IF (PE%Master) CALL message(io8, .TRUE., .TRUE., mesg)
        CALL FxrChiGen(Core, FmInfo%Fxr, FmInfo, GroupInfo, PE, nTracerCntl)
        CALL HomogenizeXS_Cusping(Core, FmInfo, cuGeometry%superPin, FmInfo%Fxr, cuCMFD%PinXS, FmInfo%phis, cuGeometry%ng, cuGeometry%nxyc, &
          cuDevice%myzb, cuDevice%myze, PE%nCMFDThread, nTracerCntl%lxslib, nTracerCntl%lscat1, .FALSE.)
        CALL HomogenizeKinParam(Core, FmInfo, TranInfo, TranCntl, GroupInfo, cuGeometry%superPin, cuCMFD%PinXS,&
          nxy, myzb, myze, ng, nprec, lxslib)
        CALL cuSetCMFDSourceOperator(cuCMFD, cuDevice, GroupInfo%lUpscat)
        CALL cuCMFDPsiUpdt(cuCMFD, cuDevice)
        n = cuGeometry%nxyc * cuDevice%nzCMFD
        ierr = cublasDscal_v2(cuDevice%myblasHandle, n, inv_eigv0, cuCMFD%psi8, 1)
        CALL cuTranReactivityUpdt_Adj(GroupInfo, TranInfo, TranCntl, l3dim, .FALSE.)
      END IF

      DO WHILE(.TRUE.) ! PKE&TH Iteration
        IF(lampconv) EXIT
        IF(PE%Master) CALL ShowTransientFlag(io8, 1, tth_prev, tth)
        CALL Transient_Corrector(TranInfo, TranCntl, nTracerCntl, tth_prev, tth, dtpke)
        fmult = TranInfo%amp / ampsave(ithstep)
        WRITE(mesg, '(6x,3(a11,es11.4),a)'), 'Begin Amp: ', ampd,  'End Amp: ', TranInfo%amp, 'Ratio: ', TranInfo%amp/ampd*100., ' [%]'
        IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
        CALL cuUpdtPowerLevel(TranInfo, .FALSE.)
        TranInfo%PowerLevel = fmult * TranInfo%PowerLevel
        nTracerCntl%PowerLevel = TranInfo%PowerLevel; ThInfo%PowLv = nTracerCntl%PowerLevel
        WRITE(mesg, '(6x,3(a11,es11.4),a)'), 'Begin Pow: ', powd,  'End Pow: ', TranInfo%PowerLevel, 'Ratio: ', TranInfo%PowerLevel/powd*100., ' [%]'
        IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
        !-2.TH Update
        pkethconv = .FALSE.
        IF(PE%Master) CALL ShowTransientFlag(io8, 2, tth_prev, tth)
        CALL TransientTH(Core, CmInfo, FmInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE, dtth)
        IF(ThInfo%TDopChg .LT. 5.E-4) THEN
          pkethconv = .TRUE.
        ELSE
          CALL THFeedBack(Core, CmInfo, FmInfo, ThInfo, nTRACERCntl, PE)
#ifdef MPI_ENV
          CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif
          IF (nTracerCntl%lrestrmt) THEN
            CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, Eigv, GroupInfo, nTracerCntl, PE)
            CALL KinParamGen(Core, FmInfo, TranInfo, ThInfo, GroupInfo, .TRUE., nTracerCntl, PE)
          ENDIF
          IF (nTracerCntl%lMacro) CALL SetCoreMacXs(Core, FmInfo)
          WRITE(mesg, '(a)') 'Cell Homogenization...'
          IF (PE%Master) CALL message(io8, .TRUE., .TRUE., mesg)
          CALL FxrChiGen(Core, FmInfo%Fxr, FmInfo, GroupInfo, PE, nTracerCntl)
          CALL HomogenizeXS_Cusping(Core, FmInfo, cuGeometry%superPin, FmInfo%Fxr, cuCMFD%PinXS, FmInfo%phis, cuGeometry%ng, cuGeometry%nxyc, &
            cuDevice%myzb, cuDevice%myze, PE%nCMFDThread, nTracerCntl%lxslib, nTracerCntl%lscat1, .FALSE.)
          CALL HomogenizeKinParam(Core, FmInfo, TranInfo, TranCntl, GroupInfo, cuGeometry%superPin, cuCMFD%PinXS,&
            nxy, myzb, myze, ng, nprec, lxslib)
          CALL cuSetCMFDSourceOperator(cuCMFD, cuDevice, GroupInfo%lUpscat)
          CALL cuCMFDPsiUpdt(cuCMFD, cuDevice)
          n = cuGeometry%nxyc * cuDevice%nzCMFD
          ierr = cublasDscal_v2(cuDevice%myblasHandle, n, inv_eigv0, cuCMFD%psi8, 1)

          !--Save PKE Parameters
          rhod = TranInfo%delrho_dynamic;  betad = TranInfo%corebeta_dynamic; gammad = TranInfo%lifetime_dynamic
          CALL cuTranReactivityUpdt_Adj(GroupInfo, TranInfo, TranCntl, l3dim, .FALSE.)
          alpha = (TranInfo%delrho_dynamic * 1.E-5 - TranInfo%corebeta_dynamic) / TranInfo%lifetime_dynamic
          alphad = (rhod * 1.E-5 - betad) / gammad
          pketherr = (alpha - alphad) * dtpke
          IF(abs(pketherr) .LT. 1.E-5) pkethconv = .TRUE.
          WRITE(mesg, '(1x,a27,x,es12.4,a15,L3)'), 'ESTIMATED PKE LTE from TH: ', pketherr, ' Convergence - ', pkethconv
          IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
        END IF
        IF (pkethconv) THEN
          ampsave(ithstep) = TranInfo%amp
          CALL cuSavePKEParameters(TranInfo, TranCntl)
          EXIT
        ELSE
          !--Recover amplitude
          TranInfo%amp = ampd
        END IF
      END DO

      IF(ithstep .EQ. 1) THEN
        CALL SaveBaseTH(Core, THInfo)
      END IF
      CALL SaveTranTHsol(Core, ThInfo, TranCntl, nTracerCntl, PE)

      IF(lthfinish) THEN
        amperr = (amp_now - TranInfo%amp) / TranInfo%amp
        IF(abs(amperr) .LT. 1.E-5) lampconv = .TRUE.
        CALL UpdtTdopErr(ThInfo, TDopErr, TDopSave, Core%nxy, Core%nz)
        IF(TDopErr .GT. 0.01) THEN
          CALL CUDASubGrpFsp(Core, FmInfo, THInfo, RayInfo, GroupInfo)
#ifdef MPI_ENV
          CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif
          CALL SubGrpEffXsGen(Core, FmInfo%Fxr, THInfo, Eigv, GroupInfo, nTracerCntl, PE)
          CALL KinParamGen(Core, FmInfo, TranInfo, ThInfo, GroupInfo, .TRUE., nTracerCntl, PE)
          IF (nTracerCntl%lMacro) CALL SetCoreMacXs(Core, FmInfo)
        END IF
        IF(TDopErr .LT. 5.E-4) lthconv = .TRUE.
        EXIT
      END IF
    END DO
300 CONTINUE

    IF(TranCntl%cmfdres .LT. 1.E-6) lcmfdconv = .TRUE.
    lstepconv = lcmfdconv .AND. lthconv .AND. ItrCntl%lConv
    IF((ItrCntl%Cmfdit - ItrCntl%Cmfdit0) .LT. 5) lstepconv = .FALSE.
    WRITE(mesg, '(1x,a27,x,es12.4,a15,L3)'), 'AMPLITUDE CHANGE: ', amperr, ' CONVERGENCE - ', lampconv
    IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
    WRITE(mesg, '(1x,a27,x,es12.4,a15,L3)'), 'DOPPLER CHANGE  : ', TdopErr , ' CONVERGENCE - ', lthconv
    IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
    WRITE(mesg, '(1x,a27,x,es12.4,a15,L3)'), 'CMFD RESIDUAL   : ', TranCntl%cmfdres , ' CONVERGENCE - ', lcmfdconv
    IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
    IF(TranCntl%lMOCUpdt) THEN
      WRITE(mesg, '(1x,a27,x,es12.4,a15,L3)'), 'FM PSI ERROR    : ', TranCntl%delpsifm , ' CONVERGENCE - ', ItrCntl%lConv
      IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
    END IF
    IF(lstepconv) THEN
      CALL CorrectAmp(Core, CmInfo, FmInfo, GroupInfo, TranInfo, nTracerCntl)
      EXIT
    END IF
    !--Reflect PKE amplitude
    CALL ScaleFmFlux(FmInfo, TranInfo)
    IF(nthstep .GT. 1) TranInfo%AmpRatio = TranInfo%Amp / ampd_cmfd
    cuTranCMInfo%Expo = TranInfo%AmpRatio
    cuTranCmInfo%Expo_Alpha = TranInfo%AmpTilt
    IF(cuCntl%lAxial .AND. l3dim) THEN
      cuAxial%Expo = TranInfo%AmpRatio
      cuAxial%Expo_Alpha = TranInfo%AmpTilt
    END IF
    TranInfo%fmExpo = TranInfo%AmpRatio
    TranInfo%fmExpo_Alpha = TranInfo%AmpTilt
    IF(TranCntl%lMOCUpdt) THEN
      CALL CUDATranMOCSweep(Core, RayInfo, FmInfo, TranInfo, TranCntl)
    END IF
    IF(PE%Master) CALL ShowTransientFlag(io8, 3, tcmfd_prev, tcmfd)
    CALL CUDATranCMFDAcc(Core, CmInfo, FmInfo, TranInfo, TranCntl, .FALSE.)
    !--Save PKE Parameters
    rhod = TranInfo%delrho_dynamic;  betad = TranInfo%corebeta_dynamic; gammad = TranInfo%lifetime_dynamic
    CALL cuTranReactivityUpdt_Adj(GroupInfo, TranInfo, TranCntl, l3dim, .FALSE.)
    WRITE(mesg, '(a15,x, 3a12)'), 'PKE Parameters', 'RHO [pcm]',  'BETA [pcm]', 'GAMMA [pcm]'
    IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
    WRITE(mesg, '(a15,x, 4es12.4)'), 'CMFD UPDATED: ', TranInfo%delrho_dynamic,  TranInfo%corebeta_dynamic*1.E5, TranInfo%lifetime_dynamic*1.E5
    IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
    WRITE(mesg, '(a15,x, 4es12.4)'), 'DIFFERENCE: ', (TranInfo%delrho_dynamic-rhod),  (TranInfo%corebeta_dynamic-betad)*1.E5, (TranInfo%lifetime_dynamic-gammad)*1.E5
    IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
    CALL GetShape(TranInfo)
    ampsave(1:nthstep) = ampsave(1:nthstep) * amp_now / ampsave(nthstep)
  END DO

  !Theta method source term
  CALL cuUpdtResSrc(Core, FmInfo, cuCMFD, cuDevice, TranInfo, TranCntl, l3dim)
  !Precursor Update
  CALL cuUpdtPrec(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, nprec)
  !Calculate
  CALL cuUpdtPowerLevel(TranInfo, .FALSE.)
  nTracerCntl%PowerLevel = TranInfo%PowerLevel; ThInfo%PowLv = nTracerCntl%PowerLevel

  CALL cuSaveTranSol(FmInfo, nTracerCntl)
  CALL SaveTranTHsol(Core, THInfo, TranCntl, nTracerCntl, PE)
  CALL cuSavePKEParameters(TranInfo, TranCntl, .TRUE.)
  CALL SaveShape()
  IF(TranCntl%lSGFSPUpdt) ThInfo%TDopBase(1:Core%nz, 1:Core%nxy) = ThInfo%Tdop(1:Core%nz, 1:Core%nxy)
  IF(TranCntl%lMOCUpdt) CALL CheckMOCUpdt(Core, FmInfo, CmInfo, GroupInfo, nTracerCntl, TranCntl, PE, .TRUE.)

  IF(PE%MASTER) THEN
    IF(.NOT. nTracerCntl%lAdjoint) THEN
      WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, 1e11.4, x, a, 1e10.3, x, a, 1e12.5)') '@1', tcmfd,"s, P=", TranInfo%PowerLevel * 100, ', $= ', TranInfo%reactivity_Dynamic, ', G=', TranInfo%lifetime_Dynamic
      CALL message(io8, FALSE, TRUE, MESG)
      WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, F7.1, x, a, F8.5, x, a, 1e10.3)') '@2', tcmfd,"s, r(pcm)= ", TranInfo%delrho_Dynamic/10._8, ', k= ', TranInfo%TranEig_Dynamic/10._8, ', B= ', TranInfo%CoreBeta_Dynamic
      CALL message(io8, FALSE, TRUE, MESG)
    ELSE
      WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, 1e14.7, x, a, 1e14.7, x, a, 1e14.7)') '@3', tcmfd,"s, P=", TranInfo%PowerLevel * 100, ', $= ', TranInfo%reactivity_Dynamic, ', G=', TranInfo%lifetime_dynamic
      CALL message(io8, FALSE, TRUE, MESG)
      WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, F7.1, x, a, F8.5, x, a, 1e14.7)') '@4', tcmfd,"s, r(pcm)= ", &
        TranInfo%delrho_Dynamic/10._8, ', k= ', TranInfo%TranEig_Dynamic/10._8, ', B= ', TranInfo%CoreBeta_Dynamic
      CALL message(io8, FALSE, TRUE, MESG)
    ENDIF
    IF(nTracerCntl%lFeedback) THEN
      WRITE(mesg, '(a, 1p, e12.4, x, a, f10.3, x, a, f10.3, x, a)') '@5', tcmfd,"s, Max Tf=  ", ThInfo%TfMax/10, "C, Max TfAvg= ", ThInfo%Tfavg_max/10, "C"
      CALL message(io8, FALSE, TRUE, MESG)
      WRITE(mesg, '(a, 1p, e12.4, x, a, f10.3, x, a, f10.3, x, a)') '@6', tcmfd,"s, Max Tcl= ", ThInfo%Tfcl_max/10, "C, Max TMOA=  ", ThInfo%TModoutAvg/10, "C"
      CALL message(io8, FALSE, TRUE, MESG)
    ENDIF
    WRITE(mesg, '(A, 1p, e12.4, x, a, i8)') '@9', TranCntl%T(istep), "s, MOC = ", ItrCntl%Mocit - ItrCntl%Mocit0
    CALL message(io8, FALSE, TRUE, MESG)
    WRITE(mesg, '(A, 1p, e12.4, x, a, i8, x, a, i8)') '@10', TranCntl%T(istep), "s, MGCMFD = ", ItrCntl%Cmfdit - ItrCntl%Cmfdit0, "Inner = ", ItrCntl%InnerIt
    CALL message(io8, FALSE, TRUE, MESG)
    CALL ShowHbar2(io8)
  END IF
  IF(lfinish) THEN
    EXIT
  END IF
END DO

DEALLOCATE(PKEParamSave, TDopSave)

END SUBROUTINE

SUBROUTINE CUDATransient_Adpt_Driver()
USE PARAM
USE geom,               ONLY : Core
USE RAYS,               ONLY : RayInfo
USE Core_mod,           ONLY : FmInfo,            CmInfo,             THInfo,           GroupInfo,   &
                               eigv
USE PE_Mod,             ONLY : PE
USE CNTL,               ONLY : nTracerCntl
USE itrcntl_mod,        ONLY : ItrCntl
USE SUBGRP_MOD,         ONLY : SubGrpEffXsGen
USE TRAN_MOD,           ONLY : TranInfo,          TranCntl,                                         &
                               InitTransient,     KinParamGen,        InitPrecursor
USE Boron_mod,          ONLY : SetBoronCoolant
USE MOC_COMMON,         ONLY : SetCoreMacXs
USE CMFD_COMMON,        ONLY : SetRadialCoupling, HomogenizeKinParam,  SetCMFDPrecCoeff
USE XsPerturb_mod,      ONLY : XsPerturbation,    InitXsPerturbation,  SaveXsPerturbation,  RecoverXsPerturbation
USE files,              ONLY : io8
USE ioutil,             ONLY : message,           terminate
USE CUDA_MASTER
USE CUDA_INIT,          ONLY : AllocTransientVar, deallocTranCMFD
USE CUDA_Transient
USE PromptFeedback_mod, ONLY : UpdtFuelTemp, InitPromptFeedback
USE TH_Mod,             ONLY : TransientTH,       THFeedback,         SaveTranThSol
IMPLICIT NONE
REAL, POINTER :: TDopSave(:,:)
REAL :: ampSave, ampP, corrfactor
REAL :: shapechg, trunErr, shpCrit
REAL :: tprev, tnow
REAL :: tthprev, tth
REAL :: dtcmfd, dtth, dtpke
REAL :: iqsaa_r, iqsaa_h(2), iqsaa_f(2,2), iqsaa_j(2), denom, numer, iqsaa_a, iqsaa_rd
INTEGER :: nxy, myzb, myze, nprec, ng
INTEGER :: istep, i, iter, im, inewt
LOGICAL :: l3dim, lxslib, lout
LOGICAL :: lcmfdconv, lfinish, laccept

INTEGER :: ierr

nxy = cuGeometry%nxyc
myzb = cuDevice%myzb
myze = cuDevice%myze
ng = cuGeometry%ng
nprec = cuGeometry%nprec
l3dim = nTracerCntl%l3dim
lxslib = nTracerCntl%lxslib
lfinish = .FALSE.

shpCrit = 0.01

!nTracerCntl%lfeedback = .FALSE.
!nTracerCntl%lMacro = .FALSE.
IF(nTracerCntl%lKineticBen) THEN
  nTracerCntl%lchidgen = .FALSE.
  nTracerCntl%lchidkgen = .FALSE.
END IF
TranCntl%lchidk = nTracerCntl%lchidkgen

ALLOCATE(TDopSave(Core%nz, Core%nxy))

WRITE(mesg, '(a)') 'Initialize Transient Calculation...'
IF (PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)
TranCntl%NowStep = 1
TranCntl%T(1) = TranCntl%Tdiv_inp(1)
TranCntl%DelT(1) = TranCntl%Tdiv_inp(1)
CALL AllocTransientVar(cuDevice, PE%lCUDACMFD, l3dim, TranCntl%lchidk)
CALL cuInitTransient(Core, FmInfo, CmInfo, RayInfo, ThInfo, GroupInfo, nTracerCntl, TranInfo, TranCntl, eigv, l3dim)
IF(TranCntl%lDynamicBen) CALL InitPromptFeedback(Core, FmInfo, CmInfo, TranInfo, PE, ng)

CALL HomogenizeKinParam(Core, FmInfo, TranInfo, TranCntl, GroupInfo, cuGeometry%superPin, cuCMFD%PinXS,&
                        nxy, myzb, myze, ng, nprec, lxslib)
CALL SetCMFDPrecCoeff(TranInfo, TranCntl, cuCMFD%PinXS, cuTranCMInfo%CellOmegam, cuTranCMInfo%CellOmega0, cuTranCMInfo%CellOmegap, &
                      nxy, myzb, myze, nprec)
CALL cuUpdtResSrc(Core, FmInfo, cuCMFD, cuDevice, TranInfo, TranCntl, l3dim)

CALL cuTranReactivityUpdt_Adj(GroupInfo, TranInfo, TranCntl, l3dim, .TRUE.)
IF(TranCntl%lCorrector .OR. TranCntl%lGuess) CALL cuSavePKEParameters(TranInfo, TranCntl)
CALL GetShape(TranInfo)
CALL SaveShape()
trunErr = MathErr(TranInfo, TranCntl)
IF(nTracerCntl%lFeedback) TDopSave(1:Core%nz, 1:Core%nxy) = THInfo%Tdop(1:Core%nz, 1:Core%nxy)
ampSave = 1.
IF(PE%MASTER) THEN
  IF(.NOT. nTracerCntl%lAdjoint) THEN
    WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, 1e11.4, x, a, 1e10.3, x, a, 1e12.5)') '@1', 0.,"s, P=", TranInfo%PowerLevel * 100, ', $= ', TranInfo%reactivity_Dynamic, ', G=', TranInfo%lifetime_Dynamic
    CALL message(io8, FALSE, TRUE, MESG)
    WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, F7.1, x, a, F8.5, x, a, 1e10.3)') '@2', 0.,"s, r(pcm)= ", TranInfo%delrho_Dynamic/10._8, ', k= ', TranInfo%TranEig_Dynamic/10._8, ', B= ', TranInfo%CoreBeta_Dynamic
    CALL message(io8, FALSE, TRUE, MESG)
  ELSE
    WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, 1e14.7, x, a, 1e14.7, x, a, 1e14.7)') '@3', 0.,"s, P=", TranInfo%PowerLevel * 100, ', $= ', TranInfo%reactivity_Dynamic, ', G=', TranInfo%lifetime_dynamic
    CALL message(io8, FALSE, TRUE, MESG)
    WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, F7.1, x, a, F8.5, x, a, 1e14.7)') '@4', 0.,"s, r(pcm)= ", &
      TranInfo%delrho_Dynamic/10._8, ', k= ', TranInfo%TranEig_Dynamic/10._8, ', B= ', TranInfo%CoreBeta_Dynamic
    CALL message(io8, FALSE, TRUE, MESG)
    IF(nTracerCntl%lFeedback) THEN
      WRITE(mesg, '(a, 1p, e12.4, x, a, f10.3, x, a, f10.3, x, a)') '@5', 0.,"s, Max Tf= ", ThInfo%TfMax/10, "C, Max TfAvg= ", ThInfo%Tfavg_max/10, "C"
      CALL message(io8, FALSE, TRUE, MESG)
      WRITE(mesg, '(a, 1p, e12.4, x, a, f10.3, x, a, f10.3, x, a)') '@6', 0.,"s, Max Tcl= ", ThInfo%Tfcl_max/10, "C, Max TMOA= ", ThInfo%TModoutAvg/10, "C"
      CALL message(io8, FALSE, TRUE, MESG)
    ENDIF
  ENDIF
  WRITE(mesg, '(A)') hbar2(1:77)
  CALL message(io8, FALSE, TRUE, MESG)
ENDIF
CALL InitXsPerturbation(Core, TranCntl, PE)

tnow = 0.
dtpke = TranCntl%Tdiv_Corrector
dtcmfd = TranCntl%Tdiv_inp(1)
dtth = dtcmfd

laccept = .TRUE.
istep = 0
DO WHILE(.TRUE.)
  IF(laccept) THEN
    istep = istep + 1
    tprev = tnow
    tnow = tprev + dtcmfd
    TranCntl%nowstep = istep
    IF((TranCntl%tend-tnow) .LT. 1.E-6) THEN
      tnow = TranCntl%tend
      dtcmfd = tnow - tprev
      lfinish = .TRUE.
      IF(dtth .GT. dtcmfd) THEN
        dtth = dtcmfd
      ELSE IF(dtcmfd - dtth .LT. 1.E-6) THEN
        dtth = dtcmfd
      END IF
    END IF
    CALL SaveXsPerturbation(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, PE, tnow, tprev)
  ELSE
    tnow = tprev + dtcmfd
    IF((TranCntl%tend-tnow) .LT. 1.E-6) THEN
      tnow = TranCntl%tend
      dtcmfd = tnow - tprev
      lfinish = .TRUE.
      IF(dtth .GT. dtcmfd) THEN
        dtth = dtcmfd
      ELSE IF(dtcmfd - dtth .LT. 1.E-6) THEN
        dtth = dtcmfd
      END IF
    END IF
  END IF

  TranCntl%T(istep) = tnow
  TranCntl%DelT(istep) = dtcmfd
  ItrCntl%lConv = .FALSE.
  CALL InitIterVar(Core, FMInfo, CmInfo, GroupInfo,  .FALSE., ItrCntl, nTracerCntl, PE)
  lcmfdconv = .FALSE.

  TranCntl%lMOCUpdt = .FALSE.

  CALL CUDA_amplitude_updt(Core, CmInfo, FmInfo, GroupInfo, ThInfo, TranInfo, TranCntl, nTracerCntl, PE, TDopSave, ampSave, tprev, tnow, dtth)
  IF(TranCntl%lIQSAA) THEN
    TranCntl%IQSAA_g(1,TranCntl%IQSAA_m) = TranInfo%Amp
    TranCntl%IQSAA_g(2,TranCntl%IQSAA_m) = TranInfo%Amptilt
    TranCntl%IQSAA_x(1,TranCntl%IQSAA_m) = TranInfo%Amp
    TranCntl%IQSAA_x(2,TranCntl%IQSAA_m) = TranInfo%Amptilt
  END IF
  IF(nTracerCntl%lFeedback .AND. nTracerCntl%lrestrmt) THEN
    CALL CUDASubGrpFsp(Core, FmInfo,THInfo, RayInfo, GroupInfo)
  ENDIF
  iter = 0
  DO WHILE(.TRUE.)
    iter = iter + 1
    CALL ScaleFmFlux(FmInfo, TranInfo)
    TranInfo%AmpRatio = TranInfo%Amp / ampSave
    cuTranCMInfo%Expo = TranInfo%AmpRatio
    cuTranCmInfo%Expo_Alpha = TranInfo%AmpTilt
    IF(cuCntl%lAxial .AND. l3dim) THEN
      cuAxial%Expo = TranInfo%AmpRatio
      cuAxial%Expo_Alpha = TranInfo%AmpTilt
    END IF
    TranInfo%fmExpo = TranInfo%AmpRatio
    TranInfo%fmExpo_Alpha = TranInfo%AmpTilt

    CALL CUDATranCMFDAcc(Core, CmInfo, FmInfo, TranInfo, TranCntl, .TRUE.)
    CALL cuTranReactivityUpdt_Adj(GroupInfo, TranInfo, TranCntl, l3dim, .FALSE.)
    CALL GetShape(TranInfo)
    ampP = TranInfo%lifetime_Dynamic * TranInfo%Factor_F * TranInfo%Inv_Factor_K0
    shapechg = ShapeChange2()
    WRITE(mesg, '(2(a15,ES14.7))') 'Amplitude',TranInfo%amp, 'Shape Chg.', shapechg
    IF(PE%master)  CALL message(io8, FALSE, TRUE, MESG)
    CALL RecoverXsPerturbation(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, PE, tnow, tprev)
    !IF(shapechg .GT. shpCrit) EXIT

    CALL CUDA_amplitude_updt(Core, CmInfo, FmInfo, GroupInfo, ThInfo, TranInfo, TranCntl, nTracerCntl, PE, TDopSave, ampSave, tprev, tnow, dtth)
    IF(TranCntl%lIQSAA) THEN
      DO im = 2, TranCntl%IQSAA_m
        TranCntl%IQSAA_g(:,im-1) = TranCntl%IQSAA_g(:,im)
      END DO
      TranCntl%IQSAA_g(1,TranCntl%IQSAA_m) = TranInfo%Amp
      TranCntl%IQSAA_g(2,TranCntl%IQSAA_m) = TranInfo%Amptilt
      IF(TranCntl%IQSAA_m .EQ. 2) THEN
        iqsaa_f(:,1) = iqsaa_f(:,2)
        iqsaa_f(:,2) = TranCntl%IQSAA_g(:,2) - TranCntl%IQSAA_x(:,2)
        iqsaa_h = iqsaa_f(:,2)
        iqsaa_j = iqsaa_f(:,2) - iqsaa_f(:,1)
        iqsaa_r = sqrt(iqsaa_h(1)**2 + iqsaa_h(2)**2)
        iqsaa_a = 0
        DO inewt = 1, 10
          denom = iqsaa_j(1)**2 + iqsaa_j(2)**2
          numer = iqsaa_j(1)*iqsaa_h(1) + iqsaa_j(2)*iqsaa_h(2)
          iqsaa_a = iqsaa_a - numer / denom
          iqsaa_h(:) = iqsaa_f(:,2) + iqsaa_a * iqsaa_j(:)
          iqsaa_rd = iqsaa_r
          iqsaa_r = sqrt(iqsaa_h(1)**2 + iqsaa_h(2)**2)
          PRINT*, inewt, iqsaa_r, iqsaa_rd, iqsaa_a
          IF(iqsaa_r .LT. 1.E-3) EXIT
        END DO
        TranInfo%Amp = iqsaa_a * TranCntl%IQSAA_g(1,1) + (1-iqsaa_a) * TranCntl%IQSAA_g(1,2)
        TranInfo%Amptilt = iqsaa_a * TranCntl%IQSAA_g(2,1) + (1-iqsaa_a) * TranCntl%IQSAA_g(2,2)
        PRINT*, TranCntl%IQSAA_g(1,1), TranCntl%IQSAA_x(1,2)
        PRINT*, TranInfo%Amp, TranCntl%IQSAA_g(1,2)
      ELSE
        CALL terminate('high order AA is under construction')
      END IF
      DO im = 2, TranCntl%IQSAA_m
        TranCntl%IQSAA_x(:,im-1) = TranCntl%IQSAA_x(:,im)
      END DO
      TranCntl%IQSAA_x(1,TranCntl%IQSAA_m) = TranInfo%Amp
      TranCntl%IQSAA_x(2,TranCntl%IQSAA_m) = TranInfo%Amptilt
    END IF
    corrfactor = TranInfo%Amp/ampP - 1.
    WRITE(mesg, '(a15, ES11.3, 2(a8,ES14.7))') 'EPKE Correction', corrfactor, 'Amp_C:', TranInfo%Amp, 'Amp_P:', ampP
    IF (PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)
    !IF(abs(corrfactor) .LT. 1.E-5) lcmfdconv = .TRUE.
    IF(ItrCntl%lconv) lcmfdconv = .TRUE.
    IF(iter .GT. 10) lcmfdconv = .TRUE.
    IF(lcmfdconv) EXIT
  END DO

  !CALL TimeStepSelection(dtcmfd, shapechg, shpCrit, laccept)
  IF(laccept) THEN
    CALL cuCorrectNowStep(Core, CmInfo, FmInfo, GroupInfo, TranInfo, TranCntl, nTracerCntl)
    CALL cuUpdtPowerLevel(TranInfo, .FALSE.)
    nTracerCntl%PowerLevel = TranInfo%PowerLevel; ThInfo%PowLv = nTracerCntl%PowerLevel

    ampSave = TranInfo%amp
    IF(nTracerCntl%lFeedback) THEN
      TDopSave(1:Core%nz, 1:Core%nxy) = THInfo%Tdop(1:Core%nz, 1:Core%nxy)
      CALL SaveTranTHsol(Core, THInfo, TranCntl, nTracerCntl, PE)
    END IF
    CALL cuUpdtPrec(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, nprec)
    CALL cuSaveTranSol(FmInfo, nTracerCntl)
    CALL cuSavePKEParameters(TranInfo, TranCntl)
    CALL PrintTranInfo(PE, TranInfo, TranCntl, ThInfo, nTracerCntl, ItrCntl, tnow)
    CALL SaveShape()
    CALL cuUpdtResSrc(Core, FmInfo, cuCMFD, cuDevice, TranInfo, TranCntl, l3dim)
    IF(lfinish) EXIT
  END IF
END DO

END SUBROUTINE

#endif
