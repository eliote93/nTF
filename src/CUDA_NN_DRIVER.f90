!#include <defines.h>
!#include <CUDADEFINES.h>
!#ifdef __PGI
!SUBROUTINE CUDANNFSP_Driver()
!USE PARAM
!USE geom,           ONLY : Core
!USE RAYS,           ONLY : RayInfo
!USE Core_mod,       ONLY : FmInfo,      CmInfo,     GroupInfo,     eigv
!USE PE_Mod,         ONLY : PE
!USE CNTL,           ONLY : nTracerCntl
!USE TRAN_MOD,       ONLY : TranCntl
!USE CUDA_MASTER
!USE CUDA_NNFSP
!IMPLICIT NONE
!INTEGER :: i
!
!IF(nTracerCntl%lDynamicBen .OR. nTracerCntl%lKineticBen .OR. nTracerCntl%libtyp .EQ. 11) THEN
!  nTracerCntl%lchidgen = .FALSE.
!  nTracerCntl%lchidkgen = .FALSE.
!END IF
!TranCntl%lchidk = nTracerCntl%lchidkgen
!
!WRITE(mesg, '(a)') 'Initialize Neutron Noise Fixed Source Problem Calculation...'
!IF (PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)
!CALL cuInitNNFSP(Core, FmInfo, CmInfo, TranCntl, nTracerCntl, eigv)
!
!DO i = 1, TranCntl%nfreq
!  TranCntl%freq_now = TranCntl%freq_inp(i)
!  CALL SetNNSource(Core, FmInfo, TranCntl)
!
!
!
!
!
!END DO
!
!END SUBROUTINE
!
!SUBROUTINE cuNNMOCSweep(Core, RayInfo, FmInfo, eigv0)
!USE PARAM
!USE TYPEDEF,          ONLY: coreinfo_type,    RayInfo_Type,   FMInfo_TYPE, &
!                            Fxrinfo_type
!USE PE_Mod,           ONLY: PE
!USE CNTL,             ONLY: nTracerCntl
!USE CUDA_MOC
!USE CUDA_MASTER
!USE CUDA_INIT,        ONLY: AllocMOCVar,      DeallocMOCVar,  CopyRayVar,     DeleteRayVar
!USE files,            ONLY: io8
!USE ioutil,           ONLY: message
!USE timer,            ONLY: nTracer_dclock,   TimeChk
!USE omp_lib
!USE CUDAFOR
!USE OPENACC
!USE MOC_COMMON,       ONLY: SetMOCtrXS
!USE MOC_MOD,          ONLY: SetRtMacXsNM_Cusping
!USE SPH_Mod,          ONLY: calcPinSSPH
!IMPLICIT NONE
!TYPE(coreinfo_type) :: Core
!TYPE(RayInfo_Type) :: RayInfo
!TYPE(FMInfo_TYPE) :: FmInfo
!REAL :: eigv0
!
!TYPE(Fxrinfo_type), POINTER :: Fxr(:,:)
!REAL :: inv_eigv0
!REAL :: MocTimeBeg, MocTimeEnd, TimeRtBeg, TimeRtEnd
!INTEGER :: ngBlock, nPolarAngle, nPhiAngSv, nFsr, nxy, myzb, myze, ng
!INTEGER :: igb, ig, ifsr
!LOGICAL :: lxslib, lsSPH, lsSPHreg, lTrCorrection, lscat1, lrst
!!--- CUDA Host Variables --------------------------------------------------------------------------
!TYPE hostMOC_Type
!  REAL(GPU_PRECISION), ALLOCATABLE, PINNED :: PhiAngIn(:, :, :)
!  REAL(GPU_FLUX_PRECISION), ALLOCATABLE, PINNED :: phis(:, :), phim(:, :, :), jout(:, :, :, :)
!  REAL(GPU_SOURCE_PRECISION), ALLOCATABLE, PINNED :: xst(:, :), src(:, :), srcm(:, :, :)
!END TYPE
!
!TYPE(hostMOC_Type), POINTER :: hostMOC(:)
!REAL, POINTER :: xst(:, :), src(:, :), srcm(:, :, :)
!!--------------------------------------------------------------------------------------------------
!
!MocTimeBeg = nTracer_dclock(FALSE, FALSE)
!
!Fxr => FmInfo%Fxr
!
!myzb = PE%myzb
!myze = PE%myze
!nxy = cuGeometry%nxy
!nFsr = cuGeometry%nfsr
!ng = cuGeometry%ng
!nPolarAngle = RayInfo%nPolarAngle
!nPhiAngSv = RayInfo%nPhiAngSv
!
!lxslib = nTracerCntl%lXsLib
!lsSPH = nTracerCntl%lSSPH
!lsSPHreg = nTracerCntl%lSSPHreg
!lscat1 = nTracerCntl%lScat1
!lTrCorrection = nTracerCntl%lTrCorrection
!lrst = nTracerCntl%lRST
!
!IF (lsSPHreg) THEN
!  WRITE(mesg, '(a24)') 'Calculating Pin SSPH...'
!  IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)
!  CALL calcPinSSPH(Core, Fxr, PE)
!END IF
!
!itrcntl%mocit = itrcntl%mocit + 1
!WRITE(mesg, '(a32,I5,a3)') 'Performing Ray Tracing for NNFSP', itrcntl%mocit, '...'
!IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)
!
!CALL omp_set_nested(.TRUE.)
!CALL omp_set_max_active_levels(2)
!
!TimeRtBeg = nTracer_dclock(FALSE, FALSE)
!
!CALL acc_set_device_num(MPI_CUDA_SHARED_RANK, acc_device_nvidia)
!CALL omp_set_num_threads(cuCntl%nMOCHelper)
!
!CALL AllocMOCVar(Core, RayInfo, cuMOC, cuDevice, 2*ng, nTracerCntl%lScat1, cuCntl%lAsync, .NOT. cuCntl%lAsync)
!
!IF(.NOT. nTracerCntl%lScat1) THEN
!  ngBlock = (2*ng-1) / P0_BLOCK_SIZE + 1
!  ALLOCATE(hostMOC(ngBlock))
!  DO igb = 1, ngBlock
!    ALLOCATE(hostMOC(igb)%phis(P0_BLOCK_SIZE, nFsr))
!    ALLOCATE(hostMOC(igb)%PhiAngIn(nPolarAngle, P0_BLOCK_SIZE, nPhiAngSv))
!    ALLOCATE(hostMOC(igb)%jout(3,P0_BLOCK_SIZE, 4, nxy))
!    ALLOCATE(hostMOC(igb)%xst(P0_BLOCK_SIZE, nFsr))
!    ALLOCATE(hostMOC(igb)%src(P0_BLOCK_SIZE, nFsr))
!  END DO
!
!  ALLOCATE(xst(ng, nFsr))
!  ALLOCATE(src(P0_BLOCK_SIZE, nFsr))
!
!  DO iz = cuDevice%myzb, cuDevice%myze
!    IF(.NOT. cuDevice%lRayStatic) CALL CopyRayVar(cuGeometry%AxType)
!    IF(nTracerCntl%lMacro) THEN
!      CALL SetMOCTrXS(Core, Fxr, xst, iz, lsSPH, lsSPHreg, lTrCorrection)
!    ELSE
!      CALL SetRtMacXsNM_Cusping(Core, FmInfo, Fxr, xst, iz, ng, lxslib, lTrCorrection, lrst, lssph, lsSPHreg, PE)
!    END IF
!  END DO
!
!ELSE
!END IF
!
!NULLIFY(Fxr)
!
!MocTimeEnd = nTracer_dclock(FALSE, FALSE)
!TimeChk%MocTime = TimeChk%MocTime + MocTimeEnd - MocTimeBeg
!
!END SUBROUTINE
!
!#endif
