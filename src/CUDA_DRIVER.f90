#include <defines.h>
#include <CUDADEFINES.h>

#ifdef __PGI

SUBROUTINE CUDAInitialize(Core, RayInfo, FmInfo)
USE PARAM
USE TYPEDEF,		ONLY : CoreInfo_Type,		RayInfo_Type,		 FmInfo_Type
USE CUDA_MASTER
USE CUDA_INIT
USE CNTL,			ONLY : nTracerCntl
USE FILES,			ONLY : io8
USE IOUTIL,			ONLY : message
USE CORE_MOD,		ONLY : GroupInfo,           GcGroupInfo
USE PE_MOD,			ONLY : PE
USE CUDAFOR
USE CMFD_COMMON,    ONLY : AllocHomoXSVar,      AllocPinXS
USE CUBLAS,			ONLY : cublasCreate,		cublasSetStream
USE CUSPARSE,		ONLY : cusparseCreate,		cusparseSetStream
USE OPENACC
USE OMP_LIB
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo
TYPE(FmInfo_Type) :: FmInfo

TYPE(cudaDeviceProp) :: cudaProperty

INTEGER :: i, j, ierr
INTEGER :: nDevice, totalDevice, nz
INTEGER(8) :: totalByte, allocByte, freeByte(2)
CHARACTER(10) :: totalMBChar, allocMBChar, freeMBChar
REAL :: totalMB, allocMB, freeMB

nz = PE%myze - PE%myzb + 1

WRITE(mesg, '(a)') 'Searching for CUDA Capable GPUs...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)

CALL SetMPIEnv()

nDevice = acc_get_num_devices(acc_device_nvidia)
nDevice = min(1, nDevice)
cuCntl%nMOCHelper = PE%nThread
cuCntl%nCMFDHelper = PE%nCMFDThread
PE%nThread = cuCntl%nMOCHelper
PE%nCMFDThread = cuCntl%nCMFDHelper

CALL MPI_ALLREDUCE(nDevice, totalDevice, 1, MPI_INTEGER, MPI_SUM, MPI_CUDA_COMM, ierr)

WRITE(mesg, '(2x, i2, 1x, a)') totalDevice, 'CUDA Enabled GPUs Detected...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)

IF (totalDevice .NE. NUM_CUDA_PROC) THEN
  PE%lCUDA = .FALSE.;
  WRITE(mesg, '(a)') 'GPU Acceleration Option is Disabled...'
  IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)
  RETURN
ENDIF

WRITE(mesg, '(a)') 'Preparing GPU Environment...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)

WRITE(mesg, '(2x, a)') 'Checking Axial Heterogeneity...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)

CALL AxialCheck(Core)
CALL AssignMOCPlane(Core, cuDevice)

WRITE(mesg, '(2x, a)') 'Generating CUDA Rays...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)

CALL CUDAFastRay1DGen(Core, RayInfo)

WRITE(mesg, '(2x, a)') 'Setting CUDA Geometries...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)

CALL SetGeometry(Core, RayInfo, FmInfo, GroupInfo, cuDevice)
IF(cuCntl%lGcCMFD) CALL SetGcInfo(GroupInfo, GcGroupInfo)

CALL acc_set_device_num(MPI_CUDA_SHARED_RANK, acc_device_nvidia)

!--- Maximize L1 Cache
ierr = cudaDeviceSetCacheConfig(cudaFuncCachePreferL1)

ierr = cudaMemGetInfo(freeByte(1), totalByte)

WRITE(mesg, '(2x, a)') 'Setting Constant Data...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)

CALL SetConstVar(Core, RayInfo, GroupInfo)

WRITE(mesg, '(2x, a)') 'Setting Device Control Variables...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)

! ierr = cudaStreamCreate(cuDevice%myStream)
cuDevice%myStream = 0
CALL acc_set_cuda_stream(0, cuDevice%myStream)

CALL SetControlVar(Core, RayInfo, cuDevice)

IF (PE%lCUDACMFD) THEN

  WRITE(mesg, '(2x, a)') 'Allocating CMFD Variables...'
  IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)

  ierr = cublasCreate(cuDevice%myblasHandle)
  ierr = cusparseCreate(cuDevice%mySparseHandle)
  ierr = cublasSetStream(cuDevice%myblasHandle, cuDevice%myStream)
  ierr = cusparseSetStream(cuDevice%mySparseHandle, cuDevice%myStream)

  cuCMFD%ng = cuGeometry%ng
  cuGcCMFD%ng = cuGeometry%ngc

  CALL AllocHomoXSVar(Core, cuGeometry%ng)

  CALL AllocCMFDVar(cuCMFD, cuDevice, GroupInfo, .FALSE.)
  IF (cuCntl%lGcCMFD) CALL AllocCMFDVar(cuGcCMFD, cuDevice, GcGroupInfo, .TRUE.)

  IF (nTracerCntl%l3dim) THEN

    WRITE(mesg, '(2x, a)') 'Allocating Axial Solver Variables...'
    IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)

    CALL AllocAxialVar(cuAxial, cuDevice)

  ENDIF

ENDIF

WRITE(mesg, '(2x, a)') 'Copying Core Data...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)

CALL CopyCoreVar(Core, cuDevice)

IF (cuDevice%lRayStatic) THEN

  WRITE(mesg, '(2x, a)') 'Copying Ray Data...'
  IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)

  CALL CopyRayVar(cuDevice%myAxType)

ENDIF

WRITE(mesg, '(2x, a)') 'Copying Device Control Data...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)

CALL CopyControlVar()

!$ACC WAIT

ierr = cudaMemGetInfo(freeByte(2), totalByte)

allocByte = freeByte(1) - freeByte(2)

DO i = 0, NUM_CUDA_PROC - 1
  IF (i .EQ. MPI_CUDA_RANK) THEN
    totalMB = DBLE(totalByte) / 1024.0 ** 2
    allocMB = DBLE(allocByte) / 1024.0 ** 2
    freeMB = DBLE(freeByte(2)) / 1024.0 ** 2
    ierr = cudaGetDeviceProperties(cudaProperty, MPI_CUDA_SHARED_RANK)
    WRITE(totalMBChar, '(f8.2)') totalMB; totalMBChar = ADJUSTL(totalMBChar)
    WRITE(allocMBChar, '(f8.2)') allocMB; allocMBChar = ADJUSTL(allocMBChar)
    WRITE(freeMBChar, '(f8.2)') freeMB; freeMBChar = ADJUSTL(freeMBChar)
    WRITE(*, '(15x, a)'), TRIM(cudaProperty%name) // ' : ' // TRIM(totalMBChar) // 'MB Total, ' //                  &
                          TRIM(allocMBChar) // 'MB Allocated, ' // TRIM(freeMBChar) // 'MB Free'
  ENDIF
  CALL MPI_BARRIER(MPI_CUDA_COMM, ierr)
ENDDO

END SUBROUTINE

SUBROUTINE CUDASubGrpFsp(Core, FmInfo, THInfo, RayInfo, GroupInfo)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_type,   RayInfo_Type,       FmInfo_Type,        FxrInfo_type,                   &
                           THInfo_Type,     GroupInfo_Type
USE CNTL,           ONLY : nTracerCntl
USE xslib_mod,      ONLY : mlgdata,         mlgdata0
USE SUBGRP_MOD,     ONLY : UpdtCoreIsoInfo
USE TH_MOD,         ONLY : Cal_RefFuelTemp
USE PE_MOD,         ONLY : PE
USE FILES,          ONLY : io8
USE SubGrpFspNM
USE CUDA_MASTER
USE CUDA_INIT,      ONLY : AllocMOCVar,     DeallocMOCVar,      CopyRayVar,         DeleteRayVar
USE CUDA_MOC,       ONLY : CUDARayTrace
USE IOUTIL,         ONLY : message
USE TIMER,          ONLY : nTracer_dclock,          TimeChk
USE OMP_LIB
USE OPENACC
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(THInfo_Type) :: THInfo
TYPE(RayInfo_Type) :: RayInfo
TYPE(GroupInfo_Type) :: GroupInfo

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
INTEGER :: ig, iz, ilv, iter, itermax, itersum, ierr
INTEGER :: ng, nlv, nFsr, nFxr, nPhiAngSv, nPolarAngle
INTEGER :: gb, ge
REAL :: err, errmax, Tbeg, Tend, rtTime
REAL, POINTER :: phis(:, :), phisd(:, :), PhiAngIn(:, :, :)
REAL, POINTER :: Siglp(:, :), xst(:, :), src(:, :), jout(:, :, :, :)
LOGICAL :: lCLD, lAIC

IF (.NOT. nTracerCntl%lxslib) RETURN

Tbeg = nTracer_dclock(.FALSE., .FALSE.)

Fxr => FmInfo%Fxr
nFxr = Core%nCoreFxr
nFsr = Core%nCoreFsr
nPolarAngle = RayInfo%nPolarAngle
nPhiAngSv = RayInfo%nPhiAngSv
gb = GroupInfo%nofg + 1
ge = GroupInfo%nofg + GroupInfo%norg

iter = 0; itermax = 100
IF (any(cuGeometry%RadBC .EQ. Void)) itermax = 1

CALL Cal_RefFuelTemp(THInfo%RefFuelTemp, Core, Fxr, nTracerCntl, PE)
CALL UpdtCoreIsoInfo(Core, Fxr, PE)

WRITE(mesg,'(11x, a, f10.2, a)') "Reference Fuel Temperature", THInfo%RefFuelTemp(0), "C"
IF (PE%MASTER) CALL message(io8, .TRUE., .TRUE., mesg)
WRITE(mesg,'(a)') 'Solving Subgroup FSP...'
IF (PE%MASTER) CALL message(io8, .TRUE., .TRUE., mesg)

CALL omp_set_nested(.TRUE.)
CALL omp_set_max_active_levels(2)

CALL acc_set_device_num(MPI_CUDA_SHARED_RANK, acc_device_nvidia)
CALL omp_set_num_threads(cuCntl%nMOCHelper)

rtTime = 0.0

DO iz = cuDevice%myzb, cuDevice%myze

  IF (.NOT. Core%lFuelPlane(iz)) CYCLE

  ng = mlgdata(iz)%f_nmaclv * GroupInfo%norg

  ALLOCATE(phis(ng, nFsr)); ALLOCATE(phisd(ng, nFsr)); ALLOCATE(PhiAngIn(nPolarAngle, ng, nPhiAngSv))
  ALLOCATE(Siglp(ng, nFxr)); ALLOCATE(xst(ng, nFsr)); ALLOCATE(src(ng, nFsr))
  phis = 1.0; PhiAngIn = 1.0; PhiAngIn(:, :, 1) = 0.0

  CALL AllocMOCVar(Core, RayInfo, cuMOC, cuDevice, ng, .FALSE., .FALSE., .FALSE.)
  IF (.NOT. cuDevice%lRayStatic) CALL CopyRayVar(cuGeometry%AxType(iz))

  CALL UpdtFnAdj_NM(Core, Fxr, iz, gb, ge)
  CALL SetPlnLsigP_MLG_NM(Core, Fxr, Siglp, xst, iz, gb, ge)
  CALL SetSubGrpSrc_NM(Core, Fxr, Siglp, xst, src, iz, 1, ng)
  DO iter = 1, itermax
    cuMOC%src = src; cuMOC%xst = xst
    CALL CopyFlux(phis, phisd, nFsr, ng)
    CALL CUDARayTrace(cuMOC, cuDevice, phis, jout, PhiAngIn, .FALSE., iz, 1, ng, rtTime)
    CALL EquipXSGen_MLG_NM(Core, Fxr, Siglp, xst, phis, iz, gb, ge)
    CALL UpdtFtAdj_NM(Core, Fxr, iz, gb, ge)
    err = SubGrpFspErr_NM(phis, phisd, nFsr, ng)
    IF (err .LT. 1.0D-03) EXIT
    IF (iter .EQ. itermax) EXIT
    CALL SetPlnLsigP_MLG_NM(Core, Fxr, Siglp, xst, iz, gb, ge)
    CALL SetSubGrpSrc_NM(Core, Fxr, Siglp, xst, src, iz, 1, ng)
  ENDDO

  DEALLOCATE(phis, phisd, PhiAngIn)
  DEALLOCATE(Siglp, xst, src)

  CALL DeallocMOCVar(cuMOC, .FALSE., .FALSE., .FALSE.)
  IF (.NOT. cuDevice%lRayStatic) CALL DeleteRayVar(cuGeometry%AxType(iz))

ENDDO

CALL MPI_REDUCE(err, errmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_CUDA_COMM, ierr)
CALL MPI_REDUCE(iter, itersum, 1, MPI_INTEGER, MPI_SUM, 0, MPI_CUDA_COMM, ierr)

WRITE(mesg,'(a, i9, 1p, E20.5)') 'Subgroup FSP (Fuel) ', itersum, errmax
IF (PE%MASTER) CALL message(io8, .TRUE., .TRUE., mesg)

DO iz = cuDevice%myzb, cuDevice%myze

  IF (.NOT. Core%lFuelPlane(iz)) CYCLE
  IF (.NOT. Core%lCladPlane(iz)) CYCLE

  lCLD = .TRUE.; lAIC = .FALSE.
  ng = mlgdata0%c_nmaclv1G

  ALLOCATE(phis(ng, nFsr)); ALLOCATE(phisd(ng, nFsr)); ALLOCATE(PhiAngIn(nPolarAngle, ng, nPhiAngSv))
  ALLOCATE(Siglp(ng, nFxr)); ALLOCATE(xst(ng, nFsr)); ALLOCATE(src(ng, nFsr))
  phis = 1.0; PhiAngIn = 1.0; PhiAngIn(:, :, 1) = 0.0

  CALL AllocMOCVar(Core, RayInfo, cuMOC, cuDevice, ng, .FALSE., .FALSE., .FALSE.)
  IF (.NOT. cuDevice%lRayStatic) CALL CopyRayVar(cuGeometry%AxType(iz))

  CALL SetPlnLsigP_1gMLG_NM(Core, Fxr, Siglp, xst, iz, lCLD, lAIC)
  CALL SetSubGrpSrc_NM(Core, Fxr, Siglp, xst, src, iz, 1, ng)
  DO iter = 1, itermax
    cuMOC%src = src; cuMOC%xst = xst
    CALL CopyFlux(phis, phisd, nFsr, ng)
    CALL CUDARayTrace(cuMOC, cuDevice, phis, jout, PhiAngIn, .FALSE., iz, 1, ng, rtTime)
    err = SubGrpFspErr_NM(phis, phisd, nFsr, ng)
    IF (err .LT. 1.0D-03) EXIT
    IF (iter .EQ. itermax) EXIT
  ENDDO
  CALL EquipXSGen_1gMLG_NM(Core, Fxr, Siglp, phis, xst, iz, ng, lCLD, lAIC)

  DEALLOCATE(phis, phisd, PhiAngIn)
  DEALLOCATE(Siglp, xst, src)

  CALL DeallocMOCVar(cuMOC, .FALSE., .FALSE., .FALSE.)
  IF (.NOT. cuDevice%lRayStatic) CALL DeleteRayVar(cuGeometry%AxType(iz))

ENDDO

CALL MPI_REDUCE(err, errmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_CUDA_COMM, ierr)
CALL MPI_REDUCE(iter, itersum, 1, MPI_INTEGER, MPI_SUM, 0, MPI_CUDA_COMM, ierr)

WRITE(mesg,'(a, i9, 1p, E20.5)') 'Subgroup FSP (Clad) ', itersum, errmax
IF (PE%MASTER) CALL message(io8, .TRUE., .TRUE., mesg)

DO iz = cuDevice%myzb, cuDevice%myze

  IF (.NOT. Core%lFuelPlane(iz)) CYCLE
  IF (.NOT. Core%lAICPlane(iz)) CYCLE

  lCLD = .FALSE.; lAIC = .TRUE.
  ng = mlgdata(iz)%f_nmaclv1G

  ALLOCATE(phis(ng, nFsr)); ALLOCATE(phisd(ng, nFsr)); ALLOCATE(PhiAngIn(nPolarAngle, ng, nPhiAngSv))
  ALLOCATE(Siglp(ng, nFxr)); ALLOCATE(xst(ng, nFsr)); ALLOCATE(src(ng, nFsr))
  phis = 1.0; PhiAngIn = 1.0; PhiAngIn(:, :, 1) = 0.0

  CALL AllocMOCVar(Core, RayInfo, cuMOC, cuDevice, ng, .FALSE., .FALSE., .FALSE.)
  IF (.NOT. cuDevice%lRayStatic) CALL CopyRayVar(cuGeometry%AxType(iz))

  CALL SetPlnLsigP_1gMLG_NM(Core, Fxr, Siglp, xst, iz, lCLD, lAIC)
  CALL SetSubGrpSrc_NM(Core, Fxr, Siglp, xst, src, iz, 1, ng)
  DO iter = 1, itermax
    cuMOC%src = src; cuMOC%xst = xst
    CALL CopyFlux(phis, phisd, nFsr, ng)
    CALL CUDARayTrace(cuMOC, cuDevice, phis, jout, PhiAngIn, .FALSE., iz, 1, ng, rtTime)
    err = SubGrpFspErr_NM(phis, phisd, nFsr, ng)
    IF (err .LT. 1.0D-03) EXIT
    IF (iter .EQ. itermax) EXIT
  ENDDO
  CALL EquipXSGen_1gMLG_NM(Core, Fxr, Siglp, phis, xst, iz, ng, lCLD, lAIC)

  DEALLOCATE(phis, phisd, PhiAngIn)
  DEALLOCATE(Siglp, xst, src)

  CALL DeallocMOCVar(cuMOC, .FALSE., .FALSE., .FALSE.)
  IF (.NOT. cuDevice%lRayStatic) CALL DeleteRayVar(cuGeometry%AxType(iz))

ENDDO

CALL MPI_REDUCE(err, errmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_CUDA_COMM, ierr)
CALL MPI_REDUCE(iter, itersum, 1, MPI_INTEGER, MPI_SUM, 0, MPI_CUDA_COMM, ierr)

WRITE(mesg,'(a, i9, 1p, E20.5)') 'Subgroup FSP (AIC)  ', itersum, errmax
IF (PE%MASTER) CALL message(io8, .TRUE., .TRUE., mesg)

nTracerCntl%lSubGrpSweep = .TRUE.
Tend = nTracer_dclock(.FALSE., .FALSE.)
TimeChk%SubGrpTime = TimeChk%SubGrpTime + (Tend - Tbeg)
TimeChk%NetRTSubGrpTime = TimeChk%NetRTSubGrpTime + rtTime

IF (PE%MASTER) CALL message(io8, .TRUE., .TRUE., 'Finished Subgroup FSP...')

END SUBROUTINE

SUBROUTINE CUDAMOCSweep(Core, RayInfo, FmInfo, eigv)
USE PARAM
USE TYPEDEF,		ONLY : CoreInfo_Type,       RayInfo_Type,       FmInfo_Type,        FxrInfo_Type
USE PE_MOD,         ONLY : PE
USE CNTL,			ONLY : nTracerCntl
USE ITRCNTL_MOD,	ONLY : ItrCntl
USE CORE_MOD,		ONLY : GroupInfo
USE MOC_MOD,		ONLY : SetRtMacXsNM,	    SetRtSrcNM,         SetRtP1SrcNM,       PsiUpdate,                  &
                           UpdateEigv,          MocResidual,        PsiErr,             PseudoAbsorptionNM,         &
                           PowerUpdate
USE CUDA_MASTER
USE CUDA_INIT,      ONLY : AllocMOCVar,         DeallocMOCVar,      CopyRayVar,         DeleteRayVar
USE CUDA_MOC
USE SUbGrp_Mod,     ONLY : FxrChiGen
USE IOUTIL,         ONLY : message
USE Timer,          ONLY : nTracer_dclock,      TimeChk
USE FILES,          ONLY : io8
#ifdef MPI_ENV
USE MPICOMM_MOD,    ONLY : BCAST,               MPI_SYNC,           MPI_MAX_REAL
#endif
USE MOC_COMMON,     ONLY : SetMOCPsi,           SetMOCtrXS,         SetMOCSource,       SetMOCPower,                &
                           CopyXS
USE SPH_mod,        ONLY : ssphfnm,             calcPinSSPH
USE XSLIB_MOD,      ONLY : igresb,              igrese
USE OMP_LIB
USE CUDAFOR
USE OPENACC
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo
TYPE(FmInfo_Type) :: FmInfo
REAL :: eigv

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
REAL, POINTER :: psi(:, :), psid(:, :), AxSrc(:, :, :), AxPXS(:, :, :)
REAL :: eigerr, fiserr, peigv, reserr, errdat(3)
REAL :: eigconv, psiconv
REAL :: rtTime, TimeRtBeg, TimeRtEnd, TimeRtElapsed, MocTimeBeg, MocTimeEnd
INTEGER :: ig, igb, ifsr, iray, ipin, isurf, iz, i, ierr, InIter, nInIter
INTEGER :: myzb, myze
INTEGER :: nFsr, nxy, nPhiAngSv, nPolarAngle, nGroupInfo, nMoment
INTEGER :: GrpBeg, GrpEnd, gb, ge, ng, ngBlock, ngLocal
LOGICAL :: lJout, lxslib, l3dim, lscat1, lTrCorrection, lRST, lsSPH, lsSPHreg

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

!--------------------------------------------------------------------------------------------------

MOCTimeBeg = nTracer_dclock(FALSE, FALSE)

Fxr => FmInfo%Fxr
psi => FmInfo%psi; psid => FmInfo%psid
AxSrc => FmInfo%AxSrc; AxPXS => FmInfo%AxPXS

myzb = PE%myzb; myze = PE%myze
nFsr = cuGeometry%nFsr; nxy = cuGeometry%nxy
nPolarAngle = RayInfo%nPolarAngle; nPhiAngSv = RayInfo%nPhiAngSv; ng = cuGeometry%ng
nInIter = 2; nGroupInfo = 2; IF (.NOT. GroupInfo%lUpScat) nGroupInfo = 1

IF (nTracerCntl%ScatOd .EQ. 1) nMoment = 2
IF (nTracerCntl%ScatOd .EQ. 2) nMoment = 5
IF (nTracerCntl%ScatOd .EQ. 3) nMoment = 9

lxslib = nTracerCntl%lXslib; l3dim = nTracerCntl%l3dim; lscat1 = nTracerCntl%lscat1
lRST = nTracerCntl%lRST; lsSPH = nTracerCntl%lsSPH; lsSPHreg = nTracerCntl%lsSPHreg
lTrCorrection = nTracerCntl%lTrCorrection

psiconv = itrcntl%psiconv; eigconv = itrcntl%eigconv

IF (lsSPHreg) CALL calcPinSSPH(Core, Fxr, PE)

psid(:, myzb : myze) = psi(:, myzb : myze)
IF (nTracerCntl%lMacro) THEN
  CALL SetMOCPsi(Core, FmInfo%phis, psi)
ELSE
  CALL PsiUpdate(Core, Fxr, FmInfo%phis, psi, myzb, myze, ng, lxslib, GroupInfo)
ENDIF
CALL FxrChiGen(Core, Fxr, FmInfo, GroupInfo, PE, nTracerCntl)

itrcntl%mocit = itrcntl%mocit + 1
WRITE(mesg, '(a22,I5,a3)') 'Performing Ray Tracing', itrcntl%mocit, '...'
IF (PE%MASTER) CALL message(io8, .TRUE., .TRUE., mesg)

CALL omp_set_nested(.TRUE.)
CALL omp_set_max_active_levels(2)

TimeRtBeg = nTracer_dclock(.FALSE., .FALSE.)

CALL acc_set_device_num(MPI_CUDA_SHARED_RANK, acc_device_nvidia)
CALL omp_set_num_threads(cuCntl%nMOCHelper)

CALL AllocMOCVar(Core, RayInfo, cuMOC, cuDevice, ng, nTracerCntl%lScat1, cuCntl%lAsync, .NOT. cuCntl%lAsync)

IF (.NOT. nTracerCntl%lScat1) THEN

  IF (.NOT. cuCntl%lAsync) THEN

    ALLOCATE(phis(ng, nFsr))
    ALLOCATE(xst(ng, nFsr))
    ALLOCATE(Jout(3, ng, 4, nxy))
    ALLOCATE(PhiAngIn(nPolarAngle, ng, nPhiAngSv))

    rtTime = 0.0

    DO iz = cuDevice%myzb, cuDevice%myze
      IF (.NOT. cuDevice%lRayStatic) CALL CopyRayVar(cuGeometry%AxType(iz))
      CALL CopyXS(Core, FmInfo, iz)
      DO ig = 1, ng
        phis(ig, :) = FmInfo%phis(:, iz, ig)
        PhiAngIn(:, ig, :) = FmInfo%PhiAngIn(:, :, iz, ig)
      ENDDO
      IF (nTracerCntl%lMacro .AND. .NOT. lsSPH) THEN
        CALL SetMOCtrXS(Core, xst, iz)
      ELSE
        CALL SetRtMacXsNM(Core, Fxr(:, iz), xst, iz, ng, lxslib, lTrCorrection, lRST, lsSPH, lsSPHreg, PE)
      ENDIF
      CALL PseudoAbsorptionNM(Core, Fxr(:, iz), AxPXS, xst, iz, ng, GroupInfo, l3dim)
      cuMOC%phis = phis; cuMOC%psi = psi(:, iz); cuMOC%xst = xst
      DO i = 1, nGroupInfo
        GrpBeg = 1; GrpEnd = ng
        IF (i .GT. 1) THEN
          GrpBeg = GroupInfo%UpScatRange(1); GrpEnd = GroupInfo%UpScatRange(2)
        ENDIF
        DO InIter = 1, nInIter
          lJout = .FALSE.
          IF (InIter .EQ. nInIter) lJout = .TRUE.
          CALL CUDASourceUpdate(iz, GrpBeg, GrpEnd, eigv)
          CALL CUDARayTrace(cuMOC, cuDevice, phis, Jout, PhiAngIn, lJout, iz, GrpBeg, GrpEnd, rtTime)
          IF (lsSPH) THEN
            IF (GrpBeg .LE. igrese) THEN
              gb = max(igresb, GrpBeg); ge = min(igrese, GrpEnd)
              phis(gb : ge, :) = phis(gb : ge, :) * ssphfnm(gb : ge, :, iz)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      DO ig = 1, ng
        FmInfo%phis(:, iz, ig) = phis(ig, :)
        FmInfo%PhiAngIn(:, :, iz, ig) = PhiAngIn(:, ig, :)
        FmInfo%RadJout(:, :, :, iz, ig) = Jout(:, ig, :, :)
      ENDDO
      IF (.NOT. cuDevice%lRayStatic) CALL DeleteRayVar(cuGeometry%AxType(iz))
    ENDDO

    DEALLOCATE(phis)
    DEALLOCATE(xst)
    DEALLOCATE(Jout)
    DEALLOCATE(PhiAngIn)

  ELSE

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
    ALLOCATE(Jout(3, ng, 4, nxy))
    ALLOCATE(xst(ng, nFsr))
    ALLOCATE(src(P0_BLOCK_SIZE, nFsr))

    DO iz = cuDevice%myzb, cuDevice%myze
      IF (.NOT. cuDevice%lRayStatic) CALL CopyRayVar(cuGeometry%AxType(iz))
      !$OMP PARALLEL DO COLLAPSE(2)
      DO ig = 1, ng
        DO ifsr = 1, nFsr
          phis(ig, ifsr) = FmInfo%phis(ifsr, iz, ig)
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
      IF (nTracerCntl%lMacro .AND. .NOT. lsSPH) THEN
        CALL SetMOCtrXS(Core, xst, iz)
      ELSE
        CALL SetRtMacXsNM(Core, Fxr(:, iz), xst, iz, ng, lxslib, lTrCorrection, lRST, lsSPH, lsSPHreg, PE)
      ENDIF
      CALL PseudoAbsorptionNM(Core, Fxr(:, iz), AxPXS, xst, iz, ng, GroupInfo, l3dim)
      DO i = 1, nGroupInfo
        GrpBeg = 1; GrpEnd = ng
        IF (i .GT. 1) THEN
          GrpBeg = GroupInfo%UpScatRange(1); GrpEnd = GroupInfo%UpScatRange(2)
        ENDIF
        ngBlock = (GrpEnd - GrpBeg) / P0_BLOCK_SIZE + 1
        DO InIter = 1, nInIter
          IF (InIter .EQ. nInIter) lJout = .TRUE.
          DO igb = 1, ngBlock
            gb = (igb - 1) * P0_BLOCK_SIZE + GrpBeg
            ge = min(igb * P0_BLOCK_SIZE + GrpBeg - 1, ng)
            ngLocal = ge - gb + 1
            IF (nTracerCntl%lMacro) THEN
              CALL SetMOCSource(Core, Fxr(:, iz), src, phis, psi, AxSrc, xst, eigv, iz, gb, ge, gb - 1)
            ELSE
              CALL SetRtSrcNM(Core, Fxr(:, iz), src, phis, psi, AxSrc, xst, eigv, iz,                               &
                              gb, ge, ng, GroupInfo, l3dim, lxslib, lscat1, .FALSE., PE, gb - 1)
            ENDIF
            !$OMP PARALLEL DO COLLAPSE(2)
            DO ifsr = 1, nFsr
              DO ig = 1, ngLocal
                hostMOC(igb)%xst(ig, ifsr) = xst(ig + gb - 1, ifsr)
                hostMOC(igb)%src(ig, ifsr) = src(ig, ifsr)
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
    DEALLOCATE(Jout)
    DEALLOCATE(xst)
    DEALLOCATE(src)

    ngBlock = (ng - 1) / P0_BLOCK_SIZE + 1
    DO igb = 1, ngBlock
      DEALLOCATE(hostMOC(igb)%phis)
      DEALLOCATE(hostMOC(igb)%PhiAngIn)
      DEALLOCATE(hostMOC(igb)%jout)
      DEALLOCATE(hostMOC(igb)%xst)
      DEALLOCATE(hostMOC(igb)%src)
    ENDDO
    DEALLOCATE(hostMOC)

  ENDIF

ELSE

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
  ALLOCATE(Jout(3, ng, 4, nxy))
  ALLOCATE(xst(ng, nFsr))
  ALLOCATE(src(PN_BLOCK_SIZE, nFsr))
  ALLOCATE(srcm(nMoment, PN_BLOCK_SIZE, nFsr))

  DO iz = cuDevice%myzb, cuDevice%myze
    IF (.NOT. cuDevice%lRayStatic) CALL CopyRayVar(cuGeometry%AxType(iz))
    !$OMP PARALLEL DO COLLAPSE(2)
    DO ig = 1, ng
      DO ifsr = 1, nFsr
        phis(ig, ifsr) = FmInfo%phis(ifsr, iz, ig)
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    phim => FmInfo%phim(:, :, :, iz)
    CALL SetRtMacXsNM(Core, Fxr(:, iz), xst, iz, ng, lxslib, lTrCorrection, lRST, lsSPH, lsSPHreg, PE)
    CALL PseudoAbsorptionNM(Core, Fxr(:, iz), AxPXS, xst, iz, ng, GroupInfo, l3dim)
    DO i = 1, nGroupInfo
      GrpBeg = 1; GrpEnd = ng
      IF (i .GT. 1) THEN
        GrpBeg = GroupInfo%UpScatRange(1); GrpEnd = GroupInfo%UpScatRange(2)
      ENDIF
      ngBlock = (GrpEnd - GrpBeg) / PN_BLOCK_SIZE + 1
      DO InIter = 1, nInIter
        IF (InIter .EQ. nInIter) lJout = .TRUE.
        DO igb = 1, ngBlock
          gb = (igb - 1) * PN_BLOCK_SIZE + GrpBeg
          ge = min(igb * PN_BLOCK_SIZE + GrpBeg - 1, ng)
          ngLocal = ge - gb + 1
          CALL SetRtSrcNM(Core, Fxr(:, iz), src, phis, psi, AxSrc, xst, eigv, iz,                                   &
                          gb, ge, ng, GroupInfo, l3dim, lxslib, lscat1, .FALSE., PE, gb - 1)
          CALL SetRtP1SrcNM(Core, Fxr(:, iz), srcm, phim, xst, iz, gb, ge, ng, GroupInfo,                           &
                            lxslib, nTracerCntl%ScatOd, PE, gb - 1)
          !$OMP PARALLEL DO COLLAPSE(2)
          DO ifsr = 1, nFsr
            DO ig = 1, ngLocal
              hostMOC(igb)%xst(ig, ifsr) = xst(ig + gb - 1, ifsr)
              hostMOC(igb)%src(ig, ifsr) = src(ig, ifsr)
              hostMOC(igb)%srcm(:, ig, ifsr) = srcm(:, ig, ifsr)
            ENDDO
          ENDDO
          !$OMP END PARALLEL DO
          hostMOC(igb)%PhiAngIn = 0.0
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
  DEALLOCATE(Jout)
  DEALLOCATE(xst)
  DEALLOCATE(src)
  DEALLOCATE(srcm)

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

TimeRtEnd = nTracer_dclock(.FALSE., .FALSE.)
TimeRtElapsed = TimeRtEnd - TimeRtBeg

#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_RTMASTER_COMM)
#endif

CALL MPI_MAX_REAL(TimeRtElapsed, PE%MPI_RTMASTER_COMM, .TRUE.)
write(mesg, '(3x, a, f10.3, 2x, a)') 'CUDA MOC Sweep Finished in ', TimeRtElapsed, 'Sec'
IF (PE%MASTER) CALL message(io8, .TRUE., .TRUE., mesg)

psid(:, myzb : myze) = psi(:, myzb : myze)
IF (nTracerCntl%lMacro) THEN
  CALL SetMOCPsi(Core, FmInfo%phis, psi)
ELSE
  CALL PsiUpdate(Core, Fxr, FmInfo%phis, psi, myzb, myze, ng, lxslib, GroupInfo)
ENDIF

CALL UpdateEigv(Core, psi, psid, eigv, peigv, myzb, myze, PE)
eigerr = abs((eigv - peigv)) / eigv; fiserr = zero; reserr = zero;
fiserr = PsiErr(Core, Psi, PsiD, myzb, myze, PE)

#ifdef MPI_ENV
errdat = (/fiserr, eigerr, reserr/)
CALL BCAST(errdat, 3, PE%MPI_RT_COMM)
fiserr = errdat(1); eigerr = errdat(2); reserr = errdat(3)
#endif
WRITE(mesg ,'(A5,I7,F15.6, 2(1pE12.3), A12)') 'RT', itrcntl%mocit, eigv, eigerr, fiserr, '   *********'
IF (PE%MASTER) CALL message(io8, .TRUE., .TRUE., mesg)

IF (nTracerCntl%lMacro) THEN
  CALL SetMOCPower(Core, FmInfo%phis, FmInfo%Power)
ELSE
  CALL PowerUpdate(Core, Fxr, FmInfo%phis, FmInfo%Power, myzb, myze, ng, lxslib, GroupInfo, PE)
ENDIF

ItrCntl%lConv = fiserr .LT. psiconv .AND. eigerr .LT. eigconv

MocTimeEnd = nTracer_dclock(.FALSE., .FALSE.)
TimeChk%MocTime = TimeChk%MocTime + MocTimeEnd - MocTimeBeg
TimeChk%MocRtTime = TimeChk%MocRtTime + rtTime

END SUBROUTINE

SUBROUTINE CUDACmfdAcc(Core, CmInfo, FmInfo, eigv)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,       CmInfo_Type,        FmInfo_Type,                                    &
                           FxrInfo_Type
USE CORE_MOD,       ONLY : GroupInfo
USE PE_MOD,         ONLY : PE
USE SUBGRP_MOD,     ONLY : FxrChiGen
USE IOUTIL,         ONLY : message
USE FILES,          ONLY : io8
USE CNTL,           ONLY : nTracerCntl
USE ITRCNTL_MOD,    ONLY : ItrCntl
USE TIMER,          ONLY : nTracer_dclock,      TimeChk
USE CMFD_COMMON,    ONLY : HomogenizeXS,        SetRadialCoupling
USE CUDA_INIT,      ONLY : DeallocCMFDVar
USE CUDA_AXIAL,     ONLY : cuAxialSolver,       cuSetAxialDtil,     cuSetAxialSourceOperator
USE CUDA_CMFD
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(CMInfo_Type) :: CmInfo
TYPE(FMInfo_Type) :: FmInfo
REAL :: eigv

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
REAL, POINTER :: phis(:, :, :), phim(:, :, :, :), phic(:, :, :), phia(:, :, :, :)
REAL, POINTER :: Jout(:, :, :, :, :), AxSrc(:, :, :), AxPXS(:, :, :)
REAL :: CmfdTimeBeg, CmfdTimeEnd, CmfdInitBeg, CmfdInitEnd
REAL :: outTol = 0.1, outErr, outRes, outRes0
INTEGER :: ng, nxy, nInIter
INTEGER :: myzb, myze
INTEGER :: outIter
INTEGER :: ierr, iz
LOGICAL :: lxslib, l3dim, lscat1, lDhat
LOGICAL :: loutConv

CmfdTimeBeg = nTracer_dclock(.FALSE., .FALSE.)

l3dim = nTracerCntl%l3dim
lxslib = nTracerCntl%lxslib
lscat1 = nTracerCntl%lscat1
Fxr => FmInfo%Fxr
phis => FmInfo%phis
phim => FmInfo%phim
phic => CmInfo%phic
Jout => CmInfo%RadJout
AxSrc => FmInfo%AxSrc
AxPXS => FmInfo%AxPXS

lDhat = ItrCntl%CMFDIt .NE. ItrCntl%CMFDIt0
ng = cuGeometry%ng; nxy = cuGeometry%nxyc; nInIter = 1
myzb = cuDevice%myzb; myze = cuDevice%myze

CALL omp_set_num_threads(PE%nCMFDThread)
CALL omp_set_nested(.TRUE.)
CALL omp_set_max_active_levels(2)

CmfdInitBeg = nTracer_dclock(.FALSE., .FALSE.)

WRITE(mesg, '(a)') 'Cell Homogenization...'
IF (PE%Master) CALL message(io8, .TRUE., .TRUE., mesg)
CALL FxrChiGen(Core, Fxr, FmInfo, GroupInfo, PE, nTracerCntl)

CALL HomogenizeXS(Core, cuGeometry%superPin, Fxr, cuCMFD%PinXS, phis, ng, nxy, myzb, myze, lxslib, lscat1, .FALSE.)
CALL SetRadialCoupling(cuGeometry%superPin, cuCMFD%PinXS, Jout, ng, nxy, myzb, myze, lDhat)

CALL acc_set_device_num(MPI_CUDA_SHARED_RANK, acc_device_nvidia)
CALL omp_set_num_threads(cuCntl%nCMFDHelper)

CALL cuSetCMFDPhis(cuCMFD, cuDevice, lFirstCMFD)

IF (l3dim) CALL cuSetAxialDtil(cuCMFD, cuAxial, cuDevice)

CmfdInitEnd = nTracer_dclock(.FALSE., .FALSE.)
TimeChk%CmfdInitTime = TimeChk%CmfdInitTime + (CmfdInitEnd - CmfdInitBeg)

IF (cuCntl%lAxial) THEN
  IF (l3dim .AND. .NOT. lFirstCMFD) THEN
    CALL cuGetNeighborFlux(cuCMFD, cuDevice)
    CALL cuSetAxialSourceOperator(Core, FmInfo, GroupInfo, cuCMFD, cuAxial, cuDevice)
    CALL cuAxialSolver(cuCMFD, cuAxial, cuDevice, eigv)
  ENDIF
ENDIF

IF (.NOT. cuCntl%lNatural) CALL cuReorderFlux(cuCMFD, cuDevice, 1)
CALL cuCopyFlux(cuCMFD, cuDevice, 1)

CmfdInitBeg = nTracer_dclock(.FALSE., .FALSE.)

WRITE(mesg, '(a)') 'Linear System Construction...'
IF (PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)

IF (cuCntl%lNatural) THEN
  IF (cuCntl%CMFDSolver .EQ. 1) CALL cuSetNaturalBiCGSystem(cuCMFD, cuDevice, l3dim, .FALSE., 0.0)
  IF (cuCntl%CMFDSolver .EQ. 2) CALL cuSetNaturalJacobiSystem(cuCMFD, cuDevice, l3dim)
ELSE
  CALL cuSetRedblackSORSystem(cuCMFD, cuDevice, l3dim)
ENDIF

CALL cuSetCMFDSourceOperator(cuCMFD, cuDevice, GroupInfo%lUpscat)

CmfdInitEnd = nTracer_dclock(.FALSE., .FALSE.)
TimeChk%CmfdInitTime = TimeChk%CmfdInitTime + (CmfdInitEnd - CmfdInitBeg)

WRITE(mesg, '(a)') 'Performing CMFD Acceleration...'
IF (PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)

CALL cuCMFDPsiUpdt(cuCMFD, cuDevice)
CALL cuCMFDSrcUpdt(cuCMFD, cuDevice, eigv, GroupInfo%lUpscat)

outIter = 0; loutConv = FALSE
DO WHILE (.NOT. loutConv)

  CALL cuInnerSolve(cuCMFD, cuDevice, 0.01)
  CALL cuCMFDPsiUpdt(cuCMFD, cuDevice)
  CALL cuCMFDEigUpdt(cuCMFD, cuDevice, eigv)
  CALL cuCMFDSrcUpdt(cuCMFD, cuDevice, eigv, GroupInfo%lUpscat)

  outRes = cuCMFDPsiErr(cuCMFD, cuDevice)
  IF (outIter .EQ. 0) outRes0 = outRes
  outErr = outRes / outRes0
  outIter = outIter + 1
  loutConv = (outErr .LE. outTol) .AND. (outIter .GE. 7)

  ItrCntl%CMFDIt = ItrCntl%CMFDIt + 1
  IF (PE%MASTER) WRITE(mesg, '(a9, i9, f22.6, 3x, f10.5, 1p, e15.3)') 'MGOUTER', ItrCntl%CMFDIt, eigv, outErr, outRes
  IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)

  IF (mod(outIter, 5) .NE. 0) CYCLE
  IF (loutConv) CYCLE

  IF (cuCntl%lAxial) THEN
    IF (l3dim .AND. .NOT. lFirstCMFD) THEN
      CALL cuCopyFlux(cuCMFD, cuDevice, 2)
      IF (.NOT. cuCntl%lNatural) CALL cuReorderFlux(cuCMFD, cuDevice, 2)
      CALL cuGetNeighborFlux(cuCMFD, cuDevice)
      CALL cuAxialSolver(cuCMFD, cuAxial, cuDevice, eigv)
    ENDIF
  ENDIF

  !$OMP PARALLEL NUM_THREADS(2)
  !$OMP MASTER
  IF (cuCntl%lGcCMFD) THEN
    CALL omp_set_num_threads(cuCntl%nCMFDHelper - 1)
    CALL cuCopyFlux(cuCMFD, cuDevice, 2)
    IF (.NOT. cuCntl%lNatural) CALL cuReorderFlux(cuCMFD, cuDevice, 2)
    CALL cuGetNeighborFlux(cuCMFD, cuDevice)
    CALL CUDAGcCmfdAcc(cuCMFD, cuGcCMFD, cuDevice, eigv)
    IF (.NOT. cuCntl%lNatural) CALL cuReorderFlux(cuCMFD, cuDevice, 1)
    CALL cuCopyFlux(cuCMFD, cuDevice, 1)
    CALL cuCMFDPsiUpdt(cuCMFD, cuDevice)
    CALL cuCMFDSrcUpdt(cuCMFD, cuDevice, eigv, GroupInfo%lUpscat)
  ENDIF
  CmfdInitBeg = nTracer_dclock(.FALSE., .FALSE.)
  !$OMP END MASTER
  !$OMP SINGLE
  IF (cuCntl%lAxial) THEN
    IF (l3dim .AND. .NOT. lFirstCMFD) THEN
      CALL omp_set_num_threads(1)
      CALL acc_set_device_num(MPI_CUDA_SHARED_RANK, acc_device_nvidia)
      WRITE(mesg, '(a)') 'Linear System Construction...'
      IF (PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)
      IF (cuCntl%lNatural) THEN
        IF (cuCntl%CMFDSolver .EQ. 1) CALL cuSetNaturalBiCGSystem(cuCMFD, cuDevice, l3dim, .FALSE., 0.0)
        IF (cuCntl%CMFDSolver .EQ. 2) CALL cuSetNaturalJacobiSystem(cuCMFD, cuDevice, l3dim)
      ELSE
        CALL cuSetRedblackSORSystem(cuCMFD, cuDevice, l3dim)
      ENDIF
    ENDIF
  ENDIF
  !$OMP END SINGLE
  !$OMP END PARALLEL
  CmfdInitEnd = nTracer_dclock(.FALSE., .FALSE.)
  TimeChk%CmfdInitTime = TimeChk%CmfdInitTime + (CmfdInitEnd - CmfdInitBeg)

ENDDO

CALL cuCopyFlux(cuCMFD, cuDevice, 2)
IF (.NOT. cuCntl%lNatural) CALL cuReorderFlux(cuCMFD, cuDevice, 2)
CALL cuSetMOCPhis(Core, cuCMFD, cuDevice, phis, phim, phic, lScat1)

IF (l3dim) THEN
  CALL cuGetNeighborFlux(cuCMFD, cuDevice)
  CALL cuSetAxialSrc(cuCMFD, cuDevice, AxSrc, AxPXS, phic)
ENDIF

CALL DeallocCMFDVar(cuCMFD)

IF (lFirstCMFD) lFirstCMFD = .FALSE.

CmfdTimeEnd = nTracer_dclock(.FALSE., .FALSE.)
TimeChk%CmfdTime = TimeChk%CmfdTime + (CmfdTimeEnd - CmfdTimeBeg)

END SUBROUTINE

SUBROUTINE CUDAGcCmfdAcc(cuCMFD, cuGcCMFD, cuDevice, Keff)
USE PARAM
USE PE_MOD,         ONLY : PE
USE CORE_MOD,       ONLY : GcGroupInfo
USE IOUTIL,         ONLY : message
USE FILES,          ONLY : io8
USE ITRCNTL_MOD,    ONLY : ItrCntl
USE CUDA_INIT,      ONLY : DeallocCMFDVar
USE CUDA_AXIAL,     ONLY : cuSetAxialGcDtil,    cuSetAxialGcDhat
USE CUDA_CMFD
USE CUDA_GCCMFD
IMPLICIT NONE

TYPE(cuCMFD_Type) :: cuCMFD, cuGcCMFD
TYPE(cuDevice_Type) :: cuDevice
REAL :: Keff

REAL :: outTol = 0.1, outErr, outRes, outRes0
REAL :: eigv, seigv, tol
INTEGER :: outIter, ierr
LOGICAL :: l3dim, loutConv, lNatural, lDcpl

l3dim = cuGeometry%l3dim
lNatural = cuCntl%lNatural
lDcpl = cuCntl%lDcpl

WRITE(mesg, '(a)') 'Group Condensing...'
IF (PE%Master) CALL message(io8, .TRUE., .TRUE., mesg)

IF (cuCntl%lShift) THEN
  cuCntl%lNatural = .TRUE.
  cuCntl%lDcpl = .FALSE.
  seigv = 1.0 / (Keff + cuCntl%Shift)
  eigv = 1.0 / (1.0 / Keff - seigv)
  tol = 0.1
  WRITE(mesg, '(a27, f4.2)') 'Wielandt Shift Parameter : ', cuCntl%Shift
  IF (PE%Master) CALL message(io8, TRUE, TRUE, mesg)
ELSE
  seigv = 0.0
  eigv = Keff
  tol = 0.01
ENDIF

CALL cuHomogenizeGcXS(cuCMFD, cuGcCMFD, cuDevice)
CALL cuSetRadialGcCoupling(cuCMFD, cuGcCMFD, cuDevice)

CALL cuSetCMFDPhis(cuGcCMFD, cuDevice, .TRUE.)

IF (l3dim) THEN
  CALL cuGetNeighborFlux(cuGcCMFD, cuDevice)
  CALL cuSetAxialGcDtil(cuGcCMFD, cuDevice)
  CALL cuSetAxialGcDhat(cuCMFD, cuGcCMFD, cuDevice)
ENDIF

IF (.NOT. cuCntl%lNatural) CALL cuReorderFlux(cuGcCMFD, cuDevice, 1)
CALL cuCopyFlux(cuGcCMFD, cuDevice, 1)

IF (cuCntl%lNatural) THEN
  IF (cuCntl%lShift) THEN
    CALL cuSetNaturalBiCGSystem(cuGcCMFD, cuDevice, l3dim, .TRUE., seigv)
  ELSE
    IF (cuCntl%CMFDSolver .EQ. 1) CALL cuSetNaturalBiCGSystem(cuGcCMFD, cuDevice, l3dim, .FALSE., 0.0)
    IF (cuCntl%CMFDSolver .EQ. 2) CALL cuSetNaturalJacobiSystem(cuGcCMFD, cuDevice, l3dim)
  ENDIF
ELSE
  CALL cuSetRedblackSORSystem(cuGcCMFD, cuDevice, l3dim)
ENDIF

CALL cuSetCMFDSourceOperator(cuGcCMFD, cuDevice, GcGroupInfo%lUpscat)
CALL cuCMFDPsiUpdt(cuGcCMFD, cuDevice)
CALL cuCMFDSrcUpdt(cuGcCMFD, cuDevice, eigv, GcGroupInfo%lUpscat)

outIter = 0; loutConv = FALSE
DO WHILE (.NOT. loutConv)

  CALL cuInnerSolve(cuGcCMFD, cuDevice, tol)
  CALL cuCMFDPsiUpdt(cuGcCMFD, cuDevice)
  CALL cuCMFDEigUpdt(cuGcCMFD, cuDevice, eigv)
  CALL cuCMFDSrcUpdt(cuGcCMFD, cuDevice, eigv, GcGroupInfo%lUpscat)

  outRes = cuCMFDPsiErr(cuGcCMFD, cuDevice)
  IF (outIter .EQ. 0) outRes0 = outRes
  outErr = outRes / outRes0
  outIter = outIter + 1
  loutConv = outErr .LE. outTol
  Keff = 1.0 / (1.0 / eigv + seigv)

  ItrCntl%GcCMFDIt = ItrCntl%GcCMFDIt + 1
  IF (PE%MASTER) WRITE(mesg, '(a9, i9, f22.6, 3x, f10.5, 1p, e15.3)') 'CGOUTER', ItrCntl%GcCMFDIt, Keff, outErr, outRes
  IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)

ENDDO

WRITE(mesg, '(a)') 'Group Reconstruction...'
IF (PE%Master) CALL message(io8, TRUE, TRUE, mesg)

CALL cuCopyFlux(cuGcCMFD, cuDevice, 2)
IF (.NOT. cuCntl%lNatural) CALL cuReorderFlux(cuGcCMFD, cuDevice, 2)
CALL cuGcReconstruction(cuCMFD, cuGcCMFD, cuDevice)

CALL DeallocCMFDVar(cuGcCMFD)

cuCntl%lNatural = lNatural
cuCntl%lDcpl = lDcpl

END SUBROUTINE

#endif
