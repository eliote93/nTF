#include <CUDADEFINES.h>
#include <defines.h>

#ifdef __PGI
MODULE CUDA_Transient

USE CUDA_MASTER
USE CUDA_UTIL
USE CUDA_SYSTEM
USE CUDA_SOLVER
IMPLICIT NONE
REAL(8), ALLOCATABLE, DEVICE :: psidbg(:), phidbg(:,:)
REAL(8), ALLOCATABLE :: TranDtil(:,:,:,:), TranDhat(:,:,:,:)
LOGICAL :: lalloc = .TRUE.
REAL(8), ALLOCATABLE, DEVICE :: shape_prev(:,:), shape_now(:,:)
REAL(8), ALLOCATABLE, DEVICE :: rtilda_prev(:,:), rtilda_now(:,:)
REAL(8) :: amp_prev, amp_now
CONTAINS

SUBROUTINE cuInitTransient(Core, FmInfo, CmInfo, RayInfo, ThInfo, GroupInfo, nTracerCntl, TranInfo, TranCntl, eigv, l3dim)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,     FmInfo_Type,      Pin_Type,     Cell_Type,   &
                             TranInfo_Type,     TranCntl_Type,    CmInfo_Type,  RayInfo_Type,&
                             ThInfo_Type,       GroupInfo_Type
USE CNTL,             ONLY : nTracerCntl_Type
USE PE_Mod,           ONLY : PE
USE TH_MOD,           ONLY : ThVar,             ThOpt, CoolantTHSave, FuelTHSave, fxrtempsave
USE TRAN_MOD,         ONLY : FluxNormTransient, SettimeStep,      InitFmPrecursor, SetSamplingTimeStep
USE MOC_MOD,          ONLY : PsiUpdate
USE MOC_COMMON,       ONLY : setMocPsi
USE BasicOperation,   ONLY : CP_CA,             CP_VA,            MULTI_CA
USE BenchXs,          ONLY : DnpLambdaBen,      ChidBen
USE TranMacXsLib_Mod, ONLY : InitChidLib,       InitLambdaLib,    initchidklib
USE CUDA_AXMOC,       ONLY : cuSetAxPsiNowStep
USE CUDA_PWDIST,      ONLY : cuInithPrec
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(RayInfo_Type) :: RayInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
REAL :: eigv
LOGICAL :: l3dim

TYPE(superPin_Type), POINTER :: superPin(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(dim3) :: Blocks, Threads
REAL(8), ALLOCATABLE :: psi(:, :)
REAL(8), ALLOCATABLE :: pPrec(:, :, :)
REAL(8), ALLOCATABLE, DEVICE :: axPrec(:, :)
REAL, POINTER :: hz(:), hzfm(:)
REAL, ALLOCATABLE :: precsum(:)
INTEGER, POINTER :: PinMap(:), fmRange(:, :)
REAL :: vol, volsum, tmpvolsum
REAL :: psisum, inv_psisum
REAL(GPU_NODAL_PRECISION) :: inveigv_ax
INTEGER(KIND = cuda_stream_kind) :: myStream
INTEGER :: FsrIdxSt, nLocalFsr
INTEGER :: nxy, nzCMFD, myzb, myze, myzbf, myzef, nprec, ng, nzSub, nfsr, nfuelcell, nz
INTEGER :: ixy, ixy_map, iz, izf, ipin, icel, ifsr, iprec, ig, ifxr
INTEGER :: i, j
INTEGER :: n, ierr

nxy = cuGeometry%nxyc
nprec = cuGeometry%nprec
ng = cuGeometry%ng
nzCMFD = cuDevice%nzCMFD
myzb = cuDevice%myzb
myze = cuDevice%myze
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
myStream = cuDevice%myStream

nfsr = Core%nCorefsr

Pin => Core%Pin
Cell => Core%CellInfo
superPin => cuGeometry%superPin
pinMap => cuGeometry%pinMap
fmRange => cuGeometry%fmRange
hz => cuGeometry%hz
hzfm => cuGeometry%hzfm

volsum = 0.
!$OMP PARALLEL PRIVATE(iz, ixy_map)
DO izf = myzbf, myzef
  iz = cuGeometry%planeMap(izf)
  !$OMP DO SCHEDULE(GUIDED) REDUCTION(+:volsum)
  DO ixy = 1, nxy
    ixy_map = pinMap(ixy)
    IF(superPin(ixy_map)%lFuel(iz)) THEN
      volsum = volsum + cuGeometry%PinVolFm(ixy_map, izf)
    END IF
  END DO
  !$OMP END DO
END DO
!$OMP END PARALLEL
#ifdef MPI_ENV
CALL MPI_ALLREDUCE(volsum, tmpvolsum, 1, MPI_DOUBLE, MPI_SUM, MPI_CUDA_COMM, ierr)
volsum = tmpvolsum
#endif
cuTranCMInfo%Inv_FuelVol = 1._8 / volsum

nfuelcell = 0
!$OMP PARALLEL PRIVATE(icel)
DO iz = myzb, myze
  !$OMP DO SCHEDULE(GUIDED) REDUCTION(+:nfuelcell)
  DO ixy = 1, Core%nxy
    icel = Pin(ixy)%Cell(iz)
    IF(Cell(icel)%lFuel) THEN
      nfuelcell = nfuelcell + 1
    END IF
  END DO
  !$OMP END DO
END DO
!$OMP END PARALLEL
TranInfo%nfuelcell = nfuelcell
#ifdef MPI_ENV
CALL MPI_ALLREDUCE(TranInfo%nfuelcell, nfuelcell, 1, MPI_INTEGER, MPI_SUM, MPI_CUDA_COMM, ierr)
TranInfo%nfuelcell = nfuelcell
#endif

IF(.NOT. nTracerCntl%ladjoint) THEN
  ALLOCATE(cuCMFD_adj%phis8(ng, nxy*nzCMFD))
  cuCMFD_adj%phis8 = 1.
END IF
!-------------------------------------------------------------------------------------------------------
!CALL FluxNormTransient(Core, RayInfo, FmInfo, CmInfo, TranInfo, ng, nTracerCntl, PE)
!CALL SetUnitPowerLevel(Core, CmInfo, TranInfo, nTracerCntl, PE, ng)
CALL cuUpdtPowerLevel(TranInfo, .TRUE., cuCntl%lPwDist, nTracerCntl%lDcyHeat)
CALL SetTimeStep(TranCntl)
IF(TranCntl%lNNSampling) CALL SetSamplingTimeStep(TranCntl)
TranInfo%eigv0 = eigv
TranInfo%PowerLevel0 = nTracerCntl%PowerLevel
TranInfo%PowerLevel = nTracerCntl%PowerLevel
TranInfo%PowerLeveld = nTracerCntl%PowerLevel

! Eigenvalue Normalization
FmInfo%psid(:, myzb : myze) = FmInfo%psi(:, myzb : myze)
IF (nTracerCntl%lMacro) THEN
  CALL SetMOCPsi(Core, FmInfo%phis, FmInfo%psi)
ELSE
  CALL PsiUpdate(Core, FmInfo%Fxr, FmInfo%phis, FmInfo%psi, myzb, myze, ng, nTracerCntl%lxslib, GroupInfo)
ENDIF
CALL MULTI_CA(1._8/eigv, FmInfo%Psi(1:nfsr, myzb:myze), nfsr, myze-myzb+1)

CALL CP_VA(FmInfo%TranPhi(1:nFsr, myzb:myze, 1:ng), FmInfo%Phis(1:nFsr, myzb:myze, 1:ng), nFsr, myze-myzb+1, ng)
CALL CP_VA(FmInfo%TranPsi(1:nFsr, myzb:myze), FmInfo%Psi(1:nFsr, myzb:myze), nFsr, myze-myzb+1)
CALL CP_VA(FmInfo%TranPsid(1:nFsr, myzb:myze), FmInfo%Psi(1:nFsr, myzb:myze), nFsr, myze-myzb+1)
CALL CP_VA(FmInfo%TranPower(1:nFsr, myzb:myze), FmInfo%Power(1:nFsr, myzb:myze), nFsr, myze-myzb+1)

!Exponential Transformation
CALL CP_CA(TranInfo%FmExpo(1:nFsr, myzb:myze, 1:ng), 1._8, nFsr, myze - myzb + 1, ng)
CALL CP_CA(TranInfo%FmExpo_Alpha(1:nFsr, myzb:myze, 1:ng), 0._8, nFsr, myze - myzb + 1, ng)

!Init Transient TH
IF(nTracerCntl%lFeedback) THEN
  CALL InitTransientTH(THInfo, ThVar, Core%nxy, Core%nz)
ENDIF

CALL KinParamGen(Core, FmInfo, TranInfo,  ThInfo, GroupInfo, .TRUE., nTracerCntl, PE)
IF(nTRACERCntl%lXsLib) THEN
  CALL InitLambdaLib(TranInfo%Lambda, nTracerCntl%refdcy_del, nTracerCntl%llibdcy_del, nTracerCntl%lfitbeta)
  IF(nTracerCntl%lchidkgen) THEN
    CALL InitChidkLib(Core, FmInfo, GroupInfo, PE, nTracerCntl, TranInfo%chid, TranInfo%Chidk, ng, nprec)
  ELSE
    CALL InitChidLib(Core, FmInfo, GroupInfo, PE, nTracerCntl, TranInfo%Chid, ng)
  END IF
  !IF(nTracerCntl%lchidkgen) THEN
  !  TranInfo%chid = 0.
  !  DO iprec = 1, nprec
  !    DO ig = 1, ng
  !      TranInfo%chid(ig) = TranInfo%chid(ig) + TranInfo%chid(ig, iprec)
  !    END DO
  !  END DO
  !  TranInfo%chid = TranInfo%chid / sum(TranInfo%chid)
  !END IF
ELSE
  !Chid, Lambda
   CALL ChidBen(TranInfo%Chid)
   CALL DnpLambdaBen(TranInfo%Lambda)
ENDIF

DO i = 1, nprec
  TranInfo%InvLambda(i) = 1._8 / TranInfo%Lambda(i)
ENDDO

DO iz = myzb, myze
  DO ifxr = 1, Core%nCoreFxr
    FmInfo%Fxr(ifxr, iz)%imix0 = FmInfo%Fxr(ifxr, iz)%imix
  ENDDO
ENDDO

CALL InitFmPrecursor(Core, FmInfo, TranInfo, nTracerCntl, PE, ng, nprec)
!-------------------------------------------------------------------------------------------------------

n = ng * nxy * nzCMFD
ierr = cublasDcopy_v2(cuDevice%myblasHandle, n, cuCMFD%phis8, 1, cuTranCMInfo%TranPhi, 1)

n = nxy * nzCMFD
ierr = cublasDscal_v2(cuDevice%myblasHandle, n, 1._8/eigv, cuCMFD%psi8, 1)
ierr = cublasDcopy_v2(cuDevice%myblasHandle, n, cuCMFD%psi8, 1, cuTranCMInfo%TranPsi, 1)
ierr = cublasDcopy_v2(cuDevice%myblasHandle, n, cuCMFD%psi8, 1, cuTranCMInfo%TranPsid, 1)

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO iz = myzb, myze
  DO ixy = 1, nxy
    DO ig = 1, ng
      cuTranCMInfo%PhiC(ig, ixy, iz) = cuCMFD%h_phic8(ig, ixy, iz)
      cuTranCMInfo%TranPhiC(ig, ixy, iz) = cuTranCMInfo%PhiC(ig, ixy, iz)
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

ALLOCATE(psi(nxy, myzbf:myzef), pPrec(nxy, myzbf:myzef, nprec))
ierr = cudaMemcpy(psi, cuCMFD%psi8, n, cudaMemcpyDeviceToHost)

ALLOCATE(precsum(nprec))
!!$OMP PARALLEL PRIVATE(ixy_map, precsum, ipin, FsrIdxSt, icel, nLocalFsr, ifsr, vol, psisum, inv_psisum)
!!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO iz = myzb, myze
  DO ixy = 1, nxy
    ixy_map = pinMap(ixy)
    precsum = 0.
    DO j = 1, superPin(ixy_map)%nxy
      ipin = superPin(ixy_map)%pin(j)
      FsrIdxSt = Pin(ipin)%FsrIdxSt
      icel = Pin(ipin)%Cell(iz)
      nLocalFsr = Cell(icel)%nFsr

      DO i = 1, nLocalFsr
        ifsr = FsrIdxSt + i - 1
        vol = Cell(icel)%vol(i) * Core%hz(iz)
        DO iprec = 1, nprec
          precsum(iprec) = precsum(iprec) + FmInfo%Prec(iprec, ifsr, iz) * vol
        END DO
      END DO
    END DO
    psisum = 0.
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      psisum = psisum + psi(ixy, izf)
      DO iprec = 1, nprec
        cuTranCMInfo%Prec(iprec, ixy, izf) = precsum(iprec) * psi(ixy, izf)
      END DO
    END DO
    IF(psisum .GT. 0) THEN
      inv_psisum = 1./psisum
      DO izf = fmRange(iz, 1), fmRange(iz, 2)
        DO iprec = 1, nprec
          cuTranCMInfo%Prec(iprec, ixy, izf) = cuTranCMInfo%Prec(iprec, ixy, izf) * inv_psisum
        END DO
      END DO
    ELSE
      DO izf = fmRange(iz, 1), fmRange(iz, 2)
        DO iprec = 1, nprec
          cuTranCMInfo%Prec(iprec, ixy, izf) = 0.
        END DO
      END DO
    END IF
  END DO
END DO
!!$OMP END DO
!!$OMP END PARALLEL
DEALLOCATE(precsum)
DEALLOCATE(psi)

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO iprec = 1, nprec
  DO izf = myzbf, myzef
    DO ixy = 1, nxy
      pPrec(ixy, izf, iprec) = cuTranCMInfo%Prec(iprec, ixy, izf)
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

n = nxy * nzCMFD * nprec
ierr = cudaMemcpy(cuTranCMInfo%dPrec, pPrec, n, cudaMemcpyHostToDevice)


CALL cuSetTransientConst(TranInfo, TranCntl)
IF(cuCntl%lAxial .AND. l3dim) THEN
  !$OMP PARALLEL PRIVATE(ixy_map)
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
  DO iprec = 1, nprec
    DO izf = myzbf, myzef
      DO ixy = 1, nxy
        ixy_map = pinMap(ixy)
        pPrec(ixy, izf, iprec) = pPrec(ixy, izf, iprec) / cuGeometry%PinVolFm(ixy_map, izf)
      END DO
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
  ALLOCATE(axPrec(nxy * nzCMFD, nprec))
  ierr = cudaMemcpy(axPrec, pPrec, n, cudaMemcpyHostToDevice)
END IF
DEALLOCATE(pPrec)

IF(cuCntl%lPwDist) THEN
  cuPwDist%lTran = .TRUE.
  IF(nTracerCntl%lDcyHeat) THEN
    CALL cuInithPrec(cuPwDist)
  END IF
END IF

cuTranCMInfo%Expo = 1.
cuTranCMInfo%Expo_Alpha = 0.

IF(cuCntl%lAxial .AND. l3dim) THEN
  CALL cuSetAxPsiNowStep(1./eigv)
  nzSub = cuDevice%nzSub
  n = nxy * ng * nzSub
  ierr = cublasCopy(cuDevice%myblasHandle, n, cuAxial%phiCoeff(:,:,:,0), 1, cuAxial%TranPhi, 1)
  n = nxy * nzSub
  ierr = cublasCopy(cuDevice%myblasHandle, n, cuAxial%PsiCoeff(:,:,0), 1, cuAxial%TranPsi, 1)
  ierr = cublasCopy(cuDevice%myblasHandle, n, cuAxial%PsiCoeff(:,:,0), 1, cuAxial%TranPsid, 1)

  Threads = dim3(cuDevice%cuThreadPerBlock, 1, 1)
  Blocks = dim3(nxy/ Threads%x + 1, 1, 1)
  !$ACC HOST_DATA USE_DEVICE(cuDevice)
  CALL cuCopyPrecKernel <<<Blocks, Threads, 0, myStream>>>              &
                        (cuDevice, axPrec, cuAxial%Prec, cuAxial%PsiCoeff)
  !$ACC END HOST_DATA

  DEALLOCATE(axPrec)
  cuAxial%Expo = 1.
  cuAxial%Expo_Alpha = 0.
END IF

cuCntl%PCRtype = TranCntl%PCRtype

IF(TranCntl%ladptt) THEN
  ALLOCATE(shape_prev(ng, nxy*nzCMFD), shape_now(ng, nxy*nzCMFD))
  ALLOCATE(rtilda_prev(ng, nxy*nzCMFD), rtilda_now(ng, nxy*nzCMFD))
  rtilda_now = 0.
  ALLOCATE(CoolantTHSave(Core%nxy), FuelTHSave(Core%nxy))
  nz = Core%nz
  DO ixy = 1, Core%nxy
    ALLOCATE(CoolantTHSave(ixy)%Tcool(1:nz))
    ALLOCATE(CoolantTHSave(ixy)%DenCool(1:nz))
    ALLOCATE(CoolantTHSave(ixy)%hCool(1:nz))
    ALLOCATE(CoolantTHSave(ixy)%rhou(0:nz))
    ALLOCATE(CoolantTHSave(ixy)%rhohu(0:nz))
    ALLOCATE(CoolantTHSave(ixy)%u(0:nz))
    ALLOCATE(CoolantTHSave(ixy)%qeff(1:nz))
    ALLOCATE(FuelTHSave(ixy)%qvol(1:nz))
    ALLOCATE(FuelTHSave(ixy)%Tfuel(1:ThOpt%nrpellet+5, 1:nz))
    ALLOCATE(CoolantTHSave(ixy)%Tcoold(1:nz))
    ALLOCATE(CoolantTHSave(ixy)%DenCoold(1:nz))
    ALLOCATE(CoolantTHSave(ixy)%hCoold(1:nz))
    ALLOCATE(CoolantTHSave(ixy)%rhoud(0:nz))
    ALLOCATE(CoolantTHSave(ixy)%rhohud(0:nz))
    ALLOCATE(CoolantTHSave(ixy)%ud(0:nz))
    ALLOCATE(CoolantTHSave(ixy)%qeffd(1:nz))
    ALLOCATE(FuelTHSave(ixy)%qvold(1:nz))
    ALLOCATE(FuelTHSave(ixy)%Tfueld(1:ThOpt%nrpellet+5, 1:nz))
  END DO
  ALLOCATE(FxrTempSave(Core%nCoreFxr, myzb:myze, 10))
END IF

END SUBROUTINE

SUBROUTINE cuSetTransientConst(TranInfo, TranCntl)
USE TYPEDEF,          ONLY : TranInfo_Type,     TranCntl_Type
USE CUDA_CONST
IMPLICIT NONE
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl

INTEGER :: ig, iprec

nprec = cuGeometry%nprec
theta = TranCntl%theta
DO ig = 1, ng
  chid(ig) = TranInfo%chid(ig)
END DO
IF(TranCntl%lchidk) THEN
  DO iprec = 1, nprec
    DO ig = 1, ng
      chidk(ig, iprec) = TranInfo%chidk(ig, iprec)
    END DO
  END DO
END IF


END SUBROUTINE

SUBROUTINE cuUpdtResSrc(Core, FmInfo, cuCMFD, cuDevice, TranInfo, TranCntl, l3dim)
USE TYPEDEF,          ONLY : CoreInfo_Type,     FmInfo_Type,      TranInfo_Type,      TranCntl_Type
USE CUDAFOR
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(cuCMFD_Type) :: cuCMFD
TYPE(cuDevice_Type) :: cuDevice
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
LOGICAL :: l3dim

TYPE(superPin_Type), POINTER :: superPin(:)
REAL(8), ALLOCATABLE :: ResSrc(:, :, :)
REAL(GPU_NODAL_PRECISION), ALLOCATABLE :: ResSrca(:, :, :)
REAL, POINTER :: chz(:), cPinVol(:, :)
INTEGER, POINTER :: pinMap(:), fmRange(:, :)
TYPE(dim3) :: Blocks, Threads
REAL :: ResSum(cuCMFD%ng)
REAL :: invArea, pinarea, invVol, resVal
INTEGER(KIND = cuda_stream_kind) :: myStream
INTEGER :: bottomRange(2), topRange(2)
INTEGER :: nr, nc, nnz
INTEGER :: ng, nxy, myzbf, myzef, nzCMFD, myzb, myze, ncel, nzSub, nprec
INTEGER :: i, j, ixy, ixy_map, iz, ipin, izf, ig, iprec
INTEGER :: n, ierr

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nprec = cuGeometry%nprec
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
nzCMFD = cuDevice%nzCMFD
myzb = cuDevice%myzb
myze = cuDevice%myze
ncel = cuGeometry%nxyc * nzCMFD
myStream = cuDevice%myStream

Threads = dim3(ng, cuDevice%cuMaxThreadPerBlock/ng, 1)
Blocks = dim3(ncel / Threads%y + 1, 1, 1)

bottomRange = cuCMFD%bottomRange
topRange = cuCMFD%topRange
nc = cuCMFD%M%nc

CALL cuSetNaturalTranBiCGSystem(cuCMFD, cuDevice, TranInfo, TranCntl, l3dim, .FALSE., .FALSE.)
IF(TranCntl%lchidk) THEN
  DO iprec = 1, nprec
    CALL cuSetPrecSrcOperator(iprec, TranInfo%lambda(iprec))
    CALL cuCMFDPrecSrcKUpdt(TranInfo, TranCntl, iprec)
  END DO
  CALL cuComputeDelayedSrc_chidk <<< Blocks, Threads, 0, myStream>>>                             &
    (cuCMFD%src8, cuTranCMInfo%PrecSrcK, nzCMFD)
ELSE
  CALL cuSetPrecSrcOperator()
  CALL cuCMFDPrecSrcUpdt(TranInfo, TranCntl)
  CALL cuComputeDelayedSrc <<< Blocks, Threads, 0, myStream>>>                             &
    (cuCMFD%src8, cuTranCMInfo%PrecSrc, nzCMFD)
END IF
CALL cuMatMul3D(cuCMFD%M, cuCMFD%offDiag8, cuCMFD%phis8, cuTranCMInfo%ResSrc, bottomRange, topRange,                               &
                  cuDevice%mySparseHandle, cuDevice%myStream, MPI_CUDA_COMM)
CALL cuVectorOp('-', nc, cuCMFD%src8, cuTranCMInfo%ResSrc, cuTranCMInfo%ResSrc, cuDevice%myStream)

superPin => cuGeometry%superPin
pinMap => cuGeometry%PinMap
fmRange => cuGeometry%fmRange
chz => Core%hz
cPinVol => Core%PinVol

ALLOCATE(ResSrc(ng, nxy, myzbf:myzef))

n = ng * nxy * nzCMFD
ierr = cudaMemcpy(ResSrc, cuTranCMInfo%ResSrc, n, cudaMemcpyDeviceToHost)

!$OMP PARALLEL PRIVATE(ixy_map, ResSum, invArea, ipin, pinarea)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO iz = myzb, myze
  DO ixy = 1, nxy
    ixy_map = PinMap(ixy)
    ResSum = 0.
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      ResSum = ResSum + ResSrc(:, ixy, izf)
    END DO
    !print*, ressum, ixy, iz
    !STOP
    invArea = 1./superPin(ixy_map)%Area
    DO j = 1, superPin(ixy_map)%nxy
      ipin = superPin(ixy_map)%pin(j)
      pinarea = cPinVol(ipin, iz) / chz(iz)
      FmInfo%ResSrc(ipin, iz, :) = ResSum * pinarea * invArea
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

IF(cuCntl%lAxial .AND. l3dim) THEN
  nzSub = cuDevice%nzSub
  ALLOCATE(ResSrca(nxy, ng, nzSub))
  !$OMP PARALLEL PRIVATE(ixy_map, invVol, ResVal)
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
  DO iz = myzbf, myzef
    DO ipin = 1, nxy
      ixy_map = PinMap(ipin)
      invVol = 1./ cuGeometry%PinVolFm(ixy_map, iz)
      DO ig = 1, ng
        ResVal = ResSrc(ig, ipin, iz) * invVol
        DO izf = cuDevice%fmSubrange(iz, 1), cuDevice%fmSubRange(iz, 2)
          ResSrca(ipin, ig, izf) = ResVal
        END DO
      END DO
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
  n = nxy * ng * nzSub
  ierr = cudaMemcpy(cuAxial%ResSrc, ResSrca, n, cudaMemcpyHostToDevice)
  DEALLOCATE(ResSrca)
END IF

DEALLOCATE(ResSrc)
NULLIFY(superPin, pinMap, fmRange, chz, cPinVol)
CALL destroyCsr(cuCMFD%M)

END SUBROUTINE

SUBROUTINE cuCMFDPrecSrcUpdt(TranInfo, TranCntl)
USE TYPEDEF,          ONLY : TranInfo_Type,     TranCntl_Type
IMPLICIT NONE
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl

REAL(8), ALLOCATABLE :: kappa(:)
REAL :: delt
INTEGER :: nowstep
INTEGER :: nprec, nxy, myzbf, myzef, nzCMFD
INTEGER :: iprec, ixy, izf
INTEGER :: n, ierr

nowstep = TranCntl%nowstep
delt = TranCntl%delt(nowstep)

nprec = cuGeometry%nprec
nxy = cuGeometry%nxyc
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
nzCMFD = cuDevice%nzCMFD

ALLOCATE(kappa(nprec))
DO iprec = 1, nprec
  kappa(iprec) = exp(-delt * TranInfo%lambda(iprec))
END DO

n = nxy * nzCMFD
CALL cuInitArray(n, cuTranCMInfo%PrecSrc, cuDevice%myStream)
DO iprec = 1, nprec
  CALL cuSumDecayedPrecursor(cuTranCMInfo%dPrec(:,iprec), TranInfo%lambda(iprec) * kappa(iprec), n)
END DO

CALL cuVectorOp('*+', n, cuTranCMInfo%Omegalm, cuTranCMInfo%TranPsid, &
                cuTranCMInfo%PrecSrc, cuDevice%myStream)
CALL cuVectorOp('*+', n, cuTranCMInfo%Omegal0, cuTranCMInfo%TranPsi, &
                cuTranCMInfo%PrecSrc, cuDevice%myStream)

DEALLOCATE(kappa)

END SUBROUTINE

SUBROUTINE cuCMFDPrecSrcKUpdt(TranInfo, TranCntl, iprec)
USE TYPEDEF,          ONLY : TranInfo_Type,     TranCntl_Type
IMPLICIT NONE
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl

REAL(8) :: kappa
REAL :: delt
INTEGER :: nowstep
INTEGER :: nprec, nxy, myzbf, myzef, nzCMFD
INTEGER :: iprec, ixy, izf
INTEGER :: n, ierr

REAL, allocatable :: dummy(:)

nowstep = TranCntl%nowstep
delt = TranCntl%delt(nowstep)

nprec = cuGeometry%nprec
nxy = cuGeometry%nxyc
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
nzCMFD = cuDevice%nzCMFD

kappa = exp(-delt * TranInfo%lambda(iprec))

n = nxy * nzCMFD
CALL cuInitArray(n, cuTranCMInfo%PrecSrcK(:,iprec), cuDevice%myStream)
ierr = cublasDaxpy_v2(cuDevice%myblasHandle, n, TranInfo%lambda(iprec) * kappa, cuTranCMInfo%dPrec(:,iprec), 1, cuTranCMInfo%PrecSrcK(:,iprec), 1)

!ALLOCATE(dummy(n))
!dummy = cuTranCMINfo%precsrck(:,iprec)
!print*, dummy
CALL cuVectorOp('*+', n, cuTranCMInfo%Omegalm, cuTranCMInfo%TranPsid, &
                cuTranCMInfo%PrecSrcK(:,iprec) , cuDevice%myStream)
CALL cuVectorOp('*+', n, cuTranCMInfo%Omegal0, cuTranCMInfo%TranPsi, &
                cuTranCMInfo%PrecSrcK(:,iprec), cuDevice%myStream)


END SUBROUTINE

SUBROUTINE cuSumDecayedPrecursor(Prec, lkappa, n)
IMPLICIT NONE
REAL(8), DEVICE :: Prec(*)
REAL(8) :: lkappa
INTEGER :: n

INTEGER :: ierr

ierr = cublasDaxpy_v2(cuDevice%myblasHandle, n, lkappa, Prec, 1, cuTranCMInfo%PrecSrc, 1)

END SUBROUTINE

SUBROUTINE cuSetPrecSrcOperator(iprec, lambda)
IMPLICIT NONE
INTEGER, OPTIONAL :: iprec
REAL, OPTIONAl :: lambda

REAL(8), ALLOCATABLE :: Omegalm(:, :), Omegal0(:, :)
INTEGER, POINTER :: planeMap(:), PinMap(:)
INTEGER :: nxy, nzCMFD, myzbf, myzef
INTEGER :: ipin,  ipin_map, izf, iz
INTEGER :: n, ierr
LOGICAL :: lchidk

lchidk = PRESENT(iprec)

nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef

pinMap => cuGeometry%pinMap
planeMap => cuCMFD%planeMap

ALLOCATE(Omegalm(nxy, myzbf:myzef))
ALLOCATE(Omegal0(nxy, myzbf:myzef))

!$OMP PARALLEL PRIVATE(ipin_map)
DO izf = myzbf, myzef
  iz = planeMap(izf)
  !$OMP DO SCHEDULE(GUIDED)
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    IF(lchidk) THEN
      Omegalm(ipin, izf) = cuTranCMInfo%CellOmegam(iprec, ipin_map, iz) * lambda
    ELSE
      Omegalm(ipin, izf) = cuTranCMInfo%CellOmegam(0, ipin_map, iz)
    END IF
  END DO
  !$OMP END DO
END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(ipin_map)
DO izf = myzbf, myzef
  iz = planeMap(izf)
  !$OMP DO SCHEDULE(GUIDED)
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    IF(lchidk) THEN
      Omegal0(ipin, izf) = cuTranCMInfo%Cellomega0(iprec, ipin_map, iz) * lambda
    ELSE
      Omegal0(ipin, izf) = cuTranCMInfo%Cellomega0(0, ipin_map, iz)
    END IF
  END DO
  !$OMP END DO
END DO
!$OMP END PARALLEL

n = nxy * nzCMFD
ierr = cudaMemcpy(cuTranCMInfo%Omegalm, Omegalm, n, cudaMemcpyHostToDevice)
ierr = cudaMemcpy(cuTranCMInfo%Omegal0, Omegal0, n, cudaMemcpyHostToDevice)

NULLIFY(PinMap, PlaneMap)
DEALLOCATE(Omegalm, Omegal0)

END SUBROUTINE

SUBROUTINE cuSetInvVelo()
IMPLICIT NONE

REAL(8), ALLOCATABLE :: pinvVel(:, :, :)
REAL, POINTER :: PinVolFm(:,:)
INTEGER, POINTER :: PinMap(:), PlaneMap(:)
REAL :: vol
INTEGER :: ng, nxy, nzCMFD, myzbf, myzef
INTEGER :: ig, ipin, ipin_map, iz, izf
INTEGER :: n, ierr

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
nzCMFD = cuDevice%nzCMFD

pinMap => cuGeometry%pinMap
planeMap => cuCMFD%planeMap
PinVolFm => cuGeometry%PinVolFm

ALLOCATE(pInvVel(ng, nxy, myzbf:myzef))

!$OMP PARALLEL PRIVATE(ipin_map, vol)
DO izf = myzbf, myzef
  iz = PlaneMap(izf)
  !$OMP DO SCHEDULE(GUIDED)
  DO ipin = 1, nxy
    ipin_map = PinMap(ipin)
    vol = PinVolFm(ipin_map, izf)
    DO ig = 1, ng
      pInvVel(ig, ipin, izf) = vol /cuCMFD%PinXS(ipin_map, iz)%velo(ig)
    END DO
  END DO
  !$OMP END DO
END DO
!$OMP END PARALLEL
NULLIFY(PinMap, PlaneMap, PinVolFm)

n = nxy * nzCMFD * ng
ierr = cudaMemcpy(cuTranCMInfo%VolInvVel, pInvVel, n, cudaMemcpyHostToDevice)
DEALLOCATE(pInvVel)

END SUBROUTINE

SUBROUTINE cuTranCMFDSrcUpdt(TranCntl)
USE TYPEDEF,          ONLY : TranCntl_Type
USE CUDAFOR
IMPLICIT NONE
TYPE(TranCntl_Type) :: TranCntl

TYPE(dim3) :: Blocks, Threads
REAL(8) :: delt
INTEGER(KIND = cuda_stream_kind) :: myStream
INTEGER :: ng, nzCMFD, ncel
INTEGER :: ierr, isync

ng = cuGeometry%ng
myStream = cuDevice%myStream
nzCMFD = cuDevice%nzCMFD
ncel = cuGeometry%nxyc * nzCMFD

CALL cuSetInvVelo()

delt = TranCntl%delt(TranCntl%nowstep)

!Threads = dim3(ng, cuDevice%cuThreadPerBlock/ng, 1)
Threads = dim3(ng, 2, 1)
Blocks = dim3(ncel / Threads%y + 1, 1, 1)

IF(TranCntl%lchidk) THEN
  CALL cuComputeTranCMFDSrc_chidk <<< Blocks, Threads, 0, myStream>>>                             &
    (cuCMFD%src8, cuTranCMInfo%ResSrc, cuTranCMInfo%PrecSrcK, cuTranCMInfo%TranPhi,       &
    cuTranCMInfo%VolInvVel, cuTranCMInfo%expo_alpha, cuTranCMInfo%expo, delt, nzCMFD)
ELSE
  CALL cuComputeTranCMFDSrc <<< Blocks, Threads, 0, myStream>>> &
       (cuCMFD%src8, cuTranCMInfo%ResSrc, cuTranCMInfo%PrecSrc, cuTranCMInfo%TranPhi, &
       cuTranCMInfo%VolInvVel, cuTranCMInfo%expo_alpha, cuTranCMInfo%expo, delt, nzCMFD)
END IF
isync = cudaDeviceSynchronize()
ierr = cudaGetLastError()
if (ierr.NE.0) print*, cudaGetErrorString(ierr)

END SUBROUTINE

SUBROUTINE cuExpTrsfUpdt(Core, TranInfo, TranCntl, l3dim)
USE TYPEDEF,          ONLY : CoreInfo_Type,    TranInfo_Type,     TranCntl_Type,     Pin_Type,  &
                             Cell_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
LOGICAL :: l3dim

TYPE(superPin_Type), POINTER :: superPin(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
REAL(8), ALLOCATABLE :: Expo(:, :, :), Expo_alpha(:, :, :)
REAL(GPU_NODAL_PRECISION), ALLOCATABLE :: axExpo(:, :, :), axExpo_alpha(:, :, :)
REAL, POINTER :: FmExpo_Alpha(:, :, :), FmExpo(:, :, :)
INTEGER, POINTER :: PinMap(:)
INTEGER, POINTER :: FmRange(:, :)
REAL :: delt, deltn
REAL :: ratio, alphaval, expoval
INTEGER :: nowstep
INTEGER :: ng, nxy, myzb, myze, myzbf, myzef, nzCMFD
INTEGER :: FsrIdxSt, nLocalFsr
INTEGER :: iz, izf, ig, ixy, ixy_map, i, j, ipin
INTEGER :: icel, ifsr
INTEGER :: n, ierr

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
myzb = cuDevice%myzb
myze = cuDevice%myze
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
nzCMFD = cuDevice%nzCMFD

superPin => cuGeometry%superPin
Pin => Core%Pin
Cell => Core%CellInfo
PinMap => cuGeometry%PinMap
fmRange => cuGeometry%fmRange

FmExpo_Alpha => TranInfo%FmExpo_Alpha; FmExpo => TranInfO%FmExpo

nowstep = TranCntl%nowstep
delt = TranCntl%delt(nowstep)
deltn = delt

ALLOCATE(Expo(ng, nxy, myzbf:myzef), Expo_alpha(ng, nxy, myzbf:myzef))
IF(cuCntl%lAxial .AND. l3dim) THEN
  ALLOCATE(axExpo(nxy, ng, myzb:myze), axExpo_Alpha(nxy, ng, myzb:myze))
END IF
!$OMP PARALLEL PRIVATE(ixy_map, ratio, alphaval, expoval, ipin, FsrIdxSt, icel, nLocalFsr, ifsr)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO iz = myzb, myze
  DO ixy = 1, nxy
    ixy_map = PinMap(ixy)
    DO ig = 1, ng
      IF(cuCMFD%h_phic8(ig, ixy, iz) .GT. 0 .AND. cuTranCMInfo%TranPhiC(ig, ixy, iz) .GT. 0) THEN
        ratio = cuCMFD%h_phic8(ig, ixy, iz) / cuTranCMInfo%TranPhiC(ig, ixy, iz)
        !ratio = min(ratio, 10.0_8)
        alphaval = log(abs(ratio)) / delt
        expoval = exp(alphaval * deltn)
      ELSE
        alphaval = 0.
        expoval = 1.
      END IF
      DO izf = fmRange(iz,1), fmRange(iz,2)
        Expo_alpha(ig, ixy, izf) = alphaval
        Expo(ig, ixy, izf) = expoval
      END DO
      DO j = 1, superPin(ixy_map)%nxy
        ipin = superPin(ixy_map)%pin(j)
        FsrIdxSt = Pin(ipin)%FsrIdxSt
        icel = Pin(ipin)%Cell(iz)
        nLocalFsr = Cell(icel)%nFsr
        DO i = 1, nLocalFsr
          ifsr = FsrIdxSt + i - 1
          fmExpo_Alpha(ifsr, iz, ig) = alphaval
          fmExpo(ifsr, iz, ig) = expoval
        END DO
      END DO
      IF(cuCntl%lAxial .AND. l3dim) THEN
        axExpo(ixy, ig, iz) = expoval
        axExpo_Alpha(ixy, ig, iz) = alphaval
      END IF
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

n = nxy * nzCMFD * ng
ierr = cudaMemcpy(cuTranCMInfo%Expo, Expo, n, cudaMemcpyHostToDevice)
ierr = cudaMemcpy(cuTranCMInfo%Expo_alpha, Expo_alpha, n, cudaMemcpyHostToDevice)
IF(cuCntl%lAxial .AND. l3dim) THEN
  n = nxy * ng * (myze - myzb + 1)
  ierr = cudaMemcpy(cuAxial%Expo, axExpo, n, cudaMemcpyHostToDevice)
  ierr = cudaMemcpy(cuAxial%Expo_alpha, axExpo_alpha, n, cudaMemcpyHostToDevice)
END IF

NULLIFY(superPin, Pin, Cell, PinMap, fmRange)
NULLIFY(fmExpo_Alpha, fmExpo)
DEALLOCATE(Expo, Expo_alpha)
IF(cuCntl%lAxial .AND. l3dim) THEN
  DEALLOCATE(axExpo, axExpo_Alpha)
END IF

END SUBROUTINE

SUBROUTINE cuTranReactivityUpdt(GroupInfo, TranInfo, l3dim, lFirst)
USE TYPEDEF,          ONLY : GroupInfo_Type,      TranInfo_Type
USE CUDA_CMFD
USE CUDA_INIT,        ONLY : DeallocCMFDVar
IMPLICIT NONE
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranInfo_Type) :: TranInfo
LOGICAL :: l3dim, lFirst

TYPE(dim3) :: Blocks, Threads
REAL(8), ALLOCATABLE, DEVICE :: migration(:), InvVelo(:), Precursor(:,:), One(:)
REAL(8), ALLOCATABLE :: recv(:,:,:), psi(:,:), rInvVelo(:,:,:)
INTEGER, POINTER :: PlaneMap(:), PinMap(:)
REAL :: psisum, resphi, betaavg, lifetime
REAL, POINTER :: precsum(:)
REAL :: Buf0(4), Buf(4)
REAL :: locpsi
INTEGER :: bottomRange(2), topRange(2)
INTEGER :: n, nxy, nzCMFD, ng, ncel
INTEGER :: myzbf, myzef
INTEGER :: ixy, izf, ig, iz, ixy_map, iprec
INTEGER :: ierr

REAL :: invvelsum(8), phivolsum(8)
REAL :: invvelsum0(8), phivolsum0(8)
integer :: gidx

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
ncel = nxy * nzCMFD
n = ng * nxy * nzCMFD
bottomRange = cuCMFD%bottomRange
topRange = cuCMFD%topRange
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
!eigv0 = TranInfo%eigv0
!inv_eigv0 = 1./eigv0

ALLOCATE(migration(n), invVelo(n))
CALL cuSetNaturalBiCGSystem(cuCMFD, cuDevice, l3dim, .FALSE., 0.0)
CALL cuSetInvVelo()
CALL cuMatMul3D(cuCMFD%M, cuCMFD%offDiag8, cuCMFD%phis8, migration, bottomRange, topRange,                               &
                  cuDevice%mySparseHandle, cuDevice%myStream, MPI_CUDA_COMM)
CALL cuVectorOp('*', n, cuCMFD%phis8, cuTranCMInfo%VolInvVel, InvVelo, cuDevice%myStream)

ALLOCATE(recv(ng, nxy, myzbf:myzef), rInvVelo(ng, nxy, myzbf:myzef))
ALLOCATE(psi(nxy, myzbf:myzef))

ierr = cudaMemcpy(recv, migration, n, cudaMemcpyDeviceToHost)
ierr = cudaMemcpy(psi, cuCMFD%psi8, nxy*nzCMFD, cudaMemcpyDeviceToHost)
ierr = cudaMemcpy(rInvVelo, InvVelo, n, cudaMemcpyDeviceToHost)

DEALLOCATE(migration, invVelo)

PlaneMap => cuCMFD%PlaneMap
PinMap => cuGeometry%PinMap
psisum = 0
resphi = 0
betaavg = 0
lifetime = 0
invvelsum = 0.
phivolsum = 0.
!print*, GroupInfo%InvGcStruct(1:47)
!!$OMP PARALLEL PRIVATE(iz, ixy_map, locpsi)
DO izf = myzbf, myzef
  iz = PlaneMap(izf)
!  !$OMP DO SCHEDULE(GUIDED) REDUCTION(+:psisum, resphi, betaavg, lifetime, invvelsum, phivolsum)
  DO ixy = 1, nxy
    ixy_map = PinMap(ixy)
    DO ig = 1, ng
      resphi = resphi - recv(ig, ixy, izf) * psi(ixy, izf) / cuGeometry%pinvolfm(cuGeometry%pinmap(ixy), izf)
      lifetime = lifetime + rInvVelo(ig, ixy, izf) * psi(ixy, izf)/ cuGeometry%pinvolfm(cuGeometry%pinmap(ixy), izf)
      gidx = GroupInfo%INvGCStruct(ig)
      invvelsum(gidx) = invvelsum(gidx) + rInvVelo(ig, ixy, izf)
      phivolsum(gidx) = PhiVolSum(gidx) + cuCMFD%h_phis8(ig, ixy, izf) * cuGeometry%pinvolfm(cuGeometry%pinmap(ixy), izf)
    END DO
    locpsi = psi(ixy,izf)
    resphi = resphi + locpsi * locpsi/ cuGeometry%pinvolfm(cuGeometry%pinmap(ixy), izf)
    psisum = psisum + locpsi * locpsi/ cuGeometry%pinvolfm(cuGeometry%pinmap(ixy), izf)
    betaavg = betaavg + cuCMFD%PinXS(ixy_map, iz)%betat * locpsi * locpsi/ cuGeometry%pinvolfm(cuGeometry%pinmap(ixy), izf)
    !write(905,'(i3, es14.6)'), ixy, cuCMFD%PinXS(ixy_map, iz)%betat
  END DO
!  !$OMP END DO
END DO
!!$OMP END PARALLEL

#ifdef MPI_ENV
Buf0 =(/resphi, lifetime, betaavg, psisum/)
CALL MPI_ALLREDUCE(Buf0, Buf, 4, MPI_DOUBLE, MPI_SUM, MPI_CUDA_COMM, ierr)
resphi = Buf(1); lifetime = Buf(2); betaavg = Buf(3); psisum = Buf(4)
CALL MPI_ALLREDUCE(invvelsum, invvelsum0, 8, MPI_DOUBLE, MPI_SUM, MPI_CUDA_COMM, ierr)
CALL MPI_ALLREDUCE(phivolsum, phivolsum0, 8, MPI_DOUBLE, MPI_SUM, MPI_CUDA_COMM, ierr)
#endif

PRINT'(2es14.6)', invvelsum0(1) / phivolsum0(1) , invvelsum0(2) / phivolsum0(2)
PRINT'(2es14.6)', invvelsum0(3) / phivolsum0(3) , invvelsum0(4) / phivolsum0(4)
PRINT'(2es14.6)', invvelsum0(5) / phivolsum0(5) , invvelsum0(6) / phivolsum0(6)
PRINT'(2es14.6)', invvelsum0(7) / phivolsum0(7) , invvelsum0(8) / phivolsum0(8)
!stop

resphi = resphi / psisum
betaavg = betaavg / psisum
TranInfo%lifetime = lifetime / psisum
TranInfo%CoreBeta = betaavg
TranInfo%delRho = resphi * 1.e+5
TranInfo%Reactivity = resphi / betaavg
TranInfo%TranEig = 1./ (1. - resphi)
TranInfo%Factor_F = psisum

IF(lFirst) THEN
  TranInfo%Inv_Factor_F0 = 1./ psisum
  TranInfo%Inv_Factor_K0 = 1./ lifetime
  TranInfo%Inv_lifetime0 = 1./ TranInfo%lifetime
  TranInfo%Amp = 1.
END IF

CALL DeallocCMFDVar(cuCMFD)
DEALLOCATE(psi, recv, rInvVelo)
NULLIFY(PlaneMap)

END SUBROUTINE

SUBROUTINE cuTranReactivityUpdt_Adj(GroupInfo, TranInfo, TranCntl, l3dim, lFirst)
USE TYPEDEF,          ONLY : GroupInfo_Type, TranInfo_Type, TranCntl_Type
USE CUDA_CMFD
USE CUDA_INIT,        ONLY : DeallocCMFDVar
IMPLICIT NONE
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranInfo_Type) ::  TranInfo
TYPE(TranCntl_Type) :: TranCntl
LOGICAL :: l3dim, lFirst

TYPE(dim3) :: Blocks, Threads
REAL(8), ALLOCATABLE, DEVICE :: InvVel(:), FisSrc(:), Migration(:), Residual(:), BetaWeight(:, :)
REAL(8), ALLOCATABLE, DEVICE :: Beta(:)
REAL(8), ALLOCATABLE, DEVICE :: BetakWeight(:,:,:), Betak(:,:)
REAL(8), POINTER, DEVICE :: csrVal(:)
INTEGER, POINTER, DEVICE :: csrRowPtr(:), csrColIdx(:)
REAL, POINTER :: BetaHost(:, :), BetakHost(:,:,:)
REAL :: factor_F, avg_res, avg_beta, avg_invVel
INTEGER, POINTER :: PlaneMap(:), PinMap(:)
INTEGER :: bottomRange(2), topRange(2)
INTEGER :: n, ng, nxy, nzCMFD, myzbf, myzef, ncel, nprec
INTEGER :: ixy, ixy_map, iz, izf, iprec
INTEGER :: nr, nc, nnz
INTEGER :: ierr

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
ncel = nxy * nzCMFD
nprec = TranInfo%nprec
n = ng * nxy * nzCMFD
bottomRange = cuCMFD%bottomRange
topRange = cuCMFD%topRange
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef

PlaneMap => cuCMFD%PlaneMap
PinMap => cuGeometry%PinMap

ALLOCATE(InvVel(n), FisSrc(n), Migration(n), Residual(n))
ALLOCATE(Beta(nxy * nzCMFD), BetaHost(nxy, myzbf:myzef), BetaWeight(ng, nxy*nzCMFD))
ALLOCATE(Betak(nxy * nzCMFD, nprec), BetakHost(nxy, myzbf:myzef, nprec), BetakWeight(ng, nxy*nzCMFD, nprec))

CALL cuSetNaturalBiCGSystem(cuCMFD, cuDevice, l3dim, .FALSE., 0.0)
CALL cuSetCMFDSourceOperator(cuCMFD, cuDevice, GroupInfo%lUpscat)
CALL cuSetInvVelo()

!$OMP PARALLEL PRIVATE(iz, ixy_map)
DO izf = myzbf, myzef
  iz = PlaneMap(izf)
  !$OMP DO SCHEDULE(GUIDED)
  DO ixy = 1, nxy
    ixy_map = PinMap(ixy)
    BetaHost(ixy, izf) = cuCMFD%PinXS(ixy_map, iz)%betat
  END DO
  !$OMP END DO
END DO
!$OMP END PARALLEL
ierr = cudaMemcpy(Beta, BetaHost, nxy*nzCMFD, cudaMemcpyHostToDevice)
DEALLOCATE(BetaHost)

!$OMP PARALLEL PRIVATE(iz, ixy_map)
DO iprec = 1, nprec
  DO izf = myzbf, myzef
    iz = PlaneMap(izf)
    !$OMP DO SCHEDULE(GUIDED)
    DO ixy = 1, nxy
      ixy_map = PinMap(ixy)
      BetakHost(ixy, izf, iprec) = cuCMFD%PinXS(ixy_map, iz)%beta(iprec)
    END DO
    !$OMP END DO
  END DO
END DO
!$OMP END PARALLEL
ierr = cudaMemcpy(Betak, BetakHost, nxy*nzCMFD*nprec, cudaMemcpyHostToDevice)
DEALLOCATE(BetakHost)

Threads = dim3(ng, cuDevice%cuMaxThreadPerBlock/ng, 1)
Blocks = dim3(ncel / Threads%y + 1, 1, 1)

CALL cuComputeBetaAverage <<< Blocks, Threads, 0, cuDevice%myStream>>>                             &
     (BetaWeight, Beta, cuCMFD%psi8, nzCMFD)

DO iprec = 1, nprec
  IF(TranCntl%lchidk) THEN
    CALL cuComputeBetaAverage_iprec <<< Blocks, Threads, 0, cuDevice%myStream>>>                             &
      (BetakWeight(:,:,iprec), Betak(:,iprec), cuCMFD%psi8, nzCMFD, iprec)
  ELSE
    CALL cuComputeBetaAverage <<< Blocks, Threads, 0, cuDevice%myStream>>>                             &
      (BetakWeight(:,:,iprec), Betak(:,iprec), cuCMFD%psi8, nzCMFD)
  END IF
END DO

CALL cuVectorOp('*', n, cuCMFD%phis8, cuTranCMInfo%VolInvVel, InvVel, cuDevice%myStream)
IF(cuDevice%lFuel) THEN
  csrVal => cuCMFD%Chi%d_csrVal
  csrRowPtr => cuCMFD%Chi%d_csrRowPtr
  csrColIdx => cuCMFD%Chi%d_csrColIdx
  nr = cuCMFD%Chi%nr
  nc = cuCMFD%Chi%nc
  nnz = cuCMFD%Chi%nnz
  ierr = cusparseDcsrmv(cuDevice%mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, 1.0_8,                &
    cuCMFD%Chi%descr, csrVal, csrRowPtr, csrColIdx, cuCMFD%psi8, 0.0_8, FisSrc)
ELSE
  FisSrc = 0.
END IF
CALL cuMatMul3D(cuCMFD%M, cuCMFD%offDiag8, cuCMFD%phis8, Migration, bottomRange, topRange,                               &
                  cuDevice%mySparseHandle, cuDevice%myStream, MPI_CUDA_COMM)
CALL cuVectorOp('-', n, FisSrc, Migration, Residual, cuDevice%myStream)

factor_F = dotMulti(cuCMFD_adj%phis8, FisSrc, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
avg_invVel = dotMulti(cuCMFD_adj%phis8, invVel, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
avg_res = dotMulti(cuCMFD_adj%phis8, Residual, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
avg_beta = dotMulti(cuCMFD_adj%phis8, BetaWeight, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)

DO iprec = 1, nprec
  TranInfo%corebetak(iprec) = dotMulti(cuCMFD_adj%phis8, BetakWeight(:,:,iprec), n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
  TranInfo%corebetak(iprec) = TranInfo%corebetak(iprec) / factor_F
END DO

avg_res = avg_res / factor_F
TranInfo%lifetime_Dynamic = avg_invVel / factor_F
TranInfo%delRho_Dynamic = avg_res * 1.e+5
TranInfo%TranEig_Dynamic = 1./(1. - avg_res)
TranInfo%CoreBeta_Dynamic = avg_beta / factor_F
TranInfo%Reactivity_Dynamic = avg_res / TranInfo%CoreBeta_Dynamic
TranInfo%Factor_F = factor_F
IF(TranCntl%lchidk) TranInfo%CoreBeta_Dynamic = SUM(TranInfo%corebetak(1:6))

IF(lFirst) THEN
  TranInfo%Inv_Factor_F0 = 1./ factor_F
  TranInfo%Inv_Factor_K0 = 1./ avg_invVel
  TranInfo%Inv_lifetime0 = 1./ TranInfo%lifetime_Dynamic
  TranInfo%Amp = 1.
END IF

CALL DeallocCMFDVar(cuCMFD)
DEALLOCATE(Beta, BetaWeight)
DEALLOCATE(Betak, BetakWeight)
DEALLOCATE(InvVel, FisSrc, Residual, Migration)

! PRINT '(6es14.6)', TranInfo%corebetak
! PRINT '(3es14.6)', SUM(TranInfo%corebetak(1:6)), TranInfo%coreBeta_Dynamic

END SUBROUTINE

SUBROUTINE cuSavePKEParameters(TranInfo, Trancntl, lthstep)
USE TYPEDEF,        ONLY : TranInfo_Type,   TranCntl_Type
USE PARAM
USE PE_Mod,           ONLY : PE
USE files,            ONLY : io8
USE ioutil,           ONLY : message,           terminate
IMPLICIT NONE
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
LOGICAL, OPTIONAL :: lthstep

TYPE(dim3) :: Blocks, Threads
REAL(8), ALLOCATABLE, DEVICE :: Precursor(:, :)
INTEGER :: ng, nxy, ncel, nzCMFD, n
INTEGER :: iprec

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
ncel = nxy * nzCMFD
n = ng * nxy * nzCMFD

TranInfo%Prev_delrho = TranInfo%delrho_Dynamic
TranInfo%Prev_corebetat = TranInfo%corebeta_Dynamic
TranInfo%Prev_lifetime = TranInfo%lifetime_Dynamic
TranInfo%Prev_factor_F = TranInfo%factor_F

DO iprec = 1, TranInfo%nprec
  TranInfo%Prev_corebeta(iprec) = TranInfo%coreBetak(iprec)
END DO

IF(present(lthstep)) THEN
  DO iprec = 1, TranInfo%nprec
    TranInfo%Prev_corePrec(iprec) = TranInfo%Prev_corePrec(iprec) * TranInfo%precRatio(iprec)
  END DO
ELSE
  Threads = dim3(ng, cuDevice%cuMaxThreadPerBlock/ng, 1)
  Blocks = dim3(ncel / Threads%y + 1, 1, 1)
  ALLOCATE(Precursor(ng, nxy* nzCMFD))
  DO iprec = 1, TranInfo%nprec
    IF(TranCntl%lchidk) THEN
      CALL cuComputeDelayedSrc_iprec <<< Blocks, Threads, 0, cuDevice%myStream>>>                             &
        (Precursor, cuTranCMInfo%dPrec(:, iprec), nzCMFD, iprec)
    ELSE
      CALL cuComputeDelayedSrc <<< Blocks, Threads, 0, cuDevice%myStream>>>                             &
        (Precursor, cuTranCMInfo%dPrec(:, iprec), nzCMFD)
    END IF
    TranInfo%Prev_corePrec(iprec) = dotMulti(cuCMFD_adj%phis8, Precursor, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
  END DO
  DEALLOCATE(Precursor)
  DO iprec = 1, TranInfo%nprec
    TranInfo%Prev_corePrec(iprec) = TranInfo%Prev_corePrec(iprec) * TranInfo%Inv_Factor_F0
  END DO
  WRITE(mesg, '(a10,x, 6es13.5)'), 'PsiAvg:', TranInfo%Factor_F, 1./TranInfo%Inv_Factor_F0
  IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
  WRITE(mesg, '(a10,x, 3es13.5)'), 'CorePrec:', TranInfo%Prev_CorePrec(1:3)
  IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
  WRITE(mesg, '(11x, 3es13.5)'), TranInfo%Prev_CorePrec(4:6)
  IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
END IF


END SUBROUTINE

SUBROUTINE cuUpdtPowerLevel(TranInfo, lFirst, lPwDist, lfpowersave)
USE TYPEDEF,          ONLY : TranInfo_Type
USE CNTL,             ONLY : nTracerCntl
USE CUDA_PWDIST,      ONLY : cuUpdtPinPw
IMPLICIT NONE
TYPE(TranInfo_Type) :: TranInfo
LOGICAL :: lFirst
LOGICAL, OPTIONAL :: lPwDist, lfpowersave

REAL, ALLOCATABLE, DEVICE :: XSkf(:,:)
REAL, POINTER :: XSkfHost(:,:,:)
REAL, POINTER :: PinVolFm(:,:)
REAL :: locVol
REAL :: PwSum, Plevel, AvgPw
INTEGER, POINTER :: PlaneMap(:), PinMap(:)
INTEGER :: myzbf, myzef, nxy, ng, nzCMFD
INTEGER :: ig, ixy, izf
INTEGER :: iz, ixy_map
INTEGER :: ierr, n

IF(PRESENT(lPwDist) .AND. lPwDist) THEN
  CALL cuUpdtPinPw(cuPwDist, nTracerCntl%lDcyHeat .AND. (.NOT. lFirst), lfpowersave, AvgPw)
  IF(lFirst) THEN
    TranInfo%UnitPowerLevel0 = nTracerCntl%PowerLevel / AvgPw
    TranInfo%PwSum0 = PwSum
  END IF
  Plevel = AvgPw * TranInfo%UnitPowerLevel0
  TranInfo%PowerLevel = Plevel
  RETURN
END IF

nxy = cuGeometry%nxyc
ng = cuGeometry%ng
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
nzCMFD = cuDevice%nzCMFD

PlaneMap => cuCMFD%PlaneMap
PinMap => cuGeometry%PinMap
PinVolFm => cuGeometry%PinVolFm

ALLOCATE(XSkfHost(ng, nxy, myzbf:myzef))
!$OMP PARALLEL PRIVATE(iz, ixy_map, locVol)
DO izf = myzbf, myzef
  iz = PlaneMap(izf)
 !$OMP DO SCHEDULE(GUIDED)
  DO ixy = 1, nxy
    ixy_map = PinMap(ixy)
    locVol = PinVolFm(ixy_map, izf)
    DO ig = 1, ng
      XSkfHost(ig, ixy, izf) = cuCMFD%PinXS(ixy_map, iz)%xskf(ig) * locVol
    END DO
  END DO
 !$OMP END DO
END DO
!$OMP END PARALLEL

ALLOCATE(XSkf(ng, nxy*nzCMFD))
n = ng * nxy * nzCMFD
ierr = cudaMemcpy(XSkf, XSkfHost, n, cudaMemcpyHostToDevice)
DEALLOCATE(XskfHost)

PwSum = dotMulti(cuCMFD%phis8, XSkf, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
DEALLOCATE(XSkf)

AvgPw = PwSum * cuTranCMInfo%Inv_FuelVol
IF(lFirst) THEN
  TranInfo%UnitPowerLevel0 = nTracerCntl%PowerLevel / AvgPw
  TranInfo%PwSum0 = PwSum
END IF
Plevel = AvgPw * TranInfo%UnitPowerLevel0
TranInfo%PowerLevel = Plevel
NULLIFY(PlaneMap, PinMap, PinVolFm)

END SUBROUTINE

SUBROUTINE cuSaveTranSol(FmInfo, nTracerCntl)
USE TYPEDEF,          ONLY : FmInfo_Type
USE CNTL,             ONLY : nTracerCntl_Type
USE ioutil,           ONLY : terminate
IMPLICIT NONE
TYPE(FmInfo_Type) :: FmInfo
TYPE(nTracerCntl_Type) :: nTracerCntl

INTEGER :: nxy, nfsr, ng, nzCMFD, nzSub
INTEGER :: myzb, myze, myzbf, myzef
INTEGER :: izf, ixy, ifsr, iz, ig
INTEGER :: n, ierr

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nfsr = cuGeometry%nfsr
nzCMFD = cuDevice%nzCMFD
myzb = cuDevice%myzb
myze = cuDevice%myze
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
!MOC Variables
FmInfo%TranPhi(1:nfsr, myzb:myze, 1:ng) = FmInfo%Phis(1:nfsr, myzb:myze, 1:ng)
FmInfo%TranPsid(1:nfsr, myzb:myze) = FmInfo%TranPsi(1:nfsr, myzb:myze)
FmInfo%TranPsi(1:nfsr, myzb:myze) = FmInfo%Psi(1:nfsr, myzb:myze)
FmInfo%TranPower(1:nfsr, myzb:myze) = FmInfo%Power(1:nfsr, myzb:myze)

!CMFD Variables
n = ng * nxy * nzCMFD
ierr = cublasDcopy_v2(cuDevice%myblasHandle, n, cuCMFD%phis8, 1, cuTranCMInfo%TranPhi, 1)

n = nxy * nzCMFD
ierr = cublasDcopy_v2(cuDevice%myblasHandle, n, cuTranCMInfo%TranPsi, 1, cuTranCMInfo%TranPsid, 1)
ierr = cublasDcopy_v2(cuDevice%myblasHandle, n, cuCMFD%psi8, 1, cuTranCMInfo%TranPsi, 1)

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO iz = myzb, myze
  DO ixy = 1, nxy
    DO ig = 1, ng
      !cuTranCMInfo%TranPhiC(ig, ixy, iz) = cuTranCMInfo%PhiC(ig, ixy, iz)
      !cuTranCMInfo%PhiC(ig, ixy, iz) = cuCMFD%h_phic8(ig, ixy, iz)
      cuTranCMInfo%TranPhiC(ig, ixy, iz) = cuCMFD%h_phic8(ig, ixy, iz)
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

IF(cuCntl%lAxial .AND. nTracerCntl%l3dim) THEN
  nzSub = cuDevice%nzSub
  n = nxy * ng * nzSub
  ierr = cublasCopy(cuDevice%myblasHandle, n, cuAxial%phiCoeff(:,:,:,0), 1, cuAxial%TranPhi, 1)
  n = nxy * nzSub
  ierr = cublasCopy(cuDevice%myblasHandle, n, cuAxial%TranPsi, 1, cuAxial%TranPsid, 1)
  ierr = cublasCopy(cuDevice%myblasHandle, n, cuAxial%PsiCoeff(:,:,0), 1, cuAxial%TranPsi, 1)
END IF


END SUBROUTINE

SUBROUTINE cuUpdtPsiNowStep(Core, FmInfo, GroupInfo, nTracerCntl, eigv0)
USE TYPEDEF,          ONLY : CoreInfo_Type,    FmInfo_Type,      GroupInfo_Type
USE CNTL,             ONLY : nTracerCntl_Type
USE MOC_COMMON,       ONLY : setMocPsi
USE MOC_MOD,          ONLY : psiUpdate
USE CUDA_CMFD,        ONLY : cuCMFDPsiUpdt
USE ioutil,           ONLY : terminate
USE CUDA_AXMOC,       ONLY : cuSetAxPsiNowStep
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
REAL :: eigv0

REAL :: inv_eigv0
INTEGER :: myzb, myze, ng
INTEGER :: n, ierr
LOGICAL :: lxslib

myzb = cuDevice%myzb
myze = cuDevice%myze
ng = cuGeometry%ng

inv_eigv0 = 1./ eigv0
lxslib = nTracerCntl%lxslib

FmInfo%psid(:, myzb : myze) = FmInfo%psi(:, myzb : myze)
IF (nTracerCntl%lMacro) THEN
  CALL SetMOCPsi(Core, FmInfo%phis, FmInfo%psi)
ELSE
  CALL PsiUpdate(Core, FmInfo%Fxr, FmInfo%phis, FmInfo%psi, myzb, myze, ng, lxslib, GroupInfo)
ENDIF
FmInfo%psi = inv_eigv0 * FmInfo%psi

CALL cuSetCMFDSourceOperator(cuCMFD, cuDevice, GroupInfo%lUpscat)
!CALL cuSetTranCMFDSourceOperator(cuCMFD, cuDevice, TranInfo, TranCntl, GroupInfo%lUpscat)
CALL cuCMFDPsiUpdt(cuCMFD, cuDevice)
n = cuGeometry%nxyc * cuDevice%nzCMFD
ierr = cublasDscal_v2(cuDevice%myblasHandle, n, inv_eigv0, cuCMFD%psi8, 1)

IF(cuCntl%lAxial .AND. nTracerCntl%l3dim) THEN
  CALL cuSetAxPsiNowStep(inv_eigv0)
END IF

CALL destroyCsr(cuCMFD%F)
CALL destroyCsr(cuCMFD%Chi)

END SUBROUTINE

SUBROUTINE cuCorrectNowStep(Core, CmInfo, FmInfo, GroupInfo, TranInfo, TranCntl, nTracerCntl)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type, CmInfo_Type, FmInfo_Type, GroupInfo_Type, TranInfo_Type, &
                             TranCntl_Type
USE CNTL,             ONLY : nTracerCntl_Type
USE CUDA_CMFD,        ONLY : cuCopyFlux, cuSetMOCPhis
USE PE_Mod,           ONLY : PE
USE files,            ONLY : io8
USE ioutil,           ONLY : message,           terminate
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl

REAL :: amp, ampK, Know, Kerr
REAL :: corrfactor
INTEGER :: nfsr, myzb, myze
INTEGER :: ng, nxy, nzCMFD
INTEGER :: iz, ifsr, ixy
INTEGER :: n, ierr

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
nfsr = Core%nCoreFsr
myzb = cuDevice%myzb
myze = cuDevice%myze

ampK = TranInfo%lifetime_Dynamic * TranInfo%Factor_F
IF(.NOT. TranCntl%lIQS) TranInfo%AmpPredictor = ampK * TranInfo%Inv_Factor_K0
corrfactor = TranInfo%Amp / TranInfo%AmpPredictor
WRITE(mesg, '(a15, ES11.3, 2(a8,ES14.7))') 'EPKE Correction', corrfactor-1., 'Amp_C:', TranInfo%Amp, 'Amp_P:',TranInfo%AmpPredictor
IF (PE%MASTER .AND. TranCntl%ladptt) CALL message(io8, TRUE, TRUE, mesg)
Know = ampK / TranInfo%amp
Kerr = Know * TranInfo%Inv_Factor_K0 - 1.
WRITE(mesg, '(a13, ES13.3, 2(a8,ES14.7))') 'K Convergence', Kerr, 'K_now:', Know, 'K_0:',1./TranInfo%Inv_Factor_K0
IF (PE%MASTER .AND. TranCntl%ladptt) CALL message(io8, .TRUE., .TRUE., mesg)
IF(TranCntl%lIQS) THEN
  TranInfo%AmpPredictor = TranInfo%Amp
END IF

n = ng * nxy * nzCMFD
ierr = cublasDscal_v2(cuDevice%myblasHandle, n, corrfactor, cuCMFD%phis8, 1)
CALL cuCopyFlux(cuCMFD, cuDevice, 2)

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO iz = myzb, myze
  DO ifsr = 1, nfsr
    FmInfo%phis(ifsr, iz, :) = FmInfo%phis(ifsr, iz, :) * corrfactor
    FmInfo%phim(:, :, ifsr, iz) = FmInfo%phim(:, :, ifsr, iz) * corrfactor
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO iz = myzb, myze
  DO ixy = 1, Core%nxy
    CmInfo%phic(ixy, iz, :) = CmInfo%phic(ixy, iz, :) * corrfactor
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

CALL cuUpdtPsiNowStep(Core, FmInfo, GroupInfo, nTracerCntl, TranInfo%eigv0)

END SUBROUTINE

SUBROUTINE cuCorrectorGuess(GroupInfo, TranInfo, TranCntl, nTracerCntl, l3dim)
USE TYPEDEF,          ONLY : GroupInfo_Type,    TranInfo_Type,    TranCntl_Type
USE CNTL,             ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
LOGICAL :: l3dim

REAL :: amp
REAL :: fmult
INTEGER :: ng, nxy, myzb, myze
INTEGER :: ig, ipin, iz

ng = cuCMFD%ng
nxy = cuGeometry%nxyc
myzb = cuDevice%myzb
myze = cuDevice%myze

CALL cuTranReactivityUpdt_Adj(GroupInfo, TranInfo, TranCntl, l3dim, .FALSE.)
! Save Amplitude
IF(.NOT. TranCntl%lCorrector) TranInfo%amp = TranInfo%lifetime_Dynamic * TranInfo%Factor_F * TranInfo%Inv_Factor_K0
amp = TranInfo%amp
CALL Transient_Corrector(TranInfo, TranCntl, nTracerCntl)
!CALL Transient_Corrector(TranInfo, TranCntl, nTracerCntl, TranCntl%T(TranCntl%nowstep-1), TranCntl%T(TranCntl%nowstep), TranCntl%Tdiv_corrector)

fmult = TranInfo%amp / amp
!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO iz = myzb, myze
  DO ipin = 1, nxy
    DO ig = 1, ng
      cuCMFD%PinXS(ipin, iz)%Phi(ig) = fmult * cuCMFD%PinXS(ipin, iz)%Phi(ig)
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

! recover Amplitude
TranInfo%amp = amp

! recover PKE parameters
TranInfo%delrho_Dynamic = TranInfo%Prev_delrho
TranInfo%corebeta_Dynamic = TranInfo%Prev_corebetat
TranInfo%lifetime_Dynamic = TranInfo%Prev_lifetime
TranInfo%factor_F = TranInfo%Prev_factor_F

END SUBROUTINE

SUBROUTINE cuCorrectorIter(GroupInfo, TranInfo, TranCntl, nTracerCntl, l3dim)
USE PARAM
USE TYPEDEF,          ONLY : GroupInfo_Type,    TranInfo_Type,    TranCntl_Type
USE CNTL,             ONLY : nTracerCntl_Type
USE PE_Mod,           ONLY : PE
USE files,            ONLY : io8
USE ioutil,           ONLY : message,           terminate
USE timer,            ONLY : nTracer_dclock,    TimeChk
IMPLICIT NONE
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
LOGICAL :: l3dim

REAL :: amp, ampK, Know, Kerr
REAL :: fmult
INTEGER :: n, ierr

REAL :: TimeBeg, TimeEnd
REAL, SAVE :: Timee = 0

CALL cuTranReactivityUpdt_Adj(GroupInfo, TranInfo, TranCntl, l3dim, .FALSE.)
! Save Amplitude
amp = TranInfo%amp
CALL Transient_Corrector(TranInfo, TranCntl, nTracerCntl)
!CALL Transient_Corrector(TranInfo, TranCntl, nTracerCntl, TranCntl%T(TranCntl%nowstep-1), TranCntl%T(TranCntl%nowstep), TranCntl%Tdiv_corrector)
ampK = TranInfo%lifetime_Dynamic * TranInfo%Factor_F
IF(.NOT. TranCntl%lIQS) TranInfo%AmpPredictor = ampK * TranInfo%Inv_Factor_K0
fmult = TranInfo%amp / TranInfo%ampPredictor
!print*, traninfo%amp, traninfo%amppredictor
WRITE(mesg, '(a15, ES11.3, 2(a8,ES11.3))') 'EPKE Correction', fmult-1., 'Amp_C:', TranInfo%Amp, 'Amp_P:',TranInfo%AmpPredictor
IF (PE%MASTER) CALL message(io8, .TRUE., .TRUE., mesg)
Know = ampK / TranInfo%amp
Kerr = Know * TranInfo%Inv_Factor_K0 - 1.
WRITE(mesg, '(a13, ES13.3, 2(a8,ES11.3))') 'K Convergence', Kerr, 'K_now:', Know, 'K_0:',1./TranInfo%Inv_Factor_K0
IF (PE%MASTER) CALL message(io8, .TRUE., .TRUE., mesg)
IF(TranCntl%lIQS) THEN
  TranInfo%AmpPredictor = TranInfo%Amp
END IF

n = cuGeometry%nxyc * cuDevice%nzCMFD
ierr = cublasDscal_v2(cuDevice%myblasHandle, n, fmult, cuCMFD%psi8, 1)
n = n * cuGeometry%ng
ierr = cublasDscal_v2(cuDevice%myblasHandle, n, fmult, cuCMFD%phis8, 1)

TimeBeg = nTracer_dclock(FALSE, FALSE)
IF(.NOT. TranCntl%lExpTrsf) THEN
  cuTranCMInfo%Expo = TranInfo%AmpRatio
  cuTranCmInfo%Expo_Alpha = TranInfo%AmpTilt
  IF(cuCntl%lAxial .AND. l3dim) THEN
    cuAxial%Expo = TranInfo%AmpRatio
    cuAxial%Expo_Alpha = TranInfo%AmpTilt
  END IF
  TranInfo%fmExpo = TranInfo%AmpRatio
  TranInfo%fmExpo_Alpha = TranInfo%AmpTilt
END IF
TimeEnd = nTracer_dclock(FALSE, FALSE)
Timee = Timee + TimeEnd - TimeBeg
WRITE(mesg, '(a13, ES13.3, 2(a8,ES11.3))') 'Corrector Update Time  : ', Timee
IF (PE%MASTER) CALL message(io8, .TRUE., .TRUE., mesg)

! recover Amplitude
IF(TranCntl%lIQS) THEN
  TranInfo%ampPredictor = TranInfo%amp
END IF
TranInfo%amp = amp

! recover PKE parameters
TranInfo%delrho_Dynamic = TranInfo%Prev_delrho
TranInfo%corebeta_Dynamic = TranInfo%Prev_corebetat
TranInfo%lifetime_Dynamic = TranInfo%Prev_lifetime
TranInfo%factor_F = TranInfo%Prev_factor_F

END SUBROUTINE


SUBROUTINE cuUpdtPrec(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, nprec)
USE TYPEDEF,          ONLY : CoreInfo_Type,     FmInfo_Type,      TranInfo_Type,      TranCntl_Type,  &
                             FxrInfo_Type,      Pin_Type,         Cell_Type
USE CNTL,             ONLY : nTracerCntl_Type
USE ioutil,           ONLY : terminate
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
INTEGER :: nprec

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
REAL, POINTER :: Psi(:, :), TranPsi(:, :), TranPsid(:, :)
REAL, POINTER :: Omegam(:, :, :), Omega0(:, :, :), Omegap(:, :, :)
REAL :: delt, lambda(nprec), kappa(nprec)
REAL(GPU_NODAL_PRECISION) :: akappa
INTEGER :: nowStep
INTEGER :: FxrIdxSt, FsrIdxSt
INTEGER :: myzb, myze, nxy, nxyc, nzCMFD, nzSub
INTEGER :: ixy, iz, i, j, iprec, ifxr, ifsr, icel
INTEGER :: n, ierr

nowstep = TranCntl%nowstep
delt = TranCntl%delt(nowstep)

Fxr => FmInfo%Fxr
Pin => Core%Pin
Cell => Core%CellInfo

DO iprec = 1, nprec
  lambda(iprec) = TranInfo%lambda(iprec)
  kappa(iprec) = exp(-delt * lambda(iprec))
END DO

myzb = cuDevice%myzb
myze = cuDevice%myze
nzCMFD = cuDevice%nzCMFD

! Fm Precursor
nxy = Core%nxy
Psi => FmInfo%Psi; TranPsi => FmInfo%TranPsi; TranPsid => FmInfo%TranPsid
Omegam => TranInfo%FxrOmegam; Omega0 => TranInfo%FxrOmega0; Omegap => TranInfo%FxrOmegap
!$OMP PARALLEL PRIVATE(FxrIdxSt, FsrIdxSt, icel, ifxr, ifsr)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO iz = myzb, myze
  DO ixy = 1, nxy
    FxrIdxSt = Pin(ixy)%FxrIdxSt
    FsrIdxSt = Pin(ixy)%FsrIdxSt
    icel = Pin(ixy)%Cell(iz)
    DO j = 1, Cell(icel)%nFxr
      ifxr = FxrIdxSt + j - 1
      IF(.NOT. Fxr(ifxr, iz)%lfuel) CYCLE
      DO i = 1, Cell(icel)%nFsrInFxr(j)
        ifsr = FsrIdxSt + Cell(icel)%MapFxr2FsrIdx(i, j) - 1
        DO iprec = 1, nprec
          FmInfo%Prec(iprec, ifsr, iz) = kappa(iprec) * FmInfo%Prec(iprec, ifsr, iz)
          FmInfo%Prec(iprec, ifsr, iz) = FmInfo%Prec(iprec, ifsr, iz) + Omegam(iprec, ifxr, iz) * TranPsid(ifsr, iz)
          FmInfo%Prec(iprec, ifsr, iz) = FmInfo%Prec(iprec, ifsr, iz) + Omega0(iprec, ifxr, iz) * TranPsi(ifsr, iz)
          FmInfo%Prec(iprec, ifsr, iz) = FmInfo%Prec(iprec, ifsr, iz) + Omegap(iprec, ifxr, iz) * Psi(ifsr, iz)
        END DO
      END DO
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
NULLIFY(Psi, TranPsi, TranPsid)
NULLIFY(Omegam, Omega0, Omegap)

! Cm Precursor
nxyc = cuGeometry%nxyc
n = nxyc * nzCMFD
CALL CopyCMFDPrecCoeff(nTracerCntl%l3dim)
DO iprec = 1, nprec
  ierr = cublasDscal_v2(cuDevice%myblasHandle, n, kappa(iprec), cuTranCMInfo%dPrec(:, iprec), 1)
  CALL cuVectorOp('*+', n, cuTranCMInfo%dOmegam(:, iprec), cuTranCMInfo%TranPsid, &
                  cuTranCMInfo%dPrec(:, iprec), cuDevice%myStream)
  CALL cuVectorOp('*+', n, cuTranCMInfo%dOmega0(:, iprec), cuTranCMInfo%TranPsi, &
                  cuTranCMInfo%dPrec(:, iprec), cuDevice%myStream)
  CALL cuVectorOp('*+', n, cuTranCMInfo%dOmegap(:, iprec), cuCMFD%Psi8, &
                  cuTranCMInfo%dPrec(:, iprec), cuDevice%myStream)
END DO

IF(cuCntl%lAxial .AND. nTracerCntl%l3dim) THEN
  nzSub = cuDevice%nzSub
  n = nxyc * nzSub
  DO iprec = 1, nprec
    akappa = kappa(iprec)
    ierr = cublasScal(cuDevice%myblasHandle, n, akappa, cuAxial%Prec(:, :, iprec), 1)
    CALL cuVectorOp('*+', n, cuAxial%Omegam(:, :, iprec), cuAxial%TranPsid, &
                    cuAxial%Prec(:, :, iprec), cuDevice%myStream)
    CALL cuVectorOp('*+', n, cuAxial%Omega0(:, :, iprec), cuAxial%TranPsi, &
                    cuAxial%Prec(:, :, iprec), cuDevice%myStream)
    CALL cuVectorOp('*+', n, cuAxial%Omegap(:, :, iprec), cuAxial%PsiCoeff(:, :, 0), &
                    cuAxial%Prec(:, :, iprec), cuDevice%myStream)
  END DO
  !CALL terminate('Undercondtruction - Transient Axial')
END IF
END SUBROUTINE

SUBROUTINE CopyCMFDPrecCoeff(l3dim)
IMPLICIT NONE
LOGICAL :: l3dim

REAL, ALLOCATABLE :: Omegam(:, :, :), Omega0(:, :, :), Omegap(:, :, :)
REAL(GPU_NODAL_PRECISION), ALLOCATABLE :: Omegama(:, :, :), Omega0a(:, :, :), Omegapa(:, :, :)
INTEGER, POINTER :: fmRange(:, :), PinMap(:)
INTEGER :: nxy, nzCMFD, nprec, nzSub
INTEGER :: myzbf, myzef, myzb, myze, izfb, izfe
INTEGER :: izf, ixy, iprec, iz, ixy_map
INTEGER :: n, ierr

nxy = cuGeometry%nxyc
nprec = cuGeometry%nprec
nzCMFD = cuDevice%nzCMFD
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef
myzb = cuDevice%myzb
myze = cuDevice%myze
fmRange => cuGeometry%fmRange
PinMap => cuGeometry%PinMap

ALLOCATE(Omegam(nxy, myzbf:myzef, nprec))
ALLOCATE(Omega0(nxy, myzbf:myzef, nprec))
ALLOCATE(Omegap(nxy, myzbf:myzef, nprec))

!$OMP PARALLEL PRIVATE(ixy_map)
DO iz = myzb, myze
  izfb = fmRange(iz, 1)
  izfe = fmRange(iz, 2)
  !$OMP DO SCHEDULE(GUIDED)
  DO ixy = 1, nxy
    ixy_map = PinMap(ixy)
    DO iprec = 1, nprec
      DO izf = izfb, izfe
        Omegam(ixy, izf, iprec) = cuTranCMInfo%CellOmegam(iprec, ixy_map, iz)
        Omega0(ixy, izf, iprec) = cuTranCMInfo%CellOmega0(iprec, ixy_map, iz)
        Omegap(ixy, izf, iprec) = cuTranCMInfo%CellOmegap(iprec, ixy_map, iz)
      END DO
    END DO
  END DO
  !$OMP END DO
END DO
!$OMP END PARALLEL

n = nprec * nxy * nzCMFD
ierr = cudaMemcpy(cuTranCMInfo%dOmegam, Omegam, n, cudaMemcpyHostToDevice)
ierr = cudaMemcpy(cuTranCMInfo%dOmega0, Omega0, n, cudaMemcpyHostToDevice)
ierr = cudaMemcpy(cuTranCMInfo%dOmegap, Omegap, n, cudaMemcpyHostToDevice)

DEALLOCATE(Omegam, Omega0, Omegap)

IF(cuCntl%lAxial .AND. l3dim) THEN
  nzSub = cuDevice%nzSub
  ALLOCATE(Omegama(nxy, nzSub, nprec))
  ALLOCATE(Omega0a(nxy, nzSub, nprec))
  ALLOCATE(Omegapa(nxy, nzSub, nprec))
  !!$OMP PARALLEL PRIVATE(ixy_map)
  DO iz = myzb, myze
    izfb = cuDevice%cmSubRange(iz, 1)
    izfe = cuDevice%cmSubRange(iz, 2)
    !!$OMP DO SCHEDULE(GUIDED)
    DO ixy = 1, nxy
      ixy_map = PinMap(ixy)
      DO iprec = 1, nprec
        DO izf = izfb, izfe
          Omegama(ixy, izf, iprec) = cuTranCMInfo%CellOmegam(iprec, ixy_map, iz)
          Omega0a(ixy, izf, iprec) = cuTranCMInfo%CellOmega0(iprec, ixy_map, iz)
          Omegapa(ixy, izf, iprec) = cuTranCMInfo%CellOmegap(iprec, ixy_map, iz)
        END DO
      END DO
    END DO
    !!$OMP END DO
  END DO
  !!$OMP END PARALLEL
  n = nprec * nxy * nzSub
  ierr = cudaMemcpy(cuAxial%Omegam, Omegama, n, cudaMemcpyHostToDevice)
  ierr = cudaMemcpy(cuAxial%Omega0, Omega0a, n, cudaMemcpyHostToDevice)
  ierr = cudaMemcpy(cuAxial%Omegap, Omegapa, n, cudaMemcpyHostToDevice)
  DEALLOCATE(Omegama, Omega0a, Omegapa)
END IF

NULLIFY(fmRange, PinMap)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuComputeTranCMFDSrc_chidk(trSrc, resSrc, precSrcK, TranPhi, VolInvVel, expo_alpha, expo, delt, nzCMFD)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE
REAL(8), DEVICE :: trSrc(:,:), resSrc(:,:), precSrcK(:,:), TranPhi(:,:)
REAL(8), DEVICE :: VolInvVel(:,:), expo_alpha(:,:), expo(:,:)
REAL(8), VALUE :: delt
INTEGER, VALUE :: nzCMFD

REAL(8) :: thetah, prevSrc
INTEGER :: ncel
INTEGER :: ig, icel, iprec

ncel = nxyc * nzCMFD

ig = threadIdx%x
icel = threadIdx%y + (blockIdx%x - 1) * blockDim%y
IF(icel .GT. ncel) RETURN

thetah = 1./ theta - 1.

trSrc(ig, icel) = 0.
DO iprec = 1, nprec
  trSrc(ig, icel) =  trSrc(ig,icel) + chidk(ig, iprec) * PrecSrcK(icel, iprec)
END DO
prevSrc = Thetah * (ResSrc(ig, icel) - expo_alpha(ig, icel) * VolInvVel(ig, icel) * TranPhi(ig, icel))
prevSrc = prevSrc + VolInvVel(ig, icel) * TranPhi(ig, icel) / (delt * theta)
!prevSrc =  VolInvVel(ig, icel) / (delt * theta)
prevSrc = prevSrc * Expo(ig, icel)
trSrc(ig, icel) = trSrc(ig, icel) + prevSrc

END SUBROUTINE



ATTRIBUTES(GLOBAL) SUBROUTINE cuComputeBetaAverage(BetaWeight, Beta, Psi, nzCMFD)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE
REAL(8), DEVICE :: BetaWeight(:,:), Beta(:), Psi(:)
INTEGER, VALUE :: nzCMFD

INTEGER :: ncel
INTEGER :: ig, icel

ncel = nxyc * nzCMFD

ig = threadIdx%x
icel = threadIdx%y + (blockIdx%x - 1) * blockDim%y
IF(icel .GT. ncel) RETURN

BetaWeight(ig, icel) = chid(ig) * Beta(icel) * Psi(icel)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuComputeBetaAverage_iprec(BetaWeight, Beta, Psi, nzCMFD, iprec)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE
REAL(8), DEVICE :: BetaWeight(:,:), Beta(:), Psi(:)
INTEGER, VALUE :: nzCMFD
INTEGER, VALUE :: iprec

INTEGER :: ncel
INTEGER :: ig, icel

ncel = nxyc * nzCMFD

ig = threadIdx%x
icel = threadIdx%y + (blockIdx%x - 1) * blockDim%y
IF(icel .GT. ncel) RETURN

BetaWeight(ig, icel) = chidk(ig, iprec) * Beta(icel) * Psi(icel)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuComputeDelayedSrc(trSrc, precSrc, nzCMFD)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE
REAL(8), DEVICE :: trSrc(:,:), precSrc(*)
INTEGER, VALUE :: nzCMFD

INTEGER :: ncel
INTEGER :: ig, icel

ncel = nxyc * nzCMFD

ig = threadIdx%x
icel = threadIdx%y + (blockIdx%x - 1) * blockDim%y
IF(icel .GT. ncel) RETURN

trSrc(ig, icel) = chid(ig) * PrecSrc(icel)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuComputeDelayedSrc_iprec(trSrc, precSrc, nzCMFD, iprec)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE
REAL(8), DEVICE :: trSrc(:,:), precSrc(*)
INTEGER, VALUE :: nzCMFD
INTEGER, VALUE ::  iprec

INTEGER :: ncel
INTEGER :: ig, icel

ncel = nxyc * nzCMFD

ig = threadIdx%x
icel = threadIdx%y + (blockIdx%x - 1) * blockDim%y
IF(icel .GT. ncel) RETURN

trSrc(ig, icel) = chidk(ig, iprec) * PrecSrc(icel)

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuComputeDelayedSrc_chidk(trSrc, precSrcK, nzCMFD)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE
REAL(8), DEVICE :: trSrc(:,:), precSrcK(:, :)
INTEGER, VALUE :: nzCMFD

INTEGER :: ncel
INTEGER :: ig, icel, iprec

ncel = nxyc * nzCMFD

ig = threadIdx%x
icel = threadIdx%y + (blockIdx%x - 1) * blockDim%y
IF(icel .GT. ncel) RETURN

trSrc(ig, icel) = 0.
DO iprec = 1, 6
  trSrc(ig, icel) = chidk(ig, iprec) * PrecSrcK(icel, iprec)
END DO

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuCopyPrecKernel(cuDevice, CMPrec, SubPrec, SubPsi)
USE CUDA_CONST
USE CUDAFOR
IMPLICIT NONE
TYPE(cuDevice_Type), INTENT(IN) :: cuDevice
REAL(8), DEVICE, INTENT(IN) :: CMPrec(:, :)
REAL(GPU_NODAL_PRECISION), DEVICE, INTENT(IN) :: SubPsi(:, :, 0:)
REAL(GPU_NODAL_PRECISION), DEVICE :: subPrec(:, :, :)

REAL(GPU_NODAL_PRECISION) :: psisum, fmult
INTEGER :: ipin, iz, izf, iprec, icel

ipin = (blockIdx%x - 1) * blockDim%x + threadIdx%x

IF(ipin .GT. nxyc) RETURN

DO iz = cuDevice%myzbf, cuDevice%myzef
  icel = (iz-cuDevice%myzbf) * nxyc + ipin
  psisum = 0.
  DO izf = cuDevice%fmSubRange(iz, 1), cuDevice%fmSubrange(iz, 2)
    psisum = psisum + SubPsi(ipin, izf, 0)
  END DO
  psisum = psisum / REAL((cuDevice%fmSubRange(iz, 2) - cuDevice%fmSubRange(iz,1) + 1))

  IF(psisum .GT. 0) THEN
    DO izf = cuDevice%fmSubRange(iz, 1), cuDevice%fmSubrange(iz, 2)
      fmult = SubPsi(ipin, izf, 0) / psisum
      DO iprec = 1, nprec
        SubPrec(ipin, izf, iprec) = fmult * CMPrec(icel, iprec)
      END DO
    END DO
  ELSE
    DO izf = cuDevice%fmSubRange(iz, 1), cuDevice%fmSubrange(iz, 2)
      DO iprec = 1, nprec
        SubPrec(ipin, izf, iprec) = 0.0
      END DO
    END DO
  END IF
END DO

END SUBROUTINE

FUNCTION cuTranSrcNorm(eigv0) RESULT(snorm)
USE CUDA_CMFD,        ONLY : cuCMFDPsiUpdt
IMPLICIT NONE
REAL :: snorm
REAL :: eigv0

REAL(8), ALLOCATABLE, DEVICE :: b(:)
REAL(8), POINTER, DEVICE :: csrVal(:)
INTEGER, POINTER, DEVICE :: csrRowPtr(:), csrColIdx(:)
REAL :: inv_eigv0
INTEGER :: nr, nc, nnz
INTEGER :: n, ierr

CALL cuCMFDPsiUpdt(cuCMFD, cuDevice)
n = cuGeometry%nxyc * cuDevice%nzCMFD
inv_eigv0 = 1./eigv0
ierr = cublasDscal_v2(cuDevice%myblasHandle, n, inv_eigv0, cuCMFD%psi8, 1)

n = n * cuGeometry%ng
ALLOCATE(b(n))
ierr = cublasDcopy_v2(cuDevice%myblasHandle, n, cuCMFD%src8, 1, b, 1)

csrVal => cuCMFD%S%d_csrVal
csrRowPtr => cuCMFD%S%d_csrRowPtr
csrColIdx => cuCMFD%S%d_csrColIdx
nr = cuCMFD%S%nr
nc = cuCMFD%S%nc
nnz = cuCMFD%S%nnz

ierr = cusparseDcsrmv(cuDevice%mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, 1.0_8,                &
                      cuCMFD%S%descr, csrVal, csrRowPtr, csrColIdx, cuCMFD%phis8, 1.0_8, b)

IF(cuDevice%lFuel) THEN
  csrVal => cuCMFD%Chi%d_csrVal
  csrRowPtr => cuCMFD%Chi%d_csrRowPtr
  csrColIdx => cuCMFD%Chi%d_csrColIdx
  nr = cuCMFD%Chi%nr
  nc = cuCMFD%Chi%nc
  nnz = cuCMFD%Chi%nnz
  ierr = cusparseDcsrmv(cuDevice%mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, 1.0_8,                &
    cuCMFD%Chi%descr, csrVal, csrRowPtr, csrColIdx, cuCMFD%psi8, 1.0_8, b)
END IF

snorm = normMulti(b, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
NULLIFY(csrVal, csrRowPtr, csrColIdx)
DEALLOCATE(b)

END FUNCTION

SUBROUTINE cuCalcFisRate(fisRate, nxy, myzb, myze, ng)
USE TYPEDEF,          ONLY : PinXs_Type
IMPLICIT NONE
REAL :: FisRate(nxy, myzb:myze)
INTEGER :: nxy, myzb, myze, ng

REAL, POINTER :: hz(:), hzfm(:)
INTEGER, POINTER :: pinMap(:), fmRange(:,:)
REAL :: hmult, frsum
INTEGER :: nxyc
INTEGER :: ixy, ixy_map, ig, ipin, iz, izf
INTEGER :: i

PinMap => cuGeometry%PinMap
fmRange => cuGeometry%fmRange
hz => cuGeometry%hz
hzfm => cuGeometry%hzfm
nxyc = cuGeometry%nxyc

!$OMP PARALLEL PRIVATE(ixy_map, frsum , hmult, ipin)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO iz = myzb, myze
  DO ixy = 1, nxyc
    ixy_map = pinMap(ixy)
    frsum = 0.
    DO izf = fmRange(iz, 1), fmRange(iz, 2)
      hmult = hzfm(izf) / hz(iz)
      DO ig = 1, ng
        frsum = frsum + cuCMFD%PinXS(ixy_map, iz)%xskf(ig) * cuCMFD%h_phis8(ig, ixy, izf) * hmult
      END DO
      DO i = 1, cuGeometry%superPin(ixy_map)%nxy
        ipin = cuGeometry%superPin(ixy_map)%pin(i)
        fisRate(ipin, iz) = frsum
      END DO
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE cuTranReac_dbg(GroupInfo, TranInfo, l3dim, lsave)
USE PARAM
USE TYPEDEF,          ONLY : GroupInfo_Type, TranInfo_Type
USE CUDA_CMFD
USE CUDA_INIT,        ONLY : DeallocCMFDVar
USE PE_Mod,           ONLY : PE
USE files,            ONLY : io8
USE ioutil,           ONLY : message,           terminate
IMPLICIT NONE
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranInfo_Type) ::  TranInfo
LOGICAL :: l3dim, lsave

TYPE(dim3) :: Blocks, Threads
REAL(8), ALLOCATABLE, DEVICE :: InvVel(:), FisSrc(:), Migration(:), Residual(:), BetaWeight(:, :)
REAL(8), ALLOCATABLE, DEVICE :: Beta(:)
REAL(8), ALLOCATABLE, DEVICE :: BetakWeight(:,:,:), Betak(:,:)
REAL(8), POINTER, DEVICE :: csrVal(:)
INTEGER, POINTER, DEVICE :: csrRowPtr(:), csrColIdx(:)
REAL, POINTER :: BetaHost(:, :), BetakHost(:,:,:)
REAL :: factor_F, avg_res, avg_beta, avg_invVel
INTEGER, POINTER :: PlaneMap(:), PinMap(:)
INTEGER :: bottomRange(2), topRange(2)
INTEGER :: n, ng, nxy, nzCMFD, myzbf, myzef, ncel, nprec
INTEGER :: ixy, ixy_map, iz, izf, iprec
INTEGER :: nr, nc, nnz
INTEGER :: ierr

REAL :: lifetime, delrho, traneig, corebeta, reactivity

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
ncel = nxy * nzCMFD
nprec = TranInfo%nprec
n = ng * nxy * nzCMFD
bottomRange = cuCMFD%bottomRange
topRange = cuCMFD%topRange
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef

PlaneMap => cuCMFD%PlaneMap
PinMap => cuGeometry%PinMap

ALLOCATE(InvVel(n), FisSrc(n), Migration(n), Residual(n))
ALLOCATE(Beta(nxy * nzCMFD), BetaHost(nxy, myzbf:myzef), BetaWeight(ng, nxy*nzCMFD))
ALLOCATE(Betak(nxy * nzCMFD, nprec), BetakHost(nxy, myzbf:myzef, nprec), BetakWeight(ng, nxy*nzCMFD, nprec))

IF(lsave) THEN
  ALLOCATE(psidbg(nxy*nzCMFD), phidbg(ng, nxy*nzCMFD))
  ierr = cublasDcopy_v2(cuDevice%myblasHandle, n, cuCMFD%phis8, 1, phidbg, 1)
END IF

CALL cuSetNaturalBiCGSystem(cuCMFD, cuDevice, l3dim, .FALSE., 0.0)
CALL cuSetCMFDSourceOperator(cuCMFD, cuDevice, GroupInfo%lUpscat)
CALL cuSetInvVelo()

IF (cuDevice%lFuel) THEN
  ierr = cusparseDcsrmv(cuDevice%mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, cuCMFD%F%nr, cuCMFD%F%nc, cuCMFD%F%nnz, 1.0_8,                &
    cuCMFD%F%descr, cuCMFD%F%d_csrVal, cuCMFD%F%d_csrRowPtr, cuCMFD%F%d_csrColIdx, phidbg, 0.0_8, psidbg)
  ierr = cublasDscal_v2(cuDevice%myblasHandle, nxy*nzCMFD, 1./TranInfo%eigv0, psidbg, 1)
ELSE
  psidbg = 0.
END IF

!$OMP PARALLEL PRIVATE(iz, ixy_map)
DO izf = myzbf, myzef
  iz = PlaneMap(izf)
  !$OMP DO SCHEDULE(GUIDED)
  DO ixy = 1, nxy
    ixy_map = PinMap(ixy)
    BetaHost(ixy, izf) = cuCMFD%PinXS(ixy_map, iz)%betat
  END DO
  !$OMP END DO
END DO
!$OMP END PARALLEL
ierr = cudaMemcpy(Beta, BetaHost, nxy*nzCMFD, cudaMemcpyHostToDevice)
DEALLOCATE(BetaHost)

Threads = dim3(ng, cuDevice%cuMaxThreadPerBlock/ng, 1)
Blocks = dim3(ncel / Threads%y + 1, 1, 1)

CALL cuComputeBetaAverage <<< Blocks, Threads, 0, cuDevice%myStream>>>                             &
     (BetaWeight, Beta, psidbg, nzCMFD)

CALL cuVectorOp('*', n, phidbg, cuTranCMInfo%VolInvVel, InvVel, cuDevice%myStream)
IF(cuDevice%lFuel) THEN
  csrVal => cuCMFD%Chi%d_csrVal
  csrRowPtr => cuCMFD%Chi%d_csrRowPtr
  csrColIdx => cuCMFD%Chi%d_csrColIdx
  nr = cuCMFD%Chi%nr
  nc = cuCMFD%Chi%nc
  nnz = cuCMFD%Chi%nnz
  ierr = cusparseDcsrmv(cuDevice%mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, 1.0_8,                &
    cuCMFD%Chi%descr, csrVal, csrRowPtr, csrColIdx, psidbg, 0.0_8, FisSrc)
ELSE
  FisSrc = 0.
END IF
CALL cuMatMul3D(cuCMFD%M, cuCMFD%offDiag8, phidbg, Migration, bottomRange, topRange,                               &
                  cuDevice%mySparseHandle, cuDevice%myStream, MPI_CUDA_COMM)
CALL cuVectorOp('-', n, FisSrc, Migration, Residual, cuDevice%myStream)

factor_F = dotMulti(cuCMFD_adj%phis8, FisSrc, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
avg_invVel = dotMulti(cuCMFD_adj%phis8, invVel, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
avg_res = dotMulti(cuCMFD_adj%phis8, Residual, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
avg_beta = dotMulti(cuCMFD_adj%phis8, BetaWeight, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)

avg_res = avg_res / factor_F
lifetime = avg_invVel / factor_F
delRho = avg_res * 1.e+5
TranEig = 1./(1. - avg_res)
CoreBeta = avg_beta / factor_F
Reactivity = avg_res / CoreBeta
Factor_F = factor_F

WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, 1e14.7, x, a, 1e14.7, x, a, 1e14.7)') '@13', 1.,"s, P=", TranInfo%PowerLevel * 100, ', $= ', reactivity, ', G=', lifetime
CALL message(io8, TRUE, TRUE, MESG)
WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, F7.1, x, a, F8.5, x, a, 1e14.7)') '@14',1. ,"s, r(pcm)= ", &
 delrho/10._8, ', k= ', TranEig/10._8, ', B= ', CoreBeta
CALL message(io8, TRUE, TRUE, MESG)
WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, 1e14.7, x, a, 1e14.7, x, a, 1e14.7)') '@15', 1.,"s, P=", TranInfo%PowerLevel * 100, ', $= ', avg_res, ', F=', factor_F
CALL message(io8, TRUE, TRUE, MESG)


CALL DeallocCMFDVar(cuCMFD)
DEALLOCATE(Beta, BetaWeight)
DEALLOCATE(Betak, BetakWeight)
DEALLOCATE(InvVel, FisSrc, Residual, Migration)

  END SUBROUTINE

SUBROUTINE cuTransient_Monitor(GroupInfo, TranInfo, l3dim)
USE PARAM
USE TYPEDEF,          ONLY : GroupInfo_Type, TranInfo_Type
!USE CUDA_CMFD
USE CUDA_INIT,        ONLY : DeallocCMFDVar
USE PE_Mod,           ONLY : PE
USE files,            ONLY : io8
USE ioutil,           ONLY : message,           terminate
IMPLICIT NONE
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranInfo_Type) ::  TranInfo
LOGICAL :: l3dim

TYPE(dim3) :: Blocks, Threads
TYPE(CSR_DOUBLE) :: L, A, S, F
TYPE(superPin_Type), POINTER :: Pin(:)
INTEGER, POINTER :: planeMap(:), pinMap(:), pinMapRev(:)
REAL(8), ALLOCATABLE, DEVICE :: ScatSrc(:), FisSrc(:), Leakage(:), Absorption(:), BetaWeight(:, :)
REAL(8), ALLOCATABLE, DEVICE :: Beta(:)
REAL(8), ALLOCATABLE, DEVICE :: BetakWeight(:,:,:), Betak(:,:)
REAL(8), POINTER, DEVICE :: csrVal(:)
INTEGER, POINTER, DEVICE :: csrRowPtr(:), csrColIdx(:)
REAL(8), ALLOCATABLE :: offDiag(:, :, :)
REAL(8), ALLOCATABLE :: offDiag8(:, :, :)
REAL, POINTER :: BetaHost(:, :), BetakHost(:,:,:)
REAL, POINTER :: PinVolFm(:, :), hzfm(:)
REAL :: factor_F, avg_res, avg_beta, avg_invVel
REAL :: avg_fis, avg_leak, avg_scat, avg_abs
INTEGER :: bottomRange(2), topRange(2)
INTEGER :: n, ng, nxy, nzCMFD, myzbf, myzef, ncel, nprec
INTEGER :: ixy, ixy_map, iz, izf, iprec, ir, ibd, isurf, ineighpin, dz, ic
INTEGER :: ig, jg, gb, ge
INTEGER :: nr, nc, nnz
INTEGER :: ierr
INTEGER, PARAMETER :: DOWN = 5, UP = 6, SELF = 7
INTEGER :: surf_(7) = (/ UP, NORTH, WEST, SELF, EAST, SOUTH, DOWN /)
REAL :: diagVal(7), Dtil, Dhat, val

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
ncel = nxy * nzCMFD
nprec = TranInfo%nprec
n = ng * nxy * nzCMFD
bottomRange = cuCMFD%bottomRange
topRange = cuCMFD%topRange
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef

Pin => cuGeometry%superPin
PinVolFm => cuGeometry%PinVolFm
hzfm => cuGeometry%hzfm
PlaneMap => cuCMFD%PlaneMap
PinMap => cuGeometry%PinMap
pinMapRev => cuGeometry%pinMapRev

ALLOCATE(ScatSrc(n), FisSrc(n), Leakage(n), Absorption(n))
ALLOCATE(Beta(nxy * nzCMFD), BetaHost(nxy, myzbf:myzef), BetaWeight(ng, nxy*nzCMFD))
ALLOCATE(Betak(nxy * nzCMFD, nprec), BetakHost(nxy, myzbf:myzef, nprec), BetakWeight(ng, nxy*nzCMFD, nprec))

CALL createCSR(L, 7 * ng * nxy * nzCMFD, ng * nxy * nzCMFD, ng * nxy * nzCMFD)
CALL createCSR(A, ng * nxy * nzCMFD, ng * nxy * nzCMFD, ng * nxy * nzCMFD)
CALL createCSR(S, ng * ng * nxy * nzCMFD, ng * nxy * nzCMFD, ng * nxy * nzCMFD)
CALL createCSR(F, ng * ng * nxy * nzCMFD, ng * nxy * nzCMFD, ng * nxy * nzCMFD)
CALL cuSetInvVelo()

DO izf = myzbf, myzef
  iz = planeMap(izf)
  DO ixy = 1, nxy
    ixy_map = pinMap(ixy)
    DO ig = 1, ng
      ir = ig + (ixy-1) * ng + (izf - myzbf) * ng * nxy
      diagVal = 0.0
      DO ibd = 1, 4
        Dtil = cuCMFD%PinXS(ixy_map, iz)%Dtil(ibd, ig)
        Dhat = cuCMFD%PinXS(ixy_map, iz)%Dhat(ibd, ig)
        diagVal(ibd) = - (Dtil + Dhat) * hzfm(izf)
        diagVal(SELF) = diagVal(SELF) + (Dtil - Dhat) * hzfm(izf)
      END DO
      IF(l3dim) THEN
        DO ibd = 5,6
          Dtil = cuCMFD%AxDtil(ibd-4, ig, ixy, izf)
          Dhat = cuCMFD%AxDhat(ibd-4, ig, ixy, izf)
          diagVal(ibd) = -(Dtil + Dhat) * PinVolFm(ixy_map, izf) / hzfm(izf)
          diagVal(SELF) = diagVal(SELF) + (Dtil - Dhat) * PinVolFm(ixy_map, izf) / hzfm(izf)
        END DO
      END IF
      DO ibd = 1, 7
        isurf = surf_(ibd)
        SELECT CASE(isurf)
        CASE(SELF)
          ineighpin = ixy
          dz = 0
        CASE (UP)
          ineighpin = ixy
          dz = 1
        CASE (DOWN)
          ineighpin = ixy
          dz = -1
        CASE (NORTH, WEST, EAST, SOUTH)
          ineighpin = pin(ixy_map)%neighidx(isurf)
          ineighpin = pinmaprev(ineighpin)
          IF(ineighpin .LE. 0) diagVal(isurf) = 0.
          IF(ineighpin .EQ. ixy) THEN
            Dtil = cuCMFD%PinXS(ixy_map, iz)%Dtil(isurf, ig)
            Dhat = cuCMFD%PinXS(ixy_map, iz)%Dhat(isurf, ig)
            diagVal(SELF) = diagVal(SELF) - (Dtil +Dhat) * hzfm(izf)
            diagVal(isurf) = 0.
          END IF
          dz = 0
        END SELECT
        ic = ir + (ineighpin - ixy) * ng + dz * ng * nxy
        val = diagVal(isurf)
        CALL pushCsr(L, val, ir, ic)
      END DO
      val = PinVolFm(ixy_map,  izf) * cuCMFD%PinXS(ixy_map, iz)%xsr(ig)
      CALL pushCsr(A, val, ir, ir)

      gb = cuCMFD%PinXS(ixy_map, iz)%XSs(ig)%ib
      ge = cuCMFD%PinXS(ixy_map, iz)%XSs(ig)%ie
      DO jg = 1, ng
        IF(jg .EQ. ig) CYCLE
        ic = ir + (jg-ig)
        val = - cuCMFD%PinXS(ixy_map, iz)%Xsnf(jg) * cuCMFD%PinXS(ixy_map, iz)%Chi(ig) * PinVolFm(ixy_map, izf)
        CALL pushCsr(F, val, ir, ic)
        IF(jg.GE. gb .AND. jg .LE. ge) THEN
          val = - cuCMFD%PinXS(ixy_map, iz)%xss(ig)%from(jg) * PinVolFm(ixy_map, izf)
          CALL pushCsr(S, val, ir, ic)
        END IF
      END DO
    END DO
  END DO
END DO

CALL finalizeSortCsr(L, .TRUE.)
CALL finalizeSortCsr(A, .TRUE.)
CALL finalizeSortCsr(S, .TRUE.)
CALL finalizeSortCsr(F, .TRUE.)

IF (l3dim) THEN

  ALLOCATE(offDiag(ng, nxy, 2), offDiag8(ng, nxy, 2))
  offDiag = 0.0; offDiag8 = 0.0

  IF (.NOT. cuDevice%lBottom) THEN
    !$OMP PARALLEL PRIVATE(Dtil, Dhat, ixy_map, val)
    !$OMP DO SCHEDULE(GUIDED)
    DO ixy = 1, nxy
      ixy_map = pinMap(ixy)
      DO ig = 1, ng
        Dtil = cuCMFD%AxDtil(bottom, ig, ixy, myzbf)
        Dhat = cuCMFD%AxDhat(bottom, ig, ixy, myzbf)
        val = - (Dtil + Dhat) * PinVolFm(ixy_map, myzbf) / hzfm(myzbf)
        offDiag(ig, ixy, bottom) = val
        offDiag8(ig, ixy, bottom) = val
      ENDDO
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL
  ENDIF
  IF (.NOT. cuDevice%lTop) THEN
    !$OMP PARALLEL PRIVATE(Dtil, Dhat, ixy_map, val)
    !$OMP DO SCHEDULE(GUIDED)
    DO ixy = 1, nxy
      ixy_map = pinMap(ixy)
      DO ig = 1, ng
        Dtil = cuCMFD%AxDtil(top, ig, ixy, myzef)
        Dhat = cuCMFD%AxDhat(top, ig, ixy, myzef)
        val = - (Dtil + Dhat) * PinVolFm(ixy_map, myzef) / hzfm(myzef)
        offDiag(ig, ixy, top) = val
        offDiag8(ig, ixy, top) = val
      ENDDO
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL
  ENDIF

  ierr = cudaMemcpy(cuCMFD%offDiag, offDiag, ng * nxy * 2, cudaMemcpyHostToDevice)
  ierr = cudaMemcpy(cuCMFD%offDiag8, offDiag8, ng * nxy * 2, cudaMemcpyHostToDevice)

  DEALLOCATE(offDiag, offDiag8)

ENDIF

Threads = dim3(ng, cuDevice%cuMaxThreadPerBlock/ng, 1)
Blocks = dim3(ncel / Threads%y + 1, 1, 1)

CALL cuMatMul3D(L, cuCMFD%offDiag8, cuCMFD%phis8, Leakage, bottomRange, topRange,                               &
                  cuDevice%mySparseHandle, cuDevice%myStream, MPI_CUDA_COMM)
csrVal => S%d_csrVal
csrRowPtr => S%d_csrRowPtr
csrColIdx => S%d_csrColIdx
nr = S%nr
nc = S%nc
nnz = S%nnz
ierr = cusparseDcsrmv(cuDevice%mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, 1.0_8,                &
  S%descr, csrVal, csrRowPtr, csrColIdx, cuCMFD%phis8, 0.0_8, ScatSrc)
csrVal => F%d_csrVal
csrRowPtr => F%d_csrRowPtr
csrColIdx => F%d_csrColIdx
nr = F%nr
nc = F%nc
nnz = F%nnz
ierr = cusparseDcsrmv(cuDevice%mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, 1.0_8,                &
  F%descr, csrVal, csrRowPtr, csrColIdx, cuCMFD%phis8, 0.0_8, FisSrc)
csrVal => A%d_csrVal
csrRowPtr => A%d_csrRowPtr
csrColIdx => A%d_csrColIdx
nr = A%nr
nc = A%nc
nnz = A%nnz
ierr = cusparseDcsrmv(cuDevice%mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, 1.0_8,                &
  A%descr, csrVal, csrRowPtr, csrColIdx, cuCMFD%phis8, 0.0_8, Absorption)

avg_leak = dotMulti(cuCMFD_adj%phis8, Leakage, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
avg_scat = dotMulti(cuCMFD_adj%phis8, ScatSrc, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
avg_abs = dotMulti(cuCMFD_adj%phis8, Absorption, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
avg_fis = dotMulti(cuCMFD_adj%phis8, FisSrc, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)

CALL DeallocCMFDVar(cuCMFD)
DEALLOCATE(Beta, BetaWeight)
DEALLOCATE(Betak, BetakWeight)
DEALLOCATE(ScatSrc, FisSrc, Leakage, Absorption)

WRITE(mesg, '(a, x, a, 1p, 1e14.6, x, a, 1p, 1e14.6, x, a, 1e14.6, x, a, 1e14.6)') '@13', "L=", avg_leak, ', A= ', avg_abs, ', F=', avg_fis, ', S=', avg_scat
IF(PE%master) CALL message(io8, TRUE, TRUE, MESG)
END SUBROUTINE

SUBROUTINE cuCondiMOC_Monitor(Core)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type
USE PE_Mod,           ONLY : PE
USE files,            ONLY : io8
USE ioutil,           ONLY : message,           terminate
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core

TYPE(CSR_DOUBLE), POINTER :: L(:), tranL(:)
TYPE(superPin_Type), POINTER :: Pin(:)
REAL, POINTER, DEVICE :: Leakage(:,:), tranLeakage(:,:), Ldiff(:,:)
REAL, POINTER :: hLeakage(:,:), hLdiff(:,:)
REAL, POINTER :: PinVolFm(:,:), hzfm(:)
REAL(8), POINTER, DEVICE :: csrVal(:)
INTEGER, POINTER, DEVICE :: csrRowPtr(:), csrColIdx(:)
INTEGER, POINTER :: PlaneMap(:), PinMap(:), PinMapRev(:)
REAL :: Dtil, Dhat
REAL :: val, maxval, sumval, sumleak, totleak
REAL :: realn, totn
REAL :: fmult
INTEGER :: bottomRange(2), topRange(2)
INTEGER :: ng, nxy, nzCMFD, myzb, myze, myzbf, myzef, nsurf
INTEGER :: n, nr, nc, nnz, FsrIdxSt, nLocalFsr
INTEGER :: i, j, ig, ixy, ixy_map, iz, izf, ibd, ipin, icel, ifsr
INTEGER :: ir, ic, ineigh
INTEGER :: ierr
CHARACTER(20) :: DirName(4) = (/'SOUTH', 'WEST', 'NORTH', 'EAST'/)

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
nsurf = 4

bottomRange = cuCMFD%bottomRange
topRange = cuCMFD%topRange
myzb = cuDevice%myzb
myze = cuDevice%myze
myzbf = cuDevice%myzbf
myzef = cuDevice%myzef

Pin => cuGeometry%superPin
PinVolFm => cuGeometry%PinVolFm
hzfm => cuGeometry%hzfm
PlaneMap => cuCMFD%PlaneMap
PinMap => cuGeometry%PinMap
PinMapRev => cuGeometry%pinMapRev

IF(lalloc) THEN
  ALLOCATE(TranDhat(4, ng, nxy, myzb:myze), TranDtil(4, ng, nxy, myzb:myze))
  lalloc = .FALSE.
  GO TO 111
END IF

ALLOCATE(Leakage(ng * nxy * nzCMFD, nsurf))
ALLOCATE(tranLeakage(ng * nxy * nzCMFD, nsurf))
ALLOCATE(Ldiff(ng * nxy * nzCMFD, nsurf))
ALLOCATE(hLeakage(ng * nxy * nzCMFD, nsurf))
ALLOCATE(hLdiff(ng * nxy * nzCMFD, nsurf))
ALLOCATE(L(nsurf), tranL(nsurf))
DO ibd = 1, nsurf
  CALL createCSR(L(ibd), 2*ng * nxy * nzCMFD, ng * nxy * nzCMFD, ng * nxy * nzCMFD)
  CALL createCSR(tranL(ibd), 2*ng * nxy * nzCMFD, ng * nxy * nzCMFD, ng * nxy * nzCMFD)
END DO

DO izf = myzbf, myzef
  iz = planeMap(izf)
  DO ixy = 1, nxy
    ixy_map = pinMap(ixy)
    DO ig = 1, ng
      ir = ig + (ixy-1) * ng + (izf - myzbf) * ng * nxy
      DO ibd = 1, 4
        Dtil = cuCMFD%PinXS(ixy_map, iz)%Dtil(ibd, ig) / Pin(ixy_map)%BDLength(ibd)
        Dhat = cuCMFD%PinXS(ixy_map, iz)%Dhat(ibd, ig) / Pin(ixy_map)%BDLength(ibd)
        val = (Dtil - Dhat)
        CALL pushCSR(L(ibd), val, ir, ir)
        val = -(Dtil + Dhat)
        ineigh = Pin(ixy_map)%neighidx(ibd)
        ineigh = PinMapRev(ineigh)
        IF(ineigh .EQ. 0) val = 0
        ic = ir + (ineigh - ixy) * ng
        CALL pushCSR(L(ibd), val, ir, ic)
      END DO
      DO ibd = 1, 4
        Dtil = TranDtil(ibd, ig, ixy_map, iz) / Pin(ixy_map)%BDLength(ibd)
        Dhat = TranDhat(ibd, ig, ixy_map, iz) / Pin(ixy_map)%BDLength(ibd)
        val = (Dtil - Dhat)
        CALL pushCSR(tranL(ibd), val, ir, ir)
        val = -(Dtil + Dhat)
        ineigh = Pin(ixy_map)%neighidx(ibd)
        ineigh = PinMapRev(ineigh)
        IF(ineigh .EQ. 0) val = 0
        ic = ir + (ineigh - ixy) * ng
        CALL pushCSR(tranL(ibd), val, ir, ic)
      END DO
    END DO
  END DO
END DO

DO ibd = 1, nsurf
  CALL finalizeSortCSR(L(ibd), .TRUE.)
  CALL finalizeSortCSR(tranL(ibd), .TRUE.)
END DO

DO ibd = 1, nsurf
  nr = L(ibd)%nr
  nc = L(ibd)%nc
  nnz = L(ibd)%nnz
  csrVal => L(ibd)%d_csrVal
  csrRowPtr => L(ibd)%d_csrRowPtr
  csrColIdx => L(ibd)%d_csrColIdx
  ierr = cuSparseDcsrmv(cuDevice%mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, 1.0_8, &
                        L(ibd)%descr, csrVal, csrRowPtr, csrColIdx, cuCMFD%phis8, 0.0_8, Leakage(:,ibd))
END DO

DO ibd = 1, nsurf
  nr = tranL(ibd)%nr
  nc = tranL(ibd)%nc
  nnz = tranL(ibd)%nnz
  csrVal => tranL(ibd)%d_csrVal
  csrRowPtr => tranL(ibd)%d_csrRowPtr
  csrColIdx => tranL(ibd)%d_csrColIdx
  ierr = cuSparseDcsrmv(cuDevice%mySparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, nr, nc, nnz, 1.0_8, &
                        tranL(ibd)%descr, csrVal, csrRowPtr, csrColIdx, cuCMFD%phis8, 0.0, tranLeakage(:,ibd))
END DO

n = nsurf * ng * nxy * nzCMFD
CALL cuVectorOp('-', n, Leakage, tranLeakage, Ldiff, cuDevice%myStream)

hLeakage = Leakage
hLdiff = Ldiff
n = ng * nxy * nzCMFD
WRITE(mesg, '(a)') 'Leakage Difference Summary...'
IF(PE%master) CALL message(io8, TRUE, TRUE, MESG)
realn = n
CALL MPI_ALLREDUCE(realn, totn, 1, MPI_DOUBLE, MPI_SUM, MPI_CUDA_COMM, ierr)
DO ibd = 1, nsurf
  maxval = 0
  sumval = 0
  sumleak = 0
  DO i = 1, n
    val = abs(hLdiff(i, ibd))
    sumval = sumval + val
    sumleak = sumleak + abs(hLeakage(i, ibd))
    IF(val .GT. maxval) THEN
      maxval = val
    END IF
  END DO
  CALL MPI_ALLREDUCE(sumleak, totleak, 1, MPI_DOUBLE, MPI_SUM, MPI_CUDA_COMM, ierr)
  sumleak = totleak / totn
  CALL MPI_ALLREDUCE(sumval, val, 1, MPI_DOUBLE, MPI_SUM, MPI_CUDA_COMM, ierr)
  sumval = val / totleak
  CALL MPI_ALLREDUCE(maxval, val, 1, MPI_DOUBLE, MPI_MAX, MPI_CUDA_COMM, ierr)
  maxval = val / sumleak
  WRITE(mesg, '(1x,a6,a1,x,a10,es12.5,2x,a9,es12.5,2x,a9,es12.5)')  DirName(ibd), ':', 'Avg Leak =',sumleak, 'Avg Err =', sumval, 'Max Err =', maxval
  IF(PE%master) CALL message(io8, FALSE, TRUE, MESG)
END DO

DO ibd = 1, nsurf
  CALL DestroyCsr(L(ibd))
  CALL DestroyCSr(tranL(ibd))
END DO

DEALLOCATE(L, tranL)
DEALLOCATE(Leakage, TranLeakage, Ldiff, hLeakage, hLdiff)

!DO ig = 1, ng
!  DO iz = myzb, myze
!    DO ixy = 1, nxy
!      fmult = cuCMFD%h_phic8(ig, ixy, iz) / cuTranCMInfo%TranPhiC(ig, ixy, iz)
!      DO j = 1, Pin(ixy_map)%nxy
!        ipin = Pin(ixy_map)%pin(j)
!        FsrIdxSt = Core%Pin(ipin)%FsrIdxSt
!        icel = Core%Pin(ipin)%Cell(iz)
!        nLocalFsr = Core%CellInfo(icel)%nFsr
!        DO i = 1, nLocalFsr
!          ifsr = FsrIdxSt + i - 1
!          phis(ifsr, iz, :) = phis(ifsr, iz, :) * fmult
!        ENDDO
!        phic(ipin, iz, :) = phisum
!      ENDDO
!      cuCMFD%h_phic8(:, ixy, iz) = phisum
!    ENDDO
!  ENDDO
!END DO

111 CONTINUE

DO iz = myzb, myze
  DO ixy = 1, nxy
    TranDtil(:,:,ixy, iz) = cuCMFD%PinXS(ixy, iz)%Dtil(:,:)
    TranDhat(:,:,ixy, iz) = cuCMFD%PinXS(ixy, iz)%Dhat(:,:)
  END DO
END DO
NULLIFY(Pin, PinVolFm, hzfm, PlaneMap, PinMap, PinMapRev)


END SUBROUTINE

SUBROUTINE Transient_Corrector(TranInfo, TranCntl, nTracerCntl, Tbeg_inp, Tend_inp, Tdiv_inp)
USE PARAM
USE TYPEDEF,        ONLY : TranInfo_Type,     TranCntl_Type
USE CNTL,           ONLY : nTracerCntl_Type
USE TRAN_MOD,       ONLY : SetCorrectorPrecParam
USE PE_Mod,           ONLY : PE
USE files,            ONLY : io8
USE ioutil,           ONLY : message,           terminate
IMPLICIT NONE
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
REAL, OPTIONAL :: Tbeg_inp, Tend_inp, Tdiv_inp

REAL, POINTER :: T(:), delt(:)
REAL, POINTER :: omegam(:), omega0(:), omegap(:)
REAL, POINTER :: kappa(:), prec(:)
REAL, POINTER :: betak(:), betak_beg(:), betak_end(:)
REAL :: Tbeg, Tend, Tdiv, Tnow
REAL :: rho, beta, gamma, factor_F
REAL :: rho_beg, beta_beg, gamma_beg, factor_F_beg
REAL :: rho_end, beta_end, gamma_end, factor_F_end
REAL :: ressrc, precsrc, omegalm, omegal0, omegalp
REAL :: wt, wtbar
REAL :: theta
REAL :: denom, Invthdelt, coeff1
REAL :: nowdelt, prevdelt
REAL :: amp0
REAL :: alphamax, alphatilt, alphabeg, alphaend, tmod, p2dot, fmult
INTEGER :: nowstep, nstep
INTEGER :: nprec, norder
INTEGER :: istep, iprec

nowstep = TranCntl%nowstep
nprec = TranInfo%nprec


IF(PRESENT(Tbeg_inp)) THEN
  Tbeg = Tbeg_inp
  Tend = Tend_inp
  Tdiv = Tdiv_inp
ELSE
  Tend = TranCntl%T(nowstep)
  Tbeg = TranCntl%T(nowstep-1)
  Tdiv = TranCntl%Tdiv_corrector
END IF

ALLOCATE(betak(nprec), betak_beg(nprec), betak_end(nprec))

rho_beg = TranInfo%Prev_delrho * 1.E-5_8
beta_beg = TranInfo%Prev_corebetat
gamma_beg = TranInfo%Prev_lifetime
factor_F_beg = TranInfo%Prev_Factor_F / TranInfo%Amp
DO iprec = 1, nprec
  betak_beg(iprec) = TranInfo%Prev_corebeta(iprec)
END DO
rho_end = TranInfo%delrho_Dynamic * 1.E-5_8
beta_end = TranInfo%corebeta_Dynamic
gamma_end = TranInfo%lifetime_Dynamic
factor_F_end = TranInfo%factor_F / TranInfo%AmpPredictor
DO iprec = 1, nprec
  betak_end(iprec) = TranInfo%coreBetak(iprec)
END DO

theta = 1.
!theta = TranCntl%theta

WRITE(mesg, '(6x,a14,x, 3a12)'), 'PKE Parameters', 'RHO [pcm]',  'BETA [pcm]', 'GAMMA [pcm]'
IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
WRITE(mesg, '(6x,a15,x, 4es12.4)'), 'BEGIN: ', rho_beg*1.E5,  beta_beg*1.E5, gamma_beg*1.E5
IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
WRITE(mesg, '(6x,a15,x, 4es12.4)'), 'END  : ', rho_end*1.E5,  beta_end*1.E5, gamma_end*1.E5
IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)

ALLOCATE(omegam(nprec), omega0(nprec), omegap(nprec))
ALLOCATE(kappa(nprec), prec(nprec))

DO iprec = 1, nprec
  prec(iprec) = TranInfo%Prev_corePrec(iprec)
END DO

precsrc = 0.
DO iprec = 1, nprec
  precsrc = precsrc + TranInfo%lambda(iprec) * prec(iprec)
END DO
precsrc = precsrc * TranInfo%Inv_lifetime0
ressrc = (rho_beg - beta_beg) * TranInfo%Amp / gamma_beg + precsrc
TranInfo%TranAmp = TranInfo%Amp
amp0 = TranInfo%Amp

!--Adjust time-step size
IF(PRESENT(Tbeg_inp)) THEN
  alphabeg = (rho_beg - beta_beg) / gamma_beg
  alphaend = (rho_end - beta_end) / gamma_end
  alphatilt = (alphaend - alphabeg) / (Tend - Tbeg)
  IF(abs(alphabeg) .GT. abs(alphaend)) THEN
    alphamax = alphabeg
  ELSE
    alphamax = alphaend
  END IF
  p2dot = alphatilt + alphamax * (alphamax + precsrc / amp0)
  tmod = sqrt(1.E-5/abs(p2dot))
  fmult = tmod / Tdiv
  IF(fmult .GT. 1.5) fmult = 1.5
  IF(fmult .LT. 0.667) fmult = 0.667
  WRITE(mesg, '(6x,a15,x,es12.4,a8,es12.4,a3)'), 'PKE TIME-STEP: ', Tdiv*1000, ' ms --> ', fmult*Tdiv*1000, ' ms'
  IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
  Tdiv  = fmult * Tdiv
  Tdiv_inp = Tdiv
END IF

nstep = (Tend - Tbeg) / Tdiv + 1
ALLOCATE(T(nstep), delt(nstep))
nstep = 0
Tnow = Tbeg
DO WHILE(.TRUE.)
  nstep = nstep + 1
  Tnow = Tnow + Tdiv
  T(nstep) = Tnow
  delt(nstep) = Tdiv
  IF(Tnow .GT. Tend .OR. abs(Tend - Tnow) .LT. 1.E-6_8) THEN
    T(nstep) = Tend
    delt(nstep) = Tend - T(nstep-1)
    EXIT
  END IF
END DO


!WRITE(mesg, '(a10, 5es13.5)'), 'Amplitude:', TranInfo%Amp,  ressrc
!IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
DO istep = 1, nstep
  Tnow = T(istep)

  !-- interpolation of kinetic parameters
  wt = (Tnow - Tbeg) / (Tend - Tbeg)
  wtbar = 1. - wt
  rho = wtbar * rho_beg + wt * rho_end
  beta = wtbar * beta_beg  + wt * beta_end
  gamma = wtbar * gamma_beg + wt * gamma_end
  factor_F =  1./ (gamma * TranInfo%Inv_Factor_K0)
  DO iprec = 1, nprec
    betak(iprec) = wtbar * betak_beg(iprec) + wt * betak_end(iprec)
  END DO

  !-- Set Precursor Source
  nowdelt = delt(istep)
  IF(istep .EQ. 1 .OR. abs(TranCntl%theta - 0.5) .GT. epsm6) THEN
    prevdelt = nowdelt
    norder = 1
  ELSE
    prevdelt = delt(istep-1)
    norder = 2
  END IF
  CALL SetCorrectorPrecParam(TranInfo, omegam, omega0, omegap, betak, factor_F, nowdelt, prevdelt, nprec, norder)
  DO iprec = 1, nprec
    kappa(iprec) = exp(-delt(istep) * TranInfo%lambda(iprec))
  END DO
  precsrc = 0.
  omegalm = 0.
  omegal0 = 0.
  omegalp = 0.
  DO iprec = 1, nprec
    precsrc = precsrc + TranInfo%lambda(iprec) * kappa(iprec) * prec(iprec)
    omegalm = omegalm + TranInfo%lambda(iprec) * omegam(iprec)
    omegal0 = omegal0 + TranInfo%lambda(iprec) * omega0(iprec)
    omegalp = omegalp + TranInfo%lambda(iprec) * omegap(iprec)
  END DO
  precsrc = precsrc + omegalm * TranInfo%TranAmpd
  precsrc = precsrc + omegal0 * TranInfo%TranAmp
  precsrc = precsrc * TranInfo%Inv_lifetime0
  omegalp = omegalp * TranInfo%Inv_lifetime0

  !-- Calculate Amplitude
  Invthdelt = 1./ (theta * delt(istep))
  coeff1 = (rho - beta) / gamma
  denom = Invthdelt - coeff1 - omegalp
  TranInfo%Amp = Invthdelt * TranInfo%Amp + (1.-theta) / theta * ressrc + precsrc
  TranInfo%Amp = TranInfo%Amp / denom

  !-- Update ResSrc
  ressrc = (coeff1 + omegalp) * TranInfo%Amp + precsrc

  !-- Update Precursor
  DO iprec = 1, nprec
    prec(iprec) = kappa(iprec) * prec(iprec)
    prec(iprec) = prec(iprec) + omegam(iprec) * TranInfo%TranAmpd
    prec(iprec) = prec(iprec) + omega0(iprec) * TranInfo%TranAmp
    prec(iprec) = prec(iprec) + omegap(iprec) * TranInfo%Amp
  END DO

  !WRITE(mesg, '(a10, 5es13.5)'), 'Amplitude:', TranInfo%Amp, ressrc, tnow, nowdelt, prevdelt
  !IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)

  TranInfo%TranAmpd = TranInfo%TranAmp
  TranInfo%TranAmp = TranInfo%Amp
END DO
TranInfo%AmpTilt = ressrc / TranInfo%Amp
TranInfo%AmpRatio = TranInfo%Amp / amp0

DO iprec = 1, nprec
  TranInfo%PrecRatio(iprec) = prec(iprec) / TranInfo%Prev_CorePrec(iprec)
END DO

DEALLOCATE(omegam, omega0, omegap)
DEALLOCATE(kappa, prec)
DEALLOCATE(betak, betak_beg, betak_end)
DEALLOCATE(T, delt)

END SUBROUTINE

SUBROUTINE GetShape(TranInfo)
USE TYPEDEF,        ONLY : TranInfo_Type
IMPLICIT NONE
TYPE(TranInfo_Type) :: TranInfo

REAL :: invamp
INTEGER :: ng, nxy, nzCMFD
INTEGER :: n, ierr

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
n = ng * nxy * nzCMFD
amp_now = TranInfo%lifetime_Dynamic * TranInfo%Factor_F * TranInfo%Inv_Factor_K0
invamp = 1./amp_now

ierr = cublasDcopy_v2(cuDevice%myblasHandle, n, cuCMFD%phis8, 1, shape_now, 1)
ierr = cublasDscal_v2(cuDevice%myblasHandle, n, invamp, shape_now, 1)

END SUBROUTINE

SUBROUTINE SaveShape()
IMPLICIT NONE

INTEGER :: ng, nxy, nzCMFD
INTEGER :: n, ierr

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
n = ng * nxy * nzCMFD

ierr = cublasDcopy_v2(cuDevice%myblasHandle, n, shape_now, 1, shape_prev, 1)
amp_prev = amp_now

END SUBROUTINE

FUNCTION ShapeChange()
IMPLICIT NONE
REAL :: ShapeChange
REAL(8), POINTER, DEVICE :: diffvec(:)
REAL :: shapeavg
INTEGER :: ng, nxy, nzCMFD
INTEGER :: n, ierr

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
n = ng * nxy * nzCMFD

ALLOCATE(diffvec(n))
CALL cuVectorOp('-', n, shape_now, shape_prev, diffvec, cuDevice%myStream)
ShapeChange = dotMulti(cuCMFD_adj%phis8, diffvec, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
shapeavg = dotMulti(cuCMFD_adj%phis8, shape_now, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
ShapeChange = ShapeChange / shapeavg
DEALLOCATE(diffvec)
END FUNCTION

FUNCTION ShapeChange2()
USE CUDA_MASTER
IMPLICIT NONE
REAL :: ShapeChange2
REAL(8), POINTER, DEVICE :: diffvec(:)
REAL :: shapeavg
INTEGER :: ng, nxy, nzCMFD
INTEGER :: n, ierr

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
n = ng * nxy * nzCMFD

ALLOCATE(diffvec(n))
CALL cuVectorOp('-', n, shape_now, shape_prev, diffvec, cuDevice%myStream)
ShapeChange2 = normMulti(diffvec, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
shapeavg = normMulti(shape_now, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
ShapeChange2 = ShapeChange2 / shapeavg
DEALLOCATE(diffvec)
END FUNCTION


SUBROUTINE InterpolateShape(tnow, tbeg, tend, amp)
IMPLICIT NONE
REAL :: tnow, tbeg, tend, amp

REAL :: wt, wtbar
INTEGER :: ng, nxy, nzCMFD
INTEGER :: n, ierr

wt = (tend - tnow) / (tend - tbeg)
wtbar = 1._8 - wt

wt = wt * amp
wtbar = wtbar * amp

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
n = ng * nxy * nzCMFD
ierr = cublasDcopy_v2(cuDevice%myblasHandle, n, shape_prev, 1, cuCMFD%Phis8, 1)
ierr = cublasDscal_v2(cuDevice%myblasHandle, n, wt, cuCMFD%phis8, 1)
ierr = cublasDaxpy_v2(cuDevice%myblasHandle, n, wtbar,  shape_now, 1, cuCMFD%phis8, 1)

 END SUBROUTINE

SUBROUTINE CorrectAmp(Core, CmInfo, FmInfo, GroupInfo, TranInfo, nTracerCntl)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type, CmInfo_Type, FmInfo_Type, GroupInfo_Type, TranInfo_Type
USE CNTL,             ONLY : nTracerCntl_Type
USE CUDA_CMFD,        ONLY : cuCopyFlux, cuSetMOCPhis
USE PE_Mod,           ONLY : PE
USE files,            ONLY : io8
USE ioutil,           ONLY : message,           terminate
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(nTracerCntl_Type) :: nTracerCntl

REAL :: fmult
INTEGER :: ng, nxy, nzCMFD, myzb, myze, nfsr
INTEGER :: iz, ixy, ifsr
INTEGER :: n, ierr

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzCMFD = cuDevice%nzCMFD
nfsr = Core%nCoreFsr
n = ng * nxy * nzCMFD

myzb = cuDevice%myzb
myze = cuDevice%myze

fmult = TranInfo%amp / amp_now
ierr = cublasDscal_v2(cuDevice%myblasHandle, n, fmult, cuCMFD%phis8, 1)
WRITE(mesg, '(1x, a27, ES12.4)') 'PKE CORRECTION  : ', fmult-1.
IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)

CALL cuCopyFlux(cuCMFD, cuDevice, 2)

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO iz = myzb, myze
  DO ifsr = 1, nfsr
    FmInfo%phis(ifsr, iz, :) = FmInfo%phis(ifsr, iz, :) * fmult
    FmInfo%phim(:, :, ifsr, iz) = FmInfo%phim(:, :, ifsr, iz) * fmult
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO iz = myzb, myze
  DO ixy = 1, Core%nxy
    CmInfo%phic(ixy, iz, :) = CmInfo%phic(ixy, iz, :) * fmult
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

CALL cuUpdtPsiNowStep(Core, FmInfo, GroupInfo, nTracerCntl, TranInfo%eigv0)
END SUBROUTINE

SUBROUTINE ScaleFmFlux(FmInfo, TranInfo)
USE TYPEDEF,        ONLY : FmInfo_Type,     TranInfo_Type
IMPLICIT NONE
TYPE(FmInfo_Type) :: FmInfo
TYPE(TranInfo_Type) :: TranInfo

REAL :: fmult
INTEGER :: nfsr, myzb, myze
INTEGER :: ifsr, iz

nfsr = cuGeometry%nfsr
myzb = cuDevice%myzb
myze = cuDevice%myze
fmult = TranInfo%amp / amp_now

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO iz = myzb, myze
  DO ifsr = 1, nfsr
    FmInfo%phis(ifsr, iz, :) = FmInfo%phis(ifsr, iz, :) * fmult
    FmInfo%phim(:, :, ifsr, iz) = FmInfo%phim(:, :, ifsr, iz) * fmult
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE CUDA_amplitude_updt(Core, CmInfo, FmInfo, GroupInfo, ThInfo, TranInfo, TranCntl, nTracerCntl, PE, tdop0, amp0, tprev, tnow, dtth0)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,    CmInfo_Type,     FmInfo_Type,     GroupInfo_Type,    &
                           ThInfo_Type,      TranInfo_Type,   TranCntl_Type,   PE_Type
USE CNTL,           ONLY : nTracerCntl_Type
USE Core_mod,       ONLY : eigv
USE files,          ONLY : io8
USE ioutil,         ONLY : message,           terminate,          ShowHbar2
USE SUBGRP_MOD,     ONLY : SubGrpEffXsGen,    FxrChiGen
USE MOC_COMMON,     ONLY : SetCoreMacXs
USE CMFD_COMMON,    ONLY : HomogenizeXS_Cusping,  SetRadialCoupling, HomogenizeKinParam,  SetCMFDPrecCoeff
USE CUDA_CMFD,      ONLY : cuCMFDPsiUpdt,     cuSetCoarsePhis,      cuCopyFlux
USE XsPerturb_mod,  ONLY : XsPerturbation,   InitXsPerturbation,    SaveXsPerturbation,   RecoverXsPerturbation
USE TH_Mod,         ONLY : TransientTH,      THFeedBack,      UpdtTdopErr
USE Boron_mod,      ONLY : SetBoronCoolant
USE TRAN_MOD,       ONLY : SavePrevPKEParameters,   RecoverPrevPKEParameters
USE CUDA_MASTER
#ifdef MPI_ENV
USE MPIComm_Mod,    ONLY : BCAST,               MPI_SYNC,           MPI_MAX_REAL
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
REAL, POINTER :: tdop0(:,:)
REAL :: amp0, tprev, tnow, dtth0

REAL, POINTER :: PKEparamSave(:)
REAL :: tdopPrev(Core%nz, Core%nxy)
REAL :: amp0th, inv_eigv0
REAL :: tdopchg
REAL :: dtth, tth, tthd
REAL :: rtemp
INTEGER :: nthstep, nparam
INTEGER :: ith, iter
INTEGER :: nxy, myzb, myze, ng, nprec, n, ierr
LOGICAL :: lconv, laccept, lxslib, lprevUpdt, lsaved

nxy = cuGeometry%nxyc
myzb = cuDevice%myzb
myze = cuDevice%myze
ng = cuGeometry%ng
nprec = cuGeometry%nprec
lxslib = nTracerCntl%lxslib
inv_eigv0 = 1./ TranInfo%eigv0
nparam = 4 + 2*nprec

lsaved = .FALSE.
DO WHILE(.TRUE.)
  rtemp = FLOOR((tnow-tprev) * 1.E+5_8/ dtth0) * 1.E-5_8
  nthstep = CEILING(rtemp)
  dtth = (tnow-tprev) / REAL(nthstep)
  IF(nthstep .GT. 1) THEN
    lprevUpdt = .TRUE.
  ELSE
    lprevUpdt = .FALSE.
  END IF
  IF(lprevUpdt .AND. .NOT. lsaved) THEN
    ALLOCATE(PKEparamSave(nparam))
    CALL SavePrevPKEParameters(TranInfo, PKEparamSave, nparam)
    CALL SaveBaseTH(Core, THInfo)
    lsaved = .TRUE.
  END IF
  laccept = .TRUE.
  TranInfo%amp = amp0
  IF(nTracerCntl%lFeedback) tdopPrev(1:Core%nz, 1:Core%nxy) = tdop0(1:Core%nz, 1:Core%nxy)
  tth = tprev
  DO ith = 1, nthstep
    lconv = .FALSE.
    tthd = tth
    tth = tth + dtth
    amp0th = TranInfo%amp
    CALL XsPerturbation(Core, FmInfo, TranInfo, TranCntl, nTracerCntl, PE, tth, tthd)
    CALL SetBoronCoolant(Core, FmInfo%Fxr, nTracerCntl%boronppm, PE%myzb, PE%myze)
    CALL InterpolateShape(tth, tprev, tnow, amp0th)
    IF(nTracerCntl%lFeedback) THEN
      CALL cuUpdtPowerLevel(TranInfo, .FALSE.)
      nTracerCntl%PowerLevel = TranInfo%PowerLevel; ThInfo%PowLv = nTracerCntl%PowerLevel

      CALL cuCopyFlux(cuCMFD, cuDevice, 2)
      CALL cuSetCoarsePhis(Core, cuCMFD, cuDevice, CmInfo%Phic)
      WRITE(mesg, '(a, a,i3,a,i3,a,x,e11.4,a,e11.4,a)') 'Performing T/H Calculation...', '(', ith, '/', nthstep, ')', tthd, ' s, ', tth, ' s'
      IF(PE%Master) CALL MESSAGE(io8,TRUE,TRUE,mesg)
      CALL TransientTH(Core, CMInfo, FmInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE, dtth)
    END IF
    iter = 0
    DO WHILE(.TRUE.)
      iter = iter + 1
      IF(nTracerCntl%lFeedback) THEN
        CALL THFeedBack(Core,  CMInfo, FmInfo, ThInfo, nTracerCntl, PE)
#ifdef MPI_ENV
        CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif
      END IF
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
      CALL cuTranReactivityUpdt_Adj(GroupInfo, TranInfo, TranCntl, nTracerCntl%l3dim, .FALSE.)

      CALL Transient_Corrector(TranInfo, TranCntl,  nTracerCntl, tthd, tth, TranCntl%Tdiv_corrector)

      IF(nTracerCntl%lFeedback) THEN
        CALL cuUpdtPowerLevel(TranInfo, .FALSE.)
        TranInfo%PowerLevel = TranInfo%PowerLevel * TranInfo%amp / amp0th
        nTracerCntl%PowerLevel = TranInfo%PowerLevel; ThInfo%PowLv = nTracerCntl%PowerLevel

        WRITE(mesg, '(a, a,i3,a,i3,a,x,e11.4,a,e11.4,a)') 'Performing T/H Calculation...', '(', ith, '/', nthstep, ')', tthd, ' s, ', tth, ' s'
        IF(PE%Master) CALL MESSAGE(io8,TRUE,TRUE, mesg)
        CALL TransientTH(Core, CMInfo, FmInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE, dtth)
        IF(ThInfo%TdopChg .LT. 1.E-5) THEN
          lconv = .TRUE.
          CALL UpdtTdopErr(ThInfo, tdopchg, tdopPrev, Core%nxy, Core%nz)
        ELSE
          TranInfo%amp = amp0th
        END IF
      ELSE
        lconv = .TRUE.
      END IF
      IF(lconv) THEN
        EXIT
      END IF
    END DO
    IF(ith .NE. nthstep) THEN
      CALL cuSavePKEParameters(TranInfo, TranCntl, .TRUE.)
      CALL SaveTranTHsol(Core, ThInfo, TranCntl, nTracerCntl, PE)
    END IF
  END DO
  IF(lprevUpdt) THEN
    CALL RecoverPrevPKEParameters(TranInfo, PKEparamSave, nparam)
    CALL RecoverBaseTH(Core, THInfo)
  END IF
  IF(laccept) EXIT
END DO
IF(lsaved) DEALLOCATE(PKEparamSave)
IF(nTracerCntl%lFeedback) THEN
  TranCntl%dtth = dtth0
  TranCntl%nthstep = nthstep
END IF
PRINT*, TranInfo%amp, amp0th
END SUBROUTINE

SUBROUTINE PrintTranInfo(PE, TranInfo, TranCntl, ThInfo, nTracerCntl, ItrCntl, tnow)
USE PARAM
USE TYPEDEF,        ONLY : PE_Type, TranInfo_Type, ThInfo_Type, TranCntl_Type
USE CNTL,           ONLY : nTracerCntl_Type
USE itrcntl_mod,    ONLY : ItrCntl_Type
USE files,          ONLY : io8
USE ioutil,         ONLY : message,           terminate,          ShowHbar2
IMPLICIT NONE
TYPE(PE_Type) :: PE
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_TYpe) :: TranCntl
TYPE(ThInfo_Type) :: ThInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(ItrCntl_Type) :: ItrCntl
REAL :: tnow

IF(PE%MASTER) THEN
  IF(.NOT. nTracerCntl%lAdjoint) THEN
    WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, 1e11.4, x, a, 1e10.3, x, a, 1e12.5)') '@1', tnow,"s, P=", TranInfo%PowerLevel * 100, ', $= ', TranInfo%reactivity_Dynamic, ', G=', TranInfo%lifetime_Dynamic
    CALL message(io8, FALSE, TRUE, MESG)
    WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, F7.1, x, a, F8.5, x, a, 1e10.3)') '@2', tnow,"s, r(pcm)= ", TranInfo%delrho_Dynamic/10._8, ', k= ', TranInfo%TranEig_Dynamic/10._8, ', B= ', TranInfo%CoreBeta_Dynamic
    CALL message(io8, FALSE, TRUE, MESG)
  ELSE
    WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, 1e14.7, x, a, 1e14.7, x, a, 1e14.7)') '@3', tnow,"s, P=", TranInfo%PowerLevel * 100, ', $= ', TranInfo%reactivity_Dynamic, ', G=', TranInfo%lifetime_dynamic
    CALL message(io8, FALSE, TRUE, MESG)
    WRITE(mesg, '(a, 1p, e12.4, x, a, 1p, F7.1, x, a, F8.5, x, a, 1e14.7)') '@4', tnow,"s, r(pcm)= ", &
      TranInfo%delrho_Dynamic/10._8, ', k= ', TranInfo%TranEig_Dynamic/10._8, ', B= ', TranInfo%CoreBeta_Dynamic
    CALL message(io8, FALSE, TRUE, MESG)
  ENDIF
  IF(nTracerCntl%lFeedback) THEN
    WRITE(mesg, '(a, 1p, e12.4, x, a, f10.3, x, a, f10.3, x, a)') '@5', tnow,"s, Max Tf=  ", ThInfo%TfMax/10, "C, Max TfAvg= ", ThInfo%Tfavg_max/10, "C"
    CALL message(io8, FALSE, TRUE, MESG)
    WRITE(mesg, '(a, 1p, e12.4, x, a, f10.3, x, a, f10.3, x, a)') '@6', tnow,"s, Max Tcl= ", ThInfo%Tfcl_max/10, "C, Max TMOA=  ", ThInfo%TModoutAvg/10, "C"
    CALL message(io8, FALSE, TRUE, MESG)
  ENDIF
  WRITE(mesg, '(A, 1p, e12.4, x, a, i8)') '@9', tnow, "s, MOC = ", ItrCntl%Mocit - ItrCntl%Mocit0
  CALL message(io8, FALSE, TRUE, MESG)
  WRITE(mesg, '(A, 1p, e11.4, x, a, i8, x, a, i8)') '@10', tnow, "s, MGCMFD = ", ItrCntl%Cmfdit - ItrCntl%Cmfdit0, "Inner = ", ItrCntl%InnerIt
  CALL message(io8, FALSE, TRUE, MESG)
  IF(TranCntl%ladptt) THEN
    WRITE(mesg, '(A, 1p, e11.4, x, a, e11.4, x, a)') '@11', tnow, "s, DelT CMFD = ", TranCntl%delT(TranCntl%nowstep), "s"
    CALL message(io8, FALSE, TRUE, MESG)
    IF(nTracerCntl%lFeedback) THEN
      WRITE(mesg, '(A, 1p, e11.4, x, a, e11.4, x, a, i8)') '@12', tnow, "s, DelT TH =   ", TranCntl%dtth, "s, TH nstep = ", TranCntl%nthstep
      CALL message(io8, FALSE, TRUE, MESG)
    END IF
  END IF
  CALL ShowHbar2(io8)
END IF


END SUBROUTINE

FUNCTION MathErr(TranInfo, TranCntl)
USE TYPEDEF,        ONLY : TranInfo_Type,     TranCntl_Type
IMPLICIT NONE
REAL :: MathErr
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl

REAL(8), POINTER, DEVICE :: diffvec(:)
REAL :: fmult, delt, shapeavg
INTEGER :: ng, nxy, nzCMFD
INTEGER :: n, ierr

ng = cuGeometry%ng
nxy = cuGeometry%nxyc
nzCMFD = cuGeometry%nzCMFD
n = ng * nxy * nzCMFD
ALlOCATE(diffvec(n))

ierr = cublasDCopy_v2(cuDevice%myblasHandle, n, rtilda_now, 1, rtilda_prev, 1)
!rtilda_now = cuTranCMInfo%ResSrc
MathErr = normMulti(cuTranCMInfo%ResSrc, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)

CALL cuVectorOp('/', n, cuTranCMInfo%ResSrc, cuTranCMInfo%volInvVel, rtilda_now, cuDevice%myStream)
MathErr = normMulti(rtilda_prev, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
CALL cuVectorOp('*-', n, cuTranCMInfo%expo_alpha, cuCMFD%phis8, rtilda_now, cuDevice%myStream)

delt = TranCntl%delT(TranCntl%nowstep)
fmult = 1./TranInfo%amp
ierr = cublasDscal_v2(cuDevice%myblasHandle, n, fmult, rtilda_now, 1)
CALL cuVectorOp('-', n, rtilda_now, rtilda_prev, diffvec, cuDevice%myStream)
MathErr = normMulti(diffvec, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
MathErr = normMulti(diffvec, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
shapeavg = normMulti(shape_now, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
MathErr = MathErr * delt / shapeavg
DEALLOCATE(diffvec)

MathErr = normMulti(rtilda_now, n, MPI_CUDA_COMM, cuDevice%myblasHandle, .TRUE.)
END FUNCTION

SUBROUTINE TimeStepSelection(dt, err, errcrit, laccept)
IMPLICIT NONE
REAL :: dt
REAL :: err, errcrit
LOGICAL :: laccept

REAL :: safef = 0.90
REAL :: lowlim = 0.7
REAL :: highlim = 1.3
REAL :: corrfactor

corrfactor = safef * errcrit / err
IF(corrfactor .GT. highlim) THEN
  corrfactor = highlim
ELSE IF(corrfactor .LT. lowlim) THEN
  corrfactor = lowlim
END IF

dt = dt * corrfactor

IF(err .LE. errcrit) THEN
  laccept = .TRUE.
ELSE
  laccept = .FALSE.
END IF

END SUBROUTINE

END MODULE


#endif
