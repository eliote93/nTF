#include <defines.h>
MODULE PromptFeedback_mod
IMPLICIT NONE

!REAL, PRIVATE, PARAMETER :: zeroPow = 1.E-6 ! W/cm3
REAL, PRIVATE, PARAMETER :: epsilon = 3.204E-11 ! J/fission
REAL, PRIVATE, PARAMETER :: alpha = 3.83E-11 ! K*cm3
REAL, PRIVATE :: zeroNormFactor
REAL, PRIVATE :: alphaEff

REAL, POINTER, PRIVATE :: fisRate0(:,:)
REAL, POINTER, PRIVATE :: fisRate(:,:)
REAL, POINTER, PRIVATE :: prevT(:,:)
REAL, POINTER, PRIVATE :: fvolratio(:,:)

  CONTAINS
  SUBROUTINE InitPromptFeedback(CoreInfo, FmInfo, CmInfo, TranInfo, PE, ng)
  USE TYPEDEF,              ONLY : CoreInfo_Type,  FmInfo_Type,    CmInfo_Type,  TranInfo_Type,  &
                                   PE_Type
  USE BenchXs,              ONLY : xsnfDynBen
  IMPLICIT NONE
  TYPE(CoreInfo_Type) :: CoreInfo
  TYPE(FmInfo_Type) :: FmInfo
  TYPE(CmInfo_Type) :: CmInfo
  TYPE(TranInfo_Type) :: TranInfo
  TYPE(PE_Type) :: PE
  INTEGER :: ng

  REAL, POINTER :: XsMacNf(:)
  REAL :: totVol, PwSum0
  REAL :: PwLevel
  REAL :: celvol, fvol
  INTEGER :: nxy, myzb, myze
  INTEGER :: FxrIdxSt, nLocalFxr
  INTEGER :: ixy, iz, icel, ifxr, i, imix

  totVol = CoreInfo%totVol
  PwSum0 = TranInfo%PwSum0
  nxy = CoreInfo%nxy
  myzb = PE%myzb
  myze = PE%myze

  ALLOCATE(fisRate0(nxy, myzb:myze), fisRate(nxy, myzb:myze))
  ALLOCATE(fvolratio(nxy, myzb:myze))
  ALLOCATE(prevT(nxy, myzb:myze))
  ALLOCATE(XsMacNf(ng))
  fvolratio = 0.

  PwLevel = epsilon * PwSum0 / totVol
  zeroNormFactor = TranInfo%InitPow / PwLevel
  alphaEff = alpha * zeroNormFactor

  print'(4es14.6)', alphaEff, TranInfo%InitPow, PwLevel, totVol

  CALL UpdtFisRate(CmInfo, PE, nxy, myzb, myze, ng, .TRUE.)
  DO iz = myzb, myze
    DO ixy = 1, nxy
      IF(.NOT. CoreInfo%Pin(ixy)%lfuel) CYCLE
      icel = CoreInfo%Pin(ixy)%Cell(iz)
      FxrIdxSt = CoreInfo%Pin(ixy)%FxrIdxSt
      nLocalFxr = CoreInfo%CellInfo(icel)%nFxr
      celvol = 0.
      fvol = 0.
      DO i = 1, nLocalFxr
        ifxr = FxrIdxSt + i - 1
        imix = FmInfo%Fxr(ifxr, iz)%imix
        celvol = celvol + FmInfo%Fxr(ifxr, iz)%area
        CALL xsnfDynBen(imix, TranInfo%fuelTemp(ixy, iz), 1, ng, XsMacNf)
        IF(SUM(XsMacNF) .GT. 0) THEN
          fvol = fvol + FmInfo%Fxr(ifxr, iz)%area
        END IF
      END DO
      fvolratio(ixy, iz) = celvol / fvol
    END DO
  END DO
  DEALLOCATE(XsMacNf)

  END SUBROUTINE

  SUBROUTINE UpdtFisRate(CmInfo, PE, nxy, myzb, myze, ng, lFirst)
  USE TYPEDEF,          ONLY : CmInfo_Type,       PE_Type
  USE TRAN_MOD,         ONLY : CalcFisRate
#ifdef __PGI
  USE CUDA_Transient,   ONLY : cuCalcFisRate
#endif
  IMPLICIT NONE
  TYPE(CmInfo_Type) :: CmInfo
  TYPE(PE_Type) :: PE
  INTEGER :: nxy, myzb, myze, ng
  LOGICAL :: lFirst

#ifdef __PGI
  IF(PE%lCUDA) THEN
    IF(lFirst) THEN
      CALL cuCalcFisRate(fisRate0, nxy, myzb, myze, ng)
    ELSE
      CALL cuCalcFisRate(fisRate, nxy, myzb, myze, ng)
    END IF
    GOTO 1
  END IF
#endif
  IF(lFirst) THEN
    CALL CalcFisRate(fisRate0, CmInfo, nxy, myzb, myze, ng)
  ELSE
    CALL CalcFisRate(fisRate, CmInfo, nxy, myzb, myze, ng)
  END IF
1 CONTINUE

  END SUBROUTINE

  SUBROUTINE UpdtFuelTemp(CoreInfo, CmInfo, TranInfo, TranCntl, PE, ng, lprevupdt)
  USE PARAM
  USE TYPEDEF,              ONLY : CoreInfo_Type,  CmInfo_Type,  TranInfo_Type,  &
                                   TranCntl_Type,  PE_Type
  USE files,                ONLY : io8
  USE ioutil,               ONLY : message
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  TYPE(CoreInfo_Type) :: CoreInfo
  TYPE(CmInfo_Type) :: CmInfo
  TYPE(TranInfo_Type) :: TranInfo
  TYPE(TranCntl_Type) :: TranCntl
  TYPE(PE_Type) :: PE
  INTEGER :: ng
  LOGICAL :: lprevupdt

  REAL :: delt, frdiff
  REAL :: tfmax, tfavg, tmpavg, tmpmax
  INTEGER :: nowstep
  INTEGER :: nxy, myzb, myze
  INTEGER :: ixy, iz
  INTEGER :: ierr

  REAL :: tfhfp, w_t, lambda_t
  LOGICAL :: lHFPmodel
  logical, save :: lfirst
  data lfirst /.true./

  nxy = CoreInfo%nxy
  myzb = PE%myzb
  myze = PE%myze

  nowstep = TranCntl%nowstep
  delt = Trancntl%delt(nowstep)

  lHFPmodel = .TRUE.
  lambda_t = 1
  w_t = exp(-lambda_t*TranCntl%T(nowstep))

  !IF(lprevupdt) THEN
  IF(lfirst) THEN
    prevT(:,myzb:myze) = TranInfo%fuelTemp(:,myzb:myze)
    lfirst = .FALSE.
  END IF

  CALL UpdtFisRate(CmInfo, PE, nxy, myzb, myze, ng, .FALSE.)

  DO iz = myzb, myze
    IF(.NOT. CoreInfo%lFuelPlane(iz)) CYCLE
    !$OMP PARALLEL PRIVATE(frdiff, tfhfp)
    !$OMP DO SCHEDULE(GUIDED)
    DO ixy = 1, nxy
      IF(.NOT. CoreInfo%Pin(ixy)%lfuel) CYCLE
      frdiff = fisRate(ixy, iz) - fisRate0(ixy, iz)
      TranInfo%fuelTemp(ixy, iz) = prevT(ixy, iz) + delt * alphaEff * frdiff * fvolratio(ixy, iz)
      IF(lHFPmodel) THEN
        tfhfp = 300 + 600 * TranInfo%PowerLevel
        IF(.not. lfirst) prevT(ixy,iz) = TranInfo%fuelTemp(ixy,iz)
        TranInfo%fuelTemp(ixy,iz) = w_t * TranInfo%fuelTemp(ixy,iz) + (1-w_t)*tfhfp
      END IF
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
  END DO

  tfmax = 0.
  tfavg = 0.
  DO iz = myzb, myze
    IF(.NOT. CoreInfo%lFuelPlane(iz)) CYCLE
    DO ixy = 1, nxy
      IF(.NOT. CoreInfo%Pin(ixy)%lfuel) CYCLE
      tfmax = max(TranInfo%fuelTemp(ixy, iz), tfmax)
      tfavg = tfavg + TranInfo%fuelTemp(ixy, iz)
    END DO
  END DO
#ifdef MPI_ENV
  CALL MPI_ALLREDUCE(tfavg, tmpavg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, PE%MPI_CMFD_COMM, ierr)
  tfavg = tmpavg
  CALL MPI_ALLREDUCE(tfmax, tmpmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, PE%MPI_CMFD_COMM, ierr)
  tfmax = tmpmax
#endif
  tfavg = tfavg / TranInfo%nfuelcell

  WRITE(mesg,601) tfmax, tfavg
  IF(PE%master) CALL message(io8,false,TRUE,mesg)
601 format(11x,"Max Tf=",f8.2, " K,  Avg Tf",f8.2," K")

  END SUBROUTINE


END MODULE
