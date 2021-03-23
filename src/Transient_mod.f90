#include <defines.h>
MODULE TRAN_MOD
USE PARAM
USE TYPEDEF,    ONLY : TranInfo_TYPE,       XsChange_TYPE,     TranCntl_TYPE, XSNoise_Type, XsCntlRod_Type
IMPLICIT NONE


TYPE(TranCntl_Type), TARGET :: TranCntl
TYPE(XsChange_Type), TARGET :: XsChange(100)
TYPE(XSNoise_Type), TARGET :: XsNoise(100)
TYPE(XsCntlRod_Type), TARGET :: XsCntlRod(100)

TYPE(TranInfo_Type) :: TranInfo


INTERFACE


SUBROUTINE TransientFsp_Driver(Core, FmInfo, CmInfo, TranInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,           ONLY : CoreInfo_Type,           FmInfo_Type,           CmInfo_Type,         &
                              TranInfo_Type,           GroupInfo_Type,        TranCntl_Type,       &
                              ThInfo_Type,             PE_Type
USE CNTL,              ONLY : nTracerCntl_Type

IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
END SUBROUTINE



SUBROUTINE InitTransient(Core, RayInfo, FmInfo, CmInfo, ThInfo, eigv, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,     FmInfo_Type,      CmInfo_Type,      TranInfo_Type,     &
                           RayInfo_Type,      GroupInfo_Type,   TranCntl_Type,    PE_TYPE,           &
                           ThInfo_Type
USE CNTL,           ONLY : nTracerCntl_Type
IMPLICIT NONE


TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE

REAL :: eigv
END SUBROUTINE

SUBROUTINE FluxNormTransient(Core, RayInfo, FmInfo, CmInfo, TranInfo, ng, nTracerCntl, PE)
!Flux level normalization
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,     FmInfo_Type,       CmInfo_Type,     &
                           RayInfo_Type,      TranInfo_Type,     PE_TYPE
USE CNTL,           ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(RayInfo_Type) :: RayInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
INTEGER :: ng

END SUBROUTINE

SUBROUTINE InitPrecursor(Core, FmInfo, CmInfo, TranInfo, GroupInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,     ONLY : CoreInfo_Type,       FmInfo_Type,           CmInfo_Type,        &
                        TranInfo_Type,       PE_Type,               GroupInfo_Type
USE CNTL,        ONLY : nTracerCntl_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

END SUBROUTINE
SUBROUTINE InitFmPrecursor(Core, FmInfo, TranInfo, nTracerCntl, PE, ng, nprec)
USE PARAM
USE TYPEDEF,           ONLY : CoreInfo_Type,       FmInfo_Type,         TranInfo_Type,       PE_Type,   &
                              FxrInfo_Type,        Pin_Type,            Cell_Type
USE CNTL,              ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
INTEGER :: ng, nprec

END SUBROUTINE

SUBROUTINE SetTimeStep(TranCntl)
USE PARAM
USE TYPEDEF,       ONLY : TranCntl_Type
IMPLICIT NONE

TYPE(TranCntl_Type) :: TranCntl

END SUBROUTINE

SUBROUTINE SetSamplingTimeStep(TranCntl)
USE PARAM
USE TYPEDEF,       ONLY : TranCntl_Type
IMPLICIT NONE

TYPE(TranCntl_Type) :: TranCntl

END SUBROUTINE

SUBROUTINE UpdtPowerLevel(Core, FmInfo, CmInfo, TranInfo, ng, nTracerCntl, PE)
! Flux level normalization
! normalize flux such that average flux in fuel region be unity
! then update fission source and moments accordingly
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,     FmInfo_Type,       CmInfo_Type,     &
                           TranInfo_Type,     PE_TYPE,                            &
                           PinXs_Type
USE CNTL,           ONLY : nTracerCntl_Type
USE BasicOperation, ONLY : MULTI_CA
#ifdef MPI_ENV
USE MPIComm_Mod,    ONLY : REDUCE
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
INTEGER :: ng

END SUBROUTINE

SUBROUTINE UpdtPrec(Core, FmInfo, CmInfo, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,             ONLY : CoreInfo_Type,          FmInfo_Type,           CmInfo_Type,         &
                                TranInfo_Type,          GroupInfo_Type,        TranCntl_Type,       &
                                PE_Type,                                                            &
                                FxrInfo_Type,           Pin_Type,               Cell_Type
USE CNTL,                ONLY : nTracerCntl_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

END SUBROUTINE

SUBROUTINE CmFmInfoSync(Core, FmInfo, CmInfo, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,             ONLY : CoreInfo_Type,          FmInfo_Type,           CmInfo_Type,         &
                                TranInfo_Type,          GroupInfo_Type,        TranCntl_Type,       &
                                PE_Type,                                                            &
                                FxrInfo_Type,           Pin_Type,               Cell_Type
USE CNTL,                ONLY : nTracerCntl_Type
USE MOC_MOD,             ONLY : PsiUpdate
USE BasicOperation,      ONLY : MULTI_CA
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

END SUBROUTINE

SUBROUTINE SaveTranSol(Core, FmInfo, CmInfo, TranInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,             ONLY : CoreInfo_Type,          FmInfo_Type,        CmInfo_Type,       &
                                TranInfo_Type,          TranCntl_Type,      PE_Type,           &
                                ThInfo_Type,            GroupInfo_Type
USE CNTL,                ONLY : nTracerCntl_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(THInfo_Type) :: ThInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

END SUBROUTINE

SUBROUTINE UpdtResSrc(Core, FmInfo, CmInfo, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,             ONLY : CoreInfo_Type,          FmInfo_Type,         CmInfo_Type,       &
                                TranInfo_Type,          GroupInfo_Type,      TranCntl_Type,     &
                                PE_Type
USE CNTL,                ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCNtl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE

END SUBROUTINE
!
!SUBROUTINE XsPerturbation(TranInfo, TranCntl, nTracerCntl)
!USE PARAM
!USE TYPEDEF,             ONLY : TranInfo_Type,          TranCntl_Type
!USE CNTL,                ONLY : nTracerCntl_Type
!IMPLICIT NONE
!TYPE(TranInfo_Type) :: TranInfo
!TYPE(TranCntl_Type) :: TranCntl
!TYPE(nTracerCntl_Type) :: nTracerCntl
!END SUBROUTINE


SUBROUTINE KinParamGen(Core, FmInfo, TranInfo, ThInfo,  GroupInfo,  lBetaUpdt, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,    FmInfo_Type,      ThInfo_Type,          &
                             TranInfo_Type,    GroupInfo_Type,   PE_Type,              &
                             FXRInfo_Type,     PinXs_Type,       Pin_Type,             &
                             PinInfo_Type,     Cell_Type
USE CNTL,             ONLY : nTracerCntl_Type
USE BenchXs,          ONLY : XsBaseBen,        DnpBetaBen,        NeutVeloBen
USE BasicOperation,   ONLY : CP_VA,            CP_CA,            MULTI_VA
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(ThInfo_Type) :: ThInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) ::  GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
LOGICAL :: lBetaUpdt
!
END SUBROUTINE


SUBROUTINE SolExpExtpltn(Core, FmInfo, CmInfo, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,                 ONLY : CoreInfo_Type,          FmInfo_Type,            CmInfo_Type,       &
                                    TranInfo_Type,          TranCntl_Type,          PE_Type,           &
                                    GroupInfo_Type
USE CNTL,                    ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

END SUBROUTINE

SUBROUTINE UpdtExpTrsf(Core, FmInfo, CmInfo, TranInfo, GroupInfo, lupdt, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,                 ONLY : CoreInfo_Type,          FmInfo_Type,            CmInfo_Type,       &
                                    TranInfo_Type,          TranCntl_Type,          PE_Type,           &
                                    GroupInfo_Type
USE CNTL,                    ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
LOGICAL :: lupdt

ENDSUBROUTINE

SUBROUTINE CheckCondMOC(Core, FmInfo, CmInfo, TranInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,            ONLY : CoreInfo_Type,         FmInfo_Type,         CmInfo_Type,    &
                               TranInfo_Type,         GroupInfo_Type,      TranCntl_Type,  &
                               ThInfo_Type,           PE_Type
USE CNTL,               ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(THInfo_Type) :: THInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

END SUBROUTINE

SUBROUTINE CheckSGFSP(Core, ThInfo, TranCntl, PE)
USE PARAM
USE TYPEDEF,            ONLY : CoreInfo_Type,   ThInfo_Type,    TranCntl_Type,    PE_TYPE
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(ThInfo_Type) :: ThInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(PE_Type) :: PE
END SUBROUTINE

SUBROUTINE CheckMOCUpdt(Core, FmInfo, CmInfo, GroupInfo, nTracerCntl, TranCntl, PE, lSave)
USE PARAM
USE TYPEDEF,            ONLY : CoreInfo_Type,     FmInfo_Type,         CmInfo_Type,          TranCntl_Type,  &
                               Pin_Type,          PinXS_Type,          FxrInfo_TYpe,         Cell_Type,      &
                               PE_Type,           GroupInfo_Type
USE CNTL,               ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(TranCntl_Type) :: TranCntl
TYPE(PE_Type) :: PE
LOGICAL :: lSave
END SUBROUTINE

SUBROUTINE UpdtBaseXsCondiMOC(Core, FmInfo, CmInfo, TranInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,            ONLY : CoreInfo_Type,         FmInfo_Type,         CmInfo_Type,    &
                               TranInfo_Type,         GroupInfo_Type,      TranCntl_Type,  &
                               ThInfo_Type,           PE_Type
USE CNTL,               ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(THInfo_Type) :: THInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

END SUBROUTINE

SUBROUTINE TranReactivityUpdt(CmInfo, eigv, TranCntl, PE, TranInfo)
USE PARAM
USE TYPEDEF,   ONLY : CmInfo_Type,  PE_TYPE, TranCntl_Type, TranInfo_Type
IMPLICIT NONE

TYPE(CmInfo_Type) :: CmInfo
REAL :: eigv
TYPE(PE_TYPE) :: PE
TYPE(TranCntl_Type) :: TranCntl
TYPE(TranInfo_Type) :: TranInfo
END SUBROUTINE

FUNCTION UpdtHomXsTr1g(PinXS, PhiC, ng)
USE PARAM
USE TYPEDEF,       ONLY : PinXS_Type
IMPLICIT NONE
TYPE(PinXS_Type) :: PinXS
REAL :: PhiC(ng)
REAL :: UpdtHomXsTr1g
INTEGER :: ng

END FUNCTION

FUNCTION UpdtSigA1g(Fxr, PhiFxr, ng, GroupInfo, nTracerCntl)
USE PARAM
USE TYPEDEF,              ONLY : FxrInfo_Type,         ThInfo_Type,        GroupInfo_Type
USE CNTL,                 ONLY : nTracerCntl_Type
IMPLICIT NONE

REAL :: UpdtSigA1g
TYPE(FxrInfo_Type) :: Fxr
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
REAL :: PhiFxr(ng)
INTEGER :: ng


END FUNCTION


FUNCTION UpdtThermalSigA(Fxr, PhiFxr, ng, GroupInfo, nTracerCntl)
USE PARAM
USE TYPEDEF,              ONLY : FxrInfo_Type,         ThInfo_Type,        GroupInfo_Type,       &
                                 XsMAc_Type
USE CNTL,                 ONLY : nTracerCntl_Type
USE MacXsLib_mod,         ONLY : MacXsBase
USE XsUtil_mod,           ONLY : GetXsMacDat,          ReturnXsMacDat
IMPLICIT NONE

REAL :: UpdtThermalSigA

TYPE(FxrInfo_Type) :: Fxr
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
REAL :: PhiFxr(ng)
INTEGER :: ng

END FUNCTION

FUNCTION UpdtSigA4g(Fxr, PhiFxr, ng, GroupInfo, nTracerCntl)
USE PARAM
USE TYPEDEF,              ONLY : FxrInfo_Type,         ThInfo_Type,        GroupInfo_Type,       &
                                 XsMAc_Type
USE CNTL,                 ONLY : nTracerCntl_Type
USE MacXsLib_mod,         ONLY : MacXsBase
USE XsUtil_mod,           ONLY : GetXsMacDat,          ReturnXsMacDat
IMPLICIT NONE

REAL :: UpdtSigA4g(4)

TYPE(FxrInfo_Type) :: Fxr
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
REAL :: PhiFxr(ng)
INTEGER :: ng

END FUNCTION

SUBROUTINE SetCorrectorPrecParam(TranInfo, omegam, omega0, omegap, beta, factor_F, delt, deltp, nprec, norder)
USE TYPEDEF,        ONLY : TranInfo_Type
IMPLICIT NONE
TYPE(TranInfo_Type) :: TranInfo
REAL :: omegam(nprec), omega0(nprec), omegap(nprec), beta(nprec)
REAL :: factor_F
REAL :: delt, deltp
INTEGER :: nprec, norder

END SUBROUTINE

SUBROUTINE CalcFisRate(FisRate, CmInfo, nxy, myzb, myze, ng)
USE TYPEDEF,          ONLY : CmInfo_Type,       PinXs_Type
IMPLICIT NONE
REAL :: FisRate(nxy, myzb:myze)
TYPE(CmInfo_Type) :: CmInfo
INTEGER :: nxy, myzb, myze, ng
END SUBROUTINE

END INTERFACE

  CONTAINS
  SUBROUTINE SavePrevPKEParameters(TranInfo, paramsave, nparam)
  USE TYPEDEF,      ONLY : TranInfo_Type
  IMPLICIT NONE
  TYPE(TranInfo_Type) :: TranInfo
  REAL :: paramsave(nparam)
  INTEGER, INTENT(IN) :: nparam

  INTEGER :: nprec
  INTEGER :: ibeg, iend

  nprec = TranInfo%nprec

  paramsave(1) = TranInfo%Prev_delrho
  paramsave(2) = TranInfo%Prev_corebetat
  paramsave(3) = TranInfo%Prev_lifetime
  paramsave(4) = TranInfo%Prev_factor_F
  ibeg = 5; iend = ibeg + nprec - 1
  paramsave(ibeg:iend) = TranInfo%Prev_corebeta(1:nprec)
  ibeg = iend + 1; iend = ibeg + nprec - 1
  paramsave(ibeg:iend) = TranInfo%Prev_coreprec(1:nprec)

  END SUBROUTINE

  SUBROUTINE RecoverPrevPKEParameters(TranInfo, paramsave, nparam)
  USE TYPEDEF,      ONLY : TranInfo_Type
  IMPLICIT NONE
  TYPE(TranInfo_Type) :: TranInfo
  REAL :: paramsave(nparam)
  INTEGER, INTENT(IN) :: nparam

  INTEGER :: nprec
  INTEGER :: ibeg, iend

  nprec = TranInfo%nprec

  TranInfo%Prev_delrho = paramsave(1)
  TranInfo%Prev_corebetat = paramsave(2)
  TranInfo%Prev_lifetime = paramsave(3)
  TranInfo%Prev_factor_F = paramsave(4)
  ibeg = 5; iend = ibeg + nprec - 1
  TranInfo%Prev_corebeta(1:nprec) = paramsave(ibeg:iend)
  ibeg = iend + 1; iend = ibeg + nprec - 1
  TranInfo%Prev_coreprec(1:nprec) = paramsave(ibeg:iend)

  END SUBROUTINE

  SUBROUTINE ShowTransientFlag(io, calclevel, tbeg, tend)
  USE PARAM,      ONLY : mesg
  IMPLICIT NONE
  INTEGER :: io
  INTEGER :: calclevel
  CHARACTER*17 :: CARD
  REAL :: tbeg, tend

  CHARACTER*12 :: hbar

  SELECT CASE(calclevel)
  CASE(0)
    card = 'MATERIAL UPDATE  '
  CASE(1)
    card = 'PKE UPDATE       '
  CASE(2)
    card = 'TH UPDATE        '
  CASE(3)
    card = 'CMFD UPDATE      '
  END SELECT


  hbar = '------------'
  WRITE(mesg,'(a12, i1,a2,a17,a,es12.4,a6,es12.4,a2,a12)') hbar, calclevel, ': ', card, '-', tbeg, ' s to ', tend, ' s', hbar

  print '(a)', mesg
  write(io, '(a)') mesg

  END SUBROUTINE

  SUBROUTINE EDIT_NNSAMPLING(Core, TranCntl, CmInfo, PE, ng)
  USE TYPEDEF,    ONLY : coreinfo_type, TranCntl_TYPE, CMInfo_Type, PE_TYPE, &
                         pin_type, cell_type
#ifdef MPI_ENV
  USE MPIComm_Mod,  ONLY : SENDRECV, MPI_SYNC
#endif
  IMPLICIT NONE
  TYPE(coreinfo_type) :: Core
  TYPE(TranCntl_TYPE) :: TranCntl
  TYPE(CMInfo_Type) :: CmInfo
  TYPE(PE_TYPE) :: PE
  INTEGER :: ng

  TYPE(pin_type), POINTER :: Pin(:)
  TYPE(cell_type), POINTER :: Cell(:)
  REAL, POINTER :: PhiC3D(:,:,:)
  INTEGER :: nxy, nz, myzb, myze
  INTEGER :: ixy, iz, icel, ig, iz1, iz2
  INTEGER :: nowstep, iSstep
  INTEGER :: Sio
  INTEGER :: iAsy, iAsyTyp
  INTEGER :: i

  nowstep = TranCntl%nowstep
  iSstep = TranCntl%iSstep
  IF(iSstep .GT. TranCntl%nSstep) RETURN
  IF(ABS(TranCntl%Ssteps(iSstep) - TranCntl%T(nowstep)) .GT. 1.E-5) THEN
    IF(TranCntl%T(nowstep) .GT. TranCntl%Ssteps(iSstep)) THEN
      STOP 'Noise Sampling Step is Skipped'
    ELSE
      RETURN
    END IF
  END IF

  Sio = TranCntl%Sio
  IF(iSstep .EQ. 1) THEN
    OPEN(Sio, FILE = 'NNSAMPLING', STATUS = 'replace')
  END IF
  nxy = Core%nxy
  nz = Core%nz
  Pin => Core%Pin
  Cell => Core%CellInfo
  myzb = PE%myzb
  myze = PE%myze

  IF(pe%master) THEN
  ALLOCATE(PhiC3D(ng,nxy,nz))
  ELSE
  ALLOCATE(PhiC3D(ng,nxy,myzb:myze))
  END IF
  PhiC3D = 1.
  !$OMP PARALLEL
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
  DO iz = myzb, myze
    DO ixy = 1, nxy
      DO ig = 1, ng
        phiC3D(ig,ixy,iz) = CmInfo%PhiC(ixy,iz,ig)
      END DO
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL

#ifdef MPI_ENV
  IF(PE%master) THEN
    DO i = 1, PE%nCmfdProc-1
      iz1 = PE%AxDomRange(1,i)
      iz2 = PE%AxDomRange(2,i)
      CALL SENDRECV(PhiC3D(1:ng,1:nxy,iz1:iz2), ng, nxy, iz2-iz1+1,i,.FALSE., PE%MPI_CMFD_COMM)
    END DO
  ELSE
      iz1 = myzb
      iz2 = myze
      CALL SENDRECV(PhiC3D(1:ng,1:nxy,iz1:iz2), ng, nxy, iz2-iz1+1, 0, .TRUE., PE%MPI_CMFD_COMM)
  END IF
  CALL MPI_SYNC(PE%MPI_CMFD_COMM)
#endif

  IF(PE%master) THEN
    WRITE(Sio, '(a15, i5)') 'SAMPLING_STEP',  iSstep
    DO iz = 1, nz
      IF(.not. core%lfuelplane(iz)) CYCLE
      DO ixy = 1, nxy
        icel = Pin(ixy)%Cell(iz)
        iAsy = Pin(ixy)%iasy
        iAsyTyp = Core%Asy(iAsy)%AsyType
        IF(.NOT. Core%AsyInfo(iAsyTyp)%lfuel) CYCLE
        WRITE(Sio,'(3i5, 2es14.6)', advance = 'NO') Pin(ixy)%ix, Pin(ixy)%iy, iz, Cell(icel)%Geom%lx, Cell(icel)%Geom%ly
        DO ig = 1, ng
          WRITE(Sio, '(es14.6)', advance = 'NO') PhiC3D(ig,ixy,iz)
        END DO
        WRITE(Sio, *)
      END DO
    END DO
  END IF

  DEALLOCATE(PhiC3D)
  NULLIFY(Pin, Cell)
  TranCntl%iSstep = iSstep + 1

  END SUBROUTINE

END MODULE
