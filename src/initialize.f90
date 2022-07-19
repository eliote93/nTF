#include <defines.h>
SUBROUTINE initialize()

USE Geom,          ONLY : Core
USE Core_Mod,      ONLY : FmInfo, CmInfo, GroupInfo
USE RAYS,          ONLY : RayInfo
USE FILES,         ONLY : io8
USE ioutil,        ONLY : message
USE CNTL,          ONLY : nTracerCntl
USE VTK_Mod,       ONLY : ProcessVTK2DPlnModel
USE PE_MOD,        ONLY : PE
USE MPIConfig_Mod, ONLY : SetMPIEnv
USE LpShf_mod,     ONLY : LP_Shuffling
USE CNTLROD_mod,   ONLY : InitCntlRodConf, SetCrBankPosition
USE FXRVAR_MOD,    ONLY : lFXRVAR
USE XS_COMMON,     ONLY : AlignFxrs, AlignIsodata, AlignMLGdata, AlignResVarPin, AlignResIsoData
USE MOC_COMMON,    ONLY : AllocCoreMacXs
USE MPIComm_Mod,   ONLY : MPI_SYNC
USE MOC_MOD,       ONLY : initRT
#ifdef __INTEL_MKL
USE MKL_INIT
#endif

IMPLICIT NONE
! ----------------------------------------------------

IF (nTracerCntl%lHex) CALL HexMain

! Parallel Environment
#ifdef MPI_ENV
CALL SetMPIEnv(PE)
IF (PE%lidle) RETURN
#endif

! Geo. & Ray
IF (nTracerCntl%lHex) THEN
  CALL HexInit
  
  ! Control Rod Information
  IF (PE%RTMASTER .AND. nTracerCntl%lCrInfo) CALL InitCntlRodConf(Core, FmInfo, CMInfo, GroupInfo, nTracerCntl, PE)
ELSE
  IF (nTracerCntl%lnTIGRst) CALL CopyBaseCell

  CALL SetGeometries

  ! Control Rod Information
  IF (PE%RTMASTER .AND. nTracerCntl%lCrInfo) CALL InitCntlRodConf(Core, FmInfo, CMInfo, GroupInfo, nTracerCntl, PE)
  IF (nTracerCntl%lDcpl) CALL initDcpl

  ! Modular ray information
  CALL MakeRays
END IF

! Allocate Problem Dependent Memory
CALL AllocPDM

! Loading Patern
IF (nTracerCntl%LPShf .OR. nTracerCntl%lCooling) CALL LP_Shuffling(Core, FmInfo, GroupInfo, nTRACERCntl, PE)

! Caterogies Information of Resonanase Isotope and Eq. XS
IF (nTracerCntl%lXsLib) THEN
  mesg = 'Allocating XS...'
  IF (PE%master) CALL message(io8, TRUE, TRUE, mesg)
  
  IF (PE%RTmaster)        CALL AllocXsEQ
  IF (nTracerCntl%lDcpl)  CALL AllocDcplXsEq
  IF (nTracerCntl%lScat1) CALL P1ScatXsTreatment
END IF

! T/H
IF (PE%RTmaster .AND. nTracerCntl%lFeedBack .OR. lFXRVAR) CALL InitTH

! FXR Variation
IF (lFXRVAR .AND. .NOT.nTracerCntl%lBranch ) CALL FXRVARIATION(Core%ncorefxr, FmInfo, PE)

! Flux Variables
IF (PE%RTmaster) CALL SetCorePower(Core, nTracerCntl)
IF (PE%RTmaster) CALL InitFluxVariables

! RT
IF (PE%RTmaster) CALL initRT(RayInfo, Core, nTracerCntl, PE)

! Pre-generated Macro Cross Section
IF (nTracerCntl%lMacro) CALL AllocCoreMacXs(Core)

! 3D Solver Employing Intel MKL
#ifdef __INTEL_MKL
mesg = 'Allocating MKL Env...'
IF (PE%master) CALL message(io8, TRUE, TRUE, mesg)
IF (PE%lMKL.AND.mklCntl%lEnter) CALL SetMKLEnv(Core, FmInfo, RayInfo)
#endif

! GPU Acceleration
#ifdef __PGI
IF (PE%lCUDA .AND. nTracerCntl%lXsAlign) THEN
  CALL AlignIsodata(GroupInfo)
  CALL AlignFxrs(FmInfo%Fxr, GroupInfo)
  CALL AlignMLGdata
  CALL AlignResVarPin
  CALL AlignResIsoData(GroupInfo)
END IF

IF (PE%lCUDA) CALL CUDAInitialize(Core, RayInfo, FmInfo)
#endif

! Gamma Transport
#ifdef __GAMMA_TRANSPORT
IF (nTracerCntl%lGamma) CALL PrepareGamma(Core, RayInfo)
#endif

! VTK
IF(nTracerCntl%OutpCntl%VisMod .EQ. 4) THEN
  CALL ProcessVTK2DPlnMOdel(Core, FmInfo, GroupInfo, nTracerCntl, PE)
  CALL MPI_SYNC(PE%MPI_CMFD_COMM)
  
  STOP
END IF

! Detail Output
IF (PE%master .AND. nTracerCntl%lDetailOutput) CALL InputEdit(io8)
! ----------------------------------------------------

END SUBROUTINE initialize