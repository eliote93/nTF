#include <defines.h>
SUBROUTINE initialize()
USE PARAM
USE Geom,           ONLY : Core
USE Core_Mod,       ONLY : FmInfo,      CmInfo,                   GroupInfo
USE RAYS,           ONLY : RayInfo
USE FILES,          ONLY : IO8
USE CNTL,           ONLY : nTracerCntl
USE DcplCore_Mod,   ONLY : DcplInfo
USE VTK_Mod,        ONLY : ProcessVTK
USE PE_MOD,         ONLY : PE,           DcplPE
USE MPIConfig_Mod,  ONLY : SetMPIEnv, SetGeomPEVariables
USE LpShf_mod,      ONLY : LP_Shuffling
USE CNTLROD_mod,    ONLY : InitCntlRodConf,                                           &
                           SetCrBankPosition
USE SubChGen_mod,   ONLY : ConfigSubChGeom  !---BYS edit : off SubChGen Mod ??
USE FXRVAR_MOD,     ONLY : lFXRVAR
#ifdef __INTEL_MKL
USE MKL_INIT
#endif
USE MOC_COMMON,     ONLY : AllocCoreMacXs
IMPLICIT NONE
INTEGER :: i
!USE CNTL, ONLY : lXsLib

IF (nTracerCntl%lHex) CALL HexMain()

!Set Up Parallel Enviorment
#ifdef MPI_ENV
CALL SetMPIEnv(PE)
IF(PE%lidle) RETURN
#endif
!CALL SetGeomPEVariables(PE)

IF (nTracerCntl%lHex) THEN
  CALL HexInit()
ELSE
  IF (nTracerCntl%lnTIGRst) CALL CopyBaseCell()
  
  CALL SetGeometries()
  
  IF(nTracerCntl%lDcpl) CALL initDcpl()
  
  !Generate modular ray information
  CALL MakeRays()
END IF


!Allocate Problem Dependent Memory
CALL AllocPDM()

!Loading Patern
IF(nTracerCntl%LPShf .OR. nTracerCntl%lCooling) CALL LP_Shuffling(Core, FmInfo, GroupInfo,  nTRACERCntl, PE)

!Control Rod Position
!Set Control Rod Information
IF(PE%RTMASTER) CALL InitCntlRodConf(Core, FmInfo, CMInfo, GroupInfo, nTracerCntl, PE)

!IF(nTracerCntl%lCrInfo) CALL SetCrBankPosition(Core, FmInfo, CmInfo, GroupInfo, nTracerCntl, PE)

!Initialize and Set up Caterogies information of resonanase isotope and Eq. XS
IF(nTracerCntl%lXsLib) THEN
  IF(PE%RTmaster) CALL AllocXsEQ()
  IF(nTracerCntl%lDcpl) CALL AllocDcplXsEq()
  IF(nTracerCntl%lScat1) CALL P1ScatXsTreatment()
ENDIF

!CALL ConfigSubChGeom(Core, nTracerCntl, PE)
IF(PE%RTmaster .AND. nTracerCntl%lFeedBack .OR. lFXRVAR) CALL InitTH()

! ---BYS edit !!!!!!!!
IF( lFXRVAR .AND. .NOT. nTracerCntl%lBranch ) CALL FXRVARIATION(Core%ncorefxr,FmInfo,PE)

!Initialize the flux variables
IF(PE%RTmaster) CALL SetCorePower(Core, nTracerCntl)
IF(PE%RTmaster) CALL InitFluxVariables()
IF(nTRACERCNTL%lReadXeEQ) THEN
  !IF(PE%RTmaster) CALL ReadXeND(Core, FmInfo, PE)
ENDIF

!--- CNJ Edit : Pre-generated Macro Cross Section
IF (nTracerCntl%lMacro) CALL AllocCoreMacXs(Core)

!--- CNJ Edit : 3D Solver Employing Intel MKL
#ifdef __INTEL_MKL
IF (PE%lMKL) CALL SetMKLEnv(Core, FmInfo, RayInfo)
#endif

!--- CNJ Edit : GPU Acceleration
#ifdef __PGI
IF (PE%lCUDA) CALL CUDAInitialize(Core, RayInfo, FmInfo)
#endif

!--- CNJ Edit : Gamma Transport Initialization
#ifdef __GAMMA_TRANSPORT
IF (nTracerCntl%lGamma) CALL PrepareGamma(Core, RayInfo)
#endif

!Input Edit - writing output file
!CALL ProcessVTK(Core, FmInfo, GroupInfo%ng, PE)
IF(PE%master .AND. nTracerCntl%lDetailOutput) CALL InputEdit(io8)

END SUBROUTINE
