#include <defines.h>
MODULE InitMC_MOD
    
INTEGER :: ninact, nact, nht       !MC temporal variables
LOGICAL :: lMCcmfd, lfdb(2),lMCcmfdset          !MC temporal variables
INTEGER :: idx, period, skip, flag !MC temporal variables
CONTAINS

SUBROUTINE initializeMC()
USE PARAM
USE Geom,           ONLY : Core, AsyInfo, CellInfo, PinInfo
USE Core_Mod,       ONLY : FmInfo,      CmInfo,                   GroupInfo
USE RAYS,           ONLY : RayInfo
USE FILES,          ONLY : IO8
USE CNTL,           ONLY : nTracerCntl
USE DcplCore_Mod,   ONLY : DcplInfo
USE VTK_Mod,        ONLY : ProcessVTK
USE PE_MOD,         ONLY : PE,           DcplPE
USE MPIConfig_Mod,  ONLY : SetMPIEnv
USE LpShf_mod,      ONLY : LP_Shuffling
USE CNTLROD_mod,    ONLY : InitCntlRodConf,                                           &
                           SetCrBankPosition
USE SubChGen_mod,   ONLY : ConfigSubChGeom  !---BYS edit : off SubChGen Mod ??
USE Files,     ONLY : io8
USE IOUTIL,    ONLY : message
IMPLICIT NONE
INTEGER :: i
!USE CNTL, ONLY : lXsLib

!Set Up Parallel Enviorment
#ifdef MPI_ENV
CALL SetMPIEnv(PE)
IF(PE%lidle) RETURN
#endif
!CALL SetGeomPEVariables(PE)

CALL SetGeometries()

IF(nTracerCntl%lDcpl) CALL initDcpl()


! MC Debugging
!nTracerCntl.TRSolver=2 ! set to MC
  IF( nTracerCntl%lrestrmt .AND. nTracerCntl%lxslib )THEN
  IF(PE%Master) THEN
    mesg = 'Initialize for MOC variables to using nTRACER mic-library in MC'
    CALL message(io8,TRUE, TRUE, mesg)
  ENDIF
    !Generate modular ray information
    CALL MakeRays()
    !Allocate Problem Dependent Memory
    CALL AllocPDM()
    CALL FXRVARIATION(Core%ncorefxr,FmInfo)
    IF(PE%RTmaster) CALL AllocXsEQ()
    IF(nTracerCntl%lDcpl) CALL AllocDcplXsEq()
    IF(nTracerCntl%lScat1) CALL P1ScatXsTreatment()
    IF(PE%RTmaster .AND. nTracerCntl%lFeedBack) CALL InitTH()
    IF(PE%RTmaster) CALL SetCorePower(Core, nTracerCntl)
    IF(PE%RTmaster) CALL InitFluxVariables()
  ENDIF  
ENDSUBROUTINE
ENDMODULE
