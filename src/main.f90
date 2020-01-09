!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include <defines.h>
PROGRAM NTRACER
USE param
USE timer,            ONLY : nTracer_reset_dclock, nTracer_dclock, TimeChk
USE CNTL,             ONLY : nTracerCntl
USE PE_MOD,           ONLY : PE
USE MPIConfig_Mod,    ONLY : PEInitialize
USE GEOM,             ONLY : Core, ng
USE RAYS,             ONLY : RayInfo
USE Core_mod,         ONLY : CMInfo, FmInfo, THInfo, GroupInfo, eigv, peigv, xkconv, xkcrit
USE Restart_mod,      ONLY : WriteRestartFile
USE MCP_Util,         ONLY : lMCP_Restart
USE MCModular,        ONLY : MCSimulation, InitMCModular
USE InitMC_MOD,       ONLY : initializeMC
IMPLICIT NONE

REAL :: tinit
INTEGER :: lproblem
LOGICAL :: lfirst

CALL nTracer_reset_dclock()
TimeChk%TotalTime = nTracer_dclock(TRUE, FALSE)
TimeChk%TotalTime = nTracer_dclock(FALSE, FALSE)
#ifdef MPI_ENV
CALL PEInitialize()
#endif
CALL PreProcess()   !READ INPUT

lfirst = TRUE

DO WHILE(nTracerCntl%morecase)
  CALL readinput
  lproblem = nTracerCntl%lproblem
  IF (nTracerCntl%TRSolver .EQ. 2) THEN ! Transport solution from MC 
    CALL InitializeMC()
    CALL InitMCModular()
    CALL MCSimulation()
    CALL OutputEditMC()
  ELSE
    IF (lfirst) THEN
      CALL initialize()      !Prepare the calculation
    ENDIF
    SELECT CASE(lproblem)
      CASE(lsseigv)
        CALL sseig()
      CASE(lDcplsseigv)
        CALL DcplSseig()
      CASE(ldepletion)
        CALL Depletion_Driver()
      CASE(lTransient)
        CALL Transient_Driver()
      CASE(lCrCspGen)
        CALL CrCspGen_Driver()
      CASE(lXenonDynamics)
        CALL XeDyn_Driver()
      CASE(lBranch)
        CALL Branch_Driver()
      CASE(lEFTsearch)  ! Effective Fuel Temperature search
        CALL EFT_Driver()
    END SELECT
    IF (nTracerCntl%lWriteRst) CALL WriteRestartFile(Core, FmInfo, .FALSE., 0, nTracerCntl, PE)
    lfirst = FALSE
    nTracerCntl%CaseNo = nTracerCntl%CaseNo + 1
  ENDIF
ENDDO
TimeChk%TotalTime = nTracer_dclock(FALSE, FALSE) - TimeChk%TotalTime
CALL finalize()

END PROGRAM