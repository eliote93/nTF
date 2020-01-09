#include <defines.h>
SUBROUTINE CrCspGen_Driver()
USE PARAM
USE TYPEDEF,         ONLY : 
USE GEOM,            ONLY : Core
USE Core_mod,        ONLY : CmInfo,        FmInfo,          GroupInfo
USE itrcntl_mod,     ONLY : ItrCntl
USE FILES,           ONLY : io8
USE IOUTIL,          ONLY : message
USE CNTL,            ONLY : nTracerCntl
USE PE_Mod,          ONLY : PE
USE CNTLROD_Mod,     ONLY : GetBankName
USE CrCspGen_MOD,    ONLY : CrCspGenCntl,                                      &
                            InitCrCspGen,  CrCspGenCmfd,    MovingCRPosition,  &
                            SaveCrCspFlux, WriteCrCspFile_ascii
#ifdef MPI_ENV
USE MPIComm_mod,    ONLY : MPI_SYNC
#endif                          

IMPLICIT NONE

INTEGER :: icrpos, iCrBeg, iCrEnd, ibank
REAL :: eigv

CHARACTER(10) :: BankName

LOGICAL :: MASTER

MASTER = PE%MASTER

IF(MASTER) THEN
  CALL Message(io8, TRUE, TRUE, 'Control Rod Cusping Function Generation ...')
ENDIF


iCrBeg = CrCspGenCntl%CrPln(1); iCrEnd = CrCspGenCntl%CrPln(2)
ibank = CrCspGenCntl%BankID; BankName = GetBankName(ibank)
CALL InitCrCspGen(core, GroupInfo, PE)
eigv=1
DO icrpos = iCrBeg, iCrEnd
  IF(MASTER) THEN
  	CALL Message(io8, TRUE, TRUE, 'Set Control Rod Position for Cusping Function Generation...')
    WRITE(mesg, '(5x, A10,A2,x, A, I5, x,A2,x, I5)') BankName,':', 'Plane', icrpos, 'to', iCrEnd
    CALL message(io8, FALSE, TRUE, mesg)
  	WRITE(MESG, '(A)')  
  ENDIF
  !PRINT *, 'CR POSITION', icrpos
  CALL MovingCRPosition(CmInfo, icrpos, PE)
  CALL CrCspGenCmfd(Core, CmInfo, eigv, PE, nTracerCntl, ItrCntl)
  IF(MASTER) CALL Message(io8, TRUE, TRUE, 'Save Control Rod Region Flux...')
  CALL SaveCrCspFlux(Core, CmInfo, Eigv, iBank, iCrPos, nTracerCntl, PE)
  IF(MASTER) THEN
    WRITE(mesg, '(A)') hbar2(1:77)
    CALL message(io8, FALSE, TRUE, MESG)
  ENDIF
ENDDO

IF(MASTER) THEN
  WRITE(mesg, '(A)') 'Writing Control Rod Cusping File...'
  CALL message(io8, FALSE, TRUE, MESG)
ENDIF
CALL WriteCrCspFile_ascii(Core, PE)
CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
END SUBROUTINE
