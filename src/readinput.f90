#include <defines.h>
SUBROUTINE readinput
use param
use files,          only : filename,       InputFileIdx,  io5,            io8
use ioutil,         only : terminate,      toupper,       openfile,       IFnumeric,   &
                           nfields,        message,       ShowHbar2
use inputcards ,    only : oneline,        probe,         mxcard,         nblock,      &
                           FindBlockId,    FindCardId,    cards,          blocks
!USE CNTL,           ONLY : MASTER
USE PE_Mod,         ONLY : PE
USE ReadInpCard
USE CNTL,           ONLY : nTracerCntl
USE PE_MOD,         ONLY : PE
USE Geom,           ONLY : Core

#ifdef MPI_ENV
USE MPIComm_Mod,    ONLY : MPIWaitTurn,  MPI_SYNC
#endif
USE MCP_Util,         ONLY : lMCP_Restart
USE MCP_Restart_Main, ONLY : MCP_Initialize
implicit none
character(15)   :: cardname, blockname, astring
INTEGER         :: indev
INTEGER         :: idblock,idcard
LOGICAL :: master
master = PE%master

IF (master) THEN
  IF(nTracerCntl%CaseNo .GT. 0) CALL ShowHbar2(io8)

  WRITE(io8,'(//45x,a,a/a)') 'Echo of Input Deck ',trim(filename(5)),hbar1
ENDIF
#ifndef MPI_ASYNC_READ
#ifdef MPI_ENV
!Waiting while other CPU reads
CALL MPIWaitTurn(PE%MPI_COMM, PE%myrank, PE%nproc, .FALSE.)
#endif
#endif

IF (nTracerCntl%lnTIGRst) CALL ProcMatBin !--- 180531 JSR Edit

indev = io5
DO WHILE(TRUE)
  READ(io5,'(a512)',END=1000) oneline
  IF(probe .eq. DOT) nTracerCntl%morecase = .FALSE.
  IF(probe.eq.BANG) cycle;     IF(probe.eq.POUND) cycle
  IF(oneline.eq.BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe.eq.DOT) exit;       IF(probe.eq.SLASH) exit
  READ(oneline,*) blockname; CALL toupper(blockname)
  IF(BLOCKNAME .EQ. 'DEPL') nTracerCntl%lproblem = ldepletion
  IF(master) CALL message(io8,FALSE,FALSE,oneline)
  SELECT CASE(BLOCKNAME)
    CASE('STATE');        CALL ReadStateCard(indev, io8)
    CASE('GEOM')
      IF (nTracerCntl%lHex) THEN
        CALL HexReadInp(indev)
      ELSE
        CALL ReadGeomCard(indev, io8)
      END IF
    CASE('OPTION');       CALL ReadOptionCard(indev, io8)
    CASE('MATERIAL');     CALL ReadMaterialCard(indev, io8)
    CASE('TH');           CALL ReadTHCard(indev, io8)
    CASE('TRAN');         CALL ReadTranCard(indev, io8)
    CASE('DEPL');         CALL ReadDeplCard(indev, io8)
    CASE('XSEC');         CALL ReadXsecCard(indev, io8)
    CASE('DECOUPLE');     CALL ReadDcplCard(indev, io8)
    CASE('PARALLEL');     CALL ReadParallelCard(indev, io8)
    CASE('EDIT');         CALL ReadEditCard(indev, io8)
    CASE('LP_SHF');       CALL ReadLpShfCard(indev, io8)
    CASE('VISUAL');       CALL ReadVisualCard(indev, io8)
    CASE('CNTLROD');      CALL ReadCntlRodCard(indev, io8)
    CASE('CUSPING');      CALL ReadCuspingCard(indev,io8)
    CASE('XEDYN');        CALL ReadXeDynCard(indev, io8)
    CASE('MCP_RESTART');  CALL ReadMcpRestartCard(indev, io8)
    CASE('SUBCH_OP');     CALL ReadSubchOptionCard(indev, io8)
    CASE('MKL');          CALL ReadMKLCard(indev, io8)   !--- CNJ Edit : Intel MKL Option Parsing
    CASE('CUDA');         CALL ReadCUDACard(indev, io8)  !--- CNJ Edit : CUDA Option Parsing
    CASE('GCGEN');        CALL ReadGCCard(indev, io8)
    CASE('NTIG_RESTART'); CALL ReadNTIGCard(indev, io8) ! --- 180827 JSR
  END SELECT
ENDDO
1000 continue

!IF(nTracerCntl%ScatOd .GT. 0) THEN
!  IF(PE%nThread .GT. core%NXYA*core%NASYCELL) THEN
!    PE%nThread = core%NXYA*core%NASYCELL !nThread Bound for Pn Scat
!  ENDIF
!ENDIF

IF (MASTER) THEN
  write(mesg,'(a)') hbar2
  CALL message(io8,FALSE, FALSE, mesg)
  mesg = 'Reading Input from '//trim(filename(5))//'...'
  CALL message(io8, TRUE, TRUE, mesg)
    IF( nTracerCntl%TRSolver .EQ. 1 )THEN
    IF( nTracerCntl%lNodeMajor )THEN
        IF( nTracerCntl%lDomainDcmp )THEN
            IF( nTracerCntl%lLinSrcCASMO )THEN
                IF( nTracerCntl%lHybrid )THEN
                    mesg =  'SCHEME : NODE MAJOR / DCMP : T / LinSrc : T(Hybrid)'
                ELSE
                    mesg =  'SCHEME : NODE MAJOR / DCMP : T / LinSrc : T'
                ENDIF
            ELSE
                mesg =  'SCHEME : NODE MAJOR / DCMP : T / LinSrc : F'
            ENDIF
        ELSE
            IF( nTracerCntl%lLinSrcCASMO )THEN
                mesg =  'SCHEME : NODE MAJOR / DCMP : F / LinSrc : T'
            ELSE
                mesg =  'SCHEME : NODE MAJOR / DCMP : F / LinSrc : F'
            ENDIF
        ENDIF
        IF( nTracerCntl%ScatOd .GT. 0 )THEN
            IF ( .NOT. PE%lCUDA) THEN
                mesg =  'SCHEME : NODE MAJOR / DCMP : T / Pn AFSS : T'
            ELSE
                mesg =  'SCHEME : NODE MAJOR / DCMP : F / Pn AFSS : T'
            ENDIF
        ENDIF
    ELSE !GROUP MAJOR
        IF( nTracerCntl%lAFSS )THEN
            mesg =  'SCHEME : GROUP MAJOR / AFSS : T'
        ELSE
            mesg =  'SCHEME : GROUP MAJOR / AFSS : F'
        ENDIF
    ENDIF
    ELSE
        mesg =  'SCHEME : multigroup MC'
    ENDIF
    CALL message(io8,TRUE, TRUE, mesg)
ENDIF
#ifndef MPI_ASYNC_READ
#ifdef MPI_ENV
!waiting for other processor reading input
CALL MPIWaitTurn(PE%MPI_COMM, PE%myrank, PE%nproc, .TRUE.)
#endif
#else
CALL MPI_SYNC(PE%MPI_COMM)
#endif
!CLOSE(io5)

IF(nTracerCntl%lISOTXS) CALL PostProc_ISOTXS

IF(lMCP_Restart) CALL MCP_Initialize()

END SUBROUTINE

SUBROUTINE CheckMKLOptions()
USE PARAM
USE PE_MOD
USE FILES
USE IOUTIL
IMPLICIT NONE

#ifndef __INTEL_MKL
mesg = 'MKL Block Disabled: MKL Not Found'
CALL message(io8, TRUE, TRUE, mesg)
PE%lMKL = FALSE
RETURN
#else

#endif

END SUBROUTINE
