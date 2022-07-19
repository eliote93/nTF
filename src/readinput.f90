#include <defines.h>
SUBROUTINE readinput

USE ReadInpCard
USE param,      ONLY : TRUE, FALSE, mesg, hbar2
USE files,      ONLY : filename, InputFileIdx, io5, io8
USE ioutil,     ONLY : terminate, toupper, openfile, IFnumeric, nfields, message, ShowHbar2
USE inputcards, ONLY : oneline, probe, mxcard, nblock, FindBlockId, FindCardId, cards, blocks
USE CNTL,       ONLY : nTracerCntl
USE PE_MOD,     ONLY : PE
USE Geom,       ONLY : Core
USE RAYS,       ONLY : RayInfo

USE MCP_Util,         ONLY : lMCP_Restart
USE MCP_Restart_Main, ONLY : MCP_Initialize

#ifdef MPI_ENV
USE MPIComm_Mod,    ONLY : MPIWaitTurn, MPI_SYNC
#endif

IMPLICIT NONE

CHARACTER(15) :: cardname, blockname, astring, dumc
INTEGER :: indev, idblock,idcard
REAL :: del
LOGICAL :: master
! ----------------------------------------------------

master = PE%master

IF (master) THEN
  IF (nTracerCntl%CaseNo .GT. 0) CALL ShowHbar2(io8)
  
  WRITE (io8, '(//45X, A,A/A)') 'Echo of Input Deck ', trim(filename(5)), hbar1
END IF

#ifndef MPI_ASYNC_READ
#ifdef MPI_ENV
CALL MPIWaitTurn(PE%MPI_COMM, PE%myrank, PE%nproc, FALSE) ! Waiting while other CPU reads
#endif
#endif

IF (nTracerCntl%lnTIGRst) CALL ProcMatBin !--- 180531 JSR Edit
! ----------------------------------------------------
indev = io5

DO WHILE (TRUE)
  READ (io5,'(a512)',END = 1000) oneline
  
  IF (probe .EQ. DOT) nTracerCntl%morecase = FALSE
  
  IF (probe .EQ. BANG)     CYCLE
  IF (probe .EQ. POUND)    CYCLE
  IF (oneline .EQ. BLANK)  CYCLE
  IF (IFnumeric(oneline))  CYCLE
  
  IF (probe .EQ. DOT)   EXIT
  IF (probe .EQ. SLASH) EXIT
  
  IF (master) CALL message(io8, FALSE, FALSE, oneline)
  READ (oneline,*) blockname
  CALL toupper(blockname)
  
  IF (BLOCKNAME .EQ. 'DEPL') nTracerCntl%lproblem = ldepletion
  
  SELECT CASE(BLOCKNAME)
    CASE('STATE');        CALL ReadStateCard(indev, io8)
    CASE('GEOM')
      IF (nTracerCntl%lHex) THEN
        CALL HexReadInp(indev)
      ELSE
        CALL ReadGeomCard(indev, io8)
      END IF
    CASE('OPTION');       CALL ReadOptionCard     (indev, io8)
    CASE('MATERIAL');     CALL ReadMaterialCard   (indev, io8)
    CASE('TH');           CALL ReadTHCard         (indev, io8)
    CASE('TRAN');         CALL ReadTranCard       (indev, io8)
    CASE('DEPL');         CALL ReadDeplCard       (indev, io8)
    CASE('XSEC');         CALL ReadXsecCard       (indev, io8)
    CASE('DECOUPLE');     CALL ReadDcplCard       (indev, io8)
    CASE('PARALLEL');     CALL ReadParallelCard   (indev, io8)
    CASE('EDIT');         CALL ReadEditCard       (indev, io8)
    CASE('LP_SHF');       CALL ReadLpShfCard      (indev, io8)
    CASE('VISUAL');       CALL ReadVisualCard     (indev, io8)
    CASE('CNTLROD');      CALL ReadCntlRodCard    (indev, io8)
    CASE('CUSPING');      CALL ReadCuspingCard    (indev, io8)
    CASE('XEDYN');        CALL ReadXeDynCard      (indev, io8)
    CASE('MCP_RESTART');  CALL ReadMcpRestartCard (indev, io8)
    CASE('SUBCH_OP');     CALL ReadSubchOptionCard(indev, io8)
    CASE('MKL');          CALL ReadMKLCard        (indev, io8) !--- CNJ Edit : Intel MKL Option Parsing
    CASE('CUDA');         CALL ReadCUDACard       (indev, io8) !--- CNJ Edit : CUDA Option Parsing
    CASE('GCGEN');        CALL ReadGCCard         (indev, io8)
    CASE('NTIG_RESTART'); CALL ReadNTIGCard       (indev, io8) ! --- 180827 JSR
  END SELECT
END DO

1000 CONTINUE
! ----------------------------------------------------
IF (MASTER) THEN 
  WRITE (mesg, '(A)') hbar2 
  CALL message(io8, FALSE, FALSE, mesg) 
  
  mesg = 'Reading Input from ' // trim(filename(5)) // '...'
  CALL message(io8, TRUE, TRUE, mesg)
  
  IF (PE%nproc .GT. 1) THEN
    WRITE (dumc, '(I2)') PE%nproc
    mesg = 'Prc. # : ' // trim(dumc)
    CALL message(io8, TRUE, TRUE, mesg)
  END IF
  
  IF (PE%nThread .GT. 1) THEN
    WRITE (dumc, '(I2)') PE%nThread
    mesg = 'Thr. # : ' // trim(dumc)
    CALL message(io8, TRUE, TRUE, mesg)
  END IF
  
  del = RayInfo%Del
  
  IF (int(1000*del) .EQ. 10*int(100*del)) THEN
    WRITE (dumc, '(F4.2, X, I2, X, I1)') del, RayInfo%nAziAngle/2, RayInfo%nPolarAngle
  ELSE
    WRITE (dumc, '(F5.3, X, I2, X, I1)') del, RayInfo%nAziAngle/2, RayInfo%nPolarAngle
  END IF
  
  mesg = 'Ray    : ( ' // trim(dumc) // ' )'
  CALL message(io8, TRUE, TRUE, mesg)
  
  IF (nTracerCntl%TRSolver .NE. 1) THEN 
    mesg = 'SCHEME : multigroup MC' 
    CALL message(io8, TRUE, TRUE, mesg) 
  ELSE
    IF (nTracerCntl%lFeedback) THEN 
      mesg = 'T/H    : T' 
    ELSE
      mesg = 'T/H    : F' 
    END IF  
    CALL message(io8, TRUE, TRUE, mesg) 
    
    IF (nTracerCntl%lNodeMajor) THEN 
      IF (nTracerCntl%ScatOd .GT. 0) THEN 
        mesg = 'SCHEME : Node Major Pn' 
      ELSE 
        mesg = 'SCHEME : Node Major Trc' 
      END IF 
      CALL message(io8, TRUE, TRUE, mesg) 
    ELSE 
      IF (nTracerCntl%ScatOd .GT. 0) THEN 
        mesg = 'SCHEME : Group Major Pn' 
      ELSE 
        mesg = 'SCHEME : Group Major Trc' 
      END IF 
      CALL message(io8, TRUE, TRUE, mesg) 
    END IF 
     
    IF (nTracerCntl%lLinSrcCASMO) THEN 
      IF (nTracerCntl%lHybrid) THEN 
        mesg = 'LINSRC : T (Hybrid)' 
      ELSE 
        mesg = 'LINSRC : T' 
      END IF 
      CALL message(io8, TRUE, TRUE, mesg) 
    ELSE 
      mesg = 'LINSRC : F' 
      CALL message(io8, TRUE, TRUE, mesg) 
    END IF 
     
    IF (nTracerCntl%lDomainDcmp) THEN 
      mesg = 'DCMP   : T' 
      CALL message(io8, TRUE, TRUE, mesg) 
    ELSE 
      mesg = 'DCMP   : F' 
      CALL message(io8, TRUE, TRUE, mesg) 
    END IF 
     
    IF (nTracerCntl%lAFSS) THEN 
      mesg = 'AFSS   : T' 
      CALL message(io8, TRUE, TRUE, mesg) 
    ELSE 
      mesg = 'AFSS   : F' 
      CALL message(io8, TRUE, TRUE, mesg) 
    END IF 
  END IF 
END IF 
! ----------------------------------------------------
#ifndef MPI_ASYNC_READ
#ifdef MPI_ENV
CALL MPIWaitTurn(PE%MPI_COMM, PE%myrank, PE%nproc, TRUE) ! Waiting for other processor reading input
#endif
#else
CALL MPI_SYNC(PE%MPI_COMM)
#endif

IF (nTracerCntl%lISOTXS) THEN 
  CALL PostProc_ISOTXS
  IF (nTracerCntl%lGamma) CALL PostProc_GAMISO
END IF

IF (lMCP_Restart) CALL MCP_Initialize
! ----------------------------------------------------

END SUBROUTINE readinput