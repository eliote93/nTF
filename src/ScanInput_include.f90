#include <defines.h>
SUBROUTINE ScanInput_Include()
use param
use files,          only : filename,    InputFileIdx,  SteamTblIdx,                 &
                           TempFileIdx, io5,           io6,            io8,         &
                           io_quick,    CaseId,        TempInclIdx,    lIncFile      
use ioutil,         only : terminate,   toupper,       openfile,       IFnumeric,   &
                           nfields,     getfn
use inputcards ,    only : oneline,     probe,         mxcard,         nblock,      &
                           FindBlockId, FindCardId,    cards,          blocks
use geom,           only : LEDGE,        core,          nCellX0,       nAsyType0,   &
                           CellPitch,    AsyPitch,      nz,            ng,          &
                           nPinType0,    lgap,          lRot,          lCbd,        &
                           nCellType0,   nCellType,     lCoreAng,      nPinType,    &
                           nAsyType,       nCellX,      nPreC,         Albedo,      &
                           nAsyGapType,  lhgap
USE CNTLROD_mod,    ONLY : AddNCrConf
USE CrCsp_mod,      ONLY : SetnCrCspDat
USE PE_MOD,         ONLY : PE
USE Core_Mod,       ONLY : GroupInfo                           
USE Material_Mod,   ONLY : Mixture,      nMixType
USE Nuclidmap_mod,  ONLY : nuclidmap
USE TRAN_MOD,       ONLY : TranCntl
USE cntl,           ONLY : nTracerCntl
USE itrcntl_mod,    ONLY : ConvItrCntl
USE MCP_Util,       ONLY : nMCP,         MCP_Path,     lMCP_restart

#ifdef MPI_ENV
USE MPIComm_Mod,    ONLY : MPIWaitTurn, MPI_SYNC
#endif        
IMPLICIT NONE

CHARACTER(15)   :: cardname, blockname, astring, optfield
INTEGER         :: nCardInp(mxcard, nblock)
INTEGER         :: indev
INTEGER         :: idblock,idcard
INTEGER         :: iline                    !Number of Line of Input file
INTEGER         :: nMxCardInp               !Maximum number of card in the one block
INTEGER         :: nLineField
INTEGER         :: iasy, icell, ipin              !Assembly, Cell, Pin Number Idx
INTEGER         :: i, j, k,  NB
REAL            :: REALTemp(100)
INTEGER         :: MCP_temp
CHARACTER(512)  :: MCP_oneline
CHARACTER(512)  :: Include_Path
CHARACTER(7)    :: Include_char
INTEGER :: index, lenblank

!Input File Open
indev = io5
CALL OpenFile(io5, TRUE, FALSE, FALSE, filename(InputFIleIdx))
!Preminlary Card and Block Serch
nCardInp = 0
i = 0
DO WHILE (TRUE)
  read(indev,'(a512)',END=5) oneline
  i = i + 1
  IF(probe.eq.BANG) cycle;     IF(probe.eq.POUND) cycle
  IF(oneline.eq.BLANK0) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe.eq.DOT) exit;       IF(probe.eq.SLASH) exit
  read(oneline,*) cardname
  CALL toupper(cardname)
  IF(probe .ne. BLANK) THEN   !BLOCK Name
    blockname = cardname
    idblock = FindBlockId(blockname)
    IF( idblock .eq. 0) then
      mesg='BLOCK '//trim(blockname)//' not available'
      CALL terminate(mesg)
    ENDIF
  ELSE                        !Case Name which is belong to the certain group
    idcard = FindCardId(idblock,cardname)
    IF( idcard .eq. 0) then
      mesg='CARD '//trim(cardname)//' not available'
      CALL terminate(mesg)
    ENDIF
    nCardInp(idcard,idblock) = nCardInp(idcard,idblock) + 1
  ENDIF
ENDDO
!Unexpected 
5 mesg='CASE Input Should be ended with DOt(.) or slash(/)'
IF(probe.ne.DOT .AND. probe.ne.SLASH) CALL terminate(mesg)
rewind(indev)

!Determine Maximum  number of cards in the one block
nMxCardInp = MAXVAL(nCardInp)

Filename(TempInclIdx) = trim(CaseId)//trim('.inptmp')
!IF(PE%SLAVE) RETURN
!CALL OpenFIle(86, FALSE, FALSE, FALSE, filename(TempInclIdx))
REWIND(INDEV)


DO
  read(indev,'(a512)', END=6) oneline
  IF(oneline.eq.BLANK) cycle;
  lenblank = 0  
  DO index = 1, 30
    IF(oneline(index:index) .EQ. ' ') THEN
      lenblank = lenblank + 1
    ELSE
      EXIT
    END IF
  END DO
  IF (lenblank .LT. 30) THEN
    Include_char(1:7) = oneline(lenblank+1:lenblank+8)
  END IF
  CALL toupper(Include_Char)
  IF(Include_Char .NE. 'INCLUDE') WRITE(86, '(A)') TRIM(oneline)
 
  IF(probe.eq.BANG) cycle;     IF(probe.eq.POUND) cycle
  IF(IFnumeric(oneline)) cycle
  IF(probe.eq.DOT) exit;

  read(oneline,*) cardname; CALL toupper(cardname)
  nLineField = nfields(oneline)-1
  SELECT CASE(cardname)
    CASE('CASEID') 
      READ(oneline,*) ASTRING, CaseId  
      Filename(TempInclIdx) = trim(CaseID)//trim('.inptmp')
      IF(PE%master) THEN
        CALL OpenFIle(86, FALSE, FALSE, FALSE, filename(TempInclIdx))
      END IF
    CASE('INCLUDE')
      
      IF(PE%master) THEN
        CALL GETFN(oneline, 2, Include_Path)
        OPEN(UNIT = 87, File = Include_Path, STATUS='OLD', ACTION = 'READ', POSITION = 'REWIND')
        DO
          READ(87, '(a256)', END=1000) Oneline
          WRITE(87, '(A)') TRIM(Oneline)
        END DO
        1000 CONTINUE
        CLOSE(87)
      END IF
  ENDSELECT
ENDDO
6 CONTINUE

filename(InputFIleIdx) = filename(TempInclIdx)

CLOSE(io5)
IF(PE%master) CLOSE(86)

END SUBROUTINE