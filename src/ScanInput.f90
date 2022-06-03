#include <defines.h>
SUBROUTINE ScanInput()
use param
use files,          only : filename,    InputFileIdx,  XsFileIdx,      SteamTblIdx, &
                           TempFileIdx, io5,           io6,            io8,         &
                           io_quick,    CaseId,        TempInclIdx,    lIncFile,    &
                           rlbIdx,      slbIdx,        DeplFileIdx,    phlIdx, PXSIDX, pmatrxIdx
use ioutil,         only : terminate,   toupper,       openfile,       IFnumeric,   &
                           nfields,     getfn,         icolfield,      expdast, nhexdat
use inputcards ,    only : oneline,     probe,         mxcard,         nblock,      &
                           FindBlockId, FindCardId,    cards,          blocks
use geom,           only : LEDGE,        core,          nCellX0,       nAsyType0,   &
                           CellPitch,    AsyPitch,      nz,            ng,          &
                           nPinType0,    lgap,          lRot,          lCbd,        &
                           nCellType0,   nCellType,     lCoreAng,      nPinType,    &
                           nAsyType,       nCellX,      nPreC,         Albedo,      &
                           nAsyGapType,   lhgap, ngthrm,                            & 
                           nBasecell,   nBasecell0,                                 & ! --- 180625 JSR  
                           !--- CNJ Edit : Flexible Cell & Pin Numbering
                           CellTypeMap,  PinTypeMap,                                &
                           nGapType,     nGapPinType,   nVssTyp,                    & ! KSC Edited - 170807 
                           nCrCell
USE CNTLROD_mod,    ONLY : AddNCrConf
USE CrCsp_mod,      ONLY : SetnCrCspDat
USE PE_MOD,         ONLY : PE
USE Core_Mod,       ONLY : GroupInfo
USE Material_Mod,   ONLY : Mixture,      nMixType
USE ReadXsec
USE BenchXs,        ONLY : scatod
USE TRAN_MOD,       ONLY : TranCntl,     TranInfo
USE cntl,           ONLY : nTracerCntl
USE itrcntl_mod,    ONLY : ConvItrCntl
USE MCP_Util,       ONLY : nMCP,         MCP_Path,     lMCP_restart,    nMCP_Plane
USE SubChCoupling_mod,    ONLY: is_coupled
USE DEPL_MOD,      ONLY : DeplCntl
USE DeplLib_MOD,   ONLY : DeplLib
USE PointXSRT_MOD,       ONLY : ReadPWXS
USE HexData,       ONLY : hLgc
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
INTEGER         :: i, j, k,  NB, ist
INTEGER         :: maxCellID, maxPinID      !--- CNJ Edit : Flexible Cell & Pin Numbering
INTEGER         :: matfid = 28, ibasecell        !--- 180531 JSR Edit
REAL            :: REALTemp(100)
INTEGER         :: MCP_temp
CHARACTER(512)  :: MCP_oneline
CHARACTER*512   :: aline, bline
INTEGER :: xsod
!INTEGER,pointer :: iCard(:,:)

!Input File Open
#ifndef MPI_ASYNC_READ
#ifdef MPI_ENV
CALL MPIWaitTurn(PE%MPI_COMM, PE%myrank, PE%nproc, .FALSE.)
#endif
#endif

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

!--- JSR Edit : nTIG Restart
DO WHILE (TRUE)
  read(indev,'(a512)',END=5) oneline
  i = i + 1
  IF(probe.eq.BANG) cycle;     IF(probe.eq.POUND) cycle
  IF(oneline.eq.BLANK0) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe.eq.DOT) exit;       IF(probe.eq.SLASH) exit
  read(oneline,*) cardname
  CALL toupper(cardname)
  IF (cardname .EQ. 'NTIG') THEN
    READ(oneline, *) astring, nTracerCntl%lnTIGRst; EXIT
  ENDIF
ENDDO

REWIND(indev)

IF (nTracerCntl%lnTIGRst) THEN
  OPEN(matfid, file = "MATERIAL.bin", form = 'unformatted')
  READ(matfid) nMixType
  CLOSE(matfid)
ENDIF

!Determine Maximum  number of cards in the one block
nMxCardInp = MAXVAL(nCardInp)

nAsyType0 = 0; nCellType0 = 0; nPinType0 = 0; nAsyGapType = 0
nGapType = 0; nGappinType = 0; nVssTyp = 0; nPinType = 0; nCellType = 0 ! KSC EDIT 17/12/11
nBasecell0 = 0; nBasecell = 0 ! --- 180625 JSR
maxCellID = 0; maxPinID = 0   !--- CNJ Edit : Flexible Cell & Pin Numbering
nCrCell = 0

DO while(TRUE)
  read(indev,'(a512)') oneline
  IF(probe.eq.BANG) cycle;     IF(probe.eq.POUND) cycle
  IF(oneline.eq.BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe.eq.DOT) exit;       IF(probe.eq.SLASH) cycle
  read(oneline,*) cardname; CALL toupper(cardname)
  nLineField = nfields(oneline)-1
  SELECT CASE(cardname)
    CASE('CASEID')
      READ(oneline,*) ASTRING, CaseId
    CASE('RESTART')

    CASE('MIXTURE')
      READ(oneline,*) ASTRING, i
      nMixType = Max(i, nMixType)
    CASE('TRISO')

    CASE('HEX')
      READ(oneline,*) ASTRING, nTracerCntl%lHex

    CASE('NPINS')
      READ(oneline,*) ASTRING, nCellX0
      IF(nLineField .GT. 1) READ(oneline, *) ASTRING, nCellX0, Core%nAsyGT
      Core%nASyCell = nCellX0*nCellX0
    CASE('PITCH')
      READ(oneline,*) ASTRING, CellPitch
      IF(nLineField .gt. 1) READ(oneline,*) ASTRING, CellPitch, AsyPitch
    CASE('AX_MESH')
      ist = icolfield(oneline, 2)
      aline = oneline(ist:256)
      CALL expdast(aline, bline)
      nz = nhexdat(bline)
      CORE%nz = nz; PE%nz = nz
    CASE('CELL')
      read(oneline,*) astring, icell
      !--- CNJ Edit : Flexible Cell & Pin Numbering
      nCellType0 = max(icell, nCellType0)
!      maxCellID = max(icell, maxCellID)
    CASE('BASECELL')
      nBasecell0 = nBasecell0 + 1
      ! --- 180626 JSR
      READ(oneline,*) ASTRING, ibasecell
    CASE('CAD_CELL')
      read(oneline,*) astring, icell
      !--- CNJ Edit : Flexible Cell & Pin Numbering
      nCellType0 = max(icell, nCellType0)
!      maxCellID = max(icell, maxCellID)
    CASE('GAP_CELL')
      read(oneline,*) astring, icell
      nCellType0 = max(icell, nCellType0)
    CASE('GAP')
      read(oneline,*) astring, lgap
    CASE('HGAP')
      read(oneline,*) astring, lhgap
    CASE('PIN')
      read(oneline,*) astring, ipin
      !--- CNJ Edit : Flexible Cell & Pin Numbering
      nPinType0 = max(ipin, nPinType0)
!      maxPinID = max(ipin, maxPinID)
    CASE('GAP_PIN')
      read(oneline,*) astring, ipin
      !--- CNJ Edit : Flexible Cell & Pin Numbering
      nPinType0 = max(ipin, nPinType0)
!      maxPinID = max(ipin, maxPinID)
    CASE('ASSEMBLY')
      read(oneline,*) astring, iasy
      nAsyType0 = max(iasy,nAsyType0)
    CASE('GAP_ASY')
      read(oneline, *) astring, iasy
      nAsyGapType = max(nAsyGapType, iasy)
    CASE('ALBEDO')
      read(oneline,*) astring,albeDO(1),(REALTemp(k),k=2,nLineField)
      nb = 4
      DO i= 2, nb
        albeDO(i) = albeDO(1)
        IF(nLineField .gt. 3) albeDO(i) = REALTemp(i)
      ENDDO
      albeDO(nb+1) = REALTemp(nLineField)
      albeDO(nb+2) = REALTemp(nLineField)
      IF((nLineField .eq. 3) .or. (nLineField .eq. (nb+2))) albeDO(nb+1) = REALTemp(nLineField-1)
      IF(lCoreAng .ne. 360) then ! 90deg symmetry
        albeDO(nb-1) = 0
        albeDO(nb-2) = 0
      ENDIF
    CASE('RAD_CONF')
      core%nxa=0
      core%nya=0
      core%nxya=0
      read(oneline,*) cardname, lCoreAng
      IF(nLineField .ne. 1) then
        read(oneline,*) cardname, lCoreAng, astring
        CALL TOUPPER(astring)
        IF(astring .eq. 'EDGE') then
          lEdge = TRUE
        ELSEIF(astring .eq. 'CENT') then
          lEdge = FALSE
        ELSEIF(astring .eq. 'CHECKER') then
          lCbd = .TRUE.
        ENDIF
      ENDIF
      IF(nLineField .GE. 3) THEN
        read(oneline,*) cardname, lCoreAng, astring, optfield
        CALL TOUPPER(optfield)
        IF(optfield .EQ. 'ROT' .AND. lCoreAng .EQ. 90) lRot = .TRUE.
      ENDIF

      DO WHILE(TRUE)
        read(indev,'(a256)') oneline
        IF(.not. IFnumeric(oneline)) exit
        nLineField = nfields(oneline)
        core%nxa = max(core%nxa, nLineField)
        core%nya = core%nya + 1
        core%nxya = core%nxya + nLineField
      ENDDO
      BACKSPACE(indev)
      IF(lCbd .AND. lRot) CALL Terminate('Error : Checkerboard and Rotational Configuration are not compitable')
      IF(lCbd .AND.  (core%nxya .NE. core%nxa * core%nxa)) CALL Terminate('Dummy Assembly in not allowed in Checkerboard Configuraiton')
    CASE('VESSEL')
      read(oneline,*) astring, icell
      nVssTyp = max(iCell, nVssTyp)
      hLgc%lvss = TRUE
    CASE('VYGORODKA')
      hLgc%lvyg = TRUE
    CASE('CR_CELL')
      CALL AddNCrConf(1)
      nCellType0 = nCellType0 + 1
      nCrCell = nCrCell + 1
    CASE('CR_ASYCONF')
      CALL AddNCrConf(2)
    CASE('CR_BANK')
      CALL AddNCrConf(3)
    CASE('CSP_FILE')
      READ(oneline,*) ASTRING, i
      CALL SetnCrCspDat(i)
    CASE('LIB_TYPE')
      ! nTRACERCntl%libtyp   1:benchmark 2:mlb 3:plc 4:isotxs
      IF(nLineField .eq. 1) then
        read(oneline,*) astring, nTracerCntl%libtyp
      ELSEIF(nLineField .eq. 2) then
        read(oneline,*) astring, nTracerCntl%libtyp, nTracerCntl%XsOrder
      ELSEIF(nLineField .eq. 3) then
        read(oneline,*) astring, nTracerCntl%libtyp, nTracerCntl%XsOrder, nTracerCntl%lTrCorrection
      ENDIF

      IF(nTracerCntl%libtyp .eq. 1) nTracerCntl%lBenchXs = .TRUE.   ! nTRACER benchmark format, FULL Scattering matrix
      IF(nTracerCntl%libtyp .eq. 11) nTracerCntl%lBenchXs = .TRUE.   ! nTRACER benchmark format, NEACRP-l-335
      IF(nTracerCntl%libtyp .eq. 4) nTracerCntl%lISOTXS = .TRUE.    ! ISOTXS Format from EXUS-F and MCC-3
      nTracerCntl%lXsLib = .not. nTracerCntl%lBenchXs

      !IF(nLineField .EQ. 2) THEN
      !  read(oneline,*) astring, i, j
      !  nTracerCntl%ScatOd = j
      !  IF( j .GT. 0) nTracerCntl%lscat1 = .TRUE.
      !ENDIF
    !--- CNJ Edit : Read Gamma Trasnport Order
    CASE('SCAT_ORDER')
        IF(nLineField .EQ. 1) THEN
            read(oneline,*) astring, nTracerCntl%ScatOd
            IF(nTracerCntl%ScatOd .GT. 0)THEN
                nTracerCntl%lscat1 = .TRUE.
                nTracerCntl%lTrcorrection =.FALSE.
            ENDIF
            scatod=nTracerCntl%ScatOd
        ELSEIF(nLineField .EQ. 2) THEN
            read(oneline,*) astring, nTracerCntl%ScatOd, nTracerCntl%GammaScatOd
            IF(nTracerCntl%ScatOd .GT. 0)THEN
                nTracerCntl%lscat1 = .TRUE.
                nTracerCntl%lTrcorrection =.FALSE.
            ENDIF
            IF(nTracerCntl%GammaScatOd .GT. 0)THEN
                nTracerCntl%lGammaScat1 = .TRUE.
            ENDIF
            scatod=nTracerCntl%ScatOd
        ENDIF
    CASE('BCR_OPT')
        read(oneline,*) astring, nTracerCntl%lBCRopt
    CASE('FILE')
      CALL getfn(oneline,2,filename(XsFileIdx))
#ifdef __GAMMA_TRANSPORT
    CASE('PHL')
      CALL getfn(oneline,2,filename(phlIdx))
      nTracerCntl%lGamma = .TRUE.
    CASE('PMATRX')
      CALL getfn(oneline,2,filename(pmatrxIdx)) 
#endif
    CASE('RLB')
      CALL getfn(oneline,2,filename(rlbIdx))
      nTracerCntl%lrestrmt=.TRUE.
    CASE('SLB')
      CALL getfn(oneline,2,filename(slbIdx))      
      nTracerCntl%lSSPH=.TRUE.
    CASE('DEPLFILE')
      CALL getfn(oneline,2,filename(DeplFileIdx))
      DeplLib%Filename = filename(DeplFileIdx)
      DeplCntl%lDeplFile = .TRUE.
    CASE('GROUP_SPEC')
      read(oneline,*) astring, ng !, ngthrm
    CASE('GC_SPEC')
      read(oneline, *) astring, nTracerCntl%nGC
      nTracerCntl%lGcCmfdGrp = .TRUE.
    CASE('GC_CMFD')
      read(oneline, *) astring, nTracerCntl%lGcCmfd
    CASE('LINSRC')
      READ(oneline, *) astring, nTracerCntl%lLinSrc
    CASE('RSTFILE')

    CASE('STEAM_TBL')
      CALL GetFn(oneline, 2, filename(SteamTblIdx))
    CASE('BORON')
      READ(oneline, *) astring, nTracerCntl%BoronPPM
      IF(nTracerCntl%BoronPPM .GT. 0._8) nTracerCntl%lInitBoron = .TRUE.
    CASE('XENON')
      READ(oneline, *) astring, optfield
      nTracerCntl%lXeDyn = .TRUE.
    CASE('XS_CHANGE')
      TranCntl%nChange = TranCntl%nChange + 1    
    CASE('XS_NOISE')
      TranCntl%nNoise = TranCntl%nNoise + 1
    CASE('XS_CNTLROD')
      nTracerCntl%nCntlRod = nTracerCntl%nCntlRod + 1
    CASE('KIN_BENCH')
      READ(oneline, *) astring, TranCntl%lKineticBen
      nTracerCntl%lKineticBen = TranCntl%lKineticBen
    CASE('DYN_BENCH')
      IF(nLineField .Eq. 1) THEN 
        READ(oneline, *) astring, TranCntl%lDynamicBen, TranInfo%InitTemp
      ELSE
        READ(oneline, *) astring, TranCntl%lDynamicBen, TranInfo%InitTemp, TranInfo%Initpow
      END IF
      nTracerCntl%lDynamicBen = TranCntl%lDynamicBen
    CASE('AF_SRC')
      READ(oneline, *) astring, nTracerCntl%lAfSrc, TranCntl%lMethod
    CASE('COBRA_TF')
      READ(oneline, *) astring, is_coupled
    CASE('ESCOT')
      READ(oneline, *) astring, is_coupled
    CASE('NMAXOUTER')
      READ(oneline, *) astring, ConvItrCntl%OuterMax
    CASE('EIG_CONV')
      READ(oneline, *) astring, ConvItrCntl%eigconv
    CASE('PSI_CONV')
      READ(oneline, *) astring, ConvItrCntl%psiconv
    CASE('RES_CONV')
      READ(oneline, *) astring, ConvItrCntl%resconv
    CASE('DECUSP_CONV')
      READ(oneline, *) astring, ConvItrCntl%decuspconv
    CASE('NMCP')
      READ(oneline, *) astring, nMCP
      ALLOCATE(MCP_Path(1:nMCP))
    CASE('NMCP_PLANE')
      READ(oneline, *) astring, nMCP_Plane
    CASE('MCP_CALL')
      READ(oneline, *) astring, lMCP_restart
    CASE('ED')
      nTracerCntl%lED = .true.
    CASE('MCXSFILE')
      READ(oneline, *) astring, nTracerCntl%MCXSfile  
      nTracerCntl%lMCXS = .true.  
      call readMCXSfile(io_quick,nTracerCntl%MCXSfile)
    CASE('BURNUP')
      DeplCNTL%nBurnUpStep = nLineField
      READ(ONELINE, *) ASTRING, REALTemp(1:nLineField)

      DO WHILE(TRUE)
        read(indev,'(a256)') oneline
        IF(.NOT. IFnumeric(oneline)) EXIT
        nLineField = nfields(oneline)
        DeplCNTL%nBurnUpStep = DeplCNTL%nBurnUpStep + nLineField
      ENDDO
      Backspace(indev);
    CASE('PSM') ! Point-wise Slowing-down Method
      nTRACERCntl%lPSM = .TRUE.
      nTracerCntl%lrestrmt=.TRUE.
      nTRACERCntl%lpointrt = .TRUE.
      CALL getfn(oneline,2,filename(pxsIdx))
    CASE('OTFRIF')
      nTRACERCntl%lOTFRIF  = .TRUE.
      nTRACERCntl%lRIF     = .TRUE.
      nTRACERCntl%lrestrmt = .TRUE.
      nTRACERCntl%lpointrt = .TRUE.
      CALL getfn(oneline,2,filename(pxsIdx))
  ENDSELECT
ENDDO

IF (nTracerCntl%lHex) THEN
  REWIND(indev)
  DO while(TRUE)
    read(indev,'(a512)') oneline
    IF(probe.eq.BANG) cycle;     IF(probe.eq.POUND) cycle
    IF(oneline.eq.BLANK) cycle;  IF(IFnumeric(oneline)) cycle
    IF(probe.eq.DOT) exit;       IF(probe.eq.SLASH) exit
    read(oneline,*) cardname; CALL toupper(cardname)
    nLineField = nfields(oneline)-1
    SELECT CASE(cardname)
      CASE('CELL')
        read(oneline,*) astring, icell
        nCellType  = max(icell, nCellType) ! KSC EDIT 17/12/11
      CASE('GAP_CELL')
        read(oneline,*) astring, icell
        nGapType   = max(iCell, nGapType) ! KSC EDIT 17/12/11
      CASE('PIN')
        read(oneline,*) astring, ipin
        nPinType  = max(ipin, nPinType) ! KSC EDIT 17/12/11
      CASE('GAP_PIN')
        read(oneline,*) astring, ipin
        nGapPinType = max(ipin, nGapPinType) ! KSC EDIT 17/12/11
    ENDSELECT
  ENDDO
ENDIF

! REWIND(indev)
! 
! !--- CNJ Edit : Flexible Cell & Pin Numbering
! nAsyType0 = 0; nCellType0 = 0; nPinType0 = 0
! ALLOCATE(CellTypeMap(maxCellID), PinTypeMap(maxPinID))
! DO while(TRUE)
!   read(indev,'(a512)') oneline
!   IF(probe.eq.BANG) cycle;     IF(probe.eq.POUND) cycle
!   IF(oneline.eq.BLANK) cycle;  IF(IFnumeric(oneline)) cycle
!   IF(probe.eq.DOT) exit;       IF(probe.eq.SLASH) exit
!   read(oneline,*) cardname; CALL toupper(cardname)
!   nLineField = nfields(oneline)-1
!   SELECT CASE(cardname)
!     CASE('CELL')
!       read(oneline,*) astring, icell
!       nCellType0 = nCellType0 + 1
!       CellTypeMap(icell) = nCellType0
!     CASE('GAP_CELL')
!       read(oneline,*) astring, icell
!       nCellType0 = nCellType0 + 1
!       CellTypeMap(icell) = nCellType0
!     CASE('PIN')
!       read(oneline,*) astring, ipin
!       nPinType0 = nPinType0 + 1
!       PinTypeMap(ipin) = nPinType0
!     CASE('GAP_PIN')
!       read(oneline,*) astring, ipin
!       nPinType0 = nPinType0 + 1
!       PinTypeMap(ipin) = nPinType0
!   ENDSELECT
! ENDDO

IF (.NOT. nTracerCntl%lHex) THEN ! KSC EDIT - 17/12/11
  !Cell and Pincell Number
  nCellX = nCellX0
  
  !nCellType0 = nCellType0 + 1 !Cell For Baffle
  !nTracerCntl%StructCell = nCellType0
  IF( lgap .AND. lhgap )THEN
      CALL Terminate('Error : Both two GAP and HGAP option cannot be used at the same time')
  ENDIF
  
  IF(lgap) then
    nBasecell0 = nBasecell0 + 3 ! --- 180724 JSR
    nCellType0 = nCellType0 + 3
    nPinType0 = nPinType0 + 3
    nCellX = nCellX0 + 2
    nTracerCntl%lGap = .TRUE.
    IF(nAsyGapType .GT. 0) nTracerCntl%lUsrGap = .TRUE.
  !--- BYS edit start 16/05/23/
  ELSEIF(lhgap)THEN
    nCellType0 = nCellType0 *9
    nPinType0 = nPinType0 *9
    nCellX = nCellX0 + 2
    nTracerCntl%lhGap = .TRUE.
  !--- BYS edit end 16/05/23/
  ENDIF
  
  nBasecell = nBasecell0
  nCellType = nCellType0
  nPinType = nPinType0
  nAsyType = nAsyType0
  
  IF(nz .EQ. 1 .and. albedo(nb+1) .EQ. ZERO .and. albedo(nb+2) .EQ. ZERO) THEN
    nTracerCntl%l3dim = FALSE
    nTracerCntl%nDim = 2
  ENDIF
  
  IF(.not. lEdge) then
    nBasecell = nBasecell*4 ! --- 180724 JSR
    nCellType = nCellType0*4
    nPinType = nPinType0*4
    nAsyType = nAsyType*4
  ENDIF
ENDIF

!Read XS file
IF(nTracerCntl%lxslib) CALL ReadDeplFile(io_quick, FileName(DeplFileIdx), DeplLib(1), GroupInfo%ntiso_depl)
IF(nTracerCntl%libtyp .eq.0) THEN !Binary File, HELIOS format
  !CALL Terminate('Error : ''Lib_Type'' card should not be 0 anymore.')
  CALL ReadXSL(io_quick, filename(XsFileIdx), GroupInfo%ng, GroupInfo%nofg, GroupInfo%norg, GroupInfo%nTiso)
  CALL setXeDynVar(DeplLib(1))
  CALL setKappa(DeplLib(1))
  ng = GroupInfo%ng
ELSEIF (nTracerCntl%libtyp .eq. 1) THEN !Benchmark XS
  CALL openfile(io_quick, TRUE, FALSE, FALSE, filename(XsFileIdx))
  xsod=nTracerCntl%XsOrder
  IF( xsod .EQ. -1 )THEN  ! Xs order not given
    xsod=nTracerCntl%ScatOd
  ENDIF
  CALL ReadBenchXs(io_quick, io8, ng, nPrec, TranCntl%nChange, XsOd)  !--- r544 edit > r563 updated
  close(io_quick)
ELSEIF (nTracerCntl%libtyp .eq. 2) THEN !Default, by PHS
  CALL readMLX(io_quick, filename(XsFileIdx), GroupInfo%ng, GroupInfo%nofg,   &
                  GroupInfo%norg, GroupInfo%nTiso, GroupInfo%ngg)
  CALL setXeDynVar(DeplLib(1))
  CALL setKappa(DeplLib(1))  
  IF (nTracerCntl%lrestrmt) THEN
      CALL readRLB(io_quick, filename(rlbIdx))
  ELSE
      nTracerCntl%lRIF=.FALSE.; nTracerCntl%lRST=.FALSE.
  ENDIF
  IF (nTracerCntl%lsSPH) CALL readSLB(io_quick, filename(slbIdx), nTracerCntl%ScatOd)
  ng = GroupInfo%ng
#ifdef __GAMMA_TRANSPORT
  IF (nTracerCntl%lGamma) CALL ReadPHL(io_quick,filename(phlIdx),GroupInfo%NELE,GroupInfo%nTiso,GroupInfo%ngg)
#endif
  IF (nTRACERCntl%lPSM.OR.nTRACERCntl%lOTFRIF) THEN
    CALL ReadPWXS(io_quick, filename(pxsIdx))
  END IF
ELSEIF (nTracerCntl%libtyp .eq. 3) THEN !PLC, by LCH
  nTracerCntl%lRIF=.FALSE.; nTracerCntl%lrestrmt = .TRUE.
  CALL ReadPLC(io_quick, filename(XsFileIdx), GroupInfo%ng, GroupInfo%nofg, GroupInfo%norg, GroupInfo%nTiso)
  CALL setXeDynVar(DeplLib(1))
  CALL setKappa(DeplLib(1))
  ng = GroupInfo%ng
ELSEIF (nTracerCntl%libtyp .eq. 11) THEN !NEACRP-l-335 benchmark problem
  CALL openfile(io_quick, TRUE, FALSE, FALSE, filename(XsFileIdx))
  xsod=nTracerCntl%XsOrder
  IF( xsod .EQ. -1 )THEN  ! Xs order not given
    xsod=nTracerCntl%ScatOd
  ENDIF
  CALL ReadBenchXs_NEACRP(io_quick, ng, nPrec, TranCntl%nChange, XsOd)  !--- r544 edit > r563 updated
  close(io_quick)
ELSEIF (nTRACERCntl%libtyp .EQ. 4) THEN ! ISOTXS, by JSU
  nTracerCntl%lRIF=.FALSE.; nTracerCntl%lrestrmt = .FALSE.
  CALL readISOTXS(io_quick, filename(XsFileIdx), GroupInfo%ng, GroupInfo%nTiso, XsOd)
  IF(nTracerCntl%ScatOd > XsOd) CALL terminate('SCAT_ORDER is Larger than the order in XS Lib. !')
  nTracerCntl%lGcCMFD = .FALSE.
#ifdef __GAMMA_TRANSPORT
  IF (nTracerCntl%lGamma) THEN 
    CALL ReadGAMISO(io_quick,filename(phlIdx),GroupInfo%ngg,GroupInfo%nele)
    CALL ReadPMATRX(io_quick,filename(pmatrxIdx))
  END IF
#endif
ELSE
  CALL Terminate('Error : ''Lib_Type'' card should be given between 0 - 4.')
ENDIF

#ifndef MPI_ASYNC_READ
#ifdef MPI_ENV
CALL MPIWaitTurn(PE%MPI_COMM, PE%myrank, PE%nproc, .TRUE.)
#endif

!COpy Input File

CALL CopyInputFile(io5, io8, CaseID)   !--- CNJ Edit : io0 --> io8
CLOSE(io5)
#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_COMM)
#endif
CALL OpenFIle(io5, TRUE, FALSE, FALSE, filename(TempFileIdx))
#else
REWIND(io5)
#endif

#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_COMM)
#endif
!
END SUBROUTINE

SUBROUTINE CopyInputFile(InDev, OutDev, CaseID)
USE PARAM
USE ioutil,         only : terminate,   toupper,       openfile
USE files,          only : filename,    TempFIleIdx
USE inputcards ,    only : oneline,     probe
USE PE_MOD,         ONLY : PE
IMPLICIT NONE
INTEGER :: InDev, OutDev
CHARACTER(80) :: CaseID

Filename(TempFIleIdx) = trim(CaseId)//trim('.tmp')
IF(PE%SLAVE) RETURN
CALL OpenFIle(OutDev, FALSE, FALSE, FALSE, filename(TempFIleIdx))
REWIND(INDEV)
DO
  READ(indev,'(a512)',END=6) oneline
  WRITE(OutDev,'(a512)') oneline
  IF(probe.eq.DOT) exit
ENDDO
6 CONTINUE
CLOSE(OutDev)
REWIND(INDEV)
END SUBROUTINE

SUBROUTINE readMCXSfile(indev,file)
USE ioutil,   ONLY : openfile  
use CNTL,     ONLY : nTracerCntl
IMPLICIT NONE
integer,intent(in) :: indev
character*512,intent(in) :: file

integer :: nc,ic,ig,nfxr,igi,igf,igx,iasyx,iasyy,ix,iy,ifxr1,ifxr2,iz,niso,iso

igi=10; igf=25

call openfile(indev,.TRUE.,.FALSE.,.FALSE.,file)
read(indev,*) nc
nTracerCntl%nmcxs=nc
allocate(nTracerCntl%mcxs(nc))
do ic=1,nc  
    read(indev,*) iz,iasyx,iasyy,ix,iy,ifxr1,ifxr2,niso
    nTracerCntl%mcxs(ic)%iz=iz
    nTracerCntl%mcxs(ic)%iasyx=iasyx
    nTracerCntl%mcxs(ic)%iasyy=iasyy
    nTracerCntl%mcxs(ic)%ix=ix
    nTracerCntl%mcxs(ic)%iy=iy
    nTracerCntl%mcxs(ic)%ifxr1=ifxr1
    nTracerCntl%mcxs(ic)%ifxr2=ifxr2
    nTracerCntl%mcxs(ic)%niso=niso 
    allocate(nTracerCntl%mcxs(ic)%idiso(niso))     
    nfxr=ifxr2-ifxr1+1
    allocate(nTracerCntl%mcxs(ic)%isoxsa(igi:igf,ifxr1:ifxr2,niso),nTracerCntl%mcxs(ic)%isoxsnf(igi:igf,ifxr1:ifxr2,niso))
    do iso=1,niso
        read(indev,*) nTracerCntl%mcxs(ic)%idiso(iso)
        do ig=igi,igf
            read(indev,*) igx,nTracerCntl%mcxs(ic)%isoxsa(ig,ifxr1:ifxr2,iso)
            read(indev,*) igx,nTracerCntl%mcxs(ic)%isoxsnf(ig,ifxr1:ifxr2,iso)
        enddo
    enddo
enddo
close(indev)
END SUBROUTINE