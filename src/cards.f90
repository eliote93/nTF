module inputcards
implicit none
!character(1) DOT,BANG,BLANK,SLASH,AST
!parameter (DOT='.',BANG='!',BLANK=' ',SLASH='/',AST='*')
INTEGER, parameter :: mxncol=512, mxnfield=128
INTEGER, parameter :: mxcard=30, nblock=25
character*1 probe,sc(mxncol),sc80(80)
character*512  oneline
!character*10 cardname,blockname


logical IFfile
equivalence (probe,oneline)
equivalence (sc(1),oneline)

INTEGER ::  ncard(nblock)
data ncard/ 1,  2, 20, 29, 3,            & !CASEID    STATE      XSEC      OPTION     MATERIAL
           25, 20,  6, 12, 8,            & !GEOM      TH         TRAN      DEPL       BRANCH
           24,  7,  5, 2,  9,            & !EDIT      DECOUPLE   PARALLEL  VISUAL     LP_SHF  !--- BYS edit : # of EDIT card
            9,  4,  3, 4,  6,            & !CNTLROD   CUSPING    XEDYN     CONV       MCP_RESTART
            4, 13, 10, 14, 3/              !SUBCH_OP  MKL        CUDA      GCGEN      NTIG_RESTART
character(15):: blocks(nblock),cards(mxcard,nblock)

data blocks /'CASEID'    ,'STATE'     ,'XSEC'      ,'OPTION'    ,'MATERIAL'        &
            ,'GEOM'      ,'TH'        ,'TRAN'      ,'DEPL'      ,'BRANCH'          &
            ,'EDIT'      ,'DECOUPLE'  ,'PARALLEL'  ,'VISUAL'    ,'LP_SHF'          &
            ,'CNTLROD'   ,'CUSPING'   ,'XEDYN'     ,'CONV'      ,'MCP_RESTART'     &
            ,'SUBCH_OP'  ,'MKL'       ,'CUDA'      ,'GCGEN'     ,'NTIG_RESTART'    /

data cards/   'RESTART'   ,29*' '                                                  &
!  2  STATE   2                                                                    &
      ,'CORE_POWER' ,'TH_COND'   ,28*' '                                           &
!  3  XSEC   20                                                                    &
      ,'LIB_TYPE'   ,'GROUP_SPEC','FILE'      ,'BASE_MICRO','DNEUT_CHI'            &
      ,'NEUT_VELO'  ,'DNP_NGRP'  ,'DNP_BETA'  ,'DNP_LAMBDA','DUMMY'                &
      ,'BUCKLING'   ,'RLB'       ,'SLB'       ,'GC_SPEC'   ,'DEPLFILE'             &
      ,'GFILE'      ,'RT'        ,'ED'        ,'MCXSFILE'  ,'RESO_OPT', 10*''      &
!  4  OPTION  29                                                                   &
      ,'FEEDBACK'   ,'CMFD'      ,'RAY'       ,'MOC'       ,'ITER_LIM'             &
      ,'GRP_ALBEDO' ,'MPI_GROUPS','GC_CMFD'   ,'SEARCH'    ,'LKGSPLIT'             &
      ,'BORON'      ,'XENON'     ,'NODAL'     ,'CRITICAL'  ,'LINSRC'               &
      ,'FRR'        ,'SCAT_ORDER','BCR_OPT'   ,'MOCUR'     ,'AXREF_FDM'            &
      ,'MC'         ,'MCCMFD'    ,'CMFDSET'   ,'DCMP'      ,'NODE_MAJOR'           &
      ,'AFSS'       ,'NODALOPT'  ,'MULTIGRID' ,'GRIDSTR'   ,1*''                   &
!  5  MATERIAL 3                                                                   &
      ,'MIXTURE'    ,'TRISO'    ,'INCLUDE'    ,27*''                               &
!  6  GEOM    24                                                                   &
      ,'HEX'        ,'NPINS'     ,'PITCH'     ,'AX_MESH'   ,'ALBEDO'               &
      ,'STRUCT'     ,'CELL'      ,'GAP'       ,'PIN'       ,'ASSEMBLY'             &
      ,'RAD_CONF'   ,'CAD_CELL'  ,'GAP_CELL'  ,'GAP_PIN'   ,'GAP_ASY'              &
      ,'RING_STRUCT','BAFFLE'    ,'BARREL'    ,'INCLUDE'   ,'HGAP'                 &
      ,'VESSEL'     ,'VYGORODKA' ,'HEXOPT'    ,'NTIG'      ,'BASECELL', 5*' '      &
!  7  TH    20                                                                     &
      ,'PIN_DIM'    ,'NRING_COND','EFF_DOPLT' ,'KCOND_FUEL','RHOCP_FUEL'           &
      ,'KCOND_CLAD' ,'RHOCP_CLAD','STEAM_TBL' ,'SIMPLE_TH' ,'CHGRP'                &
      ,'MATRA'      ,'FDM_COND'  ,'IHP'       ,'MODT_FILE' ,'FUELT_FILE'           &
      ,'COBRA_TF'   ,'NPIN_TH'   ,'PITCH_TH'  ,'HGAP_TH'   ,'ESCOT'                &
      ,10*' '                                                                      &
!  8  TRAN  6                                                                      &
      ,'TIME_STEP'  ,'EXPO_OPT'  ,'THETA'     ,'COND_RT'   ,'TIME_EDIT'            &
      ,'XS_CHANGE'  ,24*''                                                         &
!  9  DEPL   10                                                                    &
      ,'BU_TYPE'    ,'BURNUP'    ,'OPT_DEP'   ,'DEPLFILE'  ,'PC_OPT'               &
      ,'B10DEPL'    ,'COREFOLLOW','CORE_STATE','FLUX_INIT' ,'EFFXSDEPL'            &
      ,'NSUBSTEP'   ,'SOLVER'    ,18*''                                            &
! 10  BRANCH  8                                                                    &
      ,'BORONWORTH' ,'MTC'       ,'FTC'       ,'MAT'       ,'PPM'                  &
      ,'TM'         ,'TF'        ,'TMAT'      ,22*''                               &
! 11  EDIT   24                                                                    &
      ,'ISOTOPE'    ,'RSTFILE'   ,'DEPL_RST'  ,'FLUX'      ,'PINXS'                &
      ,'CSPOUT'     ,'EFFXS'     ,'BOUTP'     ,'FLUX_BOUTP','BSTEP_BOUTP'          &
      ,'ISOTOPE_BOUTP','GC_OPT'  ,'TMOD'      ,'TFUEL'     ,'RHO'                  & !--- BYS edit : options for GC_opt and FXR variation
      ,'FSRXS'      ,'BKLG'      ,'RAYGEN'    ,'DETAIL'    ,'NODETIME'             &
      ,'FXRMGMAC'   ,'SSPHOUT'   ,'EFFMAT0'   ,'PHIM'      ,6*''                   &
! 12  DECOUPLE 7                                                                   &
      ,'REF_PLN'    ,'REF_TEMP'  ,'PLN_MAP'   ,'XSFTN'     ,'DCPLITR'              &
      ,'DCPLCMFDITR','DCPLFBITR' ,23*''                                            &
! 13  PARALLEL  5                                                                  &
      ,'MOC_TRD'    ,'NODAL_TRD', 'CMFD_TRD'  ,'DEPL_TRD', 'AX_DCP', 25*''         &
! 14  VISUAL    2                                                                  &
      ,'VISMOD'     ,'FLUXOUT'   ,28*''                                            &
! 15  LP_SHF    7                                                                  &
      ,'CYCLEID'    ,'RSTLIB'   ,'PUL'       ,'SHF'       ,'CYCLE'                 &
      ,'RST'        ,'['        ,'PLN_MAP'   ,'RMV_BP'    ,21*''                   &
! 16  CNTLROD   9                                                                  &
      ,'CR_CELL'    ,'CR_ASYCONF','CR_BANK'  ,'CR_POS'     ,'CR_POSCHG'            &
      ,'CSP_FILE'   ,'CSP_MAP'   ,'CRMV_DOM' ,'CR_POSDAT'  ,21*''                  &
! 17                                                                               &
      ,'PINXS_FILE' ,'XS_MAP'    ,'CSP_BANK' ,'CSP_PLN'    ,26*''                  &
! 18                                                                               &
      ,'TIME_STEP'  ,'TIME_UNIT' ,'CORE_STATE',27*''                               &
! 19                                                                               &
      ,'NMAXOUTER'  ,'EIG_CONV'  ,'RES_CONV' ,'PSI_CONV'    ,26*''                 &
! 20                                                                               &
      ,'NMCP'       ,'MCP_LIST'  ,'MCP_FILE' ,'MCP_CALL'    ,'PLANE'               &
      ,'NMCP_PLANE' ,24*''                                                         &
! 21                                                                               &
      ,'CONTROL'    ,'CMAXOUTER'  ,'BORONTRACK', 'REL_CONV',26*''                  &
!--- CNJ Edit : Intel MKL Options                                                  &
! 22  MKL                                                                          &
      ,'DAVIDSON'   ,'BILU'       ,'AXTYPE'     ,'DCPL'     ,'SHIFT'               &
      ,'PCMFD'      ,'DIRECT'     ,'CHEBYSHEV'  ,'POLAR_XS' ,'GCCMFD'              &
      ,'GCSTR'      ,'SUPERPIN'   ,'SUBPLN'     ,17*''                             &
!--- CNJ Edit : CUDA Options                                                       &
! 23  CUDA                                                                         &
      ,'CU_CMFD'    ,'CU_GCCMFD'  ,'CU_GCSTR'   ,'CU_DCPL'   ,'CU_RBSOR'           &
      ,'CU_AXTYPE'  ,'CU_SPAI'    ,'CU_AAJ'     ,'CU_SUBPLN' ,'CU_SHIFT'           &
      ,20*''                                                                       &
! 24  Group Constants                                                              &
      ,'ASY_GCL'    ,'PIN_GCL'    ,'PDQ_GCL'    ,'NTIG_RST' ,'SPEC'                &
      ,'NFEWG'      ,'TMOD'       ,'TFUEL'      ,'RHO'      ,'BORON'               &
      ,'EFT'        ,'UNICOOL'    ,'EFT_BORON'  ,'EMT'      ,16*''                 &
! 25  NTIG_RESTART    3                                                            &
      ,'HM_MASS'    ,'FRESH_ASM'  ,'PUL'        , 27*''                            /
CONTAINS
FUNCTION FindBlockId(blockname)
implicit none
INTEGER                   :: FindBlockId
character(15), intent(in) ::blockname
INTEGER :: i

DO i = 1, nblock
  IF(blockname .eq. blocks(i)) exit
ENDDO
FindBlockId = i
IF(i .gt. nblock) FindBlockId=0
END FUNCTION

FUNCTION FindCardId(BlockId,cardname)
implicit none
INTEGER                   :: FindCardId
INTEGER, intent(in)       :: blockId
character(15)             :: cardname
INTEGER :: i
DO i = 1, ncard(blockid)
  IF(cardname.eq.cards(i,BlockId)) exit
ENDDO

FindCardId = i
IF(i .gt. ncard(blockid)) FindCardId=0
END FUNCTION

END module
