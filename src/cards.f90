module inputcards
implicit none
!character(1) DOT,BANG,BLANK,SLASH,AST
!parameter (DOT='.',BANG='!',BLANK=' ',SLASH='/',AST='*')
INTEGER, parameter :: mxncol=512, mxnfield=128
INTEGER, parameter :: mxcard=40, nblock=25
character*1 probe,sc(mxncol),sc80(80)
character*512  oneline
character*1024  longline
!character*10 cardname,blockname


logical IFfile
equivalence (probe,oneline)
equivalence (sc(1),oneline)

INTEGER ::  ncard(nblock)
data ncard/ 1,  2, 24, 31, 3,            & !CASEID    STATE      XSEC      OPTION     MATERIAL
           25, 25, 25, 12, 8,            & !GEOM      TH         TRAN      DEPL       BRANCH
           26,  7,  5, 2,  9,            & !EDIT      DECOUPLE   PARALLEL  VISUAL     LP_SHF  !--- BYS edit : # of EDIT card
           10,  4,  3, 5,  6,            & !CNTLROD   CUSPING    XEDYN     CONV       MCP_RESTART
            6, 15, 14, 14, 3/              !SUBCH_OP  MKL        CUDA      GCGEN      NTIG_RESTART
character(15):: blocks(nblock),cards(mxcard,nblock)

data blocks /'CASEID'    ,'STATE'     ,'XSEC'      ,'OPTION'    ,'MATERIAL'        &
            ,'GEOM'      ,'TH'        ,'TRAN'      ,'DEPL'      ,'BRANCH'          &
            ,'EDIT'      ,'DECOUPLE'  ,'PARALLEL'  ,'VISUAL'    ,'LP_SHF'          &
            ,'CNTLROD'   ,'CUSPING'   ,'XEDYN'     ,'CONV'      ,'MCP_RESTART'     &
            ,'SUBCH_OP'  ,'MKL'       ,'CUDA'      ,'GCGEN'     ,'NTIG_RESTART'    /

data cards/   'RESTART'   ,39*' '                                                  &
!  2  STATE   2                                                                    &
      ,'CORE_POWER' ,'TH_COND'   ,38*' '                                           &
!  3  XSEC   24                                                                    &
      ,'LIB_TYPE'   ,'GROUP_SPEC','FILE'      ,'BASE_MICRO','DNEUT_CHI'            &
      ,'NEUT_VELO'  ,'DNP_NGRP'  ,'DNP_BETA'  ,'DNP_LAMBDA','DUMMY'                &
      ,'BUCKLING'   ,'RLB'       ,'SLB'       ,'GC_SPEC'   ,'DEPLFILE'             &
      ,'PHL'        ,'RT'        ,'ED'        ,'MCXSFILE'  ,'RESO_OPT'             &
      ,'PMATRX'     ,'PSM'       ,'OTFRIF'    ,'RIF'       ,16*''                  &
!  4  OPTION  31                                                                   &
      ,'FEEDBACK'   ,'CMFD'      ,'RAY'       ,'MOC'       ,'ITER_LIM'             &
      ,'GRP_ALBEDO' ,'MPI_GROUPS','GC_CMFD'   ,'SEARCH'    ,'LKGSPLIT'             &
      ,'BORON'      ,'XENON'     ,'NODAL'     ,'CRITICAL'  ,'LINSRC'               &
      ,'FRR'        ,'SCAT_ORDER','BCR_OPT'   ,'MOCUR'     ,'AXREF_FDM'            &
      ,'MC'         ,'MCCMFD'    ,'CMFDSET'   ,'DCMP'      ,'NODE_MAJOR'           &
      ,'AFSS'       ,'NODALOPT'  ,'MULTIGRID' ,'GRIDSTR'   ,'SUBPLN'               &
      ,'POWERMODE'  , 9*''                                                         &
!  5  MATERIAL 3                                                                   &
      ,'MIXTURE'    ,'TRISO'    ,'INCLUDE'    ,37*''                               &
!  6  GEOM    25                                                                   &
      ,'HEX'        ,'NPINS'     ,'PITCH'     ,'AX_MESH'   ,'ALBEDO'               &
      ,'STRUCT'     ,'CELL'      ,'GAP'       ,'PIN'       ,'ASSEMBLY'             &
      ,'RAD_CONF'   ,'CAD_CELL'  ,'GAP_CELL'  ,'GAP_PIN'   ,'GAP_ASY'              &
      ,'RING_STRUCT','BAFFLE'    ,'BARREL'    ,'INCLUDE'   ,'HGAP'                 &
      ,'VESSEL'     ,'VYGORODKA' ,'HEXOPT'    ,'NTIG'      ,'BASECELL', 15*' '     &
!  7  TH    25                                                                     &
      ,'PIN_DIM'    ,'NRING_COND','EFF_DOPLT' ,'KCOND_FUEL','RHOCP_FUEL'           &
      ,'KCOND_CLAD' ,'RHOCP_CLAD','STEAM_TBL' ,'SIMPLE_TH' ,'CHGRP'                &
      ,'MATRA'      ,'FDM_COND'  ,'IHP'       ,'MODT_FILE' ,'FUELT_FILE'           &
      ,'COBRA_TF'   ,'NPIN_TH'   ,'PITCH_TH'  ,'HGAP_TH'   ,'ESCOT'                &
      ,'WRAPPER'    ,'URELX_FUEL','FRAC_DC'   ,'CH_CONF'   ,'AA_STH'               &
      , 15*' '                                                                        &
!  8  TRAN  25                                                                     &
      ,'TIME_STEP'  ,'EXPO_OPT'  ,'THETA'     ,'COND_RT'   ,'TIME_EDIT'            &
      ,'XS_CHANGE'  ,'KIN_BENCH' ,'AF_SRC'    ,'USERTHETA' ,'NMAXOUTER'            &
      ,'CONV_CMFD'  ,'NMAXCMFD'  ,'STEPREAC'  ,'METHOD'    ,'RES_CONV'             &
      ,'PSI_CONV'   ,'CORRECTOR' ,'ADJOINT'   ,'DYN_BENCH' ,'XS_NOISE'             &
      ,'XS_CNTLROD' ,'WORKING'   ,'DCY_HEAT'  ,'NOISE_SAMPLING', 'NNFSP', 15*''    &
!  9  DEPL   12                                                                    &
      ,'BU_TYPE'    ,'BURNUP'    ,'OPT_DEP'   ,'DEPLFILE'  ,'PC_OPT'               &
      ,'B10DEPL'    ,'COREFOLLOW','CORE_STATE','FLUX_INIT' ,'EFFXSDEPL'            &
      ,'NSUBSTEP'   ,'SOLVER'    ,28*''                                            &
! 10  BRANCH  8                                                                    &
      ,'BORONWORTH' ,'MTC'       ,'FTC'       ,'MAT'       ,'PPM'                  &
      ,'TM'         ,'TF'        ,'TMAT'      ,32*''                               &
! 11  EDIT   26                                                                    &
      ,'ISOTOPE'    ,'RSTFILE'   ,'DEPL_RST'  ,'FLUX'      ,'PINXS'                &
      ,'CSPOUT'     ,'EFFXS'     ,'BOUTP'     ,'FLUX_BOUTP','BSTEP_BOUTP'          &
      ,'ISOTOPE_BOUTP','GC_OPT'  ,'TMOD'      ,'TFUEL'     ,'RHO'                  & !--- BYS edit : options for GC_opt and FXR variation
      ,'FSRXS'      ,'BKLG'      ,'RAYGEN'    ,'DETAIL'    ,'NODETIME'             &
      ,'FXRMGMAC'   ,'SSPHOUT'   ,'EFFMAT0'   ,'PHIM'      ,'GEFFXS'               &
      ,'KERMA'      ,14*''                                                         &
! 12  DECOUPLE 7                                                                   &
      ,'REF_PLN'    ,'REF_TEMP'  ,'PLN_MAP'   ,'XSFTN'     ,'DCPLITR'              &
      ,'DCPLCMFDITR','DCPLFBITR' ,33*''                                            &
! 13  PARALLEL  5                                                                  &
      ,'MOC_TRD'    ,'NODAL_TRD', 'CMFD_TRD'  ,'DEPL_TRD', 'AX_DCP', 35*''         &
! 14  VISUAL    2                                                                  &
      ,'VISMOD'     ,'FLUXOUT'   ,38*''                                            &
! 15  LP_SHF    9                                                                  &
      ,'CYCLEID'    ,'RSTLIB'   ,'PUL'       ,'SHF'       ,'CYCLE'                 &
      ,'RST'        ,'['        ,'PLN_MAP'   ,'RMV_BP'    ,31*''                   &
! 16  CNTLROD   10                                                                 &
      ,'CR_CELL'    ,'CR_ASYCONF','CR_BANK'  ,'CR_POS'     ,'CR_POSCHG'            &
      ,'CSP_FILE'   ,'CSP_MAP'   ,'CRMV_DOM' ,'CR_POSDAT'  ,'CR_DECUSP'            &
      ,30*''                                                                       &
! 17  CUSPING   4                                                                  &
      ,'PINXS_FILE' ,'XS_MAP'    ,'CSP_BANK' ,'CSP_PLN'    ,36*''                  &
! 18  XEDYN     3                                                                  &
      ,'TIME_STEP'  ,'TIME_UNIT' ,'CORE_STATE',37*''                               &
! 19  CONV      4                                                                  &
      ,'NMAXOUTER'  ,'EIG_CONV'  ,'RES_CONV' ,'PSI_CONV'    ,'DECUSP_CONV'         &
      ,35*''                                                                       &
! 20  MCP_RESTART  6                                                               &
      ,'NMCP'       ,'MCP_LIST'  ,'MCP_FILE' ,'MCP_CALL'    ,'PLANE'               &
      ,'NMCP_PLANE' ,34*''                                                         &
! 21  SUBCH_OP  6                                                                  &
      ,'CONTROL'    ,'CMAXOUTER'  ,'BORONTRACK', 'REL_CONV' ,'COURANT'             &
      ,'SBCH_OUT'   ,34*''                                                         &
!--- CNJ Edit : Intel MKL Options                                                  &
! 22  MKL       14                                                                 &
      ,'DAVIDSON'   ,'BILU'       ,'AXTYPE'     ,'DCPL'     ,'SHIFT'               &
      ,'PCMFD'      ,'DIRECT'     ,'CHEBYSHEV'  ,'POLAR_XS' ,'GCCMFD'              &
      ,'GCSTR'      ,'SUPERPIN'   ,'SUBPLN'     ,'DEPL'     ,'CMFD'                &
      ,25*''                                                                       &
!--- CNJ Edit : CUDA Options                                                       &
! 23  CUDA      14                                                                 &
      ,'CU_CMFD'    ,'CU_GCCMFD'  ,'CU_GCSTR'   ,'CU_DCPL'   ,'CU_RBSOR'           &
      ,'CU_AXTYPE'  ,'CU_SPAI'    ,'CU_AAJ'     ,'CU_SUBPLN' ,'CU_SHIFT'           &
      ,'CU_DEPL'    ,'CU_XS'      ,'CU_SUPERPIN','CU_PWDIST' ,26*''                &
! 24  Group Constants  14                                                          &
      ,'ASY_GCL'    ,'PIN_GCL'    ,'PDQ_GCL'    ,'NTIG_RST' ,'SPEC'                &
      ,'NFEWG'      ,'TMOD'       ,'TFUEL'      ,'RHO'      ,'BORON'               &
      ,'EFT'        ,'UNICOOL'    ,'EFT_BORON'  ,'EMT'      ,26*''                 &
! 25  NTIG_RESTART    3                                                            &
      ,'HM_MASS'    ,'FRESH_ASM'  ,'PUL'        , 37*''                            /
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
