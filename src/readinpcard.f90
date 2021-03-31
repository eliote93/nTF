#include <defines.h>
module readinpcard
use param
use files,          only : filename,    InputFileIdx,                  ModTFileIdx, FuelTFileIdx, &
                           io5,         io8
use ioutil,         only : terminate,   toupper,       IFnumeric,      nfields,     icolfield,    &
                           fndchara,    GetFn,         message
use inputcards ,    only : oneline,     probe,         mxcard,         nblock,      FindBlockId,  &
                           FindCardId,  cards,         blocks, longline
use allocs
implicit none
character(15),private   :: cardname, blockname, astring
character(512),private :: dataline
INTEGER,private :: nLineField

contains

SUBROUTINE ReadStateCard(InDev, OutDev)
USE PARAM
USE CNTL, ONLY : nTracerCntl
USE PE_MOD, ONLY : PE
IMPLICIT NONE
INTEGER :: InDev, OutDev
INTEGER :: idcard
INTEGER,parameter :: idblock = 2
INTEGER :: IntTemp(100)
INTEGER :: i, j, k
REAL :: REALTemp(100)

!REAL::w
!CHARACTER(10) :: warg

LOGICAL :: MASTER

!call getarg(2,warg)
!READ(warg, *) w
!print *, w
master = PE%master
DO while(TRUE)
  read(indev,'(a256)') oneline
  IF(master) CALL message(io8,FALSE,FALSE,oneline)
  IF(probe.eq.BANG) CYCLE;     IF(probe.eq.POUND) CYCLE
  IF(oneline.eq.BLANK) CYCLE;  IF(IFnumeric(oneline)) CYCLE
  IF(probe.eq.DOT) EXIT;       IF(probe.eq.SLASH) EXIT
  READ(oneline,*) cardname; CALL toupper(cardname)
  idcard = FindCardId(idblock,cardname); IF(idcard .eq. 0) EXIT
  nLineField = nfields(oneline)-1
  SELECT CASE(idcard)
    CASE(1) !'CORE_POWER'
      READ(oneline, *) astring, nTracerCntl%PowerLevel
      nTracerCntl%PowerLevel = nTracerCntl%PowerLevel * 0.01_8
    CASE(2) !TH_COND
      READ(oneline, *) astring, (REALTemp(j), j = 1, nLineField)
      nTracerCntl%PowerFA = RealTemp(1); nTracerCntl%TempInlet = RealTemp(2)
      nTracerCntl%Pexit = RealTemp(3); nTracerCntl%fmdotfa = RealTemp(4)
      nTracerCntl%PowerFA = nTracerCntl%PowerFA * 1.0E6_8  !MW => W
      nTracerCntl%Pexit = nTracerCntl%Pexit * 1.0E6_8      !MPA => Pa

      nTracerCntl%PowerFA = nTracerCntl%PowerFA
      nTracerCntl%fmdotfa = nTracerCntl%fmdotfa
  END SELECT
ENDDO
Backspace(indev); IF(Master) Backspace(outdev)
END SUBROUTINE

SUBROUTINE ReadMaterialCard(indev, outdev)
USE PE_MOD,    ONLY : PE
USE CNTL, ONLY : nTracerCntl
IMPLICIT NONE
INTEGER           :: indev,outdev
INTEGER           :: idcard
INTEGER,parameter :: idblock = 5
INTEGER :: IntTemp(100)
INTEGER :: i,k
REAL :: REALTemp(100)
LOGICAL :: Master

Master = PE%master
DO while(TRUE)
  read(indev,'(a256)') oneline
  IF(Master) CALL message(io8,FALSE,FALSE,oneline)
  IF(probe.eq.BANG) CYCLE;     IF(probe.eq.POUND) CYCLE
  IF(oneline.eq.BLANK) CYCLE;  IF(IFnumeric(oneline)) CYCLE
  IF(probe.eq.DOT) EXIT;       IF(probe.eq.SLASH) EXIT
  READ(oneline,*) cardname; CALL toupper(cardname)
  idcard = FindCardId(idblock,cardname); IF(idcard .eq. 0) EXIT
  nLineField = nfields(oneline)-1
  SELECT CASE(idcard)
    CASE(1) !'CORE_POWER'
      i=icolfield(oneline,2)
      dataline = oneline(i:256)
      IF(nTracerCntl%lISOTXS) THEN ! for fast reactor...
        CALL ProcMat_ISOTXS(dataline, indev, outdev)
      ELSE
        CALL ProcMat(dataline, indev, outdev)
      END IF
    CASE(2)
      CONTINUE
  END SELECT
ENDDO
Backspace(indev); IF(Master) Backspace(outdev)
END SUBROUTINE

SUBROUTINE ReadGeomCard(indev,outdev)
use geom,      only : nz,        hz,           HzInv,                                                 &
                      lCoreAng,  CellPitch,    AsyPitch,    lgap,                                     &
                      nCellX0,   nCellType0,  nCellType,   CellInfo,      nPinType0,    nPinType,     &
                      PinInfo
USE CNTL,      ONLY : nTracerCntl
USE PE_MOD,         ONLY : PE
implicit none
INTEGER           :: indev,outdev
INTEGER           :: idcard
INTEGER,parameter :: idblock = 6
INTEGER :: IntTemp(100)
INTEGER :: i,k
REAL :: REALTemp(100)
LOGICAL :: Master
LOGICAL :: lCellInit, lBaseCellInit
Master = PE%master

lCellInit = .FALSE.
lBaseCellInit = .FALSE.

DO while(TRUE)
  read(indev,'(a)') oneline
  IF(Master) CALL message(io8,FALSE,FALSE,oneline)
  IF(probe.eq.BANG) cycle;     IF(probe.eq.POUND) cycle
  IF(oneline.eq.BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe.eq.DOT) exit;       IF(probe.eq.SLASH) exit
  read(oneline,*) cardname; CALL toupper(cardname)
  idcard = FindCardId(idblock,cardname); IF(idcard .eq. 0) exit
  nLineField = nfields(oneline)-1

  select case(idcard)
    case(1)   !HEX

    case(2)   !NPINS
      !Already read in the scaninput
    case(3)   !PITCH
      !Already read in the scaninput
      READ(oneline,*) ASTRING, CellPitch
      AsyPitch = CellPitch*nCellX0
      IF(nLineField .gt. 1) READ(oneline,*) ASTRING, CellPitch, AsyPitch
    case(4)   !AX_MESH
      CALL dmalloc(hz,nz); CALL dmalloc(HzInv, nz)
      read(oneline,*) astring, (hz(k),k=1,nz)
      HzInv = 1./Hz
      CONTINUE
    case(5)   !ALBEDO

    case(6)   !STRUCT

    case(7)   !CELL
      i=icolfield(oneline,2)
      dataline = oneline(i:512) ! BYS edit 160708
      !--- JSR Edit : nTIG Restart
      IF (nTracerCntl%lnTIGRst) THEN
        CALL ReadNInit_Cell_nTIG(dataline, FALSE, FALSE, lCellInit)
      ELSE
        CALL ReadNInit_Cell(dataline, FALSE, FALSE, lCellInit)
      ENDIF
    case(8)   !GAP
      i=icolfield(oneline,2)
      dataline = oneline(i:512)
      !--- JSR Edit : nTIG Restart
      IF (nTracerCntl%lnTIGRst) THEN
        CALL Init_GapCell_nTIG(DataLine, lBaseCellInit)
      ELSE
        CALL Init_GapCell(DataLine, lCellInit)
      ENDIF
    case(9)   !PIN
      i=icolfield(oneline,2)
      dataline = oneline(i:512) ! BYS edit 160708
      CALL ReadNinit_Pin(dataline, FALSE)
    case(10)  !ASEMBLY
      IF(lGap) CALL Init_GapPin
      i=icolfield(oneline,2)
      dataline = oneline(i:512)
      CALL ReadNInit_ASSEMBLY(indev,dataline)
    case(11)  !RAD_CONF
      CALL ReadNInit_Core(indev)
    case(12)  !CAD_CELL
      i=icolfield(oneline,2)
      dataline = oneline(i:512)
      CALL ReadNInit_Cell(dataline, TRUE, FALSE, lCellInit)
    CASE(13)
      i=icolfield(oneline,2)
      dataline = oneline(i:512)
      !CALL ReadNInit_UsrGapCell(Dataline)
      CALL ReadNInit_Cell(dataline, FALSE, TRUE, lCellInit)
      !CALL TERMINATE('Underconstruction : Baffle Option')
    CASE(14) !ASY_GAP
      i=icolfield(oneline,2)
      dataline = oneline(i:512)
      CALL ReadNinit_Pin(dataline, TRUE)
    CASE(15) !
      i=icolfield(oneline,2)
      dataline = oneline(i:512)
      CALL ReadNInit_AsyGapConf(indev, dataline)
    CASE(16) ! RING_STRUCT
      i=icolfield(oneline,2)
      dataline = oneline(i:512)
      CALL ReadNInit_RingStruct(dataline)
    CASE(17) ! BAFFLE

    CASE(18) ! BARREL
      i=icolfield(oneline,2)
      dataline = oneline(i:512)
      CALL ReadNInit_Barrel(dataline)
    case(20)   !GAP
      i=icolfield(oneline,2)
      dataline = oneline(i:512)
      !CALL  Init_HGapCell(DataLine, lCellInit)
    !--- JSR Edit : nTIG Restart
    case(25)   !BASECELL
      i=icolfield(oneline,2)
      dataline = oneline(i:512) ! BYS edit 160708
      CALL ReadNInit_BaseCell_nTIG(dataline, FALSE, FALSE, lBaseCellInit)
  END select
ENDDO
Backspace(indev); IF(Master) Backspace(io8)

END SUBROUTINE

SUBROUTINE ReadTHCard(InDev, OutDev)
USE PARAM
USE Geom,           ONLY : Core, CellPitch, AsyPitch
USE TH_Mod,         ONLY : ThOpt, ThVar
USE CNTL,           ONLY : nTracerCntl
USE PE_MOD,         ONLY : PE
USE ioutil,         ONLY: terminate
USE SubChCoupling_mod,    ONLY: is_coupled, CodeName
USE Anderson_Acceleration_SIMPLETH, ONLY : m_AA
IMPLICIT NONE

INTEGER :: InDev, OutDev
INTEGER           :: idcard
INTEGER,parameter :: idblock = 7
INTEGER :: IntTemp(100)
INTEGER :: i,k,nb, nfn, AA_STH
REAL :: REALTemp(100)
LOGICAL :: Master

Master = PE%master
ThVar%ChannelPitch = CellPitch * epsm2
ThVar%AsyPitch = AsyPitch * epsm2

ThVar%nAsych=REAL(Core%nAsyCell - Core%nAsyGT)       !Number of Fuel Channel
ThVar%nAsyGT = REAL(Core%nAsyGT)                     !Number of Guide Tube
ThVar%nzth = Core%nz                           !Number of Axial Node
nfn = 0
DO while(TRUE)
  read(indev,'(a256)') oneline
  IF(master) CALL message(io8,FALSE,FALSE,oneline)
  IF(probe .EQ. BANG) cycle;     IF(probe .EQ. POUND) cycle
  IF(oneline .EQ. BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe .EQ. DOT) exit;       IF(probe.eq.SLASH) exit
  read(oneline,*) cardname; CALL toupper(cardname)
  idcard = FindCardId(idblock,cardname); IF(idcard .eq. 0) exit
  nLineField = nfields(oneline)-1

  !AA_STH = 1
  !m_AA = AA_STH
  !ThOpt%AA_STH = .true.
  SELECT CASE(idcard)
    CASE(1)  ! PIN_DIM
      READ(oneline, *) astring, (THOpt%PinDim(i), i = 1, nLineField)
      IF (nLineField .EQ. 4) THEN
        THOpt%PinDim(5) = THOpt%PinDim(4) - THOpt%PinDim(3)
      endif
    CASE(2)  ! NRING_COND
      READ(oneline, *) astring, THOpt%nrpellet
    CASE(3)  ! EFF_DOPLT

    CASE(4)  ! KCOND_FUEL
      !READ(oneline, *) astring, (THOpt%kFUelCorrelation(i), i=0, 5)
      READ(oneline, *) astring, THOpt%kFuelModel
    CASE(5)  ! RHOCP_FUEL
      !READ(oneline, *) astring, (THOpt%CpFuelCorrelation(i), i=0, 3)
      READ(oneline, *) astring, THOpt%CpFuelModel, THOpt%RhoFuel
    CASE(6)  ! KCOND_CLAD
      !READ(oneline, *) astring, (THOpt%KCladCorrelation(i), i=0, 3)
      READ(oneline, *) astring, THOpt%KCladModel
    CASE(7)  ! RHOCP_CLAD
      !READ(oneline, *) astring, (THOpt%CpCladCorrelation(i), i=0, 3)
      READ(oneline, *) astring, THOpt%CpCladModel, THOpt%RhoClad
    CASE(8)  ! STEAM_TBL

    CASE(9)  ! SIMPLE_TH
      READ(oneline, *) astring, nTracerCntl%lSimpleTH
    CASE(10) !ChGrp
      READ(oneline, *) astring, nTracerCntl%ThCh_mod
    CASE(11) !MATRA
      READ(oneline, *) astring, nTracerCntl%lMatra
    CASE(12) !FDM_COND
      READ(oneline, *) astring, nTracerCntl%lFuelCondFDM
    CASE(13) !IHP
      READ(oneline, *) astring, nTracerCntl%LIHP
    CASE(14) ! Moderator
      !READ(oneline, *) astring,
      CALL GetFn(oneline, 2, filename(ModTFileIdx))
      nfn = nfn + 1
    CASE(15) ! Fuel Temperature
      !READ(oneline, *) astring,
      CALL GetFn(oneline, 2, filename(FuelTFileIdx))
      nfn = nfn + 1
    case(16) ! COBRA_TF
      READ(oneline,*) astring, is_coupled
      CodeName='CTF'
    case(17)
      read(oneline, *) astring, ThVar%nAsych, ThVar%nAsyGT
    case(18)
      IF(nLineField .EQ. 2) THEN
        read(oneline, *) astring, ThVar%ChannelPitch, ThVar%AsyPitch
      ELSE
        read(oneline, *) astring, ThVar%ChannelPitch, ThVar%AsyPitch, ThVar%HAct
        ThVar%lhact = .TRUE.
      END IF

      ThVar%ChannelPitch = ThVar%ChannelPitch* epsm2
      ThVar%AsyPitch=ThVar%AsyPitch*epsm2
    case(19) ! Hgap_TH
      !READ(oneline,*) astring, ThOpt%hgap
      READ(oneline,*) astring, ThOpt%hGapModel
    case(20) ! ESCOT
      READ(oneline,*) astring, is_coupled
      CodeName='ESCOT'
    case(23) ! FRAC_DC
      READ(oneline,*) astring, ThVar%FracDC
      !ThVar%FracDF = 1._8 - ThVar%FracDC
    case(24) ! ch_conf
      nTracerCntl%lthchconf = .TRUE.
      CALL ReadCh_conf(indev, ThVar)
    case(25) !AA_STH
      READ(oneline,*) astring, AA_STH
      if (AA_STH < 1) then
        ThOpt%AA_STH = .false.
        m_AA = AA_STH
      else
        ThOpt%AA_STH = .true.
        m_AA = AA_STH
      endif
  END SELECT
ENDDO

IF(nfn .EQ. 2) nTracerCntl%lUserDefTemp = .TRUE.
IF(nTracerCntl%lIHP) nTracerCntl%lFuelCondFDM = .TRUE.
Backspace(indev); IF(Master) Backspace(io8)
END SUBROUTINE

SUBROUTINE ReadCh_conf(indev, ThVar)
USE PARAM
USE TYPEDEF
USE PE_Mod,         ONLY : PE
USE files,          ONLY : io5, io8
USE inputcards,     ONLY : oneline, probe
USE ioutil,         ONLY : toupper, IFnumeric, nfields, message, terminate
USE geom,           ONLY : Core
IMPLICIT NONE
INTEGER, INTENT(IN) :: indev
TYPE(ThVar_Type) :: ThVar

INTEGER :: ja
INTEGER :: i, j, jfr, jto
INTEGER :: nFieldsLine

READ(oneline, *) astring, ThVar%nChType
ALLOCATE(ThVar%THCh(ThVar%nChType))
DO i = 1, ThVar%nChType
  READ(indev, '(a256)') oneline
  IF(PE%Master) CALL message(io8,FALSE,FALSE,oneline)
  READ(oneline, *)  astring, ThVar%THCh(i)%nAsyCh, ThVar%THCh(i)%nAsyGT, ThVar%THCh(i)%ChannelPitch, ThVar%THCh(i)%AsyPitch, ThVar%THCh(i)%hact
  ThVar%THCh(i)%ChannelPitch = ThVar%THCh(i)%ChannelPitch * epsm2
  ThVar%THCh(i)%AsyPitch = ThVar%THCh(i)%AsyPitch * epsm2
  ThVar%THCh(i)%hact = ThVar%THCh(i)%hact * epsm2
END DO

ALLOCATE(Core%ThChMap(Core%nxya))

ja = 0; jfr = 1
DO WHILE(TRUE)
  READ(indev, '(a256)') oneline
  IF(PE%Master) CALL message(io8,FALSE,FALSE,oneline)
  IF(probe.eq.BANG) cycle;     IF(probe.eq.POUND) cycle
  ja = ja + 1
  nFieldsLine = nfields(oneline)
  jto = jfr + nFieldsLine-1
  READ(oneline, *) (Core%THChMap(j), j = jfr, jto)
  jfr = jto + 1
  IF(ja .EQ. Core%nya) EXIT
END DO

END SUBROUTINE

SUBROUTINE ReadTranCard(indev, outdev)
USE PARAM
USE CNTL,         ONLY : nTracerCntl
USE Tran_mod,     ONLY : TranCntl,     XsCHANGE,      XsNoise,    XsCntlRod
USE PE_MOD,        ONLY : PE

INTEGER           :: indev,outdev
INTEGER           :: idcard
INTEGER,parameter :: idblock = 8

INTEGER :: n, m, l, k
INTEGER :: i, ipos(100)
CHARACTER(10) :: optfield
LOGICAL :: Master

CHARACTER(15) :: working_card

n = 0
l = 0
k = 0
master = PE%master
TranCntl%Tstep_inp = 0._8
TranCntl%Tdiv_inp = 0._8
DO while(TRUE)
  read(indev,'(a256)') oneline
  IF(master) CALL message(io8,FALSE,FALSE,oneline)
  IF(probe .EQ. BANG) cycle;     IF(probe .EQ. POUND) cycle
  IF(oneline .EQ. BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe .EQ. DOT) exit;       IF(probe.eq.SLASH) exit
  read(oneline,*) cardname; CALL toupper(cardname)
  idcard = FindCardId(idblock,cardname); IF(idcard .eq. 0) exit
  nLineField = nfields(oneline)-1

  SELECT CASE(idcard)
    CASE(1)  !Time Step
      CALL FndChara(oneline, ipos, m, SLASH)
      m = m + 1; TranCntl%Tstep_inp(0) = m
      !IF(MOD(nLineField, 2) .NE. 0) CALL TERMINATE('Not Proper TIME_STEP INPUT')
      i = 1
      READ(oneline, *) astring, TranCntl%Tstep_inp(i), TranCntl%Tdiv_inp(i)
      DO i = 2, m
        READ(oneline(ipos(i-1)+1:256), *) TranCntl%Tstep_inp(i), TranCntl%Tdiv_inp(i)
      ENDDO
      !READ(oneline, *) astring, (TranCntl%Tstep_inp(i), TranCntl%Tdiv_inp(i), i = 1, nLineField)
    CASE(2)  !EXPO_OPT
    CASE(3)  !THETA
    CASE(4)  !COND_RT
      READ(oneline, *) astring, TranCntl%lCondiMOC
    CASE(5)  !TIME_EDIT
      TranCntl%nTWriteOut = nLineField
      READ(oneline, *) astring, TranCntl%TWriteOut(1:nLineField)
    CASE(6)  !XS_CHANGE
      n = n + 1
      XsChange(n)%lUse = .TRUE.
      !READ(oneline, *) astring, XsChange(n)%iso0, XsChange(n)%iso1, XsChange(n)%tbeg, XsChange(n)%tend
      CALL ReadXsChange(Oneline, XsChange(n), TranCntl%lCusping)
    CASE(7)  !KIN_BENCH
    CASE(8)  !Af_Src
    CASE(9) !USERTHETA
      READ(oneline, *) astring, TranCntl%Theta
    CASE(10) !nMAXouter
      READ(oneline, *) astring, TranCntl%nMaxOuter
    CASE(11) !CONV_CMFD
      READ(oneline, *) astring, TranCntl%cmfd_res_conv
    CASE(12) !nMAXCMFD
      READ(oneline, *) astring, TranCntl%nMaxCmfd
    CASE(13) !lStepFunc
      !READ(oneline, *) astring, TranCntl%lStepFunc
      READ(oneline, *) astring, TranCntl%lStepImplicit
    CASE(14) !Method
      READ(oneline, *) astring, TranCntl%TD
      IF(TranCntl%TD .EQ. 1) THEN !CN
        TranCntl%lExpTrsf = .FALSE.; TranCntl%lAdpTheta = .FALSE.
        TranCntl%MOC_BDF = .FALSE.;  TranCntl%CMFD_BDF = .FALSE.
        TranCntl%lSCM = .FALSE.;     TranCntl%lAM = .FALSE.
      ELSEIF(TranCntl%TD .EQ. 2) THEN !CNET
        TranCntl%lExpTrsf = .TRUE.; TranCntl%lAdpTheta = .FALSE.
        TranCntl%MOC_BDF = .FALSE.;  TranCntl%CMFD_BDF = .FALSE.
        TranCntl%lSCM = .FALSE.;     TranCntl%lAM = .FALSE.
      ELSEIF(TranCntl%TD .EQ. 3) THEN !AT
        TranCntl%lExpTrsf = .FALSE.; TranCntl%lAdpTheta = .TRUE.
        TranCntl%MOC_BDF = .FALSE.;  TranCntl%CMFD_BDF = .FALSE.
        TranCntl%lSCM = .FALSE.;     TranCntl%AdpThetaMethod = .TRUE.
        TranCntl%lAM = .FALSE.
        TranCntl%AdpThetaStep = 1
      ELSEIF(TranCntl%TD .EQ. 4) THEN !BDF
        TranCntl%lExpTrsf = .FALSE.; TranCntl%lAdpTheta = .FALSE.
        TranCntl%MOC_BDF = .TRUE.;  TranCntl%CMFD_BDF = .TRUE.
        TranCntl%lSCM = .FALSE.;     TranCntl%lAM = .FALSE.
      ELSEIF(TranCntl%TD .EQ. 5) THEN !SCM
        TranCntl%lExpTrsf = .FALSE.; TranCntl%lAdpTheta = .FALSE.
        TranCntl%MOC_BDF = .FALSE.;  TranCntl%CMFD_BDF = .FALSE.
        TranCntl%lSCM = .TRUE.;     TranCntl%lAM = .FALSE.
      ELSEIF(TranCntl%TD .EQ. 6) THEN !SCM_Prec
        TranCntl%lExpTrsf = .FALSE.; TranCntl%lAdpTheta = .FALSE.
        TranCntl%MOC_BDF = .FALSE.;  TranCntl%CMFD_BDF = .FALSE.
        TranCntl%lSCM = .TRUE.;     TranCntl%lAM = .FALSE.
        TranCntl%lSCM_Prec = .TRUE.
      ELSEIF(TranCntl%TD .EQ. 7) THEN !AM3
        TranCntl%lExpTrsf = .FALSE.; TranCntl%lAdpTheta = .FALSE.
        TranCntl%MOC_BDF = .TRUE.;  TranCntl%CMFD_BDF = .TRUE.
        TranCntl%lSCM = .FALSE.;     TranCntl%lAM = .TRUE.
      ELSEIF(TranCntl%TD .EQ. 8) THEN !AM3ET
        TranCntl%lExpTrsf = .TRUE.; TranCntl%lAdpTheta = .FALSE.
        TranCntl%MOC_BDF = .TRUE.;  TranCntl%CMFD_BDF = .TRUE.
        TranCntl%lSCM = .FALSE.;     TranCntl%lAM = .TRUE.
      ELSE
        TranCntl%TD = 2
      ENDIF
    CASE(15) !RES_CONV
      READ(oneline, *) astring, TranCntl%res_conv
    CASE(16) !PSI_CONV
      READ(oneline, *) astring, TranCntl%psi_conv
    CASE(17) !CORRECTOR
      READ(oneline, *) astring, TranCntl%lCorrector, TranCntl%Tdiv_corrector
    CASE(18) !ADJOINT
      IF(nLineField .EQ. 1) THEN
        READ(oneline, *) astring, nTracerCntl%lAdjoint
      ELSE
        READ(oneline, *) astring, nTracerCntl%lAdjoint, nTracerCntl%lCMFDAdj
      END IF
    CASE(19) !Dyn_Bench
    CASE(20) !NOISE
      l = l+1
      CALL ReadXsNoise(Oneline, XsNoise(l))
    CASE(21) !XS_CNTLROD: Options for Control Rod cusping
      k = k + 1
      CALL ReadXsCntlRod(Oneline, XsCntlRod(k), TranCntl%lCusping)
    CASE(22) !WORKING: Options for a debugging
      DO WHILE(TRUE)
        READ(indev, '(a256)') oneline
        IF(probe .EQ. DOT) exit;       IF(probe.eq.SLASH) exit
        IF(master) CALL message(io8,FALSE,FALSE,oneline)
        nLineField = nfields(oneline)-1
        read(oneline,*) astring, working_card; CALL toupper(working_card)
        IF(working_card .EQ. 'PKE_GUESS') THEN
          IF(nLineField .EQ. 2) THEN
            READ(oneline, *) astring, astring, TranCntl%lGuess
          ELSE
            READ(oneline, *) astring, astring, TranCntl%lGuess, TranCntl%Tdiv_corrector
          END IF
        ELSE IF(working_card .EQ. 'TRAN_PCR') THEN
          READ(oneline, *) astring, astring, TranCntl%PCRtype
        ELSE IF(working_card .EQ. 'PCQS_ITER') THEN
          READ(oneline, *) astring, astring, TranCntl%lPCQSIter, TranCntl%lIQS
        ELSE IF(working_card .EQ. 'ADAPTIVET') THEN
          READ(oneline, *) astring, astring, TranCntl%lAdptT, TranCntl%Tend, TranCntl%Tdiv_inp(1)
        ELSE IF(working_card .EQ. 'IQSAA') THEN
          READ(oneline, *) astring, astring, TranCntl%lIQSAA, TranCntl%IQSAA_m
          TranCntl%IQSAA_m = TranCntl%IQSAA_m+1
          ALLOCATE(TranCntl%IQSAA_x(2, TranCntl%IQSAA_m), TranCntl%IQSAA_g(2, TranCntl%IQSAA_m))
        ELSE IF(working_card .EQ. 'FIXTMP_RW') THEN
          READ(oneline, *) astring, astring, TranCntl%lfixtmprw
        ELSE IF(working_card .EQ. 'WORKING_END') THEN
          EXIT
        END IF
      END DO
    CASE(23) ! DCY_HEAT: Power Calculation with Decay Heat
      READ(oneline,*) astring, nTracerCntl%lDcyHeat
    CASE(24) ! NOISE_SAMPLING : SAMPLING For Noise Analysis
      TranCntl%lNNSampling = .TRUE.
      READ(oneline,*) astring, TranCntl%SBeg, TranCntl%SEnd, TranCntl%SPeriod
    CASE(25) ! NNFSP: Neutron Noise Fixed Source Problem
      TranCntl%nfreq = nLineField
      READ(oneline,*) astring, (TranCntl%freq_inp(i), i = 1, nLineField)
      nTracerCntl%lproblem = lNNFSP
    END SELECT
ENDDO
IF(nTracerCntl%lProblem .NE. lNNFSP) nTracerCntl%lProblem = lTransient
nTracerCntl%lTranOn = .TRUE.
TranCntl%nchange = n
DO i = 1, n
  IF(XsChange(i)%lCusping) THEN
    nTracerCntl%lDecusp = .TRUE.
    EXIT
  END IF
END DO
Backspace(indev); IF(Master) Backspace(outdev)
END SUBROUTINE

SUBROUTINE ReadOptionCard(indev,outdev)
USE PARAM
USE RAYS,   ONLY : RayInfo
USE GEOM,   ONLY : ng, Core
USE CNTL,   ONLY : nTracerCntl
USE ItrCNTL_mod,      ONLY : ItrCntl
USE GEOM,   ONLY : nSubPlane, maxhzfm
USE PE_MOD, ONLY : PE
USE IOUTIL, ONLY : GetFn
USE FILES,  ONLY : filename,       XeEqFileIdx
USE InitMC_MOD
IMPLICIT NONE
INTEGER           :: indev,outdev
INTEGER           :: idcard
INTEGER,parameter :: idblock = 4
INTEGER :: IntTemp(100)
INTEGER :: i,k,nb,opt
REAL :: REALTemp(100)
CHARACTER(10) :: optfield
CHARACTER(100) :: charTemp
LOGICAL :: Master

master = PE%master
DO while(TRUE)
  read(indev,'(a256)') oneline
  IF(master) CALL message(io8,FALSE,FALSE,oneline)
  IF(probe .EQ. BANG) cycle;     IF(probe .EQ. POUND) cycle
  IF(oneline .EQ. BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe .EQ. DOT) exit;       IF(probe.eq.SLASH) exit
  read(oneline,*) cardname; CALL toupper(cardname)
  idcard = FindCardId(idblock,cardname); IF(idcard .eq. 0) exit
  nLineField = nfields(oneline)-1

  SELECT CASE(idcard)
    CASE(1)    !FEEDBACK
      READ(oneline, *) astring, nTracerCntl%lFeedBack
    CASE(2)    !CMFD
      READ(oneline, *) astring, ItrCntl%CmfdItrCntl%nIterMax, ItrCntl%CmfdItrCntl%nIterMin, itrcntl%CMFDItrCntl%ncmfdpnodal, ItrCntl%CmfdItrCntl%ConvCrit
    CASE(3)    !RAY
      READ(oneline, *) astring, RealTemp(1), (IntTemp(i), i = 1, nLineField - 1)
      IF(nLineField .GE. 1 .AND. RealTemp(1) .ne. 0)  RayInfo%Del = RealTemp(1)
      IF(nLineField .GE. 2) RayInfo%nAziAngle = IntTemp(1) * 2
      IF(nLineField .GE. 3) RayInfo%nPolarAngle = IntTemp(2)
      CONTINUE
    CASE(4)    !CONV_CRIT
      READ(oneline, *) astring, ItrCntl%MocItrCntl%nitermax
    CASE(5)    !ITER_LIM

    CASE(6)   !--- CNJ Edit : Group-wise Albedo (BYS Request)
      READ(oneline, *) astring, charTemp
      nTracerCntl%lGroupAlbedo = TRUE
    CASE(7)    !MPI_GROUPS

    CASE(9)    !SEARCH
      READ(oneline, *) astring, nTracerCntl%lSearch
      !IF(nTracerCntl%lSearch) nTracerCntl%lBoronSearch
      OptField = 'BORON'
      IF(nTracerCntl%lSearch) THEN
        IF(nLinefield .GT. 1) THEN
          READ(oneline, *) astring, nTracerCntl%lSearch, nTracerCntl%Target_eigv
        ENDIF
        IF(nLinefield .GT. 2) THEN
          READ(oneline, *) astring, nTracerCntl%lSearch, nTracerCntl%Target_eigv, optfield
        ENDIF
        !Select Reactivity Search
        CALL TOUPPER(OptField)
        SELECT CASE(OptField)
          CASE('BORON')
            nTracerCntl%lBoronSearch = .TRUE.; nTracerCntl%lCrSearch = .FALSE.
          CASE('CR')
            nTracerCntl%lBoronSearch = .FALSE.; nTracerCntl%lCrSearch = .TRUE.
          CASE DEFAULT
            CALL TERMINATE('Wrong Input Arguments in SEARCH Card')
        END SELECT
      ENDIF
    CASE(10) ! LkgSplit
      nTracerCntl%LkgSplitLv=0
      nTracerCntl%lDhom=FALSE
      READ(oneline, *) astring, nTracerCntl%LkgSplitLv
      IF(nLineFIeld .EQ. 2) THEN
        READ(oneline, *) astring, nTracerCntl%LkgSplitLv, nTracerCntl%lDhom
      ENDIF
    CASE(11)   !BORON
      READ(oneline, *) astring, nTracerCntl%BoronPPM
      IF(nTracerCntl%BoronPPM .GT. 0._8) nTracerCntl%lInitBoron = .TRUE.
    CASE(12)   !XENON
      nTracerCntl%lXeDyn = .TRUE.
      READ(oneline, *) astring, optfield
      CALL TOUPPER(optfield)
      IF(OptField .EQ. 'EQ') THEN
        nTracerCntl%lEqXe = .TRUE.
        nTracerCntl%lTrXe = .FALSE.
      ELSEIF(OptField .EQ. 'TR') THEN
        nTracerCntl%lEqXe = .FALSE.
        nTracerCntl%lTrXe = .TRUE.
      ENDIF
      IF(nLineField .EQ. 2) THEN
        READ(oneline, *) astring, optfield, nTracerCntl%lWriteXeEq
      ELSEIF(nLineField .EQ. 3) THEN
        READ(oneline, *) astring, optfield, nTracerCntl%lWriteXeEq
        CALL GetFn(oneline, 4, FileName(XeEqFileIdx))
        nTracerCntl%lReadXeEq = .TRUE.
      ENDIF
    CASE(13)   !NODAL
      READ(oneline, *) astring, nTracerCntl%AxSolver
      IF(nLineField .EQ. 2) THEN
        READ(oneline, *) astring, nTracerCntl%AxSolver, nSubPlane
        IF(nSubPlane .GT. 1) nTracerCntl%lSubPlane = TRUE
      ENDIF

    CASE(14)   !CRITICAL
      READ(oneline, *) astring, nTracerCntl%lCritSpec
    CASE(15)   !LINSRC
        nTracerCntl%lLinSrc=FALSE
        nTracerCntl%lLinSrcCASMO=FALSE
        IF(nLineFIeld .EQ. 1) THEN
          READ(oneline, *) astring, nTracerCntl%lLinSrc
        ELSEIF(nLineFIeld .EQ. 2) THEN
          READ(oneline, *) astring, nTracerCntl%lLinSrc, nTracerCntl%lHybrid
        ENDIF
    CASE(16)   !Fast Ray Reconstrunction
      READ(oneline, *) astring, nTracerCntl%lFastMOC
      nTracerCntl%FastMocLv = 0
      IF(nTracerCntl%lFastMOC) nTracerCntl%FastMocLv = 1
    CASE(17)   !SCAT_ORDER
    CASE(18)   !BCR_OPT
    CASE(19)   !MOCUR
       READ(oneline, *) astring, nTracerCntl%lMOCUR
       nTracerCntl%lOptUR = .TRUE.
       nTracerCntl%UserUR = 0.75
       IF(nLineFIeld .EQ. 2) THEN
         READ(oneline, *) astring, nTracerCntl%lMOCUR, nTracerCntl%lOptUR
       ELSEIF(nLineField .GE. 3) THEN
         READ(oneline, *) astring, nTracerCntl%lMOCUR, nTracerCntl%lOptUR, nTracerCntl%UserUR
       ENDIF
    CASE(20) ! AXREF_FDM
       READ(oneline, *) astring, nTracerCntl%lAxRefFDM

    CASE(21) ! MC
      nTracerCntl%TRSolver=2 ! Set Monte Carlo
      READ(oneline, *) astring, nht, ninact, nact
    CASE(22) ! MCCMFD
      READ(oneline, *) astring, lMCcmfd, lfdb
    CASE(23) ! CMFDSET
      lMCcmfdset=.TRUE.
      READ(oneline, *) astring, idx, period, skip, flag
    CASE(24) ! DCMP
      nTracerCntl%lDomainDcmp=.FALSE.
      READ(oneline, *) astring, nTracerCntl%lDomainDcmp
    CASE(25) ! NODE_MAJOR
      nTracerCntl%lNodeMajor=.FALSE. ! group major
      READ(oneline, *) astring, nTracerCntl%lNodeMajor
    CASE(26) ! AFSS - angular flux save scheme - in group major
      nTracerCntl%lAFSS=.FALSE.
      READ(oneline, *) astring, nTracerCntl%lAFSS
    CASE(27) ! NODALOPT
      READ(oneline, *) astring, nTracerCntl%lnonFuelpinFDM
      IF(nLineFIeld .EQ. 2) THEN
        READ(oneline, *) astring, nTracerCntl%lnonFuelpinFDM, nTracerCntl%lResetResErr
      ENDIF
    !--- CNJ Edit : Angular Multigrid Ray Tracing
    CASE(28) ! MULTIGRID
      READ(oneline, *) astring, nTracerCntl%lMultigrid
      IF (nLineField .EQ. 2) THEN
        READ(oneline, *) astring, nTracerCntl%lMultigrid, nTracerCntl%MultigridLV
      ENDIF
    CASE(29) ! GRIDSTR
      READ(oneline, *) astring, nTracerCntl%gridStr(1 : nLineField - 1)
      nTracerCntl%gridNum = nLineField - 1
    CASE(30) ! SUBPLN
      nTracerCntl%lSubPlane = TRUE
      READ(oneline, *) astring, maxhzfm
    CASE(31) ! POWERMODE  JSU EDIT 20190819
      IF (nLineField .EQ. 1) THEN
        READ(oneline, * ) astring, nTRACERCntl%pmode
      ELSE IF(nLineField .eq. 2) THEN
        READ(oneline, * ) astring, nTRACERCntl%pmode, nTRACERCntl%pout_mode
      END IF

  END SELECT
ENDDO
!--- Scheme case
IF (nTracerCntl%lNodeMajor) THEN
  IF (nTracerCntl%lLinSrc) THEN
    nTracerCntl%lLinSrcCASMO = TRUE
    nTracerCntl%lLinSrc = FALSE
  ENDIF
!  IF (.NOT. PE%lCUDA) THEN
!    IF (nTracerCntl%ScatOd .GT. 0 ) THEN
!      IF (.NOT. nTracerCntl%lDomainDcmp) THEN
!        nTracerCntl%lDomainDcmp = TRUE
!        IF (PE%Master) THEN
!          mesg = 'Turning DCMP scheme on : DCMP must be followed to perform Pn NM'
!          CALL message(io8, TRUE, TRUE, mesg)
!        ENDIF
!      ENDIF
!    ENDIF
!  ENDIF
ENDIF

IF (nTracerCntl%lDomainDcmp) THEN
  IF (.NOT. nTracerCntl%lNodeMajor) THEN
    nTracerCntl%lNodeMajor = TRUE
    IF (PE%Master) THEN
      mesg = 'Turning NODE MAJOR scheme on : NM must be followed to perform DCMP'
      CALL message(io8, TRUE, TRUE, mesg)
    ENDIF
  ENDIF
ELSE
  IF (nTracerCntl%lHybrid) THEN
    nTracerCntl%lDomainDcmp = TRUE
    IF (PE%Master) THEN
      mesg = 'Turning DCMP scheme on : DCMP must be followed to perform Hybrid CASMO LinSrc'
      CALL message(io8, TRUE, TRUE, mesg)
    ENDIF
    iF (.NOT. nTracerCntl%lNodeMajor) THEN
      nTracerCntl%lNodeMajor = TRUE
      IF (PE%Master) THEN
        mesg = 'Turning NODE MAJOR scheme on : NM must be followed to perform DCMP'
        CALL message(io8, TRUE, TRUE, mesg)
      ENDIF
    ENDIF
  ENDIF
ENDIF

!--- CNJ Edit : Group-wise Albedo (BYS Request)
ALLOCATE(Core%groupAlbedo(4, ng)); Core%groupAlbedo = 1.0
IF (nTracerCntl%lGroupAlbedo) THEN
  OPEN(80, FILE = TRIM(charTemp))
  DO i = 1, ng
    READ(80, *) Core%groupAlbedo(:, i)
  ENDDO
ENDIF

Backspace(indev); IF(Master) Backspace(outdev)
END SUBROUTINE

SUBROUTINE ReadDeplCard(InDev, OutDev)
USE PARAM
!USE TH_Mod,        ONLY : ThOpt
USE DeplType,       ONLY : DeplCntl_Type,   DeplLib_Type
USE DEPL_MOD,       ONLY : DeplCntl
USE DeplLib_MOD,    ONLY : DeplLib
USE FILES,          ONLY : FILENAME,        DeplFileIdx
USE IOUTIL,         ONLY : GetFn,           Terminate,      TOUPPER
USE CNTL,           ONLY : nTracerCntl
USE ALLOCS
USE PE_MOD,         ONLY : PE
IMPLICIT NONE

INTEGER :: InDev, OutDev
INTEGER           :: idcard
INTEGER,parameter :: idblock = 9
INTEGER :: IntTemp(100)
INTEGER :: i1, i2, id, i, k, nspt, ipos(100)
REAL :: REALTemp(2000)
!LOGICAL :: LDeplFile
CHARACTER(256) :: CharTemp
LOGICAL :: Master


Master = PE%master
DO while(TRUE)
  read(indev,'(a256)') oneline
  If(Master) CALL message(io8,FALSE,FALSE,oneline)
  IF(probe .EQ. BANG) cycle;     IF(probe .EQ. POUND) cycle
  IF(oneline .EQ. BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe .EQ. DOT) exit;       IF(probe.eq.SLASH) exit
  read(oneline,*) cardname; CALL toupper(cardname)
  idcard = FindCardId(idblock,cardname); IF(idcard .eq. 0) exit
  nLineField = nfields(oneline)-1
  SELECT CASE(idcard)
  !'BU_TYPE'   ,'BURNUP'    ,'OPT_DEP'    ,12*' '
    CASE(1) !BU_TYPE
      READ(ONELINE, *) ASTRING, DeplCNTL%BurnUpType
    CASE(2) ! BURNUP
      DeplCNTL%nBurnUpStep = nLineField
      READ(ONELINE, *) ASTRING, REALTemp(1:nLineField)

      DO WHILE(TRUE)
        read(indev,'(a256)') oneline
        IF(.NOT. IFnumeric(oneline)) EXIT
        If(Master) CALL message(io8,FALSE,FALSE,oneline)
        nLineField = nfields(oneline)
        i1 = DeplCNTL%nBurnUpStep + 1
        i2 = DeplCNTL%nBurnUpStep + nLineField
        DeplCNTL%nBurnUpStep = i2
        READ(ONELINE, *) REALTEMP(i1:i2)
      ENDDO
      Backspace(indev);
      CALL Dmalloc0(DeplCntl%T_efpd, 0, DeplCNTL%nBurnUpStep)
      CALL Dmalloc0(DeplCntl%T_mwdkghm, 0, DeplCNTL%nBurnUpStep)
      IF(DeplCNTL%BurnUpType == 1) THEN
        DeplCntl%T_efpd(1:DeplCNTL%nBurnUpStep)=REALTEMP(1:DeplCNTL%nBurnUpStep)
!        READ(ONELINE, *) ASTRING, DeplCntl%T_efpd(1:nLineField)
      ELSE
        DeplCntl%T_mwdkghm(1:DeplCNTL%nBurnUpStep)=REALTEMP(1:DeplCNTL%nBurnUpStep)
!        READ(ONELINE, *) ASTRING, DeplCntl%T_mwdkghm(1:nLineField)
      ENDIF
    CASE(3) !OPT_DEP

    CASE(4) !DEPLFILE
      !LDeplFile = TRUE
      !READ(ONELINE, *) ASTRING, filename(DeplFileIdx)
      !CALL getfn(oneline,2,filename(DeplFileIdx))
      !DeplLib%Filename = filename(DeplFileIdx)
      !DeplCntl%lDeplFile = .TRUE.
    CASE(5) !PC_OPT
      READ(ONELINE, *) ASTRING, DeplCntl%PC_OPT
    CASE(6)  !B10DEPL
      READ(ONELINE, *) ASTRING, DeplCntl%lB10Depl, DeplCntl%B10DeplMod, DeplCntl%vratio
      IF(.NOT. DeplCntl%lB10Depl) DeplCntl%B10DeplMod = 0
      IF(DeplCntl%vratio .LE.0.0001) DeplCntl%vratio = 0.053957
    CASE(7)  !Core_Foloow
      READ(ONELINE, *) ASTRING, DeplCntl%lCoreFollow
      IF(nLineField .GT. 1) THEN
        READ(oneline, *) ASTRING, DeplCntl%lCoreFollow, CharTemp
        !IF(CharTemp .EQ. 'STEP') DeplCntl%CoreState%LinStateChg = .FALSE.
        !IF(CharTemp .EQ. 'LIN') DeplCntl%CoreState%LinStateChg = .TRUE.
        DeplCntl%CoreState%StepStateChg = .NOT. DeplCntl%CoreState%LinStateChg
      ENDIF
    CASE(8)  !Core State
      CALL ProcCoreState(oneline)
    CASE(9)  !Flux_INIT
      READ(ONELINE, *) ASTRING, nTracerCntl%lDeplVarInit
    CASE(10)  !EFFXSDEPL
      nTracerCntl%OutpCntl%nDeplRegXsOut = nTracerCntl%OutpCntl%nDeplRegXsOut + 1
      i = nTracerCntl%OutpCntl%nDeplRegXsOut
      k = nLineField - 11
      IF(k .GT. 100) CALL TERMINATE('Max. Isotope is 100 : EFFXSDEPL')
      IF(k .LT. 1) CALL TERMINATE('Not Enough Input Argmuents : EFFXSDEPL')
      CALL FndChara(oneline, ipos, nspt, SLASH)
      READ(ONELINE, *) astring, nTracerCntl%OutpCntl%DeplRegXsOutList(1, i)
      READ(ONELINE(ipos(1)+1:), *) nTracerCntl%OutpCntl%DeplRegXsOutList(2:3, i)
      READ(ONELINE(ipos(2)+1:), *) nTracerCntl%OutpCntl%DeplRegXsOutList(4:5, i)
      READ(ONELINE(ipos(3)+1:), *) nTracerCntl%OutpCntl%DeplRegXsOutList(6:7, i)
      READ(ONELINE(ipos(4)+1:), *) nTracerCntl%OutpCntl%DeplIsoXsOutList(1:k, i)
      nTracerCntl%OutpCntl%DeplIsoXsOutList(0,i) = k
    CASE(11) !SubStep
      READ(ONELINE, *) astring,DeplCntl%nSubStep
    CASE(12) !Solver
      READ(ONELINE, *) astring,DeplCntl%Solver
    END SELECT
END DO
IF(.NOT. DeplCntl%LDeplFile) CALL Terminate('ReadInpCard : Can Not Find Depletion Chain File')
Backspace(indev); IF(Master) Backspace(outdev)

IF(nTracerCntl%lReadXeEq) nTracerCntl%lReadXeEq = .FALSE.
IF(DeplCntl%lCoreFollow) THEN
  DeplCntl%nBurnUpStep = DeplCntl%CoreState%nBurnupStep
  CALL Dmalloc0(DeplCntl%T_efpd, 0, DeplCNTL%nBurnUpStep)
  CALL Dmalloc0(DeplCntl%T_mwdkghm, 0, DeplCNTL%nBurnUpStep)
  DeplCntl%T_efpd(1:DeplCNTL%nBurnUpStep) = DeplCntl%CoreState%T_efpd(1:DeplCNTL%nBurnUpStep)
  DeplCntl%T_mwdkghm(1:DeplCNTL%nBurnUpStep) = DeplCntl%CoreState%T_mwdkghm(1:DeplCNTL%nBurnUpStep)
ENDIF

CONTAINS

SUBROUTINE ProcCoreState(Oneline0)
USE IOUTIL,       ONLY : fndchara,         terminate
IMPLICIT NONE
CHARACTER(256) :: Oneline0

REAL :: tbu, PowLv, FlowLV
REAL :: boron_keff
INTEGER :: id0
LOGICAL :: lBoronSearch
CHARACTER(50) :: astring0, boronmode

INTEGER :: ipos(100), nspt

READ(oneline0, *) astring0, id0, tbu
CALL fndchara(oneline0, ipos, nspt, SLASH)
IF(nspt .LT. 2) CALL TERMINATE('ReadInpCard : Error In CORE_STATE CARD')
! id burnup / powLv, flowLV / boron mode
DeplCntl%CoreState%nBurnupStep = max(DeplCntl%CoreState%nBurnupStep, id0)
IF(DeplCNTL%BurnUpType == 1) THEN
  DeplCntl%CoreState%T_efpd(id0) = tbu
ELSE
  DeplCntl%CoreState%T_mwdkghm(id0) = tbu
ENDIF

READ(oneline0(ipos(1)+1:ipos(2)-1), *) PowLv
FlowLV = 100.
DeplCntl%CoreState%RelPow(id0) = PowLv / 100.0_8; DeplCntl%CoreState%FlowRate(id0) = FlowLv / 100.0_8
READ(oneline0(ipos(2)+1:256), *) lBoronSearch, BoronMode
DeplCntl%CoreState%lBoronSearch(id0) = lBoronSearch
IF(lBoronSearch) THEN
  READ(BoronMode, *) boron_keff
  DeplCntl%CoreState%Target_keff(id0) = boron_keff
ELSE
  IF(boronmode .EQ. '-') THEN
    DeplCntl%CoreState%BoronPPM(id0) = -100._8
  ELSE
    IF(ifnumeric(Boronmode)) THEN
      READ(BoronMode, *) boron_keff
      DeplCntl%CoreState%BoronPPM(id0) = boron_keff
    ENDIF
  ENDIF
!  IF(ifnumeric(Boronmode)) THEN
!    READ(BoronMode, *) boron_keff
!    DeplCntl%CoreState%BoronPPM(id0) = boron_keff
!  ELSE
!    CALL TOUPPER(boronmode)
!    IF(boronmode .EQ. '-') DeplCntl%CoreState%BoronPPM(id0) = -100._8
!  ENDIF
ENDIF

END SUBROUTINE

END SUBROUTINE

SUBROUTINE ReadXsecCard(indev, outdev)
USE CNTL,          ONLY : nTracerCntl
USE DEPL_MOD,      ONLY : DeplCntl
USE DeplLib_MOD,   ONLY : DeplLib
USE PE_MOD,        ONLY : PE
USE FILES,         ONLY : FILENAME,        DeplFileIdx
IMPLICIT NONE
INTEGER           :: indev,outdev
INTEGER           :: idcard
INTEGER,PARAMETER :: idblock = 3
INTEGER :: IntTemp(100)
INTEGER :: i,k, idum1, idum2
REAL :: REALTemp(100)
CHARACTER(10) :: dumc, optfield
LOGICAL :: Master

Master = PE%Master
DO while(TRUE)
  read(indev,'(a256)') oneline
  IF(Master) CALL message(io8,FALSE,FALSE,oneline)
  IF(probe.eq.BANG) CYCLE;     IF(probe.eq.POUND) CYCLE
  IF(oneline.eq.BLANK) CYCLE;  IF(IFnumeric(oneline)) CYCLE
  IF(probe.eq.DOT) EXIT;       IF(probe.eq.SLASH) EXIT
  READ(oneline,*) cardname; CALL toupper(cardname)
  idcard = FindCardId(idblock,cardname); IF(idcard .eq. 0) EXIT
  nLineField = nfields(oneline)-1
  SELECT CASE(idcard)
    CASE(1) !LIB_TYPE

    CASE(2)  !GROUP_SPEC

    CASE(3)  !FILE

    CASE(4)  !BASE_MICRO

    CASE(5)  !DNEUT_CHI
      READ(oneline,*) astring, nTracerCntl%lchidgen, nTracerCntl%lchidkgen

    CASE(6)  !NEUT_VELO
      READ(oneline,*) astring, nTracerCntl%lfixvel
    CASE(7)  !DNP_NGRP

    CASE(8)  !DNP_BETA
      READ(oneline,*) astring, nTracerCntl%lfitbeta
    CASE(9)  !DNP_LAMBDA
      READ(oneline,*) astring, nTracerCntl%llibdcy_del, nTracerCntl%refdcy_del
    CASE(10) !DUMMY

    CASE(11) !BUCKLING
       READ(oneline, *) astring, nTracerCntl%Bsq
       IF(nTracerCntl%Bsq .GT. 0._8) nTracerCntl%lBsq = .TRUE.
    CASE(12) !RLB
    CASE(13) !SLB
    CASE(14) !GC_SPEC
      IF(nLineField .GT. 1) THEN
        IF(nLineField .NE. nTracerCntl%NGC+1) THEN
          nTracerCntl%lGCCMFD = .FALSE.
          CYCLE
        ENDIF
        !READ(oneline, *) astring, idum1, nTracerCntl%GcStructMod

        IF(nTracerCntl%GcStructMod .EQ. 1) THEN
          READ(oneline, *) astring, idum1, nTracerCntl%GcStructMod, dumc,(nTracerCntl%ELB(i), i = 1, nTracerCntl%NGC-1)
        ELSEIF(nTracerCntl%GcStructMod .EQ. 2) THEN
          READ(oneline, *) astring, idum1, nTracerCntl%GcStructMod, dumc,(nTracerCntl%GLB(i), i = 1, nTracerCntl%NGC-1)
        ELSE

        ENDIF

      ENDIF
    CASE(15) ! Deplfile
      CALL getfn(oneline,2,filename(DeplFileIdx))
      DeplLib%Filename = filename(DeplFileIdx)
      DeplCntl%lDeplFile = .TRUE.
    CASE(16) ! PHL ** photoatomic data..

    CASE(17) ! RT (stands for Resonance Treatment.)
      IF(nTracerCntl%libtyp == 2) THEN
        SELECT CASE(nLineField)
        CASE(1)
          READ(oneline, *) astring, nTracerCntl%lMLG
        CASE(2)
          READ(oneline, *) astring, nTracerCntl%lMLG , nTracerCntl%lRIF
        CASE(3)
          READ(oneline, *) astring, nTracerCntl%lMLG , nTracerCntl%lRIF, nTracerCntl%lRST
        CASE(4)
          READ(oneline, *) astring, nTracerCntl%lMLG , nTracerCntl%lRIF, nTracerCntl%lRST, nTracerCntl%nMLG
        CASE(5)
          READ(oneline, *) astring, nTracerCntl%lMLG , nTracerCntl%lRIF, nTracerCntl%lRST, nTracerCntl%nMLG, nTracerCntl%lCAT
        CASE(6)
          READ(oneline, *) astring, nTracerCntl%lMLG , nTracerCntl%lRIF, nTracerCntl%lRST, nTracerCntl%nMLG, nTracerCntl%lCAT, nTracerCntl%l4Lv
        END SELECT
      END IF
    CASE(20) ! RESO_OPT for Old XS Lib.
      IF(nTracerCntl%libtyp == 0 .or. nTracerCntl%libtyp == 3) THEN
        nTracerCntl%lrestrmt = .TRUE.
        nTracerCntl%lCAT = .TRUE.
        nTracerCntl%l4Lv = .TRUE.
        nTracerCntl%lMLG = .TRUE.
        nTracerCntl%lRIF = .FALSE.
        nTracerCntl%lRST = .FALSE.
        SELECT CASE(nLineField)
        CASE(1)
          READ(oneline, *) astring, nTracerCntl%lrestrmt
		    CASE(2)
          READ(oneline, *) astring, nTracerCntl%lrestrmt, nTracerCntl%lMLG
        CASE(3)
          READ(oneline, *) astring, nTracerCntl%lrestrmt, nTracerCntl%lMLG, nTracerCntl%l4Lv
        CASE(4)
          READ(oneline, *) astring, nTracerCntl%lrestrmt, nTracerCntl%lMLG, nTracerCntl%l4Lv, nTracerCntl%lCAT
        END SELECT
      END IF
    CASE(24) ! RIF
      READ(ONELINE, *) astring, optfield
      CALL TOUPPER(optfield)
      IF (optfield.EQ.'FXR') THEN
        nTRACERCntl%lRIFFXR = .TRUE.
        nTRACERCntl%lRIF = .TRUE.
      ELSEIF (optfield.EQ.'BON') THEN
        nTRACERCntl%lRIF = .FALSE.
      END IF
  END SELECT
ENDDO
nTRACERCntl%lRIFFXR = nTRACERCntl%lRIFFXR.AND.nTracerCntl%lrestrmt
Backspace(indev); IF(Master) Backspace(outdev)
END SUBROUTINE


SUBROUTINE ReadDcplCard(indev, outdev)
USE CNTL,          ONLY : nTracerCntl
USE GEOM,          ONLY : Core, nz
USE DcplCore_Mod,  ONLY : DcplInfo
USE PE_MOD,        ONLY : PE
IMPLICIT NONE
INTEGER           :: indev,outdev
INTEGER           :: idcard
INTEGER,PARAMETER :: idblock = 12
INTEGER :: IntTemp(100)
INTEGER :: i,k
INTEGER :: npln
REAL :: REALTemp(100)
LOGICAL :: Master

nTracerCntl%lDcpl = .TRUE.
nTracerCntl%lProblem = ldcplsseigv
master = PE%master


DO while(TRUE)
  read(indev,'(a256)') oneline
  IF(Master) CALL message(io8,FALSE,FALSE,oneline)
  IF(probe.eq.BANG) CYCLE;     IF(probe.eq.POUND) CYCLE
  IF(oneline.eq.BLANK) CYCLE;  IF(IFnumeric(oneline)) CYCLE
  IF(probe.eq.DOT) EXIT;       IF(probe.eq.SLASH) EXIT
  READ(oneline,*) cardname; CALL toupper(cardname)
  idcard = FindCardId(idblock,cardname); IF(idcard .eq. 0) EXIT
  nLineField = nfields(oneline)-1
  SELECT CASE(idcard)
    CASE(1)  !REF_PLN
      DcplInfo%nRefPln = nLineField
      CALL Dmalloc(DcplInfo%RefPln, nLineField)
      READ(oneline, *) astring, (DcplInfo%RefPln(i), i = 1, nLineField)
    CASE(2)  !REF_TEMP
      DcplInfo%lFeedBack = .TRUE.
    CASE(3)  !PLN_MAP
      nPln = nLineField
      IF(npln .NE. nz) CALL terminate('ReadInpCard : Dimmension Mismatch of PLN_MAP')
      CALL DMALLOC(DcplInfo%Pln_Map, nPln)
      READ(oneline, *) astring, (DcplInfo%Pln_Map(i), i = 1,nPln)
    CASE(4)
      READ(oneline, *) astring, nTracerCntl%lXsFtn
      IF(nLineField .GT. 1) THEN
        READ(oneline, *) astring, nTracerCntl%lXsFtn, nTracerCntl%XsFtnMod
      ENDIF
  END SELECT
ENDDO
IF(.NOT. ASSOCIATED(DcplInfo%Ref_Temp)) THEN
  CALL Dmalloc(DcplInfo%Ref_Temp, DcplInfo%nRefTemp)
   DcplInfo%Ref_Temp(1) = nTracerCntl%TempInlet
ENDIF
IF(nTracerCntl%lFeedBack)  DcplInfo%nRefTemp  = 3
IF(.NOT. nTracerCntl%lFeedBack) THEN
  DcplInfo%nRefTemp  = 3
  DEALLOCATE(DcplInfo%RefPln)
  DcplInfo%nRefPln = npln
  CALL Dmalloc(DcplInfo%RefPln, nPln)
  DO i = 1, npln
    DcplInfo%RefPln(i) = i;    DcplInfo%Pln_Map(i) = i
  ENDDO
ENDIF
Backspace(indev); IF(Master) Backspace(outdev)
END SUBROUTINE

SUBROUTINE ReadParallelCard(indev,outdev)
USE PARAM
USE CNTL,   ONLY : nTracerCntl
USE PE_MOD, ONLY : PE
IMPLICIT NONE
INTEGER           :: indev,outdev
INTEGER           :: idcard
INTEGER,parameter :: idblock = 13
INTEGER :: IntTemp(100)
INTEGER :: i,k,nb
REAL :: REALTemp(100)
LOGICAL :: Master

master = PE%master
DO while(TRUE)
  read(indev,'(a256)') oneline
  IF(master) CALL message(io8,FALSE,FALSE,oneline)
  IF(probe .EQ. BANG) cycle;     IF(probe .EQ. POUND) cycle
  IF(oneline .EQ. BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe .EQ. DOT) exit;       IF(probe.eq.SLASH) exit
  read(oneline,*) cardname; CALL toupper(cardname)
  idcard = FindCardId(idblock,cardname); IF(idcard .eq. 0) exit
  nLineField = nfields(oneline)-1
  PE%lThread = .TRUE.
  SELECT CASE(idcard)
    CASE(1)    !MOC_TRD
      READ(oneline, *) astring, PE%nThread
    CASE(2)    !NODAL_TRD
      READ(oneline, *) astring,PE%nAxThread
    CASE(3)    !CMFD_TRD
      READ(oneline, *) astring,PE%nCmfdThread
    CASE(4)    !DEPL_TRD
      READ(oneline, *) astring,PE%nDeplThread
    CASE(5)    !AX_DCP
      PE%lUsrAxDcp = .TRUE.
      PE%UsrAxDecomp(0) = nLineField
      READ(oneline, *) astring, PE%UsrAxDecomp(1:nLineField)
  END SELECT
ENDDO
Backspace(indev); IF(Master) Backspace(outdev)
END SUBROUTINE

SUBROUTINE ReadVisualCard(indev,outdev)
USE PARAM
USE CNTL,   ONLY : nTracerCntl
USE PE_MOD, ONLY : PE
IMPLICIT NONE
INTEGER           :: indev,outdev
INTEGER           :: idcard
INTEGER,parameter :: idblock = 14
INTEGER :: IntTemp(100)
INTEGER :: i,k,nb
REAL :: REALTemp(100)
LOGICAL :: Master

master = PE%master
DO while(TRUE)
  read(indev,'(a256)') oneline
  IF(master) CALL message(io8,FALSE,FALSE,oneline)
  IF(probe .EQ. BANG) cycle;     IF(probe .EQ. POUND) cycle
  IF(oneline .EQ. BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe .EQ. DOT) exit;       IF(probe.eq.SLASH) exit
  read(oneline,*) cardname; CALL toupper(cardname)
  idcard = FindCardId(idblock,cardname); IF(idcard .eq. 0) exit
  nLineField = nfields(oneline)-1
  SELECT CASE(idcard)
    CASE(1)    !VIS_MOD
      READ(oneline, *) astring, nTracerCntl%OutpCntl%VisMod
      !0 no VTK output
      !1 3-D output
      !2 2-D plane output
      !3 3-D homogenized ouptut
    CASE(2)    !FLUX_OUT
  END SELECT
ENDDO
Backspace(indev); IF(Master) Backspace(outdev)
END SUBROUTINE


SUBROUTINE ReadEditCard(indev, outdev)
USE PARAM
use geom,   ONLY : nz
USE CNTL,   ONLY : nTracerCntl
USE PE_MOD, ONLY : PE
USE Files,  ONLY : caseid,            filename,         RstFileIdx,        &
                   DeplRstFileIdx
USE IOUTIL, ONLY : Getfn
USE ALLOCS
USE FXRVAR_MOD !--- BYS edit
USE B1_Mod, ONLY : b1k, lb1k
USE XSLIB_MOD, ONLY : noghel
IMPLICIT NONE
INTEGER           :: indev,outdev
INTEGER           :: idcard
INTEGER,parameter :: idblock = 11
INTEGER :: IntTemp(100), ipos1(100), ipos2(100)
INTEGER :: i,k,nb, nspt1, nspt2
REAL :: REALTemp(100)
INTEGER :: nspt, ipos(100)
LOGICAL :: Master
CHARACTER(80) :: CharTemp(10)


master = PE%master
DO while(TRUE)
  read(indev,'(a256)') oneline
  IF(master) CALL message(io8,FALSE,FALSE,oneline)
  IF(probe .EQ. BANG) cycle;     IF(probe .EQ. POUND) cycle
  IF(oneline .EQ. BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe .EQ. DOT) exit;       IF(probe.eq.SLASH) exit
  read(oneline,*) cardname; CALL toupper(cardname)
  idcard = FindCardId(idblock,cardname); IF(idcard .eq. 0) exit
  nLineField = nfields(oneline)-1

  SELECT CASE(idcard)
    CASE(1)    !ISOTOPE
      READ(oneline, *) astring, CharTemp(1)
      CALL TOUPPER(CharTemp(1))
      IF(CharTemp(1) .EQ. 'MAJOR') THEN
        nTRACERCntl%outpCntl%IsoOutList(0) = 18
      ELSEIF(CharTemp(1) .EQ. 'ALL') THEN
        nTRACERCntl%outpCntl%IsoOutList(0) = -1
      ELSE
        nTRACERCntl%outpCntl%IsoOutList(0) = nLineField
        READ(oneline, *) astring, nTRACERCntl%OutpCntl%IsoOutList(1:nLineField)
      ENDIF
    CASE(2)    !RSTFILE
      READ(oneline, *) astring, nTRACERCntl%lWriteRst
      IF(nTRACERCntl%lWriteRst) filename(RstFileIdx) = trim(caseid)//'.RST'
      IF(nLineField .GT. 1) THEN
          !filename(RstFileIdx)
        CALL GETFN(oneline, 3, filename(RstFileIdx))
      ENDIF
    CASE(3)
      READ(oneline, *) astring, nTRACERCntl%lWriteDeplRst
      IF(nTRACERCntl%lWriteRst) filename(DeplRstFileIdx) = trim(caseid)//'.DeplRST'
      IF(nLineField .GT. 1) THEN
        CALL GETFN(oneline, 3, filename(RstFileIdx))
      ENDIF
    CASE(4) !Flux
      CALL fndchara(oneline, ipos1, nspt1, '[')
      CALL fndchara(oneline, ipos2, nspt2, ']')
      IF(nspt1 .NE. nspt2) CALL TERMINATE('MISMATCH at FLUX CARD')
      k = nTRACERCntl%OutpCntl%FluxOutList(1, 0)
      DO i = 1, nspt1
        k = k + 1;
        READ(oneline(ipos1(i)+1:ipos2(i)-1), *) IntTemp(1:4)
        nTracerCntl%outpCntl%FluxOutList(1:4, k) = IntTemp(1:4)
      ENDDO
      nTRACERCntl%outpCntl%FluxOutList(1, 0) = k
    CASE(5) !PINXS
      IF(nTracerCntl%OutpCntl%npinxsout .EQ. 0) THEN
        CALL DMALLOC(nTracerCntl%OutpCntl%PinXsOutList, nz)
        ALLOCATE(nTracerCntl%OutpCntl%PinXsFn(nz))
      ENDIF
      nTracerCntl%OutpCntl%npinxsout = nTracerCntl%OutpCntl%npinxsout + 1
      k = nTracerCntl%OutpCntl%npinxsout
      READ(ONELINE, *) astring, nTracerCntl%OutpCntl%PinXsOutList(k)
      CALL GetFn(oneline, 3, nTracerCntl%OutpCntl%PinXsFn(k))
      CONTINUE
    CASE(6) !CSPOUT
      READ(ONELINE, *) astring, nTracerCntl%OutpCntl%lCspGenOut
      IF(nTracerCntl%OutpCntl%lCspGenOut) THEN
        IF(nLineField .LT. 2) CALL TERMINATE('Not Enough Input Argument: CSPOUT')
        READ(ONELINE, *, ERR=11001) astring, nTracerCntl%OutpCntl%lCspGenOut, nTracerCntl%OutpCntl%CspGenBankId
      ENDIF
    CASE(7) !EFFXS
      BACKSPACE(indev)
      READ(indev, '(A1024)') longline
      nTracerCntl%OutpCntl%nRegXsOut = nTracerCntl%OutpCntl%nRegXsOut + 1
      i = nTracerCntl%OutpCntl%nRegXsOut
      nLineField = nfields(longline)-1
      k = nLineField - 5
      IF(k .LT. 1) CALL TERMINATE('Not Enough Input Argmuents : GEFFXS')
      CALL FndChara(longline, ipos, nspt, SLASH)
      IF (nspt.eq.2) THEN ! ASM
          IF(k .GT. 200) CALL TERMINATE('Max. Isotope is 200 : GEFFXS')
          nTracerCntl%OutpCntl%RegXsOutASM(i)=.TRUE.
          nTracerCntl%OutpCntl%RegXsOutPIN(i)=.FALSE.
          nTracerCntl%OutpCntl%RegXsOutFXR(i)=.FALSE.
          READ(longline, *, ERR=11002) astring, nTracerCntl%OutpCntl%RegXsOutList(1, i)
          READ(longline(ipos(1)+1:), *,  ERR=11002) nTracerCntl%OutpCntl%RegXsOutList(2:3, i)
          READ(longline(ipos(2)+1:), *,  ERR=11002) nTracerCntl%OutpCntl%IsoXsOutList(1:k, i)
      ELSEIF (nspt.eq.3) THEN ! PIN
          k=k-3
          IF(k .GT. 200) CALL TERMINATE('Max. Isotope is 200 : GEFFXS')
          nTracerCntl%OutpCntl%RegXsOutASM(i)=.FALSE.
          nTracerCntl%OutpCntl%RegXsOutPIN(i)=.TRUE.
          nTracerCntl%OutpCntl%RegXsOutFXR(i)=.FALSE.
          READ(longline, *, ERR=11002) astring, nTracerCntl%OutpCntl%RegXsOutList(1, i)
          READ(longline(ipos(1)+1:), *,  ERR=11002) nTracerCntl%OutpCntl%RegXsOutList(2:3, i)
          READ(longline(ipos(2)+1:), *,  ERR=11002) nTracerCntl%OutpCntl%RegXsOutList(4:5, i)
          READ(longline(ipos(3)+1:), *,  ERR=11002) nTracerCntl%OutpCntl%IsoXsOutList(1:k, i)
      ELSEIF (nspt.eq.4) THEN ! FXR
          k=k-6
          IF(k .GT. 200) CALL TERMINATE('Max. Isotope is 200 : GEFFXS')
          nTracerCntl%OutpCntl%RegXsOutASM(i)=.FALSE.
          nTracerCntl%OutpCntl%RegXsOutPIN(i)=.FALSE.
          nTracerCntl%OutpCntl%RegXsOutFXR(i)=.TRUE.
          READ(longline, *, ERR=11002) astring, nTracerCntl%OutpCntl%RegXsOutList(1, i)
          READ(longline(ipos(1)+1:), *,  ERR=11002) nTracerCntl%OutpCntl%RegXsOutList(2:3, i)
          READ(longline(ipos(2)+1:), *,  ERR=11002) nTracerCntl%OutpCntl%RegXsOutList(4:5, i)
          READ(longline(ipos(3)+1:), *,  ERR=11002) nTracerCntl%OutpCntl%RegXsOutList(6:7, i)
          READ(longline(ipos(4)+1:), *,  ERR=11002) nTracerCntl%OutpCntl%IsoXsOutList(1:k, i)
      ELSE
          CALL TERMINATE('Not Enough Input Argmuents : GEFFXS')
      ENDIF
      nTracerCntl%OutpCntl%IsoXsOutList(0,i) = k
    CASE(8)  !BINOUTP
      READ(ONELINE, *) astring, (Chartemp(i), i= 1, nLineField)
      DO i = 1, nLineField
        CALL TOUPPER(Chartemp(i))
        SELECT CASE(CharTemp(i))
          CASE('FLUX')
            nTRACERCntl%OutpCntl%lFlux_BOutp = .TRUE.
            nTRACERCntl%OutpCntl%lBoutp = nTRACERCntl%OutpCntl%lBoutp .OR. nTRACERCntl%OutpCntl%lFlux_BOutp
          CASE('DEPL')
            nTRACERCntl%OutpCntl%lDepl_BOutp = .TRUE.
            nTRACERCntl%OutpCntl%lBoutp = nTRACERCntl%OutpCntl%lDepl_BOutp .OR. nTRACERCntl%OutpCntl%lFlux_BOutp
          CASE('TH')
            nTRACERCntl%OutpCntl%lTh_BOutp = .TRUE.
            nTRACERCntl%OutpCntl%lBoutp = nTRACERCntl%OutpCntl%lTh_BOutp .OR. nTRACERCntl%OutpCntl%lFlux_BOutp
        ENDSELECT
      ENDDO

    CASE(9) !FLUX_EDIT
      READ(ONELINE, *) astring, nTRACERCntl%OutpCntl%FluxEdit
    CASE(10) !BSTEP_BOUTP

    CASE(11) !ISOTOPE_BOUTP
      READ(oneline, *) astring, CharTemp(1)
      CALL TOUPPER(CharTemp(1))
      IF(CharTemp(1) .EQ. 'MAJOR') THEN
        nTRACERCntl%outpCntl%IsoBoutp(0) = 18
      ELSEIF(CharTemp(1) .EQ. 'ALL') THEN
        nTRACERCntl%outpCntl%IsoBoutp(0) = -1
      ELSE
        nTRACERCntl%outpCntl%IsoOutList(0) = nLineField
        READ(oneline, *) astring, nTRACERCntl%OutpCntl%IsoOutList(1:nLineField)
      ENDIF
    CASE(12)  ! ---BYS edit / gc_opt GroupConst_option
      nTracerCntl%GC_spec=2  !Solution Spectrum as Default
      nTracerCntl%nGCgrp=2   !two-group GC as default
      nTracerCntl%lPinwiseGC=.FALSE.
      SELECT CASE(nLineField)
          CASE(1)
            READ(oneline, *) astring, nTracerCntl%lGC
          CASE(2)
            READ(oneline, *) astring, nTracerCntl%lGC , nTracerCntl%GC_spec
          CASE(3)
           READ(oneline, *) astring, nTracerCntl%lGC , nTracerCntl%GC_spec, nTracerCntl%lPinwiseGC
          CASE(4)
           READ(oneline, *) astring, nTracerCntl%lGC , nTracerCntl%GC_spec, nTracerCntl%lPinwiseGC, nTracerCntl%nGCgrp
          CASE(5)
           READ(oneline, *) astring, nTracerCntl%lGC , nTracerCntl%GC_spec, nTracerCntl%lPinwiseGC, nTracerCntl%nGCgrp, nTracerCntl%lGCgapHom
      ENDSELECT
    CASE(13)  ! ---BYS edit / TMOD : moderator temperature variation calc with given delta value
      READ(oneline, *) astring, tmod(1)
      ltmod=.TRUE.
    CASE(14)  ! ---BYS edit / TFUEL : fuel temperature variation calc with given delta value
      READ(oneline, *) astring, tfuel(1)
      ltfuel=.TRUE.
    CASE(15)  ! ---BYS edit / RHO : moderator density variation calc with given modified density
      READ(oneline, *) astring, rho(1)
      lrho=.TRUE.
    CASE(16)  ! ---BYS edit / CEA XS generation in FXR wise
      nTracerCntl%lFSRXS=.FALSE.
      nTracerCntl%lFSRPhi=.FALSE.
      SELECT CASE(nLineField)
          CASE(1)
            READ(oneline, *) astring, nTracerCntl%lFSRXS
          CASE(2)
            READ(oneline, *) astring, nTracerCntl%lFSRXS, nTracerCntl%lFSRPhi
      ENDSELECT

    CASE(17)  ! ---BYS edit / BKLG : buckling for k.NE.1 search option
      READ(oneline, *) astring, b1k
      lb1k=.TRUE.
    CASE(18)  ! --- BYS edit / Ray structure Gen/Read
      nTracerCntl%lRayGen=.TRUE.  ! Default : True,  if false-> read RayData
      nTracerCntl%RayGenOpt=0     ! Default : 0(no print out), if 1, printout, 2, print out then terminate
      SELECT CASE(nLineField)
          CASE(1)
            READ(oneline, *) astring, nTracerCntl%lRayGen
          CASE(2)
            READ(oneline, *) astring, nTracerCntl%lRayGen, nTracerCntl%RayGenOpt
            IF(.NOT. nTracerCntl%lRayGen) nTracerCntl%RayGenOpt=0
      ENDSELECT
    CASE(19)  ! ---BYS edit / BKLG : buckling for k.NE.1 search option
      READ(oneline, *) astring, nTracerCntl%lDetailOutput
    CASE(20)
      !READ(oneline, *) astring, nTracerCntl%lnodetime ! DELETED : KSC

    CASE(21)  !FXRMGMAC
      nTracerCntl%OutpCntl%nRegMGXsOut = nTracerCntl%OutpCntl%nRegMGXsOut + 1
      i = nTracerCntl%OutpCntl%nRegMGXsOut
      k = nLineField - 11
      IF(k .LT. 1) CALL TERMINATE('Not Enough Input Argmuents : FXRMGMAC')
      CALL FndChara(oneline, ipos, nspt, SLASH)
      READ(ONELINE, *, ERR=11003) astring, nTracerCntl%OutpCntl%RegMGXsOutList(1, i)
      READ(ONELINE(ipos(1)+1:), *,  ERR=11003) nTracerCntl%OutpCntl%RegMGXsOutList(2:3, i)
      READ(ONELINE(ipos(2)+1:), *,  ERR=11003) nTracerCntl%OutpCntl%RegMGXsOutList(4:5, i)
      READ(ONELINE(ipos(3)+1:), *,  ERR=11003) nTracerCntl%OutpCntl%RegMGXsOutList(6:7, i)
      READ(ONELINE(ipos(4)+1:), *,  ERR=11003) nTracerCntl%OutpCntl%MGBdry(1:k, i)
      nTracerCntl%OutpCntl%MGBdry(0, i)=0
      nTracerCntl%OutpCntl%MGBdry(k+1, i)=noghel
      nTracerCntl%OutpCntl%nMG(i)=k+1
    CASE(22)  ! SSPHOUT
      nTracerCntl%OutpCntl%lssphout = .TRUE.
    CASE(23)  ! --- PHS edit / 0th scattering matrix out
      nTracerCntl%OutpCntl%nRegMAT0Out = nTracerCntl%OutpCntl%nRegMAT0Out + 1
      i = nTracerCntl%OutpCntl%nRegMAT0Out
      k = nLineField - 5
      IF(k .LT. 1) CALL TERMINATE('Not Enough Input Argmuents : EFFMAT0')
      CALL FndChara(oneline, ipos, nspt, SLASH)
      IF (nspt.eq.2) THEN
          IF(k .GT. 200) CALL TERMINATE('Max. Isotope is 200 : EFFMAT0')
          nTracerCntl%OutpCntl%RegMAT0OutASM(i)=.TRUE.
          nTracerCntl%OutpCntl%RegMAT0OutPIN(i)=.FALSE.
          nTracerCntl%OutpCntl%RegMAT0OutFXR(i)=.FALSE.
          READ(ONELINE, *, ERR=11004) astring, nTracerCntl%OutpCntl%RegMAT0OutList(1, i)
          READ(ONELINE(ipos(1)+1:), *,  ERR=11004) nTracerCntl%OutpCntl%RegMAT0OutList(2:3, i)
          READ(ONELINE(ipos(2)+1:), *,  ERR=11004) nTracerCntl%OutpCntl%IsoMAT0OutList(1:k, i)
      ELSEIF (nspt.eq.3) THEN
          k=k-3
          IF(k .GT. 200) CALL TERMINATE('Max. Isotope is 200 : EFFMAT0')
          nTracerCntl%OutpCntl%RegMAT0OutASM(i)=.FALSE.
          nTracerCntl%OutpCntl%RegMAT0OutPIN(i)=.TRUE.
          nTracerCntl%OutpCntl%RegMAT0OutFXR(i)=.FALSE.
          READ(ONELINE, *, ERR=11004) astring, nTracerCntl%OutpCntl%RegMAT0OutList(1, i)
          READ(ONELINE(ipos(1)+1:), *,  ERR=11004) nTracerCntl%OutpCntl%RegMAT0OutList(2:3, i)
          READ(ONELINE(ipos(2)+1:), *,  ERR=11004) nTracerCntl%OutpCntl%RegMAT0OutList(4:5, i)
          READ(ONELINE(ipos(3)+1:), *,  ERR=11004) nTracerCntl%OutpCntl%IsoMAT0OutList(1:k, i)
      ELSEIF (nspt.eq.4) THEN
          k=k-6
          IF(k .GT. 200) CALL TERMINATE('Max. Isotope is 200 : EFFMAT0')
          nTracerCntl%OutpCntl%RegMAT0OutASM(i)=.FALSE.
          nTracerCntl%OutpCntl%RegMAT0OutPIN(i)=.FALSE.
          nTracerCntl%OutpCntl%RegMAT0OutFXR(i)=.TRUE.
          READ(ONELINE, *, ERR=11004) astring, nTracerCntl%OutpCntl%RegMAT0OutList(1, i)
          READ(ONELINE(ipos(1)+1:), *,  ERR=11004) nTracerCntl%OutpCntl%RegMAT0OutList(2:3, i)
          READ(ONELINE(ipos(2)+1:), *,  ERR=11004) nTracerCntl%OutpCntl%RegMAT0OutList(4:5, i)
          READ(ONELINE(ipos(3)+1:), *,  ERR=11004) nTracerCntl%OutpCntl%RegMAT0OutList(6:7, i)
          READ(ONELINE(ipos(4)+1:), *,  ERR=11004) nTracerCntl%OutpCntl%IsoMAT0OutList(1:k, i)
      ELSE
          CALL TERMINATE('Not Enough Input Argmuents : EFFMAT0')
      ENDIF
      nTracerCntl%OutpCntl%IsoMAT0OutList(0,i) = k
    CASE(24) !PHIM
      nTracerCntl%OutpCntl%nRegPhimOut = nTracerCntl%OutpCntl%nRegPhimOut + 1
      i = nTracerCntl%OutpCntl%nRegPhimOut
      k = nLineField - 5
      IF(k .LT. 1) CALL TERMINATE('Not Enough Input Argmuents : EFFXS')
      CALL FndChara(oneline, ipos, nspt, SLASH)
      IF (nspt.eq.2) THEN
          nTracerCntl%OutpCntl%RegPhimOutASM(i)=.TRUE.
          nTracerCntl%OutpCntl%RegPhimOutPIN(i)=.FALSE.
          nTracerCntl%OutpCntl%RegPhimOutFXR(i)=.FALSE.
          READ(ONELINE, *, ERR=11005) astring, nTracerCntl%OutpCntl%RegPhimOutList(1, i)
          READ(ONELINE(ipos(1)+1:), *,  ERR=11005) nTracerCntl%OutpCntl%RegPhimOutList(2:3, i)
          READ(ONELINE(ipos(2)+1:), *,  ERR=11005) nTracerCntl%OutpCntl%PhimOrdOutList(1:k, i)
      ELSEIF (nspt.eq.3) THEN
          k=k-3
          nTracerCntl%OutpCntl%RegPhimOutASM(i)=.FALSE.
          nTracerCntl%OutpCntl%RegPhimOutPIN(i)=.TRUE.
          nTracerCntl%OutpCntl%RegPhimOutFXR(i)=.FALSE.
          READ(ONELINE, *, ERR=11005) astring, nTracerCntl%OutpCntl%RegPhimOutList(1, i)
          READ(ONELINE(ipos(1)+1:), *,  ERR=11005) nTracerCntl%OutpCntl%RegPhimOutList(2:3, i)
          READ(ONELINE(ipos(2)+1:), *,  ERR=11005) nTracerCntl%OutpCntl%RegPhimOutList(4:5, i)
          READ(ONELINE(ipos(3)+1:), *,  ERR=11005) nTracerCntl%OutpCntl%PhimOrdOutList(1:k, i)
      ELSEIF (nspt.eq.4) THEN
          k=k-6
          nTracerCntl%OutpCntl%RegPhimOutASM(i)=.FALSE.
          nTracerCntl%OutpCntl%RegPhimOutPIN(i)=.FALSE.
          nTracerCntl%OutpCntl%RegPhimOutFXR(i)=.TRUE.
          READ(ONELINE, *, ERR=11005) astring, nTracerCntl%OutpCntl%RegPhimOutList(1, i)
          READ(ONELINE(ipos(1)+1:), *,  ERR=11005) nTracerCntl%OutpCntl%RegPhimOutList(2:3, i)
          READ(ONELINE(ipos(2)+1:), *,  ERR=11005) nTracerCntl%OutpCntl%RegPhimOutList(4:5, i)
          READ(ONELINE(ipos(3)+1:), *,  ERR=11005) nTracerCntl%OutpCntl%RegPhimOutList(6:7, i)
          READ(ONELINE(ipos(4)+1:), *,  ERR=11005) nTracerCntl%OutpCntl%PhimOrdOutList(1:k, i)
      ELSE
          CALL TERMINATE('Not Enough Input Argmuents : PHIM')
      ENDIF
      nTracerCntl%OutpCntl%PhimOrdOutList(0,i) = k
    CASE(25) !GEFFXS
      nTracerCntl%OutpCntl%nRegPhXsOut = nTracerCntl%OutpCntl%nRegPhXsOut + 1
      i = nTracerCntl%OutpCntl%nRegPhXsOut
      k = nLineField - 5
      IF(k .LT. 1) CALL TERMINATE('Not Enough Input Argmuents : GEFFXS')
      CALL FndChara(oneline, ipos, nspt, SLASH)
      IF (nspt.eq.2) THEN
          IF(k .GT. 200) CALL TERMINATE('Max. Isotope is 200 : GEFFXS')
          nTracerCntl%OutpCntl%RegPhXsOutASM(i)=.TRUE.
          nTracerCntl%OutpCntl%RegPhXsOutPIN(i)=.FALSE.
          nTracerCntl%OutpCntl%RegPhXsOutFXR(i)=.FALSE.
          READ(ONELINE, *, ERR=11006) astring, nTracerCntl%OutpCntl%RegPhXsOutList(1, i)
          READ(ONELINE(ipos(1)+1:), *,  ERR=11006) nTracerCntl%OutpCntl%RegPhXsOutList(2:3, i)
          READ(ONELINE(ipos(2)+1:), *,  ERR=11006) nTracerCntl%OutpCntl%IsoPhXsOutList(1:k, i)
      ELSEIF (nspt.eq.3) THEN
          k=k-3
          IF(k .GT. 200) CALL TERMINATE('Max. Isotope is 200 : GEFFXS')
          nTracerCntl%OutpCntl%RegPhXsOutASM(i)=.FALSE.
          nTracerCntl%OutpCntl%RegPhXsOutPIN(i)=.TRUE.
          nTracerCntl%OutpCntl%RegPhXsOutFXR(i)=.FALSE.
          READ(ONELINE, *, ERR=11006) astring, nTracerCntl%OutpCntl%RegPhXsOutList(1, i)
          READ(ONELINE(ipos(1)+1:), *,  ERR=11006) nTracerCntl%OutpCntl%RegPhXsOutList(2:3, i)
          READ(ONELINE(ipos(2)+1:), *,  ERR=11006) nTracerCntl%OutpCntl%RegPhXsOutList(4:5, i)
          READ(ONELINE(ipos(3)+1:), *,  ERR=11006) nTracerCntl%OutpCntl%IsoPhXsOutList(1:k, i)
      ELSEIF (nspt.eq.4) THEN
          k=k-6
          IF(k .GT. 200) CALL TERMINATE('Max. Isotope is 200 : GEFFXS')
          nTracerCntl%OutpCntl%RegPhXsOutASM(i)=.FALSE.
          nTracerCntl%OutpCntl%RegPhXsOutPIN(i)=.FALSE.
          nTracerCntl%OutpCntl%RegPhXsOutFXR(i)=.TRUE.
      READ(ONELINE, *, ERR=11006) astring, nTracerCntl%OutpCntl%RegPhXsOutList(1, i)
      READ(ONELINE(ipos(1)+1:), *,  ERR=11006) nTracerCntl%OutpCntl%RegPhXsOutList(2:3, i)
      READ(ONELINE(ipos(2)+1:), *,  ERR=11006) nTracerCntl%OutpCntl%RegPhXsOutList(4:5, i)
      READ(ONELINE(ipos(3)+1:), *,  ERR=11006) nTracerCntl%OutpCntl%RegPhXsOutList(6:7, i)
      READ(ONELINE(ipos(4)+1:), *,  ERR=11006) nTracerCntl%OutpCntl%IsoPhXsOutList(1:k, i)
      ELSE
          CALL TERMINATE('Not Enough Input Argmuents : GEFFXS')
      ENDIF
      nTracerCntl%OutpCntl%IsoPhXsOutList(0,i) = k
    CASE(26) ! KERMA
      nTracerCntl%OutpCntl%nRegKERMAOut = nTracerCntl%OutpCntl%nRegKERMAOut + 1
      i = nTracerCntl%OutpCntl%nRegKERMAOut
      k = nLineField - 5
      IF(k .LT. 1) CALL TERMINATE('Not Enough Input Argmuents : KERMA')
      CALL FndChara(oneline, ipos, nspt, SLASH)
      IF (nspt.eq.2) THEN
          IF(k .GT. 200) CALL TERMINATE('Max. Isotope is 200 : KERMA')
          nTracerCntl%OutpCntl%RegKERMAOutASM(i)=.TRUE.
          nTracerCntl%OutpCntl%RegKERMAOutPIN(i)=.FALSE.
          nTracerCntl%OutpCntl%RegKERMAOutFXR(i)=.FALSE.
          READ(ONELINE, *, ERR=11007) astring, nTracerCntl%OutpCntl%RegKERMAOutList(1, i)
          READ(ONELINE(ipos(1)+1:), *,  ERR=11007) nTracerCntl%OutpCntl%RegKERMAOutList(2:3, i)
          READ(ONELINE(ipos(2)+1:), *,  ERR=11007) nTracerCntl%OutpCntl%IsoKERMAOutList(1:k, i)
      ELSEIF (nspt.eq.3) THEN
          k=k-3
          IF(k .GT. 200) CALL TERMINATE('Max. Isotope is 200 : KERMA')
          nTracerCntl%OutpCntl%RegKERMAOutASM(i)=.FALSE.
          nTracerCntl%OutpCntl%RegKERMAOutPIN(i)=.TRUE.
          nTracerCntl%OutpCntl%RegKERMAOutFXR(i)=.FALSE.
          READ(ONELINE, *, ERR=11007) astring, nTracerCntl%OutpCntl%RegKERMAOutList(1, i)
          READ(ONELINE(ipos(1)+1:), *,  ERR=11007) nTracerCntl%OutpCntl%RegKERMAOutList(2:3, i)
          READ(ONELINE(ipos(2)+1:), *,  ERR=11007) nTracerCntl%OutpCntl%RegKERMAOutList(4:5, i)
          READ(ONELINE(ipos(3)+1:), *,  ERR=11007) nTracerCntl%OutpCntl%IsoKERMAOutList(1:k, i)
      ELSEIF (nspt.eq.4) THEN
          k=k-6
          IF(k .GT. 200) CALL TERMINATE('Max. Isotope is 200 : KERMA')
          nTracerCntl%OutpCntl%RegKERMAOutASM(i)=.FALSE.
          nTracerCntl%OutpCntl%RegKERMAOutPIN(i)=.FALSE.
          nTracerCntl%OutpCntl%RegKERMAOutFXR(i)=.TRUE.
      READ(ONELINE, *, ERR=11007) astring, nTracerCntl%OutpCntl%RegKERMAOutList(1, i)
      READ(ONELINE(ipos(1)+1:), *,  ERR=11007) nTracerCntl%OutpCntl%RegKERMAOutList(2:3, i)
      READ(ONELINE(ipos(2)+1:), *,  ERR=11007) nTracerCntl%OutpCntl%RegKERMAOutList(4:5, i)
      READ(ONELINE(ipos(3)+1:), *,  ERR=11007) nTracerCntl%OutpCntl%RegKERMAOutList(6:7, i)
      READ(ONELINE(ipos(4)+1:), *,  ERR=11007) nTracerCntl%OutpCntl%IsoKERMAOutList(1:k, i)
      ELSE
          CALL TERMINATE('Not Enough Input Argmuents : KERMA')
      ENDIF
      nTracerCntl%OutpCntl%IsoKERMAOutList(0,i) = k
  END SELECT
ENDDO
Backspace(indev); IF(Master) Backspace(outdev)
RETURN
11001 CALL TERMINATE('Input Argument Error : CSPOUT')
11002 CALL TERMINATE('Input Argument Error : EFFXS')
11003 CALL TERMINATE('Input Argument Error : FXRMGMAC')
11004 CALL TERMINATE('Input Argument Error : EFFMAT0')
11005 CALL TERMINATE('Input Argument Error : PHIM')
11006 CALL TERMINATE('Input Argument Error : GEFFXS')
11007 CALL TERMINATE('Input Argument Error : KERMA')
END SUBROUTINE

SUBROUTINE ReadLpShfCard(indev, outdev)
USE PARAM
USE PE_MOD, ONLY : PE
USE CNTL,   ONLY : nTracerCntl
USE GEOM,   ONLY : Core
USE LpShf_mod, ONLY : LpShf
IMPLICIT NONE
INTEGER           :: indev,outdev
INTEGER           :: idcard
INTEGER,parameter :: idblock = 15
INTEGER :: IntTemp(100)
INTEGER :: i,k,nb
REAL :: REALTemp(100)
LOGICAL :: Master

INTEGER :: ixya, icyc

master = PE%master
nTRACERCntl%LpShf = .TRUE.
ALLOCATE(LpShf%SHF(Core%nxya))
DO while(TRUE)
  read(indev,'(a256)') oneline
  IF(master) CALL message(io8,FALSE,FALSE,oneline)
  IF(probe .EQ. BANG) cycle;     IF(probe .EQ. POUND) cycle
  IF(oneline .EQ. BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe .EQ. DOT) exit;       IF(probe.eq.SLASH) exit
  read(oneline,*) cardname; CALL toupper(cardname)
  idcard = FindCardId(idblock,cardname); IF(idcard .eq. 0) exit
  nLineField = nfields(oneline)-1
  SELECT CASE(idcard)
    CASE(1) !CYCLE ID
      READ(oneline , *) astring, LpShf%cycleid
    CASE(2) !RSTLIB     !Restart file
      CALL ReadRstLib(indev, outdev)
    CASE(3) !PUL
      CALL ReadPUL(indev, outdev)
    CASE(4) !SHF
      CALL ReadShfInfo(indev, outdev)
    CASE(5) !CYCLE
      CALL TERMINATE('Improper location of CYCLE CARD')
    CASE(6)
      READ(oneline, *) astring, nTracerCntl%lRstCal
      IF(nTracerCntl%lRstCal) THEN
        READ(oneline, *) astring, nTracerCntl%lRstCal, LpShf%RstCycleId
      ENDIF
    CASE(7) !
      !CALL TERMINATE('Wrong input line staring with [ ')
    CASE(8) !PLN_MAP
      READ(oneline, *) astring, (LpShf%Pln_MAP(i), i = 1, Core%nz)
      LpShf%lPln_Map = .TRUE.
    CASE(9) !RMV_BP
      CALL ReadRmvBP(oneline)
    END SELECT
ENDDO
BACKSPACE(indev); IF(Master) BACKSPACE(outdev)
!LpShf Check
IF(nTracerCntl%lRstCal) RETURN
DO ixya = 1, Core%nxya
  IF(LpShf%Shf(ixya)%lFreshFuel) CYCLE
  IF(LpShf%Shf(ixya)%lNoneFuel) CYCLE
  icyc = LpShf%Shf(ixya)%cycleid
  IF(.NOT. LpShf%lDataExist(icyc)) THEN
    CALL TERMINATE('Required Restart File is not specified')
  ENDIF
ENDDO
END SUBROUTINE


SUBROUTINE ReadRmvBP(datline)
USE PARAM
USE CNTL,      ONLY : nTracerCntl
USE LpShf_mod, ONLY : RmvBP, nRmvBP
IMPLICIT NONE
CHARACTER(256) :: datline
INTEGER :: ipos(100), nspt
INTEGER :: i
CALL fndchara(datline, ipos, nspt, SLASH)
DO i = 1, nspt-1
  READ(datline(ipos(i)+1:256), *) RmvBP(1:2, i)
  nRmvBP = nRmvBP + 1
ENDDO
END SUBROUTINE

SUBROUTINE ReadShfInfo(indev, outdev)
USE PARAM
USE PE_MOD, ONLY : PE
USE CNTL,   ONLY : nTracerCntl
USE LpShf_mod, ONLY : LpShf,                                &
                      ShfInfo,    Xidx2No
USE geom,    ONLY : CORE
IMPLICIT NONE
INTEGER           :: indev,outdev
INTEGER           :: idcard
INTEGER,parameter :: idblock = 15
INTEGER :: ipos1(100), nspt1, ipos2(100), nspt2
INTEGER :: id
INTEGER :: i, j, k
LOGICAL :: Master

CHARACTER(80) :: field
CHARACTER(10) :: dat1
INTEGER :: ixya

master = PE%master

ixya = 0

nLineField = nfields(oneline)-1
IF(nLineField .GT. 1) THEN
  READ(oneline, *) astring, dat1(1:1), LpShf%IdxSt(2)
  CALL TOUPPER(dat1(1:1))
  LpShf%IdxSt(1) = Xidx2No(dat1(1:1))
ENDIF


nTRACERCntl%LpShf = .TRUE.
DO while(TRUE)
  read(indev,'(a256)') oneline
  IF(master) CALL message(io8,FALSE,FALSE,oneline)
  IF(probe .EQ. BANG) cycle;     IF(probe .EQ. POUND) cycle
  IF(oneline .EQ. BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe .EQ. DOT) exit;       IF(probe.eq.SLASH) exit
  read(oneline,*) cardname; CALL toupper(cardname)
  idcard = FindCardId(idblock,cardname);
  IF(idcard .NE. 7) exit
  nLineField = nfields(oneline)-1
  CALL fndchara(oneline, ipos1, nspt1, '[')
  CALL fndchara(oneline, ipos2, nspt2, ']')
  IF(nspt1 .NE. nspt2) CALL TERMINATE('Wrong Shf Card Input')
  DO j = 1, nspt1
    ixya = ixya + 1
    field = ''
    field = oneline(ipos1(j)+1:ipos2(j)-1)
    READ(field, *) dat1
    CALL TOUPPER(dat1)
    IF(dat1 .EQ.'F') THEN !Fresh Fuel
      LpShf%Shf(ixya)%lFreshFuel = .TRUE.
      LpShf%Shf(ixya)%cycleid = 0
    ELSEIF(dat1 .EQ. 'N') THEN ! Non Fuel
      LpShf%Shf(ixya)%lFreshFuel = .FALSE.
      LpShf%Shf(ixya)%LNoneFuel = .TRUE.
      LpShf%Shf(ixya)%cycleid = 0
    ELSE
      CALL ShfInfo(field, LpShf%Shf(ixya), LpShf%IdxSt)
    ENDIF
  ENDDO
  CONTINUE
ENDDO
BACKSPACE(indev); IF(Master) BACKSPACE(outdev)
END SUBROUTINE


SUBROUTINE ReadRstLib(indev, outdev)
USE PARAM
USE PE_MOD, ONLY : PE
USE CNTL,   ONLY : nTracerCntl
USE LpShf_mod, ONLY : LpShf
IMPLICIT NONE
INTEGER           :: indev,outdev
INTEGER           :: idcard
INTEGER,parameter :: idblock = 15
INTEGER :: id
LOGICAL :: Master

master = PE%master
DO while(TRUE)
  read(indev,'(a256)') oneline
  IF(master) CALL message(io8,FALSE,FALSE,oneline)
  IF(probe .EQ. BANG) cycle;     IF(probe .EQ. POUND) cycle
  IF(oneline .EQ. BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe .EQ. DOT) exit;       IF(probe.eq.SLASH) exit
  read(oneline,*) cardname; CALL toupper(cardname)
  idcard = FindCardId(idblock,cardname); IF(idcard .NE. 5) exit
  nLineField = nfields(oneline)-1
  read(oneline, *) astring, id
  CALL GETFN(oneline, 3, LpShf%Rstfiles(id))
  LpSHf%lDataExist(id) = .TRUE.
  CONTINUE
  LpShf%nrstfile = max(LpShf%nrstfile, id)
ENDDO
Backspace(indev); IF(Master) Backspace(outdev)
END SUBROUTINE

SUBROUTINE ReadPUL(indev, outdev)
USE PARAM
USE PE_MOD, ONLY : PE
USE CNTL,   ONLY : nTracerCntl
USE LpShf_mod, ONLY : LpShf
IMPLICIT NONE
INTEGER           :: indev,outdev
INTEGER           :: idcard
INTEGER,parameter :: idblock = 15
INTEGER :: id, idum
LOGICAL :: Master

master = PE%master
LpShf%lPUL = .TRUE.
DO while(TRUE)
  read(indev,'(a256)') oneline
  IF(master) CALL message(io8,FALSE,FALSE,oneline)
  IF(probe .EQ. BANG) cycle;     IF(probe .EQ. POUND) cycle
  IF(oneline .EQ. BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe .EQ. DOT) exit;       IF(probe.eq.SLASH) exit
  read(oneline,*) cardname; CALL toupper(cardname)
  idcard = FindCardId(idblock,cardname); IF(idcard .NE. 5) exit
  nLineField = nfields(oneline)-1
  read(oneline, *) astring, id
  read(oneline, *) astring, idum, LpShf%PUL(id)
  CONTINUE
ENDDO
Backspace(indev); IF(Master) Backspace(outdev)
END SUBROUTINE

SUBROUTINE ReadXsChange(Dataline, XSChange, CntlCusping)
USE PARAM
USE TYPEDEF,           ONLY : XsChange_Type
USE IOUTIL,            ONLY : FndChara
USE GEOM,              ONLY : Core
USE ALLOCS
IMPLICIT NONE
CHARACTER(512) :: Dataline
TYPE(XsChange_Type) :: XsChange
LOGICAL :: CntlCusping

CHARACTER(256) :: Dataline2
INTEGER :: ipos(100), ipos2(100)
INTEGER :: nspt, nspt2
INTEGER :: ixya
INTEGER :: nxa, nya, i1, i2
INTEGER :: i
LOGICAL :: lpln, lasy, lcusping

READ(oneline, *) astring, XsChange%iso0, XsChange%iso1, XsChange%tbeg, XsChange%tend
IF(XsChange%tbeg .EQ. XsChange%tend) XsChange%lStepFunc = .TRUE.
CALL fndchara(dataline,ipos,nspt,SLASH)

lpln = .FALSE.
lasy = .FALSE.
lcusping = .FALSE.
IF(nspt .GT. 1) THEN
  DO i = ipos(1)+1, ipos(2)-1
    IF(dataline(i:i) .NE. '') THEN
      lpln = .TRUE.
      EXIT
    END IF
  END DO
  IF(nspt .GT. 2) THEN
    DO i = ipos(2)+1, ipos(3)-1
      IF(dataline(i:i) .NE. '') THEN
        lasy = .TRUE.
        EXIT
      END IF
    END DO
    IF(nspt .GT. 3) THEN
      DO i = ipos(3)+1, ipos(4)-1
        IF(dataline(i:i) .NE. '') THEN
          lcusping = .TRUE.
          EXIT
        END IF
      END DO
    END IF
  END IF
END IF
IF(lpln .OR. lasy) THEN
  XsChange%lSupl = .TRUE.
END IF

IF(lpln) THEN
  Dataline2 = ''; Dataline2 = Dataline(ipos(1)+1:ipos(2)-1)
  XsChange%field1 = Dataline2
  READ(Dataline2, *) i1, i2
  XsChange%izbeg = i1; XsChange%izend = i2
ELSE
  XsChange%izbeg = 1;   XsChange%izend = Core%nz
END IF

IF(lasy) THEN
  Dataline2 = ''; Dataline2 = Dataline(ipos(2)+1:ipos(3)-1)
  XsChange%field2 = Dataline2

  CALL FndChara(XsChange%field2, ipos2, nspt2, '[')

  XsChange%nasy = nspt2
  CALL DMALLOC(XsChange%AsyList, nspt2)
  DO ixya = 1, nspt2
    READ(XsChange%field2(ipos2(ixya)+1:256), *) i1, i2
    XsChange%AsyList(ixya) = Core%CoreIdx(i1, i2)
  ENDDO
ELSE
  XsChange%nasy = Core%nxya
  CALL DMALLOC(XsChange%AsyList, Core%nxya)
  DO ixya = 1, Core%nxya
    XsChange%AsyList(ixya) = ixya
  ENDDO
END IF

IF(lcusping) THEN
  Dataline2 = ''; Dataline2 = Dataline(ipos(3)+1:ipos(4)-1)
  READ(Dataline2, *) XsChange%lCusping, XsChange%lCuspingDirection
ELSE
  XsChange%lCusping = .FALSE.
END IF
IF(XsChange%lCusping) CntlCusping = .TRUE.

END SUBROUTINE

SUBROUTINE ReadXsCntlRod(Dataline, XSCntlRod, CntlCusping)
USE PARAM
USE TYPEDEF,           ONLY : XsCntlRod_Type
USE IOUTIL,            ONLY : FndChara
USE GEOM,              ONLY : Core
USE ALLOCS
IMPLICIT NONE
CHARACTER(512) :: Dataline
TYPE(XsCntlRod_Type) :: XsCntlRod
LOGICAL :: CntlCusping

CHARACTER(256) :: Dataline2
INTEGER :: ipos(100), ipos2(100)
INTEGER :: nspt, nspt2
INTEGER :: ixya
INTEGER :: nxa, nya, i1, i2
INTEGER :: i
LOGICAL :: lpln, lasy, lcusping

READ(oneline, *) astring, XsCntlRod%iso0, XsCntlRod%iso1, XsCntlRod%wt
CALL fndchara(dataline,ipos,nspt,SLASH)

lpln = .FALSE.
lasy = .FALSE.
lcusping = .FALSE.
IF(nspt .GT. 1) THEN
  DO i = ipos(1)+1, ipos(2)-1
    IF(dataline(i:i) .NE. '') THEN
      lpln = .TRUE.
      EXIT
    END IF
  END DO
  IF(nspt .GT. 2) THEN
    DO i = ipos(2)+1, ipos(3)-1
      IF(dataline(i:i) .NE. '') THEN
        lasy = .TRUE.
        EXIT
      END IF
    END DO
    IF(nspt .GT. 3) THEN
      DO i = ipos(3)+1, ipos(4)-1
        IF(dataline(i:i) .NE. '') THEN
          lcusping = .TRUE.
          EXIT
        END IF
      END DO
    END IF
  END IF
END IF

IF(lpln) THEN
  Dataline2 = ''; Dataline2 = Dataline(ipos(1)+1:ipos(2)-1)
  XsCntlRod%field1 = Dataline2
  READ(Dataline2, *) i1, i2
  XsCntlRod%izbeg = i1; XsCntlRod%izend = i2
ELSE
  XsCntlRod%izbeg = 1;   XsCntlRod%izend = Core%nz
END IF

IF(lasy) THEN
  Dataline2 = ''; Dataline2 = Dataline(ipos(2)+1:ipos(3)-1)
  XsCntlRod%field2 = Dataline2

  CALL FndChara(XsCntlRod%field2, ipos2, nspt2, '[')

  XsCntlRod%nasy = nspt2
  CALL DMALLOC(XsCntlRod%AsyList, nspt2)
  DO ixya = 1, nspt2
    READ(XsCntlRod%field2(ipos2(ixya)+1:256), *) i1, i2
    XsCntlRod%AsyList(ixya) = Core%CoreIdx(i1, i2)
  ENDDO
ELSE
  XsCntlRod%nasy = Core%nxya
  CALL DMALLOC(XsCntlRod%AsyList, Core%nxya)
  DO ixya = 1, Core%nxya
    XsCntlRod%AsyList(ixya) = ixya
  ENDDO
END IF

IF(lcusping) THEN
  Dataline2 = ''; Dataline2 = Dataline(ipos(3)+1:ipos(4)-1)
  READ(Dataline2, *) XsCntlRod%lCusping, XsCntlRod%lCuspingDirection
ELSE
  XsCntlRod%lCusping = .FALSE.
END IF
IF(XsCntlRod%lCusping) CntlCusping = .TRUE.

END SUBROUTINE

SUBROUTINE ReadXsNoise(Dataline, XsNoise)
USE PARAM
USE TYPEDEF,        ONLY : XsNoise_Type
USE ioutil,         ONLY : FndChara
USE geom,           ONLY : Core
IMPLICIT NONE
CHARACTER(512) :: Dataline
TYPE(XsNoise_Type) :: XsNoise


CHARACTER(512) :: Dataline2
CHARACTER(1) :: xstype
INTEGER :: ipos(100)
INTEGER :: nspt, nLinefield

READ(oneline, *) astring, XsNoise%iso0, xstype
CALL toupper(xstype)
SELECT CASE(xstype)
CASE('C') ! Capture
  XsNoise%itype = 1
CASE('F') ! Fission
  XsNoise%itype = 2
CASE('S') ! Scattering
  XsNoise%itype = 3
END SELECT
CALL fndChara(dataline, ipos, nspt, SLASH)

dataline2 = ' '
dataline2 = dataline(ipos(1)+1:ipos(2)-1)
READ(dataline2, *) XsNoise%amp, XsNoise%freq, XsNoise%phase

dataline2 = ' '
dataline2 = dataline(ipos(2)+1:ipos(3)-1)
READ(dataline2, *) XsNoise%ixa, XsNoise%iya
XsNoise%ixya = Core%CoreIdx(XsNoise%ixa, XsNoise%iya)

dataline2 = ' '
dataline2 = dataline(ipos(3)+1:ipos(4)-1)
READ(dataline2, *) XsNoise%ix, XsNoise%iy

dataline2 = ' '
dataline2 = dataline(ipos(4)+1:ipos(5)-1)
READ(dataline2, *) XsNoise%izbeg, XsNoise%izend

END SUBROUTINE

SUBROUTINE ReadCntlRodCard(indev,outdev)
USE PARAM
use geom,        ONLY : lgap,                    nCellX0,                Core
USE CNTL,        ONLY : nTracerCntl
USE PE_MOD,      ONLY : PE
USE CNTLROD_MOD, ONLY : GetInp_CrCellMap,       GetInp_CrAsyConf,         GetInp_CrBank,     &
                        GetInp_CrBank,          GetInp_CrPosition,        GetInp_CrPosChg,   &
                        GetInp_CrMvDom,         GetInp_CrDecusp
USE CrCsp_Mod,   ONLY : GetInp_CSPFILE,         GetInp_CSPMAP
IMPLICIT NONE
INTEGER           :: indev,outdev
INTEGER           :: idcard
INTEGER,parameter :: idblock = 16
INTEGER :: IntTemp(100)
INTEGER :: i,k,nb
CHARACTER(10) :: BankName
LOGICAL :: Master

master = PE%master
DO while(TRUE)
  read(indev,'(a256)') oneline
  IF(master) CALL message(io8,FALSE,FALSE,oneline)
  IF(probe .EQ. BANG) cycle;     IF(probe .EQ. POUND) cycle
  IF(oneline .EQ. BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe .EQ. DOT) exit;       IF(probe.eq.SLASH) exit
  read(oneline,*) cardname; CALL toupper(cardname)
  idcard = FindCardId(idblock,cardname); IF(idcard .eq. 0) exit
  nLineField = nfields(oneline)-1
  SELECT CASE(idcard)
    CASE(1)    !CR_CELL
      READ(oneline, *) astring, IntTemp(1:2)
      CALL GetInp_CrCellMap(IntTemp(1), IntTemp(2))
    CASE(2)    !CR_ASYCONF
      READ(oneline, *) astring, IntTemp(1:2)
      CALL GetInp_CrAsyConf(IntTemp(1), IntTemp(2), nCellX0, lgap, indev, outdev, PE)
    CASE(3)    !CR_BANK
      READ(oneline, *) astring, IntTemp(1), BankName
      CALL GetInp_CrBank(IntTemp(1), BankName, Core%nxya, Core%nxa, Core%nya, indev, outdev, PE)
    CASE(4)    !CR_POS
      CALL GetInp_CrPosition(Oneline, nLineField, PE)
    CASE(5)    !CR_POSCHG
      CALL GetInp_CrPosChg(Oneline, nLineField, PE)
    CASE(6)    !CSP_FILE
      CALL GetInp_CSPFILE(Oneline)
    CASE(7)    !CSP_MAP
      CALL GetInp_CSPMAP(Core%nxya, Core%nxa, Core%nya, Indev, outdev, PE)
    CASE(8)    !CRMV_DOM  Control Rod Moving Domain
      CALL GetInp_CrMvDom(OneLine, PE)
    CASE(10)   !CR_DECUSP
      nTracerCntl%lDecusp = .TRUE.
      CALL GetInp_CrDeCusp(Oneline)
  END SELECT
ENDDO
Backspace(indev); IF(Master) Backspace(outdev)
nTracerCntl%lCrInfo = .TRUE.

END SUBROUTINE

SUBROUTINE ReadCuspingCard(indev,outdev)
USE PARAM
USE GEOM,           ONLY : nz
!use geom,           ONLY : lgap,                    nCellX0,                Core
USE CNTL,           ONLY : nTracerCntl
USE PE_MOD,         ONLY : PE
USE CrCspGen_MOD,   ONLY : CrCspGenCntl
USE ioutil,         ONLY : Getfn,        TERMINATE
USE ALLOCS
!USE CNTLROD_MOD, ONLY : GetInp_CrCellMap,      GetInp_CrAsyConf,         GetInp_CrBank,     &
!                        GetInp_CrBank,         GetInp_CrPosition,        GetInp_CrPosChg
IMPLICIT NONE
INTEGER           :: indev,outdev
INTEGER           :: idcard
INTEGER,parameter :: idblock = 17
INTEGER :: IntTemp(100)
INTEGER :: i,k,nb
CHARACTER(10) :: Chartemp
LOGICAL :: Master

master = PE%master
CALL DMALLOC(CrCspGenCntl%RodOutMap, nz)
CALL DMALLOC(CrCspGenCntl%RodInMap, nz)
DO while(TRUE)
  read(indev,'(a256)') oneline
  IF(master) CALL message(io8,FALSE,FALSE,oneline)
  IF(probe .EQ. BANG) cycle;     IF(probe .EQ. POUND) cycle
  IF(oneline .EQ. BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe .EQ. DOT) exit;       IF(probe.eq.SLASH) exit
  read(oneline,*) cardname; CALL toupper(cardname)
  idcard = FindCardId(idblock,cardname); IF(idcard .eq. 0) exit
  nLineField = nfields(oneline)-1
  SELECT CASE(idcard)
    CASE(1)    !PINXS_FILE
      CrCspGenCntl%nPinxsDat = CrCspGenCntl%nPinxsDat + 1
      READ(oneline, *) astring, k
      CALL GetFn(oneline, 3, CrCspGenCntl%PinXsFn(k))
    CASE(2)    !XS_MAP
      READ(oneline, *) astring, CharTemp
      CALL TOUPPER(CharTemp)
      IF(CharTemp .EQ. 'IN') THEN
        READ(oneline, *) astring, chartemp, CrCspGenCntl%RodINMap(1:nz)
      ELSEIF(CharTemp .EQ. 'OUT') THEN
        READ(oneline, *) astring, chartemp, CrCspGenCntl%RodOutMap(1:nz)
      ELSE
        CALL TERMINATE('Wrong Input CARD : XS_MAP')
      ENDIF
    CASE(3)    !CSP_BANK
      READ(oneline, *) astring, CrCspGenCntl%BankId
    CASE(4)    !CSP_PLN
      READ(oneline, *) astring, CrCspGenCntl%CrPln(1:2)
  END SELECT
ENDDO
Backspace(indev); IF(Master) Backspace(outdev)
nTracerCntl%lCrInfo = .TRUE.
nTracerCntl%lProblem = lCrCspGen
END SUBROUTINE

SUBROUTINE ReadXeDynCard(indev,outdev)
USE PARAM
USE CNTL,           ONLY : nTracerCntl
USE PE_MOD,         ONLY : PE
USE XeDyn_Mod,      ONLY : GetInp_TimeUnit,      GetInp_TimeStep,               &
                           GetInp_CoreState,     ProcessXeDynInp
USE ALLOCS
IMPLICIT NONE
INTEGER           :: indev,outdev
INTEGER           :: idcard
INTEGER,parameter :: idblock = 18
LOGICAL :: Master
DO while(TRUE)
  read(indev,'(a256)') oneline
  IF(master) CALL message(io8,FALSE,FALSE,oneline)
  IF(probe .EQ. BANG) cycle;     IF(probe .EQ. POUND) cycle
  IF(oneline .EQ. BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe .EQ. DOT) exit;       IF(probe.eq.SLASH) exit
  read(oneline,*) cardname; CALL toupper(cardname)
  idcard = FindCardId(idblock,cardname); IF(idcard .eq. 0) exit
  nLineField = nfields(oneline)-1
  SELECT CASE(idcard)
    CASE(1)    !Time Step
      CALL GetInp_TimeStep(oneline)
    CASE(2)    !Unit
      CALL GetInp_TimeUnit(oneline)
    CASE(3)    !CoreState
      CALL GetInp_CoreState(oneline)
  END SELECT
ENDDO
nTracerCntl%lProblem = lXenonDynamics
CALL ProcessXeDynInp()
Backspace(indev); IF(Master) Backspace(outdev)
END SUBROUTINE

SUBROUTINE ReadMcpRestartCard(indev, outdev)
USE PARAM
USE CNTL,           ONLY : nTracerCntl
USE PE_MOD,         ONLY : PE
USE ALLOCS
USE MCP_Util,       ONLY : lMCP_Restart
IMPLICIT NONE
INTEGER           :: indev,outdev
INTEGER           :: idcard
INTEGER,parameter :: idblock = 20
LOGICAL :: master
master = PE%master
DO while(TRUE)
  read(indev,'(a256)') oneline
  IF(master) CALL message(io8,FALSE,FALSE,oneline)
  IF(probe .EQ. BANG) cycle;     IF(probe .EQ. POUND) cycle
  IF(oneline .EQ. BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe .EQ. DOT) exit;       IF(probe.eq.SLASH) exit
  read(oneline,*) cardname; CALL toupper(cardname)
  idcard = FindCardId(idblock,cardname); IF(idcard .eq. 0) exit
  nLineField = nfields(oneline)-1
  SELECT CASE(idcard)
    CASE(1) ! nMCP
    CASE(2) ! MCP_LIST
      CALL ReadMcpLib(indev, outdev)
      !lMCP_Restart = .TRUE.
    CASE(3) ! MCP_FILE
      CALL TERMINATE('Improper location of MCP_FILE CARD')
  ENDSELECT
ENDDO
Backspace(indev); IF(Master) Backspace(outdev)
END SUBROUTINE

SUBROUTINE ReadMcpLib(indev, outdev)
USE PARAM
USE PE_MOD, ONLY : PE
USE CNTL,   ONLY : nTracerCntl
USE MCP_Util, ONLY : MCP_Path
IMPLICIT NONE
INTEGER           :: indev,outdev
INTEGER           :: idcard
INTEGER,parameter :: idblock = 20
INTEGER :: id
LOGICAL :: Master
master = PE%master
DO while(TRUE)
  read(indev,'(a256)') oneline
  IF(probe .EQ. BANG) cycle;     IF(probe .EQ. POUND) cycle
  IF(oneline .EQ. BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe .EQ. DOT) exit;       IF(probe.eq.SLASH) exit
  read(oneline,*) cardname; CALL toupper(cardname)
  idcard = FindCardId(idblock,cardname); IF(idcard .NE. 3) exit
  nLineField = nfields(oneline)-1
  read(oneline, *) astring, id
  CALL GETFN(oneline, 3, MCP_Path(id))
  CONTINUE
ENDDO
Backspace(indev); IF(Master) Backspace(outdev)
END SUBROUTINE

SUBROUTINE ReadSubchOptionCard(InDev, OutDev)
USE PARAM
USE PE_MOD,         ONLY : PE
USE SubChCoupling_mod,    ONLY: coupled_maxsteps, coupled_relconv, Courant, sbch_outop
IMPLICIT NONE

INTEGER :: InDev, OutDev
INTEGER           :: idcard, I
INTEGER,parameter :: idblock = 21
LOGICAL :: Master

Master = PE%master

DO while(TRUE)
  read(indev,'(a256)') oneline
  IF(master) CALL message(io8,FALSE,FALSE,oneline)
  IF(probe .EQ. BANG) cycle;     IF(probe .EQ. POUND) cycle
  IF(oneline .EQ. BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe .EQ. DOT) exit;       IF(probe.eq.SLASH) exit
  read(oneline,*) cardname; CALL toupper(cardname)
  idcard = FindCardId(idblock,cardname); IF(idcard .eq. 0) exit
  nLineField = nfields(oneline)-1
  SELECT CASE(idcard)
    CASE(1)  ! CONTROL

    CASE(2)  ! CMAXOUTER
      READ(oneline, *) astring, coupled_maxsteps
    CASE(3)  ! BORONTRACK

    CASE(4)  ! REL_CONV
      READ(oneline, *) astring, coupled_relconv

    CASE(5)  ! COURANT
      READ(oneline, *) astring, Courant
    CASE(6)  ! SBCH_OUT
      READ(oneline, *) astring, (sbch_outop(I), I = 1, nLineField)
  END SELECT
ENDDO

Backspace(indev); IF(Master) Backspace(io8)
END SUBROUTINE

!--- CNJ Edit : Intel MKL Option Parsing
SUBROUTINE ReadMKLCard(InDev, OutDev)
USE PARAM
USE PE_MOD,         ONLY : PE
#ifdef __INTEL_MKL
USE MKL_3D,         ONLY : mklGeom,         mklCntl,        mklDepl
#endif
IMPLICIT NONE

INTEGER :: InDev, OutDev
INTEGER           :: idcard
INTEGER,parameter :: idblock = 22
INTEGER           :: option = 1
LOGICAL :: Master

Master = PE%master
PE%lMKL = TRUE

#ifdef __INTEL_MKL
DO while(TRUE)
  read(indev,'(a256)') oneline
  IF(master) CALL message(io8,FALSE,FALSE,oneline)
  IF(probe .EQ. BANG) cycle;     IF(probe .EQ. POUND) cycle
  IF(oneline .EQ. BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe .EQ. DOT) exit;       IF(probe.eq.SLASH) exit
  read(oneline,*) cardname; CALL toupper(cardname)
  idcard = FindCardId(idblock,cardname); IF(idcard .eq. 0) exit
  nLineField = nfields(oneline)-1
  SELECT CASE(idcard)
    CASE(1)  ! DAVIDSON
      READ(oneline, *) astring, mklCntl%lDavidson
    CASE(2)  ! BILU

    CASE(3)  ! AXTYPE
      IF (nLineField .EQ. 1) THEN
        READ(oneline, *) astring, mklCntl%AxSolver
      ELSEIF (nLineField .EQ. 2) THEN
        READ(oneline, *) astring, mklCntl%AxSolver, option
      ELSEIF (nLineField .EQ. 3) THEN
        READ(oneline, *) astring, mklCntl%AxSolver, option, mklCntl%MOCHeight
      ELSEIF (nLineField .EQ. 4) THEN
        READ(oneline, *) astring, mklCntl%AxSolver, option, mklCntl%MOCHeight, mklGeom%nPolar1D
      ENDIF
      IF (mklCntl%AxSolver .EQ. 0) THEN
        mklCntl%lAxial = .FALSE.
      ELSEIF (mklCntl%AxSolver .EQ. 1) THEN
        IF (option .EQ. 1) mklCntl%lSENM = .TRUE.
      ELSEIF (mklCntl%AxSolver .EQ. 2) THEN
        IF (PE%nCMFDProc .GT. 1) STOP 'WARNING: MPI Execution is not Supported in Axial FDM'
      ELSEIF (mklCntl%AxSolver .EQ. 3) THEN
        mklCntl%lMOC = .TRUE.
        mklCntl%polyOrder = option
        IF (mklCntl%polyOrder .EQ. 1) mklCntl%lCASMO = .TRUE.
      ELSE
        STOP 'WARNING: Unknown Axial Solver Type'
      ENDIF
    CASE(4)  ! DCPL
      IF (nLineField .EQ. 1) THEN
        READ(oneline, *) astring, mklCntl%lDcpl
      ELSE
        READ(oneline, *) astring, mklCntl%lDcpl, mklCntl%DcplLv
      ENDIF
      IF (.NOT. mklCntl%lDcpl) mklCntl%DcplLv = 0
    CASE(5)  ! SHIFT
      READ(oneline, *) astring, mklCntl%lShift, mklCntl%Shift
    CASE(6)  ! PCMFD

    CASE(7)  ! DIRECT
      READ(oneline, *) astring, mklCntl%lDirect
      IF (PE%nProc .GT. 1) mklCntl%lDcpl = .FALSE.
    CASE(8)  ! CHEBYSHEV
      IF (nLineField .EQ. 1) THEN
        READ(oneline, *) astring, mklCntl%lChebyshev
        mklCntl%chebyOrder = mklCntl%maxOuter
      ELSE
        READ(oneline, *) astring, mklCntl%lChebyshev, mklCntl%chebyOrder
      ENDIF
    CASE(9)

    CASE(10) ! GCCMFD
      READ(oneline, *) astring, mklCntl%lGcCMFD
    CASE(11) ! GCSTRUCT
      mklGeom%ngc = nLineField / 2
      ALLOCATE(mklGeom%GcStruct(2, mklGeom%ngc))
      READ(oneline, *) astring, mklGeom%GcStruct
    CASE(12) ! SUPERPIN
      READ(oneline, *) astring, mklCntl%lSuperpin
    CASE(13) ! SUBPLN
      READ(oneline, *) astring, mklCntl%lSubplane, mklCntl%CMFDHeight
    CASE(14) ! DEPL
      READ(oneline, *) astring, mklCntl%lDepl, mklDepl%SysByte, mklDepl%Scale
  END SELECT
ENDDO
#endif

Backspace(indev); IF(Master) Backspace(io8)
END SUBROUTINE

!--- CNJ Edit : CUDA Option Parsing
SUBROUTINE ReadCUDACard(InDev, OutDev)
USE PARAM
USE PE_MOD,         ONLY : PE
USE CNTL,           ONLY : nTracerCntl
#ifdef __PGI
USE CUDA_MASTER,    ONLY : cuGeometry,      cuCntl,       cuDepl
#endif
IMPLICIT NONE

INTEGER :: InDev, OutDev
INTEGER           :: idcard
INTEGER,parameter :: idblock = 23
INTEGER           :: option = 1
LOGICAL :: Master

Master = PE%master
PE%lCUDA = TRUE
PE%lCUDACMFD = TRUE
nTracerCntl%lNodeMajor = TRUE
!nTracerCntl%lMacro = TRUE
nTracerCntl%lED = FALSE
nTracerCntl%lMLG = TRUE
IF (nTracerCntl%lXsAlign) THEN
 nTracerCntl%lMacro = FALSE
END IF

#ifdef __PGI
DO while(TRUE)
  read(indev,'(a256)') oneline
  IF(master) CALL message(io8,FALSE,FALSE,oneline)
  IF(probe .EQ. BANG) cycle;     IF(probe .EQ. POUND) cycle
  IF(oneline .EQ. BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe .EQ. DOT) exit;       IF(probe.eq.SLASH) exit
  read(oneline,*) cardname; CALL toupper(cardname)
  idcard = FindCardId(idblock,cardname); IF(idcard .eq. 0) exit
  nLineField = nfields(oneline)-1
  SELECT CASE(idcard)
    CASE(1)  ! CU_CMFD
      READ(oneline, *) astring, PE%lCUDACMFD
    CASE(2)  ! CU_GCCMFD
      READ(oneline, *) astring, cuCntl%lGcCMFD
    CASE(3)  ! CU_GCSTR
      cuGeometry%ngc = nLineField / 2
      ALLOCATE(cuGeometry%GcStruct(2, cuGeometry%ngc))
      READ(oneline, *) astring, cuGeometry%GcStruct
    CASE(4)  ! CU_DCPL
      IF (nLineField .EQ. 1) THEN
        READ(oneline, *) astring, cuCntl%lDcpl
        IF (cuCntl%lDcpl) cuCntl%DcplLv = 2
      ELSEIF (nLineField .EQ. 2) THEN
        READ(oneline, *) astring, cuCntl%lDcpl, cuCntl%DcplLv
        IF (.NOT. cuCntl%lDcpl) cuCntl%DcplLv = 0
      ENDIF
    CASE(5)  ! CU_RBSOR
      READ(oneline, *) astring, cuCntl%lNatural
      cuCntl%lNatural = .NOT. cuCntl%lNatural
    CASE(6)  ! CU_AXTYPE
      IF (nLineField .EQ. 1) THEN
        READ(oneline, *) astring, cuCntl%AxSolver
      ELSEIF (nLineField .EQ. 2) THEN
        READ(oneline, *) astring, cuCntl%AxSolver, option
      ELSEIF (nLineField .EQ. 3) THEN
        READ(oneline, *) astring, cuCntl%AxSolver, option, cuCntl%MOCHeight
      ELSEIF (nLineField .EQ. 4) THEN
        READ(oneline, *) astring, cuCntl%AxSolver, option, cuCntl%MOCHeight, cuGeometry%nPolar1D
      ENDIF
      IF (cuCntl%AxSolver .EQ. 0) THEN
        cuCntl%lAxial = .FALSE.
      ELSEIF (cuCntl%AxSolver .EQ. 1) THEN
        IF (option .EQ. 0) cuCntl%lSENM = .FALSE.
      ELSEIF (cuCntl%AxSolver .EQ. 3) THEN
        IF (option .EQ. 0) cuCntl%lCASMO = .FALSE.
      ENDIF
    CASE(7)  ! CU_SPAI
      READ(oneline, *) astring, cuCntl%lSPAI
    CASE(8)  ! CU_AAJ

    CASE(9)  ! CU_SUBPLN
      READ(oneline, *) astring, cuCntl%lSubplane, cuCntl%CMFDHeight
    CASE(10) ! CU_SHIFT
      READ(oneline, *) astring, cuCntl%lShift, cuCntl%Shift
    CASE(11) ! CU_DEPL
      READ(oneline, *) astring, PE%lCUDADepl, cuDepl%sysByte, cuDepl%Scale
    CASE(12) ! CU_XS
      READ(oneline, *) astring, nTracerCntl%lXsAlign
      IF (nTracerCntl%lXsAlign) THEN
        nTracerCntl%lMacro = FALSE
      ELSE
        nTracerCntl%lMacro = TRUE
      END IF
    CASE(13) ! CU_SUPERPIN
      READ(oneline, *) astring, cuCntl%lsuperpin
    CASE(14) ! CU_PWDIST
      READ(oneline, *) astring, cuCntl%lPwDist
  END SELECT
ENDDO
#endif

Backspace(indev); IF(Master) Backspace(io8)
END SUBROUTINE

SUBROUTINE ReadGCCard(indev, outdev)
USE PARAM
USE CNTL,   ONLY : nTracerCntl
USE Files,  ONLY : caseid,            filename
USE IOUTIL, ONLY : Getfn
USE PE_MOD, ONLY : PE
USE ALLOCS
USE FXRVAR_MOD !--- BYS edit
USE B1_Mod, ONLY : b1k, lb1k
IMPLICIT NONE
INTEGER           :: indev,outdev
INTEGER           :: idcard
INTEGER,parameter :: idblock = 24
INTEGER :: IntTemp(100), ipos1(100), ipos2(100)
INTEGER :: i,k,nb, nspt1, nspt2
REAL :: REALTemp(100)
INTEGER :: nspt, ipos(100)
LOGICAL :: Master
CHARACTER(80) :: CharTemp(10)


master = PE%master
nTracerCntl%GC_spec=2  !Solution Spectrum as Default
nTracerCntl%nGCgrp=2   !two-group GC as default
nTracerCntl%lPIN_GCL=.FALSE.
DO while(TRUE)
    read(indev,'(a256)') oneline
    IF(master) CALL message(io8,FALSE,FALSE,oneline)
    IF(probe .EQ. BANG) cycle;     IF(probe .EQ. POUND) cycle
    IF(oneline .EQ. BLANK) cycle;  IF(IFnumeric(oneline)) cycle
    IF(probe .EQ. DOT) exit;       IF(probe.eq.SLASH) exit
    read(oneline,*) cardname; CALL toupper(cardname)
    idcard = FindCardId(idblock,cardname); IF(idcard .eq. 0) exit
    nLineField = nfields(oneline)-1


    SELECT CASE(idcard)
        CASE(1)   ! ASY_GCL : Assembly-wise Group Constants
            nTracerCntl%lASY_MIC = .TRUE.
            SELECT CASE(nLineField)
                CASE(1)
                    READ(oneline, *) astring, nTracerCntl%lASY_GCL
                CASE(2)
                    READ(oneline, *) astring, nTracerCntl%lASY_GCL, nTracerCntl%lASY_MIC
            ENDSELECT
        CASE(2)   ! PIN_GCL : Pin-wise Group Constants for SPHINCS
            nTracerCntl%lGCgapHom= .TRUE.
            nTracerCntl%lPIN_MIC = .FALSE.
            SELECT CASE(nLineField)
                CASE(1)
                    READ(oneline, *) astring, nTracerCntl%lPIN_GCL
                CASE(2)
                    READ(oneline, *) astring, nTracerCntl%lPIN_GCL, nTracerCntl%lGCgapHom
                CASE(3)
                    READ(oneline, *) astring, nTracerCntl%lPIN_GCL, nTracerCntl%lGCgapHom, nTracerCntl%lPIN_MIC
            ENDSELECT
        CASE(3)   ! PDQ_GCL : PDQ format for nTRACER-MASTER
            READ(oneline, *) astring, nTracerCntl%lPDQ_GCL
        CASE(4)   ! NTIG_RST : nTIG Restart file, nTRACER-INPUT-GENERATOR
            READ(oneline, *) astring, nTracerCntl%lGCRST
        CASE(5)   ! SPEC : GC_GEN Spectrum / 0: B1 critical / 1: B1D-inf spec / 2: SA spec
            READ(oneline, *) astring, nTracerCntl%GC_spec
        CASE(6)   ! NFEWG : # of Few-Group / 1, 2, 4, 8 available
            READ(oneline, *) astring, nTracerCntl%nGCgrp
        CASE(7)   ! TMOD : moderator temperature variation calc with given delta value
            READ(oneline, *) astring, tmod(1:nLineField)
            ntmod=nLineField
            lFXRVAR=.TRUE.
            ltmod=.TRUE.
            lBranchRun=.TRUE.
        CASE(8)   ! TFUEL : fuel temperature variation calc with given delta value
            READ(oneline, *) astring, tfuel(1:nLineField)
            ntfuel=nLineField
            ltfuel=.TRUE.
            lBranchRun=.TRUE.
        CASE(9)   ! RHO : moderator density variation calc with given modified density
            READ(oneline, *) astring, rho(1:nLineField)
            nrho=nLineField
            lFXRVAR=.TRUE.
            lrho=.TRUE.
            lBranchRun=.TRUE.
        CASE(10)  ! BORON : moderator density variation calc with given modified density
            READ(oneline, *) astring, bboron(1:nLineField)
            IF(nTracerCntl%BoronPPM .GT. 0._8) nTracerCntl%lInitBoron = .TRUE.
            nboron=nLineField
            lFXRVAR=.TRUE.
            lboron=.TRUE.
            lBranchRun=.TRUE.
        CASE(11)  ! EFT : Effective Fuel Temperature search
            !SELECT CASE(nLineField)
            !    CASE(1)
            !        READ(oneline, *) astring, nTracerCntl%lEFTsearch
            !    CASE DEFAULT
            !ENDSELECT
            nTracerCntl%nEFTpts=nLineField-2
            READ(oneline, *) astring, nTracerCntl%lEFTsearch, nTracerCntl%crit_EFT, nTracerCntl%P_EFT(1:nLineField-2)
            IF( nTracerCntl%lEFTsearch )THEN
                nTracerCntl%lProblem=lEFTsearch
                nTracerCntl%lFeedBack=.TRUE.
                !nTracerCntl%fmdotfa=1E20 ! infinite mass flow rate
                nTracerCntl%lASY_MIC = .TRUE.
                nTracerCntl%lASY_GCL = .TRUE.
            ENDIF
        CASE(12)  ! UNICOOL : EFT search option; Set uniform coolant temperature
            READ(oneline, *) astring, nTracerCntl%EFTunicoolant
        CASE(13)  ! EFT_BORON : EFT search option; Set uniform coolant temperature
            nTracerCntl%EFT_nboron=nLineField
            READ(oneline, *) astring, nTracerCntl%EFT_boron(1:nTracerCntl%EFT_nboron)
    END SELECT
ENDDO
IF( nTracerCntl%lASY_GCL .OR. nTracerCntl%lPIN_GCL .OR. nTracerCntl%lGCRST .OR. nTracerCntl%lPDQ_GCL )THEN
    nTracerCntl%lGC=.TRUE.
ENDIF
IF( (lTmod .OR. lTfuel .OR. lRho .OR. lboron) .AND. lBranchRun)THEN
    nTracerCntl%lBranch=.TRUE.
    nTracerCntl%lProblem=lBranch ! lBranch=7
ENDIF

Backspace(indev); IF(Master) Backspace(outdev)
RETURN
END SUBROUTINE

SUBROUTINE ReadNTIGCard(indev, outdev)
USE PARAM
USE PE_MOD, ONLY : PE
USE CNTL,   ONLY : nTracerCntl
USE GEOM,   ONLY : Core
USE LpShf_mod, ONLY : LpShf
IMPLICIT NONE
INTEGER           :: indev,outdev
INTEGER           :: idcard
INTEGER,parameter :: idblock = 25
LOGICAL :: Master
INTEGER :: i, j
INTEGER :: ixya, icyc
INTEGER :: jfr, jto, nxasy
INTEGER, ALLOCATABLE :: lfreshfuel(:)

master = PE%master
DO while(TRUE)
  read(indev,'(a256)') oneline
  IF(master) CALL message(io8,FALSE,FALSE,oneline)
  IF(probe .EQ. BANG) cycle;     IF(probe .EQ. POUND) cycle
  IF(oneline .EQ. BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe .EQ. DOT) exit;       IF(probe.eq.SLASH) exit
  read(oneline,*) cardname; CALL toupper(cardname)
  idcard = FindCardId(idblock,cardname); IF(idcard .eq. 0) exit
  SELECT CASE(idcard)
    CASE(1) !HM_MASS
      READ(oneline, *), astring, Core%HM_mass0
    CASE(2) !FRESH_ASM
      nTracerCntl%lCooling = .TRUE.
      ALLOCATE(lfreshfuel(Core%nxya))
      ALLOCATE(LpShf%SHF(Core%nxya))
      jfr = 1;
      DO i = 1, Core%nya
        READ(indev,'(A)'), ONELINE
        nxasy = nfields(ONELINE)
        jto = jfr + nxasy - 1
        READ(ONELINE,*), (lfreshfuel(j), j = jfr, jto)
        DO j = jfr, jto
          IF (lfreshfuel(j) .EQ. 0) THEN
            LpShf%Shf(j)%lfreshfuel = .FALSE.
          END IF
        END DO
        jfr = jto + 1
      END DO
    CASE(3) !PUL
      nTracerCntl%lCooling = .TRUE.
      LpShf%lPUL = .TRUE.
      READ(ONELINE, *), astring, LpShf%PUL(1)
    END SELECT
ENDDO
BACKSPACE(indev); IF(Master) BACKSPACE(outdev)

END SUBROUTINE

END module
