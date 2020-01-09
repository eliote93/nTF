#include <defines.h>
MODULE CNTLROD_mod
USE PARAM
USE TYPEDEF,            ONLY : CrCell_TYPE,          CrAsyConf_TYPE,          CrBank_TYPE,           &
                               CrInfo_Type,          CrSearch_Type,           CrPosDat_Type,        &                                                         
                               CoreInfo_Type,        Cell_Type,               PE_Type
USE CNTL,               ONLY : nTracerCntl_Type
USE ALLOCS
IMPLICIT NONE
INTEGER, PRIVATE, PARAMETER :: C_CRCELL = 1
INTEGER, PRIVATE, PARAMETER :: C_CRASYCONF = 2
INTEGER, PRIVATE, PARAMETER :: C_CRBANK = 3

INTEGER, PRIVATE, PARAMETER :: C_CrStepInp = 0
INTEGER, PRIVATE, PARAMETER :: C_CrPosInp = 1

INTEGER, PRIVATE, PARAMETER :: nCrBankMax = 100
INTEGER, PRIVATE, PARAMETER :: nCrPosDatMax = 100
REAL, PRIVATE :: hzmax


!Control Rod Geometry
INTEGER, PRIVATE :: nCrBank0 = 0
INTEGER, PRIVATE :: nCrCell0 = 0
INTEGER, PRIVATE :: nCrAsyConf0 = 0
INTEGER, PRIVATE :: nCrBank = 0
INTEGER, PRIVATE :: nCrCell = 0
INTEGER, PRIVATE :: nCrAsyConf = 0
LOGICAL, PRIVATE :: lCrInfo = .FALSE.
LOGICAL, PRIVATE :: lCrIn = .FALSE.

!INTEGER, PARAMETER :: iCrOutIdx = 1
!INTEGER, PARAMETER :: iCrInIdx = 2


INTEGER, POINTER :: CRCellMAP(:, :)
TYPE(CrAsyConf_Type), POINTER, PRIVATE :: CrAsyConf(:)
TYPE(CrBank_Type), POINTER, PRIVATE :: CrBank(:)
TYPE(CrPosDat_Type), TARGET :: CrPosDat
TYPE(CrInfo_Type), POINTER, PRIVATE :: CrInfo
LOGICAL, PRIVATE :: LAlloc(3) = .FALSE.
CONTAINS

SUBROUTINE AddNCrConf(imod)
IMPLICIT NONE
INTEGER :: imod
IF(imod .EQ. C_CRCELL) nCrCell0 = nCrCell0 + 1
IF(imod .EQ. C_CRASYCONF) nCrAsyConf0 = nCrAsyConf0 + 1
IF(imod .EQ. C_CRBANK) nCrBank0 = nCrBank0 + 1
END SUBROUTINE

SUBROUTINE GetInp_CrCellMap(CellUCR, CellCR)
IMPLICIT NONE
INTEGER :: CellUCR, CellCR
IF(.NOT. lAlloc(C_CRCELL)) THEN
  lAlloc(C_CRCELL) = .TRUE.
  CALL DMALLOC(CRCellMAP, 2, nCrCell0)
ENDIF

nCrCell = nCrCell + 1
CrCellMap(1, nCrCell) = CellUCR
CrCellMap(2, nCrCell) = CellCR
END SUBROUTINE

SUBROUTINE GetInp_CrAsyConf(ConfId, iang, nxc0, lgap, indev, outdev, PE)
USE ioutil,         only : terminate,   toupper,       IFnumeric,   nfields,     message
USE inputcards ,    only : oneline,     probe
USE BasicOperation, ONLY : CP_CA,       CP_VA
IMPLICIT NONE
INTEGER :: indev, outdev
INTEGER :: ConfId, iang, nxc0h
LOGICAL :: lGap
TYPE(PE_TYPE) :: PE

INTEGER :: ndat, nxc0, nxc
INTEGER :: i,j, k, n, jfr, jto, jfr0, jto0, jx,iy
INTEGER :: inpdat(6400), AsyConf(80,80)

IF(.NOT. lAlloc(C_CRASYCONF)) THEN
  lAlloc(C_CRASYCONF) = .TRUE.
  ALLOCATE(CrAsyConf(nCrAsyConf0))
ENDIF
nxc0h = nxc0 / 2 + mod(nxc0, 2)
ndat = nxc0 * nxc0
IF(iang .EQ. 90) ndat = nxc0h * nxc0h
IF(iang .EQ. 45) ndat = nxc0h * (1+nxc0h) / 2 

jfr = 1; jto = 0

DO WHILE(.TRUE.)
  read(indev, '(a256)') oneline
  IF(PE%Master) CALL message(outdev,FALSE, FALSE,oneline)
  IF(probe.eq.BANG) cycle;     IF(probe.eq.POUND) cycle  
  n = nFields(oneline)
  jto = jfr + n -1
  read(oneline, *) (inpdat(j), j = jfr, jto)
  jfr = jfr + n 
  IF(jto .eq. ndat) exit  
ENDDO

nxc = nxc0; jfr = 1;jto = nxc;
IF(lGap) THEN
  nxc = nxc+2
  jfr = 2; jto = nxc -1
ENDIF

CALL CP_CA(AsyConf(1:nxc+2, 1:nxc+2), 0, nxc+2, nxc+2)

SELECT CASE(iAng)
  case(360)
    k=0
    DO i = jfr, jto
      DO j = jfr, jto
        k = k + 1
        AsyConf(j, i) = inpdat(k)  
      ENDDO
    ENDDO
  case(90)
    k=0
    jto0 = jto; jfr0 = jto - nxc0h +1
    DO i = jfr0, jto0
      DO j = jfr0, jto0
        k = k + 1
        AsyConf(j, i) = inpdat(k)
      ENDDO
    ENDDO
  CASE(45)
    k=0
    jto0 = jto; jfr0 = jto - nxc0h +1
    DO i = jfr0, jto0
      DO j = jfr0, i !jfr0+(i-jfr0+1)
        k = k + 1
        AsyConf(j, i) = inpdat(k)
      ENDDO
    ENDDO
    !45 -> 90 Symmetry
    k= ndat+1 
    DO i = jto0, jfr0, -1
      DO j = i, jfr0,-1
      !DO j = jto0, i-1,-1
         k=k-1
         AsyConf(i, j) = inpdat(k)
      ENDDO
    ENDDO         
  END SELECT
  
IF( iAng .NE. 360) then
  !90 Symmetry -> 360 
  DO i = jfr0, jto0
    iy = nxc - i + 1 
    DO j = jfr0, jto0
      jx = nxc - j + 1 
      AsyConf(jx, i) = AsyConf(j, i)
      AsyConf(j, iy) = AsyConf(j, i)
      AsyConf(jx, iy) = AsyConf(j, i)
    ENDDO 
  ENDDO
ENDIF
  
nCrAsyConf = nCrAsyConf + 1
CALL DMALLOC(CrAsyConf(ConfId)%CrLoc, nxc, nxc)
CALL CP_VA(CrAsyConf(ConfId)%CrLoc(1:nxc, 1:nxc), AsyConf(1:nxc, 1:nxc), nxc, nxc)

END SUBROUTINE

SUBROUTINE GetInp_CrBank(BankId, BankName, nxya, nxa, nya, indev, outdev, PE)
USE ioutil,         only : terminate,   toupper,       IFnumeric,   nfields,     message
USE inputcards ,    only : oneline,     probe
USE BasicOperation, ONLY : CP_CA,       CP_VA
IMPLICIT NONE
INTEGER :: BankId
CHARACTER(10) :: BankName
INTEGER :: nxya, nxa, nya
INTEGER :: indev, outdev
TYPE(PE_TYPE) :: PE

INTEGER :: ndat, ncrasy, nFieldsLine, ixa, iya, ixya, i
INTEGER :: inpdat(nxa*nya)

IF(.NOT. lAlloc(C_CRBANK)) THEN
  lAlloc(C_CRBANK) = .TRUE.
  ALLOCATE(CrBank(nCrBank0))
ENDIF
nCrBank = nCrBank + 1
CALL DMALLOC(CrBank(BankId)%CrAsyLoc,nxya)
ndat = 0; ixa = 0; iya = 0; ixya=0
DO WHILE(.TRUE.)
  read(indev, '(a256)') oneline
  IF(PE%Master) CALL message(outdev,FALSE, FALSE,oneline)
  IF(probe.eq.BANG) cycle;     IF(probe.eq.POUND) cycle  
  iya = iya + 1
  nFieldsLine = nfields(oneline)
  READ(oneline, *) (inpdat(i), i = 1, nFieldsLine)
  
  DO ixa = 1, nFieldsLine
    ixya = ixya + 1
    CrBank(BankID)%CrAsyLoc(ixya)=inpdat(ixa)
    IF(inpdat(ixa) .NE. 0) ndat = ndat + 1
  ENDDO
  IF(iya .EQ. nya) EXIT
ENDDO

CrBank(BankID)%BankName = BankName
CrBank(BankID)%nCrAsy = ndat
END SUBROUTINE

SUBROUTINE GetInp_CrPosChg(Oneline, nFields, PE)
USE TYPEDEF,           ONLY : PE_TYPE
USE IOUTIL,            ONLY : TERMINATE
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
INTEGER :: nFields
CHARACTER(256) :: Oneline 

CHARACTER(50) :: ErrMesg
CHARACTER(20) :: astring(2)
LOGICAL :: lStepInp
REAL :: RealDat(5)
INTEGER :: i, idat, BankId

IF(nFields .GT. 1) THEN
  READ(oneline, *) astring(1), BANKID, astring(2)
  SELECT CASE(astring(2)) 
    CASE('POS')
      lStepInp = .FALSE.
    CASE('STEP')
      lStepInp = .TRUE.
    CASE DEFAULT
      CALL TERMINATE('Wrong CR_POS Input Line')
  END SELECT

  IF(lStepInp) THEN
    
  ELSE
    READ(oneline, *) astring(1), BankId, astring(2), RealDat(1)
    CrBank(BankId)%RodPos = RealDat(1)
  ENDIF
ELSE
  READ(oneline, *) astring(1), idat
  IF(.NOT. CrPosDat%LExist(idat)) THEN
    WRITE(ErrMesg, '(A, I5)') 'CNTL Input Error - Control Rod Position Data is not defined :', idat
    CALL TERMINATE('Control Rod Position is not Defined')
  ENDIF
  DO i = 1, nCrBank
    IF(CrBank(i)%PosType .EQ. 1) THEN !STATE
      CrBank(i)%RodStep = CrPosDat%RodStep(i, idat)
    ELSE
      CrBank(i)%RodPos = CrPosDat%RodPos(i, idat)
    ENDIF
  
  ENDDO
ENDIF
CrInfo%lCrChg = .TRUE.

END SUBROUTINE 

SUBROUTINE GetInp_CrPositionDat(Oneline, nFields, PE)
USE TYPEDEF,           ONLY : PE_TYPE
USE IOUTIL,            ONLY : FndChara,     TERMINATE
IMPLICIT NONE
CHARACTER(256) :: Oneline
INTEGER :: nFIelds
TYPE(PE_Type) :: PE

CHARACTER(50) :: ErrMesg
CHARACTER(20) :: astring
INTEGER :: i, j, idat, ibank
INTEGER :: ipos(200), nspt

READ(oneline, *) astring, idat
IF(idat .GT. nCrPosDatMax) THEN
  CALL TERMINATE('CNTL Input Error - Maximum Number of CR_POSDAT 100')
ENDIF

CALL FndChara(oneline, ipos, nspt, SLASH)
IF(nCrBank .NE. nspt - 1) THEN
  WRITE(ErrMesg, '(A, I5)') 'CNTL Input Error - CR_POSDAT : ', idat
  CALL TERMINATE(ErrMesg)  
ENDIF
CrPosDat%nDat = MAX(CrPosDat%nDat, idat)
CrPosDat%lExist(idat) = .TRUE.
DO i = 1, nspt-1
  READ(oneline(ipos(i)+1:ipos(i+1)-1), *, ERR=1000) ibank 
  READ(oneline(ipos(i)+1:ipos(i+1)-1), *, ERR=1000) j, CrPosDat%PosInpDat(ibank, idat)
  CrPosDat%RodPos(ibank, idat) = CrPosDat%PosInpDat(ibank, idat)
  CrPosDat%RodStep(ibank, idat) = NINT(CrPosDat%PosInpDat(ibank, idat))
ENDDO
1000 WRITE(ErrMesg, '(A, I5)') 'CNTL Input Error - CR_POSDAT : ', idat
END SUBROUTINE

SUBROUTINE GetInp_CrPosition(Oneline, nFields, PE)
USE TYPEDEF,           ONLY : PE_TYPE
USE IOUTIL,            ONLY : TERMINATE
IMPLICIT NONE
TYPE(PE_TYPE) :: PE
INTEGER :: nFields
CHARACTER(256) :: Oneline 

CHARACTER(20) :: astring(10)
REAL :: RealDat(5)
INTEGER :: BankID

LOGICAL :: lCrSearch, lStepInp

lCrSearch = .FALSE.; lStepInp = .FALSE.

READ(oneline, *) astring(1), BANKID, astring(2)
SELECT CASE(astring(2)) 
  CASE('POS')
    lStepInp = .FALSE.
    CrBank(BankId)%PosType = 2
  CASE('STEP')
    lStepInp = .TRUE.
    CrBank(BankId)%PosType = 1
  CASE DEFAULT
    CALL TERMINATE('Wrong CR_POS Input Line')
END SELECT

IF(lStepInp) THEN
    
ELSE
  IF(nFields .GT. 4) THEN
    READ(oneline, *) astring(1), BankId, astring(2), RealDat(1), RealDat(2), lCrSearch
  ELSE
    READ(oneline, *) astring(1), BankId, astring(2), RealDat(1), RealDat(2)
  ENDIF
  CrBank(BankId)%lCrSearch = lCrSearch
  CrBank(BankId)%RodInPos = RealDat(1)
  CrBank(BankId)%RodPos = RealDat(2)
ENDIF

END SUBROUTINE 

SUBROUTINE GetInp_CrMvDom(OneLine, PE)
IMPLICIT NONE
CHARACTER(256) :: oneline
TYPE(PE_TYPE) :: PE
CHARACTER(20) :: astring

READ(oneline, *) astring, CrInfo%CrSearch%CrMv_Dom(1:2)
END SUBROUTINE

SUBROUTINE AssignCRPosDat(idat)
IMPLICIT NONE
INTEGER :: idat, i

DO i =  1, nCrBank
  IF(CrBank(i)%PosType .EQ. 1) THEN !STATE
    CrBank(i)%RodStep = CrPosDat%RodStep(i, idat)
  ELSE
    CrBank(i)%RodPos = CrPosDat%RodPos(i, idat)
  ENDIF    
ENDDO
CrInfo%lCrChg = .TRUE.
  
END SUBROUTINE

SUBROUTINE InitCntlRodConf(Core, FmInfo, CmInfo, GroupInfo, nTRACERCntl, PE)
USE TYPEDEF,         ONLY : FmInfo_Type,        CmInfo_Type,        GroupInfo_Type
USE CrCsp_MOD,       ONLY : InitCrCsp
USE ioutil,          ONLY : message
USE FILES,           ONLY : io8
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTRACERCntl_Type):: nTracerCntl
TYPE(PE_Type) :: PE

IF(PE%MASTER) THEN
  WRITE(mesg, '(A)') 'Initialize Control Rod Features...'
  CALL message(io8, .TRUE., .TRUE., MESG)
ENDIF

ALLOCATE(CrInfo)
Core%lCrInfo = .TRUE.
Core%CrInfo=> CrInfo
Core%CrInfo%lCrChg = .TRUE.
CrInfo%nCrCell = nCrCell; CrInfo%nCrAsyConf = nCrAsyConf; CrInfo%nCrBank = nCrBank;
CrInfo%CrCellMap => CrCellMap; CrInfo%CrAsyConf => CrAsyConf; CrInfo%CrBank => CrBank
CrInfo%CrPosDat => CrPosDat
hzmax = sum(core%hz)

CALL SetCrCellInfo(Core%CellInfo, Core%nCellType, Core%lEdge)
CALL SetCrAsyConf(Core)
CALL SetCrBank(Core)
CALL SetGlobalCrPinInfo(Core)
CALL SetCrFxrInfo(Core, FmInfo, GroupInfo, nTracerCntl, PE)
CALL SetCspFxrInfo(Core, FmInfo, GroupInfo, nTracerCntl, PE)

IF(CrInfo%CrSearch%CrMv_Dom(1) .EQ. 0) CrInfo%CrSearch%CrMv_Dom =(/1, Core%nz/)

IF(nTracerCntl%lCrCsp) THEN
  CALL InitCrCsp(Core, GroupInfo, nTracerCntl, PE)
ENDIF
CONTINUE
END SUBROUTINE



SUBROUTINE SetCrBankPosition(Core, FmInfo, CmInfo, GroupInfo, nTracerCntl, PE)
USE TYPEDEF,      ONLY : FmInfo_Type,        CmInfo_Type,        GroupInfo_Type,          &
                         PE_Type,                                                         &
                         FxrInfo_Type,       Pin_Type,           Cell_Type
USE CrCsp_Mod,    ONLY : SetCrCspFtn
!USE XsTypeDef,    ONLY : ChangeMixture
USE FILES,        ONLY : io8
USE IOUTIL,       ONLY : message
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTRACERCntl_Type):: nTracerCntl
TYPE(PE_Type) :: PE

REAL :: CrPos0
INTEGER :: nxy, myzb, myze, xyb, xye   !--- CNJ Edit : Domain Decomposition + MPI
INTEGER :: FxrIdxSt, nFxrLocal, nCrFxr
INTEGER :: ibank, iz, ixy, ifxr, icel, ifxrbeg, ifxrend
INTEGER :: icrincell, icroutcell, imix_crin, imix_crout
INTEGER :: i, j
LOGICAL :: lcspftn
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)

IF(PE%MASTER) CALL Message(io8, TRUE, TRUE, 'Set Control Rod Position...')

myzb = PE%myzb; myze = PE%myze; nxy = Core%nxy

!-- CNJ Edit : Domain Decomposition + MPI
xyb = PE%myPinBeg; xye = PE%myPinEnd
IF (PE%RTMASTER) THEN
  xyb = 1; xye = nxy
ENDIF

Pin => Core%Pin; CellInfo => Core%CellInfo; Fxr => FmInfo%Fxr
DO ibank = 1, nCrBank
  CALL CalCrBankLocation(CrBank(ibank)%FracCrIn, CrBank(ibank)%lCrIn, CrBank(ibank)%RodInPos,      &
                         CrBank(ibank)%RodPos, Core%hz, Core%nz, nTracerCntl%lCrCsp)
  CrPos0 = CrBank(ibank)%RodInPos + CrBank(ibank)%RodPos
  IF(CrBank(ibank)%RodPos .LE. 0.001) CrBank(ibank)%lCrFullIn = .TRUE.
  IF(PE%MASTER) THEN
    IF(.NOT. CrBank(ibank)%lCrIn) THEN
      WRITE(mesg, '(5x, A10,A2,x, A15)') CrBank(ibank)%BankName,':', 'Rod Out'
    ELSE
      WRITE(mesg, '(5x, A10,A2,x, F12.3,A3, x, A)') CrBank(ibank)%BankName,':', CrBank(ibank)%RodPos, 'cm', 'from Fuel Bottom'
    ENDIF
    CALL message(io8, FALSE, TRUE, mesg)
  ENDIF
  DO ixy = xyb, xye
    IF(.NOT. Pin(ixy)%lCrPin) CYCLE
    IF(Pin(ixy)%CrBankId .NE. iBank) CYCLE
    DO iz = myzb, myze
      icel = Pin(ixy)%Cell(iz)
      IF(.NOT. CellInfo(icel)%lCrCell) CYCLE
      IF(ABS(CrBank(ibank)%FracCrIn(iz) - Pin(ixy)%FracCrIn(iz)) .LT. epsm2) CYCLE
      nFxrLocal = CellInfo(icel)%nFxr; FxrIdxSt = Pin(ixy)%FxrIdxSt
      nCrFxr = CellInfo(icel)%CrCell%nCrFxr
      IF(CellInfo(icel)%CrCell%lCrIn) THEN
        iCrInCell = icel; iCrOutCell = CellInfo(icel)%CrCell%CellCrIdx
      ELSE
        iCrOutCell = CellInfo(icel)%CrCell%CellCrIdx; iCrInCell = icel
      ENDIF
      
      lcspftn = .FALSE.
      DO i = 1, nCrFxr 
        j = CellInfo(icel)%CrCell%CrFxrIdx(i)
        ifxr = FxrIdxSt + j - 1
        imix_crin = CellInfo(icel)%CrCell%imix_Crin(i)
        imix_crout = CellInfo(icel)%CrCell%imix_Crout(i)
        IF(nTracerCntl%lXsLib) THEN
          Fxr(ifxr, iz)%lCrCspFtn = .FALSE.; Fxr(ifxr, iz)%lMixtureMixing = .FALSE.
          IF(CrBank(ibank)%FracCrIn(iz) .GT. 0.9999_8) THEN !out -> in
            CALL ReplaceMixture(Fxr(ifxr, iz), imix_crin)
          ELSEIF(CrBank(ibank)%FracCrIn(iz) .LT. 0.0001_8) THEN !in-> out
            CALL ReplaceMixture(Fxr(ifxr, iz), imix_crout)
          ELSE
            CALL CuspingMixture(Fxr(ifxr, iz), imix_crin, imix_crout, CrBank(ibank)%FracCrIn(iz))
            IF(.NOT. nTracerCntl%lCrCspFtn) Fxr(ifxr, iz)%lCrCspFtn = .FALSE.
            lcspftn = .TRUE.
            !Fxr(ifxr, iz)%lCrCspFtn = .FALSE.
            !CALL SetCrCspFtn(Fxr(ifxr:ifxr, iz), 1, CrPos0, CrBank(ibank)%FracCrIn(iz), Pin(ixy)%iasy, iz, nTracerCntl, PE)
          ENDIF
        ELSE
          IF(CrBank(ibank)%FracCrIn(iz) .GT. 0.5_8) THEN !out -> in
            Fxr(ifxr, iz)%imix = iMix_CrIn
          ELSE
            Fxr(ifxr, iz)%imix = iMix_CrOut
          ENDIF                    
        ENDIF
      ENDDO
      IF(nTracerCntl%lXsLib .AND. nTracerCntl%lCrCspFtn .AND. lCspFtn) THEN
        ifxrbeg = FxrIdxst; ifxrend = FxrIdxSt + nFxrLocal -1
        CALL SetCrCspFtn(Fxr(ifxrbeg:ifxrend, iz), nFxrLocal, CrPos0, CrBank(ibank)%FracCrIn(iz), Pin(ixy)%iasy, iz, nTracerCntl, PE)
      ENDIF
      Pin(ixy)%FracCrIn(iz) = CrBank(ibank)%FracCrIn(iz)
    ENDDO
  ENDDO
ENDDO

CrInfo%lCrChg = .FALSE.
END SUBROUTINE





SUBROUTINE ReplaceMixture(myFxr, imix)
USE TYPEDEF,               ONLY : FxrInfo_Type
USE Material_Mod,          ONLY : Mixture
USE BasicOperation,        ONLY : CP_VA
USE ALLOCS
IMPLICIT NONE
TYPE(FxrInfo_Type) :: myFxr
INTEGER :: imix
INTEGER :: niso

niso = Mixture(imix)%niso
IF(niso .GT. myFxr%ndim) THEN
  DEALLOCATE(myFxr%Pnum, myFxr%idiso)
  CALL Dmalloc(myFxr%pnum, niso + 4); CALL Dmalloc(myFxr%IdIso, niso + 4)
  myFxr%ndim = niso + 4  
ENDIF

myFxr%lMixtureMixing = .FALSE.; myFxr%lCrCspFtn = .FALSE.
myFxr%imix = imix
myFxr%lRes = Mixture(imix)%lRes; myFxr%lh2o = Mixture(imix)%lh2o

myFxr%niso = Mixture(imix)%niso
CALL CP_VA(myFxr%idiso(1:niso), Mixture(imix)%idiso(1:niso), niso)
CALL CP_VA(myFxr%pnum(1:niso), Mixture(imix)%pnum(1:niso), niso)

END SUBROUTINE

SUBROUTINE CuspingMixture(myFxr, iso_crin, iso_crout, FracCrIn)
USE TYPEDEF,          ONLY : FxrInfo_Type,       CspFxr_Type
USE Material_Mod,     ONLY : Mixture
USE BasicOperation,   ONLY : CP_VA
USE ALLOCS
IMPLICIT NONE
TYPE(FxrInfo_Type) :: myFxr
INTEGER :: iso_crin, iso_crout
REAL :: FracCrIn

TYPE(CspFxr_TYPE), POINTER :: CspFxr

INTEGER :: Idiso_new(1024)
REAL :: pnum_new(1024)
INTEGER :: niso_new

REAL :: wt1, wt2
INTEGER :: iso1, iso2, id
LOGICAL :: lh2o_1, lh2o_2
INTEGER :: i, j

iso1 = iso_crin; iso2 = iso_crout
wt1 = FracCrIn; wt2 = 1._8 - wt1
niso_new = 0

niso_new = mixture(iso1)%niso
CALL CP_VA(IdIso_new(1:niso_new), mixture(iso1)%IdIso(1:niso_new), niso_new)
DO i = 1, niso_new
  pnum_new(i) = wt1 * mixture(iso1)%pnum(i)
ENDDO

DO i = 1, Mixture(iso2)%niso
  id = Mixture(iso2)%idiso(i)
  DO j = 1, niso_new
    IF(id .EQ. IdIso_new(j)) EXIT
  ENDDO
  IF(j .EQ. niso_new + 1) THEN
    niso_new = niso_new + 1;     IdIso_new(j) = id
    pnum_new(j) = 0
  ENDIF
  pnum_new(j) = pnum_new(j) + wt2 * Mixture(iso2)%pnum(i)
ENDDO

IF(niso_new .GT. myFxr%ndim) THEN
  DEALLOCATE(myFxr%pnum, myFxr%Idiso)
  CALL Dmalloc(myFxr%pnum, niso_new + 4); CALL Dmalloc(myFxr%IdIso, niso_new + 4)
  myFxr%ndim = niso_new + 4
ENDIF

myFxr%niso = niso_new
DO i = 1, niso_new
  myFxr%idiso(i) = idiso_new(i); myFxr%pnum(i) = pnum_new(i)
ENDDO

lh2o_1 = Mixture(iso1)%lh2o .AND. Mixture(iso2)%lh2o; lh2o_2 = Mixture(iso1)%lh2o .OR. Mixture(iso2)%lh2o
myFxr%lMixtureMixing = .TRUE.

myFxr%lH2o = .FALSE.
IF(lh2o_1 .AND. lh2o_2) THEN
  myFxr%h2ofrac = Mixture(iso1)%h2ofrac0 * wt1 + Mixture(iso2)%h2ofrac0 * wt2
  myFxr%lh2o = .TRUE.
ELSEIF(.NOT. lh2o_1 .AND. lh2o_2) THEN
  IF(Mixture(iso1)%lh2o) myFxr%h2ofrac = Mixture(iso1)%h2ofrac0 * wt1
  IF(Mixture(iso2)%lh2o) myFxr%h2ofrac = Mixture(iso2)%h2ofrac0 * wt2
  myFxr%lh2o = .TRUE.
ENDIF
myFxr%lRes = Mixture(iso1)%lRes .OR. Mixture(iso1)%lRes
myFxr%lMixtureMixing = .TRUE.; myFxr%lCrCspFtn = .FALSE.

IF(.NOT. myFxr%lCrFxr) RETURN
myFxr%lCrCspFtn = .TRUE.
CspFxr => myFxr%CspFxr

!Mixture Map
CspFxr%niso(1) = mixture(iso1)%niso; CspFxr%niso(2) = mixture(iso2)%niso
CspFxr%vwt = (/wt1, wt2/)

CspFxr%lH2O(1) = Mixture(iso1)%lh2o; CspFxr%lH2O(2) = Mixture(iso2)%lh2o
CspFxr%H2oFrac = (/ Mixture(iso1)%h2ofrac0, Mixture(iso2)%h2ofrac0 /)
DO i = 1, mixture(iso1)%niso
  CspFxr%pnum(i, 1) = mixture(iso1)%pnum(i)
  CspFxr%isolist(i, 1) = mixture(iso1)%idiso(i)
  DO j = 1, myFxr%niso
    IF(CspFxr%isolist(i, 1) .EQ. myFXR%idiso(j)) EXIT
  ENDDO
  CspFxr%isomap(i, 1) = j
ENDDO
DO i = 1, mixture(iso2)%niso
  CspFxr%pnum(i, 2) = mixture(iso2)%pnum(i)
  CspFxr%isolist(i, 2) = mixture(iso2)%idiso(i)
  DO j = 1, myFxr%niso
    IF(CspFxr%isolist(i, 2) .EQ. myFXR%idiso(j)) EXIT
  ENDDO
  CspFxr%isomap(i, 2) = j
ENDDO
END SUBROUTINE

SUBROUTINE MixingMixture(myFxr, iso1, iso2, wt0, lComplete)
USE TYPEDEF,          ONLY : FxrInfo_Type
USE Material_Mod,     ONLY : Mixture
USE BasicOperation,   ONLY : CP_VA
USE ALLOCS
IMPLICIT NONE
TYPE(FxrInfo_Type) :: myFxr
INTEGER :: iso1, iso2
REAL :: wt0
LOGICAL :: lComplete

INTEGER :: Idiso_new(1024)
REAL :: pnum_new(1024)
INTEGER :: niso_new

INTEGER :: i, j
INTEGER :: id

REAL :: wtbar

LOGICAL :: lh2o_1, lh2o_2

wtbar = 1._8 - wt0
niso_new = 0
IF(.NOT. lComplete .OR. (wt0 .NE. 1._8)) THEN
  niso_new = mixture(iso1)%niso
  CALL CP_VA(IdIso_new(1:niso_new), mixture(iso1)%IdIso(1:niso_new), niso_new)
  !CALL CP_VA(pnum_new(1:niso_new), myFxr%pnum(1:niso_new), niso_new)
  DO i = 1, niso_new
    pnum_new(i) = wtbar * mixture(iso1)%pnum(i)
  ENDDO
ENDIF

DO i = 1, Mixture(iso2)%niso
  id = Mixture(iso2)%idiso(i)
  DO j = 1, niso_new
    IF(id .EQ. IdIso_new(j)) EXIT
  ENDDO
  IF(j .EQ. niso_new + 1) THEN
    niso_new = niso_new + 1;     IdIso_new(j) = id
    pnum_new(j) = 0
  ENDIF
  pnum_new(j) = pnum_new(j) + wt0 * Mixture(iso2)%pnum(i)
ENDDO

IF(niso_new .GT. myFxr%ndim) THEN
  DEALLOCATE(myFxr%pnum, myFxr%Idiso)
  CALL Dmalloc(myFxr%pnum, niso_new + 4); CALL Dmalloc(myFxr%IdIso, niso_new + 4)
  myFxr%ndim = niso_new + 4
ENDIF

myFxr%niso = niso_new
DO i = 1, niso_new
  myFxr%idiso(i) = idiso_new(i); myFxr%pnum(i) = pnum_new(i)
ENDDO

lh2o_1 = Mixture(iso1)%lh2o .AND. Mixture(iso2)%lh2o; lh2o_2 = Mixture(iso1)%lh2o .OR. Mixture(iso2)%lh2o

myFxr%lMixtureMixing = .TRUE.

IF(lh2o_1 .AND. lh2o_2) THEN
  myFxr%h2ofrac = Mixture(iso1)%h2ofrac0 * (1-wt0) + Mixture(iso2)%h2ofrac0 * wt0
ELSEIF(.NOT. lh2o_1 .AND. lh2o_2) THEN
  IF(Mixture(iso1)%lh2o) myFxr%h2ofrac = Mixture(iso1)%h2ofrac0 * (1-wt0)
  IF(Mixture(iso2)%lh2o) myFxr%h2ofrac = Mixture(iso2)%h2ofrac0 * wt0
  myFxr%lh2o = .TRUE.
ENDIF

myFxr%lRes = Mixture(iso1)%lRes .OR. Mixture(iso1)%lRes

IF(lComplete) THEN
  myFxr%lh2o = Mixture(iso2)%lh2o
  myFxr%imix = iso2
  myFxr%lMixtureMixing = .FALSE.
  myFxr%lRes = Mixture(iso2)%lRes
ENDIF
END SUBROUTINE

SUBROUTINE CalCrBankLocation(FracCrIn, lCrIn, RodInPos, RodPos, hz, nz, lCusping)
REAL :: FracCrIn(nz), hz(nz)
REAL :: RodInPos, RodPos
INTEGER :: nz
LOGICAL :: lCrIn
LOGICAL :: lCusping

REAL :: z(0:nz)
REAL :: CrTip

INTEGER :: iz
REAL :: FlagCrIn
z(0) =0
DO iz = 1, nz
  z(iz) = z(iz-1) + hz(iz)  
ENDDO

CrTip = RodInPos + RodPos
DO iz = nz, 1, -1
  FracCrIn(iz) = (z(iz) - CrTip) / hz(iz)
  FracCrIn(iz) = max(0._8, FracCrIn(iz))
  FracCrIn(iz) = min(1._8, FracCrIn(iz))
  IF(.NOT. lCusping) THEN
    IF(FracCrIn(iz) .LT. 0.5_8) THEN
      FracCrIn(iz) = 0._8
    ELSE
      FracCrIn(iz) = 1._8
    ENDIF
  ENDIF
ENDDO
FlagCrIn = sum(FracCrIn(1:nz))
lCrIn = .FALSE.
IF(FlagCrIn .GT. 1.e-6) lCrIn = .TRUE.
END SUBROUTINE
SUBROUTINE SetCspFxrInfo(Core, FmInfo, GroupInfo, nTracerCntl, PE)
USE TYPEDEF,          ONLY : FmInfo_Type,           FxrInfo_Type,           GroupInfo_Type,        &
                             PE_TYPE,                                                              &
                             CspFxr_Type
USE Material_Mod,      ONLY : Mixture
USE BasicOperation,    ONLY : CP_CA
USE ALLOCS
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTRACERCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE

TYPE(CspFxr_Type), POINTER :: CspFxr
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
INTEGER :: ng, ntiso, nfxr, nfsr, nxy
INTEGER :: myzb, myze
INTEGER :: iz, ifxr

IF(.NOT. nTracerCntl%lCrCsp) RETURN

ntiso = GroupInfo%ntiso
myzb = PE%myzb; myze = PE%myze
nfsr = Core%nCorefsr; nfxr = Core%nCorefxr; ng = GroupInfo%ng

Fxr => FmInfo%Fxr

DO iz = myzb, myze
!DO iz = CrInfo%CrMv_Dom(1), CrInfo%CrMv_Dom(2)
  DO ifxr = 1, nfxr
    IF(.NOT. Fxr(ifxr, iz)%lCrFxr) CYCLE
    ALLOCATE(Fxr(ifxr, iz)%CspFxr)
    CspFxr => Fxr(ifxr, iz)%CspFxr
    CALL DMALLOC(CspFxr%pnum, ntiso, 2)
    CALL DMALLOC(CspFxr%isolist, ntiso, 2)
    CALL DMALLOC(CspFxr%fcsp, ng, 2)
    CALL DMALLOC(CspFxr%isomap, ntiso, 2)
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE SetCrFxrInfo(Core, FmInfo, GroupInfo, nTracerCntl, PE)
USE TYPEDEF,           ONLY : FmInfo_Type,           FxrInfo_Type,           GroupInfo_Type,        &
                              Cell_Type,             Pin_Type,                                                    &
                              PE_TYPE,               Mixture_Type
USE Material_Mod,      ONLY : Mixture
USE BasicOperation,    ONLY : CP_CA
USE XsLib_Mod,         only : nlvflxmax, nreshel
USE ALLOCS
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

INTEGER, PARAMETER :: MaxFxr = 100, MaxFsr = 1000

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Pin_Type), POINTER :: Pin(:)
LOGICAL, POINTER :: AllocFxr(:, :)

INTEGER :: nCatMg, nSubFlx, iResGrpBeg, iResGrpEnd
INTEGER :: myzb, myze, xyb, xye, nCellType   !--- CNJ Edit : Domain Decomposition + MPI

LOGICAL :: lRes1, lRes2
INTEGER :: i, j
INTEGER :: iz, ixy, ifxr, icel, ireg, itype, icrtype, imix1, imix2
INTEGER :: FsrIdxSt, FxrIdxst, nxy 
myzb = PE%myzb; myze = PE%myze
nCellType = Core%nCellType; nxy = Core%nxy

!-- CNJ Edit : Domain Decomposition + MPI
xyb = PE%myPinBeg; xye = PE%myPinEnd
IF (PE%RTMASTER) THEN
  xyb = 1; xye = nxy
ENDIF

iResGrpBeg = GroupInfo%nofg + 1;
iResGrpEnd = GroupInfo%nofg + GroupInfo%norg
Fxr => FmInfo%Fxr; CellInfo => Core%CellInfo; Pin => Core%Pin

ALLOCATE(AllocFxr(MaxFxr, nCellType))
DO itype = 1, nCellType
  AllocFxr(:, itype) = .FALSE.
  IF(.NOT. CellInfo(itype)%lCrCell) CYCLE
  icrtype=CellInfo(itype)%CrCell%CellCrIdx
  DO i = 1, CellInfo(itype)%CrCell%nCrFxr
    j = CellInfo(itype)%CrCell%CrFxrIdx(i)
    !ifxr = FxrIdxst + j - 1
    !Fxr(ifxr, iz)%lCrFxr = .TRUE.
    ireg = CellInfo(itype)%MapFxr2FsrIdx(1, j)
    imix1 = CellInfo(itype)%iReg(ireg); imix2 = CellInfo(icrtype)%iReg(ireg)
    lRes1 = Mixture(imix1)%lres; lRes2 = Mixture(imix2)%lres
    If(lRes2 .AND. .NOT. lRes1) THEN
      AllocFxr(j, itype) = .TRUE.
    ENDIF
  ENDDO
ENDDO

! Why resonance data in FMInfo%Fxr are redundantly allocated ?????
! -> annotated.
!DO iz = myzb, myze
!  DO ixy = xyb, xye
!    FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxst = Pin(ixy)%FxrIdxst 
!    icel = Pin(ixy)%Cell(iz)
!    IF(.NOT. CellInfo(icel)%lCrCell) CYCLE
!    DO j = 1, CellInfo(icel)%nFxr
!      IF(.NOT. AllocFxr(j, icel)) CYCLE
!      
!      ifxr = FxrIdxst + j - 1
!      CALL Dmalloc0(Fxr(ifxr, iz)%XsEq, 1, nlvflxmax, 1, nreshel, iResGrpBeg, iResGrpEnd)   
!      Fxr(ifxr, iz)%XsEq = 0._8
!      CALL Dmalloc0(Fxr(ifxr, iz)%fresoa, iResGrpBeg, iResGrpEnd)
!      CALL Dmalloc0(Fxr(ifxr, iz)%fresoS, iResGrpBeg, iResGrpEnd)
!      CALL Dmalloc0(Fxr(ifxr, iz)%fresostr, iResGrpBeg, iResGrpEnd)
!      CALL Dmalloc0(Fxr(ifxr, iz)%fresoF, iResGrpBeg, iResGrpEnd)
!      CALL CP_CA(Fxr(ifxr, iz)%fresoa(iResGrpBeg:iResGrpEnd), 1._8, iResGrpEnd-iResGrpBeg+1)
!      CALL CP_CA(Fxr(ifxr, iz)%fresoS(iResGrpBeg:iResGrpEnd), 1._8, iResGrpEnd-iResGrpBeg+1)   
!      CALL CP_CA(Fxr(ifxr, iz)%fresoStr(iResGrpBeg:iResGrpEnd), 1._8, iResGrpEnd-iResGrpBeg+1)   
!      CALL CP_CA(Fxr(ifxr, iz)%fresoF(iResGrpBeg:iResGrpEnd), 1._8, iResGrpEnd-iResGrpBeg+1)          
!      
!      ALLOCATE(Fxr(ifxr, iz)%fresoAIso(Fxr(ifxr, iz)%ndim,iResGrpBeg:iResGrpEnd)) 
!      ALLOCATE(Fxr(ifxr, iz)%fresoFIso(Fxr(ifxr, iz)%ndim,iResGrpBeg:iResGrpEnd)) 
!      Fxr(ifxr, iz)%fresoAIso=1._8
!      Fxr(ifxr, iz)%fresoFIso=1._8
!      
!    ENDDO
!  ENDDO
!ENDDO

DO iz = myzb, myze
  DO ixy = xyb, xye
    FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxst = Pin(ixy)%FxrIdxst 
    icel = Pin(ixy)%Cell(iz)
    IF(.NOT. CellInfo(icel)%lCrCell) CYCLE
    icrtype=CellInfo(icel)%CrCell%CellCrIdx
    DO i = 1, CellInfo(icel)%CrCell%nCrFxr
      j = CellInfo(icel)%CrCell%CrFxrIdx(i)
      ifxr = FxrIdxSt + j -1
      Fxr(ifxr, iz)%lCrFxr = .TRUE.
    ENDDO
  ENDDO
ENDDO


DEALLOCATE(AllocFxr)
NULLIFY(CellInfo, Fxr)
END SUBROUTINE

SUBROUTINE SetGlobalCrPinInfo(Core)
USE TYPEDEF,            ONLY : Asy_Type,    Pin_TYPE,          Cell_Type
USE IOUTIL,             ONLY : TERMINATE
USE ALLOCS
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core

TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
INTEGER :: nxy, nxya, nCrAsy, nCrPin
INTEGER :: i, iz, ixy, ixya, ibank, icel

Asy => Core%Asy; Pin => Core%Pin; CellInfo => Core%CellInfo
nxy = Core%nxy; nxya = Core%nxya

nCrAsy = 0; nCrPin = 0
DO ixya = 1, nxya
  IF(.NOT. Core%Asy(ixya)%lCrAsy) CYCLE
  nCrAsy = nCrAsy + 1; nCrPin = nCrPin + Asy(ixya)%nCrPin
  DO i = 1, Asy(ixya)%nCrPin
    ixy = Asy(ixya)%CrPinIdx(i); ibank = Asy(ixya)%CrPinBankid(i)
    IF(Pin(ixy)%lCrPin) CALL TERMINATE('Overlapping of Control Rod')
    Pin(ixy)%lCrPin = .TRUE.
    Pin(ixy)%CrBankId = ibank
    CALL DMALLOC(Pin(ixy)%FracCrIn, Core%nz)
    DO iz = 1, Core%nz
      icel = Pin(ixy)%Cell(iz)
      IF(.NOT. CellInfo(icel)%lCrCell) CYCLE
      IF(CellInfo(icel)%CrCell%lCrIn) THEN
        Pin(ixy)%FracCrIn = 1._8
      ELSE
        Pin(ixy)%FracCrIn = 0._8
      ENDIF
    ENDDO
  ENDDO
ENDDO

CrInfo%nCrAsy = nCrAsy; CrInfo%nCrPin = nCrpin
NULLIFY(Asy, Pin, CellInfo)
END SUBROUTINE

SUBROUTINE SetCrBank(Core)
USE TYPEDEF,          ONLY : Asy_Type,    AsyInfo_Type
USE BasicOperation,   ONLY : CP_VA
USE allocs
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core

INTEGER :: i, j, ix, iy, ixy, iasy, iasytype, icrconf, icrpin
INTEGER :: nxc0, nyc0, nxc, nyc, nxyc, nxya, nCrPin

INTEGER, POINTER :: CrPinLoc(:, :)

TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)

nxc0 = Core%nxc; nyc0 = Core%nyc; nxyc = core%nxyc
nxya = Core%nxya
Asy => Core%Asy; AsyInfo => Core%AsyInfo

ALLOCATE(CrPinLoc(nxc0, nyc0))

DO i = 1, nCrBank
  CALL DMALLOC(CrBank(i)%FracCrIn, Core%nz)
  DO iasy = 1, Core%nxya
    icrconf = CrBank(i)%CrAsyLoc(iasy)
    IF(icrconf .EQ. 0) CYCLE
    iasytype = Asy(iasy)%AsyType
    !icrconf = CrBank(i)%
    nxc = AsyInfo(iasytype)%nx
    nyc = AsyInfo(iasytype)%ny
    IF(nxc .EQ. nxc0 .AND. nyc .EQ. nyc0) THEN
      nCrPin = CrAsyConf(icrconf)%nCrPin
      CALL CP_VA(CrPinLoc(1:nxc0, 1:nyc0), CrAsyConf(icrconf)%CrLoc(1:nxc0, 1:nyc0), nxc0, nyc0)
    ELSE
      CALL SetPartialAsyCrConf(CrPinLoc, nCrpin, nxc, nyc, CrAsyConf(icrconf)%CrLoc, nxc0, nyc0)
    ENDIF 
    IF(nCrPin .EQ. 0) CYCLE
    IF(.NOT. Asy(iasy)%lCrAsy) THEN
      CALL DMALLOC(Asy(iasy)%CrPinIdx, nxyc)
      CALL DMALLOC(Asy(iasy)%CrPinBankId, nxyc)
      Asy(iasy)%lCrAsy = .TRUE.
      iCrPin = 0
    ELSE
      icrPin = Asy(iasy)%nCrPin
    ENDIF
    DO iy = 1, nyc
      DO ix = 1, nxc
        IF(CrPinLoc(ix, iy) .EQ. 0) CYCLE
        iCrPin = iCrPin + 1
        ixy = AsyInfo(iasytype)%Pin2DIdx(ix, iy); ixy = Asy(iasy)%GlobalPinIdx(ixy)
        Asy(iasy)%CrPinIdx(iCrPin) = ixy   !Pin Index
        Asy(iasy)%CrPinBankId(iCrPin) = i  !Bank ID
      ENDDO
    ENDDO
    Asy(iasy)%nCrPin = iCrPin
    CONTINUE
  ENDDO
ENDDO

NULLIFY(Asy, AsyInfo)
DEALLOCATE(CrPinLoc)

END SUBROUTINE

SUBROUTINE SetPartialAsyCrConf(CrLoc, nCrPin, nx, ny, CrLocBase, nx0, ny0)
INTEGER :: CrLoc(nx0,ny0)
INTEGER :: CrLocBase(nx0, ny0)
INTEGER :: nx, ny, nx0, ny0
INTEGER :: nCrPin, ixbeg, iybeg
INTEGER :: ix, iy, ix0, iy0
ixbeg = nx0 - nx + 1
iybeg = ny0 - ny + 1
iy = 0
DO iy0 = iybeg, ny0
  iy = iy + 1; ix = 0
  DO ix0 = ixbeg, nx0
    ix = ix + 1
    CrLoc(ix, iy) = CrLocBase(ix0, iy0)
  ENDDO
ENDDO
nCrPin = 0
DO iy = 1, ny
  DO ix = 1, nx
    IF(CrLoc(ix, iy) .EQ. 0) CYCLE
    nCrPin = nCrPin + 1
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE SetCrAsyConf(Core)
USE ALLOCS
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
INTEGER :: nxc, nyc, ncrPin
INTEGER :: i, ix, iy

nxc = Core%nxc; nyc = Core%nyc

DO i = 1, nCrAsyConf
  ncrPin = 0
  DO iy = 1, nyc
    DO ix = 1, nxc
      IF(CrAsyConf(i)%CrLoc(ix, iy) .EQ. 0) CYCLE
      ncrPin = ncrPin + 1
    ENDDO
  ENDDO
  CrAsyConf(i)%nCrPin = ncrPin
  CALL DMALLOC(CrAsyConf(i)%CrPinIdx, ncrPin)
  
  ncrPin = 0
  DO iy = 1, nyc
    DO ix = 1, nxc
      IF(CrAsyConf(i)%CrLoc(ix, iy) .EQ. 0) CYCLE
      ncrPin = ncrPin + 1
      CrAsyConf(i)%CrPinIdx(ncrPin) = nxc * (iy - 1) + ix
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE SetCrCellInfo(CellInfo, nCellType, lEDGE)
USE TYPEDEF,            ONLY : Cell_Type
USE IOUTIL,             ONLY : TERMINATE
IMPLICIT NONE
TYPE(Cell_Type) :: CellInfo(nCellType)
INTEGER :: nCellType
LOGICAL :: lEdge

INTEGER, PARAMETER :: MaxFxr = 100, MaxFsr = 1000

INTEGER :: i, j, ifsr, ifxr, imix1, imix2
INTEGER :: icrout, icrin, icrout1, icrin1
INTEGER :: nFxr, nFsr, nCrFxr, nCrFsr
INTEGER :: CrFxrIdx(MaxFxr), CrFsrIdx(MaxFsr)
INTEGER :: imix_crin(MaxFxr), imix_crout(MaxFxr)

LOGICAL :: lCent, lCCell

DO i = 1, nCrCell
  icrout=CrCellMap(1, i); icrin=CrCellMap(2, i)
  CellInfo(icrout)%lCrCell = .TRUE.
  CellInfo(icrin)%lCrCell = .TRUE.
  ALLOCATE(CellInfo(icrout)%CrCell)
  ALLOCATE(CellInfo(icrin)%CrCell)
  CellInfo(icrout)%CrCell%CellCrIdx = icrin
  CellInfo(icrin)%CrCell%CellCrIdx = icrout
  CellInfo(icrout)%CrCell%lCrout = .TRUE.;CellInfo(icrout)%CrCell%lCrIn = .FALSE.
  CellInfo(icrin)%CrCell%lCrout = .TRUE.;CellInfo(icrin)%CrCell%lCrIn = .TRUE.
  
  IF(lEDGE) CYCLE
  IF(Cellinfo(iCrout)%lCCell) CYCLE 
  DO j = 1, 3
    icrout1 = CellInfo(icrout)%EdgeCellIdx(j)
    icrin1 = CellInfo(icrin)%EdgeCellIdx(j)
    
    CellInfo(icrout1)%lCrCell = .TRUE.
    CellInfo(icrin1)%lCrCell = .TRUE.
    ALLOCATE(CellInfo(icrout1)%CrCell)
    ALLOCATE(CellInfo(icrin1)%CrCell)
    CellInfo(icrout1)%CrCell%CellCrIdx = icrin1
    CellInfo(icrin1)%CrCell%CellCrIdx = icrout1
    CellInfo(icrout1)%CrCell%lCrout = .TRUE.;CellInfo(icrout1)%CrCell%lCrIn = .FALSE.
    CellInfo(icrin1)%CrCell%lCrout = .TRUE.;CellInfo(icrin1)%CrCell%lCrIn = .TRUE.    
  ENDDO
ENDDO
DO i = 1, nCellType
  IF(.NOT. CellInfo(i)%lCrCell) CYCLE
  IF(CellInfo(i)%CrCell%LCrIn) CYCLE
  lCent = CellInfo(i)%lCentX .or. CellInfo(i)%lCentY .or. CellInfo(i)%lCentXY
  lCCell = CellInfo(i)%lCCell
  IF(lCent .AND. lCCell) CYCLE
  icrout = i; icrin = CellInfo(i)%CrCell%CellCrIdx
  IF(CellInfo(icrout)%nFxr .NE. CellInfo(icrin)%nFxr) CALL TERMINATE('Wrong Control Rod Cell Assignment')
  IF(CellInfo(icrout)%nFsr .NE. CellInfo(icrin)%nFsr) CALL TERMINATE('Wrong Control Rod Cell Assignment')
  nFxr = CellInfo(icrout)%nFxr; nFsr = CellInfo(icrout)%nFsr
  nCrFxr = 0; nCrFsr = 0
  DO ifxr = 1, nFXR
    IF(CellInfo(icrout)%nFsrInFxr(ifxr) .NE. CellInfo(icrin)%nFsrInFxr(ifxr)) THEN
      CALL TERMINATE('Wrong Control Rod Cell Assignment')
    ENDIF
    iFsr = CellInfo(icrout)%MapFxr2FsrIdx(1, ifxr)
    imix1 = CellInfo(icrout)%iReg(ifsr)
    imix2 = CellInfo(icrin)%iReg(ifsr)
    IF(imix1 .EQ. imix2) CYCLE
    nCrFxr = nCrFxr + 1
    CrFxrIdx(nCrFxr) = iFxr;
    imix_crout(nCrFxr) = imix1; imix_crin(nCrFxr) = imix2; 
    DO j = 1, CellInfo(icrout)%nFsrInFxr(ifxr)
      iFsr =  CellInfo(icrout)%MapFxr2FsrIdx(j, ifxr)
      IF(CellInfo(icrin)%MapFxr2FsrIdx(j, ifxr) .NE. iFsr) CALL TERMINATE('Wrong Control Rod Cell Assignment')
      nCrFsr = nCrFsr + 1
      CrFsrIdx(nCrFsr) = iFsr
    ENDDO
  ENDDO
  CONTINUE
  CALL DMALLOC(CellInfo(icrin)%CrCell%CrFxrIdx, nCrFxr)
  CALL DMALLOC(CellInfo(icrin)%CrCell%CrFsrIdx, nCrFsr)
  CALL DMALLOC(CellInfo(icrin)%CrCell%imix_crin, nCrFxr)
  CALL DMALLOC(CellInfo(icrin)%CrCell%imix_crout, nCrFxr)
  
  CALL DMALLOC(CellInfo(icrout)%CrCell%CrFxrIdx, nCrFxr)
  CALL DMALLOC(CellInfo(icrout)%CrCell%CrFsrIdx, nCrFsr)
  CALL DMALLOC(CellInfo(icrout)%CrCell%imix_crin, nCrFxr)
  CALL DMALLOC(CellInfo(icrout)%CrCell%imix_crout, nCrFxr)
  
  CellInfo(icrin)%CrCell%nCrFxr = nCrFxr; CellInfo(icrin)%CrCell%nCrFsr = nCrFsr  
  CellInfo(icrin)%CrCell%CrFxrIdx(1:nCrFxr) = CrFxrIdx(1:nCrFxr)
  CellInfo(icrin)%CrCell%CrFsrIdx(1:nCrFsr) = CrFxrIdx(1:nCrFsr)
  CellInfo(icrin)%CrCell%imix_crin(1:nCrFxr) = imix_crin(1:nCrFxr)
  CellInfo(icrin)%CrCell%imix_crout(1:nCrFxr) = imix_crout(1:nCrFxr)

  CellInfo(icrOut)%CrCell%nCrFxr = nCrFxr; CellInfo(icrOut)%CrCell%nCrFsr = nCrFsr
  CellInfo(icrout)%CrCell%CrFxrIdx(1:nCrFxr) = CrFxrIdx(1:nCrFxr)
  CellInfo(icrout)%CrCell%CrFsrIdx(1:nCrFsr) = CrFxrIdx(1:nCrFsr)
  CellInfo(icrout)%CrCell%imix_crin(1:nCrFxr) = imix_crin(1:nCrFxr)
  CellInfo(icrout)%CrCell%imix_crout(1:nCrFxr) = imix_crout(1:nCrFxr)
ENDDO

END SUBROUTINE

SUBROUTINE CalCrFluxDip(io, Core, BankId, iz, FmInfo, CmInfo, GroupInfo, nTracerCntl, PE)
USE TYPEDEF,              ONLY : FmInfo_Type,        CmInfo_Type,      GroupInfo_Type,      &
                                 Pin_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
INTEGER :: BankId, iz, io

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
REAL, POINTER :: Phis(:, :, :)
REAL :: FluxDip(500)

INTEGER :: myzb, myze, xyb, xye, nxy, ng   !--- CNJ Edit : Domain Decomposition + MPI
INTEGER :: ig, ixy, icel, ifxr, ifsr
INTEGER :: FxrIdxSt, FsrIdxSt, nFxrLocal, nFsrInFxr, nCrFxr
INTEGER :: i, j, k, l
REAL :: Fluxsum1(500), FluxSum2, Volsum1(500), volsum2
ng = GroupInfo%ng
myzb = PE%myzb; myze = PE%myze; nxy = Core%nxy

!-- CNJ Edit : Domain Decomposition + MPI
xyb = PE%myPinBeg; xye = PE%myPinEnd
IF (PE%RTMASTER) THEN
  xyb = 1; xye = nxy
ENDIF

Pin => Core%Pin; CellInfo => Core%CellInfo
Phis => FmInfo%Phis

IF( iz .GT. myze .OR. iz .LT. myzb) RETURN
DO ig = 1, ng
  Fluxsum1= 0; FluxSum2=0
  volsum1=0; volsum2=0
  nCRFXR = 0
  DO ixy = xyb, xye
    IF(.NOT. Pin(ixy)%lCrPin) CYCLE
    IF(Pin(ixy)%CrBankId .NE. Bankid) CYCLE
    icel = Pin(ixy)%Cell(iz)
    IF(.NOT. CellInfo(icel)%lCrCell) CYCLE
    nFxrLocal = CellInfo(icel)%nFxr; FxrIdxSt = Pin(ixy)%FxrIdxSt
    FsrIdxSt = Pin(ixy)%FsrIdxSt
    nCrFxr = CellInfo(icel)%CrCell%nCrFxr
    
    DO i = 1, nCrFxr 
      j = CellInfo(icel)%CrCell%CrFxrIdx(i)
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
      ifxr = FxrIdxSt + j - 1        
      DO k = 1, nFsrInFxr
        l = CellInfo(icel)%MapFxr2FsrIdx(k, j)
        ifsr = FsrIdxSt + l - 1
        Fluxsum1(i) = FluxSum1(i) +  Phis(ifsr, iz, ig) * CellInfo(icel)%vol(l)
        volsum1(i) = volsum1(i) + CellInfo(icel)%vol(l)
      ENDDO
    ENDDO
    
    DO j = 1, nFxrLocal
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
      ifxr = FxrIdxSt + j - 1        
      DO k = 1, nFsrInFxr
        l = CellInfo(icel)%MapFxr2FsrIdx(k, j)
        ifsr = FsrIdxSt + l - 1
        Fluxsum2 = FluxSum2 +  Phis(ifsr, iz, ig) * CellInfo(icel)%vol(l)
        volsum2 = volsum2 + CellInfo(icel)%vol(l)
      ENDDO
    ENDDO
  ENDDO
  DO i = 1, nCrFxr
    Fluxsum1(i) = Fluxsum1(i) / volsum1(i)
  ENDDO
  fluxsum2= fluxsum2 / volsum2
  IF(ig .EQ. 1) WRITE(io) nCRFXR
  WRITE(io) ig, (fluxsum1(i)/fluxsum2, i = 1, nCrFxr)
  !fluxsum2 = fluxsum2 / volsum2
  !fluxdip(ig) = fluxsum1 / fluxsum2
ENDDO
!WRITE(999, '(5e20.7)') (fluxdip(ig), ig = 1, ng)
END SUBROUTINE


FUNCTION GetBankName(BankId)
IMPLICIT NONE
CHARACTER(10) GetBankName
INTEGER :: BankID

IF(BankID .GT. nCrBank ) THEN
  GetBankName=''
  RETURN
ENDIF
GetBankName = CrBank(BankID)%BankName
END FUNCTION


SUBROUTINE CntlRodSearch(Core, FmInfo, CmInfo, eigv, GroupInfo, GcGroupInfo, nTracerCntl, PE)
USE TYPEDEF,          ONLY : CoreInfo_Type,       FmInfo_Type,       CmInfo_Type,   &
                             GroupInfo_Type,      PE_Type
USE CNTL,             ONLY : nTracerCntl_Type
USE itrcntl_mod,      ONLY : ItrCntl_Type
USE IOUTIL,            ONLY : message
USE FILES,            ONLY : io8
IMPLICIT NONE 
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(GroupInfo_Type) :: GcGroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl

TYPE(PE_Type) :: PE

TYPE(ItrCntl_TYpe), SAVE :: WorthItrCntl
REAL :: eigv

REAL :: eigv_est, deigv
LOGICAL :: lCrUpdt

IF(CrInfo%CrSearch%lSearchConv) RETURN
IF(PE%MASTER) THEN
  WRITE(mesg,'(a)') 'Control Rod Position Search ...'
  CALL message(io8, TRUE, TRUE, mesg)      
ENDIF
eigv_est = eigv

IF(.NOT. nTracerCntl%lGcCMFD) THEN
  WorthItrCntl%GcCmfdItrCntl%conv_rescrit=0.01
  WorthItrCntl%GcCmfdItrCntl%conv_eigvcrit=1.e-5
  WorthItrCntl%GcCmfdItrCntl%lLogOut = .FALSE.
  IF(PE%MASTER) THEN
    WRITE(mesg,'(a)') 'Performing CGCMFD for CR Worth Estimation ...'
    CALL message(io8, TRUE, TRUE, mesg)      
  ENDIF
  CALL GcCmfdAcc(Core, CmInfo, Eigv_est, .FALSE., GroupInfo, GcGroupInfo, nTracerCntl, WorthItrCntl, PE)   
  IF(PE%MASTER) THEN
    WRITE(mesg,'(a20, F10.6)') 'Estimated k-eff ', eigv_est
    CALL message(io8, TRUE, TRUE, mesg)      
  ENDIF  
ENDIF
deigv = abs(eigv - nTracerCntl%target_eigv)

IF(CrInfo%CrSearch%iter .EQ. 0) THEN
  lCrUpdt =.TRUE.
ELSE
  CALL CalCrWorth(Eigv_est, lCrUpdt, PE)
ENDIF

IF(deigv < 0.00010_8 .AND. lCrUpdt) THEN
  CrInfo%CrSearch%lSearchConv = .TRUE.  
  
ENDIF

IF(lCrUpdt) THEN
  CALL UpdtCrSearchPosition(nTracerCntl%Target_eigv, eigv_est, .FALSE., PE%MASTER)
  CALL SetCrBankPosition(Core, FmInfo, CmInfo, GroupInfo, nTracerCntl, PE)
ENDIF

END SUBROUTINE

SUBROUTINE CalCrWorth(Eigv, lCrUpdt, PE)
USE IOUTIL,         ONLY : message
USE FILES,          ONLY : io8
IMPLICIT NONE
REAL :: Eigv
LOGICAL :: lCrUpdt
TYPE(PE_TYPE) :: PE
REAL :: CrMv
REAL :: Worthd, Worth
REAL :: reigvd, reigv
REAL :: rdif

CrInfo%CrSearch%eigv_iter = CrInfo%CrSearch%eigv_iter + 1
Worthd = CrInfo%CrSearch%CrWorth
CrMv = CrInfo%CrSearch%CrMv
reigv = 1._8 / eigv; reigvd = 1._8 / CrInfo%CrSearch%eigvd
Worth = 1.0_8 * (-reigvd + reigv) * 1.e+5 / (CrMv)

lCrUpdt = .FALSE.
rdif = abs((Worth-Worthd) / Worth )
!PRINT *, rdif, Worth, Worthd
IF(rdif < 0.10_8) lCrUpdt = .TRUE.
IF(CrInfo%CrSearch%eigv_iter .GT. 3) lCrUpdt = .TRUE.
IF(lCrUpdt) CrInfo%CrSearch%eigv_iter = 0
IF(PE%MASTER) THEN
  WRITE(mesg,'(11x, a22, F9.4, A7, F9.2, A, L3)') 'Estimated CR Worth :', Worth, 'pcm/cm -', rdif*100, '%', lCrUpdt
  CALL message(io8, FALSE, TRUE, mesg)        
ENDIF

CrInfo%CrSearch%CrWorth = Worth
END SUBROUTINE

SUBROUTINE UpdtCrSearchPosition(Target_eigv, eigv, lreset, MASTER)
IMPLICIT NONE
REAL :: Target_eigv, eigv
LOGICAL :: lReset, MASTER

REAL :: reigv, reigvd
REAL :: Worth, Worth0, CrMvd, CrMv


IF(.NOT. CrInfo%CrSearch%linit) CALL InitSearchCrPosition(CrInfo%CrSearch, MASTER)

CrMv = CrInfo%CrSearch%CrMv
CrMvd =  CrInfo%CrSearch%CrMv
reigv = 1._8 / eigv; reigvd = 1._8 / CrInfo%CrSearch%eigvd
CrInfo%CrSearch%iter = CrInfo%CrSearch%iter + 1
Worth = CrInfo%CrSearch%CrWorth
IF(CrInfo%CrSearch%iter .GT. 1) THEN
  !PRINT *, 'eigv:',reigvd, 'eigv:',reigv, 'crmvd:', CrMv
  Worth = 1.0_8 * (-reigvd + reigv) * 1.e+5 / (CrMv)
  Worth0 = worth
  Worth = min(Worth, 100.0_8)
  Worth = max(Worth, 1.0_8)
ENDIF
CrMv = (1._8 / target_eigv - reigv) * 1.e+5_8 / Worth
CrMv = min(40._8, CrMv)
CrMv = max(-40._8, CrMv)
CALL ChangeCrSearchBankPos(CrInfo%CrSearch, CrMv)

CrInfo%CrSearch%CrMv = CrMv
CrInfo%CrSearch%CrMvd = CrInfo%CrSearch%CrMv
CrInfo%CrSearch%eigvd = eigv
CrInfo%CrSearch%CrWorthd = CrInfo%CrSearch%CrWorth
CrInfo%CrSearch%CrWorth = Worth
END SUBROUTINE

SUBROUTINE ChangeCrSearchBankPos(CrSearch, CrMv)
IMPLICIT NONE
TYPE(CrSearch_Type) :: CrSearch
REAL :: CrMv


REAL :: RodPos(nCrBankMax), RodPos0(nCrBankMax)
LOGICAL :: lCrSearch(nCrBankMax)
REAL :: MaxPos, MaxPos0, Del
INTEGER :: ibank

CrSearch%lFullIn = .FALSE.; CrSearch%lFullOut = .FALSE.

MaxPos = 10000._8
MaxPos0 = 10000._8
DO ibank = 1, nCrBank
  RodPos(ibank) = CrBank(ibank)%RodPos
  RodPos0(ibank) = RodPos(ibank) + CrBank(ibank)%RodInPos
  lCrSearch(ibank) = CrBank(ibank)%lCrSearch
ENDDO

DO ibank = 1, nCrBank
  IF(.NOT. lCrSearch(ibank)) CYCLE
  RodPos(ibank) = RodPos(ibank) - CrMv
  RodPos0(ibank) = RodPos0(ibank) - CrMv
  MaxPos = MIN(MaxPos, RodPos(ibank))
  MaxPos0 = MIN(MaxPos0, RodPos0(ibank))
ENDDO

Del = 0 
IF(MaxPos .LT. 0._8) THEN
  Del = MaxPos
  CrSearch%lFullIn = .TRUE.
ELSEIF(MaxPos0 .GT. hzmax) THEN
  Del = hzmax - MaxPos0
  CrSearch%lFullOut = .TRUE.
ENDIF

DO ibank = 1, nCrBank
  IF(.NOT. lCrSearch(ibank)) CYCLE
  CrBank(ibank)%RodPos = RodPos(ibank) + Del
  PRINT *, CrBank(ibank)%RodPos
ENDDO
!CrMv = CrMv  + Del
END SUBROUTINE

SUBROUTINE InitSearchCrPosition(CrSearch, MASTER)
IMPLICIT NONE
TYPE(CrSearch_Type) :: CrSearch
LOGICAL :: MASTER

INTEGER :: iBank
LOGICAL :: lCrIn
INTEGER :: MaxBankIdx
REAL :: MaxPos, Pos, Del

lCrIn = .FALSE.
MaxPos = 100000._8
DO ibank = 1, nCrBank
  IF(.NOT. CrBank(ibank)%lCrSearch) CYCLE
  Pos = CrBank(ibank)%RodInPos + CrBank(ibank)%RodPos
  IF(CrBank(ibank)%lCrIn)  lCrIn = .TRUE.
  IF(Pos .LT. MaxPos) THEN
    MaxPos = Pos
    MaxBankIdx = ibank
  ENDIF
ENDDO

IF(.NOT. lCrIn) THEN
  Del = abs(MaxPos-hzMax)
  DO ibank = 1, nCrBank
    IF(.NOT. CrBank(ibank)%lCrSearch) CYCLE
    CrBank(ibank)%RodPos = CrBank(ibank)%RodPos - Del
  ENDDO
ENDIF

CrSearch%lInit = .TRUE.
END SUBROUTINE

END MODULE
