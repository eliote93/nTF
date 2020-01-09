#include <defines.h>
SUBROUTINE DcplSseig()
USE PARAM

USE GEOM,            ONLY : Core, ng
USE RAYS,            ONLY : RayInfo
USE Core_mod,        ONLY : GroupInfo        ,CmInfo            ,THInfo                 &
                           ,Eigv             ,FmInfo            ,GcGroupInfo
USE DcplCore_Mod,    ONLY : DcplInfo         ,DcplFmInfo        ,DcplCmInfo             &
                           ,DcplThInfo       ,DcplFxr
USE PE_MOD,          ONLY : PE               ,DcplPE
USE CNTL,            ONLY : nTracerCntl     ,DcplControl
USE SubGrp_mod,      ONLY : SubGrpFsp
USE DcplXsGen_Mod,   ONLY : RefPlnXsGeneration, DcplConvChk

USE Th_Mod,          ONLY : SteadyStateTH, THFeedBack
USE ItrCNTL_mod,     ONLY : ItrCntl, DcplItrCntl
USE FILES,           ONLY : io8
USE IOUTIL,          ONLY : message
IMPLICIT NONE
INTEGER :: iRefPln, iRefTemp
INTEGER :: nRefpln, nRefTemp
INTEGER :: iter
LOGICAL :: lSubGrp, lPlnMocConv
LOGICAL :: Master

IF(nTracerCntl%lFeedBack) THEN
  CALL DcplSseigTH()
  RETURN
ENDIF



nRefTemp = DcplInfo%nRefTemp; nRefPln = DcplInfo%nRefPln
DcplInfo%DcplFmInfo => DcplFmInfo; DcplInfo%DcplCmInfo => DcplCmInfo

Master = PE%Master

WRITE(mesg, '(A)') hbar1(1:77)
IF(Master) CALL message(io8, FALSE, TRUE, MESG)

DO iter = 1, 2
  lSubGrp = FALSE
  IF(iter .EQ. 1) lSubGrp = FALSE
  CALL RefPlnXsGeneration(Core, DcplInfo, DcplThInfo, GroupInfo, DcplControl, DcplItrCntl, PE, DcplPE, lSubGrp, .TRUE.)
  CALL DcplCmfd3d(Core, CmInfo, FmInfo, DcplInfo, Eigv, ng, GroupInfo, GcGroupInfo, nTracerCntl, ItrCntl, PE)
ENDDO
!CALL RefPlnXsGeneration(Core, DcplInfo, DcplThInfo, GroupInfo, DcplControl, DcplItrCntl, DcplPE, FALSE)
!CALL DcplCmfd3d(Core, CmInfo, DcplInfo, Eigv, ng, GroupInfo, nTracerCntl, ItrCntl, PE)

WRITE(mesg, '(A)') hbar2(1:77)
IF(Master) CALL message(io8, FALSE, TRUE, MESG)

Write(mesg, '(A27, 5X, F8.5)')  'k-eff =', Eigv
IF(Master) CALL message(io8, FALSE, TRUE, MESG)
END SUBROUTINE

SUBROUTINE DcplSseig_woThFtn()
USE PARAM

USE GEOM,            ONLY : Core, ng
USE RAYS,            ONLY : RayInfo
USE Core_mod,        ONLY : GroupInfo        ,CmInfo            ,FmInfo                 &
                           ,THInfo           ,Eigv              ,GcGroupInfo
USE DcplCore_Mod,    ONLY : DcplInfo         ,DcplFmInfo        ,DcplCmInfo             &
                           ,DcplThInfo       ,DcplFxr
USE PE_MOD,          ONLY : PE               ,DcplPE
USE CNTL,            ONLY : nTracerCntl      ,DcplControl
USE SubGrp_mod,      ONLY : SubGrpFsp
USE DcplXsGen_Mod,   ONLY : RefPlnXsGeneration ,DcplConvChk

USE DcplTH_Mod,      ONLY : DcplThUpdate
USE ItrCNTL_mod,     ONLY : ItrCntl, DcplItrCntl
USE FILES,           ONLY : io8
USE IOUTIL,          ONLY : message
USE MPIDcpl_Mod,     ONLY : MPI_DcplMessage
IMPLICIT NONE
INTEGER :: iRefPln, iRefTemp
INTEGER :: nRefpln, nRefTemp
INTEGER :: iter
INTEGER :: myrank
LOGICAL :: lSubGrp, lPlnMocConv, lParallel


nRefTemp = DcplInfo%nRefTemp; nRefPln = DcplInfo%nRefPln
DcplInfo%DcplFmInfo => DcplFmInfo; DcplInfo%DcplCmInfo => DcplCmInfo
DcplInfo%Fxr => FmInfo%Fxr
!Parallel Variables
lParallel = PE%lDcplParallel; myrank = PE%myCmfdrank
lParallel = TRUE
!Initial TH-Sweep
IF(nTracerCntl%lFeedBack) THEN
  CALL DcplThUpdate(Core, CmInfo, FmInfo, ThInfo, DcplInfo, DcplThInfo, Eigv, ng, TRUE, GroupInfo, nTracerCntl, DcplItrCntl, PE)
  CALL DcplCmfd3d(Core, CmInfo, FmInfo, DcplInfo, Eigv, ng, GroupInfo, GcGroupInfo, nTracerCntl, ItrCntl, PE)
  !CALL DcplThUpdate(Core, CmInfo, FmInfo, ThInfo, DcplInfo, Eigv, ng, nTracerCntl, PE)
ENDIF
nTracerCntl%lDcplCmfdXsGen = .FALSE.

WRITE(mesg, '(A)') hbar1(1:77); 
IF(.NOT. lParallel) CALL message(io8, FALSE, TRUE, MESG)
IF(lParallel) CALL MPI_DcplMessage(mesg, myrank, FALSE, TRUE, FALSE)
DO iter = 1, 20
  lSubGrp = FALSE
  IF(iter .EQ. 1) lSubGrp = TRUE
  IF(nTracerCntl%lFeedBack) THEN
    CALL DcplThUpdate(Core, CmInfo, FmInfo, ThInfo, DcplInfo, DcplThInfo, Eigv, ng, TRUE, GroupInfo, nTracerCntl, DcplItrCntl, PE)
  ENDIF
  CALL RefPlnXsGeneration(Core, DcplInfo, DcplThInfo, GroupInfo, DcplControl, DcplItrCntl, PE, DcplPE, lSubGrp, .TRUE.)
  lPlnMocConv =  DcplConvChk(DcplItrCntl(1:nRefPln), nRefPln)
  PRINT *, '===Convgence Check :', lPlnMocConv  
  CALL DcplCmfd3d(Core, CmInfo, FmInfo, DcplInfo, Eigv, ng, GroupInfo, GcGroupInfo, nTracerCntl, ItrCntl, PE)
  PAUSE
ENDDO
!CALL RefPlnXsGeneration(Core, DcplInfo, DcplThInfo, GroupInfo, DcplControl, DcplItrCntl, DcplPE, FALSE)
!CALL DcplCmfd3d(Core, CmInfo, DcplInfo, Eigv, ng, GroupInfo, nTracerCntl, ItrCntl, PE)

WRITE(mesg, '(A)') hbar2(1:77)
CALL message(io8, FALSE, TRUE, MESG)

Write(mesg, '(A27, 5X, F8.5)')  'k-eff =', Eigv
CALL message(io8, FALSE, TRUE, MESG)
END SUBROUTINE

SUBROUTINE DcplSseigTH()
USE PARAM

USE GEOM,            ONLY : Core, ng
USE RAYS,            ONLY : RayInfo
USE Core_mod,        ONLY : GroupInfo        ,CmInfo            ,FmInfo                 &
                           ,THInfo           ,Eigv              ,GcGroupInfo
USE DcplCore_Mod,    ONLY : DcplInfo         ,DcplFmInfo        ,DcplCmInfo             &
                           ,DcplThInfo       ,DcplFxr
USE PE_MOD,          ONLY : PE               ,DcplPE
USE CNTL,            ONLY : nTracerCntl      ,DcplControl
USE SubGrp_mod,      ONLY : SubGrpFsp
USE DcplXsGen_Mod,   ONLY : RefPlnXsGeneration ,DcplConvChk     ,ThXsGeneration
USE DcplTH_Mod,      ONLY : DcplThUpdate
USE ItrCNTL_mod,     ONLY : ItrCntl, DcplItrCntl
USE FILES,           ONLY : io8
USE IOUTIL,          ONLY : message
IMPLICIT NONE
INTEGER :: iRefPln, iRefTemp
INTEGER :: nRefpln, nRefTemp
INTEGER :: iter, iterTH
LOGICAL :: lSubGrp, lPlnMocConv, lXsFtn, lUpdt
LOGICAL :: lThConv, lMOC, lQuickUpdt
LOGICAL :: Master

nRefTemp = DcplInfo%nRefTemp; nRefPln = DcplInfo%nRefPln
DcplInfo%DcplFmInfo => DcplFmInfo; DcplInfo%DcplCmInfo => DcplCmInfo

Master = PE%Master

lXsFtn = nTracerCntl%lXsFtn

!Initial TH-Sweep
IF(nTracerCntl%lFeedBack) THEN
  CALL DcplThUpdate(Core, CmInfo, FmInfo, ThInfo, DcplInfo, DcplThInfo, Eigv, ng, TRUE, GroupInfo, nTracerCntl, DcplItrCntl, PE)
  CALL DcplCmfd3d(Core, CmInfo, FmInfo, DcplInfo, Eigv, ng, GroupInfo, GcGroupInfo, nTracerCntl, ItrCntl, PE)
ENDIF

iterTH = 0
!iterquickTH = 0
lThConv = FALSE
ItrCntl%nThCondiGen = 2
lMOC = .TRUE.
!#define mod2
DO iter = 1, 12
  IF(lXsftn) THEN
    lUpdt = FALSE; lQuickUpdt = FALSE
    IF(iter .GT. 2 .AND. mod(iter, 3) .EQ. 0) lUpdt = .TRUE.
    IF(iterTH .GE. ItrCntl%nThCondiGen) lUpdt = .FALSE.  
#ifdef mod1   
    IF(iterTH .GE. 2) THEN
      lMOC =.FALSE. 
      lQuickUpdt = .TRUE.
    ENDIF
#endif
#ifdef mod2
    IF(iterTH .GE. 1 .AND. .NOT. lUpdt) THEN
        lMOC =.FALSE. 
        lQuickUpdt = .TRUE.      
    ENDIF
#endif
    
    CALL DcplThUpdate(Core, CmInfo, FmInfo, ThInfo, DcplInfo, DcplThInfo, Eigv, ng, lUpdt, GroupInfo, nTracerCntl, DcplItrCntl, PE)
    
    IF(lUpdt) THEN
      WRITE(mesg, '(A)') hbar1(1:77); IF(Master) CALL message(io8, FALSE, TRUE, MESG)  
      iterTH = iterTH + 1
      lSubGrp = FALSE; IF(iterTH .EQ. 2) lSubGrp = TRUE
      CALL ThXsGeneration(Core, DcplInfo, DcplThInfo, GroupInfo, DcplControl, DcplItrCntl, PE, DcplPE, lSubGrp, lMOC)
      nTracerCntl%lDcplCmfdXsGen = .FALSE. 
      WRITE(mesg, '(A)') hbar1(1:77); IF(Master)  CALL message(io8, FALSE, TRUE, MESG)
    ENDIF
    IF(lQuickUpdt) THEN
      lSubGrp = .FALSE.
      CALL ThXsGeneration(Core, DcplInfo, DcplThInfo, GroupInfo, DcplControl, DcplItrCntl, PE, DcplPE, lSubGrp, .FALSE.)
    ENDIF
  ELSEIF(.NOT. lXsFtn) THEN
    lSubGrp = FALSE; IF(iter .EQ. 1) lsubGrp = FALSE
    lUpdt = TRUE
    CALL DcplThUpdate(Core, CmInfo, FmInfo, ThInfo, DcplInfo, DcplThInfo, Eigv, ng, TRUE, GroupInfo, nTracerCntl, DcplItrCntl, PE)
    CALL RefPlnXsGeneration(Core, DcplInfo, DcplThInfo, GroupInfo, DcplControl, DcplItrCntl, PE, DcplPE, lSubGrp, lMOC)
    nTracerCntl%lDcplCmfdXsGen = .FALSE. 
    !lPlnMocConv =  DcplConvChk(DcplItrCntl(1:nRefPln), nRefPln)
  ENDIF
  IF(ThInfo%TdopChg .LT. epsm4) lThConv = TRUE
  CALL DcplCmfd3d(Core, CmInfo, FmInfo, DcplInfo, Eigv, ng, GroupInfo, GcGroupInfo, nTracerCntl, ItrCntl, PE)
  IF(iter .LT. 9) CYCLE
  IF(lThConv .AND. ItrCntl%lConv) THEN
    IF(lXsFtn .AND. iterTH .GE. ItrCntl%nThCondiGen) EXIT
    IF(.NOT. lXsFtn) EXIT
  ENDIF
ENDDO

WRITE(mesg, '(A)') hbar2(1:77)
IF(Master)  CALL message(io8, FALSE, TRUE, MESG)

Write(mesg, '(A27, 5X, F8.5)')  'k-eff =', Eigv
IF(Master) CALL message(io8, FALSE, TRUE, MESG)
END SUBROUTINE
