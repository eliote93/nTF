#include <defines.h>
SUBROUTINE ThXsGeneration(Core, DcplInfo, THInfo, GroupInfo, DcplCntl, DcplItrCntl, PE, DcplPE, lSubGrp0, lMOC)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type        ,DcplInfo_Type        ,GroupInfo_Type    &
                          ,ThInfo_Type          ,PE_Type                                 &
                          ,FmInfo_Type          ,CmInfo_Type          ,FxrInfo_Type
USE CNTL,            ONLY : nTracerCntl_Type
USE itrcntl_mod,     ONLY : ItrCntl_TYPE                                                &
                           ,DcplXsGenCntl_Type
USE RAYS,            ONLY : RayInfo
USE DcplTh_Mod,      ONLY : SetThXsGenCondition ,SaveXsGenInfo
USE DcplXsGen_Mod,   ONLY : DcplPlnMOC
USE DcplMicXsFtn_mod,ONLY : UpdtMicroXsFtn
USE SubGrp_mod,      ONLY : SubGrpFsp                           
USE FILES,           ONLY : io8
USE IOUTIL,          ONLY : message
USE BasicOperation,  ONLY : CP_VA
#ifdef MPI_ENV
USE MPIDcpl_Mod,     ONLY : Idle_DcplPlnMOC     ,MPI_DcplMessage
USE MPIComm_mod,     ONLY : MPI_SeqMessage
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(DcplInfo_Type) :: DcplInfo
TYPE(THInfo_Type), POINTER :: THInfo(:)
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: DcplCntl
TYPE(ItrCntl_Type) :: DcplItrCntl(100)
TYPE(PE_TYPE) :: PE
TYPE(PE_TYPE) :: DcplPE(:)
LOGICAL :: lSubGrp0, lMOC

TYPE(CmInfo_Type), POINTER :: DcplCmInfo(:, :)
TYPE(FmInfo_Type), POINTER :: DcplFmInfo(:, :)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)

INTEGER :: iRefPln, iRefTemp
INTEGER :: nRefpln, nRefTemp, ng
INTEGER :: myRefPlnBeg, myRefPlnEnd
INTEGER :: nfsr, nxy
LOGICAL :: lSubGrp, lXsLib, lFuelPln, lMPI

INTEGER :: i, j, k
INTEGER :: ipln, imode
INTEGER :: nXsGenFbIterMax(3)
#ifdef MPI_ENV
INTEGER :: myrank, nproc, comm
LOGICAL :: Master 
#endif
REAL :: Frac(5)

!DATA Frac /0.0_8, 0.2_8, 0,2_8, 0.8_8, 0.8_8/
!DATA FRAC /0._8, 2*0.3_8, 2*0.8_8/
DATA Frac / 0._8, 0.3_8, 0.8_8, 2*0._8/
lXslib = DcplCntl%lXsLib; 
lSubGrp = lXsLib .AND. lSubGrp0

ng = GroupInfo%ng; nfsr = Core%nCoreFsr; nxy = Core%nxy
nRefTemp = DcplInfo%nRefTemp; nRefPln = DcplInfo%nRefPln
#ifndef MPI_ENV
myRefPlnBeg = 1; myRefPlnEnd = nRefPln
lMPI = .FALSE.
#else
myRefPlnBeg = PE%myRefPlnBeg; myRefPlnEnd = PE%myRefPlnEnd
lMPI = PE%lDcplParallel; myrank = PE%myCMFDrank
nproc = PE%nCmfdProc; Comm = PE%MPI_CMFD_Comm
Master = PE%Master
#endif

DcplCmInfo => DcplInfo%DcplCmInfo; DcplFmInfo => DcplInfo%DcplFmInfo

!Iteration
DO iRefPln = myRefPlnBeg, myRefPlnEnd
  iRefTemp = 1;
  WRITE(mesg, '(A, I5)') 'Decoupled Planar MOC - Ref. Plane No.', iRefPln
#ifndef MPI_ENV  
  CALL message(io8, TRUE, TRUE, mesg)
#else
  IF(.NOT. lMPI) CALL message(io8, TRUE, TRUE, mesg)
  IF(lMPI) THEN
    CALL MPI_SeqMessage(mesg, io8, TRUE, TRUE, myrank, nproc, comm)
    CALL MPI_DcplMessage(mesg, myrank, TRUE, TRUE, FALSE)
  ENDIF
#endif
  IF(lSubGrp) CALL SubGrpFsp(Core, DcplFmInfo(iRefTemp, iRefPln)%Fxr,  &
                             THInfo(iRefPln), RayInfo, GroupInfo, DcplCntl, DcplPE(iRefPln))
  !
  DcplItrCntl(iRefPln)%DcplXsGenCntl%iRefPln = iRefPln; DcplItrCntl(iRefPln)%DcplXsGenCntl%iRefTemp = iRefTemp
  ipln = DcplPE(iRefPln)%myzb; lFuelPln = Core%lFuelPlane(ipln)
  nXsGenFbIterMax = DcplItrCntl(iRefPln)%DcplXsGenCntl%nXsGenFbIterMax(1:3)
  DO imode = 1, 3
  !!CALL SetThXsGenCondition(Core, DcplFmInfo(iRefTemp, iRefPln)%Fxr(:, ipln), ThInfo(iRefPln), ipln, 3, 0.3_8)
    iRefTemp = imode
#ifndef MPI_ENV
    IF(.NOT. lFuelPln .AND. imode .EQ. 2) CYCLE
#else
#endif
    DcplItrCntl(iRefPln)%DcplXsGenCntl%nXsGenIterMax = nXsGenFbIterMax(imode)
    IF(imode .NE. 1 ) THEN
      DcplInfo%eigv(iRefTemp, iRefPln) = DcplInfo%eigv(1, iRefPln)
      !IF(lsubgrp) THEN
        DO i = 1, ng 
          CALL CP_VA(DcplFmInfo(iRefTemp, iRefPln)%phis(1:nfsr, ipln:ipln, 1:ng),    &
                    DcplFmInfo( 1, iRefPln)%phis(1:nfsr, ipln:ipln, 1:ng), nfsr, 1, ng)
        ENDDO
      !ENDIF
    ENDIF
    WRITE(mesg, '(A, I5)') 'Decoupled Planar MOC - T-H Condition Case :', imode
#ifndef MPI_ENV
    CALL message(io8, TRUE, TRUE, mesg)
#else
    IF(Master) CALL message(io8, TRUE, TRUE, mesg)
    IF(lMPI) CALL MPI_DcplMessage(mesg, myrank, TRUE, FALSE, FALSE)
    IF(.NOT. lFuelPln .AND. imode .EQ. 2) THEN
      IF(lMPI) CALL Idle_DcplPlnMOC(ng, GroupInfo, DcplItrCntl(iRefPln), DcplPE(iRefPln), lMOC)
      CYCLE
    ENDIF
#endif
    CALL CP_VA(DcplInfo%DcplFmInfo(iRefTemp, iRefPln)%AxPXS(1:nxy, ipln:ipln, 1:ng)           &
                   ,DcplInfo%DcplFmInfo(1, iRefPln)%AxPXS(1:nxy, ipln:ipln, 1:ng), nxy, 1, ng) 
    CALL SetThXsGenCondition(Core, DcplFmInfo(iRefTemp, iRefPln)%Fxr, ThInfo(iRefPln),           &
                             iRefPln, ipln, imode, Frac(imode), DcplCntl%BoronPPM)    
    CALL DcplPlnMOC(Core, DcplFmInfo(iRefTemp, iRefPln), DcplCmInfo(iRefTemp, iRefPln),                   &
                    ThInfo(iRefPln), RayInfo, GroupInfo, DcplInfo%eigv(iRefTemp, iRefPln),        &
                     ng, DcplCntl, DcplItrCntl(iRefPln), DcplPE(iRefPln), lMOC)
    IF(DcplCntl%XsFtnMod .EQ. 2) CALL UpdtMicroXsFtn(Core, DcplFmInfo(iRefTemp, iRefPln), DcplInfo%MicXsFtn(1:nxy, iRefPln),   &
                        iRefPln, nxy, imode, GroupInfo, DcplPE(iRefPln))
    CALL SaveXsGenInfo(Core, DcplFmInfo(iRefTemp, iRefPln), DcplCmInfo(iRefTemp, iRefPln),   &
                       imode, DcplPE(iRefPln))

  ENDDO
#ifdef MPI_ENV  
  IF(lMPI) CALL MPI_DcplMessage(mesg, myrank, TRUE, FALSE, TRUE)
#endif  
ENDDO
END SUBROUTINE

SUBROUTINE SetThXsGenCondition(Core, Fxr, ThInfo, iRefPln, ipln, imod, frac, BoronPPM)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type        ,FxrInfo_Type        ,THInfo_Type
USE DcplTh_Mod,       ONLY : DcplFuelTChg         ,DcplModTChg
USE Boron_mod,        ONLY : SetBoronCoolant
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(ThInfo_Type) :: ThInfo
INTEGER :: ipln, imod, iRefPln
REAL :: frac
REAL :: BoronPPM

SELECT CASE(imod)
  CASE(1)
    !Do Nothing
  CASE(2)  !Fuel Temperature Increase
    CALL DcplFuelTChg(Core, Fxr, ThInfo, iRefPln, ipln, frac)
  CASE(3)  !Fuel Temperature Decrease
    CALL DcplModTChg(Core, Fxr, ThInfo, iRefPln, ipln, 2, frac)
    IF(BoronPPM .GT. 0._8) CALL SetBoronCoolant(Core, Fxr, BoronPPM, ipln, ipln)
  CASE(4)  !Modeartor Temperature Increase
    !CALL DcplFuelTChg(Core, Fxr(:), ThInfo, ipln, -frac)
  CASE(5)  !Moderator Temperature Decrease
    CALL DcplModTChg(Core, Fxr, ThInfo, iRefPln, ipln, 1, -frac)
END SELECT

END SUBROUTINE

SUBROUTINE SaveXsGenInfo(Core, FmInfo, CmInfo, imode, DcplPE)
USE PARAM
USE TYPEDEF,     ONLY : CoreInfo_Type         ,FmInfo_Type         ,PE_Type       &
                       ,FxrInfo_Type          ,Pin_Type            ,Cell_Type     &
                       ,CmInfo_Type           ,PinXs_Type
USE BasicOperation, ONLY : CP_VA
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(PE_Type) :: DcplPE
INTEGER :: imode

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(PinXs_Type), POINTER :: PinXs(:, :)

INTEGER :: ipln ,ixy, icel, ifxr
INTEGER :: nxy, nCoreFsr, nCoreFxr
INTEGER :: niso, nFsrInFxr, FsrIdxSt, FxrIdxSt, nlocalFxr
INTEGER :: j

Fxr => FmInfo%Fxr;    Pin => Core%Pin
CellInfo => Core%CellInfo
PinXs => CmInfo%PinXS
nxy = Core%nxy

ipln = DcplPE%myzb
DO ixy = 1, nxy
  FxrIdxSt = Pin(ixy)%FxrIdxSt
  icel = Pin(ixy)%Cell(ipln)
  nlocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nlocalFxr
    ifxr = FxrIdxSt + j - 1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    niso = Fxr(ifxr, ipln)%nIso
    Fxr(ifxr, ipln)%Dcpl_Temp(imode) = Fxr(ifxr, ipln)%Temp
    Fxr(ifxr, ipln)%Dcpl_Pnum(1 : niso, imode) = Fxr(ifxr, ipln)%pnum(1 : niso)
    Fxr(ifxr, ipln)%Dcpl_FresoA(:, imode) = Fxr(ifxr, ipln)%FresoA(:)
    Fxr(ifxr, ipln)%Dcpl_FresoS(:, imode) = Fxr(ifxr, ipln)%FresoS(:)
    Fxr(ifxr, ipln)%Dcpl_FresoF(:, imode) = Fxr(ifxr, ipln)%FresoF(:)
  ENDDO
ENDDO

IF(imode .NE. 1) THEN
  DO ixy = 1, nxy
    FxrIdxSt = Pin(ixy)%FxrIdxSt
    icel = Pin(ixy)%Cell(ipln)
    nlocalFxr = CellInfo(icel)%nFxr
    DO j = 1, nlocalFxr
      ifxr = FxrIdxSt + j - 1
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
      niso = Fxr(ifxr, ipln)%nIso
      Fxr(ifxr, ipln)%Temp = Fxr(ifxr, ipln)%Dcpl_Temp(1)
      Fxr(ifxr, ipln)%pnum(1 : niso) = Fxr(ifxr, ipln)%Dcpl_Pnum(1 : niso, 1)
      Fxr(ifxr, ipln)%FresoA(:) = Fxr(ifxr, ipln)%Dcpl_FresoA(:, 1)
      Fxr(ifxr, ipln)%FresoS(:) = Fxr(ifxr, ipln)%Dcpl_FresoS(:, 1)
      Fxr(ifxr, ipln)%FresoF(:) = Fxr(ifxr, ipln)%Dcpl_FresoF(:, 1)
    ENDDO
  ENDDO
ENDIF
#ifdef DcplMemSave
DO ixy = 1, nxy
  PinXS(ixy, ipln)%Dcpl_Dhat(:, :, imode) = PinXS(ixy, ipln)%Dhat(:, :)
ENDDO
#endif
ENDSUBROUTINE
