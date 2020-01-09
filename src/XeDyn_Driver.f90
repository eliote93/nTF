#include <defines.h>

SUBROUTINE XeDyn_Driver()
USE PARAM
USE TYPEDEF,        ONLY : XeDynInfo_Type,    XeDynState_Type,                         &
                           FmInfo_Type,       CmInfo_Type,               PE_Type,           &
                           CoreInfo_Type,     FxrInfo_Type
USE PE_Mod,         ONLY : PE
USE CNTL,           ONLY : nTracerCntl
USE geom,           ONLY : Core
USE Core_mod,       ONLY : FmInfo,             CmInfo,            THInfo,              &
                           GroupInfo,                                                  &
                            xkconv,             xkcrit
USE DEPL_MOD,       ONLY : DeplCntl
USE XeDyn_Mod,      ONLY : XeDynInfo,                                                  &
                           FinalXe,           InitXe,                                   &
                           InitXeDynInfo,     SetXeDynCoreState, UpdtXeDyn
USE CritSpec_mod,    ONLY : XsKinf
USE ioutil,         ONLY : message
USE files,          ONLY : io8
IMPLICIT NONE
INTEGER :: myzb, myze, nxy, nfxr
INTEGER :: istep, istate
LOGICAL :: lPerturb
LOGICAL :: MASTER, SLAVE
nFxr = Core%nCoreFxr
myzb = PE%myzb; myze = PE%myze
nxy = Core%nxy 

Master = PE%Master; Slave = PE%Slave

WRITE(mesg, '(A)') 'Performing Xe Dynamics Calculation ...'
IF(Master) CALL message(io8, TRUE, TRUE, mesg)
IF(.NOT. nTracerCntl%lXeDyn) CALL InitXeDynInfo(FmInfo%Fxr, nFxr, myzb, myze)
IF(XeDynInfo%lPerturb(0)) THEN
  istep = 0
  CALL SetXeDynCoreState(istep, XeDynInfo, nTracerCntl, PE)
ELSE 
  CALL UpdtXeDyn(Core, FmInfo, THInfo, GroupInfo, nTracerCntl, DeplCntl, PE, Mode = InitXe)
ENDIF
DO istep = 1, XeDynInfo%nTimeStep
  Core%DelT = XeDynInfo%TSec(istep) - XeDynInfo%TSec(istep-1)
  IF(XeDynInfo%lPerturb(istep)) CALL SetXeDynCoreState(istep, XeDynInfo, nTracerCntl, PE)
  CALL SSEIG()
  CALL UpdtXeDyn(Core, FmInfo, THInfo, GroupInfo, nTracerCntl, DeplCntl, PE, Mode = FinalXe)
  xkcrit = XsKinf(Core, FmInfo, GroupInfo, .TRUE., nTracerCntl%lCritSpec, PE)

  CALL XeDynBanner(istep)
  
  IF(PE%lCmfdGrp) CALL OutputEdit()
ENDDO



nTracerCntl%lXeDyn = .TRUE.; nTracerCntl%lTrXe = .TRUE.
nTracerCntl%lEqXe = .FALSE.
END SUBROUTINE

SUBROUTINE SetXeDynCoreState(istep, XeDynInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,         ONLY : XeDynInfo_Type,    XeDynState_Type, PE_Type
USE CNTL,            ONLY : nTracerCntl_Type
USE geom,            ONLY : Core
USE Core_mod,        ONLY : FmInfo,             GroupInfo,            THInfo,           &                
                            xkconv,             xkcrit
USE XeDyn_Mod,       ONLY : UpdtXeDyn,                                                  &
                            InitXe,             UpdtXe,               FinalXe
USE DEPL_MOD,        ONLY : DeplCntl
USE CritSpec_mod,    ONLY : XsKinf
USE FILES,           ONLY : io8
USE IOUTIL,          ONLY : message

IMPLICIT NONE
TYPE(XeDynInfo_Type) :: XeDynInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
INTEGER :: istep

INTEGER :: istate
CHARACTER(100) :: boronppm

istate = XeDynInfo%StateMapping(istep)


IF(XeDynInfo%CoreState(istate)%lBoronSearch) THEN
  nTracerCntl%lBoronSearch = .TRUE.
  nTracerCntl%Target_eigv = XeDynInfo%CoreState(istate)%Target_keff
ELSEIF(.NOT. XeDynInfo%CoreState(istate)%lPrevBoroN) THEN
  nTracerCntl%lBoronSearch = .FALSE.
  nTracerCntl%BoronPPM = XeDynInfo%CoreState(istate)%BoronPPM
ENDIF

nTracerCntl%PowerLevel = XeDynInfo%CoreState(istate)%PowLv
ThInfo%PowLv = XeDynInfo%CoreState(istate)%PowLv

IF(nTracerCntl%lBoronSearch) THEN
  WRITE(boronppm, '(x, A, x, F8.5)')  'SEARCH - ', nTracerCntl%Target_eigv
ELSE
  WRITE(boronppm, '(F8.2 x, A)') nTracerCntl%BoronPpm, 'ppm'  
ENDIF


WRITE(mesg, '(A)') 'Update Steady-State Solution for Core State Change ...'
IF(PE%MASTER) CALL message(io8, TRUE, TRUE, mesg)
nTracerCntl%lXeDyn = .FALSE.
CALL SSEIG()
xkcrit = XsKinf(Core, FmInfo, GroupInfo, .TRUE., nTracerCntl%lCritSpec, PE)
nTracerCntl%lXeDyn = .TRUE.
nTracerCntl%lTrXe = .TRUE.
CALL UpdtXeDyn(Core, FmInfo, THInfo, GroupInfo, nTracerCntl, DeplCntl, PE, Mode = InitXe)
WRITE(mesg,'(A, (F7.3, A), A,(F8.3,A),2(A,F8.5),A,F8.2,A)') 'Time = ', XeDynInfo%Th(istep), ' h ', &
        'Burnup = ',DeplCntl%NowBurnUp(2),  ' MWD/kgHM','  keff = ',xkconv,'  kinf(BK) = ', &
        xkcrit,'  ppm = ',nTracerCntl%boronppm,' ppm'
IF(PE%Master) CALL message(io8, FALSE, TRUE, mesg)
!5 + 9

WRITE(mesg, '(A)') hbar1(1:77)
IF(PE%Master) CALL message(io8, FALSE, TRUE, MESG)

END SUBROUTINE

SUBROUTINE XeDynBanner(istep)
USE PARAM
USE TYPEDEF,         ONLY : CoreInfo_Type,      FmInfo_Type,        PE_Type
USE DeplType,        ONLY : DeplCntl_Type
USE geom,            ONLY : Core
USE Core_mod,        ONLY : FmInfo,                                            &
                            xkconv,             xkcrit
USE XeDyn_Mod,       ONLY : XeDynInfo
USE CNTL,            ONLY : nTracerCntl
USE PE_Mod,          ONLY : PE
USE Depl_Mod,        ONLY : DeplCntl,            GetCoreHmMass
USE FILES,           ONLY : io8
USE IOUTIL,          ONLY : message
IMPLICIT NONE
INTEGER :: istep
WRITE(mesg,'(A, (F7.3, A), A,(F8.3,A),2(A,F8.5),A,F8.2,A)') 'Time = ', XeDynInfo%Th(istep), ' h ', &
        'Burnup = ',DeplCntl%NowBurnUp(2),  ' MWD/kgHM','  keff = ',xkconv,'  kinf(BK) = ', &
        xkcrit,'  ppm = ',nTracerCntl%boronppm,' ppm'
IF(PE%Master) CALL message(io8, FALSE, TRUE, mesg)
WRITE(mesg, '(A)') hbar2(1:77)
IF(PE%Master) CALL message(io8, FALSE, TRUE, MESG)
END SUBROUTINE