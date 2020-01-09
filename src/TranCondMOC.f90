#include <defines.h>
  
SUBROUTINE CheckCondMOC(Core, FmInfo, CmInfo, TranInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,            ONLY : CoreInfo_Type,         FmInfo_Type,         CmInfo_Type,    &
                               TranInfo_Type,         GroupInfo_Type,      TranCntl_Type,  &
                               ThInfo_Type,           PE_Type,             PinXS_Type,     &
                               FxrInfo_Type,          Cell_Type,           Pin_Type
USE CNTL,               ONLY : nTracerCntl_Type
USE TRAN_MOD,           ONLY : UpdtHomXsTr1g,         UpdtSigA1g,          UpdtThermalSigA,  &
                               UpdtSigA4g
USE MOC_Mod,            ONLY : FxrAvgPhi
#ifdef MPI_ENV
USE MpiComm_mod,        ONLY : REDUCE,                MPI_MAX_REAL
#endif
USE BasicOperation,     ONLY : CP_VA
USE IOUTIL,             ONLY : message
USE FILES,              ONLY : io8
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(THInfo_Type) :: THInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(PinXS_TYPE), POINTER :: PinXS(:, :)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(Cell_Type), POINTER :: CellInfo(:)
REAL, POINTER :: PhiC(:, :, :), Phis(:, :, :)
REAL :: PhiG(ngmax)

INTEGER :: ng, nxy, myzb, myze, nz
INTEGER :: FxrIdxSt, nLocalFxr
INTEGER :: ixy, iz, ifxr, icel
INTEGER :: i, ig

REAL :: XsChg, TempChg, XsChgCrit, XsChgMax
REAL :: Xs4gChg(4), Xs4gChgMax(4)
REAL :: PhiFxr(1000)
INTEGER :: MocUpdt, MocUpdtG

ng = GroupInfo%ng; myzb = PE%myzb; myze = PE%myze
nz = PE%nz; nxy = Core%nxy
PinXS => CmInfo%PinXS; PhiC => CmInfo%PhiC
Fxr => FmInfo%Fxr; Phis => FmInfo%Phis
CellInfo => Core%CellInfo; Pin => Core%Pin
XsChgCrit = 0.0005_8

TranCntl%lMocUpdt = .FALSE.
TranCntl%lSGFSPUpdt = .FALSE.

!:PRINT *, 'Hesre'
IF(TranCntl%NowStep .EQ. 0) THEN
  TranInfo%Tfcl_ref_SG = ThInfo%Tfmax
  TranInfo%Tfcl_ref_MOC = ThInfo%Tfmax  
  !RETURN
ENDIF
  
TranCntl%it_woMOC = TranCntl%it_woMOC + 1
TranCntl%it_woSG = TranCntl%it_woSG + 1
IF(TranCntl%lXsPerturb) THEN
  TranCntl%lMocUpdt = .TRUE.
  TranCntl%lSGFSPUpdt = .TRUE.
ENDIF
#define CondiSGFSP
#define FxrXsMonitor
#ifdef TempMonitor
IF(abs(TranInfo%Tfcl_ref_MOC - THinfo%Tfcl_max) .GT. 7._8) TranCntl%lMOCUpdt = .TRUE.
IF(mod(TranCntl%it_woMOC, 10) .EQ. 0) TranCntl%lMocUpdt = .TRUE.
!TranCntl%lMocUpdt = .TRUE. 

IF(nTracerCntl%lFeedback) THEN
  IF(abs(TranInfo%Tfcl_ref_SG-ThInfo%TfCl_max) .GT. 7._8)  TranCntl%lSgFspUpdt = .TRUE. 
  !IF(TranCntl%lSgFspUpdt) TranInfo%TfCl_max0 = ThInfo%TfCl_max
  IF(mod(TranCntl%it_woSG, 10) .EQ. 0) TranCntl%lSgFspUpdt = .TRUE.
ENDIF


!IF(nTracerCntl%lFeedback)
#endif
#ifdef HomXsMonitor

IF(TranCntl%NowStep .EQ. 0) THEN
  DO iz = myzb, myze
    DO ixy = 1, nxy
      CALL CP_VA(Phig(1:ng), PhiC(ixy, iz, 1:ng), ng)
      PinXs(ixy, iz)%xstr1g0 = UpdtHomXsTr1g(PinXs(ixy, iz), PhiG, ng)
                               
    ENDDO
  ENDDO
  IF(nTracerCntl%lXsLib) TranInfo%RefTemp0(1:nz) = ThInfo%RefFuelTemp(1:nz)
  RETURN
ENDIF

MocUpdt = 0; XsChgMax=0
DO iz = myzb, myze
  DO ixy = 1, nxy
    CALL CP_VA(Phig(1:ng), PhiC(ixy, iz, 1:ng), ng)
    PinXs(ixy, iz)%xstr1g = UpdtHomXsTr1g(PinXs(ixy, iz), PhiG, ng)    
    XsChg = abs(PinXs(ixy, iz)%xstr1g - PinXs(ixy, iz)%xstr1g0) / abs(PinXs(ixy, iz)%xstr1g)
    XsChgMax = MAX(XsChgMax, XsChg)
    IF(XsChg  .GT. XsChgCrit) MocUpdt = 1
  ENDDO
ENDDO

#ifdef MPI_ENV
CALL MPI_MAX_REAL(XsChgMax, PE%MPI_CMFD_COMM, .TRUE.)
CALL REDUCE(MocUpdt, MocUpdtg, PE%MPI_CMFD_COMM, .TRUE.)
MocUpdt = MocUpdtg
#endif
IF(MocUpdt .GT. 0) THEN
  TranCntl%lMocUpdt = .TRUE.
ENDIF

TranCntl%lSgFspUpdt = TranCntl%lMocUpdt

IF(TranCntl%lMocUpdt) THEN
  !TranCntl%lMocUpdt = .TRUE.
  DO iz = myzb, myze
    DO ixy = 1, nxy
      PinXs(ixy, iz)%xstr1g0 = PinXs(ixy, iz)%xstr1g    
    ENDDO
  ENDDO
ENDIF

IF(nTracerCntl%lXsLib) THEN
  TranCntl%lSgFspUpdt = .FALSE.  
  DO iz = 1, nz
    TranInfo%RefTemp(iz) = ThInfo%RefFuelTemp(iz)
  ENDDO
  IF(nTracerCntl%lFeedback) THEN
   ! IF(abs(TranInfo%Tfcl_ref_SG-ThInfo%TfCl_max) .GT. 7._8)  TranCntl%lSgFspUpdt = .TRUE. 
    !IF(TranCntl%lSgFspUpdt) TranInfo%TfCl_max0 = ThInfo%TfCl_max
    !IF(mod(TranCntl%it_woSG, 5) .EQ. 0) TranCntl%lSgFspUpdt = .TRUE.
  ENDIF
ENDIF
#endif
#ifdef FxrXsMonitor
!#define thermalXS
IF(TranCntl%NowStep .EQ. 0) THEN
  DO iz = myzb, myze
    DO ixy = 1, nxy
      icel = Pin(ixy)%Cell(iz)
      nLocalFxr = CellInfo(icel)%nFxr; FxrIdxSt = Pin(ixy)%FxrIdxSt
      DO i = 1, nLocalFxr
        ifxr = FxrIdxSt + i - 1
       
        PhiFxr(1:ng) = FxrAvgPhi(Core, Fxr, Phis, ixy, i, iz, ng, PE)
#ifdef thermalXS 
        Fxr(ifxr, iz)%siga0 = UpdtThermalSigA(Fxr(ifxr, iz), PhiFxr, ng, GroupInfo, nTracerCntl)
#else
        Fxr(ifxr, iz)%siga0 = UpdtSigA1g(Fxr(ifxr, iz), PhiFxr, ng, GroupInfo, nTracerCntl)
#endif
        Fxr(ifxr, iz)%siga4g0 = UpdtSigA4g(Fxr(ifxr, iz), PhiFxr, ng, GroupInfo, nTracerCntl)
      ENDDO
    ENDDO
  ENDDO
  IF(nTracerCntl%lXsLib) TranInfo%RefTemp0(1:nz) = ThInfo%RefFuelTemp(1:nz)
ENDIF

MocUpdt = 0; XsChgMax=0
Xs4gChgMax = 0 
DO iz = myzb, myze
  DO ixy = 1, nxy
    icel = Pin(ixy)%Cell(iz)
    nLocalFxr = CellInfo(icel)%nFxr; FxrIdxSt = Pin(ixy)%FxrIdxSt
    DO i = 1, nLocalFxr
      ifxr = FxrIdxSt + i - 1
      PhiFxr(1:ng) = FxrAvgPhi(Core, Fxr, Phis, ixy, i, iz, ng, PE)
#ifdef thermalXS 
        Fxr(ifxr, iz)%siga = UpdtThermalSigA(Fxr(ifxr, iz), PhiFxr, ng, GroupInfo, nTracerCntl)
#else
        Fxr(ifxr, iz)%siga = UpdtSigA1g(Fxr(ifxr, iz), PhiFxr, ng, GroupInfo, nTracerCntl)
#endif
      Fxr(ifxr, iz)%siga4g = UpdtSigA4g(Fxr(ifxr, iz), PhiFxr, ng, GroupInfo, nTracerCntl)
      XsChg = abs(Fxr(ifxr, iz)%siga - Fxr(ifxr, iz)%siga0) / abs(Fxr(ifxr, iz)%siga0)
      XsChgMax = MAX(XsChgMax, XsChg)
      IF(XsChg  .GT. XsChgCrit) MocUpdt = 1
      DO ig = 1, 4
        Xs4gChg(ig) = abs(Fxr(ifxr, iz)%siga4g(ig) - Fxr(ifxr, iz)%siga4g0(ig)) / abs(Fxr(ifxr, iz)%siga4g0(ig))
        Xs4gChgMax(ig) = MAX(Xs4gChgMax(ig), Xs4gChg(ig))
      ENDDO
    ENDDO
  ENDDO
ENDDO

#ifdef MPI_ENV
CALL MPI_MAX_REAL(XsChgMax, PE%MPI_CMFD_COMM, .TRUE.)
CALL REDUCE(MocUpdt, MocUpdtg, PE%MPI_CMFD_COMM, .TRUE.)
MocUpdt = MocUpdtg
DO ig = 1, 4
  CALL MPI_MAX_REAL(Xs4gChgMax(ig), PE%MPI_CMFD_COMM, .TRUE.)
ENDDO
#endif
IF(MocUpdt .GT. 0) THEN
  TranCntl%lMocUpdt = .TRUE.
ENDIF
!TranCntl%lMocUpdt = .TRUE.
IF(abs(TranInfo%Tfcl_ref_MOC - THinfo%Tfcl_max) .GT. 10._8) TranCntl%lMOCUpdt = .TRUE.

#ifndef CondiSGFSP
TranCntl%lSgFspUpdt = TranCntl%lMocUpdt
#else
IF(abs(TranInfo%Tfcl_ref_SG-ThInfo%TfCl_max) .GT. 10._8)  TranCntl%lSgFspUpdt = .TRUE.
#endif
IF(TranCntl%lMocUpdt) THEN
  DO iz = myzb, myze
    DO ifxr = 1, Core%nCoreFxr
      Fxr(ifxr, iz)%siga0 = Fxr(ifxr, iz)%siga
      Fxr(ifxr, iz)%siga4g0 = Fxr(ifxr, iz)%siga4g
    ENDDO
  ENDDO
ENDIF

#endif

IF(TranCntl%lMocUpdt) THEN
  TranCntl%it_woMOC = 0
  TranInfo%TfCl_ref_MOC = ThInfo%TfCl_max
ELSE
  !TranCntl%it_woMOC = TranCntl%it_woMOC + 1
  !TranInfo%TfCl_ref_SG = ThInfo%TfCl_max
ENDIF

IF(TranCntl%lSgfspUpdt) THEN
  TranCntl%it_woSG = 0
  TranInfo%TfCl_ref_SG = ThInfo%TfCl_max
ELSE

ENDIF
IF(PE%MASTER) THEN
!  PRINT '(I5, 2x, A, 3L, F12.6)', TranCntl%NowStep, 'XS Perturb/ MOC / SGFSP', TranCntl%lXsPerturb, TranCntl%lMocUpdt, TranCntl%lSGFSPUpdt, XsChgMax
  WRITE(mesg, '(I5, 2x, A, 3L, F12.6)') TranCntl%NowStep, 'XS Perturb/ MOC / SGFSP', TranCntl%lXsPerturb, TranCntl%lMocUpdt, TranCntl%lSGFSPUpdt, XsChgMax
  CALL message(io8, TRUE, TRUE, mesg)
  !PRINT *, XsChgMax
  WRITE(100, '(I7, 4e15.5)') TranCntl%NowStep, Xs4gChgMax(1:4)
ENDIF

NULLIFY(PinXS, PhiC) 

END SUBROUTINE

SUBROUTINE UpdtBaseXsCondiMOC(Core, FmInfo, CmInfo, TranInfo, ThInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,            ONLY : CoreInfo_Type,         FmInfo_Type,         CmInfo_Type,    &
                               TranInfo_Type,         GroupInfo_Type,      TranCntl_Type,  &
                               ThInfo_Type,           PE_Type,             PinXS_Type,     &
                               FxrInfo_Type,          Cell_Type,           Pin_Type
USE CNTL,               ONLY : nTracerCntl_Type
USE TRAN_MOD,           ONLY : UpdtHomXsTr1g,         UpdtSigA1g,          UpdtThermalSigA,  &
                               UpdtSigA4g
USE MOC_Mod,            ONLY : FxrAvgPhi
#ifdef MPI_ENV
USE MpiComm_mod,        ONLY : REDUCE,                MPI_MAX_REAL
#endif
USE BasicOperation,     ONLY : CP_VA
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(THInfo_Type) :: THInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(PinXS_TYPE), POINTER :: PinXS(:, :)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(Cell_Type), POINTER :: CellInfo(:)
REAL, POINTER :: PhiC(:, :, :), Phis(:, :, :)
REAL :: PhiG(ngmax)

INTEGER :: ng, nxy, myzb, myze, nz
INTEGER :: FxrIdxSt, nLocalFxr
INTEGER :: ixy, iz, ifxr, icel
INTEGER :: i, ig


REAL :: Xs4gChg(4), Xs4gChgMax(4)
REAL :: PhiFxr(1000)

ng = GroupInfo%ng; myzb = PE%myzb; myze = PE%myze
nz = PE%nz; nxy = Core%nxy
PinXS => CmInfo%PinXS; PhiC => CmInfo%PhiC
Fxr => FmInfo%Fxr; Phis => FmInfo%Phis
CellInfo => Core%CellInfo; Pin => Core%Pin

!#define thermalXS
DO iz = myzb, myze
  DO ixy = 1, nxy
    icel = Pin(ixy)%Cell(iz)
    nLocalFxr = CellInfo(icel)%nFxr; FxrIdxSt = Pin(ixy)%FxrIdxSt
    DO i = 1, nLocalFxr
      ifxr = FxrIdxSt + i - 1
       
      PhiFxr(1:ng) = FxrAvgPhi(Core, Fxr, Phis, ixy, i, iz, ng, PE)
#ifdef thermalXS 
      Fxr(ifxr, iz)%siga0 = UpdtThermalSigA(Fxr(ifxr, iz), PhiFxr, ng, GroupInfo, nTracerCntl)
#else
      Fxr(ifxr, iz)%siga0 = UpdtSigA1g(Fxr(ifxr, iz), PhiFxr, ng, GroupInfo, nTracerCntl)
#endif
      Fxr(ifxr, iz)%siga4g0 = UpdtSigA4g(Fxr(ifxr, iz), PhiFxr, ng, GroupInfo, nTracerCntl)
    ENDDO
  ENDDO
ENDDO
IF(PE%MASTER) PRINT *, 'Update Base XS for Conditional MOC'
END SUBROUTINE
SUBROUTINE CheckPhiShape(Core, CmInfo, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,           CmInfo_Type,          TranInfo_Type,        &
                             GroupInfo_Type,          TranCntl_Type,                              &
                             PE_Type
USE CNTL,             ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE

END SUBROUTINE


FUNCTION UpdtHomXsTr1g(PinXS, PhiC, ng)
USE PARAM
USE TYPEDEF,       ONLY : PinXS_Type
IMPLICIT NONE
TYPE(PinXS_Type) :: PinXS
REAL :: PhiC(ng)
REAL :: UpdtHomXsTr1g
INTEGER :: ng

REAL :: RRsum, PhiSum
INTEGER :: ig

RRSum = 0; PhiSum = 0
DO ig = 1, ng
  RRsum = RRsum + PinXS%xstr(ig) * PhiC(ig)
  PhiSum = PhiSum + PhiC(ig)
ENDDO
UpdtHomXsTr1g = RRsum / PhiSum

END FUNCTION


FUNCTION UpdtSigA1g(Fxr, PhiFxr, ng, GroupInfo, nTracerCntl)
USE PARAM
USE TYPEDEF,              ONLY : FxrInfo_Type,         ThInfo_Type,        GroupInfo_Type,       &
                                 XsMAc_Type
USE CNTL,                 ONLY : nTracerCntl_Type
USE MacXsLib_mod,         ONLY : MacXsBase
USE XsUtil_mod,           ONLY : GetXsMacDat,          ReturnXsMacDat
IMPLICIT NONE

REAL :: UpdtSigA1g

TYPE(FxrInfo_Type) :: Fxr
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
REAL :: PhiFxr(ng)
INTEGER :: ng

INTEGER :: ig
INTEGER :: iresogrp1, iresogrp2
REAL :: PhiSum, XsSum


TYPE(XsMac_Type), POINTER :: XsMac

iresoGrp1 = GroupInfo%nofg + 1; iresoGrp2 = GroupInfo%nofg + GroupInfo%norg
CALL GetXsMacDat(XsMac, ng, .FALSE.)
CALL MacXsBase(XSMac, Fxr, 1, ng, ng, 1._8, FALSE, TRUE, TRUE)

IF(Fxr%lres) THEN
  DO ig = iresogrp1, iresogrp2
    XsMac%XsMacA(ig) = XsMac%XsMacA(ig) * Fxr%fresoa(ig)  
    !XsMac%XsMacS(ig) = XsMac%XsMacS(ig) * Fxr%fresos(ig)  
  ENDDO
ENDIF
XsSum = 0; PhiSum = 0;

DO ig = 1, ng
  XsSum = XsSum + PhiFxr(ig) * XsMac%XsMacA(ig)
  PhiSum = PhiSum + PhiFxr(ig)
ENDDO
UpdtSigA1g = XsSum / PhiSum
CALL ReturnXsMacDat(XsMac)
END FUNCTION

FUNCTION UpdtThermalSigA(Fxr, PhiFxr, ng, GroupInfo, nTracerCntl)
USE PARAM
USE TYPEDEF,              ONLY : FxrInfo_Type,         ThInfo_Type,        GroupInfo_Type,       &
                                 XsMAc_Type
USE CNTL,                 ONLY : nTracerCntl_Type
USE MacXsLib_mod,         ONLY : MacXsBase
USE XsUtil_mod,           ONLY : GetXsMacDat,          ReturnXsMacDat
IMPLICIT NONE

REAL :: UpdtThermalSigA

TYPE(FxrInfo_Type) :: Fxr
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
REAL :: PhiFxr(ng)
INTEGER :: ng

INTEGER :: ig, igbeg, igend
INTEGER :: iresogrp1, iresogrp2
REAL :: PhiSum, XsSum


TYPE(XsMac_Type), POINTER :: XsMac

iresoGrp1 = GroupInfo%nofg + 1; iresoGrp2 = GroupInfo%nofg + GroupInfo%norg
CALL GetXsMacDat(XsMac, ng, .FALSE.)
CALL MacXsBase(XSMac, Fxr, 1, ng, ng, 1._8, FALSE, TRUE, TRUE)

IF(Fxr%lres) THEN
  DO ig = iresogrp1, iresogrp2
    XsMac%XsMacA(ig) = XsMac%XsMacA(ig) * Fxr%fresoa(ig)  
    !XsMac%XsMacS(ig) = XsMac%XsMacS(ig) * Fxr%fresos(ig)  
  ENDDO
ENDIF
XsSum = 0; PhiSum = 0;


igbeg = 1; igend = ng;
IF(ng .EQ. 47) igbeg = 36
IF(ng .EQ. 190) igend = 160
DO ig = igbeg, igend
  XsSum = XsSum + PhiFxr(ig) * XsMac%XsMacA(ig)
  PhiSum = PhiSum + PhiFxr(ig)
ENDDO
UpdtThermalSigA = XsSum / PhiSum
CALL ReturnXsMacDat(XsMac)
END FUNCTION

FUNCTION UpdtSigA4g(Fxr, PhiFxr, ng, GroupInfo, nTracerCntl)
USE PARAM
USE TYPEDEF,              ONLY : FxrInfo_Type,         ThInfo_Type,        GroupInfo_Type,       &
                                 XsMAc_Type
USE CNTL,                 ONLY : nTracerCntl_Type
USE MacXsLib_mod,         ONLY : MacXsBase
USE XsUtil_mod,           ONLY : GetXsMacDat,          ReturnXsMacDat
IMPLICIT NONE

REAL :: UpdtSigA4g(4)

TYPE(FxrInfo_Type) :: Fxr
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
REAL :: PhiFxr(ng)
INTEGER :: ng

INTEGER :: ig0, ig, igbeg, igend
INTEGER :: iresogrp1, iresogrp2
INTEGER :: groupspec(2,4)
REAL :: PhiSum, XsSum


TYPE(XsMac_Type), POINTER :: XsMac

iresoGrp1 = GroupInfo%nofg + 1; iresoGrp2 = GroupInfo%nofg + GroupInfo%norg
CALL GetXsMacDat(XsMac, ng, .FALSE.)
CALL MacXsBase(XSMac, Fxr, 1, ng, ng, 1._8, FALSE, TRUE, TRUE)

IF(Fxr%lres) THEN
  DO ig = iresogrp1, iresogrp2
    XsMac%XsMacA(ig) = XsMac%XsMacA(ig) * Fxr%fresoa(ig)  
    !XsMac%XsMacS(ig) = XsMac%XsMacS(ig) * Fxr%fresos(ig)  
  ENDDO
ENDIF

GroupSpec(:,1)=(/1, iresogrp1-1/)
GroupSpec(:,2)=(/iresogrp1, iresogrp2/)
GroupSpec(:,3)=(/iresogrp2+1, 35/)
GroupSpec(:,4)=(/36, ng/)
DO ig0 = 1, 4
  XsSum = 0; PhiSum = 0;
  igbeg = GroupSpec(1, ig0); igend = GroupSpec(2, ig0)
  DO ig = igbeg, igend
    XsSum = XsSum + PhiFxr(ig) * XsMac%XsMacA(ig)
    PhiSum = PhiSum + PhiFxr(ig)
  ENDDO
  UpdtSiga4g(ig0) = XsSum / PhiSum
ENDDO
!UpdtThermalSigA = XsSum / PhiSum
CALL ReturnXsMacDat(XsMac)
END FUNCTION
