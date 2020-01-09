
SUBROUTINE NULLTransient_Driver()
USE PARAM
USE GEOM,          ONLY : Core, ng
USE RAYS,          ONLY : RAYINFO
USE Core_mod,      ONLY : CMInfo,            FmInfo,               THInfo,        &
                          GroupInfo,         eigv
USE PE_MOD,        ONLY : PE
USE CNTL,          ONLY : nTracerCntl
USE ItrCNTL_mod,   ONLY : ItrCntl
USE FILES,         ONLY : io8
USE IOUTIL,        ONLY : message
USE Timer,         ONLY : nTracer_dclock, TimeChk
#ifdef MPI_ENV
USE MPIComm_mod,   ONLY : MPI_SYNC
#endif

USE TRAN_MOD,      ONLY : TranInfo,          TranCntl,                            &
                          InitTransient,     InitPrecursor  
USE TRANCMFD_MOD,  ONLY : 
IMPLICIT NONE 

INTEGER :: istep

!CALL SSEIG()
CALL InitTransient(Core, RayInfo, FmInfo, CmInfo, ThInfo, eigv, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)
!Calculate homgenized chid and beta
CALL HomKineticParamGen(Core, FmInfo, CmInfo, TranInfo, GroupInfo, nTracerCntl, PE)
!Calculate the Precursor Quadratic approx. Coefficient
TranCntl%NowStep = 1
CALL SetCmfdPrecParam(Core, CmInfo, TranInfo, TranCntl, nTracerCntl, PE)

CALL InitNullTransient(Core, FmInfo, CMInfo, GroupInfo, nTracerCntl, PE) 

CALL InitTransient(Core, RayInfo, FmInfo, CmInfo, ThINfo, eigv, TranInfo, GroupInfo, TranCntl, nTracerCntl, PE)
CALL InitPrecursor(Core, FmInfo, CmInfo, TranInfo, GroupInfo, nTracerCntl, PE)

TranCntl%NowStep = 1
WRITE(88, *) 0, TranInfo%PowerLevel
DO istep = 1, 2500
  CALL TranCmfd_Driver(Core, FmInfo, CmInfo, TranInfo, THInfo, GroupInfo, TranCntl, nTracerCntl, ItrCntl, PE)
  CALL UpdtPowerLevel(Core, FmInfo, CmInfo, TranInfo, ng, nTracerCntl, PE)
  WRITE(88, '(I5, F10.5, F10.4)') istep, TranInfo%PowerLevel * 100, TranInfo%Reactivity
  CONTINUE
ENDDO
END SUBROUTINE

SUBROUTINE InitNullTransient(Core, FmInfo, CMInfo, GroupInfo, nTracerCntl, PE) 
USE PARAM
USE TYPEDEF,   ONLY : CoreInfo_Type,      FmInfo_Type,      CmInfo_Type,     &
                      GroupInfo_Type,     PE_Type
USE CNTL,      ONLY : nTRACERCntl_Type
USE BasicOperation, ONLY : CP_CA, MULTI_CA
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE

INTEGER :: myzb, myze, myzbf, myzef
INTEGER :: nxy, nfsr, ng, nprec

nFsr = Core%nCoreFsr; nxy = Core%nxy
myzb = PE%myzb; myze = PE%myze
myzbf = PE%myzbf; myzef = PE%myzef
ng = GroupInfo%ng

nprec = GroupInfo%nprec
!CALL MULTI_CA(0.5_8, FmInfo%Prec(1:nprec, 1:nFsr, myzb:myze), nprec, nfsr, myze - myzb + 1)
!CALL MULTI_CA(0.5_8, CmInfo%PrecCm(1:nprec, 1:nxy, myzb:myze), nprec, nxy, myze - myzb + 1)
!CALL MULTI_CA(0.5_8, CmInfo%PrecFm(1:nprec, 1:nxy, myzbf:myzef), nprec, nxy, myzef - myzbf + 1)

!CALL MULTI_CA(0.5_8, FmInfo%Phis(1:nFsr, myzb:myze, 1:ng), nfsr, myze - myzb + 1, ng)

CALL MULTI_CA(0.5_8, FmInfo%TranPhi(1:nFsr, myzb:myze, 1:ng), nfsr, myze - myzb + 1, ng)
CALL MULTI_CA(0.5_8, FmInfo%TranPsi(1:nFsr, myzb:myze), nfsr, myze - myzb + 1)
!CALL MULTI_CA(0.5_8, FmInfo%TranPsid(1:nFsr, myzb:myze), nfsr, myze - myzb + 1)

!CALL MULTI_CA(0.5_8, CmInfo%PhiC(1:nxy, myzb:myze, 1:ng), nxy, myze - myzb + 1, ng)
CALL MULTI_CA(0.5_8, CmInfo%TranPhiCm(1:nxy, myzb:myze, 1:ng), nxy, myze - myzb + 1, ng)
CALL MULTI_CA(0.5_8, CmInfo%TranPsiCm(1:nxy, myzb:myze), nxy, myze - myzb + 1)
!CALL MULTI_CA(0.5_8, CmInfo%TranPsiCmd(1:nxy, myzb:myze), nxy, myze - myzb + 1)
!CALL MULTI_CA(0.5_8, CmInfo%PhiFm(1:nxy, myzb:myze, 1:ng), nxy, myze - myzb + 1, ng)
!CALL MULTI_CA(0.5_8, CmInfo%TranPhiFm(1:nxy, myzbf:myzef, 1:ng), nxy, myzef - myzbf + 1, ng)
!CALL MULTI_CA(0.5_8, CmInfo%TranPsiFm(1:nxy, myzbf:myzef), nxy, myzef - myzbf + 1)
!CALL MULTI_CA(0.5_8, CmInfo%TranPsiFmd(1:nxy, myzbf:myzef), nxy, myzef - myzbf + 1)

END SUBROUTINE