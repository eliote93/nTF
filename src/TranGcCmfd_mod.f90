Module TranGcCmfd_mod
INTERFACE
SUBROUTINE TranGcCmfdAcc(Core, CmInfo, TranInfo, Eigv, lGcFdk, GroupInfo, GcGroupInfo, TranCntl, nTracerCntl, ItrCntl, PE)
USE PARAM
USE TYPEDEF,       ONLY : CoreInfo_Type,              CmInfo_Type,         GroupInfo_Type,      &
                          TranInfo_Type,              TranCntl_Type,                            &
                          PE_TYPE,                    PinXS_Type
USE Cntl,           ONLY : nTracerCntl_Type
USE ItrCntl_mod,    ONLY : ItrCntl_TYPE
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_TYpe) :: GroupInfo, GcGroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(ItrCntl_TYPE) :: ItrCntl
TYPE(PE_TYPE) :: PE
REAL :: Eigv
LOGICAL :: lGcFdk 
END SUBROUTINE

SUBROUTINE MakeGcKineticParam(Core, CmInfo, TranInfo, GroupInfo, GcGroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,         CmInfo_Type,            TranInfo_Type,         &
                             GroupInfo_Type,        PE_TYPE,                TranCntl_Type,         &
                             PinXS_Type
USE CNTL,             ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo,     GcGroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
END SUBROUTINE

SUBROUTINE SetTranGcCmfdLinearSystem()

END SUBROUTINE

SUBROUTINE SetTranCmfd2GLinearSystem()

END SUBROUTINE

FUNCTION TranGcResidualError(phi, psi, TranSrc,  PE, constsrc)
USE PARAM
USE TYPEDEF,   ONLY : PE_TYPE
USE CMFD_MOD, ONLY : nxy,         myzbf,       myzef,               &
                     src,                                           &
                     SubPlaneMap, PinVol,      PinVolFm,            &
                     hzfm,        PinNeighIdx
IMPLICIT NONE
REAL :: TranGcResidualError
REAL, POINTER :: phi(:, :, :), psi(:, :), TranSrc(:, :, :)
TYPE(PE_TYPE) :: PE
REAL, OPTIONAL :: ConstSrc

END FUNCTION
END INTERFACE
  
END MODULE