MODULE SUBGRP_MOD
!USE PARAM
!USE TYPEDEF,          ONLY : CoreInfo_Type,   FxrInfo_Type,   PE_TYPE
!USE CNTL,             ONLY : nTracerCntl_Type
IMPLICIT NONE
INTERFACE
  SUBROUTINE SubGrpFsp(Core, Fxr, THInfo, RayInfo,  GroupInfo, nTracerCntl, PE)
    USE PARAM
    USE TYPEDEF,      ONLY : CoreInfo_Type,    FxrInfo_Type,   THInfo_Type,  RayInfo_Type,  GroupInfo_Type,   PE_Type     
    USE CNTL,         ONLY : nTracerCntl_Type
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FxrInfo_Type),POINTER :: Fxr(:, :)
    TYPE(THInfo_Type) :: THInfo
    TYPE(RayInfo_Type) :: RayInfo
    TYPE(GroupInfo_Type) :: GroupInfo
    TYPE(nTracerCntl_Type) :: nTracerCntl
    TYPE(PE_TYPE) :: PE
  END SUBROUTINE
  SUBROUTINE SubGrpFsp_CAT(Core, Fxr, THInfo, RayInfo, GroupInfo, nTracerCntl, PE)
    USE PARAM
    USE TYPEDEF,  ONLY : CoreInfo_Type,   FxrInfo_Type,   GroupInfo_Type,    RayInfo_Type,    PE_TYPE, THInfo_Type
    USE CNTL,     ONLY : nTracerCntl_Type
    IMPLICIT NONE                            
    TYPE(CoreInfo_Type) :: Core
    TYPE(FxrInfo_Type),POINTER :: Fxr(:, :)
    TYPE(THInfo_Type) :: THInfo
    TYPE(RayInfo_Type) :: RayInfo
    TYPE(GroupInfo_Type) :: GroupInfo  
    TYPE(nTracerCntl_Type) :: nTracerCntl
    TYPE(PE_TYPE) :: PE
  END SUBROUTINE  
  SUBROUTINE SubGrpFsp_ISO(Core, Fxr, THInfo, RayInfo, GroupInfo, nTracerCntl, PE)
    USE PARAM
    USE TYPEDEF,  ONLY : CoreInfo_Type,   FxrInfo_Type,   GroupInfo_Type,    RayInfo_Type,    PE_TYPE, THInfo_Type
    USE CNTL,     ONLY : nTracerCntl_Type
    IMPLICIT NONE                            
    TYPE(CoreInfo_Type) :: Core
    TYPE(FxrInfo_Type),POINTER :: Fxr(:, :)
    TYPE(THInfo_Type) :: THInfo
    TYPE(RayInfo_Type) :: RayInfo
    TYPE(GroupInfo_Type) :: GroupInfo  
    TYPE(nTracerCntl_Type) :: nTracerCntl
    TYPE(PE_TYPE) :: PE
  END SUBROUTINE  
    
  SUBROUTINE SetPlnLsigP(Siglp, Sigtr, lv, irc, Core, Fxr, ilv, iz, ig, PE) 
    USE TYPEDEF, ONLY : CoreInfo_Type, FxrInfo_Type, PE_TYPE
    IMPLICIT NONE               
    TYPE(CoreInfo_Type) :: Core
    TYPE(FxrInfo_Type),POINTER :: Fxr(:, :)
    TYPE(PE_TYPE) :: PE
    INTEGER,INTENT(IN) :: irc, ilv, iz, ig
    REAL,INTENT(IN) :: lv
    REAL,POINTER :: Siglp(:), Sigtr(:)
  END SUBROUTINE
  SUBROUTINE SetPlnLsigP_cat(Siglp, Sigtr, lv, icat, Core, Fxr, ilv, iz, ig, PE) 
    USE TYPEDEF,     ONLY : CoreInfo_Type, FxrInfo_Type, PE_Type 
    IMPLICIT NONE               
    TYPE(PE_Type) :: PE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FxrInfo_Type),POINTER :: Fxr(:, :), myFxr
    INTEGER,INTENT(IN) :: icat, ilv, iz, ig
    REAL,INTENT(IN) :: lv
    REAL,POINTER :: Siglp(:), Sigtr(:)
  END SUBROUTINE
  SUBROUTINE SetPlnLsigP_1gMLG(Siglp, Sigtr, lv, Core, Fxr, iz, lCLD, lAIC, PE) 
    USE TYPEDEF, ONLY : CoreInfo_Type, FxrInfo_Type, PE_TYPE 
    IMPLICIT NONE               
    TYPE(CoreInfo_Type) :: Core
    TYPE(FxrInfo_Type),POINTER :: Fxr(:, :)
    TYPE(PE_TYPE) :: PE
    INTEGER,INTENT(IN) :: iz
    REAL,INTENT(IN) :: lv
    REAL,POINTER :: Siglp(:), Sigtr(:)
    LOGICAL,INTENT(IN) :: lCLD, lAIC
  END SUBROUTINE  
  SUBROUTINE SetPlnLsigP_Dancoff(Siglp, Sigtr, Core, Fxr, iz, PE) 
    USE TYPEDEF, ONLY : CoreInfo_Type, FxrInfo_Type, PE_TYPE 
    IMPLICIT NONE               
    TYPE(CoreInfo_Type) :: Core
    TYPE(FxrInfo_Type),POINTER :: Fxr(:, :)
    TYPE(PE_TYPE) :: PE
    INTEGER,INTENT(IN) :: iz
    REAL,POINTER :: Siglp(:), Sigtr(:)
  END SUBROUTINE  
  SUBROUTINE SetPlnLsigP_DancoffAIC(Siglp, Sigtr, Core, Fxr, iz, PE) 
    USE TYPEDEF, ONLY : CoreInfo_Type, FxrInfo_Type, PE_TYPE 
    IMPLICIT NONE               
    TYPE(CoreInfo_Type) :: Core
    TYPE(FxrInfo_Type),POINTER :: Fxr(:, :)
    TYPE(PE_TYPE) :: PE
    INTEGER,INTENT(IN) :: iz
    REAL,POINTER :: Siglp(:), Sigtr(:)
  END SUBROUTINE  
  SUBROUTINE SetPlnLsigP_MLG(Siglp, Sigtr, lv, Core, Fxr, ilv, iz, ig, PE) 
    USE TYPEDEF, ONLY : CoreInfo_Type, FxrInfo_Type, PE_TYPE 
    IMPLICIT NONE               
    TYPE(CoreInfo_Type) :: Core
    TYPE(FxrInfo_Type),POINTER :: Fxr(:, :)
    TYPE(PE_TYPE) :: PE
    INTEGER,INTENT(IN) :: ilv, iz, ig
    REAL,INTENT(IN) :: lv
    REAL,POINTER :: Siglp(:), Sigtr(:)
  END SUBROUTINE
  
  SUBROUTINE SetSubGrpSrc1g(src, SigLamPot, xstr1g, Core, Fxr, iz, PE)
  USE TYPEDEF, ONLY : CoreInfo_Type, FxrInfo_Type, PE_TYPE 
  IMPLICIT NONE
  TYPE(CoreInfo_Type) :: Core
  TYPE(FxrInfo_Type), POINTER  :: Fxr(:, :)
  TYPE(PE_TYPE) :: PE
  INTEGER,INTENT(IN) :: iz
  REAL,POINTER,INTENT(IN) :: xstr1g(:), SigLamPot(:)
  REAL,POINTER :: src(:)
  END SUBROUTINE
  
  SUBROUTINE EquipXSGen(phis1g, SigLamPot, xstr1g, ilv, irc, Core, Fxr, iz, ig, PE)
    USE Typedef,  ONLY : CoreInfo_Type, FxrInfo_Type, PE_TYPE
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
    TYPE(PE_TYPE) :: PE
    INTEGER,INTENT(IN) :: iz, ilv, irc, ig
    REAL, POINTER,INTENT(IN) :: phis1g(:),SigLamPot(:),xstr1g(:)
  END SUBROUTINE
  SUBROUTINE EquipXSGen_1gMLG(phis1g, SigLamPot, xstr1g, ilv, Core, Fxr, iz, lCLD, lAIC, PE)
    USE Typedef,  ONLY : CoreInfo_Type, FxrInfo_Type, PE_TYPE
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
    TYPE(PE_TYPE) :: PE
    INTEGER,INTENT(IN) :: iz, ilv
    REAL,POINTER,INTENT(IN) :: phis1g(:),SigLamPot(:), xstr1g(:)
    LOGICAL,INTENT(IN) :: lCLD, lAIC
  END SUBROUTINE
  SUBROUTINE EquipXSGen_MLG(phis1g, SigLamPot, xstr1g, ilv, Core, Fxr, iz, ig, PE)
    USE Typedef,  ONLY : CoreInfo_Type, FxrInfo_Type, PE_TYPE
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
    TYPE(PE_TYPE) :: PE
    INTEGER,INTENT(IN) :: iz, ilv, ig
    REAL,POINTER,INTENT(IN) :: phis1g(:),SigLamPot(:), xstr1g(:)
  END SUBROUTINE

  FUNCTION SubGrpFspErr(phis1g, phis1gd, nfsr, PE)
    USE Typedef,  ONLY : PE_TYPE
    USE PARAM
    IMPLICIT NONE
    TYPE(PE_TYPE) :: PE
    INTEGER :: nfsr
    REAL :: SubGrpFspErr
    REAL, POINTER :: phis1g(:), phis1gd(:)
  END FUNCTION
  
  SUBROUTINE SubGrpEffXsGen(Core, Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
    !Effective Xs Generation
    USE PARAM
    USE TYPEDEF,      ONLY : CoreInfo_Type,       FxrInfo_Type,          GroupInfo_Type,     &
                             PE_Type,             THInfo_Type
    USE CNTL,           ONLY : nTracerCntl_Type
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
    TYPE(GroupInfo_TYPE) :: GroupInfo
    TYPE(THInfo_Type) :: THInfo
    TYPE(nTracerCntl_Type) :: nTracerCntl
    TYPE(PE_Type) :: PE
    REAL :: eigv

  END SUBROUTINE  
  
  SUBROUTINE FxrChiGen(Core, Fxr, FmInfo, GroupInfo, PE, nTracerCntl)
    USE TYPEDEF,        ONLY :   CoreInfo_Type,     FxrInfo_Type,  FmInfo_Type,            &
                                 GroupInfo_Type,    PE_Type   
    USE CNTL,           ONLY : nTracerCntl_type
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
    TYPE(FmInfo_Type) :: FmInfo
    TYPE(GroupInfo_Type) :: GroupInfo
    TYPE(PE_Type) :: PE
    TYPE(nTracerCntl_Type) :: nTracerCntl
  END SUBROUTINE
  
  SUBROUTINE CalcDancoff(Core, Fxr, RayInfo, THInfo, nTracerCntl, PE)
    USE PARAM
    USE TYPEDEF,  ONLY : CoreInfo_Type,   FxrInfo_Type,   RayInfo_Type,  THInfo_Type,   PE_TYPE
    USE CNTL,     ONLY : nTracerCntl_Type
    IMPLICIT NONE                            
    TYPE(CoreInfo_Type) :: Core
    TYPE(FxrInfo_Type),POINTER :: Fxr(:, :)
    TYPE(RayInfo_Type) :: RayInfo
    TYPE(THInfo_Type) :: THInfo
    TYPE(nTracerCntl_Type) :: nTracerCntl
    TYPE(PE_TYPE) :: PE
  END SUBROUTINE
  SUBROUTINE Set_Dancoff(Siglp, Sigtr, phi, Core, Fxr, iz, nitertot, PE) 
    USE TYPEDEF,     ONLY : CoreInfo_Type, FxrInfo_Type, PE_Type 
    IMPLICIT NONE               
    TYPE(PE_TYPE) :: PE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FxrInfo_Type),POINTER :: Fxr(:, :)
    INTEGER,INTENT(IN) :: iz
    INTEGER :: nitertot
    REAL,POINTER :: Siglp(:), Sigtr(:), phi(:)
  END SUBROUTINE
  SUBROUTINE Set_DancoffAIC(Siglp, Sigtr, phi, Core, Fxr, iz, nitertot, PE) 
    USE TYPEDEF,     ONLY : CoreInfo_Type, FxrInfo_Type, PE_Type 
    IMPLICIT NONE               
    TYPE(PE_TYPE) :: PE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FxrInfo_Type),POINTER :: Fxr(:, :)
    INTEGER,INTENT(IN) :: iz
    INTEGER :: nitertot
    REAL,POINTER :: Siglp(:), Sigtr(:), phi(:)
  END SUBROUTINE
  SUBROUTINE CalcEscXSCP(Core, Fxr, GroupInfo, THInfo, nTracerCntl, PE)
    USE TYPEDEF,        ONLY : coreinfo_type, Fxrinfo_type, GroupInfo_Type, PE_TYPE, THInfo_Type
    USE CNTL,           ONLY : nTracerCntl_Type
    IMPLICIT NONE 
    TYPE(CoreInfo_Type) :: Core
    TYPE(FxrInfo_Type),POINTER :: Fxr(:, :)
    TYPE(THInfo_Type) :: THInfo
    TYPE(GroupInfo_Type) :: GroupInfo
    TYPE(nTracerCntl_Type) :: nTracerCntl
    TYPE(PE_TYPE) :: PE    
  END SUBROUTINE
  SUBROUTINE UpdtNDAF(id, ilv, irc, Core, Fxr, iz, ig, PE)
    USE Typedef,      ONLY : CoreInfo_Type, FxrInfo_Type, PE_Type
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
    TYPE(PE_TYPE) :: PE
    INTEGER,INTENT(IN) :: id, ilv, irc, iz, ig
  END SUBROUTINE
  SUBROUTINE UpdtNDAF_CAT(ilv, icat, Core, Fxr, iz, ig, PE)
    USE Typedef,      ONLY : CoreInfo_Type, FxrInfo_Type, PE_Type
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
    TYPE(PE_TYPE) :: PE
    INTEGER,INTENT(IN) :: ilv, icat, iz, ig
  END SUBROUTINE 
  SUBROUTINE UpdtFnAdj(Core, Fxr, ig, iz, PE) 
    USE Typedef,      ONLY : CoreInfo_Type, FxrInfo_Type, PE_Type
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)  
    TYPE(PE_TYPE) :: PE
    INTEGER,INTENT(IN) :: iz, ig
  END SUBROUTINE
  SUBROUTINE UpdtFtAdj(Core, Fxr, ilv, ig, iz, PE) 
    USE Typedef,      ONLY : CoreInfo_Type, FxrInfo_Type, PE_Type
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)  
    TYPE(PE_TYPE) :: PE
    INTEGER,INTENT(IN) :: ilv, ig, iz
  END SUBROUTINE
  
  SUBROUTINE UpdtCoreIsoInfo(Core, Fxr, PE)
    USE TYPEDEF,      ONLY : CoreInfo_Type,    FxrInfo_Type,   PE_Type   
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(FxrInfo_Type),POINTER :: Fxr(:, :)
    TYPE(PE_TYPE) :: PE
  END SUBROUTINE

  SUBROUTINE UpdtFnAdj_pin(FnAdj,Fxr,TempAvgsq,ig1,ig2,ridx,idx,nlocalFxr,iz) 
    USE TYPEDEF,        ONLY : Fxrinfo_type
    IMPLICIT NONE
    TYPE(FxrInfo_Type),POINTER :: Fxr(:, :)
    INTEGER :: nlocalFxr
    INTEGER :: ridx(nlocalFxr),idx(nlocalFxr),ig1,ig2,iz
    REAL :: FnAdj(nlocalFxr,ig1:ig2),TempAvgsq
  END SUBROUTINE
  SUBROUTINE GetSigtrlp_pin(LocalSiglp,LocalSiglpC,Fxr,ig1,ig2,ridx,idx,FxrIdxSt,nlocalFxr,iz)    
    USE TYPEDEF,        ONLY : Fxrinfo_type
    IMPLICIT NONE
    TYPE(FxrInfo_Type),POINTER :: Fxr(:, :)
    INTEGER :: nlocalFxr
    INTEGER :: ridx(nlocalFxr),idx(nlocalFxr),ig1,ig2,iz,FxrIdxSt
    REAL :: LocalSiglp(nlocalFxr,ig1:ig2)
    REAL :: LocalSiglpC(nlocalFxr)
  END SUBROUTINE
  SUBROUTINE UpdtFtAdj_pin(FtAdj,LocalSiglp,Fxr,ilv,ig,ridx,idx,FxrIdxSt,nlocalFxr,TempAvgsq,iz) 
    USE TYPEDEF,        ONLY : Fxrinfo_type
    IMPLICIT NONE
    TYPE(FxrInfo_Type),POINTER :: Fxr(:, :)
    INTEGER :: nlocalFxr
    INTEGER :: ridx(nlocalFxr),idx(nlocalFxr),ilv,FxrIdxSt,iz,ig
    REAL :: LocalSiglp(nlocalFxr),FtAdj(nlocalFxr),TempAvgsq
  END SUBROUTINE
  SUBROUTINE Reset1Dpingeom(RP,eR,Core,icel,nlocalfxr)  
    USE TYPEDEF,  ONLY : CoreInfo_Type,ResVarPin_Type
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    integer :: nlocalfxr,icel
    real :: eR
    TYPE(ResVarPin_Type), POINTER :: RP
  END SUBROUTINE
END INTERFACE
END MODULE