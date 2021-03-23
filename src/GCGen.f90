#include <defines.h>
SUBROUTINE GCGen(Core, FmInfo, THInfo, CmInfo, GroupInfo, nTracerCntl, PE, ng)
USE PARAM
USE TYPEDEF,         ONLY : CoreInfo_Type,FmInfo_Type, THInfo_Type, CMInfo_Type,Cell_Type,Pin_Type, &
                            FxrInfo_Type, PinXs_Type, XsMac_Type, GroupInfo_Type,PE_TYPE, Powerdist_TYPE, &
                            Asy_Type, Asyinfo_type
USE BasicOperation,  ONLY : CP_VA, CP_CA, MULTI_VA, MULTI_CA
USE files,           ONLY : caseid
USE XSLIB_MOD
USE CNTL,            ONLY : nTracerCntl_Type
USE CritSpec_mod,    ONLY : GetDiffusionCoeff
USE Core_mod,        ONLY : eigv
USE GC_Mod,          ONLY : iasy, iz, nxya, nxy, nx, ny, nz, hz, lSA, Asy, AsyInfo, iasytype, lPinwise, lPin, ipin, npinlist, ngrp, lgapHom
USE GC_mod,          ONLY : lPDQ_GCL, lASY_GCL, lPIN_GCL, lASY_MIC, lPIN_MIC, lB1
USE FILES,            ONLY : io8
USE IOUTIL,           ONLY : message
IMPLICIT NONE

!--- input variables
TYPE(CoreInfo_Type) :: Core
TYPE(FMInfo_Type) :: FmInfo
TYPE(THInfo_Type) :: THInfo
TYPE(CMInfo_Type) :: CMInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
INTEGER :: ng
INTEGER :: iB1, nB1
!---end of input variables


!----------INITIALIZE-----------
ngrp=nTracerCntl%nGCgrp
lGapHom=(Core%lGap).AND.(nTracerCntl%lGCgapHom)
!--- Allocates and zero initialize
CALL GCInit(Core, CmInfo, ng)
WRITE(mesg, '(A)') 'Generating Group Constants ... INIT'
CALL message(io8, TRUE, TRUE, mesg)

lPinWise=nTracerCntl%lPIN_GCL
lPDQ_GCL=nTracerCntl%lPDQ_GCL
lASY_GCL=nTracerCntl%lASY_GCL
lASY_MIC=nTracerCntl%lASY_MIC
lPIN_GCL=nTracerCntl%lPIN_GCL
lPIN_MIC=nTracerCntl%lPIN_MIC
!--- Group Boundary Search
CALL GroupBoundary(ng)
CALL GCPwrNorm(Core, FmInfo, THInfo, CmInfo, nTracerCntl, GroupInfo, PE, ng)
WRITE(mesg, '(A)') 'Generating Group Constants ... Pnorm'
CALL message(io8, TRUE, TRUE, mesg)

!--- 
DO iasy = 1, nxya
    iasytype=Asy(iasy)%AsyType
    nxy=AsyInfo(iasytype)%nxy
    nx=AsyInfo(iasytype)%nx
    ny=AsyInfo(iasytype)%ny
    DO iz = PE%myzb, PE%myze
        hz=Core%HZ(iz)
        nB1=1
        IF( lSA .AND. nTracerCntl%lXsLib ) nB1=2
        DO iB1 = 1, nB1
        !--- POWER and Burn-up Exposure calculation
        CALL GCBU(FmInfo)
        CALL GCPWR(Core, CmInfo, nTracerCntl, PE, ng)    
        WRITE(mesg, '(A)') 'Generating Group Constants ... Set BU and PWR'
        CALL message(io8, TRUE, TRUE, mesg)
        
        !--- Isotope index search sweep
        CALL GCIso(Core, FmInfo, nTracerCntl, ng)
        !CALL GCIso_new(Core, FmInfo, nTracerCntl, ng)
        !CALL GCIso_fast(Core, FmInfo, nTracerCntl, ng)
        WRITE(mesg, '(A)') 'Generating Group Constants ... Set IsoList'
        CALL message(io8, TRUE, TRUE, mesg)
        IF( nTracerCntl%lGCRST .AND. iB1 .EQ. 1 )THEN
            CALL GCRst(Core, FmInfo, nTracerCntl, ng)
            WRITE(mesg, '(A)') 'Generating Group Constants ... Print GCRst'
            CALL message(io8, TRUE, TRUE, mesg)
        ENDIF
        !--- GET flux spectrum
        IF( lSA )THEN
            IF( iB1 .EQ. 1 )THEN
                lB1=.FALSE.
                CALL GCSpec(Core, FmInfo, GroupInfo, nTracerCntl, PE, ng, nTracerCntl%gc_spec)
            ELSE
                lB1=.TRUE.
                CALL GCSpec(Core, FmInfo, GroupInfo, nTracerCntl, PE, ng, 0)
            ENDIF                
        ENDIF
        
        !---------- HOMOGENIZATION --------------------------------------------------
        CALL GCHom_FAST(Core, THInfo, FmInfo, GroupInfo, nTracerCntl, PE, ng)
        !CALL GCHom_SA(Core, THInfo, FmInfo, GroupInfo, nTracerCntl, PE, ng)
        !CALL GCHom_pin(Core, THInfo, FmInfo, GroupInfo, nTracerCntl, PE, ng)
        
        WRITE(mesg, '(A)') 'Generating Group Constants ... Hom'
        CALL message(io8, TRUE, TRUE, mesg)
        !--- B1 diffusion correction        
        CALL GCB1D(ng)
        !--- Macroscopic to Microscopic XS
        CALL GCMac2Mic(ng)
        WRITE(mesg, '(A)') 'Generating Group Constants ... Mac2Mic'
        CALL message(io8, TRUE, TRUE, mesg)

        !---------- 2G CONDENSATION -----------------------------------------------
        CALL GC2GFlux(ng)
        !--- Diffusion coefficient & Transport XS
        CALL GC2GDiffCoeff(ng) ! << GCCon 안에도 있고 (TA용) 정리해서 GCCon 내부에서 호출하도록 수정해야함
        !--- Microscopic 2G XS
        CALL GCCon(ng, FALSE)
        WRITE(mesg, '(A)') 'Generating Group Constants ... Con'
        CALL message(io8, TRUE, TRUE, mesg)
        
        !--- Assembly Discontinuity Factor
        CALL GCSurf(Core, FmInfo, ng)
        CALL GCADF(Core, ng)
        CALL GCFluxNorm(ng)
        CALL GCJNorm(ng)
        WRITE(mesg, '(A)') 'Generating Group Constants ... ADFs and Norm'
        CALL message(io8, TRUE, TRUE, mesg)
        
        !--------- PDQ file edit-----------------------------------------
        CALL GCEDIT_SA(ng)
        WRITE(mesg, '(A)') 'Generating Group Constants ... Edit'
        CALL message(io8, TRUE, TRUE, mesg)
        CALL GCCLR()            
        WRITE(mesg, '(A)') 'Generating Group Constants ... Clear'
        CALL message(io8, TRUE, TRUE, mesg)
        ENDDO !--- B1 loop
    ENDDO
    IF( lPinwise )THEN
        lPin=.TRUE.
        lB1=.FALSE.
        !CALL SetPinList()
        IF( lGapHom )THEN
            IF( Asy(iasy)%lCentX )THEN
                npinlist=(nx-2)*(ny-1)
            ELSEIF( Asy(iasy)%lCentY )THEN
                npinlist=(nx-1)*(ny-2)
            ELSEIF( Asy(iasy)%lCentXY )THEN
                npinlist=(nx-1)*(ny-1)
            ELSE
                npinlist=(nx-2)*(ny-2)
            ENDIF   
            !npinlist=(SQRT(REAL(nxy))-2)**2
        ELSE
            npinlist=nxy
        ENDIF
        DO iz = PE%myzb, PE%myze
            hz=Core%HZ(iz)
            DO ipin = 1, npinlist
                CALL SetPinList(Asy(iasy)%lCentX,Asy(iasy)%lCentY,Asy(iasy)%lCentXY)
                !--- POWER and Burn-up Exposure calculation
                CALL GCPWR(Core, CmInfo, nTracerCntl, PE, ng)    
                
                !--- Isotope index search sweep
                CALL GCIso(Core, FmInfo, nTracerCntl, ng)
                !--- GET flux spectrum
                IF( lSA )THEN
                    CALL GCSpec(Core, FmInfo, GroupInfo, nTracerCntl, PE, ng, nTracerCntl%gc_spec)
                ENDIF
                
                !---------- HOMOGENIZATION --------------------------------------------------
                CALL GCHom_pin(Core, THInfo, FmInfo, GroupInfo, nTracerCntl, PE, ng)
                
                !--- B1 diffusion correction        
                CALL GCB1D(ng)
                !--- Macroscopic to Microscopic XS
                CALL GCMac2Mic(ng)

                !---------- 2G CONDENSATION -----------------------------------------------
                CALL GC2GFlux(ng)
                !--- Diffusion coefficient & Transport XS
                CALL GC2GDiffCoeff(ng) ! << GCCon 안에도 있고 (TA용) 정리해서 GCCon 내부에서 호출하도록 수정해야함
                !--- Microscopic 2G XS
                CALL GCCon(ng, TRUE)
                
                !--- Assembly Discontinuity Factor
                CALL GCSurf_Pin(Core, FmInfo, ng)
                CALL GCADF_Pin(Core, ng)
                CALL GCFluxNorm(ng)
                CALL GCJNorm(ng)
                
                !--------- PDQ file edit-----------------------------------------
                CALL GCEDIT_SA(ng)
                CALL GCCLR()
            ENDDO !--end of Pin
        ENDDO !--end of Z        
    ENDIF !--end of pinwise 
    lPin=.FALSE.
ENDDO
CALL GCFin()


END SUBROUTINE