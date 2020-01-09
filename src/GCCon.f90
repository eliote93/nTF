SUBROUTINE GroupBoundary(ng)
    !--- Determination of Macro Group 
    USE XSLIB_MOD
    USE GC_mod,      ONLY : nglist, GrpLBList, ngrp, GrpIdx
    USE CNTL,        ONLY : nTracerCntl     !-- EDIT JSU 12/07/2018 
    IMPLICIT NONE
    INTEGER :: ng
    
    REAL :: EStruct(ng+1),Eavg(ng)
    INTEGER :: ig, gidx, iGrpSet
    REAL :: ELB, temp
    
    iGrpSet=0
    SELECT CASE(ngrp)
    CASE(2)
        iGrpSet=1
    CASE(4)
        iGrpSet=2
    CASE(8)
        iGrpSet=3
    CASE(16)
        iGrpSet=4
    ENDSELECT
    IF( iGrpSet .EQ. 0 )THEN
        write(*,*) '    None Excisting Grp Boundary Info ... BYS edit 14/02/11'
    ENDIF    
    
    gidx=1    
    IF (nTracerCntl%libtyp .eq. 2 .OR. nTracerCntl%libtyp .eq. 3)  THEN
        EStruct(1:ng+1) = enbhel(1:ng+1)
        DO ig = 1, ng
          Eavg(ig) = (EStruct(ig) + EStruct(ig+1))/2._8
        ENDDO
        ELB = GrpLBList(gidx, iGrpSet)
        DO ig = 1, ng-1
            temp = (ELB-Eavg(ig))*(ELB-Eavg(ig+1))
            GrpIdx(ig)=gidx
            IF(temp .LE. 0._8) THEN
                gidx=gidx+1
                ELB = GrpLBList(gidx, iGrpSet)
            ENDIF
        ENDDO
        GrpIdx(ng)=ngrp
    ELSE
        IF (ng .eq. 7) THEN
            DO ig = 1, 4
                GrpIdx(ig)=1
            ENDDO
            DO ig = 5, 7
                GrpIdx(ig)=2
            ENDDO
        ELSEIF (ng .eq. 2) THEN
            GrpIdx(1)=1
            GrpIdx(2)=2
        ENDIF
    ENDIF    
    CONTINUE
    
    
ENDSUBROUTINE

SUBROUTINE GC2GFlux(ng)
    USE GC_mod,     ONLY : PhiG, PhiVol, GrpIdx
    IMPLICIT NONE
    INTEGER :: ng
    
    INTEGER :: ig
    INTEGER ::  gidx
    
    PhiG=0
    DO ig = 1, ng
        gidx=GrpIdx(ig)
        PhiG(gidx)=PhiG(gidx)+PhiVol(ig)  !instead of fphi, tphi
    ENDDO    
    
    
ENDSUBROUTINE

SUBROUTINE GC2GDiffCoeff(ng)
    USE GC_mod
    IMPLICIT NONE
    INTEGER :: ng
    
    INTEGER :: ig
    INTEGER ::  gidx
    
    !--- Diffusion coeff.
    DO ig = 1, ng    
        gidx=GrpIdx(ig)
        isoMacXs2g(0,1,gidx)=isoMacXs2g(0,1,gidx)+Dng(ig)*PhiVol(ig)        
    ENDDO
    DO gidx = 1, ngrp
        isoMacXs2g(0,1,gidx)=isoMacXs2g(0,1,gidx)/PhiG(gidx)
    ENDDO
ENDSUBROUTINE


SUBROUTINE GCCon(ng)
    USE GC_mod
    IMPLICIT NONE
    INTEGER :: ng
    
    INTEGER :: ig, ig2
    INTEGER ::  gidx, gidx2
    INTEGER :: iso
    !D, a, r, f, nf, kf
    !ALLOCATE(isoMacXs(nisotot,6,ng), isoMacSm(nisotot,ng,ng),isoMacXs2g(nisotot,6,2)) 
    !tr, a, r, f, nu, k
    !ALLOCATE(isoMicXs(nisotot,6,ng), isoMicSm(nisotot,ng,ng),isoMicXs2g(nisotot,6,2))

    DO iso = 1, nisotot
        DO ig = 1, ng    
            gidx=GrpIdx(ig)
            isoMacXs2g(iso,0,gidx)=isoMacXs2g(iso,0,gidx)+ isoMacXs(iso,0,ig)*PhiVol(ig) ! total
            isoMacXs2g(iso,1,gidx)=isoMacXs2g(iso,1,gidx)+ 1/isoMacXs(iso,1,ig)*PhiVol(ig) ! 1/tr
            isoMacXs2g(iso,2,gidx)=isoMacXs2g(iso,2,gidx)+ isoMacXs(iso,2,ig)*PhiVol(ig) ! a
            isoMacXs2g(iso,3,gidx)=isoMacXs2g(iso,3,gidx)+ isoMacXs(iso,2,ig)*PhiVol(ig) ! r= !a! + s_out
            isoMacXs2g(iso,4,gidx)=isoMacXs2g(iso,4,gidx)+ isoMacXs(iso,4,ig)*PhiVol(ig) ! f
            isoMacXs2g(iso,5,gidx)=isoMacXs2g(iso,5,gidx)+ isoMacXs(iso,5,ig)*PhiVol(ig) ! nf
            isoMacXs2g(iso,6,gidx)=isoMacXs2g(iso,6,gidx)+ isoMacXs(iso,6,ig)*PhiVol(ig) ! kf
            DO ig2 = 1, ng
                gidx2=GrpIdx(ig2)                
                isoMacSm2g(iso,gidx,gidx2)  =isoMacSm2g(iso,gidx,gidx2)  +isoMacSm(iso,ig,ig2)  *PhiVol(ig)
                isoMacSm2gP1(iso,gidx,gidx2)=isoMacSm2gP1(iso,gidx,gidx2)+isoMacSmP1(iso,ig,ig2)*PhiVol(ig)
                IF( gidx .NE. gidx2 )THEN
                    isoMacXs2g(iso,3,gidx)=isoMacXs2g(iso,3,gidx)+ isoMacSm(iso,ig,ig2)*PhiVol(ig) ! r= !a! + s_out
                ENDIF                
            ENDDO    
        ENDDO
    ENDDO
    DO iso = 1, nisotot
        DO gidx = 1, ngrp
            rphi=1.0/PhiG(gidx)
        
            isoMacXs2g(iso,0,gidx)=isoMacXs2g(iso,0,gidx)*rphi  ! total
            isoMicXs2g(iso,0,gidx)=isoMacXs2g(iso,0,gidx)/isoNumden(iso)
            isoMacXs2g(0,0,gidx) = isoMacXs2g(0,0,gidx) + isoMacXs2g(iso,0,gidx)
            
            isoMacXs2g(iso,1,gidx)=isoMacXs2g(iso,1,gidx)*rphi ! 1/tr_iso_G = 3DG
            isoMacXs2g(iso,1,gidx)=1/(isoMacXs2g(iso,1,gidx)) ! 1/tr = DG -> tr_iso_G 
            isoMicXs2g(iso,1,gidx)=isoMacXs2g(iso,1,gidx)/isoNumden(iso)  ! sig_tr
            
            isoMacXs2g(iso,2,gidx)=isoMacXs2g(iso,2,gidx)*rphi  ! a
            isoMicXs2g(iso,2,gidx)=isoMacXs2g(iso,2,gidx)/isoNumden(iso)
            isoMacXs2g(0,2,gidx) = isoMacXs2g(0,2,gidx) + isoMacXs2g(iso,2,gidx)
            
            isoMacXs2g(iso,3,gidx)=isoMacXs2g(iso,3,gidx)*rphi  ! r
            isoMicXs2g(iso,3,gidx)=isoMacXs2g(iso,3,gidx)/isoNumden(iso)
            isoMacXs2g(0,3,gidx) = isoMacXs2g(0,3,gidx) + isoMacXs2g(iso,3,gidx)
            
            IF( isoMacXs2g(iso,4,gidx) .NE. 0 )THEN
                isoMacXs2g(iso,4,gidx)=isoMacXs2g(iso,4,gidx)*rphi  ! f
                isoMicXs2g(iso,4,gidx)=isoMacXs2g(iso,4,gidx)/isoNumden(iso) !f
                isoMacXs2g(0,4,gidx) = isoMacXs2g(0,4,gidx) + isoMacXs2g(iso,4,gidx)
                
                isoMacXs2g(iso,5,gidx)=isoMacXs2g(iso,5,gidx)*rphi  ! isoMac nuf
                isoMicXs2g(iso,5,gidx)=isoMacXs2g(iso,5,gidx)/isoMacXs2g(iso,4,gidx) ! micro nu
                isoMacXs2g(0,5,gidx) = isoMacXs2g(0,5,gidx) + isoMacXs2g(iso,5,gidx) ! Macro nuf
                
                isoMacXs2g(iso,6,gidx)=isoMacXs2g(iso,6,gidx)*rphi  ! kf
                isoMicXs2g(iso,6,gidx)=isoMacXs2g(iso,6,gidx)/isoMacXs2g(iso,4,gidx) ! kappa
                isoMacXs2g(0,6,gidx) = isoMacXs2g(0,6,gidx) + isoMacXs2g(iso,6,gidx)
            ELSE
                isoMacXs2g(iso,4,gidx)=0 ! f
                isoMicXs2g(iso,4,gidx)=0 ! f
                isoMacXs2g(iso,5,gidx)=0 ! nuf
                isoMicXs2g(iso,5,gidx)=0 ! nu
                isoMacXs2g(iso,6,gidx)=0 ! kf
                isoMicXs2g(iso,6,gidx)=0 ! kappa
            ENDIF
            DO gidx2 = 1, ngrp
                isoMacSm2g(iso,gidx,gidx2)=isoMacSm2g(iso,gidx,gidx2)*rphi
                IF( isoNumden(iso) .EQ. 0_8 )THEN
                    isoMicSm2g(iso,gidx,gidx2)=0_8
                ELSE
                    isoMicSm2g(iso,gidx,gidx2)=isoMacSm2g(iso,gidx,gidx2)/isoNumden(iso)  ! debugged in 17/03/10 r563d
                ENDIF                
                isoMacSm2g(0,gidx,gidx2)=isoMacSm2g(0,gidx,gidx2)+isoMacSm2g(iso,gidx,gidx2)
            ENDDO
        ENDDO
    ENDDO  
    ChiG=0
    DO ig = 1, ng
        gidx=GrpIdx(ig)
        ChiG(gidx)=ChiG(gidx)+isoMacXs(0,6,ig)
    ENDDO
ENDSUBROUTINE