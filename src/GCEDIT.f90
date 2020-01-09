SUBROUTINE GCEDIT_SA(ng)
    USE GC_mod
    USE GCpin_mod
    USE Depl_Mod,  ONLY : DeplCntl  !15/03/10 r534 : added for GC gen while depl_calc
    USE FXRVAR_MOD
    USE EFT_MOD
    IMPLICIT NONE
    INTEGER :: ng
    
    INTEGER :: i, j, ix, iy, ixy, nnxy
    INTEGER :: ig, gidx
    
    INTEGER :: iso
    INTEGER :: io2, io3
    CHARACTER(256) :: fn, fn2
    
    io2=41
    io3=42
    WRITE(fn,'(A)') TRIM(caseid)
    IF( lB1 ) WRITE(fn,'(A,A)') 'B1C_', TRIM(caseid)
    IF( lBranchRun )THEN
        SELECT CASe(CurrentVarId)
        CASE(branch_base)
            WRITE(fn,'(A,A,i1)') TRIM(fn), '_base'
        CASE(branch_tmod)
            WRITE(fn,'(A,A,i1)') TRIM(fn), '_tmod',CurrentCaseId
        CASE(branch_tfuel)
            WRITE(fn,'(A,A,i1)') TRIM(fn), '_tfuel',CurrentCaseId
        CASE(branch_rho)
            WRITE(fn,'(A,A,i1)') TRIM(fn), '_rho',CurrentCaseId
        CASE(branch_boron)
            WRITE(fn,'(A,A,i1)') TRIM(fn), '_boron',CurrentCaseId
        ENDSELECT
    ENDIF
    IF( lEFT_GCGEN )THEN
        SELECT CASE(calcid)
        CASE(1) ! reference
            WRITE(fn,'(A,A)') TRIM(fn), '_REF'
        CASE(2) ! Uniform Tmod
            WRITE(fn,'(A,A)') TRIM(fn), '_UMT'
        CASE(3) ! Uniform Tfuel
            WRITE(fn,'(A,A)') TRIM(fn), '_UFT'
        ENDSELECT
        WRITE(fn,'(A,A,i1,A,i1)') TRIM(fn), '_', iboroncase,'_',ipowercase        
    ENDIF
    IF( .NOT. lSA )THEN
        WRITE(fn,'(A,A)') TRIM(fn),'_'
        IF( iasy .LT. 10 )THEN
            WRITE(fn,'(A,A,i1)') TRIM(fn),'0', iasy
        ELSE
            WRITE(fn,'(A,i2)') TRIM(fn), iasy
        ENDIF
    ENDIF
    IF( lZdir )THEN
        WRITE(fn,'(A,A)') TRIM(fn),'_z'
        IF( iz .LT. 10 )THEN
            WRITE(fn,'(A,A,i1)') TRIM(fn),'0', iz
        ELSE
            WRITE(fn,'(A,i2)') TRIM(fn), iz
        ENDIF
    ENDIF
    IF( DeplCntl%nburnupstep .GE. 1 )THEN
        WRITE(fn,'(A,A)') TRIM(fn),'_d'
        IF( DeplCntl%nowstep .LT. 10 )THEN
            WRITE(fn,'(A,A,i1)') TRIM(fn),'00',DeplCntl%nowstep
        ELSEIF( DeplCntl%nowstep .LT. 100 )THEN
            WRITE(fn,'(A,A,i2)') TRIM(fn),'0',DeplCntl%nowstep
        ELSE
            WRITE(fn,'(A,i3)') TRIM(fn),DeplCntl%nowstep
        ENDIF
    ENDIF        
    
    !--- Old Source 17.12.05. BYS_EDIT
    !IF( lPin )THEN        
    !    WRITE(fn,'(A,A)') TRIM(fn),'_p'
    !    IF( ipin .LT. 10 )THEN
    !        WRITE(fn,'(A,A,i1)') TRIM(fn),'00', ipin
    !    ELSEIF( ipin .LT. 100 )THEN
    !        WRITE(fn,'(A,A,i2)') TRIM(fn),'0', ipin
    !    ELSE
    !        WRITE(fn,'(A,i3)') TRIM(fn), ipin
    !    ENDIF
    !ENDIF    
    IF( .NOT. lPin .AND. lPDQ_GCL )THEN
        WRITE(fn2,'(A,A)') TRIM(fn),'.PDQ'
        OPEN(unit=io2, file = fn2, status = 'replace')
        
        WRITE(io2,'(a)') '! 2G_PinFlux'
        WRITE(io2,'(i4,i4)') nx, 360
        DO gidx = 1, ngrp
            WRITE(io2,'(a4,i1)') '! g=', gidx
            DO i = 1, ny
                WRITE(io2,'(t2,1p,100e12.5)') PinPhiVolg(gidx,:,i)
            ENDDO
        ENDDO
        
        WRITE(io2,'(a12)') '! Power_Dist'
        DO iy = 1, ny
            WRITE(io2,'(t2,1p,100e12.5)') pwr(:,iy)
        ENDDO
        WRITE(io2,'(a13)') '! BU-Exposure'
        DO iy = 1, ny
            WRITE(io2,'(t2,1p,100e12.5)') buexp(:,iy)
        ENDDO
        WRITE(io2,'(a)') '! 2G_ADF (SWNE) '
        !avgADF2g=1 !adf off
        !WRITE(io2,'(t2,100f12.5)') avgADF2g(1), avgADF2g(2)
        DO gidx = 1,ngrp
            WRITE(io2,'(t2,i3,100f12.5)') gidx, ADF2g(1:4,gidx) !--- 14/10/06 edit for non-symmetric assembly : SWNE
        ENDDO
        
        WRITE(io2,'(a)') '! 2G_SurFlx'
        avgBPhiG(:)=PhiG(:)
        WRITE(io2,'(t2,1p,100e12.5)') avgBPhiG(:)
        WRITE(io2,'(a)') '! 2G_ConFlx'
        DO i = 1, 3
            conPhiG(i,:)=PhiG(:)
        ENDDO
        DO gidx = 1, ngrp
            WRITE(io2,'(t2,1p,100e12.5)') conPhiG(:,gidx)
        ENDDO
   
        WRITE(io2,'(a)') '! 2G_flux'
        WRITE(io2,'(t2,1p,100e12.5)') PhiG(:) !fphi, tphi    
        
        !D, a, r, f, nf, kf
        !ALLOCATE(isoMacXs(nisotot,6,ng), isoMacSm(nisotot,ng,ng),isoMacXs2g(nisotot,6,2)) 
        !tr, a, r, f, nu, k
        !ALLOCATE(isoMicXs(nisotot,6,ng), isoMicSm(nisotot,ng,ng),isoMicXs2g(nisotot,6,2))
        WRITE(io2,'(a)') '! 2G_MAC'
        DO gidx = 1, ngrp
            WRITE(io2,'(t2,1p,100e12.5)') isoMacXs2g(0,1:6,gidx), isoMacXs2g(0,0,gidx), isoMacXs2g(0,3,gidx)-isoMacXs2g(0,2,gidx)
        ENDDO
        
        WRITE(io2,'(a)') '! k_inf'
        WRITE(io2,'(t2,100f12.9)') kinf
        WRITE(io2,'(a)') '! Msquare'
        WRITE(io2,'(t2,1p,100e12.5)') Msq
        WRITE(io2,'(a)') '! Bsquare'
        WRITE(io2,'(t2,1p,100e12.5)') Bsq
        
        WRITE(io2,'(a)') '! 2G_Mic'
        WRITE(io2,'(i4)') nisotot
        DO iso = 1, nisotot
            WRITE(io2,'(i5,1p,e12.5)') isoname(iso),isonumden(iso)
            DO gidx = 1, ngrp
                WRITE(io2,'(t2,1p,100e12.5)') isoMicXs2g(iso,1:6,gidx),isoMicXs2g(iso,0,gidx),isoMicXs2g(iso,3,gidx)-isoMicXs2g(iso,2,gidx) 
            ENDDO
        ENDDO
        
        WRITE(io2,'(a)') '.'  !---END of official PDQ EDIT LINE -----------------------------------------------
        CLOSE(io2)
    ENDIF
    
    !--- AGL library --- old FBX
    IF( .NOT.lPin .AND. lASY_GCL )THEN
        WRITE(fn2,'(A,A)') TRIM(fn),'.AGL'
        OPEN(unit=io2, file = fn2, status = 'replace')
        
        
        WRITE(io2,'(t2,100f15.7)') kinf
        DO gidx = 1, ngrp
            WRITE(io2,'(t2,i15,100f15.7)') gidx, ADF2g(1:4,gidx) !--- 14/10/06 edit for non-symmetric assembly : SWNE
        ENDDO
        
        WRITE(io2,'(t2,1p,100e15.7)') Bsq
        
        WRITE(io2,'(t2,1p,100e15.7)') PhiG(:)
        
        DO gidx = 1, ngrp
            WRITE(io2,'(t2,1p,100e15.7)') isoMacXs2g(0,1:6,gidx),isoMacXs2g(0,0,gidx), chig(gidx)
        ENDDO
        !--- 16/03/14 few-group Scat MAtrix
        DO gidx = 1, ngrp
            WRITE(io2,'(t2,1p,100e15.7)') isoMacSm2g(0,gidx,:)
        ENDDO
        !--- 16/03/14 few-group Scat MAtrix end
        
        DO ig = 1, ng
            WRITE(io2,'(t2,1p,100e15.7)') Phi(ig)
        ENDDO
        
        !--- mG XS
        DO ig = 1, ng
            WRITE(io2,'(t2,1p,100e15.7)') isoMacXs(0,1:6,ig),isoMacXs(0,0,ig) !D, a, r, f, nf, kf 
        ENDDO
        DO ig = 1, ng
            WRITE(io2,'(t2,1p,100e15.7)') isoMacSm(0,ig,:)
        ENDDO
        DO i = 0, 2
            DO gidx = 1, ngrp            
                WRITE(io2,'(t2,1p,100e15.7)') SurfJG(i, :, gidx)
            ENDDO
        ENDDO
        !--
        DO gidx = 1, ngrp
            WRITE(io2,'(t2,i15,100f15.5)') gidx, ADF2gSurf(1:4, gidx) !--- 15/06/12 edit for No extrapolation of surf Flux : SWNE
        ENDDO
        DO gidx = 1, ngrp
            WRITE(io2,'(t2,1p,100e15.7)') SurfPhiG(:, gidx)
        ENDDO
        DO ig = 1, ng
            WRITE(io2,'(t2,1p,100e15.7)') isoMacSmP1(0,ig,:)
        ENDDO
        
        DO iy = 1, ny
            WRITE(io2,'(t2,1p,100e15.7)') pwr(:,iy)
        ENDDO
        
        DO gidx = 1, ngrp
            DO i = 1, ny
                WRITE(io2,'(t2,1p,100e15.7)') PinPhiVolg(gidx,:,i)
            ENDDO
        ENDDO
        
        DO iy = 1, ny
            WRITE(io2,'(t2,1p,100e15.7)') buexp(:,iy)
        ENDDO
        
        i = 0
        DO ig = 1, ng            
            WRITE(io2,'(t2,1p,100e15.7)') SurfJ(i, :, ig)
        ENDDO
        
        IF( lASY_MIC )THEN
            !---     Micro-Scopic XS ---
            WRITE(io2,'(i16)') nisotot
            DO iso = 1, nisotot
                IF(isoname(iso).EQ.92238)THEN
                    WRITE(io2,'(i16,1p,3e15.7)') isoname(iso),isonumden(iso),u238n2n/u238phi/isonumden(iso),u238n2n/u238phi/isonumden(iso)/(isoMicXs2g(iso,2,1)*PhiG(1)+isoMicXs2g(iso,2,2)*phig(2))*(phig(1)+phig(2))
                    u238n2n=0
                    u238phi=0
                ELSE
                    WRITE(io2,'(i16,1p,e15.7)') isoname(iso),isonumden(iso)
                ENDIF
                IF(isonumden(iso).NE.0.0_8)THEN
                    DO gidx = 1, ngrp
                        WRITE(io2,'(t2,1p,100e15.7)') isoMicXs2g(iso,1:6,gidx),isoMicXs2g(iso,0,gidx)
                    ENDDO
                ELSE
                    DO gidx = 1, ngrp
                        WRITE(io2,'(t2,1p,100e15.7)') 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8
                    ENDDO
                ENDIF
            ENDDO
!#define mic47Gxs    
#ifdef mic47Gxs    
            WRITE(io2,'(i4)') nisotot
            DO iso = 1, nisotot
                WRITE(io2,'(i15,1p,e15.7)') isoname(iso),isonumden(iso)
                DO ig = 1, ng
                    WRITE(io2,'(t2,1p,100e15.7)') isoMicXs(iso,1:6,ig)
                ENDDO
            ENDDO
#endif    
        ENDIF
        CLOSE(io2)
    ENDIF
    
    !--- FBX library ---
    IF( lPin .AND. lPIN_GCL )THEN
        IF( iPin .EQ. 1 )THEN
            WRITE(fn2,'(A,A)') TRIM(fn),'.PGL'
            OPEN(unit=io2, file = fn2, status = 'replace')
            WRITE(io2,'(t2,100f15.7)') kinf
            WRITE(io2,'(t2,i15)') npinlist
        ENDIF
        WRITE(io2,'(t2,i15)') ipin
        
        WRITE(io2,'(t2,1p,100e15.7)') PhiG(:)
        DO gidx = 1, ngrp
            WRITE(io2,'(t2,1p,100e15.7)') isoMacXs2g(0,1:6,gidx),isoMacXs2g(0,0,gidx), chig(gidx)
        ENDDO
        !--- 16/03/14 few-group Scat MAtrix
        DO gidx = 1, ngrp
            WRITE(io2,'(t2,1p,100e15.7)') isoMacSm2g(0,gidx,:)
        ENDDO
        !--- 16/03/14 few-group Scat MAtrix end
        DO ig = 1, ng
            WRITE(io2,'(t2,1p,100e15.7)') Phi(ig)
        ENDDO
        !--- mG XS
        DO ig = 1, ng
            WRITE(io2,'(t2,1p,100e15.7)') isoMacXs(0,1:6,ig),isoMacXs(0,0,ig) !D, a, r, f, nf, kf 
        ENDDO
        DO ig = 1, ng
            WRITE(io2,'(t2,1p,100e15.7)') isoMacSm(0,ig,:)
        ENDDO
        DO i = 0, 2
            DO gidx = 1, ngrp            
                WRITE(io2,'(t2,1p,100e15.7)') SurfJG(i, :, gidx)
            ENDDO
        ENDDO
        !--
        DO gidx = 1, ngrp
            WRITE(io2,'(t2,i15,100f15.5)') gidx, ADF2gSurf(1:4, gidx) !--- 15/06/12 edit for No extrapolation of surf Flux : SWNE
        ENDDO
        DO gidx = 1, ngrp
            WRITE(io2,'(t2,1p,100e15.7)') SurfPhiG(:, gidx)
        ENDDO
        i = 0
        DO ig = 1, ng            
            WRITE(io2,'(t2,1p,100e15.7)') SurfJ(i, :, ig)
        ENDDO
        
        IF( lPIN_MIC )THEN
            !---     Micro-Scopic XS ---
            WRITE(io2,'(i16)') nisotot
            DO iso = 1, nisotot
                WRITE(io2,'(i16,1p,e15.7)') isoname(iso),isonumden(iso)
                DO gidx = 1, ngrp
                    WRITE(io2,'(t2,1p,100e15.7)') isoMicXs2g(iso,1:6,gidx),isoMicXs2g(iso,0,gidx)
                ENDDO
            ENDDO
#ifdef mic47Gxs    
            WRITE(io2,'(i4)') nisotot
            DO iso = 1, nisotot
                WRITE(io2,'(i15,1p,e15.7)') isoname(iso),isonumden(iso)
                DO ig = 1, ng
                    WRITE(io2,'(t2,1p,100e15.7)') isoMicXs(iso,1:6,ig)
                ENDDO
            ENDDO
#endif    
        ENDIF
        IF( ipin .EQ. npinlist ) CLOSE(io2)
    ENDIF
    
ENDSUBROUTINE


  
SUBROUTINE GCCLR()
    USE GC_mod
    IMPLICIT NONE
    checkiso=.FALSE. ! 0=not Used. 1=used
    nisotot=0; isoNumden=0; 
    h2oNumden=0; uo2Numden=0;
    XSt=0.0; !Sm1=0.0; 
    nCoreFxr = 0
    PhiVol=0; PinPhiVol=0; Phi=0;
    Dng=0;    Bsq=0;  Msq=0
    u238n2n=0.;u238phi=0.0
    DEALLOCATE(isoMicXs, isoMicSm, isoMicXs2g, isoMicSm2g) !tr, a, r, f, nu, k
    DEALLOCATE(isoMacXs, isoMacSm, isoMacXs2g, isoMacSm2g) !D, a, r, f, nf, kf
    DEALLOCATE(isoMacSmP1, isoMacSm2gP1, isoMicSmP1, isoMicSm2gP1)
    IF( .NOT. lPin )THEN
    DEALLOCATE(PinPhiVolG)
    ENDIF
ENDSUBROUTINE
