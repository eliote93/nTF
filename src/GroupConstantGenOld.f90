#include <defines.h>
SUBROUTINE GroupConstGenMac(Core, FmInfo, CmInfo, GroupInfo, nTracerCntl, PE, ng)
USE PARAM
USE TYPEDEF,         ONLY : CoreInfo_Type,FmInfo_Type, CMInfo_Type, PinXs_Type , GroupInfo_Type,PE_TYPE
USE BasicOperation, ONLY : CP_VA, CP_CA, MULTI_VA, MULTI_CA
USE files,           ONLY : caseid
USE XSLIB_MOD
USE CNTL,         ONLY : nTracerCntl_Type
USE CritSpec_mod,     ONLY : GetDiffusionCoeff
USE GroupConst_Mod
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(CMInfo_Type) :: CMInfo
TYPE(FMInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
INTEGER, INTENT(IN) :: ng

REAL :: cellbphi(4,core%nxy,ng), cellbJ(4,core%nxy,ng), cellavgphi(core%nxy,ng)
REAL :: cellbphi2g(4,core%nxy,2), cellbJ2g(4,core%nxy,2), cellavgphi2g(core%nxy,2)
LOGICAL :: noString

TYPE(PinXS_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: PhiC(:, :, :)
REAL :: Phi(ng,Core%nxy), Phi2g(2,Core%nxy)
REAL :: Bsq, Msq, kinf, keff
REAL :: f_ratio(ng)
REAL :: phicrit(ng), Dng(ng)
REAL :: critphisum, infphisum, Rcritphisum, Rinfphisum

INTEGER :: npins, pidx(2,100), gidx
REAL :: ADF(4,ng), ADFden,ADF2g(4,2), phiavg(2), Rnxy
REAL :: avgADF2g(2), avgbphi2g(2)
REAL :: conphi2g(3,2) !corner, center, avg/ ig

INTEGER :: i, j, k
INTEGER :: ig, ig2, iz, ixy, igb, ige
INTEGER :: myzb, myze, nxy, nx, ny
INTEGER :: io, io2, io3

REAL :: RR(ng,0:ng+5), scat(ng,ng)
REAL :: phisum(ng), psisum(ng), volsum, Rphisum, Rpsisum 
REAL :: localphi, localphi2, localpsi, localpsi2, totalpsi, totalphi
CHARACTER(256) :: fn,fn2,fn3

REAL :: FRR(0:7),TRR(0:7), ORR(6)
REAL :: fphi,tphi,tpsi,fpsi

REAL :: EStruct(ng+1),Eavg(ng)
REAL :: ELB, temp
INTEGER :: EIdx

!----------INITIALIZE-----------
iz=1
PinXS => CMInfo%PinXS; PHIC => CmInfo%PhiC
nxy = Core%nxy; nx = Core%nx; ny = Core%ny;
SELECT CASE(nTracerCntl%gc_spec-3) !--- flux spectrum selection
CASE(0) !---critical spectrum
    CALL GetDiffusionCoeff(Core, FmInfo, GroupInfo, nTracerCntl, PE, Dng, phicrit, Bsq, Msq, kinf, keff)    
    critphisum=0; infphisum=0;
    DO ig = 1, ng
        f_ratio(ig)=0
        DO ixy = 1, nxy
            localPhi = Core%PinVOl(ixy, iz) * PhiC(ixy, iz, ig)
            f_ratio(ig)=f_ratio(ig)+localphi
        ENDDO
        infphisum=infphisum+f_ratio(ig)
        critphisum=critphisum+phicrit(ig)
    ENDDO
    Rinfphisum=one/infphisum
    Rcritphisum=one/critphisum
    CALL MULTI_CA(Rinfphisum, f_ratio(:), ng)
    CALL MULTI_CA(Rcritphisum, phicrit(:), ng)
    DO ig = 1, ng
        f_ratio(ig)=phicrit(ig)/f_ratio(ig)
        DO ixy = 1, nxy
            Phi(ig, ixy)=PhiC(ixy, iz, ig)*f_ratio(ig)
        ENDDO
    ENDDO
CASE(1) !---Infinit medium spectrum
    CALL GetDiffusionCoeff(Core, FmInfo, GroupInfo, nTracerCntl, PE, Dng, phicrit, Bsq, Msq, kinf, keff)    
    
    DO ig = 1, ng
        DO ixy = 1, nxy
            Phi(ig,ixy) = PhiC(ixy, iz, ig)
        ENDDO
    ENDDO
CASE(2) !---Infinit medium spectrum (no critical)
    Dng=0
    Bsq=0
    DO ig = 1, ng
        DO ixy = 1, nxy
            Phi(ig,ixy) = PhiC(ixy, iz, ig)
        ENDDO
    ENDDO
ENDSELECT

CONTINUE


!----------HOMOGENIZATION-------------------------------------------------------

volsum=0
totalpsi=0
DO ixy = 1, nxy
    volsum = volsum + Core%PinVol(ixy,iz)
ENDDO
DO ig = 1, ng
    CALL CP_CA(RR(ig,0:5+ng), ZERO, ng + 6)
ENDDO
DO ig = 1, ng
    phisum(ig)=0
    psisum(ig)=0
    DO ixy = 1, nxy
        localPhi = Core%PinVOl(ixy, iz) * Phi(ig,ixy)
        localpsi = localphi * PinXS(ixy,iz)%XSNF(ig)          !-- nsig_F * phi_g*V
        phisum(ig) = phisum(ig) + localphi      !--- g-th group flux that integrated over the total area 
        psisum(ig) = psisum(ig) + localpsi      !--- volume integrated fission source from g-th group
        RR(ig,0) = RR(ig,0) + localphi * PinXS(ixy,iz)%XST(ig)
        RR(ig,1) = RR(ig,1) + localphi * PinXS(ixy,iz)%XSTR(ig)
        RR(ig,2) = RR(ig,2) + localphi * PinXS(ixy,iz)%XSR(ig)
        RR(ig,3) = RR(ig,3) + localphi * PinXS(ixy,iz)%XSA(ig)
        RR(ig,4) = RR(ig,4) + localphi * PinXS(ixy,iz)%XSNF(ig)
        igb = PinXS(ixy,iz)%XSS(ig)%IB
        ige = PinXS(ixy,iz)%XSS(ig)%IE
        DO ig2 = igb, ige
        !DO ig2 = 1, ng
            localPhi2 = Core%PinVOl(ixy, iz) * Phi(ig2,ixy)
            RR(ig2,4 + ig) = RR(ig2,4 + ig) + localphi2 * PinXS(ixy,iz)%XSS(ig)%FROM(ig2)
        ENDDO
        RR(ig,4+ig)=RR(ig,4+ig)+localphi*PinXS(ixy,iz)%XSS(ig)%Withingroupscat
        localpsi2=0
        DO ig2 = 1, ng
            localpsi2=localpsi2+Phi(ig2,ixy) * PinXS(ixy,iz)%XSNF(ig2) * Core%PinVol(ixy,iz)             !-- nsig_F * phi_g*V
        ENDDO
        RR(ig,5+ng) = RR(ig,5+ng) + localpsi2 * PinXS(ixy,iz)%CHI(ig)
    ENDDO
ENDDO
DO ig = 1, ng
    IF( phisum(ig) .EQ. zero )THEN
        Rphisum=0
    ELSE
        Rphisum = one/phisum(ig)
    ENDIF
    CALL MULTI_CA(Rphisum, RR(ig,0:4+ng), 5+ng)
    totalpsi=totalpsi + psisum(ig)
ENDDO

IF (totalpsi .EQ. zero) THEN
    Rpsisum=0
ELSE
    Rpsisum = one/totalpsi
ENDIF
DO ig = 1, ng
    RR(ig,5+ng) = RR(ig,5+ng)*Rpsisum
ENDDO
!DO ig = 1,ng
!    IF (RR(ig,3) .eq. zero) THEN
!        RR(ig,3)=RR(ig,2)
!        DO ig2=1,ng
!            IF( ig2 .NE. ig )THEN
!            RR(ig,3)=RR(ig,3)-RR(ig,4+ig2) !--sig_a
!            ENDIF
!        ENDDO
!    ENDIF
!ENDDO
DO ig = 1,ng
    RR(ig,3)=RR(ig,0)
    DO ig2=1,ng
        IF( ig2 .NE. ig )THEN
        RR(ig,3)=RR(ig,3)-RR(ig,4+ig2) !--sig_a
        ENDIF
    ENDDO
ENDDO
!----------Assembly Discontinuity Factor
#ifdef adf
DO ig = 1, ng
    ADFden=0
    DO ixy = 1, nxy
        ADFden = ADFden + cellphibdry(1,ixy,2,ig)
        Cellavgphi(ixy,ig) = cellphibdry(1,ixy,2,ig)
        DO i = 1, 4
            CellbPhi(i,ixy,ig) = cellphibdry(i,ixy,1,ig)
            CellbJ(i,ixy,ig) = cellphibdry(i,ixy,3,ig)
        ENDDO
    ENDDO
    ADFden=ADFden/nxy
    DO i = 1, 4
        ADF(i, ig)=0
        SELECT CASE(i)
        CASE(1) !---SOUTH
            npins=nx
            DO j = 1, npins
                pidx(1,j)=nxy-nx+j
                pidx(2,j)=pidx(1,j)-nx
            ENDDO
        CASE(2) !---WEST
            npins=ny
            DO j = 1, npins
                pidx(1,j)=(j-1)*nx+1
                pidx(2,j)=pidx(1,j)+1
            ENDDO
        CASE(3) !---NORTH
            npins=nx
            DO j = 1, npins
                pidx(1,j)=j
                pidx(2,j)=pidx(1,j)+nx
            ENDDO
        CASE(4) !---EAST
            npins=ny
            DO j = 1, npins
                pidx(1,j)=j*nx
                pidx(2,j)=pidx(1,j)-1
            ENDDO
        ENDSELECT
        IF(.TRUE.)THEN
            DO j = 1, npins
                ADF(i,ig)=ADF(i,ig) + cellphibdry(i,j,1,ig) !surf/pin/1/g
            ENDDO
            ADF(i,ig)=ADF(i,ig)/(ADFden*npins)
        ELSE
            DO j = 1, npins
                ADF(i, ig) = ADF(i, ig) + ( 3*Phi(ig,pidx(1,j)) - Phi(ig,pidx(2,j)) )*0.5            
            ENDDO
            ADF(i, ig)=ADF(i, ig)*volsum/(npins*phisum(ig))
        ENDIF
    ENDDO    
ENDDO
#endif
io=42
WRITE(fn,'(A,A,A)') 'H_',TRIM(caseid),'.xslib'
OPEN(unit=io, file = fn, status = 'replace')
WRITE(io,'(a)') ' base_micro 1'
DO ig = 1, ng
    IF( Dng(ig) .NE. 0 )THEN
        WRITE(io,'(t2,1p,7e20.12)') RR(ig,1),RR(ig,3),RR(ig,4),RR(ig,4),RR(ig,5+ng), Dng(ig), RR(ig,0)
    ELSE
        WRITE(io,'(t2,1p,7e20.12)') RR(ig,1),RR(ig,3),RR(ig,4),RR(ig,4),RR(ig,5+ng), 1/3/RR(ig,1), RR(ig,0)
    ENDIF        
ENDDO
DO ig = 1, ng
    WRITE(io,'(t2,1p,100e20.12)') RR(ig,5:4+ng) 
ENDDO
DO ig = 1, ng
    WRITE(io,'(t2,1p,100e20.12)') phisum(ig)
ENDDO
io=43
WRITE(fn,'(A,A,A)') 'HET_',TRIM(caseid),'.xslib'
OPEN(unit=io, file = fn, status = 'replace')
DO ixy = 1, nxy
    scat=zero
    DO ig = 1, ng
        igb = PinXS(ixy,iz)%XSS(ig)%IB
        ige = PinXS(ixy,iz)%XSS(ig)%IE
        scat(ig,ig) = PinXS(ixy,iz)%XSS(ig)%WITHINGROUPSCAT
        DO ig2 = igb, ige
            scat(ig2,ig) = scat(ig2,ig) + PinXS(ixy,iz)%XSS(ig)%FROM(ig2)        
        ENDDO
    ENDDO
    !WRITE(io,'(a12,i3)') ' base_micro ', ixy
    WRITE(io,'(i3)') ixy
    DO ig = 1, ng
        WRITE(io,'(t2,1p,5e20.12)') PinXS(ixy,iz)%XSTR(ig),PinXS(ixy,iz)%XSA(ig),PinXS(ixy,iz)%XSNF(ig),PinXS(ixy,iz)%XSNF(ig),PinXS(ixy,iz)%CHI(ig)
    ENDDO
    DO ig = 1, ng
        WRITE(io,'(t2,1p,100e20.12)') scat(ig,:)
    ENDDO
ENDDO
DO ig = 1, ng
    WRITE(io,'(t2,1p,1000e20.12)') Phi(ig,:)
ENDDO

!-------------- 2group condensing-----------------------------------------------
IF (ng .eq. 47) THEN
EStruct(1:ng) = enbhel(1:ng)
EStruct(ng+1) = 1.0E-4_8
DO ig = 1, ng
  Eavg(ig) = (EStruct(ig) + EStruct(ig+1))/2._8
ENDDO
ELB = 6.2506E-01
DO ig = 1, ng-1
    temp = (ELB-Eavg(ig))*(ELB-Eavg(ig+1))
    IF(temp .LE. 0._8) THEN
        EIdx=ig !---1~EIdx-1 : fast , EIdx~ng : thermal / EIdx=35 / 1~34, 35~47
        EXIT
    ENDIF
ENDDO
    EIdx=EIdx+1 !---old bug
ELSEIF (ng .eq. 7) THEN
    EIdx=4
ELSEIF (ng .eq. 2) THEN
    EIdx=2
ENDIF


CALL CP_CA(FRR(:), ZERO, 8)
CALL CP_CA(TRR(:), ZERO, 8)
fphi=0; tphi=0;

DO ig = 1, ng    
    IF (ig .LT. EIdx) THEN
        fphi=fphi+phisum(ig)
        FRR(0) = FRR(0) + RR(ig,0) * phisum(ig)  ! Sig_t
        FRR(1) = FRR(1) + RR(ig,1) * phisum(ig)  ! Sig_tr
        FRR(2) = FRR(2) + RR(ig,3) * phisum(ig)  ! Sig_a
        FRR(3) = FRR(3) + RR(ig,4) * phisum(ig)  ! nSig_f
        FRR(4) = FRR(4) + 1/Dng(ig)  * phisum(ig)  ! Diff_coeff
        !---CRR(:,5~6) = scattering XS
        DO ig2 = 1, ng
            IF (ig2 .LT. Eidx) THEN
                FRR(5) = FRR(5) + RR(ig,4+ig2) * phisum(ig)
            ELSE
                FRR(6) = FRR(6) + RR(ig,4+ig2) * phisum(ig)
            ENDIF
        ENDDO
        FRR(7) = FRR(7) + RR(ig,ng+5)  ! X (chi)
    ELSE
        tphi=tphi+phisum(ig)
        TRR(0) = TRR(0) + RR(ig,0) * phisum(ig)  ! Sig_t
        TRR(1) = TRR(1) + RR(ig,1) * phisum(ig)  ! Sig_tr
        TRR(2) = TRR(2) + RR(ig,3) * phisum(ig)  ! Sig_a
        TRR(3) = TRR(3) + RR(ig,4) * phisum(ig)  ! nSig_f
        TRR(4) = TRR(4) + 1/Dng(ig)  * phisum(ig)  ! kSig_f
        !---CRR(:,5~6) = scattering XS
        DO ig2 = 1, ng
            IF (ig2 .LT. EIdx) THEN
                TRR(5) = TRR(5) + RR(ig,4+ig2) * phisum(ig)
            ELSE
                TRR(6) = TRR(6) + RR(ig,4+ig2) * phisum(ig)
            ENDIF
        ENDDO
        TRR(7) = TRR(7) + RR(ig,ng+5)  ! X (chi)
    ENDIF
ENDDO
Rphisum = one/fphi
CALL MULTI_CA(Rphisum, FRR(0:6), 7)
Rphisum = one/tphi
CALL MULTI_CA(Rphisum, TRR(0:6), 7)

!----------Assembly Discontinuity Factor in 2-Groups
phiavg=0
phi2g=0
Cellavgphi2g=0
CellbPhi2g=0
CellbJ2g=0
#ifdef adf
DO ig = 1, ng
    IF (ig .LT. EIdx) THEN
        gidx=1
    ELSE
        gidx=2
    ENDIF
    DO ixy = 1, nxy
        IF( .TRUE. )THEN
            phiavg(gidx) = phiavg(gidx) + cellphibdry(1,ixy,2,ig)
        ELSE
            phiavg(gidx)=phiavg(gidx) + Phi(ig,ixy)
        ENDIF
        phi2g(gidx,ixy) = phi2g(gidx,ixy) + Phi(ig,ixy)
        !--- ADF verfication work
        Cellavgphi2g(ixy,gidx) = Cellavgphi2g(ixy,gidx) + Cellavgphi(ixy,ig)
        DO i = 1, 4
            CellbPhi2g(i,ixy,gidx) = CellbPhi2g(i,ixy,gidx) + CellbPhi(i,ixy,ig)
            CellbJ2g(i,ixy,gidx) = CellbJ2g(i,ixy,gidx) + CellbJ(i,ixy,ig)
        ENDDO
    ENDDO
ENDDO

Rnxy=one/nxy
CALL MULTI_CA(Rnxy, phiavg(1:2), 2)
avgADF2g=0
DO i = 1, 4
    SELECT CASE(i)
    CASE(1) !---SOUTH
        npins=nx
        DO j = 1, npins
            pidx(1,j)=nxy-nx+j
            pidx(2,j)=pidx(1,j)-nx
        ENDDO
    CASE(2) !---WEST
        npins=ny
        DO j = 1, npins
            pidx(1,j)=(j-1)*nx+1
            pidx(2,j)=pidx(1,j)+1
        ENDDO
    CASE(3) !---NORTH
        npins=nx
        DO j = 1, npins
            pidx(1,j)=j
            pidx(2,j)=pidx(1,j)+nx
        ENDDO
    CASE(4) !---EAST
        npins=ny
        DO j = 1, npins
            pidx(1,j)=j*nx
            pidx(2,j)=pidx(1,j)-1
        ENDDO
    ENDSELECT        
    DO ig = 1, ng
        IF (ig .LT. EIdx) THEN
            gidx=1
        ELSE
            gidx=2
        ENDIF
        DO j = 1, npins
            IF( .TRUE. )THEN
                ADF2g(i, gidx) = ADF2g(i, gidx) + cellphibdry(i,pidx(1,j),1,ig)
            ELSE
                ADF2g(i, gidx) = ADF2g(i, gidx) + ( 3*Phi(ig,pidx(1,j)) - Phi(ig,pidx(2,j)) )*0.5            
            ENDIF
        ENDDO
    ENDDO    
    DO ig = 1, 2
        avgbphi2g(ig)=avgbphi2g(ig)+ADF2g(i, ig)/(4*npins)
        ADF2g(i, ig)=ADF2g(i, ig)/(npins*phiavg(ig))
        avgADF2g(ig)=avgADF2g(ig)+ADF2g(i,ig)*0.25
    ENDDO
ENDDO

!Corner flux

DO ig = 1, 2
    conphi2g(1,ig)=cellphibdry(3,1,1,ig)
    conphi2g(2,ig)=cellphibdry(3,nx/2,1,ig)
    conphi2g(3,ig)=(conphi2g(1,ig)+conphi2g(2,ig))/2
ENDDO
#endif
    


io2=44
noString=.TRUE.
noString=.FALSE.
WRITE(fn2,'(A,A,A)') 'H_',TRIM(caseid),'_2G_MAC.xslib'
OPEN(unit=io2, file = fn2, status = 'replace')
IF(.NOT. noString) WRITE(io2,'(a)') ' MacroXS Sig_t/sig_a/nusig_f/ksig_f/Chi/D /total'
WRITE(io2,'(t2,1p,6e20.12)') FRR(1),FRR(2),FRR(3),FRR(3),FRR(7),1/FRR(4), FRR(0)
WRITE(io2,'(t2,1p,6e20.12)') TRR(1),TRR(2),TRR(3),TRR(3),TRR(7),1/TRR(4), TRR(0)
WRITE(io2,'(t2,1p,100e20.12)') FRR(5), FRR(6) 
WRITE(io2,'(t2,1p,100e20.12)') TRR(5), TRR(6) 
IF(.NOT. noString) WRITE(io2,'(a)') ' Fast/ Thermal flux'
WRITE(io2,'(t2,1p,100e20.12)') fphi, tphi
IF(.NOT. noString) WRITE(io2,'(a)') ' 2G_ADF Fast/ Thermal'
WRITE(io2,'(t2,1p,100e20.12)') avgADF2g(1), avgADF2g(2)
IF(.NOT. noString) WRITE(io2,'(a)') ' 2G_SurFlx Fast/ Thermal'
WRITE(io2,'(t2,1p,100e20.12)') avgbphi2g(1), avgbphi2g(2)
IF(.NOT. noString) WRITE(io2,'(a)') ' 2G_ConFlx Fast(con-con-avg)/ Thermal(con-con-avg)'
WRITE(io2,'(t2,1p,100e20.12)') conphi2g(1,1), conphi2g(2,1), conphi2g(3,1)
WRITE(io2,'(t2,1p,100e20.12)') conphi2g(1,2), conphi2g(2,2), conphi2g(3,2)

!IF(.NOT. noString) WRITE(io2,'(a)') ' 2G_ADF : SF,T  WF,T  NF,T  EF,T '
!WRITE(io2,'(t2,1p,100e20.12)') ADF2g(1,1), ADF2g(1,2), ADF2g(2,1), ADF2g(2,2), ADF2g(3,1), ADF2g(3,2), ADF2g(4,1), ADF2g(4,2)
!WRITE(io2,'(a)') ' Bsquare '
!WRITE(io2,'(t2,1p,100e20.12)') Bsq
WRITE(io2,'(a)') ' '
IF(.NOT. noString) WRITE(io2,'(a)') ' 2group flux distribution'
DO ig = 1, 2
    IF(.NOT. noString) WRITE(io2,'(a4,i1)') 'g = ', ig
    DO i = 1, ny
        j=(i-1)*nx+1
        WRITE(io2,'(t2,1p,100e20.12)') cellavgphi2g(j:j+nx-1,ig)
    ENDDO
ENDDO
WRITE(io2,'(a)') ' '
IF(.NOT. noString) WRITE(io2,'(a)') ' 2group boundary flux distribution '
DO k = 1, 4 !direction
    IF(.NOT. noString) WRITE(io2,'(a4,i1)') 'd = ', k
    DO ig = 1, 2
        IF(.NOT. noString) WRITE(io2,'(a6,i1)') '  g = ', ig
        DO i = 1, ny
            j=(i-1)*nx+1
            WRITE(io2,'(t2,1p,100e20.12)') cellbphi2g(k,j:j+nx-1,ig)
        ENDDO
    ENDDO
ENDDO
WRITE(io2,'(a)') ' '
IF(.NOT. noString) WRITE(io2,'(a)') ' 2group boundary current distribution '
DO k = 1, 4 !direction
    IF(.NOT. noString) WRITE(io2,'(a4,i1)') 'd = ', k
    DO ig = 1, 2
        IF(.NOT. noString) WRITE(io2,'(a6,i1)') '  g = ', ig
        DO i = 1, ny
            j=(i-1)*nx+1
            WRITE(io2,'(t2,1p,100e20.12)') cellbJ2g(k,j:j+nx-1,ig)
        ENDDO
    ENDDO
ENDDO
WRITE(io2,'(a)') ' '
IF(.NOT. noString) WRITE(io2,'(a)') ' 2group flux distribution (normalized)'
DO ig = 1, 2
    IF(.NOT. noString) WRITE(io2,'(a4,i1)') 'g = ', ig
    DO i = 1, ny
        j=(i-1)*nx+1
        WRITE(io2,'(t2,1p,100e20.12)') phi2g(ig,j:j+nx-1)
    ENDDO
ENDDO

#ifdef nostr
io2=45
noString=.TRUE.
WRITE(fn2,'(A,A,A)') 'H_',TRIM(caseid),'_2G_noStr.xslib'
OPEN(unit=io2, file = fn2, status = 'replace')
IF(.NOT. noString) WRITE(io2,'(a)') ' base_micro 1'
WRITE(io2,'(t2,1p,6e20.12)') FRR(1),FRR(2),FRR(3),FRR(3),FRR(7),FRR(4)
WRITE(io2,'(t2,1p,6e20.12)') TRR(1),TRR(2),TRR(3),TRR(3),TRR(7),TRR(4)
WRITE(io2,'(t2,1p,100e20.12)') FRR(5), FRR(6) 
WRITE(io2,'(t2,1p,100e20.12)') TRR(5), TRR(6) 
IF(.NOT. noString) WRITE(io2,'(a)') ' 2G_ADF : SF,T  WF,T  NF,T  EF,T '
WRITE(io2,'(t2,1p,100e20.12)') ADF2g(1,1), ADF2g(1,2), ADF2g(2,1), ADF2g(2,2), ADF2g(3,1), ADF2g(3,2), ADF2g(4,1), ADF2g(4,2)
!WRITE(io2,'(a)') ' Bsquare '
!WRITE(io2,'(t2,1p,100e20.12)') Bsq
WRITE(io2,'(a)') ' '
IF(.NOT. noString) WRITE(io2,'(a)') ' 2group flux distribution'
DO ig = 1, 2
    IF(.NOT. noString) WRITE(io2,'(a4,i1)') 'g = ', ig
    DO i = 1, ny
        j=(i-1)*nx+1
        WRITE(io2,'(t2,1p,100e20.12)') cellavgphi2g(j:j+nx-1,ig)
    ENDDO
ENDDO
WRITE(io2,'(a)') ' '
IF(.NOT. noString) WRITE(io2,'(a)') ' 2group boundary flux distribution '
DO k = 1, 4 !direction
    IF(.NOT. noString) WRITE(io2,'(a4,i1)') 'd = ', k
    DO ig = 1, 2
        IF(.NOT. noString) WRITE(io2,'(a6,i1)') '  g = ', ig
        DO i = 1, ny
            j=(i-1)*nx+1
            WRITE(io2,'(t2,1p,100e20.12)') cellbphi2g(k,j:j+nx-1,ig)
        ENDDO
    ENDDO
ENDDO
WRITE(io2,'(a)') ' '
IF(.NOT. noString) WRITE(io2,'(a)') ' 2group boundary current distribution '
DO k = 1, 4 !direction
    IF(.NOT. noString) WRITE(io2,'(a4,i1)') 'd = ', k
    DO ig = 1, 2
        IF(.NOT. noString) WRITE(io2,'(a6,i1)') '  g = ', ig
        DO i = 1, ny
            j=(i-1)*nx+1
            WRITE(io2,'(t2,1p,100e20.12)') cellbJ2g(k,j:j+nx-1,ig)
        ENDDO
    ENDDO
ENDDO
WRITE(io2,'(a)') ' '
IF(.NOT. noString) WRITE(io2,'(a)') ' 2group flux distribution (normalized)'
DO ig = 1, 2
    IF(.NOT. noString) WRITE(io2,'(a4,i1)') 'g = ', ig
    DO i = 1, ny
        j=(i-1)*nx+1
        WRITE(io2,'(t2,1p,100e20.12)') phi2g(ig,j:j+nx-1)
    ENDDO
ENDDO



!-------------- 1 group condensing-----------------------------------------------
IF( 1 .EQ. 2 )THEN
CALL CP_CA(ORR(:), ZERO, 6)
totalphi=tphi+fphi
ORR(1)=FRR(1)*fphi+TRR(1)*tphi
ORR(2)=FRR(2)*fphi+TRR(2)*tphi
ORR(3)=FRR(3)*fphi+TRR(3)*tphi
ORR(4)=FRR(4)*fphi+TRR(4)*tphi
ORR(5)=(FRR(5)+FRR(6))*fphi+(TRR(5)+TRR(6))*tphi
ORR(6)=FRR(7)+TRR(7)

Rphisum = one/totalphi
CALL MULTI_CA(Rphisum, ORR(1:5), 5)

io3=46

WRITE(fn3,'(A,A,A)') 'H_',TRIM(caseid),'_1G.xslib'
OPEN(unit=io3, file = fn3, status = 'replace')

WRITE(io3,'(a)') ' base_micro 1'
WRITE(io3,'(t2,1p,6e20.12)') ORR(1),ORR(2),ORR(3),ORR(3),ORR(6),ORR(4)
WRITE(io3,'(t2,1p,100e20.12)') ORR(5)
!WRITE(io3,'(t2,1p,2e20.12)') fphi/tphi, tphi/tphi 
ENDIF
#endif

END SUBROUTINE
