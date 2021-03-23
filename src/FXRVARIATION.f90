SUBROUTINE Branch_driver
    USE PARAM
    USE Geom,           ONLY : Core
    USE Core_Mod,       ONLY : FmInfo
    USE PE_MOD,         ONLY : PE       
    USE FILES,         ONLY : io8
    USE IOUTIL,        ONLY : message,           ShowHbar1,            ShowHbar2
    USE CNTL,          ONLY : nTracerCntl
    USE FXRVAR_MOD
    IMPLICIT NONE
    
    INTEGER :: icase
    LOGICAL :: Master, Slave
    LOGICAL, SAVE :: lFirst
    Data lFirst /.TRUE./
    Master = PE%Master; Slave = PE%Slave
    
    CurrentVarId=branch_base
    CurrentCaseId=0
    !--- Initialize Coolant density
    CALL FXRVARIATION_BRANCH(Core%ncorefxr,FmInfo,PE,branch_tmod,0,.FALSE.)
    nTracerCntl%lInitBoron = .TRUE.
    !--- Base condition   
    CALL SSEIG()
    
    nTracerCntl%CaseNo = nTracerCntl%CaseNo + 1
    IF( lTmod )THEN
        IF(MASTER) CALL ShowHbar2(io8)
        WRITE(mesg, '(A)') 'Performing T_mod variation ...'
        IF(Master) CALL message(io8, TRUE, TRUE, mesg)
        CurrentVarId=branch_tmod
        DO icase = 1, ntmod
            CurrentCaseId=icase
            IF(Master) CALL ShowHbar1(io8)
            WRITE(mesg, '(A,i2,a,i2)') 'T_mod variation ... CaseId : ',icase, ' /', ntmod
            IF(Master) CALL message(io8, TRUE, TRUE, mesg)
            CALL FXRVARIATION_BRANCH(Core%ncorefxr,FmInfo,PE,branch_tmod,icase,.FALSE.)
            CALL SSEIG()
            CALL FXRVARIATION_BRANCH(Core%ncorefxr,FmInfo,PE,branch_tmod,icase,.TRUE.)
        ENDDO
        nTracerCntl%CaseNo = nTracerCntl%CaseNo + 1
    ENDIF
    IF( lTFuel )THEN
        IF(MASTER) CALL ShowHbar2(io8)
        WRITE(mesg, '(A)') 'Performing T_fuel variation ...'
        IF(Master) CALL message(io8, TRUE, TRUE, mesg)
        CurrentVarId=branch_tfuel
        DO icase = 1, ntfuel
            CurrentCaseId=icase
            IF(Master) CALL ShowHbar1(io8)
            WRITE(mesg, '(A,i2,a,i2)') 'T_fuel variation ... CaseId : ',icase, ' /', ntfuel
            IF(Master) CALL message(io8, TRUE, TRUE, mesg)
            CALL FXRVARIATION_BRANCH(Core%ncorefxr,FmInfo,PE,branch_tfuel,icase,.FALSE.)
            CALL SSEIG()
            CALL FXRVARIATION_BRANCH(Core%ncorefxr,FmInfo,PE,branch_tfuel,icase,.TRUE.)
        ENDDO
        nTracerCntl%CaseNo = nTracerCntl%CaseNo + 1
    ENDIF
    IF( lboron )THEN
        IF(MASTER) CALL ShowHbar2(io8)
        WRITE(mesg, '(A)') 'Performing Boron variation ...'
        IF(Master) CALL message(io8, TRUE, TRUE, mesg)
        CurrentVarId=branch_boron
        bboron(0)=nTracerCntl%BoronPPM
        DO icase = 1, nboron
            CurrentCaseId=icase
            IF(Master) CALL ShowHbar1(io8)
            WRITE(mesg, '(A,i2,a,i2)') 'Boron variation ... CaseId : ',icase, ' /', nboron
            IF(Master) CALL message(io8, TRUE, TRUE, mesg)
            nTracerCntl%BoronPPM=bboron(icase)
            nTracerCntl%lInitBoron = .TRUE.
            CALL SSEIG()
        ENDDO
        nTracerCntl%BoronPPM=bboron(0)
        nTracerCntl%CaseNo = nTracerCntl%CaseNo + 1
    ENDIF
    IF( lRho )THEN
        IF(MASTER) CALL ShowHbar2(io8)
        WRITE(mesg, '(A)') 'Performing Rho_mod variation ...'
        IF(Master) CALL message(io8, TRUE, TRUE, mesg)
        CurrentVarId=branch_rho
        DO icase = 1, nrho
            CurrentCaseId=icase
            IF(Master) CALL ShowHbar1(io8)
            WRITE(mesg, '(A,i2,a,i2)') 'Rho_mod variation ... CaseId : ',icase, ' /', nrho
            IF(Master) CALL message(io8, TRUE, TRUE, mesg)
            CALL FXRVARIATION_BRANCH(Core%ncorefxr,FmInfo,PE,branch_rho,icase,.FALSE.)
            CALL SSEIG()
        ENDDO
        CALL FXRVARIATION_BRANCH(Core%ncorefxr,FmInfo,PE,branch_rho,0,.TRUE.)
        nTracerCntl%CaseNo = nTracerCntl%CaseNo + 1
    ENDIF
    IF(MASTER) CALL ShowHbar2(io8)
    CurrentVarId=branch_base
    CurrentCaseId=0
    
ENDSUBROUTINE
SUBROUTINE FXRVAR_INIT
    USE FXRVAR_MOD
    IMPLICIT NONE
    lFXRvar=.FALSE.
    ltmod=.FALSE.
    ltfuel=.FALSE.
    lrho=.FALSE.
    tmod=0
    tfuel=0
    rho=0
ENDSUBROUTINE
SUBROUTINE FXRVARIATION_BRANCH(nFXR, FmInfo,PE,ivar,icase,lrollback)
    USE PARAM
    USE FXRVAR_MOD
    USE NuclidMap_mod, ONLY : AtomicWeight
    USE SteamTBL_mod, ONLY : steamtbl
    USE TYPEDEF, ONLY : FmInfo_Type, FxrInfo_TYPE, PE_TYPE
    USE CNTL, ONLY : nTracerCntl
    USE FILES,            ONLY : io8
    USE IOUTIL,           ONLY : message
    IMPLICIT NONE
    TYPE(FMInfo_Type) :: FmInfo
    TYPE(FxrInfo_Type), POINTER :: FXR(:,:)
    TYPE(PE_Type) :: PE
    INTEGER :: ivar, icase
    LOGICAL :: lrollback
    INTEGER :: ifxr, nFXR, iz
    INTEGER :: hidx, oidx, idx, chk
    REAL :: rate, aw, ndenh2o
    REAL :: pexit
    REAL :: wt, wh, hin, wrho, wvin, wxin, wbetain, wkapain, wcpin
    
    DO iz= PE%myzb, PE%myze
    FXR=>Fminfo%FXR
    SELECT CASE(ivar)
        CASE(branch_tmod)
            PEXIT = nTracerCntl%PExit
            DO ifxr = 1, nFXR
                IF( FXR(ifxr,iz)%lh2o )THEN
                    IF( .NOT. lRollback )THEN
                        FXR(ifxr,iz)%temp=FXR(ifxr,iz)%temp+tmod(icase)
                    ELSE
                        FXR(ifxr,iz)%temp=FXR(ifxr,iz)%temp-tmod(icase)
                    ENDIF
                    wt=FXR(ifxr,iz)%temp
                    CALL steamtbl(TRUE,pexit,wt,wh,wrho,wvin,wxin,wbetain,wkapain,wcpin)
                    DO idx = 1, FXR(ifxr,iz)%ndim
                        SELECT CASE(FXR(ifxr,iz)%idiso(idx))
                        CASE(1001)
                            hidx=idx
                        CASE(8016)
                            oidx=idx
                        ENDSELECT
                    ENDDO
                    wrho=wrho/1000
                    aw=2*AtomicWeight(1001)+AtomicWeight(8016)      ! 170911 BYS JJH EDIT
                    ndenh2o=wrho* AVOGADRO / aw                     ! *0.03344 !(0.03344=6.02/18)
                    FXR(ifxr,iz)%pnum(oidx)=ndenh2o
                    FXR(ifxr,iz)%pnum(hidx)=2*ndenh2o
                ENDIF
            ENDDO
            IF(PE%MASTER .AND. .NOT. lRollBack)THEN
                WRITE(mesg,'(a,f7.5,a)') '  FXR_VARIATION : Applying Water gram density  = ', wrho, ' g/cc'
                CALL message(io8, TRUE, TRUE, mesg)
            ENDIF
        CASE(branch_tfuel)
            DO ifxr = 1, nFXR
                IF( FXR(ifxr,iz)%lfuel )THEN
                    IF( .NOT. lRollback )THEN
                        FXR(ifxr,iz)%temp=FXR(ifxr,iz)%temp+tfuel(icase)
                    ELSE
                        FXR(ifxr,iz)%temp=FXR(ifxr,iz)%temp-tfuel(icase)
                    ENDIF
                ENDIF
            ENDDO
        CASE(branch_rho)
            IF( .NOT. lRollback )THEN
                DO ifxr = 1, nFXR
                    IF( FXR(ifxr,iz)%lh2o )THEN
                        DO idx = 1, FXR(ifxr,iz)%ndim
                            SELECT CASE(FXR(ifxr,iz)%idiso(idx))
                            CASE(1001)
                                hidx=idx
                            CASE(8016)
                                oidx=idx
                            ENDSELECT
                        ENDDO
                        aw=2*AtomicWeight(1001)+AtomicWeight(8016)      ! 170911 BYS JJH EDIT
                        !rate=rho/FXR(ifxr,iz)%pnum(oidx)* AVOGADRO / aw !*0.03344 !(0.03344=6.02/18)
                        !FXR(ifxr,iz)%pnum(oidx)=FXR(ifxr,iz)%pnum(oidx)*rate
                        !FXR(ifxr,iz)%pnum(hidx)=FXR(ifxr,iz)%pnum(hidx)*rate
                        ndenh2o=rho(icase)* AVOGADRO / aw                     ! *0.03344 !(0.03344=6.02/18)
                        FXR(ifxr,iz)%pnum(oidx)=ndenh2o
                        FXR(ifxr,iz)%pnum(hidx)=2*ndenh2o
                    ENDIF
                ENDDO
            ELSE
                DO ifxr = 1, nFXR
                    IF( FXR(ifxr,iz)%lh2o )THEN
                        wt=FXR(ifxr,iz)%temp
                        CALL steamtbl(TRUE,pexit,wt,wh,wrho,wvin,wxin,wbetain,wkapain,wcpin)
                        DO idx = 1, FXR(ifxr,iz)%ndim
                            SELECT CASE(FXR(ifxr,iz)%idiso(idx))
                            CASE(1001)
                                hidx=idx
                            CASE(8016)
                                oidx=idx
                            ENDSELECT
                        ENDDO
                        wrho=wrho/1000
                        aw=2*AtomicWeight(1001)+AtomicWeight(8016)      ! 170911 BYS JJH EDIT
                        ndenh2o=wrho* AVOGADRO / aw                     ! *0.03344 !(0.03344=6.02/18)
                        FXR(ifxr,iz)%pnum(oidx)=ndenh2o
                        FXR(ifxr,iz)%pnum(hidx)=2*ndenh2o
                    ENDIF
                ENDDO
            ENDIF
    ENDSELECT
    ENDDO
ENDSUBROUTINE
    
SUBROUTINE FXRVARIATION(nFXR, FmInfo,PE)
    USE PARAM
    USE FXRVAR_MOD
    USE NuclidMap_mod, ONLY : AtomicWeight
    USE SteamTBL_mod, ONLY : steamtbl
    USE TYPEDEF, ONLY : FmInfo_Type, FxrInfo_TYPE, PE_TYPE
    USE CNTL, ONLY : nTracerCntl
    USE FILES,            ONLY : io8
    USE IOUTIL,           ONLY : message
    IMPLICIT NONE
    TYPE(FMInfo_Type) :: FmInfo
    TYPE(FxrInfo_Type), POINTER :: FXR(:,:)
    TYPE(PE_Type) :: PE
    INTEGER :: ifxr, nFXR, iz
    INTEGER :: hidx, oidx, idx, chk
    REAL :: rate, aw, ndenh2o
    REAL :: pexit
    REAL :: wt, wh, hin, wrho, wvin, wxin, wbetain, wkapain, wcpin
    
    DO iz= PE%myzb, PE%myze
    FXR=>Fminfo%FXR
    IF(ltfuel)THEN
        DO ifxr = 1, nFXR
            IF( FXR(ifxr,iz)%lfuel )THEN
                FXR(ifxr,iz)%temp=FXR(ifxr,iz)%temp+tfuel(1)
            ENDIF
        ENDDO
    ENDIF
    IF(ltmod)THEN
        PEXIT = nTracerCntl%PExit
        DO ifxr = 1, nFXR
            IF( FXR(ifxr,iz)%lh2o )THEN
                FXR(ifxr,iz)%temp=FXR(ifxr,iz)%temp+tmod(1)
                wt=FXR(ifxr,iz)%temp
                CALL steamtbl(TRUE,pexit,wt,wh,wrho,wvin,wxin,wbetain,wkapain,wcpin)
                DO idx = 1, FXR(ifxr,iz)%ndim
                    SELECT CASE(FXR(ifxr,iz)%idiso(idx))
                    CASE(1001)
                        hidx=idx
                    CASE(8016)
                        oidx=idx
                    ENDSELECT
                ENDDO
                wrho=wrho/1000
                aw=2*AtomicWeight(1001)+AtomicWeight(8016)      ! 170911 BYS JJH EDIT
                ndenh2o=wrho* AVOGADRO / aw                     ! *0.03344 !(0.03344=6.02/18)
                FXR(ifxr,iz)%pnum(oidx)=ndenh2o
                FXR(ifxr,iz)%pnum(hidx)=2*ndenh2o
            ENDIF
        ENDDO
        IF(PE%MASTER) THEN
            WRITE(mesg,'(a,f7.5,a)') '  FXR_VARIATION : Applying Water gram density  = ', wrho, ' g/cc'
            CALL message(io8, TRUE, TRUE, mesg)
        ENDIF
    ENDIF
    IF(lrho)THEN
        DO ifxr = 1, nFXR
            IF( FXR(ifxr,iz)%lh2o )THEN
                DO idx = 1, FXR(ifxr,iz)%ndim
                    SELECT CASE(FXR(ifxr,iz)%idiso(idx))
                    CASE(1001)
                        hidx=idx
                    CASE(8016)
                        oidx=idx
                    ENDSELECT
                ENDDO
                aw=2*AtomicWeight(1001)+AtomicWeight(8016)      ! 170911 BYS JJH EDIT
                !rate=rho/FXR(ifxr,iz)%pnum(oidx)* AVOGADRO / aw !*0.03344 !(0.03344=6.02/18)
                !FXR(ifxr,iz)%pnum(oidx)=FXR(ifxr,iz)%pnum(oidx)*rate
                !FXR(ifxr,iz)%pnum(hidx)=FXR(ifxr,iz)%pnum(hidx)*rate
                ndenh2o=rho(1)* AVOGADRO / aw                     ! *0.03344 !(0.03344=6.02/18)
                FXR(ifxr,iz)%pnum(oidx)=ndenh2o
                FXR(ifxr,iz)%pnum(hidx)=2*ndenh2o
            ENDIF
        ENDDO
    ENDIF
    ENDDO
ENDSUBROUTINE
    
MODULE EFT_MOD
    USE TYPEDEF, ONLY : FuelTH_Type
    INTEGER :: nEFTP=100
    REAL :: P_eft(100)
    REAL :: k_ref(100) 
    REAL :: surfT(100)
    REAL :: meanT(100)    
    REAL :: centT(100) 
    REAL :: rho_mod(100)
    REAL :: T_mod(100)
    REAL :: T_emt(100)    
    INTEGER :: ntest(100)
    
    REAL :: T_eft(100,100)  ! power, test
    REAL :: k_eft(100,100)  ! power, test
    REAL :: e_eft(100,100)  ! power, test
    REAL :: modeff(100)
    REAL :: w1(100,100) ! power, test
    REAL :: w2(100,100) ! power, test
    TYPE(FuelTH_Type), POINTER :: FuelTH0(:)
    LOGICAL :: lEFT_GCGEN
    INTEGER :: iboroncase, ipowercase, calcid
ENDMODULE

SUBROUTINE EFT_driver
    USE PARAM
    USE Geom,           ONLY : Core
    USE Core_Mod,       ONLY : FmInfo, eigv, THinfo, CMInfo, GroupInfo
    USE TH_mod,         ONLY : THVar
    USE PE_MOD,         ONLY : PE       
    USE FILES,          ONLY : io8, io18, caseid, LOCALFN
    USE IOUTIL,         ONLY : message,           ShowHbar1,            ShowHbar2, openfile
    USE CNTL,           ONLY : nTracerCntl
    USE itrcntl_mod,    ONLY : ItrCntl
    USE EFT_MOD    
    INTEGER :: icase = 0, itest = 0
    LOGICAL :: Master, Slave
    LOGICAL, SAVE :: lFirst
    Data lFirst /.TRUE./
    LOGICAL :: lfind
    REAL :: T1, T2, e1, e2, T, modT, tmod, e, avgrho, constrho, mdot, inletT
    REAL :: dum1, dum2, dum3
    
    Master = PE%Master; Slave = PE%Slave
    IF( nTracerCntl%EFT_nboron .EQ. 0 )THEN
        nTracerCntl%EFT_nboron=1
        nTracerCntl%EFT_boron(1)=nTracerCntl%BoronPPM
    ENDIF
    
    IF(PE%master) THEN
        localfn=trim(caseid)//'.EFT'
        CALL openfile(io18,FALSE,FALSE,FALSE, localfn)
    ENDIF
        
    DO iboron = 1, nTracerCntl%EFT_nboron
        nTracerCntl%BoronPPM=nTracerCntl%EFT_boron(iboron)
        IF(MASTER) CALL ShowHbar2(io8)
        WRITE(mesg, '(A,f6.1,A)') 'EFT search at Boron : ', nTracerCntl%BoronPPM, ' ppm'
        IF(Master) CALL message(io8, TRUE, TRUE, mesg)
        
    DO icase = 1, nTracerCntl%nEFTpts
        IF(MASTER) CALL ShowHbar2(io8)
        WRITE(mesg, '(A,f6.1,A)') 'Performing Reference case for Eff. Fuel Temperature search at P = ', nTracerCntl%P_EFT(icase), ' %'
        IF(Master) CALL message(io8, TRUE, TRUE, mesg)
        !--- reference
        P_eft(icase)=nTracerCntl%P_EFT(icase)
        nTracerCntl%PowerLevel=P_eft(icase)* 0.01_8
        nTracerCntl%lFeedBack=.TRUE.
        lEFT_GCGEN=.TRUE.;  iboroncase=iboron;  ipowercase=icase;  calcid=1
        IF( .NOT. nTracerCntl%lMATRA )THEN
            ItrCntl%cmfdit=0
            CALL InitTH()
            CALL SteadyStateTH(Core, CmInfo, FmInfo, ThInfo, Eigv, ng, GroupInfo, nTracerCntl, PE)        
            CALL GET_T(centT(icase),meanT(icase),surfT(icase),modT,Core,THInfo,THVar)
            CALL GET_RHO(avgrho,Core,THInfo,THVar)
            rho_mod(icase)=avgrho
            T_mod(icase)=modT
        ENDIF
        CALL InitTH()
        CALL SSEIG()
        IF( nTracerCntl%lMATRA )THEN
            CALL GET_T(centT(icase),meanT(icase),surfT(icase),modT,Core,THInfo,THVar)
            CALL GET_RHO(avgrho,Core,THInfo,THVar)
            rho_mod(icase)=avgrho
            T_mod(icase)=modT
        ENDIF
        lEFT_GCGEN=.FALSE.; 
        nTracerCntl%lFeedBack=.FALSE.
        k_ref(icase)=eigv
        IF(Master) CALL ShowHbar1(io8)
        WRITE(mesg, '(A,4f7.2)') 'Reference case : Cent-Mean-Surf-Mod T.(C) = ',centT(icase),meanT(icase),surfT(icase),modT
        IF(Master) CALL message(io8, TRUE, TRUE, mesg)
        WRITE(mesg, '(A,1f12.5)') '               : Volume average mod. density (g/cc) = ',avgrho
        IF(Master) CALL message(io8, TRUE, TRUE, mesg)
        CALL STL(constrho, modT)
        WRITE(mesg, '(A,1f12.5,a,f7.1)') '               : Rho(avgT) mod. density (g/cc)      = ',constrho, ' , r(pcm) = ', (constrho/avgrho-1)*100000
        IF(Master) CALL message(io8, TRUE, TRUE, mesg)        
        T_emt(icase)=T_mod(icase)
        !--- CONST COOLANT
        IF( nTracerCntl%EFTUniCoolant .EQ. 1 )THEN !old approach
            IF( .NOT. nTracerCntl%lMATRA )THEN
                IF(Master) CALL ShowHbar1(io8)
                WRITE(mesg, '(A,f6.2,A)') 'Setting Infinite Coolant flowrate with T_in =', T_mod(icase), ' C'
                IF(Master) CALL message(io8, TRUE, TRUE, mesg)
            ELSE  
                CALL GET_T(dum1,dum2,dum3,modT,Core,THInfo,THVar)
                WRITE(mesg, '(A,f6.2,A)') 'Setting Infinite Coolant flowrate from MATRA with T_in =', T_mod(icase), ' C'
                IF(Master) CALL message(io8, TRUE, TRUE, mesg)
            ENDIF
            mdot=nTracerCntl%fMdotFA    
            inletT=nTracerCntl%TempInlet  
            nTracerCntl%fMdotFA=1000*mdot
            nTracerCntl%TempInlet=T_mod(icase)
            nTracerCntl%lFeedBack=.TRUE.
            lEFT_GCGEN=.TRUE.;  iboroncase=iboron;  ipowercase=icase;  calcid=2
            CALL InitTH()
            CALL SSEIG()
            lEFT_GCGEN=.FALSE.; 
            nTracerCntl%lFeedBack=.FALSE.    
            nTracerCntl%fMdotFA=mdot
            nTracerCntl%TempInlet=inletT
            modeff(icase)=-(eigv-k_ref(icase))*100000
            WRITE(mesg, '(A,f7.1,A,f5.1)') 'Moderator Temperature Effect at P = ', P_eft(icase), ' %, del_k = ',modeff(icase)
            IF(Master) CALL message(io8, TRUE, TRUE, mesg) 
            IF(Master) CALL message(io18, FALSE, TRUE, mesg)
            k_ref(icase)=eigv  
        ELSEIF( nTracerCntl%EFTUniCoolant .EQ. 2 )THEN ! 
            IF( .NOT. nTracerCntl%lMATRA )THEN
                IF(Master) CALL ShowHbar1(io8)
                WRITE(mesg, '(A,f7.2,A)') 'Setting uniform coolant temperature and density T_mod = ', T_emt(icase), ' C'
                IF(Master) CALL message(io8, TRUE, TRUE, mesg)
                CALL SET_ModTemp(Core%ncorefxr,FmInfo,PE,T_emt(icase))
            ELSE
                CALL GET_RHO_MATRA(avgrho,Core,THInfo,THVar)
                IF(Master) CALL ShowHbar1(io8)
                WRITE(mesg, '(A,f12.5,A)') 'Setting uniform coolant density for MATRA rho = ', avgrho, ' g/cc'
                IF(Master) CALL message(io8, TRUE, TRUE, mesg)
                CALL SET_Mod_MATRA(Core%ncorefxr,FmInfo,PE,T_emt(icase),avgrho)
            ENDIF
            lEFT_GCGEN=.TRUE.;  iboroncase=iboron;  ipowercase=icase;  calcid=2
            CALL SSEIG()
            lEFT_GCGEN=.FALSE.; 
            modeff(icase)=-(eigv-k_ref(icase))*100000
            WRITE(mesg, '(A,f7.1,A,f5.1)') 'Moderator Temperature Effect at P = ', P_eft(icase), ' %, del_k = ',modeff(icase)
            IF(Master) CALL message(io8, TRUE, TRUE, mesg) 
            k_ref(icase)=eigv
        ENDIF
        IF( iboron .EQ. 1 )THEN
        lfind = .FALSE.
        itest=0
        DO WHILE( .NOT. lfind )
            itest=itest+1
            IF(MASTER) CALL ShowHbar1(io8)
            WRITE(mesg, '(A,f7.1,A,i3)') 'Performing test case for power ', P_eft(icase), ' %, Test No. : ',itest
            IF(Master) CALL message(io8, TRUE, TRUE, mesg)
            IF( itest .EQ. 1 )THEN
                T_eft(icase,itest)=(1.0-0.78)*surfT(icase)+0.78*meanT(icase)
            ELSEIF( itest .EQ. 2 )THEN
                T_eft(icase,itest)=(1.0-0.64)*centT(icase)+0.64*surfT(icase)
            ELSE
                T1=sqrt(T_eft(icase,itest-2)+CKELVIN)
                T2=sqrt(T_eft(icase,itest-1)+CKELVIN)
                e1=e_eft(icase,itest-2)
                e2=e_eft(icase,itest-1)
                T=T2-e2*(T2-t1)/(e2-e1)
                T_eft(icase,itest)=T**2-CKELVIN
            ENDIF
            CALL GET_W_EFT(centT(icase),meanT(icase),surfT(icase),T_eft(icase,itest),w1(icase,itest),w2(icase,itest))
            WRITE(mesg, '(A,f7.2,A,f7.1,A)') 'Fuel T.(C) = ', T_eft(icase,itest), ' C for P = ', P_eft(icase), ' %'
            IF(Master) CALL message(io8, TRUE, TRUE, mesg)            
            CALL SET_FuelTemp(Core%ncorefxr,FmInfo,PE,T_eft(icase,itest))
            CALL SET_ModTemp(Core%ncorefxr,FmInfo,PE,T_emt(icase))
            calcid=3   
            IF( itest .GE. 3 ) lEFT_GCGEN=.TRUE.
            CALL SSEIG()
            lEFT_GCGEN=.FALSE.;
            k_eft(icase,itest)=eigv
            e_eft(icase,itest)=abs(eigv-k_ref(icase))*100000
            WRITE(mesg, '(A,f7.2,A,f7.1,A)') ' > Error at Fuel T.(C)', T_eft(icase,itest), ' C, E = ', (eigv-k_ref(icase))*100000, ' pcm'
            IF(Master) CALL message(io8, TRUE, TRUE, mesg)            
            IF( e_eft(icase,itest) .LT. nTracerCntl%crit_eft .OR. itest .EQ. nEFTP )THEN
                lfind=.TRUE.     
                ntest(icase)=itest
                !lEFT_GCGEN=.TRUE.;  iboroncase=iboron;  ipowercase=icase;  calcid=3                
                !IF(Master) CALL GCGen(Core, FmInfo, THInfo, CmInfo, GroupInfo, nTracerCntl, PE, ng)
                !lEFT_GCGEN=.FALSE.;
            ENDIF   
        ENDDO        
        ENDIF
    ENDDO
    IF(Master) CALL ShowHbar2(io8)
    IF(Master) CALL ShowHbar2(io18)
    WRITE(mesg, '(A,f6.1,A)') '              Summary of EFT search at Boron : ', nTracerCntl%BoronPPM, ' ppm'
    IF(Master) CALL message(io8, FALSE, TRUE, mesg)
    IF(Master) CALL message(io18, FALSE, TRUE, mesg)
    IF( iboron .EQ. 1 )THEN
    DO icase = 1, nTracerCntl%nEFTpts
        itest=ntest(icase)
        WRITE(mesg, '(A,f6.1,A,f7.2,f8.2,3f7.2,f8.5)')  ' P :', P_eft(icase), ' %, EFT-CMS-Mod T(C) =',T_eft(icase,itest),centT(icase),meanT(icase),surfT(icase),T_mod(icase),rho_mod(icase)
        IF(Master) CALL message(io8, FALSE, TRUE, mesg)
        IF(Master) CALL message(io18, FALSE, TRUE, mesg)
    ENDDO
    IF(Master) CALL ShowHbar2(io8)
    IF(Master) CALL ShowHbar2(io18)
    DO icase = 1, nTracerCntl%nEFTpts
        itest=ntest(icase)
        WRITE(mesg, '(A,f6.1,A,f7.5,A,f6.3,A,f6.3,A,f3.1,A,i2)') ' P :', P_eft(icase), ' %, k-eff = ', k_ref(icase),' w_s = ', w1(icase,itest), ' , w_n = ', w2(icase,itest), ' , E = ', e_eft(icase,itest),' Try : ',itest
        IF(Master) CALL message(io8, FALSE, TRUE, mesg)
        IF(Master) CALL message(io18, FALSE, TRUE, mesg)
    ENDDO
    IF(Master) CALL ShowHbar2(io8)
    IF(Master) CALL ShowHbar2(io18)
    ENDIF
    DO icase = 1, nTracerCntl%nEFTpts
        WRITE(mesg, '(A,f7.1,A,f6.1)') 'Moderator Temperature Effect at P = ', P_eft(icase), ' %, del_k = ',modeff(icase)
        IF(Master) CALL message(io18, FALSE, TRUE, mesg) 
    ENDDO
    IF(Master) CALL ShowHbar2(io8)
    IF(Master) CALL ShowHbar2(io18)
    ENDDO ! END OF BORON LOOP
    IF(Master) close(io18)
ENDSUBROUTINE
SUBROUTINE GET_W_EFT(centT, meanT, surfT, Teft, w1, w2)
    REAL :: centT, meanT, surfT, Teft, w1, w2
    w1 = (Teft-surfT)/(meanT-surfT)
    w2 = (Teft-centT)/(surfT-centT)
ENDSUBROUTINE
SUBROUTINE SET_FuelTemp(nFXR, FmInfo, PE, Tfuel)
    USE PARAM
    USE TYPEDEF, ONLY : FmInfo_Type, FxrInfo_TYPE, PE_TYPE
    USE CNTL, ONLY : nTracerCntl
    IMPLICIT NONE
    TYPE(FMInfo_Type) :: FmInfo
    TYPE(FxrInfo_Type), POINTER :: FXR(:,:)
    TYPE(PE_Type) :: PE
    REAL :: Tfuel
    INTEGER :: ifxr, nFXR, iz    
    DO iz= PE%myzb, PE%myze
        FXR=>Fminfo%FXR
        DO ifxr = 1, nFXR
            IF( FXR(ifxr,iz)%lfuel )THEN
                FXR(ifxr,iz)%temp=tfuel+CKELVIN
            ENDIF
        ENDDO
    ENDDO
ENDSUBROUTINE
SUBROUTINE STL(wrho, tmod)
    USE PARAM
    USE CNTL, ONLY : nTracerCntl
    USE SteamTBL_mod, ONLY : steamtbl
    REAL :: pexit, aw, ndenh2o
    REAL :: wt, wh, hin, wrho, wvin, wxin, wbetain, wkapain, wcpin
    PEXIT = nTracerCntl%PExit         
    wt=tmod+CKELVIN
    CALL steamtbl(TRUE,pexit,wt,wh,wrho,wvin,wxin,wbetain,wkapain,wcpin)
    wrho=wrho/1000
ENDSUBROUTINE   

SUBROUTINE SET_ModTemp(nFXR, FmInfo, PE, Tmod)
    USE PARAM
    USE NuclidMap_mod, ONLY : AtomicWeight
    USE SteamTBL_mod, ONLY : steamtbl
    USE TYPEDEF, ONLY : FmInfo_Type, FxrInfo_TYPE, PE_TYPE
    USE CNTL, ONLY : nTracerCntl
    USE FILES,            ONLY : io8
    USE IOUTIL,           ONLY : message
    IMPLICIT NONE
    TYPE(FMInfo_Type) :: FmInfo
    TYPE(FxrInfo_Type), POINTER :: FXR(:,:)
    TYPE(PE_Type) :: PE
    REAL :: Tmod
    INTEGER :: ifxr, nFXR, iz, idx, hidx, oidx
    REAL :: pexit, aw, ndenh2o
    REAL :: wt, wh, hin, wrho, wvin, wxin, wbetain, wkapain, wcpin
    PEXIT = nTracerCntl%PExit
    DO iz= PE%myzb, PE%myze
        FXR=>Fminfo%FXR
        DO ifxr = 1, nFXR
            IF( FXR(ifxr,iz)%lh2o )THEN
                FXR(ifxr,iz)%temp=tmod+CKELVIN                
                wt=FXR(ifxr,iz)%temp
                CALL steamtbl(TRUE,pexit,wt,wh,wrho,wvin,wxin,wbetain,wkapain,wcpin)
                DO idx = 1, FXR(ifxr,iz)%ndim
                    SELECT CASE(FXR(ifxr,iz)%idiso(idx))
                    CASE(1001)
                        hidx=idx
                    CASE(8016)
                        oidx=idx
                    ENDSELECT
                ENDDO
                wrho=wrho/1000
                aw=2*AtomicWeight(1001)+AtomicWeight(8016)      ! 170911 BYS JJH EDIT
                ndenh2o=wrho* AVOGADRO / aw                     ! *0.03344 !(0.03344=6.02/18)
                FXR(ifxr,iz)%pnum(oidx)=ndenh2o
                FXR(ifxr,iz)%pnum(hidx)=2*ndenh2o                
            ENDIF
        ENDDO
    ENDDO    
    IF(PE%MASTER) THEN
        WRITE(mesg,'(a,f7.5,a)') '  FXR_VARIATION : Applying Water gram density  = ', wrho, ' g/cc'
        CALL message(io8, TRUE, TRUE, mesg)
    ENDIF
ENDSUBROUTINE
SUBROUTINE SAVE_FuelT(Core, THInfo, THVar)
    USE PARAM
    USE TYPEDEF, ONLY : CoreInfo_Type,      ThInfo_Type ,      CoolantTh_Type,  &
                        FuelTH_Type,        THVar_Type,                         &
                        AsyInfo_Type,       Asy_Type
    USE EFT_MOD, ONLY : FuelTH0
    USE ALLOCS
    IMPLICIT NONE
    REAL :: centT, meanT, surfT, modT
    REAL :: centTz, meanTz, surfTz, modTz
    
    TYPE(CoreInfo_Type) :: Core
    TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
    TYPE(Asy_Type), POINTER :: Asy(:)
    
    TYPE(ThInfo_Type) :: ThInfo
    TYPE(FuelTH_Type), POINTER :: FuelTH(:)
    TYPE(THVar_Type) :: THVar
    
    REAL, POINTER :: hz(:)
    REAL :: hz_tot
    INTEGER :: iz, ixy, ixya, ixa, iya, iasy, AsyType, ir, i
    INTEGER :: nz, nxy, nxya, nxa, nya, nfuelpin, nrp
    INTEGER :: navg
    
    nz = Core%nz; nxy = Core%nxy
    nxya = Core%nxya;
    AsyInfo => Core%AsyInfo; Asy => Core%Asy
    hz=>Core%hz
    FuelTH => THInfo%FuelTH
    nrp = THVar%npr
    ALLOCATE(FuelTH0(nxy))
    DO i = 1, nxy
        IF(.NOT. Core%Pin(i)%lFuel)CYCLE
        CALL Dmalloc(FuelTh0(i)%tfuel, nrp + 5, nz)
    ENDDO
    DO iz = 1, nz
        DO iasy = 1, nxya
            AsyType = Asy(iasy)%AsyType
            IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
            nxya = AsyInfo(AsyType)%nxy
            DO ixya = 1, nxya
                ixy = Asy(iasy)%GlobalPinIdx(ixya)
                IF(.NOT. FuelTh(ixy)%lfuel) CYCLE
                DO ir = 1, nrp+5
                    FuelTh0(ixy)%Tfuel(ir, iz)=FuelTh(ixy)%Tfuel(ir, iz)
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDSUBROUTINE    
SUBROUTINE SET_FuelT(Core, THInfo, THVar)
    USE PARAM
    USE TYPEDEF, ONLY : CoreInfo_Type,      ThInfo_Type ,      CoolantTh_Type,  &
                        FuelTH_Type,        THVar_Type,                         &
                        AsyInfo_Type,       Asy_Type
    USE EFT_MOD, ONLY : FuelTH0
    IMPLICIT NONE
    REAL :: centT, meanT, surfT, modT
    REAL :: centTz, meanTz, surfTz, modTz
    
    TYPE(CoreInfo_Type) :: Core
    TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
    TYPE(Asy_Type), POINTER :: Asy(:)
    
    TYPE(ThInfo_Type) :: ThInfo
    TYPE(FuelTH_Type), POINTER :: FuelTH(:)
    TYPE(THVar_Type) :: THVar
    
    REAL, POINTER :: hz(:)
    REAL :: hz_tot
    INTEGER :: iz, ixy, ixya, ixa, iya, iasy, AsyType, ir
    INTEGER :: nz, nxy, nxya, nxa, nya, nfuelpin, nrp
    INTEGER :: navg
    
    nz = Core%nz; nxy = Core%nxy
    nxya = Core%nxya;
    AsyInfo => Core%AsyInfo; Asy => Core%Asy
    hz=>Core%hz
    FuelTH => THInfo%FuelTH
    nrp = THVar%npr
    DO iz = 1, nz
        DO iasy = 1, nxya
            AsyType = Asy(iasy)%AsyType
            IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
            nxya = AsyInfo(AsyType)%nxy
            DO ixya = 1, nxya
                ixy = Asy(iasy)%GlobalPinIdx(ixya)
                IF(.NOT. FuelTh(ixy)%lfuel) CYCLE
                DO ir = 1, nrp+5
                    FuelTh(ixy)%Tfuel(ir, iz)=FuelTh0(ixy)%Tfuel(ir, iz)
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDSUBROUTINE    
               
SUBROUTINE GET_RHO(avgrho, Core, THInfo, THVar)
    USE PARAM
    USE SteamTBL_mod, ONLY : steamtbl
    USE TYPEDEF, ONLY : CoreInfo_Type,      ThInfo_Type ,      CoolantTh_Type,  &
                        FuelTH_Type,        THVar_Type,                         &
                        AsyInfo_Type,       Asy_Type
    USE CNTL, ONLY : nTracerCntl
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
    TYPE(Asy_Type), POINTER :: Asy(:)
    
    TYPE(ThInfo_Type) :: ThInfo
    TYPE(FuelTH_Type), POINTER :: FuelTH(:)
    TYPE(CoolantTH_Type), POINTER :: CoolantTH(:)
    TYPE(THVar_Type) :: THVar
    
    REAL, POINTER :: hz(:)
    REAL :: hz_tot, avgrho, totarea, agt, af
    INTEGER :: iz, ixy, ixya, ixa, iya, iasy, AsyType
    INTEGER :: nz, nxy, nxya, nxa, nya, nfuelpin
    INTEGER :: navg
    REAL :: pexit, aw, ndenh2o
    REAL :: wt, wh, hin, wrho, wvin, wxin, wbetain, wkapain, wcpin
    
    PEXIT = nTracerCntl%PExit
    nz = Core%nz; nxy = Core%nxy
    nxya = Core%nxya;
    AsyInfo => Core%AsyInfo; Asy => Core%Asy
    hz=>Core%hz
    FuelTH => THInfo%FuelTH
    CoolantTH => THInfo%CoolantTH
    agt=THVar%acfgt
    af=THVar%acf
    nfuelpin=0;
    avgrho=0.0_8;
    DO iz = 1, nz
        DO iasy = 1, nxya
            AsyType = Asy(iasy)%AsyType
            IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
            nxya = AsyInfo(AsyType)%nxy
            DO ixya = 1, nxya
                ixy = Asy(iasy)%GlobalPinIdx(ixya)
                !IF(.NOT. CoolantTH(ixy)%lfuel)THEN
                !    avgrho=avgrho+agt*CoolantTH(ixy)%DenCool(iz)*hz(iz)
                !    totarea=totarea+agt*hz(iz)
                !ELSE
                !    avgrho=avgrho+ af*CoolantTH(ixy)%DenCool(iz)*hz(iz)
                !    totarea=totarea+af*hz(iz)
                !ENDIF
                IF(.NOT. FuelTH(ixy)%lfuel)CYCLE
                wt=FuelTh(ixy)%Tcool(iz)+CKELVIN
                CALL steamtbl(TRUE,pexit,wt,wh,wrho,wvin,wxin,wbetain,wkapain,wcpin)
                nfuelpin=nfuelpin+1
                avgrho=avgrho+ wrho                   
            ENDDO
        ENDDO
    ENDDO
    !avgrho=avgrho/totarea
    avgrho=avgrho/nfuelpin/1000
    NULLIFY(AsyInfo)
    NULLIFY(Asy)
END SUBROUTINE 
SUBROUTINE GET_RHO_MATRA(avgrho, Core, THInfo, THVar)
    USE PARAM
    USE SteamTBL_mod, ONLY : steamtbl
    USE TYPEDEF, ONLY : CoreInfo_Type,      ThInfo_Type ,      CoolantTh_Type,  &
                        FuelTH_Type,        THVar_Type,                         &
                        AsyInfo_Type,       Asy_Type
    USE CNTL, ONLY : nTracerCntl
    IMPLICIT NONE
    TYPE(CoreInfo_Type) :: Core
    TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
    TYPE(Asy_Type), POINTER :: Asy(:)
    
    TYPE(ThInfo_Type) :: ThInfo
    TYPE(FuelTH_Type), POINTER :: FuelTH(:)
    TYPE(CoolantTH_Type), POINTER :: CoolantTH(:)
    TYPE(THVar_Type) :: THVar
    
    REAL, POINTER :: hz(:)
    REAL :: hz_tot, avgrho, totarea, agt, af
    INTEGER :: iz, ixy, ixya, ixa, iya, iasy, AsyType
    INTEGER :: nz, nxy, nxya, nxa, nya, nfuelpin
    INTEGER :: navg
    REAL :: pexit, aw, ndenh2o
    REAL :: wt, wh, hin, wrho, wvin, wxin, wbetain, wkapain, wcpin
    
    PEXIT = nTracerCntl%PExit
    nz = Core%nz; nxy = Core%nxy
    nxya = Core%nxya;
    AsyInfo => Core%AsyInfo; Asy => Core%Asy
    hz=>Core%hz
    FuelTH => THInfo%FuelTH
    CoolantTH => THInfo%CoolantTH
    agt=THVar%acfgt
    af=THVar%acf
    nfuelpin=0;
    avgrho=0.0_8;
    DO iz = 1, nz
        DO iasy = 1, nxya
            AsyType = Asy(iasy)%AsyType
            IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
            nxya = AsyInfo(AsyType)%nxy
            DO ixya = 1, nxya
                ixy = Asy(iasy)%GlobalPinIdx(ixya)
                IF(.NOT. FuelTH(ixy)%lfuel)CYCLE
                wrho=THInfo%Dencool(iz,ixy)
                nfuelpin=nfuelpin+1
                avgrho=avgrho+ wrho
            ENDDO
        ENDDO
    ENDDO
    avgrho=avgrho/nfuelpin/1000
    NULLIFY(AsyInfo)
    NULLIFY(Asy)
END SUBROUTINE 

SUBROUTINE SET_Mod_MATRA(nFXR, FmInfo, PE, Tmod, wrho)
    USE PARAM
    USE NuclidMap_mod, ONLY : AtomicWeight
    USE SteamTBL_mod, ONLY : steamtbl
    USE TYPEDEF, ONLY : FmInfo_Type, FxrInfo_TYPE, PE_TYPE
    USE CNTL, ONLY : nTracerCntl
    USE FILES,            ONLY : io8
    USE IOUTIL,           ONLY : message
    IMPLICIT NONE
    TYPE(FMInfo_Type) :: FmInfo
    TYPE(FxrInfo_Type), POINTER :: FXR(:,:)
    TYPE(PE_Type) :: PE
    REAL :: Tmod
    INTEGER :: ifxr, nFXR, iz, idx, hidx, oidx
    REAL :: pexit, aw, ndenh2o
    REAL :: wt, wh, hin, wrho, wvin, wxin, wbetain, wkapain, wcpin
    PEXIT = nTracerCntl%PExit
    DO iz= PE%myzb, PE%myze
        FXR=>Fminfo%FXR
        DO ifxr = 1, nFXR
            IF( FXR(ifxr,iz)%lh2o )THEN
                FXR(ifxr,iz)%temp=tmod+CKELVIN                
                DO idx = 1, FXR(ifxr,iz)%ndim
                    SELECT CASE(FXR(ifxr,iz)%idiso(idx))
                    CASE(1001)
                        hidx=idx
                    CASE(8016)
                        oidx=idx
                    ENDSELECT
                ENDDO
                aw=2*AtomicWeight(1001)+AtomicWeight(8016)      ! 170911 BYS JJH EDIT
                ndenh2o=wrho* AVOGADRO / aw                     ! *0.03344 !(0.03344=6.02/18)
                FXR(ifxr,iz)%pnum(oidx)=ndenh2o
                FXR(ifxr,iz)%pnum(hidx)=2*ndenh2o                
            ENDIF
        ENDDO
    ENDDO    
    IF(PE%MASTER) THEN
        WRITE(mesg,'(a,f7.5,a)') '  FXR_VARIATION : Applying Water gram density  = ', wrho, ' g/cc'
        CALL message(io8, TRUE, TRUE, mesg)
    ENDIF
ENDSUBROUTINE

SUBROUTINE GET_T(centT, meanT, surfT, modT, Core, THInfo, THVar)
    USE PARAM
    USE TYPEDEF, ONLY : CoreInfo_Type,      ThInfo_Type ,      CoolantTh_Type,  &
                        FuelTH_Type,        THVar_Type,                         &
                        AsyInfo_Type,       Asy_Type
    IMPLICIT NONE
    REAL :: centT, meanT, surfT, modT
    REAL :: centTz, meanTz, surfTz, modTz
    
    TYPE(CoreInfo_Type) :: Core
    TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
    TYPE(Asy_Type), POINTER :: Asy(:)
    
    TYPE(ThInfo_Type) :: ThInfo
    TYPE(FuelTH_Type), POINTER :: FuelTH(:)
    TYPE(THVar_Type) :: THVar
    
    REAL, POINTER :: hz(:)
    REAL :: hz_tot
    INTEGER :: iz, ixy, ixya, ixa, iya, iasy, AsyType
    INTEGER :: nz, nxy, nxya, nxa, nya, nfuelpin
    INTEGER :: navg
    
    nz = Core%nz; nxy = Core%nxy
    nxya = Core%nxya;
    AsyInfo => Core%AsyInfo; Asy => Core%Asy
    hz=>Core%hz
    FuelTH => THInfo%FuelTH
    navg = THVar%npr5
    centT=0._8;
    meanT=0._8;
    surfT=0._8;
    modT=0._8;
    hz_tot=0._8;
    DO iz = 1, nz
        centTz=0._8;
        meanTz=0._8;
        surfTz=0._8;
        modTz=0._8;
        nfuelpin=0
        DO iasy = 1, nxya
            AsyType = Asy(iasy)%AsyType
            IF(.NOT. AsyInfo(AsyType)%lFuel) CYCLE
            nxya = AsyInfo(AsyType)%nxy
            DO ixya = 1, nxya
                ixy = Asy(iasy)%GlobalPinIdx(ixya)
                IF(.NOT. FuelTh(ixy)%lfuel) CYCLE
                nfuelpin=nfuelpin+1
                centTz =centTz +FuelTh(ixy)%Tfuel(1, iz)  !cent
                meanTz =meanTz +FuelTh(ixy)%tfuel(navg, iz) ! avg
                surfTz =surfTz +FuelTh(ixy)%Tfuel(ThVar%npr1, iz) ! surf
                modTz =modTz +FuelTh(ixy)%Tcool(iz) ! surf
                !WRITE(123123,'(2f17.5)')  FuelTh(ixy)%Tfuel(1, iz), FuelTh(ixy)%Tcool(iz) ! surf
            ENDDO
        ENDDO
        centTz=centTz/nfuelpin
        meanTz=meanTz/nfuelpin
        surfTz=surfTz/nfuelpin
        modTz=modTz/nfuelpin
        centT=centT+centTz*hz(iz)
        meanT=meanT+meanTz*hz(iz)
        surfT=surfT+surfTz*hz(iz)
        modT=modT+modTz*hz(iz)
        hz_tot=hz_tot+hz(iz)
    ENDDO
    centT=centT/hz_tot
    meanT=meanT/hz_tot
    surfT=surfT/hz_tot
    modT=modT/hz_tot
    
    NULLIFY(AsyInfo)
    NULLIFY(Asy)
    NULLIFY(FuelTH)
END SUBROUTINE

SUBROUTINE EFT_driver_old
    USE PARAM
    USE Geom,           ONLY : Core
    USE Core_Mod,       ONLY : FmInfo, eigv, THinfo
    USE TH_mod,         ONLY : THVar
    USE PE_MOD,         ONLY : PE       
    USE FILES,          ONLY : io8
    USE IOUTIL,         ONLY : message,           ShowHbar1,            ShowHbar2
    USE CNTL,           ONLY : nTracerCntl
    USE EFT_MOD
    INTEGER :: icase = 0, itest = 0
    LOGICAL :: Master, Slave
    LOGICAL, SAVE :: lFirst
    Data lFirst /.TRUE./
    LOGICAL :: lfind
    REAL :: T1, T2, e1, e2, T, modT
    
    Master = PE%Master; Slave = PE%Slave
    
    DO icase = 1, nTracerCntl%nEFTpts
        IF(MASTER) CALL ShowHbar2(io8)
        WRITE(mesg, '(A)') 'Performing Reference case for Eff. Fuel Temperature search ...'
        IF(Master) CALL message(io8, TRUE, TRUE, mesg)
        !--- reference
        P_eft(icase)=nTracerCntl%P_EFT(icase)
        nTracerCntl%PowerLevel=P_eft(icase)* 0.01_8
        nTracerCntl%lFeedBack=.TRUE.
        CALL InitTH()
        CALL SSEIG()
        nTracerCntl%lFeedBack=.FALSE.
        k_ref(icase)=eigv
        CALL GET_T(centT(icase),meanT(icase),surfT(icase),modT,Core,THInfo,THVar)
        IF(Master) CALL ShowHbar1(io8)
        WRITE(mesg, '(A,4f7.2)') 'Reference case : Cent-Mean-Surf-Mod T.(C) = ',centT(icase),meanT(icase),surfT(icase),modT
        IF(Master) CALL message(io8, TRUE, TRUE, mesg)
        lfind = .FALSE.
        itest=0
        DO WHILE( .NOT. lfind )
            itest=itest+1
            IF(MASTER) CALL ShowHbar1(io8)
            WRITE(mesg, '(A,f7.1,A,i3)') 'Performing test case for power ', P_eft(icase), ' %, Test No. : ',itest
            IF(Master) CALL message(io8, TRUE, TRUE, mesg)
            IF( itest .EQ. 1 )THEN
                T_eft(icase,itest)=(1.0-0.84)*surfT(icase)+0.84*meanT(icase)
            ELSEIF( itest .EQ. 2 )THEN
                T_eft(icase,itest)=(1.0-0.62)*centT(icase)+0.62*surfT(icase)
            ELSE
                T1=sqrt(T_eft(icase,itest-2)+CKELVIN)
                T2=sqrt(T_eft(icase,itest-1)+CKELVIN)
                e1=e_eft(icase,itest-2)
                e2=e_eft(icase,itest-1)
                T=T2-e2*(T2-t1)/(e2-e1)
                T_eft(icase,itest)=T**2-CKELVIN
            ENDIF
            CALL GET_W_EFT(centT(icase),meanT(icase),surfT(icase),T_eft(icase,itest),w1(icase,itest),w2(icase,itest))
            WRITE(mesg, '(A,f7.2,A,f5.3,A,f5.3)') 'Fuel T.(C) = ', T_eft(icase,itest), ' C, w_Std = ', w1(icase,itest), ' , w_NEA = ', w2(icase,itest)
            IF(Master) CALL message(io8, TRUE, TRUE, mesg)            
            CALL SET_FuelTemp(Core%ncorefxr,FmInfo,PE,T_eft(icase,itest))
            CALL SSEIG()
            k_eft(icase,itest)=eigv
            e_eft(icase,itest)=abs(eigv-k_ref(icase))*100000
            WRITE(mesg, '(A,f7.2,A,f5.3,A,f5.3,A,f4.1)') ' > Fuel T.(C)', T_eft(icase,itest), ' C, w_Std = ', w1(icase,itest), ' , w_NEA = ', w2(icase,itest), ' , E = ', e_eft(icase,itest)
            IF(Master) CALL message(io8, TRUE, TRUE, mesg)            
            IF( e_eft(icase,itest) .LT. nTracerCntl%crit_eft .OR. itest .EQ. nEFTP )THEN
                lfind=.TRUE.     
                ntest(icase)=itest
            ENDIF   
        ENDDO        
    ENDDO
    IF(Master) CALL ShowHbar2(io8)
    DO icase = 1, nTracerCntl%nEFTpts
        itest=ntest(icase)
        WRITE(mesg, '(A,f6.1,A,f7.5,A,f5.3,A,f5.3,A,f3.1)') 'P :', P_eft(icase), ' %, k-eff = ', k_ref(icase),' w_s = ', w1(icase,itest), ' , w_n = ', w2(icase,itest), ' , E = ', e_eft(icase,itest)
        IF(Master) CALL message(io8, TRUE, TRUE, mesg)
    ENDDO
ENDSUBROUTINE