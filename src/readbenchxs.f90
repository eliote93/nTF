SUBROUTINE ReadBenchXs(indev,ioutdev, ng,  nPrec, nXsChange, ScatOd)
!Reading Benchmakr XS
USE param
USE cntl,           only : nTracerCntl
!USE PE_Mod,         ONLY : Master
USE files,          only : filename,    InputFileIdx,  io5,            io8
USE inputcards ,    only : oneline,     probe,         mxcard,         nblock !,      &
                         !  FindBlockId, FindCardId,    cards,          blocks
use ioutil,         only : terminate,   toupper,       openfile,       IFnumeric,   &
                           nfields,     message
use BenchXs,        only : nxsltype,    ngben,         MacXsBen,       lKinParam,   &
                           nxschangeben,  benchxstype,                                &
                           DNP_BETA,    DNP_Lambda,    Neut_Velo,      Dneut_Chi
IMPLICIT NONE

INTEGER, intent(in) :: indev,ioutdev           !Input and Output Device Unit number
INTEGER :: nPreC                               !Number of Delayed Neutron Precursor Group
INTEGER :: ng                                  !Number of Energy Group
INTEGER :: nXsChange                           !Number of XS change for transient
!LOGICAL :: lscat1                              !P1 Scattering Option Flag
INTEGER :: ScatOd



INTEGER :: nLineField
INTEGER :: ixsl, iso                           !Benchmark XS index
INTEGER :: i
INTEGER :: ig, igs
real :: sumxss, sumchi
CHARACTER(15)   :: cardname,  astring
CHARACTER(1024) :: XsOneLine

INTEGER :: infinp
LOGICAL :: linflow = .FALSE.
LOGICAL :: lmod = .FALSE.
INTEGER :: nmod, fiso, iiso
INTEGER, ALLOCATABLE :: modiso(:)

IF(linflow) THEN 
    infinp = 802
    OPEN(infinp, file='inflowinput.inp',status='OLD')
    read(infinp,*)
    read(infinp,*) nmod
    ALLOCATE(modiso(nmod))
    read(infinp,*)
    read(infinp,*) modiso
    read(infinp,*) 
    read(infinp,*) fiso
END IF 

benchxstype=1
nGben = ng; nXsChangeBen = nXsChange
nxsltype = 0
!SCAN Benchmark XS files
DO while(TRUE)
  read(indev, '(a256)', END = 199) oneline
  IF(probe.eq.BANG) cycle;     IF(probe.eq.POUND) cycle
  IF(oneline.eq.BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe .ne. BLANK) CALL terminate('First Column should be blank:'//trim(oneline))
  read(oneline,*) cardname
  CALL toupper(cardname)
  nLineField = nfields(oneline) - 1
  IF(cardname .eq. 'BASE_MICRO') then
    read(oneline,*) astring, ixsl;
    nXslType = max(nXslType,ixsl)
  ENDIF

  IF(cardname.eq.'DNP_NGRP') THEN
    read(oneline,*) astring,nprec   !Precurror Group
    lKinParam = .TRUE.
  ENDIF
ENDDO
199 continue
rewind(indev)
!Allocate Benchmark XS
CALL alloc_BenchXs(ngben, nXslType + nXsChange, scatod)
IF(lKinParam) CALL Alloc_BenchKinParam(nPrec, ng)
!Read Benchnark XS
DO while(TRUE)
  read(indev, '(a256)', END = 299) oneline
  IF(probe.eq.BANG) cycle;     IF(probe.eq.POUND) cycle          !Cycle Condition(go to next line)
  IF(oneline.eq.BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe.eq.DOT) exit;       IF(probe.eq.SLASH) exit           !Exit Condition
  read(oneline,*) cardname
  CALL toupper(cardname)
  nLineField = nfields(oneline) - 1
  SELECT CASE(CardName)
    CASE('BASE_MICRO')
      READ(oneline,*) astring, iso
      MacXsBen(iso)%lempty = FALSE
      !Reading Group-Wise XS - Total, absorption, nu-fs, kaffa-fs, chi
      DO ig = 1, ngben
        READ(indev, '(a256)') oneline
        IF(.not. ifnumeric(oneline)) THEN
          WRITE(mesg,'(2(a,i3),a)') 'Character Found in BASE_MICRO ',iso
          CALL Terminate(mesg)
        ENDIF
        READ(oneline, *) MacXsBen(iso)%xst(ig), MacXsBen(iso)%xsa(ig), MacXsBen(iso)%xsnf(ig), MacXsBen(iso)%xskf(ig), MacXsBen(iso)%chi(ig)
        IF(MacXsBen(iso)%xsnf(ig) .gt. epsm4) MacXsBen(iso)%lFuel = TRUE
      ENDDO
      ! Scattering XS
      DO ig = 1, ngben
        READ(indev, '(a1024)') XsOneLine
        IF(.not. ifnumeric(oneline)) THEN
          WRITE(mesg,'(2(a,i3),a)') 'Character Found in BASE_MICRO ',iso
          CALL Terminate(mesg)
        ENDIF
        read(XsOneLine,*) (MacXsBen(iso)%xss0(ig,igs), igs = 1, ngben)
      ENDDO
      !IF(lscat1) THEN
  IF(scatod .ge. 1)THEN
        ! HigherOrder Scattering XS
        DO ig = 1, ngben
            READ(indev, '(a1024)') XsOneLine
            IF(.not. ifnumeric(oneline)) THEN
            WRITE(mesg,'(2(a,i3),a)') 'Character Found in BASE_MICRO ',iso
            CALL Terminate(mesg)
            ENDIF
            read(XsOneLine,*) (MacXsBen(iso)%xss1(ig,igs), igs = 1, ngben)
        ENDDO
      IF(scatod .ge. 2)THEN
        DO ig = 1, ngben
            READ(indev, '(a1024)') XsOneLine
            IF(.not. ifnumeric(oneline)) THEN
            WRITE(mesg,'(2(a,i3),a)') 'Character Found in BASE_MICRO ',iso
            CALL Terminate(mesg)
            ENDIF
            read(XsOneLine,*) (MacXsBen(iso)%xss2(ig,igs), igs = 1, ngben)
        ENDDO
        IF(scatod .ge. 3)THEN
            DO ig = 1, ngben
                READ(indev, '(a1024)') XsOneLine
                IF(.not. ifnumeric(oneline)) THEN
                WRITE(mesg,'(2(a,i3),a)') 'Character Found in BASE_MICRO ',iso
                CALL Terminate(mesg)
                ENDIF
                read(XsOneLine,*) (MacXsBen(iso)%xss3(ig,igs), igs = 1, ngben)
            ENDDO
        ENDIF
      ENDIF
  ENDIF


      !ENDIF
      !CHI verifciation
      sumchi = sum(MacXsBen(iso)%chi(1:ngben))
      IF(sumchi .GT. epsm3) THEN
        sumchi = 1._8/sumchi
        MacXsBen(iso)%chi(1:ngben)  = sumchi*MacXsBen(iso)%chi(1:ngben)
      endif
      !Transport XS
      IF( ScatOd .EQ. 0 )THEN
        DO ig = 1 ,ng
        !sumXss = sum(MacXsBen(iso)%Xss(ig,1:ngben))
        !MacXsBen(iso)%Xss(ig,ig) = MacXsBen(iso)%xss(ig,ig) - sumXss- MacXsBen(iso)%xsa(ig)+MacXsBen(iso)%xst(ig)
        !MacXsBen(iso)%Xstr(ig) = MacXsBen(iso)%xst(ig)
        !MacXsBen(iso)%xst(ig) = sumXss+ MacXsBen(iso)%xsa(ig)

        !MacXsBen(iso)%Xstr(ig) = MacXsBen(iso)%xst(ig)  ! 170202 off



            sumXss = sum(MacXsBen(iso)%Xss(ig,1:ngben))
            MacXsBen(iso)%Xst(ig) = MacXsBen(iso)%xsa(ig) + sumXss
            MacXsBen(iso)%Xs0sum(ig) = sumXss
            MacXsBen(iso)%xstr(ig)=MacXsBen(iso)%xst(ig)
        ENDDO
      ELSE
        
        DO ig = 1 ,ng
            !macXsben(iso)%xstr(ig)=macxsben(iso)%xst(ig)
            !
            sumXss = sum(MacXsBen(iso)%Xss(ig,1:ngben))
            MacXsBen(iso)%Xst(ig) = MacXsBen(iso)%xsa(ig) + sumXss
            MacXsBen(iso)%Xs0sum(ig) = sumXss
            sumXss = sum(MacXsBen(iso)%Xss1(ig,1:ngben))
            MacXsBen(iso)%Xs1sum(ig) = sumXss
            
        END DO      
        DO iiso = 1, nmod
            IF(iso.eq.modiso(iiso)) lmod = .TRUE.
        END DO 
        IF(linflow .AND. lmod) THEN 
            CALL inflow_BenchXS(iso, fiso)
            lmod = .FALSE.
        ELSE
            DO ig = 1, ng
                MacXsBen(iso)%Xstr(ig) = MacXsBen(iso)%xst(ig)  - MacXsBen(iso)%Xs1sum(ig)
            END DO 
            !-- Transport corrected
            !MacXsBen(iso)%Xss(ig,ig)=MacXsBen(iso)%Xss(ig,ig)-sumXss
        END IF 
        
      ENDIF

    CASE('DNP_LAMBDA')
      READ(oneline, *) astring, DNP_Lambda(1:nprec)
    CASE('DNP_BETA')
      READ(oneline, *) astring, DNP_Beta(1:nprec)
    CASE('NEUT_VELO')
      READ(oneline, *) astring, Neut_Velo(1:ng)
    CASE('DNEUT_CHI')
      READ(oneline, *) astring, Dneut_chi(1:ng)
  END SELECT
ENDDO

299 continue

END SUBROUTINE

SUBROUTINE ReadBenchXs_CEA(indev,ioutdev, ng,  nPrec, nXsChange, ScatOd)
!Reading Benchmakr XS
USE param
USE cntl,           only : nTracerCntl
!USE PE_Mod,         ONLY : Master
USE files,          only : filename,    InputFileIdx,  io5,            io8
USE inputcards ,    only : oneline,     probe,         mxcard,         nblock !,      &
                         !  FindBlockId, FindCardId,    cards,          blocks
use ioutil,         only : terminate,   toupper,       openfile,       IFnumeric,   &
                           nfields,     message
use BenchXs,        only : nxsltype,    ngben,         MacXsBen,       lKinParam,   &
                           nxschangeben, benchxstype,                                            &
                           DNP_BETA,    DNP_Lambda,    Neut_Velo,      Dneut_Chi
IMPLICIT NONE

INTEGER, intent(in) :: indev,ioutdev           !Input and Output Device Unit number
INTEGER :: nPreC                               !Number of Delayed Neutron Precursor Group
INTEGER :: ng                                  !Number of Energy Group
INTEGER :: nXsChange                           !Number of XS change for transient
!LOGICAL :: lscat1                              !P1 Scattering Option Flag
INTEGER :: ScatOd



INTEGER :: nLineField
INTEGER :: ixsl, iso                           !Benchmark XS index
INTEGER :: i, iscat, idx
INTEGER :: ig, igs, nscat
INTEGER :: prof(12)
REAL :: scat1d(5), tempsm(ng,ng)
real :: sumxss, sumchi, scat
CHARACTER(10)   :: cardname,  astring
CHARACTER(1024) :: XsOneLine
INTEGER :: infinp
LOGICAL :: linflow = .TRUE.
LOGICAL :: lmod = .FALSE.
INTEGER :: nmod, fiso, iiso
INTEGER, ALLOCATABLE :: modiso(:)

IF(linflow) THEN 
    infinp = 802
    OPEN(infinp, file='inflowinput.inp',status='OLD')
    read(infinp,*)
    read(infinp,*) nmod
    ALLOCATE(modiso(nmod))
    read(infinp,*)
    read(infinp,*) modiso
    read(infinp,*) 
    read(infinp,*) fiso
END IF 


benchxstype=4
nGben = ng; nXsChangeBen = nXsChange
nxsltype = 0
!SCAN Benchmark XS files
DO while(TRUE)
  read(indev, '(a256)', END = 199) oneline
  IF(probe.eq.BANG) cycle;     IF(probe.eq.POUND) cycle
  IF(oneline.eq.BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe .ne. BLANK) CALL terminate('First Column should be blank:'//trim(oneline))
  read(oneline,*) cardname
  CALL toupper(cardname)
  nLineField = nfields(oneline) - 1
  IF(cardname .eq. 'BASE_MICRO') then
    read(oneline,*) astring, ixsl;
    nXslType = max(nXslType,ixsl)
  ENDIF

  IF(cardname.eq.'DNP_NGRP') THEN
    read(oneline,*) astring,nprec   !Precurror Group
    lKinParam = .TRUE.
  ENDIF
ENDDO
199 continue
rewind(indev)
!Allocate Benchmark XS
CALL alloc_BenchXs(ngben, nXslType + nXsChange, scatod)
IF(lKinParam) CALL Alloc_BenchKinParam(nPrec, ng)
!Read Benchnark XS
DO while(TRUE)
  read(indev, '(a256)', END = 299) oneline
  IF(probe.eq.BANG) cycle;     IF(probe.eq.POUND) cycle          !Cycle Condition(go to next line)
  IF(oneline.eq.BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe.eq.DOT) exit;       IF(probe.eq.SLASH) exit           !Exit Condition
  read(oneline,*) cardname
  CALL toupper(cardname)
  nLineField = nfields(oneline) - 1
  SELECT CASE(CardName)
    CASE('BASE_MICRO')
      READ(oneline,*) astring, iso
      MacXsBen(iso)%lempty = FALSE
      !Reading Group-Wise XS - Total, absorption, nu-fs, kaffa-fs, chi
      DO ig = 1, ngben
        READ(indev, '(a256)') oneline
        IF(.not. ifnumeric(oneline)) THEN
          WRITE(mesg,'(2(a,i3),a)') 'Character Found in BASE_MICRO ',iso
          CALL Terminate(mesg)
        ENDIF
        READ(oneline, *) MacXsBen(iso)%xst(ig), MacXsBen(iso)%xsa(ig), MacXsBen(iso)%xsnf(ig), MacXsBen(iso)%xskf(ig), MacXsBen(iso)%chi(ig)
        !WRITE(*,'(a)') oneline, 'ori'
        !WRITE(*, '(5e17.5)') MacXsBen(iso)%xst(ig), MacXsBen(iso)%xsa(ig), MacXsBen(iso)%xsnf(ig), MacXsBen(iso)%xskf(ig), MacXsBen(iso)%chi(ig)
        IF(MacXsBen(iso)%xsnf(ig) .gt. epsm4) MacXsBen(iso)%lFuel = TRUE
      ENDDO
      MacXsBen(iso)%lCR = FALSE
      IF(MacXsBen(iso)%xsa(ng) .gt. 1._8 .AND. .NOT.MacXsBen(iso)%lFuel) MacXsBen(iso)%lCR = TRUE
      ! Scattering XS
      idx=0
      READ(indev, '(a1024)') XsOneLine
      read(XsOneline,*) prof(1:nfields(XsOneLine))
      DO ig = 1, ngben
        IF( idx.EQ.12 )THEN
          READ(indev, '(a1024)') XsOneLine
          read(XsOneline,*) prof(1:nfields(XsOneLine))
          idx=0
        ENDIF
        idx=idx+1;
        DO iscat = 1, scatod+1
          MacXSBen(iso)%PnSM(iscat,ig)%ib=prof(idx)
        ENDDO
        idx=idx+1;
        DO iscat = 1, scatod+1
          MacXSBen(iso)%PnSM(iscat,ig)%ie=prof(idx)
        ENDDO
        DO iscat = 1, scatod+1
          ALLOCATE(MacXSBen(iso)%PnSM(iscat,ig)%from(MacXSBen(iso)%PnSM(iscat,ig)%ib:MacXSBen(iso)%PnSM(iscat,ig)%ie))
          MacXSBen(iso)%PnSM(iscat,ig)%from(MacXSBen(iso)%PnSM(iscat,ig)%ib:MacXSBen(iso)%PnSM(iscat,ig)%ie)=0._8
        ENDDO
      ENDDO
      DO iscat = 1, scatod+1
        READ(indev, '(a1024)') XsOneLine
        READ(XsOneLine,*) scat1d(1:nfields(XsOneLine))
        idx=0
        DO ig = 1, ngben
          MacXsBen(iso)%PnSm(iscat,ig)%withingroupscat=0._8
          DO igs = MacXsBen(iso)%PnSm(iscat,ig)%ib, MacXsBen(iso)%PnSm(iscat,ig)%ie
            IF( idx .EQ. 5)THEN
              READ(indev, '(a1024)') XsOneLine
              READ(XsOneLine,*) scat1d(1:nfields(XsOneLine))
              idx=0
            ENDIF
            idx=idx+1
            MacXsBen(iso)%PnSm(iscat,ig)%from(igs)=scat1d(idx)
          ENDDO
          IF( (ig.GE.MacXsBen(iso)%PnSm(iscat,ig)%ib).AND.(ig.LE.MacXsBen(iso)%PnSm(iscat,ig)%ie) )THEN
            MacXsBen(iso)%PnSm(iscat,ig)%withingroupscat=scat1d(idx)
          ENDIF
        ENDDO
      ENDDO

      !ENDIF
      !CHI verifciation
      sumchi = sum(MacXsBen(iso)%chi(1:ngben))
      IF(sumchi .GT. epsm3) THEN
        sumchi = 1._8/sumchi
        MacXsBen(iso)%chi(1:ngben)  = sumchi*MacXsBen(iso)%chi(1:ngben)
      endif

      !Transport XS  -- recent
      IF( ScatOd .EQ. 0 )THEN
          DO ig = 1, ng
              DO igs = 1, ng
                  tempsm(igs,ig)=0._8
              ENDDO
              DO igs = MacXsBen(iso)%PnSm(1,ig)%ib, MacXsBen(iso)%PnSm(1,ig)%ie
                  tempsm(igs,ig)=MacXsBen(iso)%PnSm(1,ig)%from(igs)
              ENDDO
          ENDDO
        DO ig = 1 ,ng
            sumXss=0._8
            DO igs = 1, ng
                sumXss=sumXss+tempsm(ig,igs)
            ENDDO
            MacXsBen(iso)%xst(ig)=sumXss+ MacXsBen(iso)%xsa(ig)
            MacXsBen(iso)%xstr(ig)=MacXsBen(iso)%xst(ig)
            MacXsBen(iso)%Xs0sum(ig) = sumXss
        ENDDO
      ELSE ! not P0
          DO ig = 1, ng
              DO igs = 1, ng
                  tempsm(igs,ig)=0._8
              ENDDO
              DO igs = MacXsBen(iso)%PnSm(1,ig)%ib, MacXsBen(iso)%PnSm(1,ig)%ie
                  tempsm(igs,ig)=MacXsBen(iso)%PnSm(1,ig)%from(igs)
              ENDDO
          ENDDO

        DO ig = 1 ,ng
            sumXss=0._8
            DO igs = 1, ng
                sumXss=sumXss+tempsm(ig,igs)
            ENDDO
            MacXsBen(iso)%Xst(ig) = MacXsBen(iso)%xsa(ig) + sumXss
            MacXsBen(iso)%Xs0sum(ig) = sumXss
        ENDDO

          DO ig = 1, ng
              DO igs = 1, ng
                  tempsm(igs,ig)=0._8
              ENDDO
              DO igs = MacXsBen(iso)%PnSm(2,ig)%ib, MacXsBen(iso)%PnSm(2,ig)%ie
                  tempsm(igs,ig)=MacXsBen(iso)%PnSm(2,ig)%from(igs)
              ENDDO
          ENDDO

        DO ig = 1, ng
            sumXss=0
            DO igs = 1, ng
                sumXss=sumXss+tempsm(ig,igs)
            ENDDO
            MacXsBen(iso)%Xstr(ig) = MacXsBen(iso)%xst(ig) - sumXss  ! 17/1/18 TR
            MacXsBen(iso)%Xs1sum(ig) = sumXss
            IF( MacXsBen(iso)%Xstr(ig) .LE. 0._8 ) WRITE(*,'(a,i5,i3,7e13.5)') 'Negative Sig_TR at', iso, ig, MacXsBen(iso)%xst(ig), MacXsBen(iso)%xsa(ig), MacXsBen(iso)%Xs0sum(ig), sumXss, MacXsBen(iso)%xstr(ig)
            !-- Transport corrected
            !IF( (ig.GE.MacXsBen(iso)%PnSm(1,ig)%ib).AND.(ig.LE.MacXsBen(iso)%PnSm(1,ig)%ie) )THEN
            !  MacXsBen(iso)%PnSm(1,ig)%withingroupscat=MacXsBen(iso)%PnSm(1,ig)%withingroupscat-MacXsBen(iso)%Xs1sum(ig)
            !  MacXsBen(iso)%PnSm(1,ig)%from(ig)=MacXsBen(iso)%PnSm(1,ig)%from(ig)-MacXsBen(iso)%Xs1sum(ig)
            !ENDIF
        ENDDO
        
        DO iiso = 1, nmod
            IF(iso.eq.modiso(iiso)) lmod = .TRUE.
        END DO 
        IF(linflow .AND. lmod) THEN 
            !CALL inflow_BenchXS_CEA(iso, fiso)            
            
            lmod = .FALSE.
        ELSE
            DO ig = 1, ng
                MacXsBen(iso)%Xstr(ig) = MacXsBen(iso)%xst(ig)  - MacXsBen(iso)%Xs1sum(ig)
            END DO 
            !-- Transport corrected
            !MacXsBen(iso)%Xss(ig,ig)=MacXsBen(iso)%Xss(ig,ig)-sumXss
        END IF 

      ENDIF


    CASE('DNP_LAMBDA')
      READ(oneline, *) astring, DNP_Lambda(1:nprec)
    CASE('DNP_BETA')
      READ(oneline, *) astring, DNP_Beta(1:nprec)
    CASE('NEUT_VELO')
      READ(oneline, *) astring, Neut_Velo(1:ng)
    CASE('DNEUT_CHI')
      READ(oneline, *) astring, Dneut_chi(1:ng)
  END SELECT
END DO

299 continue

DO iiso = 1, nmod
    CALL inflow_BenchXS_CEA(modiso(iiso), fiso)
END DO 
IF(linflow) deallocate(modiso)

END SUBROUTINE


SUBROUTINE alloc_BenchXS(ngben, nXslType, scatod)
!Allocate BenchmarkXS
use param
use allocs
use BenchXS,   only : MacXsBen, benchxstype
INTEGER, intent(in) :: ngben, nXslType
!logical, intent(in) :: lscat1
INTEGER :: scatod
INTEGER :: ixsl, i
allocate(MacXsBen(nXslType))

DO ixsl = 1, nXslType
  MacXsBen(ixsl)%lempty = TRUE;     MacXsBen(ixsl)%lFuel = FALSE
  CALL Dmalloc(MacXsBen(ixsl)%xst, ngben);  CALL Dmalloc(MacXsBen(ixsl)%xstr, ngben);    !Allocate
  CALL Dmalloc(MacXsBen(ixsl)%xskf, ngben); CALL Dmalloc(MacXsBen(ixsl)%xsnf, ngben);
  CALL Dmalloc(MacXsBen(ixsl)%chi, ngben);  CALL Dmalloc(MacXsBen(ixsl)%xsa, ngben);
  CALL Dmalloc(MacXsBen(ixsl)%xs0sum, ngben);
  CALL Dmalloc(MacXsBen(ixsl)%xs1sum, ngben);

  if(benchxstype.eq.1)THEN
    CALL Dmalloc(MacXsBen(ixsl)%xss, ngben, ngben);
    MacXsBen(ixsl)%Xss0 => MacXsBen(ixsl)%Xss
    IF(scatod .ge. 1)THEN
        CALL Dmalloc(MacXsBen(ixsl)%xss1, ngben, ngben)
        IF(scatod .ge. 2)THEN
            CALL Dmalloc(MacXsBen(ixsl)%xss2, ngben, ngben)
            IF(scatod .ge. 3)THEN
                CALL Dmalloc(MacXsBen(ixsl)%xss3, ngben, ngben)
            ENDIF
        ENDIF
    ENDIF
  ELSE
    CALL alloc_BenchXSSM(ngben, nXSltype, scatod)
  ENDIF
  !Transient Variable
ENDDO
END SUBROUTINE

SUBROUTINE alloc_BenchXSSM(ngben, nXslType, scatod)
!Allocate BenchmarkXS
use param
use allocs
use BenchXS,   only : MacXsBen
INTEGER, intent(in) :: ngben, nXslType
!logical, intent(in) :: lscat1
INTEGER :: scatod
INTEGER :: ixsl, i

DO ixsl = 1, nXslType
    ALLOCATE(MacXsBen(ixsl)%PnSM(scatod+1,ngben))
ENDDO

END SUBROUTINE

SUBROUTINE Alloc_BenchKinParam(nProc, ng)
USE PARAM
USE ALLOCS
USE BENCHXS,     ONLY : lKinParam,        DNP_NGRP,     DNP_BETA,    &
                        DNP_Lambda,       Neut_velo,    Dneut_Chi
IMPLICIT NONE
INTEGER :: nProc, ng
DNP_NGRP = nProc
CALL DMALLOC(DNP_BETA, nProc)
CALL DMALLOC(DNP_Lambda, nProc)
CALL DMALLOC(Neut_Velo, ng)
CALL DMALLOC(Dneut_Chi, ng)

END SUBROUTINE

    
SUBROUTINE inflow_BenchXS( modiso, fiso)
USE BenchXs, ONLY : MacXSBen, ngben
USE InflowMathUtil
IMPLICIT NONE 
INTEGER ::  modiso, fiso
REAL :: bckl, tol, err
REAL, ALLOCATABLE :: phi(:,:) , phiOld(:,:), L(:), U(:), A(:,:,:), uchi(:), b(:,:), diffPhi(:,:)
REAL :: trcorr

INTEGER :: norder, niter
INTEGER :: io, ig, jg, iiso

bckl = 1.e-4
tol=1.e-12
    
norder=1

allocate(phi(ngben,-1:norder+1),phiOld(ngben,0:norder),diffPhi(ngben,0:norder))
allocate(L(0:norder), U(0:norder),A(ngben,ngben,0:norder))
allocate(uchi(ngben))

! Lower triangular matrix of A 
L=0.
do io=1,norder
    L(io)=real(io,8)/(2.*real(io,8)+1.)
    if(mod(io,2).eq.1) L(io)=-L(io)
end do 
    
! Upper triangular matrix of A 
U=0.
do io=0,norder-1
    U(io)=(real(io,8)+1.)/(2.*real(io,8)+1.)
    if(mod(io,2).eq.1) U(io)=-U(io)
end do 
L=L*bckl
U=u*bckl

uchi = MacXSben(fiso)%chi

allocate(b(ngben,0:norder))
b=0.
b(1:ngben,0)=uchi

! Set A in Linear System Ax=b
A=0. 
phi=0.
phi(:,0)=uchi

DO ig = 1, ngben
    A(ig, ig,0) = MacXSben(modiso)%xst(ig)   
    
    DO jg = 1, ngben
        A(ig, jg,0) = A(ig, jg,0) - MacXSben(modiso)%xss(jg,ig)
    END DO 
    
    A(ig, ig,1) = MacXSben(modiso)%xst(ig)  
    
    DO jg = 1, ngben
        A(ig, jg,1) = A(ig, jg,1) - MacXSben(modiso)%xss1(jg,ig)
    END DO 
    
END DO

DO io=0,norder
    CALL inverseA(A(:,:,io),ngben)
END DO 

err=1.
phiold(:,0:norder)=phi(:,0:norder)
niter=1
do while(err.ge.tol)
    do io=0,norder
        call matvecprod(A(:,:,io),b(:,io)-L(io)*phi(:,io-1)-U(io)*phi(:,io+1),phi(:,io),ngben)
    end do
            
    diffPhi=(phi(:,0:norder)-phiold(:,0:norder))/phiold(:,0:norder)
    err=rmserr(diffPhi,1,ngben,0,norder)
    phiold(:,0:norder)=phi(:,0:norder)
    !write(*,'(2x,a7,i4,a11,es15.5)') 'Iter : ',niter,' / Error : ',err
            
    niter=niter+1
end do 

! Transport Correction
   
DO ig=1,ngben    
    trcorr = 0.
    DO jg = 1, ngben
        trcorr = trcorr - MacXsben(iiso)%xss1(jg,ig) * phi(jg,1)        
    END DO        
    trcorr = trcorr / phi(ig,1)
    MacXsben(iiso)%xstr(ig) = MacXsben(iiso)%xst(ig) + trcorr        
END DO    



END SUBROUTINE


SUBROUTINE inflow_BenchXS_CEA(modiso, fiso)

USE BenchXs, ONLY : MacXSBen, ngben
USE InflowMathUtil
IMPLICIT NONE 
INTEGER :: modiso, fiso
REAL :: bckl, tol, err
REAL, ALLOCATABLE :: phi(:,:) , phiOld(:,:), L(:), U(:), A(:,:,:), uchi(:), b(:,:), diffPhi(:,:)
REAL :: trcorr

INTEGER :: norder, niter
INTEGER :: io, ig, jg

bckl = 1.e-4
tol=1.e-12
    
norder=1

allocate(phi(ngben,-1:norder+1),phiOld(ngben,0:norder),diffPhi(ngben,0:norder))
allocate(L(0:norder), U(0:norder),A(ngben,ngben,0:norder))
allocate(uchi(ngben))

! Lower triangular matrix of A 
L=0.
do io=1,norder
    L(io)=real(io,8)/(2.*real(io,8)+1.)
    if(mod(io,2).eq.1) L(io)=-L(io)
end do 
    
! Upper triangular matrix of A 
U=0.
do io=0,norder-1
    U(io)=(real(io,8)+1.)/(2.*real(io,8)+1.)
    if(mod(io,2).eq.1) U(io)=-U(io)
end do 
L=L*bckl
U=u*bckl

uchi = MacXSben(fiso)%chi

allocate(b(ngben,0:norder))
b=0.
b(1:ngben,0)=uchi

! Set A in Linear System Ax=b
A=0. 
phi=0.
phi(:,0)=uchi

DO ig = 1, ngben
    A(ig, ig,0) = MacXSben(modiso)%xst(ig)   
    
    
    DO jg = MacXSben(modiso)%PnSm(1,ig)%ib , MacXSben(modiso)%PnSm(1,ig)%ie
        A(ig, jg,0) = A(ig, jg,0) - MacXSben(modiso)%PnSm(1,ig)%from(jg)
    END DO 
    
    A(ig, ig,1) = MacXSben(modiso)%xst(ig)  
    
    DO jg = MacXSben(modiso)%PnSm(2,ig)%ib , MacXSben(modiso)%PnSm(2,ig)%ie
        A(ig, jg,1) = A(ig, jg,1) - MacXSben(modiso)%PnSm(2,ig)%from(jg)
    END DO 
    
END DO

DO io=0,norder
    CALL inverseA(A(:,:,io),ngben)
END DO 

err=1.
phiold(:,0:norder)=phi(:,0:norder)
niter=1
do while(err.ge.tol)
    do io=0,norder
        call matvecprod(A(:,:,io),b(:,io)-L(io)*phi(:,io-1)-U(io)*phi(:,io+1),phi(:,io),ngben)
    end do
            
    diffPhi=(phi(:,0:norder)-phiold(:,0:norder))/phiold(:,0:norder)
    err=rmserr(diffPhi,1,ngben,0,norder)
    phiold(:,0:norder)=phi(:,0:norder)
    !write(*,'(2x,a7,i4,a11,es15.5)') 'Iter : ',niter,' / Error : ',err
            
    niter=niter+1
end do 

! Transport Correction
    DO ig=1,ngben    
        trcorr = 0.
        DO jg = MacXSben(modiso)%PnSm(2,ig)%ib , MacXSben(modiso)%PnSm(2,ig)%ie
            trcorr = trcorr -  MacXSben(modiso)%PnSm(2,ig)%from(jg) * phi(jg,1)        
        END DO        
        trcorr = trcorr / phi(ig,1)
        MacXsben(modiso)%xstr(ig) = MacXsben(modiso)%xst(ig) + trcorr        
    END DO    



END SUBROUTINE