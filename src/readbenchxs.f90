SUBROUTINE ReadBenchXs(indev,ioutdev, ng,  nPrec, nXsChange, ScatOd)
!Reading Benchmark XS
USE param
USE cntl,           only : nTracerCntl
USE TRAN_MOD,       ONLY : TranCntl
!USE PE_Mod,         ONLY : Master
USE files,          only : filename,    InputFileIdx,  io5,            io8
USE inputcards ,    only : oneline,     probe,         mxcard,         nblock !,      &
                         !  FindBlockId, FindCardId,    cards,          blocks
use ioutil,         only : terminate,   toupper,       openfile,       IFnumeric,     &
                           nfields,     message
use BenchXs,        only : nxsltype,    ngben,         MacXsBen,       lKinParam,     &
                           nxschangeben,  benchxstype,                                &
                           DNP_BETA,    DNP_Lambda,    Neut_Velo,      Dneut_Chi,     &
                           KinParamXsBen
IMPLICIT NONE

INTEGER, intent(in) :: indev,ioutdev           !Input and Output Device Unit number
INTEGER :: nPreC                               !Number of Delayed Neutron Precursor Group
INTEGER :: ng                                  !Number of Energy Group
INTEGER :: nXsChange                           !Number of XS change for transient
INTEGER :: nXsNoise
!LOGICAL :: lscat1                              !P1 Scattering Option Flag
INTEGER :: ScatOd



INTEGER :: nLineField
INTEGER :: ixsl, iso                           !Benchmark XS index
INTEGER :: i
INTEGER :: ig, igs
real :: sumxss, sumchi
REAL :: temp, temp2
CHARACTER(15)   :: cardname,  astring
CHARACTER(1024) :: XsOneLine

INTEGER :: infinp
LOGICAL :: linflow = .FALSE.
LOGICAL :: lmod = .FALSE.
INTEGER :: nmod, fiso, iiso
INTEGER, ALLOCATABLE :: modiso(:)


IF(nTracerCntl%lDynamicBen) THEN
  CALL ReadBenchXS_Dynamic(indev, ng, nprec, nXsChange, ScatOd) 
  RETURN
END IF

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

nXsNoise = TranCntl%nNoise
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
CALL alloc_BenchXs(ngben, nXslType + nXsChange + nXsNoise, scatod)
IF(lKinParam) CALL Alloc_BenchKinParam(nPrec, ng, nXslType + nXsChange + nXsNoise)
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

  IF(TranCntl%lKineticBen .AND. lKinParam) THEN
    READ(indev, '(a1024)') XsOneline
    READ(Xsoneline, *)  (KinParamXsBen(iso)%Beta(igs), igs = 1, nprec)
    READ(indev, '(a1024)') XsOneline
    READ(Xsoneline, *)  (KinParamXsBen(iso)%Velo(igs), igs = 1, ngben)
    DO ig = 1, ngben
      READ(indev, '(a1024)') XsOneline
      READ(Xsoneline, *) (KinParamXsBen(iso)%ChiDg(ig, igs), igs = 1, nprec)
    ENDDO
  END IF

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

IF(lKinParam .AND. TranCntl%lKineticBen) THEN
  DO iso = 1, nXslType
    !KinParamXsBen(iso)%Beta(:) = 1E-4
    DO ig = 1, ngben
      temp = 0._8
      temp2 = 0._8
      DO igs = 1, nprec
        temp = temp + KinParamXsben(iso)%ChiDg(ig, igs) * KinParamXsben(iso)%Beta(igs)
        temp2 = temp2 + KinParamXsben(iso)%Beta(igs)
      ENDDO
      KinParamXsben(iso)%ChiD(ig) = temp / temp2
    ENDDO
  ENDDO
  Dneut_chi = KinParamXsBen(2)%ChiD
  !Dneut_chi = KinParamXsBen)%ChiD
END IF

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
IF(lKinParam) CALL Alloc_BenchKinParam(nPrec, ng, nXslType)
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

SUBROUTINE Alloc_BenchKinParam(nProc, ng, nXslType)
USE PARAM
USE TRAN_MOD,    ONLY : TranCntl
USE ALLOCS
USE BENCHXS,     ONLY : lKinParam,        DNP_NGRP,     DNP_BETA,    &
                        DNP_Lambda,       Neut_velo,    Dneut_Chi,   &
                        KinParamXsBen
IMPLICIT NONE
INTEGER :: nProc, ng, nXslType

INTEGER :: ixsl

DNP_NGRP = nProc
CALL DMALLOC(DNP_BETA, nProc)
CALL DMALLOC(DNP_Lambda, nProc)
CALL DMALLOC(Neut_Velo, ng)
CALL DMALLOC(Dneut_Chi, ng)

IF(TranCntl%lKineticBen) THEN
  ALLOCATE(KinParamXsBen(nXslType))
  DO ixsl = 1, nXslType
    CALL Dmalloc(KinParamXsBen(ixsl)%Beta, nProc)
    CALL Dmalloc(KinParamXsBen(ixsl)%Velo, ng)
    CALL Dmalloc(KinParamXsBen(ixsl)%ChiD, ng)
    CALL Dmalloc(KinParamXsBen(ixsl)%ChiDg, ng, nProc)
  END DO
END IF

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

SUBROUTINE ReadBenchXS_Dynamic(indev, ng, nprec, nXsChange, ScatOd)
USE PARAM
USE CNTL,           ONLY : nTracerCntl
USE TRAN_MOD,       ONLY : TranCntl
USE files,          ONLY : filename,    InputFileIdx,  io5,            io8
USE inputcards ,    ONLY : oneline,     probe,         mxcard,         nblock 
USE ioutil,         ONLY : terminate,   toupper,       openfile,       IFnumeric,     &
                           nfields,     message
USE BenchXs,        ONLY : nxsltype,    lKinParam,     ngben,          nXsChangeBen,  &
                           benchxstype, DynMacXsBen,   DNP_Lambda,     Dneut_Chi ,    &
                           alloc_DynBenchXs
IMPLICIT NONE
INTEGER, INTENT(IN) :: indev
INTEGER :: nprec, ng, nXsChange, ScatOd

REAL :: sumchi, sumxss
INTEGER, PARAMETER :: nominal_xsl = 2, nominal_temp = 2
INTEGER :: nXsNoise
INTEGER :: nLineField, libod
INTEGER :: ixsl, itemp, ig, igs, iprec, iod
CHARACTER(15) :: cardname, astring, astring1, astring2
CHARACTER(1024) :: XsOneLine

nXsNoise = TranCntl%nNoise
benchxstype = 1
ngben = ng 
nXsChangeBen = nXsChange
nxsltype = 0
libod = 0
DO WHILE(TRUE)
  READ(indev, '(a256)', END = 199) oneline  
  IF(probe.eq.BANG) CYCLE;     IF(probe.eq.POUND) CYCLE
  IF(oneline.eq.BLANK) CYCLE;  IF(IFnumeric(oneline)) CYCLE
  IF(probe .ne. BLANK) CALL terminate('First Column should be blank:'//trim(oneline))
  READ(oneline, *) cardname
  CALL toupper(cardname)
  nLineField = nfields(oneline) - 1
  IF(cardname .EQ. 'BASE_MICRO') THEN
    READ(oneline, *) astring, ixsl
    nXslType = max(nXslType, ixsl)
  END IF
  IF(cardname .EQ. 'DNP_NGRP') THEN
    READ(oneline, *) astring, nprec
    lKinParam = .TRUE.
  END IF
  IF(cardname .EQ. 'SCATOD') THEN
    READ(oneline, *) astring, iod
    libod = max(libod, iod) 
  END IF
END DO
199 CONTINUE
REWIND(indev)
ALLOCATE(DynMacXsBen(nXslType+nXsChange+nXsNoise))
CALL Alloc_BenchKinParam(nprec, ng, nxsltype+nXsChange+nXsNoise)

DO WHILE(.TRUE.)
  READ(indev, '(a256)', END = 299) oneline
  IF(probe.eq.BANG) cycle;     IF(probe.eq.POUND) cycle          !Cycle Condition(go to next line)
  IF(oneline.eq.BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe.eq.DOT) exit;       IF(probe.eq.SLASH) exit           !Exit Condition
  read(oneline,*) cardname
  CALL toupper(cardname)
  nLineField = nfields(oneline) - 1
  SELECT CASE(cardname)
  CASE('BASE_MICRO')
    READ(oneline, *) astring, ixsl
    READ(indev, '(a256)') oneline
    READ(oneline, *) astring, DynMacXsBen(ixsl)%ntemp
    CALL alloc_DynBenchXS(ixsl, DynMacXsBen(ixsl)%ntemp, nprec)
    READ(indev, '(a256)') oneline
    READ(oneline, *) astring, DynMacXsBen(ixsl)%temp
    READ(indev, '(a256)') oneline
    DO ig = 1, ngben
      DO itemp = 1, DynMacXsBen(ixsl)%ntemp
        READ(indev, '(a256)') oneline
        READ(oneline, *) astring, astring1, DynMacXsBen(ixsl)%xst(itemp, ig), DynMacXsBen(ixsl)%xsa(itemp, ig), DynMacXsBen(ixsl)%xsnf(itemp, ig), &
                         DynMacXsBen(ixsl)%xskf(itemp, ig), DynMacXsBen(ixsl)%chi(itemp, ig)
        IF(DynMacXSBen(ixsl)%xsnf(itemp, ig) .GT. epsm4) DynMacXsBen(ixsl)%lFuel = .TRUE.
      END DO
    END DO
    DO iod = 0, libod
      READ(indev, '(a256)') oneline !--scattering
      READ(indev, '(a256)') oneline !--scattering
      IF(iod .EQ. 0) THEN
        DO ig = 1, ngben
          DO itemp = 1, DynMacXsBen(ixsl)%ntemp
            READ(indev, '(a1024)') XsOneline
            READ(XsOneLine, *) astring, astring1, (DynMacXsBen(ixsl)%xss0(itemp, igs, ig), igs = 1, ngben)
          END DO
        END DO
      ELSE IF(iod .EQ. 1 .AND. scatod .GE. 1) THEN
        DO ig = 1, ngben
          DO itemp = 1, DynMacXsBen(ixsl)%ntemp
            READ(indev, '(a1024)') XsOneline
            READ(XsOneLine, *) astring, astring1, (DynMacXsBen(ixsl)%xss1(itemp, igs, ig), igs = 1, ngben)
          END DO
        END DO
        !DynMacXsBen(ixsl)%xss1 = DynMacXsBen(ixsl)%xss1 * 0.5
      ELSE IF(iod .EQ. 2 .AND. scatod .GE. 2) THEN
        DO ig = 1, ngben
          DO itemp = 1, DynMacXsBen(ixsl)%ntemp
            READ(indev, '(a1024)') XsOneline
            READ(XsOneLine, *) astring, astring1, (DynMacXsBen(ixsl)%xss2(itemp, igs, ig), igs = 1, ngben)
          END DO
        END DO
      ELSE IF(iod .EQ. 3 .AND. scatod .GE. 3) THEN
        DO ig = 1, ngben
          DO itemp = 1, DynMacXsBen(ixsl)%ntemp
            READ(indev, '(a1024)') XsOneline
            READ(XsOneLine, *) astring, astring1, (DynMacXsBen(ixsl)%xss3(itemp, igs, ig), igs = 1, ngben)
          END DO
        END DO
      ELSE
        DO ig = 1, ngben
          DO itemp = 1, DynMacXsBen(ixsl)%ntemp
            READ(indev, '(a1024)') XsOneline
          END DO
        END DO
      END IF
    END DO
    !READ(indev, '(a256)') oneline !--scattering
    !DO ig = 1, ngben
    !  DO itemp = 1, DynMacXsBen(ixsl)%ntemp
    !    READ(indev, '(a1024)') XsOneline
    !    READ(XsOneLine, *) astring, astring1, (DynMacXsBen(ixsl)%xss0(itemp, igs, ig), igs = 1, ngben)
    !  END DO
    !END DO
    !IF(scatod .ge. 1) THEN
    !  DO ig = 1, ngben
    !    DO itemp = 1, DynMacXsBen(ixsl)%ntemp
    !      READ(indev, '(a1024)') XsOneline
    !      READ(XsOneLine, *) astring, astring1, (DynMacXsBen(ixsl)%xss1(itemp, igs, ig), igs = 1, ngben)
    !    END DO
    !  END DO
    !  IF(scatod .ge. 2) THEN
    !    DO ig = 1, ngben
    !      DO itemp = 1, DynMacXsBen(ixsl)%ntemp
    !        READ(indev, '(a1024)') XsOneline
    !        READ(XsOneLine, *) astring, astring1, (DynMacXsBen(ixsl)%xss2(itemp, igs, ig), igs = 1, ngben)
    !      END DO
    !    END DO
    !    IF(scatod .ge. 3) THEN
    !      DO ig = 1, ngben
    !        DO itemp = 1, DynMacXsBen(ixsl)%ntemp
    !          READ(indev, '(a1024)') XsOneline
    !          READ(XsOneLine, *) astring, astring1, (DynMacXsBen(ixsl)%xss3(itemp, igs, ig), igs = 1, ngben)
    !        END DO
    !      END DO
    !    END IF
    !  END IF
    !END IF
    IF(lKinParam) THEN
      READ(indev, '(a256)') oneline !--kinetic parameters
      DO itemp = 1, DynMacXsBen(ixsl)%ntemp
        READ(indev, '(a1024)') XsOneline
        READ(XsOneLine, *) astring, astring1, (DynMacXsBen(ixsl)%Beta(itemp, iprec), iprec = 1, nprec)
      END DO
      DO itemp = 1, DynMacXsBen(ixsl)%ntemp
        READ(indev, '(a1024)') XsOneline
        READ(XsOneLine, *) astring, astring1, (DynMacXsBen(ixsl)%Lambda(itemp, iprec), iprec = 1, nprec)
      END DO
      DO itemp = 1, DynMacXsBen(ixsl)%ntemp
        READ(indev, '(a1024)') XsOneline
        READ(XsOneLine, *) astring, astring1, (DynMacXsBen(ixsl)%Velo(itemp, ig), ig = 1, ngben)
      END DO
      DO itemp = 1, DynMacXsBen(ixsl)%ntemp
        READ(indev, '(a1024)') XsOneline
        READ(XsOneLine, *) astring, astring1, (DynMacXsBen(ixsl)%Chid(itemp, ig), ig = 1, ngben)
      END DO
    END IF
    
    IF(.NOT. DynMacXsBen(ixsl)%lFuel) THEN
      DynMacXsBen(ixsl)%lCR = TRUE
      DO itemp = 1, DynMacXsBen(ixsl)%ntemp
        IF(DynMacXsBen(ixsl)%xsa(itemp, ng) .LE. 1._8) THEN
          DynMacXsBen(ixsl)%lCR = FALSE
          EXIT
        END IF
      END DO
    ELSE
      DynMacXsBen(ixsl)%lCR = FALSE
    END IF
    !--chi verification
    DO itemp = 1, DynMacXsBen(ixsl)%ntemp
      sumchi = SUM(DynMacXsBen(ixsl)%chi(itemp, 1:ngben))
      IF(sumchi .GT. epsm3) THEN
        sumchi = 1. / sumchi
        DynMacXsBen(ixsl)%chi(itemp, 1:ngben) = sumchi * DynMacXsBen(ixsl)%chi(itemp, 1:ngben)
      END IF
    END DO

    !-- Transport XS
    IF(scatod .EQ. 0) THEN
      DO ig = 1, ngben
        DO itemp = 1, DynMacXsBen(ixsl)%ntemp
          sumXss = sum(DynMacXsBen(ixsl)%xss(itemp, ig, 1:ngben))
          DynMacXsBen(ixsl)%xst(itemp, ig) = DynMacXsBen(ixsl)%xsa(itemp, ig) + sumXss
          DynMacXsBen(ixsl)%xs0sum(itemp, ig) = sumXss
          DynMacXsBen(ixsl)%xstr(itemp, ig) = DynMacXsBen(ixsl)%xst(itemp, ig)
        END DO
      END DO
    ELSE
      DO ig = 1, ngben
        DO itemp = 1, DynMacXsBen(ixsl)%ntemp
          sumXss = sum(DynMacXsBen(ixsl)%xss(itemp, ig, 1:ngben))
          DynMacXsBen(ixsl)%xst(itemp, ig) = DynMacXsBen(ixsl)%xsa(itemp, ig) + sumXss
          DynMacXsBen(ixsl)%xs0sum(itemp, ig) = sumXss
          sumXss = sum(DynMacXsBen(ixsl)%xss1(itemp, ig, 1:ngben))
          DynMacXsBen(ixsl)%xs1sum(itemp, ig) = sumXss
        END DO
      END DO
      DO ig = 1, ngben
        DO itemp  = 1, DynMacXsBen(ixsl)%ntemp
          DynMacXsBen(ixsl)%xstr(itemp, ig) = DynMacXsBen(ixsl)%xst(itemp, ig) ! - DynMacXsBen(ixsl)%xs1sum(itemp, ig) 
        END DO
      END DO
    END IF
  END SELECT

END DO  
299 CONTINUE

IF(lKinParam) THEN
  Dneut_chi = DynMacXsBen(nominal_xsl)%Chid(nominal_temp, :)
  DNP_lambda = DynMacXsBen(nominal_xsl)%Lambda(nominal_temp, :)
END IF

END SUBROUTINE

SUBROUTINE ReadBenchXs_NEACRP(indev, ng, nprec, nxschange, scatod)
USE PARAM
USE CNTL,           ONLY : nTracerCntl
USE TRAN_MOD,       ONLY : TranCntl
USE files,          ONLY : filename,    InputFileIdx,  io5,            io8
USE inputcards ,    ONLY : oneline,     probe,         mxcard,         nblock 
USE ioutil,         ONLY : terminate,   toupper,       openfile,       IFnumeric,     &
                           nfields,     message
USE BenchXs,        ONLY : nxsltype,    lKinParam,     ngben,          nXsChangeBen,  &
                           benchxstype, NeaCrpXsBen,   MacXsBen,       nxsltype_total,&
                           DNP_Lambda,  DNP_BETA,      Neut_Velo,      Dneut_Chi,     &
                           NeaCrpCABen, SetBoronCoolant_NEACRP
IMPLICIT NONE
INTEGER, INTENT(IN) :: indev
INTEGER :: ng, nprec, nxschange, scatod

TYPE comp_rodded_type
  INTEGER :: icomp
  INTEGER :: nrodded
  INTEGER :: basecomp
  INTEGER, POINTER :: ca_idx(:)
  REAL, POINTER :: ca_wt(:)
END TYPE

TYPE(comp_rodded_type), POINTER :: comp_rodded(:)
REAL :: wt
INTEGER :: nLineField
INTEGER :: nXsNoise
INTEGER :: ncomp, nca, ncomp_rodded
INTEGER :: icomp, nrodded, ica, irod
INTEGER :: basecomp, i, j, ig
CHARACTER(15) :: cardname, astring

nXsNoise = TranCntl%nNoise
benchxstype = 1
ngben = ng
nXsChangeBen = nXsChange

!Scan library file
ncomp = 0
ncomp_rodded = 0
nca = 0
DO WHILE(TRUE)
  READ(indev, '(a256)', END = 199) oneline
  IF(probe.eq.BANG) CYCLE;     IF(probe.eq.POUND) CYCLE
  IF(oneline.eq.BLANK) CYCLE;  IF(IFnumeric(oneline)) CYCLE
  IF(probe .ne. BLANK) CALL terminate('First Column should be blank:'//trim(oneline))  
  READ(oneline,*) cardname
  CALL toupper(cardname)
  nLineField = nfields(oneline) - 1
  IF(cardname .EQ. 'COMP') THEN 
    READ(oneline, *) astring, icomp, nrodded
    ncomp = max(ncomp, icomp)
    IF(nrodded .GT. 0) ncomp_rodded = ncomp_rodded + 1 
  END IF
  IF(cardname .EQ. 'CA') THEN
    nca = nca + 1
  END IF
  IF(cardname .EQ.'DNP_NGRP') THEN
    read(oneline,*) astring,nprec   
    lKinParam = .TRUE.
  ENDIF
END DO
199 CONTINUE
REWIND(indev)

nxsltype = ncomp
nxsltype_total = ncomp + nXsChange + nXsNoise
!Allocate benchmark XS
CALL alloc_BenchXS(ngben, ncomp + nXsChange + nXsNoise, scatod)
IF(lKinParam) CALL Alloc_BenchKinParam(nPrec, ng, ncomp + nXsChange + nXsNoise)
ALLOCATE(NeaCrpXsBen(ncomp + nXsChange + nXsNoise))
ALLOCATE(NeaCrpCaBen(nca))
ALLOCATE(comp_rodded(ncomp_rodded))

!Read library file
irod = 0
DO WHILE(TRUE)
  READ(indev, '(a256)', END = 299) oneline
  IF(probe .EQ. BANG) CYCLE;     IF(probe .EQ. POUND) CYCLE          !Cycle Condition(go to next line)
  IF(oneline .EQ. BLANK) CYCLE;  IF(IFnumeric(oneline)) CYCLE
  IF(probe .EQ. DOT) EXIT;       IF(probe .EQ. SLASH) EXIT           !Exit Condition
  READ(oneline,*) cardname
  CALL toupper(cardname)
  nLineField = nfields(oneline) - 1
  SELECT CASE(cardname) 
  CASE('COMP')
    READ(oneline, *) astring, icomp, nrodded
    IF(nrodded .EQ. 0) THEN
      READ(indev, '(a256)') oneline
      READ(oneline, *) NeaCrpXsBen(icomp)%xs0(1:5)
      READ(indev, '(a256)') oneline
      READ(oneline, *) NeaCrpXsBen(icomp)%xs0(6:9), astring
      READ(indev, '(a256)') oneline
      READ(oneline, *) NeaCrpXsBen(icomp)%gradxs_c(1:5)
      READ(indev, '(a256)') oneline
      READ(oneline, *) NeaCrpXsBen(icomp)%gradxs_c(6:9), NeaCrpXsBen(icomp)%c0
      READ(indev, '(a256)') oneline
      READ(oneline, *) NeaCrpXsBen(icomp)%gradxs_tm(1:5)
      READ(indev, '(a256)') oneline
      READ(oneline, *) NeaCrpXsBen(icomp)%gradxs_tm(6:9), NeaCrpXsBen(icomp)%tm0
      NeaCrpXsBen(icomp)%tm0 = NeaCrpXsBen(icomp)%tm0 + CKELVIN
      READ(indev, '(a256)') oneline
      READ(oneline, *) NeaCrpXsBen(icomp)%gradxs_rho(1:5)
      READ(indev, '(a256)') oneline
      READ(oneline, *) NeaCrpXsBen(icomp)%gradxs_rho(6:9), NeaCrpXsBen(icomp)%rho0
      READ(indev, '(a256)') oneline
      READ(oneline, *) NeaCrpXsBen(icomp)%gradxs_td(1:5)
      READ(indev, '(a256)') oneline
      READ(oneline, *) NeaCrpXsBen(icomp)%gradxs_td(6:9), NeaCrpXsBen(icomp)%td0
      NeaCrpXsBen(icomp)%td0 = NeaCrpXsBen(icomp)%td0 + CKELVIN
    ELSE
      irod = irod + 1
      comp_rodded(irod)%icomp = icomp
      comp_rodded(irod)%nrodded = nrodded
      ALLOCATE(comp_rodded(irod)%ca_idx(nrodded))
      ALLOCATE(comp_rodded(irod)%ca_wt(nrodded))
      READ(indev, '(a256)') oneline
      READ(oneline, *) comp_rodded(irod)%basecomp, (comp_rodded(irod)%ca_idx(ica), comp_rodded(irod)%ca_wt(ica), ica = 1, nrodded)
    END IF
  CASE('CA')
    READ(oneline, *) astring, ica
    NeaCrpCABen(ica)%lCA = .TRUE.
    READ(indev, '(a256)') oneline
    READ(oneline, *) NeaCrpCABen(ica)%xs0(1:5)
    READ(indev, '(a256)') oneline
    READ(oneline, *) NeaCrpCABen(ica)%xs0(6:9), astring
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
299 CONTINUE

DO irod = 1, ncomp_rodded
  icomp = comp_rodded(irod)%icomp
  basecomp = comp_rodded(irod)%basecomp 
  NeaCrpXsBen(icomp)%c0 = NeaCrpXsBen(basecomp)%c0
  NeaCrpXsBen(icomp)%tm0 = NeaCrpXsBen(basecomp)%tm0
  NeaCrpXsBen(icomp)%rho0 = NeaCrpXsBen(basecomp)%rho0
  NeaCrpXsBen(icomp)%td0 = NeaCrpXsBen(basecomp)%td0
  NeaCrpXsBen(icomp)%xs0 = NeaCrpXsBen(basecomp)%xs0
  NeaCrpXsBen(icomp)%gradxs_c = NeaCrpXsBen(basecomp)%gradxs_c
  NeaCrpXsBen(icomp)%gradxs_tm = NeaCrpXsBen(basecomp)%gradxs_tm
  NeaCrpXsBen(icomp)%gradxs_rho = NeaCrpXsBen(basecomp)%gradxs_rho
  NeaCrpXsBen(icomp)%gradxs_td = NeaCrpXsBen(basecomp)%gradxs_td
  DO i = 1, comp_rodded(irod)%nrodded
    ica = comp_rodded(irod)%ca_idx(i)
    wt = comp_rodded(irod)%ca_wt(i)
    DO j = 1, 9
      NeaCrpXsBen(icomp)%xs0(j) = NeaCrpXsBen(icomp)%xs0(j) + wt * NeaCrpCABen(ica)%xs0(j)
    END DO 
  END DO
END DO

DO irod = 1, ncomp_rodded
  DEALLOCATE(comp_rodded(irod)%ca_idx)
  DEALLOCATE(comp_rodded(irod)%ca_wt)
END DO 
DEALLOCATE(comp_rodded)

DO icomp = 1, ncomp
  MacXsBen(icomp)%lempty = FALSE
  MacXsBen(icomp)%xstr(1) = NeaCrpXsBen(icomp)%xs0(1)
  MacXsBen(icomp)%xstr(2) = NeaCrpXsBen(icomp)%xs0(6)
  MacXsBen(icomp)%xsa(1) = NeaCrpXsBen(icomp)%xs0(3)
  MacXsBen(icomp)%xsa(2) = NeaCrpXsBen(icomp)%xs0(7)
  MacXsBen(icomp)%xsnf(1) = NeaCrpXsBen(icomp)%xs0(4)
  MacXsBen(icomp)%xsnf(2) = NeaCrpXsBen(icomp)%xs0(8)
  MacXsBen(icomp)%xskf(1) = NeaCrpXsBen(icomp)%xs0(5)
  MacXsBen(icomp)%xskf(2) = NeaCrpXsBen(icomp)%xs0(9)
  MacXsBen(icomp)%xss0(1,2) = NeaCrpXsBen(icomp)%xs0(2)
  MacXsBen(icomp)%xss0(1,1) = NeaCrpXsBen(icomp)%xs0(1) - NeaCrpXsBen(icomp)%xs0(2)- NeaCrpXsBen(icomp)%xs0(3)
  MacXsBen(icomp)%xss0(2,2) = NeaCrpXsBen(icomp)%xs0(6) - NeaCrpXsBen(icomp)%xs0(7)
  MacXsBen(icomp)%xss0(2,1) = 0.
  MacXsBen(icomp)%xst(1) = MacXsBen(icomp)%xstr(1)
  MacXsBen(icomp)%xst(2) = MacXsBen(icomp)%xstr(2)
  DO ig = 1, 2
    MacXsBen(icomp)%Xs0sum(ig) = sum(MacXsBen(icomp)%Xss(ig,1:ngben))
  END DO 
  MacXsBen(icomp)%chi(1) = 1.
  MacXsBen(icomp)%chi(0) = 0.
  IF(sum(MacXsBen(icomp)%xsnf(:)) .GT. epsm4) MacXsBen(icomp)%lfuel = .TRUE.
  MacXsBen(icomp)%boronppm = NeaCrpXsBen(icomp)%c0
END DO
CALL SetBoronCoolant_NEACRP(nTracerCntl%boronppm)

END SUBROUTINE

