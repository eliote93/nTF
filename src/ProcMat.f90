#include <defines.h>
SUBROUTINE ProcMat(dataline, indev, outdev)
USE PARAM
USE TYPEDEF,       ONLY : Mixture_Type
USE Material_Mod,  ONLY : Mixture,    nMixType
USE CNTL,          ONLY : nTracerCntl
USE PE_MOD,        ONLY : PE
USE IOUTIL,        ONLY : toupper,   IFnumeric,    nfields,   fndchara,    nfieldto, &
                          CharSkip,  terminate,    IfNumeric_except
USE NuclidMap_mod, ONLY : lExistIsotope, AtomicWeight                          
USE XsLib_Mod,     ONLY : mapnucl,   ldiso
USE Boron_mod,     ONLY : b10frac0
USE ALLOCS

IMPLICIT NONE
character(256),intent(in) :: dataline
character(512) :: dataline0, ReadLine
character(512) :: MixtureInfo(200),dum
INTEGER :: indev, outdev

TYPE(Mixture_Type) :: bufMixture

INTEGER :: ipos(5)
INTEGER :: ndata, nspt, imix, niso, idiso, nsphiso, id
INTEGER :: i, j, k, n, m, jfr, jend
REAL :: fsum, fchk, aw
REAL :: natB, fwt
LOGICAL :: lRes,lAg,lIn,lCd
LOGICAL, SAVE :: lfirst
DATA lfirst /.TRUE./

IF(lfirst) THEN
  ALLOCATE(Mixture(nMixType))
  lfirst = FALSE
ENDIF

dataline0 = dataline
ndata = nfieldto(dataline, SLASH)
READ(dataline, *) imix

IF(ndata .LE. 3 .OR. ndata .GE. 6) THEN
  WRITE(MESG, '(a, i3)') '***  INPUT ERROR at MIXTURE ', imix
ENDIF

CALL fndchara(dataline, ipos, nspt, SLASH)

!Mixture information 
Mixture(imix)%lempty = .FALSE.
Mixture(imix)%lMox = .FALSE.
Mixture(imix)%lAIC = .FALSE.
READ(dataline, *) i, Mixture(imix)%name, j, Mixture(imix)%dens
Mixture(imix)%temp = nTracerCntl%TempInlet + ckelvin           !Mixture Temperature
IF(ndata .eq. 5) THEN
  READ(dataline, *) i, Mixture(imix)%name, j, Mixture(imix)%dens, Mixture(imix)%temp
  Mixture(imix)%temp = Mixture(imix)%temp + CKELVIN
ENDIF


Mixture(imix)%deplopt=j;
IF(j .ge. 1) THEN
  IF(j .EQ. 2) Mixture(imix)%lfuel = TRUE        !Fuel
  IF(j .GE. 2) Mixture(imix)%ldepl = TRUE        !Depletable
ENDIF

write(MixtureInfo(1), *) ''
CALL CharSkip(dataline, MixtureInfo(1), SLASH, 1)

n = 1
ndata = 0
ndata = nfields(MixtureInfo(1))
DO WHILE(TRUE)
  n = n +1
  READ(indev, '(a256)') MixtureInfo(n)
  IF(PE%Master) WRITE(outdev, '(A)')  trim(MixtureInfo(n))
  IF(.NOT. IfNumeric_except(MixtureInfo(n),35)) EXIT  ! IF isotope number starts with '#' (which is for SPH factor), it is regarded as just isotope number.
  m = nfields(MixtureInfo(n))
  ndata = ndata + nfields(MixtureInfo(n))
ENDDO
BACKSPACE(indev); if(PE%Master) BACKSPACE(outdev)
IF(ndata .GT. 1 .AND. MOD(ndata, 2) .NE. 0) THEN
  WRITE(mesg, '(a, i5)') '(Isotope, fraction) pair in error for mixture : ', imix
  CALL terminate(mesg)
ENDIF

niso = ndata/2; Mixture(imix)%niso = niso
CALL dmalloc(Mixture(imix)%idiso, niso)
CALL dmalloc(Mixture(imix)%fweig, niso)
CALL dmalloc(Mixture(imix)%pnum, niso)
CALL dmalloc(Mixture(imix)%lSSPH, niso)
jfr = 0
jend = 0
lAg=.false.; lIn=.false.; lCd=.false.
DO i = 1, n-1
  k = nfields(MixtureInfo(i))/2
  jfr = jend + 1; jend = jfr + k-1
  dum=adjustl(MixtureInfo(i))
  IF (ichar(dum(1:1)).eq.35) THEN
      Mixture(imix)%lSSPH(jfr:jend)=.TRUE.
      Read(dum(2:), *) (Mixture(imix)%idiso(j), Mixture(imix)%fweig(j), j = jfr, jend)
  ELSE
      Mixture(imix)%lSSPH(jfr:jend)=.FALSE.
      Read(dum, *) (Mixture(imix)%idiso(j), Mixture(imix)%fweig(j), j = jfr, jend)
  ENDIF
  DO j=jfr, jend
     idiso = Mixture(imix)%idiso(j)
     IF (idiso.eq.92238) THEN
         Mixture(imix)%lSSPH(j)=.TRUE.
     ELSE
         idiso = idiso/1000
         if (idiso.eq.47) lAg=.true.
         if (idiso.eq.48) lIn=.true.
         if (idiso.eq.49) lCd=.true.
     ENDIF
  ENDDO  
ENDDO
if (lAg.and.lIn.and.lCd) Mixture(imix)%lAIC=.true.
CONTINUE
!Isotope Existence Check
fsum = 0
DO i = 1, niso
  idiso = Mixture(imix)%idiso(i)
  IF(.NOT. lExistIsotope(idiso)) THEN
    WRITE(mesg, *) 'This nuclide is not included in the library',Mixture(imix)%idiso(i)
    CALL Terminate(mesg)
  ENDIF
  fsum = fsum + Mixture(imix)%fweig(i)
ENDDO

IF(Mixture(imix)%name .EQ. 'MOD') THEN
  Mixture(imix)%lh2o = TRUE
ENDIF

IF(Mixture(imix)%name .EQ. 'CLD') THEN
  Mixture(imix)%lCLD = TRUE
ENDIF

IF(Mixture(imix)%name .EQ. 'MOX') THEN
  Mixture(imix)%lMOX = TRUE
ENDIF

IF(Mixture(imix)%name .EQ. 'AIC') THEN
  Mixture(imix)%lAIC = TRUE
ENDIF
fchk = abs(fsum - 100._8)

IF(fchk .lt. 1e-10_8) THEN !RATIO
  DO k = 1, niso
    idiso = Mixture(imix)%idiso(K)
    aw = AtomicWeight(idiso)
    Mixture(imix)%pnum(k) = Mixture(imix)%dens * Mixture(imix)%fweig(k) * AVOGADRO / aw / 100.0_8
  ENDDO
ELSEIF(fsum .lt. 1._8) THEN !Number Density
  DO k = 1, niso
    Mixture(imix)%pnum(k) = Mixture(imix)%fweig(k)
  ENDDO
ELSE ! Wrong Ratio Case
    WRITE(mesg, *) 'Ratio of this Mixture is not correct', imix
    CALL Terminate(mesg)
ENDIF
!Gd
DO i = 1, niso
  idiso = Mixture(imix)%idiso(i)
  idiso = idiso/100
  IF(idiso == 641) Mixture(imix)%lGd = .TRUE.    
ENDDO
!Nat. Boron Conversion
IF (Mixture(imix)%lh2o) THEN
  DO i = 1, niso
    idiso = Mixture(imix)%idiso(i)
    IF (idiso .EQ. 5000) THEN
      ALLOCATE(bufMixture%idiso(niso),bufMixture%pnum(niso),bufMixture%fweig(niso),bufMixture%lSSPH(niso))
      bufMixture%idiso = Mixture(imix)%idiso
      bufMixture%pnum = Mixture(imix)%pnum
      bufMixture%fweig = Mixture(imix)%fweig
      bufMixture%lSSPH = Mixture(imix)%lSSPH
      DEALLOCATE(Mixture(imix)%idiso,Mixture(imix)%pnum,Mixture(imix)%fweig,Mixture(imix)%lSSPH)
      
      Mixture(imix)%niso = niso+1
      ALLOCATE(Mixture(imix)%idiso(niso+1),Mixture(imix)%pnum(niso+1),Mixture(imix)%fweig(niso+1),Mixture(imix)%lSSPH(niso+1))
      DO k = 1, i-1
        Mixture(imix)%idiso(k) = bufMixture%idiso(k)
        Mixture(imix)%pnum(k) = bufMixture%pnum(k)
        Mixture(imix)%fweig(k) = bufMixture%fweig(k)
        Mixture(imix)%lSSPH(k) = bufMixture%lSSPH(k)
      END DO
      DO k = i+2, niso+1
        Mixture(imix)%idiso(k) = bufMixture%idiso(k-1)
        Mixture(imix)%pnum(k) = bufMixture%pnum(k-1)
        Mixture(imix)%fweig(k) = bufMixture%fweig(k-1)
        Mixture(imix)%lSSPH(k) = bufMixture%lSSPH(k-1)
      END DO
      Mixture(imix)%idiso(i) = 5010; Mixture(imix)%idiso(i+1) = 5011;
      natB = bufMixture%pnum(i)
      Mixture(imix)%pnum(i) = natB*b10frac0; Mixture(imix)%pnum(i+1) = natB*(1.-b10frac0)
      natB = bufMixture%fweig(i); fwt = AtomicWeight(5011)*(1.-b10frac0)/AtomicWeight(5010)/b10frac0
      Mixture(imix)%fweig(i) = natB/(1.+fwt); Mixture(imix)%fweig(i+1) = natB*fwt/(1.+fwt)
      Mixture(imix)%lSSPH(i) = bufMixture%lSSPH(i); Mixture(imix)%lSSPH(i+1) = bufMixture%lSSPH(i);
      
      DEALLOCATE(bufMixture%idiso, bufMixture%pnum, bufMixture%fweig, bufMixture%lSSPH)
      EXIT
    END IF
  END DO  
END IF
  
! Check Resonance
DO i=1, niso
  idiso = Mixture(imix)%idiso(i)
  id = mapnucl(idiso)  
  IF (ldiso(id)%lreso) Mixture(imix)%lres=.TRUE.
ENDDO

END SUBROUTINE

!--- JSR Edit : nTIG Restart

SUBROUTINE ProcMatBin
USE PARAM
USE TYPEDEF,       ONLY : Mixture_Type
USE Material_Mod,  ONLY : Mixture,    nMixType
USE CNTL,          ONLY : nTracerCntl
USE PE_MOD,        ONLY : PE
USE IOUTIL,        ONLY : toupper,   IFnumeric,    nfields,   fndchara,    nfieldto, &
                          CharSkip,  terminate,    IfNumeric_except
USE NuclidMap_mod, ONLY : lExistIsotope, AtomicWeight                          
USE XsLib_Mod,     ONLY : mapnucl,   ldiso
USE ALLOCS

IMPLICIT NONE
character(512) :: dataline0, ReadLine
character(512) :: MixtureInfo(200),dum
INTEGER :: indev, outdev

INTEGER :: ipos(5)
INTEGER :: matfid = 28, Nmix, nameid, imixture, matmat= 35
INTEGER :: ndata, nspt, imix, niso, idiso, nsphiso, id
INTEGER :: i, j, k, n, m, jfr, jend
REAL :: fsum, fchk, aw
LOGICAL :: lRes 
LOGICAL, SAVE :: lfirst
DATA lfirst /.TRUE./

IF(lfirst) THEN
  ALLOCATE(Mixture(nMixType))
  lfirst = FALSE
ENDIF

!Mixture information 
OPEN(matfid, file = "MATERIAL.bin", form = 'unformatted')
READ(matfid), Nmix

ALLOCATE(Mixture(Nmix))

DO imix = 1, Nmix
  Mixture(imix)%lempty = .FALSE.
  Mixture(imix)%lMox = .FALSE.
  READ(matfid), Niso,  imixture
  Mixture(imix)%niso = niso
  READ(matfid), nameid, j, Mixture(imix)%dens, Mixture(imix)%temp  ! 180530 %nameid
  Mixture(imix)%temp = Mixture(imix)%temp + CKELVIN
  ! ENDIF
  Mixture(imix)%name = "undefined"
  !------- 180530 %name
  SELECT CASE(nameid)
  CASE(1)
      Mixture(imix)%name = "UO2"
  CASE(2)
      Mixture(imix)%name = "MOD"
  CASE(3)
      Mixture(imix)%name = "CLD"
  CASE(4)
      Mixture(imix)%name = "MOX"
  END SELECT
  ! ------------------
  Mixture(imix)%deplopt=j;
  IF(j .ge. 1) THEN
    IF(j .EQ. 2) Mixture(imix)%lfuel = TRUE        !Fuel
    IF(j .GE. 2) Mixture(imix)%ldepl = TRUE        !Depletable
  ENDIF

  CALL dmalloc(Mixture(imix)%idiso, niso)
  CALL dmalloc(Mixture(imix)%fweig, niso)
  CALL dmalloc(Mixture(imix)%pnum, niso)
  CALL dmalloc(Mixture(imix)%lSSPH, niso)

  !--------------- 180530 edit
  DO i = 1, Niso
      READ(matfid), Mixture(imix)%idiso(i), Mixture(imix)%fweig(i)
      IF (Mixture(imix)%idiso(i).eq.92238) Mixture(imix)%lSSPH(i)=.TRUE.
  END DO
  !--------------------------------------
  
  !Isotope Existence Check
  fsum = 0
  DO i = 1, niso
    idiso = Mixture(imix)%idiso(i)
    IF(.NOT. lExistIsotope(idiso)) THEN
      WRITE(mesg, *) 'This nuclide is not included in the library',Mixture(imix)%idiso(i)
      CALL Terminate(mesg)
    ENDIF
    fsum = fsum + Mixture(imix)%fweig(i)
  ENDDO
  
  IF(Mixture(imix)%name .EQ. 'MOD') THEN
    Mixture(imix)%lh2o = TRUE
  ENDIF
  
  IF(Mixture(imix)%name .EQ. 'CLD') THEN
    Mixture(imix)%lCLD = TRUE
  ENDIF
  
  IF(Mixture(imix)%name .EQ. 'MOX') THEN
    Mixture(imix)%lMOX = TRUE
  ENDIF
  
  fchk = abs(fsum - 100._8)
  
  IF(fchk .lt. 1e-10_8) THEN !RATIO
    DO k = 1, niso
      idiso = Mixture(imix)%idiso(K)
      aw = AtomicWeight(idiso)
      Mixture(imix)%pnum(k) = Mixture(imix)%dens * Mixture(imix)%fweig(k) * AVOGADRO / aw / 100.0_8
    ENDDO
  ELSEIF(fsum .lt. 1._8) THEN !Number Density
    DO k = 1, niso
      Mixture(imix)%pnum(k) = Mixture(imix)%fweig(k)
    ENDDO
  ELSE ! Wrong Ratio Case
      WRITE(mesg, *) 'Ratio of this Mixture is not correct', imix
      CALL Terminate(mesg)
  ENDIF
  !Gd
  DO i = 1, niso
    idiso = Mixture(imix)%idiso(i)
    idiso = idiso/1000
    IF(idiso == 64) Mixture(imix)%lGd = .TRUE.    
  ENDDO
  
  ! Check Resonance
  DO i=1, niso
    idiso = Mixture(imix)%idiso(i)
    id = mapnucl(idiso)  
    IF (ldiso(id)%lreso) Mixture(imix)%lres=.TRUE.
  ENDDO
END DO
END SUBROUTINE