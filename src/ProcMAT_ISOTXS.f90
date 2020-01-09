#include <defines.h>
SUBROUTINE ProcMat_ISOTXS(dataline, indev, outdev)
USE PARAM
USE TYPEDEF,       ONLY : Mixture_Type
USE Material_Mod,  ONLY : Mixture,    nMixType
USE CNTL,          ONLY : nTracerCntl
USE PE_MOD,        ONLY : PE
USE IOUTIL,        ONLY : toupper,   IFnumeric,    nfields,   fndchara,    nfieldto, &
                          CharSkip,  terminate
USE NuclidMap_mod, ONLY : lExistIsotope, AtomicWeight                          
USE XSLIB_MOD,     ONLY : ldiso,nelthel,mapnucl
USE ALLOCS

IMPLICIT NONE
character(256),intent(in) :: dataline
character(512) :: dataline0, ReadLine
character(512) :: MixtureInfo(200)
character(6),pointer :: aid(:)
INTEGER :: indev, outdev

INTEGER :: ipos(5)
INTEGER :: ndata, nspt, imix, niso, idiso
INTEGER :: i, j, k, n, m, jfr, jend
REAL :: fsum, fchk, aw
LOGICAL :: lRes, lExist
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
READ(dataline, *) i, Mixture(imix)%name, j, Mixture(imix)%dens
Mixture(imix)%temp = nTracerCntl%TempInlet + ckelvin           !Mixture Temperature
IF(ndata .eq. 5) THEN
  READ(dataline, *) i, Mixture(imix)%name, j, Mixture(imix)%dens, Mixture(imix)%temp
  Mixture(imix)%temp = Mixture(imix)%temp + CKELVIN
ENDIF


IF(j .ge. 1) THEN
!  Mixture(imix)%lres = TRUE
!  nTracerCntl%lRes = TRUE
  IF(j .EQ. 2) Mixture(imix)%lfuel = TRUE        !Fuel
  IF(j .GE. 2) Mixture(imix)%ldepl = TRUE        !Depletable
ENDIF
Mixture(imix)%lres = .FALSE.



write(MixtureInfo(1), *) ''
CALL CharSkip(dataline, MixtureInfo(1), SLASH, 1)

n = 1
ndata = 0
ndata = nfields(MixtureInfo(1))
DO WHILE(TRUE)
  n = n +1
  READ(indev, '(a256)') MixtureInfo(n)
  IF(PE%Master) WRITE(outdev, '(a256)')  MixtureInfo(n)
  IF(.NOT. IfNumeric(MixtureInfo(n))) EXIT
  m = nfields(MixtureInfo(n))
  ndata = ndata + nfields(MixtureInfo(n))
ENDDO
BACKSPACE(indev); if(PE%Master) BACKSPACE(outdev)
IF(ndata .GT. 1 .AND. MOD(ndata, 3) .NE. 0) THEN
  WRITE(mesg, '(a, i5)') '(Isotope, fraction) pair in error for mixture : ', imix
  CALL terminate(mesg)
ENDIF

niso = ndata/3; Mixture(imix)%niso = niso
CALL dmalloc(Mixture(imix)%idiso, niso)
CALL dmalloc(Mixture(imix)%fweig, niso)
CALL dmalloc(Mixture(imix)%pnum, niso)
ALLOCATE(aid(niso))
jfr = 0
jend = 0
DO i = 1, n-1
  k = nfields(MixtureInfo(i))/3
  jfr = jend + 1; jend = jfr + k-1
  Read(MixtureInfo(i), *) (Mixture(imix)%idiso(j), aid(j), Mixture(imix)%fweig(j), j = jfr, jend)
ENDDO
CONTINUE
!Isotope Existence Check
fsum = 0
DO i = 1, niso
  lExist = .FALSE.
  DO j = 1, nelthel
    IF(aid(i) /= ldiso(j)%aid) CYCLE
    lExist = .TRUE.
    EXIT
  ENDDO
  IF(.NOT. lExist) THEN
    WRITE(mesg, *) 'This nuclide is not included in the library',aid(i)
    CALL Terminate(mesg)
  ENDIF
  fsum = fsum + Mixture(imix)%fweig(i)
  idiso = Mixture(imix)%idiso(i)
  idiso = idiso - mod(idiso,100) + mod(idiso,10)*10 + 1
  DO WHILE(.TRUE.)
    IF(.NOT. lExistIsotope(idiso)) EXIT
    IF(aid(i) == ldiso(mapnucl(idiso))%aid) EXIT
    idiso = idiso + 1
  ENDDO
  ldiso(j)%nid = idiso
  Mixture(imix)%idiso(i) = idiso
  mapnucl(idiso) = j
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

END SUBROUTINE