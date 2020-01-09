SUBROUTINE P1ScatXsTreatment()
USE PARAM
USE XSLIB_MOD
USE XsUtil_mod
USE ALLOCS
IMPLICIT NONE
INTEGER :: ig, id,p1idx, it
INTEGER :: ntemp

DO id = 1, nelthel
  
  p1idx=mapnuclp13(ldiso(id)%NID)
  IF(p1idx .NE. 0) CYCLE
  ntemp =  ldiso(id)%ntemp
  DO ig = 1, noghel
    DO it = 1, ntemp
      ldiso(id)%sigs(ig, it) =ldiso(id)%sigstr(ig, it)
      ldiso(id)%sigss(ig,it) = ldiso(id)%sm(ig,it)%from(ig)
    ENDDO
  ENDDO
ENDDO
#ifdef decativeP12
DO id = 1, nelrhel
  p1idx=mapnuclp13(ldiso(id)%NID)
  IF(p1idx .eq. 0) CYCLE
  ntemp = ldiso(id)%np1temp
  DO it = 1 , ntemp
    DO ig = 1, noghel
      ldiso(id)%smp1(ig, it)%from(:) = 0
      ldiso(id)%smp2(ig, it)%from(:) = 0
      ldiso(id)%smp3(ig, it)%from(:) = 0
    ENDDO
  ENDDO
   ntemp = ldiso(id)%ntemp
  DO ig = 1, noghel
    ldiso(id)%sigs(ig, 1:ntemp) =ldiso(id)%sigstr(ig, 1:ntemp)
  ENDDO
ENDDO
#endif
END SUBROUTINE