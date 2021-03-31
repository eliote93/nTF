SUBROUTINE ApproxExp(PolarAng, npr)
  
USE PARAM
USE ALLOCS
USE TYPEDEF, ONLY : PolarAngle_Type
USE MOC_MOD, ONLY : expa, expb, expa_p, expb_p

IMPLICIT NONE

TYPE(PolarAngle_Type) :: PolarAng(npr)

INTEGER :: npr, ipr, i

REAL :: ttt, xval1, xval2, dx
REAL :: yval1(100), yval2(100), rsinvpol(100)
! ----------------------------------------------------

DO ipr = 1, npr
  rsinvpol(ipr) = one/PolarAng(ipr)%sinv
END DO

dx    = epsm3
xval1 = -40.0_8

DO ipr = 1, npr
  expa(-40000, ipr) = 0; expa_p(ipr, -40000) = 0
  expb(-40000, ipr) = 1; expb_p(ipr, -40000) = 1
  
  ttt = rsinvpol(ipr) * xval1
  
  IF (ttt .LT. -200.0_8) THEN
    yval1(ipr) = 1
  ELSE
    yval1(ipr) = 1 - exp(ttt)
  END IF
END DO

DO i = -39999, 0
  xval2 = xval1 + dx
  
  DO ipr = 1, npr
    ttt = rsinvpol(ipr) * xval2
    
    IF (ttt .LT. -100) THEN
      yval2(ipr) = 1
    ELSE
      yval2(ipr) = 1 - exp(ttt)
    END IF
    
    expa(i, ipr) = 1000._8 * (yval2(ipr) - yval1(ipr)); expa_p(ipr, i) = 1000._8 * (yval2(ipr) - yval1(ipr))
    expb(i, ipr) = yval1(ipr) - expa(i, ipr) * xval1;   expb_p(ipr, i) = yval1(ipr) - expa_p(ipr, i) * xval1
    
    yval1(ipr) = yval2(ipr)
  END DO
  
  xval1 = xval2
END DO

CONTINUE

DO ipr = 1, npr
  DO i = -39999,0
    expa(i, ipr) = dx * expa(i, ipr); expa_p(ipr, i) = dx * expa_p(ipr, i)
  END DO
END DO 
! ----------------------------------------------------

END SUBROUTINE ApproxExp