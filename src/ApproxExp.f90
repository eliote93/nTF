SUBROUTINE ApproxExp(PolarAng, npol)
  
USE PARAM,   ONLY : ONE, EPSM3
USE TYPEDEF, ONLY : PolarAngle_Type
USE MOC_MOD, ONLY : ExpA, ExpB

IMPLICIT NONE

TYPE(PolarAngle_Type) :: PolarAng(npol)

INTEGER :: npol, ipol, idx

REAL :: ttt, xval1, xval2, dx
REAL :: yval1(100), yval2(100), rsinvpol(100)
! ----------------------------------------------------

DO ipol = 1, npol
  rsinvpol(ipol) = ONE / PolarAng(ipol)%sinv
END DO

dx    = EPSM3
xval1 = -40.0_8

DO ipol = 1, npol
  ExpA(-40000, ipol) = 0
  ExpB(-40000, ipol) = 1
  
  ttt = rsinvpol(ipol) * xval1
  
  IF (ttt .LT. -200.0_8) THEN
    yval1(ipol) = 1
  ELSE
    yval1(ipol) = 1 - exp(ttt)
  END IF
END DO

DO idx = -39999, 0
  xval2 = xval1 + dx
  
  DO ipol = 1, npol
    ttt = rsinvpol(ipol) * xval2
    
    IF (ttt .LT. -100) THEN
      yval2(ipol) = 1
    ELSE
      yval2(ipol) = 1 - exp(ttt)
    END IF
    
    ExpA(idx, ipol) = 1000._8 * (yval2(ipol) - yval1(ipol))
    ExpB(idx, ipol) = yval1(ipol) - ExpA(idx, ipol) * xval1
    
    yval1(ipol) = yval2(ipol)
  END DO
  
  xval1 = xval2
END DO

CONTINUE

DO ipol = 1, npol
  DO idx = -39999,0
    ExpA(idx, ipol) = dx * ExpA(idx, ipol)
  END DO
END DO 
! ----------------------------------------------------

END SUBROUTINE ApproxExp