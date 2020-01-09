!Subroutine ApproxExp(PolarAng, npr, Expa, Expb)
Subroutine ApproxExp(PolarAng, npr)
USE PARAM
USE TYPEDEF, ONLY : PolarAngle_Type
USE ALLOCS
USE MOC_MOD, ONLY : expa, expb, expa_p, expb_p
IMPLICIT NONE
TYPE(PolarAngle_Type) :: PolarAng(npr)
!REAL  :: Expa(:, :), Expb(:, :)
INTEGER :: npr
INTEGER :: ipr,i
REAL :: ttt, xval1,xval2, dx
REAL :: yval1(100), yval2(100), rsinvpol(100)

!call dmalloc0(expa,-40000,0,1,npr)
!call dmalloc0(expb,-40000,0,1,npr)

DO ipr = 1, npr
  rsinvpol(ipr) = one/PolarAng(ipr)%sinv
ENDDO

dx=epsm3
xval1=-40.0_8
do ipr=1,npr
  expa(-40000,ipr)=0; expa_p(ipr,-40000)=0;
  expb(-40000,ipr)=1; expb_p(ipr,-40000)=1;
  ttt=rsinvpol(ipr)*xval1
  if(ttt.lt.-200.0_8) then
    yval1(ipr)=1
  else
    yval1(ipr)=1-exp(ttt)
  endif  
enddo
do i=-39999,0
  xval2=xval1+dx
  do ipr=1,npr
    ttt=rsinvpol(ipr)*xval2
    if(ttt.lt.-100) then
      yval2(ipr)=1
    else
      yval2(ipr)=1-exp(ttt)
    endif
    expa(i,ipr)=1000._8*(yval2(ipr)-yval1(ipr)); expa_p(ipr,i)=1000._8*(yval2(ipr)-yval1(ipr))
    expb(i,ipr)=yval1(ipr)-expa(i,ipr)*xval1; expb_p(ipr,i)=yval1(ipr)-expa_p(ipr,i)*xval1
    yval1(ipr)=yval2(ipr)
  enddo
  xval1=xval2
enddo
CONTINUE
DO ipr = 1, npr
  DO i = -39999,0
    expa(i,ipr) = dx * expa(i,ipr); expa_p(ipr,i) = dx * expa_p(ipr,i)
  ENDDO
ENDDO 
end subroutine 