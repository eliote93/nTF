MODULE HighOrderSC
  USE MCDefine
  USE TYPEDEF, ONLY : BenchXs_type
  USE Cntl, ONLY : nTracerCntl
  IMPLICIT NONE
  PUBLIC :: WeightCount, NegWeightCount
  INTEGER :: WeightCount = 0
  INTEGER :: NegWeightCount = 0
  INTEGER :: hscmg
  
  CONTAINS

  SUBROUTINE SmpDirAnIsoGaussian(nst, dir, xs)
    USE RNG
    USE BenchXs
    !USE MCBasics
    USE ioutil, ONLY : terminate
!    USE RootFinding
    IMPLICIT NONE
    TYPE(nstat), INTENT(IN) :: nst
    REAL(8), INTENT(OUT) :: dir(3)
    TYPE(XsecSet), INTENT(IN) ::xs
    REAL(8) :: rn, sum
    REAL(8), POINTER, DIMENSION(:) :: cdf
    REAL(8) :: rt1, rt2, rt3, ww1, ww2, ww3, w1, w2, w3
    INTEGER :: i, ig, jg, ng
    REAL(8) :: mu, ms, alpha
    REAL(8) :: tmpdirx, tmpdiry, tmpdirz

    ig=nst%g
    sum=0._8
    cdf=>xs%SctCDF(:,ig)
    rn=GetRNS(nst%seed)
    ng=hscmg
    DO i=1, ng
      sum=sum+cdf(i)
      IF (rn .lt. sum) EXIT
    END DO
    jg=i

    IF(nTracerCntl%scatod .eq. 1) THEN
      mu=xs%rt1(ig,jg)
    ELSEIF(nTracerCntl%scatod .eq. 2) THEN
      CALL TERMINATE('Input Error in option card. Check SCAT_ORDER')
    ELSEIF(nTracerCntl%scatod .eq. 3) THEN
        rn=GetRNS(nst%seed)
        rt1=xs%rt1(ig,jg)
        rt2=xs%rt2(ig,jg)
        w1= xs%w1(ig,jg)
        w2= xs%w2(ig,jg)
        IF (rn .lt. w1/(w1+w2)) THEN
            mu=rt1
        ELSE
            mu=rt2
        END IF
    ELSEIF(nTracerCntl%scatod .eq. 5) THEN
    END IF

    alpha=2._8*pi*GetRNS(nst%seed)
    tmpdirx=dir(1)
    tmpdiry=dir(2)
    tmpdirz=dir(3)
    dir(1)=mu*tmpdirx+SQRT(1._8-mu**2._8)*(tmpdirx*tmpdirz*cos(alpha)-tmpdiry*sin(alpha))/SQRT(1._8-tmpdirz**2._8)
    dir(2)=mu*tmpdiry+SQRT(1._8-mu**2._8)*(tmpdiry*tmpdirz*cos(alpha)+tmpdirx*sin(alpha))/SQRT(1._8-tmpdirz**2._8)
    dir(3)=mu*tmpdirz-SQRT(1._8-mu**2._8)*SQRT(1._8-tmpdirz**2._8)*cos(alpha)

  END SUBROUTINE

  SUBROUTINE SetAnIsoGaussian(xs, scatod)
    USE BenchXs, ONLY : ngben
    !USE MCBasics
    USE Core_mod, ONLY : GroupInfo
    USE ioutil, ONLY : terminate
    IMPLICIT NONE
    TYPE(XsecSet) ::xs
    REAL(8) :: rn, sum
    REAL(8) :: c0, c1, c2, c3, c4, c5, M0, M1, M2, M3, M4, M5
    REAL(8), POINTER, DIMENSION(:) :: cdf
    REAL(8) :: mu1, a10, N1, sigma1, L1, L2, mu2, a20, a21, L3, mu3, N2, sigma2, a30, a31, a32
    REAL(8) :: rt1, rt2, rt3, ww1, ww2, ww3, w1, w2, w3
    REAL(8) :: Mcheck0, Mcheck1, Mcheck2, Mcheck3
    INTEGER :: i, ig, jg, ng, scatod
    REAL(8) :: mu, ms, alpha
    REAL(8) :: tmpdirx, tmpdiry, tmpdirz
    REAL(8) :: x_old, fx_new, x_new, numerator, denominator

    ng=GroupInfo%ng  
    IF( ng .EQ. 0 )THEN
        ng=ngben
    ENDIF    
    hscmg=ng
    IF( Scatod .EQ. 1)THEN
        ALLOCATE(xs%rt1(ng, ng))
        DO ig=1, ng ! from
            DO jg= 1, ng ! to
                c0=1/2._8
                c1=3*xs%xss1(ig,jg)/(2._8*xs%xss0(ig,jg))
                M0=1._8
                M1=2._8*c1/3._8

                mu1=M1; a10=-M1
                rt1=M1
                mu=rt1
                xs%rt1(ig,jg)=rt1
            ENDDO
        ENDDO
    ELSEIF( ScatOd .EQ. 3)THEN
        ALLOCATE(xs%rt1(ng, ng), xs%rt2(ng, ng), xs%w1(ng, ng), xs%w2(ng, ng))
        DO ig=1, ng ! from
            DO jg= 1, ng ! to
                c0=1/2._8 - 5*xs%xss2(ig,jg)/(4._8*xs%xss0(ig,jg))
                c1=3._8*xs%xss1(ig,jg)/(2._8*xs%xss0(ig,jg)) - 21._8*xs%xss3(ig,jg)/(4._8*xs%xss0(ig,jg))
                c2=15._8*xs%xss2(ig,jg)/(4._8*xs%xss0(ig,jg))
                c3=35._8*xs%xss3(ig,jg)/(4._8*xs%xss0(ig,jg))
                M0=2._8*c2/3._8+2._8*c0; M1=2._8*c3/5._8+2._8*c1/3._8; M2=2._8*c2/5._8+2._8*c0/3._8; M3=2._8*c3/7._8+2._8*c1/5._8;
                
                mu1=M1; a10=-M1; N1=a10*M1+M2; sigma1=N1; L1=M1
                
                L2=a10*M2+M3
                mu2=L2/N1-sigma1*L1/N1
                a20=-a10*mu2-sigma1
                a21=a10-mu2
                rt1=(-a21-SQRT(a21**2-4*a20))/2._8
                rt2=(-a21+SQRT(a21**2-4*a20))/2._8
                ww1=1._8 + (rt1**2._8+2._8*a10*rt1+a10**2._8)/N1; w1=1._8/ww1
                ww2=1._8 + (rt2**2._8+2._8*a10*rt2+a10**2._8)/N1; w2=1._8/ww2
                xs%rt1(ig,jg)=rt1
                xs%rt2(ig,jg)=rt2
                xs%w1(ig,jg)=w1
                xs%w2(ig,jg)=w2
            ENDDO
        ENDDO  
    ELSEIF( ScatOd .NE. 0 )THEN        
        CALL TERMINATE('Input Error in option card. Check SCAT_ORDER')
    ENDIF
    DEALLOCATE(xs%xss1,xs%xss2,xs%xss3)      
        
END SUBROUTINE

  SUBROUTINE SetAnIsoGaussian_MAC(xs, scatod)
    USE BenchXs, ONLY : ngben
    !USE MCBasics
    USE Core_mod, ONLY : GroupInfo
    USE ioutil, ONLY : terminate
    IMPLICIT NONE
    TYPE(BenchXs_type) ::xs
    REAL(8) :: rn, sum
    REAL(8) :: c0, c1, c2, c3, c4, c5, M0, M1, M2, M3, M4, M5
    REAL(8), POINTER, DIMENSION(:) :: cdf
    REAL(8) :: mu1, a10, N1, sigma1, L1, L2, mu2, a20, a21, L3, mu3, N2, sigma2, a30, a31, a32
    REAL(8) :: rt1, rt2, rt3, ww1, ww2, ww3, w1, w2, w3
    REAL(8) :: Mcheck0, Mcheck1, Mcheck2, Mcheck3
    INTEGER :: i, ig, jg, ng, scatod
    REAL(8) :: mu, ms, alpha
    REAL(8) :: tmpdirx, tmpdiry, tmpdirz
    REAL(8) :: x_old, fx_new, x_new, numerator, denominator

    ng=GroupInfo%ng  
    IF( ng .EQ. 0 )THEN
        ng=ngben
    ENDIF    
    hscmg=ng
    IF( Scatod .EQ. 1)THEN
        ALLOCATE(xs%rt1(ng, ng))
        DO ig=1, ng ! from
            DO jg= 1, ng ! to
                c0=1/2._8
                c1=3*xs%xss1(ig,jg)/(2._8*xs%xss0(ig,jg))
                M0=1._8
                M1=2._8*c1/3._8

                mu1=M1; a10=-M1
                rt1=M1
                mu=rt1
                xs%rt1(ig,jg)=rt1
            ENDDO
        ENDDO
        DEALLOCATE(xs%xss1)      
    ELSEIF( ScatOd .EQ. 3)THEN
        ALLOCATE(xs%rt1(ng, ng), xs%rt2(ng, ng), xs%w1(ng, ng), xs%w2(ng, ng))
        DO ig=1, ng ! from
            DO jg= 1, ng ! to
                c0=1/2._8 - 5*xs%xss2(ig,jg)/(4._8*xs%xss0(ig,jg))
                c1=3._8*xs%xss1(ig,jg)/(2._8*xs%xss0(ig,jg)) - 21._8*xs%xss3(ig,jg)/(4._8*xs%xss0(ig,jg))
                c2=15._8*xs%xss2(ig,jg)/(4._8*xs%xss0(ig,jg))
                c3=35._8*xs%xss3(ig,jg)/(4._8*xs%xss0(ig,jg))
                M0=2._8*c2/3._8+2._8*c0; M1=2._8*c3/5._8+2._8*c1/3._8; M2=2._8*c2/5._8+2._8*c0/3._8; M3=2._8*c3/7._8+2._8*c1/5._8;
                
                mu1=M1; a10=-M1; N1=a10*M1+M2; sigma1=N1; L1=M1
                
                L2=a10*M2+M3
                mu2=L2/N1-sigma1*L1/N1
                a20=-a10*mu2-sigma1
                a21=a10-mu2
                rt1=(-a21-SQRT(a21**2-4*a20))/2._8
                rt2=(-a21+SQRT(a21**2-4*a20))/2._8
                ww1=1._8 + (rt1**2._8+2._8*a10*rt1+a10**2._8)/N1; w1=1._8/ww1
                ww2=1._8 + (rt2**2._8+2._8*a10*rt2+a10**2._8)/N1; w2=1._8/ww2
                xs%rt1(ig,jg)=rt1
                xs%rt2(ig,jg)=rt2
                xs%w1(ig,jg)=w1
                xs%w2(ig,jg)=w2
            ENDDO
        ENDDO  
        DEALLOCATE(xs%xss1,xs%xss2,xs%xss3)      
    ELSEIF( ScatOd .NE. 0 )THEN        
        CALL TERMINATE('Input Error in option card. Check SCAT_ORDER')
    ENDIF
        
END SUBROUTINE
!
!  SUBROUTINE SmpDirAnIsoGaussian_OLD(nst, dir, xs)
!    USE RNG
!    USE BenchXs
!    !USE MCBasics
!    USE ioutil, ONLY : terminate
!!    USE RootFinding
!    IMPLICIT NONE
!    TYPE(nstat), INTENT(IN) :: nst
!    REAL(8), INTENT(OUT) :: dir(3)
!    TYPE(XsecSet), INTENT(IN) ::xs
!    REAL(8) :: rn, sum
!    REAL(8) :: c0, c1, c2, c3, c4, c5, M0, M1, M2, M3, M4, M5
!    REAL(8), POINTER, DIMENSION(:) :: cdf
!    REAL(8) :: mu1, a10, N1, sigma1, L1, L2, mu2, a20, a21, L3, mu3, N2, sigma2, a30, a31, a32
!    REAL(8) :: rt1, rt2, rt3, ww1, ww2, ww3, w1, w2, w3
!    REAL(8) :: Mcheck0, Mcheck1, Mcheck2, Mcheck3
!    INTEGER :: i, ig, jg
!    REAL(8) :: mu, ms, alpha
!    REAL(8) :: tmpdirx, tmpdiry, tmpdirz
!    REAL(8) :: x_old, fx_new, x_new, numerator, denominator
!
!    ig=nst%g
!    sum=0._8
!    cdf=>xs%SctCDF(:,ig)
!    rn=GetRNS(nst%seed)
!    DO i=1, ngben
!      sum=sum+cdf(i)
!      IF (rn .lt. sum) EXIT
!    END DO
!    jg=i
!
!    IF(nTracerCntl%scatod .eq. 1) THEN
!      c0=1/2._8
!      c1=3*MacXsBen(xs%idx).xss1(ig,jg)/(2._8*MacXsBen(xs%idx).xss0(ig,jg))
!      M0=1._8
!      M1=2._8*c1/3._8
!
!      mu1=M1; a10=-M1
!      rt1=M1
!      mu=rt1
!    ELSEIF(nTracerCntl%scatod .eq. 2) THEN
!      CALL TERMINATE('Input Error in option card. Check SCAT_ORDER')
!    ELSEIF(nTracerCntl%scatod .eq. 3) THEN
!      c0=1/2._8 - 5*MacXsBen(xs%idx).xss2(ig,jg)/(4._8*MacXsBen(xs%idx).xss0(ig,jg))
!      c1=3._8*MacXsBen(xs%idx).xss1(ig,jg)/(2._8*MacXsBen(xs%idx).xss0(ig,jg)) - 21._8*MacXsBen(xs%idx).xss3(ig,jg)/(4._8*MacXsBen(xs%idx).xss0(ig,jg))
!      c2=15._8*MacXsBen(xs%idx).xss2(ig,jg)/(4._8*MacXsBen(xs%idx).xss0(ig,jg))
!      c3=35._8*MacXsBen(xs%idx).xss3(ig,jg)/(4._8*MacXsBen(xs%idx).xss0(ig,jg))
!      M0=2._8*c2/3._8+2._8*c0; M1=2._8*c3/5._8+2._8*c1/3._8; M2=2._8*c2/5._8+2._8*c0/3._8; M3=2._8*c3/7._8+2._8*c1/5._8;
!      
!      mu1=M1; a10=-M1; N1=a10*M1+M2; sigma1=N1; L1=M1
!      
!      L2=a10*M2+M3
!      mu2=L2/N1-sigma1*L1/N1
!      a20=-a10*mu2-sigma1
!      a21=a10-mu2
!      
!      rt1=(-a21-SQRT(a21**2-4*a20))/2._8
!      rt2=(-a21+SQRT(a21**2-4*a20))/2._8
!      ww1=1._8 + (rt1**2._8+2._8*a10*rt1+a10**2._8)/N1; w1=1._8/ww1
!      ww2=1._8 + (rt2**2._8+2._8*a10*rt2+a10**2._8)/N1; w2=1._8/ww2
!      Mcheck0=w1+w2; Mcheck1=rt1*w1+rt2*w2; Mcheck2=rt1**2._8*w1+rt2**2._8*w2; Mcheck3=rt1**3._8*w1+rt2**3._8*w2;
!      
!      IF (w1 .lt. -1 .or. w1 .gt. 1) THEN
!        NegWeightCount=NegWeightCount+1
!      ELSEIF (w2 .lt. -1 .or. w2 .gt. 1) THEN
!        NegWeightCount=NegWeightCount+1
!      ELSE
!        WeightCount=WeightCount+1
!      END IF
!      
!      rn=GetRNS(nst%seed)
!      IF (rn .lt. w1/(w1+w2)) THEN
!        mu=rt1
!      ELSE
!        mu=rt2
!      END IF
!    ELSEIF(nTracerCntl%scatod .eq. 5) THEN
!!      c0=1/2._8 - 5*MacXsBen(xs%idx).xss2(ig,jg)/(4*MacXsBen(xs%idx).xss0(ig,jg)) + 27._8*MacXsBen(xs%idx).xss4(ig,jg)/(16._8*MacXsBen(xs%idx).xss0(ig,jg))
!!      c1=3*MacXsBen(xs%idx).xss1(ig,jg)/(2*MacXsBen(xs%idx).xss0(ig,jg)) - 21*MacXsBen(xs%idx).xss3(ig,jg)/(4*MacXsBen(xs%idx).xss0(ig,jg)) + 11._8*15._8*MacXsBen(xs%idx).xss5(ig,jg)/(16._8*MacXsBen(xs%idx).xss0(ig,jg))
!!      c2=15*MacXsBen(xs%idx).xss2(ig,jg)/(4*MacXsBen(xs%idx).xss0(ig,jg)) - 270._8*MacXsBen(xs%idx).xss4(ig,jg)/(16._8*MacXsBen(xs%idx).xss0(ig,jg))
!!      c3=35*MacXsBen(xs%idx).xss3(ig,jg)/(4*MacXsBen(xs%idx).xss0(ig,jg)) - 770._8*MacXsBen(xs%idx).xss5(ig,jg)/(16._8*MacXsBen(xs%idx).xss0(ig,jg))
!!      c4=9._8*35._8*MacXsBen(xs%idx).xss4(ig,jg)/(16._8*MacXsBen(xs%idx).xss0(ig,jg))
!!      c5=11._8*63._8*MacXsBen(xs%idx).xss5(ig,jg)/(16._8*MacXsBen(xs%idx).xss0(ig,jg))
!!  
!!      M0=2._8*c4/5._8+2._8*c2/3._8+2._8*c0
!!      M1=2._8*c5/7._8+2._8*c3/5._8+2._8*c1/3._8; M2=2._8*c4/7._8+2._8*c2/5._8+2._8*c0/3._8; 
!!      M3=2._8*c5/9._8+2._8*c3/7._8+2._8*c1/5._8; M4=2._8*c4/9._8+2._8*c2/7._8+2._8*c0/5._8; 
!!      M5=2._8*c5/11._8+2._8*c3/9._8+2._8*c1/7._8; 
!!  
!!      mu1=M1; a10=-M1; N1=a10*M1+M2; sigma1=N1; L1=M1;
!!      L2=a10*M2+M3; mu2=L2/N1-sigma1*L1/N1;
!!      a20=-a10*mu2-sigma1; a21=a10-mu2;
!!  
!!      L3 = a20*M3 + a21*M4 + M5; N2 = a20*M2 + a21*M3 + M4; sigma2 = N2/N1; mu3 = L3/N2 - sigma2*L2/N2;
!!      a30 = -mu3*a20 - sigma2*a10
!!      a31 = a20 - mu3*a21 - sigma2
!!      a32 = a21 - mu3
!!      x_old=1._8; fx_new=1._8;
!!      DO WHILE ( ABS(fx_new) .gt. 1.0E-10 ) 
!!        numerator = PolyValue(x_old,a32,a31,a30)
!!        denominator = Deriv(x_old,a32,a31)
!!        x_new = x_old - numerator/denominator
!!        fx_new = PolyValue(x_new,a32,a31,a30)
!!        x_old = x_new
!!      END DO
!!      rt1 = x_new
!!      
!!      x_old=1._8; fx_new=1._8;
!!      DO WHILE ( ABS(fx_new) .gt. 1.0E-10 ) 
!!        numerator = PolyValue2(x_old,a32,a31,a30,rt1)
!!        denominator = Deriv2(x_old,a32,a31,a30,rt1)
!!        x_new = x_old - numerator/denominator
!!        fx_new = PolyValue2(x_new,a32,a31,a30,rt1)
!!        x_old = x_new
!!      END DO
!!      rt2 = x_new
!!      
!!      x_old=1._8; fx_new=1._8;
!!      DO WHILE ( ABS(fx_new) .gt. 1.0E-10 ) 
!!        numerator = PolyValue3(x_old,a32,a31,a30,rt1,rt2)
!!        denominator = Deriv3(x_old,a32,a31,a30,rt1,rt2)
!!        x_new = x_old - numerator/denominator
!!        fx_new = PolyValue3(x_new,a32,a31,a30,rt1,rt2)
!!        x_old = x_new
!!      END DO
!!      rt3 = x_new
!!      
!!      ww1 = 1._8 + (rt1+a10)**2._8/N1 + (rt1**2._8 + a21*rt1 + a20)**2._8/N2;
!!      w1 = 1._8/ww1
!!      ww2 = 1._8 + (rt2+a10)**2._8/N1 + (rt2**2._8 + a21*rt2 + a20)**2._8/N2;
!!      w2 = 1._8/ww2
!!      ww3 = 1._8 + (rt3+a10)**2._8/N1 + (rt3**2._8 + a21*rt3 + a20)**2._8/N2;
!!      w3 = 1._8/ww3
!!      
!!      CALL RANDOM_NUMBER(rn)
!!      IF( rn .lt. w1 ) THEN
!!        mu = rt1
!!      ELSEIF( rn .gt. w1 .AND. rn .lt. w1+w2 ) THEN
!!        mu = rt2
!!      ELSE
!!        mu = rt3
!!      END IF
!    END IF
!
!    alpha=2._8*pi*GetRNS(nst%seed)
!    tmpdirx=dir(1)
!    tmpdiry=dir(2)
!    tmpdirz=dir(3)
!    dir(1)=mu*tmpdirx+SQRT(1._8-mu**2._8)*(tmpdirx*tmpdirz*cos(alpha)-tmpdiry*sin(alpha))/SQRT(1._8-tmpdirz**2._8)
!    dir(2)=mu*tmpdiry+SQRT(1._8-mu**2._8)*(tmpdiry*tmpdirz*cos(alpha)+tmpdirx*sin(alpha))/SQRT(1._8-tmpdirz**2._8)
!    dir(3)=mu*tmpdirz-SQRT(1._8-mu**2._8)*SQRT(1._8-tmpdirz**2._8)*cos(alpha)
!
!  END SUBROUTINE
!  SUBROUTINE SmpDirAnIsoGaussian_OLDmic(nst, dir, xs)
!    USE RNG
!    USE BenchXs
!    !USE MCBasics
!    USE ioutil, ONLY : terminate
!!    USE RootFinding
!    IMPLICIT NONE
!    TYPE(nstat), INTENT(IN) :: nst
!    REAL(8), INTENT(OUT) :: dir(3)
!    TYPE(XsecSet), INTENT(IN) ::xs
!    REAL(8) :: rn, sum
!    REAL(8) :: c0, c1, c2, c3, c4, c5, M0, M1, M2, M3, M4, M5
!    REAL(8), POINTER, DIMENSION(:) :: cdf
!    REAL(8) :: mu1, a10, N1, sigma1, L1, L2, mu2, a20, a21, L3, mu3, N2, sigma2, a30, a31, a32
!    REAL(8) :: rt1, rt2, rt3, ww1, ww2, ww3, w1, w2, w3
!    REAL(8) :: Mcheck0, Mcheck1, Mcheck2, Mcheck3
!    INTEGER :: i, ig, jg
!    REAL(8) :: mu, ms, alpha
!    REAL(8) :: tmpdirx, tmpdiry, tmpdirz
!    REAL(8) :: x_old, fx_new, x_new, numerator, denominator
!
!    ig=nst%g  !from
!    sum=0._8
!    cdf=>xs%SctCDF(:,ig)
!    rn=GetRNS(nst%seed)
!    DO i=1, ngben
!      sum=sum+cdf(i)
!      IF (rn .lt. sum) EXIT
!    END DO
!    jg=i  ! to
!
!    IF(nTracerCntl%scatod .eq. 1) THEN
!      c0=1/2._8
!      c1=3*xs%xss1(ig,jg)/(2._8*xs%xss0(ig,jg))
!      M0=1._8
!      M1=2._8*c1/3._8
!
!      mu1=M1; a10=-M1
!      rt1=M1
!      mu=rt1
!    ELSEIF(nTracerCntl%scatod .eq. 2) THEN
!      CALL TERMINATE('Input Error in option card. Check SCAT_ORDER')
!    ELSEIF(nTracerCntl%scatod .eq. 3) THEN
!      c0=1/2._8 - 5*xs%xss2(ig,jg)/(4._8*xs%xss0(ig,jg))
!      c1=3._8*xs%xss1(ig,jg)/(2._8*xs%xss0(ig,jg)) - 21._8*xs%xss3(ig,jg)/(4._8*xs%xss0(ig,jg))
!      c2=15._8*xs%xss2(ig,jg)/(4._8*xs%xss0(ig,jg))
!      c3=35._8*xs%xss3(ig,jg)/(4._8*xs%xss0(ig,jg))
!      M0=2._8*c2/3._8+2._8*c0; M1=2._8*c3/5._8+2._8*c1/3._8; M2=2._8*c2/5._8+2._8*c0/3._8; M3=2._8*c3/7._8+2._8*c1/5._8;
!      
!      mu1=M1; a10=-M1; N1=a10*M1+M2; sigma1=N1; L1=M1
!      
!      L2=a10*M2+M3
!      mu2=L2/N1-sigma1*L1/N1
!      a20=-a10*mu2-sigma1
!      a21=a10-mu2
!      
!      rt1=(-a21-SQRT(a21**2-4*a20))/2._8
!      rt2=(-a21+SQRT(a21**2-4*a20))/2._8
!      ww1=1._8 + (rt1**2._8+2._8*a10*rt1+a10**2._8)/N1; w1=1._8/ww1
!      ww2=1._8 + (rt2**2._8+2._8*a10*rt2+a10**2._8)/N1; w2=1._8/ww2
!      Mcheck0=w1+w2; Mcheck1=rt1*w1+rt2*w2; Mcheck2=rt1**2._8*w1+rt2**2._8*w2; Mcheck3=rt1**3._8*w1+rt2**3._8*w2;
!      
!      IF (w1 .lt. -1 .or. w1 .gt. 1) THEN
!        NegWeightCount=NegWeightCount+1
!      ELSEIF (w2 .lt. -1 .or. w2 .gt. 1) THEN
!        NegWeightCount=NegWeightCount+1
!      ELSE
!        WeightCount=WeightCount+1
!      END IF
!      
!      rn=GetRNS(nst%seed)
!      IF (rn .lt. w1/(w1+w2)) THEN
!        mu=rt1
!      ELSE
!        mu=rt2
!      END IF
!    ELSEIF(nTracerCntl%scatod .eq. 5) THEN
!!      c0=1/2._8 - 5*MacXsBen(xs%idx).xss2(ig,jg)/(4*MacXsBen(xs%idx).xss0(ig,jg)) + 27._8*MacXsBen(xs%idx).xss4(ig,jg)/(16._8*MacXsBen(xs%idx).xss0(ig,jg))
!!      c1=3*MacXsBen(xs%idx).xss1(ig,jg)/(2*MacXsBen(xs%idx).xss0(ig,jg)) - 21*MacXsBen(xs%idx).xss3(ig,jg)/(4*MacXsBen(xs%idx).xss0(ig,jg)) + 11._8*15._8*MacXsBen(xs%idx).xss5(ig,jg)/(16._8*MacXsBen(xs%idx).xss0(ig,jg))
!!      c2=15*MacXsBen(xs%idx).xss2(ig,jg)/(4*MacXsBen(xs%idx).xss0(ig,jg)) - 270._8*MacXsBen(xs%idx).xss4(ig,jg)/(16._8*MacXsBen(xs%idx).xss0(ig,jg))
!!      c3=35*MacXsBen(xs%idx).xss3(ig,jg)/(4*MacXsBen(xs%idx).xss0(ig,jg)) - 770._8*MacXsBen(xs%idx).xss5(ig,jg)/(16._8*MacXsBen(xs%idx).xss0(ig,jg))
!!      c4=9._8*35._8*MacXsBen(xs%idx).xss4(ig,jg)/(16._8*MacXsBen(xs%idx).xss0(ig,jg))
!!      c5=11._8*63._8*MacXsBen(xs%idx).xss5(ig,jg)/(16._8*MacXsBen(xs%idx).xss0(ig,jg))
!!  
!!      M0=2._8*c4/5._8+2._8*c2/3._8+2._8*c0
!!      M1=2._8*c5/7._8+2._8*c3/5._8+2._8*c1/3._8; M2=2._8*c4/7._8+2._8*c2/5._8+2._8*c0/3._8; 
!!      M3=2._8*c5/9._8+2._8*c3/7._8+2._8*c1/5._8; M4=2._8*c4/9._8+2._8*c2/7._8+2._8*c0/5._8; 
!!      M5=2._8*c5/11._8+2._8*c3/9._8+2._8*c1/7._8; 
!!  
!!      mu1=M1; a10=-M1; N1=a10*M1+M2; sigma1=N1; L1=M1;
!!      L2=a10*M2+M3; mu2=L2/N1-sigma1*L1/N1;
!!      a20=-a10*mu2-sigma1; a21=a10-mu2;
!!  
!!      L3 = a20*M3 + a21*M4 + M5; N2 = a20*M2 + a21*M3 + M4; sigma2 = N2/N1; mu3 = L3/N2 - sigma2*L2/N2;
!!      a30 = -mu3*a20 - sigma2*a10
!!      a31 = a20 - mu3*a21 - sigma2
!!      a32 = a21 - mu3
!!      x_old=1._8; fx_new=1._8;
!!      DO WHILE ( ABS(fx_new) .gt. 1.0E-10 ) 
!!        numerator = PolyValue(x_old,a32,a31,a30)
!!        denominator = Deriv(x_old,a32,a31)
!!        x_new = x_old - numerator/denominator
!!        fx_new = PolyValue(x_new,a32,a31,a30)
!!        x_old = x_new
!!      END DO
!!      rt1 = x_new
!!      
!!      x_old=1._8; fx_new=1._8;
!!      DO WHILE ( ABS(fx_new) .gt. 1.0E-10 ) 
!!        numerator = PolyValue2(x_old,a32,a31,a30,rt1)
!!        denominator = Deriv2(x_old,a32,a31,a30,rt1)
!!        x_new = x_old - numerator/denominator
!!        fx_new = PolyValue2(x_new,a32,a31,a30,rt1)
!!        x_old = x_new
!!      END DO
!!      rt2 = x_new
!!      
!!      x_old=1._8; fx_new=1._8;
!!      DO WHILE ( ABS(fx_new) .gt. 1.0E-10 ) 
!!        numerator = PolyValue3(x_old,a32,a31,a30,rt1,rt2)
!!        denominator = Deriv3(x_old,a32,a31,a30,rt1,rt2)
!!        x_new = x_old - numerator/denominator
!!        fx_new = PolyValue3(x_new,a32,a31,a30,rt1,rt2)
!!        x_old = x_new
!!      END DO
!!      rt3 = x_new
!!      
!!      ww1 = 1._8 + (rt1+a10)**2._8/N1 + (rt1**2._8 + a21*rt1 + a20)**2._8/N2;
!!      w1 = 1._8/ww1
!!      ww2 = 1._8 + (rt2+a10)**2._8/N1 + (rt2**2._8 + a21*rt2 + a20)**2._8/N2;
!!      w2 = 1._8/ww2
!!      ww3 = 1._8 + (rt3+a10)**2._8/N1 + (rt3**2._8 + a21*rt3 + a20)**2._8/N2;
!!      w3 = 1._8/ww3
!!      
!!      CALL RANDOM_NUMBER(rn)
!!      IF( rn .lt. w1 ) THEN
!!        mu = rt1
!!      ELSEIF( rn .gt. w1 .AND. rn .lt. w1+w2 ) THEN
!!        mu = rt2
!!      ELSE
!!        mu = rt3
!!      END IF
!    END IF
!
!    alpha=2._8*pi*GetRNS(nst%seed)
!    tmpdirx=dir(1)
!    tmpdiry=dir(2)
!    tmpdirz=dir(3)
!    dir(1)=mu*tmpdirx+SQRT(1._8-mu**2._8)*(tmpdirx*tmpdirz*cos(alpha)-tmpdiry*sin(alpha))/SQRT(1._8-tmpdirz**2._8)
!    dir(2)=mu*tmpdiry+SQRT(1._8-mu**2._8)*(tmpdiry*tmpdirz*cos(alpha)+tmpdirx*sin(alpha))/SQRT(1._8-tmpdirz**2._8)
!    dir(3)=mu*tmpdirz-SQRT(1._8-mu**2._8)*SQRT(1._8-tmpdirz**2._8)*cos(alpha)
!
!  END SUBROUTINE
END MODULE