#define exactki3
module CP_mod
    implicit none
    
    integer,parameter :: ngauss=5
    REAL,parameter :: pi=3.141592653589793_8, hpi=1.570796326794897_8 , ki0=0.7853981633974483_8
    REAL :: Quad(ngauss),Wgt(ngauss)
    
    contains
    
    subroutine runCP(phi,DelRk,vol,Qsurfvol,invsurf4,X,sigt,q,nr)
    ! white boundary condition without self-scattering XS
        implicit none
        integer,intent(in) :: nr
        real,intent(in) :: DelRk(nr),sigt(nr),q(nr),vol(nr),Qsurfvol(nr),invsurf4,X(nr,nr,ngauss)
        real,intent(out) :: phi(nr)        
        
        INTEGER :: k,k2
        REAL :: Pij(nr,nr),gam(nr),b(nr),A(nr,nr),Xik(nr,nr),Y(nr),xk(nr)
        REAL :: sigtvol(nr),invsigt(nr),sigrvol(nr),invsigtvol(nr)
        
        !vol(1)=pi*rad(1)*rad(1)
        !sumvol=vol(1)
        !DO i=2,nr
        !    vol(i)=pi*rad(i)*rad(i)-sumvol
        !    sumvol=sumvol+vol(i)
        !ENDDO 
        !quartersurf=hpi*rad(nr)
        
        invsigt=1._8/sigt
        sigtvol=vol*sigt
        invsigtvol=1._8/sigtvol
        sigrvol=sigtvol   ! zero scattering XS
        !DO i=1,nr
        !    Qsurfvol(i)=vol(i)*quartersurf
        !ENDDO
        
        call CalcPijGam(Pij,gam,X,DelRk,invsurf4,sigt,sigtvol,nr)
        !CALL ConstCPA(A,Pij,sigtvol,nr)  
        !CALL LUfac(A,nr)
        DO k=1,nr
            !CALL ConstCPB(Pij,b,invsigt,k,nr)   ! Construct # of collisions in each region due to unit source in region, k.
            !CALL SolveSys(A,b,Xik(1:nr,k),nr)   ! Obtain Xik(1:nr,k), all fluxes of each region, due to unit source in region, k.
            !! getXiK
            DO k2=1,nr
                Xik(k2,k)=Pij(k2,k)*invsigtvol(k2)*invsigt(k)
            ENDDO
        ENDDO
        !CALL SolveSys(A,gam,Y,nr)               ! Obtain Y, all fluxes of each region, due to one neutron coming from the outer surface.
        !! get Y
        DO k=1,nr
            Y(k)=gam(k)*invsigtvol(k)
        ENDDO
        CALL ApplyAlbedo(Xik,xk,Y,Qsurfvol,sigrvol,nr)  
        CALL GetPhi(phi,Xik,q,nr)
    
    end subroutine
    
    SUBROUTINE CalcPijGam(Pij,gam,X,DelRk,invsurf4,sigt,sigtvol,nr)
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: nr
        REAL,INTENT(IN) :: DelRk(nr),invsurf4,sigt(nr),sigtvol(nr),X(nr,nr,ngauss)
        REAL,INTENT(OUT) :: Pij(nr,nr),gam(nr)
        INTEGER :: m,i,k,n,j
        REAL :: sigt2(nr),Y,Y22,SUMX,a,b,c,d,e,f,sumpij
        REAL :: Sij(0:nr,0:nr),Sijm(1:nr,1:nr),TAUN(1:nr,1:nr),TAUP(1:nr,1:nr)
        LOGICAL,SAVE :: lfirst=.true.
        
        if (lfirst) then
            CALL GaussJacobiQuad(ngauss,Quad,Wgt)
            lfirst=.false.
        endif
            
        !Rk=r    
        sigt2 = 2._8*sigt        
        !DelRk(1)=Rk(1)
        !Rk22(1)=Rk(1)*Rk(1)
        !DO i=2,nr
        !    DelRk(i)=Rk(i)-Rk(i-1)
        !    Rk22(i)=Rk(i)*Rk(i)
        !ENDDO  
        !
        !invsurf4=1._8/(hpi*Rk(nr))
        
        !DO m=1,ngauss
        !    DO k=1,nr
        !        Y=Rk(k)-DelRk(k)*Quad(m); Y22=Y*Y
        !        X(k,k,m)=dsqrt(RK22(k)-Y22)
        !        SUMX=X(k,k,m)
        !        DO n=k+1,nr
        !            X(n,k,m)=dsqrt(RK22(n)-Y22)-SUMX
        !            SUMX=SUMX+X(n,k,m)
        !        ENDDO
        !    ENDDO
        !ENDDO   
        !Do i=1,nr
        !    DelRk(i)=DelRk(i)*4._8
        !ENDDO
        
        Pij=0._8; gam=0._8
        Sij=0._8
        DO m=1,ngauss
            Sijm=0._8
            DO k=1,nr
                TAUN=0._8; TAUP=0._8;
                a=sigt2(k)*X(k,k,m); TAUP(k,k)=a
                Sijm(k,k)=Sijm(k,k)+(ki_3(a)-ki0)*DelRk(k)
                DO n=k+1,nr
                    b=sigt(n)*X(n,k,m);  TAUN(n-1,n)=b
                    c=b+TAUP(n-1,n-1);   TAUP(n-1,n)=c
                    Sijm(n-1,n)=Sijm(n-1,n)+(ki_3(c)-ki_3(b))*DelRk(k)
                    d=b+c;    TAUP(n,n)=d
                    Sijm(n,n)=Sijm(n,n)+(ki_3(d)-ki0)*DelRk(k)
                    DO j=k,n-2
                        e=TAUN(j,n-1)+b;  TAUN(j,n)=e
                        f=e+TAUP(j,j);    TAUP(j,n)=f
                        Sijm(j,n)=Sijm(j,n)+(Ki_3(f)-ki_3(e))*DelRk(k)
                    ENDDO
                ENDDO           
            ENDDO
            DO n=1,nr
                DO j=1,n
                    Sij(j,n)=Sij(j,n)+Sijm(j,n)*WGt(m)
                ENDDO
            ENDDO
        ENDDO
        DO n=1,nr
            DO j=1,n-1
                Sij(n,j)=Sij(j,n)        
            ENDDO
        ENDDO
        DO n=1,nr
            Pij(n,n)=sigtvol(n)
            gam(n)=Pij(n,n)
        ENDDO
        DO n=1,nr
            DO j=1,n-1
                Pij(j,n)=Pij(j,n)+(Sij(j,n)+Sij(j-1,n-1)-Sij(j,n-1)-Sij(j-1,n))
                Pij(n,j)=Pij(j,n)
            ENDDO
            Pij(n,n)=Pij(n,n)+(Sij(n,n)+Sij(n-1,n-1)-Sij(n,n-1)-Sij(n-1,n))
        ENDDO        
        DO n=1,nr
            sumpij=sum(Pij(:,n)) ! originally sum(Pij(:,n)), but using reciprocity relation for optimizing indexing for fast programming
            gam(n)=gam(n)-sumpij
            gam(n)=gam(n)*invsurf4
        ENDDO   
    END SUBROUTINE  
    
    subroutine runCP_(phi,rad,sigt,q,nr)
    ! white boundary condition without self-scattering XS
        implicit none
        integer,intent(in) :: nr
        real,intent(in) :: rad(nr),sigt(nr),q(nr)
        real,intent(out) :: phi(nr)        
        
        INTEGER :: k,i,k2
        REAL :: sumvol,vol(nr)
        REAL :: Pij(nr,nr),gam(nr),b(nr),A(nr,nr),Xik(nr,nr),Y(nr),xk(nr)
        REAL :: sigtvol(nr),invsigt(nr),quartersurf,quartersurfvol(nr),sigrvol(nr),invsigtvol(nr)
        
        vol(1)=pi*rad(1)*rad(1)
        sumvol=vol(1)
        DO i=2,nr
            vol(i)=pi*rad(i)*rad(i)-sumvol
            sumvol=sumvol+vol(i)
        ENDDO 
        quartersurf=hpi*rad(nr)
        
        invsigt=1._8/sigt
        sigtvol=vol*sigt
        invsigtvol=1._8/sigtvol
        sigrvol=sigtvol   ! zero scattering XS
        DO i=1,nr
            quartersurfvol(i)=vol(i)*quartersurf
        ENDDO
        
        call CalcPijGam_(Pij,gam,rad,sigt,sigtvol,nr)
        !CALL ConstCPA(A,Pij,sigtvol,nr)  
        !CALL LUfac(A,nr) 
        DO k=1,nr
            !CALL ConstCPB(Pij,b,invsigt,k,nr)   ! Construct # of collisions in each region due to unit source in region, k.
            !CALL SolveSys(A,b,Xik(1:nr,k),nr)   ! Obtain Xik(1:nr,k), all fluxes of each region, due to unit source in region, k.
            !! getXiK
            DO k2=1,nr
                Xik(k2,k)=Pij(k2,k)*invsigtvol(k2)*invsigt(k)
            ENDDO
        ENDDO
        !CALL SolveSys(A,gam,Y,nr)               ! Obtain Y, all fluxes of each region, due to one neutron coming from the outer surface.
        !! get Y
        DO k=1,nr
            Y(k)=gam(k)*invsigtvol(k)
        ENDDO
        CALL ApplyAlbedo(Xik,xk,Y,quartersurfvol,sigrvol,nr)  
        CALL GetPhi(phi,Xik,q,nr)
    
    end subroutine
    
    SUBROUTINE CalcPijGam_(Pij,gam,r,sigt,sigtvol,nr)
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: nr
        REAL,INTENT(IN) :: r(nr),sigt(nr),sigtvol(nr)
        REAL,INTENT(OUT) :: Pij(nr,nr),gam(nr)
        INTEGER :: m,i,k,n,j
        REAL :: Rk(nr),DelRk(nr),Rk22(nr),fourinvsurf
        REAL :: sigt2(nr),Y,Y22,SUMX,X(nr,nr,ngauss),a,b,c,d,e,f,sumpij
        REAL :: Sij(0:nr,0:nr),Sijm(1:nr,1:nr),TAUN(1:nr,1:nr),TAUP(1:nr,1:nr)
        LOGICAL,SAVE :: lfirst=.true.
        
        if (lfirst) then
            CALL GaussJacobiQuad(ngauss,Quad,Wgt)
            lfirst=.false.
        endif
            
        Rk=r    
        sigt2 = 2._8*sigt        
        DelRk(1)=Rk(1)
        Rk22(1)=Rk(1)*Rk(1)
        DO i=2,nr
            DelRk(i)=Rk(i)-Rk(i-1)
            Rk22(i)=Rk(i)*Rk(i)
        ENDDO  
    
        fourinvsurf=1._8/(hpi*Rk(nr))
        
        DO m=1,ngauss
            DO k=1,nr
                Y=Rk(k)-DelRk(k)*Quad(m); Y22=Y*Y
                X(k,k,m)=dsqrt(RK22(k)-Y22)
                SUMX=X(k,k,m)
                DO n=k+1,nr
                    X(n,k,m)=dsqrt(RK22(n)-Y22)-SUMX
                    SUMX=SUMX+X(n,k,m)
                ENDDO
            ENDDO
        ENDDO   
        Do i=1,nr
            DelRk(i)=DelRk(i)*4._8
        ENDDO
        
        Pij=0._8; gam=0._8
        Sij=0._8
        DO m=1,ngauss
            Sijm=0._8
            DO k=1,nr
                TAUN=0._8; TAUP=0._8;
                a=sigt2(k)*X(k,k,m); TAUP(k,k)=a
                Sijm(k,k)=Sijm(k,k)+(ki_3(a)-ki0)*DelRk(k)
                DO n=k+1,nr
                    b=sigt(n)*X(n,k,m);  TAUN(n-1,n)=b
                    c=b+TAUP(n-1,n-1);   TAUP(n-1,n)=c
                    Sijm(n-1,n)=Sijm(n-1,n)+(ki_3(c)-ki_3(b))*DelRk(k)
                    d=b+c;    TAUP(n,n)=d
                    Sijm(n,n)=Sijm(n,n)+(ki_3(d)-ki0)*DelRk(k)
                    DO j=k,n-2
                        e=TAUN(j,n-1)+b;  TAUN(j,n)=e
                        f=e+TAUP(j,j);    TAUP(j,n)=f
                        Sijm(j,n)=Sijm(j,n)+(Ki_3(f)-ki_3(e))*DelRk(k)
                    ENDDO
                ENDDO           
            ENDDO
            DO n=1,nr
                DO j=1,n
                    Sij(j,n)=Sij(j,n)+Sijm(j,n)*WGt(m)
                ENDDO
            ENDDO
        ENDDO
        DO n=1,nr
            DO j=1,n-1
                Sij(n,j)=Sij(j,n)        
            ENDDO
        ENDDO
        DO n=1,nr
            Pij(n,n)=sigtvol(n)
            gam(n)=Pij(n,n)
        ENDDO
        DO n=1,nr
            DO j=1,n-1
                Pij(j,n)=Pij(j,n)+(Sij(j,n)+Sij(j-1,n-1)-Sij(j,n-1)-Sij(j-1,n))
                Pij(n,j)=Pij(j,n)
            ENDDO
            Pij(n,n)=Pij(n,n)+(Sij(n,n)+Sij(n-1,n-1)-Sij(n,n-1)-Sij(n-1,n))
        ENDDO        
        DO n=1,nr
            sumpij=sum(Pij(:,n)) ! originally sum(Pij(:,n)), but using reciprocity relation for optimizing indexing for fast programming
            gam(n)=gam(n)-sumpij
            gam(n)=gam(n)*fourinvsurf
        ENDDO   
    END SUBROUTINE  
    SUBROUTINE GaussJacobiQuad(N,P22,W)
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: N
        REAL,INTENT(OUT) :: P22(N),W(N)
        INTEGER :: i
        REAL :: GJC(6,6)
        DATA GJC/                                                         &
            .0         ,.6666666667,.3550510257,.2123405382,.1397598643,.0985350858, & 
            .5000000000,.0         ,.8449489743,.5905331356,.4164095676,.3045357266, &    
            .1819586183,.3180413817,.0         ,.9114120405,.7231569864,.5620251898, &  
            .0698269799,.2292411064,.2009319137,.0         ,.9428958039,.8019865821, &  
            .0311809710,.1298475476,.2034645680,.1355069134,.0         ,.9601901429, &  
            .0157479145,.0739088701,.1463869871,.1671746381,.0967815902,.0       /  
        DO i=1,N
            P22(i)=GJC(N+1,i)*GJC(N+1,i)
            W(i)=GJC(i,N+1)
        ENDDO
    END SUBROUTINE     
#ifdef exactki3
    FUNCTION KI_3(x)
!     computes Bickley-Naylor functions
      IMPLICIT NONE
      ! arguments
      DOUBLE PRECISION :: KI_3
      DOUBLE PRECISION, INTENT(IN) :: x 
!**
      DOUBLE PRECISION, DIMENSION(0:7), PARAMETER :: ki0=     &
           (/1.D+20,1.570796326794897D0,1.D0,0.785398163397448D0, &
           2.D0/3.D0,0.589048622548086D0,8.D0/15.D0,0.490873852123405D0/)                                    
      DOUBLE PRECISION, PARAMETER :: gamma=0.57721566490153D0
      DOUBLE PRECISION, DIMENSION(21), PARAMETER :: a= &
               (/0.D0,0.D0,1.3845825179D-02,1.9748007111d-01,8.3482645579D-01, & 
                1.2570304451D0,5.8904861765d-01,0.D0,0.D0,0.0050175543121D0, &    
                0.086457801977D0,0.4541310603D0,0.87229048537D0,0.53333333163D0,0.D0, &    
                9.4163370495D-04,0.023256193951D0,0.1917908614D0,0.6605888802D0, &     
                0.96563978851D0,0.49087385208D0/)
      DOUBLE PRECISION, DIMENSION(21), PARAMETER :: b= &
               (/0.D0,1.104745766d-02,0.19208156163D0,1.1126047386D0,2.6888166891D0, &     
                2.765768548D0,1.D0,0.D0,4.0034500346D-03,0.083493681282D0, &    
                0.58951057779D0,1.7577448589D0,2.2400104946D0,1.D0,7.513155597D-04 , &     
                0.021654858879D0,0.2238978073D0,1.0572470136D0,2.4117190402D0, &     
                2.5536829558D0,1.D0/)                                                                                                                                                  
      DOUBLE PRECISION, DIMENSION(3) :: za,zb 
      DOUBLE PRECISION, DIMENSION(0:7) :: ki
      DOUBLE PRECISION               :: z,c,x1,es
      INTEGER                        :: i,i1,j,j2                                                      
!----------------------------------------------------------------------- 
      IF(x <= 1.0D0)THEN                                                    
         IF(x <= 0.D0)THEN
            ki=ki0
         ELSE
            z=0.5D0*x; c=-(gamma+DLOG(z)); z=z*z 
            ki(3)=ki0(2)-x*(ki0(1)/2.D0-x*((c+3.0D0/2.D0+1.0D0/3.D0)/6.D0+z* &
                  ((c+4.D0/3.D0+1.D0/4.D0+1.D0/5.D0)/60.D0+(z/4.D0)* &
                  ((c+1.7D0+1.D0/6.D0+1.D0/7.D0)/210.D0+(z/9.D0)*((c+11.D0/6.D0+1.D0/7.D0+1.D0/8.D0+1.D0/ &
                   9.D0)/504.D0+(z/15840.D0)*(c+25.D0/12.D0+1.D0/9.D0+1.D0/10.D0+1.D0/11.D0))))))
            ki(3)=ki0(3)-x*ki(3)
         ENDIF
      ELSE
        x1=1.D0/x; es=DEXP(-x)*DSQRT(1.D0+x); i1=1                                                                   
        DO i=1,3                                                                                                                          
            za(i)=a(i1); zb(i)=b(i1)                                                            
            DO j=1,6                                                             
               j2=i1+j; za(i)=za(i)*x+a(j2); zb(i)=zb(i)*x+b(j2)
            ENDDO
            i1=i1+7     
        ENDDO
        ki(5)=za(1)/zb(1)*es                                                  
        ki(6)=za(2)/zb(2)*es                                           
        ki(7)=za(3)/zb(3)*es                                                  
        ki(4)=ki(6)+x1*(ki(7)+5.D0*(ki(7)-ki(5)))                          
        ki(3)=ki(5)+x1*(ki(6)+4.D0*(ki(6)-ki(4)))     
      ENDIF
      KI_3 = ki(3)                                                     
    END FUNCTION KI_3
#else
    FUNCTION KI_3(x)
      REAL :: XX,X,KI_3
      INTEGER :: I
      X   =ABS(XX)                                                            
      IF(X.GT.0.99999)  GO TO 16                                              
      I  =INT(20.0*X+1.00001)                                                 
      GO TO (1,2,3,4,5,6,7,8,9,10,11,11,12,12,13,13,14,14,15,15),I            
!                          ** RANGE 0.00-0.05 **                              
    1 KI_3=(.7266088*X-.9990226)*X+.7853961                                    
      RETURN                                                                  
!                          ** RANGE 0.05-0.10 **                              
    2 KI_3=(.6466375*X-.9912340)*X+.7852024                                    
      RETURN                                                                  
!                          ** RANGE 0.10-0.15 **                              
    3 KI_3=(.5856605*X-.9791293)*X+.7845986                                    
      RETURN                                                                  
!                          ** RANGE 0.15-0.20 **                              
    4 KI_3=(.5346648*X-.9638914)*X+.7834577                                    
      RETURN                                                                  
!                          ** RANGE 0.20-0.25 **                              
    5 KI_3=(.4907827*X-.9463843)*X+.7817094                                    
      RETURN                                                                  
!                          ** RANGE 0.25-0.30 **                              
    6 KI_3=(.4521752*X-.9271152)*X+.7793031                                    
      RETURN                                                                  
!                          ** RANGE 0.30-0.35 **                              
    7 KI_3=(.4177388*X-.9064822)*X+.7762107                                    
      RETURN                                                                  
!                          ** RANGE 0.35-0.40 **                              
    8 KI_3=(.3869945*X-.8849865)*X+.7724519                                    
      RETURN                                                                  
!                          ** RANGE 0.40-0.45 **                              
    9 KI_3=(.3590753*X-.8626685)*X+.7679903                                    
      RETURN                                                                  
!                          ** RANGE 0.45-0.50 **                              
   10 KI_3=(.3338676*X-.8400133)*X+.7628988                                    
      RETURN                                                                  
!                          ** RANGE 0.50-0.60 **                              
   11 KI_3=(.2998569*X-.8054172)*X+.7540982                                    
      RETURN                                                                  
!                          ** RANGE 0.60-0.70 **                              
   12 KI_3=(.2609154*X-.7587821)*X+.7401279                                    
      RETURN                                                                  
!                          ** RANGE 0.70-0.80 **                              
   13 KI_3=(.2278226*X-.7125290)*X+.7239594                                    
      RETURN                                                                  
!                          ** RANGE 0.80-0.90 **                              
   14 KI_3=(.1994999*X-.6672761)*X+.7058777                                    
      RETURN                                                                  
!                          ** RANGE 0.90-1.00 **                              
   15 KI_3=(.1751248*X-.6234536)*X+.6861762                                    
      RETURN                                                                  
                                                                             
   16 IF(X.GT.9.0) GO TO 160                                                  
      I   =INT(2.5*(X-0.99998)) +1                                            
      GO TO (17,18,19,20,21,22,23,24,25,26,27,27,28,28,29,29,30,30,31,31,31),I                              
  160 KI_3=0.0                                                                 
     RETURN                                                                  
!                           ** RANGE 1.0-1.4 **                               
   17 KI_3=((-.05337485*X+.3203223)*X-.7538355)*X+.7247294                     
      RETURN                                                                  
!                           ** RANGE 1.4-1.8 **                               
   18 KI_3=((-.03146833*X+.2295280)*X-.6279752)*X+.6663720                     
      RETURN                                                                  
!                           ** RANGE 1.8-2.2 **                               
   19 KI_3=((-.01906198*X+.1631667)*X-.5094124)*X+.5956163                     
      RETURN                                                                  
!                           ** RANGE 2.2-2.6 **                               
   20 KI_3=((-.01174752*X+.1152418)*X-.4046007)*X+.5191031                     
      RETURN                                                                  
!                           ** RANGE 2.6-3.0 **                               
   21 KI_3=((-.007328415*X+.08097913)*X-.3159648)*X+.4425954                   
      RETURN                                                                  
!                           ** RANGE 3.0-3.4 **                               
   22 KI_3=((-.004617254*X+.05669960)*X-.2434341)*X+.3703178                   
      RETURN                                                                  
!                           ** RANGE 3.4-3.8 **                               
   23 KI_3=(.007923547*X-.07158569)*X+.1684022                                 
      RETURN                                                                  
!                           ** RANGE 3.8-4.2 **                               
   24 KI_3=(.005095111*X-.05016344)*X+.1278307                                 
      RETURN                                                                  
!                           ** RANGE 4.2-4.6 **                               
   25 KI_3=(.003286040*X-.03501524)*X+.09611422                                
      RETURN                                                                  
!                           ** RANGE 4.6-5.0 **                               
   26 KI_3=(.002126242*X-.02437465)*X+.07170491                                
      RETURN                                                                  
!                           ** RANGE 5.0-5.8 **                               
   27 KI_3=(.001123687*X-.01425519)*X+.04616317                                
      RETURN                                                                  
!                           ** RANGE 5.8-6.6 **                               
   28 KI_3=(4.762937E-4*X-6.810124E-3)*X+.02475115                             
      RETURN                                                                  
!                           ** RANGE 6.6-7.4 **                               
   29 KI_3=(2.031843E-4*X-3.232035E-3)*X+.01302864                             
      RETURN                                                                  
!                           ** RANGE 7.4-8.2 **                               
   30 KI_3=(8.701440E-5*X-1.524126E-3)*X+6.749972E-3                           
      RETURN                                                                  
!                           ** RANGE 8.2-9.0 **                               
   31 KI_3=(3.742673E-5*X-7.157367E-4)*X+3.454768E-3                           
                                                                             
      RETURN                                         
    END FUNCTION KI_3
#endif
    SUBROUTINE ConstCPA(A,Pij,sigtvol,nr)
        IMPLICIT NONE    
        INTEGER,INTENT(IN) :: nr
        REAL,INTENT(IN) :: Pij(nr,nr),sigtvol(nr)
        REAL,INTENT(OUT) :: A(nr,nr)
        INTEGER :: i,j
        REAL :: c
        DO i=1,nr
          DO j=1,nr  
            A(j,i)=0._8
          ENDDO
          A(i,i)=sigtvol(i)
        ENDDO
    END SUBROUTINE 
    
    SUBROUTINE ConstCPB(Pij,b,invsigt,k,nr)
        IMPLICIT NONE    
        INTEGER,INTENT(IN) :: k,nr
        REAL,INTENT(IN) :: Pij(nr,nr),invsigt(nr)
        REAL,INTENT(OUT) :: b(nr)
        INTEGER :: i    
        REAL :: invsigtk
        invsigtk=invsigt(k)  
        DO i=1,nr
            b(i)=Pij(k,i)*invsigtk
        ENDDO
    END SUBROUTINE
    
    SUBROUTINE SolveSys(A,b,Sol,N)    
        IMPLICIT NONE    
        INTEGER,INTENT(IN) :: N
        REAL,INTENT(IN) :: A(N,N),b(N)
        REAL,INTENT(OUT) :: Sol(N)    
        REAL :: y(N),sum
        INTEGER :: i,j
    
        Sol=0._8
        y(1)=b(1)
        DO i=2,N
            sum=0._8
            DO j=1,i-1
                sum=sum+A(i,j)*y(j)
            ENDDO
            y(i)=b(i)-sum
        ENDDO
        Sol(N)=y(N)/A(N,N)
        DO i=N-1,1,-1
            sum=0._8
            DO j=i+1,N
                sum=sum+A(i,j)*Sol(j)
            ENDDO
            Sol(i)=(y(i)-sum)/A(i,i)
        ENDDO    
    END SUBROUTINE
    
    SUBROUTINE LUfac(A,N)
        IMPLICIT NONE        
        INTEGER,INTENT(IN) :: N
        REAL,INTENT(INOUT) :: A(N,N)    
        INTEGER :: k,i,j
        REAL :: fmult,rpivot
    
        DO k=1,N-1
            rpivot=1._8/A(k,k)
            DO i=k+1,N
                fmult=rpivot*A(i,k)
                A(i,k)=fmult
                DO j=k+1,N
                    A(i,j)=A(i,j)-fmult*A(k,j)
                ENDDO
            ENDDO
        ENDDO    
    END SUBROUTINE
    
    SUBROUTINE ApplyAlbedo(Xik,xk,Y,quartersurfvol,sigrvol,nr)
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: nr
        REAL,INTENT(IN) :: quartersurfvol(nr),sigrvol(nr)
        REAL,INTENT(INOUT) :: Xik(nr,nr),xk(nr),Y(nr)
        REAL :: denom,totgam,Lgam(nr),xkk
        INTEGER :: i,k
    
        totgam = 0._8
        DO k=1,nr
            xk(k)=quartersurfvol(k)*Y(k)        ! Obtain xk, # of neutrons reaching to the outer surface from region, k, without collision.
            Lgam(k)=sigrvol(k)*Y(k)             ! Obtain Lgam, absorption partial blackness.
            totgam = totgam + Lgam(k)           ! Obtain totgam, absorption partial blackness.
        ENDDO
        
        denom=1._8/totgam
        DO k=1,nr
            Y(k)=Y(k)*denom
        ENDDO
        DO k=1,nr
            xkk=xk(k)
            DO i=1,nr
                Xik(i,k)=Xik(i,k)+xkk*Y(i)
            ENDDO
        ENDDO    
    END SUBROUTINE
    
    SUBROUTINE GetPhi(phi,Xik,q,nr)    
        IMPLICIT NONE    
        INTEGER,INTENT(IN) :: nr
        REAL,INTENT(IN) :: Xik(nr,nr),q(nr)    
        REAL,INTENT(OUT) :: phi(nr)
        INTEGER :: i,k
        REAL :: qk
        DO i=1,nr
            phi(i)=0._8
        ENDDO
        DO k=1,nr
            qk=q(k)
            DO i=1,nr
                phi(i) = phi(i) + qk*Xik(i,k)
            ENDDO
        ENDDO
    END SUBROUTINE
    
end module