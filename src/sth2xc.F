      SUBROUTINE sth2x3(a,s,it,err) 
!   compute water thermodynamic properties as a function of temperature 
!   and pressure
      USE SteamTBL_mod, ONLY : nt,np,ns,ns2,klp,klp2,llp,nt5,jpl, mxst
      IMPLICIT real (8)(a-h,o-z) 
      REAL(8) a(mxst),s(26) 
      LOGICAL err 
!      INCLUDE 'machds.h' 
!                                                                       
      REAL(8) k(9) 
      LOGICAL s1,s2 
      INTEGER iunp(2) 
      EQUIVALENCE(unp,iunp(1)) 
      DATA k/-7.691234564d0,-2.608023696d1,-1.681706546d2,6.423285504d1,
     +-1.189646225d2,4.167117320d0,2.097506760d1,1.d9,6.0d0/            
      DATA crt/647.3d0/,crp/22120000.d0/ 
      IF(s(1).lt.a(1).or.s(1).gt.5000.d0) GOTO 1001 
      IF(s(2).le.0.0d0.or.s(2).gt.a(jpl)) GOTO 1001 
      unp=s(23) 
      ip=iunp(1) 
      jp=iunp(2) 
      IF(ip.le.0.or.ip.gt.nt)ip=1 
      IF(jp.le.0.or.jp.ge.np)jp=1 
      IF(s(1).ge.a(nt)) GOTO 46 
!   set indexes in temperature and pressure tables for saturation       
!   computations                                                        
   11 IF(s(1).ge.a(ip)) GOTO 10 
      ip=ip-1 
      GOTO 11 
   10 IF(s(1).lt.a(ip+1)) GOTO 12 
      ip=ip+1 
      GOTO 10 
   12 IF(ip.ge.ns) GOTO 44 
      s1=.false. 
      s2=.false. 
      fr=s(1)/crt 
      fr1=1.0d0-fr 
      s(10)=crp*exp((((((k(5)*fr1+k(4))*fr1+k(3))*fr1+k(2))*fr1+k(1))*  
     +fr1)/(((k(7)*fr1+k(6))*fr1+1.0d0)*fr)-fr1/(k(8)*fr1*fr1+k(9)))    
      jpp=jp+nt 
   13 IF(s(10).ge.a(jpp)) GOTO 14 
      jpp=jpp-1 
      IF(jpp.gt.nt) GOTO 13 
      s2=.true. 
      GOTO 15 
   14 IF(s(10).lt.a(jpp+1)) GOTO 15 
      jpp=jpp+1 
      IF(jpp.lt.jpl) GOTO 14 
      s1=.true. 
   15 kp2=klp2+(jpp-nt) * 13 
      kp=klp+ip * 13 
      IF(s2.or.a(jpp).le.a(kp)) GOTO 16 
      pa=a(jpp) 
      ta=a(kp2) 
      ia=kp2 
      GOTO 17 
   16 ta=a(ip) 
      pa=a(kp) 
      ia=kp 
   17 IF(s1.or.a(jpp+1).ge.a(kp+13)) GOTO 18 
      pb=a(jpp+1) 
      tb=a(kp2+13) 
      ib=kp2+13 
      GOTO 19 
   18 tb=a(ip+1) 
      pb=a(kp+13) 
      ib=kp+13 
   19 fr1=s(1)-ta 
!   compare input pressure to saturation pressure to determine vapor or 
!   liquid                                                              
      fr=fr1/(tb-ta) 
      IF(s(2).le.s(10)) GOTO 43 
      it=1 
      s(9)=0.0d0 
      GOTO 50 
   46 ip=nt 
   44 it=4 
      GOTO 45 
   43 it=3 
   45 s(9)=1.0d0 
   50 jpq=jp+nt 
!   search for single phase indexes.                                    
   51 IF(s(2).ge.a(jpq)) GOTO 53 
      jpq=jpq-1 
      IF(jpq.gt.nt) GOTO 51 
      GOTO 90 
   53 IF(s(2).lt.a(jpq+1)) GOTO 54 
      jpq=jpq+1 
      IF(jpq.lt.jpl) GOTO 53 
      GOTO 1001 
   54 jp=jpq-nt 
      lpp=llp+jp*nt5+ip * 6 
      lqq=lpp+nt5 
      kp2=klp2+jp * 13 
      frn=s(1)-a(ip) 
      IF(it-3)60,70,73 
!   liquid phase.                                                       
   60 IF(a(jpq).ge.s(10)) GOTO 61 
      px=s(10) 
      s2=.false. 
      GOTO 62 
   61 px=a(jpq) 
      s2=.true. 
      IF(jp.gt.ns2) GOTO 63 
      IF(a(kp2).gt.a(ip+1)) GOTO 63 
      frc=frn/(a(kp2)-a(ip)) 
      frc2=s(1)-a(kp2) 
      ic=kp2+1 
      GOTO 62 
   63 frc=frn/(a(ip+1)-a(ip)) 
      frc2=s(1)-a(ip+1) 
      ic=lpp+6 
      frd=frc 
      frd2=frc2 
      GOTO 65 
   62 IF(jp.ge.ns2) GOTO 67 
      IF(a(kp2+13).gt.a(ip+1)) GOTO 67 
      frd=frn/(a(kp2+13)-a(ip)) 
      frd2=s(1)-a(kp2+13) 
      id=kp2+14 
      GOTO 66 
   67 frd=frn/(a(ip+1)-a(ip)) 
      frd2=s(1)-a(ip+1) 
   65 id=lqq+6 
   66 IF(s2) GOTO 69 
      hfg1=a(ia+8)-a(ia+2)+pa*(a(ia+7)-a(ia+1)) 
      hfg2=a(ib+8)-a(ib+2)+pb*(a(ib+7)-a(ib+1)) 
      dpdt1=hfg1/(ta*(a(ia+7)-a(ia+1))) 
      dpdt2=hfg2/(tb*(a(ib+7)-a(ib+1))) 
      f1=a(ia+1)*(a(ia+3)-a(ia+4)*dpdt1) 
      f2=a(ib+1)*(a(ib+3)-a(ib+4)*dpdt2) 
      d1=f1*(tb-ta) 
      d2=f2*(tb-ta) 
      c0=a(ia+1) 
      c1=d1 
      c2=3.d0*(a(ib+1)-a(ia+1))-d2-2.d0*d1 
      c3=d2+d1-2.d0*(a(ib+1)-a(ia+1)) 
      s(11)=c0+fr*(c1+fr*(c2+fr*c3)) 
      s(13)=a(ia+2)+(a(ib+2)-a(ia+2))*fr 
      s(17)=a(ia+3)+(a(ib+3)-a(ia+3))*fr 
      s(19)=a(ia+4)+(a(ib+4)-a(ia+4))*fr 
      s(21)=a(ia+5)+(a(ib+5)-a(ia+5))*fr 
      s(25)=a(ia+6)+(a(ib+6)-a(ia+6))*fr 
      GOTO 56 
   69 d1=a(lpp+2)*a(lpp)*(frn-frc2) 
      d2=a(ic+2)*a(ic)*(frn-frc2) 
      c0=a(lpp) 
      c1=d1 
      c2=3.d0*(a(ic)-a(lpp))-d2-2.d0*d1 
      c3=d2+d1-2.d0*(a(ic)-a(lpp)) 
      s(11)=c0+frc*(c1+frc*(c2+frc*c3)) 
      s(13)=a(lpp+1)+(a(ic+1)-a(lpp+1))*frc 
      s(17)=a(lpp+2)+(a(ic+2)-a(lpp+2))*frc 
      s(19)=a(lpp+3)+(a(ic+3)-a(lpp+3))*frc 
      s(21)=a(lpp+4)+(a(ic+4)-a(lpp+4))*frc 
      s(25)=a(lpp+5)+(a(ic+5)-a(lpp+5))*frc 
   56 d1=a(lqq+2)*a(lqq)*(frn-frd2) 
      d2=a(id+2)*a(id)*(frn-frd2) 
      c0=a(lqq) 
      c1=d1 
      c2=3.d0*(a(id)-a(lqq))-d2-2.d0*d1 
      c3=d2+d1-2.d0*(a(id)-a(lqq)) 
      s(12)=c0+frd*(c1+frd*(c2+frd*c3)) 
      s(14)=a(lqq+1)+(a(id+1)-a(lqq+1))*frd 
      s(18)=a(lqq+2)+(a(id+2)-a(lqq+2))*frd 
      s(20)=a(lqq+3)+(a(id+3)-a(lqq+3))*frd 
      s(22)=a(lqq+4)+(a(id+4)-a(lqq+4))*frd 
      s(26)=a(lqq+5)+(a(id+5)-a(lqq+5))*frd 
      IF(s(11).gt.s(12)) GOTO 83 
      s(3)=s(11) 
      fr=0.0d0 
      GOTO 84 
   83 s(15)=s(19)*s(11) 
      s(16)=s(20)*s(12) 
      fr1=s(16)-s(15) 
      IF(abs(fr1).lt.1.0d-10) GOTO 81 
      fr=s(11)+s(12)-(a(jpq+1)-px)*s(15)*s(16)/fr1 
      fr1=sqrt(fr*fr-4.0d0*s(11)*s(12)*(s(16)*(1.0d0-s(19)*(s(2)-px))-s(
     +15)*(1.0d0-s(20)*(s(2)-a(jpq+1))))/fr1)                           
      s(3)=0.5d0*(fr+fr1) 
      IF(s(3).gt.s(11))s(3)=0.5d0*(fr-fr1) 
      IF(s(3).ge.s(12)) GOTO 82 
   81 fr=(s(2)-px)/(a(jpq+1)-px) 
      s(3)=s(11)*(1.0d0-fr)+s(12)*fr 
   82 fr=(s(3)-s(11))/(s(12)-s(11)) 
   84 s(4)=s(13)+(s(14)-s(13))*fr 
      s(5)=s(4)+s(2)*s(3) 
      s(6)=s(17)+(s(18)-s(17))*fr 
      s(7)=s(19)+(s(20)-s(19))*fr 
      s(8)=s(21)+(s(22)-s(21))*fr 
      s(24)=s(25)+(s(26)-s(25))*fr 
      GOTO 20 
!   vapor phase.                                                        
   70 s1=.false. 
      IF(a(jpq+1).le.s(10)) GOTO 71 
      frc=s(10) 
      hfg1=a(ia+8)-a(ia+2)+pa*(a(ia+7)-a(ia+1)) 
      hfg2=a(ib+8)-a(ib+2)+pb*(a(ib+7)-a(ib+1)) 
      dpdt1=hfg1/(ta*(a(ia+7)-a(ia+1))) 
      dpdt2=hfg2/(tb*(a(ib+7)-a(ib+1))) 
      f1=a(ia+7)*(a(ia+9)-a(ia+10)*dpdt1) 
      f2=a(ib+7)*(a(ib+9)-a(ib+10)*dpdt2) 
      d1=f1*(tb-ta) 
      d2=f2*(tb-ta) 
      c0=a(ia+7) 
      c1=d1 
      c2=3.d0*(a(ib+7)-a(ia+7))-d2-2.d0*d1 
      c3=d2+d1-2.d0*(a(ib+7)-a(ia+7)) 
      s(12)=c0+fr*(c1+fr*(c2+fr*c3)) 
      s(14)=a(ia+8)+(a(ib+8)-a(ia+8))*fr 
      s(18)=a(ia+9)+fr*tb/s(1)*(a(ib+9)-a(ia+9)) 
      s(20)=a(ia+10)+(s(10)-pa)/(pb-pa)*pb/s(10)*(a(ib+10)-a(ia+10)) 
      s(22)=a(ia+11)+(a(ib+11)-a(ia+11))*fr 
      s(26)=a(ia+12)+(a(ib+12)-a(ia+12))*fr 
      GOTO 72 
   71 IF(a(kp2+13).lt.a(ip)) GOTO 73 
      frd=(s(1)-a(kp2+13))/(a(ip+1)-a(kp2+13)) 
      ic=kp2+20 
      GOTO 74 
   73 IF(ip.eq.nt) GOTO 80 
      frd=frn/(a(ip+1)-a(ip)) 
      s1=.true. 
      ic=lqq 
   74 frc=a(jpq+1) 
      c0=a(ic) 
      d1=a(ic)*a(ic+2)*(a(ip+1)-a(ip)) 
      d2=a(lqq+8)*a(lqq+6)*(a(ip+1)-a(ip)) 
      c1=d1 
      c2=3.d0*(a(lqq+6)-a(ic))-d2-2.d0*d1 
      c3=d2+d1-2.d0*(a(lqq+6)-a(ic)) 
      s(12)=c0+frd*(c1+frd*(c2+frd*c3)) 
      s(14)=a(ic+1)+(a(lqq+7)-a(ic+1))*frd 
      s(18)=a(ic+2)+frd*a(ip+1)/s(1)*(a(lqq+8)-a(ic+2)) 
      s(20)=a(ic+3)+(a(lqq+9)-a(ic+3))*frd 
      s(22)=a(ic+4)+(a(lqq+10)-a(ic+4))*frd 
      s(26)=a(ic+5)+(a(lqq+11)-a(ic+5))*frd 
   72 IF(s1) GOTO 75 
      IF(a(kp2).lt.a(ip)) GOTO 77 
      frd=(s(1)-a(kp2))/(a(ip+1)-a(kp2)) 
      ia=kp2+7 
      GOTO 76 
   77 frd=frn/(a(ip+1)-a(ip)) 
   75 ia=lpp 
   76 c0=a(ia) 
      d1=a(ia)*a(ia+2)*(a(ip+1)-a(ip)) 
      d2=a(lpp+8)*a(lpp+6)*(a(ip+1)-a(ip)) 
      c1=d1 
      c2=3.d0*(a(lpp+6)-a(ia))-d2-2.d0*d1 
      c3=d2+d1-2.d0*(a(lpp+6)-a(ia)) 
      s(11)=c0+frd*(c1+frd*(c2+frd*c3)) 
      s(13)=a(ia+1)+(a(lpp+7)-a(ia+1))*frd 
      s(17)=a(ia+2)+frd*a(ip+1)/s(1)*(a(lpp+8)-a(ia+2)) 
      s(19)=a(ia+3)+(a(lpp+9)-a(ia+3))*frd 
      s(21)=a(ia+4)+(a(lpp+10)-a(ia+4))*frd 
      s(25)=a(ia+5)+(a(lpp+11)-a(ia+5))*frd 
      fr=s(12)*(frc-a(jpq)) 
      s(3)=s(11)*fr/(fr+(s(11)-s(12))*(s(2)-a(jpq))) 
      fr=(s(3)-s(11))/(s(12)-s(11)) 
      frn=fr*s(12)/s(3) 
      s(4)=s(13)+(s(14)-s(13))*fr 
      s(5)=s(4)+s(2)*s(3) 
      s(6)=s(17)+frn*(s(18)-s(17)) 
      s(7)=s(19)+(s(20)-s(19))*fr 
      s(8)=s(21)+frn*(s(22)-s(21)) 
      s(24)=s(25)+(s(26)-s(25))*fr 
   20 iunp(1)=ip 
      iunp(2)=jp 
      s(23)=unp 
      err=.false. 
      RETURN 
!   vapor phase, temperature greater than highest table temperature.    
   80 fr=a(lqq)*(a(jpq+1)-a(jpq)) 
      s(3)=a(lpp)*fr/(fr+(a(lpp)-a(lqq))*(s(2)-a(jpq))) 
      fr=(s(3)-a(lpp))/(a(lqq)-a(lpp)) 
      frc=fr*a(lqq)/s(3) 
      s(5)=a(lpp+1)+(a(lqq+1)-a(lpp+1))*fr+s(2)*s(3) 
      s(8)=a(lpp+4)+frc*(a(lqq+4)-a(lpp+4)) 
      frd=s(1)/a(nt) 
      s(3)=s(3)*frd 
      s(5)=s(5)+s(8)*frn 
      s(4)=s(5)-s(2)*s(3) 
      s(6)=(a(lpp+2)+frc*(a(lqq+2)-a(lpp+2)))/frd 
      s(7)=a(lpp+3)+(a(lqq+3)-a(lpp+3))*fr 
      s(24)=a(lpp+5)+(a(lqq+5)-a(lpp+5))*fr 
      cv=s(8)-a(nt)*s(6)*s(6)*s(3)/s(7) 
      s(24)=s(24)+cv*log(frd**(s(8)/cv)) 
      GOTO 20 
!   vapor phase, pressure less than lowest table pressure               
   90 IF(it.eq.1) GOTO 1001 
      IF(s(1).lt.a(klp2+13)) GOTO 92 
      lpp=llp+nt5+ip * 6 
      IF(ip.eq.nt) GOTO 95 
      IF(a(ip).lt.a(klp2+13)) GOTO 93 
      ia=ip 
      lqq=lpp 
      GOTO 91 
   93 ia=klp2+13 
      lqq=ia+7 
   91 fr=(s(1)-a(ia))/(a(ip+1)-a(ia)) 
      s(3)=(fr*a(lpp+6)/a(ip+1)+(1.0d0-fr)*a(lqq)/a(ia))*a(jpq+1)*s(1)/ 
     +s(2)                                                              
      s(4)=a(lqq+1)+(a(lpp+7)-a(lqq+1))*fr 
      s(6)=a(lqq+2)+(a(lpp+8)-a(lqq+2))*fr*a(ip+1)/s(1) 
      s(8)=a(lqq+4)+(a(lpp+10)-a(lqq+4))*fr 
      s(24)=a(lqq+5)+(a(lpp+11)-a(lqq+5))*fr 
      ren=s(2)*s(3)/s(1) 
      s(24)=s(24)-ren*log(s(2)/a(jpq+1)) 
   94 s(5)=s(4)+s(2)*s(3) 
      s(7)=1.0d0/s(2) 
      jp=1 
      GOTO 20 
   92 s(3)=(fr*pb*a(ib+7)/tb+(1.0d0-fr)*pa*a(ia+7)/ta)*s(1)/s(2) 
      s(4)=a(ia+8)+(a(ib+8)-a(ia+8))*fr 
      s(6)=1.0d0/s(1) 
      s(8)=a(ia+11)+(a(ib+11)-a(ia+11))*fr 
      s(24)=a(ia+12)+(a(ib+12)-a(ia+12))*fr 
      ren=s(2)*s(3)/s(1) 
      s(24)=s(24)-ren*log(s(2)/s(10)) 
      GOTO 94 
   95 frd=s(1)/a(nt) 
      frc=a(nt+1)*a(lpp) 
      s(3)=frc*frd/s(2) 
      s(8)=a(lpp+4) 
      s(5)=a(lpp+1)+frc+s(8)*(s(1)-a(nt)) 
      s(4)=s(5)-s(2)*s(3) 
      s(6)=a(lpp+2)/frd 
      s(7)=a(lpp+3)*a(nt+1)/s(2) 
      ren=s(2)*s(3)/s(1) 
      s(24)=a(lpp+5)+s(8)*log(s(1)/a(nt))-ren*log(s(2)/a(jpq+1)) 
      GOTO 20 
 1001 err=.true. 
      RETURN 
      END SUBROUTINE sth2x3                         

      SUBROUTINE sth2x5(a,s,it,err) 
!   compute water thermodynamic properties as a function of enthalpy 
!   and pressure                                                        
      USE SteamTBL_mod, ONLY : nt,np,ns,ns2,klp,klp2,llp,nt5,jpl, mxst
      IMPLICIT real (8)(a-h,o-z) 
!
      LOGICAL err 
!     INCLUDE 'machds.h' 
!                                                                       
      dimension a(mxst),s(26)
      real(4)  k(9),c(9),b(6),cc(3) 
      LOGICAL s1,s2,s3 
!
      INTEGER iunp(2) 
      EQUIVALENCE(unp,iunp(1)) 
      DATA b/6669.352222d0,-4658.899d0,1376.536722d0,-201.9126167d0,    
     +14.82832111d0,-.4337434056d0/                                     
      DATA c/274.9043833d0,13.66254889d0,1.176781611d0,-.189693d0,      
     +8.74535666d-2,-1.7405325d-2,2.147682333d-3,-1.383432444d-4,       
     +3.800086611d-6/                                                   
      DATA crt/647.3d0/,crp/22120000.d0/,ctrp/273.16d0/ 
      DATA k/-7.691234564d0,-26.08023696d0,-168.1706546d0,6.423285504d1,
     +-1.189646225d2,4.167117320d0,2.097506760d1,1.d9,6.0d0/            
      DATA pxxx/1378.951459d0/,pxxy/3102640.782d0/,pxx1/1.450377377d-3/,
     +pxx2/1.450377377d-4/                                              
      DATA cc/0.84488898d0,2.9056480d0,219.74589d0/ 
!                                                                       
      s3=.false. 
!   check for valid input                                               
   23 IF(s(2).le.0.0d0.or.s(2).gt.a(jpl)) GOTO 1001 
      unp=s(23) 
      ip=iunp(1) 
      jp=iunp(2) 
      IF(ip.le.0.or.ip.ge.nt)ip=1 
      IF(jp.le.0.or.jp.ge.np)jp=1 
      jpp=jp+nt 
      s1=.false. 
!   set indexes in temperature and pressure tables for saturation       
!   computations                                                        
   11 IF(s(2).ge.a(jpp)) GOTO 10 
      jpp=jpp-1 
      IF(jpp.gt.nt) GOTO 11 
      jpp=jpp+1 
      jp=1 
      s1=.true. 
      IF(s(2).lt.a(klp+13)) GOTO 44 
      GOTO 12 
   10 IF(s(2).lt.a(jpp+1)) GOTO 12 
      jpp=jpp+1 
      GOTO 10 
   12 jp=jpp-nt 
      IF(s(2).ge.crp) GOTO 44 
      IF(s3) GOTO 15 
      IF(s(2).lt.pxxx) GOTO 18 
      IF(s(2).gt.pxxy) GOTO 17 
      fr=log(pxx1*s(2)) 
      s(10)=(((((((fr*c(9)+c(8))*fr+c(7))*fr+c(6))*fr+c(5))*fr+c(4))*fr+
     +c(3))*fr+c(2))*fr+c(1)                                            
      GOTO 27 
   17 fr=log(pxx2*s(2)) 
      s(10)=((((fr*b(6)+b(5))*fr+b(4))*fr+b(3))*fr+b(2))*fr+b(1) 
      GOTO 27 
   18 fr=log(s(2)) 
      s(10)=(cc(1)*fr+cc(2))*fr+cc(3) 
   27 fr=s(10)/crt 
      fr1=1.0d0-fr 
      d1=((((k(5)*fr1+k(4))*fr1+k(3))*fr1+k(2))*fr1+k(1))*fr1 
      d2=(((5.0d0*k(5)*fr1+4.0d0*k(4))*fr1+3.0d0*k(3))*fr1+2.0d0*k(2))* 
     +fr1+k(1)                                                          
      c2=k(7)*fr1 
      c1=(c2+k(6))*fr1+1.0d0 
      c2=2.0d0*c2+k(6) 
      f2=k(8)*fr1 
      f1=1.0d0/(f2*fr1+k(9)) 
      f2=2.0d0*f2 
      hfg1=1.0d0/(fr*c1) 
      hfg2=fr1*f1 
      pa=crp*exp(d1*hfg1-hfg2) 
      s(10)=max(s(10)+(s(2)-pa)*crt/(pa*((d1*hfg1*(fr*c2-c1)-d2)*hfg1+( 
     +1.0d0-hfg2*f2)*f1)),ctrp)                                         
   15 ic=ip 
   16 IF(s(10).ge.a(ic)) GOTO 13 
      ic=ic-1 
      IF(ic.gt.0) GOTO 16 
      ic=1 
      GOTO 14 
   13 IF(s(10).lt.a(ic+1)) GOTO 14 
      ic=ic+1 
      IF(ic.lt.ns) GOTO 13 
      ic=ic-1 
   14 kp=klp+ic * 13 
      kp2=klp2+jp * 13 
      IF(s1.or.a(jpp).le.a(kp)) GOTO 19 
      pa=a(jpp) 
      ta=a(kp2) 
      ia=kp2 
      GOTO 20 
   19 ta=a(ic) 
      pa=a(kp) 
      ia=kp 
   20 IF(a(jpp+1).ge.a(kp+13)) GOTO 21 
      pb=a(jpp+1) 
      tb=a(kp2+13) 
      ib=kp2+13 
      GOTO 22 
   21 pb=a(kp+13) 
      tb=a(ic+1) 
      ib=kp+13 
   22 fr1=s(10)-ta 
      fr=fr1/(tb-ta) 
!   compute vsubf and vsubg to determine liquid, two phase, or vapor    
!   state                                                               
      hfg1=a(ia+8)-a(ia+2)+pa*(a(ia+7)-a(ia+1)) 
      hfg2=a(ib+8)-a(ib+2)+pb*(a(ib+7)-a(ib+1)) 
      dpdt1=hfg1/(ta*(a(ia+7)-a(ia+1))) 
      dpdt2=hfg2/(tb*(a(ib+7)-a(ib+1))) 
      f1=a(ia+1)*(a(ia+3)-a(ia+4)*dpdt1) 
      f2=a(ib+1)*(a(ib+3)-a(ib+4)*dpdt2) 
      c1=f1*(tb-ta) 
      d2=f2*(tb-ta) 
      c0=a(ia+1) 
      c2=3.d0*(a(ib+1)-c0)-d2-2.d0*c1 
      c3=d2+c1-2.d0*(a(ib+1)-c0) 
      s(11)=c0+fr*(c1+fr*(c2+fr*c3)) 
   24 s(13)=a(ia+2)+(a(ib+2)-a(ia+2))*fr 
      s(15)=s(13)+s(2)*s(11) 
      IF(s(5).le.s(15)) GOTO 41 
      f1=a(ia+7)*(a(ia+9)-a(ia+10)*dpdt1) 
      f2=a(ib+7)*(a(ib+9)-a(ib+10)*dpdt2) 
      c1=f1*(tb-ta) 
      d2=f2*(tb-ta) 
      c0=a(ia+7) 
      c2=3.d0*(a(ib+7)-c0)-d2-2.d0*c1 
      c3=d2+c1-2.d0*(a(ib+7)-c0) 
      s(12)=c0+fr*(c1+fr*(c2+fr*c3)) 
      s(14)=a(ia+8)+(a(ib+8)-a(ia+8))*fr 
      s(16)=s(14)+s(2)*s(12) 
      IF(s(5).ge.s(16)) GOTO 43 
!   two phase fluid.                                                    
      it=2 
      s(17)=a(ia+3)+(a(ib+3)-a(ia+3))*fr 
      s(18)=a(ia+9)+fr*tb/s(10)*(a(ib+9)-a(ia+9)) 
      s(19)=a(ia+4)+(a(ib+4)-a(ia+4))*fr 
      s(20)=a(ia+10)+(s(2)-pa)/(pb-pa)*pb/s(2)*(a(ib+10)-a(ia+10)) 
      s(21)=a(ia+5)+(a(ib+5)-a(ia+5))*fr 
      s(22)=a(ia+11)+(a(ib+11)-a(ia+11))*fr 
      s(25)=a(ia+6)+(a(ib+6)-a(ia+6))*fr 
      s(26)=a(ia+12)+(a(ib+12)-a(ia+12))*fr 
      s(9)=(s(5)-s(15))/(s(16)-s(15)) 
      fr=1.0d0-s(9) 
      s(1)=s(10) 
      s(3)=fr*s(11)+s(9)*s(12) 
      s(4)=fr*s(13)+s(9)*s(14) 
      s(24)=fr*s(25)+s(9)*s(26) 
      ip=ic 
   25 iunp(1)=ip 
      iunp(2)=jp 
      s(23)=unp 
      err=.false. 
      RETURN 
   44 it=4 
      GOTO 50 
!   single phase fluid, search for single phase indexes.                
   41 it=1 
      IF(s1) GOTO 1001 
      GOTO 50 
   43 it=3 
   50 lpp=llp+jp*nt5+ip * 6 
   51 ht=a(lpp+1)+a(jpp)*a(lpp) 
      IF(s(5).ge.ht) GOTO 52 
      lpp=lpp-6 
      ip=ip-1 
      IF(ip.gt.0) GOTO 51 
      lpp=lpp+6 
      ip=ip+1 
      GOTO 54 
   52 ht=a(lpp+7)+a(jpp)*a(lpp+6) 
      IF(s(5).le.ht) GOTO 54 
      lpp=lpp+6 
      ip=ip+1 
      IF(ip.lt.nt) GOTO 52 
      IF(s1) GOTO 95 
      lpp=lpp-6 
      ip=ip-1 
      GOTO 53 
   54 IF(s1) GOTO 96 
   53 kp=klp+ip * 13 
      lqq=lpp+nt5 
      frn=s(2)-a(jpp) 
      frc2=s(2)-a(jpp+1) 
      s3=.false. 
      IF(it-3)58,70,48 
   48 IF(ip.ge.ns) GOTO 157 
      GOTO 57 
!   liquid phase.                                                       
   58 ht=s(15) 
      ut=s(13) 
      vt=s(11) 
   60 IF(a(ip).lt.s(10)) GOTO 57 
      ip=ip-1 
      IF(ip.le.0) GOTO 1001 
      kp=kp-13 
      lpp=lpp-6 
      lqq=lqq-6 
      GOTO 60 
   57 s1=.false. 
      IF(it.eq.4) GOTO 61 
      IF(a(ip+1).le.s(10)) GOTO 61 
      ta=s(10) 
      s2=.false. 
      GOTO 62 
   61 ta=a(ip+1) 
      s2=.true. 
      IF(ip+1.gt.ns) GOTO 63 
      IF(a(kp+13).le.a(jpp)) GOTO 63 
      frc3=s(2)-a(kp+13) 
      px=a(kp+13) 
      frc=a(jpp+1)-a(kp+13) 
      ic=kp+14 
      GOTO 64 
   63 frc3=frn 
      frc=a(jpp+1)-a(jpp) 
      px=a(jpp) 
      ic=lpp+6 
      frd3=frc3 
      frd=frc 
      s1=.true. 
   64 pr=(s(2)-px)/(a(jpp+1)-px) 
      c0=1.d0/a(ic) 
      r2=1.d0/a(lqq+6) 
      c1=c0*a(ic+3)*(a(jpp+1)-px) 
      d2=r2*a(lqq+9)*(a(jpp+1)-px) 
      c2=3.d0*(r2-c0)-d2-2.d0*c1 
      c3=d2+c1-2.d0*(r2-c0) 
      s11i=c0+pr*(c1+pr*(c2+pr*c3)) 
      s(11)=1.d0/s11i 
   82 fr1=(s(11)-a(ic))/(a(lqq+6)-a(ic)) 
      s(13)=a(ic+1)+(a(lqq+7)-a(ic+1))*fr1 
      s(15)=s(13)+s(2)*s(11) 
   62 IF(s1) GOTO 65 
      IF(ip.gt.ns) GOTO 67 
      IF(a(kp).le.a(jpp)) GOTO 67 
      frd3=s(2)-a(kp) 
      frd=a(jpp+1)-a(kp) 
      px=a(kp) 
      id=kp+1 
      GOTO 66 
   67 frd3=frn 
      frd=a(jpp+1)-a(jpp) 
      px=a(jpp) 
      s1=.true. 
   65 id=lpp 
   66 IF(s3) GOTO 68 
      pr=(s(2)-px)/(a(jpp+1)-px) 
      c0=1.d0/a(id) 
      r2=1.d0/a(lqq) 
      c1=c0*a(id+3)*(a(jpp+1)-px) 
      d2=r2*a(lqq+3)*(a(jpp+1)-px) 
      c2=3.d0*(r2-c0)-d2-2.d0*c1 
      c3=d2+c1-2.d0*(r2-c0) 
      s12i=c0+pr*(c1+pr*(c2+pr*c3)) 
      s(12)=1.d0/s12i 
  182 fr2=(s(12)-a(id))/(a(lqq)-a(id)) 
      s(14)=a(id+1)+(a(lqq+1)-a(id+1))*fr2 
      s(16)=s(14)+s(2)*s(12) 
      IF(s(16).le.s(5)) GOTO 68 
      s2=.true. 
      s(11)=s(12) 
      s(13)=s(14) 
      s(15)=s(16) 
      fr1=fr2 
      ip=ip-1 
      ta=a(ip+1) 
      IF(ip.le.0) GOTO 1001 
      kp=kp-13 
      lqq=lqq-6 
      lpp=lpp-6 
      ic=id 
      GOTO 62 
   68 IF(s(15).ge.s(5)) GOTO 59 
      s(12)=s(11) 
      s(11)=vt 
      s(14)=s(13) 
      s(13)=ut 
      s(16)=s(15) 
      s(15)=ht 
      fr2=fr1 
      lqq=lqq+6 
      lpp=lpp+6 
      kp=kp+13 
      ip=ip+1 
      IF(ip.lt.ns) GOTO 158 
      s3=.false. 
      GOTO 157 
  158 id=ic 
      s3=.true. 
      GOTO 57 
   59 IF(s2) GOTO 69 
      s(17)=a(ia+3)+(a(ib+3)-a(ia+3))*fr 
      s(19)=a(ia+4)+(a(ib+4)-a(ia+4))*fr 
      s(21)=a(ia+5)+(a(ib+5)-a(ia+5))*fr 
      s(25)=a(ia+6)+(a(ib+6)-a(ia+6))*fr 
      GOTO 56 
   69 s(17)=a(ic+2)+(a(lqq+8)-a(ic+2))*fr1 
      s(19)=a(ic+3)+(a(lqq+9)-a(ic+3))*fr1 
      s(21)=a(ic+4)+(a(lqq+10)-a(ic+4))*fr1 
      s(25)=a(ic+5)+(a(lqq+11)-a(ic+5))*fr1 
   56 s(18)=a(id+2)+(a(lqq+2)-a(id+2))*fr2 
      s(20)=a(id+3)+(a(lqq+3)-a(id+3))*fr2 
      s(22)=a(id+4)+(a(lqq+4)-a(id+4))*fr2 
      s(26)=a(id+5)+(a(lqq+5)-a(id+5))*fr2 
      fr=(s(5)-s(16))/(s(15)-s(16)) 
      frd=ta-a(ip) 
      fr2=frd*fr 
      s(1)=a(ip)+fr2 
      fr1=fr2/frd 
      tr=(s(1)-a(ip))/frd 
      c0=s(12) 
      c1=s(12)*s(18)*frd 
      d2=s(11)*s(17)*frd 
      c2=3.d0*(s(11)-s(12))-d2-2.d0*c1 
      c3=d2+c1-2.d0*(s(11)-s(12)) 
      s(3)=c0+tr*(c1+tr*(c2+tr*c3)) 
      s(4)=s(5)-s(2)*s(3) 
      s(6)=(c1+tr*(2.d0*c2+3.d0*tr*c3))/(frd*s(3)) 
      s(7)=s(20)+(s(19)-s(20))*fr1 
      s(8)=s(22)+(s(21)-s(22))*fr1 
      s(24)=s(26)+(s(25)-s(26))*fr1 
      s(9)=0.d0 
      GOTO 25 
   70 ht=s(16) 
      ut=s(14) 
      vt=s(12) 
!   vapor phase.                                                        
  160 IF(a(ip+1).gt.s(10)) GOTO 157 
      ip=ip+1 
      kp=kp+13 
      lpp=lpp+6 
      lqq=lqq+6 
      GOTO 160 
  157 s1=.false. 
      IF(it.eq.4) GOTO 161 
      IF(a(ip).ge.s(10)) GOTO 161 
      ta=s(10) 
      s2=.false. 
      GOTO 162 
  161 ta=a(ip) 
      s2=.true. 
      IF(ip.ge.ns) GOTO 163 
      IF(a(kp).ge.a(jpp+1)) GOTO 163 
      frc3=a(kp)-a(jpp) 
      px=a(kp) 
      ic=kp+7 
      GOTO 164 
  163 frc3=a(jpp+1)-a(jpp) 
      px=a(jpp+1) 
      ic=lqq 
      frd3=frc3 
      s1=.true. 
  164 fr1=a(ic)*frc3 
      pr=(s(2)-a(jpp))/(px-a(jpp)) 
      c0=1.d0/a(lpp) 
      r2=1.d0/a(ic) 
      c1=c0*a(lpp+3)*(px-a(jpp)) 
      d2=r2*a(ic+3)*(px-a(jpp)) 
      c2=3.d0*(r2-c0)-d2-2.d0*c1 
      c3=d2+c1-2.d0*(r2-c0) 
      s12i=c0+pr*(c1+pr*(c2+pr*c3)) 
      s(12)=1.d0/s12i 
      frc2=(s(12)-a(lpp))/(a(ic)-a(lpp)) 
      frc=frc2*a(ic)/s(12) 
      s(14)=a(lpp+1)+(a(ic+1)-a(lpp+1))*frc2 
      s(16)=s(14)+s(2)*s(12) 
  162 IF(s1) GOTO 165 
      IF(ip+1.gt.ns) GOTO 167 
      IF(a(kp+13).ge.a(jpp+1)) GOTO 167 
      frd3=a(kp+13)-a(jpp) 
      px=a(kp+13) 
      id=kp+20 
      GOTO 166 
  167 frd3=a(jpp+1)-a(jpp) 
      px=a(jpp+1) 
      s1=.true. 
  165 id=lqq+6 
  166 IF(s3) GOTO 168 
      fr1=a(id)*frd3 
      pr=(s(2)-a(jpp))/(px-a(jpp)) 
      c0=1.d0/a(lpp+6) 
      r2=1.d0/a(id) 
      c1=c0*a(lpp+9)*(px-a(jpp)) 
      d2=r2*a(id+3)*(px-a(jpp)) 
      c2=3.d0*(r2-c0)-d2-2.d0*c1 
      c3=d2+c1-2.d0*(r2-c0) 
      s11i=c0+pr*(c1+pr*(c2+pr*c3)) 
      s(11)=1.d0/s11i 
      frd2=(s(11)-a(lpp+6))/(a(id)-a(lpp+6)) 
      frd=frd2*a(id)/s(11) 
      s(13)=a(lpp+7)+(a(id+1)-a(lpp+7))*frd2 
      s(15)=s(13)+s(2)*s(11) 
      IF(s(15).ge.s(5)) GOTO 168 
      s2=.true. 
      ip=ip+1 
      lqq=lqq+6 
      lpp=lpp+6 
      IF(ip.eq.nt) GOTO 80 
      ta=a(ip) 
      kp=kp+13 
      ic=id 
      s(12)=s(11) 
      s(14)=s(13) 
      s(16)=s(15) 
      frc2=frd2 
      frc=frd 
      GOTO 162 
  168 IF(s(16).le.s(5)) GOTO 159 
      s3=.true. 
      s(11)=s(12) 
      s(12)=vt 
      s(13)=s(14) 
      s(14)=ut 
      s(15)=s(16) 
      s(16)=ht 
      frd2=frc2 
      frd=frc 
      ip=ip-1 
      IF(ip.le.0) GOTO 1001 
      kp=kp-13 
      lpp=lpp-6 
      lqq=lqq-6 
      IF(ip.lt.ns) GOTO 57 
      GOTO 157 
  159 IF(s2) GOTO 169 
      s(18)=a(ia+9)+fr*tb/s(10)*(a(ib+9)-a(ia+9)) 
      s(20)=a(ia+10)+(s(2)-pa)/(pb-pa)*pb/s(2)*(a(ib+10)-a(ia+10)) 
      s(22)=a(ia+11)+(a(ib+11)-a(ia+11))*fr 
      s(26)=a(ia+12)+(a(ib+12)-a(ia+12))*fr 
      GOTO 156 
  169 s(18)=a(lpp+2)+frc*(a(ic+2)-a(lpp+2)) 
      s(20)=a(lpp+3)+(a(ic+3)-a(lpp+3))*frc2 
      s(22)=a(lpp+4)+frc*(a(ic+4)-a(lpp+4)) 
      s(26)=a(lpp+5)+(a(ic+5)-a(lpp+5))*frc2 
  156 s(17)=a(lpp+8)+frd*(a(id+2)-a(lpp+8)) 
      s(19)=a(lpp+9)+(a(id+3)-a(lpp+9))*frd2 
      s(21)=a(lpp+10)+frd*(a(id+4)-a(lpp+10)) 
      s(25)=a(lpp+11)+(a(id+5)-a(lpp+11))*frd2 
      fr=(s(5)-s(16))/(s(15)-s(16)) 
      frd=a(ip+1)-ta 
      fr2=frd*fr 
      s(1)=ta+fr2 
      fr1=fr2/frd 
      tr=(s(1)-ta)/frd 
      c0=s(12) 
      c1=s(12)*s(18)*frd 
      d2=s(11)*s(17)*frd 
      c2=3.d0*(s(11)-s(12))-d2-2.d0*c1 
      c3=d2+c1-2.d0*(s(11)-s(12)) 
      s(3)=c0+tr*(c1+tr*(c2+tr*c3)) 
      s(4)=s(5)-s(2)*s(3) 
      s(6)=(c1+tr*(2.d0*c2+3.d0*tr*c3))/(frd*s(3)) 
      s(7)=s(20)+(s(19)-s(20))*fr1 
      s(8)=s(22)+(s(21)-s(22))*fr1 
      s(24)=s(26)+(s(25)-s(26))*fr1 
   99 s(9)=1.0d0 
      GOTO 25 
!   vapor phase, temperature greater than highest table temperature.    
   80 fr=a(lqq)*(a(jpp+1)-a(jpp)) 
      s(3)=a(lpp)*fr/(fr+(a(lpp)-a(lqq))*frn) 
      fr=(s(3)-a(lpp))/(a(lqq)-a(lpp)) 
      frc=fr*a(lqq)/s(3) 
      ht=a(lpp+1)+(a(lqq+1)-a(lpp+1))*fr+s(2)*s(3) 
      s(8)=a(lpp+4)+frc*(a(lqq+4)-a(lpp+4)) 
      s(1)=(s(5)-ht+s(8)*a(nt))/s(8) 
      frd=s(1)/a(nt) 
      s(3)=s(3)*frd 
      s(4)=s(5)-s(2)*s(3) 
      s(6)=(a(lpp+2)+frc*(a(lqq+2)-a(lpp+2)))/frd 
      s(7)=a(lpp+3)+(a(lqq+3)-a(lpp+3))*fr 
      s(24)=a(lpp+5)+(a(lqq+5)-a(lpp+5))*fr 
      cv=s(8)-a(nt)*s(6)*s(6)*s(3)/s(7) 
      s(24)=s(24)+cv*log(frd**(s(8)/cv)) 
      GOTO 99 
!   vapor phase, pressure less than lowest table pressure               
   96 ht=a(klp2+21)+a(nt+1)*a(klp2+20) 
      IF(s(5).lt.ht) GOTO 90 
      IF(a(ip).lt.a(klp2+13)) GOTO 93 
      ia=ip 
      lqq=lpp 
      ht=a(lqq+1)+a(nt+1)*a(lqq) 
      GOTO 91 
   93 ia=klp2+13 
      lqq=ia+7 
   91 fr=(s(5)-ht)/(a(lpp+7)+a(nt+1)*a(lpp+6)-ht) 
      fr1=a(ip+1)-a(ia) 
      fr2=fr1*fr 
      s(1)=a(ia)+fr2 
      fr1=fr2/fr1 
      s(3)=(fr1*a(lpp+6)/a(ip+1)+(1.0d0-fr1)*a(lqq)/a(ia))*a(nt+1)*s(1)/
     +s(2)                                                              
      s(6)=a(lqq+2)+(a(lpp+8)-a(lqq+2))*fr*a(ip+1)/s(1) 
      s(8)=a(lqq+4)+(a(lpp+10)-a(lqq+4))*fr 
      s(24)=a(lqq+5)+(a(lpp+11)-a(lqq+5))*fr 
      ren=s(2)*s(3)/s(1) 
      s(24)=s(24)-ren*log(s(2)/a(nt+1)) 
   94 s(7)=1.0d0/s(2) 
   98 s(4)=s(5)-s(2)*s(3) 
      GOTO 99 
   95 s(8)=a(lpp+4) 
      s(1)=(s(5)-a(lpp+1)-a(nt+1)*a(lpp)+s(8)*a(nt))/s(8) 
      frd=s(1)/a(nt) 
      frc=a(nt+1)*a(lpp) 
      s(3)=frc*frd/s(2) 
      s(6)=a(lpp+2)/frd 
      s(7)=a(lpp+3)*a(nt+1)/s(2) 
      ren=s(2)*s(3)/s(1) 
      IF(s(1).le.0.0d0) GOTO 1001 
      s(24)=a(lpp+5)+s(8)*log(s(1)/a(nt))-ren*log(s(2)/a(nt+1)) 
      GOTO 98 
   90 ht=a(klp+21)+a(klp+13)*a(klp+20) 
      IF(s(5).lt.ht) GOTO 1001 
      ip=1 
      kp=klp+13 
  202 frd=a(kp+21)+a(kp+13)*a(kp+20) 
      IF(s(5).le.frd) GOTO 201 
      ip=ip+1 
      kp=kp+13 
      ht=frd 
      GOTO 202 
  201 fr=(s(5)-ht)/(frd-ht) 
      s(1)=a(ip)+fr*(a(ip+1)-a(ip)) 
      s(3)=(fr*a(kp+13)*a(kp+20)/a(ip+1)+(1.0d0-fr)*a(kp)*a(kp+7)/a(ip))
     +*s(1)/s(2)                                                        
      s(6)=1.0d0/s(1) 
      s(8)=a(kp+11)+(a(kp+24)-a(kp+11))*fr 
      ren=s(2)*s(3)/s(1) 
      IF(s(1).le.0.0d0) GOTO 1001 
      s(24)=a(kp+12)+s(8)*log(s(1)/a(ip))-ren*log(s(2)/a(kp)) 
      GOTO 94 
 1001 err=.true. 
      RETURN 
      END subroutine
