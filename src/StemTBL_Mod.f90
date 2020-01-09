MODULE SteamTBL_mod
USE PARAM
IMPLICIT NONE
INTEGER, PARAMETER :: mxst = 16140
REAL(8) :: stdata(mxst)
INTEGER :: nt,np,nst,nsp,it3bp,it4bp,it5bp,nprpnt,it3p0
REAL :: ttrip, ptrip, vtrip, tcrit, pcrit, vcrit, tmin, pmin, tmax, pmax
INTEGER :: ns,ns2,klp,klp2,llp,nt5,jpl 
real(8) state(26) !current state 
real(8) stt,stp,stv,stu,sth,stbeta,stkapa,stcp,stx,stpsat,stvl,stvv,stul,stuv,sthl, &
        sthv,stbetal,stbetav,stkapal,stkapv, stcpl,stcpv,stlv,sts,stsl, stsv, stkapav
character(80) STeamTBL_record(2)       
EQUIVALENCE(ns, nt); EQUIVALENCE(nsp, ns2)
EQUIVALENCE(it3bp, klp); EQUIVALENCE(it4bp, klp2)
EQUIVALENCE(it5bp, llp); EQUIVALENCE(nprpnt, nt5)
EQUIVALENCE(it3p0, jpl)


equivalence(state( 1),stt)    !temperature (K)
equivalence(state( 2),stp)    !pressure    (Pa)
equivalence(state( 3),stv)    !specific volume (m3/Kg)
equivalence(state( 4),stu)    !internal energy (J)
equivalence(state( 5),sth)    !specific enthalpy (J/Kg)
equivalence(state( 6),stbeta) !thermal exp. eoef. 
equivalence(state( 7),stkapa) !compressibility (1/Pa)
equivalence(state( 8),stcp)   !Cp (J/Kg)
equivalence(state( 9),stx)    !quality
equivalence(state(10),stpsat) !saturation pressure (Pa)
equivalence(state(11),stvl)   
equivalence(state(12),stvv)
equivalence(state(13),stul)
equivalence(state(14),stuv)
equivalence(state(15),sthl)
equivalence(state(16),sthv)
equivalence(state(17),stbetal)
equivalence(state(18),stbetav)
equivalence(state(19),stkapal)
equivalence(state(20),stkapav)
equivalence(state(21),stcpl)
equivalence(state(22),stcpv)
equivalence(state(23),stlv)   !two-phase indicator
equivalence(state(24),sts)    !specific entropy (J/Kg-K)
equivalence(state(25),stsl)
equivalence(state(26),stsv)

INTERFACE
  SUBROUTINE sth2x3(a,s,it,err) 
  INTEGER, PARAMETER :: mxst = 16140
  REAL(8) a(mxst),s(26) 
  INTEGER :: it
  LOGICAL err 
  END SUBROUTINE

  SUBROUTINE sth2x5(a,s,it,err) 
  INTEGER, PARAMETER :: mxst = 16140
  REAL(8) a(mxst),s(26) 
  INTEGER :: it
  LOGICAL err 
  END SUBROUTINE
    
END INTERFACE
CONTAINS

SUBROUTINE steamtbl(iftemp,wp,wt,wh,wrho,wv,wx,wbeta,wkapa,wcp)
USE PARAM
LOGICAL :: iftemp, err
REAL :: wp,wt,wh,wrho,wv,wx,wbeta,wkapa,wcp
INTEGER :: isterr
stp=wp
if(iftemp) then
  stt=wt
  call sth2x3(stdata,state,isterr,err)  ! function of temperature and pressure
  wh=sth
else
  sth=wh
  call sth2x5(stdata,state,isterr,err)  ! function of enthalpy and pressure
  wt=stt
endif        
!
wv=stv
wx=stx
wbeta=stbeta
wkapa=stkapa
wcp=stcp
wrho=1/wv
END SUBROUTINE

SUBROUTINE ReadSteamTbl(InDev)
INTEGER :: InDev
INTEGER :: nuse
nuse = mxst
CALL StRead(InDev, nuse, STeamTBL_record, stdata)
CONTINUE
END SUBROUTINE

SUBROUTINE stread(n,nuse,record,a) 
! 
!      stread  - read and initialize steam tables from thermodynamic 
!                properties file data 
!      Calling sequence: 
!                call  stread (ip1,ip2,cp3,rp4)
!      Parameters: 
! 
!                ip1 = n      = FORTRAN unit number from which 
!                               thermodynamic properties file data is 
!                               read (input) 
!                ip2 = nuse   = number of words available in rp4 for 
!                               storage of steam tables (input) 
!                             = number of words of rp4 actually needed f
!                               steam tables (output) 
!                             = -1 if error detected during steam table 
!                                read or initialization (output) 
!                cp3 = record = character array into which information 
!                               about the steam tables and generating 
!                               program is placed (output) 
!                rp4 = a      = array into which steam tables are read 
!                               (output) 
!      I/O units: 
!                ip1 (input);  see above 
!                * (default output) 
!      This routine adapted from sth2xi routine written by R. J. Wagner 
!      for light water steam tables 
REAL(8) a(mxst) 
INTEGER n,nuse 
CHARACTER(80) record(*) 
!     INCLUDE 'stcom.h'  !hj 3apr99
!     common /decart_sth2xc/nt,np,ns, ns2,  klp,klp2, llp,  nt5,   jpl !hj 3apr99
INTEGER i,ios,ntot 
!--rewind thermodynamic properties file 
REWIND n 
!--get thermodynamic properties file title, and information about the 
!--generating program 
READ(n,end=10,err=20,iostat=ios)record(1) 
READ(n,end=10,err=20,iostat=ios)record(2) 
!--get triple point and critical point data, minimum and maximum 
!--temperatures and pressures, table statistics, and table pointers 
READ(n,end=10,err=20,iostat=ios)ttrip,ptrip,vtrip,tcrit,pcrit,                &  
     vcrit,tmin,pmin,tmax,pmax,nt,np,nst,nsp,it3bp,it4bp,it5bp,nprpnt,        &
     it3p0  
!--get number of words in steam tables 
READ(n,end=10,err=20,iostat=ios)ntot 
!--check number of words in steam tables against number of words        
!--available for steam tables storage 
IF(ntot.gt.nuse) GOTO 30
nuse=ntot 
!--get steam tables 
READ(n,end=10,err=20,iostat=ios)(a(i),i=1,ntot) 
GOTO 50 
!--premature end of data encountered 
10 WRITE(*,1001) 
GOTO 40 
!--error reading steam table data  
20 WRITE(*,1002)ios 
GOTO 40 
!--insufficient space 
30 WRITE(*,1003) 
40 nuse=-1 
!--done 
50 RETURN 
1001 FORMAT  ('0***** end of data encountered reading thermodynamic ', 'property file') 
1002 FORMAT  ('0***** read error encountered reading thermodynamic ', 'property file, iostat =',i4) 
1003 FORMAT  ('0***** insufficient space furnished for thermodynamic ', 'property file') 
END subroutine
      
END MODULE
