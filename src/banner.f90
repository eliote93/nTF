SUBROUTINE Banner(OutDev)
USE PARAM
use files,          only : CaseId
#ifdef __PGI
#ifdef __WIN64
USE DFPORT,         only : HOSTNM, GETLOG
#endif
#endif
#ifdef __INTEL_COMPILER
USE IFPORT,         only : HOSTNM
#endif
USE IOUTIL,         only : message
IMPLICIT NONE

INTEGER, intent(in) :: OutDeV
INTEGER :: idatea(3), itimea(3)
INTEGER :: i, j, k
CHARACTER(256) :: BannerData(10)
CHARACTER(80) :: USERNAME, HOSTNAME
CHARACTER(8) :: ADATE
CHARACTER(10) :: ATIME
CHARACTER(2) :: yy,dd,hh,mm,ss
CHARACTER(4) :: mon(12)
CHARACTER(30) :: version
#ifdef __PGI
#ifdef __linux
CHARACTER(80) :: getlog
INTEGER :: hostnm
#endif
#endif

DATA mon /'Jan.','Feb.','Mar.','Apr.','May ','June','July','Aug.','Sep.','Oct.','Nov.','Dec.'/
DATA BannerData &
/"   b.             8 8888888 8888888888 8 888888888o.            .8.           ,o888888o.    8 8888888888   8 888888888o.     ", &
"   888o.          8       8 8888       8 8888    `88.          .888.         8888     `88.  8 8888         8 8888    `88.    ", &
"   Y88888o.       8       8 8888       8 8888     `88         :88888.     ,8 8888       `8. 8 8888         8 8888     `88    ", &
"   .`Y888888o.    8       8 8888       8 8888     ,88        . `88888.    88 8888           8 8888         8 8888     ,88    ", &
"   8o. `Y888888o. 8       8 8888       8 8888.   ,88'       .8. `88888.   88 8888           8 888888888888 8 8888.   ,88'    ", &
"   8`Y8o. `Y88888o8       8 8888       8 888888888P'       .8`8. `88888.  88 8888           8 8888         8 888888888P'     ", &
"   8   `Y8o. `Y8888       8 8888       8 8888`8b          .8' `8. `88888. 88 8888           8 8888         8 8888`8b         ", &
"   8      `Y8o. `Y8       8 8888       8 8888 `8b.       .8'   `8. `88888.`8 8888       .8' 8 8888         8 8888 `8b.       ", &
"   8         `Y8o.`       8 8888       8 8888   `8b.    .888888888. `88888.  8888     ,88'  8 8888         8 8888   `8b.     ", &
"   8            `Yo       8 8888       8 8888     `88. .8'       `8. `88888.  `8888888P'    8 888888888888 8 8888     `88.   "/

!Get HOST Machine Information
#ifdef __PGI
#ifdef __WIN64
! i = hostnm(hostname); CALL getlog(username)
hostname = 'hostname'; username = 'username'
#endif
#ifdef __linux
i = hostnm(hostname); username = getlog()
#endif
#else
i = hostnm(hostname); CALL getlog(username)
#endif

!Get Time
if(len_trim(username).eq.80) username='remote_shell'
call date_and_time(adate,atime)
read(adate(3:4),'(i2)') idatea(3) !year
idatea(3)=idatea(3)+2000
read(adate(5:6),'(i2)') idatea(2) !month
read(adate(7:8),'(i2)') idatea(1) !day
hh=atime(1:2)    !year
mm=atime(3:4)    !min
ss=atime(5:6)    !sec

!PRINT Banner
WRITE(OutDev, '(a126)') 'nTRACERnTRACERnTRACERnTRACERnTRACERnTRACERnTRACERnTRACERnTRACERnTRACERnTRACERnTRACERnTRACERnTRACERnTRACERnTRACERnTRACERnTRACER'
WRITE(OutDev, '(a256)'); WRITE(OutDev, '(a256)'); WRITE(OutDev, '(a256)')
DO i = 1, 10
  WRITE(OutDev, '(x,a256)') BannerData(i)
ENDDO

WRITE(OutDev, '(a256)'); WRITE(OutDev, '(a256)'); WRITE(OutDev, '(a256)')
WRITE(version, '(a)') 'v2.50 Beta - 19/01/07'
WRITE(OutDev, '(56x,A)') TRIM(version)
WRITE(OutDev, '(a256)')
WRITE(OutDev, '(a126)') 'nTRACERnTRACERnTRACERnTRACERnTRACERnTRACERnTRACERnTRACERnTRACERnTRACERnTRACERnTRACERnTRACERnTRACERnTRACERnTRACERnTRACERnTRACER'
WRITE(OutDev, '(a256)');WRITE(OutDev, '(a256)')
!
!Executing Case Print
PRINT '(25x,a8,a,25x)', 'nTRACER ', TRIM(version)
PRINT '(a80)', hbar2
PRINT 600, trim(caseid)
WRITE(OutDev, 600) trim(caseid)
write(mesg,601) "by ",trim(username)," on ",trim(hostname)," at ",hh,mm,ss,mon(idatea(2)),idatea(1),idatea(3),"..."
call message(OutDev,.false.,.TRUE.,mesg)
WRITE(OutDev, '(a256)'); WRITE(OutDev, '(a256)'); WRITE(OutDev, '(a256)')

600 format(13x,"Executing Case ",a)
601 FORMAT (13x,5a,a2,':',a2,':',a2,', ',a4,1x,i2,', ',i4,a)


END SUBROUTINE