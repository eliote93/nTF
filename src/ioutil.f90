module ioutil
#ifdef __PGI
#ifdef __WIN64
USE DFPORT
#endif
#endif
#ifdef __INTEL_COMPILER
USE IFPORT
#endif
contains

SUBROUTINE terminate(errmesg)
!
use param
use inputcards
use files
!
character errmesg*(*)
character*6 stars
data stars/'#####'/
print '(a,a)',stars,errmesg
!write(io8,'(a,a)') stars,errmesg
!write(io8,*)
#define DBG
#ifdef DBG
stop '***** Abnormal Termination Due to an Input Error *****'
#endif
stop '***** Abnormal Termination Due to an Input Error *****'
END SUBROUTINE

SUBROUTINE openlf(indev,aline)
use param
use inputcards
implicit none
INTEGER :: indev
character*512  aline
character*256 filename
!
CALL getfn(aline,2,filename)
CALL openfile(indev,TRUE,FALSE, FALSE, filename)
!
return
END SUBROUTINE
!

SUBROUTINE toupper(aa)
use param
! convert lowercase string to uppercase
INTEGER, parameter :: INDXA=97,IDNXZ=122
character aa*(*)
!
lenaa=len_trim(aa)
i=1
DO while (aa(i:i).ne.' ' .and. i.le.lenaa)
   ia=ichar(aa(i:i))
   IF(ia.ge.INDXA) aa(i:i)=char(ia-32)
   i=i+1
   IF(i.gt.lenaa) return
ENDDO
END SUBROUTINE

!
FUNCTION nfields(aline)
!
use param
use inputcards
use cntl
!
!
character*512  aline
logical nonblankd,nonblank,multidata
!
nonblankd=FALSE
multidata=FALSE
!oneline=aline
ncol=len_trim(aline)
n=0
DO i=ncol,1,-1
   IF(aline(i:i).eq.' ' .or. ichar(aline(i:i)).eq.9) then !ichar(tab)=9
      nonblank=.false.
   else
      IF(aline(i:i).eq.'*') then
         multidata=.true.
         multidcol=i
      ENDIF
      nonblank=.true.
      IF(aline(i:i).eq.'!') then
         n=-1
         nonblankd=.true.
      ENDIF
   ENDIF
   IF((.not.nonblankd.and.nonblank) .or. (nonblankd.and..not.nonblank)) then
      n=n+1
      IF(multidata.and.(nonblankd.and. .not.nonblank)) then
         read(aline(i+1:multidcol-1),*) nmult
         n=n+(nmult-1)*2
         multidata=.false.
      ENDIF
   ENDIF
   nonblankd=nonblank
ENDDO
IF(mod(n,2).ne.0) then
  nfields=n/2+1
else
  nfields=n/2
ENDIF
!
return
END FUNCTION

SUBROUTINE expdast(aline, bline)

USE param, ONLY : TRUE, FALSE, BLANK, AST, BANG

IMPLICIT NONE

CHARACTER*512 :: aline, bline, cline, multline

INTEGER :: nhexdat, icol, jcol, kcol, ncol, nmult, nline, ist, jst, jed, lgh1, lgh2, lgh3
CHARACTER*1 :: chr
! ----------------------------------------------------

bline = aline
ncol  = len_trim(bline)
icol = 0

DO
  icol = icol + 1
  
  IF (icol .GT. ncol) EXIT
  
  chr = bline(icol:icol)
  
  IF (chr .EQ. BANG) EXIT
  IF (chr .NE. AST) CYCLE
  
  ! SET : Number of Multiple
  DO jcol = icol-1, 1, -1 ! Non-Blank
    chr = bline(jcol:jcol)
    IF (chr .NE. BLANK) EXIT
  END DO
  
  IF (jcol .EQ. 1) THEN
    kcol = 0
  ELSE
    DO kcol = jcol-1, 1, -1 ! Blank
      chr = bline(kcol:kcol)
      IF (chr .EQ. BLANK) EXIT
    END DO
  END IF
  
  READ (bline(kcol+1:jcol), *) nmult
  ist = kcol+1
  
  ! SET : Multiplier
  DO jcol = icol+1, ncol ! Non-Blank
    chr = bline(jcol:jcol)
    IF (chr .NE. BLANK) EXIT
  END DO
  
  IF (chr .EQ. '(') THEN
    DO kcol = jcol+1, ncol ! End of Paranthesis
      chr = bline(kcol:kcol)
      IF (chr .EQ. ')') EXIT
    END DO
    
    IF (kcol .GT. ncol) CALL terminate("NHEXDAT : PARANTHESIS")
    
    multline = bline(jcol+1:kcol-1)
    cline = bline(kcol+1:ncol) ! Posterior
  ELSE
    DO kcol = jcol+1, ncol ! Blank
      chr = bline(kcol:kcol)
      IF (chr .EQ. BLANK) EXIT
    END DO
    
    multline = bline(jcol:kcol-1)
    cline = bline(kcol:ncol) ! Posterior
  END IF
  
  nline = len_trim(multline)
  
  ! CHK : Size
  lgh1 = ist - 1 ! Anterior
  lgh2 = nmult*(nline + 1) ! Multline
  lgh3 = ncol - kcol + 1 ! Posterior
  IF (lgh1 + lgh2 + lgh3 .GT. 512) CALL terminate("NHEXDAT : LGH")
  
  ! MULT
  bline(ist:) = BLANK
  
  DO jcol = 1, nmult
    jst = ist + (nline+1) * (jcol-1)
    jed = jst + nline-1
    
    bline(jst:jed) = multline
  END DO
  
  bline(jed+1:) = cline ! Posterior
  ncol = len_trim(bline)
END DO
! ----------------------------------------------------

END SUBROUTINE expdast

FUNCTION nhexdat(aline)

USE param, ONLY : BLANK, BANG

IMPLICIT NONE

INTEGER :: nhexdat
CHARACTER*512 :: aline
CHARACTER*1   :: chr

INTEGER :: icol
LOGICAL :: lblkprv, lblknow
! ----------------------------------------------------

lblkprv = aline(1:1) .EQ. BLANK

nhexdat = 0
IF (.NOT.lblkprv) nhexdat = 1

DO icol = 2, len_trim(aline)
  chr = aline(icol:icol)
  
  IF (chr .EQ. BANG) EXIT
  
  lblknow = chr .EQ. BLANK
  
  IF (lblkprv .AND. .NOT.lblknow) nhexdat = nhexdat + 1
  
  lblkprv = lblknow
END DO
! ----------------------------------------------------

END FUNCTION nhexdat

SUBROUTINE multiline(idin,iDOut,ndataf)
use param
use inputcards
use files

nfread=0
nlread=0
DO while(TRUE)
   read(idin,'(a512)',END=100) oneline
   write(iDOut,'(a)') trim(oneline)
   nfread=nfread+nfields(oneline)
   nlread=nlread+1
   IF(nfread.ge.ndataf) exit
ENDDO
DO i=1,nlread
   backspace(idin)
ENDDO
!
return
!
100 continue
CALL terminate('Fewer Data than Required')
!
END SUBROUTINE

!
SUBROUTINE intexpand(n,ival,IField,itarget,ntarget)
!
use param
!
! expand list of INTEGERs containing minus signs and assigns ival
! to target array
dimension IField(*),itarget(*)
iprev=IField(1)
DO i=1,n
   IF(IField(i).gt.0) then
      IF(IField(i).le.ntarget) itarget(IField(i))=ival
      iprev=IField(i)
   elseIF(IField(i).lt.0) then
      DO j=iprev+1,-IField(i)
         IF(j.le.ntarget) itarget(j)=ival
      ENDDO
   else
      mesg='Value 0 Not Allowed'
      CALL terminate(mesg)
   ENDIF
ENDDO
!
return
END SUBROUTINE
!
SUBROUTINE openfile(indev,IFread,IFbin, IFappend, filename)
use param
!
logical IFex,IFread,IFbin, IFappend
character*256 filename
character(256) filename0
INTEGER :: istat

IF(filename(1:1) .eq. '$') THEN
#ifdef __GFORTRAN__
  CALL GET_ENVIRONMENT_VARIABLE(NAME = filename(2:80), VALUE = filename0, STATUS = istat)
  IF(istat .NE. 0) THEN
    mesg='Env. variable does not exist - '//filename
    CALL TERMINATE(mesg)
  ENDIF
#else
  istat = GETENVQQ(filename(2:80), filename0)
  IF(istat .EQ. 0) THEN
    mesg='Env. variable does not exist - '//filename
    CALL TERMINATE(mesg)
  ENDIF
#endif
ELSE
  filename0 = filename
ENDIF
!
IF(IFread) then
   inquire(file=filename0,exist=IFex)
   IF(IFex) then !input file
      IF(IFbin) then
!        IF(IFappend) THEN
          open(indev,file=filename0,status='old',form='unformatted')
!        ELSE
!          open(indev,file=filename,status='old',form='unformatted', ACCESS = 'APPEND')
!        ENDIF
      else
!        IF(IFappend) THEN
!          open(indev,file=filename,status='old', ACCESS = 'APPEND')
!        ELSE
          open(indev,file=filename0,status='old')
!        ENDIF
      ENDIF
   else
      mesg='File DOes not exist - '//filename0
      CALL terminate(mesg)
   ENDIF
else       !open output file
   IF(IFbin) then
      IF(IFappend) THEN
        open(indev,file=filename0,status='old',form='unformatted', ACCESS = 'APPEND')
      ELSE
        open(indev,file=filename0,status='unknown',form='unformatted')
      ENDIF
   else
     IF(IFappend) THEN
       open(indev,file=filename0,status='old', ACCESS = 'APPEND')
     ELSE
       open(indev,file=filename0,status='unknown')
     ENDIF
   ENDIF
ENDIF
!
IF(.NOT. IFappend) rewind(indev)
!
return
END SUBROUTINE

!
SUBROUTINE getfn(aline,i,str)
!
use param
!
! get i-th string on a line of strings
!
character*(*) aline,str
!
ncol=len(aline)
icharjm1=ichar(aline(1:1))  !tab
IF(aline(1:1).eq.' ' .or. icharjm1.eq.9) then
   jf=0
else
   jf=1
ENDIF
jm1=1
DO j=2,ncol
   icharj=ichar(aline(j:j))
   IF((aline(jm1:jm1).eq.' ' .or. icharjm1.eq.9) .and. (aline(j:j).ne.' ' .and. icharj.ne.9)) then
      jf=jf+1
      jbeg=j
   ENDIF
   IF(jf.eq.i) then
      IF(aline(j:j).eq.' ' .or. icharj.eq.9) then
         str=aline(jbeg:j-1)
         return
      elseIF(j.eq.ncol) then
         str=aline(jbeg:j)
         return
      ENDIF
   ENDIF
   jm1=j
   icharjm1=icharj
ENDDO
!
write(mesg,'("File name not present at",i3,"-th field.")') i
CALL terminate(mesg)
!
return
END SUBROUTINE

!
SUBROUTINE cpfile(fn)
!
use param
use inputcards
use files
!
  character*(*) fn
  INTEGER iEND
  probe=BLANK
  indev=io5
  open(io14,file=fn,status='unknown')
  DO while (probe.ne.BANG)
     read(indev,'(a512)',END=100) oneline
     DO iEND=512,1,-1
        IF (oneline(iEND:iEND).ne.BLANK) goto 200
     ENDDO
200        continue
     IF(probe.eq.BANG .or. probe.eq.DOT .or. probe.eq.SLASH) then
        write(io8,'(a)') oneline(1:iEND)
        IF(.not.(probe.eq.DOT .or. probe.eq.SLASH)) then
           read(indev,'(a512)',END=100) oneline
           DO iEND=512,1,-1
              IF (oneline(iEND:iEND).ne.BLANK) goto 210
           ENDDO
210              continue
           write(io8,'(a)') oneline(1:iEND)
        ENDIF
        backspace(indev)
        go to 100
     else
        write(io14,'(a)') oneline(1:iEND)
        write(io8,'(a)') oneline(1:iEND)
     ENDIF
  ENDDO
100     continue
  close(io14)
  return
END SUBROUTINE
!

SUBROUTINE fndchara(aline,ipos,nchar,onec)
character(256) aline
character(1) onec
dimension ipos(*)
nchar=0
nstr=len_trim(aline)
DO i=1,nstr
  IF(aline(i:i).eq.'!') exit
  IF(aline(i:i).eq.onec) then
    nchar=nchar+1
    ipos(nchar)=i
  ENDIF
ENDDO
ipos(nchar+1)=nstr+1
return
END SUBROUTINE
SUBROUTINE fndchara512(aline,ipos,nchar,onec)
character(512) aline
character(1) onec
dimension ipos(*)
nchar=0
nstr=len_trim(aline)
DO i=1,nstr
  IF(aline(i:i).eq.'!') exit
  IF(aline(i:i).eq.onec) then
    nchar=nchar+1
    ipos(nchar)=i
  ENDIF
ENDDO
ipos(nchar+1)=nstr+1
return
END SUBROUTINE
!
SUBROUTINE fndchar(aline,ipos,onec)
character(256) aline
character(1) onec
DO i=1,len_trim(aline)
  IF(aline(i:i).eq.onec) goto 999
ENDDO
999 ipos=i
return
END SUBROUTINE

subroutine charskip(aline1,aline2,onec,nskip)
character(256) aline1,aline2
character(1) onec
integer ipos(10)
call fndchara(aline1,ipos,nchar,onec)
iend=len_trim(aline1)
aline2=aline1(ipos(nskip)+1:iend)
return
end subroutine

SUBROUTINE message(indev,IFcpu,IFdisp,amesg)
!
use param
use timer
use cntl
use files
USE PE_MOD, ONLY : PE
!
character  amesg*(*)
character*1 sc(128)
character*10 form
character*12 atime
INTEGER :: timevals(8),BaseDate(3)
REAL :: firsttm, addtm
INTEGER :: iday
logical IFcpu,IFdisp,IFfdisp
logical :: first=TRUE

data tprinted/0./
data  iday /0/
data addtm /0./
save tprinted,first,firsttm, addtm, iday, BaseDate

REAL :: tlap
INTEGER :: ilo
!
#ifdef __PGI
INTERFACE
FUNCTION setvbuf3f(lu, typ, size)
INTEGER setvbuf3f, lu, typ, size
END FUNCTION
END INTERFACE
ilo = setvbuf3f(indev, 2, 0)
#endif
!
IFfdisp = IFdisp
IF(.NOT. prtscrn) then
  IFfdisp=FALSE
ENDIF

IF(first) then
  CALL date_and_time(values=timevals)
  BaseDate(1:3) = timevals(1:3)
  firsttm=chglong(timevals)
  first=FALSE
ENDIF
!
IF(IFcpu) then
   CALL date_and_time(values=timevals)
   IF(timevals(3) .NE. BaseDate(3)) THEN
      IF(iday == 0) THEN
        addtm = 86400._8 - firsttm
        firsttm = 0
        iday = iday + 1
        BaseDate(1:3) = timevals(1:3)
      ELSE
        addtm = addtm + 86400._8
        firsttm = 0
        iday = iday + 1
        BaseDate(1:3) = timevals(1:3)
      ENDIF
   ENDIF

   tm=chglong(timevals)-firsttm + addtm
   tm=tm-mod(tm,0.001)
   itm=tm
   mss=mod(tm,1.0)*1000
   idd = itm / 86400
   itm = itm - idd * 86400
   ihh=itm/3600
   itm=itm-ihh*3600
   imm=itm/60
   itm=itm-imm*60
   iss=itm
   write(atime,'(2(i2.2,":"),i2.2,".",i3.3)') ihh,imm,iss,mss
   if(idd .GT. 0) write(atime,'((i3,":"),2(i2.2,":"),i2.2)') idd, ihh,imm,iss
   !write(atime,'(2(i2.2,":"),i2.2)') ihh,imm,iss
   IF(IFfdisp) print 600, atime,trim(amesg)
   write(indev,600) atime,trim(amesg)
else
   IF(IFfdisp) print 601,trim(amesg)
   write(indev,601) trim(amesg)
ENDIF
600 format(a,1x,a)
601 format(a)
END SUBROUTINE


FUNCTION IFnumeric(aline)
use param
use inputcards
character aline*(*)
logical :: IFnumeric

IFnumeric=FALSE
oneline=aline
DO i=1,mxncol
   iascii=ichar(sc(i))
   IF(sc(i).ne.BLANK .and. iascii.ne.9 ) then  !determine IF the first character is numeric
      IF((iascii-48)*(iascii-57) .le. 0) IFnumeric=TRUE
      IF(sc(i) .EQ. '+') IFnumeric=TRUE
      IF(sc(i) .EQ. '-') IFnumeric=TRUE
      return
   ENDIF
ENDDO

return
END FUNCTION

FUNCTION IfNumeric_except(aline,iasciiX)
use param
use inputcards
character aline*(*)
logical :: IfNumeric_except
integer :: iasciiX

IfNumeric_except=FALSE
oneline=aline
DO i=1,mxncol
   iascii=ichar(sc(i))
   IF(sc(i).ne.BLANK .and. iascii.ne.9 ) then  !determine IF the first character is numeric
      IF(((iascii-48)*(iascii-57) .le. 0) .or. (iascii.eq.iasciiX)) IfNumeric_except=TRUE
      IF(sc(i) .EQ. '+') IfNumeric_except=TRUE
      IF(sc(i) .EQ. '-') IfNumeric_except=TRUE
      return
   ENDIF
ENDDO

return
END FUNCTION

SUBROUTINE getIFile(inpfile)
character*(*) inpfile
character(256) aarg
INTEGER i2
narg=iargc()
DO i=0,narg
  i2=i
  CALL getarg(i2,aarg)
  lenaa=len_trim(aarg)
  DO j=lenaa,6,-1
    CALL toupper(aarg(j-6:j))
    IF(aarg(j-6:j).eq.'NTRACER') go to 100
  ENDDO
ENDDO
100 continue
IF(i.ge.narg) then
  inpfile=' '
else
  i2=i+1
  CALL getarg(i2,inpfile)
  IF(inpfile(1:1).eq.'-') inpfile=' '
ENDIF
!
return
END SUBROUTINE

SUBROUTINE skiplines(idin,iDOut)
!
use param
!
! skip lines until meeting numeric lines
use inputcards
use files
!logical IFnumeric
!
DO while(TRUE)
   read(idin,'(a512)',END=100) oneline
   IF(IFnumeric(oneline)) then
      backspace(idin)
      return
   ENDIF
   write(iDOut,'(a)') trim(oneline)
ENDDO
100 continue
CALL terminate("Can't Find Numeric Line...")
END SUBROUTINE
!
SUBROUTINE read1more(indev,onelinet,ndataf)
!
use param
use inputcards
use files
!
character*512  onelinet
!character*7 form7

logical IFreadmore
!
IFreadmore=TRUE
ndataf=0
!
DO while (IFreadmore)
   read(indev,'(a512)',END=100) oneline
   IF(oneline.ne.BLANK .and. probe.ne.BANG) then
      ndataf=nfields(oneline)
      IFreadmore=FALSE
   ENDIF
   write(io8,form7(oneline,ncol)) oneline
ENDDO
onelinet=oneline
return
!
100 continue
mesg='END of file reached during read'
CALL terminate(mesg)
!
END SUBROUTINE
!
FUNCTION form7(a,ncol)

use param
! determines the form
use inputcards
!
character*1 a(*)
character*7 form7
DO i=mxncol,1,-1
   IF(a(i).ne.BLANK) go to 50
ENDDO
i=1
50 ncol=i
write(form7,'("(a",i3,")")') ncol
return
END FUNCTION

FUNCTION icolfield(aline,ithfield)
character(256)  aline
logical nonblankd,nonblank,multidata
nonblankd=.false.
multidata=.false.
ncol=len_trim(aline)
n=0
IF(aline(1:1).ne.' ' .and. ichar(aline(1:1)).ne.9) then !ichar(tab)=9
  IField=1
else
  IField=0
ENDIF
DO i=2,256
   IF((aline(i-1:i-1).eq.' ' .or. ichar(aline(i-1:i-1)).eq.9).and. &
     (aline(i:i).ne.' ' .and. ichar(aline(i:i)).ne.9)) then
      IField=IField+1
   ENDIF
   IF(IField.eq.ithfield) exit
ENDDO
icolfield=min(i,256)
return
END FUNCTION

FUNCTION nfieldto(oneline,charto)
character(1) oneline(256),onec,charto
n=1
ipblank=0
IF(oneline(1).eq.' ') ipblank=1
DO i=2,256
  IF(oneline(i).eq.charto) goto 100
  iblank=0
  IF(oneline(i).eq.' ') iblank=1
  IF(ipblank.ne.iblank) n=n+1
  ipblank=iblank
ENDDO
100 continue
nfieldto=n/2
return
END FUNCTION

FUNCTION nfieldto512(oneline,charto)
character(1) oneline(512),onec,charto
n=1
ipblank=0 
IF(oneline(1).eq.' ') ipblank=1 
DO i=2,512
  IF(oneline(i).eq.charto) goto 100
  iblank=0
  IF(oneline(i).eq.' ') iblank=1
  IF(ipblank.ne.iblank) n=n+1
  ipblank=iblank
ENDDO
100 continue
nfieldto512=n/2
return
END FUNCTION

SUBROUTINE PrintReal1DarrayTo2Darray(io, array1D, nxy, nx, ny, FormatData)
USE PARAM
IMPLICIT NONE
REAL :: Array1D(nxy)
INTEGER :: io, nxy, nx, ny
INTEGER :: i, j,jbeg, jend
CHARACTER(120) :: FormatData
DO i = 1, ny
  jbeg = nx * (i - 1) + 1; jend = nx * i
  !WRITE(io, '(7x, 200F8.4)') (Array1D(j), j = jbeg, jend)
  WRITE(io, FormatData) (Array1D(j), j = jbeg, jend)
ENDDO
END SUBROUTINE

SUBROUTINE PrintInt1DarrayTo2Darray(io, array1D, nxy, nx, ny, FormatData)
USE PARAM
IMPLICIT NONE
INTEGER :: Array1D(nxy)
INTEGER :: io, nxy, nx, ny
INTEGER :: i, j,jbeg, jend
CHARACTER(120) :: FormatData
DO i = 1, ny
  jbeg = nx * (i - 1) + 1; jend = nx * i
 !WRITE(io, '(7x, 200I5)') (Array1D(j), j = jbeg, jend)
  WRITE(io, FormatData) (Array1D(j), j = jbeg, jend)
ENDDO
END SUBROUTINE


SUBROUTINE GoToFileEnd(io)
USE PARAM
IMPLICIT NONE
INTEGER :: io
DO
  read(io, *, end=2000)
ENDDO
2000 continue
END SUBROUTINE

SUBROUTINE FnTrim(fn)
CHARACTER(256) :: fn0, fn
INTEGER :: i, j
j=0
fn0=fn
fn=''
DO i = 1, 256
  IF(fn0(i:i) .EQ. '') CYCLE
  IF(fn0(i:i) .EQ. ' ') CYCLE
  IF(fn0(i:i) .EQ. '    ') CYCLE
  j=j+1
  Fn(j:j) = fn0(i:i)
ENDDO

END SUBROUTINE

FUNCTION CreateDir(dir) RESULT(STATUS)
CHARACTER(*) :: dir
LOGICAL :: STATUS
INTEGER :: ISTATUS

STATUS = .FALSE.
INQUIRE(FILE=dir, EXIST=STATUS)

#ifdef __PGI
#ifdef __WIN64
IF(.NOT. STATUS) STATUS = MAKEDIRQQ(dir)
#endif
#ifdef __linux
IF(.NOT. STATUS) THEN
  ISTATUS = system("mkdir " // dir)
  IF(ISTATUS .EQ. 0) STATUS = .TRUE.
ENDIF
#endif
#endif
#ifdef __GFORTRAN__
IF(.NOT. STATUS) THEN
  ISTATUS = system("mkdir " // dir)
  IF(ISTATUS .EQ. 0) STATUS = .TRUE.
ENDIF
#endif
#ifdef __INTEL_COMPILER
IF(.NOT. STATUS) STATUS = MAKEDIRQQ(dir)
#endif

END FUNCTION

FUNCTION ChangeDir(dir) RESULT(STATUS)
CHARACTER(*) :: dir
INTEGER :: STATUS_INT
LOGICAL :: STATUS

STATUS = .FALSE.
STATUS_INT = chdir(dir)
IF (STATUS_INT .EQ. 0) STATUS = .TRUE.

END FUNCTION

FUNCTION PWD() RESULT(dir)
CHARACTER(1024) :: dir
INTEGER :: length

#ifdef __GFORTRAN__
length = GETCWD(dir)
#else
length = GETDRIVEDIRQQ(dir)
#endif

END FUNCTION

FUNCTION GetOutputDir() RESULT(OutputDir)
USE FILES,        ONLY : caseid,      workingdir
CHARACTER(1024) :: outputDir
CHARACTER(1024) :: nowdir
LOGICAL :: STATUS
nowdir = PWD()
STATUS = ChangeDir(workingdir)
STATUS = ChangeDir(caseid)
outputDir = PWD()
STATUS = ChangeDir(nowdir)
END FUNCTION

SUBROUTINE ShowHbar1(io)
USE PARAM
IMPLICIT NONE
INTEGER :: io
WRITE(mesg, '(A)') hbar1(1:77)
CALL message(io, FALSE, TRUE, MESG)
END SUBROUTINE
SUBROUTINE ShowHbar2(io)
USE PARAM
IMPLICIT NONE
INTEGER :: io
WRITE(mesg, '(A)') hbar2(1:77)
CALL message(io, FALSE, TRUE, MESG)
END SUBROUTINE

! This is a simple function to search for an available unit.
! LUN_MIN and LUN_MAX define the range of possible LUNs to check.
! The UNIT value is returned by the function, and also by the optional
! argument. This allows the function to be used directly in an OPEN
! statement, and optionally save the result in a local variable.
! If no units are available, -1 is returned.
integer function newunit(unit)
  implicit none
  integer, intent(out), optional :: unit
! local
  integer, parameter :: LUN_MIN=10, LUN_MAX=1000
  logical :: opened
  integer :: lun
! begin
  newunit=-1
  do lun=LUN_MIN,LUN_MAX
    inquire(unit=lun,opened=opened)
    if (.not. opened) then
      newunit=lun
      exit
    end if
  end do
  if (present(unit)) unit=newunit
end function newunit


END module
