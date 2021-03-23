 module timer
INTEGER, parameter :: size=100
REAL, dimension(:) :: starttime(size)=0
INTEGER :: index = 0
TYPE TimeChk_Type
  REAL :: TotalTime = 0
  REAL :: MocTime = 0
  REAL :: MocRtTime = 0
  REAL :: CmfdTime = 0
  REAL :: CmfdInitTime = 0    !--- CNJ Edit : CPU Overhead Measurement in CUDA CMFD Solver
  REAL :: AxialNodalTime = 0
  REAL :: SubGrpTime = 0
  REAL :: PSMTime = 0
  REAL :: PSMfresoTime = 0
  REAL :: DancoffTime = 0
  REAL :: CPTime = 0
  REAL :: NetRTSubGrpTime = 0
  REAL :: SubGrpGenEFFXSTime = 0
  INTEGER :: SubGrpRTniter = 0
  INTEGER :: SubGrpFSPNum = 0
  REAL :: MocEffTime = 0
  REAL :: AxBTime = 0
  REAL :: AxNSolverTime = 0   !--- CNJ Edit : Comm Overhead Measurement in CUDA Axial Solver
  REAL :: CommTime = 0
  REAL :: DeplTime = 0
  REAL :: ThTime = 0
  REAL :: DeplSetTime = 0, DeplSysTime = 0, DeplSolTime = 0, DeplPostTime = 0.
  REAL :: DeplBase1GTime = 0, DeplBaseMGTime = 0, DeplcuSysTime = 0
  REAL :: XSsubTime = 0!, XSefriTime = 0, XSefmcTime = 0, XSTime = 0
  REAL :: cuXSPreTime = 0, cuXSMainTime = 0, cuHomDevTime = 0, cuHomHostTime = 0
END TYPE
TYPE(TimeChk_Type) :: TimeChk
REAL,save :: tinit,tth,ttotal,txsec,tppr

REAL,SAVE :: sysclock0 = 0, rnr = 0
INTEGER(8), SAVE :: ncmax
contains
SUBROUTINE timeron()
  INTEGER :: timevals(8)
  IF(index.ge.size) return
  index = index+1
!            CALL cpu_time(starttime(index))
  CALL date_and_time(values = timevals)
  starttime(index)=chglong(timevals)
END SUBROUTINE

SUBROUTINE timeroff(totalelapsed)
  INTEGER :: timevals(8)
  REAL ::  elapsed, totalelapsed
  IF(index.le.0) return
!            CALL cpu_time(elapsed)
  CALL date_and_time(values=timevals)
  elapsed=chglong(timevals)
  totalelapsed = totalelapsed+elapsed - starttime(index)
  index = index-1
END SUBROUTINE

SUBROUTINE reset()
starttime=0
index=0
END SUBROUTINE

FUNCTION chglong(timevals)
INTEGER :: timevals(8)
REAL :: chglong
chglong=timevals(5)*3600.0+timevals(6)*60.0+timevals(7)+timevals(8)*0.001
END FUNCTION

function nTracer_dclock(reset,onomp)
USE OMP_LIB
double precision nTracer_dclock !,sysclock0 !,nrn
!common /nTracer_dclock_sysclock/sysclock0,rnr,ncmax
integer(8),save :: ncd=0,nccyc=0
INTEGER(8) :: ncint
logical reset,onomp
REAL :: nc
IF(reset) CALL nTracer_reset_dclock
call system_clock(ncint)
if(ncint.lt.ncd) nccyc=nccyc+1
ncd=ncint
nc=ncint+nccyc*ncmax
nTracer_dclock=nc*rnr-sysclock0

nTracer_dclock=omp_get_wtime()
return
end function

subroutine nTracer_reset_dclock

!double precision sysclock0 !,rnr
INTEGER :: ncint
REAL :: nrint
!common /nTracer_sysclock/sysclock0,rnr,ncmax
call system_clock(ncint,nrint,ncmax)
rnr=nrint
rnr=1/rnr
sysclock0=DBLE(ncint)*rnr
!      write(mesg,'(a,f15.3,a,2(a,f15.3))') " Timer Resolution",rnr*1000000," microsec"
!     +, " Time Left",ncmax*rnr-sysclock0
!      call decart_message(io,.FALSE.,.true.,mesg)
return
end subroutine

END module
