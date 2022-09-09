#include <defines.h>
#include <DefDBG.h>
SUBROUTINE Finalize()
USE PARAM
USE TYPEDEF,       ONLY : PE_TYPE
USE GEOM,          ONLY : CORE
USE CORE_MOD,      ONLY : eigv
USE IOUTIL,        ONLY : message
USE FILES,         ONLY : io5,          io8,        filename,         TempFileIdx, TempInclIdx, lIncFile
USE TIMER,         ONLY : TimeChk
USE itrcntl_mod,   ONLY : ItrCntl
USE CNTL,          ONLY : nTracerCntl
USE PE_MOD,        ONLY : PE
USE SubChCoupling_mod,    ONLY : is_coupled, CodeName
USE ESCOTTh_mod,          ONLY : Finalize_ESCOT
#ifdef MPI_ENV
USE MPIConfig_MOD, ONLY : PeFinalize
USE MPICOMM_MOD,   ONLY : MPI_SYNC,     MPI_MAX_REAL
#endif
#ifdef __PGI
#ifdef __WIN64
USE DFLIB,         ONLY : DELFILESQQ
#endif
#endif
#ifdef __INTEL_COMPILER
USE IFPORT,        ONLY : DELFILESQQ
#endif
#ifdef __INTEL_MKL
USE MKL_3D,         ONLY : mklCntl
#endif
IMPLICIT NONE
CHARACTER(8) :: adate
CHARACTER(10) :: atime
character(2) year,date,hour,min,sec
LOGICAL :: status
#ifdef __PGI
#ifdef __linux
INTEGER :: UNLINK
INTEGER :: istatus
#endif
#endif

IF (is_coupled) THEN
  SELECT CASE (TRIM(CodeName))
  CASE("ESCOT")
    CALL Finalize_ESCOT
  END SELECT
END IF

CALL MPI_MAX_REAL(TimeChk%AxBTime, PE%MPI_NTRACER_COMM, .FALSE.)
CALL MPI_MAX_REAL(TimeChk%CmfdInitTime, PE%MPI_NTRACER_COMM, .FALSE.)
CALL MPI_MAX_REAL(TimeChk%AxNSolverTime, PE%MPI_NTRACER_COMM, .FALSE.)

111 format (A27, 3X, F10.2, A5,A1 A5, I7)
222 format (A27, 3X, F10.2, A5,A1 A5, I7,A1,2x,A10,I4)
333 format (A27, 3X, F10.2, A5,A1 A5, I8,A1,2x,A12,I8)

601 FORMAT (a,a2,':',a2,':',a2,'.')
IF(PE%MASTER) THEN
  Write(mesg, 111)  'Total Time =', TimeChk%TotalTime, 'Sec'
  CALL message(io8, FALSE, TRUE, MESG)
  Write(mesg, 111)  'tst 1 =', TimeChk%tst1, 'Sec' ! DEBUG
  CALL message(io8, FALSE, TRUE, MESG)
  Write(mesg, 111)  'tst 2 =', TimeChk%tst2, 'Sec' ! DEBUG
  CALL message(io8, FALSE, TRUE, MESG)
  IF (nTracerCntl%lED) THEN
    Write(mesg, 333)  'Dancoff. FSP =', TimeChk%DancoffTime, 'Sec', ',', '  N =', nTracerCntl%nCP_er
    CALL message(io8, FALSE, TRUE, MESG)
    Write(mesg, 333)  '1-D CP =', TimeChk%CPTime, 'Sec', ',', '  N =', nTracerCntl%nCP, ',', 'Total # CP =',nTracerCntl%nCP_er+nTracerCntl%nCP
    CALL message(io8, FALSE, TRUE, MESG)
  ELSE
    Write(mesg, 111)  'Subgrp. FSP =', TimeChk%SubGrpTime, 'Sec'
    CALL message(io8, FALSE, TRUE, MESG)
    Write(mesg, 111)  'Subgrp. RT. =', TimeChk%NetRTSubGrpTime, 'Sec', ',', 'N =', TimeChk%SubGrpRTniter
    CALL message(io8, FALSE, TRUE, MESG)
  ENDIF
  Write(mesg, 222)  'MOC MG =', TimeChk%MocTime, 'Sec', ',', 'N =', ItrCntl%MocIt,',', '# THREAD =', PE%nThread
  CALL message(io8, FALSE, TRUE, MESG)
#ifdef __PGI
  Write(mesg, 111)  'MOC RT. =', TimeChk%MocRtTime, 'Sec'
  CALL message(io8, FALSE, TRUE, MESG)
#endif
  Write(mesg, 111)  'CMFD MG =', TimeChk%CmfdTime, 'Sec', ',', 'N =', ItrCntl%Cmfdit
  CALL message(io8, FALSE, TRUE, MESG)
#ifdef __INTEL_MKL
  IF (mklCntl%lGcCMFD) Write(mesg, 111)  'CMFD CG =', 0., 'Sec', ',', 'N =', ItrCntl%GcCMFDIt
  CALL message(io8, FALSE, TRUE, MESG)
#else
  IF (nTracerCntl%lGcCMFD) Write(mesg, 111)  'CMFD CG =', 0., 'Sec', ',', 'N =', ItrCntl%GcCMFDIt
  CALL message(io8, FALSE, TRUE, MESG)
#endif
  IF (PE%lMKL .OR. PE%lCUDACMFD) THEN
    Write(mesg, 111)  'CMFD Init =', TimeChk%CmfdInitTime, 'Sec', ','
    CALL message(io8, FALSE, TRUE, MESG)
  ENDIF
  Write(mesg, 111)  'Ax NODAL =', TimeChk%AxialNodalTime, 'Sec'
  CALL message(io8, FALSE, TRUE, MESG)
  IF (PE%lMKL .OR. PE%lCUDACMFD) THEN
    Write(mesg, 111)  'Ax Kernel =', TimeChk%AxNSolverTime, 'Sec'
    CALL message(io8, FALSE, TRUE, MESG)
  ENDIF
  Write(mesg, 111)  'Ax=b =', TimeChk%AxBTime, 'Sec'
  CALL message(io8, FALSE, TRUE, MESG)
  Write(mesg, 111)  'DEPL =', TimeChk%DeplTime, 'Sec'
  CALL message(io8, FALSE, TRUE, MESG)
  IF (PE%lMKL .OR. PE%lCUDADepl) THEN
#ifdef DEPL_TCHK
    Write(mesg, 111)  'Preset DEPL =', TimeChk%DeplSetTime, 'Sec'
    CALL message(io8, FALSE, TRUE, MESG)
    Write(mesg, 111)  'Sys. Setup DEPL =', TimeChk%DeplSysTime, 'Sec'
    CALL message(io8, FALSE, TRUE, MESG)
    Write(mesg, 111)  'Sol DEPL =', TimeChk%DeplSolTime, 'Sec'
    CALL message(io8, FALSE, TRUE, MESG)
    Write(mesg, 111)  'Post DEPL =', TimeChk%DeplPostTime, 'Sec'
    CALL message(io8, FALSE, TRUE, MESG)
#endif
#ifdef DEPL_SUB_TCHK
    IF (.NOT. PE%lMKL) Write(mesg, 111)  'CUDA Copy DEPL =', TimeChk%DeplcuSysTime, 'Sec'
    IF (.NOT. PE%lMKL) CALL message(io8, FALSE, TRUE, MESG)
#endif
  END IF
#ifdef XS_TCHK
!  Write(mesg, 111)  'XS Gen =', TimeChk%XSTime, 'Sec'
!  CALL message(io8, FALSE, TRUE, MESG)
  Write(mesg, 111)  'Sub Eff Gen =', TimeChk%XSsubTime, 'Sec'
  CALL message(io8, FALSE, TRUE, MESG)
#ifdef XSsub_TCHK
  Write(mesg, 111)  'cuXS Sub Siglp Gen =', TimeChk%cuXSPreTime, 'Sec'
  CALL message(io8, FALSE, TRUE, MESG)
  Write(mesg, 111)  'cuXS Sub Main Part =', TimeChk%cuXSMainTime, 'Sec'
  CALL message(io8, FALSE, TRUE, MESG)
!  Write(mesg, 111)  'RIF Eff Gen =', TimeChk%XSefriTime, 'Sec'
!  CALL message(io8, FALSE, TRUE, MESG)
!  Write(mesg, 111)  'Mac Eff Gen =', TimeChk%XSefmcTime, 'Sec'
!  CALL message(io8, FALSE, TRUE, MESG)
  WRITE(mesg, 111)  'cuXS Hom Dev. =', TimeChk%cuHomDevTime, 'Sec'
  CALL message(io8, FALSE, TRUE, MESG)
  WRITE(mesg, 111)  'cuXS Hom Hst. =', TimeChk%cuHomHostTime, 'Sec'
  CALL message(io8, FALSE, TRUE, MESG)
#endif
#endif
  write(mesg, 111)  'Sub Eff Gen =', TimeChk%SubGrpGenEFFXSTime, 'Sec'
  CALL message(io8, FALSE, TRUE, MESG)
  Write(mesg, 111)  'T-H =', TimeChk%ThTime, 'Sec'
  CALL message(io8, FALSE, TRUE, MESG)
  Write(mesg, 111)  'Comm  =', TimeChk%CommTime, 'Sec'
  CALL message(io8, FALSE, TRUE, MESG)
  !Finalize Message
  CALL message(io8, FALSE, TRUE, '')
  call date_and_time(adate,atime)
  hour=atime(1:2)    !year
  min=atime(3:4)     !min
  sec=atime(5:6)     !sec
  write(mesg,601) "Job Finished at ", hour , min, sec
  CALL message(io8, TRUE, TRUE, mesg)
  WRITE(io8, *)
  WRITE(io8, '(2x, A)') 'Devoloped by'
  WRITE(io8, '(7x, A)') "  ______    __    __   __    __   _______    _______    __        "
  WRITE(io8, '(7x, A)') " /      \  /  \  /  | /  |  /  | /       \  /       \  /  |       "
  WRITE(io8, '(7x, A)') "/$$$$$$  | $$  \ $$ | $$ |  $$ | $$$$$$$  | $$$$$$$  | $$ |       "
  WRITE(io8, '(7x, A)') "$$ \__$$/  $$$  \$$ | $$ |  $$ | $$ |__$$ | $$ |__$$ | $$ |       "
  WRITE(io8, '(7x, A)') "$$      \  $$$$  $$ | $$ |  $$ | $$    $$<  $$    $$/  $$ |       "
  WRITE(io8, '(7x, A)') " $$$$$$  | $$ $$ $$ | $$ |  $$ | $$$$$$$  | $$$$$$$/   $$ |       "
  WRITE(io8, '(7x, A)') "/  \__$$ | $$ |$$$$ | $$ \__$$ | $$ |  $$ | $$ |       $$ |_____  "
  WRITE(io8, '(7x, A)') "$$    $$/  $$ | $$$ | $$    $$/  $$ |  $$ | $$ |       $$       | "
  WRITE(io8, '(7x, A)') " $$$$$$/   $$/   $$/   $$$$$$/   $$/   $$/  $$/        $$$$$$$$/  "
  WRITE(io8, *)
  CLOSE(io8)
ENDIF

IF(PE%MASTER) close(io5)
#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
#endif
IF(PE%MASTER) THEN
  CLOSE(io5)
#ifdef __PGI
#ifdef __linux
  istatus = UNLINK(filename(TempFileIdx))
  IF (.NOT. lIncFile) THEN
    istatus = UNLINK(filename(TempInclIdx))
    istatus = UNLINK('fort.86')
  ENDIF
#endif
#ifdef __WIN64
  status = DELFILESQQ(filename(TempFileIdx))
  IF (.NOT.lIncFile) THEN
    status = DELFILESQQ(filename(TempInclIdx))
    status = DELFILESQQ('fort.86')
  ENDIF
#endif
#endif
#ifdef __GFORTRAN__
  CALL UNLINK(filename(TempFileIdx))
  IF (.NOT. lIncFile) THEN
    CALL UNLINK(filename(TempInclIdx))
    CALL UNLINK('fort.86')
  ENDIF
#endif
#ifdef __INTEL_COMPILER
  status = DELFILESQQ(filename(TempFileIdx))
  IF (.NOT.lIncFile) THEN
    status = DELFILESQQ(filename(TempInclIdx))
    status = DELFILESQQ('fort.86')
  ENDIF
#endif
ENDIF
#ifdef MPI_ENV
CALL PeFinalize()
#endif
END SUBROUTINE
