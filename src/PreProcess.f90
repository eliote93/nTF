#include <defines.h>
SUBROUTINE PREPROCESS()
USE PARAM
USE FILES,              ONLY : LOCALFN,        FILENAME,          InputFileIdx,           &
                               caseid,         WorkingDir,                                &
                               io8
USE PE_Mod,             ONLY : PE
!USE CNTL, ONLY : MASTER
USE IOUTIL,             ONLY : GETIFILE,       openfile,          PWD,                    &
                               CreateDir,      GetOutputDir,      terminate
#ifdef MPI_ENV
USE MPIComm_Mod,        ONLY : GetMPIFile,		 MPI_SYNC
#endif
USE MCP_Util,           ONLY : lMCP_restart
IMPLICIT NONE

#ifdef MPI_ENV
INCLUDE 'mpif.h'
#endif

INTEGER :: ICHAR1ST

#ifdef __GFORTRAN__
CALL gfortModuleInit
#endif

!Get Working Dir
WorkingDir = PWD()

LOCALFN = ' '
IF(PE%MASTER) CALL GETIFILE(LOCALFN)
#ifdef MPI_ENV
CALL GetMPIFile(LocalFn, PE%MASTER, PE%MPI_COMM)
#endif
ICHAR1ST = ICHAR(LOCALFN(1 : 1))
IF(ICHAR1ST .EQ. 0 .OR. ICHAR1ST .EQ. 32) LOCALFN = 'nTRACER.INP'
FILENAME(InputFileIdx) = LOCALFN
!IF(PE%MASTER) CALL ScanInput_Include
IF (lMCP_restart) CALL ScanInput_Include
#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_COMM)
#endif
CALL ScanInput
!OPEN OUTPut Files
localfn=trim(caseid)//'.out'
IF(PE%master) THEN
  CALL openfile(io8,FALSE,FALSE,FALSE, localfn)
  CALL Banner(io8)
ENDIF
!Create Outuput Dir
!IF(PE%MASTER) THEN
!  mkdir = CreateDir(caseid)
!  IF(.NOT. mkdir) CALL TERMINATE('Fail to Create Output Directory')
!  OutputDir = GetOutputDir()
!ENDIF

END SUBROUTINE

#ifdef __GFORTRAN__
SUBROUTINE gfortModuleInit
USE PARAM
USE TH_Mod,         ONLY : ThOpt
USE itrcntl_mod,    ONLY : ItrCntl
USE CNTL,           ONLY : nTracerCntl

! conductivity, w/m-C, t in K
ThOPT%kFUelCorrelation = (/ 1.05_8, 0._8, 0._8, 0._8, 2150._8, 73.15_8 /)
ThOPT%KCladCorrelation = (/ 7.51_8, 2.09e-2_8, -1.45e-5_8, 7.67e-9_8 /)
!  volumetric heat capacity, in J/m^3-C, t in K
ThOpt%CpFuelCorrelation = (/ 1668768.6_8, 3123.6716_8, -2.4584262_8, 0.000658459_8 /) !rhof=10.282
ThOpt%CpCladCorrelation = (/ 1666764.0_8, 757.284_8, 0., 0. /) !rhoclad=6.6 g/cc

ItrCntl%DcplItrData(0:3, 1) = (/ 1, 3, 2, 2 /)
ItrCntl%DcplItrData(1:2, 2) = (/ 200, 40 /)
ItrCntl%DcplItrData(1:2, 3) = (/ 20, 10 /)
ItrCntl%DcplItrData(1:2, 4) = (/ 2, 2 /)
ItrCntl%DcplItrData(1:3, 5) = (/ 15, 30, 3 /)

nTracerCntl%OutpCntl%IsoOutList(1:18) = (/                                    &
                     90232,      92233,      92234,      92235,      92236,   &
                     92238,      93237,      93239,      94233,      94239,   &
                     94240,      94241,      94242,      95241,      95242,   &
                     95242,      96242,      96244                           /)

nTracerCntl%OutpCntl%IsoBoutp(1:18) = (/                                      &
                     90232,      92233,      92234,      92235,      92236,   &
                     92238,      93237,      93239,      94233,      94239,   &
                     94240,      94241,      94242,      95241,      95242,   &
                     95242,      96242,      96244                           /)

END SUBROUTINE
#endif
