#include <defines.h>
SUBROUTINE OutputEditMC
USE PARAM
USE TYPEDEF,   ONLY : CoreInfo_Type,    PowerDist_Type
USE GEOM,      ONLY : Core,             ng
USE Core_mod,  ONLY : FmInfo,           CmInfo,          eigv,        &
                      ThInfo,           GroupInfo
USE TH_mod,    ONLY : THVar
USE CNTL,      ONLY : nTracerCntl
USE FILES,     ONLY : IO8
USE PE_mod,    ONLY : PE
USE ioutil,    ONLY : message,          PrintReal1DarrayTo2Darray
!USE VTK_Mod,        ONLY : ProcessVTK, ProcessVTK3D, ProcessVTK3D_HOM
#ifdef MPI_ENV
USE MpiComm_mod, ONLY : BCAST, MPI_SYNC
#endif
IMPLICIT NONE
TYPE(PowerDist_Type) :: PowerDist
INTEGER :: io
LOGICAL :: Master

io = io8
WRITE(io, '(a)') '========================================================================'
WRITE(io, '(5x, a)') 'Results'
WRITE(io, '(a)') '========================================================================'
WRITE(io, *)

END SUBROUTINE