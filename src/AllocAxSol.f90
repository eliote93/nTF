#include <defines.h>
SUBROUTINE AllocAxSol(AxSolver)
USE PARAM
USE TYPEDEF,      ONLY : CoreInfo_TYPE
USE GEOM,         ONLY : Core,         ng
USE CORE_MOD,     ONLY : AxSrc,       AxPXS,          &
                         CmInfo,       FmInfo,    GroupInfo     
USE PE_MOD,       ONLY : PE
USE MPIAxSolver_Mod, ONLY : InitMPIAxNSolver
USE P1SENM_MOD,   ONLY : AllocP1Solver
USE Sp3Senm_Mod,  ONLY : AllocSP3Solver
USE Allocs
IMPLICIT NONE
INTEGER :: AxSolver
INTEGER :: nxy, myzbf, myzef
INTEGER :: mynxybeg, mynxyend

nxy = Core%nxy
myzbf = PE%myzbf; myzef = PE%myzef
!ALLOCATE(AxFlx(myzbf:myzef, nxy))
CALL Dmalloc0(AxSrc, 1, nxy, PE%myzb, PE%myze, 1, ng)
CALL Dmalloc0(AxPxs, 1, nxy, PE%myzb, PE%myze, 1, ng)
CALL InitMPIAxNSolver(Core, GroupInfo, ng, AxSolver, PE)
!#ifdef MPI_ENV
!IF(PE%nCMFDproc .EQ. 1) THEN
!#endif
!  SELECT CASE(AxSolver)
!    CASE(lP1SENM)
!      CALL AllocP1Solver(AxFlx, 1, nxy, myzbf, myzef, ng, TRUE)
!    CASE(lP3SENM)
!      CALL AllocSP3Solver(AxFlx, 1, nxy, myzbf, myzef, ng, TRUE)
!  END SELECT
!#ifdef MPI_ENV
!ELSE
!  !mynxybeg = PE%mynxyBeg; mynxyEnd = PE%mynxyEnd
!  !ALLOCATE(MpiAxFlx(myzbf:myzef, myNxyBeg:myNxyEnd))
!  SELECT CASE(AxSolver)
!    CASE(lP1SENM)
!      CALL AllocP1Solver(AxFlx, 1, nxy, myzbf, myzef, ng, FALSE)
!      !CALL AllocP1Solver(MpiAxFlx, Core, PE, MyNxyBeg, MyNxyEnd, ng, TRUE)
!    CASE(lP3SENM)
!      CALL AllocSP3Solver(AxFlx, 1, nxy, myzbf, myzef, ng, FALSE)
!      !CALL AllocSP3Solver(MpiAxFlx, Core, PE, MyNxyBeg, MyNxyEnd, ng, TRUE)
!  END SELECT
!ENDIF
!CmInfo%MpiAxFlx => MpiAxFlx
!#endif
!CmInfo%AxFlx => AxFlx
CmInfo%AxSrc => AxSrc
CmInfo%AxPXS => AxPXS
FmInfo%AxSrc => AxSrc
FmInfo%AxPXS => AxPXS
END SUBROUTINE