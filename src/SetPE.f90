#include <defines.h>
SUBROUTINE MPI_Array_OUT(x, n, myrank, id)
IMPLICIT NONE
REAL :: x(n)
CHARACTER(*) :: id
INTEGER :: n, myrank 
INTEGER :: i, io

io = 100 + myrank
WRITE(io, '(A20)') id
DO i = 1, n
  WRITE(io, '(1pe20.5)') x(i)
ENDDO
END SUBROUTINE
!SUBROUTINE PEInitialize()
!USE PARAM
!#ifdef MPI_ENV
!USE MPI
!#endif
!USE TypeDef,     ONLY : PE_TYPE
!USE PE_Mod,      ONLY : PE
!
!IMPLICIT NONE
!INTEGER :: ierr
!
!#ifdef MPI_ENV
!CALL MPI_INIT(iErr)
!CALL MPI_COMM_RANK(MPI_COMM_WORLD, PE%myrank, ierr)
!CALL MPI_COMM_SIZE(MPI_COMM_WORLD, PE%nproc, ierr)
!PE%MPI_COMM = MPI_COMM_WORLD
!
!PE%master = TRUE
!IF(PE%myrank .NE. 0) PE%master = FALSE
!PE%SLAVE = .NOT. PE%master
!#else
!
!#endif
!ENDSUBROUTINE
!
!SUBROUTINE Pe_Finlize(PE)
!USE PARAM
!USE TypeDef,    ONLY : PE_TYPE
!#ifdef MPI_ENV
!USE MPI
!#endif
!IMPLICIT NONE
!TYPE(PE_TYPE) :: PE
!INTEGER :: ierr
!
!END SUBROUTINE
!
!SUBROUTINE SetPEVariables()
!USE PARAM
!USE TYPEDEF, ONLY : PE_TYPE
!USE PE_Mod,  ONLY : PE
!USE GEOM,    ONLY : CORE, nz, nzfm, nSubPlane
!#ifdef MPI_ENV
!PE%myzb = 1; PE%myze = nz
!PE%myzbf = 1; PE%myzef = nz*nSubPlane
!#else
!PE%myzb = 1; PE%myze = nz
!PE%myzbf = 1; PE%myzef = nz*nSubPlane
!#endif
!END SUBROUTINE

