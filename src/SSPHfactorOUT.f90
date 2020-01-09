#include <defines.h>

SUBROUTINE ProcessSSPHfactor(Core,nTracerCntl)
USE PARAM
USE TYPEDEF
USE CNTL,             ONLY : nTracerCntl_Type
USE FILES,            ONLY : io16,                caseid
USE ioutil
USE SPH_mod,          ONLY : ssphf
USE XSLIB_MOD,        ONLY : igresb,igrese
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Pin_Type), POINTER :: PIN(:)

INTEGER :: ig,iz,ipin,nxy,FsrIdxSt,FxrIdxSt,icel,nlocalFxr,j,ifsrst,k
CHARACTER(256) :: fn

IF(.NOT. nTracerCntl%lSSPH) RETURN

Pin => Core%Pin; CellInfo => Core%CellInfo; nxy = Core%nxy

fn=trim(caseid)//'.sphout'
CALL openfile(io16,FALSE,FALSE,FALSE, fn)

WRITE(io16, '(a3, a3, a6, a6, a8)') 'ig','iz','iasy','ipin','ifxrst' 

DO ig=igresb,igrese
    DO iz=1,Core%nz
        DO ipin = 1, nxy
            FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
            icel = Pin(ipin)%Cell(iz)
            nlocalFxr = CellInfo(icel)%nFxr 
            WRITE(io16, '(i3, i3, i6, i6, i8)', ADVANCE='NO') ig,iz,Pin(ipin)%iasy,ipin,FxrIdxSt
            DO j = nLocalFxr,1,-1
                ifsrst = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(1, j) - 1
                WRITE(io16, '(ES18.8)', ADVANCE='NO') ssphf(ifsrst,iz,ig)
            ENDDO !End of Fxr Sweep
            WRITE(io16,*)
        ENDDO !End of Pin Sweep
    ENDDO !End of Plane Sweep
ENDDO

close(io16)

END SUBROUTINE 