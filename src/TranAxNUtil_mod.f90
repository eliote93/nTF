MODULE TranAxNUtil_Mod
USE PARAM
USE TYPEDEF,        ONLY : CmInfo_Type,          AxFlx_Type,          PE_Type,      &
                           TranCntl_Type


USE CNTL,           ONLY : nTRACERCntl_Type
USE ALLOCS
IMPLICIT NONE
INTEGER, PRIVATE :: AxSolver, ng, nprec
INTEGER, PRIVATE :: mynxybeg, mynxyend, myzb, myze, myzbf, myzef, nz, nzfm

CONTAINS

SUBROUTINE InitTranAxNUtil(AxSolver0, ng0, nprec0, PE)
USE PARAM
IMPLICIT NONE
INTEGER :: AxSolver0, ng0, nprec0
TYPE(PE_TYPE) :: PE
AxSolver = AxSolver0; ng = ng0; nprec = nprec0
mynxybeg = PE%mynxybeg; mynxyend = PE%mynxyend
myzb = PE%myzb; myze = PE%myze
myzbf = PE%myzbf; myzef = PE%myzef
nzfm = PE%nzfm; nz = PE%nz
END SUBROUTINE

SUBROUTINE AllocTranAxNodeVar(Flx, ng0, nprec0)
IMPLICIT NONE
TYPE(AxFlx_Type) :: Flx
INTEGER :: ng0, nprec0
CALL DMALLOC0(Flx%TranPsi, 0, 4)
CALL DMALLOC0(Flx%TranPsid, 0, 4)
CALL DMALLOC0(Flx%Prec, 0, 4, 1, nprec0)
CALL DMALLOC0(Flx%TranSrc, 0, 4, 1, ng0)
END SUBROUTINE

SUBROUTINE InitAxNodePrecursor(Flx, InvLamBeta, reigv)
TYPE(AxFlx_Type) :: Flx
REAL :: InvLamBeta(nprec), reigv
INTEGER :: iprec, i

!reigv = 1._8 / eigv
DO iprec = 1, nprec
  DO i = 0, 4
    Flx%prec(i, iprec) = InvLamBeta(iprec) * reigv * Flx%psi(i)
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE UpdtAxNodePrecursor(Flx, kappa, omegam, omega0, omegap, reigv)
IMPLICIT NONE
TYPE(AxFlx_Type) :: Flx
REAL :: reigv
REAL :: kappa(nprec), omegam(0:nprec), omega0(0:nprec), omegap(0:nprec)
INTEGER :: i, iprec

DO iprec = 1, nprec
  DO i = 0, 4
    Flx%prec(i, iprec) = kappa(iprec) * Flx%prec(i, iprec)
    Flx%prec(i, iprec) = Flx%prec(i, iprec) + omegam(iprec) * Flx%TranPsid(i)
    Flx%prec(i, iprec) = Flx%prec(i, iprec) + omega0(iprec) * Flx%TranPsi(i)
    Flx%prec(i, iprec) = Flx%prec(i, iprec) + omegap(iprec) * Flx%Psi(i) * reigv
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE SaveTranAxNodeSol(Flx, reigv)
!Save Transient solution 
IMPLICIT NONE
TYPE(AxFlx_Type) :: Flx
REAL :: reigv
Flx%TranPsid(0:4) = Flx%TranPsi(0:4)
Flx%TranPsi(0:4) = Flx%Psi(0:4) * reigv

END SUBROUTINE

FUNCTION UpdtAxNodePrecSrc(Flx , lambda, kappa, omegam, omega0)
TYPE(AxFlx_Type) :: Flx
REAL :: lambda(nprec), kappa(nprec), omegam, omega0
REAL :: UpdtAxNodePrecSrc(0:4)
INTEGER :: i, iprec

UpdtAxNodePrecSrc = 0
DO iprec = 1, nprec
  DO i = 0, 4
    UpdtAxNodePrecSrc(i) = UpdtAxNodePrecSrc(i) + lambda(iprec) * kappa(iprec) * Flx%Prec(i, iprec)
  ENDDO
ENDDO

DO i = 0, 4
  UpdtAxNodePrecSrc(i) = UpdtAxNodePrecSrc(i) + Omegam * Flx%TranPsid(i) + Omega0 * Flx%TranPsi(i)
ENDDO
END FUNCTION

SUBROUTINE UpdtAxNodeTranSrc2d(AxFlx, TranSrc2nd, BC)
USE TYPEDEF,      ONLY : AxFlx_TYPE
USE Tlkg_mod,     ONLY : Convtlkg2nd
IMPLICIT NONE
TYPE(AxFlx_Type) :: AxFlx(nzfm)
REAL :: TranSrc2nd(nzfm, ng)
REAL :: Shape2nd(0:2)
INTEGER :: BC(2)

REAL :: NodeAvg(0:nzfm+1)
INTEGER :: iz, ig

DO ig = 1, ng
  NodeAvg = 0
  DO iz = 1, nzfm
    NodeAvg(iz) = TranSrc2nd(iz, ig)
  ENDDO
  IF(BC(1) .EQ. RefCell) NodeAvg(0) = NodeAvg(1)
  IF(BC(2) .EQ. RefCell) NodeAvg(nzfm+1) = NodeAvg(nzfm)

  DO iz = 1, nzfm
    Shape2nd = ConvTlkg2nd(NodeAvg(iz-1:iz+1))
    AxFlx(iz)%TranSrc(0:2, ig) = AxFlx(iz)%TranSrc(0:2, ig) + Shape2nd(0:2)
  ENDDO
ENDDO

END SUBROUTINE

END MODULE
