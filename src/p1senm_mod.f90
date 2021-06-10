#include <defines.h>
MODULE P1SENM_MOD
USE TYPEDEF, ONLY : PINXS_TYPE, SCATMAT, AXFLX_TYPE, AxGeom_Type
IMPLICIT NONE
TYPE AxSolItrCntl_TYPE
  INTEGER :: GroupItr = 10
  INTEGER :: NodalItr = 2
END TYPE

REAL,PARAMETER, PRIVATE :: pi=3.1415926535897932384626433832795_8
REAL,PARAMETER, PRIVATE :: tor=1.0E-6
REAL,PARAMETER, PRIVATE :: infinite=1.0e+9
TYPE(AxSolItrCntl_TYPE), PRIVATE :: AxSolItrCntl
INTEGER, PRIVATE :: ng
INTEGER, PRIVATE :: myzbf, myzef, myzb, myze, nz, BC(2)
REAL, POINTER, PRIVATE :: hz(:)
LOGICAL, PRIVATE :: lTran = .FALSE.

REAL, PARAMETER :: kmax = 200._8
CONTAINS

SUBROUTINE SetTranP1SenmEnv(TranCntl)
USE TYPEDEF,   ONLY : TranCntl_Type
IMPLICIT NONE
TYPE(TranCntl_Type) :: TranCntl
lTran = .TRUE.

END SUBROUTINE

SUBROUTINE SetP1SenmEnv(Geom)
TYPE(AxGeom_Type) :: Geom
hz => Geom%H; !COMP => Geom%Comp
BC = Geom%BC; ng = Geom%ng
nz = Geom%nmesh
myzbf = Geom%myzbf; myzef = Geom%myzef
END SUBROUTINE
!SUBROUTINE P1SENM(FLX, XS, eigv, hz, comp, ng0, BC, PE, nz)

SUBROUTINE P1SENM(Flx, XS, eigv, PE)
USE TYPEDEF, ONLY : PE_TYPE
!USE AxSolver_mod, ONLY : AxGeom_TYpe
IMPLICIT NONE
TYPE(AXFLX_TYPE) :: FLX(myzbf:myzef)
TYPE(PINXS_TYPE) :: XS(myzbf:myzef)
!REAL, POINTER :: hz(:)
!INTEGER, POINTER :: comp(:)
TYPE(PE_TYPE) :: PE
REAL :: eigv
!INTEGER :: BC(2), nz
INTEGER :: iz, icL, icR, it

DO it = 1, AxSolItrCntl%NodalItr
  icL =1
  if(myzbf .eq. 1) CALL lbd_P1SENM2n(Flx(1), eigv, xs(icL), Hz(1), BC(1))
  DO iz = myzbf, myzef - 1, 1
    icL = iz; icR = iz + 1
    CALL P1Senm_2n(FLX(iz), FLX(iz + 1), eigv, XS(icL), XS(icR), Hz(iz : iz + 1))
  ENDDO
  icR = myzbf
  IF(myzef .eq. nz) CALL rbd_P1SENM2n(Flx(nz), eigv, xs(icR), Hz(nz), BC(2))
ENDDO

DO iz =  myzbf, myzef
  FLX(iz)%PDHAT = FLX(iz)%DHAT
ENDDO

END SUBROUTINE

SUBROUTINE P1senm_2n(fluxL, fluxR, eigv, sigL, SigR, h)
IMPLICIT NONE
TYPE(AXFLX_TYPE), intent(inout) :: fluxL, fluxR
REAL ::  h(2), eigv
TYPE(PINXS_TYPE) :: SigL,  SigR
!type(cmxs_type), intent(in) :: sig(2)
REAL :: tlkg(0:2, ng, 2)
REAL :: Jn(ng), phis(ng)
REAL :: q(0:4, 2), c(0:4, 2)
REAL :: A(2), B(2), k(2), beta(2), D(2)
REAL :: sinhk(2), coshk(2)
REAL :: BetaL, BetaR, myphi, NeighPhi, Jnet, Jfdm, Dtil, DHat,PDHAT, DEL, ALPHA
INTEGER :: ig, i, j
tlkg(:, :, 1) = FluxL%Tlkg(:, 1, :)
tlkg(:, :, 2) = FluxR%Tlkg(:, 1, :)

DO j = 1, AxSolItrCntl%GroupItr
  fluxL%psi = set_psishape(fluxL%phi(:,  1,  :),  SigL%XSNF)
  fluxR%psi = set_psishape(fluxR%phi(:,  1,  :),  SigR%XSNF)
  DO ig = 1 ,ng
    q(:, 1) = set_qshape(fluxL, tlkg(:, ig, 1), eigv, SigL%chi(ig), SigL%XSS(ig), ig)
    q(:, 2) = set_qshape(fluxR, tlkg(:, ig, 2), eigv, SigR%chi(ig), SigR%XSS(ig), ig)  
    
    k(1) = 0.5_8 * h(1) * sqrt(SigL%XSR(ig)/SigL%XSD(ig)); k(2) = 0.5_8 * h(2) * sqrt(SigR%XSR(ig)/SigR%XSD(ig))
    k(1) = min(k(1), kmax); k(2) = min(k(2), kmax)
    beta(1)=SigL%XSD(ig)/h(1); beta(2)=SigR%XSD(ig)/h(2)
    coshk(1)=cosh(k(1)); coshk(2)=cosh(k(2))
    sinhk(1)=sinh(k(1)); sinhk(2)=sinh(k(2))
    c(:, 1)=set_ptclsol(q(:, 1), k(1), sigL%XSR(ig)); c(:, 2)=set_ptclsol(q(:, 2), k(2), SigR%XSR(ig))
    
    B(1) = k(1) / sinhk(1) * (fluxL%phi(0, 1, ig) - c(0, 1))
    B(2) = k(2) / sinhk(2) * (fluxR%phi(0, 1, ig) - c(0, 2))       

    A = set_Acoeff(sinhk, coshk, B, beta, k, c)
    fluxL%phi(:, 1, ig) = set_phishape(fluxL%phi(0, 1, ig), c(:, 1), A(1), B(1), k(1)) 
    fluxR%phi(:, 1, ig) = set_phishape(fluxR%phi(0, 1, ig), c(:, 2), A(2), B(2), k(2)) 
    
    Jn(ig) = Jright(A(1), B(1), c(:, 1), k(1), beta(1))
    FluxL%JOut(1, 2, ig) = Jn(ig); FluxR%JOut(1, 1, ig) = -Jn(ig)
    !phis(ig)=phisright(A(1), B(1), c(:, 1), k(1), beta(1))    
  ENDDO
ENDDO
!--- Negative Flux Correction ST
!DO ig = 1 ,ng
!  IF(FluxL%phi(0, 1, ig) .LT. 0)THEN
!      FluxL%phi(0, 1, ig)=1E-9
!       WRITE(97,'(a,3i)') 'Negative P1 Correction ', ig
!       WRITE(*,'(a,3i)') 'Negative P1 Correction ', ig
!  ENDIF
!  IF(FluxR%phi(0, 1, ig) .LT. 0)THEN
!      FluxR%phi(0, 1, ig)=1E-9
!       WRITE(97,'(a,3i)') 'Negative P1 Correction ', ig
!       WRITE(*,'(a,3i)') 'Negative P1 Correction ', ig
!  ENDIF
!ENDDO
!Generate Dhat & Dtil
DO ig = 1, ng
  betaL = SigL%XSD(ig)/H(1); betaR = SigR%XSD(ig)/H(2);
  Dtil = 2._8 * BetaL * BetaR / (betaR + betaL)
  myphi = FluxL%phi(0, 1, ig); NeighPhi = FluxR%phi(0, 1, ig)
  JFDM = Dtil*(myPhi - NeighPhi)
  Jnet = FluxL%Jout(1, 2, ig)
  Dhat = (JFDM - Jnet) / (myphi + NeighPhi)
  PDHAT = FluxL%Pdhat(2,IG)
  DEL = ABS(DHAT-PDHAT)
  IF(DEL .GT. 10.*DTIL) THEN
    ALPHA = DTIL/(DEL-DTIL)
    DHAT = PDHAT + ALPHA*(DHAT - PDHAT)
  ENDIF   
  FluxL%Dhat(2, ig) = Dhat; FluxL%Dtil(2, ig) = Dtil
  FluxR%Dhat(1, ig) = -Dhat; FluxR%Dtil(1, ig) = Dtil
ENDDO
CONTINUE
END SUBROUTINE

SUBROUTINE lbd_P1SENM2n(flux, cm_keff, sig, h, ibc)
IMPLICIT NONE
TYPE(AXFLX_TYPE),  INTENT(inout) :: flux
REAL, intent(in) ::  h, cm_keff
INTEGER, intent(in) :: ibc
TYPE(PinXs_Type), intent(in) :: sig
REAL :: tlkg(0:2, ng)
REAL  :: Jn(ng), phis(ng)
REAL  :: q(0:4), c(0:4)
REAL  :: A, B, k, beta, D
REAL  :: albedo(-1:0), alpha
REAL  :: sinhk, coshk
INTEGER :: ig, j
REAL :: myphi, NeighPhi, Jnet, Jfdm, Dtil, DHat

albedo(-1) = 0
albedo(0) = 0.5 !1.E+30
alpha = albedo(ibc)
!call ax_psishape(flux, sig%nu_fs(:))
tlkg(:, :) = Flux%tlkg(:, 1, :)
DO j=1, AxSolItrCntl%GroupItr
  flux%psi = set_psishape(flux%phi(:, 1, :), sig%XSNF)
  DO ig = 1, ng
    !Set Up Source Term
    q = set_qshape(flux, tlkg(:, ig), cm_keff, Sig%chi(ig), Sig%XSS(ig), ig)
    k = 0.5_8 * h * sqrt(sig%XSR(ig)/sig%XSD(ig))
    k = min(k,kmax)
    !k = min(k, 20.0_8)
    sinhk = sinh(k);coshk = cosh(k); beta = sig%XSD(ig)/h
    c(0:4) = set_ptclsol(q(0:4), k, sig%XSR(ig))
    B = k/sinhk * (flux%phi(0, 1, ig) - c(0))
    A = B * (alpha * coshk + 2 * sinhk * beta * k)                      &
        - 2 * beta * (c(1) - 3 * c(2) + 6 * c(3) - 10 * c(4))           &
        + alpha * (c(0) - c(1) + c(2) - c(3) + c(4))
    A = A / (alpha * sinhk + 2 * coshk * beta * k)
    flux%phi(:, 1, ig) = set_phishape(flux%phi(0, 1, ig), c, A, B, k)
    Jn(ig) = Jleft(A, B, c, k, beta)
    Flux%JOUT(1, 1, ig) = -Jn(ig)
!    phis(ig) = phisleft(A, B, c, k, beta)
  ENDDO
ENDDO

!--- Negative Flux Correction ST
!DO ig = 1 ,ng
!  IF(Flux%phi(0, 1, ig) .LT. 0)THEN
!      Flux%phi(0, 1, ig)=1E-9
!       WRITE(97,'(a,3i)') 'Negative P1 Correction LB', ig
!       WRITE(*,'(a,3i)') 'Negative P1 Correction LB', ig
!  ENDIF
!ENDDO
DO ig = 1, ng
  beta = sig%XSD(ig)/h
  Dtil = beta * alpha / (beta + 0.5_8 * alpha)
  myphi = Flux%phi(0, 1, ig); NeighPhi = 0
  JFDM = Dtil * (myphi - NeighPhi); JNet = Flux%JOUT(1, 1, ig)
  Dhat = (JFDM - JNet) / (myphi + NeighPhi)
  Flux%Dhat(1, ig) = Dhat; Flux%Dtil(1, ig) = Dtil
ENDDO
ENDSUBROUTINE

SUBROUTINE rbd_P1SENM2n(flux,cm_keff,sig,h,ibc)
IMPLICIT NONE
type(AXFLX_TYPE),intent(inout) :: flux
REAL,intent(in) :: h, cm_keff
INTEGER, intent(in) :: ibc
type(PinXs_Type), intent(in) :: sig
REAL :: tlkg(0:2, ng)
REAL :: Jn(ng), phis(ng)
REAL:: q(0:4), c(0:4)
REAL:: A, B, k, beta,D
REAL:: albedo(-1:0), alpha
REAL:: sinhk, coshk
REAL :: myphi, NeighPhi, Jnet, Jfdm, Dtil, DHat
INTEGER:: ig,j
albedo(-1)=0
albedo(0)=0.5 !1.E+30

alpha=albedo(ibc)
tlkg(:, :) = Flux%tlkg(:, 1, :)

DO j=1, AxSolItrCntl%GroupItr
  flux%psi = set_psishape(flux%phi, sig%XSNF)
  DO ig = 1, ng
    q = set_qshape(flux, tlkg(:, ig), cm_keff, sig%chi(ig), sig%XSS(ig), ig)
    k = 0.5*h*sqrt(sig%XSR(ig)/sig%XSD(ig))
    k = min(k,kmax)
    !k = min(k, 20.0_8)
    sinhk = sinh(k);coshk = cosh(k)
    beta = sig%XSD(ig)/h
    c(0:4) = set_ptclsol(q(0:4), k, sig%XSR(ig))
    B = k/sinhk*(flux%phi(0, 1, ig)-c(0))
    A = B*(alpha*coshk+2*sinhk*beta*k)+2*beta*(c(1)+3*c(2)+6*c(3)+10*c(4))+alpha*(c(0)+c(1)+c(2)+c(3)+c(4))
    A = -A/(alpha*sinhk+2*coshk*beta*k)
    flux%phi(:,  1,  ig) = set_phishape(flux%phi(0, 1, ig), c, A, B, k)
    Jn(ig) = Jright(A, B, c, k, beta)
    Flux%JOUT(1,  2,  ig)  =  Jn(ig) 
    !phis(ig) = phisright(A, B, c, k, beta)
  ENDDO
  CONTINUE
ENDDO
!DO ig = 1 ,ng
!  IF(Flux%phi(0, 1, ig) .LT. 0)THEN
!      Flux%phi(0, 1, ig)=1E-9
!       WRITE(97,'(a,3i)') 'Negative P1 Correction RB', ig
!       WRITE(*,'(a,3i)') 'Negative P1 Correction RB', ig
!  ENDIF
!ENDDO
DO ig = 1, ng
  beta = sig%XSD(ig)/h
  Dtil = beta * alpha / (beta + 0.5_8 * alpha)
  myphi = Flux%phi(0, 1, ig); NeighPhi = 0
  JFDM = Dtil * (myphi - NeighPhi); JNet = Flux%JOUT(1, 2, ig)
  Dhat = (JFDM - JNet) / (myphi + NeighPhi)
  Flux%Dhat(2, ig) = Dhat; Flux%Dtil(2, ig) = Dtil  
ENDDO
ENDSUBROUTINE


FUNCTION set_psishape(phi,sig)
IMPLICIT NONE
INTEGER :: ig
REAL,INTENT(in) :: phi(0:4, ng), sig(ng)
REAL :: set_psishape(0 : 4), psi(0 : 4)
PSI=0
Set_PsiShape=0
DO ig = 1, ng
  Psi = Psi + Phi(:, ig) * sig(ig)
ENDDO
Set_PsiShape(:) = Psi(:)
END FUNCTION

FUNCTION set_qshape(flux, tlkg, cm_keff, chi, sig, ig)
IMPLICIT NONE
REAL :: set_qshape(0:4)
REAL, INTENT(in) :: tlkg(0:2), cm_keff
TYPE(AXFLX_TYPE), INTENT(inout) :: flux
REAL, INTENT(in) :: chi
TYPE(SCATMAT) :: sig
INTEGER :: ig
INTEGER :: ig0, igb, ige

set_qshape=0
set_qshape=chi*flux%psi/cm_keff

igb = sig%ib; ige = sig%ie
DO ig0= igb, ige
  set_qshape = set_qshape+flux%phi(:, 1, ig0) * sig%from(ig0)
ENDDO

set_qshape(0:2)=set_qshape(0:2)-tlkg(0:2)
IF(lTran) set_qshape(0:4) = set_qshape(0:4) + Flux%TranSrc(0:4, ig)
END FUNCTION

FUNCTION set_phishape(phiavg, c, A, B, k)
IMPLICIT NONE
REAL, intent(in) :: phiavg, c(0:4), A, B, k
REAL :: phi(0:4), set_phishape(0:4)
REAL :: sinhk, coshk, kinv
sinhk = sinh(k)
coshk = cosh(k)
kinv = 1/k
phi(0) = c(0) + sinhk * kinv * B
phi(1) = c(1) + 3._8 * kinv*(coshk - sinhk * kinv) * A
phi(2) = c(2) - 5._8 * (3._8 * coshk * kinv**2-(1._8 + 3._8 * kinv**2) * sinhk * kinv) * B
phi(3) = c(3) + 7._8 * kinv * ((1._8 + 15._8*kinv**2) * coshk - (6._8 + 15._8 * kinv**2) * sinhk * kinv) * A
phi(4) = c(4) - 9._8 * (5._8 * kinv**2 * (2._8 + 21._8 * kinv**2) * coshk-(1._8 + 45._8 * kinv**2 + 105._8 * kinv**4) * sinhk * kinv) * B
set_phishape = phi
END FUNCTION


FUNCTION Jleft(A , B, c , k, beta)
REAL, INTENT(in) :: A, B, k, beta, c(0:4)
REAL :: Jleft
Jleft = 0
Jleft = -2._8 * beta * (cosh(k) * A * k -sinh(k) * B * k + c(1) - 3._8 * c(2) + 6._8 * c(3) - 10._8 * c(4))
END FUNCTION

FUNCTION Jright(A, B, c, k, beta)
implicit none
REAL, intent(in) :: A, B, k, beta, c(0:4)
REAL :: Jright
Jright = 0._8
Jright = -2._8*beta*(cosh(k)*A*k+sinh(k)*B*k+c(1)+3._8*c(2)+6*c(3)+10._8*c(4))
END FUNCTION

FUNCTION phisleft(A, B, c, k, beta)
IMPLICIT NONE
REAL, INTENT(in) :: A, B, k, beta, c(0:4)
REAL :: phisleft
phisleft = -sinh(k)*A+cosh(k)*B+c(0)-c(1)+c(2)-c(3)+c(4)
END FUNCTION

function phisright(A, B, c, k, beta)
IMPLICIT NONE
REAL, INTENT(in) :: A, B, k, beta, c(0:4)
REAL :: phisright
phisright = sinh(k)*A+cosh(k)*B+c(0)+c(1)+c(2)+c(3)+c(4)
END FUNCTION

FUNCTION set_ptclsol(q, k, sigrmv)
IMPLICIT NONE
REAL, INTENT(IN) :: q(0:4), k, sigrmv
REAL :: set_ptclsol(0:4), c(0:4), kinv, siginv

siginv = 1._8 / sigrmv
kinv = 1._8 / k

c(0) = siginv * (q(0) + kinv**2*(3._8*q(2) + 10._8 * q(4)) + kinv**4 * 105._8 * q(4))
c(1) = siginv * (q(1) + 15._8 * kinv**2._8 * q(3))
c(2) = siginv * (q(2) + 35._8 * kinv**2._8 * q(4))
c(3) = siginv * q(3)
c(4) = siginv * q(4)
set_ptclsol = c
END FUNCTION

FUNCTION set_Acoeff(sinhk, coshk, B, beta, k, c)
IMPLICIT NONE
REAL, intent(in) :: sinhk(2), coshk(2), B(2), beta(2), k(2), c(0:4, 2)
REAL :: gamma(2), A(2), set_Acoeff(2), mat(2, 2), detinv, temp
INTEGER :: i

gamma(1) = B(2)*coshk(2)-B(1)*coshk(1)
DO i = 0, 4
  gamma(1) = gamma(1)+(-1)**(i)*c(i, 2)-c(i, 1)
ENDDO
gamma(2) = -2._8*beta(1)*(B(1)*k(1)*sinhk(1)+c(1, 1)+3._8*c(2, 1)+6._8*c(3, 1)+10._8*c(4, 1))
gamma(2) = gamma(2)-2._8*beta(2)*(B(2)*k(2)*sinhk(2)-c(1, 2)+3*c(2, 2)-6._8*c(3, 2)+10._8*c(4, 2))

mat(1, 1) = sinhk(1); mat(1, 2) = sinhk(2)
mat(2, 1) = 2._8*beta(1)*k(1)*coshk(1); mat(2, 2) = -2._8*beta(2)*k(2)*coshk(2)

detinv = 1._8/(mat(1, 1)*mat(2, 2)-mat(1, 2)*mat(2, 1))
temp = mat(1, 1); mat(1, 1) = mat(2, 2)*detinv
mat(2, 2) = temp*detinv
mat(1, 2) = -mat(1, 2)*detinv
mat(2, 1) = -mat(2, 1)*detinv
A(1) = mat(1, 1)*gamma(1)+mat(1, 2)*gamma(2)
A(2) = mat(2, 1)*gamma(1)+mat(2, 2)*gamma(2)
set_Acoeff = A
END FUNCTION

SUBROUTINE GetSubNodePhi(SubPhi, A, B, k, C)
USE PARAM
IMPLICIT NONE
REAL :: SubPhi(-1:1)
REAL :: A, B, k, C(0:4)
REAL :: delx,x1, x2

DelX = 2._8/3._8
x1 = -1._8/3._8; x2 = -1._8; 
SUbPhi(-1) = IntegrateSol(A, B, k, C, x1, x2) / DelX

x1 = 1._8/3._8; x2 = -1._8/3._8; 
SUbPhi(0) = IntegrateSol(A, B, k, C, x1, x2) / DelX

x1 = 1._8; x2 = 1._8/3._8; 
SUbPhi(1) = IntegrateSol(A, B, k, C, x1, x2) / DelX


END SUBROUTINE

FUNCTION IntegrateSol(A, B, k, C, x1, x2)
USE PARAM
IMPLICIT NONE
REAL :: IntegrateSol
REAL :: A, B, C(0:4), k, x1, x2
IntegrateSol = IntegrateSol0(A, B, k, C, x1) - IntegrateSol0(A, B, k, C, x2)
END FUNCTION

FUNCTION IntegrateSol0(A, B, k, C, x)
IMPLICIT NONE
REAL :: IntegrateSol0
REAL :: A, B, C(0:4), k, x
REAL :: y

y = A * Cosh(k * x) / k + B * Sinh(k * x) / k
y = y + x * c(0);
y = y + c(1) * (0.5_8 * x *x)
y = y + c(2) * 0.5_8 * (-x + x**3)
y = y + c(3) * (-0.75_8 * (x**2) + 0.625 * (x**4))
y = y + c(4) * (0.375 * x - 1.25*(x**3) + 0.875 * (x**5))
IntegrateSol0 = y
END FUNCTION

SUBROUTINE AllocP1AxFlxType(AxFlx, ng0, luse)
USE PARAM
USE ALLOCS
IMPLICIT NONE
TYPE(AxFlx_Type) :: AxFlx
INTEGER :: ng0
LOGICAL :: luse       

IF(Luse) THEN
  AxFlx%luse = .TRUE.
  CALL Dmalloc0(AxFlx%Phi, 0, 4, 1, 1, 1, ng0); CALL Dmalloc0(AxFlx%Psi, 0, 4)
  CALL Dmalloc0(AxFlx%Tlkg, 0, 2, 1, 1, 1, ng0); CALL Dmalloc0(AxFlx%Jout, 1, 1, 1, 2, 1, ng0)
  CALL Dmalloc0(AxFlx%LkgShape, 0, 4, 1, ng0); 
  !CALL Dmalloc0(AxFlx%SubNodePhi, -1, 1, 1, ng0); CALL Dmalloc0(AxFlx%SubNodeLkg, -1, 1, 1, ng0)
ENDIF

CALL Dmalloc0(AxFlx%Dhat, 1, 2, 1, ng0); CALL Dmalloc0(AxFlx%Dtil, 1, 2, 1, ng0)
CALL Dmalloc0(AxFlx%PDhat, 1, 2, 1, ng0);
END SUBROUTINE



SUBROUTINE AllocP1Solver(AxFlx, ixybeg, ixyend, izbeg, izend, ng0, lUse)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type, PE_TYPE
IMPLICIT NONE
TYPE(AxFlx_Type), POINTER :: AxFlx(:, :)
!TYPE(PE_TYPE) :: PE
INTEGER :: ixybeg, ixyend, izbeg, izend, ng0
INTEGER :: ixy, iz
LOGICAL :: lUse

!nxy = Core%nxy
!myzbf = PE%myzbf; myzef = PE%myzef
DO ixy = ixybeg, ixyend
  DO iz = izbeg, izend
    CALL AllocP1AxFlxType(AxFlx(iz, ixy), ng0, lUse)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE P1senm_2n_backup(fluxL, fluxR, eigv, sigL, SigR, h)
IMPLICIT NONE
TYPE(AXFLX_TYPE), intent(inout) :: fluxL, fluxR
REAL ::  h(2), eigv
TYPE(PINXS_TYPE) :: SigL,  SigR
!type(cmxs_type), intent(in) :: sig(2)
REAL :: tlkg(0:2, ng, 2)
REAL :: Jn(ng), phis(ng)
REAL :: q(0:4, 2), c(0:4, 2)
REAL :: A(2), B(2), k(2), beta(2), D(2)
REAL :: sinhk(2), coshk(2)
REAL :: BetaL, BetaR, myphi, NeighPhi, Jnet, Jfdm, Dtil, DHat,PDHAT, DEL, ALPHA
INTEGER :: ig, i, j
tlkg(:, :, 1) = FluxL%Tlkg(:, 1, :)
tlkg(:, :, 2) = FluxR%Tlkg(:, 1, :)

DO j = 1, AxSolItrCntl%GroupItr
  fluxL%psi = set_psishape(fluxL%phi(:,  1,  :),  SigL%XSNF)
  fluxR%psi = set_psishape(fluxR%phi(:,  1,  :),  SigR%XSNF)
  DO ig = 1 ,ng
    q(:, 1) = set_qshape(fluxL, tlkg(:, ig, 1), eigv, SigL%chi(ig), SigL%XSS(ig), ig)
    q(:, 2) = set_qshape(fluxR, tlkg(:, ig, 2), eigv, SigR%chi(ig), SigR%XSS(ig), ig)  
    
    k(1) = 0.5_8 * h(1) * sqrt(SigL%XSR(ig)/SigL%XSD(ig)); k(2) = 0.5_8 * h(2) * sqrt(SigR%XSR(ig)/SigR%XSD(ig))
    k(1) = min(k(1), kmax); k(2) = min(k(2), kmax)
    beta(1)=SigL%XSD(ig)/h(1); beta(2)=SigR%XSD(ig)/h(2)
    coshk(1)=cosh(k(1)); coshk(2)=cosh(k(2))
    sinhk(1)=sinh(k(1)); sinhk(2)=sinh(k(2))
    c(:, 1)=set_ptclsol(q(:, 1), k(1), sigL%XSR(ig)); c(:, 2)=set_ptclsol(q(:, 2), k(2), SigR%XSR(ig))
    
    B(1) = k(1) / sinhk(1) * (fluxL%phi(0, 1, ig) - c(0, 1))
    B(2) = k(2) / sinhk(2) * (fluxR%phi(0, 1, ig) - c(0, 2))       

    A = set_Acoeff(sinhk, coshk, B, beta, k, c)
    fluxL%phi(:, 1, ig) = set_phishape(fluxL%phi(0, 1, ig), c(:, 1), A(1), B(1), k(1)) 
    fluxR%phi(:, 1, ig) = set_phishape(fluxR%phi(0, 1, ig), c(:, 2), A(2), B(2), k(2)) 
    
    Jn(ig) = Jright(A(1), B(1), c(:, 1), k(1), beta(1))
    FluxL%JOut(1, 2, ig) = Jn(ig); FluxR%JOut(1, 1, ig) = -Jn(ig)
    !phis(ig)=phisright(A(1), B(1), c(:, 1), k(1), beta(1))    
  ENDDO
ENDDO
ENDSUBROUTINE
!

END MODULE