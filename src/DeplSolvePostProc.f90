#include <Depletion.h>
#include <DefDBG.h>
! Solve System
!   : Predictor/Corrector, Gd QD routines
! Post Process
!   : Burnup for each fxr, Final ND afte PC, Gd Post-correction, [Xe,Sm] Eq/Tr
#if defined __INTEL_MKL || __PGI
#ifdef __INTEL_MKL
SUBROUTINE SolveDeplSys(DeplSysBundle, DeplFxrBundle, lCorrector, lGd, SolverTyp, Nths)
USE HPDeplType
USE CSRMATRIX
USE OMP_LIB
IMPLICIT NONE

INTERFACE
  SUBROUTINE SolveMatExp(DeplMat, pnum0, pnum1, SolverTyp)
  USE CSRMATRIX
  IMPLICIT NONE
  TYPE(CSR_DOUBLE) :: DeplMat
  REAL(8),POINTER :: pnum0(:), pnum1(:)
  INTEGER :: SolverTyp
  END SUBROUTINE SolveMatExp
END INTERFACE

TYPE(DeplSysBundle_Type) :: DeplSysBundle
TYPE(DeplFxrBundle_Type) :: DeplFxrBundle
LOGICAL :: lCorrector, lGd
INTEGER :: Nths, SolverTyp

TYPE(DeplFxr_Type), POINTER :: Fxrs(:)

REAL(8), POINTER :: pnums(:)
REAL(8), POINTER :: pnums1(:)
REAL(8), POINTER :: pnum0(:), pnum1(:)
TYPE(CSR_DOUBLE), POINTER :: DeplMats(:)

INTEGER, POINTER :: IdTrue(:)

INTEGER :: NofIso
INTEGER :: Nsys, ifxrbeg
INTEGER :: i, j, k, l


NofIso = DeplSysBundle%NofIso; Nsys = DeplSysBundle%Nsys; ifxrbeg = DeplSysBundle%ifxrbeg
pnums => DeplSysBundle%pnums; pnums1 => DeplSysBundle%pnums1
DeplMats => DeplSysBundle%DeplMats; Fxrs => DeplFxrBundle%FxrBundle
pnums1 = 0._8
call OMP_SET_NUM_THREADS(Nths);
#ifdef ITER_CHK
call OMP_SET_NUM_THREADS(1);
#endif
!$OMP PARALLEL DO PRIVATE(pnum0, pnum1,j,k,l) SCHEDULE(GUIDED)
DO i = 1, Nsys
  j = (i-1)*NofIso+1; k = i*NofIso; l = ifxrbeg+i-1
  pnum0 => pnums(j:k); pnum1 => pnums1(j:k)
  CALL SolveMatExp(DeplMats(i),pnum0,pnum1,SolverTyp)
  NULLIFY(pnum0,pnum1)
END DO
!$OMP END PARALLEL DO

IF (lGd) THEN
  IdTrue => DeplFxrBundle%IdTrueGd
ELSE
  IdTrue => DeplFxrBundle%IdTrueDepl
END IF
IF (lCorrector) THEN
  !$OMP PARALLEL DO PRIVATE(pnum1,k) SCHEDULE(GUIDED)
  DO i = 1, Nsys
    k = IdTrue(ifxrbeg+i-1)
    IF (Fxrs(k)%lDepl) Fxrs(k)%pnum_cor(1:NofIso) = pnums1((i-1)*NofIso+1:i*NofIso)
  END DO
  !$OMP END PARALLEL DO
ELSE
  !$OMP PARALLEL DO PRIVATE(pnum1,k) SCHEDULE(GUIDED)
  DO i = 1, Nsys
    k = IdTrue(ifxrbeg+i-1)
    IF (Fxrs(k)%lDepl) Fxrs(k)%pnum_pre(1:NofIso) = pnums1((i-1)*NofIso+1:i*NofIso)
  END DO
  !$OMP END PARALLEL DO
END IF
NULLIFY(pnums, pnums1, DeplMats, Fxrs)
END SUBROUTINE

SUBROUTINE SolveMatExp(DeplMat, pnum0, pnum1, SolverTyp)
USE CSRMATRIX
USE MatExponential, ONLY : MatExpKrylov_CSR, MatExpCRAM_CSR, MatExpCRAM_Iter, iLU0, Diag, GS
IMPLICIT NONE
TYPE(CSR_DOUBLE) :: DeplMat
REAL(8),POINTER :: pnum0(:), pnum1(:)
INTEGER :: SolverTyp
TYPE(CSR_DOUBLE) :: Dummy
#ifdef ITER_CHK
  CALL MatExpKrylov_CSR(DeplMat, pnum0, pnum1, 1)
  CALL MatExpCRAM_Iter(.FALSE., DeplMat, Dummy, pnum0, pnum1, GS, numthread=1)
  CALL MatExpCRAM_Iter(.FALSE., DeplMat, Dummy, pnum0, pnum1, Diag, numthread=1)
#endif
IF(SolverTyp .EQ. 2) THEN
  !CALL MatExpCRAM_CSR(.FALSE., DeplMat, Dummy, pnum0, pnum1, numthread=1)
  CALL MatExpCRAM_Iter(.FALSE., DeplMat, Dummy, pnum0, pnum1, iLU0, numthread=1)
ELSE
  CALL MatExpKrylov_CSR(DeplMat, pnum0, pnum1, 1)
END IF

END SUBROUTINE SolveMatExp
#endif

SUBROUTINE UpdatePnumDepl(DeplFxrBundle, lCorrector, lSavePre)
USE HPDeplType
IMPLICIT NONE
TYPE(DeplFxrBundle_Type) :: DeplFxrBundle
LOGICAL :: lCorrector
LOGICAL :: lSavePre

TYPE(DeplFxr_Type), POINTER :: Fxrs(:)

INTEGER :: nTrueDepl
INTEGER :: i, j, k

Fxrs => DeplFxrBundle%FxrBundle; nTrueDepl = DeplFxrBundle%nTrueDepl
IF(lCorrector) THEN
  DO i = 1, nTrueDepl
    k = DeplFxrBundle%IdTrueDepl(i)
    Fxrs(k)%pnum_depl(:) = 0.5*(Fxrs(k)%pnum_pre(:)+Fxrs(k)%pnum_cor(:))
  END DO
  NULLIFY(Fxrs)
ELSEIF(lSavePre) THEN
  DO i = 1, nTrueDepl
    k = DeplFxrBundle%IdTrueDepl(i)
    Fxrs(k)%pnum_depl(:) = Fxrs(k)%pnum_pre(:)
  END DO
  NULLIFY(Fxrs)
END IF
  END SUBROUTINE

SUBROUTINE UpdatePnumDeplGd(DeplFxrBundle, lCorrector)
USE HPDeplType
IMPLICIT NONE
TYPE(DeplFxrBundle_Type) :: DeplFxrBundle
LOGICAL :: lCorrector

TYPE(DeplFxr_Type), POINTER :: aFxr
TYPE(GdFxrInfo_Type), POINTER :: GdFxrs(:)

INTEGER :: nTrueGd
INTEGER :: i, j, k
INTEGER, POINTER :: IdIsoGd(:)

nTrueGd = DeplFxrBundle%nTrueGd
IdIsoGd => DeplFxrBundle%IdIsoGd; GdFxrs => DeplFxrBundle%GdFxrBundle

IF (lCorrector) THEN
  DO i = 1, nTrueGd
    aFxr => GdFxrs(i)%aFxr
    aFxr%pnum_depl(IdIsoGd(:)) = aFxr%pnum_cor(IdIsoGd(:))!*0.3+aFxr%pnum_pre(IdIsoGd(:))*0.7
  END DO
END IF
END SUBROUTINE UpdatePnumDeplGd

SUBROUTINE UpdatePnumSS(DeplFxrBundle, lCorrector)
USE HPDeplType
IMPLICIT NONE
TYPE(DeplFxrBundle_Type) :: DeplFxrBundle
LOGICAL :: lCorrector

TYPE(DeplFxr_Type), POINTER :: Fxrs(:)
INTEGER, POINTER :: IdIsoEig(:)

INTEGER :: nTrueDepl, NisoEig
INTEGER :: i, j, k, l

Fxrs => DeplFxrBundle%FxrBundle; nTrueDepl = DeplFxrBundle%nTrueDepl

DO i = 1, nTrueDepl
  k = DeplFxrBundle%IdTrueDepl(i)
  IdIsoEig => Fxrs(k)%IdIsoEig; NisoEig = Fxrs(k)%NisoEig
  DO j = 1, NisoEig
    l = IdIsoEig(j)
    IF (l.eq.0) CYCLE
    IF(lCorrector) THEN
      Fxrs(k)%pnum_sseig(j) = Fxrs(k)%pnum_depl(l)
    ELSE
      Fxrs(k)%pnum_sseig(j) = Fxrs(k)%pnum_pre(l)
    END IF
  END DO
  NULLIFY(IdIsoEig)
END DO
NULLIFY(Fxrs)
END SUBROUTINE

SUBROUTINE UpdateBU(DeplFxrBundle, lCorrector)
USE HPDeplType
IMPLICIT NONE
TYPE(DeplFxrBundle_Type) :: DeplFxrBundle
LOGICAL :: lCorrector

TYPE(DeplFxr_Type), POINTER :: Fxrs(:)
INTEGER, POINTER :: IdIsoEig(:)

INTEGER :: nTrueDepl, NisoEig
INTEGER :: i, j, k, l
REAL(8) :: power, Hmkg0, delT, bupin
REAL(8), PARAMETER :: barnPcm2 = 1.e+24
REAL(8), PARAMETER :: MWperW = 1.e-6
REAL(8), PARAMETER :: DperS = 1.15740740740740741E-5

Fxrs => DeplFxrBundle%FxrBundle; nTrueDepl = DeplFxrBundle%nTrueDepl
delT = DeplFxrBundle%delT;

DO i = 1, nTrueDepl
  k = DeplFxrBundle%IdTrueDepl(i)
  Hmkg0 = Fxrs(k)%Hmkg0
  IdIsoEig => Fxrs(k)%IdIsoEig
  NisoEig = Fxrs(k)%NisoEig
  power = 0.;
  DO j = 1, NisoEig
    l = IdIsoEig(j)
    power = power + Fxrs(k)%pnum_sseig(j)*Fxrs(k)%Kappa(l)*Fxrs(k)%xs1g(RctIdFis,l)
  END DO
  bupin = Fxrs(k)%burnup0+power/Hmkg0*Fxrs(k)%NormFlux1g*Fxrs(k)%Vol*barnPcm2*MWperW*delT*DperS
  IF (lCorrector) THEN
    Fxrs(k)%burnup = 0.5*(Fxrs(k)%burnup+bupin)
    Fxrs(k)%burnup0 = Fxrs(k)%burnup
  ELSE
    Fxrs(k)%burnup = bupin
  END IF

  NULLIFY(IdIsoEig)
END DO
NULLIFY(Fxrs)
  END SUBROUTINE

SUBROUTINE UpdateEqXe(DeplLib, DeplFxrBundle)
USE HPDeplType
IMPLICIT NONE
TYPE(DeplLib_Type) :: DeplLib
TYPE(DeplFxrBundle_Type) :: DeplFxrBundle

TYPE(DeplFxr_Type), POINTER :: Fxrs(:)
INTEGER, POINTER :: IdIsoEig(:)

INTEGER :: nTrueDepl
INTEGER :: i, j, k, l
INTEGER :: NofIso, NofDec, NofRct, NofFis, NofYld, FisRctId
INTEGER :: NisoEig
REAL(8) :: delT
REAL(8), POINTER :: pnum_sseig(:), pnum_depl(:)

INTEGER, PARAMETER :: IdI135 = 531350, IdXe = 541350
INTEGER :: idepl_I135, idepl_Xe
INTEGER :: iSS_I135, iSS_Xe
INTEGER :: iyld_I135, iyld_Xe
INTEGER :: idecXe, irctXe

REAL(8), POINTER :: DecFrac(:,:), FisFrac(:,:), xs1g(:,:), FisXs(:), Fispnum(:)
INTEGER, POINTER :: IdAftRct(:,:), IdAftDec(:,:)
INTEGER, POINTER :: IdFis(:), IdYld(:)
REAL(8) :: I135Rate, XeRate, FisRate, phi
REAL(8) :: N_I135, N_Xe

INTERFACE
   INTEGER FUNCTION IsoidSrch(IsoArr, NofIso, Idiso)
   USE HPDeplType
   INTEGER, INTENT(IN) :: NofIso, Idiso
   INTEGER, INTENT(IN) :: IsoArr(*)
   END FUNCTION
END INTERFACE

Fxrs => DeplFxrBundle%FxrBundle; nTrueDepl = DeplFxrBundle%nTrueDepl
NofIso = DeplLib%NofIso; NofRct = DeplLib%NofRct;
NofDec = DeplLib%NofDec; NofFis = DeplLib%NofFis;
NofYld = DeplLib%NofYld; FisRctId = DeplLib%FisRctId
!delT = DeplFxrBundle%delT;

ALLOCATE(FisXs(NofFis), Fispnum(NofFis))

DecFrac => DeplLib%DecFrac; FisFrac => DeplLib%FisFrac;
IdAftRct => DeplLib%IdAftRct; IdAftDec => DeplLib%IdAftDec;
IdFis => DeplLib%IdFis; IdYld => DeplLib%IdYld

idepl_I135 = 0; idepl_Xe = 0;

idepl_I135 = IsoidSrch(DeplLib%IdIso, NofIso, IdI135)
idepl_Xe = IsoidSrch(DeplLib%IdIso, NofIso, IdXe)
iyld_I135 = IsoidSrch(IdYld, NofYld, idepl_I135)
iyld_Xe = IsoidSrch(IdYld, NofYld, idepl_Xe)
idecXe = IsoidSrch(IdAftDec(:,idepl_I135), NofDec, idepl_Xe)
irctXe = IsoidSrch(IdAftRct(:,idepl_I135), NofRct, idepl_Xe)

DO i = 1, nTrueDepl
  iSS_I135 = 0; iSS_Xe = 0;
  I135Rate = 0.; XeRate = 0.; FisRate = 0.;
  N_I135 = 0.; N_Xe = 0.;
  k = DeplFxrBundle%IdTrueDepl(i)
  xs1g => Fxrs(k)%xs1g;
  pnum_sseig => Fxrs(k)%pnum_sseig; pnum_depl=> Fxrs(k)%pnum_pre
  phi = Fxrs(k)%NormFlux1g
  FisXs(1:NofFis) = xs1g(FisRctId, IdFis(1:NofFis))
  Fispnum(1:NofFis) = pnum_depl(IdFis(1:NofFis))
  IdIsoEig => Fxrs(k)%IdIsoEig
  NisoEig = Fxrs(k)%NisoEig
  iSS_I135 = IsoidSrch(IdIsoEig, NisoEig, idepl_I135)
  iSS_Xe = IsoidSrch(IdIsoEig, NisoEig, idepl_Xe)

! I-135
  I135Rate = SUM(DecFrac(:,idepl_I135))!*delT
  !I135Rate = I135Rate+SUM(xs1g(:,idepl_I135))*phi!*delT
  FisRate = SUM(FisXs(1:NofFis)*FisFrac(1:NofFis, iyld_I135)*Fispnum(1:NofFis))*phi!*delT
  N_I135 = FisRate/I135Rate;
! Xe-135
  XeRate = SUM(DecFrac(:,idepl_Xe))!*delT
  XeRate = XeRate+SUM(xs1g(:,idepl_Xe))*phi!*delT
  IF(idecXe .GT. 0) I135Rate = DecFrac(idecXe,idepl_Xe)!*delT
  IF(irctXe .GT. 0) I135Rate = I135Rate + xs1g(irctXe,idepl_Xe)*phi!*delT
  FisRate = FisRate+SUM(FisXs(1:NofFis)*FisFrac(1:NofFis, iyld_Xe)*Fispnum(1:NofFis))*phi!*delT
  !N_Xe = (FisRate+I135Rate*N_I135)/XeRate
  N_Xe = (FisRate)/XeRate
! Overwrite
  IF (iSS_I135 .ne. 0) pnum_sseig(iSS_I135) = N_I135;
  IF (iSS_Xe .ne. 0) pnum_sseig(iSS_Xe) = N_Xe;
  !pnum_depl(idepl_I135) = N_I135;
  !pnum_depl(idepl_Xe) = N_Xe;
  NULLIFY(IdIsoEig)
END DO
DEALLOCATE(FisXS, Fispnum)
NULLIFY(Fxrs)
END SUBROUTINE

SUBROUTINE SolveAnalyticGd(DeplLib, DeplFxrBundle, pnumbuf, Nths, nSubStep)
USE HPDeplType
USE OMP_LIB
IMPLICIT NONE
TYPE(DeplLib_Type) :: DeplLib
TYPE(DeplFxrBundle_Type) :: DeplFxrBundle
REAL :: pnumbuf(:)
INTEGER :: Nths, nSubStep

REAL :: delT
REAL :: yld(7), rem(7), expt(7), invDrem(7,6), pnum0(7), bufk, bufl, phi

INTEGER :: ifxr, igd, nfxr
INTEGER :: i, j, k, l
INTEGER, POINTER :: IdIsoGd(:), IdTrueGd(:)

nfxr = DeplFxrBundle%nTrueGd
IdIsoGd => DeplFxrBundle%IdIsoGd
IdTrueGd => DeplFxrBundle%IdTrueGd

delT = DeplFxrBundle%delT / REAL(nSubStep)

call OMP_SET_NUM_THREADS(Nths);

!$OMP PARALLEL DO PRIVATE(igd, yld,rem,expt,invDrem,pnum0,i,j,k,l,bufk,bufl,phi) SCHEDULE (GUIDED)
DO ifxr = 1, nfxr
  igd = IdTrueGd(ifxr)
  pnum0 = pnumbuf(1+(ifxr-1)*7:ifxr*7)
  phi = DeplFxrBundle%FxrBundle(igd)%NormFlux1g
  DO i = 1, 7
    rem(i) = SUM(DeplFxrBundle%FxrBundle(igd)%xs1g(:, IdIsoGd(i)))*phi + &
            SUM(DeplLib%DecFrac(:,IdIsoGd(i)))
    expt(i) = exp(-rem(i)*delT)
  END DO
  DO i = 3, 6
    yld(i) = DeplFxrBundle%FxrBundle(igd)%xs1g(RctIdCAP, IdIsoGd(i-1))*phi
  END DO
  DO i = 2, 5
    DO j = i+1, 6
      invDrem(j,i) = 1./(rem(j)-rem(i))
    END DO
  END DO

  pnumbuf(1+(ifxr-1)*7) = pnum0(1)*expt(1); pnumbuf(7+(ifxr-1)*7) = pnum0(7)*expt(7);
  DO i = 2, 6
    pnumbuf(i+(ifxr-1)*7) = pnum0(i)*expt(i); bufk = 0.;
    DO j = 2, i-1
      bufk = 0.
      DO k = j, i-1
        bufl = yld(k+1)*(expt(k)-expt(i))*invDrem(i,k)
        DO l = j, k-1
          bufl = bufl*yld(l+1)*invDrem(k,l)
        END DO
        DO l = k+1, i-1
          bufl = bufl*yld(l+1)*invDrem(l,k)
        END DO
      END DO
      bufk = bufk+bufl*pnum0(j)
    END DO
    pnumbuf(i+(ifxr-1)*7) = pnumbuf(i+(ifxr-1)*7)+bufk
  END DO
  DO i = 1, 7
    IF (.NOT.(pnumbuf(i+(ifxr-1)*7).GT.0 .OR. pnumbuf(i+(ifxr-1)*7).LE.0)) THEN
      WRITE(*, '(A,I6,A,I6,A)') "NaN during GdQuad at ", ifxr,"th fxr, ",IdIsoGd(i)," isotope"
      WRITE(*, '(A,ES12.5, A, 7ES12.5)') "delT =", delT, " pnum0 :", pnum0(1:7)
      WRITE(*, '(A,7ES12.5)') "Yld :", yld(1:7)
      WRITE(*, '(A,7ES12.5)') "Rem :", rem(1:7)
      STOP
    END IF
  END DO
!  IF (ifxr.EQ.1) THEN
!    WRITE(*,'(3ES12.5)') delT, expt(1), rem(1)
!    WRITE(*,'(A, 2ES12.5)') '152', pnum0(1), pnumbuf(1+(ifxr-1)*7)
!    WRITE(*,'(A, 2ES12.5)') '154',pnum0(2), pnumbuf(2+(ifxr-1)*7)
!    WRITE(*,'(A, 2ES12.5)') '155',pnum0(3), pnumbuf(3+(ifxr-1)*7)
!    WRITE(*,'(A, 2ES12.5)') '156',pnum0(4), pnumbuf(4+(ifxr-1)*7)
!    WRITE(*,'(A, 2ES12.5)') '157',pnum0(5), pnumbuf(5+(ifxr-1)*7)
!    WRITE(*,'(A, 2ES12.5)') '158',pnum0(6), pnumbuf(6+(ifxr-1)*7)
!    WRITE(*,'(A, 2ES12.5)') '160',pnum0(7), pnumbuf(7+(ifxr-1)*7)
!  END IF
END DO
!$OMP END PARALLEL DO


END SUBROUTINE
#endif
