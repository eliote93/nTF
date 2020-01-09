MODULE B1_Mod
USE PARAM
USE TYPEDEF,        ONLY : XsMac_Type
USE ALLOCS
USE BasicOperation, ONLY : CP_CA
USE LU_MOD,         ONLY : LU_solver, LU_factorize
IMPLICIT NONE
REAL, POINTER, PRIVATE :: D0(:, :), DInv0(:, :), PHI0(:)
INTEGER, PRIVATE :: ng
LOGICAL :: lP1Mod = .FALSE.
LOGICAL :: lAlloc = .FALSE.
LOGICAL :: lTHonly= .TRUE. !14/04/25 for Thermal group only
INTEGER :: bg=35 !boundary group =35 / 1~34 / 35~37
LOGICAL :: lb1k =.FALSE.
REAL :: b1k =1.0
CONTAINS

SUBROUTINE B1Calculation(XsMac, Phi, Bsq, ng0, lP1mod0)
TYPE(XsMac_Type) :: XsMac
REAL :: Phi(ng0), Bsq
REAL :: kinf, keff, kinf_Msq
LOGICAL :: lP1Mod0
INTEGER :: ng0, i 

ng = ng0
lP1mod = lP1mod0

IF(.NOT. lAlloc) THEN
  CALL Dmalloc(D0, ng, ng)
  CALL Dmalloc(DInv0, ng, ng)
  CALL Dmalloc(PHI0, ng)
  lAlloc = .TRUE.
ENDIF

Bsq = 1.0E-8; kinf = 1._8
!lP1mod = .TRUE.
lP1mod = .FALSE. !BYS edit 14/05/01
CALL Set_Dinv(Dinv0, XsMac, Bsq)
CALL Cal_D(D0, Dinv0)
CALL Updt_B1Phi(D0, Bsq, XsMac, Phi0, Kinf)

bsq=0.01_8
IF(.NOT. lP1Mod) THEN
  CALL Set_Dinv(Dinv0, XsMac, Bsq)
  CALL Cal_D(D0, Dinv0)
ENDIF
CALL Updt_B1Phi(D0, Bsq, XsMac, Phi0, Keff)

CALL cal_kinf_Msq(Bsq, keff, kinf, kinf_Msq)
DO i = 1, 100
  CALL search_Bsq(Bsq, keff, kinf_Msq)
  IF(.NOT. lP1Mod) THEN
    CALL Set_Dinv(Dinv0, XsMac, Bsq)
    CALL Cal_D(D0, Dinv0)  
  ENDIF
  CALL Updt_B1Phi(D0, Bsq, XsMac, Phi0, Keff)
  IF(abs(keff-1._8)<1.0E-6) exit 
  CALL cal_kinf_Msq(Bsq, keff, kinf, kinf_Msq)
ENDDO
!WRITE(88, *), i, Bsq, keff
Phi(1:ng) = Phi0(1:ng)
!DO i = 1, ng
!  WRITE(78, *) i, PHI0(i)
!ENDDO
!
!STOP
CONTINUE
END SUBROUTINE


!--- BYSedit
SUBROUTINE B1Calculation_D(XsMac, Phi, Bsq, ng0, lP1mod0, Dng,kinf, keff, kinf_Msq)
    IMPLICIT NONE
    TYPE(XsMac_Type) :: XsMac
    REAL :: Phi(ng0), Bsq
    REAL :: kinf, keff, kinf_Msq
    LOGICAL :: lP1Mod0
    INTEGER :: ng0, i ,g, gg
    REAL :: Dng(ng0)
    ng = ng0
    lP1mod = lP1mod0
    lThonly = .True.
    
    IF(.NOT. lAlloc) THEN
      CALL Dmalloc(D0, ng, ng)
      CALL Dmalloc(DInv0, ng, ng)
      CALL Dmalloc(PHI0, ng)
      lAlloc = .TRUE.
    ENDIF
    
    Bsq = 1.0E-8; kinf = 1._8
    !lP1mod = .TRUE.
    lP1mod = .FALSE. !BYS edit 14/05/01
    CALL Set_Dinv(Dinv0, XsMac, Bsq)
    CALL Cal_D(D0, Dinv0)
    CALL Updt_B1Phi(D0, Bsq, XsMac, Phi0, Kinf)
    
    bsq=0.01_8
    IF(.NOT. lP1Mod) THEN
      CALL Set_Dinv(Dinv0, XsMac, Bsq)
      CALL Cal_D(D0, Dinv0)
    ENDIF
    CALL Updt_B1Phi(D0, Bsq, XsMac, Phi0, Keff)
    
    CALL cal_kinf_Msq(Bsq, keff, kinf, kinf_Msq)
    DO i = 1, 100
      CALL search_Bsq_k(Bsq, keff, kinf_Msq)
      !CALL search_Bsq(Bsq, keff, kinf_Msq)
      IF(.NOT. lP1Mod) THEN
        CALL Set_Dinv(Dinv0, XsMac, Bsq)
        CALL Cal_D(D0, Dinv0)  
      ENDIF
      CALL Updt_B1Phi(D0, Bsq, XsMac, Phi0, Keff)
      IF(abs(keff-b1k)<1.0E-6) exit 
      !IF(abs(keff-1._8)<1.0E-6) exit 
      CALL cal_kinf_Msq(Bsq, keff, kinf, kinf_Msq)
    ENDDO
    !WRITE(88, *), i, Bsq, keff
    Phi(1:ng) = Phi0(1:ng)
    
    Dng=zero
    DO g = 1, ng
        DO gg = 1, ng
            Dng(g)=Dng(g)+D0(g,gg)*phi(gg)
        ENDDO
        Dng(g)=Dng(g)/phi(g)
    ENDDO
    
    
    CONTINUE
END SUBROUTINE
!--- BYSedit end

SUBROUTINE search_Bsq_k(Bsq, keff, kinf_Msq)
    IMPLICIT NONE
    REAL :: Bsq,keff,kinf_Msq,B
    
    Bsq = Bsq + Kinf_Msq * (1._8 - b1k / keff)
    !Bsq = Bsq + Kinf_Msq * (1._8 - 1._8 / keff)
ENDSUBROUTINE

SUBROUTINE search_Bsq(Bsq, keff, kinf_Msq)
    IMPLICIT NONE
    REAL :: Bsq,keff,kinf_Msq,B
    
    Bsq = Bsq + Kinf_Msq * (1._8 - 1._8 / keff)
ENDSUBROUTINE

SUBROUTINE cal_kinf_Msq(Bsq,keff,kinf,kinf_Msq)
    IMPLICIT NONE
    REAL :: Bsq,keff,kinf
    REAL :: kinf_Msq
    kinf_Msq= Bsq / (1._8 / keff - 1._8 / kinf)
END SUBROUTINE


SUBROUTINE Updt_B1Phi(D, Bsq, XS, phi, kinf)
    IMPLICIT NONE
    REAL :: D(ng, ng), Phi(ng)
    REAL :: Bsq, kinf
    TYPE(XsMac_Type) :: Xs
    
    REAL :: A(ng, ng), LU(ng, ng), lmntvec(ng), B(ng)
    INTEGER :: pvec(ng)
    INTEGER :: i,j
    
    DO j = 1, NG
      DO i = 1, NG
        A(i, j) = bsq * D(i, j) - XS%XsMacSm(j, i)
      ENDDO
      A(j, j) = A(j, j) + XS%XsMacT(j)
      b(j) = XS%chi(j)
    ENDDO
    
    CALL LU_factorize(A, LU, pvec, ng)
    CALL LU_solver(LU, pvec, b, phi, ng)
    
    Kinf = 0
    
    DO i = 1, ng
      kinf = kinf + PHI(i) * XS%XsMacNf(i)
    ENDDO
END SUBROUTINE

SUBROUTINE Cal_D(D, Dinv)
    IMPLICIT NONE
    REAL :: D(ng, ng), Dinv(ng, ng)
    REAL :: LU(ng, ng), lmntVec(ng)
    INTEGER :: PVEC(ng)
    INTEGER :: I
    
    LmntVec = 0._8
    CALL LU_FACTORIZE(DInv, LU, PVEC, NG)
    DO I = 1, ng
      LmntVec(i) = 1._8
      D(:, i) = 0._8
      CALL LU_SOLVER(LU, PVEC, lmntVec, D(:, i), ng)
      lmntVec(I) = 0._8
    ENDDO
    
END SUBROUTINE

SUBROUTINE Set_Dinv(Dinv, Xs, Bsq)
    IMPLICIT NONE
    REAL :: Dinv(ng, ng), Bsq
    TYPE(XsMac_Type) :: Xs
    REAL :: alpha
    INTEGER :: i, j
    
    DO j = 1, ng
      alpha = 1._8
      !IF(.NOT. lP1Mod) alpha = falpha(Bsq, Xs%XsMacT(j))
      IF(.NOT. lP1Mod) THEN !---BYS edit 14/04/25
          alpha = falpha(Bsq, Xs%XsMacT(j))
      ENDIF
      DO i = 1, ng
        DInv(i, j) = -3._8 * XS%XsMacP1Sm(j, i)
      ENDDO 
      Dinv(j, j) = Dinv(j, j) + 3._8 * Alpha * XS%XsMacT(j)
    ENDDO
    
END SUBROUTINE

FUNCTION Falpha(Bsq, Sigt)
    IMPLICIT NONE
    REAL :: Bsq,sigt,falpha
    REAL :: a00,B
    
    IF(Bsq>0) THEN
      B=sqrt(Bsq)
      a00 = atan( B / sigt) / B
    ELSE
      B=sqrt(-Bsq)
      a00 = log((sigt+B)/(sigt-B))/(2._8 * B)
    ENDIF
    falpha =a00 * Bsq / sigt / (3._8 * (1._8 - a00 * sigt))
    !---BYS edit 14/04/30 
    ! if alpha == 1 >> B1-P1 hybrid form
    falpha=1._8
    !---BYS edit end
    END FUNCTION

END MODULE
