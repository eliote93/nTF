MODULE DirectSENM_Mod
USE SP3DEF,        ONLY : AxFLX_TYPE,      MODALFLX_TYPE,          &
                          PinXS_TYPE,      AxGEOM_TYPE
USE SENMOP_MOD
USE MAT2x2OP
IMPLICIT NONE
INTEGER,PARAMETER,PRIVATE :: VoidCell = 0, RefCell = -1
CONTAINS

SUBROUTINE SetRHS(B, FLX, MFLX, XS, GEOM, NMESH, IG, NG)
REAL :: B(2, 2*NMESH)
TYPE(AxFLX_TYPE)  :: FLX(NMESH)
TYPE(MODALFLX_TYPE) :: MFLX(NMESH)
TYPE(PinXS_TYPE) :: XS(NMESH)
TYPE(AxGEOM_TYPE) :: GEOM
INTEGER :: NMESH, IG, NG

REAL,DIMENSION(2,NMESH) :: K,KSQ,PHIAVG,SFLX
REAL,DIMENSION(4,NMESH) :: SHK,CHK,SM,KM,MM,QT
REAL :: H(NMESH)

REAL :: TEMP(2), PSOL(0:4, 2, 2)
INTEGER :: I, IROW

H = GEOM%H

DO I=1,NMESH
  KSQ(1, I) = FLX(I)%KSQ(1, IG);             KSQ(2, I) = FLX(I)%KSQ(2, IG)
  K(1,I) = SQRT(FLX(I)%KSQ(1,IG));           K(2,I) = SQRT(FLX(I)%KSQ(2,IG))
  KM(:,I) = KMAT(K(1,I),K(2,I));             SM(:,I) = FLX(I)%S(:,IG)
  MM(:,I) = MMAT(XS(I)%XSD(IG),XS(I)%XSD2(IG),H(I))  
ENDDO

!
IROW = 1
TEMP = SetBoundaryRHS(KM(:,1), SM(:,1), MM(:,1), MFLX(1)%PSOL(:, :, IG), GEOM%BC(1), 1)
B(:, IROW) = TEMP
DO I = 1, NMESH-1
  PSOL(:,:, 1) = MFLX(I)%PSOL(:, :, IG); PSOL(:,:,2) = MFLX(I+1)%PSOL(:, :, IG)
  IROW = IROW + 1; TEMP = SetInterfaceMmtRHS(SM(:,I:I+1), PSOL)
  B(:, IROW) = TEMP  
  IROW = IROW + 1; TEMP = SetInterfaceCurrentRHS(MM(:,I:I+1), SM(:,I:I+1), PSOL)
  B(:, IROW) = TEMP
ENDDO
IROW = IROW + 1;
TEMP = SetBoundaryRHS(KM(:,NMESH), SM(:,NMESH), MM(:,NMESH), MFLX(NMESH)%PSOL(:, :, IG), GEOM%BC(2), 2)
B(:, IROW) = TEMP
END SUBROUTINE

FUNCTION SetBoundaryRHS(KM, SM, MM, C, BC, IDX)
USE  SENMOP_MOD,   ONLY : KMAT
IMPLICIT NONE
REAL :: SetBoundaryRHS(2)
REAL,DIMENSION(4) :: KM,SM,MM
REAL :: C(0:4, 2)
!TYPE(AxFLX_TYPE) :: FLX
!TYPE(ModalFLX_TYPE) :: MFLX
INTEGER :: IDX, BC
REAL :: H
REAL :: ABD(4, -1:0),LP_S(0:4,2),DLP_S(0:4,2)

REAL :: MPHIS(2),MJS(2), PPHIS(2), PJS(2), K(2)

INTEGER :: I

DATA ABD(:,RefCell) / 4*0.0/
DATA ABD(:,VoidCell) / 0.5_8,   0.625_8, -0.125_8, 0.625_8/

DATA LP_S(:,1) / 1._8, -1._8,  1._8,  -1._8,  1._8/
!LEGENDRE POLYNOMIAL AT X= 1
DATA LP_S(:,2) / 1._8,  1._8,  1._8,   1._8,  1._8/
!DERIVATIVES OF LEGENDRE POLYNOMIAL AT X=-1
DATA DLP_S(:,1) / 0._8,  1._8, -3._8,   6._8, -10._8/
!DERIVATIVES OF LEGENDRE POLYNOMIAL AT X=1
DATA DLP_S(:,2) / 0._8,  1._8,  3._8,   6._8,  10._8/

MPHIS = 0; MJS = 0
DO I = 0, 4
  MPHIS = MPHIS + LP_S(I, IDX) * C(I, :)
  MJS = MJS + DLP_S(I, IDX) * C(I, :)
ENDDO
PPHIS = SUBMATVECOP(SM, MPHIS)
PJS = -SUBMATVECOP(SUBMATOP(MM,SM), MJS) 
SetBoundaryRHS = PJS +(-1._8)**(IDX+1) * SUBMATVECOP(ABD(:, BC), PPHIS)
END FUNCTION

FUNCTION SetInterfaceMmtRHS(SM, C)
IMPLICIT NONE
INTEGER, PARAMETER :: L=1,R=2
REAL,DIMENSION(4,2) :: SM
REAL :: C(0:4,2,2)
REAL :: SetInterfaceMmtRHS(2)

REAL :: MPHIS(2,2), PPHIS(2,2), TEMP(2)
MPHIS(:, L) = C(0,:,L) + C(1,:,L) + C(2,:,L) + C(3,:,L) + C(4,:,L)
PPHIS(:, L) = SubMatVecOP(SM(:,L), MPHIS(:,L))
MPHIS(:, R) = C(0,:,R) - C(1,:,R) + C(2,:,R) - C(3,:,R) + C(4,:,R)
PPHIS(:, R) = SubMatVecOP(SM(:,R), MPHIS(:,R))
SetInterfaceMmtRHS = PPHIS(:, L) - PPHIS(:, R)
END FUNCTION

FUNCTION SetInterfaceCurrentRHS(MM, SM, C)
IMPLICIT NONE
INTEGER, PARAMETER :: L=1,R=2
REAL,DIMENSION(4,2) :: SM, MM
REAL :: C(0:4,2,2)
REAL :: SetInterfaceCurrentRHS(2)
REAL :: MJS(2,2), PJS(2,2)

MJS(:,L) = C(1,:,L)+3._8*C(2,:,L)+6._8*C(3,:,L)+10._8*C(4,:,L)
PJS(:,L) = -SUBMATVECOP(SUBMATOP(MM(:,L),SM(:,L)), MJS(:,L))

MJS(:,R) = C(1,:,R)-3._8*C(2,:,R)+6._8*C(3,:,R)-10._8*C(4,:,R)
PJS(:,R) = -SUBMATVECOP(SUBMATOP(MM(:,R),SM(:,R)), MJS(:,R))

SetInterfaceCurrentRHS = PJS(:, L) - PJS(:, R)

END FUNCTION

SUBROUTINE SetSENMMAT(SENMMAT, FLX, XS, GEOM, NMESH, NG)
USE SP3DEF,   ONLY : SENMMAT_TYPE,     PinXS_TYPE,       AxGEOM_TYPE
!USE SenmMatSolver_mod
USE SenmMatSolver_mod, ONLY : SenmMat_LU,   UpdtMatElmt
IMPLICIT NONE
INTEGER :: NMESH, NG
TYPE(SENMMAT_TYPE) :: SENMMAT(NG)
TYPE(AxFlx_TYPE) :: FLX(NMESH)
TYPE(PinXS_TYPE) :: XS(NMESH)
TYPE(AxGEOM_TYPE) :: GEOM

REAL, DIMENSION(4, NG) :: AL, AR, BL, BR
REAL :: BC(2)
INTEGER :: IG, I, ICOL, IROW

BC = GEOM%BC
IROW = 1
CALL BcCondition(AL, BL, FLX(1), XS(1), Geom%H(1), 1, Geom%BC(1), NG)
CALL UpdtMatElmt(SENMMAT, AL, IROW, 1, ng); CALL UpdtMatElmt(SENMMAT, BL, IROW, 2, ng)

DO I = 1, NMESH-1
  IROW = IROW + 1;
  CALL InterfaceMmtCondition(AL, AR, BL, BR, FLX(I:I+1), XS(I:I+1), Geom%H(I:I+1), NG)
  ICOL = 2*I -1;       CALL UpdtMatElmt(SENMMAT, AL, IROW, ICOL, NG)
  ICOL = 2*I;          CALL UpdtMatElmt(SENMMAT, BL, IROW, ICOL, NG)
  ICOL = 2*I +1;       CALL UpdtMatElmt(SENMMAT, AR, IROW, ICOL, NG)
  ICOL = 2*I +2;       CALL UpdtMatElmt(SENMMAT, BR, IROW, ICOL, NG)
  !Apply Interface CUrrent C
  IROW = IROW + 1
  CALL InterfaceCurrentCondition(AL, AR, BL, BR, FLX(I:I+1), XS(I:I+1), Geom%H(I:I+1), NG)
  ICOL = 2*I -1;       CALL UpdtMatElmt(SENMMAT, AL, IROW, ICOL, NG)
  ICOL = 2*I;          CALL UpdtMatElmt(SENMMAT, BL, IROW, ICOL, NG)
  ICOL = 2*I +1;       CALL UpdtMatElmt(SENMMAT, AR, IROW, ICOL, NG)
  ICOL = 2*I +2;       CALL UpdtMatElmt(SENMMAT, BR, IROW, ICOL, NG)
ENDDO
IROW = IROW + 1
!CALL BcCondition(AL, BL, FLX(1), XS(1), Geom%H(1), 1, Geom%BC(1), NG)
CALL BcCondition(AR, BR, FLX(NMESH), XS(NMESH), Geom%H(NMESH), 2, Geom%BC(2), NG)
ICOL = 2 * NMESH - 1;  CALL UpdtMatElmt(SENMMAT, AR, IROW, ICOL, NG)
ICOL = 2 * NMESH;      CALL UpdtMatElmt(SENMMAT, BR, IROW, ICOL, NG)

!LU Factorization
DO IG = 1, NG
  CALL SenmMat_LU(SenmMat(IG), 2*nmesh)
ENDDO
END SUBROUTINE


SUBROUTINE BCCondition(A, B, FLX, XS, H, IDX, BC, NG)
IMPLICIT NONE
TYPE(AxFLX_TYPE) :: FLX
TYPE(PinXS_TYPE) :: XS
REAL :: H
INTEGER :: IDX, BC, NG

REAL :: ABD(4, -1:0)

REAL :: A(4, NG), B(4, NG)
REAL,DIMENSION(2) :: KSQ, K, PHIAVG
REAL,DIMENSION(4),TARGET :: SHK, CHK, SM, KM, MM
REAL,DIMENSION(:),POINTER :: SHKG,CHKG,SMG,KMG,MMG
REAL :: X
INTEGER :: IG

DATA ABD(:,RefCell) /4*0._8/
DATA ABD(:,VoidCell) /0.5_8, 0.625_8, -0.125_8, 0.625_8/

X = (-1._8)**IDX

DO IG = 1, NG
  KSQ(1) = FLX%KSQ(1,IG);  KSQ(2) = FLX%KSQ(2,IG)
  K(1) = SQRT(FLX%KSQ(1,IG));   K(2) = SQRT(FLX%KSQ(2,IG))
  KM(:) = KMAT(K(1),K(2));   SM(:) = FLX%S(:,IG)
  SHK = SINHMAT(K(1),K(2),1._8);   CHK = COSHMAT(K(1),K(2),1._8)
  MM = MMAT(XS%XSD(IG),XS%XSD2(IG),H)
  A(:, IG) = MULTI_MATOP(MM, SM, KM, SHK) + MULTI_MATOP(ABD(:, BC), SM, CHK)
  A(:, IG) = (-1)**IDX * A(:, IG)
  B(:, IG) = MULTI_MATOP(ABD(:, BC), SM, SHK) + MULTI_MATOP(MM, SM, KM, CHK)
ENDDO
END SUBROUTINE

SUBROUTINE InterfaceMmtCondition(A1, A2, B1, B2, FLX, XS, H, NG)
IMPLICIT NONE
REAL, DIMENSION(4, NG) :: A1, A2, B1, B2
TYPE(AxFLX_TYPE) :: FLX(2)
TYPE(PinXS_TYPE) :: XS(2)
REAL :: H(2)
INTEGER :: NG

REAL,DIMENSION(0:4,2,2) :: Q,QH,C,MFLX4TH,FLX4TH
REAL,DIMENSION(2,2) :: K,KSQ,SFLX
REAL,DIMENSION(4,2) ::  SHK, CHK,SM,KM,MM

INTEGER,PARAMETER :: L=1, R=2
INTEGER :: IG, I
!ASSIGN VARIABLE
DO IG = 1, NG
  DO I=1,2
    KSQ(1,I) = FLX(I)%KSQ(1,IG); KSQ(2,I) = FLX(I)%KSQ(2,IG)
    K(1,I) = SQRT(FLX(I)%KSQ(1,IG)); K(2,I) = SQRT(FLX(I)%KSQ(2,IG))
    KM(:,I) = KMAT(K(1,I),K(2,I)); SM(:,I) = FLX(I)%S(:,IG)
    SHK(:,I) = SINHMAT(K(1,I),K(2,I),1._8); CHK(:,I) = COSHMAT(K(1,I),K(2,I),1._8)
    MM(:,I) = MMAT(XS(I)%XSD(IG),XS(I)%XSD2(IG),H(I))
  ENDDO
  A1(:, IG) = 0; A2(:, IG) = 0; 
  B1(:, IG) = 0; B2(:, IG) = 0;
  A1(:, IG) = - SUBMATOP(SM(:, L), CHK(:, L)); A2(:, IG) = SUBMATOP(SM(:, R), CHK(:, R))
  B1(:, IG) = - SUBMATOP(SM(:, L), SHK(:, L)); B2(:, IG) = -SUBMATOP(SM(:, R), SHK(:, R))
ENDDO
END SUBROUTINE

SUBROUTINE InterfaceCurrentCondition(A1, A2, B1, B2, FLX, XS, H, NG)
IMPLICIT NONE
REAL, DIMENSION(4, NG) :: A1, A2, B1, B2
TYPE(AxFLX_TYPE) :: FLX(2)
TYPE(PinXS_TYPE) :: XS(2)
REAL :: H(2)
INTEGER :: NG

REAL,DIMENSION(0:4,2,2) :: Q,QH,C,MFLX4TH,FLX4TH
REAL,DIMENSION(2,2) :: K,KSQ,SFLX
REAL,DIMENSION(4,2) ::  SHK, CHK,SM,KM,MM

INTEGER,PARAMETER :: L=1, R=2
INTEGER :: IG, I
!ASSIGN VARIABLE
DO IG = 1, NG
  DO I=1,2
    KSQ(1,I) = FLX(I)%KSQ(1,IG); KSQ(2,I) = FLX(I)%KSQ(2,IG)
    K(1,I) = SQRT(FLX(I)%KSQ(1,IG)); K(2,I) = SQRT(FLX(I)%KSQ(2,IG))
    KM(:,I) = KMAT(K(1,I),K(2,I)); SM(:,I) = FLX(I)%S(:,IG)
    SHK(:,I) = SINHMAT(K(1,I),K(2,I),1._8); CHK(:,I) = COSHMAT(K(1,I),K(2,I),1._8)
    MM(:,I) = MMAT(XS(I)%XSD(IG),XS(I)%XSD2(IG),H(I))
  ENDDO
  A1(:, IG) = 0; A2(:, IG) = 0; 
  B1(:, IG) = 0; B2(:, IG) = 0;
  !A1(:, IG) = - SUBMATOP(SM(:, L), CHK(:, L)); A2(:, IG) = SUBMATOP(SM(:, R), CHK(:, R))
  !B1(:, IG) = - SUBMATOP(SM(:, L), SHK(:, L)); B2(:, IG) = -SUBMATOP(SM(:, R), SHK(:, R))
  A1(:, IG) =  MULTI_MATOP(MM(:, L), SM(:, L), KM(:, L), SHK(:, L))
  A2(:, IG) =  MULTI_MATOP(MM(:, R), SM(:, R), KM(:, R), SHK(:, R))
  B1(:, IG) =  MULTI_MATOP(MM(:, L), SM(:, L), KM(:, L), CHK(:, L))
  B2(:, IG) = -MULTI_MATOP(MM(:, R), SM(:, R), KM(:, R), CHK(:, R))
ENDDO
END SUBROUTINE


SUBROUTINE UpdtAxialPsi(PSI, Flx, XS, GEOM, NMESH, lPsiLvUpdt)
USE SENMOP_MOD,    ONLY : UPDT_PSI
IMPLICIT NONE
TYPE(AxFLX_TYPE) :: FLX(NMESH)
TYPE(PinXS_TYPE) :: XS(NMESH)
TYPE(AxGEOM_TYPE) :: GEOM
REAL :: PSI(NMESH)
INTEGER :: NMESH, NG
LOGICAL :: lPsiLvUpdt
INTEGER :: I

NG = GEOM%NG
IF(.NOT. lPsiLvUpdt) THEN
  DO I = 1, NMESH
    FLX(I)%PSI = UPDT_PSISHAPE(FLX(I)%PHI(:, 1, :), XS(I), NG)
    PSI(I) = FLX(I)%PSI(0)
  ENDDO
ELSE
  DO I = 1, NMESH
    FLX(I)%PSI = UPDT_PSI(FLX(I)%PHI(:, 1, :), XS(I), NG)
    PSI(I) = FLX(I)%PSI(0)
  ENDDO
ENDIF

END SUBROUTINE

SUBROUTINE UpdtAxialQ(Q, Flx, MFLX, XS, GEOM, KEFF, NMESH, IG, NG, lTran)
USE SENMOP_MOD,    ONLY : UPDT_Q
USE MAT2x2OP
IMPLICIT NONE
TYPE(AxFLX_TYPE) :: FLX(NMESH)
TYPE(MODALFLX_TYPE) :: MFLX(NMESH)
TYPE(PinXS_TYPE) :: XS(NMESH)
TYPE(AxGEOM_TYPE) :: GEOM
LOGICAL :: lTran
INTEGER :: NMESH, IG, NG

REAL :: KEFF
REAL :: Q(0:4, 2, NMESH)

REAL :: QT(4), Q0(0:4, 2)
INTEGER :: I, IOD

  DO I = 1, NMESH
    Q(:, :, I) = UPDT_Q(FLX(I), XS(I), KEFF, IG, NG, lTran)
    QT =  FLX(I)%QT(:,IG)
    DO IOD = 0, 4
      MFLX(I)%QH(IOD, :) = SUBMATVECOP(QT, Q(IOD, :, I))
    ENDDO
  ENDDO

END SUBROUTINE

SUBROUTINE UpdtAxialPSOL(FLX, MFLX, XS, GEOM, NMESH, IG, NG)
USE SENMOP_MOD,    ONLY : PARTICULAR_SOL
IMPLICIT NONE
TYPE(AxFLX_TYPE) :: FLX(NMESH)
TYPE(MODALFLX_TYPE) :: MFLX(NMESH)
TYPE(PinXS_TYPE) :: XS(NMESH)
TYPE(AxGEOM_TYPE) :: GEOM
!REAL :: QH(0:4, 2, NMESH), PSOL(0:4, 2, NMESH)
INTEGER :: NMESH, IG, NG

REAL, DIMENSION(2) :: KSQ
INTEGER :: I, J
DO I = 1, NMESH
  KSQ(1) = FLX(I)%KSQ(1,IG); KSQ(2) = FLX(I)%KSQ(2,IG)
  MFLX(I)%PSOL(:, :, IG) = PARTICULAR_SOL(MFLX(I)%QH(:, :), KSQ)
ENDDO
END SUBROUTINE

SUBROUTINE UpdtFlux(Coeff, FLX, MFLX, XS, GEOM, NMESH, IG, NG)
IMPLICIT NONE
REAL :: COEFF(2, 2*NMESH)
TYPE(AxFlx_TYPE) :: FLX(NMESH)
TYPE(MODALFLX_TYPE) :: MFLX(NMESH)
TYPE(PINXS_TYPE) :: XS(NMESH)
TYPE(AxGEOM_TYPE) :: GEOM
INTEGER :: NMESH, IG, NG
INTEGER :: I, J
REAL :: SM(4)
!REAL :: Vec(2), MAT(4), INVKM(4), KM(4), SHK(4), CHK(4),SM(4)
REAL :: MFLX4TH(0:4,2),FLX4TH(0:4,2)
REAL :: KSQ(2), K(2), A(2), B(2)

DO I = 1, NMESH
  KSQ(1) = FLX(I)%KSQ(1,IG); KSQ(2) = FLX(I)%KSQ(2,IG)
  SM(:) = FLX(I)%S(:, IG)
  MFLX(I)%A(:, IG) = Coeff(:, 2*I-1)
  MFLX(I)%B(:, IG) = Coeff(:, 2*I)
  A(:) = MFLX(I)%A(:, IG); B(:) = MFLX(I)%B(:, IG)
  !4-TH order flux
  MFLX4TH(:, 1) = FLX_4THOD_EXPANSION(A(1), B(1), MFLX(I)%PSOL(:, 1, IG), KSQ(1))
  MFLX4TH(:, 2) = FLX_4THOD_EXPANSION(A(2), B(2), MFLX(I)%PSOL(:, 2, IG), KSQ(2))
  DO J = 0, 4
    FLX4TH(J, :) = SUBMATVECOP(SM, MFLX4TH(J, :))
  ENDDO
  FLX(I)%PHI(:, :, IG)=FLX4TH(:, :)
  !  WRITE(197,'(3i,e)') ig, I, geom%ix, FLX(I)%PHI(0,1,IG)
  !IF(FLX(I)%PHI(0,1,IG) .LT. 0 )THEN
  !     WRITE(198,'(a,3i)') 'Negative AxFLux at-1', ig, I, geom%ix
  !     !WRITE(*,'(a,3i)') 'Negative AxFLux at', ig, I, geom%ix
  !    !FLX(I)%PHI(0,1,IG)=1E-4
  !    !FLX(I)%PHI(1:4,1,IG)=0._8      
  !ENDIf
  !IF(FLX(I)%PHI(0,2,IG) .LT. 0 )THEN
  !     WRITE(198,'(a,3i)') 'Negative AxFLux at-2', ig, I, geom%ix
  !     WRITE(*,'(a,3i)') 'Negative AxFLux at', ig, I, geom%ix
  !    !FLX(I)%PHI(0,2,IG)=1E-4
  !    !FLX(I)%PHI(1:4,2,IG)=0._8      
  !ENDIf
  
  !FLX(I)%PHI(1:4, :, IG)=FLX4TH(1:4, :)
  CONTINUE
ENDDO
END SUBROUTINE

SUBROUTINE UpdtCurrent(FLX, MFLX, XS, GEOM, NMESH, NG)
IMPLICIT NONE
TYPE(AxFlx_TYPE) :: FLX(NMESH)
TYPE(MODALFLX_TYPE) :: MFLX(NMESH)
TYPE(PinXS_TYPE) :: XS(NMESH)
TYPE(AxGEOM_TYPE) :: GEOM
INTEGER :: NMESH, NG

INTEGER :: BC(2)
INTEGER :: I, IG

REAL :: H(NMESH)
REAL :: jnet(2), A(2), B(2), PSOL(0:4, 2)
REAL :: K(2), KSQ(2), KM(4, NMESH, NG), SM(4, NMESH, NG), MM(4, NMESH, NG)
BC = GEOM%BC; H = GEOM%H

DO i = 1, NMESH
  DO ig = 1, NG
    KSQ(1) = FLX(I)%KSQ(1, IG); KSQ(2) = FLX(I)%KSQ(2, IG)
    K(1) = SQRT(KSQ(1)); K(2) = SQRT(KSQ(2))
    KM(:, I, IG) = KMAT(K(1), K(2)); SM(:, I, IG) = FLX(I)%S(:, IG)
    MM(:, I, IG) = MMAT(XS(I)%XSD(IG),XS(I)%XSD2(IG),H(I))
  ENDDO
ENDDO


!Left Boundary
IF(BC(1) .EQ. REFCELL) THEN
  FLX(1)%JOUT(:, 1, :) = 0
ELSE
  DO ig = 1, ng
    !Obtain Current
    A = MFLX(1)%A(:, IG); B = MFLX(1)%B(:, IG)
    PSOL = MFLX(1)%PSOL(:, :, IG)
    JNET = -CAL_J(A, B, PSOL, KM(:, 1, IG), SM(:, 1, IG), MM(:, 1, IG),-1._8)
    FLX(1)%JOUT(:, 1, IG) = JNET    
  ENDDO
ENDIF

DO I = 1, NMESH - 1
  DO ig = 1, ng
    !Obtain Current
    A = MFLX(I)%A(:, IG); B = MFLX(I)%B(:, IG)
    PSOL = MFLX(I)%PSOL(:, :, IG)
    JNET = CAL_J(A, B, PSOL, KM(:, I, IG), SM(:, I, IG), MM(:, I, IG), 1._8)
    FLX(I)%JOUT(:, 2, IG) = JNET; FLX(I+1)%JOUT(:, 1, IG) = -JNET
  ENDDO
ENDDO

!Right Boundary
IF(BC(2) .EQ. REFCELL) THEN
  FLX(NMESH)%JOUT(:, 2, :) = 0
ELSE
  DO ig = 1, ng
    !Obtain Current
    A = MFLX(NMESH)%A(:, IG); B = MFLX(NMESH)%B(:, IG)
    PSOL = MFLX(NMESH)%PSOL(:, :, IG)
    JNET = CAL_J(A, B, PSOL, KM(:, NMESH, IG), SM(:, NMESH, IG), MM(:, NMESH, IG), 1._8)
    FLX(NMESH)%JOUT(:, 2, IG) = JNET    
  ENDDO
ENDIF
CONTINUE

END SUBROUTINE

SUBROUTINE UpdtP1Dhat(FLX, MFLX, XS, GEOM, NMESH, NG)
IMPLICIT NONE
TYPE(AxFlx_TYPE) :: FLX(NMESH)
TYPE(MODALFLX_TYPE) :: MFLX(NMESH)
TYPE(PinXS_TYPE) :: XS(NMESH)
TYPE(AxGEOM_TYPE) :: GEOM
INTEGER :: NMESH, NG

INTEGER :: BC(2)
INTEGER :: I, IG

REAL :: H(NMESH)
REAL :: beta1, beta2, DTIL, DHAT, phi1, phi2, jfdm
REAL :: jnet(2)
REAL :: partialDhat(2)   !--- CNJ Edit : p-CMFD Acceleration
REAL :: PDHAT, DEL, ALPHA, crit
LOGICAL :: lDhatCor
BC = GEOM%BC; H = GEOM%H

lDhatCor=.TRUE.
! lDhatCor=.FALSE.
! ALPHA=0.01_8;
! crit=0.1_8;
crit=10._8;
!Left Boundary
IF(BC(1) .EQ. REFCELL) THEN
  FLX(1)%DTIL(1, :) = 0; FLX(1)%DHAT(1, :) = 0
ELSE
  beta2 = 0.25_8
  DO ig = 1, ng
    !Obtain Current
    JNET = FLX(1)%JOUT(:, 1, IG)
    beta1 = XS(1)%XSD(IG)/H(1); phi1 = FLX(1)%PHI(0, 1, IG)
    DTIL = 2._8 * beta1 * beta2 / (beta1 + beta2)
    JFDM = DTIL*PHI1
    DHAT = (JFDM-JNET(1))/PHI1
    PDHAT = FLX(1)%PDHAT(1, IG)
    DEL = ABS(DHAT-PDHAT)
    IF(DEL .GT. crit*DTIL .AND. lDhatCor) THEN
      ALPHA = DTIL/(DEL-DTIL)
      DHAT = PDHAT + ALPHA*(DHAT - PDHAT)
    ENDIF        
    FLX(1)%DTIL(1, IG) = DTIL; FLX(1)%DHAT(1, IG) = DHAT
  ENDDO
ENDIF

DO I = 1, NMESH - 1
  DO ig = 1, ng
    !Obtain Current
    JNET = FLX(I)%JOUT(:, 2, IG)
    beta1 = XS(I)%XSD(IG)/H(I);     phi1 = FLX(I)%PHI(0, 1, IG)
    beta2 = XS(I+1)%XSD(IG)/H(I+1); phi2 = FLX(I+1)%PHI(0, 1, IG)
    DTIL = 2._8 * beta1 * beta2 / (beta1 + beta2)
    JFDM = DTIL * (PHI1 - PHI2)
    DHAT = (JFDM-JNET(1))/(PHI1+PHI2)
    PDHAT = FLX(I)%PDHAT(2, IG)
    DEL = ABS(DHAT-PDHAT)
    IF(DEL .GT. crit*DTIL .AND. lDhatCor) THEN
      ALPHA = DTIL/(DEL-DTIL)
      DHAT = PDHAT + ALPHA*(DHAT - PDHAT)
    ENDIF        
    FLX(I)%DTIL(2, IG) = DTIL; FLX(I)%DHAT(2, IG) = DHAT
    FLX(I+1)%DTIL(1, IG) = DTIL; FLX(I+1)%DHAT(1, IG) = -DHAT
  ENDDO
ENDDO

!Right Boundary
IF(BC(2) .EQ. REFCELL) THEN
  FLX(NMESH)%DTIL(2, :) = 0; FLX(NMESH)%DHAT(2, :) = 0
ELSE
  beta2 = 0.25_8
  DO ig = 1, ng
    !Obtain Current
    JNET = FLX(NMESH)%JOUT(:, 2, IG)
    beta1 = XS(NMESH)%XSD(IG)/H(NMESH); phi1 = FLX(NMESH)%PHI(0, 1, IG)
    DTIL = 2._8 * beta1 * beta2 / (beta1 + beta2)
    JFDM = DTIL*PHI1
    DHAT = (JFDM-JNET(1))/PHI1
    PDHAT = FLX(NMESH)%PDHAT(2, IG)
    DEL = ABS(DHAT-PDHAT)
    IF(DEL .GT. crit*DTIL .AND. lDhatCor) THEN
      ALPHA = DTIL/(DEL-DTIL)
      DHAT = PDHAT + ALPHA*(DHAT - PDHAT)
    ENDIF        
    FLX(NMESH)%DTIL(2, IG) = DTIL; FLX(NMESH)%DHAT(2, IG) = DHAT
  ENDDO
ENDIF

END SUBROUTINE

SUBROUTINE PrintFlux(FLX, MFLX, GEOM, NMESH, IG1, IG2, NG)
IMPLICIT NONE
TYPE(AxFlx_TYPE) :: FLX(NMESH)
TYPE(MODALFLX_TYPE) :: MFLX(NMESH)
TYPE(PINXS_TYPE) :: XS(NMESH)
TYPE(AxGEOM_TYPE) :: GEOM
INTEGER :: NMESH, IG1, IG2, NG

REAL :: KSQ(2, NG), K(2, NG), A(2, NG), B(2, NG)
REAL :: KM(4, NG), SM(4, NG), PSOL(0:4, 2, NG)
REAL :: X, xh, xh0, dx, MMT(2, NG)
REAL :: activeH, psisum
INTEGER :: i, j, IG

activeH = 0
DO i = 1, NMESH
  IF(FLX(I)%PSI(0) .LT. 1.0E-5) CYCLE
  activeH = activeH + Geom%H(i)
  psisum = psisum + FLX(I)%PSI(0) * Geom%H(i)
ENDDO
psisum = psisum / activeH; psisum = 1._8 /psisum
 
xh0 = 0
open(unit = 99, file ='MMT.out', STATUS='REPLACE')
WRITE(99, '(15x, 500(I18,15x))')(IG, IG= ig1,ig2)
DO i = 1, NMESH
  DO IG = IG1, IG2
    KSQ(1, IG) = FLX(I)%KSQ(1,IG); KSQ(2, IG) = FLX(I)%KSQ(2,IG)
    K(1, IG) = SQRT(KSQ(1, IG)); K(2, IG) = SQRT(KSQ(2, IG))
    KM(:, IG) = KMAT(K(1, IG), K(2, IG)); SM(:, IG) = FLX(I)%S(:, IG)
    A(:, IG) = MFLX(I)%A(:, IG);  B(:, IG) = MFLX(I)%B(:, IG)
    PSOL(:, :, IG) = MFLX(I)%PSOL(:, :, IG)
    dx = 2._8/100
    x= -1._8 - dx
  ENDDO
  DO j = 1, 100
    x= x + dx
    xh = xh0 + 0.5*(x+1)*Geom%H(i)
    DO IG = IG1, IG2
      MMT(1:2, IG) = CAL_SP3FLX(A(:, IG), B(:, IG), PSOL(:, :, IG), KM(:, IG), SM(:, IG), x)
    ENDDO
    MMT = MMT * PSISUM
    write(99, '(I5, f10.5, 500(2e15.5, 2x, A))') I, xh, (MMT(1:2, IG), '|', IG = ig1, ig2)
  ENDDO
  xh0 = xh0 + GEOM%H(I)
ENDDO
!write(99, '(A)')'=============================================================================================='
CLOSE(99)
END SUBROUTINE

END MODULE

