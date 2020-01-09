SUBROUTINE GdFxrBurnUp(DeplVars, DeplLib, DeplCntl)
USE PARAM
USE DeplType,      ONLY : DeplVars_Type,      DeplLib_Type,        DeplCntl_Type,          &
                          MatExp_Type
USE Depl_mod,      ONLY : MatExp,             AllocMatExpVec,      InitDeplXS,             &
                          TranDepXS2DeplLib,  CopyIsoNumVector,    MakeDeplMat
USE MatExp_mod,    ONLY : MatExpSolver,       Mat_Type
USE GdDepl_mod,    ONLY : IsoList,            IsoLoc,                                       &
                          UpdateGdXs
USE BasicOperation,ONLY : CP_VA
IMPLICIT NONE

TYPE(DeplVars_Type) :: DeplVars
TYPE(DeplLib_Type) :: DeplLib
TYPE(DeplCntl_Type) :: DeplCntl

REAL, POINTER :: IsoNum_RK(:, :)
REAL :: N155(1:4)
INTEGER :: Gd155Idx

REAL, POINTER :: IsoNum(:), BurnUpXs(:, :)
REAL :: BurnUpTime, Phi1g, num_rk
INTEGER :: nIsoDep, nsubstep
INTEGER :: i, j
INTEGER :: Tid

Tid = DeplVars%tid

nIsoDep = DeplLib%nIsoDep
MatExp(Tid)%Mat => DeplVars%Dmat; MatExp(Tid)%nIsoDepl = nIsoDep
MatExp(Tid)%Solver = DeplCntl%Solver
IsoNum => DeplVars%IsoNum; BurnUpXs => DeplVars%BurnUpXS

IF(.NOT. MatExp(Tid)%lAllocVec) THEN
  CALL AllocMatExpVec(MatExp(Tid), DeplLib%nIsoDep)
ENDIF

ALLOCATE(IsoNum_RK(nIsoDep, 0:4))
    
Gd155Idx = IsoLoc(3)
nsubstep = DeplCntl%nSubStep
BurnUpTime = DeplCntl%Tsec / nsubstep
DO j = 1, nsubstep
  CALL CP_VA(IsoNum_RK(1:nIsoDep, 0), IsoNum(1:nIsoDep), nIsoDep)
  !1st step
  N155(1) = IsoNum(Gd155Idx)
  CALL  UpdateGdXs(DeplVars, N155(1))

  CALL InitDeplXS(DeplLib)
  CALL TranDepXS2DeplLib(DeplLib, BurnUpXs, IsoNum_RK(:, 0), nIsoDep)
  CALL CopyIsoNumVector(MatExp(Tid)%Viso0, IsoNum_RK(:, 0), nIsoDep, epsm30)
  phi1g = DeplVars%phi1g * 1.0E-24_8
  CALL MakeDeplMat(MatExp(Tid)%Mat, DeplLib, Phi1g, BurnUpTime)
  CALL MatExpSolver(MatExp(Tid))
  CALL Cp_VA(IsoNum_RK(1:nIsoDep, 1),  MatExp(Tid)%Viso(1:nIsoDep), nIsoDep)

  N155(2) = N155(1) + 0.5_8 * (IsoNum_RK(Gd155Idx, 1) - N155(1))

  CALL  UpdateGdXs(DeplVars, N155(2))
  CALL InitDeplXS(DeplLib)
  CALL TranDepXS2DeplLib(DeplLib, BurnUpXs, IsoNum_RK(:, 0), nIsoDep)
  CALL CopyIsoNumVector(MatExp(Tid)%Viso0, IsoNum_RK(:, 0), nIsoDep, epsm30)
  phi1g = DeplVars%phi1g * 1.0E-24_8
  CALL MakeDeplMat(MatExp(Tid)%Mat, DeplLib, Phi1g, BurnUpTime)
  CALL MatExpSolver(MatExp(Tid))
  CALL Cp_VA(IsoNum_RK(1:nIsoDep, 2),  MatExp(Tid)%Viso(1:nIsoDep), nIsoDep)

  N155(3) = N155(1) + 0.5_8 * (IsoNum_RK(Gd155Idx, 2) - N155(1))
  CALL  UpdateGdXs(DeplVars, N155(3))
  CALL InitDeplXS(DeplLib)
  CALL TranDepXS2DeplLib(DeplLib, BurnUpXs, IsoNum_RK(:, 0), nIsoDep)
  CALL CopyIsoNumVector(MatExp(Tid)%Viso0, IsoNum_RK(:, 0), nIsoDep, epsm30)
  phi1g = DeplVars%phi1g * 1.0E-24_8
  CALL MakeDeplMat(MatExp(Tid)%Mat, DeplLib, Phi1g, BurnUpTime)
  CALL MatExpSolver(MatExp(Tid))
  CALL Cp_VA(IsoNum_RK(1:nIsoDep, 3),  MatExp(Tid)%Viso(1:nIsoDep), nIsoDep)

  N155(4) = IsoNum_RK(Gd155Idx, 3)
  CALL  UpdateGdXs(DeplVars, N155(4))
  CALL InitDeplXS(DeplLib)
  CALL TranDepXS2DeplLib(DeplLib, BurnUpXs, IsoNum_RK(:, 0), nIsoDep)
  CALL CopyIsoNumVector(MatExp(Tid)%Viso0, IsoNum_RK(:, 0), nIsoDep, epsm30)
  phi1g = DeplVars%phi1g * 1.0E-24_8
  CALL MakeDeplMat(MatExp(Tid)%Mat, DeplLib, Phi1g, BurnUpTime)
  CALL MatExpSolver(MatExp(Tid))
  CALL Cp_VA(IsoNum_RK(1:nIsoDep, 4),  MatExp(Tid)%Viso(1:nIsoDep), nIsoDep)

  DO i = 1, nisodep
    num_RK = IsoNum_RK(i, 1) + 2.0_8 * IsoNum_RK(i, 2) + 2.0_8 * IsoNum_RK(i, 3) + IsoNum_RK(i, 4)
    num_RK = num_RK / 6.0_8
    IF(num_RK .LT. 1.e-40_8) num_RK =1.e-40_8
    IsoNum(i) = num_RK
  ENDDO
ENDDO
DEALLOCATE(IsoNum_RK)
END SUBROUTINE

SUBROUTINE GdXsFunction(DeplVars, DeplXs, Fxr, lCorrectStep)
USE PARAM
USE TypeDef,     ONLY : FxrInfo_Type
USE DeplType,    ONLY : DeplXs_Type,      DeplVars_Type
USE MAT2x2OP,    ONLY : DIR4x4
IMPLICIT NONE

TYPE(DeplVars_Type) :: DeplVars
TYPE(DeplXs_Type) :: DeplXS
TYPE(FxrInfo_Type) :: Fxr
LOGICAL :: lCorrectStep

REAL :: XsFtn(0:2, 64152:64160)
REAL :: RR(-1:1, 64152:64160), N155(-1:1)
REAL :: A(16), B(4), X(4) 
INTEGER :: i, j, id

!Init. RR
DO i = -1, 1
  DO j = 64152, 64160
    RR(i,j) = 0.
  ENDDO
ENDDO

DO i = -1, 0
  DO j = 1, Fxr%DeplXs1g(i)%niso
    id = Fxr%DeplXs1g(i)%idiso(j)
    IF(id .LT. 64152 .OR. id .GT. 64160) CYCLE
    RR(i, id) =  Fxr%DeplXs1g(i)%xsa(j)
!    RR(i, id) =  Fxr%DeplXs1g(i)%xsf(j)
!    RR(i, id) =  Fxr%DeplXs1g(i)%xsn2n(j)
!    RR(i, id) =  Fxr%DeplXs1g(i)%xsn3n(j)
    RR(i, id) = RR(i, id) * Fxr%DeplXs1g(i)%phi1g * 1.0E-24_8
    IF(id == 64155) N155(i) = Fxr%DeplXs1g(i)%n155
  ENDDO
ENDDO
i=1
DO j = 1, Fxr%niso_depl
  id = Fxr%idiso(j)
  IF(id .LT. 64152 .OR. id .GT. 64160) CYCLE
  RR(i, id) =  DeplXs%xsa(j)
!  RR(i, 2, id) =  DeplXs%xsf(j)
!  RR(i, 3, id) =  DeplXs%xsn2n(j)
!  RR(i, 4, id) =  DeplXs%xsn3n(j)
  RR(i, id) = RR(i, id) * DeplXs%phi1g * 1.0E-24_8
  IF(id == 64155) N155(i) = Fxr%pnum(j)
ENDDO

!#define GdXS_LINEAR
#ifdef GdXS_LINEAR
DO id = 64152, 64160
  IF(id .EQ. 64153) CYCLE
  IF(id .EQ. 64159) CYCLE
  A = 0; A(16) = 1
  A(1) = 1; A(2) = N155(0)
  A(5) = 1; A(6) = N155(1)
  B(1) = RR(0, id); B(2) = RR(1, id)
    
  x(1) = (A(6)*B(1) - A(2)*B(2))/(A(1)*A(6)-A(2)*A(5))
  x(2) = (-A(5)*B(1) + A(1)*B(2))/(A(1)*A(6)-A(2)*A(5))
    
  x(1:2) = x(1:2) / (DeplXs%phi1g * 1.0E-24_8)
  XsFtn(0:1, id) = x(1:2)
  XsFtn(2, id) = 0.
ENDDO
#else
DO id = 64152, 64160
  IF(id .EQ. 64153) CYCLE
  IF(id .EQ. 64159) CYCLE
  A = 0; A(16) = 1
  A(1) = 1; A(2) = N155(-1);  A(3) = N155(-1) * N155(-1)
  A(5) = 1; A(6) = N155( 0);  A(7) = N155( 0) * N155( 0)
  A(9) = 1; A(10) = N155( 1); A(11) = N155( 1) * N155( 1)
  B(1) = RR(-1, id); B(2) = RR(0, id); B(3) = RR(1, id); B(4) = 1
  CALL Dir4x4(A, b, x)
  x(1:3) = x(1:3) / (DeplXs%phi1g * 1.0E-24_8)
  XsFtn(0:2, id) = x(1:3)
ENDDO
#endif
!!
!DO id = 64152, 64160
!  IF(id .EQ. 64153) CYCLE
!  IF(id .EQ. 64159) CYCLE
!  A = 0; A(16) = 1; A(11) = 1
!  IF(lCorrectStep) THEN 
!    A(1) = 1; A(2) = N155(0)
!    A(5) = 1; A(6) = N155(1)
!    B(1) = RR(0, id); B(2) = RR(1, id); B(3) = 1; B(4) = 1
!  ELSE
!    A(1) = 1; A(2) = N155(0)
!    A(5) = 1; A(6) = N155(1)
!    B(1) = RR(0, id); B(2) = RR(1, id); B(3) = 1; B(4) = 1  
!  ENDIF
!  CALL Dir4x4(A, b, x)
!  X(3:4) = 0
!  x(1:3) = x(1:3) / (DeplXs%phi1g * 1.0E-24_8)
!  XsFtn(0:2, id) = x(1:3)
!ENDDO

DeplVars%GdXsFtn(0:2, 64152:64160) = XsFtn(0:2, 64152:64160)
END SUBROUTINE

SUBROUTINE UpdateGdXs(DeplVars, N155)
USE PARAM
USE DeplType,    ONLY : DeplVars_Type

USE GdDepl_mod,    ONLY : IsoList, IsoLoc
IMPLICIT NONE
TYPE(DeplVars_Type) :: DeplVars
REAL :: N155

REAL :: XsFtn(0:2, 64152:64160)
REAL :: xs1g, N1550, f

INTEGER :: i, id

!  DO i = 1, 7
!    id = IsoList(i); id_lib = iposiso(id) 
!    IsoLoc(i) = DeplVars%MapXs2Dep(id_lib)
!  ENDDO
XsFtn = DeplVars%GdXsFtn
#ifdef GdXS_LINEAR
DO i = 1, 7
  id = IsoList(i)
  xs1g = XsFtn(0, id) + N155 * XsFtn(1, id)
  f = 1._8 / DeplVars%BurnUpXs(1, IsoLoc(i))
  DeplVars%BurnUpXs(1, IsoLoc(i)) = xs1g
  f = xs1g * f
  DeplVars%BurnUpXs(2:4, IsoLoc(i)) = DeplVars%BurnUpXs(2:4, IsoLoc(i)) * f
ENDDO
#else
DO i = 1, 7
  id = IsoList(i)
  xs1g = XsFtn(0, id) + N155 * XsFtn(1, id) + N155 * N155 * XsFtn(2,id)
  f = 1._8 / DeplVars%BurnUpXs(1, IsoLoc(i))
  DeplVars%BurnUpXs(1, IsoLoc(i)) = xs1g
  f = xs1g * f
  DeplVars%BurnUpXs(2:4, IsoLoc(i)) = DeplVars%BurnUpXs(2:4, IsoLoc(i)) * f
ENDDO
#endif

END SUBROUTINE

SUBROUTINE UpdateDeplGdFxrInfo(Fxr, DeplVars, GroupInfo, lCorrectStep, lXeDyn)
USE PARAM
USE TYPEDEF,            ONLY : FxrInfo_Type,       GroupInfo_Type
USE DeplType,           ONLY : DeplVars_Type
USE BasicOperation,     ONLY : MULTI_CA,  CP_CA
USE nuclidmap_mod,      ONLY : iposiso,   PndCrit, nuclide
USE GdDepl_mod,         ONLY : IsoList,           IsoLoc
USE XsUtil_mod,         ONLY : SetXeDynEnv
IMPLICIT NONE
TYPE(FxrInfo_Type) :: Fxr
TYPE(DeplVars_Type) :: DeplVars
TYPE(GroupInfo_Type) :: GroupInfo
LOGICAL :: lCorrectStep, lXeDyn

INTEGER, POINTER :: MapXs2Dep(:), MapDep2Xs(:), idiso(:)
REAL, POINTER :: IsoNum(:), pnum(:)
REAL :: IsoNum_Gd(7)

INTEGER :: i, j, IdXs, IdDep
INTEGER :: niso, nisoDepl, ntiso
INTEGER :: niso_wolib, idiso_wolib(100)
REAL ::  pnum_wolib(100)

MapXs2Dep => DeplVars%MapXs2Dep
MapDep2Xs => DeplVars%MapDep2Xs
IsoNum => DeplVars%IsoNum
nisodepl = DeplVars%nIsoDepl
IdIso => Fxr%IdIso
pnum => Fxr%pnum

IF(lCorrectStep) THEN
  DO i = 1, 7
    j = IsoLoc(i)
    IsoNum_Gd(i) = IsoNum(j)
  ENDDO
  
  CALL MULTI_CA(0.5_8, IsoNum(1:nIsoDepl), nIsoDepl)
  IF(Fxr%l_pnum_all) THEN ! Depletion TimeStep Bug Fixed
    niso = Fxr%nIso_depl
    DO i = 1, nIsoDepl
      IsoNum(i) = IsoNum(i) + 0.5_8 * Fxr%pnum_all(i)
    END DO ! of i
  ELSE
    niso = Fxr%nIso_depl
    DO i = 1, niso
      IdXs = iposiso(IdIso(i))
      IdDep = MapXs2Dep(IdXs)
      IF(IdDep .EQ. 0) CYCLE
      IsoNum(IdDep) = IsoNum(IdDep) + 0.5_8 * pnum(i)
    ENDDO
  END IF
  
  DO i = 1, 7
    j = IsoLoc(i)
    IsoNum(j) = IsoNum_Gd(i)
  ENDDO
ENDIF

niso_wolib = 0 !Depletion TimeStep Bug Fixed
DO i = 1, Fxr%nIso
  IdXs = iPosIso(IdIso(i))
  IdDep = MapXs2Dep(IdXs)
  IF(IdDep /= 0) CYCLE
  nIso_wolib = nIso_wolib + 1
  idiso_wolib(niso_wolib) = IdIso(i)
  pnum_wolib(niso_wolib) = pnum(i)
END DO ! of i

ntiso = GroupInfo%ntiso

CALL CP_CA(IdIso, 0, ntiso)
CALL CP_CA(pnum, 0._8, ntiso)
niso = 0
DO i = 1, nIsoDepl
  j = MapDep2Xs(i)
  IF(j .EQ. 0) CYCLE
  IF(abs(IsoNum(i)) .LT. epsm20) CYCLE
  IF(IsoNum(i) .LT. pndcrit(j)) CYCLE
  niso = niso + 1
  IdXs = nuclide(j)
  pnum(niso) = IsoNum(i)
  idiso(niso) = IdXs
ENDDO

DO i = 1, nIsoDepl
  Fxr%pnum_all(i) = IsoNum(i)
END DO ! of i

IF(niso_wolib > 0) THEN
  DO i = 1, niso_wolib
    niso = niso + 1
    pnum(niso) = pnum_wolib(i)
    idiso(niso) = idiso_wolib(i)
  END DO ! of i
END IF

Fxr%niso = niso
Fxr%niso_depl = niso

IF(lXeDyn) CALL SetXeDynEnv(IdIso, pnum, Fxr%niso, Fxr%niso_depl, ntiso)
NULLIFY(MapXs2Dep, IsoNum)
NULLIFY(pnum, idiso)
END SUBROUTINE

SUBROUTINE UpdateGdDeplMat(DMat, DeplLib, PHI, BurnUpTime)
USE PARAM
USE DEPLTYPE,      ONLY : Mat_TYPE,         DeplLib_Type,                    &
                          AtomKind,         STATEKIND,        ISOTOPEKIND,   &
                          FisYield_Type      
USE BasicOperation, ONLY : CP_CA,           MULTI_CA
IMPLICIT NONE
TYPE(Mat_Type) :: DMat
TYPE(DeplLib_Type) :: DeplLib
REAL :: PHI, BurnUpTime

TYPE(AtomKind), POINTER :: Lib1(:)
TYPE(StateKind), POINTER :: MyStat, MyStatFr
TYPE(IsotopeKind), POINTER :: MyIso
TYPE(FisYield_Type), POINTER :: FisYield(:)


REAL, POINTER :: DIAG(:), OffDiag(:,:)              !Diag and Off Diag Elements for Sparse Matrix
INTEGER, POINTER :: nlmnt(:), lmntIdx(:,:)          !#of element for row, non-zero element index
INTEGER :: nmaxoffdiag 

REAL :: Y

INTEGER :: nIsoDep, nAtomDep
INTEGER :: I, J, K, M, JB, JE
INTEGER :: IM, IP, IFR, ITO, IZFR, IAFR

Lib1 => DeplLib%AtomLib1
FisYield => DeplLib%FisYield

nIsoDep = DeplLib%nIsoDep; nAtomDep = DeplLib%nAtomDep

Diag => DMAT%DIAG; OffDiag => DMAT%OffDiag
nLmnt => DMAT%nLmnt; lmntIdx => DMAT%LmntIdx
nMaxOffDiag = DMAT%nMaxOffDiag

CALL CP_CA(Diag, 0._8, nIsoDep)
CALL CP_CA(OffDiag, 0._8, nMaxOffDiag, nIsoDep)

I = 64
DO J = Lib1(I)%IB, Lib1(I)%IE
  MyIso => Lib1(I)%A(J)
  DO K = 0, MyIso%NSTATE
    IF(MyIso%STAT(K)%IMAT .EQ. 0) CYCLE
    MyStat => MyIso%STAT(K)
    IM = MyStat%IMAT
    Diag(IM) = -PHI * MyStat%XS(0) - MyStat%Rambda   ! LOSS
    DO M = 1, MyStat%NTO1    !  SOURCE FROM DECAY CHAIN
      IFR = MyStat%ITO1(1, M); ITO = MyStat%ITO1(2, M); IP = MyStat%ITO1(3, M)
      OffDiag(IP, ITO) = OffDiag(IP, ITO) + MyStat%FRAC(IFR)
    ENDDO

    DO M = 1, MyStat%NTO2   !  SOURCE FROM NUCLEAR REACTION
      IFR = MyStat%ITO2(1, M); ITO = MyStat%ITO2(2, M); IP = MyStat%ITO2(3, M)
      OffDiag(IP, ITO) = OffDiag(IP, ITO) + PHI * MyStat%XS(IFR)
    ENDDO
    
    IF(MySTAT%IGRP .EQ. 3) THEN
      DO M = 1, MyStat%NFR3  !   SOURCE FROM FISSION
        IFR = MyStat%IFR3(1, M); IP = MyStat%IFR3(2, M); Y = MyStat%Y(IFR)
        IZFR = FisYield(IFR)%AtomNum; IAFR = FisYield(IFR)%AtomWeight
        MyStatFr => Lib1(IZFR)%A(IAFR)%STAT(0)
        OffDiag(IP, IM) = OffDiag(IP, IM) + Y * PHI * MyStatFr%XS(4)
        !IZFRIAFR
      ENDDO
    ENDIF
    Diag(IM) = Diag(IM) * BurnUpTime
    M = nLmnt(I)
    OffDiag(1:M, IM) = OffDiag(1:M, IM) * BurnUpTime 
  ENDDO
ENDDO

IF(ASSOCIATED(MyStat)) NULLIFY(MyStat)
IF(ASSOCIATED(MyStatFr)) NULLIFY(MyStatFr)


NULLIFY(FisYield);  NULLIFY(Lib1)
NULLIFY(Diag); NULLIFY(OffDiag)
NULLIFY(nLmnt); NULLIFY(LmntIdx)
END SUBROUTINE

!SUBROUTINE PCPout(Fxr, nFxr)
!USE PARAM
!USE TYPEDEF,    ONLY : FxrInfo_Type
!TYPE(FxrInfo_Type) :: Fxr(nFxr)
!INTEGER :: nFxr
!REAL :: dat155(1000), dat157(1000)
!id =0
!dat155=0; dat157=0
!DO i = 1, nFxr
!  IF(.NOT. Fxr(i)%lGd) CYCLE
!  id = id + 1
!  dat155(id) = Fxr(i)%DeplPCP%numdat(1, 64155)
!  dat157(id) = Fxr(i)%DeplPCP%numdat(1, 64157)
!  id = id + 1
!  dat155(id) = Fxr(i)%DeplPCP%numdat(2, 64155)
!  dat157(id) = Fxr(i)%DeplPCP%numdat(2, 64157)
!  id = id + 1
!  dat155(id) = Fxr(i)%DeplPCP%numdat(3, 64155)
!  dat157(id) = Fxr(i)%DeplPCP%numdat(3, 64157)
!  id = id + 1
!  dat155(id) = Fxr(i)%DeplPCP%f(64155)
!  dat157(id) = Fxr(i)%DeplPCP%numdat(3, 64157)
!!  dat155(1, id) = Fxr(i)%DeplPCP%f_r(64155)
!!  dat155(2, id) = Fxr(i)%DeplPCP%f_c(64155)
!!  dat157(1, id) = Fxr(i)%DeplPCP%f_r(64157)
!!  dat157(2, id) = Fxr(i)%DeplPCP%f_c(64157)
!ENDDO
!1234 format(1000(4e15.6, 2x, '/', 2x))
!write(155, 1234) dat155(1:id)
!write(157, 1234) dat157(1:id)
!!write(156, '(1000f10.6)') dat155(2, 1:id)
!!write(157, '(1000f10.6)') dat157(1, 1:id)
!!write(158, '(1000f10.6)') dat157(2, 1:id)
!END SUBROUTINE