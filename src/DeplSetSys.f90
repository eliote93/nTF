! Get mass of heavy metals : SetHmkg0
! Get the time step of depletion cal. : GetDelT
! Set a bunch of depletion systems : SetDeplSys
SUBROUTINE SetHmkg0(DeplFxrBundle, DeplLib)
  USE HPDeplType
  IMPLICIT NONE
  TYPE(DeplFxrBundle_Type) :: DeplFxrBundle
  TYPE(DeplLib_Type) :: DeplLib
  TYPE(DeplFxr_Type), POINTER :: Fxrs(:)
  INTEGER :: nfxr, nisoeig
  INTEGER :: i, j, k
  REAL(8), PARAMETER :: invAvgdr = 1.660539040e-24 ! 1./6.022140857e+23
  REAL(8), PARAMETER :: kgperg = 1.e-3
  REAL(8), PARAMETER :: barnPcm2 = 1.e+24
  REAL(8) :: Hmkg0, Hmkg0_fxr

  nfxr = DeplFxrBundle%nfxr
  Fxrs => DeplFxrBundle%FxrBundle
  Hmkg0 = 0._8
  DO i = 1, nfxr
    IF (.NOT. Fxrs(i)%lFuel) CYCLE
    nisoeig = Fxrs(i)%NisoEig
    DO j = 1, nisoeig
      k = Fxrs(i)%IdIsoEig(j)
      IF (.NOT. DeplLib%lActinide(k)) CYCLE
      Hmkg0_fxr = Fxrs(i)%pnum_sseig(j)*barnPcm2*Fxrs(i)%Vol*dble(MOD(DeplLib%Idiso(k),10000)/10)*invAvgdr*kgperg
      Fxrs(i)%Hmkg0 = Fxrs(i)%Hmkg0+Hmkg0_fxr
      Hmkg0 = Hmkg0 + Hmkg0_fxr
    END DO
  END DO
  DeplFxrBundle%Hmkg0 = Hmkg0
END SUBROUTINE

SUBROUTINE GetDelT(DeplFxrBundle, BurnUP, lDelBU)
  USE HPDeplType
  IMPLICIT NONE
  TYPE(DeplFxrBundle_Type) :: DeplFxrBundle
  REAL(8) :: BurnUP  ! lDelBU = .TRUE. -> BurnUP(MWD/kg) / lDelBU = .FALSE. -> BurnUP(days)
  LOGICAL :: lDelBU

  IF (lDelBU) THEN
    DeplFxrBundle%delT = BurnUP*DeplFxrBundle%Hmkg0/DeplFxrBundle%Power*86400._8
  ELSE
    DeplFxrBundle%delT = BurnUP*86400._8
  END IF
END SUBROUTINE

INTEGER FUNCTION CalSizeSysBun(DeplLib, ByteSys, Scale, NsysMax, IsCRAM)
  USE HPDeplType
#ifdef __PGI
  USE CUDA_Workspace, ONLY : Grid, Bloc, lGS, lILP
#endif
  IMPLICIT NONE
  TYPE(DeplLib_Type) :: DeplLib
  INTEGER :: ByteSys, Scale, NsysMax
  LOGICAL :: IsCRAM
  INTEGER, PARAMETER :: Precision = 8

  INTEGER :: NNZ, NR
  REAL(8) :: SizeSys, BunSize

  BunSize = DBLE(ByteSys)*10.**DBLE(Scale)
  NNZ = SIZE(DeplLib%YldMapColIdx); NR = SIZE(DeplLib%YldMapRowPtr)
  SizeSys = DBLE(NNZ*Precision);
#ifdef __PGI
  SizeSys = DBLE(Precision*(NNZ*2+NR*3))
  IF (lGS) THEN
    IF (.NOT. lILP) BunSize = BunSize-DBLE(Precision*(2+1)*Grid*Bloc)
    IF (lILP) BunSize = BunSize-DBLE(Precision*(2*8+1)*Grid*Bloc)
  ELSE
    BunSize = BunSize-DBLE(Precision*2*9*Grid*Bloc)
  END IF
#else
  IF (IsCRAM) THEN
    SizeSys = SizeSys*2.;
  END IF
#endif
  CalSizeSysBun = MIN(INT(BunSize/SizeSys), NsysMax)
END FUNCTION

SUBROUTINE SetDeplSys(DeplLib, DeplFxrBundle, Nsys, ifxrbeg, DeplSysBundle, lGd, Nths, nSubStp)
USE HPDeplType
USE CSRMATRIX
IMPLICIT NONE
TYPE(DeplLib_Type) :: DeplLib
TYPE(DeplFxrBundle_Type) :: DeplFxrBundle
INTEGER :: Nsys, ifxrbeg
TYPE(DeplSysBundle_Type) :: DeplSysBundle
LOGICAL :: lGd
INTEGER :: Nths
INTEGER, OPTIONAL :: nSubStp

TYPE(DeplFxr_Type), POINTER :: Fxrs(:)
TYPE(CSR_DOUBLE), POINTER :: DeplMats(:)
REAL(8), POINTER :: pnums(:), pnums0(:)
REAL(8), POINTER :: xs1g(:,:)
REAL(8), POINTER :: DecFrac(:,:)
REAL(8), POINTER :: FisFrac(:,:)
REAL(8), ALLOCATABLE :: FisXs(:)
INTEGER, POINTER :: IdAftRct(:,:), IdAftDec(:,:)
INTEGER, POINTER :: IdFis(:), IdYld(:)
INTEGER, POINTER :: SubYldRct(:), SubYldRctAct(:), SubYldDec(:)
LOGICAL, POINTER :: lActinide(:)

INTEGER, POINTER :: IdTrue(:)

INTEGER :: ifxr, i ,j, k, l, ir, ic, icbeg, icend
INTEGER :: NofIso, NofFis, NofYld, NofDec, NofRct, FisRctId
REAL(8) :: LossRate, YldRate, phi, delT

REAL(8), POINTER :: csrval(:)
INTEGER, POINTER :: csrColIdx(:), csrRowPtr(:)
INTEGER :: nnz, nr, nc

REAL(8) :: SubStep

SubStep = 1.
IF (PRESENT(nSubStp)) SubStep = nSubStp

DeplSysBundle%Nsys = Nsys; DeplSysBundle%ifxrbeg = ifxrbeg
Fxrs => DeplFxrBundle%FxrBundle
delT = DeplFxrBundle%delT/SubStep
!IF (SubStep .GT. 1.) print*, SubStep, delT, DeplFxrBundle%delT

NofIso = DeplLib%NofIso; NofFis = DeplLib%NofFis; NofYld = DeplLib%NofYld
NofDec = DeplLib%NofDec; NofRct = DeplLib%NofRct; FisRctId = DeplLib%FisRctId

DeplSysBundle%NofIso = NofIso

DecFrac => DeplLib%DecFrac; FisFrac => DeplLib%FisFrac
IdAftRct=>DeplLib%IdAftRct; IdAftDec=>DeplLib%IdAftDec
IdFis => DeplLib%IdFis; IdYld => DeplLib%IdYld
SubYldRct => DeplLib%SubYldRct; SubYldRctAct => DeplLib%SubYldRctAct; SubYldDec => DeplLib%SubYldDec
lActinide => DeplLib%lActinide

ALLOCATE(DeplSysBundle%DeplMats(Nsys), DeplSysBundle%pnums(NofIso*Nsys))
ALLOCATE(DeplSysBundle%pnums1(NofIso*Nsys))
DeplMats => DeplSysBundle%DeplMats
pnums => DeplSysBundle%pnums
pnums = 0.;

nr = NofIso; nc = NofIso; nnz = DeplLib%YldMapRowPtr(nr+1)-1

IF (lGD) THEN
  IdTrue => DeplFxrBundle%IdTrueGd
ELSE
  IdTrue => DeplFxrBundle%IdTrueDepl
END IF

!ALLOCATE(FisXs(NofFis))
CALL OMP_SET_NUM_THREADS(Nths)
!$OMP PARALLEL DO PRIVATE(xs1g, phi, pnums0, csrval, csrRowPtr, csrColIdx, j, k, l, icbeg, icend, ifxr, ir, ic, YldRate, LossRate, FisXs) SCHEDULE(GUIDED)
DO i = 1, Nsys
  ifxr = IdTrue(i-1+ifxrbeg)
  ALLOCATE(FisXs(NofFis))
  xs1g => Fxrs(ifxr)%xs1g
  phi = Fxrs(ifxr)%NormFlux1g

  pnums0=>Fxrs(ifxr)%pnum_depl
  pnums(NofIso*(i-1)+1:NofIso*i) = pnums0(:)
  NULLIFY(pnums0)
  CALL createCSR(DeplMats(i), nnz, nr, nc)
  csrval=>DeplMats(i)%csrVal; csrColIdx=>DeplMats(i)%csrColIdx
  csrRowPtr=>DeplMats(i)%csrRowPtr
  csrval = 0.
  csrColIdx(:) = DeplLib%YldMapColIdx(:)
  csrRowPtr(:) = DeplLib%YldMapRowPtr(:)

  ! Decay
  DO j = 1, NofIso
    LossRate=0.
    ic = j
    DO k = 1, NofDec
      ir = IdAftDec(k,j)
      YldRate = DecFrac(k,j)*delT
      LossRate = LossRate+YldRate
      IF (ir .EQ. 0) CYCLE
      icbeg = csrRowPtr(ir); icend = csrRowPtr(ir+1)-1
      DO l = icbeg, icend
        IF (csrColIdx(l) .NE. ic) CYCLE
        csrval(l) = csrval(l)+YldRate
        EXIT
      END DO
    END DO
    icbeg = csrRowPtr(j); icend = csrRowPtr(j+1)-1
    DO l = icbeg, icend
      IF (csrColIdx(l) .NE. j) CYCLE
      csrval(l) = csrval(l)-LossRate
      !if(i.EQ.1)print*, j, csrColIdx(l)
    END DO
  END DO

  ! Reaction
  DO j = 1, NofIso
    LossRate = 0.
    ic = j
    DO k = 1, NofRct
      ir = IdAftRct(k,j)
      YldRate = xs1g(k,j)*phi*delT
      LossRate = LossRate+YldRate
      If (ir .EQ. 0) CYCLE
      icbeg = csrRowPtr(ir); icend = csrRowPtr(ir+1)-1
      DO l = icbeg, icend
        IF (csrColIdx(l) .NE. ic) CYCLE
        csrval(l) = csrval(l)+YldRate
        EXIT
      END DO
    END DO
    icbeg = csrRowPtr(j); icend = csrRowPtr(j+1)-1
    DO l = icbeg, icend
      IF (csrColIdx(l) .NE. j) CYCLE
      csrval(l) = csrval(l)-LossRate
      !if(i.eq.1)print*, j, csrColIdx(l)
    END DO
  END DO

  ! Fission
  DO j = 1, NofFis
    FisXs(j) = xs1g(FisRctId, IdFis(j))
  END DO
  DO j = 1, NofYld
    ir = IdYld(j)
    DO k = 1, NofFis
      YldRate = FisFrac(k,j)*FisXs(k)*phi*delT
      ic = IdFis(k)
      icbeg = csrRowPtr(ir); icend = csrRowPtr(ir+1)-1
      DO l = icbeg, icend
        IF (csrColIdx(l) .NE. ic) CYCLE
        csrval(l) = csrval(l)+YldRate
        EXIT
      END DO
    END DO
  END DO

  ! SubYld
  DO j = 1, NofRct
    DO k = 1, NofIso
      ic = k
      IF (lActinide(k)) THEN
        ir = SubYldRctAct(j)
      ELSE
        ir = subYldRct(j)
      END IF
      IF(ir .EQ. 0) CYCLE
      icbeg = csrRowPtr(ir); icend = csrRowPtr(ir+1)-1
      YldRate = xs1g(j,k)*phi*delT
      DO l = icbeg, icend
        IF (csrColIdx(l) .NE. ic) CYCLE
        csrval(l) = csrval(l)+YldRate
        EXIT
      END DO
    END DO
  END DO
  DO j = 1, NofDec
    ir = SubYldDec(j)
    IF (ir .EQ. 0) CYCLE
    DO k = 1, NofIso
      ic = k
      icbeg = csrRowPtr(ir); icend = csrRowPtr(ir+1)-1
      YldRate = DecFrac(j,k)*delT
      DO l = icbeg, icend
        IF (csrColIdx(l) .NE. ic) CYCLE
        csrval(l) = csrval(l)+YldRate
      END DO
    END DO
  ENDDO
  DeplMats(i)%nnz = nnz; DeplMats(i)%lFinalized = .TRUE.
  !IF(i.EQ.9) CALL printCSR(DeplMats(i), 'CSROUT.txt',671)
  !CALL finalizeCsr(DeplMats(i), .FALSE.)
  NULLIFY(xs1g, csrval, csrRowPtr, csrColIdx)
  DEALLOCATE(FisXs)
END DO
!$OMP END PARALLEL DO

NULLIFY(lActinide, FisFrac, DecFrac, SubYldRct, SubYldRctAct, SubYldDec)
NULLIFY(IdAftRct, IdAftDec, IdFis, IdYld)
NULLIFY(pnums, DeplMats, Fxrs)

END SUBROUTINE

SUBROUTINE DestroySys(DeplSysBundle)
  USE HPDeplType
  IMPLICIT NONE
  TYPE(DeplSysBundle_Type) :: DeplSysBundle
  INTEGER :: Nsys, i

  Nsys = DeplSysBundle%Nsys
  DO i = 1, Nsys
    CALL destroyCSR(DeplSysBundle%DeplMats(i))
  END DO
  DEALLOCATE(DeplSysBundle%DeplMats)
  DEALLOCATE(DeplSysBundle%pnums, DeplSysBundle%pnums1)
END SUBROUTINE

#ifdef __PGI
SUBROUTINE SetDeplSys_woCSRT(DeplLib, DeplFxrBundle, Nsys, ifxrbeg, DeplSysBundle, lGd, Nths, nSubStp)
USE HPDeplType
USE CSRMATRIX
IMPLICIT NONE
TYPE(DeplLib_Type) :: DeplLib
TYPE(DeplFxrBundle_Type) :: DeplFxrBundle
INTEGER :: Nsys, ifxrbeg
TYPE(DeplSysBundle_Type) :: DeplSysBundle
LOGICAL :: lGd
INTEGER :: Nths
INTEGER, OPTIONAL :: nSubStp

TYPE(DeplFxr_Type), POINTER :: Fxrs(:)
REAL(8), POINTER :: pnums(:), pnums0(:)
REAL(8), POINTER :: xs1g(:,:)
REAL(8), POINTER :: DecFrac(:,:)
REAL(8), POINTER :: FisFrac(:,:)
REAL(8), ALLOCATABLE :: FisXs(:)
INTEGER, POINTER :: IdAftRct(:,:), IdAftDec(:,:)
INTEGER, POINTER :: IdFis(:), IdYld(:)
INTEGER, POINTER :: SubYldRct(:), SubYldRctAct(:), SubYldDec(:)
LOGICAL, POINTER :: lActinide(:)

INTEGER, POINTER :: IdTrue(:)

INTEGER :: ifxr, i ,j, k, l, ir, ic, icbeg, icend
INTEGER :: NofIso, NofFis, NofYld, NofDec, NofRct, FisRctId
REAL(8) :: LossRate, YldRate, phi, delT

!COMPLEX(8), POINTER :: Diag(:,:)
!REAL(8), POINTER :: OffDiag(:,:)
!INTEGER, POINTER :: csrColIdx(:), csrRowPtr(:)
INTEGER :: nnz, nndz, nr, nc

REAL(8) :: SubStep

SubStep = 1.
IF (PRESENT(nSubStp)) SubStep = nSubStp

DeplSysBundle%Nsys = Nsys; DeplSysBundle%ifxrbeg = ifxrbeg
Fxrs => DeplFxrBundle%FxrBundle
delT = DeplFxrBundle%delT/SubStep
!IF (SubStep .GT. 1.) print*, SubStep, delT, DeplFxrBundle%delT

NofIso = DeplLib%NofIso; NofFis = DeplLib%NofFis; NofYld = DeplLib%NofYld
NofDec = DeplLib%NofDec; NofRct = DeplLib%NofRct; FisRctId = DeplLib%FisRctId

DeplSysBundle%NofIso = NofIso

DecFrac => DeplLib%DecFrac; FisFrac => DeplLib%FisFrac
IdAftRct=>DeplLib%IdAftRct; IdAftDec=>DeplLib%IdAftDec
IdFis => DeplLib%IdFis; IdYld => DeplLib%IdYld
SubYldRct => DeplLib%SubYldRct; SubYldRctAct => DeplLib%SubYldRctAct; SubYldDec => DeplLib%SubYldDec
lActinide => DeplLib%lActinide

!csrColIdx => DeplLib%YldMapColIdx; csrRowPtr => DeplLib%YldMapRowPtr
nr = NofIso; nc = NofIso; nnz = DeplLib%YldMapRowPtr(nr+1)-1
nndz = nnz - nr; ! number of non-diagonal-zero

!ALLOCATE(DeplSysBundle%DeplMats(Nsys))
ALLOCATE(DeplSysBundle%Diag(Nsys,nr), DeplSysBundle%OffDiag(Nsys,nndz))
ALLOCATE(DeplSysBundle%pnums(NofIso*Nsys), DeplSysBundle%pnums1(NofIso*Nsys))
!DeplMats => DeplSysBundle%DeplMats
pnums => DeplSysBundle%pnums
pnums = 0.;
!Diag => DeplSysBundle%Diag; OffDiag => DeplSysBundle%OffDiag;
DeplSysBundle%Diag = 0.; DeplSysBundle%OffDiag = 0.;

IF (lGD) THEN
  IdTrue => DeplFxrBundle%IdTrueGd
ELSE
  IdTrue => DeplFxrBundle%IdTrueDepl
END IF

!ALLOCATE(FisXs(NofFis))
CALL OMP_SET_NUM_THREADS(Nths)
!$OMP PARALLEL DO PRIVATE(xs1g, phi, pnums0, j,k,l,icbeg,icend,ifxr,ir, ic, YldRate, LossRate, FisXs) SCHEDULE(GUIDED)
DO i = 1, Nsys
  ifxr = IdTrue(i-1+ifxrbeg)
  ALLOCATE(FisXs(NofFis))
  xs1g => Fxrs(ifxr)%xs1g
  phi = Fxrs(ifxr)%NormFlux1g

  pnums0=>Fxrs(ifxr)%pnum_depl
  pnums(NofIso*(i-1)+1:NofIso*i) = pnums0(:)
  NULLIFY(pnums0)
!  CALL createCSR(DeplMats(i), nnz, nr, nc)
!  csrColIdx(:) = DeplLib%YldMapColIdx(:)
!  csrRowPtr(:) = DeplLib%YldMapRowPtr(:)

  ! Decay
  DO j = 1, NofIso
    LossRate=0.
    ic = j
    DO k = 1, NofDec
      ir = IdAftDec(k,j)
      YldRate = DecFrac(k,j)*delT
      LossRate = LossRate+YldRate
      IF (ir .EQ. 0) CYCLE
      l = DeplLib%IzDec(k,j)
      IF (l .LE. 0) CYCLE
      l = l-ir+(ir/ic)
!      print*, 'DECAY', k*1000+j, l
      DeplSysBundle%OffDiag(i,l) = DeplSysBundle%OffDiag(i,l)+YldRate
!      icbeg = DeplLib%YldMapRowPtr(ir); icend = DeplLib%YldMapRowPtr(ir+1)-1
!      DO l = icbeg, icend
!        IF (DeplLib%YldMapColIdx(l) .EQ. ic) THEN
!          ir = l - ir + (ir/ic)
!          DeplSysBundle%OffDiag(i,ir) = DeplSysBundle%OffDiag(i,ir)+YldRate
!          EXIT
!        END IF
!      END DO
    END DO
    DeplSysBundle%Diag(i,j) = DeplSysBundle%Diag(i,j)-LossRate
  END DO

  ! Reaction
  DO j = 1, NofIso
    LossRate = 0.
    ic = j
    DO k = 1, NofRct
      ir = IdAftRct(k,j)
      YldRate = xs1g(k,j)*phi*delT
      LossRate = LossRate+YldRate
      If (ir .EQ. 0) CYCLE
      l = DeplLib%IzRct(k,j)
      IF (l .LE. 0) CYCLE
      l = l-ir+(ir/ic)
!      print*, 'React', k*1000+j, l
      DeplSysBundle%OffDiag(i,l) = DeplSysBundle%OffDiag(i,l)+YldRate
!      icbeg = DeplLib%YldMapRowPtr(ir); icend = DeplLib%YldMapRowPtr(ir+1)-1
!      DO l = icbeg, icend
!        IF (DeplLib%YldMapColIdx(l) .EQ. ic) THEN
!          ir = l - ir + (ir/ic)
!          DeplSysBundle%OffDiag(i,ir) = DeplSysBundle%OffDiag(i,ir)+YldRate
!          EXIT
!        END IF
!      END DO
    END DO
    DeplSysBundle%Diag(i,j) = DeplSysBundle%Diag(i,j)-LossRate
  END DO

  ! Fission
  DO j = 1, NofFis
    FisXs(j) = xs1g(FisRctId, IdFis(j))
  END DO
  DO j = 1, NofYld
    ir = IdYld(j)
    DO k = 1, NofFis
      YldRate = FisFrac(k,j)*FisXs(k)*phi*delT
      ic = IdFis(k)
      l = DeplLib%IzFY(k,j)
      IF (l.LE.0) CYCLE
      l = l-ir+(ir/ic)
!      print*, 'Fission', k*1000+j, l
      DeplSysBundle%OffDiag(i,l) = DeplSysBundle%OffDiag(i,l)+YldRate
!      icbeg = DeplLib%YldMapRowPtr(ir); icend = DeplLib%YldMapRowPtr(ir+1)-1
!      DO l = icbeg, icend
!        IF (DeplLib%YldMapColIdx(l) .EQ. ic) THEN
!          ic = l - ir + (ir/ic)
!          DeplSysBundle%OffDiag(i,ic) = DeplSysBundle%OffDiag(i,ic)+YldRate
!          EXIT
!        END IF
!      END DO
    END DO
  END DO

  ! SubYld
  DO j = 1, NofRct
    DO k = 1, NofIso
      ic = k
      IF (lActinide(k)) THEN
        ir = SubYldRctAct(j)
        l = DeplLib%IzSYRAct(j,k)
      ELSE
        ir = subYldRct(j)
        l = DeplLib%IzSYR(j,k)
      END IF
      IF (l.LE.0) CYCLE
      IF(ir .EQ. 0) CYCLE
      l = l-ir+(ir/ic)
!      print*, 'SubYld', j*1000+k, l
      YldRate = xs1g(j,k)*phi*delT
      DeplSysBundle%OffDiag(i,l) = DeplSysBundle%OffDiag(i,l)+YldRate
!      icbeg = DeplLib%YldMapRowPtr(ir); icend = DeplLib%YldMapRowPtr(ir+1)-1
!      DO l = icbeg, icend
!        IF (DeplLib%YldMapColIdx(l) .EQ. ic) THEN
!          ic = l - ir + (ir/ic)
!          DeplSysBundle%OffDiag(i,ic) = DeplSysBundle%OffDiag(i,ic)+YldRate
!          EXIT
!        END IF
!      END DO
    END DO
  END DO
  DO j = 1, NofDec
    ir = SubYldDec(j)
    IF (ir .EQ. 0) CYCLE
    DO k = 1, NofIso
      ic = k
      YldRate = DecFrac(j,k)*delT
      l = DeplLib%IzSYD(j,k)
      IF (l.LE.0) CYCLE
      l=l-ir+(ir/ic)
!      print*, 'SYDec', j*1000+k, l
      DeplSysBundle%OffDiag(i,l) = DeplSysBundle%OffDiag(i,l)+YldRate
!      icbeg = DeplLib%YldMapRowPtr(ir); icend = DeplLib%YldMapRowPtr(ir+1)-1
!      DO l = icbeg, icend
!        IF (DeplLib%YldMapColIdx(l) .EQ. ic) THEN
!          ic = l - ir + (ir/ic)
!          DeplSysBundle%OffDiag(i,ic) = DeplSysBundle%OffDiag(i,ic)+YldRate
!          EXIT
!        END IF
!      END DO
    END DO
  ENDDO
!  DeplMats(i)%nnz = nnz; DeplMats(i)%lFinalized = .TRUE.
  !IF(i.EQ.9) CALL printCSR(DeplMats(i), 'CSROUT.txt',671)
  !CALL finalizeCsr(DeplMats(i), .FALSE.)
!  NULLIFY(xs1g, csrval, csrRowPtr, csrColIdx)
  DEALLOCATE(FisXs)
END DO
!$OMP END PARALLEL DO

NULLIFY(lActinide, FisFrac, DecFrac, SubYldRct, SubYldRctAct, SubYldDec)
NULLIFY(IdAftRct, IdAftDec, IdFis, IdYld)
NULLIFY(pnums, Fxrs)

END SUBROUTINE

SUBROUTINE DestroySys_wPoint(DeplSysBundle, Nvec, Solvec, lcsrT)
  USE CSRMATRIX
  USE HPDeplType
  IMPLICIT NONE
  TYPE(DeplSysBundle_Type) :: DeplSysBundle
  REAL(8), POINTER :: Nvec(:), Solvec(:)
  INTEGER :: Nsys, i
  LOGICAL :: lcsrT
  IF (lcsrT) THEN
    Nsys = DeplSysBundle%Nsys
    DO i = 1, Nsys
      CALL destroyCSR(DeplSysBundle%DeplMats(i))
    END DO
    DEALLOCATE(DeplSysBundle%DeplMats)
  ELSE
    DEALLOCATE(DeplSysBundle%Diag, DeplSysBundle%OffDiag)
  END IF
  Nvec => DeplSysBundle%pnums
  Solvec => DeplSysBundle%pnums1
END SUBROUTINE

SUBROUTINE DestroySysnVec(DeplSysBundle, lcsrT)
  USE CSRMATRIX
  USE HPDeplType
  IMPLICIT NONE
  TYPE(DeplSysBundle_Type) :: DeplSysBundle
  INTEGER :: Nsys, i
  LOGICAL :: lcsrT
  IF (lcsrT) THEN
    Nsys = DeplSysBundle%Nsys
    DO i = 1, Nsys
      CALL destroyCSR(DeplSysBundle%DeplMats(i))
    END DO
    DEALLOCATE(DeplSysBundle%DeplMats)
  ELSE
    DEALLOCATE(DeplSysBundle%Diag, DeplSysBundle%OffDiag)
  END IF
  DEALLOCATE(DeplSysBundle%pnums1, DeplSysBundle%pnums)
END SUBROUTINE

SUBROUTINE CopySolVec(DeplSysBundle, NofIso, Solvec, lFM, stream)
  USE HPDeplType
  USE CUDAFOR
  IMPLICIT NONE
  TYPE(DeplSysBundle_Type) :: DeplSysBundle
  INTEGER :: NofIso
  REAL(8) :: Solvec(:)
  LOGICAL :: lFM
  INTEGER(KIND = cuda_stream_kind) :: stream
  INTEGER :: Nsys, i, j, k
  Nsys = DeplSysBundle%Nsys
  i = cudaStreamSynchronize(stream)
  IF (lFM) THEN
    Solvec = DeplSysBundle%pnums1
  ELSE
    DO i = 1, Nsys
      k = (i-1)*NofIso
      DO j = 1, NofIso
        !WRITE(601,*) '(', j+k, DeplSysBundle%pnums1(j+k), ')'
          Solvec(j+k)=DeplSysBundle%pnums1((j-1)*Nsys+i)
!        END IF
      END DO
    END DO
  END IF
END SUBROUTINE

SUBROUTINE DestroyVecs_wCopy(DeplFxrBundle, ifxrbeg, Nsys, NofIso, Nvec, Solvec, lCorrector, lGd, lFM, stream)
  USE HPDeplType
  USE CUDAFOR
  IMPLICIT NONE
  TYPE(DeplFxrBundle_Type) :: DeplFxrBundle
  INTEGER :: ifxrbeg, Nsys, NofIso
  REAL(8), POINTER :: Nvec(:), Solvec(:)
  LOGICAL :: lCorrector, lGd, lFM
  INTEGER(KIND = cuda_stream_kind) :: stream

  INTEGER :: i, j, k
  INTEGER, POINTER :: IdTrue(:)
  TYPE(DeplFxr_Type), POINTER :: Fxrs(:)
  REAL(8), POINTER :: pnum_fxr(:)

  Fxrs => DeplFxrBundle%FxrBundle
  IF (lGd) THEN
    IdTrue => DeplFxrBundle%IdTrueGd
  ELSE
    IdTrue => DeplFxrBundle%IdTrueDepl
  END IF

  i = cudaStreamSynchronize(stream)
  DO  i = 1, Nsys
    k = IdTrue(i+ifxrbeg-1)
    IF (lCorrector) THEN
      pnum_fxr => Fxrs(k)%pnum_cor
    ELSE
      pnum_fxr => Fxrs(k)%pnum_pre
    END IF

    IF (lFM) THEN
      DO j = 1, NofIso
        pnum_fxr(j) = Solvec(j+(i-1)*NofIso)
      END DO
    ELSE
      DO j = 1, NofIso
        pnum_fxr(j) = Solvec((j-1)*Nsys+i)
      END DO
    END IF
    !print*, pnum_fxr(:)
!    IF (.NOT. lCorrector) WRITE(101,'(611E13.5)') pnum_fxr(:)
  END DO
  DEALLOCATE(Nvec, Solvec)
END SUBROUTINE
#endif
