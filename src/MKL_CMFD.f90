#include <defines.h>
#ifdef __INTEL_MKL
MODULE MKL_CMFD

USE MKL_3D

IMPLICIT NONE

LOGICAL :: lFirstCMFD = TRUE

PRIVATE
PUBLIC :: MKL_CmfdPower, MKL_CmfdDavidson

CONTAINS
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE MKL_CmfdPower(CoreInfo, CmInfo, FmInfo, eigv)

USE MKL_POWER
USE MKL_CHEBYSHEV
USE PARAM,       ONLY : TRUE, FALSE, mesg
USE TYPEDEF,     ONLY : CoreInfo_Type, CmInfo_Type, FmInfo_Type, FxrInfo_Type, PinXs_Type
USE CORE_MOD,    ONLY : GroupInfo
USE PE_MOD,      ONLY : PE
USE SUBGRP_MOD,  ONLY : FxrChiGen
USE IOUTIL,      ONLY : message
USE FILES,       ONLY : io8
USE CNTL,        ONLY : nTracerCntl
USE ITRCNTL_MOD, ONLY : ItrCntl
USE TIMER,       ONLY : nTracer_dclock, TimeChk
USE CMFD_COMMON, ONLY : HomogenizeXS, SetRadialCoupling
USE MKL_HOMOXS,  ONLY : HomogenizePnXS
USE MKL_AXIAL,   ONLY : MKL_AxialSolver, SetAxialDtil
USE HexCmfd,     ONLY : HexSetMocPhiIn

IMPLICIT NONE

TYPE (CoreInfo_Type) :: CoreInfo
TYPE (CmInfo_Type)   :: CmInfo
TYPE (FmInfo_Type)   :: FmInfo

REAL :: eigv
! ----------------------------------------------------
TYPE (FxrInfo_type), POINTER, DIMENSION(:,:) :: Fxr
TYPE (PinXS_Type),   POINTER, DIMENSION(:,:) :: PinXS

REAL, POINTER, DIMENSION(:,:,:)     :: phis, phic, AxSrc, AxPXS
REAL, POINTER, DIMENSION(:,:,:,:)   :: phim
REAL, POINTER, DIMENSION(:,:,:,:,:) :: Jout

REAL :: CmfdTimeBeg, CmfdTimeEnd, CmfdInitBeg, CmfdInitEnd, outTol, outErr, outRes, outRes0

INTEGER :: outIter, outMin, outMax, ig, iter, InIter, GrpBeg, GrpEnd, ng, nxy, nInIter, myzb, myze
INTEGER :: nGroupInfo = 2
LOGICAL :: lxslib, lscat1, lDhat, lAxRefFDM, l3dim, lGcCMFD, loutConv
LOGICAL, SAVE :: lChebyshev
! ----------------------------------------------------

CmfdTImeBeg = nTracer_dclock(FALSE, FALSE)

CALL omp_set_num_threads(PE%nCMFDThread)

Fxr   => FmInfo%Fxr
phis  => FmInfo%phis
phim  => FmInfo%phim
AxSrc => FmInfo%AxSrc
AxPXS => FmInfo%AxPXS

phic => CmInfo%phic
Jout => Cminfo%RadJout

PinXS    => mklCMFD%PinXS
outTol    = mklCntl%outerConv
outMin    = mklCntl%minOuter
outMax    = mklCntl%maxOuter
lGcCmfd   = mklCntl%lGcCMFD
lAxRefFDM = mklCntl%lAxRefFDM

myzb = mklGeom%myzb
myze = mklGeom%myze
ng   = mklGeom%ng
nxy  = mklGeom%nxy

lxslib = nTracerCntl%lxslib
lscat1 = nTracerCntl%lScat1
l3dim  = nTracerCntl%l3dim

nInIter  = 5 ! NOTICE
loutConv = FALSE
outIter  = 0

lDhat = TRUE
IF (ItrCntl%CMFDIt .EQ. ItrCntl%CMFDIt0) lDhat = FALSE

IF (lFirstCMFD) lChebyshev = mklCntl%lChebyshev
mklCntl%lChebyshev = lChebyshev .AND. .NOT. lFirstCMFD

mklCntl%lPardiso = mklCntl%lDirect .AND. (mklCntl%lDcpl .OR. PE%nCMFDproc .EQ. 1)
! ----------------------------------------------------
CmfdInitBeg = nTracer_dclock(FALSE, FALSE)

WRITE (mesg, '(A)') 'Cell Homogenization...'
IF (PE%Master) CALL message(io8, TRUE, TRUE, mesg)

CALL FxrChiGen(CoreInfo, Fxr, FmInfo, GroupInfo, nTracerCntl, myzb, myze)
CALL HomogenizeXS(CoreInfo, mklGeom%superPin, Fxr, PinXS, phis, ng, nxy, myzb, myze, lxslib, lscat1, FALSE)
CALL SetRadialCoupling(mklGeom%superPin, mklCMFD%superJout, PinXS, Jout, ng, nxy, myzb, myze, lDhat)

IF (lAxRefFDM) CALL SetReflectorDhatZero(PinXS)
IF (l3dim) CALL SetAxialDtil(mklCMFD, mklAxial)
IF (l3dim .AND. lxslib) CALL HomogenizePnXS(CoreInfo, FmInfo, GroupInfo)

CALL SetCMFDPhis(mklCMFD, PinXS, lFirstCMFD)

CmfdInitEnd = nTracer_dclock(FALSE, FALSE)
TimeChk%CmfdInitTime = TimeChk%CmfdInitTime + (CmfdInitEnd - CmfdInitBeg)
! ----------------------------------------------------
IF (mklCntl%lAxial .AND. l3dim .AND. .NOT.lFirstCMFD) THEN
  CALL GetNeighborFlux(mklCMFD)
  CALL MKL_AxialSolver(CoreInfo, CmInfo, PinXS, eigv)
END IF
! ----------------------------------------------------
CmfdInitBeg = nTracer_dclock(FALSE, FALSE)

WRITE (mesg, '(A)') 'Linear System Construction...'
IF (PE%Master) CALL message(io8, TRUE, TRUE, mesg)

IF (mklCntl%lChebyshev) CALL ChebyshevReset

CALL SetSourceOperator(mklCMFD, PinXS)

IF (mklCntl%lJacobi) THEN
  CALL SetCsrJacobiSystem(mklCMFD, PInXS, l3dim, 0.0)
ELSE
  CALL SetCsrBiCGSystem(mklCMFD, PinXS, l3dim, lFirstCMFD, 0.0)
END IF

CmfdInitEnd = nTracer_dclock(FALSE, FALSE)
TimeChk%CmfdInitTime = TimeChk%CmfdInitTime + (CmfdInitEnd - CmfdInitBeg)
! ----------------------------------------------------
WRITE (mesg, '(a)') 'Performing CMFD Calculation...'
IF (PE%Master) CALL message(io8, TRUE, TRUE, mesg)

CALL CMFDPsiUpdt(mklCMFD)
! ----------------------------------------------------
DO WHILE (.NOT. loutConv)
  ! Iter. : LS
  DO iter = 1, nGroupInfo
    IF (iter .GT. 1) THEN
      GrpBeg = GroupInfo%UpScatRange(1); GrpEnd = GroupInfo%UpScatRange(2)
    ELSE
      GrpBeg = 1; GrpEnd = ng
    END IF
    
    DO InIter = 1, nInIter
      IF (mklCntl%lDcpl) CALL GetNeighborFlux(mklCMFD)
      
      IF (mklCntl%lPardiso) THEN
        DO ig = GrpBeg, GrpEnd
          CALL CMFDSrcUpdt(mklCMFD, ig, eigv)
          CALL pardisoSolve(mklCMFD%M(ig), mklCMFD%phis(:, :, ig), mklCMFD%src(:, ig))
        END DO
      ELSE
        DO ig = GrpBeg, GrpEnd
          CALL CMFDSrcUpdt(mklCMFD, ig, eigv)
          
          IF (mklCntl%lJacobi) THEN
            CALL AAJ(mklCMFD, mklCMFD%Jacobi(ig), mklCMFD%phis(:, :, ig), mklCMFD%src(:, ig))
          ELSE
            CALL BiCGSTAB(mklCMFD, ig)
          END IF
        END DO
      END IF
    END DO
  END DO
  
  ! Updt
  CALL CMFDPsiUpdt(mklCMFD)
  CALL CMFDEigUpdt(mklCMFD, eigv)
  
  IF (mklCntl%lDcpl) CALL GetNeighborFlux(mklCMFD)
  
  ! CHK : Cnv.
  outRes = CMFDResidual(mklCMFD, eigv)
  
  IF (outIter .EQ. 0) outRes0 = outRes
  
  outErr   = outRes / outRes0
  outIter  = outIter + 1
  loutConv = (outErr .LE. outTol) .AND. (outIter .GE. outMin)
  loutConv = loutConv .OR. (outIter .GE. outMax)
  
  ItrCntl%CMFDIt = ItrCntl%CMFDIt + 1
  
  WRITE (mesg, '(a9, i9, f22.6, 3x, f10.5, 1p, e15.3)') 'MGOUTER', ItrCntl%CMFDIt, eigv, outErr, outRes
  IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
  IF (eigv.GT.2. .OR. eigv.LT.0. .OR. eigv.NE.eigv .OR. outErr.NE.outErr .OR. outRes.NE.outRes) THEN
    CALL finalize
    STOP
  END IF
  
  ! ChevyShev
  IF (lChebyshev) THEN
    CALL ChebyshevAcc(mklCMFD, mklCntl%lChebyshev)
    
    IF (mklCntl%lChebyshev) CALL CMFDPsiUpdt(mklCMFD)
  END IF
  
  !IF (mod(outIter, 5) .NE. 0) CYCLE
  IF (outIter .NE. 5) CYCLE ! Only Once
  
  ! Ax.
  IF (mklCntl%lAxial .AND. l3dim .AND. .NOT.lFirstCMFD) THEN
    CALL GetNeighborFlux(mklCMFD)
    CALL MKL_AxialSolver(CoreInfo, CmInfo, PinXS, eigv)
    
    CmfdInitBeg = nTracer_dclock(FALSE, FALSE)
    
    WRITE (mesg, '(a)') 'Linear System Construction...'
    IF (PE%Master) CALL message(io8, TRUE, TRUE, mesg)
    
    CALL SetCsrBiCGSystem(mklCMFD, PinXS, l3dim, FALSE, 0.0)
    
    CmfdInitEnd = nTracer_dclock(FALSE, FALSE)
    TimeChk%CmfdInitTime = TimeChk%CmfdInitTime + (CmfdInitEnd - CmfdInitBeg)
  END IF
  
  ! GC
  IF (lGcCMFD) THEN
    IF (mklCntl%lChebyshev) CALL ChebyshevReset
    
    CALL GetNeighborFlux(mklCMFD)
    CALL MKL_GcCmfdPower(CoreInfo, CmInfo, eigv)
  END IF
END DO
! ----------------------------------------------------
IF (mklCntl%lPardiso) THEN
  DO ig = 1, ng
    CALL pardisoDelete(mklCMFD%M(ig))
  END DO
END IF

CALL SetMOCPhis(CoreInfo, PinXS, phis, phic)

CALL HexSetMOCPhiIn(CoreInfo, mklGeom%superPin, mklCMFD%superJout, PinXS, ng, nxy, myzb, myze, ItrCntl, nTracerCntl)

IF (l3dim) THEN
  CALL GetNeighborFlux(mklCMFD)
  CALL SetAxialSrc(AxSrc, AxPXS, phic)
END IF

IF (lFirstCMFD) lFirstCMFD = FALSE

CmfdTImeEnd = nTracer_dclock(FALSE, FALSE)
TimeChk%CmfdTime = TimeChk%CmfdTime + (CmfdTimeEnd - CmfdTimeBeg)
! ----------------------------------------------------
NULLIFY (Fxr)
NULLIFY (PinXS)
NULLIFY (phis)
NULLIFY (phic)
NULLIFY (AxSrc)
NULLIFY (AxPXS)
NULLIFY (phim)
NULLIFY (Jout)
! ----------------------------------------------------

END SUBROUTINE MKL_CmfdPower
! ------------------------------------------------------------------------------------------------------------
! ======================================================================= !
!        CMFD Acceleration Employing Generalized Davidson's Method        !
!  http://www.netlib.org/utk/people/JackDongarra/etemplates/node138.html  !
! ======================================================================= !
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE MKL_CmfdDavidson(CoreInfo, CmInfo, FmInfo, eigv)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,       CmInfo_Type,        FmInfo_Type,                                    &
                           FxrInfo_Type,        PinXs_Type
USE CMFD_MOD,       ONLY : HomBuckling
USE CORE_MOD,       ONLY : GroupInfo
USE PE_MOD,         ONLY : PE
USE SUBGRP_MOD,     ONLY : FxrChiGen
USE IOUTIL,         ONLY : message
USE FILES,          ONLY : io8
USE CNTL,           ONLY : nTracerCntl
USE ITRCNTL_MOD,    ONLY : ItrCntl
USE TIMER,          ONLY : nTracer_dclock,      TimeChk
USE CMFD_COMMON,    ONLY : HomogenizeXS,        SetRadialCoupling
USE MKL_AXIAL,      ONLY : MKL_AxialSolver,     SetAxialDtil
USE MKL_DAVIDSON
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(FmInfo_Type) :: FmInfo
REAL :: eigv

TYPE(mklDavidson_Type), POINTER :: Davidson
TYPE(FxrInfo_type), POINTER :: Fxr(:, :)
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: phis(:, :, :), phic(:, :, :), phim(:,:,:,:)
REAL, POINTER :: Jout(:, :, :, :, :), AxSrc(:, :, :), AxPXS(:, :, :)
REAL :: CmfdTimeBeg, CmfdTimeEnd
REAL :: resTol = 1.0D-07, relTol = 1.0D-04
REAL :: outErr, outRes, outRes0
INTEGER :: ng, nxy
INTEGER :: myzb, myze
INTEGER :: outIter
LOGICAL :: lxslib, lscat1, lDhat, l3dim
LOGICAL :: loutConv

CmfdTImeBeg = nTracer_dclock(FALSE, FALSE)

CALL omp_set_num_threads(PE%nCMFDThread)

Davidson => mklCMFD%Davidson
Fxr => FmInfo%Fxr
phis => FmInfo%phis
phic => CmInfo%phic
Jout => Cminfo%RadJout
AxSrc => FmInfo%AxSrc
AxPXS => FmInfo%AxPXS
PinXS => mklCMFD%PinXS

myzb = mklGeom%myzb; myze = mklGeom%myze
ng = mklGeom%ng; nxy = mklGeom%nxy
lxslib = nTracerCntl%lxslib; lscat1 = nTracerCntl%lScat1; l3dim = nTracerCntl%l3dim
loutConv = FALSE; outIter = 0

lDhat = TRUE
IF (ItrCntl%Cmfdit .EQ. ItrCntl%Cmfdit0) lDhat = FALSE

WRITE(mesg,'(a)') 'Cell Homogenization...'
IF (PE%Master) CALL message(io8, TRUE, TRUE, mesg)
CALL FxrChiGen(CoreInfo, Fxr, FmInfo, GroupInfo, nTracerCntl, myzb, myze)
CALL HomogenizeXS(CoreInfo, mklGeom%superPin, Fxr, PinXS, phis, ng, nxy, myzb, myze, lxslib, lscat1, FALSE)
#ifdef Buckling
IF (nTracerCntl%lBsq) THEN
  CALL HomBuckling(CoreInfo, Fxr, phis, PinXS, myzb, myze, ng, nTracerCntl%bsq, lxsLib)
END IF
#endif
CALL SetRadialCoupling(mklGeom%superPin, mklCMFD%superJout, PinXS, Jout, ng, nxy, myzb, myze, lDhat)
IF (l3dim) CALL SetAxialDtil(mklCMFD, mklAxial)

WRITE(mesg, '(a)') 'Performing CMFD Calculation...'
IF (PE%Master) CALL message(io8, TRUE, TRUE, mesg)

IF (l3dim) CALL SetSourceOperator(mklCMFD, PinXS)
CALL SetDavidsonOperator(PinXS, lDhat, l3dim)

CALL SetCMFDPhis(mklCMFD, PinXS, lFirstCMFD)
CALL ReorderFlux(1)

IF (mklCntl%lAxial) THEN
  IF (l3dim .AND. .NOT. lFirstCMFD) THEN
    CALL GetNeighborFlux(mklCMFD)
    CALL MKL_AxialSolver(CoreInfo, CmInfo, PinXS, eigv)
  END IF
END IF

DO WHILE (.NOT. loutConv)
  outIter = outIter + 1
  ItrCntl%CMFDIt = ItrCntl%CMFDIt + 1
  CALL RayleighRitz(Davidson, outIter, eigv)
  outRes = GetResidual(Davidson, eigv)
  IF (outIter .EQ. 1) outRes0 = outRes
  outErr = outRes / outRes0
  loutConv = (outRes .LE. resTol) .OR. (outErr .LE. relTol)
  IF (mod(outIter, 5) .EQ. 0 .OR. loutConv) THEN
    WRITE(mesg, '(a10, i8, i8, f14.6, 3x, f10.5, 1p, e15.3)') 'DAVIDSON', ItrCntl%CMFDIt, outIter, eigv, outErr, outRes
    IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
  END IF
  CALL SolveDavidsonEq(Davidson, eigv)
END DO

CALL ReorderFlux(2)
CALL SetMOCPhis(CoreInfo, PinXS, phis, phic)
CALL SetMOCPhim(CoreInfo, PinXS, phim)

IF (l3dim) THEN
  CALL GetNeighborFlux(mklCMFD)
  CALL SetAxialSrc(AxSrc, AxPXS, phic)
END IF

IF (lFirstCMFD) lFirstCMFD = FALSE

CmfdTImeEnd = nTracer_dclock(FALSE, FALSE)
TimeChk%CmfdTime = TimeChk%CmfdTime + (CmfdTimeEnd - CmfdTimeBeg)

END SUBROUTINE MKL_CmfdDavidson
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE MKL_GcCmfdPower(CoreInfo, CmInfo, Keff)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,       CmInfo_Type,        PinXs_Type
USE CORE_MOD,       ONLY : GroupInfo,           GcGroupInfo
USE PE_MOD,         ONLY : PE
USE IOUTIL,         ONLY : message
USE FILES,          ONLY : io8
USE CNTL,           ONLY : nTracerCntl
USE ITRCNTL_MOD,    ONLY : ItrCntl
USE MKL_HOMOXS,     ONLY : HomogenizeGcXS,      SetRadialGcCoupling
USE MKL_AXIAL,      ONLY : SetAxialGcDtil,      SetAxialGcDhat
USE MKL_POWER
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(CmInfo_Type) :: CmInfo
REAL :: Keff

TYPE(PinXS_Type), POINTER :: PinXS(:, :), GcPinXS(:, :)
REAL :: outTol, outErr, outRes, outRes0
REAL :: eigv, seigv
INTEGER :: outIter, outMin, outMax
INTEGER :: ig, iter
INTEGER :: GrpBeg, GrpEnd, nGroupInfo
INTEGER :: ng, ngc
LOGICAL :: l3dim, loutConv, lDcpl

CALL omp_set_num_threads(PE%nCMFDThread)

PinXS => mklCMFD%PinXS
GcPinXS => mklGcCMFD%PinXS

outTol = mklCntl%outerConv
outMin = mklCntl%minOuter
outMax = mklCntl%maxOuter

l3dim = nTracerCntl%l3dim; loutConv = FALSE; lDcpl = mklCntl%lDcpl
ng = mklGeom%ng; ngc = mklGeom%ngc; nGroupInfo = 1; outIter = 0
IF (GcGroupInfo%lUpScat) nGroupInfo = 2

WRITE(mesg, '(a)') 'Group Condensing...'
IF (PE%Master) CALL message(io8, TRUE, TRUE, mesg)

IF (mklCntl%lShift) THEN
  mklCntl%lDcpl = FALSE
  seigv = 1.0 / (Keff + mklCntl%Shift)
  eigv = 1.0 / (1.0 / Keff - seigv)
  WRITE(mesg, '(a27, f4.2)') 'Wielandt Shift Parameter : ', mklCntl%Shift
  IF (PE%Master) CALL message(io8, TRUE, TRUE, mesg)
ELSE
  seigv = 0.0
  eigv = Keff
END IF

CALL HomogenizeGcXS(CoreInfo, PinXS, GcPinXS)
CALL SetCMFDPhis(mklGcCMFD, GcPinXS, TRUE)

CALL SetRadialGcCoupling(PinXS, GcPinXS)
IF (l3dim) THEN
  CALL GetNeighborFlux(mklGcCMFD)
  CALL SetAxialGcDtil(GcPinXS)
  CALL SetAxialGcDhat()
END IF

CALL SetSourceOperator(mklGcCMFD, GcPinXS)

IF (mklCntl%lJacobi) THEN
  CALL SetCsrJacobiSystem(mklGcCMFD, GcPInXS, l3dim, seigv)
ELSE
  CALL SetCsrBiCGSystem(mklGcCMFD, GcPinXS, l3dim, TRUE, seigv)
END IF

CALL CMFDPsiUpdt(mklGcCMFD)

DO WHILE (.NOT. loutConv)
  DO iter = 1, nGroupInfo
    GrpBeg = 1; GrpEnd = ngc
    IF (iter .GT. 1) THEN
      GrpBeg = GcGroupInfo%UpScatRange(1); GrpEnd = GcGroupInfo%UpScatRange(2)
    END IF
    IF (mklCntl%lDcpl) CALL GetNeighborFlux(mklGcCMFD)
    IF (mklCntl%lPardiso) THEN
      DO ig = 1, ngc
        CALL CMFDSrcUpdt(mklGcCMFD, ig, eigv)
        CALL pardisoSolve(mklGcCMFD%M(ig), mklGcCMFD%phis(:, :, ig), mklGcCMFD%src(:, ig))
      END DO
    ELSE
      DO ig = 1, ngc
        CALL CMFDSrcUpdt(mklGcCMFD, ig, eigv)
        IF (mklCntl%lJacobi) THEN
          CALL AAJ(mklGcCMFD, mklGcCMFD%Jacobi(ig), mklGcCMFD%phis(:, :, ig), mklGcCMFD%src(:, ig))
        ELSE
          CALL BiCGSTAB(mklGcCMFD, ig)
        END IF
      END DO
    END IF
  END DO
  CALL CMFDPsiUpdt(mklGcCMFD)
  CALL CMFDEigUpdt(mklGcCMFD, eigv)
  IF (mklCntl%lDcpl) CALL GetNeighborFlux(mklGcCMFD)
!  outRes = CMFDPsiErr(mklGcCMFD)
  outRes = CMFDResidual(mklGcCMFD, eigv)
  IF (outIter .EQ. 0) outRes0 = outRes
  outErr = outRes / outRes0
  outIter = outIter + 1
  loutConv = (outErr .LE. outTol) .AND. (outIter .GE. outMin)
  loutConv = loutConv .OR. (outIter .GE. outMax)
  ItrCntl%GcCMFDIt = ItrCntl%GcCMFDIt + 1
  Keff = 1.0 / (1.0 / eigv + seigv)
  WRITE(mesg, '(a9, i9, f22.6, 3x, f10.5, 1p, e15.3)') 'CGOUTER', ItrCntl%GcCMFDIt, Keff, outErr, outRes
  IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
  IF (Keff.GT.2. .OR. Keff.LT.0. .OR. Keff.NE.Keff .OR. outErr.NE.outErr .OR. outRes.NE.outRes) THEN
    CALL finalize
    STOP
  END IF
END DO

WRITE(mesg, '(a)') 'Group Reconstruction...'
IF (PE%Master) CALL message(io8, TRUE, TRUE, mesg)
CALL GcReconstruction(GcPinXS)

mklCntl%lDcpl = lDcpl

END SUBROUTINE MKL_GcCmfdPower
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetReflectorDhatZero(PinXS)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :)
INTEGER :: ipin, iz
INTEGER :: nxy, myzb, myze

nxy = mklGeom%nxy
myzb = mklGeom%myzb
myze = mklGeom%myze

DO iz = myzb, myze
  IF (.NOT. mklGeom%lRefPlane(iz)) CYCLE
  DO ipin = 1, nxy
    PinXS(ipin, iz)%Dhat = 0.0
  END DO
END DO

END SUBROUTINE SetReflectorDhatZero
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetSourceOperator(CMFD, PinXS)
USE PARAM
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
TYPE(PinXS_Type), POINTER :: PinXS(:, :)

INTEGER, POINTER :: pinMap(:), planeMap(:)
INTEGER :: ir, iz, izf, ig, igf, igt, ipin, ipin_map
INTEGER :: gb, ge
INTEGER :: ng, nxy, nzCMFD
REAL, POINTER :: PinVolFm(:, :)
REAL :: val

ng = CMFD%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
planeMap => CMFD%planeMap
pinMap => mklGeom%pinMap
PinVolFm => mklGeom%PinVolFm

CMFD%S = 0.0

!$OMP PARALLEL PRIVATE(iz, ir, ipin_map, val)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO igt = 1, ng
  DO igf = 1, ng
    DO izf = 1, nzCMFD
      iz = planeMap(izf)
      DO ipin = 1, nxy
        ir = ipin + (izf - 1) * nxy
        ipin_map = pinMap(ipin)
        IF (PinXS(ipin_map, iz)%XSs(igt)%ib .GT. igf) CYCLE
        IF (PinXS(ipin_map, iz)%XSs(igt)%ie .LT. igf) CYCLE
        IF (igt .EQ. igf) THEN
          val = PinXS(ipin_map, iz)%XSs(igt)%self * PinVolFm(ipin, izf)
        ELSE
          val = PinXS(ipin_map, iz)%XSs(igt)%from(igf) * PinVolFm(ipin, izf)
        END IF
        CMFD%S(ir, igf, igt) = val
      END DO
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

CMFD%F = 0.0

!$OMP PARALLEL PRIVATE(iz, ir, ipin_map)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ig = 1, ng
  DO izf = 1, nzCMFD
    iz = planeMap(izf)
    DO ipin = 1, nxy
      ir = ipin + (izf - 1) * nxy
      ipin_map = pinMap(ipin)
      CMFD%F(ir, ig) = PinXS(ipin_map, iz)%XSnf(ig) * PinVolFm(ipin, izf)
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

CMFD%Chi = 0.0

!$OMP PARALLEL PRIVATE(iz, ir, ipin_map)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO ig = 1, ng
  DO izf = 1, nzCMFD
    iz = planeMap(izf)
    DO ipin = 1, nxy
      ir = ipin + (izf - 1) * nxy
      ipin_map = pinMap(ipin)
      CMFD%Chi(ir, ig) = PinXS(ipin_map, iz)%Chi(ig)
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE SetSourceOperator
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetCMFDPhis(CMFD, PinXS, lCopy)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
LOGICAL :: lCopy

INTEGER :: ig, iz, idx, izf, ipin, ipin_map
INTEGER :: ng, nxy, nzCMFD
INTEGER :: myzb, myze
INTEGER, POINTER :: pinMap(:), planeMap(:), fmRange(:, :)
REAL :: fmult
REAL, POINTER :: hzfm(:), hz(:)

ng = CMFD%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
myzb = mklGeom%myzb
myze = mklGeom%myze
pinMap => mklGeom%pinMap
planeMap => CMFD%planeMap
fmRange => mklGeom%fmRange
hzfm => mklGeom%hzfm
hz => mklGeom%hz

IF (lCopy) THEN
  !$OMP PARALLEL PRIVATE(iz, ipin_map)
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
  DO ig = 1, ng
    DO izf = 1, nzCMFD
      iz = planeMap(izf)
      DO ipin = 1, nxy
        ipin_map = pinMap(ipin)
        CMFD%phis(ipin, izf, ig) = PinXS(ipin_map, iz)%Phi(ig)
      END DO
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
ELSE
  !$OMP PARALLEL PRIVATE(ipin_map, fmult)
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
  DO ig = 1, ng
    DO ipin = 1, nxy
      ipin_map = pinMap(ipin)
      DO iz = myzb, myze
        fmult = PinXS(ipin_map, iz)%Phi(ig) / CMFD%phic(ipin, iz, ig)
        DO izf = fmRange(iz, 1), fmRange(iz, 2)
          CMFD%phis(ipin, izf, ig) = CMFD%phis(ipin, izf, ig) * fmult
        END DO
      END DO
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
END IF

END SUBROUTINE SetCMFDPhis
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE GcReconstruction(GcPinXS)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: GcPinXS(:, :)

INTEGER :: ng, ngc, nxy, nzCMFD
INTEGER :: ig, igc, ipin, ipin_map, iz
INTEGER, POINTER :: pinMap(:)
REAL :: fmult

ng = mklGeom%ng
ngc = mklGeom%ngc
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
pinMap => mklGeom%pinMap

!$OMP PARALLEL PRIVATE(ipin_map, igc, fmult)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO iz = 1, nzCMFD
    DO ipin = 1, nxy
      ipin_map = pinMap(ipin)
      igc = mklGeom%GcStructInv(ig)
      fmult = mklGcCMFD%phis(ipin, iz, igc) / GcPinXS(ipin_map, iz)%Phi(igc)
      mklCMFD%phis(ipin, iz, ig) = mklCMFD%phis(ipin, iz, ig) * fmult
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

CALL dcopy(nxy * nzCMFD, mklGcCMFD%psi, 1, mklCMFD%psi, 1)

END SUBROUTINE GcReconstruction
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetMOCPhis(CoreInfo, PinXS, phis, phic)
USE TYPEDEF,        ONLY : CoreInfo_Type,       PinXS_Type,         Pin_Type,       Cell_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: phis(:, :, :), phic(:, :, :)

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(superPin_Type), POINTER :: superPin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
INTEGER :: i, j, ig, iz, izf, ifsr, icel, ixy, ixy_map, ipin
INTEGER :: FsrIdxSt
INTEGER :: ng, nxy, nzCMFD, nLocalFsr
INTEGER :: myzb, myze
INTEGER, POINTER :: pinMap(:), fmRange(:, :)
REAL :: fmult, phisum
REAL, POINTER :: hzfm(:), hz(:)

Pin => CoreInfo%Pin
Cell => CoreInfo%CellInfo
superPin => mklGeom%superPin
pinMap => mklGeom%pinMap
ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
myzb = mklGeom%myzb
myze = mklGeom%myze
fmRange => mklGeom%fmRange
hzfm => mklGeom%hzfm
hz => mklGeom%hz

!$OMP PARALLEL PRIVATE(FsrIdxSt, ipin, icel, ixy_map, nLocalFsr, phisum, fmult, ifsr)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO iz = myzb, myze
    DO ixy = 1, nxy
      phisum = 0.0
      ixy_map = pinMap(ixy)
      DO izf = fmRange(iz, 1), fmRange(iz, 2)
        phisum = phisum + mklCMFD%phis(ixy, izf, ig) * (hzfm(izf) / hz(iz))
      END DO
      fmult = phisum / PinXS(ixy_map, iz)%Phi(ig)
      DO j = 1, superPin(ixy_map)%nxy
        ipin = superPin(ixy_map)%pin(j)
        FsrIdxSt = Pin(ipin)%FsrIdxSt
        icel = Pin(ipin)%Cell(iz)
        nLocalFsr = Cell(icel)%nFsr
        DO i = 1, nLocalFsr
          ifsr = FsrIdxSt + i - 1
          phis(ifsr, iz, ig) = phis(ifsr, iz, ig) * fmult
        END DO
        phic(ipin, iz, ig) = phisum
      END DO
      mklCMFD%phic(ixy, iz, ig) = phisum
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE SetMOCPhis
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetMOCPhim(CoreInfo, PinXS, phim)

USE TYPEDEF, ONLY : CoreInfo_Type, PinXS_Type, Pin_Type, Cell_Type
USE CNTL,    ONLY : nTracerCntl

IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo

TYPE(PinXS_Type), POINTER, DIMENSION(:,:) :: PinXS

REAL, POINTER, DIMENSION(:,:,:,:) :: phim

TYPE(Pin_Type),      POINTER, DIMENSION(:) :: Pin
TYPE(superPin_Type), POINTER, DIMENSION(:) :: superPin
TYPE(Cell_Type),     POINTER, DIMENSION(:) :: Cell

INTEGER :: ig, iz, izf, ifsr, jfsr, icel, ixy, ixy_map, ipin, jpin
INTEGER :: FsrIdxSt, ng, nxy, nLocalFsr, myzb, myze

INTEGER, POINTER, DIMENSION(:) :: pinMap 

REAL :: fmult
! ----------------------------------------------------

IF (.NOT. nTracerCntl%lScat1) RETURN

Pin  => CoreInfo%Pin
Cell => CoreInfo%CellInfo

superPin => mklGeom%superPin
pinMap   => mklGeom%pinMap
ng        = mklGeom%ng
nxy       = mklGeom%nxy
myzb      = mklGeom%myzb
myze      = mklGeom%myze

!$OMP PARALLEL PRIVATE(FsrIdxSt, ipin, jpin, icel, ixy_map, nLocalFsr, fmult, ifsr, jfsr)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO iz = myzb, myze
    DO ixy = 1, nxy
      ixy_map = pinMap(ixy)
      fmult   = mklCMFD%phic(ixy, iz, ig) / PinXS(ixy_map, iz)%Phi(ig)
      
      DO ipin = 1, superPin(ixy_map)%nxy
        jpin      = superPin(ixy_map)%pin(ipin)
        FsrIdxSt  = Pin(jpin)%FsrIdxSt
        icel      = Pin(jpin)%Cell(iz)
        nLocalFsr = Cell(icel)%nFsr
        
        DO ifsr = 1, nLocalFsr
          jfsr = FsrIdxSt + ifsr - 1
          
          phim(:, ig, jfsr, iz) = phim(:, ig, jfsr, iz) * fmult
        END DO
      END DO
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
! ----------------------------------------------------

END SUBROUTINE SetMOCPhim
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetAxialSrc(AxSrc, AxPXS, phic)
USE CNTL,           ONLY : nTracerCntl
IMPLICIT NONE

REAL, POINTER :: AxSrc(:, :, :), AxPXS(:, :, :), phic(:, :, :)

TYPE(superPin_Type), POINTER :: Pin(:)
INTEGER :: i, ig, iz, izf, ipin, ixy, ixy_map
INTEGER :: ng, nxy
INTEGER :: myzb, myze
INTEGER, POINTER :: pinMap(:), fmRange(:, :)
REAL :: myphi, neighphi, Dtil, Dhat, Jnet
REAL, POINTER :: hz(:)

Pin => mklGeom%superPin
ng = mklGeom%ng
nxy = mklGeom%nxy
myzb = mklGeom%myzb
myze = mklGeom%myze
pinMap => mklGeom%pinMap
fmRange => mklGeom%fmRange
hz => mklGeom%hz

!$OMP PARALLEL PRIVATE(ixy_map, ipin, izf, myphi, neighphi, Dtil, Dhat, Jnet)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO iz = myzb, myze
    DO ixy = 1, nxy
      ixy_map = pinMap(ixy)
      !--- Axial Source from Bottom
      izf = fmRange(iz, bottom)
      myphi = mklCMFD%phis(ixy, izf, ig)
      IF (iz .EQ. myzb) THEN
        neighphi = mklCMFD%neighphis(ixy, ig, bottom)
      ELSE
        neighphi = mklCMFD%phis(ixy, izf - 1, ig)
      END IF
      Dtil = mklCMFD%AxDtil(bottom, ixy, izf, ig)
      Dhat = mklCMFD%AxDhat(bottom, ixy, izf, ig)
      Jnet = - Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
      DO i = 1, Pin(ixy_map)%nxy
        ipin = Pin(ixy_map)%pin(i)
        AxSrc(ipin, iz, ig) = Jnet
      END DO
      !--- Axial Source from Top
      izf = fmRange(iz, top)
      myphi = mklCMFD%phis(ixy, izf, ig)
      IF (iz .EQ. myze) THEN
        neighphi = mklCMFD%neighphis(ixy, ig, top)
      ELSE
        neighphi = mklCMFD%phis(ixy, izf + 1, ig)
      END IF
      Dtil = mklCMFD%AxDtil(top, ixy, izf, ig)
      Dhat = mklCMFD%AxDhat(top, ixy, izf, ig)
      Jnet = - Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
      DO i = 1, Pin(ixy_map)%nxy
        ipin = Pin(ixy_map)%pin(i)
        AxSrc(ipin, iz, ig) = AxSrc(ipin, iz, ig) + Jnet
      END DO
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(ipin)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO iz = myzb, myze
    DO ixy = 1, nxy
      DO i = 1, Pin(ixy)%nxy
        ipin = Pin(ixy)%pin(i)
        AxSrc(ipin, iz, ig) = AxSrc(ipin, iz, ig) / hz(iz)
      END DO
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

AxPXS(:, myzb : myze, :) = 0.0

IF (nTracerCntl%LkgSplitLv .EQ. 0) THEN
  !$OMP PARALLEL PRIVATE(ipin)
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
  DO ig = 1, ng
    DO iz = myzb, myze
      DO ixy = 1, nxy
        DO i = 1, Pin(ixy)%nxy
          ipin = Pin(ixy)%pin(i)
          IF (AxSrc(ipin, iz, ig) .LT. 0.0) CYCLE
          IF (phic(ipin, iz, ig) .LT. 0.0) CYCLE
          AxPXS(ipin, iz, ig) = AxSrc(ipin, iz, ig) / phic(ipin, iz, ig)
        END DO
      END DO
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
END IF

END SUBROUTINE SetAxialSrc
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE GetNeighborFlux(CMFD)

IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
INTEGER :: ig, ng, nxy, nzCMFD
REAL, POINTER :: myphis(:, :, :), neighphis(:, :, :)

ng = CMFD%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

ALLOCATE(myphis(nxy, ng, 2))
neighphis => CMFD%neighphis

DO ig = 1, ng
  CALL dcopy(nxy, CMFD%phis(:, 1, ig), 1, myphis(:, ig, bottom), 1)
  CALL dcopy(nxy, CMFD%phis(:, nzCMFD, ig), 1, myphis(:, ig, top) , 1)
END DO

CALL InitFastComm()
CALL GetNeighborFast(ng * nxy, myphis(:, :, bottom), neighphis(:, :, top), bottom)
CALL GetNeighborFast(ng * nxy, myphis(:, :, top), neighphis(:, :, bottom), top)
CALL FinalizeFastComm()

IF (mklGeom%lBottom) THEN
  IF (mklGeom%AxBC(bottom) .EQ. VoidCell) neighphis(:, :, bottom) = 0.0
  IF (mklGeom%AxBC(bottom) .EQ. RefCell) CALL dcopy(ng * nxy, myphis(:, :, bottom), 1, neighphis(:, :, bottom), 1)
END IF

IF (mklGeom%lTop) THEN
  IF (mklGeom%AxBC(top) .EQ. VoidCell) neighphis(:, :, top) = 0.0
  IF (mklGeom%AxBC(top) .EQ. RefCell) CALL dcopy(ng * nxy, myphis(:, :, top), 1, neighphis(:, :, top), 1)
END IF

DEALLOCATE(myphis)

END SUBROUTINE GetNeighborFlux
! ------------------------------------------------------------------------------------------------------------

END MODULE MKL_CMFD
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE MKL_ReorderPinXS(CoreInfo, CmInfo)
USE TYPEDEF,        ONLY : CoreInfo_Type,       CmInfo_Type,         PinXS_Type
USE CMFD_COMMON
USE MKL_3D
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(CmInfo_Type) :: CmInfo

TYPE(PinXS_Type), POINTER :: PinXS(:, :)
TYPE(superPin_Type), POINTER :: Pin(:)
REAL, POINTER :: PinVol(:, :)
INTEGER :: ng, nxy, myzb, myze
INTEGER :: i, ipin, iFuelPin, ixy, iz

PinXS => CmInfo%PinXS
Pin => mklGeom%superPin
PinVol => CoreInfo%PinVol
ng = mklGeom%ng
nxy = mklGeom%nxy
myzb = mklGeom%myzb
myze = mklGeom%myze

DO iz = myzb, myze
  DO ixy = 1, nxy
    DO i = 1, Pin(ixy)%nxy
      ipin = Pin(ixy)%pin(i)
      CALL CopyPinXS(mklCMFD%PinXS(ixy, iz), PinXS(ipin, iz), ng)
    END DO
  END DO
END DO

DO iz = myzb, myze
  DO ixy = 1, nxy
    IF (.NOT. Pin(ixy)%lFuel(iz)) CYCLE
    iFuelPin = Pin(ixy)%iFuelPin
    PinXS(iFuelPin, iz)%XSnf = PinXS(iFuelPin, iz)%XSnf * PinVol(iFuelPin, iz)
    PinXS(iFuelPin, iz)%XSkf = PinXS(iFuelPin, iz)%XSkf * PinVol(iFuelPin, iz)
    DO i = 1, Pin(ixy)%nxy
      ipin = Pin(ixy)%pin(i); IF (ipin .EQ. iFuelPin) CYCLE
      PinXS(iFuelPin, iz)%XSnf = PinXS(iFuelPin, iz)%XSnf + PinXS(ipin, iz)%XSnf * PinVol(ipin, iz)
      PinXS(iFuelPin, iz)%XSkf = PinXS(iFuelPin, iz)%XSkf + PinXS(ipin, iz)%XSkf * PinVol(ipin, iz)
      PinXS(ipin, iz)%XSnf = 0.0; PinXS(ipin, iz)%XSkf = 0.0
    END DO
    PinXS(iFuelPin, iz)%XSnf = PinXS(iFuelPin, iz)%XSnf / PinVol(iFuelPin, iz)
    PinXS(iFuelPin, iz)%XSkf = PinXS(iFuelPin, iz)%XSkf / PinVol(iFuelPin, iz)
  END DO
END DO

END SUBROUTINE MKL_ReorderPinXS
#endif