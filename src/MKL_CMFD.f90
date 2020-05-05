#include <defines.h>
!--- CNJ Edit : 3D CMFD Acceleration Modules with Intel MKL 
#ifdef __INTEL_MKL

MODULE MKL_HOMOXS

USE MKL_3D
IMPLICIT NONE

CONTAINS

SUBROUTINE HomogenizePnXS(CoreInfo, FmInfo, GroupInfo)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,       FmInfo_Type,        GroupInfo_Type,     FxrInfo_Type,               &
                           Pin_Type,            XsMac_Type
USE MacXsLib_Mod,   ONLY : MacP1XsScatMatrix,   MacP2XsScatMatrix,  MacP3XsScatMatrix
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo

TYPE(superPin_Type), POINTER :: superPin(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(XsMac_Type) :: XsMac
INTEGER :: ig, ifxr, ixy, ixy_map, ipin, iz, izf
INTEGER :: ng, nxy, nzCMFD
INTEGER, POINTER :: pinMap(:), planeMap(:)

Pin => CoreInfo%Pin
Fxr => FmInfo%Fxr
superPin => mklGeom%superPin
ng = mklGeom%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
planeMap => mklGeom%planeMap
pinMap => mklGeom%pinMap

DO izf = 1, nzCMFD
  iz = planeMap(izf)
  DO ixy = 1, nxy
    IF (.NOT. mklGeom%lH2OCell(izf, ixy)) CYCLE
    ixy_map = pinMap(ixy)
    ipin = superPin(ixy_map)%pin(1)
    ifxr = Pin(ipin)%FxrIdxSt
    CALL MacP1XsScatMatrix(XsMac, Fxr(ifxr, iz), 1, ng, ng, GroupInfo)
!    CALL MacP2XsScatMatrix(XsMac, Fxr(ifxr, iz), 1, ng, ng, GroupInfo)
!    CALL MacP3XsScatMatrix(XsMac, Fxr(ifxr, iz), 1, ng, ng, GroupInfo)
    mklAxial%SmP1(:, :, izf, ixy) = XsMac%XsMacP1Sm
!    mklAxial%SmP2(:, :, izf, ixy) = XsMac%XsMacP2Sm
!    mklAxial%SmP3(:, :, izf, ixy) = XsMac%XsMacP3Sm
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE HomogenizeGcXS(CoreInfo, PinXS, GcPinXS)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,       PinXS_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(PinXS_Type), POINTER :: PinXS(:, :), GcPinXS(:, :)

REAL :: localphis(mklGeom%ng)
INTEGER :: nxy, nzCMFD
INTEGER :: ipin, ipin_map, iz, izf
INTEGER, POINTER :: pinMap(:), planeMap(:)

nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
pinMap => mklGeom%pinMap
planeMap => mklGeom%planeMap

DO izf = 1, nzCMFD
  iz = planeMap(izf)
  !$OMP PARALLEL PRIVATE(localphis, ipin_map)
  !$OMP DO SCHEDULE(GUIDED)
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    localphis = mklCMFD%phis(ipin, izf, :)
    CALL HomogenizeCellGcXS(PinXS(ipin_map, iz), GcPinXS(ipin_map, izf), localphis)
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDDO

END SUBROUTINE

SUBROUTINE HomogenizeCellGcXS(PinXS, GcPinXS, phis)
USE PARAM
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type) :: PinXS, GcPinXS
REAL :: phis(:)

REAL :: RR(3), RRS(mklGeom%ngc, mklGeom%ngc)
REAL :: phisum, chisum
INTEGER :: ng, ngc
INTEGER :: igc, ig, igb, ige
INTEGER :: igs, igsb, igse, ig0

ng = mklGeom%ng
ngc = mklGeom%ngc

DO igc = 1, ngc
  igb = mklGeom%GcStruct(1, igc); ige = mklGeom%GcStruct(2, igc)
  RR = 0; phisum = 0; chisum = 0
  DO ig = igb, ige
    RR(1) = RR(1) + PinXS%XStr(ig) * phis(ig)
    RR(2) = RR(2) + PinXS%XSnf(ig) * phis(ig)
    RR(3) = RR(3) + PinXS%XSD(ig) * phis(ig)
    phisum = phisum + phis(ig)
    chisum = chisum + PinXS%Chi(ig)
  ENDDO
  RR = RR / phisum
  GcPinXS%XStr(igc) = RR(1)
  GcPinXS%XSnf(igc) = RR(2)
  GcPinXS%XSD(igc) = RR(3)
  GcPinXS%Chi(igc) = chisum
  GcPinXS%Phi(igc) = phisum
ENDDO

RRS = 0
DO ig = 1, ng
  igc = mklGeom%GcStructInv(ig)
  igsb = PinXS%XSs(ig)%ib; igse = PinXS%XSs(ig)%ie
  ig0 = igc
  RRS(ig0, igc) = RRS(ig0, igc) + PinXS%XSs(ig)%WithInGroupScat * phis(ig)
  DO igs = igsb, igse
    ig0 = mklGeom%GcStructInv(igs)
    RRS(ig0, igc) = RRS(ig0, igc) + PinXS%XSs(ig)%from(igs) * phis(igs)
  ENDDO
ENDDO

DO igc = 1, ngc
  RRS(igc, :) = RRS(igc, :) / GcPinXS%Phi(igc)
ENDDO

DO ig = 1, ngc
  igsb = GcPinXS%XSs(ig)%ib; igse = GcPinXS%XSs(ig)%ie
  DO igs = igsb, igse
    GcPinXS%XSs(ig)%from(igs) = RRS(igs, ig)
  ENDDO
  GcPinXS%XSs(ig)%WithInGroupScat = GcPinXS%XSs(ig)%from(ig)
  GcPinXS%XSs(ig)%from(ig) = 0.0
ENDDO

DO ig = 1, ngc
  GcPinXS%XSr(ig) = GcPinXS%XStr(ig) - GcPinXS%XSs(ig)%WithInGroupScat
ENDDO

END SUBROUTINE

SUBROUTINE SetRadialGcCoupling(PinXS, GcPinXS)
USE PARAM
USE TYPEDEF,        ONLY : PinXS_Type
USE CNTL,           ONLY : nTracerCntl
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :), GcPinXS(:, :)

TYPE(superPin_Type), POINTER :: Pin(:)
INTEGER, POINTER :: pinMap(:), pinMapRev(:), planeMap(:)
INTEGER :: ng, ngc, nxy, nzCMFD
INTEGER :: ig, igc, ipin, ipin_map, ineighpin, iz, izf, ibd, inbd
REAL :: Dtil, Dhat, pDhat(2), atil, myphi, neighphi, mybeta, neighbeta, albedo, jfdm, surfphifdm, smy
REAL, POINTER :: Jnet(:, :, :, :), Jpart(:, :, :, :, :)

IF (nTracerCntl%lHex) THEN
  CALL HexSetRadialGcCoupling(PinXs, GcPinXs)
  
  RETURN
END IF

Pin => mklGeom%superPin
ng = mklGeom%ng
ngc = mklGeom%ngc
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
pinMap => mklGeom%pinMap
pinMapRev => mklGeom%pinMapRev
planeMap => mklGeom%planeMap

ALLOCATE(Jnet(4, nxy, nzCMFD, ngc)); Jnet = 0.0

!--- Condense Currents

DO igc = 1, ngc
  DO ig = mklGeom%GcStruct(1, igc), mklGeom%GcStruct(2, igc)
    DO izf = 1, nzCMFD
      iz = planeMap(izf)
      !$OMP PARALLEL PRIVATE(ipin_map, ineighpin, myphi, neighphi, Dtil, Dhat, pDhat, jfdm)
      !$OMP DO SCHEDULE(GUIDED)
      DO ipin = 1, nxy
        ipin_map = pinMap(ipin)
        myphi = mklCMFD%phis(ipin, izf, ig)
        DO ibd = 1, 4
          ineighpin = Pin(ipin_map)%Neighidx(ibd)
          ineighpin = pinMapRev(ineighpin)
          IF (ineighpin .LE. 0) THEN
            neighphi = 0.0
          ELSE
            neighphi = mklCMFD%phis(ineighpin, izf, ig)
          ENDIF
          Dtil = PinXS(ipin_map, iz)%Dtil(ibd, ig)
          Dhat = PinXS(ipin_map, iz)%Dhat(ibd, ig)
          jfdm = - Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
          Jnet(ibd, ipin, izf, igc) = Jnet(ibd, ipin, izf, igc) + jfdm
        ENDDO
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL
    ENDDO
  ENDDO
ENDDO

!--- Compute Coupling Coefficients

!$OMP PARALLEL PRIVATE(ipin_map, ineighpin, inbd, myphi, neighphi, mybeta, neighbeta, Dtil, Dhat, pDhat, jfdm, albedo, smy)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO igc = 1, ngc
  DO izf = 1, nzCMFD
    DO ipin = 1, nxy
      ipin_map = pinMap(ipin)
      DO ibd = 1, 4
        ineighpin = Pin(ipin_map)%NeighIdx(ibd)
        smy = Pin(ipin_map)%BdLength(ibd)
        myphi = GcPinXS(ipin_map, izf)%Phi(igc)
        mybeta = GcPinXS(ipin_map, izf)%XSD(igc) / Pin(ipin_map)%Center2SurfaceL(ibd)
        IF (ineighpin .GT. 0) THEN
          inbd = Pin(ineighpin)%NeighSurfIdx(ibd)
          neighphi = GcPinXS(ineighpin, izf)%Phi(igc)
          neighbeta = GcPinXS(ineighpin, izf)%XSD(igc) / Pin(ineighpin)%Center2SurfaceL(inbd)
          Dtil = mybeta * neighbeta / (mybeta + neighbeta) * smy
          jfdm = - Dtil * (neighphi - myphi)
          Dhat = - (Jnet(ibd, ipin, izf, igc) - jfdm) / (myphi + neighphi)
        ELSE
          IF (ineighpin .EQ. VoidCell) THEN
            neighbeta = 0.5; neighphi = 0.0; albedo = 0.5
          ELSEIF (ineighpin .EQ. RefCell) THEN
            neighbeta = 0.0; neighphi = myphi; albedo = 0.0
          ENDIF
          Dtil = mybeta * neighbeta / (mybeta + neighbeta) * smy
          jfdm = - Dtil * (neighphi - myphi)
          Dhat = - (Jnet(ibd, ipin, izf, igc) - jfdm) / (myphi + neighphi)
        ENDIF 
        GcPinXS(ipin_map, izf)%Dtil(ibd, igc) = Dtil
        GcPinXS(ipin_map, izf)%Dhat(ibd, igc) = Dhat
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

DEALLOCATE(Jnet)

END SUBROUTINE

SUBROUTINE HexSetRadialGcCoupling(PinXS, GcPinXS)

USE PARAM
USE allocs
USE geom,    ONLY : ncbd
USE TYPEDEF, ONLY : PinXS_Type

IMPLICIT NONE

TYPE(PinXS_Type),    POINTER :: PinXS(:, :), GcPinXS(:, :)
TYPE(superPin_Type), POINTER :: Pin(:)

INTEGER, POINTER :: pinMap(:), pinMapRev(:), planeMap(:)
INTEGER :: ng, ngc, nxy, nzCMFD
INTEGER :: ig, igc, ipin, ipin_map, ineighpin, iz, izf, ibd, jbd, iNgh, jNgh
REAL :: Dtil, Dhat, myphi, neighphi, mybeta, neighbeta, albedo, jfdm, smy
REAL, POINTER :: Jnet(:, :, :, :)
! ----------------------------------------------------

ng     = mklGeom%ng
ngc    = mklGeom%ngc
nxy    = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

Pin       => mklGeom%superPin
pinMap    => mklGeom%pinMap
pinMapRev => mklGeom%pinMapRev
planeMap  => mklGeom%planeMap

CALL dmalloc(Jnet, ncbd, nxy, nzCMFD, ngc)
! ----------------------------------------------------
!               01. CONDENSE : current
! ----------------------------------------------------
DO igc = 1, ngc
  DO ig = mklGeom%GcStruct(1, igc), mklGeom%GcStruct(2, igc)
    DO izf = 1, nzCMFD
      iz = planeMap(izf)
      !$OMP PARALLEL PRIVATE(ipin_map, ineighpin, myphi, neighphi, Dtil, Dhat, jfdm, iNgh)
      !$OMP DO SCHEDULE(GUIDED)
      DO ipin = 1, nxy
        ipin_map = pinMap(ipin)
        myphi    = mklCMFD%phis(ipin, izf, ig)
        
        DO iNgh = 1, Pin(ipin_map)%nNgh
          ineighpin = Pin(ipin_map)%Neighidx(iNgh)
          ineighpin = pinMapRev(ineighpin)
          
          IF (ineighpin .LE. 0) THEN
            neighphi = ZERO
          ELSE
            neighphi = mklCMFD%phis(ineighpin, izf, ig)
          END IF
          
          Dtil = PinXS(ipin_map, iz)%Dtil(iNgh, ig)
          Dhat = PinXS(ipin_map, iz)%Dhat(iNgh, ig)
          jfdm = - Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
          
          Jnet(iNgh, ipin, izf, igc) = Jnet(iNgh, ipin, izf, igc) + jfdm
        END DO
      END DO
      !$OMP END DO
      !$OMP END PARALLEL
    END DO
  END DO
END DO
! ----------------------------------------------------
!               02. CAL : dhat, dtil
! ----------------------------------------------------
!$OMP PARALLEL PRIVATE(ipin_map, ineighpin, jbd, myphi, neighphi, mybeta, neighbeta, Dtil, Dhat, jfdm, albedo, smy, iNgh, jNgh)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO igc = 1, ngc
  DO izf = 1, nzCMFD
    DO ipin = 1, nxy
      ipin_map = pinMap(ipin)
      
      DO iNgh = 1, Pin(ipin_map)%nNgh
        ibd       = Pin(ipin_map)%NghBd(iNgh)
        ineighpin = Pin(ipin_map)%NeighIdx(iNgh)
        jNgh      = Pin(ipin_map)%NeighSurfIdx(iNgh)
        smy       = Pin(ipin_map)%BdLength(iNgh)
        
        myphi  = GcPinXS(ipin_map, izf)%Phi(igc)
        mybeta = GcPinXS(ipin_map, izf)%XSD(igc) / Pin(ipin_map)%Center2SurfaceL(ibd)
        
        IF (ineighpin .GT. 0) THEN
          jbd = Pin(ineighpin)%NghBd(jNgh)
          
          neighphi  = GcPinXS(ineighpin, izf)%Phi(igc)
          neighbeta = GcPinXS(ineighpin, izf)%XSD(igc) / Pin(ineighpin)%Center2SurfaceL(jbd)
          
          Dtil = mybeta * neighbeta / (mybeta + neighbeta) * smy
          jfdm = - Dtil * (neighphi - myphi)
          Dhat = - (Jnet(iNgh, ipin, izf, igc) - jfdm) / (myphi + neighphi)
        ELSE
          IF (ineighpin .EQ. VoidCell) THEN
            neighbeta = HALF
            neighphi  = ZERO
            albedo    = HALF
          ELSE IF (ineighpin .EQ. RefCell) THEN
            neighbeta = ZERO
            neighphi  = myphi
            albedo    = ZERO
          END IF
          Dtil = mybeta * neighbeta / (mybeta + neighbeta) * smy
          jfdm = - Dtil * (neighphi - myphi)
          Dhat = - (Jnet(iNgh, ipin, izf, igc) - jfdm) / (myphi + neighphi)
        END IF
        
        GcPinXS(ipin_map, izf)%Dtil(iNgh, igc) = Dtil
        GcPinXS(ipin_map, izf)%Dhat(iNgh, igc) = Dhat
      END DO
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

DEALLOCATE(Jnet)
! ----------------------------------------------------

END SUBROUTINE HexSetRadialGcCoupling

END MODULE
    
!--- Driver & General Routines --------------------------------------------------------------------
    
MODULE MKL_CMFD

USE MKL_3D
IMPLICIT NONE

LOGICAL :: lFirstCMFD = TRUE

PRIVATE
PUBLIC :: MKL_CmfdPower, MKL_CmfdDavidson

CONTAINS

!--- Public Routines ------------------------------------------------------------------------------

SUBROUTINE MKL_CmfdPower(CoreInfo, CmInfo, FmInfo, eigv)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,       CmInfo_Type,        FmInfo_Type,                                    &
                           FxrInfo_Type,        PinXs_Type
USE CORE_MOD,       ONLY : GroupInfo
USE PE_MOD,         ONLY : PE
USE SUBGRP_MOD,     ONLY : FxrChiGen
USE IOUTIL,         ONLY : message
USE FILES,          ONLY : io8
USE CNTL,           ONLY : nTracerCntl
USE ITRCNTL_MOD,    ONLY : ItrCntl
USE TIMER,          ONLY : nTracer_dclock,      TimeChk
USE CMFD_COMMON,    ONLY : HomogenizeXS,        SetRadialCoupling
USE MKL_HOMOXS,     ONLY : HomogenizePnXS
USE MKL_AXIAL,      ONLY : MKL_AxialSolver,     SetAxialDtil
USE MKL_POWER
USE MKL_CHEBYSHEV
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(FmInfo_Type) :: FmInfo
REAL :: eigv

TYPE(FxrInfo_type), POINTER :: Fxr(:, :)
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: phis(:, :, :), phim(:, :, :, :), phic(:, :, :)
REAL, POINTER :: Jout(:, :, :, :, :), AxSrc(:, :, :), AxPXS(:, :, :)
REAL :: CmfdTimeBeg, CmfdTimeEnd, CmfdInitBeg, CmfdInitEnd
REAL :: outTol, outErr, outRes, outRes0
INTEGER :: outIter, outMin, outMax
INTEGER :: ig, iter, InIter
INTEGER :: GrpBeg, GrpEnd, nGroupInfo = 2
INTEGER :: ng, nxy, nInIter
INTEGER :: myzb, myze
LOGICAL :: lxslib, lscat1, lDhat, lAxRefFDM, l3dim, lGcCMFD, loutConv
LOGICAL, SAVE :: lChebyshev

CmfdTImeBeg = nTracer_dclock(FALSE, FALSE)

CALL omp_set_num_threads(PE%nCMFDThread)

Fxr => FmInfo%Fxr
phis => FmInfo%phis
phim => FmInfo%phim
phic => CmInfo%phic
Jout => Cminfo%RadJout
AxSrc => FmInfo%AxSrc
AxPXS => FmInfo%AxPXS
PinXS => mklCMFD%PinXS

outTol = mklCntl%outerConv
outMin = mklCntl%minOuter
outMax = mklCntl%maxOuter
lGcCmfd = mklCntl%lGcCMFD
lAxRefFDM = mklCntl%lAxRefFDM

myzb = mklGeom%myzb; myze = mklGeom%myze
ng = mklGeom%ng; nxy = mklGeom%nxy; nInIter = 1
lxslib = nTracerCntl%lxslib; lscat1 = nTracerCntl%lScat1; l3dim = nTracerCntl%l3dim
loutConv = FALSE; outIter = 0

lDhat = TRUE
IF (ItrCntl%CMFDIt .EQ. ItrCntl%CMFDIt0) lDhat = FALSE

IF (lFirstCMFD) lChebyshev = mklCntl%lChebyshev
mklCntl%lChebyshev = lChebyshev .AND. .NOT. lFirstCMFD
mklCntl%lPardiso = mklCntl%lDirect .AND. (mklCntl%lDcpl .OR. PE%nCMFDproc .EQ. 1)

CmfdInitBeg = nTracer_dclock(.FALSE., .FALSE.)

WRITE(mesg,'(a)') 'Cell Homogenization...'
IF (PE%Master) CALL message(io8, TRUE, TRUE, mesg)
CALL FxrChiGen(CoreInfo, Fxr, FmInfo, GroupInfo, PE, nTracerCntl)
CALL HomogenizeXS(CoreInfo, mklGeom%superPin, Fxr, PinXS, phis, ng, nxy, myzb, myze, lxslib, lscat1, FALSE)
CALL SetRadialCoupling(mklGeom%superPin, PinXS, Jout, ng, nxy, myzb, myze, lDhat)
IF (lAxRefFDM) CALL SetReflectorDhatZero(PinXS)
IF (l3dim) CALL SetAxialDtil(PinXS)
IF (l3dim .AND. lxslib) CALL HomogenizePnXS(CoreInfo, FmInfo, GroupInfo)

CALL SetCMFDPhis(mklCMFD, PinXS, lFirstCMFD)

CmfdInitEnd = nTracer_dclock(.FALSE., .FALSE.)
TimeChk%CmfdInitTime = TimeChk%CmfdInitTime + (CmfdInitEnd - CmfdInitBeg)

IF (mklCntl%lAxial) THEN
  IF (l3dim .AND. .NOT. lFirstCMFD) THEN
    CALL GetNeighborFlux(mklCMFD)
    CALL MKL_AxialSolver(CoreInfo, CmInfo, PinXS, eigv)
  ENDIF
ENDIF

CmfdInitBeg = nTracer_dclock(.FALSE., .FALSE.)

WRITE(mesg, '(a)') 'Linear System Construction...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)

IF (mklCntl%lChebyshev) CALL ChebyshevReset()

CALL SetSourceOperator(mklCMFD, PinXS)

IF (mklCntl%lJacobi) THEN
  CALL SetCsrJacobiSystem(mklCMFD, PInXS, l3dim, 0.0)
ELSE
  CALL SetCsrBiCGSystem(mklCMFD, PinXS, l3dim, lFirstCMFD, 0.0)
ENDIF

CmfdInitEnd = nTracer_dclock(.FALSE., .FALSE.)
TimeChk%CmfdInitTime = TimeChk%CmfdInitTime + (CmfdInitEnd - CmfdInitBeg)

WRITE(mesg, '(a)') 'Performing CMFD Calculation...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)

CALL CMFDPsiUpdt(mklCMFD)

DO WHILE (.NOT. loutConv)
  DO iter = 1, nGroupInfo
    GrpBeg = 1; GrpEnd = ng
    IF (iter .GT. 1) THEN
      GrpBeg = GroupInfo%UpScatRange(1); GrpEnd = GroupInfo%UpScatRange(2)
    ENDIF
    DO InIter = 1, nInIter
      IF (mklCntl%lDcpl) CALL GetNeighborFlux(mklCMFD)
      IF (mklCntl%lPardiso) THEN
        DO ig = GrpBeg, GrpEnd
          CALL CMFDSrcUpdt(mklCMFD, ig, eigv)
          CALL pardisoSolve(mklCMFD%M(ig), mklCMFD%phis(:, :, ig), mklCMFD%src(:, ig))
        ENDDO
      ELSE
        DO ig = GrpBeg, GrpEnd
          CALL CMFDSrcUpdt(mklCMFD, ig, eigv)
          IF (mklCntl%lJacobi) THEN
            CALL AAJ(mklCMFD, mklCMFD%Jacobi(ig), mklCMFD%phis(:, :, ig), mklCMFD%src(:, ig))
          ELSE
            CALL BiCGSTAB(mklCMFD, ig)
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDDO
  CALL CMFDPsiUpdt(mklCMFD)
  CALL CMFDEigUpdt(mklCMFD, eigv)
  IF (mklCntl%lDcpl) CALL GetNeighborFlux(mklCMFD)
!  outRes = CMFDPsiErr(mklCMFD)
  outRes = CMFDResidual(mklCMFD, eigv)
  IF (outIter .EQ. 0) outRes0 = outRes
  outErr = outRes / outRes0
  outIter = outIter + 1
  loutConv = (outErr .LE. outTol) .AND. (outIter .GE. outMin)
  loutConv = loutConv .OR. (outIter .GE. outMax)
  ItrCntl%CMFDIt = ItrCntl%CMFDIt + 1
  WRITE(mesg, '(a9, i9, f22.6, 3x, f10.5, 1p, e15.3)') 'MGOUTER', ItrCntl%CMFDIt, eigv, outErr, outRes
  IF (PE%MASTER) CALL message(io8, FALSE, TRUE, mesg)
  IF (lChebyshev) THEN
    CALL ChebyshevAcc(mklCMFD, mklCntl%lChebyshev)
    IF (mklCntl%lChebyshev) CALL CMFDPsiUpdt(mklCMFD)
  ENDIF
  IF (mod(outIter, 5) .NE. 0) CYCLE
  IF (mklCntl%lAxial .AND. .NOT. lFirstCMFD) THEN
    IF (l3dim) THEN
      CALL GetNeighborFlux(mklCMFD)
      CALL MKL_AxialSolver(CoreInfo, CmInfo, PinXS, eigv)
      CmfdInitBeg = nTracer_dclock(.FALSE., .FALSE.)
      WRITE(mesg, '(a)') 'Linear System Construction...'
      IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)
      CALL SetCsrBiCGSystem(mklCMFD, PinXS, l3dim, FALSE, 0.0)
      CmfdInitEnd = nTracer_dclock(.FALSE., .FALSE.)
      TimeChk%CmfdInitTime = TimeChk%CmfdInitTime + (CmfdInitEnd - CmfdInitBeg)
    ENDIF
  ENDIF
  IF (lGcCMFD) THEN
    IF (mklCntl%lChebyshev) CALL ChebyshevReset()
    CALL GetNeighborFlux(mklCMFD)
    CALL MKL_GcCmfdPower(CoreInfo, CmInfo, eigv)
!    CALL CMFDPsiUpdt(mklCMFD)
  ENDIF
ENDDO

IF (mklCntl%lPardiso) THEN
  DO ig = 1, ng
    CALL pardisoDelete(mklCMFD%M(ig))
  ENDDO
ENDIF

CALL SetMOCPhis(CoreInfo, PinXS, phis, phic)

IF (l3dim) THEN
  CALL GetNeighborFlux(mklCMFD)
  CALL SetAxialSrc(AxSrc, AxPXS, phic)
ENDIF

IF (lFirstCMFD) lFirstCMFD = FALSE

CmfdTImeEnd = nTracer_dclock(FALSE, FALSE)
TimeChk%CmfdTime = TimeChk%CmfdTime + (CmfdTimeEnd - CmfdTimeBeg)

END SUBROUTINE

! ======================================================================= !
!        CMFD Acceleration Employing Generalized Davidson's Method        !
!  http://www.netlib.org/utk/people/JackDongarra/etemplates/node138.html  !
! ======================================================================= !

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
REAL, POINTER :: phis(:, :, :), phic(:, :, :)
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
CALL FxrChiGen(CoreInfo, Fxr, FmInfo, GroupInfo, PE, nTracerCntl)
CALL HomogenizeXS(CoreInfo, mklGeom%superPin, Fxr, PinXS, phis, ng, nxy, myzb, myze, lxslib, lscat1, FALSE)
#ifdef Buckling
IF(nTracerCntl%lBsq) THEN
  CALL HomBuckling(CoreInfo, Fxr, phis, PinXS, myzb, myze, ng, nTracerCntl%bsq, lxsLib)
ENDIF
#endif
CALL SetRadialCoupling(mklGeom%superPin, PinXS, Jout, ng, nxy, myzb, myze, lDhat)
IF (l3dim) CALL SetAxialDtil(PinXS)

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
  ENDIF
ENDIF

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
  ENDIF
  CALL SolveDavidsonEq(Davidson, eigv)
ENDDO

CALL ReorderFlux(2)
CALL SetMOCPhis(CoreInfo, PinXS, phis, phic)

IF (l3dim) THEN
  CALL GetNeighborFlux(mklCMFD)
  CALL SetAxialSrc(AxSrc, AxPXS, phic)
ENDIF

IF (lFirstCMFD) lFirstCMFD = FALSE

CmfdTImeEnd = nTracer_dclock(FALSE, FALSE)
TimeChk%CmfdTime = TimeChk%CmfdTime + (CmfdTimeEnd - CmfdTimeBeg)

END SUBROUTINE

!--- Private Routines -----------------------------------------------------------------------------

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
  mklCntl%lDcpl = .FALSE.
  seigv = 1.0 / (Keff + mklCntl%Shift)
  eigv = 1.0 / (1.0 / Keff - seigv)
  WRITE(mesg, '(a27, f4.2)') 'Wielandt Shift Parameter : ', mklCntl%Shift
  IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)
ELSE
  seigv = 0.0
  eigv = Keff
ENDIF

CALL HomogenizeGcXS(CoreInfo, PinXS, GcPinXS)
CALL SetCMFDPhis(mklGcCMFD, GcPinXS, TRUE)

CALL SetRadialGcCoupling(PinXS, GcPinXS)
IF (l3dim) THEN
  CALL GetNeighborFlux(mklGcCMFD)
  CALL SetAxialGcDtil(GcPinXS)
  CALL SetAxialGcDhat()
ENDIF

CALL SetSourceOperator(mklGcCMFD, GcPinXS)

IF (mklCntl%lJacobi) THEN
  CALL SetCsrJacobiSystem(mklGcCMFD, GcPInXS, l3dim, seigv)
ELSE
  CALL SetCsrBiCGSystem(mklGcCMFD, GcPinXS, l3dim, TRUE, seigv)
ENDIF

CALL CMFDPsiUpdt(mklGcCMFD)

DO WHILE (.NOT. loutConv)
  DO iter = 1, nGroupInfo
    GrpBeg = 1; GrpEnd = ngc
    IF (iter .GT. 1) THEN
      GrpBeg = GcGroupInfo%UpScatRange(1); GrpEnd = GcGroupInfo%UpScatRange(2)
    ENDIF
    IF (mklCntl%lDcpl) CALL GetNeighborFlux(mklGcCMFD)
    IF (mklCntl%lPardiso) THEN
      DO ig = 1, ngc
        CALL CMFDSrcUpdt(mklGcCMFD, ig, eigv)
        CALL pardisoSolve(mklGcCMFD%M(ig), mklGcCMFD%phis(:, :, ig), mklGcCMFD%src(:, ig))
      ENDDO
    ELSE
      DO ig = 1, ngc
        CALL CMFDSrcUpdt(mklGcCMFD, ig, eigv)
        IF (mklCntl%lJacobi) THEN
          CALL AAJ(mklGcCMFD, mklGcCMFD%Jacobi(ig), mklGcCMFD%phis(:, :, ig), mklGcCMFD%src(:, ig))
        ELSE
          CALL BiCGSTAB(mklGcCMFD, ig)
        ENDIF
      ENDDO
    ENDIF
  ENDDO
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
ENDDO

WRITE(mesg, '(a)') 'Group Reconstruction...'
IF (PE%Master) CALL message(io8, TRUE, TRUE, mesg)
CALL GcReconstruction(GcPinXS)

mklCntl%lDcpl = lDcpl

END SUBROUTINE

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
  ENDDO
ENDDO

END SUBROUTINE

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
        ENDIF
        CMFD%S(ir, igf, igt) = val
      ENDDO
    ENDDO
  ENDDO
ENDDO
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
    ENDDO
  ENDDO
ENDDO
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
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

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
      ENDDO
    ENDDO
  ENDDO
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
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDIF

END SUBROUTINE

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
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

CALL dcopy(nxy * nzCMFD, mklGcCMFD%psi, 1, mklCMFD%psi, 1)

END SUBROUTINE

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
      ENDDO
      fmult = phisum / PinXS(ixy_map, iz)%Phi(ig)
      DO j = 1, superPin(ixy_map)%nxy
        ipin = superPin(ixy_map)%pin(j)
        FsrIdxSt = Pin(ipin)%FsrIdxSt
        icel = Pin(ipin)%Cell(iz)
        nLocalFsr = Cell(icel)%nFsr
        DO i = 1, nLocalFsr
          ifsr = FsrIdxSt + i - 1
          phis(ifsr, iz, ig) = phis(ifsr, iz, ig) * fmult
        ENDDO
        phic(ipin, iz, ig) = phisum
      ENDDO
      mklCMFD%phic(ixy, iz, ig) = phisum
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

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
      ENDIF
      Dtil = mklCMFD%AxDtil(bottom, ixy, izf, ig)
      Dhat = mklCMFD%AxDhat(bottom, ixy, izf, ig)
      Jnet = - Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
      DO i = 1, Pin(ixy_map)%nxy
        ipin = Pin(ixy_map)%pin(i)
        AxSrc(ipin, iz, ig) = Jnet
      ENDDO
      !--- Axial Source from Top
      izf = fmRange(iz, top)
      myphi = mklCMFD%phis(ixy, izf, ig)
      IF (iz .EQ. myze) THEN
        neighphi = mklCMFD%neighphis(ixy, ig, top)
      ELSE
        neighphi = mklCMFD%phis(ixy, izf + 1, ig)
      ENDIF
      Dtil = mklCMFD%AxDtil(top, ixy, izf, ig)
      Dhat = mklCMFD%AxDhat(top, ixy, izf, ig)
      Jnet = - Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
      DO i = 1, Pin(ixy_map)%nxy
        ipin = Pin(ixy_map)%pin(i)
        AxSrc(ipin, iz, ig) = AxSrc(ipin, iz, ig) + Jnet
      ENDDO
    ENDDO
  ENDDO
ENDDO
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
      ENDDO
    ENDDO
  ENDDO
ENDDO
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
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDIF

END SUBROUTINE

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
ENDDO

CALL InitFastComm()
CALL GetNeighborFast(ng * nxy, myphis(:, :, bottom), neighphis(:, :, top), bottom)
CALL GetNeighborFast(ng * nxy, myphis(:, :, top), neighphis(:, :, bottom), top)
CALL FinalizeFastComm()

IF (mklGeom%lBottom) THEN
  IF (mklGeom%AxBC(bottom) .EQ. VoidCell) neighphis(:, :, bottom) = 0.0
  IF (mklGeom%AxBC(bottom) .EQ. RefCell) CALL dcopy(ng * nxy, myphis(:, :, bottom), 1, neighphis(:, :, bottom), 1)
ENDIF

IF (mklGeom%lTop) THEN
  IF (mklGeom%AxBC(top) .EQ. VoidCell) neighphis(:, :, top) = 0.0
  IF (mklGeom%AxBC(top) .EQ. RefCell) CALL dcopy(ng * nxy, myphis(:, :, top), 1, neighphis(:, :, top), 1)
ENDIF

DEALLOCATE(myphis)

END SUBROUTINE

END MODULE

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
    ENDDO
  ENDDO
ENDDO

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
    ENDDO
    PinXS(iFuelPin, iz)%XSnf = PinXS(iFuelPin, iz)%XSnf / PinVol(iFuelPin, iz)
    PinXS(iFuelPin, iz)%XSkf = PinXS(iFuelPin, iz)%XSkf / PinVol(iFuelPin, iz)
  ENDDO
ENDDO

END SUBROUTINE

#endif