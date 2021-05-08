#include <defines.h>
!--- CNJ Edit : Driver Routine for Gamma CMFD Calculation
#ifdef __INTEL_MKL
#ifdef __GAMMA_TRANSPORT
MODULE GAMMA_HOMOXS

USE MKL_3D
USE GammaTYPEDEF
USE GamXSUtil
IMPLICIT NONE

TYPE(GamMacXs_Type), POINTER, PRIVATE :: XsMac(:, :)

CONTAINS

SUBROUTINE AllocHomoXSVar(CoreInfo, ngg)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,       Cell_Type
USE CMFD_MOD,       ONLY : CMFDPinXS
USE PE_MOD,         ONLY : PE

IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
INTEGER :: ngg

TYPE(Cell_Type), POINTER :: Cell(:)
INTEGER :: ifxr, icel, ipin, tid
INTEGER :: nCellType, nFxrMax

Cell => CoreInfo%CellInfo
nCellType = CoreInfo%nCellType

nFxrMax = 0
DO icel = 1, nCellType
  nFxrMax = max(nFxrMax, Cell(icel)%nFxr)
ENDDO
nFxrMax = nFxrMax + 3

CALL omp_set_num_threads(PE%nCMFDThread)

ALLOCATE(XsMac(nFxrMax, PE%nCMFDThread))

!$OMP PARALLEL PRIVATE(tid)
tid = omp_get_thread_num() + 1
DO ifxr = 1, nFxrMax
  XsMac(ifxr, tid)%ngg = ngg
  CALL AllocGamMacXs(XsMac(ifxr, tid))
ENDDO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE AllocPinXS(PinXS, GroupInfo, nxy, myzb, myze)
USE PARAM
USE TYPEDEF,        ONLY : GroupInfo_Type
IMPLICIT NONE

TYPE(GPinXS_Type), POINTER :: PinXS(:, :)
TYPE(GroupInfo_Type) :: GroupInfo
INTEGER :: ng, ngg, nxy, myzb, myze

TYPE(GPinXS_Type), POINTER :: myPinXS
INTEGER :: ig, igg, igb, ige, ipin, iz

ng = GroupInfo%ng
ngg = GroupInfo%ngg

ALLOCATE(PinXS(nxy, myzb : myze))

DO iz = myzb, myze
  DO ipin = 1, nxy
    myPinXS => PinXS(ipin, iz)
    ALLOCATE(myPinXS%Dtil(4, ngg)); myPinXS%Dtil = 0.0
    ALLOCATE(myPinXS%Dhat(4, ngg)); myPinXS%Dhat = 0.0
    ALLOCATE(myPinXS%XSD(ngg)); myPinXS%XSD = 0.0
    ALLOCATE(myPinXS%XSt(ngg)); myPinXS%XSt = 0.0
    ALLOCATE(myPinXS%XStr(ngg)); myPinXS%XStr = 0.0
    ALLOCATE(myPinXS%XSr(ngg)); myPinXS%XSr = 0.0
    ALLOCATE(myPinXS%XSa(ngg)); myPinXS%XSa = 0.0
    ALLOCATE(myPinXS%NPhi(ng)); myPinXS%NPhi = 0.0
    ALLOCATE(myPinXS%GPhi(ngg)); myPinXS%GPhi = 0.0
    ALLOCATE(myPinXS%XSp(ng, ngg)); myPinXS%XSp = 0.0
    ALLOCATE(myPinXS%XSs(ngg))
    DO igg = 1, ngg
      igb = GroupInfo%InScatRange_ph(1, igg); ige = GroupInfo%InScatRange_ph(2, igg)
      ALLOCATE(myPinXS%XSs(igg)%from(igb : ige)); myPinXS%XSs(igg)%from = 0.0
      myPinXS%XSs(igg)%ib = igb; myPinXS%XSs(igg)%ie = ige
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE HomogenizeXS(CoreInfo, superPin, Fxr, PinXS, phis, gphis, ng, ngg, nxy, myzb, myze, lscat1)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,       FxrInfo_Type,       PinXS_Type,                                     &
                           Pin_Type,            Cell_Type
USE CORE_MOD,       ONLY : GroupInfo
USE CNTL,           ONLY : nTracerCntl
USE GamXsLib_Mod,   ONLY : GamXsBase,           GamScatMatrix,      GamProdMatrix
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(superPin_Type), POINTER :: superPin(:)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(GPinXS_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: phis(:, :, :), gphis(:, :, :)
INTEGER :: ng, ngg, nxy, myzb, myze
LOGICAL :: lscat1

TYPE(FxrInfo_Type), POINTER :: myFxr
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)

INTEGER :: i, j, ixy, icel, ipin, iz, ifxr, ifxr_global, itype, tid
INTEGER :: FsrIdxSt, FxrIdxSt
INTEGER :: norg, nLocalFxr, nFsrInFxr
INTEGER :: irgb, irge

Pin => CoreInfo%Pin
Cell => CoreInfo%CellInfo


!$OMP PARALLEL PRIVATE(tid, ipin, icel, ifxr, ifxr_global, itype, FsrIdxSt, FxrIdxSt, nLocalFxr, myFxr)
tid = omp_get_thread_num() + 1

DO iz = myzb, myze
  !$OMP DO SCHEDULE(GUIDED)
  DO ixy = 1, nxy
    nLocalFxr = 0
    DO j = 1, superPin(ixy)%nxy
      ipin = superPin(ixy)%pin(j)
      icel = Pin(ipin)%Cell(iz)
      FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
      DO i = 1, Cell(icel)%nFxr
        ifxr = nLocalFxr + i; ifxr_global = FxrIdxSt + i - 1; nFsrInFxr = Cell(icel)%nFsrInFxr(i)
        myFxr => Fxr(ifxr_global, iz)
!         print*, ifxr, j, 'gamxs', myfxr%niso
        CALL GamXsBase(XsMac(ifxr, tid), myFxr, 1, ngg, ngg, FALSE, TRUE)
!         print*, i, j, 'gamprod', XsMac(ifxr, tid)%lProdAlloc
        CALL GamProdMatrix(XsMac(ifxr, tid), myFxr, 1, ngg, ng, ngg, GroupInfo, FALSE)
!  print*, Xsmac(20,1)%lProdAlloc
        ! print*, i, j, 'gamsm'
         CALL GamScatMatrix(XsMac(ifxr, tid), myFxr, 1, ngg, ngg, FALSE, TRUE)
!  print*, Xsmac(20,1)%lProdAlloc
        XsMac(ifxr, tid)%XsMacTr = XsMac(ifxr, tid)%XsMacA + XsMac(ifxr, tid)%XsMacStr
        XsMac(ifxr, tid)%XsMacT = XsMac(ifxr, tid)%XsMacA + XsMac(ifxr, tid)%XsMacS
      ENDDO
      nLocalFxr = nLocalFxr + Cell(icel)%nFxr
    ENDDO
!    print*, 'homcellxs', ixy, nxy
    CALL HomogenizeCellXS(CoreInfo, superPin(ixy), PinXS(ixy, iz), XsMac(1 : nLocalFxr, tid), phis, gphis, iz, ng, ngg, lscat1)
  ENDDO
  !$OMP END DO
ENDDO

!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE HomogenizeCellXS(CoreInfo, superPin, PinXS, XsMac, phis, gphis, iz, ng, ngg, lscat1)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,       Pin_Type,           Cell_Type,                                      &
                           PinXS_Type,          XsMac_Type
USE CORE_MOD,       ONLY : GroupInfo
USE CNTL,           ONLY : nTracerCntl
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(superPin_Type) :: superPin
TYPE(GPinXS_Type) :: PinXS
TYPE(GamMacXs_Type) :: XsMac(:)
REAL, POINTER :: phis(:, :, :), gphis(:, :, :)
INTEGER :: iz, ng, ngg
LOGICAL :: lscat1

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
INTEGER :: i, j, ixy, ifxr, ipin, icel, ig, igg, ig2, igb, ige, ireg, ifsr
INTEGER :: nxy, nFxr, nFsrInFxr, nLocalFxr, FsrIdxSt
REAL :: localphi, localgphi, vol, phisum, gphisum, volsum
REAL :: RR(3), RRS(ngg), RRP(ngg)
LOGICAL :: lFuel

Pin => CoreInfo%Pin
Cell => CoreInfo%CellInfo
nxy = superPin%nxy

DO igg = 1, ngg
  PinXS%XSs(igg)%from = 0.0
ENDDO

DO igg = 1, ngg
  igb = GroupInfo%OutScatRange_Ph(1, igg)
  ige = GroupInfo%OutScatRange_Ph(2, igg)
  RR = 0.0; RRS = 0.0; gphisum = 0.0; volsum = 0.0
  nLocalFxr = 0
  DO ixy = 1, nxy
    ipin = superPin%pin(ixy)
    icel = Pin(ipin)%Cell(iz)
    nFxr = Cell(icel)%nFxr
    FsrIdxSt = Pin(ipin)%FsrIdxSt
    DO i = 1, nFxr
      ifxr = nLocalFxr + i
      nFsrInFxr = Cell(icel)%nFsrInFxr(i)
      localgphi = 0.0
      DO j = 1, nFsrInFxr
        ireg = Cell(icel)%MapFxr2FsrIdx(j, i)
        ifsr = FsrIdxSt + ireg - 1
        vol = Cell(icel)%vol(ireg)
        IF (lscat1) THEN
          localgphi = localgphi + gphis(ifsr, igg, iz) * vol
        ELSE
          localgphi = localgphi + gphis(igg, ifsr, iz) * vol
        ENDIF
        volsum = volsum + vol
      ENDDO
      gphisum = gphisum + localgphi
      RR(1) = RR(1) + localgphi * XsMac(ifxr)%XsMacT(igg)
      RR(2) = RR(2) + localgphi * XsMac(ifxr)%XsMacTr(igg)
      RR(3) = RR(3) + localgphi * XsMac(ifxr)%XsMacA(igg)
      DO ig2 = igb, ige
        RRS(ig2) = RRS(ig2) + localgphi * XsMac(ifxr)%XsMacSm(igg, ig2)
      ENDDO
    ENDDO
    nLocalFxr = nLocalFxr + nFxr
  ENDDO
  PinXS%GPhi(igg) = gphisum / volsum
  RR = RR / gphisum
  RRS = RRS / gphisum
  PinXS%XSt(igg) = RR(1)
  PinXS%XStr(igg) = RR(2)
  PinXS%XSa(igg) = RR(3)
  DO ig2 = igb, ige
    IF ((igg - PinXS%XSs(ig2)%ib) * (igg - PinXS%XSs(ig2)%ie) .GT. 0) CYCLE
    PinXS%XSs(ig2)%from(igg) = RRS(ig2)
  ENDDO
  PinXS%XSr(igg) = RR(3) + sum(RRS) - RRS(igg)
  PinXS%XSs(igg)%self = PinXS%XSs(igg)%from(igg)
  PinXS%XSs(igg)%from(igg) = 0.0
  PinXs%XSD(igg) = 1.0 / 3.0 / PinXS%XStr(igg)
ENDDO

PinXS%XSP = 0.0

DO ig = 1, ng
  volsum = 0.0; phisum = 0.0; RRP = 0.0
  nLocalFxr = 0
  DO ixy = 1, nxy
    ipin = superPin%pin(ixy)
    icel = Pin(ipin)%Cell(iz)
    nFxr = Cell(icel)%nFxr
    FsrIdxSt = Pin(ipin)%FsrIdxSt
    DO i = 1, nFxr
      ifxr = nLocalFxr + i
      nFsrInFxr = Cell(icel)%nFsrInFxr(i)
      localphi = 0.0
      DO j = 1, nFsrInFxr
        ireg = Cell(icel)%MapFxr2FsrIdx(j, i)
        ifsr = FsrIdxSt + ireg - 1
        vol = Cell(icel)%vol(ireg)
        localphi = phis(ifsr, iz, ig) * vol
        phisum = phisum + localphi
        volsum = volsum + vol
        DO igg = 1, ngg
          RRP(igg) = RRP(igg) + localphi * XsMac(ifxr)%GProdTot(ig, igg)
        ENDDO
      ENDDO
    ENDDO
    nLocalFxr = nLocalFxr + nFxr
  ENDDO
  PinXS%NPhi(ig) = phisum / volsum
  RRP = RRP / phisum
  PinXS%XSP(ig, :) = RRP
ENDDO

END SUBROUTINE

SUBROUTINE HomogenizePnXS(CoreInfo, FmInfo)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,       FmInfo_Type,        FxrInfo_Type,       Pin_Type
USE GamXsLib_Mod,   ONLY : GamP1XsScatMatrix,   GamP2XsScatMatrix,  GamP3XsScatMatrix
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(FmInfo_Type) :: FmInfo

TYPE(superPin_Type), POINTER :: superPin(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(GamMacXs_Type) :: XsMac
INTEGER :: ig, ifxr, ixy, ixy_map, ipin, iz, izf
INTEGER :: ngg, nxy, nzCMFD
INTEGER, POINTER :: pinMap(:), planeMap(:)

Pin => CoreInfo%Pin
Fxr => FmInfo%Fxr
superPin => mklGeom%superPin
ngg = mklGeom%ngg
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
    CALL GamP1XsScatMatrix(XsMac, Fxr(ifxr, iz), 1, ngg, ngg)
!    CALL GamP2XsScatMatrix(XsMac, Fxr(ifxr, iz), 1, ngg, ngg)
!    CALL GamP3XsScatMatrix(XsMac, Fxr(ifxr, iz), 1, ngg, ngg)
    mklAxial%SmP1(:, :, izf, ixy) = XsMac%MacGSM1
!    mklAxial%SmP2(:, :, izf, ixy) = XsMac%MacGSM2
!    mklAxial%SmP3(:, :, izf, ixy) = XsMac%MacGSM3
  ENDDO
ENDDO

END SUBROUTINE

END MODULE

MODULE GAMMA_CMFD

USE MKL_3D
USE GAMMA_HOMOXS
USE GAMMA_AXIAL
IMPLICIT NONE

LOGICAL :: lFirstCMFD = TRUE

PRIVATE
PUBLIC :: GammaCMFDAcc

CONTAINS

SUBROUTINE GammaCMFDAcc(CoreInfo, CmInfo, FmInfo)
USE PARAM
USE TYPEDEF,            ONLY : CoreInfo_Type,       CmInfo_Type,        FmInfo_Type,        FxrInfo_Type
USE GammaCore_MOD,      ONLY : GPhiAngIn,           GPhis,              GPhic,              GJout,                  &
                               GAxSrc,              GAxPXS
USE Core_mod,           ONLY : GroupInfo
USE PE_MOD,             ONLY : PE
USE CNTL,               ONLY : nTracerCntl
USE TIMER,              ONLY : nTracer_dclock,      TimeChk
USE MKL_POWER,          ONLY : CMFDResidual
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(FmInfo_Type) :: FmInfo

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(GPinXS_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: phis(:, :, :)
REAL :: CmfdTimeBeg, CmfdTimeEnd, CmfdInitBeg, CmfdInitEnd
INTEGER :: ig, iter
INTEGER :: myzb, myze
INTEGER :: ng, ngg, nxy
INTEGER :: GrpBeg, GrpEnd, nGroupInfo = 2
LOGICAL :: lDhat, lScat1, l3dim

CmfdTImeBeg = nTracer_dclock(FALSE, FALSE)

Fxr => FmInfo%Fxr
Phis => FmInfo%Phis
PinXS => mklGammaCMFD%GPinXS
myzb = mklGeom%myzb; myze = mklGeom%myze
ng = mklGeom%ng; ngg = mklGeom%ngg; nxy = mklGeom%nxy
lScat1 = nTracerCntl%lGammaScat1; l3dim = nTracerCntl%l3dim
lDhat = .NOT. lFirstCMFD
IF (.NOT. GroupInfo%lUpScat_ph) nGroupInfo = 1

CmfdInitBeg = nTracer_dclock(.FALSE., .FALSE.)

CALL HomogenizeXS(CoreInfo, mklGeom%superPin, Fxr, PinXS, phis, gphis, ng, ngg, nxy, myzb, myze, lScat1)
CALL SetRadialCoupling(mklGeom%superPin, PinXS, gJout, ngg, nxy, myzb, myze, lDhat, lScat1)
IF (l3dim) THEN
  CALL SetAxialDtil(mklGammaCMFD, mklGammaAxial)
  CALL HomogenizePnXS(CoreInfo, FmInfo)
ENDIF

CALL SetCMFDPhis(mklGammaCMFD, PinXS, lFirstCMFD)

CmfdInitEnd = nTracer_dclock(.FALSE., .FALSE.)
TimeChk%CmfdInitTime = TimeChk%CmfdInitTime + (CmfdInitEnd - CmfdInitBeg)

IF (mklCntl%lAxial) THEN
  IF (l3dim .AND. .NOT. lFirstCMFD) THEN
    CALL GetNeighborFlux(mklGammaCMFD)
    CALL GammaAxialSolver()
  ENDIF
ENDIF

CALL SetCsrBiCGSystem(mklGammaCMFD, l3dim, lFirstCMFD)

DO iter = 1, nGroupInfo
  GrpBeg = 1; GrpEnd = ngg
  IF (iter .GT. 1) THEN
    GrpBeg = GroupInfo%UpScatRange_Ph(1); GrpEnd = GroupInfo%UpScatRange_Ph(2)
  ENDIF
  DO ig = GrpBeg, GrpEnd
!    print*, 'srcupdt'
    CALL CMFDSrcUpdt(mklGammaCMFD, ig)
!    print*, 'srcupdt fin', mklGeom%l3dim
    CALL BiCGSTAB(mklGammaCMFD, ig)
!    print*, 'bicgstagb fin', mklGeom%l3dim
  ENDDO
ENDDO

CALL SetMOCPhis(CoreInfo, PinXS, gphis, gphic, lScat1)

IF (l3dim) THEN
  CALL GetNeighborFlux(mklGammaCMFD)
  CALL SetAxialSrc(GAxSrc, GAxPXS, gphic)
ENDIF


IF (lFirstCMFD) lFirstCMFD = FALSE

CmfdTImeEnd = nTracer_dclock(FALSE, FALSE)
TimeChk%CmfdTime = TimeChk%CmfdTime + (CmfdTimeEnd - CmfdTimeBeg)

END SUBROUTINE

SUBROUTINE SetCMFDPhis(CMFD, PinXS, lCopy)

IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
TYPE(GPinXS_Type), POINTER :: PinXS(:, :)
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
        CMFD%phis(ipin, izf, ig) = PinXS(ipin_map, iz)%GPhi(ig)
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
        fmult = PinXS(ipin_map, iz)%GPhi(ig) / CMFD%phic(ipin, iz, ig)
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

SUBROUTINE SetRadialCoupling(Pin, PinXS, Jout, ng, nxy, myzb, myze, lDhat, lScat1)
USE PARAM
USE cntl,     ONLY : nTracerCntl
IMPLICIT NONE

TYPE(superPin_Type), POINTER :: Pin(:)
TYPE(GPinXS_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: Jout(:, :, :, :, :)
INTEGER :: ng, nxy, myzb, myze
LOGICAL :: lDhat, lScat1

INTEGER :: ig, ipin, ineighpin, iz, ibd, inbd
REAL :: Dtil, Dhat, myphi, neighphi, mybeta, neighbeta, jnet, jfdm, smy
REAL, POINTER :: superJout(:, :, :, :, :)

IF(nTracerCntl%lHex) THEN
  CALL HexSetRadialCoupling(Pin, PinXS, Jout, ng, nxy, myzb, myze, lDhat, lScat1)
  RETURN
END IF


ALLOCATE(superJout(3, 4, nxy, myzb : myze, ng))

CALL superPinCurrent(Pin, Jout, superJout, ng, nxy, myzb, myze, lScat1)

!$OMP PARALLEL PRIVATE(ineighpin, inbd, myphi, neighphi, mybeta, neighbeta, Dtil, Dhat, jnet, jfdm, smy)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO iz = myzb, myze
    DO ipin = 1, nxy
      DO ibd = 1, 4
        ineighpin = Pin(ipin)%NeighIdx(ibd)
        inbd = Pin(ipin)%NeighSurfIdx(ibd)
        smy = Pin(ipin)%BdLength(ibd)
        myphi = PinXS(ipin, iz)%GPhi(ig)
        mybeta = PinXS(ipin, iz)%XSD(ig) / Pin(ipin)%Center2SurfaceL(ibd)
        IF (ineighpin .GT. 0) THEN
          neighphi = PinXS(ineighpin, iz)%GPhi(ig)
          neighbeta = PinXS(ineighpin, iz)%XSD(ig) / Pin(ineighpin)%Center2SurfaceL(inbd)
          Dtil = mybeta * neighbeta / (mybeta + neighbeta) * smy
          jfdm = - Dtil * (neighphi - myphi)
          jnet = superJout(2, ibd, ipin, iz, ig) - superJout(1, ibd, ipin, iz, ig)
          Dhat = - (jnet - jfdm) / (myphi + neighphi)
        ELSE
          IF (ineighpin .EQ. VoidCell) THEN
            neighbeta = 0.5; neighphi = 0.0
          ELSEIF (ineighpin .EQ. RefCell) THEN
            neighbeta = 0.0; neighphi = myphi
          ENDIF
          Dtil = mybeta * neighbeta / (mybeta + neighbeta) * smy
          jfdm = - Dtil * (neighphi - myphi)
          jnet = superJout(2, ibd, ipin, iz, ig) - superJout(1, ibd, ipin, iz, ig)
          Dhat = - (jnet - jfdm) / (myphi + neighphi)
        ENDIF
        PinXS(ipin, iz)%Dtil(ibd, ig) = Dtil
        IF (lDhat) PinXS(ipin, iz)%Dhat(ibd, ig) = Dhat
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

DEALLOCATE(superJout)

END SUBROUTINE

SUBROUTINE HexSetRadialCoupling(Pin, PinXS, Jout, ng, nxy, myzb, myze, lDhat, lScat1)
USE PARAM
USE geom,    ONLY : ncbd
IMPLICIT NONE

TYPE(superPin_Type), POINTER :: Pin(:)
TYPE(GPinXS_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: Jout(:, :, :, :, :)
INTEGER :: ng, nxy, myzb, myze
LOGICAL :: lDhat, lScat1

INTEGER :: ig, ipin, ineighpin, iz, iNgh, ibd, jNgh, jbd
REAL :: Dtil, Dhat, myphi, neighphi, mybeta, neighbeta, jnet, jfdm, smy
REAL, POINTER :: superJout(:, :, :, :, :)

ALLOCATE(superJout(3, ncbd, nxy, myzb : myze, ng)) ! # of Ngh is fixed as 15, artibrary #

CALL HexsuperPinCurrent(Pin, Jout, superJout, ng, nxy, myzb, myze, lScat1)

!$OMP PARALLEL PRIVATE(ineighpin, iNgh, jNgh, ibd, jbd, myphi, neighphi, mybeta, neighbeta, Dtil, Dhat, jnet, jfdm, smy)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO iz = myzb, myze
    DO ipin = 1, nxy
      DO iNgh = 1, Pin(ipin)%nNgh
        ibd = Pin(ipin)%NghBd(iNgh)
        ineighpin = Pin(ipin)%NeighIdx(iNgh)
        jNgh = Pin(ipin)%NeighSurfIdx(iNgh)
        smy = Pin(ipin)%NghLgh(iNgh)
        myphi = PinXS(ipin, iz)%GPhi(ig)
        mybeta = PinXS(ipin, iz)%XSD(ig) / Pin(ipin)%Center2SurfaceL(ibd)
        IF (ineighpin .GT. 0) THEN
          jbd = Pin(ineighpin)%NghBd(jNgh)
          neighphi = PinXS(ineighpin, iz)%GPhi(ig)
          neighbeta = PinXS(ineighpin, iz)%XSD(ig) / Pin(ineighpin)%Center2SurfaceL(jbd)
          Dtil = mybeta * neighbeta / (mybeta + neighbeta) * smy
          jfdm = - Dtil * (neighphi - myphi)
          jnet = superJout(2, iNgh, ipin, iz, ig) - superJout(1, iNgh, ipin, iz, ig)
          Dhat = - (jnet - jfdm) / (myphi + neighphi)
        ELSE
          IF (ineighpin .EQ. VoidCell) THEN
            neighbeta = 0.5; neighphi = 0.0
          ELSEIF (ineighpin .EQ. RefCell) THEN
            neighbeta = 0.0; neighphi = myphi
          ENDIF
          Dtil = mybeta * neighbeta / (mybeta + neighbeta) * smy
          jfdm = - Dtil * (neighphi - myphi)
          jnet = superJout(2, iNgh, ipin, iz, ig) - superJout(1, iNgh, ipin, iz, ig)
          Dhat = - (jnet - jfdm) / (myphi + neighphi)
        ENDIF
        PinXS(ipin, iz)%Dtil(iNgh, ig) = Dtil
        IF (lDhat) PinXS(ipin, iz)%Dhat(iNgh, ig) = Dhat
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

DEALLOCATE(superJout)

END SUBROUTINE HexSetRadialCoupling

SUBROUTINE superPinCurrent(Pin, Jout, superJout, ng, nxy, myzb, myze, lScat1)

IMPLICIT NONE

TYPE(superPin_Type), POINTER :: Pin(:)
REAL, POINTER :: Jout(:, :, :, :, :), superJout(:, :, :, :, :)
INTEGER :: ig, ix, iy, iz, ixy, ipin
INTEGER :: ng, nx, ny, nxy
INTEGER :: myzb, myze
LOGICAL :: lScat1

superJout = 0.0

IF (lScat1) THEN
  !$OMP PARALLEL PRIVATE(nx, ny, ipin)
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
  DO ig = 1, ng
    DO iz = myzb, myze
      DO ixy = 1, nxy
        nx = Pin(ixy)%nx
        ny = Pin(ixy)%ny
        DO ix = 1, nx
          ipin = Pin(ixy)%pin2D(ix, 1)
          superJout(:, NORTH, ixy, iz, ig) = superJout(:, NORTH, ixy, iz, ig) + Jout(:, NORTH, ipin, ig, iz)
          ipin = Pin(ixy)%pin2D(ix, ny)
          superJout(:, SOUTH, ixy, iz, ig) = superJout(:, SOUTH, ixy, iz, ig) + Jout(:, SOUTH, ipin, ig, iz)
        ENDDO
        DO iy = 1, ny
          ipin = Pin(ixy)%pin2D(1, iy)
          superJout(:, WEST, ixy, iz, ig) = superJout(:, WEST, ixy, iz, ig) + Jout(:, WEST, ipin, ig, iz)
          ipin = Pin(ixy)%pin2D(nx, iy)
          superJout(:, EAST, ixy, iz, ig) = superJout(:, EAST, ixy, iz, ig) + Jout(:, EAST, ipin, ig, iz)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ELSE
  !$OMP PARALLEL PRIVATE(nx, ny, ipin)
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
  DO ig = 1, ng
    DO iz = myzb, myze
      DO ixy = 1, nxy
        nx = Pin(ixy)%nx
        ny = Pin(ixy)%ny
        DO ix = 1, nx
          ipin = Pin(ixy)%pin2D(ix, 1)
          superJout(:, NORTH, ixy, iz, ig) = superJout(:, NORTH, ixy, iz, ig) + Jout(:, ig, NORTH, ipin, iz)
          ipin = Pin(ixy)%pin2D(ix, ny)
          superJout(:, SOUTH, ixy, iz, ig) = superJout(:, SOUTH, ixy, iz, ig) + Jout(:, ig, SOUTH, ipin, iz)
        ENDDO
        DO iy = 1, ny
          ipin = Pin(ixy)%pin2D(1, iy)
          superJout(:, WEST, ixy, iz, ig) = superJout(:, WEST, ixy, iz, ig) + Jout(:, ig, WEST, ipin, iz)
          ipin = Pin(ixy)%pin2D(nx, iy)
          superJout(:, EAST, ixy, iz, ig) = superJout(:, EAST, ixy, iz, ig) + Jout(:, ig, EAST, ipin, iz)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDIF

END SUBROUTINE

SUBROUTINE HexSuperPinCurrent(Pin, Jout, superJout, ng, nxy, myzb, myze, lScat1)
IMPLICIT NONE

TYPE(superPin_Type), POINTER :: Pin(:)
REAL, POINTER :: Jout(:, :, :, :, :), superJout(:, :, :, :, :)
INTEGER :: ig, ix, iy, iz, ixy, jxy, ipin, iNgh, iBndy, jBndy
INTEGER :: ng, nxy
INTEGER :: myzb, myze
LOGICAL :: lScat1
REAL    :: ratio
! ----------------------------------------------------

superJout = 0.0

!$OMP PARALLEL PRIVATE(ig, iz, ixy, iBndy, jxy, iPin, jBndy)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO iz = myzb, myze
    DO ixy = 1, nxy
      DO iNgh = 1, Pin(ixy)%nNgh
        iBndy = Pin(ixy)%NghBd(iNgh)
        ratio = Pin(ixy)%NghLgh(iNgh) / Pin(ixy)%BdLength(iBndy)

        DO jxy = 1, Pin(ixy)%nBdmPin(iBndy)
          iPin  = Pin(ixy)%BdMPidx(jxy, iBndy) ! MOC Pin
          jBndy = Pin(ixy)%BdMPsuf(jxy, iBndy) ! MOC Suf

          superJout(:, iNgh, ixy, iz, ig) = superJout(:, iNgh, ixy, iz, ig) &
                                               + Jout(:, jBndy, iPin, iz, ig) * ratio
        END DO
      END DO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

IF(lScat1) THEN
  !$OMP PARALLEL PRIVATE(iBndy, iPin, jBndy, ratio)
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
  DO ig = 1, ng
    DO iz = myzb, myze
      DO ixy = 1, nxy
        DO iNgh = 1, Pin(ixy)%nNgh
          iBndy = Pin(ixy)%NghBd(iNgh)
          ratio = Pin(ixy)%NghLgh(iNgh) / Pin(ixy)%BdLength(iBndy)
          DO jxy = 1, Pin(ixy)%nBdmPin(iBndy)
            iPin = Pin(ixy)%BdMPidx(jxy, iBndy)
            jBndy = Pin(ixy)%BdMPsuf(jxy, iBndy)
            superJout(:, iNgh, ixy, iz, ig) = superJout(:, iNgh, ixy, iz, ig) + Jout(:, jBndy, iPin, ig, iz) * ratio
          END DO
        END DO
      END DO
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
ELSE
  !$OMP PARALLEL PRIVATE(iBndy, iPin, jBndy, ratio)
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
  DO ig = 1, ng
    DO iz = myzb, myze
      DO ixy = 1, nxy
        DO iNgh = 1, Pin(ixy)%nNgh
          iBndy = Pin(ixy)%NghBd(iNgh)
          ratio = Pin(ixy)%NghLgh(iNgh) / Pin(ixy)%BdLength(iBndy)
          DO jxy = 1, Pin(ixy)%nBdmPin(iBndy)
            iPin = Pin(ixy)%BdMPidx(jxy, iBndy)
            jBndy = Pin(ixy)%BdMPsuf(jxy, iBndy)
            superJout(:, iNgh, ixy, iz, ig) = superJout(:, iNgh, ixy, iz, ig) + Jout(:, ig, jBndy, iPin, iz) * ratio
          END DO
        END DO
      END DO
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
END IF

! ----------------------------------------------------

END SUBROUTINE HexSuperPinCurrent


SUBROUTINE SetCsrBiCGSystem(CMFD, l3dim, lPrecond)
USE PARAM
USE geom,           ONLY : ncbd
USE MKL_BILU,       ONLY : MKL_PrepareILU
IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
LOGICAL :: l3dim, lPrecond

TYPE(GPinXS_Type), POINTER :: PinXS(:, :)
TYPE(superPin_Type), POINTER :: Pin(:)
INTEGER, POINTER :: pinMap(:), pinMapRev(:), planeMap(:)
INTEGER :: DOWN, UP, SELF
INTEGER :: ir, ic, iz, izf, ig, igf, igt, ibd, isurf, ipin, ipin_map, ineighpin
INTEGER :: dz, gb, ge
INTEGER :: ng, nxy, nzCMFD, nbd
REAL, POINTER :: PinVolFm(:, :), hzfm(:)
REAL :: diagVal(ncbd+3), Dtil, Dhat, val

DOWN = ncbd+1
UP   = ncbd+2
SELF = ncbd+3

PinXS => CMFD%GPinXS
Pin => mklGeom%superPin
ng = CMFD%ng
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
hzfm => mklGeom%hzfm
pinMap => mklGeom%pinMap
pinMapRev => mklGeom%pinMapRev
planeMap => CMFD%planeMap
PinVolFm => mklGeom%PinVolFm

!--- Set Group Major Diffusion Operator

!$OMP PARALLEL PRIVATE(diagVal, iz, ipin_map, Dtil, Dhat, isurf, ineighpin, dz, ir, ic, val)
!$OMP DO SCHEDULE(DYNAMIC)
DO ig = 1, ng
  CALL createCsr(CMFD%M(ig), (ncbd+3) * nxy * nzCMFD, nxy * nzCMFD, nxy * nzCMFD)
  DO izf = 1, nzCMFD
    iz = planeMap(izf)
    DO ipin = 1, nxy
      ir = ipin + (izf - 1) * nxy
      ipin_map = pinMap(ipin)
      diagVal = 0.0
      DO ibd = 1, pin(ipin_map)%nNgh
        Dtil = PinXS(ipin_map, iz)%Dtil(ibd, ig)
        Dhat = PinXS(ipin_map, iz)%Dhat(ibd, ig)
        diagVal(ibd) = - (Dtil + Dhat) * hzfm(izf)
        diagVal(SELF) = diagVal(SELF) + (Dtil - Dhat) * hzfm(izf)
      ENDDO
      IF (l3dim) THEN
        DO ibd = DOWN, UP
          Dtil = CMFD%AxDtil(ibd - ncbd, ipin, izf, ig)
          Dhat = CMFD%AxDhat(ibd - ncbd, ipin, izf, ig)
          diagVal(ibd) = - (Dtil + Dhat) * PinVolFm(ipin, izf) / hzfm(izf)
          diagVal(SELF) = diagVal(SELF) + (Dtil - Dhat) * PinVolFm(ipin, izf) / hzfm(izf)
        ENDDO
      ENDIF
      diagVal(SELF) = diagVal(SELF) + PinVolFm(ipin, izf) * PinXS(ipin_map, iz)%XSr(ig)
      DO ibd = 1, ncbd+3
        isurf = ibd
        IF (isurf .EQ. SELF) THEN
          ineighpin = ipin
          dz = 0
        ELSEIF (isurf .EQ. UP) THEN
          ineighpin = ipin
          dz = 1
        ELSEIF (isurf .EQ. DOWN) THEN
          ineighpin = ipin
          dz = -1
        ELSE
          ineighpin = Pin(ipin_map)%NeighIdx(isurf)
          ineighpin = pinMapRev(ineighpin)
          IF (ineighpin .LE. 0) diagVal(isurf) = 0.0
          IF (ineighpin .EQ. ipin) THEN
            Dtil = PinXS(ipin_map, iz)%Dtil(isurf, ig)
            Dhat = PinXS(ipin_map, iz)%Dhat(isurf, ig)
            diagVal(SELF) = diagVal(SELF) - (Dtil + Dhat) * hzfm(izf)
            diagVal(isurf) = 0.0
          ENDIF
          dz = 0
        END IF
        ic = ir + (ineighpin - ipin) + dz * nxy
        val = diagVal(isurf)
        CALL pushCsr(CMFD%M(ig), val, ir, ic)
      ENDDO
    ENDDO
  ENDDO
  CALL finalizeSortCsr(CMFD%M(ig), FALSE)
ENDDO
!$OMP END DO
!$OMP END PARALLEL

!--- Factorization

IF (lPrecond) THEN
  !$OMP PARALLEL DO SCHEDULE(DYNAMIC)
  DO ig = 1, ng
    CALL MKL_PrepareILU(CMFD%M(ig), CMFD%ILU(ig))
  ENDDO
  !$OMP END PARALLEL DO
ENDIF

!--- Set Axial Off-diagonals for MPI

IF (l3dim) THEN
  IF (.NOT. mklGeom%lBottom) THEN
    !$OMP PARALLEL PRIVATE(Dtil, Dhat, val)
    !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
    DO ig = 1, ng
      DO ipin = 1, nxy
        Dtil = CMFD%AxDtil(bottom, ipin, 1, ig)
        Dhat = CMFD%AxDhat(bottom, ipin, 1, ig)
        val = - (Dtil + Dhat) * PinVolFm(ipin, 1) / hzfm(1)
        CMFD%AxOffDiag(ipin, bottom, ig) = val
      ENDDO
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL
  ENDIF
  IF (.NOT. mklGeom%lTop) THEN
    !$OMP PARALLEL PRIVATE(Dtil, Dhat, val)
    !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
    DO ig = 1, ng
      DO ipin = 1, nxy
        Dtil = CMFD%AxDtil(top, ipin, nzCMFD, ig)
        Dhat = CMFD%AxDhat(top, ipin, nzCMFD, ig)
        val = - (Dtil + Dhat) * PinVolFm(ipin, nzCMFD) / hzfm(nzCMFD)
        CMFD%AxOffDiag(ipin, top, ig) = val
      ENDDO
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL
  ENDIF
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

SUBROUTINE CMFDSrcUpdt(CMFD, igg)

IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
TYPE(GPinXS_Type), POINTER :: PinXS(:, :)
INTEGER :: idx, ipin, iz, izf, ig, igg, igf, gb, ge
INTEGER :: ng, ngg, nxy, nzCMFD
INTEGER, POINTER :: planeMap(:)
REAL, POINTER :: PinVolFm(:, :)

PinXS => CMFD%GPinXS
planeMap => CMFD%planeMap
PinVolFm => mklGeom%PinVolFm
ng = mklGeom%ng
ngg = mklGeom%ngg
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD

CMFD%src(:, igg) = 0.0

DO ig = 1, ng
  DO izf = 1, nzCMFD
    iz = planeMap(izf)
    !$OMP PARALLEL DO PRIVATE(idx)
    DO ipin = 1, nxy
      idx = ipin + (izf - 1) * nxy
      CMFD%src(idx, igg) = CMFD%src(idx, igg) + PinXS(ipin, iz)%XSP(ig, igg) * PinXS(ipin, iz)%NPhi(ig)
    ENDDO
    !$OMP END PARALLEL DO
  ENDDO
ENDDO

gb = CMFD%InScatRange(1, igg)
ge = CMFD%InScatRange(2, igg)

DO igf = gb, ge
  DO izf = 1, nzCMFD
    iz = planeMap(izf)
    !$OMP PARALLEL DO PRIVATE(idx)
    DO ipin = 1, nxy
      idx = ipin + (izf - 1) * nxy
      CMFD%src(idx, igg) = CMFD%src(idx, igg) + PinXS(ipin, iz)%XSS(igg)%from(igf) * CMFD%phis(ipin, izf, igf)
    ENDDO
    !$OMP END PARALLEL DO
  ENDDO
ENDDO

!$OMP PARALLEL DO PRIVATE(idx) COLLAPSE(2)
DO izf = 1, nzCMFD
  DO ipin = 1, nxy
    idx = ipin + (izf - 1) * nxy
    CMFD%src(idx, igg) = CMFD%src(idx, igg) * PinVolFm(ipin, izf)
  ENDDO
ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE

SUBROUTINE SetMOCPhis(CoreInfo, PinXS, phis, phic, lScat1)
USE TYPEDEF,        ONLY : CoreInfo_Type,       Pin_Type,       Cell_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(GPinXS_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: phis(:, :, :), phic(:, :, :)
LOGICAL :: lScat1

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(superPin_Type), POINTER :: superPin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
INTEGER :: i, j, igg, iz, izf, ifsr, icel, ixy, ixy_map, ipin
INTEGER :: FsrIdxSt
INTEGER :: ngg, nxy, nzCMFD, nLocalFsr
INTEGER :: myzb, myze
INTEGER, POINTER :: pinMap(:), fmRange(:, :)
REAL :: fmult, phisum
REAL, POINTER :: hzfm(:), hz(:)

Pin => CoreInfo%Pin
Cell => CoreInfo%CellInfo
superPin => mklGeom%superPin
pinMap => mklGeom%pinMap
ngg = mklGeom%ngg
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
myzb = mklGeom%myzb
myze = mklGeom%myze
fmRange => mklGeom%fmRange
hzfm => mklGeom%hzfm
hz => mklGeom%hz

!$OMP PARALLEL PRIVATE(FsrIdxSt, ipin, icel, ixy_map, nLocalFsr, phisum, fmult, ifsr)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO igg = 1, ngg
  DO iz = myzb, myze
    DO ixy = 1, nxy
      phisum = 0.0
      ixy_map = pinMap(ixy)
      DO izf = fmRange(iz, 1), fmRange(iz, 2)
        phisum = phisum + mklGammaCMFD%phis(ixy, izf, igg) * (hzfm(izf) / hz(iz))
      ENDDO
      fmult = phisum / PinXS(ixy_map, iz)%GPhi(igg)
      DO j = 1, superPin(ixy_map)%nxy
        ipin = superPin(ixy_map)%pin(j)
        FsrIdxSt = Pin(ipin)%FsrIdxSt
        icel = Pin(ipin)%Cell(iz)
        nLocalFsr = Cell(icel)%nFsr
        DO i = 1, nLocalFsr
          ifsr = FsrIdxSt + i - 1
          IF (lScat1) THEN
            phis(ifsr, igg, iz) = phis(ifsr, igg, iz) * fmult
          ELSE
            phis(igg, ifsr, iz) = phis(igg, ifsr, iz) * fmult
          ENDIF
        ENDDO
        phic(ipin, igg, iz) = phisum
      ENDDO
      mklGammaCMFD%phic(ixy, iz, igg) = phisum
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
INTEGER :: i, igg, iz, izf, ipin, ixy, ixy_map
INTEGER :: ngg, nxy
INTEGER :: myzb, myze
INTEGER, POINTER :: pinMap(:), fmRange(:, :)
REAL :: myphi, neighphi, Dtil, Dhat, Jnet
REAL, POINTER :: hz(:)

Pin => mklGeom%superPin
ngg = mklGeom%ngg
nxy = mklGeom%nxy
myzb = mklGeom%myzb
myze = mklGeom%myze
pinMap => mklGeom%pinMap
fmRange => mklGeom%fmRange
hz => mklGeom%hz

!$OMP PARALLEL PRIVATE(ixy_map, ipin, izf, myphi, neighphi, Dtil, Dhat, Jnet)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO igg = 1, ngg
  DO iz = myzb, myze
    DO ixy = 1, nxy
      ixy_map = pinMap(ixy)
      !--- Axial Source from Bottom
      izf = fmRange(iz, bottom)
      myphi = mklGammaCMFD%phis(ixy, izf, igg)
      IF (iz .EQ. myzb) THEN
        neighphi = mklGammaCMFD%neighphis(ixy, igg, bottom)
      ELSE
        neighphi = mklGammaCMFD%phis(ixy, izf - 1, igg)
      ENDIF
      Dtil = mklGammaCMFD%AxDtil(bottom, ixy, izf, igg)
      Dhat = mklGammaCMFD%AxDhat(bottom, ixy, izf, igg)
      Jnet = - Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
      DO i = 1, Pin(ixy_map)%nxy
        ipin = Pin(ixy_map)%pin(i)
        AxSrc(ipin, igg, iz) = Jnet
      ENDDO
      !--- Axial Source from Top
      izf = fmRange(iz, top)
      myphi = mklGammaCMFD%phis(ixy, izf, igg)
      IF (iz .EQ. myze) THEN
        neighphi = mklGammaCMFD%neighphis(ixy, igg, top)
      ELSE
        neighphi = mklGammaCMFD%phis(ixy, izf + 1, igg)
      ENDIF
      Dtil = mklGammaCMFD%AxDtil(top, ixy, izf, igg)
      Dhat = mklGammaCMFD%AxDhat(top, ixy, izf, igg)
      Jnet = - Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
      DO i = 1, Pin(ixy_map)%nxy
        ipin = Pin(ixy_map)%pin(i)
        AxSrc(ipin, igg, iz) = AxSrc(ipin, igg, iz) + Jnet
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(ipin)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO iz = myzb, myze
  DO igg = 1, ngg
    DO ixy = 1, nxy
      DO i = 1, Pin(ixy)%nxy
        ipin = Pin(ixy)%pin(i)
        AxSrc(ipin, igg, iz) = AxSrc(ipin, igg, iz) / hz(iz)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

AxPXS = 0.0

IF (nTracerCntl%LkgSplitLv .EQ. 0) THEN
  !$OMP PARALLEL PRIVATE(ipin)
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
  DO iz = myzb, myze
    DO igg = 1, ngg
      DO ixy = 1, nxy
        DO i = 1, Pin(ixy)%nxy
          ipin = Pin(ixy)%pin(i)
          IF (AxSrc(ipin, igg, iz) .LT. 0.0) CYCLE
          IF (phic(ipin, igg, iz) .LT. 0.0) CYCLE
          AxPXS(ipin, igg, iz) = AxSrc(ipin, igg, iz) / phic(ipin, igg, iz)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDIF

END SUBROUTINE

SUBROUTINE BiCGSTAB(CMFD, ig)
USE PE_MOD,         ONLY : PE
USE TIMER,          ONLY : nTracer_dclock,      TimeChk
USE MKL_BILU,       ONLY : MKL_SolveILU
IMPLICIT NONE

TYPE(mklCMFD_Type) :: CMFD
INTEGER :: ig

TYPE(CSR_DOUBLE), POINTER :: M, ILU
INTEGER :: i, iter, itermin, itermax
INTEGER :: nr, nc, nnz
INTEGER :: nxy, nzCMFD
REAL :: Tbeg, Tend
REAL :: err, tol
REAL :: rho0, rho1, w0, w1, norm0, norm1, dot0, dot1, alpha, beta
REAL, POINTER :: x(:, :), b(:)
REAL, POINTER :: h(:), s(:), t(:), y(:), z(:)
REAL, POINTER :: r0(:), r1(:), rhat(:), v0(:), v1(:), p0(:), p1(:)
REAL, POINTER :: csrVal(:)
INTEGER, POINTER :: csrRowPtr(:), csrColIdx(:)
INTEGER :: bottomRange(2), topRange(2)
LOGICAL :: lConv

Tbeg = nTracer_dclock(FALSE, FALSE)

M => CMFD%M(ig)
ILU => CMFD%ILU(ig)
x => CMFD%phis(:, :, ig)
b => CMFD%src(:, ig)
csrVal => M%csrVal
csrRowPtr => M%csrRowPtr
csrColIdx => M%csrColIdx

nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
bottomRange = CMFD%bottomRange
topRange = CMFD%topRange

itermin = 5
itermax = 500
tol = 1.0E-06

nr = M%nr
nc = M%nc
nnz = M%nnz

ALLOCATE(h(nc), s(nc), t(nc), y(nc), z(nc))
ALLOCATE(r0(nc), r1(nc), rhat(nc))
ALLOCATE(v0(nc), v1(nc))
ALLOCATE(p0(nc), p1(nc))

!--- 1 ----------------------
CALL MatMulMPIFast(M, CMFD%AxOffDiag(:, :, ig), x, r1, bottomRange, topRange)
CALL vdsub(nc, b, r1, r1)
!--- 2 ----------------------
CALL dcopy(nc, r1, 1, rhat, 1)
!--- 3 ----------------------
norm0 = normMPI(r1, nc, PE%MPI_CMFD_COMM)
rho1 = 1.0; w1 = 1.0; alpha = 1.0
!--- 4 ----------------------
v1 = 0.0; p1 = 0.0
!--- 5 ----------------------
lConv = FALSE; iter = 0
DO WHILE (.NOT. lConv)
  !--- 5.0 ------------------
  iter = iter + 1
  rho0 = rho1; w0 = w1
  CALL dcopy(nc, r1, 1, r0, 1)
  CALL dcopy(nc, p1, 1, p0, 1)
  CALL dcopy(nc, v1, 1, v0, 1)
  !--- 5.1 ------------------
  rho1 = dotMPI(rhat, r0, nc, PE%MPI_CMFD_COMM)
  !--- 5.2 ------------------
  beta = (rho1 / rho0) * (alpha / w0)
  !--- 5.3 ------------------
  CALL dcopy(nc, r0, 1, p1, 1)
  CALL daxpy(nc, -w0, v0, 1, p0, 1)
  CALL daxpy(nc, beta, p0, 1, p1, 1)
  !--- 5.4 ------------------
  CALL MKL_SolveILU(ILU, y, p1)
  !--- 5.5 ------------------
  CALL MatMulMPIFast(M, CMFD%AxOffDiag(:, :, ig), y, v1, bottomRange, topRange)
  !--- 5.6 ------------------
  dot1 = dotMPI(rhat, v1, nc, PE%MPI_CMFD_COMM)
  alpha = rho1 / dot1
  !--- 5.7 ------------------
  CALL dcopy(nc, x, 1, h, 1)
  CALL daxpy(nc, alpha, y, 1, h, 1)
  !--- 5.9 ------------------
  CALL dcopy(nc, r0, 1, s, 1)
  CALL daxpy(nc, -alpha, v1, 1, s, 1)
  !--- 5.10 -----------------
  CALL MKL_SolveILU(ILU, z, s)
  !--- 5.11 -----------------
  CALL MatMulMPIFast(M, CMFD%AxOffDiag(:, :, ig), z, t, bottomRange, topRange)
  !--- 5.12 -----------------
  dot0 = dotMPI(t, s, nc, PE%MPI_CMFD_COMM)
  dot1 = dotMPI(t, t, nc, PE%MPI_CMFD_COMM)
  w1 = dot0 / dot1
  !--- 5.13 -----------------
  CALL dcopy(nc, h, 1, x, 1)
  CALL daxpy(nc, w1, z, 1, x, 1)
  !--- 5.15 -----------------
  CALL dcopy(nc, s, 1, r1, 1)
  CALL daxpy(nc, -w1, t, 1, r1, 1)
  !--- Convergence ----------
  IF (mklCntl%lDcpl) THEN
    norm1 = dnrm2(nc, r1, 1)
  ELSE
    norm1 = normMPI(r1, nc, PE%MPI_CMFD_COMM)
  ENDIF
  err = norm1 / norm0; lConv = (err .LE. tol) .AND. (iter .GE. itermin)
  lConv = lConv .OR. (iter .GE. itermax)
ENDDO
print*, norm1, norm0, iter

DEALLOCATE(h, s, t, y, z)
DEALLOCATE(r0, r1, rhat)
DEALLOCATE(v0, v1)
DEALLOCATE(p0, p1)

Tend = nTracer_dclock(FALSE, FALSE)
TimeChk%AxBTime = TimeChk%AxBTime + (Tend - Tbeg)

END SUBROUTINE

END MODULE
#endif
#endif
