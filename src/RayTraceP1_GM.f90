#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceP1_GM(RayInfo, CoreInfo, phis1g, phim1g, PhiAngIn1g, xst1g, src1g, srcm1g, jout1g, iz, ljout, ScatOd, FastMocLv)

USE TIMER
USE ALLOCS
USE OMP_LIB
USE PARAM,   ONLY : FALSE, ZERO, ONE, FORWARD, BACKWARD
USE TYPEDEF, ONLY : RayInfo_Type, CoreInfo_type, Pin_Type, Cell_Type
USE Moc_Mod, ONLY : TrackingDat, Comp, SrcAng1g1, SrcAng1g2, wtang, mwt
USE PE_MOD,  ONLY : PE
USE CNTL,    ONLY : nTracerCntl

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:)     :: phis1g, xst1g, src1g
REAL, POINTER, DIMENSION(:,:)   :: PhiAngIn1g, srcm1g, phim1g
REAL, POINTER, DIMENSION(:,:,:) :: jout1g

INTEGER :: iz, ScatOd
LOGICAL :: ljout, lAFSS
INTEGER, OPTIONAL :: FastMocLv
! ----------------------------------------------------
TYPE (Cell_Type), POINTER, DIMENSION(:) :: Cell
TYPE (Pin_Type),  POINTER, DIMENSION(:) :: Pin

INTEGER :: nAziAng, nPolarAng, nFsr, nxy, nThr
INTEGER :: ithr, FsrIdxSt, icel, iazi, ipol, iod, iRotRay, ifsr, jfsr, ixy, krot
REAL :: wttmp, tmpsrc, ONETHREE, ONEFIVE, ONESEVEN
! ----------------------------------------------------

nAziAng   = RayInfo%nAziAngle
nPolarAng = RayInfo%nPolarAngle

nFsr = CoreInfo%nCoreFsr
nxy  = CoreInfo%nxy

nThr = PE%nThread

lAFSS = nTracerCntl%lAFSS
! ----------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ithr, iazi, ifsr, ipol, tmpsrc)
ithr = omp_get_thread_num() + 1

DO iazi = 1, nAziAng
  DO ifsr = PE%myOmpFsrBeg(ithr), PE%myOmpFsrEnd(ithr)
    DO ipol = 1, nPolarAng
      SrcAng1g1(ipol, ifsr, iazi) = src1g(ifsr)
      SrcAng1g2(ipol, ifsr, iazi) = src1g(ifsr)
      
      tmpsrc = comp(1, ipol, iazi) * srcm1g(1, ifsr) + comp(2, ipol, iazi) * srcm1g(2, ifsr)
      
      SrcAng1g1(ipol, ifsr, iazi) = SrcAng1g1(ipol, ifsr, iazi) + tmpsrc
      SrcAng1g2(ipol, ifsr, iazi) = SrcAng1g2(ipol, ifsr, iazi) - tmpsrc
      
      IF (scatod .LT. 2) CYCLE
      
      tmpsrc = comp(3, ipol, iazi) * srcm1g(3, ifsr) + comp(4, ipol, iazi) * srcm1g(4, ifsr) + comp(5, ipol, iazi) * srcm1g(5, ifsr)
      
      SrcAng1g1(ipol, ifsr, iazi) = SrcAng1g1(ipol, ifsr, iazi) +  tmpsrc
      SrcAng1g2(ipol, ifsr, iazi) = SrcAng1g2(ipol, ifsr, iazi) +  tmpsrc
      
      IF (scatod .LT. 3) CYCLE
      
      tmpsrc = comp(6, ipol, iazi) * srcm1g(6, ifsr) + comp(7, ipol, iazi) * srcm1g(7, ifsr) + comp(8, ipol, iazi) * srcm1g(8, ifsr) + comp(9, ipol, iazi) * srcm1g(9, ifsr)
      
      SrcAng1g1(ipol, ifsr, iazi) = SrcAng1g1(ipol, ifsr, iazi) + tmpsrc
      SrcAng1g2(ipol, ifsr, iazi) = SrcAng1g2(ipol, ifsr, iazi) - tmpsrc
     END DO
  END DO
END DO
!$OMP END PARALLEL
! ----------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ithr, krot, iRotRay)
ithr = omp_get_thread_num() + 1

TrackingDat(ithr)%phis1g = ZERO
TrackingDat(ithr)%phim1g = ZERO
IF (ljout) TrackingDat(ithr)%jout1g = ZERO
IF (lAFSS) TrackingDat(ithr)%phia1g = ZERO

TrackingDat(ithr)%PhiAngIn1g => PhiAngIn1g
TrackingDat(ithr)%src1g      => src1g
TrackingDat(ithr)%xst1g      => xst1g
TrackingDat(ithr)%SrcAng1g1  => SrcAng1g1
TrackingDat(ithr)%SrcAng1g2  => SrcAng1g2

DO krot = 1, 2
  IF (nTracerCntl%lHex) THEN
    !$OMP DO SCHEDULE(GUIDED)
    DO iRotRay = 1, RayInfo%nRotRay
      CALL HexTrackRotRayP1_GM(RayInfo, CoreInfo, TrackingDat(ithr), lJout, iRotRay, iz, krot, ScatOd, FastMocLv, lAFSS)
    END DO
    !$OMP END DO NOWAIT
  ELSE
    !$OMP DO SCHEDULE(GUIDED)
    DO iRotRay = 1, RayInfo%nRotRay
      CALL RecTrackRotRayP1_GM(RayInfo, CoreInfo, TrackingDat(ithr), lJout, iRotRay, iz, krot, ScatOd, FastMocLv, lAFSS)
    END DO
    !$OMP END DO NOWAIT
  END IF
END DO
!$OMP END PARALLEL
! ----------------------------------------------------
IF (.NOT. lAFSS) THEN
  phis1g = ZERO
  phim1g = ZERO

  DO ithr = 1, nThr
    phis1g = phis1g + TrackingDat(ithr)%phis1g
    phim1g = phim1g + TrackingDat(ithr)%phim1g
  END DO
ELSE
  DO ithr = 1, nThr
    DO iazi = 1, RayInfo%nAziAngle
      DO ipol = 1, RayInfo%nPolarAngle
        DO ifsr = 1, nFsr
          phis1g(ifsr) = phis1g(ifsr) + wtang(ipol, iAzi) * (TrackingDat(ithr)%phia1g(FORWARD, ipol, iazi, ifsr) + TrackingDat(ithr)%phia1g(BACKWARD, ipol, iazi, ifsr))
          
          phim1g(1:2, ifsr) = phim1g(1:2, ifsr) + mwt(1:2, ipol, iazi) * (TrackingDat(ithr)%phia1g(FORWARD, ipol, iazi, ifsr) - TrackingDat(ithr)%phia1g(BACKWARD, ipol, iazi, ifsr))
          
          IF (scatod .LT. 2) CYCLE
          
          phim1g(3:5, ifsr) = phim1g(3:5, ifsr) + mwt(3:5, ipol, iazi) * (TrackingDat(ithr)%phia1g(FORWARD, ipol, iazi, ifsr) + TrackingDat(ithr)%phia1g(BACKWARD, ipol, iazi, ifsr))
          
          IF (scatod .LT. 3) CYCLE
          
          phim1g(6:9, ifsr) = phim1g(6:9, ifsr) + mwt(6:9, ipol, iazi) * (TrackingDat(ithr)%phia1g(FORWARD, ipol, iazi, ifsr) - TrackingDat(ithr)%phia1g(BACKWARD, ipol, iazi, ifsr))
        END DO
      END DO
    END DO
  END DO
  !phia1g = ZERO
  !
  !DO ithr = 1, nThr
  !  phia1g = phia1g + TrackingDat(ithr)%phia1g
  !END DO
END IF

IF (ljout) THEN
  jout1g = ZERO
  
  DO ithr = 1, nThr
    jout1g = jout1g + TrackingDat(ithr)%jout1g
  END DO
END IF
! ----------------------------------------------------
ONETHREE = ONE / 3._8
ONEFIVE  = ONE / 5._8
ONESEVEN = ONE / 7.

Cell => CoreInfo%CellInfo
Pin  => CoreInfo%Pin

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ithr, ixy, FsrIdxSt, icel, ifsr, jfsr, wttmp)
ithr = omp_get_thread_num() + 1

DO ixy = PE%myOmpNxyBeg(ithr), PE%myOmpNxyEnd(ithr)
  FsrIdxSt = Pin(ixy)%FsrIdxSt
  icel     = Pin(ixy)%Cell(iz)
  
  DO ifsr = 1, Cell(icel)%nFsr
    jfsr  = FsrIdxSt + ifsr - 1
    wttmp = ONE / xst1g(jfsr) / Cell(icel)%vol(ifsr)
    
    phis1g(jfsr) = phis1g(jfsr) * wttmp + src1g(jfsr)
    
    phim1g(1:2, jfsr) = phim1g(1:2, jfsr) * wttmp + srcm1g(1:2, jfsr) * ONETHREE
    
    IF (scatod .LT. 2) CYCLE
    
    phim1g(3:5, jfsr) = phim1g(3:5, jfsr) * wttmp + srcm1g(3:5, jfsr) * ONEFIVE
    
    IF (scatod .LT. 3) CYCLE
    
    phim1g(6:9, jfsr) = phim1g(6:9, jfsr) * wttmp + srcm1g(6:9, jfsr) * ONESEVEN
  END DO
END DO
!$OMP END PARALLEL

NULLIFY (Cell)
NULLIFY (Pin)
! ----------------------------------------------------

END SUBROUTINE RayTraceP1_GM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayP1_GM(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, krot, ScatOd, FastMocLv, lAFSS)

USE PARAM,   ONLY : TRUE, FALSE, ZERO, FORWARD, BACKWARD
USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type, Pin_Type, Asy_Type, AsyInfo_Type, PinInfo_Type, Cell_Type, AziAngleInfo_Type, PolarAngle_Type, ModRayInfo_type, AsyRayInfo_type,  &
                    CoreRayInfo_Type, RotRayInfo_Type, CellRayInfo_type, FastCoreRayDat_Type, TrackingDat_Type, FastRaySegDat_Type
USE Moc_Mod, ONLY : nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay

IMPLICIT NONE

TYPE (RayInfo_Type),     INTENT(INOUT) :: RayInfo
TYPE (CoreInfo_Type),    INTENT(INOUT) :: CoreInfo
TYPE (TrackingDat_Type), INTENT(INOUT) :: TrackingDat

LOGICAL, INTENT(IN) :: ljout, lAFSS
INTEGER, INTENT(IN) :: irotray, iz, krot, ScatOd, FastMocLv
! ----------------------------------------------------
TYPE (Pin_Type),            POINTER, DIMENSION(:) :: Pin
TYPE (Asy_Type),            POINTER, DIMENSION(:) :: Asy
TYPE (PinInfo_Type),        POINTER, DIMENSION(:) :: PinInfo
TYPE (Cell_Type),           POINTER, DIMENSION(:) :: Cell
TYPE (AziAngleInfo_Type),   POINTER, DIMENSION(:) :: AziAng
TYPE (PolarAngle_Type),     POINTER, DIMENSION(:) :: PolarAng
TYPE (AsyRayInfo_type),     POINTER, DIMENSION(:) :: AsyRay
TYPE (CoreRayInfo_Type),    POINTER, DIMENSION(:) :: CoreRay
TYPE (RotRayInfo_Type),     POINTER, DIMENSION(:) :: RotRay

TYPE (FastCoreRayDat_Type), POINTER :: FastRay
TYPE (CellRayInfo_Type),    POINTER :: CellRay
TYPE (CellRayInfo_Type),    POINTER :: CellRay1D
TYPE (FastRaySegDat_Type),  POINTER :: FastRaySeg

INTEGER, POINTER, DIMENSION(:)   :: LocalFsrIdx
INTEGER, POINTER, DIMENSION(:,:) :: FsrIdx, ExpAppIdx

REAL, POINTER, DIMENSION(:)       :: LenSeg, phis1g, xst1g
REAL, POINTER, DIMENSION(:,:)     :: OptLenList, PhiAngOut1g, PhiAngIn1g, ExpA, ExpB, wtang, phim1g
REAL, POINTER, DIMENSION(:,:,:)   :: ExpAppPolar, jout1g, wtsurf, mwt, mwt2, SrcAng1g1, SrcAng1g2
REAL, POINTER, DIMENSION(:,:,:,:) :: phia1g

INTEGER :: mp(2)
INTEGER :: iazi, ipol, iCoreRay, iasyray, iceray, irayseg, itype, idir, nRotRay, nCoreRay, nAsyRay, nPinRay, nRaySeg, FsrIdxSt, nPolarAng, nAziAng, nPhiAngSv, od
INTEGER :: ipin, icel, iasy, ireg, isurf, irsegidx, icellrayidx, PhiAnginSvIdx, PhiAngOutSvIdx, nFsr, nxy, i, j, k, l, m, jbeg, jend, jinc, ir, ir1, iod, ibcel

INTEGER, DIMENSION(nMaxCellRay, nMaxCoreRay, 2) :: CellRayIdxSt, SurfIdx
INTEGER, DIMENSION(nMaxCellRay, nMaxCoreRay) :: PinIdx
INTEGER, DIMENSION(nMaxCoreRay) :: nTotRaySeg, nTotCellRay

REAL :: tau, phiobd, phid, wt, wt2(4)
REAL, ALLOCATABLE :: phiobdPolar(:)

LOGICAL :: lFast

DATA mp /2, 1/
! ----------------------------------------------------

lFast = FALSE
IF (FastMocLv .GT. 0) lFast = TRUE

AziAng   => RayInfo%AziAngle
PolarAng => RayInfo%PolarAngle
AsyRay   => RayInfo%AsyRay
CoreRay  => RayInfo%CoreRay
RotRay   => RayInfo%RotRay
nAziAng   = RayInfo%nAziAngle
nPolarAng = RayInfo%nPolarAngle
nPhiAngSv = RayInfo%nPhiAngSv

Asy     => CoreInfo%Asy
Pin     => CoreInfo%Pin
PinInfo => CoreInfo%Pininfo
Cell    => CoreInfo%CellInfo
nFsr     = CoreInfo%nCoreFsr
nxy      = CoreInfo%nxy

FsrIdx      => TrackingDat%FsrIdx
ExpAppIdx   => TrackingDat%ExpAppIdx
OptLenList  => TrackingDat%OptLenList
ExpAppPolar => TrackingDat%ExpAppPolar
Phis1g      => TrackingDat%phis1g
phim1g      => TrackingDat%phim1g
xst1g       => TrackingDat%xst1g
jout1g      => TrackingDat%jout1g
PhiAngOut1g => TrackingDat%PhiAngOut1g
PhiAngIn1g  => TrackingDat%phiAngIn1g
ExpA        => TrackingDat%ExpA
ExpB        => TrackingDat%ExpB
wtang       => TrackingDat%wtang
wtsurf      => TrackingDat%wtsurf
mwt         => TrackingDat%mwt
mwt2        => TrackingDat%mwt2
SrcAng1g1   => TrackingDat%SrcAng1g1
SrcAng1g2   => TrackingDat%SrcAng1g2

IF (lAFSS) phia1g => TrackingDat%phia1g
! ----------------------------------------------------
OD = 2
IF (ScatOd .EQ. 2) OD = 5
IF (ScatOd .EQ. 3) OD = 9

i = iRotRay
nCoreRay = RotRay(irotRay)%nRay

IF (.NOT. lFast) THEN
  DO j = 1, nCoreRay
    irsegidx    = 0
    icellrayidx = 0
    
    iCoreRay = RotRay(iRotRay)%RayIdx(j)
    nAsyRay  = CoreRay(iCoreRay)%nRay
    
    DO k = 1, nAsyRay
      iasyray = CoreRay(iCoreRay)%AsyRayIdx(k)
      iasy    = CoreRay(iCoreRay)%AsyIdx(k)
      
      IF (iasy .EQ. 0) CYCLE ! Skip Dummy Assembly
      
      nPinRay = AsyRay(iAsyRay)%nCellRay
      itype   = Asy(iasy)%PartialAsyFlag
      
      DO l = 1, nPinRay
        ipin   = AsyRay(iAsyRay)%PinIdx(l) ! Local Pin Idx(within Assembly)
        iceray = AsyRay(iAsyRay)%PinRayIdx(l)
        
        ipin = Asy(iAsy)%GlobalPinIdx(ipin) ! Global Pin Index
        
        icel     = Pin(ipin)%Cell(iz)
        FsrIdxSt = Pin(ipin)%FsrIdxSt
        
        ibcel    = Cell(icel)%basecellstr
        CellRay => Cell(ibcel)%CellRay(iceray)
        
        icellrayidx = icellrayidx + 1
        
        PinIdx      (icellrayidx, j)    = ipin
        CellRayIdxSt(icellrayidx, j, 2) = irsegidx + 1
        
        nRaySeg      = CellRay%nSeg
        LocalFsrIdx => CellRay%LocalFsrIdx
        LenSeg      => CellRay%LenSeg
        
        DO iRaySeg = 1, nRaySeg
          ireg = FsrIdxSt + CellRay%LocalFsrIdx(iRaySeg) - 1
          
          tau = -LenSeg(iRaySeg) * xst1g(ireg)
          tau = -CellRay%LenSeg(iRaySeg) * xst1g(ireg)
          
          irsegidx = irsegidx + 1
          
          FsrIdx    (irsegidx, j) = ireg
          OptLenList(irsegidx, j) = tau
          ExpAppIdx (irsegidx, j) = max(INT(tau), -40000)
          ExpAppIdx (irsegidx, j) = min(0, ExpAppIdx(irsegidx, j))
        END DO
        CellRayIdxSt(icellrayidx, j, 1) = irsegidx
        SurfIdx     (icellRayIdx, j, 1) = AsyRay(iAsyRay)%PinRaySurf(2, l) ! OutSurface
        SurfIdx     (icellRayIdx, j, 2) = AsyRay(iAsyRay)%PinRaySurf(1, l) ! Insurface
      END DO
    END DO
    
    nTotRaySeg(j) = irsegidx
    nTotCellRay(j) = icellRayIdx
  END DO
ELSE IF (FastMocLV .EQ. 1) THEN
  FastRay   => RayInfo%FastCoreRayDat(i, iz)
  CellRay1D => RayInfo%CellRay1D
  
  DO j = 1, nCoreRay
    irsegidx = 0
    
    nTotRaySeg (j) = FastRay%nTotRaySeg(j)
    nTotCellRay(j) = FastRay%nTotCellRay(j)
    
    DO l = 1, FastRay%nTotCellRay(j)
      PinIdx      (l, j)    = FastRay%PinIdx      (l, j)
      CellRayIdxSt(l, j, 1) = FastRay%CellRayIdxSt(l, j, 1)
      CellRayIdxSt(l, j, 2) = FastRay%CellRayIdxSt(l, j, 2)
      SurfIdx     (l, j, 1) = FastRay%SurfIdx     (l, j, 1)
      SurfIdx     (l, j, 2) = FastRay%SurfIdx     (l, j, 2)
      
      ipin = PinIdx(l, j)
      
      icel     = Pin(ipin)%Cell(iz)
      FsrIdxSt = Pin(ipin)%FsrIdxSt
      
      DO k = FastRay%Ray1DIdx(1, l, j), FastRay%Ray1DIdx(2, l, j)
        irsegidx = irsegidx + 1
        
        ireg = FsrIdxSt + CellRay1D%LocalFsrIdx(K) - 1
        
        FsrIdx(irsegidx, j) = ireg
        
        tau = -CellRay1D%LenSeg(k) * xst1g(ireg)
        
        OptLenList(irsegidx, j) = tau
        ExpAppIdx (irsegidx, j) = max(INT(tau), -40000)
        ExpAppIdx (irsegidx, j) = min(0, ExpAppIdx(irsegidx, j))
      END DO
    END DO
  END DO
ELSE IF (FastMocLv .EQ. 2) THEN
  FastRay   => RayInfo%FastCoreRayDat(i, iz)
  CellRay1D => RayInfo%CellRay1D
  
  DO j = 1, nCoreRay
    nTotRaySeg (j) = FastRay%nTotRaySeg (j)
    nTotCellRay(j) = FastRay%nTotCellRay(j)
    
    DO l = 1, FastRay%nTotCellRay(j)
      PinIdx      (l, j) = FastRay%PinIdx        (l, j)
      CellRayIdxSt(l, j, 1) = FastRay%CellRayIdxSt(l, j, 1)
      CellRayIdxSt(l, j, 2) = FastRay%CellRayIdxSt(l, j, 2)
      SurfIdx     (l, j, 1) = FastRay%SurfIdx     (l, j, 1)
      SurfIdx    (l, j, 2) = FastRay%SurfIdx      (l, j, 2)
    END DO
    
    FastRaySeg => RayInfo%FastCoreRayDat(i, iz)%RaySeg(j)
    
    DO l = 1, FastRay%nTotRaySeg(j)
      ireg = FastRaySeg%FsrIdx(l)
      
      FsrIdx(l, j) = ireg
      
      tau = -FastRaySeg%LenSeg(l) * xst1g(ireg)
      
      OptLenList(l, j) = tau
      ExpAppIdx (l, j) = max(INT(tau), -40000)
      ExpAppIdx (l, j) = min(0, ExpAppIdx(l, j))
    END DO
  END DO
END IF
! ----------------------------------------------------
DO j = 1, nCoreRay
  DO l = 1, nTotRaySeg(j)
    DO ipol = 1, nPolarAng
      ExpAppPolar(ipol, l, j) = ExpA(ExpAppIdx(l, j), ipol) * optlenlist(l, j) + ExpB(ExpAppIdx(l, j), ipol)
    END DO
  END DO
END DO
! ----------------------------------------------------
ALLOCATE (phiobdPolar(nPolarAng))

PhiAnginSvIdx  = RayInfo%PhiAngInSvIdx(iRotRay ,krot)
PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay, krot)
phiobdPolar(:) = PhiAngIn1g(:,PhiAnginSvIdx)

jinc = 1; jbeg = 1; jend = nCoreRay
IF (krot .EQ. 2) THEN ! Backward Sweep
  jinc = -1; jbeg = nCoreRay; jend = 1
END IF
! ----------------------------------------------------
DO j = jbeg, jend, jinc
  iazi = CoreRay(RotRay(irotray)%RayIdx(j))%iang
  
  idir = RotRay(i)%DIR(j)
  IF (krot .EQ. 2) idir = mp(idir) ! Reverse the sweep direction
  
  nRaySeg = nTotRaySeg(j)
  
  IF (idir .EQ. 1) THEN ! Forward Sweep
    PhiAngOut1g(:, 1) = phiobdPolar(:)
    
    DO ir = 1, nRaySeg
      ireg = FsrIdx(ir, j)
      
      DO ipol = 1, nPolarAng
        wt = wtang(ipol, iazi)
        
        phid = (PhiAngOut1g(ipol, ir) - SrcAng1g1(ipol, ireg, iazi)) * ExpAppPolar(ipol, ir, j)
        
        PhiAngOut1g(ipol, ir+1) = PhiAngOut1g(ipol, ir) - phid
        
        IF (lAFSS) THEN
          phia1g(FORWARD, ipol, iazi, ireg) = phia1g(FORWARD, ipol, iazi, ireg) + phid
        ELSE
          phis1g(ireg) = phis1g(ireg) + wt*phid
          
          DO iod = 1, od
            phim1g(iod, ireg) = phim1g(iod, ireg) + mwt(iod, ipol, iazi) * phid
          END DO
        END IF
      END DO
    END DO
    
    phiobdPolar(:) = PhiAngOut1g(:, nRaySeg+1)
    
    ! Surface
    IF (ljout) THEN
      DO ir = 1, nTotCellRay(j)
        DO ipol = 1, nPolarANg
          wt = wtang(ipol, iazi)
          
          wt2(1:4) = wtsurf(ipol, iazi, 1:4)
          
          icel = PinIdx(ir, j); isurf = SurfIdx(ir, j, 1)
          
          Jout1g(2, isurf, icel) = Jout1g(2, isurf, icel) + PhiAngOut1g(ipol, CellRayIdxSt(ir, j, 1)+1) * wt
          Jout1g(3, isurf, icel) = Jout1g(3, isurf, icel) + PhiAngOut1g(ipol, CellRayIdxSt(ir, j, 1)+1) * wt2(isurf)
          
          isurf = SurfIdx(ir, j, 2)
          
          Jout1g(1, isurf, icel) = Jout1g(1, isurf, icel) + PhiAngOut1g(ipol, CellRayIdxSt(ir, j, 2)) * wt
          Jout1g(3, isurf, icel) = Jout1g(3, isurf, icel) + PhiAngOut1g(ipol, CellRayIdxSt(ir, j, 2)) * wt2(isurf)
        END DO
      END DO
    END IF
  ELSE
    PhiAngOut1g(:, nRaySeg+2) = phiobdPolar(:)
    
    ir = nRaySeg + 1
    
    DO ir1 = 1, nRaySEg
      ir   = ir - 1
      ireg = FsrIdx(ir, j)
      
      DO ipol = 1, nPolarAng
        wt = wtang(ipol, iazi)
        
        phid = (PhiAngOut1g(ipol, ir + 2) - SrcAng1g2(ipol, ireg, iazi)) * ExpAppPolar(ipol, ir, j)
        
        PhiAngOut1g(ipol, ir+1) = PhiAngOut1g(ipol, ir + 2) - phid
        
        IF (lAFSS) THEN
          phia1g(BACKWARD, ipol, iazi, ireg) = phia1g(BACKWARD, ipol, iazi, ireg) + phid
        ELSE
          phis1g(ireg) = phis1g(ireg) + wt * phid
          
          DO iod = 1, od
            phim1g(iod, ireg) = phim1g(iod, ireg) + mwt2(iod, ipol, iazi) * phid
          END DO
        END IF
      END DO
    END DO
    
    phiobdPolar(:) = PhiAngOut1g(:, 2)
    
    IF (lJout) THEN
      DO ir = 1, nTotCellRay(j)
        DO ipol = 1, nPolarAng
          wt = wtang(ipol, iazi)
          
          wt2(1:4) = wtsurf(ipol, iazi, 1:4)
          
          icel = PinIdx(ir, j); isurf = SurfIdx(ir, j, 2)
          
          Jout1g(2, isurf, icel) = Jout1g(2, isurf, icel) + PhiAngOut1g(ipol, CellRayIdxSt(ir, j, 2)+1) * wt
          Jout1g(3, isurf, icel) = Jout1g(3, isurf, icel) + PhiAngOut1g(ipol, CellRayIdxSt(ir, j, 2)+1) * wt2(isurf)
          
          isurf = SurfIdx(ir, j, 1)
          
          Jout1g(1, isurf, icel) = Jout1g(1, isurf, icel) + PhiAngOut1g(ipol, CellRayIdxSt(ir, j, 1)+2) * wt
          Jout1g(3, isurf, icel) = Jout1g(3, isurf, icel) + PhiAngOut1g(ipol, CellRayIdxSt(ir, j, 1)+2) * wt2(isurf)
        END DO
      END DO
    END IF
  END IF
END DO

PhiAngIn1g(:, PhiAngOutSvIdx) = phiobdPolar(:)
! ----------------------------------------------------
DEALLOCATE(phiobdPolar)

! Ray
NULLIFY (AziAng)
NULLIFY (PolarAng)
NULLIFY (AsyRay)
NULLIFY (CoreRay)
NULLIFY (RotRay)
NULLIFY (CellRay)

! Geo.
NULLIFY (Asy)
NULLIFY (Pin)
NULLIFY (PinInfo)
NULLIFY (Cell)

! Local
NULLIFY (FsrIdx)
NULLIFY (ExpAppIdx)
NULLIFY (OptLenList)
NULLIFY (ExpAppPolar)
NULLIFY (LenSeg)
NULLIFY (LocalFsrIdx)
NULLIFY (Phis1g)
NULLIFY (phim1g)
NULLIFY (xst1g)
NULLIFY (jout1g)
NULLIFY (PhiAngOut1g)
NULLIFY (PhiAngIn1g)
NULLIFY (ExpA)
NULLIFY (ExpB)
NULLIFY (wtang)
NULLIFY (wtsurf)
NULLIFY (mwt)
NULLIFY (mwt2)
NULLIFY (SrcAng1g1)
NULLIFY (SrcAng1g2)

IF (lAFSS) NULLIFY (phia1g)
! ----------------------------------------------------

END SUBROUTINE RecTrackRotRayP1_GM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayP1_GM(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, krot, ScatOd, FastMocLv, lAFSS)

USE param,   ONLY : FORWARD, BACKWARD
USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, TrackingDat_Type, Pin_Type
USE Moc_Mod, ONLY : nMaxCellRay, nMaxCoreRay
USE HexType, ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexCoreRay, Type_HexRotRay
USE HexData, ONLY : hAsy, haRay, hcRay, hRotRay, hAsyTypInfo

IMPLICIT NONE

TYPE (RayInfo_Type),     INTENT(INOUT) :: RayInfo
TYPE (CoreInfo_Type),    INTENT(INOUT) :: CoreInfo
TYPE (TrackingDat_Type), INTENT(INOUT) :: TrackingDat

LOGICAL, INTENT(IN) :: ljout, lAFSS
INTEGER, INTENT(IN) :: irotray, iz, krot, ScatOd
INTEGER, INTENT(IN) :: FastMocLv
! ----------------------------------------------------
INTEGER :: iazi, ipol, iaRay, jaRay, iAsy, iSurf, icRay, jbeg, jend, jinc, ihpRay, ifsr, iRaySeg, iRaySeg1, iCel, iReg
INTEGER :: PhiAnginSvIdx, PhiAngOutSvIdx, iGeoTyp, iAsyTyp, jhPin, icBss, jcBss, jcRay, nod, iod, irsegidx, icellrayidx
INTEGER :: nCoreRay, nAsyRay, nPolarAng, nRaySeg

INTEGER :: CellRayIdxSt(nMaxCellRay, nMaxCoreRay, 2)
INTEGER :: PinIdx      (nMaxCellRay, nMaxCoreRay)
INTEGER :: SurfIdx     (nMaxCellRay, nMaxCoreRay, 2)
INTEGER :: nTotRaySeg  (nMaxCoreRay)
INTEGER :: nTotCellRay (nMaxCoreRay)

INTEGER, POINTER, DIMENSION(:,:) :: FsrIdx, ExpAppIdx

REAL :: tau, phid, wt

REAL :: locphiout(RayInfo%nPolarAngle)

REAL, POINTER, DIMENSION(:)       :: phis1g, xst1g
REAL, POINTER, DIMENSION(:,:)     :: ExpA, ExpB, wtang, phim1g, OptLenList, PhiAngOut1g
REAL, POINTER, DIMENSION(:,:,:)   :: jout1g, mwt, mwt2, SrcAng1g1, SrcAng1g2, ExpAppPolar
REAL, POINTER, DIMENSION(:,:,:,:) :: phia1g

TYPE (Pin_Type), POINTER, DIMENSION(:) :: Pin

TYPE (Type_HexAsyRay),  POINTER :: haRay_Loc
TYPE (Type_HexCelRay),  POINTER :: CelRay_Loc
TYPE (Type_HexRotRay),  POINTER :: hRotRay_Loc
! ----------------------------------------------------

nPolarAng      = RayInfo%nPolarAngle
PhiAngInSvIdx  = RayInfo%PhiAngInSvIdx (iRotRay, krot)
PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay, krot)

hRotRay_Loc => hRotRay(iRotRay)
nCoreRay     = hRotRay_Loc%ncRay

Pin => CoreInfo%Pin

phis1g         => TrackingDat%phis1g
phim1g         => TrackingDat%phim1g
xst1g          => TrackingDat%xst1g
jout1g         => TrackingDat%jout1g
ExpA           => TrackingDat%ExpA
ExpB           => TrackingDat%ExpB
wtang          => TrackingDat%wtang
mwt            => TrackingDat%mwt
mwt2           => TrackingDat%mwt2
SrcAng1g1      => TrackingDat%SrcAng1g1
SrcAng1g2      => TrackingDat%SrcAng1g2
locphiout       = TrackingDat%PhiAngIn1g(:, PhiAnginSvIdx)
FsrIdx         => TrackingDat%FsrIdx
ExpAppIdx      => TrackingDat%ExpAppIdx
OptLenList     => TrackingDat%OptLenList
ExpAppPolar    => TrackingDat%ExpAppPolar
PhiAngOut1g => TrackingDat%PhiAngOut1g

IF (lAFSS) phia1g => TrackingDat%phia1g

nod = 2
IF (ScatOd .EQ. 2) nod = 5
IF (ScatOd .EQ. 3) nod = 9
! ----------------------------------------------------
DO icRay = 1, nCoreRay
  jcRay   = abs(hRotRay_Loc%cRayIdx(icRay))
  nAsyRay = hcRay(jcRay)%nmRay
  
  irSegIdx    = 0
  iCellRayIdx = 0
  
  DO iaRay = 1, nAsyRay
    jaRay   = hcRay(jcRay)%mRayIdx(iaRay)
    iAsy    = hcRay(jcRay)%AsyIdx (iaRay)
    iAsyTyp = hAsy(iAsy)%AsyTyp
    iGeoTyp = hAsy(iAsy)%GeoTyp
    icBss   = hAsyTypInfo(iAsyTyp)%iBss
    
    haRay_Loc => haRay(iGeoTyp, icBss, jaRay)
    
    DO ihpRay = 1, haRay_Loc%nhpRay
      iCellRayIdx = iCellRayIdx + 1
      
      jhPin = haRay_Loc%CelRay(ihpRay)%hPinIdx
      jhPin = hAsy(iAsy)%PinIdxSt + hAsyTypInfo(iAsyTyp)%PinLocIdx(iGeoTyp, jhPin) - 1
      jcBss = Pin(jhPin)%hCelGeo(iz)
      
      CelRay_Loc => haRay(iGeoTyp, jcBss, jaRay)%CelRay(ihpRay)
      
      ! Start of Cell
      CellRayIdxSt(iCellRayIdx, icRay, 2) = irSegIdx + 1
      SurfIdx     (iCellRayIdx, icRay, 2) = CelRay_Loc%hSufIdx(1) ! y : Small
      
      DO iRaySeg = 1, CelRay_Loc%nSegRay
        irSegIdx = irSegIdx + 1
        iReg     =  CelRay_Loc%MshIdx(iRaySeg) + Pin(jhPin)%FsrIdxSt - 1
        tau      = -CelRay_Loc%SegLgh(iRaySeg) * xst1g(iReg) ! Optimum Length
        
        FsrIdx    (irSegIdx, icRay) = iReg
        OptLenList(irSegIdx, icRay) = tau
        ExpAppIdx (irSegIdx, icRay) = max(INT(tau), -40000)
        ExpAppIdx (irSegIdx, icRay) = min(0, ExpAppIdx(irSegIdx, icRay))
      END DO
      
      ! End of Cel
      CellRayIdxSt(iCellRayIdx, icRay, 1) = irSegIdx
      SurfIdx     (iCellRayIdx, icRay, 1) = CelRay_Loc%hSufIdx(2) ! y : Big
      PinIdx      (iCellRayIdx, icRay)    = jhPin
    END DO
  END DO
  
  nTotCellRay(icRay) = iCellRayIdx
  nTotRaySeg (icRay) = irSegIdx
END DO
! ----------------------------------------------------
DO icRay = 1, nCoreRay
  DO iRaySeg = 1, nTotRaySeg(icRay)
    DO iPol = 1, nPolarAng
      ExpAppPolar(iPol, iRaySeg, icRay) = ExpA(ExpAppIdx(iRaySeg, icRay), iPol) * optlenlist(iRaySeg, icRay) + ExpB(ExpAppIdx(iRaySeg, icRay), iPol)
    END DO
  END DO
END DO
! ----------------------------------------------------
IF (krot .EQ. 1) THEN
  jbeg = 1; jend = nCoreRay; jinc = 1
ELSE
  jbeg = nCoreRay; jend = 1; jinc = -1
END IF

DO icRay = jbeg, jend, jinc
  jcRay = hRotRay_Loc%cRayIdx(icRay)
  iAzi  = hcRay(abs(hRotRay_Loc%cRayIdx(icRay)))%AzmIdx
  
  IF (krot .EQ. 2) jcRay = -jcRay !Reverse the Sweep Direction
  
  nRaySeg = nTotRaySeg(icRay)
  ! ----------------------------------------------------
  IF (jcRay .GT. 0) THEN
    PhiAngOut1g(:, 1) = locphiout(:)
    
    DO iRaySeg = 1, nRaySeg
      iReg = FsrIdx(iRaySeg, icRay)
      
      DO iPol = 1, nPolarAng
        wt   = wtang(iPol, iAzi)
        phid = (PhiAngOut1g(iPol, iRaySeg) - SrcAng1g1(iPol, iReg, iAzi)) * ExpAppPolar(iPol, iRaySeg, icRay)
        
        PhiAngOut1g(ipol, iRaySeg+1) = PhiAngOut1g(iPol, iRaySeg) - phid
        
        IF (lAFSS) THEN
          phia1g(FORWARD, iPol, iAzi, iReg) = phia1g(FORWARD, iPol, iAzi, iReg) + phid
        ELSE
          phis1g(iReg) = phis1g(iReg) + wt * phid
          
          DO iod = 1, nod
            phim1g(iod, iReg) = phim1g(iod, iReg) + mwt(iod, iPol, iAzi) * phid
          END DO
        END IF
      END DO
    END DO
    
    locphiout(:) = PhiAngOut1g(:, nRaySeg+1)
    
    ! Surface
    IF (ljout) THEN
      DO iRaySeg = 1, nTotCellRay(icRay)
        DO iPol = 1, nPolarAng
          wt    = wtang(iPol, iAzi)
          iCel  = PinIdx(iRaySeg, icRay)
          iSurf = SurfIdx(iRaySeg, icRay, 1)
          
          Jout1g(2, iSurf, iCel) = Jout1g(2, iSurf, iCel) + wt * PhiAngOut1g(iPol, CellRayIdxSt(iRaySeg, icRay, 1)+1)
          
          isurf = SurfIdx(iRaySeg, icRay, 2)
          
          Jout1g(1, iSurf, iCel) = Jout1g(1, iSurf, iCel) + wt * PhiAngOut1g(iPol, CellRayIdxSt(iRaySeg, icRay, 2))
        END DO
      END DO
    END IF
  ! ----------------------------------------------------
  ELSE
    PhiAngOut1g(:, nRaySeg+2) = locphiout(:)
    
    iRaySeg = nRaySeg + 1
    
    DO iRayseg1 = 1, nRaySeg
      iRaySeg = iRaySeg - 1
      iReg    = FsrIdx(iRaySeg, icRay)
      
      DO iPol = 1, nPolarAng
        wt   = wtang(iPol, iAzi)
        phid = (PhiAngOut1g(iPol, iRaySeg + 2) - SrcAng1g2(iPol, iReg, iAzi)) * ExpAppPolar(iPol, iRaySeg, icRay)
        
        PhiAngOut1g(iPol, iRaySeg+1) = PhiAngOut1g(iPol, iRaySeg + 2) - phid
        
        IF (lAFSS) THEN
          phia1g(BACKWARD, iPol, iAzi, iReg) = phia1g(BACKWARD, iPol, iAzi, iReg) + phid
        ELSE
          phis1g(iReg) = phis1g(iReg) + wt * phid
          
          DO iod = 1, nod
            phim1g(iod, iReg) = phim1g(iod, iReg) + mwt2(iod, iPol, iAzi) * phid
          END DO
        END IF
      END DO
    END DO
    
    locphiout(:) = PhiAngOut1g(:, 2)
    
    ! Surface 
    IF (lJout) THEN
      DO iRaySeg = 1, nTotCellRay(icRay)
        DO iPol = 1, nPolarAng
          wt    = wtang(iPol, iAzi)
          iCel  = PinIdx(iRaySeg, icRay)
          iSurf = SurfIdx(iRaySeg, icRay, 2)
          
          Jout1g(2, iSurf, iCel) = Jout1g(2, iSurf, iCel) + wt * PhiAngOut1g(iPol, CellRayIdxSt(iRaySeg, icRay, 2)+1)
          
          isurf = SurfIdx(iRaySeg, icRay, 1)
          
          Jout1g(1, iSurf, iCel) = Jout1g(1, iSurf, iCel) + wt * PhiAngOut1g(iPol, CellRayIdxSt(iRaySeg, icRay, 1)+2)
        END DO
      END DO
    END IF
  END IF
END DO

TrackingDat%PhiAngIn1g(:, PhiAngOutSvIdx) = locphiout
! ----------------------------------------------------
! Loc.
NULLIFY (wtang)
NULLIFY (phis1g)
NULLIFY (xst1g)
NULLIFY (ExpA)
NULLIFY (ExpB)
NULLIFY (Pin)
NULLIFY (FsrIdx)
NULLIFY (ExpAppIdx)
NULLIFY (OptLenList)
NULLIFY (jout1g)
NULLIFY (ExpAppPolar)
NULLIFY (PhiAngOut1g)

IF (lAFSS) NULLIFY (phia1g)

! P1
NULLIFY (phim1g)
NULLIFY (mwt)
NULLIFY (mwt2)
NULLIFY (SrcAng1g1)
NULLIFY (SrcAng1g2)

! Hex.
NULLIFY (hRotRay_Loc)
NULLIFY (CelRay_Loc)
NULLIFY (hRotRay_Loc)
! ----------------------------------------------------

END SUBROUTINE HexTrackRotRayP1_GM
! ------------------------------------------------------------------------------------------------------------