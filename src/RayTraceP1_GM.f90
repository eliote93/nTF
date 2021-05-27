#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceP1GM_OMP(RayInfo, CoreInfo, phis, phim, PhiAngIn, xst, src, srcm, jout, iz, ljout, ScatOd, FastMocLv)

USE TIMER
USE ALLOCS
USE OMP_LIB
USE PARAM,   ONLY : FALSE, ZERO, ONE
USE TYPEDEF, ONLY : RayInfo_Type, CoreInfo_type, Pin_Type, Cell_Type
USE Moc_Mod, ONLY : TrackingDat, Comp, SrcAng1, SrcAng2
USE PE_MOD,  ONLY : PE
USE CNTL,    ONLY : nTracerCntl

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:)     :: phis, xst, src
REAL, POINTER, DIMENSION(:,:)   :: PhiAngIn, srcm, phim
REAL, POINTER, DIMENSION(:,:,:) :: jout

INTEGER :: iz, ScatOd
LOGICAL :: ljout
INTEGER, OPTIONAL :: FastMocLv
! ----------------------------------------------------
TYPE (Cell_Type), POINTER, DIMENSION(:) :: Cell
TYPE (Pin_Type),  POINTER, DIMENSION(:) :: Pin

INTEGER :: nAziAng, nPolarAng, nFsr, nxy, nThread
INTEGER :: ithr, FsrIdxSt, icel, iazi, ipol, iod, iRotRay, ifsr, jfsr, ipin, krot
REAL :: wttmp, tmpsrc, ONETHREE, ONEFIVE, ONESEVEN
! ----------------------------------------------------

nAziAng   = RayInfo%nAziAngle
nPolarAng = RayInfo%nPolarAngle

nFsr = CoreInfo%nCoreFsr
nxy  = CoreInfo%nxy

nThread = PE%nThread
! ----------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ithr, iazi, ifsr, ipol, tmpsrc)
ithr = omp_get_thread_num() + 1

DO iazi = 1, nAziAng
  DO ifsr = PE%myOmpFsrBeg(ithr), PE%myOmpFsrEnd(ithr)
    DO ipol = 1, nPolarAng
      SrcAng1(ipol, ifsr, iazi) = src(ifsr)
      SrcAng2(ipol, ifsr, iazi) = src(ifsr)
      
      tmpsrc = comp(1, ipol, iazi) * srcm(1, ifsr) + comp(2, ipol, iazi) * srcm(2, ifsr)
      
      SrcAng1(ipol, ifsr, iazi) = SrcAng1(ipol, ifsr, iazi) + tmpsrc
      SrcAng2(ipol, ifsr, iazi) = SrcAng2(ipol, ifsr, iazi) - tmpsrc
      
      IF (scatod .LT. 2) CYCLE
      
      tmpsrc = comp(3, ipol, iazi) * srcm(3, ifsr) + comp(4, ipol, iazi) * srcm(4, ifsr) + comp(5, ipol, iazi) * srcm(5, ifsr)
      
      SrcAng1(ipol, ifsr, iazi) = SrcAng1(ipol, ifsr, iazi) +  tmpsrc
      SrcAng2(ipol, ifsr, iazi) = SrcAng2(ipol, ifsr, iazi) +  tmpsrc
      
      IF (scatod .LT. 3) CYCLE
      
      tmpsrc = comp(6, ipol, iazi) * srcm(6, ifsr) + comp(7, ipol, iazi) * srcm(7, ifsr) + comp(8, ipol, iazi) * srcm(8, ifsr) + comp(9, ipol, iazi) * srcm(9, ifsr)
      
      SrcAng1(ipol, ifsr, iazi) = SrcAng1(ipol, ifsr, iazi) + tmpsrc
      SrcAng2(ipol, ifsr, iazi) = SrcAng2(ipol, ifsr, iazi) - tmpsrc
     END DO
  END DO
END DO
!$OMP END PARALLEL
! ----------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ithr, krot, iRotRay)
ithr = omp_get_thread_num() + 1

TrackingDat(ithr)%phis = ZERO
TrackingDat(ithr)%phim = ZERO
IF (ljout) TrackingDat(ithr)%jout = ZERO

TrackingDat(ithr)%PhiAngIn => PhiAngIn
TrackingDat(ithr)%src      => src
TrackingDat(ithr)%xst      => xst
TrackingDat(ithr)%srcm     => srcm
TrackingDat(ithr)%SrcAng1  => SrcAng1
TrackingDat(ithr)%SrcAng2  => SrcAng2

DO krot = 1, 2
  IF (nTracerCntl%lHex) THEN
    !$OMP DO SCHEDULE(GUIDED)
    DO iRotRay = 1, RayInfo%nRotRay
      CALL HexTrackRotRayP1GM_OMP(RayInfo, CoreInfo, TrackingDat(ithr), lJout, iRotRay, iz, krot, ScatOd, FastMocLv)
    END DO
    !$OMP END DO NOWAIT
  ELSE
    !$OMP DO SCHEDULE(GUIDED)
    DO iRotRay = 1, RayInfo%nRotRay
      CALL RecTrackRotRayP1GM_OMP(RayInfo, CoreInfo, TrackingDat(ithr), lJout, iRotRay, iz, krot, ScatOd, FastMocLv)
    END DO
    !$OMP END DO NOWAIT
  END IF
END DO
!$OMP END PARALLEL
! ----------------------------------------------------
phis = ZERO
phim = ZERO
IF (ljout) jout = ZERO

DO ithr = 1, nThread
  phis = phis + TrackingDat(ithr)%phis
  phim = phim + TrackingDat(ithr)%phim
  
  IF (.NOT. ljout) CYCLE
  
  jout = jout + TrackingDat(ithr)%jout
END DO
! ----------------------------------------------------
ONETHREE = ONE / 3._8
ONEFIVE  = ONE / 5._8
ONESEVEN = ONE / 7.

Cell => CoreInfo%CellInfo
Pin  => CoreInfo%Pin

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ithr, ipin, FsrIdxSt, icel, ifsr, jfsr, wttmp)
ithr = omp_get_thread_num() + 1

DO ipin = PE%myOmpNxyBeg(ithr), PE%myOmpNxyEnd(ithr)
  FsrIdxSt = Pin(ipin)%FsrIdxSt
  icel     = Pin(ipin)%Cell(iz)
  
  DO ifsr = 1, Cell(icel)%nFsr
    jfsr  = FsrIdxSt + ifsr - 1
    wttmp = ONE / xst(jfsr) / Cell(icel)%vol(ifsr)
    
    phis(jfsr) = phis(jfsr) * wttmp + src(jfsr)
    
    phim(1:2, jfsr) = phim(1:2, jfsr) * wttmp + srcm(1:2, jfsr) * ONETHREE
    
    IF (scatod .LT. 2) CYCLE
    
    phim(3:5, jfsr) = phim(3:5, jfsr) * wttmp + srcm(3:5, jfsr) * ONEFIVE
    
    IF (scatod .LT. 3) CYCLE
    
    phim(6:9, jfsr) = phim(6:9, jfsr) * wttmp + srcm(6:9, jfsr) * ONESEVEN
  END DO
END DO
!$OMP END PARALLEL

NULLIFY (Cell)
NULLIFY (Pin)
! ----------------------------------------------------

END SUBROUTINE RayTraceP1GM_OMP

! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayP1GM_OMP(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, krot, ScatOd, FastMocLv)

USE ALLOCS
USE PARAM,   ONLY : TRUE, FALSE, ZERO
USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type, Pin_Type, Asy_Type, AsyInfo_Type, PinInfo_Type, Cell_Type, AziAngleInfo_Type, PolarAngle_Type, ModRayInfo_type, AsyRayInfo_type,  &
                    CoreRayInfo_Type, RotRayInfo_Type, CellRayInfo_type, FastCoreRayDat_Type, TrackingDat_Type, FastRaySegDat_Type
USE Moc_Mod, ONLY : nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay

IMPLICIT NONE

TYPE (RayInfo_Type),     INTENT(INOUT) :: RayInfo
TYPE (CoreInfo_Type),    INTENT(INOUT) :: CoreInfo
TYPE (TrackingDat_Type), INTENT(INOUT) :: TrackingDat

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz, krot, ScatOd, FastMocLv

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

REAL, POINTER, DIMENSION(:)     :: LenSeg, phis, src, xst
REAL, POINTER, DIMENSION(:,:)   :: OptLenList, PhiAngOutPolar, PhiAngIn, EXPA, EXPB, wtang, phim
REAL, POINTER, DIMENSION(:,:,:) :: ExpAppPolar, jout, wtsurf, mwt, mwt2, SrcAng1, SrcAng2

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

FsrIdx         => TrackingDat%FsrIdx
ExpAppIdx      => TrackingDat%ExpAppIdx
OptLenList     => TrackingDat%OptLenList
ExpAppPolar    => TrackingDat%ExpAppPolar
Phis           => TrackingDat%phis
src            => TrackingDat%src
xst            => TrackingDat%xst
jout           => TrackingDat%jout
PhiAngOutPolar => TrackingDat%PhiAngOutPolar
PhiAngIn       => TrackingDat%phiAngIn
EXPA           => TrackingDat%EXPA
EXPB           => TrackingDat%EXPB
wtang          => TrackingDat%wtang
wtsurf         => TrackingDat%wtsurf
mwt            => TrackingDat%mwt
mwt2           => TrackingDat%mwt2
SrcAng1        => TrackingDat%SrcAng1
SrcAng2        => TrackingDat%SrcAng2
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
          
          tau = -LenSeg(iRaySeg) * xst(ireg)   !
          tau = -CellRay%LenSeg(iRaySeg) * xst(ireg)   !
          
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
        
        tau = -CellRay1D%LenSeg(k) * xst(ireg)
        
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
      
      tau = -FastRaySeg%LenSeg(l) * xst(ireg)
      
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
      ExpAppPolar(ipol, l, j) = expa(ExpAppIdx(l, j), ipol) * optlenlist(l, j) + expb(ExpAppIdx(l, j), ipol)
    END DO
  END DO
END DO
! ----------------------------------------------------
ALLOCATE (phiobdPolar(nPolarAng))

PhiAnginSvIdx  = RayInfo%PhiAngInSvIdx(iRotRay ,krot)
PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay, krot)
phiobdPolar(:) = PhiAngIn(:,PhiAnginSvIdx)

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
    PhiAngOutPolar(:, 1) = phiobdPolar(:)
    
    DO ir = 1, nRaySeg
      ireg = FsrIdx(ir, j)
      
      DO ipol = 1, nPolarAng
        wt = wtang(ipol, iazi)
        
        phid = (PhiAngOutPolar(ipol, ir) - SrcAng1(ipol, ireg, iazi)) * ExpAppPolar(ipol, ir, j)
        
        PhiAngOutPolar(ipol, ir+1) = PhiAngOutPolar(ipol, ir) - phid
        
        phis(ireg) = phis(ireg) + wt*phid
        
        DO iod = 1, od
          phim(iod, ireg) = phim(iod, ireg) + mwt(iod, ipol, iazi) * phid
        END DO
      END DO
    END DO
    
    phiobdPolar(:) = PhiAngOutPolar(:, nRaySeg+1)
    
    ! Surface
    IF (ljout) THEN
      DO ir = 1, nTotCellRay(j)
        DO ipol = 1, nPolarANg
          wt = wtang(ipol, iazi)
          
          wt2(1:4) = wtsurf(ipol, iazi, 1:4)
          
          icel = PinIdx(ir, j); isurf = SurfIdx(ir, j, 1)
          
          Jout(2, isurf, icel) = Jout(2, isurf, icel) + PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 1)+1) * wt
          Jout(3, isurf, icel) = Jout(3, isurf, icel) + PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 1)+1) * wt2(isurf)
          
          isurf = SurfIdx(ir, j, 2)
          
          Jout(1, isurf, icel) = Jout(1, isurf, icel) + PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 2)) * wt
          Jout(3, isurf, icel) = Jout(3, isurf, icel) + PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 2)) * wt2(isurf)
        END DO
      END DO
    END IF
  ELSE
    PhiAngOutPolar(:, nRaySeg+2) = phiobdPolar(:)
    
    ir = nRaySeg + 1
    
    DO ir1 = 1, nRaySEg
      ir   = ir - 1
      ireg = FsrIdx(ir, j)
      
      DO ipol = 1, nPolarAng
        wt = wtang(ipol, iazi)
        
        phid = (PhiAngOutPolar(ipol, ir + 2) - SrcAng2(ipol, ireg, iazi)) * ExpAppPolar(ipol, ir, j)
        
        PhiAngOutPolar(ipol, ir+1) = PhiAngOutPolar(ipol, ir + 2) - phid
        
        phis(ireg) = phis(ireg) + wt * phid
        
        DO iod = 1, od
          phim(iod, ireg) = phim(iod, ireg) + mwt2(iod, ipol, iazi) * phid
        END DO
      END DO
    END DO
    
    phiobdPolar(:) = PhiAngOutPolar(:, 2)
    
    IF (lJout) THEN
      DO ir = 1, nTotCellRay(j)
        DO ipol = 1, nPolarAng
          wt = wtang(ipol, iazi)
          
          wt2(1:4) = wtsurf(ipol, iazi, 1:4)
          
          icel = PinIdx(ir, j); isurf = SurfIdx(ir, j, 2)
          
          Jout(2, isurf, icel) = Jout(2, isurf, icel) + PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 2)+1) * wt
          Jout(3, isurf, icel) = Jout(3, isurf, icel) + PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 2)+1) * wt2(isurf)
          
          isurf = SurfIdx(ir, j, 1)
          
          Jout(1, isurf, icel) = Jout(1, isurf, icel) + PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 1)+2) * wt
          Jout(3, isurf, icel) = Jout(3, isurf, icel) + PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 1)+2) * wt2(isurf)
        END DO
      END DO
    END IF
  END IF
END DO

PhiAngIn(:, PhiAngOutSvIdx) = phiobdPolar(:)
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
NULLIFY (Phis)
NULLIFY (src)
NULLIFY (xst)
NULLIFY (jout)
NULLIFY (PhiAngOutPolar)
NULLIFY (PhiAngIn)
NULLIFY (EXPA)
NULLIFY (EXPB)
NULLIFY (wtang)
NULLIFY (wtsurf)
NULLIFY (mwt)
NULLIFY (mwt2)
NULLIFY (SrcAng1)
NULLIFY (SrcAng2)
! ----------------------------------------------------

END SUBROUTINE RecTrackRotRayP1GM_OMP
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayP1GM_OMP(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, krot, ScatOd, FastMocLv)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, TrackingDat_Type, Pin_Type
USE Moc_Mod, ONLY : nMaxCellRay, nMaxCoreRay
USE HexData, ONLY : hAsy
USE HexType, ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexCoreRay, Type_HexRotRay
USE HexData, ONLY : haRay, hcRay, hRotRay, hAsyTypInfo

IMPLICIT NONE

TYPE (RayInfo_Type),     INTENT(INOUT) :: RayInfo
TYPE (CoreInfo_Type),    INTENT(INOUT) :: CoreInfo
TYPE (TrackingDat_Type), INTENT(INOUT) :: TrackingDat

LOGICAL, INTENT(IN) :: ljout
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

REAL, POINTER, DIMENSION(:)     :: phis, src, xst
REAL, POINTER, DIMENSION(:,:)   :: EXPA, EXPB, wtang, phim, OptLenList, PhiAngOutPolar
REAL, POINTER, DIMENSION(:,:,:) :: jout, mwt, mwt2, SrcAng1, SrcAng2, ExpAppPolar

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

phis           => TrackingDat%phis
phim           => TrackingDat%phim
src            => TrackingDat%src
xst            => TrackingDat%xst
jout           => TrackingDat%jout
EXPA           => TrackingDat%EXPA
EXPB           => TrackingDat%EXPB
wtang          => TrackingDat%wtang
mwt            => TrackingDat%mwt
mwt2           => TrackingDat%mwt2
SrcAng1        => TrackingDat%SrcAng1
SrcAng2        => TrackingDat%SrcAng2
locphiout       = TrackingDat%PhiAngIn(:, PhiAnginSvIdx)
FsrIdx         => TrackingDat%FsrIdx
ExpAppIdx      => TrackingDat%ExpAppIdx
OptLenList     => TrackingDat%OptLenList
ExpAppPolar    => TrackingDat%ExpAppPolar
PhiAngOutPolar => TrackingDat%PhiAngOutPolar

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
        tau      = -CelRay_Loc%SegLgh(iRaySeg) * XsT(iReg) ! Optimum Length
        
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
      ExpAppPolar(iPol, iRaySeg, icRay) = expa(ExpAppIdx(iRaySeg, icRay), iPol) * optlenlist(iRaySeg, icRay) &
                                        + expb(ExpAppIdx(iRaySeg, icRay), iPol)
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
  IF (jcRay > 0) THEN
    PhiAngOutPolar(:, 1) = locphiout(:)
    
    DO iRaySeg = 1, nRaySeg
      iReg = FsrIdx(iRaySeg, icRay)
      
      DO iPol = 1, nPolarAng
        wt   = wtang(iPol, iAzi)
        phid = (PhiAngOutPolar(iPol, iRaySeg) - SrcAng1(iPol, iReg, iAzi)) * ExpAppPolar(iPol, iRaySeg, icRay)
        
        PhiAngOutPolar(ipol, iRaySeg+1) = PhiAngOutPolar(iPol, iRaySeg) - phid
        
        phis(iReg) = phis(iReg) + wt * phid
        
        DO iod = 1, nod
          phim(iod, iReg) = phim(iod, iReg) + mwt(iod, iPol, iAzi) * phid
        END DO
      END DO
    END DO
    
    locphiout(:) = PhiAngOutPolar(:, nRaySeg+1)
    
    ! Surface
    IF (ljout) THEN
      DO iRaySeg = 1, nTotCellRay(icRay)
        DO iPol = 1, nPolarAng
          wt    = wtang(iPol, iAzi)
          iCel  = PinIdx(iRaySeg, icRay)
          iSurf = SurfIdx(iRaySeg, icRay, 1)
          
          Jout(2, iSurf, iCel) = Jout(2, iSurf, iCel) + wt * PhiAngOutPolar(iPol, CellRayIdxSt(iRaySeg, icRay, 1)+1)
          !Jout(3, iSurf, iCel) = Jout(3, iSurf, iCel) + wt2(isurf) * PhiAngOut(CellRayIdxSt(iRaySeg, icRay, 1)+1)
          
          isurf = SurfIdx(iRaySeg, icRay, 2)
          
          Jout(1, iSurf, iCel) = Jout(1, iSurf, iCel) + wt * PhiAngOutPolar(iPol, CellRayIdxSt(iRaySeg, icRay, 2))
          !Jout(3, iSurf, iCel) = Jout(3, iSurf, iCel) + wt2(iSurf) * PhiAngOut(CellRayIdxSt(iRaySeg, icRay, 2))
        END DO
      END DO
    END IF
  ! ----------------------------------------------------
  ELSE
    PhiAngOutPolar(:, nRaySeg+2) = locphiout(:)
    
    iRaySeg = nRaySeg + 1
    
    DO iRayseg1 = 1, nRaySeg
      iRaySeg = iRaySeg - 1
      iReg    = FsrIdx(iRaySeg, icRay)
      
      DO iPol = 1, nPolarAng
        wt   = wtang(iPol, iAzi)
        phid = (PhiAngOutPolar(iPol, iRaySeg + 2) - SrcAng2(iPol, iReg, iAzi)) * ExpAppPolar(iPol, iRaySeg, icRay)
        
        PhiAngOutPolar(iPol, iRaySeg+1) = PhiAngOutPolar(iPol, iRaySeg + 2) - phid
        
        phis(iReg) = phis(iReg) + wt * phid
        
        DO iod = 1, nod
          phim(iod, iReg) = phim(iod, iReg) + mwt2(iod, iPol, iAzi) * phid
        END DO
      END DO
    END DO
    
    locphiout(:) = PhiAngOutPolar(:, 2)
    
    ! Surface 
    IF (lJout) THEN
      DO iRaySeg = 1, nTotCellRay(icRay)
        DO iPol = 1, nPolarAng
          wt    = wtang(iPol, iAzi)
          iCel  = PinIdx(iRaySeg, icRay)
          iSurf = SurfIdx(iRaySeg, icRay, 2)
          
          Jout(2, iSurf, iCel) = Jout(2, iSurf, iCel) + wt * PhiAngOutPolar(iPol, CellRayIdxSt(iRaySeg, icRay, 2)+1)
          !Jout(3, iSurf, iCel) = Jout(3, iSurf, iCel) + wt2(iSurf) * PhiAngOut(CellRayIdxSt(iRaySeg, icRay, 2)+1)
          
          isurf = SurfIdx(iRaySeg, icRay, 1)
          
          Jout(1, iSurf, iCel) = Jout(1, iSurf, iCel) + wt * PhiAngOutPolar(iPol, CellRayIdxSt(iRaySeg, icRay, 1)+2)
          !Jout(3, iSurf, iCel) = Jout(3, iSurf, iCel) + wt2(iSurf) * PhiAngOut(CellRayIdxSt(iRaySeg, icRay, 1) + 2)
        END DO
      END DO
    END IF
  END IF
END DO

TrackingDat%PhiAngIn(:, PhiAngOutSvIdx) = locphiout
! ----------------------------------------------------
! Loc.
NULLIFY (wtang)
NULLIFY (phis)
NULLIFY (src)
NULLIFY (xst)
NULLIFY (EXPA)
NULLIFY (EXPB)
NULLIFY (Pin)
NULLIFY (FsrIdx)
NULLIFY (ExpAppIdx)
NULLIFY (OptLenList)
NULLIFY (jout)
NULLIFY (ExpAppPolar)
NULLIFY (PhiAngOutPolar)

! P1
NULLIFY (phim)
NULLIFY (mwt)
NULLIFY (mwt2)
NULLIFY (SrcAng1)
NULLIFY (SrcAng2)

! Hex.
NULLIFY (hRotRay_Loc)
NULLIFY (CelRay_Loc)
NULLIFY (hRotRay_Loc)
! ----------------------------------------------------

END SUBROUTINE HexTrackRotRayP1GM_OMP
! ------------------------------------------------------------------------------------------------------------