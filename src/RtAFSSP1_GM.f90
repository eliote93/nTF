SUBROUTINE RayTraceP1GM_AFSS(RayInfo, CoreInfo, phis, phim, PhiAngIn, xst, src, srcm, jout, iz, ljout, ScatOd, FastMocLv)

USE OMP_LIB
USE PARAM,   ONLY : TRUE, FALSE, ZERO, ONE
USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type, Pin_Type, Cell_Type
USE Moc_Mod, ONLY : TrackingDat, wtang, comp, mwt, mwt2, SrcAng1, SrcAng2, nOmpAng, OmpRayBeg, OmpRayEnd, OmpRayBegBd, OmpRayEndBd, OmpMap
USE cntl,    ONLY : nTracerCntl
USE IOUTIL,  ONLY : terminate
USE PE_MOD,  ONLY : PE

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:)     :: phis, xst, src
REAL, POINTER, DIMENSION(:,:)   :: PhiAngIn, phim, srcm
REAL, POINTER, DIMENSION(:,:,:) :: jout

INTEGER :: iz, ScatOd
LOGICAL :: ljout
INTEGER, OPTIONAL :: FastMocLv
! ----------------------------------------------------
TYPE (Cell_Type), POINTER, DIMENSION(:) :: Cell
TYPE (Pin_Type),  POINTER, DIMENSION(:) :: Pin

INTEGER :: nAziAng, nPolarAng, nxy, nthr, iRotRay, ithr, FsrIdxSt, icel, ireg, iazi, ipol, OmpAng, ixy, iDir, iFSR, jFSR
REAL :: wttmp, tempsrc, ONETHREE, ONEFIVE, ONESEVEN
! ----------------------------------------------------

nAziAng   = RayInfo%nAziAngle
nPolarAng = RayInfo%nPolarAngle

nxy  = CoreInfo%nxy
nthr = PE%nthread
! ----------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ithr, iazi, iFSR, ipol, tempsrc)
ithr = 1
!$ ithr = omp_get_thread_num()+1

IF (ScatOd .EQ. 1) THEN
  DO iazi = 1, nAziAng
    DO iFSR = PE%myOmpFsrBeg(ithr), PE%myOmpFsrEnd(ithr)
      DO ipol = 1, nPolarAng
        SrcAng1(ipol, iFSR, iazi) = src(iFSR)
        SrcAng2(ipol, iFSR, iazi) = src(iFSR)
        
        tempsrc = comp(1, ipol, iazi) * srcm(1, iFSR) + comp(2, ipol, iazi) * srcm(2, iFSR)
        
        SrcAng1(ipol, iFSR, iazi) = SrcAng1(ipol, iFSR, iazi) + tempsrc
        SrcAng2(ipol, iFSR, iazi) = SrcAng2(ipol, iFSR, iazi) - tempsrc
       END DO
    END DO
  END DO
ELSE IF (ScatOd .EQ. 2) THEN
  DO iazi = 1, nAziAng
    DO iFSR = PE%myOmpFsrBeg(ithr), PE%myOmpFsrEnd(ithr)
      DO ipol = 1, nPolarAng
        SrcAng1(ipol, iFSR, iazi) = src(iFSR)
        SrcAng2(ipol, iFSR, iazi) = src(iFSR)
        
        tempsrc = comp(1, ipol, iazi) * srcm(1, iFSR) + comp(2, ipol, iazi) * srcm(2, iFSR)
        
        SrcAng1(ipol, iFSR, iazi) = SrcAng1(ipol, iFSR, iazi) + tempsrc
        SrcAng2(ipol, iFSR, iazi) = SrcAng2(ipol, iFSR, iazi) - tempsrc
        
        tempsrc = comp(3, ipol, iazi) * srcm(3, iFSR) + comp(4, ipol, iazi) * srcm(4, iFSR) + comp(5, ipol, iazi) * srcm(5, iFSR)
        
        SrcAng1(ipol, iFSR, iazi) = SrcAng1(ipol, iFSR, iazi) +  tempsrc
        SrcAng2(ipol, iFSR, iazi) = SrcAng2(ipol, iFSR, iazi) +  tempsrc
      END DO
    END DO
  END DO
ELSE IF (ScatOd .EQ. 3) THEN
  DO iazi = 1, nAziAng
    DO iFSR = PE%myOmpFsrBeg(ithr), PE%myOmpFsrEnd(ithr)
      DO ipol = 1, nPolarAng
        SrcAng1(ipol, iFSR, iazi) = src(iFSR)
        SrcAng2(ipol, iFSR, iazi) = src(iFSR)
        
        tempsrc = comp(1, ipol, iazi) * srcm(1, iFSR) + comp(2, ipol, iazi) * srcm(2, iFSR) + comp(6, ipol, iazi) * srcm(6, iFSR) + comp(7, ipol, iazi) * srcm(7, iFSR) + comp(8, ipol, iazi) * srcm(8, iFSR) + comp(9, ipol, iazi) * srcm(9, iFSR)
        
        SrcAng1(ipol, iFSR, iazi) = SrcAng1(ipol, iFSR, iazi) + tempsrc
        SrcAng2(ipol, iFSR, iazi) = SrcAng2(ipol, iFSR, iazi) - tempsrc
        
        tempsrc = comp(3, ipol, iazi) * srcm(3, iFSR) + comp(4, ipol, iazi) * srcm(4, iFSR) + comp(5, ipol, iazi) * srcm(5, iFSR)
        
        SrcAng1(ipol, iFSR, iazi) = SrcAng1(ipol, iFSR, iazi) +  tempsrc
        SrcAng2(ipol, iFSR, iazi) = SrcAng2(ipol, iFSR, iazi) +  tempsrc
      END DO
    END DO
  END DO
END IF
!$OMP END PARALLEL
! ----------------------------------------------------
!$ call omp_set_dynamic(FALSE)
!$ call omp_set_num_threads(nthr)

DO ithr = 1, nthr
  TrackingDat(ithr)%PhiAngIn => PhiAngIn
  TrackingDat(ithr)%src      => src
  TrackingDat(ithr)%xst      => xst
  TrackingDat(ithr)%SrcAng1  => SrcAng1
  TrackingDat(ithr)%SrcAng2  => SrcAng2
  
  DO iazi = 1, nOmpAng
    TrackingDat(ithr)%phi1a(:,:,iazi) = ZERO
    TrackingDat(ithr)%phi2a(:,:,iazi) = ZERO
  END DO
  
  DO ixy = 1, nxy
    TrackingDat(ithr)%jout(:, :, ixy) = zero
  END DO
END DO

phis = ZERO
phim = ZERO
IF (ljout) Jout = ZERO
! ----------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ithr, iRotRay, iDir, iazi, ipol, OmpAng, ixy)
ithr = 1
!$  ithr = omp_get_thread_num()+1
!$OMP BARRIER
IF (nTracerCntl%lScatBd) THEN
  DO iRotRay = OmpRayBeg(ithr), OmpRayEnd(ithr)
    CALL RecTrackRotRayP1GM_AFSS(RayInfo, CoreInfo, TrackingDat(ithr), lJout, iRotRay, iz, FastMocLv)
  END DO
ELSE
  DO iDir = 1, 2
    DO iRotRay = OmpRayBegBd(iDir, ithr), OmpRayEndBd(iDir, ithr)
      CALL RecTrackRotRayP1GM_AFSS(RayInfo, CoreInfo, TrackingDat(ithr), lJout, iRotRay, iz, FastMocLv)
    END DO
  END DO
END IF
!$OMP BARRIER
IF (ScatOd .EQ. 1) THEN
  DO ithr = 1, nthr
    DO iazi = 1, nOmpAng
      DO iFSR = PE%myOmpFsrBeg(ithr), PE%myOmpFsrEnd(ithr)
        DO ipol = 1, nPolarAng
          OmpAng = OmpMap(ithr, iazi)
          
          phis(iFSR) = phis(iFSR) + wtang(ipol, OmpAng) * (TrackingDat(ithr)%phi1a(ipol, iFSR, iazi) + TrackingDat(ithr)%phi2a(ipol, iFSR, iazi))
          
          phim(1:2, iFSR) = phim(1:2, iFSR) + mwt(1:2, ipol, OmpAng) * (TrackingDat(ithr)%phi1a(ipol, iFSR, iazi) - TrackingDat(ithr)%phi2a(ipol, iFSR, iazi))
        END DO
      END DO
    END DO
  END DO
ELSE IF (ScatOd .EQ. 2) THEN
  DO ithr = 1, nthr
    DO iazi = 1, nOmpAng
      DO iFSR = PE%myOmpFsrBeg(ithr), PE%myOmpFsrEnd(ithr)
        DO ipol = 1, nPolarAng
          OmpAng = OmpMap(ithr, iazi)
          
          phis(iFSR) = phis(iFSR) + wtang(ipol, OmpAng) * (TrackingDat(ithr)%phi1a(ipol, iFSR, iazi) + TrackingDat(ithr)%phi2a(ipol, iFSR, iazi))
          
          phim(1:2, iFSR) = phim(1:2, iFSR) + mwt(1:2, ipol, OmpAng) * (TrackingDat(ithr)%phi1a(ipol, iFSR, iazi) - TrackingDat(ithr)%phi2a(ipol, iFSR, iazi))
          phim(3:5, iFSR) = phim(3:5, iFSR) + mwt(3:5, ipol, OmpAng) * (TrackingDat(ithr)%phi1a(ipol, iFSR, iazi) + TrackingDat(ithr)%phi2a(ipol, iFSR, iazi))
        END DO
      END DO
    END DO
  END DO
ELSE IF (ScatOd .EQ. 3) THEN
  DO ithr = 1, nthr
    DO iazi = 1, nOmpAng
      DO iFSR = PE%myOmpFsrBeg(ithr), PE%myOmpFsrEnd(ithr)
        DO ipol = 1, nPolarAng
          OmpAng = OmpMap(ithr, iazi)
          
          phis(iFSR) = phis(iFSR) + wtang(ipol, OmpAng) * (TrackingDat(ithr)%phi1a(ipol, iFSR, iazi) + TrackingDat(ithr)%phi2a(ipol, iFSR, iazi))
          
          phim(1:2, iFSR) = phim(1:2, iFSR) + mwt(1:2, ipol, OmpAng) * (TrackingDat(ithr)%phi1a(ipol, iFSR, iazi) - TrackingDat(ithr)%phi2a(ipol, iFSR, iazi))
          phim(3:5, iFSR) = phim(3:5, iFSR) + mwt(3:5, ipol, OmpAng) * (TrackingDat(ithr)%phi1a(ipol, iFSR, iazi) + TrackingDat(ithr)%phi2a(ipol, iFSR, iazi))
          phim(6:9, iFSR) = phim(6:9, iFSR) + mwt(6:9, ipol, OmpAng) * (TrackingDat(ithr)%phi1a(ipol, iFSR, iazi) - TrackingDat(ithr)%phi2a(ipol, iFSR, iazi))
        END DO
      END DO
    END DO
  END DO
END IF

IF (ljout) THEN
  DO ithr = 1, nthr
    DO ixy = PE%myOmpNxyBeg(ithr), PE%myOmpNxyEnd(ithr)
      jout(:, :, ixy) = jout(:, :, ixy) + TrackingDat(ithr)%jout(:, :, ixy)
    END DO
  END DO
END IF
!$OMP END PARALLEL
! ----------------------------------------------------
ONETHREE = ONE / 3._8
ONEFIVE  = ONE / 5._8
ONESEVEN = ONE / 7._8

Cell => CoreInfo%CellInfo
Pin  => CoreInfo%Pin

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ithr, ixy, FsrIdxSt, icel, iFSR, jFSR, wttmp)
ithr = 1
!$ ithr = omp_get_thread_num()+1

IF (ScatOd .EQ. 1) THEN
  DO ixy = PE%myOmpNxyBeg(ithr), PE%myOmpNxyEnd(ithr)
    FsrIdxSt = Pin(ixy)%FsrIdxSt
    icel     = Pin(ixy)%Cell(iz)
    
    DO iFSR = 1, Cell(icel)%nFsr
      jFSR = FsrIdxSt + iFSR - 1
      
      wttmp = 1._8 / (xst(jFSR) * Cell(icel)%vol(iFSR))
      
      phis(jFSR) = phis(jFSR) * wttmp + src(jFSR)
      
      phim(1:2, jFSR) = phim(1:2, jFSR) * wttmp + srcm(1:2, jFSR) * ONETHREE
    END DO
  END DO
ELSE IF (ScatOd .EQ. 2) THEN
  DO ixy = PE%myOmpNxyBeg(ithr), PE%myOmpNxyEnd(ithr)
    FsrIdxSt = Pin(ixy)%FsrIdxSt
    icel     = Pin(ixy)%Cell(iz)
    
    DO iFSR = 1, Cell(icel)%nFsr
      jFSR = FsrIdxSt + iFSR - 1
      
      wttmp = 1._8 / (xst(jFSR) * Cell(icel)%vol(iFSR))
      
      phis(jFSR) = phis(jFSR) * wttmp + src(jFSR)
      
      phim(1:2, jFSR) = phim(1:2, jFSR) * wttmp + srcm(1:2, jFSR) * ONETHREE
      phim(3:5, jFSR) = phim(3:5, jFSR) * wttmp + srcm(3:5, jFSR) * ONEFIVE
    END DO
  END DO
ELSE IF (ScatOd .EQ. 3) THEN
  DO ixy = PE%myOmpNxyBeg(ithr), PE%myOmpNxyEnd(ithr)
    FsrIdxSt = Pin(ixy)%FsrIdxSt
    icel     = Pin(ixy)%Cell(iz)
    
    DO iFSR = 1, Cell(icel)%nFsr
      jFSR = FsrIdxSt + iFSR - 1
      
      wttmp = 1._8 / (xst(jFSR) * Cell(icel)%vol(iFSR))
      
      phis(jFSR) = phis(jFSR) * wttmp + src(jFSR)
      
      phim(1:2, jFSR) = phim(1:2, jFSR) * wttmp + srcm(1:2, jFSR) * ONETHREE
      phim(3:5, jFSR) = phim(3:5, jFSR) * wttmp + srcm(3:5, jFSR) * ONEFIVE
      phim(6:9, jFSR) = phim(6:9, jFSR) * wttmp + srcm(6:9, jFSR) * ONESEVEN
    END DO
  END DO
END IF
!$OMP END PARALLEL

NULLIFY (Cell)
NULLIFY (Pin)
! ----------------------------------------------------

END SUBROUTINE RayTraceP1GM_AFSS
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayP1GM_AFSS(RayInfo, CoreInfo, TrackingDat, ljout, iRotRay, iz, FastMocLv)

USE PARAM,   ONLY : TRUE, FALSE
USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type, Pin_Type, Asy_Type, Cell_Type, AsyRayInfo_type,  &
                    CoreRayInfo_Type, RotRayInfo_Type, CellRayInfo_type, FastCoreRayDat_Type, TrackingDat_Type, FastRaySegDat_Type
USE Moc_Mod, ONLY : nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay

IMPLICIT NONE

TYPE (RayInfo_Type),     INTENT(INOUT) :: RayInfo
TYPE (CoreInfo_Type),    INTENT(INOUT) :: CoreInfo
TYPE (TrackingDat_Type), INTENT(INOUT) :: TrackingDat

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: iRotRay, iz, FastMocLv
! ----------------------------------------------------
TYPE (Pin_Type),         POINTER, DIMENSION(:) :: Pin
TYPE (Asy_Type),         POINTER, DIMENSION(:) :: Asy
TYPE (Cell_Type),        POINTER, DIMENSION(:) :: Cell
TYPE (AsyRayInfo_type),  POINTER, DIMENSION(:) :: AsyRay
TYPE (CoreRayInfo_Type), POINTER, DIMENSION(:) :: CoreRay
TYPE (RotRayInfo_Type),  POINTER, DIMENSION(:) :: RotRay

TYPE (FastCoreRayDat_Type), POINTER :: FastRay
TYPE (CellRayInfo_Type),    POINTER :: CellRay
TYPE (CellRayInfo_Type),    POINTER :: CellRay1D
TYPE (FastRaySegDat_Type),  POINTER :: FastRaySeg

INTEGER, POINTER, DIMENSION(:)   :: LocalFsrIdx
INTEGER, POINTER, DIMENSION(:,:) :: FsrIdx, ExpAppIdx

REAL, POINTER, DIMENSION(:)     :: LenSeg, phis, src, xst
REAL, POINTER, DIMENSION(:,:)   :: OptLenList, PhiAngOutPolar, PhiAngIn, ExpA, ExpB, wtang
REAL, POINTER, DIMENSION(:,:,:) :: ExpAppPolar, jout, wtsurf, phi1a, phi2a, SrcAng1, SrcAng2

INTEGER :: mp(2)
INTEGER :: iazi, ipol, iCoreRay, jCoreRay, iAsyRay, jAsyRay, iceray, irayseg, irot, idir, nCoreRay, nAsyRay, nPinRay, nRaySeg, FsrIdxSt, nPolarAng, nAziAng, iFSR
INTEGER :: ipin, icel, ibcel, iasy, ireg, isurf, irsegidx, icellrayidx, PhiAnginSvIdx, PhiAngOutSvIdx, jbeg, jend, jinc, ir, ir1, iOmpAzi, iPinRay

INTEGER, DIMENSION(nMaxCoreRay)                 :: nTotRaySeg, nTotCellRay
INTEGER, DIMENSION(nMaxCellRay, nMaxCoreRay)    :: PinIdx
INTEGER, DIMENSION(nMaxCellRay, nMaxCoreRay, 2) :: CellRayIdxSt, SurfIdx

REAL :: tau, phiobd, phid, wt
REAL, ALLOCATABLE :: phiobdPolar(:)

LOGICAL :: lFast

DATA mp /2, 1/
! ----------------------------------------------------

lFast = FALSE
IF (FastMocLv .GT. 0) lFast = TRUE

AsyRay   => RayInfo%AsyRay
CoreRay  => RayInfo%CoreRay
RotRay   => RayInfo%RotRay
nAziAng   = RayInfo%nAziAngle
nPolarAng = RayInfo%nPolarAngle

Asy  => CoreInfo%Asy
Pin  => CoreInfo%Pin
Cell => CoreInfo%CellInfo

FsrIdx         => TrackingDat%FsrIdx
ExpAppIdx      => TrackingDat%ExpAppIdx
OptLenList     => TrackingDat%OptLenList
ExpAppPolar    => TrackingDat%ExpAppPolar
Phis           => TrackingDat%phis
src            => TrackingDat%src
phi1a          => TrackingDat%phi1a
phi2a          => TrackingDat%phi2a
xst            => TrackingDat%xst
jout           => TrackingDat%jout
PhiAngOutPolar => TrackingDat%PhiAngOutPolar
PhiAngIn       => TrackingDat%phiAngIn
ExpA           => TrackingDat%ExpA
ExpB           => TrackingDat%ExpB
Wtang          => TrackingDat%wtang
wtsurf         => TrackingDat%wtsurf
SrcAng1        => TrackingDat%SrcAng1
SrcAng2        => TrackingDat%SrcAng2
! ----------------------------------------------------
nCoreRay = RotRay(irotRay)%nRay

IF (.NOT. lFast) THEN
  DO iCoreRay = 1, nCoreRay
    irsegidx    = 0
    icellrayidx = 0
    
    jCoreRay = RotRay(iRotRay)%RayIdx(iCoreRay)
    nAsyRay  = CoreRay(jCoreRay)%nRay
    
    DO iAsyRay = 1, nAsyRay
      jAsyRay = CoreRay(jCoreRay)%AsyRayIdx(iAsyRay)
      iasy    = CoreRay(jCoreRay)%AsyIdx   (iAsyRay)
      
      IF (iasy .EQ. 0) CYCLE ! Skip Dummy Assembly
      
      nPinRay = AsyRay(jAsyRay)%nCellRay
      
      DO iPinRay = 1, nPinRay
        ipin   = AsyRay(jAsyRay)%PinIdx   (iPinRay) ! Local Pin Idx(within Assembly)
        iceray = AsyRay(jAsyRay)%PinRayIdx(iPinRay)
        
        ipin = Asy(iAsy)%GlobalPinIdx(ipin) ! Global Pin Index
        
        icel     = Pin(ipin)%Cell(iz)
        FsrIdxSt = Pin(ipin)%FsrIdxSt
        
        ibcel    = Cell(icel)%BaseCellStr
        CellRay => Cell(ibcel)%CellRay(iceray)
        
        icellrayidx = icellrayidx + 1
        
        PinIdx      (icellrayidx, iCoreRay)    = ipin
        CellRayIdxSt(icellrayidx, iCoreRay, 2) = irsegidx + 1
        
        nRaySeg      = CellRay%nSeg
        LocalFsrIdx => CellRay%LocalFsrIdx
        LenSeg      => CellRay%LenSeg
        
        DO iRaySeg = 1, nRaySeg
          ireg = FsrIdxSt + CellRay%LocalFsrIdx(iRaySeg) - 1
          
          tau = -LenSeg(iRaySeg) * xst(ireg)   !
          tau = -CellRay%LenSeg(iRaySeg) * xst(ireg)   !
          
          irsegidx = irsegidx + 1
          
          FsrIdx    (irsegidx, iCoreRay) = ireg
          OptLenList(irsegidx, iCoreRay) = tau
          ExpAppIdx (irsegidx, iCoreRay) = max(INT(tau), -40000)
          ExpAppIdx (irsegidx, iCoreRay) = min(0, ExpAppIdx(irsegidx, iCoreRay))
        END DO
        
        CellRayIdxSt(icellrayidx, iCoreRay, 1) = irsegidx
        SurfIdx     (icellRayIdx, iCoreRay, 1) = AsyRay(jAsyRay)%PinRaySurf(2, iPinRay) ! OutSurface
        SurfIdx     (icellRayIdx, iCoreRay, 2) = AsyRay(jAsyRay)%PinRaySurf(1, iPinRay) ! Insurface
      END DO
    END DO
    
    nTotRaySeg (iCoreRay) = irsegidx
    nTotCellRay(iCoreRay) = icellRayIdx
  END DO
ELSE IF (FastMocLV .EQ. 1) THEN
  FastRay   => RayInfo%FastCoreRayDat(iRotRay, iz)
  CellRay1D => RayInfo%CellRay1D
  
  DO iCoreRay = 1, nCoreRay
    irsegidx = 0
    
    nTotRaySeg (iCoreRay) = FastRay%nTotRaySeg (iCoreRay)
    nTotCellRay(iCoreRay) = FastRay%nTotCellRay(iCoreRay)
    
    DO iPinRay = 1, FastRay%nTotCellRay(iCoreRay)
      PinIdx      (iPinRay, iCoreRay)    = FastRay%PinIdx      (iPinRay, iCoreRay)
      CellRayIdxSt(iPinRay, iCoreRay, 1) = FastRay%CellRayIdxSt(iPinRay, iCoreRay, 1)
      CellRayIdxSt(iPinRay, iCoreRay, 2) = FastRay%CellRayIdxSt(iPinRay, iCoreRay, 2)
      SurfIdx     (iPinRay, iCoreRay, 1) = FastRay%SurfIdx     (iPinRay, iCoreRay, 1)
      SurfIdx     (iPinRay, iCoreRay, 2) = FastRay%SurfIdx     (iPinRay, iCoreRay, 2)
      
      ipin = PinIdx(iPinRay, iCoreRay)
      
      icel     = Pin(ipin)%Cell(iz)
      FsrIdxSt = Pin(ipin)%FsrIdxSt
      
      DO iFSR = FastRay%Ray1DIdx(1, iPinRay, iCoreRay), FastRay%Ray1DIdx(2, iPinRay, iCoreRay)
        irsegidx = irsegidx + 1
        ireg     = FsrIdxSt + CellRay1D%LocalFsrIdx(iFSR) - 1
        
        FsrIdx(irsegidx, iCoreRay) = ireg
        
        tau = -CellRay1D%LenSeg(iFSR) * xst(ireg)
        
        OptLenList(irsegidx, iCoreRay) = tau
        ExpAppIdx (irsegidx, iCoreRay) = max(INT(tau), -40000)
        ExpAppIdx (irsegidx, iCoreRay) = min(0, ExpAppIdx(irsegidx, iCoreRay))
      END DO
    END DO
  END DO
ELSE IF (FastMocLv .EQ. 2) THEN
  FastRay   => RayInfo%FastCoreRayDat(iRotRay, iz)
  CellRay1D => RayInfo%CellRay1D
  
  DO iCoreRay = 1, nCoreRay
    nTotRaySeg (iCoreRay) = FastRay%nTotRaySeg (iCoreRay)
    nTotCellRay(iCoreRay) = FastRay%nTotCellRay(iCoreRay)
    
    DO iPinRay = 1, FastRay%nTotCellRay(iCoreRay)
      PinIdx      (iPinRay, iCoreRay)    = FastRay%PinIdx      (iPinRay, iCoreRay)
      CellRayIdxSt(iPinRay, iCoreRay, 1) = FastRay%CellRayIdxSt(iPinRay, iCoreRay, 1)
      CellRayIdxSt(iPinRay, iCoreRay, 2) = FastRay%CellRayIdxSt(iPinRay, iCoreRay, 2)
      SurfIdx     (iPinRay, iCoreRay, 1) = FastRay%SurfIdx     (iPinRay, iCoreRay, 1)
      SurfIdx     (iPinRay, iCoreRay, 2) = FastRay%SurfIdx     (iPinRay, iCoreRay, 2)
    END DO
    
    FastRaySeg => RayInfo%FastCoreRayDat(iRotRay, iz)%RaySeg(iCoreRay)
    
    DO iFSR = 1, FastRay%nTotRaySeg(iCoreRay)
      ireg = FastRaySeg%FsrIdx(iFSR)
      
      FsrIdx(iFSR, iCoreRay) = ireg
      
      tau = -FastRaySeg%LenSeg(iFSR) * xst(ireg)
      
      OptLenList(iFSR, iCoreRay) = tau
      ExpAppIdx (iFSR, iCoreRay) = max(INT(tau), -40000)
      ExpAppIdx (iFSR, iCoreRay) = min(0, ExpAppIdx(iFSR, iCoreRay))
    END DO
  END DO
END IF
! ----------------------------------------------------
DO iCoreRay = 1, nCoreRay
  DO iFSR = 1, nTotRaySeg(iCoreRay)
    DO ipol = 1, nPolarAng
      ExpAppPolar(ipol, iFSR, iCoreRay) = ExpA(ExpAppIdx(iFSR, iCoreRay), ipol) * optlenlist(iFSR, iCoreRay) + ExpB(ExpAppIdx(iFSR, iCoreRay), ipol)
      CONTINUE
    END DO
  END DO
END DO
! ----------------------------------------------------
ALLOCATE(phiobdPolar(1:nPolarAng))

DO irot = 1, 2
  PhiAnginSvIdx  = RayInfo%PhiAngInSvIdx (iRotRay ,irot)
  PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay ,irot)
  
  phiobdPolar(:) = PhiAngIn(:,PhiAnginSvIdx)

  jinc = 1; jbeg = 1; jend = nCoreRay
  IF (irot .eq. 2) THEN ! Backward Sweep
    jinc = -1; jbeg = nCoreRay; jend = 1
  END IF
  
  DO iCoreRay = jbeg, jend, jinc
    iazi    = CoreRay(RotRay(irotray)%RayIdx(iCoreRay))%iang
    iompazi = RotRay(irotray)%OmpRayIdx(iCoreRay)
    
    idir = RotRay(iRotRay)%DIR(iCoreRay)
    IF (irot .eq. 2) idir = mp(idir) ! Reverse the sweep direction
    
    nRaySeg = nTotRaySeg(iCoreRay)

    IF (idir .eq. 1) THEN ! Forward Sweep
      PhiAngOutPolar(:, 1) = phiobdPolar(:)
      
      DO ir = 1, nRaySeg
        ireg = FsrIdx(ir, iCoreRay)
        
        DO ipol = 1, nPolarAng
          phid = (PhiAngOutPolar(ipol, ir) - SrcAng1(ipol, ireg, iazi)) * ExpAppPolar(ipol, ir, iCoreRay)
          
          PhiAngOutPolar(ipol, ir+1) = PhiAngOutPolar(ipol, ir) - phid
          
          phi1a(ipol, ireg, iompazi) =  phi1a(ipol, ireg, iompazi) + phid
        END DO
      END DO
      
      phiobdPolar(:) = PhiAngOutPolar(:, nRaySeg+1)
      
      ! Surface
      IF (ljout) THEN
        DO ir = 1, nTotCellRay(iCoreRay)
          DO ipol = 1, nPolarANg
            wt = wtang(ipol, iazi)
            
            icel  = PinIdx(ir, iCoreRay)
            isurf = SurfIdx(ir, iCoreRay, 1)
            
            Jout(2, isurf, icel) = Jout(2, isurf, icel) + wt * PhiAngOutPolar(ipol, CellRayIdxSt(ir, iCoreRay, 1)+1)
            
            isurf = SurfIdx(ir, iCoreRay, 2)
            
            Jout(1, isurf, icel) = Jout(1, isurf, icel) + wt * PhiAngOutPolar(ipol, CellRayIdxSt(ir, iCoreRay, 2))
          END DO
        END DO
      END IF
    ELSE
      PhiAngOutPolar(:, nRaySeg+2) = phiobdPolar(:)
      
      ir = nRaySeg + 1
      
      DO ir1 = 1, nRaySeg
        ir   = ir - 1
        ireg = FsrIdx(ir, iCoreRay)
        
        DO ipol = 1, nPolarAng
          phid = (PhiAngOutPolar(ipol, ir + 2) - SrcAng2(ipol, ireg, iazi)) * ExpAppPolar(ipol, ir, iCoreRay)
          
          PhiAngOutPolar(ipol, ir+1) = PhiAngOutPolar(ipol, ir + 2) - phid
          
          phi2a(ipol, ireg, iompazi) =  phi2a(ipol, ireg, iompazi) + phid
        END DO
      END DO
      
      phiobdPolar(:) = PhiAngOutPolar(:, 2)
      
      IF (lJout) THEN
        DO ir = 1, nTotCellRay(iCoreRay)
          DO ipol = 1, nPolarAng
            wt = wtang(ipol, iazi)
            
            icel  = PinIdx (ir, iCoreRay)
            isurf = SurfIdx(ir, iCoreRay, 2)
            
            Jout(2, isurf, icel) = Jout(2, isurf, icel) + wt * PhiAngOutPolar(ipol, CellRayIdxSt(ir, iCoreRay, 2)+1)
            
            isurf = SurfIdx(ir, iCoreRay, 1)
            
            Jout(1, isurf, icel) = Jout(1, isurf, icel) + wt * PhiAngOutPolar(ipol, CellRayIdxSt(ir, iCoreRay, 1)+2)
          END DO
        END DO
      END IF
    END IF
  END DO
  
  PhiAngIn(:,PhiAngOutSvIdx) = phiobdPolar(:)
END DO
! ----------------------------------------------------
DEALLOCATE(phiobdPolar)

! Ray
NULLIFY (AsyRay)
NULLIFY (CoreRay)
NULLIFY (RotRay)
NULLIFY (CellRay)

! Geo.
NULLIFY (Asy)
NULLIFY (Pin)
NULLIFY (Cell)

! Loc.
NULLIFY (FsrIdx)
NULLIFY (ExpAppIdx)
NULLIFY (OptLenList)
NULLIFY (LenSeg)
NULLIFY (LocalFsrIdx)
NULLIFY (Phis)
NULLIFY (src)
NULLIFY (xst)
NULLIFY (jout)
NULLIFY (PhiAngOutPolar)
NULLIFY (PhiAngIn)
NULLIFY (ExpA)
NULLIFY (ExpB)
NULLIFY (wtang)
NULLIFY (wtsurf)
NULLIFY (Phi1a)
NULLIFY (Phi2a)
NULLIFY (SrcAng1)
NULLIFY (SrcAng2)
! ----------------------------------------------------

END SUBROUTINE RecTrackRotRayP1GM_AFSS
! ------------------------------------------------------------------------------------------------------------