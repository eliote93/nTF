! ------------------------------------------------------------------------------------------------------------
SUBROUTINE InitAFSS(RayInfo, CoreInfo, PE, nTracerCntl)

USE allocs
USE TYPEDEF, ONLY : RayInfo_Type, CoreInfo_Type, PE_TYPE, RotRayInfo_Type
USE CNTL,    ONLY : nTracerCntl_Type
USE IOUTIL,  ONLY : terminate
USE moc_mod, ONLY : OmpRayBeg, OmpRayEnd, OmpRayBegBd, OmpRayEndBd, OmpMap, trackingdat, nOmpAng

IMPLICIT NONE

TYPE (RayInfo_Type)     :: RayInfo
TYPE (CoreInfo_Type)    :: CoreInfo
TYPE (nTracerCntl_Type) :: nTracerCntl
TYPE (PE_TYPE)          :: PE
! ----------------------------------------------------
INTEGER :: nAziAng, nPolarAng, nRotray, nFSR, nThread, OmpTemp, ithr, nAzi, iRotRay, iazi, iRay, jazi

TYPE(RotRayInfo_Type), POINTER :: LocRotRay, RotRay(:)

INTEGER, POINTER, DIMENSION(:) :: nAziRotRay, nAziRotRay0, nSegRotRay
! ----------------------------------------------------

IF (.NOT. nTracerCntl%lAFSS) RETURN

nAziAng   = RayInfo%nAziAngle
nPolarAng = RayInfo%nPolarAngle
nRotRay   = RayInfo%nRotRay
RotRay   => RayInfo%RotRay

nFSR = CoreInfo%nCoreFsr

nThread = PE%nThread
! ----------------------------------------------------
IF (nTracerCntl%lScatBd) THEN
  ! Basic
  CALL dmalloc(nAziRotRay, nAziAng/2)
  CALL dmalloc(nSegRotRay, nAziAng/2)
  
  DO iRotRay = 1, nRotRay
    iazi = RotRay(iRotRay)%Ang1
    
    nAziRotRay(iazi) = nAziRotRay(iazi) + 1
    nSegRotRay(iazi) = nSegRotRay(iazi) + RotRay(iRotRay)%nSeg
  END DO
  
  CALL dmalloc(nAziRotRay0,   nThread)
  CALL dmalloc0(OmpRayBeg, 0, nThread)
  CALL dmalloc0(OmpRayEnd, 0, nThread)
  
  OmpRayBeg(0) = 1
  OmpRayEnd(0) = 0
  
  nAzi    = nAziAng / 2 / nThread
  OmpTemp = 0
  jazi    = 0
  
  DO ithr = 1, nThread
    OmpRayBeg(ithr) = OmpRayBeg(ithr-1) + OmpTemp
    
    DO iazi = 1, nAzi
      jazi = jazi + 1
      
      nAziRotRay0(ithr) = nAziRotRay0(ithr) + nAziRotRay(jazi)
    END DO
    
    OmpTemp = nAziRotRay0(ithr)
    OmpRayEnd(ithr) = OmpRayEnd(ithr-1) + OmpTemp
  END DO
  
  ! Parallel Ang
  DO iRotRay = 1, nRotRay
    LocRotRay => RotRay(iRotRay)
    
    LocRotRay%OmpAng1 = MOD(LocRotRay%Ang1, nAzi)
    
    IF (LocRotRay%OmpAng1 .EQ. 0) LocRotRay%OmpAng1 = nAzi
    
    LocRotRay%OmpAng2 = LocRotRay%OmpAng1 + nAzi
  END DO
  
  DO iRotRay = 1, nRotRay
    LocRotRay => RotRay(iRotRay)
    
    DO iRay = 1, LocRotRay%nRay, 2
      LocRotRay%OmpRayIdx(iRay)   = LocRotRay%OmpAng1
      LocRotRay%OmpRayIdx(iRay+1) = LocRotRay%OmpAng2
    END DO
  END DO
  
  CALL dmalloc(OmpMap, nThread, nAzi*2)
  
  DO ithr = 1, nThread
    DO iazi = 1, nAzi
      OmpTemp = OmpTemp + 1
      
      OmpMap(ithr, iazi)      = OmpTemp
      OmpMap(ithr, iazi+nAzi) = nAziAng - OmpTemp + 1
    END DO
  END DO
! ----------------------------------------------------
ELSE
  CALL dmalloc(nAziRotRay, nAziAng)
  CALL dmalloc(nSegRotRay, nAziAng)

  DO iRotRay = 1, nRotRay
    iazi = RotRay(iRotRay)%Ang1
    
    nAziRotRay(iazi) = nAziRotRay(iazi) + 1
    nSegRotRay(iazi) = nSegRotRay(iazi) + RotRay(iRotRay)%nSeg
  END DO
  
  CALL dmalloc(nAziRotRay0,   2*nThread)
  CALL dmalloc0(OmpRayBeg, 0, 2*nThread)
  CALL dmalloc0(OmpRayEnd, 0, 2*nThread)
  
  CALL dmalloc0(OmpRayBegBd, 1, 2, 0, nThread)
  CALL dmalloc0(OmpRayEndBd, 1, 2, 0, nThread)
  
  OmpRayBegBd(1,0) = 1
  OmpRayEndBd(1,0) = 0
  
  nAzi    = nAziAng / nThread / 2
  OmpTemp = 0
  jazi    = 0
  
  DO ithr = 1, nThread
    OmpRayBegBd(1, ithr) = OmpRayBegBd(1, ithr-1) + OmpTemp
    
    DO iazi = 1, nAzi
      jazi = jazi + 1
      
      nAziRotRay0(ithr) = nAziRotRay0(ithr) + nAziRotRay(jazi)
    END DO
    
    OmpTemp = nAziRotRay0(ithr)
    
    OmpRayEndBd(1, ithr) = OmpRayEndBd(1, ithr-1) + OmpTemp
  END DO
  
  nAziRotRay0 = 0
  jazi        = nAziAng + 1
  OmpTemp     = 0
  
  OmpRayBegBd(2, 0) = nRotRay + 1
  OmpRayEndBd(2, 0) = nRotRay
  
  DO ithr = 1, nThread
    OmpRayEndBd(2, ithr) = OmpRayEndBd(2, ithr-1) - OmpTemp
    
    DO iazi = 1, nAzi
      jazi = jazi - 1
      
      nAziRotRay0(ithr) = nAziRotRay0(ithr) + nAziRotRay(jazi)
    END DO
    
    OmpTemp = nAziRotRay0(ithr)
    
    OmpRayBegBd(2, ithr) = OmpRayBegBd(2, ithr-1) - OmpTemp
  END DO
  
  ! Parallel Ang
  DO iRotRay = 1, nRotRay
    LocRotRay => RotRay(iRotRay)
    
    IF (LocRotRay%nRay .GT. 1) THEN
      IF (LocRotRay%Ang1 .LE. nAziAng/2) THEN
        LocRotRay%OmpAng1 = MOD(LocRotRay%Ang1, nAzi)
        
        IF (LocRotRay%OmpAng1 .EQ. 0) LocRotRay%OmpAng1 = nAzi
        
        LocRotRay%OmpAng2 = LocRotRay%OmpAng1 + nAzi
      ELSE IF (LocRotRay%Ang1 .GT. nAziAng/2) THEN
        LocRotRay%OmpAng1 = MOD(nAziAng+1 - LocRotRay%Ang1, nAzi) + nAzi
        
        IF (LocRotRay%OmpAng1 .EQ. nAzi) LocRotRay%OmpAng1 = LocRotRay%OmpAng1 + nAzi
        
        LocRotRay%OmpAng2 = LocRotRay%OmpAng1 - nAzi
      END IF
    ELSE IF (LocRotRay%nRay .EQ. 1) THEN
      IF (LocRotRay%Ang1 .LE. nAziAng/2) THEN
        LocRotRay%OmpAng1 = MOD(LocRotRay%Ang1, nAzi)
        
        IF (LocRotRay%OmpAng1 .EQ. 0) LocRotRay%OmpAng1 = nAzi
        
        LocRotRay%OmpAng2 = 0
      ELSE IF (LocRotRay%Ang1 .GT. nAziAng/2) THEN
        
        LocRotRay%OmpAng1 = MOD(nAziAng+1-LocRotRay%Ang1, nAzi) + nAzi
        
        IF (LocRotRay%OmpAng1 .EQ. nAzi) LocRotRay%OmpAng1 = LocRotRay%OmpAng1 + nAzi
        
        LocRotRay%OmpAng2 = 0
      END IF
    END IF
  END DO
  
  DO iRotRay = 1, nRotRay
    LocRotRay => RotRay(iRotRay)
    
    IF (LocRotRay%nRay .GT. 1) THEN
      IF (MOD(LocRotRay%nRay, 2) .EQ. 0) THEN
        DO iRay = 1, LocRotRay%nRay, 2
          LocRotRay%OmpRayIdx(iRay)   = LocRotRay%OmpAng1
          LocRotRay%OmpRayIdx(iRay+1) = LocRotRay%OmpAng2
        END DO
      ELSE
        DO iRay = 1, LocRotRay%nRay-1, 2
          LocRotRay%OmpRayIdx(iRay)   = LocRotRay%OmpAng1
          LocRotRay%OmpRayIdx(iRay+1) = LocRotRay%OmpAng2
        END DO
        
        LocRotRay%OmpRayIdx(LocRotRay%nRay) = LocRotRay%OmpAng1
      END IF
    ELSE
      LocRotRay%OmpRayIdx(1) = LocRotRay%OmpAng1
    END IF
  END DO
  
  CALL dmalloc(OmpMap, nThread, nAzi*2)
  
  DO ithr = 1, nThread
    DO iazi = 1, nAzi
      OmpTemp = OmpTemp + 1
      
      OmpMap(ithr, iazi)      = OmpTemp
      OmpMap(ithr, iazi+nAzi) = nAziAng - OmpTemp + 1
    END DO
  END DO
END IF
! ----------------------------------------------------
IF (MOD(nAziANg / 2, nThread) .NE. 0) CALL terminate('WRONG_MOC_TRD')

nOmpAng = nAziAng / nThread

DO ithr = 1, nThread
  CALL dmalloc(TrackingDat(ithr)%phi1a, nPolarAng, nFsr, nOmpAng)
  CALL dmalloc(TrackingDat(ithr)%phi2a, nPolarAng, nFsr, nOmpAng)
END DO
! ----------------------------------------------------
! Ray
NULLIFY (LocRotRay)
NULLIFY (RotRay)

! Local
NULLIFY (nAziRotRay)
NULLIFY (nAziRotRay0)
NULLIFY (nSegRotRay)
! ----------------------------------------------------

END SUBROUTINE InitAFSS
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceGM_AFSS(RayInfo, CoreInfo, phis, PhiAngIn, xst, src, jout, iz, ljout, FastMocLv)

USE OMP_LIB
USE PARAM,   ONLY : FALSE, ZERO
USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type, Pin_Type, Cell_Type
USE Moc_Mod, ONLY : TrackingDat, nOmpAng, OmpRayBeg, OmpRayEnd, OmpRayBegBd, OmpRayEndBd, OmpMap, wtang
USE cntl,    ONLY : nTracerCntl
USE PE_MOD,  ONLY : PE

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:)     :: phis, xst, src
REAL, POINTER, DIMENSION(:,:)   :: PhiAngIn
REAL, POINTER, DIMENSION(:,:,:) :: jout

INTEGER :: iz
LOGICAL :: ljout
INTEGER, OPTIONAL :: FastMocLv
! ----------------------------------------------------
TYPE (Cell_Type), POINTER, DIMENSION(:) :: Cell
TYPE (Pin_Type),  POINTER, DIMENSION(:) :: Pin

INTEGER :: nPolarAng, nxy, nThread, iRotRay, ithr, FsrIdxSt, icel, iazi, ipol, OmpAng, ixy, iDir, iFSR, jFSR
! ----------------------------------------------------

nxy = CoreInfo%nxy
! ----------------------------------------------------
!$ call omp_set_dynamic(FALSE)
!$ call omp_set_num_threads(nThread)

DO ithr = 1, nThread
  TrackingDat(ithr)%PhiAngIn => PhiAngIn
  TrackingDat(ithr)%src      => src
  TrackingDat(ithr)%xst      => xst
  
  DO iazi = 1, nOmpAng
    TrackingDat(ithr)%phi1a(:,:,iazi) = ZERO
    TrackingDat(ithr)%phi2a(:,:,iazi) = ZERO
  END DO
  
  DO ixy = 1, nxy
    TrackingDat(ithr)%jout(:, :, ixy) = ZERO
  END DO
END DO

phis = ZERO
IF (ljout) Jout = ZERO
! ----------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ithr, iDir, iazi, iRotRay, iFSR, ipol, OmpAng)
ithr = 1
!$ ithr = omp_get_thread_num()+1
!$OMP BARRIER
IF (nTracerCntl%lScatBd) THEN
  DO iRotRay = OmpRayBeg(ithr), OmpRayEnd(ithr)
    CALL RecTrackRotRayGM_AFSS(RayInfo, CoreInfo, TrackingDat(ithr), ljout, iRotRay, iz, FastMocLv)
  END DO
ELSE
  DO iDir = 1, 2
    DO iRotRay = OmpRayBegBd(iDir, ithr), OmpRayEndBd(iDir, ithr)
      CALL RecTrackRotRayGM_AFSS(RayInfo, CoreInfo, TrackingDat(ithr), ljout, iRotRay, iz, FastMocLv)
    END DO
  END DO
END IF
!$OMP BARRIER
DO ithr = 1, nThread
  DO iazi = 1, nOmpAng
   DO iFSR = PE%myOmpFsrBeg(ithr), PE%myOmpFsrEnd(ithr)
      DO ipol = 1, nPolarAng
        OmpAng = OmpMap(ithr, iazi)
        
        phis(iFSR) = phis(iFSR) + wtang(ipol, OmpAng) * (TrackingDat(ithr)%phi1a(ipol, iFSR, iazi) + TrackingDat(ithr)%phi2a(ipol, iFSR, iazi))
      END DO
    END DO
  END DO
END DO

IF (ljout) THEN
  DO ithr = 1, nThread
    DO ixy = PE%myOmpNxyBeg(ithr), PE%myOmpNxyEnd(ithr)
      jout(:, :, ixy) = jout(:, :, ixy) + TrackingDat(ithr)%jout(:, :, ixy)
    END DO
  END DO
END IF
!$OMP END PARALLEL
! ----------------------------------------------------
Cell => CoreInfo%CellInfo
Pin  => CoreInfo%Pin

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ixy, FsrIdxSt, icel, iFSR, jFSR)
!$OMP DO
DO ixy = 1, nxy
  FsrIdxSt = Pin(ixy)%FsrIdxSt
  icel     = Pin(ixy)%Cell(iz)
  
  DO iFSR = 1, Cell(icel)%nFsr
    jFSR = FsrIdxSt + iFSR - 1
    
    phis(jFSR) = phis(jFSR) / xst(jFSR) / Cell(icel)%vol(iFSR) + src(jFSR)
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

NULLIFY (Cell)
NULLIFY (Pin)
! ----------------------------------------------------

END SUBROUTINE RayTraceGM_AFSS
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayGM_AFSS(RayInfo, CoreInfo, TrackingDat, ljout, iRotRay, iz, FastMocLv)

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
TYPE (Pin_Type),          POINTER, DIMENSION(:) :: Pin
TYPE (Asy_Type),          POINTER, DIMENSION(:) :: Asy
TYPE (Cell_Type),         POINTER, DIMENSION(:) :: Cell
TYPE (AsyRayInfo_type),   POINTER, DIMENSION(:) :: AsyRay
TYPE (CoreRayInfo_Type),  POINTER, DIMENSION(:) :: CoreRay
TYPE (RotRayInfo_Type),   POINTER, DIMENSION(:) :: RotRay

TYPE (FastCoreRayDat_Type), POINTER :: FastRay
TYPE (CellRayInfo_Type),    POINTER :: CellRay
TYPE (CellRayInfo_Type),    POINTER :: CellRay1D
TYPE (FastRaySegDat_Type),  POINTER :: FastRaySeg

INTEGER, POINTER, DIMENSION(:)   :: LocalFsrIdx
INTEGER, POINTER, DIMENSION(:,:) :: FsrIdx,  ExpAppIdx

REAL, POINTER, DIMENSION(:)     :: LenSeg, phis, src, xst
REAL, POINTER, DIMENSION(:,:)   :: OptLenList, PhiAngOutPolar, PhiAngIn, EXPA, EXPB, wtang
REAL, POINTER, DIMENSION(:,:,:) :: ExpAppPolar, jout, wtsurf, phi1a, phi2a

INTEGER :: mp(2)
INTEGER :: iazi, ipol, iCoreRay, jCoreRay, iAsyRay, jAsyRay, iceray, irayseg, irot, idir, nCoreRay, nAsyRay, nPinRay, nRaySeg, FsrIdxSt, nPolarAng, nAziAng
INTEGER :: ipin, icel, iasy, ireg, isurf, irsegidx, icellrayidx, PhiAnginSvIdx, PhiAngOutSvIdx, jbeg, jend, jinc, ir, ir1, iOmpAzi, ibcel, iFSR, iPinRay

INTEGER, DIMENSION(nMaxCoreRay)                 :: nTotRaySeg, nTotCellRay
INTEGER, DIMENSION(nMaxCellRay, nMaxCoreRay)    :: PinIdx
INTEGER, DIMENSION(nMaxCellRay, nMaxCoreRay, 2) :: CellRayIdxSt, SurfIdx

REAL :: tau, phiobd, phid, wt, wt2(4)

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
EXPA           => TrackingDat%EXPA
EXPB           => TrackingDat%EXPB
Wtang          => TrackingDat%wtang
wtsurf         => TrackingDat%wtsurf
! ----------------------------------------------------
nCoreRay = RotRay(iRotRay)%nRay

IF (.NOT. lFast) THEN
  DO iCoreRay = 1, nCoreRay
    irsegidx    = 0
    icellrayidx = 0
    
    jCoreRay = RotRay(iRotRay)%RayIdx(iCoreRay)
    nAsyRay  = CoreRay(jCoreRay)%nRay
    
    DO iAsyRay = 1, nAsyRay
      jAsyRay = CoreRay(jCoreRay)%AsyRayIdx(iAsyRay)
      iasy    = CoreRay(jCoreRay)%AsyIdx   (iAsyRay)
      
      IF (iasy .EQ. 0)  CYCLE ! Skip Dummy Assembly
      
      nPinRay = AsyRay(jAsyRay)%nCellRay
      
      DO iPinRay = 1, nPinRay
        ipin   = AsyRay(jAsyRay)%PinIdx   (iPinRay) ! Local Pin Idx(within Assembly)
        iceray = AsyRay(jAsyRay)%PinRayIdx(iPinRay)
        
        ipin = Asy(iAsy)%GlobalPinIdx(ipin) ! Global Pin Index
        
        icel     = Pin(ipin)%Cell(iz)
        FsrIdxSt = Pin(ipin)%FsrIdxSt
        
        ibcel    = Cell(icel)%basecellstr
        CellRay => Cell(ibcel)%CellRay(iceray)
        
        icellrayidx = icellrayidx + 1
        
        PinIdx      (icellrayidx, iCoreRay)    = ipin
        CellRayIdxSt(icellrayidx, iCoreRay, 2) = irsegidx + 1
        
        nRaySeg      = CellRay%nSeg
        LocalFsrIdx => CellRay%LocalFsrIdx
        LenSeg      => CellRay%LenSeg
        
        DO iRaySeg = 1, nRaySeg
          ireg = FsrIdxSt + CellRay%LocalFsrIdx(iRaySeg) - 1
          tau  = -LenSeg(iRaySeg)         * xst(ireg)   !
          tau  = -CellRay%LenSeg(iRaySeg) * xst(ireg)   !
          
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
ELSEIF (FastMocLV .EQ. 1) THEN
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
         
         ireg = FsrIdxSt + CellRay1D%LocalFsrIdx(iFSR) - 1
         
         FsrIdx(irsegidx, iCoreRay) = ireg
         
         tau = -CellRay1D%LenSeg(iFSR) * xst(ireg)
         
         OptLenList(irsegidx, iCoreRay) = tau
         ExpAppIdx (irsegidx, iCoreRay) = max(INT(tau), -40000)
         ExpAppIdx (irsegidx, iCoreRay) = min(0, ExpAppIdx(irsegidx, iCoreRay))
       END DO
     END DO
   END DO
ELSEIF (FastMocLv .EQ. 2) THEN
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
    DO ipol = 1, nPolarANg
      ExpAppPolar(ipol, iFSR, iCoreRay) = expa(ExpAppIdx(iFSR, iCoreRay), ipol) * optlenlist(iFSR, iCoreRay) + expb(ExpAppIdx(iFSR, iCoreRay), ipol)
      CONTINUE
    END DO
  END DO
END DO

ALLOCATE(phiobdPolar(1:nPolarAng))
! ----------------------------------------------------
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
    
    IF (lJout) wt2(1:4) = wtsurf(ipol, iazi, 1:4)

    nRaySeg = nTotRaySeg(iCoreRay)

    IF (idir .EQ. 1) THEN  ! Forward Sweep
      PhiAngOutPolar(:, 1) = phiobdPolar(:)
      
      DO ir = 1, nRaySeg
        ireg = FsrIdx(ir, iCoreRay)
        
        DO ipol = 1, nPolarAng
          phid = (PhiAngOutPolar(ipol, ir) - src(ireg)) * ExpAppPolar(ipol, ir, iCoreRay)
          
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
            
            Jout(2, isurf, icel) = Jout(2, isurf, icel) + PhiAngOutPolar(ipol, CellRayIdxSt(ir, iCoreRay, 1)+1) * wt
            Jout(3, isurf, icel) = Jout(3, isurf, icel) + PhiAngOutPolar(ipol, CellRayIdxSt(ir, iCoreRay, 1)+1) * wt2(isurf)
            
            isurf = SurfIdx(ir, iCoreRay, 2)
            
            Jout(1, isurf, icel) = Jout(1, isurf, icel) + PhiAngOutPolar(ipol, CellRayIdxSt(ir, iCoreRay, 2)) * wt
            Jout(3, isurf, icel) = Jout(3, isurf, icel) + PhiAngOutPolar(ipol, CellRayIdxSt(ir, iCoreRay, 2)) * wt2(isurf)
          END DO
        END DO
      END IF
    ELSE
      PhiAngOutPolar(:, nRaySeg+2) = phiobdPolar(:)
      
      ir = nRaySeg + 1
      
      DO ir1 = 1, nRaySeg
        ir = ir - 1
        
        ireg = FsrIdx(ir, iCoreRay)
        
        DO ipol = 1, nPolarAng
          phid = (PhiAngOutPolar(ipol, ir + 2) - src(ireg)) * ExpAppPolar(ipol, ir, iCoreRay)
          
          PhiAngOutPolar(ipol, ir+1) = PhiAngOutPolar(ipol, ir + 2) - phid
          
          phi2a(ipol, ireg, iompazi) =  phi2a(ipol, ireg, iompazi) + phid
        END DO
      END DO
      
      phiobdPolar(:) = PhiAngOutPolar(:, 2)
      
      IF (lJout) THEN
        DO ir = 1, nTotCellRay(iCoreRay)
          DO ipol = 1, nPolarAng
            wt = wtang(ipol, iazi)
            
            icel  = PinIdx(ir, iCoreRay)
            isurf = SurfIdx(ir, iCoreRay, 2)
            
            Jout(2, isurf, icel) = Jout(2, isurf, icel) + PhiAngOutPolar(ipol, CellRayIdxSt(ir, iCoreRay, 2)+1) * wt
            Jout(3, isurf, icel) = Jout(3, isurf, icel) + PhiAngOutPolar(ipol, CellRayIdxSt(ir, iCoreRay, 2)+1) * wt2(isurf)
            
            isurf = SurfIdx(ir, iCoreRay, 1)
            
            Jout(1, isurf, icel) = Jout(1, isurf, icel) + PhiAngOutPolar(ipol, CellRayIdxSt(ir, iCoreRay, 1)+2) * wt
            Jout(3, isurf, icel) = Jout(3, isurf, icel) + PhiAngOutPolar(ipol, CellRayIdxSt(ir, iCoreRay, 1)+2) * wt2(isurf)
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
NULLIFY (EXPA)
NULLIFY (EXPB)
NULLIFY (wtang)
NULLIFY (wtsurf)
NULLIFY (Phi1a)
NULLIFY (Phi2a)
! ----------------------------------------------------

END SUBROUTINE RecTrackRotRayGM_AFSS
! ------------------------------------------------------------------------------------------------------------
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

INTEGER :: nAziAng, nPolarAng, nxy, nThread, iRotRay, ithr, FsrIdxSt, icel, ireg, iazi, ipol, OmpAng, ixy, iDir, iFSR, jFSR
REAL :: wttmp, tempsrc, ONETHREE, ONEFIVE, ONESEVEN
! ----------------------------------------------------

nAziAng   = RayInfo%nAziAngle
nPolarAng = RayInfo%nPolarAngle

nxy = CoreInfo%nxy
! ----------------------------------------------------
!$ call omp_set_dynamic(FALSE)
!$ call omp_set_num_threads(nThread)

DO ithr = 1, nThread
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
  DO ithr = 1, nThread
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
  DO ithr = 1, nThread
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
  DO ithr = 1, nThread
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
  DO ithr = 1, nThread
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
REAL, POINTER, DIMENSION(:,:)   :: OptLenList, PhiAngOutPolar, PhiAngIn, EXPA, EXPB, wtang
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
EXPA           => TrackingDat%EXPA
EXPB           => TrackingDat%EXPB
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
      ExpAppPolar(ipol, iFSR, iCoreRay) = expa(ExpAppIdx(iFSR, iCoreRay), ipol) * optlenlist(iFSR, iCoreRay) + expb(ExpAppIdx(iFSR, iCoreRay), ipol)
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
NULLIFY (EXPA)
NULLIFY (EXPB)
NULLIFY (wtang)
NULLIFY (wtsurf)
NULLIFY (Phi1a)
NULLIFY (Phi2a)
NULLIFY (SrcAng1)
NULLIFY (SrcAng2)
! ----------------------------------------------------

END SUBROUTINE RecTrackRotRayP1GM_AFSS
! ------------------------------------------------------------------------------------------------------------