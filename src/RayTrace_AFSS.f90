SUBROUTINE RayTraceGM_AFSS(RayInfo, CoreInfo, phis, PhiAngIn, xst, src, jout, iz, ljout, FastMocLv)

USE TIMER
USE ALLOCS
USE OMP_LIB
USE PARAM,   ONLY : TRUE, FALSE, ZERO
USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type, TrackingDat_Type, Pin_Type, Cell_Type, AziAngleInfo_Type, PolarAngle_Type
USE Moc_Mod, ONLY : nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay, Expa, Expb, TrackingDat, ApproxExp, wtang, wtsurf, Phia1g
USE geom,    ONLY : nbd
USE cntl,    ONLY : nTracerCntl
USE PE_MOD,  ONLY : PE
USE IOUTIL,  ONLY : terminate

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:)     :: phis, xst, src
REAL, POINTER, DIMENSION(:,:)   :: PhiAngIn
REAL, POINTER, DIMENSION(:,:,:) :: jout
INTEGER :: iz
LOGICAL :: ljout
INTEGER, OPTIONAL :: FastMocLv

LOGICAL, SAVE :: lfirst
DATA lfirst /TRUE/

TYPE (AziAngleInfo_Type), POINTER, DIMENSION(:) :: AziAng
TYPE (PolarAngle_Type),   POINTER, DIMENSION(:) :: PolarAng
TYPE (Cell_Type),         POINTER, DIMENSION(:) :: Cell
TYPE (Pin_Type),          POINTER, DIMENSION(:) :: Pin

INTEGER :: nAziAng, nPolarAng, nPhiAngSv, nRotray, nFsr, nxy, nThread
INTEGER :: iRotRay, ithr, FsrIdxSt, icel, ireg, iazi, ipol, i, j, k, l, m
REAL :: wttmp

LOGICAL, SAVE :: lScatBd

INTEGER, SAVE :: OmpTemp, OmpAng, nOmpAng

INTEGER, POINTER, SAVE, DIMENSION(:)   :: nAziRotRay, nAziRotRay0, nSegRotRay, OmpRayBeg, OmpRayEnd, OmpRayList
INTEGER, POINTER, SAVE, DIMENSION(:,:) :: OmpRayBegBd, OmpRayEndBd, OmpMap, OmpRayNum
! ----------------------------------------------------

AziAng   => RayInfo%AziAngle
PolarAng => RayInfo%PolarAngle
nAziAng   = RayInfo%nAziAngle
nPolarAng = RayInfo%nPolarAngle
nPhiAngSv = RayInfo%nPhiAngSv
nRotRay   = RayInfo%nRotRay

nFsr = CoreInfo%nCoreFsr
nxy  = CoreInfo%nxy

nThread = PE%nThread
! ----------------------------------------------------
IF (lfirst) THEN
  lScatBd = nTracerCntl%lScatBd
  
  IF (lScatBd) THEN
    CALL dmalloc(nAziRotRay,  nAziAng/2)
    CALL dmalloc(nAziRotRay0, nThread)
    CALL dmalloc(nSegRotRay,  nAziAng/2)
        
    DO i = 1, nRotRay
      j = RayInfo%RotRay(i)%Ang1
      
      nAziRotRay(j) = nAziRotRay(j) + 1
      nSegRotRay(j) = nSegRotRay(j) + RayInfo%RotRay(i)%nSeg
    END DO
    
    CALL dmalloc0(OmpRayBeg, 0, nThread)
    CALL dmalloc0(OmpRayEnd, 0, nThread)
    
    k = nAziAng/2/nThread
    nAziRotRay0(:) = 0
    OmpRayBeg(0) = 1
    OmpTemp = 0
    OmpRayEnd(0) = 0
    l = 0
    
    DO i = 1, nThread
      OmpRayBeg(i) = OmpRayBeg(i-1) + OmpTemp
      
      DO j = 1, k
        l = l + 1
        
        nAziRotRay0(i) = nAziRotRay0(i) + nAziRotRay(l)
      END DO
      
      OmpTemp = nAziRotRay0(i)
      OmpRayEnd(i) = OmpRayEnd(i-1) + OmpTemp
    END DO
    
    ! Parallel Ang
    DO i = 1, nRotRay
      RayInfo%RotRay(i)%OmpAng1 = MOD(RayInfo%RotRay(i)%Ang1, k)
      
      IF (RayInfo%RotRay(i)%OmpAng1 .EQ. 0) RayInfo%RotRay(i)%OmpAng1 = k
      
      RayInfo%RotRay(i)%OmpAng2 = RayInfo%RotRay(i)%OmpAng1 + k
    END DO
    
    DO i = 1, nRotRay
      DO j = 1, RayInfo%RotRay(i)%nRay, 2
        RayInfo%RotRay(i)%OmpRayIdx(j)   = RayInfo%RotRay(i)%OmpAng1
        RayInfo%RotRay(i)%OmpRayIdx(j+1) = RayInfo%RotRay(i)%OmpAng2
      END DO
    END DO
    
    CALL dmalloc(OmpMap, nThread, k*2)
    
    DO i = 1, nThread
      DO j = 1, k
        OmpTemp = OmpTemp + 1
        
        OmpMap(i,j)   = OmpTemp
        OmpMap(i,j+k) = nAziAng - OmpTemp + 1
      END DO
    END DO
    nOmpAng = nAziAng/nThread
    
  ELSE
    CALL dmalloc(nAziRotRay,  nAziAng)
    CALL dmalloc(nAziRotRay0, 2*nThread)
    CALL dmalloc(nSegRotRay,  nAziAng)
  
    DO i = 1, nRotRay
      j = RayInfo%RotRay(i)%Ang1
      
      nAziRotRay(j) = nAziRotRay(j) + 1
      nSegRotRay(j) = nSegRotRay(j) + RayInfo%RotRay(i)%nSeg
    END DO
    
    CALL dmalloc0(OmpRayBeg, 0, 2*nThread)
    CALL dmalloc0(OmpRayEnd, 0, 2*nThread)
    
    CALL dmalloc0(OmpRayBegBd, 1, 2, 0, nThread)
    CALL dmalloc0(OmpRayEndBd, 1, 2, 0, nThread)
    
    k = nAziAng/nThread/2
    nAziRotRay0(:) = 0
    OmpRayBegBd(1,0) = 1
    OmpTemp = 0
    OmpRayEndBd(1,0) = 0
    l = 0
    
    DO i = 1, nThread
      OmpRayBegBd(1,i) = OmpRayBegBd(1,i-1) + OmpTemp
      
      DO j = 1, k
        l = l + 1
        
        nAziRotRay0(i) = nAziRotRay0(i) + nAziRotRay(l)
      END DO
      OmpTemp = nAziRotRay0(i)
      OmpRayEndBd(1,i) = OmpRayEndBd(1,i-1) + OmpTemp
    END DO
    
    nAziRotRay0(:) = 0
    l = nAziAng + 1
    OmpRayBegBd(2,0) = nRotRay + 1
    OmpTemp = 0
    OmpRayEndBd(2,0) = nRotRay
    
    DO i = 1, nThread
      OmpRayEndBd(2,i) = OmpRayEndBd(2,i-1) - OmpTemp
      
      DO j = 1, k
        l = l - 1
        
        nAziRotRay0(i) = nAziRotRay0(i) + nAziRotRay(l)
      END DO
      
      OmpTemp = nAziRotRay0(i)
      OmpRayBegBd(2,i) = OmpRayBegBd(2,i-1) - OmpTemp
    END DO
    
    ! Parallel Ang
    DO i = 1, nRotRay
      IF (RayInfo%RotRay(i)%nRay .GT. 1) THEN
        IF (RayInfo%RotRay(i)%Ang1 .LE. nAziAng/2) THEN
          RayInfo%RotRay(i)%OmpAng1 = MOD(RayInfo%RotRay(i)%Ang1, k)
          
          IF (RayInfo%RotRay(i)%OmpAng1 .EQ. 0) RayInfo%RotRay(i)%OmpAng1 = k
          
          RayInfo%RotRay(i)%OmpAng2 = RayInfo%RotRay(i)%OmpAng1 + k
        ELSE IF (RayInfo%RotRay(i)%Ang1 .GT. nAziAng/2) THEN
          RayInfo%RotRay(i)%OmpAng1 = MOD(nAziAng+1-RayInfo%RotRay(i)%Ang1, k) + k
          
          IF (RayInfo%RotRay(i)%OmpAng1 .EQ. k) RayInfo%RotRay(i)%OmpAng1 = RayInfo%RotRay(i)%OmpAng1 + k
          
          RayInfo%RotRay(i)%OmpAng2 = RayInfo%RotRay(i)%OmpAng1 - k
        END IF
      ELSE IF (RayInfo%RotRay(i)%nRay .EQ. 1) THEN
        IF (RayInfo%RotRay(i)%Ang1 .LE. nAziAng/2) THEN
          RayInfo%RotRay(i)%OmpAng1 = MOD(RayInfo%RotRay(i)%Ang1, k)
          
          IF (RayInfo%RotRay(i)%OmpAng1 .EQ. 0) RayInfo%RotRay(i)%OmpAng1 = k
          
          RayInfo%RotRay(i)%OmpAng2 = 0
        ELSE IF (RayInfo%RotRay(i)%Ang1 .GT. nAziAng/2) THEN
          
          RayInfo%RotRay(i)%OmpAng1 = MOD(nAziAng+1-RayInfo%RotRay(i)%Ang1, k) + k
          
          IF (RayInfo%RotRay(i)%OmpAng1 .EQ. k) RayInfo%RotRay(i)%OmpAng1 = RayInfo%RotRay(i)%OmpAng1 + k
          
          RayInfo%RotRay(i)%OmpAng2 = 0
        END IF
      END IF
    END DO
    
    DO i = 1, nRotRay
      IF (RayInfo%RotRay(i)%nRay .GT. 1) THEN
        IF (MOD(RayInfo%RotRay(i)%nRay, 2) .EQ. 0) THEN
          DO j = 1, RayInfo%RotRay(i)%nRay, 2
            RayInfo%RotRay(i)%OmpRayIdx(j)   = RayInfo%RotRay(i)%OmpAng1
            RayInfo%RotRay(i)%OmpRayIdx(j+1) = RayInfo%RotRay(i)%OmpAng2
          END DO
        ELSE
          DO j = 1, RayInfo%RotRay(i)%nRay-1, 2
            RayInfo%RotRay(i)%OmpRayIdx(j)   = RayInfo%RotRay(i)%OmpAng1
            RayInfo%RotRay(i)%OmpRayIdx(j+1) = RayInfo%RotRay(i)%OmpAng2
          END DO
          
          RayInfo%RotRay(i)%OmpRayIdx(RayInfo%RotRay(i)%nRay) = RayInfo%RotRay(i)%OmpAng1
        END IF
      ELSE
        RayInfo%RotRay(i)%OmpRayIdx(1) = RayInfo%RotRay(i)%OmpAng1
      END IF
    END DO
    
    CALL dmalloc(OmpMap, nThread, k*2)
    
    DO i = 1, nThread
      DO j = 1, k
        OmpTemp = OmpTemp + 1
        
        OmpMap(i,j)   = OmpTemp
        OmpMap(i,j+k) = nAziAng - OmpTemp + 1
      END DO
    END DO
    
    CALL dmalloc(OmpRayList, nThread)
    CALL dmalloc(OmpRayNum, nThread, nRotRay)
    
    nOmpAng = nAziAng/nThread
    
    DO i = 1, nThread
      DO j = 1, 2
        OmpRayList(i) = OmpRayList(i) + OmpRayEndBd(j,i) - OmpRayBegBd(j,i) + 1
      END DO
    END DO
    
    DO i = 1, nThread
      m = 0
      
      DO j = 1, 2
        l = 0
        
        DO k = 1, OmpRayEndBd(j,i) - OmpRayBegBd(j,i) + 1
          m = m + 1
          
          OmpRayNum(i,m) = OmpRayBegBd(j,i) + l
          
          l = l + 1
        END DO
      END DO
    END DO
  END IF
  
  IF (MOD(nAziANg/2, nThread) .NE. 0) CALL terminate('WRONG_MOC_TRD')
END IF
! ----------------------------------------------------
IF (lfirst) THEN
  lFirst = FALSE
  
  CALL ApproxExp(RayInfo%PolarAngle, nPolarAng)
  
  DO ithr = 1, nThread
    IF (TrackingDat(ithr)%lAlloc) CYCLE
    
    CALL dmalloc(TrackingDat(ithr)%FsrIdx,         nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%ExpAppIdx,      nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%OptLenList,     nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%ExpAppPolar,    nPolarAng, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%PhiAngOutPolar, nPolarAng, nMaxRaySeg+2)
    CALL dmalloc(TrackingDat(ithr)%phi1a,          nPolarAng, nFsr, nOmpAng)
    CALL dmalloc(TrackingDat(ithr)%phi2a,          nPolarAng, nFsr, nOmpAng)
    CALL dmalloc(TrackingDat(ithr)%Jout,           3, nbd, nxy)
    
    TrackingDat(ithr)%Expa => Expa
    TrackingDat(ithr)%Expb => Expb
    TrackingDat(ithr)%lAlloc = TRUE
  END DO
  
  CALL dmalloc(wtang, nPolarAng, nAziAng)
  
  DO ipol = 1, nPolarAng
    wttmp = PolarAng(ipol)%weight * PolarAng(ipol)%sinv
    
    DO iazi = 1, nAziAng
      wtang(ipol, iazi) = wttmp  * AziAng(iazi)%weight * AziAng(iazi)%del
      
      wtsurf(ipol, iazi, 1) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv)
      wtsurf(ipol, iazi, 3) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv)
      wtsurf(ipol, iazi, 2) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv)
      wtsurf(ipol, iazi, 4) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv)
    END DO
  END DO
END IF
! ----------------------------------------------------
!$  call omp_set_dynamic(FALSE)
!$  call omp_set_num_threads(nThread)

DO ithr = 1, nThread
  TrackingDat(ithr)%PhiAngIn => PhiAngIn
  TrackingDat(ithr)%src      => src
  TrackingDat(ithr)%xst      => xst
  TrackingDat(ithr)%wtang    => wtang
  TrackingDat(ithr)%wtsurf   => wtsurf
END DO

phis =  ZERO
IF (ljout) Jout = ZERO
! ----------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(irotray, ithr, i, j, k, l, iazi, ipol, OmpAng)
ithr = 1
!$  ithr = omp_get_thread_num()+1
j = ithr

DO iazi = 1, nOmpAng
  TrackingDat(j)%phi1a(:,:,iazi) = ZERO
  TrackingDat(j)%phi2a(:,:,iazi) = ZERO
END DO

DO i = 1, nxy
  TrackingDat(j)%jout(:, :, i) = ZERO
END DO
!$OMP BARRIER
IF (lScatBd) THEN
  DO i = OmpRayBeg(ithr), OmpRayEnd(ithr)
    irotray = i
    CALL RecTrackRotRayGM_AFSS(RayInfo, CoreInfo, TrackingDat(ithr), ljout, irotray, iz, FastMocLv)
  END DO
ELSE
  DO i = 1, 2
    DO j = OmpRayBegBd(i, ithr), OmpRayEndBd(i, ithr)
      irotray = j
      CALL RecTrackRotRayGM_AFSS(RayInfo, CoreInfo, TrackingDat(ithr), ljout, irotray, iz, FastMocLv)
    END DO
  END DO
END IF
!$OMP BARRIER
DO j = 1, nThread
  DO iazi = 1, nOmpAng
   DO i = PE%myOmpFsrBeg(ithr), PE%myOmpFsrEnd(ithr)
      DO ipol = 1, nPolarAng
        OmpAng = OmpMap(j, iazi)
        
        phis(i) = phis(i) + wtang(ipol, OmpAng) * (TrackingDat(j)%phi1a(ipol, i, iazi) + TrackingDat(j)%phi2a(ipol, i, iazi))
      END DO
    END DO
  END DO
END DO

IF (ljout) THEN
  DO j = 1, nThread
    DO i = PE%myOmpNxyBeg(ithr), PE%myOmpNxyEnd(ithr)
      jout(:, :, i) = jout(:, :, i) + TrackingDat(j)%jout(:, :, i)
    END DO
  END DO
END IF
!$OMP END PARALLEL
! ----------------------------------------------------
Cell => CoreInfo%CellInfo
Pin  => CoreInfo%Pin

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, j, FsrIdxSt, icel, ireg)
!$OMP DO
DO l = 1, nxy
  FsrIdxSt = Pin(l)%FsrIdxSt
  icel = Pin(l)%Cell(iz)
  
  DO j = 1, Cell(icel)%nFsr
    ireg = FsrIdxSt + j - 1
    
    phis(ireg) = phis(ireg)/xst(ireg)/Cell(icel)%vol(j) + src(ireg)
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
! ----------------------------------------------------
NULLIFY (AziAng)
NULLIFY (PolarAng)
NULLIFY (Cell)
NULLIFY (Pin)
! ----------------------------------------------------

END SUBROUTINE RayTraceGM_AFSS
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayGM_AFSS(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, FastMocLv)

USE PARAM,   ONLY : TRUE, FALSE, ZERO
USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type, Pin_Type, Asy_Type, AsyInfo_Type, PinInfo_Type, Cell_Type, AziAngleInfo_Type, PolarAngle_Type, ModRayInfo_type, AsyRayInfo_type,  &
                    CoreRayInfo_Type, RotRayInfo_Type, CellRayInfo_type, FastCoreRayDat_Type, TrackingDat_Type, FastRaySegDat_Type
USE Moc_Mod, ONLY : nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay

IMPLICIT NONE

TYPE (RayInfo_Type),     INTENT(INOUT) :: RayInfo
TYPE (CoreInfo_Type),    INTENT(INOUT) :: CoreInfo
TYPE (TrackingDat_Type), INTENT(INOUT) :: TrackingDat

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz, FastMocLv

TYPE (Pin_Type),          POINTER, DIMENSION(:) :: Pin
TYPE (Asy_Type),          POINTER, DIMENSION(:) :: Asy
TYPE (PinInfo_Type),      POINTER, DIMENSION(:) :: PinInfo
TYPE (Cell_Type),         POINTER, DIMENSION(:) :: Cell
TYPE (AziAngleInfo_Type), POINTER, DIMENSION(:) :: AziAng
TYPE (PolarAngle_Type),   POINTER, DIMENSION(:) :: PolarAng
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
INTEGER :: iazi, ipol, iCoreRay, iasyray, iceray, irayseg, irot, itype, idir, nRotRay, nCoreRay, nAsyRay, nPinRay, nRaySeg, FsrIdxSt, nPolarAng, nAziAng, nPhiAngSv
INTEGER :: ipin, icel, iasy, ireg, isurf, irsegidx, icellrayidx, PhiAnginSvIdx, PhiAngOutSvIdx, nFsr, nxy, i, j, k, l, m, jbeg, jend, jinc, ir, ir1, iOmpAzi, ibcel

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
      
      IF (iasy .EQ. 0)  CYCLE ! Skip Dummy Assembly
      
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
          tau  = -LenSeg(iRaySeg) * xst(ireg)   !
          tau  = -CellRay%LenSeg(iRaySeg) * xst(ireg)   !
          
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
ELSEIF (FastMocLV .EQ. 1) THEN
   FastRay   => RayInfo%FastCoreRayDat(i, iz)
   CellRay1D => RayInfo%CellRay1D
   
   DO j = 1, nCoreRay
     irsegidx = 0
     
     nTotRaySeg (j) = FastRay%nTotRaySeg (j)
     nTotCellRay(j) = FastRay%nTotCellRay(j)
     
     DO l = 1, FastRay%nTotCellRay(j)
       PinIdx      (l, j) = FastRay%PinIdx         (l, j)
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
ELSEIF (FastMocLv .EQ. 2) THEN
  FastRay   => RayInfo%FastCoreRayDat(i, iz)
  CellRay1D => RayInfo%CellRay1D
  
  DO j = 1, nCoreRay
    nTotRaySeg (j) = FastRay%nTotRaySeg (j)
    nTotCellRay(j) = FastRay%nTotCellRay(j)
    
    DO l = 1, FastRay%nTotCellRay(j)
      PinIdx      (l, j) = FastRay%PinIdx         (l, j)
      CellRayIdxSt(l, j, 1) = FastRay%CellRayIdxSt(l, j, 1)
      CellRayIdxSt(l, j, 2) = FastRay%CellRayIdxSt(l, j, 2)
      SurfIdx     (l, j, 1) = FastRay%SurfIdx     (l, j, 1)
      SurfIdx     (l, j, 2) = FastRay%SurfIdx     (l, j, 2)
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
    DO ipol = 1, nPolarANg
      ExpAppPolar(ipol, l, j) = expa(ExpAppIdx(l, j), ipol) * optlenlist(l, j) + expb(ExpAppIdx(l, j), ipol)
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
  
  DO j = jbeg, jend, jinc
    iazi    = CoreRay(RotRay(irotray)%RayIdx(j))%iang
    iompazi = RotRay(irotray)%OmpRayIdx(j)
    
    idir = RotRay(i)%DIR(j)
    IF (irot .eq. 2) idir = mp(idir) ! Reverse the sweep direction
    
    IF (lJout) wt2(1:4) = wtsurf(ipol, iazi, 1:4)

    nRaySeg = nTotRaySeg(j)

    IF (idir .EQ. 1) THEN  ! Forward Sweep
      PhiAngOutPolar(:, 1) = phiobdPolar(:)
      
      DO ir = 1, nRaySeg
        ireg = FsrIdx(ir, j)
        
        DO ipol = 1, nPolarAng
          phid = (PhiAngOutPolar(ipol, ir) - src(ireg)) * ExpAppPolar(ipol, ir, j)
          
          PhiAngOutPolar(ipol, ir+1) = PhiAngOutPolar(ipol, ir) - phid
          
          phi1a(ipol, ireg, iompazi) =  phi1a(ipol, ireg, iompazi) + phid
        END DO
      END DO
      
      phiobdPolar(:) = PhiAngOutPolar(:, nRaySeg+1)
      
      ! Surface
      IF (ljout) THEN
        DO ir = 1, nTotCellRay(j)
          DO ipol = 1, nPolarANg
            wt = wtang(ipol, iazi)
            
            icel  = PinIdx(ir, j)
            isurf = SurfIdx(ir, j, 1)
            
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
      
      DO ir1 = 1, nRaySeg
        ir = ir - 1
        
        ireg = FsrIdx(ir, j)
        
        DO ipol = 1, nPolarAng
          phid = (PhiAngOutPolar(ipol, ir + 2) - src(ireg)) * ExpAppPolar(ipol, ir, j)
          
          PhiAngOutPolar(ipol, ir+1) = PhiAngOutPolar(ipol, ir + 2) - phid
          
          phi2a(ipol, ireg, iompazi) =  phi2a(ipol, ireg, iompazi) + phid
        END DO
      END DO
      
      phiobdPolar(:) = PhiAngOutPolar(:, 2)
      
      IF (lJout) THEN
        DO ir = 1, nTotCellRay(j)
          DO ipol = 1, nPolarAng
            wt = wtang(ipol, iazi)
            
            icel   = PinIdx(ir, j)
            isurf = SurfIdx(ir, j, 2)
            
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
  
  PhiAngIn(:,PhiAngOutSvIdx) = phiobdPolar(:)
END DO
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

USE TIMER
USE ALLOCS
USE OMP_LIB
USE PARAM,   ONLY : TRUE, FALSE, ZERO, ONE
USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type, TrackingDat_Type, Pin_Type, Cell_Type, AziAngleInfo_Type, PolarAngle_Type
USE Moc_Mod, ONLY : nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay, Expa, Expb, TrackingDat, ApproxExp, wtang, SrcAng1, SrcAng2, comp, mwt
USE geom,    ONLY : nbd
USE cntl,    ONLY : nTracerCntl
USE PE_MOD,  ONLY : PE
USE IOUTIL,  ONLY : terminate

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:)     :: phis, xst, src
REAL, POINTER, DIMENSION(:,:)   :: PhiAngIn, phim, srcm
REAL, POINTER, DIMENSION(:,:,:) :: jout

INTEGER :: iz, ScatOd
LOGICAL :: ljout
INTEGER, OPTIONAL :: FastMocLv

LOGICAL, SAVE :: lfirst
DATA lfirst /TRUE/

TYPE (AziAngleInfo_Type), POINTER, DIMENSION(:) :: AziAng
TYPE (PolarAngle_Type),   POINTER, DIMENSION(:) :: PolarAng
TYPE (Cell_Type),         POINTER, DIMENSION(:) :: Cell
TYPE (Pin_Type),          POINTER, DIMENSION(:) :: Pin

INTEGER :: nAziAng, nPolarAng, nPhiAngSv, nRotray, nFsr, nxy, nThread, od, iod, iRotRay, tid, FsrIdxSt, icel, ireg, iazi, ipol, i, j, l, k, m
REAL :: wttmp, wtcos, wtpolar, wtsin2, tempsrc, ONETHREE, ONEFIVE, ONESEVEN

REAL, ALLOCATABLE :: wtp(:)

LOGICAL, SAVE :: lScatBd

INTEGER, SAVE :: OmpTemp, OmpAng, nOmpAng

INTEGER, POINTER, SAVE, DIMENSION(:)   :: nAziRotRay, nAziRotRay0, nSegRotRay, OmpRayBeg, OmpRayEnd, OmpRayList
INTEGER, POINTER, SAVE, DIMENSION(:,:) :: OmpRayBegBd, OmpRayEndBd, OmpMap, OmpRayNum
! ----------------------------------------------------

AziAng   => RayInfo%AziAngle
PolarAng => RayInfo%PolarAngle
nAziAng   = RayInfo%nAziAngle
nPolarAng = RayInfo%nPolarAngle
nPhiAngSv = RayInfo%nPhiAngSv
nRotRay   = RayInfo%nRotRay

nFsr = CoreInfo%nCoreFsr
nxy  = CoreInfo%nxy

nThread = PE%nThread
! ----------------------------------------------------
OD = 2

IF (ScatOd .EQ. 2) OD = 5
IF (ScatOd .EQ. 3) OD = 9

ALLOCATE(wtp(1:OD))
! ----------------------------------------------------
IF (lfirst) THEN
  lScatBd = nTracerCntl%lScatBd
  
  IF (lScatBd) THEN
    CALL dmalloc(nAziRotRay, nAziAng/2)
    CALL dmalloc(nAziRotRay0, nThread)
    CALL dmalloc(nSegRotRay, nAziAng/2)
        
    DO i = 1, nRotRay
      j = RayInfo%RotRay(i)%Ang1
      
      nAziRotRay(j) = nAziRotRay(j) + 1
      nSegRotRay(j) = nSegRotRay(j) + RayInfo%RotRay(i)%nSeg
    END DO
    
    CALL dmalloc0(OmpRayBeg, 0, nThread)
    CALL dmalloc0(OmpRayEnd, 0, nThread)
    
    k = nAziAng/2/nThread
    nAziRotRay0(:) = 0
    OmpRayBeg(0) = 1
    OmpTemp = 0
    OmpRayEnd(0) = 0
    l = 0
    
    DO i = 1, nThread
      OmpRayBeg(i) = OmpRayBeg(i-1) + OmpTemp
      
      DO j = 1, k
        l = l + 1
        
        nAziRotRay0(i) = nAziRotRay0(i) + nAziRotRay(l)
      END DO
      
      OmpTemp = nAziRotRay0(i)
      
      OmpRayEnd(i) = OmpRayEnd(i-1) + OmpTemp
    END DO
    
    ! Parallel Ang
    DO i = 1, nRotRay
      RayInfo%RotRay(i)%OmpAng1 = MOD(RayInfo%RotRay(i)%Ang1, k)
      
      IF (RayInfo%RotRay(i)%OmpAng1 .EQ. 0) RayInfo%RotRay(i)%OmpAng1 = k
      
      RayInfo%RotRay(i)%OmpAng2 = RayInfo%RotRay(i)%OmpAng1 + k
    END DO
    
    DO i = 1, nRotRay
      DO j = 1, RayInfo%RotRay(i)%nRay, 2
        RayInfo%RotRay(i)%OmpRayIdx(j)   = RayInfo%RotRay(i)%OmpAng1
        RayInfo%RotRay(i)%OmpRayIdx(j+1) = RayInfo%RotRay(i)%OmpAng2
      END DO
    END DO
    
    CALL dmalloc(OmpMap, nThread, k*2)
    
    DO i = 1, nThread
      DO j = 1, k
        OmpTemp = OmpTemp + 1
        OmpMap(i,j) = OmpTemp
        OmpMap(i,j+k) = nAziAng - OmpTemp + 1
      END DO
    END DO
    
    nOmpAng = nAziAng/nThread
  ELSE
    CALL dmalloc(nAziRotRay, nAziAng)
    CALL dmalloc(nAziRotRay0, 2*nThread)
    CALL dmalloc(nSegRotRay, nAziAng)
        
    DO i = 1, nRotRay
      j = RayInfo%RotRay(i)%Ang1
      
      nAziRotRay(j) = nAziRotRay(j) + 1
      nSegRotRay(j) = nSegRotRay(j) + RayInfo%RotRay(i)%nSeg
    END DO
    
    CALL dmalloc0(OmpRayBeg, 0, 2*nThread)
    CALL dmalloc0(OmpRayEnd, 0, 2*nThread)
    
    CALL dmalloc0(OmpRayBegBd, 1, 2, 0, nThread)
    CALL dmalloc0(OmpRayEndBd, 1, 2, 0 ,nThread)
    
    k = nAziAng/nThread/2
    nAziRotRay0(:) = 0
    OmpRayBegBd(1,0) = 1
    OmpTemp = 0
    OmpRayEndBd(1,0) = 0
    l = 0
    
    DO i = 1, nThread
      OmpRayBegBd(1,i) = OmpRayBegBd(1,i-1) + OmpTemp
      
      DO j = 1, k
        l = l + 1
        
        nAziRotRay0(i) = nAziRotRay0(i) + nAziRotRay(l)
      END DO
      
      OmpTemp = nAziRotRay0(i)
      OmpRayEndBd(1,i) = OmpRayEndBd(1,i-1) + OmpTemp
    END DO
    
    nAziRotRay0(:) = 0
    l = nAziAng + 1
    OmpRayBegBd(2,0) = nRotRay + 1
    OmpTemp = 0
    OmpRayEndBd(2,0) = nRotRay
    
    DO i = 1, nThread
      OmpRayEndBd(2,i) = OmpRayEndBd(2,i-1) - OmpTemp
      
      DO j = 1, k
        l = l - 1
        
        nAziRotRay0(i) = nAziRotRay0(i) + nAziRotRay(l)
      END DO
      OmpTemp = nAziRotRay0(i)
      OmpRayBegBd(2,i) = OmpRayBegBd(2,i-1) - OmpTemp
    END DO
    
    ! Parallel Ang
    DO i = 1, nRotRay
      IF (RayInfo%RotRay(i)%nRay .GT. 1) THEN
        IF (RayInfo%RotRay(i)%Ang1 .LE. nAziAng/2) THEN
          RayInfo%RotRay(i)%OmpAng1 = MOD(RayInfo%RotRay(i)%Ang1, k)
          
          IF (RayInfo%RotRay(i)%OmpAng1 .EQ. 0) RayInfo%RotRay(i)%OmpAng1 = k
          
          RayInfo%RotRay(i)%OmpAng2 = RayInfo%RotRay(i)%OmpAng1 + k
        ELSE IF (RayInfo%RotRay(i)%Ang1 .GT. nAziAng/2) THEN
          RayInfo%RotRay(i)%OmpAng1 = MOD(nAziAng+1-RayInfo%RotRay(i)%Ang1, k) + k
          
          IF (RayInfo%RotRay(i)%OmpAng1 .EQ. k) RayInfo%RotRay(i)%OmpAng1 = RayInfo%RotRay(i)%OmpAng1 + k
          
          RayInfo%RotRay(i)%OmpAng2 = RayInfo%RotRay(i)%OmpAng1 - k
        END IF
      ELSE IF (RayInfo%RotRay(i)%nRay .EQ. 1) THEN
        IF (RayInfo%RotRay(i)%Ang1 .LE. nAziAng/2) THEN
          RayInfo%RotRay(i)%OmpAng1 = MOD(RayInfo%RotRay(i)%Ang1, k)
          
          IF (RayInfo%RotRay(i)%OmpAng1 .EQ. 0) RayInfo%RotRay(i)%OmpAng1 = k
          
          RayInfo%RotRay(i)%OmpAng2 = 0
        ELSE IF (RayInfo%RotRay(i)%Ang1 .GT. nAziAng/2) THEN
          RayInfo%RotRay(i)%OmpAng1 = MOD(nAziAng+1-RayInfo%RotRay(i)%Ang1, k) + k
          
          IF (RayInfo%RotRay(i)%OmpAng1 .EQ. k) RayInfo%RotRay(i)%OmpAng1 = RayInfo%RotRay(i)%OmpAng1 + k
          
          RayInfo%RotRay(i)%OmpAng2 = 0
        END IF
      END IF
    END DO
    
    DO i = 1, nRotRay
      IF (RayInfo%RotRay(i)%nRay .GT. 1) THEN
        IF (MOD(RayInfo%RotRay(i)%nRay, 2) .EQ. 0) THEN
          DO j = 1, RayInfo%RotRay(i)%nRay, 2
            RayInfo%RotRay(i)%OmpRayIdx(j)   = RayInfo%RotRay(i)%OmpAng1
            RayInfo%RotRay(i)%OmpRayIdx(j+1) = RayInfo%RotRay(i)%OmpAng2
          END DO
        ELSE
          DO j = 1, RayInfo%RotRay(i)%nRay-1, 2
            RayInfo%RotRay(i)%OmpRayIdx(j)   = RayInfo%RotRay(i)%OmpAng1
            RayInfo%RotRay(i)%OmpRayIdx(j+1) = RayInfo%RotRay(i)%OmpAng2
          END DO
          RayInfo%RotRay(i)%OmpRayIdx(RayInfo%RotRay(i)%nRay) = RayInfo%RotRay(i)%OmpAng1
        END IF
      ELSE
        RayInfo%RotRay(i)%OmpRayIdx(1) = RayInfo%RotRay(i)%OmpAng1
      END IF
    END DO
    
    CALL dmalloc(OmpMap, nThread, k*2)
    
    DO i = 1, nThread
      DO j = 1, k
        OmpTemp = OmpTemp + 1
        
        OmpMap(i,j)   = OmpTemp
        OmpMap(i,j+k) = nAziAng - OmpTemp + 1
      END DO
    END DO
    
    CALL dmalloc(OmpRayList, nThread)
    CALL dmalloc(OmpRayNum, nThread, nRotRay)
    
    nOmpAng = nAziAng/nThread
    OmpRayList(:) = 0; OmpRayNum(:,:) = 0
    
    DO i = 1, nThread
      DO j = 1, 2
        OmpRayList(i) = OmpRayList(i) + OmpRayEndBd(j,i) - OmpRayBegBd(j,i) + 1
      END DO
    END DO
    
    DO i = 1, nThread
      m = 0
      
      DO j = 1, 2
        l = 0
        
        DO k = 1, OmpRayEndBd(j,i) - OmpRayBegBd(j,i) + 1
          m = m + 1
          
          OmpRayNum(i,m) = OmpRayBegBd(j,i) + l
          
          l = l + 1
        END DO
      END DO
    END DO
  END IF
  
  IF (MOD(nAziANg/2, nThread) .NE. 0) CALL terminate('WRONG_MOC_TRD')
END IF
! ----------------------------------------------------
IF (lfirst) THEN
  
  lFirst = FALSE
  
  CALL ApproxExp(RayInfo%PolarAngle, nPolarAng)
  
  DO tid = 1, nThread
     IF (TrackingDat(tid)%lAllocP1) CYCLE
     
     CALL dmalloc(TrackingDat(tid)%phi1a, nPolarAng, nFsr, nOmpAng)
     CALL dmalloc(TrackingDat(tid)%phi2a, nPolarAng, nFsr, nOmpAng)
     
     TrackingDat(tid)%lAllocP1 = TRUE
     
     IF (TrackingDat(tid)%lAlloc) CYCLE
     
     CALL dmalloc(TrackingDat(tid)%FsrIdx,         nMaxRaySeg, nMaxCoreRay)
     CALL dmalloc(TrackingDat(tid)%ExpAppIdx,      nMaxRaySeg, nMaxCoreRay)
     CALL dmalloc(TrackingDat(tid)%OptLenList,     nMaxRaySeg, nMaxCoreRay)
     CALL dmalloc(TrackingDat(tid)%ExpAppPolar,    nPolarAng, nMaxRaySeg, nMaxCoreRay)
     CALL dmalloc(TrackingDat(tid)%PhiAngOutPolar, nPolarAng, nMaxRaySeg+2)
     CALL dmalloc(TrackingDat(tid)%Jout,           3, nbd, nxy)
     
     TrackingDat(tid)%lAlloc = TRUE
  END DO
  
  CALL dmalloc(SrcAng1, nPolarAng, nFsr, nAziAng)
  CALL dmalloc(SrcAng2, nPolarAng, nFsr, nAziAng)
  CALL dmalloc(Comp,Od, nPolarAng, nAziAng)
  CALL dmalloc(mwt, Od, nPolarAng, nAziAng)
  CALL dmalloc(wtang,   nPolarAng, nAziAng)
  
  DO ipol = 1, nPolarAng
    wttmp = PolarAng(ipol)%weight * PolarAng(ipol)%sinv
    
    DO iazi = 1, nAziAng
      wtang(ipol, iazi) = wttmp  * AziAng(iazi)%weight * AziAng(iazi)%del
    END DO
  END DO
  
  DO ipol = 1, nPolarAng
    wttmp = PolarAng(ipol)%sinv
    
    DO iazi = 1, nAziAng
      comp(1, ipol, iazi) = wttmp * AziAng(iazi)%cosv
      comp(2, ipol, iazi) = wttmp * AziAng(iazi)%sinv
      
      mwt(1:2, ipol, iazi) = comp(1:2, ipol, iazi) * wtang(ipol, iazi)
    END DO
  END DO
  
  IF (ScatOd .GE. 2) THEN
    DO ipol = 1, nPolarAng
      wttmp   = PolarAng(ipol)%sinv
      wtsin2  = PolarAng(ipol)%sinv * PolarAng(ipol)%sinv
      wtcos   = PolarAng(ipol)%cosv
      wtpolar = 1.5_8 * PolarAng(ipol)%cosv * PolarAng(ipol)%cosv - 0.5_8
      
      DO iazi = 1, nAziAng
        Comp(3, ipol, iazi) = wtpolar
        Comp(4, ipol, iazi) = wtsin2 * (1._8-2._8*AziAng(iazi)%sinv*AziAng(iazi)%sinv)
        Comp(5, ipol, iazi) = wtsin2 * (2._8 * AziAng(iazi)%sinv * AziAng(iazi)%cosv)

        mwt(3,   ipol, iazi) = comp(3, ipol, iazi) *  wtang(ipol, iazi)
        mwt(4:5, ipol, iazi) = 0.75_8 * comp(4:5, ipol, iazi) *  wtang(ipol, iazi)
      END DO
    END DO
  END IF
  
  IF (ScatOd .EQ. 3) THEN
    DO ipol = 1, nPolarAng
      wttmp = PolarAng(ipol)%sinv
      DO iazi = 1, nAziAng
        Comp(6, ipol, iazi) = (5._8 * PolarAng(ipol)%cosv * PolarAng(ipol)%cosv - 1._8) * wttmp * AziAng(iazi)%cosv
        Comp(7, ipol, iazi) = (5._8 * PolarAng(ipol)%cosv * PolarAng(ipol)%cosv - 1._8) * wttmp * AziAng(iazi)%sinv
        Comp(8, ipol, iazi) = (wttmp ** 3._8) * (4._8 * (AziAng(iazi)%cosv ** 3._8) - 3._8 * AziAng(iazi)%cosv)
        Comp(9, ipol, iazi) = (wttmp ** 3._8) * (- 4._8 * (AziAng(iazi)%sinv ** 3._8) + 3._8 * AziAng(iazi)%sinv)

        mwt(6:7, ipol, iazi) = 0.375_8 * comp(6:7, ipol, iazi) * wtang(ipol, iazi)
        mwt(8:9, ipol, iazi) = 0.625_8 * comp(8:9, ipol, iazi) * wtang(ipol, iazi)
      END DO
    END DO
  END IF
END IF
! ----------------------------------------------------
!$  call omp_set_dynamic(FALSE)
!$  call omp_set_num_threads(nThread)

DO tid = 1, nThread
  TrackingDat(tid)%Expa     => Expa
  TrackingDat(tid)%Expb     => Expb
  TrackingDat(tid)%PhiAngIn => PhiAngIn
  TrackingDat(tid)%src      => src
  TrackingDat(tid)%xst      => xst
  TrackingDat(tid)%wtang    => WtAng
  TrackingDat(tid)%SrcAng1  => SrcAng1
  TrackingDat(tid)%SrcAng2  => SrcAng2
END DO

phis = ZERO
phim = ZERO
IF (ljout) Jout = ZERO
! ----------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(irotray, tid, iazi, ipol, i, j, tempsrc, k, l)
tid = 1

!$  tid = omp_get_thread_num()+1

IF (ScatOd .EQ. 1) THEN
  DO iazi = 1, nAziAng
    DO i = PE%myOmpFsrBeg(tid), PE%myOmpFsrEnd(tid)
      DO ipol = 1, nPolarAng
        SrcAng1(ipol, i, iazi) = src(i)
        SrcAng2(ipol, i, iazi) = src(i)
        
        tempsrc = comp(1, ipol, iazi) * srcm(1, i) + comp(2, ipol, iazi) * srcm(2, i)
        
        SrcAng1(ipol, i, iazi) = SrcAng1(ipol, i, iazi) + tempsrc
        SrcAng2(ipol, i, iazi) = SrcAng2(ipol, i, iazi) - tempsrc
       END DO
    END DO
  END DO
ELSE IF (ScatOd .EQ. 2) THEN
  DO iazi = 1, nAziAng
    DO i = PE%myOmpFsrBeg(tid), PE%myOmpFsrEnd(tid)
      DO ipol = 1, nPolarAng
        SrcAng1(ipol, i, iazi) = src(i)
        SrcAng2(ipol, i, iazi) = src(i)
        
        tempsrc = comp(1, ipol, iazi) * srcm(1, i) + comp(2, ipol, iazi) * srcm(2, i)
        
        SrcAng1(ipol, i, iazi) = SrcAng1(ipol, i, iazi) + tempsrc
        SrcAng2(ipol, i, iazi) = SrcAng2(ipol, i, iazi) - tempsrc
        
        tempsrc = comp(3, ipol, iazi) * srcm(3, i) + comp(4, ipol, iazi) * srcm(4, i) + comp(5, ipol, iazi) * srcm(5, i)
        
        SrcAng1(ipol, i, iazi) = SrcAng1(ipol, i, iazi) +  tempsrc
        SrcAng2(ipol, i, iazi) = SrcAng2(ipol, i, iazi) +  tempsrc
      END DO
    END DO
  END DO
ELSE IF (ScatOd .EQ. 3) THEN
  DO iazi = 1, nAziAng
    DO i = PE%myOmpFsrBeg(tid), PE%myOmpFsrEnd(tid)
      DO ipol = 1, nPolarAng
        SrcAng1(ipol, i, iazi) = src(i)
        SrcAng2(ipol, i, iazi) = src(i)
        
        tempsrc = comp(1, ipol, iazi) * srcm(1, i) + comp(2, ipol, iazi) * srcm(2, i) + comp(6, ipol, iazi) * srcm(6, i) + comp(7, ipol, iazi) * srcm(7, i) + comp(8, ipol, iazi) * srcm(8, i) + comp(9, ipol, iazi) * srcm(9, i)
        
        SrcAng1(ipol, i, iazi) = SrcAng1(ipol, i, iazi) + tempsrc
        SrcAng2(ipol, i, iazi) = SrcAng2(ipol, i, iazi) - tempsrc
        
        tempsrc = comp(3, ipol, iazi) * srcm(3, i) + comp(4, ipol, iazi) * srcm(4, i) + comp(5, ipol, iazi) * srcm(5, i)
        
        SrcAng1(ipol, i, iazi) = SrcAng1(ipol, i, iazi) +  tempsrc
        SrcAng2(ipol, i, iazi) = SrcAng2(ipol, i, iazi) +  tempsrc
      END DO
    END DO
  END DO
END IF
!$OMP END PARALLEL
! ----------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(irotray, tid, i, j, k, l, iazi, ipol, OmpAng)
tid = 1
!$  tid = omp_get_thread_num()+1
j = tid

DO iazi = 1, nOmpAng
  TrackingDat(j)%phi1a(:,:,iazi) = ZERO
  TrackingDat(j)%phi2a(:,:,iazi) = ZERO
END DO

DO i = 1, nxy
  TrackingDat(j)%jout(:, :, i) = zero
END DO
!$OMP BARRIER
IF (lScatBd) THEN
  DO i = OmpRayBeg(tid), OmpRayEnd(tid)
    irotray = i
    
    CALL RecTrackRotRayP1GM_AFSS(RayInfo, CoreInfo, TrackingDat(tid), lJout, irotray, iz, FastMocLv)
  END DO
ELSE
  DO i = 1, 2
    DO j = OmpRayBegBd(i, tid), OmpRayEndBd(i, tid)
      irotray = j
      
      CALL RecTrackRotRayP1GM_AFSS(RayInfo, CoreInfo, TrackingDat(tid), lJout, irotray, iz, FastMocLv)
    END DO
  END DO
END IF
!$OMP BARRIER
IF (ScatOd .EQ. 1) THEN
  DO j = 1, nThread
    DO iazi = 1, nOmpAng
      DO i = PE%myOmpFsrBeg(tid), PE%myOmpFsrEnd(tid)
        DO ipol = 1, nPolarAng
          OmpAng = OmpMap(j, iazi)
          
          phis(i) = phis(i) + wtang(ipol, OmpAng) * (TrackingDat(j)%phi1a(ipol, i, iazi) + TrackingDat(j)%phi2a(ipol, i, iazi))
          
          phim(1:2, i) = phim(1:2, i) + mwt(1:2, ipol, OmpAng) * (TrackingDat(j)%phi1a(ipol, i, iazi) - TrackingDat(j)%phi2a(ipol, i, iazi))
        END DO
      END DO
    END DO
  END DO
ELSE IF (ScatOd .EQ. 2) THEN
  DO j = 1, nThread
    DO iazi = 1, nOmpAng
      DO i = PE%myOmpFsrBeg(tid), PE%myOmpFsrEnd(tid)
        DO ipol = 1, nPolarAng
          OmpAng = OmpMap(j, iazi)
          
          phis(i) = phis(i) + wtang(ipol, OmpAng) * (TrackingDat(j)%phi1a(ipol, i, iazi) + TrackingDat(j)%phi2a(ipol, i, iazi))
          
          phim(1:2, i) = phim(1:2, i) + mwt(1:2, ipol, OmpAng) * (TrackingDat(j)%phi1a(ipol, i, iazi) - TrackingDat(j)%phi2a(ipol, i, iazi))
          phim(3:5, i) = phim(3:5, i) + mwt(3:5, ipol, OmpAng) * (TrackingDat(j)%phi1a(ipol, i, iazi) + TrackingDat(j)%phi2a(ipol, i, iazi))
        END DO
      END DO
    END DO
  END DO
ELSE IF (ScatOd .EQ. 3) THEN
  DO j = 1, nThread
    DO iazi = 1, nOmpAng
      DO i = PE%myOmpFsrBeg(tid), PE%myOmpFsrEnd(tid)
        DO ipol = 1, nPolarAng
          OmpAng = OmpMap(j, iazi)
          
          phis(i) = phis(i) + wtang(ipol, OmpAng) * (TrackingDat(j)%phi1a(ipol, i, iazi) + TrackingDat(j)%phi2a(ipol, i, iazi))
          
          phim(1:2, i) = phim(1:2, i) + mwt(1:2, ipol, OmpAng) * (TrackingDat(j)%phi1a(ipol, i, iazi) - TrackingDat(j)%phi2a(ipol, i, iazi))
          phim(3:5, i) = phim(3:5, i) + mwt(3:5, ipol, OmpAng) * (TrackingDat(j)%phi1a(ipol, i, iazi) + TrackingDat(j)%phi2a(ipol, i, iazi))
          phim(6:9, i) = phim(6:9, i) + mwt(6:9, ipol, OmpAng) * (TrackingDat(j)%phi1a(ipol, i, iazi) - TrackingDat(j)%phi2a(ipol, i, iazi))
        END DO
      END DO
    END DO
  END DO
END IF

IF (ljout) THEN
  DO j = 1, nThread
    DO i = PE%myOmpNxyBeg(tid), PE%myOmpNxyEnd(tid)
      jout(:, :, i) = jout(:, :, i) + TrackingDat(j)%jout(:, :, i)
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

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, j, FsrIdxSt, icel, ireg, wttmp, tid)
tid = 1

!$  tid = omp_get_thread_num()+1

IF (ScatOd .EQ. 1) THEN
  DO l = PE%myOmpNxyBeg(tid), PE%myOmpNxyEnd(tid)
    FsrIdxSt = Pin(l)%FsrIdxSt
    icel     = Pin(l)%Cell(iz);
    
    DO j = 1, Cell(icel)%nFsr
      ireg = FsrIdxSt + j - 1
      
      wttmp = 1._8/(xst(ireg)*Cell(icel)%vol(j))
      
      phis(ireg) = phis(ireg) * wttmp + src(ireg)
      
      phim(1:2, ireg) = phim(1:2, ireg) * wttmp + srcm(1:2, ireg) * ONETHREE
    END DO
  END DO
ELSE IF (ScatOd .EQ. 2) THEN
  DO l = PE%myOmpNxyBeg(tid), PE%myOmpNxyEnd(tid)
    FsrIdxSt = Pin(l)%FsrIdxSt
    icel     = Pin(l)%Cell(iz);
    
    DO j = 1, Cell(icel)%nFsr
      ireg = FsrIdxSt + j - 1
      
      wttmp = 1._8/(xst(ireg)*Cell(icel)%vol(j))
      
      phis(ireg) = phis(ireg) * wttmp + src(ireg)
      
      phim(1:2, ireg) = phim(1:2, ireg) * wttmp + srcm(1:2, ireg) * ONETHREE
      phim(3:5, ireg) = phim(3:5, ireg) * wttmp + srcm(3:5, ireg) * ONEFIVE
    END DO
  END DO
ELSE IF (ScatOd .EQ. 3) THEN
  DO l = PE%myOmpNxyBeg(tid), PE%myOmpNxyEnd(tid)
    FsrIdxSt = Pin(l)%FsrIdxSt
    icel     = Pin(l)%Cell(iz)
    
    DO j = 1, Cell(icel)%nFsr
      ireg = FsrIdxSt + j - 1
      
      wttmp = 1._8/(xst(ireg)*Cell(icel)%vol(j))
      
      phis(ireg) = phis(ireg) * wttmp + src(ireg)
      
      phim(1:2, ireg) = phim(1:2, ireg) * wttmp + srcm(1:2, ireg) * ONETHREE
      phim(3:5, ireg) = phim(3:5, ireg) * wttmp + srcm(3:5, ireg) * ONEFIVE
      phim(6:9, ireg) = phim(6:9, ireg) * wttmp + srcm(6:9, ireg) * ONESEVEN
    END DO
  END DO
END IF
!$OMP END PARALLEL
! ----------------------------------------------------
DEALLOCATE(wtp)

NULLIFY (AziAng)
NULLIFY (PolarAng)
NULLIFY (Cell)
NULLIFY (Pin)
! ----------------------------------------------------

END SUBROUTINE RayTraceP1GM_AFSS
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayP1GM_AFSS(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, FastMocLv)

USE PARAM,   ONLY : TRUE, FALSE, ZERO
USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type, Pin_Type, Asy_Type, AsyInfo_Type, PinInfo_Type, Cell_Type, AziAngleInfo_Type, PolarAngle_Type, ModRayInfo_type,  AsyRayInfo_type,  &
                    CoreRayInfo_Type, RotRayInfo_Type, CellRayInfo_type, FastCoreRayDat_Type, TrackingDat_Type, FastRaySegDat_Type
USE Moc_Mod, ONLY : nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay

IMPLICIT NONE

TYPE (RayInfo_Type),     INTENT(INOUT) :: RayInfo
TYPE (CoreInfo_Type),    INTENT(INOUT) :: CoreInfo
TYPE (TrackingDat_Type), INTENT(INOUT) :: TrackingDat

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz, FastMocLv

TYPE (Pin_Type),          POINTER, DIMENSION(:) :: Pin
TYPE (Asy_Type),          POINTER, DIMENSION(:) :: Asy
TYPE (PinInfo_Type),      POINTER, DIMENSION(:) :: PinInfo
TYPE (Cell_Type),         POINTER, DIMENSION(:) :: Cell
TYPE (AziAngleInfo_Type), POINTER, DIMENSION(:) :: AziAng
TYPE (PolarAngle_Type),   POINTER, DIMENSION(:) :: PolarAng
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
REAL, POINTER, DIMENSION(:,:,:) :: ExpAppPolar, jout, wtsurf, phi1a, phi2a, SrcAng1, SrcAng2

INTEGER :: mp(2)
INTEGER :: iazi, ipol, iCoreRay, iasyray, iceray, irayseg, irot, itype, idir, nRotRay, nCoreRay, nAsyRay, nPinRay, nRaySeg, FsrIdxSt, nPolarAng, nAziAng, nPhiAngSv
INTEGER :: ipin, icel, ibcel, iasy, ireg, isurf, irsegidx, icellrayidx, PhiAnginSvIdx, PhiAngOutSvIdx, nFsr, nxy, i, j, k, l, m, jbeg, jend, jinc, ir, ir1, iOmpAzi

INTEGER, DIMENSION(nMaxCellRay, nMaxCoreRay, 2) :: CellRayIdxSt, SurfIdx
INTEGER, DIMENSION(nMaxCellRay, nMaxCoreRay) :: PinIdx
INTEGER, DIMENSION(nMaxCoreRay) :: nTotRaySeg, nTotCellRay

REAL :: tau, phiobd, phid, wt
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
        
        ibcel    = Cell(icel)%BaseCellStr
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
    
    nTotRaySeg (j) = FastRay%nTotRaySeg (j)
    nTotCellRay(j) = FastRay%nTotCellRay(j)
    
    DO l = 1, FastRay%nTotCellRay(j)
      PinIdx      (l, j) = FastRay%PinIdx         (l, j)
      CellRayIdxSt(l, j, 1) = FastRay%CellRayIdxSt(l, j, 1)
      CellRayIdxSt(l, j, 2) = FastRay%CellRayIdxSt(l, j, 2)
      SurfIdx     (l, j, 1) = FastRay%SurfIdx     (l, j, 1)
      SurfIdx     (l, j, 2) = FastRay%SurfIdx     (l, j, 2)
      
      ipin = PinIdx(l, j)
      
      icel     = Pin(ipin)%Cell(iz)
      FsrIdxSt = Pin(ipin)%FsrIdxSt
      
      DO k = FastRay%Ray1DIdx(1, l, j), FastRay%Ray1DIdx(2, l, j)
        irsegidx = irsegidx + 1
        ireg     = FsrIdxSt + CellRay1D%LocalFsrIdx(K) - 1
        
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
      PinIdx      (l, j)    = FastRay%PinIdx      (l, j)
      CellRayIdxSt(l, j, 1) = FastRay%CellRayIdxSt(l, j, 1)
      CellRayIdxSt(l, j, 2) = FastRay%CellRayIdxSt(l, j, 2)
      SurfIdx     (l, j, 1) = FastRay%SurfIdx     (l, j, 1)
      SurfIdx     (l, j, 2) = FastRay%SurfIdx     (l, j, 2)
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
    DO ipol = 1, nPolarANg
      ExpAppPolar(ipol, l, j) = expa(ExpAppIdx(l, j), ipol) * optlenlist(l, j) + expb(ExpAppIdx(l, j), ipol)
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
  
  DO j = jbeg, jend, jinc
    iazi    = CoreRay(RotRay(irotray)%RayIdx(j))%iang
    iompazi = RotRay(irotray)%OmpRayIdx(j)
    
    idir = RotRay(i)%DIR(j)
    IF (irot .eq. 2) idir = mp(idir) ! Reverse the sweep direction
    
    nRaySeg = nTotRaySeg(j)

    IF (idir .eq. 1) THEN ! Forward Sweep
      PhiAngOutPolar(:, 1) = phiobdPolar(:)
      
      DO ir = 1, nRaySeg
        ireg = FsrIdx(ir, j)
        
        DO ipol = 1, nPolarAng
          phid = (PhiAngOutPolar(ipol, ir) - SrcAng1(ipol, ireg, iazi)) * ExpAppPolar(ipol, ir, j)
          
          PhiAngOutPolar(ipol, ir+1) = PhiAngOutPolar(ipol, ir) - phid
          
          phi1a(ipol, ireg, iompazi) =  phi1a(ipol, ireg, iompazi) + phid
        END DO
      END DO
      
      phiobdPolar(:) = PhiAngOutPolar(:, nRaySeg+1)
      
      ! Surface
      IF (ljout) THEN
        DO ir = 1, nTotCellRay(j)
          DO ipol = 1, nPolarANg
            wt = wtang(ipol, iazi)
            
            icel  = PinIdx(ir, j)
            isurf = SurfIdx(ir, j, 1)
            
            Jout(2, isurf, icel) = Jout(2, isurf, icel) + wt * PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 1)+1)
            
            isurf = SurfIdx(ir, j, 2)
            
            Jout(1, isurf, icel) = Jout(1, isurf, icel) + wt * PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 2))
          END DO
        END DO
      END IF
    ELSE
      PhiAngOutPolar(:, nRaySeg+2) = phiobdPolar(:)
      
      ir = nRaySeg + 1
      
      DO ir1 = 1, nRaySeg
        ir   = ir - 1
        ireg = FsrIdx(ir, j)
        
        DO ipol = 1, nPolarAng
          phid = (PhiAngOutPolar(ipol, ir + 2) - SrcAng2(ipol, ireg, iazi)) * ExpAppPolar(ipol, ir, j)
          
          PhiAngOutPolar(ipol, ir+1) = PhiAngOutPolar(ipol, ir + 2) - phid
          
          phi2a(ipol, ireg, iompazi) =  phi2a(ipol, ireg, iompazi) + phid
        END DO
      END DO
      
      phiobdPolar(:) = PhiAngOutPolar(:, 2)
      
      IF (lJout) THEN
        DO ir = 1, nTotCellRay(j)
          DO ipol = 1, nPolarAng
            wt = wtang(ipol, iazi)
            
            icel = PinIdx(ir, j)
            isurf = SurfIdx(ir, j, 2)
            
            Jout(2, isurf, icel) = Jout(2, isurf, icel) + wt * PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 2)+1)
            
            isurf = SurfIdx(ir, j, 1)
            
            Jout(1, isurf, icel) = Jout(1, isurf, icel) + wt * PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 1)+2)
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