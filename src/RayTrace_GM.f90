#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceGM_OMP(RayInfo, CoreInfo, phis, PhiAngIn, xst, src, jout, iz, ljout, FastMocLv)

USE TIMER
USE ALLOCS
USE OMP_LIB
USE PARAM,   ONLY : TRUE, FALSE, ZERO
USE TYPEDEF, ONLY : RayInfo_Type, CoreInfo_type, Pin_Type, Cell_Type, AziAngleInfo_Type, PolarAngle_Type
USE Moc_Mod, ONLY : nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay, Expa, Expb, ApproxExp, TrackingDat, wtang
USE geom,    ONLY : nbd
USE PE_MOD,  ONLY : PE
USE CNTL,    ONLY : nTracerCntl

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
LOGICAL, SAVE :: lfirst
DATA lfirst /TRUE/

TYPE (AziAngleInfo_Type), POINTER, DIMENSION(:) :: AziAng
TYPE (PolarAngle_Type),   POINTER, DIMENSION(:) :: PolarAng
TYPE (Cell_Type),         POINTER, DIMENSION(:) :: Cell
TYPE (Pin_Type),          POINTER, DIMENSION(:) :: Pin

INTEGER :: nAziAng, nPolarAng, nPhiAngSv, nRotray, nFsr, nxy, nThread
INTEGER :: ithr, FsrIdxSt, icel, ireg, iazi, ipol, iray, ifsr, jfsr, ipin, krot

REAL :: wttmp
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
  lFirst = FALSE
  
  CALL ApproxExp(RayInfo%PolarAngle, nPolarAng)
  
  DO ithr = 1, nThread
    IF (TrackingDat(ithr)%lAlloc) CYCLE
    
    CALL dmalloc(TrackingDat(ithr)%FsrIdx,         nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%ExpAppIdx,      nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%OptLenList,     nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%ExpAppPolar,    nPolarAng,  nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%PhiAngOutPolar, nPolarAng,  nMaxRaySeg + 2)
    CALL dmalloc(TrackingDat(ithr)%Phis,           nFsr)
    CALL dmalloc(TrackingDat(ithr)%Jout,           3, nbd, nxy)
    
    TrackingDat(ithr)%lAlloc = TRUE
  END DO
  
  CALL dmalloc(wtang, nPolarAng, nAziAng)
  
  DO ipol = 1, nPolarAng
    wttmp = PolarAng(ipol)%weight * PolarAng(ipol)%sinv
    
    DO iazi = 1, nAziAng
      wtang(ipol, iazi) = wttmp * AziAng(iazi)%weight * AziAng(iazi)%del
    END DO
  END DO
END IF
! ----------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ithr, iray)
ithr = omp_get_thread_num() + 1

TrackingDat(ithr)%phis = ZERO
TrackingDat(ithr)%jout = ZERO

TrackingDat(ithr)%Expa     => Expa
TrackingDat(ithr)%Expb     => Expb
TrackingDat(ithr)%src      => src
TrackingDat(ithr)%xst      => xst
TrackingDat(ithr)%PhiAngIn => PhiAngIn
TrackingDat(ithr)%wtang    => wtang

DO krot = 1, 2
  IF (nTracerCntl%lHex) THEN
    !$OMP DO SCHEDULE(GUIDED)
    DO iray = 1, nRotRay
      CALL HexTrackRotRayGM_OMP(RayInfo, CoreInfo, TrackingDat(ithr), ljout, iray, iz, krot, FastMocLv)
    END DO
    !$OMP END DO NOWAIT
  ELSE
    !$OMP DO SCHEDULE(GUIDED)
    DO iray = 1, nRotRay
      CALL RecTrackRotRayGM_OMP(RayInfo, CoreInfo, TrackingDat(ithr), ljout, iray, iz, krot, FastMocLv)
    END DO
    !$OMP END DO NOWAIT
  END IF
END DO
!$OMP END PARALLEL
! ----------------------------------------------------
phis = ZERO

IF (ljout) jout = ZERO

DO ithr = 1, nThread
  phis = phis + TrackingDat(ithr)%phis
  
  IF (.NOT. ljout) CYCLE
  
  jout = jout + TrackingDat(ithr)%jout
END DO
! ----------------------------------------------------
Cell => CoreInfo%CellInfo
Pin  => CoreInfo%Pin

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ipin, FsrIdxSt, icel, ifsr, jfsr)
!$OMP DO
DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt
  icel     = Pin(ipin)%Cell(iz)
  
  DO ifsr = 1, Cell(icel)%nFsr
    jfsr = FsrIdxSt + ifsr - 1
    
    phis(jfsr) = phis(jfsr) / xst(jfsr) / Cell(icel)%vol(ifsr) + src(jfsr)
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

END SUBROUTINE RayTraceGM_OMP
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceGM_AFSS(RayInfo, CoreInfo, phis, PhiAngIn, xst, src, jout, iz, ljout, FastMocLv)

USE TIMER
USE ALLOCS
USE OMP_LIB
USE PARAM,   ONLY : TRUE, FALSE, ZERO
USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type, TrackingDat_Type, Pin_Type, Cell_Type, AziAngleInfo_Type, PolarAngle_Type
USE Moc_Mod, ONLY : nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay, Expa, Expb, TrackingDat, ApproxExp, wtang, Phia1g
USE geom,    ONLY : nbd
USE cntl,    ONLY : nTracerCntl
USE PE_MOD,  ONLY : PE
USE ioutil,  ONLY : terminate

IMPLICIT NONE

TYPE (RayInfo_Type) :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:)     :: phis, xst, src
REAL, POINTER, DIMENSION(:,:)   :: PhiAngIn
REAL, POINTER, DIMENSION(:,:,:) :: jout

INTEGER :: iz
LOGICAL :: ljout
INTEGER, OPTIONAL :: FastMocLv
! ----------------------------------------------------
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

INTEGER, POINTER, SAVE, DIMENSION(:)   :: nAziRotRay, nAziRotRay0, nSegRotRay, OmpRayBeg, OmpRayEnd, OmpRayList
INTEGER, POINTER, SAVE, DIMENSION(:,:) :: OmpRayBegBd, OmpRayEndBd, OmpMap, OmpRayNum

INTEGER, SAVE :: OmpTemp, OmpAng, nOmpAng
! ----------------------------------------------------

! Basic
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
    
    k = nAziAng / 2 / nThread
    l = 0
    
    OmpTemp = 0
    
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
    
    OmpTemp = 0
    
    DO i = 1, nThread
      DO j = 1, k
        OmpTemp = OmpTemp + 1
        
        OmpMap(i, j)   = OmpTemp
        OmpMap(i, j+k) = nAziAng - OmpTemp + 1
      END DO
    END DO
      
    nOmpAng = nAziAng / nThread
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
    
    k = nAziAng / nThread / 2
    l = 0
    OmpRayBegBd(1, 0) = 1
    OmpTemp = 0
    
    DO i = 1, nThread
      OmpRayBegBd(1, i) = OmpRayBegBd(1, i-1) + OmpTemp
      
      DO j = 1, k
        l = l + 1
        
        nAziRotRay0(i) = nAziRotRay0(i) + nAziRotRay(l)
      END DO
      
      OmpTemp = nAziRotRay0(i)
      
      OmpRayEndBd(1, i) = OmpRayEndBd(1, i-1) + OmpTemp
    END DO
    
    l = nAziAng + 1
    
    OmpRayBegBd(2, 0) = nRotRay + 1
    OmpRayEndBd(2, 0) = nRotRay
    
    OmpTemp = 0
    
    DO i = 1, nThread
      OmpRayEndBd(2, i) = OmpRayEndBd(2, i-1) - OmpTemp
      
      DO j = 1, k
        l = l - 1
        
        nAziRotRay0(i) = nAziRotRay0(i) + nAziRotRay(l)
      END DO
      
      OmpTemp = nAziRotRay0(i)
      
      OmpRayBegBd(2, i) = OmpRayBegBd(2, i-1) - OmpTemp
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
      ELSE IF(RayInfo%RotRay(i)%nRay .EQ. 1) THEN
        IF(RayInfo%RotRay(i)%Ang1 .LE. nAziAng/2) THEN
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
    
    OmpTemp = 0
    
    DO i = 1, nThread
      DO j = 1, k
        OmpTemp = OmpTemp + 1
        
        OmpMap(i, j)   = OmpTemp
        OmpMap(i, j+k) = nAziAng - OmpTemp + 1
      END DO
    END DO
    
    CALL dmalloc(OmpRayList, nThread)
    CALL dmalloc(OmpRayNum,  nThread, nRotRay)
    
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
          
          OmpRayNum(i,m) = OmpRayBegBd(j, i) + l
          
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
    
    CALL Dmalloc(TrackingDat(ithr)%FsrIdx,     nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(ithr)%ExpAppIdx,  nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(ithr)%OptLenList, nMaxRaySeg, nMaxCoreRay)
    
    CALL Dmalloc(TrackingDat(ithr)%ExpAppPolar,    nPolarAng, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(ithr)%PhiAngOutPolar, nPolarAng, nMaxRaySeg+2)
    CALL Dmalloc(TrackingDat(ithr)%phi1a,          nPolarAng, nFsr, nOmpAng)
    CALL Dmalloc(TrackingDat(ithr)%phi2a,          nPolarAng, nFsr, nOmpAng)
    
    CALL Dmalloc(TrackingDat(ithr)%Jout, 3, nbd, nxy)
    
    TrackingDat(ithr)%Expa => Expa
    TrackingDat(ithr)%Expb => Expb
    
    TrackingDat(ithr)%lAlloc = TRUE
  END DO
  
  CALL dmalloc(wtang, nPolarAng, nAziAng)
  
  DO ipol = 1, nPolarAng
    wttmp = PolarAng(ipol)%weight * PolarAng(ipol)%sinv
    
    DO iazi = 1, nAziAng
      wtang(ipol, iazi) = wttmp  * AziAng(iazi)%weight * AziAng(iazi)%del
    END DO
  END DO
END IF
! ----------------------------------------------------
!$ call omp_set_dynamic(FALSE)
!$ call omp_set_num_threads(nThread)

DO ithr = 1, nThread
  TrackingDat(ithr)%PhiAngIn => PhiAngIn(:, :)
  TrackingDat(ithr)%src      => src
  TrackingDat(ithr)%xst      => xst
  TrackingDat(ithr)%wtang    => WtAng
END DO
! ----------------------------------------------------
phis = ZERO

IF (ljout) Jout = ZERO

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(irotray, ithr, i, j, k, l, iazi, ipol, OmpAng)
ithr = 1
!$ ithr = omp_get_thread_num()+1
j = ithr

TrackingDat(j)%phi1a = ZERO
TrackingDat(j)%phi2a = ZERO
TrackingDat(j)%jout  = ZERO
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
  icel     = Pin(l)%Cell(iz)
  
  DO j = 1, Cell(icel)%nFsr
    ireg = FsrIdxSt + j - 1
    
    phis(ireg) = phis(ireg) / xst(ireg) / Cell(icel)%vol(j) + src(ireg)
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
! ----------------------------------------------------

END SUBROUTINE RayTraceGM_AFSS
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayGM_OMP(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, krot, FastMocLv)

USE PARAM,   ONLY : FALSE, TRUE
USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type, Pin_Type, Asy_Type, AsyInfo_Type, PinInfo_Type, Cell_Type, AziAngleInfo_Type, PolarAngle_Type, ModRayInfo_type,  AsyRayInfo_type, &
                    CoreRayInfo_Type, RotRayInfo_Type, CellRayInfo_type, FastCoreRayDat_Type, TrackingDat_Type, FastRaySegDat_Type
USE Moc_Mod, ONLY : nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay

IMPLICIT NONE

TYPE(RayInfo_Type),     INTENT(INOUT) :: RayInfo
TYPE(CoreInfo_Type),    INTENT(INOUT) :: CoreInfo
TYPE(TrackingDat_Type), INTENT(INOUT) :: TrackingDat

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz, krot
INTEGER, INTENT(IN) :: FastMocLv
! ----------------------------------------------------
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
REAL, POINTER, DIMENSION(:,:,:) :: ExpAppPolar, jout
REAL, ALLOCATABLE, DIMENSION(:) :: phiobdPolar

INTEGER :: mp(2)
DATA mp /2, 1/

INTEGER :: iazi, ipol, iCoreRay, iasyray, iceray, irayseg, itype, idir, ipin, icel, iasy, ireg, isurf, irsegidx, icellrayidx, PhiAnginSvIdx, PhiAngOutSvIdx
INTEGER :: i, j, k, l, m, jbeg, jend, jinc, ir, ir1, ibcel
INTEGER :: nRotRay, nCoreRay, nAsyRay, nPinRay, nRaySeg, FsrIdxSt, nPolarAng, nAziAng, nPhiAngSv, nFsr, nxy

INTEGER, DIMENSION(nMaxCellRay, nMaxCoreRay, 2) :: CellRayIdxSt
INTEGER, DIMENSION(nMaxCellRay, nMaxCoreRay)    :: PinIdx
INTEGER, DIMENSION(nMaxCellRay, nMaxCoreRay, 2) :: SurfIdx

INTEGER, DIMENSION(nMaxCoreRay) :: nTotRaySeg
INTEGER, DIMENSION(nMaxCoreRay) :: nTotCellRay

REAL :: tau, phiobd, phid, wt
REAL :: wt2(4), wtang2(100, 100, 4)

LOGICAL :: lFast
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
! ----------------------------------------------------
IF (lJout)THEN
  DO ipol = 1, nPolarAng
    DO iazi = 1, nAziAng
      wtang2(ipol, iazi, 1) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv)
      wtang2(ipol, iazi, 3) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv)
      wtang2(ipol, iazi, 2) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv)
      wtang2(ipol, iazi, 4) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv)
    END DO
  END DO
END IF
! ----------------------------------------------------
i        = iRotRay
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
      
      IF (iasy .EQ. 0)  CYCLE
      nPinRay = AsyRay(iAsyRay)%nCellRay
      itype   = Asy(iasy)%PartialAsyFlag
      
      DO l = 1, nPinRay
        ipin   = AsyRay(iAsyRay)%PinIdx(l)
        iceray = AsyRay(iAsyRay)%PinRayIdx(l)
        ipin   = Asy(iAsy)%GlobalPinIdx(ipin)
        
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
          
          tau = - LenSeg(iRaySeg) * xst(ireg)   !
          tau = - CellRay%LenSeg(iRaySeg) * xst(ireg)   !
          
          irsegidx = irsegidx + 1
          
          FsrIdx    (irsegidx, j) = ireg
          OptLenList(irsegidx, j) = tau
          ExpAppIdx (irsegidx, j) = max(INT(tau), -40000)
          ExpAppIdx (irsegidx, j) = min(0, ExpAppIdx(irsegidx, j))
        END DO
        
        CellRayIdxSt(icellrayidx, j, 1) = irsegidx
        
        SurfIdx(icellRayIdx, j, 1) = AsyRay(iAsyRay)%PinRaySurf(2, l) !OutSurface
        SurfIdx(icellRayIdx, j, 2) = AsyRay(iAsyRay)%PinRaySurf(1, l) !Insurface
      END DO
    END DO
    
    nTotRaySeg (j) = irsegidx
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
       
       ipin     = PinIdx(l, j)
       icel     = Pin(ipin)%Cell(iz)
       FsrIdxSt = Pin(ipin)%FsrIdxSt
       
       DO k = FastRay%Ray1DIdx(1, l, j), FastRay%Ray1DIdx(2, l, j)
         irsegidx = irsegidx + 1
         ireg     = FsrIdxSt + CellRay1D%LocalFsrIdx(K) - 1
         
         FsrIdx(irsegidx, j) = ireg
         
         tau = - CellRay1D%LenSeg(k) * xst(ireg)
         
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
      
      tau = - FastRaySeg%LenSeg(l) * xst(ireg)
      
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
    END DO
  END DO
END DO
! ----------------------------------------------------
ALLOCATE (phiobdPolar(1:nPolarAng))

PhiAnginSvIdx  = RayInfo%PhiAngInSvIdx (iRotRay, krot)
PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay ,krot)

phiobdPolar = PhiAngIn(:, PhiAnginSvIdx)

jinc = 1; jbeg = 1; jend = nCoreRay
IF (krot .EQ. 2) THEN ! Backward Sweep
  jinc = -1; jbeg = nCoreRay; jend = 1
END IF

DO j = jbeg, jend, jinc
  iazi = CoreRay(RotRay(irotray)%RayIdx(j))%iang
  
  idir = RotRay(i)%DIR(j)
  IF (krot .EQ. 2) idir = mp(idir)  !Reverse the sweep direction
  
  IF (lJout) wt2(1:4) = wtang2(ipol, iazi, 1:4)
  
  nRaySeg = nTotRaySeg(j)
  
  IF (idir .EQ. 1) THEN ! Forward Sweep
    PhiAngOutPolar(:, 1) = phiobdPolar(:)
    
    DO ir = 1, nRaySeg
      ireg = FsrIdx(ir, j)
      
      DO ipol = 1, nPolarAng
        wt   = wtang(ipol, iazi)
        phid = (PhiAngOutPolar(ipol, ir) - src(ireg)) * ExpAppPolar(ipol, ir, j)
        
        PhiAngOutPolar(ipol, ir+1) = PhiAngOutPolar(ipol, ir) - phid
        
        phis(ireg) = phis(ireg) + wt*phid
      END DO
    END DO
    
    phiobdPolar(:) = PhiAngOutPolar(:, nRaySeg+1)
    
    ! Surface
    IF (ljout) THEN
      DO ir = 1, nTotCellRay(j)
        DO ipol = 1, nPolarANg
          wt    = wtang(ipol, iazi)
          icel  = PinIdx (ir, j)
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
      ir   = ir - 1
      ireg = FsrIdx(ir, j)
      
      DO ipol = 1, nPolarAng
        wt   = wtang(ipol, iazi)
        phid = (PhiAngOutPolar(ipol, ir + 2) - src(ireg)) * ExpAppPolar(ipol, ir, j)
        
        PhiAngOutPolar(ipol, ir+1) = PhiAngOutPolar(ipol, ir + 2) - phid
        
        phis(ireg) = phis(ireg) + wt * phid
      END DO
    END DO
    
    phiobdPolar(:) = PhiAngOutPolar(:, 2)
    
    IF (lJout) THEN
      DO ir = 1, nTotCellRay(j)
        DO ipol = 1, nPolarAng
          wt    = wtang(ipol, iazi)
          icel  = PinIdx (ir, j)
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

PhiAngIn(:, PhiAngOutSvIdx) = phiobdPolar(:)
! ----------------------------------------------------
DEALLOCATE (phiobdPolar)

NULLIFY (AziAng)
NULLIFY (PolarAng)
NULLIFY (AsyRay)
NULLIFY (CoreRay)
NULLIFY (RotRay)
NULLIFY (CellRay)
NULLIFY (Asy)
NULLIFY (Pin)
NULLIFY (PinInfo)
NULLIFY (Cell)
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
! ----------------------------------------------------

END SUBROUTINE RecTrackRotRayGM_OMP
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayGM_AFSS(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, FastMocLv)

USE PARAM,   ONLY : TRUE, FALSE
USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type, Pin_Type, Asy_Type, AsyInfo_Type, PinInfo_Type, Cell_Type, AziAngleInfo_Type, PolarAngle_Type, ModRayInfo_type,  AsyRayInfo_type,  &
                    CoreRayInfo_Type, RotRayInfo_Type, CellRayInfo_type, FastCoreRayDat_Type, TrackingDat_Type, FastRaySegDat_Type
USE Moc_Mod, ONLY : nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay

IMPLICIT NONE

TYPE (RayInfo_Type),     INTENT(INOUT) :: RayInfo
TYPE (CoreInfo_Type),    INTENT(INOUT) :: CoreInfo
TYPE (TrackingDat_Type), INTENT(INOUT) :: TrackingDat

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz
INTEGER, INTENT(IN) :: FastMocLv
! ----------------------------------------------------
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

INTEGER :: mp(2)
DATA mp /2, 1/

INTEGER :: iazi, ipol, iCoreRay, iasyray, iceray, irayseg, irot, itype, idir, i, j, k, l, m, jbeg, jend, jinc, ir, ir1, iOmpAzi, ibcel
INTEGER :: ipin, icel, iasy, ireg, isurf, irsegidx, icellrayidx, PhiAnginSvIdx, PhiAngOutSvIdx
INTEGER :: nRotRay, nCoreRay, nAsyRay, nPinRay, nRaySeg, FsrIdxSt, nPolarAng, nAziAng, nPhiAngSv, nFsr, nxy

INTEGER, DIMENSION(nMaxCellRay, nMaxCoreRay, 2) :: CellRayIdxSt
INTEGER, DIMENSION(nMaxCellRay, nMaxCoreRay)    :: PinIdx
INTEGER, DIMENSION(nMaxCellRay, nMaxCoreRay, 2) :: SurfIdx

INTEGER, DIMENSION(nMaxCoreRay) :: nTotRaySeg
INTEGER, DIMENSION(nMaxCoreRay) :: nTotCellRay

REAL :: tau, phiobd, phid, wt
REAL :: wt2(4), wtang2(100,100,4)

REAL, POINTER, DIMENSION(:)     :: LenSeg, phis, src, xst, PhiAngOut
REAL, POINTER, DIMENSION(:,:)   :: OptLenList, ExpApp, PhiAngIn, PhiAngOutPolar, EXPA, EXPB, wtang
REAL, POINTER, DIMENSION(:,:,:) :: ExpAppPolar, jout, Phi1a, Phi2a
REAL, ALLOCATABLE, DIMENSION(:) :: phiobdPolar

LOGICAL :: lFast
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
ExpApp         => TrackingDat%ExpApp
ExpAppPolar    => TrackingDat%ExpAppPolar
Phis           => TrackingDat%phis
src            => TrackingDat%src
phi1a          => TrackingDat%phi1a
phi2a          => TrackingDat%phi2a
xst            => TrackingDat%xst
jout           => TrackingDat%jout
PhiAngOut      => TrackingDat%PhiAngOut
PhiAngOutPolar => TrackingDat%PhiAngOutPolar
PhiAngIn       => TrackingDat%phiAngIn
EXPA           => TrackingDat%EXPA
EXPB           => TrackingDat%EXPB
Wtang          => TrackingDat%wtang
! ----------------------------------------------------
IF (lJout) THEN
  DO ipol = 1, nPolarAng
    DO iazi = 1, nAziAng
      wtang2(ipol, iazi, 1) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv)
      wtang2(ipol, iazi, 3) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv)
      wtang2(ipol, iazi, 2) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv)
      wtang2(ipol, iazi, 4) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv)
    END DO
  END DO
END IF
! ----------------------------------------------------
i        = iRotRay
nCoreRay = RotRay(irotRay)%nRay

IF (.NOT. lFast) THEN
  DO j = 1, nCoreRay
    irsegidx    = 0
    icellrayidx = 0
    iCoreRay    = RotRay(iRotRay)%RayIdx(j)
    nAsyRay     = CoreRay(iCoreRay)%nRay
    
    DO k = 1, nAsyRay
      iasyray = CoreRay(iCoreRay)%AsyRayIdx(k)
      iasy    = CoreRay(iCoreRay)%AsyIdx(k)
      
      IF (iasy .EQ. 0) CYCLE
      
      nPinRay = AsyRay(iAsyRay)%nCellRay
      itype   = Asy(iasy)%PartialAsyFlag
      
      DO l = 1, nPinRay
        ipin   = AsyRay(iAsyRay)%PinIdx(l)
        iceray = AsyRay(iAsyRay)%PinRayIdx(l)
        ipin   = Asy(iAsy)%GlobalPinIdx(ipin)
        
        icel     = Pin(ipin)%Cell(iz)
        FsrIdxSt = Pin(ipin)%FsrIdxSt
        
        ibcel = Cell(icel)%basecellstr
        CellRay => Cell(ibcel)%CellRay(iceray)
        
        icellrayidx = icellrayidx + 1
        
        PinIdx      (icellrayidx, j)    = ipin
        CellRayIdxSt(icellrayidx, j, 2) = irsegidx + 1
        
        nRaySeg      = CellRay%nSeg
        LocalFsrIdx => CellRay%LocalFsrIdx
        LenSeg      => CellRay%LenSeg
        
        DO iRaySeg = 1, nRaySeg
          ireg = FsrIdxSt + CellRay%LocalFsrIdx(iRaySeg) - 1
          
          tau = -LenSeg(iRaySeg) * xst(ireg)
          tau = - CellRay%LenSeg(iRaySeg) * xst(ireg)
          
          irsegidx = irsegidx + 1
          
          FsrIdx    (irsegidx, j) = ireg
          OptLenList(irsegidx, j) = tau
          ExpAppIdx (irsegidx, j) = max(INT(tau), -40000)
          ExpAppIdx (irsegidx, j) = min(0, ExpAppIdx(irsegidx, j))
        END DO
        
        CellRayIdxSt(icellrayidx, j, 1) = irsegidx
        SurfIdx     (icellRayIdx, j, 1) = AsyRay(iAsyRay)%PinRaySurf(2, l) !OutSurface
        SurfIdx     (icellRayIdx, j, 2) = AsyRay(iAsyRay)%PinRaySurf(1, l) !Insurface
      END DO
    END DO
    
    nTotRaySeg (j) = irsegidx
    nTotCellRay(j) = icellRayIdx
  END DO
ELSE IF(FastMocLV .EQ. 1) THEN
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
       
       ipin     = PinIdx(l, j)
       icel     = Pin(ipin)%Cell(iz)
       FsrIdxSt = Pin(ipin)%FsrIdxSt
       
       DO k = FastRay%Ray1DIdx(1, l, j), FastRay%Ray1DIdx(2, l, j)
         irsegidx = irsegidx + 1
         ireg     = FsrIdxSt + CellRay1D%LocalFsrIdx(K) - 1
         
         FsrIdx(irsegidx, j) = ireg
         
         tau = - CellRay1D%LenSeg(k) * xst(ireg)
         
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
      
      tau = - FastRaySeg%LenSeg(l) * xst(ireg)
      
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
    END DO
  END DO
END DO
! ----------------------------------------------------
ALLOCATE (phiobdPolar(1:nPolarAng))

DO irot = 1, 2
  PhiAnginSvIdx  = RayInfo%PhiAngInSvIdx (iRotRay, irot)
  PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay, irot)
  
  phiobdPolar(:) = PhiAngIn(:, PhiAnginSvIdx)

  jinc = 1; jbeg = 1; jend = nCoreRay
  IF (irot .EQ. 2) THEN !Backward Sweep
    jinc = -1; jbeg = nCoreRay; jend = 1
  END IF
  
  DO j = jbeg, jend, jinc
    iazi    = CoreRay(RotRay(irotray)%RayIdx(j))%iang
    iompazi = RotRay(irotray)%OmpRayIdx(j)
    
    idir = RotRay(i)%DIR(j)
    IF(irot .EQ. 2) idir = mp(idir)
    
    IF (lJout) wt2(1:4) = wtang2(ipol, iazi, 1:4)
    
    nRaySeg = nTotRaySeg(j)

    IF (idir .EQ. 1) THEN ! Forward Sweep
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
            wt    = wtang(ipol, iazi)
            icel  = PinIdx (ir, j)
            isurf = SurfIdx(ir, j, 1)
            
            Jout(2, isurf, icel) = Jout(2, isurf, icel) + PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 1)+1) * wt
            Jout(3, isurf, icel) = Jout(3, isurf, icel) + PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 1)+1) * wt2(isurf)
            
            isurf = SurfIdx(ir, j, 2)
            
            Jout(1, isurf, icel) = Jout(1, isurf, icel) + PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 2)) * wt
            Jout(3, isurf, icel) = Jout(3, isurf, icel) + PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 2)) * wt2(isurf)
          END DO
        END DO
      END IF
    ELSE ! Backward Sweep
      PhiAngOutPolar(:, nRaySeg+2) = phiobdPolar(:)
      
      ir = nRaySeg + 1
      
      DO ir1 = 1, nRaySeg
        ir   = ir - 1
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
            wt    = wtang(ipol, iazi)
            icel  = PinIdx (ir, j)
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
  
  PhiAngIn(:, PhiAngOutSvIdx) = phiobdPolar(:)
END DO
! ----------------------------------------------------
DEALLOCATE (phiobdPolar)

NULLIFY (AziAng)
NULLIFY (PolarAng)
NULLIFY (AsyRay)
NULLIFY (CoreRay)
NULLIFY (RotRay)
NULLIFY (CellRay)
NULLIFY (Asy)
NULLIFY (Pin)
NULLIFY (PinInfo)
NULLIFY (Cell)
NULLIFY (FsrIdx)
NULLIFY (ExpAppIdx)
NULLIFY (OptLenList)
NULLIFY (ExpApp)
NULLIFY (LenSeg)
NULLIFY (LocalFsrIdx)
NULLIFY (Phis)
NULLIFY (src)
NULLIFY (Phi1a)
NULLIFY (Phi2a)
NULLIFY (xst)
NULLIFY (jout)
NULLIFY (PhiAngOut)
NULLIFY (PhiAngIn)
NULLIFY (EXPA)
NULLIFY (EXPB)
NULLIFY (wtang)
! ----------------------------------------------------

END SUBROUTINE RecTrackRotRayGM_AFSS
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayGM_OMP(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, krot, FastMocLv)

USE TYPEDEF,  ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, TrackingDat_Type
USE Moc_Mod,  ONLY : nMaxCellRay, nMaxCoreRay
USE HexData,  ONLY : hAsy
USE HexType,  ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexCoreRay, Type_HexRotRay
USE HexData,  ONLY : haRay, hcRay, hRotRay, hAsyTypInfo

IMPLICIT NONE

TYPE (RayInfo_Type),     INTENT(INOUT) :: RayInfo
TYPE (CoreInfo_Type),    INTENT(INOUT) :: CoreInfo
TYPE (TrackingDat_Type), INTENT(INOUT) :: TrackingDat

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz, krot
INTEGER, INTENT(IN) :: FastMocLv
! ----------------------------------------------------
INTEGER :: iazi, ipol, iAsyRay, iAsy, iSurf, PhiAnginSvIdx, PhiAngOutSvIdx, ifsr, irsegidx, icellrayidx
INTEGER :: icRay, imRay, jbeg, jend, jinc, ihpRay, iRaySeg, iGeoTyp, iAsyTyp, jhPin, icBss, jcBss, jcRay, iReg, iCel, iRaySeg1
INTEGER :: nCoreRay, nAsyRay, nPolarAng, nRaySeg

INTEGER :: CellRayIdxSt(nMaxCellRay, nMaxCoreRay, 2)
INTEGER :: PinIdx      (nMaxCellRay, nMaxCoreRay)
INTEGER :: SurfIdx     (nMaxCellRay, nMaxCoreRay, 2)
INTEGER :: nTotRaySeg  (nMaxCoreRay)
INTEGER :: nTotCellRay (nMaxCoreRay)

INTEGER, POINTER, DIMENSION(:,:) :: FsrIdx, ExpAppIdx

REAL :: tau, phid, wt

REAL, DIMENSION(RayInfo%nPolarAngle) :: locphiout

REAL, POINTER, DIMENSION(:)     :: phis, src, xst
REAL, POINTER, DIMENSION(:,:)   :: EXPA, EXPB, wtang, OptLenList, PhiAngOutPolar
REAL, POINTER, DIMENSION(:,:,:) :: jout, ExpAppPolar

TYPE (Pin_Type), POINTER, DIMENSION(:) :: Pin

TYPE (Type_HexAsyRay),  POINTER :: haRay_Loc
TYPE (Type_HexCelRay),  POINTER :: CelRay_Loc
TYPE (Type_HexRotRay),  POINTER :: hRotRay_Loc
! ----------------------------------------------------

nPolarAng      = RayInfo%nPolarAngle
PhiAngInSvIdx  = RayInfo%PhiAngInSvIdx (iRotRay, krot)
PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay, krot)

Pin => CoreInfo%Pin

phis           => TrackingDat%phis
src            => TrackingDat%src
xst            => TrackingDat%xst
jout           => TrackingDat%jout
EXPA           => TrackingDat%EXPA
EXPB           => TrackingDat%EXPB
wtang          => TrackingDat%wtang
ExpAppPolar    => TrackingDat%ExpAppPolar
PhiAngOutPolar => TrackingDat%PhiAngOutPolar
FsrIdx         => TrackingDat%FsrIdx
ExpAppIdx      => TrackingDat%ExpAppIdx
OptLenList     => TrackingDat%OptLenList
locphiout      =  TrackingDat%PhiAngIn(:, PhiAnginSvIdx)

hRotRay_Loc => hRotRay(iRotRay)
nCoreRay     = hRotRay_Loc%ncRay
! ----------------------------------------------------
DO icRay = 1, nCoreRay
  jcRay   = abs(hRotRay_Loc%cRayIdx(icRay))
  nAsyRay = hcRay(jcRay)%nmRay
  
  irSegIdx    = 0
  iCellRayIdx = 0
  
  DO imRay = 1, nAsyRay
    iAsyRay = hcRay(jcRay)%mRayIdx(imRay)
    iAsy    = hcRay(jcRay)%AsyIdx(imRay)
    iAsyTyp = hAsy(iAsy)%AsyTyp
    iGeoTyp = hAsy(iAsy)%GeoTyp
    icBss   = hAsyTypInfo(iAsyTyp)%iBss
    
    haRay_Loc => haRay(iGeoTyp, icBss, iAsyRay)
    
    DO ihpRay = 1, haRay_Loc%nhpRay
      iCellRayIdx = iCellRayIdx + 1
      
      jhPin = haRay_Loc%CelRay(ihpRay)%hPinIdx
      jhPin = hAsy(iAsy)%PinIdxSt + hAsyTypInfo(iAsyTyp)%PinLocIdx(iGeoTyp, jhPin) - 1
      jcBss = Pin(jhPin)%hCelGeo(iz)
      
      CelRay_Loc => haRay(iGeoTyp, jcBss, iAsyRay)%CelRay(ihpRay)
      
      ! Start of Cell
      CellRayIdxSt(iCellRayIdx, icRay, 2) = irSegIdx + 1
      SurfIdx     (iCellRayIdx, icRay, 2) = CelRay_Loc%SurfIdx(1) ! y : Small
      
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
      SurfIdx     (iCellRayIdx, icRay, 1) = CelRay_Loc%SurfIdx(2) ! y : Big
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
      ExpAppPolar(iPol, iRaySeg, icRay) = EXPA(ExpAppIdx(iRaySeg, icRay), iPol) * optlenlist(iRaySeg, icRay) &
                                        + EXPB(ExpAppIdx(iRaySeg, icRay), iPol)
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
  iAzi  = hcRay(abs(jcRay))%AzmIdx
  
  IF (krot .EQ. 2) jcRay = -jcRay !Reverse the Sweep Direction
  
  nRaySeg = nTotRaySeg(icRay)
  ! --------------------------------------------------
  IF (jcRay > 0) THEN
    PhiAngOutPolar(:, 1) = locphiout(:)
    
    DO iRaySeg = 1, nRaySeg
      iReg = FsrIdx(iRaySeg, icRay)
      
      DO iPol = 1, nPolarAng
        wt   = wtang(iPol, iAzi)
        phid = (PhiAngOutPolar(iPol, iRaySeg) - src(iReg)) * ExpAppPolar(iPol, iRaySeg, icRay)
        
        PhiAngOutPolar(ipol, iRaySeg+1) = PhiAngOutPolar(iPol, iRaySeg) - phid
        
        phis(iReg) = phis(iReg) + wt * phid
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
        phid = (PhiAngOutPolar(iPol, iRaySeg + 2) - src(iReg)) * ExpAppPolar(iPol, iRaySeg, icRay)
        
        PhiAngOutPolar(iPol, iRaySeg+1) = PhiAngOutPolar(iPol, iRaySeg + 2) - phid
        
        phis(iReg) = phis(iReg) + wt * phid
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
NULLIFY (phis)
NULLIFY (src)
NULLIFY (xst)
NULLIFY (EXPA)
NULLIFY (EXPB)
NULLIFY (Pin)
NULLIFY (FsrIdx)
NULLIFY (ExpAppIdx)
NULLIFY (wtang)
NULLIFY (OptLenList)
NULLIFY (jout)
NULLIFY (ExpAppPolar)
NULLIFY (PhiAngOutPolar)
NULLIFY (haRay_Loc)
NULLIFY (CelRay_Loc)
NULLIFY (hRotRay_Loc)
! ----------------------------------------------------

END SUBROUTINE HexTrackRotRayGM_OMP
! ------------------------------------------------------------------------------------------------------------