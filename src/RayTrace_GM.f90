#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTrace_GM(RayInfo, CoreInfo, phis1g, PhiAngIn1g, xst1g, src1g, jout1g, iz, ljout, FastMocLv)

USE OMP_LIB
USE PARAM,   ONLY : ZERO
USE TYPEDEF, ONLY : RayInfo_Type, CoreInfo_type, Pin_Type, Cell_Type
USE Moc_Mod, ONLY : RecTrackRotRay_GM, HexTrackRotRay_GM, TrackingDat, wtang
USE PE_MOD,  ONLY : PE
USE geom,    ONLY : nbd
USE CNTL,    ONLY : nTracerCntl

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:)     :: phis1g, xst1g, src1g
REAL, POINTER, DIMENSION(:,:)   :: PhiAngIn1g
REAL, POINTER, DIMENSION(:,:,:) :: jout1g

INTEGER :: iz
LOGICAL :: ljout
INTEGER, OPTIONAL :: FastMocLv
! ----------------------------------------------------
TYPE (Pin_Type),  POINTER, DIMENSION(:) :: Pin
TYPE (Cell_Type), POINTER, DIMENSION(:) :: Cell

INTEGER :: ithr, iRotRay, krot, ixy, icel, ifsr, jfsr, iazi, ipol, ibd, nThr, nxy, FsrIdxSt, nFsr, nAzi, nPol
REAL :: phia1gp
LOGICAL :: lAFSS
! ----------------------------------------------------

nFsr = CoreInfo%nCoreFsr
nxy  = CoreInfo%nxy

lAFSS = nTracerCntl%lAFSS

nAzi = RayInfo%nAziAngle
nPol = RayInfo%nPolarAngle

nthr = PE%nThread
CALL omp_set_num_threads(nthr)
! ----------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ithr, krot, iRotRay)
ithr = omp_get_thread_num() + 1

TrackingDat(ithr)%phis1g = ZERO
IF (ljout) TrackingDat(ithr)%jout1g  = ZERO
IF (lAFSS) TrackingDat(ithr)%phia1g1 = ZERO

TrackingDat(ithr)%src1g      => src1g
TrackingDat(ithr)%xst1g      => xst1g
TrackingDat(ithr)%PhiAngIn1g => PhiAngIn1g

IF (nTracerCntl%lHex) THEN
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
  DO krot = 1, 2
    DO iRotRay = 1, RayInfo%nRotRay
      CALL HexTrackRotRay_GM(RayInfo, CoreInfo, TrackingDat(ithr), ljout, iRotRay, iz, krot, lAFSS)
    END DO
  END DO
  !$OMP END DO
ELSE
  !$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
  DO krot = 1, 2
    DO iRotRay = 1, RayInfo%nRotRay
      CALL RecTrackRotRay_GM(RayInfo, CoreInfo, TrackingDat(ithr), ljout, iRotRay, iz, krot, FastMocLv)
    END DO
  END DO
  !$OMP END DO
END IF
!$OMP END PARALLEL
! ----------------------------------------------------
phis1g = ZERO

IF (lAFSS) THEN
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ifsr, iazi, ipol, phia1gp, ithr)
  !$OMP DO SCHEDULE(GUIDED)
  DO ifsr = 1, nFsr
    DO iazi = 1, nAzi
      DO ipol = 1, nPol
        phia1gp = ZERO
        
        DO ithr = 1, nthr
          phia1gp = phia1gp + TrackingDat(ithr)%phia1g1(ipol, ifsr, iazi)
        END DO
        
        phis1g(ifsr) = phis1g(ifsr) + wtang(ipol, iazi) * phia1gp
      END DO
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
ELSE
  DO ithr = 1, nthr
    DO ifsr = 1, nfsr
      phis1g(ifsr) = phis1g(ifsr) + TrackingDat(ithr)%phis1g(ifsr)
    END DO
  END DO
END IF

IF (ljout) THEN
  jout1g = ZERO
  
  DO ithr = 1, nthr
    DO ixy = 1, nxy
      DO ibd = 1, nbd
        Jout1g(:, ibd, ixy) = Jout1g(:, ibd, ixy) + TrackingDat(ithr)%Jout1g(:, ibd, ixy)
      END DO
    END DO
  END DO
END IF
! ----------------------------------------------------
Pin  => CoreInfo%Pin
Cell => CoreInfo%CellInfo

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ixy, FsrIdxSt, icel, ifsr, jfsr)
!$OMP DO SCHEDULE(GUIDED)
DO ixy = 1, nxy
  FsrIdxSt = Pin(ixy)%FsrIdxSt
  icel     = Pin(ixy)%Cell(iz)
  
  DO ifsr = 1, Cell(icel)%nFsr
    jfsr = FsrIdxSt + ifsr - 1
    
    phis1g(jfsr) = phis1g(jfsr) / xst1g(jfsr) / Cell(icel)%vol(ifsr) + src1g(jfsr)
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

NULLIFY (Pin)
NULLIFY (Cell)
! ----------------------------------------------------

END SUBROUTINE RayTrace_GM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRay_GM(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, krot, FastMocLv)

USE PARAM,   ONLY : TRUE, FALSE, ZERO
USE TYPEDEF, ONLY : RayInfo_Type, coreinfo_type, Pin_Type, Asy_Type, AsyInfo_Type, PinInfo_Type, Cell_Type, AziAngleInfo_Type, PolarAngle_Type, ModRayInfo_type, AsyRayInfo_type,  &
                    CoreRayInfo_Type, RotRayInfo_Type, CellRayInfo_type, FastCoreRayDat_Type, TrackingDat_Type, FastRaySegDat_Type
USE Moc_Mod, ONLY : nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay

IMPLICIT NONE

TYPE(RayInfo_Type),     INTENT(INOUT) :: RayInfo
TYPE(CoreInfo_Type),    INTENT(INOUT) :: CoreInfo
TYPE(TrackingDat_Type), INTENT(INOUT) :: TrackingDat

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz, krot, FastMocLv

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

REAL, POINTER, DIMENSION(:)     :: LenSeg, phis1g, src1g, xst1g
REAL, POINTER, DIMENSION(:,:)   :: OptLenList, PhiAngOut1g, PhiAngIn1g, ExpA, ExpB, wtang
REAL, POINTER, DIMENSION(:,:,:) :: ExpAppPolar, jout1g, wtsurf

INTEGER :: mp(2)
INTEGER :: iazi, ipol, iCoreRay, iasyray, iceray, irayseg, itype, idir, nRotRay, nCoreRay, nAsyRay, nPinRay, nRaySeg, FsrIdxSt, nPolarAng, nAziAng, nPhiAngSv
INTEGER :: ipin, icel, iasy, ifsr, isurf, irsegidx, icellrayidx, PhiAnginSvIdx, PhiAngOutSvIdx, nFsr, nxy, i, j, k, l, m, jbeg, jend, jinc, ir, ir1, ibcel

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
src1g       => TrackingDat%src1g
xst1g       => TrackingDat%xst1g
jout1g      => TrackingDat%jout1g
PhiAngOut1g => TrackingDat%PhiAngOut1g
PhiAngIn1g  => TrackingDat%phiAngIn1g
ExpA        => TrackingDat%ExpA
ExpB        => TrackingDat%ExpB
wtang       => TrackingDat%wtang
wtsurf      => TrackingDat%wtsurf
! ----------------------------------------------------
i = iRotRay
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
      
      IF (iasy .EQ. 0)  CYCLE ! Skip Dummy Assembly
      
      nPinRay = AsyRay(iAsyRay)%nCellRay
      itype   = Asy(iasy)%PartialAsyFlag
      
      DO l = 1, nPinRay   !Pin Ray Sweep
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
          ifsr = FsrIdxSt + CellRay%LocalFsrIdx(iRaySeg) - 1
          
          tau = -LenSeg(iRaySeg) * xst1g(ifsr)
          tau = -CellRay%LenSeg(iRaySeg) * xst1g(ifsr)
          
          irsegidx = irsegidx + 1
          
          FsrIdx    (irsegidx, j) = ifsr
          OptLenList(irsegidx, j) = tau
          ExpAppIdx (irsegidx, j) = max(INT(tau), -40000)
          ExpAppIdx (irsegidx, j) = min(0, ExpAppIdx(irsegidx, j))
        END DO
        CellRayIdxSt(icellrayidx, j, 1) = irsegidx
        SurfIdx     (icellRayIdx, j, 1) = AsyRay(iAsyRay)%PinRaySurf(2, l) ! OutSurface
        SurfIdx     (icellRayIdx, j, 2) = AsyRay(iAsyRay)%PinRaySurf(1, l) ! Insurface
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
         ifsr     = FsrIdxSt + CellRay1D%LocalFsrIdx(K) - 1
         
         FsrIdx(irsegidx, j) = ifsr
         
         tau = -CellRay1D%LenSeg(k) * xst1g(ifsr)
         
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
      ifsr = FastRaySeg%FsrIdx(l)
      
      FsrIdx(l, j) = ifsr
      
      tau = -FastRaySeg%LenSeg(l) * xst1g(ifsr)
      
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
      ExpAppPolar(ipol, l, j) = ExpA(ExpAppIdx(l, j), ipol) * optlenlist(l, j) + ExpB(ExpAppIdx(l, j), ipol)
    END DO
  END DO
END DO

ALLOCATE (phiobdPolar(1:nPolarAng))

PhiAnginSvIdx  = RayInfo%PhiAngInSvIdx (iRotRay, krot)
PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay ,krot)

phiobdPolar = PhiAngIn1g(:, PhiAnginSvIdx)

jinc = 1; jbeg = 1; jend = nCoreRay
IF (krot .EQ. 2) THEN ! Backward Sweep
  jinc = -1; jbeg = nCoreRay; jend = 1
END IF
! ----------------------------------------------------
DO j = jbeg, jend, jinc
  idir = RotRay(i)%DIR(j)
  iazi = CoreRay(RotRay(irotray)%RayIdx(j))%iang
  
  IF (krot .eq. 2) idir = mp(idir) ! Reverse the sweep direction
  
  IF (lJout) wt2(1:4) = wtsurf(ipol, iazi, 1:4)
    
  nRaySeg = nTotRaySeg(j)
  
  IF (idir .EQ. 1) THEN ! Forward
    PhiAngOut1g(:, 1) = phiobdPolar(:)
    
    ! Iter. : FSR
    DO ir = 1, nRaySeg
      ifsr = FsrIdx(ir, j)
      
      DO ipol = 1, nPolarAng
        wt = wtang(ipol, iazi)
        
        phid = (PhiAngOut1g(ipol, ir) - src1g(ifsr)) * ExpAppPolar(ipol, ir, j)
        
        PhiAngOut1g(ipol, ir+1) = PhiAngOut1g(ipol, ir) - phid
        
        phis1g(ifsr) = phis1g(ifsr) + wt*phid
      END DO
    END DO
    
    phiobdPolar(:) = PhiAngOut1g(:, nRaySeg+1)
    
    ! Surf.
    IF (ljout) THEN
      DO ir = 1, nTotCellRay(j)
        DO ipol = 1, nPolarANg
          wt = wtang(ipol, iazi)
          
          icel  = PinIdx (ir, j)
          isurf = SurfIdx(ir, j, 1)
          
          Jout1g(2, isurf, icel) = Jout1g(2, isurf, icel) + PhiAngOut1g(ipol, CellRayIdxSt(ir, j, 1)+1) * wt
          Jout1g(3, isurf, icel) = Jout1g(3, isurf, icel) + PhiAngOut1g(ipol, CellRayIdxSt(ir, j, 1)+1) * wt2(isurf)
          
          isurf = SurfIdx(ir, j, 2)
          
          Jout1g(1, isurf, icel) = Jout1g(1, isurf, icel) + PhiAngOut1g(ipol, CellRayIdxSt(ir, j, 2)) * wt
          Jout1g(3, isurf, icel) = Jout1g(3, isurf, icel) + PhiAngOut1g(ipol, CellRayIdxSt(ir, j, 2)) * wt2(isurf)
        END DO
      END DO
    END IF
  ELSE ! Backward
    PhiAngOut1g(:, nRaySeg+2) = phiobdPolar(:)
    
    ir = nRaySeg + 1
    
    ! Iter. : FSR
    DO ir1 = 1, nRaySeg
      ir   = ir - 1
      ifsr = FsrIdx(ir, j)
      
      DO ipol = 1, nPolarAng
        wt = wtang(ipol, iazi)
        
        phid = (PhiAngOut1g(ipol, ir + 2) - src1g(ifsr)) * ExpAppPolar(ipol, ir, j)
        
        PhiAngOut1g(ipol, ir+1) = PhiAngOut1g(ipol, ir + 2) - phid
        
        phis1g(ifsr) = phis1g(ifsr) + wt * phid
      END DO
    END DO
    
    phiobdPolar(:) = PhiAngOut1g(:, 2)
    
    ! Surf.
    IF (lJout) THEN
      DO ir = 1, nTotCellRay(j)
        DO ipol = 1, nPolarAng
          wt = wtang(ipol, iazi)
          
          icel  = PinIdx (ir, j)
          isurf = SurfIdx(ir, j, 2)
          
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

PhiAngIn1g(:,PhiAngOutSvIdx) = phiobdPolar(:)

DEALLOCATE (phiobdPolar)
! ----------------------------------------------------
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
NULLIFY (src1g)
NULLIFY (xst1g)
NULLIFY (jout1g)
NULLIFY (PhiAngOut1g)
NULLIFY (PhiAngIn1g)
NULLIFY (ExpA)
NULLIFY (ExpB)
NULLIFY (wtang)
NULLIFY (wtsurf)
! ----------------------------------------------------

END SUBROUTINE RecTrackRotRay_GM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRay_GM(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, krot, lAFSS)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, TrackingDat_Type, Pin_Type
USE HexType, ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexCoreRay, Type_HexRotRay
USE HexData, ONLY : hAsy, haRay, hcRay, hRotRay, hAsyTypInfo

IMPLICIT NONE

TYPE (RayInfo_Type)     :: RayInfo
TYPE (CoreInfo_Type)    :: CoreInfo
TYPE (TrackingDat_Type) :: TrackingDat

LOGICAL, INTENT(IN) :: ljout, lAFSS
INTEGER, INTENT(IN) :: irotray, iz, krot
! ----------------------------------------------------
INTEGER :: iazi, ipol, icRay, jcRay, iaRay, jaRay, iRaySeg, ihpRay, iAsy, ifsr, iSurf, jbeg, jend, jinc, ig, iGeoTyp, iAsyTyp, jhPin, icBss, jcBss
INTEGER :: nCoreRay, nAsyRay, nPolarAng, PhiAnginSvIdx, PhiAngOutSvIdx, ExpAppIdx
INTEGER :: iast, iaed, iainc, ipst, iped, ipinc, isfst, isfed, isgst, isged, isginc

REAL :: phid, tau, locsrc, ExpApp

REAL :: wtazi(10)
REAL, DIMENSION(RayInfo%nPolarAngle) :: PhiAngOut

REAL, POINTER, DIMENSION(:)     :: phis1g, src1g, xst1g
REAL, POINTER, DIMENSION(:,:)   :: PhiAngIn1g, ExpA, ExpB, wtang
REAL, POINTER, DIMENSION(:,:,:) :: Jout1g, phia1g

TYPE (Pin_Type), POINTER, DIMENSION(:) :: Pin

TYPE (Type_HexAsyRay), POINTER :: haRay_Loc
TYPE (Type_HexCelRay), POINTER :: CelRay_Loc
TYPE (Type_HexRotRay), POINTER :: hRotRay_Loc
! ----------------------------------------------------

! Ray
nPolarAng      = RayInfo%nPolarAngle
PhiAngInSvIdx  = RayInfo%PhiAngInSvIdx (iRotRay, krot)
PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay, krot)

hRotRay_Loc => hRotRay(iRotRay)
nCoreRay     = hRotRay_Loc%ncRay

! Geo.
Pin => CoreInfo%Pin

! Loc.
phis1g     => TrackingDat%phis1g
src1g      => TrackingDat%src1g
xst1g      => TrackingDat%xst1g
Jout1g     => TrackingDat%Jout1g
PhiAngIn1g => TrackingDat%PhiAngIn1g
wtang      => TrackingDat%wtang
ExpA       => TrackingDat%ExpA
ExpB       => TrackingDat%ExpB

IF (lAFSS) phia1g => TrackingDat%phia1g1

! Iter.
PhiAngOut(1:nPolarAng) = PhiAngIn1g(1:nPolarAng, PhiAnginSvIdx)

IF (krot .EQ. 1) THEN
  jbeg = 1; jend = nCoreRay; jinc = 1
ELSE
  jend = 1; jbeg = nCoreRay; jinc = -1
END IF
! ----------------------------------------------------
DO icRay = jbeg, jend, jinc
  jcRay   = hRotRay_Loc%cRayIdx(icRay)
  nAsyRay = hcRay(abs(jcRay))%nmRay
  iazi    = hcRay(abs(jcRay))%AzmIdx
  
  IF (krot .EQ. 2) jcRay = -jcRay ! Reverse the Sweep Direction
  
  DO ipol = 1, nPolarAng
    wtazi(ipol) = wtang(ipol, iazi)
  END DO
  ! --------------------------------------------------
  IF (jcRay .GT. 0) THEN
    iast = 1; iaed = nAsyRay; iainc = 1
  ELSE
    iaed = 1; iast = nAsyRay; iainc = -1
  END IF
  
  DO iaRay = iast, iaed, iainc
    jaRay = hcRay(abs(jcRay))%mRayIdx(iaRay)
    iAsy  = hcRay(abs(jcRay))%AsyIdx (iaRay)
    
    IF (iAsy .EQ. 0) CYCLE
    
    iAsyTyp = hAsy(iAsy)%AsyTyp
    iGeoTyp = hAsy(iAsy)%GeoTyp
    icBss   = hAsyTypInfo(iAsyTyp)%iBss
    
    haRay_Loc => haRay(iGeoTyp, icBss, jaRay)
    
    IF (jcRay .GT. 0) THEN
      ipst = 1; iped = haRay_Loc%nhpRay; ipinc = 1;  isfst = 1; isfed = 2
    ELSE
      iped = 1; ipst = haRay_Loc%nhpRay; ipinc = -1; isfed = 1; isfst = 2
    END IF
    
    DO ihpRay = ipst, iped, ipinc
      jhPin = haRay_Loc%CelRay(ihpRay)%hPinIdx
      jhPin = hAsy(iAsy)%PinIdxSt + hAsyTypInfo(iAsyTyp)%PinLocIdx(iGeoTyp, jhPin) - 1
      jcBss = Pin(jhPin)%hCelGeo(iz)
      
      CelRay_Loc => haRay(iGeoTyp, jcBss, jaRay)%CelRay(ihpRay)
      
      ! Surface : In-coming
      IF (lJout) THEN
        iSurf = CelRay_Loc%hSufIdx(isfst)
        
        DO ipol = 1, nPolarAng
          Jout1g(1, iSurf, jhPin) = Jout1g(1, isurf, jhPin) + PhiAngOut(ipol) * wtazi(ipol)
        END DO
      END IF
      
      ! Iter. : FSR
      IF (jcRay .GT. 0) THEN
        isgst = 1; isged = CelRay_Loc%nSegRay; isginc = 1
      ELSE
        isged = 1; isgst = CelRay_Loc%nSegRay; isginc = -1
      END IF
      
      DO iRaySeg = isgst, isged, isginc
        ifsr = CelRay_Loc%MshIdx(iRaySeg) + Pin(jhPin)%FsrIdxSt - 1
        
        tau = -CelRay_Loc%SegLgh(iRaySeg) * xst1g(ifsr) ! Optimum Length
        
        ExpAppIdx = max(INT(tau), -40000)
        ExpAppIdx = min(0, ExpAppIdx)
        
        locsrc = src1g(ifsr)
        
        DO ipol = 1, nPolarAng
          ExpApp = ExpA(ExpAppIdx, ipol) * tau + ExpB(ExpAppIdx, ipol)
          
          phid = (PhiAngOut(ipol) - locsrc) * ExpApp
          
          PhiAngOut(ipol) = PhiAngOut(ipol) - phid
          
          IF (lAFSS) THEN
            phia1g(ipol, ifsr, iazi) = phia1g(ipol, ifsr, iazi) + phid
          ELSE
            phis1g(ifsr) = phis1g(ifsr) + wtazi(ipol) * phid
          END IF
        END DO
      END DO
      
      ! Surface : Out-going
      IF (lJout) THEN
        isurf = CelRay_Loc%hSufIdx(isfed)
        
        DO ipol = 1, nPolarAng
          Jout1g(2, iSurf, jhPin) = Jout1g(2, iSurf, jhPin) + PhiAngOut(ipol) * wtazi(ipol)
        END DO
      END IF
    END DO
  END DO
END DO

PhiAngIn1g(1:nPolarAng, PhiAngOutSvIdx) = PhiAngOut(1:nPolarAng)
! ----------------------------------------------------
! Loc.
NULLIFY (phis1g)
NULLIFY (src1g)
NULLIFY (xst1g)
NULLIFY (Jout1g)
NULLIFY (PhiAngIn1g)
NULLIFY (wtang)
NULLIFY (ExpA)
NULLIFY (ExpB)
NULLIFY (Pin)

IF (lAFSS) NULLIFY (phia1g)

! Hex.
NULLIFY (hRotRay_Loc)
NULLIFY (CelRay_Loc)
NULLIFY (hRotRay_Loc)
! ----------------------------------------------------

END SUBROUTINE HexTrackRotRay_GM
! ------------------------------------------------------------------------------------------------------------