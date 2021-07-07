#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceGM_OMP(RayInfo, CoreInfo, phis, PhiAngIn, xst, src, jout, iz, ljout, FastMocLv)

USE OMP_LIB
USE PARAM,   ONLY : ZERO
USE TYPEDEF, ONLY : RayInfo_Type, CoreInfo_type, Pin_Type, Cell_Type
USE Moc_Mod, ONLY : TrackingDat
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
TYPE (Cell_Type), POINTER, DIMENSION(:) :: Cell
TYPE (Pin_Type),  POINTER, DIMENSION(:) :: Pin

INTEGER :: nFsr, nxy, nThread
INTEGER :: ithr, FsrIdxSt, icel, iRotRay, ifsr, jfsr, ixy, krot
! ----------------------------------------------------

nFsr = CoreInfo%nCoreFsr
nxy  = CoreInfo%nxy

nThread = PE%nThread
! ----------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ithr, krot, iRotRay)
ithr = omp_get_thread_num() + 1

TrackingDat(ithr)%phis = ZERO
IF (ljout) TrackingDat(ithr)%jout = ZERO

TrackingDat(ithr)%src      => src
TrackingDat(ithr)%xst      => xst
TrackingDat(ithr)%PhiAngIn => PhiAngIn

DO krot = 1, 2
  IF (nTracerCntl%lHex) THEN
    !$OMP DO SCHEDULE(GUIDED)
    DO iRotRay = 1, RayInfo%nRotRay
      CALL HexTrackRotRayGM_OMP(RayInfo, CoreInfo, TrackingDat(ithr), ljout, iRotRay, iz, krot, FastMocLv)
    END DO
    !$OMP END DO NOWAIT
  ELSE
    !$OMP DO SCHEDULE(GUIDED)
    DO iRotRay = 1, RayInfo%nRotRay
      CALL RecTrackRotRayGM_OMP(RayInfo, CoreInfo, TrackingDat(ithr), ljout, iRotRay, iz, krot, FastMocLv)
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

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ixy, FsrIdxSt, icel, ifsr, jfsr)
!$OMP DO
DO ixy = 1, nxy
  FsrIdxSt = Pin(ixy)%FsrIdxSt
  icel     = Pin(ixy)%Cell(iz)
  
  DO ifsr = 1, Cell(icel)%nFsr
    jfsr = FsrIdxSt + ifsr - 1
    
    phis(jfsr) = phis(jfsr) / xst(jfsr) / Cell(icel)%vol(ifsr) + src(jfsr)
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

NULLIFY (Cell)
NULLIFY (Pin)
! ----------------------------------------------------

END SUBROUTINE RayTraceGM_OMP
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayGM_OMP(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, krot, FastMocLv)

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

REAL, POINTER, DIMENSION(:)     :: LenSeg, phis, src, xst
REAL, POINTER, DIMENSION(:,:)   :: OptLenList, PhiAngOutPolar, PhiAngIn, ExpA, ExpB, wtang
REAL, POINTER, DIMENSION(:,:,:) :: ExpAppPolar, jout, wtsurf

INTEGER :: mp(2)
INTEGER :: iazi, ipol, iCoreRay, iasyray, iceray, irayseg, itype, idir, nRotRay, nCoreRay, nAsyRay, nPinRay, nRaySeg, FsrIdxSt, nPolarAng, nAziAng, nPhiAngSv
INTEGER :: ipin, icel, iasy, ireg, isurf, irsegidx, icellrayidx, PhiAnginSvIdx, PhiAngOutSvIdx, nFsr, nxy, i, j, k, l, m, jbeg, jend, jinc, ir, ir1, ibcel

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
ExpA           => TrackingDat%ExpA
ExpB           => TrackingDat%ExpB
wtang          => TrackingDat%wtang
wtsurf         => TrackingDat%wtsurf
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
          ireg = FsrIdxSt + CellRay%LocalFsrIdx(iRaySeg) - 1
          
          tau = -LenSeg(iRaySeg) * xst(ireg)
          tau = -CellRay%LenSeg(iRaySeg) * xst(ireg)
          
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
      ExpAppPolar(ipol, l, j) = ExpA(ExpAppIdx(l, j), ipol) * optlenlist(l, j) + ExpB(ExpAppIdx(l, j), ipol)
    END DO
  END DO
END DO

ALLOCATE (phiobdPolar(1:nPolarAng))

PhiAnginSvIdx  = RayInfo%PhiAngInSvIdx (iRotRay, krot)
PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay ,krot)

phiobdPolar = PhiAngIn(:, PhiAnginSvIdx)

jinc = 1; jbeg = 1; jend = nCoreRay
IF (krot .EQ. 2) THEN ! Backward Sweep
  jinc = -1; jbeg = nCoreRay; jend = 1
END IF
! ----------------------------------------------------
DO j = jbeg, jend, jinc
  idir = RotRay(i)%DIR(j)
  iazi = CoreRay(RotRay(irotray)%RayIdx(j))%iang
  
  IF (krot .eq. 2) idir = mp(idir) ! Reverse the sweep direction
  
  IF (lJout) THEN
    wt2(1:4) = wtsurf(ipol, iazi, 1:4)
  END IF
  
  nRaySeg = nTotRaySeg(j)
  
  IF (idir .EQ. 1) THEN ! Forward Sweep
    PhiAngOutPolar(:, 1) = phiobdPolar(:)
    
    DO ir = 1, nRaySeg
      ireg = FsrIdx(ir, j)
      
      DO ipol = 1, nPolarAng
        wt = wtang(ipol, iazi)
        
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
          wt = wtang(ipol, iazi)
          
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
        wt = wtang(ipol, iazi)
        
        phid = (PhiAngOutPolar(ipol, ir + 2) - src(ireg)) * ExpAppPolar(ipol, ir, j)
        
        PhiAngOutPolar(ipol, ir+1) = PhiAngOutPolar(ipol, ir + 2) - phid
        
        phis(ireg) = phis(ireg) + wt * phid
      END DO
    END DO
    
    phiobdPolar(:) = PhiAngOutPolar(:, 2)
    
    IF (lJout) THEN
      DO ir = 1, nTotCellRay(j)
        DO ipol = 1, nPolarAng
          wt = wtang(ipol, iazi)
          
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

PhiAngIn(:,PhiAngOutSvIdx) = phiobdPolar(:)

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
! ----------------------------------------------------

END SUBROUTINE RecTrackRotRayGM_OMP
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayGM_OMP(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, krot, FastMocLv)

USE TYPEDEF,  ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, TrackingDat_Type
USE Moc_Mod,  ONLY : nMaxCellRay, nMaxCoreRay
USE HexType,  ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexCoreRay, Type_HexRotRay
USE HexData,  ONLY : hAsy, haRay, hcRay, hRotRay, hAsyTypInfo

IMPLICIT NONE

TYPE (RayInfo_Type),     INTENT(INOUT) :: RayInfo
TYPE (CoreInfo_Type),    INTENT(INOUT) :: CoreInfo
TYPE (TrackingDat_Type), INTENT(INOUT) :: TrackingDat

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz, krot
INTEGER, INTENT(IN) :: FastMocLv
! ----------------------------------------------------
INTEGER :: iazi, ipol, iaRay, jaRay, iAsy, iSurf, PhiAnginSvIdx, PhiAngOutSvIdx, ifsr, irsegidx, icellrayidx
INTEGER :: icRay, jbeg, jend, jinc, ihpRay, iRaySeg, iGeoTyp, iAsyTyp, jhPin, icBss, jcBss, jcRay, iReg, iCel, iRaySeg1
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
REAL, POINTER, DIMENSION(:,:)   :: ExpA, ExpB, wtang, OptLenList, PhiAngOutPolar
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
ExpA           => TrackingDat%ExpA
ExpB           => TrackingDat%ExpB
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
      ExpAppPolar(iPol, iRaySeg, icRay) = ExpA(ExpAppIdx(iRaySeg, icRay), iPol) * optlenlist(iRaySeg, icRay) + ExpB(ExpAppIdx(iRaySeg, icRay), iPol)
    END DO
  END DO
END DO
! ----------------------------------------------------
IF (krot .EQ. 1) THEN
  jbeg = 1; jend = nCoreRay; jinc = 1
ELSE
  jend = 1; jbeg = nCoreRay; jinc = -1
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
! Loc.
NULLIFY (phis)
NULLIFY (src)
NULLIFY (xst)
NULLIFY (ExpA)
NULLIFY (ExpB)
NULLIFY (Pin)
NULLIFY (FsrIdx)
NULLIFY (ExpAppIdx)
NULLIFY (wtang)
NULLIFY (OptLenList)
NULLIFY (jout)
NULLIFY (ExpAppPolar)
NULLIFY (PhiAngOutPolar)

! Hex.
NULLIFY (haRay_Loc)
NULLIFY (CelRay_Loc)
NULLIFY (hRotRay_Loc)
! ----------------------------------------------------

END SUBROUTINE HexTrackRotRayGM_OMP
! ------------------------------------------------------------------------------------------------------------