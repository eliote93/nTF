#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTrace_NM(RayInfo, CoreInfo, phisNM, PhiAngInNM, xstNM, srcNM, joutNM, iz, gb, ge, ljout)

USE OMP_LIB
USE PARAM,   ONLY : ZERO
USE TYPEDEF, ONLY : RayInfo_Type, CoreInfo_type, Pin_Type, Cell_Type
USE Moc_Mod, ONLY : TrackingDat
USE geom,    ONLY : nbd
USE PE_MOD,  ONLY : PE
USE CNTL,    ONLY : nTracerCntl

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

INTEGER :: iz, gb, ge
LOGICAL :: ljout

REAL, POINTER, DIMENSION(:,:)     :: phisNM, xstNM, srcNM
REAL, POINTER, DIMENSION(:,:,:)   :: PhiAngInNM
REAL, POINTER, DIMENSION(:,:,:,:) :: joutNM
! ----------------------------------------------------
TYPE (Cell_Type), POINTER, DIMENSION(:) :: Cell
TYPE (Pin_Type),  POINTER, DIMENSION(:) :: Pin

INTEGER :: nFsr, nxy, nthr, iRotRay, krot, ithr, FsrIdxSt, icel, ibd, ifsr, jfsr, ig, ixy
! ----------------------------------------------------

nFsr = CoreInfo%nCoreFsr
nxy  = CoreInfo%nxy

nthr = PE%nThread
CALL omp_set_num_threads(nthr)
! ----------------------------------------------------
!$OMP PARALLEL PRIVATE(ithr, krot, iRotRay)
ithr = omp_get_thread_num() + 1

TrackingDat(ithr)%phisNM = ZERO
IF (ljout) TrackingDat(ithr)%joutNM = ZERO

TrackingDat(ithr)%srcNM      => srcNM
TrackingDat(ithr)%xstNM      => xstNM
TrackingDat(ithr)%PhiAngInNM => PhiAngInNM

!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO krot = 1, 2
  DO iRotRay = 1, RayInfo%nRotRay
    IF (nTracerCntl%lHex) THEN
      CALL HexTrackRotRay_NM(RayInfo, CoreInfo, TrackingDat(ithr), ljout, iRotRay, iz, krot, gb, ge)
    ELSE
      CALL RecTrackRotRay_NM(RayInfo, CoreInfo, TrackingDat(ithr), ljout, iRotRay, iz, krot, gb, ge)
    END IF
  END DO
END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
! ----------------------------------------------------
phisNM(gb:ge, :) = ZERO
IF (ljout) joutNM(:, gb:ge, :, :) = ZERO

! Iter. is Necessary to avoid Stack Over-flow
DO ithr = 1, nthr
  DO ig = gb, ge
    DO ifsr = 1, nfsr
      phisNM(ig, ifsr) = phisNM(ig, ifsr) + TrackingDat(ithr)%phisNM(ig, ifsr)
    END DO
  END DO
  
  IF (.NOT. ljout) CYCLE
  
  DO ixy = 1, nxy
    DO ibd = 1, nbd
      DO ig = gb, ge
        joutNM(:, ig, ibd, ixy) = joutNM(:, ig, ibd, ixy) + TrackingDat(ithr)%joutNM(:, ig, ibd, ixy)
      END DO
    END DO
  END DO
END DO
! ----------------------------------------------------
Cell => CoreInfo%CellInfo
Pin  => CoreInfo%Pin

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ixy, FsrIdxSt, icel, ifsr, jfsr, ig)
!$OMP DO SCHEDULE(GUIDED)
DO ixy = 1, nxy
  FsrIdxSt = Pin(ixy)%FsrIdxSt
  icel     = Pin(ixy)%Cell(iz)
  
  DO ifsr = 1, Cell(icel)%nFsr
    jfsr = FsrIdxSt + ifsr - 1
    
    DO ig = gb, ge
      phisNM(ig, jfsr) = phisNM(ig, jfsr) / (xstNM(ig, jfsr) * Cell(icel)%vol(ifsr)) + srcNM(ig, jfsr)
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

NULLIFY (Cell)
NULLIFY (Pin)
! ----------------------------------------------------

END SUBROUTINE RayTrace_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRay_NM(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, krot, gb, ge)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, Asy_Type, PinInfo_Type, Cell_Type, AsyRayInfo_type, CoreRayInfo_Type, RotRayInfo_Type, CellRayInfo_type, TrackingDat_Type

IMPLICIT NONE

TYPE(RayInfo_Type)     :: RayInfo
TYPE(CoreInfo_Type)    :: CoreInfo
TYPE(TrackingDat_Type) :: TrackingDat

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz, krot, gb, ge

TYPE(Pin_Type),         POINTER, DIMENSION(:) :: Pin
TYPE(Asy_Type),         POINTER, DIMENSION(:) :: Asy
TYPE(PinInfo_Type),     POINTER, DIMENSION(:) :: PinInfo
TYPE(Cell_Type),        POINTER, DIMENSION(:) :: Cell
TYPE(AsyRayInfo_type),  POINTER, DIMENSION(:) :: AsyRay
TYPE(CoreRayInfo_Type), POINTER, DIMENSION(:) :: CoreRay
TYPE(RotRayInfo_Type),  POINTER, DIMENSION(:) :: RotRay

TYPE(CellRayInfo_Type), POINTER :: CellRay

INTEGER :: iazi, ipol, icray, jcray, iaray, jaray, ipray, iceray, irayseg, idir, ipin, icel, iasy, ifsr, isurf, jbeg, jend, jinc, ig, ibcel
INTEGER :: nRotRay, nCoreRay, nAsyRay, nPinRay, nRaySeg, nPolarAng, FsrIdxSt, PhiAnginSvIdx, PhiAngOutSvIdx, ExpAppIdx
INTEGER :: iast, iaed, iainc, ipst, iped, ipinc, isfst, isfed, isgst, isged, isginc

INTEGER :: mp(2)

INTEGER, POINTER, DIMENSION(:) :: LocalFsrIdx

REAL :: PhiAngOut(RayInfo%nPolarAngle, gb:ge)
REAL :: phid, tau, locsrc, ExpApp
REAL :: wt(10), wt2(10, 4)

REAL, POINTER, DIMENSION(:)       :: LenSeg
REAL, POINTER, DIMENSION(:,:)     :: phisNM, srcNM, xstNM, ExpA, ExpB, wtang
REAL, POINTER, DIMENSION(:,:,:)   :: PhiAngInNM, wtsurf
REAL, POINTER, DIMENSION(:,:,:,:) :: joutNM

DATA mp /2, 1/
! ----------------------------------------------------

AsyRay   => RayInfo%AsyRay
CoreRay  => RayInfo%CoreRay
RotRay   => RayInfo%RotRay
nPolarAng = RayInfo%nPolarAngle

Asy     => CoreInfo%Asy
Pin     => CoreInfo%Pin
PinInfo => CoreInfo%Pininfo
Cell    => CoreInfo%CellInfo

wtang      => TrackingDat%wtang
wtsurf     => TrackingDat%wtsurf
phisNM     => TrackingDat%phisNM
srcNM      => TrackingDat%srcNM
xstNM      => TrackingDat%xstNM
joutNM     => TrackingDat%joutNM
PhiAngInNM => TrackingDat%PhiAngInNM
ExpA       => TrackingDat%ExpA
ExpB       => TrackingDat%ExpB

nCoreRay = RotRay(irotRay)%nRay

PhiAngInSvIdx  = RayInfo%PhiAngInSvIdx (iRotRay, krot)
PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay, krot)

PhiAngOut(1:nPolarAng, gb:ge) = PhiAngInNM(1:nPolarAng, gb:ge, PhiAnginSvIdx)

IF (krot .EQ. 1) THEN
  jbeg = 1; jend = nCoreRay; jinc = 1
ELSE
  jend = 1; jbeg = nCoreRay; jinc = -1
END IF
! ----------------------------------------------------
DO icray = jbeg, jend, jinc
  jcray   = RotRay(iRotRay)%RayIdx(icray)
  nAsyRay = CoreRay(jcray)%nRay
  iazi    = CoreRay(jcray)%iang
  
  idir = RotRay(iRotRay)%dir(icray)
  IF (krot .eq. 2) idir = mp(idir)
  
  DO ipol = 1, nPolarAng
    wt(ipol) = wtang(ipol, iazi)
    
    IF (lJout) wt2(ipol, 1:4) = wtsurf(ipol, iazi, 1:4)
  END DO
  ! --------------------------------------------------
  IF (idir .EQ. 1) THEN
    iast = 1; iaed = nAsyRay; iainc = 1
  ELSE
    iaed = 1; iast = nAsyRay; iainc = -1
  END IF
  
  DO iaray = iast, iaed, iainc
    jaray = CoreRay(jcray)%AsyRayIdx(iaray)
    iasy  = CoreRay(jcray)%AsyIdx(iaray)
    
    IF (iAsy .EQ. 0) CYCLE
    
    nPinRay = AsyRay(jaray)%nCellRay
    
    IF (idir .EQ. 1) THEN
      ipst = 1; iped = nPinRay; ipinc = 1;  isfst = 1; isfed = 2
    ELSE
      iped = 1; ipst = nPinRay; ipinc = -1; isfed = 1; isfst = 2
    END IF
    
    DO ipray = ipst, iped, ipinc
      ipin   = AsyRay(jaray)%PinIdx(ipray)
      iceray = AsyRay(jaray)%PinRayIdx(ipray)
      ipin   = Asy(iAsy)%GlobalPinIdx(ipin)
      
      icel     = Pin(ipin)%Cell(iz)
      FsrIdxSt = Pin(ipin)%FsrIdxSt
      
      ibcel    = Cell(icel)%basecellstr
      CellRay => Cell(ibcel)%CellRay(iceray)
      
      nRaySeg      = CellRay%nSeg
      LocalFsrIdx => CellRay%LocalFsrIdx
      LenSeg      => CellRay%LenSeg
      
      ! Surf. : In-coming
      IF (lJout) THEN
        isurf = AsyRay(jaray)%PinRaySurf(isfst, ipray)
        
        DO ig = gb, ge
          DO ipol = 1, nPolarAng
            JoutNM(1, ig, isurf, ipin) = JoutNM(1, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
            JoutNM(3, ig, isurf, ipin) = JoutNM(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
          END DO
        END DO
      END IF
      
      ! Iter. : FSR
      IF (idir .EQ. 1) THEN
        isgst = 1; isged = nRaySeg; isginc = 1
      ELSE
        isged = 1; isgst = nRaySeg; isginc = -1
      END IF
      
      DO iRaySeg = isgst, isged, isginc
        ifsr = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
        
        DO ig = gb, ge
          tau = -LenSeg(iRaySeg) * xstNM(ig, ifsr)
          
          ExpAppIdx = max(INT(tau), -40000)
          ExpAppIdx = min(0, ExpAppIdx)
          
          locsrc = srcNM(ig, ifsr)
          
          DO ipol = 1, nPolarAng
            ExpApp = ExpA(ExpAppIdx, ipol) * tau + ExpB(ExpAppIdx, ipol)
            
            phid = (PhiAngOut(ipol, ig) - locsrc) * ExpApp
            
            PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
            
            phisNM(ig, ifsr) = phisNM(ig, ifsr) + wt(ipol) * phid
          END DO
        END DO
      END DO
      
      ! Surf. : Out-going
      IF (lJout) THEN
        isurf = AsyRay(jaray)%PinRaySurf(isfed, ipray)
        
        DO ig = gb, ge
          DO ipol = 1, nPolarAng
            JoutNM(2, ig, isurf, ipin) = JoutNM(2, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
            JoutNM(3, ig, isurf, ipin) = JoutNM(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
          END DO
        END DO
      END IF
    END DO
  END DO
END DO

PhiAngInNM(:, gb:ge, PhiAngOutSvIdx) = PhiAngOut
! ----------------------------------------------------
! Geo.
NULLIFY (Asy)
NULLIFY (Pin)
NULLIFY (PinInfo)
NULLIFY (Cell)

! Ray.
NULLIFY (AsyRay)
NULLIFY (CoreRay)
NULLIFY (RotRay)
NULLIFY (CellRay)

! Local
NULLIFY (LocalFsrIdx)
NULLIFY (LenSeg)
NULLIFY (phisNM)
NULLIFY (srcNM)
NULLIFY (xstNM)
NULLIFY (joutNM)
NULLIFY (PhiAngInNM)
NULLIFY (ExpA)
NULLIFY (ExpB)
NULLIFY (wtang)
NULLIFY (wtsurf)
! ----------------------------------------------------

END SUBROUTINE RecTrackRotRay_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRay_NM(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, krot, gb, ge)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, TrackingDat_Type, Pin_Type
USE HexType, ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexCoreRay, Type_HexRotRay
USE HexData, ONLY : hAsy, haRay, hcRay, hRotRay, hAsyTypInfo

IMPLICIT NONE

TYPE (RayInfo_Type)     :: RayInfo
TYPE (CoreInfo_Type)    :: CoreInfo
TYPE (TrackingDat_Type) :: TrackingDat

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz, krot, gb, ge
! ----------------------------------------------------
INTEGER :: iazi, ipol, icRay, jcRay, iaRay, jaRay, iRaySeg, ihpRay, iAsy, ifsr, iSurf, jbeg, jend, jinc, ig, iGeoTyp, iAsyTyp, jhPin, icBss, jcBss
INTEGER :: nCoreRay, nAsyRay, nPolarAng, PhiAnginSvIdx, PhiAngOutSvIdx, ExpAppIdx
INTEGER :: iast, iaed, iainc, ipst, iped, ipinc, isfst, isfed, isgst, isged, isginc

REAL :: phid, tau, locsrc, ExpApp

REAL :: wtazi(10)
REAL, DIMENSION(RayInfo%nPolarAngle, gb:ge) :: PhiAngOut

REAL, POINTER, DIMENSION(:,:)     :: phisNM, srcNM, xstNM, ExpA, ExpB, wtang
REAL, POINTER, DIMENSION(:,:,:)   :: PhiAngInNM
REAL, POINTER, DIMENSION(:,:,:,:) :: joutNM

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
phisNM     => TrackingDat%phisNM
srcNM      => TrackingDat%srcNM
xstNM      => TrackingDat%xstNM
joutNM     => TrackingDat%joutNM
PhiAngInNM => TrackingDat%PhiAngInNM
wtang      => TrackingDat%wtang
ExpA       => TrackingDat%ExpA
ExpB       => TrackingDat%ExpB

! Iter.
PhiAngOut(1:nPolarAng, gb:ge) = PhiAngInNM(1:nPolarAng, gb:ge, PhiAnginSvIdx)

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
          DO ig = gb, ge
            JoutNM(1, ig, iSurf, jhPin) = JoutNM(1, ig, isurf, jhPin) + PhiAngOut(ipol, ig) * wtazi(ipol)
          END DO
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
        
        DO ig = gb, ge
          tau = -CelRay_Loc%SegLgh(iRaySeg) * xstNM(ig, ifsr) ! Optimum Length
          
          ExpAppIdx = max(INT(tau), -40000)
          ExpAppIdx = min(0, ExpAppIdx)
          
          locsrc = srcNM(ig, ifsr)
          
          DO ipol = 1, nPolarAng
            ExpApp = ExpA(ExpAppIdx, ipol) * tau + ExpB(ExpAppIdx, ipol)
            
            phid = (PhiAngOut(ipol, ig) - locsrc) * ExpApp
            
            PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
            
            phisNM(ig, ifsr) = phisNM(ig, ifsr) + wtazi(ipol) * phid
          END DO
        END DO
      END DO
      
      ! Surface : Out-going
      IF (lJout) THEN
        isurf = CelRay_Loc%hSufIdx(isfed)
        
        DO ipol = 1, nPolarAng
          DO ig = gb, ge
            JoutNM(2, ig, iSurf, jhPin) = JoutNM(2, ig, iSurf, jhPin) + PhiAngOut(ipol, ig) * wtazi(ipol)
          END DO
        END DO
      END IF
    END DO
  END DO
END DO

PhiAngInNM(1:nPolarAng, gb:ge, PhiAngOutSvIdx) = PhiAngOut(1:nPolarAng, gb:ge)
! ----------------------------------------------------
! Loc.
NULLIFY (phisNM)
NULLIFY (srcNM)
NULLIFY (xstNM)
NULLIFY (joutNM)
NULLIFY (PhiAngInNM)
NULLIFY (wtang)
NULLIFY (ExpA)
NULLIFY (ExpB)
NULLIFY (Pin)

! Hex.
NULLIFY (hRotRay_Loc)
NULLIFY (CelRay_Loc)
NULLIFY (hRotRay_Loc)
! ----------------------------------------------------

END SUBROUTINE HexTrackRotRay_NM
! ------------------------------------------------------------------------------------------------------------