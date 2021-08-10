#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceP1_NM(RayInfo, CoreInfo, phisNM, phimNM, PhiAngInNM, xstNM, srcNM, srcmNM, joutNM, iz, gb, ge, ljout)

USE TIMER
USE OMP_LIB
USE PARAM,   ONLY : FALSE, ZERO, ONE, RTHREE, RFIVE, RSEVEN
USE TYPEDEF, ONLY : RayInfo_Type, CoreInfo_type, Pin_Type, Cell_Type
USE Moc_Mod, ONLY : TrackingDat, SrcAngNM1, SrcAngNM2, comp
USE geom,    ONLY : ng, nbd
USE PE_MOD,  ONLY : PE
USE CNTL,    ONLY : nTracerCntl

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

INTEGER :: iz, gb, ge
LOGICAL :: ljout

REAL, POINTER, DIMENSION(:,:)     :: phisNM, xstNM, srcNM
REAL, POINTER, DIMENSION(:,:,:)   :: phimNM, PhiAngInNM, srcmNM
REAL, POINTER, DIMENSION(:,:,:,:) :: joutNM
! ----------------------------------------------------
TYPE (Cell_Type), POINTER, DIMENSION(:) :: Cell
TYPE (Pin_Type),  POINTER, DIMENSION(:) :: Pin

INTEGER :: nAziAng, nPolarAng, nFsr, nxy, nthr, scatod, nod
INTEGER :: iRotRay, ipol, iazi, krot, ithr, FsrIdxSt, icel, ibd, ig, ifsr, jfsr, iod, ixy
REAL :: wttmp, srctmp
! ----------------------------------------------------

nFsr = CoreInfo%nCoreFsr
nxy  = CoreInfo%nxy

nthr = PE%nThread

ScatOd = nTracerCntl%ScatOd

nAziAng   = RayInfo%nAziAngle
nPolarAng = RayInfo%nPolarAngle

nod = 2
IF (ScatOd .EQ. 2) nod = 5
IF (ScatOd .EQ. 3) nod = 9
! ----------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ithr, iazi, Ifsr, ipol, ig, srctmp)
ithr = omp_get_thread_num() + 1

DO iazi = 1, nAziAng
  DO ifsr = PE%myOmpFsrBeg(ithr), PE%myOmpFsrEnd(ithr)
    DO ipol = 1, nPolarAng
      DO ig = gb, ge
        SrcAngNM1(ig, ipol, ifsr, iazi) = srcNM(ig, ifsr)
        SrcAngNM2(ig, ipol, ifsr, iazi) = srcNM(ig, ifsr)
        
        srctmp = comp(1, ipol, iazi) * srcmNM(1, ig, ifsr) + comp(2, ipol, iazi) * srcmNM(2, ig, ifsr)
        
        SrcAngNM1(ig, ipol, ifsr, iazi) = SrcAngNM1(ig, ipol, ifsr, iazi) + srctmp
        SrcAngNM2(ig, ipol, ifsr, iazi) = SrcAngNM2(ig, ipol, ifsr, iazi) - srctmp
        
        IF (ScatOd .LT. 2) CYCLE
        
        srctmp = comp(3, ipol, iazi) * srcmNM(3, ig, ifsr) + comp(4, ipol, iazi) * srcmNM(4, ig, ifsr) + comp(5, ipol, iazi) * srcmNM(5, ig, ifsr)
        
        SrcAngNM1(ig, ipol, ifsr, iazi) = SrcAngNM1(ig, ipol, ifsr, iazi) + srctmp
        SrcAngNM2(ig, ipol, ifsr, iazi) = SrcAngNM2(ig, ipol, ifsr, iazi) + srctmp
        
        IF (ScatOd .LT. 3) CYCLE
        
        srctmp = comp(6, ipol, iazi) * srcmNM(6, ig, ifsr) + comp(7, ipol, iazi) * srcmNM(7, ig, ifsr) + comp(8, ipol, iazi) * srcmNM(8, ig, ifsr) + comp(9, ipol, iazi) * srcmNM(9, ig, ifsr)
        
        SrcAngNM1(ig, ipol, ifsr, iazi) = SrcAngNM1(ig, ipol, ifsr, iazi) + srctmp
        SrcAngNM2(ig, ipol, ifsr, iazi) = SrcAngNM2(ig, ipol, ifsr, iazi) - srctmp
      END DO
     END DO
  END DO
END DO
!$OMP END PARALLEL
! ----------------------------------------------------
!$OMP PARALLEL PRIVATE(ithr, krot, iRotRay)
ithr = omp_get_thread_num() + 1

TrackingDat(ithr)%phisNM = ZERO
TrackingDat(ithr)%phimNM = ZERO
IF (ljout) TrackingDat(ithr)%joutNM = ZERO

TrackingDat(ithr)%srcNM      => srcNM
TrackingDat(ithr)%xstNM      => xstNM
TrackingDat(ithr)%PhiAngInNM => PhiAngInNM
TrackingDat(ithr)%SrcAngNM1  => SrcAngNM1
TrackingDat(ithr)%SrcAngNM2  => SrcAngNM2

DO krot = 1, 2
  IF (nTracerCntl%lHex) THEN
    !$OMP DO SCHEDULE(GUIDED)
    DO iRotRay = 1, RayInfo%nRotRay
      CALL HexTrackRotRayP1_NM(RayInfo, CoreInfo, TrackingDat(ithr), ljout, iRotRay, iz, krot, gb, ge, scatod)
    END DO
    !$OMP END DO NOWAIT
  ELSE
    !$OMP DO SCHEDULE(GUIDED)
    DO iRotRay = 1, RayInfo%nRotRay
      CALL RecTrackRotRayP1_NM(RayInfo, CoreInfo, TrackingDat(ithr), ljout, iRotRay, iz, krot, gb, ge, scatod)
    END DO
    !$OMP END DO NOWAIT
  END IF
END DO
!$OMP END PARALLEL
! ----------------------------------------------------
phisNM(gb:ge, :)    = ZERO
phimNM(:, gb:ge, :) = ZERO
IF (ljout) joutNM(:, gb:ge, :, :) = ZERO

DO ithr = 1, nthr
  DO ifsr = 1, nFsr
    DO ig = gb, ge
      phisNM(ig, ifsr) = phisNM(ig, ifsr) + TrackingDat(ithr)%phisNM(ig, ifsr)
      
      DO iod = 1, nod
        phimNM(iod, ig, ifsr) = phimNM(iod, ig, ifsr) + TrackingDat(ithr)%phimNM(iod, ig, ifsr)
      END DO
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

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ixy, FsrIdxSt, icel, ifsr, jfsr, ig, wttmp)
!$OMP DO SCHEDULE(GUIDED)
DO ixy = 1, nxy
  FsrIdxSt = Pin(ixy)%FsrIdxSt
  icel     = Pin(ixy)%Cell(iz)
  
  DO ifsr = 1, Cell(icel)%nFsr
    jfsr = FsrIdxSt + ifsr - 1
    
    DO ig = gb, ge
      wttmp = ONE / (xstNM(ig, jfsr) * Cell(icel)%vol(ifsr))
      
      phisNM(ig, jfsr) = phisNM(ig, jfsr) * wttmp + srcNM(ig, jfsr)
      
      phimNM(1:2, ig, jfsr) = phimNM(1:2, ig, jfsr) * wttmp + srcmNM(1:2, ig, jfsr) * RTHREE
      
      IF (ScatOd .LT. 2) CYCLE
      
      phimNM(3:5, ig, jfsr) = phimNM(3:5, ig, jfsr) * wttmp + srcmNM(3:5, ig, jfsr) * RFIVE
      
      IF (ScatOd .LT. 3) CYCLE
      
      phimNM(6:9, ig, jfsr) = phimNM(6:9, ig, jfsr) * wttmp + srcmNM(6:9, ig, jfsr) * RSEVEN
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

NULLIFY (Cell)
NULLIFY (Pin)
! ----------------------------------------------------

END SUBROUTINE RayTraceP1_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayP1_NM(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, krot, gb, ge, scatod)

USE PARAM,   ONLY : TRUE, FALSE, ZERO
USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, Asy_Type, PinInfo_Type, Cell_Type, AsyRayInfo_type, CoreRayInfo_Type, RotRayInfo_Type, CellRayInfo_type, TrackingDat_Type

IMPLICIT NONE

TYPE(RayInfo_Type)     :: RayInfo
TYPE(CoreInfo_Type)    :: CoreInfo
TYPE(TrackingDat_Type) :: TrackingDat

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz, krot, gb, ge, scatod

TYPE(Pin_Type),         POINTER, DIMENSION(:) :: Pin
TYPE(Asy_Type),         POINTER, DIMENSION(:) :: Asy
TYPE(PinInfo_Type),     POINTER, DIMENSION(:) :: PinInfo
TYPE(Cell_Type),        POINTER, DIMENSION(:) :: Cell
TYPE(AsyRayInfo_type),  POINTER, DIMENSION(:) :: AsyRay
TYPE(CoreRayInfo_Type), POINTER, DIMENSION(:) :: CoreRay
TYPE(RotRayInfo_Type),  POINTER, DIMENSION(:) :: RotRay

TYPE(CellRayInfo_Type), POINTER :: CellRay

INTEGER :: iazi, ipol, icray, jcray, iaray, jaray, ipray, iceray, irayseg, idir, ipin, icel, iasy, ifsr, isurf, jbeg, jend, jinc, ig, ibcel, iod, nod
INTEGER :: nRotRay, nCoreRay, nAsyRay, nPinRay, nRaySeg, nPolarAng, FsrIdxSt, PhiAnginSvIdx, PhiAngOutSvIdx, ExpAppIdx
INTEGER :: iast, iaed, iainc, ipst, iped, ipinc, isfst, isfed, isgst, isged, isginc

INTEGER :: mp(2)

INTEGER, POINTER, DIMENSION(:) :: LocalFsrIdx

REAL :: PhiAngOut(RayInfo%nPolarAngle, gb : ge)
REAL :: phid, tau, ExpApp

REAL :: wt(10), wt2(10, 4)

REAL, POINTER, DIMENSION(:)       :: LenSeg
REAL, POINTER, DIMENSION(:,:)     :: phisNM, xstNM, ExpA, ExpB, wtang
REAL, POINTER, DIMENSION(:,:,:)   :: phimNM, PhiAngInNM, wtsurf, mwt, mwt2
REAL, POINTER, DIMENSION(:,:,:,:) :: joutNM, SrcAngNM1, SrcAngNM2

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

phisNM      => TrackingDat%phisNM
phimNM      => TrackingDat%phimNM
xstNM       => TrackingDat%xstNM
joutNM      => TrackingDat%joutNM
PhiAngInNM  => TrackingDat%PhiAngInNM
ExpA        => TrackingDat%ExpA
ExpB        => TrackingDat%ExpB
wtang       => TrackingDat%wtang
wtsurf      => TrackingDat%wtsurf
mwt         => TrackingDat%mwt
mwt2        => TrackingDat%mwt2
SrcAngNM1   => TrackingDat%SrcAngNM1
SrcAngNM2   => TrackingDat%SrcAngNM2

nCoreRay = RotRay(irotRay)%nRay

PhiAngInSvIdx  = RayInfo%PhiAngInSvIdx (iRotRay, krot)
PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay, krot)

PhiAngOut = PhiAngInNM(:, gb:ge, PhiAnginSvIdx)

nod = 2
IF(ScatOd .EQ. 2) nod = 5
IF(ScatOd .EQ. 3) nod = 9

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
    DO iaray = 1, nAsyRay
      jaray = CoreRay(jcray)%AsyRayIdx(iaray)
      iasy  = CoreRay(jcray)%AsyIdx(iaray)
      
      IF (iAsy .EQ. 0) CYCLE
      
      nPinRay = AsyRay(jaray)%nCellRay
      
      DO ipray = 1, nPinRay
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
        
        IF (lJout) THEN
          isurf = AsyRay(jaray)%PinRaySurf(1, ipray)
          
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              JoutNM(1, ig, isurf, ipin) = JoutNM(1, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
              JoutNM(3, ig, isurf, ipin) = JoutNM(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
            END DO
          END DO
        END IF
        
        DO iRaySeg = 1, nRaySeg
          ifsr = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
          
          DO ig = gb, ge
            tau = -LenSeg(iRaySeg) * xstNM(ig, ifsr)
            
            ExpAppIdx = max(INT(tau), -40000)
            ExpAppIdx = min(0, ExpAppIdx)
            
            DO ipol = 1, nPolarAng
              ExpApp = ExpA(ExpAppIdx, ipol) * tau + ExpB(ExpAppIdx, ipol)
              phid   = (PhiAngOut(ipol, ig) - SrcAngNM1(ig, ipol, ifsr, iazi)) * ExpApp ! NOTICE : 1
              
              PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
              
              phisNM(ig, ifsr) = phisNM(ig, ifsr) + wt(ipol) * phid
              
              DO iod = 1, nod
                phimNM(iod, ig, ifsr) = phimNM(iod, ig, ifsr) + mwt(iod, ipol, iazi) * phid ! NOTICE : 1
              END DO
            END DO
          END DO
        END DO
        
        IF (lJout) THEN
          isurf = AsyRay(jaray)%PinRaySurf(2, ipray)
          
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              JoutNM(2, ig, isurf, ipin) = JoutNM(2, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
              JoutNM(3, ig, isurf, ipin) = JoutNM(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
            END DO
          END DO
        END IF
      END DO
    END DO
  ! --------------------------------------------------
  ELSE
    DO iaray = nAsyRay, 1, -1
      jaray = CoreRay(jcray)%AsyRayIdx(iaray)
      iasy  = CoreRay(jcray)%AsyIdx(iaray)
      
      IF (iAsy .EQ. 0) CYCLE
      
      nPinRay = AsyRay(jaray)%nCellRay
      
      DO ipray = nPinRay, 1, -1
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
        
        IF (lJout) THEN
          isurf = AsyRay(jaray)%PinRaySurf(2, ipray)
          
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              JoutNM(1, ig, isurf, ipin) = JoutNM(1, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
              JoutNM(3, ig, isurf, ipin) = JoutNM(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
            END DO
          END DO
        END IF
        
        DO iRaySeg = nRaySeg, 1, -1
          ifsr = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
          
          DO ig = gb, ge
            tau = -LenSeg(iRaySeg) * xstNM(ig, ifsr)
            
            ExpAppIdx = max(INT(tau), -40000)
            ExpAppIdx = min(0, ExpAppIdx)
            
            DO ipol = 1, nPolarAng
              ExpApp = ExpA(ExpAppIdx, ipol) * tau + ExpB(ExpAppIdx, ipol)
              
              phid = (PhiAngOut(ipol, ig) - SrcAngNM2(ig, ipol, ifsr, iazi)) * ExpApp ! NOTICE : 2
              
              PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
              
              phisNM(ig, ifsr) = phisNM(ig, ifsr) + wt(ipol) * phid
              
              DO iod = 1, nod
                phimNM(iod, ig, ifsr) = phimNM(iod, ig, ifsr) + mwt2(iod, ipol, iazi) * phid ! NOTICE : 2
              END DO
            END DO
          END DO
        END DO
        
        IF (lJout) THEN
          isurf = AsyRay(jaray)%PinRaySurf(1, ipray)
          
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              JoutNM(2, ig, isurf, ipin) = JoutNM(2, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
              JoutNM(3, ig, isurf, ipin) = JoutNM(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
            END DO
          END DO
        END IF
      END DO
    END DO
  END IF
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
NULLIFY (PhisNM)
NULLIFY (xstNM)
NULLIFY (joutNM)
NULLIFY (PhiAngInNM)
NULLIFY (ExpA)
NULLIFY (ExpB)
NULLIFY (wtang)
NULLIFY (wtsurf)

! P1
NULLIFY (phimNM)
NULLIFY (mwt)
NULLIFY (mwt2)
NULLIFY (SrcAngNM1)
NULLIFY (SrcAngNM2)
! ----------------------------------------------------

END SUBROUTINE RecTrackRotRayP1_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayP1_NM(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, krot, gb, ge, ScatOd)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, TrackingDat_Type, Pin_Type
USE HexType, ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexCoreRay, Type_HexRotRay
USE HexData, ONLY : hAsy, haRay, hcRay, hRotRay, hAsyTypInfo

IMPLICIT NONE

TYPE (RayInfo_Type)     :: RayInfo
TYPE (CoreInfo_Type)    :: CoreInfo
TYPE (TrackingDat_Type) :: TrackingDat

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz, krot, gb, ge, ScatOd
! ----------------------------------------------------
INTEGER :: iAzi, iPol, icRay, jcRay, iaRay, jaRay, iRaySeg, ihpRay, iAsy, ifsr, iSurf, jbeg, jend, jinc, ig, iGeoTyp, iAsyTyp, jhPin, icBss, jcBss, iod
INTEGER :: nCoreRay, nAsyRay, nPolarAng, PhiAnginSvIdx, PhiAngOutSvIdx, ExpAppIdx, nod
INTEGER :: iast, iaed, iainc, ipst, iped, ipinc, isfst, isfed, isgst, isged, isginc

REAL :: phid, tau, ExpApp

REAL :: wtazi(10)
REAL, DIMENSION(RayInfo%nPolarAngle, gb:ge) :: PhiAngOut

REAL, POINTER, DIMENSION(:,:)     :: phisNM, xstNM, ExpA, ExpB, wtang
REAL, POINTER, DIMENSION(:,:,:)   :: PhiAngInNM, phimNM, LocMwt
REAL, POINTER, DIMENSION(:,:,:,:) :: joutNM, LocSrc

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

! P1
phimNM => TrackingDat%phimNM

SELECT CASE (ScatOd)
CASE (1); nOd = 2
CASE (2); nOd = 5
CASE (3); nOd = 9
END SELECT
! ----------------------------------------------------
DO icRay = jbeg, jend, jinc
  jcRay   = hRotRay_Loc%cRayIdx(icRay)
  nAsyRay = hcRay(abs(jcRay))%nmRay
  iAzi    = hcRay(abs(jcRay))%AzmIdx
  
  IF (krot .EQ. 2) jcRay = -jcRay ! Reverse the Sweep Direction
  
  DO ipol = 1, nPolarAng
    wtazi(ipol) = wtang(ipol, iazi)
  END DO
  ! --------------------------------------------------
  IF (jcRay .GT. 0) THEN
    iast = 1; iaed = nAsyRay; iainc = 1;  LocSrc => TrackingDat%SrcAngNM1; LocMwt => TrackingDat%mwt
  ELSE
    iaed = 1; iast = nAsyRay; iainc = -1; LocSrc => TrackingDat%SrcAngNM2; LocMwt => TrackingDat%mwt2
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
          
          DO ipol = 1, nPolarAng
            ExpApp = ExpA(ExpAppIdx, ipol) * tau + ExpB(ExpAppIdx, ipol)
            
            phid = (PhiAngOut(ipol, ig) - LocSrc(ig, ipol, ifsr, iazi)) * ExpApp ! NOTICE
            
            PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
            
            phisNM(ig, ifsr) = phisNM(ig, ifsr) + wtazi(ipol) * phid
            
            DO iod = 1, nod
              phimNM(iod, ig, ifsr) = phimNM(iod, ig, ifsr) + LocMwt(iod, ipol, iazi) * phid ! NOTICE
            END DO
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
NULLIFY (xstNM)
NULLIFY (joutNM)
NULLIFY (PhiAngInNM)
NULLIFY (wtang)
NULLIFY (ExpA)
NULLIFY (ExpB)
NULLIFY (Pin)

! Hex
NULLIFY (hRotRay_Loc)
NULLIFY (CelRay_Loc)
NULLIFY (hRotRay_Loc)

! P1
NULLIFY (phimNM)
NULLIFY (LocMwt)
NULLIFY (LocSrc)
! ----------------------------------------------------

END SUBROUTINE HexTrackRotRayP1_NM
! ------------------------------------------------------------------------------------------------------------