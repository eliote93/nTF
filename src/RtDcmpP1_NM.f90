#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmpP1_NM(RayInfo, CoreInfo, phisNM, phimNM, PhiAngInNM, xstNM, srcNM, srcmNM, MocJoutNM, iz, gb, ge, lJout)

USE OMP_LIB
USE PARAM,       ONLY : ZERO, ONE, RTHREE, RFIVE, RSEVEN
USE TYPEDEF,     ONLY : RayInfo_Type, Coreinfo_type, Asy_Type, Pin_Type, Cell_Type
USE Moc_Mod,     ONLY : TrackingDat, DcmpPhiAngInNg, DcmpPhiAngOutNg, DcmpColorAsy, RtDcmpPnThr_NM, DcmpGatherBndyFluxNg, DcmpScatterBndyFluxNg, DcmpLinkBndyFluxNg
USE PE_MOD,      ONLY : PE
USE CNTL,        ONLY : nTracerCntl
USE itrcntl_mod, ONLY : itrcntl
USE HexData,     ONLY : hLgc

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:,:)     :: phisNM, xstNM, srcNM
REAL, POINTER, DIMENSION(:,:,:)   :: phimNM, PhiAngInNM, srcmNM
REAL, POINTER, DIMENSION(:,:,:,:) :: MocJoutNM

INTEGER :: iz, gb, ge
LOGICAL :: lJout
! ----------------------------------------------------
TYPE (Pin_Type),  POINTER, DIMENSION(:) :: Pin
TYPE (Cell_Type), POINTER, DIMENSION(:) :: Cell

INTEGER :: ithr, nThr, ScatOd, iAsy, jAsy, iit, icolor, jcolor, ncolor, ixy, nxy, FsrIdxSt, icel, ifsr, jfsr, ig
LOGICAL :: lHex
REAL :: wttmp

INTEGER, PARAMETER :: AuxRec(2, 0:1) = [2, 1,  1, 2]
INTEGER, PARAMETER :: AuxHex(3, 0:2) = [3, 1, 2,  1, 2, 3,  2, 3, 1]
! ----------------------------------------------------

nxy   = CoreInfo%nxy
Cell => CoreInfo%CellInfo
Pin  => CoreInfo%Pin

lHex   = nTracerCntl%lHex
ScatOd = nTracerCntl%ScatOd

IF (lHex) THEN
  ncolor = 3; iit = mod(itrcntl%mocit, 3)
  
  IF (hLgc%l060) ncolor = 1
ELSE
  ncolor = 2; iit = mod(itrcntl%mocit, 2)
END IF

nthr = PE%nthread
CALL OMP_SET_NUM_THREADS(nThr)
! ----------------------------------------------------
DO ithr = 1, nThr
  TrackingDat(ithr)%srcNM => srcNM
  TrackingDat(ithr)%xstNM => xstNM
END DO

DcmpPhiAngOutNg(:, gb:ge, :, :, :) = ZERO

phisNM   (gb:ge, :) = ZERO
phimNM(:, gb:ge, :) = ZERO
IF (lJout)  MocJoutNM(:, gb:ge, :, :) = ZERO
! ----------------------------------------------------
DO icolor = 1, ncolor
  IF (lHex) THEN
    jcolor = AuxHex(icolor, iit)
    
    IF (hLgc%l060) jcolor = icolor
  ELSE
    jcolor = AuxRec(icolor, iit)
  END IF
  
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) CALL DcmpScatterBndyFluxNg(RayInfo, PhiAngInNM, DcmpPhiAngInNg)
#endif
  
  DO ithr = 1, nthr
    TrackingDat(ithr)%PhiAngInNM      => PhiAngInNM
    TrackingDat(ithr)%DcmpPhiAngInNg  => DcmpPhiAngInNg
    TrackingDat(ithr)%DcmpPhiAngOutNg => DcmpPhiAngOutNg
  END DO
  
  !$OMP PARALLEL PRIVATE(ithr, iAsy)
  ithr = 1
  !$ ithr = omp_get_thread_num()+1
  !$OMP DO SCHEDULE(GUIDED)
  DO iAsy = 1, DcmpColorAsy(0, jcolor)
    jAsy = DcmpColorAsy(iAsy, jcolor)
    
    CALL RtDcmpPnThr_NM(RayInfo, CoreInfo, TrackingDat(ithr), phisNM, phimNM, srcmNM, MocjoutNM, jAsy, iz, gb, ge, ScatOd, lJout, lHex)
  END DO
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL
  
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) CALL DcmpGatherBndyFluxNg(RayInfo, DcmpPhiAngOutNg)
#endif
  
  IF (PE%RTMASTER) CALL DcmpLinkBndyFluxNg(CoreInfo, RayInfo, PhiAngInNM, DcmpPhiAngInNg, DcmpPhiAngOutNg, gb, ge, jcolor)
END DO
! ----------------------------------------------------
!$OMP PARALLEL PRIVATE(ixy, FsrIdxSt, icel, ifsr, jfsr, ig, wttmp)
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
!$OMP END DO NOWAIT
!$OMP END PARALLEL

NULLIFY (Pin)
NULLIFY (Cell)
! ----------------------------------------------------

END SUBROUTINE RayTraceDcmpP1_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RtDcmpPnThr_NM(RayInfo, CoreInfo, TrackingLoc, phisNM, phimNM, srcmNM, MocjoutNM, jAsy, iz, gb, ge, ScatOd, lJout, lHex)

USE ALLOCS
USE PARAM,   ONLY : TRUE, ZERO
USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type, Pin_Type, Cell_Type, DcmpAsyRayInfo_Type
USE Moc_Mod, ONLY : Comp, DcmpPhiAngInNg, DcmpPhiAngOutNg, RecTrackRotRayPn_Dcmp, HexTrackRotRayPn_Dcmp
USE PE_MOD,  ONLY : PE
USE geom,    ONLY : nbd
USE CNTL,    ONLY : nTracerCntl
USE HexData, ONLY : hAsy

IMPLICIT NONE

TYPE (RayInfo_Type)     :: RayInfo
TYPE (CoreInfo_Type)    :: CoreInfo
TYPE (TrackingDat_Type) :: TrackingLoc

REAL, POINTER, DIMENSION(:,:)     :: phisNM
REAL, POINTER, DIMENSION(:,:,:)   :: phimNM, srcmNM
REAL, POINTER, DIMENSION(:,:,:,:) :: MocJoutNM

INTEGER :: iz, jAsy, gb, ge, ScatOd
LOGICAL :: lJout, lHex
! ----------------------------------------------------
TYPE (Cell_Type), POINTER, DIMENSION(:) :: Cell
TYPE (Pin_Type),  POINTER, DIMENSION(:) :: Pin

TYPE (DcmpAsyRayInfo_Type), POINTER, DIMENSION(:,:) :: DcmpAsyRay

REAL, POINTER, DIMENSION(:,:)     :: srcNM
!REAL, POINTER, DIMENSION(:,:,:,:) :: phiaNM, SrcAngNM
REAL, POINTER, DIMENSION(:,:,:,:) :: SrcAngNM1, SrcAngNM2

INTEGER, POINTER, DIMENSION(:) :: DcmpAsyRayCount

INTEGER :: nAziAng, nPolarAng, nAsyRay, FsrSt, FsrEd, PinSt, PinEd, nOd, iOd
INTEGER :: icel, ixy, iazi, ipol, iDcmpAsyRay, imRay, ig, ifsr, jfsr, iAsyRay, krot, ibd
REAL :: srctmp
! ----------------------------------------------------

Cell => CoreInfo%CellInfo
Pin  => CoreInfo%Pin

nAziAng          = RayInfo%nAziAngle
nPolarAng        = RayInfo%nPolarAngle
DcmpAsyRay      => RayInfo%DcmpAsyRay
DcmpAsyRayCount => RayInfo%DcmpAsyRayCount

IF (lHex) THEN
  PinSt = hAsy(jAsy)%PinIdxSt
  PinEd = hAsy(jAsy)%PinIdxSt + hAsy(jAsy)%nTotPin - 1
ELSE
  PinSt = CoreInfo%Asy(jAsy)%GlobalPinIdx(1)
  PinEd = CoreInfo%Asy(jAsy)%GlobalPinIdx(CoreInfo%AsyInfo(CoreInfo%Asy(jAsy)%AsyType)%nxy)
END IF

FsrSt = Pin(PinSt)%FsrIdxSt
FsrEd = Pin(PinEd)%FsrIdxSt + Cell(Pin(PinEd)%Cell(iz))%nFsr - 1

SELECT CASE (ScatOd)
  CASE (1); nod = 2
  CASE (2); nod = 5
  CASE (3); nod = 9
END SELECT
! ----------------------------------------------------
CALL dmalloc0(TrackingLoc%phisNM, gb, ge, FsrSt, FsrEd)
IF (ljout) CALL dmalloc0(TrackingLoc%JoutNM, 1, 3, gb, ge, 1, nbd, PinSt, PinEd)

!CALL dmalloc0(TrackingLoc%phiaNM,   1, nPolarAng, gb, ge, FsrSt, FsrEd, 1, 4)
!CALL dmalloc0(TrackingLoc%SrcAngNM, 1, nPolarAng, gb, ge, FsrSt, FsrEd, 1, 4)

CALL dmalloc0(TrackingLoc%phimNM, 1, nOd, gb, ge, FsrSt, FsrEd)
CALL dmalloc0(TrackingLoc%SrcAngNM1, gb, ge, 1, nPolarAng, FsrSt, FsrEd, 1, nAziAng)
CALL dmalloc0(TrackingLoc%SrcAngNM2, gb, ge, 1, nPolarAng, FsrSt, FsrEd, 1, nAziAng)

srcNM    => TrackingLoc%srcNM
!SrcAngNM => TrackingLoc%SrcAngNM
SrcAngNM1 => TrackingLoc%SrcAngNM1
SrcAngNM2 => TrackingLoc%SrcAngNM2
! ----------------------------------------------------
DO iazi = 1, nAziAng
  DO ifsr = FsrSt, FsrEd
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
! ----------------------------------------------------
IF (nTracerCntl%lHex) THEN
  DO iAsyRay = 1, DcmpAsyRayCount(jAsy)
    DO krot = 1, 2
      CALL HexTrackRotRayDcmpP1_NM(RayInfo, CoreInfo, TrackingLoc, DcmpAsyRay(iAsyRay, jAsy), ljout, iz, gb, ge, krot, ScatOd)
    END DO
  END DO
ELSE
  DO iAsyRay = 1, DcmpAsyRayCount(jAsy)
    DO krot = 1, 2
      !CALL RecTrackRotRayDcmpP1_NM(RayInfo, CoreInfo, TrackingLoc, DcmpAsyRay(iAsyRay, jAsy), ljout, iz, gb, ge, krot, ScatOd)
    END DO
  END DO
END IF
! ----------------------------------------------------
DO ifsr = FsrSt, FsrEd
  DO ig = gb, ge
    phisNM(ig, ifsr) = phisNM(ig, ifsr) + TrackingLoc%phisNM(ig, ifsr)
    
    DO iod = 1, nod
      phimNM(iod, ig, ifsr) = phimNM(iod, ig, ifsr) + TrackingLoc%phimNM(iod, ig, ifsr)
    END DO
  END DO
END DO

IF (ljout) THEN
  DO ixy = PinSt, PinEd
    DO ibd = 1, nbd
      DO ig = gb, ge
        MocjoutNM(:, ig, ibd, ixy) = MocjoutNM(:, ig, ibd, ixy) + TrackingLoc%joutNM(:, ig, ibd, ixy)
      END DO
    END DO
  END DO
END IF
! ----------------------------------------------------
DEALLOCATE (TrackingLoc%phisNM)
!DEALLOCATE (TrackingLoc%phiaNM)
!DEALLOCATE (TrackingLoc%SrcAngNM)
DEALLOCATE (TrackingLoc%phimNM)
DEALLOCATE (TrackingLoc%SrcAngNM1)
DEALLOCATE (TrackingLoc%SrcAngNM2)
IF (ljout) DEALLOCATE (TrackingLoc%JoutNM)

NULLIFY (Cell)
NULLIFY (Pin)
NULLIFY (srcNM)
!NULLIFY (phiaNM)
!NULLIFY (SrcAngNM)
NULLIFY (SrcAngNM1)
NULLIFY (SrcAngNM2)
NULLIFY (DcmpAsyRay)
NULLIFY (DcmpAsyRayCount)
! ----------------------------------------------------

END SUBROUTINE RtDcmpPnThr_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayDcmpP1_NM(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, ljout, iz, gb, ge, krot, ScatOd)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, TrackingDat_Type, DcmpAsyRayInfo_Type, AziAngleInfo_Type
USE HexType, ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexCoreRay, Type_HexRotRay
USE HexData, ONLY : hAsy, haRay, hAsyTypInfo

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: iz, gb, ge, krot, ScatOd
! ----------------------------------------------------
TYPE (Pin_Type),          POINTER, DIMENSION(:) :: Pin
TYPE (AziAngleInfo_Type), POINTER, DIMENSION(:) :: AziAng
  
TYPE (Type_HexAsyRay), POINTER :: haRay_Loc
TYPE (Type_HexCelRay), POINTER :: CelRay_Loc

REAL, POINTER, DIMENSION(:,:)       :: phisNM, srcNM, xstNM, EXPA, EXPB, wtang, hwt
REAL, POINTER, DIMENSION(:,:,:)     :: PhiAngInNM, phimNM, LocMwt
REAL, POINTER, DIMENSION(:,:,:,:)   :: joutNM, LocSrc
REAL, POINTER, DIMENSION(:,:,:,:,:) :: DcmpPhiAngInNg, DcmpPhiAngOutNg

INTEGER, POINTER, DIMENSION(:) :: AsyRayList, DirList, AziList

INTEGER :: mp(2)
INTEGER :: iazi, ipol, irotray, iasyray, iray, irayseg, idir, icel, iasy, iFSR, isurf, ig, jbeg, jend, jinc, imray, ipray
INTEGER :: ipst, iped, ipinc, isfst, isfed, isgst, isged, isginc, iAsyTyp, iGeoTyp, icBss, jhPin, jcBss, iOd, nOd
INTEGER :: nAsyRay, nPinRay, nRaySeg, FsrIdxSt, nPolarAng, PhiAnginSvIdx, ExpAppIdx

REAL :: wtazi(10), locwt(10), loccs, locsn

REAL, DIMENSION(RayInfo%nPolarAngle, gb:ge) :: PhiAngOut
REAL :: phid, tau, ExpApp, wtsurf

DATA mp /2, 1/
! ----------------------------------------------------

! Ray
nPolarAng     = RayInfo%nPolarAngle
AziAng       => RayInfo%AziAngle
PhiAngInSvIdx = RayInfo%PhiAngInSvIdx(iRotRay, krot)

! Loc.
phisNM     => TrackingDat%phisNM
srcNM      => TrackingDat%srcNM
xstNM      => TrackingDat%xstNM
joutNM     => TrackingDat%joutNM
PhiAngInNM => TrackingDat%PhiAngInNM
wtang      => TrackingDat%wtang
EXPA       => TrackingDat%EXPA
EXPB       => TrackingDat%EXPB
hwt        => TrackingDat%hwt

! Dcmp.
nAsyRay     = DcmpAsyRay%nAsyRay
iRotRay     = DcmpAsyRay%iRotRay
iAsy        = DcmpAsyRay%iAsy
iRay        = DcmpAsyRay%iRay
AsyRayList => DcmpAsyRay%AsyRayList
DirList    => DcmpAsyRay%DirList
AziList    => DcmpAsyRay%AziList

DcmpPhiAngInNg  => TrackingDat%DcmpPhiAngInNg
DcmpPhiAngOutNg => TrackingDat%DcmpPhiAngOutNg

! Geo.
Pin => CoreInfo%Pin

iAsyTyp = hAsy(iAsy)%AsyTyp
iGeoTyp = hAsy(iAsy)%GeoTyp
icBss   = hAsyTypInfo(iAsyTyp)%iBss

! Iter.
IF (DcmpAsyRay%lRotRayBeg(krot)) THEN
  PhiAngOut(1:nPolarAng, gb:ge) = PhiAngInNM    (1:nPolarAng, gb:ge, PhiAnginSvIdx)
ELSE
  PhiAngOut(1:nPolarAng, gb:ge) = DcmpPhiAngInNg(1:nPolarAng, gb:ge, krot, iRay, iAsy)
END IF

IF (krot .EQ. 1) THEN
  jbeg = 1; jend = nAsyRay; jinc = 1
ELSE
  jend = 1; jbeg = nAsyRay; jinc = -1
END IF

! P1
phimNM => TrackingDat%phimNM

SELECT CASE (ScatOd)
CASE (1); nOd = 2
CASE (2); nOd = 5
CASE (3); nOd = 9
END SELECT
! ----------------------------------------------------
DO imray = jbeg, jend, jinc
  iAsyRay = AsyRayList(imray)
  iazi    = AziList   (imray)
  idir    = DirList   (imray)
  IF (krot .EQ. 2) idir = mp(idir)
  
  DO ipol = 1, nPolarAng
    wtazi(ipol) = wtang(ipol, iazi)
    locwt(ipol) = hwt  (ipol, iazi)
  END DO
  
  loccs = AziAng(iazi)%cosv
  locsn = AziAng(iazi)%sinv
  
  haRay_Loc => haRay(iGeoTyp, icBss, iAsyRay)
  
  nPinRay = haRay_Loc%nhpRay
  
  IF (idir .EQ. 1) THEN
    ipst = 1; iped = nPinRay; ipinc = 1;  isfst = 1; isfed = 2; LocSrc => TrackingDat%SrcAngNM1; LocMwt => TrackingDat%mwt
  ELSE
    iped = 1; ipst = nPinRay; ipinc = -1; isfed = 1; isfst = 2; LocSrc => TrackingDat%SrcAngNM2; LocMwt => TrackingDat%mwt2
  END IF
  ! --------------------------------------------------
  DO ipray = ipst, iped, ipinc
    jhPin = haRay_Loc%CelRay(ipRay)%hPinIdx
    jhPin = hAsy(iAsy)%PinIdxSt + hAsyTypInfo(iAsyTyp)%PinLocIdx(iGeoTyp, jhPin) - 1
    jcBss = Pin(jhPin)%hCelGeo(iz)
    
    CelRay_Loc => haRay(iGeoTyp, jcBss, iAsyRay)%CelRay(ipRay)
    
    ! Surface : In-coming
    IF (lJout) THEN
      iSurf = CelRay_Loc%hSufIdx(isfst)
      
      DO ipol = 1, nPolarAng
        wtsurf = locwt(ipol) / abs(loccs * CelRay_Loc%hsn(isfst) - locsn * CelRay_Loc%hcs(isfst))
        
        DO ig = gb, ge
          joutNM(1, ig, isurf, jhpin) = joutNM(1, ig, isurf, jhpin) + PhiAngOut(ipol, ig) * wtazi(ipol)
          joutNM(3, ig, isurf, jhpin) = joutNM(3, ig, isurf, jhpin) + PhiAngOut(ipol, ig) * wtsurf
        END DO
      END DO
    END IF
    
    ! Iter. : FSR
    nRaySeg = CelRay_Loc%nSegRay
    
    IF (idir .EQ. 1) THEN
      isgst = 1; isged = nRaySeg; isginc = 1
    ELSE
      isged = 1; isgst = nRaySeg; isginc = -1
    END IF
    
    DO iRaySeg = isgst, isged, isginc
      iFSR = CelRay_Loc%MshIdx(iRaySeg) + Pin(jhPin)%FsrIdxSt - 1
      
      DO ig = gb, ge
        tau = -CelRay_Loc%SegLgh(iRaySeg) * xstNM(ig, iFSR) ! Optimum Length
        
        ExpAppIdx = max(INT(tau), -40000)
        ExpAppIdx = min(0, ExpAppIdx)
        
        DO ipol = 1, nPolarAng
          ExpApp = ExpA(ExpAppIdx, ipol) * tau + ExpB(ExpAppIdx, ipol)
          
          phid = (PhiAngOut(ipol, ig) - LocSrc(ig, ipol, ifsr, iazi)) * ExpApp
          
          PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
          
          phisNM(ig, iFSR) = phisNM(ig, iFSR) + wtazi(ipol) * phid
          
          DO iod = 1, nod
            phimNM(iod, ig, ifsr) = phimNM(iod, ig, ifsr) + LocMwt(iod, ipol, iazi) * phid ! NOTICE
          END DO
        END DO
      END DO
    END DO
    
    ! Surface : Out-going
    IF (lJout) THEN
      isurf = CelRay_Loc%hSufIdx(isfed) ! y : Big
      
      DO ipol = 1, nPolarAng
        wtsurf = locwt(ipol) / abs(loccs * CelRay_Loc%hsn(isfed) - locsn * CelRay_Loc%hcs(isfed))
        
        DO ig = gb, ge
          joutNM(2, ig, isurf, jhpin) = joutNM(2, ig, isurf, jhpin) + PhiAngOut(ipol, ig) * wtazi(ipol)
          joutNM(3, ig, isurf, jhpin) = joutNM(3, ig, isurf, jhpin) + PhiAngOut(ipol, ig) * wtsurf
        END DO
      END DO
    END IF
  END DO
END DO

DcmpPhiAngOutNg(1:nPolarAng, gb:ge, krot, iRay, iAsy) = PhiAngOut(1:nPolarAng, gb:ge)
! ----------------------------------------------------
! Ray
NULLIFY (AziAng)
NULLIFY (haRay_Loc)
NULLIFY (CelRay_Loc)

! Geo.
NULLIFY (Pin)

! Loc.
NULLIFY (phisNM)
NULLIFY (srcNM)
NULLIFY (xstNM)
NULLIFY (joutNM)
NULLIFY (PhiAngInNM)
NULLIFY (wtang)
NULLIFY (EXPA)
NULLIFY (EXPB)
NULLIFY (hwt)

! Dcmp.
NULLIFY (AsyRayList)
NULLIFY (DirList)
NULLIFY (AziList)
NULLIFY (DcmpPhiAngInNg)
NULLIFY (DcmpPhiAngOutNg)

! P1
NULLIFY (phimNM)
NULLIFY (LocMwt)
NULLIFY (LocSrc)
! ----------------------------------------------------

END SUBROUTINE HexTrackRotRayDcmpP1_NM
! ------------------------------------------------------------------------------------------------------------