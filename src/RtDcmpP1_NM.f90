#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmpP1_NM(RayInfo, CoreInfo, phisNM, phimNM, PhiAngInNM, xstNM, srcNM, srcmNM, MocJoutNM, iz, gb, ge, lJout)

USE OMP_LIB
USE PARAM,       ONLY : ZERO, RTHREE, RFIVE, RSEVEN
USE TYPEDEF,     ONLY : RayInfo_Type, Coreinfo_type, Asy_Type, Pin_Type, Cell_Type
USE Moc_Mod,     ONLY : TrackingDat, DcmpPhiAngInNg, DcmpPhiAngOutNg, DcmpColorAsy, RtDcmpPnThr_NM, DcmpGatherBndyFlux, DcmpScatterBndyFlux, DcmpLinkBndyFlux
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
TYPE (Asy_Type),  POINTER, DIMENSION(:) :: Asy
TYPE (Pin_Type),  POINTER, DIMENSION(:) :: Pin
TYPE (Cell_Type), POINTER, DIMENSION(:) :: Cell

INTEGER :: ithr, nThr, ScatOd, iAsy, jAsy, iit, icolor, jcolor, ncolor, ixy, nxy, FsrIdxSt, icel, iFSR, jFSR, ig
LOGICAL :: lHex
REAL :: wttmp

INTEGER, PARAMETER :: AuxRec(2, 0:1) = [2, 1,  1, 2]
INTEGER, PARAMETER :: AuxHex(3, 0:2) = [3, 1, 2,  1, 2, 3,  2, 3, 1]
! ----------------------------------------------------

nxy   = CoreInfo%nxy
Asy  => CoreInfo%Asy
Cell => CoreInfo%CellInfo
Pin  => CoreInfo%Pin

ScatOd = nTracerCntl%ScatOd
lHex   = nTracerCntl%lHex

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

phisNM(   gb:ge, :) = ZERO
phimNM(:, gb:ge, :) = ZERO
IF (lJout) MocJoutNM(:, gb:ge, :, :) = ZERO
! ----------------------------------------------------
DO icolor = 1, ncolor
  IF (.NOT. nTracerCntl%lHex) THEN
    jcolor = AuxRec(icolor, iit)
    
    IF (hLgc%l060) jcolor = icolor
  ELSE
    jcolor = AuxHex(icolor, iit)
  END IF
  
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) CALL DcmpScatterBndyFlux(RayInfo, PhiAngInNM, DcmpPhiAngInNg)
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
    
    CALL RtDcmpPnThr_NM(RayInfo, CoreInfo, phisNM, phimNM, srcmNM, MocjoutNM, iz, jAsy, gb, ge, ScatOd, lJout)
  END DO
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL
  
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) CALL DcmpGatherBndyFlux(RayInfo, DcmpPhiAngOutNg)
#endif
  
  IF (PE%RTMASTER) CALL DcmpLinkBndyFlux(CoreInfo, RayInfo, PhiAngInNM, DcmpPhiAngInNg, DcmpPhiAngOutNg, gb, ge, jcolor)
END DO
! ----------------------------------------------------
!$OMP PARALLEL PRIVATE(ixy, FsrIdxSt, icel, iFSR, jFSR, ig, wttmp)
!$OMP DO SCHEDULE(GUIDED)
DO ixy = 1, nxy
  FsrIdxSt = Pin(ixy)%FsrIdxSt
  icel     = Pin(ixy)%Cell(iz)
  
  DO iFSR = 1, Cell(icel)%nFsr
    jFSR = FsrIdxSt + iFSR - 1
    
    DO ig = gb, ge
      wttmp = 1._8 / (xstNM(ig, jFSR) * Cell(icel)%vol(iFSR))
      
      phisNM(ig, jFSR) = phisNM(ig, jFSR) * wttmp + srcNM(ig, jFSR)
      
      phimNM(1:2, ig, jFSR) = phimNM(1:2, ig, jFSR) * wttmp + srcmNM(1:2, ig, jFSR) * RTHREE
      
      IF (ScatOd .LT. 2) CYCLE
      
      phimNM(3:5, ig, jFSR) = phimNM(3:5, ig, jFSR) * wttmp + srcmNM(3:5, ig, jFSR) * RFIVE
      
      IF (ScatOd .LT. 3) CYCLE
      
      phimNM(6:9, ig, jFSR) = phimNM(6:9, ig, jFSR) * wttmp + srcmNM(6:9, ig, jFSR) * RSEVEN
    END DO
  END DO  
END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

NULLIFY (Asy)
NULLIFY (Pin)
NULLIFY (Cell)
! ----------------------------------------------------

END SUBROUTINE RayTraceDcmpP1_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RtDcmpPnThr_NM(RayInfo, CoreInfo, phisNM, phimNM, srcmNM, MocJoutNM, iz, jAsy, gb, ge, ScatOd, lJout)

USE ALLOCS
USE OMP_LIB
USE PARAM,   ONLY : TRUE, ZERO
USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type, Pin_Type, Cell_Type, Asy_Type, AsyInfo_Type, DcmpAsyRayInfo_Type
USE Moc_Mod, ONLY : TrackingDat, Comp, mwt, AziMap, DcmpPhiAngInNg, DcmpPhiAngOutNg, RecTrackRotRayPn_Dcmp, HexTrackRotRayPn_Dcmp, wtang
USE PE_MOD,  ONLY : PE
USE geom,    ONLY : nbd
USE CNTL,    ONLY : nTracerCntl
USE HexData, ONLY : hAsy

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:,:)     :: phisNM
REAL, POINTER, DIMENSION(:,:,:)   :: phimNM, srcmNM
REAL, POINTER, DIMENSION(:,:,:,:) :: MocJoutNM

INTEGER :: iz, jAsy, gb, ge, ScatOd
LOGICAL :: ljout
! ----------------------------------------------------
TYPE (Asy_Type),            POINTER, DIMENSION(:)   :: Asy
TYPE (AsyInfo_Type),        POINTER, DIMENSION(:)   :: AsyInfo
TYPE (Cell_Type),           POINTER, DIMENSION(:)   :: Cell
TYPE (Pin_Type),            POINTER, DIMENSION(:)   :: Pin
TYPE (DcmpAsyRayInfo_Type), POINTER, DIMENSION(:,:) :: DcmpAsyRay

REAL, POINTER, DIMENSION(:,:)     :: srcNM
REAL, POINTER, DIMENSION(:,:,:,:) :: phiaNM, SrcAngNM

INTEGER, POINTER, DIMENSION(:) :: DcmpAsyRayCount

INTEGER :: nAziAng, nPolarAng, nAsyRay, FsrSt, FsrEd, PinSt, PinEd, AziIdx
INTEGER :: ithr, icel, ixy, iazi, ipol, iDcmpAsyRay, imRay, ig, iFSR, jFSR, iAsyRay, krot, ibd, jAzi
REAL :: tmpsrc
! ----------------------------------------------------

Asy     => CoreInfo%Asy
AsyInfo => CoreInfo%AsyInfo
Cell    => CoreInfo%CellInfo
Pin     => CoreInfo%Pin

nAziAng          = RayInfo%nAziAngle
nPolarAng        = RayInfo%nPolarAngle
DcmpAsyRay      => RayInfo%DcmpAsyRay
DcmpAsyRayCount => RayInfo%DcmpAsyRayCount

IF (.NOT. nTracerCntl%lHex) THEN
  PinSt = Asy(jAsy)%GlobalPinIdx(1)
  PinEd = Asy(jAsy)%GlobalPinIdx(AsyInfo(Asy(jAsy)%AsyType)%nxy)
ELSE
  PinSt = hAsy(jAsy)%PinIdxSt
  PinEd = hAsy(jAsy)%PinIdxSt + hAsy(jAsy)%nTotPin - 1
END IF

FsrSt = Pin(PinSt)%FsrIdxSt
FsrEd = Pin(PinEd)%FsrIdxSt + Cell(Pin(PinEd)%Cell(iz))%nFsr - 1
! ----------------------------------------------------
ithr = omp_get_thread_num() + 1

CALL dmalloc0(TrackingDat(ithr)%phisNM, gb, ge, FsrSt, FsrEd)
IF (ljout) CALL dmalloc0(TrackingDat(ithr)%JoutNM, 1, 3, gb, ge, 1, nbd, PinSt, PinEd)

CALL dmalloc0(TrackingDat(ithr)%phiaNM,   1, nPolarAng, gb, ge, FsrSt, FsrEd, 1, 4)
CALL dmalloc0(TrackingDat(ithr)%SrcAngNM, 1, nPolarAng, gb, ge, FsrSt, FsrEd, 1, 4)

srcNM    => TrackingDat(ithr)%srcNM
SrcAngNM => TrackingDat(ithr)%SrcAngNM
! ----------------------------------------------------
DO iAzi = 1, nAziAng / 2
  ! SET : Src.
  DO jAzi = 1, 2
    AziIdx = (nAziAng - 2 * iAzi + 1) * (jAzi - 1) + iAzi
    
    DO iFSR = FsrSt, FsrEd
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngNM(ipol, ig, iFSR, AziMap(AziIdx, 1)) = srcNM(ig, iFSR)
          SrcAngNM(ipol, ig, iFSR, AziMap(AziIdx, 2)) = srcNM(ig, iFSR)
          
          tmpsrc = Comp(1, ipol, AziIdx) * srcmNM(1, ig, iFSR) + Comp(2, ipol, AziIdx) * srcmNM(2, ig, iFSR)
          
          SrcAngNM(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngNM(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tmpsrc
          SrcAngNM(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngNM(ipol, ig, iFSR, AziMap(AziIdx, 2)) - tmpsrc
          
          IF (ScatOd .LT. 2) CYCLE
          
          tmpsrc = Comp(3, ipol, AziIdx) * srcmNM(3, ig, iFSR) + Comp(4, ipol, AziIdx) * srcmNM(4, ig, iFSR) + Comp(5, ipol, AziIdx) * srcmNM(5, ig, iFSR)
          
          SrcAngNM(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngNM(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tmpsrc
          SrcAngNM(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngNM(ipol, ig, iFSR, AziMap(AziIdx, 2)) + tmpsrc
          
          IF (ScatOd .LT. 3) CYCLE
          
          tmpsrc = Comp(6, ipol, AziIdx) * srcmNM(6, ig, iFSR) + Comp(7, ipol, AziIdx) * srcmNM(7, ig, iFSR) &
                 + Comp(8, ipol, AziIdx) * srcmNM(8, ig, iFSR) + Comp(9, ipol, AziIdx) * srcmNM(9, ig, iFSR)
          
          SrcAngNM(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngNM(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tmpsrc
          SrcAngNM(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngNM(ipol, ig, iFSR, AziMap(AziIdx, 2)) - tmpsrc
        END DO
      END DO
    END DO
  END DO
END DO
! ----------------------------------------------------

! Need to Go Back to the Previous Version
! To Achieve Angle Decomposition

IF (nTracerCntl%lHex) THEN
  DO iAsyRay = 1, DcmpAsyRayCount(jAsy)
    DO krot = 1, 2
      CALL HexTrackRotRayPn_Dcmp(RayInfo, CoreInfo, TrackingDat(ithr), DcmpAsyRay(iAsyRay, jAsy), ljout, iz, gb, ge, krot)
    END DO
  END DO
ELSE
  DO iAsyRay = 1, DcmpAsyRayCount(jAsy)
    DO krot = 1, 2
      CALL RecTrackRotRayPn_Dcmp(RayInfo, CoreInfo, TrackingDat(ithr), DcmpAsyRay(iAsyRay, jAsy), ljout, iz, gb, ge, krot)
    END DO
  END DO
END IF
! ----------------------------------------------------
phiaNM => TrackingDat(ithr)%phiaNM

DO iAzi = 1, nAziAng / 2
  DO jAzi = 1, 2
    AziIdx = (nAziAng - 2 * iAzi + 1) * (jAzi - 1) + iAzi
    
    DO iFSR = FsrSt, FsrEd  
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          phisNM(ig, iFSR) = phisNM(ig, iFSR) + wtang(ipol, AziIdx) * (phiaNM(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phiaNM(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          
          phimNM(1:2, ig, iFSR) = phimNM(1:2, ig, iFSR) + mwt(1:2, ipol, AziIdx) * (phiaNM(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phiaNM(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          
          IF (ScatOd .LT. 2) CYCLE
          
          phimNM(3:5, ig, iFSR) = phimNM(3:5, ig, iFSR) + mwt(3:5, ipol, AziIdx) * (phiaNM(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phiaNM(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          
          IF (ScatOd .LT. 3) CYCLE
          
          phimNM(6:9, ig, iFSR) = phimNM(6:9, ig, iFSR) + mwt(6:9, ipol, AziIdx) * (phiaNM(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phiaNM(ipol, ig, iFSR, AziMap(AziIdx, 2)))
        END DO
      END DO
    END DO
  END DO
END DO
! ----------------------------------------------------
IF (ljout) THEN
  DO ixy = PinSt, PinEd
    DO ibd = 1, nbd
      DO ig = gb, ge
        MocjoutNM(:, ig, ibd, ixy) = MocjoutNM(:, ig, ibd, ixy) + TrackingDat(ithr)%joutNM(:, ig, ibd, ixy)
      END DO
    END DO
  END DO
END IF
! ----------------------------------------------------
DEALLOCATE (TrackingDat(ithr)%phisNM)
DEALLOCATE (TrackingDat(ithr)%phiaNM)
DEALLOCATE (TrackingDat(ithr)%SrcAngNM)
IF (ljout) DEALLOCATE (TrackingDat(ithr)%JoutNM)

NULLIFY (Asy)
NULLIFY (AsyInfo)
NULLIFY (Cell)
NULLIFY (Pin)
NULLIFY (srcNM)
NULLIFY (phiaNM)
NULLIFY (SrcAngNM)
NULLIFY (DcmpAsyRay)
NULLIFY (DcmpAsyRayCount)
! ----------------------------------------------------

END SUBROUTINE RtDcmpPnThr_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayPn_Dcmp(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, lJout, iz, gb, ge, krot)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, Asy_Type, Cell_Type, AsyRayInfo_type, CoreRayInfo_Type, CellRayInfo_type, TrackingDat_Type, DcmpAsyRayInfo_Type

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

LOGICAL :: lJout
INTEGER :: iz, gb, ge, krot
! ----------------------------------------------------
TYPE (Pin_Type),          POINTER, DIMENSION(:) :: Pin
TYPE (Asy_Type),          POINTER, DIMENSION(:) :: Asy
TYPE (Cell_Type),         POINTER, DIMENSION(:) :: Cell
TYPE (AsyRayInfo_type),   POINTER, DIMENSION(:) :: AsyRay
TYPE (CoreRayInfo_Type),  POINTER, DIMENSION(:) :: CoreRay

TYPE (CellRayInfo_Type), POINTER :: CellRay

REAL, POINTER, DIMENSION(:)         :: LenSeg
REAL, POINTER, DIMENSION(:,:)       :: xstNM, EXPA, EXPB, wtang
REAL, POINTER, DIMENSION(:,:,:)     :: PhiAngInNM, wtsurf
REAL, POINTER, DIMENSION(:,:,:,:)   :: SrcAngNM, phiaNM, JoutNM
REAL, POINTER, DIMENSION(:,:,:,:,:) :: DcmpPhiAngInNg, DcmpPhiAngOutNg

INTEGER, POINTER, DIMENSION(:)   :: LocalFsrIdx, AsyRayList, DirList, AziList
INTEGER, POINTER, DIMENSION(:,:) :: AziMap

INTEGER :: mp(2), AziSvIdx(2)
INTEGER :: iAzi, ipol, iRotRay, iAsyRay, jAsyRay, iRay, iRaySeg, idir, jbeg, jend, jinc
INTEGER :: ipin, icel, iasy, ireg, isurf, ibcel, iceray, ig, ir, iPinRay
INTEGER :: nCoreRay, nAsyRay, nDcmpRay, nPinRay, nRaySeg, FsrIdxSt, nPolarAng, PhiAnginSvIdx, ExpAppIdx

REAL :: wt(10), wt2(10, 4)
REAL :: PhiAngOut(RayInfo%nPolarAngle, gb:ge)
REAL :: phid, tau, ExpApp

DATA mp /2, 1/
! ----------------------------------------------------

! Ray
nPolarAng = RayInfo%nPolarAngle
AsyRay   => RayInfo%AsyRay
CoreRay  => RayInfo%CoreRay

! Geo.
Asy  => CoreInfo%Asy
Pin  => CoreInfo%Pin
Cell => CoreInfo%CellInfo

! Tracking Dat
xstNM         => TrackingDat%xstNM
PhiAngInNM    => TrackingDat%PhiAngInNM
DcmpPhiAngInNg  => TrackingDat%DcmpPhiAngInNg
DcmpPhiAngOutNg => TrackingDat%DcmpPhiAngOutNg
EXPA          => TrackingDat%EXPA
EXPB          => TrackingDat%EXPB
wtang         => TrackingDat%wtang
wtsurf        => TrackingDat%wtsurf
AziMap        => TrackingDat%AziMap
phiaNM        => TrackingDat%phiaNM
SrcAngNM      => TrackingDat%SrcAngNM
JoutNM        => TrackingDat%JoutNM

! Dcmp.
nAsyRay     = DcmpAsyRay%nAsyRay
iRotRay     = DcmpAsyRay%iRotRay
iAsy        = DcmpAsyRay%iAsy
iRay        = DcmpAsyRay%iRay
AsyRayList => DcmpAsyRay%AsyRayList
DirList    => DcmpAsyRay%DirList
AziList    => DcmpAsyRay%AziList
! ----------------------------------------------------
PhiAngInSvIdx = RayInfo%PhiAngInSvIdx(iRotRay, krot)

IF (DcmpAsyRay%lRotRayBeg(krot)) THEN
  PhiAngOut = PhiAngInNM(1:nPolarAng, gb:ge, PhiAnginSvIdx)
ELSE
  PhiAngOut = DcmpPhiAngInNg(1:nPolarAng, gb:ge, krot, iRay, iAsy)
END IF

IF (krot .EQ. 2) THEN
  jbeg = nAsyRay; jend = 1; jinc = -1
ELSE
  jbeg = 1; jend = nAsyRay; jinc = 1
END IF

DO iAsyRay = jbeg, jend, jinc
  jAsyRay = AsyRayList(iAsyRay)
  nPinRay = AsyRay(jAsyRay)%nCellRay
  iazi    = AziList(iAsyRay)
  idir    = DirList(iAsyRay)
  IF (krot .eq. 2) idir = mp(idir)
  
  AziSvIdx = AziMap(iAzi, :)
  
  DO ipol = 1, nPolarAng
    wt(ipol) = wtang(ipol, iazi)
    
    IF (lJout) wt2(ipol, 1:4) = wtsurf(ipol, iazi, 1:4)
  END DO
  
  IF (idir .EQ. 1) THEN
    DO iPinRay = 1, nPinRay
      ipin   = AsyRay(jAsyRay)%PinIdx   (iPinRay)
      iceray = AsyRay(jAsyRay)%PinRayIdx(iPinRay)
      
      ipin = Asy(iAsy)%GlobalPinIdx(ipin)   
      
      icel     = Pin(ipin)%Cell(iz)
      FsrIdxSt = Pin(ipin)%FsrIdxSt
      
      ibcel    = Cell(icel)%basecellstr
      CellRay => Cell(ibcel)%CellRay(iceray)
      nRaySeg  = CellRay%nSeg
      
      LocalFsrIdx => CellRay%LocalFsrIdx
      LenSeg      => CellRay%LenSeg
      
      IF (lJout) THEN
        isurf = AsyRay(jAsyRay)%PinRaySurf(1, iPinRay)
        
        DO ig = gb, ge
          DO ipol = 1, nPolarAng
            joutNM(1, ig, isurf, ipin) = joutNM(1, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
            joutNM(3, ig, isurf, ipin) = joutNM(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
          END DO
        END DO
      END IF
      
      DO iRaySeg = 1, nRaySeg
        ireg = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
        
        DO ig = gb, ge
          tau = - LenSeg(iRaySeg) * xstNM(ig, ireg)
          
          ExpAppIdx = max(INT(tau), -40000)
          ExpAppIdx = min(0, ExpAppIdx)
          
          DO ipol = 1, nPolarAng
            ExpApp = ExpA(ExpAppIdx, ipol) * tau + EXPB(ExpAppIdx, ipol)
            
            phid = (PhiAngOut(ipol, ig) - SrcAngNM(ipol, ig, ireg, AziSvIdx(idir))) * ExpApp
            
            PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
            
            phiaNM(ipol, ig, ireg, AziSvIdx(idir)) = phiaNM(ipol, ig, ireg, AziSvIdx(idir)) + phid
          END DO
        END DO
      END DO
      
      IF (lJout) THEN
        isurf = AsyRay(jAsyRay)%PinRaySurf(2, iPinRay)
        
        DO ig = gb, ge
          DO ipol = 1, nPolarAng
            joutNM(2, ig, isurf, ipin) = joutNM(2, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
            joutNM(3, ig, isurf, ipin) = joutNM(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
          END DO
        END DO
      END IF
    END DO
  ELSE
    DO iPinRay = nPinRay, 1, -1
      ipin   = AsyRay(jAsyRay)%PinIdx   (iPinRay)
      iceray = AsyRay(jAsyRay)%PinRayIdx(iPinRay)
      
      ipin = Asy(iAsy)%GlobalPinIdx(ipin)   
      
      icel     = Pin(ipin)%Cell(iz)             
      FsrIdxSt = Pin(ipin)%FsrIdxSt
      
      ibcel    = Cell(icel)%basecellstr
      CellRay => Cell(ibcel)%CellRay(iceray)
      nRaySeg  = CellRay%nSeg
      
      LocalFsrIdx => CellRay%LocalFsrIdx
      LenSeg      => CellRay%LenSeg
      
      IF (lJout) THEN
        isurf = AsyRay(jAsyRay)%PinRaySurf(2, iPinRay)
        
        DO ig = gb, ge
          DO ipol = 1, nPolarAng
            joutNM(1, ig, isurf, ipin) = joutNM(1, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
            joutNM(3, ig, isurf, ipin) = joutNM(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
          END DO
        END DO
      END IF
      
      DO iRaySeg = nRaySeg, 1, -1
        ireg = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
        
        DO ig = gb, ge
          tau = - LenSeg(iRaySeg) * xstNM(ig, ireg)
          
          ExpAppIdx = max(INT(tau), -40000)
          ExpAppIdx = min(0, ExpAppIdx)
          
          DO ipol = 1, nPolarAng
            ExpApp = EXPA(ExpAppIdx, ipol) * tau + EXPB(ExpAppIdx, ipol)
            
            phid = (PhiAngOut(ipol, ig) - SrcAngNM(ipol, ig, ireg, AziSvIdx(idir))) * ExpApp
            
            PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
            
            phiaNM(ipol, ig, ireg, AziSvIdx(idir)) = phiaNM(ipol, ig, ireg, AziSvIdx(idir)) + phid
          END DO
        END DO
      END DO
      
      IF (lJout) THEN
        isurf = AsyRay(jAsyRay)%PinRaySurf(1, iPinRay)
        
        DO ig = gb, ge
          DO ipol = 1, nPolarAng
            joutNM(2, ig, isurf, ipin) = joutNM(2, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
            joutNM(3, ig, isurf, ipin) = joutNM(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
          END DO
        END DO
      END IF
    END DO
  END IF
END DO

DcmpPhiAngOutNg(1:nPolarAng, gb:ge, krot, iRay, iAsy) = PhiAngOut
! ----------------------------------------------------
! Geo.
NULLIFY (Pin)
NULLIFY (Asy)
NULLIFY (Cell)

! Ray
NULLIFY (AsyRay)
NULLIFY (CellRay)

! Loc.
NULLIFY (LenSeg)
NULLIFY (SrcAngNM)
NULLIFY (phiaNM)
NULLIFY (xstNM)
NULLIFY (EXPA)
NULLIFY (EXPB)
NULLIFY (wtang)
NULLIFY (wtsurf)
NULLIFY (PhiAngInNM)
NULLIFY (LocalFsrIdx)
NULLIFY (JoutNM)

! Dcmp.
NULLIFY (DcmpPhiAngInNg)
NULLIFY (DcmpPhiAngOutNg)
NULLIFY (AsyRayList)
NULLIFY (DirList)
NULLIFY (AziList)
NULLIFY (AziMap)
! ----------------------------------------------------

END SUBROUTINE RecTrackRotRayPn_Dcmp
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayPn_Dcmp(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, lJout, iz, gb, ge, krot)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, AziAngleInfo_type, TrackingDat_Type, DcmpAsyRayInfo_Type
USE HexType, ONLY : Type_HexAsyRay, Type_HexCelRay
USE HexData, ONLY : hAsy, hAsyTypInfo, haRay

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

LOGICAL :: lJout
INTEGER :: iz, gb, ge, krot
! ----------------------------------------------------
TYPE (Pin_Type),          POINTER, DIMENSION(:) :: Pin
TYPE (AziAngleInfo_type), POINTER, DIMENSION(:) :: AziAng

TYPE (Type_HexAsyRay), POINTER :: haRay_Loc
TYPE (Type_HexCelRay), POINTER :: CelRay_Loc

REAL, POINTER, DIMENSION(:)         :: LenSeg
REAL, POINTER, DIMENSION(:,:)       :: xstNM, EXPA, EXPB, wtang, hwt
REAL, POINTER, DIMENSION(:,:,:)     :: PhiAngInNM
REAL, POINTER, DIMENSION(:,:,:,:)   :: SrcAngNM, phiaNM, JoutNM
REAL, POINTER, DIMENSION(:,:,:,:,:) :: DcmpPhiAngInNg, DcmpPhiAngOutNg

INTEGER, POINTER, DIMENSION(:)   :: LocalFsrIdx, AsyRayList, DirList, AziList
INTEGER, POINTER, DIMENSION(:,:) :: AziMap

INTEGER :: mp(2)
INTEGER :: iAzi, ipol, iRotRay, imRay, iAsyRay, iRay, iRaySeg, idir, jbeg, jend, jinc, ipst, iped, ipinc, isfst, isfed, isgst, isged, isginc
INTEGER :: ipin, icel, iasy, iFSR, isurf, ibcel, iceray, ig, ir, ipRay, iAsyTyp, iGeoTyp, icBss, jhPin, jcBss, jazi
INTEGER :: nAsyRay, nDcmpRay, nPinRay, nRaySeg, FsrIdxSt, nPolarAng, PhiAnginSvIdx, ExpAppIdx

REAL :: wtazi(10), locwt(10), loccs, locsn, wtsurf
REAL :: PhiAngOut(RayInfo%nPolarAngle, gb:ge)
REAL :: phid, tau, ExpApp

DATA mp /2, 1/
! ----------------------------------------------------

! Ray
nPolarAng = RayInfo%nPolarAngle
AziAng   => RayInfo%AziAngle

! Geo.
Pin  => CoreInfo%Pin

! Tracking Dat
xstNM         => TrackingDat%xstNM
PhiAngInNM    => TrackingDat%PhiAngInNM
DcmpPhiAngInNg  => TrackingDat%DcmpPhiAngInNg
DcmpPhiAngOutNg => TrackingDat%DcmpPhiAngOutNg
EXPA          => TrackingDat%EXPA
EXPB          => TrackingDat%EXPB
wtang         => TrackingDat%wtang
hwt           => TrackingDat%hwt
AziMap        => TrackingDat%AziMap
phiaNM        => TrackingDat%phiaNM
SrcAngNM      => TrackingDat%SrcAngNM
JoutNM        => TrackingDat%JoutNM

! Dcmp.
nAsyRay     = DcmpAsyRay%nAsyRay
iRotRay     = DcmpAsyRay%iRotRay
iAsy        = DcmpAsyRay%iAsy
iRay        = DcmpAsyRay%iRay
AsyRayList => DcmpAsyRay%AsyRayList
DirList    => DcmpAsyRay%DirList
AziList    => DcmpAsyRay%AziList

iAsyTyp = hAsy(iAsy)%AsyTyp
iGeoTyp = hAsy(iAsy)%GeoTyp
icBss   = hAsyTypInfo(iAsyTyp)%iBss

PhiAngInSvIdx = RayInfo%PhiAngInSvIdx(iRotRay, krot)

IF (DcmpAsyRay%lRotRayBeg(krot)) THEN
  PhiAngOut = PhiAngInNM  (1:nPolarAng, gb:ge, PhiAnginSvIdx)
ELSE
  PhiAngOut = DcmpPhiAngInNg(1:nPolarAng, gb:ge, krot, iRay, iAsy)
END IF

IF (krot .EQ. 1) THEN
  jbeg = 1; jend = nAsyRay; jinc = 1
ELSE
  jend = 1; jbeg = nAsyRay; jinc = -1
END IF
! ----------------------------------------------------
DO imRay = jbeg, jend, jinc
  iAsyRay = AsyRayList(imRay)
  iazi    = AziList   (imRay)
  idir    = DirList   (imRay)
  IF (krot .eq. 2) idir = mp(idir)
  
  jazi = AziMap(iAzi, idir)
    
  DO ipol = 1, nPolarAng
    wtazi(ipol) = wtang(ipol, iazi)
    
    IF (lJout) locwt(ipol) = hwt(ipol, iazi)
  END DO
  
  loccs = AziAng(iazi)%cosv
  locsn = AziAng(iazi)%sinv
  
  haRay_Loc => haRay(iGeoTyp, icBss, iAsyRay)
  
  nPinRay = haRay_Loc%nhpRay
  
  IF (idir .EQ. 1) THEN
    ipst = 1; iped = nPinRay; ipinc =  1; isfst = 1; isfed = 2
  ELSE
    iped = 1; ipst = nPinRay; ipinc = -1; isfed = 1; isfst = 2
  END IF
  
  DO ipRay = ipst, iped, ipinc
    jhPin = haRay_Loc%CelRay(ipRay)%hPinIdx
    jhPin = hAsy(iAsy)%PinIdxSt + hAsyTypInfo(iAsyTyp)%PinLocIdx(iGeoTyp, jhPin) - 1
    jcBss = Pin(jhPin)%hCelGeo(iz)
    
    CelRay_Loc => haRay(iGeoTyp, jcBss, iAsyRay)%CelRay(ipRay)
    
    ! Surf. : In-coming
    IF (lJout) THEN
      iSurf = CelRay_Loc%hSufIdx(isfst)
      
      DO ipol = 1, nPolarAng
        wtsurf = locwt(ipol) / abs(loccs * CelRay_Loc%hsn(isfst) - locsn * CelRay_Loc%hcs(isfst))
        
        DO ig = gb, ge
          joutNM(1, ig, isurf, ipin) = joutNM(1, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wtazi(ipol)
          joutNM(3, ig, isurf, ipin) = joutNM(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wtsurf
        END DO
      END DO
    END IF
    
    ! Iter. : FSR
    nRaySeg = CelRay_Loc%nSegRay
    
    IF (idir .EQ. 1) THEN
      isgst = 1; isged = nRaySeg; isginc =  1
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
          ExpApp = ExpA(ExpAppIdx, ipol) * tau + EXPB(ExpAppIdx, ipol)
          
          phid = (PhiAngOut(ipol, ig) - SrcAngNM(ipol, ig, iFSR, jazi)) * ExpApp
          
          PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
          
          phiaNM(ipol, ig, iFSR, jazi) = phiaNM(ipol, ig, iFSR, jazi) + phid
        END DO
      END DO
    END DO
    
    ! Surf. : Out-going
    IF (lJout) THEN
      isurf = CelRay_Loc%hSufIdx(isfed)
      
      DO ipol = 1, nPolarAng
        wtsurf = locwt(ipol) / abs(loccs * CelRay_Loc%hsn(isfst) - locsn * CelRay_Loc%hcs(isfst))
        
        DO ig = gb, ge
          joutNM(2, ig, isurf, ipin) = joutNM(2, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wtazi(ipol)
          joutNM(3, ig, isurf, ipin) = joutNM(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wtsurf
        END DO
      END DO
    END IF
  END DO
END DO

DcmpPhiAngOutNg(1:nPolarAng, gb:ge, krot, iRay, iAsy) = PhiAngOut(1:nPolarAng, gb:ge)
! ----------------------------------------------------
! Geo.
NULLIFY (Pin)

! Ray
NULLIFY (AziAng)
NULLIFY (haRay_Loc)
NULLIFY (CelRay_Loc)

! Loc.
NULLIFY (LenSeg)
NULLIFY (SrcAngNM)
NULLIFY (phiaNM)
NULLIFY (xstNM)
NULLIFY (EXPA)
NULLIFY (EXPB)
NULLIFY (wtang)
NULLIFY (hwt)
NULLIFY (PhiAngInNM)
NULLIFY (LocalFsrIdx)
NULLIFY (JoutNM)

! Dcmp.
NULLIFY (DcmpPhiAngInNg)
NULLIFY (DcmpPhiAngOutNg)
NULLIFY (AsyRayList)
NULLIFY (DirList)
NULLIFY (AziList)
NULLIFY (AziMap)
! ----------------------------------------------------

END SUBROUTINE HexTrackRotRayPn_Dcmp
! ------------------------------------------------------------------------------------------------------------