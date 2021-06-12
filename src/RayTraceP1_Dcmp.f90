#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceP1_Dcmp(RayInfo, CoreInfo, phisNM, phimNM, PhiAngInNM, xstNM, srcNM, srcmNM, MocJoutNM, iz, gb, ge, lJout)

USE OMP_LIB
USE PARAM,       ONLY : ZERO, RED, BLACK, GREEN
USE TYPEDEF,     ONLY : RayInfo_Type, Coreinfo_type, Asy_Type, AsyInfo_Type
USE Moc_Mod,     ONLY : TrackingDat, DcmpPhiAngIn, DcmpPhiAngOut, &
                        RayTraceDcmp_OMP, RayTraceDcmp_Pn, RayTraceDcmp_LSCASMO, DcmpGatherBoundaryFlux, DcmpScatterBoundaryFlux, DcmpLinkBoundaryFlux
USE Core_mod,    ONLY : phisSlope, srcSlope
USE PE_MOD,      ONLY : PE
USE CNTL,        ONLY : nTracerCntl
USE itrcntl_mod, ONLY : itrcntl

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:,:)     :: phisNM, xstNM, srcNM
REAL, POINTER, DIMENSION(:,:,:)   :: phimNM, PhiAngInNM, srcmNM
REAL, POINTER, DIMENSION(:,:,:,:) :: MocJoutNM

INTEGER :: iz, gb, ge
LOGICAL :: lJout
! ----------------------------------------------------
TYPE (AsyInfo_Type), POINTER, DIMENSION(:) :: AsyInfo
TYPE (Asy_Type),     POINTER, DIMENSION(:) :: Asy

INTEGER :: color, startColor, endColor, colorInc, ithr, nThr, ScatOd, iAsy
! ----------------------------------------------------

AsyInfo => CoreInfo%AsyInfo
Asy     => CoreInfo%Asy

ScatOd = nTracerCntl%ScatOd

IF (.NOT. nTracerCntl%lHex) THEN
  IF (mod(itrcntl%mocit, 2) .EQ. 0) THEN
    startColor = BLACK
    endColor   = RED
    colorInc   = RED - BLACK
  ELSE
    startColor = RED
    endColor   = BLACK
    colorInc   = BLACK - RED
  END IF
ELSE
  IF (mod(itrcntl%mocit, 2) .EQ. 0) THEN
    startColor = GREEN
    endColor   = RED
    colorInc   = (RED - GREEN)/2
  ELSE
    startColor = RED
    endColor   = GREEN
    colorInc   = (GREEN - RED)/2
  END IF
END IF

nthr = PE%nthread
CALL OMP_SET_NUM_THREADS(nThr)
! ----------------------------------------------------
DO ithr = 1, nThr
  TrackingDat(ithr)%srcNM => srcNM
  TrackingDat(ithr)%xstNM => xstNM
END DO

DcmpPhiAngOut(:, gb:ge, :, :, :) = ZERO
! ----------------------------------------------------
DO color = startColor, endColor, colorInc
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) CALL DcmpScatterBoundaryFlux(RayInfo, PhiAngInNM, DcmpPhiAngIn)
#endif

  !$OMP PARALLEL PRIVATE(ithr, iAsy)
  ithr = 1
  !$ ithr = omp_get_thread_num()+1
  
  TrackingDat(ithr)%PhiAngInNM    => PhiAngInNM
  TrackingDat(ithr)%DcmpPhiAngIn  => DcmpPhiAngIn
  TrackingDat(ithr)%DcmpPhiAngOut => DcmpPhiAngOut
  !$OMP BARRIER
  !$OMP DO SCHEDULE(GUIDED)
  DO iAsy = PE%myAsyBeg, PE%myAsyEnd
    IF (Asy(iAsy)%color .NE. color) CYCLE
    
    CALL RayTraceDcmp_Pn(RayInfo, CoreInfo, phisNM, phimNM, srcmNM, MocjoutNM, iz, iAsy, gb, ge, ScatOd, lJout)
  END DO
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL
  
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) CALL DcmpGatherBoundaryFlux(RayInfo, DcmpPhiAngOut)
#endif
  
  IF (PE%RTMASTER) CALL DcmpLinkBoundaryFlux(CoreInfo, RayInfo, PhiAngInNM, DcmpPhiAngIn, DcmpPhiAngOut, gb, ge, color)
END DO

NULLIFY (Asy)
NULLIFY (AsyInfo)
! ----------------------------------------------------

END SUBROUTINE RayTraceP1_Dcmp
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmp_Pn(RayInfo, CoreInfo, phisNM, phimNM, srcmNM, joutNM, iz, iAsy, gb, ge, ScatOd, lJout)

USE ALLOCS
USE OMP_LIB
USE PARAM,   ONLY : TRUE, ZERO, RTHREE, RFIVE, RSEVEN
USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type, Pin_Type, Cell_Type, Asy_Type, AsyInfo_Type, DcmpAsyRayInfo_Type
USE Moc_Mod, ONLY : TrackingDat, Comp, mwt, AziMap, DcmpPhiAngIn, DcmpPhiAngOut, RecTrackRotRayPn_Dcmp, wtang
USE PE_MOD,  ONLY : PE
USE CNTL,    ONLY : nTracerCntl
USE HexData, ONLY : hAsy

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:,:)     :: phisNM
REAL, POINTER, DIMENSION(:,:,:)   :: phimNM, srcmNM
REAL, POINTER, DIMENSION(:,:,:,:) :: joutNM

INTEGER :: iz, iAsy, gb, ge, ScatOd
LOGICAL :: ljout
! ----------------------------------------------------
TYPE (Asy_Type),            POINTER, DIMENSION(:)   :: Asy
TYPE (AsyInfo_Type),        POINTER, DIMENSION(:)   :: AsyInfo
TYPE (Cell_Type),           POINTER, DIMENSION(:)   :: Cell
TYPE (Pin_Type),            POINTER, DIMENSION(:)   :: Pin
TYPE (DcmpAsyRayInfo_Type), POINTER, DIMENSION(:,:) :: DcmpAsyRay

REAL, POINTER, DIMENSION(:,:)     :: xstNM, srcNM
REAL, POINTER, DIMENSION(:,:,:,:) :: phianm, SrcAngnm

INTEGER, POINTER, DIMENSION(:) :: DcmpAsyRayCount

INTEGER :: nAziAng, nPolarAng, nAsyRay, nxy, ithr, FsrIdxSt, FsrIdxEnd, PinIdxSt, PinIdxEd, AziIdx
INTEGER :: icel, ipin, iazi, ipol, iDcmpAsyRay, imRay, ig, iFSR, jFSR, iAsyRay, krot
REAL :: wttemp, tempsrc
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
  nxy = AsyInfo(Asy(iAsy)%AsyType)%nxy
  
  PinIdxSt = Asy(iAsy)%GlobalPinIdx(1)
  PinIdxEd = Asy(iAsy)%GlobalPinIdx(nxy)
ELSE
  PinIdxSt = hAsy(iAsy)%PinIdxSt
  PinIdxEd = hAsy(iAsy)%PinIdxSt + hAsy(iAsy)%nTotPin - 1
END IF

FsrIdxSt  = Pin(PinIdxSt)%FsrIdxSt
FsrIdxEnd = Pin(PinIdxEd)%FsrIdxSt + Cell(Pin(PinIdxEd)%Cell(iz))%nFsr - 1
! ----------------------------------------------------
ithr = omp_get_thread_num() + 1

CALL dmalloc0(TrackingDat(ithr)%phianm,   1, nPolarAng, gb, ge, FsrIdxSt, FsrIdxEnd, 1, 4) ! Can be Huge Time-consuming
CALL dmalloc0(TrackingDat(ithr)%SrcAngnm, 1, nPolarAng, gb, ge, FsrIdxSt, FsrIdxEnd, 1, 4)

xstNM    => TrackingDat(ithr)%xstNM
srcNM    => TrackingDat(ithr)%srcNM
SrcAngnm => TrackingDat(ithr)%SrcAngnm

phisNM(   gb:ge, FsrIdxSt:FsrIdxEnd) = ZERO
phimNM(:, gb:ge, FsrIdxSt:FsrIdxEnd) = ZERO

IF (lJout) joutNM(:, gb:ge, :, PinIdxSt:PinIdxEd) = ZERO
! ----------------------------------------------------
DO iAzi = 1, nAziAng / 2
  IF (ScatOd .EQ. 1) THEN
    DO iFSR = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = srcNM(ig, iFSR)
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = srcNM(ig, iFSR)
          
          tempsrc = Comp(1, ipol, AziIdx) * srcmNM(1, ig, iFSR) + Comp(2, ipol, AziIdx) * srcmNM(2, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) - tempsrc
        END DO
      END DO
      
      AziIdx = nAziAng - iAzi + 1
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = srcNM(ig, iFSR)
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = srcNM(ig, iFSR)
          
          tempsrc = Comp(1, ipol, AziIdx) * srcmNM(1, ig, iFSR) + Comp(2, ipol, AziIdx) * srcmNM(2, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) - tempsrc
        END DO
      END DO
    END DO
  ELSE IF (ScatOd .EQ. 2) THEN 
    DO iFSR = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = srcNM(ig, iFSR)
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = srcNM(ig, iFSR)
          
          tempsrc = Comp(1, ipol, AziIdx) * srcmNM(1, ig, iFSR) + Comp(2, ipol, AziIdx) * srcmNM(2, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) - tempsrc
          
          tempsrc = Comp(3, ipol, AziIdx) * srcmNM(3, ig, iFSR) + Comp(4, ipol, AziIdx) * srcmNM(4, ig, iFSR) + Comp(5, ipol, AziIdx) * srcmNM(5, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) + tempsrc
        END DO
      END DO
      
      AziIdx = nAziAng - iAzi + 1
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = srcNM(ig, iFSR)
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = srcNM(ig, iFSR)
          
          tempsrc = Comp(1, ipol, AziIdx) * srcmNM(1, ig, iFSR) + Comp(2, ipol, AziIdx) * srcmNM(2, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) - tempsrc
          
          tempsrc = Comp(3, ipol, AziIdx) * srcmNM(3, ig, iFSR) + Comp(4, ipol, AziIdx) * srcmNM(4, ig, iFSR) + Comp(5, ipol, AziIdx) * srcmNM(5, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) + tempsrc
        END DO
      END DO
    END DO
  ELSEIF (ScatOd .EQ. 3) THEN
    DO iFSR = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = srcNM(ig, iFSR)
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = srcNM(ig, iFSR)
          
          tempsrc = Comp(1, ipol, AziIdx) * srcmNM(1, ig, iFSR) + Comp(2, ipol, AziIdx) * srcmNM(2, ig, iFSR) &
                  + Comp(6, ipol, AziIdx) * srcmNM(6, ig, iFSR) + Comp(7, ipol, AziIdx) * srcmNM(7, ig, iFSR) &
                  + Comp(8, ipol, AziIdx) * srcmNM(8, ig, iFSR) + Comp(9, ipol, AziIdx) * srcmNM(9, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) - tempsrc
          
          tempsrc = Comp(3, ipol, AziIdx) * srcmNM(3, ig, iFSR) + Comp(4, ipol, AziIdx) * srcmNM(4, ig, iFSR) + Comp(5, ipol, AziIdx) * srcmNM(5, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) + tempsrc
        END DO
      END DO
      
      AziIdx = nAziAng - iAzi + 1
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = srcNM(ig, iFSR)
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = srcNM(ig, iFSR)
          
          tempsrc = Comp(1, ipol, AziIdx) * srcmNM(1, ig, iFSR) + Comp(2, ipol, AziIdx) * srcmNM(2, ig, iFSR) &
                  + Comp(6, ipol, AziIdx) * srcmNM(6, ig, iFSR) + Comp(7, ipol, AziIdx) * srcmNM(7, ig, iFSR) &
                  + Comp(8, ipol, AziIdx) * srcmNM(8, ig, iFSR) + Comp(9, ipol, AziIdx) * srcmNM(9, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) - tempsrc
          
          tempsrc = Comp(3, ipol, AziIdx) * srcmNM(3, ig, iFSR) + Comp(4, ipol, AziIdx) * srcmNM(4, ig, iFSR) + Comp(5, ipol, AziIdx) * srcmNM(5, ig, iFSR)
          
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) = SrcAngnm(ipol, ig, iFSR, AziMap(AziIdx, 2)) + tempsrc
        END DO
      END DO
    END DO
  END IF
END DO
! ----------------------------------------------------
IF (nTracerCntl%lHex) THEN
  DO iAsyRay = 1, DcmpAsyRayCount(iAsy)
    DO krot = 1, 2
      CALL HexTrackRotRayPn_Dcmp(RayInfo, CoreInfo, TrackingDat(ithr), DcmpAsyRay(iAsyRay, iAsy), joutNM, ljout, iz, gb, ge, krot)
    END DO
  END DO
ELSE
  DO iAsyRay = 1, DcmpAsyRayCount(iAsy)
    DO krot = 1, 2
      CALL RecTrackRotRayPn_Dcmp(RayInfo, CoreInfo, TrackingDat(ithr), DcmpAsyRay(iAsyRay, iAsy), joutNM, ljout, iz, gb, ge, krot)
    END DO
  END DO
END IF

!DO iAzi = 1, nAziAng / 2
!  DO imRay = 1, DcmpAsyAziList(0, iAzi, iAsy)
!    iDcmpAsyRay = DcmpAsyAziList(imRay, iAzi, iAsy)
!    
!    CALL RecTrackRotRayPn_Dcmp(RayInfo, CoreInfo, TrackingDat(ithr), DcmpAsyRay(iDcmpAsyRay, iAsy), joutNM, lJout, iz, gb, ge)
!  END DO
!END DO
! ----------------------------------------------------
phianm => TrackingDat(ithr)%phianm

DO iAzi = 1, nAziAng / 2
  IF (ScatOd .EQ. 1) THEN
    DO iFSR = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          phisNM(ig, iFSR) = phisNM(ig, iFSR) + wtang(ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          
          phimNM(1:2, ig, iFSR) = phimNM(1:2, ig, iFSR) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
        END DO
      END DO
      
      AziIdx = nAziAng - iAzi + 1
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          phisNM(ig, iFSR) = phisNM(ig, iFSR) + wtang(ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          
          phimNM(1:2, ig, iFSR) = phimNM(1:2, ig, iFSR) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
        END DO
      END DO
    END DO
  ELSEIF (ScatOd .EQ. 2) THEN
    DO iFSR = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          phisNM(ig, iFSR) = phisNM(ig, iFSR) + wtang(ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          
          phimNM(1:2, ig, iFSR) = phimNM(1:2, ig, iFSR) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          phimNM(3:5, ig, iFSR) = phimNM(3:5, ig, iFSR) + mwt(3:5, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
        END DO
      END DO
      
      AziIdx = nAziAng - iAzi + 1
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          phisNM(ig, iFSR) = phisNM(ig, iFSR) + wtang(ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          
          phimNM(1:2, ig, iFSR) = phimNM(1:2, ig, iFSR) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          phimNM(3:5, ig, iFSR) = phimNM(3:5, ig, iFSR) + mwt(3:5, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
        END DO
      END DO
    END DO
  ELSEIF (ScatOd .EQ. 3) THEN
    DO iFSR = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          phisNM(ig, iFSR) = phisNM(ig, iFSR) + wtang(ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          
          phimNM(1:2, ig, iFSR) = phimNM(1:2, ig, iFSR) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          phimNM(3:5, ig, iFSR) = phimNM(3:5, ig, iFSR) + mwt(3:5, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          phimNM(6:9, ig, iFSR) = phimNM(6:9, ig, iFSR) + mwt(6:9, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
        END DO
      END DO
      
      AziIdx = nAziAng - iAzi + 1
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          phisNM(ig, iFSR) = phisNM(ig, iFSR) + wtang(ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          
          phimNM(1:2, ig, iFSR) = phimNM(1:2, ig, iFSR) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          phimNM(3:5, ig, iFSR) = phimNM(3:5, ig, iFSR) + mwt(3:5, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) + phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
          phimNM(6:9, ig, iFSR) = phimNM(6:9, ig, iFSR) + mwt(6:9, ipol, AziIdx) * (phianm(ipol, ig, iFSR, AziMap(AziIdx, 1)) - phianm(ipol, ig, iFSR, AziMap(AziIdx, 2)))
        END DO
      END DO
    END DO
  END IF
END DO

DEALLOCATE (TrackingDat(ithr)%phianm)
DEALLOCATE (TrackingDat(ithr)%SrcAngnm)
! ----------------------------------------------------
IF (ScatOd .EQ. 1) THEN
  DO iPin = PinIdxSt, PinIdxEd
    FsrIdxSt = Pin(iPin)%FsrIdxSt
    icel     = Pin(iPin)%Cell(iz)
    
    DO iFSR = 1, Cell(icel)%nFsr
      jFSR = FsrIdxSt + iFSR - 1
      
      DO ig = gb, ge
        wttemp = 1._8 / (xstNM(ig, jFSR) * Cell(icel)%vol(iFSR))
        
        phisNM(ig, jFSR) = phisNM(ig, jFSR) * wttemp + srcNM(ig, jFSR)
        
        phimNM(1:2, ig, jFSR) = phimNM(1:2, ig, jFSR) * wttemp + srcmNM(1:2, ig, jFSR) * RTHREE
      END DO
    END DO  
  END DO
ELSE IF (ScatOd .EQ. 2) THEN
  DO iPin = PinIdxSt, PinIdxEd
    FsrIdxSt = Pin(iPin)%FsrIdxSt
    icel     = Pin(iPin)%Cell(iz)
    
    DO iFSR = 1, Cell(icel)%nFsr
      jFSR = FsrIdxSt + iFSR - 1
      
      DO ig = gb, ge
        wttemp = 1._8 / (xstNM(ig, jFSR) * Cell(icel)%vol(iFSR))
        
        phisNM(ig, jFSR) = phisNM(ig, jFSR) * wttemp + srcNM(ig, jFSR)
        
        phimNM(1:2, ig, jFSR) = phimNM(1:2, ig, jFSR) * wttemp + srcmNM(1:2, ig, jFSR) * rthree
        phimNM(3:5, ig, jFSR) = phimNM(3:5, ig, jFSR) * wttemp + srcmNM(3:5, ig, jFSR) * rfive
      END DO
    END DO  
  END DO
ELSE IF (ScatOd .EQ. 3) THEN
  DO iPin = PinIdxSt, PinIdxEd
    FsrIdxSt = Pin(iPin)%FsrIdxSt
    icel     = Pin(iPin)%Cell(iz)
    
    DO iFSR = 1, Cell(icel)%nFsr
      jFSR = FsrIdxSt + iFSR - 1
      
      DO ig = gb, ge
        wttemp = 1._8 / (xstNM(ig, jFSR) * Cell(icel)%vol(iFSR))
        
        phisNM(ig, jFSR) = phisNM(ig, jFSR) * wttemp + srcNM(ig, jFSR)
        
        phimNM(1:2, ig, jFSR) = phimNM(1:2, ig, jFSR) * wttemp + srcmNM(1:2, ig, jFSR) * rthree
        phimNM(3:5, ig, jFSR) = phimNM(3:5, ig, jFSR) * wttemp + srcmNM(3:5, ig, jFSR) * rfive
        phimNM(6:9, ig, jFSR) = phimNM(6:9, ig, jFSR) * wttemp + srcmNM(6:9, ig, jFSR) * rseven
      END DO
    END DO  
  END DO
END IF
! ----------------------------------------------------
! Geo.
NULLIFY (Asy)
NULLIFY (AsyInfo)
NULLIFY (Cell)
NULLIFY (Pin)

! Loc.
NULLIFY (xstNM)
NULLIFY (srcNM)
NULLIFY (phianm)
NULLIFY (SrcAngnm)

! Dcmp.
NULLIFY (DcmpAsyRay)
NULLIFY (DcmpAsyRayCount)
! ----------------------------------------------------

END SUBROUTINE RayTraceDcmp_Pn
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayPn_Dcmp(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, joutNM, lJout, iz, gb, ge, krot)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, Asy_Type, Cell_Type, AsyRayInfo_type, CoreRayInfo_Type, CellRayInfo_type, TrackingDat_Type, DcmpAsyRayInfo_Type

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

REAL, POINTER, DIMENSION(:,:,:,:) :: joutNM

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
REAL, POINTER, DIMENSION(:,:,:,:)   :: SrcAngNM, phiaNM
REAL, POINTER, DIMENSION(:,:,:,:,:) :: DcmpPhiAngIn, DcmpPhiAngOut

INTEGER, POINTER, DIMENSION(:)   :: LocalFsrIdx, AsyRayList, DirList, AziList
INTEGER, POINTER, DIMENSION(:,:) :: AziMap

INTEGER :: mp(2), AziSvIdx(2)
INTEGER :: iAzi, ipol, iRotRay, iAsyRay, jAsyRay, iRay, iRaySeg, idir, jbeg, jend, jinc
INTEGER :: ipin, icel, iasy, ireg, isurf, ibcel, iceray, ig, ir, iPinRay
INTEGER :: nCoreRay, nAsyRay, nDcmpRay, nPinRay, nRaySeg, FsrIdxSt, nPolarAng, PhiAnginSvIdx, ExpAppIdx

REAL :: wt(10), wt2(10, 4)
REAL :: PhiAngOut(RayInfo%nPolarAngle, gb : ge)
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
DcmpPhiAngIn  => TrackingDat%DcmpPhiAngIn
DcmpPhiAngOut => TrackingDat%DcmpPhiAngOut
EXPA          => TrackingDat%EXPA
EXPB          => TrackingDat%EXPB
wtang         => TrackingDat%wtang
wtsurf        => TrackingDat%wtsurf
AziMap        => TrackingDat%AziMap
phiaNM        => TrackingDat%phianm
SrcAngNM      => TrackingDat%SrcAngnm

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
  PhiAngOut = DcmpPhiAngIn(1:nPolarAng, gb:ge, krot, iRay, iAsy)
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

DcmpPhiAngOut(1:nPolarAng, gb:ge, krot, iRay, iAsy) = PhiAngOut
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
NULLIFY (srcAngNM)
NULLIFY (phiaNM)
NULLIFY (xstNM)
NULLIFY (EXPA)
NULLIFY (EXPB)
NULLIFY (wtang)
NULLIFY (wtsurf)
NULLIFY (PhiAngInNM)
NULLIFY (LocalFsrIdx)

! Dcmp.
NULLIFY (DcmpPhiAngIn)
NULLIFY (DcmpPhiAngOut)
NULLIFY (AsyRayList)
NULLIFY (DirList)
NULLIFY (AziList)
NULLIFY (AziMap)
! ----------------------------------------------------

END SUBROUTINE RecTrackRotRayPn_Dcmp
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayPn_Dcmp(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, joutNM, lJout, iz, gb, ge, krot)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, AziAngleInfo_type, TrackingDat_Type, DcmpAsyRayInfo_Type
USE HexType, ONLY : Type_HexAsyRay, Type_HexCelRay
USE HexData, ONLY : hAsy, hAsyTypInfo, haRay

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

REAL, POINTER, DIMENSION(:,:,:,:) :: joutNM

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
REAL, POINTER, DIMENSION(:,:,:,:)   :: SrcAngNM, phiaNM
REAL, POINTER, DIMENSION(:,:,:,:,:) :: DcmpPhiAngIn, DcmpPhiAngOut

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
DcmpPhiAngIn  => TrackingDat%DcmpPhiAngIn
DcmpPhiAngOut => TrackingDat%DcmpPhiAngOut
EXPA          => TrackingDat%EXPA
EXPB          => TrackingDat%EXPB
wtang         => TrackingDat%wtang
hwt           => TrackingDat%hwt
AziMap        => TrackingDat%AziMap
phiaNM        => TrackingDat%phianm
SrcAngNM      => TrackingDat%SrcAngnm

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
  PhiAngOut = DcmpPhiAngIn(1:nPolarAng, gb:ge, krot, iRay, iAsy)
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

DcmpPhiAngOut(1:nPolarAng, gb:ge, krot, iRay, iAsy) = PhiAngOut(1:nPolarAng, gb:ge)
! ----------------------------------------------------
! Geo.
NULLIFY (Pin)

! Ray
NULLIFY (AziAng)
NULLIFY (haRay_Loc)
NULLIFY (CelRay_Loc)

! Loc.
NULLIFY (LenSeg)
NULLIFY (srcAngNM)
NULLIFY (phiaNM)
NULLIFY (xstNM)
NULLIFY (EXPA)
NULLIFY (EXPB)
NULLIFY (wtang)
NULLIFY (hwt)
NULLIFY (PhiAngInNM)
NULLIFY (LocalFsrIdx)

! Dcmp.
NULLIFY (DcmpPhiAngIn)
NULLIFY (DcmpPhiAngOut)
NULLIFY (AsyRayList)
NULLIFY (DirList)
NULLIFY (AziList)
NULLIFY (AziMap)
! ----------------------------------------------------

END SUBROUTINE HexTrackRotRayPn_Dcmp
! ------------------------------------------------------------------------------------------------------------