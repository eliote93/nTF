#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTrace_Dcmp(RayInfo, CoreInfo, iz, gb, ge, lJout, lSubGrp, lScat1, lLinSrcCASMO, lHybrid)

USE ALLOCS
USE OMP_LIB
USE PARAM,       ONLY : TRUE, FALSE, ZERO, RED, BLACK
USE TYPEDEF,     ONLY : RayInfo_Type, Coreinfo_type, PolarAngle_Type, AziAngleInfo_Type, Pin_Type, Asy_Type, AsyInfo_Type
USE Moc_Mod,     ONLY : nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay, Expa_p, Expb_p, ApproxExp, TrackingDat, wtang, Comp, mwt, AziMap, phisnm, phimnm, srcnm, srcmnm, &
                        xstnm, MocJoutnm, PhiAngInnm, DcmpPhiAngIn, DcmpPhiAngOut, RayTraceDcmp_NM, RayTraceDcmp_Pn, RayTraceDcmp_LSCASMO, DcmpGatherBoundaryFlux, DcmpScatterBoundaryFlux, &
                        DcmpLinkBoundaryFlux
USE Core_mod,    ONLY : phisSlope, srcSlope
USE PE_MOD,      ONLY : PE
USE GEOM,        ONLY : ng
USE CNTL,        ONLY : nTracerCntl
USE itrcntl_mod, ONLY : itrcntl

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

INTEGER :: iz, gb, ge
LOGICAL :: lJout, lSubGrp, lScat1, lLinSrcCASMO, lHybrid

TYPE (PolarAngle_Type),   POINTER, DIMENSION(:) :: PolarAng
TYPE (AziAngleInfo_Type), POINTER, DIMENSION(:) :: AziAng
TYPE (AsyInfo_Type),      POINTER, DIMENSION(:) :: AsyInfo
TYPE (Asy_Type),          POINTER, DIMENSION(:) :: Asy

LOGICAL, SAVE :: lFirst
DATA lFirst /TRUE/

INTEGER :: color, tid, nThread, nPolarAngle, nAziAngle, ScatOd, AsyType, iAsy, iAzi, ipol, od, startColor, endColor, colorInc
REAL :: wtcos, wtpolar, wtsin2, wttemp
! ----------------------------------------------------

PolarAng   => RayInfo%PolarAngle
AziAng     => RayInfo%AziAngle
nPolarAngle = RayInfo%nPolarAngle
nAziAngle   = RayInfo%nAziAngle

AsyInfo => CoreInfo%AsyInfo
Asy     => CoreInfo%Asy

nThread = PE%nThread
ScatOd  = nTracerCntl%ScatOd

IF (mod(itrcntl%mocit, 2) .EQ. 0) THEN
  startColor = BLACK
  endColor   = RED
  colorInc   = RED - BLACK
ELSE
  startColor = RED
  endColor   = BLACK
  colorInc   = BLACK - RED
END IF

CALL OMP_SET_NUM_THREADS(nThread)
! ----------------------------------------------------
IF (lFirst) THEN
  lFirst = FALSE
  
  CALL ApproxExp(RayInfo%PolarAngle, nPolarAngle)
  
  DO tid = 1, nThread
    IF (TrackingDat(tid)%lAllocNM) CYCLE
    
    TrackingDat(tid)%Expa          => Expa_p
    TrackingDat(tid)%Expb          => Expb_p
    TrackingDat(tid)%srcnm         => srcnm
    TrackingDat(tid)%xstnm         => xstnm
    TrackingDat(tid)%PhiAngInnm    => PhiAngInnm
    TrackingDat(tid)%DcmpPhiAngIn  => DcmpPhiAngIn
    TrackingDat(tid)%DcmpPhiAngOut => DcmpPhiAngOut
    
    TrackingDat(tid)%lAllocNM = TRUE
  END DO
  
  IF (lLinSrcCASMO) THEN
    DO tid = 1, nThread
      IF (TrackingDat(tid)%lAllocLinSrc) CYCLE
      
      CALL dmalloc(TrackingDat(tid)%FsrIdx,                        nMaxRaySeg, nMaxCoreRay)
      CALL dmalloc(TrackingDat(tid)%ExpAppIdxnm,               ng, nMaxRaySeg, nMaxCoreRay)
      CALL dmalloc(TrackingDat(tid)%OptLenListnm,              ng, nMaxRaySeg, nMaxCoreRay)
      CALL dmalloc(TrackingDat(tid)%ExpAppnm,     nPolarAngle, ng, nMaxRaySeg, nMaxCoreRay)
      CALL dmalloc(TrackingDat(tid)%PhiAngOutnm,  nPolarAngle, ng, nMaxRaySeg + 2)
      CALL dmalloc(TrackingDat(tid)%cmOptLen,     nPolarAngle, ng, nMaxRaySeg, nMaxCoreRay)
      CALL dmalloc(TrackingDat(tid)%cmOptLenInv,  nPolarAngle, ng, nMaxRaySeg, nMaxCoreRay)
      CALL dmalloc(TrackingDat(tid)%q0,                        ng, nMaxRaySeg, nMaxCoreRay)
      CALL dmalloc(TrackingDat(tid)%q1,           nPolarAngle, ng, nMaxRaySeg, nMaxCoreRay)
      CALL dmalloc(TrackingDat(tid)%E1,           nPolarAngle, ng, nMaxRaySeg, nMaxCoreRay)
      CALL dmalloc(TrackingDat(tid)%E3,           nPolarAngle, ng, nMaxRaySeg, nMaxCoreRay)
      CALL dmalloc(TrackingDat(tid)%R1,           nPolarAngle, ng, nMaxRaySeg, nMaxCoreRay)
      CALL dmalloc(TrackingDat(tid)%R3,           nPolarAngle, ng, nMaxRaySeg, nMaxCoreRay)
      CALL dmalloc(TrackingDat(tid)%x0,                         2, nMaxRaySeg, nMaxCoreRay)
      CALL dmalloc(TrackingDat(tid)%y0,                         2, nMaxRaySeg, nMaxCoreRay)
      
      TrackingDat(tid)%lAllocLinSrc = TRUE
    END DO
  END IF
  ! ----------------------------------------------------
  IF (lScat1) THEN
    IF (ScatOd .EQ. 1) od = 2
    IF (ScatOd .EQ. 2) od = 5
    IF (ScatOd .EQ. 3) od = 9
    
    CALL dmalloc(Comp, od, nPolarAngle, nAziAngle)
    CALL dmalloc(mwt,  od, nPolarAngle, nAziAngle)
    CALL dmalloc(wtang,    nPolarAngle, nAziAngle)
    
    DO iAzi = 1, nAziAngle / 2
      AziMap(iAzi, 1) = 1
      AziMap(iAzi, 2) = 2
      
      AziMap(nAziAngle - iAzi + 1, 1) = 3
      AziMap(nAziAngle - iAzi + 1, 2) = 4
    END DO
    
    DO ipol = 1, nPolarAngle
      wttemp = PolarAng(ipol)%weight * PolarAng(ipol)%sinv
      
      DO iazi = 1, nAziAngle
        wtang(ipol, iazi) = wttemp * AziAng(iazi)%weight * AziAng(iazi)%del
      END DO
    END DO
    
    DO ipol = 1, nPolarAngle
      wttemp = PolarAng(ipol)%sinv
      
      DO iazi = 1, nAziAngle
        Comp(1, ipol, iazi) = wttemp * AziAng(iazi)%cosv
        Comp(2, ipol, iazi) = wttemp * AziAng(iazi)%sinv
        
        mwt(1:2, ipol, iazi) = Comp(1:2, ipol, iazi) * wtang(ipol, iazi)      
      END DO
    END DO
    
    IF (ScatOd .GE. 2) THEN
      DO ipol = 1, nPolarAngle
        wttemp = PolarAng(ipol)%sinv
        wtsin2 = PolarAng(ipol)%sinv * PolarAng(ipol)%sinv
        wtcos  = PolarAng(ipol)%cosv
        
        wtpolar = 1.5_8 * PolarAng(ipol)%cosv * PolarAng(ipol)%cosv - 0.5_8    
        
        DO iazi = 1, nAziAngle
          Comp(3, ipol, iazi) = wtpolar
          Comp(4, ipol, iazi) = wtsin2 * (1._8 - 2._8 * AziAng(iazi)%sinv * AziAng(iazi)%sinv)
          Comp(5, ipol, iazi) = wtsin2 * (2._8 * AziAng(iazi)%sinv * AziAng(iazi)%cosv)
          
          mwt(3,   ipol, iazi) = Comp(3, ipol, iazi) * wtang(ipol, iazi)
          mwt(4:5, ipol, iazi) = 0.75_8 * Comp(4:5, ipol, iazi) * wtang(ipol, iazi)              
        END DO
      END DO 
    END IF
    
    IF (ScatOd .EQ. 3) THEN
      DO ipol = 1, nPolarAngle
        wttemp = PolarAng(ipol)%sinv
        
        DO iazi = 1, nAziAngle
          Comp(6, ipol, iazi) = (5._8 * PolarAng(ipol)%cosv * PolarAng(ipol)%cosv - 1._8) * wttemp * AziAng(iazi)%cosv
          Comp(7, ipol, iazi) = (5._8 * PolarAng(ipol)%cosv * PolarAng(ipol)%cosv - 1._8) * wttemp * AziAng(iazi)%sinv
          Comp(8, ipol, iazi) = (wttemp ** 3._8) * (4._8 * (AziAng(iazi)%cosv ** 3._8) - 3._8 * AziAng(iazi)%cosv)
          Comp(9, ipol, iazi) = (wttemp ** 3._8) * (-4._8 * (AziAng(iazi)%sinv ** 3._8) + 3._8 * AziAng(iazi)%sinv)
          
          mwt(6:7, ipol, iazi) = 0.375_8 * Comp(6:7, ipol, iazi) * wtang(ipol, iazi)
          mwt(8:9, ipol, iazi) = 0.625_8 * Comp(8:9, ipol, iazi) * wtang(ipol, iazi)
        END DO
      END DO
    END IF
    
    DO tid = 1, nThread
      TrackingDat(tid)%wtang  => wtang
      TrackingDat(tid)%AziMap => AziMap
    END DO
  END IF
END IF
! ----------------------------------------------------
IF (lLinSrcCASMO) THEN
  DO tid = 1, nThread
    TrackingDat(tid)%srcSlope => srcSlope(:, :, :, iz)
  END DO
END IF

DcmpPhiAngOut(:, gb:ge, :, :, :) = ZERO
! ----------------------------------------------------
DO color = startColor, endColor, colorInc
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) CALL DcmpScatterBoundaryFlux(RayInfo, PhiAngInnm, DcmpPhiAngIn)
#endif

  !$OMP PARALLEL DO PRIVATE(iAsy, AsyType) SCHEDULE(DYNAMIC)
  DO iAsy = PE%myAsyBeg, PE%myAsyEnd
    IF (Asy(iAsy)%color .NE. color) CYCLE
    
    AsyType = Asy(iAsy)%AsyType
    
    IF (.NOT.lScat1 .OR. lSubGrp) THEN
      IF (AsyInfo(AsyType)%lFuel) THEN
        IF (lLinSrcCASMO) THEN
          IF (lHybrid) THEN
            CALL RayTraceDcmp_NM(RayInfo, CoreInfo, phisnm, PhiAngInnm, xstnm, srcnm, MocJoutnm, iz, iAsy, gb, ge, ljout)
          ELSE
            CALL RayTraceDcmp_LSCASMO(RayInfo, CoreInfo, phisnm, phisSlope, PhiAngInnm, srcnm, srcSlope, xstnm, MocJoutnm, iz, iAsy, gb, ge, ljout)
          END IF
        ELSE
          CALL RayTraceDcmp_NM(RayInfo, CoreInfo, phisnm, PhiAngInnm, xstnm, srcnm, MocJoutnm, iz, iAsy, gb, ge, ljout)
        END IF
      ELSE
        IF (lLinSrcCASMO) THEN
          CALL RayTraceDcmp_LSCASMO(RayInfo, CoreInfo, phisnm, phisSlope, PhiAngInnm, srcnm, srcSlope, xstnm, MocJoutnm, iz, iAsy, gb, ge, ljout)
        ELSE
          CALL RayTraceDcmp_NM(RayInfo, CoreInfo, phisnm, PhiAngInnm, xstnm, srcnm, MocJoutnm, iz, iAsy, gb, ge, ljout)
        END IF
      END IF
    ELSE
      CALL RayTraceDcmp_Pn(RayInfo, CoreInfo, phisnm, phimnm, PhiAngInnm, xstnm, srcnm, srcmnm, MocJoutnm, iz, iAsy, gb, ge, ScatOd, lJout)
    END IF
  END DO
  !$OMP END PARALLEL DO
  
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) CALL DcmpGatherBoundaryFlux(RayInfo, DcmpPhiAngOut)
#endif

  IF (PE%RTMASTER) CALL DcmpLinkBoundaryFlux(CoreInfo, RayInfo, PhiAngInnm, DcmpPhiAngIn, DcmpPhiAngOut, gb, ge, color)
END DO
! ----------------------------------------------------

END SUBROUTINE RayTrace_Dcmp
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmp_NM(RayInfo, CoreInfo, phisnm, PhiAngInnm, xstnm, srcnm, joutnm, iz, iasy, gb, ge, ljout)

USE TIMER
USE ALLOCS
USE OMP_LIB
USE PARAM,   ONLY : ZERO
USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, AsyInfo_Type, Asy_Type, Pin_Type, Cell_Type, TrackingDat_Type, DcmpAsyRayInfo_Type
USE Moc_Mod, ONLY : nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay, Expa, Expb, TrackingDat, ApproxExp, RecTrackRotRayNM_Dcmp
USE geom,    ONLY : ng
USE PE_MOD,  ONLY : PE
USE CNTL,    ONLY : nTracerCntl

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:,:)     :: phisnm, xstnm, srcnm
REAL, POINTER, DIMENSION(:,:,:)   :: PhiAngInnm
REAL, POINTER, DIMENSION(:,:,:,:) :: joutnm

INTEGER :: iz, iasy, gb, ge
LOGICAL :: ljout

TYPE (AsyInfo_Type),        POINTER, DIMENSION(:)   :: AsyInfo
TYPE (Asy_Type),            POINTER, DIMENSION(:)   :: Asy
TYPE (Cell_Type),           POINTER, DIMENSION(:)   :: Cell
TYPE (Pin_Type),            POINTER, DIMENSION(:)   :: Pin
TYPE (DcmpAsyRayInfo_Type), POINTER, DIMENSION(:,:) :: DcmpAsyRay

INTEGER, POINTER, DIMENSION(:) :: DcmpAsyRayCount

INTEGER :: AsyType, nxy, nThread, iAsyRay, iAzi, tid, icel, ig, ipin, ifsr, jfsr, krot
INTEGER :: FsrIdxSt, FsrIdxEnd, PinIdxSt, PinIdxEnd
! ----------------------------------------------------

AsyInfo => CoreInfo%AsyInfo
Asy     => CoreInfo%Asy
Cell    => CoreInfo%CellInfo
Pin     => CoreInfo%Pin

nxy     = AsyInfo(AsyType)%nxy
nThread = PE%nThread

AsyType   = Asy(iAsy)%AsyType
PinIdxSt  = Asy(iAsy)%GlobalPinIdx(1)
PinIdxEnd = Asy(iAsy)%GlobalPinIdx(nxy)

FsrIdxSt  = Pin(PinIdxSt)%FsrIdxSt
FsrIdxEnd = Pin(PinIdxEnd)%FsrIdxSt + Cell(Pin(PinIdxEnd)%Cell(iz))%nFsr - 1

DcmpAsyRay      => RayInfo%DcmpAsyRay
DcmpAsyRayCount => RayInfo%DcmpAsyRayCount

phisnm(gb:ge, FsrIdxSt:FsrIdxEnd) = ZERO

IF (ljout) joutnm(:, gb:ge, :, PinIdxSt:PinIdxEnd) = ZERO
! ----------------------------------------------------
tid = omp_get_thread_num() + 1

IF (nTracerCntl%lHex) THEN
  DO iAsyRay = 1, DcmpAsyRayCount(iAsy)
    DO krot = 1, 2
      CALL HexTrackRotRayNM_Dcmp(RayInfo, CoreInfo, TrackingDat(tid), phisnm, joutnm, ljout, DcmpAsyRay(iAsyRay, iAsy), iz, gb, ge, krot)
    END DO
  END DO
ELSE
  DO iAsyRay = 1, DcmpAsyRayCount(iAsy)
    DO krot = 1, 2
      CALL RecTrackRotRayNM_Dcmp(RayInfo, CoreInfo, TrackingDat(tid), phisnm, joutnm, ljout, DcmpAsyRay(iAsyRay, iAsy), iz, gb, ge, krot)
    END DO
  END DO
END IF

DO ipin = PinIdxSt, PinIdxEnd
  FsrIdxSt = Pin(ipin)%FsrIdxSt
  icel     = Pin(ipin)%Cell(iz)
  
  DO ifsr = 1, Cell(icel)%nFsr
    jfsr = FsrIdxSt + ifsr - 1
    
    DO ig = gb, ge
      phisnm(ig, jfsr) = phisnm(ig, jfsr) / (xstnm(ig, jfsr) * Cell(icel)%vol(ifsr)) + srcnm(ig, jfsr)
    END DO
  END DO  
END DO
! ----------------------------------------------------

END SUBROUTINE RayTraceDcmp_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayNM_Dcmp(RayInfo, CoreInfo, TrackingDat, phis, jout, ljout, DcmpAsyRay, iz, gb, ge, krot)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, Asy_Type, PinInfo_Type, Cell_Type, AziAngleInfo_Type, PolarAngle_Type, AsyRayInfo_type, RotRayInfo_Type, &
                    CellRayInfo_type, TrackingDat_Type, DcmpAsyRayInfo_Type

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

REAL, POINTER, DIMENSION(:,:)     :: phis
REAL, POINTER, DIMENSION(:,:,:,:) :: jout

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: iz, gb, ge, krot

TYPE (Pin_Type),          POINTER, DIMENSION(:) :: Pin
TYPE (Asy_Type),          POINTER, DIMENSION(:) :: Asy
TYPE (Cell_Type),         POINTER, DIMENSION(:) :: Cell
TYPE (AziAngleInfo_Type), POINTER, DIMENSION(:) :: AziAng
TYPE (PolarAngle_Type),   POINTER, DIMENSION(:) :: PolarAng
TYPE (AsyRayInfo_type),   POINTER, DIMENSION(:) :: AsyRay

TYPE (CellRayInfo_Type),  POINTER :: CellRay

REAL, POINTER, DIMENSION(:)         :: LenSeg
REAL, POINTER, DIMENSION(:,:)       :: src, xst, EXPA, EXPB
REAL, POINTER, DIMENSION(:,:,:)     :: PhiAngIn
REAL, POINTER, DIMENSION(:,:,:,:,:) :: DcmpPhiAngIn, DcmpPhiAngOut

INTEGER, POINTER, DIMENSION(:) :: LocalFsrIdx, AsyRayList, DirList, AziList

INTEGER :: mp(2)
INTEGER :: iazi, ipol, irotray, iasyray, iceray, iray, irayseg, idir, ipin, icel, ibcel, iasy, ireg, isurf, ig, jbeg, jend, jinc, imray, ipray
INTEGER :: ipst, iped, ipinc, isfst, isfed, isgst, isged, isginc
INTEGER :: nAsyRay, nPinRay, nRaySeg, FsrIdxSt, nAziAng, nPolarAng, PhiAnginSvIdx, ExpAppIdx

REAL :: wtang(10, 100), wtang2(10, 100, 4), wt(10), wt2(10, 4)

REAL, DIMENSION(RayInfo%nPolarAngle, gb:ge) :: PhiAngOut
REAL :: wttemp, phid, tau, ExpApp

DATA mp /2, 1/
! ----------------------------------------------------

! Ray Info.
AziAng   => RayInfo%AziAngle
PolarAng => RayInfo%PolarAngle
AsyRay   => RayInfo%AsyRay
nPolarAng = RayInfo%nPolarAngle
nAziAng   = RayInfo%nAziAngle

! Geo.
Asy  => CoreInfo%Asy
Pin  => CoreInfo%Pin
Cell => CoreInfo%CellInfo

! Tracking Dat Pointing
src           => TrackingDat%srcnm
xst           => TrackingDat%xstnm
PhiAngIn      => TrackingDat%phiAngInnm
DcmpPhiAngIn  => TrackingDat%DcmpPhiAngIn
DcmpPhiAngOut => TrackingDat%DcmpPhiAngOut
EXPA          => TrackingDat%EXPA
EXPB          => TrackingDat%EXPB

! Dcmp. Ray
nAsyRay     = DcmpAsyRay%nAsyRay
iRotRay     = DcmpAsyRay%iRotRay
AsyRayList => DcmpAsyRay%AsyRayList
DirList    => DcmpAsyRay%DirList
AziList    => DcmpAsyRay%AziList
iAsy        = DcmpAsyRay%iAsy
iRay        = DcmpAsyRay%iRay

! Wgh.
DO ipol = 1, nPolarAng
  wttemp = PolarAng(ipol)%weight * PolarAng(ipol)%sinv
  
  DO iazi = 1, nAziAng
    wtang(ipol, iazi) = wttemp * AziAng(iazi)%weight * AziAng(iazi)%del
  END DO
END DO

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

! Ray. B.C.
PhiAngInSvIdx = RayInfo%PhiAngInSvIdx(iRotRay, krot)

IF (DcmpAsyRay%lRotRayBeg(krot)) THEN
  PhiAngOut = PhiAngIn(1:nPolarAng, gb:ge, PhiAnginSvIdx)
ELSE
  PhiAngOut = DcmpPhiAngIn(1:nPolarAng, gb:ge, krot, iRay, iAsy)
END IF
  
IF (krot .EQ. 2) THEN
  jbeg = nAsyRay; jend = 1; jinc = -1
ELSE
  jbeg = 1; jend = nAsyRay; jinc = 1
END IF
! ----------------------------------------------------  
DO imray = jbeg, jend, jinc
  iAsyRay = AsyRayList(imray)
  nPinRay = AsyRay(iAsyRay)%nCellRay
  iazi    = AziList(imray)
  idir    = DirList(imray)
  IF (krot .EQ. 2) idir = mp(idir)
  
  DO ipol = 1, nPolarAng
    wt(ipol) = wtang(ipol, iazi)
    
    IF (lJout) wt2(ipol, 1:4) = wtang2(ipol, iazi, 1:4)
  END DO
  
  IF (idir .EQ. 1) THEN
    ipst = 1; iped = nPinRay; ipinc = 1;  isfst = 1; isfed = 2;
  ELSE
    ipst = nPinRay; iped = 1; ipinc = -1; isfst = 2; isfed = 1;
  END IF
  
  DO ipray = ipst, iped, ipinc
    ipin   = AsyRay(iAsyRay)%PinIdx(ipray)
    iceray = AsyRay(iAsyRay)%PinRayIdx(ipray)
    
    ipin = Asy(iAsy)%GlobalPinIdx(ipin)   
    
    icel     = Pin(ipin)%Cell(iz)
    FsrIdxSt = Pin(ipin)%FsrIdxSt
    
    ibcel    = Cell(icel)%basecellstr
    CellRay => Cell(ibcel)%CellRay(iceray)
    
    nRaySeg      = CellRay%nSeg
    LocalFsrIdx => CellRay%LocalFsrIdx
    LenSeg      => CellRay%LenSeg
    
    IF (lJout) THEN
      isurf = AsyRay(iAsyRay)%PinRaySurf(isfst, ipray)
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          Jout(1, ig, isurf, ipin) = Jout(1, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
          Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
        END DO
      END DO
    END IF
    
    IF (idir .EQ. 1) THEN
      isgst = 1; isged = nRaySeg; isginc = 1
    ELSE
      isgst = nRaySeg; isged = 1; isginc = -1
    END IF
    
    DO iRaySeg = isgst, isged, isginc
      ireg = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
      
      DO ig = gb, ge
        tau = - LenSeg(iRaySeg) * xst(ig, ireg)
        
        ExpAppIdx = max(INT(tau), -40000)
        ExpAppIdx = min(0, ExpAppIdx)
        
        DO ipol = 1, nPolarAng
          ExpApp = EXPA(ipol, ExpAppIdx) * tau + EXPB(ipol, ExpAppIdx)
          
          phid = (PhiAngOut(ipol, ig) - src(ig, ireg)) * ExpApp
          
          PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
          
          phis(ig, ireg) = phis(ig, ireg) + wt(ipol) * phid
        END DO
      END DO
    END DO
    
    IF (lJout) THEN
      isurf = AsyRay(iAsyRay)%PinRaySurf(isfed, ipray)
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          Jout(2, ig, isurf, ipin) = Jout(2, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
          Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
        END DO
      END DO
    END IF
  END DO
END DO

DcmpPhiAngOut(1:nPolarAng, gb:ge, krot, iRay, iAsy) = PhiAngOut
! ----------------------------------------------------

END SUBROUTINE RecTrackRotRayNM_Dcmp
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayNM_Dcmp(RayInfo, CoreInfo, TrackingDat, phis, jout, ljout, DcmpAsyRay, iz, gb, ge, krot)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, PinInfo_Type, AziAngleInfo_Type, PolarAngle_Type, RotRayInfo_Type, TrackingDat_Type, DcmpAsyRayInfo_Type
USE HexData, ONLY : hAsy
USE HexType, ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexCoreRay, Type_HexRotRay
USE HexData, ONLY : haRay, hAsyTypInfo

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

REAL, POINTER, DIMENSION(:,:)     :: phis
REAL, POINTER, DIMENSION(:,:,:,:) :: jout

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: iz, gb, ge, krot

TYPE (Pin_Type),          POINTER, DIMENSION(:) :: Pin
TYPE (AziAngleInfo_Type), POINTER, DIMENSION(:) :: AziAng
TYPE (PolarAngle_Type),   POINTER, DIMENSION(:) :: PolarAng

TYPE (Type_HexAsyRay), POINTER :: haRay_Loc
TYPE (Type_HexCelRay), POINTER :: CelRay_Loc

REAL, POINTER, DIMENSION(:,:)       :: src, xst, EXPA, EXPB
REAL, POINTER, DIMENSION(:,:,:)     :: PhiAngIn
REAL, POINTER, DIMENSION(:,:,:,:,:) :: DcmpPhiAngIn, DcmpPhiAngOut

INTEGER, POINTER, DIMENSION(:) :: AsyRayList, DirList, AziList

INTEGER :: mp(2)
INTEGER :: iazi, ipol, irotray, iasyray, iray, irayseg, idir, ipin, icel, iasy, ireg, isurf, ig, jbeg, jend, jinc, imray, ipray
INTEGER :: ipst, iped, ipinc, isfst, isfed, isgst, isged, isginc, iAsyTyp, iGeoTyp, icBss, jhPin, jcBss
INTEGER :: nAsyRay, nPinRay, nRaySeg, FsrIdxSt, nAziAng, nPolarAng, PhiAnginSvIdx, ExpAppIdx

REAL :: wtang(10, 100), wtang2(10, 100, 4), wt(10), wt2(10, 4)

REAL, DIMENSION(RayInfo%nPolarAngle, gb:ge) :: PhiAngOut
REAL :: wttemp, phid, tau, ExpApp

DATA mp /2, 1/
! ----------------------------------------------------

! Ray Info.
AziAng   => RayInfo%AziAngle
PolarAng => RayInfo%PolarAngle
nPolarAng = RayInfo%nPolarAngle
nAziAng   = RayInfo%nAziAngle

! Geo.
Pin => CoreInfo%Pin

! Tracking Dat Pointing
src           => TrackingDat%srcnm
xst           => TrackingDat%xstnm
PhiAngIn      => TrackingDat%phiAngInnm
DcmpPhiAngIn  => TrackingDat%DcmpPhiAngIn
DcmpPhiAngOut => TrackingDat%DcmpPhiAngOut
EXPA          => TrackingDat%EXPA
EXPB          => TrackingDat%EXPB

! Dcmp. Ray
iRotRay     = DcmpAsyRay%iRotRay
AsyRayList => DcmpAsyRay%AsyRayList
DirList    => DcmpAsyRay%DirList
AziList    => DcmpAsyRay%AziList
iAsy        = DcmpAsyRay%iAsy
iRay        = DcmpAsyRay%iRay

! Wgh.
DO ipol = 1, nPolarAng
  wttemp = PolarAng(ipol)%weight * PolarAng(ipol)%sinv
  
  DO iazi = 1, nAziAng
    wtang(ipol, iazi) = wttemp * AziAng(iazi)%weight * AziAng(iazi)%del
  END DO
END DO

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

! Ray. B.C.
PhiAngInSvIdx = RayInfo%PhiAngInSvIdx(iRotRay, krot)

IF (DcmpAsyRay%lRotRayBeg(krot)) THEN
  PhiAngOut = PhiAngIn(1:nPolarAng, gb:ge, PhiAnginSvIdx)
ELSE
  PhiAngOut = DcmpPhiAngIn(1:nPolarAng, gb:ge, krot, iRay, iAsy)
END IF
  
IF (krot .EQ. 2) THEN
  jbeg = nAsyRay; jend = 1; jinc = -1
ELSE
  jbeg = 1; jend = nAsyRay; jinc = 1
END IF
! ----------------------------------------------------  
DO imray = jbeg, jend, jinc
  iAsyRay = AsyRayList(imray)
  iazi    = AziList(imray)
  idir    = DirList(imray)
  IF (krot .EQ. 2) idir = mp(idir)
  
  DO ipol = 1, nPolarAng
    wt(ipol) = wtang(ipol, iazi)
    
    IF (lJout) wt2(ipol, 1:4) = wtang2(ipol, iazi, 1:4)
  END DO
  
  iAsyTyp = hAsy(iAsy)%AsyTyp
  iGeoTyp = hAsy(iAsy)%GeoTyp
  icBss   = hAsyTypInfo(iAsyTyp)%iBss
  
  haRay_Loc => haRay(iGeoTyp, icBss, iAsyRay)
  
  nPinRay = haRay_Loc%nhpRay
  
  IF (idir .EQ. 1) THEN
    ipst = 1; iped = nPinRay; ipinc = 1;  isfst = 1; isfed = 2;
  ELSE
    ipst = nPinRay; iped = 1; ipinc = -1; isfst = 2; isfed = 1;
  END IF
  
  DO ipray = ipst, iped, ipinc
    jhPin = haRay_Loc%CelRay(ipRay)%hPinIdx
    jhPin = hAsy(iAsy)%PinIdxSt + hAsyTypInfo(iAsyTyp)%PinLocIdx(iGeoTyp, jhPin) - 1
    jcBss = Pin(jhPin)%hCelGeo(iz)
    
    CelRay_Loc => haRay(iGeoTyp, jcBss, iAsyRay)%CelRay(ipRay)
    
    IF (lJout) THEN
      iSurf = CelRay_Loc%SurfIdx(isfst)
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          Jout(1, ig, isurf, ipin) = Jout(1, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
          Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
        END DO
      END DO
    END IF
    
    nRaySeg = CelRay_Loc%nSegRay
    
    IF (idir .EQ. 1) THEN
      isgst = 1; isged = nRaySeg; isginc = 1
    ELSE
      isgst = nRaySeg; isged = 1; isginc = -1
    END IF
    
    DO iRaySeg = isgst, isged, isginc
      ireg = CelRay_Loc%MshIdx(iRaySeg) + Pin(jhPin)%FsrIdxSt - 1
      
      DO ig = gb, ge
        tau = -CelRay_Loc%SegLgh(iRaySeg) * xst(ig, ireg) ! Optimum Length
        
        ExpAppIdx = max(INT(tau), -40000)
        ExpAppIdx = min(0, ExpAppIdx)
        
        DO ipol = 1, nPolarAng
          ExpApp = EXPA(ipol, ExpAppIdx) * tau + EXPB(ipol, ExpAppIdx)
          
          phid = (PhiAngOut(ipol, ig) - src(ig, ireg)) * ExpApp
          
          PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
          
          phis(ig, ireg) = phis(ig, ireg) + wt(ipol) * phid
        END DO
      END DO
    END DO
    
    IF (lJout) THEN
      isurf = CelRay_Loc%SurfIdx(isfed) ! y : Big
      
      DO ig = gb, ge
        DO ipol = 1, nPolarAng
          Jout(2, ig, isurf, ipin) = Jout(2, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
          Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
        END DO
      END DO
    END IF
  END DO
END DO

DcmpPhiAngOut(1:nPolarAng, gb:ge, krot, iRay, iAsy) = PhiAngOut
! ----------------------------------------------------

END SUBROUTINE HexTrackRotRayNM_Dcmp
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmp_LSCASMO(RayInfo, CoreInfo, phisnm, phisSlope, PhiAngInnm, srcnm, srcSlope,          &
                                xstnm, joutnm, iz, iasy, gb, ge, ljout)
USE PARAM
USE TYPEDEF,    ONLY :  RayInfo_Type,       Coreinfo_type,      Pin_Type,           Asy_Type,               &
                        AsyInfo_Type,       PinInfo_Type,       Cell_Type,          DcmpAsyRayInfo_Type
USE Moc_Mod,    ONLY :  nMaxRaySeg,         nMaxCellRay,        nMaxAsyRay,         nMaxCoreRay,            &
                        Expa,               Expb,               ApproxExp,          TrackingDat,            &
                        DcmpPhiAngIn,       DcmpPhiAngOut,      TrackRotRayLSDcmp_CASMO
USE geom,       ONLY :  ng
USE PE_Mod,     ONLY :  PE
USE ALLOCS
USE OMP_LIB
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phisnm(:, :), PhiAngInnm(:, :, :)
REAL, POINTER :: srcnm(:, :), xstnm(:, :), joutnm(:, :, :, :)
REAL, POINTER :: phisSlope(:, :, :, :), srcSlope(:, :, :, :)
INTEGER :: iz, iasy, gb, ge
LOGICAL :: ljout

TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(DcmpAsyRayInfo_Type), POINTER :: DcmpAsyRay(:, :)
INTEGER, POINTER :: DcmpAsyRayCount(:)
REAL, POINTER :: FsrMxx(:), FsrMyy(:), FsrMxy(:)

INTEGER :: AsyType

INTEGER :: nxy, nFsr, nThread
INTEGER :: iAsyRay
INTEGER :: tid
INTEGER :: FsrIdxSt, FsrIdxEnd
INTEGER :: PinIdxSt, PinIdxEnd
INTEGER :: icel, ireg
INTEGER :: i, j, l, g
REAL :: detinv

REAL, POINTER :: phimx(:, :, :), phimy(:, :, :)

AsyInfo => CoreInfo%AsyInfo; Asy => CoreInfo%Asy
Cell => CoreInfo%CellInfo; Pin => CoreInfo%Pin
FsrMxx => CoreInfo%CoreFsrMxx(:, iz)
FsrMyy => CoreInfo%CoreFsrMyy(:, iz)
FsrMxy => CoreInfo%CoreFsrMxy(:, iz)

AsyType = Asy(iAsy)%AsyType
nxy = AsyInfo(AsyType)%nxy
nThread = PE%nThread

PinIdxSt = Asy(iAsy)%GlobalPinIdx(1)
PinIdxEnd = Asy(iAsy)%GlobalPinIdx(nxy)
FsrIdxSt = Pin(PinIdxSt)%FsrIdxSt
FsrIdxEnd = Pin(PinIdxEnd)%FsrIdxSt + Cell(Pin(PinIdxEnd)%Cell(iz))%nFsr - 1

DcmpAsyRay => RayInfo%DcmpAsyRay
DcmpAsyRayCount => RayInfo%DcmpAsyRayCount

ALLOCATE(phimx(2, gb : ge, FsrIdxSt : FsrIdxEnd), phimy(2, gb : ge, FsrIdxSt : FsrIdxEnd))
phimx(:, :, :) = zero; phimy(:, :, :) = zero

phisnm(gb : ge, FsrIdxSt : FsrIdxEnd) = zero
IF(ljout) joutnm(:, gb : ge, :, PinIdxSt : PinIdxEnd) = zero

tid = omp_get_thread_num() + 1
DO iAsyRay = 1, DcmpAsyRayCount(iAsy)
  CALL TrackRotRayLSDcmp_CASMO(RayInfo, CoreInfo, TrackingDat(tid), phisnm, phimx, phimy, joutnm,   &
                               DcmpAsyRay(iAsyRay, iAsy), ljout, iz, gb, ge)
END DO
    
DO i = FsrIdxSt, FsrIdxEnd
  DO g = gb, ge
    phimx(1, g, i) = phimx(1, g, i) / xstnm(g, i)
    phimy(1, g, i) = phimy(1, g, i) / xstnm(g, i)
  END DO
END DO

DO l = PinIdxSt, PinIdxEnd
  icel = Pin(l)%Cell(iz);
  DO j = 1, Cell(icel)%nFsr
    ireg = Pin(l)%FsrIdxSt + j - 1
    DO g = gb, ge
      phisnm(g, ireg) = phisnm(g, ireg) / (xstnm(g, ireg) * Cell(icel)%vol(j)) + srcnm(g, ireg)
      phimx(:, g, ireg) = phimx(:, g, ireg) / (xstnm(g, ireg) * Cell(icel)%vol(j))
      phimy(:, g, ireg) = phimy(:, g, ireg) / (xstnm(g, ireg) * Cell(icel)%vol(j))
    END DO
  END DO  
END DO

DO i = FsrIdxSt, FsrIdxEnd
  detinv = 1.0 / (FsrMxx(i) * FsrMyy(i) - FsrMxy(i) * FsrMxy(i))
  DO g = gb, ge
    phisSlope(1, g, i, iz) = detinv * (FsrMyy(i) * SUM(phimx(:, g, i)) - FsrMxy(i) * SUM(phimy(:, g, i)))
    phisSlope(2, g, i, iz) = detinv * (FsrMxx(i) * SUM(phimy(:, g, i)) - FsrMxy(i) * SUM(phimx(:, g, i)))
  END DO
END DO

DEALLOCATE(phimx, phimy)

NULLIFY(Cell); NULLIFY(Pin)
NULLIFY(FsrMxx); NULLIFY(FsrMyy); NULLIFY(FsrMxy)

END SUBROUTINE RayTraceDcmp_LSCASMO
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE TrackRotRayLSDcmp_CASMO(RayInfo, CoreInfo, TrackingDat, phisC, phimx, phimy, jout, DcmpAsyRay, ljout, iz, gb, ge)
USE PARAM
USE TYPEDEF,    ONLY :  RayInfo_Type,       Coreinfo_type,      Pin_Type,           Asy_Type,               &
                        AsyInfo_Type,       PinInfo_Type,       Cell_Type,          AziAngleInfo_Type,      &
                        PolarAngle_Type,    AsyRayInfo_type,    CoreRayInfo_Type,   DcmpAsyRayInfo_Type,    &
                        RotRayInfo_Type,    CellRayInfo_type,   TrackingDat_Type
USE Moc_Mod,    ONLY :  nMaxCellRay,        nMaxCoreRay
USE BasicOperation, ONLY : CP_CA, CP_VA
IMPLICIT NONE

TYPE(RayInfo_Type), INTENT(INOUT) :: RayInfo
TYPE(CoreInfo_Type), INTENT(INOUT) :: CoreInfo
TYPE(TrackingDat_Type), INTENT(INOUT) :: TrackingDat
TYPE(DcmpAsyRayInfo_Type) :: DcmpAsyRay
REAL, POINTER :: phisC(:, :), phimx(:, :, :), phimy(:, :, :), jout(:, :, :, :)
LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: iz, gb, ge

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(PinInfo_Type), POINTER :: PinInfo(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(AziAngleInfo_Type), POINTER :: AziAng(:)
TYPE(PolarAngle_Type), POINTER :: PolarAng(:)
TYPE(AsyRayInfo_type), POINTER :: AsyRay(:)
TYPE(CoreRayInfo_Type), POINTER :: CoreRay(:)
TYPE(RotRayInfo_Type), POINTER :: RotRay(:)
TYPE(CellRayInfo_Type), POINTER :: CellRay

REAL, POINTER :: LenSeg(:)
INTEGER, POINTER :: LocalFsrIdx(:)
INTEGER, POINTER :: FsrIdx(:, :),  ExpAppIdx(:, :, :)
REAL, POINTER :: OptLenList(:, :, :)
REAL, POINTER :: srcC(:, :), srcSlope(:, :, :), xst(:, :)
REAL, POINTER :: PhiAngOut(:, :, :), PhiAngIn(:, :, :)
REAL, POINTER :: DcmpPhiAngOut(:, :, :, :, :), DcmpPhiAngIn(:, :, :, :, :)
REAL, POINTER :: EXPA(:, :), EXPB(:, :)
REAL, POINTER :: E1(:, :, :, :), E3(:, :, :, :), R1(:, :, :, :), R3(:, :, :, :)
REAL, POINTER :: cmOptLen(:, :, :, :), cmOptLenInv(:, :, :, :)
REAL, POINTER :: q0(:, :, :), q1(:, :, :, :)
REAL, POINTER :: x0(:, :, :), y0(:, :, :)
REAL, POINTER :: FsrCentroid(:, :)
INTEGER, POINTER :: AsyRayList(:), DirList(:), AziList(:)

INTEGER :: mp(2)
INTEGER :: iazi, ipol, ig, iRotRay, iCoreRay, iasyray, iceray, irayseg, irot, itype, idir, iray
INTEGER :: nRotRay, nCoreRay, nAsyRay, nPinRay, nRaySeg, FsrIdxSt
INTEGER :: nPolarAng, nAziAng, nPhiAngSv
INTEGER :: ipin, icel, iasy, ireg, isurf
INTEGER :: irsegidx, icellrayidx
INTEGER :: nFsr, nxy
INTEGER :: i, j, k, l, m, jbeg, jend, jinc, ir, ir1
INTEGER :: ibcel

INTEGER :: CellRayIdxSt(DcmpAsyRay%nMaxCellRay, DcmpAsyRay%nAsyRay, 2)
INTEGER :: PinIdx(DcmpAsyRay%nMaxCellRay, DcmpAsyRay%nAsyRay)
INTEGER :: SurfIdx(DcmpAsyRay%nMaxCellRay, DcmpAsyRay%nAsyRay, 2)
INTEGER :: PhiAnginSvIdx(2)
INTEGER :: nTotRaySeg(DcmpAsyRay%nAsyRay), nTotCellRay(DcmpAsyRay%nAsyRay)

REAL :: wtang(10, 100), wttemp, wtang2(10, 100, 4)
REAL :: tau, phiobd(10, gb : ge), phid, phim, wt(10), wt2(10, 4)
REAL :: ax(10), ay(10), polsininv(10)
REAL :: segCenter(2), segSlope(2, 10)

DATA mp /2, 1/

!Ray Info Pointing
AziAng => RayInfo%AziAngle;
PolarAng => RayInfo%PolarAngle;
AsyRay => RayInfo%AsyRay
CoreRay => RayInfo%CoreRay
RotRay => RayInfo%RotRay

!Geometry Info Pointing
Asy => CoreInfo%Asy
Pin => CoreInfo%Pin
PinInfo => CoreInfo%Pininfo
Cell => CoreInfo%CellInfo
FsrCentroid => CoreInfo%CoreFsrCentroid(:, :, iz)

!Tracking Dat Pointing
FsrIdx => TrackingDat%FsrIdx
ExpAppIdx => TrackingDat%ExpAppIdxnm
OptLenList => TrackingDat%OptLenListnm
cmOptLen => TrackingDat%cmOptLen; cmOptLenInv => TrackingDat%cmOptLenInv
srcC => TrackingDat%srcnm; srcSlope => TrackingDat%srcSlope
xst => TrackingDat%xstnm
PhiAngOut => TrackingDat%PhiAngOutnm
PhiAngIn => TrackingDat%phiAngInnm
DcmpPhiAngIn => TrackingDat%DcmpPhiAngIn
DcmpPhiAngOut => TrackingDat%DcmpPhiAngOut
EXPA => TrackingDat%EXPA
EXPB => TrackingDat%EXPB
E1 => TrackingDat%E1; E3 => TrackingDat%E3
R1 => TrackingDat%R1; R3 => TrackingDat%R3
q0 => TrackingDat%q0; q1 => TrackingDat%q1
x0 => TrackingDat%x0; y0 => TrackingDat%y0

nAziAng = RayInfo%nAziAngle; nPolarAng = RayInfo%nPolarAngle
nPhiAngSv = RayInfo%nPhiAngSv
nFsr = CoreInfo%nCoreFsr; nxy = CoreInfo%nxy

DO ipol = 1, nPolarAng
  wttemp = PolarAng(ipol)%weight * PolarAng(ipol)%sinv
  polsininv(ipol) = 1._8 / PolarAng(ipol)%sinv
  DO iazi = 1, nAziAng
    wtang(ipol, iazi) = wttemp * AziAng(iazi)%weight * AziAng(iazi)%del
  END DO
END DO
IF(lJout)THEN
  DO ipol = 1, nPolarAng
    DO iazi = 1, nAziAng
      wtang2(ipol, iazi, 1) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv)
      wtang2(ipol, iazi, 3) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv)
      wtang2(ipol, iazi, 2) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv)
      wtang2(ipol, iazi, 4) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv)
    END DO
  END DO
END IF

nAsyRay = DcmpAsyRay%nAsyRay
iRotRay = DcmpAsyRay%iRotRay
AsyRayList => DcmpAsyRay%AsyRayList
DirList => DcmpAsyRay%DirList
AziList => DcmpAsyRay%AziList
iAsy = DcmpAsyRay%iAsy; iRay = DcmpAsyRay%iRay
PhiAngInSvIdx = RayInfo%PhiAngInSvIdx(iRotRay, 1 : 2)

DO j = 1, nAsyRay
  irsegidx = 0
  iAsyRay = AsyRayList(j); iAzi = AziList(j)
  nPinRay = AsyRay(iAsyRay)%nCellRay
  DO ipol = 1, nPolarAng
    segSlope(:, ipol) = (/ AziAng(iazi)%cosv * PolarAng(ipol)%sinv, AziAng(iazi)%sinv * PolarAng(ipol)%sinv /)
  END DO
  DO l = 1, nPinRay
    ipin = AsyRay(iAsyRay)%PinIdx(l)
    iceray = AsyRay(iAsyRay)%PinRayIdx(l)
    ipin = Asy(iAsy)%GlobalPinIdx(ipin)
    icel = Pin(ipin)%Cell(iz)
    FsrIdxSt = Pin(ipin)%FsrIdxSt
    ibcel = Cell(icel)%basecellstr
    CellRay => Cell(ibcel)%CellRay(iceray)
    PinIdx(l, j) = ipin
    CellRayIdxSt(l, j, 2) = irsegidx + 1
    nRaySeg = CellRay%nSeg
    LocalFsrIdx => CellRay%LocalFsrIdx
    LenSeg => CellRay%LenSeg
    DO iRaySeg = 1, nRaySeg
      ireg = FsrIdxSt + CellRay%LocalFsrIdx(iRaySeg) - 1
      irsegidx = irsegidx + 1
      FsrIdx(irsegidx, j) = ireg
      segCenter(:) = half * (CellRay%pts(:, iRaySeg) + CellRay%pts(:, iRaySeg + 1)) - FsrCentroid(:, ireg)
      x0(1, irsegidx, j) = segCenter(1) - LenSeg(iRaySeg) * 0.0005_8 * AziAng(iazi)%cosv
      x0(2, irsegidx, j) = segCenter(1) + LenSeg(iRaySeg) * 0.0005_8 * AziAng(iazi)%cosv
      y0(1, irsegidx, j) = segCenter(2) - LenSeg(iRaySeg) * 0.0005_8 * AziAng(iazi)%sinv
      y0(2, irsegidx, j) = segCenter(2) + LenSeg(iRaySeg) * 0.0005_8 * AziAng(iazi)%sinv
      DO ig = gb, ge
        tau = - LenSeg(iRaySeg) * xst(ig, ireg)
        OptLenList(ig, irsegidx, j) = tau
        q0(ig, irsegidx, j) = srcC(ig, ireg) + sum(srcSlope(1 : 2, ig, ireg) * segCenter(:))
        DO ipol = 1, nPolarAng
          cmOptLen(ipol, ig, irsegidx, j) = - tau * polsininv(ipol) * 0.001_8
          cmOptLenInv(ipol, ig, irsegidx, j) = 1._8 / cmOptLen(ipol, ig, irsegidx, j)
          q1(ipol, ig, irsegidx, j) = half * sum(srcSlope(3 : 4, ig, ireg) * segSlope(:, ipol))
        END DO
        ExpAppIdx(ig, irsegidx, j) = max(INT(tau), -40000)
        ExpAppIdx(ig, irsegidx, j) = min(0, ExpAppIdx(ig, irsegidx, j))
      END DO
    END DO
    CellRayIdxSt(l, j, 1) = irsegidx
    SurfIdx(l, j, 1) = AsyRay(iAsyRay)%PinRaySurf(2, l) !OutSurface
    SurfIdx(l, j, 2) = AsyRay(iAsyRay)%PinRaySurf(1, l) !Insurface
  END DO
  nTotRaySeg(j) = irsegidx
  nTotCellRay(j) = nPinRay
END DO

DO ipol = 1, nPolarAng
  DO j = 1, nAsyRay
    DO l = 1, nTotRaySeg(j)
      DO ig = gb, ge
        E1(ipol, ig, l, j) = EXPA(ExpAppIdx(ig, l, j), ipol) * OptLenList(ig, l, j) + EXPB(ExpAppIdx(ig, l, j), ipol)
        E3(ipol, ig, l, j) = 2._8 * (cmOptLen(ipol, ig, l, j) - E1(ipol, ig, l, j)) - cmOptLen(ipol, ig, l, j) * E1(ipol, ig, l, j)
        R1(ipol, ig, l, j) = 1._8 + half * cmOptLen(ipol, ig, l, j) - (1._8 + cmOptLenInv(ipol, ig, l, j)) * E1(ipol, ig, l, j)
        R3(ipol, ig, l, j) = rsix * cmOptLen(ipol, ig, l, j) - 2._8 * cmOptLenInv(ipol, ig, l, j) - 2._8  &
                             + (1._8 + cmOptLenInv(ipol, ig, l, j)) * (1._8 + 2._8 * cmOptLenInv(ipol, ig, l, j)) * E1(ipol, ig, l, j)
      END DO
    END DO
  END DO
END DO

DO irot = 1, 2
  IF (DcmpAsyRay%lRotRayBeg(irot)) THEN
    phiobd(1 : nPolarAng, gb : ge) = PhiAngIn(1 : nPolarAng, gb : ge, PhiAnginSvIdx(irot))
  ELSE
    phiobd(1 : nPolarAng, gb : ge) = DcmpPhiAngIn(1 : nPolarAng, gb : ge, irot, iRay, iAsy)
  END IF
  jinc = 1; jbeg = 1; jend = nAsyRay
  IF(irot .eq. 2) THEN
    jinc = -1; jbeg = nAsyRay; jend = 1
  END IF
  DO j = jbeg, jend, jinc
    idir = DirList(j); iazi = AziList(j)
    nRaySeg = nTotRaySeg(j)
    IF(irot .eq. 2) idir = mp(idir)
    DO ipol = 1, nPolarAng
      wt(ipol) = wtang(ipol, iazi)
      ax(ipol) = wt(ipol) * AziAng(iazi)%cosv * PolarAng(ipol)%sinv
      ay(ipol) = wt(ipol) * AziAng(iazi)%sinv * PolarAng(ipol)%sinv
    END DO
    IF(lJout) THEN
      wt2(1 : nPolarAng, :) = wtang2(1 : nPolarAng, iazi, :)
    END IF
    nRaySeg = nTotRaySeg(j)
    IF(idir .eq. 1) THEN
      PhiAngOut(1 : nPolarAng, gb : ge, 1) = phiobd(1 : nPolarAng, gb : ge)
      DO ir = 1, nRaySeg
        ireg = FsrIdx(ir, j)
        DO ig = gb, ge
          DO ipol = 1, nPolarAng          
            phid = (q0(ig, ir, j) - PhiAngOut(ipol, ig, ir)) * E1(ipol, ig, ir, j) + q1(ipol, ig, ir, j) * E3(ipol, ig, ir, j)
            PhiAngOut(ipol, ig, ir + 1) = PhiAngOut(ipol, ig, ir) + phid
            phisC(ig, ireg) = phisC(ig, ireg) - wt(ipol) * phid
            phim = PhiAngOut(ipol, ig, ir) * (cmOptLen(ipol, ig, ir, j) * half) + (q0(ig, ir, j) - PhiAngOut(ipol, ig, ir)) * R1(ipol, ig, ir, j)  &
                   + q1(ipol, ig, ir, j) * cmOptLen(ipol, ig, ir, j) * R3(ipol, ig, ir, j)
            phimx(1, ig, ireg) = phimx(1, ig, ireg) + ax(ipol) * phim * cmOptLen(ipol, ig, ir, j)
            phimy(1, ig, ireg) = phimy(1, ig, ireg) + ay(ipol) * phim * cmOptLen(ipol, ig, ir, j)
            phimx(2, ig, ireg) = phimx(2, ig, ireg) + wt(ipol) * x0(idir, ir, j) * (-phid + q0(ig, ir, j) * cmOptLen(ipol, ig, ir, j))
            phimy(2, ig, ireg) = phimy(2, ig, ireg) + wt(ipol) * y0(idir, ir, j) * (-phid + q0(ig, ir, j) * cmOptLen(ipol, ig, ir, j))
          END DO
        END DO
      END DO
      phiobd(1 : nPolarAng, gb : ge) = PhiAngOut(1 : nPolarAng, gb : ge, nRaySeg + 1)
      IF(ljout) THEN
        DO ir = 1, nTotCellRay(j)
          icel = PinIdx(ir, j); isurf = SurfIdx(ir, j, 1)
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(2, ig, isurf, icel) = Jout(2, ig, isurf, icel) + wt(ipol) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 1) + 1)
              Jout(3, ig, isurf, icel) = Jout(3, ig, isurf, icel) + wt2(ipol, isurf) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 1) + 1)
            END DO
          END DO
          isurf = SurfIdx(ir, j, 2)
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(1, ig, isurf, icel) = Jout(1, ig, isurf, icel) + wt(ipol) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 2))
              Jout(3, ig, isurf, icel) = Jout(3, ig, isurf, icel) + wt2(ipol, isurf) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 2))
            END DO
          END DO
        END DO
      END IF
    ELSE
      PhiAngOut(1 : nPolarAng, gb : ge, nRaySeg + 2) = phiobd(1 : nPolarAng, gb : ge)
      ir = nRaySeg + 1
      DO ir1 = 1, nRaySeg
        ir = ir - 1
        ireg = FsrIdx(ir, j)
        DO ig = gb, ge
          DO ipol = 1, nPolarAng
            phid = (q0(ig, ir, j) - PhiAngOut(ipol, ig, ir + 2)) * E1(ipol, ig, ir, j) - q1(ipol, ig, ir, j) * E3(ipol, ig, ir, j)
            PhiAngOut(ipol, ig, ir + 1) = PhiAngOut(ipol, ig, ir + 2) + phid
            phisC(ig, ireg) = phisC(ig, ireg) - wt(ipol) * phid
            phim = PhiAngOut(ipol, ig, ir + 2) * (cmOptLen(ipol, ig, ir, j) * half) + (q0(ig, ir, j) - PhiAngOut(ipol, ig, ir + 2)) * R1(ipol, ig, ir, j)  &
                   - q1(ipol, ig, ir, j) * cmOptLen(ipol, ig, ir, j) * R3(ipol, ig, ir, j)
            phimx(1, ig, ireg) = phimx(1, ig, ireg) - ax(ipol) * phim * cmOptLen(ipol, ig, ir, j)
            phimy(1, ig, ireg) = phimy(1, ig, ireg) - ay(ipol) * phim * cmOptLen(ipol, ig, ir, j)
            phimx(2, ig, ireg) = phimx(2, ig, ireg) + wt(ipol) * x0(idir, ir, j) * (-phid + q0(ig, ir, j) * cmOptLen(ipol, ig, ir, j))
            phimy(2, ig, ireg) = phimy(2, ig, ireg) + wt(ipol) * y0(idir, ir, j) * (-phid + q0(ig, ir, j) * cmOptLen(ipol, ig, ir, j))
          END DO
        END DO
      END DO
      phiobd(1 : nPolarAng, gb : ge) = PhiAngOut(1 : nPolarAng, gb : ge, 2)
      IF(lJout) THEN
        DO ir = 1, nTotCellRay(j)
          icel = PinIdx(ir, j); isurf = SurfIdx(ir, j, 2)
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(2, ig, isurf, icel) = Jout(2, ig, isurf, icel) + wt(ipol) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 2) + 1)
              Jout(3, ig, isurf, icel) = Jout(3, ig, isurf, icel) + wt2(ipol, isurf) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 2) + 1)
            END DO
          END DO
          isurf = SurfIdx(ir, j, 1)
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(1, ig, isurf, icel) = Jout(1, ig, isurf, icel) + wt(ipol) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 1) + 2)
              Jout(3, ig, isurf, icel) = Jout(3, ig, isurf, icel) + wt2(ipol, isurf) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 1) + 2)
            END DO
          END DO
        END DO
      END IF
    END IF
  END DO
  DcmpPhiAngOut(1 : nPolarAng, gb : ge, irot, iRay, iAsy) = phiobd(1 : nPolarAng, gb : ge)
END DO

NULLIFY(AziAng); NULLIFY(PolarAng)
NULLIFY(AsyRay); NULLIFY(CoreRay)
NULLIFY(RotRay); NULLIFY(CellRay)

!Geometry Info Pointing
NULLIFY(Asy); NULLIFY(Pin)
NULLIFY(PinInfo); NULLIFY(Cell)
NULLIFY(FsrCentroid)

!Tracking Dat Pointing
NULLIFY(FsrIdx); NULLIFY(ExpAppIdx)
NULLIFY(OptLenList); NULLIFY(cmOptLen); NULLIFY(cmOptLenInv)
NULLIFY(LenSeg); NULLIFY(LocalFsrIdx)
NULLIFY(srcC); NULLIFY(srcSlope)
NULLIFY(xst)
NULLIFY(PhiAngOut); NULLIFY(PhiAngIn)
NULLIFY(EXPA); NULLIFY(EXPB)
NULLIFY(E1); NULLIFY(E3); NULLIFY(R1); NULLIFY(R3)
NULLIFY(q0); NULLIFY(q1); NULLIFY(x0); NULLIFY(y0)

END SUBROUTINE TrackRotRayLSDcmp_CASMO
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmp_Pn(RayInfo, CoreInfo, phisnm, phimnm, PhiAngInnm, xstnm, srcnm, srcmnm, joutnm, iz, iAsy, gb, ge, ScatOd, lJout)

USE TIMER
USE ALLOCS
USE OMP_LIB
USE PARAM,   ONLY : TRUE, ZERO, RTHREE, RFIVE, RSEVEN
USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, TrackingDat_Type, Pin_Type, Cell_Type, PolarAngle_Type, AziAngleInfo_Type, Asy_Type, AsyInfo_Type, DcmpAsyRayInfo_Type
USE GEOM,    ONLY : ng
USE Moc_Mod, ONLY : nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay, Expa, Expb, TrackingDat, ApproxExp, wtang, Comp, mwt, AziMap, DcmpPhiAngIn, DcmpPhiAngOut, TrackRotRayPn_Dcmp
USE PE_MOD,  ONLY : PE

USE BasicOperation, ONLY : CP_CA, CP_VA

IMPLICIT NONE

TYPE (RayInfo_Type)  :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:,:)     :: phisnm, xstnm, srcnm
REAL, POINTER, DIMENSION(:,:,:)   :: phimnm, PhiAngInnm, srcmnm
REAL, POINTER, DIMENSION(:,:,:,:) :: joutnm

INTEGER :: iz, iAsy, gb, ge, ScatOd
LOGICAL :: ljout

LOGICAL, SAVE :: lfirst
DATA lfirst /TRUE/

TYPE (AziAngleInfo_Type),   POINTER, DIMENSION(:)   :: AziAng
TYPE (PolarAngle_Type),     POINTER, DIMENSION(:)   :: PolarAng
TYPE (Asy_Type),            POINTER, DIMENSION(:)   :: Asy
TYPE (AsyInfo_Type),        POINTER, DIMENSION(:)   :: AsyInfo
TYPE (Cell_Type),           POINTER, DIMENSION(:)   :: Cell
TYPE (Pin_Type),            POINTER, DIMENSION(:)   :: Pin
TYPE (DcmpAsyRayInfo_Type), POINTER, DIMENSION(:,:) :: DcmpAsyRay

REAL, POINTER, DIMENSION(:,:,:,:) :: phianm, SrcAngnm

INTEGER, POINTER, DIMENSION(:,:,:) :: DcmpAsyAziList

INTEGER :: nAziAng, nPolarAng, nAsyRay, nAsy, nxy, tid, FsrIdxSt, FsrIdxEnd, PinIdxSt, PinIdxEnd, AziIdx
INTEGER :: icel, ireg, ipin, iazi, ipol, iDcmpAsyRay, iCoreRay, i, j, l, k, m, g
REAL :: wttemp, wtcos, wtpolar, wtsin2, tempsrc
! ----------------------------------------------------

Asy     => CoreInfo%Asy
AsyInfo => CoreInfo%AsyInfo
Cell    => CoreInfo%CellInfo
Pin     => CoreInfo%Pin
nAsy     = CoreInfo%nxya

AziAng         => RayInfo%AziAngle
PolarAng       => RayInfo%PolarAngle
DcmpAsyRay     => RayInfo%DcmpAsyRay
DcmpAsyAziList => RayInfo%DcmpAsyAziList
nAziAng         = RayInfo%nAziAngle
nPolarAng       = RayInfo%nPolarAngle

tid = omp_get_thread_num() + 1

nxy = AsyInfo(Asy(iAsy)%AsyType)%nxy

PinIdxSt  = Asy(iAsy)%GlobalPinIdx(1)
PinIdxEnd = Asy(iAsy)%GlobalPinIdx(nxy)

FsrIdxSt  = Pin(Asy(iAsy)%GlobalPinIdx(1))%FsrIdxSt
FsrIdxEnd = Pin(Asy(iAsy)%GlobalPinIdx(nxy))%FsrIdxSt + Cell(Pin(Asy(iAsy)%GlobalPinIdx(nxy))%Cell(iz))%nFsr - 1
   
CALL dmalloc0(TrackingDat(tid)%phianm,   1, nPolarAng, gb, ge, FsrIdxSt, FsrIdxEnd, 1, 4)
CALL dmalloc0(TrackingDat(tid)%SrcAngnm, 1, nPolarAng, gb, ge, FsrIdxSt, FsrIdxEnd, 1, 4)

phianm   => TrackingDat(tid)%phianm
SrcAngnm => TrackingDat(tid)%SrcAngnm

phisnm(   gb:ge, FsrIdxSt:FsrIdxEnd) = ZERO
phimnm(:, gb:ge, FsrIdxSt:FsrIdxEnd) = ZERO

IF (lJout) joutnm(:, gb:ge, :, PinIdxSt:PinIdxEnd) = ZERO
! ----------------------------------------------------
DO iAzi = 1, nAziAng / 2
  IF (ScatOd .EQ. 1) THEN
    DO ireg = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      
      DO g = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = srcnm(g, ireg)
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = srcnm(g, ireg)
          
          tempsrc = Comp(1, ipol, AziIdx) * srcmnm(1, g, ireg) + Comp(2, ipol, AziIdx) * srcmnm(2, g, ireg)
          
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) - tempsrc
        END DO
      END DO
      
      AziIdx = nAziAng - iAzi + 1
      
      DO g = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = srcnm(g, ireg)
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = srcnm(g, ireg)
          
          tempsrc = Comp(1, ipol, AziIdx) * srcmnm(1, g, ireg) + Comp(2, ipol, AziIdx) * srcmnm(2, g, ireg)
          
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) - tempsrc
        END DO
      END DO
    END DO
  ELSE IF (ScatOd .EQ. 2) THEN 
    DO ireg = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      
      DO g = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = srcnm(g, ireg)
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = srcnm(g, ireg)
          
          tempsrc = Comp(1, ipol, AziIdx) * srcmnm(1, g, ireg) + Comp(2, ipol, AziIdx) * srcmnm(2, g, ireg)
          
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) - tempsrc
          
          tempsrc = Comp(3, ipol, AziIdx) * srcmnm(3, g, ireg) + Comp(4, ipol, AziIdx) * srcmnm(4, g, ireg) + Comp(5, ipol, AziIdx) * srcmnm(5, g, ireg)
          
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) + tempsrc
        END DO
      END DO
      
      AziIdx = nAziAng - iAzi + 1
      
      DO g = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = srcnm(g, ireg)
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = srcnm(g, ireg)
          
          tempsrc = Comp(1, ipol, AziIdx) * srcmnm(1, g, ireg) + Comp(2, ipol, AziIdx) * srcmnm(2, g, ireg)
          
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) - tempsrc
          
          tempsrc = Comp(3, ipol, AziIdx) * srcmnm(3, g, ireg) + Comp(4, ipol, AziIdx) * srcmnm(4, g, ireg) + Comp(5, ipol, AziIdx) * srcmnm(5, g, ireg)
          
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) + tempsrc
        END DO
      END DO
    END DO
  ELSEIF (ScatOd .EQ. 3) THEN
    DO ireg = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      
      DO g = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = srcnm(g, ireg)
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = srcnm(g, ireg)
          
          tempsrc = Comp(1, ipol, AziIdx) * srcmnm(1, g, ireg) + Comp(2, ipol, AziIdx) * srcmnm(2, g, ireg) &
                  + Comp(6, ipol, AziIdx) * srcmnm(6, g, ireg) + Comp(7, ipol, AziIdx) * srcmnm(7, g, ireg) &
                  + Comp(8, ipol, AziIdx) * srcmnm(8, g, ireg) + Comp(9, ipol, AziIdx) * srcmnm(9, g, ireg)
          
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) - tempsrc
          
          tempsrc = Comp(3, ipol, AziIdx) * srcmnm(3, g, ireg) + Comp(4, ipol, AziIdx) * srcmnm(4, g, ireg) &
                  + Comp(5, ipol, AziIdx) * srcmnm(5, g, ireg)
          
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) + tempsrc
        END DO
      END DO
      
      AziIdx = nAziAng - iAzi + 1
      
      DO g = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = srcnm(g, ireg)
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = srcnm(g, ireg)
          
          tempsrc = Comp(1, ipol, AziIdx) * srcmnm(1, g, ireg) + Comp(2, ipol, AziIdx) * srcmnm(2, g, ireg) &
                  + Comp(6, ipol, AziIdx) * srcmnm(6, g, ireg) + Comp(7, ipol, AziIdx) * srcmnm(7, g, ireg) &
                  + Comp(8, ipol, AziIdx) * srcmnm(8, g, ireg) + Comp(9, ipol, AziIdx) * srcmnm(9, g, ireg)
          
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) - tempsrc
          
          tempsrc = Comp(3, ipol, AziIdx) * srcmnm(3, g, ireg) + Comp(4, ipol, AziIdx) * srcmnm(4, g, ireg) &
                  + Comp(5, ipol, AziIdx) * srcmnm(5, g, ireg)
          
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) + tempsrc
        END DO
      END DO
    END DO
  END IF

  phianm = ZERO
  
  DO i = 1, DcmpAsyAziList(0, iAzi, iAsy)
    iDcmpAsyRay = DcmpAsyAziList(i, iAzi, iAsy)
    
    CALL TrackRotRayPn_Dcmp(RayInfo, CoreInfo, TrackingDat(tid), DcmpAsyRay(iDcmpAsyRay, iAsy), Joutnm, lJout, iz, gb, ge)
  END DO
  
  IF (ScatOd .EQ. 1) THEN
    DO ireg = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      
      DO g = gb, ge
        DO ipol = 1, nPolarAng            
          phisnm(g, ireg) = phisnm(g, ireg) + wtang(ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) + phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
          
          phimnm(1:2, g, ireg) = phimnm(1:2, g, ireg) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) - phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
        END DO
      END DO
      
      AziIdx = nAziAng - iAzi + 1
      
      DO g = gb, ge
        DO ipol = 1, nPolarAng            
          phisnm(g, ireg) = phisnm(g, ireg) + wtang(ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) + phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
          
          phimnm(1:2, g, ireg) = phimnm(1:2, g, ireg) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) - phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
        END DO
      END DO
    END DO
  ELSEIF (ScatOd .EQ. 2) THEN
    DO ireg = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      
      DO g = gb, ge
        DO ipol = 1, nPolarAng            
          phisnm(g, ireg) = phisnm(g, ireg) + wtang(ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) + phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
          
          phimnm(1:2, g, ireg) = phimnm(1:2, g, ireg) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) - phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
          phimnm(3:5, g, ireg) = phimnm(3:5, g, ireg) + mwt(3:5, ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) + phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
        END DO
      END DO
      
      AziIdx = nAziAng - iAzi + 1
      
      DO g = gb, ge
        DO ipol = 1, nPolarAng            
          phisnm(g, ireg) = phisnm(g, ireg) + wtang(ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) + phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
          
          phimnm(1:2, g, ireg) = phimnm(1:2, g, ireg) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) - phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
          phimnm(3:5, g, ireg) = phimnm(3:5, g, ireg) + mwt(3:5, ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) + phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
        END DO
      END DO
    END DO
  ELSEIF (ScatOd .EQ. 3) THEN
    DO ireg = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      
      DO g = gb, ge
        DO ipol = 1, nPolarAng            
          phisnm(g, ireg) = phisnm(g, ireg) + wtang(ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) + phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
          
          phimnm(1:2, g, ireg) = phimnm(1:2, g, ireg) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) - phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
          phimnm(3:5, g, ireg) = phimnm(3:5, g, ireg) + mwt(3:5, ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) + phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
          phimnm(6:9, g, ireg) = phimnm(6:9, g, ireg) + mwt(6:9, ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) - phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
        END DO
      END DO
      
      AziIdx = nAziAng - iAzi + 1
      
      DO g = gb, ge
        DO ipol = 1, nPolarAng            
          phisnm(g, ireg) = phisnm(g, ireg) + wtang(ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) + phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
          
          phimnm(1:2, g, ireg) = phimnm(1:2, g, ireg) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) - phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
          phimnm(3:5, g, ireg) = phimnm(3:5, g, ireg) + mwt(3:5, ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) + phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
          phimnm(6:9, g, ireg) = phimnm(6:9, g, ireg) + mwt(6:9, ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) - phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
        END DO
      END DO
    END DO
  END IF
END DO
    
DEALLOCATE (TrackingDat(tid)%phianm)
DEALLOCATE (TrackingDat(tid)%SrcAngnm)
! ----------------------------------------------------
IF (ScatOd .EQ. 1) THEN
  !$OMP PARALLEL DO PRIVATE(FsrIdxSt, ireg, icel, wttemp) SCHEDULE(DYNAMIC)
  DO l = PinIdxSt, PinIdxEnd
    FsrIdxSt = Pin(l)%FsrIdxSt
    icel     = Pin(l)%Cell(iz)
    
    DO j = 1, Cell(icel)%nFsr
      ireg = FsrIdxSt + j - 1
      
      DO g = gb, ge
        wttemp = 1._8 / (xstnm(g, ireg) * Cell(icel)%vol(j))
        
        phisnm(g, ireg) = phisnm(g, ireg) * wttemp + srcnm(g, ireg)
        
        phimnm(1:2, g, ireg) = phimnm(1:2, g, ireg) * wttemp + srcmnm(1:2, g, ireg) * rthree
      END DO
    END DO  
  END DO
  !$OMP END PARALLEL DO
ELSEIF (ScatOd .EQ. 2) THEN
  !$OMP PARALLEL DO PRIVATE(FsrIdxSt, ireg, icel, wttemp) SCHEDULE(DYNAMIC)
  DO l = PinIdxSt, PinIdxEnd
    FsrIdxSt = Pin(l)%FsrIdxSt
    icel     = Pin(l)%Cell(iz)
    DO j = 1, Cell(icel)%nFsr
      ireg = FsrIdxSt + j - 1
      
      DO g = gb, ge
        wttemp = 1._8 / (xstnm(g, ireg) * Cell(icel)%vol(j))
        
        phisnm(g, ireg) = phisnm(g, ireg) * wttemp + srcnm(g, ireg)
        
        phimnm(1:2, g, ireg) = phimnm(1:2, g, ireg) * wttemp + srcmnm(1:2, g, ireg) * rthree
        phimnm(3:5, g, ireg) = phimnm(3:5, g, ireg) * wttemp + srcmnm(3:5, g, ireg) * rfive
      END DO
    END DO  
  END DO
  !$OMP END PARALLEL DO
ELSEIF (ScatOd .EQ. 3) THEN
  !$OMP PARALLEL DO PRIVATE(FsrIdxSt, ireg, icel, wttemp) SCHEDULE(DYNAMIC)
  DO l = PinIdxSt, PinIdxEnd
    FsrIdxSt = Pin(l)%FsrIdxSt
    icel     = Pin(l)%Cell(iz)
    
    DO j = 1, Cell(icel)%nFsr
      ireg = FsrIdxSt + j - 1
      
      DO g = gb, ge
        wttemp = 1._8 / (xstnm(g, ireg) * Cell(icel)%vol(j))
        
        phisnm(g, ireg) = phisnm(g, ireg) * wttemp + srcnm(g, ireg)
        
        phimnm(1:2, g, ireg) = phimnm(1:2, g, ireg) * wttemp + srcmnm(1:2, g, ireg) * rthree
        phimnm(3:5, g, ireg) = phimnm(3:5, g, ireg) * wttemp + srcmnm(3:5, g, ireg) * rfive
        phimnm(6:9, g, ireg) = phimnm(6:9, g, ireg) * wttemp + srcmnm(6:9, g, ireg) * rseven
      END DO
    END DO  
  END DO
  !$OMP END PARALLEL DO
END IF
! ----------------------------------------------------

END SUBROUTINE RayTraceDcmp_Pn
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE TrackRotRayPn_Dcmp(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, Jout, lJout, iz, gb, ge)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, Asy_Type, AsyInfo_Type, PinInfo_Type, Cell_Type, AziAngleInfo_Type, PolarAngle_Type, ModRayInfo_type, &
                    AsyRayInfo_type, CoreRayInfo_Type, RotRayInfo_Type, CellRayInfo_type, TrackingDat_Type, DcmpAsyRayInfo_Type
USE Moc_Mod, ONLY : nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay

USE BasicOperation, ONLY : CP_CA, CP_VA

IMPLICIT NONE

TYPE (RayInfo_Type)        :: RayInfo
TYPE (CoreInfo_Type)       :: CoreInfo
TYPE (TrackingDat_Type)    :: TrackingDat
TYPE (DcmpAsyRayInfo_Type) :: DcmpAsyRay

REAL, POINTER, DIMENSION(:,:,:,:) :: Jout

LOGICAL :: lJout
INTEGER :: iz, gb, ge, ScatOd

TYPE (Pin_Type),          POINTER, DIMENSION(:) :: Pin
TYPE (Asy_Type),          POINTER, DIMENSION(:) :: Asy
TYPE (PinInfo_Type),      POINTER, DIMENSION(:) :: PinInfo
TYPE (Cell_Type),         POINTER, DIMENSION(:) :: Cell
TYPE (AziAngleInfo_Type), POINTER, DIMENSION(:) :: AziAng
TYPE (PolarAngle_Type),   POINTER, DIMENSION(:) :: PolarAng
TYPE (AsyRayInfo_type),   POINTER, DIMENSION(:) :: AsyRay
TYPE (CoreRayInfo_Type),  POINTER, DIMENSION(:) :: CoreRay
TYPE (RotRayInfo_Type),   POINTER, DIMENSION(:) :: RotRay

TYPE (CellRayInfo_Type), POINTER :: CellRay

REAL, POINTER, DIMENSION(:)         :: LenSeg
REAL, POINTER, DIMENSION(:,:)       :: xst, EXPA, EXPB, wtang
REAL, POINTER, DIMENSION(:,:,:)     :: PhiAngIn
REAL, POINTER, DIMENSION(:,:,:,:)   :: SrcAng, phia
REAL, POINTER, DIMENSION(:,:,:,:,:) :: DcmpPhiAngIn, DcmpPhiAngOut

INTEGER, POINTER, DIMENSION(:)   :: LocalFsrIdx, AsyRayList, DirList, AziList
INTEGER, POINTER, DIMENSION(:,:) :: AziMap

INTEGER :: mp(2), AziSvIdx(2)
INTEGER :: iAzi, ipol, iRotRay, iCoreRay, iAsyRay, iRay, iRaySeg, irot, idir, jbeg, jend, jinc
INTEGER :: ipin, icel, iasy, ireg, isurf, ibcel, iceray, i, j, k, l, m, ig, ir
INTEGER :: nCoreRay, nAsyRay, nDcmpRay, nPinRay, nRaySeg, FsrIdxSt, nPolarAng, nAziAng, PhiAnginSvIdx, ExpAppIdx

REAL :: wtang2(10, 100, 4), wt(10), wt2(10, 4)
REAL :: PhiAngOut(RayInfo%nPolarAngle, gb : ge)
REAL :: phid, tau, ExpApp

DATA mp /2, 1/
! ----------------------------------------------------

! Ray Info Pointing
AziAng   => RayInfo%AziAngle
PolarAng => RayInfo%PolarAngle
AsyRay   => RayInfo%AsyRay
CoreRay  => RayInfo%CoreRay
RotRay   => RayInfo%RotRay
nAziAng   = RayInfo%nAziAngle
nPolarAng = RayInfo%nPolarAngle

! Geometry Info Pointing
Asy     => CoreInfo%Asy
Pin     => CoreInfo%Pin
PinInfo => CoreInfo%Pininfo
Cell    => CoreInfo%CellInfo

! Tracking Dat Pointing
phia          => TrackingDat%phianm
SrcAng        => TrackingDat%SrcAngnm
xst           => TrackingDat%xstnm
PhiAngIn      => TrackingDat%PhiAngInnm
DcmpPhiAngIn  => TrackingDat%DcmpPhiAngIn
DcmpPhiAngOut => TrackingDat%DcmpPhiAngOut
EXPA          => TrackingDat%EXPA
EXPB          => TrackingDat%EXPB
wtang         => TrackingDat%wtang
AziMap        => TrackingDat%AziMap

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

nAsyRay     = DcmpAsyRay%nAsyRay
iRotRay     = DcmpAsyRay%iRotRay
AsyRayList => DcmpAsyRay%AsyRayList
DirList    => DcmpAsyRay%DirList
AziList    => DcmpAsyRay%AziList
iAsy        = DcmpAsyRay%iAsy
iRay        = DcmpAsyRay%iRay
! ----------------------------------------------------
DO irot = 1, 2  
  PhiAngInSvIdx = RayInfo%PhiAngInSvIdx(iRotRay, irot)
  
  IF (DcmpAsyRay%lRotRayBeg(irot)) THEN
    PhiAngOut = PhiAngIn(1:nPolarAng, gb:ge, PhiAnginSvIdx)
  ELSE
    PhiAngOut = DcmpPhiAngIn(1:nPolarAng, gb:ge, irot, iRay, iAsy)
  END IF
  
  jbeg = 1; jend = nAsyRay; jinc = 1
  IF (irot .EQ. 2) THEN
    jbeg = nAsyRay; jend = 1; jinc = -1
  END IF
  
  DO j = jbeg, jend, jinc
    iAsyRay = AsyRayList(j)
    nPinRay = AsyRay(iAsyRay)%nCellRay
    idir    = DirList(j)
    iazi    = AziList(j)
    IF (irot .eq. 2) idir = mp(idir)
    
    AziSvIdx = AziMap(iAzi, :)
    
    DO ipol = 1, nPolarAng
      wt(ipol) = wtang(ipol, iazi)
      
      IF (lJout) wt2(ipol, 1:4) = wtang2(ipol, iazi, 1:4)
    END DO
    
    IF (idir .EQ. 1) THEN
      DO l = 1, nPinRay
        ipin   = AsyRay(iAsyRay)%PinIdx(l)
        iceray = AsyRay(iAsyRay)%PinRayIdx(l)
        
        ipin = Asy(iAsy)%GlobalPinIdx(ipin)   
        
        icel     = Pin(ipin)%Cell(iz)
        FsrIdxSt = Pin(ipin)%FsrIdxSt
        
        ibcel    = Cell(icel)%basecellstr
        CellRay => Cell(ibcel)%CellRay(iceray)
        nRaySeg  = CellRay%nSeg
        
        LocalFsrIdx => CellRay%LocalFsrIdx
        LenSeg      => CellRay%LenSeg
        
        IF (lJout) THEN
          isurf = AsyRay(iAsyRay)%PinRaySurf(1, l)
          
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(1, ig, isurf, ipin) = Jout(1, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
              Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
            END DO
          END DO
        END IF
        
        DO iRaySeg = 1, nRaySeg
          ireg = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
          
          DO ig = gb, ge
            tau = - LenSeg(iRaySeg) * xst(ig, ireg)
            
            ExpAppIdx = max(INT(tau), -40000)
            ExpAppIdx = min(0, ExpAppIdx)
            
            DO ipol = 1, nPolarAng
              ExpApp = EXPA(ipol, ExpAppIdx) * tau + EXPB(ipol, ExpAppIdx)
              
              phid = (PhiAngOut(ipol, ig) - SrcAng(ipol, ig, ireg, AziSvIdx(idir))) * ExpApp
              
              PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
              
              phia(ipol, ig, ireg, AziSvIdx(idir)) = phia(ipol, ig, ireg, AziSvIdx(idir)) + phid
            END DO
          END DO
        END DO
        
        IF (lJout) THEN
          isurf = AsyRay(iAsyRay)%PinRaySurf(2, l)
          
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(2, ig, isurf, ipin) = Jout(2, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
              Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
            END DO
          END DO
        END IF
      END DO
    ELSE
      DO l = nPinRay, 1, -1
        ipin   = AsyRay(iAsyRay)%PinIdx(l)
        iceray = AsyRay(iAsyRay)%PinRayIdx(l)
        
        ipin = Asy(iAsy)%GlobalPinIdx(ipin)   
        
        icel     = Pin(ipin)%Cell(iz)             
        FsrIdxSt = Pin(ipin)%FsrIdxSt
        
        ibcel    = Cell(icel)%basecellstr
        CellRay => Cell(ibcel)%CellRay(iceray)
        nRaySeg  = CellRay%nSeg
        
        LocalFsrIdx => CellRay%LocalFsrIdx
        LenSeg      => CellRay%LenSeg
        
        IF (lJout) THEN
          isurf = AsyRay(iAsyRay)%PinRaySurf(2, l)
          
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(1, ig, isurf, ipin) = Jout(1, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
              Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
            END DO
          END DO
        END IF
        
        DO iRaySeg = nRaySeg, 1, -1
          ireg = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
          
          DO ig = gb, ge
            tau = - LenSeg(iRaySeg) * xst(ig, ireg)
            
            ExpAppIdx = max(INT(tau), -40000)
            ExpAppIdx = min(0, ExpAppIdx)
            
            DO ipol = 1, nPolarAng
              ExpApp = EXPA(ipol, ExpAppIdx) * tau + EXPB(ipol, ExpAppIdx)
              
              phid = (PhiAngOut(ipol, ig) - SrcAng(ipol, ig, ireg, AziSvIdx(idir))) * ExpApp
              
              PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
              
              phia(ipol, ig, ireg, AziSvIdx(idir)) = phia(ipol, ig, ireg, AziSvIdx(idir)) + phid
            END DO
          END DO
        END DO
        
        IF (lJout) THEN
          isurf = AsyRay(iAsyRay)%PinRaySurf(1, l)
          
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(2, ig, isurf, ipin) = Jout(2, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
              Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
            END DO
          END DO
        END IF
      END DO
    END IF
  END DO
  
  DcmpPhiAngOut(1:nPolarAng, gb:ge, irot, iRay, iAsy) = PhiAngOut
END DO
! ----------------------------------------------------

END SUBROUTINE TrackRotRayPn_Dcmp
! ------------------------------------------------------------------------------------------------------------