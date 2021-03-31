#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTrace_Dcmp(RayInfo, CoreInfo, iz, gb, ge, lJout, lSubGrp, lScat1, lLinSrcCASMO, lHybrid)

USE PARAM
USE TYPEDEF,  ONLY : RayInfo_Type,      Coreinfo_type,   PolarAngle_Type,    AziAngleInfo_Type,                     &
                     Pin_Type,          Asy_Type,        AsyInfo_Type
USE Moc_Mod,  ONLY : nMaxRaySeg,        nMaxCellRay,     nMaxAsyRay,         nMaxCoreRay,                           &
                     Expa_p,            Expb_p,          ApproxExp,          TrackingDat,                           &
                     wtang,             Comp,            mwt,                AziMap,                                &
                     phisnm,            phimnm,          srcnm,              srcmnm,                                &
                     xstnm,             MocJoutnm,       PhiAngInnm,         DcmpPhiAngIn,                          &
                     DcmpPhiAngOut,     RayTraceDcmp_NM, RayTraceDcmp_Pn,    RayTraceDcmp_LSCASMO,                  &
                     DcmpGatherBoundaryFlux,             DcmpScatterBoundaryFlux,                                   &
                     DcmpLinkBoundaryFlux
USE Core_mod, ONLY : phisSlope,         srcSlope
USE PE_MOD,   ONLY : PE
USE GEOM,     ONLY : ng
USE CNTL,     ONLY : nTracerCntl
USE itrcntl_mod,    ONLY : itrcntl
USE ALLOCS
USE OMP_LIB
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
INTEGER :: iz, gb, ge
LOGICAL :: lJout, lSubGrp, lScat1, lLinSrcCASMO, lHybrid

TYPE(PolarAngle_Type), POINTER :: PolarAng(:)
TYPE(AziAngleInfo_Type), POINTER :: AziAng(:)
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)

LOGICAL, SAVE :: lFirst
DATA lFirst /.TRUE./

INTEGER :: color, tid
INTEGER :: nThread, nPolarAngle, nAziAngle, ScatOd
INTEGER :: AsyType, iAsy, iAzi
INTEGER :: ipol, od
INTEGER :: startColor, endColor, colorInc
REAL :: wtcos, wtpolar, wtsin2, wttemp

PolarAng => RayInfo%PolarAngle
AziAng => RayInfo%AziAngle
AsyInfo => CoreInfo%AsyInfo
Asy => CoreInfo%Asy

nThread = PE%nThread
nPolarAngle = RayInfo%nPolarAngle
nAziAngle = RayInfo%nAziAngle
ScatOd = nTracerCntl%ScatOd

startColor = red
endColor = black
colorInc = black - red
IF (mod(itrcntl%mocit, 2) .EQ. 0) THEN
  startColor = black
  endColor = red
  colorInc = red - black
ENDIF

CALL OMP_SET_NUM_THREADS(nThread)

IF (lFirst) THEN
  lFirst = FALSE
  CALL ApproxExp(RayInfo%PolarAngle, nPolarAngle)
  DO tid = 1, nThread
    IF (TrackingDat(tid)%lAllocNM) CYCLE
    TrackingDat(tid)%Expa => Expa_p; TrackingDat(tid)%Expb => Expb_p
    TrackingDat(tid)%srcnm => srcnm; TrackingDat(tid)%xstnm => xstnm
    TrackingDat(tid)%PhiAngInnm => PhiAngInnm
    TrackingDat(tid)%DcmpPhiAngIn => DcmpPhiAngIn
    TrackingDat(tid)%DcmpPhiAngOut => DcmpPhiAngOut
    TrackingDat(tid)%lAllocNM = TRUE
  ENDDO
  IF (lLinSrcCASMO) THEN
    DO tid = 1, nThread
      IF (TrackingDat(tid)%lAllocLinSrc) CYCLE
      CALL Dmalloc(TrackingDat(tid)%FsrIdx, nMaxRaySeg, nMaxCoreRay)
      CALL Dmalloc(TrackingDat(tid)%ExpAppIdxnm, ng, nMaxRaySeg, nMaxCoreRay)
      CALL Dmalloc(TrackingDat(tid)%OptLenListnm, ng, nMaxRaySeg, nMaxCoreRay)
      CALL Dmalloc(TrackingDat(tid)%ExpAppnm, nPolarAngle, ng, nMaxRaySeg, nMaxCoreRay)
      CALL Dmalloc(TrackingDat(tid)%PhiAngOutnm, nPolarAngle, ng, nMaxRaySeg + 2)
      CALL Dmalloc(TrackingDat(tid)%cmOptLen, nPolarAngle, ng, nMaxRaySeg, nMaxCoreRay)
      CALL Dmalloc(TrackingDat(tid)%cmOptLenInv, nPolarAngle, ng, nMaxRaySeg, nMaxCoreRay)
      CALL Dmalloc(TrackingDat(tid)%q0, ng, nMaxRaySeg, nMaxCoreRay)
      CALL Dmalloc(TrackingDat(tid)%q1, nPolarAngle, ng, nMaxRaySeg, nMaxCoreRay)
      CALL Dmalloc(TrackingDat(tid)%E1, nPolarAngle, ng, nMaxRaySeg, nMaxCoreRay)
      CALL Dmalloc(TrackingDat(tid)%E3, nPolarAngle, ng, nMaxRaySeg, nMaxCoreRay)
      CALL Dmalloc(TrackingDat(tid)%R1, nPolarAngle, ng, nMaxRaySeg, nMaxCoreRay)
      CALL Dmalloc(TrackingDat(tid)%R3, nPolarAngle, ng, nMaxRaySeg, nMaxCoreRay)
      CALL Dmalloc(TrackingDat(tid)%x0, 2, nMaxRaySeg, nMaxCoreRay)
      CALL Dmalloc(TrackingDat(tid)%y0, 2, nMaxRaySeg, nMaxCoreRay)
      TrackingDat(tid)%lAllocLinSrc = TRUE
    ENDDO
  ENDIF
  IF (lScat1) THEN
    IF (ScatOd .EQ. 1) od = 2
    IF (ScatOd .EQ. 2) od = 5
    IF (ScatOd .EQ. 3) od = 9
    ALLOCATE(Comp(od, nPolarAngle, nAziAngle))
    ALLOCATE(mwt(od, nPolarAngle, nAziAngle))
    ALLOCATE(wtang(nPolarAngle, nAziAngle))
    DO iAzi = 1, nAziAngle / 2
      AziMap(iAzi, 1) = 1
      AziMap(iAzi, 2) = 2
      AziMap(nAziAngle - iAzi + 1, 1) = 3
      AziMap(nAziAngle - iAzi + 1, 2) = 4
    ENDDO
    DO ipol = 1, nPolarAngle
      wttemp = PolarAng(ipol)%weight * PolarAng(ipol)%sinv
      DO iazi = 1, nAziAngle
        wtang(ipol, iazi) = wttemp * AziAng(iazi)%weight * AziAng(iazi)%del
      ENDDO
    ENDDO
    DO ipol = 1, nPolarAngle
      wttemp = PolarAng(ipol)%sinv
      DO iazi = 1, nAziAngle
        Comp(1, ipol, iazi) = wttemp * AziAng(iazi)%cosv
        Comp(2, ipol, iazi) = wttemp * AziAng(iazi)%sinv
        mwt(1:2, ipol, iazi) = Comp(1:2, ipol, iazi) * wtang(ipol, iazi)      
      ENDDO
    ENDDO  
    IF (ScatOd .GE. 2) THEN
      DO ipol = 1, nPolarAngle
        wttemp = PolarAng(ipol)%sinv
        wtsin2 = PolarAng(ipol)%sinv * PolarAng(ipol)%sinv
        wtcos = PolarAng(ipol)%cosv
        wtpolar = 1.5_8 * PolarAng(ipol)%cosv * PolarAng(ipol)%cosv - 0.5_8    
        DO iazi = 1, nAziAngle
          Comp(3, ipol, iazi) = wtpolar
          Comp(4, ipol, iazi) = wtsin2 * (1._8 - 2._8 * AziAng(iazi)%sinv * AziAng(iazi)%sinv)
          Comp(5, ipol, iazi) = wtsin2 * (2._8 * AziAng(iazi)%sinv * AziAng(iazi)%cosv)
          mwt(3, ipol, iazi) = Comp(3, ipol, iazi) * wtang(ipol, iazi)
          mwt(4:5, ipol, iazi) = 0.75_8 * Comp(4:5, ipol, iazi) * wtang(ipol, iazi)              
        ENDDO
      ENDDO 
    ENDIF
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
        ENDDO
      ENDDO
    ENDIF  
    DO tid = 1, nThread
      TrackingDat(tid)%wtang => wtang
      TrackingDat(tid)%AziMap => AziMap
    ENDDO
  ENDIF
ENDIF

DO tid = 1, nThread
  IF (lLinSrcCASMO) TrackingDat(tid)%srcSlope => srcSlope(:, :, :, iz)
ENDDO

DcmpPhiAngOut(:, gb : ge, :, :, :) = zero

DO color = startColor, endColor, colorInc
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) THEN
    CALL DcmpScatterBoundaryFlux(RayInfo, PhiAngInnm, DcmpPhiAngIn)
  ENDIF
#endif
  !$OMP PARALLEL DO PRIVATE(iAsy, AsyType) SCHEDULE(DYNAMIC)
  DO iAsy = PE%myAsyBeg, PE%myAsyEnd
    IF (Asy(iAsy)%color .NE. color) CYCLE
    AsyType = Asy(iAsy)%AsyType
    IF (.NOT. lScat1 .OR. lSubGrp) THEN
      IF (AsyInfo(AsyType)%lFuel) THEN
        IF (lLinSrcCASMO) THEN
          IF (lHybrid) THEN
            CALL RayTraceDcmp_NM(RayInfo, CoreInfo, phisnm, PhiAngInnm, xstnm, srcnm,                               &
                                 MocJoutnm, iz, iAsy, gb, ge, ljout)
          ELSE
            CALL RayTraceDcmp_LSCASMO(RayInfo, CoreInfo, phisnm, phisSlope, PhiAngInnm, srcnm, srcSlope,            &
                                      xstnm, MocJoutnm, iz, iAsy, gb, ge, ljout)
          ENDIF
        ELSE
          CALL RayTraceDcmp_NM(RayInfo, CoreInfo, phisnm, PhiAngInnm, xstnm, srcnm,                                 &
                               MocJoutnm, iz, iAsy, gb, ge, ljout)
        ENDIF
      ELSE
        IF (lLinSrcCASMO) THEN
          CALL RayTraceDcmp_LSCASMO(RayInfo, CoreInfo, phisnm, phisSlope, PhiAngInnm, srcnm, srcSlope,              &
                                    xstnm, MocJoutnm, iz, iAsy, gb, ge, ljout)
        ELSE
          CALL RayTraceDcmp_NM(RayInfo, CoreInfo, phisnm, PhiAngInnm, xstnm, srcnm,                                 &
                               MocJoutnm, iz, iAsy, gb, ge, ljout)
        ENDIF
      ENDIF
    ELSE
      CALL RayTraceDcmp_Pn(RayInfo, CoreInfo, phisnm, phimnm, PhiAngInnm, xstnm, srcnm, srcmnm,                     &
                           MocJoutnm, iz, iAsy, gb, ge, ScatOd, lJout)
    ENDIF
  ENDDO
  !$OMP END PARALLEL DO
#ifdef MPI_ENV
  IF (PE%nRTProc .GT. 1) THEN
    CALL DcmpGatherBoundaryFlux(RayInfo, DcmpPhiAngOut)
  ENDIF
#endif
  IF (PE%RTMASTER) THEN
    CALL DcmpLinkBoundaryFlux(CoreInfo, RayInfo, PhiAngInnm, DcmpPhiAngIn, DcmpPhiAngOut, gb, ge, color)
  ENDIF
ENDDO
      
END SUBROUTINE RayTrace_Dcmp
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceDcmp_NM(RayInfo, CoreInfo, phisnm, PhiAngInnm, xstnm, srcnm,               &
                           joutnm, iz, iasy, gb, ge, ljout)
USE PARAM
USE TYPEDEF, ONLY :   RayInfo_Type,      Coreinfo_type,   AsyInfo_Type,                       &
                      Asy_Type,          Pin_Type,        Cell_Type,                          &       
                      TrackingDat_Type,  DcmpAsyRayInfo_Type
USE Moc_Mod, ONLY :   nMaxRaySeg,        nMaxCellRay,     nMaxAsyRay,       nMaxCoreRay,      &
                      Expa,              Expb,            TrackingDat,                        &
                      ApproxExp,         TrackRotRayNM_Dcmp
USE TIMER
USE geom,    ONLY :   ng
USE ALLOCS
USE PE_MOD,  ONLY :   PE
USE OMP_LIB
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phisnm(:, :), PhiAngInnm(:, :, :), xstnm(:, :), srcnm(:, :), joutnm(:, :, :, :)
INTEGER :: iz, iasy, gb, ge
LOGICAL :: ljout

TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(DcmpAsyRayInfo_Type), POINTER :: DcmpAsyRay(:, :)
INTEGER, POINTER :: DcmpAsyRayCount(:)

INTEGER :: AsyType

INTEGER :: nxy, nThread
INTEGER :: iAsyRay, iAzi
INTEGER :: tid
INTEGER :: FsrIdxSt, FsrIdxEnd
INTEGER :: PinIdxSt, PinIdxEnd
INTEGER :: icel, ireg
INTEGER :: i, j, l, g

AsyInfo => CoreInfo%AsyInfo; Asy => CoreInfo%Asy
Cell => CoreInfo%CellInfo; Pin => CoreInfo%Pin

AsyType = Asy(iAsy)%AsyType
nxy = AsyInfo(AsyType)%nxy
nThread = PE%nThread

PinIdxSt = Asy(iAsy)%GlobalPinIdx(1)
PinIdxEnd = Asy(iAsy)%GlobalPinIdx(nxy)
FsrIdxSt = Pin(PinIdxSt)%FsrIdxSt
FsrIdxEnd = Pin(PinIdxEnd)%FsrIdxSt + Cell(Pin(PinIdxEnd)%Cell(iz))%nFsr - 1

DcmpAsyRay => RayInfo%DcmpAsyRay
DcmpAsyRayCount => RayInfo%DcmpAsyRayCount

phisnm(gb : ge, FsrIdxSt : FsrIdxEnd) = zero
IF(ljout) joutnm(:, gb : ge, :, PinIdxSt : PinIdxEnd) = zero
    
tid = omp_get_thread_num() + 1
DO iAsyRay = 1, DcmpAsyRayCount(iAsy)
  CALL TrackRotRayNM_Dcmp(RayInfo, CoreInfo, TrackingDat(tid), phisnm, joutnm,   &
                          ljout, DcmpAsyRay(iAsyRay, iAsy), iz, gb, ge)
ENDDO

DO l = PinIdxSt, PinIdxEnd
  FsrIdxSt = Pin(l)%FsrIdxSt; icel = Pin(l)%Cell(iz);
  DO j = 1, Cell(icel)%nFsr
    ireg = FsrIdxSt + j - 1
    DO g = gb, ge
      phisnm(g, ireg) = phisnm(g, ireg) / (xstnm(g, ireg) * Cell(icel)%vol(j)) + srcnm(g, ireg)
    ENDDO
  ENDDO  
ENDDO

END SUBROUTINE RayTraceDcmp_NM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE TrackRotRayNM_Dcmp(RayInfo, CoreInfo, TrackingDat, phis, jout, ljout, DcmpAsyRay, iz, gb, ge)
USE PARAM
USE TYPEDEF,  ONLY :  RayInfo_Type,       Coreinfo_type,                                                &
                      Pin_Type,           Asy_Type,          PinInfo_Type,        Cell_Type,            &
                      AziAngleInfo_Type,  PolarAngle_Type,   AsyRayInfo_type,     RotRayInfo_Type,      &
                      CellRayInfo_type,   TrackingDat_Type,  DcmpAsyRayInfo_Type
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
TYPE(TrackingDat_Type) :: TrackingDat
TYPE(DcmpAsyRayInfo_Type) :: DcmpAsyRay
REAL, POINTER :: phis(:, :), jout(:, :, :, :)
LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: iz, gb, ge

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(AziAngleInfo_Type), POINTER :: AziAng(:)
TYPE(PolarAngle_Type), POINTER :: PolarAng(:)
TYPE(AsyRayInfo_type), POINTER :: AsyRay(:)
TYPE(CellRayInfo_Type), POINTER :: CellRay

REAL, POINTER :: LenSeg(:)
INTEGER, POINTER :: LocalFsrIdx(:)
REAL, POINTER :: src(:, :), xst(:, :)
REAL, POINTER :: PhiAngIn(:, :, :)
REAL, POINTER :: DcmpPhiAngIn(:, :, :, :, :), DcmpPhiAngOut(:, :, :, :, :)
REAL, POINTER :: EXPA(:, :), EXPB(:, :)
INTEGER, POINTER :: AsyRayList(:), DirList(:), AziList(:)

INTEGER :: mp(2)
INTEGER :: iazi, ipol, irotray, iasyray, iceray, iray, irayseg, irot, idir
INTEGER :: nAsyRay, nPinRay, nRaySeg, FsrIdxSt
INTEGER :: nAziAng, nPolarAng
INTEGER :: ipin, icel, ibcel, iasy, ireg, isurf
INTEGER :: PhiAnginSvIdx
INTEGER :: i, j, k, l, m, ig, ir, ir1
INTEGER :: jbeg, jend, jinc

REAL :: wttemp, wtang(10, 100), wtang2(10, 100, 4), wt(10), wt2(10, 4)

REAL :: PhiAngOut(RayInfo%nPolarAngle, gb : ge)
REAL :: phid, tau, ExpApp
INTEGER :: ExpAppIdx

DATA mp /2, 1/

!Ray Info Pointing
AziAng => RayInfo%AziAngle
PolarAng => RayInfo%PolarAngle
AsyRay => RayInfo%AsyRay

!Geometry Info Pointing
Asy => CoreInfo%Asy
Pin => CoreInfo%Pin
Cell => CoreInfo%CellInfo

!Tracking Dat Pointing
src => TrackingDat%srcnm
xst => TrackingDat%xstnm
PhiAngIn => TrackingDat%phiAngInnm
DcmpPhiAngIn => TrackingDat%DcmpPhiAngIn
DcmpPhiAngOut => TrackingDat%DcmpPhiAngOut
EXPA => TrackingDat%EXPA
EXPB => TrackingDat%EXPB

nPolarAng = RayInfo%nPolarAngle
nAziAng = RayInfo%nAziAngle

DO ipol = 1, nPolarAng
  wttemp = PolarAng(ipol)%weight * PolarAng(ipol)%sinv
  DO iazi = 1, nAziAng
    wtang(ipol, iazi) = wttemp * AziAng(iazi)%weight * AziAng(iazi)%del
  ENDDO
ENDDO

IF (lJout) THEN
  DO ipol = 1, nPolarAng
    DO iazi = 1, nAziAng
      wtang2(ipol, iazi, 1) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv)
      wtang2(ipol, iazi, 3) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv)
      wtang2(ipol, iazi, 2) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv)
      wtang2(ipol, iazi, 4) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv)
    ENDDO
  ENDDO
ENDIF

nAsyRay = DcmpAsyRay%nAsyRay
iRotRay = DcmpAsyRay%iRotRay
AsyRayList => DcmpAsyRay%AsyRayList
DirList => DcmpAsyRay%DirList
AziList => DcmpAsyRay%AziList
iAsy = DcmpAsyRay%iAsy; iRay = DcmpAsyRay%iRay

DO irot = 1, 2  
  PhiAngInSvIdx = RayInfo%PhiAngInSvIdx(iRotRay, irot)
  IF (DcmpAsyRay%lRotRayBeg(irot)) THEN
    PhiAngOut = PhiAngIn(1 : nPolarAng, gb : ge, PhiAnginSvIdx)
  ELSE
    PhiAngOut = DcmpPhiAngIn(1 : nPolarAng, gb : ge, irot, iRay, iAsy)
  ENDIF
  jbeg = 1; jend = nAsyRay; jinc = 1
  IF (irot .EQ. 2) THEN
    jbeg = nAsyRay; jend = 1; jinc = -1
  ENDIF
  DO j = jbeg, jend, jinc
    iAsyRay = AsyRayList(j)
    nPinRay = AsyRay(iAsyRay)%nCellRay
    idir = DirList(j); iazi = AziList(j)
    IF (irot .eq. 2) idir = mp(idir)
    DO ipol = 1, nPolarAng
      wt(ipol) = wtang(ipol, iazi)
      IF (lJout) wt2(ipol, 1 : 4) = wtang2(ipol, iazi, 1 : 4)
    ENDDO
    IF (idir .EQ. 1) THEN
      DO l = 1, nPinRay
        ipin = AsyRay(iAsyRay)%PinIdx(l)
        iceray = AsyRay(iAsyRay)%PinRayIdx(l)
        ipin = Asy(iAsy)%GlobalPinIdx(ipin)   
        icel = Pin(ipin)%Cell(iz)
        FsrIdxSt = Pin(ipin)%FsrIdxSt
        ibcel = Cell(icel)%basecellstr
        CellRay => Cell(ibcel)%CellRay(iceray)
        nRaySeg = CellRay%nSeg
        LocalFsrIdx => CellRay%LocalFsrIdx
        LenSeg => CellRay%LenSeg
        IF (lJout) THEN
          isurf = AsyRay(iAsyRay)%PinRaySurf(1, l)
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(1, ig, isurf, ipin) = Jout(1, ig, isurf, ipin) + wt(ipol) * PhiAngOut(ipol, ig)
              Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol, ig)
            ENDDO
          ENDDO
        ENDIF
        DO iRaySeg = 1, nRaySeg
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
            ENDDO
          ENDDO
        ENDDO
        IF (lJout) THEN
          isurf = AsyRay(iAsyRay)%PinRaySurf(2, l)
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(2, ig, isurf, ipin) = Jout(2, ig, isurf, ipin) + wt(ipol) * PhiAngOut(ipol, ig)
              Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol, ig)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ELSE
      DO l = nPinRay, 1, -1
        ipin = AsyRay(iAsyRay)%PinIdx(l)
        iceray = AsyRay(iAsyRay)%PinRayIdx(l)
        ipin = Asy(iAsy)%GlobalPinIdx(ipin)   
        icel = Pin(ipin)%Cell(iz)             
        FsrIdxSt = Pin(ipin)%FsrIdxSt
        ibcel = Cell(icel)%basecellstr
        CellRay => Cell(ibcel)%CellRay(iceray)
        nRaySeg = CellRay%nSeg
        LocalFsrIdx => CellRay%LocalFsrIdx
        LenSeg => CellRay%LenSeg
        IF (lJout) THEN
          isurf = AsyRay(iAsyRay)%PinRaySurf(2, l)
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(1, ig, isurf, ipin) = Jout(1, ig, isurf, ipin) + wt(ipol) * PhiAngOut(ipol, ig)
              Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol, ig)
            ENDDO
          ENDDO
        ENDIF
        DO iRaySeg = nRaySeg, 1, -1
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
            ENDDO
          ENDDO
        ENDDO
        IF (lJout) THEN
          isurf = AsyRay(iAsyRay)%PinRaySurf(1, l)
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(2, ig, isurf, ipin) = Jout(2, ig, isurf, ipin) + wt(ipol) * PhiAngOut(ipol, ig)
              Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol, ig)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  DcmpPhiAngOut(1 : nPolarAng, gb : ge, irot, iRay, iAsy) = PhiAngOut
ENDDO

END SUBROUTINE TrackRotRayNM_Dcmp
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
ENDDO
    
DO i = FsrIdxSt, FsrIdxEnd
  DO g = gb, ge
    phimx(1, g, i) = phimx(1, g, i) / xstnm(g, i)
    phimy(1, g, i) = phimy(1, g, i) / xstnm(g, i)
  ENDDO
ENDDO

DO l = PinIdxSt, PinIdxEnd
  icel = Pin(l)%Cell(iz);
  DO j = 1, Cell(icel)%nFsr
    ireg = Pin(l)%FsrIdxSt + j - 1
    DO g = gb, ge
      phisnm(g, ireg) = phisnm(g, ireg) / (xstnm(g, ireg) * Cell(icel)%vol(j)) + srcnm(g, ireg)
      phimx(:, g, ireg) = phimx(:, g, ireg) / (xstnm(g, ireg) * Cell(icel)%vol(j))
      phimy(:, g, ireg) = phimy(:, g, ireg) / (xstnm(g, ireg) * Cell(icel)%vol(j))
    ENDDO
  ENDDO  
ENDDO

DO i = FsrIdxSt, FsrIdxEnd
  detinv = 1.0 / (FsrMxx(i) * FsrMyy(i) - FsrMxy(i) * FsrMxy(i))
  DO g = gb, ge
    phisSlope(1, g, i, iz) = detinv * (FsrMyy(i) * SUM(phimx(:, g, i)) - FsrMxy(i) * SUM(phimy(:, g, i)))
    phisSlope(2, g, i, iz) = detinv * (FsrMxx(i) * SUM(phimy(:, g, i)) - FsrMxy(i) * SUM(phimx(:, g, i)))
  ENDDO
ENDDO

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
  ENDDO
ENDDO
IF(lJout)THEN
  DO ipol = 1, nPolarAng
    DO iazi = 1, nAziAng
      wtang2(ipol, iazi, 1) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv)
      wtang2(ipol, iazi, 3) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv)
      wtang2(ipol, iazi, 2) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv)
      wtang2(ipol, iazi, 4) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv)
    ENDDO
  ENDDO
ENDIF

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
  ENDDO
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
        ENDDO
        ExpAppIdx(ig, irsegidx, j) = max(INT(tau), -40000)
        ExpAppIdx(ig, irsegidx, j) = min(0, ExpAppIdx(ig, irsegidx, j))
      ENDDO
    ENDDO
    CellRayIdxSt(l, j, 1) = irsegidx
    SurfIdx(l, j, 1) = AsyRay(iAsyRay)%PinRaySurf(2, l) !OutSurface
    SurfIdx(l, j, 2) = AsyRay(iAsyRay)%PinRaySurf(1, l) !Insurface
  ENDDO
  nTotRaySeg(j) = irsegidx
  nTotCellRay(j) = nPinRay
ENDDO

DO ipol = 1, nPolarAng
  DO j = 1, nAsyRay
    DO l = 1, nTotRaySeg(j)
      DO ig = gb, ge
        E1(ipol, ig, l, j) = EXPA(ExpAppIdx(ig, l, j), ipol) * OptLenList(ig, l, j) + EXPB(ExpAppIdx(ig, l, j), ipol)
        E3(ipol, ig, l, j) = 2._8 * (cmOptLen(ipol, ig, l, j) - E1(ipol, ig, l, j)) - cmOptLen(ipol, ig, l, j) * E1(ipol, ig, l, j)
        R1(ipol, ig, l, j) = 1._8 + half * cmOptLen(ipol, ig, l, j) - (1._8 + cmOptLenInv(ipol, ig, l, j)) * E1(ipol, ig, l, j)
        R3(ipol, ig, l, j) = rsix * cmOptLen(ipol, ig, l, j) - 2._8 * cmOptLenInv(ipol, ig, l, j) - 2._8  &
                             + (1._8 + cmOptLenInv(ipol, ig, l, j)) * (1._8 + 2._8 * cmOptLenInv(ipol, ig, l, j)) * E1(ipol, ig, l, j)
      ENDDO
    ENDDO
  ENDDO
ENDDO

DO irot = 1, 2
  IF (DcmpAsyRay%lRotRayBeg(irot)) THEN
    phiobd(1 : nPolarAng, gb : ge) = PhiAngIn(1 : nPolarAng, gb : ge, PhiAnginSvIdx(irot))
  ELSE
    phiobd(1 : nPolarAng, gb : ge) = DcmpPhiAngIn(1 : nPolarAng, gb : ge, irot, iRay, iAsy)
  ENDIF
  jinc = 1; jbeg = 1; jend = nAsyRay
  IF(irot .eq. 2) THEN
    jinc = -1; jbeg = nAsyRay; jend = 1
  ENDIF
  DO j = jbeg, jend, jinc
    idir = DirList(j); iazi = AziList(j)
    nRaySeg = nTotRaySeg(j)
    IF(irot .eq. 2) idir = mp(idir)
    DO ipol = 1, nPolarAng
      wt(ipol) = wtang(ipol, iazi)
      ax(ipol) = wt(ipol) * AziAng(iazi)%cosv * PolarAng(ipol)%sinv
      ay(ipol) = wt(ipol) * AziAng(iazi)%sinv * PolarAng(ipol)%sinv
    ENDDO
    IF(lJout) THEN
      wt2(1 : nPolarAng, :) = wtang2(1 : nPolarAng, iazi, :)
    ENDIF
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
          ENDDO
        ENDDO
      ENDDO
      phiobd(1 : nPolarAng, gb : ge) = PhiAngOut(1 : nPolarAng, gb : ge, nRaySeg + 1)
      IF(ljout) THEN
        DO ir = 1, nTotCellRay(j)
          icel = PinIdx(ir, j); isurf = SurfIdx(ir, j, 1)
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(2, ig, isurf, icel) = Jout(2, ig, isurf, icel) + wt(ipol) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 1) + 1)
              Jout(3, ig, isurf, icel) = Jout(3, ig, isurf, icel) + wt2(ipol, isurf) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 1) + 1)
            ENDDO
          ENDDO
          isurf = SurfIdx(ir, j, 2)
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(1, ig, isurf, icel) = Jout(1, ig, isurf, icel) + wt(ipol) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 2))
              Jout(3, ig, isurf, icel) = Jout(3, ig, isurf, icel) + wt2(ipol, isurf) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 2))
            ENDDO
          ENDDO
        ENDDO
      ENDIF
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
          ENDDO
        ENDDO
      ENDDO
      phiobd(1 : nPolarAng, gb : ge) = PhiAngOut(1 : nPolarAng, gb : ge, 2)
      IF(lJout) THEN
        DO ir = 1, nTotCellRay(j)
          icel = PinIdx(ir, j); isurf = SurfIdx(ir, j, 2)
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(2, ig, isurf, icel) = Jout(2, ig, isurf, icel) + wt(ipol) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 2) + 1)
              Jout(3, ig, isurf, icel) = Jout(3, ig, isurf, icel) + wt2(ipol, isurf) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 2) + 1)
            ENDDO
          ENDDO
          isurf = SurfIdx(ir, j, 1)
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(1, ig, isurf, icel) = Jout(1, ig, isurf, icel) + wt(ipol) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 1) + 2)
              Jout(3, ig, isurf, icel) = Jout(3, ig, isurf, icel) + wt2(ipol, isurf) * PhiAngOut(ipol, ig, CellRayIdxSt(ir, j, 1) + 2)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  ENDDO
  DcmpPhiAngOut(1 : nPolarAng, gb : ge, irot, iRay, iAsy) = phiobd(1 : nPolarAng, gb : ge)
ENDDO

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

USE PARAM
USE TYPEDEF, ONLY :   RayInfo_Type,      Coreinfo_type,     TrackingDat_Type,                                           &
                      Pin_Type,          Cell_Type,         PolarAngle_Type,        AziAngleInfo_Type,                  &
                      Asy_Type,          AsyInfo_Type,      DcmpAsyRayInfo_Type
USE GEOM,    ONLY :   ng
USE Moc_Mod, ONLY :   nMaxRaySeg,        nMaxCellRay,       nMaxAsyRay,             nMaxCoreRay,                        &
                      Expa,              Expb,              TrackingDat,            ApproxExp,                          &
                      wtang,             Comp,              mwt,                    AziMap,                             &
                      DcmpPhiAngIn,      DcmpPhiAngOut,     TrackRotRayPn_Dcmp
USE TIMER
USE BasicOperation, ONLY : CP_CA, CP_VA
USE ALLOCS
USE PE_MOD,  ONLY :   PE
USE OMP_LIB
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phisnm(:, :), phimnm(:, :, :), PhiAngInnm(:, :, :), xstnm(:, :)
REAL, POINTER :: srcnm(:, :), srcmnm(:, :, :), joutnm(:, :, :, :)
INTEGER :: iz, iAsy, gb, ge
LOGICAL :: ljout
INTEGER :: ScatOd

LOGICAL, SAVE :: lfirst
DATA lfirst /.TRUE./

!Pointing Variable
TYPE(AziAngleInfo_Type), POINTER :: AziAng(:)
TYPE(PolarAngle_Type), POINTER :: PolarAng(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(DcmpAsyRayInfo_Type), POINTER :: DcmpAsyRay(:, :)
REAL, POINTER :: phianm(:, :, :, :), SrcAngnm(:, :, :, :)
INTEGER, POINTER :: DcmpAsyAziList(:, :, :)

INTEGER :: nAziAng, nPolarAng
INTEGER :: nAsyRay, nAsy, nxy

INTEGER :: tid
INTEGER :: FsrIdxSt, FsrIdxEnd, PinIdxSt, PinIdxEnd
INTEGER :: icel, ireg, ipin, iazi, ipol, iDcmpAsyRay, iCoreRay
INTEGER :: AziIdx
REAL :: wttemp, wtcos, wtpolar, wtsin2, tempsrc
INTEGER :: i, j, l, k, m, g

Asy => CoreInfo%Asy
AsyInfo => CoreInfo%AsyInfo
Cell => CoreInfo%CellInfo
Pin => CoreInfo%Pin
nAsy = CoreInfo%nxya

AziAng => RayInfo%AziAngle
PolarAng => RayInfo%PolarAngle
DcmpAsyRay => RayInfo%DcmpAsyRay
DcmpAsyAziList => RayInfo%DcmpAsyAziList
nAziAng = RayInfo%nAziAngle
nPolarAng = RayInfo%nPolarAngle

tid = omp_get_thread_num() + 1

nxy = AsyInfo(Asy(iAsy)%AsyType)%nxy
PinIdxSt = Asy(iAsy)%GlobalPinIdx(1)
PinIdxEnd = Asy(iAsy)%GlobalPinIdx(nxy)
FsrIdxSt = Pin(Asy(iAsy)%GlobalPinIdx(1))%FsrIdxSt
FsrIdxEnd = Pin(Asy(iAsy)%GlobalPinIdx(nxy))%FsrIdxSt + Cell(Pin(Asy(iAsy)%GlobalPinIdx(nxy))%Cell(iz))%nFsr - 1
   
CALL Dmalloc0(TrackingDat(tid)%phianm, 1, nPolarAng, gb, ge, FsrIdxSt, FsrIdxEnd, 1, 4)
CALL Dmalloc0(TrackingDat(tid)%SrcAngnm, 1, nPolarAng, gb, ge, FsrIdxSt, FsrIdxEnd, 1, 4)
phianm => TrackingDat(tid)%phianm
SrcAngnm => TrackingDat(tid)%SrcAngnm

phisnm(gb : ge, FsrIdxSt : FsrIdxEnd) = zero
phimnm(:, gb : ge, FsrIdxSt : FsrIdxEnd) = zero
IF (lJout) joutnm(:, gb : ge, :, PinIdxSt : PinIdxEnd) = zero

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
        ENDDO
      ENDDO
      AziIdx = nAziAng - iAzi + 1
      DO g = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = srcnm(g, ireg)
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = srcnm(g, ireg)
          tempsrc = Comp(1, ipol, AziIdx) * srcmnm(1, g, ireg) + Comp(2, ipol, AziIdx) * srcmnm(2, g, ireg)
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) - tempsrc
        ENDDO
      ENDDO
    ENDDO
  ELSEIF (ScatOd .EQ. 2) THEN 
    DO ireg = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      DO g = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = srcnm(g, ireg)
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = srcnm(g, ireg)
          tempsrc = Comp(1, ipol, AziIdx) * srcmnm(1, g, ireg) + Comp(2, ipol, AziIdx) * srcmnm(2, g, ireg)
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) - tempsrc
          tempsrc = Comp(3, ipol, AziIdx) * srcmnm(3, g, ireg) + Comp(4, ipol, AziIdx) * srcmnm(4, g, ireg)                        &
                    + Comp(5, ipol, AziIdx) * srcmnm(5, g, ireg)
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) + tempsrc
        ENDDO
      ENDDO
      AziIdx = nAziAng - iAzi + 1
      DO g = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = srcnm(g, ireg)
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = srcnm(g, ireg)
          tempsrc = Comp(1, ipol, AziIdx) * srcmnm(1, g, ireg) + Comp(2, ipol, AziIdx) * srcmnm(2, g, ireg)
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) - tempsrc
          tempsrc = Comp(3, ipol, AziIdx) * srcmnm(3, g, ireg) + Comp(4, ipol, AziIdx) * srcmnm(4, g, ireg)                        &
                    + Comp(5, ipol, AziIdx) * srcmnm(5, g, ireg)
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) + tempsrc
        ENDDO
      ENDDO
    ENDDO
  ELSEIF (ScatOd .EQ. 3) THEN
    DO ireg = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      DO g = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = srcnm(g, ireg)
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = srcnm(g, ireg)
          tempsrc = Comp(1, ipol, AziIdx) * srcmnm(1, g, ireg) + Comp(2, ipol, AziIdx) * srcmnm(2, g, ireg)                        &
                    + Comp(6, ipol, AziIdx) * srcmnm(6, g, ireg) + Comp(7, ipol, AziIdx) * srcmnm(7, g, ireg)                      &
                    + Comp(8, ipol, AziIdx) * srcmnm(8, g, ireg) + Comp(9, ipol, AziIdx) * srcmnm(9, g, ireg)
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) - tempsrc
          tempsrc = Comp(3, ipol, AziIdx) * srcmnm(3, g, ireg) + Comp(4, ipol, AziIdx) * srcmnm(4, g, ireg)                        &
                    + Comp(5, ipol, AziIdx) * srcmnm(5, g, ireg)
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) + tempsrc
        ENDDO
      ENDDO
      AziIdx = nAziAng - iAzi + 1
      DO g = gb, ge
        DO ipol = 1, nPolarAng
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = srcnm(g, ireg)
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = srcnm(g, ireg)
          tempsrc = Comp(1, ipol, AziIdx) * srcmnm(1, g, ireg) + Comp(2, ipol, AziIdx) * srcmnm(2, g, ireg)                        &
                    + Comp(6, ipol, AziIdx) * srcmnm(6, g, ireg) + Comp(7, ipol, AziIdx) * srcmnm(7, g, ireg)                      &
                    + Comp(8, ipol, AziIdx) * srcmnm(8, g, ireg) + Comp(9, ipol, AziIdx) * srcmnm(9, g, ireg)
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) - tempsrc
          tempsrc = Comp(3, ipol, AziIdx) * srcmnm(3, g, ireg) + Comp(4, ipol, AziIdx) * srcmnm(4, g, ireg)                        &
                    + Comp(5, ipol, AziIdx) * srcmnm(5, g, ireg)
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 1)) + tempsrc
          SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) = SrcAngnm(ipol, g, ireg, AziMap(AziIdx, 2)) + tempsrc
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  phianm = zero
  DO i = 1, DcmpAsyAziList(0, iAzi, iAsy)
    iDcmpAsyRay = DcmpAsyAziList(i, iAzi, iAsy)
    CALL TrackRotRayPn_Dcmp(RayInfo, CoreInfo, TrackingDat(tid), DcmpAsyRay(iDcmpAsyRay, iAsy), Joutnm, lJout, iz, gb, ge)
  ENDDO
        
  IF (ScatOd .EQ. 1) THEN
    DO ireg = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      DO g = gb, ge
        DO ipol = 1, nPolarAng            
          phisnm(g, ireg) = phisnm(g, ireg) + wtang(ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) + phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
          phimnm(1:2, g, ireg) = phimnm(1:2, g, ireg) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) - phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
        ENDDO
      ENDDO
      AziIdx = nAziAng - iAzi + 1
      DO g = gb, ge
        DO ipol = 1, nPolarAng            
          phisnm(g, ireg) = phisnm(g, ireg) + wtang(ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) + phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
          phimnm(1:2, g, ireg) = phimnm(1:2, g, ireg) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) - phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
        ENDDO
      ENDDO
    ENDDO
  ELSEIF (ScatOd .EQ. 2) THEN
    DO ireg = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      DO g = gb, ge
        DO ipol = 1, nPolarAng            
          phisnm(g, ireg) = phisnm(g, ireg) + wtang(ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) + phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
          phimnm(1:2, g, ireg) = phimnm(1:2, g, ireg) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) - phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
          phimnm(3:5, g, ireg) = phimnm(3:5, g, ireg) + mwt(3:5, ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) + phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
        ENDDO
      ENDDO
      AziIdx = nAziAng - iAzi + 1
      DO g = gb, ge
        DO ipol = 1, nPolarAng            
          phisnm(g, ireg) = phisnm(g, ireg) + wtang(ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) + phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
          phimnm(1:2, g, ireg) = phimnm(1:2, g, ireg) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) - phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
          phimnm(3:5, g, ireg) = phimnm(3:5, g, ireg) + mwt(3:5, ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) + phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
        ENDDO
      ENDDO
    ENDDO
  ELSEIF (ScatOd .EQ. 3) THEN
    DO ireg = FsrIdxSt, FsrIdxEnd
      AziIdx = iAzi
      DO g = gb, ge
        DO ipol = 1, nPolarAng            
          phisnm(g, ireg) = phisnm(g, ireg) + wtang(ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) + phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
          phimnm(1:2, g, ireg) = phimnm(1:2, g, ireg) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) - phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
          phimnm(3:5, g, ireg) = phimnm(3:5, g, ireg) + mwt(3:5, ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) + phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
          phimnm(6:9, g, ireg) = phimnm(6:9, g, ireg) + mwt(6:9, ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) - phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
        ENDDO
      ENDDO
      AziIdx = nAziAng - iAzi + 1
      DO g = gb, ge
        DO ipol = 1, nPolarAng            
          phisnm(g, ireg) = phisnm(g, ireg) + wtang(ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) + phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
          phimnm(1:2, g, ireg) = phimnm(1:2, g, ireg) + mwt(1:2, ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) - phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
          phimnm(3:5, g, ireg) = phimnm(3:5, g, ireg) + mwt(3:5, ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) + phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
          phimnm(6:9, g, ireg) = phimnm(6:9, g, ireg) + mwt(6:9, ipol, AziIdx) * (phianm(ipol, g, ireg, AziMap(AziIdx, 1)) - phianm(ipol, g, ireg, AziMap(AziIdx, 2)))
        ENDDO
      ENDDO
    ENDDO
  ENDIF
    
ENDDO
    
DEALLOCATE(TrackingDat(tid)%phianm)
DEALLOCATE(TrackingDat(tid)%SrcAngnm)

IF (ScatOd .EQ. 1) THEN
  !$OMP PARALLEL DO PRIVATE(FsrIdxSt, ireg, icel, wttemp) SCHEDULE(DYNAMIC)
  DO l = PinIdxSt, PinIdxEnd
    FsrIdxSt = Pin(l)%FsrIdxSt; icel = Pin(l)%Cell(iz);
    DO j = 1, Cell(icel)%nFsr
      ireg = FsrIdxSt + j - 1
      DO g = gb, ge
        wttemp = 1._8 / (xstnm(g, ireg) * Cell(icel)%vol(j))
        phisnm(g, ireg) = phisnm(g, ireg) * wttemp + srcnm(g, ireg)
        phimnm(1:2, g, ireg) = phimnm(1:2, g, ireg) * wttemp + srcmnm(1:2, g, ireg) * rthree
      ENDDO
    ENDDO  
  ENDDO
  !$OMP END PARALLEL DO
ELSEIF (ScatOd .EQ. 2) THEN
  !$OMP PARALLEL DO PRIVATE(FsrIdxSt, ireg, icel, wttemp) SCHEDULE(DYNAMIC)
  DO l = PinIdxSt, PinIdxEnd
    FsrIdxSt = Pin(l)%FsrIdxSt; icel = Pin(l)%Cell(iz);
    DO j = 1, Cell(icel)%nFsr
      ireg = FsrIdxSt + j - 1
      DO g = gb, ge
        wttemp = 1._8 / (xstnm(g, ireg) * Cell(icel)%vol(j))
        phisnm(g, ireg) = phisnm(g, ireg) * wttemp + srcnm(g, ireg)
        phimnm(1:2, g, ireg) = phimnm(1:2, g, ireg) * wttemp + srcmnm(1:2, g, ireg) * rthree
        phimnm(3:5, g, ireg) = phimnm(3:5, g, ireg) * wttemp + srcmnm(3:5, g, ireg) * rfive
      ENDDO
    ENDDO  
  ENDDO
  !$OMP END PARALLEL DO
ELSEIF (ScatOd .EQ. 3) THEN
  !$OMP PARALLEL DO PRIVATE(FsrIdxSt, ireg, icel, wttemp) SCHEDULE(DYNAMIC)
  DO l = PinIdxSt, PinIdxEnd
    FsrIdxSt = Pin(l)%FsrIdxSt; icel = Pin(l)%Cell(iz);
    DO j = 1, Cell(icel)%nFsr
      ireg = FsrIdxSt + j - 1
      DO g = gb, ge
        wttemp = 1._8 / (xstnm(g, ireg) * Cell(icel)%vol(j))
        phisnm(g, ireg) = phisnm(g, ireg) * wttemp + srcnm(g, ireg)
        phimnm(1:2, g, ireg) = phimnm(1:2, g, ireg) * wttemp + srcmnm(1:2, g, ireg) * rthree
        phimnm(3:5, g, ireg) = phimnm(3:5, g, ireg) * wttemp + srcmnm(3:5, g, ireg) * rfive
        phimnm(6:9, g, ireg) = phimnm(6:9, g, ireg) * wttemp + srcmnm(6:9, g, ireg) * rseven
      ENDDO
    ENDDO  
  ENDDO
  !$OMP END PARALLEL DO
ENDIF

END SUBROUTINE RayTraceDcmp_Pn
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE TrackRotRayPn_Dcmp(RayInfo, CoreInfo, TrackingDat, DcmpAsyRay, Jout, lJout, iz, gb, ge)
USE PARAM
USE TYPEDEF,  ONLY :  RayInfo_Type,         Coreinfo_type,      Pin_Type,            Asy_Type,                           &
                      AsyInfo_Type,         PinInfo_Type,       Cell_Type,           AziAngleInfo_Type,                  &
                      PolarAngle_Type,      ModRayInfo_type,    AsyRayInfo_type,     CoreRayInfo_Type,                   &
                      RotRayInfo_Type,      CellRayInfo_type,   TrackingDat_Type,    DcmpAsyRayInfo_Type
USE Moc_Mod, ONLY :   nMaxRaySeg,           nMaxCellRay,        nMaxAsyRay,          nMaxCoreRay
USE BasicOperation, ONLY : CP_CA, CP_VA
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
TYPE(TrackingDat_Type) :: TrackingDat
TYPE(DcmpAsyRayInfo_Type) :: DcmpAsyRay
REAL, POINTER :: Jout(:, :, :, :)
LOGICAL :: lJout
INTEGER :: iz, gb, ge, ScatOd

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
REAL, POINTER :: phia(:, :, :, :)
REAL, POINTER :: SrcAng(:, :, :, :), xst(:, :)
REAL, POINTER :: PhiAngIn(:, :, :)
REAL, POINTER :: DcmpPhiAngIn(:, :, :, :, :), DcmpPhiAngOut(:, :, :, :, :)
REAL, POINTER :: EXPA(:, :), EXPB(:, :)
REAL, POINTER :: wtang(:, :)
INTEGER, POINTER :: AziMap(:, :)
INTEGER, POINTER :: AsyRayList(:), DirList(:), AziList(:)

INTEGER :: mp(2)
INTEGER :: iAzi, ipol, iRotRay, iCoreRay, iAsyRay, iRay, iRaySeg, irot, idir
INTEGER :: nCoreRay, nAsyRay, nDcmpRay, nPinRay, nRaySeg, FsrIdxSt
INTEGER :: nPolarAng, nAziAng
INTEGER :: ipin, icel, iasy, ireg, isurf, ibcel, iceray
INTEGER :: i, j, k, l, m, ig, ir
INTEGER :: jbeg, jend, jinc
INTEGER :: PhiAnginSvIdx, AziSvIdx(2)

REAL :: wtang2(10, 100, 4), wt(10), wt2(10, 4)

REAL :: PhiAngOut(RayInfo%nPolarAngle, gb : ge)
REAL :: phid, tau, ExpApp
INTEGER :: ExpAppIdx

DATA mp /2, 1/

!Ray Info Pointing
AziAng => RayInfo%AziAngle
PolarAng => RayInfo%PolarAngle
AsyRay => RayInfo%AsyRay
CoreRay => RayInfo%CoreRay
RotRay => RayInfo%RotRay

!Geometry Info Pointing
Asy => CoreInfo%Asy
Pin => CoreInfo%Pin
PinInfo => CoreInfo%Pininfo
Cell => CoreInfo%CellInfo

!Tracking Dat Pointing
phia => TrackingDat%phianm
SrcAng => TrackingDat%SrcAngnm
xst => TrackingDat%xstnm
PhiAngIn => TrackingDat%PhiAngInnm
DcmpPhiAngIn => TrackingDat%DcmpPhiAngIn
DcmpPhiAngOut => TrackingDat%DcmpPhiAngOut
EXPA => TrackingDat%EXPA
EXPB => TrackingDat%EXPB

wtang => TrackingDat%wtang; AziMap => TrackingDat%AziMap
nAziAng = RayInfo%nAziAngle; nPolarAng = RayInfo%nPolarAngle

IF (lJout) THEN
  DO ipol = 1, nPolarAng
    DO iazi = 1, nAziAng
      wtang2(ipol, iazi, 1) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv)
      wtang2(ipol, iazi, 3) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv)
      wtang2(ipol, iazi, 2) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv)
      wtang2(ipol, iazi, 4) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv)
    ENDDO
  ENDDO
ENDIF

nAsyRay = DcmpAsyRay%nAsyRay
iRotRay = DcmpAsyRay%iRotRay
AsyRayList => DcmpAsyRay%AsyRayList
DirList => DcmpAsyRay%DirList
AziList => DcmpAsyRay%AziList
iAsy = DcmpAsyRay%iAsy; iRay = DcmpAsyRay%iRay

DO irot = 1, 2  
  PhiAngInSvIdx = RayInfo%PhiAngInSvIdx(iRotRay, irot)
  IF (DcmpAsyRay%lRotRayBeg(irot)) THEN
    PhiAngOut = PhiAngIn(1 : nPolarAng, gb : ge, PhiAnginSvIdx)
  ELSE
    PhiAngOut = DcmpPhiAngIn(1 : nPolarAng, gb : ge, irot, iRay, iAsy)
  ENDIF
  jbeg = 1; jend = nAsyRay; jinc = 1
  IF (irot .EQ. 2) THEN
    jbeg = nAsyRay; jend = 1; jinc = -1
  ENDIF
  DO j = jbeg, jend, jinc
    iAsyRay = AsyRayList(j)
    nPinRay = AsyRay(iAsyRay)%nCellRay
    idir = DirList(j); iazi = AziList(j)
    IF (irot .eq. 2) idir = mp(idir)
    AziSvIdx = AziMap(iAzi, :)
    DO ipol = 1, nPolarAng
      wt(ipol) = wtang(ipol, iazi)
      IF (lJout) wt2(ipol, 1 : 4) = wtang2(ipol, iazi, 1 : 4)
    ENDDO
    IF (idir .EQ. 1) THEN
      DO l = 1, nPinRay
        ipin = AsyRay(iAsyRay)%PinIdx(l)
        iceray = AsyRay(iAsyRay)%PinRayIdx(l)
        ipin = Asy(iAsy)%GlobalPinIdx(ipin)   
        icel = Pin(ipin)%Cell(iz)
        FsrIdxSt = Pin(ipin)%FsrIdxSt
        ibcel = Cell(icel)%basecellstr
        CellRay => Cell(ibcel)%CellRay(iceray)
        nRaySeg = CellRay%nSeg
        LocalFsrIdx => CellRay%LocalFsrIdx
        LenSeg => CellRay%LenSeg
        IF (lJout) THEN
          isurf = AsyRay(iAsyRay)%PinRaySurf(1, l)
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(1, ig, isurf, ipin) = Jout(1, ig, isurf, ipin) + wt(ipol) * PhiAngOut(ipol, ig)
              Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol, ig)
            ENDDO
          ENDDO
        ENDIF
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
            ENDDO
          ENDDO
        ENDDO
        IF (lJout) THEN
          isurf = AsyRay(iAsyRay)%PinRaySurf(2, l)
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(2, ig, isurf, ipin) = Jout(2, ig, isurf, ipin) + wt(ipol) * PhiAngOut(ipol, ig)
              Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol, ig)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ELSE
      DO l = nPinRay, 1, -1
        ipin = AsyRay(iAsyRay)%PinIdx(l)
        iceray = AsyRay(iAsyRay)%PinRayIdx(l)
        ipin = Asy(iAsy)%GlobalPinIdx(ipin)   
        icel = Pin(ipin)%Cell(iz)             
        FsrIdxSt = Pin(ipin)%FsrIdxSt
        ibcel = Cell(icel)%basecellstr
        CellRay => Cell(ibcel)%CellRay(iceray)
        nRaySeg = CellRay%nSeg
        LocalFsrIdx => CellRay%LocalFsrIdx
        LenSeg => CellRay%LenSeg
        IF (lJout) THEN
          isurf = AsyRay(iAsyRay)%PinRaySurf(2, l)
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(1, ig, isurf, ipin) = Jout(1, ig, isurf, ipin) + wt(ipol) * PhiAngOut(ipol, ig)
              Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol, ig)
            ENDDO
          ENDDO
        ENDIF
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
            ENDDO
          ENDDO
        ENDDO
        IF (lJout) THEN
          isurf = AsyRay(iAsyRay)%PinRaySurf(1, l)
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(2, ig, isurf, ipin) = Jout(2, ig, isurf, ipin) + wt(ipol) * PhiAngOut(ipol, ig)
              Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol, ig)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  DcmpPhiAngOut(1 : nPolarAng, gb : ge, irot, iRay, iAsy) = PhiAngOut
ENDDO

END SUBROUTINE TrackRotRayPn_Dcmp
! ------------------------------------------------------------------------------------------------------------