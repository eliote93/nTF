!--- CNJ Edit : CASMO Linear Source
SUBROUTINE GetFSRCentroid(RayInfo, CoreInfo, iz, lHybrid)
USE PARAM
USE TYPEDEF,    ONLY :  RayInfo_Type,       CoreInfo_Type,      Pin_Type,       Cell_Type
USE MOC_MOD,    ONLY :  TrackingDat,        nMaxRaySeg,         nMaxCoreRay
USE ALLOCS
USE PE_Mod,     ONLY :  PE
USE OMP_LIB
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
INTEGER :: iz
LOGICAL :: lHybrid

TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(Pin_Type), POINTER :: Pin(:)
REAL, POINTER :: CoreFsrCentroid(:, :, :), CoreFsrMxx(:, :), CoreFsrMyy(:, :), CoreFsrMxy(:, :)

INTEGER :: i, j, l, tid, FsrIdxSt, icel, ireg
INTEGER :: nFsr, nThread, nRotRay, nAziAng, nPolarAng, nxy, myzb, myze

nFsr = CoreInfo%nCoreFsr; nxy = CoreInfo%nxy
nRotRay = RayInfo%nRotRay
nAziAng = RayInfo%nAziAngle; nPolarAng = RayInfo%nPolarAngle
nThread = PE%nThread
myzb = PE%myzb; myze = PE%myze

IF (iz .EQ. myzb) THEN
  CALL Dmalloc0(CoreInfo%CoreFsrCentroid, 1, 2, 1, nFsr, myzb, myze)
  CALL Dmalloc0(CoreInfo%CoreFsrMxx, 1, nFsr, myzb, myze)
  CALL Dmalloc0(CoreInfo%CoreFsrMyy, 1, nFsr, myzb, myze)
  CALL Dmalloc0(CoreInfo%CoreFsrMxy, 1, nFsr, myzb, myze)
  DO i = 1, nThread
    CALL Dmalloc(TrackingDat(i)%FsrIdx, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(i)%dCentroid, 2, nPolarAng, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(i)%FsrCentroid, 2, nFsr)
    CALL Dmalloc(TrackingDat(i)%dMxx, nPolarAng, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(i)%FsrMxx, nFsr)
    CALL Dmalloc(TrackingDat(i)%dMyy, nPolarAng, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(i)%FsrMyy, nFsr)
    CALL Dmalloc(TrackingDat(i)%dMxy, nPolarAng, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(i)%FsrMxy, nFsr)
  ENDDO
ENDIF

CoreFsrCentroid => CoreInfo%CoreFsrCentroid
CoreFsrMxx => CoreInfo%CoreFsrMxx
CoreFsrMyy => CoreInfo%CoreFsrMyy
CoreFsrMxy => CoreInfo%CoreFsrMxy

CALL omp_set_dynamic(.FALSE.)
CALL omp_set_num_threads(nThread) 

!!!! Centroid Sweep !!!!

!$OMP PARALLEL PRIVATE(tid)
tid = OMP_GET_THREAD_NUM() + 1
!$OMP DO SCHEDULE(GUIDED)
DO i = 1, nRotRay
  CALL IntegrateRay(RayInfo, CoreInfo, TrackingDat(tid), i, iz, 0, lHybrid)
ENDDO
!$OMP END DO
!$OMP END PARALLEL

DO j = 1, nThread
  DO i = 1, nFsr
    CoreFsrCentroid(:, i, iz) = CoreFsrCentroid(:, i, iz) + TrackingDat(j)%FsrCentroid(:, i)
  ENDDO
ENDDO

Cell => CoreInfo%CellInfo; Pin => CoreInfo%Pin
DO l = 1, nxy
  FsrIdxSt = Pin(l)%FsrIdxSt; icel = Pin(l)%Cell(iz);
  DO j = 1, Cell(icel)%nFsr
    ireg = FsrIdxSt + j - 1
    CoreFsrCentroid(:, ireg, iz) = CoreFsrCentroid(:, ireg, iz) / Cell(icel)%vol(j)
  ENDDO
ENDDO

!!!! Moment Sweep !!!!

!$OMP PARALLEL PRIVATE(tid)
tid = OMP_GET_THREAD_NUM() + 1
!$OMP DO SCHEDULE(GUIDED)
DO i = 1, nRotRay
  CALL IntegrateRay(RayInfo, CoreInfo, TrackingDat(tid), i, iz, 1, lHybrid)
ENDDO
!$OMP END DO
!$OMP END PARALLEL

DO j = 1, nThread
  DO i = 1, nFsr
    CoreFsrMxx(i, iz) = CoreFsrMxx(i, iz) + TrackingDat(j)%FsrMxx(i)
    CoreFsrMyy(i, iz) = CoreFsrMyy(i, iz) + TrackingDat(j)%FsrMyy(i)
    CoreFsrMxy(i, iz) = CoreFsrMxy(i, iz) + TrackingDat(j)%FsrMxy(i)
  ENDDO
ENDDO

Cell => CoreInfo%CellInfo; Pin => CoreInfo%Pin
DO l = 1, nxy
  FsrIdxSt = Pin(l)%FsrIdxSt; icel = Pin(l)%Cell(iz);
  DO j = 1, Cell(icel)%nFsr
    ireg = FsrIdxSt + j - 1
    CoreFsrMxx(ireg, iz) = CoreFsrMxx(ireg, iz) / Cell(icel)%vol(j)
    CoreFsrMyy(ireg, iz) = CoreFsrMyy(ireg, iz) / Cell(icel)%vol(j)
    CoreFsrMxy(ireg, iz) = CoreFsrMxy(ireg, iz) / Cell(icel)%vol(j)
  ENDDO
ENDDO

!!!! Deallocate !!!!

IF (iz .EQ. myze) THEN
  DO i = 1, nThread
    DEALLOCATE(TrackingDat(i)%FsrIdx)
    DEALLOCATE(TrackingDat(i)%dCentroid, TrackingDat(i)%FsrCentroid)
    DEALLOCATE(TrackingDat(i)%dMxx, TrackingDat(i)%FsrMxx)
    DEALLOCATE(TrackingDat(i)%dMyy, TrackingDat(i)%FsrMyy)
    DEALLOCATE(TrackingDat(i)%dMxy, TrackingDat(i)%FsrMxy)
  ENDDO
ENDIF

END SUBROUTINE
    
SUBROUTINE IntegrateRay(RayInfo, CoreInfo, TrackingDat, irotray, iz, order, lHybrid)
USE PARAM
USE TYPEDEF,    ONLY :  RayInfo_Type,       CoreInfo_Type,      Pin_Type,           Cell_Type,              &
                        Asy_Type,           AsyInfo_Type,       AsyRayInfo_Type,    AziAngleInfo_Type,      &
                        PolarAngle_Type,    CoreRayInfo_Type,   RotRayInfo_Type,                            &
                        TrackingDat_Type,   CellRayInfo_Type
USE Moc_Mod,    ONLY :  nMaxRaySeg,         nMaxCoreRay
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
TYPE(TrackingDat_Type) :: TrackingDat
INTEGER :: irotray, iz, order
LOGICAL :: lHybrid

REAL, POINTER :: CoreFsrCentroid(:, :)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(AziAngleInfo_Type), POINTER :: AziAng(:)
TYPE(PolarAngle_Type), POINTER :: PolarAng(:)
TYPE(AsyRayInfo_type), POINTER :: AsyRay(:)
TYPE(CoreRayInfo_Type), POINTER :: CoreRay(:)
TYPE(RotRayInfo_Type), POINTER :: RotRay(:)
TYPE(CellRayInfo_Type), POINTER :: CellRay

REAL, POINTER :: LenSeg(:)
REAL, POINTER :: dCentroid(:, :, :, :), dMxx(:, :, :), dMyy(:, :, :), dMxy(:, :, :)
REAL, POINTER :: FsrCentroid(:, :), FsrMxx(:), FsrMyy(:), FsrMxy(:)
INTEGER, POINTER :: LocalFsrIdx(:), FsrIdx(:, :)

INTEGER :: j, k, l
INTEGER :: irsegidx, icoreray, iasyray, irayseg, iasy, itype, iceray, ipin, icel, ibcel, ireg, iazi, ipol, ir
INTEGER :: nAziAng, nPolarAng, nCoreRay, nAsyRay, nPinRay, nRaySeg, nTotRaySeg(nMaxCoreRay)
INTEGER :: FsrIdxSt
LOGICAL :: lFuel
REAL :: cmLenSeg, cmLenSeg1, cmLenSeg2, cmLenSeg3
REAL :: wttemp, wtang(10, 100)
REAL :: ax, ay, x0(2), y0(2), segCenter(2)

nAziAng = RayInfo%nAziAngle
nPolarAng = RayInfo%nPolarAngle

AziAng => RayInfo%AziAngle
PolarAng => RayInfo%PolarAngle
AsyRay => RayInfo%AsyRay
CoreRay => RayInfo%CoreRay
RotRay => RayInfo%RotRay

CoreFsrCentroid => CoreInfo%CoreFsrCentroid(:, :, iz)
Asy => CoreInfo%Asy
AsyInfo => CoreInfo%AsyInfo
Pin => CoreInfo%Pin
Cell => CoreInfo%CellInfo

FsrIdx => TrackingDat%FsrIdx
dCentroid => TrackingDat%dCentroid
FsrCentroid => TrackingDat%FsrCentroid
dMxx => TrackingDat%dMxx
FsrMxx => TrackingDat%FsrMxx
dMyy => TrackingDat%dMyy
FsrMyy => TrackingDat%FsrMyy
dMxy => TrackingDat%dMxy
FsrMxy => TrackingDat%FsrMxy

DO ipol = 1, nPolarAng
  wttemp = PolarAng(ipol)%weight * PolarAng(ipol)%sinv
  DO iazi = 1, nAziAng
    wtang(ipol, iazi) = wttemp * AziAng(iazi)%weight * AziAng(iazi)%del
  ENDDO
ENDDO

IF (order .EQ. 0) THEN
    
  nCoreRay = RotRay(iRotRay)%nRay
  DO j = 1, nCoreRay
    irsegidx = 0
    iCoreRay = RotRay(iRotRay)%RayIdx(j)
    nAsyRay = CoreRay(iCoreRay)%nRay
    iazi = CoreRay(iCoreRay)%iang
    DO k = 1, nAsyRay
      iasyray = CoreRay(iCoreRay)%AsyRayIdx(k)
      iasy = CoreRay(iCoreRay)%AsyIdx(k)
      IF (iasy .EQ. 0) CYCLE
      IF (lHybrid) THEN
        itype = Asy(iasy)%AsyType
        lFuel = AsyInfo(itype)%lFuel
        IF (lFuel) CYCLE
      ENDIF
      nPinRay = AsyRay(iAsyRay)%nCellRay
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
        DO iRaySeg = 1, nRaySeg
          ireg = FsrIdxSt + CellRay%LocalFsrIdx(iRaySeg) - 1 
          irsegidx = irsegidx + 1
          FsrIdx(irsegidx, j) = ireg
          cmLenSeg = LenSeg(iRaySeg) * 0.001_8
          segCenter(:) = CellRay%pts(:, iRaySeg) + CellRay%pts(:, iRaySeg + 1)
          DO ipol = 1, nPolarAng
            cmLenSeg1 = cmLenSeg / PolarAng(ipol)%sinv
            dCentroid(:, ipol, irsegidx, j) = segCenter(:) * cmLenSeg1 * wtang(ipol, iazi)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    nTotRaySeg(j) = irsegidx
  ENDDO

  DO j = 1, nCoreRay
    nRaySeg = nTotRaySeg(j)
    DO ir = 1, nRaySeg
      ireg = FsrIdx(ir, j)
      DO ipol = 1, nPolarAng
        FsrCentroid(:, ireg) = FsrCentroid(:, ireg) + dCentroid(:, ipol, ir, j)
      ENDDO
    ENDDO
  ENDDO

ELSEIF (order .EQ. 1) THEN
 
  nCoreRay = RotRay(iRotRay)%nRay
  DO j = 1, nCoreRay
    irsegidx = 0
    iCoreRay = RotRay(iRotRay)%RayIdx(j)
    nAsyRay = CoreRay(iCoreRay)%nRay
    iazi = CoreRay(iCoreRay)%iang
    DO k = 1, nAsyRay
      iasyray = CoreRay(iCoreRay)%AsyRayIdx(k)
      iasy = CoreRay(iCoreRay)%AsyIdx(k)
      IF (iasy .EQ. 0) CYCLE
      IF (lHybrid) THEN
        itype = Asy(iasy)%AsyType
        lFuel = AsyInfo(itype)%lFuel
        IF (lFuel) CYCLE
      ENDIF
      nPinRay = AsyRay(iAsyRay)%nCellRay
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
        DO iRaySeg = 1, nRaySeg
          ireg = FsrIdxSt + CellRay%LocalFsrIdx(iRaySeg) - 1 
          irsegidx = irsegidx + 1
          FsrIdx(irsegidx, j) = ireg
          cmLenSeg = LenSeg(iRaySeg) * 0.001_8
          x0(1) = CellRay%pts(1, iRaySeg) + CellRay%pts(1, iRaySeg + 1) - cmLenSeg * AziAng(iazi)%cosv
          x0(2) = CellRay%pts(1, iRaySeg) + CellRay%pts(1, iRaySeg + 1) + cmLenSeg * AziAng(iazi)%cosv
          y0(1) = CellRay%pts(2, iRaySeg) + CellRay%pts(2, iRaySeg + 1) - cmLenSeg * AziAng(iazi)%sinv
          y0(2) = CellRay%pts(2, iRaySeg) + CellRay%pts(2, iRaySeg + 1) + cmLenSeg * AziAng(iazi)%sinv
          x0(:) = x0(:) * half; y0(:) = y0(:) * half
          x0(:) = x0(:) - CoreFsrCentroid(1, ireg); y0(:) = y0(:) - CoreFsrCentroid(2, ireg)
          DO ipol = 1, nPolarAng
            ax = AziAng(iazi)%cosv * PolarAng(ipol)%sinv; ay = AziAng(iazi)%sinv * PolarAng(ipol)%sinv
            cmLenSeg1 = cmLenSeg / PolarAng(ipol)%sinv; cmLenSeg2 = half * cmLenSeg1 ** 2; cmLenSeg3 = rthree * cmLenSeg1 ** 3
            dMxx(ipol, irsegidx, j) = (cmLenSeg1 * (x0(1) ** 2 + x0(2) ** 2)                                                 &
                                      + 2._8 * cmLenSeg2 * (x0(1) - x0(2)) * ax                                              &
                                      + 2._8 * cmLenSeg3 * ax ** 2) * wtang(ipol, iazi)
            dMyy(ipol, irsegidx, j) = (cmLenSeg1 * (y0(1) ** 2 + y0(2) ** 2)                                                 &
                                      + 2._8 * cmLenSeg2 * (y0(1) - y0(2)) * ay                                              &
                                      + 2._8 * cmLenSeg3 * ay ** 2) * wtang(ipol, iazi)
            dMxy(ipol, irsegidx, j) = (cmLenSeg1 * (x0(1) * y0(1) + x0(2) * y0(2))                                           &
                                      + cmLenSeg2 * ((x0(1) - x0(2)) * ay + (y0(1) - y0(2)) * ax)                            &
                                      + 2._8 * cmLenSeg3 * ax * ay) * wtang(ipol, iazi)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    nTotRaySeg(j) = irsegidx
  ENDDO

  DO j = 1, nCoreRay
    nRaySeg = nTotRaySeg(j)
    DO ir = 1, nRaySeg
      ireg = FsrIdx(ir, j)
      DO ipol = 1, nPolarAng
        FsrMxx(ireg) = FsrMxx(ireg) + dMxx(ipol, ir, j)
        FsrMyy(ireg) = FsrMyy(ireg) + dMyy(ipol, ir, j)
        FsrMxy(ireg) = FsrMxy(ireg) + dMxy(ipol, ir, j)
      ENDDO
    ENDDO
  ENDDO

ENDIF

END SUBROUTINE