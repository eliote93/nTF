#include <defines.h>
SUBROUTINE RayTrace_OMP(RayInfo, CoreInfo, phis, PhiAngIn, xst, src, jout, iz, ljout, FastMocLv, lAFSS)
USE PARAM
USE TYPEDEF, ONLY :   RayInfo_Type,      CoreInfo_type,   Pin_Type,         Cell_Type,        &
                      AziAngleInfo_Type, PolarAngle_Type
USE Moc_Mod, ONLY :   nMaxRaySeg,        nMaxCellRay,     nMaxAsyRay,       nMaxCoreRay,      &
                      Expa,              Expb,            ApproxExp,        TrackingDat,      &
                      wtang
USE geom,    ONLY : nbd
USE TIMER
USE BasicOperation, ONLY : CP_CA, CP_VA
USE ALLOCS
USE PE_MOD,  ONLY :   PE
USE OMP_LIB
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phis(:), PhiAngIn(:, :), xst(:), src(:), jout(:, :, :)
INTEGER :: iz
LOGICAL :: ljout
INTEGER, OPTIONAL :: FastMocLv
LOGICAL, OPTIONAL :: lAFSS

LOGICAL, SAVE :: lfirst
DATA lfirst /.TRUE./

TYPE(AziAngleInfo_Type), POINTER :: AziAng(:)
TYPE(PolarAngle_Type), POINTER :: PolarAng(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(Pin_Type), POINTER :: Pin(:)

INTEGER :: nAziAng, nPolarAng, nPhiAngSv
INTEGER :: nRotray
INTEGER :: nFsr, nxy
INTEGER :: nThread

INTEGER :: tid        !Thread Id
INTEGER :: FsrIdxSt, icel, ireg, iazi, ipol
INTEGER :: i, j, l
REAL :: wttemp

LOGICAL :: lJout_OMP(100)

!Get Essential Variable
AziAng => RayInfo%AziAngle; PolarAng => RayInfo%PolarAngle;
nAziAng = RayInfo%nAziAngle; nPolarAng = RayInfo%nPolarAngle
nPhiAngSv = RayInfo%nPhiAngSv
nRotRay = RayInfo%nRotRay
nFsr = CoreInfo%nCoreFsr; nxy = CoreInfo%nxy
nThread = PE%nThread

!Allocate Static Memery
IF(lfirst) THEN
  lFirst = FALSE
  CALL ApproxExp(RayInfo%PolarAngle, nPolarAng)
  DO tid = 1, nThread
    IF(TrackingDat(tid)%lAlloc) CYCLE
    CALL Dmalloc(TrackingDat(tid)%FsrIdx, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(tid)%ExpAppIdx, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(tid)%OptLenList, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(tid)%ExpAppPolar, nPolarAng, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(tid)%PhiAngOutPolar, nPolarAng, nMaxRaySeg + 2)
    CALL Dmalloc(TrackingDat(tid)%Phis, nFsr)
    !CALL Dmalloc(TrackingDat(tid)%Jout, 2, nbd, nxy)
    CALL Dmalloc(TrackingDat(tid)%Jout, 3, nbd, nxy) !---BYS edit / 150612 Surface flux (1:in 2:out 3:surfphi)
    TrackingDat(tid)%lAlloc = .TRUE.
  ENDDO
  ALLOCATE(wtang(nPolarAng, nAziAng))
  DO ipol = 1, nPolarAng
    wttemp = PolarAng(ipol)%weight * PolarAng(ipol)%sinv
    DO iazi = 1, nAziAng
      wtang(ipol, iazi) = wttemp * AziAng(iazi)%weight * AziAng(iazi)%del
    ENDDO
  ENDDO
ENDIF

!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(nThread)

DO tid = 1, nThread
  TrackingDat(tid)%Expa => Expa; TrackingDat(tid)%Expb => Expb
  TrackingDat(tid)%PhiAngIn => PhiAngIn(:, :)
  TrackingDat(tid)%src => src; TrackingDat(tid)%xst => xst
  TrackingDat(tid)%wtang => wtang
ENDDO

!Pointing
!Flux and current set to zero
CALL CP_CA(phis, ZERO, nFsr)
!IF(ljout)  CALL CP_CA(Jout, ZERO, 2, 4, CoreInfo%nxy) !---BYS edit / 150612 Surface flux
IF(ljout)  CALL CP_CA(Jout, ZERO, 3, nbd, CoreInfo%nxy) !---BYS edit / 150612 Surface flux
lJout_OMP = lJout
!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(tid, i, j)
!$  tid = omp_get_thread_num()+1
!$OMP MASTER
DO j = 1, nThread
  DO i = 1, nFsr
    TrackingDat (j)%phis(i) = zero
  ENDDO
ENDDO
DO j = 1, nThread
  DO i = 1, nxy
     TrackingDat(j)%jout(:, :, i) = zero
     !TrackingDat(j)%Phisurf(:, i) = zero !---BYS edit / 150612 Surface flux
  ENDDO
ENDDO
!$OMP END MASTER

!$OMP BARRIER

!$OMP DO SCHEDULE(GUIDED)
DO i = 1, nRotRay
  CALL TrackRotRay(RayInfo, CoreInfo, TrackingDat(tid), lJout_OMP(tid), i, iz, FastMocLv)
END DO
!$OMP END DO

!$OMP END PARALLEL
DO j = 1, nThread
  DO i = 1, nFsr
    phis(i) = phis(i) + TrackingDat(j)%Phis(i)
  ENDDO
ENDDO
  IF(ljout) THEN
    DO j = 1, nThread
      DO i = 1, nxy
          jout(:, :, i) = jout(:, :, i) + TrackingDat(j)%jout(:, :, i)
          !phisurf(:, i) = phisurf(:, i) + TrackingDat(j)%phisurf(:, i)    !---BYS edit / 150612 Surface flux
      ENDDO
    ENDDO
  ENDIF

Cell => CoreInfo%CellInfo; Pin => CoreInfo%Pin
!$OMP PARALLEL DEFAULT(SHARED)  &
!$OMP PRIVATE(l, j, FsrIdxSt, icel, ireg)
!$OMP DO
DO l = 1, nxy
  FsrIdxSt = Pin(l)%FsrIdxSt; icel = Pin(l)%Cell(iz);
  DO j = 1, Cell(icel)%nFsr
    ireg = FsrIdxSt + j - 1
    phis(ireg) = phis(ireg)/xst(ireg)/Cell(icel)%vol(j) + src(ireg)
  ENDDO
ENDDO
!$OMP END DO

!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE RayTrace_OMP_AFSS(RayInfo, CoreInfo, phis, PhiAngIn, xst, src, jout, iz, ljout, FastMocLv, lAFSS)
USE PARAM
USE TYPEDEF, ONLY :   RayInfo_Type,      coreinfo_type,                                       &
                      TrackingDat_Type,                                                       &
                      Pin_Type,          Cell_Type,                                           &
                      AziAngleInfo_Type, PolarAngle_Type
USE Moc_Mod, ONLY :   nMaxRaySeg,        nMaxCellRay,     nMaxAsyRay,       nMaxCoreRay,      &
                      Expa,              Expb,            TrackingDat,                        &
                      ApproxExp,         wtang,           Phia1g,           Phia2g
USE geom,           ONLY : nbd
USE cntl,           ONLY : nTracerCntl
USE TIMER
USE BasicOperation, ONLY : CP_CA, CP_VA
USE ALLOCS
USE PE_MOD,  ONLY :   PE
USE OMP_LIB
IMPLICIT NONE

!INTEGER, PARAMETER :: nThread = 4
TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phis(:), PhiAngIn(:, :), xst(:), src(:), jout(:, :, :)
INTEGER :: iz
LOGICAL :: ljout
INTEGER, OPTIONAL :: FastMocLv
LOGICAL, OPTIONAL :: lAFSS

LOGICAL, SAVE :: lfirst
DATA lfirst /.TRUE./

!Pointing Variable
TYPE(AziAngleInfo_Type), POINTER :: AziAng(:)
TYPE(PolarAngle_Type), POINTER :: PolarAng(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(Pin_Type), POINTER :: Pin(:)

INTEGER :: nAziAng, nPolarAng, nPhiAngSv
INTEGER :: nRotray
INTEGER :: nFsr, nxy
INTEGER :: nThread

INTEGER :: iRotRay
INTEGER :: tid        !Thread Id
INTEGER :: FsrIdxSt, icel, ireg, iazi, ipol
INTEGER :: i, j, k, l, m
REAL :: wttemp

LOGICAL :: lJout_OMP(100)

LOGICAL, SAVE :: lScatBd
INTEGER, POINTER, SAVE :: nAziRotRay(:), nAziRotRay0(:), nSegRotRay(:), OmpRayBeg(:), OmpRayEnd(:)
INTEGER, POINTER, SAVE :: OmpRayBegBd(:,:), OmpRayEndBd(:,:)
INTEGER, SAVE :: OmpTemp, OmpAng
INTEGER, SAVE :: nOmpAng
INTEGER, POINTER, SAVE :: OmpMap(:,:)
INTEGER, POINTER, SAVE :: OmpRayList(:)
INTEGER, POINTER, SAVE :: OmpRayNum(:,:)

!Get Essential Variable
AziAng => RayInfo%AziAngle; PolarAng => RayInfo%PolarAngle;
nAziAng = RayInfo%nAziAngle; nPolarAng = RayInfo%nPolarAngle
nPhiAngSv = RayInfo%nPhiAngSv
nRotRay = RayInfo%nRotRay
nFsr = CoreInfo%nCoreFsr; nxy = CoreInfo%nxy
nThread = PE%nThread

IF(lfirst) THEN
lScatBd = nTracerCntl%lScatBd
IF(lScatBd) THEN
  ALLOCATE(nAziRotRay(1:nAziAng/2)); ALLOCATE(nAziRotRay0(1:nThread)); ALLOCATE(nSegRotRay(1:nAziAng/2))
  nAziRotRay(:) = 0; nSegRotRay(:) = 0
  DO i = 1, nRotRay
    j = RayInfo%RotRay(i)%Ang1; nAziRotRay(j) = nAziRotRay(j) + 1; nSegRotRay(j) = nSegRotRay(j) + RayInfo%RotRay(i)%nSeg
  ENDDO
  ALLOCATE(OmpRayBeg(0:nThread), OmpRayEnd(0:nThread))
  k = nAziAng/2/nThread; nAziRotRay0(:) = 0; OmpRayBeg(0) = 1; OmpTemp = 0; OmpRayEnd(0) = 0; l = 0
  DO i = 1, nThread
    OmpRayBeg(i) = OmpRayBeg(i-1) + OmpTemp
    DO j = 1, k
      l = l + 1; nAziRotRay0(i) = nAziRotRay0(i) + nAziRotRay(l)
    ENDDO
    OmpTemp = nAziRotRay0(i); OmpRayEnd(i) = OmpRayEnd(i-1) + OmpTemp
  ENDDO
  !Parallel Ang
  DO i = 1, nRotRay
    RayInfo%RotRay(i)%OmpAng1 = MOD(RayInfo%RotRay(i)%Ang1, k)
    IF(RayInfo%RotRay(i)%OmpAng1 .EQ. 0) RayInfo%RotRay(i)%OmpAng1 = k
    RayInfo%RotRay(i)%OmpAng2 = RayInfo%RotRay(i)%OmpAng1 + k
  ENDDO
  DO i = 1, nRotRay
    DO j = 1, RayInfo%RotRay(i)%nRay, 2
      RayInfo%RotRay(i)%OmpRayIdx(j) = RayInfo%RotRay(i)%OmpAng1; RayInfo%RotRay(i)%OmpRayIdx(j+1) = RayInfo%RotRay(i)%OmpAng2
    ENDDO
  ENDDO
  ALLOCATE(OmpMap(nThread, k*2)); OmpTemp = 0
  DO i = 1, nThread
    DO j = 1, k
      OmpTemp = OmpTemp + 1; OmpMap(i,j) = OmpTemp; OmpMap(i,j+k) = nAziAng - OmpTemp + 1
    ENDDO
  ENDDO
  nOmpAng = nAziAng/nThread
ELSE
  ALLOCATE(nAziRotRay(1:nAziAng)); ALLOCATE(nAziRotRay0(1:2*nThread)); ALLOCATE(nSegRotRay(1:nAziAng))
  nAziRotRay(:) = 0; nSegRotRay(:) = 0
  DO i = 1, nRotRay
    j = RayInfo%RotRay(i)%Ang1; nAziRotRay(j) = nAziRotRay(j) + 1; nSegRotRay(j) = nSegRotRay(j) + RayInfo%RotRay(i)%nSeg
  ENDDO
  ALLOCATE(OmpRayBeg(0:2*nThread), OmpRayEnd(0:2*nThread)); ALLOCATE(OmpRayBegBd(2, 0:nThread), OmpRayEndBd(2, 0:nThread))
  k = nAziAng/nThread/2; nAziRotRay0(:) = 0; OmpRayBegBd(1,0) = 1; OmpTemp = 0; OmpRayEndBd(1,0) = 0; l = 0
  DO i = 1, nThread
    OmpRayBegBd(1,i) = OmpRayBegBd(1,i-1) + OmpTemp
    DO j = 1, k
      l = l + 1; nAziRotRay0(i) = nAziRotRay0(i) + nAziRotRay(l)
    ENDDO
    OmpTemp = nAziRotRay0(i); OmpRayEndBd(1,i) = OmpRayEndBd(1,i-1) + OmpTemp
  ENDDO
  nAziRotRay0(:) = 0; l = nAziAng + 1; OmpRayBegBd(2,0) = nRotRay + 1; OmpTemp = 0; OmpRayEndBd(2,0) = nRotRay;
  DO i = 1, nThread
    OmpRayEndBd(2,i) = OmpRayEndBd(2,i-1) - OmpTemp
    DO j = 1, k
      l = l - 1; nAziRotRay0(i) = nAziRotRay0(i) + nAziRotRay(l)
    ENDDO
    OmpTemp = nAziRotRay0(i); OmpRayBegBd(2,i) = OmpRayBegBd(2,i-1) - OmpTemp
  ENDDO
  !Parallel Ang
  DO i = 1, nRotRay
    IF(RayInfo%RotRay(i)%nRay .GT. 1) THEN
      IF(RayInfo%RotRay(i)%Ang1 .LE. nAziAng/2) THEN
        RayInfo%RotRay(i)%OmpAng1 = MOD(RayInfo%RotRay(i)%Ang1, k)
        IF(RayInfo%RotRay(i)%OmpAng1 .EQ. 0) RayInfo%RotRay(i)%OmpAng1 = k
        RayInfo%RotRay(i)%OmpAng2 = RayInfo%RotRay(i)%OmpAng1 + k
      ELSEIF(RayInfo%RotRay(i)%Ang1 .GT. nAziAng/2) THEN
        RayInfo%RotRay(i)%OmpAng1 = MOD(nAziAng+1-RayInfo%RotRay(i)%Ang1, k) + k
        IF(RayInfo%RotRay(i)%OmpAng1 .EQ. k) RayInfo%RotRay(i)%OmpAng1 = RayInfo%RotRay(i)%OmpAng1 + k
        RayInfo%RotRay(i)%OmpAng2 = RayInfo%RotRay(i)%OmpAng1 - k
      ENDIF
    ELSEIF(RayInfo%RotRay(i)%nRay .EQ. 1) THEN
      IF(RayInfo%RotRay(i)%Ang1 .LE. nAziAng/2) THEN
        RayInfo%RotRay(i)%OmpAng1 = MOD(RayInfo%RotRay(i)%Ang1, k)
        IF(RayInfo%RotRay(i)%OmpAng1 .EQ. 0) RayInfo%RotRay(i)%OmpAng1 = k
        RayInfo%RotRay(i)%OmpAng2 = 0
      ELSEIF(RayInfo%RotRay(i)%Ang1 .GT. nAziAng/2) THEN
        RayInfo%RotRay(i)%OmpAng1 = MOD(nAziAng+1-RayInfo%RotRay(i)%Ang1, k) + k
        IF(RayInfo%RotRay(i)%OmpAng1 .EQ. k) RayInfo%RotRay(i)%OmpAng1 = RayInfo%RotRay(i)%OmpAng1 + k
        RayInfo%RotRay(i)%OmpAng2 = 0
      ENDIF
    ENDIF
  ENDDO
  DO i = 1, nRotRay
    IF(RayInfo%RotRay(i)%nRay .GT. 1) THEN
      IF(MOD(RayInfo%RotRay(i)%nRay, 2) .EQ. 0) THEN
        DO j = 1, RayInfo%RotRay(i)%nRay, 2
          RayInfo%RotRay(i)%OmpRayIdx(j) = RayInfo%RotRay(i)%OmpAng1; RayInfo%RotRay(i)%OmpRayIdx(j+1) = RayInfo%RotRay(i)%OmpAng2
        ENDDO
      ELSE
        DO j = 1, RayInfo%RotRay(i)%nRay-1, 2
          RayInfo%RotRay(i)%OmpRayIdx(j) = RayInfo%RotRay(i)%OmpAng1; RayInfo%RotRay(i)%OmpRayIdx(j+1) = RayInfo%RotRay(i)%OmpAng2
        ENDDO
        RayInfo%RotRay(i)%OmpRayIdx(RayInfo%RotRay(i)%nRay) = RayInfo%RotRay(i)%OmpAng1
      ENDIF
    ELSE
      RayInfo%RotRay(i)%OmpRayIdx(1) = RayInfo%RotRay(i)%OmpAng1
    ENDIF
  ENDDO
  ALLOCATE(OmpMap(nThread, k*2)); OmpTemp = 0
  DO i = 1, nThread
    DO j = 1, k
      OmpTemp = OmpTemp + 1; OmpMap(i,j) = OmpTemp; OmpMap(i,j+k) = nAziAng - OmpTemp + 1
    ENDDO
  ENDDO
  ALLOCATE(OmpRayList(nThread)); ALLOCATE(OmpRayNum(nThread, nRotRay))
  nOmpAng = nAziAng/nThread; OmpRayList(:) = 0; OmpRayNum(:,:) = 0
  DO i = 1, nThread
    DO j = 1, 2
      OmpRayList(i) = OmpRayList(i) + OmpRayEndBd(j,i) - OmpRayBegBd(j,i) + 1
    ENDDO
  ENDDO
  DO i = 1, nThread
    m = 0
    DO j = 1, 2
      l = 0
      DO k = 1, OmpRayEndBd(j,i) - OmpRayBegBd(j,i) + 1
        m = m + 1; OmpRayNum(i,m) = OmpRayBegBd(j,i) + l; l = l + 1
      ENDDO
    ENDDO
  ENDDO
ENDIF
IF(MOD(nAziANg/2, nThread) .NE. 0) THEN
  PRINT*, 'WRONG_MOC_TRD'
  STOP
ENDIF
ENDIF

!Allocate Static Memery
IF(lfirst) THEN
  lFirst = FALSE
  CALL ApproxExp(RayInfo%PolarAngle, nPolarAng)
  DO tid = 1, nThread
    IF(TrackingDat(tid)%lAlloc) CYCLE
    CALL Dmalloc(TrackingDat(tid)%FsrIdx, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(tid)%ExpAppIdx, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(tid)%OptLenList, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(tid)%ExpAppPolar, nPolarAng, nMaxRaySeg, nMaxCoreRay)
    CALL Dmalloc(TrackingDat(tid)%PhiAngOutPolar, nPolarAng, nMaxRaySeg+2)
    CALL Dmalloc(TrackingDat(tid)%phi1a, nPolarAng, nFsr, nOmpAng)
    CALL Dmalloc(TrackingDat(tid)%phi2a, nPolarAng, nFsr, nOmpAng)
    CALL Dmalloc(TrackingDat(tid)%Jout, 3, nbd, nxy)
    TrackingDat(tid)%Expa => Expa; TrackingDat(tid)%Expb => Expb
    TrackingDat(tid)%lAlloc = .TRUE.
  ENDDO
  ALLOCATE(wtang(nPolarAng, nAziAng))
  DO ipol = 1, nPolarAng
    wttemp = PolarAng(ipol)%weight * PolarAng(ipol)%sinv
    DO iazi = 1, nAziAng
      wtang(ipol, iazi) = wttemp  * AziAng(iazi)%weight * AziAng(iazi)%del
    ENDDO
  ENDDO
ENDIF

!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(nThread)

DO tid = 1, nThread
  TrackingDat(tid)%PhiAngIn => PhiAngIn(:, :)
  TrackingDat(tid)%src => src; TrackingDat(tid)%xst => xst
  TrackingDat(tid)%wtang => WtAng;
ENDDO

!Pointing
!Flux and current set to zero
CALL CP_CA(phis, ZERO, nFsr)
IF(ljout) CALL CP_CA(Jout, ZERO, 3, nbd, CoreInfo%nxy)
lJout_OMP = lJout
!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(irotray, tid, i, j, k, l, iazi, ipol, OmpAng)
tid = 1
!$  tid = omp_get_thread_num()+1
j = tid

DO iazi = 1, nOmpAng
  TrackingDat(j)%phi1a(:,:,iazi) = zero
  TrackingDat(j)%phi2a(:,:,iazi) = zero
ENDDO

DO i = 1, nxy
  TrackingDat(j)%jout(:, :, i) = zero
ENDDO

!$OMP BARRIER

IF(lScatBd) THEN
  DO i = OmpRayBeg(tid), OmpRayEnd(tid)
    irotray = i
    CALL TrackRotRay_AFSS(RayInfo, CoreInfo, TrackingDat(tid), lJout_OMP(tid), irotray, iz, FastMocLv)
  ENDDO
ELSE
  DO i = 1, 2
    DO j = OmpRayBegBd(i, tid), OmpRayEndBd(i, tid)
      irotray = j
      CALL TrackRotRay_AFSS(RayInfo, CoreInfo, TrackingDat(tid), lJout_OMP(tid), irotray, iz, FastMocLv)
    ENDDO
  ENDDO
ENDIF

!$OMP BARRIER

DO j = 1, nThread
  DO iazi = 1, nOmpAng
   DO i = PE%myOmpFsrBeg(tid), PE%myOmpFsrEnd(tid)
      DO ipol = 1, nPolarAng
        OmpAng = OmpMap(j, iazi)
        phis(i) = phis(i) + wtang(ipol, OmpAng) * (TrackingDat(j)%phi1a(ipol, i, iazi) + TrackingDat(j)%phi2a(ipol, i, iazi))
      ENDDO
    ENDDO
  ENDDO
ENDDO

IF(ljout) THEN
  DO j = 1, nThread
    DO i = PE%myOmpNxyBeg(tid), PE%myOmpNxyEnd(tid)
      jout(:, :, i) = jout(:, :, i) + TrackingDat(j)%jout(:, :, i)
    ENDDO
  ENDDO
ENDIF

!$OMP END PARALLEL

Cell => CoreInfo%CellInfo; Pin => CoreInfo%Pin

!$OMP PARALLEL DEFAULT(SHARED)  &
!$OMP PRIVATE(l, j, FsrIdxSt, icel, ireg)
!$OMP DO
DO l = 1, nxy
  FsrIdxSt = Pin(l)%FsrIdxSt; icel = Pin(l)%Cell(iz);
  DO j = 1, Cell(icel)%nFsr
    ireg = FsrIdxSt + j - 1
    phis(ireg) = phis(ireg)/xst(ireg)/Cell(icel)%vol(j) + src(ireg)
  ENDDO
ENDDO
!$OMP END DO

!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE TrackRotRay(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, FastMocLv)
USE PARAM
USE TYPEDEF,  ONLY :  RayInfo_Type,      coreinfo_type,                                       &
                      Pin_Type,          Asy_Type,        AsyInfo_Type,     PinInfo_Type,     &
                      Cell_Type,                                                              &
                      AziAngleInfo_Type, PolarAngle_Type, ModRayInfo_type,  AsyRayInfo_type,  &
                      CoreRayInfo_Type,  RotRayInfo_Type, CellRayInfo_type, FastCoreRayDat_Type,&
                      TrackingDat_Type,  FastRaySegDat_Type
USE Moc_Mod, ONLY :   nMaxRaySeg,        nMaxCellRay,     nMaxAsyRay,       nMaxCoreRay
USE BasicOperation, ONLY : CP_CA, CP_VA
USE CNTL, ONLY: nTracerCntl
IMPLICIT NONE

TYPE(RayInfo_Type), INTENT(INOUT) :: RayInfo
TYPE(CoreInfo_Type), INTENT(INOUT) :: CoreInfo
TYPE(TrackingDat_Type), INTENT(INOUT) :: TrackingDat
LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz
INTEGER, INTENT(IN) :: FastMocLv

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Asy_Type), POINTER :: Asy(:)
!TYPE(AsyInfo_Type), POINTER :: AsyInfo
TYPE(PinInfo_Type), POINTER :: PinInfo(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(AziAngleInfo_Type), POINTER :: AziAng(:)
TYPE(PolarAngle_Type), POINTER :: PolarAng(:)
!TYPE(ModRayInfo_type), POINTER :: ModRay
TYPE(AsyRayInfo_type), POINTER :: AsyRay(:)
TYPE(CoreRayInfo_Type), POINTER :: CoreRay(:)
TYPE(RotRayInfo_Type), POINTER :: RotRay(:)
TYPE(FastCoreRayDat_Type), POINTER :: FastRay
TYPE(CellRayInfo_Type), POINTER :: CellRay
TYPE(CellRayInfo_Type), POINTER :: CellRay1D
TYPE(FastRaySegDat_Type), POINTER :: FastRaySeg

REAL, POINTER :: LenSeg(:)
INTEGER, POINTER :: LocalFsrIdx(:)
INTEGER, POINTER :: FsrIdx(:, :),  ExpAppIdx(:, :)
REAL, POINTER :: OptLenList(:, :), ExpAppPolar(:,:,:)
REAL, POINTER :: phis(:), src(:), xst(:), jout(:, :, :)
REAL, POINTER :: PhiAngOutPolar(:, :), PhiAngIn(:, :)
REAL, POINTER :: EXPA(:, :), EXPB(:, :)
REAL, POINTER :: wtang(:, :)

INTEGER :: mp(2)
INTEGER :: iazi, ipol, iCoreRay, iasyray, iceray, irayseg, irot, itype, idir     !Ray related index
INTEGER :: nRotRay, nCoreRay, nAsyRay, nPinRay, nRaySeg, FsrIdxSt
INTEGER :: nPolarAng, nAziAng, nPhiAngSv
INTEGER :: ipin, icel, iasy, ireg, isurf                                                  !Geometries index
INTEGER :: irsegidx, icellrayidx, PhiAnginSvIdx, PhiAngOutSvIdx
INTEGER :: nFsr, nxy
INTEGER :: i, j, k, l, m, jbeg, jend, jinc, ir, ir1
INTEGER :: ibcel

INTEGER :: CellRayIdxSt(nMaxCellRay, nMaxCoreRay, 2)
INTEGER :: PinIdx(nMaxCellRay, nMaxCoreRay)
INTEGER :: SurfIdx(nMaxCellRay, nMaxCoreRay, 2)
INTEGER :: nTotRaySeg(nMaxCoreRay)
INTEGER :: nTotCellRay(nMaxCoreRay)

REAL :: wtang2(100, 100, 4)
REAL :: tau, phiobd, phid, wt, wt2(4)
REAL, ALLOCATABLE :: phiobdPolar(:)

LOGICAL :: lFast

DATA mp /2, 1/

IF (nTracerCntl%lHex) THEN
  CALL HexTrackRotRay(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, FastMocLv)
  RETURN
END IF

lFast = FALSE
IF(FastMocLv .GT. 0) lFast = .TRUE.

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

!Tracking Dat Pointing
FsrIdx => TrackingDat%FsrIdx
ExpAppIdx => TrackingDat%ExpAppIdx
OptLenList => TrackingDat%OptLenList
ExpAppPolar => TrackingDat%ExpAppPolar
Phis => TrackingDat%phis; src => TrackingDat%src
xst => TrackingDat%xst; jout => TrackingDat%jout
PhiAngOutPolar => TrackingDat%PhiAngOutPolar
PhiAngIn => TrackingDat%phiAngIn
EXPA => TrackingDat%EXPA
EXPB => TrackingDat%EXPB
wtang => TrackingDat%wtang

nAziAng = RayInfo%nAziAngle; nPolarAng = RayInfo%nPolarAngle
nPhiAngSv = RayInfo%nPhiAngSv
nFsr = CoreInfo%nCoreFsr; nxy = CoreInfo%nxy

IF(lJout)THEN
    DO ipol = 1, nPolarAng
        DO iazi = 1, nAziAng !(SWNE) ????
            wtang2(ipol, iazi, 1) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv)
            wtang2(ipol, iazi, 3) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv)
            wtang2(ipol, iazi, 2) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv)
            wtang2(ipol, iazi, 4) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv)
        ENDDO
    ENDDO
ENDIF

i = iRotRay
nCoreRay = RotRay(irotRay)%nRay
IF(.NOT. lFast) THEN
  DO j = 1, nCoreRay    !Core Ray Sweep
    irsegidx = 0;   icellrayidx = 0
    iCoreRay = RotRay(iRotRay)%RayIdx(j)
    nAsyRay = CoreRay(iCoreRay)%nRay
    DO k = 1, nAsyRay    !Assembly Ray Sweep
      iasyray = CoreRay(iCoreRay)%AsyRayIdx(k)
      iasy = CoreRay(iCoreRay)%AsyIdx(k)
      IF(iasy .EQ. 0)  CYCLE   !Skip Dummy Assembly
      nPinRay = AsyRay(iAsyRay)%nCellRay
      itype = Asy(iasy)%PartialAsyFlag
      DO l = 1, nPinRay   !Pin Ray Sweep
        ipin = AsyRay(iAsyRay)%PinIdx(l)      !Local Pin Idx(within Assembly)
        iceray = AsyRay(iAsyRay)%PinRayIdx(l) !Cell Ray Index
        ipin = Asy(iAsy)%GlobalPinIdx(ipin)   !Global Pin Index
        icel = Pin(ipin)%Cell(iz)             !Cell Type
        FsrIdxSt = Pin(ipin)%FsrIdxSt
        ibcel=Cell(icel)%basecellstr
        CellRay => Cell(ibcel)%CellRay(iceray) !Pointing Cell Ray
        icellrayidx = icellrayidx + 1
        PinIdx(icellrayidx, j) = ipin
        CellRayIdxSt(icellrayidx, j, 2) = irsegidx + 1
        nRaySeg = CellRay%nSeg
        LocalFsrIdx => CellRay%LocalFsrIdx
        LenSeg => CellRay%LenSeg
        DO iRaySeg = 1, nRaySeg
          !ireg = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
          ireg = FsrIdxSt + CellRay%LocalFsrIdx(iRaySeg) - 1
          tau = - LenSeg(iRaySeg) * xst(ireg)   !

             !Global Fsr Index
          !tau = -1000._8 * CellRay%LenSeg(iRaySeg) * xst(ireg)   !
          tau = - CellRay%LenSeg(iRaySeg) * xst(ireg)   !
          irsegidx = irsegidx + 1
          FsrIdx(irsegidx, j) = ireg
          OptLenList(irsegidx, j) = tau
          ExpAppIdx(irsegidx, j) = max(INT(tau), -40000)
          ExpAppIdx(irsegidx, j) = min(0, ExpAppIdx(irsegidx, j))
          !OptLenList(irsegidx, j) = OptLenList(irsegidx, j)/1000._8
          continue
        ENDDO   !End of Ray Segments Sweep, irayseg
        CellRayIdxSt(icellrayidx, j, 1) = irsegidx
        SurfIdx(icellRayIdx, j, 1) = AsyRay(iAsyRay)%PinRaySurf(2, l) !OutSurface
        SurfIdx(icellRayIdx, j, 2) = AsyRay(iAsyRay)%PinRaySurf(1, l) !Insurface
      ENDDO    !End of Pin Ray Seep, l
    ENDDO    !End of Asy Ray Swwep, k
    nTotRaySeg(j) = irsegidx
    nTotCellRay(j) = icellRayIdx
  ENDDO    !End of Core Ray Sweep, j
ELSEIF(FastMocLV .EQ. 1) THEN
   FastRay => RayInfo%FastCoreRayDat(i, iz)
   CellRay1D => RayInfo%CellRay1D
   DO j = 1, nCoreRay
     irsegidx = 0
     nTotRaySeg(j) = FastRay%nTotRaySeg(j)
     nTotCellRay(j) = FastRay%nTotCellRay(j)
     DO l = 1, FastRay%nTotCellRay(j)
       PinIdx(l, j) = FastRay%PinIdx(l, j)
       CellRayIdxSt(l, j, 1) = FastRay%CellRayIdxSt(l, j, 1)
       CellRayIdxSt(l, j, 2) = FastRay%CellRayIdxSt(l, j, 2)
       SurfIdx(l, j, 1) = FastRay%SurfIdx(l, j, 1)
       SurfIdx(l, j, 2) = FastRay%SurfIdx(l, j, 2)
       ipin =  PinIdx(l, j); icel = Pin(ipin)%Cell(iz)
       FsrIdxSt = Pin(ipin)%FsrIdxSt
       DO k = FastRay%Ray1DIdx(1, l, j), FastRay%Ray1DIdx(2, l, j)
         irsegidx = irsegidx + 1
         ireg = FsrIdxSt + CellRay1D%LocalFsrIdx(K) - 1
         FsrIdx(irsegidx, j) = ireg
         tau = - CellRay1D%LenSeg(k) * xst(ireg)
         OptLenList(irsegidx, j) = tau
         ExpAppIdx(irsegidx, j) = max(INT(tau), -40000)
         ExpAppIdx(irsegidx, j) = min(0, ExpAppIdx(irsegidx, j))
       ENDDO
     ENDDO
   ENDDO
ELSEIF(FastMocLv .EQ. 2) THEN
  FastRay => RayInfo%FastCoreRayDat(i, iz)
  CellRay1D => RayInfo%CellRay1D
  DO j = 1, nCoreRay
    nTotRaySeg(j) = FastRay%nTotRaySeg(j)
    nTotCellRay(j) = FastRay%nTotCellRay(j)
    DO l = 1, FastRay%nTotCellRay(j)
      PinIdx(l, j) = FastRay%PinIdx(l, j)
      CellRayIdxSt(l, j, 1) = FastRay%CellRayIdxSt(l, j, 1)
      CellRayIdxSt(l, j, 2) = FastRay%CellRayIdxSt(l, j, 2)
      SurfIdx(l, j, 1) = FastRay%SurfIdx(l, j, 1)
      SurfIdx(l, j, 2) = FastRay%SurfIdx(l, j, 2)
    ENDDO
    FastRaySeg => RayInfo%FastCoreRayDat(i, iz)%RaySeg(j)
    DO l = 1, FastRay%nTotRaySeg(j)
      ireg = FastRaySeg%FsrIdx(l)
      FsrIdx(l, j) = ireg
      tau = - FastRaySeg%LenSeg(l) * xst(ireg)
      OptLenList(l, j) = tau
      ExpAppIdx(l, j) = max(INT(tau), -40000)
      ExpAppIdx(l, j) = min(0, ExpAppIdx(l, j))
    ENDDO
  ENDDO
ENDIF

!Approximate 1-exp
DO j = 1, nCoreRay
  DO l = 1, nTotRaySeg(j)
    DO ipol = 1, nPolarANg
      ExpAppPolar(ipol, l, j) = expa(ExpAppIdx(l, j), ipol)*optlenlist(l, j) + expb(ExpAppIdx(l, j), ipol)
    ENDDO
  ENDDO
ENDDO

ALLOCATE(phiobdPolar(1:nPolarAng))
DO irot = 1, 2
  PhiAnginSvIdx = RayInfo%PhiAngInSvIdx(iRotRay ,irot)
  PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay ,irot)
  phiobdPolar = PhiAngIn(:, PhiAnginSvIdx)
  jinc = 1; jbeg = 1; jend = nCoreRay
  IF(irot .eq. 2) THEN !Backward Sweep
    jinc = -1; jbeg = nCoreRay; jend = 1
  ENDIF
  DO j = jbeg, jend, jinc
    idir = RotRay(i)%DIR(j); iazi = CoreRay(RotRay(irotray)%RayIdx(j))%iang
    IF(irot .eq. 2) idir = mp(idir)  !Reverse the sweep direction
    IF (lJout) THEN
      wt2(1:4) = wtang2(ipol, iazi, 1:4)
    ENDIF
    nRaySeg = nTotRaySeg(j)
    IF(idir .eq. 1) THEN  !Forward Sweep
      PhiAngOutPolar(:, 1) = phiobdPolar(:)
      DO ir = 1, nRaySeg
        ireg = FsrIdx(ir, j)
        DO ipol = 1, nPolarAng
          wt = wtang(ipol, iazi)
          phid = (PhiAngOutPolar(ipol, ir) - src(ireg)) * ExpAppPolar(ipol, ir, j)
          PhiAngOutPolar(ipol, ir+1) = PhiAngOutPolar(ipol, ir) - phid
          phis(ireg) = phis(ireg) + wt*phid
        ENDDO
      ENDDO
      phiobdPolar(:) = PhiAngOutPolar(:, nRaySeg+1)
      !Surface
      IF (ljout) THEN
        DO ir = 1, nTotCellRay(j)
          DO ipol = 1, nPolarANg
            wt = wtang(ipol, iazi)
            icel = PinIdx(ir, j); isurf = SurfIdx(ir, j, 1)
            Jout(2, isurf, icel) = Jout(2, isurf, icel) + wt * PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 1)+1)
            Jout(3, isurf, icel) = Jout(3, isurf, icel) + wt2(isurf) * PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 1)+1)
            isurf = SurfIdx(ir, j, 2)
            Jout(1, isurf, icel) = Jout(1, isurf, icel) + wt * PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 2))
            Jout(3, isurf, icel) = Jout(3, isurf, icel) + wt2(isurf) * PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 2))
          ENDDO
        ENDDO
      ENDIF
    ELSE
      PhiAngOutPolar(:, nRaySeg+2) = phiobdPolar(:)
      ir = nRaySeg + 1
      DO ir1 = 1, nRaySeg
        ir = ir - 1
        ireg = FsrIdx(ir, j)
        DO ipol = 1, nPolarAng
          wt = wtang(ipol, iazi)
          phid = (PhiAngOutPolar(ipol, ir + 2) - src(ireg)) * ExpAppPolar(ipol, ir, j)
          PhiAngOutPolar(ipol, ir+1) = PhiAngOutPolar(ipol, ir + 2) - phid
          phis(ireg) = phis(ireg) + wt * phid
        ENDDO
      ENDDO
      phiobdPolar(:) = PhiAngOutPolar(:, 2)
      IF (lJout) THEN
        DO ir = 1, nTotCellRay(j)
          DO ipol = 1, nPolarAng
            wt = wtang(ipol, iazi)
            icel = PinIdx(ir, j); isurf = SurfIdx(ir, j, 2)
            Jout(2, isurf, icel) = Jout(2, isurf, icel) + wt * PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 2)+1)
            Jout(3, isurf, icel) = Jout(3, isurf, icel) + wt2(isurf) * PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 2)+1)
            isurf = SurfIdx(ir, j, 1)
            Jout(1, isurf, icel) = Jout(1, isurf, icel) + wt * PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 1)+2)
            Jout(3, isurf, icel) = Jout(3, isurf, icel) + wt2(isurf) * PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 1)+2)
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  ENDDO !End of CoreRay Sweep
  PhiAngIn(:,PhiAngOutSvIdx) = phiobdPolar(:)
ENDDO !Backward and forward Sweep

DEALLOCATE(phiobdPolar)

NULLIFY(AziAng); NULLIFY(PolarAng)
NULLIFY(AsyRay); NULLIFY(CoreRay)
NULLIFY(RotRay); NULLIFY(CellRay)

!Geometry Info Pointing
NULLIFY(Asy); NULLIFY(Pin)
NULLIFY(PinInfo); NULLIFY(Cell)

!Tracking Dat Pointing
NULLIFY(FsrIdx); NULLIFY(ExpAppIdx)
NULLIFY(OptLenList); NULLIFY(ExpAppPolar)
NULLIFY(LenSeg); NULLIFY(LocalFsrIdx)
NULLIFY(Phis); NULLIFY(src)
NULLIFY(xst); NULLIFY(jout)
NULLIFY(PhiAngOutPolar); NULLIFY(PhiAngIn)
NULLIFY(EXPA); NULLIFY(EXPB)
NULLIFY(wtang)

END SUBROUTINE

SUBROUTINE TrackRotRay_AFSS(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, FastMocLv)
USE PARAM
USE TYPEDEF,  ONLY :  RayInfo_Type,      coreinfo_type,                                       &
                      Pin_Type,          Asy_Type,        AsyInfo_Type,     PinInfo_Type,     &
                      Cell_Type,                                                              &
                      AziAngleInfo_Type, PolarAngle_Type, ModRayInfo_type,  AsyRayInfo_type,  &
                      CoreRayInfo_Type,  RotRayInfo_Type, CellRayInfo_type, FastCoreRayDat_Type,&
                      TrackingDat_Type,  FastRaySegDat_Type
USE Moc_Mod, ONLY :   nMaxRaySeg,        nMaxCellRay,     nMaxAsyRay,       nMaxCoreRay
USE BasicOperation, ONLY : CP_CA, CP_VA
IMPLICIT NONE

TYPE(RayInfo_Type), INTENT(INOUT) :: RayInfo
TYPE(CoreInfo_Type), INTENT(INOUT) :: CoreInfo
TYPE(TrackingDat_Type), INTENT(INOUT) :: TrackingDat
LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz
INTEGER, INTENT(IN) :: FastMocLv

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Asy_Type), POINTER :: Asy(:)
!TYPE(AsyInfo_Type), POINTER :: AsyInfo
TYPE(PinInfo_Type), POINTER :: PinInfo(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(AziAngleInfo_Type), POINTER :: AziAng(:)
TYPE(PolarAngle_Type), POINTER :: PolarAng(:)
!TYPE(ModRayInfo_type), POINTER :: ModRay
TYPE(AsyRayInfo_type), POINTER :: AsyRay(:)
TYPE(CoreRayInfo_Type), POINTER :: CoreRay(:)
TYPE(RotRayInfo_Type), POINTER :: RotRay(:)
TYPE(FastCoreRayDat_Type), POINTER :: FastRay
TYPE(CellRayInfo_Type), POINTER :: CellRay
TYPE(CellRayInfo_Type), POINTER :: CellRay1D
TYPE(FastRaySegDat_Type), POINTER :: FastRaySeg

REAL, POINTER :: LenSeg(:)
INTEGER, POINTER :: LocalFsrIdx(:)
INTEGER, POINTER :: FsrIdx(:, :),  ExpAppIdx(:, :)
REAL, POINTER :: OptLenList(:, :), ExpApp(:, :), ExpAppPolar(:,:,:)
REAL, POINTER :: phis(:), src(:), xst(:), jout(:, :, :)
REAL, POINTER :: PhiAngOut(:), PhiAngIn(:, :)
REAL, POINTER :: PhiAngOutPolar(:,:)

REAL, POINTER :: EXPA(:, :), EXPB(:, :)

INTEGER :: mp(2)
INTEGER :: iazi, ipol, iCoreRay, iasyray, iceray, irayseg, irot, itype, idir     !Ray related index
INTEGER :: nRotRay, nCoreRay, nAsyRay, nPinRay, nRaySeg, FsrIdxSt
INTEGER :: nPolarAng, nAziAng, nPhiAngSv
INTEGER :: ipin, icel, iasy, ireg, isurf                                                  !Geometries index
INTEGER :: irsegidx, icellrayidx, PhiAnginSvIdx, PhiAngOutSvIdx
INTEGER :: nFsr, nxy
INTEGER :: i, j, k, l, m, jbeg, jend, jinc, ir, ir1
INTEGER :: iOmpAzi
INTEGER :: ibcel

INTEGER :: CellRayIdxSt(nMaxCellRay, nMaxCoreRay, 2)
INTEGER :: PinIdx(nMaxCellRay, nMaxCoreRay)
INTEGER :: SurfIdx(nMaxCellRay, nMaxCoreRay, 2)
INTEGER :: nTotRaySeg(nMaxCoreRay)
INTEGER :: nTotCellRay(nMaxCoreRay)

REAL, POINTER :: wtang(:, :)
REAL, POINTER :: Phi1a(:,:,:), Phi2a(:,:,:)

REAL :: wtang2(100,100,4)
REAL :: tau, phiobd, phid, wt, wt2(4)
REAL, ALLOCATABLE :: phiobdPolar(:)

LOGICAL :: lFast

DATA mp /2, 1/

lFast = FALSE
IF(FastMocLv .GT. 0) lFast = .TRUE.

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

!Tracking Dat Pointing
FsrIdx => TrackingDat%FsrIdx
ExpAppIdx => TrackingDat%ExpAppIdx
OptLenList => TrackingDat%OptLenList
ExpApp => TrackingDat%ExpApp
ExpAppPolar => TrackingDat%ExpAppPolar
Phis => TrackingDat%phis; src => TrackingDat%src
phi1a => TrackingDat%phi1a; phi2a => TrackingDat%phi2a
xst => TrackingDat%xst; jout => TrackingDat%jout
PhiAngOut => TrackingDat%PhiAngOut
PhiAngOutPolar => TrackingDat%PhiAngOutPolar
PhiAngIn => TrackingDat%phiAngIn
EXPA => TrackingDat%EXPA
EXPB => TrackingDat%EXPB

Wtang => TrackingDat%wtang
nAziAng = RayInfo%nAziAngle; nPolarAng = RayInfo%nPolarAngle
nPhiAngSv = RayInfo%nPhiAngSv
nFsr = CoreInfo%nCoreFsr; nxy = CoreInfo%nxy

IF(lJout)THEN
  DO ipol = 1, nPolarAng
    DO iazi = 1, nAziAng !(SWNE) ½º¿þÀÎ?
      wtang2(ipol, iazi, 1) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv)
      wtang2(ipol, iazi, 3) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv)
      wtang2(ipol, iazi, 2) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv)
      wtang2(ipol, iazi, 4) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv)
    ENDDO
  ENDDO
ENDIF

i = iRotRay
nCoreRay = RotRay(irotRay)%nRay
IF(.NOT. lFast) THEN
  DO j = 1, nCoreRay    !Core Ray Sweep
    irsegidx = 0;   icellrayidx = 0
    iCoreRay = RotRay(iRotRay)%RayIdx(j)
    nAsyRay = CoreRay(iCoreRay)%nRay
    DO k = 1, nAsyRay    !Assembly Ray Sweep
      iasyray = CoreRay(iCoreRay)%AsyRayIdx(k)
      iasy = CoreRay(iCoreRay)%AsyIdx(k)
      IF(iasy .EQ. 0)  CYCLE   !Skip Dummy Assembly
      nPinRay = AsyRay(iAsyRay)%nCellRay
      itype = Asy(iasy)%PartialAsyFlag
      DO l = 1, nPinRay   !Pin Ray Sweep
        ipin = AsyRay(iAsyRay)%PinIdx(l)      !Local Pin Idx(within Assembly)
        iceray = AsyRay(iAsyRay)%PinRayIdx(l) !Cell Ray Index
        ipin = Asy(iAsy)%GlobalPinIdx(ipin)   !Global Pin Index
        icel = Pin(ipin)%Cell(iz)             !Cell Type
        FsrIdxSt = Pin(ipin)%FsrIdxSt
        ibcel=Cell(icel)%basecellstr
        CellRay => Cell(ibcel)%CellRay(iceray) !Pointing Cell Ray
        icellrayidx = icellrayidx + 1
        PinIdx(icellrayidx, j) = ipin
        CellRayIdxSt(icellrayidx, j, 2) = irsegidx + 1
        nRaySeg = CellRay%nSeg
        LocalFsrIdx => CellRay%LocalFsrIdx
        LenSeg => CellRay%LenSeg
        DO iRaySeg = 1, nRaySeg
          !ireg = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
          ireg = FsrIdxSt + CellRay%LocalFsrIdx(iRaySeg) - 1
          tau = - LenSeg(iRaySeg) * xst(ireg)   !

             !Global Fsr Index
          !tau = -1000._8 * CellRay%LenSeg(iRaySeg) * xst(ireg)   !
          tau = - CellRay%LenSeg(iRaySeg) * xst(ireg)   !
          irsegidx = irsegidx + 1
          FsrIdx(irsegidx, j) = ireg
          OptLenList(irsegidx, j) = tau
          ExpAppIdx(irsegidx, j) = max(INT(tau), -40000)
          ExpAppIdx(irsegidx, j) = min(0, ExpAppIdx(irsegidx, j))
          !OptLenList(irsegidx, j) = OptLenList(irsegidx, j)/1000._8
          continue
        ENDDO   !End of Ray Segments Sweep, irayseg
        CellRayIdxSt(icellrayidx, j, 1) = irsegidx
        SurfIdx(icellRayIdx, j, 1) = AsyRay(iAsyRay)%PinRaySurf(2, l) !OutSurface
        SurfIdx(icellRayIdx, j, 2) = AsyRay(iAsyRay)%PinRaySurf(1, l) !Insurface
      ENDDO    !End of Pin Ray Seep, l
    ENDDO    !End of Asy Ray Swwep, k
    nTotRaySeg(j) = irsegidx
    nTotCellRay(j) = icellRayIdx
  ENDDO    !End of Core Ray Sweep, j
ELSEIF(FastMocLV .EQ. 1) THEN
   FastRay => RayInfo%FastCoreRayDat(i, iz)
   CellRay1D => RayInfo%CellRay1D
   DO j = 1, nCoreRay
     irsegidx = 0
     nTotRaySeg(j) = FastRay%nTotRaySeg(j)
     nTotCellRay(j) = FastRay%nTotCellRay(j)
     DO l = 1, FastRay%nTotCellRay(j)
       PinIdx(l, j) = FastRay%PinIdx(l, j)
       CellRayIdxSt(l, j, 1) = FastRay%CellRayIdxSt(l, j, 1)
       CellRayIdxSt(l, j, 2) = FastRay%CellRayIdxSt(l, j, 2)
       SurfIdx(l, j, 1) = FastRay%SurfIdx(l, j, 1)
       SurfIdx(l, j, 2) = FastRay%SurfIdx(l, j, 2)
       ipin =  PinIdx(l, j); icel = Pin(ipin)%Cell(iz)
       FsrIdxSt = Pin(ipin)%FsrIdxSt
       DO k = FastRay%Ray1DIdx(1, l, j), FastRay%Ray1DIdx(2, l, j)
         irsegidx = irsegidx + 1
         ireg = FsrIdxSt + CellRay1D%LocalFsrIdx(K) - 1
         FsrIdx(irsegidx, j) = ireg
         tau = - CellRay1D%LenSeg(k) * xst(ireg)
         OptLenList(irsegidx, j) = tau
         ExpAppIdx(irsegidx, j) = max(INT(tau), -40000)
         ExpAppIdx(irsegidx, j) = min(0, ExpAppIdx(irsegidx, j))
       ENDDO
     ENDDO
   ENDDO
ELSEIF(FastMocLv .EQ. 2) THEN
  FastRay => RayInfo%FastCoreRayDat(i, iz)
  CellRay1D => RayInfo%CellRay1D
  DO j = 1, nCoreRay
    nTotRaySeg(j) = FastRay%nTotRaySeg(j)
    nTotCellRay(j) = FastRay%nTotCellRay(j)
    DO l = 1, FastRay%nTotCellRay(j)
      PinIdx(l, j) = FastRay%PinIdx(l, j)
      CellRayIdxSt(l, j, 1) = FastRay%CellRayIdxSt(l, j, 1)
      CellRayIdxSt(l, j, 2) = FastRay%CellRayIdxSt(l, j, 2)
      SurfIdx(l, j, 1) = FastRay%SurfIdx(l, j, 1)
      SurfIdx(l, j, 2) = FastRay%SurfIdx(l, j, 2)
    ENDDO
    FastRaySeg => RayInfo%FastCoreRayDat(i, iz)%RaySeg(j)
    DO l = 1, FastRay%nTotRaySeg(j)
      ireg = FastRaySeg%FsrIdx(l)
      FsrIdx(l, j) = ireg
      tau = - FastRaySeg%LenSeg(l) * xst(ireg)
      OptLenList(l, j) = tau
      ExpAppIdx(l, j) = max(INT(tau), -40000)
      ExpAppIdx(l, j) = min(0, ExpAppIdx(l, j))
    ENDDO
  ENDDO
ENDIF

!Approximate 1-exp
DO j = 1, nCoreRay
  DO l = 1, nTotRaySeg(j)
    DO ipol = 1, nPolarANg
      ExpAppPolar(ipol, l, j) = expa(ExpAppIdx(l, j), ipol)*optlenlist(l, j) + expb(ExpAppIdx(l, j), ipol)
      CONTINUE
    ENDDO
  ENDDO
ENDDO

ALLOCATE(phiobdPolar(1:nPolarAng))

DO irot = 1, 2
  PhiAnginSvIdx = RayInfo%PhiAngInSvIdx(iRotRay ,irot)
  PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay ,irot)
  phiobdPolar(:) = PhiAngIn(:,PhiAnginSvIdx)

  jinc = 1; jbeg = 1; jend = nCoreRay
  IF(irot .eq. 2) THEN !Backward Sweep
    jinc = -1; jbeg = nCoreRay; jend = 1
  ENDIF
  DO j = jbeg, jend, jinc
    idir = RotRay(i)%DIR(j); iazi = CoreRay(RotRay(irotray)%RayIdx(j))%iang; iompazi = RotRay(irotray)%OmpRayIdx(j)
    IF(irot .eq. 2) idir = mp(idir)  !Reverse the sweep direction
    !wt = wtang(ipol, iazi)

      IF(lJout)THEN
          wt2(1:4)=wtang2(ipol,iazi,1:4)
      ENDIF
    nRaySeg = nTotRaySeg(j)

    IF(idir .eq. 1) THEN  !Forward Sweep
      PhiAngOutPolar(:, 1) = phiobdPolar(:)
      DO ir = 1, nRaySeg
        ireg = FsrIdx(ir, j)
        DO ipol = 1, nPolarAng
          phid = (PhiAngOutPolar(ipol, ir) - src(ireg)) * ExpAppPolar(ipol, ir, j)
          PhiAngOutPolar(ipol, ir+1) = PhiAngOutPolar(ipol, ir) - phid
          phi1a(ipol, ireg, iompazi) =  phi1a(ipol, ireg, iompazi) + phid
        ENDDO
      ENDDO
      phiobdPolar(:) = PhiAngOutPolar(:, nRaySeg+1)
      !Surface
      IF(ljout) THEN
        DO ir = 1, nTotCellRay(j)
          DO ipol = 1, nPolarANg
            wt = wtang(ipol, iazi)
            icel = PinIdx(ir, j); isurf = SurfIdx(ir, j, 1)
            Jout(2, isurf, icel) = Jout(2, isurf, icel) + wt * PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 1)+1)
            Jout(3, isurf, icel) = Jout(3, isurf, icel) + wt2(isurf) * PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 1)+1)
            isurf = SurfIdx(ir, j, 2)
            Jout(1, isurf, icel) = Jout(1, isurf, icel) + wt * PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 2))
            Jout(3, isurf, icel) = Jout(3, isurf, icel) + wt2(isurf) * PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 2))
          ENDDO
        ENDDO
      ENDIF
    ELSE
      PhiAngOutPolar(:, nRaySeg+2) = phiobdPolar(:)
      ir = nRaySeg + 1
      DO ir1 = 1, nRaySeg
        ir = ir - 1
        ireg = FsrIdx(ir, j)
        DO ipol = 1, nPolarAng
          phid = (PhiAngOutPolar(ipol, ir + 2) - src(ireg)) * ExpAppPolar(ipol, ir, j)
          PhiAngOutPolar(ipol, ir+1) = PhiAngOutPolar(ipol, ir + 2) - phid
          phi2a(ipol, ireg, iompazi) =  phi2a(ipol, ireg, iompazi) + phid
        ENDDO
      ENDDO
      phiobdPolar(:) = PhiAngOutPolar(:, 2)
      IF(lJout) THEN
        DO ir = 1, nTotCellRay(j)
          DO ipol = 1, nPolarAng
            wt = wtang(ipol, iazi)
            icel = PinIdx(ir, j); isurf = SurfIdx(ir, j, 2)
            Jout(2, isurf, icel) = Jout(2, isurf, icel) + wt * PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 2)+1)
            Jout(3, isurf, icel) = Jout(3, isurf, icel) + wt2(isurf) * PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 2)+1)
            isurf = SurfIdx(ir, j, 1)
            Jout(1, isurf, icel) = Jout(1, isurf, icel) + wt * PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 1)+2)
            Jout(3, isurf, icel) = Jout(3, isurf, icel) + wt2(isurf) * PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 1)+2)
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  ENDDO !End of CoreRay Sweep
  PhiAngIn(:,PhiAngOutSvIdx) = phiobdPolar(:)
ENDDO !Backward and forward Sweep


DEALLOCATE(phiobdPolar)
NULLIFY(AziAng); NULLIFY(PolarAng)
NULLIFY(AsyRay); NULLIFY(CoreRay)
NULLIFY(RotRay); NULLIFY(CellRay)

!Geometry Info Pointing
NULLIFY(Asy); NULLIFY(Pin)
NULLIFY(PinInfo); NULLIFY(Cell)

!Tracking Dat Pointing
NULLIFY(FsrIdx); NULLIFY(ExpAppIdx)
NULLIFY(OptLenList); NULLIFY(ExpApp)
NULLIFY(LenSeg); NULLIFY(LocalFsrIdx)
NULLIFY(Phis); NULLIFY(src)
NULLIFY(Phi1a); NULLIFY(Phi2a)
NULLIFY(xst); NULLIFY(jout)
NULLIFY(PhiAngOut); NULLIFY(PhiAngIn)
NULLIFY(EXPA); NULLIFY(EXPB)
NULLIFY(wtang)

END SUBROUTINE
! ------------------------------------------------------------------------------------------------------------
! CNJ Edit : Node Majors, Domain Decomposition, Angular Multigrid Ray Tracing
SUBROUTINE RayTraceNM_OMP(RayInfo, CoreInfo, phisnm, PhiAngInnm, xstnm, srcnm, joutnm, iz, gb, ge, ljout, lDomainDcmp, FastMocLv)

USE PARAM
USE TIMER
USE ALLOCS
USE OMP_LIB
USE TYPEDEF,     ONLY : RayInfo_Type, CoreInfo_type, Pin_Type, Cell_Type, MultigridInfo_Type
USE ITRCNTL_MOD, ONLY : ItrCntl
USE Moc_Mod,     ONLY : TrackingDat
USE geom,        ONLY : ng, nbd
USE PE_MOD,      ONLY : PE
USE CNTL,        ONLY : nTracerCntl

IMPLICIT NONE

TYPE(RayInfo_Type)  :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo

INTEGER :: iz, gb, ge, FastMocLv
LOGICAL :: ljout, lDomainDcmp

REAL, POINTER, DIMENSION(:,:)     :: phisnm, xstnm, srcnm
REAL, POINTER, DIMENSION(:,:,:)   :: PhiAngInnm
REAL, POINTER, DIMENSION(:,:,:,:) :: joutnm

LOGICAL, SAVE :: lfirst
DATA lfirst /TRUE/

TYPE(Cell_Type),          POINTER, DIMENSION(:) :: Cell
TYPE(Pin_Type),           POINTER, DIMENSION(:) :: Pin
TYPE(MultigridInfo_Type), POINTER, DIMENSION(:) :: MultigridInfo

INTEGER :: nAziAng, nPolarAng, nPhiAngSv, nRotray, nFsr, nAsy, nxy, nThread
INTEGER :: iRotRay, ipol, iazi, irot, AziIdx, ithr, FsrIdxSt, ipin, icel, ibd, ilv, iray, ifsr, jfsr, ig
! ----------------------------------------------------

nAsy = CoreInfo%nxya
nFsr = CoreInfo%nCoreFsr
nxy  = CoreInfo%nxy

MultigridInfo => RayInfo%MultigridInfo
nPhiAngSv      = RayInfo%nPhiAngSv
nRotRay        = RayInfo%nRotRay

nThread = PE%nThread

ilv = nTracerCntl%gridStr(ItrCntl%mocit)
IF (ItrCntl%mocit .GT. nTracerCntl%gridNum) ilv = nTracerCntl%MultigridLV

nAziAng   = MultigridInfo(ilv)%nAzi
nPolarAng = MultigridInfo(ilv)%nPolar

CALL omp_set_num_threads(nThread)
! ----------------------------------------------------
IF (lfirst) THEN
  lFirst = FALSE

  DO ithr = 1, nThread
    IF (TrackingDat(ithr)%lAllocNM) CYCLE
    
    CALL Dmalloc(TrackingDat(ithr)%phisnm, ng, nFsr)
    
    IF (.NOT. lDomainDcmp) CALL Dmalloc(TrackingDat(ithr)%Joutnm, 3, ng, nbd, nxy)
    
    TrackingDat(ithr)%lAllocNM = TRUE
  END DO
END IF
! ----------------------------------------------------
!$OMP PARALLEL PRIVATE(ithr, irot, iazi, iRotRay, iray, AziIdx)
ithr = omp_get_thread_num() + 1

TrackingDat(ithr)%phisnm      = ZERO
TrackingDat(ithr)%EXPA       => MultigridInfo(ilv)%EXPA
TrackingDat(ithr)%EXPB       => MultigridInfo(ilv)%EXPB
TrackingDat(ithr)%srcnm      => srcnm
TrackingDat(ithr)%xstnm      => xstnm
TrackingDat(ithr)%PhiAngInnm => PhiAngInnm

IF (ljout) TrackingDat(ithr)%joutnm = ZERO

DO irot = 1, 2
  DO iazi = 1, nAziAng / 2
    AziIdx = MultigridInfo(ilv)%AziList(iazi)
    
    !$OMP DO SCHEDULE(GUIDED)
    DO iray = 1, RayInfo%RotRayAziList(0, Aziidx)
      iRotRay = RayInfo%RotRayAziList(iray, AziIdx)
      
      CALL TrackRotRayNM(RayInfo, CoreInfo, TrackingDat(ithr), ljout, iRotRay, iz, ilv, irot, gb, ge)
    END DO
    !$OMP END DO NOWAIT
  END DO
END DO
!$OMP END PARALLEL
! ----------------------------------------------------
phisnm(gb:ge, :) = ZERO
IF(ljout) joutnm(:, gb:ge, :, :) = ZERO

DO ithr = 1, nThread
  DO ifsr = 1, nFsr
    DO ig = gb, ge
      phisnm(ig, ifsr) = phisnm(ig, ifsr) + TrackingDat(ithr)%phisnm(ig, ifsr)
    END DO
  END DO
  
  IF (.NOT. ljout) CYCLE
  
  DO ipin = 1, nxy
    DO ibd = 1, nbd
      DO ig = gb, ge
        joutnm(:, ig, ibd, ipin) = joutnm(:, ig, ibd, ipin) + TrackingDat(ithr)%joutnm(:, ig, ibd, ipin)
      END DO
    END DO
  END DO
END DO
! ----------------------------------------------------
Cell => CoreInfo%CellInfo
Pin  => CoreInfo%Pin

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ipin, FsrIdxSt, icel, ifsr, jfsr, ig)
!$OMP DO SCHEDULE(GUIDED)
DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt
  icel     = Pin(ipin)%Cell(iz);
  
  DO ifsr = 1, Cell(icel)%nFsr
    jfsr = FsrIdxSt + ifsr - 1
    
    DO ig = gb, ge
      phisnm(ig, jfsr) = phisnm(ig, jfsr) / (xstnm(ig, jfsr) * Cell(icel)%vol(ifsr)) + srcnm(ig, jfsr)
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
! ----------------------------------------------------
NULLIFY (Cell)
NULLIFY (Pin)
NULLIFY (MultigridInfo)
! ----------------------------------------------------

END SUBROUTINE RayTraceNM_OMP
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE TrackRotRayNM(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, ilv, irot, gb, ge)

USE PARAM
USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, Asy_Type, PinInfo_Type, Cell_Type, AsyRayInfo_type, CoreRayInfo_Type, RotRayInfo_Type, CellRayInfo_type, TrackingDat_Type
USE CNTL,    ONLY : nTracerCntl

IMPLICIT NONE

TYPE(RayInfo_Type)     :: RayInfo
TYPE(CoreInfo_Type)    :: CoreInfo
TYPE(TrackingDat_Type) :: TrackingDat

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz, ilv, irot, gb, ge

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

INTEGER :: mp(2)

INTEGER, POINTER, DIMENSION(:) :: LocalFsrIdx

REAL :: PhiAngOut(RayInfo%nPolarAngle, gb:ge)
REAL :: phid, tau, locsrc, ExpApp
REAL :: wt(10), wt2(10, 4)

REAL, POINTER, DIMENSION(:)       :: LenSeg
REAL, POINTER, DIMENSION(:,:)     :: phis, src, xst, expa, expb, wtang, wtazi
REAL, POINTER, DIMENSION(:,:,:)   :: PhiAngIn, wtsurf
REAL, POINTER, DIMENSION(:,:,:,:) :: jout

DATA mp /2, 1/
! ----------------------------------------------------

IF (nTracerCntl%lHex) THEN
  CALL HexTrackRotRayNM(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, ilv, irot, gb, ge)
  RETURN
END IF

AsyRay  => RayInfo%AsyRay
CoreRay => RayInfo%CoreRay
RotRay  => RayInfo%RotRay

wtang    => RayInfo%MultigridInfo(ilv)%wtang
wtsurf   => RayInfo%MultigridInfo(ilv)%wtsurf
nPolarAng = RayInfo%MultigridInfo(ilv)%nPolar

Asy     => CoreInfo%Asy
Pin     => CoreInfo%Pin
PinInfo => CoreInfo%Pininfo
Cell    => CoreInfo%CellInfo

phis     => TrackingDat%phisnm
src      => TrackingDat%srcnm
xst      => TrackingDat%xstnm
jout     => TrackingDat%joutnm
PhiAngIn => TrackingDat%phiAngInnm
EXPA     => TrackingDat%EXPA
EXPB     => TrackingDat%EXPB

nCoreRay = RotRay(irotRay)%nRay

PhiAngInSvIdx  = RayInfo%PhiAngInSvIdx(iRotRay, irot)
PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay, irot)

PhiAngOut = PhiAngIn(:, gb:ge, PhiAnginSvIdx)

IF (irot .EQ. 1) THEN
  jbeg = 1; jend = nCoreRay; jinc = 1
ELSE
  jbeg = nCoreRay; jend = 1; jinc = -1
END IF
! ----------------------------------------------------
DO icray = jbeg, jend, jinc
  jcray   = RotRay(iRotRay)%RayIdx(icray)
  nAsyRay = CoreRay(jcray)%nRay
  iazi    = CoreRay(jcray)%iang
  
  idir = RotRay(iRotRay)%dir(icray)
  IF (irot .eq. 2) idir = mp(idir)
  
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
              Jout(1, ig, isurf, ipin) = Jout(1, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
              Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
            END DO
          END DO
        END IF
        
        DO iRaySeg = 1, nRaySeg
          ifsr = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
          
          DO ig = gb, ge
            tau = -LenSeg(iRaySeg) * xst(ig, ifsr)
            
            ExpAppIdx = max(INT(tau), -40000)
            ExpAppIdx = min(0, ExpAppIdx)
            
            locsrc = src(ig, ifsr)
            
            DO ipol = 1, nPolarAng
              ExpApp = EXPA(ipol, ExpAppIdx) * tau + EXPB(ipol, ExpAppIdx)
              
              phid = (PhiAngOut(ipol, ig) - locsrc) * ExpApp
              
              PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
              
              phis(ig, ifsr) = phis(ig, ifsr) + wt(ipol) * phid
            END DO
          END DO
        END DO
        
        IF (lJout) THEN
          isurf = AsyRay(jaray)%PinRaySurf(2, ipray)
          
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(2, ig, isurf, ipin) = Jout(2, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
              Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
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
              Jout(1, ig, isurf, ipin) = Jout(1, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
              Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
            END DO
          END DO
        END IF
        
        DO iRaySeg = nRaySeg, 1, -1
          ifsr = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
          
          DO ig = gb, ge
            tau = -LenSeg(iRaySeg) * xst(ig, ifsr)
            
            ExpAppIdx = max(INT(tau), -40000)
            ExpAppIdx = min(0, ExpAppIdx)
            
            locsrc = src(ig, ifsr)
            
            DO ipol = 1, nPolarAng
              ExpApp = EXPA(ipol, ExpAppIdx) * tau + EXPB(ipol, ExpAppIdx)
              
              phid = (PhiAngOut(ipol, ig) - locsrc) * ExpApp
              
              PhiAngOut(ipol, ig) = PhiAngOut(ipol, ig) - phid
              
              phis(ig, ifsr) = phis(ig, ifsr) + wt(ipol) * phid
            END DO
          END DO
        END DO
        
        IF (lJout) THEN
          isurf = AsyRay(jaray)%PinRaySurf(1, ipray)
          
          DO ig = gb, ge
            DO ipol = 1, nPolarAng
              Jout(2, ig, isurf, ipin) = Jout(2, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt(ipol)
              Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + PhiAngOut(ipol, ig) * wt2(ipol, isurf)
            END DO
          END DO
        END IF
      END DO
    END DO
  END IF
END DO

PhiAngIn(:, gb:ge, PhiAngOutSvIdx) = PhiAngOut
! ----------------------------------------------------

END SUBROUTINE TrackRotRayNM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRay(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, FastMocLv)

USE PARAM
USE TYPEDEF,  ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, PolarAngle_Type, TrackingDat_Type
USE Moc_Mod,  ONLY : nMaxCellRay, nMaxCoreRay
USE HexData,  ONLY : hAsy
USE HexType,  ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexCoreRay, Type_HexRotRay
USE HexData,  ONLY : haRay, hCelBss, hCel, hLgc, hcRay, hRotRay, hAsyTypInfo

IMPLICIT NONE

TYPE(RayInfo_Type), INTENT(INOUT) :: RayInfo
TYPE(CoreInfo_Type), INTENT(INOUT) :: CoreInfo
TYPE(TrackingDat_Type), INTENT(INOUT) :: TrackingDat
LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz
INTEGER, INTENT(IN) :: FastMocLv

INTEGER, POINTER :: FsrIdx(:, :),  ExpAppIdx(:, :)
REAL,    POINTER :: OptLenList(:, :), ExpAppPolar(:,:,:)
REAL,    POINTER :: phis(:), src(:), xst(:), jout(:, :, :)
REAL,    POINTER :: PhiAngOutPolar(:, :), PhiAngIn(:, :)
REAL,    POINTER :: EXPA(:, :), EXPB(:, :)
REAL,    POINTER :: wtang(:, :), phiobdPolar(:)

INTEGER :: iAzi, iPol, iAsyRay, iRot     !Ray related index
INTEGER :: nRotRay, nCoreRay, nAsyRay, nRaySeg, nPolarAng
INTEGER :: iCel, iAsy, iReg, iSurf
INTEGER :: irsegidx, icellrayidx, PhiAnginSvIdx, PhiAngOutSvIdx
INTEGER :: icRay, imRay, jbeg, jend, jinc, ihpRay, iRaySeg, iRaySeg1
INTEGER :: iGeoTyp, iAsyTyp, jhPin, icBss, jcBss, jcRay

INTEGER :: CellRayIdxSt(nMaxCellRay, nMaxCoreRay, 2)
INTEGER :: PinIdx      (nMaxCellRay, nMaxCoreRay)
INTEGER :: SurfIdx     (nMaxCellRay, nMaxCoreRay, 2)
INTEGER :: nTotRaySeg  (nMaxCoreRay)
INTEGER :: nTotCellRay (nMaxCoreRay)

REAL :: tau, phiobd, phid, wt

TYPE(Pin_Type),        POINTER :: Pin(:)
TYPE(PolarAngle_Type), POINTER :: PolarAng(:)
TYPE(Type_HexAsyRay),  POINTER :: haRay_Loc
TYPE(Type_HexCelRay),  POINTER :: CelRay_Loc
TYPE(Type_HexRotRay),  POINTER :: hRotRay_Loc
! ----------------------------------------------------

PolarAng => RayInfo%PolarAngle
nPolarAng = RayInfo%nPolarAngle

hRotRay_Loc => hRotRay(iRotRay)
nCoreRay     = hRotRay_Loc%ncRay

Pin => CoreInfo%Pin

FsrIdx     => TrackingDat%FsrIdx
ExpAppIdx  => TrackingDat%ExpAppIdx
OptLenList => TrackingDat%OptLenList
Phis       => TrackingDat%phis
src        => TrackingDat%src
xst        => TrackingDat%xst
Jout       => TrackingDat%jout
PhiAngIn   => TrackingDat%phiAngIn
EXPA       => TrackingDat%EXPA
EXPB       => TrackingDat%EXPB
wtang      => TrackingDat%wtang

ExpAppPolar    => TrackingDat%ExpAppPolar
PhiAngOutPolar => TrackingDat%PhiAngOutPolar
! ----------------------------------------------------
DO icRay = 1, nCoreRay
  jcRay   = abs(hRotRay_Loc%cRayIdx(icRay))
  nAsyRay = hcRay(jcRay)%nmRay
  
  irSegIdx    = 0
  iCellRayIdx = 0
  
  DO imRay = 1, nAsyRay
    iAsyRay = hcRay(jcRay)%mRayIdx(imRay)
    iAsy    = hcRay(jcRay)%AsyIdx(imRay)
    iAsyTyp = hAsy(iAsy)%AsyTyp
    iGeoTyp = hAsy(iAsy)%GeoTyp
    icBss   = hAsyTypInfo(iAsyTyp)%iBss
    
    haRay_Loc => haRay(iGeoTyp, icBss, iAsyRay)
    
    DO ihpRay = 1, haRay_Loc%nhpRay
      iCellRayIdx = iCellRayIdx + 1
      
      jhPin = haRay_Loc%CelRay(ihpRay)%hPinIdx
      jhPin = hAsy(iAsy)%PinIdxSt + hAsyTypInfo(iAsyTyp)%PinLocIdx(iGeoTyp, jhPin) - 1
      jcBss = Pin(jhPin)%hCelGeo(iz)
      
      CelRay_Loc => haRay(iGeoTyp, jcBss, iAsyRay)%CelRay(ihpRay)
      
      ! Start of Cell
      CellRayIdxSt(iCellRayIdx, icRay, 2) = irSegIdx + 1
      SurfIdx     (iCellRayIdx, icRay, 2) = CelRay_Loc%SurfIdx(1) ! y : Small
      
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
      SurfIdx     (iCellRayIdx, icRay, 1) = CelRay_Loc%SurfIdx(2) ! y : Big
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
      ExpAppPolar(iPol, iRaySeg, icRay) = expa(ExpAppIdx(iRaySeg, icRay), iPol) * optlenlist(iRaySeg, icRay) &
                                        + expb(ExpAppIdx(iRaySeg, icRay), iPol)
    END DO
  END DO
END DO

ALLOCATE(phiobdPolar(1:nPolarAng))
! ----------------------------------------------------
DO iRot = 1, 2
  PhiAnginSvIdx  = RayInfo%PhiAngInSvIdx (iRotRay, iRot)
  PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay, iRot)
  phiobdPolar    = PhiAngIn(:, PhiAnginSvIdx)
  
  jinc = 1
  jbeg = 1
  jend = nCoreRay
  
  IF(iRot .eq. 2) THEN !Backward Sweep
    jinc = -1
    jbeg = nCoreRay
    jend = 1
  ENDIF
  
  DO icRay = jbeg, jend, jinc
    jcRay = hRotRay_Loc%cRayIdx(icRay)
    iAzi  = hcRay(abs(jcRay))%AzmIdx
    
    IF(iRot .eq. 2) jcRay = -jcRay !Reverse the Sweep Direction
    
    nRaySeg = nTotRaySeg(icRay)
    ! ----------------------------------------------------! ----------------------------
    !      1. Forward
    ! ----------------------------
    IF(jcRay > 0) THEN
      PhiAngOutPolar(:, 1) = phiobdPolar(:)
      
      DO iRaySeg = 1, nRaySeg
        iReg = FsrIdx(iRaySeg, icRay)
        
        DO iPol = 1, nPolarAng
          wt   = wtang(iPol, iAzi)
          phid = (PhiAngOutPolar(iPol, iRaySeg) - src(iReg)) * ExpAppPolar(iPol, iRaySeg, icRay)
          
          PhiAngOutPolar(ipol, iRaySeg+1) = PhiAngOutPolar(iPol, iRaySeg) - phid
          
          phis(iReg) = phis(iReg) + wt * phid
        END DO
      END DO
      
      phiobdPolar(:) = PhiAngOutPolar(:, nRaySeg+1)
      
      ! Surface
      IF(ljout) THEN
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
          ENDDO
        END DO
      ENDIF
    ! ----------------------------------------------------
    ELSE
      PhiAngOutPolar(:, nRaySeg+2) = phiobdPolar(:)
      
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
      ENDDO
      
      phiobdPolar(:) = PhiAngOutPolar(:, 2)
      
      ! Surface
      IF(lJout) THEN
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
        ENDDO
      ENDIF
    ENDIF
  ENDDO !End of CoreRay Sweep
  
  PhiAngIn(:,PhiAngOutSvIdx) = phiobdPolar(:)
ENDDO !Backward and forward Sweep
! ----------------------------------------------------

DEALLOCATE (phiobdPolar)

NULLIFY (PolarAng)

!Tracking Dat Pointing
NULLIFY (FsrIdx, ExpAppIdx, OptLenList, ExpAppPolar, Phis, src, xst, PhiAngOutpolar, PhiAngIn, EXPA, EXPB, Pin)

! Hex
NULLIFY (haRay_Loc, CelRay_Loc, hRotRay_Loc)
! ----------------------------------------------------

END SUBROUTINE HexTrackRotRay
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayNM(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, ilv, irot, gb, ge)

USE PARAM
USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, PolarAngle_Type, TrackingDat_Type, Pin_Type
USE Moc_Mod, ONLY : nMaxCellRay, nMaxCoreRay
USE HexData, ONLY : hAsy
USE HexType, ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexCoreRay, Type_HexRotRay
USE HexData, ONLY : haRay, hCelBss, hCel, hLgc, hcRay, hRotRay, hAsyTypInfo

IMPLICIT NONE

TYPE(RayInfo_Type)     :: RayInfo
TYPE(CoreInfo_Type)    :: CoreInfo
TYPE(TrackingDat_Type) :: TrackingDat
LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz, ilv, irot, gb, ge

TYPE(Pin_Type), POINTER, DIMENSION(:) :: Pin

TYPE(Type_HexAsyRay), POINTER :: haRay_Loc
TYPE(Type_HexCelRay), POINTER :: CelRay_Loc
TYPE(Type_HexRotRay), POINTER :: hRotRay_Loc

INTEGER :: iAzi, iPol, icRay, jcRay, iAsyRay, iRaySeg, imRay, ihpRay, iAsy, ifsr, iSurf, jbeg, jend, jinc, ig, iGeoTyp, iAsyTyp, jhPin, icBss, jcBss
INTEGER :: nCoreRay, nAsyRay, nPolarAng, PhiAnginSvIdx, PhiAngOutSvIdx, ExpAppIdx

REAL :: phid, tau, locsrc, ExpApp

REAL :: wt(10)
REAL :: PhiAngOut(RayInfo%nPolarAngle, gb:ge)

REAL, POINTER, DIMENSION(:,:)     :: phis, src, xst, expa, expb, wtang
REAL, POINTER, DIMENSION(:,:,:)   :: PhiAngIn
REAL, POINTER, DIMENSION(:,:,:,:) :: jout
! ----------------------------------------------------

wtang    => RayInfo%MultigridInfo(ilv)%wtang
!wtsurf  => RayInfo%MultigridInfo(ilv)%wtsurf
nPolarAng = RayInfo%MultigridInfo(ilv)%nPolar

hRotRay_Loc => hRotRay(iRotRay)
nCoreRay     = hRotRay_Loc%ncRay

Pin => CoreInfo%Pin

Phis     => TrackingDat%phisnm
src      => TrackingDat%srcnm
xst      => TrackingDat%xstnm
jout     => TrackingDat%joutnm
PhiAngIn => TrackingDat%phiAngInnm
EXPA     => TrackingDat%EXPA
EXPB     => TrackingDat%EXPB

PhiAngInSvIdx  = RayInfo%PhiAngInSvIdx (iRotRay, irot)
PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay, irot)

PhiAngOut = PhiAngIn(:, gb:ge, PhiAnginSvIdx)

IF (iRot .EQ. 1) THEN
  jbeg = 1; jend = nCoreRay; jinc = 1
ELSE
  jbeg = nCoreRay; jend = 1; jinc = -1
END IF
! ----------------------------------------------------
DO icRay = jbeg, jend, jinc
  jcRay   = hRotRay_Loc%cRayIdx(icRay)
  nAsyRay = hcRay(abs(jcRay))%nmRay
  iAzi    = hcRay(abs(jcRay))%AzmIdx
  
  IF (iRot .eq. 2) jcRay = -jcRay !Reverse the Sweep Direction
  
  DO iPol = 1, nPolarAng
    wt(ipol) = wtang(iPol, iAzi)
    
    !IF (lJout) wt2(ipol, 1 : 4) = wtsurf(ipol, iazi, 1 : 4)
  END DO
  ! --------------------------------------------------
  IF(jcRay > 0) THEN
    DO imRay = 1, nAsyRay
      iAsyRay = hcRay(abs(jcRay))%mRayIdx(imRay)
      iAsy    = hcRay(abs(jcRay))%AsyIdx(imRay)
      
      IF (iAsy .EQ. 0) CYCLE
      
      iAsyTyp = hAsy(iAsy)%AsyTyp
      iGeoTyp = hAsy(iAsy)%GeoTyp
      icBss   = hAsyTypInfo(iAsyTyp)%iBss
      
      haRay_Loc => haRay(iGeoTyp, icBss, iAsyRay)
      
      DO ihpRay = 1, haRay_Loc%nhpRay
        jhPin = haRay_Loc%CelRay(ihpRay)%hPinIdx
        jhPin = hAsy(iAsy)%PinIdxSt + hAsyTypInfo(iAsyTyp)%PinLocIdx(iGeoTyp, jhPin) - 1
        jcBss = Pin(jhPin)%hCelGeo(iz)
        
        CelRay_Loc => haRay(iGeoTyp, jcBss, iAsyRay)%CelRay(ihpRay)
        
        IF (lJout) THEN
          iSurf = CelRay_Loc%SurfIdx(1) ! y : small
          
          DO ig = gb, ge
            DO iPol = 1, nPolarAng
              Jout(1, ig, iSurf, jhPin) = Jout(1, ig, isurf, jhPin) + wt(ipol) * PhiAngOut(iPol, ig)
              !Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol, ig)
            END DO
          END DO
        END IF
        
        DO iRaySeg = 1, CelRay_Loc%nSegRay
          ifsr = CelRay_Loc%MshIdx(iRaySeg) + Pin(jhPin)%FsrIdxSt - 1
          
          DO ig = gb, ge
            tau = -CelRay_Loc%SegLgh(iRaySeg) * xst(ig, ifsr) ! Optimum Length
            
            ExpAppIdx = max(INT(tau), -40000)
            ExpAppIdx = min(0, ExpAppIdx)
            
            locsrc = src(ig, ifsr)
            
            DO iPol = 1, nPolarAng
              ExpApp = EXPA(iPol, ExpAppIdx) * tau + EXPB(iPol, ExpAppIdx)
              
              phid = (PhiAngOut(iPol, ig) - locsrc) * ExpApp
              
              PhiAngOut(iPol, ig) = PhiAngOut(iPol, ig) - phid
              
              phis(ig, ifsr) = phis(ig, ifsr) + wt(iPol) * phid
            END DO
          END DO
        END DO
        
        IF (lJout) THEN
          isurf = CelRay_Loc%SurfIdx(2) ! y : Big
          
          DO ig = gb, ge
            DO iPol = 1, nPolarAng
              Jout(2, ig, iSurf, jhPin) = Jout(2, ig, iSurf, jhPin) + wt(iPol) * PhiAngOut(iPol, ig)
              !Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol, ig)
            END DO
          END DO
        END IF
      END DO
    END DO
  ! --------------------------------------------------
  ELSE
    DO imRay = nAsyRay, 1, -1
      iAsyRay = hcRay(abs(jcRay))%mRayIdx(imRay)
      iAsy    = hcRay(abs(jcRay))%AsyIdx(imRay)
      
      IF (iAsy .EQ. 0) CYCLE
      
      iAsyTyp = hAsy(iAsy)%AsyTyp
      iGeoTyp = hAsy(iAsy)%GeoTyp
      icBss   = hAsyTypInfo(iAsyTyp)%iBss
      
      haRay_Loc => haRay(iGeoTyp, icBss, iAsyRay)
      
      DO ihpRay = haRay_Loc%nhpRay, 1, -1
        jhPin = haRay_Loc%CelRay(ihpRay)%hPinIdx
        jhPin = hAsy(iAsy)%PinIdxSt + hAsyTypInfo(iAsyTyp)%PinLocIdx(iGeoTyp, jhPin) - 1
        jcBss = Pin(jhPin)%hCelGeo(iz)
        
        CelRay_Loc => haRay(iGeoTyp, jcBss, iAsyRay)%CelRay(ihpRay)
        
        IF (lJout) THEN
          iSurf = CelRay_Loc%SurfIdx(2) ! y : Big
          
          DO ig = gb, ge
            DO iPol = 1, nPolarAng
              Jout(1, ig, iSurf, jhPin) = Jout(1, ig, isurf, jhPin) + wt(ipol) * PhiAngOut(iPol, ig)
              !Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol, ig)
            END DO
          END DO
        END IF
        
        DO iRaySeg = CelRay_Loc%nSegRay, 1, -1
          ifsr = CelRay_Loc%MshIdx(iRaySeg) + Pin(jhPin)%FsrIdxSt - 1
          
          DO ig = gb, ge
            tau = -CelRay_Loc%SegLgh(iRaySeg) * XsT(ig, ifsr) ! Optimum Length
            
            ExpAppIdx = max(INT(tau), -40000)
            ExpAppIdx = min(0, ExpAppIdx)
            
            locsrc = src(ig, ifsr)
            
            DO iPol = 1, nPolarAng
              ExpApp = EXPA(iPol, ExpAppIdx) * tau + EXPB(iPol, ExpAppIdx)
              
              phid = (PhiAngOut(iPol, ig) - locsrc) * ExpApp
              
              PhiAngOut(iPol, ig) = PhiAngOut(iPol, ig) - phid
              
              phis(ig, ifsr) = phis(ig, ifsr) + wt(iPol) * phid
            END DO
          END DO
        END DO
        
        IF (lJout) THEN
          isurf = CelRay_Loc%SurfIdx(1) ! y : small
          
          DO ig = gb, ge
            DO iPol = 1, nPolarAng
              Jout(2, ig, iSurf, jhPin) = Jout(2, ig, iSurf, jhPin) + wt(iPol) * PhiAngOut(iPol, ig)
              !Jout(3, ig, isurf, ipin) = Jout(3, ig, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol, ig)
            END DO
          END DO
        END IF
      END DO
    END DO
  END IF
END DO

PhiAngIn(:, gb:ge, PhiAngOutSvIdx) = PhiAngOut
! ----------------------------------------------------

END SUBROUTINE HexTrackRotRayNM
! ------------------------------------------------------------------------------------------------------------