#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceP1GM_OMP(RayInfo, CoreInfo, phis, phim, PhiAngIn, xst, src, srcm, jout, iz, ljout, ScatOd, FastMocLv)

USE TIMER
USE ALLOCS
USE OMP_LIB
USE PARAM,   ONLY : TRUE, FALSE, ZERO, ONE
USE TYPEDEF, ONLY : RayInfo_Type, CoreInfo_type, Pin_Type, Cell_Type, AziAngleInfo_Type, PolarAngle_Type
USE Moc_Mod, ONLY : nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay, Expa, Expb, ApproxExp, TrackingDat, wtang, SrcAng, Comp, SrcAng1, SrcAng2, mwt, mwt2
USE geom,    ONLY : nbd
USE PE_MOD,  ONLY : PE
USE CNTL,    ONLY : nTracerCntl

USE BasicOperation, ONLY : CP_CA, CP_VA

IMPLICIT NONE

TYPE (RayInfo_Type) :: RayInfo
TYPE (CoreInfo_Type) :: CoreInfo

REAL, POINTER, DIMENSION(:)     :: phis, xst, src
REAL, POINTER, DIMENSION(:,:)   :: PhiAngIn, srcm, phim
REAL, POINTER, DIMENSION(:,:,:) :: jout

INTEGER :: iz, ScatOd
LOGICAL :: ljout
INTEGER, OPTIONAL :: FastMocLv
! ----------------------------------------------------
LOGICAL, SAVE :: lfirst
DATA lfirst /TRUE/

TYPE (AziAngleInfo_Type), POINTER, DIMENSION(:) :: AziAng
TYPE (PolarAngle_Type),   POINTER, DIMENSION(:) :: PolarAng
TYPE (Cell_Type),         POINTER, DIMENSION(:) :: Cell
TYPE (Pin_Type),          POINTER, DIMENSION(:) :: Pin

INTEGER :: nAziAng, nPolarAng, nPhiAngSv, nRotray, nFsr, nxy, nThread, nod
INTEGER :: ithr, FsrIdxSt, icel, iazi, ipol, iod, iray, ifsr, jfsr, ipin, krot
REAL :: wttmp, wtcos, wtpolar, wtsin2, tmpsrc, ONETHREE, ONEFIVE, ONESEVEN
! ----------------------------------------------------

AziAng   => RayInfo%AziAngle
PolarAng => RayInfo%PolarAngle
nAziAng   = RayInfo%nAziAngle
nPolarAng = RayInfo%nPolarAngle
nPhiAngSv = RayInfo%nPhiAngSv
nRotRay   = RayInfo%nRotRay

nFsr = CoreInfo%nCoreFsr
nxy  = CoreInfo%nxy

nThread = PE%nThread

nod = 2
IF (ScatOd .EQ. 2) nod = 5
IF (ScatOd .EQ. 3) nod = 9
! ----------------------------------------------------
IF (lfirst) THEN
  lFirst = FALSE
  
  CALL ApproxExp(RayInfo%PolarAngle, nPolarAng)
  
  DO ithr = 1, nThread
    IF (TrackingDat(ithr)%lAlloc) CYCLE
    
    CALL dmalloc(TrackingDat(ithr)%FsrIdx,         nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%ExpAppIdx,      nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%OptLenList,     nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%ExpAppPolar,    nPolarAng,  nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%PhiAngOutPolar, nPolarAng,  nMaxRaySeg + 2)
    CALL dmalloc(TrackingDat(ithr)%Phis,           nFsr)
    CALL dmalloc(TrackingDat(ithr)%Jout,           3, nbd, nxy)
    
    TrackingDat(ithr)%Expa  => Expa
    TrackingDat(ithr)%Expb  => Expb
    TrackingDat(ithr)%lAlloc = TRUE
  END DO
  
  CALL dmalloc(wtang, nPolarAng, nAziAng)
  
  DO ipol = 1, nPolarAng
    wttmp = PolarAng(ipol)%weight * PolarAng(ipol)%sinv
    
    DO iazi = 1, nAziAng
      wtang(ipol, iazi) = wttmp * AziAng(iazi)%weight * AziAng(iazi)%del
    END DO
  END DO
  
  ! P1
  DO ithr = 1, nThread
    IF (TrackingDat(ithr)%lAllocP1) CYCLE
    
    CALL Dmalloc(TrackingDat(ithr)%Phim, nod, nFsr)
    
    TrackingDat(ithr)%lAllocP1 = TRUE
  END DO
  
  CALL dmalloc(SrcAng1, nPolarAng, nFsr, nAziAng)
  CALL dmalloc(SrcAng2, nPolarAng, nFsr, nAziAng)
  CALL dmalloc(Comp,    nod, nPolarAng, nAziAng)
  CALL dmalloc(mwt,     nod, nPolarAng, nAziAng)
  CALL dmalloc(mwt2,    nod, nPolarAng, nAziAng)
  
  DO ipol = 1, nPolarAng
    wttmp   = PolarAng(ipol)%sinv
    wtsin2  = PolarAng(ipol)%sinv * PolarAng(ipol)%sinv
    wtcos   = PolarAng(ipol)%cosv
    wtpolar = 1.5_8 * PolarAng(ipol)%cosv * PolarAng(ipol)%cosv - 0.5_8
    
    DO iazi = 1, nAziAng
      comp(1, ipol, iazi) = wttmp * AziAng(iazi)%cosv
      comp(2, ipol, iazi) = wttmp * AziAng(iazi)%sinv
      
      mwt (1:2, ipol, iazi) = comp(1:2, ipol, iazi) * wtang(ipol, iazi)
      mwt2(1:2, ipol, iazi) = -mwt(1:2, ipol, iazi)
            
      IF (scatod .LT. 2) CYCLE
      
      Comp(3, ipol, iazi) = wtpolar
      Comp(4, ipol, iazi) = wtsin2 * (1._8-2._8*AziAng(iazi)%sinv*AziAng(iazi)%sinv)
      Comp(5, ipol, iazi) = wtsin2 * (2._8 * AziAng(iazi)%sinv * AziAng(iazi)%cosv)
      
      mwt (3,   ipol, iazi) = comp(3, ipol, iazi) *  wtang(ipol, iazi)
      mwt (4:5, ipol, iazi) = 0.75_8 * comp(4:5, ipol, iazi) *  wtang(ipol, iazi)
      mwt2(3:5, ipol, iazi) = mwt(3:5, ipol, iazi)
      
      IF (scatod .LT. 3) CYCLE
      
      Comp(6, ipol, iazi) = (5._8 * PolarAng(ipol)%cosv * PolarAng(ipol)%cosv - 1._8) * wttmp * AziAng(iazi)%cosv
      Comp(7, ipol, iazi) = (5._8 * PolarAng(ipol)%cosv * PolarAng(ipol)%cosv - 1._8) * wttmp * AziAng(iazi)%sinv
      Comp(8, ipol, iazi) = (wttmp ** 3._8) * (4._8 * (AziAng(iazi)%cosv ** 3._8) - 3._8 * AziAng(iazi)%cosv)
      Comp(9, ipol, iazi) = (wttmp ** 3._8) * (- 4._8 * (AziAng(iazi)%sinv ** 3._8) + 3._8 * AziAng(iazi)%sinv)
      
      mwt (6:7, ipol, iazi) = 0.375_8 * comp(6:7, ipol, iazi) * wtang(ipol, iazi)
      mwt (8:9, ipol, iazi) = 0.625_8 * comp(8:9, ipol, iazi) * wtang(ipol, iazi)
      mwt2(6:9, ipol, iazi) = -mwt(6:9, ipol, iazi)
    END DO
  END DO
END IF
! ----------------------------------------------------
!$call omp_set_dynamic(FALSE)
!$call omp_set_num_threads(nThread)

DO ithr = 1, nThread
  TrackingDat(ithr)%PhiAngIn => PhiAngIn
  TrackingDat(ithr)%src      => src
  TrackingDat(ithr)%xst      => xst
  TrackingDat(ithr)%wtang    => wtang
  TrackingDat(ithr)%srcm     => srcm
  TrackingDat(ithr)%SrcAng1  => SrcAng1
  TrackingDat(ithr)%SrcAng2  => SrcAng2
  TrackingDat(ithr)%mwt      => mwt
  TrackingDat(ithr)%mwt2     => mwt2
END DO
! ----------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ithr, iazi, ifsr, ipol, tmpsrc)
ithr = omp_get_thread_num() + 1

DO iazi = 1, nAziAng
  DO ifsr = PE%myOmpFsrBeg(ithr), PE%myOmpFsrEnd(ithr)
    DO ipol = 1, nPolarAng
      SrcAng1(ipol, ifsr, iazi) = src(ifsr)
      SrcAng2(ipol, ifsr, iazi) = src(ifsr)
      
      tmpsrc = comp(1, ipol, iazi) * srcm(1, ifsr) + comp(2, ipol, iazi) * srcm(2, ifsr)
      
      SrcAng1(ipol, ifsr, iazi) = SrcAng1(ipol, ifsr, iazi) + tmpsrc
      SrcAng2(ipol, ifsr, iazi) = SrcAng2(ipol, ifsr, iazi) - tmpsrc
      
      IF (scatod .LT. 2) CYCLE
      
      tmpsrc = comp(3, ipol, iazi) * srcm(3, ifsr) + comp(4, ipol, iazi) * srcm(4, ifsr) + comp(5, ipol, iazi) * srcm(5, ifsr)
      
      SrcAng1(ipol, ifsr, iazi) = SrcAng1(ipol, ifsr, iazi) +  tmpsrc
      SrcAng2(ipol, ifsr, iazi) = SrcAng2(ipol, ifsr, iazi) +  tmpsrc
      
      IF (scatod .LT. 3) CYCLE
      
      tmpsrc = comp(6, ipol, iazi) * srcm(6, ifsr) + comp(7, ipol, iazi) * srcm(7, ifsr) + comp(8, ipol, iazi) * srcm(8, ifsr) + comp(9, ipol, iazi) * srcm(9, ifsr)
      
      SrcAng1(ipol, ifsr, iazi) = SrcAng1(ipol, ifsr, iazi) + tmpsrc
      SrcAng2(ipol, ifsr, iazi) = SrcAng2(ipol, ifsr, iazi) - tmpsrc
     END DO
  END DO
END DO
!$OMP END PARALLEL
! ----------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ithr, krot, iray)
ithr = omp_get_thread_num() + 1

TrackingDat(ithr)%phis = ZERO
TrackingDat(ithr)%phim = ZERO
TrackingDat(ithr)%jout = ZERO

DO krot = 1, 2
  IF (nTracerCntl%lHex) THEN
    !$OMP DO SCHEDULE(GUIDED)
    DO iray = 1, nRotRay
      CALL HexTrackRotRayP1GM_OMP(RayInfo, CoreInfo, TrackingDat(ithr), lJout, iray, iz, krot, ScatOd, FastMocLv)
    END DO
    !$OMP END DO NOWAIT
  ELSE
    !$OMP DO SCHEDULE(GUIDED)
    DO iray = 1, nRotRay
      CALL RecTrackRotRayP1GM_OMP(RayInfo, CoreInfo, TrackingDat(ithr), lJout, iray, iz, krot, ScatOd, FastMocLv)
    END DO
    !$OMP END DO NOWAIT
  END IF
END DO
!$OMP END PARALLEL
! ----------------------------------------------------
phis = ZERO
phim = ZERO

IF (ljout) jout = ZERO

DO ithr = 1, nThread
  phis = phis + TrackingDat(ithr)%phis
  phim = phim + TrackingDat(ithr)%phim
  
  IF (.NOT. ljout) CYCLE
  
  jout = jout + TrackingDat(ithr)%jout
END DO
! ----------------------------------------------------
ONETHREE = ONE / 3._8
ONEFIVE  = ONE / 5._8
ONESEVEN = ONE / 7.

Cell => CoreInfo%CellInfo
Pin  => CoreInfo%Pin

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ithr, ipin, FsrIdxSt, icel, ifsr, jfsr, wttmp)
ithr = omp_get_thread_num() + 1

DO ipin = PE%myOmpNxyBeg(ithr), PE%myOmpNxyEnd(ithr)
  FsrIdxSt = Pin(ipin)%FsrIdxSt
  icel     = Pin(ipin)%Cell(iz)
  
  DO ifsr = 1, Cell(icel)%nFsr
    jfsr  = FsrIdxSt + ifsr - 1
    wttmp = ONE / xst(jfsr) / Cell(icel)%vol(ifsr)
    
    phis(jfsr) = phis(jfsr) * wttmp + src(jfsr)
    
    phim(1:2, jfsr) = phim(1:2, jfsr) * wttmp + srcm(1:2, jfsr) * ONETHREE
    
    IF (scatod .LT. 2) CYCLE
    
    phim(3:5, jfsr) = phim(3:5, jfsr) * wttmp + srcm(3:5, jfsr) * ONEFIVE
    
    IF (scatod .LT. 3) CYCLE
    
    phim(6:9, jfsr) = phim(6:9, jfsr) * wttmp + srcm(6:9, jfsr) * ONESEVEN
  END DO
END DO
!$OMP END PARALLEL
! ----------------------------------------------------
NULLIFY (AziAng)
NULLIFY (PolarAng)
NULLIFY (Cell)
NULLIFY (Pin)
! ----------------------------------------------------

END SUBROUTINE RayTraceP1GM_OMP
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceP1GM_AFSS(RayInfo, CoreInfo, phis, phim, PhiAngIn, xst, src, srcm, jout, iz, ljout, ScatOd, FastMocLv)

USE PARAM
USE TYPEDEF, ONLY :   RayInfo_Type,      coreinfo_type,                                       &
                      TrackingDat_Type,                                                       &
                      Pin_Type,          Cell_Type,                                           &
                      AziAngleInfo_Type, PolarAngle_Type
USE Moc_Mod, ONLY :   nMaxRaySeg,        nMaxCellRay,     nMaxAsyRay,       nMaxCoreRay,      &
                      Expa,              Expb,            TrackingDat,                        &
                      ApproxExp,         wtang, &
                      SrcAng1,           SrcAng2,         comp,             mwt
USE geom,           ONLY : nbd
USE cntl,           ONLY : nTracerCntl
USE TIMER
USE BasicOperation, ONLY : CP_CA, CP_VA
USE ALLOCS
USE PE_MOD,  ONLY :   PE
USE OMP_LIB
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phis(:), PhiAngIn(:, :), xst(:), src(:), jout(:, :, :)
REAL, POINTER :: phim(:,:), srcm(:,:)
INTEGER :: iz
LOGICAL :: ljout
INTEGER :: ScatOd
INTEGER, OPTIONAL :: FastMocLv

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
INTEGER :: od, iod

INTEGER :: iRotRay
INTEGER :: tid        !Thread Id
INTEGER :: FsrIdxSt, icel, ireg, iazi, ipol
REAL :: wttmp, wtcos, wtpolar, wtsin2, tempsrc
REAL :: ONETHREE, ONEFIVE, ONESEVEN
INTEGER :: i, j, l, k, m
REAL, ALLOCATABLE :: wtp(:)

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

OD = 2
IF(ScatOd .EQ. 2) OD = 5
IF(ScatOd .EQ. 3) OD = 9
ALLOCATE(wtp(1:OD))

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
     IF(TrackingDat(tid)%lAllocP1) CYCLE
     CALL Dmalloc(TrackingDat(tid)%phi1a, nPolarAng, nFsr, nOmpAng)
     CALL Dmalloc(TrackingDat(tid)%phi2a, nPolarAng, nFsr, nOmpAng)
     TrackingDat(tid)%lAllocP1 = .TRUE.
     IF(TrackingDat(tid)%lAlloc) CYCLE
     CALL Dmalloc(TrackingDat(tid)%FsrIdx, nMaxRaySeg, nMaxCoreRay)
     CALL Dmalloc(TrackingDat(tid)%ExpAppIdx, nMaxRaySeg, nMaxCoreRay)
     CALL Dmalloc(TrackingDat(tid)%OptLenList, nMaxRaySeg, nMaxCoreRay)
     CALL Dmalloc(TrackingDat(tid)%ExpAppPolar, nPolarAng, nMaxRaySeg, nMaxCoreRay)
     CALL Dmalloc(TrackingDat(tid)%PhiAngOutPolar, nPolarAng, nMaxRaySeg+2)
     CALL Dmalloc(TrackingDat(tid)%Jout, 3, nbd, nxy)
     TrackingDat(tid)%lAlloc = .TRUE.
  ENDDO
  ALLOCATE(SrcAng1(nPolarAng, nFsr, nAziAng)); ALLOCATE(SrcAng2(nPolarAng, nFsr, nAziAng))
  ALLOCATE(Comp(Od, nPolarAng, nAziAng))
  ALLOCATE(mwt(Od, nPolarAng, nAziAng))
  ALLOCATE(wtang(nPolarAng, nAziAng))
  DO ipol = 1, nPolarAng
    wttmp = PolarAng(ipol)%weight * PolarAng(ipol)%sinv
    DO iazi = 1, nAziAng
      wtang(ipol, iazi) = wttmp  * AziAng(iazi)%weight * AziAng(iazi)%del
    ENDDO
  ENDDO
  DO ipol = 1, nPolarAng
    wttmp =  PolarAng(ipol)%sinv
    DO iazi = 1, nAziAng
      comp(1, ipol, iazi) = wttmp * AziAng(iazi)%cosv
      comp(2, ipol, iazi) = wttmp * AziAng(iazi)%sinv
      mwt(1:2, ipol, iazi) = comp(1:2, ipol, iazi) *  wtang(ipol, iazi)
    ENDDO
  ENDDO
  IF(ScatOd .GE. 2) THEN
    DO ipol = 1, nPolarAng
      wttmp =  PolarAng(ipol)%sinv
      wtsin2 =  PolarAng(ipol)%sinv * PolarAng(ipol)%sinv
      wtcos  =  PolarAng(ipol)%cosv
      wtpolar =  1.5_8 * PolarAng(ipol)%cosv * PolarAng(ipol)%cosv - 0.5_8
      DO iazi = 1, nAziAng
        Comp(3, ipol, iazi) = wtpolar
        Comp(4, ipol, iazi) = wtsin2 * (1._8-2._8*AziAng(iazi)%sinv*AziAng(iazi)%sinv)
        Comp(5, ipol, iazi) = wtsin2 * (2._8 * AziAng(iazi)%sinv * AziAng(iazi)%cosv)

        mwt(3, ipol, iazi) = comp(3, ipol, iazi) *  wtang(ipol, iazi)
        mwt(4:5, ipol, iazi) = 0.75_8 * comp(4:5, ipol, iazi) *  wtang(ipol, iazi)
      ENDDO
    ENDDO
  ENDIF
  IF(ScatOd .EQ. 3) THEN
    DO ipol = 1, nPolarAng
      wttmp =  PolarAng(ipol)%sinv
      DO iazi = 1, nAziAng
        Comp(6, ipol, iazi) = (5._8 * PolarAng(ipol)%cosv * PolarAng(ipol)%cosv - 1._8) * wttmp * AziAng(iazi)%cosv
        Comp(7, ipol, iazi) = (5._8 * PolarAng(ipol)%cosv * PolarAng(ipol)%cosv - 1._8) * wttmp * AziAng(iazi)%sinv
        Comp(8, ipol, iazi) = (wttmp ** 3._8) * (4._8 * (AziAng(iazi)%cosv ** 3._8) - 3._8 * AziAng(iazi)%cosv)
        Comp(9, ipol, iazi) = (wttmp ** 3._8) * (- 4._8 * (AziAng(iazi)%sinv ** 3._8) + 3._8 * AziAng(iazi)%sinv)

        mwt(6:7, ipol, iazi) = 0.375_8 * comp(6:7, ipol, iazi) * wtang(ipol, iazi)
        mwt(8:9, ipol, iazi) = 0.625_8 * comp(8:9, ipol, iazi) * wtang(ipol, iazi)
      ENDDO
    ENDDO
  ENDIF
ENDIF

!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(nThread)

DO tid = 1, nThread
  TrackingDat(tid)%Expa => Expa; TrackingDat(tid)%Expb => Expb
  TrackingDat(tid)%PhiAngIn => PhiAngIn(:, :)
  TrackingDat(tid)%src => src; TrackingDat(tid)%xst => xst
  TrackingDat(tid)%wtang => WtAng;
  TrackingDat(tid)%SrcAng1 => SrcAng1; TrackingDat(tid)%SrcAng2 => SrcAng2
ENDDO

!Pointing
!Flux and current set to zero
CALL CP_CA(phis, ZERO, nFsr)
CALL CP_CA(phim, ZERO, Od, nFsr)
IF(ljout) CALL CP_CA(Jout, ZERO, 3, nbd, CoreInfo%nxy)

!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(irotray, tid, iazi, ipol, i, j, tempsrc, k, l)
tid = 1
!$  tid = omp_get_thread_num()+1
IF(ScatOd .EQ. 1) THEN
  DO iazi = 1, nAziAng
    DO i = PE%myOmpFsrBeg(tid), PE%myOmpFsrEnd(tid)
      DO ipol = 1, nPolarAng
        SrcAng1(ipol, i, iazi) = src(i); SrcAng2(ipol, i, iazi) = src(i)
        tempsrc = comp(1, ipol, iazi) * srcm(1, i) + comp(2, ipol, iazi) * srcm(2, i)
        SrcAng1(ipol, i, iazi) = SrcAng1(ipol, i, iazi) + tempsrc
        SrcAng2(ipol, i, iazi) = SrcAng2(ipol, i, iazi) - tempsrc
       ENDDO
    ENDDO
  ENDDO
ELSEIF (ScatOd .EQ. 2) THEN
  DO iazi = 1, nAziAng
    DO i = PE%myOmpFsrBeg(tid), PE%myOmpFsrEnd(tid)
      DO ipol = 1, nPolarAng
        SrcAng1(ipol, i, iazi) = src(i); SrcAng2(ipol, i, iazi) = src(i)
        tempsrc = comp(1, ipol, iazi) * srcm(1, i) + comp(2, ipol, iazi) * srcm(2, i)
        SrcAng1(ipol, i, iazi) = SrcAng1(ipol, i, iazi) + tempsrc
        SrcAng2(ipol, i, iazi) = SrcAng2(ipol, i, iazi) - tempsrc
        tempsrc = comp(3, ipol, iazi) * srcm(3, i) + comp(4, ipol, iazi) * srcm(4, i) + comp(5, ipol, iazi) * srcm(5, i)
        SrcAng1(ipol, i, iazi) = SrcAng1(ipol, i, iazi) +  tempsrc
        SrcAng2(ipol, i, iazi) = SrcAng2(ipol, i, iazi) +  tempsrc
      ENDDO
    ENDDO
  ENDDO
ELSEIF (ScatOd .EQ. 3) THEN
  DO iazi = 1, nAziAng
    DO i = PE%myOmpFsrBeg(tid), PE%myOmpFsrEnd(tid)
      DO ipol = 1, nPolarAng
        SrcAng1(ipol, i, iazi) = src(i); SrcAng2(ipol, i, iazi) = src(i)
        tempsrc = comp(1, ipol, iazi) * srcm(1, i) + comp(2, ipol, iazi) * srcm(2, i) + comp(6, ipol, iazi) * srcm(6, i) + comp(7, ipol, iazi) * srcm(7, i) + comp(8, ipol, iazi) * srcm(8, i) + comp(9, ipol, iazi) * srcm(9, i)
        SrcAng1(ipol, i, iazi) = SrcAng1(ipol, i, iazi) + tempsrc
        SrcAng2(ipol, i, iazi) = SrcAng2(ipol, i, iazi) - tempsrc
        tempsrc = comp(3, ipol, iazi) * srcm(3, i) + comp(4, ipol, iazi) * srcm(4, i) + comp(5, ipol, iazi) * srcm(5, i)
        SrcAng1(ipol, i, iazi) = SrcAng1(ipol, i, iazi) +  tempsrc
        SrcAng2(ipol, i, iazi) = SrcAng2(ipol, i, iazi) +  tempsrc
      ENDDO
    ENDDO
  ENDDO
ENDIF
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(irotray, tid, i, j, k, l, iazi, ipol, OmpAng)
tid = 1
!$  tid = omp_get_thread_num()+1
j = tid

DO iazi = 1, nOmpAng
  TrackingDat(j)%phi1a(:,:,iazi) = zero
  TrackingDat(j)%phi2a(:,:,iazi) = zero
END DO

DO i = 1, nxy
  TrackingDat(j)%jout(:, :, i) = zero
END DO

!$OMP BARRIER

IF(lScatBd) THEN
  DO i = OmpRayBeg(tid), OmpRayEnd(tid)
    irotray = i
    CALL RecTrackRotRayP1GM_AFSS(RayInfo, CoreInfo, TrackingDat(tid), lJout, irotray, iz, FastMocLv)
  END DO
ELSE
  DO i = 1, 2
    DO j = OmpRayBegBd(i, tid), OmpRayEndBd(i, tid)
      irotray = j
      CALL RecTrackRotRayP1GM_AFSS(RayInfo, CoreInfo, TrackingDat(tid), lJout, irotray, iz, FastMocLv)
    END DO
  END DO
END IF

!$OMP BARRIER

IF(ScatOd .EQ. 1) THEN
  DO j = 1, nThread
    DO iazi = 1, nOmpAng
      DO i = PE%myOmpFsrBeg(tid), PE%myOmpFsrEnd(tid)
        DO ipol = 1, nPolarAng
          OmpAng = OmpMap(j, iazi)
          phis(i) = phis(i) + wtang(ipol, OmpAng) * (TrackingDat(j)%phi1a(ipol, i, iazi) + TrackingDat(j)%phi2a(ipol, i, iazi))
          phim(1:2, i) = phim(1:2, i) + mwt(1:2, ipol, OmpAng) * (TrackingDat(j)%phi1a(ipol, i, iazi) - TrackingDat(j)%phi2a(ipol, i, iazi))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSEIF(ScatOd .EQ. 2) THEN
  DO j = 1, nThread
    DO iazi = 1, nOmpAng
      DO i = PE%myOmpFsrBeg(tid), PE%myOmpFsrEnd(tid)
        DO ipol = 1, nPolarAng
          OmpAng = OmpMap(j, iazi)
          phis(i) = phis(i) + wtang(ipol, OmpAng) * (TrackingDat(j)%phi1a(ipol, i, iazi) + TrackingDat(j)%phi2a(ipol, i, iazi))
          phim(1:2, i) = phim(1:2, i) + mwt(1:2, ipol, OmpAng) * (TrackingDat(j)%phi1a(ipol, i, iazi) - TrackingDat(j)%phi2a(ipol, i, iazi))
          phim(3:5, i) = phim(3:5, i) + mwt(3:5, ipol, OmpAng) * (TrackingDat(j)%phi1a(ipol, i, iazi) + TrackingDat(j)%phi2a(ipol, i, iazi))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSEIF(ScatOd .EQ. 3) THEN
  DO j = 1, nThread
    DO iazi = 1, nOmpAng
      DO i = PE%myOmpFsrBeg(tid), PE%myOmpFsrEnd(tid)
        DO ipol = 1, nPolarAng
          OmpAng = OmpMap(j, iazi)
          phis(i) = phis(i) + wtang(ipol, OmpAng) * (TrackingDat(j)%phi1a(ipol, i, iazi) + TrackingDat(j)%phi2a(ipol, i, iazi))
          phim(1:2, i) = phim(1:2, i) + mwt(1:2, ipol, OmpAng) * (TrackingDat(j)%phi1a(ipol, i, iazi) - TrackingDat(j)%phi2a(ipol, i, iazi))
          phim(3:5, i) = phim(3:5, i) + mwt(3:5, ipol, OmpAng) * (TrackingDat(j)%phi1a(ipol, i, iazi) + TrackingDat(j)%phi2a(ipol, i, iazi))
          phim(6:9, i) = phim(6:9, i) + mwt(6:9, ipol, OmpAng) * (TrackingDat(j)%phi1a(ipol, i, iazi) - TrackingDat(j)%phi2a(ipol, i, iazi))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF

IF(ljout) THEN
  DO j = 1, nThread
    DO i = PE%myOmpNxyBeg(tid), PE%myOmpNxyEnd(tid)
      jout(:, :, i) = jout(:, :, i) + TrackingDat(j)%jout(:, :, i)
    ENDDO
  ENDDO
ENDIF

!$OMP END PARALLEL
ONETHREE = one/3._8; ONEFIVE = one/5._8; ONESEVEN = one/7._8
Cell => CoreInfo%CellInfo; Pin => CoreInfo%Pin

!$OMP PARALLEL DEFAULT(SHARED)  &
!$OMP PRIVATE(l, j, FsrIdxSt, icel, ireg, wttmp, tid)
tid = 1
!$  tid = omp_get_thread_num()+1
IF(ScatOd .EQ. 1) THEN
  DO l = PE%myOmpNxyBeg(tid), PE%myOmpNxyEnd(tid)
    FsrIdxSt = Pin(l)%FsrIdxSt; icel = Pin(l)%Cell(iz);
    DO j = 1, Cell(icel)%nFsr
      ireg = FsrIdxSt + j - 1
      wttmp = 1._8/(xst(ireg)*Cell(icel)%vol(j))
      phis(ireg) = phis(ireg) * wttmp + src(ireg)
      phim(1:2, ireg) = phim(1:2, ireg) * wttmp + srcm(1:2, ireg) * ONETHREE
    ENDDO
  ENDDO
ELSEIF (ScatOd .EQ. 2) THEN
  DO l = PE%myOmpNxyBeg(tid), PE%myOmpNxyEnd(tid)
    FsrIdxSt = Pin(l)%FsrIdxSt; icel = Pin(l)%Cell(iz);
    DO j = 1, Cell(icel)%nFsr
      ireg = FsrIdxSt + j - 1
      wttmp = 1._8/(xst(ireg)*Cell(icel)%vol(j))
      phis(ireg) = phis(ireg) * wttmp + src(ireg)
      phim(1:2, ireg) = phim(1:2, ireg) * wttmp + srcm(1:2, ireg) * ONETHREE
      phim(3:5, ireg) = phim(3:5, ireg) * wttmp + srcm(3:5, ireg) * ONEFIVE
    ENDDO
  ENDDO
ELSEIF (ScatOd .EQ. 3) THEN
  DO l = PE%myOmpNxyBeg(tid), PE%myOmpNxyEnd(tid)
    FsrIdxSt = Pin(l)%FsrIdxSt; icel = Pin(l)%Cell(iz);
    DO j = 1, Cell(icel)%nFsr
      ireg = FsrIdxSt + j - 1
      wttmp = 1._8/(xst(ireg)*Cell(icel)%vol(j))
      phis(ireg) = phis(ireg) * wttmp + src(ireg)
      phim(1:2, ireg) = phim(1:2, ireg) * wttmp + srcm(1:2, ireg) * ONETHREE
      phim(3:5, ireg) = phim(3:5, ireg) * wttmp + srcm(3:5, ireg) * ONEFIVE
      phim(6:9, ireg) = phim(6:9, ireg) * wttmp + srcm(6:9, ireg) * ONESEVEN
    ENDDO
  ENDDO
ENDIF
!$OMP END PARALLEL

DEALLOCATE(wtp)

END SUBROUTINE RayTraceP1GM_AFSS
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayP1GM_OMP(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, krot, ScatOd, FastMocLv)

USE PARAM
USE TYPEDEF,  ONLY :  RayInfo_Type,      coreinfo_type,                                       &
                      Pin_Type,          Asy_Type,        AsyInfo_Type,     PinInfo_Type,     &
                      Cell_Type,                                                              &
                      AziAngleInfo_Type, PolarAngle_Type, ModRayInfo_type,  AsyRayInfo_type,  &
                      CoreRayInfo_Type,  RotRayInfo_Type, CellRayInfo_type, FastCoreRayDat_Type,&
                      TrackingDat_Type,  FastRaySegDat_Type
USE Moc_Mod, ONLY :   nMaxRaySeg,        nMaxCellRay,     nMaxAsyRay,       nMaxCoreRay
USE BasicOperation, ONLY : CP_CA, CP_VA
USE ALLOCS
IMPLICIT NONE

TYPE(RayInfo_Type), INTENT(INOUT) :: RayInfo
TYPE(CoreInfo_Type), INTENT(INOUT) :: CoreInfo
TYPE(TrackingDat_Type), INTENT(INOUT) :: TrackingDat
LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz, krot, ScatOd
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
REAL, POINTER :: phim(:,:)
REAL, POINTER :: PhiAngOutPolar(:, :), PhiAngIn(:, :)
REAL, POINTER :: EXPA(:, :), EXPB(:, :)

INTEGER :: mp(2)
INTEGER :: iazi, ipol, iCoreRay, iasyray, iceray, irayseg, itype, idir     !Ray related index
INTEGER :: nRotRay, nCoreRay, nAsyRay, nPinRay, nRaySeg, FsrIdxSt
INTEGER :: nPolarAng, nAziAng, nPhiAngSv
INTEGER :: ipin, icel, iasy, ireg, isurf                                                  !Geometries index
INTEGER :: irsegidx, icellrayidx, PhiAnginSvIdx, PhiAngOutSvIdx
INTEGER :: nFsr, nxy
INTEGER :: i, j, k, l, m, jbeg, jend, jinc, ir, ir1, iod, ibcel

INTEGER :: CellRayIdxSt(nMaxCellRay, nMaxCoreRay, 2)
INTEGER :: PinIdx(nMaxCellRay, nMaxCoreRay)
INTEGER :: SurfIdx(nMaxCellRay, nMaxCoreRay, 2)
INTEGER :: nTotRaySeg(nMaxCoreRay)
INTEGER :: nTotCellRay(nMaxCoreRay)

REAL, POINTER :: wtang(:, :), mwt(:, :, :), mwt2(:,:,:)
REAL, POINTER :: SrcAng1(:,:,:), SrcAng2(:,:,:)
INTEGER :: od
REAL :: tau, phiobd, phid, wt, wt2(4), wtang2(100,100,4)
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
ExpAppPolar => TrackingDat%ExpAppPolar
Phis => TrackingDat%phis; src => TrackingDat%src
Phim => TrackingDat%phim
xst => TrackingDat%xst; jout => TrackingDat%jout
PhiAngOutPolar => TrackingDat%PhiAngOutPolar
PhiAngIn => TrackingDat%phiAngIn
EXPA => TrackingDat%EXPA
EXPB => TrackingDat%EXPB
wtang => TrackingDat%wtang
mwt => TrackingDat%mwt
mwt2 => TrackingDat%mwt2
SrcAng1 => TrackingDat%SrcAng1
SrcAng2 => TrackingDat%SrcAng2

nAziAng = RayInfo%nAziAngle; nPolarAng = RayInfo%nPolarAngle
nPhiAngSv = RayInfo%nPhiAngSv
nFsr = CoreInfo%nCoreFsr; nxy = CoreInfo%nxy

OD = 2
IF(ScatOd .EQ. 2) OD = 5
IF(ScatOd .EQ. 3) OD = 9

IF(lJout)THEN
    DO ipol = 1, nPolarAng
        DO iazi = 1, nAziAng !(SWNE) 스웨인?
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

DO j = 1, nCoreRay
  DO l = 1, nTotRaySeg(j)
    DO ipol = 1, nPolarAng
      ExpAppPolar(ipol, l, j) = expa(ExpAppIdx(l, j), ipol)*optlenlist(l, j) + expb(ExpAppIdx(l, j), ipol)
    ENDDO
  ENDDO
ENDDO

ALLOCATE(phiobdPolar(nPolarAng))

PhiAnginSvIdx = RayInfo%PhiAngInSvIdx(iRotRay ,krot)
PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay, krot)
phiobdPolar(:) = PhiAngIn(:,PhiAnginSvIdx)
jinc = 1; jbeg = 1; jend = nCoreRay
IF(krot .EQ. 2) THEN !Backward Sweep
  jinc = -1; jbeg = nCoreRay; jend = 1
ENDIF
DO j = jbeg, jend, jinc
  idir = RotRay(i)%DIR(j); iazi = CoreRay(RotRay(irotray)%RayIdx(j))%iang
  IF(krot .EQ. 2) idir = mp(idir)  !Reverse the sweep direction
  !wt = wtang(ipol, iazi)
  nRaySeg = nTotRaySeg(j)
  
  IF(idir .eq. 1) THEN  !Forward Sweep
    PhiAngOutPolar(:, 1) = phiobdPolar(:)
    DO ir = 1, nRaySeg
      ireg = FsrIdx(ir, j)
      DO ipol = 1, nPolarAng
        wt = wtang(ipol, iazi)
        phid = (PhiAngOutPolar(ipol, ir) - SrcAng1(ipol, ireg, iazi)) * ExpAppPolar(ipol, ir, j)
        PhiAngOutPolar(ipol, ir+1) = PhiAngOutPolar(ipol, ir) - phid
        phis(ireg) = phis(ireg) + wt*phid
        DO iod = 1, od
          phim(iod, ireg) = phim(iod, ireg) + mwt(iod, ipol, iazi) * phid
        ENDDO
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
    DO ir1 = 1, nRaySEg
      ir = ir - 1
      ireg = FsrIdx(ir, j)
      DO ipol = 1, nPolarAng
        wt = wtang(ipol, iazi)
        phid = (PhiAngOutPolar(ipol, ir + 2) - SrcAng2(ipol, ireg, iazi)) * ExpAppPolar(ipol, ir, j)
        PhiAngOutPolar(ipol, ir+1) = PhiAngOutPolar(ipol, ir + 2) - phid
        phis(ireg) = phis(ireg) + wt * phid
        DO iod = 1, od
          phim(iod, ireg) = phim(iod, ireg) + mwt2(iod, ipol, iazi) * phid
        ENDDO
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
NULLIFY(Phis); NULLIFY(src); NULLIFY(Phim)
NULLIFY(xst); NULLIFY(jout)
NULLIFY(PhiAngOutPolar); NULLIFY(PhiAngIn)
NULLIFY(EXPA); NULLIFY(EXPB)
NULLIFY(wtang); NULLIFY(mwt); NULLIFY(mwt2)
NULLIFY(SrcAng1); NULLIFY(SrcAng2)

END SUBROUTINE RecTrackRotRayP1GM_OMP
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayP1GM_AFSS(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, FastMocLv)

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
INTEGER :: ipin, icel, ibcel, iasy, ireg, isurf                                  !Geometries index
INTEGER :: irsegidx, icellrayidx, PhiAnginSvIdx, PhiAngOutSvIdx
INTEGER :: nFsr, nxy
INTEGER :: i, j, k, l, m, jbeg, jend, jinc, ir, ir1
INTEGER :: iOmpAzi

INTEGER :: CellRayIdxSt(nMaxCellRay, nMaxCoreRay, 2)
INTEGER :: PinIdx(nMaxCellRay, nMaxCoreRay)
INTEGER :: SurfIdx(nMaxCellRay, nMaxCoreRay, 2)
INTEGER :: nTotRaySeg(nMaxCoreRay)
INTEGER :: nTotCellRay(nMaxCoreRay)

REAL, POINTER :: wtang(:, :)
REAL, POINTER :: Phi1a(:,:,:), Phi2a(:,:,:)

REAL :: wtang2(100,100,4)
REAL :: tau, phiobd, phid, wt
REAL, ALLOCATABLE :: phiobdPolar(:)
REAL, POINTER :: SrcAng1(:,:,:), SrcAng2(:,:,:)

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

SrcAng1 => TrackingDat%SrcAng1
SrcAng2 => TrackingDat%SrcAng2

IF(lJout)THEN
  DO ipol = 1, nPolarAng
    DO iazi = 1, nAziAng !(SWNE) 스웨인?
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
        ibcel = Cell(icel)%BaseCellStr
        FsrIdxSt = Pin(ipin)%FsrIdxSt
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
    nRaySeg = nTotRaySeg(j)

    IF(idir .eq. 1) THEN  !Forward Sweep
      PhiAngOutPolar(:, 1) = phiobdPolar(:)
      DO ir = 1, nRaySeg
        ireg = FsrIdx(ir, j)
        DO ipol = 1, nPolarAng
          phid = (PhiAngOutPolar(ipol, ir) - SrcAng1(ipol, ireg, iazi)) * ExpAppPolar(ipol, ir, j)
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
            isurf = SurfIdx(ir, j, 2)
            Jout(1, isurf, icel) = Jout(1, isurf, icel) + wt * PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 2))
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
          phid = (PhiAngOutPolar(ipol, ir + 2) - SrcAng2(ipol, ireg, iazi)) * ExpAppPolar(ipol, ir, j)
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
            isurf = SurfIdx(ir, j, 1)
            Jout(1, isurf, icel) = Jout(1, isurf, icel) + wt * PhiAngOutPolar(ipol, CellRayIdxSt(ir, j, 1)+2)
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
NULLIFY(SrcAng1); NULLIFY(SrcAng2)

END SUBROUTINE RecTrackRotRayP1GM_AFSS
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RayTraceP1GM_MGD(RayInfo, CoreInfo, phis, phim, PhiAngIn, xst, src, srcm, jout, iz, ScatOd, lJout)

USE PARAM
USE TYPEDEF,        ONLY : RayInfo_Type,    CoreInfo_type,      Pin_Type,           Cell_Type,                      &
                           MultigridInfo_Type
USE CNTL,           ONLY : nTracerCntl
USE ITRCNTL_MOD,    ONLY : ItrCntl
USE MOC_MOD,        ONLY : AziMap,          TrackingDat
USE ALLOCS
USE PE_MOD,         ONLY : PE
USE OMP_LIB
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phis(:), phim(:, :), PhiAngIn(:, :), xst(:)
REAL, POINTER :: src(:), srcm(:, :), jout(:, :, :)
INTEGER :: iz
LOGICAL :: ljout
INTEGER :: ScatOd

LOGICAL, SAVE :: lfirst
DATA lfirst /.TRUE./

TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(MultigridInfo_Type), POINTER :: MultigridInfo(:)
REAL, POINTER :: phia(:, :, :), SrcAng(:, :, :)
REAL, POINTER :: wtang(:, :), mwt(:, :, :), Comp(:, :, :)

INTEGER :: nAziAng, nPolarAng
INTEGER :: nxy, nFsr, nThread, nAzipThread, nThreadpAzi
INTEGER :: tid, od, icel, ireg, ipin, iazi, ipol, iRotRay
INTEGER :: AziIdx, FsrIdxSt
INTEGER :: ilv, i, j, l, k, m
INTEGER, POINTER :: AziList(:, :)
REAL :: wttmp, tempsrc

Cell => CoreInfo%CellInfo
Pin => CoreInfo%Pin
MultigridInfo => RayInfo%MultigridInfo
nxy = CoreInfo%nxy
nFsr = CoreInfo%nCoreFsr
nThread = PE%nThread

ilv = nTracerCntl%gridStr(ItrCntl%mocit)
IF (ItrCntl%mocit .GT. nTracerCntl%gridNum) ilv = nTracerCntl%MultigridLV

nAziAng = MultigridInfo(ilv)%nAzi
nPolarAng = MultigridInfo(ilv)%nPolar
wtang => MultigridInfo(ilv)%wtang
mwt => MultigridInfo(ilv)%mwt
Comp => MultigridInfo(ilv)%Comp

IF (ScatOd .EQ. 1) od = 2
IF (ScatOd .EQ. 2) od = 5
IF (ScatOd .EQ. 3) od = 9

phis = zero
phim = zero
IF (lJout) jout = zero

DO tid = 1, nThread
  TrackingDat(tid)%EXPA => MultigridInfo(ilv)%EXPA
  TrackingDat(tid)%EXPB => MultigridInfo(ilv)%EXPB
  TrackingDat(tid)%xst => xst
  TrackingDat(tid)%PhiAngIn => PhiAngIn
ENDDO

CALL omp_set_num_threads(nThread)

IF (lfirst) THEN
  lFirst = FALSE
  DO tid = 1, nThread
    IF (.NOT. TrackingDat(tid)%lAllocP1) THEN
      CALL Dmalloc(TrackingDat(tid)%phia, RayInfo%nPolarAngle, nFsr, 4)
      CALL Dmalloc(TrackingDat(tid)%SrcAng, RayInfo%nPolarAngle, nFsr, 4)
      CALL Dmalloc(TrackingDat(tid)%Jout, 3, 4, nxy)
      TrackingDat(tid)%lAllocP1 = .TRUE.
    ENDIF
  ENDDO
ENDIF

nAzipThread = max(1, nAziAng / 2 / nThread)
nThreadpAzi = max(1, nThread / nAziAng * 2)

ALLOCATE(AziList(nThread, nAzipThread))

IF (nThreadpAzi .GT. 1) THEN
  DO i = 1, nThreadpAzi
    iAzi = 0
    DO tid = i, nThread, nThreadpAzi
      iAzi = iAzi + 1
      AziList(tid, 1) = iAzi
    ENDDO
  ENDDO
ELSE
  iAzi = 0
  DO tid = 1, nThread
    DO i = 1, nAzipThread
      iAzi = iAzi + 1
      AziList(tid, i) = iAzi
    ENDDO
  ENDDO
ENDIF

DO i = 1, nAzipThread

  !$OMP PARALLEL PRIVATE(iRotRay, iAzi, AziIdx, SrcAng, tempsrc, tid)
  tid = omp_get_thread_num() + 1; iAzi = AziList(tid, i)
  SrcAng => TrackingDat(tid)%SrcAng

  IF (ScatOd .EQ. 1) THEN
    DO ireg = 1, nFsr
      AziIdx = MultigridInfo(ilv)%AziList(iAzi)
      DO ipol = 1, nPolarAng
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = src(ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = src(ireg)
        tempsrc = Comp(1, ipol, AziIdx) * srcm(1, ireg) + Comp(2, ipol, AziIdx) * srcm(2, ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = SrcAng(ipol, ireg, AziMap(AziIdx, 1)) + tempsrc
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = SrcAng(ipol, ireg, AziMap(AziIdx, 2)) - tempsrc
      ENDDO
      AziIdx = MultigridInfo(ilv)%AziList(nAziAng - iAzi + 1)
      DO ipol = 1, nPolarAng
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = src(ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = src(ireg)
        tempsrc = Comp(1, ipol, AziIdx) * srcm(1, ireg) + Comp(2, ipol, AziIdx) * srcm(2, ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = SrcAng(ipol, ireg, AziMap(AziIdx, 1)) + tempsrc
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = SrcAng(ipol, ireg, AziMap(AziIdx, 2)) - tempsrc
      ENDDO
    ENDDO
  ELSEIF (ScatOd .EQ. 2) THEN
    DO ireg = 1, nFsr
      AziIdx = MultigridInfo(ilv)%AziList(iAzi)
      DO ipol = 1, nPolarAng
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = src(ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = src(ireg)
        tempsrc = Comp(1, ipol, AziIdx) * srcm(1, ireg) + Comp(2, ipol, AziIdx) * srcm(2, ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = SrcAng(ipol, ireg, AziMap(AziIdx, 1)) + tempsrc
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = SrcAng(ipol, ireg, AziMap(AziIdx, 2)) - tempsrc
        tempsrc = Comp(3, ipol, AziIdx) * srcm(3, ireg) + Comp(4, ipol, AziIdx) * srcm(4, ireg)                     &
                  + Comp(5, ipol, AziIdx) * srcm(5, ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = SrcAng(ipol, ireg, AziMap(AziIdx, 1)) + tempsrc
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = SrcAng(ipol, ireg, AziMap(AziIdx, 2)) + tempsrc
      ENDDO
      AziIdx = MultigridInfo(ilv)%AziList(nAziAng - iAzi + 1)
      DO ipol = 1, nPolarAng
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = src(ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = src(ireg)
        tempsrc = Comp(1, ipol, AziIdx) * srcm(1, ireg) + Comp(2, ipol, AziIdx) * srcm(2, ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = SrcAng(ipol, ireg, AziMap(AziIdx, 1)) + tempsrc
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = SrcAng(ipol, ireg, AziMap(AziIdx, 2)) - tempsrc
        tempsrc = Comp(3, ipol, AziIdx) * srcm(3, ireg) + Comp(4, ipol, AziIdx) * srcm(4, ireg)                     &
                  + Comp(5, ipol, AziIdx) * srcm(5, ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = SrcAng(ipol, ireg, AziMap(AziIdx, 1)) + tempsrc
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = SrcAng(ipol, ireg, AziMap(AziIdx, 2)) + tempsrc
      ENDDO
    ENDDO
  ELSEIF (ScatOd .EQ. 3) THEN
    DO ireg = 1, nFsr
      AziIdx = MultigridInfo(ilv)%AziList(iAzi)
      DO ipol = 1, nPolarAng
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = src(ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = src(ireg)
        tempsrc = Comp(1, ipol, AziIdx) * srcm(1, ireg) + Comp(2, ipol, AziIdx) * srcm(2, ireg)                     &
                  + Comp(6, ipol, AziIdx) * srcm(6, ireg) + Comp(7, ipol, AziIdx) * srcm(7, ireg)                   &
                  + Comp(8, ipol, AziIdx) * srcm(8, ireg) + Comp(9, ipol, AziIdx) * srcm(9, ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = SrcAng(ipol, ireg, AziMap(AziIdx, 1)) + tempsrc
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = SrcAng(ipol, ireg, AziMap(AziIdx, 2)) - tempsrc
        tempsrc = Comp(3, ipol, AziIdx) * srcm(3, ireg) + Comp(4, ipol, AziIdx) * srcm(4, ireg)                     &
                  + Comp(5, ipol, AziIdx) * srcm(5, ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = SrcAng(ipol, ireg, AziMap(AziIdx, 1)) + tempsrc
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = SrcAng(ipol, ireg, AziMap(AziIdx, 2)) + tempsrc
      ENDDO
      AziIdx = MultigridInfo(ilv)%AziList(nAziAng - iAzi + 1)
      DO ipol = 1, nPolarAng
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = src(ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = src(ireg)
        tempsrc = Comp(1, ipol, AziIdx) * srcm(1, ireg) + Comp(2, ipol, AziIdx) * srcm(2, ireg)                     &
                  + Comp(6, ipol, AziIdx) * srcm(6, ireg) + Comp(7, ipol, AziIdx) * srcm(7, ireg)                   &
                  + Comp(8, ipol, AziIdx) * srcm(8, ireg) + Comp(9, ipol, AziIdx) * srcm(9, ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = SrcAng(ipol, ireg, AziMap(AziIdx, 1)) + tempsrc
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = SrcAng(ipol, ireg, AziMap(AziIdx, 2)) - tempsrc
        tempsrc = Comp(3, ipol, AziIdx) * srcm(3, ireg) + Comp(4, ipol, AziIdx) * srcm(4, ireg)                     &
                  + Comp(5, ipol, AziIdx) * srcm(5, ireg)
        SrcAng(ipol, ireg, AziMap(AziIdx, 1)) = SrcAng(ipol, ireg, AziMap(AziIdx, 1)) + tempsrc
        SrcAng(ipol, ireg, AziMap(AziIdx, 2)) = SrcAng(ipol, ireg, AziMap(AziIdx, 2)) + tempsrc
      ENDDO
    ENDDO
  ENDIF

  TrackingDat(tid)%phia = zero
  IF (lJout) TrackingDat(tid)%Jout = zero

  AziIdx = MultigridInfo(ilv)%AziList(iAzi)
  DO j = mod(tid - 1, nThreadpAzi) + 1, RayInfo%RotRayAziList(0, AziIdx), nThreadpAzi
    iRotRay = RayInfo%RotRayAziList(j, AziIdx)
    CALL RecTrackRotRayP1GM_MGD(RayInfo, CoreInfo, TrackingDat(tid), ljout, iRotRay, iz, ilv)
  ENDDO

  !$OMP END PARALLEL

  DO tid = 1, nThread

    iAzi = AziList(tid, i)
    phia => TrackingDat(tid)%phia

    IF (ScatOd .EQ. 1) THEN
    !$OMP PARALLEL DO PRIVATE(AziIdx) SCHEDULE(GUIDED)
      DO ireg = 1, nFsr
        AziIdx = MultigridInfo(ilv)%AziList(iAzi)
        DO ipol = 1, nPolarAng
          phis(ireg) = phis(ireg) + wtang(ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) + phia(ipol, ireg, AziMap(AziIdx, 2)))
          phim(1:2, ireg) = phim(1:2, ireg) + mwt(1:2, ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) - phia(ipol, ireg, AziMap(AziIdx, 2)))
        ENDDO
        AziIdx = MultigridInfo(ilv)%AziList(nAziAng - iAzi + 1)
        DO ipol = 1, nPolarAng
          phis(ireg) = phis(ireg) + wtang(ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) + phia(ipol, ireg, AziMap(AziIdx, 2)))
          phim(1:2, ireg) = phim(1:2, ireg) + mwt(1:2, ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) - phia(ipol, ireg, AziMap(AziIdx, 2)))
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
    ELSEIF (ScatOd .EQ. 2) THEN
      !$OMP PARALLEL DO PRIVATE(AziIdx) SCHEDULE(GUIDED)
      DO ireg = 1, nFsr
        AziIdx = MultigridInfo(ilv)%AziList(iAzi)
        DO ipol = 1, nPolarAng
          phis(ireg) = phis(ireg) + wtang(ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) + phia(ipol, ireg, AziMap(AziIdx, 2)))
          phim(1:2, ireg) = phim(1:2, ireg) + mwt(1:2, ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) - phia(ipol, ireg, AziMap(AziIdx, 2)))
          phim(3:5, ireg) = phim(3:5, ireg) + mwt(3:5, ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) + phia(ipol, ireg, AziMap(AziIdx, 2)))
        ENDDO
        AziIdx = MultigridInfo(ilv)%AziList(nAziAng - iAzi + 1)
        DO ipol = 1, nPolarAng
          phis(ireg) = phis(ireg) + wtang(ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) + phia(ipol, ireg, AziMap(AziIdx, 2)))
          phim(1:2, ireg) = phim(1:2, ireg) + mwt(1:2, ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) - phia(ipol, ireg, AziMap(AziIdx, 2)))
          phim(3:5, ireg) = phim(3:5, ireg) + mwt(3:5, ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) + phia(ipol, ireg, AziMap(AziIdx, 2)))
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
    ELSEIF (ScatOd .EQ. 3) THEN
      !$OMP PARALLEL DO PRIVATE(AziIdx) SCHEDULE(GUIDED)
      DO ireg = 1, nFsr
        AziIdx = MultigridInfo(ilv)%AziList(iAzi)
        DO ipol = 1, nPolarAng
          phis(ireg) = phis(ireg) + wtang(ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) + phia(ipol, ireg, AziMap(AziIdx, 2)))
          phim(1:2, ireg) = phim(1:2, ireg) + mwt(1:2, ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) - phia(ipol, ireg, AziMap(AziIdx, 2)))
          phim(3:5, ireg) = phim(3:5, ireg) + mwt(3:5, ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) + phia(ipol, ireg, AziMap(AziIdx, 2)))
          phim(6:9, ireg) = phim(6:9, ireg) + mwt(6:9, ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) - phia(ipol, ireg, AziMap(AziIdx, 2)))
        ENDDO
        AziIdx = MultigridInfo(ilv)%AziList(nAziAng - iAzi + 1)
        DO ipol = 1, nPolarAng
          phis(ireg) = phis(ireg) + wtang(ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) + phia(ipol, ireg, AziMap(AziIdx, 2)))
          phim(1:2, ireg) = phim(1:2, ireg) + mwt(1:2, ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) - phia(ipol, ireg, AziMap(AziIdx, 2)))
          phim(3:5, ireg) = phim(3:5, ireg) + mwt(3:5, ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) + phia(ipol, ireg, AziMap(AziIdx, 2)))
          phim(6:9, ireg) = phim(6:9, ireg) + mwt(6:9, ipol, AziIdx) * (phia(ipol, ireg, AziMap(AziIdx, 1)) - phia(ipol, ireg, AziMap(AziIdx, 2)))
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
    ENDIF

    IF (lJout) Jout = Jout + TrackingDat(tid)%Jout

  ENDDO

ENDDO

IF (ScatOd .EQ. 1) THEN
  !$OMP PARALLEL DO PRIVATE(FsrIdxSt, ireg, icel, wttmp) SCHEDULE(GUIDED)
  DO l = 1, nxy
    FsrIdxSt = Pin(l)%FsrIdxSt; icel = Pin(l)%Cell(iz);
    DO j = 1, Cell(icel)%nFsr
      ireg = FsrIdxSt + j - 1
      wttmp = 1._8 / (xst(ireg) * Cell(icel)%vol(j))
      phis(ireg) = phis(ireg) * wttmp + src(ireg)
      phim(1:2, ireg) = phim(1:2, ireg) * wttmp + srcm(1:2, ireg) * rthree
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
ELSEIF (ScatOd .EQ. 2) THEN
  !$OMP PARALLEL DO PRIVATE(FsrIdxSt, ireg, icel, wttmp) SCHEDULE(GUIDED)
  DO l = 1, nxy
    FsrIdxSt = Pin(l)%FsrIdxSt; icel = Pin(l)%Cell(iz);
    DO j = 1, Cell(icel)%nFsr
      ireg = FsrIdxSt + j - 1
      wttmp = 1._8 / (xst(ireg) * Cell(icel)%vol(j))
      phis(ireg) = phis(ireg) * wttmp + src(ireg)
      phim(1:2, ireg) = phim(1:2, ireg) * wttmp + srcm(1:2, ireg) * rthree
      phim(3:5, ireg) = phim(3:5, ireg) * wttmp + srcm(3:5, ireg) * rfive
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
ELSEIF (ScatOd .EQ. 3) THEN
  !$OMP PARALLEL DO PRIVATE(FsrIdxSt, ireg, icel, wttmp) SCHEDULE(GUIDED)
  DO l = 1, nxy
    FsrIdxSt = Pin(l)%FsrIdxSt; icel = Pin(l)%Cell(iz);
    DO j = 1, Cell(icel)%nFsr
      ireg = FsrIdxSt + j - 1
      wttmp = 1._8 / (xst(ireg) * Cell(icel)%vol(j))
      phis(ireg) = phis(ireg) * wttmp + src(ireg)
      phim(1:2, ireg) = phim(1:2, ireg) * wttmp + srcm(1:2, ireg) * rthree
      phim(3:5, ireg) = phim(3:5, ireg) * wttmp + srcm(3:5, ireg) * rfive
      phim(6:9, ireg) = phim(6:9, ireg) * wttmp + srcm(6:9, ireg) * rseven
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
ENDIF

END SUBROUTINE RayTraceP1GM_MGD
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE RecTrackRotRayP1GM_MGD(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, ilv)

USE PARAM
USE TYPEDEF,        ONLY : RayInfo_Type,        Coreinfo_type,      Pin_Type,           Asy_Type,                   &
                           PinInfo_Type,        Cell_Type,          AsyRayInfo_type,    CoreRayInfo_Type,           &
                           RotRayInfo_Type,     CellRayInfo_type,   TrackingDat_Type
USE MOC_MOD,        ONLY : AziMap
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
TYPE(TrackingDat_Type) :: TrackingDat
LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz, ilv

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(PinInfo_Type), POINTER :: PinInfo(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(AsyRayInfo_type), POINTER :: AsyRay(:)
TYPE(CoreRayInfo_Type), POINTER :: CoreRay(:)
TYPE(RotRayInfo_Type), POINTER :: RotRay(:)
TYPE(CellRayInfo_Type), POINTER :: CellRay

INTEGER, POINTER :: LocalFsrIdx(:)
REAL, POINTER :: LenSeg(:)
REAL, POINTER :: phia(:, :, :), SrcAng(:, :, :), xst(:), jout(:, :, :)
REAL, POINTER :: PhiAngIn(:, :)
REAL, POINTER :: EXPA(:, :), EXPB(:, :)
REAL, POINTER :: wtang(:, :), wtsurf(:, :, :)

INTEGER :: mp(2)
INTEGER :: iazi, ipol, iCoreRay, iasyray, iceray, irayseg, idir
INTEGER :: nRotRay, nCoreRay, nAsyRay, nPinRay, nRaySeg, nPolarAng, FsrIdxSt
INTEGER :: ipin, icel, iasy, ireg, isurf
INTEGER :: PhiAnginSvIdx, PhiAngOutSvIdx, AziSvIdx(2)
INTEGER :: i, j, k, l, m, jbeg, jend, jinc, irot
INTEGER :: ibcel

REAL :: wttmp
REAL :: wt(10), wt2(10, 4)

REAL :: PhiAngOut(RayInfo%nPolarAngle)
REAL :: phid, tau, ExpApp
INTEGER :: ExpAppIdx

LOGICAL :: lFast

DATA mp /2, 1/

AsyRay => RayInfo%AsyRay
CoreRay => RayInfo%CoreRay
RotRay => RayInfo%RotRay
wtang => RayInfo%MultigridInfo(ilv)%wtang
wtsurf => RayInfo%MultigridInfo(ilv)%wtsurf
nPolarAng = RayInfo%MultigridInfo(ilv)%nPolar

Asy => CoreInfo%Asy
Pin => CoreInfo%Pin
PinInfo => CoreInfo%Pininfo
Cell => CoreInfo%CellInfo

phia => TrackingDat%phia; SrcAng => TrackingDat%SrcAng
xst => TrackingDat%xst; jout => TrackingDat%Jout
PhiAngIn => TrackingDat%PhiAngIn
EXPA => TrackingDat%EXPA
EXPB => TrackingDat%EXPB

nCoreRay = RotRay(irotRay)%nRay

DO irot = 1, 2
  PhiAngInSvIdx = RayInfo%PhiAngInSvIdx(iRotRay, irot)
  PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay, irot)
  PhiAngOut = PhiAngIn(:, PhiAnginSvIdx)
  jbeg = 1; jend = nCoreRay; jinc = 1
  IF (irot .EQ. 2) THEN
    jbeg = nCoreRay; jend = 1; jinc = -1
  ENDIF
  DO j = jbeg, jend, jinc
    iCoreRay = RotRay(iRotRay)%RayIdx(j)
    nAsyRay = CoreRay(iCoreRay)%nRay
    idir = RotRay(iRotRay)%dir(j); iazi = CoreRay(iCoreRay)%iang
    AziSvIdx = AziMap(iazi, :)
    IF (irot .eq. 2) idir = mp(idir)
    DO ipol = 1, nPolarAng
      wt(ipol) = wtang(ipol, iazi)
      IF (lJout) wt2(ipol, 1 : 4) = wtsurf(ipol, iazi, 1 : 4)
    ENDDO
    IF (idir .EQ. 1) THEN
      DO k = 1, nAsyRay
        iasyray = CoreRay(iCoreRay)%AsyRayIdx(k)
        iasy = CoreRay(iCoreRay)%AsyIdx(k)
        IF (iAsy .EQ. 0) CYCLE
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
          IF (lJout) THEN
            isurf = AsyRay(iAsyRay)%PinRaySurf(1, l)
            DO ipol = 1, nPolarAng
              Jout(1, isurf, ipin) = Jout(1, isurf, ipin) + wt(ipol) * PhiAngOut(ipol)
              Jout(3, isurf, ipin) = Jout(3, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol)
            ENDDO
          ENDIF
          DO iRaySeg = 1, nRaySeg
            ireg = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
            tau = - LenSeg(iRaySeg) * xst(ireg)
            ExpAppIdx = max(INT(tau), -40000)
            ExpAppIdx = min(0, ExpAppIdx)
            DO ipol = 1, nPolarAng
              ExpApp = EXPA(ipol, ExpAppIdx) * tau + EXPB(ipol, ExpAppIdx)
              phid = (PhiAngOut(ipol) - SrcAng(ipol, ireg, AziSvIdx(idir))) * ExpApp
              PhiAngOut(ipol) = PhiAngOut(ipol) - phid
              phia(ipol, ireg, AziSvIdx(idir)) = phia(ipol, ireg, AziSvIdx(idir)) + phid
            ENDDO
          ENDDO
          IF (lJout) THEN
            isurf = AsyRay(iAsyRay)%PinRaySurf(2, l)
            DO ipol = 1, nPolarAng
              Jout(2, isurf, ipin) = Jout(2, isurf, ipin) + wt(ipol) * PhiAngOut(ipol)
              Jout(3, isurf, ipin) = Jout(3, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    ELSE
      DO k = nAsyRay, 1, -1
        iasyray = CoreRay(iCoreRay)%AsyRayIdx(k)
        iasy = CoreRay(iCoreRay)%AsyIdx(k)
        IF (iAsy .EQ. 0) CYCLE
        nPinRay = AsyRay(iAsyRay)%nCellRay
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
            DO ipol = 1, nPolarAng
              Jout(1, isurf, ipin) = Jout(1, isurf, ipin) + wt(ipol) * PhiAngOut(ipol)
              Jout(3, isurf, ipin) = Jout(3, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol)
            ENDDO
          ENDIF
          DO iRaySeg = nRaySeg, 1, -1
            ireg = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
            tau = - LenSeg(iRaySeg) * xst(ireg)
            ExpAppIdx = max(INT(tau), -40000)
            ExpAppIdx = min(0, ExpAppIdx)
            DO ipol = 1, nPolarAng
              ExpApp = EXPA(ipol, ExpAppIdx) * tau + EXPB(ipol, ExpAppIdx)
              phid = (PhiAngOut(ipol) - SrcAng(ipol, ireg, AziSvIdx(idir))) * ExpApp
              PhiAngOut(ipol) = PhiAngOut(ipol) - phid
              phia(ipol, ireg, AziSvIdx(idir)) = phia(ipol, ireg, AziSvIdx(idir)) + phid
            ENDDO
          ENDDO
          IF (lJout) THEN
            isurf = AsyRay(iAsyRay)%PinRaySurf(1, l)
            DO ipol = 1, nPolarAng
              Jout(2, isurf, ipin) = Jout(2, isurf, ipin) + wt(ipol) * PhiAngOut(ipol)
              Jout(3, isurf, ipin) = Jout(3, isurf, ipin) + wt2(ipol, isurf) * PhiAngOut(ipol)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDDO
  PhiAngIn(:, PhiAngOutSvIdx) = PhiAngOut
ENDDO

END SUBROUTINE RecTrackRotRayP1GM_MGD
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTrackRotRayP1GM_OMP(RayInfo, CoreInfo, TrackingDat, ljout, irotray, iz, krot, ScatOd, FastMocLv)

USE TYPEDEF, ONLY : RayInfo_Type, Coreinfo_type, Pin_Type, TrackingDat_Type, Pin_Type
USE Moc_Mod, ONLY : nMaxCellRay, nMaxCoreRay
USE HexData, ONLY : hAsy
USE HexType, ONLY : Type_HexAsyRay, Type_HexCelRay, Type_HexCoreRay, Type_HexRotRay
USE HexData, ONLY : haRay, hcRay, hRotRay, hAsyTypInfo

IMPLICIT NONE

TYPE (RayInfo_Type),     INTENT(INOUT) :: RayInfo
TYPE (CoreInfo_Type),    INTENT(INOUT) :: CoreInfo
TYPE (TrackingDat_Type), INTENT(INOUT) :: TrackingDat

LOGICAL, INTENT(IN) :: ljout
INTEGER, INTENT(IN) :: irotray, iz, krot, ScatOd
INTEGER, INTENT(IN) :: FastMocLv
! ----------------------------------------------------
INTEGER :: iazi, ipol, iaRay, jaRay, iAsy, iSurf, icRay, jbeg, jend, jinc, ihpRay, ifsr, iRaySeg, iRaySeg1, iCel, iReg
INTEGER :: PhiAnginSvIdx, PhiAngOutSvIdx, iGeoTyp, iAsyTyp, jhPin, icBss, jcBss, jcRay, nod, iod, irsegidx, icellrayidx
INTEGER :: nCoreRay, nAsyRay, nPolarAng, nRaySeg

INTEGER :: CellRayIdxSt(nMaxCellRay, nMaxCoreRay, 2)
INTEGER :: PinIdx      (nMaxCellRay, nMaxCoreRay)
INTEGER :: SurfIdx     (nMaxCellRay, nMaxCoreRay, 2)
INTEGER :: nTotRaySeg  (nMaxCoreRay)
INTEGER :: nTotCellRay (nMaxCoreRay)

INTEGER, POINTER, DIMENSION(:,:) :: FsrIdx, ExpAppIdx

REAL :: tau, phid, wt

REAL :: locphiout(RayInfo%nPolarAngle)

REAL, POINTER, DIMENSION(:)     :: phis, src, xst
REAL, POINTER, DIMENSION(:,:)   :: EXPA, EXPB, wtang, phim, OptLenList, PhiAngOutPolar
REAL, POINTER, DIMENSION(:,:,:) :: jout, mwt, mwt2, SrcAng1, SrcAng2, ExpAppPolar

TYPE (Pin_Type), POINTER, DIMENSION(:) :: Pin

TYPE (Type_HexAsyRay),  POINTER :: haRay_Loc
TYPE (Type_HexCelRay),  POINTER :: CelRay_Loc
TYPE (Type_HexRotRay),  POINTER :: hRotRay_Loc
! ----------------------------------------------------

nPolarAng      = RayInfo%nPolarAngle
PhiAngInSvIdx  = RayInfo%PhiAngInSvIdx (iRotRay, krot)
PhiAngOutSvIdx = RayInfo%PhiangOutSvIdx(iRotRay, krot)

hRotRay_Loc => hRotRay(iRotRay)
nCoreRay     = hRotRay_Loc%ncRay

Pin => CoreInfo%Pin

phis           => TrackingDat%phis
phim           => TrackingDat%phim
src            => TrackingDat%src
xst            => TrackingDat%xst
jout           => TrackingDat%jout
EXPA           => TrackingDat%EXPA
EXPB           => TrackingDat%EXPB
wtang          => TrackingDat%wtang
mwt            => TrackingDat%mwt
mwt2           => TrackingDat%mwt2
SrcAng1        => TrackingDat%SrcAng1
SrcAng2        => TrackingDat%SrcAng2
locphiout       = TrackingDat%PhiAngIn(:, PhiAnginSvIdx)
FsrIdx         => TrackingDat%FsrIdx
ExpAppIdx      => TrackingDat%ExpAppIdx
OptLenList     => TrackingDat%OptLenList
ExpAppPolar    => TrackingDat%ExpAppPolar
PhiAngOutPolar => TrackingDat%PhiAngOutPolar

nod = 2
IF (ScatOd .EQ. 2) nod = 5
IF (ScatOd .EQ. 3) nod = 9
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
      ExpAppPolar(iPol, iRaySeg, icRay) = expa(ExpAppIdx(iRaySeg, icRay), iPol) * optlenlist(iRaySeg, icRay) &
                                        + expb(ExpAppIdx(iRaySeg, icRay), iPol)
    END DO
  END DO
END DO
! ----------------------------------------------------
IF (krot .EQ. 1) THEN
  jbeg = 1; jend = nCoreRay; jinc = 1
ELSE
  jbeg = nCoreRay; jend = 1; jinc = -1
END IF

DO icRay = jbeg, jend, jinc
  jcRay = hRotRay_Loc%cRayIdx(icRay)
  iAzi  = hcRay(abs(hRotRay_Loc%cRayIdx(icRay)))%AzmIdx
  
  IF (krot .EQ. 2) jcRay = -jcRay !Reverse the Sweep Direction
  
  nRaySeg = nTotRaySeg(icRay)
  ! ----------------------------------------------------
  IF (jcRay > 0) THEN
    PhiAngOutPolar(:, 1) = locphiout(:)
    
    DO iRaySeg = 1, nRaySeg
      iReg = FsrIdx(iRaySeg, icRay)
      
      DO iPol = 1, nPolarAng
        wt   = wtang(iPol, iAzi)
        phid = (PhiAngOutPolar(iPol, iRaySeg) - SrcAng1(iPol, iReg, iAzi)) * ExpAppPolar(iPol, iRaySeg, icRay)
        
        PhiAngOutPolar(ipol, iRaySeg+1) = PhiAngOutPolar(iPol, iRaySeg) - phid
        
        phis(iReg) = phis(iReg) + wt * phid
        
        DO iod = 1, nod
          phim(iod, iReg) = phim(iod, iReg) + mwt(iod, iPol, iAzi) * phid
        END DO
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
    ENDIF
  ! ----------------------------------------------------
  ELSE
    PhiAngOutPolar(:, nRaySeg+2) = locphiout(:)
    
    iRaySeg = nRaySeg + 1
    
    DO iRayseg1 = 1, nRaySeg
      iRaySeg = iRaySeg - 1
      iReg    = FsrIdx(iRaySeg, icRay)
      
      DO iPol = 1, nPolarAng
        wt   = wtang(iPol, iAzi)
        phid = (PhiAngOutPolar(iPol, iRaySeg + 2) - SrcAng2(iPol, iReg, iAzi)) * ExpAppPolar(iPol, iRaySeg, icRay)
        
        PhiAngOutPolar(iPol, iRaySeg+1) = PhiAngOutPolar(iPol, iRaySeg + 2) - phid
        
        phis(iReg) = phis(iReg) + wt * phid
        
        DO iod = 1, nod
          phim(iod, iReg) = phim(iod, iReg) + mwt2(iod, iPol, iAzi) * phid
        END DO
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
NULLIFY (wtang)
NULLIFY (phis)
NULLIFY (src)
NULLIFY (xst)
NULLIFY (EXPA)
NULLIFY (EXPB)
NULLIFY (Pin)
NULLIFY (FsrIdx)
NULLIFY (ExpAppIdx)
NULLIFY (OptLenList)
NULLIFY (jout)
NULLIFY (ExpAppPolar)
NULLIFY (PhiAngOutPolar)
NULLIFY (hRotRay_Loc)
NULLIFY (CelRay_Loc)
NULLIFY (hRotRay_Loc)
NULLIFY (phim)
NULLIFY (mwt)
NULLIFY (mwt2)
NULLIFY (SrcAng1)
NULLIFY (SrcAng2)
! ----------------------------------------------------

END SUBROUTINE HexTrackRotRayP1GM_OMP
! ------------------------------------------------------------------------------------------------------------