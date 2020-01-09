#include <defines.h>
SUBROUTINE RayTrace(RayInfo, CoreInfo, phis, PhiAngIn, xst, src, jout, iz, ljout, FastMocLv, lAFSS)
USE PARAM
USE TYPEDEF,  ONLY :  RayInfo_Type,      coreinfo_type,                                       &
                      Pin_Type,          Asy_Type,        AsyInfo_Type,     PinInfo_Type,     &
                      Cell_Type,                                                              &
                      AziAngleInfo_Type, PolarAngle_Type, ModRayInfo_type,  AsyRayInfo_type,  &
                      CoreRayInfo_Type,  RotRayInfo_Type, CellRayInfo_type
USE Moc_Mod, ONLY :   nMaxRaySeg,        nMaxCellRay,     nMaxAsyRay,       nMaxCoreRay,      &
                      Expa,              Expb,                                                &
                      ApproxExp,         RayTrace_OMP,    RayTrace_OMP_AFSS
USE PE_MOD,  ONLY :   PE
#ifdef MPI_ENV
USE MPICOMM_MOD, ONLY : REDUCE,          BCAST
#endif
USE BasicOperation, ONLY : CP_CA, CP_VA
USE ALLOCS
IMPLICIT NONE
!Input Argumentsww
TYPE(RayInfo_Type) :: RayInfo
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phis(:), PhiAngIn(:, :), xst(:), src(:), jout(:, :, :)
!REAL :: phis(:), PhiAngIn(:, :), xst(:), src(:), jout(:, :, :)
INTEGER :: iz
LOGICAL :: ljout
INTEGER, OPTIONAL :: FastMocLv
LOGICAL, OPTIONAL :: lAFSS

INTEGER :: FastMOCLv0
LOGICAL :: lAFSS0

FastMocLv0 = 0
lAFSS0 = .FALSE.
IF(Present(FastMocLv)) FastMocLv0 = FastMocLv
IF(Present(lAFSS)) lAFSS0 = lAFSS

  IF (.NOT. lAFSS0) THEN
    CALL RayTrace_OMP(RayInfo, CoreInfo, phis, PhiAngIn, xst, src, jout, iz, ljout, FastMocLv0, lAFSS0)
  ELSE
    CALL RayTrace_OMP_AFSS(RayInfo, CoreInfo, phis, PhiAngIn, xst, src, jout, iz, ljout, FastMocLv0, lAFSS0)
  ENDIF

END SUBROUTINE 

!--- CNJ Edit : Domain Decomposition    
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
      
END SUBROUTINE 
    
SUBROUTINE MPI_PhiAngInReduce(RayInfo, PhiAngIn1g, n1, n2, PE)
USE PARAM
USE TYPEDEF,     ONLY : RayInfo_Type, PE_TYPE
#ifdef MPI_ENV
USE MPICOMM_MOD, ONLY : REDUCE, MPI_SYNC
#endif
USE BasicOperation, ONLY : CP_CA, CP_VA
IMPLICIT NONE
TYPE(RayInfo_Type) :: RayInfo
TYPE(PE_TYPE) :: PE
INTEGER :: n1, n2
REAL :: PhiAngIn1g(n1, n2)

!REAL :: Buf0(n1,n2), buf1(n1,n2)
!INTEGER :: Lists(n2)
REAL, POINTER :: Buf0(:, :), buf1(:, :)
INTEGER, POINTER :: Lists(:)

INTEGER :: i, j, k1, k2
INTEGER :: n
#ifdef MPI_ENV
!IF(PE%nRTproc .EQ. 1) RETURN
ALLOCATE(Buf0(n1,n2), buf1(n1, n2))
ALLOCATE(Lists(n2))
CALL CP_CA(Lists(1:n2), 0, n2)
CALL CP_CA(buf0(1:n1,1:n2), 0._8, n1, n2)
CALL CP_CA(buf1(1:n1,1:n2), 0._8, n1, n2)
!Lists = 0
!Buf0 = 0  
!Buf1 = 0
DO i =  PE%myRayBeg, PE%myRayEnd 
  k1 = RayInfo%PhiAngInSvIdx(i,1);   k2 = RayInfo%PhiangOutSvIdx(i,1)
  IF(k1 .NE. 1) Buf0(:, k1) = PhiAngIn1g(:, k1)
  IF(k2 .NE. 1) Buf0(:, k1) = PhiAngIn1g(:, k2)
  k1 = RayInfo%PhiAngInSvIdx(i,2);   k2 = RayInfo%PhiangOutSvIdx(i,2)
  IF(k1 .NE. 1) Buf0(:, k1) = PhiAngIn1g(:, k1)
  IF(k2 .NE. 1) Buf0(:, k1) = PhiAngIn1g(:, k2)  
ENDDO
CALL REDUCE(Buf0(:,:), Buf1(:,:), n1, n2, PE%MPI_RT_COMM, FALSE)
IF(PE%RTMASTER) THEN
  DO i = 1, RayInfo%nRotRay
    DO j = 1, 2
      k1 = RayInfo%PhiAngInSvIdx(i, j);   k2 = RayInfo%PhiangOutSvIdx(i, j)
      IF(k1 .NE.1) PhiAngIn1g(:,k1) = buf1(:,k1)
      IF(k2 .NE.1) PhiAngIn1g(:,k2) = buf1(:,k2)
    ENDDO
  ENDDO
ENDIF
DEALLOCATE(Lists, buf0, buf1)
#endif
!CALL MPI_SYNC(PE%MPI_RT_COMM)
END SUBROUTINE