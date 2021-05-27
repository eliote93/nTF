SUBROUTINE initRT(RayInfo, CoreInfo, nTracerCntl, PE)

USE OMP_LIB
USE allocs
USE PARAM,   ONLY : TRUE
USE TYPEDEF, ONLY : RayInfo_Type, CoreInfo_type, PE_TYPE, AziAngleInfo_Type, PolarAngle_Type
USE Cntl,    ONLY : nTracerCntl_Type
USE MOC_MOD, ONLY : trackingdat, ApproxExp, nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay, wtang, wtsurf, Comp, mwt, mwt2, SrcAng1, SrcAng2, SrcAngnm1, SrcAngnm2, EXPA, EXPB
USE geom,    ONLY : nbd, ng

IMPLICIT NONE

TYPE (RayInfo_Type)     :: RayInfo
TYPE (CoreInfo_Type)    :: CoreInfo
TYPE (nTracerCntl_Type) :: nTracerCntl
TYPE (PE_TYPE)          :: PE

INTEGER :: nFsr, nxy, ithr, scatod, nod, nPolarAng, nAziAng, ipol, iazi, nthr
REAL :: wttmp, wtsin2, wtcos, wtpolar
LOGICAL :: lscat1

TYPE (AziAngleInfo_Type), POINTER, DIMENSION(:) :: AziAng
TYPE (PolarAngle_Type),   POINTER, DIMENSION(:) :: PolarAng
! ----------------------------------------------------

AziAng   => RayInfo%AziAngle
PolarAng => RayInfo%PolarAngle
nPolarAng = RayInfo%nPolarAngle
nAziAng   = RayInfo%nAziAngle

CALL omp_set_num_threads(PE%nThread)
nthr = PE%nThread

nFsr = CoreInfo%nCoreFsr
nxy  = CoreInfo%nxy

scatod = nTracerCntl%scatod
lscat1 = nTracerCntl%lscat1

IF (nTracerCntl%lscat1) THEN
  SELECT CASE (scatod)
    CASE (1); nod = 2
    CASE (2); nod = 5
    CASE (3); nod = 9
  END SELECT
END IF
! ----------------------------------------------------
! EXP
CALL ApproxExp(PolarAng, nPolarAng)

! Wgt.
CALL dmalloc(wtang,  nPolarAng, nAziAng)
CALL dmalloc(wtsurf, nPolarAng, nAziAng, 4)

DO ipol = 1, nPolarAng
  wttmp = PolarAng(ipol)%weight * PolarAng(ipol)%sinv
  
  DO iazi = 1, nAziAng
    wtang(ipol, iazi) = wttmp * AziAng(iazi)%weight * AziAng(iazi)%del
    
    wtsurf(ipol, iazi, 1) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv)
    wtsurf(ipol, iazi, 3) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv)
    wtsurf(ipol, iazi, 2) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv)
    wtsurf(ipol, iazi, 4) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv)
  END DO
END DO
! ----------------------------------------------------
! Basic
DO ithr = 1, nThr
  CALL dmalloc(TrackingDat(ithr)%FsrIdx,         nMaxRaySeg, nMaxCoreRay)
  CALL dmalloc(TrackingDat(ithr)%ExpAppIdx,      nMaxRaySeg, nMaxCoreRay)
  CALL dmalloc(TrackingDat(ithr)%OptLenList,     nMaxRaySeg, nMaxCoreRay)
  CALL dmalloc(TrackingDat(ithr)%ExpAppPolar,    nPolarAng,  nMaxRaySeg, nMaxCoreRay)
  CALL dmalloc(TrackingDat(ithr)%PhiAngOutPolar, nPolarAng,  nMaxRaySeg + 2)
  
  TrackingDat(ithr)%Expa   => Expa
  TrackingDat(ithr)%Expb   => Expb
  TrackingDat(ithr)%wtang  => wtang
  TrackingDat(ithr)%wtsurf => wtsurf
END DO

! Basic : GM vs. NM
IF (.NOT. nTracerCntl%lNodeMajor) THEN
  DO ithr = 1, nThr
    CALL dmalloc(TrackingDat(ithr)%Phis, nFsr)
    CALL dmalloc(TrackingDat(ithr)%Jout, 3, nbd, nxy)
    
    TrackingDat(ithr)%lAlloc = TRUE
  END DO
ELSE
  DO ithr = 1, nThr
    CALL dmalloc(TrackingDat(ithr)%phisnm,    ng, nFsr)
    CALL dmalloc(TrackingDat(ithr)%Joutnm, 3, ng, nbd, nxy)
    
    TrackingDat(ithr)%lAllocNM = TRUE
  END DO
END IF
! ----------------------------------------------------
IF (.NOT. lscat1) GO TO 1000

! P1 : mwt, mwt2
CALL dmalloc(Comp, nod, nPolarAng, nAziAng)
CALL dmalloc(mwt,  nod, nPolarAng, nAziAng)
CALL dmalloc(mwt2, nod, nPolarAng, nAziAng)

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
    Comp(4, ipol, iazi) = wtsin2 * (1._8 - 2._8 * AziAng(iazi)%sinv * AziAng(iazi)%sinv)
    Comp(5, ipol, iazi) = wtsin2 * (       2._8 * AziAng(iazi)%sinv * AziAng(iazi)%cosv)
    
    mwt (3,   ipol, iazi) = comp(3, ipol, iazi) *  wtang(ipol, iazi)
    mwt (4:5, ipol, iazi) = 0.75_8 * comp(4:5, ipol, iazi) *  wtang(ipol, iazi)
    mwt2(3:5, ipol, iazi) = mwt(3:5, ipol, iazi)
    
    IF (scatod .LT. 3) CYCLE
    
    Comp(6, ipol, iazi) = (5._8 * PolarAng(ipol)%cosv * PolarAng(ipol)%cosv - 1._8) * wttmp * AziAng(iazi)%cosv
    Comp(7, ipol, iazi) = (5._8 * PolarAng(ipol)%cosv * PolarAng(ipol)%cosv - 1._8) * wttmp * AziAng(iazi)%sinv
    Comp(8, ipol, iazi) = (wttmp ** 3._8) * ( 4._8 * (AziAng(iazi)%cosv ** 3._8) - 3._8 * AziAng(iazi)%cosv)
    Comp(9, ipol, iazi) = (wttmp ** 3._8) * (-4._8 * (AziAng(iazi)%sinv ** 3._8) + 3._8 * AziAng(iazi)%sinv)
    
    mwt (6:7, ipol, iazi) = 0.375_8 * comp(6:7, ipol, iazi) * wtang(ipol, iazi)
    mwt (8:9, ipol, iazi) = 0.625_8 * comp(8:9, ipol, iazi) * wtang(ipol, iazi)
    mwt2(6:9, ipol, iazi) = -mwt(6:9, ipol, iazi)
  END DO
END DO

DO ithr = 1, nThr
  TrackingDat(ithr)%mwt  => mwt
  TrackingDat(ithr)%mwt2 => mwt2
  
  TrackingDat(ithr)%lAllocP1 = TRUE
END DO
! ----------------------------------------------------
! P1 : GM vs. NM
IF (.NOT. nTracerCntl%lNodeMajor) THEN
  DO ithr = 1, nThr
    CALL dmalloc(SrcAng1, nPolarAng, nFsr, nAziAng)
    CALL dmalloc(SrcAng2, nPolarAng, nFsr, nAziAng)
    
    TrackingDat(ithr)%SrcAng1  => SrcAng1
    TrackingDat(ithr)%SrcAng2  => SrcAng2
    
    CALL dmalloc(TrackingDat(ithr)%phim, nod, nFsr)
  END DO
! NM
ELSE
  DO ithr = 1, nThr
    CALL dmalloc(SrcAngnm1, ng, nPolarAng, nFsr, nAziAng)
    CALL dmalloc(SrcAngnm2, ng, nPolarAng, nFsr, nAziAng)
        
    TrackingDat(ithr)%SrcAngnm1  => SrcAngnm1
    TrackingDat(ithr)%SrcAngnm2  => SrcAngnm2
    
    CALL dmalloc(TrackingDat(ithr)%phimnm, nod, ng, nFsr)
  END DO
END IF
! ----------------------------------------------------
1000 CONTINUE

NULLIFY (AziAng)
NULLIFY (PolarAng)
! ----------------------------------------------------

END SUBROUTINE initRT