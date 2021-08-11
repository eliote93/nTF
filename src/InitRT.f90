SUBROUTINE initRT(RayInfo, CoreInfo, nTracerCntl, PE)

USE OMP_LIB
USE allocs
USE PARAM,   ONLY : mesg, TRUE
USE TYPEDEF, ONLY : RayInfo_Type, CoreInfo_type, PE_TYPE, AziAngleInfo_Type, PolarAngle_Type
USE Cntl,    ONLY : nTracerCntl_Type
USE MOC_MOD, ONLY : trackingdat, ApproxExp, nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay, wtang, wtsurf, Comp, mwt, mwt2, SrcAng1g1, SrcAng1g2, SrcAngNg1, SrcAngNg2, EXPA, EXPB, hwt, DcmpAsyClr
USE geom,    ONLY : nbd, ng
USE files,   ONLY : io8
USE ioutil,  ONLY : message, terminate

IMPLICIT NONE

TYPE (RayInfo_Type)     :: RayInfo
TYPE (CoreInfo_Type)    :: CoreInfo
TYPE (nTracerCntl_Type) :: nTracerCntl
TYPE (PE_TYPE)          :: PE

INTEGER :: nFsr, nxy, ithr, scatod, nod, nPolarAng, nAziAng, ipol, iazi, nthr
REAL :: wttmp, wtsin2, wtcos, wtpolar
LOGICAL :: lscat1, ldcmp, lLinSrcCASMO

TYPE (AziAngleInfo_Type), POINTER, DIMENSION(:) :: AziAng
TYPE (PolarAngle_Type),   POINTER, DIMENSION(:) :: PolarAng
! ----------------------------------------------------

mesg = 'Allocating RT Variables...'
IF (PE%master) CALL message(io8, TRUE, TRUE, mesg)

AziAng   => RayInfo%AziAngle
PolarAng => RayInfo%PolarAngle
nPolarAng = RayInfo%nPolarAngle
nAziAng   = RayInfo%nAziAngle

nthr = PE%nThread
CALL omp_set_num_threads(nThr)

nFsr = CoreInfo%nCoreFsr
nxy  = CoreInfo%nxy

scatod       = nTracerCntl%scatod
lscat1       = nTracerCntl%lscat1
ldcmp        = nTracerCntl%lDomainDcmp
lLinSrcCASMO = nTracerCntl%lLinSrcCASMO

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
    
    wtsurf(ipol, iazi, 1) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv) ! S
    wtsurf(ipol, iazi, 3) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv) ! N
    wtsurf(ipol, iazi, 2) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv) ! W
    wtsurf(ipol, iazi, 4) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv) ! E
  END DO
END DO

! Hex.
IF (nTracerCntl%lHex) THEN
  CALL dmalloc(hwt, nPolarAng, nAziAng)
  
  DO ipol = 1, nPolarAng
    DO iazi = 1, nAziAng
      hwt(ipol, iazi) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del
    END DO
  END DO
END IF
! ----------------------------------------------------
! Basic
DO ithr = 1, nThr
  TrackingDat(ithr)%ExpA   => ExpA
  TrackingDat(ithr)%ExpB   => ExpB
  TrackingDat(ithr)%wtang  => wtang
  TrackingDat(ithr)%wtsurf => wtsurf
  
  IF (.NOT. nTracerCntl%lHex) CYCLE
  
  TrackingDat(ithr)%hwt => hwt
END DO

! Basic : GM vs. NM
IF (.NOT. nTracerCntl%lNodeMajor) THEN
  IF (.NOT. ldcmp) THEN
    DO ithr = 1, nThr
      CALL dmalloc(TrackingDat(ithr)%phis1g, nFsr)
      CALL dmalloc(TrackingDat(ithr)%Jout1g, 3, nbd, nxy)
    END DO
  END IF
  
  DO ithr = 1, nThr
    CALL dmalloc(TrackingDat(ithr)%FsrIdx,                  nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%ExpAppIdx,               nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%OptLenList,              nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%ExpAppPolar, nPolarAng,  nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%PhiAngOut1g, nPolarAng,  nMaxRaySeg + 2)
    
    TrackingDat(ithr)%lAlloc = TRUE
  END DO
ELSE IF (.NOT. ldcmp) THEN
  DO ithr = 1, nThr
    CALL dmalloc(TrackingDat(ithr)%phisNg,    ng, nFsr)
    CALL dmalloc(TrackingDat(ithr)%JoutNg, 3, ng, nbd, nxy)
    
    TrackingDat(ithr)%lAllocNM = TRUE
  END DO
END IF

! AFSS
IF (nTracerCntl%lAFSS) THEN
  DO ithr = 1, nthr
    IF (.NOT. nTracerCntl%lNodeMajor) THEN
      CALL dmalloc(TrackingDat(ithr)%phia1g, 2,     nPolarAng, nAziAng, nFsr)
    ELSE
      CALL dmalloc(TrackingDat(ithr)%phiaNg, 2, ng, nPolarAng, nAziAng, nFsr)
    END IF
  END DO
END IF

! Dcmp. & Lin Src CASMO
IF (ldcmp .AND. lLinSrcCASMO) THEN
  DO ithr = 1, nThr
    CALL dmalloc(TrackingDat(ithr)%FsrIdx,                      nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%ExpAppIdxnm,             ng, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%OptLenListNg,            ng, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%PhiAngOutNg,  nPolarAng, ng, nMaxRaySeg + 2)
    CALL dmalloc(TrackingDat(ithr)%cmOptLen,     nPolarAng, ng, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%cmOptLenInv,  nPolarAng, ng, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%q0,                      ng, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%q1,           nPolarAng, ng, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%E1,           nPolarAng, ng, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%E3,           nPolarAng, ng, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%R1,           nPolarAng, ng, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%R3,           nPolarAng, ng, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%x0,                       2, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%y0,                       2, nMaxRaySeg, nMaxCoreRay)
    
    TrackingDat(ithr)%lAllocLinSrc = TRUE
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
IF (.NOT. ldcmp) THEN
  IF (.NOT. nTracerCntl%lNodeMajor) THEN
    CALL dmalloc(SrcAng1g1, nPolarAng, nFsr, nAziAng)
    CALL dmalloc(SrcAng1g2, nPolarAng, nFsr, nAziAng)
      
    DO ithr = 1, nThr
      TrackingDat(ithr)%SrcAng1g1 => SrcAng1g1
      TrackingDat(ithr)%SrcAng1g2 => SrcAng1g2
      
      CALL dmalloc(TrackingDat(ithr)%phim1g, nod, nFsr)
    END DO
  ELSE
    CALL dmalloc(SrcAngNg1, ng, nPolarAng, nFsr, nAziAng)
    CALL dmalloc(SrcAngNg2, ng, nPolarAng, nFsr, nAziAng)
    
    DO ithr = 1, nThr
      TrackingDat(ithr)%SrcAngNg1 => SrcAngNg1
      TrackingDat(ithr)%SrcAngNg2 => SrcAngNg2
      
      CALL dmalloc(TrackingDat(ithr)%phimNg, nod, ng, nFsr)
    END DO
  END IF
END IF
! ----------------------------------------------------
1000 CONTINUE

NULLIFY (AziAng)
NULLIFY (PolarAng)
! ----------------------------------------------------

END SUBROUTINE initRT