SUBROUTINE initRT(RayInfo, CoreInfo, nTracerCntl, PE)

USE OMP_LIB
USE allocs
USE PARAM,   ONLY : mesg, TRUE
USE TYPEDEF, ONLY : RayInfo_Type, CoreInfo_type, PE_TYPE, AziAngleInfo_Type, PolarAngle_Type
USE Cntl,    ONLY : nTracerCntl_Type
USE MOC_MOD, ONLY : ApproxExp, nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay, wtang, wtsurf, Comp, mwt, mwt2, EXPA, EXPB, hwt, DcmpAsyClr, AziRotRay, OmpRotRay, OmpAzi, &
                    trackingdat, SrcAng1g1, SrcAng1g2, SrcAngNg1, SrcAngNg2, phia1g1, phia1g2
USE geom,    ONLY : nbd, ng
USE files,   ONLY : io8
USE ioutil,  ONLY : message, terminate
USE HexData, ONLY : hLgc, hRotRay, hcRay

IMPLICIT NONE

TYPE (RayInfo_Type)     :: RayInfo
TYPE (CoreInfo_Type)    :: CoreInfo
TYPE (nTracerCntl_Type) :: nTracerCntl
TYPE (PE_TYPE)          :: PE

INTEGER :: nFsr, nxy, ithr, scatod, nod, nPol, nAzi, ipol, iazi, nthr, iray, ntmp1, ntmp2, itmp, nRotRay
REAL :: wttmp, wtsin2, wtcos, wtpolar
LOGICAL :: lscat1, ldcmp, lLinSrcCASMO, lGM, lHex

TYPE (AziAngleInfo_Type), POINTER, DIMENSION(:) :: AziAng
TYPE (PolarAngle_Type),   POINTER, DIMENSION(:) :: PolarAng
! ----------------------------------------------------

mesg = 'Allocating RT Variables...'
IF (PE%master) CALL message(io8, TRUE, TRUE, mesg)

AziAng   => RayInfo%AziAngle
PolarAng => RayInfo%PolarAngle
nPol      = RayInfo%nPolarAngle
nAzi      = RayInfo%nAziAngle
nRotRay   = RayInfo%nRotRay

nthr = PE%nThread
CALL omp_set_num_threads(nThr)

nFsr = CoreInfo%nCoreFsr
nxy  = CoreInfo%nxy

lGM          = .NOT. nTracerCntl%lNodeMajor
lHex         = nTracerCntl%lHex
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
CALL ApproxExp(PolarAng, nPol)

! Wgt.
CALL dmalloc(wtang,  nPol, nAzi)
CALL dmalloc(wtsurf, nPol, nAzi, 4)

DO ipol = 1, nPol
  wttmp = PolarAng(ipol)%weight * PolarAng(ipol)%sinv
  
  DO iazi = 1, nAzi
    wtang(ipol, iazi) = wttmp * AziAng(iazi)%weight * AziAng(iazi)%del
    
    wtsurf(ipol, iazi, 1) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv) ! S
    wtsurf(ipol, iazi, 3) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%sinv) ! N
    wtsurf(ipol, iazi, 2) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv) ! W
    wtsurf(ipol, iazi, 4) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del / ABS(AziAng(iazi)%cosv) ! E
  END DO
END DO

! Hex.
IF (lHex) THEN
  CALL dmalloc(hwt, nPol, nAzi)
  
  DO ipol = 1, nPol
    DO iazi = 1, nAzi
      hwt(ipol, iazi) = PolarAng(ipol)%weight * AziAng(iazi)%weight * AziAng(iazi)%del
    END DO
  END DO
END IF
! ----------------------------------------------------
! Omp Rot Ray, DEBUG
ntmp1 = nRotRay / nthr
ntmp2 = nRotRay - ntmp1 * nthr

iray = 1

DO itmp = 1, ntmp2
  OmpRotRay(1, itmp) = iray
  OmpRotRay(2, itmp) = iray + ntmp1
  
  iray = iray + ntmp1 + 1
END DO

DO itmp = ntmp2 + 1, nthr
  OmpRotRay(1, itmp) = iray
  OmpRotRay(2, itmp) = iray + ntmp1 - 1
  
  iray = iray + ntmp1
END DO

! Omp Azi
ntmp1 = nAzi / nthr
ntmp2 = nAzi - ntmp1 * nthr

OmpAzi(0, 1:ntmp2)      = ntmp1 + 1
OmpAzi(0, ntmp2+1:nthr) = ntmp1

iazi = 0

DO itmp = 1, ntmp1
  DO ithr = 1, nthr
    iazi = iazi + 1
    
    OmpAzi(itmp, ithr) = iazi
  END DO
END DO

DO itmp = 1, ntmp2
  iazi = iazi + 1
  
  OmpAzi(ntmp1+1, itmp) = iazi
END DO
! ----------------------------------------------------
! Basic
DO ithr = 1, nThr
  TrackingDat(ithr)%ExpA   => ExpA
  TrackingDat(ithr)%ExpB   => ExpB
  TrackingDat(ithr)%wtang  => wtang
  TrackingDat(ithr)%wtsurf => wtsurf
  
  IF (.NOT. lHex) CYCLE
  
  TrackingDat(ithr)%hwt => hwt
END DO

! Basic : GM vs. NM
IF (lGM) THEN
  IF (.NOT. ldcmp) THEN
    DO ithr = 1, nThr
      CALL dmalloc(TrackingDat(ithr)%phis1g, nFsr)
      CALL dmalloc(TrackingDat(ithr)%Jout1g, 3, nbd, nxy)
    END DO
  END IF
  
  DO ithr = 1, nThr
    CALL dmalloc(TrackingDat(ithr)%FsrIdx,            nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%ExpAppIdx,         nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%OptLenList,        nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%ExpAppPolar, nPol, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%PhiAngOut1g, nPol, nMaxRaySeg + 2)
    
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
IF (nTracerCntl%lAFSS .AND. .NOT.ldcmp .AND. lGM) THEN
  !CALL dmalloc(phia1g1, nPol, nFsr, nAzi) ! pol 2 fsr
  !CALL dmalloc(phia1g2, nPol, nFsr, nAzi)
  CALL dmalloc(phia1g1, nFsr, nPol, nAzi) ! fsr 2 pol
  CALL dmalloc(phia1g2, nFsr, nPol, nAzi)
  
  IF (hLgc%lNoRef) THEN
    ! SET : Azi Rot Ray
    CALL dmalloc0(AziRotRay, 0, nRotRay, 1, nAzi)
    
    DO iray = 1, nRotRay
      iazi = hcRay(hRotRay(iray)%cRayIdx(1))%AzmIdx
      
      AziRotRay(0, iazi) = AziRotRay(0, iazi) + 1
      
      AziRotRay(AziRotRay(0, iazi), iazi) = iray
    END DO
    
    ! ALLOC
    DO ithr = 1, ithr
      !CALL dmalloc(TrackingDat(ithr)%phia1g1, nPol, nFsr, OmpAzi(0, ithr)) ! pol 2 fsr
      !CALL dmalloc(TrackingDat(ithr)%phia1g2, nPol, nFsr, OmpAzi(0, ithr))
      CALL dmalloc(TrackingDat(ithr)%phia1g1, nFsr, nPol, OmpAzi(0, ithr)) ! fsr 2 pol
      CALL dmalloc(TrackingDat(ithr)%phia1g2, nFsr, nPol, OmpAzi(0, ithr))
    END DO
  ELSE
    DO ithr = 1, nthr
      !CALL dmalloc(TrackingDat(ithr)%phia1g1, nPol, nFsr, nAzi) ! pol 2 fsr
      !CALL dmalloc(TrackingDat(ithr)%phia1g2, nPol, nFsr, nAzi)
      CALL dmalloc(TrackingDat(ithr)%phia1g1, nFsr, nPol, nAzi) ! fsr 2 pol
      CALL dmalloc(TrackingDat(ithr)%phia1g2, nFsr, nPol, nAzi)
    END DO
  END IF
END IF

! Dcmp. & Lin Src CASMO
IF (ldcmp .AND. lLinSrcCASMO) THEN
  DO ithr = 1, nThr
    CALL dmalloc(TrackingDat(ithr)%FsrIdx,                 nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%ExpAppIdxnm,        ng, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%OptLenListNg,       ng, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%PhiAngOutNg,  nPol, ng, nMaxRaySeg + 2)
    CALL dmalloc(TrackingDat(ithr)%cmOptLen,     nPol, ng, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%cmOptLenInv,  nPol, ng, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%q0,                 ng, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%q1,           nPol, ng, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%E1,           nPol, ng, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%E3,           nPol, ng, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%R1,           nPol, ng, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%R3,           nPol, ng, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%x0,                  2, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%y0,                  2, nMaxRaySeg, nMaxCoreRay)
    
    TrackingDat(ithr)%lAllocLinSrc = TRUE
  END DO
END IF
! ----------------------------------------------------
IF (.NOT. lscat1) GO TO 1000

! P1 : mwt, mwt2
CALL dmalloc(Comp, nod, nPol, nAzi)
CALL dmalloc(mwt,  nod, nPol, nAzi)
CALL dmalloc(mwt2, nod, nPol, nAzi)

DO ipol = 1, nPol
  wttmp   = PolarAng(ipol)%sinv
  wtsin2  = PolarAng(ipol)%sinv * PolarAng(ipol)%sinv
  wtcos   = PolarAng(ipol)%cosv
  wtpolar = 1.5_8 * PolarAng(ipol)%cosv * PolarAng(ipol)%cosv - 0.5_8
  
  DO iazi = 1, nAzi
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

DO ipol = 1, nPol
  DO iazi = 1, nAzi
    comp(1:2, ipol, iazi) = comp(1:2, ipol, iazi) * 3._8
    
    IF (scatod .GE. 2) comp(3:5, ipol, iazi) = comp(3:5, ipol, iazi) * 5._8
    IF (scatod .EQ. 3) comp(6:9, ipol, iazi) = comp(6:9, ipol, iazi) * 7._8
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
  IF (lGM) THEN
    CALL dmalloc(SrcAng1g1, nPol, nFsr, nAzi)
    CALL dmalloc(SrcAng1g2, nPol, nFsr, nAzi)
      
    DO ithr = 1, nThr
      TrackingDat(ithr)%SrcAng1g1 => SrcAng1g1
      TrackingDat(ithr)%SrcAng1g2 => SrcAng1g2
      
      CALL dmalloc(TrackingDat(ithr)%phim1g, nod, nFsr)
    END DO
  ELSE
    CALL dmalloc(SrcAngNg1, ng, nPol, nFsr, nAzi)
    CALL dmalloc(SrcAngNg2, ng, nPol, nFsr, nAzi)
    
    DO ithr = 1, nThr
      TrackingDat(ithr)%SrcAngNg1 => SrcAngNg1
      TrackingDat(ithr)%SrcAngNg2 => SrcAngNg2
      
      CALL dmalloc(TrackingDat(ithr)%phimNg, nod, nFsr, ng)
    END DO
  END IF
END IF
! ----------------------------------------------------
1000 CONTINUE

NULLIFY (AziAng)
NULLIFY (PolarAng)
! ----------------------------------------------------

END SUBROUTINE initRT