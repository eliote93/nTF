SUBROUTINE initRT(RayInfo, CoreInfo, nTracerCntl, PE)

USE OMP_LIB
USE allocs
USE PARAM,   ONLY : mesg, TRUE
USE TYPEDEF, ONLY : RayInfo_Type, CoreInfo_type, PE_TYPE, AziAngleInfo_Type, PolarAngle_Type
USE Cntl,    ONLY : nTracerCntl_Type
USE MOC_MOD, ONLY : trackingdat, ApproxExp, nMaxRaySeg, nMaxCellRay, nMaxAsyRay, nMaxCoreRay, wtang, wtsurf, Comp, mwt, mwt2, SrcAng1, SrcAng2, SrcAngnm1, SrcAngnm2, EXPA, EXPB, hwt, &
                    AziMap, DcmpAsyClr, nOmpAng
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
      CALL dmalloc(TrackingDat(ithr)%Phis, nFsr)
      CALL dmalloc(TrackingDat(ithr)%Jout, 3, nbd, nxy)
    END DO
  END IF
  
  DO ithr = 1, nThr
    CALL dmalloc(TrackingDat(ithr)%FsrIdx,                     nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%ExpAppIdx,                  nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%OptLenList,                 nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%ExpAppPolar,    nPolarAng,  nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%PhiAngOutPolar, nPolarAng,  nMaxRaySeg + 2)
    
    TrackingDat(ithr)%lAlloc = TRUE
  END DO
ELSE IF (.NOT. ldcmp) THEN
  DO ithr = 1, nThr
    CALL dmalloc(TrackingDat(ithr)%phisnm,    ng, nFsr)
    CALL dmalloc(TrackingDat(ithr)%Joutnm, 3, ng, nbd, nxy)
    
    TrackingDat(ithr)%lAllocNM = TRUE
  END DO
END IF

! AFSS
IF (nTracerCntl%lAFSS) THEN
  IF (MOD(nAziAng / 2, nthr) .NE. 0) CALL terminate('WRONG_MOC_TRD')
  
  nOmpAng = nAziAng / nthr
  
  DO ithr = 1, nthr
    !CALL dmalloc(TrackingDat(ithr)%phi1a, nPolarAng, nFsr, nOmpAng)
    !CALL dmalloc(TrackingDat(ithr)%phi2a, nPolarAng, nFsr, nOmpAng)
    
    CALL dmalloc(TrackingDat(ithr)%phi1a, nPolarAng, nFsr, nAziAng)
    CALL dmalloc(TrackingDat(ithr)%phi2a, nPolarAng, nFsr, nAziAng)
  END DO
  
  !CALL initAFSS(RayInfo, CoreInfo, nTracerCntl, PE)
END IF

! Dcmp. & Lin Src CASMO
IF (ldcmp .AND. lLinSrcCASMO) THEN
  DO ithr = 1, nThr
    CALL dmalloc(TrackingDat(ithr)%FsrIdx,                      nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%ExpAppIdxnm,             ng, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%OptLenListnm,            ng, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%ExpAppnm,     nPolarAng, ng, nMaxRaySeg, nMaxCoreRay)
    CALL dmalloc(TrackingDat(ithr)%PhiAngOutnm,  nPolarAng, ng, nMaxRaySeg + 2)
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
    CALL dmalloc(SrcAng1, nPolarAng, nFsr, nAziAng)
    CALL dmalloc(SrcAng2, nPolarAng, nFsr, nAziAng)
      
    DO ithr = 1, nThr
      TrackingDat(ithr)%SrcAng1 => SrcAng1
      TrackingDat(ithr)%SrcAng2 => SrcAng2
      
      CALL dmalloc(TrackingDat(ithr)%phim, nod, nFsr)
    END DO
  ELSE
    CALL dmalloc(SrcAngnm1, ng, nPolarAng, nFsr, nAziAng)
    CALL dmalloc(SrcAngnm2, ng, nPolarAng, nFsr, nAziAng)
    
    DO ithr = 1, nThr
      TrackingDat(ithr)%SrcAngnm1 => SrcAngnm1
      TrackingDat(ithr)%SrcAngnm2 => SrcAngnm2
      
      CALL dmalloc(TrackingDat(ithr)%phimnm, nod, ng, nFsr)
    END DO
  END IF
END IF

! P1 : Dcmp.
IF (ldcmp) THEN
  DO iAzi = 1, nAziAng / 2
    AziMap(iAzi, 1) = 1
    AziMap(iAzi, 2) = 2
    
    AziMap(nAziAng - iAzi + 1, 1) = 3
    AziMap(nAziAng - iAzi + 1, 2) = 4
  END DO
  
  DO ithr = 1, nThr
    TrackingDat(ithr)%AziMap => AziMap
  END DO
END IF
! ----------------------------------------------------
1000 CONTINUE

NULLIFY (AziAng)
NULLIFY (PolarAng)
! ----------------------------------------------------

END SUBROUTINE initRT
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE InitAFSS(RayInfo, CoreInfo, PE, nTracerCntl)

USE allocs
USE TYPEDEF, ONLY : RayInfo_Type, CoreInfo_Type, PE_TYPE, RotRayInfo_Type
USE CNTL,    ONLY : nTracerCntl_Type
USE IOUTIL,  ONLY : terminate
USE moc_mod, ONLY : OmpRayBeg, OmpRayEnd, OmpRayBegBd, OmpRayEndBd, OmpMap, nOmpAng, trackingdat

IMPLICIT NONE

TYPE (RayInfo_Type)     :: RayInfo
TYPE (CoreInfo_Type)    :: CoreInfo
TYPE (nTracerCntl_Type) :: nTracerCntl
TYPE (PE_TYPE)          :: PE
! ----------------------------------------------------
INTEGER :: nAziAng, nPolarAng, nRotray, nFSR, nthr, OmpTmp, ithr, nAzi, iRotRay, iazi, iRay, jazi

TYPE(RotRayInfo_Type), POINTER :: LocRotRay, RotRay(:)

INTEGER, POINTER, DIMENSION(:) :: nAziRotRay, nAziRotRay0
! ----------------------------------------------------

nAziAng   = RayInfo%nAziAngle
nPolarAng = RayInfo%nPolarAngle
nRotRay   = RayInfo%nRotRay
RotRay   => RayInfo%RotRay

nFSR = CoreInfo%nCoreFsr

nthr = PE%nThread
! ----------------------------------------------------
IF (nTracerCntl%lScatBd) THEN
  ! Rot. Ray for Azm. Ang.
  CALL dmalloc(nAziRotRay, nAziAng/2) ! # of Rot. Rays at [0 ~ pi]
  
  DO iRotRay = 1, nRotRay
    iazi = RotRay(iRotRay)%Ang1
    
    nAziRotRay(iazi) = nAziRotRay(iazi) + 1
  END DO
  
  ! Rot. Ray for Thr.
  CALL dmalloc(nAziRotRay0,   nthr)
  CALL dmalloc0(OmpRayBeg, 0, nthr)
  CALL dmalloc0(OmpRayEnd, 0, nthr)
  
  OmpRayBeg(0) = 1
  OmpRayEnd(0) = 0
  
  nAzi   = nAziAng / 2 / nthr
  OmpTmp = 0
  jazi   = 0
  
  DO ithr = 1, nthr
    OmpRayBeg(ithr) = OmpRayBeg(ithr-1) + OmpTmp
    
    DO iazi = 1, nAzi
      jazi = jazi + 1
      
      nAziRotRay0(ithr) = nAziRotRay0(ithr) + nAziRotRay(jazi)
    END DO
    
    OmpTmp = nAziRotRay0(ithr)
    OmpRayEnd(ithr) = OmpRayEnd(ithr-1) + OmpTmp
  END DO
  
  ! Parallel Ang.
  DO iRotRay = 1, nRotRay
    LocRotRay => RotRay(iRotRay)
    
    LocRotRay%OmpAng1 = MOD(LocRotRay%Ang1, nAzi)
    
    IF (LocRotRay%OmpAng1 .EQ. 0) LocRotRay%OmpAng1 = nAzi
    
    LocRotRay%OmpAng2 = LocRotRay%OmpAng1 + nAzi
  END DO
  
  DO iRotRay = 1, nRotRay
    LocRotRay => RotRay(iRotRay)
    
    DO iRay = 1, LocRotRay%nRay, 2
      LocRotRay%OmpRayIdx(iRay)   = LocRotRay%OmpAng1
      LocRotRay%OmpRayIdx(iRay+1) = LocRotRay%OmpAng2
    END DO
  END DO
  
  CALL dmalloc(OmpMap, nthr, nAzi*2)
  
  DO ithr = 1, nthr
    DO iazi = 1, nAzi
      OmpTmp = OmpTmp + 1
      
      OmpMap(ithr, iazi)      = OmpTmp
      OmpMap(ithr, iazi+nAzi) = nAziAng - OmpTmp + 1
    END DO
  END DO
! ----------------------------------------------------
ELSE
  ! Rot. Ray for Azm. Ang.
  CALL dmalloc(nAziRotRay, nAziAng) ! # of Rot. Rays at [0 ~ 2*pi]
  
  DO iRotRay = 1, nRotRay
    iazi = RotRay(iRotRay)%Ang1
    
    nAziRotRay(iazi) = nAziRotRay(iazi) + 1
  END DO
  
  ! Rot. Ray for Thr. with Forward Dir.
  CALL dmalloc(nAziRotRay0,   2*nthr)
    
  CALL dmalloc0(OmpRayBegBd, 1, 2, 0, nthr)
  CALL dmalloc0(OmpRayEndBd, 1, 2, 0, nthr)
  
  OmpRayBegBd(1, 0) = 1 ! Range of Rot. Ray for each Thr.
  OmpRayEndBd(1, 0) = 0
  
  nAzi   = nAziAng / nthr / 2 ! # of Azm Ang. for 1 Thr.
  OmpTmp = 0
  jazi   = 0
  
  DO ithr = 1, nthr
    OmpRayBegBd(1, ithr) = OmpRayBegBd(1, ithr-1) + OmpTmp
    
    DO iazi = 1, nAzi
      jazi = jazi + 1
      
      nAziRotRay0(ithr) = nAziRotRay0(ithr) + nAziRotRay(jazi) ! # of Rot. Rays for each Thr.
    END DO
    
    OmpTmp = nAziRotRay0(ithr)
    
    OmpRayEndBd(1, ithr) = OmpRayEndBd(1, ithr-1) + OmpTmp
  END DO
  
  ! Rot. Ray for Thr. with Backward Dir.
  nAziRotRay0 = 0
  jazi        = nAziAng + 1
  OmpTmp      = 0
  
  OmpRayBegBd(2, 0) = nRotRay + 1
  OmpRayEndBd(2, 0) = nRotRay
  
  DO ithr = 1, nthr
    OmpRayEndBd(2, ithr) = OmpRayEndBd(2, ithr-1) - OmpTmp
    
    DO iazi = 1, nAzi
      jazi = jazi - 1
      
      nAziRotRay0(ithr) = nAziRotRay0(ithr) + nAziRotRay(jazi)
    END DO
    
    OmpTmp = nAziRotRay0(ithr)
    
    OmpRayBegBd(2, ithr) = OmpRayBegBd(2, ithr-1) - OmpTmp
  END DO
  
  ! Parallel Ang. for Core Rays in Rot. Ray
  DO iRotRay = 1, nRotRay
    LocRotRay => RotRay(iRotRay)
    
    IF (LocRotRay%nRay .GT. 1) THEN ! Including Reflection
      IF (LocRotRay%Ang1 .LE. nAziAng/2) THEN
        LocRotRay%OmpAng1 = MOD(LocRotRay%Ang1, nAzi)
        
        IF (LocRotRay%OmpAng1 .EQ. 0) LocRotRay%OmpAng1 = nAzi
        
        LocRotRay%OmpAng2 = LocRotRay%OmpAng1 + nAzi
      ELSE IF (LocRotRay%Ang1 .GT. nAziAng/2) THEN
        LocRotRay%OmpAng1 = MOD(nAziAng+1 - LocRotRay%Ang1, nAzi) + nAzi
        
        IF (LocRotRay%OmpAng1 .EQ. nAzi) LocRotRay%OmpAng1 = LocRotRay%OmpAng1 + nAzi
        
        LocRotRay%OmpAng2 = LocRotRay%OmpAng1 - nAzi
      END IF
    ELSE IF (LocRotRay%nRay .EQ. 1) THEN ! No Reflection
      IF (LocRotRay%Ang1 .LE. nAziAng/2) THEN
        LocRotRay%OmpAng1 = MOD(LocRotRay%Ang1, nAzi)
        
        IF (LocRotRay%OmpAng1 .EQ. 0) LocRotRay%OmpAng1 = nAzi
        
        LocRotRay%OmpAng2 = 0
      ELSE IF (LocRotRay%Ang1 .GT. nAziAng/2) THEN
        
        LocRotRay%OmpAng1 = MOD(nAziAng+1-LocRotRay%Ang1, nAzi) + nAzi
        
        IF (LocRotRay%OmpAng1 .EQ. nAzi) LocRotRay%OmpAng1 = LocRotRay%OmpAng1 + nAzi
        
        LocRotRay%OmpAng2 = 0
      END IF
    END IF
  END DO
  
  DO iRotRay = 1, nRotRay
    LocRotRay => RotRay(iRotRay)
    
    IF (LocRotRay%nRay .GT. 1) THEN ! Including Reflection
      IF (MOD(LocRotRay%nRay, 2) .EQ. 0) THEN
        DO iRay = 1, LocRotRay%nRay, 2 ! NOTICE : increase with 2
          LocRotRay%OmpRayIdx(iRay)   = LocRotRay%OmpAng1
          LocRotRay%OmpRayIdx(iRay+1) = LocRotRay%OmpAng2 ! Azm. Ang. after Reflection ??
        END DO
      ELSE
        DO iRay = 1, LocRotRay%nRay-1, 2
          LocRotRay%OmpRayIdx(iRay)   = LocRotRay%OmpAng1
          LocRotRay%OmpRayIdx(iRay+1) = LocRotRay%OmpAng2
        END DO
        
        LocRotRay%OmpRayIdx(LocRotRay%nRay) = LocRotRay%OmpAng1
      END IF
    ELSE ! No Reflection
      LocRotRay%OmpRayIdx(1) = LocRotRay%OmpAng1
    END IF
  END DO
  
  CALL dmalloc(OmpMap, nthr, nAzi*2)
  
  DO ithr = 1, nthr
    DO iazi = 1, nAzi
      OmpTmp = OmpTmp + 1
      
      OmpMap(ithr, iazi)      = OmpTmp ! Numeric # of Azm. Ang
      OmpMap(ithr, iazi+nAzi) = nAziAng - OmpTmp + 1 ! Reflcetion to x-axis
    END DO
  END DO
END IF
! ----------------------------------------------------
NULLIFY (LocRotRay)
NULLIFY (RotRay)

DEALLOCATE (nAziRotRay)
DEALLOCATE (nAziRotRay0)
! ----------------------------------------------------

END SUBROUTINE InitAFSS
! ------------------------------------------------------------------------------------------------------------