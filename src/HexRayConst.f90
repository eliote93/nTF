MODULE HexRayConst

USE ioutil, ONLY : terminate

IMPLICIT NONE

CONTAINS

! ------------------------------------------------------------------------------------------------------------
!                                     01. HEX SET : Core Ray
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetCoreRay()

USE PARAM,   ONLY : TRUE, FALSE, ZERO
USE geom,    ONLY : nZ
USE Moc_Mod, ONLY : nMaxCellRay, nMaxRaySeg
USE HexType, ONLY : Type_HexCoreRay, Type_HexAsyRay
USE HexData, ONLY : nhAsy, ncRay, Asy2Dto1DMap, hAsy, haRay, hcRay, hAsyTypInfo, &
                    nAzmAng, NumMray, AngMray, hmRay, hPinInfo, RodPin, GapPin, hCel, gCel, gCelBss

IMPLICIT NONE

INTEGER :: iAng, iAsy, imRay, jmRay, icRay, iz, iAsyTyp, iGeoTyp, icBss, jcBss, ihPin, jhPin, iCel
INTEGER :: iAsy_Loc, imRay_Loc, iaX_Loc, iaY_Loc
INTEGER :: iAsy_Nxt, imRay_Nxt, iaX_Nxt, iaY_Nxt
INTEGER :: iDir, tDir, iNxt, nSeg, tNumCel, tNumSeg, icRayNxt(-1:1)

INTEGER, PARAMETER :: nMaxModRay = 1000 ! ARBITRARY

INTEGER :: cRayLst(2, -nMaxModRay:nMaxModRay) ! (mRay Idx/Asy Idx, imRay)

LOGICAL, POINTER :: mRayLst(:, :) ! (imRay, iAsy)

TYPE(Type_HexCoreRay), POINTER :: hcRay_Loc(:) ! (icRay)
TYPE(Type_HexAsyRay),  POINTER :: haRay_Loc
! ----------------------------------------------------

ALLOCATE (mRayLst(NumMRay(0), nhAsy))
! ----------------------------------------------------
!               01. SET : mRay Lst
! ----------------------------------------------------
mRayLst = FALSE

DO iAsy = 1, nhAsy
  iGeoTyp = hAsy(iAsy)%GeoTyp
  
  DO imRay = 1, NumMRay(0)
    IF (haRay(iGeoTyp, 1, imRay)%nhpRay < 1) CYCLE
    
    mRayLst(imRay, iAsy) = TRUE
  END DO
END DO
! ----------------------------------------------------
!               02. FIND : Nxt mRay
! ----------------------------------------------------
ALLOCATE(hcRay_Loc (nhAsy * NumMRay(0)))

ncRay = 0

DO iAng = 1, nAzmAng
  DO iAsy = 1, nhAsy
    DO imRay = AngMray(1, iAng), AngMray(2, iAng)
      IF(.NOT. mRayLst(imRay, iAsy)) CYCLE
      
      ncRay = ncRay + 1
      
      mRayLst(imRay, iAsy) = FALSE
      ! ----------------------------
      !      1. GO : Nxt mRay
      ! ----------------------------
      cRayLst  = 0
      icRayNxt = 0
      
      cRayLst(1, 0) = imRay
      cRayLst(2, 0) = iAsy
      
      DO tDir = 1, 2
        iDir = 2*tDir - 3 ! Negative : y¢Ù, Positive : y¢Ö
        
        imRay_Loc = imRay
        iAsy_Loc  = iAsy
        iaX_Loc   = hAsy(iAsy)%iaX
        iaY_Loc   = hAsy(iAsy)%iaY
        
        DO iNxt = 1, nMaxModRay
          iaX_Nxt   = hmRay(imRay_Loc)%NxtAsy_Mov(1, iDir) + iaX_Loc
          iaY_Nxt   = hmRay(imRay_Loc)%NxtAsy_Mov(2, iDir) + iaY_Loc
          imRay_Nxt = hmRay(imRay_Loc)%NxtMray_Mov(iDir)
          iAsy_Nxt  = Asy2Dto1DMap(iaX_Nxt, iaY_Nxt)
          
          IF(iAsy_Nxt .EQ. 0) EXIT
          IF(.NOT. mRayLst(imRay_Nxt, iAsy_Nxt)) EXIT
          
          mRayLst(imRay_Nxt, iAsy_Nxt) = FALSE
          
          imRay_Loc = imRay_Nxt
          iAsy_Loc  = iAsy_Nxt
          iaX_Loc   = iaX_Nxt
          iaY_Loc   = iaY_Nxt
          
          icRayNxt(iDir) = icRayNxt(iDir) + iDir
          
          cRayLst(1, icRayNxt(iDir)) = imRay_Loc
          cRayLst(2, icRayNxt(iDir)) = iAsy_Loc
        END DO
        
        IF (iNxt > nMaxModRay) CALL terminate("SET CORE RAY")
      END DO
      ! ----------------------------
      !      2. COLLECT
      ! ----------------------------
      hcRay_Loc(ncRay)%nmRay  = icRayNxt(1) - icRayNxt(-1) + 1
      hcRay_Loc(ncRay)%AzmIdx = iAng
      
      ALLOCATE(hcRay_Loc(ncRay)%mRayIdx (hcRay_Loc(ncRay)%nmRay))
      ALLOCATE(hcRay_Loc(ncRay)% AsyIdx (hcRay_Loc(ncRay)%nmRay))
      
      jmRay = 0
      
      DO iNxt = icRayNxt(-1), icRayNxt(1)
        jmRay = jmRay + 1
        
        hcRay_Loc(ncRay)%mRayIdx(jmRay) = cRayLst(1, iNxt)
        hcRay_Loc(ncRay)% AsyIdx(jmRay) = cRayLst(2, iNxt)
      END DO
    END DO
  END DO
END DO

ALLOCATE (hcRay(ncRay))
! ----------------------------------------------------
!               03. CP : hcRay
! ----------------------------------------------------
!$OMP PARALLEL PRIVATE(icRay, imRay)
!$OMP DO SCHEDULE(GUIDED)
DO icRay = 1, ncRay
  hcRay(icRay)%nmRay  = hcRay_Loc(icRay)%nmRay
  hcRay(icRay)%AzmIdx = hcRay_Loc(icRay)%AzmIdx
  
  ALLOCATE (hcRay(icRay)%mRayIdx (hcRay(icRay)%nmRay))
  ALLOCATE (hcRay(icRay)% AsyIdx (hcRay(icRay)%nmRay))
  
  DO imRay = 1, hcRay_Loc(icRay)%nmRay
    hcRay(icRay)%mRayIdx(imRay) = hcRay_Loc(icRay)%mRayIdx(imRay)
    hcRay(icRay)% AsyIdx(imRay) = hcRay_Loc(icRay)% AsyIdx(imRay)
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
! ----------------------------------------------------
!               04. SET : Max Cel, Seg
! ----------------------------------------------------
nMaxCellRay = 0
nMaxRaySeg  = 0

!!$OMP PARALLEL PRIVATE(icRay, imRay, jmRay, iAsy, iAsyTyp, iGeoTyp, icBss, haRay_Loc, iCel, ihPin, jhPin, nSeg, iz, jcBss, tNumCel, tNumSeg)
!!$OMP DO REDUCTION(MAX:nMaxCellRay, nMaxRaySeg) SCHEDULE(GUIDED)
DO icRay = 1, ncRay
  tNumCel = 0
  tNumSeg = 0
  
  DO imRay = 1, hcRay_Loc(icRay)%nmRay
    jmRay   = hcRay(icRay)%mRayIdx(imRay)
    iAsy    = hcRay(icRay)%AsyIdx (imRay)
    iAsyTyp = hAsy(iAsy)%AsyTyp
    iGeoTyp = hAsy(iAsy)%GeoTyp
    icBss   = hAsyTypInfo(hAsy(iAsy)%AsyTyp)%iBss
    
    haRay_Loc => haRay(iGeoTyp, icBss, jmRay)
    
    tNumCel = tNumCel + haRay_Loc%nhpRay
    
    DO iCel = 1, haRay_Loc%nhpRay
      ihPin = haRay_Loc%CelRay(iCel)%hPinIdx
      jhPin = hAsyTypInfo(iAsyTyp)%PinLocIdx(iGeoTyp, ihPin) + hAsy(iAsy)%PinIdxSt - 1 ! Global Pin Idx
      nSeg  = 0
      
      IF (ihPin > hAsyTypInfo(iAsyTyp)%nRodPin(1)) THEN
        DO iz = 1, nZ
          jcBss = gCelBss(gCel(GapPin(hPinInfo(jhPin)%PinTyp)%iCel(iz))%igBss)%icBss ! Gap Pin
          nSeg  = max(nSeg, haRay(iGeoTyp, jcBss, jmRay)%CelRay(iCel)%nSegRay)
        END DO
      ELSE
        DO iz = 1, nZ
          jcBss = hCel(RodPin(hPinInfo(jhPin)%PinTyp)%iCel(iz))%icBss ! Rod Pin
          nSeg  = max(nSeg, haRay(iGeoTyp, jcBss, jmRay)%CelRay(iCel)%nSegRay)
        END DO
      END IF
      
      tNumSeg = tNumSeg + nSeg
    END DO
  END DO
  
  nMaxCellRay = max(nMaxCellRay, tNumCel)
  nMaxRaySeg  = max(nMaxRaySeg,  tNumSeg)
END DO
!!$OMP END DO
!!$OMP END PARALLEL

NULLIFY (hcRay_Loc, haRay_Loc)
! ----------------------------------------------------

END SUBROUTINE HexSetCoreRay
! ------------------------------------------------------------------------------------------------------------
!                                     02. SET : Rot Ray
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetRotRay()

USE PARAM,   ONLY : ZERO, FALSE, TRUE
USE HexData, ONLY : ncRay, nRotRay, hcRay, hAsy, hRotRay, hLgc, hmRay, Asy2Dto1Dmap, Asy1Dto2DMap
USE HexUtil, ONLY : SetSgn_INT
USE Moc_Mod, ONLY : nMaxCoreRay

IMPLICIT NONE

INTEGER :: icRay, jcRay, imRay, jmRay, iAsy, jAsy, iGeo, jxAsy, jyAsy, iDir, jDir, tDir
INTEGER :: iNum, jNum, iRotRay
LOGICAL :: lErr

INTEGER, PARAMETER :: nMaxVacRotRay = 10 ! ARBITRARY

INTEGER, POINTER :: RayConn(:, :)
LOGICAL, POINTER :: lcRayUse(:)
! ----------------------------------------------------
TYPE Type_TmpRotRayLst
  INTEGER :: ncRay = 1
  
  INTEGER, POINTER :: RotRay(:)
END TYPE Type_TmpRotRayLst

TYPE(Type_TmpRotRayLst), POINTER :: tLst(:)
! ----------------------------------------------------

ALLOCATE (RayConn (-1:1, ncRay)); RayConn  = 0
ALLOCATE (lcRayUse      (ncRay)); lcRayUse = FALSE
ALLOCATE (tLst          (ncRay))
! ----------------------------------------------------
!               01. CONNECT : cRay
! ----------------------------------------------------
lErr = FALSE

!$OMP PARALLEL PRIVATE(icRay, tDir, iDir, iNum, imRay, iAsy, iGeo, jmRay, jDir, jxAsy, jyAsy, jAsy, jcRay, jNum)
!$OMP DO SCHEDULE (GUIDED)
DO icRay = 1, ncRay
  DO tDir = 1, 2
    iDir = 2*tDir - 3 ! Negative : y¢Ù, Positive : y¢Ö
    iNum = ((1 - iDir)/2) * 1 &                ! y¢Ù
         + ((iDir + 1)/2) * hcRay(icRay)%nmRay ! y¢Ö
    
    imRay = hcRay(icRay)%mRayIdx(iNum)
    iAsy  = hcRay(icRay)% AsyIdx(iNum)
    iGeo  = hAsy(iAsy)%GeoTyp
    
    jmRay = hmRay(imRay)%NxtMray_Ref(iDir, iGeo)
    jDir  = hmRay(imRay)% NxtDir_Ref(iDir, iGeo)
    
    IF (jDir .EQ. 0) CYCLE
    
    IF (hLgc%lAzmRot) iAsy = hAsy(iAsy)%RotNgh
    
    jxAsy = Asy1Dto2Dmap(1, iAsy) + hmRay(imRay)%NxtAsy_Ref(1, iDir, iGeo)
    jyAsy = Asy1Dto2Dmap(2, iAsy) + hmRay(imRay)%NxtAsy_Ref(2, iDir, iGeo)
    jAsy  = Asy2Dto1Dmap(jxAsy, jyAsy)
    
    IF (jAsy .EQ. 0) CYCLE
    
    DO jcRay = 1, ncRay
      IF (icRay .EQ. jcRay) CYCLE
      
      jNum = ((1 - jDir)/2) * hcRay(jcRay)%nmRay &  ! y¢Ù
           + ((jDir + 1)/2) * 1                     ! y¢Ö
      
      IF (hcRay(jcRay)% AsyIdx(jNum) .NE. jAsy) CYCLE
      IF (hcRay(jcRay)%mRayIdx(jNum) .NE. jmRay) CYCLE
      
      EXIT
    END DO
    
    IF (RayConn(iDir, icRay) .NE. 0) lErr = TRUE
    
    RayConn(iDir, icRay) = jcRay * jDir ! Negative : y¢Ù, Positive : y¢Ö
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

IF (lErr) CALL terminate("CONNECT CORE RAY")
! ----------------------------------------------------
!               02. SET : Rot Ray
! ----------------------------------------------------
nRotRay = 0
! ----------------------------
!      1. RAD BC : REF
! ----------------------------
IF (hLgc%lRadRef) THEN
  DO icRay = 1, ncRay
    IF(lcRayUse(icRay)) CYCLE
    
    nRotRay = nRotRay + 1
    jcRay   = icRay
    jDir    = 1
    
    ALLOCATE (tLst(nRotRay)%RotRay (ncRay))
    
    DO
      lcRayUse(abs(jcRay)) = TRUE
      
      tLst(nRotRay)%RotRay(tLst(nRotRay)%ncRay) = jcRay ! Negative : y¢Ù, Positive : y¢Ö
      
      jcRay = RayConn(jDir, abs(jcRay))
      
      IF(abs(jcRay) .EQ. icRay) EXIT ! Return to Start mRay
      
      jDir = SetSgn_INT(jcRay) ! Negative : y¢Ù, Positive : y¢Ö
      
      tLst(nRotRay)%ncRay = tLst(nRotRay)%ncRay + 1
    END DO
  END DO
! ----------------------------
!      2. RAD BC : VAC
! ----------------------------
ELSE
  DO icRay = 1, ncRay
    DO tDir = 1, 2
      iDir = 2*tDir - 3
      
      IF(lcRayUse(icRay)) CYCLE
      IF(RayConn(iDir, icRay) .NE. 0) CYCLE ! Start from Vacuum
      
      nRotRay = nRotRay + 1
      jDir    = -iDir ! NOTICE : Minus
      jcRay   = icRay * jDir
      
      ALLOCATE (tLst(nRotRay)%RotRay (nMaxVacRotRay))
      
      DO
        lcRayUse(abs(jcRay)) = TRUE
        
        tLst(nRotRay)%RotRay(tLst(nRotRay)%ncRay) = jcRay ! Negative : y¢Ù, Positive : y¢Ö
        
        IF(RayConn(jDir, abs(jcRay)) .EQ. 0) EXIT ! End to Vacuum
        
        jcRay = RayConn(jDir, abs(jcRay))
        jDir  = SetSgn_INT(jcRay) ! Negative : y¢Ù, Positive : y¢Ö
        
        tLst(nRotRay)%ncRay = tLst(nRotRay)%ncRay + 1
        
        IF (tLst(nRotRay)%ncRay > nMaxVacRotRay) CALL terminate("SET : ROT RAY - VAC")
      END DO
    END DO
  END DO
END IF
! ----------------------------------------------------
!               03. CP : Rot Ray
! ----------------------------------------------------
ALLOCATE(hRotRay (nRotRay))

!$OMP PARALLEL PRIVATE(iRotRay)
!$OMP DO SCHEDULE(GUIDED)
DO iRotRay = 1, nRotRay
  hRotRay(iRotRay)%ncRay = tLst(iRotRay)%ncRay
  
  ALLOCATE(hRotRay(iRotRay)%cRayIdx (tLst(iRotRay)%ncRay))
  
  ! Negative : y¢Ù, Positive : y¢Ö
  hRotRay(iRotRay)%cRayIdx(1:tLst(iRotRay)%ncRay) = tLst(iRotRay)%RotRay(1:tLst(iRotRay)%ncRay)
END DO
!$OMP END DO
!$OMP END PARALLEL

nMaxCoreRay = 0

DO iRotRay = 1, nRotRay
  nMaxCoreRay = max(nMaxCoreRay, hRotRay(iRotRay)%ncRay)
END DO

NULLIFY (RayConn, lcRayUse, tLst)
! ----------------------------------------------------

END SUBROUTINE HexSetRotRay
! ------------------------------------------------------------------------------------------------------------

END MODULE HexRayConst