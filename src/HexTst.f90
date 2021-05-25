#include <defines.h>
MODULE HexTst

USE HexUtil, ONLY : HexChkRange_INT, HexChkEqual_INT
USE PARAM,   ONLY : TRUE, FALSE, ZERO, HALF
USE ioutil,  ONLY : terminate

IMPLICIT NONE

CONTAINS

! ------------------------------------------------------------------------------------------------------------
!                                     01. TST : Hc Pin
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTstHcPin()

USE geom,    ONLY : ncbd
USE HexData, ONLY : nhcPin, hcPin

IMPLICIT NONE

INTEGER :: iPin

INTEGER, POINTER :: tnNgh(:), tnmPin(:), tmpIdx(:, :), tnBdmPin(:, :), tBdMPidx(:, :, :), tBdMPsuf(:, :, :)
INTEGER, POINTER :: tNghPin(:, :), tNghSuf(:, :)
REAL,    POINTER :: tNghLgh(:, :)
! ----------------------------------------------------

ALLOCATE (tnNgh          (nhcPin));    tnNgh = -1
ALLOCATE (tnmPin         (nhcPin));   tnmPin = -1
ALLOCATE (tmpIdx      (3, nhcPin));   tmpIdx = -1
ALLOCATE (tnBdmPin    (6, nhcPin)); tnBdmPin = -1
ALLOCATE (tBdMPidx (2, 6, nhcPin)); tBdMPidx = -1
ALLOCATE (tBdMPsuf (2, 6, nhcPin)); tBdMPsuf = -1
ALLOCATE (tNghPin  (ncbd, nhcPin));  tNghPin = -1
ALLOCATE (tNghSuf  (ncbd, nhcPin));  tNghSuf = -1
ALLOCATE (tNghLgh  (ncbd, nhcPin));  tNghLgh = ZERO

DO iPin = 1, nhcPin
  tnNgh         (iPin) = hcPin(iPin)%nNgh
  tnmPin        (iPin) = hcPin(iPin)%nmPin
  tmpIdx     (:, iPin) = hcPin(iPin)%mpIdx
  tnBdmPin   (:, iPin) = hcPin(iPin)%nBdmPin
  tBdMPidx(:, :, iPin) = hcPin(iPin)%BdMPidx
  tBdMPsuf(:, :, iPin) = hcPin(iPin)%BdMPsuf
  tNghPin    (:, iPin) = hcPin(iPin)%NghPin
  tNghSuf    (:, iPin) = hcPin(iPin)%NghSuf
  tNghLgh    (:, iPin) = hcPin(iPin)%NghLgh
END DO

NULLIFY (tnNgh, tnmPin, tmpIdx, tnBdmPin, tBdMPidx, tBdMPsuf, tNghPin, tNghSuf, tNghLgh)

CALL terminate("END OF TEST = CMFD CONST")
! ----------------------------------------------------

END SUBROUTINE HexTstHcPin
! ------------------------------------------------------------------------------------------------------------
!                                     02. TST : Asy Ray Seg Num
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTstAsyRaySegNum()

USE HexData, ONLY : ncBss, NumMray, haRay

IMPLICIT NONE

INTEGER :: icBss, imRay, iGeo, ihpRay

INTEGER, POINTER :: nSeg(:)
! ----------------------------------------------------

ALLOCATE (nSeg (ncBss)); nSeg = 0

DO icBss = 1, ncBss
  DO imRay = 1, NumMray(0)
    DO iGeo = 1, 7
      DO ihpRay = 1, haRay(iGeo, icBss, imRay)%nhpRay
        nSeg(icBss) = nSeg(icBss) + haRay(iGeo, icBss, imRay)%CelRay(ihpRay)%nSegRay
      END DO
    END DO
  END DO
END DO

PRINT *, nSeg(1)

DEALLOCATE (nSeg)

CALL terminate("END OF TEST = ASY RAY SEG NUM")
! ----------------------------------------------------

END SUBROUTINE HexTstAsyRaySegNum
! ------------------------------------------------------------------------------------------------------------
!                                     03. PRINT : Local Pin Typ
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexPrintPinTyp()

use HexData, ONLY : nHexPin, hAsy, hAsyTypInfo, hPinInfo, mpTypNumNgh

IMPLICIT NONE

INTEGER :: io, iAsy, iaTyp, iGeo, iPin, jPin, ivTyp, nPt, iPt

REAL :: Pts(2, 6)
! ----------------------------------------------------

io = 11

OPEN (UNIT = io, file = "TST Pin Typ")

WRITE(io, '(A, X, I6)') 'Pin Num', nHexPin

DO iPin = 1, nHexPin
  iAsy  = hPinInfo(iPin)%AsyIdx
  iaTyp = hPinInfo(iPin)%AsyTyp
  jPin  = hPinInfo(iPin)%OrdInAsy01
  iGeo  = hAsy(iAsy)%GeoTyp
  ivTyp = hAsyTypInfo(iaTyp)%PinVtxTyp(iGeo, jPin)
  nPt   = mpTypNumNgh(ivTyp)
  
  DO iPt = 1, nPt
    Pts(1:2, iPt) = hAsyTypInfo(iaTyp)%mpVtx(1:2, iPt, iGeo, jPin) + hAsy(iAsy)%Cnt(1:2)
  END DO
  
  IF (hPinInfo(iPin)%lRod) THEN
    WRITE (io, "(X, I6, X, I2, X, I2)") iPin,  hPinInfo(iPin)%PinTyp, nPt
  ELSE
    WRITE (io, "(X, I6, X, I2, X, I2)") iPin, -hPinInfo(iPin)%PinTyp, nPt
  END IF
  
  WRITE (io, "(X, 100(F12.5))") Pts(1:2, 1:nPt)
END DO

CLOSE (io)

CALL terminate("END OF TEST = PIN TYP")
! ----------------------------------------------------

END SUBROUTINE HexPrintPinTyp
! ------------------------------------------------------------------------------------------------------------
!                                     04. CHK : Eff Xs Out
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexChkEffXs()

USE CNTL,    ONLY : nTracerCntl
USE geom,    ONLY : nZ, Core, Pin, CellInfo
USE HexData, ONLY : nAsyCore, nhAsy, hCore, hAsy, hAsyTypInfo

IMPLICIT NONE

INTEGER :: i, iz, ixa, iya, ixc, iyc, xSt, xEd
INTEGER :: ixya, iaTyp, iaGeo, ixyc, ixy, iCel
INTEGER :: nPin, nRng, nFXR
! ----------------------------------------------------

IF (.NOT. nTracerCntl%lXsLib) CALL terminate("EFF XS : XS LIB")

DO i = 1, nTracerCntl%OutpCntl%nRegXsOut
  iz  = nTracerCntl%OutpCntl%RegXsOutList(1, i)
  ixa = nTracerCntl%OutpCntl%RegXsOutList(2, i)
  iya = nTracerCntl%OutpCntl%RegXsOutList(3, i)
  ixc = nTracerCntl%OutpCntl%RegXsOutList(4, i)
  iyc = nTracerCntl%OutpCntl%RegXsOutList(5, i)
  xSt = nTracerCntl%OutpCntl%RegXsOutList(6, i)
  xEd = nTracerCntl%OutpCntl%RegXsOutList(7, i)
  
  CALL HexChkRange_INT( iz, 1,       nZ, "EFF XS : IZ")
  CALL HexChkRange_INT(ixa, 1, Core%nya, "EFF XS : IXA")
  CALL HexChkRange_INT(iya, 1, Core%nya, "EFF XS : IYA")
  
  ixya = hCore(ixa, iya)
  
  CALL HexChkRange_INT(ixya, 1, nhAsy, "EFF XS : IXYA")
  
  iaTyp = hAsy(ixya)%AsyTyp
  iaGeo = hAsy(ixya)%GeoTyp
  !nPin  = hAsy(ixya)%nPin
  nRng  = 2 * nPin - 1
  
  CALL HexChkRange_INT(ixc, 1, nRng, "EFF XS : IXC")
  CALL HexChkRange_INT(iyc, 1, nRng, "EFF XS : IYC")
  
  ixyc = hAsyTypInfo(iaTyp)%Pin2Dto1Dmap(ixc, iyc)
  ixyc = hAsyTypInfo(iaTyp)%PinLocIdx(iaGeo, ixyc)
  
  CALL HexChkRange_INT(ixyc, 1, hAsy(ixya)%nTotPin, "EFF XS : IXYC")
  
  ixy  = hAsy(ixya)%PinIdxSt + ixyc - 1
  iCel = Pin(ixy)%Cell(iz)
  nFXR = CellInfo(iCel)%nFXR
  
  CALL HexChkRange_INT(xSt, 1, nFXR, "EFF XS : FXR BOUNDARY")
  CALL HexChkRange_INT(xEd, 1, nFXR, "EFF XS : FXR BOUNDARY")
END DO
! ----------------------------------------------------

END SUBROUTINE HexChkEffXs
! ------------------------------------------------------------------------------------------------------------
!                                     05. TST : hPin Info
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTsthPinInfo()

USE HexData, ONLY : nHexPin, hPinInfo

IMPLICIT NONE

INTEGER, POINTER :: tPinTyp(:), tAsyIdx(:), tAsyTyp(:), tOrdInAsy01(:), tix(:), tiy(:), tVtxTyp(:)
INTEGER, POINTER :: tnSct(:), tFsrIdxSt(:), tFxrIdxSt(:)
REAL, POINTER :: tWt(:), tCnt(:,:), tVol(:), tVolFm(:)
LOGICAL, POINTER :: tlInn(:), tlBndy(:), tlRod(:), tlGap(:)

INTEGER :: iPin
! ----------------------------------------------------

ALLOCATE(tPinTyp     (nHexPin))
ALLOCATE(tAsyIdx     (nHexPin))
ALLOCATE(tAsyTyp     (nHexPin))
ALLOCATE(tOrdInAsy01 (nHexPin))
ALLOCATE(tix         (nHexPin))
ALLOCATE(tiy         (nHexPin))
ALLOCATE(tVtxTyp     (nHexPin))
ALLOCATE(tnSct       (nHexPin))
ALLOCATE(tFsrIdxSt   (nHexPin))
ALLOCATE(tFxrIdxSt   (nHexPin))
ALLOCATE(tWt         (nHexPin))
ALLOCATE(tCnt     (2, nHexPin))
ALLOCATE(tVol        (nHexPin))
ALLOCATE(tVolFm      (nHexPin))
ALLOCATE(tlInn       (nHexPin))
ALLOCATE(tlBndy      (nHexPin))
ALLOCATE(tlRod       (nHexPin))
ALLOCATE(tlGap       (nHexPin))

DO iPin = 1, nHexPin
  tPinTyp     (iPin) = hPinInfo(iPin)%PinTyp
  tAsyIdx     (iPin) = hPinInfo(iPin)%AsyIdx
  tAsyTyp     (iPin) = hPinInfo(iPin)%AsyTyp
  tVtxTyp     (iPin) = hPinInfo(iPin)%VtxTyp
  tOrdInAsy01 (iPin) = hPinInfo(iPin)%OrdInAsy01
  tix         (iPin) = hPinInfo(iPin)%ix
  tiy         (iPin) = hPinInfo(iPin)%iy
  tnSct       (iPin) = hPinInfo(iPin)%nSct
  tFsrIdxSt   (iPin) = hPinInfo(iPin)%FsrIdxSt
  tFxrIdxSt   (iPin) = hPinInfo(iPin)%FxrIdxSt
  tWt         (iPin) = hPinInfo(iPin)%Wt
  tCnt   (1:2, iPin) = hPinInfo(iPin)%Cnt(1:2)
  tVol        (iPin) = hPinInfo(iPin)%Vol(1)
  tVolFm      (iPin) = hPinInfo(iPin)%VolFm(1)
  tlInn       (iPin) = hPinInfo(iPin)%lInn
  tlBndy      (iPin) = hPinInfo(iPin)%lBndy
  tlRod       (iPin) = hPinInfo(iPin)%lRod
  tlGap       (iPin) = hPinInfo(iPin)%lGap
END DO

NULLIFY(tPinTyp, tAsyIdx, tAsyTyp, tOrdInAsy01, tix, tiy, tVtxTyp, tnSct, tFsrIdxSt, tFxrIdxSt)
NULLIFY(tWt, tCnt, tVol, tVolFm, tlInn, tlBndy, tlRod, tlGap)

CALL terminate("END OF TEST = HPIN INFO")
! ----------------------------------------------------

END SUBROUTINE HexTsthPinInfo
! ------------------------------------------------------------------------------------------------------------
!                                     06. TST : hmRay
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTsthmRay()

USE HexData, ONLY : hmRay, NumMray

IMPLICIT NONE

INTEGER :: imRay

INTEGER, POINTER :: tAzmIdx(:)

REAL, POINTER :: tEq(:, :)
REAL, POINTER :: tPt(:, :, :)

INTEGER, POINTER :: tNxtAsy_Mov(:, :, :)
INTEGER, POINTER :: tNxtmRay_Mov(:, :)

INTEGER, POINTER :: tNxtAsy_Ref(:, :, :, :)
INTEGER, POINTER :: tNxtmRay_Ref(:, :, :)
INTEGER, POINTER :: tNxtDir_Ref(:, :, :)
! ----------------------------------------------------

ALLOCATE (tAzmIdx                 (NumMray(0)))
ALLOCATE (tEq                  (3, NumMray(0)))
ALLOCATE (tPt               (2, 2, NumMray(0)))
ALLOCATE (tNxtAsy_Mov    (2, -1:1, NumMray(0)))
ALLOCATE (tNxtmRay_Mov      (-1:1, NumMray(0)))
ALLOCATE (tNxtAsy_Ref (2, -1:1, 7, NumMray(0)))
ALLOCATE (tNxtmRay_Ref   (-1:1, 7, NumMray(0)))
ALLOCATE (tNxtDir_Ref    (-1:1, 7, NumMray(0)))

DO imRay = 1, NumMray(0)
  tAzmIdx              (imRay) = hmRay(imRay)%AzmIdx
  tEq               (:, imRay) = hmRay(imRay)%Eq               (:)
  tPt            (:, :, imRay) = hmRay(imRay)%Pt            (:, :)
  tNxtAsy_Mov    (:, :, imRay) = hmRay(imRay)%NxtAsy_Mov    (:, :)
  tNxtmRay_Mov      (:, imRay) = hmRay(imRay)%NxtmRay_Mov      (:)
  tNxtAsy_Ref (:, :, :, imRay) = hmRay(imRay)%NxtAsy_Ref (:, :, :)
  tNxtmRay_Ref   (:, :, imRay) = hmRay(imRay)%NxtmRay_Ref   (:, :)
  tNxtDir_Ref    (:, :, imRay) = hmRay(imRay)%NxtDir_Ref    (:, :)
END DO

NULLIFY (tAzmIdx, tEq, tPt, tNxtAsy_Mov, tNxtmRay_Mov, tNxtAsy_Ref, tNxtmRay_Ref, tNxtDir_Ref)

CALL terminate("END OF TEST = HMRAY")
! ----------------------------------------------------

END SUBROUTINE HexTsthmRay
! ------------------------------------------------------------------------------------------------------------
!                                     07. TST : hpRay
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTsthpRay(RayEqn, nPin, hpRay)

USE HexType, ONLY : Type_HexPinRay

IMPLICIT NONE

INTEGER, POINTER :: tPinIdx(:)
INTEGER, POINTER :: tSurfIdx(:, :)
REAL,    POINTER :: tPinPt (:, :, :)

INTEGER :: nPin, iPin
REAL    :: RayEqn(3)

TYPE(Type_HexPinRay) :: hpRay(:)
! ----------------------------------------------------

ALLOCATE (tPinIdx      (nPin))
ALLOCATE (tSurfIdx  (2, nPin))
ALLOCATE (tPinPt (2, 2, nPin))

DO iPin = 1, nPin
  tPinIdx          (iPin) = hpRay(iPin)%PinIdx
  tSurfIdx    (1:2, iPin) = hpRay(iPin)%SurfIdx(1:2)
  tPinPt (1:2, 1:2, iPin) = hpRay(iPin)%PinPt (1:2, 1:2)
END DO

NULLIFY (tPinIdx, tSurfIdx, tPinPt)

CALL terminate("END OF TEST = HPRAY")
! ----------------------------------------------------

END SUBROUTINE HexTsthpRay
! ------------------------------------------------------------------------------------------------------------
!                                     08. TST : haRay
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTsthaRay(iGeo, icBss, imRay)

USE HexType, ONLY : Type_HexAsyRay
USE HexData, ONLY : haRay, hmRay

IMPLICIT NONE

INTEGER :: iGeo, icBss, imRay
REAL    :: RayEqn(3)

INTEGER :: ihpRay, ihsRay, nhsRay

INTEGER, POINTER :: MshIdx(:)
!REAL,    POINTER :: SegPts(:, :)
REAL,    POINTER :: SegLgh(:)
TYPE(Type_HexAsyRay), POINTER :: haRay_Loc
! ----------------------------------------------------

haRay_Loc  => haRay(iGeo, icBss, imRay)
RayEqn(1:3) = hmRay(imRay)%Eq

IF (haRay_Loc%nhpRay .LT. 1) CALL terminate("INVALID GEO TYP & MRAY")

nhsRay = 0

DO ihpRay = 1, haRay_Loc%nhpRay
  nhsRay = nhsRay + haRay_Loc%CelRay(ihpRay)%nSegRay
END DO

ALLOCATE (MshIdx    (nhsRay))
ALLOCATE (SegLgh    (nhsRay))
!ALLOCATE (SegPts (2, nhsRay + 1))

nhsRay = 0

DO ihpRay = 1, haRay_Loc%nhpRay
  DO ihsRay = 1, haRay_Loc%CelRay(ihpRay)%nSegRay
    nhsRay = nhsRay + 1

    MshIdx     (nhsRay) = haRay_Loc%CelRay(ihpRay)%MshIdx(ihsRay)
    SegLgh     (nhsRay) = haRay_Loc%CelRay(ihpRay)%SegLgh(ihsRay)
    !SegPts(1:2, nhsRay) = haRay_Loc%CelRay(ihpRay)%SegPts(1:2, ihsRay)
  END DO
END DO

!SegPts(1:2, nhsRay + 1) = haRay_Loc%CelRay(haRay_Loc%nhpRay)%SegPts(1:2, ihsRay)

NULLIFY (MshIdx, SegLgh, haRay_Loc)
!NULLIFY (SegPts)

CALL terminate("END OF TEST = HARAY")
! ----------------------------------------------------

END SUBROUTINE HexTsthaRay
! ------------------------------------------------------------------------------------------------------------
!                                     09. TST : hcRay
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTsthcRay()

USE PARAM,   ONLY : ZERO
USE HexType, ONLY : Type_HexCoreRay
USE HexData, ONLY : nAzmAng, hcRay, ncRay, hAsy, hmRay

IMPLICIT NONE

INTEGER :: io, iAng, icRay, jcRay, imRay, iAsy, nMaxRay, nmRay
REAL    :: aCnt(2), Pts(2, 2)

INTEGER, POINTER :: cRayRng(:, :) ! (St/Ed, iAng)

TYPE(Type_HexCoreRay), POINTER :: hcRay_Loc
! ----------------------------------------------------

io = 11

OPEN (UNIT = io, file = "TST Core Ray Pts")
! ----------------------------------------------------
!               01. SET : cRay Range
! ----------------------------------------------------
ALLOCATE (cRayRng (2, nAzmAng))

nMaxRay = 0

DO iAng = 1, nAzmAng
  DO icRay = 1, ncRay
    IF (hcRay(icRay)%AzmIdx .EQ. iAng) EXIT
  END DO

  DO jcRay = icRay + 1, ncRay
    IF (hcRay(jcRay)%AzmIdx .NE. iAng) EXIT
  END DO

  cRayRng(1, iAng) = icRay
  cRayRng(2, iAng) = jcRay - 1

  nMaxRay = max(nMaxRay, jcRay - icRay)
END DO
! ----------------------------------------------------
!               02. PLOT
! ----------------------------------------------------
DO iAng = 1, nAzmAng
  DO icRay = cRayRng(1, iAng), cRayRng(2, iAng)
    hcRay_Loc => hcRay(icRay)

    jcRay = icRay + 1 - cRayRng(1, iAng)
    imRay = hcRay_Loc%mRayIdx(1)
    iAsy  = hcRay_Loc% AsyIdx(1)
    aCnt  = hAsy(iAsy)%Cnt

    Pts(1:2, 1) = hmRay(imRay)%Pt(1:2, 1) + aCnt(1:2)

    nmRay = hcRay_Loc%nmRay
    imRay = hcRay_Loc%mRayIdx(nmRay)
    iAsy  = hcRay_Loc% AsyIdx(nmRay)
    aCnt  = hAsy(iAsy)%Cnt

    Pts(1:2, 2) = hmRay(imRay)%Pt(1:2, 2) + aCnt(1:2)

    WRITE(io, '(2(X, I5), 4(X, F12.5))') iAng, jcRay, Pts(1:2, 1), Pts(1:2, 2)
  END DO
END DO

NULLIFY (cRayRng, hcRay_Loc)

CLOSE(io)

CALL terminate("END OF TEST = HCRAY")
! ----------------------------------------------------

END SUBROUTINE HexTsthcRay
! ------------------------------------------------------------------------------------------------------------
!                                     10. TST : Rot Ray Pts
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTsthRotRay()

USE HexType, ONLY : Type_HexRotRay, Type_HexCoreRay
USE HexData, ONLY : hRotRay, nRotRay, hcRay, hmRay, hAsy

IMPLICIT NONE

INTEGER :: io, iRotRay, icRay, iDir, imRay, jmRay, iAsy
REAL    :: Pts(2, 2)

TYPE(Type_HexRotRay),  POINTER :: hrRay_Loc
TYPE(Type_HexCoreRay), POINTER :: hcRay_Loc
! ----------------------------------------------------

io = 11

OPEN (UNIT = io, file = "TST Rot Ray")

DO iRotRay = 1, nRotRay
  hrRay_Loc => hRotRay(iRotRay)

  WRITE(io, '(4(X, I5))') iRotRay, hrRay_Loc%ncRay, &
                         hcRay(abs(hrRay_Loc%cRayIdx(1)))%AzmIdx, &
                         hcRay(abs(hrRay_Loc%cRayIdx(hrRay_Loc%ncRay)))%AzmIdx

  DO icRay = 1, hrRay_Loc%ncRay
    hcRay_Loc => hcRay(abs(hrRay_Loc%cRayIdx(icRay)))

    DO iDir = 1, 2
      imRay = (2 - iDir) * 1 &             ! y¢Ù
            + (iDir - 1) * hcRay_Loc%nmRay ! y¢Ö

      iAsy  = hcRay_Loc% AsyIdx(imRay)
      jmRay = hcRay_Loc%mRayIdx(imRay)

      Pts(1:2, iDir) = hmRay(jmRay)%Pt(1:2, iDir) + hAsy(iAsy)%Cnt(1:2)
    END DO

    WRITE(io, '(2(X, I5), 4(X, F12.5))') icRay, hrRay_Loc%cRayIdx(icRay), Pts(1:2, 1:2)
  END DO

  WRITE(io, *)
END DO

NULLIFY (hrRay_Loc, hcRay_Loc)

CLOSE(io)

CALL terminate("END OF TEST = ROT RAY PTS")
! ----------------------------------------------------

END SUBROUTINE HexTsthRotRay
! ------------------------------------------------------------------------------------------------------------
!                                     11. TST : Pin Vtx
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTstPinVtx(iaTyp)

USE HexType, ONLY : Type_HexAsyTypInfo
USE HexData, ONLY : hAsyTypInfo, nGeoTyp, mpTypNumNgh

IMPLICIT NONE

INTEGER :: iaTyp, io, iGeo, iPin, ivTyp, nBndy

TYPE(Type_HexAsyTypInfo), POINTER :: aInf_Loc
! ----------------------------------------------------

io = 11

OPEN (UNIT = io, file = "TST PIN VTX")

aInf_Loc => hAsyTypInfo(iaTyp)

DO iGeo = 1, nGeoTyp
  WRITE(io, *) iGeo, aInf_Loc%nTotPin(1)

  DO iPin = 1, aInf_Loc%nTotPin(1)
    IF (.NOT. aInf_Loc%lGeoPin(iGeo, iPin)) THEN
      WRITE(io, *) iPin, 0
    ELSE
      ivTyp = aInf_Loc%PinVtxTyp(iGeo, iPin)
      nBndy = mpTypNumNgh(ivTyp)

      WRITE(io, '(2(X, I5), 14(X, F12.5))') iPin, nBndy, aInf_Loc%mpVtx(1:2, 1:nBndy+1, iGeo, iPin)
    END IF
  END DO
END DO

CLOSE(io)

NULLIFY (aInf_Loc)

CALL terminate("END OF TEST = PIN VTX")
! ----------------------------------------------------

END SUBROUTINE HexTstPinVtx
! ------------------------------------------------------------------------------------------------------------
!                                     12. TST : Super-Pin Vtx
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexTstSpVtx(iaTyp)

USE HexType, ONLY : Type_HexAsyTypInfo
USE HexData, ONLY : hAsyTypInfo, nGeoTyp, spTypNumNgh

IMPLICIT NONE

INTEGER :: iaTyp, io, iGeo, iPin, ivTyp, nBndy

TYPE(Type_HexAsyTypInfo), POINTER :: aInf_Loc
! ----------------------------------------------------

io = 11

OPEN (UNIT = io, file = "TST PIN VTX")

aInf_Loc => hAsyTypInfo(iaTyp)

DO iGeo = 1, nGeoTyp
  WRITE(io, *) iGeo, aInf_Loc%nRodPin(1)

  DO iPin = 1, aInf_Loc%nRodPin(1)
    IF (.NOT. aInf_Loc%lGeoPin(iGeo, iPin)) THEN
      WRITE(io, *) iPin, 0
    ELSE
      ivTyp = aInf_Loc%PinVtxTyp(iGeo, iPin)
      nBndy = spTypNumNgh(ivTyp)

      WRITE(io, '(2(X, I5), 14(X, F12.5))') iPin, nBndy, aInf_Loc%spVtx(1:2, 1:nBndy+1, iGeo, iPin)
    END IF
  END DO
END DO

CLOSE(io)

NULLIFY (aInf_Loc)

CALL terminate("END OF TEST = SUPER PIN VTX")
! ----------------------------------------------------

END SUBROUTINE HexTstSpVtx
! ------------------------------------------------------------------------------------------------------------
!                                     13. PRT : Flx
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexPrtFlx()

USE Core_mod, ONLY : FmInfo
USE geom,     ONLY : Core, nz, ng

IMPLICIT NONE

INTEGER :: io, iFsr, ig, iz
INTEGER, SAVE :: itr
LOGICAL, SAVE :: lFst
! ----------------------------------------------------

io  = 11
itr = itr + 1
ig  = 2 ! Fixed

IF (lFst) THEN
  lFst = FALSE
  itr  = 1

  OPEN (UNIT = io, file = "TST FLX")
END IF

WRITE (io, *) itr
WRITE (*,  *) itr
WRITE (io, *)

WRITE (io, '(X, A8, X, 100(I17, X))') 'FSR / iz', (iz, iz = 1, nz)

DO iFsr = 1, Core%nCoreFsr
  WRITE (io, '(X, I8, X, 100(E17.5, X))') iFsr, FmInfo%Phis(iFsr, 1:nz, ig)
END DO

WRITE (io, *)
! ----------------------------------------------------

END SUBROUTINE HexPrtFlx
! ------------------------------------------------------------------------------------------------------------
!                                     14. PRT : mkl Phi C
! ------------------------------------------------------------------------------------------------------------
#ifdef __INTEL_MKL
SUBROUTINE HexPrtMklPinXsPhi()

USE MKL_3D

IMPLICIT NONE

INTEGER :: io, ixy, ig, iz, ixy_map
INTEGER, SAVE :: itr
LOGICAL, SAVE :: lFst
! ----------------------------------------------------

io  = 11
itr = itr + 1
ig  = 2 ! Fixed

IF (lFst) THEN
  lFst = FALSE
  itr  = 1

  OPEN (UNIT = io, file = "TST MKL PIN XS PHI")
END IF

WRITE (io, *) itr
WRITE (*,  *) itr
WRITE (io, *)
WRITE (io, '(X, A8, X, 100(I17, X))') 'Pin / iz', (iz, iz = 1, mklGeom%nzCMFD)

DO ixy = 1, mklGeom%nxy
  ixy_map = mklGeom%pinMap(ixy)

  WRITE (io, '(X, I8, X, 100(E17.5, X))') ixy, (mklCMFD%PinXS(ixy_map, iz)%Phi(ig), iz = 1, mklGeom%nzCMFD)
END DO

WRITE (io, *)
! ----------------------------------------------------

END SUBROUTINE HexPrtMklPinXsPhi
! ------------------------------------------------------------------------------------------------------------
!                                     15. PRT : mkl Phi S
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexPrtMklPhiS()

USE MKL_3D

IMPLICIT NONE

INTEGER :: io, ixy, ig, iz, ixy_map
INTEGER, SAVE :: itr
LOGICAL, SAVE :: lFst
! ----------------------------------------------------

io  = 11
itr = itr + 1
ig  = 2 ! Fixed

IF (lFst) THEN
  lFst = FALSE
  itr  = 1

  OPEN (UNIT = io, file = "TST MKL PHI S")
END IF

WRITE (io, *) itr
WRITE (*,  *) itr
WRITE (io, *)
WRITE (io, '(X, A8, X, 100(I17, X))') 'Pin / iz', (iz, iz = 1, mklGeom%nzCMFD)

DO ixy = 1, mklGeom%nxy
  ixy_map = mklGeom%pinMap(ixy)

  WRITE (io, '(X, I8, X, 100(E17.5, X))') ixy, (mklCMFD%phis(ixy_map, iz, ig), iz = 1, mklGeom%nzCMFD)
END DO

WRITE (io, *)
! ----------------------------------------------------

END SUBROUTINE HexPrtMklPhiS
! ------------------------------------------------------------------------------------------------------------
!                                     15. PRT : mkl Phi S
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexPrtCSRval(csrDouble)

USE MKL_3D

IMPLICIT NONE

TYPE(CSR_DOUBLE), POINTER :: csrDouble(:)

INTEGER :: io, iCol, iRow, idxNN, NN, ig
INTEGER, SAVE :: itr
LOGICAL, SAVE :: lFst
! ----------------------------------------------------

io  = 11
itr = itr + 1

IF (lFst) THEN
  lFst = FALSE
  itr  = 1

  OPEN (UNIT = io, file = "TST CSR VAL")
END IF

WRITE (io, *) itr
WRITE (*,  *) itr
WRITE (io, *)

DO ig = 1, 47
  NN    = 10
  idxNN = csrDouble(ig)%nnz / NN

  WRITE (io, '(X, A8, X, 100(I17, X))') 'ig',  ig
  WRITE (io, '(X, A8, X, 100(I17, X))') 'idx', (iRow, iRow = 1, NN)

  DO iCol = 1, idxNN
    WRITE (io, '(X, I8, X, 100(E17.5, X))') NN*(iCol-1), (csrDouble(ig)%csrval(NN*(iCol-1)+iRow), iRow = 1, NN)
  END DO

  WRITE (io, '(X, I8, X, 100(E17.5, X))') idxNN*NN, (csrDouble(ig)%csrval(idxNN*NN+iRow), iRow = 1, csrDouble(ig)%nnz - idxNN*NN)
  WRITE (io, *)
END DO
! ----------------------------------------------------

END SUBROUTINE HexPrtCSRval
#endif
! ------------------------------------------------------------------------------------------------------------
!                                     15. PRT : mkl Phi S
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexPrtArray_REAL_2D(aa, Array, d1, d2)

IMPLICIT NONE

REAL :: Array(:, :)
INTEGER :: d1, d2
CHARACTER*256 :: aa

INTEGER :: io, iCol, iRow, idxNN, NN, ig
INTEGER, SAVE :: itr
LOGICAL, SAVE :: lFst
! ----------------------------------------------------

io  = 11
itr = itr + 1

IF (lFst) THEN
  lFst = FALSE
  itr  = 1

  OPEN (UNIT = io, file = "TST ARRAY VAL")
END IF

WRITE (io, *) itr
WRITE (*,  *) itr
WRITE (io, *)

WRITE (io, '(X, A8, X, 1000(I17, X))') aa, (iRow, iRow = 1, d1)

DO iCol = 1, d2
  WRITE (io, '(X, I8, X, 1000(E17.5, X))') iCol, (Array(iRow, iCol), iRow = 1, d1)
END DO

WRITE (io, *)
! ----------------------------------------------------

END SUBROUTINE HexPrtArray_REAL_2D
! ------------------------------------------------------------------------------------------------------------
!                                     15. PRT : PHI C
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexPrtPhiC(myRank, myzb, myze)

USE GEOM,     ONLY : Core, ng, nz
USE Core_mod, ONLY : CmInfo
USE HexData,  ONLY : hAsy, hPinInfo
USE PE_MOD,   ONLY : PE

IMPLICIT NONE

INTEGER :: myzb, myze, myRank
INTEGER :: iAsy, iz, ixy, io, AsyType, jxy, iPin
CHARACTER*256 :: fName
! ----------------------------------------------------

io = 11 + myRank

WRITE(fName, "(A11, I2)"), "TST PHI C -", myRank

OPEN (UNIT = io, file = fName)

DO iz = myzb, myze   !Axial Sweep
  IF(.NOT. Core%lFuelPlane(iz)) CYCLE
  
  DO iAsy = 1, Core%nxya
    AsyType = Core%Asy(iasy)%AsyType
    
    IF(.NOT. Core%AsyInfo(AsyType)%lFuel) CYCLE
    
    DO ixy = 1, hAsy(iAsy)%nTotPin  !Pin Sweep
      jxy  = hAsy(iAsy)%PinIdxSt - 1 + ixy
      iPin = hPinInfo(jxy)%OrdInAsy01
      
      WRITE (io, '(2(X, I4), 100(X, E17.5))') iz, jxy, CmInfo%PhiC(jxy, iz, 1:ng)
    END DO 
  END DO
END DO

CLOSE (io)
! ----------------------------------------------------

END SUBROUTINE HexPrtPhiC
! ------------------------------------------------------------------------------------------------------------
!                                     15. PRT : PHI C
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexPrtXsKF(myRank, myzb, myze)

USE GEOM,     ONLY : Core, ng, nz
USE Core_mod, ONLY : CmInfo
USE HexData,  ONLY : hAsy, hPinInfo
USE PE_MOD,   ONLY : PE

IMPLICIT NONE

INTEGER :: myzb, myze, myRank
INTEGER :: iAsy, iz, ixy, io, AsyType, jxy, iPin
CHARACTER*256 :: fName
! ----------------------------------------------------

io = 11 + myRank

WRITE(fName, "(A11, I2)"), "TST XS KF -", myRank

OPEN (UNIT = io, file = fName)

DO iz = 1, nz   !Axial Sweep
  IF(.NOT. Core%lFuelPlane(iz)) CYCLE
  
  DO iAsy = 1, Core%nxya
    AsyType = Core%Asy(iasy)%AsyType
    
    IF(.NOT. Core%AsyInfo(AsyType)%lFuel) CYCLE
    
    DO ixy = 1, hAsy(iAsy)%nTotPin  !Pin Sweep
      jxy  = hAsy(iAsy)%PinIdxSt - 1 + ixy
      iPin = hPinInfo(jxy)%OrdInAsy01
      
      WRITE (io, '(2(X, I4), 100(X, E17.5))') iz, jxy, CmInfo%PinXs(jxy, iz)%XsKF(1:ng)
    END DO 
  END DO
END DO

CLOSE (io)
! ----------------------------------------------------

END SUBROUTINE HexPrtXsKF
! ------------------------------------------------------------------------------------------------------------

END MODULE HexTst