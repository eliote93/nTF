#include <defines.h>
MODULE HexCmfdConst

IMPLICIT NONE

CONTAINS

! ------------------------------------------------------------------------------------------------------------
!                                     01. SET : Super Pin
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetHcPin()

USE HexData, ONLY : hLgc

USE HexType, ONLY : Type_HexAsyTypInfo, Type_HexCmfdPin
USE HexData, ONLY : nhAsy, hAsy, hAsyTypInfo, nHexPin, hcPin, hPinInfo, spTypNumNgh, nhcPin
USE HexUtil, ONLY : SetEqn, FindCnt

IMPLICIT NONE
! ----------------------------------------------------

! CASE : Sng Cel
IF (hLgc%lSngCel) THEN
  CALL HexSetHcPin_Sng

  RETURN
END IF

! SET : hcPin Self
IF (hLgc%lspCMFD) THEN
  CALL HexSetHcPinSlf_SP
ELSE
  CALL HexSetHcPinSlf_MP
END IF

!  SET : hcPin Ngh
CALL HexSetHcPinNgh
! ----------------------------------------------------

END SUBROUTINE HexSetHcPin
! ------------------------------------------------------------------------------------------------------------
!                                     02. SET : Super Pin Self - SP
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetHcPinSlf_SP()

USE allocs
USE HexType, ONLY : Type_HexAsyTypInfo, Type_HexCmfdPin
USE HexData, ONLY : nhAsy, hAsy, hAsyTypInfo, nHexPin, hcPin, hPinInfo, spTypNumNgh, nhcPin
USE HexUtil, ONLY : SetEqn, FindCnt

IMPLICIT NONE

INTEGER :: iAsy, iGeo, iPin, jPin, kPin, tPin, ivTyp, iBndy, nPin, nBndy
REAL    :: Cnt(2)

INTEGER, POINTER :: cpSlfMPnum(:, :), cpSlfMPidx(:, :, :), PinMap(:) ! Global to Super-Pin
INTEGER, POINTER :: cpSufMPnum(:, :, :), cpSufMPidx(:, :, :, :), cpSufMPsuf(:, :, :, :)

TYPE(Type_HexAsyTypInfo), POINTER :: aInf_Loc
TYPE(Type_HexCmfdPin),    POINTER :: cPin_Loc
! ----------------------------------------------------

CALL dmalloc(PinMap, nHexPin)

nPin = 0

DO iAsy = 1, nhAsy
  nhcPin = nhcPin + hAsy(iAsy)%nRodPin

  DO iPin = 1, hAsy(iAsy)%nRodPin
    nPin = nPin + 1
    jPin = iPin + hAsy(iAsy)%PinIdxSt - 1

    PinMap(jPin) = nPin
  END DO
END DO

ALLOCATE (hcPin (nhcPin))
! ----------------------------------------------------
DO iAsy = 1, nhAsy
  iGeo = hAsy(iAsy)%GeoTyp

  aInf_Loc => hAsyTypInfo(hAsy(iAsy)%AsyTyp)

  cpSufMPnum => aInf_Loc%cpSufMPnum
  cpSlfMPnum => aInf_Loc%cpSlfMPnum
  cpSlfMPidx => aInf_Loc%cpSlfMPidx
  cpSufMPidx => aInf_Loc%cpSufMPidx
  cpSufMPsuf => aInf_Loc%cpSufMPsuf

  !$OMP PARALLEL PRIVATE(iPin, jPin, kPin, ivTyp, nBndy, cPin_Loc, Cnt, iBndy, tPin)
  !$OMP DO SCHEDULE(GUIDED)
  DO iPin = 1, hAsy(iAsy)%nRodPin
    jPin  = iPin + hAsy(iAsy)%PinIdxSt - 1  ! Global Pin Idx
    kPin  = hPinInfo(jPin)%OrdInAsy01       ! Pin Idx in "hAsyTypInfo"
    ivTyp = aInf_Loc%PinVtxTyp(iGeo, kPin)  ! Vtx Typ
    nBndy = spTypNumNgh(ivTyp)

    cPin_Loc => hcPin(PinMap(jPin))
    ! ----------------------------
    !      1. CP
    ! ----------------------------
    cPin_Loc%nBndy = nBndy
    cPin_Loc%aIdx  = iAsy

    cPin_Loc%BdPts(1, 1:7) = aInf_Loc%spVtx(1, 1:7, iGeo, kPin) + hAsy(iAsy)%Cnt(1)
    cPin_Loc%BdPts(2, 1:7) = aInf_Loc%spVtx(2, 1:7, iGeo, kPin) + hAsy(iAsy)%Cnt(2)

    Cnt(1:2) = FindCnt(nBndy, cPin_Loc%BdPts(1:2, 1:nBndy))

    DO iBndy = 1, nBndy
      cPin_Loc%BdLgh(iBndy) = aInf_Loc%spBndyLgh(iBndy, ivTyp)
      cPin_Loc%BdC2B(iBndy) = aInf_Loc%spBndyC2B(iBndy, ivTyp)

      cPin_Loc%BdEqn(1:3, iBndy) = SetEqn(cPin_Loc%BdPts(1:2, iBndy), cPin_Loc%BdPts(1:2, iBndy+1), Cnt)
    END DO
    ! ----------------------------
    !      2. MP
    ! ----------------------------
    cPin_Loc%Area = aInf_Loc%spTypAre(ivTyp)

    cPin_Loc%nmPin   = cpSlfMPnum   (iGeo, kPin)
    cPin_Loc%nBdmPin = cpSufMPnum(:, iGeo, kPin)

    DO jPin = 1, cPin_Loc%nmPin
      tPin = cpSlfMPidx(jPin, iGeo, kPin)

      cPin_Loc%mpIdx(jPin) = aInf_Loc%PinLocIdx(iGeo, tPin) + hAsy(iAsy)%PinIdxSt - 1
    END DO

    DO iBndy = 1, cPin_Loc%nBndy
      DO jPin = 1, cpSufMPnum(iBndy, iGeo, kPin)
        tPin = cpSufMPidx(jPin, iBndy, iGeo, kPin)

        cPin_Loc%BdMPsuf(jPin, iBndy) = cpSufMPsuf(jPin, iBndy, iGeo, kPin)
        cPin_Loc%BdMPidx(jPin, iBndy) = aInf_Loc%PinLocIdx(iGeo, tPin) + hAsy(iAsy)%PinIdxSt - 1 ! Global Pin Idx
      END DO
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
END DO

NULLIFY (PinMap, cpSufMPnum, cpSlfMPnum, cpSlfMPidx, cpSufMPidx, cpSufMPsuf, aInf_Loc, cPin_Loc)
! ----------------------------------------------------

END SUBROUTINE HexSetHcPinSlf_SP
! ------------------------------------------------------------------------------------------------------------
!                                     03. SET : Super Pin Self - MP
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetHcPinSlf_MP()

USE allocs
USE HexType, ONLY : Type_HexAsyTypInfo, Type_HexCmfdPin
USE HexData, ONLY : nhAsy, hAsy, hAsyTypInfo, nHexPin, hcPin, hPinInfo, mpTypNumNgh, nhcPin
USE HexUtil, ONLY : SetEqn, FindCnt

IMPLICIT NONE

INTEGER :: iAsy, iGeo, iPin, jPin, kPin, tPin, ivTyp, iBndy, nBndy
REAL    :: Cnt(2)

INTEGER, POINTER :: cpSlfMPnum(:, :), cpSlfMPidx(:, :, :), PinMap(:) ! Global to Super-Pin
INTEGER, POINTER :: cpSufMPnum(:, :, :), cpSufMPidx(:, :, :, :), cpSufMPsuf(:, :, :, :)

TYPE(Type_HexAsyTypInfo), POINTER :: aInf_Loc
TYPE(Type_HexCmfdPin),    POINTER :: cPin_Loc
! ----------------------------------------------------

CALL dmalloc(PinMap, nHexPin)

nhcPin = nHexPin

DO iPin = 1, nhcPin
  PinMap(iPin) = iPin
END DO

ALLOCATE (hcPin (nhcPin))
! ----------------------------------------------------
DO iAsy = 1, nhAsy
  iGeo = hAsy(iAsy)%GeoTyp

  aInf_Loc => hAsyTypInfo(hAsy(iAsy)%AsyTyp)

  cpSufMPnum => aInf_Loc%cpSufMPnum
  cpSlfMPnum => aInf_Loc%cpSlfMPnum
  cpSlfMPidx => aInf_Loc%cpSlfMPidx
  cpSufMPidx => aInf_Loc%cpSufMPidx
  cpSufMPsuf => aInf_Loc%cpSufMPsuf

  !$OMP PARALLEL PRIVATE(iPin, jPin, kPin, ivTyp, nBndy, cPin_Loc, Cnt, iBndy, tPin)
  !$OMP DO SCHEDULE(GUIDED)
  DO iPin = 1, hAsy(iAsy)%nTotPin
    jPin  = iPin + hAsy(iAsy)%PinIdxSt - 1  ! Global Pin Idx
    kPin  = hPinInfo(jPin)%OrdInAsy01       ! Pin Idx in "hAsyTypInfo"
    ivTyp = aInf_Loc%PinVtxTyp(iGeo, kPin)  ! Vtx Typ
    nBndy = mpTypNumNgh(ivTyp)

    cPin_Loc => hcPin(PinMap(jPin))
    ! ----------------------------
    !      1. CP
    ! ----------------------------
    cPin_Loc%nBndy = nBndy
    cPin_Loc%aIdx  = iAsy

    cPin_Loc%BdPts(1, 1:7) = aInf_Loc%mpVtx(1, 1:7, iGeo, kPin) + hAsy(iAsy)%Cnt(1)
    cPin_Loc%BdPts(2, 1:7) = aInf_Loc%mpVtx(2, 1:7, iGeo, kPin) + hAsy(iAsy)%Cnt(2)

    Cnt(1:2) = FindCnt(nBndy, cPin_Loc%BdPts(1:2, 1:nBndy))

    DO iBndy = 1, nBndy
      cPin_Loc%BdLgh(iBndy) = aInf_Loc%mpBndyLgh(iBndy, ivTyp)
      cPin_Loc%BdC2B(iBndy) = aInf_Loc%mpBndyC2B(iBndy, ivTyp)

      cPin_Loc%BdEqn(1:3, iBndy) = SetEqn(cPin_Loc%BdPts(1:2, iBndy), cPin_Loc%BdPts(1:2, iBndy+1), Cnt)
    END DO
    ! ----------------------------
    !      2. MP
    ! ----------------------------
    cPin_Loc%Area = aInf_Loc%mpTypAre(ivTyp)

    cPin_Loc%nmPin   = cpSlfMPnum   (iGeo, kPin)
    cPin_Loc%nBdmPin = cpSufMPnum(:, iGeo, kPin)

    DO jPin = 1, cPin_Loc%nmPin
      tPin = cpSlfMPidx(jPin, iGeo, kPin)

      cPin_Loc%mpIdx(jPin) = aInf_Loc%PinLocIdx(iGeo, tPin) + hAsy(iAsy)%PinIdxSt - 1
    END DO

    DO iBndy = 1, cPin_Loc%nBndy
      DO jPin = 1, cpSufMPnum(iBndy, iGeo, kPin)
        tPin = cpSufMPidx(jPin, iBndy, iGeo, kPin)

        cPin_Loc%BdMPsuf(jPin, iBndy) = cpSufMPsuf(jPin, iBndy, iGeo, kPin)
        cPin_Loc%BdMPidx(jPin, iBndy) = aInf_Loc%PinLocIdx(iGeo, tPin) + hAsy(iAsy)%PinIdxSt - 1 ! Global Pin Idx
      END DO
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
END DO

NULLIFY (PinMap, cpSufMPnum, cpSlfMPnum, cpSlfMPidx, cpSufMPidx, cpSufMPsuf, aInf_Loc, cPin_Loc)
! ----------------------------------------------------

END SUBROUTINE HexSetHcPinSlf_MP
! ------------------------------------------------------------------------------------------------------------
!                                     04. SET : hcPin Ngh
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetHcPinNgh()

USE PARAM,   ONLY : ZERO, VoidCell, RefCell, TRUE, FALSE, HALF
USE HexType, ONLY : Type_HexCmfdPin
USE HexData, ONLY : hcPin, nhcPin, hLgc, cBndyEq, hAsy, hAsyTypInfo, Sq3
USE HexUtil, ONLY : SetEqn, ChkPtEqn, ChkSamePts, SetPtsPri

IMPLICIT NONE

INTEGER :: iPin, jPin, iBndy, jBndy, nNgh, iPri
LOGICAL :: lChk, lChk01, lChk02, lChk03, lChk04, lChk05, lNgh
REAL    :: Cnt(2), PtsSlf(2, 2), PtsNgh(2, 2), EqnSlf(3), Lgh, EqnTst(3)

TYPE(Type_HexCmfdPin), POINTER :: cPin_Loc, jPin_Loc
! ----------------------------------------------------

Cnt    = ZERO
Lgh    = hAsyTypInfo(hAsy(1)%AsyTyp)%pPch
EqnTst = SetEqn([Lgh, ZERO], [Lgh * HALF, -Lgh * HALF * Sq3], Cnt)

!$OMP PARALLEL PRIVATE(cPin_Loc, nNgh, iBndy, EqnSlf, PtsSlf, iPri, lNgh, lChk, lChk01, lChk02, lChk03, lChk04, &
!$OMP                  jPin_Loc, jBndy, PtsNgh, Lgh)
!$OMP DO SCHEDULE(GUIDED)
DO iPin = 1, nhcPin
  cPin_Loc => hcPin(iPin)
  
  nNgh = 0
  ! ----------------------------
  DO iBndy = 1, cPin_Loc%nBndy
    EqnSlf(1:3)    = cPin_Loc%BdEqn(1:3, iBndy)
    PtsSlf(1:2, 1) = cPin_Loc%BdPts(1:2, iBndy)
    PtsSlf(1:2, 2) = cPin_Loc%BdPts(1:2, iBndy + 1)
    
    lChk01 = ChkPtEqn(PtsSlf(1:2, 1), cBndyEq(1:3, 1)) .AND. ChkPtEqn(PtsSlf(1:2, 2), cBndyEq(1:3, 1))
    lChk02 = ChkPtEqn(PtsSlf(1:2, 1), cBndyEq(1:3, 2)) .AND. ChkPtEqn(PtsSlf(1:2, 2), cBndyEq(1:3, 2))
    lChk   = hLgc%l060 .AND. (lChk01 .OR. lChk02)
    ! ----------------------------
    !      1. Core REF bndy
    ! ----------------------------
    IF (lChk .AND. hLgc%lAzmRef) THEN
      nNgh = nNgh + 1
      
      cPin_Loc%NghPin(nNgh) = RefCell
      cPin_Loc%NghBd (nNgh) = iBndy
      cPin_Loc%NghLgh(nNgh) = cPin_Loc%BdLgh(iBndy)
      
      CYCLE
    END IF
    
    IF (lChk .AND. hLgc%lAzmRot) THEN
      IF (ChkPtEqn(PtsSlf(1:2, 1), EqnTst) .OR. ChkPtEqn(PtsSlf(1:2, 2), EqnTst)) THEN
        nNgh = nNgh + 1
        
        cPin_Loc%NghPin(nNgh) = RefCell
        cPin_Loc%NghBd (nNgh) = iBndy
        cPin_Loc%NghLgh(nNgh) = cPin_Loc%BdLgh(iBndy)
        
        CYCLE
      END IF
      
      PtsSlf = CalRotPt(EqnSlf, PtsSlf)
      EqnSlf = SetEqn(PtsSlf(1:2, 1), PtsSlf(1:2, 2), Cnt) ! Cnt is meaningless
    END IF
    
    iPri = SetPtsPri(PtsSlf(1:2, 1:2))
    lNgh = FALSE
    ! ----------------------------
    !      2. CP - CP
    ! ----------------------------
    DO jPin = 1, nhcPin
      IF (iPin .EQ. jPin) CYCLE
      
      jPin_Loc => hcPin(jPin)
      
      DO jBndy = 1, jPin_Loc%nBndy
        PtsNgh(1:2, 1) = jPin_Loc%BdPts(1:2, jBndy)
        PtsNgh(1:2, 2) = jPin_Loc%BdPts(1:2, jBndy+1)
        
        lChk01 = ChkPtEqn(PtsNgh(1:2, 1), EqnSlf(1:3))
        lChk02 = ChkPtEqn(PtsNgh(1:2, 2), EqnSlf(1:3))
        lChk03 = lChk01 .AND. lChk02
        
        IF (.NOT. lChk03) CYCLE
        
        lChk01 = ChkSamePts(PtsSlf(1:2, 1), PtsNgh(1:2, 1))
        lChk02 = ChkSamePts(PtsSlf(1:2, 2), PtsNgh(1:2, 2))
        lChk03 = lChk01 .AND. lChk02
        
        lChk01 = ChkSamePts(PtsSlf(1:2, 1), PtsNgh(1:2, 2))
        lChk02 = ChkSamePts(PtsSlf(1:2, 2), PtsNgh(1:2, 1))
        lChk04 = lChk01 .AND. lChk02
        
        ! Exact Ngh
        IF (lChk03 .OR. lChk04) THEN
          nNgh = nNgh + 1
          lNgh = TRUE
          
          cPin_Loc%NghPin(nNgh) = jPin
          cPin_Loc%NghBd (nNgh) = iBndy
          cPin_Loc%NghSuf(nNgh) = jBndy
          cPin_Loc%NghLgh(nNgh) = cPin_Loc%BdLgh(iBndy)
          
          EXIT
        END IF
        
        Lgh = CalNghSegLgh(PtsSlf, PtsNgh, iPri)
        
        IF (Lgh < 1E-3) CYCLE ! ARTIBRARY
        
        ! Segment Ngh
        nNgh = nNgh + 1
        lNgh = TRUE
        
        cPin_Loc%NghPin(nNgh) = jPin
        cPin_Loc%NghBd (nNgh) = iBndy
        cPin_Loc%NghSuf(nNgh) = jBndy
        cPin_Loc%NghLgh(nNgh) = Lgh
      END DO
    END DO
    
    IF (lNgh) CYCLE
    ! ----------------------------
    !      3. Core VAC bndy
    ! ----------------------------
    nNgh = nNgh + 1
    
    cPin_Loc%NghPin(nNgh) = VoidCell
    cPin_Loc%NghBd (nNgh) = iBndy
    cPin_Loc%NghLgh(nNgh) = cPin_Loc%BdLgh(iBndy)
    
    IF (hLgc%iSym > 3) cPin_Loc%NghPin(nNgh) = RefCell ! Sng Asy / Cel
    IF (hLgc%lRadRef)  cPin_Loc%NghPin(nNgh) = RefCell ! REF on Radial direction
    IF (lChk)          cPin_Loc%NghPin(nNgh) = RefCell ! Ngh Pin = Self Pin with Azm ROT
  END DO
  ! ----------------------------
  
  cPin_Loc%nNgh = nNgh
END DO
!$OMP END DO
!$OMP END PARALLEL

NULLIFY (cPin_Loc, jPin_Loc)
! ----------------------------------------------------

END SUBROUTINE HexSetHcPinNgh
! ------------------------------------------------------------------------------------------------------------
!                                     05. CAL : Rot Pt
! ------------------------------------------------------------------------------------------------------------
FUNCTION CalRotPt(Eqn, Pts)

USE PARAM,   ONLY : ZERO
USE HexData, ONLY : PI_3
USE HexUtil, ONLY : ChkSameVal, RotPt

IMPLICIT NONE

REAL    :: Eqn(3), Pts(2, 2), Theta, CalRotPt(2, 2)
LOGICAL :: lChk
! ----------------------------------------------------

lChk = ChkSameVal(Eqn(1), ZERO)

IF (lChk) THEN
  Theta = -PI_3
ELSE
  Theta =  PI_3
END IF

CalRotPt(1:2, 1) = RotPt(Pts(1:2, 1), Theta)
CalRotPt(1:2, 2) = RotPt(Pts(1:2, 2), Theta)
! ----------------------------------------------------

END FUNCTION CalRotPt
! ------------------------------------------------------------------------------------------------------------
!                                     06. CAL : Ngh Seg Lgh
! ------------------------------------------------------------------------------------------------------------
FUNCTION CalNghSegLgh(PtsSlf, PtsNgh, iPri)

USE PARAM,   ONLY : ZERO
USE HexUtil, ONLY : FindPtLgh

IMPLICIT NONE

REAL :: CalNghSegLgh, PtsSlf(2, 2), PtsNgh(2, 2), Bndy(2, 2), Ngh(2, 2), Seg(2, 2)
INTEGER :: iPri
! ----------------------------------------------------

CalNghSegLgh = ZERO

IF (PtsSlf(iPri, 1) < PtsSlf(iPri, 2)) THEN
  Bndy(1:2, 1) = PtsSlf(1:2, 1)
  Bndy(1:2, 2) = PtsSlf(1:2, 2)
ELSE
  Bndy(1:2, 1) = PtsSlf(1:2, 2)
  Bndy(1:2, 2) = PtsSlf(1:2, 1)
END IF

IF (PtsNgh(iPri, 1) < PtsNgh(iPri, 2)) THEN
  Ngh(1:2, 1) = PtsNgh(1:2, 1)
  Ngh(1:2, 2) = PtsNgh(1:2, 2)
ELSE
  Ngh(1:2, 1) = PtsNgh(1:2, 2)
  Ngh(1:2, 2) = PtsNgh(1:2, 1)
END IF

IF (Ngh(iPri, 1) < Bndy(iPri, 1)) THEN
  Seg(1:2, 1) = Bndy(1:2, 1)
ELSE
  Seg(1:2, 1) = Ngh(1:2, 1)
END IF

IF (Ngh(iPri, 2) > Bndy(iPri, 2)) THEN
  Seg(1:2, 2) = Bndy(1:2, 2)
ELSE
  Seg(1:2, 2) = Ngh(1:2, 2)
END IF

IF (Seg(iPri, 1) > Seg(iPri, 2)) RETURN

CalNghSegLgh = FindPtLgh(Seg(1:2, 1), Seg(1:2, 2))
! ----------------------------------------------------

END FUNCTION CalNghSegLgh
! ------------------------------------------------------------------------------------------------------------
!                                     07. CnP : Super Pin
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexCnPSuperPin(CoreInfo, superPin, nxy, lSuperpin)

USE allocs
USE PARAM,          ONLY : TRUE, FALSE
USE TYPEDEF,        ONLY : CoreInfo_Type, Cell_Type, Pin_Type
USE HexType,        ONLY : Type_HexCmfdPin
USE HexData,        ONLY : nhcPin, hcPin
USE TYPEDEF_COMMON, ONLY : superPin_Type

IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo

TYPE(superPin_Type), POINTER :: superPin(:)
TYPE(Cell_Type),     POINTER :: Cell(:)
TYPE(Pin_Type),      POINTER :: Pin(:)

INTEGER :: nxy, ixy, jxy, iz, iPin, iCel
LOGICAL :: lSuperpin

TYPE(superPin_Type),   POINTER :: sPin_Loc
TYPE(Type_HexCmfdPin), POINTER :: cPin_Loc
! ----------------------------------------------------

nxy       = nhcPin
lSuperpin = TRUE

ALLOCATE (superPin (nxy))
! ----------------------------------------------------
!               01. CP
! ----------------------------------------------------
DO ixy = 1, nxy
  sPin_Loc => superPin(ixy)
  cPin_Loc => hcPin(ixy)
  
  ! Garbage
  sPin_Loc%nx = 0
  sPin_Loc%ny = 0
  sPin_Loc%ix = 0
  sPin_Loc%iy = 0
  
  ! CnP
  sPin_Loc%nxy = cPin_Loc%nmPin
  
  CALL dmalloc(sPin_Loc%pin,             3)
  CALL dmalloc(sPin_Loc%neighidx,       15)
  CALL dmalloc(sPin_Loc%neighsurfidx,   15)
  CALL dmalloc(sPin_Loc%BDLength,        6)
  CALL dmalloc(sPin_Loc%Center2SurfaceL, 6)
  
  ! SELF
  sPin_Loc%pin  = cPin_Loc%mpIdx  ! MOC Pin
  sPin_Loc%Area = cPin_Loc%Area
  
  sPin_Loc%BdLength(1:6) = cPin_Loc%BdLgh(1:6)
  
  sPin_Loc%Center2SurfaceL(1:6) = cPin_Loc%BdC2B(1:6)
  
  ! MP
  sPin_Loc%nBdmPin = cPin_Loc%nBdmPin
  sPin_Loc%BdMPidx = cPin_Loc%BdMPidx ! MOC Pin
  sPin_Loc%BdMPsuf = cPin_Loc%BdMPsuf ! MOC Pin
  
  ! CnP
  sPin_Loc%nNgh         = cPin_Loc%nNgh
  sPin_Loc%NeighIdx     = cPin_Loc%NghPin
  sPin_Loc%NeighSurfIdx = cPin_Loc%NghSuf
  sPin_Loc%NghBd        = cPin_Loc%NghBd
  sPin_Loc%NghLgh       = cPin_Loc%NghLgh
END DO
! ----------------------------------------------------
!               02. SET : Fuel Dat
! ----------------------------------------------------
Cell => CoreInfo%CellInfo
Pin  => CoreInfo%Pin

DO ixy = 1, nxy
  CALL dmalloc(superPin(ixy)%lFuel, CoreInfo%nz)
  
  DO iz = 1, CoreInfo%nz
    DO jxy = 1, superPin(ixy)%nxy
      iPin = superPin(ixy)%pin(jxy)
      iCel = Pin(iPin)%Cell(iz)

      IF (Cell(iCel)%lFuel) THEN
        superPin(ixy)%lFuel(iz) = TRUE

        EXIT
      END IF
    END DO
  END DO

  DO jxy = 1, superPin(ixy)%nxy
    iPin = superPin(ixy)%pin(jxy)

    IF (Pin(iPin)%lFuel) THEN
      superPin(ixy)%iFuelPin = ipin
      EXIT
    END IF
  END DO
END DO

NULLIFY (sPin_Loc, cPin_Loc, Pin, Cell)
! ----------------------------------------------------

END SUBROUTINE HexCnPSuperPin
! ------------------------------------------------------------------------------------------------------------
!                                     08. SET : Global Pin
! ------------------------------------------------------------------------------------------------------------
#ifdef __INTEL_MKL
SUBROUTINE HexSetGlobalPin(CoreInfo, FmInfo)

USE allocs
USE MKL_3D
USE PARAM,   ONLY : TRUE, FALSE
USE HexData, ONLY : nhcPin
USE TYPEDEF, ONLY : CoreInfo_Type, FmInfo_Type, FxrInfo_Type, Pin_Type, Cell_Type

IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(FmInfo_Type)   :: FmInfo

TYPE(superPin_Type), POINTER :: superPin(:)
TYPE(FxrInfo_Type),  POINTER :: Fxr(:, :)
TYPE(Pin_Type),      POINTER :: Pin(:)
TYPE(Cell_Type),     POINTER :: Cell(:)

INTEGER :: izf, iz, iPin, jPin, kPin, iCel, iFXR, jFXR, nzCMFD

REAL, POINTER :: hzfm(:)
! ----------------------------------------------------

nzCMFD    = mklGeom%nzCMFD
hzfm     => mklGeom%hzfm
superPin => mklGeom%superPin
Fxr      => FmInfo%Fxr
Pin      => CoreInfo%Pin
Cell     => CoreInfo%CellInfo

CALL dmalloc (mklGeom%pinMap,           nhcPin)
CALL dmalloc0(mklGeom%pinMapRev, -1,    nhcPin)
CALL dmalloc (mklGeom%PinVolFm,         nhcPin, nzCMFD)
CALL dmalloc (mklGeom%lRefPin,          nhcPin)
CALL dmalloc (mklGeom%lRefCell, nzCMFD, nhcPin)
CALL dmalloc (mklGeom%lH2OCell, nzCMFD, nhcPin)

! Pin Map
mklGeom%pinMapRev(-1) = -1
mklGeom%pinMapRev (0) =  0

DO iPin = 1, nhcPin
  mklGeom%pinMap   (iPin) = iPin
  mklGeom%pinMapRev(iPin) = iPin
END DO

! Pin Vol Fm
DO izf = 1, nzCMFD
  iz = mklGeom%planeMap(izf)

  DO iPin = 1, nhcPin
    mklGeom%PinVolFm(iPin, izf) = superPin(iPin)%Area * hzfm(izf)
  END DO
END DO

! Ref Dat
DO iPin = 1, nhcPin
  mklGeom%lRefPin(iPin) = .NOT. ANY(superPin(iPin)%lFuel)

  DO izf = 1, nzCMFD
    iz = mklGeom%planeMap(izf)

    mklGeom%lRefCell(izf, iPin) = .NOT. superPin(iPin)%lFuel(iz)
  END DO
END DO

! H2O Dat
DO iPin = 1, nhcPin
  DO izf = 1, nzCMFD
    iz = mklGeom%planeMap(izf)

    mklGeom%lH2OCell(izf, iPin) = TRUE

    DO jPin = 1, superPin(iPin)%nxy
      kPin = superPin(iPin)%pin(jPin)
      iCel = Pin(kPin)%Cell(iz)

      DO iFXR = 1, Cell(iCel)%nFxr
        jfxr = Pin(kPin)%FxrIdxSt + iFXR - 1

        IF (.NOT. Fxr(jfxr, iz)%lH2O) mklGeom%lH2OCell(izf, iPin) = FALSE
      END DO
    END DO
  END DO
END DO
! ----------------------------------------------------

NULLIFY (superPin, hzfm, Pin, Cell, FXR)

END SUBROUTINE HexSetGlobalPin
#endif
! ------------------------------------------------------------------------------------------------------------
!                                     09. SET : CMFD Sng Pin
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetHcPin_Sng()

USE PARAM,   ONLY : HALF
USE geom,    ONLY : Core, Pin
USE HexData, ONLY : hCelBss

IMPLICIT NONE
! ----------------------------------------------------

Core%Pin => Pin

ALLOCATE (Core%Pin(1)%NeighIdx        (6))
ALLOCATE (Core%Pin(1)%NeighSurfIdx    (6))
ALLOCATE (Core%Pin(1)%Center2SurfaceL (6))
ALLOCATE (Core%Pin(1)%BdLength        (6))

Core%Pin(1)%NeighIdx        = -1
Core%Pin(1)%NeighSurfIdx    = -1
Core%Pin(1)%Center2SurfaceL = hCelBss(1)%pF2F * HALF ! # of hCelBss must be 1 for Sng Cel
Core%Pin(1)%BdLength        = hCelBss(1)%pPch
! ----------------------------------------------------

END SUBROUTINE HexSetHcPin_Sng
! ------------------------------------------------------------------------------------------------------------

END MODULE HexCMFDConst
