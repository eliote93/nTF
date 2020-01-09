SUBROUTINE HexConnCelBss()


USE PARAM,   ONLY : epsm7, TRUE, FALSE
USE geom,    ONLY : nGapType, nAsyType0
USE ioutil,  ONLY : terminate
USE HexType, ONLY : Type_HexGapCelBss
USE HexData, ONLY : hCelBss, gCelBss, gCel, ncBss, ngBss, hAsyTypInfo, hLgc
USE HexUtil, ONLY : ChkSameVal

IMPLICIT NONE

INTEGER :: iBss, jBss, iAsyTyp
LOGICAL :: l01, l02, l03

TYPE(Type_HexGapCelBss) :: gCelBss_Tmp(100)
! ----------------------------------------------------

hAsyTypInfo(1)%iBss = 1
hCelBss(1)%iaTyp    = 1

IF (hLgc%lSngCel) RETURN

! ----------------------------------------------------
!               01. CONN : Rod and Gap Cel Basis
! ----------------------------------------------------
DO iBss = 1, ncBss
  DO jBss = 1, ngBss
    l01 = ChkSameVal(hCelBss(iBss)%aiF2F,     gCelBss(jBss)%aiF2F)
    l02 = ChkSameVal(hCelBss(iBss)% pF2F,     gCelBss(jBss)% pF2F)
    l03 =            hCelBss(iBss)% nPin .EQ. gCelBss(jBss)% nPin
    
    IF (.NOT.(l01 .AND. l02 .AND. l03)) CYCLE
    
    hCelBss(iBss)%igBss = jBss
    gCelBss(jBss)%icBss = iBss
  END DO
END DO
! ----------------------------------------------------
!               02. CONN : Asy Cel Basis
! ----------------------------------------------------
DO iBss = 1, ncBss
  DO iAsyTyp = 1, nAsyType0
    l01 = ChkSameVal(hAsyTypInfo(iAsyTyp)%aiF2F,     hCelBss(iBss)%aiF2F)
    l02 = ChkSameVal(hAsyTypInfo(iAsyTyp)% pF2F,     hCelBss(iBss)% pF2F)
    l03 =            hAsyTypInfo(iAsyTyp)% nPin .EQ. hCelBss(iBss)% nPin
    
    IF (.NOT.(l01 .AND. l02 .AND. l03)) CYCLE
    
    hAsyTypInfo(iAsyTyp)%iBss = iBss
    hCelBss(iBss)%iaTyp       = iAsyTyp
  END DO
END DO

DO iBss = 1, ncBss
  IF (hCelBss(iBss)%igBss .EQ. 0) CALL terminate("UNPAIRED ROD CEL BASIS EXISTS")
  IF (hCelBss(iBss)%iaTyp .EQ. 0) CALL terminate("UNPAIRED ROD CEL BASIS EXISTS")
END DO

DO iAsyTyp = 1, nAsyType0
  IF (hAsyTypInfo(iAsyTyp)%iBss .EQ. 0) CALL terminate("UNPAIRED ROD CEL BASIS EXISTS")
END DO
! ----------------------------------------------------

END SUBROUTINE HexConnCelBss
! ------------------------------------------------------------------------------------------------------------