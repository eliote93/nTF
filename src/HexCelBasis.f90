SUBROUTINE HexConnCelBss()

USE PARAM,   ONLY : epsm7, TRUE, FALSE
USE geom,    ONLY : nAsyType0
USE ioutil,  ONLY : terminate
USE HexType, ONLY : Type_HexGapCelBss
USE HexData, ONLY : hCelBss, gCelBss, gCel, ncBss, ngBss, hAsyTypInfo, hLgc
USE HexUtil, ONLY : ChkSameVal

IMPLICIT NONE

INTEGER :: iBss, jBss, iaTyp
LOGICAL :: l01, l02, l03

TYPE(Type_HexGapCelBss) :: gCelBss_Tmp(100)
! ----------------------------------------------------

! ----------------------------------------------------
!               01. CASE : Sng Cel
! ----------------------------------------------------
IF (hLgc%lSngCel) THEN
  DO iaTyp = 1, nAsyType0
    IF (hAsyTypInfo(iaTyp)%luse) EXIT
  END DO
  
  hAsyTypInfo(iaTyp)%iBss = 1 ! # of hCel Bss must be 1 for Sng Cel
  
  hCelBss(1)%iaTyp = iaTyp
  
  RETURN
END IF
! ----------------------------------------------------
!               02. CONN : Rod and Gap Cel Basis
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
!               03. CONN : Asy Cel Basis
! ----------------------------------------------------
DO iBss = 1, ncBss
  DO iaTyp = 1, nAsyType0
    IF (.NOT. hAsyTypInfo(iaTyp)%luse) CYCLE
    
    l01 = ChkSameVal(hAsyTypInfo(iaTyp)%aiF2F,     hCelBss(iBss)%aiF2F)
    l02 = ChkSameVal(hAsyTypInfo(iaTyp)% pF2F,     hCelBss(iBss)% pF2F)
    l03 =            hAsyTypInfo(iaTyp)% nPin .EQ. hCelBss(iBss)% nPin
    
    IF (.NOT.(l01 .AND. l02 .AND. l03)) CYCLE
    
    hAsyTypInfo(iaTyp)%iBss = iBss
    
    hCelBss(iBss)%iaTyp = iaTyp
  END DO
END DO

DO iBss = 1, ncBss
  IF (hCelBss(iBss)%igBss .EQ. 0) CALL terminate("UNPAIRED ROD CEL BASIS EXISTS")
  IF (hCelBss(iBss)%iaTyp .EQ. 0) CALL terminate("UNPAIRED ROD CEL BASIS EXISTS")
END DO

DO iaTyp = 1, nAsyType0
  IF (.NOT. hAsyTypInfo(iaTyp)%luse) CYCLE
  
  IF (hAsyTypInfo(iaTyp)%iBss .EQ. 0) CALL terminate("UNPAIRED ROD CEL BASIS EXISTS")
END DO
! ----------------------------------------------------

END SUBROUTINE HexConnCelBss
! ------------------------------------------------------------------------------------------------------------