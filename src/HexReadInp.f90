SUBROUTINE HexReadInp(indev)

USE param,      ONLY : ONE, BLANK, DOT
USE geom,       ONLY : nZ, hz, hzInv
USE ioutil,     ONLY : toupper, ifnumeric, icolfield
USE RAYS,       ONLY : RayInfo
USE inputcards, ONLY : oneline, probe, FindCardId
USE PE_MOD,     ONLY : PE

USE HexInpCard
USE HexChk
USE HexData,      ONLY : aoF2F, aoPch, nAzmAng, nPolAng, Del_Inp, Sq3Inv
USE HexCelBssGap, ONLY : HexSetGapCelBss
USE HexCelBssRod, ONLY : HexSetRodCelBss

IMPLICIT NONE

CHARACTER(15)  :: cardname, astring
CHARACTER(256) :: dataline

INTEGER, PARAMETER :: idblock = 6

INTEGER :: ist, iz, indev, idcard, nLineField
LOGICAL :: Master
! ----------------------------------------------------

RayInfo%nAziAngle = (RayInfo%nAziAngle / 6) * 6

nAzmAng = RayInfo%nAziAngle
nPolAng = RayInfo%nPolarAngle
Del_Inp = RayInfo%Del
Master  = PE%master

CALL HexInitInp
! ----------------------------------------------------
DO WHILE (TRUE)
  READ (indev, '(A256)') oneline
  
  IF (Master) CALL message(io8, FALSE, FALSE, oneline) ! ECHO
  
  IF (probe   .EQ. BANG)  CYCLE
  IF (probe   .EQ. POUND) CYCLE
  IF (oneline .EQ. BLANK) CYCLE
  IF (IFnumeric(oneline)) CYCLE
  
  IF (probe .EQ. DOT)   EXIT
  IF (probe .EQ. SLASH) EXIT
  
  READ (oneline, *) cardname
  CALL toupper(cardname)
  
  idcard     = FindCardId(idblock, cardname)
  nLineField = nfields(oneline)-1
  
  IF (idcard .EQ. 0) EXIT
  
  ist = icolfield(oneline, 2)
  dataline = oneline(ist:256)
  ! --------------------------------------------------
  SELECT CASE (idcard)
    ! PITCH
    CASE (3)
      IF (nLineField .NE. 1) CALL terminate("PITCH")
      
      READ (oneline, *) ASTRING, aoF2F
      
      aoPch = aoF2F * Sq3Inv
    
    ! AX_MESH
    CASE (4)
      CALL dmalloc(hz,    nz)
      CALL dmalloc(HzInv, nz)
      
      READ (oneline, *) astring, (hz(iz), iz = 1, nz)
      
      HzInv = ONE / Hz
      
    ! Other
    CASE (5);  CALL HexRead_Albedo (dataline)
    CASE (7);  CALL HexRead_Cel    (dataline)
    CASE (9);  CALL HexRead_Pin    (dataline)
    CASE (10); CALL HexRead_Asy    (indev, dataline)
    CASE (11); CALL HexRead_Core   (indev, dataline)
    CASE (13); CALL HexRead_GapCel (dataline)
    CASE (14); CALL HexRead_GapPin (dataline)
    CASE (21); CALL HexRead_Vss    (dataline)
    CASE (22); CALL HexRead_Vyg    (dataline)
    CASE (23); CALL HexRead_Opt    (dataline)
  END SELECT
END DO
! ----------------------------------------------------
BACKSPACE (indev)
BACKSPACE (io8)

CALL HexChkInp
CALL HexSetRodCelBss
CALL HexSetGapCelBss
CALL HexConnCelBss
CALL HexChkCelBss
! ----------------------------------------------------

END SUBROUTINE HexReadInp