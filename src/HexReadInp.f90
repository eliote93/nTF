SUBROUTINE HexReadInp(indev)

USE param,      ONLY : ONE, BLANK, DOT
USE geom,       ONLY : nZ, hz, hzInv
USE ioutil,     ONLY : toupper, ifnumeric, icolfield
USE RAYS,       ONLY : RayInfo
USE inputcards, ONLY : oneline, probe, FindCardId
USE PE_MOD,     ONLY : PE
USE HexData,    ONLY : aoF2F, aoPch, nAzmAng, nPolAng, Del_Inp, Sq3Inv
USE HexInpCard

IMPLICIT NONE

CHARACTER(15)  :: cardname, astring
CHARACTER(256) :: dataline

INTEGER, PARAMETER :: idblock = 6

INTEGER :: i, k, indev, idcard, nLineField
LOGICAL :: Master
! ----------------------------------------------------

RayInfo%nAziAngle = (RayInfo%nAziAngle / 6) * 6

nAzmAng = RayInfo%nAziAngle
nPolAng = RayInfo%nPolarAngle
Del_Inp = RayInfo%Del
Master  = PE%master

CALL HexInitInp
! ----------------------------------------------------
DO WHILE(TRUE)
  READ (indev,'(a256)') oneline
  
  IF (Master) CALL message(io8, FALSE, FALSE, oneline) ! ECHO
  ! ----------------------------------------------------
  !               01. CHK : Valid Input
  ! ----------------------------------------------------
  IF(probe.eq.BANG) cycle;     IF(probe.eq.POUND) cycle
  IF(oneline.eq.BLANK) cycle;  IF(IFnumeric(oneline)) cycle
  IF(probe.eq.DOT) exit;       IF(probe.eq.SLASH) exit
  
  READ (oneline, *) cardname
  CALL toupper(cardname)
  
  idcard     = FindCardId(idblock, cardname)
  nLineField = nfields(oneline)-1
  
  IF(idcard .eq. 0) EXIT
  ! ----------------------------------------------------
  !               02. READ
  ! ----------------------------------------------------  
  SELECT CASE(idcard)
    ! ----------------------------
    !      3. PITCH
    ! ----------------------------
    CASE(3)
      IF(nLineField .ne. 1) CALL terminate("PITCH")
      
      READ (oneline,*) ASTRING, aoF2F
      
      aoPch = aoF2F * Sq3Inv
    ! ----------------------------
    !      4. AX_MESH
    ! ----------------------------
    CASE(4)
      CALL dmalloc(hz,    nz)
      CALL dmalloc(HzInv, nz)
      
      READ (oneline, *) astring, (hz(k), k = 1, nz)
      HzInv = 1. / Hz
    ! ----------------------------
    !      5. ALBEDO
    ! ----------------------------
    CASE(5)
      i = icolfield(oneline,2)
      dataline = oneline(i:256)
      
      CALL HexRead_Albedo(dataline)
    ! ----------------------------
    !      7. CEL
    ! ----------------------------
    CASE(7)
      i = icolfield(oneline,2)
      dataline = oneline(i:256)
      
      CALL HexRead_Cel(dataline)
    ! ----------------------------
    !      8. PIN
    ! ----------------------------
    CASE(9)
      i=icolfield(oneline,2)
      dataline = oneline(i:256)
      
      CALL HexRead_Pin(dataline)
    ! ----------------------------
    !      9. ASSEMBLY
    ! ----------------------------
    CASE(10)
      i=icolfield(oneline,2)
      dataline = oneline(i:256)
      
      CALL HexRead_Asy(indev,dataline)
    ! ----------------------------
    !     11. RAD_CONF
    ! ----------------------------
    CASE(11)
      i=icolfield(oneline,2)
      dataline = oneline(i:256)
      
      CALL HexRead_Core(indev,dataline)
    ! ----------------------------
    !     13. GAP CEL
    ! ----------------------------
    CASE(13)
      i=icolfield(oneline,2)
      dataline = oneline(i:256)
      
      CALL HexRead_GapCel(dataline)
    ! ----------------------------
    !     14. GAP Pin
    ! ----------------------------
    CASE(14)
      i=icolfield(oneline,2)
      dataline = oneline(i:256)
      
      CALL HexRead_GapPin(dataline)
    ! ----------------------------
    !     21. Vessel
    ! ----------------------------
    CASE(21)
      i=icolfield(oneline,2)
      dataline = oneline(i:256)
      
      CALL HexRead_Vss(dataline)
    ! ----------------------------
    !     22. Vygordka
    ! ----------------------------
    CASE(22)
      i=icolfield(oneline,2)
      dataline = oneline(i:256)
      
      CALL HexRead_Vyg(dataline)
    ! ----------------------------
    !     23. Hex Opt
    ! ----------------------------
    CASE(23)
      i=icolfield(oneline,2)
      dataline = oneline(i:256)
      
      CALL HexRead_Opt(dataline)
  END SELECT
END DO

BACKSPACE(indev)
! ----------------------------------------------------
!               03. SET : Cel Geo Basis
! ----------------------------------------------------
CALL HexChkInp
CALL HexSetRodCelBss
CALL HexSetGapCelBss
CALL HexConnCelBss
CALL HexChkCelBss
! ----------------------------------------------------

END SUBROUTINE HexReadInp