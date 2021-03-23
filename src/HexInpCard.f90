MODULE HexInpCard

USE allocs
USE param,      ONLY : TRUE, FALSE, BANG, POUND, SLASH
USE files,      ONLY : io8
USE ioutil,     ONLY : nfields, fndchara, message, terminate
USE inputcards, ONLY : oneline, probe
USE HexData,    ONLY : hLgc

IMPLICIT NONE

CONTAINS

! ------------------------------------------------------------------------------------------------------------
!                                     01. HEX INIT : INP
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexInitInp()

USE allocs
USE geom,    ONLY : nZ, nCellType, nPinType, nGapType, nGapPinType, nAsyType0, nVssTyp
USE HexData, ONLY : RodPin, GapPin, hCel, gCel, hAsy, hAsyTypInfo, hVss

IMPLICIT NONE

INTEGER :: iPin
! ----------------------------------------------------

IF (hLgc%lvyg) THEN
  nGapPinType = nGapPinType + 2
  nGapType    = nGapType    + 2
END IF

ALLOCATE (hCel        (nCellType    * (nVssTyp + 1)))
ALLOCATE (gCel        (nGapType     * (nVssTyp + 1)))
ALLOCATE (RodPin      (nPinType     * (nVssTyp + 1)))
ALLOCATE (GapPin      (nGapPinType  * (nVssTyp + 1)))
ALLOCATE (hAsyTypInfo (nAsyType0))
ALLOCATE (hVss        (nVssTyp))

DO iPin = 1, nPinType
  CALL dmalloc(RodPin(iPin)%iCel, nZ)
    
  RodPin(iPin)%lRod = TRUE
END DO

DO iPin = 1, nGapPinType
  CALL dmalloc(GapPin(iPin)%iCel, nZ)
  
  GapPin(iPin)%lGap = TRUE
END DO
! ----------------------------------------------------

END SUBROUTINE HexInitInp
! ------------------------------------------------------------------------------------------------------------
!                                     02. HEX READ : Cel
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexRead_Cel(dataline0)

USE cntl,    ONLY : nTracerCntl
USE ioutil,  ONLY : nfieldto
USE HexType, ONLY : Type_HexRodCel, nMaxFXR
USE HexData, ONLY : hCel, aoF2F, Sq3Inv

IMPLICIT NONE

CHARACTER(256),INTENT(IN) :: dataline0
CHARACTER(256) :: dataline

INTEGER :: i, iCel, nSpt, nDataField, nData, mData, tData
INTEGER :: ipos(100), nDiv(100)

REAL :: RR(100)

TYPE(Type_HexRodCel), POINTER :: hCel_Loc
! ----------------------------------------------------

dataline = dataline0

READ (dataline,*) icel
hCel_Loc => hCel(iCel)

nDataField = len_trim(dataline)

CALL fndchara(dataline, ipos, nspt, SLASH)

IF (nSpt .NE. 5)          CALL terminate("WRONG [CELL] - ROD CEL SLASH")
IF (hCel_Loc%nFXR .NE. 0) CALL terminate("WRONG [CELL] - OVERLAPPED CEL INPUT")

nData = nfieldto(dataline, SLASH) ! Number of FXR including Mod
hCel_Loc%nFXR = nData

mData = nfieldto(dataline(ipos(2)+1:nDataField), SLASH)
tData = nfieldto(dataline(ipos(1)+1:nDataField), SLASH)

IF (nData .NE. mData)        CALL terminate("WRONG [CELL] - # of FXRS")
IF (nData .NE. tData)        CALL terminate("WRONG [CELL] - # of FXRS")
IF (hCel_Loc%nFXR > nMaxFXR) CALL terminate("WRONG [CELL] - NMAXFXR")

READ (dataline(ipos(5)+1:nDataField),*) hCel_Loc%aiF2F
READ (dataline(ipos(4)+1:nDataField),*) hCel_Loc%pF2F
READ (dataline(ipos(3)+1:nDataField),*) hCel_Loc%nPin
READ (dataline(ipos(2)+1:nDataField),*) (hCel_Loc%xDiv(nData - i + 1),i=1,nData)
READ (dataline(ipos(1)+1:nDataField),*) (hCel_Loc%xMix(nData - i + 1),i=1,nData)
READ (dataline,*) (RR(i),i=1,ndata)

IF (hCel_Loc%xDiv(1).NE.1 .AND. nTracerCntl%lxslib) THEN
  hCel_Loc%xDiv(1) = 1
  
  WRITE (*,*) "LAST # of FXR DIVISION MUST BE 1"
END IF

DO i = 1, nData - 1
  hCel_Loc%xRad(nData - i + 1) = RR(i+1)
END DO

hCel_Loc%xRad(1) = hCel_Loc%pF2F * Sq3Inv

NULLIFY (hCel_Loc)
! ----------------------------------------------------

END SUBROUTINE HexRead_Cel
! ------------------------------------------------------------------------------------------------------------
!                                     03. HEX READ : Gap Cel
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexRead_GapCel(dataline0)

USE param,   ONLY : HALF
USE ioutil,  ONLY : nfieldto
USE HexType, ONLY : Type_HexGapCel, nMaxFXR
USE HexData, ONLY : gCel, aoF2F

IMPLICIT NONE

CHARACTER(256), INTENT(IN) :: dataline0
CHARACTER(256) :: dataline

INTEGER :: i, iCel
INTEGER :: nSpt, nDataField, nData
INTEGER :: ipos(100)

REAL :: RR(100), ii(7)

TYPE(Type_HexGapCel), POINTER :: gCel_Loc
! ----------------------------------------------------

dataline = dataline0

READ (dataline,*) iCel

nDataField = len_trim(dataline)
CALL fndchara(dataline, ipos, nSpt, SLASH)

gCel_Loc => gCel(iCel)

IF (gCel_Loc%nFXR .NE. 0) CALL terminate("WRONG [GAP_CELL] - OVERLAPPED GAP CEL INPUT")
! ----------------------------------------------------
!               01. CASE : One FXR
! ----------------------------------------------------
SELECT CASE (nSpt)
CASE (4)
  gCel_Loc%nFXR = 1
  
  READ (dataline(ipos(4)+1:nDataField),*) gCel_Loc%aiF2F
  READ (dataline(ipos(3)+1:nDataField),*) gCel_Loc%pF2F
  READ (dataline(ipos(2)+1:nDataField),*) gCel_Loc%nPin
  READ (dataline(ipos(1)+1:nDataField),*) gCel_Loc%xDiv(1)
  READ (dataline,*) (ii(i),i=1,2)
  
  gCel_Loc%xMix(1) = ii(2)
  gCel_Loc%xHgt(1) = (aoF2F - gCel_Loc%aiF2F) * HALF
! ----------------------------------------------------
!               02. CASE : Multi FXR
! ----------------------------------------------------
CASE (5)
  ndata = nfieldto(dataline, SLASH) ! Number of FXR including Mod
  gCel_Loc%nFXR = ndata
  
  IF (gCel_Loc%nFXR > nMaxFXR) CALL terminate("WRONG [GAP_CELL] - NMAXFXR")
  
  READ (dataline(ipos(5)+1:nDataField),*) gCel_Loc%aiF2F
  READ (dataline(ipos(4)+1:nDataField),*) gCel_Loc%pF2F
  READ (dataline(ipos(3)+1:nDataField),*) gCel_Loc%nPin
  
  READ (dataline(ipos(2)+1:nDataField),*) (gCel_Loc%xDiv(ndata - i + 1),i=1,ndata)
  READ (dataline(ipos(1)+1:nDataField),*) (gCel_Loc%xMix(ndata - i + 1),i=1,ndata)
  READ (dataline,*) (RR(i),i=1,ndata)
  
  DO i = 1, ndata - 1
    gCel_Loc%xHgt(ndata - i + 1) = RR(i+1)
  END DO
  
  gCel_Loc%xHgt(1) = (aoF2F - gCel_Loc%aiF2F) * HALF
CASE DEFAULT
  CALL terminate("WRONG [GAP_CELL] - OVERLAPPED PIN INPUT")
END SELECT

NULLIFY (gCel_Loc)
! ----------------------------------------------------

END SUBROUTINE HexRead_GapCel
! ------------------------------------------------------------------------------------------------------------
!                                     04. HEX READ : Pin
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexRead_Pin(dataline0)

USE geom,    ONLY : nz, nCellType
USE HexData, ONLY : RodPin

IMPLICIT NONE

CHARACTER(256), INTENT(IN) :: dataline0
CHARACTER(256) :: dataline

INTEGER :: iz, iPin
INTEGER :: ii(300)
! ----------------------------------------------------

dataline = dataline0

READ (dataline,*) iPin
READ (dataline,*) (ii(iz),iz = 1, nZ + 1)

IF (RodPin(iPin)%iCel(1) .NE. 0) CALL terminate("WRONG [PIN] - OVERLAPPED GAP PIN INPUT")

DO iz = 1, nZ
  RodPin(iPin)%iCel(iz) = ii(iz + 1)
END DO
! ----------------------------------------------------

END SUBROUTINE HexRead_Pin
! ------------------------------------------------------------------------------------------------------------
!                                     05. HEX READ : Gap Pin
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexRead_GapPin(dataline0)

USE geom,    ONLY : nz
USE HexData, ONLY : GapPin

IMPLICIT NONE

CHARACTER(256), INTENT(IN) :: dataline0
CHARACTER(256) :: dataline

INTEGER :: iz, iPin
INTEGER :: ii(300)
! ----------------------------------------------------

dataline = dataline0

READ (dataline,*) iPin
READ (dataline,*) (ii(iz),iz = 1, nZ + 1)

IF (GapPin(iPin)%iCel(1) .NE. 0) CALL terminate("WRONG [GAP_PIN] - OVERLAPPED GAP PIN INPUT")

DO iz = 1, nZ
  GapPin(iPin)%iCel(iz) = ii(iz + 1)
END DO
! ----------------------------------------------------

END SUBROUTINE HexRead_GapPin
! ------------------------------------------------------------------------------------------------------------
!                                     06. HEX READ : Asy
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexRead_Asy(indev, dataline0)

USE geom,    ONLY : nZ, nPinType
USE HexType, ONLY : Type_HexAsyTypInfo
USE HexData, ONLY : hAsyTypInfo, RodPin, hCel, aoF2F, Sq3Inv
USE PE_MOD,  ONLY : PE

IMPLICIT NONE

INTEGER, INTENT(IN) :: indev 
CHARACTER(256), INTENT(IN) :: dataline0
CHARACTER(256) :: dataline

INTEGER :: nSpt, nDataField, nData
INTEGER :: ipos(100), ii(100)
 
INTEGER :: i, j, k, n, iAsy, jfr, jto, Azm
INTEGER :: iCel, iCase

LOGICAL :: Master

TYPE(Type_HexAsyTypInfo), POINTER :: aInf_Loc
! ----------------------------------------------------

Master   = PE%master
dataline = dataline0
! ----------------------------------------------------
!               01. READ : Basic Data
! ----------------------------------------------------
READ (dataline,*) iAsy, Azm

IF (Azm .NE. 360) CALL terminate("WRONG [ASSEMBLY] - HEX ASY AZM")

aInf_Loc => hAsyTypInfo(iAsy)

nDataField = len_trim(dataline)
CALL fndchara(dataline, ipos, nSpt, SLASH)

IF (aInf_Loc%nPin .NE. 0) CALL terminate("WRONG [ASSEMBLY] - OVERLAPPED ASY INPUT")
! ----------------------------------------------------
!               02. CASE : Sng Cel
! ----------------------------------------------------
IF (nSpt .EQ. 1) THEN
  READ (dataline(ipos(1)+1:nDataField),*) iCase
  
  READ (indev,'(a256)') oneline
  IF (Master) CALL message(io8, FALSE, FALSE, oneline) ! ECHO
  READ (oneline, *) ii(1)
  
  IF (iCase .NE. 1) CALL terminate("WRONG [ASSEMBLY] - ASY TYPE INPUT")
  
  aInf_Loc%nPin = 1
  
  CALL dmalloc(aInf_Loc%PinIdx, 1, 1)
  
  aInf_Loc%PinIdx = ii(1)
  aInf_Loc%pF2F   = hCel(ii(1))%pF2F
  aInf_Loc%pPch   = hCel(ii(1))%pF2F * Sq3Inv
  aInf_Loc%aiF2F  = hCel(ii(1))%pF2F
    
  RETURN
END IF
! ----------------------------------------------------
!               03. READ : Asy Data
! ----------------------------------------------------
IF (nSpt .NE. 4) CALL terminate("WRONG [ASSEMBLY] - ASY INPUT")

READ (dataline(ipos(4)+1:nDataField),*) aInf_Loc%aiF2F
READ (dataline(ipos(3)+1:nDataField),*) aInf_Loc%pF2F
READ (dataline(ipos(2)+1:nDataField),*) aInf_Loc%nPin
READ (dataline(ipos(1)+1:nDataField),*) aInf_Loc%gTyp

aInf_Loc%pPCH  = aInf_Loc%pF2F  * Sq3Inv
aInf_Loc%aiPch = aInf_Loc%aiF2F * Sq3Inv
! ----------------------------------------------------
!               04. READ : Pin Data
! ----------------------------------------------------
CALL dmalloc(aInf_Loc%PinIdx, 2*aInf_Loc%nPin-1, 2*aInf_Loc%nPin-1)

jfr = 0
jto = aInf_Loc%nPin

DO i = 1, 2 * aInf_Loc%nPin - 1
  READ (indev, '(a256)') oneline
  IF (Master) CALL message(io8, FALSE, FALSE, oneline) ! ECHO
  IF (probe.EQ.BANG .OR. probe.EQ.POUND) CYCLE
  
  n = nFields(oneline)
  READ (oneline, *) (ii(j), j = 1, n)
  
  IF (n .EQ. 1) THEN
    aInf_Loc%PinIdx = ii(1)
    
    EXIT
  END IF
  
  IF (n .NE. jto) CALL terminate("WRONG [ASSEMBLY] - ASY INPUT MATRIX")
  
  DO j = 1, jto
    k = j + jfr
    
    aInf_Loc%PinIdx(k, i) = ii(j)
  END DO
  
  IF (i < aInf_Loc%nPin) THEN
    jto = jto + 1
  ELSE
    jfr = jfr + 1
    jto = jto - 1
  END IF
END DO

NULLIFY (aInf_Loc)
! ----------------------------------------------------

END SUBROUTINE HexRead_Asy
! ------------------------------------------------------------------------------------------------------------
!                                     07. HEX READ : Core
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexRead_Core(indev, dataline0)

USE allocs
USE geom,    ONLY : Core, Albedo, nAsyType0, nZ, nPinType, nGapPinType
USE ioutil,  ONLY : toupper
USE HexUtil, ONLY : HexChkRange_INT
USE HexData, ONLY : hCore, nAsyCore, nhAsy, Asy2Dto1DMap, Asy1Dto2DMap, hAsyTypInfo, nVA, &
                    RodPin, GapPin, hCel, gCel
USE PE_MOD,  ONLY : PE

IMPLICIT NONE

INTEGER, INTENT(IN) :: indev 
CHARACTER(256), INTENT(IN) :: dataline0
CHARACTER(256) :: dataline, dataline2

INTEGER :: i, j, k, n, iAng, jfr, jto, nAsy, nAsyTot
INTEGER :: iAsy, jAsy, ixPin, iyPin, nPin, iz, jPin
INTEGER :: ii(100)

LOGICAL :: Master
LOGICAL, SAVE :: lfirst = TRUE
! ----------------------------------------------------

IF (.NOT. lfirst) CALL terminate ("WRONG [RAD_CONF] - INPUTTED MORE THAN ONE ")

lfirst   = FALSE
dataline = dataline0
Master   = PE%master

ii   = 0
nVA  = 0
nAsy = 0
n    = nfields(oneline) - 1

IF (n .EQ. 1) THEN
  READ (dataline,*) iAng
  
  hLgc%lAzmRef = TRUE
ELSE
  READ (dataline,*) iAng, dataline2
  
  CALL toupper(dataline2)
  
  hLgc%lAzmRef = dataline2 .NE. 'ROT'
  hLgc%lAzmRot = dataline2 .EQ. 'ROT'
END IF

CALL dmalloc(hCore, Core%nya, Core%nya)
! ----------------------------------------------------
!               CASE : Sng Asy
! ----------------------------------------------------
IF (Core%nya .EQ. 1) THEN
  hLgc%lSngAsy = TRUE
  hLgc%l360    = TRUE
  hLgc%iSym    = 4
  
  CALL dmalloc0(Asy2Dto1Dmap, 0, 2, 0, 2)
  CALL dmalloc (Asy1Dto2Dmap, 2, 1)
    
  Asy2Dto1Dmap(1,1) = 1
  Asy1Dto2Dmap(:,1) = 1
  
  READ (indev,'(a256)') oneline
  IF (Master) CALL message(io8, FALSE, FALSE, oneline) ! ECHO
  READ (oneline, *) ii(1)
  
  hCore(1, 1) = ii(1)
  
  nAsy = 1
  
  IF (hAsyTypInfo(hCore(1, 1))%nPin .EQ. 1) THEN
    hLgc%lSngAsy = FALSE
    hLgc%lSngCel = TRUE
    hLgc%iSym    = 5
  END IF
! ----------------------------------------------------
!               CASE : 360
! ----------------------------------------------------
ELSE IF (iAng .EQ. 360) THEN
  nAsyCore   = (Core%nya + 1)/2
  nAsyTot    = 3 * nAsyCore * (nAsyCore - 1) + 1
  hLgc%l360  = TRUE
  hLgc%iSym  = 3
  
  CALL dmalloc0(Asy2Dto1Dmap, 0, 2 * nAsyCore, 0, 2 * nAsyCore)
  CALL dmalloc (Asy1Dto2Dmap, 2, nAsyTot)
    
  jfr = 0
  jto = nAsyCore

  DO i = 1, Core%nya
    READ (indev,'(a256)') oneline
    IF (Master) CALL message(io8, FALSE, FALSE, oneline) ! ECHO
    IF (probe.EQ.BANG .OR. probe.EQ.POUND) CYCLE
    
    n = nFields(oneline)
    READ (oneline, *) (ii(j), j = 1, n)
    
    DO j = 1, jto
      k = j + jfr
      hCore(k, i) = ii(j)
      
      nAsy = nAsy + 1
      
      Asy2Dto1DMap(k, i)    = nAsy
      Asy1Dto2DMap(1, nAsy) = k
      Asy1Dto2DMap(2, nAsy) = i
      
      IF (ii(j) .NE. 0) CYCLE
      
      nVA  = nVA  + 1
      nAsy = nAsy - 1 ! DON'T COUNT VOID ASY
            
      Asy2Dto1DMap(k, i) = 0
    END DO
    
    IF (i < nAsyCore) THEN
      jto = jto + 1
    ELSE
      jto = jto - 1
      jfr = jfr + 1
    END IF
  END DO
! ----------------------------------------------------
!               CASE : 120
! ----------------------------------------------------
ELSE IF (iAng .EQ. 120) THEN
  nAsyCore   = Core%nya
  nAsyTot    = nAsyCore * nAsyCore
  hLgc%l120  = TRUE
  hLgc%iSym  = 2
  
  CALL dmalloc0(Asy2Dto1Dmap, 0, nAsyCore+1, 0, nAsyCore+1)
  CALL dmalloc (Asy1Dto2Dmap, 2, nAsyTot)
    
  jto = nAsyCore
  
  DO i = 1, Core%nya
    READ (indev,'(a256)') oneline
    IF (Master) CALL message(io8, FALSE, FALSE, oneline) ! ECHO
    IF (probe.EQ.BANG .OR. probe.EQ.POUND) CYCLE
    
    n = nFields(oneline)
    READ (oneline, *) (ii(j), j = 1, n)
    
    DO j = 1, jto
      hCore(j, i) = ii(j)
      
      nAsy = nAsy + 1
      
      Asy2Dto1DMap(j, i)    = nAsy
      Asy1Dto2DMap(1, nAsy) = j
      Asy1Dto2DMap(2, nAsy) = i
      
      IF (ii(j) .NE. 0) CYCLE
      
      nVA  = nVA  + 1
      nAsy = nAsy - 1 ! DON'T COUNT VOID ASY
            
      Asy2Dto1DMap(k, i) = 0
    END DO
  END DO
! ----------------------------------------------------
!               CASE : 060
! ----------------------------------------------------
ELSE IF (iAng .EQ. 60) THEN
  nAsyCore   = Core%nya
  nAsyTot    = nAsyCore * (nAsyCore + 1) / 2
  hLgc%l060  = TRUE
  hLgc%iSym  = 1
  
  CALL dmalloc0(Asy2Dto1Dmap, 0, nAsyCore+1, 0, nAsyCore+1)
  CALL dmalloc (Asy1Dto2Dmap, 2, nAsyTot)
    
  jfr = 0
  jto = nAsyCore
  
  DO i = 1, Core%nya
    READ (indev,'(a256)') oneline
    IF (Master) CALL message(io8, FALSE, FALSE, oneline) ! ECHO
    IF (probe.EQ.BANG .OR. probe.EQ.POUND) CYCLE
    
    n = nFields(oneline)
    READ (oneline, *) (ii(j), j = 1, n)
    
    DO j = 1, jto
      k = j + jfr
      hCore(k, i) = ii(j)
      
      nAsy = nAsy + 1
      
      Asy2Dto1DMap(k, i)    = nAsy
      Asy1Dto2DMap(1, nAsy) = k
      Asy1Dto2DMap(2, nAsy) = i
      
      IF (ii(j) .NE. 0) CYCLE
      
      nVA  = nVA  + 1
      nAsy = nAsy - 1 ! DON'T COUNT VOID ASY
            
      Asy2Dto1DMap(k, i) = 0
    END DO
    
    jfr = jfr + 1
    jto = jto - 1
  END DO
ELSE
  CALL terminate("WRONG [RAD_CONF] - CORE SYMMETRY")
END IF

nhAsy = nAsy
! ----------------------------------------------------
!               03. SET : luse
! ----------------------------------------------------
DO iAsy = 1, nhAsy
  jAsy = hCore(Asy1Dto2DMap(1, iAsy), Asy1Dto2DMap(2, iAsy))
  
  CALL HexChkRange_INT(jAsy, 1, nAsyType0, "WRONG [RAD_CONF] - CORE ASY INPUT")
  
  hAsyTypInfo(jAsy)%luse = TRUE
  
  nPin = hAsyTypInfo(jAsy)%nPin
  
  DO ixPin = 1, 2 * nPin - 1
    DO iyPin = 1, 2 * nPin - 1
      jPin = hAsyTypInfo(jAsy)%PinIdx(ixPin, iyPin)
      
      IF ((jPin < 1).OR.(jPin > nPinType)) CYCLE
      
      RodPin(jPin)%luse = TRUE
      
      DO iz = 1, nZ
        hCel(RodPin(jPin)%iCel(iz))%luse = TRUE
      END DO
    END DO
  END DO
  
  IF (hLgc%lSngCel) CYCLE
  
  jPin = hAsyTypInfo(jAsy)%gTyp
  
  IF ((jPin < 1).OR.(jPin > nGapPinType)) CYCLE
  
  GapPin(jPin)%luse = TRUE
  
  DO iz = 1, nZ
    gCel(GapPin(jPin)%iCel(iz))%luse = TRUE
  END DO
END DO
! ----------------------------------------------------

END SUBROUTINE HexRead_Core
! ------------------------------------------------------------------------------------------------------------
!                                     08. HEX READ : Albedo
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexRead_Albedo(dataline0)

USE PARAM, ONLY : EPSM5

IMPLICIT NONE

CHARACTER(256), INTENT(IN) :: dataline0
CHARACTER(256) :: dataline

INTEGER :: k

REAL :: Albedo(3)
! ----------------------------------------------------

dataline = dataline0

READ (dataline,*) (Albedo(k),k=1,3)

IF (abs(Albedo(1)) < epsm5) THEN
  hLgc%lRadRef = TRUE
ELSE
  hLgc%lRadVac = TRUE
END IF

IF (abs(Albedo(2)) < epsm5) THEN
  hLgc%lAxRef(1) = TRUE
ELSE
  hLgc%lAxVac(1) = TRUE
END IF

IF (abs(Albedo(3)) < epsm5) THEN
  hLgc%lAxRef(2) = TRUE
ELSE
  hLgc%lAxVac(2) = TRUE
END IF
! ----------------------------------------------------

END SUBROUTINE HexRead_Albedo
! ------------------------------------------------------------------------------------------------------------
!                                     09. HEX READ : Vss
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexRead_Vss(dataline0)

USE PARAM,   ONLY : ZERO
USE HexData, ONLY : hVss

IMPLICIT NONE

CHARACTER(256), INTENT(IN) :: dataline0
CHARACTER(256) :: dataline

INTEGER :: nSpt, nDataField
INTEGER :: iVss, Tmp
INTEGER :: ipos(100)
! ----------------------------------------------------

dataline = dataline0

READ (dataline,*) iVss

nDataField = len_trim(dataline)
CALL fndchara(dataline, ipos, nSpt, SLASH)

SELECT CASE (nSpt)
CASE (2)
  READ (dataline(ipos(2)+1:nDataField), *) hVss(iVss)%zSt, hVss(iVss)%zEd
  READ (dataline(ipos(1)+1:nDataField), *) hVss(iVss)%Mat
  READ (dataline, *) Tmp, hVss(iVss)%Rad(1), hVss(iVss)%Rad(2)
  
  hVss(iVss)%Cnt = ZERO
CASE (3)
  READ (dataline(ipos(3)+1:nDataField), *) hVss(iVss)%zSt, hVss(iVss)%zEd
  READ (dataline(ipos(2)+1:nDataField), *) hVss(iVss)%Mat
  READ (dataline(ipos(1)+1:nDataField), *) hVss(iVss)%Rad(1), hVss(iVss)%Rad(2)
  
  READ (dataline, *) Tmp, hVss(iVss)%Cnt(1), hVss(iVss)%Cnt(2)
CASE DEFAULT
  CALL terminate("WRONG [VSS] - VESSEL SLASH")
END SELECT
! ----------------------------------------------------

END SUBROUTINE HexRead_Vss
! ------------------------------------------------------------------------------------------------------------
!                                     10. HEX READ : Vyg
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexRead_Vyg(dataline0)

USE HexData, ONLY : vAsyTyp, vRefTyp, vMat1, vMat2, vFxr, vzSt, vzEd

IMPLICIT NONE

CHARACTER(256), INTENT(IN) :: dataline0
CHARACTER(256) :: dataline

INTEGER :: nSpt, nDataField, Tmp
INTEGER :: ipos(100)
! ----------------------------------------------------

dataline  = dataline0

nDataField = len_trim(dataline)
CALL fndchara(dataline, ipos, nSpt, SLASH)

IF (nSpt .NE. 2) CALL terminate("WRONG [VYG] - WRONG VYGORODKA INPUT")

READ (dataline(ipos(2)+1:nDataField), *) vzSt, vzEd
READ (dataline(ipos(1)+1:nDataField), *) vMat1, vMat2, vFXR
READ (dataline, *) vAsyTyp, vRefTyp
! ----------------------------------------------------

END SUBROUTINE HexRead_Vyg
! ------------------------------------------------------------------------------------------------------------
!                                     11. HEX READ : Opt
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexRead_Opt(dataline0)

USE HexData,     ONLY : hLgc, nInnMOCItr
USE CNTL,        ONLY : nTracerCntl
USE ItrCNTL_mod, ONLY : ItrCntl

#ifdef __INTEL_MKL
USE MKL_3D,  ONLY : mklcntl
#endif

IMPLICIT NONE

CHARACTER(256), INTENT(IN) :: dataline0
CHARACTER(256) :: dataline

INTEGER :: nInnCMFDmax = 100
INTEGER :: nOutCMFDmax = 100

REAL :: InnCMFDConvCrit = 0.01
REAL :: OutCMFDConvCrit = 0.1

LOGICAL :: lCmfd = TRUE
! ----------------------------------------------------

dataline = dataline0

READ (dataline,*) lCmfd, hLgc%lspCMFD, nInnMOCItr, nInnCMFDmax, InnCMFDConvCrit, nOutCMFDmax, OutCMFDConvCrit

!mklcntl%maxOuter  = nOutCmfdmax
!mklcntl%maxInner  = nInnCmfdmax
!mklcntl%outerConv = OutCMFDConvCrit
!mklcntl%innerConv = InnCMFDConvCrit

ItrCntl%InSolverItrCntl%ninmax   = nInnCMFDmax
ItrCntl%InSolverItrCntl%convcrit = InnCMFDConvCrit

nTracerCntl%lCMFD = lCMFD
! ----------------------------------------------------

END SUBROUTINE HexRead_Opt
! ------------------------------------------------------------------------------------------------------------
!                                     13. HEX READ : Corner Stiffener
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexRead_CrnStff(dataline0)

USE PARAM,   ONLY : TRUE
USE ioutil,  ONLY : nfields
USE HexData, ONLY : csMat, csWdt, csLgh, hLgc

IMPLICIT NONE

CHARACTER(256), INTENT(IN) :: dataline0
CHARACTER(512) :: dataline
! ----------------------------------------------------

dataline = dataline0

IF (nfields(dataline) .NE. 3) CALL terminate("WRONG [CRN_STFF]")

READ (dataline, *) csMat, csWdt, csLgh

hLgc%lcrnstff = TRUE
! ----------------------------------------------------

END SUBROUTINE HexRead_CrnStff
! ------------------------------------------------------------------------------------------------------------

END MODULE HexInpCard