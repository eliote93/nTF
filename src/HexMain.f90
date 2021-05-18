#include <defines.h>
SUBROUTINE HexMain

USE PARAM,  ONLY : TRUE, FALSE, MESG
USE geom,   ONLY : nbd, ncbd, nz
USE CNTL,   ONLY : nTracerCntl
USE PE_MOD, ONLY : PE
USE IOUTIL, ONLY : message, terminate
USE FILES,  ONLY : io8

USE HexRayBasic
USE HexRayConst
USE HexCmfdConst
USE HexCnP
USE HexTst

USE HexData,     ONLY : hLgc, ncBss, ngBss, NumMray, haRay, nGeoTyp
USE HexGeoConst, ONLY : HexSetVyg
USE HexPinConst, ONLY : HexSetVss

IMPLICIT NONE

INTEGER :: icBss
! ----------------------------------------------------

nbd = 6

IF (ngBss .EQ. 1) ncbd = 6
IF (ngBss .NE. 1) ncbd = 10

nTracerCntl%MultigridLV = 6

IF (nz .EQ. 1) nTracerCntl%l3dim = FALSE

CALL ConvertXs

WRITE(MESG, '(A)') '-------------------------------------------------------------------'
IF(PE%Master) CALL message(io8, TRUE, TRUE, MESG)
! ----------------------------------------------------
! Geo
WRITE(MESG, '(A)') 'HEX : Set Geom ...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, MESG)

CALL HexSetGeo

! CMFD map
WRITE(MESG, '(A)') '      Set CMFD Map ...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, MESG)

CALL HexSetHcPin

! Vyg & Vss
CALL HexSetVyg
CALL HexSetVss ! Vss must follow Vyg

! CnP
CALL HexCnPnT
! ----------------------------------------------------
! Ray Basic Data
WRITE(MESG, '(A)') '      Set Ray Basic Data ...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, MESG)

CALL HexSetRayParam
CALL HexSetModRay
CALL HexSetModRayNxt

! Asy Ray Base
WRITE(MESG, '(A)') '      Set Asy Ray Base ...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, MESG)

ALLOCATE (haRay (nGeoTyp, ncBss, NumMray(0)))

DO icBss = 1, ncBss
  CALL HexSetAsyRay(icBss)
END DO

! Ray
WRITE(MESG, '(A)') '      Set Core Ray ...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, MESG)

CALL HexSetCoreRay

WRITE(MESG, '(A)') '      Connect Core Ray ...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, MESG)

CALL HexSetRotRay

WRITE(MESG, '(A)') '      Convert Hex Ray ...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, MESG)

CALL ConvertRay

WRITE(MESG, '(A)') '-------------------------------------------------------------------'
IF(PE%Master) CALL message(io8, TRUE, TRUE, MESG)
! ----------------------------------------------------
! Tst
!CALL HexTstHcPin
!CALL HexTsthPinInfo
!CALL HexPrintPinTyp
!CALL HexTsthmRay
!CALL HexTsthaRay(1, 1, 53)
!CALL HexTstAsyRaySegNum
!CALL HexTsthcRay
!CALL HexTsthRotRay
! ----------------------------------------------------

END SUBROUTINE HexMain