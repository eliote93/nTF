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

USE HexData, ONLY : hLgc, ncBss, NumMray, haRay, nGeoTyp
USE HexCP,   ONLY : ConvertRay, ConvertXs
USE HexTst,  ONLY : HexTstHcPin, HexTstAsyRaySegNum, HexTsthmRay, HexTsthaRay, HexTsthcRay, HexTsthRotRay

IMPLICIT NONE

INTEGER :: icBss
! ----------------------------------------------------

nbd  = 6
ncbd = 15

nTracerCntl%MultigridLV = 6

IF (nz .EQ. 1) nTracerCntl%l3dim = .FALSE.

CALL ConvertXs

WRITE(MESG, '(A)') '------------------------------------------------------------------------'
IF(PE%Master) CALL message(io8, TRUE, TRUE, MESG)
! ----------------------------------------------------
!               01. SET : Geo
! ----------------------------------------------------
WRITE(MESG, '(A)') 'HEX : Set Geom ...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, MESG)

CALL HexSetGeo
! ----------------------------------------------------
!               02. SET : CMFD map
! ----------------------------------------------------
WRITE(MESG, '(A)') '      Set CMFD Map ...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, MESG)

CALL HexSetHcPin
!CALL HexTstHcPin
! ----------------------------------------------------
!               03. SET : Ray Basic Data
! ----------------------------------------------------
WRITE(MESG, '(A)') '      Set Ray Basic Data ...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, MESG)

CALL HexSetRayParam
CALL HexSetModRay
CALL HexSetModRayNxt
!CALL HexTsthmRay
! ----------------------------------------------------
!               04. SET : Asy Ray Base
! ----------------------------------------------------
WRITE(MESG, '(A)') '      Set Asy Ray Base ...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, MESG)

ALLOCATE (haRay (nGeoTyp, ncBss, NumMray(0)))

DO icBss = 1, ncBss
  CALL HexSetAsyRay(icBss)
END DO

!CALL HexTsthaRay(1, 1, 53)
!CALL HexTstAsyRaySegNum
! ----------------------------------------------------
!               04. SET : Ray
! ----------------------------------------------------
WRITE(MESG, '(A)') '      Set Core Ray ...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, MESG)

CALL HexSetCoreRay
!CALL HexTsthcRay

WRITE(MESG, '(A)') '      Connect Core Ray ...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, MESG)

CALL HexSetRotRay
!CALL HexTsthRotRay

WRITE(MESG, '(A)') '      Convert Hex Ray ...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, MESG)

CALL ConvertRay

WRITE(MESG, '(A)') '------------------------------------------------------------------------'
IF(PE%Master) CALL message(io8, TRUE, TRUE, MESG)
! ----------------------------------------------------

END SUBROUTINE HexMain