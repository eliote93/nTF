MODULE HexNgh

IMPLICIT NONE

CONTAINS
! ------------------------------------------------------------------------------------------------------------
!                                     01. HEX SET : Asy Ngh
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetAsyNgh(iAsy)

USE HexType, ONLY : Type_HexAsy
USE HexData, ONLY : Asy2Dto1DMap, hAsy, hLgc

IMPLICIT NONE

INTEGER :: iAsy, iBndy, jaX, jaY

INTEGER, PARAMETER :: xNxt(6) = [ 0, 1, 1, 0, -1, -1] ! 1 = NE, 2 = EE, ... (Clock-Wise)
INTEGER, PARAMETER :: yNxt(6) = [-1, 0, 1, 1,  0, -1]

TYPE(Type_HexAsy), POINTER :: hAsy_Loc
! ----------------------------------------------------

hAsy_Loc => hAsy(iAsy)

DO iBndy = 1, 6
  jaX = hAsy_Loc%iaX + xNxt(iBndy)
  jaY = hAsy_Loc%iaY + yNxt(iBndy)
  
  hAsy_Loc%NghIdx(iBndy) = Asy2Dto1DMap(jaX, jaY) ! 0 : VOID
END DO

NULLIFY(hAsy_Loc)
! ----------------------------------------------------

END SUBROUTINE HexSetAsyNgh
! ------------------------------------------------------------------------------------------------------------
!                                     02. HEX SET : Asy Rot Ngh
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetAsyRotNgh(iAsy)

USE PARAM,   ONLY : HALF
USE HexData, ONLY : hAsy, Asy2Dto1Dmap, Asy1Dto2Dmap, hLgc

IMPLICIT NONE

INTEGER :: iAsy, ix, iy
REAL    :: cx
! ----------------------------------------------------

IF (.NOT. hLgc%lAzmRot) RETURN

ix = Asy1Dto2Dmap(1, iAsy)
iy = Asy1Dto2Dmap(2, iAsy)

IF ((iy .NE. 1) .AND. (ix .NE. iy)) RETURN

cx = real(ix + 1) * HALF
iy = iy + int(2. * (cx - real(iy)))

hAsy(iAsy)%RotNgh = Asy2Dto1Dmap(ix, iy)
! ----------------------------------------------------

END SUBROUTINE HexSetAsyRotNgh
! ------------------------------------------------------------------------------------------------------------
!                                     03. HEX SET : Asy Typ Pin Ngh
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetAsyTypPinNgh(iaTyp)

USE HexData, ONLY : hLgc

IMPLICIT NONE

INTEGER :: iaTyp
! ----------------------------------------------------

IF (hLgc%lSngCel) THEN
  RETURN
ELSE IF (hLgc%lspCMFD) THEN
  CALL HexSetAsyTypSPNgh(iaTyp)
ELSE
  CALL HexSetAsyTypMPNgh(iaTyp)
END IF
! ----------------------------------------------------

END SUBROUTINE HexSetAsyTypPinNgh
! ------------------------------------------------------------------------------------------------------------

END MODULE HexNgh