SUBROUTINE HexSetAsyRay_SngCel()

USE PARAM,   ONLY : TRUE, MESG
USE PE_MOD,  ONLY : PE
USE IOUTIL,  ONLY : message
USE FILES,   ONLY : io8
USE HexType, ONLY : Type_HexRayPinInfo
USE HexData, ONLY : NumMray, hmRay, RayCel, RayPinInfo, haRay, AsyVtx, AsyEqn
USE HexUtil, ONLY : ChkEqnTwoVtx, SolveLineEqn

USE HexAsyRayConst
USE HexAsyRayIntSct

IMPLICIT NONE

INTEGER :: imRay, iBndy, nPt
LOGICAL :: lSol
REAL :: RayEqn(3), Sol(2), PinPt(2, 2)

TYPE(Type_HexRayPinInfo), POINTER :: rInf_Loc
! ----------------------------------------------------

WRITE(MESG, '(9X, I2, X, A3)') 1, '...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, MESG)
! ----------------------------------------------------
!               01. INIT
! ----------------------------------------------------
ALLOCATE (RayPinInfo (1)); rInf_Loc => RayPinInfo(1)

rInf_Loc%nBndy  = 6
rInf_Loc%VtxTyp = 3

rInf_Loc%Vtx(1:2, 1:7, 1) = AsyVtx(1:2, 1:7)
rInf_Loc%Eqn(1:3, 1:6, 1) = AsyEqn(1:3, 1:6)

ALLOCATE (RayCel (1))

CALL HexSetRayRodCel_12(RayCel, 1)
! ----------------------------------------------------
!               02. CALC : Msh Int Sct
! ----------------------------------------------------
!$OMP PARALLEL PRIVATE(imRay, RayEqn, nPt, iBndy, lSol, PinPt, Sol)
!$OMP DO SCHEDULE(GUIDED)
DO imRay = 1, NumMray(0)
  haRay(1, 1, imRay)%nhpRay = 1
  
  ALLOCATE (haRay(1, 1, imRay)%CelRay (1))
  
  haRay(1, 1, imRay)%CelRay(1)%hPinIdx = 1
  haRay(1, 1, imRay)%CelRay(1)%hSufIdx = 1 ! TRASH
  
  RayEqn = hmRay(imRay)%Eq ! Origin : Cel Cnt
  nPt    = 0
  
  DO iBndy = 1, rInf_Loc%nBndy(1)
    lSol = ChkEqnTwoVtx(RayEqn, rInf_Loc%Vtx(1:2, iBndy, 1), rInf_Loc%Vtx(1:2, iBndy + 1, 1))
    
    IF (.NOT. lSol) CYCLE
    
    nPt = nPt + 1
    
    CALL SolveLineEqn(RayEqn, rInf_Loc%Eqn(1:3, iBndy, 1), Sol, lSol)
    
    PinPt(1:2, nPt) = Sol(1:2)
  END DO
  
  CALL HexSetRayIntSct_Msh(haRay(1, 1, imRay)%CelRay(1), RayEqn, 1, 1, PinPt)
END DO
!$OMP END DO
!$OMP END PARALLEL

NULLIFY (rInf_Loc)
! ----------------------------------------------------

END SUBROUTINE HexSetAsyRay_SngCel