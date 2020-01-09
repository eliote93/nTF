MODULE TYPEDEF_COMMON
    
IMPLICIT NONE

TYPE superPin_Type
  INTEGER :: nx, ny, nxy
  INTEGER :: ix, iy, iFuelPin
  INTEGER :: NeighIdx(15), NeighSurfIdx(15)
  INTEGER, POINTER :: pin(:), pin2D(:, :)
  REAL :: BdLength(6), Center2SurfaceL(6), Area
  LOGICAL, POINTER :: lFuel(:)
  
  ! HEX
  INTEGER :: nNgh = 4, nBdmPin(6)
  INTEGER :: BdMPidx(2, 6), BdMPsuf(2, 6)
  INTEGER :: NghBd(15)
  REAL    :: NghLgh(15)
END TYPE

TYPE CoreXsMac_Type
  REAL, POINTER :: XSt(:, :)
  REAL, POINTER :: XStr(:, :)
  REAL, POINTER :: XSa(:, :)
  REAL, POINTER :: XSnf(:, :)
  REAL, POINTER :: XSkf(:, :)
  REAL, POINTER :: XSsm(:, :), XSsmP1(:, :), XSsmP2(:, :), XSsmP3(:, :)
END TYPE

END MODULE