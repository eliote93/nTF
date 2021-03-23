SUBROUTINE SetThCell(CellInfo, THOpt, nTracerCntl)
USE PARAM
USE TYPEDEF,              ONLY : Cell_Type,        THCell_Type
USE Cntl,                 ONLY : nTracerCntl_Type
USE TH_Mod,               ONLY : THOpt_Type
USE Material_Mod,         ONLY : Mixture
USE ALLOCS
IMPLICIT NONE
TYPE(Cell_Type) :: CellInfo
TYPE(ThOpt_Type) :: ThOpt
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(THCell_Type), POINTER :: THCell
REAL :: FsrR(2), ThR(0:100), R(0:100)
REAL :: rpallet, delR, AreaFxr, AreaPallet, frac
INTEGER :: i, j, k, ifsr, imix
INTEGER :: nCondRing, nrpellet, nFxr
INTEGER :: Reg(2)


ALLOCATE(CellInfo%THCell)
THCell => CellInfo%THCell
IF(CellInfo%lRect) THEN
  THCell%FuelReg = (/1, -1/); THCell%CldReg = (/1, -1/)
  THCell%CoolReg = (/1, -1/)
  IF(CellInfo%lFuel)  THCell%FuelReg = (/1,1/)
  IF(.NOT. CellInfo%lFuel) THCell%CoolReg = (/1,1/)
  RETURN
ENDIF

nFxr = CellInfo%nFxr
THCell%FuelReg = (/1000, 0/); THCell%CoolReg = (/1000, 0/)
THCell%CldReg = (/1000, 0/)
!Region Information
DO i = 1, nFxr
  ifsr = CellInfo%MapFxr2FsrIdx(1, i); imix = CellInfo%iReg(ifsr)
  IF(Mixture(imix)%lh2o) THEN
    ThCell%CoolReg(1) = min(ThCell%CoolReg(1), i)
    ThCell%CoolReg(2) = max(ThCell%CoolReg(2), i)
  ENDIF
  IF(Mixture(imix)%lfuel) THEN
    ThCell%FuelReg(1) = min(ThCell%FuelReg(1), i)
    ThCell%FuelReg(2) = max(ThCell%FuelReg(2), i)
  ENDIF
  IF(Mixture(imix)%lCLD) THEN
    ThCell%CldReg(1) = min(ThCell%CldReg(1), i)
    ThCell%CldReg(2) = max(ThCell%CldReg(2), i)
  ENDIF  
ENDDO

IF(.NOT. CellInfo%lfuel) THEN
  ThCell%FuelReg = (/nFxr+1, nFxr/)
END IF

rpallet = ThOpt%PinDim(1)*0.1; nrpellet = ThOpt%nrpellet
delR = rPallet / nrpellet
ThR = 0
IF(.NOT. nTracerCntl%lSimpleTh) THEN
  DO i = 1, nrpellet
    ThR(i) = delR*I
  ENDDO
ELSE
  DO i = 1, nrpellet
    ThR(i) = rPallet * SQRT(1._8/nrpellet*i)
  ENDDO
ENDIF
! Neutronic <-> TH node Mapping
CALL DMALLOC(ThCell%FuelMapping, 2, nFxr)
CALL DMALLOC(ThCell%Frac, nrpellet, nFxr)
CALL DMALLOC(ThCell%CondMapping, 2, nrpellet)
CALL DMALLOC(ThCell%InvFrac, nFXR, nrpellet)
FsrR = 0
AreaPallet = rpallet * rpallet
DO i = ThCell%FuelReg(2), ThCell%FuelReg(1), -1 
  FsrR(2) = CellInfo%Geom%circle(3, i-1)
  IF(i .NE. nFxr) FsrR(1) = CellInfo%Geom%circle(3, i)
  Reg = (/1000, -1/)
  DO j = 1, nrpellet
    IF(ThR(j) .GT. FsrR(1)) Reg(1) = MIN(Reg(1), j)
    IF(ThR(j) .LT. FsrR(2)) Reg(2) = Max(Reg(2), j)
  ENDDO
  R(Reg(1):Reg(2)) = ThR(Reg(1):Reg(2))
  R(Reg(1)-1) = FsrR(1);
  IF(Reg(2) .NE. nrpellet .AND. ABS(ThR(Reg(2))-FsrR(2)) .GT. 1.E-4) THEN
    Reg(2)=Reg(2) + 1
    R(Reg(2)) = FsrR(2)
  ENDIF
  nCondRing = Reg(2) - Reg(1) + 1
  AreaFxr = FsrR(2) * FsrR(2) - FsrR(1) * FsrR(1)
  
  DO j = Reg(1), Reg(2)
    frac = (R(j) * R(j) - R(j-1) * R(j-1)) / AreaFxr
    ThCell%Frac(j, i) = frac
    ThCell%InvFrac(i,j) = frac*AreaFxr/AreaPallet
  ENDDO
  ThCell%FuelMapping(:, i) = Reg
ENDDO


DO i = 1, nrpellet
  Reg = (/1000, -1/)
  DO j = ThCell%FuelReg(2), ThCell%FuelReg(1), -1
    IF(ThCell%Frac(i, j) .EQ. 0._8) CYCLE
    Reg(1) = MIN(Reg(1), j); Reg(2) = Max(Reg(2), j)
  ENDDO
  ThCell%CondMapping(:, i) = Reg
ENDDO

NULLIFY(THCell)
END SUBROUTINE