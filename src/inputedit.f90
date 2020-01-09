SUBROUTINE inputEdit(io)
USE PARAM
USE GEOM,        ONLY : CORE, AsyGeom
USE Rays,        ONLY : RayInfo
USE CORE_MOD,    ONLY : eigv
USE IOUTIL,      ONLY : message
USE CNTL,        ONLY : nTracerCntl
IMPLICIT NONE
INTEGER :: io
!USE FILES,       ONLY : io8
!USE TIMER,       ONLY : TimeChk
!USE itrcntl_mod, ONLY : ItrCntl
CALL message(io, TRUE, TRUE, 'Generating Input Edits...')
WRITE(io, '(a)') hbar2
!Cell Info Out
CALL CellInfoOut(io, Core%CellInfo, Core%nCellType)
!Pin Info Out
CALL PinInfoOut(io, Core%PinInfo, Core%nPinType, Core%nz)
!
CALL AssemblyOut(io, Core%AsyInfo, AsyGeom, Core%nasytype)
!
CALL LoadingPattern(io, Core)
!
CALL RayInfoOut(io, RayInfo)

IF(nTracerCntl%lXsLib) CALL PrintLibEStructure(io)
WRITE(io, '(a)') hbar2
END SUBROUTINE

SUBROUTINE AsyTypeOut(AsyInfo, nas, io)
USE PARAM
USE TYPEDEF, ONLY : AsyInfo_Type
TYPE(AsyInfo_type) :: AsyInfo(nas)
INTEGER, INTENT(IN) :: nas
INTEGER, INTENT(in) :: io

INTEGER :: ias, ix, iy
END SUBROUTINE

SUBROUTINE CellInfoOut(IO, CellInfo, ncelltype)
USE PARAM
USE TYPEDEF, ONLY : cell_type
IMPLICIT NONE

TYPE(Cell_Type) :: CellInfo(nCellType)
CHARACTER(12) :: PartialCell
CHARACTER(7)  :: USED
INTEGER :: io, nCellType, nDivAzi

INTEGER :: i, j
101 FORMAT (5x, A4, 2x, A5, 2x, A5, 2x, A12, 2x, A7, 3x, A35)
102 FORMAT (    I9,    I7,    I7, 2x, A12, 2x, A7, 3x, 100(F8.4, '(',I4 ,')' ))
103 FORMAT (5x, A4, 2x, A5, 2x, A5, 2x, A12, 2x, A7, 3x, A25, 2x, A8, 2x, A8)
104 FORMAT (    I9,     I7,     I7, 2x, A12, 2x, A7, 3x, I7)
105 FORMAT (5x, A4, 2x, A5, 2x, A5, 2x, A7, 3x, A25, 2x, A8, 2x, A8)
106 FORMAT (5x, I4, 2x, I5, 2x, I5, 2x, A7, 3x, I17, 10x, F8.4, 2x, F8.4)

WRITE(io, *)
WRITE(io, '(2x, A13)') '= Cell Type ='
WRITE(io, '(3x, A13, 2x, F8.5)') 'Cell Pitch = ', CellInfo(1)%Geom%lx
WRITE(io, '(3x, A16)') 'Cylindrical Cell'
WRITE(io, 101) 'TYPE', '# FSR', '# FXR', 'Partial Cell' , 'USE','RADIUS(Composition or Mixture)'
DO i = 1, nCellType
  IF(CellInfo(i)%lRECT) CYCLE
  WRITE(USED, '(A7)') 'TRUE'
  IF(.NOT. CellInfo(i)%luse) WRITE(USED, '(A7)') 'FALSE'
  WRITE(PartialCell, '(A12)') 'FALSE'
  IF(CellInfo(i)%lCentX) WRITE(PartialCell, '(A12)') 'TRUE (1)'
  IF(CellInfo(i)%lCentY) WRITE(PartialCell, '(A12)') 'TRUE (2)'
  IF(CellInfo(i)%lCentXY) WRITE(PartialCell, '(A12)') 'TRUE (3)'
  nDivAzi = CellInfo(i)%nDivAzi
  WRITE(io ,102) i, CellInfo(i)%nFsr, CellInfo(i)%nFxr, PartialCell, USED,                &
                (CellInfo(i)%Geom%circle(3, j), CellInfo(i)%iReg(nDivAzi*(j-1) + 1), j = CellInfo(i)%Geom%ncircle, 1, -1)
ENDDO

WRITE(io, *)
WRITE(io, '(3x, A16)') 'Rectangular Cell'
WRITE(io, 103) 'TYPE', '# FSR', '# FXR', 'Partial Cell', 'USE', 'Composition or Mixture'
DO i = 1, nCellType
  IF(.NOT. CellInfo(i)%lRECT) CYCLE
  IF(CellInfo(i)%lGap) CYCLE
  WRITE(USED, '(A7)') 'TRUE'
  IF(.NOT. CellInfo(i)%luse) WRITE(USED, '(A7)') 'FALSE'
  WRITE(PartialCell, '(A12)') 'FALSE'
  IF(.NOT. CellInfo(i)%luse) WRITE(USED, '(A7)') 'FALSE'
  IF(CellInfo(i)%lCentX) WRITE(PartialCell, '(A12)') 'TRUE (1)'
  IF(CellInfo(i)%lCentY) WRITE(PartialCell, '(A12)') 'TRUE (2)'
  IF(CellInfo(i)%lCentXY) WRITE(PartialCell, '(A12)') 'TRUE (3)'
  WRITE(io ,104) i, CellInfo(i)%nFsr, CellInfo(i)%nFxr, PartialCell, USED, CellInfo(i)%iReg(1)
ENDDO

WRITE(io, *)
WRITE(io, '(3x, A8)') 'Gap Cell'
WRITE(io, 105) 'TYPE', '# FSR', '# FXR', 'USE', 'Composition or Mixture', 'Pitch x', 'Pitch y'
DO i = 1, nCellType
  IF(.NOT. CellInfo(i)%lGap) CYCLE
  WRITE(USED, '(A7)') 'TRUE'
  IF(.NOT. CellInfo(i)%luse) WRITE(USED, '(A7)') 'FALSE'
  WRITE(PartialCell, '(A12)') 'FALSE'
  IF(CellInfo(i)%lCentX) WRITE(PartialCell, '(A12)') 'TRUE (1)'
  IF(CellInfo(i)%lCentY) WRITE(PartialCell, '(A12)') 'TRUE (2)'
  IF(CellInfo(i)%lCentXY) WRITE(PartialCell, '(A12)') 'TRUE (3)'
  WRITE(io ,106) i, CellInfo(i)%nFsr, CellInfo(i)%nFxr, USED, CellInfo(i)%iReg(1), CellInfo(i)%Geom%lx, CellInfo(i)%Geom%ly
ENDDO
WRITE(io, *)
END SUBROUTINE

SUBROUTINE PinInfoOut(io, PinInfo, nPinType, nz)
USE PARAM
USE TYPEDEF, ONLY : pininfo_type
IMPLICIT NONE
TYPE(PinInfo_Type), INTENT(in) :: PinInfo(nPinType)
INTEGER, INTENT(in) :: io, nPinType, nz

CHARACTER(12) :: PartialPin
CHARACTER(7)  :: USED
INTEGER :: i, j
201 FORMAT (5x, A4, 2x, A12, 2x, A7, 2x, A12)
202 FORMAT (    I9, 2x, A12, 2x, A7, 2x, 100I5)
WRITE(io, *)
WRITE(io, '(2x, A12)') '= Pin Type ='
!WRITE(io, '(3x, A18, 2x, 200F8.5)') 'Axial Mesh Size : ',
WRITE(io, 201) 'TYPE', 'Partial Pin', 'USE', 'Cell Type'
DO i = 1, nPinType
  WRITE(USED, '(A7)') 'TRUE'
  IF(.NOT. PinInfo(i)%luse) WRITE(USED, '(A7)') 'FALSE'
  WRITE(PartialPin, '(A12)') 'FALSE'
  IF(PinInfo(i)%lCentX) WRITE(PartialPin, '(A12)') 'TRUE (1)'
  IF(PinInfo(i)%lCentY) WRITE(PartialPin, '(A12)') 'TRUE (2)'
  IF(PinInfo(i)%lCentXY) WRITE(PartialPin, '(A12)') 'TRUE (3)'
  WRITE(io, 202) i, PartialPin, USED, (pinInfo(i)%Cell(j), j = 1, nz)
ENDDO
WRITE(io, *)
END SUBROUTINE

SUBROUTINE AssemblyOut(io, AsyInfo, AsyGeom, nasytype)
USE PARAM
USE TYPEDEF, ONLY : AsyInfo_type, BasicGeom
USE IOUTIL, ONLY : PrintInt1DarrayTo2Darray
IMPLICIT NONE
TYPE(AsyInfo_Type) :: AsyInfo(nasytype)
TYPE(BasicGeom) :: AsyGeom(0:3)
CHARACTER(120) :: FormatData
INTEGER, INTENT(in) :: io, nasytype
INTEGER :: i, j, k
INTEGER :: ix, iy, ixy, itype

WRITE(io, *)
WRITE(io, '(2x, A18)') '= Assembly Type ='

WRITE(formatdata,'(a)') '(7x, 200I5)'
DO i = 1, nAsyType
  IF(.NOT. AsyInfo(i)%luse) CYCLE
  WRITE(io, '(5x, A7, I5)') 'TYPE = ', i
  itype = 0
  IF(AsyInfo(i)%lCentX) itype = 1
  IF(AsyInfo(i)%lCentY) itype = 2
  IF(AsyInfo(i)%lCentXY) itype = 3
  WRITE(io, '(7x, A32, F8.5, A3, 3x, F8.5)') 'Assembly Pitch(x-dir / y-dir) = ', AsyGeom(itype)%lx, '/', AsyGeom(itype)%ly
  WRITE(io, '(7x, A19, I5, A3, 3x, I5)') 'No. of Pins/Side = ', AsyInfo(i)%nx, '/', AsyInfo(i)%ny
  WRITE(io, *)
  CALL PRINTInt1DArrayTo2DArray(io, AsyInfo(i)%pin, AsyInfo(i)%nxy, AsyInfo(i)%nx, AsyInfo(i)%ny, formatdata)
  IF(i .NE. nAsyType) WRITE(io, *)
ENDDO
WRITE(io, *)

END SUBROUTINE

SUBROUTINE LoadingPattern(io, Core)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type, Asy_Type, Pin_Type
USE IOUTIL, ONLY : PrintInt1DarrayTo2Darray
IMPLICIT NONE
INTEGER :: io
TYPE(CoreInfo_Type) :: Core
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(Pin_Type), POINTER :: Pin(:)
INTEGER :: nxa, nya, nxya
INTEGER :: i, j, k, ix, iy, ixy
INTEGER :: PrintArray(1000)
CHARACTER(7) :: NeighCell(4, 200)
CHARACTER(120) :: FormatData
Pin => Core%Pin
Asy => Core%Asy
nxa = Core%nxa; nya = Core%nya; nxya = Core%nxya

WRITE(io, *)
WRITE(io, '(2x, A24)') '= Core Loading Pattern ='
WRITE(io, *)
WRITE(io, '(7x, A15, 2I6)') 'N. Column/Row =', nxa, nya
WRITE(io, *)
!Core Loading Pattern
WRITE(formatdata,'(a)') '(7x, 200I5)'
!CALL PrintInt1DarrayTo2Darray(io, Core%CoreMap, nxya, nxa, nya, formatdata)
DO iy = 1, nya
  k = 0
  DO ix = 1, nxa
    IF(Core%CoreIdx(ix, iy) .EQ. 0) CYCLE
    k = k + 1; i = Core%CoreIdx(ix, iy)
    PrintArray(k) = Core%CoreMap(i)
  ENDDO
  WRITE(io, '(7x, 200I5)') (PrintArray(i), i = 1, k)
ENDDO
WRITE(io, *)
!Cell Numbering Index Print Out
WRITE(io, '(7x, A)') 'Cell Numbering'
DO i = 1, nxya
  WRITE(io, *)
  WRITE(io, '(7x, A17, i7)') 'Assembly Number = ' , i
  j = Core%CoreMap(i)
  CALL PrintInt1DarrayTo2Darray(io, Asy(i)%GlobalPinIdx, Core%AsyInfo(j)%nxy, Core%AsyInfo(j)%nx, Core%AsyInfo(j)%ny, formatdata)
ENDDO

WRITE(io, *)
WRITE(io, '(7x, A)') 'Neighbor Cells'
DO i = 1, nxya
  WRITE(io, '(7x, A17, i7)') 'Assembly Number = ' , i
  j = Core%CoreMap(i)
  DO iy = 1, Core%AsyInfo(j)%ny
    DO ix = 1, Core%AsyInfo(j)%nx
      ixy = Asy(i)%GlobalPinIdx(Core%AsyInfo(j)%pin2DIdx(ix, iy))
      DO k = 1, 4
        IF(Pin(ixy)%NeighIdx(k) .GT. 0) THEN
          WRITE(NeighCell(k, ix), '(I7)') Pin(ixy)%NeighIdx(k)
        ELSE
          WRITE(NeighCell(k, ix), '(A7)') '-'
        ENDIF
      ENDDO
    ENDDO

    DO k = 1, 4
      WRITE(io, '(7x, 200(A7, 2x))') (NeighCell(k, ix), ix = 1, Core%AsyInfo(j)%nx)
    ENDDO
    WRITE(io, *)
  ENDDO
  IF(Core%AsyInfo(j)%ny .EQ. 0) WRITE(io, *)
ENDDO
NULLIFY(Pin)
NULLIFY(Asy)
END SUBROUTINE

SUBROUTINE RayInfoOut(io, RayInfo)
USE PARAM
USE TYPEDEF, ONLY : RayInfo_Type, AziAngleInfo_Type, PolarAngle_Type

IMPLICIT NONE
TYPE(RayInfo_Type) :: RayInfo
INTEGER :: io

TYPE(AziAngleInfo_Type), POINTER :: AziAngle(:)
TYPE(PolarAngle_Type), POINTER :: PolarAngle(:)
INTEGER :: nAziAng, nPolarAng
INTEGER :: i, j, k

nAziAng = RayInfo%nAziAngle
AziAngle => RayInfo%AziAngle
nPolarAng = RayInfo%nPolarAngle
PolarAngle => RayInfo%PolarAngle

WRITE(io, *)
WRITE(io, '(2x, A)') '=Tracking Ray Information='
WRITE(io, *)
WRITE(io, '(5x, A)') 'Azimuthal Angles'
WRITE(io, '(7x, A)') '   ID/ rays/ start/   Angle/         Sin/         Cos/         Tan/        dels/          wt'
301 format (7x, I4, '/', I6, '/', I6, '/', F8.2, 5('/', e12.4))
302 format (7x, I4, 3('/', e12.4))
DO i = 1, nAziAng
  WRITE(io, 301) i, AziAngle(i)%nModRay, AziAngle(i)%ModRayIdxSt, AziAngle(i)%ang*180._8/pi,                 &
              AziAngle(i)%sinv, AziAngle(i)%cosv, AziAngle(i)%tanv, AziAngle(i)%del, AziAngle(i)%weight
ENDDO
WRITE(io, *)
WRITE(io, '(5x, A)') 'Polar Angles'
WRITE(io, '(5x, A)') '    ID/         Sin/         Cos/          wt'

DO i = 1, nPolarAng
  WRITE(io, 302) i, PolarAngle(i)%sinv, PolarAngle(i)%cosv, PolarAngle(i)%weight
ENDDO
WRITE(io, *)

NULLIFY(AziAngle)
NULLIFY(PolarAngle)
END SUBROUTINE

SUBROUTINE PrintLibEStructure(io)
USE XSLIB_MOD,  ONLY : noghel,       enbhel,    LibInfo
IMPLICIT NONE
INTEGER :: IO
REAL, PARAMETER :: e0 = 1.0E-4_8
INTEGER :: ig

WRITE(IO, *)
WRITE(io, '(2x, A)') '=Multi-Group Library Information='
WRITE(IO, *)
WRITE(IO, '(5x, A)') LibInfo(2)
WRITE(IO, *)
WRITE(IO, '(5x, A)') '- Energy Group Structures -'
write(io,*) '       Group    Higher Bound     Lower Bound'
DO ig = 1, noghel-1
  WRITE(IO,'(10x, i3, 1p2e15.5)') ig, enbhel(ig), enbhel(ig+1)
ENDDO
WRITE(IO, '(10x, i3, 1p2e15.5)') noghel, enbhel(ig), e0
WRITE(IO, *)
END SUBROUTINE

!SUBROUTINE PrintInt1DarrayTo2Darray(io, array1D, nxy, nx, ny)
!USE PARAM
!IMPLICIT NONE
!INTEGER :: Array1D(nxy)
!INTEGER :: io, nxy, nx, ny
!INTEGER :: i, j,jbeg, jend
!DO i = 1, ny
!  jbeg = nx * (i - 1) + 1; jend = nx * i
!  WRITE(io, '(7x, 200I5)') (Array1D(j), j = jbeg, jend)
!ENDDO
!END SUBROUTINE
