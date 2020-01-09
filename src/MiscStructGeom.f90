#include <defines.h>
SUBROUTINE SetCoreMiscStruct()
USE PARAM
USE CNTL,       ONLY : nTracerCntl
USE geom,       ONLY : Core
USE Core_mod,   ONLY : FmInfo
USE PE_Mod,     ONLY : PE
IMPLICIT NONE

return
IF(nTracerCntl%lAutoBarrel) CALL SetCoreAnnularRingStruct()
!IF(nTracerCntl%lInitBoron) CALL SetBoronCoolant(Core, FmInfo%Fxr, nTracerCntl%BoronPPM , PE%myzb, PE%myze)

END SUBROUTINE
  
SUBROUTINE SetCoreAnnularRingStruct()
USE PARAM
USE TYPEDEF,       ONLY : MiscStruct_Type
USE GEOM,          ONLY : MiscStruct,                                             &
                          lCoreAng,         lEdge,                                &
                          nCellX,           nxy,            nx,          ny,      &
                          Asy,              Pin,            CellInfo
USE Core_mod,      ONLY : GroupInfo
USE CNTL,          ONLY : nTracerCntl
USE PE_Mod,        ONLY : PE
USE IOUTIL,        ONLY : Terminate
USE BasicOperation,ONLY : CP_CA
IMPLICIT NONE
LOGICAL :: lBarrel
LOGICAL :: lRing

!Cell Cordinate Information
INTEGER, POINTER :: PinIdx(:, :)
REAL, POINTER :: Cor_X(:), Cor_Y(:)
REAL :: Vx(2), Vy(2)
INTEGER :: ix, iy, ixy

!Ring Information
REAL :: CX, CY
REAL :: rin, rout
REAL :: VolFrac

INTEGER, PARAMETER :: OutCell = 0
INTEGER, PARAMETER :: InCell = 1
INTEGER, PARAMETER :: BoundCell_In = 2
INTEGER, PARAMETER :: BoundCell_Out = 3
INTEGER, PARAMETER :: BoundCell_InOut = 4

!Ring Location
INTEGER, POINTER :: RingCell(:, :)
INTEGER :: iringbeg, iringend
INTEGER :: ir, RingCellType

IF(.NOT. nTracerCntl%lMiscStruct) RETURN
lBarrel = nTracerCntl%lAutoBarrel; lRing = nTracerCntl%lAutoRingStruct

iringbeg = 1; iringend = 0
IF(lBarrel) iringbeg = 0
IF(lRing) iringend = MiscStruct%nring

IF(lCoreAng .EQ. 360 .AND. nx .NE. ny) THEN
  IF(lBarrel) CALL Terminate('BARREL CARD is not available : nx is not equal to ny')
  IF(lRing) CALL Terminate('BARREL CARD is not available : nx is not equal to ny')
ENDIF

ALLOCATE(PinIdx(nx, ny), Cor_X(0:nx), Cor_Y(0:ny))
ALLOCATE(RingCell(nx, ny))

CALL SetPinIdx(PinIdx, nx, ny, nxy)
!Set Cell Cordinate Info
CALL SetCelCorInfo(Cor_X, Cor_Y, PinIdx, nx, ny) 
!Set Ring Center
SELECT CASE(lCoreAng)
  CASE(360)
    CX = 0.5_8*(Cor_X(0) + Cor_X(nx))
    CY = 0.5_8*(Cor_Y(0) + Cor_Y(ny))
  CASE(90)
    CX= 0; CY = 0
END SELECT

CALL CP_CA(RingCell, 0, nx, ny)
!CALL CP_CA(RingVolFrac(1:nx, 1:ny), 0._8, nx, ny)
DO ir = iringbeg, iringend
  rin = MiscStruct%rad_ring(1, ir); rout = MiscStruct%rad_ring(2, ir)
  DO ixy = 1, nxy
    ix = Pin(ixy)%ix; iy = Pin(ixy)%iy
    Vx(1:2) = Cor_X(ix-1:ix); Vy(1:2) = Cor_Y(iy-1:iy)
    CALL CheckRingCell(Vx, Vy, Rin, Rout, Cx, Cy, RingCellType)
    IF(RingCellType .EQ. OutCell) CYCLE
    !RingCell(ix, iy) = RingCellType
    IF(RingCellType .EQ. InCell) THEN
      RingCell(ix, iy) = ir+1
      !RingVolFrac(ix, iy) = 1
    ELSEIF(RingCellType .EQ. BoundCell_InOut) THEN
      RingCell(ix, iy) = ir+1
    ELSE  
      CALL CalRingVolFrac(VolFrac, Vx, Vy, Rin, Rout, Cx, Cy, RingCellType)
      IF(VolFrac .GT. 0.5_8) RingCell(ix, iy) = ir+1
    ENDIF
  ENDDO
  IF(MiscStruct%lVolCor_ring(ir)) THEN
    CAll AnnualVolCorrection(ir, RingCell, nx, ny)
  ELSE
    MiscStruct%VolCor_ring(ir) = 1._8
  ENDIF
  
ENDDO

CALL SetAnnualMaterial(RingCell, nx, ny)
!STOP
END SUBROUTINE




SUBROUTINE SetAnnualMaterial(RingCell, nx0, ny0)
USE PARAM
USE TYPEDEF,      ONLY : MiscStruct_Type,  Asy_Type,      Pin_Type,    Cell_Type, &
                         FxrInfo_Type
USE GEOM,         ONLY : MiscStruct,                                              &
                         lCoreAng,         lEdge,                                 &
                         nCellX,           nxy,            nx,          ny,       &
                         Asy,              Pin,            CellInfo
USE Core_mod,       ONLY : Fxr,              GroupInfo
USE Material_mod,   ONLY : Mixture
USE CNTL,           ONLY : nTRACERCntl
USE PE_Mod,         ONLY : PE
USE allocs
USE BasicOperation, ONLY : CP_VA,         MULTI_CA
IMPLICIT NONE
INTEGER, INTENT(IN) :: nx0, ny0
INTEGER, INTENT(IN) :: RingCell(nx0, ny0)
!INTEGER :: ixy ,iz, imix, ir, izbeg, izend

INTEGER :: ir, iringbeg, iringend
LOGICAL :: lBarrel, lRing

INTEGER :: iz, izbeg, izend
INTEGER :: myzb, myze

INTEGER :: ix, iy, ixy, ifxr, ifxr0, icel, imix
INTEGER :: FxrIdxSt, nLocalFxr

INTEGER :: ntiso, niso

LOGICAL :: lRes, lXsLib

lBarrel = nTracerCntl%lAutoBarrel; lRing = nTracerCntl%lAutoRingStruct
iringbeg = 1; iringend = 0
IF(lBarrel) iringbeg = 0
IF(lRing) iringend = MiscStruct%nring

lRes = nTracerCntl%lRestrmt
lXsLib = nTracerCntl%lXsLib

myzb = PE%myzb; myze = PE%myze
ntiso = GroupInfo%ntiso
DO ir = iringbeg, iringend
  izbeg = MiscStruct%ring_plnbeg(ir); izend = MiscStruct%ring_plnend(ir)
  DO iz = izbeg, izend
    IF(iz .LT. myzb .OR. iz .GT. myze) CYCLE
    DO ixy = 1, nxy
      ix = Pin(ixy)%ix; iy = Pin(ixy)%iy
      IF(RIngCell(ix, iy) .NE. ir + 1) CYCLE
      icel = Pin(ixy)%Cell(iz); FxrIdxSt = Pin(ixy)%FxrIdxSt
      nLocalFxr = CellInfo(icel)%nFxr
      DO ifxr0 = 1, nLocalFxr
        ifxr = FxrIdxSt + ifxr0 - 1
        imix = MiscStruct%mix_ring(ir)
        Fxr(ifxr, iz)%imix = imix
        
        IF(lXsLib) THEN
          IF(lRes) Fxr(ifxr, iz)%lRes = Mixture(imix)%lRes
          niso =  Mixture(imix)%niso
          Fxr(ifxr, iz)%temp = Mixture(imix)%temp
          Fxr(ifxr, iz)%lh2o = Mixture(imix)%lh2o
          IF(Fxr(ifxr, iz)%ndim .LT. Min(ntiso,niso+4)) THEN
            Fxr(ifxr, iz)%ndim = Min(ntiso,niso+4)  
            DEALLOCATE(Fxr(ifxr,iz)%pnum, Fxr(ifxr, iz)%idiso)
            ALLOCATE(Fxr(ifxr,iz)%pnum(Fxr(ifxr, iz)%ndim))
            ALLOCATE(Fxr(ifxr,iz)%idiso(Fxr(ifxr, iz)%ndim))
          ELSE
            Fxr(ifxr, iz)%ndim = Min(ntiso,niso+4)  
          ENDIF
          Fxr(ifxr, iz)%ldepl = Mixture(imix)%ldepl
          CALL CP_VA(Fxr(ifxr, iz)%idiso(1:niso), Mixture(imix)%idiso(1:niso), niso)
          CALL CP_VA(Fxr(ifxr, iz)%pnum(1:niso), Mixture(imix)%pnum(1:niso), niso)
          IF(MiscStruct%lVolCor_ring(ir)) THEN
            CALL MULTI_CA(MiscStruct%VolCor_ring(ir), Fxr(ifxr, iz)%pnum(1:niso), niso)
          ENDIF
        ENDIF
        
      ENDDO
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE AnnualVolCorrection(ir, RingCell, nx0, ny0)
USE PARAM
USE TYPEDEF,       ONLY : Cell_Type,     Pin_Type
USE GEOM,          ONLY : nxy,            ny,               nx,         &
                          lCoreAng,                                    &
                          MiscStruct,     CellInfo,      Pin
IMPLICIT NONE
INTEGER :: RingCell(nx0, ny0)
INTEGER :: ir
INTEGER :: nx0, ny0
INTEGER :: ixy, ix, iy, icel
REAL :: CelVol
REAL :: Vol0, Vol1
Vol0 = (MiscStruct%rad_ring(2, ir)**2 - MiscStruct%rad_ring(1, ir)**2 ) * pi
IF(lCoreAng .EQ. 90) Vol0 = 0.25_8 * Vol0
Vol1 = 0
DO ixy = 1, nxy
  ix = Pin(ixy)%ix; iy = Pin(ixy)%iy
  IF(RingCell(ix, iy) .NE. ir+1) CYCLE
  icel = PIn(ixy)%Cell(1)
  CelVol = CellInfo(icel)%Geom%lx * CellInfo(icel)%Geom%ly
  Vol1 = Vol1 + CelVol
ENDDO
MiscStruct%VolCor_ring(ir) = Vol0 / Vol1
END SUBROUTINE

SUBROUTINE CalRingVolFrac(VolFrac, Vx, Vy, Rin, Rout, Cx, Cy, RingCellTyp)
USE PARAM
USE GeomTreatment,          ONLY : Line_Circle_Intersection
IMPLICIT NONE
REAL, INTENT(IN) :: Vx(2), Vy(2), Rin, Rout, Cx, Cy
INTEGER, INTENT(In) :: RingCellTyp
REAL, INTENT(OUT) :: VolFrac

INTEGER, PARAMETER :: OutCell = 0
INTEGER, PARAMETER :: InCell = 1
INTEGER, PARAMETER :: BoundCell_In = 2
INTEGER, PARAMETER :: BoundCell_Out = 3
INTEGER, PARAMETER :: BoundCell_InOut = 4

INTEGER, PARAMETER :: Out = 0
INTEGER, PARAMETER :: In = 1
INTEGER, PARAMETER :: Intsect = 2

INTEGER :: ir, iv, iv1, iv2, is, ipt
INTEGER :: npt(2), npt0, m, n
LOGICAL :: lInts(4, 2)
REAL :: IntsPt(2, 4, 2)

INTEGER :: VInOut(4,2), nVIn(2)

INTEGER :: Vneigh(2,4), V2L(2, 4), L2V(2, 4)
LOGICAL :: IntsectCrit(0:1,0:1)
REAL :: vertex(2, 4), line(3, 4), circle(3, 2), Sol(2, 2)
REAL :: r0(2)
REAL :: Dist, r, Vol0
REAL :: Vfrac(2)

DATA Vneigh / 2, 3, 1, 3, 1, 4, 2, 3 / 
DATA V2L    / 3, 2, 3, 4, 2, 1, 4, 1 /
DATA L2V    / 3, 4, 1, 3, 1, 2, 2, 4 /
DATA IntSectCrit / .FALSE., .TRUE., .TRUE., .FALSE./ 
!Vertex and Line Index
!   1    (3)    2
!    ***********
!    *         *
! (2)*         *(4)
!    *         *
!    ***********
!   3    (1)    4
vertex = reshape((/Vx(1), Vy(1), Vx(2), Vy(1), Vx(1), Vy(2), Vx(2), Vy(2)/), (/2,4/))
r0 = (/Rin, Rout/)
Line(:, 1) = (/0._8, 1._8, -Vy(2)/)
Line(:, 2) = (/1._8, 0._8, -Vx(1)/)
Line(:, 3) = (/0._8, 1._8, -Vy(1)/)
Line(:, 4) = (/1._8, 0._8, -Vx(2)/)

Vol0 = abs(Vx(1) -Vx(2)) * abs(Vy(1) - Vy(2))

Circle(:, 1) = (/CX, CY, Rin/)
Circle(:, 2) = (/CX, CY, Rout/)

VInOut = OUT; nVIn = 0
DO ir = 1, 2
  r = R0(ir)
  DO iv = 1, 4
    dist = sqrt((vertex(1, iv)-CX)**2+(vertex(2, iv)-CY)**2)
    IF(dist .LT. r) THEN
      VInOut(iv, ir) = IN
      nVIn(ir) = nVIn(ir) + 1
    ENDIF
  ENDDO
ENDDO


Vfrac = 1; IntsPt = 0; lInts = .FALSE.; npt = 0
DO ir = 1, 2
  r = R0(ir)
  IF(ir .EQ. 1 .AND. RingCellTyp .EQ. BoundCell_out) THEN
    Vfrac(1) = 1.
    CYCLE
  ENDIF  
  IF(ir.EQ. 2 .AND. RingCellTyp .EQ. BoundCell_In) THEN
    Vfrac(2) = 0.
    CYCLE
  ENDIF
  DO is = 1, 4
    iv1 = L2V(1, is); iv2 = L2V(2, is)
    IF(.NOT. IntsectCrit(VInOut(iv1, ir), VInOut(iv2, ir))) CYCLE
    Sol = Line_Circle_Intersection(Line(:, is), Circle(:, ir), npt0)
    m = 0
    DO ipt = 1, npt0
      dist = (Vertex(1, iv1)-Sol(1, ipt)) * (Vertex(1, iv2)-Sol(1, ipt))
      IF(dist .GT. EPSM5) CYCLE
      dist = (Vertex(2, iv1)-Sol(2, ipt)) * (Vertex(2, iv2)-Sol(2, ipt))
      IF(dist .GT. EPSM5) CYCLE
      m = m + 1
      EXit
    ENDDO
    npt(ir) = npt(ir) + 1
    lInts(is, ir) = .TRUE.
    IntsPt(:, is, ir) = Sol(:, ipt)
  ENDDO
  SELECT CASE(nVin(ir))
    CASE(1)
      Vfrac(ir) = VolumeCase1(lInts(:, ir), IntsPt(:, :, ir), VInOut(:, ir))
    CASE(2)
      Vfrac(ir) = VolumeCase2(lInts(:, ir), IntsPt(:, :, ir), VInOut(:, ir))
    CASE(3)
      Vfrac(ir) = VolumeCase3(lInts(:, ir), IntsPt(:, :, ir), VInOut(:, ir))
  END SELECT
  Vfrac(ir) = Vfrac(ir) / Vol0
ENDDO

!INTEGER, PARAMETER :: OutCell = 0
!INTEGER, PARAMETER :: InCell = 1
!INTEGER, PARAMETER :: BoundCell_In = 2
!INTEGER, PARAMETER :: BoundCell_Out = 3
!INTEGER, PARAMETER :: BoundCell_InOut = 4
SELECT CASE(RingCellTyp)
CASE(OutCell)
  VolFrac = 0
CASE(InCell)
  VolFrac = 1
CASE(BoundCell_In)
  VolFrac = 1._8 - Vfrac(1)
CASE(BoundCell_out)
  VolFrac = Vfrac(2)
CASE(BoundCell_InOut)
  VolFrac = Vfrac(2) - Vfrac(1)
END SELECT

CONTAINS
FUNCTION VolumeCase1(lInts0, IntsPt0, VInOut0)
IMPLICIT NONE
LOGICAL, INTENT(IN) :: lInts0(4)
REAL, INTENT(IN) :: IntsPt0(2, 4)
INTEGER, INTENT(IN) :: VInOut0(4)
REAL :: VolumeCase1
REAL :: l(2)
INTEGER :: is0, is1, iv0
DO iv0 = 1, 4
  IF(VInOut0(iv0) .EQ. OUT) CYCLE
  DO is1 = 1, 2
    is0 = V2L(is1, iv0)
    l(is1) = (IntsPt0(1, is0) - Vertex(1, iv0))**2 + (IntsPt0(2, is0) - Vertex(2, iv0))**2 
    l(is1) = SQRT(l(is1))
  ENDDO
  EXIT
ENDDO
VolumeCase1 = 0.5_8 * l(1) * l(2)

END FUNCTION

FUNCTION VolumeCase2(lInts0, IntsPt0, VInOut0)
IMPLICIT NONE
LOGICAL, INTENT(IN) :: lInts0(4)
REAL, INTENT(IN) :: IntsPt0(2, 4)
INTEGER, INTENT(IN) :: VInOut0(4)
REAL :: VolumeCase2
INTEGER :: is0, is1, iv0, iv1, iv2, i0
REAL :: l(0:2)

DO is0 = 1, 4
  iv1 = L2V(1, is0); iv2 = L2V(2, is0)
  IF(VInOut0(iv1) .EQ. IN .AND. VInOut0(iv2) .EQ. IN) THEN
    l(0) = (Vertex(1, iv1) - Vertex(1, iv2))**2 + (Vertex(2, iv1) - Vertex(2, iv2)) **2
    l(0) = SQRT(l(0))
    EXIT
  ENDIF
ENDDO

i0 = 0
DO iv0 = 1, 4
  IF(VInOut0(iv0) .EQ. OUT) CYCLE
  DO is1 = 1, 2
    is0 = V2L(is1, iv0)
    IF(.NOT. lInts0(is0)) CYCLE
    i0 = i0 + 1
    l(i0) = (IntsPt0(1, is0) - Vertex(1, iv0))**2 + (IntsPt0(2, is0) - Vertex(2, iv0))**2 
    l(i0) = SQRT(l(i0))
  ENDDO
ENDDO
VolumeCase2 = 0.5_8 * l(0) * (l(1) + l(2))
CONTINUE
END FUNCTION

FUNCTION VolumeCase3(lInts0, IntsPt0, VInOut0)
IMPLICIT NONE
LOGICAL, INTENT(IN) :: lInts0(4)
REAL, INTENT(IN) :: IntsPt0(2, 4)
INTEGER, INTENT(IN) :: VInOut0(4)
REAL :: VolumeCase3
INTEGER :: VInout1(4)
INTEGER :: InvInOut(0:1)


REAL :: lx, ly, vol0
INTEGER :: iv0

lx = Vertex(1, 2) - Vertex(1, 1)
ly = Vertex(2, 1) - Vertex(2, 3)
vol0  = lx * ly

InvInOut(0) = 1; InvInout(1) = 0
DO iv0 = 1, 4
  VInOut1(iv0) = InvInOut(VInOut0(iv0))
ENDDO
VolumeCase3 = VolumeCase1(lInts0, IntsPt0, VInOut1)
VolumeCase3 = Vol0 - VolumeCase3
CONTINUE

END FUNCTION


END SUBROUTINE



SUBROUTINE CheckRingCell(Vx, Vy, Rin, Rout, Cx, Cy, CellType)
REAL, INTENT(IN) :: Vx(2), Vy(2), Rin, Rout, Cx, Cy
INTEGER, INTENT(OUT) :: CellType

INTEGER, PARAMETER :: OutCell = 0
INTEGER, PARAMETER :: InCell = 1
INTEGER, PARAMETER :: BoundCell_In = 2
INTEGER, PARAMETER :: BoundCell_Out = 3
INTEGER, PARAMETER :: BoundCell_InOut = 4

INTEGER, PARAMETER :: Out = 0
INTEGER, PARAMETER :: In = 1
INTEGER, PARAMETER :: Intsect = 2

INTEGER, PARAMETER :: Err = -1

REAL :: vertex(2, 4)
INTEGER :: RingStat(2), RingStatDat(0:2, 0:2)

REAL :: r0(2)
REAL :: r, dist
INTEGER :: ir, iv
INTEGER :: nv
DATA RingStatDat  /OutCell,          InCell,             BoundCell_Out, &
                   Err,              OutCell,            Err,           &
                   Err,              BOundCell_In,       BoundCell_InOut/

vertex = reshape((/Vx(1), Vy(1), Vx(2), Vy(1), Vx(1), Vy(2), Vx(2), Vy(2)/), (/2,4/))
r0 = (/Rin, Rout/)
DO ir = 1, 2
  r = r0(ir)
  nv = 0; RingStat(ir) = Intsect
  DO iv = 1, 4
    dist = sqrt((vertex(1, iv)-CX)**2+(vertex(2, iv)-CY)**2)
    IF(dist .LT. r) nv = nv + 1
  ENDDO
  IF(nv .EQ. 4) RingStat(ir) = In
  IF(nv .EQ. 0) RingStat(ir) = Out
ENDDO

CellType = RingStatDat(RingStat(2), RingStat(1))
CONTINUE
END SUBROUTINE

SUBROUTINE SetPinIdx(PinIdx, nx, ny, nxy)
USE PARAM
USE TYPEDEF,        ONLY : Pin_Type
USE GEOM,           ONLY : Pin
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
INTEGER, INTENT(IN) :: nx, ny, nxy
INTEGER, INTENT(INOUT) :: PinIdx(nx, ny)
INTEGER :: ixy, ix, iy

CALL CP_CA(PinIdx, 0, nx, ny)
DO ixy = 1, nxy
  ix = Pin(ixy)%ix; iy = Pin(ixy)%iy
  PinIdx(ix, iy) = ixy
ENDDO
END SUBROUTINE

SUBROUTINE SetCelCorInfo(Cor_X, Cor_Y, PinIdx, nx, ny)
USE PARAM
USE TYPEDEF,       ONLY : Cell_Type,        Pin_Type
USE GEOM,          ONLY : CellInfo,         Pin
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE

INTEGER, INTENT(IN) :: PinIdx(nx, ny)
INTEGER, INTENT(IN) :: nx, ny
REAL, INTENT(INOUT) :: Cor_X(0:nx), Cor_Y(0:ny)

INTEGER :: ix, iy, ixy, icel

CALL CP_CA(Cor_X(0:nx), 0._8, nx+1)
CALL CP_CA(Cor_Y(0:ny), 0._8, ny+1)

DO ix = 1, nx
  DO iy = 1, ny
    IF(PinIdx(ix, iy) .NE. 0) EXIT
  ENDDO
  ixy = PinIdx(ix, iy)
  icel = Pin(ixy)%Cell(1)
  Cor_X(ix) = Cor_X(ix-1) + CellInfo(icel)%Geom%lx
  
ENDDO

DO iy = 1, ny
  DO ix = 1, nx
    IF(PinIdx(ix, iy) .NE. 0) EXIT
  ENDDO
  ixy = PinIdx(ix, iy)
  icel = Pin(ixy)%Cell(1)
  Cor_Y(iy) = Cor_Y(iy-1) - CellInfo(icel)%Geom%ly
ENDDO

END SUBROUTINE

!Subroutines for configuring the structure geometry = Baffle & Barrel

!SUBROUTINE InitStructCell(No, pitch, matid)
!USE PARAM
!IMPLICIT NONE
!
!REAL :: Pitch
!INTEGER :: No, matid
!
!CHARACTER(256) :: dataline
!REAL :: del
!INTEGER :: i
!
!Del = Pitch / 7
!!1111 format('(i3, 6F10.5, A, 6F10.5, A, 49I3, A, A5, A, A5)')
!!WRITE(dataline, 1111) No, (Del, i=1,6), '/', (Del, i=1,6), '/', (matid, i=1,49),'/', '7*1','/', '7*1'
!
!1001 format(i3, 6F10.5, 1x, A, 6F10.5, 1x, A, 1x, A3, I1,1x, A, A5, 1x, A, A5)
!1002 format(i3, 6F10.5, 1x, A, 6F10.5, 1x, A, 1x, A3, I2,1x, A, A5, 1x, A, A5)
!1003 format(i3, 6F10.5, 1x, A, 6F10.5, 1x, A, 1x, A3, I3,1x, A, A5, 1x, A, A5)
!IF(matid .LT. 10) THEN
!  WRITE(dataline, 1001) No, (Del, i=1,6), '/', (Del, i=1,6), '/', '49*',matid, '/', '7*1','/', '7*1'
!ELSEIF(matid .LT. 100) THEN
!  WRITE(dataline, 1002) No, (Del, i=1,6), '/', (Del, i=1,6), '/', '49*',matid, '/', '7*1','/', '7*1'
!ELSE 
!  WRITE(dataline, 1003) No, (Del, i=1,6), '/', (Del, i=1,6), '/', '49*',matid, '/', '7*1','/', '7*1'
!ENDIF
!CALL ReadNInit_Cell(dataline, FALSE, FALSE)
!END SUBROUTINE

!SUBROUTINE SetBaffleCell()
!USE PARAM
!USE ALLOCS
!!USE CNTL, ONLY : master
!USE GEOM, ONLY : nAsyType,   nPinType,    nCellType,    nPinType0,     nCellType0,         &  
!                 Core,       AsyInfo,     PinInfo,      CellInfo,      Pin,                &
!                 Asy,        hz,          nz,           nzfm,                              &
!                 hzInv,      HzFmInv,                                                      &
!                 nCellX0,    nCellX,                                                       &
!                 AsyPitch,   lEdge,       lGap,                                            &
!                 hzfm,       SubPlaneMap, SubPlaneRange, nSubPlane
!USE PE_MOD,         ONLY : PE            
!USE CNTL,           ONLY : nTracerCntl
!IMPLICIT NONE
!INTEGER :: iasytype, ixy, ibd, ineigh, iz, icel, istructcell
!INTEGER :: nxy, nbd, nFsr, nFxr
!
!nxy = Core%nxy; nbd = 4
!istructcell = nTracerCntl%StructCell
!nFsr = CellInfo(iStructCell)%nFsr; nFxr = CellInfo(iStructCell)%nFxr
!DO ixy = 1, nxy
!  !IF(.NOT. Pin(ixy)%lRadRef) CYCLE
!  iasytype=Pin(ixy)%AsyType
!  IF(AsyInfo(iasyType)%lfuel) CYCLE
!  DO ibd = 1, 4
!    ineigh = Pin(ixy)%NeighIdx(ibd)
!    IF(ineigh .LT. 1) CYCLE
!    !IF(.NOT. Pin(ineigh)%lfuel) CYCLE
!    iasytype=Pin(ineigh)%AsyType
!    IF(.NOT. AsyInfo(iasyType)%lfuel) CYCLE
!    Pin(ixy)%lBaffle = .TRUE.
!    DO iz = 1, nz
!      icel = Pin(ixy)%cell(iz)
!      IF(.NOT. Core%lFuelPlane(iz)) CYCLE
!      IF(.NOT. lGap) THEN
!        Pin(ixy)%Cell(iz) = istructcell
!        Pin(ixy)%nFsrMax = MAX(Pin(ixy)%nFsrMax, nFSR)
!        Pin(ixy)%nFxrMax = MAX(Pin(ixy)%nFxrMax, nFXR)
!      ENDIF
!    ENDDO
!  ENDDO
!ENDDO
!END SUBROUTINE 