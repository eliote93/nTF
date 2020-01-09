SUBROUTINE MakeCellRay1DArray(RayInfo, CellRayBase, CellInfo, nCell)
USE PARAM
USE TYPEDEF,     ONLY : CellRayBase_Type, CellRayInfo_Type, RayInfo_Type, &
                        Cell_Type
USE allocs
IMPLICIT NONE

TYPE(RayInfo_Type) :: RayInfo
TYPE(CellRayBase_Type) :: CellRayBase(nCell)
TYPE(CellRayInfo_Type), POINTER :: CellRay1D
Type(Cell_Type) :: CellInfo(ncell)
INTEGER :: ncell

INTEGER :: icel, iseg, iline, idx, iazi, ibeg, iend
INTEGER :: nline_tot, nseg_tot, nline, nseg
INTEGER :: ibcel

ALLOCATE(RayInfo%CellRay1D)
CellRay1D => RayInfo%CellRay1D
nline_tot = 0; nseg_tot = 0
DO icel = 1, ncell
  ibcel=CellRayBase(icel)%CellType   !--- BYS edit : same base cell ray structure
  IF(.NOT. CellInfo(icel)%luse) CYCLE
  nline_tot = nline_tot + CellRayBase(ibcel)%nline
  nline = CellRayBase(ibcel)%nline
  DO iline = 1, CellRayBase(ibcel)%nline
    nseg_tot = nseg_tot + CellRayBase(ibcel)%CellRay(iline)%nseg
  ENDDO
ENDDO

CellRay1D%nseg = nseg_Tot
CALL DMALLOC(CellRay1D%LocalFsrIdx, nseg_Tot)
CALL DMALLOC(CellRay1D%LenSeg, nseg_Tot)
idx = 0

DO iazi = 1, RayInfo%nAziAngle
  DO icel = 1, ncell
    ibcel=CellRayBase(icel)%CellType
    IF(.NOT. CellInfo(icel)%luse) CYCLE
    nline = CellRayBase(ibcel)%nline
    ibeg =CellRayBase(ibcel)%AziAngIdx(iazi); iend = CellRayBase(ibcel)%AziAngIdx(iazi)+CellRayBase(ibcel)%nlines(iazi)-1
    
    DO iline = ibeg, iend 
      nseg = CellRayBase(ibcel)%CellRay(iline)%nseg
      CellRayBase(ibcel)%CellRay(iline)%idxst= idx  + 1
      DO iseg = 1, nseg
        idx = idx + 1
        CellRay1D%LenSeg(idx) = CellRayBase(ibcel)%CellRay(iline)%lenseg(iseg)
        CellRay1D%LocalFsrIdx(idx) = CellRayBase(ibcel)%CellRay(iline)%LocalFsrIdx(iseg)
      ENDDO
      CellRayBase(ibcel)%CellRay(iline)%idxend = idx
    ENDDO
  ENDDO
ENDDO

!RayInfo%CellRay1D => CellRay1D
END SUBROUTINE
  
SUBROUTINE SetFastMOC(RayInfo, Core, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,      ONLY : RayInfo_Type,      coreinfo_type,   FastCoreRayDat_Type,                &
                         FastRaySegDat_Type,                                                     & 
                         Pin_Type,          Asy_Type,        AsyInfo_Type,     PinInfo_Type,     &
                         Cell_Type,                                                              &
                         AziAngleInfo_Type, PolarAngle_Type, ModRayInfo_type,  AsyRayInfo_type,  &
                         CoreRayInfo_Type,  RotRayInfo_Type, CellRayInfo_type,                   &
                         PE_TYPE
USE Moc_Mod,      ONLY : nMaxRaySeg,        nMaxCellRay,     nMaxAsyRay,       nMaxCoreRay
USE CNTL,         ONLY : nTracerCntl_Type
USE ALLOCS
IMPLICIT NONE
TYPE(RayInfo_Type), INTENT(INOUT) :: RayInfo
TYPE(CoreInfo_Type), INTENT(INOUT) :: Core
!TYPE(FastCoreRayDat_Type), INTENT(INOUT), POINTER :: FastCoreRayDat(:, :)
TYPE(nTracerCntl_Type), INTENT(INOUT) :: nTracerCntl
TYPE(PE_Type) :: PE

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Asy_Type), POINTER :: Asy(:)
!TYPE(AsyInfo_Type), POINTER :: AsyInfo
TYPE(PinInfo_Type), POINTER :: PinInfo(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(AziAngleInfo_Type), POINTER :: AziAng(:)
TYPE(PolarAngle_Type), POINTER :: PolarAng(:)
!TYPE(ModRayInfo_type), POINTER :: ModRay
TYPE(AsyRayInfo_type), POINTER :: AsyRay(:)
TYPE(CoreRayInfo_Type), POINTER :: CoreRay(:)
TYPE(RotRayInfo_Type), POINTER :: RotRay(:)
TYPE(CellRayInfo_Type), POINTER :: CellRay


TYPE(FastCoreRayDat_Type), POINTER :: FastCoreRayDat(:, :)
REAL, POINTER :: LenSeg(:)
INTEGER, POINTER :: LocalFsrIdx(:)
INTEGER, POINTER :: FsrIdx(:, :)

INTEGER :: i, j, k, l, m
INTEGER :: iasy, ipin, icel, iz, ireg
INTEGER :: iCoreRay, iRotRay, iasyray, iceray
INTEGER :: irsegidx, icellrayidx, itype, irayseg
INTEGER :: FsrIdxSt 
INTEGER :: nCoreRay, nAsyRay, nPinRay, nRaySeg
INTEGER :: ibcel

INTEGER, POINTER :: CellRayIdxSt(:, :, :), PinIdx(:, :), SurfIdx(:, :, :), Ray1DIdx(:, :, :), CeRayIdx(:, :)
INTEGER, POINTER :: nTotRaySeg(:), nTotCellRay(:)

!nTracerCntl%FastMocLv = 1

ALLOCATE(CellRayIdxSt(nMaxCellRay, nMaxCoreRay, 2))
ALLOCATE(PinIdx(nMaxCellRay, nMaxCoreRay))
ALLOCATE(SurfIdx(nMaxCellRay, nMaxCoreRay, 2))
ALLOCATE(nTotRaySeg(nMaxCoreRay))
ALLOCATE(nTotCellRay(nMaxCoreRay))
ALLOCATE(Ray1DIdx(2, nMaxCellRay, nMaxCoreRay))
ALLOCATE(CeRayIDx(nMaxCellRay, nMaxCoreRay))

ALLOCATE(RayInfo%FastCoreRayDat(RayInfo%nRotRay, PE%myzb:PE%myze))
FastCoreRayDat => RayInfo%FastCoreRayDat

!Ray Info Pointing
AziAng => RayInfo%AziAngle;
PolarAng => RayInfo%PolarAngle;
AsyRay => RayInfo%AsyRay
CoreRay => RayInfo%CoreRay
RotRay => RayInfo%RotRay

!Geometry Info Pointing
Asy => Core%Asy
Pin => Core%Pin
PinInfo => Core%Pininfo
Cell => Core%CellInfo
!ALLOCATE(FastCoreRayDat(RayInfo%nRotRay, PE%myzb:PE%myze))
DO iz = PE%myzb, PE%myze
  DO i = 1, RayInfo%nRotRay
    irotRay= i
    nCoreRay = RotRay(irotRay)%nRay
    nTotCellRay = 0; nTotRaySeg =0
    DO j = 1, nCoreRay    !Core Ray Sweep
      irsegidx = 0;   icellrayidx = 0
      iCoreRay = RotRay(iRotRay)%RayIdx(j)
      nAsyRay = CoreRay(iCoreRay)%nRay
      DO k = 1, nAsyRay    !Assembly Ray Sweep
        iasyray = CoreRay(iCoreRay)%AsyRayIdx(k)
        iasy = CoreRay(iCoreRay)%AsyIdx(k)
        IF(iasy .EQ. 0)  CYCLE   !Skip Dummy Assembly
        nPinRay = AsyRay(iAsyRay)%nCellRay
        itype = Asy(iasy)%PartialAsyFlag
        DO l = 1, nPinRay   !Pin Ray Sweep
          ipin = AsyRay(iAsyRay)%PinIdx(l)      !Local Pin Idx(within Assembly)
          iceray = AsyRay(iAsyRay)%PinRayIdx(l) !Cell Ray Index
          ipin = Asy(iAsy)%GlobalPinIdx(ipin)   !Global Pin Index
          icel = Pin(ipin)%Cell(iz)             !Cell Type
          FsrIdxSt = Pin(ipin)%FsrIdxSt
          ibcel=Cell(icel)%basecellstr
          CellRay => Cell(ibcel)%CellRay(iceray) !Pointing Cell Ray
          !CeRayIdx(icellrayidx, j) = iceray
          icellrayidx = icellrayidx + 1
          PinIdx(icellrayidx, j) = ipin
          CellRayIdxSt(icellrayidx, j, 2) = irsegidx + 1
          CeRayIdx(icellrayidx, j) = iceray
          nRaySeg = CellRay%nSeg
          LocalFsrIdx => CellRay%LocalFsrIdx
          LenSeg => CellRay%LenSeg
          ray1DIdx(1, icellrayidx, j) = CellRay%Idxst
          ray1DIdx(2, icellrayidx, j) = CellRay%Idxend
          DO iRaySeg = 1, nRaySeg
            !ireg = FsrIdxSt + LocalFsrIdx(iRaySeg) - 1
           ! ireg = FsrIdxSt + CellRay%LocalFsrIdx(iRaySeg) - 1 
            irsegidx = irsegidx + 1
          ENDDO   !End of Ray Segments Sweep, irayseg
          CellRayIdxSt(icellrayidx, j, 1) = irsegidx
          SurfIdx(icellRayIdx, j, 1) = AsyRay(iAsyRay)%PinRaySurf(2, l) !OutSurface
          SurfIdx(icellRayIdx, j, 2) = AsyRay(iAsyRay)%PinRaySurf(1, l) !Insurface
        ENDDO    !End of Pin Ray Seep, l
      ENDDO    !End of Asy Ray Swwep, k
      nTotRaySeg(j) = irsegidx
      nTotCellRay(j) = icellRayIdx
    ENDDO    !End of Core Ray Sweep, j 
    CALL DMALLOC(FastCoreRayDat(i,iz)%PinIdx, maxval(nTotCellRay),nCoreRay)
    CALL DMALLOC(FastCoreRayDat(i,iz)%CellRayIdxSt, maxval(nTotCellRay),nCoreRay, 2)
    CALL DMALLOC(FastCoreRayDat(i,iz)%SurfIdx, maxval(nTotCellRay),nCoreRay, 2)
    CALL DMALLOC(FastCoreRayDat(i,iz)%Ray1DIdx, 3, maxval(nTotCellRay),nCoreRay)
    CALL DMALLOC(FastCoreRayDat(i,iz)%nTotRaySeg, nCoreRay)
    CALL DMALLOC(FastCoreRayDat(i,iz)%nTotCellRay, nCoreRay)
    
    IF(nTracerCntl%FastMocLv .EQ. 2) THEN
      ALLOCATE(FastCoreRAyDat(i, iz)%RaySeg(nCoreRay))
    ENDIF
  
    
    DO j =1, nCoreRay
      FastCoreRayDat(i,iz)%nTotRaySeg(j) = nTotRaySeg(j)
      FastCoreRayDat(i,iz)%nTotCellRay(j) = nTotCellRay(j)
      DO l = 1, nTotCellRay(j)
        FastCoreRayDat(i,iz)%PinIdx(l, j) = PinIDx(l,j)
        FastCoreRayDat(i,iz)%CellRayIdxSt(l, j, 1) = CellRayIdxSt(l,j ,1)
        FastCoreRayDat(i,iz)%CellRayIdxSt(l, j, 2) = CellRayIdxSt(l,j ,2)
        FastCoreRayDat(i,iz)%SurfIdx(l, j, 1) = SurfIdx(l,j ,1)
        FastCoreRayDat(i,iz)%SurfIdx(l, j, 2) = SurfIdx(l,j ,2)
        FastCoreRayDat(i,iz)%Ray1DIdx(1:2, l, j) = Ray1DIdx(1:2, l,j)
        FastCoreRayDat(i,iz)%Ray1DIdx(3, l, j) = CeRayIdx(l,j)
      ENDDO
    ENDDO
    
    IF(nTracerCntl%FastMocLv .NE. 2) CYCLE
    
    DO j = 1, nCoreRay
      CALL Dmalloc(FastCoreRayDat(i, iz)%RaySeg(j)%FsrIdx, nTotRaySeg(j))
      CALL Dmalloc(FastCoreRayDat(i, iz)%RaySeg(j)%Lenseg, nTotRaySeg(j))
      iRaySeg = 0
      DO l = 1, nTotCellRay(j)
        ipin = PinIDx(l,j)
        FsrIdxSt = Pin(ipin)%FsrIdxSt
        DO m = Ray1DIdx(1, l, j), Ray1DIdx(2, l, j)
          iRaySeg = iRaySeg + 1
          ireg = FsrIdxSt + RayInfo%CellRay1D%LocalFsrIdx(m) - 1
          FastCoreRayDat(i, iz)%RaySeg(j)%FsrIdx(irayseg) = ireg
          FastCoreRayDat(i, iz)%RaySeg(j)%Lenseg(irayseg) = RayInfo%CellRay1D%LenSeg(m)
        ENDDO
      ENDDO
    ENDDO
    !CALL DMALLOC(FastCoreRayDat(i,iz)%CellIdx, maxval(nTotCellRay),nCoreRay)
  ENDDO
ENDDO


!RayInfo%FastCoreRayDat => FastCoreRayDat

END SUBROUTINE