#include <defines.h>
MODULE BinOutp_Mod
USE PARAM
USE TYPEDEF,         ONLY : CoreInfo_Type,        FmInfo_Type,         CmInfo_Type,        &
                            FxrInfo_Type,         Pin_Type,            THInfo_Type,        &
                            PinXS_Type,           Cell_Type
USE GEOM,            ONLY : Core,                 ng
USE Core_mod,        ONLY : FmInfo,               CmInfo,              eigv,               &
                            ThInfo,               GroupInfo
USE DEPL_MOD,        ONLY : DeplCntl
USE CNTL,            ONLY : nTracerCntl
USE PE_Mod,          ONLY : PE
USE FILES,           ONLY : CaseID,               io8,                 io17,               &
                            workingdir,           outputdir
USE IOUTIL,          ONLY : GETIFILE,             openfile,            PWD,                &
                            CreateDir,            GetOutputDir,        ChangeDir,          &
                            OpenFile
USE BasicOperation,  ONLY : DotProduct,           CP_CA,              CP_VA
USE MPIComm_Mod,     ONLY : BCAST,                REDUCE,             MPI_SYNC
IMPLICIT NONE

TYPE AsyOutp_Type
  SEQUENCE
  INTEGER :: ixa, iya, nx, ny, nfuel
  INTEGER :: itype
  INTEGER, POINTER :: PinIdx(:, :)
  INTEGER, POINTER :: PinType(:, :) !1: Fuel 2: Non-Fuel
END TYPE


LOGICAL, PRIVATE :: lBOUTP, lFLUX_BOUTP, lTH_BOUTP, lDEPL_BOUTP
LOGICAL, PRIVATE :: lFirst = .TRUE.
INTEGER, PARAMETER :: nMajorIsotope = 18
INTEGER, PRIVATE :: MajorIsotope(nMajorIsotope)
INTEGER, PRIVATE :: nglist(3)
REAL, PRIVATE :: GrpLBList(8,3)

INTEGER, PRIVATE :: ngc
INTEGER, PRIVATE :: Grpstruct(2,ngmax)
LOGICAL, PRIVATE :: MASTER
INTEGER, PRIVATE :: COMM, NPROC
INTEGER, PRIVATE :: myzb, myze, nz, nxy
INTEGER, PRIVATE :: CalNoId

TYPE(AsyOutp_Type), POINTER, PRIVATE :: AsyOutp(:)

TYPE(Pin_Type), POINTER, PRIVATE :: Pin(:)
TYPE(FxrInfo_Type), POINTER, PRIVATE :: Fxr(:, :)

TYPE(PinXS_Type), POINTER, PRIVATE :: PinXS(:, :)
TYPE(Cell_TYPE), POINTER, PRIVATE :: CellInfo(:)
REAL, POINTER, PRIVATE :: PHI(:, :, :), PSI(:, :), PHIC(:, :, :), PSIC(:, :)
REAL, POINTER :: PinVol(:, :)
REAL, PRIVATE :: NormPin3D, NormPin2D, NormAsy3D, NormAsy2D
REAL, PRIVATE :: NormAbsFlux
REAl, PRIVATE :: hactive

INTEGER, PRIVATE :: nIsoOut = 0
INTEGER, PRIVATE :: IsoOutList(1:1000) = 0
INTEGER, PRIVATE :: ntiso


DATA MajorIsotope /                                                          &
                    90232,      92233,      92234,      92235,      92236,   &
                    92238,      93237,      93239,      94233,      94239,   &
                    94240,      94241,      94242,      95241,      95242,   &
                    95242,      96242,      96244                            /
DATA nglist    /2, 4, 8/
DATA GrpLBList / 6.2506E-01,         0._8,       0._8,         0._8,       0._8,        0._8,         0._8,  0._8,     &
                 9.1188E+3_8, 3.9279E+0_8, 6.2506E-1_8,        0._8,       0._8,        0._8,         0._8,  0._8,     &
                 2.2313E+6_8, 8.2085E+5_8, 9.1188E+3_8, 1.3007E+2_8, 3.9279E+0_8, 6.2506E-1_8, 1.4572E-1_8, 0._8/

CONTAINS
SUBROUTINE BinaryOutput()
IMPLICIT NONE
LOGICAL :: mkdir, status

IF(lFirst) THEN
  CALL SetBinOutpCntl()
  CALL SetAsyOutp()
  IF(MASTER)  mkdir = CreateDir(caseid)
  OutputDir = GetOutputDir()
  lFirst = .FALSE.
ENDIF

CalNoId = nTracerCntl%CalNoId

CALL Normalization()

Status = ChangeDir(OutPutDir)

IF(CalNoId .EQ. 0) THEN
  IF(MASTER) CALL BinGeomEdit()
#ifdef MPI_ENV
  CALL MPI_SYNC(COMM)
#endif
ENDIF

IF(lFlux_BOUTP) CALL BinFluxEdit()

Status = ChangeDir(WorkingDir)
CONTINUE
END SUBROUTINE

SUBROUTINE BinGeomEdit()
IMPLICIT NONE
!Radial Information
CHARACTER(256) :: fn

INTEGER :: FuelPlane(nzmax), FuelPin(nzmax)
INTEGER :: ixa, iya, ixya, ix, iy, ix0, iy0, ixy, iz, ityp
INTEGER :: nx, ny

IF(MASTER) THEN
  CALL GetFileName(fn, 'GEOM', CalNoId)
  CALL Openfile(io17, .FALSE., .TRUE., .FALSE., Fn)
ENDIF

!Axial Information
DO iz = 1, nz
  FuelPlane(iz) = 0
  IF(Core%lFuelPlane(iz)) FuelPlane(iz) = 1
ENDDO
WRITE(io17) Core%nz
WRITE(io17) (REAL(Core%hz(iz), 4), iz = 1, nz)
WRITE(io17) (FuelPlane(iz), iz = 1, nz)
WRITE(io17) Core%nxa, Core%nya, Core%nxya
!Assembly Index
DO iy = 1, Core%nya
  WRITE(io17) (Core%CoreIdx(ix, iy), ix = 1, Core%nxa)
ENDDO

!Assembly Information
DO ixya = 1, Core%nxya
  nx = AsyOutp(ixya)%nx; ny = AsyOutp(ixya)%ny
  WRITE(io17) ixya, AsyOutp(ixya)%ixa, AsyOutp(ixya)%iya, AsyOutp(ixya)%nx, AsyOutp(ixya)%ny,  AsyOutp(ixya)%itype, AsyOutp(ixya)%nfuel
  DO iy = 1, ny
    WRITE(io17) (AsyOutp(ixya)%PinIdx(ix, iy), ix = 1, nx)
  ENDDO
  DO iy = 1, ny
    WRITE(io17) (AsyOutp(ixya)%PinType(ix, iy), ix = 1, nx)
  ENDDO
ENDDO

IF(MASTER) CLOSE(io17)

END SUBROUTINE

SUBROUTINE BinFluxEdit()
!USE BasicOperation
IMPLICIT NONE
!INTEGER :: Pw1D(nz), Psi1D(nz), Phi1D(nz, ngc)
INTEGER :: iz, ixya, ixa, iya, ix, iy, ixy, ig,  igb, ige
INTEGER :: nx, ny
REAL, POINTER :: AsyPw(:, :, :), AsyPsi(:, :, :), AsyPhi(:, :, :, :), AbsPhi(:, :, :, :)
REAL, POINTER :: AsyPw2D(:, :), AsyPsi2D(:, :), AsyPhi2D(:, :, :)
REAL, POINTER :: AsyBuf(:, :, :)
REAL :: f
CHARACTER(256) :: fn

IF(MASTER) THEN
  CALL GetFileName(fn, 'FLUX', CalNoId)
  CALL Openfile(io17, .FALSE., .TRUE., .FALSE., Fn)
ENDIF 

DO ixya = 1, Core%nxya
  nx = AsyOutp(ixya)%nx; ny = AsyOutp(ixya)%ny
  ALLOCATE(AsyPw(nz, nx, ny)); ALLOCATE(AsyPsi(nz, nx, ny))
  ALLOCATE(AsyPhi(nz, nx, ny, ngc))
  ALLOCATE(AsyBuf(nz, nx, ny))
  CALL CP_CA(AsyPsi, 0._8, nz, nx, ny); CALL CP_CA(AsyPw, 0._8, nz, nx, ny)
  CALL CP_CA(AsyPhi, 0._8, nz,  nx, ny, ngc)
  CALL CP_CA(AsyBuf, 0._8, nz, nx, ny)
  DO iy = 1, ny
    DO ix = 1, nx
      ixy = AsyOutp(ixya)%PinIdx(ix, iy)
      DO iz = myzb, myze
        AsyPsi(iz, ix, iy) = DotProduct(PinXs(ixy, iz)%XSNF(1:ng), PhiC(ixy, iz, 1:ng), ng) * PinVol(ixy, iz)
        AsyPw(iz, ix, iy) = DotProduct(PinXs(ixy, iz)%XSKF(1:ng), PhiC(ixy, iz, 1:ng), ng) * PinVol(ixy, iz)
        DO ig = 1, ngc
          igb = GrpStruct(1, ig); ige= GrpStruct(2, ig)
          AsyPhi(iz, ix, iy, ig) = SUM(PhiC(ixy, iz, igb:ige)) * PinVol(ixy, iz)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
#ifdef MPI_ENV
  CALL REDUCE(AsyPw, AsyBuf, nz, nx, ny, COMM, .FALSE.)
  CALL CP_VA(AsyPw, AsyBuf, nz, nx, ny)
  CALL REDUCE(AsyPsi, AsyBuf, nz, nx, ny, COMM, .FALSE.)
  CALL CP_VA(AsyPsi, AsyBuf, nz, nx, ny)
  DO ig = 1, ngc
    CALL REDUCE(AsyPhi(1:nz, 1:nx, 1:ny, ig), AsyBuf, nz, nx, ny, COMM, .FALSE.)
    CALL CP_VA(AsyPhi(1:nz, 1:nx, 1:ny, ig), AsyBuf, nz, nx, ny)
  ENDDO
#endif
  IF(.NOT. MASTER) THEN
    DEALLOCATE(AsyPw, AsyPsi, AsyPhi)
    DEALLOCATE(Asybuf)
    CALL MPI_SYNC(COMM)
    CYCLE
  ENDIF
  ALLOCATE(AsyPw2D(nx, ny))
  ALLOCATE(AsyPsi2D(nx, ny))
  ALLOCATE(AsyPhi2D(nx, ny, ngc))
  ALLOCATE(AbsPhi(nz, nx, ny, ngc))
  DO iy = 1, ny
    DO ix = 1, nx
      AsyPw2D(ix, iy) = SUM(AsyPw(1:nz, ix, iy))
      AsyPsi2D(ix, iy) = SUM(AsyPsi(1:nz, ix, iy))
      DO ig= 1, ngc
        AsyPhi2D(ix, iy, ig) = SUM(AsyPhi(1:nz, ix, iy, ig))
      ENDDO
    ENDDO
  ENDDO
  
  !Normalization
  DO iy = 1, ny
    DO ix = 1, nx
      ixy = AsyOutP(ixya)%PinIdx(ix, iy)
      f = NormPin2D / (PinVol(ixy, 1) / Core%hz(1))
      AsyPw2D(ix, iy) = AsyPw2D(ix, iy) * f
      AsyPsi2D(ix, iy) = AsyPsi2D(ix, iy) * f
      DO ig= 1, ngc
        AsyPhi2D(ix, iy, ig) = AsyPhi2D(ix, iy, ig) * f
      ENDDO
      DO iz = 1, nz
        f = NormPin3D / PinVol(ixy, iz)
        AsyPw(iz, ix, iy) = AsyPw(iz, ix, iy) * f
        AsyPsi(iz, ix, iy) = AsyPsi(iz, ix, iy) * f
        DO ig = 1, ngc
          AbsPhi(iz, ix, iy, ig) = AsyPhi(iz, ix, iy, ig) * NormAbsFlux / PinVol(ixy, iz)
          AsyPhi(iz, ix, iy, ig) = AsyPhi(iz, ix, iy, ig) * f
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  
  WRITE(io17) ixya, AsyOutp(ixya)%nx, AsyOutp(ixya)%ny,  AsyOutp(ixya)%ny, AsyOutp(ixya)%itype
  IF(AsyOutp(ixya)%itype .EQ. 1) THEN
    DO iy = 1, ny
      WRITE(io17) (REAL(AsyPw2D(ix, iy),4), ix = 1, nx)
    ENDDO
    DO iy = 1, ny
      WRITE(io17) (REAL(AsyPsi2D(ix, iy), 4), ix = 1, nx)
    ENDDO
  ENDIF
  DO ig = 1, ngc
    DO iy = 1, ny
      WRITE(io17) (REAL(AsyPhi2D(ix, iy, ig), 4), ix = 1, nx)
    ENDDO
  ENDDO
  DO iy = 1, ny
    DO ix = 1, nx
      IF(AsyOutp(ixya)%itype .EQ. 1) THEN
        WRITE(io17) (REAL(AsyPw(iz, ix, iy), 4), iz = 1, nz)
        WRITE(io17) (REAL(AsyPsi(iz, ix, iy), 4), iz = 1, nz)
      ENDIF

      DO ig = 1, ngc
        WRITE(io17) (REAL(AsyPhi(iz, ix, iy, ig), 4), iz = 1, nz)
      ENDDO
      DO ig = 1, ngc
        WRITE(io17) (REAL(AbsPhi(iz, ix, iy, ig), 4), iz = 1, nz)
      ENDDO
    ENDDO
  ENDDO
  DEALLOCATE(AsyPhi2D, AsyPw2D, AsyPsi2D)
  DEALLOCATE(AsyPw, AsyPsi, AsyPhi)
  DEALLOCATE(Asybuf)
  CALL MPI_SYNC(COMM)
ENDDO

IF(MASTER) CLOSE(io17)

END SUBROUTINE

SUBROUTINE BinDeplEdit()
USE XSLIB_MOD,      ONLY : mapnucl
IMPLICIT NONE
CHARACTER(256) :: fn
INTEGER :: i, j, k, iz, ixya, ixa, iya, ix, iy, ixy, icel
INTEGER :: FxrIdxSt, nLocalFxr, ifxr, iso, id
INTEGER :: iproc
INTEGER :: nx, ny

REAL :: volsum

REAL, POINTER :: Isotope(:, :), Isotope_Buf(:, :)

IF(MASTER) THEN
  CALL GetFileName(fn, 'FLUX', CalNoId)
  CALL Openfile(io17, .FALSE., .TRUE., .FALSE., Fn)
ENDIF

IF(MASTER) THEN
  WRITE(io17) nIsoOut
  WRITE(io17) (IsoOutList(i), i = 1, nIsoOut)
ENDIF

ALLOCATE(Isotope(0:ntiso, nz))
ALLOCATE(Isotope_buf(0:ntiso, nz))
DO ixya = 1, Core%nxya
  nx = AsyOutp(ixya)%nx; ny = AsyOutp(ixya)%ny
  !Non-Fuel Region
  IF(AsyOutp(ixya)%itype .EQ. 0) THEN
    IF(MASTER) WRITE(io17) ixya, AsyOutp(ixya)%nx, AsyOutp(ixya)%ny,  AsyOutp(ixya)%ny, AsyOutp(ixya)%itype
    CYCLE
  ENDIF
  DO iy = 1, ny
    DO ix = 1, nx
      ixy = AsyOutp(ixya)%PinIdx(ix, iy)
      FxrIdxSt = Pin(Ixy)%FxrIdxSt
      DO iz = 1, nz
        CALL CP_CA(Isotope(0:ntiso, iz), 0._8, ntiso+1)
        icel = Pin(ixy)%cell(iz)
        IF(MASTER) WRITE(io17) ix, iy, AsyOutp(ixya)%PinType(ix, iy)
        IF(.NOT. CellInfo(icel)%lfuel) CYCLE

        nLocalFxr = CellInfo(icel)%nFxr
        volsum = 0
        DO j = 1, nLocalFxr
          ifxr = FxrIdxSt + j - 1
          IF(.NOT. Fxr(ifxr, iz)%ldepl) CYCLE
          volsum = volsum + Fxr(ifxr, iz)%area
          Isotope(0, iz) = Isotope(0, iz) + Fxr(ifxr, iz)%area * Fxr(ifxr, iz)%burnup
          DO k = 1, Fxr(ifxr, iz)%niso
            iso = Fxr(ifxr, iz)%idiso(k);id = mapnucl(iso)
            Isotope(id, iz) = Isotope(id, iz) + Fxr(ifxr, iz)%area *   Fxr(ifxr,iz)%pnum(k)
          ENDDO
        ENDDO
        IF(volsum .GT. 0._8) THEN 
          volsum = 1._8 / volsum
          DO k = 0, ntiso
            Isotope(i, iz) = Isotope(i, iz) * volsum
          ENDDO
        ENDIF
      ENDDO
#ifdef MPI_ENV
      CALL REDUCE(Isotope(0:ntiso, 1:nz), Isotope_buf(0:ntiso, 1:nz), ntiso+1, nz, COMM, .FALSE.)
#endif
      IF(.NOT. MASTER) CYCLE
    ENDDO 
  ENDDO
ENDDO
IF(MASTER) CLOSE(io17)
DEALLOCATE(Isotope, Isotope_buf)
END SUBROUTINE



SUBROUTINE Normalization()
USE Depl_mod,         ONLY : FluxNormalizeFactor
IMPLICIT NONE
INTEGER :: ixy, iz, ig, ityp
REAL :: pwsum, volsum,AsyVolSum, pwsum_loc, volsum_loc, AsyVolSum_loc
REAL :: LocalPinPw

pwsum = 0; volsum = 0
DO iz = myzb, myze
  DO ixy = 1, nxy
    IF(.NOT. Pin(ixy)%lfuel) CYCLE
    LocalPinPw =  DotProduct(PinXs(ixy, iz)%XSKF(1:ng), PhiC(ixy, iz, 1:ng), ng)
    LocalPinPw = LocalPinPw * PinVol(ixy, iz)
    IF(LocalPinPw .LT. epsm10) CYCLE
    Pwsum = Pwsum + LocalPinPw
    volsum = volsum + PinVol(ixy, iz)
  ENDDO
ENDDO

hactive = 0
DO iz = 1, nz
  IF(Core%lFuelPlane(iz)) hactive = hactive + Core%hz(iz)  
ENDDO

AsyVolSum = 0
DO iz = 1, nz
  DO ixy = 1, Core%nxya
    ityp = Core%CoreMap(ixy)
    IF(Core%AsyInfo(ityp)%lFuel) THEN
      AsyVolSum = AsyVolSum + Core%AsyVol(ixy, iz)
    ENDIF
  ENDDO
ENDDO

#ifdef MPI_ENV
volsum_loc = volsum; pwsum_loc = pwsum
CALL REDUCE(volsum_loc, volsum,  PE%MPI_CMFD_COMM, .TRUE.)
CALL REDUCE(pwsum_loc, pwsum,  PE%MPI_CMFD_COMM, .TRUE.)
#endif
NormPin3D = volsum / pwsum
NormPin2D =  NormPin3D / hactive
NormAsy3D = AsyVolSum / pwsum
NormAsy2D = NormAsy3D / hactive
IF(.NOT. nTracerCntl%lXslib) NormAbsFlux = NormPin3D
IF(nTracerCntl%lXslib) THEN
  NormAbsFlux = FluxNormalizeFactor(Core, FmInfo, GroupInfo, DeplCntl%PowerCore, FALSE, TRUE, PE)
  NormAbsFlux = NormAbsFlux * nTracerCntl%PowerLevel
ENDIF


CONTINUE
END SUBROUTINE

SUBROUTINE SetAsyOutp()
IMPLICIT NONE
!Radial Information
CHARACTER(80) :: fn

INTEGER :: FuelPlane(nzmax)
INTEGER :: ixa, iya, ixya, ix, iy, ix0, iy0, ixy, iz, ityp
INTEGER :: nx, ny, nxbeg, nxend, nybeg, nyend


!Assembly Information
ALLOCATE(AsyOutp(Core%nxya))
DO ixya = 1, Core%nxya
  ityp = Core%Asy(ixya)%Asytype
  nx = Core%AsyInfo(ityp)%nx; ny = Core%AsyInfo(ityp)%ny
  IF(Core%lGap) THEN
    nx = nx - 2; ny = ny - 2
    IF(Core%AsyInfo(ityp)%lCentX .OR. Core%AsyInfo(ityp)%lCentXY) ny = ny + 1
    IF(Core%AsyInfo(ityp)%lCentY .OR. Core%AsyInfo(ityp)%lCentXY) nx = nx + 1
  ENDIF
  !WRITE(io17) ixy, Core%Asy(ixya)%ixa, Core%Asy(ixya)%iya, nx, ny
  
  AsyOutp(ixya)%ixa = Core%Asy(ixya)%ixa; AsyOutp(ixya)%iya = Core%Asy(ixya)%iya
  AsyOutp(ixya)%nx = nx; AsyOutp(ixya)%ny = ny; 
  AsyOutp(ixya)%nfuel = Core%AsyInfo(ityp)%nFuelPin
  
  AsyOutp(ixya)%itype = 0
  IF(Core%AsyInfo(ityp)%lfuel) AsyOutp(ixya)%itype = 1
  
  ALLOCATE(AsyOutp(ixya)%PinIdx(nx, ny))
  ALLOCATE(AsyOutp(ixya)%PinType(nx, ny))
  
  nxbeg = 1; nybeg  = 1
  nxend = Core%AsyInfo(ityp)%nx; nyend  = Core%AsyInfo(ityp)%ny
  IF(Core%lGap) THEN
    nxbeg = nxbeg + 1; nybeg = nybeg + 1
    nxend = nxend - 1; nyend = nyend - 1
    IF(Core%AsyInfo(ityp)%lCentX .OR. Core%AsyInfo(ityp)%lCentXY) nybeg = nybeg - 1
    IF(Core%AsyInfo(ityp)%lCentY .OR. Core%AsyInfo(ityp)%lCentXY) nxbeg = nxbeg - 1
  ENDIF
  
  iy = 0
  DO iy0 = nybeg, nyend
    iy = iy + 1; ix = 0
    DO ix0 = nxbeg, nxend
      ix = ix + 1
      ixy = Core%AsyInfo(ityp)%Pin2DIdx(ix0, iy0)
      ixy = Core%Asy(ixya)%GlobalPinIdx(ixy)
      AsyOutp(ixya)%PinIdx(ix, iy) = ixy
      AsyOutp(ixya)%PinType(ix, iy) = 0
      IF(Pin(ixy)%lFuel) AsyOutp(ixya)%PinType(ix, iy) = 1
    ENDDO
  ENDDO

ENDDO

END SUBROUTINE

SUBROUTINE SetBinOutpCntl()
!Binary Output file control
USE XSLIB_MOD,       ONLY : enbhel,     ldiso  
IMPLICIT NONE
REAL :: EStruct(ngmax),Eavg(ngmax), temp, ELB
INTEGER :: EIdx(ngmax)
INTEGER :: ig, ig0
INTEGER :: i, j, m

MASTER = PE%MASTER

ngc = nTRACERCntl%OutpCntl%FluxEdit
IF(ngc .EQ. 0) THEN
  ngc = ng
  DO ig = 1, ng
    Grpstruct(1:2, ig) = i
  ENDDO
ELSE
  EStruct(1:ng) = enbhel(1:ng)
  EStruct(ng+1) = 1.0E-4_8
  DO i = 1, ng
    Eavg(i) = (EStruct(i) + EStruct(i+1))/2._8
  ENDDO 
  DO i = 1, 3
    IF(nGC .EQ. nGList(i)) EXIT
  ENDDO
  m = i
  DO i = 1, nGC-1
    ELB = GrpLBList(i, m)
    DO j = 1, ng-1
      temp = (ELB-Eavg(j))*(ELB-Eavg(j+1))
      IF(temp .LE. 0._8) THEN
        EIdx(i)=j
        EXIT
      ENDIF
    ENDDO
    GrpStruct(1:2, 1) = (/1, EIdx(1)/)
    DO ig = 2, nGC-1
      GrpStruct(1:2, ig) = (/EIdx(ig-1)+1, Eidx(ig)/)
    ENDDO
    GrpStruct(1:2, nGC) = (/EIdx(nGC-1)+1, ng/)
  ENDDO
ENDIF

lBOUTP = nTracerCntl%OutpCntl%lBoutp
lFlux_BOutp = nTracerCntl%OutpCntl%lFlux_Boutp
lTH_BOutp = nTracerCntl%OutpCntl%lTH_Boutp
lDepl_Boutp = nTracerCntl%OutpCntl%lDepl_Boutp

COMM = PE%MPI_CMFD_COMM; NPROC = PE%nCMFDProc
nxy = Core%nxy; nz = Core%nz
myzb = PE%myzb; myze = PE%myze
IF(lDepl_Boutp) THEN
  ntiso = GroupInfo%ntiso
  IF(nTRACERCntl%OutpCntl%IsoBoutp(0) .EQ. 0) THEN
    nTRACERCntl%OutpCntl%IsoBoutp(0) = nMajorIsotope
    DO i = 1, nMajorIsotope
      nTRACERCntl%OutpCntl%IsoBoutp(i) = MajorIsotope(i)
    ENDDO
  ELSEIF(nTRACERCntl%OutpCntl%IsoBoutp(0) .EQ. -1) THEN
    nTRACERCntl%OutpCntl%IsoBoutp(0) = ntiso
    DO i = 1, ntiso
      nTRACERCntl%OutpCntl%IsoBoutp(i) = ldiso(i)%Nid
    ENDDO
  ENDIF
  nisoout = nTRACERCntl%OutpCntl%IsoBoutp(0)
  DO i = 1, nisoout
    IsoOutList(i) = nTRACERCntl%OutpCntl%IsoBoutp(i)
  ENDDO
ENDIF

Pin => Core%Pin
Fxr => FmInfo%Fxr; PinXs => CmInfo%PinXS
Phi => FmInfo%Phis; Psi => FmInfo%psi
PhiC => CmInfo%PhiC; PsiC => CmInfo%PsiC
PinVol => Core%PinVol; CellInfo => Core%CellInfo
END SUBROUTINE

SUBROUTINE GetFileName(fn, OutType, Id)
IMPLICIT NONE
INTEGER :: id
CHARACTER(*) :: OutType
CHARACTER(80) :: fn
CHARACTER(4) :: ext

ext='.nbo'
SELECT CASE(OutType)
  CASE('FLUX')
    IF(Id .LT. 10) THEN
      WRITE(fn, '(A, I1, A)') 'FLUX_00', id, ext
    ELSEIF(Id .LT. 100) THEN
      WRITE(fn, '(A, I2, A)') 'FLUX_0', id, ext 
    ELSE
      WRITE(fn, '(A, I2, A)') 'FLUX_', id, ext
    ENDIF
  CASE('DEPL')
    IF(Id .LT. 10) THEN
      WRITE(fn, '(A, I1, A)') 'DEPL_00', id, ext
    ELSEIF(Id .LT. 100) THEN
      WRITE(fn, '(A, I2, A)') 'DEPL_0', id, ext 
    ELSE
      WRITE(fn, '(A, I2, A)') 'DEPL_', id, ext
    ENDIF
  CASE('TH')
    IF(Id .LT. 10) THEN
      WRITE(fn, '(A, I1, A)') 'TH_00', id, ext
    ELSEIF(Id .LT. 100) THEN
      WRITE(fn, '(A, I2, A)') 'TH_0', id, ext 
    ELSE
      WRITE(fn, '(A, I2, A)') 'TH_', id, ext
    ENDIF
  CASE('GEOM')
    WRITE(fn, '(A,A)') 'GEOM', ext
END SELECT

END SUBROUTINE

END MODULE
