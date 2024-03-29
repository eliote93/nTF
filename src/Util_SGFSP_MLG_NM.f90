MODULE Util_SGFSP_MLG_NM

IMPLICIT NONE

CONTAINS

SUBROUTINE SetPlnLsigP_MLG_NM(Core, Fxr, Siglp, xst, iz, gb, ge)
USE TYPEDEF,        ONLY : CoreInfo_Type,   FxrInfo_Type,       Pin_Type,           Cell_Type
USE xslib_mod,      ONLY : libdata,         ldiso,              mapnucl,            mlgdata
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
REAL, POINTER :: Siglp(:, :), xst(:, :)
INTEGER :: iz, gb, ge

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
INTEGER :: ixy, nxy, nlv

nxy = Core%nxy
Pin => Core%Pin
CellInfo => Core%CellInfo
nlv = mlgdata(iz)%f_nmaclv

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED)
DO ixy = 1, nxy
 CALL SetPlnLsigP_MLG_NM_Pin(ixy)
ENDDO
!$OMP END DO
!$OMP END PARALLEL

CONTAINS

SUBROUTINE SetPlnLsigP_MLG_NM_Pin(ixy)

TYPE(libdata), POINTER :: isodata
TYPE(FxrInfo_Type), POINTER :: myFxr
INTEGER :: FxrIdxSt, FsrIdxSt, nLocalFxr, nFsrInFxr
INTEGER :: i, j, iso, id, idx, ilv, ig, ixy, ifxr, ifsr, icel
INTEGER :: niso
INTEGER, POINTER :: idiso(:)
REAL :: lv, LocalSigtr
REAL, POINTER :: pnum(:)

icel = Pin(ixy)%Cell(iz)

FxrIdxSt = Pin(ixy)%FxrIdxSt
FsrIdxSt = Pin(ixy)%FsrIdxSt
nLocalFxr = CellInfo(icel)%nFxr
DO j = 1, nLocalFxr
  ifxr = FxrIdxSt + j - 1
  myFxr => Fxr(ifxr, iz)
  niso = myFxr%niso; pnum => myFxr%pnum; idiso => myFxr%idiso
  DO ig = gb, ge
    DO ilv = 1, nlv
      idx = ilv + (ig - gb) * nlv
      lv = mlgdata(iz)%f_maclv(ilv)
      LocalSigtr = 0.0
      DO iso = 1, niso
        id = mapnucl(idiso(iso)); isodata => ldiso(id)
        IF (isodata%sigp .EQ. 0.0) THEN
          LocalSigtr = LocalSigtr + pnum(iso) * isodata%lamsigp(ig)
        ELSE
          LocalSigtr = LocalSigtr + pnum(iso) * isodata%sigp
        ENDIF
      ENDDO
      Siglp(idx, ifxr) = LocalSigtr
      IF (Fxr(ifxr, iz)%lres) THEN
        IF ((.NOT. Fxr(ifxr, iz)%lCLD) .AND. (.NOT. Fxr(ifxr, iz)%lAIC)) THEN
          LocalSigtr = LocalSigtr + lv * myFxr%FnAdj(ig) * myFxr%FtAdj(ilv, ig)
        ENDIF
      ENDIF
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
      DO i = 1, nFsrInFxr
        ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
        xst(idx, ifsr) = LocalSigtr
      ENDDO
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE

END SUBROUTINE

SUBROUTINE SetPlnLsigP_1gMLG_NM(Core, Fxr, Siglp, xst, iz, lCLD, lAIC)
USE TYPEDEF,        ONLY : CoreInfo_Type,   FxrInfo_Type,       Pin_Type,           Cell_Type
USE xslib_mod,      ONLY : libdata,         ldiso,              mapnucl,            mlgdata,                        &
                           mlgdata0
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
REAL, POINTER :: Siglp(:, :), xst(:, :)
INTEGER :: iz
LOGICAL :: lCLD, lAIC

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
INTEGER :: ixy, nxy, nlv

nxy = Core%nxy
Pin => Core%Pin
CellInfo => Core%CellInfo
IF (lCLD) nlv = mlgdata0%c_nmaclv1G
IF (lAIC) nlv = mlgdata(iz)%f_nmaclv1G

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED)
DO ixy = 1, nxy
  CALL SetPlnLsigP_1gMLG_NM_Pin(ixy)
ENDDO
!$OMP END DO
!$OMP END PARALLEL

CONTAINS

SUBROUTINE SetPlnLsigP_1gMLG_NM_Pin(ixy)

TYPE(FxrInfo_Type), POINTER :: myFxr
TYPE(libdata), POINTER :: isodata
INTEGER :: FxrIdxSt, FsrIdxSt, nLocalFxr, nFsrInFxr
INTEGER :: i, j, ilv, ixy, ifxr, ifsr, icel, iso, id
INTEGER :: niso
INTEGER, POINTER :: idiso(:)
REAL :: Localsiglp, LocalSigtr
REAL, POINTER :: pnum(:)

icel = Pin(ixy)%Cell(iz)

FxrIdxSt = Pin(ixy)%FxrIdxSt
FsrIdxSt = Pin(ixy)%FsrIdxSt
nLocalFxr = CellInfo(icel)%nFxr
DO j = 1, nLocalFxr
  ifxr = FxrIdxSt + j - 1
  myFxr => Fxr(ifxr, iz)
  niso = myFxr%niso; pnum => myFxr%pnum; idiso => myFxr%idiso
  DO ilv = 1, nlv
    LocalSiglp = 0.0
    DO iso = 1, niso
      id = mapnucl(idiso(iso)); isodata => ldiso(id)
      LocalSiglp = LocalSiglp + pnum(iso) * isodata%sigp
    ENDDO
    LocalSigtr = LocalSiglp
    IF (lCLD) THEN
      IF (Fxr(ifxr, iz)%lres .AND. Fxr(ifxr, iz)%lCLD) LocalSigtr = LocalSigtr + mlgdata0%c_maclv1G(ilv)
    ELSEIF (lAIC) THEN
      IF (Fxr(ifxr, iz)%lres .AND. Fxr(ifxr, iz)%lAIC) LocalSigtr = LocalSigtr + mlgdata(iz)%f_maclv1G(ilv)
    ENDIF
    Siglp(ilv, ifxr) = LocalSiglp
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      xst(ilv, ifsr) = LocalSigtr
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE

END SUBROUTINE

SUBROUTINE SetSubGrpSrc_NM(Core, Fxr, Siglp, xst, src, iz, gb, ge)
USE TYPEDEF,        ONLY : CoreInfo_Type,   FxrInfo_Type,       Pin_Type,           Cell_Type
USE xslib_mod,      ONLY : mlgdata
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER  :: Fxr(:, :)
REAL, POINTER :: Siglp(:, :), xst(:, :), src(:, :)
INTEGER :: iz, gb, ge

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
INTEGER :: ixy, nxy

nxy = Core%nxy
Pin => Core%Pin
CellInfo => Core%CellInfo

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED)
DO ixy = 1, nxy
  CALL SetSubGrpSrc_NM_Pin(ixy)
ENDDO
!$OMP END DO
!$OMP END PARALLEL

CONTAINS

SUBROUTINE SetSubGrpSrc_NM_Pin(ixy)

INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nFsrInFxr
INTEGER :: i, j, ig, icel, ixy, ifxr, ifsr
REAL :: xstrinv, localsrc

icel = Pin(ixy)%Cell(iz)

FsrIdxSt = Pin(ixy)%FsrIdxSt
FxrIdxSt = Pin(ixy)%FxrIdxSt
nLocalFxr = CellInfo(icel)%nFxr
DO j = 1, nLocalFxr
  ifxr = FxrIdxSt + j - 1
  nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
  DO i = 1, nFsrInFxr
    ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
    DO ig = gb, ge
      src(ig, ifsr) = Siglp(ig, ifxr) / xst(ig, ifsr)
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE

END SUBROUTINE

SUBROUTINE EquipXSGen_MLG_NM(Core, Fxr, Siglp, xst, phis, iz, gb, ge)
USE TYPEDEF,        ONLY : CoreInfo_Type,   FxrInfo_Type,       ResVarPin_Type,     Pin_Type,                       &
                           Cell_Type
USE xslib_mod,      ONLY : mlgdata
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
REAL, POINTER :: Siglp(:, :), xst(:, :), phis(:, :)
INTEGER :: iz, gb, ge

TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(ResVarPin_Type), POINTER :: ResVarPin(:, :)
INTEGER :: ixy, nlv, nxy

nxy = Core%nxy
Pin => Core%Pin
ResVarPin => Core%ResVarPin
CellInfo => Core%CellInfo
nlv = mlgdata(iz)%f_nmaclv

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED)
DO ixy = 1, nxy
  CALL EquipXSGen_MLG_NM_Pin(ixy)
ENDDO
!$OMP END DO
!$OMP END PARALLEL

CONTAINS

SUBROUTINE EquipXSGen_MLG_NM_Pin(ixy)

TYPE(FxrInfo_Type), POINTER :: myFxr
INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nFsrInFxr
INTEGER :: idx, ilv, ig, ixy, icel, ifxr, ifsr, i, j, l
REAL :: vol, volsum, phisum, xseq, maclp, maclv, avgphisum, avgvolsum, maclpavg, maclvavg

icel = Pin(ixy)%Cell(iz)
IF (.NOT. CellInfo(icel)%lres) RETURN

FsrIdxSt = Pin(ixy)%FsrIdxSt
FxrIdxSt = Pin(ixy)%FxrIdxSt
nLocalFxr = CellInfo(icel)%nFxr
DO ig = gb, ge
  DO ilv = 1, nlv
    idx = ilv + (ig - gb) * nlv
    avgphisum = 0.0; avgvolsum = 0.0
    maclvavg = 0.0; maclpavg = 0.0
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j - 1; nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
      myFxr => Fxr(ifxr, iz); ifsr = myFxr%FsrIdxSt
      IF (.NOT. myFxr%lres) CYCLE
      ifsr = myFxr%FsrIdxSt
      maclp = Siglp(idx, ifxr)
      maclv = xst(idx, ifsr) - maclp
      IF (maclv .EQ. 0.0) CYCLE
      phisum = 0.0; volsum = 0.0
      DO i = 1, nFsrInFxr
        l = Cellinfo(icel)%MapFxr2FsrIdx(i, j); vol = CellInfo(icel)%vol(l)
        ifsr = FsrIdxSt + l - 1
        phisum = phisum + phis(idx, ifsr) * vol; volsum = volsum + vol
      ENDDO
      maclpavg = maclpavg + maclp * phisum
      maclvavg = maclvavg + maclv * phisum
      avgphisum = avgphisum + phisum
      avgvolsum = avgvolsum + volsum
      phisum = phisum / volsum
      IF (abs(phisum - 1.0) .LT. 1.0E-10) THEN
        xseq = 1.0E+10
      ELSE
        xseq = - maclp + maclv * phisum / (1.0 - phisum)
      ENDIF
!      xseq = - maclp + maclv * phisum / (1.0 - phisum)
      myFxr%xseq_f_mg(ilv, ig) = xseq
    ENDDO
    IF (maclpavg .EQ. 0.0) CYCLE
    maclpavg = maclpavg / avgphisum
    maclvavg = maclvavg / avgphisum
    avgphisum = avgphisum / avgvolsum
    IF (abs(avgphisum - 1.0) .LT. 1.0E-10) THEN
      xseq = 1.0E+10
    ELSE
      xseq = - maclpavg + maclvavg * avgphisum / (1.0 - avgphisum)
    ENDIF
!    xseq = - maclpavg + maclvavg * avgphisum / (1.0 - avgphisum)
    ResVarPin(ixy, iz)%avgxseq_mg(ilv, ig) = xseq
  ENDDO
ENDDO

END SUBROUTINE

END SUBROUTINE

SUBROUTINE EquipXSGen_1gMLG_NM(Core, Fxr, Siglp, phis, xst, iz, nlv, lCLD, lAIC)
USE TYPEDEF,        ONLY : CoreInfo_Type,   FxrInfo_Type,       ResVarPin_Type,     Pin_Type,                       &
                           Cell_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
REAL, POINTER :: Siglp(:, :), phis(:, :), xst(:, :)
INTEGER :: iz, nlv
LOGICAL :: lCLD, lAIC

TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(ResVarPin_Type), POINTER :: ResVarPin(:, :)
INTEGER :: ixy, nxy

nxy = Core%nxy
Pin => Core%Pin
ResVarPin => Core%ResVarPin
CellInfo => Core%CellInfo

IF (lCLD) THEN
  !$OMP PARALLEL
  !$OMP DO SCHEDULE(GUIDED)
  DO ixy = 1, nxy
    CALL EquipXSGen_1gMLG_NM_Pin_CLD(ixy)
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDIF

IF (lAIC) THEN
  !$OMP PARALLEL
  !$OMP DO SCHEDULE(GUIDED)
  DO ixy = 1, nxy
    CALL EquipXSGen_1gMLG_NM_Pin_AIC(ixy)
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDIF

CONTAINS

SUBROUTINE EquipXSGen_1gMLG_NM_Pin_CLD(ixy)

TYPE(FxrInfo_Type), POINTER :: myFxr
INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nFsrInFxr
INTEGER :: ilv, ixy, icel, ifxr, ifsr, i, j, l
REAL :: vol, volsum, phisum, xseq, maclp, maclv

icel = Pin(ixy)%Cell(iz)
IF (.NOT. CellInfo(icel)%lres) RETURN

FsrIdxSt = Pin(ixy)%FsrIdxSt
FxrIdxSt = Pin(ixy)%FxrIdxSt
nLocalFxr = CellInfo(icel)%nFxr
DO j = 1, nLocalFxr
  ifxr = FxrIdxSt + j - 1; nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
  myFxr => Fxr(ifxr, iz)
  IF (.NOT. myFxr%lres) CYCLE
  DO ilv = 1, nlv
    ifsr = myFxr%FsrIdxSt
    maclp = Siglp(ilv, ifxr)
    maclv = xst(ilv, ifsr) - maclp
    IF (maclv .EQ. 0.0) CYCLE
    phisum = 0.0; volsum = 0.0
    DO i = 1, nFsrInFxr
      l = Cellinfo(icel)%MapFxr2FsrIdx(i, j); vol = CellInfo(icel)%vol(l)
      ifsr = FsrIdxSt + l - 1
      phisum = phisum + phis(ilv, ifsr) * vol; volsum = volsum + vol
    ENDDO
    phisum = phisum / volsum
    IF (abs(phisum - 1.0) .LT. 1.0E-10) THEN
      xseq = 1.0E+10
    ELSE
      xseq = - maclp + maclv * phisum / (1.0 - phisum)
    ENDIF
    myFxr%xseq_c_1g(ilv) = xseq
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE EquipXSGen_1gMLG_NM_Pin_AIC(ixy)

TYPE(FxrInfo_Type), POINTER :: myFxr
INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nFsrInFxr
INTEGER :: ilv, ixy, ifxr, icel, ifsr, nxy, i, j, l
REAL :: vol, volsum, phisum, xseq, maclp, maclv, avgphisum, avgvolsum, maclpavg, maclvavg

icel = Pin(ixy)%Cell(iz)
IF (.NOT. CellInfo(icel)%lres) RETURN
IF (.NOT. CellInfo(icel)%lAIC) RETURN

FsrIdxSt = Pin(ixy)%FsrIdxSt
FxrIdxSt = Pin(ixy)%FxrIdxSt
nLocalFxr = CellInfo(icel)%nFxr
DO ilv = 1, nlv
  maclpavg = 0.0; maclvavg = 0.0; avgphisum = 0.0; avgvolsum = 0.0
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1; nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    myFxr => Fxr(ifxr, iz)
    IF (.NOT. myFxr%lres) CYCLE
    ifsr = myFxr%FsrIdxSt
    maclp = Siglp(ilv, ifxr)
    maclv = xst(ilv, ifsr) - maclp
    IF (maclv .EQ. 0.0) CYCLE
    phisum = 0.0; volsum = 0.0
    DO i = 1, nFsrInFxr
      l = Cellinfo(icel)%MapFxr2FsrIdx(i, j); vol = CellInfo(icel)%vol(l)
      ifsr = FsrIdxSt + l - 1
      phisum = phisum + phis(ilv, ifsr) * vol; volsum = volsum + vol
    ENDDO
    maclpavg = maclpavg + maclp * phisum
    maclvavg = maclvavg + maclv * phisum
    avgphisum = avgphisum + phisum
    avgvolsum = avgvolsum + volsum
    phisum = phisum / volsum
    IF (abs(phisum - 1.0) .LT. 1.0E-10) THEN
      xseq = 1.0E+10
    ELSE
      xseq = - maclp + maclv * phisum / (1.0 - phisum)
    ENDIF
    myFxr%xseq_f_1g(ilv) = xseq
  ENDDO
  IF (maclpavg .EQ. 0.0) CYCLE
  maclpavg = maclpavg / avgphisum
  maclvavg = maclvavg / avgphisum
  avgphisum = avgphisum / avgvolsum
  IF (abs(avgphisum - 1.0) .LT. 1.0E-10) THEN
    xseq = 1.0E+10
  ELSE
    xseq = - maclpavg + maclvavg * avgphisum / (1.0 - avgphisum)
  ENDIF
  ResVarPin(ixy, iz)%avgxseq_1g(ilv) = xseq
ENDDO

END SUBROUTINE

END SUBROUTINE

SUBROUTINE UpdtFnAdj_NM(Core, Fxr, iz, gb, ge)
USE TYPEDEF,        ONLY : CoreInfo_Type,   FxrInfo_Type,       Cell_Type,          ResVarPin_Type,                 &
                           Pin_Type
USE xslib_mod,      ONLY : libdata,         ldiso,              mapnucl
USE XsUtil_mod,     ONLY : LineIntpol
USE TH_Mod,         ONLY : GetPinFuelTemp
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
INTEGER :: iz, gb, ge

TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(ResVarPin_Type), POINTER :: ResVarPin(:, :)
REAL, POINTER :: NresPinAvg(:, :)
INTEGER :: ig, nxy

nxy = Core%nxy
CellInfo => Core%CellInfo
Pin => Core%Pin
ResVarPin => Core%ResVarPin

ALLOCATE(NresPinAvg(nxy, gb : ge)); NresPinAvg = 0.0

!$OMP PARALLEL
!$OMP DO SCHEDULE(DYNAMIC)
DO ig = gb, ge
  CALL UpdtFnAdj_NM_Group(ig)
ENDDO
!$OMP END DO
!$OMP END PARALLEL

DEALLOCATE(NresPinAvg)

CONTAINS

SUBROUTINE UpdtFnAdj_NM_Group(ig)

TYPE(FxrInfo_Type), POINTER :: myFxr
TYPE(libdata), POINTER :: isodata
INTEGER :: FxrIdxSt, nLocalFxr
INTEGER :: i, j, jd, ig, ipin, icel, ifxr, iso
INTEGER :: niso, nRiTemp, npot
INTEGER, POINTER :: idiso(:)
REAL :: TempAvgsq, N, Nsum, Npin, areasum, pinareasum, ria
REAL, POINTER :: pnum(:)

Nsum = 0.0; areasum = 0.0

DO ipin = 1, nxy
  FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz); nLocalFxr = CellInfo(icel)%nFxr
  IF (.NOT. ResVarPin(ipin, iz)%lres) CYCLE
  IF (.NOT. CellInfo(icel)%lres) CYCLE
  IF (CellInfo(icel)%lAIC) CYCLE
  TempAvgsq = dsqrt(GetPinFuelTemp(Core, Fxr, iz, ipin))
  Npin = 0.0; pinareasum = 0.0
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1
    myFxr => Fxr(ifxr, iz)
    IF (.NOT. myFxr%lres) CYCLE
    IF (myFxr%lCLD) CYCLE
    niso = myFxr%niso; idiso => myFxr%idiso; pnum => myFxr%pnum
    N = 0.0
    DO iso = 1, niso
      jd = mapnucl(idiso(iso)); isodata => ldiso(jd)
      IF (.NOT. isodata%lreso) CYCLE
      npot = isodata%nsig0; nRiTemp = isodata%nrtemp
      ria = LineIntpol(TempAvgsq, nRiTemp, isodata%rtempsq(1 : nRiTemp), isodata%ri_a(npot, ig, 1 : nRiTemp))
      N = N + pnum(iso) * ria
    ENDDO
    myFxr%FnAdj(ig) = N
    Npin = Npin + N * myFxr%area
    pinareasum = pinareasum + myFxr%area
  ENDDO
  IF (pinareasum .GT. 0.0) NresPinAvg(ipin, ig) = Npin / pinareasum
  areasum = areasum + pinareasum
  Nsum = Nsum + Npin
ENDDO

Nsum = Nsum / areasum

DO ipin = 1, nxy
  FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz); nLocalFxr = CellInfo(icel)%nFxr
  IF (.NOT. ResVarPin(ipin, iz)%lres) CYCLE
  IF (.NOT. CellInfo(icel)%lres) CYCLE
  IF (CellInfo(icel)%lAIC) CYCLE
  ResVarPin(ipin, iz)%FnAdj(ig) = NresPinAvg(ipin, ig) / Nsum
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1
    myFxr => Fxr(ifxr, iz)
    IF (.NOT. myFxr%lres) CYCLE
    IF (myFxr%lCLD) CYCLE
    myFxr%FnAdj(ig) = myFxr%FnAdj(ig) / Nsum
  ENDDO
ENDDO

END SUBROUTINE

END SUBROUTINE

SUBROUTINE UpdtFtAdj_NM(Core, Fxr, iz, gb, ge)
USE Typedef,        ONLY : CoreInfo_Type,   FxrInfo_Type,       Cell_Type,          Pin_Type
USE xslib_mod,      ONLY : libdata,         ldiso,              mapnucl,            mapnuclres,                     &
                           mlgdata
USE XsUtil_mod,     ONLY : LineIntPol,      LineIntPol2
USE TH_Mod,         ONLY : GetPinFuelTemp
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
INTEGER :: iz, gb, ge

TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Pin_Type), POINTER :: Pin(:)
INTEGER :: ixy, nxy, nlv

nxy = Core%nxy
Pin => Core%Pin
CellInfo => Core%CellInfo
nlv = mlgdata(iz)%f_nmaclv

!$OMP PARALLEL
!$OMP DO SCHEDULE(GUIDED)
DO ixy = 1, nxy
  CALL UpdtFtAdj_NM_Pin(ixy)
ENDDO
!$OMP END DO
!$OMP END PARALLEL

CONTAINS

SUBROUTINE UpdtFtAdj_NM_Pin(ixy)

TYPE(FxrInfo_Type), POINTER :: myFxr
TYPE(libdata), POINTER :: isodata, jsodata

INTEGER :: FxrIdxSt, nLocalFxr
INTEGER :: i, j, ig, ilv, ifxr, icel, ixy, idres, jdres, id, jd, iso, jso
INTEGER :: niso, npot, nRiTemp, n1, n2, n3, n4, lb, ub
INTEGER, POINTER :: idiso(:)
REAL :: ind, siglp, TempAvgsq, Tempsq, micsigb, micsig0, sigbsq, xdat(20), ydat(20), Nreg, Navg, I_avg, I_reg
REAL, POINTER :: pnum(:)

icel = Pin(ixy)%Cell(iz)
IF (.NOT. CellInfo(icel)%lres) RETURN
IF (CellInfo(icel)%lAIC) RETURN

FxrIdxSt = Pin(ixy)%FxrIdxSt
nLocalFxr = CellInfo(icel)%nFxr
TempAvgsq = dsqrt(GetPinFuelTemp(Core, Fxr, iz, ixy))
DO j = 1, nLocalFxr
  ifxr = FxrIdxSt + j - 1
  myFxr => Fxr(ifxr, iz)
  niso = myFxr%niso; idiso => myFxr%idiso; pnum => myFxr%pnum
  Tempsq = dsqrt(myFxr%temp)
  IF (.NOT. myFxr%lres) CYCLE
  IF (myFxr%lCLD) CYCLE
  DO ig = gb, ge
    DO ilv = 1, nlv
      IF (myFxr%xseq_f_mg(ilv, ig) .EQ. 0.0) CYCLE
      siglp = 0.0
      DO jso = 1, niso
        jd = mapnucl(idiso(jso)); jsodata => ldiso(jd)
        siglp = siglp + pnum(jso) * jsodata%lamsigp(ig)
      ENDDO
      Nreg = 0.0; Navg = 0.0
      DO iso = 1, niso
        id = mapnucl(idiso(iso)); isodata => ldiso(id)
        IF (.NOT. isodata%lreso) CYCLE
        idres = mapnuclres(idiso(iso), iz)
        ind = 0.0
        DO jso = 1, niso
          jd = mapnucl(idiso(jso)); jsodata => ldiso(jd)
          jdres = mapnuclres(idiso(jso), iz)
          IF (jdres .EQ. idres) ind = ind + pnum(jso)
        ENDDO
        micsigb = (myFxr%XsEq_f_mg(ilv, ig) + siglp) / ind / myFxr%FtAdj(ilv, ig)
        micsig0 = micsigb - isodata%lamsigp(ig)
        IF (micsig0 .LE. 0.0) micsig0 = 1.0E-10
        sigbsq = dsqrt(micsig0)
        npot = isodata%nsig0
        nRiTemp = isodata%nrtemp
        n1 = 1; n3 = 1;
        DO i = 1, nRiTemp
          n2 = i;
          IF (isodata%rtempsq(i)>tempsq) EXIT
          n1 = i
          !xdat(i) = isodata%rtempsq(i)
          !ydat(i) = LineIntPol2(sigbsq, npot, isodata%sig0sq(1 : npot), isodata%ri_a(1 : npot, ig, i))
        ENDDO
        DO i = 1, nRiTemp
          n4 = i;
          IF (isodata%rtempsq(i)>TempAvgsq) EXIT
          n3 = i
        ENDDO
        lb = min(n1,n3); ub = max(n2,n4)
        DO i = n1,n2
          xdat(i) = isodata%rtempsq(i)
          ydat(i) = LineIntPol2(sigbsq, npot, isodata%sig0sq(1:npot),isodata%ri_a(1:npot,ig,i))
        END DO
        IF (n1.EQ.n3) lb = n4
        IF (n2.EQ.n4) ub = n3
        DO i = lb, ub
          xdat(i) = isodata%rtempsq(i)
          ydat(i) = LineIntPol2(sigbsq, npot, isodata%sig0sq(1:npot),isodata%ri_a(1:npot,ig,i))
        END DO
        I_reg = LineIntpol(tempsq, (n2-n1)+1, xdat(n1:n2), ydat(n1:n2))
        Nreg = Nreg + pnum(iso) * I_reg
        I_avg = LineIntpol(TempAvgsq, (n4-n3)+1, xdat(n3:n4), ydat(n3:n4))
        Navg = Navg + pnum(iso) * I_avg
      ENDDO
      myFxr%FtAdj(ilv, ig) = Nreg / Navg
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE

END SUBROUTINE

END MODULE