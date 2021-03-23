#include <defines.h>
MODULE CritSpec_mod
USE PARAM
USE TYPEDEF,        ONLY : XsMac_Type
USE MOC_MOD,        ONLY : FxrAvgPhi
USE MacXsLib_mod,   ONLY : MacXsBase,  MacXsScatMatrix,  MacP1XsScatMatrix, MacXsANf
USE BenchXs,        ONLY : XsBaseBen,  XsBaseDynBen
USE BasicOperation, ONLY : CP_CA, CP_VA, MULTI_VA, DotProduct, DIV_VA,        &
                           MULTI_VA
USE XsUtil_mod,     ONLY : AllocXsMac
IMPLICIT NONE
TYPE(XsMac_Type), PRIVATE :: XsMacBase(128)
INTEGER, PARAMETER, PRIVATE :: MaxGrp = 1000
CONTAINS

SUBROUTINE GetCriticalSpectrum(Core, FmInfo, GroupInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,      ONLY: CoreInfo_Type,     FmInfo_Type,     GroupInfo_Type,    &
                        PE_TYPE
USE B1_Mod,       ONLY : B1Calculation
USE CNTL,         ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
TYPE(XsMac_Type), SAVE :: CoreXsMac
REAL :: CoreSpec(MaxGrp), CritSpec(MaxGrp)
REAL :: Bsq, Fnorm
INTEGER :: ng
INTEGER :: ig
ng = GroupInfo%ng

CALL CoreAvgXs(CoreXsMac, CoreSpec, Core, FmInfo, GroupInfo, nTracerCntl, PE)
CALL B1Calculation(CoreXsMac, CritSpec(1:ng), Bsq, ng, .FALSE.)

Fnorm = 1._8 / SUM(CritSpec(1:ng))
CritSpec(1:ng) = Fnorm * CritSpec

CALL CP_VA(FmInfo%PhiCrit(1:ng), CritSpec(1:ng), ng)
!Spectrum Conversion Factor
DO ig = 1, ng
  FmInfo%SpecConv(ig) = CritSpec(ig) / CoreSpec(ig)
ENDDO
CONTINUE
END SUBROUTINE


!--- BYSedit
SUBROUTINE GetDiffusionCoeff(Core, FmInfo, GroupInfo, nTracerCntl, PE, Dng, phicrit, Bsq, kinf_Msq, kinf, keff)
USE PARAM
USE TYPEDEF,      ONLY: CoreInfo_Type,     FmInfo_Type,     GroupInfo_Type,    &
                        PE_TYPE
USE B1_Mod,       ONLY : B1Calculation_D
USE CNTL,         ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE
TYPE(XsMac_Type), SAVE :: CoreXsMac
REAL :: CoreSpec(MaxGrp), CritSpec(MaxGrp)
REAL :: Bsq, Fnorm,kinf, keff, kinf_Msq
INTEGER :: ng
INTEGER :: ig
REAL :: phicrit(GroupInfo%ng), Dng(GroupInfo%ng)
ng = GroupInfo%ng

CALL CoreAvgXs(CoreXsMac, CoreSpec, Core, FmInfo, GroupInfo, nTracerCntl, PE)
CALL B1Calculation_D(CoreXsMac, CritSpec(1:ng), Bsq, ng, .FALSE.,Dng,kinf, keff, kinf_Msq)
phicrit(1:ng)=CritSpec(1:ng)
CONTINUE
END SUBROUTINE
!--- BYSedit end

SUBROUTINE CoreAvgXS(XsMac, Spec, Core, FmInfo, GroupInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,       ONLY : CoreInfo_Type,                                       &
                          FmInfo_Type,        GroupInfo_Type,     PE_TYPE,     &
                          FxrInfo_Type,       Cell_Type,          Pin_Type
USE CNTL,          ONLY : nTracerCntl_Type

IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
REAL :: SPEC(MaxGrp)
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE

TYPE(FxrInfo_Type), POINTER :: Fxr(:, :), myFxr
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)

REAL, POINTER :: Phis(:, :, :)
REAL, POINTER :: PinVol(:, :)
REAL, POINTER :: Hz(:)

REAL, POINTER :: XsMacA(:), XsMacNf(:), XsMacT(:), XsMacTr(:), chi(:)
REAL, POINTER :: XsMacSm(:, :), XsMacP1Sm(:, :)
REAL :: PhiFxr(MaxGrp), PhiSum(MaxGrp)
REAL :: Area, Fis, FisSum

REAL :: SigSum, SigSumP1

INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nLocalFsr, nFsrInFxr
INTEGER :: nFsr, nFxr, nxy, nz, myzb, myze, norg, nchi
INTEGER :: ng, iResoGrpBeg, iResoGrpEnd
INTEGER :: iz, ipin, icel, ifxr, ifsr, ifsrlocal
INTEGER :: i, j, ig, ig2
LOGICAL :: lXsLib

lXsLib = nTracerCntl%lXsLib
IF(.NOT. lXsLib) RETURN

Pin => Core%Pin; CellInfo => Core%CellInfo
PinVol => Core%PinVol; Hz => Core%Hz

Fxr => FmInfo%Fxr; Phis => FmInfo%Phis

ng = GroupInfo%ng
iresoGrpBeg = GroupInfo%nofg + 1; iresoGrpEnd = GroupInfo%nofg + GroupInfo%norg
norg = GroupInfo%norg;  nchi = GroupInfo%nchi

myzb = PE%myzb; myze = PE%myze
nFsr = Core%nCoreFsr;   nFxr = Core%nCoreFxr
nxy = Core%nxy

IF(.NOT. XsMac%lalloc) THEN
  XsMac%ng = ng
  CALL AllocXsMac(XsMac)
ENDIF

XsMacA => XsMac%XsMacA; XsMacNf => XsMac%XsMacNf
XsMacT => XsMac%XsMacT; XsMacTr => XsMac%XsMacTr
XsMacSm => XsMac%XsMacSm; XsMacP1Sm => XsMac%XsMacP1Sm
Chi => XsMac%Chi

CALL CP_CA(XsMacA(1:ng), 0._8, ng); CALL CP_CA(XsMacNf(1:ng), 0._8, ng)
CALL CP_CA(XsMacSm(1:ng, 1:ng), 0._8, ng, ng); CALL CP_CA(XsMacP1Sm(1:ng, 1:ng), 0._8, ng, ng)
CALL CP_CA(CHI(1:ng), 0._8, ng)
FisSum = 0; phisum = 0
!CALL

DO iz = myzb, myze
  DO ipin = 1, nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
    icel = Pin(ipin)%Cell(iz)
    nlocalFxr = CellInfo(icel)%nFxr
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      myFxr => Fxr(ifxr, iz)
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
      PhiFxr(1:ng) = 0
      DO i = 1, nFsrInFxr
        ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(i, j)
        iFsr = FsrIdxSt + ifsrlocal - 1
        Area = CellInfo(icel)%vol(ifsrlocal)
        PhiFxr(1:ng) = PhiFxr(1:ng) + Area * Phis(iFsr, iz, 1:ng)
      ENDDO
      PhiFxr(1:ng) = PhiFxr(1:ng) * Hz(iz)

      CALL MacXsBase(XsMacBase(1), myFxr, 1, ng, ng, 1._8, FALSE, TRUE)
      CALL MacXsScatMatrix(XsMacBase(1), myFxr, 1, ng, ng, GroupInfo, TRUE)
      IF(myFxr%lres) THEN
        DO ig = iResoGrpBeg,iResoGrpEnd
          XsMacBase(1)%XsMacA(ig) = XsMacBase(1)%XsMacA(ig) * myFxr%FresoA(ig)
          XsMacBase(1)%XsMacNf(ig) = XsMacBase(1)%XsMacNf(ig) * myFxr%FresoNF(ig)
        ENDDO
      ENDIF
      CALL MacP1XsScatMatrix(XsMacBase(1), myFxr, 1, ng, ng, GroupInfo)

      !Added Reaction values
      XsMacA(1:ng) = XsMacA(1:ng) + PhiFxr(1:ng) * XsMacBase(1)%XsMacA(1:ng)
      XsMacNf(1:ng) = XsMacNf(1:ng) + PhiFxr(1:ng) * XsMacBase(1)%XsMacNf(1:ng)
      DO ig = 1, ng
        DO ig2 = 1, ng
          XsMacSm(ig2, ig) = XsMacSm(ig2, ig) + XsMacBase(1)%XsMacSm(ig2, ig) * PhiFxr(ig2)
          XsMacP1Sm(ig2, ig) = XsMacP1Sm(ig2, ig) + XsMacBase(1)%XsMacP1Sm(ig2, ig) * PhiFxr(ig2)
        ENDDO
      ENDDO
      IF(myFxr%ldepl) THEN
        Fis = DotProduct(PhiFxr(1:ng), XsMacBase(1)%XsMacNf(1:ng), ng)
        DO ig = 1, nchi
          CHI(ig) = CHI(ig) + Fis * myFxr%Chi(ig)
        ENDDO
        FisSum = FisSum + Fis
      ENDIF
      PhiSum(1:ng) = PhiSum(1:ng) + PhiFxr(1:ng)
    ENDDO
  ENDDO
ENDDO

CALL DIV_VA(XsMacA(1:ng), PhiSum(1:ng), ng)
CALL DIV_VA(XsMacNf(1:ng), PhiSum(1:ng), ng)

DO ig = 1, ng
  DO ig2 = 1, ng
    XsMacSm(ig2, ig) = XsMacSm(ig2, ig) / PhiSum(ig2)
    XsMacP1Sm(ig2, ig) = XsMacP1Sm(ig2, ig) / PhiSum(ig2)
  ENDDO
ENDDO

DO ig = 1, ng
  CHI(ig) = CHI(ig) / FisSum
ENDDO

DO ig = 1, ng
  SigSum = 0; SigSumP1 = 0
  DO ig2 = 1, ng
    SigSum = SigSum + XsMacSm(ig, ig2)
    SigSumP1 = SigSumP1 + XsMacP1Sm(ig, ig2)
  ENDDO
  XsMacT(ig) = XsMacA(ig) + SigSum
  XsMacTr(ig) = XsMacA(ig) + SigSumP1
ENDDO

!Flux Normalization
FisSum = 1._8 / Sum(PhiSum(1:ng))
Spec(1:ng) = FisSum * PhiSum(1:ng)
NULLIFY(XsMacA, XsMacNf, XsMacT, XsMacTr, CHI)
NULLIFY(XsMacSm, XsMacP1Sm)

NULLIFY(myFxr)
NULLIFY(Pin, CellInfo, Pinvol, Hz)
NULLIFY(Fxr, Phis)
END SUBROUTINE



FUNCTION XsKinf(Core, FmInfo, GroupInfo, lXsLib, lCritSpec, PE)
USE TYPEDEF,      ONLY : CoreInfo_Type,        FmInfo_TYpe,       PE_TYPE,         &
                         GroupInfo_Type,                                           &
                         FxrInfo_Type,         Cell_Type,         Pin_Type
USE OMP_LIB
#ifdef MPI_ENV
USE MPIComm_mod, ONLY : REDUCE
#endif
USE TRAN_MOD,    ONLY : TranInfo,              TranCntl
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_Type) :: PE
REAL :: XsKinf
LOGICAL :: lXsLib, lCritSpec

TYPE(FxrInfo_TYpe), POINTER :: Fxr(:, :), myFxr
TYPE(Cell_Type), POINTER :: CellInfo(:)
Type(Pin_Type), POINTER :: Pin(:)

REAL, POINTER :: Phis(:, :, :)
REAL, POINTER :: SpecConv(:)
REAL, POINTER :: hz(:)

REAL :: PhiFxr(1000), xs, xl, xsv(128), xlv(128), area
REAL, ALLOCATABLE :: xsA(:), xsNF(:)
#ifdef MPI_ENV
REAL :: Buf(2), Buf0(2)
#endif

INTEGER :: nFxr, nFsr, nxy, myzb, myze, ng, norg
INTEGER :: iz, ixy, ig, icel, itype
INTEGER :: iFxr, iFsr, ifsrlocal, nlocalFxr, nFsrInFxr, FxrIdxSt, FsrIdxSt
INTEGER :: i, j, iResoGrpBeg, iResoGrpEnd
INTEGER :: tid

nFxr = Core%nCoreFxr; nFsr = Core%nCoreFsr
nxy = Core%nxy
myzb = PE%myzb; myze = PE%myze
Pin => Core%Pin
CellInfo => Core%CellInfo; hz => Core%hz;

ng = GroupInfo%ng
Fxr => FmInfo%Fxr; Phis => FmInfo%Phis
SpecConv => FmInfo%SpecConv

IF(lxslib) THEN
  iResoGrpBeg = GroupInfo%nofg + 1
  iResoGrpEnd = GroupInfo%nofg + GroupInfo%norg
  norg = GroupInfo%norg
ENDIF
DO tid = 1, 128
  IF(.NOT. XsMacBase(tid)%lAlloc) THEN
    XsMacBase(tid)%ng = ng
    CALL AllocXsMac(XsMacBase(tid))
ENDIF
END DO
xsv = 0; xlv = 0
tid = 1

DO iz = myzb, myze
  IF(.NOT. Core%lFuelPlane(iz)) CYCLE
!  !$OMP PARALLEL DEFAULT(SHARED) &
!  !$OMP PRIVATE(ixy,FsrIdxSt,FxrIdxSt,icel,nlocalFxr,j,ifxr,nFsrInFxr,myFxr,ig,itype,PhiFxr,i,ifsrlocal,iFsr,Area,xsA,xsNF)
!  !$ tid = omp_get_thread_num()+1
  ALLOCATE(xsA(ng),xsNF(ng))
!  !$OMP DO SCHEDULE(GUIDED)
  DO ixy = 1, nxy
    IF(Pin(ixy)%lRadRef) CYCLE
    FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
    icel = Pin(ixy)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
      myFxr => Fxr(ifxr, iz)
      IF(lXsLib) Then
        CALL MacXsANf(xsA,xsNF,myFxr,1,ng,ng)
        IF(myFxr%lres) THEN
          do ig = iResoGrpBeg,iResoGrpEnd
            xsA(ig) = xsA(ig) * myFxr%FresoA(ig)
            xsNF(ig) = xsNF(ig) * myFxr%FresoNF(ig)
          enddo
        ENDIF
      ELSE
        !itype = CellInfo(icel)%iReg(ifsrlocal)
        itype = myFxr%imix
        IF(TranCntl%lDynamicBen) THEN
          CALL xsbaseDynBen(itype, TranInfo%fuelTemp(ixy, iz), 1, ng, 1, ng, .FALSE., XsMacBase(tid))
        ELSE
          CALL xsbaseBen(itype, 1, ng, 1, ng, .FALSE., XsMacBase(tid))
        END IF
        xsA = XsMacBase(tid)%XsMacA(1:ng)
        xsNF = XsMacBase(tid)%XsMacNF(1:ng)
      ENDIF
      PhiFxr = 0
      DO i = 1, nFsrInFxr
        ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(i, j)
        iFsr = FsrIdxSt + ifsrlocal - 1
        Area = CellInfo(icel)%vol(ifsrlocal)
        PhiFxr(1:ng) = PhiFxr(1:ng) + Area * Phis(iFsr, iz, 1:ng)
      ENDDO
      IF(lCritSpec) THEN
        CALL MULTI_VA(SpecConv(1:ng), PhiFxr(1:ng), ng)
      ENDIF
      xlv(tid) = xlv(tid) + DotProduct(xsA(1:ng), PhiFxr(1:ng), ng) * hz(iz)
      xsv(tid) = xsv(tid) + DotProduct(xsNF(1:ng), PhiFxr(1:ng), ng) * hz(iz)
    ENDDO
  ENDDO
!  !$OMP END DO
  DEALLOCATE(xsA,xsNF)
!  !$OMP END PARALLEL
ENDDO
xl = SUM(xlv); xs = SUM(xsv)
#ifdef MPI_ENV
Buf0 = (/xl, xs/)
CALL REDUCE(buf0, buf, 2, PE%MPI_RTMASTER_COMM, .TRUE.)
xl = Buf(1); xs = Buf(2)
#endif
XsKinf = xs / xl
NULLIFY(PIN, CellInfo, hz)
NULLIFY(Fxr, Phis, SpecConv)

END FUNCTION

!SUBROUTINE PelletPowerDist(Core, FmInfo, CmInfo, GroupInfo, nTracerCntl, PE, istep)
!USE PARAM
!USE TYPEDEF,       ONLY : CoreInfo_Type,      CmInfo_Type,                     &
!                          FmInfo_Type,        GroupInfo_Type,     PE_TYPE,     &
!                          FxrInfo_Type,       Cell_Type,          Pin_Type
!USE CNTL,          ONLY : nTracerCntl_Type
!USE Depl_mod,         ONLY : FluxNormalizeFactor
!IMPLICIT NONE
!TYPE(CoreInfo_Type) :: Core
!TYPE(FmInfo_Type) :: FmInfo
!TYPE(CmInfo_Type) :: CmInfo
!TYPE(GroupInfo_Type) :: GroupInfo
!TYPE(nTracerCntl_Type) :: nTracerCntl
!TYPE(PE_TYPE) :: PE
!
!TYPE(FxrInfo_Type), POINTER :: Fxr(:, :), myFxr
!TYPE(Pin_Type), POINTER :: Pin(:)
!TYPE(Cell_Type), POINTER :: CellInfo(:)
!
!REAL, POINTER :: Phis(:, :, :)
!REAL, POINTER :: PinVol(:, :)
!REAL, POINTER :: Hz(:)
!
!REAL, POINTER :: XsMacA(:), XsMacNf(:), XsMacT(:), XsMacTr(:), chi(:)
!REAL, POINTER :: XsMacSm(:, :), XsMacP1Sm(:, :)
!REAL :: PhiFxr(MaxGrp), PhiSum(MaxGrp)
!REAL :: Area, AreaSum, Fis, FisSum
!
!REAL :: SigSum, SigSumP1
!REAL :: RPDist(100), FastFlux,NormFactor, PinVolsum
!
!INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nLocalFsr, nFsrInFxr
!INTEGER :: nFsr, nFxr, nxy, nz, myzb, myze, norg, nchi
!INTEGER :: ng, iResoGrpBeg, iResoGrpEnd
!INTEGER :: iz, ipin, icel, ifxr, ifsr, ifsrlocal
!INTEGER :: i, j, ig, ig2,istep, igb, ige
!LOGICAL :: lXsLib
!
!lXsLib = nTracerCntl%lXsLib
!IF(.NOT. lXsLib) RETURN
!
!Pin => Core%Pin; CellInfo => Core%CellInfo
!PinVol => Core%PinVol; Hz => Core%Hz
!
!Fxr => FmInfo%Fxr; Phis => FmInfo%Phis
!
!ng = GroupInfo%ng
!iresoGrpBeg = GroupInfo%nofg + 1; iresoGrpEnd = GroupInfo%nofg + GroupInfo%norg
!norg = GroupInfo%norg;  nchi = GroupInfo%nchi
!
!myzb = PE%myzb; myze = PE%myze
!nFsr = Core%nCoreFsr;   nFxr = Core%nCoreFxr
!nxy = Core%nxy
!
!
!FisSum = 0; phisum = 0
!!CALL
!RPDist = 0
!FastFlux = 0;pinvolsum = 0
!
!NormFactor = FluxNormalizeFactor(Core, FmInfo, GroupInfo, nTracerCntl%PowerCore, FALSE, TRUE, PE)
!NormFactor = NormFactor * nTracerCntl%PowerLevel
!igb = GroupInfo%GCStruct(1, 1); ige = GroupInfo%GCStruct(2, 1)
!DO iz = myzb, myze
!  !DO ipin = 1, nxy
!  IF(.NOT. Core%lFuelPlane(iz)) CYCLE
!  ipin = 1
!  !FastFlux = FastFlux + Core%Pinvol(ipin, iz)*CmInfo%GcPhiC(ipin, iz, 1)
!    FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
!    icel = Pin(ipin)%Cell(iz)
!    nlocalFxr = CellInfo(icel)%nFxr
!    DO j = 1, nLocalFxr
!      ifxr = FxrIdxSt + j -1
!      myFxr => Fxr(ifxr, iz)
!      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
!      PhiFxr(1:ng) = 0
!      AreaSum = 0
!      DO i = 1, nFsrInFxr
!        ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(i, j)
!        iFsr = FsrIdxSt + ifsrlocal - 1
!        Area = CellInfo(icel)%vol(ifsrlocal)
!        AreaSum = AreaSum + Area
!        PhiFxr(1:ng) = PhiFxr(1:ng) + Area * Phis(iFsr, iz, 1:ng)
!      ENDDO
!      IF(myFxr%ldepl) THEN
!        DO ig = 1, 4 !groupinfo%nchi
!          FastFlux = FastFlux + PhiFxr(ig)*hz(iz)
!        ENDDO
!        PinVolsum = PinVolsum + Areasum * hz(iz)
!      ENDIF
!
!      PhiFxr(1:ng) = PhiFxr(1:ng)/AreaSum
!
!      CALL MacXsBase(XsMacBase, myFxr, 1, ng, ng, 1._8, FALSE, TRUE)
!      CALL MacXsScatMatrix(XsMacBase, myFxr, 1, ng, ng, GroupInfo, TRUE)
!      IF(myFxr%lres) THEN
!        do ig = iResoGrpBeg, iResoGrpEnd
!          XsMacBase%XsMacA(ig) = XsMacBase%XsMacA(ig) * myFxr%FresoA(ig)
!          XsMacBase%XsMacNf(ig) = XsMacBase%XsMacNf(ig) * myFxr%FresoNF(ig)
!          XsMacBase%XsMacKf(ig) = XsMacBase%XsMacKf(ig) * myFxr%FresoF(ig)
!        ENDDO
!
!      ENDIF
!      CALL MacP1XsScatMatrix(XsMacBase, myFxr, 1, ng, ng, GroupInfo)
!
!      DO ig= 1, ng
!        RPDist(j) = RPDist(j) + PhiFxr(ig) * XsMacBase%XsMackf(ig)
!      ENDDO
!
!      !Added Reaction values
!!      XsMacA(1:ng) = XsMacA(1:ng) + PhiFxr(1:ng) * XsMacBase%XsMacA(1:ng)
!!      XsMacNf(1:ng) = XsMacNf(1:ng) + PhiFxr(1:ng) * XsMacBase%XsMacNf(1:ng)
!!      DO ig = 1, ng
!!        DO ig2 = 1, ng
!!          XsMacSm(ig2, ig) = XsMacSm(ig2, ig) + XsMacBase%XsMacSm(ig2, ig) * PhiFxr(ig2)
!!          XsMacP1Sm(ig2, ig) = XsMacP1Sm(ig2, ig) + XsMacBase%XsMacP1Sm(ig2, ig) * PhiFxr(ig2)
!!        ENDDO
!!      ENDDO
!!      IF(myFxr%ldepl) THEN
!!        Fis = DotProduct(PhiFxr(1:ng), XsMacBase%XsMacNf(1:ng), ng)
!!        DO ig = 1, nchi
!!          CHI(ig) = CHI(ig) + Fis * myFxr%Chi(ig)
!!        ENDDO
!!        FisSum = FisSum + Fis
!!      ENDIF
!!      PhiSum(1:ng) = PhiSum(1:ng) + PhiFxr(1:ng)
!  ENDDO
!!  ENDDO
!ENDDO
!FastFlux =FastFlux/PinVolsum * NormFactor
!!Flux Normalization
!!WRITE(88, '(I10, 200e15.4)') istep, (RPDIst(j), j=1,nLocalFxr)
!!WRITE(89, '(I10, 200e15.4)') istep, FastFlux
!!WRITE(*, '(I10, 200e15.4)') istep, (RPDIst(j), j=1,nLocalFxr)
!NULLIFY(myFxr)
!NULLIFY(Pin, CellInfo, Pinvol, Hz)
!NULLIFY(Fxr, Phis)
!END SUBROUTINE

END MODULE
