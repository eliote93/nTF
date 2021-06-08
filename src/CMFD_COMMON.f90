MODULE CMFD_COMMON

USE PARAM
USE TYPEDEF,        ONLY : XsMac_Type
USE TYPEDEF_COMMON
USE MOC_COMMON,     ONLY : CoreXsMac,           InScatRange,        InScatIdx
USE OMP_LIB
IMPLICIT NONE

TYPE(XsMac_Type), POINTER, PRIVATE :: XsMac(:, :)

CONTAINS

SUBROUTINE AllocHomoXSVar(CoreInfo, ng)
USE PARAM
USE CNTL,           ONLY : nTracerCntl
USE TYPEDEF,        ONLY : CoreInfo_Type,       Cell_Type
USE CMFD_MOD,       ONLY : CMFDPinXS
USE PE_MOD,         ONLY : PE
USE XsUtil_mod,     ONLY : AllocXsMac
USE HexData,        ONLY : ncTyp
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
INTEGER :: ng

TYPE(Cell_Type), POINTER :: Cell(:)
INTEGER :: ifxr, icel, ipin, tid
INTEGER :: nCellType, nFxrMax

Cell => CoreInfo%CellInfo
nCellType = CoreInfo%nCellType

nFxrMax = 0

IF(nTracerCntl%lHex) THEN
  DO icel = 1, ncTyp
    IF (.NOT. Cell(icel)%luse) CYCLE

    nFxrMax = max(nFxrMax, Cell(icel)%nFxr)
  END DO

  nFxrMax = nFxrMax * 3
ELSE
  DO icel = 1, nCellType
    nFxrMax = max(nFxrMax, Cell(icel)%nFxr)
  END DO

  nFxrMax = nFxrMax + 3
ENDIF

CALL omp_set_num_threads(PE%nCMFDThread)

ALLOCATE(XsMac(nFxrMax, PE%nCMFDThread))

!$OMP PARALLEL PRIVATE(tid)
tid = omp_get_thread_num() + 1
DO ifxr = 1, nFxrMax
  XsMac(ifxr, tid)%ng = ng
  CALL AllocXsMac(XsMac(ifxr, tid))
ENDDO
!$OMP END PARALLEL

END SUBROUTINE
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE AllocPinXS(PinXS, GroupInfo, nxy, myzb, myze)

USE allocs
USE geom,    ONLY : ncbd
USE TYPEDEF, ONLY : GroupInfo_Type, PinXS_Type
USE CNTL,    ONLY : nTracerCntl

IMPLICIT NONE

TYPE(PinXS_Type), POINTER, DIMENSION(:,:) :: PinXS
TYPE(GroupInfo_Type) :: GroupInfo
INTEGER :: ng, nxy, myzb, myze

TYPE(PinXS_Type), POINTER :: myPinXS
INTEGER :: ig, igb, ige, ipin, iz
! ----------------------------------------------------

ng = GroupInfo%ng

ALLOCATE (PinXS(nxy, myzb:myze))

DO iz = myzb, myze
  DO ipin = 1, nxy
    myPinXS => PinXS(ipin, iz)
    
    CALL dmalloc(myPinXS%Dtil, ncbd, ng)
    CALL dmalloc(myPinXS%Dhat, ncbd, ng)
    CALL dmalloc(myPinXS%XSD,  ng)
    CALL dmalloc(myPinXS%XSD2, ng)
    CALL dmalloc(myPinXS%XSt,  ng)
    CALL dmalloc(myPinXS%XStr, ng)
    CALL dmalloc(myPinXS%XSr,  ng)
    CALL dmalloc(myPinXS%XSa,  ng)
    CALL dmalloc(myPinXS%XSnf, ng)
    CALL dmalloc(myPinXS%XSkf, ng)
    CALL dmalloc(myPinXS%Chi,  ng)
    CALL dmalloc(myPinXS%Phi,  ng)
    
    ALLOCATE (myPinXS%XSs (ng))
    
    DO ig = 1, ng
      igb = GroupInfo%InScatRange(1, ig)
      ige = GroupInfo%InScatRange(2, ig)
      
      CALL dmalloc0(myPinXS%XSs(ig)%from, igb, ige)
      
      myPinXS%XSs(ig)%ib = igb
      myPinXS%XSs(ig)%ie = ige
    END DO
    
    IF (.NOT. nTracerCntl%lDomainDcmp) CYCLE
    
    CALL dmalloc(myPinXS%atil, ncbd, ng)
    CALL dmalloc(myPinXS%ahat, ncbd, ng)
  END DO
END DO
! ----------------------------------------------------

END SUBROUTINE AllocPinXS
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE CopyPinXS(Input, Output, ng)
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type) :: Input, Output
INTEGER :: ig, ng

Output%Dtil = Input%Dtil
Output%Dhat = Input%Dhat
Output%XSD = Input%XSD
Output%XSD2 = Input%XSD2
Output%XSt = Input%XSt
Output%XStr = Input%XStr
Output%XSr = Input%XSr
Output%XSa = Input%XSa
Output%XSnf = Input%XSnf
Output%XSkf = Input%XSkf
Output%Chi = Input%Chi

DO ig = 1, ng
  Output%XSs(ig)%from = Input%XSs(ig)%from
ENDDO

END SUBROUTINE

SUBROUTINE HomogenizeXS(CoreInfo, superPin, Fxr, PinXS, phis, ng, nxy, myzb, myze, lxslib, lscat1, lsigt)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,       FxrInfo_Type,       PinXS_Type,                                     &
                           Pin_Type,            Cell_Type
USE CORE_MOD,       ONLY : GroupInfo
USE CNTL,           ONLY : nTracerCntl
USE BenchXs,        ONLY : xsbaseBen,           xsbaseDynBen
USE MacXsLib_Mod,   ONLY : MacXsBase,           MacXsScatMatrix
USE TRAN_MOD,       ONLY : TranInfo
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(superPin_Type), POINTER :: superPin(:)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: phis(:, :, :)
INTEGER :: ng, nxy, myzb, myze
LOGICAL :: lxslib, lscat1, lsigt

TYPE(FxrInfo_Type), POINTER :: myFxr
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)

INTEGER :: i, j, ixy, icel, ipin, iz, ifxr, ifxr_global, itype, tid
INTEGER :: FsrIdxSt, FxrIdxSt
INTEGER :: nChi, norg, nLocalFsr, nLocalFxr, nFsrInFxr
INTEGER :: irgb, irge

IF (nTracerCntl%lMacro) THEN
  CALL HomogenizeMacroXS(CoreInfo, superPin, Fxr, PinXS, phis, ng, nxy, myzb, myze, lxslib, lscat1, lsigt)
  RETURN
ENDIF

Pin => CoreInfo%Pin
Cell => CoreInfo%CellInfo
nChi = GroupInfo%nChi

IF (lxslib) THEN
  norg = GroupInfo%norg
  irgb = GroupInfo%nofg + 1
  irge = GroupInfo%nofg + GroupInfo%norg
ENDIF

!$OMP PARALLEL PRIVATE(tid, ipin, icel, ifxr, ifxr_global, itype, FsrIdxSt, FxrIdxSt, nLocalFxr, myFxr)
tid = omp_get_thread_num() + 1

DO iz = myzb, myze
  !$OMP DO SCHEDULE(GUIDED)
  DO ixy = 1, nxy
    nLocalFxr = 0
    DO j = 1, superPin(ixy)%nxy
      ipin = superPin(ixy)%pin(j)
      icel = Pin(ipin)%Cell(iz)
      FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
      DO i = 1, Cell(icel)%nFxr
        ifxr = nLocalFxr + i; ifxr_global = FxrIdxSt + i - 1; nFsrInFxr = Cell(icel)%nFsrInFxr(i)
        myFxr => Fxr(ifxr_global, iz)
        XsMac(ifxr, tid)%lFuel = FALSE
        IF (lxslib) THEN
          CALL MacXsBase(XsMac(ifxr, tid), myFxr, 1, ng, ng, 1.0D0, FALSE, FALSE, TRUE)
          CALL MacXsScatMatrix(XsMac(ifxr, tid), myFxr, 1, ng, ng, GroupInfo, FALSE, TRUE)
          IF (myFxr%lres) THEN
            XsMac(ifxr, tid)%XsMacA(irgb : irge) = myFxr%FresoA(irgb : irge) * XsMac(ifxr, tid)%XsMacA(irgb : irge)
            IF (CoreInfo%lFuelPlane(iz)) THEN
              XsMac(ifxr, tid)%XsMacNf(irgb : irge) = myFxr%fresonf(irgb : irge) * XsMac(ifxr, tid)%XsMacNf(irgb : irge)
              XsMac(ifxr, tid)%XsMacKf(irgb : irge) = myFxr%fresokf(irgb : irge) * XsMac(ifxr, tid)%XsMacKf(irgb : irge)
            ENDIF
          ENDIF
          XsMac(ifxr, tid)%XsMacTr = XsMac(ifxr, tid)%XsMacA + XsMac(ifxr, tid)%XsMacStr
          XsMac(ifxr, tid)%XsMacT = XsMac(ifxr, tid)%XsMacA + XsMac(ifxr, tid)%XsMacS
          XsMac(ifxr, tid)%Chi = 0.0
          IF (myFxr%lDepl) THEN
            XsMac(ifxr, tid)%Chi(1 : nChi) = myFxr%Chi
            XsMac(ifxr, tid)%lFuel = TRUE
          ENDIF
        ELSE
          itype = myFxr%imix
          IF(nTracerCntl%lDynamicBen) THEN
            CALL xsbaseDynBen(itype, TranInfo%fuelTemp(ipin, iz), 1, ng, 1, ng, lscat1, XsMac(ifxr, tid))
          ELSE
            CALL xsbaseBen(itype, 1, ng, 1, ng, lscat1, XsMac(ifxr, tid))
          END IF
        ENDIF
      ENDDO
      nLocalFxr = nLocalFxr + Cell(icel)%nFxr
    ENDDO
    CALL HomogenizeCellXS(CoreInfo, superPin(ixy), PinXS(ixy, iz), XsMac(1 : nLocalFxr, tid), phis, iz, ng, lxslib, lsigt)
  ENDDO
  !$OMP END DO
ENDDO

!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE HomogenizeMacroXS(CoreInfo, superPin, Fxr, PinXS, phis, ng, nxy, myzb, myze, lxslib, lscat1, lsigt)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,       FxrInfo_Type,       PinXS_Type,                                     &
                           Pin_Type,            Cell_Type
USE CNTL,           ONLY : nTracerCntl
USE CORE_MOD,       ONLY : GroupInfo
USE BenchXs,        ONLY : xsbaseBen,           xsbaseDynBen,       xsbaseben_neacrp
USE TRAN_MOD,       ONLY : TranInfo,            TranCntl
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(superPin_Type), POINTER :: superPin(:)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: phis(:, :, :)
INTEGER :: ng, nxy, myzb, myze
LOGICAL :: lxslib, lscat1, lsigt

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
INTEGER :: i, j, ig, igf, igs, ixy, iz, ipin, icel, itype, ifxr, ifxr_global, tid
INTEGER :: FxrIdxSt, nLocalFxr, nChi

Pin => CoreInfo%Pin
Cell => CoreInfo%CellInfo
nChi = GroupInfo%nChi

!$OMP PARALLEL PRIVATE(tid, ipin, icel, igs, itype, ifxr, ifxr_global, FxrIdxSt, nLocalFxr)
tid = omp_get_thread_num() + 1

DO iz = myzb, myze
  !$OMP DO SCHEDULE(GUIDED)
  DO ixy = 1, nxy
    nLocalFxr = 0
    DO j = 1, superPin(ixy)%nxy
      ipin = superPin(ixy)%pin(j)
      icel = Pin(ipin)%Cell(iz)
      FxrIdxSt = Pin(ipin)%FxrIdxSt
      DO i = 1, Cell(icel)%nFxr
        ifxr = nLocalFxr + i; ifxr_global = FxrIdxSt + i - 1
        IF (lxslib) THEN
          XsMac(ifxr, tid)%XsMacT = CoreXsMac(iz)%XSt(:, ifxr_global)
          XsMac(ifxr, tid)%XsMacTr = CoreXsMac(iz)%XStr(:, ifxr_global)
          XsMac(ifxr, tid)%XsMacA = CoreXsMac(iz)%XSa(:, ifxr_global)
          XsMac(ifxr, tid)%XsMacNf = CoreXsMac(iz)%XSnf(:, ifxr_global)
          XsMac(ifxr, tid)%XsMacKf = CoreXsMac(iz)%XSkf(:, ifxr_global)
          DO ig = 1, ng
            DO igf = InScatRange(1, ig), InScatRange(2, ig)
              igs = InScatIdx(igf, ig)
              XsMac(ifxr, tid)%XsMacSm(igf, ig) = CoreXsMac(iz)%XSsm(igs, ifxr_global)
            ENDDO
          ENDDO
          XsMac(ifxr, tid)%Chi = 0.0
          XsMac(ifxr, tid)%lFuel = FALSE
          IF (Fxr(ifxr_global, iz)%lDepl) THEN
            XsMac(ifxr, tid)%Chi(1 : nChi) = Fxr(ifxr_global, iz)%Chi
            XsMac(ifxr, tid)%lFuel = TRUE
          ENDIF
        ELSE
          itype = Fxr(ifxr_global, iz)%imix
          IF(nTracerCntl%libtyp .EQ. 11) THEN
            CALL xsbaseben_NEACRP(itype, Fxr(ifxr_global,iz)%rho, Fxr(ifxr_global,iz)%temp, Fxr(ifxr_global,iz)%DopTemp, XsMac(ifxr, tid))
          ELSE
            IF(nTracerCntl%lDynamicBen) THEN
              CALL xsbaseDynBen(itype, TranInfo%fuelTemp(ipin, iz), 1, ng, 1, ng, lscat1, XsMac(ifxr, tid))
            ELSE
              CALL xsbaseBen(itype, 1, ng, 1, ng, lscat1, XsMac(ifxr, tid))
            END IF
          ENDIF
        END IF
      ENDDO
      nLocalFxr = nLocalFxr + Cell(icel)%nFxr
    ENDDO
    CALL HomogenizeCellXS(CoreInfo, superPin(ixy), PinXS(ixy, iz), XsMac(1 : nLocalFxr, tid), phis, iz, ng, lxslib, lsigt)
  ENDDO
  !$OMP END DO
ENDDO

!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE HomogenizeXS_Cusping(CoreInfo, FmInfo, superPin, Fxr, PinXS, phis, ng, nxy, myzb, myze, nthread, lxslib, lscat1, lsigt)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,       FxrInfo_Type,       PinXS_Type,                                     &
                           Pin_Type,            Cell_Type,          FmInfo_Type
USE CORE_MOD,       ONLY : GroupInfo
USE CNTL,           ONLY : nTracerCntl
USE BenchXs,        ONLY : xsbaseBen,           MacXsBen,           xsbaseBen_Cusping,    &
                           xsbaseDynBen,        xsbaseDynBen_Cusping, DynMacXsBen,        &
                           xsbaseBen_NEACRP
USE MacXsLib_Mod,   ONLY : MacXsBase,           MacXsScatMatrix
USE TRAN_MOD,       ONLY : TranInfo
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(superPin_Type), POINTER :: superPin(:)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: phis(:, :, :)
INTEGER :: ng, nxy, myzb, myze, nthread
LOGICAL :: lxslib, lscat1, lsigt

TYPE(FxrInfo_Type), POINTER :: myFxr
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)

REAL :: phiz(ng, nthread), philz(ng, nthread), phiuz(ng, nthread)
REAL :: vol(nthread), volsum(nthread)
INTEGER :: ig, ifsrlocal, ifsr
INTEGER :: i, j, k, ixy, icel, ipin, iz, ifxr, ifxr_global, itype, tid
INTEGER :: FsrIdxSt, FxrIdxSt
INTEGER :: nChi, norg, nLocalFsr, nLocalFxr, nFsrInFxr
INTEGER :: irgb, irge
LOGICAL :: lCusping

IF (nTracerCntl%lMacro) THEN
  CALL HomogenizeMacroXS_Cusping(CoreInfo, FmInfo, superPin, Fxr, PinXS, phis, ng, nxy, myzb, myze, nthread, lxslib, lscat1, lsigt)
  RETURN
ENDIF

Pin => CoreInfo%Pin
Cell => CoreInfo%CellInfo
nChi = GroupInfo%nChi

IF (lxslib) THEN
  norg = GroupInfo%norg
  irgb = GroupInfo%nofg + 1
  irge = GroupInfo%nofg + GroupInfo%norg
ENDIF

!$OMP PARALLEL PRIVATE(tid, ipin, icel, ifxr, ifxr_global, itype, FsrIdxSt, nFsrInFxr, FxrIdxSt, nLocalFxr, myFxr, ifsrlocal, ifsr)
tid = omp_get_thread_num() + 1

DO iz = myzb, myze
  !$OMP DO SCHEDULE(GUIDED)
  DO ixy = 1, nxy
    nLocalFxr = 0
    DO j = 1, superPin(ixy)%nxy
      ipin = superPin(ixy)%pin(j)
      icel = Pin(ipin)%Cell(iz)
      FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
      DO i = 1, Cell(icel)%nFxr
        ifxr = nLocalFxr + i; ifxr_global = FxrIdxSt + i - 1; nFsrInFxr = Cell(icel)%nFsrInFxr(i)
        myFxr => Fxr(ifxr_global, iz)
        XsMac(ifxr, tid)%lFuel = FALSE
        IF (lxslib) THEN
          CALL MacXsBase(XsMac(ifxr, tid), myFxr, 1, ng, ng, 1.0D0, FALSE, FALSE, TRUE)
          CALL MacXsScatMatrix(XsMac(ifxr, tid), myFxr, 1, ng, ng, GroupInfo, FALSE, TRUE)
          IF (myFxr%lres) THEN
            XsMac(ifxr, tid)%XsMacA(irgb : irge) = myFxr%FresoA(irgb : irge) * XsMac(ifxr, tid)%XsMacA(irgb : irge)
            IF (CoreInfo%lFuelPlane(iz)) THEN
              XsMac(ifxr, tid)%XsMacNf(irgb : irge) = myFxr%fresonf(irgb : irge) * XsMac(ifxr, tid)%XsMacNf(irgb : irge)
              XsMac(ifxr, tid)%XsMacKf(irgb : irge) = myFxr%fresokf(irgb : irge) * XsMac(ifxr, tid)%XsMacKf(irgb : irge)
            ENDIF
          ENDIF
          XsMac(ifxr, tid)%XsMacTr = XsMac(ifxr, tid)%XsMacA + XsMac(ifxr, tid)%XsMacStr
          XsMac(ifxr, tid)%XsMacT = XsMac(ifxr, tid)%XsMacA + XsMac(ifxr, tid)%XsMacS
          XsMac(ifxr, tid)%Chi = 0.0
          IF (myFxr%lDepl) THEN
            XsMac(ifxr, tid)%Chi(1 : nChi) = myFxr%Chi
            XsMac(ifxr, tid)%lFuel = TRUE
          ENDIF
        ELSE
          itype = myFxr%imix
          IF(nTracerCntl%lDynamicBen) THEN
            lCusping = DynMacXsBen(itype)%lCusping
          ELSE
            lCusping = MacXsBen(itype)%lCusping
          END IF
          IF(lCusping) THEN
            phiz(:,tid) = 0.
            philz(:, tid) = 0.
            phiuz(:, tid) = 0.
            volsum(tid) = 0.
            DO ig = 1, ng
              DO k = 1, nfsrinfxr
                ifsrlocal = Cell(icel)%MapFxr2FsrIdx(k, i)
                ifsr = FsrIdxSt + ifsrlocal - 1
                vol(tid) = Cell(icel)%Vol(ifsrlocal)
                IF(ig .EQ. 1) volsum(tid) = volsum(tid) + vol(tid)
                phiz(ig, tid) = Phiz(ig, tid) + Phis(ifsr, iz, ig) * vol(tid)
                IF(iz .EQ. myzb) THEN
                  philz(ig, tid) = philz(ig, tid) + FmInfo%neighphis(ifsr, ig, BOTTOM) * vol(tid)
                ELSE
                  philz(ig, tid) = philz(ig, tid) + phis(ifsr, iz-1, ig) * vol(tid)
                END IF
                IF(iz .EQ. myze) THEN
                  phiuz(ig, tid) = phiuz(ig, tid) + FmInfo%neighphis(ifsr, ig, TOP) * vol(tid)
                ELSE
                  phiuz(ig, tid) = phiuz(ig, tid) + phis(ifsr, iz+1, ig) * vol(tid)
                END IF
              END DO
              phiz(ig, tid) = phiz(ig, tid) / volsum(tid)
              philz(ig, tid) = philz(ig, tid) / volsum(tid)
              phiuz(ig, tid) = phiuz(ig, tid) / volsum(tid)
            END DO
            IF(nTracerCntl%libtyp .EQ. 11) THEN
              CALL xsbaseben_NEACRP(itype, myFxr%rho, myFxr%temp, myFxr%DopTemp, XsMac(ifxr, tid))
            ELSE
              IF(nTracerCntl%lDynamicBen) THEN
                CALL xsbaseDynBen_Cusping(itype, TranInfo%fuelTemp(ipin, iz), 1, ng, 1, ng, lscat1, XsMac(ifxr, tid),&
                  phiz(:,tid), philz(:,tid), phiuz(:,tid), CoreInfo%hzfm(iz), CoreInfo%hzfm(iz-1), CoreInfo%hzfm(iz+1))
              ELSE
                CALL xsbaseBen_Cusping(itype, 1, ng, 1, ng, lscat1, XsMac(ifxr, tid),&
                  phiz(:,tid), philz(:,tid), phiuz(:,tid), CoreInfo%hzfm(iz), CoreInfo%hzfm(iz-1), CoreInfo%hzfm(iz+1))
              END IF
            END IF
          ELSE
            IF(nTracerCntl%libtyp .EQ. 11) THEN
              CALL xsbaseben_NEACRP(itype, myFxr%rho, myFxr%temp, myFxr%DopTemp, XsMac(ifxr, tid))
            ELSE
              IF(nTracerCntl%lDynamicBen) THEN
                CALL xsbaseDynBen(itype, TranInfo%fuelTemp(ipin, iz), 1, ng, 1, ng, lscat1, XsMac(ifxr, tid))
              ELSE
                CALL xsbaseBen(itype, 1, ng, 1, ng, lscat1, XsMac(ifxr, tid))
              END IF
            ENDIF
          END IF
        ENDIF
      ENDDO
      nLocalFxr = nLocalFxr + Cell(icel)%nFxr
    ENDDO
    CALL HomogenizeCellXS(CoreInfo, superPin(ixy), PinXS(ixy, iz), XsMac(1 : nLocalFxr, tid), phis, iz, ng, lxslib, lsigt)
  ENDDO
  !$OMP END DO
ENDDO

!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE HomogenizeMacroXS_Cusping(CoreInfo, FmInfo, superPin, Fxr, PinXS, phis, ng, nxy, myzb, myze, nthread, lxslib, lscat1, lsigt)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,       FxrInfo_Type,       PinXS_Type,                                     &
                           Pin_Type,            Cell_Type,          FmInfo_Type
USE CNTL,           ONLY : nTracerCntl
USE CORE_MOD,       ONLY : GroupInfo
USE BenchXs,        ONLY : xsbaseBen,           MacXsBen,           xsbaseBen_Cusping,  &
                           xsbaseDynBen,        xsbaseDynBen_Cusping, DynMacXsBen, xsbaseBen_NEACRP
USE TRAN_MOD,       ONLY : TranInfo,            TranCntl
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(FmInfo_Type) :: FmInfo
TYPE(superPin_Type), POINTER :: superPin(:)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: phis(:, :, :)
INTEGER :: ng, nxy, myzb, myze, nthread
LOGICAL :: lxslib, lscat1, lsigt

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
REAL :: phiz(ng, nthread), philz(ng, nthread), phiuz(ng, nthread)
REAL :: vol(nthread), volsum(nthread)
INTEGER :: nfsrinfxr, FsrIdxSt
INTEGER :: ifsrlocal, ifsr
INTEGER :: i, j, k, ig, igf, igs, ixy, iz, ipin, icel, itype, ifxr, ifxr_global, tid
INTEGER :: FxrIdxSt, nLocalFxr, nChi
LOGICAL :: lCusping

Pin => CoreInfo%Pin
Cell => CoreInfo%CellInfo
nChi = GroupInfo%nChi

!$OMP PARALLEL PRIVATE(tid, ipin, icel, igs, itype, ifxr, ifxr_global, FxrIdxSt, FsrIdxSt, nLocalFxr, nfsrinfxr, ifsrlocal, ifsr)
tid = omp_get_thread_num() + 1

DO iz = myzb, myze
  !$OMP DO SCHEDULE(GUIDED)
  DO ixy = 1, nxy
    nLocalFxr = 0
    DO j = 1, superPin(ixy)%nxy
      ipin = superPin(ixy)%pin(j)
      icel = Pin(ipin)%Cell(iz)
      FxrIdxSt = Pin(ipin)%FxrIdxSt
      FsrIdxSt = Pin(ipin)%FsrIdxSt
      DO i = 1, Cell(icel)%nFxr
        ifxr = nLocalFxr + i; ifxr_global = FxrIdxSt + i - 1
        nFsrInFxr = Cell(icel)%nFsrInFxr(i)
        IF (lxslib) THEN
          XsMac(ifxr, tid)%XsMacT = CoreXsMac(iz)%XSt(:, ifxr_global)
          XsMac(ifxr, tid)%XsMacTr = CoreXsMac(iz)%XStr(:, ifxr_global)
          XsMac(ifxr, tid)%XsMacA = CoreXsMac(iz)%XSa(:, ifxr_global)
          XsMac(ifxr, tid)%XsMacNf = CoreXsMac(iz)%XSnf(:, ifxr_global)
          XsMac(ifxr, tid)%XsMacKf = CoreXsMac(iz)%XSkf(:, ifxr_global)
          DO ig = 1, ng
            DO igf = InScatRange(1, ig), InScatRange(2, ig)
              igs = InScatIdx(igf, ig)
              XsMac(ifxr, tid)%XsMacSm(igf, ig) = CoreXsMac(iz)%XSsm(igs, ifxr_global)
            ENDDO
          ENDDO
          XsMac(ifxr, tid)%Chi = 0.0
          XsMac(ifxr, tid)%lFuel = FALSE
          IF (Fxr(ifxr_global, iz)%lDepl) THEN
            XsMac(ifxr, tid)%Chi(1 : nChi) = Fxr(ifxr_global, iz)%Chi
            XsMac(ifxr, tid)%lFuel = TRUE
          ENDIF
        ELSE
          itype = Fxr(ifxr_global, iz)%imix
          IF(TranCntl%lDynamicBen) THEN
            lCusping = DynMacXsBen(itype)%lCusping
          ELSE
            lCusping = MacXsBen(itype)%lCusping
          END IF
          IF(lCusping) THEN
            phiz(:,tid) = 0.
            philz(:, tid) = 0.
            phiuz(:, tid) = 0.
            volsum(tid) = 0.
            DO ig = 1, ng
              DO k = 1, nfsrinfxr
                ifsrlocal = Cell(icel)%MapFxr2FsrIdx(k, i)
                ifsr = FsrIdxSt + ifsrlocal - 1
                vol(tid) = Cell(icel)%Vol(ifsrlocal)
                IF(ig .EQ. 1) volsum(tid) = volsum(tid) + vol(tid)
                phiz(ig, tid) = Phiz(ig, tid) + Phis(ifsr, iz, ig) * vol(tid)
                IF(iz .EQ. myzb) THEN
                  philz(ig, tid) = philz(ig, tid) + FmInfo%neighphis(ifsr, ig, BOTTOM) * vol(tid)
                ELSE
                  philz(ig, tid) = philz(ig, tid) + phis(ifsr, iz-1, ig) * vol(tid)
                END IF
                IF(iz .EQ. myze) THEN
                  phiuz(ig, tid) = phiuz(ig, tid) + FmInfo%neighphis(ifsr, ig, TOP) * vol(tid)
                ELSE
                  phiuz(ig, tid) = phiuz(ig, tid) + phis(ifsr, iz+1, ig) * vol(tid)
                END IF
              END DO
              phiz(ig, tid) = phiz(ig, tid) / volsum(tid)
              philz(ig, tid) = philz(ig, tid) / volsum(tid)
              phiuz(ig, tid) = phiuz(ig, tid) / volsum(tid)
            END DO
            IF(TranCntl%lDynamicBen) THEN
              CALL xsbaseDynBen_Cusping(itype, TranInfo%fuelTemp(ipin, iz), 1, ng, 1, ng, lscat1, XsMac(ifxr, tid),&
                phiz(:,tid), philz(:,tid), phiuz(:,tid), CoreInfo%hzfm(iz), CoreInfo%hzfm(iz-1), CoreInfo%hzfm(iz+1))
            ELSE
              CALL xsbaseBen_Cusping(itype, 1, ng, 1, ng, lscat1, XsMac(ifxr, tid),&
                phiz(:,tid), philz(:,tid), phiuz(:,tid), CoreInfo%hzfm(iz), CoreInfo%hzfm(iz-1), CoreInfo%hzfm(iz+1))
            END IF
          ELSE
            IF(nTracerCntl%libtyp .EQ. 11) THEN
              !IF(ixy .EQ. 1) PRINT*, 'NEACRP base xs'
              CALL xsbaseben_NEACRP(itype, Fxr(ifxr_global,iz)%rho, Fxr(ifxr_global,iz)%temp, Fxr(ifxr_global,iz)%DopTemp, XsMac(ifxr, tid))
            ELSE
              IF(nTracerCntl%lDynamicBen) THEN
                CALL xsbaseDynBen(itype, TranInfo%fuelTemp(ipin, iz), 1, ng, 1, ng, lscat1, XsMac(ifxr, tid))
              ELSE
                CALL xsbaseBen(itype, 1, ng, 1, ng, lscat1, XsMac(ifxr, tid))
              END IF
            ENDIF
          END IF
        ENDIF
      ENDDO
      nLocalFxr = nLocalFxr + Cell(icel)%nFxr
    ENDDO
    CALL HomogenizeCellXS(CoreInfo, superPin(ixy), PinXS(ixy, iz), XsMac(1 : nLocalFxr, tid), phis, iz, ng, lxslib, lsigt)
  ENDDO
  !$OMP END DO
ENDDO

!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE HomogenizeCellXS(CoreInfo, superPin, PinXS, XsMac, phis, iz, ng, lxslib, lsigt)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,       Pin_Type,           Cell_Type,                                      &
                           PinXS_Type,          XsMac_Type
USE CORE_MOD,       ONLY : GroupInfo
USE CNTL,           ONLY : nTracerCntl
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(superPin_Type) :: superPin
TYPE(PinXS_Type) :: PinXS
TYPE(XsMac_Type) :: XsMac(:)
REAL, POINTER :: phis(:, :, :)
INTEGER :: iz, ng
LOGICAL :: lxslib, lsigt

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
INTEGER :: i, j, ixy, ifxr, ipin, icel, ig, ig2, igb, ige, ireg, ifsr
INTEGER :: nxy, nFxr, nFsrInFxr, nLocalFxr, FsrIdxSt
REAL :: localphi, localpsi, psi, vol, phisum, psisum, volsum
REAL :: RR(6), RRS(ng), Chi(ng)
LOGICAL :: lFuel

Pin => CoreInfo%Pin
Cell => CoreInfo%CellInfo
nxy = superPin%nxy

DO ig = 1, ng
  PinXS%XSs(ig)%from = 0.0
ENDDO

DO ig = 1, ng
  igb = GroupInfo%OutScatRange(1, ig)
  ige = GroupInfo%OutScatRange(2, ig)
  RR = 0.0; RRS = 0.0; phisum = 0.0; volsum = 0.0
  nLocalFxr = 0
  DO ixy = 1, nxy
    ipin = superPin%pin(ixy)
    icel = Pin(ipin)%Cell(iz)
    nFxr = Cell(icel)%nFxr
    FsrIdxSt = Pin(ipin)%FsrIdxSt
    DO i = 1, nFxr
      ifxr = nLocalFxr + i
      nFsrInFxr = Cell(icel)%nFsrInFxr(i)
      localphi = 0.0
      DO j = 1, nFsrInFxr
        ireg = Cell(icel)%MapFxr2FsrIdx(j, i)
        ifsr = FsrIdxSt + ireg - 1
        vol = Cell(icel)%vol(ireg)
        localphi = localphi + phis(ifsr, iz, ig) * vol
        volsum = volsum + vol
      ENDDO
      phisum = phisum + localphi
      RR(1) = RR(1) + localphi * XsMac(ifxr)%XsMacT(ig)
      RR(2) = RR(2) + localphi * XsMac(ifxr)%XsMacTr(ig)
      RR(3) = RR(3) + localphi * XsMac(ifxr)%XsMacNf(ig)
      RR(4) = RR(4) + localphi * XsMac(ifxr)%XsMacKf(ig)
      RR(5) = RR(5) + localphi * XsMac(ifxr)%XsMacA(ig)
      IF (lsigt) THEN
        RR(6) = RR(6) + localphi / XsMac(ifxr)%XsMacT(ig)
      ELSE
        RR(6) = RR(6) + localphi / XsMac(ifxr)%XsMacTr(ig)
      ENDIF
      DO ig2 = igb, ige
        RRS(ig2) = RRS(ig2) + localphi * XsMac(ifxr)%XsMacSm(ig, ig2)
      ENDDO
    ENDDO
    nLocalFxr = nLocalFxr + nFxr
  ENDDO
  PinXS%Phi(ig) = phisum / volsum
  RR = RR / phisum
  RRS = RRS / phisum
  PinXS%XSt(ig) = RR(1)
  PinXS%XStr(ig) = RR(2)
  PinXS%XSnf(ig) = RR(3)
  PinXS%XSkf(ig) = RR(4)
  PinXS%XSa(ig) = RR(5)
  PinXS%XSs(ig)%WithInGroupScat = RRS(ig)
  DO ig2 = igb, ige
    IF ((ig - PinXS%XSs(ig2)%ib) * (ig - PinXS%XSs(ig2)%ie) .GT. 0) CYCLE
    PinXS%XSs(ig2)%from(ig) = RRS(ig2)
  ENDDO
  PinXS%XSr(ig) = RR(5) + sum(RRS) - RRS(ig)
  PinXS%XSs(ig)%self = PinXS%XSs(ig)%from(ig)
  PinXS%XSs(ig)%from(ig) = 0.0
  PinXs%XSD(ig) = 1.0 / 3.0 / PinXS%XStr(ig)
  PinXS%XSD2(ig) = 3.0 / 7.0 / PinXS%XSt(ig)
  IF (lsigt) THEN
    PinXS%XStr(ig) = PinXS%XSt(ig)
    PinXs%XSD(ig) = 1.0 / 3.0 / PinXS%XSt(ig)
  ENDIF
  IF (nTracerCntl%lDhom) THEN
    PinXs%XSD(ig) = 1.0 / 3.0 * RR(6)
  ENDIF
ENDDO

lFuel = FALSE
IF (.NOT. lxslib) THEN
  nLocalFxr = 0
  DO ixy = 1, nxy
    ipin = superPin%pin(ixy)
    icel = Pin(ipin)%Cell(iz)
    nFxr = Cell(icel)%nFxr
    DO i = 1, nFxr
      ifxr = nLocalFxr + i
      IF (XsMac(ifxr)%lFuel) THEN
        lFuel = TRUE; EXIT
      ENDIF
    ENDDO
    nLocalFxr = nLocalFxr + nFxr
  ENDDO
ELSE
  DO ixy = 1, nxy
    ipin = superPin%pin(ixy)
    icel = Pin(ipin)%Cell(iz)
    IF (Cell(icel)%lFuel) THEN
      lFuel = TRUE; EXIT
    ENDIF
  ENDDO
ENDIF

Chi = 0.0

IF (lFuel) THEN
  psisum = zero
  nLocalFxr = 0
  DO ixy = 1, nxy
    ipin = superPin%pin(ixy)
    icel = Pin(ipin)%Cell(iz)
    nFxr = Cell(icel)%nFxr
    FsrIdxSt = Pin(ipin)%FsrIdxSt
    DO i = 1, nFxr
      ifxr = nLocalFxr + i
      IF (.NOT. XsMac(ifxr)%lFuel) CYCLE
      nFsrInFxr = Cell(icel)%nFsrInFxr(i)
      localpsi = 0.0
      DO j = 1, nFsrInFxr
        ireg = Cell(icel)%MapFxr2FsrIdx(j, i)
        ifsr = FsrIdxSt + ireg - 1
        vol = Cell(icel)%vol(ireg)
        DO ig = 1, ng
          psi = phis(ifsr, iz, ig) * vol * XsMac(ifxr)%XsMacNf(ig)
          psisum = psisum + psi
          localpsi = localpsi + psi
        ENDDO
      ENDDO
      DO ig = 1, ng
        Chi(ig) = Chi(ig) + XsMac(ifxr)%Chi(ig) * localpsi
      ENDDO
    ENDDO
    nLocalFxr = nLocalFxr + nFxr
  ENDDO
  Chi = Chi / psisum
ENDIF

PinXS%Chi = Chi

END SUBROUTINE

! ========================================================================================================== !
!  An Optimally Diffusive Coarse Mesh Finite Difference Method to Accelerate Neutron Transport Calculations  !
!         Ang Zhu, Michael Jarret, Yunlin Xu, Brendan Kochunas, Edward Larsen, Thomas Downar (2016)          !
! ========================================================================================================== !

FUNCTION odCMFD(h, xst) RESULT(theta)

IMPLICIT NONE

REAL :: h, xst
REAL :: tau, theta
REAL :: a(0 : 6)
INTEGER :: i

DATA a / -5.542780E-02, 8.740501E-02, -2.152599E-02, 3.145553E-03, -2.683648E-04, 1.222516E-05, -2.284879E-07 /

tau = h * xst

IF (tau .LT. 1.0) THEN
  theta = 0.0; RETURN
ELSEIF (tau .GE. 1.0 .AND. tau .LT. 14.0) THEN
  theta = 0.0
  DO i = 0, 6
    theta = theta + a(i) * tau ** i
  ENDDO
  theta = theta * h; RETURN
ELSE
  theta = 0.127 * h; RETURN
ENDIF

END FUNCTION

SUBROUTINE SetRadialCoupling(Pin, PinXS, Jout, ng, nxy, myzb, myze, lDhat)
USE PARAM
USE geom,           ONLY : ncbd
USE TYPEDEF,        ONLY : PinXS_Type
USE CNTL,           ONLY : nTracerCntl
IMPLICIT NONE

TYPE(superPin_Type), POINTER :: Pin(:)
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: Jout(:, :, :, :, :)
INTEGER :: ng, nxy, myzb, myze
LOGICAL :: lDhat

INTEGER :: ig, ipin, ineighpin, iz, ibd, inbd
REAL :: Dtil, Dhat, myphi, neighphi, mybeta, neighbeta, jnet, jfdm, smy
REAL, POINTER :: superJout(:, :, :, :, :)

IF (nTracerCntl%lHex) THEN
  CALL HexSetRadialCoupling(Pin, PinXS, Jout, ng, nxy, myzb, myze, lDhat)

  RETURN
END IF
ALLOCATE(superJout(3, ncbd, nxy, myzb : myze, ng))

CALL superPinCurrent(Pin, Jout, superJout, ng, nxy, myzb, myze)

!$OMP PARALLEL PRIVATE(ineighpin, inbd, myphi, neighphi, mybeta, neighbeta, Dtil, Dhat, jnet, jfdm, smy)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO iz = myzb, myze
    DO ipin = 1, nxy
      DO ibd = 1, 4
        ineighpin = Pin(ipin)%NeighIdx(ibd)
        inbd = Pin(ipin)%NeighSurfIdx(ibd)
        smy = Pin(ipin)%BdLength(ibd)
        myphi = PinXS(ipin, iz)%Phi(ig)
        mybeta = PinXS(ipin, iz)%XSD(ig) / Pin(ipin)%Center2SurfaceL(ibd)
        IF (ineighpin .GT. 0) THEN
          neighphi = PinXS(ineighpin, iz)%Phi(ig)
          neighbeta = PinXS(ineighpin, iz)%XSD(ig) / Pin(ineighpin)%Center2SurfaceL(inbd)
          Dtil = mybeta * neighbeta / (mybeta + neighbeta) * smy
          jfdm = - Dtil * (neighphi - myphi)
          jnet = superJout(2, ibd, ipin, iz, ig) - superJout(1, ibd, ipin, iz, ig)
          Dhat = - (jnet - jfdm) / (myphi + neighphi)
        ELSE
          IF (ineighpin .EQ. VoidCell) THEN
            neighbeta = 0.5; neighphi = 0.0
          ELSEIF (ineighpin .EQ. RefCell) THEN
            neighbeta = 0.0; neighphi = myphi
          ENDIF
          Dtil = mybeta * neighbeta / (mybeta + neighbeta) * smy
          jfdm = - Dtil * (neighphi - myphi)
          jnet = superJout(2, ibd, ipin, iz, ig) - superJout(1, ibd, ipin, iz, ig)
          Dhat = - (jnet - jfdm) / (myphi + neighphi)
        ENDIF
        PinXS(ipin, iz)%Dtil(ibd, ig) = Dtil
        IF (lDhat) PinXS(ipin, iz)%Dhat(ibd, ig) = Dhat
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

DEALLOCATE(superJout)

END SUBROUTINE
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetRadialCoupling(Pin, PinXS, Jout, ng, nxy, myzb, myze, lDhat)

USE geom,    ONLY : ncbd
USE TYPEDEF, ONLY : PinXS_Type
USE HexCmfd, ONLY : HexSuperPinCurrent
USE CNTL,    ONLY : nTracerCntl

IMPLICIT NONE

TYPE (superPin_Type), POINTER, DIMENSION(:)   :: Pin
TYPE (PinXS_Type),    POINTER, DIMENSION(:,:) :: PinXS

REAL, POINTER, DIMENSION(:,:,:,:,:) :: Jout

INTEGER :: ng, nxy, myzb, myze
LOGICAL :: lDhat
! ----------------------------------------------------
INTEGER :: ig, ipin, ineighpin, iz, iNgh, ibd, jNgh, jbd
REAL :: Dtil, Dhat, myphi, neighphi, mybeta, neighbeta, jnet, jfdm, smy, atil, surfphifdm, ahat

REAL, POINTER, DIMENSION(:,:,:,:,:) :: superJout
! ----------------------------------------------------

ALLOCATE (superJout(3, ncbd, nxy, myzb:myze, ng)) ! # of Ngh is fixed as 15, artibrary #

CALL HexsuperPinCurrent(Pin, Jout, superJout, ng, nxy, myzb, myze)

!$OMP PARALLEL PRIVATE(ineighpin, iNgh, jNgh, ibd, jbd, myphi, neighphi, mybeta, neighbeta, Dtil, Dhat, jnet, jfdm, smy)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO iz = myzb, myze
    DO ipin = 1, nxy
      DO iNgh = 1, Pin(ipin)%nNgh
        ibd       = Pin(ipin)%NghBd(iNgh)
        ineighpin = Pin(ipin)%NeighIdx(iNgh)
        jNgh      = Pin(ipin)%NeighSurfIdx(iNgh)
        smy       = Pin(ipin)%NghLgh(iNgh)
        
        myphi  = PinXS(ipin, iz)%Phi(ig)
        mybeta = PinXS(ipin, iz)%XSD(ig) / Pin(ipin)%Center2SurfaceL(ibd)
        
        IF (ineighpin .GT. 0) THEN
          jbd = Pin(ineighpin)%NghBd(jNgh)
          
          neighphi  = PinXS(ineighpin, iz)%Phi(ig)
          neighbeta = PinXS(ineighpin, iz)%XSD(ig) / Pin(ineighpin)%Center2SurfaceL(jbd)
        ELSE
          IF (ineighpin .EQ. VoidCell) THEN
            neighbeta = 0.5; neighphi = 0.0
          ELSEIF (ineighpin .EQ. RefCell) THEN
            neighbeta = 0.0; neighphi = myphi
          ENDIF
        ENDIF
        
        Dtil = mybeta * neighbeta / (mybeta + neighbeta) * smy
        jfdm = -Dtil * (neighphi - myphi)
        jnet = superJout(2, iNgh, ipin, iz, ig) - superJout(1, iNgh, ipin, iz, ig)
        Dhat = -(jnet - jfdm) / (myphi + neighphi)
        
        PinXS(ipin, iz)%Dtil(iNgh, ig) = Dtil
        
        IF (lDhat) PinXS(ipin, iz)%Dhat(iNgh, ig) = Dhat
        
        ! Dcmp.
        IF (.NOT. nTracerCntl%ldomaindcmp) CYCLE
        
        atil       = mybeta / (mybeta + neighbeta) * smy
        surfphifdm = atil * myphi + (smy - atil) * neighphi
        ahat       = (superJout(3, iNgh, ipin, iz, ig) - surfphifdm) / (myphi + neighphi)
        
        PinXS(ipin, iz)%atil(ingh, ig) = atil ! Multiplied with Surf. Lgh.
        PinXS(ipin, iz)%ahat(ingh, ig) = ahat ! Multiplied with Surf. Lgh.
      END DO
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

DEALLOCATE(superJout)
! ----------------------------------------------------

END SUBROUTINE HexSetRadialCoupling
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE superPinCurrent(Pin, Jout, superJout, ng, nxy, myzb, myze)

IMPLICIT NONE

TYPE(superPin_Type), POINTER :: Pin(:)
REAL, POINTER :: Jout(:, :, :, :, :), superJout(:, :, :, :, :)
INTEGER :: ig, ix, iy, iz, ixy, ipin
INTEGER :: ng, nx, ny, nxy
INTEGER :: myzb, myze

superJout = 0.0

!$OMP PARALLEL PRIVATE(nx, ny, ipin)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO ig = 1, ng
  DO iz = myzb, myze
    DO ixy = 1, nxy
      nx = Pin(ixy)%nx
      ny = Pin(ixy)%ny
      DO ix = 1, nx
        ipin = Pin(ixy)%pin2D(ix, 1)
        superJout(:, NORTH, ixy, iz, ig) = superJout(:, NORTH, ixy, iz, ig) + Jout(:, NORTH, ipin, iz, ig)
        ipin = Pin(ixy)%pin2D(ix, ny)
        superJout(:, SOUTH, ixy, iz, ig) = superJout(:, SOUTH, ixy, iz, ig) + Jout(:, SOUTH, ipin, iz, ig)
      ENDDO
      DO iy = 1, ny
        ipin = Pin(ixy)%pin2D(1, iy)
        superJout(:, WEST, ixy, iz, ig) = superJout(:, WEST, ixy, iz, ig) + Jout(:, WEST, ipin, iz, ig)
        ipin = Pin(ixy)%pin2D(nx, iy)
        superJout(:, EAST, ixy, iz, ig) = superJout(:, EAST, ixy, iz, ig) + Jout(:, EAST, ipin, iz, ig)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE SetSuperpin(CoreInfo, superPin, nxy, myzb, myze, lSuperpin)
USE TYPEDEF,        ONLY : CoreInfo_Type,       Cell_Type,          Pin_Type,                                       &
                           Asy_Type,            AsyInfo_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(superPin_Type), POINTER :: superPin(:)
INTEGER :: nxy, nz, myzb, myze
LOGICAL :: lSuperpin

TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(superPin_Type), POINTER :: mySuperPin(:, :)
INTEGER, POINTER :: node(:, :)
INTEGER :: nx0, ny0, nx, ny, nxAsy(100), nyAsy(100), nxyAsy, pitch
INTEGER :: idx, i, j, ix, iy, ixy, ixa, iya, ixRange(2), iyRange(2), ixbeg, iybeg
INTEGER :: iAsy, iAsyType, ipin, ipin_global, icel, iz
LOGICAL :: lGap

TYPE superPinList_Type
  INTEGER :: nxy, nx, ny
  TYPE(superPin_Type), POINTER :: superPin(:, :)
END TYPE

TYPE(superPinList_Type), POINTER :: superPinList(:)

Cell => CoreInfo%CellInfo
Pin => CoreInfo%Pin
Asy => CoreInfo%Asy
AsyInfo => CoreInfo%AsyInfo

nz = CoreInfo%nz
pitch = nint(sqrt(dble(CoreInfo%nAsyCell)))

IF (.NOT. lSuperpin) THEN
  ALLOCATE(superPin(nxy))
  DO ixy = 1, nxy
    ALLOCATE(superPin(ixy)%pin(1))
    ALLOCATE(superPin(ixy)%pin2D(1, 1))
    ALLOCATE(superPin(ixy)%NeighIdx(4))
    ALLOCATE(superPin(ixy)%NeighSurfIdx(4))
    ALLOCATE(superPin(ixy)%BdLength(4))
    ALLOCATE(superPin(ixy)%Center2SurfaceL(4))
    
    superPin(ixy)%nx = 1
    superPin(ixy)%ny = 1
    superPin(ixy)%nxy = 1
    superPin(ixy)%ix = Pin(ixy)%ix
    superPin(ixy)%iy = Pin(ixy)%iy
    superPin(ixy)%NeighIdx = Pin(ixy)%NeighIdx
    superPin(ixy)%NeighSurfIdx = Pin(ixy)%NeighSurfIdx
    superPin(ixy)%BdLength = Pin(ixy)%BdLength
    superPin(ixy)%Center2SurfaceL = Pin(ixy)%Center2SurfaceL
    superPin(ixy)%Area = Pin(ixy)%BdLength(WEST) * Pin(ixy)%BdLength(SOUTH)
    superPin(ixy)%pin = ixy
    superPin(ixy)%pin2D = ixy
    Pin(ixy)%isuperPin = ixy
  ENDDO
  DO ixy = 1, nxy
    ALLOCATE(superPin(ixy)%lFuel(myzb : myze))
    superPin(ixy)%lFuel = .FALSE.
    DO iz = myzb, myze
      DO i = 1, superPin(ixy)%nxy
        ipin = superPin(ixy)%pin(i)
        icel = Pin(ipin)%Cell(iz)
        superPin(ixy)%iFuelPin = ipin
        IF (Cell(icel)%lFuel) THEN
          superPin(ixy)%lFuel(iz) = .TRUE.
          EXIT
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  RETURN
ENDIF

nxy = 0
DO ipin = 1, CoreInfo%nxy
  lGap = .FALSE.
  DO iz = myzb, myze
    icel = Pin(ipin)%Cell(iz)
    IF (Cell(icel)%lGap) lGap = .TRUE.
  ENDDO
  IF (lGap) CYCLE
  nxy = nxy + 1
ENDDO

ALLOCATE(superPinList(CoreInfo%nxya))

DO iAsy = 1, CoreInfo%nxya
  iAsyType = Asy(iAsy)%AsyType
  IF (Asy(iAsy)%lCentX) THEN
    nxyAsy = pitch * (pitch / 2 + mod(pitch, 2))
    superPinList(iAsy)%nxy = nxyAsy
    superPinList(iAsy)%nx = pitch
    superPinList(iAsy)%ny = pitch / 2 + mod(pitch, 2)
    ALLOCATE(superPinList(iAsy)%superPin(superPinList(iAsy)%nx, superPinList(iAsy)%ny))
    mySuperPin => superPinList(iAsy)%superPin
    DO iy = 1, superPinList(iAsy)%ny
      DO ix = 1, superPinList(iAsy)%nx
        IF (ix .EQ. 1) THEN
          ixRange = (/ 1, 2 /); nx = 2
        ELSEIF (ix .EQ. superPinList(iAsy)%nx) THEN
          ixRange = (/ ix + 1, ix + 2 /); nx = 2
        ELSE
          ixRange = ix + 1; nx = 1
        ENDIF
        IF (iy .EQ. superPinList(iAsy)%ny) THEN
          iyRange = (/ iy, iy + 1 /); ny = 2
        ELSE
          iyRange = iy; ny = 1
        ENDIF
        mySuperPin(ix, iy)%nx = nx
        mySuperPin(ix, iy)%ny = ny
        mySuperPin(ix, iy)%nxy = nx * ny
        ALLOCATE(mySuperPin(ix, iy)%BdLength(4))
        mySuperPin(ix, iy)%BdLength = 0.0
        ALLOCATE(mySuperPin(ix, iy)%pin(nx * ny))
        ALLOCATE(mySuperPin(ix, iy)%pin2D(nx, ny))
        DO j = 1, ny
          DO i = 1, nx
            idx = i + (j - 1) * nx
            ipin = AsyInfo(iAsyType)%Pin2DIdx(ixRange(i), iyRange(j))
            ipin_global = Asy(iAsy)%GlobalPinIdx(ipin)
            mySuperPin(ix, iy)%pin(idx) = ipin_global
            mySuperPin(ix, iy)%pin2D(i, j) = ipin_global
            IF (i .EQ. 1) mySuperPin(ix, iy)%BdLength(WEST) = mySuperPin(ix, iy)%BdLength(WEST)                     &
                                                              + Pin(ipin_global)%BdLength(WEST)
            IF (i .EQ. nx) mySuperPin(ix, iy)%BdLength(EAST) = mySuperPin(ix, iy)%BdLength(EAST)                    &
                                                               + Pin(ipin_global)%BdLength(EAST)
            IF (j .EQ. 1) mySuperPin(ix, iy)%BdLength(NORTH) = mySuperPin(ix, iy)%BdLength(NORTH)                   &
                                                               + Pin(ipin_global)%BdLength(NORTH)
            IF (j .EQ. ny) mySuperPin(ix, iy)%BdLength(SOUTH) = mySuperPin(ix, iy)%BdLength(SOUTH)                  &
                                                                + Pin(ipin_global)%BdLength(SOUTH)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ELSEIF (Asy(iAsy)%lCentY) THEN
    nxyAsy = pitch * (pitch / 2 + mod(pitch, 2))
    superPinList(iAsy)%nxy = nxyAsy
    superPinList(iAsy)%nx = pitch / 2 + mod(pitch, 2)
    superPinList(iAsy)%ny = pitch
    ALLOCATE(superPinList(iAsy)%superPin(superPinList(iAsy)%nx, superPinList(iAsy)%ny))
    mySuperPin => superPinList(iAsy)%superPin
    DO iy = 1, superPinList(iAsy)%ny
      DO ix = 1, superPinList(iAsy)%nx
        IF (ix .EQ. superPinList(iAsy)%nx) THEN
          ixRange = (/ ix, ix + 1 /); nx = 2
        ELSE
          IxRange = ix; nx = 1
        ENDIF
        IF (iy .EQ. 1) THEN
          iyRange = (/ 1, 2 /); ny = 2
        ELSEIF (iy .EQ. superPinList(iAsy)%ny) THEN
          iyRange = (/ iy + 1, iy + 2 /); ny = 2
        ELSE
          iyRange = iy + 1; ny = 1
        ENDIF
        mySuperPin(ix, iy)%nx = nx
        mySuperPin(ix, iy)%ny = ny
        mySuperPin(ix, iy)%nxy = nx * ny
        ALLOCATE(mySuperPin(ix, iy)%BdLength(4))
        mySuperPin(ix, iy)%BdLength = 0.0
        ALLOCATE(mySuperPin(ix, iy)%pin(nx * ny))
        ALLOCATE(mySuperPin(ix, iy)%pin2D(nx, ny))
        DO j = 1, ny
          DO i = 1, nx
            idx = i + (j - 1) * nx
            ipin = AsyInfo(iAsyType)%Pin2DIdx(ixRange(i), iyRange(j))
            ipin_global = Asy(iAsy)%GlobalPinIdx(ipin)
            mySuperPin(ix, iy)%pin(idx) = ipin_global
            mySuperPin(ix, iy)%pin2D(i, j) = ipin_global
            IF (i .EQ. 1) mySuperPin(ix, iy)%BdLength(WEST) = mySuperPin(ix, iy)%BdLength(WEST)                     &
                                                              + Pin(ipin_global)%BdLength(WEST)
            IF (i .EQ. nx) mySuperPin(ix, iy)%BdLength(EAST) = mySuperPin(ix, iy)%BdLength(EAST)                    &
                                                               + Pin(ipin_global)%BdLength(EAST)
            IF (j .EQ. 1) mySuperPin(ix, iy)%BdLength(NORTH) = mySuperPin(ix, iy)%BdLength(NORTH)                   &
                                                               + Pin(ipin_global)%BdLength(NORTH)
            IF (j .EQ. ny) mySuperPin(ix, iy)%BdLength(SOUTH) = mySuperPin(ix, iy)%BdLength(SOUTH)                  &
                                                                + Pin(ipin_global)%BdLength(SOUTH)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ELSEIF (Asy(iAsy)%lCentXY) THEN
    nxyAsy = (pitch / 2 + mod(pitch, 2)) ** 2
    superPinList(iAsy)%nxy = nxyAsy
    superPinList(iAsy)%nx = pitch / 2 + mod(pitch, 2)
    superPinList(iAsy)%ny = pitch / 2 + mod(pitch, 2)
    ALLOCATE(superPinList(iAsy)%superPin(superPinList(iAsy)%nx, superPinList(iAsy)%ny))
    mySuperPin => superPinList(iAsy)%superPin
    DO iy = 1, superPinList(iAsy)%ny
      DO ix = 1, superPinList(iAsy)%nx
        IF (ix .EQ. superPinList(iAsy)%nx) THEN
          ixRange = (/ ix, ix + 1 /); nx = 2
        ELSE
          ixRange = ix; nx = 1
        ENDIF
        IF (iy .EQ. superPinList(iAsy)%ny) THEN
          iyRange = (/ iy, iy + 1 /); ny = 2
        ELSE
          iyRange = iy; ny = 1
        ENDIF
        mySuperPin(ix, iy)%nx = nx
        mySuperPin(ix, iy)%ny = ny
        mySuperPin(ix, iy)%nxy = nx * ny
        ALLOCATE(mySuperPin(ix, iy)%BdLength(4))
        mySuperPin(ix, iy)%BdLength = 0.0
        ALLOCATE(mySuperPin(ix, iy)%pin(nx * ny))
        ALLOCATE(mySuperPin(ix, iy)%pin2D(nx, ny))
        DO j = 1, ny
          DO i = 1, nx
            idx = i + (j - 1) * nx
            ipin = AsyInfo(iAsyType)%Pin2DIdx(ixRange(i), iyRange(j))
            ipin_global = Asy(iAsy)%GlobalPinIdx(ipin)
            mySuperPin(ix, iy)%pin(idx) = ipin_global
            mySuperPin(ix, iy)%pin2D(i, j) = ipin_global
            IF (i .EQ. 1) mySuperPin(ix, iy)%BdLength(WEST) = mySuperPin(ix, iy)%BdLength(WEST)                     &
                                                              + Pin(ipin_global)%BdLength(WEST)
            IF (i .EQ. nx) mySuperPin(ix, iy)%BdLength(EAST) = mySuperPin(ix, iy)%BdLength(EAST)                    &
                                                               + Pin(ipin_global)%BdLength(EAST)
            IF (j .EQ. 1) mySuperPin(ix, iy)%BdLength(NORTH) = mySuperPin(ix, iy)%BdLength(NORTH)                   &
                                                               + Pin(ipin_global)%BdLength(NORTH)
            IF (j .EQ. ny) mySuperPin(ix, iy)%BdLength(SOUTH) = mySuperPin(ix, iy)%BdLength(SOUTH)                  &
                                                                + Pin(ipin_global)%BdLength(SOUTH)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ELSE
    nxyAsy = pitch ** 2
    superPinList(iAsy)%nxy = nxyAsy
    superPinList(iAsy)%nx = pitch
    superPinList(iAsy)%ny = pitch
    ALLOCATE(superPinList(iAsy)%superPin(superPinList(iAsy)%nx, superPinList(iAsy)%ny))
    mySuperPin => superPinList(iAsy)%superPin
    DO iy = 1, superPinList(iAsy)%ny
      DO ix = 1, superPinList(iAsy)%nx
        IF (ix .EQ. 1) THEN
          ixRange = (/ 1, 2 /); nx = 2
        ELSEIF (ix .EQ. superPinList(iAsy)%nx) THEN
          ixRange = (/ ix + 1, ix + 2 /); nx = 2
        ELSE
          ixRange = ix + 1; nx = 1
        ENDIF
        IF (iy .EQ. 1) THEN
          iyRange = (/ 1, 2 /); ny = 2
        ELSEIF (iy .EQ. superPinList(iAsy)%ny) THEN
          iyRange = (/ iy + 1, iy + 2 /); ny = 2
        ELSE
          iyRange = iy + 1; ny = 1
        ENDIF
        mySuperPin(ix, iy)%nx = nx
        mySuperPin(ix, iy)%ny = ny
        mySuperPin(ix, iy)%nxy = nx * ny
        ALLOCATE(mySuperPin(ix, iy)%BdLength(4))
        mySuperPin(ix, iy)%BdLength = 0.0
        ALLOCATE(mySuperPin(ix, iy)%pin(nx * ny))
        ALLOCATE(mySuperPin(ix, iy)%pin2D(nx, ny))
        DO j = 1, ny
          DO i = 1, nx
            idx = i + (j - 1) * nx
            ipin = AsyInfo(iAsyType)%Pin2DIdx(ixRange(i), iyRange(j))
            ipin_global = Asy(iAsy)%GlobalPinIdx(ipin)
            mySuperPin(ix, iy)%pin(idx) = ipin_global
            mySuperPin(ix, iy)%pin2D(i, j) = ipin_global
            IF (i .EQ. 1) mySuperPin(ix, iy)%BdLength(WEST) = mySuperPin(ix, iy)%BdLength(WEST)                     &
                                                              + Pin(ipin_global)%BdLength(WEST)
            IF (i .EQ. nx) mySuperPin(ix, iy)%BdLength(EAST) = mySuperPin(ix, iy)%BdLength(EAST)                    &
                                                               + Pin(ipin_global)%BdLength(EAST)
            IF (j .EQ. 1) mySuperPin(ix, iy)%BdLength(NORTH) = mySuperPin(ix, iy)%BdLength(NORTH)                   &
                                                               + Pin(ipin_global)%BdLength(NORTH)
            IF (j .EQ. ny) mySuperPin(ix, iy)%BdLength(SOUTH) = mySuperPin(ix, iy)%BdLength(SOUTH)                  &
                                                                + Pin(ipin_global)%BdLength(SOUTH)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDIF
ENDDO

nxAsy = 0
DO iya = 1, CoreInfo%nya
  DO ixa = 1, CoreInfo%nxa
    iAsy = CoreInfo%CoreIdx(ixa, iya)
    IF (iAsy .EQ. 0) CYCLE
    nxAsy(ixa) = max(nxAsy(ixa), superPinList(iAsy)%nx)
  ENDDO
ENDDO

nyAsy = 0
DO ixa = 1, CoreInfo%nya
  DO iya = 1, CoreInfo%nxa
    iAsy = CoreInfo%CoreIdx(ixa, iya)
    IF (iAsy .EQ. 0) CYCLE
    nyAsy(iya) = max(nxAsy(iya), superPinList(iAsy)%ny)
  ENDDO
ENDDO

DO iya = 1, CoreInfo%nya
  nx0 = 0
  DO ixa = 1, CoreInfo%nxa
    iAsy = CoreInfo%CoreIdx(ixa, iya)
    IF (iAsy .EQ. 0) CYCLE
    nx0 = nx0 + superPinList(iAsy)%nx
  ENDDO
  nx = max(nx0, nx)
ENDDO

DO ixa = 1, CoreInfo%nxa
  ny0 = 0
  DO iya = 1, CoreInfo%nya
    iAsy = CoreInfo%CoreIdx(ixa, iya)
    IF (iAsy .EQ. 0) CYCLE
    ny0 = ny0 + superPinList(iAsy)%ny
  ENDDO
  ny = max(ny0, ny)
ENDDO

ALLOCATE(node(0 : nx + 1, 0 : ny + 1))

node = VoidCell
IF (CoreInfo%RadBC(SOUTH) .EQ. RefCell) node(1 : nx, ny + 1) = RefCell
IF (CoreInfo%RadBC(WEST) .EQ. RefCell) node(0, 1 : ny) = RefCell
IF (CoreInfo%RadBC(NORTH) .EQ. RefCell) node(1 : nx, 0) = RefCell
IF (CoreInfo%RadBC(EAST) .EQ. RefCell) node(nx + 1, 1 : ny) = RefCell

ALLOCATE(superPin(nxy))

iybeg = 0; ixy = 0
DO iya = 1, CoreInfo%nya
  ixbeg = 0
  DO ixa = 1, CoreInfo%nxa
    iAsy = CoreInfo%CoreIdx(ixa, iya)
    IF (iAsy .EQ. 0) THEN
      ixbeg = ixbeg + nxAsy(ixa); CYCLE
    ENDIF
    mySuperPin => superPinList(iAsy)%superPin
    DO iy = 1, superPinList(iAsy)%ny
      DO ix = 1, superPinList(iAsy)%nx
        ixy = ixy + 1
        node(ixbeg + ix, iybeg + iy) = ixy
        superPin(ixy)%ix = ixbeg + ix
        superPin(ixy)%iy = iybeg + iy
        superPin(ixy)%nx = mySuperPin(ix, iy)%nx
        superPin(ixy)%ny = mySuperPin(ix, iy)%ny
        superPin(ixy)%nxy = mySuperPin(ix, iy)%nxy
        ALLOCATE(superPin(ixy)%pin(mySuperPin(ix, iy)%nxy))
        ALLOCATE(superPin(ixy)%pin2D(mySuperPin(ix, iy)%nx, mySuperPin(ix, iy)%ny))
        ALLOCATE(superPin(ixy)%BdLength(4))
        superPin(ixy)%pin = mySuperPin(ix, iy)%pin
        superPin(ixy)%pin2D = mySuperPin(ix, iy)%pin2D
        superPin(ixy)%BdLength = mySuperPin(ix, iy)%BdLength
      ENDDO
    ENDDO
    ixbeg = ixbeg + superPinList(iAsy)%nx
  ENDDO
  iybeg = iybeg + nyAsy(iya)
ENDDO

IF (CoreInfo%lRot) THEN
  node(0, 1 : nx) = node(1 : nx, 1)
  node(1 : ny, 0) = node(1, 1 : ny)
ELSEIF (CoreInfo%lCbd) THEN

ENDIF

DO iy = 1, ny
  DO ix = 1, nx
    ixy = node(ix, iy)
    IF (ixy .EQ. 0) CYCLE
    ALLOCATE(superPin(ixy)%NeighIdx(4))
    ALLOCATE(superPin(ixy)%NeighSurfIdx(4))
    ALLOCATE(superPin(ixy)%Center2SurfaceL(4))
    superPin(ixy)%NeighIdx(SOUTH) = node(ix, iy + 1)
    superPin(ixy)%NeighIdx(WEST) = node(ix - 1, iy)
    superPin(ixy)%NeighIdx(NORTH) = node(ix, iy - 1)
    superPin(ixy)%NeighIdx(EAST) = node(ix + 1, iy)
    superPin(ixy)%NeighSurfIdx(SOUTH) = NORTH
    superPin(ixy)%NeighSurfIdx(WEST) = EAST
    superPin(ixy)%NeighSurfIdx(NORTH) = SOUTH
    superPin(ixy)%NeighSurfIdx(EAST) = WEST
    superPin(ixy)%Center2SurfaceL(SOUTH) = superPin(ixy)%BdLength(WEST) * 0.5
    superPin(ixy)%Center2SurfaceL(WEST) = superPin(ixy)%BdLength(SOUTH) * 0.5
    superPin(ixy)%Center2SurfaceL(NORTH) = superPin(ixy)%BdLength(EAST) * 0.5
    superPin(ixy)%Center2SurfaceL(EAST) = superPin(ixy)%BdLength(NORTH) * 0.5
    superPin(ixy)%Area = superPin(ixy)%BdLength(WEST) * superPin(ixy)%BdLength(SOUTH)
  ENDDO
ENDDO

IF (CoreInfo%lRot) THEN
  DO ix = 1, nx
    ixy = node(ix, 1)
    superPin(ixy)%NeighSurfIdx(NORTH) = WEST
  ENDDO
  DO iy = 1, ny
    ixy = node(1, iy)
    superPin(ixy)%NeighSurfIdx(WEST) = NORTH
  ENDDO
ELSEIF (CoreInfo%lCbd) THEN

ENDIF

DO ixy = 1, nxy
  ALLOCATE(superPin(ixy)%lFuel(nz))
  superPin(ixy)%lFuel = .FALSE.
  DO iz = 1, nz
    DO i = 1, superPin(ixy)%nxy
      ipin = superPin(ixy)%pin(i)
      icel = Pin(ipin)%Cell(iz)
      IF (Cell(icel)%lFuel) THEN
        superPin(ixy)%lFuel(iz) = .TRUE.
        EXIT
      ENDIF
    ENDDO
  ENDDO
  DO i = 1, superPin(ixy)%nxy
    ipin = superPin(ixy)%pin(i)
    IF (Pin(ipin)%lFuel) THEN
      superPin(ixy)%iFuelPin = ipin
      EXIT
    ENDIF
  ENDDO
ENDDO

DO ixy = 1, nxy
  DO i = 1, superPin(ixy)%nxy
    ipin = superPin(ixy)%pin(i)
    Pin(ipin)%isuperPin = ixy
  END DO
END DO

END SUBROUTINE

SUBROUTINE HomogenizeKinParam(Core, FmInfo, TranInfo, TranCntl, GroupInfo, superPin, PinXS, nxy, myzb, myze, ng, nprec, lxslib)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,     FmInfo_Type,      TranInfo_Type,      TranCntl_Type,    GroupInfo_Type,   &
                             PinXS_Type,        FxrInfo_Type,     Pin_Type,           Cell_Type
USE BenchXs,          ONLY : DnpBetaBen,        NeutVeloBen,      DnpBetaDynBen,      NeutVeloDynBen
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(superPin_Type), POINTER :: superPin(:)
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
INTEGER :: nxy, myzb, myze, ng, nprec
LOGICAL :: lxslib

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(FxrInfo_Type), POINTER :: myFxr
REAL, POINTER :: Psi(:, :), phis(:, :, :)
REAL :: beta0(nprec), beta(nprec), chid(ng)
REAL :: velo0(ng), velo(ng), phisum(ng)
REAL :: vol, psisum, locpsi, locphi, betat
INTEGER :: norg, nChi, irgb, irge
INTEGER :: FsrIdxSt, FxrIdxSt, nFsrInFxr
INTEGER :: iz, ixy, ipin, icel, ifxr, ifsrlocal, ifsr, ig
INTEGER :: i, j, k

Pin => Core%Pin
Cell => Core%CellInfo
Fxr => FmInfo%Fxr
Psi => FmInfo%Psi
phis => FmInfo%Phis

nChi = GroupInfo%nChi
IF(lxslib) THEN
  norg = GroupInfo%norg
  irgb = GroupInfo%nofg + 1
  irge = GroupInfo%nofg + GroupInfo%norg
END IF

chid(1:ng) = TranInfo%chid(1:ng)

!$OMP PARALLEL PRIVATE(beta, velo, psisum, phisum, ipin, icel, FsrIdxSt, FxrIdxSt, ifxr, nFsrInFxr, &
!$OMP                  myFxr, beta0, velo0, ifsrlocal, ifsr, vol, locpsi, locphi, betat)
DO iz = myzb, myze
  !$OMP DO SCHEDULE(GUIDED)
  DO ixy = 1, nxy
    beta = 0.
    velo = 0.
    psisum = 0.
    phisum = 0.
    DO k = 1, superPin(ixy)%nxy
      ipin = superPin(ixy)%pin(k)
      icel = Pin(ipin)%Cell(iz)
      FsrIdxSt = Pin(ipin)%FsrIdxSt
      FxrIdxSt = Pin(ipin)%FxrIdxSt

      DO j = 1, Cell(icel)%nFxr
        ifxr = FxrIdxSt + j - 1
        nFsrInFxr = Cell(icel)%nFsrInFxr(j)

        myFxr => Fxr(ifxr, iz)
        IF(lxslib) THEN
          beta0(1:nprec) = myFxr%Beta(1:nprec)
          velo0(1:ng) = myFxr%veloh(1:ng)
        ELSE
          IF(TranCntl%lDynamicBen) THEN
            CALL DnpBetaDynBen(myFxr%imix, TranInfo%fuelTemp(ipin, iz), beta0(1:nprec))
            CALL NeutVeloDynBen(myFxr%imix, TranInfo%fuelTemp(ipin, iz), velo0(1:ng))
          ELSE
            CALL DnpBetaBen(myFxr%imix, beta0(1:nprec))
            CALL NeutVeloBen(myFxr%imix, velo0(1:ng))
          END IF
        END IF

        DO i = 1, nFsrInFxr
          ifsrlocal = Cell(icel)%MapFxr2FsrIdx(i, j)
          ifsr = FsrIdxSt + ifsrlocal - 1


          vol = Cell(icel)%vol(ifsrlocal)
          locpsi = Psi(ifsr, iz) * vol
          beta(1:nprec) = beta(1:nprec) + beta0(1:nprec) * locpsi
          psisum = psisum + locpsi
          DO ig = 1, ng
            locphi = phis(ifsr, iz, ig) * vol
            velo(ig) = velo(ig) + 1._8 / velo0(ig) * locphi
            phisum(ig) = phisum(ig) + locphi
          END DO
        END DO
      END DO
    END DO
    IF(psisum .GT. 0) THEN
      beta(1:nprec) = beta(1:nprec) / psisum
    ELSE
      beta(1:nprec) = 0.
    END IF
    PinXS(ixy, iz)%Beta(1:nprec) = beta(1:nprec)
    DO ig = 1, ng
      velo(ig) = phisum(ig) / velo(ig)
      PinXS(ixy, iz)%velo(ig) = velo(ig)
    END DO

    betat = sum(beta(1:nprec))
    PinXS(ixy, iz)%betat = betat
    PinXS(ixy, iz)%chip(1:ng) = (PinXS(ixy, iz)%Chi(1:ng) - betat * chid(1:ng)) / (1._8 - betat)
  END DO
  !$OMP END DO
END DO
!$OMP END PARALLEL

NULLIFY(Pin, Cell, Fxr, Psi, phis)

END SUBROUTINE

SUBROUTINE SetCMFDPrecCoeff(TranInfo, TranCntl, PinXS, CellOmegam, CellOmega0, CellOmegap, nxy, myzb, myze, nprec)
USE PARAM
USE TYPEDEF,          ONLY : TranInfo_Type,      TranCntl_Type,      PinXS_Type
IMPLICIT NONE
TYPE(TranInfo_Type) :: TranInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: CellOmegam(:, :, :), CellOmega0(:, :, :), CellOmegap(:, :, :)
INTEGER :: nxy, myzb, myze, nprec

REAL :: lambda(nprec), invlambda(nprec)
REAL :: invldt(nprec), invldtgp1(nprec), kapbinvldt(nprec), kapbinvldt2(nprec), kappa(nprec), kappap1(nprec)
REAL :: omegam(nprec), omega0(nprec), omegap(nprec)
REAL :: delt, deltp, gamma, invgamma, invgammap1
REAL :: omegalm, omegal0, omegalp
INTEGER :: norder, nowstep
INTEGER :: i, ixy, iz

nowstep = TranCntl%nowstep
delt = TranCntl%delt(nowstep)
IF(nowstep .EQ. 1) THEN ! .OR. abs(TranCntl%theta - 0.5) .GT. epsm6) THEN
  norder = 1
  deltp = TranCntl%Delt(nowstep)
ELSE
  norder = 2
  deltp = TranCntl%Delt(nowstep-1)
END IF

IF(norder .EQ. 2) THEN
  gamma = delt / deltp
  invgamma = 1._8 / gamma
  invgammap1 = 1._8 / (gamma + 1._8)
END IF
DO i = 1, nprec
  lambda(i) = TranInfo%lambda(i)
  invlambda(i) = 1._8 / TranInfo%lambda(i)
  kappa(i) = exp(-lambda(i) * delt)
  kappap1(i) = kappa(i) + 1._8
  invldt(i) = 1._8 / (deltp * lambda(i))
  kapbinvldt(i) = (1._8 - kappa(i)) * invldt(i)
  IF(norder .EQ. 2) THEN
    invldtgp1(i) = invldt(i) * invgammap1
    kapbinvldt2(i) = (1._8 - kappa(i)) * (1._8 - 2._8 * invldt(i))
  END IF
END DO
IF(norder .EQ. 1) THEN
  DO i = 1, nprec
    omegam(i) = 0
    omega0(i) = invlambda(i) * (kapbinvldt(i) - kappa(i))
    omegap(i) = invlambda(i) * (1 - kapbinvldt(i))
  END DO
ELSE
  DO i = 1, nprec
    omegam(i) = invlambda(i) * invldtgp1(i) * (2._8 * kapbinvldt(i) - gamma * kappap1(i))
    omega0(i) = invlambda(i) * (invldt(i) *(kappap1(i) + kapbinvldt2(i) * invgamma) - kappa(i))
    omegap(i) = invlambda(i) * (1 - invldtgp1(i) * (2._8 + kapbinvldt2(i) * invgamma))
  END DO
END IF

CellOmegam = 0._8
CellOmega0 = 0._8
CellOmegap = 0._8

!$OMP PARALLEL PRIVATE(omegalm, omegal0, omegalp)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(2)
DO iz = myzb, myze
  DO ixy = 1, nxy
    omegalm = 0.
    omegal0 = 0.
    omegalp = 0.
    DO i = 1, nprec
      CellOmegam(i, ixy, iz) = PinXS(ixy, iz)%beta(i) * omegam(i)
      CellOmega0(i, ixy, iz) = PinXS(ixy, iz)%beta(i) * omega0(i)
      CellOmegap(i, ixy, iz) = PinXS(ixy, iz)%beta(i) * omegap(i)

      omegalm = omegalm + CellOmegam(i, ixy, iz) * lambda(i)
      omegal0 = omegal0 + CellOmega0(i, ixy, iz) * lambda(i)
      omegalp = omegalp + CellOmegap(i, ixy, iz) * lambda(i)
    END DO
    CellOmegam(0, ixy, iz) = omegalm
    CellOmega0(0, ixy, iz) = omegal0
    CellOmegap(0, ixy, iz) = omegalp

    PinXS(ixy, iz)%omega = omegalp
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

END MODULE
