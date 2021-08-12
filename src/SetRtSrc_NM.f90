#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtsrcNM(Core, Fxr, srcNg, phisNg, psi, AxSrc, xstNg, eigv, iz,gb, ge, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix, PE, Offset)

USE OMP_LIB
USE PARAM,        ONLY : ZERO, ONE, TRUE, FALSE, nThreadMax
USE TYPEDEF,      ONLY : coreinfo_type, Fxrinfo_type, Cell_Type, pin_Type, GroupInfo_Type, PE_TYPE, XsMac_Type
USE BenchXs,      ONLY : GetChiBen, GetChiDynBen, xssben, xssDynBen
USE MacXsLib_mod, ONLY : MacXsScatMatrix
USE CNTL,         ONLY : nTracerCntl
USE TRAN_MOD,     ONLY : TranInfo

IMPLICIT NONE

TYPE (CoreInfo_Type)  :: Core
TYPE (GroupInfo_Type) :: GroupInfo
TYPE (PE_TYPE)        :: PE

TYPE (FxrInfo_Type), DIMENSION(:) :: Fxr

REAL, POINTER, DIMENSION(:,:)   :: srcNg, psi, xstNg, phisNg
REAL, POINTER, DIMENSION(:,:,:) :: AxSrc

REAL :: eigv

INTEGER :: gb, ge, ng, iz, ifsr, ifxr
LOGICAL :: lxslib, lscat1, l3dim, lNegFix
INTEGER, OPTIONAL :: Offset
! ----------------------------------------------------
TYPE (Pin_Type),  POINTER, DIMENSION(:) :: Pin
TYPE (Cell_Type), POINTER, DIMENSION(:) :: CellInfo

TYPE (XsMac_Type), SAVE, DIMENSION(nThreadMax) :: XsMac

INTEGER :: nCoreFsr, nCoreFxr, nxy, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, nchi, xyb, xye, off
INTEGER :: ipin, icel, ifsrlocal, itype, ig, ig2, tid
INTEGER :: i, j, k

LOGICAL :: lNegSrcFix

REAL :: reigv

REAL, DIMENSION(ng) :: chi
REAL, POINTER, DIMENSION(:,:) :: xsmacsm
! ----------------------------------------------------

Pin      => Core%Pin
CellInfo => Core%CellInfo
nCoreFsr  = Core%nCoreFsr
nCoreFxr  = Core%nCoreFxr
nxy       = Core%nxy

xyb = PE%myPinBeg ! Domain Decomposition + MPI
xye = PE%myPinEnd

off = 0
IF (PRESENT(Offset)) off = Offset

reigv = one / eigv

lNegSrcFix = FALSE
IF (lNegFix) lNegSrcFix = TRUE

srcNg = zero

IF (lxsLib) THEN
  nchi = GroupInfo%nchi
ELSE
  DO i = 1, PE%nThread
    ALLOCATE (XsMac(i)%XsMacSm (ng, ng))
  END DO
END IF
! ----------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i, j, k, ifsr, ifxr, ipin, icel, ifsrlocal, itype, ig, ig2, tid, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, xsmacsm, chi)
tid = omp_get_thread_num() + 1
!$OMP DO SCHEDULE(GUIDED)
DO ipin = xyb, xye
  FsrIdxSt  = Pin(ipin)%FsrIdxSt
  FxrIdxSt  = Pin(ipin)%FxrIdxSt
  icel      = Pin(ipin)%Cell(iz)
  nlocalFxr = CellInfo(icel)%nFxr
  ! --------------------------------------------------
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    
    IF (lXsLib) Then
      DO ig = gb, ge
        CHI(ig) = 0
        
        IF (ig.LE.nchi .AND. Fxr(ifxr)%ldepl) CHI(ig) = Fxr(ifxr)%chi(ig)
      END DO
    ELSE
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1, j)
      itype     = Fxr(ifxr)%imix
      
      IF (nTracerCntl%lDynamicBen) THEN
        CHI(gb:ge) = GetChiDynBen(itype, TranInfo%fuelTemp(ipin, iz), gb, ge)
      ELSE
        CHI(gb:ge) = GetChiBen(itype, gb, ge)
      END IF
    END IF
    
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      
      DO ig = gb, ge
        srcNg(ig - off, ifsr) = reigv * chi(ig) * psi(ifsr, iz)
      END DO
    END DO
  END DO
  ! --------------------------------------------------
  DO j = 1, nLocalFxr
    ifxr      = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    
    IF (lXsLib) Then
      CALL MacXsScatMatrix(XsMac(tid), Fxr(ifxr), gb, ge, ng, GroupInfo, lscat1, TRUE)
      
      XsMacSm => XsMac(tid)%XsMacSm
#ifdef inflow
      DO ig = gb, ge
        XsMacSm(ig, ig) = XsMacSm(ig, ig) + Fxr(ifxr)%DelInflow(ig)
      END DO
#endif
    ELSE
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1,j)
      itype     = Fxr(ifxr)%imix
      XsMacSm  => XsMac(tid)%XsMacSm
      
      DO ig = gb, ge
        IF (nTracerCntl%lDynamicBen) THEN
          CALL xssDynben(itype, TranInfo%fuelTemp(ipin, iz), ig, 1, ng, XsMacsm, lscat1)
        ELSE
          CALL xssben(itype, ig, 1, ng, XsMacsm, lscat1)
        END IF
      END DO
    END IF
    
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      
      DO ig = gb, ge
        DO ig2 = GroupInfo%InScatRange(1, ig), GroupInfo%InScatRange(2, ig)
          srcNg(ig - off, ifsr) = srcNg(ig - off, ifsr) + xsmacsm(ig2, ig) * phisNg(ig2, ifsr)
        END DO
        
        IF (lNegSrcFix .AND. srcNg(ig - off, ifsr).LT.ZERO) THEN
          srcNg(ig- off, ifsr) = srcNg(ig - off, ifsr) - xsmacsm(ig, ig) * phisNg(ig, ifsr)
          xstNg(ig,      ifsr) = xstNg(ig,       ifsr) - xsmacsm(ig, ig)
        END IF
      END DO
    END DO
  END DO
  ! --------------------------------------------------
  IF (l3dim) THEN
    IF (nTracerCntl%LkgSplitLv .EQ. 0) THEN
      DO j = 1, nLocalFxr
        ifxr      = FxrIdxSt + j -1
        nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
        
        DO i = 1, nFsrInFxr
          ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
          
          DO ig = gb, ge
            IF (AxSrc(ipin, iz, ig).LT.ZERO .AND. .NOT.Fxr(ifxr)%lvoid) srcNg(ig - off, ifsr) = srcNg(ig - off, ifsr) - AxSrc(ipin, iz, ig)
          END DO
        END DO
      END DO
    ELSE
      DO j = 1, nLocalFxr
        ifxr      = FxrIdxSt + j -1
        nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
        
        DO i = 1, nFsrInFxr
          ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
          
          DO ig = gb, ge
            srcNg(ig - off, ifsr) = srcNg(ig - off, ifsr) - AxSrc(ipin, iz, ig)
          END DO
        END DO
      END DO
    END IF
  END IF
  ! --------------------------------------------------
  DO j = 1, CellInfo(icel)%nFsr
    ifsr = FsrIdxSt + j - 1
    
    srcNg(gb - off:ge - off, ifsr) = srcNg(gb - off:ge - off, ifsr) / xstNg(gb:ge, ifsr)
  END DO
  
  NULLIFY (Xsmacsm)
END DO
!$OMP END DO
!$OMP END PARALLEL
! --------------------------------------------------
IF (.NOT. lxsLib) THEN
  DO i = 1, PE%nThread
    DEALLOCATE (XsMac(i)%xsmacsm)
  END DO
END IF

NULLIFY (Pin)
NULLIFY (CellInfo)
! --------------------------------------------------

END SUBROUTINE SetRtsrcNM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtsrcNM_Cusping(Core, FmInfo, Fxr, srcNg, phisNg, psi, AxSrc, xstNg, eigv, iz, gb, ge, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix, PE, Offset)

USE PARAM
USE TYPEDEF,        ONLY : coreinfo_type,          Fxrinfo_type,          Cell_Type,      &
                           pin_Type,               GroupInfo_Type,        PE_TYPE,        &
                           XsMac_Type,             FmInfo_Type
USE BenchXs,        ONLY : GetChiBen,              xssben,                GetChiDynBen,     xssben_Cusping, &
                           MacXsBen,               xssDynBen,             xssDynBen_cusping, DynMacXsBen
USE MacXsLib_mod,   ONLY : MacXsScatMatrix
USE CNTL,           ONLY : nTracerCntl
USE OMP_LIB
USE TRAN_MOD,       ONLY : TranInfo
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(FxrInfo_Type) :: Fxr(:)
TYPE(GroupInfo_Type):: GroupInfo
TYPE(PE_TYPE) :: PE
REAL, POINTER :: srcNg(:, :), phisNg(:, :), psi(:, :), AxSrc(:, :, :), xstNg(:, :)
REAL :: eigv
INTEGER :: gb, ge, ng, iz, ifsr, ifxr
LOGICAL :: lxslib, lscat1, l3dim
LOGICAL :: lNegFix
INTEGER, OPTIONAL :: Offset   !--- CNJ Edit : GPU Acceleration

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(XsMac_Type), SAVE :: XsMac(nThreadMax)
INTEGER :: nCoreFsr, nCoreFxr, nxy, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, nchi
INTEGER :: ipin, icel, ifsrlocal, itype, ig, ig2, tid
INTEGER :: i, j, k
INTEGER :: xyb, xye   !--- CNJ Edit : Domain Decomposition + MPI
INTEGER :: off        !--- CNJ Edit : GPU Acceleration

REAL :: phiz(ng), philz(ng), phiuz(ng)
REAL :: vol, volsum
REAL :: reigv
REAL ::  chi(ng)
REAL, POINTER :: xsmacsm(:,:)

LOGICAL :: lNegSrcFix
LOGICAL :: lCusping
Pin => Core%Pin
CellInfo => Core%CellInfo
nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr
nxy = Core%nxy
xyb = PE%myPinBeg; xye = PE%myPinEnd   !--- CNJ Edit : Domain Decomposition + MPI
off = 0
IF (PRESENT(Offset)) off = Offset
reigv = one/eigv

lNegSrcFix = FALSE
IF(lNegFix) lNegSrcFix = TRUE

srcNg = zero
IF(lxsLib) nchi = GroupInfo%nchi

IF(.NOT. lxsLib) THEN
  DO i = 1, PE%nThread
    ALLOCATE(XsMac(i)%XsMacSm(ng, ng))
  ENDDO
ENDIF

!$OMP PARALLEL DEFAULT(SHARED)                                                                                      &
!$OMP PRIVATE(i, j, k, ifsr, ifxr, ipin, icel, ifsrlocal, itype, ig, ig2, tid,                                      &
!$OMP         FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, xsmacsm, chi, phiz, philz, phiuz, vol, volsum)
tid = omp_get_thread_num() + 1
!$OMP DO SCHEDULE(GUIDED)
DO ipin = xyb, xye
  FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    IF(lXsLib) Then
      DO ig = gb, ge
        IF(ig .gt. nchi) THEN
          CHI(ig) = 0
        ELSE
          CHI(ig) = 0
          IF(Fxr(ifxr)%ldepl) CHI(ig) = Fxr(ifxr)%chi(ig)
        ENDIF
      ENDDO
    ELSE
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1, j)
      itype = Fxr(ifxr)%imix
      IF(nTracerCntl%lDynamicBen) THEN
        CHI(gb : ge) = GetChiDynBen(itype, TranInfo%fuelTemp(ipin, iz), gb, ge)
      ELSE
        CHI(gb : ge) = GetChiBen(itype, gb, ge)
      END IF
    ENDIF
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      DO ig = gb, ge
        srcNg(ig - off, ifsr) = reigv * chi(ig) * psi(ifsr, iz)
      ENDDO
    ENDDO
  ENDDO
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    IF(lXsLib) Then
      CALL MacXsScatMatrix(XsMac(tid), Fxr(ifxr), gb, ge, ng, GroupInfo, lscat1, .TRUE.)
      XsMacSm => XsMac(tid)%XsMacSm
#ifdef inflow
      DO ig = gb, ge
        XsMacSm(ig, ig) = XsMacSm(ig, ig) + Fxr(ifxr)%DelInflow(ig)
      ENDDO
#endif
    ELSE
      itype = Fxr(ifxr)%imix
      XsMacSm => XsMac(tid)%XsMacSm
      IF(nTracerCntl%lDynamicBen) THEN
        lCusping = DynMacXsBen(itype)%lCusping
      ELSE
        lCusping = MacXsBen(itype)%lCusping
      END IF
      IF(lCusping) THEN
        phiz = 0.
        philz = 0.
        phiuz = 0.
        volsum = 0.
        DO ig = 1, ng
          DO i = 1, nFsrInFxr
            ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(i,j)
            ifsr = FsrIdxSt + ifsrlocal - 1
            vol = CellInfo(icel)%vol(ifsrlocal)
            IF(ig .EQ. 1) volsum = volsum + vol
            phiz(ig) = phiz(ig) + FmInfo%phis(ifsr, iz, ig) * vol
            IF(iz .EQ. PE%myzb) THEN 
              philz(ig) = philz(ig) + FmInfo%neighphis(ifsr, ig, BOTTOM) * vol
            ELSE
              philz(ig) = philz(ig) + FmInfo%phis(ifsr, iz-1, ig) * vol
            END IF
            IF(iz .EQ. PE%myze) THEN
              phiuz(ig) = phiuz(ig) + FmInfo%neighphis(ifsr, ig, TOP) * vol
            ELSE
              phiuz(ig) = phiuz(ig) + FmInfo%phis(ifsr, iz+1, ig) * vol
            END IF
          END DO
          phiz(ig) = phiz(ig) / volsum
          philz(ig) = philz(ig) / volsum
          phiuz(ig) = phiuz(ig) / volsum
        END DO
        DO ig = gb, ge
          IF(nTracerCntl%lDynamicBen) THEN
            CALL xssDynben_Cusping(itype, TranInfo%fuelTemp(ipin, iz), ig, 1, ng, XsMacsm, lscat1, phiz, philz, phiuz, &
              Core%hzfm(iz), Core%hzfm(iz-1), Core%hzfm(iz+1))
          ELSE
            CALL xssben_Cusping(itype, ig, 1, ng, XsMacsm, lscat1, phiz, philz, phiuz, &
              Core%hzfm(iz), Core%hzfm(iz-1), Core%hzfm(iz+1))
          END IF
        END DO
      ELSE
        DO ig = gb, ge
          IF(nTracerCntl%lDynamicBen) THEN
            CALL xssDynben(itype, TranInfo%fuelTemp(ipin, iz), ig, 1, ng, XsMacsm, lscat1)
          ELSE
            CALL xssben(itype, ig, 1, ng, XsMacsm, lscat1)
          END IF
        ENDDO
      END IF
    ENDIF
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      DO ig = gb, ge
        DO ig2 = GroupInfo%InScatRange(1, ig), GroupInfo%InScatRange(2, ig)
          srcNg(ig - off, ifsr) = srcNg(ig - off, ifsr) + xsmacsm(ig2, ig) * phisNg(ig2, ifsr)
        ENDDO
        IF(lNegSrcFix .AND. srcNg(ig - off, ifsr) .LT. 0._8) THEN
          srcNg(ig - off, ifsr) = srcNg(ig - off, ifsr) - xsmacsm(ig, ig) * phisNg(ig, ifsr)
          xstNg(ig, ifsr) = xstNg(ig, ifsr) - xsmacsm(ig, ig)
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  IF(l3dim) THEN
    IF (nTracerCntl%LkgSplitLv .EQ. 0) THEN
      DO j = 1, nLocalFxr
        ifxr = FxrIdxSt + j -1
        nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
        DO i = 1, nFsrInFxr
          ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
          DO ig = gb, ge
            IF(AxSrc(ipin, iz, ig) .LT. 0 .AND. .NOT. Fxr(ifxr)%lvoid) THEN
              srcNg(ig - off, ifsr) = srcNg(ig - off, ifsr) - AxSrc(ipin, iz, ig)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ELSE
      DO j = 1, nLocalFxr
        ifxr = FxrIdxSt + j -1
        nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
        DO i = 1, nFsrInFxr
          ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
          DO ig = gb, ge
            srcNg(ig - off, ifsr) = srcNg(ig - off, ifsr) - AxSrc(ipin, iz, ig)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDIF
  DO j = 1, CellInfo(icel)%nFsr
    ifsr = FsrIdxSt + j - 1
    srcNg(gb - off : ge - off, ifsr) = srcNg(gb - off : ge - off, ifsr) / xstNg(gb : ge, ifsr)
  ENDDO
  NULLIFY(Xsmacsm)
ENDDO
!$OMP END DO
!$OMP END PARALLEL

IF(.NOT. lxsLib) THEN
  DO i = 1, PE%nThread
    DEALLOCATE(Xsmac(i)%XsMacsm)
  ENDDO
ENDIF

NULLIFY(Pin)
NULLIFY(CellInfo)

END SUBROUTINE SetRtsrcNM_Cusping
! ------------------------------------------------------------------------------------------------------------