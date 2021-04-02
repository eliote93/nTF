#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtSrcNM(Core, Fxr, srcnm, phisnm, psi, AxSrc, xstnm, eigv, iz,gb, ge, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix, PE, Offset)

USE PARAM
USE OMP_LIB
USE TYPEDEF,      ONLY : coreinfo_type, Fxrinfo_type, Cell_Type, pin_Type, GroupInfo_Type, PE_TYPE, XsMac_Type
USE BenchXs,      ONLY : GetChiBen, GetChiDynBen, xssben, xssDynBen
USE MacXsLib_mod, ONLY : MacXsScatMatrix
USE CNTL,         ONLY : nTracerCntl
USE TRAN_MOD,     ONLY : TranInfo

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(GroupInfo_Type):: GroupInfo
TYPE(PE_TYPE)        :: PE

TYPE(FxrInfo_Type), DIMENSION(:) :: Fxr

REAL, POINTER, DIMENSION(:,:)   :: srcnm, phisnm, psi, xstnm
REAL, POINTER, DIMENSION(:,:,:) :: AxSrc

REAL :: eigv

INTEGER :: gb, ge, ng, iz, ifsr, ifxr
LOGICAL :: lxslib, lscat1, l3dim, lNegFix
INTEGER, OPTIONAL :: Offset
! ----------------------------------------------------

TYPE(Pin_Type),  POINTER, DIMENSION(:) :: Pin
TYPE(Cell_Type), POINTER, DIMENSION(:) :: CellInfo

TYPE(XsMac_Type), SAVE, DIMENSION(nThreadMax) :: XsMac

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

srcnm = zero

IF (lxsLib) THEN
  nchi = GroupInfo%nchi
ELSE
  DO i = 1, PE%nThread
    ALLOCATE (XsMac(i)%XsMacSm (ng, ng))
  END DO
END IF
! ----------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED)                                                                                      &
!$OMP PRIVATE(i, j, k, ifsr, ifxr, ipin, icel, ifsrlocal, itype, ig, ig2, tid,                                      &
!$OMP         FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, xsmacsm, chi)
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
        srcnm(ig - off, ifsr) = reigv * chi(ig) * psi(ifsr, iz)
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
          srcnm(ig - off, ifsr) = srcnm(ig - off, ifsr) + xsmacsm(ig2, ig) * phisnm(ig2, ifsr)
        END DO
        
        IF (lNegSrcFix .AND. srcnm(ig - off, ifsr).LT.ZERO) THEN
          srcnm(ig- off, ifsr) = srcnm(ig - off, ifsr) - xsmacsm(ig, ig) * phisnm(ig, ifsr)
          xstnm(ig,      ifsr) = xstnm(ig,       ifsr) - xsmacsm(ig, ig)
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
            IF (AxSrc(ipin, iz, ig).LT.ZERO .AND. .NOT.Fxr(ifxr)%lvoid) srcnm(ig - off, ifsr) = srcnm(ig - off, ifsr) - AxSrc(ipin, iz, ig)
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
            srcnm(ig - off, ifsr) = srcnm(ig - off, ifsr) - AxSrc(ipin, iz, ig)
          END DO
        END DO
      END DO
    END IF
  END IF
  ! --------------------------------------------------
  DO j = 1, CellInfo(icel)%nFsr
    ifsr = FsrIdxSt + j - 1
    
    srcnm(gb - off:ge - off, ifsr) = srcnm(gb - off:ge - off, ifsr) / xstnm(gb:ge, ifsr)
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

END SUBROUTINE SetRtSrcNM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtP1SrcNM(Core, Fxr, srcmnm, phimnm, xstnm, iz, gb, ge, ng, GroupInfo, lxslib, ScatOd, PE, Offset)

USE PARAM
USE OMP_LIB
USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, Cell_Type, pin_Type, GroupInfo_Type, XsMac_Type, PE_Type
USE BenchXs, ONLY : GetChiBen, xssben, xssm1ben, xssm2ben, xssm3ben, xssm1DynBen, xssm2DynBen, xssm3DynBen

USE MacXsLib_mod,   ONLY : MacP1XsScatMatrix, MacP2XsScatMatrix, MacP3XsScatMatrix
USE XsUtil_mod,     ONLY : GetXsMacDat, ReturnXsMacDat
USE BasicOperation, ONLY : CP_CA
USE TRAN_MOD,       ONLY : TranInfo, TranCntl

IMPLICIT NONE

TYPE(CoreInfo_Type)  :: Core
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_Type)        :: PE

TYPE(FxrInfo_Type), DIMENSION(:) :: Fxr

REAL :: eigv

REAL, POINTER, DIMENSION(:,:)   :: xstnm
REAL, POINTER, DIMENSION(:,:,:) :: srcmnm, phimnm

INTEGER :: gb, ge, ng, iz, ifsr, ifxr, ScatOd
LOGICAL :: lxslib, lscat1, l3dim
INTEGER, OPTIONAL :: Offset ! CNJ Edit : GPU Acceleration
! ----------------------------------------------------
TYPE(XsMac_Type), POINTER :: XsMac

TYPE(Pin_Type),  POINTER, DIMENSION(:) :: Pin
TYPE(Cell_Type), POINTER, DIMENSION(:) :: CellInfo

INTEGER :: nCoreFsr, nCoreFxr, nxy, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, nchi, xyb, xye, off
INTEGER :: ipin, icel, ifsrlocal, itype, ig, ig2
INTEGER :: i, j, k

REAL :: reigv, xstinv

REAL, DIMENSION(ng) :: chi
REAL, DIMENSION(9)  :: scatsrc
REAL, POINTER, DIMENSION(:,:) :: xsmacP1sm, xsmacP2sm, xsmacP3sm
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

reigv  = ONE / eigv
srcmnm = zero

IF(lxsLib) nchi = GroupInfo%nchi
! ----------------------------------------------------
!!$OMP PARALLEL DEFAULT(SHARED)      &
!!$OMP PRIVATE(XsMac, ifsr, ifxr, icel, ifsrlocal, itype, scatSrc, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, XsMacP1Sm, XsMacP2Sm, XsMacP3Sm)
IF (.NOT. lxsLib) THEN
  ALLOCATE (XsMacP1sm (ng, ng))
  IF(ScatOd .GE. 2) ALLOCATE (XsMacP2sm (ng, ng))
  IF(ScatOd .EQ. 3) ALLOCATE (XsMacP3sm (ng, ng))
  
  XsMacP1sm = ZERO
  
  IF (ScatOd .GE. 2) XsMacP2sm = ZERO
  IF (ScatOd .EQ. 3) XsMacP3sm = ZERO
END IF
!!$OMP DO
DO ipin = xyb, xye
  CALL GetXsMacDat(XsMac, ng, TRUE)
  
  FsrIdxSt  = Pin(ipin)%FsrIdxSt
  FxrIdxSt  = Pin(ipin)%FxrIdxSt
  icel      = Pin(ipin)%Cell(iz)
  nlocalFxr = CellInfo(icel)%nFxr
  
  DO j = 1, nLocalFxr
    ifxr       = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    
    IF(lXsLib) Then
      CALL MacP1XsScatMatrix(XsMac, Fxr(ifxr), gb, ge, ng, GroupInfo)
      
      IF (ScatOd .GE. 2) CALL MacP2XsScatMatrix(XsMac, Fxr(ifxr), gb, ge, ng, GroupInfo)
      IF (ScatOd .EQ. 3) CALL MacP3XsScatMatrix(XsMac, Fxr(ifxr), gb, ge, ng, GroupInfo)
      
      XsMacP1Sm => XsMac%XsMacP1Sm
      XsMacP2Sm => XsMac%XsMacP2Sm
      XsMacP3Sm => XsMac%XsMacP3Sm
    ELSE
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1, j)
      itype     = Fxr(ifxr)%imix
      
      IF(TranCntl%lDynamicBen) THEN
        CALL XsSm1DynBen(itype, TranInfo%fuelTemp(ipin, iz), 1, ng, 1, ng, XsMacP1Sm)
      ELSE
        CALL XsSm1Ben(itype, 1, ng, 1, ng, XsMacP1Sm)
      END IF
      
      IF (ScatOd .GE. 2) THEN
        IF (TranCntl%lDynamicBen) THEN
          CALL XsSm2DynBen(itype, TranInfo%fuelTemp(ipin, iz), 1, ng, 1, ng, XsMacP2Sm)   !!  OPTIMZIE
        ELSE
          CALL XsSm2Ben(itype, 1, ng, 1, ng, XsMacP2Sm)   !!  OPTIMZIE
        END IF
      END IF
      
      IF (ScatOd .EQ. 3) THEN
        IF (TranCntl%lDynamicBen) THEN
          CALL XsSm3DynBen(itype, TranInfo%fuelTemp(ipin, iz), 1, ng, 1, ng, XsMacP3Sm)
        ELSE
          CALL XsSm3Ben(itype, 1, ng, 1, ng, XsMacP3Sm)
        END IF
      END IF
    END IF
    
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      
      DO ig = gb, ge
        IF (ScatOd .EQ. 1) THEN
          scatSrc = ZERO
          
          DO ig2 = 1, ng
            scatSrc(1:2) = scatSrc(1:2) + XsMacP1Sm(ig2, ig) * phimnm(1:2, ig2, ifsr)
          END DO
          
          srcmnm(1:2, ig - off, ifsr) = scatSrc(1:2)
        ELSE IF (ScatOd .EQ. 2) THEN
          scatSrc = ZERO
          
          DO ig2 = 1, ng
            scatSrc(1:2) = scatSrc(1:2) + XsMacP1Sm(ig2, ig) * phimnm(1:2, ig2, ifsr)
            scatSrc(3:5) = scatSrc(3:5) + XsMacP2Sm(ig2, ig) * phimnm(3:5, ig2, ifsr)
          ENDDO
          
          srcmnm(1:5, ig - off, ifsr) = scatSrc(1:5)
        ELSE IF(ScatOd .EQ. 3) THEN
          scatSrc = ZERO
          
          DO ig2 = 1, ng
            scatSrc(1:2) = scatSrc(1:2) + XsMacP1Sm(ig2, ig) * phimnm(1:2, ig2, ifsr)
            scatSrc(3:5) = scatSrc(3:5) + XsMacP2Sm(ig2, ig) * phimnm(3:5, ig2, ifsr)
            scatSrc(6:9) = scatSrc(6:9) + XsMacP3Sm(ig2, ig) * phimnm(6:9, ig2, ifsr)
          ENDDO
          
          srcmnm(1:9, ig - off, ifsr) = scatSrc(1:9)
        END IF
      END DO
    END DO
  END DO
  
  CALL ReturnXsMacDat(XsMac)
ENDDO
!!$OMP END DO
IF (.NOT. lxsLib) DEALLOCATE (XsMacP1sm) ! modified because of crash! in benchmark XS
IF (.NOT. lxsLib .AND. ScatOd .GE. 2) DEALLOCATE (XsMacP2sm)
IF (.NOT. lxsLib .AND. ScatOd .EQ. 3) DEALLOCATE (XsMacP3sm)
!!$OMP END PARALLEL
! ----------------------------------------------------
IF (ScatOd .EQ. 1) THEN
  DO ipin = xyb, xye
    FsrIdxSt = Pin(ipin)%FsrIdxSt
    icel     = Pin(ipin)%Cell(iz);
    
    DO j = 1, CellInfo(icel)%nFsr
      ifsr = FsrIdxSt + j - 1
      
      DO ig = gb, ge
        xstinv = ONE / xstnm(ig, ifsr)
        
        srcmnm(1:2, ig - off, ifsr) = 3._8 * srcmnm(1:2, ig - off, ifsr) * xstinv
      END DO
    END DO
  END DO
ELSE IF (ScatOd .EQ. 2) THEN
  DO ipin = xyb, xye
    FsrIdxSt = Pin(ipin)%FsrIdxSt
    icel     = Pin(ipin)%Cell(iz);
    
    DO j = 1, CellInfo(icel)%nFsr
      ifsr = FsrIdxSt + j - 1
      
      DO ig = gb, ge
        xstinv = ONE / xstnm(ig, ifsr)
        
        srcmnm(1:2, ig - off, ifsr) = 3._8 * srcmnm(1:2, ig - off, ifsr) * xstinv
        srcmnm(3:5, ig - off, ifsr) = 5._8 * srcmnm(3:5, ig - off, ifsr) * xstinv
      END DO
    END DO
  END DO
ELSE IF (ScatOd .EQ. 3) THEN
  DO ipin = xyb, xye
    FsrIdxSt = Pin(ipin)%FsrIdxSt
    icel     = Pin(ipin)%Cell(iz);
    
    DO j = 1, CellInfo(icel)%nFsr
      ifsr = FsrIdxSt + j - 1
      
      DO ig = gb, ge
        xstinv = ONE / xstnm(ig, ifsr)
        
        srcmnm(1:2, ig - off, ifsr) = 3._8 * srcmnm(1:2, ig - off, ifsr) * xstinv
        srcmnm(3:5, ig - off, ifsr) = 5._8 * srcmnm(3:5, ig - off, ifsr) * xstinv
        srcmnm(6:9, ig - off, ifsr) = 7._8 * srcmnm(6:9, ig - off, ifsr) * xstinv
      END DO
    END DO
  END DO
END IF
! ----------------------------------------------------
IF (lxsLib) NULLIFY (XsMacP1Sm)
NULLIFY (Pin)
NULLIFY (CellInfo)
! ----------------------------------------------------

END SUBROUTINE SetRtP1SrcNM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtSrcNM_Cusping(Core, FmInfo, Fxr, srcnm, phisnm, psi, AxSrc, xstnm, eigv, iz, gb, ge, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix, PE, Offset)

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
REAL, POINTER :: srcnm(:, :), phisnm(:, :), psi(:, :), AxSrc(:, :, :), xstnm(:, :)
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

srcnm = zero
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
        srcnm(ig - off, ifsr) = reigv * chi(ig) * psi(ifsr, iz)
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
          srcnm(ig - off, ifsr) = srcnm(ig - off, ifsr) + xsmacsm(ig2, ig) * phisnm(ig2, ifsr)
        ENDDO
        IF(lNegSrcFix .AND. srcnm(ig - off, ifsr) .LT. 0._8) THEN
          srcnm(ig - off, ifsr) = srcnm(ig - off, ifsr) - xsmacsm(ig, ig) * phisnm(ig, ifsr)
          xstnm(ig, ifsr) = xstnm(ig, ifsr) - xsmacsm(ig, ig)
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
              srcnm(ig - off, ifsr) = srcnm(ig - off, ifsr) - AxSrc(ipin, iz, ig)
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
            srcnm(ig - off, ifsr) = srcnm(ig - off, ifsr) - AxSrc(ipin, iz, ig)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDIF
  DO j = 1, CellInfo(icel)%nFsr
    ifsr = FsrIdxSt + j - 1
    srcnm(gb - off : ge - off, ifsr) = srcnm(gb - off : ge - off, ifsr) / xstnm(gb : ge, ifsr)
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

END SUBROUTINE SetRtSrcNM_Cusping
! ------------------------------------------------------------------------------------------------------------