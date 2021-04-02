#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtSrcGM(Core, Fxr, src, phis, psi, axsrc, xstr1g, eigv, iz, ig, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix, PE)

USE PARAM
USE OMP_LIB
USE TYPEDEF,        ONLY : coreinfo_type, Fxrinfo_type, Cell_Type, pin_Type, GroupInfo_Type, PE_TYPE, XsMac_Type
USE BenchXs,        ONLY : GetChiBen, getChiDynBen, xssben, xssDynBen
USE MacXsLib_mod,   ONLY : MacXsScatMatrix
USE BasicOperation, ONLY : CP_CA
USE CNTL,           ONLY : nTracerCntl
USE TRAN_MOD,       ONLY : TranInfo

USE MOC_MOD, ONLY : phissave

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(GroupInfo_Type):: GroupInfo
TYPE(PE_TYPE)       :: PE

TYPE(FxrInfo_Type), DIMENSION(:) :: Fxr

REAL, POINTER, DIMENSION(:)     :: src, AxSrc, xstr1g
REAL, POINTER, DIMENSION(:,:)   :: psi
REAL, POINTER, DIMENSION(:,:,:) :: phis

REAL :: eigv

INTEGER :: myzb, myze, ig, ng, iz, ifsr, ifxr, fsridx
LOGICAL :: lxslib, lscat1, l3dim, lNegFix
! ----------------------------------------------------
TYPE(XsMac_Type), SAVE, DIMENSION(nTHREADMAX) :: XsMac

TYPE(Pin_Type),  POINTER, DIMENSION(:) :: Pin
TYPE(Cell_Type), POINTER, DIMENSION(:) :: CellInfo

INTEGER :: nxy, nCoreFsr, nCoreFxr, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, nchi
INTEGER :: ipin, icel, ifsrlocal, itype, ig2, tid
INTEGER :: i, j, k

LOGICAL :: lNegSrcFix

REAL :: reigv, psrc, pvol

REAL, DIMENSION(ng) ::  chi
REAL, POINTER, DIMENSION(:,:) :: xsmacsm
! ----------------------------------------------------

Pin      => Core%Pin
CellInfo => Core%CellInfo
nCoreFsr  = Core%nCoreFsr
nCoreFxr  = Core%nCoreFxr
nxy       = Core%nxy

reigv = ONE / eigv

lNegSrcFix = FALSE
IF(lNegFix) lNegSrcFix = TRUE

IF (.NOT. lxsLib) THEN
  DO i = 1, PE%nThread
    ALLOCATE (XsMac(i)%XsMacSm (ng, ng))
  END DO
ELSE
  nchi = GroupInfo%nchi
END IF
! ----------------------------------------------------
tid = 1

!$  call omp_set_dynamic(FALSE)
!$  call omp_set_num_threads(PE%nThread)

!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(i, j, k, ifsr, ifxr, ipin, icel, ifsrlocal, itype, ig2, tid, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, xsmacsm, chi, psrc, pvol, fsridx)
!$  tid = omp_get_thread_num() + 1
!$OMP DO
DO ipin = 1, nxy
  FsrIdxSt  = Pin(ipin)%FsrIdxSt
  FxrIdxSt  = Pin(ipin)%FxrIdxSt
  icel      = Pin(ipin)%Cell(iz)
  nlocalFxr = CellInfo(icel)%nFxr
  ! --------------------------------------------------
  ! Fission Source
  DO j = 1, nLocalFxr
    ifxr      = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    
    IF (lXsLib) Then
      CHI(ig:ig) = 0
      
      IF (ig .LE. nchi .AND. Fxr(ifxr)%ldepl) CHI(ig:ig) = Fxr(ifxr)%chi(ig)
    ELSE
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1,j)
      itype     = Fxr(ifxr)%imix
      
      IF(nTracerCntl%lDynamicBen) THEN
        CHI(ig:ig) = GetChiDynBen(itype, TranInfo%fuelTemp(ipin, iz), ig, ig)
      ELSE
        CHI(ig:ig) = GetChiBen(itype, ig, ig)
      END IF
    END IF

    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      
      src(ifsr) = reigv * chi(ig) * psi(ifsr, iz)
    END DO
  END DO
  ! --------------------------------------------------
  !Scattering Source Update
  DO j = 1, nLocalFxr
    ifxr      = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    
    IF(lXsLib) Then
      CALL MacXsScatMatrix(XsMac(tid), Fxr(ifxr), ig, ig, ng, GroupInfo, lscat1, TRUE)
      
      XsMacSm => XsMac(tid)%XsMacSm
#ifdef inflow
      XsMacSm(ig, ig) = XsMacSm(ig, ig) + Fxr(ifxr)%DelInflow(ig)
#endif
    ELSE
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1,j)
      itype     = Fxr(ifxr)%imix
      XsMacSm  => XsMac(tid)%XsMacSm
      
      IF (nTracerCntl%lDynamicBen) THEN
        CALL xssDynben(itype, TranInfo%fuelTemp(ipin, iz), ig, 1, ng, XsMacsm, lscat1)
      ELSE
        CALL xssben(itype, ig, 1, ng, XsMacsm, lscat1)
      END IF
    END IF
    
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      
      DO ig2 = GroupInfo%InScatRange(1, ig), GroupInfo%InScatRange(2, ig)
        src(ifsr) = src(ifsr) + xsmacsm(ig2, ig) * phis(ifsr, iz, ig2)
        !src(ifsr) = src(ifsr) + xsmacsm(ig2, ig) * phissave(ifsr, iz, ig2) DEBUG
      END DO
      
      IF (lNegSrcFix .AND. src(ifsr) .LT. ZERO) THEN
        src   (ifsr) = src   (ifsr) - xsmacsm(ig, ig) * phis(ifsr, iz, ig)
        !src   (ifsr) = src   (ifsr) - xsmacsm(ig, ig) * phissave(ifsr, iz, ig) ! DEBUG
        xstr1g(ifsr) = xstr1g(ifsr) - xsmacsm(ig, ig)
      END IF
    END DO
  END DO
  ! --------------------------------------------------
  IF (l3dim) THEN
#ifdef LkgSplit
  SELECT CASE (nTracerCntl%LkgSplitLv)
    CASE(0) ! default  : if AxSrc < 0
      DO j = 1, nLocalFxr
        ifxr      = FxrIdxSt + j -1
        nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
        
        DO i = 1, nFsrInFxr
          ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
          
          IF (AxSrc(ipin).LT.0 .AND. .NOT. Fxr(ifxr)%lvoid) src(ifsr) = src(ifsr) - AxSrc(ipin)
        END DO
      END DO
    CASE(1)
      DO j = 1, nLocalFxr
        ifxr      = FxrIdxSt + j -1
        nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
        
        DO i = 1, nFsrInFxr
          ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
          
          src(ifsr) = src(ifsr) - AxSrc(ipin)
        END DO
      END DO
    CASE (2) ! radial split by src
      !-- 1) calc. pinwise total source
      psrc = ZERO
      pvol = ZERO
      
      DO j = 1, nLocalFxr
        ifxr      = FxrIdxSt + j -1
        nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
        
        DO i = 1, nFsrInFxr
          ifsr    = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
          fsridx = Cellinfo(icel)%MapFxr2FsrIdx(i, j)
          
          !-- Source proportional
          psrc = psrc + Cellinfo(icel)%Vol(fsridx) * src(ifsr) ! pin-wise total source except axial source term
          pvol = pvol + Cellinfo(icel)%Vol(fsridx)
        END DO
      END DO
      
      psrc = psrc / pvol ! volume averaged source
      
      DO j = 1, nLocalFxr
        ifxr      = FxrIdxSt + j -1
        nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
        
        DO i = 1, nFsrInFxr
          ifsr   = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
          fsridx = Cellinfo(icel)%MapFxr2FsrIdx(i, j)
          
          src(ifsr) = src(ifsr) * (ONE - AxSrc(ipin) / psrc)
        END DO
      END DO
    CASE (3) ! radial split by Dphi
      !-- 1) calc. pinwise total source
      psrc = ZERO
      pvol = ZERO
      
      DO j = 1, nLocalFxr
        ifxr      = FxrIdxSt + j -1
        nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
        
        DO i = 1, nFsrInFxr
          ifsr   = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
          fsridx = Cellinfo(icel)%MapFxr2FsrIdx(i, j)
          
          !-- D*phi proportional
          psrc = psrc + Cellinfo(icel)%Vol(fsridx) * (ONE / 3._8 / xstr1g(ifsr) * phis(ifsr, iz, ig)) ! pin-wise total source except axial source term
          !psrc = psrc + Cellinfo(icel)%Vol(fsridx) * (ONE / 3._8 / xstr1g(ifsr) * phissave(ifsr, iz, ig)) ! DEBUG
          pvol = pvol + Cellinfo(icel)%Vol(fsridx)
        END DO
      END DO
      
      psrc = psrc / pvol ! volume averaged source
      
      DO j = 1, nLocalFxr
        ifxr      = FxrIdxSt + j -1
        nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
        
        DO i = 1, nFsrInFxr
          ifsr   = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
          fsridx = Cellinfo(icel)%MapFxr2FsrIdx(i, j)
          
          src(ifsr) = src(ifsr) - AxSrc(ipin) * (ONE / 3._8 / xstr1g(ifsr) * phis(ifsr, iz, ig)) / psrc
          !src(ifsr) = src(ifsr) - AxSrc(ipin) * (ONE / 3._8 / xstr1g(ifsr) * phissave(ifsr, iz, ig)) / psrc ! DEBUG
        END DO
      END DO
    CASE (4) ! radial split by D
      !-- 1) calc. pinwise total source
      psrc = ZERO
      pvol = ZERO
      
      DO j = 1, nLocalFxr
        ifxr      = FxrIdxSt + j -1
        nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
        
        DO i = 1, nFsrInFxr
          ifsr   = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
          fsridx = Cellinfo(icel)%MapFxr2FsrIdx(i, j)
          
          !-- D*phi proportional
          psrc = psrc + Cellinfo(icel)%Vol(fsridx) * (ONE / 3._8 / xstr1g(ifsr)) ! pin-wise total source except axial source term
          pvol = pvol + Cellinfo(icel)%Vol(fsridx)
        END DO
      END DO
      
      psrc = psrc / pvol ! volume averaged source
      
      DO j = 1, nLocalFxr
        ifxr      = FxrIdxSt + j -1
        nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
        
        DO i = 1, nFsrInFxr
          ifsr   = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
          fsridx = Cellinfo(icel)%MapFxr2FsrIdx(i, j)
          
          src(ifsr) = src(ifsr) - AxSrc(ipin) * (ONE / 3._8 / xstr1g(ifsr)) / psrc
        END DO
      END DO
  END SELECT
#endif
  END IF
! --------------------------------------------------
  DO j = 1, CellInfo(icel)%nFsr
    ifsr = FsrIdxSt + j - 1
    
    src(ifsr) = src(ifsr) / xstr1g(ifsr)
  END DO
  
  NULLIFY (Xsmacsm)
END DO
!$OMP END DO
!$OMP END PARALLEL
! ----------------------------------------------------
IF (.NOT. lxsLib) THEN
  DO i = 1, PE%nThread
    DEALLOCATE (Xsmac(i)%XsMacsm)
  END DO
END IF

NULLIFY (Pin)
NULLIFY (CellInfo)
! ----------------------------------------------------

END SUBROUTINE SetRtSrcGM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtP1SrcGM(Core, Fxr, srcm, phim, xstr1g, iz, ig, ng, GroupInfo, l3dim, lxslib, lscat1, ScatOd, PE)

USE PARAM
USE OMP_LIB
USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, Cell_Type, pin_Type, GroupInfo_Type, XsMac_Type, PE_Type
USE BenchXs, ONLY : GetChiBen, xssben, xssm1ben, xssm2ben, xssm3ben, xssm1DynBen, xssm2DynBen, xssm3DynBen
USE MacXsLib_mod,   ONLY : MacP1XsScatMatrix, MacP2XsScatMatrix, MacP3XsScatMatrix
USE XsUtil_mod,     ONLY : GetXsMacDat, ReturnXsMacDat, FreeXsMac
USE BasicOperation, ONLY : CP_CA
USE TRAN_MOD,       ONLY : TranInfo, TranCntl

USE MOC_MOD, ONLY : phimsave

IMPLICIT NONE

TYPE(CoreInfo_Type)  :: Core
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_Type)        :: PE

TYPE(FxrInfo_Type), DIMENSION(:) :: Fxr

REAL, POINTER, DIMENSION(:)       :: xstr1g
REAL, POINTER, DIMENSION(:,:)     :: srcm
REAL, POINTER, DIMENSION(:,:,:,:) :: phim

REAL :: eigv

INTEGER :: myzb, myze, ig, ng, iz, ifsr, ifxr, ScatOd
LOGICAL :: lxslib, lscat1, l3dim
! ----------------------------------------------------
TYPE(XsMac_Type), POINTER :: XsMac

TYPE(Pin_Type),  POINTER, DIMENSION(:) :: Pin
TYPE(Cell_Type), POINTER, DIMENSION(:) :: CellInfo

INTEGER :: nxy, nCoreFsr, nCoreFxr, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, nchi
INTEGER :: ipin, icel, ifsrlocal, itype, ig2, tid
INTEGER :: i, j, k

REAL :: reigv

REAL, DIMENSION(ng) ::  chi

REAL, POINTER, DIMENSION(:,:) :: XsMacP1sm, XsMacP2sm, XsMacP3sm
! ----------------------------------------------------

Pin      => Core%Pin
CellInfo => Core%CellInfo
nCoreFsr  = Core%nCoreFsr
nCoreFxr  = Core%nCoreFxr
nxy       = Core%nxy

reigv = ONE / eigv
srcm  = ZERO

IF (lxsLib) nchi = GroupInfo%nchi
! ----------------------------------------------------
tid = 1

!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(PE%nThread)
!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(XsMac, i, j, k, ifsr, ifxr, ipin, icel, ifsrlocal, itype, ig2, tid, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, XsMacP1Sm, XsMacP2Sm, XsMacP3Sm)
!$  tid = omp_get_thread_num() + 1

IF (.NOT. lxsLib) THEN
  ALLOCATE (XsMacP1sm (ng, ng))
  
  IF (ScatOd .GE. 2) ALLOCATE (XsMacP2sm (ng, ng))
  IF (ScatOd .EQ. 3) ALLOCATE (XsMacP3sm (ng, ng))
  
  XsMacP1sm = ZERO
  
  IF (ScatOd .GE. 2) XsMacP2sm = ZERO
  IF (ScatOd .EQ. 3) XsMacP3sm = ZERO
END IF

!$OMP DO
DO ipin = 1, nxy
  CALL GetXsMacDat(XsMac, ng, TRUE)
  
  FsrIdxSt  = Pin(ipin)%FsrIdxSt
  FxrIdxSt  = Pin(ipin)%FxrIdxSt
  icel      = Pin(ipin)%Cell(iz)
  nlocalFxr = CellInfo(icel)%nFxr
  
  DO j = 1, nLocalFxr
    ifxr      = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    
    IF (lXsLib) Then
      CALL MacP1XsScatMatrix(XsMac, Fxr(ifxr), ig, ig, ng, GroupInfo)
      
      IF(ScatOd .GE. 2) CALL MacP2XsScatMatrix(XsMac, Fxr(ifxr), ig, ig, ng, GroupInfo)
      IF(ScatOd .EQ. 3) CALL MacP3XsScatMatrix(XsMac, Fxr(ifxr), ig, ig, ng, GroupInfo)
      
      XsMacP1Sm => XsMac%XsMacP1Sm
      XsMacP2Sm => XsMac%XsMacP2Sm
      XsMacP3Sm => XsMac%XsMacP3Sm
    ELSE
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1,j)
      itype     = Fxr(ifxr)%imix
      
      IF (TranCntl%lDynamicBen) THEN
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
      
      IF (ScatOd .EQ. 1) THEN
        DO ig2 = 1, ng
          srcm(1, ifsr) = srcm(1, ifsr) + XsMacP1Sm(ig2, ig) * phim(1, ifsr, iz, ig2)
          srcm(2, ifsr) = srcm(2, ifsr) + XsMacP1Sm(ig2, ig) * phim(2, ifsr, iz, ig2)
          !srcm(1, ifsr) = srcm(1, ifsr) + XsMacP1Sm(ig2, ig) * phimsave(1, ifsr, iz, ig2) ! DEBUG
          !srcm(2, ifsr) = srcm(2, ifsr) + XsMacP1Sm(ig2, ig) * phimsave(2, ifsr, iz, ig) ! DEBUG
        END DO
      ELSE IF(ScatOd .EQ. 2) THEN
        DO ig2 = 1, ng
          srcm(1:2, ifsr) = srcm(1:2, ifsr) + XsMacP1Sm(ig2, ig) * phim(1:2, ifsr, iz, ig2)
          srcm(3:5, ifsr) = srcm(3:5, ifsr) + XsMacP2Sm(ig2, ig) * phim(3:5, ifsr, iz, ig2)
          !srcm(1:2, ifsr) = srcm(1:2, ifsr) + XsMacP1Sm(ig2, ig) * phimsave(1:2, ifsr, iz, ig2) ! DEBUG
          !srcm(3:5, ifsr) = srcm(3:5, ifsr) + XsMacP2Sm(ig2, ig) * phimsave(3:5, ifsr, iz, ig2) ! DEBUG
        END DO
      ELSE IF(ScatOd .EQ. 3) THEN
        DO ig2 = 1, ng
          srcm(1:2, ifsr) = srcm(1:2, ifsr) + XsMacP1Sm(ig2, ig) * phim(1:2, ifsr, iz, ig2)
          srcm(3:5, ifsr) = srcm(3:5, ifsr) + XsMacP2Sm(ig2, ig) * phim(3:5, ifsr, iz, ig2)
          srcm(6:9, ifsr) = srcm(6:9, ifsr) + XsMacP3Sm(ig2, ig) * phim(6:9, ifsr, iz, ig2)
          !srcm(1:2, ifsr) = srcm(1:2, ifsr) + XsMacP1Sm(ig2, ig) * phimsave(1:2, ifsr, iz, ig2) ! DEBUG
          !srcm(3:5, ifsr) = srcm(3:5, ifsr) + XsMacP2Sm(ig2, ig) * phimsave(3:5, ifsr, iz, ig2) ! DEBUG
          !srcm(6:9, ifsr) = srcm(6:9, ifsr) + XsMacP3Sm(ig2, ig) * phimsave(6:9, ifsr, iz, ig2) ! DEBUG
        END DO
      END IF
    END DO
  END DO
  
  CALL ReturnXsMacDat(XsMac) !Memory leak problem !! 17/01/09   big memory but stable
END DO
!$OMP END DO

IF (.NOT. lxsLib) DEALLOCATE(XsMacP1sm) ! modified because of crash! in benchmark XS
IF (.NOT. lxsLib .AND. ScatOd .GE. 2) DEALLOCATE(XsMacP2sm)
IF (.NOT. lxsLib .AND. ScatOd .EQ. 3) DEALLOCATE(XsMacP3sm)
!$OMP END PARALLEL
! ----------------------------------------------------
! Axail Source Contribution
IF (ScatOd .EQ. 1) THEN
  DO ipin = 1, nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt
    icel     = Pin(ipin)%Cell(iz)
    
    DO j = 1, CellInfo(icel)%nFsr
      ifsr = FsrIdxSt + j - 1
      
      srcm(1:2, ifsr) = 3._8 * srcm(1:2, ifsr) / xstr1g(ifsr)
    END DO
  END DO
ELSE IF (ScatOd .EQ. 2) THEN
  DO ipin = 1, nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt
    icel     = Pin(ipin)%Cell(iz)
    
    DO j = 1, CellInfo(icel)%nFsr
      ifsr = FsrIdxSt + j - 1
      
      srcm(1:2, ifsr) = 3._8 * srcm(1:2, ifsr) / xstr1g(ifsr)
      srcm(3:5, ifsr) = 5._8 * srcm(3:5, ifsr) / xstr1g(ifsr)
    END DO
  END DO
ELSE IF (ScatOd .EQ. 3) THEN
  DO ipin = 1, nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt
    icel     = Pin(ipin)%Cell(iz);
    
    DO j = 1, CellInfo(icel)%nFsr
      ifsr = FsrIdxSt + j - 1
      
      srcm(1:2, ifsr) = 3._8 * srcm(1:2, ifsr) / xstr1g(ifsr)
      srcm(3:5, ifsr) = 5._8 * srcm(3:5, ifsr) / xstr1g(ifsr)
      srcm(6:9, ifsr) = 7._8 * srcm(6:9, ifsr) / xstr1g(ifsr)
    END DO
  END DO
END IF
! ----------------------------------------------------
IF (lXslib) NULLIFY (XsMacP1Sm)

NULLIFY (Pin)
NULLIFY (CellInfo)
! ----------------------------------------------------

END SUBROUTINE SetRtP1SrcGM
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtSrcGM_Cusping(Core, FmInfo, Fxr, src, phis, psi, axsrc, xstr1g, eigv, iz, ig, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix, PE)

USE PARAM
USE TYPEDEF,        ONLY : coreinfo_type,          Fxrinfo_type,          Cell_Type,      &
                           pin_Type,               GroupInfo_Type,        PE_TYPE,        &
                           XsMac_Type,             FmInfo_Type
USE BenchXs,        ONLY : GetChiBen,              xssben,                GetChiDynBen,    xssben_Cusping, &
                           MacXsBen,               xssDynBen,             xssDynben_Cusping, DynMacXsBen
USE MacXsLib_mod,   ONLY : MacXsScatMatrix
USE BasicOperation, ONLY : CP_CA
USE CNTL,             ONLY : nTracerCntl
USE OMP_LIB
USE TRAN_MOD,       ONLY : TranInfo
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(FxrInfo_Type) :: Fxr(:)
TYPE(GroupInfo_Type):: GroupInfo
TYPE(PE_TYPE) :: PE
REAL, POINTER :: src(:), phis(:, :, :), psi(:, :), AxSrc(:), xstr1g(:)
REAL :: eigv
INTEGER :: myzb, myze, ig, ng, iz, ifsr, ifxr, fsridx
LOGICAL :: lxslib, lscat1, l3dim
LOGICAL :: lNegFix


TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(XsMac_Type), SAVE :: XsMac(nTHREADMAX)
INTEGER :: nxy, nCoreFsr, nCoreFxr, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, nchi
INTEGER :: ipin, icel, ifsrlocal, itype, ig2, tid
INTEGER :: i, j, k, jg

REAL :: phiz(ng), philz(ng), phiuz(ng)
REAL :: vol, volsum
REAL :: reigv
REAL ::  chi(ng)
REAL, POINTER :: xsmacsm(:,:)

REAL :: psrc, pvol

LOGICAL :: lNegSrcFix
LOGICAL :: lCusping
Pin => Core%Pin
CellInfo => Core%CellInfo
nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr
nxy = Core%nxy

reigv = one/eigv

lNegSrcFix = FALSE
IF(lNegFix) lNegSrcFix = TRUE


IF(.NOT. lxsLib) THEN
  DO i = 1, PE%nThread
    ALLOCATE(XsMac(i)%XsMacSm(ng, ng))
  ENDDO
ENDIF

IF(lxsLib) nchi = GroupInfo%nchi

tid = 1
!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(PE%nThread)
!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(i, j, k, ifsr, ifxr, ipin, icel, ifsrlocal, itype, ig2, tid, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, xsmacsm, chi, psrc, pvol, fsridx, &
!$OMP         vol, volsum, phiz, philz, phiuz)
!$  tid = omp_get_thread_num()+1
!$OMP DO
DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  !Fission Source
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    IF(lXsLib) Then
      IF(ig .gt. nchi) THEN
        CHI(ig:ig) = 0
      ELSE
        CHI(ig:ig) = 0
        IF(Fxr(ifxr)%ldepl) CHI(ig:ig) = Fxr(ifxr)%chi(ig)
      ENDIF
    ELSE
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1,j)
      itype = Fxr(ifxr)%imix
      IF(nTracerCntl%lDynamicBen) THEN
        CHI(ig:ig) = GetChiDynBen(itype, TranInfo%fuelTemp(ipin, iz), ig, ig)
      ELSE
        CHI(ig:ig) = GetChiBen(itype, ig, ig)
      END IF
    ENDIF

    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      src(ifsr) = reigv * chi(ig) * psi(ifsr, iz)
    ENDDO !Fsr Sweep
  ENDDO
  !Scattering Source Update
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    IF(lXsLib) Then
      CALL MacXsScatMatrix(XsMac(tid), Fxr(ifxr), ig, ig, ng, GroupInfo, lscat1, .TRUE.)
      XsMacSm => XsMac(tid)%XsMacSm
#ifdef inflow
      XsMacSm(ig, ig) = XsMacSm(ig, ig) + Fxr(ifxr)%DelInflow(ig)
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
        DO jg = 1, ng
          DO i = 1, nFsrInFxr
            ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(i,j)
            ifsr = FsrIdxSt + ifsrlocal - 1
            vol = CellInfo(icel)%vol(ifsrlocal)
            IF(jg .EQ. 1) volsum = volsum + vol
            phiz(jg) = phiz(jg) + phis(ifsr, iz, jg) * vol
            IF(iz .EQ. myzb) THEN
              philz(jg) = philz(jg) + FmInfo%neighphis(ifsr, jg, BOTTOM) * vol
            ELSE
              philz(jg) = philz(jg) + phis(ifsr, iz-1, jg) * vol
            END IF
            IF(iz .EQ. myze) THEN
              phiuz(jg) = phiuz(jg) + FmInfo%neighphis(ifsr, jg, TOP) * vol
            ELSE
              phiuz(jg) = phiuz(jg) + phis(ifsr, iz+1, jg) * vol
            END IF
          END DO
          phiz(jg) = phiz(jg) / volsum
          philz(jg) = philz(jg) / volsum
          phiuz(jg) = phiuz(jg) / volsum
        END DO
        IF(nTracerCntl%lDynamicBen) THEN
          CALL xssDynben_Cusping(itype, TranInfo%fuelTemp(ipin, iz), ig, 1, ng, XsMacsm, lscat1, phiz, philz, phiuz, &
            Core%hzfm(iz), Core%hzfm(iz-1), Core%hzfm(iz+1))
        ELSE
          CALL xssben_Cusping(itype, ig, 1, ng, XsMacsm, lscat1, phiz, philz, phiuz, &
            Core%hzfm(iz), Core%hzfm(iz-1), Core%hzfm(iz+1))
        END IF
      ELSE
        IF(nTracerCntl%lDynamicBen) THEN
          CALL xssDynben(itype, TranInfo%fuelTemp(ipin, iz), ig, 1, ng, XsMacsm, lscat1)
        ELSE
          CALL xssben(itype, ig, 1, ng, XsMacsm, lscat1)
        END IF
      END IF
    ENDIF

    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      DO ig2 = GroupInfo%InScatRange(1, ig), GroupInfo%InScatRange(2, ig)
        src(ifsr) = src(ifsr) + xsmacsm(ig2, ig)*phis(ifsr, iz, ig2)
      ENDDO
      IF(lNegSrcFix .AND. src(ifsr) .LT. 0._8) THEN
        src(ifsr) = src(ifsr) - xsmacsm(ig, ig)*phis(ifsr, iz, ig)
        xstr1g(ifsr) = xstr1g(ifsr) - xsmacsm(ig, ig)
      ENDIF
    ENDDO !Fsr Sweep

  ENDDO !End of Fxr Sweep

  IF (l3dim) THEN
#ifdef LkgSplit
  SELECT CASE(nTracerCntl%LkgSplitLv)
    CASE (4) ! radial split by D
      !-- 1) calc. pinwise total source
      psrc=0.0_8;
      pvol=0.0_8;
      DO j = 1, nLocalFxr
        ifxr = FxrIdxSt + j -1
        nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
        DO i = 1, nFsrInFxr
          ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
          fsridx=Cellinfo(icel)%MapFxr2FsrIdx(i, j)
          !-- D*phi proportional
          psrc=psrc+(1._8/3._8/xstr1g(ifsr))*Cellinfo(icel)%Vol(fsridx);    ! pin-wise total source except axial source term
          pvol=pvol+Cellinfo(icel)%Vol(fsridx)
        ENDDO
      ENDDO
      psrc=psrc/pvol ! volume averaged source
      DO j = 1, nLocalFxr
        ifxr = FxrIdxSt + j -1
        nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
        DO i = 1, nFsrInFxr
          ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
          fsridx=Cellinfo(icel)%MapFxr2FsrIdx(i, j)
          src(ifsr) = src(ifsr)- AxSrc(ipin)*(1._8/3._8/xstr1g(ifsr))/psrc
        ENDDO
      ENDDO
    CASE (3) ! radial split by Dphi
      !-- 1) calc. pinwise total source
      psrc=0.0_8;
      pvol=0.0_8;
      DO j = 1, nLocalFxr
        ifxr = FxrIdxSt + j -1
        nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
        DO i = 1, nFsrInFxr
          ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
          fsridx=Cellinfo(icel)%MapFxr2FsrIdx(i, j)
          !-- D*phi proportional
          psrc=psrc+(1._8/3._8/xstr1g(ifsr)*phis(ifsr, iz, ig))*Cellinfo(icel)%Vol(fsridx);    ! pin-wise total source except axial source term
          pvol=pvol+Cellinfo(icel)%Vol(fsridx)
        ENDDO
      ENDDO
      psrc=psrc/pvol ! volume averaged source
      DO j = 1, nLocalFxr
        ifxr = FxrIdxSt + j -1
        nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
        DO i = 1, nFsrInFxr
          ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
          fsridx=Cellinfo(icel)%MapFxr2FsrIdx(i, j)
          src(ifsr) = src(ifsr)- AxSrc(ipin)*(1._8/3._8/xstr1g(ifsr)*phis(ifsr, iz, ig))/psrc
        ENDDO
      ENDDO
    CASE (2) ! radial split by src
      !-- 1) calc. pinwise total source
      psrc=0.0_8;
      pvol=0.0_8;
      DO j = 1, nLocalFxr
        ifxr = FxrIdxSt + j -1
        nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
        DO i = 1, nFsrInFxr
          ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
          fsridx=Cellinfo(icel)%MapFxr2FsrIdx(i, j)
          !-- Source proportional
          psrc=psrc+src(ifsr)*Cellinfo(icel)%Vol(fsridx);    ! pin-wise total source except axial source term
          pvol=pvol+Cellinfo(icel)%Vol(fsridx)
        ENDDO
      ENDDO
      psrc=psrc/pvol ! volume averaged source
      DO j = 1, nLocalFxr
        ifxr = FxrIdxSt + j -1
        nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
        DO i = 1, nFsrInFxr
          ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
          fsridx=Cellinfo(icel)%MapFxr2FsrIdx(i, j)
          src(ifsr) = src(ifsr)*(1.0_8- AxSrc(ipin)/psrc)
        ENDDO
      ENDDO
    CASE(0) ! default  : if AxSrc<0
      DO j = 1, nLocalFxr
        ifxr = FxrIdxSt + j -1
        nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
        DO i = 1, nFsrInFxr
          ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
          IF(AxSrc(ipin) .LT. 0 .AND. .NOT. Fxr(ifxr)%lvoid) THEN
            src(ifsr) = src(ifsr) - AxSrc(ipin)
          ENDIF
        ENDDO
      ENDDO
    CASE(1)
      DO j = 1, nLocalFxr
        ifxr = FxrIdxSt + j -1
        nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
        DO i = 1, nFsrInFxr
          ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
          src(ifsr) = src(ifsr) - AxSrc(ipin)
      ENDDO
    ENDDO
  ENDSELECT
#endif
  ENDIF

  DO j = 1, CellInfo(icel)%nFsr
    ifsr = FsrIdxSt + j - 1
    src(ifsr) = src(ifsr)/xstr1g(ifsr)
  ENDDO
  NULLIFY(Xsmacsm)
ENDDO !End of Pin
!$OMP END DO
!$OMP END PARALLEL

IF(.NOT. lxsLib) THEN
  DO i = 1, PE%nThread
    DEALLOCATE(Xsmac(i)%XsMacsm)
  ENDDO
ENDIF

NULLIFY(Pin)
NULLIFY(CellInfo)

END SUBROUTINE SetRtSrcGM_Cusping
! ------------------------------------------------------------------------------------------------------------