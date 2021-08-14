#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtSrcGM(Core, Fxr, src1g, phis, psi, axsrc1g, xst1g, eigv, iz, ig, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix, PE)

USE OMP_LIB
USE PARAM,        ONLY : ZERO, TRUE, FALSE, ONE, nTHREADMAX
USE TYPEDEF,      ONLY : coreinfo_type, Fxrinfo_type, Cell_Type, pin_Type, GroupInfo_Type, PE_TYPE, XsMac_Type
USE BenchXs,      ONLY : GetChiBen, getChiDynBen, xssben, xssDynBen
USE MacXsLib_mod, ONLY : MacXsScatMatrix
USE CNTL,         ONLY : nTracerCntl
USE TRAN_MOD,     ONLY : TranInfo

IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(GroupInfo_Type):: GroupInfo
TYPE(PE_TYPE)       :: PE

TYPE(FxrInfo_Type), DIMENSION(:) :: Fxr

REAL, POINTER, DIMENSION(:)     :: src1g, AxSrc1g, xst1g
REAL, POINTER, DIMENSION(:,:)   :: psi
REAL, POINTER, DIMENSION(:,:,:) :: phis

REAL :: eigv

INTEGER :: myzb, myze, ig, ng, iz, ifsr, ifxr, fsridx
LOGICAL :: lxslib, lscat1, l3dim, lNegFix
! ----------------------------------------------------
TYPE (XsMac_Type), SAVE, DIMENSION(nTHREADMAX) :: XsMac

TYPE (Pin_Type),  POINTER, DIMENSION(:) :: Pin
TYPE (Cell_Type), POINTER, DIMENSION(:) :: CellInfo

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
IF (lNegFix) lNegSrcFix = TRUE

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
      CHI(ig:ig) = ZERO
      
      IF (ig.LE.nchi .AND. Fxr(ifxr)%ldepl) CHI(ig:ig) = Fxr(ifxr)%chi(ig)
    ELSE
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1,j)
      itype     = Fxr(ifxr)%imix
      
      IF (nTracerCntl%lDynamicBen) THEN
        CHI(ig:ig) = GetChiDynBen(itype, TranInfo%fuelTemp(ipin, iz), ig, ig)
      ELSE
        CHI(ig:ig) = GetChiBen(itype, ig, ig)
      END IF
    END IF

    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      
      src1g(ifsr) = reigv * chi(ig) * psi(ifsr, iz)
    END DO
  END DO
  ! --------------------------------------------------
  !Scattering Source Update
  DO j = 1, nLocalFxr
    ifxr      = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    
    IF (lXsLib) Then
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
        src1g(ifsr) = src1g(ifsr) + xsmacsm(ig2, ig) * phis(ifsr, iz, ig2)
      END DO
      
      IF (lNegSrcFix .AND. src1g(ifsr) .LT. ZERO) THEN
        src1g(ifsr) = src1g(ifsr) - xsmacsm(ig, ig) * phis(ifsr, iz, ig)
        xst1g(ifsr) = xst1g(ifsr) - xsmacsm(ig, ig)
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
          
          IF (AxSrc1g(ipin).LT.0 .AND. .NOT.Fxr(ifxr)%lvoid) src1g(ifsr) = src1g(ifsr) - AxSrc1g(ipin)
        END DO
      END DO
    CASE(1)
      DO j = 1, nLocalFxr
        ifxr      = FxrIdxSt + j -1
        nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
        
        DO i = 1, nFsrInFxr
          ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
          
          src1g(ifsr) = src1g(ifsr) - AxSrc1g(ipin)
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
          psrc = psrc + Cellinfo(icel)%Vol(fsridx) * src1g(ifsr) ! pin-wise total source except axial source term
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
          
          src1g(ifsr) = src1g(ifsr) * (ONE - AxSrc1g(ipin) / psrc)
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
          psrc = psrc + Cellinfo(icel)%Vol(fsridx) * (ONE / 3._8 / xst1g(ifsr) * phis(ifsr, iz, ig)) ! pin-wise total source except axial source term
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
          
          src1g(ifsr) = src1g(ifsr) - AxSrc1g(ipin) * (ONE / 3._8 / xst1g(ifsr) * phis(ifsr, iz, ig)) / psrc
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
          psrc = psrc + Cellinfo(icel)%Vol(fsridx) * (ONE / 3._8 / xst1g(ifsr)) ! pin-wise total source except axial source term
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
          
          src1g(ifsr) = src1g(ifsr) - AxSrc1g(ipin) * (ONE / 3._8 / xst1g(ifsr)) / psrc
        END DO
      END DO
  END SELECT
#endif
  END IF
! --------------------------------------------------
  DO j = 1, CellInfo(icel)%nFsr
    ifsr = FsrIdxSt + j - 1
    
    src1g(ifsr) = src1g(ifsr) / xst1g(ifsr)
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
SUBROUTINE SetRtSrcGM_Cusping(Core, FmInfo, Fxr, src, phis, psi, axsrc, xstr1g, eigv, iz, ig, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix, PE)

USE PARAM
USE TYPEDEF,        ONLY : coreinfo_type,          Fxrinfo_type,          Cell_Type,      &
                           pin_Type,               GroupInfo_Type,        PE_TYPE,        &
                           XsMac_Type,             FmInfo_Type
USE BenchXs,        ONLY : GetChiBen,              xssben,                GetChiDynBen,    xssben_Cusping, &
                           MacXsBen,               xssDynBen,             xssDynben_Cusping, DynMacXsBen
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