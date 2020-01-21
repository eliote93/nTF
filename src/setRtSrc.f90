#include <defines.h>
SUBROUTINE SetRtSrc(Core, Fxr, src, phis, psi, axsrc, xstr1g, &
                    eigv, iz, ig, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix, PE)
USE PARAM
USE TYPEDEF,        ONLY : coreinfo_type,          Fxrinfo_type,          Cell_Type,      &
                           pin_Type,               GroupInfo_Type,        PE_TYPE,        &
                           XsMac_Type
USE BenchXs,        ONLY : GetChiBen,              xssben
USE MacXsLib_mod,   ONLY : MacXsScatMatrix
USE BasicOperation, ONLY : CP_CA
USE CNTL,             ONLY : nTracerCntl
USE OMP_LIB
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
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
INTEGER :: i, j, k

REAL :: reigv
REAL ::  chi(ng)
REAL, POINTER :: xsmacsm(:,:)

REAL :: psrc, pvol

LOGICAL :: lNegSrcFix
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
!$OMP PRIVATE(i, j, k, ifsr, ifxr, ipin, icel, ifsrlocal, itype, ig2, tid, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, xsmacsm, chi, psrc, pvol, fsridx)
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
      !itype = CellInfo(icel)%iReg(ifsrlocal)
      itype = Fxr(ifxr)%imix
      CHI(ig:ig) = GetChiBen(itype, ig, ig)
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
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1,j)
      !itype = CellInfo(icel)%iReg(ifsrlocal)
      itype = Fxr(ifxr)%imix
      XsMacSm => XsMac(tid)%XsMacSm
      CALL xssben(itype, ig, 1, ng, XsMacsm, lscat1)
      !CHI(ig:ig) = GetChiBen(itype, ig, ig)
    ENDIF

    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      DO ig2 = GroupInfo%InScatRange(1, ig), GroupInfo%InScatRange(2, ig)
      !DO ig2 = 1, ng
        src(ifsr) = src(ifsr) + xsmacsm(ig2, ig)*phis(ifsr, iz, ig2)
      ENDDO
      IF(lNegSrcFix .AND. src(ifsr) .LT. 0._8) THEN
        src(ifsr) = src(ifsr) - xsmacsm(ig, ig)*phis(ifsr, iz, ig)
        xstr1g(ifsr) = xstr1g(ifsr) - xsmacsm(ig, ig)
      ENDIF
    ENDDO !Fsr Sweep

  ENDDO !End of Fxr Sweep

!--- Backup
!  IF(l3dim) THEN
!    DO j = 1, nLocalFxr
!      ifxr = FxrIdxSt + j -1
!      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
!      DO i = 1, nFsrInFxr
!        ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
!#ifndef LkgSplit
!        src(ifsr) = src(ifsr) - AxSrc(ipin)
!#else
!        IF(AxSrc(ipin) .LT. 0 .AND. .NOT. Fxr(ifxr)%lvoid) THEN
!          src(ifsr) = src(ifsr) - AxSrc(ipin)
!        ENDIF
!#endif
!        if(src(ifsr) .LT. 0)THEN
!            !WRITE(99,*) 'Negative Src at ig, iz, ifsr ', ig, iz, ifsr
!            !WRITE(*,*)  'Negative Src at ig, iz, ifsr ', ig, iz, ifsr
!        ENDIF
!      ENDDO
!    ENDDO
!  ENDIF
!--- Backup End

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

END SUBROUTINE

!--- CNJ Edit : Node Majors
SUBROUTINE SetRtSrcNM(Core, Fxr, srcnm, phisnm, psi, AxSrc, xstnm, eigv, iz,              &
                      gb, ge, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix, PE, Offset)
USE PARAM
USE TYPEDEF,        ONLY : coreinfo_type,          Fxrinfo_type,          Cell_Type,      &
                           pin_Type,               GroupInfo_Type,        PE_TYPE,        &
                           XsMac_Type
USE BenchXs,        ONLY : GetChiBen,              xssben
USE MacXsLib_mod,   ONLY : MacXsScatMatrix
USE CNTL,           ONLY : nTracerCntl
USE OMP_LIB
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
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

REAL :: reigv
REAL ::  chi(ng)
REAL, POINTER :: xsmacsm(:,:)

LOGICAL :: lNegSrcFix
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
!$OMP         FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, xsmacsm, chi)
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
      CHI(gb : ge) = GetChiBen(itype, gb, ge)
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
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1,j)
      itype = Fxr(ifxr)%imix
      XsMacSm => XsMac(tid)%XsMacSm
      DO ig = gb, ge
        CALL xssben(itype, ig, 1, ng, XsMacsm, lscat1)
      ENDDO
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
END SUBROUTINE

SUBROUTINE PseudoAbsorption(Core, Fxr, src, phis1g, AxPXS, xstr1g, iz, ig, ng, GroupInfo, l3dim)
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,          Fxrinfo_type,          Cell_Type,      &
                         pin_Type,               GroupInfo_Type,        XsMac_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:)
TYPE(GroupInfo_Type):: GroupInfo
REAL :: phis1g(:), AxPXS(:)
REAL, POINTER :: src(:), xstr1g(:)
REAL :: eigv
INTEGER :: myzb, myze, ig, ng, iz
LOGICAL :: l3dim

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
REAL :: pAbXs, phiavg, vol, pSrc
INTEGER :: nxy, nCoreFsr, nCoreFxr
INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nFsrInFxr
INTEGER :: i, j, ipin, icel, ifsr, ifxr

IF(.NOT. l3dim) RETURN

Pin => Core%Pin
CellInfo => Core%CellInfo
nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr
nxy = Core%nxy
DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt; icel = Pin(ipin)%Cell(iz);
  FxrIdxSt = Pin(ipin)%FxrIdxSt; nLocalFxr = CellInfo(icel)%nFxr
  pSrc=0
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      pSrc=pSrc+ xstr1g(ifsr)
      IF(.NOT. Fxr(ifxr)%lVoid) xstr1g(ifsr) = xstr1g(ifsr)+AxPXS(ipin)
    ENDDO
  ENDDO
ENDDO
NULLIFY(Pin)
NULLIFY(CellInfo)
END SUBROUTINE

!--- CNJ Edit : Node Majors
SUBROUTINE PseudoAbsorptionNM(Core, Fxr, AxPXS, xstnm, iz, ng, GroupInfo, l3dim)
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,          Fxrinfo_type,          Cell_Type,      &
                         pin_Type,               GroupInfo_Type,        XsMac_Type
USE PE_MOD,       ONLY : PE
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:)
TYPE(GroupInfo_Type):: GroupInfo
REAL, POINTER :: AxPXS(:, :, :), xstnm(:, :)
REAL :: eigv
INTEGER :: iz, ng
LOGICAL :: l3dim

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
REAL :: pAbXs, phiavg, vol
INTEGER :: nCoreFsr, nCoreFxr, nxy
INTEGER :: FsrIdxSt, FxrIdxSt, nLocalFxr, nFsrInFxr
INTEGER :: i, j, ipin, icel, ifsr, ifxr, ig
INTEGER :: xyb, xye   !--- CNJ Edit : Domain Decomposition + MPI

IF(.NOT. l3dim) RETURN

Pin => Core%Pin
CellInfo => Core%CellInfo
nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr
nxy = Core%nxy
xyb = PE%myPinBeg; xye = PE%myPinEnd   !--- CNJ Edit : Domain Decomposition + MPI

DO ipin = xyb, xye
  FsrIdxSt = Pin(ipin)%FsrIdxSt; icel = Pin(ipin)%Cell(iz);
  FxrIdxSt = Pin(ipin)%FxrIdxSt; nLocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      DO ig = 1, ng
        IF(.NOT. Fxr(ifxr)%lVoid) xstnm(ig, ifsr) = xstnm(ig, ifsr) + AxPXS(ipin, iz, ig)
      ENDDO
    ENDDO
  ENDDO
ENDDO

NULLIFY(Pin)
NULLIFY(CellInfo)
END SUBROUTINE

SUBROUTINE SetRtP1Src(Core, Fxr, srcm, phim, xstr1g, iz, ig, ng, GroupInfo, l3dim, lxslib, lscat1, ScatOd, PE)
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,          Fxrinfo_type,          Cell_Type,      &
                         pin_Type,               GroupInfo_Type,        XsMac_Type,     &
                         PE_Type
USE BenchXs,      ONLY : GetChiBen,              xssben,                                &
                         xssm1ben,               xssm2ben,              xssm3ben
USE MacXsLib_mod, ONLY : MacP1XsScatMatrix,      MacP2XsScatMatrix,     MacP3XsScatMatrix
USE XsUtil_mod,   ONLY : GetXsMacDat,            ReturnXsMacDat,        FreeXsMac
USE BasicOperation, ONLY : CP_CA
USE OMP_LIB
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:)
TYPE(GroupInfo_Type):: GroupInfo
REAL, POINTER :: srcm(:, :)
REAL, POINTER :: phim(:, :, :, :)
REAL, POINTER :: xstr1g(:)
REAL :: eigv
INTEGER :: myzb, myze, ig, ng, iz, ifsr, ifxr
LOGICAL :: lxslib, lscat1, l3dim
INTEGER :: ScatOd
TYPE(PE_Type) :: PE

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(XsMac_Type), POINTER :: XsMac
INTEGER :: nxy, nCoreFsr, nCoreFxr, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, nchi
INTEGER :: ipin, icel, ifsrlocal, itype, ig2, tid
INTEGER :: i, j, k

REAL :: reigv
REAL ::  chi(ng)
REAL, POINTER :: XsMacP1sm(:,:), XsMacP2sm(:,:), XsMacP3sm(:,:)

Pin => Core%Pin
CellInfo => Core%CellInfo
nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr
nxy = Core%nxy

reigv = one/eigv

IF (ScatOd .EQ. 1) THEN
  CALL CP_CA(srcm(:, :), zero, 2, nCoreFsr)
ELSEIF (ScatOd .EQ. 2) THEN
  CALL CP_CA(srcm(:, :), zero, 5, nCoreFsr)
ELSEIF (ScatOd .EQ. 3) THEN
  CALL CP_CA(srcm(:, :), zero, 9, nCoreFsr)
ENDIF

IF (lxsLib) nchi = GroupInfo%nchi
!Fission Source Updates
tid = 1
!Scattering Source Update
!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(PE%nThread)
!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(XsMac, i, j, k, ifsr, ifxr, ipin, icel, ifsrlocal, itype, ig2, tid, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, XsMacP1Sm, XsMacP2Sm, XsMacP3Sm)
!$  tid = omp_get_thread_num()+1
!Scattering Source Update

IF (.NOT. lxsLib) THEN
  ALLOCATE(XsMacP1sm(ng, ng))
  IF(ScatOd .GE. 2) ALLOCATE(XsMacP2sm(ng, ng))
  IF(ScatOd .EQ. 3) ALLOCATE(XsMacP3sm(ng, ng))
  CALL CP_CA(XsMacP1sm, zero, ng, ng)
  IF(ScatOd .GE. 2) CALL CP_CA(XsMacP2sm, zero, ng, ng)
  IF(ScatOd .EQ. 3) CALL CP_CA(XsMacP3sm, zero, ng, ng)
ENDIF

!$OMP DO
DO ipin = 1, nxy
  CALL GetXsMacDat(XsMac, ng, TRUE)
  FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    IF(lXsLib) Then
      CALL MacP1XsScatMatrix(XsMac, Fxr(ifxr), ig, ig, ng, GroupInfo)
      IF(ScatOd .GE. 2) CALL MacP2XsScatMatrix(XsMac, Fxr(ifxr), ig, ig, ng, GroupInfo)
      IF(ScatOd .EQ. 3) CALL MacP3XsScatMatrix(XsMac, Fxr(ifxr), ig, ig, ng, GroupInfo)
      XsMacP1Sm => XsMac%XsMacP1Sm; XsMacP2Sm => XsMac%XsMacP2Sm; XsMacP3Sm => XsMac%XsMacP3Sm
    ELSE
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1,j)
      itype = Fxr(ifxr)%imix
      CALL XsSm1Ben(itype, 1, ng, 1, ng, XsMacP1Sm)
      IF(ScatOd .GE. 2) CALL XsSm2Ben(itype, 1, ng, 1, ng, XsMacP2Sm)   !!  OPTIMZIE
      IF(ScatOd .EQ. 3) CALL XsSm3Ben(itype, 1, ng, 1, ng, XsMacP3Sm)
    ENDIF

    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      !DO ig2 = GroupInfo%InScatRange(1, ig), GroupInfo%InScatRange(2, ig)
      IF(ScatOd .EQ. 1) THEN
        DO ig2 = 1, ng
           srcm(1, ifsr) = srcm(1, ifsr) + XsMacP1Sm(ig2, ig) * phim(1, ifsr, iz, ig2)
           srcm(2, ifsr) = srcm(2, ifsr) + XsMacP1Sm(ig2, ig) * phim(2, ifsr, iz, ig2)
        ENDDO
      ELSEIF(ScatOd .EQ. 2) THEN
        DO ig2 = 1, ng
          srcm(1:2, ifsr) = srcm(1:2, ifsr) + XsMacP1Sm(ig2, ig) * phim(1:2, ifsr, iz, ig2)
          srcm(3:5, ifsr) = srcm(3:5, ifsr) + XsMacP2Sm(ig2, ig) * phim(3:5, ifsr, iz, ig2)
        ENDDO
      ELSEIF(ScatOd .EQ. 3) THEN
        DO ig2 = 1, ng
          srcm(1:2, ifsr) = srcm(1:2, ifsr) + XsMacP1Sm(ig2, ig) * phim(1:2, ifsr, iz, ig2)
          srcm(3:5, ifsr) = srcm(3:5, ifsr) + XsMacP2Sm(ig2, ig) * phim(3:5, ifsr, iz, ig2)
          srcm(6:9, ifsr) = srcm(6:9, ifsr) + XsMacP3Sm(ig2, ig) * phim(6:9, ifsr, iz, ig2)
        ENDDO
      ENDIF

    ENDDO !Fsr Sweep
  ENDDO !End of Fxr Sweep
  CALL ReturnXsMacDat(XsMac)   !Memory leak problem !! 17/01/09   big memory but stable
ENDDO !End of Pin
!$OMP END DO

IF(.NOT. lxsLib) DEALLOCATE(XsMacP1sm) ! modified because of crash! in benchmark XS
IF(.NOT. lxsLib .AND. ScatOd .GE. 2) DEALLOCATE(XsMacP2sm)
IF(.NOT. lxsLib .AND. ScatOd .EQ. 3) DEALLOCATE(XsMacP3sm)
!$OMP END PARALLEL
!Axail Source Contribution
IF(ScatOd .EQ. 1) THEN
  DO ipin = 1, nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt; icel = Pin(ipin)%Cell(iz);
    DO j = 1, CellInfo(icel)%nFsr
      ifsr = FsrIdxSt + j - 1
      srcm(1:2, ifsr) = 3._8 * srcm(1:2, ifsr) / xstr1g(ifsr)
    ENDDO
  ENDDO
ELSEIF(ScatOd .EQ. 2) THEN
  DO ipin = 1, nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt; icel = Pin(ipin)%Cell(iz);
    DO j = 1, CellInfo(icel)%nFsr
      ifsr = FsrIdxSt + j - 1
      srcm(1:2, ifsr) = 3._8 * srcm(1:2, ifsr) / xstr1g(ifsr)
      srcm(3:5, ifsr) = 5._8 * srcm(3:5, ifsr) / xstr1g(ifsr)
    ENDDO
  ENDDO
ELSEIF(ScatOd .EQ. 3) THEN
  DO ipin = 1, nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt; icel = Pin(ipin)%Cell(iz);
    DO j = 1, CellInfo(icel)%nFsr
      ifsr = FsrIdxSt + j - 1
      srcm(1:2, ifsr) = 3._8 * srcm(1:2, ifsr) / xstr1g(ifsr)
      srcm(3:5, ifsr) = 5._8 * srcm(3:5, ifsr) / xstr1g(ifsr)
      srcm(6:9, ifsr) = 7._8 * srcm(6:9, ifsr) / xstr1g(ifsr)
    ENDDO
  ENDDO
ENDIF

NULLIFY(XsMacP1Sm)
NULLIFY(Pin)
NULLIFY(CellInfo)

END SUBROUTINE

!--- CNJ Edit : Node Majors
SUBROUTINE SetRtP1SrcNM(Core, Fxr, srcmnm, phimnm, xstnm, iz, gb, ge, ng, GroupInfo, lxslib, ScatOd, PE, Offset)
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,          Fxrinfo_type,          Cell_Type,      &
                         pin_Type,               GroupInfo_Type,        XsMac_Type,     &
                         PE_Type
USE BenchXs,      ONLY : GetChiBen,              xssben
USE MacXsLib_mod, ONLY : MacP1XsScatMatrix,      MacP2XsScatMatrix,     MacP3XsScatMatrix
USE XsUtil_mod,   ONLY : GetXsMacDat,            ReturnXsMacDat
USE BasicOperation, ONLY : CP_CA
USE OMP_LIB
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:)
TYPE(GroupInfo_Type):: GroupInfo
REAL, POINTER :: srcmnm(:, :, :)
REAL, POINTER :: phimnm(:, :, :)
REAL, POINTER :: xstnm(:, :)
REAL :: eigv
INTEGER :: gb, ge, ng, iz, ifsr, ifxr
LOGICAL :: lxslib, lscat1, l3dim
INTEGER :: ScatOd
TYPE(PE_Type) :: PE
INTEGER, OPTIONAL :: Offset   !--- CNJ Edit : GPU Acceleration

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(XsMac_Type), POINTER :: XsMac
INTEGER :: nCoreFsr, nCoreFxr, nxy, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, nchi
INTEGER :: ipin, icel, ifsrlocal, itype, ig, ig2
INTEGER :: i, j, k
INTEGER :: xyb, xye   !--- CNJ Edit : Domain Decomposition + MPI
INTEGER :: off        !--- CNJ Edit : GPU Acceleration
REAL :: reigv
REAL :: chi(ng)
REAL :: scatSrc(9), xstinv
REAL, POINTER :: xsmacP1sm(:,:), xsmacP2sm(:,:), xsmacP3sm(:,:)

Pin => Core%Pin
CellInfo => Core%CellInfo
nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr
nxy = Core%nxy
xyb = PE%myPinBeg; xye = PE%myPinEnd   !--- CNJ Edit : Domain Decomposition + MPI
off = 0
IF (PRESENT(Offset)) off = Offset
reigv = one/eigv

IF(.NOT. lxsLib) ALLOCATE(xsmacP1sm(ng, ng))

srcmnm = zero

IF(.NOT. lxsLib) CALL CP_CA(XsMacP1sm, zero, ng, ng)
IF(lxsLib) nchi = GroupInfo%nchi

!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(XsMac, ifsr, ifxr, icel, ifsrlocal, itype, scatSrc, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, XsMacP1Sm, XsMacP2Sm, XsMacP3Sm)
!$OMP DO
DO ipin = xyb, xye
  CALL GetXsMacDat(XsMac, ng, TRUE)
  FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    IF(lXsLib) Then
      CALL MacP1XsScatMatrix(XsMac, Fxr(ifxr), gb, ge, ng, GroupInfo)
      IF(ScatOd .GE. 2) CALL MacP2XsScatMatrix(XsMac, Fxr(ifxr), gb, ge, ng, GroupInfo)
      IF(ScatOd .EQ. 3) CALL MacP3XsScatMatrix(XsMac, Fxr(ifxr), gb, ge, ng, GroupInfo)
      XsMacP1Sm => XsMac%XsMacP1Sm; XsMacP2Sm => XsMac%XsMacP2Sm; XsMacP3Sm => XsMac%XsMacP3Sm
    ELSE
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1, j)
      itype = Fxr(ifxr)%imix
    ENDIF
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      DO ig = gb, ge
        IF(ScatOd .EQ. 1) THEN
          scatSrc = 0
          DO ig2 = 1, ng
            scatSrc(1:2) = scatSrc(1:2) + XsMacP1Sm(ig2, ig) * phimnm(1:2, ig2, ifsr)
          ENDDO
          srcmnm(1:2, ig - off, ifsr) = scatSrc(1:2)
        ELSEIF(ScatOd .EQ. 2) THEN
          scatSrc = 0
          DO ig2 = 1, ng
            scatSrc(1:2) = scatSrc(1:2) + XsMacP1Sm(ig2, ig) * phimnm(1:2, ig2, ifsr)
            scatSrc(3:5) = scatSrc(3:5) + XsMacP2Sm(ig2, ig) * phimnm(3:5, ig2, ifsr)
          ENDDO
          srcmnm(1:5, ig - off, ifsr) = scatSrc(1:5)
        ELSEIF(ScatOd .EQ. 3) THEN
          scatSrc = 0
          DO ig2 = 1, ng
            scatSrc(1:2) = scatSrc(1:2) + XsMacP1Sm(ig2, ig) * phimnm(1:2, ig2, ifsr)
            scatSrc(3:5) = scatSrc(3:5) + XsMacP2Sm(ig2, ig) * phimnm(3:5, ig2, ifsr)
            scatSrc(6:9) = scatSrc(6:9) + XsMacP3Sm(ig2, ig) * phimnm(6:9, ig2, ifsr)
          ENDDO
          srcmnm(1:9, ig - off, ifsr) = scatSrc(1:9)
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  CALL ReturnXsMacDat(XsMac)
ENDDO
!$OMP END DO
!$OMP END PARALLEL

IF(ScatOd .EQ. 1) THEN
  DO ipin = xyb, xye
    FsrIdxSt = Pin(ipin)%FsrIdxSt; icel = Pin(ipin)%Cell(iz);
    DO j = 1, CellInfo(icel)%nFsr
      ifsr = FsrIdxSt + j - 1
      DO ig = gb, ge
        xstinv = 1._8 / xstnm(ig, ifsr)
        srcmnm(1:2, ig - off, ifsr) = 3._8 * srcmnm(1:2, ig - off, ifsr) * xstinv
      ENDDO
    ENDDO
  ENDDO
ELSEIF(ScatOd .EQ. 2) THEN
  DO ipin = xyb, xye
    FsrIdxSt = Pin(ipin)%FsrIdxSt; icel = Pin(ipin)%Cell(iz);
    DO j = 1, CellInfo(icel)%nFsr
      ifsr = FsrIdxSt + j - 1
      DO ig = gb, ge
        xstinv = 1._8 / xstnm(ig, ifsr)
        srcmnm(1:2, ig - off, ifsr) = 3._8 * srcmnm(1:2, ig - off, ifsr) * xstinv
        srcmnm(3:5, ig - off, ifsr) = 5._8 * srcmnm(3:5, ig - off, ifsr) * xstinv
      ENDDO
    ENDDO
  ENDDO
ELSEIF(ScatOd .EQ. 3) THEN
  DO ipin = xyb, xye
    FsrIdxSt = Pin(ipin)%FsrIdxSt; icel = Pin(ipin)%Cell(iz);
    DO j = 1, CellInfo(icel)%nFsr
      ifsr = FsrIdxSt + j - 1
      DO ig = gb, ge
        xstinv = 1._8 / xstnm(ig, ifsr)
        srcmnm(1:2, ig - off, ifsr) = 3._8 * srcmnm(1:2, ig - off, ifsr) * xstinv
        srcmnm(3:5, ig - off, ifsr) = 5._8 * srcmnm(3:5, ig - off, ifsr) * xstinv
        srcmnm(6:9, ig - off, ifsr) = 7._8 * srcmnm(6:9, ig - off, ifsr) * xstinv
      ENDDO
    ENDDO
  ENDDO
ENDIF

IF(.NOT. lxsLib) DEALLOCATE(XsMacP1Sm)
IF(lxsLib) NULLIFY(XsMacP1Sm)
NULLIFY(Pin)
NULLIFY(CellInfo)
END SUBROUTINE

SUBROUTINE AddConstSrc(Core, Fxr, Src, xstr1g, ConstSrc, iz, ig, ng)
USE PARAM
USE TYPEDEF,       ONLY : CoreInfo_Type,        Fxrinfo_type,                          &
                          Cell_Type,            pin_Type
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:)
REAL, POINTER :: Src(:), xstr1g(:)
REAL :: ConstSrc
INTEGER :: iz, ig, ng



TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)

INTEGER :: ipin, ifsr, icel
INTEGER :: FsrIdxSt
INTEGER :: i, j, k
Pin => Core%Pin
CellInfo => Core%CellInfo
nxy = Core%nxy
DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt; icel = Pin(ipin)%Cell(iz);
  DO j = 1, CellInfo(icel)%nFsr
    ifsr = FsrIdxSt + j - 1
    !IF(src(ifsr) .lt. 0) THEN
    !  src(ifsr) = 0
    !ENDIF
    src(ifsr) = src(ifsr) + ConstSrc/xstr1g(ifsr)
  ENDDO
ENDDO
END SUBROUTINE
