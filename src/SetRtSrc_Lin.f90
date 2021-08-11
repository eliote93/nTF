#include <defines.h>
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE LinPsiUpdate(Core, Fxr, LinPsi, LinSrcSlope, myzb, myze, ng, lxslib, GroupInfo)

USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,       Fxrinfo_type,       Cell_Type,     pin_Type, &
                         GroupInfo_Type,      XsMac_Type
USE BenchXs,       ONLY : xsnfBen
USE MacXsLib_Mod, ONLY : MacXsNf
USE BasicOperation, ONLY : CP_CA, MULTI_VA
IMPLICIT NONE
TYPE(coreinfo_type) :: CORE
TYPE(Fxrinfo_type),POINTER :: Fxr(:, :)
TYPE(GroupInfo_Type) :: GroupInfo
REAL, POINTER :: LinPsi(:, :, :)
REAL, POINTER :: LinSrcSlope(:, :, :, :)
INTEGER :: myzb, myze, ng
LOGICAL :: lXsLib


TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Fxrinfo_type),POINTER :: myFxr
TYPE(XsMac_Type), SAVE :: XsMac

INTEGER :: nxy, nCoreFsr, nCoreFxr, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr
INTEGER :: ipin, icel, ifsrlocal, ifsr, ifxr, iz, itype, ig
INTEGER :: iResoGrpBeg, iResoGrpEnd, norg
INTEGER :: i, j, k

REAL, POINTER :: xsmacnf(:)

Pin => Core%Pin
CellInfo => Core%CellInfo; nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr; nxy = Core%nxy
IF(lxslib) THEN
  iResoGrpBeg = GroupInfo%nofg + 1
  iResoGrpEnd = GroupInfo%nofg + GroupInfo%norg
  norg = GroupInfo%norg
ENDIF

IF(.NOT. lxsLib) ALLOCATE(xsmacnf(ng))

DO iz = myzb, myze
  CALL CP_CA(LinPsi(:, :, iz), zero, 2, nCoreFsr)
  DO ipin = 1, nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
    icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr  
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)    
      myFxr => Fxr(ifxr, iz)
      IF(lXsLib) Then
        CALL MacXsNf(XsMac, myFxr, 1, ng, ng, 1._8, FALSE)
        xsmacnf => XsMac%XsMacNf
        IF(myFxr%lres) THEN
          do ig = iResoGrpBeg, iResoGrpEnd
            XsMacNf(ig) = XsMacNf(ig) * myFxr%fresoNF(ig) 
          enddo
        ENDIF
      ELSE
        ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1,j)
        !itype = CellInfo(icel)%iReg(ifsrlocal)      
        itype = myFxr%imix
        CALL xsnfben(itype, 1, ng, xsmacnf)
        !CHI(ig:ig) = GetChiBen(itype, ig, ig)
      ENDIF
      
      DO i = 1, nFsrInFxr
        ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1  !Global FSR Index
        DO ig = 1, ng
          LinPsi(1, ifsr, iz) = LinPsi(1, ifsr, iz) + xsmacnf(ig) * LinSrcSlope(1, ifsr, iz, ig)
          LinPsi(2, ifsr, iz) = LinPsi(2, ifsr, iz) + xsmacnf(ig) * LinSrcSlope(2, ifsr, iz, ig)
        ENDDO
        CONTINUE
        !src(ifsr) = reigv * chi(ig) * psic(ifsr, iz)
      ENDDO !Fsr Sweep    
    ENDDO
  ENDDO
ENDDO

IF(.NOT. lxsLib) THEN
  Deallocate(xsmacnf)
ELSE  
  NULLIFY(XsMaCnf)
ENDIF
NULLIFY(Pin)
NULLIFY(CellInfo)
NULLIFY(myFxr)

END SUBROUTINE LinPsiUpdate
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtLinSrc(Core, Fxr, RayInfo, src1g, LinSrc, LinPsi, LinSrcSlope, xstr1g, eigv, iz, ig, ng, GroupInfo, l3dim, lxslib, lscat1)

USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,          Fxrinfo_type,          Cell_Type,         &
                         RayInfo_Type,           AziAngleInfo_Type,     PolarAngle_Type,   &
                         pin_Type,               GroupInfo_Type,        XsMac_Type
USE BenchXs,      ONLY : GetChiBen,              xssben
USE MacXsLib_mod, ONLY : MacXsScatMatrix
USE BasicOperation, ONLY : CP_CA,                MULTI_VA
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type) :: Fxr(:)
TYPE(RayInfo_Type) :: RayInfo
TYPE(GroupInfo_Type):: GroupInfo
REAL, POINTER :: Src1g(:)
REAL, POINTER :: Linsrc(:, :)
REAL, POINTER :: LinPsi(:, :, :)
REAL, POINTER :: LinSrcSlope(:, :, :, :)
REAL, POINTER :: xstr1g(:)
REAL :: eigv
INTEGER :: myzb, myze, ig, ng, iz, ifsr, ifxr
LOGICAL :: lxslib, lscat1, l3dim

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(AziAngleInfo_Type), POINTER :: AziAng(:)
TYPE(XsMac_Type), SAVE :: XsMac
INTEGER :: nxy, nCoreFsr, nCoreFxr, FsrIdxSt, FxrIdxSt
INTEGER :: nlocalFxr, nFsrInFxr, nchi
INTEGER :: nAziAng, nPolAng
INTEGER :: ipin, icel, ifsrlocal, itype, ig2
INTEGER :: i, j, k, l

REAL :: reigv, sinv, cosv
REAL :: src0, segmax, sgn
REAL :: psrc(2, 5000)
REAL ::  chi(ng)
REAL, POINTER :: xsmacsm(:,:)

Pin => Core%Pin
CellInfo => Core%CellInfo
AziAng => RayInfo%AziAngle(:)
nCoreFsr = Core%nCoreFsr; nCoreFxr = Core%nCoreFxr
nAziAng = RayInfo%nAziAngle; nPolAng = RayInfo%nPolarAngle
nxy = Core%nxy
reigv = one/eigv


IF(.NOT. lxsLib) ALLOCATE(xsmacsm(ng, ng))

!CALL CP_CA(src(:), zero, nCoreFxr)
IF(.NOT. lxsLib) CALL CP_CA(xsmacsm, zero, ng, ng)
IF(lxsLib) nchi = GroupInfo%nchi
!Fission Source Updates
DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    !Get CHI
    IF(lXsLib) Then
      IF(ig .gt. nchi) THEN
        CHI(ig:ig) = 0
      ELSE
        CHI(ig:ig) = 0
        IF(Fxr(ifxr)%ldepl) CHI(ig:ig) = Fxr(ifxr)%chi(ig)
      ENDIF
    ELSE
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1,j)
      itype = CellInfo(icel)%iReg(ifsrlocal)      
      CHI(ig:ig) = GetChiBen(itype, ig, ig)
    ENDIF
    
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      !psrc = reigv * chi(ig) * psi(ifsr, iz)
      psrc(1, i) = reigv * chi(ig) * LinPsi(1, ifsr, iz)
      psrc(2, i) = reigv * chi(ig) * LinPsi(2, ifsr, iz)
    ENDDO !Fsr Sweep
    
    !Scattering Source
    IF(lXsLib) Then
      CALL MacXsScatMatrix(XsMac, Fxr(ifxr), ig, ig, ng, GroupInfo, FALSE)
      XsMacSm => XsMac%XsMacSm
    ELSE
      ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1,j)
      itype = Fxr(ifxr)%imix
      !itype = CellInfo(icel)%iReg(ifsrlocal)
      CALL xssben(itype, ig, 1, ng, XsMacsm, lscat1)
      !CHI(ig:ig) = GetChiBen(itype, ig, ig)
    ENDIF
    
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      DO ig2 = GroupInfo%InScatRange(1, ig), GroupInfo%InScatRange(2, ig)
      !DO ig2 = 1, ng
        psrc(1, i) = psrc(1, i) +  xsmacsm(ig2, ig) * LinSrcSlope(1, ifsr, iz, ig2)    !X-component
        psrc(2, i) = psrc(2, i) +  xsmacsm(ig2, ig) * LinSrcSlope(2, ifsr, iz, ig2)    !Y-component
      ENDDO
    ENDDO !Fsr Sweep
    
    DO l = 1, nAziAng
      sinv = AziAng(l)%sinv; cosv = AziAng(l)%cosv
      DO i = 1, nFsrInFxr
        ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
        LinSrc(ifsr, l) = psrc(1, i) * cosv + psrc(2, i) * sinv
      ENDDO    
    ENDDO
    CALL CP_CA(psrc, zero, 2, 5000)
  ENDDO !End of Fxr Sweep
ENDDO


DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt; icel = Pin(ipin)%Cell(iz);
  DO j = 1, CellInfo(icel)%nFsr
    ifsr = FsrIdxSt + j - 1
    DO l = 1, nAziAng
      src0 = LinSrc(ifsr, l)
      segmax = CellInfo(icel)%MaxSeg(l, j)
      IF(abs(src0) .LT. 2._8 * src1g(ifsr) / segmax) CYCLE
      sgn = src0/abs(src0)
      src0 = sgn * 2._8 * src1g(ifsr) / segmax
    ENDDO
    !LinSrc(ifsr, :) = LinSrc(ifsr, :) / (2._8 * xstr1g(ifsr) * xstr1g(ifsr))
  ENDDO
ENDDO


DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt; icel = Pin(ipin)%Cell(iz);
  DO j = 1, CellInfo(icel)%nFsr
    ifsr = FsrIdxSt + j - 1
    LinSrc(ifsr, :) = LinSrc(ifsr, :) / (2._8 * xstr1g(ifsr) * xstr1g(ifsr))
  ENDDO
ENDDO
continue
IF(.NOT. lxsLib) DEALLOCATE(Xsmacsm)
IF(lxsLib) THEN
  NULLIFY(Xsmacsm)
ENDIF
NULLIFY(Pin)
NULLIFY(CellInfo)
NULLIFY(AziAng)

END SUBROUTINE SetRtLinSrc
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE LinPsiUpdate_CASMO(Core, Fxr, phisSlope, psiSlope, myzb, myze, ng, lxslib, GroupInfo)

USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,       Fxrinfo_type,       Cell_Type,     pin_Type, &
                         GroupInfo_Type,      XsMac_Type
USE BenchXs,       ONLY : xsnfBen
USE MacXsLib_Mod, ONLY : MacXsNf
USE BasicOperation, ONLY : CP_CA, MULTI_VA
IMPLICIT NONE
TYPE(coreinfo_type) :: CORE
TYPE(Fxrinfo_type),POINTER :: Fxr(:, :)
TYPE(GroupInfo_Type) :: GroupInfo
REAL, POINTER :: phisSlope(:, :, :, :)
REAL, POINTER :: psiSlope(:, :, :)
INTEGER :: myzb, myze, ng
LOGICAL :: lXsLib

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Fxrinfo_type),POINTER :: myFxr
TYPE(XsMac_Type), SAVE :: XsMac

INTEGER :: nxy, nCoreFsr, nCoreFxr, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr
INTEGER :: ipin, icel, ifsrlocal, ifsr, ifxr, iz, itype, ig
INTEGER :: iResoGrpBeg, iResoGrpEnd, norg
INTEGER :: i, j, k

REAL, POINTER :: xsmacnf(:)

Pin => Core%Pin
CellInfo => Core%CellInfo; nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr; nxy = Core%nxy
IF(lxslib) THEN
  iResoGrpBeg = GroupInfo%nofg + 1
  iResoGrpEnd = GroupInfo%nofg + GroupInfo%norg
  norg = GroupInfo%norg
ENDIF

IF(.NOT. lxsLib) ALLOCATE(xsmacnf(ng))

DO iz = myzb, myze
  psiSlope(:, :, iz) = zero
  DO ipin = 1, nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
    icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr  
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)    
      myFxr => Fxr(ifxr, iz)
      IF(lXsLib) Then
        CALL MacXsNf(XsMac, myFxr, 1, ng, ng, 1._8, FALSE, TRUE)
        xsmacnf => XsMac%XsMacNf
        IF(myFxr%lres) THEN
          do ig = iResoGrpBeg, iResoGrpEnd
            XsMacNf(ig) = XsMacNf(ig) * myFxr%fresoNF(ig)  
          enddo
        ENDIF
      ELSE
        ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1,j)
        itype = myFxr%imix
        CALL xsnfben(itype, 1, ng, xsmacnf)
      ENDIF
      
      DO i = 1, nFsrInFxr
        ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
        DO ig = 1, ng
          psiSlope(:, ig, iz) = psiSlope(:, ig, iz) + xsmacnf(ig) * phisSlope(:, ig, ifsr, iz)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO

IF(.NOT. lxsLib) Deallocate(xsmacnf)
IF(lXsLib) NULLIFY(XsMacNf)
NULLIFY(Pin, CellInfo)

END SUBROUTINE LinPsiUpdate_CASMO
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRtLinSrc_CASMO(Core, Fxr, RayInfo, phisNg, phisSlope, srcNg, srcSlope, psi, psiSlope, axsrc, xstNg, eigv, iz, gb, ge, ng, GroupInfo, l3dim, lxslib, lscat1, lNegFix)

USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,          Fxrinfo_type,          Cell_Type,              &
                         RayInfo_Type,           AziAngleInfo_Type,     PolarAngle_Type,        &
                         pin_Type,               GroupInfo_Type,        XsMac_Type
USE BenchXs,      ONLY : GetChiBen,              xssben
USE MacXsLib_mod, ONLY : MacXsScatMatrix
USE BasicOperation, ONLY : CP_CA,                MULTI_VA
USE PE_Mod,       ONLY : PE
USE OMP_LIB
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(RayInfo_Type) :: RayInfo
TYPE(GroupInfo_Type):: GroupInfo
REAL, POINTER :: phisNg(:, :), phisSlope(:, :, :, :)
REAL, POINTER :: srcNg(:, :), srcSlope(:, :, :, :), axsrc(:, :, :)
REAL, POINTER :: psi(:, :), psiSlope(:, :, :)
REAL, POINTER :: xstNg(:, :)
REAL :: eigv
INTEGER :: myzb, myze, gb, ge, ng, iz, ifsr, ifxr
LOGICAL :: lxslib, lscat1, l3dim, lNegFix

TYPE(FxrInfo_Type), POINTER :: myFxr
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(XsMac_Type), SAVE :: XsMac(nTHREADMAX)
INTEGER :: nCoreFsr, nCoreFxr, nxy, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, nchi
INTEGER :: ipin, icel, ifsrlocal, itype, ig, ig2, tid
INTEGER :: i, j, k
INTEGER :: xyb, xye   !--- CNJ Edit : Domain Decomposition + MPI

REAL :: reigv
REAL :: chi(ng)
REAL :: AxSrc1Pin(ng)
REAL, POINTER :: xsmacsm(:,:)

LOGICAL :: lscat1sum
LOGICAL :: lNegSrcFix
Pin => Core%Pin
CellInfo => Core%CellInfo
nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr
nxy = Core%nxy
xyb = PE%myPinBeg; xye = PE%myPinEnd

reigv = one/eigv
lscat1sum = lscat1

lNegSrcFix = FALSE
IF(lNegFix) lNegSrcFix = TRUE

IF(.NOT. lxsLib) THEN
  DO i = 1, PE%nThread
    ALLOCATE(XsMac(i)%XsMacSm(ng, ng))
  ENDDO
ENDIF

IF(lxsLib) nchi = GroupInfo%nchi

CALL omp_set_num_threads(PE%nThread)

!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(i, j, k, ifsr, ifxr, ipin, icel, ifsrlocal, itype, ig, ig2, tid, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr, myFxr, xsmacsm, chi, AxSrc1Pin)
tid = omp_get_thread_num() + 1
!$OMP DO SCHEDULE(GUIDED)
DO ipin = 1, nxy
  FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
  icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr
  AxSrc1Pin(gb : ge) = AxSrc(ipin, iz, gb : ge)
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j -1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)  
    myFxr => Fxr(ifxr, iz)
    DO ig = gb, ge
      IF(lXsLib) Then
        IF(ig .gt. nchi) THEN
          CHI(ig) = 0
        ELSE
          CHI(ig) = 0
          IF(myFxr%ldepl) CHI(ig) = myFxr%chi(ig)
        ENDIF
      ELSE
        ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1,j)    
        itype = Fxr(ifxr, iz)%imix
        CHI(ig:ig) = GetChiBen(itype, ig, ig)
      ENDIF
    ENDDO
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      DO ig = gb, ge
        srcNg(ig, ifsr) = reigv * chi(ig) * psi(ifsr, iz)
        srcSlope(1 : 2, ig, ifsr, iz) = reigv * chi(ig) * psiSlope(:, ifsr, iz)
        srcSlope(3 : 4, ig, ifsr, iz) = reigv * chi(ig) * psiSlope(:, ifsr, iz)
      ENDDO
    ENDDO
  ENDDO
  DO j = 1, nLocalFxr
    ifxr = FxrIdxSt + j - 1
    nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
    myFxr => FXR(ifxr, iz)
    DO ig = gb, ge
      IF(lXsLib) Then
        !CALL MacXsScatMatrix(XsMac(tid), myFxr, ig, ig, ng, GroupInfo, lscat1, lscat1sum, .TRUE.)
        CALL MacXsScatMatrix(XsMac(tid), myFxr, ig, ig, ng, GroupInfo, lscat1, .TRUE.)
        XsMacSm => XsMac(tid)%XsMacSm
#ifdef inflow      
        XsMacSm(ig, ig) = XsMacSm(ig, ig) + Fxr(ifxr, iz)%DelInflow(ig)
#endif      
      ELSE
        ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1,j)
        itype = myFxr%imix
        XsMacSm => XsMac(tid)%XsMacSm
        CALL xssben(itype, ig, 1, ng, XsMacsm, lscat1)
      ENDIF
    ENDDO
    DO i = 1, nFsrInFxr
      ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
      DO ig = gb, ge
        DO ig2 = GroupInfo%InScatRange(1, ig), GroupInfo%InScatRange(2, ig)
          srcNg(ig, ifsr) = srcNg(ig, ifsr) + xsmacsm(ig2, ig) * phisNg(ig2, ifsr)
          srcSlope(1 : 2, ig, ifsr, iz) = srcSlope(1 : 2, ig, ifsr, iz) + xsmacsm(ig2, ig) * phisSlope(:, ig2, ifsr, iz)
          srcSlope(3 : 4, ig, ifsr, iz) = srcSlope(3 : 4, ig, ifsr, iz) + xsmacsm(ig2, ig) * phisSlope(:, ig2, ifsr, iz)
        ENDDO
        IF (lNegSrcFix .AND. srcNg(ig, ifsr) .LT. 0._8) THEN
          srcNg(ig, ifsr) = srcNg(ig, ifsr) - xsmacsm(ig, ig) * phisNg(ig, ifsr)
          xstNg(ig, ifsr) = xstNg(ig, ifsr) - xsmacsm(ig, ig)
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  IF(l3dim) THEN
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)
      DO i = 1, nFsrInFxr
        ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1
        DO ig = gb, ge
#ifndef LkgSplit
          srcNg(ig, ifsr) = srcNg(ig, ifsr) - AxSrc1Pin(ig)
#else
          IF(AxSrc1Pin(ig) .LT. 0 .AND. .NOT. Fxr(ifxr, iz)%lvoid) THEN
            srcNg(ig, ifsr) = srcNg(ig, ifsr) - AxSrc1Pin(ig)
          ENDIF
#endif       
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  DO j = 1, CellInfo(icel)%nFsr
    ifsr = FsrIdxSt + j - 1
    DO ig = gb, ge
      srcNg(ig, ifsr) = srcNg(ig, ifsr) / xstNg(ig, ifsr)
      srcSlope(1 : 2, ig, ifsr, iz) = srcSlope(1 : 2, ig, ifsr, iz) / xstNg(ig, ifsr)
      srcSlope(3 : 4, ig, ifsr, iz) = srcSlope(3 : 4, ig, ifsr, iz) / (xstNg(ig, ifsr) ** 2)
    ENDDO
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

NULLIFY(myFxr)
NULLIFY(Pin)
NULLIFY(CellInfo)

END SUBROUTINE SetRtLinSrc_CASMO
! ------------------------------------------------------------------------------------------------------------