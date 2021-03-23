!#include <defines.h>
!#include <CUDADEFINES.h>
!
!#ifdef __PGI
!MODULE CUDA_NNFSP
!USE CUDA_MASTER
!USE CUDA_UTIL
!USE CUDA_SYSTEM
!USE CUDA_SOLVER
!IMPLICIT NONE
!
!TYPE FMNN_Type
!  REAL, ALLOCATABLE :: phis(:,:,:)
!  REAL, ALLOCATABLE :: Jout(:,:,:,:)
!  REAL, ALLOCATABLE :: nnsrc(:,:,:)
!  REAL, ALLOCATABLE :: psi(:,:)
!  REAL, ALLOCATABLE :: AxSrc(:,:,:), AxPxs(:,:,:)
!END TYPE
!
!TYPE CMNN_Type
!END TYPE
!
!TYPE(FMNN_Type) :: FMNN
!LOGICAL :: lxsfperturb
!
!CONTAINS
!
!SUBROUTINE cuInitNNFSP(Core, FmInfo, CmInfo, TranCntl, nTracerCntl, eigv)
!USE PARAM
!USE TYPEDEF,          ONLY: coreinfo_type,    FMInfo_TYPE,    CMInfo_Type,      TranCntl_TYPE,    &
!                            XSNoise_TYPE
!USE CNTL,             ONLY: nTracerCntl_Type
!USE MOC_COMMON
!IMPLICIT NONE
!TYPE(coreinfo_type) :: Core
!TYPE(FMInfo_TYPE) :: FMInfo
!TYPE(CMInfo_Type) :: CMInfo
!TYPE(TranCntl_TYPE) :: TranCntl
!TYPE(nTracerCntl_Type) :: nTracerCntl
!REAL :: eigv
!
!TYPE(XSNoise_TYPE), POINTER :: XsNoise(:)
!INTEGER :: nNoise
!INTEGER :: i
!
!nNoise = TranCntl%nNoise
!XsNoise => TranCntl%XsNoise
!
!lxsfperturb = .FALSE.
!DO i = 1, nNoise
!  IF(XsNoise(i)%itype .EQ. 2) THEN
!    lxsfperturb = .TRUE.
!    EXIT
!  END IF
!END DO
!
!CALL AllocNNVar()
!
!END SUBROUTINE
!
!SUBROUTINE AllocNNVar()
!IMPLICIT NONE
!INTEGER :: nfxr, nfsr, nxy, ng
!INTEGER :: myzb, myze
!
!nfxr = cuGeometry%nfxr
!nfsr = cuGeometry%nfsr
!nxy = cuGeometry%nxy
!ng = cuGeometry%ng
!
!myzb = cuDevice%myzb
!myze = cuDevice%myze
!
!ALLOCATE(FMNN%phis(2*ng, nfsr, myzb:myze))
!ALLOCATE(FMNN%Jout(3, 2*ng, 4, nxy))
!ALLOCATE(FMNN%nnsrc(2*ng, nfsr, myzb:myze))
!ALLOCATE(FMNN%psi(nfsr, myzb:myze))
!ALLOCATE(FMNN%axsrc(2*ng, nxy, myzb:myze), FMNN%AxPxs(2*ng, nxy, myzb:myze))
!
!END SUBROUTINE
!
!SUBROUTINE SetNNSource(Core, FmInfo, nTracerCntl, TranCntl)
!USE TYPEDEF,      ONLY: coreinfo_type, TranCntl_TYPE, XSNoise_TYPE, pin_type, &
!                        cell_type,     FMInfo_TYPE
!USE CNTL,         ONLY : nTracerCntl_Type
!USE ioutil,       ONLY : terminate
!USE PE_Mod,       ONLY : PE
!IMPLICIT NONE
!TYPE(coreinfo_type) :: Core
!TYPE(FMInfo_TYPE) :: FmInfo
!TYPE(nTracerCntl_Type) :: nTracerCntl
!TYPE(TranCntl_TYPE) :: TranCntl
!
!TYPE(pin_type), POINTER :: Pin(:)
!TYPE(cell_type), POINTER :: Cell(:)
!TYPE(XSNoise_TYPE), POINTER :: XsNoise(:)
!INTEGER :: nNoise
!INTEGER :: izbeg, izend, FxrIdxSt, FsrIdxSt
!INTEGER :: iasytype, loc_ixy, ixy, icel, iz
!INTEGER :: i, j
!
!nNoise = TranCntl%nNoise
!XsNoise => TranCntl%XsNoise
!
!Pin => Core%Pin
!Cell => Core%CellInfo
!
!FMNN%nnsrc = 0.
!DO i = 1, nNoise
!  IF(abs(TranCntl%freq_now - XsNoise(i)%freq) .GT. 1.e-5) CYCLE
!  iasytype = Core%Asy(XsNoise(i)%ixya)%AsyType
!  loc_ixy = Core%AsyInfo(iasytype)%pin2DIdx(XsNoise(i)%ix, XsNoise(i)%iy)
!  ixy = Core%Asy(XsNoise(i)%ixya)%GlobalPinIdx(loc_ixy)
!  FxrIdxSt = Pin(ixy)%FxrIdxSt; FsrIdxSt = Pin(ixy)%FsrIdxSt
!  izbeg = XsNoise(i)%izbeg
!  izend = XsNoise(i)%izend
!  DO iz = izbeg, izend
!    IF(iz .LT. PE%myzb) CYCLE
!    IF(iz .GT. PE%myze) CYCLE
!    icel = Pin(ixy)%Cell(iz)
!    IF(.NOT. nTracerCntl%lXslib) THEN
!      IF(TranCntl%lDynamicBen) THEN
!        CALL terminate('NNFSP with Dynamic_bench libary is under development')
!      ELSE
!        CALL SetCellNNSrc_ben(FmInfo, Cell(icel), XsNoise(i), FxrIdxSt, FsrIdxSt, iz, &
!                         nTracerCntl%lTrCorrection, nTracerCntl%lScat1)
!      END IF
!    ELSE
!      CALL terminate('NNFSP with xslib is under development')
!    END IF
!  END DO
!END DO
!
!NULLIFY(XsNoise)
!NULLIFY(Pin)
!NULLIFY(Cell)
!
!END SUBROUTINE
!
!SUBROUTINE SetCellNNSrc_ben(FmInfo, myCell, XsNoise, FxrIdxSt, FsrIdxSt, iz, lTrCorrection, lScat1)
!USE TYPEDEF,    ONLY: FmInfo_Type, cell_type, XSNoise_TYPE, Fxrinfo_type, &
!                      XsMac_Type
!USE BenchXs,    ONLY: xsbaseBen, ChidBen, DnpBetaBen, DnpLambdaBen
!IMPLICIT NONE
!TYPE(FMInfo_TYPE) :: FmInfo
!TYPE(cell_type) :: myCell
!TYPE(XSNoise_TYPE) :: XsNoise
!INTEGER :: FxrIdxSt, FsrIdxSt, iz
!LOGICAL :: lTrCorrection, lScat1
!
!TYPE(Fxrinfo_type), POINTER :: Fxr(:,:)
!TYPE(XsMac_Type) :: xsMac
!REAL, POINTER :: chid(:), beta(:), lambda(:), chieff(:)
!REAL, POINTER :: xstnn(:), xsnfnn(:)
!REAL :: xsa, xsf, xsnf, delxs, betatot, locpsi(2)
!REAL :: amp, freq, phase
!INTEGER :: ng, nprec
!INTEGER :: nLocalFxr, nFsrinFxr
!INTEGER :: ifxr, ifsr, imix, ig, igr, igi, nntype, iprec
!INTEGER :: i, j
!
!ng = cuGeometry%ng
!Fxr => FmInfo%Fxr
!nLocalFxr = myCell%nFXR
!
!amp = XsNoise%amp
!freq = XsNoise%freq
!phase = xsnoise%phase
!nntype = XsNoise%itype
!
!ALLOCATE(xstnn(2*ng))
!IF(nntype .EQ. 2) THEN
!  nprec = cuGeometry%nprec
!  ALLOCATE(chid(ng), beta(nprec), lambda(nprec))
!  ALLOCATE(xsnfnn(2*ng), chieff(2*ng))
!  CALL ChidBen(chid)
!END IF
!
!DO j = 1, nLocalFxr
!  ifxr = FxrIdxSt + j - 1
!  nFsrinFxr = myCell%nFsrInFxr(j)
!  imix = Fxr(ifxr,iz)%imix
!  CALL xsbaseBen(imix, 1, ng, 1, ng, lScat1, xsMac)
!  IF(nntype .EQ. 2) THEN
!    CALL DnpBetaBen(imix, beta)
!    CALL DnpLambdaBen(lambda)
!    betatot = sum(beta)
!  END IF
!  DO ig = 1, ng
!    igr = 2*ig-1
!    igi = 2*ig
!    SELECT CASE(nntype)
!    CASE(1) ! Capture XS
!      delxs = amp * (xsMac%xsmaca(ig) - xsMac%xsmackf(ig))
!      xstnn(igr) = 0.5 * delxs * sin(phase)
!      xstnn(igi) = 0.5 * delxs * cos(phase)
!    CASE(2) ! Fission XS
!      delxs = amp * xsMac%xsmackf(ig)
!      xstnn(igr) = 0.5 * delxs * sin(phase)
!      xstnn(igi) = 0.5 * delxs * cos(phase)
!      delxs = amp * xsMac%xsmacnf(ig)
!      xsnfnn(igr) = 0.5 * delxs * sin(phase)
!      xsnfnn(igi) = 0.5 * delxs * cos(phase)
!
!      !chieff(igr) = (1-betatot) * chi(ig)
!      chieff(igi) = 0.
!      DO iprec = 1, nprec
!        chieff(igr) = chieff(igr) + chid(ig) * (lambda(iprec) ** 2) * beta(iprec) / ((freq**2)+(lambda(iprec)**2))
!        chieff(igi) = chieff(igi) - chid(ig) * freq * lambda(iprec) * beta(iprec) / ((freq**2)+(lambda(iprec)**2))
!      END DO
!    CASE(3) ! Scattering XS
!    END SELECT
!  END DO
!
!  DO i = 1, nFsrinFxr
!    ifsr = FsrIdxSt + myCell%MapFxr2FsrIdx(i,j) - 1
!    IF(nntype .EQ. 2) THEN
!      locpsi = 0.
!      DO ig = 1, ng
!        igr = 2*ig-1
!        igi = 2*ig
!        locpsi(1) = locpsi(1) + xsnf(igr) * FmInfo%phis(ifsr, iz, ig)
!        locpsi(2) = locpsi(2) + xsnf(igi) * FmInfo%phis(ifsr, iz, ig)
!      END DO
!    END IF
!    DO ig = 1, ng
!      igr = 2*ig-1
!      igi = 2*ig
!      FMNN%nnsrc(igr, ifsr, iz) = FMNN%nnsrc(igr, ifsr, iz) - xstnn(igr) * FmInfo%phis(ifsr, iz, ig)
!      FMNN%nnsrc(igi, ifsr, iz) = FMNN%nnsrc(igi, ifsr, iz) - xstnn(igi) * FmInfo%phis(ifsr, iz, ig)
!      IF(nntype .EQ. 2) THEN
!        FMNN%nnsrc(igr, ifsr, iz) = FMNN%nnsrc(igr, ifsr, iz) + chieff(igr) * locpsi(1) - chieff(igi) * locpsi(2)
!        FMNN%nnsrc(igi, ifsr, iz) = FMNN%nnsrc(igi, ifsr, iz) + chieff(igr) * locpsi(2) + chieff(igi) * locpsi(1)
!      END IF
!    END DO
!  END DO
!END DO
!
!DEALLOCATE(xstnn)
!IF(nntype .EQ. 2) THEN
!   DEALLOCATE(chid, lambda, beta)
!   DEALLOCATE(xsnfnn, chieff)
!END IF
!NULLIFY(Fxr)
!
!END SUBROUTINE
!
!END MODULE
!#endif
