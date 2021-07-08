SUBROUTINE MOCSolUpdt(Core, FmInfo, CmInfo, myzb, myze, ng)
USE PARAM
USE TYPEDEF,  ONLY : CoreInfo_Type,     FmInfo_Type,      CmInfo_Type,   &
                     FXRInfo_Type,      PinXs_Type,                      &
                     Pin_Type,          PinInfo_Type,      Cell_Type
USE PE_Mod, ONLY : PE
USE CNTL,   ONLY : nTracerCntl
USE ITRCNTL_MOD, ONLY : ItrCntl
IMPLICIT NONE

TYPE(CoreInfo_Type), INTENT(IN) :: Core
TYPE(FmInfo_Type), INTENT(INOUT) :: FmInfo
TYPE(CmInfo_Type), INTENT(INOUT) :: CmInfo
INTEGER, INTENT(IN) :: myzb, myze, ng

TYPE(PinXs_Type), POINTER :: PinXs(:, :)
REAL, POINTER :: Phis(:, :, :), Phim(:, :, :, :), PHIC(:, : , :), Jout(:, :, :, :, :)

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(PinInfo_Type), POINTER :: PinInfo(:)
TYPE(Cell_Type), POINTER :: Cell(:)

INTEGER :: nCoreFxr, nCoreFsr, nxy
INTEGER :: nlocalFsr
INTEGER :: icel, ipin, iz, ixy ,ifxr, ifsr, itype, ifsrlocal
INTEGER :: FsrIdxSt, FxrIdxSt
INTEGER :: i, j, k, l, m, ig, ibd, ineigh, ierr
REAL :: fmult, fmultMg(ng), ur0, ur, ur_t

REAL :: dhat, dtil, jnet, jnet0
REAL :: myphi, neighphi

! INCLUDE 'mpif.h'

Pin => Core%Pin
PinInfo => Core%Pininfo;   Cell => Core%CellInfo
Phis => FmInfo%Phis; Phim => FmInfo%Phim
PinXS => CmInfo%PinXs; PhiC => CmInfo%PhiC
Jout => FmInfo%RadJout
nxy = Core%nxy
nCoreFxr = Core%nCoreFxr; nCoreFsr = Core%nCoreFsr

! ur=1._8
! DO ixy = 1, nxy
!   DO iz = myzb, myze
!     DO ig = 1, ng
!       FsrIdxSt = Pin(ixy)%FsrIdxSt;
!       icel = Pin(ixy)%Cell(iz)
!       nlocalFsr = Cell(icel)%nFsr
!       IF(PhiC(ixy, iz, ig).LT.0)THEN
!           ur0=PinXS(ixy, iz)%Phi(ig)/(PinXS(ixy, iz)%Phi(ig)-PhiC(ixy, iz, ig))
!           IF( ur0 .LT. ur )THEN
!               ur=ur0
!           ENDIF
!           IF( ur.LT. 0.1 )THEN
!               ur=0.1
!               PhiC(ixy, iz, ig)=0;
!           ENDIF
!       ENDIF
!     ENDDO
!   ENDDO
! ENDDO
! 
! CALL MPI_ALLREDUCE(ur, ur_t, 1, MPI_DOUBLE_PRECISION, MPI_MIN, PE%MPI_CMFD_COMM, ierr)
! ur = ur_t
! ur=1.0
! IF (ur .NE. 1 .AND. PE%MASTER) WRITE(*,'(a,f12.5)') '-----UnderRelaxation Factor CMFD', ur

DO ig = 1, ng
  DO iz = myzb, myze
    DO ixy = 1, nxy
      FsrIdxSt = Pin(ixy)%FsrIdxSt;
      icel = Pin(ixy)%Cell(iz)
      nlocalFsr = Cell(icel)%nFsr
!      PhiC(ixy, iz, ig)=ur*PhiC(ixy, iz, ig)+PinXS(ixy, iz)%Phi(ig)*(1._8-ur) ! under relaxation
      !CmInfo%PhiFM(ixy, iz, ig)=PhiC(ixy, iz, ig)
      fmult = PhiC(ixy, iz, ig)/PinXS(ixy, iz)%Phi(ig)
      DO j = 1, nLocalFsr
        ifsr = FsrIdxSt + j - 1
        phis(ifsr, iz, ig) = phis(ifsr, iz, ig)*fmult
      ENDDO
    ENDDO
  ENDDO
ENDDO

IF (nTracerCntl%lScat1) THEN
  DO ig = 1, ng
    DO iz = myzb, myze
      DO ixy = 1, nxy
        FsrIdxSt = Pin(ixy)%FsrIdxSt
        icel = Pin(ixy)%Cell(iz)
        nLocalFsr = Cell(icel)%nFsr
        fmult = Phic(ixy, iz, ig) / PinXS(ixy, iz)%Phi(ig)
        IF (nTracerCntl%lNodeMajor) THEN
          DO j = 1, nLocalFsr
            ifsr = FsrIdxSt + j - 1
            phim(:, ig, ifsr, iz) = phim(:, ig, ifsr, iz) * fmult
          ENDDO
        ENDIF
      ENDDO
    ENDDO
  ENDDO
ENDIF

NULLIFY(Pin, PinInfo ,Cell)

#ifdef COMPILE

IF (.NOT. nTracerCntl%lscat1) RETURN

IF (nTracerCntl%lNodeMajor) THEN
  DO iz = myzb, myze
    DO ixy = 1, nxy
      FsrIdxSt = Pin(ixy)%FsrIdxSt; 
      icel = Pin(ixy)%Cell(iz)
      nlocalFsr = Cell(icel)%nFsr
      fmultMg = PhiC(ixy, iz, :)/PinXS(ixy, iz)%Phi(:)
      DO j = 1, nLocalFsr
        ifsr = FsrIdxSt + j - 1
        DO ig = 1, ng
          phim(:, ig, ifsr, iz) = phim(:, ig, ifsr, iz) * fmultMg(ig)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSE
  DO ig = 1, ng
    DO iz = myzb, myze
      DO ixy = 1, nxy
        FsrIdxSt = Pin(ixy)%FsrIdxSt; 
        icel = Pin(ixy)%Cell(iz)
        nlocalFsr = Cell(icel)%nFsr
        fmult = PhiC(ixy, iz, ig)/PinXS(ixy, iz)%Phi(ig)
        DO j = 1, nLocalFsr
          ifsr = FsrIdxSt + j - 1
          phim(:, ifsr, iz, ig) = phim(:, ifsr, iz, ig) * fmult
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF

#endif

END SUBROUTINE

SUBROUTINE MOCSolUpdt_old(Core, FmInfo, CmInfo, myzb, myze, ng)
USE PARAM
USE TYPEDEF,  ONLY : CoreInfo_Type,     FmInfo_Type,      CmInfo_Type,   &
                     FXRInfo_Type,      PinXs_Type,                      &
                     Pin_Type,          PinInfo_Type,      Cell_Type
IMPLICIT NONE

TYPE(CoreInfo_Type), INTENT(IN) :: Core
TYPE(FmInfo_Type), INTENT(INOUT) :: FmInfo
TYPE(CmInfo_Type), INTENT(INOUT) :: CmInfo
INTEGER, INTENT(IN) :: myzb, myze, ng

TYPE(PinXs_Type), POINTER :: PinXs(:, :)
REAL, POINTER :: Phis(:, :, :), PHIC(:, : , :), Jout(:, :, :, :, :)

TYPE(Pin_Type), POINTER :: Pin(:)
!TYPE(AsyInfo_Type), POINTER :: AsyInfo
TYPE(PinInfo_Type), POINTER :: PinInfo(:)
TYPE(Cell_Type), POINTER :: Cell(:)

INTEGER :: nCoreFxr, nCoreFsr, nxy
INTEGER :: nlocalFsr
INTEGER :: icel, ipin, iz, ixy ,ifxr, ifsr, itype, ifsrlocal
INTEGER :: FsrIdxSt, FxrIdxSt
INTEGER :: i, j, k, l, m, ig, ibd, ineigh
REAL :: fmult

REAL :: dhat, dtil, jnet, jnet0
REAL :: myphi, neighphi

Pin => Core%Pin
PinInfo => Core%Pininfo;   Cell => Core%CellInfo
Phis => FmInfo%Phis;
PinXS => CmInfo%PinXs; PhiC => CmInfo%PhiC
Jout => FmInfo%RadJout
nxy = Core%nxy
nCoreFxr = Core%nCoreFxr; nCoreFsr = Core%nCoreFsr
DO ig = 1, ng
  DO iz = myzb, myze
    DO ixy = 1, nxy
      FsrIdxSt = Pin(ixy)%FsrIdxSt;
      icel = Pin(ixy)%Cell(iz)
      nlocalFsr = Cell(icel)%nFsr
      fmult = PhiC(ixy, iz, ig)/PinXS(ixy, iz)%Phi(ig)
      DO j = 1, nLocalFsr
        ifsr = FsrIdxSt + j - 1
        phis(ifsr, iz, ig) = phis(ifsr, iz, ig)*fmult
      ENDDO
    ENDDO
  ENDDO
ENDDO

NULLIFY(Pin, PinInfo ,Cell)

END SUBROUTINE

SUBROUTINE MOCSolUpdtUr(Core, FmInfo, CmInfo, myzb, myze, ng)
USE PARAM
USE TYPEDEF,  ONLY : CoreInfo_Type,     FmInfo_Type,      CmInfo_Type,   &
                     FXRInfo_Type,      PinXs_Type,                      &
                     Pin_Type,          PinInfo_Type,      Cell_Type
IMPLICIT NONE

TYPE(CoreInfo_Type), INTENT(IN) :: Core
TYPE(FmInfo_Type), INTENT(INOUT) :: FmInfo
TYPE(CmInfo_Type), INTENT(INOUT) :: CmInfo
INTEGER, INTENT(IN) :: myzb, myze, ng

TYPE(PinXs_Type), POINTER :: PinXs(:, :)
REAL, POINTER :: Phis(:, :, :), PHIC(:, : , :), Jout(:, :, :, :, :)

TYPE(Pin_Type), POINTER :: Pin(:)
!TYPE(AsyInfo_Type), POINTER :: AsyInfo
TYPE(PinInfo_Type), POINTER :: PinInfo(:)
TYPE(Cell_Type), POINTER :: Cell(:)

INTEGER :: nCoreFxr, nCoreFsr, nxy
INTEGER :: nlocalFsr
INTEGER :: icel, ipin, iz, ixy ,ifxr, ifsr, itype, ifsrlocal
INTEGER :: FsrIdxSt, FxrIdxSt
INTEGER :: i, j, k, l, m, ig, ibd, ineigh
REAL :: fmult

REAL :: dhat, dtil, jnet, jnet0
REAL :: myphi, neighphi

Pin => Core%Pin
PinInfo => Core%Pininfo;   Cell => Core%CellInfo
Phis => FmInfo%Phis;
PinXS => CmInfo%PinXs; PhiC => CmInfo%PhiC
Jout => FmInfo%RadJout
nxy = Core%nxy
nCoreFxr = Core%nCoreFxr; nCoreFsr = Core%nCoreFsr
DO ig = 1, ng
  DO iz = myzb, myze
    DO ixy = 1, nxy
      FsrIdxSt = Pin(ixy)%FsrIdxSt;
      icel = Pin(ixy)%Cell(iz)
      nlocalFsr = Cell(icel)%nFsr
      fmult = PhiC(ixy, iz, ig)/PinXS(ixy, iz)%Phi(ig)
      DO j = 1, nLocalFsr
        ifsr = FsrIdxSt + j - 1
        phis(ifsr, iz, ig) = phis(ifsr, iz, ig)*fmult
      ENDDO
    ENDDO
  ENDDO
ENDDO

DO ig = 1, ng
  DO iz = myzb, myze
     DO ixy = 1, nxy
       myphi = PhiC(ixy, iz, ig)
       DO ibd = 1, 4
         jnet0 = (Jout(2, ibd, ixy, iz, ig) - Jout(1, ibd, ixy, iz, ig))
         dtil = PinXs(ixy, iz)%dtil(ibd, ig)
         dhat = PinXs(ixy, iz)%dhat(ibd, ig)
         ineigh = Pin(ixy)%NeighIdx(ibd)
         neighphi = 0
         IF(ineigh .GT. 0) neighphi = Phic(ineigh, iz, ig)
         jnet = (dtil-dhat) * myphi -(dtil+dhat) * neighphi
         fmult = jnet / jnet0
         Jout(2, ibd, ixy, iz, ig) = Jnet
         Jout(1, ibd, ixy, iz, ig) = 0
       ENDDO
     ENDDO
  ENDDO
ENDDO

NULLIFY(Pin, PinInfo ,Cell)

END SUBROUTINE


!SUBROUTINE MOCSolUpdt(Core, PinXS, PHIC, phis, myzb, myze, ng)
!USE PARAM
!USE TYPEDEF,  ONLY : CoreInfo_Type,    FXRInfo_Type,      PinXs_Type,        &
!                     Pin_Type,         PinInfo_Type,      Cell_Type
!IMPLICIT NONE
!
!TYPE(CoreInfo_Type), INTENT(IN) :: Core
!TYPE(PinXs_Type), POINTER, INTENT(INOUT) :: PinXs(:, :)
!REAL, POINTER :: Phis(:, :, :), PHIC(:, : , :)
!INTEGER, INTENT(IN) :: myzb, myze, ng
!
!TYPE(Pin_Type), POINTER :: Pin(:)
!!TYPE(AsyInfo_Type), POINTER :: AsyInfo
!TYPE(PinInfo_Type), POINTER :: PinInfo(:)
!TYPE(Cell_Type), POINTER :: Cell(:)
!
!INTEGER :: nCoreFxr, nCoreFsr, nxy
!INTEGER :: nlocalFsr
!INTEGER :: icel, ipin, iz, ixy ,ifxr, ifsr, itype, ifsrlocal
!INTEGER :: FsrIdxSt, FxrIdxSt
!INTEGER :: i, j, k, l, m, ig
!REAL :: fmult
!Pin => Core%Pin
!PinInfo => Core%Pininfo;   Cell => Core%CellInfo
!
!nxy = Core%nxy
!nCoreFxr = Core%nCoreFxr; nCoreFsr = Core%nCoreFsr
!DO ig = 1, ng
!  DO iz = myzb, myze
!    DO ixy = 1, nxy
!      FsrIdxSt = Pin(ixy)%FsrIdxSt;
!      icel = Pin(ixy)%Cell(iz)
!      nlocalFsr = Cell(icel)%nFsr
!      fmult = PhiC(ixy, iz, ig)/PinXS(ixy, iz)%Phi(ig)
!      DO j = 1, nLocalFsr
!        ifsr = FsrIdxSt + j - 1
!        phis(ifsr, iz, ig) = phis(ifsr, iz, ig)*fmult
!      ENDDO
!    ENDDO
!  ENDDO
!ENDDO
!NULLIFY(Pin, PinInfo ,Cell)
!
!END SUBROUTINE

SUBROUTINE MOCLinSrcUpdt(Core, PinXS, PHIC, LinSrcSlope, myzb, myze, ng)
USE PARAM
USE TYPEDEF,  ONLY : CoreInfo_Type,    FXRInfo_Type,      PinXs_Type,        &
                     Pin_Type,         PinInfo_Type,      Cell_Type
IMPLICIT NONE

TYPE(CoreInfo_Type), INTENT(IN) :: Core
TYPE(PinXs_Type), POINTER, INTENT(INOUT) :: PinXs(:, :)
REAL, POINTER :: LinSrcSlope(:, :, :, :), PHIC(:, : , :)
INTEGER, INTENT(IN) :: myzb, myze, ng

TYPE(Pin_Type), POINTER :: Pin(:)
!TYPE(AsyInfo_Type), POINTER :: AsyInfo
TYPE(PinInfo_Type), POINTER :: PinInfo(:)
TYPE(Cell_Type), POINTER :: Cell(:)

INTEGER :: nCoreFxr, nCoreFsr, nxy
INTEGER :: nlocalFsr
INTEGER :: icel, ipin, iz, ixy ,ifxr, ifsr, itype, ifsrlocal
INTEGER :: FsrIdxSt, FxrIdxSt
INTEGER :: i, j, k, l, m, ig
REAL :: fmult
Pin => Core%Pin
PinInfo => Core%Pininfo;   Cell => Core%CellInfo

nxy = Core%nxy
nCoreFxr = Core%nCoreFxr; nCoreFsr = Core%nCoreFsr
DO ig = 1, ng
  DO iz = myzb, myze
    DO ixy = 1, nxy
      FsrIdxSt = Pin(ixy)%FsrIdxSt;
      icel = Pin(ixy)%Cell(iz)
      nlocalFsr = Cell(icel)%nFsr
      fmult = PhiC(ixy, iz, ig)/PinXS(ixy, iz)%Phi(ig)
      DO j = 1, nLocalFsr
        ifsr = FsrIdxSt + j - 1
        LinSrcSlope(1:2, ifsr, iz, ig) = LinSrcSlope(1:2, ifsr, iz, ig)*fmult
      ENDDO
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE MOCPhiInUpdt(Core, CmInfo, myzb, myze, ng)

USE PARAM,       ONLY : ZERO
USE TYPEDEF,     ONLY : CoreInfo_Type, Pin_Type, PinXs_Type, RayInfo4CMFD_Type, CmInfo_Type
USE CMFD_MOD,    ONLY : PinNeighIdx, CMFDPinXS, SubPlaneMap
USE CNTL,        ONLY : nTracerCntl
USE itrcntl_mod, ONLY : ItrCntl

IMPLICIT NONE

TYPE (CoreInfo_Type) :: Core
TYPE (CmInfo_Type)   :: CMInfo

INTEGER :: myzb, myze, ng

TYPE (RayInfo4Cmfd_Type), POINTER :: RayInfo4CMFD

TYPE (Pin_Type),   POINTER, DIMENSION(:)   :: Pin
TYPE (PinXs_Type), POINTER, DIMENSION(:,:) :: PinXs

REAL, POINTER, DIMENSION(:,:,:) :: PHIC, PhiFm

INTEGER, POINTER, DIMENSION(:)     :: DcmpAsyRayCount
INTEGER, POINTER, DIMENSION(:,:)   :: RotRayInOutCell, PhiangInSvIdx
INTEGER, POINTER, DIMENSION(:,:,:) :: DcmpAsyRayInSurf, DcmpAsyRayInCell

INTEGER :: iz, iz0, ig, ipol, ixy, iAsy, isv, isurf, ingh, nRotRay, nCoreRay, nPolar, nAsy, icnt, idir, iRotRay

REAL, POINTER, DIMENSION(:,:,:,:)     :: PhiAngIn
REAL, POINTER, DIMENSION(:,:,:,:,:)   :: RadJout
REAL, POINTER, DIMENSION(:,:,:,:,:,:) :: AsyPhiAngIn

REAL :: myphi, neighphi, surfphi, atil, ahat, smy, fmult
! ----------------------------------------------------

RayInfo4Cmfd => CmInfo%RayInfo4Cmfd
PhiC         => CmInfo%PhiC
PhiFm        => CmInfo%PhiFm
PinXS        => CmInfo%PinXS
RadJout      => CmInfo%RadJout

Pin => Core%Pin
nAsy = Core%nxya

nRotRay           = RayInfo4CMFD%nRotRay
nCoreRay          = RayInfo4CMFD%nCoreRay
nPolar            = RayInfo4CMFD%nPolAngle
RotRayInOutCell  => RayInfo4CMFD%RotRayInOutCell
PhiAngInSvIdx    => RayInfo4CMFD%PhiAngInSvIdx
PhiAngIn         => RayInfo4CMFD%PhiAngIn
DcmpAsyRayCount  => RayInfo4CMFD%DcmpAsyRayCount
DcmpAsyRayInSurf => RayInfo4CMFD%DcmpAsyRayInSurf
DcmpAsyRayInCell => RayInfo4CMFD%DcmpAsyRayInCell
AsyPhiAngIn      => RayInfo4CMFD%AsyPhiAngIn
! ----------------------------------------------------
DO idir = 1, 2
  DO iRotRay = 1, nRotRay
    ixy = RotRayInOutCell(iRotRay, idir)
    
    IF (ixy .EQ. 0) CYCLE
    
    isv = PhiAngInSvIdx(iRotRay, idir)
    
    DO iz = myzb, myze
      DO ipol = 1, nPolar
        DO ig = 1, ng
          fmult = PhiC(ixy, iz, ig) / PinXS(ixy, iz)%phi(ig)
          
          PhiAngIn(ipol, isv, iz, ig) = PhiAngIn(ipol, isv, iz, ig) * fmult
        END DO
      END DO
    END DO
  END DO
END DO
! ----------------------------------------------------
IF (nTracerCntl%lDomainDcmp) THEN
  DO iz = myzb, myze
    iz0 = SubPlaneMap(iz)
    
    DO iAsy = 1, nAsy
      DO icnt = 1, DcmpAsyRayCount(iAsy)
        DO idir = 1, 2
          ixy   = DcmpAsyRayInCell(idir, icnt, iAsy)
          isurf = DcmpAsyRayInSurf(idir, icnt, iAsy)
          ingh  = PinNeighIdx(isurf, ixy)
          smy   = Pin(ixy)%BdLength(isurf)
          
          DO ig = 1, ng
            myphi = PhiFm(ixy, iz, ig)
            
            neighphi = ZERO
            IF (ingh .GT. 0) neighphi = PhiFm(ingh, iz, ig)
            
            atil = CMFDPinXS(ixy, iz0)%atil(isurf, ig)
            ahat = CMFDPinXS(ixy, iz0)%ahat(isurf, ig)
            
            surfphi = atil * myphi + (smy - atil) * neighphi
            
            IF (ItrCntl%mocit .EQ. 0) THEN
              IF (nTracerCntl%lNodeMajor) THEN
                AsyPhiAngIn(:, ig, idir, icnt, iAsy, iz) = surfphi / smy
              ELSE
                AsyPhiAngIn(:, idir, icnt, iAsy, ig, iz) = surfphi / smy
              END IF
            ELSE
              surfphi = surfphi + ahat * (myphi + neighphi)
              fmult   = surfphi / RadJout(3, isurf, ixy, iz, ig)
              
              IF (nTracerCntl%lNodeMajor) THEN
                AsyPhiAngIn(:, ig, idir, icnt, iAsy, iz) = AsyPhiAngIn(:, ig, idir, icnt, iAsy, iz) * fmult
              ELSE
                AsyPhiAngIn(:, idir, icnt, iAsy, ig, iz) = AsyPhiAngIn(:, idir, icnt, iAsy, ig, iz) * fmult
              END IF
            END IF
          END DO
        END DO
      END DO
    END DO
  END DO
END IF
! ----------------------------------------------------
! Dcmp.
NULLIFY (RotRayInOutCell)
NULLIFY (PhiAngInSvIdx)
NULLIFY (PhiAngIn)
NULLIFY (DcmpAsyRayCount)
NULLIFY (DcmpAsyRayInSurf)
NULLIFY (DcmpAsyRayInCell)
NULLIFY (AsyPhiAngIn)

! Geo
NULLIFY (Pin)
NULLIFY (RayInfo4Cmfd)
NULLIFY (PhiC)
NULLIFY (PhiFm)
NULLIFY (PinXS)
NULLIFY (RadJout)
! ----------------------------------------------------

END SUBROUTINE MOCPhiInUpdt
! ------------------------------------------------------------------------------------------------------------