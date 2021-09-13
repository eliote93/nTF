#include <defines.h>
#ifdef __INTEL_MKL
MODULE MKL_HOMOXS

USE MKL_3D

IMPLICIT NONE

CONTAINS
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HomogenizePnXS(CoreInfo, FmInfo, GroupInfo)

USE PARAM,        ONLY : TRUE, FALSE
USE TYPEDEF,      ONLY : CoreInfo_Type, FmInfo_Type, GroupInfo_Type, FxrInfo_Type, Pin_Type, XsMac_Type
USE MacXsLib_Mod, ONLY : MacP1XsScatMatrix, MacP2XsScatMatrix, MacP3XsScatMatrix

IMPLICIT NONE

TYPE (CoreInfo_Type)  :: CoreInfo
TYPE (FmInfo_Type)    :: FmInfo
TYPE (GroupInfo_Type) :: GroupInfo
! ----------------------------------------------------
TYPE (superPin_Type), POINTER, DIMENSION(:)   :: superPin
TYPE (Pin_Type),      POINTER, DIMENSION(:)   :: Pin
TYPE (FxrInfo_Type),  POINTER, DIMENSION(:,:) :: Fxr

TYPE (XsMac_Type) :: XsMac

INTEGER :: ig, ifxr, ixy, ixy_map, ipin, jz, iz, ng, nxy, nzCMFD
INTEGER, POINTER, DIMENSION(:) :: pinMap, planeMap

LOGICAL :: lFreeSth
! ----------------------------------------------------

Pin => CoreInfo%Pin
Fxr => FmInfo%Fxr

ng        = mklGeom%ng
nxy       = mklGeom%nxy
nzCMFD    = mklGeom%nzCMFD
superPin => mklGeom%superPin
planeMap => mklGeom%planeMap
pinMap   => mklGeom%pinMap

lFreeSth = FALSE

DO iz = 1, nzCMFD ! Coarse Pln.
  jz = planeMap(iz)
  
  DO ixy = 1, nxy
    IF (.NOT. mklGeom%lH2OCell(iz, ixy)) CYCLE
    
    lFreeSth = TRUE
    ixy_map = pinMap(ixy)
    ipin = superPin(ixy_map)%pin(1)
    ifxr = Pin(ipin)%FxrIdxSt
    
    CALL MacP1XsScatMatrix(XsMac, Fxr(ifxr, jz), 1, ng, ng, GroupInfo)
!    CALL MacP2XsScatMatrix(XsMac, Fxr(ifxr, jz), 1, ng, ng, GroupInfo)
!    CALL MacP3XsScatMatrix(XsMac, Fxr(ifxr, jz), 1, ng, ng, GroupInfo)
    
    mklAxial%SmP1(:, :, iz, ixy) = XsMac%XsMacP1Sm
!    mklAxial%SmP2(:, :, iz, ixy) = XsMac%XsMacP2Sm
!    mklAxial%SmP3(:, :, iz, ixy) = XsMac%XsMacP3Sm
  END DO
END DO

IF (lFreeSth) THEN
  DEALLOCATE(XsMac%xsmaca, XsMac%xsmacf, XsMac%xsmackf, XsMac%xsmacnf, XsMac%xsmact,XsMac%xsmactr)
  DEALLOCATE(XsMac%xsmacp1sm, XsMac%xsmacp2sm, XsMac%xsmacp3sm)
  DEALLOCATE(XsMac%xsmacs, XsMac%xsmacstr, XsMac%xsmacsm)
END IF
! ----------------------------------------------------

END SUBROUTINE HomogenizePnXS
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HomogenizeGcXS(CoreInfo, PinXS, GcPinXS)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,       PinXS_Type
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(PinXS_Type), POINTER :: PinXS(:, :), GcPinXS(:, :)

REAL :: localphis(mklGeom%ng)
INTEGER :: nxy, nzCMFD
INTEGER :: ipin, ipin_map, iz, izf
INTEGER, POINTER :: pinMap(:), planeMap(:)

nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
pinMap => mklGeom%pinMap
planeMap => mklGeom%planeMap

DO izf = 1, nzCMFD
  iz = planeMap(izf)
  !$OMP PARALLEL PRIVATE(localphis, ipin_map)
  !$OMP DO SCHEDULE(GUIDED)
  DO ipin = 1, nxy
    ipin_map = pinMap(ipin)
    localphis = mklCMFD%phis(ipin, izf, :)
    CALL HomogenizeCellGcXS(PinXS(ipin_map, iz), GcPinXS(ipin_map, izf), localphis)
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
ENDDO

END SUBROUTINE
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HomogenizeCellGcXS(PinXS, GcPinXS, phis)
USE PARAM
USE TYPEDEF,        ONLY : PinXS_Type
IMPLICIT NONE

TYPE(PinXS_Type) :: PinXS, GcPinXS
REAL :: phis(:)

REAL :: RR(3), RRS(mklGeom%ngc, mklGeom%ngc)
REAL :: phisum, chisum
INTEGER :: ng, ngc
INTEGER :: igc, ig, igb, ige
INTEGER :: igs, igsb, igse, ig0

ng = mklGeom%ng
ngc = mklGeom%ngc

DO igc = 1, ngc
  igb = mklGeom%GcStruct(1, igc); ige = mklGeom%GcStruct(2, igc)
  RR = 0; phisum = 0; chisum = 0
  DO ig = igb, ige
    RR(1) = RR(1) + PinXS%XStr(ig) * phis(ig)
    RR(2) = RR(2) + PinXS%XSnf(ig) * phis(ig)
    RR(3) = RR(3) + PinXS%XSD(ig) * phis(ig)
    phisum = phisum + phis(ig)
    chisum = chisum + PinXS%Chi(ig)
  ENDDO
  RR = RR / phisum
  GcPinXS%XStr(igc) = RR(1)
  GcPinXS%XSnf(igc) = RR(2)
  GcPinXS%XSD(igc) = RR(3)
  GcPinXS%Chi(igc) = chisum
  GcPinXS%Phi(igc) = phisum
ENDDO

RRS = 0
DO ig = 1, ng
  igc = mklGeom%GcStructInv(ig)
  igsb = PinXS%XSs(ig)%ib; igse = PinXS%XSs(ig)%ie
  ig0 = igc
  RRS(ig0, igc) = RRS(ig0, igc) + PinXS%XSs(ig)%WithInGroupScat * phis(ig)
  DO igs = igsb, igse
    ig0 = mklGeom%GcStructInv(igs)
    RRS(ig0, igc) = RRS(ig0, igc) + PinXS%XSs(ig)%from(igs) * phis(igs)
  ENDDO
ENDDO

DO igc = 1, ngc
  RRS(igc, :) = RRS(igc, :) / GcPinXS%Phi(igc)
ENDDO

DO ig = 1, ngc
  igsb = GcPinXS%XSs(ig)%ib; igse = GcPinXS%XSs(ig)%ie
  DO igs = igsb, igse
    GcPinXS%XSs(ig)%from(igs) = RRS(igs, ig)
  ENDDO
  GcPinXS%XSs(ig)%WithInGroupScat = GcPinXS%XSs(ig)%from(ig)
  GcPinXS%XSs(ig)%from(ig) = 0.0
ENDDO

DO ig = 1, ngc
  GcPinXS%XSr(ig) = GcPinXS%XStr(ig) - GcPinXS%XSs(ig)%WithInGroupScat
ENDDO

END SUBROUTINE
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE SetRadialGcCoupling(PinXS, GcPinXS)
USE PARAM
USE TYPEDEF,        ONLY : PinXS_Type
USE CNTL,           ONLY : nTracerCntl
IMPLICIT NONE

TYPE(PinXS_Type), POINTER :: PinXS(:, :), GcPinXS(:, :)

TYPE(superPin_Type), POINTER :: Pin(:)
INTEGER, POINTER :: pinMap(:), pinMapRev(:), planeMap(:)
INTEGER :: ng, ngc, nxy, nzCMFD
INTEGER :: ig, igc, ipin, ipin_map, ineighpin, iz, izf, ibd, inbd
REAL :: Dtil, Dhat, pDhat(2), atil, myphi, neighphi, mybeta, neighbeta, albedo, jfdm, surfphifdm, smy
REAL, POINTER :: Jnet(:, :, :, :), Jpart(:, :, :, :, :)

IF (nTracerCntl%lHex) THEN
  CALL HexSetRadialGcCoupling(PinXs, GcPinXs)
  
  RETURN
END IF

Pin => mklGeom%superPin
ng = mklGeom%ng
ngc = mklGeom%ngc
nxy = mklGeom%nxy
nzCMFD = mklGeom%nzCMFD
pinMap => mklGeom%pinMap
pinMapRev => mklGeom%pinMapRev
planeMap => mklGeom%planeMap

ALLOCATE(Jnet(4, nxy, nzCMFD, ngc)); Jnet = 0.0

!--- Condense Currents

DO igc = 1, ngc
  DO ig = mklGeom%GcStruct(1, igc), mklGeom%GcStruct(2, igc)
    DO izf = 1, nzCMFD
      iz = planeMap(izf)
      !$OMP PARALLEL PRIVATE(ipin_map, ineighpin, myphi, neighphi, Dtil, Dhat, pDhat, jfdm)
      !$OMP DO SCHEDULE(GUIDED)
      DO ipin = 1, nxy
        ipin_map = pinMap(ipin)
        myphi = mklCMFD%phis(ipin, izf, ig)
        DO ibd = 1, 4
          ineighpin = Pin(ipin_map)%Neighidx(ibd)
          ineighpin = pinMapRev(ineighpin)
          IF (ineighpin .LE. 0) THEN
            neighphi = 0.0
          ELSE
            neighphi = mklCMFD%phis(ineighpin, izf, ig)
          ENDIF
          Dtil = PinXS(ipin_map, iz)%Dtil(ibd, ig)
          Dhat = PinXS(ipin_map, iz)%Dhat(ibd, ig)
          jfdm = - Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
          Jnet(ibd, ipin, izf, igc) = Jnet(ibd, ipin, izf, igc) + jfdm
        ENDDO
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL
    ENDDO
  ENDDO
ENDDO

!--- Compute Coupling Coefficients

!$OMP PARALLEL PRIVATE(ipin_map, ineighpin, inbd, myphi, neighphi, mybeta, neighbeta, Dtil, Dhat, pDhat, jfdm, albedo, smy)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO igc = 1, ngc
  DO izf = 1, nzCMFD
    DO ipin = 1, nxy
      ipin_map = pinMap(ipin)
      DO ibd = 1, 4
        ineighpin = Pin(ipin_map)%NeighIdx(ibd)
        smy = Pin(ipin_map)%BdLength(ibd)
        myphi = GcPinXS(ipin_map, izf)%Phi(igc)
        mybeta = GcPinXS(ipin_map, izf)%XSD(igc) / Pin(ipin_map)%Center2SurfaceL(ibd)
        IF (ineighpin .GT. 0) THEN
          inbd = Pin(ineighpin)%NeighSurfIdx(ibd)
          neighphi = GcPinXS(ineighpin, izf)%Phi(igc)
          neighbeta = GcPinXS(ineighpin, izf)%XSD(igc) / Pin(ineighpin)%Center2SurfaceL(inbd)
          Dtil = mybeta * neighbeta / (mybeta + neighbeta) * smy
          jfdm = - Dtil * (neighphi - myphi)
          Dhat = - (Jnet(ibd, ipin, izf, igc) - jfdm) / (myphi + neighphi)
        ELSE
          IF (ineighpin .EQ. VoidCell) THEN
            neighbeta = 0.5; neighphi = 0.0; albedo = 0.5
          ELSEIF (ineighpin .EQ. RefCell) THEN
            neighbeta = 0.0; neighphi = myphi; albedo = 0.0
          ENDIF
          Dtil = mybeta * neighbeta / (mybeta + neighbeta) * smy
          jfdm = - Dtil * (neighphi - myphi)
          Dhat = - (Jnet(ibd, ipin, izf, igc) - jfdm) / (myphi + neighphi)
        ENDIF
        GcPinXS(ipin_map, izf)%Dtil(ibd, igc) = Dtil
        GcPinXS(ipin_map, izf)%Dhat(ibd, igc) = Dhat
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

DEALLOCATE(Jnet)

END SUBROUTINE
! ------------------------------------------------------------------------------------------------------------
SUBROUTINE HexSetRadialGcCoupling(PinXS, GcPinXS)

USE allocs
USE PARAM,   ONLY : VOIDCELL, REFCELL, ZERO
USE geom,    ONLY : ncbd
USE TYPEDEF, ONLY : PinXS_Type

IMPLICIT NONE

TYPE(PinXS_Type),    POINTER, DIMENSION(:,:) :: PinXS, GcPinXS
TYPE(superPin_Type), POINTER, DIMENSION(:)   :: Pin

INTEGER, POINTER, DIMENSION(:) :: pinMap, pinMapRev, planeMap
REAL, POINTER, DIMENSION(:,:,:,:) :: Jnet

INTEGER :: ng, ngc, nxy, nzCMFD, ig, igc, ipin, ipin_map, ineighpin, iz, izf, ibd, jbd, iNgh, jNgh
REAL :: Dtil, Dhat, myphi, neighphi, mybeta, neighbeta, albedo, jfdm, smy
! ----------------------------------------------------

ng         = mklGeom%ng
ngc        = mklGeom%ngc
nxy        = mklGeom%nxy
nzCMFD     = mklGeom%nzCMFD
Pin       => mklGeom%superPin
pinMap    => mklGeom%pinMap
pinMapRev => mklGeom%pinMapRev
planeMap  => mklGeom%planeMap

CALL dmalloc(Jnet, ncbd, nxy, nzCMFD, ngc)
! ----------------------------------------------------
!               01. CONDENSE : current
! ----------------------------------------------------
DO igc = 1, ngc
  DO ig = mklGeom%GcStruct(1, igc), mklGeom%GcStruct(2, igc)
    DO izf = 1, nzCMFD
      iz = planeMap(izf)
      !$OMP PARALLEL PRIVATE(ipin_map, ineighpin, myphi, neighphi, Dtil, Dhat, jfdm, iNgh)
      !$OMP DO SCHEDULE(GUIDED)
      DO ipin = 1, nxy
        ipin_map = pinMap(ipin)
        myphi    = mklCMFD%phis(ipin, izf, ig)
        
        DO iNgh = 1, Pin(ipin_map)%nNgh
          ineighpin = Pin(ipin_map)%Neighidx(iNgh)
          ineighpin = pinMapRev(ineighpin)
          
          IF (ineighpin .LE. 0) THEN
            neighphi = ZERO
          ELSE
            neighphi = mklCMFD%phis(ineighpin, izf, ig)
          END IF
          
          Dtil = PinXS(ipin_map, iz)%Dtil(iNgh, ig)
          Dhat = PinXS(ipin_map, iz)%Dhat(iNgh, ig)
          jfdm = -Dtil * (neighphi - myphi) - Dhat * (neighphi + myphi)
          
          Jnet(iNgh, ipin, izf, igc) = Jnet(iNgh, ipin, izf, igc) + jfdm
        END DO
      END DO
      !$OMP END DO
      !$OMP END PARALLEL
    END DO
  END DO
END DO
! ----------------------------------------------------
!               02. CAL : dhat, dtil
! ----------------------------------------------------
!$OMP PARALLEL PRIVATE(ipin_map, ineighpin, ibd, jbd, myphi, neighphi, mybeta, neighbeta, Dtil, Dhat, jfdm, albedo, smy, iNgh, jNgh)
!$OMP DO SCHEDULE(GUIDED) COLLAPSE(3)
DO igc = 1, ngc
  DO izf = 1, nzCMFD
    DO ipin = 1, nxy
      ipin_map = pinMap(ipin)
      
      DO iNgh = 1, Pin(ipin_map)%nNgh
        ibd       = Pin(ipin_map)%NghBd(iNgh)
        ineighpin = Pin(ipin_map)%NeighIdx(iNgh)
        jNgh      = Pin(ipin_map)%NeighSurfIdx(iNgh)
        smy       = Pin(ipin_map)%BdLength(iNgh)
        
        myphi  = GcPinXS(ipin_map, izf)%Phi(igc)
        mybeta = GcPinXS(ipin_map, izf)%XSD(igc) / Pin(ipin_map)%Center2SurfaceL(ibd)
        
        IF (ineighpin .GT. 0) THEN
          jbd = Pin(ineighpin)%NghBd(jNgh)
          
          neighphi  = GcPinXS(ineighpin, izf)%Phi(igc)
          neighbeta = GcPinXS(ineighpin, izf)%XSD(igc) / Pin(ineighpin)%Center2SurfaceL(jbd)
          
          Dtil = mybeta * neighbeta / (mybeta + neighbeta) * smy
          jfdm = -Dtil * (neighphi - myphi)
          Dhat = -(Jnet(iNgh, ipin, izf, igc) - jfdm) / (myphi + neighphi)
        ELSE
          IF (ineighpin .EQ. VoidCell) THEN
            neighbeta = HALF
            neighphi  = ZERO
            albedo    = HALF
          ELSE IF (ineighpin .EQ. RefCell) THEN
            neighbeta = ZERO
            neighphi  = myphi
            albedo    = ZERO
          END IF
          
          Dtil = mybeta * neighbeta / (mybeta + neighbeta) * smy
          jfdm = -Dtil * (neighphi - myphi)
          Dhat = -(Jnet(iNgh, ipin, izf, igc) - jfdm) / (myphi + neighphi)
        END IF
        
        GcPinXS(ipin_map, izf)%Dtil(iNgh, igc) = Dtil
        GcPinXS(ipin_map, izf)%Dhat(iNgh, igc) = Dhat
      END DO
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

DEALLOCATE(Jnet)
! ----------------------------------------------------

END SUBROUTINE HexSetRadialGcCoupling
! ------------------------------------------------------------------------------------------------------------

END MODULE MKL_HOMOXS
#endif