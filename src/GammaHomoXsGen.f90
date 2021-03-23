#include <defines.h>
#ifdef __GAMMA_TRANSPORT
!-- JSU EDIT 20170803
SUBROUTINE GamHomoXsGen(Core, FXR, Phis, GPhis, PinXS, myzb, myze, ng, ngg, lScat1, lsigT)
! GENERATE HOMOGENIZED CROSS SECTION
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,    FXRInfo_Type,      Pin_Type,        Cell_Type
USE BasicOperation, ONLY : CP_VA,            CP_CA,             MULTI_VA
USE Th_Mod,         ONLY : GetPinFuelTemp,   GetPinModTemp,     GetPinTemp

USE GAMCMFD_mod,    ONLY : GXSMac,           GamCellXsGen
USE GamXsLib_Mod,   ONLY : GamXsBase,        GamScatMatrix,     GamProdMatrix
USE GammaTYPEDEF,   ONLY : GPINXS_TYPE
USE Core_mod,       ONLY : GroupInfo

IMPLICIT NONE

TYPE(CoreInfo_Type), INTENT(IN) :: Core
TYPE(FxrInfo_Type), POINTER , INTENT(IN):: FXR(:, :)
TYPE(GPINXS_TYPE), POINTER, INTENT(INOUT) :: PinXs(:, :)
REAL, POINTER, INTENT(IN) :: Phis(:, :, :), GPhis(:, :, :)
INTEGER, INTENT(IN) :: myzb, myze, ng, ngg
LOGICAL, INTENT(IN) :: lScat1, lsigt

! Pointing Variables
TYPE(FxrInfo_Type), POINTER :: myFXR
TYPE(Pin_Type), POINTER :: Pin(:)
!TYPE(AsyInfo_Type), POINTER :: AsyInfo
TYPE(Cell_Type), POINTER :: Cell(:)
REAL, POINTER :: gphiIn(:,:)

INTEGER :: nCoreFxr, nCoreFsr, nxy
INTEGER :: nCellType, nPinType, nlocalFxr, nlocalFsr, nFsrInFxr
INTEGER :: icel, ipin, iz, ixy ,ifxr, ifsr, itype, ifsrlocal, ig
INTEGER :: FsrIdxSt, FxrIdxSt
INTEGER :: i, j, k, l, m, igg
REAL :: XsMacsTr(ng)

Pin => Core%Pin
Cell => Core%CellInfo

nxy = Core%nxy
nCoreFxr = Core%nCoreFxr; nCoreFsr = Core%nCoreFsr


DO iz = myzb, myze
  DO ixy = 1, nxy
    ! Index Preparation
    FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
    icel = Pin(ixy)%Cell(iz)
    nlocalFxr = Cell(icel)%nFxr; nlocalFsr = Cell(icel)%nFsr
    ALLOCATE(gphiIn(nLocalFSR, ngg))
    DO j = 1, nLocalFxr
      ! Index Preparation
      ifxr = FxrIdxSt + j -1; nFsrInFxr = Cell(icel)%nFsrInFxr(j)
      myFxr => FXR(ifxr, iz)
      ! Macroscopic XS and matrices making
      CALL GamXsBase(GXSMac(j), myFxr, 1, ngg, ngg, FALSE, TRUE)
      CALL GamProdMatrix(GXSMac(j), myFxr, 1, ngg, ng, ngg, GroupInfo, FALSE)
      CALL GamScatMatrix(GXSMac(j), myFxr, 1, ngg, ngg, TRUE, TRUE, FALSE)

      GXSMac(j)%XsMacTr = GXSMac(j)%XsMacA + GXSMac(j)%XsMacStr
      GXSMac(j)%XsMacT = GXSMac(j)%XsMacA + GXSMac(j)%XsMacS
    ENDDO !Fxr sweep in the Cell
    IF (lscat1) THEN  ! Manipulate Index of gphis
      gphiIn(:, :) = GPhis(FsrIdxSt:FsrIdxSt+nlocalFsr-1, 1:ngg, iz)
    ELSE
      DO igg = 1, ngg
        gphiIn(:,igg) = GPhis(igg, FsrIdxSt:FsrIdxSt+nlocalFsr-1, iz)
      END DO
    END IF
    CALL GamCellXsGen(Cell(icel), Phis(FsrIdxSt:FsrIdxSt+nlocalFsr-1, iz, 1:ng),                    &
          gphiIn, PinXS(ixy, iz), GXSMac(1:nLocalFxr), ng, ngg, nLocalFxr, nLocalFsr,               &
          GroupInfo%OutScatRange_Ph, lsigt)
    DEALLOCATE(gphiIn)
  ENDDO
ENDDO

DO iz = myzb, myze
  DO ixy = 1, nxy
    PinXS(ixy,iz)%PinTemp = GetPinTemp(Core, Fxr, iz, ixy)
    IF(Core%lFuelPlane(iz) .AND. Pin(ixy)%lFuel) THEN
      PinXS(ixy,iz)%FuelTemp = GetPinFuelTemp(Core, Fxr, iz, ixy)
      PinXS(ixy,iz)%ModTemp = GetPinModTemp(Core, Fxr, iz, ixy)
    ENDIF
  ENDDO
ENDDO
!Finalize
NULLIFY(Pin)
NULLIFY(Cell)

END SUBROUTINE

SUBROUTINE GamCellXsGen(CellInfo, Phis, GPhis, PinXs, GXSMac, ng, ngg, nFxr, nFsr, OutScatRange,    &
                             lsigT)
! USING MACROSCOPIC XS IN FXR'S, GENERATE CELL HOMOGENIZED CROSS SECTION
USE PARAM
USE BasicOperation,   ONLY : CP_CA, CP_VA, MULTI_CA
USE CNTL,             ONLY : nTracerCntl
USE TYPEDEF,          ONLY : Cell_Type
USE GammaTYPEDEF,     ONLY : GamMacXS_TYPE, GPINXS_TYPE
IMPLICIT NONE
! INPUT VARIABLES
TYPE(Cell_Type) :: CellInfo
REAL :: phis(:, :), GPhis(:, :)
TYPE(GPINXS_TYPE) :: PinXs
TYPE(GamMacXS_TYPE) :: GXSMac(nFxr)
INTEGER :: ng, ngg
INTEGER :: nFXR, nFsr
INTEGER, POINTER :: OutScatRange(:, :)
LOGICAL :: lsigt

! LOCAL VARIABLES
REAL :: phisum, vol, volsum, localphi, Rphisum, RR(ngg + 5), scatsum !BYSedit
REAL :: gphisum, localgphi, RGphisum   ! GAMMA VARIABLES
REAL :: onethree, threeseven
INTEGER :: nFsrInFxr
INTEGER :: i, j, k, ig, ig2, igb, ige, ireg, igg, igg2
LOGICAL :: lfuel
ONETHREE = one/3._8
THREESEVEN = 3._8/7._8

DO igg = 1, ngg
  PinXS%XSS(igg)%FROM = ZERO
ENDDO

! for photo-atomic reaction data (Homogenization Using Gamma Flux)
! PHOTON REACTION XS'S, WEIGHTING FUNCTION : PHOTON FLUX  "gphis"
DO igg = 1, ngg  ! PHOTON ENERGY GROUP
  volsum = 0.
  gphisum = 0.
  igb = OutScatRange(1, igg)
  ige = OutScatRange(2, igg)
  RR = 0.   ! RR : 1-TOTAL, 2-TRANSPORT, 3-KERMA, 3+1~3+ngg-SCAT, ngg+4-ABSORPTION, ngg+5-inverse
  DO i = 1, nFxr
    nFsrInFxr = CellInfo%nFsrInFxr(i)
    DO j = 1, nFsrInFxr
      ireg = CellInfo%MapFxr2FsrIdx(j, i)
      vol = CellInfo%vol(ireg)
      volsum = volsum + vol
      localgphi = gphis(ireg, igg) * vol
      gphisum = gphisum + localgphi         ! volume-weighted summation of multigroup photon flux
      ! Reaction rate and KERMA photo-atomic reactions
      RR(1) = RR(1) + localgphi * GXSMac(i)%Xsmact(igg)    ! RR(1) : Macro total
      RR(2) = RR(2) + localgphi * GXSMac(i)%Xsmactr(igg)   ! RR(2) : Macro transport
      RR(3) = RR(3) + localgphi * GXSMac(i)%MacKERMA(igg)  ! RR(3) : Macro KERMA 
      !RR(4~3+ngg) : scattering matrix summation (scattered source)
      DO igg2 = igb, ige
        RR(3 + igg2) = RR(3 + igg2) + localgphi * GXSMac(i)%Xsmacsm(igg, igg2)    ! igg -> igg2
      ENDDO
      RR(ngg+4) = RR(ngg+4) + localgphi * GXSMac(i)%Xsmaca(igg) ! RR(ngg+4) : absorption XS
      IF( lsigT )THEN ! RR(ngg+5) : diffusion coefficient?
        RR(ngg+5) = RR(ngg+5) + localgphi / GXSMac(i)%Xsmact(igg)
      ELSE
        RR(ngg+5) = RR(ngg+5) + localgphi / GXSMac(i)%Xsmactr(igg)
      ENDIF
    ENDDO ! FSR loop
  ENDDO ! FXR loop
  PinXS%Gphi(igg) = gphisum / volsum
  RGphisum = one / gphisum
  CALL MULTI_CA(RGphisum, RR(:), ngg + 5)
  PinXS%XST(igg) = RR(1)
  PinXS%XSTR(igg) = RR(2)
  !PinXS%KERMA(igg) = RR(3)                            ! Pin-homogenized KERMA
  PinXS%XSA(igg)=RR(ngg + 4)
  PinXS%Xss(igg)%self = RR(3 + igg)        ! igg -> igg
  DO igg2 = igb, ige
    IF((igg-PinXS%Xss(igg2)%ib)*(igg-PinXS%Xss(igg2)%ie) .GT. 0) CYCLE
    PinXS%Xss(igg2)%From(igg) = RR(3+ igg2)           ! igg -> igg2
  ENDDO
  scatsum=0;
  DO igg2 = 1, ngg
    scatsum = scatsum + RR(3 + igg2)                   ! sum of igg -> igg2 (all scattering from igg)
  ENDDO !
  PinXS%XSR(igg) = RR(ngg + 4) + scatsum - RR(3 + igg) ! Removal XS -- absorption + scattering - self scattering(igg+3)
  PinXS%XSS(igg)%From(igg) = 0._8                      ! Self Scattering with removal (zero)
  PinXs%XSD(igg) = ONETHREE / PinXS%XSTR(igg)          ! Diffusion Cefficient

  IF( lsigT )THEN
    PinXS%XSTR(igg) = PinXS%XST(igg)
    PinXs%XSD(igg)  = ONETHREE/PinXS%XST(igg)
  ENDIF

  IF( nTracerCntl%lDhom )THEN  ! homogenized by D
    PinXs%XSD(igg) = ONETHREE * RR(ngg + 5)
    !PinXs%XSD(ig) = PinXs%XSD(ig)*0.90_8
  ENDIF
ENDDO  ! PHOTON GROUP LOOP (igg) end

! for photon generation reaction data (Homogenization Using Neutron Flux)
! PHOTON PRODUCTION MATRIX, WEIGHTING FUNCTION : NEUTRON FLUX "phis"
PinXS%XSP = ZERO  ! production
phisum = 0.;
DO ig = 1, ng
  volsum = 0.
  phisum = 0.
  RR = 0.  ! Production,  RR(igg) : photon generated within E group igg
  DO i = 1, nFxr
    nFsrInFxr = CellInfo%nFsrInFxr(i)
    DO j = 1, nFsrInFxr
      ireg = CellInfo%MapFxr2FsrIdx(j, i)
      vol = CellInfo%vol(ireg)
      volsum = volsum + vol
      localphi = phis(ireg, ig) * vol
      phisum = phisum + localphi
      ! Production from (ig) to (igg), contained at RR(igg)
      DO igg = 1, ngg
        RR(igg) = RR(igg) + localphi * GXSMac(i)%GProdTot(ig, igg)   ! ig (neutron) -> igg (photon)
      ENDDO
    END DO ! FSR loop
  END DO ! FXR loop
  PinXS%Nphi(ig) = phisum / volsum
  Rphisum = one / phisum
  RR = RR * Rphisum
  PinXS%XSP(ig,:) = RR(1:ngg)
END DO ! Neutron Energy Group loop

END SUBROUTINE

!
SUBROUTINE RadCouplingCoeffGen_Gamma(Core, CmfdPinXS, Jout, ng, lDhatUpdt, lScat1, PE)
USE PARAM
USE TYPEDEF,      ONLY : CoreInfo_Type,       Pin_Type,        PE_TYPE
USE GammaTYPEDEF, ONLY : GPinXs_Type
USE cntl,           ONLY : nTracerCntl
IMPLICIT NONE

TYPE(CoreINfo_Type) :: Core
TYPE(GPinXs_Type),POINTER :: CmfdPinXS(:, :)
REAL, POINTER :: Jout(:, :, :, :, :)
TYPE(PE_TYpe) :: PE
INTEGER :: ng   ! Photon Energy Group
LOGICAL :: lDhatUpdt, lScat1

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(GPinXs_Type), POINTER :: PinXS
REAL, POINTER :: hzfm(:), hz(:)
INTEGER :: myzb, myze, nxy

INTEGER :: ixy, ineigh, ig, iz, ibd, inbd
INTEGER, PARAMETER :: nbd = 4  ! RADIAL
REAL :: myphi(ng), neighphi, phisum, jfdm, jnet
REAL :: PDHAT, Del, Alpha
REAL :: Dhat, Dtil, mybeta, neighbeta, smy
REAL :: atil, ahat, surfphifdm   !--- CNJ Edit : Domain Decomposition
REAL :: NineSeven
LOGICAL :: lDhatCor

lDhatCor=.FALSE.

Pin => Core%Pin
hzfm => Core%hzfm; hz => Core%hz
myzb = PE%myzb; myze = PE%myze
nxy = Core%nxy
NineSeven = 9._8/7._8

DO iz = myzb, myze
  DO ixy = 1, nxy
    PinXS => CmfdPinXS(ixy, iz)
    myphi(1:ng) = PinXS%gphi(1:ng)
    DO ibd = 1, nbd
      ineigh = Pin(ixy)%NeighIdx(ibd)
      inbd = Pin(ixy)%NeighSurfIdx(ibd)    !The Corresponding Surface Index of
      smy = Pin(ixy)%BdLength(ibd) !* hz(iz)
      IF(ineigh .GT. 0) THEN
        DO ig = 1, ng
          neighphi = CmfdPinXs(ineigh, iz)%gphi(ig); phisum = neighphi + myphi(ig)
          mybeta = CmfdPinXS(ixy, iz)%XSD(ig) / Pin(ixy)%Center2SurfaceL(ibd)
          neighbeta = CmfdPinXS(ineigh, iz)%XSD(ig) / Pin(ineigh)%Center2SurfaceL(inbd)
          Dtil = mybeta * neighbeta/(mybeta + neighbeta) * smy
          PinXs%dtil(ibd, ig) = dtil
          IF(lDhatUpdt) THEN
            jfdm = Dtil * (myphi(ig) - neighphi)
            IF (lScat1) THEN
              jnet = (Jout(2, ibd, ixy, ig, iz) - Jout(1, ibd, ixy, ig, iz)) !* hz(iz)
            ELSE
              jnet = (Jout(2, ig, ibd, ixy, iz) - Jout(1, ig, ibd, ixy, iz)) !* hz(iz)
            ENDIF
            dhat = (jfdm - jnet) / phisum

            pDhat = PinXs%PDHAT(ibd, ig)
            Del = ABS(dhat - pdhat)
            IF(Del .GT. 10.*Dtil .AND. lDhatCor) THEN
              Alpha = Dtil/(Del-Dtil)
              Dhat = PDHAT + Alpha * (Dhat - PDhat)
            ENDIF
            PinXs%PDHAT(ibd, ig) = Dhat; PinXs%dhat(ibd, ig) = dhat;
          ENDIF
        ENDDO  !ENd of Group Sweep
      ELSE     !Boundary
        IF(ineigh .EQ. VoidCell) THEN
          neighbeta = 0.5_8
        ELSE
          neighbeta = 0
        ENDIF
        DO ig = 1, ng
          neighphi = 0; phisum = neighphi + myphi(ig)
          mybeta = CmfdPinXS(ixy, iz)%XSD(ig) / Pin(ixy)%Center2SurfaceL(ibd)
          Dtil = mybeta * neighbeta/(mybeta + neighbeta) * smy
          PinXs%dtil(ibd, ig) = dtil
          IF(lDhatUpdt) THEN
            jfdm = Dtil * (myphi(ig) - neighphi)
            IF (lScat1) THEN
              jnet = (Jout(2, ibd, ixy, ig, iz) - Jout(1, ibd, ixy, ig, iz)) !* hz(iz)
            ELSE
              jnet = (Jout(2, ig, ibd, ixy, iz) - Jout(1, ig, ibd, ixy, iz)) !* hz(iz)
            ENDIF
            dhat = (jfdm - jnet) / phisum
            PinXs%PDHAT(ibd, ig) = Dhat; PinXs%dhat(ibd, ig) = dhat;
          ENDIF

        ENDDO
      ENDIF
       IF(.NOT. core%lfuelplane(iz) .AND. nTRACErCntl%lAxRefFDM) THEN
         DO ig= 1, ng
            PinXs%PDHAT(ibd, ig) = 0; PinXs%dhat(ibd, ig) = 0;
         ENDDO
       ENDIF
    ENDDO  !End of Neighborhood Sweep
    NULLIFY(PinXs)
  ENDDO  !End of Radial Pin Sweep
ENDDO
NULLIFY(Pin)
NULLIFY(hzfm)
NULLIFY(PINXS)
END SUBROUTINE
  
#endif
