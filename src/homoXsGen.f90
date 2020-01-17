SUBROUTINE HomoXsGen(Core, FXR, Phis, PinXS, myzb, myze, ng, lXsLib, lScat1, lsigT)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,    FXRInfo_Type,      PinXs_Type,        &
                           Pin_Type,         PinInfo_Type,      Cell_Type
USE CMFD_mod,       ONLY : XsMac,                                                  &
                           HomoCellXsGen
USE Core_mod,       ONLY : GroupInfo
USE BenchXs,        ONLY : XsBaseBen
USE MacXsLib_Mod,   ONLY : MacXsBase,      MacXsScatMatrix
USE BasicOperation, ONLY : CP_VA,            CP_CA,             MULTI_VA
USE Th_Mod,         ONLY : GetPinFuelTemp,   GetPinModTemp,     GetPinTemp

IMPLICIT NONE

TYPE(CoreInfo_Type), INTENT(IN) :: Core
TYPE(FxrInfo_Type), POINTER , INTENT(IN):: FXR(:, :)
TYPE(PinXs_Type), POINTER, INTENT(INOUT) :: PinXs(:, :)
REAL, POINTER, INTENT(IN) :: Phis(:, :, :)
INTEGER, INTENT(IN) :: myzb, myze, ng
LOGICAL, INTENT(IN) :: lXsLib, lScat1, lsigt

TYPE(FxrInfo_Type), POINTER :: myFXR
TYPE(Pin_Type), POINTER :: Pin(:)
!TYPE(AsyInfo_Type), POINTER :: AsyInfo
TYPE(PinInfo_Type), POINTER :: PinInfo(:)
TYPE(Cell_Type), POINTER :: Cell(:)

INTEGER :: nCoreFxr, nCoreFsr, nxy
INTEGER :: nCellType, nPinType, nlocalFxr, nlocalFsr, nFsrInFxr
INTEGER :: icel, ipin, iz, ixy ,ifxr, ifsr, itype, ifsrlocal, ig
INTEGER :: FsrIdxSt, FxrIdxSt
INTEGER :: i, j, k, l, m
INTEGER :: iResoGrpBeg, iResoGrpEnd, norg, nChi

REAL :: XsMacsTr(ng)

Pin => Core%Pin
PinInfo => Core%Pininfo;   Cell => Core%CellInfo

nxy = Core%nxy
nCoreFxr = Core%nCoreFxr; nCoreFsr = Core%nCoreFsr

IF(lxslib) THEN
  norg = GroupInfo%norg; nChi = GroupInfo%nChi
  iResoGrpBeg = GroupInfo%nofg + 1
  iResoGrpEnd = GroupInfo%nofg + GroupInfo%norg
ENDIF
DO iz = myzb, myze
  DO ixy = 1, nxy
    FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
    !Get
    icel = Pin(ixy)%Cell(iz)
    nlocalFxr = Cell(icel)%nFxr; nlocalFsr = Cell(icel)%nFsr
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1; nFsrInFxr = Cell(icel)%nFsrInFxr(j)
      XsMac(j)%lFuel = FALSE
      IF(lXsLib) THEN
        myFxr => FXR(ifxr, iz)
        CALL MacXsBase(XsMac(j), myFxr, 1, ng, ng, 1._8, FALSE, TRUE, TRUE)
        CALL MacXsScatMatrix(XsMac(j), myFxr, 1, ng, ng, GroupInfo, FALSE, TRUE)
        !Self-Sheilding Effect
        IF(myFxr%lres) THEN
           do ig = iResoGrpBeg, iResoGrpEnd 
             XsMac(j)%XsMacA(ig) = XsMac(j)%XsMacA(ig) * myFxr%FresoA(ig)  
             XsMac(j)%XsMacNf(ig) = XsMac(j)%XsMacNf(ig) * myFxr%FresoF(ig)  
             XsMac(j)%XsMacKf(ig) = XsMac(j)%XsMacKf(ig) * myFxr%FresoF(ig)  
           enddo
        ENDIF
        XsMac(j)%XsMacTr = XsMac(j)%XsMacA + XsMac(j)%XsMacStr
        XsMac(j)%XsMacT = XsMac(j)%XsMacA + XsMac(j)%XsMacS
        !Obtaining
#ifdef inflow
        DO ig = 1, ng
          XsMac(j)%XsMacTr(ig) = XsMac(j)%XsMacTr(ig) + myFxr%DelInflow(ig)
          XsMac(j)%XsMacSm(ig, ig) = XsMac(j)%XsMacSm(ig, ig) + myFxr%DelInflow(ig)
        ENDDO
#endif
        CALL CP_CA(XsMac(j)%CHI, 0._8, ng)
        IF(myFxr%lDepl) THEN
          CALL CP_VA(XsMac(j)%CHI(1:nCHI), myFxr%CHI(1:nCHI), nCHI)
          XsMac(j)%lFuel = TRUE
        ENDIF
      ELSE
        ifsrlocal = Cell(icel)%MapFxr2FsrIdx(1,j)
        itype=Fxr(ifxr,iz)%imix
        !itype = Cell(icel)%iReg(ifsrlocal)
        !CALL xsbaseBen(itype, 1, ng, 1, ng, .FALSE., XsMac(j))  !--r562 original
        CALL xsbaseBen(itype, 1, ng, 1, ng, lscat1, XsMac(j))    !--r554 old >> r562c 17/02/05
      ENDIF
    ENDDO !Fxr sweep in the Cell
    CALL HomoCellXsGen(Cell(icel), Phis(FsrIdxSt:FsrIdxSt+nlocalFsr-1, iz, 1:ng), PinXS(ixy, iz),  &
                       XsMac(1:nLocalFxr), ng, nLocalFxr, nLocalFsr, GroupInfo%OutScatRange, lXsLib, lsigt)
  ENDDO
ENDDO

IF(lXsLib) THEN
  DO iz = myzb, myze
    DO ixy = 1, nxy
      PinXS(ixy,iz)%PinTemp = GetPinTemp(Core, Fxr, iz, ixy)
      IF(Core%lFuelPlane(iz) .AND. Pin(ixy)%lFuel) THEN
        PinXS(ixy,iz)%FuelTemp = GetPinFuelTemp(Core, Fxr, iz, ixy)
        PinXS(ixy,iz)%ModTemp = GetPinModTemp(Core, Fxr, iz, ixy)
      ENDIF
    ENDDO
  ENDDO
ENDIF
!Finalize
NULLIFY(Pin)
NULLIFY(PinInfo)
NULLIFY(Cell)

END SUBROUTINE

SUBROUTINE HomoCellXsGen(CellInfo, Phis, PinXs, XsMac, ng, nFxr, nFsr, OutScatRange, lXsLib, lsigT)
USE PARAM
USE TYPEDEF, ONLY : Cell_Type, PinXs_Type, XsMac_Type
USE BasicOperation, ONLY : CP_CA, CP_VA, MULTI_CA
USE CNTL,   ONLY : nTracerCntl
IMPLICIT NONE

TYPE(Cell_Type) :: CellInfo
REAL :: phis(:, :)
TYPE(PinXs_Type) :: PinXs
TYPE(XsMac_Type) :: XsMac(nFxr)
REAL :: phisum, vol, volsum, localphi, localpsi, Rphisum,RR(ng + 6), CHI(ng), scatsum !BYSedit
REAL :: onethree, threeseven
INTEGER, POINTER :: OutScatRange(:, :)
INTEGER :: ng, nFXR, nFsr, nFsrInFxr
LOGICAL :: lXsLib
INTEGER :: i, j, k, ig, ig2, igb, ige, ireg
LOGICAL :: lfuel, lsigt
ONETHREE = one/3._8
THREESEVEN = 3._8/7._8
DO ig = 1, ng
  PinXS%XSS(ig)%FROM = ZERO
ENDDO
DO ig = 1, ng
  phisum = 0; volsum = 0
  igb = OutScatRange(1, ig)
  ige = OutScatRange(2, ig)
  CALL CP_CA(RR, ZERO, ng + 6)     !---BYS edit // ng+4 > ng+5 // for calc. SIG_a  // ng+6 // for calc. 1/3sigtr
  DO i = 1, nFxr
    nFsrInFxr = CellInfo%nFsrInFxr(i)
    DO j = 1, nFsrInFxr
      ireg = CellInfo%MapFxr2FsrIdx(j, i)
      vol = CellInfo%vol(ireg)
      
      IF (vol < 1E-10) CYCLE
      
      volsum = volsum + vol
      localphi = phis(ireg, ig) * vol
      phisum = phisum + localphi

      RR(1) = RR(1) + localphi * XsMac(i)%Xsmact(ig)
      RR(2) = RR(2) + localphi * XsMac(i)%Xsmactr(ig)
      RR(3) = RR(3) + localphi * XsMac(i)%Xsmacnf(ig)
      RR(4) = RR(4) + localphi * XsMac(i)%Xsmackf(ig)

      DO ig2 = igb, ige
        RR(4 + ig2) = RR(4 + ig2) + localphi * XsMac(i)%Xsmacsm(ig, ig2)
      ENDDO
      !---BYS edit
      RR(ng+5) = RR(ng+5) + localphi * XsMac(i)%Xsmaca(ig)
      IF( lsigT )THEN
        RR(ng+6) = RR(ng+6) + localphi / XsMac(i)%Xsmact(ig)
      ELSE
        RR(ng+6) = RR(ng+6) + localphi / XsMac(i)%Xsmactr(ig)
      ENDIF
      !---BYS edit end
    ENDDO
  ENDDO
  PinXS%phi(ig) = phisum / volsum

  Rphisum = one/phisum
  CALL MULTI_CA(Rphisum, RR(:), ng + 6)
  PinXS%XST(ig) = RR(1); PinXS%XSTR(ig) = RR(2);
  PinXS%XSNF(ig) = RR(3); PinXS%XSKF(ig) = RR(4);
  !---BYS edit
  PinXS%XSA(ig)=RR(ng+5)
  !---BYS edit end
  PinXS%Xss(ig)%WithInGroupScat = RR(4 + ig)
  DO ig2 = igb, ige
    IF((ig-PinXS%Xss(ig2)%ib)*(ig-PinXS%Xss(ig2)%ie) .GT. 0) CYCLE
    PinXS%Xss(ig2)%From(ig) = RR(4+ ig2)
  ENDDO
  !---BYS edit
  scatsum=0;
  DO ig2 = 1, ng
      scatsum=scatsum + RR(4+ig2)
  ENDDO
  PinXS%XSR(ig) = RR(ng+5) + scatsum - RR(4+ig)
  !---BYS edit end
  PinXS%Xss(ig)%From(ig) = 0
  PinXs%XSD(ig) = ONETHREE/PinXS%XSTR(ig)
  PinXS%XSD2(ig) = THREESEVEN/PinXS%XST(ig)

  !PinXS%XSR(ig) = PinXS%XSTR(ig) - PinXS%Xss(ig)%WithInGroupScat  ! BYS edit 16/06/29 (original r543)
  !PinXS%XSR(ig) = PinXS%XST(ig) - PinXS%Xss(ig)%WithInGroupScat  ! BYS edit 16/06/29
  !PinXS%XSD2(ig) = THREESEVEN/PinXS%XSTR(ig)
  !PinXs%XSD(ig) = ONETHREE/PinXS%XST(ig)
  !PinXS%XSD2(ig) = THREESEVEN/PinXS%XST(ig)
  IF( lsigT )THEN
      PinXS%XSTR(ig)=PinXS%XST(ig)
    PinXs%XSD(ig) = ONETHREE/PinXS%XST(ig)
  ENDIF
  IF( nTracerCntl%lDhom )THEN  ! homogenized by D
    PinXs%XSD(ig) = ONETHREE*RR(ng+6)
    !PinXs%XSD(ig) = PinXs%XSD(ig)*0.90_8
  ENDIF

ENDDO

!CHI Update
phisum = zero
CALL CP_CA(CHI, ZERO, ng)
lfuel=.FALSE.
IF(.NOT. lxslib )THEN
    DO i = 1, nFXR
        IF( XsMac(i)%lfuel )THEN
            lfuel = .TRUE.
            EXIT
        ENDIF
    ENDDO
ELSE
    IF(CellInfo%lfuel) THEN
        lfuel=.TRUE.
    ENDIF
ENDIF
IF(lfuel) THEN
  DO i = 1, nFxr
    IF(.NOT. XsMac(i)%lfuel) CYCLE
    nFsrInFxr = CellInfo%nFsrInFxr(i)
    localpsi = 0
    DO j = 1, nFsrInFxr
      ireg = CellInfo%MapFxr2FsrIdx(j, i)
      vol = CellInfo%vol(ireg)
      
      IF (vol < 1E-10) CYCLE
      
      DO ig = 1, ng
        localphi = phis(ireg, ig) * vol * XsMac(i)%XsMacNf(ig) !Sum of fission neutron in the cell
        phisum = phisum + localphi
        localpsi = localpsi + localphi
      ENDDO
    ENDDO
    DO ig2 = 1, ng
      CHI(ig2) = CHI(ig2) + XsMac(i)%chi(ig2) * localpsi
    ENDDO
  ENDDO
  Rphisum = one/phisum
  CALL MULTI_CA(Rphisum, CHI(:), ng)
ENDIF
CALL CP_VA(PinXS%CHI(:), CHI(:), ng)
phisum = SUM(CHI)
END SUBROUTINE

!SUBROUTINE RadCouplingCoeffGen(Pin, Jout)
!USE PARAM
!USE TYPEDEF, ONLY : PinXs_TYPE, Pin_Type
!USE CMFD_MOD, ONLY : ng,        nxy,          myzb,         myze,    &
!                     hzfm,      hz,                                  &
!                     CmfdPinXS, PinNeighIdx

SUBROUTINE RadCouplingCoeffGen(Core, CmfdPinXS, Jout, ng, lDhatUpdt, PE)
USE PARAM
USE TYPEDEF,      ONLY : CoreInfo_Type,       PinXs_Type,        PE_TYPE,     &
                         Pin_Type
USE cntl,           ONLY : nTracerCntl
IMPLICIT NONE

TYPE(CoreINfo_Type) :: Core
TYPE(PinXs_Type),POINTER :: CmfdPinXS(:, :)
REAL, POINTER :: Jout(:, :, :, :, :)
TYPE(PE_TYpe) :: PE
INTEGER :: ng
LOGICAL :: lDhatUpdt


TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(PinXS_Type), POINTER :: PinXS
REAL, POINTER :: hzfm(:), hz(:)
INTEGER :: myzb, myze, nxy

INTEGER :: ixy, ineigh, ig, iz, ibd, inbd
INTEGER, PARAMETER :: nbd = 4
REAL :: myphi(ng), neighphi, phisum, jfdm, jnet
REAL :: PDHAT, Del, Alpha
REAL :: Dhat, Dtil, mybeta, neighbeta, smy
REAL :: atil, ahat, surfphifdm   !--- CNJ Edit : Domain Decomposition
REAL :: NineSeven
LOGICAL :: lDhatCor

! lDhatCor=.FALSE.
lDhatCor=.TRUE.

Pin => Core%Pin
hzfm => Core%hzfm; hz => Core%hz
myzb = PE%myzb; myze = PE%myze
nxy = core%nxy
NineSeven = 9._8/7._8

DO iz = myzb, myze
  DO ixy = 1, nxy
    PinXS => CmfdPinXS(ixy, iz)
    myphi(1:ng) = PinXS%phi(1:ng)
    DO ibd = 1, nbd
      ineigh = Pin(ixy)%NeighIdx(ibd)
      inbd = Pin(ixy)%NeighSurfIdx(ibd)    !The Corresponding Surface Index of
      smy = Pin(ixy)%BdLength(ibd) !* hz(iz)
      IF(ineigh .GT. 0) THEN
        DO ig = 1, ng
          neighphi = CmfdPinXs(ineigh, iz)%phi(ig); phisum = neighphi + myphi(ig)
          mybeta = CmfdPinXS(ixy, iz)%XSD(ig) / Pin(ixy)%Center2SurfaceL(ibd)
          neighbeta = CmfdPinXS(ineigh, iz)%XSD(ig) / Pin(ineigh)%Center2SurfaceL(inbd)
          Dtil = mybeta * neighbeta/(mybeta + neighbeta) * smy
          PinXs%dtil(ibd, ig) = dtil
          !--- CNJ Edit : Domain Decomposition
          atil = mybeta / (mybeta + neighbeta) * smy
          surfphifdm = atil * myphi(ig) + (smy - atil) * neighphi
          ahat = (Jout(3, ibd, ixy, iz, ig) - surfphifdm) / phisum
          PinXS%atil(ibd, ig) = atil; PinXS%ahat(ibd, ig) = ahat
          IF(lDhatUpdt) THEN
            jfdm = Dtil * (myphi(ig) - neighphi)
            jnet = (Jout(2, ibd, ixy, iz, ig) - Jout(1, ibd, ixy, iz, ig)) !* hz(iz)
            dhat = (jfdm - jnet) / phisum
            pDhat = PinXs%PDHAT(ibd, ig)
            Del = ABS(dhat - pdhat)
            IF(Del .GT. 10.*Dtil .AND. lDhatCor) THEN
              Alpha = Dtil/(Del-Dtil)
              Dhat = PDHAT + Alpha * (Dhat - PDhat)
            ENDIF
            PinXs%PDHAT(ibd, ig) = Dhat; PinXs%dhat(ibd, ig) = dhat
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
          IF (nTracerCntl%lGroupAlbedo) THEN
            neighbeta = (1.0 - Core%groupAlbedo(ibd, ig)) / (2.0 * (Core%groupAlbedo(ibd, ig) + 1.0))
          ENDIF
          Dtil = mybeta * neighbeta/(mybeta + neighbeta) * smy
          PinXs%dtil(ibd, ig) = dtil
          !--- CNJ Edit : Domain Decomposition
          atil = mybeta / (mybeta + neighbeta) * smy
          surfphifdm = atil * myphi(ig) + (smy - atil) * neighphi
          ahat = (Jout(3, ibd, ixy, iz, ig) - surfphifdm) / phisum
          PinXS%atil(ibd, ig) = atil; PinXS%ahat(ibd, ig) = ahat
          IF(lDhatUpdt) THEN
            jfdm = Dtil * (myphi(ig) - neighphi)
            jnet = (Jout(2, ibd, ixy, iz, ig) - Jout(1, ibd, ixy, iz, ig)) !* hz(iz)
            dhat = (jfdm - jnet) / phisum
            PinXs%PDHAT(ibd, ig) = Dhat; PinXs%dhat(ibd, ig) = dhat
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

SUBROUTINE UpdatePhiC(PinXS, PhiC)
USE PARAM
USE TYPEDEF,  ONLY : PinXs_Type
USE CMFD_Mod, ONLY : myzb, myze, nxy, ng
IMPLICIT NONE
TYPE(PinXs_Type), POINTER :: PinXs(:, :)
REAL, POINTER :: PhiC(:, :, :)
INTEGER :: iz, ig, ixy

DO iz = myzb, myze
  DO ixy = 1, nxy
    PhiC(ixy, iz, 1:ng) = PinXS(ixy, iz)%phi(1:ng)
     !PinXS(ixy, iz)%phi(1:ng) = PinXS(ixy, iz)%phi(1:ng)
  ENDDO
ENDDO
END SUBROUTINE
