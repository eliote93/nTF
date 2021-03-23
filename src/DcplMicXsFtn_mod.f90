MODULE DcplMicXsFtn_mod
USE PARAM
USE TYPEDEF,    ONLY : MicXsFtn_Type,     FxrInfo_Type,       PinXs_Type,         &
                       Cell_Type
IMPLICIT NONE

CONTAINS

SUBROUTINE UpdtMacXs_MicXsFtn(PinXs, MicXsFtn, Fxr, CellInfo, nFxr, GroupInfo)
USE PARAM
USE TYPEDEF,    ONLY : PinXS_TYPE,        MicXsFtn_TYPE,        FxrInfo_TYPE,          &
                       Cell_TYPE,         GroupInfo_TYPE
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
TYPE(PinXS_TYPE) :: PinXS
TYPE(MicXsFtn_TYPE) :: MicXsFtn
TYPE(FxrInfo_Type) :: Fxr(nFxr)
TYPE(Cell_Type) :: CellInfo
TYPE(GroupInfo_TYPE) :: GroupInfo
INTEGER :: nFxr

REAL :: RegNum(200), RegTemp(200)
REAL :: DELT(4), RefTemp(4)
REAL :: MacSig_T, MacSig_nf, MacSig_kf, Macsig_s(300)
REAL :: MicSig_T, MicSig_nf, MicSig_kf, Micsig_s(300)
INTEGER :: FxrTag(200)

INTEGER :: NG
INTEGER :: ig, ig1, ib, ie
INTEGER :: i

NG = GroupInfo%ng

CALL GetFxrTypeTag(FxrTag(1:nFxr), Fxr(1:nFxr), CellInfo, nFxr)
CALL GetRegNum(RegNum(1:nFxr), FxrTag(1:nFxr), Fxr(1:nFxr), nFxr)
CALL GetCellTemp(RegTemp(1:4), FxrTag(1:nFxr), Fxr(1:nFxr), nFxr)
RefTemp = (/MicXsFtn%Tfuel, MicXsFtn%TMod, MicXsFtn%TClad, MicXsFtn%Tst/)

DelT = RegTemp(1:4) - RefTemp(1:4)
DO ig = 1, ng
  MacSig_T = 0; MacSig_nf = 0
  MacSig_kf = 0; 
  CALL CP_CA(Macsig_s(1:ng), 0._8, ng)
  DO i = 1, 4
    IF(RegNum(i) .EQ. 0) CYCLE
    !Make Micro XS
    MicSig_T = MicXSFtn%sig_T(ig, i, 0)
    MicSig_nf = MicXSFtn%sig_nf(ig, i, 0)
    MicSig_kf = MicXSFtn%sig_kf(ig, i, 0)
    MicSig_T = MicSig_T + MicXSFtn%sig_T(ig, i, 1) * DelT(i)
    MicSig_nf = MicSig_nf + MicXSFtn%sig_nf(ig, i, 1) * DelT(i)
    MicSig_kf = MicSig_kf + MicXSFtn%sig_kf(ig, i, 1) * DelT(i)
    ib = GroupInfo%InScatRange(1, ig); ie = GroupInfo%InScatRange(2, ig);
    CALL CP_CA(Micsig_s(1:ng), 0._8, ng)
    DO ig1 = ib, ie
      Micsig_s(ig1) = MicXsFtn%xss(ig, i, 0)%from(ig1)
      Micsig_s(ig1) =  Micsig_s(ig1) + MicXsFtn%xss(ig, i, 1)%from(ig1) * DelT(i)
    ENDDO
    !Make XS set to Macro
    MacSig_T = MacSig_T + MicSig_T * RegNum(i)
    MacSig_nf = MacSig_nf + MicSig_nf * RegNum(i)
    MacSig_kf = MacSig_kf + MicSig_kf * RegNum(i)
    MacSig_s(ib:ie) = MacSig_s(ib:ie) + MicSig_s(ib:ie) * RegNum(i)
  ENDDO

  !Copy the XS
  PinXS%xst(ig) = MacSig_T;  PinXs%xstr(ig) = MacSig_T
  PinXS%xskf(ig) = MacSig_kf;PinXS%xsnf(ig) = MacSig_nf
  PinXS%xss(ig)%from(ib : ie) = MacSig_s(ib : ie)
  PinXS%Xss(ig)%WithInGroupScat = PinXs%Xss(ig)%from(ig)
  PinXs%Xss(ig)%from(ig) = 0
ENDDO

DO ig = 1, ng
  PinXS%xsr(ig) = PinXS%xstr(ig) - PinXs%Xss(ig)%WithInGroupScat
  PinXS%XSD(ig) = 1._8/(3._8 * PinXS%xstr(ig))
ENDDO
END SUBROUTINE

SUBROUTINE UpdtMicroXsFtn(Core, DcplFmInfo, MicXsFtn, iRefPln, nxy, imode, GroupInfo, PE)
USE PARAM
USE TYPEDEF,      ONLY : CoreInfo_Type           ,FmInfo_Type          ,GroupInfo_Type,    &
                         PE_TYPE                 ,Pin_Type
USE BasicOperation, ONLY : CP_CA                 ,CP_VA
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: DcplFmInfo
TYPE(GroupInfo_Type) :: Groupinfo
TYPE(MicXsFtn_Type) :: MicXsFtn(nxy)
TYPE(PE_Type) :: PE
INTEGER, INTENT(IN) :: imode, iRefPln, nxy

REAL, POINTER :: Phis(:, :, :)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)

!REAL :: LocalPhis(2000, 500)
REAL, POINTER:: LocalPhis(:, :)
REAL, POINTER :: sig_t(:, :), sig_nf(:, :), sig_kf(:, :), sig_s(:, :, :)

REAL :: RegNum(4), RegTemp(4)
REAL :: DelT(4)
INTEGER :: ng, nFxr, nFsr, nLocalFxr, nLocalFsr
INTEGER :: FxrTag(200)
INTEGER:: FsrIdxSt, FxrIdxSt
INTEGER :: ifxr, ifsr, ixy, ig, iz, icel
INTEGER :: i, ibeg, iend
LOGICAL :: lFuelPln

Phis => DcplFMInfo%Phis
Fxr => DcplFmInfo%Fxr
Pin => Core%Pin; CellInfo => Core%CellInfo
ng = GroupInfo%ng

nFxr = Core%nCoreFxr; nFsr = Core%nCoreFsr

ALLOCATE(LocalPhis(2000, ng))
ALLOCATE(sig_t(ng, 4), sig_nf(ng, 4), sig_kf(ng, 4))
ALLOCATE(sig_s(ng, ng, 4))

iz = iRefPln
lFuelPln = Core%lFuelPlane(iz)
IF(imode .eq. 2 .and. .NOT. lfuelPln) RETURN
DO ixy = 1, nxy
  FsrIdxSt = Pin(ixy)%FsrIdxSt; FxrIdxSt = Pin(ixy)%FxrIdxSt
  icel = Pin(ixy)%Cell(iz)
  nlocalFxr = CellInfo(icel)%nFxr; nlocalFsr = CellInfo(icel)%nFsr 
  ibeg = FsrIdxSt; iend = FsrIdxSt + nLocalFsr - 1
  
  DO ig = 1, ng
    CALL CP_VA(LocalPhis(1:nLocalFsr, ig), Phis(ibeg:iend, iz, ig), nLocalFsr)
  ENDDO
  CALL CP_CA(sig_t, 0._8, ng, 4); CALL CP_CA(sig_nf, 0._8, ng, 4)
  CALL CP_CA(sig_kf, 0._8, ng, 4); CALL CP_CA(sig_s, 0._8, ng, ng, 4)
  !SUBROUTINE GetFxrTypeTag(FxrTag, CellFxr, CellInfo, nFxr)
  ibeg = FxrIdxSt; iend = FxrIdxSt + nLocalFxr - 1
  CALL GetFxrTypeTag(FxrTag(1:nlocalFxr), Fxr(ibeg:iend, iz), CellInfo(icel), nlocalFxr)
  !CALL GetCellTemp(TSet, FxrTag, CellFxr, nfxr)
  CALL GetRegNum(RegNum(1:4), FxrTag(1:nlocalFxr), Fxr(ibeg:iend, iz), nLocalFxr)
  !SUBROUTINE GetCellTemp(TSet, FxrTag, CellFxr, nfxr)
  CALL GetCellTemp(RegTemp(1:4), FxrTag(1:nlocalFxr), Fxr(ibeg:iend, iz), nLocalFxr)
  CALL GetRegMicXS(Fxr(ibeg:iend, iz), LocalPhis(1:nLocalFsr, 1:ng), CellInfo(icel),     &
                   FxrTag(1:nFxr), Sig_T, Sig_Nf, Sig_kf, sig_s, nLocalFxr, nLocalFsr,   &
                   ng, lFuelPln, GroupInfo)
  
  IF(imode .EQ. 1) THEN
    DO i = 1, 4
      IF(RegNum(i) .EQ. 0._8) CYCLE
      MicXSFtn(ixy)%Sig_T(1:ng, i, 0) = sig_T(1:ng, i) / RegNum(i)
      MicXSFtn(ixy)%Sig_kf(1:ng, i, 0) = sig_kf(1:ng, i) / RegNum(i)
      MicXSFtn(ixy)%Sig_nf(1:ng, i, 0) = sig_nf(1:ng, i) / RegNum(i)
      DO ig = 1, ng
        ibeg = GroupInfo%InScatRange(1, ig); iend = GroupInfo%InScatRange(2, ig)
        MicXSFtn(ixy)%xss(ig, i, 0)%from(ibeg:iend) =  sig_s(ibeg:iend, ig, i) / RegNum(i)
      ENDDO
    ENDDO
    MicXsFtn(ixy)%Tfuel = RegTemp(1);  MicXsFtn(ixy)%TMod= RegTemp(2)
    MicXsFtn(ixy)%TClad= RegTemp(3);   MicXsFtn(ixy)%Tst = RegTemp(4)
  ELSE !(imode .EQ. 2) THEN
    DELT(1) = RegTemp(1) - MicXsFtn(ixy)%Tfuel
    DELT(2) = RegTemp(2) - MicXsFtn(ixy)%TMod
    DELT(3) = RegTemp(3) - MicXsFtn(ixy)%TClad
    DELT(4) = RegTemp(4) - MicXsFtn(ixy)%Tst
    DO i = 1, 4
      IF(imode .eq. 2 .and. i .EQ. 2) CYCLE
      IF(imode .eq. 2 .and. i .EQ. 4) CYCLE
      IF(imode .eq. 3 .and. i .EQ. 1) CYCLE
      IF(imode .eq. 3 .and. i .EQ. 3) CYCLE      
      IF(RegNum(i) .EQ. 0._8) CYCLE
      IF(DELT(i) .LT. 1.0E-4) CYCLE
      DELT(i) = 1._8 / DELT(i)
      MicXSFtn(ixy)%Sig_T(1:ng, i, 1) = sig_T(1:ng, i) / RegNum(i)
      MicXSFtn(ixy)%Sig_kf(1:ng, i, 1) = sig_kf(1:ng, i) / RegNum(i)
      MicXSFtn(ixy)%Sig_nf(1:ng, i, 1) = sig_nf(1:ng, i) / RegNum(i)
      DO ig = 1, ng
        ibeg = GroupInfo%InScatRange(1, ig); iend = GroupInfo%InScatRange(2, ig)
        MicXSFtn(ixy)%xss(ig, i, 1)%from(ibeg:iend) =  sig_s(ibeg:iend, ig, i) / RegNum(i)
      ENDDO

      MicXSFtn(ixy)%Sig_T(1:ng, i, 1) = MicXSFtn(ixy)%Sig_T(1:ng, i, 1) - MicXSFtn(ixy)%Sig_T(1:ng, i, 0)
      MicXSFtn(ixy)%Sig_kf(1:ng, i, 1) = MicXSFtn(ixy)%Sig_kf(1:ng, i, 1) - MicXSFtn(ixy)%Sig_kf(1:ng, i, 0)
      MicXSFtn(ixy)%Sig_nf(1:ng, i, 1) = MicXSFtn(ixy)%Sig_nf(1:ng, i, 1) - MicXSFtn(ixy)%Sig_nf(1:ng, i, 0)
      
      MicXSFtn(ixy)%Sig_T(1:ng, i, 1) = MicXSFtn(ixy)%Sig_T(1:ng, i, 1) * DELT(i)
      MicXSFtn(ixy)%Sig_kf(1:ng, i, 1) = MicXSFtn(ixy)%Sig_kf(1:ng, i, 1) * DELT(i)
      MicXSFtn(ixy)%Sig_nf(1:ng, i, 1) = MicXSFtn(ixy)%Sig_nf(1:ng, i, 1) * DELT(i)
      DO ig = 1, ng
        ibeg = GroupInfo%InScatRange(1, ig); iend = GroupInfo%InScatRange(2, ig)
        MicXSFtn(ixy)%xss(ig, i, 1)%from(ibeg:iend) = MicXSFtn(ixy)%xss(ig, i, 1)%from(ibeg:iend) - MicXSFtn(ixy)%xss(ig, i, 0)%from(ibeg:iend)      
        MicXSFtn(ixy)%xss(ig, i, 1)%from(ibeg:iend) = MicXSFtn(ixy)%xss(ig, i, 1)%from(ibeg:iend) * DELT(i)
      ENDDO
    ENDDO    
  ENDIF

ENDDO

DEALLOCATE(LocalPhis)
DEALLOCATE(sig_t, sig_nf, sig_kf, sig_s)
NULLIFY(Phis); NULLIFY(Fxr)
NULLIFY(Pin); NULLIFY(CellInfo)
END SUBROUTINE

SUBROUTINE GetRegMicXs(Fxr, Phis, CellInfo, FxrTag, sig_t, sig_nf, sig_kf, sig_s, nFxr, nFsr, ng, lFuelPln, GroupInfo)
!Get Region wise Micro XS
USE PARAM
USE TYPEDEF,        ONLY : GroupInfo_Type
USE CMFD_mod,       ONLY : XsMac           ,HomoCellXsGen
USE MacXsLib_Mod,   ONLY : MacXsBase       ,MacXsScatMatrix
USE BasicOperation, ONLY : CP_VA            ,CP_CA             ,MULTI_CA          &
                          ,MULTI_VA

IMPLICIT NONE
INTEGER, INTENT(IN) :: nFxr, nFsr, ng
INTEGER, INTENT(IN) :: FxrTag(nFxr)
TYPE(FxrInfo_Type), TARGET, INTENT(IN) :: Fxr(nFxr)
TYPE(Cell_TYPE), INTENT(IN) :: CellInfo
TYPE(GroupInfo_Type), INTENT(IN) :: GroupInfo
REAL, INTENT(IN) :: Phis(nFsr, ng)
REAL, INTENT(OUT) :: sig_t(ng, 4), sig_nf(ng, 4), sig_kf(ng, 4), sig_s(ng, ng, 4)
LOGICAL :: lFuelPln

TYPE(FxrInfo_Type), POINTER :: myFxr
REAL :: Area, localphi
REAL :: AreaSum(4), PhiSum(0:4), RR(ng+4, 4), RphiSum
INTEGER :: norg, nchi, nFsrInFxr
INTEGER :: ig, ig2, igb, ige, ireg, iresogrpbeg, iresogrpend
INTEGER :: i, j, k

norg = GroupInfo%norg; nCHI = GroupInfo%nCHi
iResoGrpBeg = GroupInfo%nofg + 1
iResoGrpEnd = GroupInfo%nofg + GroupInfo%norg 

DO j = 1, nFxr
  XsMac(j)%lFuel = FALSE
  myFxr => Fxr(j)
  CALL MacXsBase(XsMac(j), myFxr, 1, ng, ng, 1._8, FALSE, TRUE)
  CALL MacXsScatMatrix(XsMac(j), myFxr, 1, ng, ng, GroupInfo, FALSE)
  !Self-Sheilding Effect
  IF(myFxr%lres) THEN
    DO ig = iResoGrpBeg,iResoGrpEnd
      XsMac(j)%XsMacA(ig) = XsMac(j)%XsMacA(ig) * myFxr%FresoA(ig)
    ENDDO
    IF(lFuelPln) THEN
      DO ig = iResoGrpBeg,iResoGrpEnd
        XsMac(j)%XsMacNf(ig) = XsMac(j)%XsMacNf(ig) * myFxr%FresoNF(ig)
        XsMac(j)%XsMacKf(ig) = XsMac(j)%XsMacKf(ig) * myFxr%FresokF(ig)
      ENDDO  
    ENDIF
    XsMac(j)%XsMacTr = XsMac(j)%XsMacA + XsMac(j)%XsMacStr
    XsMac(j)%XsMacT = XsMac(j)%XsMacA + XsMac(j)%XsMacS
     
    CALL CP_CA(XsMac(j)%CHI, 0._8, ng)
    IF(myFxr%lDepl) THEN
      CALL CP_VA(XsMac(j)%CHI(1:nCHI), myFxr%CHI(1:nCHI), nCHI)
      XsMac(j)%lFuel = TRUE
    ENDIF
  ENDIF
  NULLIFY(myFxr)
ENDDO

DO ig = 1, ng
  AreaSum = 0; PhiSum = 0
  CALL CP_CA(RR, ZERO, ng + 4, 4)
  igb = GroupInfo%OutScatRange(1, ig)
  ige = GroupInfo%OutScatRange(2, ig)  
  DO i = 1, nFxr
    k = FxrTag(i)
    nFsrInFxr = CellInfo%nFsrInFxr(i)
    DO j = 1, nFsrInFxr
      ireg = CellInfo%MapFxr2FsrIdx(j, i)
      Area = CellInfo%vol(ireg)
      !Area = CellInfo%MapFxr2FsrIdx(j, i)
      AreaSum(k) = AreaSum(k) + Area
      localphi = Phis(ireg, ig) * Area
      PhiSum(k) = PhiSum(k) + localPhi
      PhiSum(0) = PhiSum(0) + localPhi
      RR(1, k) = RR(1, k) + localphi * XsMac(i)%Xsmact(ig)
      RR(2, k) = RR(2, k) + localphi * XsMac(i)%Xsmactr(ig)
      RR(3, k) = RR(3, k) + localphi * XsMac(i)%Xsmacnf(ig)
      RR(4, k) = RR(4, k) + localphi * XsMac(i)%Xsmackf(ig)
      
      DO ig2 = igb, ige
        RR(4 + ig2, k) = RR(4 + ig2, k) + localphi * XsMac(i)%Xsmacsm(ig, ig2)
      ENDDO
                  
    ENDDO
  ENDDO
  
  DO k = 1, 4
    IF(AreaSum(k) .LT. 1.E-4) CYCLE
    !RphiSum = 1._8 / PhiSum(k)
    RphiSum = 1._8 / PhiSum(0)
    CALL MULTI_CA(RphiSum, RR(1:ng+4, k), ng + 4)
    sig_t(ig, k) = RR(2, k); 
    sig_nf(ig, k) = RR(3, k); sig_kf(ig, k) = RR(4, k)
    DO ig2 = igb, ige
      sig_s(ig, ig2, k) = RR(4 + ig2, k) 
    ENDDO 
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE GetFxrTypeTag(FxrTag, CellFxr, CellInfo, nFxr)
USE TYPEDEF,   ONLY : THCell_Type
IMPLICIT NONE
TYPE(FxrInfo_Type), INTENT(IN) :: CellFxr(nfxr)
TYPE(Cell_Type), INTENT(IN) :: CellInfo
INTEGER, INTENT(IN) ::nfxr
INTEGER :: FxrTag(nFxr)

TYPE(THCell_Type), POINTER :: ThCell
INTEGER :: CLDReg(2), FuelReg(2), CoolReg(2)
INTEGER :: i1, i2, i, j
REAL :: TempSum(4), AreaSum(4)
FxrTag = 0
IF(CellInfo%lfuel) THEN
  ThCell => CellInfo%ThCell
  CLDReg = CellInfo%ThCell%CldReg; FuelReg = CellInfo%ThCell%FuelReg
  CoolReg = CellInfo%ThCell%CoolReg
  i1 = FuelReg(1); i2 = FuelReg(2); 
  FxrTag(i1:i2) = 1
  i1 = CLDReg(1); i2 = CLDReg(2)
  FxrTag(i1:i2) = 3
  i1 = CoolReg(1); i2 = CoolReg(2)
  FxrTag(i1:i2) = 2
  DO i = 1, nFxr
    IF(FxrTag(i) .EQ. 0) FxrTag(i) = 4
  ENDDO
  NULLIFY(ThCell)
ELSE
  DO i = 1, nFxr
    FxrTag(i) = 4
    IF(CellFxr(i)%lh2o) FxrTag(i) = 3
  ENDDO
ENDIF

END SUBROUTINE

SUBROUTINE GetRegNum(RegNum, FxrTag, CellFxr, nFxr)
IMPLICIT NONE
TYPE(FxrInfo_Type), INTENT(IN) :: CellFxr(nfxr)
INTEGER, INTENT(IN) ::nfxr
INTEGER, INTENT(IN) :: FxrTag(nFxr)
REAL, INTENT(OUT) :: RegNum(4)
INTEGER :: i1, i2, i, j
INTEGER :: niso

RegNum = 0
DO j = 1, nFxr
  i = FxrTag(j)
  niso = CellFxr(j)%niso
  RegNum(i) = RegNum(i) + sum(CellFxr(j)%pnum(1 : niso))
ENDDO
END SUBROUTINE

SUBROUTINE GetCellTemp(TSet, FxrTag, CellFxr, nfxr)
IMPLICIT NONE
TYPE(FxrInfo_Type), INTENT(IN) :: CellFxr(nfxr)
INTEGER, INTENT(IN) ::nfxr
INTEGER, INTENT(IN) :: FxrTag(nFxr)
REAL, INTENT(OUT) :: Tset(4)

INTEGER :: i1, i2, i, j
REAL :: TempSum(4), AreaSum(4)

AreaSum =0; TempSum = 0
DO j = 1, nFxr
  i = FxrTag(j)
  AreaSum(i) = AreaSum(i) + CellFxr(j)%Area
  TempSum(i) = TempSum(i) + CellFxr(j)%Area * CellFxr(j)%Temp
ENDDO

Tset = 0
DO j = 1, 4
  IF(AreaSum(j) .GT. 1.E-5) THEN
    TSet(j) = TempSum(j) / AreaSum(j)
  ENDIF
ENDDO

END SUBROUTINE

!SUBROUTINE PinXsFtn0(CellInfo, PinXS, PinXSFtnData, ng, ixy, iz)
!USE PARAM
!USE TYPEDEF,        ONLY : Cell_TYPE
!USE BasicOperation, ONLY : AD_VA         ,CP_CA          ,Multi_CA             &
!                          ,SUB_VA
!IMPLICIT NONE
!TYPE(Cell_Type) :: CellInfo
!TYPE(PinXs_Type) :: PinXS, PinXSFtnData(3)
!INTEGER :: ng, ixy ,iz
!
!REAL :: RefT0, RefFuelT0, RefModT0
!REAL :: RefT1, RefFuelT1, RefModT1
!REAL :: T, FuelT, ModT
!REAL :: RDFT, RDMT, DFT, DMT, FRAC
!REAL :: DelXs(ng, 5),DelXss(ng), DelXs0(ng, 5), DelXss0(ng,ng)
!
!INTEGER :: ig, ib, ie
!RefT0 = PinXSFtnData(1)%PinTemp
!RefFuelT0 = PinXSFtnData(1)%FuelTemp; RefModT0 = PinXSFtnData(1)%ModTemp
!RefT1 = PinXSFtnData(3)%PinTemp
!RefFuelT1 = PinXSFtnData(2)%FuelTemp; RefModT1 = PinXSFtnData(3)%ModTemp
!T = PinXs%PinTemp; FuelT = PinXs%FuelTemp; ModT =PinXs%ModTemp
!
!IF(.NOT. CellInfo%lFuel) THEN
!  RefModT0 = RefT0
!  RefModT1 = RefT1
!  ModT = T
!ENDIF
!
!DO ig = 1, ng
!  PinXs%Xss(ig)%from(ig) = PinXs%Xss(ig)%WithInGroupScat
!ENDDO
!CALL CP_CA(DelXs0, 0._8, ng, 5)
!CALL CP_CA(DelXss0, 0._8, ng, ng)
!
!DMT = (ModT - RefModT0)
!RDMT = (RefModT1 - RefModT0)
!
!DelXs(1:ng, 1) = PinXsFtnData(3)%xst - PinXsFtnData(1)%xst
!DelXs(1:ng, 2) = PinXsFtnData(3)%xstr - PinXsFtnData(1)%xstr
!DelXs(1:ng, 3) = PinXsFtnData(3)%xsnf - PinXsFtnData(1)%xsnf
!DelXs(1:ng, 4) = PinXsFtnData(3)%xskf - PinXsFtnData(1)%xskf
!FRAC = 0
!IF(abs(RDMT) .GT. epsm5) THEN
!  FRAC = DMT / RDMT
!ENDIF
!
!CALL MULTI_CA(FRAC, DelXS(1:ng, 1:4), ng, 4)
!
!!WRITE(64, '(2I8, 3F15.7, 200e20.5)') iz, ixy, ModT,RefModT0, RefModT1, FRAC, DelXS(1:ng, 1)
!!DelXS(1:ng, 1:2) = 0; FRAC = 0;
!CALL AD_VA(DelXs0(1:ng, 1:4), DelXs0(1:ng, 1:4), DelXS(1:ng, 1:4), ng, 4)
!
!DO ig = 1, ng
!  ib = PinXs%Xss(ig)%ib; ie = PinXs%Xss(ig)%ie
!  DelXss(ib:ie) = PinXsFtnData(3)%Xss(ig)%from(ib:ie) - PinXsFtnData(1)%Xss(ig)%from(ib:ie)
!  DelXss(ig) = PinXsFtnData(3)%Xss(ig)%WithInGroupScat - PinXsFtnData(1)%Xss(ig)%WithInGroupScat
!  CALL MULTI_CA(Frac, DelXss(ib:ie), ie-ib+1)
!  CALL AD_VA(DelXss0(ib:ie, ig), DelXss0(ib:ie, ig), DelXSs(ib:ie), ie-ib+1)
!ENDDO
!
!IF(CellInfo%lFuel) THEN
!  DFT = (SQRT(FuelT) - SQRT(RefFuelT0))
!  RDFT = (SQRT(RefFuelT1) - SQRT(RefFuelT0))
!  !FRAC = 0
!  IF(abs(RDFT) .GT. epsm5) THEN
!    FRAC = DFT / RDFT
!  ELSE
!    PRINT *, RefFuelT1, RefFuelT0
!    PAUSE
!  ENDIF
!  DelXs(1:ng, 1) = PinXsFtnData(2)%xst - PinXsFtnData(1)%xst
!  DelXs(1:ng, 2) = PinXsFtnData(2)%xstr - PinXsFtnData(1)%xstr
!  DelXs(1:ng, 3) = PinXsFtnData(2)%xsnf - PinXsFtnData(1)%xsnf
!  DelXs(1:ng, 4) = PinXsFtnData(2)%xskf - PinXsFtnData(1)%xskf
!  CALL MULTI_CA(FRAC, DelXS(1:ng, 1:4), ng, 4)
!  CALL AD_VA(DelXs0(1:ng, 1:4), DelXs0(1:ng, 1:4), DelXS(1:ng, 1:4), ng, 4)
!  DO ig = 1, ng
!    ib = PinXs%Xss(ig)%ib; ie = PinXs%Xss(ig)%ie
!    DelXss(ib:ie) = PinXsFtnData(2)%Xss(ig)%from(ib:ie) - PinXsFtnData(1)%Xss(ig)%from(ib:ie)
!    DelXss(ig) = PinXsFtnData(2)%Xss(ig)%WithInGroupScat - PinXsFtnData(1)%Xss(ig)%WithInGroupScat
!    CALL MULTI_CA(Frac, DelXss(ib:ie), ie-ib+1)
!    CALL AD_VA(DelXss0(ib:ie, ig), DelXss0(ib:ie, ig), DelXSs(ib:ie), ie-ib+1)
!  ENDDO
!ENDIF
!
!!UPDATE XS
!CALL AD_VA(PinXS%xst(1:ng), PinXS%xst(1:ng), DelXs0(1:ng, 1), ng)
!CALL AD_VA(PinXS%xstr(1:ng), PinXS%xstr(1:ng), DelXs0(1:ng, 2), ng)
!IF(CellInfo%lFuel) THEN
!  CALL AD_VA(PinXS%xsnf(1:ng), PinXS%xsnf(1:ng), DelXs0(1:ng, 3), ng)
!  CALL AD_VA(PinXS%xskf(1:ng), PinXS%xskf(1:ng), DelXs0(1:ng, 4), ng)
!ENDIF
!
!DO ig = 1, ng
!  ib = PinXs%Xss(ig)%ib; ie = PinXs%Xss(ig)%ie
!  CALL AD_VA(PinXs%Xss(ig)%from(ib:ie), PinXs%Xss(ig)%from(ib:ie), DelXss0(ib:ie, ig), ie - ib +1)
!  PinXS%Xss(ig)%WithInGroupScat = PinXs%Xss(ig)%from(ig)
!  PinXs%Xss(ig)%from(ig) = 0
!ENDDO
!!SET RMV XS
!DO ig = 1, ng
!  PinXS%xsr(ig) = PinXS%xstr(ig) - PinXs%Xss(ig)%WithInGroupScat
!  PinXS%XSD(ig) = 1._8/(3._8 * PinXS%xstr(ig))
!ENDDO
!
!END SUBROUTINE

END MODULE