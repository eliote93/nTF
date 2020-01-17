#include <defines.h>
SUBROUTINE PsiUpdate(Core, Fxr, phis, psi, myzb, myze, ng, lxslib, GroupInfo)
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,       Fxrinfo_type,       Cell_Type,     pin_Type, &
                         GroupInfo_Type,      XsMac_Type
USE BenchXs,       ONLY : xsnfBen
USE MacXsLib_Mod, ONLY : MacXsNf, IsoMacXsnf
USE BasicOperation, ONLY : CP_CA, MULTI_VA
IMPLICIT NONE
TYPE(coreinfo_type) :: CORE
TYPE(Fxrinfo_type),POINTER :: Fxr(:, :)
TYPE(GroupInfo_Type) :: GroupInfo
REAL, POINTER :: phis(:, :, :)
REAL, POINTER :: Psi(:, :)
INTEGER :: myzb, myze, ng
LOGICAL :: lXsLib


TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
TYPE(Fxrinfo_type),POINTER :: myFxr
TYPE(XsMac_Type), SAVE :: XsMac

INTEGER :: nxy, nCoreFsr, nCoreFxr, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr
INTEGER :: ipin, icel, ifsrlocal, ifsr, ifxr, iz, itype, ig
INTEGER :: iResoGrpBeg, iResoGrpEnd, norg
INTEGER :: i, j, k, iso

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
  CALL CP_CA(Psi(:, iz), zero, nCoreFsr)
  DO ipin = 1, nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
    icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr  
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)    
      myFxr => Fxr(ifxr, iz)
      IF (lXsLib) THEN
        CALL MacXsNf(XsMac, myFxr, 1, ng, ng, 1._8, FALSE, TRUE)
        xsmacnf => XsMac%XsMacNf
        IF(myFxr%lres) THEN
          do ig = iResoGrpBeg, iResoGrpEnd
            XsMacNf(ig) = XsMacNf(ig) * myFxr%fresoF(ig)  
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
          Psi(ifsr, iz) = Psi(ifsr, iz) + xsmacnf(ig) * phis(ifsr, iz, ig)
        ENDDO
        CONTINUE
        !src(ifsr) = reigv * chi(ig) * psic(ifsr, iz)
      ENDDO !Fsr Sweep    
    ENDDO
  ENDDO
ENDDO

IF(.NOT. lxsLib) Deallocate(xsmacnf)
IF(lXsLib) NULLIFY(XsMacNf)
NULLIFY(Pin, CellInfo)
END SUBROUTINE

SUBROUTINE PowerUpdate(Core, Fxr, phis, power, myzb, myze, ng, lxslib, GroupInfo, PE)
USE PARAM
USE TYPEDEF,      ONLY : coreinfo_type,       Fxrinfo_type,       Cell_Type,     pin_Type, &
                         GroupInfo_Type,      XsMac_Type,         PE_TYPE
USE BenchXs,       ONLY : xskfBen
USE MacXsLib_Mod, ONLY : MacXskf
USE BasicOperation, ONLY : CP_CA, MULTI_VA
#ifdef MPI_ENV
USE MPIComm_Mod, ONLY : BCAST
#endif
IMPLICIT NONE
TYPE(coreinfo_type) :: CORE
TYPE(Fxrinfo_type),POINTER :: Fxr(:, :)
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_Type) :: PE
REAL, POINTER :: phis(:, :, :)
REAL, POINTER :: Power(:, :)
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

REAL, POINTER :: xsmackf(:)

Pin => Core%Pin
CellInfo => Core%CellInfo; nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr; nxy = Core%nxy
IF(lxslib) THEN
  iResoGrpBeg = GroupInfo%nofg + 1
  iResoGrpEnd = GroupInfo%nofg + GroupInfo%norg
  norg = GroupInfo%norg
ENDIF

IF(.NOT. lxsLib) ALLOCATE(xsmackf(ng))
CALL CP_CA(Power(1:nCoreFsr, 1:Core%nz), zero, nCoreFsr, Core%nz)
DO iz = myzb, myze
  
  DO ipin = 1, nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
    icel = Pin(ipin)%Cell(iz); nlocalFxr = CellInfo(icel)%nFxr  
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      nFsrInFxr = CellInfo(icel)%nFsrInFxr(j)    
      myFxr => Fxr(ifxr, iz)
      IF(lXsLib) Then
        CALL MacXsKf(XsMac, myFxr, 1, ng, ng, 1._8, FALSE)
        xsmacKf => XsMac%XsMacKf
        IF(myFxr%lres) THEN
          do ig = iResoGrpBeg, iResoGrpEnd
            XsMackf(ig) = XsMackf(ig) * myFxr%fresoF(ig)  
          enddo
        ENDIF
      ELSE
        ifsrlocal = CellInfo(icel)%MapFxr2FsrIdx(1,j)
        !itype = CellInfo(icel)%iReg(ifsrlocal)      
        itype = myFxr%imix
        CALL xskfben(itype, 1, ng, xsmackf)
        !CHI(ig:ig) = GetChiBen(itype, ig, ig)
      ENDIF
      
      DO i = 1, nFsrInFxr
        ifsr = FsrIdxSt + Cellinfo(icel)%MapFxr2FsrIdx(i, j) - 1  !Global FSR Index
        DO ig = 1, ng
          Power(ifsr, iz) = power(ifsr, iz) + xsmackf(ig) * phis(ifsr, iz, ig)
        ENDDO
        CONTINUE
        !src(ifsr) = reigv * chi(ig) * psic(ifsr, iz)
      ENDDO !Fsr Sweep    
    ENDDO
  ENDDO
ENDDO
#ifdef MPI_ENV
DO iz = 1, Core%nz
  CALL BCAST(Power(1:nCoreFsr, iz), nCoreFsr, PE%MPI_RTMASTER_COMM, PE%AxDomList(iz))
ENDDO
#endif
IF(.NOT. lxsLib) Deallocate(xsmackf)
IF(lXsLib) NULLIFY(XsMackf)
NULLIFY(Pin, CellInfo)
END SUBROUTINE


SUBROUTINE UpdateEigv(Core, psi, psid, eigv, peigv, myzb, myze, PE)
USE PARAM
USE TYPEDEF, ONLY : coreinfo_type, Cell_Type, pin_Type, GroupInfo_Type, PE_TYPE
USE BenchXs, ONLY : xsnfBen
USE BasicOperation, ONLY : CP_CA
#ifdef MPI_ENV
USE MPICOMM_MOD, ONLY : REDUCEnBCAST, REDUCE
#endif
IMPLICIT NONE
TYPE(coreinfo_type) :: CORE
REAL, POINTER :: Psi(:, :), psid(:, :)
INTEGER :: myzb, myze, ng
REAL :: eigv, peigv
TYPE(PE_TYPE), OPTIONAL :: PE

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)

INTEGER :: nxy, nCoreFsr, nCoreFxr, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr
INTEGER :: ipin, icel, ifsrlocal, ifsr, ifxr, iz, itype, ig
INTEGER :: i, j, k
REAL :: psipsi, psipsid, temp
REAL, POINTER :: hz(:)
LOGICAL :: RTMASTER

Pin => Core%Pin
CellInfo => Core%CellInfo
hz => Core%hz
nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr
nxy = Core%nxy
psipsi = 0
psipsid = 0
peigv = eigv
RtMaster = .TRUE.
IF(PRESENT(PE)) RtMaster = PE%RTMaster
IF(RTMASTER) THEN
  DO iz = myzb, myze
    DO ipin = 1, nxy
      FsrIdxSt = Pin(ipin)%FsrIdxSt; icel = Pin(ipin)%Cell(iz);
      DO j = 1, CellInfo(icel)%nFsr
        IF (CellInfo(icel)%vol(j) < 1E-10) CYCLE
        
        ifsr = FsrIdxSt + j - 1
        psipsi = psipsi + psi(ifsr, iz) * psi(ifsr, iz) * (hz(iz) * CellInfo(icel)%vol(j))**2
        psipsid = psipsid + psi(ifsr, iz) * psid(ifsr, iz) * (hz(iz) * CellInfo(icel)%vol(j))**2
      ENDDO  
    ENDDO
  ENDDO
ENDIF
#ifdef MPI_ENV
IF(PRESENT(PE)) THEN
  CALL REDUCEnBCAST(psipsi, PE%MPI_RTMASTER_COMM, PE%MPI_RT_COMM, PE%RTMaster, TRUE, TRUE)
  CALL REDUCEnBCAST(psipsid, PE%MPI_RTMASTER_COMM, PE%MPI_RT_COMM, PE%RTMaster, TRUE, TRUE)
ENDIF
#endif
eigv = eigv*psipsi/psipsid

NULLIFY(Pin)
NULLIFY(CellInfo)
NULLIFY(hz)
END SUBROUTINE

SUBROUTINE CellPsiUpdate(CORE, psi, psic, myzb, myze)
USE PARAM
USE TYPEDEF, ONLY : coreinfo_type, Fxrinfo_type, Cell_Type, pin_Type
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
TYPE(coreinfo_type) :: CORE
REAL, POINTER :: Psi(:, :)
REAL, POINTER :: psiC(:, :)
INTEGER :: myzb, myze

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)

INTEGER :: nxy, nCoreFsr, nCoreFxr, nLocalFsr
INTEGER :: l, i, j, k, iz
INTEGER :: FsrIdxSt, icel, ireg
REAL, POINTER :: hz(:)
Pin => Core%Pin; CellInfo => Core%CellInfo
hz => Core%hz
nCoreFsr = Core%nCoreFsr; nCoreFxr = Core%nCoreFxr
nxy = Core%nxy
DO iz = myzb, myze
  DO l = 1, nxy
    FsrIdxSt = Pin(l)%FsrIdxSt; icel = Pin(l)%Cell(iz);
    nLocalFsr = CellInfo(icel)%nFsr
    psic(l, iz) = zero
    DO j = 1, nLocalFsr
      ireg = FsrIdxSt + j - 1
      
      IF (CellInfo(icel)%vol(j) < 1E-10) CYCLE
      
      psic(l, iz) =  psic(l, iz) + CellInfo(icel)%vol(j) * psi(ireg, iz)
    ENDDO  
    psic(l, iz) = psic(l, iz)*hz(iz)
  ENDDO
ENDDO
NULLIFY(Pin)
NULLIFY(CellInfo)
NULLIFY(hz)
END SUBROUTINE

Function PsiErr(Core, psi, psid, myzb, myze, PE)
USE PARAM
USE TYPEDEF, ONLY : coreinfo_type, Cell_Type, pin_Type, PE_TYPE
USE BenchXs, ONLY : xsnfBen
USE BasicOperation, ONLY : CP_CA
#ifdef MPI_ENV
USE MPICOMM_MOD, ONLY : REDUCE
#endif
IMPLICIT NONE
TYPE(coreinfo_type) :: CORE
REAL :: PsiErr
REAL, POINTER :: Psi(:, :), psid(:, :)
INTEGER :: myzb, myze
TYPE(PE_TYPE) :: PE

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)

INTEGER :: nxy, nCoreFsr, nCoreFxr, FsrIdxSt, FxrIdxSt, nlocalFxr, nFsrInFxr
INTEGER :: ipin, icel, ifsrlocal, ifsr, ifxr, iz, itype
INTEGER :: i, j, k
REAL :: psipsi, ErrSqSum, localErr, temp
REAL, POINTER :: hz(:)

Pin => Core%Pin
CellInfo => Core%CellInfo
hz => Core%hz
nCoreFsr = Core%nCoreFsr
nCoreFxr = Core%nCoreFxr
nxy = Core%nxy
ErrSqSum =0; psipsi = 0
DO iz = myzb, myze
  DO ipin = 1, nxy
    FsrIdxSt = Pin(ipin)%FsrIdxSt; icel = Pin(ipin)%Cell(iz);
    DO j = 1, CellInfo(icel)%nFsr
      IF (CellInfo(icel)%vol(j) < 1E-10) CYCLE
      
      ifsr = FsrIdxSt + j - 1
      psipsi = psipsi + psi(ifsr, iz) * psi(ifsr, iz) * (hz(iz) * CellInfo(icel)%vol(j))**2
      localErr = (psi(ifsr, iz) - psid(ifsr, iz)) ** 2
      ErrSqSum = ErrSqSum + localErr* (hz(iz) * CellInfo(icel)%vol(j))**2
    ENDDO  
  ENDDO
ENDDO
#ifdef MPI_ENV
CALL REDUCE(ErrSqSum, temp, PE%MPI_RTMASTER_COMM, .TRUE.)
ErrSqSum = temp
CALL REDUCE(psipsi, temp, PE%MPI_RTMASTER_COMM, .TRUE.)
psipsi = temp
#endif
PsiErr = SQRT(ErrSqSum/psipsi)
END FUNCTION


!SUBROUTINE ReducePhiAngIn(RayInfo, PhiAngIn1g, )
!
!END SUBROUTINE

FUNCTION MocResidual(Core, FmInfo, eigv, GroupInfo, ng, PE, nTracerCntl)
USE PARAM
USE TYPEDEF,  ONLY : CoreInfo_Type,     FmInfo_Type,     GroupInfo_Type,     &
                     Cell_Type,         PinInfo_Type,     Pin_Type,             &
                     FxrInfo_Type
USE MOC_MOD,  ONLY : tSrc,              xst1g,            AxSrc1g,              &
                     AxPxs1g,                                                   & 
                     SetRtSrc,          SetRtMacXs,       PseudoAbsorption,     &
                     AddBuckling,       AddConstSrc,                            &
                     !--- CNJ Edit : Node Majors
                     phisnm,            srcnm,            xstnm,                &
                     SetRtSrcNM,        SetRtMacXsNM,     PseudoAbsorptionNM,   &
                     AddBucklingNM
USE geom,     ONLY : nbd
USE cntl,     ONLY : nTracerCntl_Type
USE PE_MOD,   ONLY : PE_TYPE
USE BenchXs,  ONLY : GetXstrBen
USE DcplXsGen_Mod,  ONLY : DcplSetMocAxEff
USE BasicOperation, ONLY : CP_CA, CP_VA
#ifdef MPI_ENV
USE MPICOMM_MOD, ONLY : REDUCE
#endif
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
REAL :: eigv
TYPE(GroupInfo_Type) :: GroupInfo
INTEGER :: ng
TYPE(PE_TYPE) :: PE
TYPE(nTracerCntl_Type) :: nTracerCntl
REAL :: MocResidual

TYPE(FxrInfo_Type), POINTER :: FXR(:, :)
REAL, POINTER :: PHIS(:, :, :)
REAL, POINTER :: PSI(:, :)
REAL, POINTER :: MocJout(:, :, :, :, :)
REAL, POINTER :: AxSrc(:, :, :), AxPXS(:, :, :)
REAL, POINTER :: Res(:)

!POINTING Variables
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(PinInfo_Type), POINTER :: PinInfo(:)
TYPE(Cell_Type), POINTER :: Cell(:)
!
INTEGER :: nxy, nFsr, nlocalFxr, nlocalFsr, FsrIdxSt, FxrIdxSt, nFsrInFxr, myzb, myze
INTEGER :: ig, ireg, ipin, icel, ifxr, ixy, iz, ifsrlocal, itype
INTEGER :: i, j, k, l, m
LOGICAL :: lXsLib, l3dim
REAL :: SrcSum, LocalResidual, localsrc, vol, temp
REAL :: Leakage(100), Collision(100), Source(100)   !--- CNJ Edit : Node Majors

LOGICAL :: lNegFix 

l3dim = nTracerCntl%l3Dim
lXsLib = nTracerCntl%lXsLib

Fxr => FmInfo%Fxr
Phis => FmInfo%Phis
Psi => FmInfo%Psi
MocJout => FmInfo%RadJout !ninout/ 1:in 2:out 3:surfphi 16/02/11 BYS edit
AxSrc => FmInfo%AxSrc; AxPXS => FmInfo%AxPXS

Pin => Core%Pin
PinInfo => Core%Pininfo
Cell => Core%CellInfo
lXsLib = nTracerCntl%lXsLib; l3dim = nTracerCntl%l3dim
myzb = PE%myzb; myze = PE%myze
nFsr = Core%nCoreFsr; nxy = Core%nxy

MocResidual = 0; SrcSum = 0
!AxSrc => FmInfo%AxSrc; AxPXS => FmInfo%AxPXS
CALL CP_CA(AxPxs1g(1:nxy), 0._8, nxy)
CALL CP_CA(AxSrc1g(1:nxy), 0._8, nxy)
ALLOCATE(Res(1:nxy))
#ifndef ResMG
DO iz = myzb, myze
  CALL CP_CA(Res(1:nxy), 0._8, nxy)
  IF(.NOT. core%lfuelplane(iz) .AND. nTRACErCntl%lAxRefFDM) CYCLE
  !--- CNJ Edit : Node Majors
  IF (nTracerCntl%lNodeMajor) THEN
    DO ig = 1, ng
      phisnm(ig, :) = phis(:, iz, ig)
    ENDDO
    CALL SetRtMacXsNM(Core, Fxr(:, iz), xstnm, iz, ng, lxslib, nTracerCntl%lTrCorrection,   &
                      nTracerCntl%lRST, FALSE, FALSE, PE)
#ifdef LkgSplit
    CALL PseudoAbsorptionNM(Core, Fxr(:, iz), AxPXS, xstnm, iz, ng, GroupInfo, l3dim)
#endif
#ifdef Buckling
    IF (nTracerCntl%lBsq) CALL AddBucklingNM(Core, Fxr, xstnm, nTracerCntl%Bsq, iz, ng, lxslib, nTracerCntl%lRST)
#endif
    CALL SetRtSrcNM(Core, Fxr(:, iz), srcnm, phisnm, psi, AxSrc, xstnm, eigv, iz,           &
                    1, ng, ng, GroupInfo, l3dim, lxslib, nTracerCntl%lscat1, lNegFix, PE)
    DO ipin = 1, nxy
      FsrIdxSt = Pin(ipin)%FsrIdxSt
      icel = Pin(ipin)%Cell(iz)
      nlocalFsr = Cell(icel)%nFsr
      Collision = 0; Source = 0; Leakage = 0
      DO ig = 1, ng
        DO i = 1, nbd
          Leakage(ig) = Leakage(ig) + MocJout(2, i, ipin, iz, ig) - MocJout(1, i, ipin, iz, ig)
        ENDDO
      ENDDO
      DO i = 1, nlocalFsr
        ireg = FsrIdxSt + i - 1
        vol = Cell(icel)%vol(i)
        
        IF (vol < 1E-10) CYCLE
        
        DO ig = 1, ng
          Collision(ig) = Collision(ig) + vol * phisnm(ig, ireg) * xstnm(ig, ireg)
          Source(ig) = Source(ig) + vol * srcnm(ig, ireg) * xstnm(ig, ireg)
        ENDDO
      ENDDO
      DO ig = 1, ng
        LocalResidual = Source(ig) - Collision(ig) - Leakage(ig)
        SrcSum = SrcSum + Source(ig) ** 2
        MocResidual = MocResidual + LocalResidual ** 2
      ENDDO
    ENDDO
  ELSE
    DO ig = 1, ng
      !CALL CP_VA(AxPx1g)
      lNegFix = nTracerCntl%lDcplCal .AND. .NOT. Core%lFuelPlane(iz)
      IF(nTracerCntl%l3dim .AND. .NOT. nTracerCntl%lDcplCal) THEN
        CALL CP_VA(AxPxs1g(1:nxy), AxPxs(1:nxy, iz, ig), nxy)
        CALL CP_VA(AxSrc1g(1:nxy), AxSrc(1:nxy, iz, ig), nxy)
      ELSEIF(nTracerCntl%l3dim .or. nTracerCntl%lDcplCal) THEN
        CALL DcplSetMocAxEff(Core, AxPxs(:, iz, ig), phis(:, iz, ig), AxSrc1g, AxPxs1g, iz)
      ENDIF
      CALL CP_CA(xst1g(1:nFsr), 1._8, nFsr)
      IF(nTracerCntl%lDcplCal) THEN
        CALL SetRtSrc(Core, Fxr(:, iz), tsrc, phis, psi, axSrc1g, xst1g,                     &
                      eigv, iz, ig, ng, GroupInfo, TRUE, lXslib, FALSE, lNegFix, PE)    
      ELSE
        !CALL SetRtSrc(Core, Fxr(:, iz), tsrc(:), phis, psi, axSrc1g(:), xst1g(:), &
                      !eigv, iz, ig, ng, GroupInfo, l3dim, lXslib, FALSE, lNegFix, PE)
        CALL SetRtSrc(Core, Fxr(:, iz), tsrc, phis, psi, axSrc1g, xst1g, &
                      eigv, iz, ig, ng, GroupInfo, l3dim, lXslib, nTracerCntl%lscat1, lNegFix, PE)
      ENDIF
      IF(nTracerCntl%lDcplCal .AND. .NOT. Core%lFuelPlane(iz)) THEN
        CALL AddConstSrc(Core, Fxr(:, iz), tsrc, xst1g, nTracerCntl%ConstSrc, iz, ig, ng)
      ENDIF
      DO ipin = 1, nxy
        FsrIdxSt = Pin(ipin)%FsrIdxSt
        icel = Pin(ipin)%Cell(iz); nlocalFsr = Cell(icel)%nFsr
        DO j = 1, nlocalFsr
          IF (Cell(icel)%vol(j) < 1E-10) CYCLE
          
          ireg = FsrIdxSt + j - 1
          vol = Cell(icel)%vol(j)
          Res(ipin) = Res(ipin)  + vol * tsrc(ireg)
        ENDDO
      ENDDO
    ENDDO
    !Get total source term
    DO ipin = 1, nxy
      SrcSum = SrcSum + Res(ipin) * Res(ipin)
    ENDDO
    DO ig = 1, ng
      lNegFix = nTracerCntl%lDcplCal .AND. .NOT. Core%lFuelPlane(iz)
      IF(nTracerCntl%l3dim .AND. .NOT. nTracerCntl%lDcplCal) THEN
        CALL CP_VA(AxPxs1g(1:nxy), AxPxs(1:nxy, iz, ig), nxy)
        CALL CP_VA(AxSrc1g(1:nxy), AxSrc(1:nxy, iz, ig), nxy)
      ELSEIF(nTracerCntl%l3dim .or. nTracerCntl%lDcplCal) THEN
        CALL DcplSetMocAxEff(Core, AxPxs(:, iz, ig), phis(:, iz, ig), AxSrc1g, AxPxs1g, iz)
      ENDIF
    
    CALL SetRtMacXs(Core, Fxr(:, iz), xst1g, iz, ig, ng, lxslib, nTracerCntl%lTrCorrection, nTracerCntl%lRST, FALSE, FALSE, PE)
      IF(nTracerCntl%l3dim .or. nTracerCntl%lDcplCal) THEN
#ifdef LkgSplit  
        CALL PseudoAbsorption(Core, Fxr(:, iz), tsrc, phis(:, iz, ig),     &
                              AxPXS1g(:), xst1g, iz, ig, ng, GroupInfo, true)  
#endif
      ENDIF
#ifdef Buckling
    IF(nTracerCntl%lBsq) CALL AddBuckling(Core, Fxr, xst1g, nTracerCntl%Bsq, iz, ig, ng, lxslib, nTracerCntl%lRST)
#endif 
      DO ipin = 1, nxy
        FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
        icel = Pin(ipin)%Cell(iz)
        nlocalFsr = Cell(icel)%nFsr
        LocalResidual = 0; localsrc = 0
    !    !Current
        DO j = 1, nbd
          LocalResidual = LocalResidual + (MocJout(2, j, ipin, iz, ig) - MocJout(1, j, ipin, iz, ig))
        ENDDO
        
        DO j = 1, nlocalFsr
          IF (Cell(icel)%vol(j) < 1E-10) CYCLE
          
          ireg = FsrIdxSt + j - 1
          vol = Cell(icel)%vol(j)
          LocalResidual = LocalResidual + vol * phis(ireg, iz, ig) * xst1g(ireg)
        ENDDO
        Res(ipin) = Res(ipin) - LocalResidual
      ENDDO  !Pin Sweep
    ENDDO  !Group Sweep
    DO ipin = 1, nxy
       MocResidual = MocResidual + Res(ipin) * Res(ipin)
      !write(*,*) Res(ipin)*Res(ipin)
    ENDDO
  ENDIF
ENDDO
DEALLOCATE(Res)
#else
DO ig = 1, ng
  DO iz = myzb, myze
    IF(.NOT. core%lfuelplane(iz) .AND. nTRACErCntl%lAxRefFDM) CYCLE
    !CALL CP_VA(AxPx1g)
    lNegFix = nTracerCntl%lDcplCal .AND. .NOT. Core%lFuelPlane(iz)
    IF(nTracerCntl%l3dim .AND. .NOT. nTracerCntl%lDcplCal) THEN
      CALL CP_VA(AxPxs1g(1:nxy), AxPxs(1:nxy, iz, ig), nxy)
      CALL CP_VA(AxSrc1g(1:nxy), AxSrc(1:nxy, iz, ig), nxy)
    ELSEIF(nTracerCntl%l3dim .or. nTracerCntl%lDcplCal) THEN
      CALL DcplSetMocAxEff(Core, AxPxs(:, iz, ig), phis(:, iz, ig), AxSrc1g, AxPxs1g, iz)
      !CALL CP_VA(AxPxs1g(1:nxy), AxPxs(1:nxy, iz, ig), nxy)
      !CALL CP_VA(AxSrc1g(1:nxy), AxSrc(1:nxy, iz, ig), nxy)
    ENDIF
    
    CALL SetRtMacXs(Core, Fxr(:, iz), xst1g, iz, ig, ng, lxslib, nTracerCntl%lTrCorrection, nTracerCntl%lRST, FALSE, FALSE, PE)
    IF(nTracerCntl%l3dim .or. nTracerCntl%lDcplCal) THEN
#ifdef LkgSplit  
      CALL PseudoAbsorption(Core, Fxr(:, iz), tsrc, phis(:, iz, ig),     &
                            AxPXS1g(:), xst1g, iz, ig, ng, GroupInfo, true)  
#endif
    ENDIF
#ifdef Buckling
   IF(nTracerCntl%lBsq) CALL AddBuckling(Core, Fxr, xst1g, nTracerCntl%Bsq, iz, ig, ng, lxslib, nTracerCntl%lRST)
#endif  
    IF(nTracerCntl%lDcplCal) THEN
      CALL SetRtSrc(Core, Fxr(:, iz), tsrc, phis, psi, axSrc1g, xst1g,                &
                    eigv, iz, ig, ng, GroupInfo, TRUE, lXslib, FALSE, lNegFix, PE)    
    ELSE
      CALL SetRtSrc(Core, Fxr(:, iz), tsrc, phis, psi, axSrc1g, xst1g,                &
                    eigv, iz, ig, ng, GroupInfo, l3dim, lXslib, FALSE, lNegFix, PE)
    ENDIF
    IF(nTracerCntl%lDcplCal .AND. .NOT. Core%lFuelPlane(iz)) THEN
      CALL AddConstSrc(Core, Fxr(:, iz), tsrc, xst1g, nTracerCntl%ConstSrc, iz, ig, ng)
    ENDIF
    DO ipin = 1, nxy
      FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
      icel = Pin(ipin)%Cell(iz)
      nlocalFsr = Cell(icel)%nFsr
      LocalResidual = 0; localsrc = 0
      !Current
      DO j = 1, nbd
        LocalResidual = LocalResidual + (MocJout(2, j, ipin, iz, ig) - MocJout(1, j, ipin, iz, ig))
      ENDDO
      
      DO j = 1, nlocalFsr
        IF (Cell(icel)%vol(j) < 1E-10) CYCLE
        
        ireg = FsrIdxSt + j - 1
        vol = Cell(icel)%vol(j)
        
        IF (vol < 1E-10) CYCLE
        
        LocalResidual = LocalResidual + vol * phis(ireg, iz, ig) * xst1g(ireg)
        localsrc = localsrc + tsrc(ireg) * vol * xst1g(ireg)
      ENDDO
      
      LocalResidual = localsrc - LocalResidual
      SrcSum = SrcSum + localsrc * localsrc
      MocResidual = MocResidual + LocalResidual * LocalResidual
    ENDDO
  ENDDO
ENDDO
#endif
#ifdef MPI_ENV
CALL REDUCE(MocResidual, temp, PE%MPI_RTMASTER_COMM, .TRUE.)
MocResidual = temp
CALL REDUCE(SrcSum, temp, PE%MPI_RTMASTER_COMM, .TRUE.)
SrcSum = temp
#endif
!!--- BYS edit / 160223 divergence problem in 3D
!IF(MocResidual .LT. 0 )THEN 
!    write(*,*) 'MocResidual .LT. ZERO !!', MocResidual
!ELSE
!    write(*,*) 'MocResidual .GT. ZERO -- good' , MocResidual   
!ENDIF
!IF(SrcSum .LT. 0 )THEN 
!    write(*,*) 'SrcSum      .LT. ZERO !!', SrcSum
!ELSE
!    write(*,*) 'SrcSum      .GT. ZERO -- good', SrcSum
!ENDIF
!!--- BYS edit END / 160223 divergence problem in 3D
MocResidual = MocResidual / SrcSum
MocResidual = SQRT(MocResidual)


NULLIFY(Pin)
NULLIFY(PinInfo)
NULLIFY(Cell)
END FUNCTION



FUNCTION FxrAvgPhi(Core, Fxr, Phis, ipin, iLocalfxr, iz, ng, PE)
USE PARAM
USE TYPEDEF,     ONLY : CoreInfo_Type,    PE_Type,     FxrInfo_Type,               &
                        Cell_Type,    Pin_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(PE_Type) :: PE
REAL, POINTER :: phis(:, :, :)
INTEGER :: iLocalfxr, ipin, iz, ng
REAL :: FxrAvgPhi(ng)

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)

REAL :: area, areasum
INTEGER :: FsrIdxSt, FxrIdxSt, nFsrInFxr
INTEGER :: ifxr, ifsr, icell
INTEGER :: i, j

CellInfo => Core%CellInfo
Pin => Core%Pin

FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt

ifxr = FxrIdxSt + iLocalFxr - 1; icell = Pin(ipin)%Cell(iz)
nFsrInFxr = Fxr(ifxr, iz)%nFsrInFxr

FxrAvgPhi = 0; areasum = 0
DO i = 1, nFsrInFxr
  j = CellInfo(icell)%MapFxr2FsrIdx(i, iLocalFxr)
  iFsr = FsrIdxSt + j - 1
  Area = CellInfo(icell)%vol(j); AreaSum = AreaSum + Area
  FxrAvgPhi = FxrAvgPhi + Area * Phis(ifsr, iz, 1:ng)
ENDDO
FxrAvgPhi = FxrAvgPhi / AreaSum

NULLIFY(Pin, CellInfo)
END FUNCTION

SUBROUTINE FluxUnderRelaxation(Core, Phis1g, Phis, w, iz, ig, PE)
USE PARAM
USE TYPEDEF,  ONLY : CoreInfo_Type,  PE_TYPE,  Pin_Type, Cell_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
REAL, POINTER :: Phis1g(:)
REAL, POINTER :: Phis(:, :, :)
REAL :: w
INTEGER :: iz, ig
TYPE(PE_TYPE) :: PE

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: CellInfo(:)
INTEGER :: nxy
INTEGER :: ixy, ireg, icel, j
INTEGER :: FsrIdxSt
REAL :: wbar
Pin => Core%Pin
CellInfo => Core%CellInfo

wbar = 1 -w
nxy = Core%nxy
DO ixy = 1, nxy
  FsrIdxSt = Pin(ixy)%FsrIdxSt
  icel = Pin(ixy)%Cell(iz)
  DO j = 1, CellInfo(icel)%nFsr
    ireg = FsrIdxSt + j - 1
    Phis(ireg, iz, ig) = w * Phis1g(ireg) + wbar * Phis(ireg, iz, ig)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE FluxInUnderRelaxation(Core, PhiAngIn1g, PhiAngIn, w, n1, n2, iz, ig, PE)
USE PARAM
USE TYPEDEF, ONLY : CoreInfo_Type, PE_TYPE
IMPLICIT NONE

REAL, POINTER :: PhiAngIn1g(:, :)
REAL, POINTER :: PhiAngIn(:, :, :, :)
INTEGER :: n1, n2, iz, ig
REAL :: w
TYPE(PE_TYPE) :: PE

TYPE(CoreInfo_Type) :: Core
INTEGER :: i, j
REAL :: wbar

wbar = 1-w

FORALL(i=1:n2,j=1:n1)
  PhiAngIn(j, i, iz, ig) = w*PhiANgIn1g(j, i) + wbar *  PhiAngIn(j, i, iz, ig) 
END FORALL
END SUBROUTINE

SUBROUTINE CurrentUnderRelaxation(Core, Jout1g, Jout, w, iz, ig, PE)
USE PARAM
USE TYPEDEF,  ONLY : CoreInfo_Type,  PE_TYPE,  Pin_Type, Cell_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
REAL, POINTER :: Jout1g(:, :, :)
REAL, POINTER :: Jout(:, :, :, :, :)
REAL :: w
INTEGER :: iz, ig
TYPE(PE_TYPE) :: PE

INTEGER :: ixy, ibd
INTEGER :: nxy

REAL :: jnet, jnet0, wbar
nxy = Core%nxy

wbar = 1 - w

DO ixy = 1, nxy
  DO ibd = 1, 4
    jnet = Jout1g(2, ibd, ixy) - Jout1g(1, ibd, ixy)
    jnet0 = Jout(2, ibd, ixy, iz, ig) - Jout(1, ibd, ixy, iz, ig)
    jnet = w * jnet + wbar * jnet0
    Jout(2, ibd, ixy, iz, ig) = jnet
    Jout(1, ibd, ixy, iz, ig) = 0
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE MOCUnderRelaxationFactor(Core, FmInfo, CmInfo, GroupInfo, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,       FmInfo_Type,     CmInfo_Type, &
                           GroupInfo_Type,      PE_Type,                      &
                           PinXS_Type,          Cell_Type
USE Cntl,           ONLY : nTracerCntl_Type
USE FILES,          ONLY : IO8
#ifdef MPI_ENV
USE MPIComm_mod,    ONLY : REDUCE
#endif
IMPLICIT NONE

TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE

TYPE(PinXS_TYPE), POINTER :: PinXS(:, :)

INTEGER :: nxy, myzb, myze, ng
INTEGER :: ig, ig2, iz, ixy

REAL, POINTER :: C(:), ScatSum(:), w(:, :), wbuf(:, :)
REAL :: crit, f, wlocal
nxy = Core%nxy; myzb = PE%myzb; myze = PE%myze
ng = GroupInfo%ng
PinXS => CmInfo%PinXS

IF(.NOT. nTracerCntl%lMOCUR) THEN
  FORALL(ig=1:ng) FMInfo%w(ig) = 1.0
  RETURN
ENDIF 

IF(.NOT. nTracerCntl%l3Dim) THEN
  FORALL(ig=1:ng) FMInfo%w(ig) = 1.0
  IF(.NOT. nTracerCntl%lOptUR) THEN
    FORALL(ig=1:ng) FMInfo%w(ig) = nTracerCntl%UserUR
    IF(PE%MASTER) CALL PrintMOCUnderRelaxation(io8, FmInfo, GroupInfo, PE)
  ELSE
    FORALL(ig=1:ng) FMInfo%w(ig) = 1.0
  ENDIF
  RETURN
ENDIF

IF(.NOT. nTracerCntl%lOptUR) THEN
  FORALL(ig=1:ng) FMInfo%w(ig) = nTracerCntl%UserUR
  IF(PE%MASTER) CALL PrintMOCUnderRelaxation(io8, FmInfo, GroupInfo, PE)
  RETURN
ENDIF
ALLOCATE(C(ng), ScatSum(ng), w(ng, Core%nz), wbuf(ng, Core%nz))
FORALL (ig = 1:ng, iz=1:Core%nz) w(ig, iz) = 1.
DO iz = myzb, myze
  DO ixy = 1, nxy
     !Scattering
     ScatSum(1:ng) = 0
     DO ig = 1, ng
       ScatSum(ig) = ScatSum(ig) + PinXs(ixy, iz)%xss(ig)%WithInGroupScat
       DO ig2 = PinXs(ixy, iz)%xss(ig)%ib, PinXs(ixy, iz)%xss(ig)%ie
         ScatSum(ig2) = ScatSum(ig2) + PinXs(ixy, iz)%xss(ig)%from(ig2)
       ENDDO
     ENDDO
     FORALL (ig=1:ng) C(ig) = ScatSum(ig)/PinXs(ixy, iz)%xstr(ig)

     DO ig = 1, ng
       crit = PinXS(ixy, iz)%xstr(ig)* Core%hz(iz)
       f = 2._8 / SQRT(3 * c(ig))
       IF( f .LT. crit) THEN
          wlocal = 2._8 / (2 - c(ig)) 
       ELSE
          wlocal = 3._8 * crit * crit
          wlocal = wlocal / (2._8+3._8*(1._8-c(ig))*crit*crit)
       ENDIF
!       IF(ixy .EQ. 2) THEN
!          PRINT *, ig
!          PRINT *, 'c', c(ig), f, crit
!          print *, wlocal
!       ENDIF
       wlocal = max(0.4_8, wlocal)
       wlocal = min(1._8, wlocal)
       w(ig, iz) = min(w(ig, iz), wlocal)
     ENDDO
  ENDDO
ENDDO
#ifdef MPI_ENV

FORALL(ig=1:ng, iz=1:core%nz) wbuf(ig, iz) = 0
FORALL(ig=1:ng, iz=myzb:myze) wbuf(ig, iz) = w(ig, iz)
FORALL(ig=1:ng, iz=1:core%nz) w(ig, iz) = 0 
CALL REDUCE(wbuf, w, ng, Core%nz, PE%MPI_CMFD_COMM, .TRUE.)
#endif
FORALL(ig=1:ng) FmInfo%w(ig) = 1.0

DO ig = 1, ng
  DO iz = 1, Core%nz
    FmInfo%w(ig) = min(FmInfo%w(ig), w(ig, iz)) 
  ENDDO
ENDDO
!DO ig = 1, ng
!  PRINT *, ig,  FmInfo%w(ig)
!  FmInfo%w(ig) = 0.75
!ENDDO
IF(PE%MASTER) CALL PrintMOCUnderRelaxation(io8, FmInfo, GroupInfo, PE)
DEALLOCATE(w, wbuf, scatsum)
END SUBROUTINE

!--- CNJ Edit : Domain Decomposition
SUBROUTINE DcmpLinkBoundaryFlux(CoreInfo, RayInfo, PhiAngIn, DcmpPhiAngIn, DcmpPhiAngOut, gb, ge, color)
USE PARAM
USE TYPEDEF,    ONLY : CoreInfo_Type,   RayInfo_Type,    DcmpAsyRayInfo_Type
USE PE_Mod,     ONLY : PE
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(RayInfo_Type) :: RayInfo
REAL, POINTER :: PhiAngIn(:, :, :)
REAL, POINTER :: DcmpPhiAngIn(:, :, :, :, :), DcmpPhiAngOut(:, :, :, :, :)
INTEGER :: gb, ge
INTEGER, OPTIONAL :: color

TYPE(DcmpAsyRayInfo_Type), POINTER :: DcmpAsyRay(:, :)
INTEGER, POINTER :: DcmpAsyLinkInfo(:, :, :, :)
INTEGER, POINTER :: DcmpAsyRayCount(:)
INTEGER :: iAsy, iRay, irot
INTEGER :: nextAsy, nextRay

DcmpAsyRay => RayInfo%DcmpAsyRay
DcmpAsyLinkInfo => RayInfo%DcmpAsyLinkInfo
DcmpAsyRayCount => RayInfo%DcmpAsyRayCount

DO iAsy = 1, CoreInfo%nxya
  IF (PRESENT(color)) THEN
    IF (CoreInfo%Asy(iAsy)%color .NE. color) CYCLE
  ENDIF
  DO iRay = 1, DcmpAsyRayCount(iAsy)
    DO irot = 1, 2
      nextAsy = DcmpAsyLinkInfo(1, irot, iRay, iAsy)
      nextRay = DcmpAsyLinkInfo(2, irot, iRay, iAsy)
      IF (nextAsy .EQ. 0) THEN
        PhiAngIn(:, gb : ge, nextRay) = DcmpPhiAngOut(:, gb : ge, irot, iRay, iAsy)  
      ELSE
        DcmpPhiAngIn(:, gb : ge, irot, nextRay, nextAsy) = DcmpPhiAngOut(:, gb : ge, irot, iRay, iAsy)
      ENDIF
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE

!--- CNJ Edit : Domain Decomposition + MPI    
SUBROUTINE DcmpScatterXS(CoreInfo, xst)
USE PARAM
USE TYPEDEF,    ONLY : CoreInfo_Type
USE GEOM,       ONLY : ng
USE PE_MOD,     ONLY : PE
USE CNTL,       ONLY : nTracerCntl
USE ALLOCS
IMPLICIT NONE

INCLUDE 'mpif.h'

TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: xst(:, :)
REAL, POINTER :: buf_xst(:, :)

INTEGER :: nFsr, nMoment, nDat
INTEGER :: myPinBeg, myPinEnd, myFsrBeg, myFsrEnd
INTEGER :: sendcounts(0 : PE%nRTProc - 1), displs(0 : PE%nRTProc - 1)
INTEGER :: ierr

nFsr = PE%nFsr(PE%myRTRank)
IF (nTracerCntl%ScatOd .EQ. 1) nMoment = 2
IF (nTracerCntl%ScatOd .EQ. 2) nMoment = 5
IF (nTracerCntl%ScatOd .EQ. 3) nMoment = 9

CALL Dmalloc(buf_xst, ng, nFsr)
sendcounts = ng * PE%nFsr
displs = ng * PE%Fsr_displs
nDat = sendcounts(PE%myRTRank)
CALL MPI_SCATTERV(xst, sendcounts, displs, MPI_DOUBLE_PRECISION, buf_xst, nDat,                                     &
                  MPI_DOUBLE_PRECISION, 0, PE%MPI_RT_COMM, ierr)
xst(:, PE%myFsrBeg : PE%myFsrEnd) = buf_xst
DEALLOCATE(buf_xst)

END SUBROUTINE

SUBROUTINE DcmpScatterSource(CoreInfo, src, srcm)
USE PARAM
USE TYPEDEF,    ONLY : CoreInfo_Type
USE GEOM,       ONLY : ng
USE PE_MOD,     ONLY : PE
USE CNTL,       ONLY : nTracerCntl
USE ALLOCS
IMPLICIT NONE

INCLUDE 'mpif.h'

TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: src(:, :)
REAL, POINTER, OPTIONAL :: srcm(:, :, :)
REAL, POINTER :: buf_src(:, :), buf_srcm(:, :, :)

INTEGER :: nFsr, nMoment, nDat
INTEGER :: myPinBeg, myPinEnd, myFsrBeg, myFsrEnd
INTEGER :: sendcounts(0 : PE%nRTProc - 1), displs(0 : PE%nRTProc - 1)
INTEGER :: ierr

nFsr = PE%nFsr(PE%myRTRank)
IF (nTracerCntl%ScatOd .EQ. 1) nMoment = 2
IF (nTracerCntl%ScatOd .EQ. 2) nMoment = 5
IF (nTracerCntl%ScatOd .EQ. 3) nMoment = 9

CALL Dmalloc(buf_src, ng, nFsr)
sendcounts = ng * PE%nFsr
displs = ng * PE%Fsr_displs
nDat = sendcounts(PE%myRTRank)
CALL MPI_SCATTERV(src, sendcounts, displs, MPI_DOUBLE_PRECISION, buf_src, nDat,                                     &
                  MPI_DOUBLE_PRECISION, 0, PE%MPI_RT_COMM, ierr)
src(:, PE%myFsrBeg : PE%myFsrEnd) = buf_src
DEALLOCATE(buf_src)

IF (PRESENT(srcm)) THEN
  IF (nTracerCntl%lScat1) THEN 
    CALL Dmalloc(buf_srcm, nMoment, ng, nFsr)
    sendcounts = nMoment * ng * PE%nFsr
    displs = nMoment * ng * PE%Fsr_displs
    nDat = sendcounts(PE%myRTRank)
    CALL MPI_SCATTERV(srcm, sendcounts, displs, MPI_DOUBLE_PRECISION, buf_srcm, nDat,                               &
                      MPI_DOUBLE_PRECISION, 0, PE%MPI_RT_COMM, ierr)
    srcm(:, :, PE%myFsrBeg : PE%myFsrEnd) = buf_srcm
    DEALLOCATE(buf_srcm)
  ENDIF
ENDIF
    
END SUBROUTINE
    
SUBROUTINE DcmpScatterBoundaryFlux(RayInfo, PhiAngIn, DcmpPhiAngIn)
USE PARAM
USE TYPEDEF,    ONLY : RayInfo_Type
USE GEOM,       ONLY : ng
USE PE_MOD,     ONLY : PE
USE CNTL,       ONLY : nTracerCntl
USE ALLOCS
IMPLICIT NONE

INCLUDE 'mpif.h'

TYPE(RayInfo_Type) :: RayInfo
REAL, POINTER :: PhiAngIn(:, :, :), DcmpPhiAngIn(:, :, :, :, :)
REAL, POINTER :: buf_DcmpPhiAngIn(:, :, :, :, :)

INTEGER :: nPolarAngle, nModRay, nPhiAngSv, nAsy, nDat
INTEGER :: sendcounts(0 : PE%nRTProc - 1), displs(0 : PE%nRTProc - 1)
INTEGER :: ierr

nPolarAngle = RayInfo%nPolarAngle
nModRay = RayInfo%nModRay
nPhiAngSv = RayInfo%nPhiAngSv
nAsy = PE%nAsy(PE%myRTRank)

nDat = nPolarAngle * ng * nPhiAngSv
CALL MPI_BCAST(PhiAngIn, nDat, MPI_DOUBLE_PRECISION, 0, PE%MPI_RT_COMM, ierr)

CALL Dmalloc(buf_DcmpPhiAngIn, nPolarAngle, ng, 2, nModRay, nAsy)
sendcounts = nPolarAngle * ng * 2 * nModRay * PE%nAsy
displs = nPolarAngle * ng * 2 * nModRay * PE%Asy_displs
nDat = sendcounts(PE%myRTRank)
CALL MPI_SCATTERV(DcmpPhiAngIn, sendcounts, displs, MPI_DOUBLE_PRECISION, buf_DcmpPhiAngIn, nDat,                   &
                  MPI_DOUBLE_PRECISION, 0, PE%MPI_RT_COMM, ierr)
DcmpPhiAngIn(:, :, :, :, PE%myAsyBeg : PE%myAsyEnd) = buf_DcmpPhiAngIn
DEALLOCATE(buf_DcmpPhiAngIn)

END SUBROUTINE

SUBROUTINE DcmpGatherFlux(CoreInfo, phis, phim)
USE PARAM
USE TYPEDEF,    ONLY : CoreInfo_Type
USE GEOM,       ONLY : ng
USE PE_MOD,     ONLY : PE
USE CNTL,       ONLY : nTracerCntl
USE ALLOCS
IMPLICIT NONE

INCLUDE 'mpif.h'

TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: phis(:, :)
REAL, POINTER, OPTIONAL :: phim(:, :, :)
REAL, POINTER :: buf_phis(:, :), buf_phim(:, :, :)

INTEGER :: nFsr, nMoment, nDat
INTEGER :: recvcounts(0 : PE%nRTProc - 1), displs(0 : PE%nRTProc - 1)
INTEGER :: ierr

nFsr = PE%nFsr(PE%myRTRank)
IF (nTracerCntl%ScatOd .EQ. 1) nMoment = 2
IF (nTracerCntl%ScatOd .EQ. 2) nMoment = 5
IF (nTracerCntl%ScatOd .EQ. 3) nMoment = 9

CALL Dmalloc(buf_phis, ng, nFsr)
recvcounts = ng * PE%nFsr
displs = ng * PE%Fsr_displs
nDat = recvcounts(PE%myRTRank)
buf_phis = phis(:, PE%myFsrBeg : PE%myFsrEnd)
CALL MPI_GATHERV(buf_phis, nDat, MPI_DOUBLE_PRECISION, phis, recvcounts, displs,                                    &
                 MPI_DOUBLE_PRECISION, 0, PE%MPI_RT_COMM, ierr)
DEALLOCATE(buf_phis)

IF (PRESENT(phim)) THEN
  IF (nTracerCntl%lScat1) THEN 
    CALL Dmalloc(buf_phim, nMoment, ng, nFsr)
    recvcounts = nMoment * ng * PE%nFsr
    displs = nMoment * ng * PE%Fsr_displs
    nDat = recvcounts(PE%myRTRank)
    buf_phim = phim(:, :, PE%myFsrBeg : PE%myFsrEnd)
    CALL MPI_GATHERV(buf_phim, nDat, MPI_DOUBLE_PRECISION, phim, recvcounts, displs,                                &
                     MPI_DOUBLE_PRECISION, 0, PE%MPI_RT_COMM, ierr)
    DEALLOCATE(buf_phim)
  ENDIF
ENDIF

END SUBROUTINE

SUBROUTINE DcmpGatherCurrent(CoreInfo, jout)
USE PARAM
USE TYPEDEF,    ONLY : CoreInfo_Type
USE GEOM,       ONLY : ng
USE PE_MOD,     ONLY : PE
USE CNTL,       ONLY : nTracerCntl
USE ALLOCS
IMPLICIT NONE

INCLUDE 'mpif.h'

TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: jout(:, :, :, :)
REAL, POINTER :: buf_jout(:, :, :, :)

INTEGER :: nPin, nDat
INTEGER :: recvcounts(0 : PE%nRTProc - 1), displs(0 : PE%nRTProc - 1)
INTEGER :: ierr

nPin = PE%nPin(PE%myRTRank)

CALL Dmalloc(buf_jout, 3, ng, 4, nPin)
recvcounts = 3 * ng * 4 * PE%nPin
displs = 3 * ng * 4 * PE%Pin_displs
nDat = recvcounts(PE%myRTRank)
buf_jout = jout(:, :, :, PE%myPinBeg : PE%myPinEnd)
CALL MPI_GATHERV(buf_jout, nDat, MPI_DOUBLE_PRECISION, jout, recvcounts, displs,                                    &
                 MPI_DOUBLE_PRECISION, 0, PE%MPI_RT_COMM, ierr)
DEALLOCATE(buf_jout)

END SUBROUTINE
    
SUBROUTINE DcmpGatherBoundaryFlux(RayInfo, DcmpPhiAngOut)
USE PARAM
USE TYPEDEF,    ONLY : RayInfo_Type
USE GEOM,       ONLY : ng
USE PE_MOD,     ONLY : PE
USE CNTL,       ONLY : nTracerCntl
USE ALLOCS
IMPLICIT NONE

INCLUDE 'mpif.h'

TYPE(RayInfo_Type) :: RayInfo
REAL, POINTER :: DcmpPhiAngOut(:, :, :, :, :)
REAL, POINTER :: buf_DcmpPhiAngOut(:, :, :, :, :)

INTEGER :: nPolarAngle, nModRay, nAsy, nDat
INTEGER :: recvcounts(0 : PE%nRTProc - 1), displs(0 : PE%nRTProc - 1)
INTEGER :: ierr

nPolarAngle = RayInfo%nPolarAngle
nModRay = RayInfo%nModRay
nAsy = PE%nAsy(PE%myRTRank)

CALL Dmalloc(buf_DcmpPhiAngOut, nPolarAngle, ng, 2, nModRay, nAsy)
recvcounts = nPolarAngle * ng * 2 * nModRay * PE%nAsy
displs = nPolarAngle * ng * 2 * nModRay * PE%Asy_displs
nDat = recvcounts(PE%myRTRank)
buf_DcmpPhiAngOut = DcmpPhiAngOut(:, :, :, :, PE%myAsyBeg : PE%myAsyEnd)
CALL MPI_GATHERV(buf_DcmpPhiAngOut, nDat, MPI_DOUBLE_PRECISION, DcmpPhiAngOut, recvcounts, displs,                  &
                 MPI_DOUBLE_PRECISION, 0, PE%MPI_RT_COMM, ierr)
DEALLOCATE(buf_DcmpPhiAngOut)

END SUBROUTINE    