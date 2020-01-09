
#include <defines.h>
SUBROUTINE InitVariables()
USE PARAM
IMPLICIT NONE
CALL InitFluxVariables()
END SUBROUTINE

SUBROUTINE InitFluxVariables()
USE PARAM 
USE TYPEDEF,         ONLY : CoreInfo_Type,  Pin_Type,    Cell_Type,      XsMac_Type, &
                            RayInfo_Type,   FxrInfo_Type
USE CNTL,            ONLY : nTracerCntl, DcplControl
USE PE_mod,          ONLY : PE,             DcplPE            
USE GEOM,            ONLY : ng,                         &
                            Core,          Pin,         CellInfo,        hz
USE CORE_mod,        ONLY : phis,          PhiAngin,    FXR,                        &
                            psi,           psid,        psic,            psicd,     &
                            phic,          PhiFm,       Power,                      &
                            GroupInfo,     FmInfo,      CmInfo,                     &
                            nPhiAngSv,     nCoreFsr,    nCoreFxr
USE DcplCore_mod,    ONLY : DcplInfo,      DcplFmInfo,  DcplCmInfo
USE Rays,            ONLY : RayInfo                     
USE BenchXs,         ONLY : xsnfBen
USE MacXsLib_Mod,    ONLY : MacXsNf
USE XsUtil_mod,      ONLY : FreeXsMac
USE BasicOperation,  ONLY : CP_CA
USE SUbGrp_Mod,      ONLY : FxrChiGen
#ifdef MPI_ENV
USE MPICOMM_MOD,     ONLY : MPI_SYNC, REDUCEnBCAST
#endif
IMPLICIT NONE

TYPE(XsMac_Type) :: XsMac
TYPE(FxrInfo_Type), POINTER :: myFxr
INTEGER :: i, j, k
INTEGER :: myzb, myze, myzbf, myzef,nz
INTEGER :: iz, ipin, icel, ireg, ixsreg, itype, ig
INTEGER :: nxy, FsrIdxSt, FxrIdxSt, niso
REAL, POINTER :: XsMacNf(:)
REAL, POINTER :: PinVol(:, :)
REAL :: totvol, fuelvol, fuelvol0, fscell, volfsr, totfsvol, temp0(2), temp1(2)
REAL :: phiinit, Jinit
LOGICAL :: lxsLib

nxy = Core%nxy; nz = Core%nz
PinVol => Core%PinVol
myzb = PE%myzb; myze = PE%myze
myzbf = PE%myzbf; myzef = PE%myzef
lxsLib = nTracerCntl%lXsLib

IF(.not. lxsLib) Allocate(XsMacNf(ng))
fuelvol = 0; totvol = 0; totfsvol = 0; FuelVol0 = 0
DO iz = myzb, myze      !Plane Sweep
  DO ipin = 1, nxy      !Pin Sweep
    FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
    icel = Pin(ipin)%Cell(iz)
    fscell = 0
    DO j = 1, CellInfo(icel)%nFxr
      ixsreg = FxrIdxSt + j - 1
      volfsr = FXR(ixsreg, iz)%area*hz(iz)
      IF(lxsLib) THEN
        myFxr => FXR(ixsreg, iz)
        niso = myFxr%niso
        CALL MacXSNf(XsMac,niso, myFxr%temp, myFxr%idiso(1:niso), myFxr%pnum(1:niso),   &
                     1, ng, ng, 1._8, FALSE)
        XsMacNf => XsMac%XSMacNf
      ELSE
        ireg = CellInfo(icel)%MapFxr2FsrIdx(1,j)
        !itype = CellInfo(icel)%iReg(ireg)
        itype = Fxr(ixsreg, iz)%imix
        CALL XsNfBen(Itype, 1, ng, XsMacNf)           !Obtaining
      ENDIF
      fscell = fscell + sum(XsMacNf(1:ng))*volfsr
      totvol = totvol + volfsr
      IF(sum(XsMacNf(1:ng)) .GT. 0._8) THEN
        fuelvol0 = fuelvol0 + volfsr
      ENDIF
    ENDDO
    IF(fscell .gt. epsm3) fuelvol = fuelvol + PinVol(ipin, iz)
    totfsvol = totfsvol + fscell
  ENDDO
ENDDO

#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_RTmaster_COMM)
CALL REDUCEnBCAST(fuelvol, PE%MPI_CMFD_COMM, PE%MPI_CMFD_COMM, PE%lCMFDGrp, TRUE, TRUE)
CALL REDUCEnBCAST(fuelvol0, PE%MPI_CMFD_COMM, PE%MPI_CMFD_COMM, PE%lCMFDGrp, TRUE, TRUE)
CALL REDUCEnBCAST(totfsvol, PE%MPI_CMFD_COMM, PE%MPI_CMFD_COMM, PE%lCMFDGrp, TRUE, TRUE)
#endif
Core%FuelVolFm = Fuelvol0
Core%fuelvol = fuelvol
Core%TotVol = Totvol
phiinit = fuelvol/totfsvol
IF(fuelvol .EQ. 0._8 .OR. totfsvol .EQ. 0._8) phiinit = 1
Jinit = phiinit * 0.25_8
!print *, phiinit, totvol
!Scalar flux Initialize
CALL CP_CA(phis(1:nCoreFsr, myzb:myze, 1:ng), phiinit ,nCoreFsr, myze - myzb + 1, ng)
!Incomming angular flux initialization
CALL CP_CA(PhiAngin(:, :, myzb:myze, 1:ng), phiinit, RayInfo%nPolarAngle, nPhiAngSv, myze - myzb + 1, ng)
PhiAngin(:, 1, :, :) = 0 !Zeroincoming flux
CALL CP_CA(psi(1:nCoreFsr, myzb:myze), phiinit, nCoreFsr, myze - myzb + 1)
CALL CP_CA(psid(1:nCoreFsr, myzb:myze), phiinit, nCoreFsr, myze - myzb + 1)
CALL CP_CA(psic(1:nxy, myzb:myze), phiinit, nxy, myze - myzb + 1)
CALL CP_CA(psicd(1:nxy, myzb:myze), phiinit, nxy, myze - myzb + 1)
DO ig = 1, ng
  IF (nTracerCntl%lCMFD) CALL CP_CA(phic(1:nxy, myzb:myze, ig), phiinit, nxy, myze - myzb + 1)
  IF(nTracerCntl%lSubPlane) CALL CP_CA(PhiFm(1:nxy, myzbf:myzef, ig), phiinit, nxy, myzef - myzbf + 1)
ENDDO
!Power 
CALL CP_CA(Power(1:nCoreFsr, 1:nz), 1._8, nCoreFsr, nz)

CALL FxrChiGen(Core, Fxr, FmInfo, GroupInfo, PE, nTracerCntl)
!CALL CP_CA(phic(1:nxy, myzb:myze, 1:ng), phiinit, nxy, myze - myzb + 1, ng)
IF(nTracerCntl%lDcpl) THEN
#ifndef MPI_ENV
  DO i = 1, DcplInfo%nRefPln
#else
  DO i = PE%myRefPlnBeg, PE%myRefPlnEnd
#endif  
    myzb = DcplInfo%RefPln(i); myze = DcplInfo%RefPln(i)
    DO j = 1, DcplInfo%nRefTemp
      CALL CP_CA(DcplFmInfo(j, i)%Phis(1:nCoreFsr, myzb:myze, 1:ng), phiinit ,nCoreFsr, myze - myzb + 1, ng)
      CALL CP_CA(DcplFmInfo(j, i)%psi(1:nCoreFsr, myzb:myze), PhiInit, nCoreFsr, myze - myzb + 1)
      CALL CP_CA(DcplFmInfo(j, i)%psid(1:nCoreFsr, myzb:myze), PhiInit, nCoreFsr, myze - myzb + 1)
      CALL CP_CA(DcplFmInfo(j, i)%psic(1:nxy, myzb:myze), phiinit, nxy, myze - myzb + 1)
      CALL CP_CA(DcplFmInfo(j, i)%psicd(1:nxy, myzb:myze), phiinit, nxy, myze - myzb + 1)
      CALL FxrChiGen(Core, DcplFmInfo(j, i)%Fxr, DcplFmInfo(j, i), GroupInfo, DcplPE(i), nTracerCntl)
    ENDDO
  ENDDO
ENDIF
IF(lXsLib) CALL FreeXsMac(XsMac)
IF(.not. lxsLib .and. ASSOCIATED(XsMacNf)) THEN
  DEALLOCATE(XsMacNf)
ELSE
  IF(ASSOCIATED(XsMacNf)) NULLIFY(XsMacNf)
ENDIF
NULLIFY(MyFxr)
NULLIFY(XsMacNf)
NULLIFY(PinVol)
END SUBROUTINE

SUBROUTINE SetCorePower(Core, nTracerCntl)
USE PARAM
USE TYPEDEF,   ONLY : CoreInfo_Type,       Asy_Type,     AsyInfo_TYPE
USE Cntl,      ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(nTracerCntl_Type) :: nTracerCntl

TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)

INTEGER :: ixya, itype
INTEGER :: nxya
REAL :: FuelFA
REAL :: PowerCore, PowerFa 

Asy => Core%Asy; AsyInfo => Core%AsyInfo
nxya = Core%nxya
PowerFa = nTracerCntl%PowerFa
FuelFa = 0

DO ixya = 1, nxya
  itype = Asy(ixya)%AsyType
  IF(.NOT. AsyInfo(itype)%lFuel) CYCLE
  FuelFa = FuelFa + Asy(ixya)%wt
ENDDO 
PowerCore = FuelFa * PowerFA * epsm6 !W unit => MW 
Core%PowerCore = PowerCore
nTracerCntl%PowerCore = PowerCore
NULLIFY(ASY)

END SUBROUTINE

SUBROUTINE InitIterVar(Core, FMInfo, CmInfo, GroupInfo,  lDhatInit, ItrCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,      ONLY : CoreInfo_Type,     FMInfo_Type,    CmInfo_Type,    &
                         GroupInfo_Type,    PE_Type,                        &
                         PinXs_Type
USE ItrCntl_mod,  ONLY : ItrCntl_Type
USE Cntl,         ONLY : nTracerCntl_Type
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(FmInfo_Type) :: FmInfo
TYPE(CmInfo_Type) :: CmInfo
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(ItrCntl_Type) :: ItrCntl
TYPE(PE_Type) :: PE 

LOGICAL :: lDhatInit

TYPE(PinXs_Type), POINTER :: PinXS(:,:)

INTEGER :: ng, nfxr, nxy, myzb, myze
INTEGER :: ig, iz, ixy
ng = GroupInfo%ng
nxy = Core%nxy
nfxr = Core%nCoreFxr
myzb = PE%myzb; myze = PE%myze
PinXS => CmInfo%PinXS

IF(nTracerCntl%lBoronSearch) ItrCntl%eigconv = 1.0E-5_8

ItrCntl%CmfdIt0 = ItrCntl%Cmfdit
ItrCntl%GcCmfdIt0 = ItrCntl%GcCmfdit
ItrCntl%MocIt0 = ItrCntl%Mocit
ItrCntl%AxIt0 = ItrCntl%AxIt
ItrCntl%SrcIt0 = ItrCntl%SrcIt
IF (PE%RTMASTER) THEN   !--- CNJ Edit : Domain Decomposition + MPI
  IF (lDhatInit .AND. nTracerCntl%lCMFD) THEN
    PinXS => CmInfo%PinXS
    DO iz = myzb, myze
      DO ixy = 1, nxy
        PinXS(ixy, iz)%Dhat(:, :) = 0.0_8
        PinXS(ixy, iz)%pdhat(:, :) = 0.0_8
      ENDDO
    ENDDO
    NULLIFY(PinXS)
  ENDIF
ENDIF

END SUBROUTINE

SUBROUTINE InitDeplFluxVar()
USE PARAM 
USE TYPEDEF,         ONLY : CoreInfo_Type,  Pin_Type,    Cell_Type,      XsMac_Type, &
                            RayInfo_Type,   FxrInfo_Type
USE CNTL,            ONLY : nTracerCntl, DcplControl
USE PE_mod,          ONLY : PE,             DcplPE            
USE GEOM,            ONLY : ng,                         &
                            Core,          Pin,         CellInfo,        hz
USE CORE_mod,        ONLY : phis,          PhiAngin,    FXR,                        &
                            psi,           psid,        psic,            psicd,     &
                            phic,          PhiFm,       Power,                      &
                            GroupInfo,     FmInfo,      CmInfo,                     &
                            nPhiAngSv,     nCoreFsr,    nCoreFxr
USE Rays,            ONLY : RayInfo                     
USE BenchXs,         ONLY : xsnfBen
USE MacXsLib_Mod,    ONLY : MacXsNf
USE XsUtil_mod,      ONLY : FreeXsMac
USE BasicOperation,  ONLY : CP_CA
USE SUbGrp_Mod,      ONLY : FxrChiGen
USE FILES,            ONLY : io8
USE IOUTIL,           ONLY : message
#ifdef MPI_ENV
USE MPICOMM_MOD,     ONLY : MPI_SYNC, REDUCEnBCAST
#endif
IMPLICIT NONE

TYPE(XsMac_Type) :: XsMac
TYPE(FxrInfo_Type), POINTER :: myFxr
INTEGER :: i, j, k
INTEGER :: myzb, myze, myzbf, myzef,nz
INTEGER :: iz, ipin, icel, ireg, ixsreg, itype, ig
INTEGER :: nxy, FsrIdxSt, FxrIdxSt, niso
REAL, POINTER :: XsMacNf(:)
REAL, POINTER :: PinVol(:, :)
REAL :: totvol, fuelvol, fuelvol0, fscell, volfsr, totfsvol, temp0(2), temp1(2)
REAL :: phiinit, Jinit
LOGICAL :: lxsLib

nxy = Core%nxy; nz = Core%nz
PinVol => Core%PinVol
myzb = PE%myzb; myze = PE%myze
myzbf = PE%myzbf; myzef = PE%myzef
lxsLib = nTracerCntl%lXsLib

IF(PE%Master) THEN
   WRITE(mesg, '(A)') 'Initialized Flux Variables ...'
   CALL message(io8, TRUE, TRUE, mesg)
ENDIF

IF(.not. lxsLib) Allocate(XsMacNf(ng))
fuelvol = 0; totvol = 0; totfsvol = 0; FuelVol0 = 0
DO iz = myzb, myze      !Plane Sweep
  DO ipin = 1, nxy      !Pin Sweep
    FsrIdxSt = Pin(ipin)%FsrIdxSt; FxrIdxSt = Pin(ipin)%FxrIdxSt
    icel = Pin(ipin)%Cell(iz)
    fscell = 0
    DO j = 1, CellInfo(icel)%nFxr
      ixsreg = FxrIdxSt + j - 1
      volfsr = FXR(ixsreg, iz)%area*hz(iz)
      IF(lxsLib) THEN
        myFxr => FXR(ixsreg, iz)
        niso = myFxr%niso
        CALL MacXSNf(XsMac,niso, myFxr%temp, myFxr%idiso(1:niso), myFxr%pnum(1:niso),   &
                     1, ng, ng, 1._8, FALSE)
        XsMacNf => XsMac%XSMacNf
      ELSE
        ireg = CellInfo(icel)%MapFxr2FsrIdx(1,j)
        !itype = CellInfo(icel)%iReg(ireg)
        itype = Fxr(ixsreg, iz)%imix
        CALL XsNfBen(Itype, 1, ng, XsMacNf)           !Obtaining
      ENDIF
      fscell = fscell + sum(XsMacNf(1:ng))*volfsr
      totvol = totvol + volfsr
      IF(sum(XsMacNf(1:ng)) .GT. 0._8) THEN
        fuelvol0 = fuelvol0 + volfsr
      ENDIF
    ENDDO
    IF(fscell .gt. epsm3) fuelvol = fuelvol + PinVol(ipin, iz)
    totfsvol = totfsvol + fscell
  ENDDO
ENDDO

#ifdef MPI_ENV
CALL MPI_SYNC(PE%MPI_RTmaster_COMM)
CALL REDUCEnBCAST(fuelvol, PE%MPI_CMFD_COMM, PE%MPI_CMFD_COMM, PE%lCMFDGrp, TRUE, TRUE)
CALL REDUCEnBCAST(fuelvol0, PE%MPI_CMFD_COMM, PE%MPI_CMFD_COMM, PE%lCMFDGrp, TRUE, TRUE)
CALL REDUCEnBCAST(totfsvol, PE%MPI_CMFD_COMM, PE%MPI_CMFD_COMM, PE%lCMFDGrp, TRUE, TRUE)
#endif
Core%FuelVolFm = Fuelvol0
Core%fuelvol = fuelvol
Core%TotVol = Totvol
phiinit = fuelvol/totfsvol
IF(fuelvol .EQ. 0._8 .OR. totfsvol .EQ. 0._8) phiinit = 1
Jinit = phiinit * 0.25_8
!print *, phiinit, totvol
!Scalar flux Initialize
CALL CP_CA(phis(1:nCoreFsr, myzb:myze, 1:ng), phiinit ,nCoreFsr, myze - myzb + 1, ng)
!Incomming angular flux initialization
CALL CP_CA(PhiAngin(:, :, myzb:myze, 1:ng), phiinit, RayInfo%nPolarAngle, nPhiAngSv, myze - myzb + 1, ng)
PhiAngin(:, 1, :, :) = 0 !Zeroincoming flux
CALL CP_CA(psi(1:nCoreFsr, myzb:myze), phiinit, nCoreFsr, myze - myzb + 1)
CALL CP_CA(psid(1:nCoreFsr, myzb:myze), phiinit, nCoreFsr, myze - myzb + 1)
CALL CP_CA(psic(1:nxy, myzb:myze), phiinit, nxy, myze - myzb + 1)
CALL CP_CA(psicd(1:nxy, myzb:myze), phiinit, nxy, myze - myzb + 1)
DO ig = 1, ng
  CALL CP_CA(phic(1:nxy, myzb:myze, ig), phiinit, nxy, myze - myzb + 1)
  IF(nTracerCntl%lSubPlane) CALL CP_CA(PhiFm(1:nxy, myzbf:myzef, ig), phiinit, nxy, myzef - myzbf + 1)
ENDDO
!Power 
CALL CP_CA(Power(1:nCoreFsr, 1:nz), 1._8, nCoreFsr, nz)
END SUBROUTINE
