SUBROUTINE Cmfd2GAcc(Core, CmInfo, Eigv, lGcFdk, GroupInfo, GcGroupInfo, nTracerCntl, ItrCntl, PE)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,              CmInfo_Type,         GroupInfo_Type,      &
                           PE_TYPE,                    PinXS_Type
USE GcCmfd_mod,     ONLY : GenGcHomoXs,                GcRadCouplingCoeff,  GcAxCouplingCoeff,   &
                           SetGcCmfdEnv,               GcCmfdPsiUpdt,       GcCmfdEigUpdt,       &
                           GcResidualError,            MgCmfdSolUpdt
USE Cmfd2g_mod,     ONLY : Src2g,                      phic2g,                                   &
                           SetCmfd2GLinearSystem,      Cmfd2GSrcUpdt,       ResidualError2G,     &
                           WielandtUpdt,               WielandtLsShift,     NegativeFixUp2G
USE BICGSTAB2G_mod, ONLY : BiCGSTAB2G
#ifdef MPI_ENV
USE MpiComm_mod,    ONLY : MPI_SYNC
#endif
USE BiCGSTAB_mod,   ONLY : BiCGSTAB
USE BasicOperation, ONLY : CP_CA, CP_VA
USE IOUTIL,         ONLY : message
USE FILES,          ONLY : io8
USE Cntl,           ONLY : nTracerCntl_Type
USE ItrCntl_mod,    ONLY : GcCmfdItrCntl_TYPE, ItrCntl_Type
USE BiLU_Mod,       ONLY : MakeBilU2G
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(GroupInfo_TYpe) :: GroupInfo, GcGroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(ItrCntl_Type) :: ItrCntl
TYPE(PE_TYPE) :: PE

REAL :: Eigv
LOGICAL :: lGcFdk

REAL :: pEigv, PsiErr, ResErr, ResErr0, EigErr

TYPE(PinXS_Type), POINTER :: PinXs(:, :), GcPinXs(:, :)

REAL, POINTER :: GcPhiC(:, :, :), PhiFm(:, :, :)
REAL, POINTER :: GcPsiC(:, :), GcPsiCD(:, :)

INTEGER :: ngc, ng, nGroupInfo
INTEGER :: myzbf, myzef, nxy
INTEGER :: ig, iter, jsweep, GrpBeg, GrpEnd

REAL :: Eigvs, rEigvDel, reigv, reigvs, reigvsdel
INTEGER :: nIterMax, nIterMin
REAL :: Conv_ResCrit, Conv_EigvCrit
INTEGER :: ConvMod = 1
LOGICAL :: lConv, lExit, lLogOut

LOGICAL :: l3dim
LOGICAL :: CMFDMaster, CmfdSlave
LOGICAL :: lwielandt

INTEGER :: i, j, k, l

lwielandt = .TRUE.
l3dim = nTracerCntl%l3dim
CMFDMaster = PE%CMFDMaster; CMFDSlave = PE%CMFDSlave

ng = GroupInfo%ng; ngc = GcGroupInfo%ng
myzbf = PE%myzbf; myzef = PE%myzef
nxy = Core%nxy
PinXS => CmInfo%PinXS; GcPinXs => CmInfo%GcPinXS

PhiFm => CmInfo%PhiFm
GcPhiC => CMInfo%GcPhiC
GcPsiC => CMInfo%GcPsiC; GcPsiCD => CMInfo%GcPsiCD

nitermax = ItrCntl%GcCMFDItrCntl%nitermax
nitermin = ItrCntl%GcCMFDItrCntl%nitermin
conv_rescrit = ItrCntl%GcCMFDItrCntl%conv_rescrit
conv_eigvcrit = ItrCntl%GcCMFDItrCntl%conv_eigvcrit
lLogOut = ItrCntl%GcCMFDItrCntl%lLogOut
!convcrit = 0.01
!IF(ItrCntl%GcCmfdIt .EQ. 0) eigv = 1.0

CALL SetGcCmfdEnv(CMInfo, ngc)

WRITE(mesg,'(a)') 'Group Condensing (C)...'
IF(CMFDMaster .AND. lLogOut) CALL message(io8, TRUE, TRUE, mesg)

CALL GenGcHomoXS(Core,PinXs, GcPinXs, PhiFm, GcPhiC, ng, ngc, GroupInfo, GcGroupInfo, PE)
CALL GcRadCouplingCoeff(Core, PinXS, GcPinXS, PhiFm, ng, ngc, GroupInfo, PE)
IF(l3dim) CALL GcAxCouplingCoeff(Core, GcPinXS, PhiFm, GcPhiC, CmInfo%AxDtil, CmInfo%AxDhat, ng, ngc, GroupInfo, PE)

CALL SetCmfd2gLinearSystem(.TRUE., nTracerCntl%l3dim)

WRITE(mesg, '(a)') 'Performing 2GCMFD Calculation...'
IF(CMFDMaster .AND. lLogOut) CALL message(io8, TRUE, TRUE, mesg)

CALL GcCmfdPsiUpdt(GcPhiC, GcPsiC)

CALL ConvertArray2G(PhiC2G(:, :, :), GcPhiC(:, :, :), nxy, myzbf-1, myzef+1, 1)



rEigvDel = 1._8/(eigv); reigv = 1._8/eigv
ResErr = ResidualError2G(PhiC2G, GcPsiC, reigv, 0._8, PE)
#define wielandt
#ifdef wielandt
!initial veilandt shift
!eigvs = eigv + 0.001
!reigv = 1._8 / eigv;
!reigvs = 1._8 / eigvs
!reigvdel = reigv - reigvs;
!reigvsdel = reigvs
!call wielandtlsshift(cminfo%gccmfdls(1), reigv, reigvsdel, pe)
#endif
!lwielandt = .FALSE.
IF(lwielandt) THEN
  call wielandtlsshift(cminfo%gccmfdls(1), reigv, reigvdel, pe)
  CALL MakeBiLU2G(Core, CmInfo%GcCMFDLs(1), PhiC2G, 2, PE)
  call wielandtlsshift(cminfo%gccmfdls(1), reigv, -reigvdel, pe)
ELSE
  CALL MakeBiLU2G(Core, CmInfo%GcCMFDLs(1), PhiC2G, 2, PE)
ENDIF
DO iter = 1, nitermax
  !
  ItrCntl%GcCmfdIt = ItrCntl%GcCmfdIt + 1
  CALL Cmfd2GSrcUpdt(Src2g, GcPsiC, reigvdel)
  CALL BiCGSTAB2G(CmInfo%GcCmfdLs(1), PhiC2g, Src2G, itrcntl%InSolverItrCntl)
  CALL ConvertArray2G(PhiC2G(:, :, :), GcPhiC(:, :, :), nxy, myzbf-1, myzef+1, 2)

  peigv = eigv
  GcPsiCD(1:nxy, myzbf:myzef) = GcPsiC(1:nxy, myzbf:myzef)
  CALL GcCmfdPsiUpdt(GcPhiC, GcPsiC)
IF(.NOT. lwielandt) THEN
  CALL GcCmfdEigUpdt(GcPsiC, GcPsiCD, Eigv, psierr, PE)
  reigv = 1._8 / eigv; reigvs = 0
  reigvdel= reigv

ELSE
  CALL WielandtUpdt(CmInfo%GCCMFDLS(1), GcPsiC, GcPsiCd, Eigv, Eigvs, PsiErr, iter, PE)
  reigv = 1._8 / eigv; reigvs = 1._8 / eigvs
  reigvdel = reigv - reigvs;
  CALL NegativeFixUp2G(PhiC2G, GcPsiC, PE)
ENDIF

  EigErr = (eigv - peigv)/eigv
  ResErr = ResidualError2G(PhiC2G, GcPsiC, reigv, reigvs, PE)
  IF(iter .eq. 1) ResErr0 = ResErr

  IF(CMFDMaster) WRITE(mesg,'(6x,A9, I8, F17.6, 3x, F10.5, 1pe15.3)') '2GOUTER', ItrCntl%GcCmfdit, eigv, ResErr/ResErr0, ResErr
  IF(CMFDMaster .AND. lLogOut) CALL message(io8, FALSE, TRUE, mesg)
!  !Convergence Check
  lExit = .TRUE.; lConv = .FALSE.
  !IF((ResErr/ResErr0) .lt. convcrit ) lconv = TRUE
  SELECT CASE(ConvMod)
    CASE(1)
      IF((ResErr/ResErr0) .lt. conv_rescrit ) lconv = TRUE
    CASE(2)
      IF(EigErr .LT. Conv_eigvCrit) lconv = .TRUE.
  END SELECT

  lExit = lExit .AND. lConv
  IF(iter .LE. 4) lExit = .FALSE.
  IF(iter .GE. 4 .AND. ResErr .LT. 5.0E-8_8) lExit = .TRUE.
  IF(lExit) EXIT
ENDDO

IF(lGcFdk) CALL MgCmfdSolUpdt(CmInfo, GroupInfo, 0.5_8)

NULLIFY(PinXS, GcPinXS)
NULLIFY(PhiFm)
NULLIFY(GcPHiC, GcPsiC, GcPsiCD)

END SUBROUTINE

SUBROUTINE SetCmfd2GLinearSystem(lcmfd, l3dim)
USE PARAM
USE TYPEDEF,     ONLY : CoreInfo_Type,         CmfdLs_Type,         PinXS_Type
USE CMFD_MOD,    ONLY : nxy,                  myzb,                myze,                &
                        myzbf,                myzef,                                                      &
                        hzfm,                 hz,                  SubPlaneMap,      SubPlaneRange,       &
                        PinNeighIdx,          PinVol,              PinVolFm
USE GcCMFD_mod,  ONLY : GcPinXs,              GcCmfdLs,            ngc
USE BASICOPERATION, ONLY : CP_CA
IMPLICIT NONE
LOGICAL :: lcmfd, l3dim
INTEGER :: ig, ig2, iz, iz0, ixy, ibd
INTEGER, PARAMETER :: nbd = 4
INTEGER, PARAMETER :: IDX(2,2)= RESHAPE((/1, 3, 2, 4/), SHAPE(IDX))
REAL, POINTER :: Diag(:, :, :), RadOffDiag(:, :, :, :), AxOffDiag(:, :, :, :)
REAL :: Vol, neighVol, lmnt, offlmnt, dhat, dtil, alpha, dhmax

!DO ig = 1, ngc
Diag => GcCmfdLs(1)%Diag2g
RadOffDiag => GcCmfdLs(1)%RadOffDiag2g
AxOffDiag => GcCmfdLs(1)%AxOffDiag2g
diag(1:4,1:nxy,myzbf:myzef) = ZERO
RadOffDiag(1:4,1:4,1:nxy, myzbf:myzef) = ZERO
AxOffDiag(1:4,1:2,1:nxy,myzbf:myzef) = ZERO
DO iz = myzbf, myzef
  iz0 = SubPlaneMap(iz)
  DO ixy = 1, nxy
    Vol = PinVolFm(ixy, iz0)
    DO ig = 1, 2
      lmnt = Vol*(GcPinXS(ixy, iz)%xsr(ig))
      DO ibd = 1, nbd
        dtil = GcPinXs(ixy, iz)%Dtil(ibd, ig)
        dhat = GcPinXs(ixy, iz)%Dhat(ibd, ig)
        lmnt = lmnt + (dtil - dhat) * hzfm(iz)
        offlmnt = -(dtil + dhat) * hzfm(iz)
        RadOffDiag(idx(ig, ig),ibd, ixy, iz) = offlmnt
      ENDDO
      Diag(idx(ig, ig), ixy, iz) = lmnt
    ENDDO
  ENDDO
ENDDO
IF(l3dim) THEN
  DO iz = myzbf, myzef
    iz0 = SubPlaneMap(iz)
    DO ixy = 1, nxy
      vol = PinVolFm(ixy, iz0)/hzfm(iz)
      DO ig = 1, 2
        dtil = GcPinXS(ixy, iz)%axdtil(1, ig)
        dhat = GcPinXS(ixy, iz)%axdhat(1, ig)
        lmnt = (dtil - dhat)
        offlmnt = -(dtil + dhat)*vol
        AxOffDiag(IDX(ig, ig),1, ixy, iz) = OffLmnt

        dtil = GcPinXS(ixy, iz)%axdtil(2, ig)
        dhat = GcPinXS(ixy, iz)%axdhat(2, ig)
        lmnt = lmnt + (dtil - dhat)
        offlmnt = -(dtil + dhat)*vol
        AxOffDiag(IDX(ig, ig), 2, ixy, iz) = OffLmnt
        lmnt = lmnt * vol
        Diag(IDX(ig, ig), ixy, iz) = Diag(IDX(ig, ig), ixy, iz) + lmnt
      ENDDO
    ENDDO
  ENDDO
ENDIF
!Move Scattering Matrix
DO iz = myzbf, myzef
  iz0 = SubPlaneMap(iz)
  DO ixy = 1, nxy
    vol =  PinVolFm(ixy, iz0)
    DO ig = 1, 2
      DO ig2 = GcPinXs(ixy, iz)%xss(ig)%ib, GcPinXs(ixy, iz)%xss(ig)%ie
        Diag(IDX(ig, ig2), ixy, iz) = Diag(IDX(ig, ig2), ixy, iz) - vol * GcPinXs(ixy, iz)%xss(ig)%from(ig2)
      ENDDO
    ENDDO
  ENDDO
ENDDO
NULLIFY(Diag, RadOffDiag, AxOffDiag)

END SUBROUTINE
