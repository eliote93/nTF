#include <defines.h>
SUBROUTINE GcCmfdAcc(Core, CmInfo, Eigv, lGcFdk, GroupInfo, GcGroupInfo, nTracerCntl, ItrCntl, PE)
USE PARAM
USE TYPEDEF,       ONLY : CoreInfo_Type,              CmInfo_Type,         GroupInfo_Type,      &
                          PE_TYPE,                    PinXS_Type
USE GcCmfd_mod,    ONLY : GenGcHomoXs,                GcRadCouplingCoeff,  GcAxCouplingCoeff,   &
                          SetGcCmfdLinearSystem,      SetGcCmfdEnv,        GcCmfdPsiUpdt,       &
                          GcCmfdSrcUpdt,              GcCmfdEigUpdt,       GcResidualError,     &
                          MgCmfdSolUpdt
USE CMFD_mod,      ONLY : PhiC1g,                     SRC
#ifdef MPI_ENV
USE MpiComm_mod,    ONLY : MPI_SYNC
#endif
USE BiCGSTAB_mod,   ONLY : BiCGSTAB
USE BasicOperation, ONLY : CP_CA, CP_VA
USE IOUTIL,         ONLY : message
USE FILES,          ONLY : io8
USE Cntl,           ONLY : nTracerCntl_Type
USE ItrCntl_mod,    ONLY : ItrCntl_TYPE
USE BiLU_Mod,       ONLY : MakeBilU
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(GroupInfo_TYpe) :: GroupInfo, GcGroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(ItrCntl_TYPE) :: ItrCntl
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

INTEGER :: nIterMax, nIterMin
REAL :: Conv_resCrit, Conv_eigvCrit
INTEGER :: ConvMod 
LOGICAL :: lConv, lExit, lLogOut

LOGICAL :: l3dim
LOGICAL :: CMFDMaster, CmfdSlave

IF(GcGroupInfo%ng .eq. 2) THEN
  CALL Cmfd2gAcc(Core, CmInfo, Eigv, lGcFdk, GroupInfo, GcGroupInfo, nTracerCntl, ItrCntl, PE)
  RETURN
ENDIF

l3dim = nTracerCntl%l3dim
CMFDMaster = PE%CMFDMaster; CMFDSlave = PE%CMFDSlave

ng = GroupInfo%ng; ngc = GcGroupInfo%ng
myzbf = PE%myzbf; myzef = PE%myzef
nxy = Core%nxy
PinXS => CmInfo%PinXS; GcPinXs => CmInfo%GcPinXS

PhiFm => CmInfo%PhiFm
GcPhiC => CMInfo%GcPhiC
GcPsiC => CMInfo%GcPsiC; GcPsiCD => CMInfo%GcPsiCD

nitermax = ItrCntl%GcCmfdItrCntl%nitermax
nitermin = ItrCntl%GcCmfdItrCntl%nitermin
conv_rescrit = ItrCntl%GcCmfdItrCntl%conv_rescrit
Conv_eigvCrit = ItrCntl%GcCmfdItrCntl%Conv_eigvCrit
ConvMod = ItrCntl%GcCmfdItrCntl%convMod
lLogOut = ItrCntl%GcCmfdItrCntl%lLogOut
!IF(ItrCntl%GcCmfdIt .EQ. 0) eigv = 1.0

CALL SetGcCmfdEnv(CMInfo, ngc)

WRITE(mesg,'(a)') 'Group Condensing (C)...'
IF(CMFDMaster .AND. lLogOut) CALL message(io8, TRUE, TRUE, mesg)      

CALL GenGcHomoXS(Core,PinXs, GcPinXs, PhiFm, GcPhiC, ng, ngc, GroupInfo, GcGroupInfo, PE)
CALL GcRadCouplingCoeff(Core, PinXS, GcPinXS, PhiFm, ng, ngc, GroupInfo, PE)
IF(l3dim) CALL GcAxCouplingCoeff(Core, GcPinXS, PhiFm, GcPhiC, CmInfo%AxDtil, CmInfo%AxDhat, ng, ngc, GroupInfo, PE)

CALL SetGcCmfdLinearSystem(.TRUE., nTracerCntl%l3dim)
CALL MakeBiLU(Core, CmInfo%GcCMFDLS(1:ngc), GcPhiC, ngc,  PE)

WRITE(mesg, '(a)') 'Performing GCCMFD Calculation...'
IF(CMFDMaster .AND. lLogOut) CALL message(io8, TRUE, TRUE, mesg)      


CALL GcCmfdPsiUpdt(GcPhiC, GcPsiC)

lconv = FALSE; 
nGroupInfo = 2;
IF(.NOT. GcGroupInfo%lUpScat) nGroupInfo = 1
ResErr = GcResidualError(GcPhiC, GcPsiC, eigv,  PE)
DO iter = 1, nitermax
  ItrCntl%GcCmfdIt = ItrCntl%GcCmfdIt + 1
  DO jsweep = 1, nGroupInfo
    GrpBeg = 1; GrpEnd = ngc
    IF(jsweep .GT. 1) THEN
      GrpBeg = GcGroupInfo%UpScatRange(1); GrpEnd = GcGroupInfo%UpScatRange(2)
    ENDIF
    DO ig = GrpBeg, GrpEnd
      CALL GcCmfdSrcUpdt(SRC, GcPsiC, GcPhiC, eigv, ig)
      CALL CP_VA(Phic1g(1:nxy, myzbf : myzef), GcPhiC(1:nxy, myzbf : myzef, ig), nxy, myzef - myzbf + 1)  
      CALL BiCGSTAB(CMInfo%GcCmfdLs(ig), Phic1g, SRC, itrcntl%InSolverItrCntl)
#ifndef MPI_ENV
      CALL CP_VA(GcPhiC(1:nxy, myzbf : myzef, ig), Phic1g(1:nxy, myzbf : myzef), nxy, myzef - myzbf + 1)    
#else
      CALL CP_VA(GcPhiC(1:nxy, myzbf-1 : myzef+1, ig), Phic1g(1:nxy, myzbf-1 : myzef+1), nxy, myzef - myzbf + 3)    
#endif   
      CONTINUE   
    ENDDO
  ENDDO
  peigv = eigv
  CALL CP_VA(GcPsiCD(1 : nxy, myzbf : myzef), GcPsiC(1 : nxy, myzbf : myzef), nxy, myzef - myzbf + 1)    
  CALL GcCmfdPsiUpdt(GcPhiC, GcPsiC)
  CALL GcCmfdEigUpdt(GcPsiC, GcPsiCD, Eigv, psierr, PE)
  EigErr = (eigv - peigv)/eigv
  ResErr = GcResidualError(GcPhiC, GcPsiC, eigv,  PE)
  IF(iter .eq. 1) ResErr0 = ResErr
  IF(CMFDMaster) WRITE(mesg,'(6x,A9, I8, F17.6, 3x, F10.5, 1pe15.3)') 'CGOUTER', ItrCntl%GcCmfdit, eigv, ResErr/ResErr0, ResErr
  IF(CMFDMaster .AND. lLogOut) CALL message(io8, FALSE, TRUE, mesg)
  !Convergence Check
  lExit = .TRUE.; lConv = .FALSE.
  
  SELECT CASE(ConvMod)
    CASE(1)
      IF((ResErr/ResErr0) .lt. conv_rescrit ) lconv = TRUE
    CASE(2)
      IF(EigErr .LT. Conv_eigvCrit) lconv = .TRUE.
  END SELECT
  
  lExit = lExit .AND. lConv
  IF(iter .LE. nitermin) lExit = .FALSE.

  IF(lExit) EXIT   
  CONTINUE
ENDDO

IF(lGcFdk) CALL MgCmfdSolUpdt(CmInfo, GroupInfo, 0.5_8)

NULLIFY(PinXS, GcPinXS)
NULLIFY(PhiFm)
NULLIFY(GcPHiC, GcPsiC, GcPsiCD)
END SUBROUTINE

!
SUBROUTINE SetGcCmfdLinearSystem(lcmfd, l3dim)
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
INTEGER :: ig, iz, iz0, ixy, ibd
INTEGER, PARAMETER :: nbd = 4
REAL, POINTER :: Diag(:, :), RadOffDiag(:, :, :), AxOffDiag(:, :, :)
REAL :: Vol, neighVol, lmnt, offlmnt, dhat, dtil, alpha, dhmax

DO ig = 1, ngc
  Diag => GcCmfdLs(ig)%Diag
  RadOffDiag => GcCmfdLs(ig)%RadOffDiag
  AxOffDiag => GcCmfdLs(ig)%AxOffDiag
  CALL CP_CA(diag(:, myzbf:myzef), zero, nxy, myzef-myzbf + 1)
  CALL CP_CA(RadOffDiag(:, :, myzbf:myzef), zero, 4, nxy, myzef-myzbf + 1)
  CALL CP_CA(AxOffDiag(:, :, myzbf:myzef), zero, 2, nxy, myzef-myzbf + 1)    
  DO iz = myzbf, myzef
    iz0 = SubPlaneMap(iz)
    DO ixy = 1, nxy
      Vol = PinVolFm(ixy, iz0)
      lmnt = Vol*(GcPinXS(ixy, iz)%xsr(ig))
      DO ibd = 1, nbd
        dtil = GcPinXs(ixy, iz)%Dtil(ibd, ig)
        dhat = GcPinXs(ixy, iz)%Dhat(ibd, ig)
        lmnt = lmnt + (dtil - dhat) * hzfm(iz)
        offlmnt = -(dtil + dhat) * hzfm(iz)
        RadOffDiag(ibd, ixy, iz) = offlmnt          
      ENDDO
      Diag(ixy, iz) = lmnt
    ENDDO
  ENDDO
  IF(.NOT. l3dim) CYCLE
  DO iz = myzbf, myzef
    iz0 = SubPlaneMap(iz)
    DO ixy = 1, nxy
      vol = PinVolFm(ixy, iz0)/hzfm(iz)
      dtil = GcPinXS(ixy, iz)%axdtil(1, ig)
      dhat = GcPinXS(ixy, iz)%axdhat(1, ig) 
      lmnt = (dtil - dhat)
      offlmnt = -(dtil + dhat)*vol
      AxOffDiag(1, ixy, iz) = OffLmnt
   
      dtil = GcPinXS(ixy, iz)%axdtil(2, ig)
      dhat = GcPinXS(ixy, iz)%axdhat(2, ig)                   
      lmnt = lmnt + (dtil - dhat)
      offlmnt = -(dtil + dhat)*vol
      AxOffDiag(2, ixy, iz) = OffLmnt
      lmnt = lmnt * vol 
      Diag(ixy, iz) = Diag(ixy, iz) + lmnt            
    ENDDO
  ENDDO
ENDDO

NULLIFY(Diag, RadOffDiag, AxOffDiag)

END SUBROUTINE
