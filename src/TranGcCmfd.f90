#include <defines.h>
SUBROUTINE TranGcCmfdAcc(Core, CmInfo, TranInfo, Eigv, lGcFdk, GroupInfo, GcGroupInfo, TranCntl, nTracerCntl, ItrCntl, PE)
USE PARAM
USE TYPEDEF,       ONLY : CoreInfo_Type,              CmInfo_Type,         GroupInfo_Type,      &
                          TranInfo_Type,              TranCntl_Type,                            &
                          PE_TYPE,                    PinXS_Type
USE GcCmfd_mod,    ONLY : GenGcHomoXs,                GcRadCouplingCoeff,  GcAxCouplingCoeff,   &
                          SetGcCmfdLinearSystem,      SetGcCmfdEnv,        GcCmfdPsiUpdt,       &
                          GcCmfdSrcUpdt,              GcCmfdEigUpdt,       GcResidualError,     &
                          MgCmfdSolUpdt
USE CMFD_mod,      ONLY : PhiC1g,                     SRC
USE TranGcCmfd_mod,ONLY : MakeGcKineticParam,                                                   &
                          TranGcResidualError,        SetTranGcCmfdLinearSystem,                &
                          SetTranCmfd2GLinearSystem
#ifdef MPI_ENV
USE MpiComm_mod,    ONLY : MPI_SYNC
#endif
USE BiCGSTAB_mod,   ONLY : BiCGSTAB
USE BasicOperation, ONLY : CP_CA, CP_VA, AD_VA, MULTI_CA
USE IOUTIL,         ONLY : message
USE FILES,          ONLY : io8
USE Cntl,           ONLY : nTracerCntl_Type
USE ItrCntl_mod,    ONLY : ItrCntl_TYPE
USE BiLU_Mod,       ONLY : MakeBilU,                  MakeBilU2G
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_TYpe) :: GroupInfo, GcGroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(ItrCntl_TYPE) :: ItrCntl
TYPE(PE_TYPE) :: PE

REAL :: Eigv
LOGICAL :: lGcFdk 
REAL :: pEigv, PsiErr, ResErr, ResErr0, EigErr

TYPE(PinXS_Type), POINTER :: PinXs(:, :), GcPinXs(:, :)

REAL, POINTER :: GcPhiC(:, :, :), PhiFm(:, :, :)
REAL, POINTER :: GcPsiC(:, :), GcPsiCD(:, :)
REAL, POINTER :: GcTranSrc(:, :, :)

INTEGER :: ngc, ng, nGroupInfo
INTEGER :: myzbf, myzef, nxy
INTEGER :: ig, iter, jsweep, GrpBeg, GrpEnd

INTEGER :: nIterMax, nIterMin
REAL :: Conv_resCrit, Conv_eigvCrit
INTEGER :: ConvMod 
LOGICAL :: lConv, lExit, lLogOut

LOGICAL :: l3dim
LOGICAL :: CMFDMaster, CmfdSlave

l3dim = nTracerCntl%l3dim
CMFDMaster = PE%CMFDMaster; CMFDSlave = PE%CMFDSlave

ng = GroupInfo%ng; ngc = GcGroupInfo%ng
myzbf = PE%myzbf; myzef = PE%myzef
nxy = Core%nxy
PinXS => CmInfo%PinXS; GcPinXs => CmInfo%GcPinXS

PhiFm => CmInfo%PhiFm
GcPhiC => CMInfo%GcPhiC
GcPsiC => CMInfo%GcPsiC; GcPsiCD => CMInfo%GcPsiCD
GcTranSrc => CmInfo%GcTranSrc
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
CALL MakeGcKineticParam(Core, CmInfo, TranInfo, GroupInfo, GcGroupInfo, TranCntl, nTracerCntl, PE)

CALL SetGcCmfdLinearSystem(.TRUE., nTracerCntl%l3dim)
CALL SetTranGcCmfdLinearSystem()
CALL MakeBiLU(Core, CmInfo%GcCMFDLS(1:ngc), GcPhiC, ngc,  PE)




WRITE(mesg, '(a)') 'Performing GCCMFD Calculation...'
IF(CMFDMaster .AND. lLogOut) CALL message(io8, TRUE, TRUE, mesg)      


CALL GcCmfdPsiUpdt(GcPhiC, GcPsiC)
CALL MULTI_CA(1._8 / eigv, GcPsiC(1:nxy, myzbf:myzef), nxy, myzef - myzbf + 1)

lconv = FALSE; 
nGroupInfo = 2;
IF(.NOT. GcGroupInfo%lUpScat) nGroupInfo = 1
ResErr = TranGcResidualError(GcPhiC, GcPsiC, GcTranSrc, PE)
DO iter = 1, nitermax
  CALL UpdtGcTranSrc(Core, CmInfo, TranInfo, GroupInfo, GcGroupInfo, TranCntl, nTracerCntl, PE)
  ItrCntl%GcCmfdIt = ItrCntl%GcCmfdIt + 1
  DO jsweep = 1, nGroupInfo
    GrpBeg = 1; GrpEnd = ngc
    IF(jsweep .GT. 1) THEN
      GrpBeg = GcGroupInfo%UpScatRange(1); GrpEnd = GcGroupInfo%UpScatRange(2)
    ENDIF
    DO ig = GrpBeg, GrpEnd
      CALL GcCmfdSrcUpdt(SRC, GcPsiC, GcPhiC, 1._8, ig)
      CALL AD_VA(SRC(1:nxy, myzbf:myzef), Src(1:nxy, myzbf:myzef), GcTranSrc(1:nxy, myzbf:myzef, ig), nxy, myzef - myzbf + 1)
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

  CALL CP_VA(GcPsiCD(1 : nxy, myzbf : myzef), GcPsiC(1 : nxy, myzbf : myzef), nxy, myzef - myzbf + 1)    
  CALL GcCmfdPsiUpdt(GcPhiC, GcPsiC)
  CALL MULTI_CA(1._8 / eigv, GcPsiC(1:nxy, myzbf:myzef), nxy, myzef - myzbf + 1)
  !CALL CP_CA(GcTranSrc(1:nxy, myzbf:myzef, 1:ngc), 0._8, nxy, myzef-myzbf+1, ngc)
  ResErr = TranGcResidualError(GcPhiC, GcPsiC, GcTranSrc, PE)
  IF(iter .eq. 1) ResErr0 = ResErr
  IF(CMFDMaster) WRITE(mesg,'(6x,A9, I8, F17.6, 3x, F10.5, 1pe15.3)') 'CGOUTER', ItrCntl%GcCmfdit, eigv, ResErr/ResErr0, ResErr
  IF(CMFDMaster .AND. lLogOut) CALL message(io8, FALSE, TRUE, mesg)
  !!Convergence Check
  lExit = .TRUE.; lConv = .FALSE.
  !
  !SELECT CASE(ConvMod)
  !  CASE(1)
      IF((ResErr/ResErr0) .lt. conv_rescrit ) lconv = TRUE
  !  CASE(2)
  !    IF(EigErr .LT. Conv_eigvCrit) lconv = .TRUE.
  !END SELECT
  !
  lExit = lExit .AND. lConv
  IF(iter .LE. nitermin) lExit = .FALSE.
  !
  IF(lExit) EXIT   
  CONTINUE
ENDDO
IF(lGcFdk) CALL MgCmfdSolUpdt(CmInfo, GroupInfo, 0.5_8)

NULLIFY(PinXS, GcPinXS)
NULLIFY(PhiFm)
NULLIFY(GcPHiC, GcPsiC, GcPsiCD)
END SUBROUTINE
  
SUBROUTINE SetTranGcCmfdLinearSystem()
USE PARAM
USE TYPEDEF,     ONLY : CoreInfo_Type,         CmfdLs_Type,         PinXS_Type
USE CMFD_MOD,    ONLY : nxy,                  myzb,                myze,                &
                        myzbf,                myzef,                                                      &
                        hzfm,                 hz,                  SubPlaneMap,      SubPlaneRange,       &
                        PinNeighIdx,          PinVol,              PinVolFm
USE GcCMFD_mod,  ONLY : GcPinXs,              GcCmfdLs,            ngc
USE BASICOPERATION, ONLY : CP_CA
IMPLICIT NONE
INTEGER :: ig, iz, iz0, ixy, ibd
INTEGER, PARAMETER :: nbd = 4
REAL, POINTER :: Diag(:, :), RadOffDiag(:, :, :), AxOffDiag(:, :, :)
REAL :: Vol, neighVol, lmnt, offlmnt, dhat, dtil, alpha, dhmax

DO ig = 1, ngc
  Diag => GcCmfdLs(ig)%Diag
  DO iz = myzbf, myzef
    iz0 = SubPlaneMap(iz)
    DO ixy = 1, nxy
      Diag(ixy, iz) = Diag(ixy, iz) + GcPinXS(ixy, iz)%velo(ig) * PinVolFm(ixy, iz0)
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE SetTranCmfd2GLinearSystem()
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
DO iz = myzbf, myzef
  iz0 = SubPlaneMap(iz)
  DO ixy = 1, nxy
    Vol = PinVolFm(ixy, iz0)
    DO ig = 1, 2
      Diag(idx(ig, ig), ixy, iz) = Diag(idx(ig, ig), ixy, iz) + GcPinXS(ixy, iz)%velo(ig) * PinVolFm(ixy, iz0)
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE MakeGcKineticParam(Core, CmInfo, TranInfo, GroupInfo, GcGroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,         CmInfo_Type,            TranInfo_Type,         &
                             GroupInfo_Type,        PE_TYPE,                TranCntl_Type,         &
                             PinXS_Type
USE TRANCMFD_MOD,     ONLY : TrSrc,                 PrecSrc,                                       &
                             CmfdPrecSrcUpdt,       CmfdTranSrc
USE CNTL,             ONLY : nTracerCntl_Type
USE BasicOperation,   ONLY : CP_CA 
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo,     GcGroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE

TYPE(PinXS_TYPE), POINTER :: PinXS(:, :), GcPinXs(:, :)
REAL, POINTER :: Phi(:, :, :), GcPhi(:, :, :)
REAL, POINTER :: GcTranSrc(:, :, :)

!REAL, POINTER :: PrecSrc(:, :)

INTEGER, POINTER :: SubPlaneMap(:)
INTEGER :: ng, ngc  
INTEGER :: nxy, myzb, myze, myzbf, myzef

INTEGER :: igc, ig, iz, iz0, ixy

REAL :: Sum1, Sum2
REAL :: Velo, Delt, Theta, RVDELT, RVALPHA

nxy = Core%nxy; ng = GroupInfo%ng; ngc = GroupInfo%ngc
myzb = PE%myzb; myze = PE%myze; myzbf  = PE%myzbf; myzef  = PE%myzef

PinXS => CmInfo%PinXS; GcPinXs => CmInfo%GcPinXS; GcTranSrc => CmInfo%GcTranSrc
Phi => CmInfo%PhiFm; GcPhi => CmInfo%GcPhiC
SubPlaneMap => Core%SubPlaneMap
!Delt = TranCntl%Delt(TranCntl%nowstep)
!Theta = TranCntl%Theta
!!ALLOCATE(PrecSrc(1:nxy, myzbf:myzef))
!!CALL CmfdPrecSrcUpdt(PrecSrc, CmInfo%PrecFm, CmInfo%TranPsiFm, CmInfo%TranPsiFmd, TranCntl)
!CALL CP_CA(GcTranSrc(1:nxy, myzbf:myzef, 1:ngc), 0._8, nxy, myzef-myzbf+1, ngc)
!DO igc = 1, ngc
!  DO ig = GroupInfo%GCStruct(1, igc), GroupInfo%GCStruct(2, igc)
!    CALL CmfdTranSrc(TrSrc, CmInfo%PhiFm, CmInfo%TranPhiFm, CmInfo%PsiFm, PrecSrc, CmInfo%ResSrcFm, TranCntl, ig, PE%nCmfdThread)
!    DO iz = myzbf, myzef
!      DO ixy = 1, nxy
!        velo = PinXS(ixy, iz)%velo(ig)
!        RVDELT = 1._8 / (Velo * Delt * Theta) * Core%PinVolFm(ixy, iz)
!        sum1 = RVDELT * CmInfo%PhiFm(ixy, iz, ig) - TrSrc(ixy, iz)
!        GcTranSrc(ixy, iz, igc) = GcTranSrc(ixy, iz, igc) + TrSrc(ixy, iz)
!      ENDDO
!    ENDDO
!  ENDDO
!ENDDO
!
!DEALLOCATE(PrecSrc)
!CALL CmfdPrecSrcUpdt(PrecSrc, PrecFm, TranPsiFm, TranPsiFmd, TranCntl)
Delt = TranCntl%Delt(TranCntl%nowstep)
Theta = TranCntl%Theta
DO iz = myzbf, myzef
  iz0 = SubPlaneMap(iz)
  DO ixy = 1, nxy
    DO igc = 1, ngc
      Sum1 = 0; Sum2 = 0
      DO ig = GroupInfo%GCStruct(1, igc), GroupInfo%GCStruct(2, igc)
        velo = PinXS(ixy, iz0)%velo(ig)
        RVDELT = 1._8 / (Velo * Delt * Theta)
        RVALPHA = TranInfo%Expo_Alpha(ixy, iz, ig) / Velo 
        Sum1 = sum1 + (Rvdelt + RvAlpha) * Phi(ixy, iz, ig)
        Sum2 = Sum2 + Phi(ixy, iz, ig)
      ENDDO
      GcPinXS(ixy, iz)%Velo(igc) = Sum1 / Sum2
      CONTINUE
    ENDDO
  ENDDO
ENDDO

NULLIFY(PinXS, GcPinXS)
NULLIFY(GcTranSrc); NULLIFY(Phi, GcPhi); NULLIFY(SubPlaneMap)
END SUBROUTINE

SUBROUTINE UpdtGcTranSrc(Core, CmInfo, TranInfo, GroupInfo, GcGroupInfo, TranCntl, nTracerCntl, PE)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,         CmInfo_Type,            TranInfo_Type,         &
                             GroupInfo_Type,        PE_TYPE,                TranCntl_Type,         &
                             PinXS_Type
USE TRANCMFD_MOD,     ONLY : TrSrc,                                                                &
                             CmfdPrecSrcUpdt,       CmfdTranSrc
USE CNTL,             ONLY : nTracerCntl_Type
USE BasicOperation,   ONLY : CP_CA 
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(CmInfo_Type) :: CmInfo
TYPE(TranInfo_Type) :: TranInfo
TYPE(GroupInfo_Type) :: GroupInfo,     GcGroupInfo
TYPE(TranCntl_Type) :: TranCntl
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_TYPE) :: PE

TYPE(PinXS_TYPE), POINTER :: PinXS(:, :), GcPinXs(:, :)
REAL, POINTER :: Phi(:, :, :), GcPhi(:, :, :)
REAL, POINTER :: GcTranSrc(:, :, :)

REAL, POINTER :: PrecSrc(:, :)

INTEGER, POINTER :: SubPlaneMap(:)
INTEGER :: ng, ngc  
INTEGER :: nxy, myzb, myze, myzbf, myzef

INTEGER :: igc, ig, iz, iz0, ixy

REAL :: Sum1, Sum2
REAL :: Velo, Delt, Theta, RVDELT, RVALPHA

nxy = Core%nxy; ng = GroupInfo%ng; ngc = GroupInfo%ngc
myzb = PE%myzb; myze = PE%myze; myzbf  = PE%myzbf; myzef  = PE%myzef

PinXS => CmInfo%PinXS; GcPinXs => CmInfo%GcPinXS; GcTranSrc => CmInfo%GcTranSrc
Phi => CmInfo%PhiFm; GcPhi => CmInfo%GcPhiC
SubPlaneMap => Core%SubPlaneMap

ALLOCATE(PrecSrc(1:nxy, myzbf:myzef))
CALL CmfdPrecSrcUpdt(PrecSrc, CmInfo%PrecFm, CmInfo%TranPsiFm, CmInfo%TranPsiFmd, TranCntl)
CALL CP_CA(GcTranSrc(1:nxy, myzbf:myzef, 1:ngc), 0._8, nxy, myzef-myzbf+1, ngc)
DO igc = 1, ngc
  DO ig = GroupInfo%GCStruct(1, igc), GroupInfo%GCStruct(2, igc)
    CALL CmfdTranSrc(TrSrc, CmInfo%PhiFm, CmInfo%TranPhiFm, CmInfo%GcPsiC, PrecSrc, CmInfo%ResSrcFm, TranCntl, ig, PE%nCmfdThread)
    DO iz = myzbf, myzef
      DO ixy = 1, nxy
        GcTranSrc(ixy, iz, igc) = GcTranSrc(ixy, iz, igc) + TrSrc(ixy, iz)
      ENDDO
    ENDDO
  ENDDO
ENDDO
DEALLOCATE(PrecSrc)
END SUBROUTINE

FUNCTION TranGcResidualError(phi, psi, TranSrc, PE, constsrc)
USE PARAM
USE TYPEDEF,   ONLY : PE_TYPE
USE CMFD_MOD, ONLY : nxy,         myzbf,       myzef,               &
                     src,                                           &
                     SubPlaneMap, PinVol,      PinVolFm,            &
                     hzfm,        PinNeighIdx
USE GcCmfd_Mod, ONLY : ngc,       GcPinXS,     GcCmfdLs
#ifdef MPI_ENV
USE MPIComm_MOD, ONLY : REDUCE, GetNeighDat
#endif
IMPLICIT NONE
REAL :: TranGcResidualError
REAL, POINTER :: phi(:, :, :), psi(:, :), TranSrc(:, :, :)
REAL :: eigv
TYPE(PE_TYPE) :: PE
REAL, OPTIONAL :: ConstSrc
INTEGER :: ig, ig0, iz, iz0, ixy, ibd, ineigh, nbd
REAL :: vol, LMNT, tsrc, area
REAL :: myphi, neighphi, jsum, jnet, dtil, dhat
#ifdef MPI_ENV
REAL :: buf(2), buf0(2)
#endif
TranGcResidualError = 0; tsrc =0
nbd = 4
#define MatOpMod
!#define ResMG
#ifndef ResMG
#ifdef MPI_ENV
DO ig = 1, ngc
  IF(PE%nCMFDproc .GT. 1) THEN 
    CALL GetNeighDat(Phi(1:nxy, myzbf-1:myzef+1, ig), nxy, myzbf, myzef,   &
                           PE%myCmfdRank, PE%nproc, PE%MPI_CMFD_COMM)
  ENDIF
ENDDO
#endif
DO iz = myzbf, myzef
  iz0 = SubPlaneMap(iz)
  DO ixy = 1, nxy
    lmnt = 0
    DO ig = 1, ngc
      lmnt = lmnt + GcCMFDLS(ig)%Diag(ixy, iz) * Phi(ixy, iz, ig)
      DO ibd = 1, nbd
        ineigh = GcCMFDLS(ig)%NeighIdx(ibd, ixy)
        IF(ineigh .LE. 0) CYCLE
        lmnt = lmnt + GcCMFDLS(ig)%RadOffDiag(ibd, ixy, iz) * Phi(ineigh, iz, ig)        
      ENDDO
      !Axial Contribution
      lmnt = lmnt + GcCMFDLS(ig)%AxOffDiag(2, ixy, iz) * Phi(ixy, iz + 1, ig)
      lmnt = lmnt + GcCMFDLS(ig)%AxOffDiag(1, ixy, iz) * Phi(ixy, iz - 1, ig)
    ENDDO
    
    SRC(ixy, iz) = Psi(ixy, iz) / PinVolFm(ixy, iz0)
    DO ig = 1, ngc
      DO ig0 = GcPinXS(ixy, iz)%XSS(ig)%ib, GcPinXS(ixy, iz)%XSS(ig)%ie
        Src(ixy, iz) = Src(ixy, iz) + phi(ixy, iz, ig0) * GcPinXS(ixy, iz)%XSS(ig)%From(ig0)
      ENDDO
      Src(ixy, iz) = Src(ixy, iz) + TranSrc(ixy, iz, ig) / PinVolFm(ixy, iz0)
    ENDDO
    IF(PRESENT(ConstSrc)) SRC(ixy, iz) = SRC(ixy, iz) + ConstSrc * ngc
    SRC(ixy, iz) = Src(ixy, iz) * PinVolFm(ixy, iz0)
    lmnt = lmnt - SRC(ixy, iz) 
    tsrc = tsrc +  src(ixy, iz) * src(ixy, iz)
    TranGcResidualError = TranGcResidualError + lmnt * lmnt    
  ENDDO
ENDDO
#else

#endif
#ifdef MPI_ENV 
buf0  = (/TranGcResidualError, tsrc/)
CALL REDUCE(buf0, buf, 2, PE%MPI_CMFD_COMM, .TRUE.)
TranGcResidualError = buf(1); tsrc = buf(2)
#endif
TranGcResidualError = SQRT(TranGcResidualError/tsrc)

END FUNCTION
