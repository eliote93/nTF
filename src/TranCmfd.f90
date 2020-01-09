#include <defines.h>
SUBROUTINE SetTranCmfdLinearSystem(lCmfd, l3dim, Axsolver, TranCntl)
USE PARAM
USE TYPEDEF,        ONLY : CMFDLS_TYPE,      TranCntl_Type
USE CMFD_MOD,       ONLY : ng,        nxy,          myzb,         myze,                &
                           myzbf,     myzef,                                           &
                           hzfm,      hz,           SubPlaneMap,  SubPlaneRange,       &
                           CmfdPinXS, PinNeighIdx,  PinVol,       PinVolFm,            &
                           CmfdLs,    CmfdLs1g 
USE TRANCMFD_MOD,   ONLY : Expo,      Expo_Alpha  
USE BASICOPERATION, ONLY : CP_CA
IMPLICIT NONE
TYPE(TranCntl_Type) :: TranCntl
LOGICAL :: lcmfd, l3dim
INTEGER :: AxSolver
INTEGER :: ig, iz, iz0, ixy
INTEGER, PARAMETER :: nbd = 4
REAL, POINTER :: Diag(:, :)

REAL :: Delt, velo, RVDELT, rvalpha, Theta
REAL :: TrTerm
!LINEAR SYSTEM

Delt = TranCntl%Delt(TranCntl%nowstep)
Theta = TranCntl%Theta
CALL SetCmfdLinearSystem(lCmfd,  l3dim, AxSolver)

DO ig = 1, ng
  CMFDLS(ig)%AxialPlaneMap => SubPlaneMap
  Diag => CmfdLs(ig)%Diag
  DO iz = myzbf, myzef
    iz0 = SubPlaneMap(iz)
    DO ixy = 1, nxy
      velo = CmfdPinXS(ixy, iz0)%velo(ig)
      RVDELT = 1._8 / (Velo * Delt * Theta)
      !IF EXPONENTIAL TRANSPFORMATION IS NOT ACTIVE, THE ALPHA VALUE SHOULD BE ZERO
      RVALPHA = Expo_Alpha(ixy, iz, ig) / Velo 
      Diag(ixy, iz) = Diag(ixy, iz) + (Rvdelt + RVALPHA) * PinVolFm(ixy, iz0)
    ENDDO
  ENDDO
ENDDO
NULLIFY(Diag)
END SUBROUTINE

SUBROUTINE CmfdPrecUpdt(Prec, Psi, TranPsi, TranPsid, TranCntl)
USE PARAM
USE TYPEDEF,      ONLY : TranCntl_Type
USE CMFD_MOD,     ONLY : ng,            nxy,                                               &
                         myzb,          myze,            myzbf,           myzef,           &
                         CmfdPinXs,     SubPlaneMap,     PinVol,          PinVolFm
USE TRANCMFD_MOD, ONLY : nprec,         CellOmegam,      CellOmega0,      CellOmegap,      &
                         Lambda
IMPLICIT NONE
REAL, POINTER :: Prec(:, :, :), Psi(:, :), TranPsi(:, :), TranPsid(:, :)
TYPE(TranCntl_Type) :: TranCntl

INTEGER :: iz, iz0, ixy, iprec, nowstep
REAL :: DelT, reigv, kappa(nprec)

NowStep = TranCntl%NowStep
Delt = TranCntl%DelT(NowStep)

DO iprec = 1, nprec
  kappa(iprec) = exp(-delt * lambda(iprec))
ENDDO
!
!!precursor term
DO iz = myzbf, myzef
  iz0 = SubPlaneMap(iz)
  DO ixy = 1, nxy
    DO iprec = 1, nprec
      Prec(iprec, ixy, iz) = kappa(iprec) * Prec(iprec, ixy, iz)
    ENDDO
  ENDDO
ENDDO
!
!!Previous step Fission Source Contribution
DO iz = myzbf, myzef
  iz0 = SubPlaneMap(iz)
  DO ixy = 1, nxy
    DO iprec = 1, nprec
      Prec(iprec, ixy, iz) = Prec(iprec, ixy, iz) + CellOmegam(iprec, ixy, iz0) * TranPsid(ixy, iz)   &
                             + CellOmega0(iprec, ixy, iz0) * TranPsi(ixy, iz)                         &
                             + CellOmegap(iprec, ixy, iz0) * Psi(ixy, iz)
    ENDDO
  ENDDO
ENDDO
!
END SUBROUTINE

SUBROUTINE CmfdPrecSrcUpdt(PrecSrc, Prec, TranPsi, TranPsid, TranCntl)
USE PARAM
USE TYPEDEF,      ONLY : TranCntl_Type
USE CMFD_MOD,     ONLY : ng,            nxy,                                               &
                         myzb,          myze,            myzbf,           myzef,           &
                         CmfdPinXs,     SubPlaneMap,     PinVol,          PinVolFm
USE TRANCMFD_MOD, ONLY : nprec,         CellOmegam,      CellOmega0,      CellOmegap,      &
                         Lambda
IMPLICIT NONE
REAL, POINTER :: PrecSrc(:, :), Prec(:, :, :), TranPsi(:, :), TranPsid(:, :)
TYPE(TranCntl_Type) :: TranCntl

INTEGER :: iz, iz0, ixy, iprec, nowstep
REAL :: DelT, reigv, kappa(nprec)

NowStep = TranCntl%NowStep
Delt = TranCntl%DelT(NowStep)

DO iprec = 1, nprec
  kappa(iprec) = exp(-delt * lambda(iprec))
ENDDO

!precursor term
DO iz = myzbf, myzef
  iz0 = SubPlaneMap(iz)
  DO ixy = 1, nxy
    PrecSrc(ixy, iz) = 0
    DO iprec = 1, nprec
      PrecSrc(ixy, iz) = PrecSrc(ixy, iz) + lambda(iprec) * kappa(iprec) * Prec(iprec, ixy, iz)
    ENDDO
  ENDDO
ENDDO

!Previous step Fission Source Contribution
DO iz = myzbf, myzef
  iz0 = SubPlaneMap(iz)
  DO ixy = 1, nxy
    PrecSrc(ixy, iz) = PrecSrc(ixy, iz) + CellOmegam(0, ixy, iz0) * TranPsid(ixy, iz) +  CellOmega0(0, ixy, iz0) * TranPsi(ixy, iz)
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE CmfdSteadyPsiUpdt(phifm, psifm, nThread)
USE PARAM
USE CMFD_MOD, ONLY : ng,          nxy,          myzbf,       myzef, &
                     CmfdPinXs,   SubPlaneMap,  PinVol,      PinVolFm
IMPLICIT NONE
REAL, POINTER :: psifm(:, :), phifm(:, :, :)
INTEGER :: nThread

INTEGER :: ig, iz, iz0, ixy
!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(nThread) 
!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(ixy, ig)
DO iz = myzbf, myzef
  iz0 = SubPlaneMap(iz)
!$OMP DO
  DO ixy = 1, nxy
    Psifm(ixy, iz) = 0
    DO ig = 1, ng
      PsiFm(ixy, iz) = PsiFm(ixy, iz) + CmfdPinXs(ixy, iz0)%XSNF(ig)*PhiFm(ixy ,iz, ig)
    ENDDO
    Psifm(ixy, iz) = Psifm(ixy, iz)* PinVolFm(ixy, iz0)
  ENDDO
!$OMP END DO
ENDDO
!$OMP END PARALLEL
END SUBROUTINE

SUBROUTINE CmfdSteadySrcUpdt(SRC, psifm, phifm, eigv, ig, nThread)
USE PARAM
USE TYPEDEF, ONLY : scatmat
USE CMFD_MOD, ONLY : ng,        nxy,                                   &
                     myzb,      myze,         myzbf,       myzef,      &
                     CmfdPinXs, SubPlaneMap,  PinVol,      PinVolFm
USE OMP_LIB
IMPLICIT NONE
REAL, POINTER :: SRC(:, :), psifm(:, :), phifm(:, :, :)
REAL :: Eigv
TYPE(scatmat), POINTER :: XSS
REAL :: reigv, SS
INTEGER :: ig, nThread

INTEGER :: iz, iz0, ixy, ig2, igb, ige
reigv = one/eigv
!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(nThread) 
!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(ixy, XSS, SS, igb, ige, ig2)
DO iz = myzbf, myzef
  iz0 = SubPlaneMap(iz)
!$OMP DO
  DO ixy = 1, nxy
    SRC(ixy, iz) = psifm(ixy, iz) * reigv * CmfdPinXs(ixy, iz0)%Chi(ig)
    XSS =>CmfdPinXS(ixy, iz0)%Xss(ig) 
    igb = XSS%ib; ige = XSS%ie
    SS = 0
    DO ig2 = igb, ige
      SS = SS + phifm(ixy, iz, ig2) * XSS%From(ig2)
    ENDDO
    NULLIFY(XSS)
    SRC(ixy, iz) = SRC(ixy, iz) + SS * PinVolFm(ixy, iz0)
    !SRC(ixy, iz) = SRC(ixy, iz) 
  ENDDO
!$OMP END DO
ENDDO
!$OMP END PARALLEL
END SUBROUTINE

SUBROUTINE CmfdTranSrc(TranSrc, Phi, TranPhi, Psi, PrecSrc, ResSrc, TranCntl, ig, nThread)
USE PARAM
USE TYPEDEF,      ONLY : TranCntl_Type
USE CMFD_MOD,     ONLY : ng,            nxy,                                               &
                         myzb,          myze,            myzbf,           myzef,           &
                         CmfdPinXs,     SubPlaneMap,     PinVol,          PinVolFm
USE TRANCMFD_MOD, ONLY : Chid,          Expo,            Expo_Alpha
USE OMP_LIB
IMPLICIT NONE

REAL, POINTER :: TranSrc(:, :), PrecSrc(:, :), Psi(:, :), TranPhi(:, :, :), Phi(:, :, :), ResSrc(:, :, :)
TYPE(TranCntl_Type) :: TranCntl
INTEGER :: ig, nThread

INTEGER :: iz, iz0, ixy, iprec, nowstep
REAL :: DelT, rvdt, vol, Theta, thetah, RvAlpha
REAL :: chieff, chieff0, prevsrc

NowStep = TranCntl%NowStep
Delt = TranCntl%DelT(NowStep)
Theta = TranCntl%Theta
thetah = 1._8 / theta - 1
!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(nThread) 
!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(ixy, vol , rvdt, chieff0, RvAlpha, PrevSrc)
DO iz = myzbf, myzef
  iz0 = SubPlaneMap(iz)
!$OMP DO
  DO ixy = 1, nxy
    vol = PinvolFm(ixy, iz0)
    rvdt = 1._8 / (Delt * CmfdPinXs(ixy, iz)%velo(ig) * Theta)
    TranSrc(ixy, iz) = chid(ig) * PrecSrc(ixy, iz)
    chieff0 = chid(ig) * (CmfdPinXs(ixy, iz0)%omega - CmfdPinXs(ixy, iz0)%betat)
    TranSrc(ixy, iz) = TranSrc(ixy, iz) + chieff0 *Psi(ixy, iz)
    
    RvAlpha = Expo_Alpha(ixy, iz, ig) * Vol / CmfdPinXs(ixy, iz)%velo(ig)
    
    PrevSrc = Thetah * (ResSrc(ixy, iz, ig) - RvAlpha * TranPhi(ixy, iz, ig))
    PrevSrc = PrevSrc + rvdt * (TranPhi(ixy, iz, ig))*Vol
    PrevSrc = PrevSrc * Expo(ixy, iz, ig)
    TranSrc(ixy, iz) = TranSrc(ixy, iz) + PrevSrc
    CONTINUE
    !TranSrc(ixy, iz)=0
    !prevsrc = rvdt * TranPhi(ixy, iz, ig) - Expo_Alpha(ixy, iz, ig) * TranPhi(ixy, iz, ig) / CmfdPinXs(ixy, iz)%velo(ig)
    !TranSrc(ixy, iz) = TranSrc(ixy, iz) + rvdt * (TranPhi(ixy, iz, ig))*Vol
    !TranSrc(ixy, iz) = TranSrc(ixy, iz) + thetah * (ResSrc(ixy, iz, ig))
  ENDDO
!$OMP END DO
ENDDO
!$OMP END PARALLEL
END SUBROUTINE

SUBROUTINE CmfdTranSrc1g(TranSrc, Phi, TranPhi, Psi, PrecSrc, ResSrc, TranCntl)
USE PARAM
USE TYPEDEF,      ONLY : TranCntl_Type
USE CMFD_MOD,     ONLY : ng,            nxy,                                               &
                         myzb,          myze,            myzbf,           myzef,           &
                         CmfdPinXs,     SubPlaneMap,     PinVol,          PinVolFm
USE TRANCMFD_MOD, ONLY : Chid,                                                             &
                         Expo,          Expo_Alpha
IMPLICIT NONE

REAL, POINTER :: TranSrc(:, :), PrecSrc(:, :), Psi(:, :), TranPhi(:, :, :), Phi(:, :, :), ResSrc(:, :, :)
TYPE(TranCntl_Type) :: TranCntl



INTEGER :: iz, iz0, ixy, iprec, nowstep
REAL :: DelT, rvdt, vol, theta, thetah
REAL :: RValpha, PrevSrc
REAL :: chieff, chieff0
INTEGER :: ig

NowStep = TranCntl%NowStep
Delt = TranCntl%DelT(NowStep)
Theta = TranCntl%Theta
ThetaH = 1._8 / Theta - 1._8
DO iz = myzbf, myzef
  iz0 = SubPlaneMap(iz)
  DO ixy = 1, nxy
    !vol = PinvolFm(ixy, iz0)
    !rvdt = 1._8 / (Delt * CmfdPinXs(ixy, iz)%velo(ig))
    TranSrc(ixy, iz) = PrecSrc(ixy, iz)
    chieff0 = (CmfdPinXs(ixy, iz0)%omega - CmfdPinXs(ixy, iz0)%betat)
    TranSrc(ixy, iz) = TranSrc(ixy, iz) + chieff0 *Psi(ixy, iz)
    
    !TranSrc(ixy, iz) = TranSrc(ixy, iz) + rvdt * TranPhi(ixy, iz, ig)*Vol
    CONTINUE
  ENDDO
ENDDO

DO ig = 1, ng
  DO iz = myzbf, myzef
    iz0 = SubPlaneMap(iz)
    DO ixy = 1, nxy
      vol = PinvolFm(ixy, iz0)
      rvdt = 1._8 / (Delt * CmfdPinXs(ixy, iz)%velo(ig) * theta)
      
      RvAlpha = Expo_Alpha(ixy, iz, ig) * Vol / CmfdPinXs(ixy, iz)%velo(ig)

      PrevSrc = Thetah * (ResSrc(ixy, iz, ig) - RvAlpha * TranPhi(ixy, iz, ig))
      PrevSrc = PrevSrc + rvdt * (TranPhi(ixy, iz, ig))*Vol
      PrevSrc = PrevSrc * Expo(ixy, iz, ig)
      TranSrc(ixy, iz) = TranSrc(ixy, iz) + PrevSrc      
      !TranSrc(ixy, iz) = TranSrc(ixy, iz) + rvdt * (TranPhi(ixy, iz, ig))*Vol
      !TranSrc(ixy, iz) = TranSrc(ixy, iz) + ThetaH * ResSrc(ixy, iz, ig)
      CONTINUE
    ENDDO  
  ENDDO
ENDDO

END SUBROUTINE
!TranPhi, Psi, PrecSrc, TranCntl
FUNCTION TranResidualError(phifm, psifm, TranPhi, PrecSrc, ResSrc, TranCntl, PE)
USE PARAM
USE TYPEDEF,   ONLY : PE_TYPE, TranCntl_Type
USE CMFD_MOD,  ONLY : ng,          nxy,         myzbf,       myzef,  &
                     CmfdPinXs,   src,         CmfdLS,              &
                     SubPlaneMap, PinVol,      PinVolFm,            &
                     hzfm,        PinNeighIdx, CmfdSrcUpdt,         &
                     AddConstCmfdSrc 
USE TRANCMFD_MOD, ONLY : TrSrc,                                    &
                         CmfdTranSrc1g
USE OMP_LIB
#ifdef MPI_ENV
USE MPIComm_MOD,  ONLY : REDUCE, GetNeighDat
#endif
IMPLICIT NONE
REAL :: TranResidualError
REAL, POINTER :: phifm(:, :, :), psifm(:, :)
REAL, POINTER :: PrecSrc(:, :), TranPhi(:, :, :), ResSrc(:, :, :)
TYPE(PE_TYPE) :: PE
TYPE(TranCntl_Type) :: TranCntl

INTEGER :: ig, ig0, iz, iz0, ixy, ibd, ineigh, nbd, tid
INTEGER :: nThread
REAL :: vol, LMNT, tsrc, area
REAL :: myphi, neighphi, jsum, jnet, dtil, dhat
REAL :: tsrc_omp(nthreadmax), TranRes_omp(nthreadmax)
#ifdef MPI_ENV
REAL :: buf(2), buf0(2)
#endif

nbd = 4; tid = 1
nThread = PE%nCmfdThread

TranResidualError = 0; tsrc =0

tsrc_omp(1:nthread) = 0; TranRes_omp(1:nthread) = 0

#ifdef MPI_ENV
DO ig = 1, ng
  IF(PE%nCMFDproc .GT. 1) THEN 
    CALL GetNeighDat(PhiFm(1:nxy, myzbf-1:myzef+1, ig), nxy, myzbf, myzef,   &
                           PE%myCmfdRank, PE%nproc, PE%MPI_CMFD_COMM)
  ENDIF
ENDDO
#endif
CALL CmfdTranSrc1g(TrSrc, PhiFm, TranPhi, PsiFm, PrecSrc, ResSrc, TranCntl)

!$  call omp_set_dynamic(.FALSE.)
!$  call omp_set_num_threads(nThread) 
!$OMP PARALLEL DEFAULT(SHARED)      &
!$OMP PRIVATE(ixy, lmnt, ibd, ineigh, ig, ig0, tid)
!$  tid = omp_get_thread_num()+1
DO iz = myzbf, myzef
  iz0 = SubPlaneMap(iz)
!$OMP DO
  DO ixy = 1, nxy
    lmnt = 0
    DO ig = 1, ng
      lmnt = lmnt + CMFDLS(ig)%Diag(ixy, iz) * PhiFm(ixy, iz, ig)
      DO ibd = 1, nbd
        ineigh = CMFDLS(ig)%NeighIdx(ibd, ixy)
        IF(ineigh .LE. 0) CYCLE
        lmnt = lmnt + CMFDLS(ig)%RadOffDiag(ibd, ixy, iz0) * PhiFm(ineigh, iz, ig)        
      ENDDO
      !Axial Contribution
      lmnt = lmnt + CMFDLS(ig)%AxOffDiag(2, ixy, iz) * PhiFm(ixy, iz + 1, ig)
      lmnt = lmnt + CMFDLS(ig)%AxOffDiag(1, ixy, iz) * PhiFm(ixy, iz - 1, ig)
    ENDDO
    
    
    SRC(ixy, iz) = PsiFm(ixy, iz)/ PinVolFm(ixy, iz0) + TrSrc(ixy, iz) / PinVolFm(ixy, iz0)
    DO ig = 1, ng
      DO ig0 = CmfdPinXS(ixy, iz0)%XSS(ig)%ib, CmfdPinXS(ixy, iz0)%XSS(ig)%ie
        Src(ixy, iz) = Src(ixy, iz) + phifm(ixy, iz, ig0) * CmfdPinXS(ixy, iz0)%XSS(ig)%From(ig0)
      ENDDO
    ENDDO
    SRC(ixy, iz) = Src(ixy, iz) * PinVolFm(ixy, iz0)
    lmnt = lmnt - SRC(ixy, iz) 
    tsrc_omp(tid) = tsrc_omp(tid) +  src(ixy, iz) * src(ixy, iz)
    TranRes_omp(tid) = TranRes_omp(tid) + lmnt * lmnt    
  ENDDO
!$OMP END DO
  CONTINUE
ENDDO
!$OMP END PARALLEL
tsrc = sum(tsrc_omp(1:nthread))
TranResidualError = sum(TranRes_omp(1:nthread))
#ifdef MPI_ENV 
buf0  = (/TranResidualError, tsrc/)
CALL REDUCE(buf0, buf, 2, PE%MPI_CMFD_COMM, .TRUE.)
TranResidualError = buf(1); tsrc = buf(2)
#endif
TranResidualError = SQRT(TranResidualError/tsrc)

END FUNCTION

FUNCTION TranReactivityUpdt(CmInfo, eigv, TranCntl, PE)
USE PARAM
USE TYPEDEF,   ONLY : CmInfo_Type, PE_TYPE, TranCntl_Type
USE CMFD_MOD,  ONLY : ng,          nxy,         myzbf,       myzef,  &
                     CmfdPinXs,    src,         CmfdLS,              &
                     AxDhat,       AxDtil,                           &
                     SubPlaneMap, PinVol,      PinVolFm,            &
                     hzfm,        PinNeighIdx, CmfdSrcUpdt,         &
                     AddConstCmfdSrc 
USE TRANCMFD_MOD, ONLY : TrSrc,   nprec,      lambda,                &
                         CmfdTranSrc1g
#ifdef MPI_ENV
USE MPIComm_MOD,  ONLY : REDUCE, GetNeighDat
#endif
IMPLICIT NONE

TYPE(CmInfo_Type) :: CmInfo
REAL :: eigv
TYPE(PE_TYPE) :: PE
TYPE(TranCntl_Type) :: TranCntl


REAL, POINTER :: phifm(:, :, :), psifm(:, :), PrecFm(:, :, :)
REAL :: TranReactivityUpdt
REAL :: vol, area, jsum
REAL :: myphi, neighphi, jnet, Dtil, Dhat, loss, lmnt
REAL :: resphi, psil1, psil2, betaavg
INTEGER :: ig, ig0, iz, iz0, ixy, ibd, ineigh, nbd

INTEGER :: COMM
REAL :: BUF0(4), BUF(4)

PhiFm => CmInfo%PhiFm; PsiFm => CmInfo%PsiFm; PrecFm => CmInfo%PrecFm

CALL SetCmfdLinearSystem(TRUE, TRUE, 1)

COMM = PE%MPI_CMFD_COMM

resphi = 0
psil1 = 0; psil2 = 0
nbd =4 
DO iz = myzbf, myzef
  iz0 = SUbPlaneMap(iz)
  DO ixy = 1, nxy  
    vol = PinVolFm(ixy, iz0)
    area = vol / hzfm(iz)
    lmnt = 0; loss = 0
    DO ig = 1, ng
      jsum  = 0
      myphi = PhiFm(ixy, iz, ig)
     !Radial Current
      DO ibd = 1, nbd
        ineigh = PinNeighIdx(ibd, ixy)
        dtil = CMfdPinXs(ixy, iz0)%Dtil(ibd, ig)
        dhat = CMfdPinXs(ixy, iz0)%Dhat(ibd, ig)
        neighphi = 0
        IF(ineigh .GT. 0) neighphi = PhiFm(ineigh, iz, ig)
        jnet = (dtil - dhat)*myphi  -(dtil + dhat)*neighphi
        jsum = jsum + jnet
      ENDDO    
      jsum = jsum * hzfm(iz)
      Dtil = AxDtil(1, ixy, iz, ig); Dhat = AxDhat(1, ixy, iz, ig)
      !Dtil = AxFlx(iz, ixy)%Dtil(1, ig); Dhat = AxFlx(iz, ixy)%Dhat(1, ig)
      neighphi = PhiFM(ixy, iz - 1, ig); 
      jnet = (dtil - dhat)*myphi  -(dtil + dhat)*neighphi
      jsum = jsum + jnet * area        
      Dtil = AxDtil(2, ixy, iz, ig); Dhat = AxDhat(2, ixy, iz, ig)
      !Dtil = AxFlx(iz, ixy)%Dtil(2, ig); Dhat = AxFlx(iz, ixy)%Dhat(2, ig)
      neighphi = PhiFM(ixy, iz + 1, ig); 
      jnet = (dtil - dhat)*myphi  -(dtil + dhat)*neighphi
      jsum = jsum + jnet * area        
      loss = loss + jsum + vol*CmfdPinXs(ixy, iz0)%xsr(ig)*myphi
    ENDDO
    lmnt = loss
    SRC(ixy, iz) = 0
    DO ig = 1, ng
      DO ig0 = CmfdPinXS(ixy, iz0)%XSS(ig)%ib, CmfdPinXS(ixy, iz0)%XSS(ig)%ie
        Src(ixy, iz) = Src(ixy, iz) + phifm(ixy, iz, ig0) * CmfdPinXS(ixy, iz0)%XSS(ig)%From(ig0)
      ENDDO
    ENDDO      
    SRC(ixy, iz) =  SRC(ixy, iz) * vol + (1._8- CmfdPinXs(ixy, iz0)%betat) * PsiFm(ixy, iz) 
    DO ig = 1, nprec
      SRC(ixy, iz) = SRC(ixy, iz) + lambda(ig) * PrecFm(ig, ixy, iz)
    ENDDO
    resphi = resphi + (SRC(ixy, iz) - Loss)*PsiFm(ixy, iz)
    psil2 = psil2 + PsiFm(ixy, iz)*PsiFm(ixy, iz)
    psil1 = psil1 + PsiFm(ixy, iz)
    betaavg = betaavg + CmfdPinXs(ixy, iz0)%betat * PsiFm(ixy, iz)
  ENDDO
ENDDO
#ifdef MPI_ENV
buf0 = (/resphi, psil2, psil1, betaavg/)
CALL REDUCE(buf0, buf, 4, COMM, .TRUE.)
resphi = buf(1); psil2 = buf(2); psil1 = buf(3); betaavg = buf(4)
#endif
resphi = resphi / psil2
betaavg = betaavg / psil1
TranReactivityUpdt = resphi / betaavg / Eigv   
END FUNCTION