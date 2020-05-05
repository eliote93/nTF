#include <defines.h>
SUBROUTINE AllocCMFD2G()
USE PARAM
USE TYPEDEF,         ONLY : CmInfo_Type,     CmfdLs_Type,         PE_TYPE
USE Geom,            ONLY : Core,            ncbd
USE Core_mod,        ONLY : GcPinXs,          GCCMFDLS,           GcPsiC,         GcPsiCD,    &
                            GcPhiC,           GroupInfo,          GcGroupInfo,    CmInfo
USE CMFD_Mod,        ONLY : PinNeighIdx,      AllocPinXs
USE CMFD2G_Mod,      ONLY : PhiC2g,           src2g
USE BiCGSTAB2G_mod,  ONLY : AllocBiCGSTAB2g
USE CNTL,            ONLY : nTracerCntl
USE PE_MOD,          ONLY : PE
USE Allocs
IMPLICIT NONE
INTEGER :: nxy, myzbf, myzef, myzb, myze, ngc
INTEGER :: ig, iz ,ixy

nxy = Core%nxy
myzb = PE%myzb; myze = PE%myze
myzbf = PE%myzbf; myzef = PE%myzef 
nGC = nTracerCntl%nGC

CALL Dmalloc0(GcPsiC, 1, nxy, myzbf, myzef)
CALL Dmalloc0(GcPsiCD, 1, nxy, myzbf, myzef)
CALL Dmalloc0(GcPhiC, 1, nxy, myzbf-1, myzef+1, 1, ngc)

CALL Dmalloc0(PhiC2g, 1, ngc, 1, nxy, myzbf-1, myzef+1)
CALL Dmalloc0(src2g, 1, ngc, 1, nxy, myzbf-1, myzef+1)

ALLOCATE(GcPinXS(nxy,myzbf-1:myzef+1))

ALLOCATE(GCCMFDLS(1))

CALL Dmalloc0(GCCMFDLS(1)%diag2g, 1, 4, 1, nxy, myzbf, myzef)
CALL Dmalloc0(GCCMFDLS(1)%RadOffDiag2g, 1, 4, 1, ncbd, 1, nxy, myzbf, myzef)
CALL Dmalloc0(GCCMFDLS(1)%AxOffDiag2g, 1, 4, 1, 2, 1, nxy, myzbf, myzef)
GCCMFDLS(1)%NeighIdx => PinNeighIdx
GCCMFDLS(1)%myzbf = myzbf; GCCMFDLS%myzef = myzef
GCCMFDLS(1)%nxy = nxy; GCCMFDLS%nbd = ncbd
GCCMFDLS(1)%l2G = .TRUE.
CALL Dmalloc0(GCCMFDLS(1)%AxialPlaneMap, myzbf, myzef)
DO iz = myzbf, myzef
  GCCMFDLS(1)%AxialPlaneMap(iz) = iz
ENDDO

CALL SetCmfdMpiOmpEnv(Core, GCCMFDLS, 1, PE)


DO iz= myzbf, myzef
  DO ixy = 1, nxy
    CALL AllocPinXS(GcPinXs(ixy, iz), ngc, ncbd, GcGroupInfo%InScatRange)
  ENDDO
ENDDO
#ifdef MPI_ENV
  DO ixy = 1, nxy
    CALL Dmalloc(GcPinXS(ixy, myzef+1)%XSD, ngc)
    CALL Dmalloc(GcPinXS(ixy, myzbf-1)%XSD, ngc)
  ENDDO
#endif
CmInfo%GcPinXS => GcPinXS; CmInfo%GCCMFDLS => GCCMFDLS
CmInfo%GcPhiC => GcPhiC; 
CmInfo%GcPsiC => GcPsiC; CmInfo%GcPsiCD => GcPsiCD

CALL AllocBiCGSTAB2g(nxy, myzbf, myzef)

END SUBROUTINE

SUBROUTINE Cmfd2GSrcUpdt(SRC, psi, reigv)
USE PARAM
USE CMFD_MOD,   ONLY : nxy,                                   &
                       myzb,        myze,         myzbf,       myzef,      &
                       SubPlaneMap, PinVolFm
USE GcCMFD_MOD, ONLY : GcPinXS
IMPLICIT NONE
REAL, POINTER :: SRC(:, :, :), psi(:, :)
REAL :: reigv

INTEGER :: iz, iz0, ixy

DO iz = myzbf, myzef
  iz0 = SubPlaneMap(iz)
  DO ixy = 1, nxy
    SRC(1, ixy, iz) = psi(ixy, iz) * reigv * GcPinXs(ixy, iz)%Chi(1)
    SRC(2, ixy, iz) = psi(ixy, iz) * reigv * GcPinXs(ixy, iz)%Chi(2)
  ENDDO
ENDDO

END SUBROUTINE


SUBROUTINE ConvertArray2G(PhiC2G, GcPhiC, nxy, iz1, iz2, imod)
USE PARAM
IMPLICIT NONE
REAL :: PhiC2G(2, nxy, iz1:iz2), GcPhiC(nxy, iz1:iz2, 2)
INTEGER :: iz1, iz2, nxy, imod
INTEGER :: iz, ixy
IF(imod .EQ. 1) THEN
  DO iz = iz1, iz2
    DO ixy = 1, nxy
      PhiC2G(1, ixy, iz) = GcPhiC(ixy, iz, 1)
      PhiC2G(2, ixy, iz) = GcPhiC(ixy, iz, 2)
    ENDDO
  ENDDO
ELSE
  DO iz = iz1, iz2
    DO ixy = 1, nxy
      GcPhiC(ixy, iz, 1) = PhiC2G(1, ixy, iz)
      GcPhiC(ixy, iz, 2) = PhiC2G(2, ixy, iz)
    ENDDO
  ENDDO
ENDIF
END SUBROUTINE

SUBROUTINE WielandtUpdt(CMFD2GLS, PsiC, PsiCd, Eigv, Eigvs, PsiErr, iter, PE)
USE PARAM
USE TYPEDEF,    ONLY : PE_TYPE,      CMFDLS_TYPE
USE CMFD_MOD,   ONLY : nxy,         myzbf,       myzef,               &
                       SubPlaneMap, PinVol,      PinVolFm,            &
                       hzfm,        PinNeighIdx
USE GcCmfd_Mod, ONLY : ngc,       GcPinXS,                            &
                       GcCmfdEigUpdt
USE MpiCOmm_mod, ONLY : REDUCE
IMPLICIT NONE
TYPE(CmfdLS_TYPE) :: CMFD2GLS
REAL, POINTER :: PsiC(:, :), PsiCd(:, :)
REAL :: Eigv, Eigvs, PsiErr
REAL :: buf0(3), buf(3)
INTEGER :: iter
TYPE(PE_TYPE) :: PE

INTEGER :: iz, ixy
REAL :: Gamma, Gammad, Gamman
REAL :: reigv, reigvs, reigvd, reigvsd, reigvsdel, reigvdel
REAL, PARAMETER :: Eshift = 0.2
REAL, PARAMETER :: Eshift0 = 0.5
IF(iter .EQ. 1) THEN
  CALL GcCmfdEigUpdt(PsiC, PsiCD, Eigv, psierr, PE)
  Eigvs = Eigv + Eshift
  IF(psierr .LT. 1.0E-5) Eigvs = Eigv + Eshift0
  reigvs = 1 / eigvs; reigvsdel = reigvs
  CALL WielandtLsShift(Cmfd2GLs, reigv, reigvsdel, PE)
  RETURN
ENDIF 

gammad = 0; gamman = 0; psierr=0
DO iz = myzbf, myzef
  DO ixy = 1, nxy 
    gamman = gamman + PsiC(ixy, iz) * PsiC(ixy, iz)
    gammad = gammad + PsiC(ixy, iz) * PsiCd(ixy, iz)
    psierr = psierr + (PsiC(ixy, iz) - PsiCd(ixy, iz)) ** 2
  ENDDO
ENDDO

#ifdef MPI_ENV
buf0=(/gamman, gammad, psierr/)
CALL REDUCE(buf0, buf, 3, PE%MPI_CMFD_COMM, .TRUE.)
gamman = buf(1); gammad = buf(2); psierr = buf(3)
#endif

reigv = 1._8/eigv; reigvs = 1._8/eigvs
reigvd = reigv; reigvsd = reigvs

gamma=gammad/gamman
eigv = 1._8 / (reigv*gamma+(1-gamma)*reigvs)
psierr= sqrt(psierr/gamman)

Eigvs = Eigv + Eshift
IF(psierr .LT. 1.0E-8) Eigvs = Eigv + Eshift0

reigv = 1._8 / eigv; reigvs = 1._8 / eigvs
reigvsdel = reigvs - reigvsd; reigvdel = reigv - reigvs
CALL SetCmfd2gLinearSystem(.TRUE., .TRUE.)  
CALL WielandtLsShift(Cmfd2GLs, reigv, reigvs, PE)
!CALL WielandtLsShift(Cmfd2GLs, reigv, reigvsdel, PE)
END SUBROUTINE

SUBROUTINE WielandtLsShift(Cmfd2GLs, reigv, reigvsdel, PE)
USE PARAM
USE TYPEDEF,    ONLY : PE_TYPE,      CMFDLS_TYPE
USE CMFD_MOD,   ONLY : nxy,         myzbf,       myzef,               &
                       SubPlaneMap, PinVol,      PinVolFm,            &
                       hzfm,        PinNeighIdx
USE GcCmfd_Mod, ONLY : ngc,       GcPinXS
IMPLICIT NONE
TYPE(CMFDLS_TYPE) :: Cmfd2Gls
REAL :: reigv, reigvsdel

REAL, POINTER :: Diag2g(:, :, :)
TYPE(PE_TYPE) :: PE
REAL :: vol, xsnf(2), chi(2)
INTEGER :: iz, ixy, idx0, iz0
INTEGER, PARAMETER :: IDX(2,2)= RESHAPE((/1, 3, 2, 4/), SHAPE(IDX))

Diag2g => Cmfd2GLs%Diag2G
DO iz = myzbf, myzef
  iz0 = SubPlaneMap(iz)
  DO ixy = 1, nxy
    vol = PinVolFm(ixy, iz0) 
    xsnf = GcPinXS(ixy, iz)%xsnf * vol
    chi = GcPinXS(ixy, iz)%chi
    idx0=IDX(1,1)
    Diag2g(idx0, ixy, iz) = Diag2g(idx0, ixy, iz) - xsnf(1)*chi(1)* reigvsdel
    idx0=IDX(1,2)
    Diag2g(idx0, ixy, iz) = Diag2g(idx0, ixy, iz) - xsnf(2)*chi(1)* reigvsdel
    idx0=IDX(2,2)
    Diag2g(idx0, ixy, iz) = Diag2g(idx0, ixy, iz) - xsnf(2)*chi(2)* reigvsdel
    idx0=IDX(2,1)
    Diag2g(idx0, ixy, iz) = Diag2g(idx0, ixy, iz) - xsnf(1)*chi(2)* reigvsdel
  ENDDO
ENDDO
NULLIFY(Diag2g)
END SUBROUTINE

FUNCTION ResidualError2G(Phi2G, psi, reigv, reigvs, PE, constsrc)
USE PARAM
USE TYPEDEF,   ONLY : PE_TYPE
USE CMFD_MOD, ONLY : nxy,         myzbf,       myzef,               &
                     SubPlaneMap, PinVol,      PinVolFm,            &
                     hzfm,        PinNeighIdx
USE GcCmfd_Mod, ONLY : ngc,       GcPinXS,     GcCmfdLs
USE MAT2x2OP
#ifdef MPI_ENV
USE MPIComm_MOD, ONLY : REDUCE, GetNeighDat
#endif
IMPLICIT NONE
REAL :: ResidualError2G
REAL, POINTER :: phi2G(:, :, :), psi(:, :)
REAL :: reigv, reigvs
TYPE(PE_TYPE) :: PE
REAL, OPTIONAL :: ConstSrc
INTEGER, PARAMETER :: IDX(2,2)= RESHAPE((/1, 3, 2, 4/), SHAPE(IDX))
INTEGER :: ig, ig0, iz, iz0, ixy, ibd, ineigh, nbd, idx0
REAL :: vol, LMNT(2), Diag2g(4), lmnt0, tsrc, area, src, xsnf(2), chi(2)
REAL :: myphi, neighphi, jsum, jnet, dtil, dhat
#ifdef MPI_ENV
REAL :: buf(2), buf0(2)
#endif
ResidualError2G = 0; tsrc =0
nbd = 4
#define MatOpMod
!#define ResMG
#ifndef ResMG
#ifdef MPI_ENV
IF(PE%nCMFDproc .GT. 1) THEN 
  CALL GetNeighDat(Phi2g(1:2, 1:nxy, myzbf-1:myzef+1), 2, nxy, myzbf, myzef,   &
                         PE%myCmfdRank, PE%nproc, PE%MPI_CMFD_COMM)
ENDIF
#endif
DO iz = myzbf, myzef
  iz0 = SubPlaneMap(iz)
  DO ixy = 1, nxy
    vol = PinVolFm(ixy, iz0) 
    xsnf = GcPinXS(ixy, iz)%xsnf * vol
    chi = GcPinXS(ixy, iz)%chi
  
    lmnt = 0
    Diag2g = GcCMFDLS(1)%Diag2g(:, ixy, iz)
    
    idx0=IDX(1,1)
    Diag2g(idx0) = Diag2g(idx0) + xsnf(1)*chi(1)* reigvs
    idx0=IDX(1,2)
    Diag2g(idx0) = Diag2g(idx0) + xsnf(2)*chi(1)* reigvs
    idx0=IDX(2,2)
    Diag2g(idx0) = Diag2g(idx0) + xsnf(2)*chi(2)* reigvs
    idx0=IDX(2,1)
    Diag2g(idx0) = Diag2g(idx0) + xsnf(1)*chi(2)* reigvs    
    
    lmnt = SUBMATVECOP(Diag2g, Phi2G(:, ixy, iz))
    
    DO ibd = 1, nbd
      ineigh = GcCMFDLS(1)%NeighIdx(ibd, ixy)
      IF(ineigh .LE. 0) CYCLE
      lmnt = lmnt +SUBMATVECOP(GcCMFDLS(1)%RadOffDiag2g(:, ibd, ixy, iz), Phi2g(:, ineigh, iz))        
    ENDDO
    !Axial Contribution
    lmnt = lmnt + SUBMATVECOP(GcCMFDLS(1)%AxOffDiag2g(:, 2, ixy, iz), Phi2G(:, ixy, iz + 1))
    lmnt = lmnt + SUBMATVECOP(GcCMFDLS(1)%AxOffDiag2g(:, 1, ixy, iz), Phi2G(:, ixy, iz - 1))
!   
    lmnt0 = sum(lmnt)
    SRC = Psi(ixy, iz) *reigv 

    IF(PRESENT(ConstSrc)) SRC = SRC + ConstSrc * ngc
    SRC = Src !* PinVolFm(ixy, iz0)
    lmnt0 = lmnt0 - SRC
    tsrc = tsrc +  src * src
    ResidualError2G = ResidualError2G + lmnt0 * lmnt0   
  ENDDO
ENDDO
#else

#endif
#ifdef MPI_ENV 
buf0  = (/ResidualError2G, tsrc/)
CALL REDUCE(buf0, buf, 2, PE%MPI_CMFD_COMM, .TRUE.)
ResidualError2G = buf(1); tsrc = buf(2)
#endif
ResidualError2G = SQRT(ResidualError2g/tsrc)

END FUNCTION

SUBROUTINE NegativeFixUp2G(PhiC2G, PsiC2G, PE)
USE PARAM
USE TYPEDEF,   ONLY : PE_TYPE
USE CMFD_MOD,  ONLY : nxy,         myzbf,       myzef,   nzfm
USE IOUTIL,    ONLY : message
USE FILES,     ONLY : io8
#ifdef MPI_ENV
USE MPICOMM_MOD, ONLY : REDUCE
#endif
IMPLICIT NONE
REAL, POINTER :: PhiC2G(:, :, :), PsiC2G(:, :)
TYPE(PE_TYPE) :: PE

INTEGER :: ixy, iz, ig, k, m
REAL :: crit
#ifdef MPI_ENV
INTEGER :: COMM
INTEGER :: buf(2), buf0(2)
#endif

#ifdef MPI_ENV
COMM = PE%MPI_CMFD_COMM
#endif
k = 0; m = 0
DO iz = myzbf, myzef
  DO ixy = 1, nxy
    DO ig =1, 2
      m = m + 1
      IF(PhiC2G(ig, ixy, iz) .GT. 0._8) CYCLE
      k = k + 1 
    ENDDO
  ENDDO
ENDDO
#ifdef MPI_ENV
buf0 = (/k, m/)
CALL REDUCE(buf0, buf, 2, COMM, .TRUE.)
k = buf(1); m = buf(2)
#endif
crit = 0.8 * dble(m)
IF(DBLE(k) .GT. crit) THEN
  IF(PE%CMFDMASTER) CALL MESSAGE(io8, TRUE, TRUE, 'Neagtive Flux Fix UP')
  DO iz = myzbf, myzef
    DO ixy = 1, nxy
      PsiC2G(ixy, iz) = - PSIC2G(ixy, iz)
      PhiC2G(:, ixy, iz) = - PhiC2G(:, ixy, iz)
    ENDDO
  ENDDO
ENDIF 
END SUBROUTINE