#include <defines.h>
SUBROUTINE CmfdSrcUpdt(SRC, psifm, phifm, eigv, ig)
USE PARAM
USE TYPEDEF, ONLY : scatmat
USE TIMER,   only : nTracer_dclock, TimeChk
USE CMFD_MOD, ONLY : ng,        nxy,                                   &
                     myzb,      myze,         myzbf,       myzef,      &
                     CmfdPinXs, SubPlaneMap,  PinVol,      PinVolFm
IMPLICIT NONE
REAL, POINTER :: SRC(:, :), psifm(:, :), phifm(:, :, :)
REAL :: Eigv
TYPE(scatmat), POINTER :: XSS
REAL :: reigv, SS
REAL :: Tbeg, Tend
INTEGER :: ig

INTEGER :: iz, iz0, ixy, ig2, igb, ige
reigv = one/eigv

TBeg = nTracer_dclock(FALSE, FALSE)

DO iz = myzbf, myzef
  iz0 = SubPlaneMap(iz)
  DO ixy = 1, nxy
    SRC(ixy, iz) = psifm(ixy, iz) * reigv * CmfdPinXs(ixy, iz0)%Chi(ig)
    XSS =>CmfdPinXS(ixy, iz0)%Xss(ig) 
    igb = XSS%ib; ige = XSS%ie
    SS = 0
    DO ig2 = igb, ige
      SS = SS + phifm(ixy, iz, ig2) * XSS%From(ig2)
    ENDDO
    SRC(ixy, iz) = SRC(ixy, iz) + SS * PinVolFm(ixy, iz0)
    !SRC(ixy, iz) = SRC(ixy, iz) 
  ENDDO
ENDDO
NULLIFY(XSS)

Tend = nTracer_dclock(FALSE, FALSE)
TimeChk%AxBTime = TimeChk%AxBTime + (Tend - Tbeg)

END SUBROUTINE




SUBROUTINE CmfdPsiUpdt(phifm, psifm)
USE PARAM
USE CMFD_MOD, ONLY : ng,          nxy,          myzbf,       myzef, &
                     CmfdPinXs,   SubPlaneMap,  PinVol,      PinVolFm
IMPLICIT NONE
REAL, POINTER :: psifm(:, :), phifm(:, :, :)
INTEGER :: ig, iz, iz0, ixy
DO iz = myzbf, myzef
  iz0 = SubPlaneMap(iz)
  DO ixy = 1, nxy
    Psifm(ixy, iz) = 0
    DO ig = 1, ng
      PsiFm(ixy, iz) = PsiFm(ixy, iz) + CmfdPinXs(ixy, iz0)%XSNF(ig)*PhiFm(ixy ,iz, ig)
    ENDDO
    Psifm(ixy, iz) = Psifm(ixy, iz)* PinVolFm(ixy, iz0)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE CmfdEigUpdate(psifm, psifmd, eigv, psierr, PE)
USE PARAM
USE TYPEDEF,  ONLY : PE_TYPE
USE CMFD_MOD, ONLY : ng,          nxy,         myzbf,       myzef, &
                     CmfdPinXs,  SubPlaneMap,  PinVol,      PinVolFm
#ifdef MPI_ENV
USE MPIComm_mod, ONLY : REDUCE
#endif
IMPLICIT NONE
TYPE(PE_TYPE), OPTIONAL :: PE
REAL, POINTER :: psifm(:, :), psifmd(:, :)
REAL :: psierr, eigv
REAL :: psipsi, psipsid
INTEGER :: ig, iz, iz0, ixy
#ifdef MPI_ENV
REAL :: buf(3), buf0(3) 
#endif
psipsi = zero; psipsid = zero; psierr = zero
DO iz = myzbf, myzef
  DO ixy = 1, nxy
    psipsi = psipsi + Psifm(ixy, iz) * PsiFm(ixy, iz)
    psipsid = psipsid + Psifm(ixy, iz) * PsiFmD(ixy, iz)
    psierr = psierr + (Psifm(ixy, iz) - PsiFmD(ixy, iz))**2
  ENDDO
ENDDO     
#ifdef MPI_ENV
IF(Present(PE)) THEN
  buf0 = (/psipsi, psipsid, psierr/)
  CALL REDUCE(buf0, buf, 3, PE%MPI_CMFD_COMM, .TRUE.)
  psipsi = buf(1); psipsid = buf(2); psierr = buf(3)
ENDIF
#endif
eigv = eigv * psipsi/psipsid              
psierr = psierr / psipsi
psierr = SQRT(psierr)
END SUBROUTINE



FUNCTION ResidualError(phifm, psifm, eigv, ig1, ig2, PE, constsrc)
USE PARAM
USE TYPEDEF,   ONLY : PE_TYPE
USE CMFD_MOD, ONLY : ng,          nxy,         myzbf,       myzef,  &
                     CmfdPinXs,   src,         CmfdLS,              &
                     SubPlaneMap, PinVol,      PinVolFm,            &
                     hzfm,        PinNeighIdx, CmfdSrcUpdt,         &
                     AddConstCmfdSrc                
#ifdef MPI_ENV
USE MPIComm_MOD, ONLY : REDUCE, GetNeighDat
#endif
IMPLICIT NONE
REAL :: ResidualError
REAL, POINTER :: phifm(:, :, :), psifm(:, :)
REAL :: eigv
INTEGER :: ig1, ig2
TYPE(PE_TYPE) :: PE
REAL, OPTIONAL :: ConstSrc
INTEGER :: ig, ig0, iz, iz0, ixy, ibd, ineigh, nbd
REAL :: vol, LMNT, tsrc, area
REAL :: myphi, neighphi, jsum, jnet, dtil, dhat
#ifdef MPI_ENV
REAL :: buf(2), buf0(2)
#endif
ResidualError = 0; tsrc =0
nbd = 4
#define MatOpMod
!#define ResMG
#ifndef ResMG
#ifdef MPI_ENV
DO ig = 1, ng
  IF(PE%nCMFDproc .GT. 1) THEN 
    CALL GetNeighDat(PhiFm(1:nxy, myzbf-1:myzef+1, ig), nxy, myzbf, myzef,   &
                           PE%myCmfdRank, PE%nproc, PE%MPI_CMFD_COMM)
  ENDIF
ENDDO
#endif
DO iz = myzbf, myzef
  iz0 = SubPlaneMap(iz)
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
    
    SRC(ixy, iz) = PsiFm(ixy, iz) / Eigv / PinVolFm(ixy, iz0)
    DO ig = 1, ng
      DO ig0 = CmfdPinXS(ixy, iz0)%XSS(ig)%ib, CmfdPinXS(ixy, iz0)%XSS(ig)%ie
        Src(ixy, iz) = Src(ixy, iz) + phifm(ixy, iz, ig0) * CmfdPinXS(ixy, iz0)%XSS(ig)%From(ig0)
      ENDDO
    ENDDO
    IF(PRESENT(ConstSrc)) SRC(ixy, iz) = SRC(ixy, iz) + ConstSrc * ng
    SRC(ixy, iz) = Src(ixy, iz) * PinVolFm(ixy, iz0)
    lmnt = lmnt - SRC(ixy, iz) 
    tsrc = tsrc +  src(ixy, iz) * src(ixy, iz)
    ResidualError = ResidualError + lmnt * lmnt    
  ENDDO
  CONTINUE
ENDDO
#else
DO ig = ig1, ig2
  CALL CmfdSrcUpdt(Src, Psifm, phifm, eigv, ig)
#ifdef MPI_ENV
  IF(PE%nCMFDproc .GT. 1) THEN 
    CALL GetNeighDat(PhiFm(1:nxy, myzbf-1:myzef+1, ig), nxy, myzbf, myzef,   &
                           PE%myCmfdRank, PE%nproc, PE%MPI_CMFD_COMM)
  ENDIF
#endif
  IF(PRESENT(ConstSrc)) CALL AddConstCmfdSrc(Src, ConstSrc)
  DO iz = myzbf, myzef
    iz0 = SUbPlaneMap(iz)
    DO ixy = 1, nxy
#ifdef MatOpMod
      lmnt = CMFDLS(ig)%Diag(ixy, iz) * PhiFm(ixy, iz, ig)
      DO ibd = 1, nbd
        ineigh = CMFDLS(ig)%NeighIdx(ibd, ixy)
        IF(ineigh .LE. 0) CYCLE
        lmnt = lmnt + CMFDLS(ig)%RadOffDiag(ibd, ixy, iz0) * PhiFm(ineigh, iz, ig)        
      ENDDO
      !Axial Contribution
      lmnt = lmnt + CMFDLS(ig)%AxOffDiag(2, ixy, iz) * PhiFm(ixy, iz + 1, ig)
      lmnt = lmnt + CMFDLS(ig)%AxOffDiag(1, ixy, iz) * PhiFm(ixy, iz - 1, ig)
      lmnt = SRC(ixy, iz) - lmnt
#else
      vol = PinVolFm(ixy, iz0)
      area = vol / hzfm(iz)
      jsum  = 0
      myphi = PhiFm(ixy, iz, ig)
      !Radial Current
      DO ibd = 1, nbd
        ineigh = PinNeighIdx(ibd, ixy)
        IF(ineigh .EQ. RefCell) CYCLE
        dtil = CMfdPinXs(ixy, iz0)%Dtil(ibd, ig)
        dhat = CMfdPinXs(ixy, iz0)%Dhat(ibd, ig)
        neighphi = 0
        IF(ineigh .GT. 0) neighphi = PhiFm(ineigh, iz, ig)
        jnet = (dtil - dhat)*myphi  -(dtil + dhat)*neighphi
        jsum = jsum + jnet
      ENDDO
      jsum = jsum * hzfm(iz)
      
      Dtil = AxFlx(iz, ixy)%Dtil(1, ig); Dhat = AxFlx(iz, ixy)%Dhat(1, ig)
      neighphi = PhiFM(ixy, iz - 1, ig); 
      jnet = (dtil - dhat)*myphi  -(dtil + dhat)*neighphi
      jsum = jsum + jnet * area        

      Dtil = AxFlx(iz, ixy)%Dtil(2, ig); Dhat = AxFlx(iz, ixy)%Dhat(2, ig)
      neighphi = PhiFM(ixy, iz + 1, ig); 
      jnet = (dtil - dhat)*myphi  -(dtil + dhat)*neighphi
      jsum = jsum + jnet * area     
      
      
      lmnt = jsum + vol*CmfdPinXs(ixy, iz0)%xsr(ig)*myphi
      lmnt = src(ixy, iz) - lmnt
#endif      
      tsrc = tsrc +  src(ixy, iz) * src(ixy, iz)
      ResidualError = ResidualError + lmnt * lmnt
    ENDDO
  ENDDO
ENDDO
#endif
#ifdef MPI_ENV 
buf0  = (/ResidualError, tsrc/)
CALL REDUCE(buf0, buf, 2, PE%MPI_CMFD_COMM, .TRUE.)
ResidualError = buf(1); tsrc = buf(2)
#endif
ResidualError = SQRT(ResidualError/tsrc)

END FUNCTION

SUBROUTINE AddConstCmfdSrc(Src, ConstSrc)
USE PARAM
USE CMFD_MOD, ONLY : ng,          nxy,          myzbf,         myzef,   &
                     SubPlaneMap,  PinVol,      PinVolFm
IMPLICIT NONE
REAL, POINTER :: Src(:,:)
REAL :: ConstSrc
INTEGER :: iz, ixy, iz0

DO iz = myzbf, myzef
  iz0 = SubPlaneMap(iz)
  DO ixy = 1, nxy
    !IF(Src(ixy, iz) .LT. 0._8) Src(ixy, iz) = 0._8
    Src(ixy, iz) =  Src(ixy, iz) + ConstSrc * PinVolFm(ixy, iz0)
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE AddCMFDPxs(Pxs, iz1, iz2, ig1, ig2)
USE PARAM
USE CMFD_MOD, ONLY : nxy, ng, CmfdPinXs 
IMPLICIT NONE
REAL, POINTER :: PXS(:, :, :)
INTEGER :: iz1, iz2, ig1, ig2
INTEGER :: iz, ig, ixy

DO iz = iz1, iz2
  DO ixy = 1, nxy
    DO ig = 1, ng
      IF(Pxs(ixy, iz, ig) .GT. 0.) THEN
        CmfdPinXs(ixy, iz)%xsr(ig) = CmfdPinXs(ixy, iz)%xsr(ig) + Pxs(ixy, iz, ig)
      ELSE
        CmfdPinXs(ixy, iz)%xss(ig)%from(ig) = - Pxs(ixy, iz, ig)
      ENDIF
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE HomBuckling(Core, Fxr, Phis, PinXS, myzb, myze, ng, bsq, lxslib)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,    FXRInfo_Type,      PinXs_Type,        &
                           Pin_Type,         PinInfo_Type,      Cell_Type
USE CMFD_mod,       ONLY : XsMac
USE Core_mod,       ONLY : GroupInfo
USE BenchXs,        ONLY : XsBaseBen
USE MacXsLib_Mod,   ONLY : MacXsBase,        MacXsScatMatrix
USE BasicOperation, ONLY : CP_VA,            CP_CA,             MULTI_VA
IMPLICIT NONE

TYPE(CoreInfo_Type), INTENT(IN) :: Core
TYPE(Fxrinfo_type), POINTER :: Fxr(:,:)
TYPE(PinXs_Type), POINTER, INTENT(INOUT) :: PinXs(:, :)
REAL, POINTER, INTENT(IN) :: Phis(:, :, :)
INTEGER, INTENT(IN) :: myzb, myze, ng
LOGICAL, INTENT(IN) :: lXsLib
REAL :: Bsq

TYPE(FxrInfo_Type), POINTER :: myFXR
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(PinInfo_Type), POINTER :: PinInfo(:)
TYPE(Cell_Type), POINTER :: Cell(:)

INTEGER :: nCoreFxr, nCoreFsr, nxy
INTEGER :: nCellType, nPinType, nlocalFxr, nlocalFsr, nFsrInFxr
INTEGER :: icel, ipin, iz, ixy ,ifxr, ifsr, itype, ifsrlocal
INTEGER :: FsrIdxSt, FxrIdxSt
INTEGER :: i, j, k, l, ig
INTEGER :: iResoGrpBeg, iResoGrpEnd, norg, nChi

REAL :: XsD, PhiSum


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
    icel = Pin(ixy)%Cell(iz)
    nlocalFxr = Cell(icel)%nFxr; nlocalFsr = Cell(icel)%nFsr
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1; nFsrInFxr = Cell(icel)%nFsrInFxr(j)
      XsMac(j)%lFuel = FALSE
      IF(lXsLib) THEN
        myFxr => FXR(ifxr, iz)
        CALL MacXsBase(XsMac(j), myFxr, 1, ng, ng, 1._8, FALSE, TRUE)
        CALL MacXsScatMatrix(XsMac(j), myFxr, 1, ng, ng, GroupInfo, FALSE)   
        !Self-Sheilding Effect
        IF(myFxr%lres) THEN
           do ig = iResoGrpBeg, iResoGrpEnd
             XsMac(j)%XsMacA(ig) = XsMac(j)%XsMacA(ig) * myFxr%FresoA(ig)  
           enddo
        ENDIF
        XsMac(j)%XsMacTr = XsMac(j)%XsMacA + XsMac(j)%XsMacStr
        XsMac(j)%XsMacT = XsMac(j)%XsMacA + XsMac(j)%XsMacS
      ELSE
        ifsrlocal = Cell(icel)%MapFxr2FsrIdx(1,j)
        !itype = Cell(icel)%iReg(ifsrlocal)
        itype=Fxr(ifxr,iz)%imix
        CALL xsbaseBen(itype, 1, ng, 1, ng, FALSE, XsMac(j))
      ENDIF
    ENDDO !Fxr sweep in the Cell
    DO ig = 1, ng
      XsD = 0.; PhiSum = 0
      DO j = 1, nLocalFxr
        ifxr = FxrIdxSt + j -1; nFsrInFxr = Cell(icel)%nFsrInFxr(j)
        DO i = 1, nFsrInFxr
          k = Cell(icel)%MapFxr2FsrIdx(i, j)
                    
          ifsr = FsrIdxSt + k - 1
          PhiSum = PhiSum + Phis(ifsr, iz, ig) * Cell(icel)%vol(k)
          IF((XsMac(j)%xsMacTr(ig)) .LT. epsm3) CYCLE
          XsD = XsD + Phis(ifsr, iz, ig) * Cell(icel)%vol(k) / (3._8 * XsMac(j)%xsMacTr(ig))
        ENDDO
      ENDDO
      XsD = XsD / PhiSum * Bsq
      PinXS(ixy, iz)%XSR(ig) = PinXS(ixy, iz)%XSR(ig) + XsD
      !PinXS(ixy, iz)%XSTR(ig) = PinXS(ixy, iz)%XSTR(ig) + XsD
      !PinXS(ixy, iz)%XST(ig) = PinXS(ixy, iz)%XST(ig) + XsD
    ENDDO
    
  ENDDO
ENDDO

NULLIFY(myFXR)
NULLIFY(Pin)
NULLIFY(PinInfo)
NULLIFY(Cell)

END SUBROUTINE

SUBROUTINE ConvertSubPlnPhi(PhiC, PhiFm, imod)
USE PARAM
USE CMFD_Mod,      ONLY : myzb,          myze,          myzbf,        myzef,     &
                          nxy,           ng,                                     &
                          SubPlaneMap,   SubPlaneRange,nSubPlane  
IMPLICIT NONE
REAL, POINTER :: PhiC(:, :, :), PhiFm(:, :, :)
INTEGER :: imod

INTEGER :: ixy, iz, iz0, iz1, iz2, ig
REAL :: phiavg, nsubpln_inv, fmult

nSubPln_inv = 1._8 / nSubPlane
DO ig = 1, ng
  DO iz0 = myzb, myze
    iz1 = SubPlaneRange(1, iz0); iz2 = SubPlaneRange(2, iz0)
    DO ixy = 1, nxy
      phiavg = sum(PhiFm(ixy, iz1:iz2, ig)) * nSubPln_inv
      IF(imod .eq. 1) THEN
        fmult = PhiC(ixy, iz0, ig) / phiavg
        PhiFm(ixy, iz1:iz2, ig) = PhiFm(ixy, iz1:iz2, ig) * fmult
      ELSE
        PhiC(ixy, iz0, ig) = PhiAvg
      ENDIF
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE

