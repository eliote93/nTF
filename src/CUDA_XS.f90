#include <CUDADEFINES.h>
#ifdef __PGI
MODULE CUDA_XS
USE XS_COMMON
USE CUDA_CONST
USE CUDA_MASTER
USE CUDA_UTIL

IMPLICIT NONE

INTEGER, PARAMETER :: MacBase = 1, MacNF = 2, MacAll = 3
INTEGER, PARAMETER :: ConMac = 1, ConFSP = 2, ConEff = 3
INTEGER, PARAMETER :: BufCon = 1, BufHold = 2, BufErs = 3, BufCnE = 4
INTEGER, PARAMETER :: EffRIP = 1, EffFxrD = 2, EffFxrT = 3, EffClr = 4, EffAll = 5, EffBlkD = 6

INTEGER, PARAMETER :: matGen = 0, matCLD = 1, matAIC = 2

TYPE(cuCoreXsMac_Type) :: cuCoreMacXs ! For all planes in (myzb:myze)
TYPE(cuBlockFxr_Type) :: cuBlockFxr
TYPE(cuIsoData_Type) :: cuIsoData
TYPE(cuMLGLib_Type) :: cuMLGData
TYPE(cuResIsoData_Type) :: cuResIsoData
TYPE(cuResPin_Type) :: cuResPin
TYPE(cuPinInfo_Type) :: cuPinInfo
TYPE(dim3) :: blockXS, blockXSD, blockSmD, blockXST, blockSmT, trdXS, trdSm
TYPE(dim3) :: blockGF, blockCF, blockAF, blockGP, blockCP, blockAP, blockFnAdj, trdFSP
TYPE(dim3) :: blockN, trdN
TYPE(dim3) :: blockRP, blockRG, blockRC, blockRA, blockRD, blockRT, trdEffIso, trdEffTri, trdEffRes

CONTAINS

SUBROUTINE SetConstXsVar
IMPLICIT NONE
INTEGER :: i
DO i = 1, 32
  ifxrbegD(i) = ibfxrD(i)
  ifxrbegT(i) = ibfxrT(i)
END DO
DO i = 1, 64
  lFuelPlane(i) = lFPlane(i)
END DO

END SUBROUTINE

SUBROUTINE InitThreadGridXS()
IMPLICIT NONE
INTEGER :: nFxrD, nBD, nFxrT, nBT, ng, nBG, nBGlobal
nFxrD = FxrD%nFxrD; nFxrT = FxrTri%nFxrTri; ng = IsoData%ng
nBD = nFxrD/XS_BLOCK_RDIM+1; nBT = nFxrT/XS_BLOCK_RDIM+1; nBG = ng/XS_BLOCK_GDIM+1;
nBGlobal = (nFxrD+nFxrT)/XS_BLOCK_RDIM+1
blockXSD = dim3(nBG, nBD, 1); blockXST = dim3(nBG, nBT, 1);
blockSmD = dim3(nBG, ng, nBD); blockSmT = dim3(nBG, ng, nBT);
blockXS = dim3(nBG, nBGlobal, 1)
trdXS = dim3(XS_BLOCK_GDIM, XS_BLOCK_RDIM, 1);
trdSm = dim3(XS_BLOCK_GDIM, 1, XS_BLOCK_RDIM);
END SUBROUTINE

SUBROUTINE InitThreadGridFSP()
USE Core_mod, ONLY : GroupInfo
USE GEOM, ONLY : core
USE PE_Mod, ONLY : PE
IMPLICIT NONE
INTEGER :: nFxrD, nFxrT
INTEGER :: nFxr, nPin, nBF, nBP, nrg, nfmaclv, ncmaclv, nAmaclv, nBrG, nBflv, nBclv, nBAlv
nFxrD = FxrD%nFxrD; nFxrT = FxrTri%nFxrTri; nrg = GroupInfo%norg
nFXR = nFxrD+nFxrT; nPin = core%nxy*(PE%myze-PE%myzb+1)
nfmaclv = MLGLib%f_nmaclv; ncmaclv = MLGLib%c_nmaclv1G; nAmaclv = MLGLib%f_nmaclv1G

nBF = nFxr/FSP_BLOCK_RDIM+1; nBP = nPin/FSP_BLOCK_RDIM+1;
nBrG = nrg/FSP_BLOCK_GDIM+1; nBflv = (nrg*nfmaclv)/FSP_BLOCK_GDIM+1;
nBclv = ncmaclv/FSP_BLOCK_GDIM+1; nBAlv = nAmaclv/FSP_BLOCK_GDIM+1;

blockGF = dim3(nBflv,nBF,1); blockCF = dim3(nBclv,nBF,1); blockAF = dim3(nBAlv,nBF,1)
blockGP = dim3(nBflv,nBP,1); blockCP = dim3(nBclv,nBP,1); blockAP = dim3(nBAlv,nBP,1)
blockFnAdj = dim3(nBrG,nBF,1)
trdFSP = dim3(FSP_BLOCK_GDIM, FSP_BLOCK_RDIM, 1);

trdN = dim3(512,1,1); blockN = dim3(1,nrg*nfmaclv,1);
END SUBROUTINE

SUBROUTINE InitThreadGridEff()
USE Core_mod, ONLY : GroupInfo
USE GEOM, ONLY : core
USE PE_Mod, ONLY : PE
IMPLICIT NONE
INTEGER :: npin, nGenD, nCldT, nAICT, nResD, nResT, ntiso, nrg
INTEGER :: nBP, nBG, nBC, nBA, nBI, nBD, nBT, nBR

npin = core%nxy*(PE%myze-PE%myzb+1); nGenD = FxrD%nGenD; nCldT = FxrTri%nCldTri
nAICT = FxrTri%nAICTri; nResD =FxrD%nResD; nResT = FxrTri%nResTri;
ntiso = IsoData%ntiso; nrg = GroupInfo%norg

nBP = npin/EFF_BLOCK_RDIM+1; nBG = nGenD/EFF_BLOCK_RDIM+1;
nBC = nCldT/EFF_BLOCK_RDIM+1; nBA = nAICT/EFF_BLOCK_RDIM+1;
nBI = ntiso/EFF_BLOCK_IDIM+1; nBD = nResD/EFF_BLOCK_RDIM+1;
nBT = nResT/EFF_BLOCK_RDIM+1; nBR = nrg/EFF_BLOCK_GDIM+1;

blockRP = dim3(nBP,nBI,nrg); blockRG = dim3(nBG,nBI,nrg); blockRC = dim3(nBC,nBR,1)
blockRA = dim3(nBA,nBR,1); blockRD = dim3(1,nBD,1); blockRT =  dim3(1,nBT,1)

trdEffIso = dim3(EFF_BLOCK_RDIM,EFF_BLOCK_IDIM,1);
trdEffTri = dim3(EFF_BLOCK_RDIM,EFF_BLOCK_GDIM,1);
trdEffRes = dim3(nrg,EFF_BLOCK_RDIM,1)
END SUBROUTINE

SUBROUTINE AllocCuCoreMacXs(ng, ScatOrder, caseID)
USE CORE_MOD,     ONLY : GroupInfo
USE GEOM,       ONLY : Core
IMPLICIT NONE
INTEGER :: ng, ScatOrder, caseID

INTEGER :: i, j, k
INTEGER :: ngtot, ngs, nFxr, nFsr

INTEGER :: istat

ngtot = IsoData%ng; nFxr = core%ncorefxr
ngs = 0;
DO i = 1, ngtot
  DO j = GroupInfo%InScatRange(1,i), GroupInfo%InScatRange(2,i)
    ngs=ngs+1
  END DO
END DO
IF (ng.NE.ngtot) ngs = ng*ngtot

IF (caseID .NE. MacBase) THEN
  ALLOCATE(cuCoreMacXs%XSnf(ng,nfxr))
  istat = cudaMemSet(cuCoreMacXs%XSnf,0.,ng*nFxr)
END IF
IF (caseID .NE. MacNF) THEN
  ALLOCATE(cuCoreMacXs%XSt(ng,nfxr))
  ALLOCATE(cuCoreMacXs%XStr(ng,nfxr))
  ALLOCATE(cuCoreMacXs%XSS(ng,nfxr))
  ALLOCATE(cuCoreMacXs%XSStr(ng,nfxr))
  ALLOCATE(cuCoreMacXs%XSa(ng,nfxr))
  ALLOCATE(cuCoreMacXs%XSkf(ng,nfxr))
  istat = cudaMemSet(cuCoreMacXs%XSkf,0.,ng*nFxr)

  ALLOCATE(cuCoreMacXs%XSsm(ngs,nfxr))
  IF (ScatOrder<1) return
  ALLOCATE(cuCoreMacXs%XSsmP1(ngs,nfxr))
  ALLOCATE(cuCoreMacXs%XSsmP2(ngs,nfxr))
  ALLOCATE(cuCoreMacXs%XSsmP3(ngs,nfxr))
END IF

END SUBROUTINE

SUBROUTINE DestroyCuCoreMacXs(caseid, ScatOrder)
IMPLICIT NONE
INTEGER :: caseid, ScatOrder
IF (caseid.NE.MacBase) THEN
  DEALLOCATE(cuCoreMacXs%XSnf)
END IF
IF (caseid.NE.MacNF) THEN
  DEALLOCATE(cuCoreMacXs%XStr)
  DEALLOCATE(cuCoreMacXs%XSsm)
  DEALLOCATE(cuCoreMacXs%XSt)
  DEALLOCATE(cuCoreMacXs%XSS)
  DEALLOCATE(cuCoreMacXs%XSStr)
  DEALLOCATE(cuCoreMacXs%XSa)
  DEALLOCATE(cuCoreMacXs%XSkf)
  IF (ScatOrder<1) return
  DEALLOCATE(cuCoreMacXs%XSsmP1)
  DEALLOCATE(cuCoreMacXs%XSsmP2)
  DEALLOCATE(cuCoreMacXs%XSsmP3)
END IF
END SUBROUTINE

SUBROUTINE ConstructCuIsoData
IMPLICIT NONE
INTEGER :: ntiso, ng, nttemp, nttempPx, Np0, Np1, Np2, Np3

ntiso = IsoData%ntiso; ng = IsoData%ng;
nttemp = IsoData%nttemp; nttempPx = IsoData%nttempPx;
Np0 = IsoData%Np0; Np1 = IsoData%Np1; Np2 = IsoData%Np2; Np3 = IsoData%Np3

ALLOCATE(cuIsoData%ptsTemp(nttemp), cuIsoData%ptsTempPx(nttempPx))
ALLOCATE(cuIsoData%ptrTemp(ntiso+1), cuIsoData%ptrTempPx(ntiso+1))
ALLOCATE(cuIsoData%siga(ng,nttemp))
ALLOCATE(cuIsoData%sigf(ng,nttemp))
ALLOCATE(cuIsoData%sigkf(ng,nttemp))
ALLOCATE(cuIsoData%signf(ng,nttemp))
ALLOCATE(cuIsoData%sigs(ng,nttemp))
ALLOCATE(cuIsoData%sigstr(ng,nttemp))
ALLOCATE(cuIsoData%gidSmp0(2*ng,nttemp), cuIsoData%ptrSmp0(ng,nttemp))
ALLOCATE(cuIsoData%gidSmp1(2*ng,nttempPx), cuIsoData%ptrSmp1(ng,nttempPx))
ALLOCATE(cuIsoData%gidSmp2(2*ng,nttempPx), cuIsoData%ptrSmp2(ng,nttempPx))
ALLOCATE(cuIsoData%gidSmp3(2*ng,nttempPx), cuIsoData%ptrSmp3(ng,nttempPx))
ALLOCATE(cuIsoData%smp0(Np0),cuIsoData%smp1(Np1),cuIsoData%smp2(Np2),cuIsoData%smp3(Np3))
cuIsoData%ntiso = ntiso; cuIsoData%ng = ng;
cuIsoData%ptsTemp = IsoData%ptsTemp;
cuIsoData%ptsTempPx = IsoData%ptsTempPx
cuIsoData%ptrTemp = IsoData%ptrTemp; cuIsoData%ptrTempPx = IsoData%ptrTempPx
cuIsoData%siga = IsoData%siga; cuIsoData%sigf = IsoData%sigf; cuIsoData%sigkf = IsoData%sigkf
cuIsoData%signf = IsoData%signf; cuIsoData%sigs = IsoData%sigs; cuIsoData%sigstr = IsoData%sigstr
cuIsoData%gidSmp0 = IsoData%gidSmp0; cuIsoData%gidSmp1 = IsoData%gidSmp1
cuIsoData%gidSmp2 = IsoData%gidSmp2; cuIsoData%gidSmp3 = IsoData%gidSmp3
cuIsoData%ptrSmp0 = IsoData%ptrSmp0; cuIsoData%ptrSmp1 = IsoData%ptrSmp1
cuIsoData%ptrSmp2 = IsoData%ptrSmp2; cuIsoData%ptrSmp3 = IsoData%ptrSmp3
cuIsoData%smp0 = IsoData%smp0; cuIsoData%smp1 = IsoData%smp1
cuIsoData%smp2 = IsoData%smp2; cuIsoData%smp3 = IsoData%smp3
END SUBROUTINE

SUBROUTINE DestroyCuIsoData
IMPLICIT NONE
! No use. Keep it until the program stops
END SUBROUTINE

SUBROUTINE ConstructcuResIsoData
USE Core_mod, ONLY : GroupInfo
IMPLICIT NONE
INTEGER :: ntiso, nresiso, ntsig0, ntlv, ntrtemp, NTempPot, NtempLv, ntrif, ntrifRat, ntrifRx, irgb, irge
IF (.NOT. lrestrmt) RETURN
ntiso = IsoData%ntiso
irgb = GroupInfo%nofg+1; irge = GroupInfo%nofg+GroupInfo%norg-1
nresiso = ResIsoData%nresiso; ntsig0 = ResIsoData%ntsig0
ntlv = ResIsoData%ntlv; ntrtemp = ResIsoData%ntrtemp
ntempPot = ResIsoData%ntempPot; ntempLv = ResIsoData%ntempLv
cuResIsoData%nresiso = ResIsoData%nresiso; cuResIsoData%ntsig0 = ResIsoData%ntsig0
cuResIsoData%ntlv = ResIsoData%ntlv; cuResIsoData%ntrtemp = ResIsoData%ntrtemp
cuResIsoData%ntempPot = ResIsoData%ntempPot; cuResIsoData%ntempLv = ResIsoData%ntempLv
ntrif = ResIsoData%ntrif; ntrifRat = ResIsoData%ntrifRat; ntrifRx = ResIsoData%ntRifRx

ALLOCATE(cuResIsoData%IsoG2R(ntiso))
ALLOCATE(cuResIsoData%lclad(ntiso), cuResIsoData%ityp(nresiso), cuResIsoData%sigp(nresiso))
ALLOCATE(cuResIsoData%ptrRTemp(nresiso+1),cuResIsoData%ptrSig0(nresiso+1),cuResIsoData%ptrNlv(nresiso+1))
ALLOCATE(cuResIsoData%ptrTempPot(nresiso+1),cuResIsoData%ptrTempLv(nresiso+1))
ALLOCATE(cuResIsoData%lamsigp(irgb:irge,nresiso))
ALLOCATE(cuResIsoData%sig0sq(ntsig0), cuResIsoData%rtempsq(ntrtemp))
ALLOCATE(cuResIsoData%xsalog(NTempPot,irgb:irge),cuResIsoData%ri_a(NTempPot,irgb:irge))
ALLOCATE(cuResIsoData%lvabs(ntlv,irgb:irge))
ALLOCATE(cuResIsoData%lvfis(ntlv,irgb:irge))
ALLOCATE(cuResIsoData%wgtabs(NTempLv,irgb:irge),cuResIsoData%wgtfis(NTempLv,irgb:irge))

ALLOCATE(cuResIsoData%ptrRifRat(ntrif+1),cuResIsoData%ptrRifRx(ntrifRat+1))
ALLOCATE(cuResIsoData%ratlog(ntrifRat),cuResIsoData%rifabs(ntrifRx,irgb:irge),cuResIsoData%riffis(ntrifRx,irgb:irge))
ALLOCATE(cuResIsoData%rifid(ntiso,nresiso))

cuResIsoData%IsoG2R = ResIsoData%IsoG2R
cuResIsoData%lclad = ResIsoData%lclad; cuResIsoData%ityp = ResIsoData%ityp; cuResIsoData%sigp = ResIsoData%sigp
cuResIsoData%ptrRTemp = ResIsoData%ptrRTemp; cuResIsoData%ptrSig0 = ResIsoData%ptrSig0
cuResIsoData%ptrNlv = ResIsoData%ptrNlv; cuResIsoData%ptrTempPot = ResIsoData%ptrTempPot
cuResIsoData%ptrTempLv = ResIsoData%ptrTempLv
cuResIsoData%lamsigp = ResIsoData%lamsigp; cuResIsoData%sig0sq = ResIsoData%sig0sq
cuResIsoData%rtempsq = ResIsoData%rtempsq; cuResIsoData%xsalog = ResIsoData%xsalog
cuResIsoData%ri_a = ResIsoData%ri_a; cuResIsoData%lvabs = ResIsoData%lvabs
cuResIsoData%lvfis = ResIsoData%lvfis
cuResIsoData%wgtabs = ResIsoData%wgtabs; cuResIsoData%wgtfis = ResIsoData%wgtfis

cuResIsoData%ptrRifRat = ResIsoData%ptrRifRat; cuResIsoData%ptrRifRx = ResIsoData%ptrRifRx
cuResIsoData%ratlog = ResIsoData%ratlog; cuResIsoData%rifid = ResIsoData%rifid
cuResIsoData%rifabs = ResIsoData%rifabs; cuResIsoData%riffis = ResIsoData%riffis
END SUBROUTINE

SUBROUTINE DestroyResIsoData
! No use
END SUBROUTINE

SUBROUTINE ConstructMLGData
USE PE_Mod, ONLY : PE
IMPLICIT NONE
INTEGER :: myzb, myze, nMLvFMg, nMLvF1g, nMLvC1g
IF (.NOT. lMLG .OR. .NOT. lrestrmt) RETURN
myzb = PE%myzb; myze = PE%myze
nMLvFMg = MLGLib%f_nmaclv; nMLvF1g = MLGLib%f_nmaclv1G; nMLvC1g = MLGLib%c_nmaclv1G
cuMLGData%f_nmaclv = nMLvFMg; cuMLGData%f_nmaclv1G = nMLvF1g; cuMLGData%c_nmaclv1G = nMLvC1g
ALLOCATE(cuMLGData%c_maclv1G(nMLvC1g),cuMLGData%c_maclv1G_log(nMLvC1g))
ALLOCATE(cuMLGData%f_maclv(nMLvFMg,myzb:myze),cuMLGData%f_maclv_log(nMLvFMg,myzb:myze))
ALLOCATE(cuMLGData%f_maclv1G(nMLvF1g,myzb:myze),cuMLGData%f_maclv1G_log(nMLvF1g,myzb:myze))
!print*, 'What'
cuMLGData%c_maclv1G = MLGLib%c_maclv1G; cuMLGData%c_maclv1G_log = MLGLib%c_maclv1G_log
cuMLGData%f_maclv = MLGLib%f_maclv; cuMLGData%f_maclv_log = MLGLib%f_maclv_log;
cuMLGData%f_maclv1G = MLGLib%f_maclv1G; cuMLGData%f_maclv1G_log = MLGLib%f_maclv1G_log
!print*, MLGLib%f_maclv
END SUBROUTINE

SUBROUTINE DestroyMLGData
! No use
END SUBROUTINE

SUBROUTINE ConstructCuBlockFxr(caseid, caseEff)
USE PE_Mod, ONLY : PE
USE Core_mod, ONLY : GroupInfo
USE XSLIB_MOD,  ONLY : mlgdata0
IMPLICIT NONE
INTEGER :: caseid, caseEff
INTEGER :: nFxrD, nFxrTri, NtTri, nResTri, nResD, ntfsrD, ntfsrT, ng, ntiso, ibres, ieres
INTEGER :: nGenD, nCldD, nGenTri, nCldTri, nAICTri, NtRTri, nMLvC1G, nMLvF1G, nMLvFMg

INTEGER :: istat

IF (caseid.EQ.ConMac) THEN
  ntfsrD = FxrD%ntfsr; ntfsrT = FxrTri%ntfsr
  nFxrD = FxrD%nFxrD; nFxrTri = FxrTri%nFxrTri; NtTri = FxrTri%NtTri
  nResD = FxrD%nResD; nResTri = FxrTri%nResTri
  ng = IsoData%ng; ntiso = IsoData%ntiso
  ibres = GroupInfo%nofg+1; ieres = ibres-1+GroupInfo%norg
  IF (nFxrD.GT.0) THEN
    ALLOCATE(cuBlockFxr%nisoD(nFxrD))
    ALLOCATE(cuBlockFxr%mapglobalidD(nFxrD))
    ALLOCATE(cuBlockFxr%mapG2RD(nFxrD))
    ALLOCATE(cuBlockFxr%pinidD(nFxrD))
    ALLOCATE(cuBlockFxr%fsrstD(nFxrD))
    ALLOCATE(cuBlockFxr%mapfsrD(ntfsrD))
    ALLOCATE(cuBlockFxr%nfsrD(nFxrD))
    ALLOCATE(cuBlockFxr%pnumD(ntiso,nFxrD))
    ALLOCATE(cuBlockFxr%lvoidD(nFxrD))
    ALLOCATE(cuBlockFxr%MapNuclD(ntiso,nFxrD))
    ALLOCATE(cuBlockFxr%chi(GroupInfo%nchi,nFxrD))
    ALLOCATE(cuBlockFxr%wtTempD(nFxrD,ntiso),cuBlockFxr%itTempD(nFxrD,ntiso))

!    cuBlockFxr%nisoD = FxrD%niso
!    cuBlockFxr%mapglobalidD = Fxrd%mapglobalid;
!    cuBlockFxr%mapG2RD = FxrD%mapG2R;
!    cuBlockFxr%pinidD = Fxrd%ipin;
!    cuBlockFxr%fsrstD = Fxrd%fsrst; cuBlockFxr%nfsrD = Fxrd%nfsr;
!    cuBlockFxr%mapfsrD = FxrD%mapfsr; cuBlockFxr%pnumD = Fxrd%pnum;
!    cuBlockFxr%lvoidD = Fxrd%lvoid;
!    cuBlockFxr%MapNuclD = FxrD%MapNucl; cuBlockFxr%chi = FxrD%chi

    istat = cudaMemCpy(cuBlockFxr%nisoD,FxrD%niso,nFxrD)
    istat = cudaMemCpy(cuBlockFxr%mapglobalidD,FxrD%mapglobalid,nFxrD)
    istat = cudaMemCpy(cuBlockFxr%mapG2RD,FxrD%mapG2R,nFxrD)
    istat = cudaMemCpy(cuBlockFxr%pinidD,FxrD%ipin,nFxrD)
    istat = cudaMemCpy(cuBlockFxr%fsrstD,FxrD%fsrst,nFxrD)
    istat = cudaMemCpy(cuBlockFxr%nfsrD,FxrD%nfsr,nFxrD)
    istat = cudaMemCpy(cuBlockFxr%mapfsrD,FxrD%mapfsr,ntfsrD)
    istat = cudaMemCpy(cuBlockFxr%pnumD,FxrD%pnum,ntiso*nFxrD)
    istat = cudaMemCpy(cuBlockFxr%lvoidD,FxrD%lvoid,nFxrD)
    istat = cudaMemCpy(cuBlockFxr%MapNuclD,FxrD%MapNucl,ntiso*nFxrD)
    istat = cudaMemCpy(cuBlockFxr%chi,FxrD%chi,GroupInfo%nchi*nFxrD)
    istat = cudaMemCpy(cuBlockFxr%wtTempD,FxrD%wtTemp,ntiso*nFxrD)
    istat = cudaMemCpy(cuBlockFxr%itTempD,FxrD%itTemp,ntiso*nFxrD)
  END IF
  IF (nResD.GT.0) THEN
    ALLOCATE(cuBlockFxr%FresoAD(ibres:ieres, nResD))
    ALLOCATE(cuBlockFxr%FresoF(ibres:ieres,nResD), cuBlockFxr%FresoNF(ibres:ieres,nResD))
!    cuBlockFxr%FresoAD = Fxrd%fresoa; cuBlockFxr%FresoF = Fxrd%fresof; cuBlockFxr%FresoNF = Fxrd%fresonf;
    istat = cudaMemCpy(cuBlockFxr%FresoAD,Fxrd%fresoa,(ieres-ibres+1)*nResD)
    istat = cudaMemCpy(cuBlockFxr%FresoF,Fxrd%fresokf,(ieres-ibres+1)*nResD)
    istat = cudaMemCpy(cuBlockFxr%FresoNF,Fxrd%fresonf,(ieres-ibres+1)*nResD)
  END IF
  IF (nResTri.GT.0) THEN
    ALLOCATE(cuBlockFxr%FresoAT(ibres:ieres, nResTri))
!    cuBlockFxr%FresoAT = FxrTri%fresoa;
    istat = cudaMemCpy(cuBlockFxr%FresoAT,FxrTri%fresoa,(ieres-ibres+1)*nResTri)
  END IF
  IF (nFxrTri.GT.0) THEN
    ALLOCATE(cuBlockFxr%nisoT(nFxrTri), cuBlockFxr%AccNT(nFxrTri+1))
    ALLOCATE(cuBlockFxr%mapglobalidT(nFxrTri))
    ALLOCATE(cuBlockFxr%mapG2RT(nFxrTri))
    ALLOCATE(cuBlockFxr%pinidT(nFxrTri))
    ALLOCATE(cuBlockFxr%fsrstT(nFxrTri))
    ALLOCATE( cuBlockFxr%mapfsrT(ntfsrT))
    ALLOCATE(cuBlockFxr%nfsrT(nFxrTri))
    ALLOCATE(cuBlockFxr%pnumT(NtTri))
    ALLOCATE(cuBlockFxr%lvoidT(nFxrTri))
    ALLOCATE(cuBlockFxr%MapnuclT(NtTri))
    ALLOCATE(cuBlockFxr%itTempT(nFxrTri,FxrTri%nisomax))
    ALLOCATE(cuBlockFxr%wtTempT(nFxrTri,FxrTri%nisomax))

!    cuBlockFxr%nisoT = FxrTri%niso; cuBlockFxr%AccNT = FxrTri%AccNt
!    cuBlockFxr%mapglobalidT = FxrTri%mapglobalid;
!    cuBlockFxr%mapG2RT = FxrTri%mapG2R
!    cuBlockFxr%pinidT = FxrTri%ipin;
!    cuBlockFxr%fsrstT = FxrTri%fsrst; cuBlockFxr%nfsrT = FxrTri%nfsr;
!    cuBlockFxr%mapfsrT = FxrTri%mapfsr; cuBlockFxr%pnumT = FxrTri%pnum;cuBlockFxr%lvoidT = FxrTri%lvoid;
!    cuBlockFxr%MapnuclT = FxrTri%MapNucl

    istat = cudaMemCpy(cuBlockFxr%nisoT,FxrTri%niso,nFxrTri)
    istat = cudaMemCpy(cuBlockFxr%AccNT,FxrTri%AccNt,nFxrTri+1)
    istat = cudaMemCpy(cuBlockFxr%mapglobalidT,FxrTri%mapglobalid,nFxrTri)
    istat = cudaMemCpy(cuBlockFxr%mapG2RT,FxrTri%mapG2R,nFxrTri)
    istat = cudaMemCpy(cuBlockFxr%pinidT,FxrTri%ipin,nFxrTri)
    istat = cudaMemCpy(cuBlockFxr%fsrstT,FxrTri%fsrst,nFxrTri)
    istat = cudaMemCpy(cuBlockFxr%nfsrT,FxrTri%nfsr,nFxrTri)
    istat = cudaMemCpy(cuBlockFxr%mapfsrT,FxrTri%mapfsr,ntfsrT)
    istat = cudaMemCpy(cuBlockFxr%pnumT,FxrTri%pnum,NtTri)
    istat = cudaMemCpy(cuBlockFxr%lvoidT,FxrTri%lvoid,nFxrTri)
    istat = cudaMemCpy(cuBlockFxr%MapnuclT,FxrTri%MapNucl,NtTri)
    istat = cudaMemCpy(cuBlockFxr%itTempT,FxrTri%itTemp,nFxrTri*FxrTri%nisomax)
    istat = cudaMemCpy(cuBlockFxr%wtTempT,FxrTri%wtTemp,nFxrTri*FxrTri%nisomax)
  END IF
ELSE IF (caseid .EQ. ConFSP) THEN
  ntfsrD = FxrD%ntfsr; ntfsrT = FxrTri%ntfsr
  nFxrD = FxrD%nFxrD; nFxrTri = FxrTri%nFxrTri; NtTri = FxrTri%NtTri
  ng = IsoData%ng; ntiso = IsoData%ntiso
  ibres = GroupInfo%nofg+1; ieres = ibres-1+GroupInfo%norg
  nGenD = FxrD%nGenD; nCldD = FxrD%nCldD;
  nGenTri = FxrTri%nGenTri; nCldTri = FxrTri%nCldTri; nAICTri = FxrTri%nAICTri;
  NtRTri = FxrTri%NtRTri
  nMLvC1G = mlgdata0%c_nmaclv1G; nMLvF1G = mlgdata0%f_nmaclv1G; nMLvFMg = mlgdata0%f_nmaclv
  !print*, 'Fxr Sorting Info', nGenD, nGenTri, nCldD, nCldTri

  IF (nFxrD.GT.0) THEN
    ALLOCATE(cuBlockFxr%nisoD(nFxrD))
    ALLOCATE(cuBlockFxr%areaD(nFxrD))
    ALLOCATE(cuBlockFxr%tempD(nFxrD))
    ALLOCATE(cuBlockFxr%mapglobalidD(nFxrD))
    ALLOCATE(cuBlockFxr%pinidD(nFxrD))
    ALLOCATE(cuBlockFxr%fsrstD(nFxrD))
    ALLOCATE(cuBlockFxr%mapfsrD(ntfsrD))
    ALLOCATE(cuBlockFxr%nfsrD(nFxrD))
    ALLOCATE(cuBlockFxr%pnumD(ntiso,nFxrD), cuBlockFxr%indD(ntiso,nFxrD))
    ALLOCATE(cuBlockFxr%MapNuclD(ntiso,nFxrD))
    ALLOCATE(cuBlockFxr%fsrvolD(ntfsrD))

!    cuBlockFxr%nisoD = FxrD%niso; cuBlockFxr%areaD = FxrD%area; cuBlockFxr%tempD = FxrD%temp
!    cuBlockFxr%mapglobalidD = Fxrd%mapglobalid;
!    cuBlockFxr%pinidD = Fxrd%ipin;
!    cuBlockFxr%fsrstD = Fxrd%fsrst; cuBlockFxr%nfsrD = Fxrd%nfsr;
!    cuBlockFxr%mapfsrD = FxrD%mapfsr; cuBlockFxr%pnumD = Fxrd%pnum; cuBlockFxr%indD = FxrD%ind
!    cuBlockFxr%MapNuclD = FxrD%MapNucl;
!    cuBlockFxr%fsrvolD = FxrD%fsrvol

    istat = cudaMemCpy(cuBlockFxr%nisoD,FxrD%niso,nFxrD)
    istat = cudaMemCpy(cuBlockFxr%areaD,FxrD%area,nFxrD)
    istat = cudaMemCpy(cuBlockFxr%tempD,FxrD%temp,nFxrD)
    istat = cudaMemCpy(cuBlockFxr%mapglobalidD,FxrD%mapglobalid,nFxrD)
    istat = cudaMemCpy(cuBlockFxr%pinidD,FxrD%ipin,nFxrD)
    istat = cudaMemCpy(cuBlockFxr%fsrstD,FxrD%fsrst,nFxrD)
    istat = cudaMemCpy(cuBlockFxr%nfsrD,FxrD%nfsr,nFxrD)
    istat = cudaMemCpy(cuBlockFxr%mapfsrD,FxrD%mapfsr,ntfsrD)
    istat = cudaMemCpy(cuBlockFxr%pnumD,FxrD%pnum,ntiso*nFxrD)
    istat = cudaMemCpy(cuBlockFxr%indD,FxrD%ind,ntiso*nFxrD)
    istat = cudaMemCpy(cuBlockFxr%MapNuclD,FxrD%MapNucl,ntiso*nFxrD)
    istat = cudaMemCpy(cuBlockFxr%fsrvolD,FxrD%fsrvol,ntfsrD)
  END IF
  IF (nFxrTri.GT.0) THEN
    ALLOCATE(cuBlockFxr%nisoT(nFxrTri))
    ALLOCATE(cuBlockFxr%areaT(nFxrTri))
    ALLOCATE(cuBlockFxr%tempT(nFxrTri))
    ALLOCATE(cuBlockFxr%AccNT(nFxrTri+1))
    ALLOCATE(cuBlockFxr%mapglobalidT(nFxrTri))
    ALLOCATE(cuBlockFxr%pinidT(nFxrTri))
    ALLOCATE(cuBlockFxr%fsrstT(nFxrTri))
    ALLOCATE(cuBlockFxr%mapfsrT(ntfsrT))
    ALLOCATE(cuBlockFxr%nfsrT(nFxrTri))
    ALLOCATE(cuBlockFxr%pnumT(NtTri))
    ALLOCATE(cuBlockFxr%MapnuclT(NtTri))
    ALLOCATE(cuBlockFxr%fsrvolT(ntfsrT))

!    cuBlockFxr%nisoT = FxrTri%niso; cuBlockFxr%areaT = FxrTri%area; cuBlockFxr%tempT = FxrTri%temp
!    cuBlockFxr%AccNT = FxrTri%AccNt
!    cuBlockFxr%mapglobalidT = FxrTri%mapglobalid;
!    cuBlockFxr%pinidT = FxrTri%ipin;
!    cuBlockFxr%fsrstT = FxrTri%fsrst; cuBlockFxr%nfsrT = FxrTri%nfsr;
!    cuBlockFxr%mapfsrT = FxrTri%mapfsr; cuBlockFxr%pnumT = FxrTri%pnum;
!    cuBlockFxr%MapnuclT = FxrTri%MapNucl;
!    cuBlockFxr%fsrvolT = FxrTri%fsrvol

    istat = cudaMemCpy(cuBlockFxr%nisoT,FxrTri%niso,nFxrTri)
    istat = cudaMemCpy(cuBlockFxr%areaT,FxrTri%area,nFxrTri)
    istat = cudaMemCpy(cuBlockFxr%tempT,FxrTri%temp,nFxrTri)
    istat = cudaMemCpy(cuBlockFxr%AccNT,FxrTri%AccNt,nFxrTri+1)
    istat = cudaMemCpy(cuBlockFxr%mapglobalidT,FxrTri%mapglobalid,nFxrTri)
    istat = cudaMemCpy(cuBlockFxr%pinidT,FxrTri%ipin,nFxrTri)
    istat = cudaMemCpy(cuBlockFxr%fsrstT,FxrTri%fsrst,nFxrTri)
    istat = cudaMemCpy(cuBlockFxr%nfsrT,FxrTri%nfsr,nFxrTri)
    istat = cudaMemCpy(cuBlockFxr%mapfsrT,FxrTri%mapfsr,ntfsrT)
    istat = cudaMemCpy(cuBlockFxr%pnumT,FxrTri%pnum,NtTri)
    istat = cudaMemCpy(cuBlockFxr%MapnuclT,FxrTri%MapNucl,NtTri)
    istat = cudaMemCpy(cuBlockFxr%fsrvolT,FxrTri%fsrvol,ntfsrT)
  END IF

  IF (nGenD.GT.0) THEN
    ALLOCATE(cuBlockFxr%mapGR2GD(nGenD))
    ALLOCATE(cuBlockFxr%xseq_f_mgD(nMLvFMg,ibres:ieres,nGenD))
    ALLOCATE(cuBlockFxr%FnAdjD(ibres:ieres,nGenD))
    ALLOCATE(cuBlockFxr%FtAdjD(nMLvFMg,ibres:ieres,nGenD))

!    cuBlockFxr%mapGR2GD = FxrD%mapGR2G;
!    cuBlockFxr%xseq_f_mgD = FxrD%xseq_f_mg
!    cuBlockFxr%FnAdjD = FxrD%FnAdj; cuBlockFxr%FtAdjD = FxrD%Ftadj

    istat = cudaMemCpy(cuBlockFxr%mapGR2GD,FxrD%mapGR2G,nGenD)
    istat = cudaMemCpy(cuBlockFxr%xseq_f_mgD,FxrD%xseq_f_mg,nMLvFMg*(ieres-ibres+1)*nGenD)
    istat = cudaMemCpy(cuBlockFxr%FnAdjD,FxrD%FnAdj,(ieres-ibres+1)*nGenD)
    istat = cudaMemCpy(cuBlockFxr%FtAdjD,FxrD%Ftadj,nMLvFMg*(ieres-ibres+1)*nGenD)
  END IF
  IF (nGenTri.GT.0) THEN
    ALLOCATE(cuBlockFxr%mapGR2GT(nGenTri))
    ALLOCATE(cuBlockFxr%xseq_f_mgT(nMLvFMg,ibres:ieres,nGenTri))
    ALLOCATE(cuBlockFxr%FnAdjT(ibres:ieres,nGenTri))
    ALLOCATE(cuBlockFxr%FtAdjT(nMLvFMg,ibres:ieres,nGenTri))

!    cuBlockFxr%mapGR2GT = FxrTri%mapGR2G
!    cuBlockFxr%xseq_f_mgT = FxrTri%xseq_f_mg
!    cuBlockFxr%FnAdjT = FxrTri%FnAdj; cuBlockFxr%FtAdjT = FxrTri%Ftadj

    istat = cudaMemCpy(cuBlockFxr%mapGR2GT,FxrTri%mapGR2G,nGenTri)
    istat = cudaMemCpy(cuBlockFxr%xseq_f_mgT,FxrTri%xseq_f_mg,nMLvFMg*(ieres-ibres+1)*nGenTri)
    istat = cudaMemCpy(cuBlockFxr%FnAdjT,FxrTri%FnAdj,(ieres-ibres+1)*nGenTri)
    istat = cudaMemCpy(cuBlockFxr%FtAdjT,FxrTri%Ftadj,nMLvFMg*(ieres-ibres+1)*nGenTri)
  END IF
  IF (nCldD.GT.0) THEN
    ALLOCATE(cuBlockFxr%mapC2GD(nCldD))
    ALLOCATE(cuBlockFxr%xseq_c_1gD(nMLvC1G,nCldD))

!    cuBlockFxr%mapC2GD = FxrD%mapC2G;
!    cuBlockFxr%xseq_c_1gD = FxrD%xseq_c_1g

    istat = cudaMemCpy(cuBlockFxr%mapC2GD,FxrD%mapC2G,nCldD)
    istat = cudaMemCpy(cuBlockFxr%xseq_c_1gD,FxrD%xseq_c_1g,nMLvC1G*nCldD)
  END IF
  IF (nCldTri.GT.0) THEN
    ALLOCATE(cuBlockFxr%mapC2GT(nCldTri))
    ALLOCATE(cuBlockFxr%xseq_c_1gT(nMLvC1G,nCldTri))

!    cuBlockFxr%mapC2GT = FxrTri%mapC2G;
!    cuBlockFxr%xseq_c_1gT = FxrTri%xseq_c_1g

    istat = cudaMemCpy(cuBlockFxr%mapC2GT,FxrTri%mapC2G,nCldTri)
    istat = cudaMemCpy(cuBlockFxr%xseq_c_1gT,FxrTri%xseq_c_1g,nMLvC1G*nCldTri)
  END IF
  IF(nAICTri.GT.0) THEN
    ALLOCATE(cuBlockFxr%mapA2GT(nAICTri))
    ALLOCATE(cuBlockFxr%xseq_f_1gT(nMLvF1G,nAICTri))

!    cuBlockFxr%mapA2GT = FxrTri%mapA2G
!    cuBlockFxr%xseq_f_1gT = FxrTri%xseq_f_1g

    istat = cudaMemCpy(cuBlockFxr%mapA2GT,FxrTri%mapA2G,nAICTri)
    istat = cudaMemCpy(cuBlockFxr%xseq_f_1gT,FxrTri%xseq_f_1g,nMLvF1G*nAICTri)
  END IF
ELSE IF (caseid.EQ.ConEff) THEN

  nFxrD = FxrD%nFxrD; nFxrTri = FxrTri%nFxrTri; NtTri = FxrTri%NtTri
  nResD = FxrD%nResD; nResTri = FxrTri%nResTri
  ng = IsoData%ng; ntiso = IsoData%ntiso
  ibres = GroupInfo%nofg+1; ieres = ibres-1+GroupInfo%norg
  nGenD = FxrD%nGenD; nCldTri = FxrTri%nCldTri; nAICTri = FxrTri%nAICTri;
  NtRTri = FxrTri%NtRTri
  nMLvC1G = mlgdata0%c_nmaclv1G; nMLvF1G = mlgdata0%f_nmaclv1G; nMLvFMg = mlgdata0%f_nmaclv
  !print*, 'Fxr Sorting Info', nGenD, nGenTri, nCldD, nCldTri
  IF (caseEff .EQ. EffAll) THEN
    IF (nFxrD.NE.0) THEN
      ALLOCATE(cuBlockFxr%mapglobalidD(nFxrD),cuBlockFxr%mapG2RD(nFxrD))
      ALLOCATE(cuBlockFxr%pinidD(nFxrD), cuBlockFxr%nisoD(nFxrD))
      ALLOCATE(cuBlockFxr%pnumD(ntiso,nFxrD), cuBlockFxr%indD(ntiso,nFxrD))
      ALLOCATE(cuBlockFxr%MapNuclD(ntiso,nFxrD))
      ALLOCATE(cuBlockFxr%itTempD(nFxrD,ntiso),cuBlockFxr%wtTempD(nFxrD,ntiso))
      istat = cudaMemCpy(cuBlockFxr%mapG2RD,FxrD%mapG2R,nFxrD)
      istat = cudaMemCpy(cuBlockFxr%mapglobalidD,FxrD%mapglobalid,nFxrD)
      istat = cudaMemCpy(cuBlockFxr%nisoD,FxrD%niso,nFxrD)
      istat = cudaMemCpy(cuBlockFxr%pinidD,FxrD%ipin,nFxrD)
      istat = cudaMemCpy(cuBlockFxr%pnumD,FxrD%pnum,ntiso*nFxrD)
      istat = cudaMemCpy(cuBlockFxr%indD,FxrD%ind,ntiso*nFxrD)
      istat = cudaMemCpy(cuBlockFxr%MapNuclD,FxrD%MapNucl,ntiso*nFxrD)
      istat = cudaMemCpy(cuBlockFxr%itTempD,FxrD%itTemp,ntiso*nFxrD)
      istat = cudaMemCpy(cuBlockFxr%wtTempD,FxrD%wtTemp,ntiso*nFxrD)
    END IF
    IF (nFxrTri.NE.0) THEN
      ALLOCATE(cuBlockFxr%AccNT(nFxrTri+1))
      ALLOCATE(cuBlockFxr%mapglobalidT(nFxrTri),cuBlockFxr%mapG2RT(nFxrTri))
      ALLOCATE(cuBlockFxr%pinidT(nFxrTri), cuBlockFxr%tempT(nFxrTri))
      ALLOCATE(cuBlockFxr%pnumT(NtTri), cuBlockFxr%nisoT(nFxrTri))
      ALLOCATE(cuBlockFxr%MapnuclT(NtTri))
      ALLOCATE(cuBlockFxr%itTempT(nFxrTri,FxrTri%nisomax))
      ALLOCATE(cuBlockFxr%wtTempT(nFxrTri,FxrTri%nisomax))
      istat = cudaMemCpy(cuBlockFxr%AccNT,FxrTri%AccNt,nFxrTri+1)
      istat = cudaMemCpy(cuBlockFxr%mapglobalidT,FxrTri%mapglobalid,nFxrTri)
      istat = cudaMemCpy(cuBlockFxr%mapG2RT,FxrTri%mapG2R,nFxrTri)
      istat = cudaMemCpy(cuBlockFxr%pinidT,FxrTri%ipin,nFxrTri)
      istat = cudaMemCpy(cuBlockFxr%tempT,FxrTri%temp,nFxrTri)
      istat = cudaMemCpy(cuBlockFxr%nisoT,FxrTri%niso,nFxrTri)
      istat = cudaMemCpy(cuBlockFxr%pnumT,FxrTri%pnum,NtTri)
      istat = cudaMemCpy(cuBlockFxr%MapnuclT,FxrTri%MapNucl,NtTri)
      istat = cudaMemCpy(cuBlockFxr%itTempT,FxrTri%itTemp,nFxrTri*FxrTri%nisomax)
      istat = cudaMemCpy(cuBlockFxr%wtTempT,FxrTri%wtTemp,nFxrTri*FxrTri%nisomax)
    END IF
    IF (nResD.GT.0) THEN
      ALLOCATE(cuBlockFxr%mapR2GD(nResD))
      ALLOCATE(cuBlockFxr%FresoAIsoD(ntiso,ibres:ieres,nResD), cuBlockFxr%FresoAD(ibres:ieres,nResD))
      ALLOCATE(cuBlockFxr%FresoFIsoD(ntiso,ibres:ieres,nResD), cuBlockFxr%FresoF(ibres:ieres,nResD))
      ALLOCATE(cuBlockFxr%FresoNF(ibres:ieres,nResD))
      istat = cudaMemCpy(cuBlockFxr%mapR2GD,FxrD%mapR2G,nResD)
      istat = cudaMemCpy(cuBlockFxr%FresoAIsoD,Fxrd%fresoAIso,ntiso*(ieres-ibres+1)*nResD)
      istat = cudaMemCpy(cuBlockFxr%FresoFIsoD,Fxrd%fresoFIso,ntiso*(ieres-ibres+1)*nResD)
      istat = cudaMemCpy(cuBlockFxr%FresoAD,Fxrd%fresoa,(ieres-ibres+1)*nResD)
      istat = cudaMemCpy(cuBlockFxr%FresoF,Fxrd%fresokf,(ieres-ibres+1)*nResD)
      istat = cudaMemCpy(cuBlockFxr%FresoNF,Fxrd%fresonf,(ieres-ibres+1)*nResD)
    END IF
    IF (nResTri.GT.0) THEN
      ALLOCATE(cuBlockFxr%ACCNtR(nResTri+1), cuBlockFxr%mapR2GT(nResTri))
      ALLOCATE(cuBlockFxr%FresoAIsoT(NtRTri*(ieres-ibres+1)),cuBlockFxr%FresoAT(ibres:ieres,nResTri))
      istat = cudaMemCpy(cuBlockFxr%ACCNtR,FxrTri%AccNtR,nResTri+1)
      istat = cudaMemCpy(cuBlockFxr%mapR2GT,FxrTri%mapR2G,nResTri)
      istat = cudaMemCpy(cuBlockFxr%FresoAIsoT,FxrTri%fresoAIso,NtRTri*(ieres-ibres+1))
      istat = cudaMemCpy(cuBlockFxr%FresoAT,FxrTri%fresoa,(ieres-ibres+1)*nResTri)
    END IF
    IF (nGenD.GT.0) THEN
      ALLOCATE(cuBlockFxr%mapGR2GD(nGenD))
      ALLOCATE(cuBlockFxr%xseq_f_mgD(nMLvFMg,ibres:ieres,nGenD))
      ALLOCATE(cuBlockFxr%FnAdjD(ibres:ieres,nGenD))
      ALLOCATE(cuBlockFxr%FtAdjD(nMLvFMg,ibres:ieres,nGenD))
      istat = cudaMemCpy(cuBlockFxr%mapGR2GD,FxrD%mapGR2G,nGenD)
      istat = cudaMemCpy(cuBlockFxr%xseq_f_mgD,FxrD%xseq_f_mg,nMLvFMg*(ieres-ibres+1)*nGenD)
      istat = cudaMemcpy(cuBlockFxr%FnAdjD,FxrD%FnAdj,(ieres-ibres+1)*nGenD)
      istat = cudaMemCpy(cuBlockFxr%FtAdjD,FxrD%Ftadj,nMLvFMg*(ieres-ibres+1)*nGenD)
    END IF
    IF (nCldTri.GT.0) THEN
      ALLOCATE(cuBlockFxr%mapC2GT(nCldTri))
      ALLOCATE(cuBlockFxr%xseq_c_1gT(nMLvC1G,nCldTri))
      istat = cudaMemCpy(cuBlockFxr%mapC2GT,FxrTri%mapC2G,nCldTri)
      istat = cudaMemCpy(cuBlockFxr%xseq_c_1gT,FxrTri%xseq_c_1g,nMLvC1G*nCldTri)
    END IF
    IF (nAICTri.GT.0) THEN
      ALLOCATE(cuBlockFxr%mapA2GT(nAICTri))
      ALLOCATE(cuBlockFxr%xseq_f_1gT(nMLvF1G,nAICTri))
      istat = cudaMemCpy(cuBlockFxr%mapA2GT,FxrTri%mapA2G,nAICTri)
      istat = cudaMemCpy(cuBlockFxr%xseq_f_1gT,FxrTri%xseq_f_1g,nMLvF1G*nAICTri)
    END IF
  ELSE IF (caseEff.EQ.EffFxrD) THEN
    IF (nFxrD.NE.0) THEN
      ALLOCATE(cuBlockFxr%mapglobalidD(nFxrD),cuBlockFxr%mapG2RD(nFxrD))
      ALLOCATE(cuBlockFxr%pinidD(nFxrD), cuBlockFxr%nisoD(nFxrD))
      ALLOCATE(cuBlockFxr%pnumD(ntiso,nFxrD), cuBlockFxr%indD(ntiso,nFxrD))
      ALLOCATE(cuBlockFxr%MapNuclD(ntiso,nFxrD))
      ALLOCATE(cuBlockFxr%wtTempD(nFxrD,ntiso),cuBlockFxr%itTempD(nFxrD,ntiso))
      istat = cudaMemCpy(cuBlockFxr%mapG2RD,FxrD%mapG2R,nFxrD)
      istat = cudaMemCpy(cuBlockFxr%mapglobalidD,FxrD%mapglobalid,nFxrD)
      istat = cudaMemCpy(cuBlockFxr%nisoD,FxrD%niso,nFxrD)
      istat = cudaMemCpy(cuBlockFxr%pinidD,FxrD%ipin,nFxrD)
      istat = cudaMemCpy(cuBlockFxr%pnumD,FxrD%pnum,ntiso*nFxrD)
      istat = cudaMemCpy(cuBlockFxr%indD,FxrD%ind,ntiso*nFxrD)
      istat = cudaMemCpy(cuBlockFxr%MapNuclD,FxrD%MapNucl,ntiso*nFxrD)
      istat = cudaMemCpy(cuBlockFxr%itTempD,FxrD%itTemp,ntiso*nFxrD)
      istat = cudaMemCpy(cuBlockFxr%wtTempD,FxrD%wtTemp,ntiso*nFxrD)
    END IF
    IF (nResD.GT.0) THEN
      ALLOCATE(cuBlockFxr%mapR2GD(nResD))
      ALLOCATE(cuBlockFxr%FresoAIsoD(ntiso,ibres:ieres,nResD), cuBlockFxr%FresoAD(ibres:ieres,nResD))
      ALLOCATE(cuBlockFxr%FresoFIsoD(ntiso,ibres:ieres,nResD), cuBlockFxr%FresoF(ibres:ieres,nResD))
      ALLOCATE(cuBlockFxr%FresoNF(ibres:ieres,nResD))
      istat = cudaMemCpy(cuBlockFxr%mapR2GD,FxrD%mapR2G,nResD)
      istat = cudaMemCpy(cuBlockFxr%FresoAIsoD,Fxrd%fresoAIso,ntiso*(ieres-ibres+1)*nResD)
      istat = cudaMemCpy(cuBlockFxr%FresoFIsoD,Fxrd%fresoFIso,ntiso*(ieres-ibres+1)*nResD)
      istat = cudaMemCpy(cuBlockFxr%FresoAD,Fxrd%fresoa,(ieres-ibres+1)*nResD)
      istat = cudaMemCpy(cuBlockFxr%FresoF,Fxrd%fresokf,(ieres-ibres+1)*nResD)
      istat = cudaMemCpy(cuBlockFxr%FresoNF,Fxrd%fresonf,(ieres-ibres+1)*nResD)
    END IF
    IF (nGenD.GT.0) THEN
      ALLOCATE(cuBlockFxr%mapGR2GD(nGenD))
      ALLOCATE(cuBlockFxr%xseq_f_mgD(nMLvFMg,ibres:ieres,nGenD))
      ALLOCATE(cuBlockFxr%FnAdjD(ibres:ieres,nGenD))
      ALLOCATE(cuBlockFxr%FtAdjD(nMLvFMg,ibres:ieres,nGenD))
      istat = cudaMemCpy(cuBlockFxr%mapGR2GD,FxrD%mapGR2G,nGenD)
      istat = cudaMemCpy(cuBlockFxr%xseq_f_mgD,FxrD%xseq_f_mg,nMLvFMg*(ieres-ibres+1)*nGenD)
      istat = cudaMemcpy(cuBlockFxr%FnAdjD,FxrD%FnAdj,(ieres-ibres+1)*nGenD)
      istat = cudaMemCpy(cuBlockFxr%FtAdjD,FxrD%Ftadj,nMLvFMg*(ieres-ibres+1)*nGenD)
    END IF
  ELSE IF (caseEff.EQ.EffBlkD) THEN
    IF (nFxrD.NE.0) THEN
      ALLOCATE(cuBlockFxr%mapglobalidD(nFxrD),cuBlockFxr%mapG2RD(nFxrD))
      ALLOCATE(cuBlockFxr%pinidD(nFxrD), cuBlockFxr%nisoD(nFxrD))
      ALLOCATE(cuBlockFxr%pnumD(ntiso,nFxrD), cuBlockFxr%indD(ntiso,nFxrD))
      ALLOCATE(cuBlockFxr%MapNuclD(ntiso,nFxrD))
      ALLOCATE(cuBlockFxr%wtTempD(nFxrD,ntiso),cuBlockFxr%itTempD(nFxrD,ntiso))
      istat = cudaMemCpy(cuBlockFxr%mapG2RD,FxrD%mapG2R,nFxrD)
      istat = cudaMemCpy(cuBlockFxr%mapglobalidD,FxrD%mapglobalid,nFxrD)
      istat = cudaMemCpy(cuBlockFxr%nisoD,FxrD%niso,nFxrD)
      istat = cudaMemCpy(cuBlockFxr%pinidD,FxrD%ipin,nFxrD)
      istat = cudaMemCpy(cuBlockFxr%pnumD,FxrD%pnum,ntiso*nFxrD)
      istat = cudaMemCpy(cuBlockFxr%indD,FxrD%ind,ntiso*nFxrD)
      istat = cudaMemCpy(cuBlockFxr%MapNuclD,FxrD%MapNucl,ntiso*nFxrD)
      istat = cudaMemCpy(cuBlockFxr%itTempD,FxrD%itTemp,ntiso*nFxrD)
      istat = cudaMemCpy(cuBlockFxr%wtTempD,FxrD%wtTemp,ntiso*nFxrD)
    END IF
    IF (nResD.GT.0) THEN
      ALLOCATE(cuBlockFxr%mapR2GD(nResD))
      istat = cudaMemCpy(cuBlockFxr%mapR2GD,FxrD%mapR2G,nResD)
    END IF
    IF (nGenD.GT.0) THEN
      ALLOCATE(cuBlockFxr%mapGR2GD(nGenD))
      ALLOCATE(cuBlockFxr%xseq_f_mgD(nMLvFMg,ibres:ieres,nGenD))
      ALLOCATE(cuBlockFxr%FnAdjD(ibres:ieres,nGenD))
      ALLOCATE(cuBlockFxr%FtAdjD(nMLvFMg,ibres:ieres,nGenD))
      istat = cudaMemCpy(cuBlockFxr%mapGR2GD,FxrD%mapGR2G,nGenD)
      istat = cudaMemCpy(cuBlockFxr%xseq_f_mgD,FxrD%xseq_f_mg,nMLvFMg*(ieres-ibres+1)*nGenD)
      istat = cudaMemcpy(cuBlockFxr%FnAdjD,FxrD%FnAdj,(ieres-ibres+1)*nGenD)
      istat = cudaMemCpy(cuBlockFxr%FtAdjD,FxrD%Ftadj,nMLvFMg*(ieres-ibres+1)*nGenD)
    END IF
  ELSE IF (caseEff.EQ.EffFxrT) THEN
    IF (nFxrTri.NE.0) THEN
      ALLOCATE(cuBlockFxr%AccNT(nFxrTri+1))
      ALLOCATE(cuBlockFxr%mapglobalidT(nFxrTri),cuBlockFxr%mapG2RT(nFxrTri))
      ALLOCATE(cuBlockFxr%pinidT(nFxrTri), cuBlockFxr%tempT(nFxrTri))
      ALLOCATE(cuBlockFxr%pnumT(NtTri), cuBlockFxr%nisoT(nFxrTri))
      ALLOCATE(cuBlockFxr%MapnuclT(NtTri))
      ALLOCATE(cuBlockFxr%wtTempT(nFxrTri,FxrTri%nisomax))
      ALLOCATE(cuBlockFxr%itTempT(nFxrTri,FxrTri%nisomax))
      istat = cudaMemCpy(cuBlockFxr%AccNT,FxrTri%AccNt,nFxrTri+1)
      istat = cudaMemCpy(cuBlockFxr%mapglobalidT,FxrTri%mapglobalid,nFxrTri)
      istat = cudaMemCpy(cuBlockFxr%mapG2RT,FxrTri%mapG2R,nFxrTri)
      istat = cudaMemCpy(cuBlockFxr%pinidT,FxrTri%ipin,nFxrTri)
      istat = cudaMemCpy(cuBlockFxr%tempT,FxrTri%temp,nFxrTri)
      istat = cudaMemCpy(cuBlockFxr%nisoT,FxrTri%niso,nFxrTri)
      istat = cudaMemCpy(cuBlockFxr%pnumT,FxrTri%pnum,NtTri)
      istat = cudaMemCpy(cuBlockFxr%MapnuclT,FxrTri%MapNucl,NtTri)
      istat = cudaMemCpy(cuBlockFxr%itTempT,FxrTri%itTemp,nFxrTri*FxrTri%nisomax)
      istat = cudaMemCpy(cuBlockFxr%wtTempT,FxrTri%wtTemp,nFxrTri*FxrTri%nisomax)
    END IF
    IF (nResTri.GT.0) THEN
      ALLOCATE(cuBlockFxr%ACCNtR(nResTri+1), cuBlockFxr%mapR2GT(nResTri))
      ALLOCATE(cuBlockFxr%FresoAIsoT(NtRTri*(ieres-ibres+1)),cuBlockFxr%FresoAT(ibres:ieres,nResTri))
      istat = cudaMemCpy(cuBlockFxr%ACCNtR,FxrTri%AccNtR,nResTri+1)
      istat = cudaMemCpy(cuBlockFxr%mapR2GT,FxrTri%mapR2G,nResTri)
      istat = cudaMemCpy(cuBlockFxr%FresoAIsoT,FxrTri%fresoAIso,NtRTri*(ieres-ibres+1))
      istat = cudaMemCpy(cuBlockFxr%FresoAT,FxrTri%fresoa,(ieres-ibres+1)*nResTri)
    END IF
    IF (nCldTri.GT.0) THEN
      ALLOCATE(cuBlockFxr%mapC2GT(nCldTri))
      ALLOCATE(cuBlockFxr%xseq_c_1gT(nMLvC1G,nCldTri))
      istat = cudaMemCpy(cuBlockFxr%mapC2GT,FxrTri%mapC2G,nCldTri)
      istat = cudaMemCpy(cuBlockFxr%xseq_c_1gT,FxrTri%xseq_c_1g,nMLvC1G*nCldTri)
    END IF
    IF (nAICTri.GT.0) THEN
      ALLOCATE(cuBlockFxr%mapA2GT(nAICTri))
      ALLOCATE(cuBlockFxr%xseq_f_1gT(nMLvF1G,nAICTri))
      istat = cudaMemCpy(cuBlockFxr%mapA2GT,FxrTri%mapA2G,nAICTri)
      istat = cudaMemCpy(cuBlockFxr%xseq_f_1gT,FxrTri%xseq_f_1g,nMLvF1G*nAICTri)
    END IF
  END IF

END IF
END SUBROUTINE

SUBROUTINE CopyInCuBlockFxr
! Deprecated, no use
IMPLICIT NONE
!
!cuBlockFxr%pnumD = Fxrd%pnum;
!cuBlockFxr%FresoAD = Fxrd%fresoa; cuBlockFxr%FresoF = Fxrd%fresof; cuBlockFxr%FresoNF = Fxrd%fresonf;
!
!cuBlockFxr%pnumT = FxrTri%pnum;
!cuBlockFxr%FresoAT = FxrTri%fresoa;

END SUBROUTINE

SUBROUTINE DestroyCuBlockFxr(caseCon, caseEff)
USE Core_mod, ONLY : GroupInfo
USE XSLIB_MOD,  ONLY : mlgdata0
!USE Core_mod, ONLY : GroupInfo
IMPLICIT NONE
INTEGER :: caseCon, caseEff
INTEGER :: nFxrD, nFxrTri, NtTri, nResTri, nResD, ntfsrD, ntfsrT, ng, ntiso, ibres, ieres
INTEGER :: nGenD, nCldD, nGenTri, nCldTri, nAICTri, NtRTri, nMLvC1G, nMLvF1G, nMLvFMg

INTEGER :: istat
!INTEGER :: ig, ilv, i, j
IF (caseCon .EQ. ConMac) THEN
  IF (FxrD%nFxrD.GT.0) THEN
    DEALLOCATE(cuBlockFxr%mapglobalidD)
    DEALLOCATE(cuBlockFxr%mapR2GD,cuBlockFxr%mapG2RD)
    DEALLOCATE(cuBlockFxr%pinidD,cuBlockFxr%lvoidD)
    DEALLOCATE(cuBlockFxr%fsrstD,cuBlockFxr%nfsrD)
    DEALLOCATE(cuBlockFxr%mapfsrD)
    DEALLOCATE(cuBlockFxr%chi)
    DEALLOCATE(cuBlockFxr%itTempD,cuBlockFxr%wtTempD)
  END IF
  IF (FxrTri%nFxrTri.GT.0) THEN
    DEALLOCATE(cuBlockFxr%mapglobalidT)
    DEALLOCATE(cuBlockFxr%mapR2GT,cuBlockFxr%mapG2RT)
    DEALLOCATE(cuBlockFxr%pinidT,cuBlockFxr%lvoidT)
    DEALLOCATE(cuBlockFxr%fsrstT,cuBlockFxr%nfsrT)
    DEALLOCATE(cuBlockFxr%mapfsrT)
    DEALLOCATE(cuBlockFxr%itTempT,cuBlockFxr%wtTempT)
  END IF
  IF (FxrD%nFxrD.GT.0) THEN
    DEALLOCATE(cuBlockFxr%nisoD)
    DEALLOCATE(cuBlockFxr%pnumD)
    DEALLOCATE(cuBlockFxr%MapNuclD)
    IF (FxrD%nResD.GT.0) THEN
      DEALLOCATE(cuBlockFxr%FresoF, cuBlockFxr%FresoAD, cuBlockFxr%FresoNF)
    END IF
  END IF
  IF (FxrTri%nFxrTri.GT.0) THEN
    DEALLOCATE(cuBlockFxr%nisoT,cuBlockFxr%AccNT)
    DEALLOCATE(cuBlockFxr%pnumT)
    DEALLOCATE(cuBlockFxr%MapNuclT)
    IF (FxrTri%nResTri.GT.0) THEN
      DEALLOCATE(cuBlockFxr%FresoAT)
    END IF
  END IF
ELSE IF (caseCon .EQ. ConFSP) THEN
  nFxrD = FxrD%nFxrD; nFxrTri = FxrTri%nFxrTri
  nResD = FxrD%nResD; nResTri = FxrTri%nResTri
  ibres = GroupInfo%nofg+1; ieres = ibres-1+GroupInfo%norg
  nGenD = FxrD%nGenD; nCldTri = FxrTri%nCldTri; nAICTri = FxrTri%nAICTri;
  nMLvC1G = mlgdata0%c_nmaclv1G; nMLvF1G = mlgdata0%f_nmaclv1G; nMLvFMg = mlgdata0%f_nmaclv
  IF (nFxrD.GT.0) THEN
    DEALLOCATE(cuBlockFxr%nisoD)
    DEALLOCATE(cuBlockFxr%mapglobalidD)
    DEALLOCATE(cuBlockFxr%pinidD)
    DEALLOCATE(cuBlockFxr%fsrstD, cuBlockFxr%mapfsrD)
    DEALLOCATE(cuBlockFxr%nfsrD)
    DEALLOCATE(cuBlockFxr%tempD)
    DEALLOCATE(cuBlockFxr%areaD)
    DEALLOCATE(cuBlockFxr%indD, cuBlockFxr%pnumD)
    DEALLOCATE(cuBlockFxr%MapNuclD)
    DEALLOCATE(cuBlockFxr%fsrvolD)
  END IF
  IF (nFxrTri.GT.0) THEN
    DEALLOCATE(cuBlockFxr%nisoT, cuBlockFxr%ACCNtR, cuBlockFxr%AccNT)
    DEALLOCATE(cuBlockFxr%mapglobalidT)
    DEALLOCATE(cuBlockFxr%pinidT)
    DEALLOCATE(cuBlockFxr%fsrstT,cuBlockFxr%mapfsrT)
    DEALLOCATE(cuBlockFxr%nfsrT)
    DEALLOCATE(cuBlockFxr%tempT)
    DEALLOCATE(cuBlockFxr%areaT)
    DEALLOCATE(cuBlockFxr%pnumT)
    DEALLOCATE(cuBlockFxr%MapnuclT)
    DEALLOCATE(cuBlockFxr%fsrvolT)
  END IF
  IF (nGenD.GT.0) THEN
!    FxrD%xseq_f_mg = cuBlockFxr%xseq_f_mgD
!    FxrD%FnAdj = cuBlockFxr%FnAdjD; FxrD%Ftadj = cuBlockFxr%FtAdjD
    istat = cudaMemCpy(FxrD%xseq_f_mg,cuBlockFxr%xseq_f_mgD,nMLvFMg*(ieres-ibres+1)*nGenD)
    istat = cudaMemCpy(FxrD%FnAdj,cuBlockFxr%FnAdjD,(ieres-ibres+1)*nGenD)
    istat = cudaMemCpy(FxrD%Ftadj,cuBlockFxr%FtAdjD,nMLvFMg*(ieres-ibres+1)*nGenD)
    DEALLOCATE(cuBlockFxr%mapGR2GD)
    DEALLOCATE(cuBlockFxr%xseq_f_mgD)
    DEALLOCATE(cuBlockFxr%FtAdjD,cuBlockFxr%FnAdjD)
  END IF
  IF (nCldD.GT.0) THEN
!    FxrD%xseq_c_1g = cuBlockFxr%xseq_c_1gD
    istat = cudaMemCpy(FxrD%xseq_c_1g,cuBlockFxr%xseq_c_1gD,nMLvC1G*nCldD)
    DEALLOCATE(cuBlockFxr%mapC2GD)
    DEALLOCATE(cuBlockFxr%xseq_c_1gD)
  END IF
  IF (nGenTri.GT.0) THEN
!    FxrTri%xseq_f_mg = cuBlockFxr%xseq_f_mgT
!    FxrTri%FnAdj = cuBlockFxr%FnAdjT; FxrTri%Ftadj = cuBlockFxr%FtAdjT
    istat = cudaMemCpy(FxrTri%xseq_f_mg,cuBlockFxr%xseq_f_mgT,nMLvFMg*(ieres-ibres+1)*nGenTri)
    istat = cudaMemCpy(FxrTri%FnAdj,cuBlockFxr%FnAdjT,(ieres-ibres+1)*nGenTri)
    istat = cudaMemCpy(FxrTri%Ftadj,cuBlockFxr%FtAdjT,nMLvFMg*(ieres-ibres+1)*nGenTri)
    DEALLOCATE(cuBlockFxr%mapGR2GT)
    DEALLOCATE(cuBlockFxr%xseq_f_mgT)
    DEALLOCATE(cuBlockFxr%FtAdjT,cuBlockFxr%FnAdjT)
  END IF
  IF (nCldTri.GT.0) THEN
!    FxrTri%xseq_c_1g = cuBlockFxr%xseq_c_1gT
    istat = cudaMemCpy(FxrTri%xseq_c_1g,cuBlockFxr%xseq_c_1gT,nMLvC1G*nCldTri)
    DEALLOCATE(cuBlockFxr%mapC2GT)
    DEALLOCATE(cuBlockFxr%xseq_c_1gT)
  END IF
  IF (nAICTri.GT.0) THEN
!    FxrTri%xseq_f_1g = cuBlockFxr%xseq_f_1gT
    istat = cudaMemCpy(FxrTri%xseq_f_1g,cuBlockFxr%xseq_f_1gT,nMLvF1G*nAICTri)
    DEALLOCATE(cuBlockFxr%mapA2GT)
    DEALLOCATE(cuBlockFxr%xseq_f_1gT)
  END IF
ELSE IF (caseCon.EQ.ConEff) THEN
  nFxrD = FxrD%nFxrD; nFxrTri = FxrTri%nFxrTri
  nResD = FxrD%nResD; nResTri = FxrTri%nResTri
  ntiso = IsoData%ntiso
  ibres = GroupInfo%nofg+1; ieres = ibres-1+GroupInfo%norg
  nGenD = FxrD%nGenD; nCldTri = FxrTri%nCldTri; nAICTri = FxrTri%nAICTri;
  NtRTri = FxrTri%NtRTri
  IF (caseEff.EQ.EffAll) THEN
    IF (nFxrD.NE.0) THEN
      DEALLOCATE(cuBlockFxr%mapglobalidD,cuBlockFxr%mapG2RD)
      DEALLOCATE(cuBlockFxr%pinidD, cuBlockFxr%nisoD)
      DEALLOCATE(cuBlockFxr%pnumD, cuBlockFxr%indD)
      DEALLOCATE(cuBlockFxr%MapNuclD)
      DEALLOCATE(cuBlockFxr%itTempD,cuBlockFxr%wtTempD)
    END IF
    IF (nFxrTri.NE.0) THEN
      DEALLOCATE(cuBlockFxr%AccNT)
      DEALLOCATE(cuBlockFxr%mapglobalidT,cuBlockFxr%mapG2RT)
      DEALLOCATE(cuBlockFxr%pinidT, cuBlockFxr%tempT)
      DEALLOCATE(cuBlockFxr%pnumT, cuBlockFxr%nisoT)
      DEALLOCATE(cuBlockFxr%MapnuclT)
      DEALLOCATE(cuBlockFxr%itTempT,cuBlockFxr%wtTempT)
    END IF
    IF (nResD.GT.0) THEN
      istat = cudaMemCpy(Fxrd%fresoAIso,cuBlockFxr%FresoAIsoD,ntiso*(ieres-ibres+1)*nResD)
      istat = cudaMemCpy(Fxrd%fresoFIso,cuBlockFxr%FresoFIsoD,ntiso*(ieres-ibres+1)*nResD)
      istat = cudaMemCpy(Fxrd%fresoa,cuBlockFxr%FresoAD,(ieres-ibres+1)*nResD)
      istat = cudaMemCpy(Fxrd%fresokf,cuBlockFxr%FresoF,(ieres-ibres+1)*nResD)
      istat = cudaMemCpy(Fxrd%fresonf,cuBlockFxr%FresoNF,(ieres-ibres+1)*nResD)
      DEALLOCATE(cuBlockFxr%mapR2GD)
      DEALLOCATE(cuBlockFxr%FresoAIsoD,cuBlockFxr%FresoAD)
      DEALLOCATE(cuBlockFxr%FresoFIsoD,cuBlockFxr%FresoF)
      DEALLOCATE(cuBlockFxr%FresoNF)
    END IF
    IF (nResTri.GT.0) THEN
      istat = cudaMemCpy(FxrTri%fresoAIso,cuBlockFxr%FresoAIsoT,NtRTri*(ieres-ibres+1))
      istat = cudaMemCpy(FxrTri%fresoa,cuBlockFxr%FresoAT,(ieres-ibres+1)*nResTri)
      DEALLOCATE(cuBlockFxr%ACCNtR,cuBlockFxr%mapR2GT)
      DEALLOCATE(cuBlockFxr%FresoAIsoT,cuBlockFxr%FresoAT)
    END IF
    IF (nGenD.GT.0) THEN
      DEALLOCATE(cuBlockFxr%mapGR2GD)
      DEALLOCATE(cuBlockFxr%xseq_f_mgD)
      DEALLOCATE(cuBlockFxr%FnAdjD)
      DEALLOCATE(cuBlockFxr%FtAdjD)
    END IF
    IF (nCldTri.GT.0) THEN
      DEALLOCATE(cuBlockFxr%mapC2GT)
      DEALLOCATE(cuBlockFxr%xseq_c_1gT)
    END IF
    IF (nAICTri.GT.0) THEN
      DEALLOCATE(cuBlockFxr%mapA2GT)
      DEALLOCATE(cuBlockFxr%xseq_f_1gT)
    END IF
  ELSE IF(caseEff.EQ.EffFxrD) THEN
    IF (nFxrD.NE.0) THEN
      DEALLOCATE(cuBlockFxr%mapglobalidD,cuBlockFxr%mapG2RD)
      DEALLOCATE(cuBlockFxr%pinidD, cuBlockFxr%nisoD)
      DEALLOCATE(cuBlockFxr%pnumD, cuBlockFxr%indD)
      DEALLOCATE(cuBlockFxr%MapNuclD)
      DEALLOCATE(cuBlockFxr%itTempD,cuBlockFxr%wtTempD)
    END IF
    IF (nResD.GT.0) THEN
      istat = cudaMemCpy(Fxrd%fresoAIso,cuBlockFxr%FresoAIsoD,ntiso*(ieres-ibres+1)*nResD)
      istat = cudaMemCpy(Fxrd%fresoFIso,cuBlockFxr%FresoFIsoD,ntiso*(ieres-ibres+1)*nResD)
      istat = cudaMemCpy(Fxrd%fresoa,cuBlockFxr%FresoAD,(ieres-ibres+1)*nResD)
      istat = cudaMemCpy(Fxrd%fresokf,cuBlockFxr%FresoF,(ieres-ibres+1)*nResD)
      istat = cudaMemCpy(Fxrd%fresonf,cuBlockFxr%FresoNF,(ieres-ibres+1)*nResD)
      DEALLOCATE(cuBlockFxr%mapR2GD)
      DEALLOCATE(cuBlockFxr%FresoAIsoD,cuBlockFxr%FresoAD)
      DEALLOCATE(cuBlockFxr%FresoFIsoD,cuBlockFxr%FresoF)
      DEALLOCATE(cuBlockFxr%FresoNF)
    END IF
    IF (nGenD.GT.0) THEN
      DEALLOCATE(cuBlockFxr%mapGR2GD)
      DEALLOCATE(cuBlockFxr%xseq_f_mgD)
      DEALLOCATE(cuBlockFxr%FnAdjD)
      DEALLOCATE(cuBlockFxr%FtAdjD)
    END IF
  ELSE IF(caseEff.EQ.EffBlkD) THEN
    IF (nFxrD.NE.0) THEN
      DEALLOCATE(cuBlockFxr%mapglobalidD,cuBlockFxr%mapG2RD)
      DEALLOCATE(cuBlockFxr%pinidD, cuBlockFxr%nisoD)
      DEALLOCATE(cuBlockFxr%pnumD, cuBlockFxr%indD)
      DEALLOCATE(cuBlockFxr%MapNuclD)
      DEALLOCATE(cuBlockFxr%itTempD,cuBlockFxr%wtTempD)
    END IF
    IF (nResD.GT.0) THEN
      !istat = cudaMemCpy(Fxrd%fresoAIso,cuBlockFxr%FresoAIsoD,ntiso*(ieres-ibres+1)*nResD)
      !istat = cudaMemCpy(Fxrd%fresoFIso,cuBlockFxr%FresoFIsoD,ntiso*(ieres-ibres+1)*nResD)
      !istat = cudaMemCpy(Fxrd%fresoa,cuBlockFxr%FresoAD,(ieres-ibres+1)*nResD)
      !istat = cudaMemCpy(Fxrd%fresof,cuBlockFxr%FresoF,(ieres-ibres+1)*nResD)
      !istat = cudaMemCpy(Fxrd%fresonf,cuBlockFxr%FresoNF,(ieres-ibres+1)*nResD)
      DEALLOCATE(cuBlockFxr%mapR2GD)
      !DEALLOCATE(cuBlockFxr%FresoAIsoD,cuBlockFxr%FresoAD)
      !DEALLOCATE(cuBlockFxr%FresoFIsoD,cuBlockFxr%FresoF)
      !DEALLOCATE(cuBlockFxr%FresoNF)
    END IF
    IF (nGenD.GT.0) THEN
      DEALLOCATE(cuBlockFxr%mapGR2GD)
      DEALLOCATE(cuBlockFxr%xseq_f_mgD)
      DEALLOCATE(cuBlockFxr%FnAdjD)
      DEALLOCATE(cuBlockFxr%FtAdjD)
    END IF
  ELSE IF (caseEff.EQ.EffFxrT) THEN
    IF (nFxrTri.NE.0) THEN
      DEALLOCATE(cuBlockFxr%AccNT)
      DEALLOCATE(cuBlockFxr%mapglobalidT,cuBlockFxr%mapG2RT)
      DEALLOCATE(cuBlockFxr%pinidT, cuBlockFxr%tempT)
      DEALLOCATE(cuBlockFxr%pnumT, cuBlockFxr%nisoT)
      DEALLOCATE(cuBlockFxr%MapnuclT)
      DEALLOCATE(cuBlockFxr%itTempT,cuBlockFxr%wtTempT)
    END IF
    IF (nResTri.GT.0) THEN
      istat = cudaMemCpy(FxrTri%fresoAIso,cuBlockFxr%FresoAIsoT,NtRTri*(ieres-ibres+1))
      istat = cudaMemCpy(FxrTri%fresoa,cuBlockFxr%FresoAT,(ieres-ibres+1)*nResTri)
      DEALLOCATE(cuBlockFxr%ACCNtR,cuBlockFxr%mapR2GT)
      DEALLOCATE(cuBlockFxr%FresoAIsoT,cuBlockFxr%FresoAT)
    END IF
    IF (nCldTri.GT.0) THEN
      DEALLOCATE(cuBlockFxr%mapC2GT)
      DEALLOCATE(cuBlockFxr%xseq_c_1gT)
    END IF
    IF (nAICTri.GT.0) THEN
      DEALLOCATE(cuBlockFxr%mapA2GT)
      DEALLOCATE(cuBlockFxr%xseq_f_1gT)
    END IF
  END IF
END IF
END SUBROUTINE

SUBROUTINE ConstructCuResPin(caseid)
USE GEOM, ONLY : core
USE PE_Mod, ONLY : PE
USE Core_mod, ONLY : GroupInfo
USE XSLIB_MOD,  ONLY : mlgdata0, nreshel, nelthel
IMPLICIT NONE
INTEGER :: caseid
INTEGER :: nxyz, ng, ntiso, ibres, ieres
INTEGER :: nMLvC1G, nMLvF1G, nMLvFMg

INTEGER :: istat

IF (caseid.EQ.ConFSP) THEN
  nxyz = core%nxy*(PE%myze-PE%myzb+1)
  ng = IsoData%ng; ntiso = IsoData%ntiso
  ibres = GroupInfo%nofg+1; ieres = ibres-1+GroupInfo%norg
  nMLvC1G = mlgdata0%c_nmaclv1G; nMLvF1G = mlgdata0%f_nmaclv1G; nMLvFMg = mlgdata0%f_nmaclv

  ALLOCATE(cuResPin%lPinRes(nxyz), cuResPin%lCellRes(nxyz), cuResPin%lAIC(nxyz))
  ALLOCATE(cuResPin%avgxseq_1g(nMLvF1G,nxyz),cuResPin%avgxseq_mg(nMLvFMg,ibres:ieres,nxyz))
  ALLOCATE(cuResPin%FnAdj(ibres:ieres,nxyz),cuResPin%temp(nxyz))

!  cuResPin%lPinRes = ResPin%lres; cuResPin%lCellRes = GlobalPinInfo%lres; cuResPin%lAIC = GlobalPinInfo%lAIC
!  cuResPin%avgxseq_1g = ResPin%avgxseq_1g; cuResPin%avgxseq_mg = ResPin%avgxseq_mg
!  cuResPin%FnAdj = ResPin%FnAdj; cuResPin%temp = ResPin%temp

  istat = cudaMemCpy(cuResPin%lPinRes,ResPin%lres,nxyz)
  istat = cudaMemCpy(cuResPin%lCellRes,GlobalPinInfo%lres,nxyz)
  istat = cudaMemCpy(cuResPin%lAIC,GlobalPinInfo%lAIC,nxyz)
  istat = cudaMemCpy(cuResPin%avgxseq_1g,ResPin%avgxseq_1g,nMLvF1G*nxyz)
  istat = cudaMemCpy(cuResPin%avgxseq_mg,ResPin%avgxseq_mg,nMLvFMg*(ieres-ibres+1)*nxyz)
  istat = cudaMemCpy(cuResPin%FnAdj,ResPin%FnAdj,(ieres-ibres+1)*nxyz)
  istat = cudaMemCpy(cuResPin%temp,ResPin%temp,nxyz)
ELSE IF (caseid.EQ.ConEff) THEN
  nxyz = core%nxy*(PE%myze-PE%myzb+1)
  ng = IsoData%ng; ntiso = IsoData%ntiso
  ibres = GroupInfo%nofg+1; ieres = ibres-1+GroupInfo%norg
  nMLvC1G = mlgdata0%c_nmaclv1G; nMLvF1G = mlgdata0%f_nmaclv1G; nMLvFMg = mlgdata0%f_nmaclv

  ALLOCATE(cuResPin%lPinRes(nxyz), cuResPin%lCellRes(nxyz), cuResPin%lAIC(nxyz))
  ALLOCATE(cuResPin%avgxseq_1g(nMLvF1G,nxyz),cuResPin%avgxseq_mg(nMLvFMg,ibres:ieres,nxyz))
  ALLOCATE(cuResPin%FnAdj(ibres:ieres,nxyz),cuResPin%temp(nxyz))
  ALLOCATE(cuResPin%rifa(nreshel,ibres:ieres,nxyz),cuResPin%riff(nreshel,ibres:ieres,nxyz))
  ALLOCATE(cuResPin%ind(nelthel,nxyz),cuResPin%mapnucl(nelthel,nxyz),cuResPin%pnum(nelthel,nxyz))
  ALLOCATE(cuResPin%niso(nxyz))

!  cuResPin%lPinRes = ResPin%lres; cuResPin%lCellRes = GlobalPinInfo%lres; cuResPin%lAIC = GlobalPinInfo%lAIC
!  cuResPin%avgxseq_1g = ResPin%avgxseq_1g; cuResPin%avgxseq_mg = ResPin%avgxseq_mg
!  cuResPin%FnAdj = ResPin%FnAdj; cuResPin%temp = ResPin%temp
!  cuResPin%rifa = ResPin%rifa; cuResPin%riff = ResPin%riff
!  cuResPin%ind = ResPin%ind; cuResPin%mapnucl = ResPin%mapnucl; cuResPin%pnum = ResPin%pnum;
!  cuResPin%niso = ResPin%niso

  istat = cudaMemCpy(cuResPin%lPinRes,ResPin%lres,nxyz)
  istat = cudaMemCpy(cuResPin%lCellRes,GlobalPinInfo%lres,nxyz)
  istat = cudaMemCpy(cuResPin%lAIC,GlobalPinInfo%lAIC,nxyz)
  istat = cudaMemCpy(cuResPin%avgxseq_1g,ResPin%avgxseq_1g,nMLvF1G*nxyz)
  istat = cudaMemCpy(cuResPin%avgxseq_mg,ResPin%avgxseq_mg,nMLvFMg*(ieres-ibres+1)*nxyz)
  istat = cudaMemCpy(cuResPin%FnAdj,ResPin%FnAdj,(ieres-ibres+1)*nxyz)
  istat = cudaMemCpy(cuResPin%temp,ResPin%temp,nxyz)
  istat = cudaMemCpy(cuResPin%rifa,ResPin%rifa,nreshel*(ieres-ibres+1)*nxyz)
  istat = cudaMemCpy(cuResPin%riff,ResPin%riff,nreshel*(ieres-ibres+1)*nxyz)
  istat = cudaMemCpy(cuResPin%ind,ResPin%ind,nelthel*nxyz)
  istat = cudaMemCpy(cuResPin%mapnucl,ResPin%mapnucl,nelthel*nxyz)
  istat = cudaMemCpy(cuResPin%pnum,ResPin%pnum,nelthel*nxyz)
  istat = cudaMemCpy(cuResPin%niso,ResPin%niso,nxyz)
END IF
END SUBROUTINE

SUBROUTINE DestroyCuResPin(caseid, caseEff)
USE Core_mod, ONLY : GroupInfo
USE XSLIB_MOD,  ONLY : mlgdata0
USE Geom, ONLY : core
USE XSLIB_MOD,  ONLY : mlgdata0, nreshel, nelthel
USE PE_Mod, ONLY : PE
IMPLICIT NONE
INTEGER :: caseid, caseEff
INTEGER :: nxyz, ng, ntiso, ibres, ieres
INTEGER :: nMLvC1G, nMLvF1G, nMLvFMg

INTEGER :: istat

IF (caseid.EQ.ConFSP) THEN
  nxyz = core%nxy*(PE%myze-PE%myzb+1)
  ng = IsoData%ng; ntiso = IsoData%ntiso
  ibres = GroupInfo%nofg+1; ieres = ibres-1+GroupInfo%norg
  nMLvC1G = mlgdata0%c_nmaclv1G; nMLvF1G = mlgdata0%f_nmaclv1G; nMLvFMg = mlgdata0%f_nmaclv

!  ResPin%avgxseq_1g = cuResPin%avgxseq_1g
!  ResPin%avgxseq_mg = cuResPin%avgxseq_mg
!  ResPin%FnAdj = cuResPin%FnAdj

  istat = cudaMemCpy(ResPin%avgxseq_1g,cuResPin%avgxseq_1g,nMLvF1G*nxyz)
  istat = cudaMemCpy(ResPin%avgxseq_mg,cuResPin%avgxseq_mg,nMLvFMg*(ieres-ibres+1)*nxyz)
  istat = cudaMemCpy(ResPin%FnAdj,cuResPin%FnAdj,(ieres-ibres+1)*nxyz)

  DEALLOCATE(cuResPin%lPinRes, cuResPin%lCellRes, cuResPin%lAIC)
  DEALLOCATE(cuResPin%avgxseq_1g,cuResPin%avgxseq_mg)
  DEALLOCATE(cuResPin%FnAdj,cuResPin%temp)
ELSE IF (caseid.EQ.ConEff) THEN
  nxyz = core%nxy*(PE%myze-PE%myzb+1)
  ng = IsoData%ng; ntiso = IsoData%ntiso
  ibres = GroupInfo%nofg+1; ieres = ibres-1+GroupInfo%norg
  nMLvC1G = mlgdata0%c_nmaclv1G; nMLvF1G = mlgdata0%f_nmaclv1G; nMLvFMg = mlgdata0%f_nmaclv

!  ResPin%rifa = cuResPin%rifa; ResPin%riff = cuResPin%riff
  IF (caseEff.EQ.EffAll) THEN
    istat = cudaMemCpy(ResPin%rifa,cuResPin%rifa,nreshel*(ieres-ibres+1)*nxyz)
    istat = cudaMemCpy(ResPin%riff,cuResPin%riff,nreshel*(ieres-ibres+1)*nxyz)

    DEALLOCATE(cuResPin%lPinRes, cuResPin%lCellRes, cuResPin%lAIC)
    DEALLOCATE(cuResPin%avgxseq_1g,cuResPin%avgxseq_mg)
    DEALLOCATE(cuResPin%FnAdj,cuResPin%temp)
    DEALLOCATE(cuResPin%rifa,cuResPin%riff)
    DEALLOCATE(cuResPin%ind,cuResPin%mapnucl,cuResPin%pnum)
    DEALLOCATE(cuResPin%niso)
  ELSE IF (caseEff.EQ.EffRIP) THEN
    DEALLOCATE(cuResPin%lPinRes, cuResPin%lCellRes, cuResPin%lAIC)
    DEALLOCATE(cuResPin%avgxseq_1g,cuResPin%avgxseq_mg)
    DEALLOCATE(cuResPin%FnAdj)
    DEALLOCATE(cuResPin%ind,cuResPin%mapnucl,cuResPin%pnum)
    DEALLOCATE(cuResPin%niso)
  ELSE IF (caseEff.EQ.EffClr) THEN
    istat = cudaMemCpy(ResPin%rifa,cuResPin%rifa,nreshel*(ieres-ibres+1)*nxyz)
    istat = cudaMemCpy(ResPin%riff,cuResPin%riff,nreshel*(ieres-ibres+1)*nxyz)
    DEALLOCATE(cuResPin%temp)
    DEALLOCATE(cuResPin%rifa,cuResPin%riff)
  END IF
END IF
END SUBROUTINE

SUBROUTINE ConstructCuPinInfo
USE GEOM, ONLY : core
USE PE_Mod, ONLY : PE
IMPLICIT NONE
INTEGER :: nxyz
INTEGER :: istat
nxyz = core%nxy*(PE%myze-PE%myzb+1)
ALLOCATE(cuPinInfo%FxrStD(nxyz+1))
ALLOCATE(cuPinInfo%FxrStT(nxyz+1))
ALLOCATE(cuPinInfo%GenStD(nxyz+1))
ALLOCATE(cuPinInfo%GenStT(nxyz+1))
ALLOCATE(cuPinInfo%AICStT(nxyz+1))
!cuPinInfo%FxrStD = GlobalPinInfo%FxrStD
!cuPinInfo%FxrStT = GlobalPinInfo%FxrStT
!cuPinInfo%GenStD = GlobalPinInfo%GenStD
!cuPinInfo%GenStT = GlobalPinInfo%GenStT
!cuPinInfo%AICStT = GlobalPinInfo%AICStT
istat = cudaMemCpy(cuPinInfo%FxrStD,GlobalPinInfo%FxrStD,nxyz+1)
istat = cudaMemCpy(cuPinInfo%FxrStT,GlobalPinInfo%FxrStT,nxyz+1)
istat = cudaMemCpy(cuPinInfo%GenStD,GlobalPinInfo%GenStD,nxyz+1)
istat = cudaMemCpy(cuPinInfo%GenStT,GlobalPinInfo%GenStT,nxyz+1)
istat = cudaMemCpy(cuPinInfo%AICStT,GlobalPinInfo%AICStT,nxyz+1)
!print*, GlobalPinInfo%genstd
END SUBROUTINE

SUBROUTINE DestroyCuPinInfo
DEALLOCATE(cuPinInfo%FxrStD)
DEALLOCATE(cuPinInfo%FxrStT)
DEALLOCATE(cuPinInfo%GenStD)
DEALLOCATE(cuPinInfo%GenStT)
DEALLOCATE(cuPinInfo%AICStT)
END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuGetMacXS_D(niso,pnum,mapfxr,a,kf,nf,s,str,McA,MckF,McnF,McS,McStr,itTmp, wtTmp, gb, ge, iz, nFxrD, caseID)
USE CUDAFOR
IMPLICIT NONE
INTEGER, DEVICE :: niso(:),mapfxr(:),itTmp(:,:)
REAL(GPU_PNUM_PRECISION), DIMENSION(:,:), DEVICE ::pnum
REAL(GPU_XS_PRECISION), DIMENSION(:,:), DEVICE :: a,kf,nf,s,str,McA,MckF,McnF,McS,McStr,wtTmp
INTEGER, VALUE :: gb, ge, iz, nFxrD, caseID

!INTEGER :: xtid, ytid
INTEGER :: ig, jg, ifxr, jfxr
INTEGER :: iso, miso
INTEGER :: jt
REAL(8) :: sigbuf, wtN

REAL(8) :: ND, micA, micnf, mickf, mics, micstr
LOGICAL :: lskipit2
!REAL(8), SHARED :: pnumLoc(XS_BLOCK_GDIM,XS_BLOCK_RDIM)

!xtid = threadIdx%x;
jg = threadIdx%x+XS_BLOCK_GDIM*(blockIdx%x-1)
ig = jg+gb-1
IF (ig .GT. ge) RETURN

!ytid = threadIdx%y
ifxr  = threadIdx%y+XS_BLOCK_RDIM*(blockIdx%y-1)
if (ifxr .GT. nFxrD) RETURN

jfxr = mapfxr(ifxr)-(iz-1)*nFxr
IF (jfxr .LT. 1) RETURN
IF (jfxr .GT. nFxr) RETURN

micA = 0.; micnf = 0.; mickf = 0.; mics = 0.; micstr = 0.

miso = niso(ifxr);

IF (caseID.EQ.MacBase) THEN
  DO iso = 1, miso
    ND = pnum(iso,ifxr)
    jt = itTmp(ifxr,iso); wtN = wtTmp(ifxr,iso)
    lskipit2 = (wtN .GT. 0.99999)

    ! Left Point
    wtN = ND*wtN;

    sigbuf = a(ig, jt)
    micA = micA+wtN*sigbuf
    sigbuf = kf(ig, jt)
    mickf = mickf+wtN*sigbuf
    sigbuf = s(ig, jt)
    mics = mics+wtN*sigbuf
    sigbuf = str(ig, jt)
    micstr = micstr+wtN*sigbuf

    IF (lskipit2) CYCLE
    ! Right Point
    wtN = ND - wtN;
    jt = jt+1

    sigbuf = a(ig, jt)
    micA = micA+wtN*sigbuf
    sigbuf = kf(ig, jt)
    mickf = mickf+wtN*sigbuf
    sigbuf = s(ig, jt)
    mics = mics+wtN*sigbuf
    sigbuf = str(ig, jt)
    micstr = micstr+wtN*sigbuf
  END DO
END IF

IF (caseID.EQ.MacNF) THEN
  DO iso = 1, miso
    ND = pnum(iso,ifxr)
    jt = itTmp(ifxr,iso); wtN = wtTmp(ifxr,iso)
    lskipit2 = (wtN .GT. 0.99999)

    ! Left Point
    wtN = ND*wtN;

    sigbuf = nf(ig, jt)
    micnf = micnf+wtN*sigbuf

    IF (lskipit2) CYCLE
    ! Right Point
    wtN = ND - wtN;
    jt = jt+1

    sigbuf = nf(ig, jt)
    micnf = micnf+wtN*sigbuf
  END DO
END IF

IF (caseID.EQ.MacAll) THEN
  DO iso = 1, miso
    ND = pnum(iso,ifxr)
    jt = itTmp(ifxr,iso); wtN = wtTmp(ifxr,iso)
    lskipit2 = (wtN .GT. 0.99999)

    ! Left Point
    wtN = ND*wtN;

    sigbuf = a(ig, jt)
    micA = micA+wtN*sigbuf
    sigbuf = kf(ig, jt)
    mickf = mickf+wtN*sigbuf
    sigbuf = s(ig, jt)
    mics = mics+wtN*sigbuf
    sigbuf = str(ig, jt)
    micstr = micstr+wtN*sigbuf
    sigbuf = nf(ig, jt)
    micnf = micnf+wtN*sigbuf

    IF (lskipit2) CYCLE
    ! Right Point
    wtN = ND - wtN;
    jt = jt+1

    sigbuf = a(ig, jt)
    micA = micA+wtN*sigbuf
    sigbuf = kf(ig, jt)
    mickf = mickf+wtN*sigbuf
    sigbuf = s(ig, jt)
    mics = mics+wtN*sigbuf
    sigbuf = str(ig, jt)
    micstr = micstr+wtN*sigbuf
    sigbuf = nf(ig, jt)
    micnf = micnf+wtN*sigbuf
  END DO
END IF
!ifxr = mapfxr(ifxr)
IF (caseID.NE.MacNF) THEN
  Mca(jg,jfxr) = micA
  Mckf(jg,jfxr) = mickf
  McS(jg,jfxr) = mics
  McStr(jg,jfxr) = micstr
END IF
IF (caseID.NE.MacBase) Mcnf(jg,jfxr) = micnf
END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuGetMacXS_T(AccNt,pnum,mapfxr,a,s,str,McA,McS,McStr,itTmp, wtTmp, gb, ge, iz, nFxrTri)
USE CUDAFOR
IMPLICIT NONE
INTEGER, DEVICE :: AccNt(:),mapfxr(:),itTmp(:,:)
REAL(GPU_PNUM_PRECISION), DEVICE :: pnum(:)
REAL(GPU_XS_PRECISION), DIMENSION(:,:), DEVICE :: a,s,str,McA,McS,McStr,wtTmp
INTEGER, VALUE :: gb, ge, iz, nFxrTri

!INTEGER :: xtid, ytid
INTEGER :: ig, jg, ifxr, jfxr
INTEGER :: iso, miso, liso, jso
INTEGER :: jt
REAL(8) :: sigbuf, wtN

REAL(8) :: ND, micA, mics, micstr
LOGICAL :: lskipit2
!REAL(8), SHARED :: pnumLoc(XS_BLOCK_GDIM,XS_BLOCK_RDIM)

!xtid = threadIdx%x;
jg = threadIdx%x+XS_BLOCK_GDIM*(blockIdx%x-1)
ig = jg+gb-1
IF (ig .GT. ge) RETURN

!ytid = threadIdx%y
ifxr  = threadIdx%y+XS_BLOCK_RDIM*(blockIdx%y-1)
if (ifxr .GT. nFxrTri) RETURN

jfxr = mapfxr(ifxr)-(iz-1)*nFxr
IF (jfxr .LT. 1) RETURN
IF (jfxr .GT. nFxr) RETURN

micA = 0.; mics = 0.; micstr = 0.

miso = AccNt(ifxr)+1
liso = AccNt(ifxr+1);

DO iso = miso, liso
  ND = pnum(iso)
  jso = iso-miso+1
  jt = itTmp(ifxr,jso); wtN = wtTmp(ifxr,jso)
  IF (jt.EQ.0) CYCLE
  lskipit2 = (wtN .GT. 0.99999)
  ! Left Point
  wtN = ND*wtN;

  sigbuf = a(ig, jt)
  micA = micA+wtN*sigbuf
  sigbuf = s(ig, jt)
  mics = mics+wtN*sigbuf
  sigbuf = str(ig, jt)
  micstr = micstr+wtN*sigbuf

  IF (lskipit2) CYCLE
  ! Right Point
  wtN = ND - wtN;
  jt = jt+1

  sigbuf = a(ig, jt)
  micA = micA+wtN*sigbuf
  sigbuf = s(ig, jt)
  mics = mics+wtN*sigbuf
  sigbuf = str(ig, jt)
  micstr = micstr+wtN*sigbuf
END DO

!ifxr = mapfxr(ifxr)
Mca(jg,jfxr) = micA
MCS(jg,jfxr) = mics
MCStr(jg,jfxr) = micstr
END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuGetMacSm_D(niso,pnum,mapfxr,gidSmp,ptrSmp,smp,McSm,itTmp,wtTmp,gb,ge,iz,ng,nFxrD)
USE CUDAFOR
IMPLICIT NONE
INTEGER, DEVICE :: niso(:),mapfxr(:),gidSmp(:,:),ptrSmp(:,:),itTmp(:,:)
REAL(GPU_PNUM_PRECISION), DEVICE :: pnum(:,:)
REAL(GPU_XS_PRECISION), DEVICE :: smp(:),wtTmp(:,:)
REAL(GPU_SM_PRECISION), DEVICE :: McSm(:,:)
INTEGER, VALUE :: gb, ge, iz, ng, nFxrD

!INTEGER :: xtid, ytid, ztid
INTEGER :: ig, og, ug, ifxr, jfxr, igs
INTEGER :: iso, miso
INTEGER :: itg, jtg, jt

REAL(8) :: sigbuf, wtN
REAL(8) :: ND, micsm
LOGICAL :: lskipit2
!xtid = threadIdx%x;
ig = threadIdx%x+XS_BLOCK_GDIM*(blockIdx%x-1)
IF (ig .GT. ng) RETURN

!ytid = threadIdx%y;
ug = blockIdx%y
og = ug+gb-1;
IF (og .GT. ge) RETURN

!ztid = threadIdx%z
ifxr  = threadIdx%z+XS_BLOCK_RDIM*(blockIdx%z-1)
if (ifxr .GT. nFxrD) RETURN

jfxr = mapfxr(ifxr)-(iz-1)*nFxr
IF (jfxr .LT. 1) RETURN
IF (jfxr .GT. nFxr) RETURN

miso = niso(ifxr);
igs = BegGrpScat(gb);
igs = BegGrpScat(og)-igs+1; itg = InScatRange(1,og); jtg = InScatRange(2,og)

IF (itg.GT.ig) RETURN
IF (jtg.LT.ig) RETURN
igs = igs+ig-itg
micsm = 0.;
DO iso = 1, miso
  ND = pnum(iso,ifxr)
  jt = itTmp(ifxr,iso); wtN = wtTmp(ifxr,iso)
  !print*, ifxr, iso, jt, wtN
  lskipit2 = (wtN .GT. 0.99999)

  wtN = ND*wtN;

  itg = gidSmp(og*2-1,jt); jtg = gidSmp(og*2,jt)
  IF (itg.GT.ig .OR. jtg.LT.ig) THEN
    sigbuf = 0.;
  ELSE
    jtg = ptrSmp(og,jt) ! First Np0 at group og
    itg = jtg+ig-itg ! Np0 at (og,ig)
    sigbuf = smp(itg)
  END IF
  micsm = micsm+sigbuf*wtN

  IF (lskipit2) CYCLE
  wtN = ND-wtN;
  jt = jt+1

  itg = gidSmp(og*2-1,jt); jtg = gidSmp(og*2,jt)
  IF (itg.GT.ig .OR. jtg.LT.ig) THEN
    sigbuf = 0.;
  ELSE
    jtg = ptrSmp(og,jt) ! First Np0 at group og
    itg = jtg+ig-itg ! Np0 at (og,ig)
    sigbuf = smp(itg)
  END IF
  micsm = micsm+sigbuf*wtN
END DO
!ifxr = mapfxr(ifxr)
!print*, igs, ifxr, SIZE(McSm,1), SIZE(McSm,2), (igs>SIZE(McSm,1)), (ifxr>SIZE(McSm,2))
McSm(igs,jfxr) = micsm;
END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuGetMacSm_T(AccNt,pnum,mapfxr,gidSmp,ptrSmp,smp,McSm,itTmp,wtTmp,gb,ge,iz,ng,nFxrTri)
USE CUDAFOR
IMPLICIT NONE
INTEGER, DEVICE :: AccNt(:),mapfxr(:),gidSmp(:,:),ptrSmp(:,:),itTmp(:,:)
REAL(GPU_PNUM_PRECISION), DEVICE :: pnum(:)
REAL(GPU_XS_PRECISION), DEVICE :: smp(:),wtTmp(:,:)
REAL(GPU_SM_PRECISION), DEVICE :: McSm(:,:)
INTEGER, VALUE :: gb, ge, iz, ng, nFxrTri

!INTEGER :: xtid, ytid, ztid
INTEGER :: ig, og, ug, ifxr, jfxr, igs
INTEGER :: iso, miso, liso
INTEGER :: itg, jtg, jt

REAL(8) :: sigbuf, wtN
REAL(8) :: ND, micsm
LOGICAL :: lskipit2

!xtid = threadIdx%x;
ig = threadIdx%x+XS_BLOCK_GDIM*(blockIdx%x-1)
IF (ig .GT. ng) RETURN

!ytid = threadIdx%y;
ug = blockIdx%y
og = ug+gb-1
IF (og .GT. ge) RETURN

!ztid = threadIdx%z
ifxr  = threadIdx%z+XS_BLOCK_RDIM*(blockIdx%z-1)
if (ifxr .GT. nFxrTri) RETURN

jfxr = mapfxr(ifxr)-(iz-1)*nFxr
IF (jfxr.LT.1) RETURN
IF (jfxr.GT.nFxr) RETURN

miso = AccNt(ifxr)+1
liso = AccNt(ifxr+1);
igs = BegGrpScat(gb);
igs = BegGrpScat(og)-igs+1; itg = InScatRange(1,og); jtg = InScatRange(2,og)
IF (itg.GT.ig) RETURN
IF (jtg.LT.ig) RETURN
igs = igs+ig-itg
micsm = 0.;
DO iso = miso, liso
  ND = pnum(iso)
  jt = itTmp(ifxr,iso-miso+1); wtN = wtTmp(ifxr,iso-miso+1)
  IF (jt.EQ.0) CYCLE
  !print*, ifxr, (iso-miso+1), SIZE(wtTmp,1), SIZE(wtTmp,2), SIZE(itTmp,1), SIZE(itTmp,2)
  !print*, ifxr, iso-miso, jt, wtN
  lskipit2 = (wtN .GT. 0.99999)
  wtN = ND*wtN;

  itg = gidSmp(og*2-1,jt); jtg = gidSmp(og*2,jt)
  IF (itg.GT.ig .OR. jtg.LT.ig) THEN
    sigbuf = 0.;
  ELSE
    jtg = ptrSmp(og,jt) ! First Np0 at group og
    itg = jtg+ig-itg ! Np0 at (og,ig)
    !print*, itg, SIZE(smp), (itg>SIZE(smp))
    sigbuf = smp(itg)
  END IF
  micsm = micsm+sigbuf*wtN

  IF (lskipit2) CYCLE
  wtN = ND-wtN;
  jt = jt+1

  itg = gidSmp(og*2-1,jt); jtg = gidSmp(og*2,jt)
  IF (itg.GT.ig .OR. jtg.LT.ig) THEN
    sigbuf = 0.;
  ELSE
    jtg = ptrSmp(og,jt) ! First Np0 at group og
    itg = jtg+ig-itg ! Np0 at (og,ig)
    !print*, itg, SIZE(smp), (itg>SIZE(smp))
    sigbuf = smp(itg)
  END IF
  micsm = micsm+sigbuf*wtN
END DO
!ifxr = mapfxr(ifxr)
McSm(igs,jfxr) = micsm;
END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuTreatMac_D(mapfxr,mapfsr,mapG2R,frA,frnF,frF,fsrst,nfsrs,McA,McnF,MckF,McT,McTr,McS,McStr,xst,rgb,rge,gb,ge,iz,nFxrD,lreso,lxst,caseID)
USE CUDAFOR
IMPLICIT NONE
INTEGER, DEVICE :: mapfxr(:),mapfsr(:),mapG2R(:),fsrst(:),nfsrs(:)
REAL(GPU_XS_PRECISION), DIMENSION(:,:), DEVICE :: McA,McnF,MckF,McT,McTr,McS,McStr
REAL(GPU_RES_PRECISION), DIMENSION(:,:), DEVICE :: frA,frnF,frF
REAL(GPU_PRECISION), DIMENSION(:,:), DEVICE :: xst
INTEGER, VALUE :: rgb, rge, gb, ge, iz, nFxrD, caseID
LOGICAL, VALUE :: lreso, lxst

INTEGER :: jg, ig, igr, ifxr, idxfxr, resfxr, ifsr, ibfsr, iefsr, offnfsr
REAL(8) :: macA, macnbuf, mackbuf, chibuf
jg = threadIdx%x+XS_BLOCK_GDIM*(blockIdx%x-1)
ig = jg+gb-1
if (ig .GT. ge) RETURN

ifxr  = threadIdx%y+XS_BLOCK_RDIM*(blockIdx%y-1)
if (ifxr .GT. nFxrD) RETURN

idxfxr = mapfxr(ifxr)-(iz-1)*nFxr
IF (idxfxr.LT.1) RETURN
IF (idxfxr.GT.nFxr) RETURN

resfxr = mapG2R(ifxr)

igr = ig-rgb+1
IF (caseID.NE.MacNF) macA = McA(jg,idxfxr)
IF (lreso) THEN
  IF (ig .LE. rge .AND. ig .GE. rgb .AND. resfxr .NE. 0) THEN
    IF (caseID.NE.MacNF) THEN
      macA = macA*frA(igr,resfxr);
      mackbuf = Mckf(jg,idxfxr)*frF(igr,resfxr);
      McA(jg,idxfxr) = macA
      Mckf(jg,idxfxr) = mackbuf;
    END IF
    IF (caseID.NE.MacBase) THEN
      macnbuf = Mcnf(jg,idxfxr)*frnF(igr,resfxr);
      Mcnf(jg,idxfxr) = macnbuf;
    END IF
  END IF
END IF
IF (caseID.NE.MacNF) THEN
  McT(jg,idxfxr) = macA+McS(jg,idxfxr)
  macA = macA+McStr(jg,idxfxr)
  McTr(jg,idxfxr) = macA
END IF
IF (lxst) THEN
  offnfsr = (iz-1)*nfsr
  ibfsr = fsrst(ifxr)
  iefsr = nfsrs(ifxr)
  iefsr = ibfsr+iefsr-1
  DO ifsr = ibfsr, iefsr
    !print*, mapfsr(ifsr), size(xst,2)
    xst(ig, mapfsr(ifsr)-offnfsr) = macA
  END DO
END IF
!print*, 'nufis', ig, idxfxr, McnF(ig,idxfxr)
END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuTreatMac_T(mapfxr,mapfsr,mapG2R,frA,fsrst,nfsrs,McA,McT,McTr,McS,McStr,xst,rgb,rge,gb,ge,iz,nFxrTri,lreso,lxst)
USE CUDAFOR
IMPLICIT NONE
INTEGER, DEVICE :: mapfxr(:),mapfsr(:),mapG2R(:),fsrst(:),nfsrs(:)
REAL(GPU_XS_PRECISION), DIMENSION(:,:), DEVICE :: McA,McT,McTr,McS,McStr
REAL(GPU_RES_PRECISION), DIMENSION(:,:), DEVICE :: frA
REAL(GPU_PRECISION), DIMENSION(:,:), DEVICE :: xst
INTEGER, VALUE :: rgb, rge, gb, ge, iz, nFxrTri
LOGICAL, VALUE :: lreso, lxst

INTEGER :: jg, ig, igr, ifxr, idxfxr, resfxr, ifsr, ibfsr, iefsr, offnfsr
REAL(8) :: macA

jg = threadIdx%x+XS_BLOCK_GDIM*(blockIdx%x-1)
ig = gb+jg-1
if (ig .GT. ge) RETURN

ifxr  = threadIdx%y+XS_BLOCK_RDIM*(blockIdx%y-1)
if (ifxr .GT. nFxrTri) RETURN

idxfxr = mapfxr(ifxr)-(iz-1)*nFxr
IF (idxfxr.LT.1) RETURN
IF (idxfxr.GT.nFxr) RETURN

resfxr = mapG2R(ifxr)

igr = ig-rgb+1

macA = McA(jg,idxfxr)
IF (lreso) THEN
  IF (ig .LE. rge .AND. ig .GE. rgb .AND. resfxr .NE. 0) THEN
    macA = macA*frA(igr,resfxr);
    McA(jg,idxfxr) = macA
  END IF
END IF

McT(jg,idxfxr) = macA+McS(jg,idxfxr)
macA = macA+McStr(jg,idxfxr)
McTr(jg,idxfxr) = macA

IF (lxst) THEN
  offnfsr = (iz-1)*nfsr
  ibfsr = fsrst(ifxr); iefsr = nfsrs(ifxr)
  iefsr = ibfsr+iefsr-1
  DO ifsr = ibfsr, iefsr
  !  print*, mapfsr(ifsr), size(xst,2)
    xst(ig,mapfsr(ifsr)-offnfsr) = macA
    !print*, 'Tr', ig, ifsr, macA
  END DO
END IF
END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuSetMOCPsi(mapfxrD,mapfxrT,mapfsrD,mapfsrT,fsrstD,fsrstT,nfsrD,nfsrT,phis,McnF,psi,iz,ng,nFxrD,nFxrTri)
USE CUDAFOR
IMPLICIT NONE
INTEGER, DIMENSION(:), DEVICE :: mapfxrD,mapfxrT,mapfsrD,mapfsrT,fsrstD,fsrstT,nfsrD,nfsrT
REAL(GPU_XS_PRECISION), DEVICE :: McnF(:,:)
REAL(GPU_FLUX_PRECISION), DEVICE :: phis(:,:)
REAL(GPU_SOURCE_PRECISION), DEVICE :: psi(:)
INTEGER, VALUE :: iz, ng, nFxrD, nFxrTri

INTEGER :: ifxr, ig, idxfxr, offnfsr
INTEGER :: ifsr, idxfsr, ibfsr, iefsr
REAL(8) :: psibuf

ifxr = threadIdx%x+blockDim%x*(blockIdx%x-1)
offnfsr = (iz-1)*nfsr

IF (ifxr .LE. nFXRD) THEN
  idxfxr = mapfxrD(ifxr)-(iz-1)*nFxr
  IF (idxfxr.LT.1) RETURN
  IF (idxfxr.GT.nFxr) RETURN
  ibfsr = fsrstD(ifxr)
  iefsr = nfsrD(ifxr)
  iefsr = iefsr+ibfsr-1
  DO ifsr = ibfsr, iefsr
    psibuf = 0.
    idxfsr = mapfsrD(ifsr)-offnfsr
    DO ig = 1, ng
      psibuf = psibuf+phis(ig,idxfsr)*McnF(ig,idxfxr)
    END DO
    psi(idxfsr) = psibuf
  END DO
ELSE IF (ifxr .LE. (nFxrD+nFxrTri)) THEN
  ifxr = ifxr - nFxrD;
  idxfxr = mapfxrT(ifxr)-(iz-1)*nFxr
  IF (idxfxr.LT.1) RETURN
  IF (idxfxr.GT.nFxr) RETURN
  ibfsr = fsrstT(ifxr)
  iefsr = nfsrT(ifxr)
  iefsr = iefsr+ibfsr-1
  DO ifsr = ibfsr, iefsr
    psibuf = 0.
    idxfsr = mapfsrT(ifsr)-offnfsr
    psi(idxfsr) = psibuf
  END DO
ELSE
  RETURN
END IF
END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuSetSource(mapfxrD,mapfxrT,mapfsrD,mapfsrT,fsrstD,fsrstT,nfsrD,nfsrT,lvoidD,lvoidT,pinidD,pinidT,&
                                McTr,McSm,Chi,phis,psi,src,xst,AxSrc,AxPxs,reigv,gb,ge,iz,nFxrD,nFxrTri,lLkgSplit)
USE CUDAFOR
IMPLICIT NONE
INTEGER, DIMENSION(:), DEVICE :: mapfxrD,mapfxrT,mapfsrD,mapfsrT,fsrstD,fsrstT,nfsrD,nfsrT,pinidD,pinidT
LOGICAL, DEVICE :: lvoidD(:), lvoidT(:)
REAL(GPU_XS_PRECISION), DIMENSION(:,:), DEVICE :: McTr, Chi
REAL(GPU_SM_PRECISION), DIMENSION(:,:), DEVICE :: McSm
REAL(8), DIMENSION(:,:), DEVICE :: AxSrc,AxPxs
REAL(GPU_FLUX_PRECISION), DEVICE :: phis(:,:)
REAL(GPU_PRECISION), DEVICE :: xst(:,:)
REAL(GPU_SOURCE_PRECISION), DEVICE :: psi(:), src(:,:)
REAL(8), VALUE :: reigv
INTEGER, VALUE :: gb, ge, iz, nFxrD, nFxrTri
LOGICAL, VALUE :: lLkgSplit

LOGICAL :: lvoidFxr
INTEGER :: ifxr, idxfxr, ig, jg, ipin, offnFsr
INTEGER :: og, igs, ibgs, iegs, ifsr, idxfsr, ibfsr, iefsr
REAL(8) :: srcbuf
REAL(8) :: axsrcbuf, axpxsbuf, xstbuf, chibuf

jg = threadIdx%x+XS_BLOCK_GDIM*(blockIdx%x-1)
ig = jg+gb-1
IF (ig .GT. ge) RETURN
ifxr = threadIdx%y+XS_BLOCK_RDIM*(blockIdx%y-1)

offnFsr = (iz-1)*nfsr

igs = BegGrpScat(gb)
ibgs = BegGrpScat(ig)-igs+1; iegs = BegGrpScat(ig+1)-igs
chibuf = 0.;

IF (ifxr .LE. nFXRD) THEN
  idxfxr = mapfxrD(ifxr)-(iz-1)*nFxr
  IF (idxfxr.LT.1) RETURN
  IF (idxfxr.GT.nFxr) RETURN
  lvoidFxr = lvoidD(ifxr); ipin = pinidD(ifxr);
  axsrcbuf = AxSrc(ipin,jg); axpxsbuf = AxPxs(ipin,jg)
  ibfsr = fsrstD(ifxr)
  iefsr = nfsrD(ifxr)
  iefsr = iefsr+ibfsr-1
  IF(ig.LE.nchi) chibuf = Chi(ig,ifxr)
  xstbuf = McTr(jg, idxfxr)+axpxsbuf;
  DO ifsr = ibfsr, iefsr
    idxfsr = mapfsrD(ifsr)-offnFsr
    srcbuf = reigv*chibuf*psi(idxfsr)
    og = InScatRange(1,ig)
    DO igs = ibgs, iegs
      srcbuf = srcbuf + McSm(igs,idxfxr)*phis(og,idxfsr)
      og = og+1
    END DO
    IF ((.NOT. lLkgSplit) .OR. axsrcbuf .LT. 0. .AND. .NOT. lvoidFxr) THEN
      srcbuf = srcbuf - axsrcbuf
    END IF
    srcbuf = srcbuf/xstbuf
    src(ig,idxfsr) = srcbuf
    xst(ig,idxfsr) = xstbuf
  END DO
ELSE IF (ifxr .LE. (nFxrD+nFxrTri)) THEN
  ifxr = ifxr - nFxrD
  idxfxr = mapfxrT(ifxr)-(iz-1)*nFxr
  IF (idxfxr.LT.1) RETURN
  IF (idxfxr.GT.nFxr) RETURN
  lvoidFxr = lvoidT(ifxr); ipin = pinidT(ifxr);
  axsrcbuf = AxSrc(ipin,jg); axpxsbuf = AxPxs(ipin,jg)
  ibfsr = fsrstT(ifxr)
  iefsr = nfsrT(ifxr)
  iefsr = iefsr+ibfsr-1
  xstbuf = McTr(jg, idxfxr)+axpxsbuf;
  DO ifsr = ibfsr, iefsr
    idxfsr = mapfsrT(ifsr)-offnFsr
    srcbuf = 0
    og = InScatRange(1,ig)
    DO igs = ibgs, iegs
      srcbuf = srcbuf + McSm(igs,idxfxr)*phis(og,idxfsr)
      og = og+1
    END DO
    IF ((.NOT.lLkgSplit) .OR. axsrcbuf .LT. 0. .AND. .NOT. lvoidFxr) THEN
      srcbuf = srcbuf - axsrcbuf
    END IF
    srcbuf = srcbuf/xstbuf
    src(ig,idxfsr) = srcbuf
    xst(ig,idxfsr) = xstbuf
  END DO
ELSE
  RETURN
END IF

END SUBROUTINE

SUBROUTINE SetupCuMacXS(ScatOrder, lxst, iz, gb, ge, caseID)
USE Core_mod, ONLY : GroupInfo
USE PE_Mod, ONLY : PE
IMPLICIT NONE
INTEGER :: ScatOrder, iz, gb, ge, caseID
LOGICAL :: lxst
INTEGER :: rgb, rge, ierr, isync, i, j, jz

INTEGER :: istat
!INTEGER :: printval(20)

jz = iz-PE%myzb+1
rgb = GroupInfo%nofg + 1; rge = GroupInfo%norg+rgb-1
CALL cuGetMacXS_D <<< blockXSD, trdXS, 0, cuDevice%myStream >>>(cuBlockFxr%nisoD,cuBlockFxr%pnumD,cuBlockFxr%mapglobalidD, &
      cuIsoData%siga,cuIsoData%sigkf,cuIsoData%signf,cuIsoData%sigs,cuIsoData%sigstr,&
      cuCoreMacXs%XSa,cuCoreMacXs%XSkf,cuCoreMacXs%XSnf,cuCoreMacXs%XSS,cuCoreMacXs%XSStr,&
      cuBlockFxr%itTempD,cuBlockFxr%wtTempD,gb,ge,jz,FxrD%nFxrD, caseID)
isync = cudaDeviceSynchronize()
ierr = cudaGetLastError()
if (ierr.NE.0) print*, __FILE__, __LINE__, cudaGetErrorString(ierr)
IF (caseID.NE.MacNF) THEN
  CALL cuGetMacXS_T <<< blockXST, trdXS, 0, cuDevice%myStream >>>(cuBlockFxr%AccNT,cuBlockFxr%pnumT,cuBlockFxr%mapglobalidT,&
        cuIsoData%siga,cuIsoData%sigs,cuIsoData%sigstr,cuCoreMacXs%XSa,cuCoreMacXs%XSS,cuCoreMacXs%XSStr,cuBlockFxr%itTempT,cuBlockFxr%wtTempT,&
        gb,ge,jz,FxrTri%nFxrTri)
  isync = cudaDeviceSynchronize()
  ierr = cudaGetLastError()
  if (ierr.NE.0) print*, __FILE__, __LINE__, cudaGetErrorString(ierr)
  CALL cuGetMacSm_D <<< blockSmD, trdSm, 0, cuDevice%myStream >>>(cuBlockFxr%nisoD,cuBlockFxr%pnumD,cuBlockFxr%mapglobalidD, &
        cuIsoData%gidSmp0,cuIsoData%ptrSmp0,cuIsoData%smp0,cuCoreMacXs%XSsm,cuBlockFxr%itTempD,cuBlockFxr%wtTempD,&
        gb,ge,jz,IsoData%ng,FxrD%nFxrD)
  isync = cudaDeviceSynchronize()
  ierr = cudaGetLastError()
  if (ierr.NE.0) print*, __FILE__, __LINE__, cudaGetErrorString(ierr)
  CALL cuGetMacSm_T <<< blockSmT, trdSm, 0, cuDevice%myStream >>>(cuBlockFxr%AccNT,cuBlockFxr%pnumT,cuBlockFxr%mapglobalidT,&
        cuIsoData%gidSmp0,cuIsoData%ptrSmp0,cuIsoData%smp0,cuCoreMacXs%XSsm,cuBlockFxr%itTempT,cuBlockFxr%wtTempT,&
        gb,ge,jz,IsoData%ng,FxrTri%nFxrTri)
  isync = cudaDeviceSynchronize()
  ierr = cudaGetLastError()
  if (ierr.NE.0) print*, __FILE__, __LINE__, cudaGetErrorString(ierr)
END IF
CALL cuTreatMac_D <<< blockXSD, trdXS, 0, cuDevice%myStream >>>(cuBlockFxr%mapglobalidD,cuBlockFxr%mapfsrD,cuBlockFxr%mapG2RD,cuBlockFxr%FresoAD,cuBlockFxr%FresoNF,cuBlockFxr%FresoF,&
      cuBlockFxr%fsrstD,cuBlockFxr%nfsrD,cuCoreMacXs%XSa,cuCoreMacXs%XSnf,cuCoreMacXs%XSkf,cuCoreMacXs%XSt,cuCoreMacXs%XStr,cuCoreMacXs%XSS, cuCoreMacXs%XSStr,&
      cuMOC%xst,rgb,rge,gb,ge,jz,FxrD%nFxrD,lrestrmt,lxst,caseID)
isync = cudaDeviceSynchronize()
ierr = cudaGetLastError()
if (ierr.NE.0) print*, __FILE__, __LINE__, cudaGetErrorString(ierr)
IF (caseID.NE.MacNF) THEN
  CALL cuTreatMac_T <<< blockXST, trdXS, 0, cuDevice%myStream >>>(cuBlockFxr%mapglobalidT,cuBlockFxr%mapfsrT,cuBlockFxr%mapG2RT,cuBlockFxr%FresoAT,cuBlockFxr%fsrstT,cuBlockFxr%nfsrT,&
        cuCoreMacXs%XSa,cuCoreMacXs%XSt,cuCoreMacXs%XStr,cuCoreMacXs%XSS,cuCoreMacXs%XSStr,cuMOC%xst,rgb, rge, gb, ge, jz, FxrTri%nFxrTri,lrestrmt,lxst)
  isync = cudaDeviceSynchronize()
  ierr = cudaGetLastError()
  if (ierr.NE.0) print*, __FILE__, __LINE__, cudaGetErrorString(ierr)
END IF

IF (ScatOrder<1) RETURN

END SUBROUTINE

SUBROUTINE SetupCuMacnF(iz)
USE Core_mod, ONLY : GroupInfo
USE PE_Mod, ONLY : PE
IMPLICIT NONE
INTEGER :: iz
LOGICAL :: lxst = .FALSE.
INTEGER :: rgb, rge, ierr, isync, i, j, jz

INTEGER :: istat
!INTEGER :: printval(20)

jz = iz-PE%myzb+1
rgb = GroupInfo%nofg + 1; rge = GroupInfo%norg+rgb-1
CALL cuGetMacXS_D <<< blockXSD, trdXS, 0, cuDevice%myStream >>>(cuBlockFxr%nisoD,cuBlockFxr%pnumD,cuBlockFxr%mapglobalidD, &
      cuIsoData%siga,cuIsoData%sigkf,cuIsoData%signf,cuIsoData%sigs,cuIsoData%sigstr,&
      cuCoreMacXs%XSa,cuCoreMacXs%XSkf,cuCoreMacXs%XSnf,cuCoreMacXs%XSS,cuCoreMacXs%XSStr,&
      cuBlockFxr%itTempD,cuBlockFxr%wtTempD,1,IsoData%ng,jz,FxrD%nFxrD,MacNF)
!isync = cudaDeviceSynchronize()
!ierr = cudaGetLastError()
!if (ierr.NE.0) print*, __FILE__, __LINE__, cudaGetErrorString(ierr)
CALL cuTreatMac_D <<< blockXSD, trdXS, 0, cuDevice%myStream >>>(cuBlockFxr%mapglobalidD,cuBlockFxr%mapfsrD,cuBlockFxr%mapG2RD,cuBlockFxr%FresoAD,cuBlockFxr%FresoNF,cuBlockFxr%FresoF,&
      cuBlockFxr%fsrstD,cuBlockFxr%nfsrD,cuCoreMacXs%XSa,cuCoreMacXs%XSnf,cuCoreMacXs%XSkf,cuCoreMacXs%XSt,cuCoreMacXs%XStr,cuCoreMacXs%XSS, cuCoreMacXs%XSStr,&
      cuMOC%xst,rgb,rge,1,IsoData%ng,jz,FxrD%nFxrD,lrestrmt,lxst,MacNF)
!isync = cudaDeviceSynchronize()
!ierr = cudaGetLastError()
!if (ierr.NE.0) print*, __FILE__, __LINE__, cudaGetErrorString(ierr)

END SUBROUTINE

SUBROUTINE SetupCuSrc(CoreInfo, AxSrc, AxPxs, eigv, ScatOrder, iz, gb, ge, lLkgSplit)
USE TYPEDEF,      ONLY : coreinfo_type
USE PE_Mod,       ONLY : PE
IMPLICIT NONE
TYPE(CoreInfo_Type) :: CoreInfo
REAL, POINTER :: AxSrc(:,:,:), AxPxs(:,:,:)
REAL :: eigv
INTEGER :: ScatOrder, iz, gb, ge
LOGICAL :: lLkgSplit

REAL(8), DEVICE, ALLOCATABLE, SAVE :: AxSrcDev(:,:), AxPxsDev(:,:)
REAL :: reigv
INTEGER :: nxyz, ng, ierr, i, j, jz
!REAL(4), ALLOCATABLE :: srcmoc(:,:)

nxyz = CoreInfo%nxy*(PE%myze-PE%myzb+1); ng = IsoData%ng

ALLOCATE(AxSrcDev(nxyz,gb:ge), AxPxsDev(nxyz,gb:ge))
ierr = cudaMemcpy(AxSrcDev, AxSrc(1:,PE%myzb:,gb:), nxyz*(ge-gb+1), cudaMemcpyHostToDevice)
ierr = cudaMemcpy(AxPxsDev, AxPxs(1:,PE%myzb:,gb:), nxyz*(ge-gb+1), cudaMemcpyHostToDevice)

jz = iz-PE%myzb+1
reigv = 1./eigv

CALL cuSetSource <<< blockXS, trdXS, 0, cuDevice%myStream >>>(cuBlockFxr%mapglobalidD,cuBlockFxr%mapglobalidT,cuBlockFxr%mapfsrD,cuBlockFxr%mapfsrT,&
      cuBlockFxr%fsrstD,cuBlockFxr%fsrstT,cuBlockFxr%nfsrD,cuBlockFxr%nfsrT,cuBlockFxr%lvoidD,cuBlockFxr%lvoidT,cuBlockFxr%pinidD,cuBlockFxr%pinidT,&
      cuCoreMacXs%XStr,cuCoreMacXs%XSsm,cuBlockFxr%chi,cuMOC%phis,cuMOC%psi,cuMOC%src,cuMOC%xst,AxSrcDev,AxPxsDev,reigv,gb,ge,jz,FxrD%nFxrD,FxrTri%nFxrTri,&
      lLkgSplit)

DEALLOCATE(AxSrcDev, AxPxsDev)

IF (ScatOrder < 1) RETURN

END SUBROUTINE

SUBROUTINE SetupCuPsiNChi(iz)
! Sorry, but it contains only Psi setup
!USE geom, ONLY : core
USE PE_Mod, ONLY : PE
IMPLICIT NONE
INTEGER :: iz
INTEGER :: nfxr, nfxrd, nfxrtri, i, jz
TYPE(dim3) :: nblk, ntrd
!REAL(4), ALLOCATABLE :: psimoc(:)

nfxrd = FxrD%nFxrD; nfxrtri = FxrTri%nFxrTri; nfxr = nfxrd+nfxrtri
nblk = dim3(nfxr/256+1,1,1); ntrd = dim3(256,1,1)

jz = iz-PE%myzb+1

CALL cuSetMOCPsi <<< nblk, ntrd, 0, cuDevice%myStream>>> (cuBlockFxr%mapglobalidD,cuBlockFxr%mapglobalidT,cuBlockFxr%mapfsrD,cuBlockFxr%mapfsrT,&
      cuBlockFxr%fsrstD,cuBlockFxr%fsrstT,cuBlockFxr%nfsrD,cuBlockFxr%nfsrT,cuMOC%phis,cuCoreMacXs%XSnf,cuMOC%psi,jz,IsoData%ng,nFxrD,nFxrTri)

!ALLOCATE(psimoc(Core%ncorefsr)); psimoc = cuMOC%psi
!DO i = 1, Core%ncorefsr
!  print*, 'psi', i, psimoc(i)
!END DO

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuSetPlnLSigP_MLG(nisoD,pnumD,mapfxrD,mapfsrD,nfsrD,fsrstD,&
                                                AccNt,nisoT,pnumT,mapfxrT,mapfsrT,nfsrT,fsrstT,&
                                                mapnuclD,mapnuclT,&
                                                lamsigp,sigp,Siglp,xst,&
                                                rgb,rge,nlv,iz,nFxrD,nFxrTri)
USE CUDAFOR
IMPLICIT NONE
INTEGER,DEVICE,DIMENSION(:) :: nisoD,AccNt,nisoT,mapfxrD,mapfxrT,&
                              mapfsrD,mapfsrT,nfsrD,nfsrT,fsrstD,fsrstT
INTEGER,DEVICE :: mapnuclD(:,:),mapnuclT(:)
REAL(GPU_PNUM_PRECISION),DEVICE :: pnumD(:,:),pnumT(:)
REAL(GPU_XS_PRECISION),DEVICE :: lamsigp(:,:),sigp(:),Siglp(:,:)
REAL(GPU_PRECISION),DEVICE :: xst(:,:)
INTEGER,VALUE :: rgb,rge,nlv,iz,nFxrD,nFxrTri

INTEGER :: ig, idlv, ifxr, idxfxr, ifxrt, jfxr, ifsr, ibfsr, iefsr, idxfsr, offnfsr
INTEGER :: iso, miso, liso, idiso
REAL(8) :: sigtr, sigp1g, ND

idlv = threadIdx%x+FSP_BLOCK_GDIM*(blockIdx%x-1)
ig = (idlv-1)/nlv
ig = ig+rgb
IF (ig.GT.rge) RETURN
ifxr = threadIdx%y+FSP_BLOCK_RDIM*(blockIdx%y-1)

offnfsr = (iz-1)*nfsr
sigtr = 0.;
IF (ifxr.LE.nFxrD) THEN
  idxfxr = mapfxrD(ifxr)-(iz-1)*nFxr
  IF (idxfxr.LT.1) RETURN
  IF (idxfxr.GT.nFxr) RETURN
  ibfsr = fsrstD(ifxr)
  iefsr = nfsrD(ifxr)
  iefsr = ibfsr+iefsr-1
  miso = nisoD(ifxr)
  DO iso = 1, miso
    ND = pnumD(iso,ifxr)
    idiso = mapnuclD(iso,ifxr)
    sigp1g = sigp(idiso)
    IF (sigp1g.EQ.0.0) THEN
      sigp1g = lamsigp(ig,idiso)
    END IF
    sigtr = sigtr + sigp1g*ND
  END DO
  Siglp(idlv,idxfxr) = sigtr
  DO ifsr = ibfsr, iefsr
    idxfsr = mapfsrD(ifsr)-offnfsr
    xst(idlv,idxfsr) = sigtr
  END DO
ELSE IF (ifxr.LE.(nFxrD+nFxrTri)) THEN
  ifxrt = ifxr-nFxrD
  idxfxr = mapfxrT(ifxrt)-(iz-1)*nFxr
  IF (idxfxr.LT.1) RETURN
  IF (idxfxr.GT.nFxr) RETURN
  ibfsr = fsrstT(ifxrt)
  iefsr = nfsrT(ifxrt)
  iefsr = ibfsr+iefsr-1
  miso = AccNt(ifxrt)+1; liso = miso-1+nisoT(ifxrt)
  DO iso = miso,liso
    ND = pnumT(iso)
    idiso = mapnuclT(iso)
    sigp1g = sigp(idiso)
    IF (sigp1g.EQ.0.0) THEN
      sigp1g = lamsigp(ig,idiso)
    END IF
    sigtr = sigtr + sigp1g*ND
!    IF (idxfxr.EQ.2) print*, 'Trivial', iso-miso+1, sigp1g, ND
  END DO
  Siglp(idlv,idxfxr) = sigtr
  DO ifsr = ibfsr, iefsr
    idxfsr = mapfsrT(ifsr)-offnfsr
    xst(idlv,idxfsr) = sigtr
  END DO
END IF
END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuSetPlnLSigP1G_MLG(nisoD,pnumD,mapfxrD,mapfsrD,nfsrD,fsrstD,&
                                                AccNt,nisoT,pnumT,mapfxrT,mapfsrT,nfsrT,fsrstT,&
                                                mapnuclD,mapnuclT,&
                                                sigp,Siglp,xst,&
                                                nlv,iz,nFxrD,nFxrTri)
USE CUDAFOR
IMPLICIT NONE
INTEGER,DEVICE,DIMENSION(:) :: nisoD,AccNt,nisoT,mapfxrD,mapfxrT,&
                              mapfsrD,mapfsrT,nfsrD,nfsrT,fsrstD,fsrstT
INTEGER,DEVICE :: mapnuclD(:,:),mapnuclT(:)
REAL(GPU_PNUM_PRECISION),DEVICE :: pnumD(:,:),pnumT(:)
REAL(GPU_XS_PRECISION),DEVICE :: sigp(:),Siglp(:,:)
REAL(GPU_PRECISION),DEVICE :: xst(:,:)
INTEGER,VALUE :: nlv,iz,nFxrD,nFxrTri

INTEGER :: idlv, ifxr, idxfxr, ifxrt, jfxr, ifsr, idxfsr, ibfsr, iefsr, offnfsr
INTEGER :: iso, miso, liso, idiso
REAL(8) :: sigtr, sigp1g, ND

idlv = threadIdx%x+FSP_BLOCK_GDIM*(blockIdx%x-1)
IF (idlv.GT.nlv) RETURN
ifxr = threadIdx%y+FSP_BLOCK_RDIM*(blockIdx%y-1)
offnfsr = (iz-1)*nfsr
sigtr = 0.;
IF (ifxr.LE.nFxrD) THEN
  idxfxr = mapfxrD(ifxr)-(iz-1)*nFxr
  IF (idxfxr.LT.1) RETURN
  IF (idxfxr.GT.nFxr) RETURN
  ibfsr = fsrstD(ifxr)
  iefsr = nfsrD(ifxr)
  iefsr = ibfsr+iefsr-1
  miso = nisoD(ifxr)
  DO iso = 1, miso
    ND = pnumD(iso,ifxr)
    idiso = mapnuclD(iso,ifxr)
    sigp1g = sigp(idiso)
    sigtr = sigtr + sigp1g*ND
  END DO
  Siglp(idlv,idxfxr) = sigtr
  DO ifsr = ibfsr, iefsr
    idxfsr = mapfsrD(ifsr)-offnfsr
    xst(idlv,idxfsr) = sigtr
  END DO
ELSE IF (ifxr.LE.(nFxrD+nFxrTri)) THEN
  ifxrt = ifxr-nFxrD
  idxfxr = mapfxrT(ifxrt)-offnfsr
  IF (idxfxr.LT.1) RETURN
  IF (idxfxr.GT.nFxr) RETURN
  ibfsr = fsrstT(ifxrt)
  iefsr = nfsrT(ifxrt)
  iefsr = ibfsr+iefsr-1
  miso = AccNt(ifxrt)+1; liso = miso+nisoT(ifxrt)-1
  DO iso = miso,liso
    ND = pnumT(iso)
    idiso = mapnuclT(iso)
    sigp1g = sigp(idiso)
    sigtr = sigtr + sigp1g*ND
  END DO
  Siglp(idlv,idxfxr) = sigtr
  DO ifsr = ibfsr, iefsr
    idxfsr = mapfsrT(ifsr)-offnfsr
    xst(idlv,idxfsr) = sigtr
  END DO
END IF
END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuSetPlnXST_MLG(mapfxrD,mapfsrD,nfsrD,fsrstD,&
                                                mapfxrT,mapfsrT,nfsrT,fsrstT,&
                                                mapGR2GD,mapGR2GT,&
                                                f_maclv,Siglp,xst,&
                                                FnAdjD,FtAdjD,FnAdjT,FtAdjT,&
                                                iz,rgb,rge,nlv,nGenD,nGenTri)
USE CUDAFOR
IMPLICIT NONE
INTEGER,DEVICE,DIMENSION(:) :: mapfxrD,mapfxrT,&
                              mapfsrD,mapfsrT,nfsrD,nfsrT,fsrstD,fsrstT,mapGR2GD,mapGR2GT
REAL(GPU_XS_PRECISION),DEVICE :: f_maclv(:,:),Siglp(:,:)
REAL(GPU_RES_PRECISION),DEVICE :: FnAdjD(:,:),FtAdjD(:,:,:),FnAdjT(:,:),FtAdjT(:,:,:)
REAL(GPU_PRECISION),DEVICE :: xst(:,:)
INTEGER,VALUE :: iz,rgb,rge,nlv,nGenD,nGenTri

INTEGER :: ig, ilv, idlv, ifxr, idxfxr, ifxrt, jfxr, ifsr, ibfsr, iefsr, idxfsr, offnfsr
REAL(8) :: sigtr, lv

idlv = threadIdx%x+FSP_BLOCK_GDIM*(blockIdx%x-1)
ig = (idlv-1)/nlv
ilv = idlv-ig*nlv
ig = ig+rgb
IF (ig.GT.rge) RETURN
ifxr = threadIdx%y+FSP_BLOCK_RDIM*(blockIdx%y-1)
offnfsr = (iz-1)*nfsr

lv = f_maclv(ilv,iz)
IF (ifxr.GT.(nGenD+nGenTri)) RETURN
IF (ifxr.LE.nGenD) THEN
  jfxr = mapGR2GD(ifxr)
  idxfxr = mapfxrD(jfxr)-(iz-1)*nFxr
  IF (idxfxr.LT.1) RETURN
  IF (idxfxr.GT.nFxr) RETURN
  ibfsr = fsrstD(jfxr); iefsr = nfsrD(jfxr);
  iefsr = iefsr+ibfsr-1

  sigtr = lv*FnAdjD(ig,ifxr)*FtAdjD(ilv,ig,ifxr)
!  print*, idxfxr, lv, FnAdjD(ig,ifxr), FtAdjD(ilv,ig,ifxr)
  sigtr = Siglp(idlv,idxfxr) + sigtr
  DO ifsr = ibfsr,iefsr
    idxfsr = mapfsrD(ifsr)-offnfsr
    xst(idlv,idxfsr) = sigtr
  END DO
ELSE
  ifxrt = ifxr-nGenD
  jfxr = mapGR2GT(ifxrt)
  idxfxr = mapfxrT(jfxr)-(iz-1)*nFxr
  IF (idxfxr.LT.1) RETURN
  IF (idxfxr.GT.nFxr) RETURN
  ibfsr = fsrstT(jfxr); iefsr = nfsrT(jfxr);
  iefsr = ibfsr+iefsr-1

  sigtr = lv*FnAdjT(ig,ifxrt)*FtAdjT(ilv,ig,ifxrt)
  sigtr = Siglp(idlv,idxfxr) + sigtr
  DO ifsr = ibfsr,iefsr
    idxfsr = mapfsrT(ifsr)-offnfsr
    xst(idlv,idxfsr) = sigtr
  END DO
END IF
END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuSetPlnXST1G_MLG(mapfxrD,mapfsrD,nfsrD,fsrstD,&
                                                mapfxrT,mapfsrT,nfsrT,fsrstT,&
                                                mapC2GD,mapC2GT,mapA2GT,&
                                                f_maclv1G,c_maclv1G,Siglp,xst,&
                                                iz,nlv,lCLD,nCldD,nCldT,nAICT)
USE CUDAFOR
IMPLICIT NONE
INTEGER,DEVICE,DIMENSION(:) :: mapfxrD,mapfxrT,&
                              mapfsrD,mapfsrT,nfsrD,nfsrT,fsrstD,fsrstT,mapC2GD,mapC2GT,mapA2GT
REAL(GPU_XS_PRECISION),DEVICE :: f_maclv1G(:,:),c_maclv1G(:),Siglp(:,:)
REAL(GPU_PRECISION),DEVICE :: xst(:,:)
INTEGER,VALUE :: iz,nlv,nCldD,nCldT,nAICT
LOGICAL,VALUE :: lCLD

INTEGER :: idlv, ifxr, idxfxr, ifxrt, jfxr, ifsr, idxfsr, ibfsr, iefsr, offnfsr
REAL(8) :: sigtr

idlv = threadIdx%x+FSP_BLOCK_GDIM*(blockIdx%x-1)
IF (idlv.GT.nlv) RETURN
ifxr = threadIdx%y+FSP_BLOCK_RDIM*(blockIdx%y-1)
offnfsr = (iz-1)*nfsr
IF (lCLD) THEN
  sigtr = c_maclv1G(idlv)
  IF (ifxr.LE.nCldD) THEN
    jfxr = mapC2GD(ifxr)
    idxfxr = mapfxrD(jfxr)-(iz-1)*nFxr
    IF (idxfxr.LT.1) RETURN
    IF (idxfxr.GT.nFxr) RETURN
    ibfsr = fsrstD(jfxr); iefsr = nfsrD(jfxr);
    iefsr = iefsr+ibfsr-1

    sigtr = Siglp(idlv,idxfxr)+sigtr
    DO ifsr = ibfsr,iefsr
      idxfsr = mapfsrD(ifsr)-offnfsr
      xst(idlv,idxfsr) = sigtr
    END DO
  ELSE IF (ifxr.LE.(nCldD+nCldT)) THEN
    ifxrt = ifxr-nCldD
    jfxr = mapC2GT(ifxrt)
    idxfxr = mapfxrT(jfxr)-(iz-1)*nFxr
    IF (idxfxr.LT.1) RETURN
    IF (idxfxr.GT.nFxr) RETURN
    ibfsr = fsrstT(jfxr); iefsr = nfsrT(jfxr);
    iefsr = ibfsr+iefsr-1

    sigtr = Siglp(idlv,idxfxr)+sigtr
    DO ifsr = ibfsr,iefsr
      idxfsr = mapfsrT(ifsr)
      xst(idlv,idxfsr) = sigtr
    END DO
  END IF
ELSE
  sigtr = f_maclv1G(idlv,iz)
  IF (ifxr.LE.nAICT) THEN
    jfxr = mapA2GT(ifxr); idxfxr = mapfxrT(jfxr)
    ibfsr = fsrstT(jfxr); iefsr = nfsrT(jfxr);
    iefsr = iefsr+ibfsr-1

    sigtr = Siglp(idlv,idxfxr)+sigtr
    DO ifsr = ibfsr,iefsr
    idxfsr = mapfsrT(ifsr)-offnfsr
    xst(idlv,idxfsr) = sigtr
    END DO
  END IF
END IF
END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuSetSubGrpSrc(mapfxrD,mapfsrD,nfsrD,fsrstD,&
                                            mapfxrT,mapfsrT,nfsrT,fsrstT,&
                                            Siglp,src,xst,&
                                            iz,nFxrD,nFxrTri,maxLv)
USE CUDAFOR
IMPLICIT NONE
INTEGER,DEVICE,DIMENSION(:) :: mapfxrD,mapfxrT,mapfsrD,mapfsrT,nfsrD,nfsrT,fsrstD,fsrstT
REAL(GPU_XS_PRECISION),DEVICE :: Siglp(:,:)
REAL(GPU_PRECISION),DEVICE :: xst(:,:)
REAL(GPU_SOURCE_PRECISION),DEVICE :: src(:,:)
INTEGER,VALUE :: iz,nFxrD,nFxrTri,maxLv

INTEGER :: idlv, ifxr, idxfxr, ifsr, idxfsr, ibfsr, iefsr, offnfsr
REAL(8) :: sigLamP
REAL(8) :: srcPnt

idlv = threadIdx%x+FSP_BLOCK_GDIM*(blockIdx%x-1)
IF (idlv.GT.maxLv) RETURN
ifxr = threadIdx%y+FSP_BLOCK_RDIM*(blockIdx%y-1)
offnfsr = (iz-1)*nfsr
IF (ifxr.LE.nFxrD) THEN
  idxfxr = mapfxrD(ifxr)-(iz-1)*nFxr
  IF (idxfxr.LT.1) RETURN
  IF (idxfxr.GT.nFxr) RETURN
  sigLamP = Siglp(idlv,idxfxr)

  ibfsr = fsrstD(ifxr); iefsr = nfsrD(ifxr)
  iefsr = ibfsr+iefsr-1
  srcPnt = sigLamP/xst(idlv,mapfsrD(ibfsr)-offnfsr)
  DO ifsr = ibfsr,iefsr
    idxfsr = mapfsrD(ifsr)-offnfsr
    src(idlv,idxfsr) = srcPnt
  END DO
ELSE IF (ifxr.LE.(nFxrD+nFxrTri)) THEN
  ifxr = ifxr-nFxrD
  idxfxr = mapfxrT(ifxr)-(iz-1)*nFxr
  IF (idxfxr.LT.1) RETURN
  IF (idxfxr.GT.nFxr) RETURN
  sigLamP = Siglp(idlv,idxfxr)

  ibfsr = fsrstT(ifxr); iefsr = nfsrT(ifxr)
  iefsr = ibfsr+iefsr-1
  srcPnt = sigLamP/xst(idlv,mapfsrT(ibfsr)-offnfsr)
  DO ifsr = ibfsr, iefsr
    idxfsr = mapfsrT(ifsr)-offnfsr
    src(idlv,idxfsr) = srcPnt
  END DO
END IF
END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuEquipXsGen_MLG(mapfxrD,mapfsrD,nfsrD,fsrstD,volD,mapGR2GD,&
                                            mapfxrT,mapfsrT,nfsrT,fsrstT,volT,mapGR2GT,&
                                            GenStD,GenStT,lPinRes,&
                                            Siglp,xst,phis,xseqD,xseqT,xseqPin,&
                                            nPin,rgb,rge,iz,nlv)
USE CUDAFOR
IMPLICIT NONE
INTEGER,DEVICE,DIMENSION(:) :: mapfxrD,mapfxrT,mapfsrD,mapfsrT,nfsrD,nfsrT,fsrstD,fsrstT
INTEGER,DEVICE,DIMENSION(:) :: mapGR2GD,mapGR2GT,GenStD,GenStT
LOGICAL,DEVICE,DIMENSION(:) :: lPinRes
REAL(GPU_XS_PRECISION),DEVICE :: Siglp(:,:)
REAL(GPU_XS_PRECISION),DEVICE :: volD(:),volT(:)
REAL(GPU_PRECISION),DEVICE :: xst(:,:)
REAL(GPU_FLUX_PRECISION),DEVICE :: phis(:,:)
REAL(GPU_RES_PRECISION),DEVICE :: xseqD(:,:,:),xseqT(:,:,:),xseqPin(:,:,:)
INTEGER,VALUE :: nPin,rgb,rge,iz,nlv

INTEGER :: ig, ilv, idlv, ifxr, ibfxr, iefxr, jfxr, idxfxr, ifsr, ibfsr, iefsr, idxfsr, offnfsr
INTEGER :: ipin
REAL(8) :: volpnt, vol, phi, maclp, maclv, volsum, phisum, maclpsum, maclvsum

idlv = threadIdx%x+FSP_BLOCK_GDIM*(blockIdx%x-1)
ig = (idlv-1)/nlv
ilv = idlv-ig*nlv
ig = ig+rgb
IF (ig.GT.rge) RETURN
ipin = threadIdx%y+FSP_BLOCK_RDIM*(blockIdx%y-1)
IF (ipin.GT.nPin) RETURN
IF (.NOT.lPinRes(ipin)) THEN
  xseqPin(ilv,ig,ipin) = 0.;
  RETURN
END IF

offnfsr = (iz-1)*nfsr

volsum = 0.; phisum = 0.; maclpsum = 0.; maclvsum = 0.;

!Depl Fxrs
ibfxr = GenStD(ipin); iefxr = GenStD(ipin+1)-1
DO ifxr = ibfxr,iefxr
  jfxr = mapGR2GD(ifxr)
  idxfxr = mapfxrD(jfxr)-(iz-1)*nFxr
  IF (idxfxr.LT.1) RETURN
  IF (idxfxr.GT.nFxr) RETURN
  ibfsr = fsrstD(jfxr); iefsr = nfsrD(jfxr); iefsr = ibfsr+iefsr-1
  idxfsr = mapfsrD(ibfsr)-offnfsr
  maclp = Siglp(idlv,idxfxr); maclv = xst(idlv,idxfsr)-maclp
  vol = 0.; phi = 0.;
  DO ifsr = ibfsr, iefsr
    idxfsr = mapfsrD(ifsr)-offnfsr
    volpnt = volD(ifsr)
    phi = phi+phis(idlv,idxfsr)*volpnt
!    print*, phis(idlv,idxfsr)
    vol = vol+volpnt
!    IF (idlv.EQ.1) print*, ipin, ifxr-ibfxr, volpnt
  END DO
  maclpsum = maclpsum + maclp*phi; maclvsum = maclvsum + maclv*phi;
  phisum = phisum+phi; volsum = volsum+vol
  phi = phi/vol;
!  print*, phi
  xseqD(ilv,ig,ifxr) = -maclp + maclv * phi/(1.0-phi)
  IF (abs(phi-1.0).LT.1.e-10) xseqD(ilv,ig,ifxr) = 1.e+10;
END DO
!Trivial Fxrs
ibfxr = GenStT(ipin); iefxr = GenStT(ipin+1)-1
DO ifxr = ibfxr,iefxr
  jfxr = mapGR2GT(ifxr)
  idxfxr = mapfxrT(jfxr)-(iz-1)*nFxr
  IF (idxfxr.LT.1) RETURN
  IF (idxfxr.GT.nFxr) RETURN
  ibfsr = fsrstT(jfxr); iefsr = nfsrT(jfxr); iefsr = ibfsr+iefsr-1
  idxfsr = mapfsrT(ibfsr)-offnfsr
  maclp = Siglp(idlv,idxfxr); maclv = xst(idlv,idxfsr)-maclp
  vol = 0.; phi = 0.;
  DO ifsr = ibfsr, iefsr
    idxfsr = mapfsrT(ifsr)-offnfsr
    volpnt = volT(ifsr)
    phi = phi+phis(idlv,idxfsr)*volpnt
    vol = vol+volpnt
  END DO
  maclpsum = maclpsum + maclp*phi; maclvsum = maclvsum + maclv*phi;
  phisum = phisum+phi; volsum = volsum+vol
  phi = phi/vol;
  xseqT(ilv,ig,ifxr) = -maclp + maclv * phi/(1.0-phi)
  IF (abs(phi-1.0).LT.1.e-10) xseqT(ilv,ig,ifxr) = 1.e+10;
END DO

IF(maclvsum.EQ.0.) RETURN
maclpsum = maclpsum/phisum
maclvsum = maclvsum/phisum
!IF (ilv.EQ.8) print*, ig*1000+ipin*10, phisum, volsum
phisum = phisum/volsum
xseqPin(ilv,ig,ipin) = -maclpsum+maclvsum*phisum/(1.-phisum)
IF (abs(phisum-1.0).LT.1.e-10) xseqPin(ilv,ig,ipin) = 1.e+10;

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuEquipXsGen1GA_MLG(mapfxrT,mapfsrT,nfsrT,fsrstT,volT,mapA2GT,&
                                            AICStT,lCellRes,lCellAIC,&
                                            Siglp,xst,phis,xseqT,xseqPin,&
                                            iz,nPin,nlv)
USE CUDAFOR
IMPLICIT NONE
INTEGER,DEVICE,DIMENSION(:) :: mapfxrT,mapfsrT,nfsrT,fsrstT
INTEGER,DEVICE,DIMENSION(:) :: mapA2GT,AICStT
LOGICAL,DEVICE,DIMENSION(:) :: lCellRes,lCellAIC
REAL(GPU_XS_PRECISION),DEVICE :: Siglp(:,:), volT(:)
REAL(GPU_PRECISION),DEVICE :: xst(:,:)
REAL(GPU_FLUX_PRECISION),DEVICE :: phis(:,:)
REAL(GPU_RES_PRECISION),DEVICE :: xseqT(:,:),xseqPin(:,:)
INTEGER,VALUE :: iz,nPin,rgb,rge,nlv

INTEGER :: idlv, ifxr, ibfxr, iefxr, jfxr, idxfxr, ifsr, ibfsr, iefsr, idxfsr, offnfsr
INTEGER :: ipin
REAL(8) :: volpnt, vol, phi, maclp, maclv, volsum, phisum, maclpsum, maclvsum

idlv = threadIdx%x+FSP_BLOCK_GDIM*(blockIdx%x-1)
IF (idlv.GT.nlv) RETURN
ipin = threadIdx%y+FSP_BLOCK_RDIM*(blockIdx%y-1)
IF (ipin.GT.nPin) RETURN
IF (.NOT.(lCellRes(ipin).AND.lCellAIC(ipin))) THEN
  xseqPin(idlv,ipin) = 0.;
  RETURN
END IF
offnfsr = (iz-1)*nfsr

volsum = 0.; phisum = 0.; maclpsum = 0.; maclvsum = 0.;

!Trivial Fxrs
ibfxr = AICStT(ipin); iefxr = AICStT(ipin+1)-1
DO ifxr = ibfxr,iefxr
  jfxr = mapA2GT(ifxr)
  idxfxr = mapfxrT(jfxr)-(iz-1)*nFxr
  IF (idxfxr.LT.1) RETURN
  IF (idxfxr.GT.nFxr) RETURN
  ibfsr = fsrstT(jfxr); iefsr = nfsrT(jfxr); iefsr = ibfsr+iefsr-1
  idxfsr = mapfsrT(ibfsr)-offnfsr
  maclp = Siglp(idlv,idxfxr); maclv = xst(idlv,idxfsr)-maclp
  vol = 0.; phi = 0.;
  DO ifsr = ibfsr, iefsr
    idxfsr = mapfsrT(ifsr)-offnfsr
    volpnt = volT(ifsr)
    phi = phi+phis(idlv,idxfsr)*volpnt
    vol = vol+volpnt
  END DO
  maclpsum = maclpsum + maclp*phi; maclvsum = maclvsum + maclv*phi;
  phisum = phisum+phi; volsum = volsum+vol
  phi = phi/vol;
  xseqT(idlv,ifxr) = -maclp + maclv * phi/(1.0-phi)
  IF (abs(phi-1.0).LT.1.e-10) xseqT(idlv,ifxr) = 1.e+10
END DO

IF(maclvsum.EQ.0.) RETURN
maclpsum = maclpsum/phisum
maclvsum = maclvsum/phisum
phisum = phisum/volsum
xseqPin(idlv,ipin) = -maclpsum+maclvsum*phisum/(1.-phisum)
IF (abs(phisum-1.).LT.1.e-10) xseqPin(idlv,ipin) = 1.e+10

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuEquipXsGen1GC_MLG(mapfxrD,mapfsrD,nfsrD,fsrstD,volD,mapC2GD,&
                                            mapfxrT,mapfsrT,nfsrT,fsrstT,volT,mapC2GT,&
                                            ipinD,ipinT,lCellRes,&
                                            Siglp,xst,phis,xseqD,xseqT,&
                                            iz,nCldD,nCldT,nlv)
USE CUDAFOR
IMPLICIT NONE
INTEGER,DEVICE,DIMENSION(:) :: mapfxrD,mapfxrT,mapfsrD,mapfsrT,nfsrD,nfsrT,fsrstD,fsrstT
INTEGER,DEVICE,DIMENSION(:) :: mapC2GD,mapC2GT,ipinD,ipinT
LOGICAL,DEVICE,DIMENSION(:) :: lCellRes
REAL(GPU_XS_PRECISION),DEVICE :: Siglp(:,:), volD(:),volT(:)
REAL(GPU_PRECISION),DEVICE :: xst(:,:)
REAL(GPU_FLUX_PRECISION),DEVICE :: phis(:,:)
REAL(GPU_RES_PRECISION),DEVICE :: xseqD(:,:),xseqT(:,:)
INTEGER,VALUE :: iz,nCldD,nCldT,nlv

INTEGER :: idlv, ifxr, ifxrt, ibfxr, iefxr, jfxr, idxfxr, ifsr, ibfsr, iefsr, idxfsr, offnfsr
INTEGER :: ipin
REAL(8) :: volpnt, vol, phi, maclp, maclv

idlv = threadIdx%x+FSP_BLOCK_GDIM*(blockIdx%x-1)
IF (idlv.GT.nlv) RETURN
ifxr = threadIdx%y+FSP_BLOCK_RDIM*(blockIdx%y-1)
IF (ifxr.LE.nCldD) THEN
  jfxr = mapC2GD(ifxr)
  ipin = ipinD(jfxr)
ELSE IF (ifxr.LE.(nCldD+nCldT)) THEN
  ifxrt = ifxr-nCldD; jfxr = mapC2GT(ifxrt)
  ipin = ipinT(jfxr)
ELSE
  RETURN
END IF
offnfsr = (iz-1)*nfsr

IF (.NOT.lCellRes(ipin)) RETURN

IF (ifxr.LE.nCldD) THEN
  idxfxr = mapfxrD(jfxr)-(iz-1)*nFxr
  IF (idxfxr.LT.1) RETURN
  IF (idxfxr.GT.nFxr) RETURN
  ibfsr = fsrstD(jfxr); iefsr = nfsrD(jfxr); iefsr = ibfsr+iefsr-1
  idxfsr = mapfsrD(ibfsr)-offnfsr
  maclp = Siglp(idlv,idxfxr); maclv = xst(idlv,idxfsr)-maclp
  vol = 0.; phi = 0.;
  DO ifsr = ibfsr, iefsr
    idxfsr = mapfsrD(ifsr)-offnfsr
    volpnt = volD(ifsr)
    phi = phi+phis(idlv,idxfsr)*volpnt
    vol = vol+volpnt
  END DO
  phi = phi/vol;
  xseqD(idlv,ifxr) = -maclp + maclv * phi/(1.0-phi)
  IF (abs(phi-1.0).LT.1.e-10) xseqD(idlv,ifxr) = 1.e+10
ELSE
  idxfxr = mapfxrT(jfxr)-(iz-1)*nFxr
  IF (idxfxr.LT.1) RETURN
  IF (idxfxr.GT.nFxr) RETURN
  ibfsr = fsrstT(jfxr); iefsr = nfsrT(jfxr); iefsr = ibfsr+iefsr-1
  idxfsr = mapfsrT(ibfsr)-offnfsr
  maclp = Siglp(idlv,idxfxr); maclv = xst(idlv,idxfsr)-maclp
  vol = 0.; phi = 0.;
  DO ifsr = ibfsr, iefsr
    idxfsr = mapfsrT(ifsr)-offnfsr
    volpnt = volT(ifsr)
    phi = phi+phis(idlv,idxfsr)*volpnt
    vol = vol+volpnt
  END DO
  phi = phi/vol;
  xseqT(idlv,ifxrt) = -maclp + maclv * phi/(1.0-phi)
  IF (abs(phi-1.0).LT.1.e-10) xseqT(idlv,ifxrt) = 1.e+10
END IF

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuUpdtFnAdj_Sum(areaD,mapGR2GD,mapnuclD,nisoD,pnumD,&
                                            areaT,mapGR2GT,mapnuclT,AccNt,nisoT,pnumT,&
                                            TempAvg,areaPin,&
                                            GenStD,GenStT,lPinRes,lCellRes,lCellAIC,&
                                            rtempsq,ri_a,FnAdjD,FnAdjT,FnAdjPin,&
                                            ptrTempPot,ptrRTemp,IsoG2R,&
                                            iz,rgb,rge)
USE CUDAFOR
IMPLICIT NONE
INTEGER,DEVICE :: mapnuclD(:,:),mapnuclT(:),AccNt(:),nisoD(:),nisoT(:)
INTEGER,DEVICE,DIMENSION(:) :: mapGR2GD,mapGR2GT,GenStD,GenStT,ptrTempPot,ptrRTemp,IsoG2R
LOGICAL,DEVICE,DIMENSION(:) :: lPinRes,lCellRes,lCellAIC
REAL(GPU_XS_PRECISION),DEVICE,DIMENSION(:) :: areaD,areaT,areaPin,TempAvg, rtempsq
REAL(GPU_PNUM_PRECISION),DEVICE,DIMENSION(:) :: pnumT
REAL(GPU_PNUM_PRECISION),DEVICE :: pnumD(:,:)
REAL(GPU_RES_PRECISION),DEVICE :: ri_a(:,:),FnAdjD(:,:),FnAdjT(:,:),FnAdjPin(:,:)
INTEGER,VALUE :: iz,rgb,rge

INTEGER :: ig, ifxr, ibfxr, iefxr, jfxr
INTEGER :: ipin, iso, miso, liso, idxiso, nRiTemp, itmp
INTEGER :: n1, n2, ibtmp, ibT0

REAL(8) :: tempavgsq, N, Npin, areaF, areaP, ria
REAL(4) :: wgt, x2, x1

ig = threadIdx%x+FSP_BLOCK_GDIM*(blockIdx%x-1)
ig = ig+rgb-1
IF (ig.GT.rge) RETURN
ipin = threadIdx%y+FSP_BLOCK_RDIM*(blockIdx%y-1)
IF ((ipin-1)/nxy.NE.(iz-1)) RETURN
IF ((.NOT.lPinRes(ipin)).OR.(.NOT.lCellRes(ipin)).OR.lCellAIC(ipin)) THEN
  FnAdjPin(ig,ipin) = 0.
  RETURN
END IF
areaP = 0.; Npin = 0.; tempavgsq = sqrt(TempAvg(ipin-(iz-1)*nxy))
!Depl Fxrs
ibfxr = GenStD(ipin); iefxr = GenStD(ipin+1)-1
!print*, ibfxr, iefxr
DO ifxr = ibfxr,iefxr
  jfxr = mapGR2GD(ifxr);
  areaF = areaD(jfxr)

  miso = nisoD(jfxr)
  N = 0.
  DO iso = 1, miso
    idxiso = mapnuclD(iso,jfxr)
    idxiso = IsoG2R(idxiso)
    IF (idxiso.EQ.0) CYCLE

    nRiTemp = ptrRTemp(idxiso+1); ibtmp = ptrRTemp(idxiso);
    nRiTemp = nRiTemp-ibtmp
    ibT0 = ptrTempPot(idxiso+1)-nRiTemp
    n1 = 1; x1 = rtempsq(ibtmp+1)
    DO itmp = 1, nRiTemp
      n2 = itmp; x2 = rtempsq(ibtmp+itmp)
      IF (tempavgsq<x2) EXIT
      n1 = itmp; x1 = x2
    END DO
    IF (n1.EQ.n2) THEN
      ria = ri_a(ibT0+n1,ig)
    ELSE
      wgt = (x2-tempavgsq)/(x2-x1)
      ria = ri_a(ibT0+n1,ig)*wgt+ri_a(ibT0+n2,ig)*(1.-wgt)
    END IF
    N = N+pnumD(iso,jfxr)*ria
!    IF (.NOT.(N.GT.0 .OR. N.LE.0)) print*, idxiso, pnumD(iso,jfxr), ria
  END DO
  FnAdjD(ig,ifxr) = N
!  print*, N
  Npin = Npin+N*areaF
  areaP = areaP+areaF
END DO
!Trivial Fxrs
ibfxr = GenStT(ipin); iefxr = GenStT(ipin+1)-1
DO ifxr = ibfxr,iefxr
  jfxr = mapGR2GT(ifxr);
  areaF = areaT(jfxr)

  miso = AccNt(jfxr)+1; liso = miso-1+nisoT(jfxr)
  N = 0.
  DO iso = miso, liso
    idxiso = mapnuclT(iso)
    idxiso = IsoG2R(idxiso)
    IF (idxiso.EQ.0) CYCLE

    nRiTemp = ptrRTemp(idxiso+1); ibtmp = ptrRTemp(idxiso);
    nRiTemp = nRiTemp-ibtmp; ibT0 = ptrTempPot(idxiso+1)-nRiTemp
    n1 = 1; x1 = rtempsq(ibtmp+1)
    DO itmp = 1, nRiTemp
      n2 = itmp; x2 = rtempsq(ibtmp+itmp)
      IF (tempavgsq<x2) EXIT
      n1 = itmp; x1 = x2
    END DO
    IF (n1.EQ.n2) THEN
      ria = ri_a(ibT0+n1,ig)
    ELSE
      wgt = (x2-tempavgsq)/(x2-x1)
      ria = ri_a(ibT0+n1,ig)*wgt+ri_a(ibT0+n2,ig)*(1.-wgt)
    END IF
    N = N+pnumT(iso)*ria
  END DO
  FnAdjT(ig,ifxr) = N
  Npin = Npin+N*areaF
  areaP = areaP+areaF
END DO
FnAdjPin(ig,ipin) = Npin
IF (ig.EQ.rgb) areaPin(ipin-(iz-1)*nxy) = areaP
END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuUpdtFnAdj_Div(areaPin,GenStD,GenStT,lPinRes,lCellRes,lCellAIC,&
                                            FnAdjD,FnAdjT,FnAdjPin,Nsum,&
                                            iz,rgb,rge)
USE CUDAFOR
IMPLICIT NONE
INTEGER,DEVICE,DIMENSION(:) :: GenStD,GenStT
LOGICAL,DEVICE,DIMENSION(:) :: lPinRes,lCellRes,lCellAIC
REAL(GPU_XS_PRECISION),DEVICE,DIMENSION(:) :: areaPin
REAL(GPU_RES_PRECISION),DEVICE :: FnAdjD(:,:),FnAdjT(:,:),FnAdjPin(:,:),Nsum(:)
INTEGER,VALUE :: iz,rgb,rge

INTEGER :: ig, jg, ifxr, ibfxr, iefxr
INTEGER :: ipin

REAL(8) :: areaP, NsumG
REAL(4) :: temp

jg = threadIdx%x+FSP_BLOCK_GDIM*(blockIdx%x-1)
ig = jg+rgb-1
IF (ig.GT.rge) RETURN
ipin = threadIdx%y+FSP_BLOCK_RDIM*(blockIdx%y-1)
IF ((ipin-1)/nxy.NE.(iz-1)) RETURN
IF ((.NOT.lPinRes(ipin)).OR.(.NOT.lCellRes(ipin)).OR.lCellAIC(ipin)) RETURN

areaP = areaPin(ipin-(iz-1)*nxy); NsumG = Nsum(jg)
!Depl Fxrs
ibfxr = GenStD(ipin); iefxr = GenStD(ipin+1)-1
DO ifxr = ibfxr,iefxr
!  print*, ipin, ifxr, jfxr, ig
  FnAdjD(ig,ifxr) = FnAdjD(ig,ifxr)/NsumG
  temp = FnAdjD(ig,ifxr)
END DO
!Trivial Fxrs
ibfxr = GenStT(ipin); iefxr = GenStT(ipin+1)-1
DO ifxr = ibfxr,iefxr
  FnAdjT(ig,ifxr) = FnAdjT(ig,ifxr)/NsumG
END DO
FnAdjPin(ig,ipin) = FnAdjPin(ig,ipin)/NsumG/areaP

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuUpdtFtAdj(mapGR2GD,mapnuclD,nisoD,pnumD,indD,ipinD,&
                                        mapGR2GT,mapnuclT,AccNt,nisoT,pnumT,ipinT,&
                                        TempAvg,tempD,tempT,&
                                        lCellRes,lCellAIC,&
                                        rtempsq,ri_a,sig0sq,lamsigp,&
                                        ptrTempPot,ptrRTemp,ptrSig0,IsoG2R,&
                                        FtAdjD,FtAdjT,xseqD,xseqT,&
                                        rgb,rge,iz,nGenD,nGenTri,nlv)
USE CUDAFOR
IMPLICIT NONE
INTEGER,DEVICE :: mapnuclD(:,:),mapnuclT(:),AccNt(:),nisoD(:),nisoT(:),ipinD(:),ipinT(:)
INTEGER,DEVICE,DIMENSION(:) :: mapGR2GD,mapGR2GT,GenStD,GenStT,ptrTempPot,ptrRTemp,ptrSig0,IsoG2R
LOGICAL,DEVICE,DIMENSION(:) :: lPinRes,lCellRes,lCellAIC
REAL(GPU_XS_PRECISION),DEVICE,DIMENSION(:) :: sig0sq,rtempsq
REAL(GPU_XS_PRECISION),DEVICE,DIMENSION(:) :: areaD,areaT,areaPin,TempAvg,tempD,tempT
REAL(GPU_PNUM_PRECISION),DEVICE,DIMENSION(:) :: pnumT
REAL(GPU_PNUM_PRECISION),DEVICE :: pnumD(:,:),indD(:,:)
REAL(GPU_XS_PRECISION),DEVICE :: lamsigp(:,:)
REAL(GPU_RES_PRECISION),DEVICE :: ri_a(:,:),FtAdjD(:,:,:),FtAdjT(:,:,:),xseqD(:,:,:),xseqT(:,:,:)
INTEGER,VALUE :: rgb,rge,iz,nGenD,nGenTri,nlv

INTEGER :: ig, ilv, idlv, ipin, ifxr, ifxrt, jfxr
INTEGER :: iso, miso, liso, idiso
INTEGER :: nRiTemp, i, npot, n1, n2, n3, n4, ibtmp, ib0,ibT0
REAL(8) :: Nreg, Navg, siglp, sigbsq, TempAvgSq, TempSq, ND
REAL(4) :: xseq,ftadj,wgt1,wgt2,x1,x2

idlv = threadIdx%x+FSP_BLOCK_GDIM*(blockIdx%x-1)
ig = (idlv-1)/nlv
ilv = idlv-ig*nlv
ig = ig+rgb;
IF (ig.GT.rge) RETURN
ifxr = threadIdx%y+FSP_BLOCK_RDIM*(blockIdx%y-1)

IF (ifxr.LE.nGenD) THEN
  jfxr = mapGR2GD(ifxr)
  ipin = ipinD(jfxr)
ELSE IF (ifxr.LE.(nGenD+nGenTri)) THEN
  ifxrt = ifxr-nGenTri; jfxr = mapGR2GT(ifxrt)
  ipin = ipinT(jfxr)
ELSE
  RETURN
END IF
IF (.NOT.lCellRes(ipin)) RETURN
IF (lCellAIC(ipin)) RETURN
ipin = ipin - (iz-1)*nxy
IF ((ipin-1)/nxy.NE.0) RETURN

TempAvgSq = SQRT(TempAvg(ipin))

IF (ifxr.LE.nGenD) THEN
  siglp = 0.; Nreg = 0.; Navg = 0.;
  TempSq = tempD(jfxr); TempSq = SQRT(TempSq)
  xseq = xseqD(ilv,ig,ifxr); ftadj = FtAdjD(ilv,ig,ifxr)

  miso = nisoD(jfxr);
  DO iso = 1, miso
    idiso = mapnuclD(iso,jfxr)
    sigbsq = lamsigp(ig,idiso)
    siglp = siglp+pnumD(iso,jfxr)*sigbsq
  END DO
  DO iso = 1, miso
    idiso = mapnuclD(iso,jfxr)
    sigbsq = lamsigp(ig,idiso)
    ND = indD(idiso,jfxr)
    idiso = IsoG2R(idiso)
    IF (idiso .EQ. 0) CYCLE

    sigbsq = (xseq+siglp)/ND/ftadj-sigbsq
    ND = pnumD(iso,jfxr)
!    sigbsq = sigbsq-lamsigp(ig,idiso)
    IF (sigbsq.LE.0.) sigbsq = 1.e-10
    sigbsq = sqrt(sigbsq)
    npot = ptrSig0(idiso+1); nRiTemp = ptrRTemp(idiso+1)

    ib0 = ptrSig0(idiso); ibtmp = ptrRTemp(idiso); ibT0 = ptrTempPot(idiso);
    npot = npot-ib0; nRiTemp = nRiTemp - ibtmp

    !! n1, n2 -> sigbsq
    n1 = 1; x1 = sig0sq(ib0+1)
    DO i = 1, npot
      n2 = i; x2 = sig0sq(ib0+i)
      IF (sigbsq<x2) EXIT
      n1 = i; x1 = x2
    END DO
    IF (n1.EQ.n2) THEN
      wgt1 = 1.
    ELSE
      wgt1 = (x2-sigbsq)/(x2-x1)
      !ria = ri_a(ibT0+n1,ig)*wgt+ri_a(ibT0+n2)*(1.-wgt)
    END IF

    ! n3, n4 -> Tempsq
    n3 = 1; x1 = rtempsq(ibtmp+1)
    DO i = 1, nRiTemp
      n4 = i; x2 = rtempsq(ibtmp+i)
      IF (TempSq<x2) EXIT
      n3 = i; x1 = x2;
    END DO
    IF (n3.EQ.n4) THEN
      wgt2 = 1.
    ELSE
      wgt2 = (x2-TempSq)/(x2-x1)
    END IF
!    ib0 = ibT0+(n3-1)*npot
!    x1 = wgt1*(ri_a(ib0+n1,ig))+(1.-wgt1)*(ri_a(ib0+n2,ig))
!    ib0 = ibT0+(n4-1)*npot
!    x2 = wgt1*(ri_a(ib0+n1,ig))+(1.-wgt1)*(ri_a(ib0+n2,ig))
!    Nreg = Nreg+ND*(x1*wgt2+x2*(1.-wgt2))
    ib0 = ibT0+(n1-1)*nRiTemp
    x1 = wgt2*(ri_a(ib0+n3,ig))+(1.-wgt2)*(ri_a(ib0+n4,ig))
    ib0 = ibT0+(n2-1)*nRiTemp
    x2 = wgt2*(ri_a(ib0+n3,ig))+(1.-wgt2)*(ri_a(ib0+n4,ig))
    Nreg = Nreg+ND*(x1*wgt1+x2*(1.-wgt1))

    ! n3, n4 -> TempAvgsq
    n3 = 1; x1 = rtempsq(ibtmp+1)
    DO i = 1, nRiTemp
      n4 = i; x2 = rtempsq(ibtmp+i)
      IF (TempAvgSq<x2) EXIT
      n3 = i; x1 = x2;
    END DO
    IF (n3.EQ.n4) THEN
      wgt2 = 1.
    ELSE
      wgt2 = (x2-TempAvgSq)/(x2-x1)
    END IF
!    ib0 = ibT0+(n3-1)*npot
!    x1 = wgt1*(ri_a(ib0+n1,ig))+(1.-wgt1)*(ri_a(ib0+n2,ig))
!    ib0 = ibT0+(n4-1)*npot
!    x2 = wgt1*(ri_a(ib0+n1,ig))+(1.-wgt1)*(ri_a(ib0+n2,ig))
!    Navg = Navg+ND*(x1*wgt2+x2*(1.-wgt2))
    ib0 = ibT0+(n1-1)*nRiTemp
    x1 = wgt2*(ri_a(ib0+n3,ig))+(1.-wgt2)*(ri_a(ib0+n4,ig))
    ib0 = ibT0+(n2-1)*nRiTemp
    x2 = wgt2*(ri_a(ib0+n3,ig))+(1.-wgt2)*(ri_a(ib0+n4,ig))
    Navg = Navg+ND*(x1*wgt1+x2*(1.-wgt1))
  END DO
  FtAdjD(ilv,ig,ifxr) = Nreg/Navg
!  IF ((ig*10+ilv).EQ.11) print*, ilv+ig*10, jfxr, Nreg, Navg
ELSE
  siglp = 0.; Nreg = 0.; Navg = 0.;
  TempSq = tempT(jfxr); TempSq = sqrt(TempSq)
  xseq = xseqT(ilv,ig,ifxrt); ftadj = FtAdjT(ilv,ig,ifxrt)

  miso = AccNt(jfxr)+1; liso = miso-1+nisoT(jfxr)
  DO iso = miso, liso
    idiso = mapnuclT(iso)
    sigbsq = lamsigp(ig,idiso)
    siglp = siglp+pnumT(iso)*sigbsq
  END DO
  DO iso = 1, miso
    idiso = mapnuclT(iso)
    sigbsq = lamsigp(ig,idiso)
    idiso = IsoG2R(idiso)
    IF (idiso .EQ. 0) CYCLE
    ND = pnumT(iso)

    sigbsq = (xseq+siglp)/ND/ftadj-sigbsq
!    sigbsq = sigbsq-lamsigp(ig,idiso)
    IF (sigbsq.LE.0.) sigbsq = 1.e-10
    sigbsq = sqrt(sigbsq)
    npot = ptrSig0(idiso+1); nRiTemp = ptrRTemp(idiso+1)

    ib0 = ptrSig0(idiso); ibtmp = ptrRTemp(idiso); ibT0 = ptrTempPot(idiso);
    npot = npot-ib0; nRiTemp = nRiTemp-ibtmp
    n1 = 1; x1 = sig0sq(ib0+1)
    DO i = 1, npot
      n2 = i; x2 = sig0sq(ib0+i)
      IF (sigbsq<x2) EXIT
      n1 = i; x1 = x2
    END DO
    IF (n1.EQ.n2) THEN
      wgt1 = 1.
    ELSE
      wgt1 = (x2-sigbsq)/(x2-x1)
      !ria = ri_a(ibT0+n1,ig)*wgt+ri_a(ibT0+n2)*(1.-wgt)
    END IF

    n3 = 1; x1 = rtempsq(ibtmp+1)
    DO i = 1, nRiTemp
      n4 = i; x2 = rtempsq(ibtmp+i)
      IF (TempSq<x2) EXIT
      n3 = i; x1 = x2;
    END DO
    IF (n3.EQ.n4) THEN
      wgt2 = 1.
    ELSE
      wgt2 = (x2-TempSq)/(x2-x1)
    END IF
    ib0 = ibT0+(n1-1)*nRiTemp
    x1 = wgt2*(ri_a(ib0+n3,ig))+(1.-wgt2)*(ri_a(ib0+n4,ig))
    ib0 = ibT0+(n2-1)*nRiTemp
    x2 = wgt2*(ri_a(ib0+n3,ig))+(1.-wgt2)*(ri_a(ib0+n4,ig))
    Nreg = Nreg+ND*(x1*wgt1+x2*(1.-wgt1))

    n3 = 1; x1 = rtempsq(ibtmp+1)
    DO i = 1, nRiTemp
      n4 = i; x2 = rtempsq(ibtmp+i)
      IF (TempAvgSq<x2) EXIT
      n3 = i; x1 = x2;
    END DO
    IF (n3.EQ.n4) THEN
      wgt2 = 1.
    ELSE
      wgt2 = (x2-TempAvgSq)/(x2-x1)
    END IF
    ib0 = ibT0+(n1-1)*nRiTemp
    x1 = wgt2*(ri_a(ib0+n3,ig))+(1.-wgt2)*(ri_a(ib0+n4,ig))
    ib0 = ibT0+(n2-1)*nRiTemp
    x2 = wgt2*(ri_a(ib0+n3,ig))+(1.-wgt2)*(ri_a(ib0+n4,ig))
    Navg = Navg+ND*(x1*wgt1+x2*(1.-wgt1))
  END DO
  FtAdjT(ilv,ig,ifxrt) = Nreg/Navg
END IF
END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE GetNavg(tarArr, resArr, denom, yb, ye, ny)
USE CUDAFOR
IMPLICIT NONE
REAL(4), DEVICE :: tarArr(:,:), resArr(:)
INTEGER, VALUE :: yb, ye, ny
REAL(8), VALUE :: denom
REAL(4), SHARED :: bank(512)

INTEGER :: ix, iy, jy
INTEGER :: resiN, i, j
REAL(4) :: tempVal

ix = threadIdx%x; iy = blockIdx%y
jy = iy+yb-1
IF(jy.GT.ye) RETURN

resiN = ny/512

! &&& Collecting global memory data &&& !
j = ix; tempVal = 0.;
DO i = 1, resiN
  tempVal = tempVal+tarArr(jy,j)
  j = j+512
END DO
IF (j.LE.ny) THEN
  tempVal = tempVal+tarArr(jy,j)
END IF
bank(ix) = tempVal
! &&& ----------------------------- &&& !
CALL syncthreads()
! &&& Reduction &&& !
resiN = 512
DO WHILE(resiN.GT.1)
  resiN = resiN/2;
  IF (ix .GT. resiN) THEN
    RETURN
  ELSE
    bank(ix) = bank(ix)+bank(ix+resiN)
  END IF
  CALL syncthreads()
END DO
! &&& --------- &&& !

resArr(iy) = bank(ix)/denom

END SUBROUTINE

SUBROUTINE SetPlnLsigP_cuMLG(Siglp,xst,nlv,rgb,rge,iz)
USE PE_Mod, ONLY : PE
IMPLICIT NONE
REAL(GPU_XS_PRECISION), DEVICE :: Siglp(:,:)
REAL(GPU_PRECISION), DEVICE :: xst(:,:)
INTEGER :: nlv, rgb, rge, iz, jz

INTEGER :: isync, ierr

jz = (iz-PE%myzb)+1

CALL cuSetPlnLSigP_MLG<<<blockGF,trdFSP,0,cuDevice%myStream>>>(cuBlockFxr%nisoD,cuBlockFxr%pnumD,cuBlockFxr%mapglobalidD,&
    cuBlockFxr%mapfsrD,cuBlockFxr%nfsrD,cuBlockFxr%fsrstD,&
    cuBlockFxr%AccNt,cuBlockFxr%nisoT,cuBlockFxr%pnumT,&
    cuBlockFxr%mapglobalidT,cuBlockFxr%mapfsrT,cuBlockFxr%nfsrT,cuBlockFxr%fsrstT,&
    cuBlockFxr%MapNuclD,cuBlockFxr%mapnuclT,&
    cuResIsoData%lamsigp,cuResIsoData%sigp,Siglp,xst,&
    rgb,rge,nlv,jz,FxrD%nFxrD,FxrTri%nFxrTri)
!isync = cudaDeviceSynchronize()
!ierr = cudaGetLastError()
!if (ierr.NE.0) print*, cudaGetErrorString(ierr)
!if (ierr.NE.0) print*, __FILE__, __LINE__
!print*, 'SetPlnSigP'
!print*, blockGF
!print*, trdFSP

CALL cuSetPlnXST_MLG<<<blockGF,trdFSP,0,cuDevice%myStream>>>(cuBlockFxr%mapglobalidD,&
    cuBlockFxr%mapfsrD,cuBlockFxr%nfsrD,cuBlockFxr%fsrstD,&
    cuBlockFxr%mapglobalidT,cuBlockFxr%mapfsrT,cuBlockFxr%nfsrT,cuBlockFxr%fsrstT,&
    cuBlockFxr%mapGR2GD,cuBlockFxr%mapGR2GT,&
    cuMLGData%f_maclv,Siglp,xst,&
    cuBlockFxr%FnAdjD,cuBlockFxr%FtAdjD,cuBlockFxr%FnAdjT,cuBlockFxr%FtAdjT,&
    jz,rgb,rge,nlv,FxrD%nGenD,FxrTri%nGenTri)
!isync = cudaDeviceSynchronize()
!ierr = cudaGetLastError()
!if (ierr.NE.0) print*, cudaGetErrorString(ierr)
!if (ierr.NE.0) print*, __FILE__, __LINE__
!print*, 'SetPlnXST'

END SUBROUTINE

SUBROUTINE SetPlnLsigP1G_cuMLG(Siglp,xst,nlv,lCLD,iz)
USE PE_Mod, ONLY : PE
IMPLICIT NONE
REAL(GPU_XS_PRECISION), DEVICE :: Siglp(:,:)
REAL(GPU_PRECISION), DEVICE :: xst(:,:)
INTEGER :: nlv, iz,jz
LOGICAL :: lCLD
TYPE(dim3) :: block1G

!REAL(8), POINTER :: SiglpHost(:,:)
!REAL(4), POINTER :: xstHost(:,:)
!INTEGER :: nfsr, nfxr, i, ilv
INTEGER :: isync, ierr

block1G = blockAF
IF (lCLD) block1G = blockCF

jz = (iz-PE%myzb)+1

CALL cuSetPlnLSigP1G_MLG<<<block1G,trdFSP,0,cuDevice%myStream>>>(cuBlockFxr%nisoD,cuBlockFxr%pnumD,cuBlockFxr%mapglobalidD,&
    cuBlockFxr%mapfsrD,cuBlockFxr%nfsrD,cuBlockFxr%fsrstD,&
    cuBlockFxr%AccNt,cuBlockFxr%nisoT,cuBlockFxr%pnumT,&
    cuBlockFxr%mapglobalidT,cuBlockFxr%mapfsrT,cuBlockFxr%nfsrT,cuBlockFxr%fsrstT,&
    cuBlockFxr%MapNuclD,cuBlockFxr%mapnuclT,&
    cuResIsoData%sigp,Siglp,xst,&
    nlv,jz,FxrD%nFxrD,FxrTri%nFxrTri)
!isync = cudaDeviceSynchronize()
!ierr = cudaGetLastError()
!if (ierr.NE.0) print*, cudaGetErrorString(ierr)
!if (ierr.NE.0) print*, __FILE__, __LINE__
!print*, 'SetPlnSigP1G'
CALL cuSetPlnXST1G_MLG<<<block1G,trdFSP,0,cuDevice%myStream>>>(cuBlockFxr%mapglobalidD,&
    cuBlockFxr%mapfsrD,cuBlockFxr%nfsrD,cuBlockFxr%fsrstD,&
    cuBlockFxr%mapglobalidT,cuBlockFxr%mapfsrT,cuBlockFxr%nfsrT,cuBlockFxr%fsrstT,&
    cuBlockFxr%mapC2GD,cuBlockFxr%mapC2GT,cuBlockFxr%mapA2GT,&
    cuMLGData%f_maclv1G,cuMLGData%c_maclv1G,Siglp,xst,&
    jz,nlv,lCLD,FxrD%nCldD,FxrTri%nCldTri,FxrTri%nAICTri)
!isync = cudaDeviceSynchronize()
!ierr = cudaGetLastError()
!if (ierr.NE.0) print*, cudaGetErrorString(ierr)
!if (ierr.NE.0) print*, __FILE__, __LINE__
!print*, 'SetPlnXST1G'

!nfxr = FxrD%nFxrD+FxrTri%nFxrTri; nfsr = FxrD%ntfsr+FxrTri%ntfsr
!ALLOCATE(SiglpHost(nlv,nfxr))
!ALLOCATE(xstHost(nlv,nfsr))
!SiglpHost = Siglp; xstHost = xst
!print*, 'Siglp', lCLD
!DO ilv = 1, nlv
!  DO i = 1, nfxr
!    print*, ilv, i, SiglpHost(ilv,i)
!  END DO
!END DO
!
!print*, 'xst', lCLD
!DO ilv = 1, nlv
!  DO i = 1, nfsr
!    print*, ilv, i, xstHost(ilv,i)
!  END DO
!END DO

!DEALLOCATE(SiglpHost,xstHost)
END SUBROUTINE

SUBROUTINE SetSubGrpSrc_cuMLG(Siglp,src,xst,nlv,iz,matcase)
USE PE_Mod, ONLY : PE
IMPLICIT NONE
REAL(GPU_XS_PRECISION), DEVICE :: Siglp(:,:)
REAL(GPU_SOURCE_PRECISION), DEVICE :: src(:,:)
REAL(GPU_PRECISION), DEVICE :: xst(:,:)
INTEGER :: nlv,iz,matcase

INTEGER :: jz
TYPE(dim3) :: blocksrc

INTEGER :: ierr, isync

!REAL(4), POINTER :: srcHost(:,:)
!INTEGER :: ig, i, nfsr

jz = (iz-PE%myzb)+1

SELECT CASE(matcase)
CASE(matGen)
  blocksrc = blockGF
CASE(matCLD)
  blocksrc = blockCF
CASE(matAIC)
  blocksrc = blockAF
END SELECT

CALL cuSetSubGrpSrc<<<blocksrc,trdFSP,0,cuDevice%myStream>>>(cuBlockFxr%mapglobalidD,cuBlockFxr%mapfsrD,&
    cuBlockFxr%nfsrD,cuBlockFxr%fsrstD,&
    cuBlockFxr%mapglobalidT,cuBlockFxr%mapfsrT,&
    cuBlockFxr%nfsrT,cuBlockFxr%fsrstT,&
    Siglp,src,xst,&
    jz,FxrD%nFxrD,FxrTri%nFxrTri,nlv)
!isync = cudaDeviceSynchronize()
!ierr = cudaGetLastError()
!if (ierr.NE.0) print*, cudaGetErrorString(ierr)
!if (ierr.NE.0) print*, __FILE__, __LINE__
!print*, 'SetSubGrpSrc', matcase

!nfsr = FxrD%ntfsr+FxrTri%ntfsr
!!ALLOCATE(srcHost(nlv,nfsr)); srcHost = src
!print*, 'src'
!DO ig = 1, nlv
!  DO i = 1, nfsr
!!    print*, ig, i, srcHost(ig,i)
!  END DO
!END DO
END SUBROUTINE

SUBROUTINE EquipXsGen_cuMLG(Siglp,xst,phis,rgb,rge,nlv,iz)
USE PE_Mod, ONLY : PE
USE geom, ONLY : core
IMPLICIT NONE
REAL(GPU_XS_PRECISION), DEVICE :: Siglp(:,:)
REAL(GPU_FLUX_PRECISION), DEVICE :: phis(:,:)
REAL(GPU_PRECISION), DEVICE :: xst(:,:)
INTEGER :: rgb,rge,nlv,iz

INTEGER :: jz, nPin

INTEGER :: ierr, isync

!REAL(4), POINTER :: xseqHost(:,:,:), xseqPin(:,:,:)
!INTEGER :: ig, ilv, i

jz = (iz-PE%myzb)+1
nPin = core%nxy*(PE%myze-PE%myzb+1)

CALL cuEquipXsGen_MLG<<<blockGP,trdFSP,0,cuDevice%myStream>>>(cuBlockFxr%mapglobalidD,cuBlockFxr%mapfsrD,cuBlockFxr%nfsrD,&
    cuBlockFxr%fsrstD,cuBlockFxr%fsrvolD,cuBlockFxr%mapGR2GD,&
    cuBlockFxr%mapglobalidT,cuBlockFxr%mapfsrT,cuBlockFxr%nfsrT,&
    cuBlockFxr%fsrstT,cuBlockFxr%fsrvolT,cuBlockFxr%mapGR2GT,&
    cuPinInfo%GenStD,cuPinInfo%GenStT,cuResPin%lPinRes,&
    Siglp,xst,phis,&
    cuBlockFxr%xseq_f_mgD,cuBlockFxr%xseq_f_mgT,cuResPin%avgxseq_mg,&
    nPin,rgb,rge,jz,nlv)
!isync = cudaDeviceSynchronize()
!ierr = cudaGetLastError()
!if (ierr.NE.0) print*, cudaGetErrorString(ierr)
!if (ierr.NE.0) print*, __FILE__, __LINE__
!print*, 'EquipXsGen'

!ALLOCATE(xseqHost(nlv,rgb:rge,FxrD%nGenD), xseqPin(nlv,rgb:rge,nPin))
!xseqHost = cuBlockFxr%xseq_f_mgD; xseqPin = cuResPin%avgxseq_mg
!print*, 'xseq_mgF'
!DO ig = rgb, rge
!  DO i = 1, FxrD%nGenD
!    DO ilv = 1, nlv
!!      print*, ig, ilv, i, xseqHost(ilv,ig,i)
!    END DO
!  END DO
!END DO
!print*, 'xseq_mgPin'
!DO ig = rgb, rge
!  DO i = 1, nPin
!    DO ilv = 1, nlv
!      IF (xseqPin(ilv,ig,i).EQ.0) CYCLE
!      print*, ig, ilv, i, xseqPin(ilv,ig,i)
!    END DO
!  END DO
!END DO
!stop
END SUBROUTINE

SUBROUTINE EquipXsGen1G_cuMLG(Siglp,xst,phis,nlv,iz,matcase)
USE PE_Mod, ONLY : PE
USE geom, ONLY : core
IMPLICIT NONE
REAL(GPU_XS_PRECISION), DEVICE :: Siglp(:,:)
REAL(GPU_PRECISION), DEVICE :: xst(:,:)
REAL(GPU_FLUX_PRECISION), DEVICE :: phis(:,:)
INTEGER :: nlv,iz,matcase

INTEGER :: nPin,jz

INTEGER :: ierr, isync

!REAL(4), POINTER :: xseqT(:,:)
!INTEGER :: i, ilv

jz = (iz-PE%myzb)+1
nPin = core%nxy*(PE%myze-PE%myzb+1)

IF (matcase.EQ.matAIC) THEN
  CALL cuEquipXsGen1GA_MLG<<<blockAP,trdFSP,0,cuDevice%myStream>>>(cuBlockFxr%mapglobalidT,cuBlockFxr%mapfsrT,&
    cuBlockFxr%nfsrT,cuBlockFxr%fsrstT,cuBlockFxr%fsrvolT,cuBlockFxr%mapA2GT,&
    cuPinInfo%AICStT,cuResPin%lCellRes,cuResPin%lAIC,&
    Siglp,xst,phis,&
    cuBlockFxr%xseq_f_1gT,cuResPin%avgxseq_1g,&
    jz,nPin,nlv)
ELSE IF (matcase.EQ.matCLD) THEN
  CALL cuEquipXsGen1GC_MLG<<<blockCF,trdFSP,0,cuDevice%myStream>>>(cuBlockFxr%mapglobalidD,cuBlockFxr%mapfsrD,&
    cuBlockFxr%nfsrD,cuBlockFxr%fsrstD,cuBlockFxr%fsrvolD,cuBlockFxr%mapC2GD,&
    cuBlockFxr%mapglobalidT,cuBlockFxr%mapfsrT,cuBlockFxr%nfsrT,&
    cuBlockFxr%fsrstT,cuBlockFxr%fsrvolT,cuBlockFxr%mapC2GT,&
    cuBlockFxr%pinidD,cuBlockFxr%pinidT,cuResPin%lCellRes,&
    Siglp,xst,phis,&
    cuBlockFxr%xseq_c_1gD,cuBlockFxr%xseq_c_1gT,&
    jz,FxrD%nCldD,FxrTri%nCldTri,nlv)
ELSE
  print*, 'WRONG MAT SPECIFIER @ EquipXSGen1G'
END IF
!isync = cudaDeviceSynchronize()
!ierr = cudaGetLastError()
!if (ierr.NE.0) print*, cudaGetErrorString(ierr)
!if (ierr.NE.0) print*, __FILE__, __LINE__
!print*, 'EquipXsGen1G', matcase
!print*, 'xseqT'
!ALLOCATE(xseqT(nlv,FxrTri%nCldTri)); xseqT = cuBlockFxr%xseq_c_1gT
!DO ilv = 1, nlv
!  DO i = 1, FxrTri%nCldTri
!    print*, ilv, i, xseqT(ilv,i)
!  END DO
!END DO
!DEALLOCATE(xseqT)
END SUBROUTINE

SUBROUTINE UpdtFnAdj_cuMLG(iz, rgb, rge)
USE PE_Mod, ONLY : PE
USE geom, ONLY : core
IMPLICIT NONE
INTEGER :: iz,rgb,rge
REAL(GPU_XS_PRECISION), DEVICE, ALLOCATABLE :: TempAvg(:), areaPin(:)
REAL(8) :: areaSum
REAL(GPU_RES_PRECISION), DEVICE, ALLOCATABLE :: Nsum(:)
INTEGER :: ierrblas, ierr, isync
INTEGER :: i, jz, ig, nxypin, nrg

INTEGER :: istat

jz = (iz-PE%myzb)+1
nxypin = core%nxy
nrg = rge-rgb+1

ALLOCATE(TempAvg(nxypin),areaPin(nxypin),Nsum(nrg))
!TempAvg = ResPin%temp; areaPin = 0.;
istat = cudaMemCpy(TempAvg,ResPin%temp((jz-1)*nxypin+1:jz*nxypin),nxypin)
istat = cudaMemSet(areaPin,0.,nxypin)

CALL cuUpdtFnAdj_Sum<<<blockGP,trdFSP,0,cuDevice%myStream>>>(cuBlockFxr%areaD,cuBlockFxr%mapGR2GD,&
  cuBlockFxr%MapNuclD,cuBlockFxr%nisoD,cuBlockFxr%pnumD,&
  cuBlockFxr%areaT,cuBlockFxr%mapGR2GT,cuBlockFxr%mapnuclT,&
  cuBlockFxr%AccNt,cuBlockFxr%nisoT,cuBlockFxr%pnumT,&
  TempAvg,areaPin,&
  cuPinInfo%GenStD,cuPinInfo%GenStT,&
  cuResPin%lPinRes,cuResPin%lCellRes,cuResPin%lAIC,&
  cuResIsoData%rtempsq,cuResIsoData%ri_a,&
  cuBlockFxr%FnAdjD,cuBlockFxr%FnAdjT,cuResPin%FnAdj,&
  cuResIsoData%ptrTempPot,cuResIsoData%ptrRTemp,cuResIsoData%IsoG2R,&
  jz,rgb,rge)
!isync = cudaDeviceSynchronize()
!ierr = cudaGetLastError()
!if (ierr.NE.0) print*, cudaGetErrorString(ierr)
!if (ierr.NE.0) print*, __FILE__, __LINE__
!print*, 'UpdtFnAdj Sum'

! &&& Sum areaPin & FnAdjPin &&&
ierrblas = cublasDasum_v2(cuDevice%myblasHandle,nxypin,areaPin,1,areaSum)
CALL GetNavg<<<blockN,trdN,512*4,cuDevice%myStream>>>(cuResPin%FnAdj(:,(jz-1)*nxypin+1:jz*nxypin),Nsum,areaSum,rgb,rge,nxypin)
!ierrblas = cublasSasum_v2(cuDevice%myblasHandle,nPin,FnAdj1D,1,Nsum)
!print*, 'SIGNAL',areaSum,Nsum,(Nsum/areaSum)
!Nsum = Nsum/areaSum
! &&& ---------------------- &&&
CALL cuUpdtFnAdj_Div<<<blockGP,trdFSP,0,cuDevice%myStream>>>(areaPin,&
  cuPinInfo%GenStD,cuPinInfo%GenStT,&
  cuResPin%lPinRes,cuResPin%lCellRes,cuResPin%lAIC,&
  cuBlockFxr%FnAdjD,cuBlockFxr%FnAdjT,cuResPin%FnAdj,Nsum,&
  jz,rgb,rge)
!isync = cudaDeviceSynchronize()
!ierr = cudaGetLastError()
!if (ierr.NE.0) print*, cudaGetErrorString(ierr)
!if (ierr.NE.0) print*, __FILE__, __LINE__
!print*, 'UpdtFnAdj Div'

DEALLOCATE(TempAvg,areaPin,Nsum)

END SUBROUTINE

SUBROUTINE UpdtFtAdj_cuMLG(rgb,rge,nlv,iz)
USE PE_Mod, ONLY : PE
USE geom, ONLY : core
IMPLICIT NONE
INTEGER :: rgb, rge, nlv, iz
REAL(GPU_XS_PRECISION), DEVICE, ALLOCATABLE :: TempAvg(:)
!REAL(4), POINTER :: FtadjD(:,:,:)
!INTEGER :: ilv, ig, i
INTEGER :: nxypin, jz
INTEGER :: ierr, isync, istat

jz = (iz-PE%myzb)+1
nxypin = core%nxy

ALLOCATE(TempAvg(nxypin))
!TempAvg = ResPin%temp
istat = cudaMemCpy(TempAvg,ResPin%temp((jz-1)*nxypin+1:jz*nxypin),nxypin)
!print*, ResPin%temp

!print*, 'FtAdj'
CALL cuUpdtFtAdj<<<blockGF,trdFSP,0,cuDevice%myStream>>>(cuBlockFxr%mapGR2GD,cuBlockFxr%MapNuclD,&
  cuBlockFxr%nisoD,cuBlockFxr%pnumD,cuBlockFxr%indD,cuBlockFxr%pinidD,&
  cuBlockFxr%mapGR2GT,cuBlockFxr%MapNuclT,cuBlockFxr%AccNt,&
  cuBlockFxr%nisoT,cuBlockFxr%pnumT,cuBlockFxr%pinidT,&
  TempAvg,cuBlockFxr%tempD,cuBlockFxr%tempT,&
  cuResPin%lCellRes,cuResPin%lAIC,&
  cuResIsoData%rtempsq,&
  cuResIsoData%ri_a,cuResIsoData%sig0sq,cuResIsoData%lamsigp,&
  cuResIsoData%ptrTempPot,cuResIsoData%ptrRTemp,cuResIsoData%ptrSig0,cuResIsoData%IsoG2R,&
  cuBlockFxr%FtAdjD,cuBlockFxr%FtAdjT,cuBlockFxr%xseq_f_mgD,cuBlockFxr%xseq_f_mgT,&
  rgb,rge,jz,FxrD%nGenD,FxrTri%nGenTri,nlv)
!isync = cudaDeviceSynchronize()
!ierr = cudaGetLastError()
!if (ierr.NE.0) print*, cudaGetErrorString(ierr)
!if (ierr.NE.0) print*, __FILE__, __LINE__
!print*, 'UpdtFtAdj'

DEALLOCATE(TempAvg)

!ALLOCATE(FtadjD(nlv,nrg,FxrD%nGenD)); FtadjD = cuBlockFxr%FtAdjD
!DO ig = 1, nrg
!  DO i = 1, FxrD%nGenD
!    DO ilv = 1, nlv
!      print*, ig, ilv, i, FtadjD(ilv,ig,i)
!    END DO
!  END DO
!END DO
!
!stop

END SUBROUTINE

FUNCTION SubGrpFSPErr_cuMLG(phis, phisd, sizephi) RESULT(errmax)
IMPLICIT NONE
REAL(GPU_FLUX_PRECISION), DEVICE, TARGET :: phis(:,:), phisd(:,:)
REAL(GPU_FLUX_PRECISION), DEVICE, POINTER :: phis1D(:), phisd1D(:)
INTEGER :: sizephi
REAL :: errmax

INTEGER :: ierrblas, imax
REAL(4) :: neg1 = -1, errmax4

! &&& Max (phis-phisd) &&&
phis1D => phis(1:,1); phisd1D => phisd(1:,1)
! 1. phisd = phisd-phis
ierrblas = cublasSaxpy_v2(cuDevice%myblasHandle,sizephi,neg1,phis1D,1,phisd1D,1)
! 2. reduction:max
ierrblas = cublasIsamax_v2(cuDevice%myblasHandle,sizephi,phisd1D,1,imax)
errmax4 = phisd1D(imax)
! 3. phisd = phis
phisd = phis
! &&& ---------------- &&&
errmax = abs(errmax4)
END FUNCTION

ATTRIBUTES(GLOBAL) SUBROUTINE cuEffRIFPin(lclad,ityp,IsoG2R,&
                                ptrRTemp,ptrSig0,ptrNlv,ptrRifRat,ptrTempPot,ptrTempLv,ptrRifRx, &
                                lamsigp,lvabs,lvfis,wgtabs,wgtfis,f_maclv_log,f_maclv1G_log,rtempsq,ratlog, &
                                xsalog,rifabs,riffis,FnAdj,xseq_1g,lCellAIC,lPinRes,rifa,riff,xseq_mg, &
                                niso,mapnucl,rifid,pnum,ind,Temp,ntiso,nrg,nmaclv,nmaclv1G,iz)
USE CUDAFOR
IMPLICIT NONE
LOGICAL,DEVICE :: lclad(:)
INTEGER,DEVICE,DIMENSION(:) :: ityp,IsoG2R,ptrRTemp,ptrSig0,ptrNlv,ptrRifRat,ptrTempPot,ptrTempLv,ptrRifRx
REAL(8),DEVICE,DIMENSION(:,:) :: lamsigp,lvabs,lvfis,wgtabs,wgtfis, f_maclv_log, f_maclv1G_log
REAL(8),DEVICE,DIMENSION(:) :: rtempsq, ratlog
REAL(4),DEVICE,DIMENSION(:,:) :: xsalog,rifabs,riffis,FnAdj,xseq_1g
LOGICAL,DEVICE :: lCellAIC(:), lPinRes(:)
REAL(4),DEVICE,DIMENSION(:,:,:) :: rifa, riff, xseq_mg
INTEGER,DEVICE :: niso(:),mapnucl(:,:),rifid(:,:)
REAL(8),DEVICE :: pnum(:,:),ind(:,:),Temp(:)
INTEGER,VALUE :: ntiso, nrg, nmaclv, nmaclv1G, iz

INTEGER :: ipin, ig, iso, idxiso, ireso, jdxiso, jreso
INTEGER :: nisoP, nrtemp, nsig0
REAL(8) :: xut8, xl8, xr8
REAL(4) :: xut4, xl4, xr4
INTEGER :: T1, T2, Mc1, Mc2, AL1, AL2, AR1, AR2, R1, R2
REAL(8) :: wT, wMc, wR
REAL(4) :: wAL, wAR
INTEGER :: iut, ibt, ilv, jso, ibTL, LB, UB

REAL(8) :: siglpiso, ND, adjintmlg, xsa, xsf, phi, lvX, sigb, wgtX
REAL(8) :: rifaT(2), riffT(2)
LOGICAL :: IsFis

ipin = threadIdx%x+(blockIdx%x-1)*EFF_BLOCK_RDIM
IF ((ipin-1)/nxy.NE.(iz-1)) RETURN
IF (.NOT.lPinRes(ipin)) RETURN
iso = threadIdx%y+(blockIdx%y-1)*EFF_BLOCK_IDIM
IF (iso.GT.ntiso) RETURN
ig = blockIdx%z
IF (ig.GT.nrg) RETURN

nisoP = niso(ipin)
IF (iso.GT.nisoP) RETURN
idxiso = mapnucl(iso,ipin)
ireso = IsoG2R(idxiso)
IF (ireso.EQ.0) RETURN
IF (lclad(idxiso)) RETURN
ND = ind(idxiso,ipin)
IF (ND.LE.0.) RETURN

! &&& ========= AvgTemperature Weights Finding ========= &&& !
xut8 = Temp(ipin); xut8 = sqrt(xut8)
ibt = ptrRTemp(ireso); nrtemp = ptrRTemp(ireso+1)-ibt
T1 = 1; xl8 = rtempsq(ibt+1)
DO iut = 1, nrtemp
  T2 = iut; xr8 = rtempsq(ibt+iut)
  IF (xut8<xr8) EXIT
  T1 = T2; xl8 = xr8
END DO
IF (T1.EQ.T2) THEN
  wT = 1.
ELSE
  wT = (xr8-xut8)/(xr8-xl8)
END IF
! &&& ================================================== &&& !
!print*, iso, T1, wT

! &&& ========= Isotopic Potential XS Cal. ========= &&& !
siglpiso = 0.;
DO jso = 1, nisoP
  ND = pnum(jso,ipin)
  jdxiso = mapnucl(jso,ipin)
  jreso = IsoG2R(jdxiso)
  IF ((jreso.NE.0).AND.(.NOT.lclad(jdxiso))) THEN
    IF (rifid(jdxiso,ireso).EQ.0) siglpiso = siglpiso+ND*lamsigp(ig,jdxiso)
  ELSE
    siglpiso = siglpiso+ND*lamsigp(ig,jdxiso)
  END IF
!  IF (.NOT.lclad(jdxiso).AND.rifid(jdxiso,ireso).EQ.0) THEN
!    siglpiso = siglpiso+ND*lamsigp(ig,jdxiso)
!  END IF
END DO
!IF (ig.EQ.1 .AND. iso.EQ.2) print*, ig, siglpiso
! &&& =============================================== &&& !

! &&& ============== XSA Calculation =============== &&& !
LB = ptrNlv(ireso)+1; UB = ptrNlv(ireso+1);
ND = ind(idxiso,ipin); ibTL = ptrTempLv(ireso);
xsa = 0.; xsf = 0.; phi = 1.;
IF (lCellAIC(ipin)) THEN
  DO ilv = LB, UB
    lvX = lvabs(ilv,ig)
    xut8 = ND*lvX;
    xut8 = dlog(xut8)
    Mc1 = 1; xl8 = f_maclv1G_log(1,iz)
    DO iut = 1, nmaclv1G
      Mc2 = iut; xr8 = f_maclv1G_log(iut,iz)
      IF (xut8<xr8) EXIT
      Mc1 = Mc2; xl8 = xr8
    END DO
    IF (Mc1.EQ.Mc2) THEN
      wMc = 1.
    ELSE
      wMc = (xr8-xut8)/(xr8-xl8)
    END IF
    sigb = wMc*xseq_1g(Mc1,ipin)+(1.-wMc)*xseq_1g(Mc2,ipin)
    sigb = (sigb+siglpiso)/ND
    wgtX = wT*wgtabs(ibTL+T1,ig)+(1.-wT)*(wgtabs(ibTL+T2,ig))
    wgtX = wgtX*lvX
    lvX = 1./(lvX+sigb) ! phimult
    phi = phi-wgtX*lvX
    sigb = lvX*sigb ! phihom
    xsa = xsa+wgtX*sigb
!    IF (ityp(ireso).EQ.3) THEN
!      lvX = lvfis(ilv,ig); wgtX = wT*wgtfis(ibTL+T1,ig)+(1.-wT)*(wgtfis(ibTL+T2,ig))
!      xsf = xsf+lvX*wgtX*sigb
!    END IF
    ibTL = ibTL+nrtemp
  END DO
ELSE
  adjintmlg = ND/FnAdj(ig,ipin)
  DO ilv = LB, UB
    lvX = lvabs(ilv,ig)
    xut8 = adjintmlg*lvX;
    xut8 = dlog(xut8)
    Mc1 = 1; xl8 = f_maclv_log(1,iz)
    DO iut = 1, nmaclv
      Mc2 = iut; xr8 = f_maclv_log(iut,iz)
      IF (xut8<xr8) EXIT
      Mc1 = Mc2; xl8 = xr8
    END DO
    IF (Mc1.EQ.Mc2) THEN
      wMc = 1.
    ELSE
      wMc = (xr8-xut8)/(xr8-xl8)
    END IF
    sigb = wMc*xseq_mg(Mc1,ig,ipin)+(1.-wMc)*xseq_mg(Mc2,ig,ipin)
    sigb = (sigb+siglpiso)/ND
!    IF (iso .EQ. 2) print*, ilv-LB+1, sigb
    wgtX = wT*wgtabs(ibTL+T1,ig)+(1.-wT)*(wgtabs(ibTL+T2,ig))
    wgtX = wgtX*lvX
!    IF (iso.EQ.2 .AND. ig.EQ.1) print*, ilv-LB+1, wgtX
    lvX = 1./(lvX+sigb) ! phimult
    phi = phi-wgtX*lvX
    sigb = lvX*sigb ! phihom
    xsa = xsa+wgtX*sigb
!    IF (ityp(ireso).EQ.3) THEN
!      lvX = lvfis(ilv,ig); wgtX = wT*wgtfis(ibTL+T1,ig)+(1.-wT)*(wgtfis(ibTL+T2,ig))
!      xsf = xsf+lvX*wgtX*sigb
!    END IF
    ibTL = ibTL+nrtemp
  END DO
END IF
! &&& ============================================== &&& !
!IF (iso.EQ.9 .AND. ig.EQ.4) print*, iso, ig, xsa
xsa = xsa/phi
!xsf = xsf/phi
IF (xsa.LE.0.) RETURN
xsa = dlog(xsa)

! &&& =============== XSA Weights Finding ============== &&& !
nsig0 = ptrSig0(ireso+1)-ptrSig0(ireso)
xut4 = xsa
!   ### For T1
ibTL = ptrTempPot(ireso)+(T1-1)*nsig0
AL2 = nsig0; xr4 = xsalog(ibTL+nsig0,ig)
DO iut = nsig0,1,-1
  AL1 = iut; xl4 = xsalog(ibTL+iut,ig)
  IF (xut4.GE.xl4) EXIT
  AL2 = AL1; xr4 = xl4
END DO
IF (AL1 .EQ. AL2) THEN
  wAL = 1.
ELSE
  wAL = (xr4-xut4)/(xr4-xl4)
END IF
!IF (iso.EQ.9 .AND. ig.EQ.4) print*, xut4, xsalog(ibTL+AL1-1,ig), xsalog(ibTL+AL1,ig)
!   ### For T2
ibTL = ibTL+(T2-T1)*nsig0
AR2 = nsig0; xr4 = xsalog(ibTL+nsig0,ig)
DO iut = nsig0,1,-1
  AR1 = iut; xl4 = xsalog(ibTL+iut,ig)
  IF (xut4.GE.xl4) EXIT
  AR2 = AR1; xr4 = xl4
END DO
IF (AR1 .EQ. AR2) THEN
  wAR = 1.
ELSE
  wAR = (xr4-xut4)/(xr4-xl4)
END IF
!IF (iso.EQ.9 .AND. ig.EQ.4) print*, xsalog(ibTL+AR2,ig), xsalog(ibTL+AR1-1,ig), xsalog(ibTL+AR1,ig)
!IF (iso.EQ.9 .AND. ig.EQ.4) print*, AR1, AR2
! &&& ================================================== &&& !
!IF (iso.EQ.2 .AND. ig .EQ. 1) print*, AL2, AR2, wAL, wAR

! &&& ========== RIF Factor Calculation & Sum =========== &&& !
rifaT = 0.; riffT = 0.; IsFis = ityp(ireso)
DO jso = 1, nisoP
  jdxiso = mapnucl(jso,ipin)
  ibt = rifid(jdxiso,ireso)
  IF (ibt.EQ.0) CYCLE
  xut8 = ind(jdxiso,ipin)
  IF (xut8.LE.0.) CYCLE
  ! ### Weights for ND ratio
  UB = ptrRifRat(ibt+1)
  LB = ptrRifRat(ibt)+1

  xut8 = dlog(xut8/ND)
  R1 = LB; xl8 = ratlog(LB)
  DO iut = LB, UB
    R2 = iut; xr8 = ratlog(iut)
    IF (xut8<xr8) EXIT
    R1 = R2; xl8 = xr8
  END DO
  IF (R1 .EQ. R2) THEN
    wR = 1.
  ELSE
    wR = (xr8-xut8)/(xr8-xl8)
  END IF

  ! ### (R1,T1)
!  ibTL = ptrRifRx(R2)+T2*nsig0+AR2
!  print*, ibTL, size(rifabs,1), (size(rifabs,1).GT.ibTL)
  ibTL = ptrRifRx(R1)+(T1-1)*nsig0
  rifaT(1) = rifaT(1)+(rifabs(ibTL+AL1,ig)*wAL+rifabs(ibTL+AL2,ig)*(1.-wAL))*wR
  IF (IsFis) riffT(1) = riffT(1)+(riffis(ibTL+AL1,ig)*wAL+riffis(ibTL+AL2,ig)*(1.-wAL))*wR
  ! ### (R1,T2)
  ibTL = ibTL+(T2-T1)*nsig0
  rifaT(2) = rifaT(2)+(rifabs(ibTL+AR1,ig)*wAR+rifabs(ibTL+AR2,ig)*(1.-wAR))*wR - 1.
  IF (IsFis) riffT(2) = riffT(2)+(riffis(ibTL+AR1,ig)*wAR+riffis(ibTL+AR2,ig)*(1.-wAR))*wR - 1.
  ! ### (R2,T2)
  ibTL = ibTL+(R2-R1)*nsig0*nrtemp
  rifaT(2) = rifaT(2)+(rifabs(ibTL+AR1,ig)*wAR+rifabs(ibTL+AR2,ig)*(1.-wAR))*(1.-wR)
  IF (IsFis) riffT(2) = riffT(2)+(riffis(ibTL+AR1,ig)*wAR+riffis(ibTL+AR2,ig)*(1.-wAR))*(1.-wR)
  ! ### (R2,T1)
  ibTL = ibTL+(T1-T2)*nsig0
  rifaT(1) = rifaT(1)+(rifabs(ibTL+AL1,ig)*wAL+rifabs(ibTL+AL2,ig)*(1.-wAL))*(1.-wR) - 1.
  IF (IsFis) riffT(1) = riffT(1)+(riffis(ibTL+AL1,ig)*wAL+riffis(ibTL+AL2,ig)*(1.-wAL))*(1.-wR) - 1.
END DO
! &&& =================================================== &&& !

rifa(ireso,ig,ipin) = rifaT(1)*wT+rifaT(2)*(1.-wT)+1.
riff(ireso,ig,ipin) = riffT(1)*wT+riffT(2)*(1.-wT)+1.

!  IF (iso.EQ.3 .AND. ig.EQ.6 .AND. nisoP.EQ.12) THEN
!    print*, riffT(1), riffT(2), riff(ireso,ig,ipin)
!  END IF
!IF (iso.EQ.9 .AND. ig.EQ.4) print*, AL1, AR1
!IF (iso.EQ.3 .AND. ig.EQ.11) print*, rifa(ireso,ig,ipin), AL1, AR1

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuEffMacXS_GenD(ityp,IsoG2R,ptrRTemp,ptrNlv,ptrTempLv, &
                                lamsigp,lvabs,lvfis,wgtabs,wgtfis,maclv_log,rtempsq, &
                                rifa,riff,xseq_mg,mapnucl,ipinD,mapGR2G,mapG2R, &
                                pnum,ind,TempAvg,siglo_nreso,FnAdj,FtAdj,McAIso,McFIso, &
                                nGenD, ntiso, nrg, nmaclv, iz)
USE CUDAFOR
IMPLICIT NONE
INTEGER,DEVICE,DIMENSION(:) :: ityp,IsoG2R,ptrRTemp,ptrNlv,ptrTempLv
REAL(GPU_XS_PRECISION),DEVICE,DIMENSION(:,:) :: lamsigp,lvabs,lvfis,wgtabs,wgtfis,maclv_log
REAL(GPU_XS_PRECISION),DEVICE,DIMENSION(:) :: rtempsq
REAL(GPU_RES_PRECISION),DEVICE,DIMENSION(:,:,:) :: rifa, riff, xseq_mg
INTEGER,DEVICE :: mapnucl(:,:),ipinD(:),mapGR2G(:),mapG2R(:)
REAL(GPU_PNUM_PRECISION),DEVICE :: pnum(:,:),ind(:,:)
REAL(GPU_XS_PRECISION),DEVICE :: TempAvg(:),siglo_nreso(:,:)
REAL(GPU_RES_PRECISION),DEVICE :: McAIso(:,:,:), McFIso(:,:,:),FnAdj(:,:),FtAdj(:,:,:) ! stored at FxrD%fresoAIso(:,:,:),fresoFIso(:,:,:)
INTEGER,VALUE :: nGenD, ntiso, nrg, nmaclv, iz

INTEGER :: ifxrGR, ifxr, ifxrR, ipin, ig, iso, idxiso, ireso
INTEGER :: nrtemp
REAL(8) :: xut8, xl8, xr8
INTEGER :: T1, T2, Mc1, Mc2
REAL(8) :: wT, wMc
INTEGER :: iut, ibt, ilv, ibTL, LB, UB

REAL(8) :: siglpiso, ND, adjintmlg, xsa, xsf, phi, lvX, sigb, wgtX
REAL(4) :: FtIso
LOGICAL :: IsFis

ifxrGR = threadIdx%x+(blockIdx%x-1)*EFF_BLOCK_RDIM
IF (ifxrGR.GT.nGenD) RETURN
iso = threadIdx%y+(blockIdx%y-1)*EFF_BLOCK_IDIM
IF (iso.GT.ntiso) RETURN
ig = blockIdx%z
IF (ig.GT.nrg) RETURN

ifxr = mapGR2G(ifxrGR); idxiso = mapnucl(iso,ifxr)
IF (idxiso.EQ.0) RETURN
ifxrR = mapG2R(ifxr); ireso = IsoG2R(idxiso)
IF (ireso.EQ.0) RETURN
ND = ind(idxiso,ifxr)
IF (ND.LE.0.) RETURN
IsFis = (ityp(ireso).EQ.3)

ipin = ipinD(ifxr)
IF ((ipin-1)/nxy.NE.(iz-1)) RETURN
! &&& ========= AvgTemperature Weights Finding ========= &&& !
xut8 = TempAvg(ipin); xut8 = sqrt(xut8)
ibt = ptrRTemp(ireso); nrtemp = ptrRTemp(ireso+1)-ibt
T1 = 1; xl8 = rtempsq(ibt+1)
DO iut = 1, nrtemp
  T2 = iut; xr8 = rtempsq(ibt+iut)
  IF (xut8<xr8) EXIT
  T1 = T2; xl8 = xr8
END DO
IF (T1.EQ.T2) THEN
  wT = 1.
ELSE
  wT = (xr8-xut8)/(xr8-xl8)
END IF
! &&& ================================================== &&& !

! &&& ========== XSA and XSF Calculation ============ &&& !
siglpiso = siglo_nreso(ig,ifxrR)+ND*lamsigp(ig,idxiso)
adjintmlg = ND/FnAdj(ig,ifxrGR)
!IF (ifxrR.EQ.2 .AND. ig.EQ.1 .AND. iso .EQ. 2) print*, siglpiso, ND, adjintmlg
LB = ptrNlv(ireso)+1; UB = ptrNlv(ireso+1);
ibTL = ptrTempLv(ireso);
xsa = 0.; xsf = 0.; phi = 1.;
DO ilv = LB, UB
  lvX = lvabs(ilv,ig)
  xut8 = adjintmlg*lvX;
  xut8 = log(xut8)
  Mc1 = 1; xl8 = maclv_log(1,iz)
  DO iut = 1, nmaclv
    Mc2 = iut; xr8 = maclv_log(iut,iz)
    IF (xut8<xr8) EXIT
    Mc1 = Mc2; xl8 = xr8
  END DO
  IF (Mc1.EQ.Mc2) THEN
    wMc = 1.
  ELSE
    wMc = (xr8-xut8)/(xr8-xl8)
  END IF
  sigb = wMc*xseq_mg(Mc1,ig,ifxrGR)+(1.-wMc)*xseq_mg(Mc2,ig,ifxrGR)
  FtIso = wMc*FtAdj(Mc1,ig,ifxrGR)+(1.-wMc)*FtAdj(Mc2,ig,ifxrGR)
!  IF (ifxrR.EQ.2 .AND. ig.EQ.1 .AND. iso .EQ. 2) print*, 'A', sigb, FtIso
  xut8 = FtIso*ND ! Might help to conserve precision
  sigb = (sigb+siglpiso)/xut8
  wgtX = wT*wgtabs(ibTL+T1,ig)+(1.-wT)*(wgtabs(ibTL+T2,ig))
!  IF (ifxrR.EQ.2 .AND. ig.EQ.1 .AND. iso .EQ. 2) print*, 'W', wgtX, lvX
  wgtX = wgtX*lvX
!  IF (ifxrR.EQ.2 .AND. ig.EQ.1 .AND. iso .EQ. 2) print*, 'S', sigb, wgtX
  lvX = 1./(lvX+sigb) ! phimult
  phi = phi-wgtX*lvX
  sigb = FtIso*lvX*sigb ! phihom
  xsa = xsa+wgtX*sigb
!  IF (ifxrR.EQ.2 .AND. ig.EQ.1 .AND. iso .EQ. 2) print*, (ilv-LB+1), xsa, phi
  IF (IsFis) THEN
    lvX = lvfis(ilv,ig); wgtX = wT*wgtfis(ibTL+T1,ig)+(1.-wT)*(wgtfis(ibTL+T2,ig))
    xsf = xsf+lvX*wgtX*sigb
  END IF
  ibTL = ibTL+nrtemp
END DO
! &&& =============================================== &&& !

! &&& ============ Isotopic Mac. XS Gen ============= &&& !
ND = pnum(iso,ifxr); xsa = xsa/phi; xsf = xsf/phi;
McAIso(iso,ig,ifxrR) = xsa*ND*rifa(ireso,ig,ipin)
IF (IsFis) McFIso(iso,ig,ifxrR) = xsf*ND*riff(ireso,ig,ipin)
! &&& =================================================== !
END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuEffMacXS_CldT(IsoG2R,ptrRTemp,ptrSig0,ptrNlv,ptrRifRat,ptrTempPot,ptrTempLv,ptrRifRx, &
                                lamsigp,lvabs,wgtabs,maclv_log,rtempsq, ratlog,xsalog,rifabs,xseq_1g, &
                                AccNtR,AccNt,mapnucl,rifid,mapC2G,mapG2R,ipinT,pnum,Temp,siglo_nreso,McAIso, &
                                nCldT, nrg, nmaclv, iz)
USE CUDAFOR
IMPLICIT NONE
INTEGER,DEVICE,DIMENSION(:) :: IsoG2R,ptrRTemp,ptrSig0,ptrNlv,ptrRifRat,ptrTempPot,ptrTempLv,ptrRifRx
REAL(GPU_XS_PRECISION),DEVICE,DIMENSION(:,:) :: lamsigp,lvabs,wgtabs
REAL(GPU_XS_PRECISION),DEVICE,DIMENSION(:) :: rtempsq, ratlog,maclv_log
REAL(GPU_RES_PRECISION),DEVICE,DIMENSION(:,:) :: xsalog,rifabs
REAL(GPU_RES_PRECISION),DEVICE :: xseq_1g(:,:),McAIso(:)
INTEGER,DEVICE :: AccNtR(:),AccNt(:),mapnucl(:),rifid(:,:),mapC2G(:),mapG2R(:),ipinT(:)
REAL(GPU_PNUM_PRECISION),DEVICE :: pnum(:)
REAL(GPU_XS_PRECISION),DEVICE :: Temp(:),siglo_nreso(:,:)
INTEGER,VALUE :: nCldT, nrg, nmaclv, iz

INTEGER :: ifxrC, ifxr, ifxrR, ig, iso, isoRG, miso, liso, idxiso, ireso, jdxiso, jreso
INTEGER :: nisoP, nrtemp, nsig0
REAL(8) :: xut8, xl8, xr8
REAL(4) :: xut4, xl4, xr4
INTEGER :: TP1, TP2, Mc1, Mc2, AL1, AL2, AR1, AR2, R1, R2
REAL(8) :: wTP, wMc, wR
REAL(4) :: wAL, wAR
INTEGER :: iut, ibt, ilv, jso, ibTL, LB, UB

REAL(8) :: siglpiso, ND, xsa, xsf, phi, lvX, sigb, wgtX, lpNoR
REAL(8) :: rifaT(2)

ifxrC = threadIdx%x+EFF_BLOCK_RDIM*(blockIdx%x-1)
IF (ifxrC.GT.nCldT) RETURN
ig = threadIdx%y+EFF_BLOCK_GDIM*(blockIdx%y-1)
IF (ig.GT.nrg) RETURN

ifxr = mapC2G(ifxrC)
IF (ifxr.LT.ifxrbegT(iz)) RETURN
IF (ifxr.GE.ifxrbegT(iz+1)) RETURN
ifxrR = mapG2R(ifxr)

lpNoR = siglo_nreso(ig,ifxrR)
miso = AccNt(ifxr)+1; liso = AccNt(ifxr+1)
isoRG = AccNtR(ifxrR)*nrg+(ig-1)*(liso-miso+1)
DO iso = miso, liso
  isoRG = isoRG+1
  ND = pnum(iso); idxiso = mapnucl(iso)
!  IF (ig.EQ.1) print*, ifxrc, idxiso, (iso-miso+1)
  IF (idxiso.EQ.0) CYCLE
  ireso = IsoG2R(idxiso)
  IF (ireso.EQ.0) CYCLE
  ! &&& ========= AvgTemperature Weights Finding ========= &&& !
  xut8 = Temp(ifxr); xut8 = sqrt(xut8)
  ibt = ptrRTemp(ireso); nrtemp = ptrRTemp(ireso+1)-ibt
  TP1 = 1; xl8 = rtempsq(ibt+1)
  DO iut = 1, nrtemp
    TP2 = iut; xr8 = rtempsq(ibt+iut)
    IF (xut8<xr8) EXIT
    TP1 = TP2; xl8 = xr8
  END DO
  IF (TP1.EQ.TP2) THEN
    wTP = 1.
  ELSE
    wTP = (xr8-xut8)/(xr8-xl8)
  END IF
  ! &&& ================================================== &&& !

  siglpiso = lpNoR+ND*lamsigp(ig,idxiso)
!  IF (ig.EQ.1) print*, lpNoR, siglpiso
  ! &&& ========== XSA and XSF Calculation ============ &&& !
  LB = ptrNlv(ireso)+1; UB = ptrNlv(ireso+1);
  ibTL = ptrTempLv(ireso);
  xsa = 0.; phi = 1.;
  DO ilv = LB, UB
    lvX = lvabs(ilv,ig)
    xut8 = ND*lvX;
    xut8 = log(xut8)
    Mc1 = 1; xl8 = maclv_log(1)
    DO iut = 1, nmaclv
      Mc2 = iut; xr8 = maclv_log(iut)
      IF (xut8<xr8) EXIT
      Mc1 = Mc2; xl8 = xr8
    END DO
    IF (Mc1.EQ.Mc2) THEN
      wMc = 1.
    ELSE
      wMc = (xr8-xut8)/(xr8-xl8)
    END IF
    sigb = wMc*xseq_1g(Mc1,ifxrC)+(1.-wMc)*xseq_1g(Mc2,ifxrC)
    sigb = (sigb+siglpiso)/ND
!    IF (ig.EQ.1) print*, (ilv-LB+1), sigb
    wgtX = wTP*wgtabs(ibTL+TP1,ig)+(1.-wTP)*(wgtabs(ibTL+TP2,ig))
    wgtX = wgtX*lvX
!    IF (ig.EQ.1) print*, (ilv-LB+1), wTP, TP1
    lvX = 1./(lvX+sigb) ! phimult
    phi = phi-wgtX*lvX
    sigb = lvX*sigb ! phihom
    xsa = xsa+wgtX*sigb
    ibTL = ibTL+nrtemp
  END DO
  ! &&& =============================================== &&& !
  xsa = xsa/phi
!  print*, ig, xsa
  IF (xsa.LE.0.) CYCLE
  ! &&& =============== XSA Weights Finding ============== &&& !
  nsig0 = ptrSig0(ireso+1)-ptrSig0(ireso)
  xut4 = log(xsa)
  !   ### For T1
  ibTL = ptrTempPot(ireso)+(TP1-1)*nsig0
  AL2 = nsig0; xr4 = xsalog(ibTL+nsig0,ig)
  DO iut = nsig0,1,-1
    AL1 = iut; xl4 = xsalog(ibTL+iut,ig)
    IF (xut4.GE.xl4) EXIT
    AL2 = AL1; xr4 = xl4
  END DO
  IF (AL1 .EQ. AL2) THEN
    wAL = 1.
  ELSE
    wAL = (xr4-xut4)/(xr4-xl4)
  END IF
  !   ### For T2
  ibTL = ibTL+(TP2-TP1)*nsig0
  AR2 = nsig0; xr4 = xsalog(ibTL+nsig0,ig)
  DO iut = nsig0,1,-1
    AR1 = iut; xl4 = xsalog(ibTL+iut,ig)
    IF (xut4.GE.xl4) EXIT
    AR2 = AR1; xr4 = xl4
  END DO
  IF (AR1 .EQ. AR2) THEN
    wAR = 1.
  ELSE
    wAR = (xr4-xut4)/(xr4-xl4)
  END IF
  ! &&& ================================================== &&& !

  ! &&& ========== RIF Factor Calculation & Sum =========== &&& !
  rifaT = 0.;
  DO jso = miso, liso
    jdxiso = mapnucl(jso)
    ibt = rifid(jdxiso,ireso)
    IF (ibt.EQ.0) CYCLE
    xut8 = pnum(jso)
    IF (xut8.LE.0.) CYCLE
    ! ### Weights for ND ratio
    UB = ptrRifRat(ibt+1)
    LB = ptrRifRat(ibt)+1

    xut8 = log(xut8/ND)
    R1 = LB; xl8 = ratlog(LB)
    DO iut = LB, UB
      R2 = iut; xr8 = ratlog(iut)
      IF (xut8<xr8) EXIT
      R1 = R2; xl8 = xr8
    END DO
    IF (R1 .EQ. R2) THEN
      wR = 1.
    ELSE
      wR = (xr8-xut8)/(xr8-xl8)
    END IF

    ! ### (R1,T1)
    ibTL = ptrRifRx(R1)+(TP1-1)*nsig0
    rifaT(1) = rifaT(1)+(rifabs(ibTL+AL1,ig)*wAL+rifabs(ibTL+AL2,ig)*(1.-wAL))*wR
    ! ### (R1,T2)
    ibTL = ibTL+(TP2-TP1)*nsig0
    rifaT(2) = rifaT(2)+(rifabs(ibTL+AR1,ig)*wAR+rifabs(ibTL+AR2,ig)*(1.-wAR))*wR - 1.
    ! ### (R2,T2)
    ibTL = ibTL+(R2-R1)*nsig0*nrtemp
    rifaT(2) = rifaT(2)+(rifabs(ibTL+AR1,ig)*wAR+rifabs(ibTL+AR2,ig)*(1.-wAR))*(1.-wR)
    ! ### (R2,T1)
    ibTL = ibTL+(TP1-TP2)*nsig0
    rifaT(1) = rifaT(1)+(rifabs(ibTL+AL1,ig)*wAL+rifabs(ibTL+AL2,ig)*(1.-wAL))*(1.-wR) - 1.
  END DO
  ! &&& =================================================== &&& !
  McAIso(isoRG) = xsa*ND*(rifaT(1)*wTP+rifaT(2)*(1.-wTP)+1.)
END DO
END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuEffMacXS_AICT(IsoG2R,ptrRTemp,ptrNlv,ptrTempLv, &
                                lamsigp,lvabs,wgtabs,maclv_log,rtempsq,rifa,xseq_1g, &
                                AccNtR,AccNt,mapnucl,mapA2G,mapG2R,ipinT, &
                                pnum,TempAvg,siglo_nreso,McAIso, &
                                nAICT, nrg, nmaclv, iz)
USE CUDAFOR
IMPLICIT NONE
INTEGER,DEVICE,DIMENSION(:) :: IsoG2R,ptrRTemp,ptrNlv,ptrTempLv
REAL(GPU_XS_PRECISION),DEVICE,DIMENSION(:,:) :: lamsigp,lvabs,wgtabs,maclv_log
REAL(GPU_XS_PRECISION),DEVICE,DIMENSION(:) :: rtempsq
REAL(GPU_RES_PRECISION),DEVICE:: rifa(:,:,:), xseq_1g(:,:)
INTEGER,DEVICE :: AccNtR(:),AccNt(:),mapnucl(:),mapA2G(:),mapG2R(:),ipinT(:)
REAL(GPU_PNUM_PRECISION),DEVICE :: pnum(:)
REAL(GPU_XS_PRECISION),DEVICE :: TempAvg(:), siglo_nreso(:,:)
REAL(GPU_RES_PRECISION),DEVICE :: McAIso(:) ! stored at FxrTri%fresoAIso(:)
INTEGER,VALUE :: nAICT, nrg, nmaclv, iz

INTEGER :: ifxrA, ifxr, ifxrR, ipin, ig, iso, isoRG, miso, liso, idxiso, ireso
INTEGER :: nrtemp
REAL(8) :: xut8, xl8, xr8
INTEGER :: T1, T2, Mc1, Mc2
REAL(8) :: wT, wMc
INTEGER :: iut, ibt, ilv, ibTL, LB, UB

REAL(8) :: lpNoR, siglpiso, ND, xsa, phi, lvX, sigb, wgtX

ifxrA = threadIdx%x+EFF_BLOCK_RDIM*(blockIdx%x-1)
IF (ifxrA.GT.nAICT) RETURN
ig = threadIdx%y+EFF_BLOCK_GDIM*(blockIdx%y-1)
IF (ig.GT.nrg) RETURN

ifxr = mapA2G(ifxrA)
ifxrR = mapG2R(ifxr)

ipin = ipinT(ifxr)
IF ((ipin-1)/nxy.NE.(iz-1)) RETURN

lpNoR = siglo_nreso(ig,ifxrR)
miso = AccNt(ifxr)+1; liso = AccNt(ifxr+1)
isoRG = AccNtR(ifxrR)*nrg+(ig-1)*(liso-miso+1)
DO iso = miso, liso
  isoRG = isoRG+1
  ND = pnum(iso); idxiso = mapnucl(iso)
  IF (idxiso.EQ.0) CYCLE
  ireso = IsoG2R(idxiso)
  IF (ireso.EQ.0) CYCLE
  ! &&& ========= AvgTemperature Weights Finding ========= &&& !
  xut8 = TempAvg(ipin); xut8 = sqrt(xut8)
  ibt = ptrRTemp(ireso); nrtemp = ptrRTemp(ireso+1)-ibt
  T1 = 1; xl8 = rtempsq(ibt+1)
  DO iut = 1, nrtemp
    T2 = iut; xr8 = rtempsq(ibt+iut)
    IF (xut8<xr8) EXIT
    T1 = T2; xl8 = xr8
  END DO
  IF (T1.EQ.T2) THEN
    wT = 1.
  ELSE
    wT = (xr8-xut8)/(xr8-xl8)
  END IF
  ! &&& ================================================== &&& !
  siglpiso = lpNoR+ND*lamsigp(ig,idxiso)
  ! &&& ========== XSA and XSF Calculation ============ &&& !
  LB = ptrNlv(ireso)+1; UB = ptrNlv(ireso+1);
  ibTL = ptrTempLv(ireso);
  xsa = 0.; phi = 1.;
  DO ilv = LB, UB
    lvX = lvabs(ilv,ig)
    xut8 = ND*lvX;
    xut8 = log(xut8)
    Mc1 = 1; xl8 = maclv_log(1,iz)
    DO iut = 1, nmaclv
      Mc2 = iut; xr8 = maclv_log(iut,iz)
      IF (xut8<xr8) EXIT
      Mc1 = Mc2; xl8 = xr8
    END DO
    IF (Mc1.EQ.Mc2) THEN
      wMc = 1.
    ELSE
      wMc = (xr8-xut8)/(xr8-xl8)
    END IF
    sigb = wMc*xseq_1g(Mc1,ifxrA)+(1.-wMc)*xseq_1g(Mc2,ifxrA)
    sigb = (sigb+siglpiso)/ND
    wgtX = wT*wgtabs(ibTL+T1,ig)+(1.-wT)*(wgtabs(ibTL+T2,ig))
    wgtX = wgtX*lvX
    lvX = 1./(lvX+sigb) ! phimult
    phi = phi-wgtX*lvX
    sigb = lvX*sigb ! phihom
    xsa = xsa+wgtX*sigb
    ibTL = ibTL+nrtemp
  END DO
  ! &&& =============================================== &&& !
  xsa = xsa/phi
  McAIso(isoRG) = xsa*ND*rifa(ireso,ig,ipin)
END DO
END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuEffFreso_D(fresoAIso,fresoFIso,fresoA,fresoF,fresonF,&
                                ityp,siga,sigf,signf,wtTmp,niso,mapnucl,itTmp,mapR2G,IsoG2R,pnum,&
                                bg,nrg,nResD,iz)
USE CUDAFOR
IMPLICIT NONE
REAL(GPU_RES_PRECISION),DEVICE :: fresoAIso(:,:,:),fresoFIso(:,:,:),fresoA(:,:),fresoF(:,:),fresonF(:,:)
REAL(GPU_XS_PRECISION),DEVICE,DIMENSION(:,:) :: siga,sigf,signf,wtTmp
INTEGER,DEVICE :: niso(:),mapnucl(:,:),itTmp(:,:),mapR2G(:),IsoG2R(:),ityp(:)
REAL(GPU_PNUM_PRECISION),DEVICE :: pnum(:,:)
INTEGER, VALUE :: bg, nrg, nResD, iz

INTEGER :: ig, igR, ifxrR, ifxr
INTEGER :: iso, miso, ireso
INTEGER :: jt, jt2
REAL(8) :: sigbuf, sigbuf2, wtN

REAL(8) :: ND, micA, micnF, micF, effA, effF, effnF, invF
LOGICAL :: lskipit2, lreso, lfis

!xtid = threadIdx%x;
igR = threadIdx%x

!ytid = threadIdx%y
ifxrR  = threadIdx%y+EFF_BLOCK_RDIM*(blockIdx%y-1)
if (ifxrR .GT. nResD) RETURN
!IF (ig .GT. nrg) RETURN

ifxr = mapR2G(ifxrR)
IF (ifxr.LT.ifxrbegD(iz)) RETURN
IF (ifxr.GE.ifxrbegD(iz+1)) RETURN
!print*, 'I', ifxr
ig = igR+bg

micA = 0.; micnF = 0.; micF = 0.
effA = 0.; effnF = 0.; effF = 0.

miso = niso(ifxr);
!print*, 'n', miso
!jt = 1;
DO iso = 1, miso
  ireso = IsoG2R(mapnucl(iso,ifxr))
  lreso = (ireso.NE.0); lfis = (ityp(ireso).EQ.3)
  ND = pnum(iso,ifxr)
  jt = itTmp(ifxr,iso); wtN = wtTmp(ifxr,iso)
!  print*, 'wt', jt, wtN
  lskipit2 = (wtN .GT. 0.99999)
  jt2 = jt+1
  IF (lskipit2) jt2 = jt

  wtN = ND*wtN;
  IF (lreso) THEN
    ! Absorption
    sigbuf = siga(ig, jt); sigbuf2 = siga(ig,jt2)
!    print*, 's', sigbuf, sigbuf2
    sigbuf2 = wtN*sigbuf+(ND-wtN)*sigbuf2
    micA = micA+sigbuf2

    sigbuf = fresoAIso(iso,igR,ifxrR)
    effA = effA+sigbuf

    fresoAIso(iso,igR,ifxrR) = sigbuf/sigbuf2
!    IF (ifxrR.EQ.2 .AND. igR.EQ.1) print*, iso, sigbuf, sigbuf2

    ! Fission
    sigbuf = sigf(ig, jt); sigbuf2 = sigf(ig,jt2)
    sigbuf2 = wtN*sigbuf+(ND-wtN)*sigbuf2
    IF (sigbuf2 .lE. 1.E-15) CYCLE
    micF = micF+sigbuf2
    invF = 1./sigbuf2

    ! Nu-Fission
    sigbuf = signf(ig, jt); sigbuf2 = signf(ig,jt2)
    sigbuf2 = wtN*sigbuf+(ND-wtN)*sigbuf2
    micnF = micnF+sigbuf2
    sigbuf2 = invF*sigbuf2

    ! Eff Fis & nuF
    IF (lfis) THEN
      sigbuf = fresoFIso(iso,igR,ifxrR)
      effF = effF+sigbuf
!      print*, iso, sigbuf, 1./invF
      effnF = effnF+sigbuf*sigbuf2

      fresoFIso(iso,igR,ifxrR) = sigbuf*invF
    ELSE
      invF = 1./invF
      effF = effF+invF
      effnF = effnF+sigbuf2*invF
    END IF
  ELSE
    ! Absorption
    sigbuf = siga(ig, jt); sigbuf2 = siga(ig,jt2)
    sigbuf2 = wtN*sigbuf+(ND-wtN)*sigbuf2
    micA = micA+sigbuf2
    effA = effA+sigbuf2
!    IF (ifxrR.EQ.2 .AND. igR.EQ.1) print*, iso, sigbuf2, sigbuf2

    ! Fission
    sigbuf = sigf(ig, jt); sigbuf2 = sigf(ig,jt2)
    sigbuf2 = wtN*sigbuf+(ND-wtN)*sigbuf2
    micF = micF+sigbuf2
    effF = effF+sigbuf2

    ! Nu-fission
    sigbuf = signf(ig, jt); sigbuf2 = signf(ig,jt2)
    sigbuf2 = wtN*sigbuf+(ND-wtN)*sigbuf2
    micnF = micnF+sigbuf2
    effnF = effnF+sigbuf2
  END IF
END DO

fresoA(igR,ifxrR) = effA/micA
fresoF(igR,ifxrR) = effF/micF
fresonF(igR,ifxrR) = effnF/micnF
!print*, ifxrR, iso*100+ig, effnF, micnF
!IF (ifxrR.EQ.2 .AND. igR .EQ. 1)print*, igR+ifxrR*100, effA, micA

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuEffFreso_T(fresoAIso,fresoA,siga,wtTmp,&
                                AccNtR,AccNt,niso,mapnucl,itTmp,mapR2G,IsoG2R,pnum,&
                                bg,nrg,nResTri,iz)
USE CUDAFOR
IMPLICIT NONE
REAL(GPU_RES_PRECISION),DEVICE :: fresoAIso(:),fresoA(:,:)
REAL(GPU_XS_PRECISION),DEVICE,DIMENSION(:,:) :: siga,wtTmp
INTEGER,DEVICE :: AccNtR(:),AccNt(:),niso(:),mapnucl(:),itTmp(:,:),mapR2G(:),IsoG2R(:)
REAL(GPU_PNUM_PRECISION),DEVICE :: pnum(:)
INTEGER,VALUE :: bg, nrg, nResTri,iz

INTEGER :: ig, igR, ifxrR, ifxr
INTEGER :: iso, isoR, jso, miso, liso
INTEGER :: jt, jt2
REAL(8) :: sigbuf, sigbuf2, wtN

REAL(8) :: ND, micA, effA
LOGICAL :: lskipit2, lreso

!xtid = threadIdx%x;
igR = threadIdx%x

!ytid = threadIdx%y
ifxrR  = threadIdx%y+EFF_BLOCK_RDIM*(blockIdx%y-1)
if (ifxrR .GT. nResTri) RETURN
!IF (ig .GT. nrg) RETURN

ifxr = mapR2G(ifxrR)
IF (ifxr.LT.ifxrbegT(iz)) RETURN
IF (ifxr.GE.ifxrbegT(iz+1)) RETURN
ig = igR+bg

micA = 0.; effA = 0.

miso = AccNt(ifxr)
liso = AccNt(ifxr+1)
isoR = AccNtR(ifxrR)*nrg+(igR-1)*(liso-miso)
liso = niso(ifxr)+miso;
miso = miso+1
DO iso = miso, liso
  isoR = isoR+1
  ND = pnum(iso)
  jso = iso-miso+1
  jt = itTmp(ifxr,jso); wtN = wtTmp(ifxr,jso)
  lskipit2 = (wtN .GT. 0.99999)
  jt2 = jt+1
  IF (lskipit2) jt2 = jt
  lreso = (IsoG2R(mapnucl(iso)).NE.0)

  wtN = ND*wtN;
  sigbuf = siga(ig, jt); sigbuf2 = siga(ig,jt2)
  sigbuf2 = wtN*sigbuf+(ND-wtN)*sigbuf2
  micA = micA+sigbuf2
  IF (lreso) THEN
    sigbuf = fresoAIso(isoR)
    effA = effA+sigbuf
    fresoAIso(isoR) = sigbuf/sigbuf2
  ELSE
    effA = effA+sigbuf2
  END IF
END DO
fresoA(igR,ifxrR) = effA/micA
END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuEffIntegrated_GenD(ityp,IsoG2R,ptrRTemp,ptrNlv,ptrTempLv, &
                                lamsigp,lvabs,lvfis,wgtabs,wgtfis,maclv_log,rtempsq, &
                                rifa,riff,xseq_mg,mapnucl,ipinD,mapGR2G,mapG2R, &
                                pnum,ind,TempAvg,siglo_nreso,FnAdj,FtAdj, &
                                fresoAIso,fresoFIso,fresoA,fresoF,fresonF, &
                                siga,sigf,signf,wtTmp,niso,itTmp,&
                                nGenD,bg,gb,ge,ifxrb,ifxre,nmaclv,iz)
USE CUDAFOR
IMPLICIT NONE
INTEGER,DEVICE,DIMENSION(:) :: ityp,IsoG2R,ptrRTemp,ptrNlv,ptrTempLv
REAL(GPU_XS_PRECISION),DEVICE,DIMENSION(:,:) :: lamsigp,lvabs,lvfis,wgtabs,wgtfis,maclv_log
REAL(GPU_XS_PRECISION),DEVICE,DIMENSION(:) :: rtempsq
REAL(GPU_RES_PRECISION),DEVICE,DIMENSION(:,:,:) :: rifa, riff, xseq_mg
INTEGER,DEVICE :: mapnucl(:,:),ipinD(:),mapGR2G(:),mapG2R(:)
REAL(GPU_PNUM_PRECISION),DEVICE :: pnum(:,:),ind(:,:)
REAL(GPU_XS_PRECISION),DEVICE :: TempAvg(:), siglo_nreso(:,:)
REAL(GPU_RES_PRECISION),DEVICE :: FnAdj(:,:),FtAdj(:,:,:) ! stored at FxrD%fresoAIso(:,:,:),fresoFIso(:,:,:)
INTEGER,VALUE :: nGenD, nmaclv, iz
REAL(GPU_RES_PRECISION),DEVICE :: fresoAIso(:,:,:),fresoFIso(:,:,:),fresoA(:,:),fresoF(:,:),fresonF(:,:)
REAL(GPU_XS_PRECISION),DEVICE,DIMENSION(:,:) :: siga,sigf,signf,wtTmp
INTEGER,DEVICE :: niso(:),itTmp(:,:)
INTEGER, VALUE :: bg,gb,ge,ifxrb,ifxre

INTEGER :: ifxrGR, ifxr, ifxrR, ipin, iso, idxiso, ireso
INTEGER :: nrtemp
REAL(8) :: xut8, xl8, xr8
INTEGER :: T1, T2, Mc1, Mc2
REAL(8) :: wT, wMc
INTEGER :: iut, ibt, ilv, ibTL, LB, UB

REAL(8) :: siglpiso, ND, adjintmlg, xsa, xsf, phi, lvX, sigb, wgtX, FtIso
LOGICAL :: IsFis

INTEGER :: ig, igR, jgR
INTEGER :: miso
INTEGER :: jt, jt2
REAL(8) :: sigbuf, sigbuf2, wtN

REAL(8) :: micA, micnF, micF, effA, effF, effnF, invF
LOGICAL :: lskipit2, lreso, lfis

ifxrGR = threadIdx%y+(blockIdx%y-1)*EFF_BLOCK_RDIM
IF (ifxrGR.GT.nGenD) RETURN
igR = threadIdx%x
ig = igR+bg
jgR = ig-gb+1
IF (ig.GT.ge) RETURN
IF (ig.LT.gb) RETURN
ifxr = mapGR2G(ifxrGR)
IF (ifxr.LT.ifxrbegD(iz)) RETURN
IF (ifxr.GE.ifxrbegD(iz+1)) RETURN
ifxrR = mapG2R(ifxr)
IF (ifxrR.GT.ifxre) RETURN
IF (ifxrR.LT.ifxrb) RETURN

micA = 0.; micnF = 0.; micF = 0.
effA = 0.; effnF = 0.; effF = 0.

miso = niso(ifxr); ipin = ipinD(ifxr)
DO iso = 1, miso
  idxiso = mapnucl(iso,ifxr)
  ireso = IsoG2R(idxiso)
  ND = ind(idxiso,ifxr)
!  IF (idxiso.EQ.306) print*, ireso
  IsFis = (ityp(ireso).EQ.3)
  lreso = (ireso.NE.0)
  IF (lreso) THEN
    ! &&& ========= AvgTemperature Weights Finding ========= &&& !
    xut8 = TempAvg(ipin); xut8 = sqrt(xut8)
    ibt = ptrRTemp(ireso); nrtemp = ptrRTemp(ireso+1)-ibt
    T1 = 1; xl8 = rtempsq(ibt+1)
    DO iut = 1, nrtemp
      T2 = iut; xr8 = rtempsq(ibt+iut)
      IF (xut8<xr8) EXIT
      T1 = T2; xl8 = xr8
    END DO
    IF (T1.EQ.T2) THEN
      wT = 1.
    ELSE
      wT = (xr8-xut8)/(xr8-xl8)
    END IF
    ! &&& ================================================== &&& !

    ! &&& ========== XSA and XSF Calculation ============ &&& !
    siglpiso = siglo_nreso(igR,ifxrR)+ND*lamsigp(igR,idxiso)
    adjintmlg = ND/FnAdj(igR,ifxrGR)
    LB = ptrNlv(ireso)+1; UB = ptrNlv(ireso+1);
    ibTL = ptrTempLv(ireso);
    xsa = 0.; xsf = 0.; phi = 1.;
    DO ilv = LB, UB
      lvX = lvabs(ilv,igR)
      xut8 = adjintmlg*lvX;
      xut8 = log(xut8)
      Mc1 = 1; xl8 = maclv_log(1,iz)
      DO iut = 1, nmaclv
        Mc2 = iut; xr8 = maclv_log(iut,iz)
        IF (xut8<xr8) EXIT
        Mc1 = Mc2; xl8 = xr8
      END DO
      IF (Mc1.EQ.Mc2) THEN
        wMc = 1.
      ELSE
        wMc = (xr8-xut8)/(xr8-xl8)
      END IF
      sigb = wMc*xseq_mg(Mc1,igR,ifxrGR)+(1.-wMc)*xseq_mg(Mc2,igR,ifxrGR)
      FtIso = wMc*FtAdj(Mc1,igR,ifxrGR)+(1.-wMc)*FtAdj(Mc2,igR,ifxrGR)
      xut8 = FtIso*ND ! Might help to conserve precision
      sigb = (sigb+siglpiso)/xut8
      wgtX = wT*wgtabs(ibTL+T1,igR)+(1.-wT)*(wgtabs(ibTL+T2,igR))
      wgtX = wgtX*lvX
      lvX = 1./(lvX+sigb) ! phimult
      phi = phi-wgtX*lvX
      sigb = FtIso*lvX*sigb ! phihom
      xsa = xsa+wgtX*sigb
      IF (IsFis) THEN
        lvX = lvfis(ilv,igR); wgtX = wT*wgtfis(ibTL+T1,igR)+(1.-wT)*(wgtfis(ibTL+T2,igR))
        xsf = xsf+lvX*wgtX*sigb
      END IF
      ibTL = ibTL+nrtemp
    END DO
    ! &&& =============================================== &&& !

!    if (iso.EQ.3 .AND. ig.EQ.15 .AND. ifxrGR.EQ.976) then
!      print*, xsf, riff(ireso,igR,ipin)
!    end if

    ! &&& ============ Isotopic Mac. XS Gen ============= &&& !
    ND = pnum(iso,ifxr); xsa = xsa/phi; xsf = xsf/phi;
    IF (xsa.LE.0.) lreso = .FALSE.
    xsa = xsa*ND*rifa(ireso,igR,ipin)
    IF (IsFis) xsf = xsf*ND*riff(ireso,igR,ipin)
    ! &&& =================================================== !
  END IF

  ND = pnum(iso,ifxr)
  jt = itTmp(ifxr,iso); wtN = wtTmp(ifxr,iso)
  lskipit2 = (wtN .GT. 0.99999)
  jt2 = jt+1
  IF (lskipit2) jt2 = jt

  wtN = ND*wtN;
  IF (lreso) THEN
    ! Absorption
    sigbuf = siga(ig, jt); sigbuf2 = siga(ig,jt2)
    sigbuf2 = wtN*sigbuf+(ND-wtN)*sigbuf2
    micA = micA+sigbuf2

    effA = effA+xsa

    fresoAIso(iso,jgR,ifxrR-ifxrb+1) = xsa/sigbuf2

    ! Fission
    sigbuf = sigf(ig, jt); sigbuf2 = sigf(ig,jt2)
    sigbuf2 = wtN*sigbuf+(ND-wtN)*sigbuf2
    IF (sigbuf2 .lE. 1.E-15) CYCLE
    micF = micF+sigbuf2
    invF = 1./sigbuf2

    ! Nu-Fission
    sigbuf = signf(ig, jt); sigbuf2 = signf(ig,jt2)
    sigbuf2 = wtN*sigbuf+(ND-wtN)*sigbuf2
    micnF = micnF+sigbuf2
    sigbuf2 = invF*sigbuf2

    ! Eff Fis & nuF
    IF (IsFis) THEN
      effF = effF+xsf
      effnF = effnF+xsf*sigbuf2

      IF (invF.LT.1.e8) fresoFIso(iso,jgR,ifxrR-ifxrb+1) = xsf*invF
!    if (iso.EQ.3 .AND. ig.EQ.15 .AND. ifxrGR.EQ.976) then
!      print*, xsf, 1./invF
!    end if
    ELSE
      invF = 1./invF
      effF = effF+invF
      effnF = effnF+sigbuf2*invF
    END IF
  ELSE
    ! Absorption
    sigbuf = siga(ig, jt); sigbuf2 = siga(ig,jt2)
    sigbuf2 = wtN*sigbuf+(ND-wtN)*sigbuf2
    micA = micA+sigbuf2
    effA = effA+sigbuf2

    ! Fission
    sigbuf = sigf(ig, jt); sigbuf2 = sigf(ig,jt2)
    sigbuf2 = wtN*sigbuf+(ND-wtN)*sigbuf2
    micF = micF+sigbuf2
    effF = effF+sigbuf2

    ! Nu-fission
    sigbuf = signf(ig, jt); sigbuf2 = signf(ig,jt2)
    sigbuf2 = wtN*sigbuf+(ND-wtN)*sigbuf2
    micnF = micnF+sigbuf2
    effnF = effnF+sigbuf2

    fresoAIso(iso,jgR,ifxrR-ifxrb+1) = 1.
    fresoFIso(iso,jgR,ifxrR-ifxrb+1) = 1.
  END IF
END DO

fresoA(jgR,ifxrR-ifxrb+1) = effA/micA
fresoF(jgR,ifxrR-ifxrb+1) = effF/micF
fresonF(jgR,ifxrR-ifxrb+1) = effnF/micnF


END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuEffIntegrated_CldT(IsoG2R,ptrRTemp,ptrSig0,ptrNlv,ptrRifRat,ptrTempPot,ptrTempLv,ptrRifRx, &
                                lamsigp,lvabs,wgtabs,maclv_log,rtempsq, ratlog,xsalog,rifabs,xseq_1g, &
                                AccNtR,AccNt,mapnucl,rifid,mapC2G,mapG2R,ipinT,pnum,Temp,siglo_nreso, &
                                fresoAIso,fresoA,siga,wtTmp,niso,itTmp,&
                                nCldT,bg,nrg,nmaclv,iz)
USE CUDAFOR
IMPLICIT NONE
INTEGER,DEVICE,DIMENSION(:) :: IsoG2R,ptrRTemp,ptrSig0,ptrNlv,ptrRifRat,ptrTempPot,ptrTempLv,ptrRifRx
REAL(GPU_XS_PRECISION),DEVICE,DIMENSION(:,:) :: lamsigp,lvabs,wgtabs
REAL(GPU_XS_PRECISION),DEVICE,DIMENSION(:) :: rtempsq, ratlog,maclv_log
REAL(GPU_RES_PRECISION),DEVICE,DIMENSION(:,:) :: xsalog,rifabs
REAL(GPU_RES_PRECISION),DEVICE :: xseq_1g(:,:)
INTEGER,DEVICE :: AccNtR(:),AccNt(:),mapnucl(:),rifid(:,:),mapC2G(:),mapG2R(:),ipinT(:)
REAL(GPU_XS_PRECISION),DEVICE :: Temp(:),siglo_nreso(:,:)
REAL(GPU_PNUM_PRECISION),DEVICE :: pnum(:)
INTEGER,VALUE :: nCldT, nrg, nmaclv, iz
REAL(GPU_RES_PRECISION),DEVICE :: fresoAIso(:),fresoA(:,:)
REAL(GPU_XS_PRECISION),DEVICE,DIMENSION(:,:) :: siga,wtTmp
INTEGER,DEVICE :: itTmp(:,:), niso(:)
INTEGER,VALUE :: bg

INTEGER :: ifxrC, ifxr, ifxrR, igR, ig, iso, isoRG, miso, liso, idxiso, ireso, jdxiso, jreso
INTEGER :: nisoP, nrtemp, nsig0
REAL(8) :: xut8, xl8, xr8
REAL(4) :: xut4, xl4, xr4
INTEGER :: TP1, TP2, Mc1, Mc2, AL1, AL2, AR1, AR2, R1, R2
REAL(8) :: wTP, wMc, wR, wAL, wAR
INTEGER :: iut, ibt, ilv, jso, ibTL, LB, UB

REAL(8) :: siglpiso, ND, xsa, xsf, phi, lvX, sigb, wgtX, lpNoR
REAL(8) :: rifaT(2)

INTEGER :: jt, jt2
REAL(8) :: sigbuf, sigbuf2, wtN

REAL(8) :: micA, effA
LOGICAL :: lskipit2, lreso

ifxrC = threadIdx%y+EFF_BLOCK_RDIM*(blockIdx%y-1)
IF (ifxrC.GT.nCldT) RETURN
igR = threadIdx%x
IF (igR.GT.nrg) RETURN
ig = igR+bg
ifxr = mapC2G(ifxrC)
IF (ifxr.LT.ifxrbegT(iz)) RETURN
IF (ifxr.GE.ifxrbegT(iz+1)) RETURN
ifxrR = mapG2R(ifxr)

lpNoR = siglo_nreso(igR,ifxrR)
miso = AccNt(ifxr)+1; liso = AccNt(ifxr+1)
isoRG = AccNtR(ifxrR)*nrg+(igR-1)*(liso-miso+1)

micA = 0.; effA = 0.

DO iso = miso, liso
  isoRG = isoRG+1
  ND = pnum(iso); idxiso = mapnucl(iso)
  IF (idxiso.EQ.0) CYCLE
  ireso = IsoG2R(idxiso)
  lreso = (ireso.GT.0)
  IF (lreso) THEN
    ! &&& ========= AvgTemperature Weights Finding ========= &&& !
    xut8 = Temp(ifxr); xut8 = sqrt(xut8)
    ibt = ptrRTemp(ireso); nrtemp = ptrRTemp(ireso+1)-ibt
    TP1 = 1; xl8 = rtempsq(ibt+1)
    DO iut = 1, nrtemp
      TP2 = iut; xr8 = rtempsq(ibt+iut)
      IF (xut8<xr8) EXIT
      TP1 = TP2; xl8 = xr8
    END DO
    IF (TP1.EQ.TP2) THEN
      wTP = 1.
    ELSE
      wTP = (xr8-xut8)/(xr8-xl8)
    END IF
    ! &&& ================================================== &&& !

    siglpiso = lpNoR+ND*lamsigp(igR,idxiso)
    ! &&& ========== XSA and XSF Calculation ============ &&& !
    LB = ptrNlv(ireso)+1; UB = ptrNlv(ireso+1);
    ibTL = ptrTempLv(ireso);
    xsa = 0.; phi = 1.;
    DO ilv = LB, UB
      lvX = lvabs(ilv,igR)
      xut8 = ND*lvX;
      xut8 = log(xut8)
      Mc1 = 1; xl8 = maclv_log(1)
      DO iut = 1, nmaclv
        Mc2 = iut; xr8 = maclv_log(iut)
        IF (xut8<xr8) EXIT
        Mc1 = Mc2; xl8 = xr8
      END DO
      IF (Mc1.EQ.Mc2) THEN
        wMc = 1.
      ELSE
        wMc = (xr8-xut8)/(xr8-xl8)
      END IF
      sigb = wMc*xseq_1g(Mc1,ifxrC)+(1.-wMc)*xseq_1g(Mc2,ifxrC)
      sigb = (sigb+siglpiso)/ND
      wgtX = wTP*wgtabs(ibTL+TP1,igR)+(1.-wTP)*(wgtabs(ibTL+TP2,igR))
      wgtX = wgtX*lvX
      lvX = 1./(lvX+sigb) ! phimult
      phi = phi-wgtX*lvX
      sigb = lvX*sigb ! phihom
      xsa = xsa+wgtX*sigb
      ibTL = ibTL+nrtemp
    END DO
    ! &&& =============================================== &&& !
    xsa = xsa/phi
    IF (xsa.LE.0.) lreso = .FALSE.
  END IF
  IF (lreso) THEN
    ! &&& =============== XSA Weights Finding ============== &&& !
    nsig0 = ptrSig0(ireso+1)-ptrSig0(ireso)
    xut4 = log(xsa)
    !   ### For T1
    ibTL = ptrTempPot(ireso)+(TP1-1)*nsig0
    AL2 = nsig0; xr4 = xsalog(ibTL+nsig0,igR)
    DO iut = nsig0,1,-1
      AL1 = iut; xl4 = xsalog(ibTL+iut,igR)
      IF (xut4.GE.xl4) EXIT
      AL2 = AL1; xr4 = xl4
    END DO
    IF (AL1 .EQ. AL2) THEN
      wAL = 1.
    ELSE
      wAL = (xr4-xut4)/(xr4-xl4)
    END IF
    !   ### For T2
    ibTL = ibTL+(TP2-TP1)*nsig0
    AR2 = nsig0; xr4 = xsalog(ibTL+nsig0,igR)
    DO iut = nsig0,1,-1
      AR1 = iut; xl4 = xsalog(ibTL+iut,igR)
      IF (xut4.GE.xl4) EXIT
      AR2 = AR1; xr4 = xl4
    END DO
    IF (AR1 .EQ. AR2) THEN
      wAR = 1.
    ELSE
      wAR = (xr4-xut4)/(xr4-xl4)
    END IF
    ! &&& ================================================== &&& !

    ! &&& ========== RIF Factor Calculation & Sum =========== &&& !
    rifaT = 0.;
    DO jso = miso, liso
      jdxiso = mapnucl(jso)
      IF (jdxiso.EQ.0) CYCLE
      ibt = rifid(jdxiso,ireso)
      IF (ibt.EQ.0) CYCLE
      xut8 = pnum(jso)
      IF (xut8.LE.0.) CYCLE
      ! ### Weights for ND ratio
      LB = ptrRifRat(ibt)+1
      UB = ptrRifRat(ibt+1)

      xut8 = log(xut8/ND)
      R1 = LB; xl8 = ratlog(LB)
      DO iut = LB, UB
        R2 = iut; xr8 = ratlog(iut)
        IF (xut8<xr8) EXIT
        R1 = R2; xl8 = xr8
      END DO
      IF (R1 .EQ. R2) THEN
        wR = 1.
      ELSE
        wR = (xr8-xut8)/(xr8-xl8)
      END IF

      ! ### (R1,T1)
      ibTL = ptrRifRx(R1)+(TP1-1)*nsig0
      rifaT(1) = rifaT(1)+(rifabs(ibTL+AL1,igR)*wAL+rifabs(ibTL+AL2,igR)*(1.-wAL))*wR
      ! ### (R1,T2)
      ibTL = ibTL+(TP2-TP1)*nsig0
      rifaT(2) = rifaT(2)+(rifabs(ibTL+AR1,igR)*wAR+rifabs(ibTL+AR2,igR)*(1.-wAR))*wR - 1.
      ! ### (R2,T2)
      ibTL = ibTL+(R2-R1)*nsig0*nrtemp
      rifaT(2) = rifaT(2)+(rifabs(ibTL+AR1,igR)*wAR+rifabs(ibTL+AR2,igR)*(1.-wAR))*(1.-wR)
      ! ### (R2,T1)
      ibTL = ibTL+(TP1-TP2)*nsig0
      rifaT(1) = rifaT(1)+(rifabs(ibTL+AL1,igR)*wAL+rifabs(ibTL+AL2,igR)*(1.-wAL))*(1.-wR) - 1.
    END DO
!    print*, iso-miso+1, igR, ifxrC
    ! &&& =================================================== &&& !
    xsa = xsa*ND*(rifaT(1)*wTP+rifaT(2)*(1.-wTP)+1.)
  END IF
  jso = iso-miso+1
  jt = itTmp(ifxr,jso); wtN = wtTmp(ifxr,jso)
  lskipit2 = (wtN .GT. 0.99999)
  jt2 = jt+1
  IF (lskipit2) jt2 = jt

  wtN = ND*wtN;
  sigbuf = siga(ig, jt); sigbuf2 = siga(ig,jt2)
  sigbuf2 = wtN*sigbuf+(ND-wtN)*sigbuf2
  micA = micA+sigbuf2
  IF (lreso) THEN
    effA = effA+xsa
    fresoAIso(isoRG) = xsa/sigbuf2
  ELSE
    effA = effA+sigbuf2
  END IF
END DO
fresoA(igR,ifxrR) = effA/micA

END SUBROUTINE

ATTRIBUTES(GLOBAL) SUBROUTINE cuEffIntegrated_AICT(IsoG2R,ptrRTemp,ptrNlv,ptrTempLv, &
                                lamsigp,lvabs,wgtabs,maclv_log,rtempsq,rifa,xseq_1g, &
                                AccNtR,AccNt,mapnucl,mapA2G,mapG2R,ipinT, &
                                pnum,TempAvg,siglo_nreso, &
                                fresoAIso,fresoA,siga,wtTmp,&
                                niso,itTmp,&
                                nAICT,bg,nrg,nmaclv,iz)
USE CUDAFOR
IMPLICIT NONE
INTEGER,DEVICE,DIMENSION(:) :: IsoG2R,ptrRTemp,ptrNlv,ptrTempLv
REAL(GPU_XS_PRECISION),DEVICE,DIMENSION(:,:) :: lamsigp,lvabs,wgtabs,maclv_log
REAL(GPU_XS_PRECISION),DEVICE,DIMENSION(:) :: rtempsq
REAL(GPU_RES_PRECISION),DEVICE:: rifa(:,:,:), xseq_1g(:,:)
INTEGER,DEVICE :: AccNtR(:),AccNt(:),mapnucl(:),mapA2G(:),mapG2R(:),ipinT(:)
REAL(GPU_PNUM_PRECISION),DEVICE :: pnum(:)
REAL(GPU_XS_PRECISION),DEVICE :: TempAvg(:),siglo_nreso(:,:)
INTEGER,VALUE :: nAICT, nrg, nmaclv, iz

REAL(GPU_RES_PRECISION),DEVICE :: fresoAIso(:),fresoA(:,:)
REAL(GPU_XS_PRECISION),DEVICE,DIMENSION(:,:) :: siga,wtTmp
INTEGER,DEVICE :: itTmp(:,:), niso(:)
INTEGER,VALUE :: bg

INTEGER :: ifxrA, ifxr, ifxrR, ipin, ig, igR, iso, isoRG, miso, liso, idxiso, ireso
INTEGER :: nrtemp
REAL(8) :: xut8, xl8, xr8
INTEGER :: T1, T2, Mc1, Mc2
REAL(8) :: wT, wMc
INTEGER :: iut, ibt, ilv, ibTL, LB, UB

REAL(8) :: lpNoR, siglpiso, ND, xsa, phi, lvX, sigb, wgtX

INTEGER :: jso
INTEGER :: jt, jt2
REAL(8) :: sigbuf, sigbuf2, wtN

REAL(8) :: micA, effA
LOGICAL :: lskipit2, lreso

ifxrA = threadIdx%y+EFF_BLOCK_RDIM*(blockIdx%y-1)
IF (ifxrA.GT.nAICT) RETURN
igR = threadIdx%x
IF (igR.GT.nrg) RETURN
ig = bg+igR

ifxr = mapA2G(ifxrA)
IF (ifxr.LT.ifxrbegT(iz)) RETURN
IF (ifxr.GE.ifxrbegT(iz+1)) RETURN
ifxrR = mapG2R(ifxr)

ipin = ipinT(ifxr)

lpNoR = siglo_nreso(igR,ifxrR)
miso = AccNt(ifxr)+1; liso = AccNt(ifxr+1)
isoRG = AccNtR(ifxrR)*nrg+(igR-1)*(liso-miso+1)

micA = 0.; effA = 0.

DO iso = miso, liso
  isoRG = isoRG+1
  ND = pnum(iso); idxiso = mapnucl(iso)
  ireso = IsoG2R(idxiso)
  lreso = (ireso.GT.0)
  IF (lreso) THEN
    ! &&& ========= AvgTemperature Weights Finding ========= &&& !
    xut8 = TempAvg(ipin); xut8 = sqrt(xut8)
    ibt = ptrRTemp(ireso); nrtemp = ptrRTemp(ireso+1)-ibt
    T1 = 1; xl8 = rtempsq(ibt+1)
    DO iut = 1, nrtemp
      T2 = iut; xr8 = rtempsq(ibt+iut)
      IF (xut8<xr8) EXIT
      T1 = T2; xl8 = xr8
    END DO
    IF (T1.EQ.T2) THEN
      wT = 1.
    ELSE
      wT = (xr8-xut8)/(xr8-xl8)
    END IF
    ! &&& ================================================== &&& !
    siglpiso = lpNoR+ND*lamsigp(igR,idxiso)
    ! &&& ========== XSA and XSF Calculation ============ &&& !
    LB = ptrNlv(ireso)+1; UB = ptrNlv(ireso+1);
    ibTL = ptrTempLv(ireso);
    xsa = 0.; phi = 1.;
    DO ilv = LB, UB
      lvX = lvabs(ilv,igR)
      xut8 = ND*lvX;
      xut8 = log(xut8)
      Mc1 = 1; xl8 = maclv_log(1,iz)
      DO iut = 1, nmaclv
        Mc2 = iut; xr8 = maclv_log(iut,iz)
        IF (xut8<xr8) EXIT
        Mc1 = Mc2; xl8 = xr8
      END DO
      IF (Mc1.EQ.Mc2) THEN
        wMc = 1.
      ELSE
        wMc = (xr8-xut8)/(xr8-xl8)
      END IF
      sigb = wMc*xseq_1g(Mc1,ifxrA)+(1.-wMc)*xseq_1g(Mc2,ifxrA)
      sigb = (sigb+siglpiso)/ND
      wgtX = wT*wgtabs(ibTL+T1,igR)+(1.-wT)*(wgtabs(ibTL+T2,igR))
      wgtX = wgtX*lvX
      lvX = 1./(lvX+sigb) ! phimult
      phi = phi-wgtX*lvX
      sigb = lvX*sigb ! phihom
      xsa = xsa+wgtX*sigb
      ibTL = ibTL+nrtemp
    END DO
    ! &&& =============================================== &&& !
    xsa = xsa/phi
    IF (xsa.LE.0.) lreso = .FALSE.
    xsa = xsa*ND*rifa(ireso,igR,ipin)
  END IF

  ND = pnum(iso)
  jso = iso-miso+1
  jt = itTmp(ifxr,jso); wtN = wtTmp(ifxr,jso)
  lskipit2 = (wtN .GT. 0.99999)
  jt2 = jt+1
  IF (lskipit2) jt2 = jt

  wtN = ND*wtN;
  sigbuf = siga(ig, jt); sigbuf2 = siga(ig,jt2)
  sigbuf2 = wtN*sigbuf+(ND-wtN)*sigbuf2
  micA = micA+sigbuf2
  IF (lreso) THEN
    effA = effA+xsa
    fresoAIso(isoRG) = xsa/sigbuf2
  ELSE
    effA = effA+sigbuf2
  END IF
END DO
fresoA(igR,ifxrR) = effA/micA
END SUBROUTINE

SUBROUTINE EffMacGen_cuEff(iz)
USE Core_mod, ONLY : GroupInfo
USE geom, ONLY : core
USE PE_Mod, ONLY : PE
USE XSLIB_MOD, ONLY : nreshel
IMPLICIT NONE
INTEGER :: iz, jz
INTEGER :: npin,nResD,nResTri,nGenD,nCldT,nAICT,rgb,rge,nrg,ntiso
REAL, ALLOCATABLE :: siglpD(:,:), siglpT(:,:)
REAL(GPU_XS_PRECISION), DEVICE, ALLOCATABLE :: siglpD_dev(:,:), siglpT_dev(:,:)

INTEGER :: isync, ierr, istat

!!REAL(4), ALLOCATABLE :: rifa(:,:,:),riff(:,:,:)
REAL(4), ALLOCATABLE :: isoA(:,:,:),isof(:,:,:)
!!REAL(4), ALLOCATABLE :: isoA(:)
REAL(4), ALLOCATABLE :: fresoA(:,:), fresoF(:,:), fresoNF(:,:)
INTEGER :: iso, ig, ifxr, ifxrC, iRG

jz = iz-PE%myzb+1
rgb = GroupInfo%nofg+1; rge = GroupInfo%nofg+GroupInfo%norg
nrg = GroupInfo%norg; ntiso = IsoData%ntiso
npin = core%nxy*(PE%myze-PE%myzb+1)

nResD =FxrD%nResD; nGenD = FxrD%nGenD
nCldT = FxrTri%nCldTri; nAICT = FxrTri%nAICTri; nResTri = FxrTri%nResTri

!print*, 'A'
IF (FxrD%nResD.GT.0) THEN
  ALLOCATE(siglpD(rgb:rge,FxrD%nResD))
  ALLOCATE(siglpD_dev(nrg,FxrD%nResD))
END IF
IF (FxrTri%nResTri.GT.0) THEN
  ALLOCATE(siglpT(rgb:rge,FxrTri%nResTri))
  ALLOCATE(siglpT_dev(nrg,FxrTri%nResTri))
END IF
!print*, 'B'
CALL CalSigPotwNoReso(siglpD,siglpT,.FALSE.,.FALSE.)
!print*, 'C'
IF (FxrD%nResD.GT.0) THEN
!  siglpD_dev = siglpD;
  istat = cudaMemCpy(siglpD_dev,siglpD,nrg*FxrD%nResD)
  DEALLOCATE(siglpD)
END IF
IF (FxrTri%nResTri.GT.0) THEN
!  siglpT_dev = siglpT
  istat = cudaMemCpy(siglpT_dev,siglpT,nrg*FxrTri%nResTri)
  DEALLOCATE(siglpT)
END IF
!print*, 'D'

CALL cuEffRIFPin<<<blockRP,trdEffIso,0,cuDevice%mystream>>>(cuResIsoData%lclad,cuResIsoData%ityp,cuResIsoData%IsoG2R,&
                cuResIsoData%ptrRTemp,cuResIsoData%ptrSig0,cuResIsoData%ptrNlv,&
                cuResIsoData%ptrRifRat,cuResIsoData%ptrTempPot,cuResIsoData%ptrTempLv,cuResIsoData%ptrRifRx,&
                cuResIsoData%lamsigp,cuResIsoData%lvabs,cuResIsoData%lvfis,cuResIsoData%wgtabs,cuResIsoData%wgtfis,&
                cuMLGData%f_maclv_log,cuMLGData%f_maclv1G_log,cuResIsoData%rtempsq,cuResIsoData%ratlog,&
                cuResIsoData%xsalog,cuResIsoData%rifabs,cuResIsoData%riffis,cuResPin%FnAdj,cuResPin%avgxseq_1g,&
                cuResPin%lAIC,cuResPin%lPinRes,cuResPin%rifa,cuResPin%riff,cuResPin%avgxseq_mg,&
                cuResPin%niso,cuResPin%mapnucl,cuResIsoData%rifid,&
                cuResPin%pnum,cuResPin%ind,cuResPin%temp,&
                ntiso,nrg,MLGLib%f_nmaclv,MLGLib%f_nmaclv1G,jz)

!isync = cudaDeviceSynchronize()
!ierr = cudaPeekAtLastError()
!if (ierr.NE.0) print*, cudaGetErrorString(ierr)
!if (ierr.NE.0) print*, __FILE__, __LINE__

!ALLOCATE(rifa(nreshel,nrg,npin),riff(nreshel,nrg,npin))
!rifa = cuResPin%rifa; riff = cuResPin%riff;
!DO ig = 1, nrg
!  DO ipin = 1, npin
!    IF (.NOT. ResPin%lres(ipin)) CYCLE
!    DO iso = 1, ResPin%niso(ipin)
!      ierr = ResPin%mapnucl(iso,ipin)
!      isync = ResIsoData%IsoG2R(ierr)
!      IF (isync.GT.0) THEN
!        IF (ResIsoData%lclad(ierr)) cycle
!        WRITE(91,*) ipin, ig*1000+iso, rifa(isync,ig,ipin)!, riff(isync,ig,ipin)
!!        WRITE(91,*) ierr*100+ig, ipin, rifa(ierr,ig,ipin)
!      END IF
!    END DO
!  END DO
!END DO
!deallocate(rifa,riff)
!FLUSH(91)

CALL cuEffMacXS_GenD<<<blockRG,trdEffIso,0,cuDevice%mystream>>>(cuResIsoData%ityp,cuResIsoData%IsoG2R,cuResIsoData%ptrRTemp,cuResIsoData%ptrNlv,cuResIsoData%ptrTempLv, &
                    cuResIsoData%lamsigp,cuResIsoData%lvabs,cuResIsoData%lvfis,cuResIsoData%wgtabs,cuResIsoData%wgtfis,&
                    cuMLGData%f_maclv_log,cuResIsoData%rtempsq,cuResPin%rifa,cuResPin%riff,cuBlockFxr%xseq_f_mgD,&
                    cuBlockFxr%MapNuclD,cuBlockFxr%pinidD,cuBlockFxr%mapGR2GD,cuBlockFxr%mapG2RD, &
                    cuBlockFxr%pnumD,cuBlockFxr%indD,cuResPin%temp,siglpD_dev,&
                    cuBlockFxr%FnAdjD,cuBlockFxr%FtAdjD,cuBlockFxr%FresoAIsoD,cuBlockFxr%FresoFIsoD, &
                    nGenD,ntiso,nrg, MLGLib%f_nmaclv, jz)
!isync = cudaDeviceSynchronize()
!ierr = cudaGetLastError()
!if (ierr.NE.0) print*, cudaGetErrorString(ierr)
!if (ierr.NE.0) print*, __FILE__, __LINE__

!ALLOCATE(isoA(IsoData%ntiso,nrg,nGenD),isof(IsoData%ntiso,nrg,nGenD))
!isoA = cuBlockFxr%FresoAIsoD; isof = cuBlockFxr%FresoFIsoD;
!DO ig = 1, nrg
!  DO ifxr = 1, nGenD
!    DO iso = 1, FxrD%niso(FxrD%mapR2G(ifxr))
!      write(150+ig,*) ifxr,ig*100+iso, isoA(iso,ig,ifxr)
!      write(250+ig,*) ifxr,ig*100+iso, isof(iso,ig,ifxr)
!    END DO
!  END DO
!  flush(150+ig)
!  flush(250+ig)
!END DO
!deallocate(isoA,isof)

CALL cuEffMacXS_CldT<<<blockRC,trdEffTri,0,cuDevice%mystream>>>(cuResIsoData%IsoG2R,cuResIsoData%ptrRTemp,cuResIsoData%ptrSig0,cuResIsoData%ptrNlv,&
                    cuResIsoData%ptrRifRat,cuResIsoData%ptrTempPot,cuResIsoData%ptrTempLv,cuResIsoData%ptrRifRx, &
                    cuResIsoData%lamsigp,cuResIsoData%lvabs,cuResIsoData%wgtabs,&
                    cuMLGData%c_maclv1G_log,cuResIsoData%rtempsq,cuResIsoData%ratlog,cuResIsoData%xsalog,cuResIsoData%rifabs,&
                    cuBlockFxr%xseq_c_1gT,cuBlockFxr%AccNtR,cuBlockFxr%AccNt,cuBlockFxr%MapNuclT,&
                    cuResIsoData%rifid,cuBlockFxr%mapC2GT,cuBlockFxr%mapG2RT,cuBlockFxr%pinidT,&
                    cuBlockFxr%pnumT,cuBlockFxr%tempT,siglpT_dev,cuBlockFxr%FresoAIsoT,&
                    nCldT, nrg, MLGLib%c_nmaclv1G, jz)
!isync = cudaDeviceSynchronize()
!ierr = cudaGetLastError()
!if (ierr.NE.0) print*, cudaGetErrorString(ierr)
!if (ierr.NE.0) print*, __FILE__, __LINE__

!allocate(isoA(FxrTri%NtRTri*nrg))
!isoA = cuBlockFxr%FresoAIsoT
!DO ig = 1, nrg
!  DO ifxrC = 1, nCldT
!    ifxr = FxrTri%mapC2G(ifxrC)
!    iRG = FxrTri%AccNtR(FxrTri%mapG2R(ifxr))*nrg+(ig-1)*(FxrTri%AccNt(ifxr+1)-FxrTri%AccNt(ifxr))
!    DO iso = FxrTri%AccNt(ifxr)+1, FxrTri%AccNt(ifxr)+FxrTri%niso(ifxr)
!      iRG = iRG+1
!      WRITE(91,*) ifxrC, ig*100+(iso-FxrTri%AccNt(ifxr)), isoA(iRG)
!    END DO
!  END DO
!END DO
!DEALLOCATE(isoA)

CALL cuEffMacXS_AICT<<<blockRA,trdEffTri,0,cuDevice%mystream>>>(cuResIsoData%IsoG2R,cuResIsoData%ptrRTemp,cuResIsoData%ptrNlv,cuResIsoData%ptrTempLv, &
                    cuResIsoData%lamsigp,cuResIsoData%lvabs,cuResIsoData%wgtabs,&
                    cuMLGData%f_maclv1G_log,cuResIsoData%rtempsq,cuResPin%rifa,cuBlockFxr%xseq_f_1gT, &
                    cuBlockFxr%AccNtR,cuBlockFxr%AccNt,cuBlockFxr%MapNuclT,cuBlockFxr%mapA2GT,cuBlockFxr%mapG2RT,&
                    cuBlockFxr%pinidT,cuBlockFxr%pnumT,cuResPin%temp,siglpT_dev,cuBlockFxr%FresoAIsoT, &
                    nAICT, nrg, MLGLib%f_nmaclv1G, jz)
!isync = cudaDeviceSynchronize()
!ierr = cudaGetLastError()
!if (ierr.NE.0) print*, cudaGetErrorString(ierr)
!if (ierr.NE.0) print*, __FILE__, __LINE__

CALL cuEffFreso_D<<<blockRD,trdEffRes,0,cuDevice%mystream>>>(cuBlockFxr%FresoAIsoD,cuBlockFxr%FresoFIsoD,cuBlockFxr%FresoAD,cuBlockFxr%FresoF,cuBlockFxr%FresoNF,&
                  cuResIsoData%ityp,cuIsoData%siga,cuIsoData%sigf,cuIsoData%signf,cuBlockFxr%wtTempD,&
                  cuBlockFxr%nisoD,cuBlockFxr%MapNuclD,cuBlockFxr%itTempD,cuBlockFxr%mapR2GD,cuResIsoData%IsoG2R,cuBlockFxr%pnumD,&
                  (rgb-1),nrg,nResD, jz)
!isync = cudaDeviceSynchronize()
!ierr = cudaGetLastError()
!if (ierr.NE.0) print*, cudaGetErrorString(ierr)
!if (ierr.NE.0) print*, __FILE__, __LINE__
!ALLOCATE(fresoA(nrg,FxrD%nResD),fresoF(nrg,FxrD%nResD),fresoNF(nrg,FxrD%nResD))
!fresoA = cuBlockFxr%FresoAD; fresoF = cuBlockFxr%FresoF; fresoNF = cuBlockFxr%FresoNF
!DO ig = 1, nrg
!  DO ifxr = 1, nResD
!    WRITE(100+ig,*) ifxr, ig, fresoA(ig, ifxr)
!    WRITE(200+ig,*) ifxr, ig, fresoF(ig, ifxr)
!    WRITE(300+ig,*) ifxr, ig, fresoNF(ig, ifxr)
!  END DO
!  FLUSH(100+ig)
!  FLUSH(200+ig)
!  FLUSH(300+ig)
!END DO
!deallocate(fresoA)
!ALLOCATE(isoA(IsoData%ntiso,nrg,nGenD),isof(IsoData%ntiso,nrg,nGenD))
!isoA = cuBlockFxr%FresoAIsoD; isof = cuBlockFxr%FresoFIsoD;
!DO ig = 1, nrg
!  DO ifxr = 1, nGenD
!    DO iso = 1, FxrD%niso(FxrD%mapR2G(ifxr))
!      write(150+ig,*) ifxr,ig*100+iso, isoA(iso,ig,ifxr)
!      write(250+ig,*) ifxr,ig*100+iso, isof(iso,ig,ifxr)
!    END DO
!  END DO
!  flush(150+ig)
!  flush(250+ig)
!END DO
!deallocate(isoA,isof)

CALL cuEffFreso_T<<<blockRT,trdEffRes,0,cuDevice%mystream>>>(cuBlockFxr%FresoAIsoT,cuBlockFxr%FresoAT,cuIsoData%siga,cuBlockFxr%wtTempT,&
                  cuBlockFxr%AccNtR,cuBlockFxr%AccNt,cuBlockFxr%nisoT,cuBlockFxr%MapNuclT,&
                  cuBlockFxr%itTempT,cuBlockFxr%mapR2GT,cuResIsoData%IsoG2R,cuBlockFxr%pnumT,&
                  (rgb-1),nrg,nResTri, jz)
!isync = cudaDeviceSynchronize()
!ierr = cudaGetLastError()
!if (ierr.NE.0) print*, cudaGetErrorString(ierr)
!if (ierr.NE.0) print*, __FILE__, __LINE__
!ALLOCATE(fresoA(nrg,nResTri))
!fresoA = cuBlockFxr%FresoAT
!DO ig = 1, nrg
!  DO ifxr = 1, nResTri
!    WRITE(95,*) ifxr, ig, fresoA(ig, ifxr)
!  END DO
!END DO
!FLUSH(95)
!deallocate(fresoA)

!print*, 'F'
!stop

IF (FxrD%nResD.GT.0) THEN
  DEALLOCATE(siglpD_dev)
END IF
IF (FxrTri%nResTri.GT.0) THEN
  DEALLOCATE(siglpT_dev)
END IF

END SUBROUTINE

SUBROUTINE EffMacIntegrated_cuEff(iz)
USE Core_mod, ONLY : GroupInfo
USE geom, ONLY : core
USE PE_Mod, ONLY : PE
USE XSLIB_MOD, ONLY : nreshel
IMPLICIT NONE
INTEGER :: iz, jz
INTEGER :: npin,nResD,nResTri,nGenD,nCldT,nAICT,rgb,rge,nrg,ntiso
REAL, ALLOCATABLE :: siglpD(:,:), siglpT(:,:)
REAL(GPU_XS_PRECISION), DEVICE, ALLOCATABLE :: siglpD_dev(:,:), siglpT_dev(:,:)

INTEGER :: isync, ierr, istat

jz = iz-PE%myzb+1
rgb = GroupInfo%nofg+1; rge = GroupInfo%nofg+GroupInfo%norg
nrg = GroupInfo%norg; ntiso = IsoData%ntiso
npin = core%nxy*(PE%myze-PE%myzb+1)

nResD =FxrD%nResD; nGenD = FxrD%nGenD
nCldT = FxrTri%nCldTri; nAICT = FxrTri%nAICTri; nResTri = FxrTri%nResTri

IF (FxrD%nResD.GT.0) THEN
  ALLOCATE(siglpD(rgb:rge,FxrD%nResD))
  ALLOCATE(siglpD_dev(nrg,FxrD%nResD))
END IF
IF (FxrTri%nResTri.GT.0) THEN
  ALLOCATE(siglpT(rgb:rge,FxrTri%nResTri))
  ALLOCATE(siglpT_dev(nrg,FxrTri%nResTri))
END IF
CALL CalSigPotwNoReso(siglpD,siglpT,.FALSE.,.FALSE.)
IF (FxrD%nResD.GT.0) THEN
  siglpD_dev = siglpD;
  DEALLOCATE(siglpD)
END IF
IF (FxrTri%nResTri.GT.0) THEN
  siglpT_dev = siglpT
  DEALLOCATE(siglpT)
END IF

CALL cuEffRIFPin<<<blockRP,trdEffIso,0,cuDevice%mystream>>>(cuResIsoData%lclad,cuResIsoData%ityp,cuResIsoData%IsoG2R,&
                cuResIsoData%ptrRTemp,cuResIsoData%ptrSig0,cuResIsoData%ptrNlv,&
                cuResIsoData%ptrRifRat,cuResIsoData%ptrTempPot,cuResIsoData%ptrTempLv,cuResIsoData%ptrRifRx,&
                cuResIsoData%lamsigp,cuResIsoData%lvabs,cuResIsoData%lvfis,cuResIsoData%wgtabs,cuResIsoData%wgtfis,&
                cuMLGData%f_maclv_log,cuMLGData%f_maclv1G_log,cuResIsoData%rtempsq,cuResIsoData%ratlog,&
                cuResIsoData%xsalog,cuResIsoData%rifabs,cuResIsoData%riffis,cuResPin%FnAdj,cuResPin%avgxseq_1g,&
                cuResPin%lAIC,cuResPin%lPinRes,cuResPin%rifa,cuResPin%riff,cuResPin%avgxseq_mg,&
                cuResPin%niso,cuResPin%mapnucl,cuResIsoData%rifid,&
                cuResPin%pnum,cuResPin%ind,cuResPin%temp,&
                ntiso,nrg,MLGLib%f_nmaclv,MLGLib%f_nmaclv1G,jz)
!isync = cudaDeviceSynchronize()
!ierr = cudaPeekAtLastError()
!if (ierr.NE.0) print*, cudaGetErrorString(ierr)
!if (ierr.NE.0) print*, __FILE__, __LINE__

CALL cuEffIntegrated_GenD<<<blockRD,trdEffRes,0,cuDevice%myStream>>>(cuResIsoData%ityp,cuResIsoData%IsoG2R,&
                                cuResIsoData%ptrRTemp,cuResIsoData%ptrNlv,cuResIsoData%ptrTempLv, &
                                cuResIsoData%lamsigp,cuResIsoData%lvabs,cuResIsoData%lvfis,&
                                cuResIsoData%wgtabs,cuResIsoData%wgtfis,cuMLGData%f_maclv_log,cuResIsoData%rtempsq, &
                                cuResPin%rifa,cuResPin%riff,cuBlockFxr%xseq_f_mgD,cuBlockFxr%MapNuclD,&
                                cuBlockFxr%pinidD,cuBlockFxr%mapGR2GD,cuBlockFxr%mapG2RD, &
                                cuBlockFxr%pnumD,cuBlockFxr%indD,cuResPin%temp,siglpD_dev,&
                                cuBlockFxr%FnAdjD,cuBlockFxr%FtAdjD, &
                                cuBlockFxr%FresoAIsoD,cuBlockFxr%FresoFIsoD,cuBlockFxr%FresoAD,&
                                cuBlockFxr%FresoF,cuBlockFxr%FresoNF, &
                                cuIsoData%siga,cuIsoData%sigf,cuIsoData%signf,cuBlockFxr%wtTempD,&
                                cuBlockFxr%nisoD,cuBlockFxr%itTempD,&
                                nGenD,(rgb-1),rgb,rge,1,nResD,MLGLib%f_nmaclv,jz)
!isync = cudaDeviceSynchronize()
!ierr = cudaPeekAtLastError()
!if (ierr.NE.0) print*, cudaGetErrorString(ierr)
!if (ierr.NE.0) print*, __FILE__, __LINE__
!print*, 'GenD'

CALL cuEffIntegrated_CldT<<<blockRT,trdEffRes,0,cuDevice%myStream>>>(cuResIsoData%IsoG2R,cuResIsoData%ptrRTemp,&
                                cuResIsoData%ptrSig0,cuResIsoData%ptrNlv,cuResIsoData%ptrRifRat,&
                                cuResIsoData%ptrTempPot,cuResIsoData%ptrTempLv,cuResIsoData%ptrRifRx, &
                                cuResIsoData%lamsigp,cuResIsoData%lvabs,cuResIsoData%wgtabs,&
                                cuMLGData%c_maclv1G_log,cuResIsoData%rtempsq,cuResIsoData%ratlog,&
                                cuResIsoData%xsalog,cuResIsoData%rifabs,cuBlockFxr%xseq_c_1gT, &
                                cuBlockFxr%AccNtR,cuBlockFxr%AccNt,cuBlockFxr%MapNuclT,&
                                cuResIsoData%rifid,cuBlockFxr%mapC2GT,cuBlockFxr%mapG2RT,cuBlockFxr%pinidT,&
                                cuBlockFxr%pnumT,cuBlockFxr%tempT,siglpT_dev, &
                                cuBlockFxr%FresoAIsoT,cuBlockFxr%FresoAT,&
                                cuIsoData%siga,cuBlockFxr%wtTempT,cuBlockFxr%nisoT,cuBlockFxr%itTempT,&
                                nCldT,(rgb-1),nrg,MLGLib%c_nmaclv1G,jz)
!isync = cudaDeviceSynchronize()
!ierr = cudaPeekAtLastError()
!if (ierr.NE.0) print*, cudaGetErrorString(ierr)
!if (ierr.NE.0) print*, __FILE__, __LINE__
!!print*, 'CLDT'

CALL cuEffIntegrated_AICT<<<blockRT,trdEffRes,0,cuDevice%myStream>>>(cuResIsoData%IsoG2R,cuResIsoData%ptrRTemp,&
                                cuResIsoData%ptrNlv,cuResIsoData%ptrTempLv,cuResIsoData%lamsigp, &
                                cuResIsoData%lvabs,cuResIsoData%wgtabs,cuMLGData%f_maclv1G_log,&
                                cuResIsoData%rtempsq,cuResPin%rifa,cuBlockFxr%xseq_f_1gT, &
                                cuBlockFxr%AccNtR,cuBlockFxr%AccNt,cuBlockFxr%mapnuclT,&
                                cuBlockFxr%mapA2GT,cuBlockFxr%mapG2RT,cuBlockFxr%pinidT, &
                                cuBlockFxr%pnumT,cuResPin%temp,siglpT_dev, &
                                cuBlockFxr%FresoAIsoT,cuBlockFxr%FresoAT, &
                                cuIsoData%siga,cuBlockFxr%wtTempT,cuBlockFxr%nisoT,cuBlockFxr%itTempT,&
                                nAICT,(rgb-1),nrg,MLGLib%f_nmaclv1G,jz)
!isync = cudaDeviceSynchronize()
!ierr = cudaPeekAtLastError()
!if (ierr.NE.0) print*, cudaGetErrorString(ierr)
!if (ierr.NE.0) print*, __FILE__, __LINE__
!!print*, 'AICT'

IF (FxrD%nResD.GT.0) THEN
  DEALLOCATE(siglpD_dev)
END IF
IF (FxrTri%nResTri.GT.0) THEN
  DEALLOCATE(siglpT_dev)
END IF

END SUBROUTINE

SUBROUTINE EffMacIntegrated_RIP(iz)
USE Core_mod, ONLY : GroupInfo
USE PE_Mod, ONLY : PE
IMPLICIT NONE
INTEGER :: iz, jz
INTEGER :: nrg,ntiso

jz = iz-PE%myzb+1
nrg = GroupInfo%norg; ntiso = IsoData%ntiso

CALL cuEffRIFPin<<<blockRP,trdEffIso,0,cuDevice%mystream>>>(cuResIsoData%lclad,cuResIsoData%ityp,cuResIsoData%IsoG2R,&
                cuResIsoData%ptrRTemp,cuResIsoData%ptrSig0,cuResIsoData%ptrNlv,&
                cuResIsoData%ptrRifRat,cuResIsoData%ptrTempPot,cuResIsoData%ptrTempLv,cuResIsoData%ptrRifRx,&
                cuResIsoData%lamsigp,cuResIsoData%lvabs,cuResIsoData%lvfis,cuResIsoData%wgtabs,cuResIsoData%wgtfis,&
                cuMLGData%f_maclv_log,cuMLGData%f_maclv1G_log,cuResIsoData%rtempsq,cuResIsoData%ratlog,&
                cuResIsoData%xsalog,cuResIsoData%rifabs,cuResIsoData%riffis,cuResPin%FnAdj,cuResPin%avgxseq_1g,&
                cuResPin%lAIC,cuResPin%lPinRes,cuResPin%rifa,cuResPin%riff,cuResPin%avgxseq_mg,&
                cuResPin%niso,cuResPin%mapnucl,cuResIsoData%rifid,&
                cuResPin%pnum,cuResPin%ind,cuResPin%temp,&
                ntiso,nrg,MLGLib%f_nmaclv,MLGLib%f_nmaclv1G,jz)
END SUBROUTINE

SUBROUTINE EffMacIntegrated_GenD_GB(iz)
USE Core_mod, ONLY : GroupInfo
USE PE_Mod, ONLY : PE
IMPLICIT NONE
INTEGER :: iz, jz
INTEGER :: nResD,nGenD,rgb,rge,nrg,ntiso,gb,ge
REAL, ALLOCATABLE :: siglpD(:,:), dummy(:,:)
REAL(GPU_XS_PRECISION), DEVICE, ALLOCATABLE :: siglpD_dev(:,:)
REAL(GPU_RES_PRECISION), DEVICE, ALLOCATABLE :: FresoAIso(:,:,:), FresoFIso(:,:,:)
REAL(GPU_RES_PRECISION), DEVICE, ALLOCATABLE :: FresoA(:,:), FresoF(:,:), FresoNF(:,:)
!REAL(GPU_RES_PRECISION), PINNED, ALLOCATABLE :: FresoAIsoHost(:,:,:), FresoFIsoHost(:,:,:)
!REAL(GPU_RES_PRECISION), PINNED, ALLOCATABLE :: FresoAHost(:,:), FresoFHost(:,:), FresoNFHost(:,:)
REAL(GPU_RES_PRECISION), ALLOCATABLE :: FresoAIsoHost(:,:,:), FresoFIsoHost(:,:,:)
REAL(GPU_RES_PRECISION), ALLOCATABLE :: FresoAHost(:,:), FresoFHost(:,:), FresoNFHost(:,:)

INTEGER :: igb, ngBlock, mg

INTEGER :: istat

jz = iz-PE%myzb+1
rgb = GroupInfo%nofg+1; rge = GroupInfo%nofg+GroupInfo%norg
nrg = GroupInfo%norg; ntiso = IsoData%ntiso

nResD =FxrD%nResD; nGenD = FxrD%nGenD

IF (FxrD%nResD.GT.0) THEN
  ALLOCATE(siglpD(rgb:rge,FxrD%nResD))
  ALLOCATE(siglpD_dev(nrg,FxrD%nResD))
END IF
CALL CalSigPotwNoReso(siglpD,dummy,.FALSE.,.TRUE.)
IF (FxrD%nResD.GT.0) THEN
  siglpD_dev = siglpD;
  DEALLOCATE(siglpD)
END IF

ngBlock = (nrg-1)/SG_BLOCK_SIZE+1

DO igb = 1, ngBlock
  gb = (igb-1)*SG_BLOCK_SIZE+rgb; ge = min(nrg, igb*SG_BLOCK_SIZE)+rgb-1
  !gb = rgb; ge = rge
  mg = ge-gb+1
  ALLOCATE(FresoAIsoHost(ntiso,gb:ge,nResD))
  ALLOCATE(FresoFIsoHost(ntiso,gb:ge,nResD))
  ALLOCATE(FresoAHost(gb:ge,nResD))
  ALLOCATE(FresoFHost(gb:ge,nResD))
  ALLOCATE(FresoNFHost(gb:ge,nResD))
  ALLOCATE(FresoAIso(ntiso,gb:ge,nResD))
  ALLOCATE(FresoFIso(ntiso,gb:ge,nResD))
  ALLOCATE(FresoA(gb:ge,nResD))
  ALLOCATE(FresoF(gb:ge,nResD))
  ALLOCATE(FresoNF(gb:ge,nResD))

  FresoAHost(gb:ge,:) = FxrD%fresoa(gb:ge,:)
  FresoFHost(gb:ge,:) = FxrD%fresof(gb:ge,:)
  FresoNFHost(gb:ge,:) = FxrD%fresonf(gb:ge,:)
  FresoAIsoHost(:,gb:ge,:) = FxrD%fresoAIso(:,gb:ge,:)
  FresoFIsoHost(:,gb:ge,:) = FxrD%fresoFIso(:,gb:ge,:)

  istat = cudaMemcpy(FresoAIso, FresoAIsoHost, ntiso*mg*nResD, cudaMemcpyHostToDevice)
  istat = cudaMemcpy(FresoFIso, FresoFIsoHost, ntiso*mg*nResD, cudaMemcpyHostToDevice)
  istat = cudaMemcpy(FresoA, FresoAHost,mg*nResD, cudaMemcpyHostToDevice)
  istat = cudaMemcpy(FresoF, FresoFHost,mg*nResD, cudaMemcpyHostToDevice)
  istat = cudaMemcpy(FresoNF, FresoNFHost,mg*nResD, cudaMemcpyHostToDevice)

  CALL cuEffIntegrated_GenD<<<blockRD,trdEffRes,0,cuDevice%myStream>>>(cuResIsoData%ityp,cuResIsoData%IsoG2R,&
                                  cuResIsoData%ptrRTemp,cuResIsoData%ptrNlv,cuResIsoData%ptrTempLv, &
                                  cuResIsoData%lamsigp,cuResIsoData%lvabs,cuResIsoData%lvfis,&
                                  cuResIsoData%wgtabs,cuResIsoData%wgtfis,cuMLGData%f_maclv_log,cuResIsoData%rtempsq, &
                                  cuResPin%rifa,cuResPin%riff,cuBlockFxr%xseq_f_mgD,cuBlockFxr%MapNuclD,&
                                  cuBlockFxr%pinidD,cuBlockFxr%mapGR2GD,cuBlockFxr%mapG2RD, &
                                  cuBlockFxr%pnumD,cuBlockFxr%indD,cuResPin%temp,siglpD_dev,&
                                  cuBlockFxr%FnAdjD,cuBlockFxr%FtAdjD, &
                                  FresoAIso, FresoFIso, FresoA, FresoF, FresoNF, &
  !                                cuBlockFxr%FresoAIsoD,cuBlockFxr%FresoFIsoD,cuBlockFxr%FresoAD,&
  !                                cuBlockFxr%FresoF,cuBlockFxr%FresoNF, &
                                  cuIsoData%siga,cuIsoData%sigf,cuIsoData%signf,cuBlockFxr%wtTempD,&
                                  cuBlockFxr%nisoD,cuBlockFxr%itTempD,&
                                  nGenD,(rgb-1),gb,ge,1,nResD,MLGLib%f_nmaclv,jz)

  istat = cudaMemcpy(FresoAIsoHost,FresoAIso, ntiso*mg*nResD, cudaMemcpyDeviceToHost)
  istat = cudaMemcpy(FresoFIsoHost,FresoFIso, ntiso*mg*nResD, cudaMemcpyDeviceToHost)
  istat = cudaMemcpy(FresoAHost,FresoA, mg*nResD, cudaMemcpyDeviceToHost)
  istat = cudaMemcpy(FresoFHost,FresoF, mg*nResD, cudaMemcpyDeviceToHost)
  istat = cudaMemcpy(FresoNFHost,FresoNF, mg*nResD, cudaMemcpyDeviceToHost)

  DEALLOCATE(FresoA,FresoF,FresoNF,FresoAIso,FresoFIso)
  FxrD%fresoa(gb:ge,:) = FresoAHost(gb:ge,:)
  FxrD%fresokf(gb:ge,:) = FresoFHost(gb:ge,:)
  FxrD%fresonf(gb:ge,:) = FresoNFHost(gb:ge,:)
  FxrD%fresoAIso(:,gb:ge,:) = FresoAIsoHost(:,gb:ge,:)
  FxrD%fresoFIso(:,gb:ge,:) = FresoFIsoHost(:,gb:ge,:)
  DEALLOCATE(FresoAHost,FresoFHost,FresoNFHost,FresoAIsoHost,FresoFIsoHost)
END DO

IF (FxrD%nResD.GT.0) THEN
  DEALLOCATE(siglpD_dev)
END IF

END SUBROUTINE

SUBROUTINE EffMacIntegrated_GenD_RB(iz)
USE Core_mod, ONLY : GroupInfo
USE PE_Mod, ONLY : PE
USE timer, ONLY : TimeChk, nTracer_dclock
USE OMP_LIB
IMPLICIT NONE
INTEGER :: iz, jz
INTEGER :: nResD,nGenD,rgb,rge,nrg,ntiso
REAL(GPU_XS_PRECISION), PINNED, ALLOCATABLE :: siglpD(:,:), dummy(:,:)
REAL(GPU_XS_PRECISION), DEVICE, ALLOCATABLE :: siglpD_dev(:,:)
REAL(GPU_RES_PRECISION), DEVICE, ALLOCATABLE :: FresoAIso(:,:,:), FresoFIso(:,:,:)
REAL(GPU_RES_PRECISION), DEVICE, ALLOCATABLE :: FresoA(:,:), FresoF(:,:), FresoNF(:,:)
!REAL(GPU_RES_PRECISION), PINNED, ALLOCATABLE :: FresoAIsoHost(:,:,:), FresoFIsoHost(:,:,:)
!REAL(GPU_RES_PRECISION), PINNED, ALLOCATABLE :: FresoAHost(:,:), FresoFHost(:,:), FresoNFHost(:,:)
TYPE BlockFreso_type
  REAL(GPU_RES_PRECISION), PINNED, ALLOCATABLE :: FresoAIsoPIN(:,:,:), FresoFIsoPIN(:,:,:)
  REAL(GPU_RES_PRECISION), PINNED, ALLOCATABLE :: FresoAPIN(:,:), FresoFPIN(:,:), FresoNFPIN(:,:)
END TYPE
TYPE (BlockFreso_type), ALLOCATABLE :: FresoHost(:)
INTEGER :: irb, nrBlock, mr0, mr, rbeg, rend, ireg

LOGICAL :: lAsync = .TRUE.
!LOGICAL :: lAsync = .FALSE.

INTEGER :: istat

REAL :: tbo, teo, telap

jz = iz-PE%myzb+1
rgb = GroupInfo%nofg+1; rge = GroupInfo%nofg+GroupInfo%norg
nrg = GroupInfo%norg; ntiso = IsoData%ntiso

nResD =FxrD%nResD; nGenD = FxrD%nGenD

tbo = nTracer_dclock(.FALSE.,.FALSE.)

IF (FxrD%nResD.GT.0) THEN
  ALLOCATE(siglpD(rgb:rge,FxrD%nResD))
  ALLOCATE(siglpD_dev(nrg,FxrD%nResD))
END IF
CALL CalSigPotwNoReso(siglpD,dummy,.FALSE.,.TRUE.)
IF (FxrD%nResD.GT.0) THEN
  istat = cudaMemcpyAsync(siglpD_dev,siglpD,nrg*nResD,cudaMemcpyHostToDevice,cuDevice%myStream)
END IF

teo = nTracer_dclock(.FALSE.,.FALSE.)
telap = teo-tbo
!CALL MPI_MAX_REAL(telap, PE%MPI_NTRACER_COMM, .FALSE.)
TimeChk%cuXSPreTime = TimeChk%cuXSPreTime+telap
tbo = teo

nrBlock = XS_NUM_BLOCK
mr0 = nResD/XS_NUM_BLOCK+1
ALLOCATE(FresoAIso(ntiso,rgb:rge,mr0))
ALLOCATE(FresoFIso(ntiso,rgb:rge,mr0))
ALLOCATE(FresoA(rgb:rge,mr0))
ALLOCATE(FresoF(rgb:rge,mr0))
ALLOCATE(FresoNF(rgb:rge,mr0))
IF (lAsync) THEN
  ALLOCATE(FresoHost(nrBlock))
  DO irb = 1, nrBlock
    ALLOCATE(FresoHost(irb)%FresoAIsoPIN(ntiso,rgb:rge,mr0))
    ALLOCATE(FresoHost(irb)%FresoFIsoPIN(ntiso,rgb:rge,mr0))
    ALLOCATE(FresoHost(irb)%FresoAPIN(rgb:rge,mr0))
    ALLOCATE(FresoHost(irb)%FresoFPIN(rgb:rge,mr0))
    ALLOCATE(FresoHost(irb)%FresoNFPIN(rgb:rge,mr0))
  END DO
END IF
!teo = omp_get_wtime()
!print*, __FILE__, __LINE__, teo-tbo

!$ call omp_set_num_threads(PE%nThread)

DO irb = 1, nrBlock
  rbeg = (irb-1)*mr0+1; rend = min(nResD, irb*mr0)
  !gb = rgb; ge = rge
  mr = rend-rbeg+1

  IF (lAsync) THEN
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    DO ireg = 1, mr
    FresoHost(irb)%FresoAIsoPIN(:,:,ireg) = FxrD%fresoAIso(:,:,rbeg-1+ireg)
    FresoHost(irb)%FresoFIsoPIN(:,:,ireg) = FxrD%fresoFIso(:,:,rbeg-1+ireg)
    FresoHost(irb)%FresoAPIN(:,ireg) = FxrD%fresoa(:,rbeg-1+ireg)
    FresoHost(irb)%FresoFPIN(:,ireg) = FxrD%fresof(:,rbeg-1+ireg)
    FresoHost(irb)%FresoNFPIN(:,ireg) = FxrD%fresonf(:,rbeg-1+ireg)
    END DO
    !$OMP END PARALLEL DO
!    teo = omp_get_wtime()
!    print*, __FILE__, __LINE__, teo-tbo

    istat = cudaMemcpyAsync(FresoAIso, FresoHost(irb)%FresoAIsoPIN, ntiso*nrg*mr, cudaMemcpyHostToDevice, cuDevice%myStream)
    istat = cudaMemcpyAsync(FresoFIso, FresoHost(irb)%FresoFIsoPIN, ntiso*nrg*mr, cudaMemcpyHostToDevice, cuDevice%myStream)
    istat = cudaMemcpyAsync(FresoA, FresoHost(irb)%FresoAPIN,nrg*mr, cudaMemcpyHostToDevice, cuDevice%myStream)
    istat = cudaMemcpyAsync(FresoF, FresoHost(irb)%FresoFPIN,nrg*mr, cudaMemcpyHostToDevice, cuDevice%myStream)
    istat = cudaMemcpyAsync(FresoNF, FresoHost(irb)%FresoNFPIN,nrg*mr, cudaMemcpyHostToDevice, cuDevice%myStream)
  ELSE
    istat = cudaMemcpy(FresoAIso, FxrD%fresoAIso(1:,1:,rbeg:), ntiso*nrg*mr, cudaMemcpyHostToDevice)
    istat = cudaMemcpy(FresoFIso, FxrD%fresoFIso(1:,1:,rbeg:), ntiso*nrg*mr, cudaMemcpyHostToDevice)
    istat = cudaMemcpy(FresoA, FxrD%fresoa(1:,rbeg:),nrg*mr, cudaMemcpyHostToDevice)
    istat = cudaMemcpy(FresoF, FxrD%fresof(1:,rbeg:),nrg*mr, cudaMemcpyHostToDevice)
    istat = cudaMemcpy(FresoNF, FxrD%fresonf(1:,rbeg:),nrg*mr, cudaMemcpyHostToDevice)
!    teo = omp_get_wtime()
!    print*, __FILE__, __LINE__, teo-tbo
  END IF


  CALL cuEffIntegrated_GenD<<<blockRD,trdEffRes,0,cuDevice%myStream>>>(cuResIsoData%ityp,cuResIsoData%IsoG2R,&
                                  cuResIsoData%ptrRTemp,cuResIsoData%ptrNlv,cuResIsoData%ptrTempLv, &
                                  cuResIsoData%lamsigp,cuResIsoData%lvabs,cuResIsoData%lvfis,&
                                  cuResIsoData%wgtabs,cuResIsoData%wgtfis,cuMLGData%f_maclv_log,cuResIsoData%rtempsq, &
                                  cuResPin%rifa,cuResPin%riff,cuBlockFxr%xseq_f_mgD,cuBlockFxr%MapNuclD,&
                                  cuBlockFxr%pinidD,cuBlockFxr%mapGR2GD,cuBlockFxr%mapG2RD, &
                                  cuBlockFxr%pnumD,cuBlockFxr%indD,cuResPin%temp,siglpD_dev,&
                                  cuBlockFxr%FnAdjD,cuBlockFxr%FtAdjD, &
                                  FresoAIso, FresoFIso, FresoA, FresoF, FresoNF, &
  !                                cuBlockFxr%FresoAIsoD,cuBlockFxr%FresoFIsoD,cuBlockFxr%FresoAD,&
  !                                cuBlockFxr%FresoF,cuBlockFxr%FresoNF, &
                                  cuIsoData%siga,cuIsoData%sigf,cuIsoData%signf,cuBlockFxr%wtTempD,&
                                  cuBlockFxr%nisoD,cuBlockFxr%itTempD,&
                                  nGenD,(rgb-1),rgb,rge,rbeg,rend,MLGLib%f_nmaclv,jz)
  IF (lAsync) THEN
    istat = cudaMemcpyAsync(FresoHost(irb)%FresoAIsoPIN,FresoAIso, ntiso*nrg*mr, cudaMemcpyDeviceToHost, cuDevice%myStream)
    istat = cudaMemcpyAsync(FresoHost(irb)%FresoFIsoPIN,FresoFIso, ntiso*nrg*mr, cudaMemcpyDeviceToHost, cuDevice%myStream)
    istat = cudaMemcpyAsync(FresoHost(irb)%FresoAPIN,FresoA, nrg*mr, cudaMemcpyDeviceToHost, cuDevice%myStream)
    istat = cudaMemcpyAsync(FresoHost(irb)%FresoFPIN,FresoF, nrg*mr, cudaMemcpyDeviceToHost, cuDevice%myStream)
    istat = cudaMemcpyAsync(FresoHost(irb)%FresoNFPIN,FresoNF, nrg*mr, cudaMemcpyDeviceToHost, cuDevice%myStream)
  ELSE
    istat = cudaMemcpy(FxrD%fresoAIso(1:,1:,rbeg:),FresoAIso, ntiso*nrg*mr, cudaMemcpyDeviceToHost)
    istat = cudaMemcpy(FxrD%fresoFIso(1:,1:,rbeg:),FresoFIso, ntiso*nrg*mr, cudaMemcpyDeviceToHost)
    istat = cudaMemcpy(FxrD%fresoa(1:,rbeg:),FresoA, nrg*mr, cudaMemcpyDeviceToHost)
    istat = cudaMemcpy(FxrD%fresokf(1:,rbeg:),FresoF, nrg*mr, cudaMemcpyDeviceToHost)
    istat = cudaMemcpy(FxrD%fresonf(1:,rbeg:),FresoNF, nrg*mr, cudaMemcpyDeviceToHost)
  END IF
!  teo = omp_get_wtime()
!  print*, __FILE__, __LINE__, teo-tbo
END DO

IF (lAsync) THEN
  istat = cudaStreamSynchronize(cuDevice%myStream)
  DO irb = 1, nrBlock
    rbeg = (irb-1)*mr0+1; rend = min(nResD, irb*mr0)
    !gb = rgb; ge = rge
    mr = rend-rbeg+1
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    DO ireg = 1, mr
    FxrD%fresoAIso(:,:,rbeg-1+ireg)=FresoHost(irb)%FresoAIsoPIN(:,:,ireg)
    FxrD%fresoFIso(:,:,rbeg-1+ireg)=FresoHost(irb)%FresoFIsoPIN(:,:,ireg)
    FxrD%fresoa(:,rbeg-1+ireg) = FresoHost(irb)%FresoAPIN(:,ireg)
    FxrD%fresokf(:,rbeg-1+ireg) = FresoHost(irb)%FresoFPIN(:,ireg)
    FxrD%fresonf(:,rbeg-1+ireg) = FresoHost(irb)%FresoNFPIN(:,ireg)
    END DO
    !$OMP END PARALLEL DO
  END DO
!  teo = omp_get_wtime()
!  print*, __FILE__, __LINE__, teo-tbo
END IF

teo = nTracer_dclock(.FALSE.,.FALSE.)
telap = teo-tbo
!CALL MPI_MAX_REAL(telap, PE%MPI_NTRACER_COMM, .FALSE.)
TimeChk%cuXSMainTime = TimeChk%cuXSMainTime+telap

DEALLOCATE(FresoAIso,FresoFIso,FresoA,FresoF,FresoNF)
IF (lAsync) THEN
  DO irb = 1, nrBlock
    DEALLOCATE(FresoHost(irb)%FresoAIsoPIN)
    DEALLOCATE(FresoHost(irb)%FresoFIsoPIN)
    DEALLOCATE(FresoHost(irb)%FresoAPIN)
    DEALLOCATE(FresoHost(irb)%FresoFPIN)
    DEALLOCATE(FresoHost(irb)%FresoNFPIN)
  END DO
  DEALLOCATE(FresoHost)
END IF

IF (FxrD%nResD.GT.0) THEN
  DEALLOCATE(siglpD)
  DEALLOCATE(siglpD_dev)
END IF
!teo = omp_get_wtime()
!print*, __FILE__, __LINE__, teo-tbo

END SUBROUTINE

SUBROUTINE EffMacIntegrated_CLDnAICT(iz)
USE Core_mod, ONLY : GroupInfo
USE PE_Mod, ONLY : PE
IMPLICIT NONE
INTEGER :: iz, jz
INTEGER :: nResTri,nCldT,nAICT,rgb,rge,nrg,ntiso
REAL, ALLOCATABLE :: dummy(:,:), siglpT(:,:)
REAL(GPU_XS_PRECISION), DEVICE, ALLOCATABLE :: siglpT_dev(:,:)

INTEGER :: istat

jz = iz-PE%myzb+1
rgb = GroupInfo%nofg+1; rge = GroupInfo%nofg+GroupInfo%norg
nrg = GroupInfo%norg; ntiso = IsoData%ntiso

nCldT = FxrTri%nCldTri; nAICT = FxrTri%nAICTri; nResTri = FxrTri%nResTri

IF (FxrTri%nResTri.GT.0) THEN
  ALLOCATE(siglpT(rgb:rge,FxrTri%nResTri))
  ALLOCATE(siglpT_dev(nrg,FxrTri%nResTri))
END IF
CALL CalSigPotwNoReso(dummy,siglpT,.TRUE.,.FALSE.)
IF (FxrTri%nResTri.GT.0) THEN
  siglpT_dev = siglpT
  DEALLOCATE(siglpT)
END IF

CALL cuEffIntegrated_CldT<<<blockRT,trdEffRes,0,cuDevice%myStream>>>(cuResIsoData%IsoG2R,cuResIsoData%ptrRTemp,&
                                cuResIsoData%ptrSig0,cuResIsoData%ptrNlv,cuResIsoData%ptrRifRat,&
                                cuResIsoData%ptrTempPot,cuResIsoData%ptrTempLv,cuResIsoData%ptrRifRx, &
                                cuResIsoData%lamsigp,cuResIsoData%lvabs,cuResIsoData%wgtabs,&
                                cuMLGData%c_maclv1G_log,cuResIsoData%rtempsq,cuResIsoData%ratlog,&
                                cuResIsoData%xsalog,cuResIsoData%rifabs,cuBlockFxr%xseq_c_1gT, &
                                cuBlockFxr%AccNtR,cuBlockFxr%AccNt,cuBlockFxr%MapNuclT,&
                                cuResIsoData%rifid,cuBlockFxr%mapC2GT,cuBlockFxr%mapG2RT,cuBlockFxr%pinidT,&
                                cuBlockFxr%pnumT,cuBlockFxr%tempT,siglpT_dev, &
                                cuBlockFxr%FresoAIsoT,cuBlockFxr%FresoAT,&
                                cuIsoData%siga,cuBlockFxr%wtTempT,cuBlockFxr%nisoT,cuBlockFxr%itTempT,&
                                nCldT,(rgb-1),nrg,MLGLib%c_nmaclv1G,jz)

CALL cuEffIntegrated_AICT<<<blockRT,trdEffRes,0,cuDevice%myStream>>>(cuResIsoData%IsoG2R,cuResIsoData%ptrRTemp,&
                                cuResIsoData%ptrNlv,cuResIsoData%ptrTempLv,cuResIsoData%lamsigp, &
                                cuResIsoData%lvabs,cuResIsoData%wgtabs,cuMLGData%f_maclv1G_log,&
                                cuResIsoData%rtempsq,cuResPin%rifa,cuBlockFxr%xseq_f_1gT, &
                                cuBlockFxr%AccNtR,cuBlockFxr%AccNt,cuBlockFxr%mapnuclT,&
                                cuBlockFxr%mapA2GT,cuBlockFxr%mapG2RT,cuBlockFxr%pinidT, &
                                cuBlockFxr%pnumT,cuResPin%temp,siglpT_dev, &
                                cuBlockFxr%FresoAIsoT,cuBlockFxr%FresoAT, &
                                cuIsoData%siga,cuBlockFxr%wtTempT,cuBlockFxr%nisoT,cuBlockFxr%itTempT,&
                                nAICT,(rgb-1),nrg,MLGLib%f_nmaclv1G,jz)

IF (FxrTri%nResTri.GT.0) THEN
  DEALLOCATE(siglpT_dev)
END IF

END SUBROUTINE

END MODULE

MODULE AuxilDriver
USE TYPEDEF,        ONLY : XsMac_Type
USE TYPEDEF_COMMON, ONLY : CoreXsMac_Type
USE OMP_LIB
IMPLICIT NONE

LOGICAL, PRIVATE :: lInitMacXs = .FALSE.

TYPE(XsMac_Type), POINTER, PRIVATE :: XsMac(:,:)
TYPE(CoreXsMac_Type) :: CoreXsMacHost
INTEGER, ALLOCATABLE :: cuInScatRange(:,:), cuInScatIdx(:,:)

CONTAINS

SUBROUTINE InitCoreMacXS(CoreInfo)
USE TYPEDEF,        ONLY : CoreInfo_Type, cell_type
USE PE_MOD,         ONLY : PE
USE CNTL,           ONLY : nTracerCntl
USE CORE_MOD,       ONLY : GroupInfo
USE XsUtil_mod,     ONLY : AllocXsMac
USE OMP_LIB
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
INTEGER :: ng, ngs, nFxr
INTEGER :: iz, ig, igs, tid
INTEGER :: myzb, myze

INTEGER :: ifxr, icel, ncelltype, nfxrmax
TYPE(cell_type), POINTER :: Cell(:)

ng = GroupInfo%ng
nFxr = CoreInfo%nCoreFxr
myzb = PE%myzb
myze = PE%myze

ALLOCATE(cuInScatRange(2, ng), cuInScatIdx(ng, ng))

cuInScatRange = GroupInfo%InScatRange
cuInScatIdx = 0

ngs = 0;
DO ig = 1, ng
  DO igs = cuInScatRange(1, ig), cuInScatRange(2, ig)
    ngs = ngs + 1
    cuInScatIdx(igs, ig) = ngs
  ENDDO
ENDDO

Cell => CoreInfo%CellInfo
ncelltype = CoreInfo%nCellType

nfxrmax = 0
DO icel = 1, nCellType
  nFxrMax = max(nFxrMax, Cell(icel)%nFxr)
ENDDO
nFxrMax = nFxrMax + 3

CALL omp_set_num_threads(PE%nCMFDThread)

ALLOCATE(XsMac(nFxrMax, PE%nCMFDThread))

!$OMP PARALLEL PRIVATE(tid)
tid = omp_get_thread_num() + 1
DO ifxr = 1, nFxrMax
  XsMac(ifxr, tid)%ng = ng
  CALL AllocXsMac(XsMac(ifxr, tid))
ENDDO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE AllocCoreMacXs(CoreInfo)
USE TYPEDEF,        ONLY : CoreInfo_Type
USE PE_MOD,         ONLY : PE
USE CNTL,           ONLY : nTracerCntl
USE CORE_MOD,       ONLY : GroupInfo
USE XsUtil_mod,     ONLY : AllocXsMac
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
INTEGER :: ng, ngs, nFxr
INTEGER :: iz, ig, igs, tid
INTEGER :: myzb, myze

ng = GroupInfo%ng
nFxr = CoreInfo%nCoreFxr
myzb = PE%myzb
myze = PE%myze
ngs = cuInScatIdx(cuInScatRange(2,ng),ng)
nFxr = nFxr*(myze-myzb+1)

!ALLOCATE(cuCoreXsMac(myzb : myze))

!DO iz = myzb, myze
  ALLOCATE(CoreXsMacHost%XSt(ng, nFxr))
  ALLOCATE(CoreXsMacHost%XStr(ng, nFxr))
  ALLOCATE(CoreXsMacHost%XSa(ng, nFxr))
  ALLOCATE(CoreXsMacHost%XSnf(ng, nFxr))
  ALLOCATE(CoreXsMacHost%XSkf(ng, nFxr))
  ALLOCATE(CoreXsMacHost%XSsm(ngs, nFxr))
  IF (nTracerCntl%ScatOd .GE. 1) ALLOCATE(CoreXsMacHost%XSsmP1(ngs, nFxr))
  IF (nTracerCntl%ScatOd .GE. 2) ALLOCATE(CoreXsMacHost%XSsmP2(ngs, nFxr))
  IF (nTracerCntl%ScatOd .EQ. 3) ALLOCATE(CoreXsMacHost%XSsmP3(ngs, nFxr))
!ENDDO

END SUBROUTINE

SUBROUTINE FreeCoreMacXS
USE CNTL,           ONLY : nTracerCntl
USE PE_MOD,         ONLY : PE
USE XsUtil_mod,     ONLY : FreeXsMac
IMPLICIT NONE
INTEGER :: i
!DO i = PE%myzb, PE%myze
  DEALLOCATE(CoreXsMacHost%XSt,CoreXsMacHost%XStr,CoreXsMacHost%XSa,CoreXsMacHost%XSnf)
  DEALLOCATE(CoreXsMacHost%XSkf,CoreXsMacHost%XSsm)
  IF (nTracerCntl%ScatOd .GE. 1) DEALLOCATE(CoreXsMacHost%XSsmP1)
  IF (nTracerCntl%ScatOd .GE. 2) DEALLOCATE(CoreXsMacHost%XSsmP2)
  IF (nTracerCntl%ScatOd .EQ. 3) DEALLOCATE(CoreXsMacHost%XSsmP3)
!END DO
!DEALLOCATE(cuCoreXsMac)
END SUBROUTINE

SUBROUTINE HomogenizeXS_cuMAC(CoreInfo, superPin, Fxr, PinXS, phis, ng, nxy, myzb, myze, lxslib, lscat1, lsigt)
USE PARAM
USE PE_Mod,         ONLY : PE
USE TYPEDEF,        ONLY : CoreInfo_Type,       FxrInfo_Type,       PinXS_Type,                                     &
                           Pin_Type,            Cell_Type,          XsMac_Type
USE CORE_MOD,       ONLY : GroupInfo
USE CMFD_COMMON
USE BenchXs,        ONLY : xsbaseBen
USE CNTL,           ONLY : nTracerCntl
USE timer,          ONLY : TimeChk, nTracer_dclock
USE MPIComm_Mod,    ONLY : MPI_MAX_REAL
USE CUDA_MASTER
USE CUDA_XS,        ONLY : ConMac,  MacAll, cuCoreMacXs, &
                            AllocCuCoreMacXs, ConstructCuBlockFxr, SetupCuMacXS, &
                            DestroyCuBlockFxr, DestroyCuCoreMacXs, SetupCuSrc, &
                            SetupCuPsiNChi, &
                            BufCon, BufErs, BufHold, BufCnE
IMPLICIT NONE

TYPE(CoreInfo_Type) :: CoreInfo
TYPE(superPin_Type), POINTER :: superPin(:)
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
REAL, POINTER :: phis(:, :, :)
INTEGER :: ng, nxy, myzb, myze
LOGICAL :: lxslib, lscat1, lsigt

TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
INTEGER :: i, j, ig, igf, igs, ixy, iz, ipin, icel, itype, ifxr, ifxr_global, tid
INTEGER :: FxrIdxSt, nLocalFxr, nChi, nfxr, ngs
INTEGER :: ierr, isync

INTEGER :: ngBlock, gb, ge, igb, mg, mgs, gbs, ges
TYPE BlockCoreXS_TYPE
  REAL(GPU_XS_PRECISION), PINNED, ALLOCATABLE, DIMENSION(:,:) :: XSt, XStr, XSa, XSnf, XSkf, XSsm
END TYPE
TYPE(BlockCoreXS_TYPE), ALLOCATABLE :: BlockCoreXS(:)
!REAL(GPU_XS_PRECISION), ALLOCATABLE, DIMENSION(:,:) :: XSt, XStr, XSa, XSnf, XSkf, XSsm

REAL :: tbo, teo

Pin => CoreInfo%Pin
Cell => CoreInfo%CellInfo
nChi = GroupInfo%nChi

nfxr = CoreInfo%nCoreFXR

IF (.NOT. lInitMacXs) THEN
  CALL InitCoreMacXS(CoreInfo)
  lInitMacXs = .TRUE.
END IF

!$ call omp_set_num_threads(PE%nThread)

tbo = nTracer_dclock(.FALSE.,.FALSE.)

!ierr = cudaMemGetInfo(FB,TB)
!print*, (TB-FB)/1024/1024
!print*, __FILE__, __LINE__
CALL ConstructCuBlockFxr(ConMac,0)
CALL AllocCoreMacXs(CoreInfo)
ngs = cuInScatIdx(cuInScatRange(2,ng),ng)
DO iz = myzb, myze
  IF (cuCntl%lGrpBlock) THEN
    ngBlock = (ng-1)/XS_BLOCK_SIZE+1
    ALLOCATE(BlockCoreXS(ngBlock))
    DO igb = 1, ngBlock
      ALLOCATE(BlockCoreXS(igb)%XSt(XS_BLOCK_SIZE,nfxr))
      ALLOCATE(BlockCoreXS(igb)%XStr(XS_BLOCK_SIZE,nfxr))
      ALLOCATE(BlockCoreXS(igb)%XSa(XS_BLOCK_SIZE,nfxr))
      ALLOCATE(BlockCoreXS(igb)%XSnf(XS_BLOCK_SIZE,nfxr))
      ALLOCATE(BlockCoreXS(igb)%XSkf(XS_BLOCK_SIZE,nfxr))
      ALLOCATE(BlockCoreXS(igb)%XSsm(XS_BLOCK_SIZE*ng,nfxr))
    END DO
    CALL AllocCuCoreMacXs(XS_BLOCK_SIZE,nTracerCntl%ScatOd,MacAll)
    DO igb = 1, ngBlock
      gb = (igb-1)*XS_BLOCK_SIZE+1
      ge = min(igb*XS_BLOCK_SIZE,ng)
      mg = XS_BLOCK_SIZE
      gbs = cuInScatIdx(cuInScatRange(1,gb),gb)
      ges = cuInScatIdx(cuInScatRange(2,ge),ge)
      mgs = mg*ng
!      ALLOCATE(XSt(mg,nfxr))
!      ALLOCATE(XStr(mg,nfxr))
!      ALLOCATE(XSa(mg,nfxr))
!      ALLOCATE(XSnf(mg,nfxr))
!      ALLOCATE(XSkf(mg,nfxr))
!      ALLOCATE(XSsm(mgs,nfxr))
!      CALL AllocCuCoreMacXs(mg, nTracerCntl%ScatOd, MacAll)
      CALL SetupCuMacXS(nTracerCntl%ScatOd, .FALSE., iz, gb, ge, MacAll)

      ierr = cudaMemcpyAsync(BlockCoreXS(igb)%XSt, cuCoreMacXs%XSt, mg*nfxr, cudaMemcpyDeviceToHost, cuDevice%myStream)
      ierr = cudaMemcpyAsync(BlockCoreXS(igb)%XStr, cuCoreMacXs%XStr, mg*nfxr, cudaMemcpyDeviceToHost, cuDevice%myStream)
      ierr = cudaMemcpyAsync(BlockCoreXS(igb)%XSa, cuCoreMacXs%XSa, mg*nfxr, cudaMemcpyDeviceToHost, cuDevice%myStream)
      ierr = cudaMemcpyAsync(BlockCoreXS(igb)%XSnf, cuCoreMacXs%XSnf, mg*nfxr, cudaMemcpyDeviceToHost, cuDevice%myStream)
      ierr = cudaMemcpyAsync(BlockCoreXS(igb)%XSkf, cuCoreMacXs%XSkf, mg*nfxr, cudaMemcpyDeviceToHost, cuDevice%myStream)
      ierr = cudaMemcpyAsync(BlockCoreXS(igb)%XSsm, cuCoreMacXs%XSsm, mgs*nfxr, cudaMemcpyDeviceToHost, cuDevice%myStream)
    END DO
    DO igb = 1, ngBlock
      gb = (igb-1)*XS_BLOCK_SIZE+1
      ge = min(igb*XS_BLOCK_SIZE,ng)
      mg = ge-gb+1
      gbs = cuInScatIdx(cuInScatRange(1,gb),gb)
      ges = cuInScatIdx(cuInScatRange(2,ge),ge)
      mgs = ges-gbs+1
      IF (igb.EQ.1) ierr = cudaStreamSynchronize(cuDevice%myStream)
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      DO i = 1, nfxr
      CoreXsMacHost%XSt(gb:ge, (iz-myzb)*nfxr+i) = BlockCoreXS(igb)%XSt(1:mg, i)
      CoreXsMacHost%XStr(gb:ge, (iz-myzb)*nfxr+i) = BlockCoreXS(igb)%XStr(1:mg, i)
      CoreXsMacHost%XSa(gb:ge, (iz-myzb)*nfxr+i) = BlockCoreXS(igb)%XSa(1:mg, i)
      CoreXsMacHost%XSnf(gb:ge, (iz-myzb)*nfxr+i) = BlockCoreXS(igb)%XSnf(1:mg, i)
      CoreXsMacHost%XSkf(gb:ge, (iz-myzb)*nfxr+i) = BlockCoreXS(igb)%XSkf(1:mg, i)
      CoreXsMacHost%XSsm(gbs:ges,(iz-myzb)*nfxr+i) = BlockCoreXS(igb)%XSsm(1:mgs, i)
      END DO
      !$OMP END PARALLEL DO

      DEALLOCATE(BlockCoreXS(igb)%XSt,BlockCoreXS(igb)%XStr,BlockCoreXS(igb)%XSa)
      DEALLOCATE(BlockCoreXS(igb)%XSnf,BlockCoreXS(igb)%XSkf,BlockCoreXS(igb)%XSsm)
    END DO
    DEALLOCATE(BlockCoreXS)
    CALL DestroyCuCoreMacXs(MacAll, nTracerCntl%ScatOd)
  ELSE
    CALL AllocCuCoreMacXs(ng, nTracerCntl%ScatOd, MacAll)
    CALL SetupCuMacXS(nTracerCntl%ScatOd, .FALSE., iz, 1, ng, MacAll)
    ierr = cudaMemcpy(CoreXsMacHost%XSt(1:,(iz-myzb)*nfxr+1:), cuCoreMacXs%XSt, ng*nfxr, cudaMemcpyDeviceToHost)
    ierr = cudaMemcpy(CoreXsMacHost%XStr(1:,(iz-myzb)*nfxr+1:), cuCoreMacXs%XStr, ng*nfxr, cudaMemcpyDeviceToHost)
    ierr = cudaMemcpy(CoreXsMacHost%XSa(1:,(iz-myzb)*nfxr+1:), cuCoreMacXs%XSa, ng*nfxr, cudaMemcpyDeviceToHost)
    ierr = cudaMemcpy(CoreXsMacHost%XSnf(1:,(iz-myzb)*nfxr+1:), cuCoreMacXs%XSnf, ng*nfxr, cudaMemcpyDeviceToHost)
    ierr = cudaMemcpy(CoreXsMacHost%XSkf(1:,(iz-myzb)*nfxr+1:), cuCoreMacXs%XSkf, ng*nfxr, cudaMemcpyDeviceToHost)
    ierr = cudaMemcpy(CoreXsMacHost%XSsm(1:,(iz-myzb)*nfxr+1:), cuCoreMacXs%XSsm, ngs*nfxr, cudaMemcpyDeviceToHost)
    CALL DestroyCuCoreMacXs(MacAll,nTracerCntl%ScatOd)
  END IF
END DO

teo = nTracer_dclock(.FALSE.,.FALSE.)
TimeChk%cuHomDevTime = TimeChk%cuHomDevTime+(teo-tbo)
CALL MPI_MAX_REAL(TimeChk%cuHomDevTime,PE%MPI_NTRACER_COMM,.TRUE.)
tbo = teo;

!$ call omp_set_num_threads(PE%nCMFDThread)

!$OMP PARALLEL PRIVATE(tid, ipin, icel, igs, itype, ifxr, ifxr_global, FxrIdxSt, nLocalFxr)
tid = omp_get_thread_num() + 1

DO iz = myzb, myze
  !$OMP DO SCHEDULE(GUIDED)
  DO ixy = 1, nxy
    nLocalFxr = 0
    DO j = 1, superPin(ixy)%nxy
      ipin = superPin(ixy)%pin(j)
      icel = Pin(ipin)%Cell(iz)
      FxrIdxSt = Pin(ipin)%FxrIdxSt
      DO i = 1, Cell(icel)%nFxr
        ifxr = nLocalFxr + i; ifxr_global = FxrIdxSt + i - 1 +(iz-myzb)*nfxr
        IF (lxslib) THEN
          XsMac(ifxr, tid)%XsMacT = CoreXsMacHost%XSt(:, ifxr_global)
          XsMac(ifxr, tid)%XsMacTr = CoreXsMacHost%XStr(:, ifxr_global)
          XsMac(ifxr, tid)%XsMacA = CoreXsMacHost%XSa(:, ifxr_global)
          XsMac(ifxr, tid)%XsMacNf = CoreXsMacHost%XSnf(:, ifxr_global)
          XsMac(ifxr, tid)%XsMacKf = CoreXsMacHost%XSkf(:, ifxr_global)
          DO ig = 1, ng
            DO igf = cuInScatRange(1, ig), cuInScatRange(2, ig)
              igs = cuInScatIdx(igf, ig)
              XsMac(ifxr, tid)%XsMacSm(igf, ig) = CoreXsMacHost%XSsm(igs, ifxr_global)
            ENDDO
          ENDDO
          XsMac(ifxr, tid)%Chi = 0.0
          XsMac(ifxr, tid)%lFuel = FALSE
          IF (Fxr(ifxr_global, iz)%lDepl) THEN
            XsMac(ifxr, tid)%Chi(1 : nChi) = Fxr(ifxr_global, iz)%Chi
            XsMac(ifxr, tid)%lFuel = TRUE
          ENDIF
        ELSE
          itype = Fxr(ifxr_global, iz)%imix
          CALL xsbaseBen(itype, 1, ng, 1, ng, lscat1, XsMac(ifxr, tid))
        ENDIF
      ENDDO
      nLocalFxr = nLocalFxr + Cell(icel)%nFxr
    ENDDO
    CALL HomogenizeCellXS(CoreInfo, superPin(ixy), PinXS(ixy, iz), XsMac(1 : nLocalFxr, tid), phis, iz, ng, lxslib, lsigt)
  ENDDO
  !$OMP END DO
ENDDO
!$OMP END PARALLEL

teo = nTracer_dclock(.FALSE.,.FALSE.)
TimeChk%cuHomHostTime = TimeChk%cuHomHostTime+(teo-tbo)
CALL MPI_MAX_REAL(TimeChk%cuHomHostTime,PE%MPI_NTRACER_COMM,.TRUE.)

!print*, __FILE__, __LINE__
CALL DestroyCuBlockFxr(ConMac,0)
CALL FreeCoreMacXS
!print*, __FILE__, __LINE__

END SUBROUTINE

SUBROUTINE CUDASubGrpEffXSGen(Core, Fxr, nTracerCntl, PE)
  USE PARAM
  USE MPIComm_mod,      ONLY : MPI_SYNC, MPI_MAX_REAL
  USE TYPEDEF,      ONLY : CoreInfo_Type, FxrInfo_Type,PE_Type
  USE CNTL,         ONLY : nTracerCntl_Type
!  USE XSLIB_MOD,    ONLY : mapnucl
  USE SUBGRP_MOD,   ONLY : UpdtCoreIsoInfo
  USE IOUTIL,       ONLY : message
  USE files,        ONLY : io8
  USE XS_COMMON,    ONLY : UpdateResPinMapNucl!, FxrD
  USE CUDA_MASTER,  ONLY : cuCntl
  USE CUDA_XS,      ONLY : EffMacGen_cuEff, ConEff, ConstructCuBlockFxr, ConstructCuResPin, &
                          DestroyCuBlockFxr, DestroyCuResPin, EffMacIntegrated_cuEff,&
                          EffMacIntegrated_RIP, EffMacIntegrated_CLDnAICT,&
                          EffMacIntegrated_GenD_RB, EffMacIntegrated_GenD_GB,&
                          EffRIP, EffFxrD, EffFxrT, EffClr, EffAll, EffBlkD
  USE Timer,        ONLY : nTracer_dclock, TimeChk
  !USE CUDAFOR
  IMPLICIT NONE
  TYPE(CoreInfo_Type) :: Core
  TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
  TYPE(nTracerCntl_Type) :: nTracerCntl
  TYPE(PE_Type) :: PE

  INTEGER :: iz
  REAL :: tbo, teo

  INTEGER(8) :: totalByte, freeByte
  INTEGER :: ierr

  LOGICAL, PARAMETER :: lRegionDcmp = .TRUE.
!  LOGICAL, PARAMETER :: lRegionDcmp = .FALSE.

  CALL MPI_SYNC(PE%MPI_NTRACER_COMM)

  IF (.not.nTracerCntl%lrestrmt) RETURN

  WRITE(mesg, '(A)') 'Update Subgroup Effective XSs...'
  IF(PE%master) CALL message(io8, .TRUE., .TRUE., mesg)

  tbo = nTracer_dclock(FALSE,FALSE)

!  print*, __FILE__, __LINE__, PE%myzb
  CALL UpdtCoreIsoInfo(Core, Fxr, PE)
!  print*, __FILE__, __LINE__, PE%myzb
  CALL UpdateResPinMapNucl
!  print*, __FILE__, __LINE__, PE%myzb
!  CALL ConstructCuBlockFxr(ConEff,EffAll)
!  print*, __FILE__, __LINE__, PE%myzb
!  CALL ConstructCuResPin(ConEff)
!  print*, __FILE__, __LINE__, PE%myzb
  IF (cuCntl%lGrpBlock) THEN
    DO iz = PE%myzb,PE%myze
      IF (Core%lFuelPlane(iz)) THEN
        CALL ConstructCuResPin(ConEff)
        CALL EffMacIntegrated_RIP(iz)
        CALL DestroyCuResPin(ConEff,EffRIP)
        CALL ConstructCuBlockFxr(ConEff,EffBlkD)
        IF (lRegionDcmp) THEN
          CALL EffMacIntegrated_GenD_RB(iz)
        ELSE
          CALL EffMacIntegrated_GenD_GB(iz)
        END IF
        CALL DestroyCuBlockFxr(ConEff,EffBlkD)
        CALL ConstructCuBlockFxr(ConEff,EffFxrT)
        CALL EffMacIntegrated_CLDnAICT(iz)
        CALL DestroyCuBlockFxr(ConEff,EffFxrT)
        CALL DestroyCuResPin(ConEff,EffClr)
      END IF
    END DO
  ELSE
    DO iz = PE%myzb,PE%myze
      IF (Core%lFuelPlane(iz)) THEN
        CALL ConstructCuBlockFxr(coneff,EffAll)
        CALL ConstructCuResPin(ConEff)
        CALL EffMacIntegrated_cuEff(iz)
        CALL DestroyCuBlockFxr(ConEff,EffAll)
        CALL DestroyCuResPin(ConEff,EffAll)
      END IF
    END DO
  END IF

!  print*, __FILE__, __LINE__, PE%myzb

!  CALL DestroyCuBlockFxr(ConEff,EffAll)
!  print*, __FILE__, __LINE__, PE%myzb
!  CALL DestroyCuResPin(ConEff,EffAll)
!  print*, __FILE__, __LINE__, PE%myzb

  CALL MPI_SYNC(PE%MPI_NTRACER_COMM)

  teo = nTracer_dclock(FALSE,FALSE)
  TimeChk%XSsubTime = TimeChk%XSsubTime+(teo-tbo)

  CALL MPI_MAX_REAL(TimeChk%cuXSPreTime, PE%MPI_NTRACER_COMM, .TRUE.)
  CALL MPI_MAX_REAL(TimeChk%cuXSMainTime, PE%MPI_NTRACER_COMM, .TRUE.)
END SUBROUTINE

END MODULE
#endif

