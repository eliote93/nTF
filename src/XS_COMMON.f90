MODULE XS_COMMON
USE TYPEDEF_COMMON, ONLY : BlockFxrDepl_Type, BlockFxrTri_Type, BlockIsodata_Type, &
                        BlockMLGLib_Type, BlockResVarPin_Type, BlockResIsoData_Type, &
                        GlobalPinInfo_Type
IMPLICIT NONE

LOGICAL :: lrestrmt, lED, lMLG, lCAT, lRST, lISO, lInitBoron = .FALSE.
INTEGER, PARAMETER :: THFeed = 1, PPMFeed = 2, DEPLFeed = 3, XeDynFeed = 4

TYPE(BlockFxrDepl_Type) :: FxrD
TYPE(BlockFxrTri_Type) :: FxrTri

TYPE(BlockIsodata_Type) :: IsoData
TYPE(BlockMLGLib_Type) :: MLGLib
TYPE(BlockResVarPin_Type) :: ResPin
TYPE(BlockResIsoData_Type) :: ResIsoData

TYPE(GlobalPinInfo_Type) :: GlobalPinInfo

INTEGER :: ibfxrD(32), ibfxrT(32)
LOGICAL :: lFPlane(64)

CONTAINS

SUBROUTINE AlignIsodata(GroupInfo)
USE XSLIB_MOD,  ONLY : ldiso
USE TYPEDEF,    ONLY : GroupInfo_Type
IMPLICIT NONE
TYPE(GroupInfo_Type) :: GroupInfo
INTEGER :: ntiso, ng
INTEGER :: nttemp, nttempPx, Np0, Np1, Np2, Np3
INTEGER :: ntemp, np1temp, ib, ie, rng
INTEGER :: iso, it, ig, jg

REAL, POINTER :: ptsTemp(:), ptsTempPx(:)
INTEGER, POINTER :: ptrTemp(:), ptrTempPx(:)

REAL, DIMENSION(:,:), POINTER :: siga, sigf, sigkf, signf, sigs, sigstr
INTEGER, DIMENSION(:,:), POINTER :: gidSmp0, gidSmp1, gidSmp2, gidSmp3
INTEGER, DIMENSION(:,:), POINTER :: ptrSmp0, ptrSmp1, ptrSmp2, ptrSmp3
REAL, POINTER :: smp0(:), smp1(:), smp2(:), smp3(:)

INTEGER :: nBytes

ntiso = GroupInfo%ntiso; ng = GroupInfo%ng
ALLOCATE(ptrTemp(ntiso+1), ptrTempPx(ntiso+1))

nBytes = ntiso+ntiso+2

nttemp = 0; nttempPx = 0;
Np0 = 0; Np1 = 0; Np2 = 0; Np3 = 0;
DO iso = 1, ntiso
  ptrTemp(iso) = nttemp; ptrTempPx(iso) = nttempPx;
  ntemp = ldiso(iso)%ntemp; np1temp = ldiso(iso)%np1temp
  nttemp = nttemp + ntemp
  nttempPx = nttempPx + np1temp
  DO it = 1, ntemp
    DO ig = 1, ng
      ib = ldiso(iso)%sm(ig,it)%ib; ie = ldiso(iso)%sm(ig,it)%ie
      Np0 = Np0+ie-ib+1
    END DO
  END DO
  DO it = 1, np1temp
    DO ig = 1, ng
      ib = ldiso(iso)%smp1(ig,it)%ib; ie = ldiso(iso)%smp1(ig,it)%ie
      Np1 = Np1+ie-ib+1
      ib = ldiso(iso)%smp2(ig,it)%ib; ie = ldiso(iso)%smp2(ig,it)%ie
      Np2 = Np2+ie-ib+1
      ib = ldiso(iso)%smp3(ig,it)%ib; ie = ldiso(iso)%smp3(ig,it)%ie
      Np3 = Np3+ie-ib+1
    END DO
  END DO
END DO
ptrTemp(iso) = nttemp; ptrTempPx(iso) = nttempPx;

ALLOCATE(ptsTemp(nttemp), ptsTempPx(nttempPx));
ALLOCATE(siga(ng,nttemp),sigkf(ng,nttemp),sigf(ng,nttemp),signf(ng,nttemp),sigs(ng,nttemp),sigstr(ng,nttemp))

ALLOCATE(gidSmp0(2*ng,nttemp),gidSmp1(2*ng,nttempPx),gidSmp2(2*ng,nttempPx),gidSmp3(2*ng,nttempPx))
ALLOCATE(ptrSmp0(ng,nttemp),ptrSmp1(ng,nttempPx),ptrSmp2(ng,nttempPx),ptrSmp3(ng,nttempPx))
ALLOCATE(smp0(Np0), smp1(Np1), smp2(Np2), smp3(Np3))

!nBytes=nBytes+ntiso+ng*nttemp*5+2*ng*(nttemp+nttempPx*3)+ng*(nttemp+nttempPx*3)+Np0+Np1+Np2+Np3
!print*, nBytes*8

nttemp = 0; nttempPx = 0;
Np0 = 0; Np1 = 0; Np2 = 0; Np3 = 0;
DO iso = 1, ntiso
  ntemp = ldiso(iso)%ntemp
  np1temp = ldiso(iso)%np1temp
  siga(1:ng,nttemp+1:nttemp+ntemp) = ldiso(iso)%siga
  if (ldiso(iso)%ifis.ne.0) then
    sigf(1:ng,nttemp+1:nttemp+ntemp) = ldiso(iso)%sigf
    sigkf(1:ng,nttemp+1:nttemp+ntemp) = ldiso(iso)%sigf*ldiso(iso)%kappa
    signf(1:ng,nttemp+1:nttemp+ntemp) = ldiso(iso)%signf
  else
    sigf(1:ng,nttemp+1:nttemp+ntemp) = 0.
    sigkf(1:ng,nttemp+1:nttemp+ntemp) = 0.
    signf(1:ng,nttemp+1:nttemp+ntemp) = 0.
  endif
  sigs(1:ng,nttemp+1:nttemp+ntemp) = ldiso(iso)%sigs
  sigstr(1:ng,nttemp+1:nttemp+ntemp) = ldiso(iso)%sigstr
  DO it = 1, ntemp
    nttemp = nttemp+1
    ptsTemp(nttemp) = ldiso(iso)%temp(it)
    DO ig = 1, ng
      ptrSmp0(ig,nttemp) = Np0+1
      ib = ldiso(iso)%sm(ig,it)%ib; ie = ldiso(iso)%sm(ig,it)%ie
      rng = ie-ib+1
      gidSmp0(ig*2-1,nttemp) = ib; gidSmp0(ig*2,nttemp) = ie
      DO jg = ib, ie
        Np0 = Np0+1
        smp0(Np0) = ldiso(iso)%sm(ig,it)%from(jg)
      END DO
    END DO
  END DO
  DO it = 1, np1temp
    nttempPx = nttempPx+1
    ptsTempPx(nttempPx) = ldiso(iso)%p1temp(it)
    DO ig = 1, ng
      ptrSmp1(ig,nttempPx) = Np1+1
      ptrSmp2(ig,nttempPx) = Np2+1
      ptrSmp3(ig,nttempPx) = Np3+1
      ib = ldiso(iso)%smp1(ig,it)%ib; ie = ldiso(iso)%smp1(ig,it)%ie
      rng = ie-ib+1
      gidSmp1(ig*2-1,nttempPx) = ib; gidSmp1(ig*2,nttempPx) = ie
      DO jg = ib,ie
        Np1 = Np1+1
        smp1(Np1) = ldiso(iso)%smp1(ig,it)%from(jg)
      END DO

      ib = ldiso(iso)%smp2(ig,it)%ib; ie = ldiso(iso)%smp2(ig,it)%ie
      rng = ie-ib+1
      gidSmp2(ig*2-1,nttempPx) = ib; gidSmp2(ig*2,nttempPx) = ie
      DO jg = ib,ie
        Np2 = Np2+1
        smp2(Np2) = ldiso(iso)%smp2(ig,it)%from(jg)
      END DO

      ib = ldiso(iso)%smp3(ig,it)%ib; ie = ldiso(iso)%smp3(ig,it)%ie
      rng = ie-ib+1
      gidSmp3(ig*2-1,nttempPx) = ib; gidSmp3(ig*2,nttempPx) = ie
      DO jg = ib,ie
        Np3 = Np3+1
        smp3(Np3) = ldiso(iso)%smp3(ig,it)%from(jg)
      END DO
    END DO
  END DO
END DO
IsoData%ntiso = ntiso; IsoData%ng = ng;
IsoData%nttemp = nttemp; IsoData%nttempPx = nttempPx;
IsoData%Np0 = Np0; IsoData%Np1 = Np1; IsoData%Np2 = Np2; IsoData%Np3 = Np3;
IsoData%ptrTemp => ptrTemp; IsoData%ptrTempPx => ptrTempPx
IsoData%ptsTemp => ptsTemp; IsoData%ptsTempPx => ptsTempPx
IsoData%siga => siga; IsoData%sigf => sigf;IsoData%sigs => sigs;
IsoData%sigkf => sigkf; IsoData%signf => signf; IsoData%sigstr => sigstr;
IsoData%gidSmp0 => gidSmp0; IsoData%gidSmp1 => gidSmp1;
IsoData%gidSmp2 => gidSmp2; IsoData%gidSmp3 => gidSmp3;
IsoData%ptrSmp0 => ptrSmp0; IsoData%ptrSmp1 => ptrSmp1;
IsoData%ptrSmp2 => ptrSmp2; IsoData%ptrSmp3 => ptrSmp3;
IsoData%smp0 => smp0; IsoData%smp1 => smp1;
IsoData%smp2 => smp2; IsoData%smp3 => smp3;

END SUbROUTINE

SUBROUTINE AlignResIsoData(GroupInfo)
USE XSLIB_MOD,  ONLY : ldiso, mapnucl, nreshel
USE TYPEDEF,    ONLY : GroupInfo_Type
IMPLICIT NONE
TYPE(GroupInfo_Type) :: GroupInfo
INTEGER :: irgb, irge
INTEGER :: ntsig0, ntrtemp, ntlv, ntrif, ntrifRat, NTempPot, NTempLv, nresiso, ntrifRx
INTEGER :: nsig0, nrtemp, nlv, nrif, nrat
INTEGER :: i,j,k,l,ig,iT,ntiso
IF (.NOT. lrestrmt) RETURN

irgb = GroupInfo%nofg+1; irge = GroupInfo%nofg+GroupInfo%norg
ntiso = IsoData%ntiso
ALLOCATE(ResIsoData%IsoG2R(ntiso))

ntsig0 = 0; ntrtemp = 0; ntlv = 0; nresiso = 0
NTempPot = 0; NTempLv = 0
ntrifRat = 0; ntrifRx = 0; ntrif = 0;
ResIsoData%IsoG2R = 0
DO i = 1, ntiso
  IF (ldiso(i)%lreso) THEN
    j = ldiso(i)%nid; k = mapnucl(j-500)
    IF (mod(j,1000).GT.500 .AND. k.NE.0) CYCLE
    nresiso = nresiso + 1
    ResIsoData%IsoG2R(i) = nresiso
    nsig0 = ldiso(i)%nsig0; nrtemp = ldiso(i)%nrtemp; nlv = ldiso(i)%nlv
    ntsig0 = ntsig0+nsig0; ntrtemp = ntrtemp+nrtemp; ntlv = ntlv+nlv
    NTempPot = NTempPot+nsig0*nrtemp
    NTempLv = NTempLv+nlv*nrtemp

    nrif = ldiso(i)%nrif
    ntrif = ntrif+nrif
    DO j = 1, nrif
      nrat = ldiso(i)%rif(j)%nrat
      ntrifRat = ntrifRat+nrat
      ntrifRx = ntrifRx+nrat*nsig0*nrtemp
    END DO
  END IF
END DO
DO i = 1, ntiso
  IF (ldiso(i)%lreso) THEN
    j = ldiso(i)%nid; k = mapnucl(j-500)
    IF (mod(j,1000).GT.500 .AND. k.NE.0) THEN
      ResIsoData%IsoG2R(i) = ResIsoData%IsoG2R(k)
!      print*, j, mapnucl(j-500), ResIsoData%IsoG2R(i)
    END IF
  END IF
END DO
ResIsoData%ntrtemp = ntrtemp; ResIsoData%ntsig0 = ntsig0; ResIsoData%ntlv = ntlv
ResIsoData%nresiso = nresiso; ResIsoData%ntempPot = NTempPot; ResIsoData%ntempLv = NTempLv
ResIsoData%ntrif = ntrif; ResIsoData%ntrifRat = ntrifRat; ResIsoData%ntRifRx = ntrifRx
ALLOCATE(ResIsoData%lclad(ntiso), ResIsoData%ityp(nresiso), ResIsoData%sigp(ntiso))
ALLOCATE(ResIsoData%ptrRTemp(nresiso+1),ResIsoData%ptrSig0(nresiso+1),ResIsoData%ptrNlv(nresiso+1))
ALLOCATE(ResIsoData%ptrTempPot(nresiso+1),ResIsoData%ptrTempLv(nresiso+1))
ALLOCATE(ResIsoData%lamsigp(irgb:irge,ntiso))
ALLOCATE(ResIsoData%sig0sq(ntsig0), ResIsoData%rtempsq(ntrtemp))
ALLOCATE(ResIsoData%xsalog(NTempPot,irgb:irge),ResIsoData%ri_a(NTempPot,irgb:irge))
ALLOCATE(ResIsoData%lvabs(ntlv,irgb:irge),ResIsoData%lvabs_log(ntlv,irgb:irge))
ALLOCATE(ResIsoData%lvfis(ntlv,irgb:irge))
ALLOCATE(ResIsoData%wgtabs(NTempLv,irgb:irge),ResIsoData%wgtfis(NTempLv,irgb:irge))

ALLOCATE(ResIsoData%ptrRifRat(ntrif+1),ResIsoData%ptrRifRx(ntrifRat+1))
ALLOCATE(ResIsoData%ratlog(ntrifRat),ResIsoData%rifabs(ntrifRx,irgb:irge),ResIsoData%riffis(ntrifRx,irgb:irge))
ALLOCATE(ResIsoData%rifid(ntiso,nresiso))
ResIsoData%rifid = 0;

ResIsoData%ptrRTemp(nresiso+1) = ntrtemp
ResIsoData%ptrSig0(nresiso+1) = ntsig0
ResIsoData%ptrNlv(nresiso+1) = ntlv
ResIsoData%ptrTempPot(nresiso+1) = NTempPot
ResIsoData%ptrTempLv(nresiso+1) = NTempLv
ResIsoData%ptrRifRat(ntrif+1) = ntrifRat
ResIsoData%ptrRifRx(ntrifRat+1) = ntrifRx

!print*, nreshel, nresiso

ntsig0 = 0; ntrtemp = 0; ntlv = 0; nresiso = 0
NTempPot = 0; NTempLv = 0
ntrifRat = 0; ntrifRx = 0; ntrif = 0;
DO i = 1, ntiso
  ResIsoData%sigp(i) = ldiso(i)%sigp
  ResIsoData%lamsigp(:,i) = ldiso(i)%lamsigp(:)
  ResIsoData%lclad(i) = ldiso(i)%lclad
  IF (ldiso(i)%lreso) THEN
    j = ldiso(i)%nid; k = mapnucl(j-500)
    IF (mod(j,1000).GT.500 .AND. k.NE.0) CYCLE
    nresiso = nresiso + 1
    ResIsoData%ityp(nresiso) = ldiso(i)%ityp

    ResIsoData%ptrRTemp(nresiso) = ntrtemp
    ResIsoData%ptrSig0(nresiso) = ntsig0
    ResIsoData%ptrNlv(nresiso) = ntlv
    ResIsoData%ptrTempPot(nresiso) = NTempPot
    ResIsoData%ptrTempLv(nresiso) = NTempLv


    nsig0 = ldiso(i)%nsig0; nrtemp = ldiso(i)%nrtemp; nlv = ldiso(i)%nlv
    DO j = 1, nsig0
      ntsig0 = ntsig0+1
      ResIsoData%sig0sq(ntsig0) = ldiso(i)%sig0sq(j)
      DO k = 1, nrtemp
        NTempPot = NTempPot+1
        DO ig = irgb,irge
          ResIsoData%xsalog(NTempPot-k-(j-1)*nrtemp+(k-1)*nsig0+j,ig) = ldiso(i)%xsalog(j,ig,k) ! Transposing, -k+(k-1)*nsig0+j
          ResIsoData%ri_a(NTempPot,ig) = ldiso(i)%ri_a(j,ig,k)
        END DO
      END DO
    END DO
    DO j = 1, nlv
      ntlv = ntlv + 1
      DO ig = irgb,irge
        ResIsoData%lvabs(ntlv,ig) = ldiso(i)%lvabs(j,ig)
        ResIsoData%lvabs_log(ntlv,ig) = ldiso(i)%lvabs_log(j,ig)
        if (ldiso(i)%ityp.NE.3) CYCLE
        ResIsoData%lvfis(ntlv,ig) = ldiso(i)%lvfis(j,ig)
      END DO
      DO k = 1, nrtemp
        NTempLv = NTempLv + 1
        DO ig = irgb,irge
          ResIsoData%wgtabs(NTempLv,ig)=ldiso(i)%wgtabs(j,k,ig)
          if (ldiso(i)%ityp.NE.3) CYCLE
          ResIsoData%wgtfis(NTempLv,ig)=ldiso(i)%wgtfis(j,k,ig)
        END DO
      END DO
    END DO
    DO j = 1, nrtemp
      ntrtemp = ntrtemp+1
      ResIsoData%rtempsq(ntrtemp) = ldiso(i)%rtempsq(j)
    END DO
    DO j = 1, ntiso
      k = ldiso(i)%rifid(j)
      IF (k.GT.0) ResIsoData%rifid(j,nresiso) = ntrif+k
    END DO
    nrif = ldiso(i)%nrif
    DO j = 1, nrif
      ntrif = ntrif+1
      ResIsoData%ptrRifRat(ntrif) = ntrifRat

      nrat = ldiso(i)%rif(j)%nrat
      ResIsoData%ratlog(ntrifRat+1:ntrifRat+nrat) = ldiso(i)%rif(j)%ratlog
      DO k = 1,nrat
        ntrifRat = ntrifRat+1
        ResIsoData%ptrRifRx(ntrifRat) = ntrifRx
        DO ig = irgb,irge
          DO iT = 1, nrtemp
            DO l = 1, nsig0
              ResIsoData%rifabs(ntrifRx+(iT-1)*nsig0+l,ig)=ldiso(i)%rif(j)%abs(k,l,ig,iT)
              IF (ResIsoData%ityp(nresiso).EQ.3) THEN
                ResIsoData%riffis(ntrifRx+(iT-1)*nsig0+l,ig)=ldiso(i)%rif(j)%fis(k,l,ig,iT)
              ELSE
                ResIsoData%riffis(ntrifRx+(iT-1)*nsig0+l,ig)=0
              END IF
            END DO
          END DO
        END DO
        ntrifRx = ntrifRx+nsig0*nrtemp;
      END DO
    END DO
  END IF
END DO

END SUBROUTINE

SUBROUTINE AlignMLGdata()
USE GEOM, ONLY : Core
USE PE_MOD, ONLY : PE
USE XSLIB_MOD, ONLY : mlgdata, mlgdata0
IMPLICIT NONE
INTEGER :: myzb, myze, iz
INTEGER :: nmlgF, nmlgF1G, nmlgC1G
IF (.NOT. lrestrmt .AND. .NOT. lMLG) RETURN

myzb = PE%myzb; myze = PE%myze

nmlgF = mlgdata0%f_nmaclv; nmlgF1G = mlgdata0%f_nmaclv1G
nmlgC1G = mlgdata0%c_nmaclv1G

MLGLib%f_nmaclv = nmlgF; MLGLib%f_nmaclv1G = nmlgF1G;
MLGLib%c_nmaclv1G = nmlgC1G;
ALLOCATE(MLGLib%f_maclv(nmlgF,myzb:myze),MLGLib%f_maclv_log(nmlgF,myzb:myze))
ALLOCATE(MLGLib%f_maclv1G(nmlgF1G,myzb:myze),MLGLib%f_maclv1G_log(nmlgF1G,myzb:myze))
ALLOCATE(MLGLib%c_maclv1G(nmlgC1G),MLGLib%c_maclv1G_log(nmlgC1G))

MLGLib%c_maclv1G = mlgdata0%c_maclv1G
MLGLib%c_maclv1G_log = mlgdata0%c_maclv1G_log
DO iz = myzb, myze
  IF (core%lfuelplane(iz)) THEN
    MLGLib%f_maclv(:,iz) = mlgdata(iz)%f_maclv
    MLGLib%f_maclv_log(:,iz) = mlgdata(iz)%f_maclv_log
  END IF
  IF (core%laicplane(iz)) THEN
    MLGLib%f_maclv1G(:,iz) = mlgdata(iz)%f_maclv1G
    MLGLib%f_maclv1G_log(:,iz) = mlgdata(iz)%f_maclv1G_log
  END IF
END DO

END SUBROUTINE

SUBROUTINE AlignFxrs(Fxr, GroupInfo)
USE TYPEDEF,    ONLY : Fxrinfo_type, GroupInfo_Type
USE CNTL,       ONLY : nTracerCntl
USE GEOM,       ONLY : Core
USE PE_Mod,     ONLY : PE
USE XSLIB_MOD,  ONLY : mlgdata0, nlvflxmax, nreshel, nCat, mapnucl, mapnuclp13, mapnuclRes, mapfis, idnp1hel

IMPLICIT NONE

TYPE(Fxrinfo_type), POINTER :: Fxr(:,:)
TYPE(GroupInfo_Type) :: GroupInfo

INTEGER :: nFxr, nxy, iRgBeg, iRgEnd, myzb, myze, ntiso, nchi, nisomaxTri, nxyz
INTEGER :: nFxrD, nResD, nCldD, nGenD, ntfsrD
INTEGER :: nFxrTri, nResTri, nAICTri, nCldTri, nGenTri, ntfsrTri
INTEGER :: NtTri, NtRTri

INTEGER, DIMENSION(:,:), POINTER :: DMapNucl, DMapNuclp13, DMapFis, DMapNuclRes
LOGICAL, DIMENSION(:), POINTER :: lfuel, lCLD, lres, lGD, lh2o, lvoid

INTEGER, DIMENSION(:), POINTER :: niso, ndim, ipin, fsrst, nfsr, mapfsr
REAL, DIMENSION(:), POINTER :: temp, area, xstilde, h2ofrac

INTEGER, DIMENSION(:,:), POINTER :: idiso
REAL, DIMENSION(:,:), POINTER :: pnum, chi
REAL(4), DIMENSION(:,:), POINTER :: xseq_f_1g, xseq_c_1g, FnAdj, &
  fresoa, fresof, fresos, fresostr, fresonf, fresokf
REAL(4), POINTER :: xseq_f_mg(:,:,:), Ftadj(:,:,:), xseq(:,:,:,:), NDAF(:,:,:,:), &
  fresoAIso(:,:,:), fresoFIso(:,:,:), fresoSIso(:,:,:), fresoSSIso(:,:,:), fresoS1Iso(:,:,:)


INTEGER, DIMENSION(:), POINTER :: TMapNucl, TMapNuclp13, TMapNuclRes
LOGICAL, DIMENSION(:), POINTER :: tlAIC, tlCLD, tlres, tlh2o, tlvoid

INTEGER, DIMENSION(:), POINTER :: tniso, tndim, tipin, tfsrst, tnfsr, tmapfsr
REAL, DIMENSION(:), POINTER :: ttemp, tarea, txstilde, th2ofrac

INTEGER, DIMENSION(:), POINTER :: tidiso, AccNt, AccNtR, tDirectMapnucl
REAL, DIMENSION(:), POINTER :: tpnum
REAL(4), DIMENSION(:,:), POINTER :: txseq_f_1g, txseq_c_1g, tFnAdj, &
  tfresoa, tfresos, tfresostr
REAL(4), DIMENSION(:), POINTER :: tfresoAIso, tfresoSIso, tfresoSSIso, tfresoS1Iso
REAL(4), POINTER :: txseq_f_mg(:,:,:), txseq_c_mg(:,:,:), tFtadj(:,:,:), txseq(:,:,:,:), tNDAF(:,:,:,:)

INTEGER :: iz, ifxr, iso, j, k, ifsr, mdim, RngResG, p1idx, jdxfxr, icel, FsrIdxSt, FxrIdxSt, ig
LOGICAL :: lFpln
TYPE(Fxrinfo_type), POINTER :: myFxr

lED = nTracerCntl%lED; lMLG = nTracerCntl%lMLG; lRST = nTracerCntl%lRST;
lCAT = nTracerCntl%lCAT; lrestrmt = nTracerCntl%lrestrmt
lISO = (.NOT.lMLG).AND.(.NOT.lCAT)
nFxr = Core%nCoreFXR; nxy = Core%nxy
ntiso = GroupInfo%ntiso
iRgBeg = GroupInfo%nofg+1
iRgEnd = GroupInfo%nofg+GroupInfo%norg
RngResG = iRgEnd-iRgBeg+1
nchi = GroupInfo%nchi
myzb = PE%myzb; myze = PE%myze
nxyz = nxy*(myze-myzb+1)

ntfsrD = 0; ntfsrTri = 0;
nFxrD = 0 ; nResD = 0; nCldD = 0; nGenD = 0;
nFxrTri = 0; nResTri = 0; nCldTri = 0; nGenTri = 0; nAICTri = 0;
NtTri = 0; NtRTri = 0;
nisomaxTri = 0;
DO iz = myzb, myze
  ibfxrD(iz-myzb+1) = nFxrD+1
  ibfxrT(iz-myzb+1) = nFxrTri+1
  lfPln = Core%lfuelplane(iz)
  lFPlane(iz) = lFpln
  DO ifxr = 1, nFxr
    IF (Fxr(ifxr,iz)%niso.EQ.0) CYCLE
    IF (Fxr(ifxr,iz)%ldepl) THEN
      nFxrD = nFxrD+1
      ntfsrD = ntfsrD + Fxr(ifxr,iz)%nFsrInFxr
      IF (Fxr(ifxr,iz)%lres .AND. lFpln .AND. lrestrmt) THEN
        nResD = nResD+1
        IF (Fxr(ifxr,iz)%lCLD) nCldD = nCldD+1
      END IF
    ELSE
      nFxrTri = nFxrTri+1
      ntfsrTri = ntfsrTri+Fxr(ifxr,iz)%nFsrInFxr
      NtTri = NtTri+Fxr(ifxr,iz)%ndim
      nisomaxTri = max(nisomaxTri,Fxr(ifxr,iz)%ndim)
      IF (Fxr(ifxr,iz)%lres .AND. lFpln .AND. lrestrmt) THEN
        nResTri = nResTri+1
        NtRTri = NtRTri+Fxr(ifxr,iz)%ndim
        IF (Fxr(ifxr,iz)%lCLD) nCldTri = nCldTri+1
        IF (Fxr(ifxr,iz)%lAIC) nAICTri = nAICTri+1
      END IF
    END IF
  END DO
END DO
DO iz = (myze-myzb+1),31
  ibfxrD(iz+1) = nFxrD+1
  ibfxrT(iz+1) = nFxrTri+1
END DO

nGenD = nResD - nCldD; nGenTri = nResTri -nCldTri - nAICTri;

FxrD%ntfsr = ntfsrD; FxrTri%ntfsr = ntfsrTri;
FxrD%nFxrD = nFxrD; FxrD%nResD= nResD; FxrD%nCldD = nCldD; FxrD%nGenD = nGenD
FxrTri%nFxrTri = nFxrTri; FxrTri%NtRTri = NtRTri; FxrTri%NtTri = NtTri; FxrTri%nisomax = nisomaxTri
FxrTri%nResTri = nResTri; FxrTri%nCldTri = nCldTri;
FxrTri%nAICTri = nAICTri; FxrTri%nGenTri = nGenTri;

ALLOCATE(GlobalPinInfo%FxrStD(nxyz+1),GlobalPinInfo%FxrStT(nxyz+1))
ALLOCATE(GlobalPinInfo%GenStD(nxyz+1),GlobalPinInfo%GenStT(nxyz+1))
ALLOCATE(GlobalPinInfo%AICStT(nxyz+1))
! Common
IF (nFxrD.GT.0) THEN
  ALLOCATE(FxrD%mapG2R(nFxrD))
  FxrD%mapG2R = 0

  ALLOCATE(FxrD%mapfxrid(2*nFxrD),FxrD%mapglobalid(nFxrD))
  ALLOCATE(mapfsr(ntfsrD),FxrD%fsrvol(ntfsrD))
  ALLOCATE(lfuel(nFxrD), lCLD(nFxrD), lres(nFxrD), lGd(nFxrD), lh2o(nFxrD), lvoid(nFxrD))
  ALLOCATE(niso(nFxrD), ndim(nFxrD), ipin(nFxrD), fsrst(nFxrD), nfsr(nFxrD))
  ALLOCATE(temp(nFxrD), area(nFxrD), h2ofrac(nFxrD))
  ALLOCATE(FxrD%itTemp(nFxrD, ntiso), FxrD%wtTemp(nFXRD,ntiso))
  ALLOCATE(FxrD%itTempPx(nFxrD, ntiso), FxrD%wtTempPx(nFXRD,ntiso))
  ALLOCATE(idiso(ntiso,nFxrD), pnum(ntiso,nFxrD), chi(nchi,nFxrD))
  ALLOCATE(DMapNucl(ntiso,nFxrD),DMapNuclp13(ntiso,nFxrD),DMapFis(ntiso,nFxrD))
  ALLOCATE(FxrD%ind(ntiso,nFxrD))
  FxrD%lfuel=>lfuel; FxrD%lCLD=>lCLD; FxrD%lres=>lres; FxrD%lGd=>lGd;
  FxrD%lh2o=>lh2o; FxrD%lvoid => lvoid;
  FxrD%niso=>niso; FxrD%ndim=>ndim; FxrD%ipin=>ipin; FxrD%fsrst=>fsrst; FxrD%nfsr=>nfsr;
  FxrD%temp=>temp; FxrD%area=>area; FxrD%h2ofrac=>h2ofrac
  FxrD%idiso=>idiso; FxrD%pnum=>pnum; FxrD%chi=>chi;
  FxrD%Mapnucl=>DMapNucl; FxrD%MapNuclp13=>DMapNuclp13; FxrD%MapFis=>DMapFis
  FxrD%mapfsr => mapfsr
END IF
IF (nResD.GT.0) THEN
  ALLOCATE(FxrD%mapR2G(nResD),DMapNuclRes(ntiso,nResD))
  ALLOCATE(fresoa(iRgBeg:iRgEnd,nResD),fresof(iRgBeg:iRgEnd,nResD),fresonf(iRgBeg:iRgEnd,nResD),fresokf(iRgBeg:iRgEnd,nResD))
  ALLOCATE(fresoAIso(ntiso,iRgBeg:iRgEnd,nResD),fresoFIso(ntiso,iRgBeg:iRgEnd,nResD))
  FxrD%MapNuclRes=>DMapNuclRes
  FxrD%fresoa=>fresoa
  FxrD%fresof=>fresof
  FxrD%fresonf=>fresonf
  FxrD%fresokf=>fresokf
  FxrD%fresoAIso=>fresoAIso
  FxrD%fresoFIso=>fresoFIso
  IF (lRST) THEN
    ALLOCATE(fresos(iRgBeg:iRgEnd,nResD),fresostr(iRgBeg:iRgEnd,nResD))
    ALLOCATE(fresoSIso(ntiso,iRgBeg:iRgEnd,nResD))
    ALLOCATE(fresoSSIso(ntiso,iRgBeg:iRgEnd,nResD))
    ALLOCATE(fresoS1Iso(ntiso,iRgBeg:iRgEnd,nResD))
    FxrD%fresos=>fresos
    FxrD%fresostr=>fresostr
    FxrD%fresoSIso=>fresoSIso
    FxrD%fresoSSIso=>fresoSSIso
    FxrD%fresoS1Iso=>fresoS1Iso
  END IF
END IF
IF (nFxrTri.GT.0) THEN
  ALLOCATE(FxrTri%mapG2R(nFxrTri))
  FxrTri%mapG2R = 0;

  ALLOCATE(FxrTri%mapfxrid(2*nFxrTri),FxrTri%mapglobalid(nFxrTri))
  ALLOCATE(tmapfsr(ntfsrTri),FxrTri%fsrvol(ntfsrTri))
  ALLOCATE(AccNt(nFxrTri+1)); AccNt(nFxrTri+1) = NtTri
  ALLOCATE(tlAIC(nFxrTri), tlCLD(nFxrTri), tlres(nFxrTri), tlh2o(nFxrTri), tlvoid(nFxrTri))
  ALLOCATE(tniso(nFxrTri), tndim(nFxrTri), tipin(nFxrTri), tfsrst(nFxrTri), tnfsr(nFxrTri))
  ALLOCATE(ttemp(nFxrTri), tarea(nFxrTri), th2ofrac(nFxrTri))
  ALLOCATE(FxrTri%itTemp(nFxrTri, nisomaxTri), FxrTri%wtTemp(nFxrTri, nisomaxTri))
  ALLOCATE(tidiso(NtTri), tpnum(NtTri))
  ALLOCATE(TMapNucl(NtTri), TMapNuclp13(NtTri))
  FxrTri%lCLD=>tlCLD; FxrTri%lAIC=>tlAIC; FxrTri%lres=>tlres; FxrTri%lh2o=>tlh2o; FxrTri%lvoid=>tlvoid;
  FxrTri%niso=>tniso; FxrTri%ndim=>tndim; FxrTri%ipin=>tipin; FxrTri%fsrst=>tfsrst; FxrTri%nfsr=>tnfsr;
  FxrTri%temp=>ttemp; FxrTri%area=>tarea; FxrTri%h2ofrac=>th2ofrac
  FxrTri%idiso=>tidiso; FxrTri%pnum=>tpnum;
  FxrTri%MapNucl => TMapNucl; FxrTri%MapNuclp13=>TMapNuclp13;
  FxrTri%mapfsr => tmapfsr; FxrTri%AccNt => AccNt
END IF
IF (nResTri.GT.0) THEN
  ALLOCATE(AccNtR(nResTri+1)); AccNtR(nResTri+1) = NtRTri
  ALLOCATE(FxrTri%mapR2G(nResTri), TMapNuclRes(NtRTri))
  FxrTri%MapNuclRes => TMapNuclRes
  FxrTri%AccNtR=>AccNtR
  ALLOCATE(tfresoa(iRgBeg:iRgEnd,nResTri))
  ALLOCATE(tfresoAIso(NtRTri*(iRgEnd-iRgBeg+1)))
  FxrTri%fresoa=>tfresoa
  FxrTri%fresoAIso=>tfresoAIso
  IF (lRST) THEN
    ALLOCATE(tfresos(iRgBeg:iRgEnd,nResTri),tfresostr(iRgBeg:iRgEnd,nResTri))
    ALLOCATE(tfresoSIso(NtRTri*RngResG))
    ALLOCATE(tfresoSSIso(NtRTri*RngResG))
    ALLOCATE(tfresoS1Iso(NtRTri*RngResG))
    FxrTri%fresos=>tfresos
    FxrTri%fresostr=>tfresostr
    FxrTri%fresoSIso=>tfresoSIso; FxrTri%fresoSSIso=>tfresoSSIso; FxrTri%fresoS1Iso=>tfresoS1Iso
  END IF
END IF

IF (lrestrmt) THEN
  IF (lMLG) THEN
    IF (nFxrD.GT.0) THEN
      ALLOCATE(FxrD%mapG2C(nFxrD),FxrD%mapG2GR(nFxrD))
      FxrD%mapG2C = 0; FxrD%mapG2GR = 0;
    END IF
    IF (nCldD.GT.0) THEN
      ALLOCATE(xseq_c_1g(mlgdata0%c_nmaclv1G,nCldD))
      FxrD%xseq_c_1g=>xseq_c_1g;
      ALLOCATE(FxrD%mapC2G(nCldD))
    END IF
    IF (nGenD.GT.0) THEN
      ALLOCATE(xseq_f_mg(mlgdata0%f_nmaclv,iRgBeg:iRgEnd,nGenD))
      ALLOCATE(FnAdj(iRgBeg:iRgEnd,nGenD),FtAdj(mlgdata0%f_nmaclv,iRgBeg:iRgEnd,nGenD))
      FxrD%xseq_f_mg=>xseq_f_mg
      FxrD%FnAdj=>FnAdj
      FxrD%Ftadj=>FtAdj
      ALLOCATE(FxrD%mapGR2G(nGenD))
    END IF
    IF (nFxrTri.GT.0) THEN
      ALLOCATE(FxrTri%mapG2C(nFxrTri),FxrTri%mapG2A(nFxrTri),FxrTri%mapG2GR(nFxrTri))
      FxrTri%mapG2GR = 0; FxrTri%mapG2A = 0; FxrTri%mapG2C = 0
    END IF
    IF (nCldTri.GT.0) THEN
      ALLOCATE(txseq_c_1g(mlgdata0%c_nmaclv1G,nCldTri))
      FxrTri%xseq_c_1g=>txseq_c_1g;
      ALLOCATE(FxrTri%mapC2G(nCldTri))
    END IF
    IF (nAICTri.GT.0) THEN
      ALLOCATE(txseq_f_1g(mlgdata0%f_nmaclv1G,nAICTri))
      FxrTri%xseq_f_1g=>txseq_f_1g
      ALLOCATE(FxrTri%mapA2G(nAICTri))
    END IF
    IF (nGenTri.GT.0) THEN
      ALLOCATE(txseq_f_mg(mlgdata0%f_nmaclv,iRgBeg:iRgEnd,nGenTri))
      ALLOCATE(tFnAdj(iRgBeg:iRgEnd,nGenTri), tFtAdj(mlgdata0%f_nmaclv,iRgBeg:iRgEnd,nGenTri))
      FxrTri%xseq_f_mg=>txseq_f_mg
      FxrTri%FnAdj=>tFnAdj; FxrTri%Ftadj=>tFtAdj
      ALLOCATE(FxrTri%mapGR2G(nGenTri))
    END IF
  ELSE IF (lCAT) THEN
    ALLOCATE(xseq(nlvflxmax,nCat,iRgBeg:iRgEnd,nResD),NDAF(nlvflxmax,nCat,iRgBeg:iRgEnd,nResD))
    FxrD%xseq=>xseq
    FxrD%NDAF=>NDAF

    ALLOCATE(txseq(nlvflxmax,nCat,iRgBeg:iRgEnd,nResTri),tNDAF(nlvflxmax,nCat,iRgBeg:iRgEnd,nResTri))
    FxrTri%xseq=>txseq
    FxrTri%NDAF=>tNDAF
  ELSE
    ALLOCATE(xseq(nlvflxmax,nreshel,iRgBeg:iRgEnd,nResD),NDAF(nlvflxmax,nreshel,iRgBeg:iRgEnd,nResD))
    FxrD%xseq=>xseq
    FxrD%NDAF=>NDAF

    ALLOCATE(txseq(nlvflxmax,nreshel,iRgBeg:iRgEnd,nResTri),tNDAF(nlvflxmax,nreshel,iRgBeg:iRgEnd,nResTri))
    FxrTri%xseq=>txseq
    FxrTri%NDAF=>tNDAF
  END IF
END IF

ntfsrD = 0; ntfsrTri = 0;
nFxrD = 0 ; nResD = 0; nCldD = 0; nGenD = 0;
nFxrTri = 0; nResTri = 0; nCldTri = 0; nGenTri = 0; nAICTri = 0;
NtTri = 0; NtRTri = 0;
GlobalPinInfo%FxrStD = 0; GlobalPinInfo%FxrStT = 0
GlobalPinInfo%GenStD = 0; GlobalPinInfo%GenStT = 0
GlobalPinInfo%AICStT = 0;
!return
DO iz = myzb, myze
  lFpln = Core%lfuelplane(iz)
  DO ifxr = 1, nFxr
    myFxr => Fxr(ifxr,iz);
    IF (myFxr%niso.EQ.0) CYCLE
    FxrIdxSt = Core%Pin(myFxr%ipin)%FxrIdxSt; FsrIdxSt = Core%pin(myFxr%ipin)%FsrIdxSt;
    icel = Core%pin(myFxr%ipin)%Cell(iz)
    jdxfxr = ifxr - FxrIdxSt + 1;
    IF (Fxr(ifxr,iz)%ldepl) THEN
      nFxrD = nFxrD+1;
      FxrD%mapfxrid(nFxrD*2-1) = ifxr;
      FxrD%mapfxrid(nFxrD*2) = iz;
      FxrD%mapglobalid(nFxrD) = ifxr+nFxr*(iz-myzb)

      DO j = 1, myFxr%nFsrInFxr
        k = core%CellInfo(icel)%MapFxr2FsrIdx(j,jdxfxr)
        ifsr = FsrIdxSt+ k - 1
        mapfsr(ntfsrD+j) = ifsr
        FxrD%fsrvol(ntfsrD+j) = core%CellInfo(icel)%vol(k)
!        print*, ifsr
      END DO

      lfuel(nFxrD) = myFxr%lfuel; lCLD(nFxrD) = myFxr%lCLD; lres(nFxrD) = myFxr%lres;
      lGD(nFxrD) = myFxr%lGD; lh2o(nFxrD) = myFxr%lh2o; lvoid(nFxrD) = myFxr%lVoid
      niso(nFxrD) = myFxr%niso; ndim(nFxrD) = myFxr%ndim; ipin(nFxrD) = myFxr%ipin+(iz-myzb)*nxy;
      fsrst(nFxrD) = ntfsrD+1; nfsr(nFxrD) = myFxr%nFsrInFxr;
      temp(nFxrD) = myFxr%temp; area(nFxrD) = myFxr%area; h2ofrac(nFxrD) = myFxr%h2ofrac
      ntfsrD  = ntfsrD+myFxr%nFsrInFxr

      GlobalPinInfo%FxrStD(ipin(nFxrD)+1) = GlobalPinInfo%FxrStD(ipin(nFxrD)+1)+1

      idiso(:,nFxrD) = myFxr%idiso(:); pnum(:,nFxrD) = myFxr%pnum(:);
      DEALLOCATE(myFxr%idiso, myFxr%pnum)
      myFxr%idiso(1:ntiso) => idiso(1:ntiso,nFxrD); myFxr%pnum(1:ntiso) => pnum(1:ntiso,nFxrD)
      DO iso = 1, niso(nFxrD)
        DMapNucl(iso,nFxrD) = mapnucl(idiso(iso,nFxrD))
        p1idx = mapnuclp13(idiso(iso,nFxrD))
        IF (p1idx.GT.0) p1idx = idnp1hel(p1idx)
        DMapNuclp13(iso,nFxrD) = p1idx
      END DO

      if (lfuel(nFxrD)) THEN
        chi(:,nFxrD) = myFxr%chi(:);
        DEALLOCATE(myFxr%chi)
        myFxr%chi(1:nchi) => chi(1:nchi,nFxrD)
        DO iso = 1, niso(nFxrD)
          DMapFis(iso,nFxrD) = mapfis(idiso(iso,nFxrD))
        END DO
      END IF

      IF (Fxr(ifxr,iz)%lres .AND. lFpln .AND. lrestrmt) THEN
        nResD = nResD+1
        FxrD%mapR2G(nResD) = nFxrD; FxrD%mapG2R(nFxrD) = nResD
        DO iso = 1, niso(nFxrD)
          DMapNuclRes(iso,nResD) = mapnuclRes(idiso(iso,nFxrD),iz)
        END DO
        DO iso = 1, niso(nFxrD)
          k = idiso(iso,nFxrD)
          k = mapnucl(k)
          FxrD%ind(k,nFxrD) = 0.
          DO j = 1, ntiso
            IF (DMapNuclRes(iso,nResD).EQ.DMapNuclRes(j,nResD)) THEN
              FxrD%ind(k,nFxrD) = FxrD%ind(k,nFxrD)+pnum(j,nfxrD)
            END IF
          END DO
        END DO

        fresoa(:,nResD) = myFxr%fresoa(:); fresof(:,nResD) = myFxr%fresof(:);
        fresonf(:,nResD) = myFxr%fresonf(:); fresokf(:,nResD) = myFxr%fresokf(:);
        fresoAIso(:,:,nResD) = myFxr%fresoAIso(:,:); fresoFIso(:,:,nResD) = myFxr%fresoFIso(:,:)

        DEALLOCATE(myFxr%fresoa,myFxr%fresof,myFxr%fresonf,myFxr%fresokf)
        DEALLOCATE(myFxr%fresoAIso, myFxr%fresoFIso)

        myFxr%fresoa(iRgBeg:iRgEnd) => fresoa(iRgBeg:iRgEnd,nResD)
        myFxr%fresof(iRgBeg:iRgEnd) => fresof(iRgBeg:iRgEnd,nResD)
        myFxr%fresonf(iRgBeg:iRgEnd) => fresonf(iRgBeg:iRgEnd,nResD)
        myFxr%fresokf(iRgBeg:iRgEnd) => fresokf(iRgBeg:iRgEnd,nResD)
        myFxr%fresoAIso(1:,iRgBeg:) => fresoAIso(1:ntiso,iRgBeg:iRgEnd,nResD)
        myFxr%fresoFIso(1:,iRgBeg:) => fresoFIso(1:ntiso,iRgBeg:iRgEnd,nResD)

        IF (lRST) THEN
          fresos(:,nResD) = myFxr%fresos(:); fresostr(:,nResD) = myFxr%fresostr(:);
          fresoSIso(:,:,nResD) = myFxr%fresoSIso(:,:)
          fresoSSIso(:,:,nResD) = myFxr%fresoSSIso(:,:)
          fresoS1Iso(:,:,nResD) = myFxr%fresoS1Iso(:,:)
          DEALLOCATE(myFxr%fresos, myFxr%fresostr)
          DEALLOCATE(myFxr%fresoSIso, myFxr%fresoSSIso, myFxr%fresoS1Iso)
          myFxr%fresos(iRgBeg:iRgEnd)=>fresos(iRgBeg:iRgEnd,nResD)
          myFxr%fresostr(iRgBeg:iRgEnd)=>fresostr(iRgBeg:iRgEnd,nResD)
          myFxr%fresoSIso(1:,iRgBeg:)=>fresoSIso(1:ntiso,iRgBeg:iRgEnd,nResD)
          myFxr%fresoSSIso(1:,iRgBeg:)=>fresoSSIso(1:ntiso,iRgBeg:iRgEnd,nResD)
          myFxr%fresoS1Iso(1:,iRgBeg:)=>fresoS1Iso(1:ntiso,iRgBeg:iRgEnd,nResD)
        END IF

        IF (.NOT.lMLG) THEN
          IF (lCat) THEN
            xseq(:,:,:,nResD) = myFxr%xseq(:,:,:); NDAF(:,:,:,nResD) = myFxr%NDAF(:,:,:)
            DEALLOCATE(myFxr%xseq,myFxr%NDAF)
            myFxr%xseq(1:,1:,iRgBeg:)=>xseq(1:nlvflxmax,1:nCat,iRgBeg:iRgEnd,nResD);
            myFxr%NDAF(1:,1:,iRgBeg:)=>NDAF(1:nlvflxmax,1:nCat,iRgBeg:iRgEnd,nResD)
          ELSE
            xseq(:,:,:,nResD) = myFxr%xseq(:,:,:); NDAF(:,:,:,nResD) = myFxr%NDAF(:,:,:)
            DEALLOCATE(myFxr%xseq,myFxr%NDAF)
            myFxr%xseq(1:,1:,iRgBeg:)=>xseq(1:nlvflxmax,1:nreshel,iRgBeg:iRgEnd,nResD);
            myFxr%NDAF(1:,1:,iRgBeg:)=>NDAF(1:nlvflxmax,1:nreshel,iRgBeg:iRgEnd,nResD)
          END IF
        END IF

        IF (Fxr(ifxr,iz)%lCLD) THEN
          IF (lMLG) THEN
            nCldD = nCldD+1
            FxrD%mapC2G(nCldD) = nFxrD; FxrD%mapG2C(nFxrD) = nCldD
            xseq_c_1g(:,nCldD)=myFxr%xseq_c_1g
            DEALLOCATE(myFxr%xseq_c_1g)
            myFxr%xseq_c_1g(1:mlgdata0%c_nmaclv1G)=>xseq_c_1g(1:mlgdata0%c_nmaclv1G,nCldD);
          END IF
        ELSE
          IF (lMLG) THEN
            nGenD = nGenD+1
            FxrD%mapGR2G(nGenD) = nFxrD; FxrD%mapG2GR(nFxrD) = nGenD
            xseq_f_mg(:,:,nGenD)=myFxr%xseq_f_mg
            FnAdj(:,nGenD) = myFxr%FnAdj(:)
            FtAdj(:,:,nGenD) = myFxr%FtAdj(:,:)
            DEALLOCATE(myFxr%xseq_f_mg, myFxr%FnAdj, myFxr%FtAdj)
            myFxr%xseq_f_mg(1:,iRgBeg:)=>xseq_f_mg(1:mlgdata0%f_nmaclv,iRgBeg:iRgEnd,nGenD)
            myFxr%FnAdj(iRgBeg:iRgEnd)=>FnAdj(iRgBeg:iRgEnd,nGenD)
            myFxr%FtAdj(1:,iRgBeg:)=>Ftadj(1:mlgdata0%f_nmaclv,iRgBeg:iRgEnd,nGenD)

            GlobalPinInfo%GenStD(ipin(nFxrD)+1) = GlobalPinInfo%GenStD(ipin(nFxrD)+1)+1
          END IF
        END IF
      END IF
    ELSE
      nFxrTri = nFxrTri+1;
      AccNt(nFxrTri) = NtTri
      FxrTri%mapfxrid(nFxrTri*2-1) = ifxr;
      FxrTri%mapfxrid(nFxrTri*2) = iz;
      FxrTri%mapglobalid(nFxrTri) = ifxr+nFxr*(iz-myzb)
      mdim = Fxr(ifxr,iz)%ndim

      DO j = 1, myFxr%nFsrInFxr
        k = core%CellInfo(icel)%MapFxr2FsrIdx(j, jdxfxr)
        ifsr = FsrIdxSt + k - 1
        tmapfsr(ntfsrTri+j) = ifsr
        FxrTri%fsrvol(ntfsrTri+j) = core%CellInfo(icel)%vol(k)
!        print*, ifsr
      END DO

      tlAIC(nFxrTri) = myFxr%lAIC; tlCLD(nFxrTri) = myFxr%lCLD; tlres(nFxrTri) = myFxr%lres;
      tlh2o(nFxrTri) = myFxr%lh2o; tlvoid(nFxrTri) = myFxr%lVoid
      tniso(nFxrTri) = myFxr%niso; tndim(nFxrTri) = myFxr%ndim; tipin(nFxrTri) = myFxr%ipin+(iz-myzb)*nxy;
      tfsrst(nFxrTri) = ntfsrTri+1; tnfsr(nFxrTri) = myFxr%nFsrInFxr;
      ttemp(nFxrTri) = myFxr%temp; tarea(nFxrTri) = myFxr%area; th2ofrac(nFxrTri) = myFxr%h2ofrac
      ntfsrTri = ntfsrTri+myFxr%nFsrInFxr

      GlobalPinInfo%FxrStT(tipin(nFxrTri)+1) = GlobalPinInfo%FxrStT(tipin(nFxrTri)+1)+1

      tidiso(NtTri+1:NtTri+mdim) = myFxr%idiso(:)
      tpnum(NtTri+1:NtTri+mdim) = myFxr%pnum(:);
      DEALLOCATE(myFxr%idiso, myFxr%pnum)
      myFxr%idiso(1:mdim) => tidiso(NtTri+1:NtTri+mdim); myFxr%pnum(1:mdim) => tpnum(NtTri+1:NtTri+mdim)

      DO iso = NtTri+1, NtTri+mdim
        IF (tidiso(iso).EQ.0) CYCLE
        TMapNucl(iso) = mapnucl(tidiso(iso))
        p1idx = mapnuclp13(tidiso(iso))
        IF (p1idx.GT.0) p1idx = idnp1hel(p1idx)
        TMapNuclp13(iso) = p1idx
      END DO

      NtTri = NtTri + mdim
      IF (Fxr(ifxr,iz)%lres .AND. lFpln .AND. lrestrmt) THEN
        nResTri = nResTri+1
        FxrTri%mapR2G(nResTri) = nFxrTri; FxrTri%mapG2R(nFxrTri) = nResTri
        AccNtR(nResTri) = NtRTri
        DO iso = NtRTri+1, NtRTri+mdim
          TMapNuclRes(iso) = mapnuclRes(tidiso(iso),iz)
        END DO

        tfresoa(:,nResTri) = myFxr%fresoa(:)
        DO ig = iRgBeg, iRgEnd
          tfresoAIso(NtRTri*RngResG+(ig-iRgBeg)*mdim+1:NtRTri*RngResG+(ig-iRgBeg+1)*mdim) = myFxr%fresoAIso(:,ig)
        END DO
        DEALLOCATE(myFxr%fresoa, myFxr%fresoAIso)
        myFxr%fresoa(iRgBeg:iRgEnd) => tfresoa(iRgBeg:iRgEnd,nResTri)
        myFxr%fresoAIso(1:mdim,iRgBeg:iRgEnd) => tfresoAIso(NtRTri*RngResG+1:(NtRTri+mdim)*RngResG)

        IF (lRST) THEN
          tfresos(:,nResTri) = myFxr%fresos(:); tfresostr(:,nResTri) = myFxr%fresostr(:);
          DO iso = 1, mdim
            tfresoSIso((NtTri+iso-1)*RngResG+1:(NtRTri+iso)*RngResG) = myFxr%fresoSIso(iso,:)
            tfresoSSIso((NtTri+iso-1)*RngResG+1:(NtRTri+iso)*RngResG) = myFxr%fresoSSIso(iso,:)
            tfresoS1Iso((NtTri+iso-1)*RngResG+1:(NtRTri+iso)*RngResG) = myFxr%fresoS1Iso(iso,:)
          END DO
          DEALLOCATE(myFxr%fresos, myFxr%fresostr)
          DEALLOCATE(myFxr%fresoSIso, myFxr%fresoSSIso, myFxr%fresoS1Iso)
          myFxr%fresos(iRgBeg:iRgEnd)=>tfresos(iRgBeg:iRgEnd,nResTri)
          myFxr%fresostr(iRgBeg:iRgEnd)=>tfresostr(iRgBeg:iRgEnd,nResTri)
          myFxr%fresoSIso(1:mdim,iRgBeg:iRgEnd)=>tfresoSIso(NtRTri*RngResG+1:(NtRTri+mdim)*RngResG)
          myFxr%fresoSSIso(1:mdim,iRgBeg:iRgEnd)=>tfresoSSIso(NtRTri*RngResG+1:(NtRTri+mdim)*RngResG)
          myFxr%fresoS1Iso(1:mdim,iRgBeg:iRgEnd)=>tfresoS1Iso(NtRTri*RngResG+1:(NtRTri+mdim)*RngResG)
        END IF

        NtRTri = NtRTri+mdim
        IF (.NOT.lMLG) THEN
          IF (lCat) THEN
            txseq(:,:,:,nResTri) = myFxr%xseq(:,:,:); tNDAF(:,:,:,nResTri) = myFxr%NDAF(:,:,:)
            DEALLOCATE(myFxr%xseq,myFxr%NDAF)
            myFxr%xseq(1:,1:,iRgBeg:)=>txseq(1:nlvflxmax,1:nCat,iRgBeg:iRgEnd,nResTri);
            myFxr%NDAF(1:,1:,iRgBeg:)=>tNDAF(1:nlvflxmax,1:nCat,iRgBeg:iRgEnd,nResTri)
          ELSE
            txseq(:,:,:,nResTri) = myFxr%xseq(:,:,:); tNDAF(:,:,:,nResTri) = myFxr%NDAF(:,:,:)
            DEALLOCATE(myFxr%xseq,myFxr%NDAF)
            myFxr%xseq(1:,1:,iRgBeg:)=>txseq(1:nlvflxmax,1:nreshel,iRgBeg:iRgEnd,nResTri);
            myFxr%NDAF(1:,1:,iRgBeg:)=>tNDAF(1:nlvflxmax,1:nreshel,iRgBeg:iRgEnd,nResTri)
          END IF
        END IF

        IF (Fxr(ifxr,iz)%lCLD) THEN
          IF (lMLG) THEN
            nCldTri = nCldTri+1
            FxrTri%mapC2G(nCldTri) = nFxrTri; FxrTri%mapG2C(nFxrTri) = nCldTri
            txseq_c_1g(:,nCldTri)=myFxr%xseq_c_1g(:)
            DEALLOCATE(myFxr%xseq_c_1g)
            myFxr%xseq_c_1g(1:mlgdata0%c_nmaclv1G)=>txseq_c_1g(1:mlgdata0%c_nmaclv1G,nCldTri)
          END IF
        ELSE IF (Fxr(ifxr,iz)%lAIC) THEN
          IF (lMLG) THEN
            nAICTri = nAICTri+1
            FxrTri%mapA2G(nAICTri) = nFxrTri; FxrTri%mapG2A(nFxrTri) = nAICTri
            txseq_f_1g(:,nAICTri)=myFxr%xseq_f_1g(:)
            DEALLOCATE(myFxr%xseq_f_1g)
            myFxr%xseq_f_1g(1:mlgdata0%f_nmaclv1G)=>txseq_f_1g(1:mlgdata0%f_nmaclv1G,nAICTri)
            GlobalPinInfo%AICStT(tipin(nFxrTri)+1) = GlobalPinInfo%AICStT(tipin(nFxrTri)+1)+1
          END IF
        ELSE
          IF (lMLG) THEN
            nGenTri = nGenTri;
            FxrTri%mapGR2G(nGenTri) = nFxrTri; FxrTri%mapG2GR(nFxrTri) = nGenTri
            txseq_f_mg(:,:,nGenTri) = myFxr%xseq_f_mg(:,:)
            tFnAdj(:,nGenTri) = myFxr%FnAdj(:)
            tFtadj(:,:,nGenTri) = myFxr%FtAdj(:,:)
            DEALLOCATE(myFxr%xseq_f_mg, myFxr%FnAdj, myFxr%FtAdj)
            myFxr%xseq_f_mg(1:,iRgBeg:)=>txseq_f_mg(1:mlgdata0%f_nmaclv,iRgBeg:iRgEnd,nGenTri)
            myFxr%FnAdj(iRgBeg:iRgEnd)=>tFnAdj(iRgBeg:iRgEnd,nGenTri)
            myFxr%FtAdj(1:,iRgBeg:)=>tFtadj(1:mlgdata0%f_nmaclv,iRgBeg:iRgEnd,nGenTri)

            GlobalPinInfo%GenStT(tipin(nFxrTri)+1) = GlobalPinInfo%GenStT(tipin(nFxrTri)+1)+1
          END IF
        END IF
      END IF
    END IF
  END DO
END DO

GlobalPinInfo%FxrStD(1) = 1
GlobalPinInfo%FxrStT(1) = 1
GlobalPinInfo%GenStD(1) = 1
GlobalPinInfo%GenStT(1) = 1
GlobalPinInfo%AICStT(1) = 1
DO j = 1, nxyz
  GlobalPinInfo%FxrStD(j+1) = GlobalPinInfo%FxrStD(j)+GlobalPinInfo%FxrStD(j+1)
  GlobalPinInfo%FxrStT(j+1) = GlobalPinInfo%FxrStT(j)+GlobalPinInfo%FxrStT(j+1)
  GlobalPinInfo%GenStD(j+1) = GlobalPinInfo%GenStD(j)+GlobalPinInfo%GenStD(j+1)
  GlobalPinInfo%GenStT(j+1) = GlobalPinInfo%GenStT(j)+GlobalPinInfo%GenStT(j+1)
  GlobalPinInfo%AICStT(j+1) = GlobalPinInfo%AICStT(j)+GlobalPinInfo%AICStT(j+1)
END DO

END SUBROUTINE

SUBROUTINE AlignResVarPin
USE GEOM, ONLY : Core
USE PE_Mod, ONLY : PE
USE XSLIB_MOD, ONLY : nreshel, nelthel, mapnucl
USE Core_mod, ONLY : GroupInfo
IMPLICIT NONE
INTEGER :: myzb, myze, nxy, nxyz
INTEGER :: irgb, irge, f_nmaclv, f_nmaclv1G
INTEGER :: i, j, ipin, iso
IF (.NOT. lrestrmt) RETURN

myzb = PE%myzb; myze = PE%myze; nxy = Core%nxy; nxyz = nxy*(myze-myzb+1);
irgb = GroupInfo%nofg+1; irge = GroupInfo%nofg+GroupInfo%norg
ALLOCATE(ResPin%igt(nxyz),ResPin%lres(nxyz),ResPin%lresA(nxyz),ResPin%lresC(nxyz))
ALLOCATE(ResPin%idiso(nelthel,nxyz),ResPin%pnum(nelthel,nxyz),ResPin%niso(nxyz))
ALLOCATE(ResPin%avgxseq_mg(MLGLib%f_nmaclv,irgb:irge,nxyz))
ALLOCATE(ResPin%avgxseq_1g(MLGLib%f_nmaclv1G,nxyz))
ALLOCATE(ResPin%FnAdj(irgb:irge,nxyz))
ALLOCATE(ResPin%rifa(nreshel,irgb:irge,nxyz), ResPin%riff(nreshel,irgb:irge,nxyz))
ALLOCATE(ResPin%temp(nxyz))
ALLOCATE(ResPin%mapnucl(nelthel,nxyz))
ALLOCATE(ResPin%ind(nelthel,nxyz))
ALLOCATE(GlobalPinInfo%lAIC(nxyz), GlobalPinInfo%lres(nxyz))

ipin = 0;
DO i = myzb, myze
  DO j = 1, nxy
    ipin = ipin+1
    ResPin%igt(ipin) = Core%ResVarPin(j,i)%igt
    ResPin%lres(ipin) = Core%ResVarPin(j,i)%lres
    ResPin%lresA(ipin) = Core%ResVarPin(j,i)%lresA
    ResPin%lresC(ipin) = Core%ResVarPin(j,i)%lresC

    GlobalPinInfo%lAIC(ipin) = Core%CellInfo(Core%Pin(j)%Cell(i))%lAIC
    GlobalPinInfo%lres(ipin) = Core%CellInfo(Core%Pin(j)%Cell(i))%lres

    IF (.NOT. core%lFuelPlane(i)) CYCLE
    IF (.NOT. ResPin%lres(ipin)) CYCLE

    ResPin%idiso(:,ipin) = Core%ResVarPin(j,i)%idiso
    DEALLOCATE(Core%ResVarPin(j,i)%idiso)
    Core%ResVarPin(j,i)%idiso => ResPin%idiso(:,ipin)

    ResPin%pnum(:,ipin) = Core%ResVarPin(j,i)%pnum
    DEALLOCATE(Core%ResVarPin(j,i)%pnum)
    Core%ResVarPin(j,i)%pnum => ResPin%pnum(:,ipin)

    ResPin%niso(ipin) = Core%ResVarPin(j,i)%niso

    ResPin%avgxseq_mg(:,:,ipin) = Core%ResVarPin(j,i)%avgxseq_mg
    DEALLOCATE(Core%ResVarPin(j,i)%avgxseq_mg)
    Core%ResVarPin(j,i)%avgxseq_mg(1:,irgb:) => ResPin%avgxseq_mg(:,:,ipin)

    ResPin%avgxseq_1g(:,ipin) = Core%ResVarPin(j,i)%avgxseq_1g
    DEALLOCATE(Core%ResVarPin(j,i)%avgxseq_1g)
    Core%ResVarPin(j,i)%avgxseq_1g => ResPin%avgxseq_1g(:,ipin)

    ResPin%FnAdj(:,ipin) = Core%ResVarPin(j,i)%FnAdj
    DEALLOCATE(Core%ResVarPin(j,i)%FnAdj)
    Core%ResVarPin(j,i)%FnAdj(irgb:) => ResPin%FnAdj(:,ipin)

    ResPin%rifa(:,:,ipin) = Core%ResVarPin(j,i)%rifa
    DEALLOCATE(Core%ResVarPin(j,i)%rifa)
    Core%ResVarPin(j,i)%rifa(1:,irgb:) => ResPin%rifa(:,:,ipin)

    ResPin%riff(:,:,ipin) = Core%ResVarPin(j,i)%riff
    DEALLOCATE(Core%ResVarPin(j,i)%riff)
    Core%ResVarPin(j,i)%riff(1:,irgb:) => ResPin%riff(:,:,ipin)

    DO iso = 1, nelthel
      IF (iso.GT.ResPin%niso(ipin)) THEN
        ResPin%mapnucl(iso,ipin) = 0
      ELSE
        ResPin%mapnucl(iso,ipin) = mapnucl(ResPin%idiso(iso,ipin))
      END IF
    END DO
  END DO
END DO

END SUBROUTINE

SUBROUTINE AlignDeplFxrMem(Fxr, GroupInfo)
USE TYPEDEF,    ONLY : Fxrinfo_type, GroupInfo_Type
USE GEOM,       ONLY : Core
USE PE_Mod,     ONLY : PE
USE DeplLib_MOD, ONLY : GetnDeplIso

IMPLICIT NONE

TYPE(Fxrinfo_type), POINTER :: Fxr(:,:)
TYPE(GroupInfo_Type) :: GroupInfo

TYPE(Fxrinfo_type), POINTER :: myFxr
INTEGER :: ntiso, nisoDepl, nFxrD, myzb, myze, nFxr
INTEGER :: iz, ifxr, iso

ntiso = GroupInfo%ntiso; nisoDepl = GetnDeplIso(); nFxr = Core%nCoreFxr
myzb = PE%myzb; myze = PE%myze;
nFxrD = FxrD%nFxrD;

IF (nFxrD.LE.0) RETURN

ALLOCATE(FxrD%niso_depl(nFxrD), FxrD%niso_past(nFxrD))
ALLOCATE(FxrD%idiso_past(ntiso,nFxrD), FxrD%pnum_past(ntiso,nFxrD))
ALLOCATE(FxrD%pnum_all(nisoDepl,nFxrD), FxrD%pnum_past_all(nisoDepl,nFxrD))
ALLOCATE(FxrD%l_pnum_all(nFxrD))
FxrD%l_pnum_all=.FALSE.

nFxrD = 0;
DO iz = myzb, myze
  DO ifxr = 1, nFxr
    myFxr => Fxr(ifxr,iz)
    IF (myFxr%niso.EQ.0) CYCLE
    IF (myFxr%ldepl) THEN
      nFxrD = nFxrD+1;
      DEALLOCATE(myFxr%idiso_past,myFxr%pnum_all,myFxr%pnum_past,myFxr%pnum_past_all)
      myFxr%idiso_past=>FxrD%idiso_past(1:ntiso,nFxrD)
      myFxr%pnum_past=>FxrD%pnum_past(1:ntiso,nFxrD)

      myFxr%pnum_all=>FxrD%pnum_all(1:nisoDepl,nFxrD)
      myFxr%pnum_past_all=>FxrD%pnum_past_all(1:nisoDepl,nFxrD)
    END IF
  END DO
END DO

END SUBROUTINE

SUBROUTINE UpdateBlockFxr(Fxr, GroupInfo, caseid)
USE TYPEDEF, ONLY : Fxrinfo_type, GroupInfo_Type
USE XSLIB_MOD,  ONLY : mapnucl, mapnuclp13, mapnuclRes, mapfis, idnp1hel
USE OMP_LIB
IMPLICIT NONE
TYPE(Fxrinfo_type), POINTER :: Fxr(:,:)
TYPE(GroupInfo_Type) :: GroupInfo
INTEGER :: caseid

TYPE(Fxrinfo_type), POINTER :: myFxr
INTEGER :: i, j, k, m, idiso
LOGICAL :: lremap

SELECT CASE(caseid)
  CASE(DEPLFeed)
    !$OMP PARALLEL PRIVATE(myFxr, i, j, m)
    !$OMP DO SCHEDULE(GUIDED)
    DO i = 1, FxrD%nFxrD
      myFxr => Fxr(FxrD%mapfxrid(i*2-1), FxrD%mapfxrid(i*2))
      FxrD%niso(i) = myFxr%niso;
    END DO
    !$OMP END DO
    !$OMP DO SCHEDULE(GUIDED)
    DO i = 1, FxrD%nFxrD
      myFxr => Fxr(FxrD%mapfxrid(i*2-1), FxrD%mapfxrid(i*2))
      FxrD%niso_depl(i) = myFxr%niso_depl;
    END DO
    !$OMP END DO
    !$OMP DO SCHEDULE(GUIDED)
    DO i = 1, FxrD%nFxrD
      myFxr => Fxr(FxrD%mapfxrid(i*2-1), FxrD%mapfxrid(i*2))
      FxrD%niso_past(i) = myFxr%niso_past;
    END DO
    !$OMP END DO
    !$OMP DO SCHEDULE(GUIDED)
    DO i = 1, FxrD%nFxrD
      myFxr => Fxr(FxrD%mapfxrid(i*2-1), FxrD%mapfxrid(i*2))
      DO j = 1, FxrD%niso(i)
        FxrD%MapNucl(j,i) = mapnucl(FxrD%idiso(j,i))
      END DO
    END DO
    !$OMP END DO
    !$OMP DO SCHEDULE(GUIDED)
    DO i = 1, FxrD%nFxrD
      myFxr => Fxr(FxrD%mapfxrid(i*2-1), FxrD%mapfxrid(i*2))
      IF (myFxr%lfuel) THEN
        DO j = 1, FxrD%niso(i)
          FxrD%MapFis(j,i) = mapfis(FxrD%idiso(j,i))
        END DO
      END IF
    END DO
    !$OMP END DO
    !$OMP DO SCHEDULE(GUIDED)
    DO i = 1, FxrD%nFxrD
      myFxr => Fxr(FxrD%mapfxrid(i*2-1), FxrD%mapfxrid(i*2))
      DO j = 1, FxrD%niso(i)
        m = mapnuclp13(FxrD%idiso(j,i))
        IF(m.ne.0) FxrD%MapNuclp13(j,i) = idnp1hel(m)
      END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
    IF (lrestrmt) THEN
      !$OMP PARALLEL PRIVATE(myFxr, i, j, m, idiso)
      !$OMP DO SCHEDULE(GUIDED)
      DO k = 1, FxrD%nResD
        i = FxrD%mapR2G(k)
        IF (i.LT.1) CYCLE
        myFxr => Fxr(FxrD%mapfxrid(i*2-1), FxrD%mapfxrid(i*2))
        DO j = 1, FxrD%niso(i)
          IF (myFxr%lres) FxrD%MapNuclRes(j,k) = mapnuclRes(FxrD%idiso(j,i),FxrD%mapfxrid(i*2))
        END DO
        DO j = 1, FxrD%niso(i)
          idiso = FxrD%MapNucl(j,i)
          FxrD%ind(idiso,i) = 0.
          DO m = 1, FxrD%niso(i)
            IF (FxrD%MapNuclRes(j,k).EQ.FxrD%MapNuclRes(m,k)) FxrD%ind(idiso,i) = FxrD%ind(idiso,i)+FxrD%pnum(m,i)
          END DO
        END DO
      END DO
      !$OMP END DO
      !$OMP END PARALLEL
    END IF
  CASE(XeDynFeed)
    !$OMP PARALLEL PRIVATE(myFxr, i, j, m, lremap, idiso)
    !$OMP DO SCHEDULE(GUIDED)
    DO i = 1, FxrD%nFxrD
      myFxr => Fxr(FxrD%mapfxrid(i*2-1), FxrD%mapfxrid(i*2))
      lremap = (FxrD%niso(i).NE.myFxr%niso)
      IF (lremap) THEN
        FxrD%niso(i) = myFxr%niso; FxrD%niso_depl(i) = myFxr%niso_depl
        DO j = 1, FxrD%niso(i)
          FxrD%MapNucl(j,i) = mapnucl(FxrD%idiso(j,i))
          IF(myFxr%lfuel) FxrD%MapFis(j,i) = mapfis(FxrD%idiso(j,i))
          m = mapnuclp13(FxrD%idiso(j,i))
          IF(m.ne.0) FxrD%MapNuclp13(j,i) = idnp1hel(m)
          IF (myFxr%lres .AND. lFPlane(FxrD%mapfxrid(i*2))) FxrD%MapNuclRes(j,FxrD%mapG2R(i)) = mapnuclRes(FxrD%idiso(j,i),FxrD%mapfxrid(i*2))
        END DO
      END IF
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
    !$OMP PARALLEL PRIVATE(i, j, m, idiso)
    !$OMP DO SCHEDULE(GUIDED)
    DO k = 1, FxrD%nResD
      i = FxrD%mapR2G(k)
      DO j = 1, FxrD%niso(i)
        idiso = FxrD%MapNucl(j,i)
        FxrD%ind(idiso,i) = 0.
        DO m = 1, FxrD%niso(i)
          IF (FxrD%MapNuclRes(j,k).EQ.FxrD%MapNuclRes(m,k)) FxrD%ind(idiso,i) = FxrD%ind(idiso,i)+FxrD%pnum(m,i)
        END DO
      END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
  CASE(THFeed)
    !$OMP PARALLEL PRIVATE(myFxr, i)
    !$OMP DO SCHEDULE(GUIDED)
    DO i = 1, FxrD%nFxrD
      myFxr => Fxr(FxrD%mapfxrid(i*2-1), FxrD%mapfxrid(i*2))
      FxrD%temp(i) = myFxr%temp
    END DO
    !$OMP END DO
    !$OMP DO SCHEDULE(GUIDED)
    DO i = 1, FxrTri%nFxrTri
      myFxr => Fxr(FxrTri%mapfxrid(i*2-1),FxrTri%mapfxrid(i*2))
      FxrTri%temp(i) = myFxr%temp
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
  CASE(PPMFeed)
!    IF (.NOT.lInitBoron) THEN
      !$OMP PARALLEL PRIVATE(myFxr, i, j, m, k, lremap)
      !$OMP DO SCHEDULE(GUIDED)
      DO i = 1, FxrTri%nFxrTri
        myFxr => Fxr(FxrTri%mapfxrid(i*2-1),FxrTri%mapfxrid(i*2))
        lremap = (FxrTri%niso(i).NE.myFxr%niso)
        IF (.NOT. lremap) CYCLE
        FxrTri%niso(i) = myFxr%niso
        DO j = 1, FxrTri%niso(i)
          k = j+FxrTri%AccNt(i)
          FxrTri%MapNucl(k) = mapnucl(FxrTri%idiso(k))
          m = mapnuclp13(FxrTri%idiso(k))
          IF (m.ne.0) FxrTri%MapNuclp13(k) = idnp1hel(m)
        END DO
        IF (myFxr%lres .AND. lFPlane(FxrTri%mapfxrid(i*2))) THEN
          DO j = 1, FxrTri%niso(i)
            k = j+FxrTri%AccNt(i)
            m = j+FxrTri%AccNtR(FxrTri%mapG2R(i))
            FxrTri%MapNuclRes(m) = mapnuclRes(FxrTri%idiso(k),FxrTri%mapfxrid(i*2))
          END DO
        END IF
      END DO
      !$OMP END DO
      !$OMP END PARALLEL
      lInitBoron = .TRUE.
 !   END IF
END SELECT

END SUBROUTINE

SUBROUTINE XsLinIntPolTemp(ScatOrder)
  USE XSLIB_MOD, ONLY : itempmap, itempmapp1
  USE OMP_LIB
  IMPLICIT NONE
  INTEGER :: ScatOrder

  INTEGER :: i, j, k
  INTEGER :: nFxr, niso, id, it, it1, ntemp, itemp
  REAL :: temp

  nFxr = FxrD%nFxrD
  !$OMP PARALLEL PRIVATE(i, j, niso, temp, it, id, it1, itemp, ntemp)
  !$OMP DO SCHEDULE(GUIDED)
  DO i = 1, nFxr
    niso = FxrD%niso(i)
    temp = FxrD%temp(i)
    it = temp
    DO j = 1, niso
      id = FxrD%MapNucl(j,i)
      it1 = itempmap(it, id)
      itemp = IsoData%ptrTemp(id)
      ntemp = IsoData%ptrTemp(id+1)- itemp
      it1 = it1+itemp
      FxrD%itTemp(i,j) = it1
      FxrD%wtTemp(i,j) = 1.
      IF (ntemp .LE. 1) CYCLE
      IF (it1 .GE. ntemp+itemp) CYCLE
      FxrD%wtTemp(i,j) = (-temp + IsoData%ptsTemp(it1+1))/(IsoData%ptsTemp(it1+1)-IsoData%ptsTemp(it1))
    END DO
    IF (ScatOrder < 1) CYCLE
    DO j = 1, niso
      id = FxrD%MapNuclp13(j,i)
      it1 = itempmapp1(it, id)
      itemp = IsoData%ptrTempPx(id)
      ntemp = IsoData%ptrTempPx(id+1)- itemp
      it1 = it1+itemp
      FxrD%itTempPx(i,j) = it1
      FxrD%wtTempPx(i,j) = 1.
      IF (ntemp .LE. 1) CYCLE
      IF (it1 .GE. ntemp+itemp) CYCLE
      FxrD%wtTempPx(i,j) = (-temp + IsoData%ptsTempPx(it1+1))/(IsoData%ptsTempPx(it1+1)-IsoData%ptsTempPx(it1))
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
  FxrTri%itTemp = 0; FxrTri%wtTemp = 0.;
  nFxr = FxrTri%nFxrTri;
  !$OMP PARALLEL PRIVATE(i, j, niso, temp, it, id, it1, itemp, ntemp, k)
  !$OMP DO SCHEDULE(GUIDED)
  DO i = 1, nFxr
    niso = FxrTri%niso(i)
    temp = FxrTri%temp(i)
    it = temp
    k = FxrTri%AccNt(i);
    DO j = 1, niso
      k = k+1
      id = FxrTri%MapNucl(k)
      it1 = itempmap(it, id)
      itemp = IsoData%ptrTemp(id)
      ntemp = IsoData%ptrTemp(id+1)- itemp
      it1 = it1+itemp
      FxrTri%itTemp(i,j) = it1
      FxrTri%wtTemp(i,j) = 1.
      IF (ntemp .LE. 1) CYCLE
      IF (it1 .GE. ntemp+itemp) CYCLE
      FxrTri%wtTemp(i,j) = (-temp + IsoData%ptsTemp(it1+1))/(IsoData%ptsTemp(it1+1)-IsoData%ptsTemp(it1))
    END DO
    IF (ScatOrder < 1) CYCLE
    k = FxrTri%AccNt(i);
    DO j = 1, niso
      k = k+1;
      id = FxrTri%MapNuclp13(k)
      it1 = itempmapp1(it, id)
      itemp = IsoData%ptrTempPx(id)
      ntemp = IsoData%ptrTempPx(id+1)- itemp
      it1 = it1+itemp
      FxrTri%itTempPx(i,j) = it1
      FxrTri%wtTempPx(i,j) = 1.
      IF (ntemp .LE. 1) CYCLE
      IF (it1 .GE. ntemp+itemp) CYCLE
      FxrTri%wtTempPx(i,j) = (-temp + IsoData%ptsTempPx(it1+1))/(IsoData%ptsTempPx(it1+1)-IsoData%ptsTempPx(it1))
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
END SUBROUTINE

SUBROUTINE UpdateResPinFuelTemp
USE GEOM, ONLY : Core
USE PE_Mod, ONLY : PE
IMPLICIT NONE
INTEGER :: i, ipin
INTEGER :: nfxr, nxyz
REAL, ALLOCATABLE :: areasum(:), tempsum(:)
REAL :: area, temp
IF (.NOT. lrestrmt) RETURN

nfxr = FxrD%nFxrD; nxyz = Core%nxy*(PE%myze-PE%myzb+1)
ALLOCATE(areasum(nxyz), tempsum(nxyz))
areasum = 0; tempsum = 0;
DO i = 1, nfxr
  IF (.NOT. FxrD%lfuel(i)) CYCLE
  ipin = FxrD%ipin(i)
  area = FxrD%area(i); temp = FxrD%temp(i)
  areasum(ipin) = areasum(ipin) + area
  tempsum(ipin) = tempsum(ipin) + area*temp
END DO
nfxr = FxrTri%nFxrTri
DO i = 1, nfxr
  IF (.NOT.FxrTri%lres(i)) CYCLE
  ipin = FxrTri%ipin(i)
  IF (.NOT.GlobalPinInfo%lAIC(ipin)) CYCLE
  area = FxrTri%area(i); temp = FxrTri%temp(i)
  areasum(ipin) = areasum(ipin) + area
  tempsum(ipin) = tempsum(ipin) + area*temp
END DO
DO i = 1, nxyz
  IF (areasum(i).GT.0.) THEN
    ResPin%temp(i) = tempsum(i)/areasum(i)
  ELSE
    ResPin%temp(i) = 0.
  END IF
END DO

DEALLOCATE(areasum, tempsum)

END SUBROUTINE

SUBROUTINE UpdateResPinMapNucl
USE GEOM, ONLY : Core
USE PE_Mod, ONLY : PE
USE XSLIB_MOD, ONLY : mapnucl, mapnuclRes
IMPLICIT NONE
INTEGER :: i,iz,k,m,l,n,iso
INTEGER :: nPin, nxy, niso

nxy = core%nxy
nPin = nxy*(PE%myze-PE%myzb+1)
!$OMP PARALLEL PRIVATE(iso,m,k,l,n,iz,niso)
!$OMP DO SCHEDULE(GUIDED)
DO i = 1, nPin
  iz = (nPin-1)/nxy+PE%myzb
  IF (.NOT.core%lFuelPlane(iz)) CYCLE
  niso = core%ResVarPin(MOD(i-1,nxy)+1,iz)%niso
  ResPin%niso(i) = niso
  DO iso = 1, niso
    ResPin%mapnucl(iso,i) = mapnucl(ResPin%idiso(iso,i))
  END DO
  DO iso = 1, niso
    n = ResPin%idiso(iso,i)
    m = ResPin%mapnucl(iso,i)
    ResPin%ind(m,i) = 0.
    DO k = 1, niso
      l = ResPin%idiso(k,i)
      IF (mapnuclRes(n,iz).EQ.mapnuclRes(l,iz)) THEN
        ResPin%ind(m,i) = ResPin%ind(m,i)+ResPin%pnum(k,i)
      END IF
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
END SUBROUTINE

SUBROUTINE CalSigPotwNoReso(siglpD, siglpT, lskipD, lskipT)
USE Core_mod, ONLY : GroupInfo
IMPLICIT NONE
REAL :: siglpD(:,:),siglpT(:,:)
LOGICAL :: lskipD, lskipT
INTEGER :: i, ifxr, iso, isoidx, ig, gb, ge
REAL :: buf

gb = GroupInfo%nofg+1; ge = GroupInfo%nofg+GroupInfo%norg

IF (.NOT. lskipD) THEN
!$OMP PARALLEL PRIVATE(ifxr,iso,isoidx,ig,buf)
!$OMP DO SCHEDULE(GUIDED)
DO i = 1, FxrD%nResD
  ifxr = FxrD%mapR2G(i)
  DO ig = gb, ge
    buf = 0.;
    DO iso = 1, FxrD%niso(ifxr)
      isoidx = FxrD%mapnucl(iso,ifxr)
      IF (ResIsoData%IsoG2R(isoidx).NE.0) CYCLE
!      IF (i.EQ.2 .AND. ig.EQ.(gb)) print*, FxrD%pnum(iso,ifxr), ResIsoData%lamsigp(ig,isoidx)
      buf = buf+FxrD%pnum(iso,ifxr)*ResIsoData%lamsigp(ig,isoidx)
    END DO
    siglpD(ig-gb+1,i) = buf
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
END IF

IF (.NOT. lskipT) THEN
!$OMP PARALLEL PRIVATE(ifxr,iso,isoidx,ig,buf)
!$OMP DO SCHEDULE(GUIDED)
DO i = 1, FxrTri%nResTri
  ifxr = FxrTri%mapR2G(i)
  DO ig = gb, ge
    buf= 0.
    DO iso = FxrTri%AccNt(ifxr)+1, FxrTri%AccNt(ifxr+1)
      isoidx = FxrTri%mapnucl(iso)
      IF (isoidx.EQ.0) CYCLE
      IF (ResIsoData%IsoG2R(isoidx).NE.0) CYCLE
!      IF (FxrTri%lCLD(ifxr) .AND. ig.EQ.gb) print*, i*10+iso-FxrTri%AccNt(ifxr), FxrTri%pnum(iso), ResIsoData%lamsigp(ig,isoidx)
      buf = buf+FxrTri%pnum(iso)*ResIsoData%lamsigp(ig,isoidx)
    END DO
    siglpT(ig-gb+1,i) = buf
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
END IF

END SUBROUTINE

END MODULE

