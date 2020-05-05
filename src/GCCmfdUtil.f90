#include <defines.h>
SUBROUTINE Init_GCCmfd()
USE PARAM
USE TYPEDEF,       ONLY : GroupInfo_Type
USE Core_mod,      ONLY : GroupInfo,            GcGroupInfo
USE GCCmfd_mod,    ONLY : PredefineGrpStruct,   UserdefineGrpStruct1,   UserdefineGrpStruct2
USE CNTL,          ONLY : nTracerCntl
IMPLICIT NONE

GroupInfo%nGC = nTracerCntl%nGC
SELECT CASE(nTracerCntl%GcStructMod) 
  CASE(0)
    CALL PreDefineGrpStruct(GroupInfo, nTracerCntl%nGC)
  CASE(1)
    CALL UserDefineGrpStruct1(GroupInfo, nTracerCntl%nGC, nTracerCntl%ELB(1:nTracerCntl%nGC-1))
  CASE(2)
    CALL UserDefineGrpStruct2(GroupInfo, nTracerCntl%nGC, nTracerCntl%GLB(1:nTracerCntl%nGC-1))
END SELECT

CALL SetGcGroupInfo(GcGroupInfo, GroupInfo)
IF(GroupInfo%nGC .NE. 2) THEN
  CALL AllocGCCMFD()
ELSE
  !CALL AllocGCCMFD()
  CALL AllocCMFD2G()
ENDIF
CONTINUE
END SUBROUTINE

SUBROUTINE AllocGCCMFD()
USE PARAM
USE TYPEDEF,         ONLY : CmInfo_Type,     CmfdLs_Type,         PE_TYPE
USE Geom,            ONLY : Core,            ncbd
USE Core_mod,        ONLY : GcPinXs,          GcCMFDLS,           GcPsiC,         GcPsiCD,    &
                            GcPhiC,           GroupInfo,          GcGroupInfo,    CmInfo
USE CMFD_Mod,        ONLY : PinNeighIdx,      AllocPinXs
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

ALLOCATE(GcPinXS(nxy,myzbf-1:myzef+1))
ALLOCATE(GcCMFDLs(nGC))
DO ig = 1, ngc
  CALL Dmalloc0(GcCMFDLs(ig)%diag, 1, nxy, myzbf, myzef)
  CALL Dmalloc0(GcCMFDLs(ig)%RadOffDiag, 1, ncbd, 1, nxy, myzbf, myzef)
  CALL Dmalloc0(GcCMFDLs(ig)%AxOffDiag, 1, 2, 1, nxy, myzbf, myzef)
  GcCMFDLS(ig)%NeighIdx => PinNeighIdx
  GcCMFDLS(ig)%myzbf = myzbf; GcCMFDLS(ig)%myzef = myzef
  GcCMFDLS(ig)%nxy = nxy; GcCMFDLS(ig)%nbd = ncbd
  CALL Dmalloc0(GcCMFDLs(ig)%AxialPlaneMap, myzbf, myzef)
  DO iz = myzbf, myzef
    GcCMFDLs(ig)%AxialPlaneMap(iz) = iz
  ENDDO
ENDDO
CALL SetCmfdMpiOmpEnv(Core, GcCMFDLS, ngc, PE)

!InScatRange(1, 1:ngc) = 1;InScatRange(2, 1:ngc) = ngc;

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
CmInfo%GcPinXS => GcPinXS; CmInfo%GcCmfdLs => GcCmfdLs
CmInfo%GcPhiC => GcPhiC; 
CmInfo%GcPsiC => GcPsiC; CmInfo%GcPsiCD => GcPsiCD
END SUBROUTINE

SUBROUTINE SetGcGroupInfo(GcGrpInfo, GrpInfo)
USE PARAM
USE TYPEDEF,    ONLY : GroupInfo_Type
USE ALLOCS
IMPLICIT NONE
TYPE(GroupInfo_Type) :: GrpInfo, GcGrpInfo
!Make InScatRange
INTEGER :: ng, ngc
INTEGER :: ig, ib, ie, igc, ibc, iec
ngc = GrpInfo%ngc; ng = GrpInfo%ng
GcGrpInfo%ng = ngc
!In&Out Scat Range
CALL Dmalloc(GcGrpInfo%InScatRange, 2, ngc)
CALL Dmalloc(GcGrpInfo%OutScatRange, 2, ngc)
GcGrpInfo%InScatRange(1, :) =  ngc; GcGrpInfo%InScatRange(2, :) = 1
GcGrpInfo%OutScatRange(1, :) =  ngc; GcGrpInfo%OutScatRange(2, :) = 1
DO ig = 1, ng
  igc = GrpInfo%InvGCStruct(ig)
  
  ib = GrpInfo%InScatRange(1, ig); ie = GrpInfo%InScatRange(2, ig);  
  ibc = GrpInfo%InvGCStruct(ib); iec = GrpInfo%InvGCStruct(ie);
  GcGrpInfo%InScatRange(1, igc) = MIN(GcGrpInfo%InScatRange(1, igc), ibc)
  GcGrpInfo%InScatRange(2, igc) = MAX(GcGrpInfo%InScatRange(2, igc), iec)

  ib = GrpInfo%OutScatRange(1, ig); ie = GrpInfo%OutScatRange(2, ig); 
  ibc = GrpInfo%InvGCStruct(ib); iec = GrpInfo%InvGCStruct(ie);
  GcGrpInfo%OutScatRange(1, igc) = MIN(GcGrpInfo%OutScatRange(1, igc), ibc)
  GcGrpInfo%OutScatRange(2, igc) = MAX(GcGrpInfo%OutscatRange(2, igc), iec)
ENDDO

!Upscattering Range
GcGrpInfo%lUpScat = .FALSE.
DO ig = 1, ngc
  IF(GcGrpInfo%InScatRange(2,ig) .GT. ig) THEN
    GcGrpInfo%lUpScat = .TRUE.  
    EXIT
  ENDIF
ENDDO
IF(GcGrpInfo%lUpScat) THEN
  GcGrpInfo%UpscatRange(1) = ngc;GcGrpInfo%UpscatRange(2) = 1
  DO ig = 1, ngc
    IF(GcGrpInfo%InScatRange(2,ig) .GT. ig) THEN
      GcGrpInfo%UpscatRange(1) = MIN(GcGrpInfo%UpscatRange(1), ig)
      GcGrpInfo%UpscatRange(2) = MAX(GcGrpInfo%UpscatRange(2), ig)
    ENDIF
  ENDDO
   GcGrpInfo%UpscatRange(2) = ngc
ENDIF

END SUBROUTINE



SUBROUTINE GenGcHomoXs(Core, PinXS, GcPinXS, phi, GcPhi, ng, ngc, GroupInfo, GcGroupInfo, PE)
USE PARAM
USE TYPEDEF,          ONLY : CoreInfo_Type,         PinXS_Type,         GroupInfo_Type,       &
                             PE_TYPE
USE GcCmfd_mod,       ONLY : MakeGcPinXS
USE BasicOperation,   ONLY : CP_VA 
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(PinXS_Type), POINTER :: PinXS(:, :), GcPinXS(:, :)
TYPE(GroupInfo_Type) :: GroupInfo,     GcGroupInfo
TYPE(PE_TYPE) :: PE
REAL, POINTER :: Phi(:, :, :), GcPhi(:, :, :)
INTEGER :: ng, ngc

REAL :: PHI0(ng)
INTEGER, POINTER :: SubPlaneMap(:)
INTEGER :: nxy, myzb, myze, myzbf, myzef
INTEGER :: ixy, iz, iz0, ig

nxy = Core%nxy
myzb = PE%myzb; myze = PE%myze
myzbf = PE%myzbf; myzef = PE%myzef

SubPlaneMap => Core%SubPlaneMap

DO iz = myzbf, myzef
  iz0 = SubPlaneMap(iz)
  DO ixy = 1, nxy
    CALL CP_VA(Phi0(1:ng), Phi(ixy, iz, 1:ng), ng)
    CALL MakeGcPinXS(PinXS(ixy, iz0), GcPinXS(ixy, iz), Phi0(1:ng), ng, ngc, GroupInfo, GcGroupInfo)
    GcPhi(ixy, iz, 1:ngc) = GcPinXS(ixy, iz)%phi(1:ngc)
  ENDDO
ENDDO

NULLIFY(SubPlaneMap)
END SUBROUTINE

SUBROUTINE MakeGcPinXS(PinXS, GcPinXS, phi, ng, ngc, GroupInfo, GcGroupInfo)
USE PARAM
USE TYPEDEF,    ONLY : PinXS_Type,     GroupInfo_Type
IMPLICIT NONE
TYPE(PinXS_Type) :: PinXS, GcPinXS
TYPE(GroupInfo_Type) :: GroupInfo,     GcGroupInfo
REAL :: Phi(ng)
INTEGER :: ng, ngc
REAL :: RR(5), RRS(ngc, ngc)
REAL :: sumphi, sumchi
INTEGER :: igc, ig, igb, ige
INTEGER :: igs, igsb, igse, ig0

DO igc = 1, ngc
  igb = GroupInfo%GCStruct(1, igc); ige = GroupInfo%GCStruct(2, igc)
  RR(1:5) = 0; sumphi =0 ; sumchi = 0
  DO ig = igb, ige
    RR(1) = RR(1) + PinXS%XST(ig) * phi(ig); RR(2) = RR(2) + PinXS%XSTR(ig) * phi(ig)
    RR(3) = RR(3) + PinXS%XSNF(ig) * phi(ig); RR(4) = RR(4) + PinXS%XSKF(ig) * phi(ig)
    RR(5) = RR(5) + PinXS%XSD(ig) * phi(ig);
    sumphi = sumphi + phi(ig)
    sumchi = sumchi + PinXS%CHI(ig)
  ENDDO
  RR = RR / sumphi
  GcPinXS%XST(igc) = RR(1); GcPinXS%XSTR(igc) = RR(2);
  GcPinXS%XSNF(igc) = RR(3); GcPinXS%XSKF(igc) = RR(4);
  GcPinXS%XSD(igc) = RR(5)
  GcPinXS%CHI(igc) = SumChi; GcPinXS%phi(igc) = sumphi  
ENDDO

!Scattering Matrix
RRS = 0; sumphi = 0
DO ig = 1, ng
  igc = GroupInfo%InvGCStruct(ig)
  igsb = PinXS%XSS(ig)%ib; igse = PinXS%XSS(ig)%ie
  ig0 = igc
  RRS(ig0, igc) = RRS(ig0, igc) + PinXS%XSS(ig)%WithInGroupScat * PHI(ig)  !ig0 -> igc
  DO igs= igsb, igse
    ig0 = GroupInfo%InvGCStruct(igs)
    !igs -> ig,  ig0 -> igc
    RRS(ig0, igc) = RRS(ig0, igc) + PinXS%XSS(ig)%from(igs) * Phi(igs)
  ENDDO
ENDDO

DO igc = 1, ngc
  RRS(igc, :) = RRS(igc, :) / GcPinXS%phi(igc) !PHI(igc)
ENDDO

DO ig = 1, ngc
  igsb = GcPinXS%Xss(ig)%ib; igse = GcPinXS%Xss(ig)%ie;
  DO igs =igsb, igse
    GcPinXS%XSS(ig)%from(igs) = RRS(igs, ig)
  ENDDO
  GcPinXS%XSS(ig)%WithInGroupScat = GcPinXS%XSS(ig)%from(ig)
  GcPinXS%XSS(ig)%from(ig) = 0
ENDDO

DO ig = 1, ngc
 !GcPinXS%XSD(ig) = 1._8 / (3._8 * GcPinXS%XSTR(ig))
  GcPinXS%XSR(ig) = GcPinXS%XSTR(ig) - GcPinXS%XSS(ig)%WithInGroupScat
ENDDO

END SUBROUTINE

SUBROUTINE GcRadCouplingCoeff(Core, PinXS, GcPinXS, Phi, ng, ngc, GroupInfo, PE)
USE PARAM
USE TYPEDEF,   ONLY : CoreInfo_Type,          PinXS_Type,          GroupInfo_Type,    &
                      Pin_Type,               PE_TYPE
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(PinXS_Type), POINTER :: PinXs(:, :), GcPinXs(:, :)
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_TYPE) :: PE
REAL, POINTER :: phi(:, :, :)
INTEGER :: ng, ngc

TYPE(PinXS_Type), POINTER :: myPinXS, neighPinXS
TYPE(PinXS_Type), POINTER :: myGcPinXS, neighGcPinXS
TYPE(Pin_Type), POINTER :: Pin(:)
INTEGER, POINTER :: SubPlaneMap(:)
INTEGER :: nxy, nbd, myzb, myze, myzbf, myzef
INTEGER :: ixy, iz, iz0, ineigh, ibd, inbd, ig, igc
REAL :: smy, Dtil, Dhat, myphi, neighPhi, phisum, mybeta, neighbeta
REAL :: Jnet(ngc), jfdm

nxy = Core%nxy; nbd = 4
myzb = PE%myzb; myze = PE%myze
myzbf = PE%myzbf; myzef = PE%myzef
Pin => Core%Pin; SubPlaneMap => Core%SubPlaneMap

DO iz = myzbf, myzef
  iz0 = SubPlaneMap(iz)
  DO ixy = 1, nxy
    myGcPinXS => GcPinXs(ixy, iz); myPinXS => PinXs(ixy, iz0)
    DO ibd = 1, nbd
      ineigh = Pin(ixy)%NeighIdx(ibd); inbd = Pin(ixy)%NeighSurfIdx(ibd)    !The Corresponding Surface Index of 
      smy = Pin(ixy)%BdLength(ibd)
      Jnet =0
      DO ig = 1, ng
        myPhi = Phi(ixy, iz, ig); neighPhi = 0
        IF(ineigh .GT. 0) neighPhi = Phi(ineigh, iz, ig)
        igc = GroupInfo%InvGCStruct(ig)
        Dhat = PinXS(ixy, iz0)%Dhat(ibd, ig); Dtil = PinXS(ixy, iz0)%Dtil(ibd, ig)
        Jnet(igc) = Jnet(igc) + (Dtil - Dhat) * MyPhi - (Dtil + Dhat) * NeighPhi
      ENDDO
      IF(ineigh .GT. 0) THEN
        DO ig = 1, ngc
          myPhi = GcPinXS(ixy, iz)%phi(ig); neighPhi = GcPinXS(ineigh, iz)%phi(ig)
          phisum = neighphi + myphi
          mybeta = GcPinXS(ixy, iz)%XSD(ig) / Pin(ixy)%Center2SurfaceL(ibd)
          neighbeta = GcPinXS(ineigh, iz)%XSD(ig) / Pin(ineigh)%Center2SurfaceL(inbd)
          Dtil = mybeta * neighbeta/(mybeta + neighbeta) * smy
          jfdm = Dtil * (myphi - neighPhi)
          Dhat = (jfdm - jnet(ig)) / Phisum
          GcPinXS(ixy, iz)%Dhat(ibd, ig) = Dhat; GcPinXS(ixy, iz)%Dtil(ibd, ig) = Dtil
        ENDDO
      ELSE
        IF(ineigh .EQ. VoidCell) THEN
          neighbeta = 0.5_8
        ELSE
          neighbeta = 0
        ENDIF
        DO ig = 1, ngc
          myPhi = GcPinXS(ixy, iz)%phi(ig); neighPhi = 0
          phisum = neighphi + myphi
          mybeta = GcPinXS(ixy, iz)%XSD(ig) / Pin(ixy)%Center2SurfaceL(ibd)
          Dtil = mybeta * neighbeta/(mybeta + neighbeta) * smy
          jfdm = Dtil * (myphi - neighphi)
          dhat = (jfdm - jnet(ig)) / phisum
          GcPinXS(ixy, iz)%Dhat(ibd, ig) = Dhat; GcPinXS(ixy, iz)%Dtil(ibd, ig) = Dtil
        ENDDO
      ENDIF  
      
    ENDDO
  ENDDO
ENDDO
NULLIFY(PIN); NULLIFY(SUbPlaneMap)

END SUBROUTINE

SUBROUTINE GcAxCouplingCoeff(Core, GcPinXS, PHI, GcPhi, AxDtil, AxDhat, ng, ngc, GroupInfo, PE)
USE PARAM
USE TYPEDEF,        ONLY : CoreInfo_Type,          AxFlx_Type,          PinXS_Type,       &
                           GroupInfo_Type,         PE_TYPE
#ifdef MPI_ENV
USE MPIComm_MOD, ONLY : REDUCE, GetNeighDat
#endif
USE BasicOperation, ONLY : CP_CA
IMPLICIT NONE
TYPE(CoreInfo_Type) :: Core
TYPE(PinXs_Type), POINTER :: GcPinXs(:, :)
TYPE(GroupInfo_Type) :: GroupInfo
TYPE(PE_Type) :: PE
REAL, POINTER :: PHI(:, :, :), GcPhi(:, :, :)
REAL, POINTER :: AxDtil(:, :, :, :), AxDhat(:, :, :, :)
INTEGER :: ng, ngc

INTEGER, POINTER :: SubPlaneMap(:)
REAL, POINTER :: hzfm(:)
REAL, POINTER :: neighXSD(:, :)
REAL :: Dtil, Dhat, mybeta, neighbeta, myphi, neighphi, phisum
REAL :: jnet(ngc), jfdm
REAL :: pDhat, alpha

INTEGER :: mp(1:2)
INTEGER :: nxy, nzfm, myzb, myze, myzbf, myzef
INTEGER :: i, iz, iz0, ixy, ineigh
INTEGER :: ig, igc
mp = (/-1,1/)

nxy = Core%nxy;nzfm = Core%nzfm
myzb =PE%myzb;   myze = PE%myze
myzbf =PE%myzbf; myzef = PE%myzef
SubPlaneMap => Core%subPlaneMap; hzfm => Core%hzfm

#ifdef MPI_ENV
IF(PE%nCMFDproc .GT. 1) THEN 
  DO ig = 1, ngc
      CALL GetNeighDat(Phi(1:nxy, myzbf-1:myzef+1, ig), nxy, myzbf, myzef,   &
                             PE%myCmfdRank, PE%nproc, PE%MPI_CMFD_COMM)
      CALL GetNeighDat(GcPhi(1:nxy, myzbf-1:myzef+1, ig), nxy, myzbf, myzef,   &
                             PE%myCmfdRank, PE%nproc, PE%MPI_CMFD_COMM)                             
  ENDDO

  ALLOCATE(neighXSD(nxy, 4))
  DO ig = 1, ngc
    CALL CP_CA(neighXSD(1:nxy, 1:4), 0._8, nxy, 4)  
    DO ixy = 1, nxy
      neighXSD(ixy, 2) = GcPinXS(ixy,  myzbf)%XSD(ig)
      neighXSD(ixy, 3) = GcPinXS(ixy,  myzef)%XSD(ig)
    ENDDO
    
    CALL GetNeighDat(neighXSD(1:nxy, 1:4), nxy, 2, 3, &
                     PE%myCmfdRank, PE%nproc, PE%MPI_CMFD_COMM)
    DO ixy = 1, nxy
      GcPinXS(ixy,  myzbf-1)%XSD(ig) = neighXSD(ixy, 1)
      GcPinXS(ixy,  myzef+1)%XSD(ig) = neighXSD(ixy, 4)
    ENDDO
  ENDDO
  DEALLOCATE(neighXSD)
ENDIF
#endif

DO ixy = 1, nxy
  DO iz = myzbf, myzef
    DO i = 1, 2
      ineigh = iz + mp(i)
      neighbeta = 0
      IF(ineigh .EQ. 0 .AND. Core%AxBc(i)  .EQ. VoidCell) neighbeta =0.5_8
      IF(ineigh .EQ. nzfm .AND. Core%AxBc(i)  .EQ. VoidCell) neighbeta =0.5_8
      !Get
      Jnet = 0
      DO ig = 1, ng
        igc = GroupInfo%InvGCStruct(ig)
        MyPhi = Phi(ixy, iz, ig); NeighPhi = Phi(ixy, ineigh, ig)
        Dtil = AxDtil(i, ixy ,iz, ig); Dhat = AxDhat(i, ixy, iz, ig)
        !Dtil = AxFlx(iz, ixy)%Dtil(i, ig); Dhat = AxFlx(iz, ixy)%Dhat(i, ig)
        Jnet(igc) = Jnet(igc) + (Dtil - Dhat) * MyPhi - (Dtil + Dhat) * NeighPhi
      ENDDO
      
      DO ig = 1, ngc
        mybeta = GcPinXS(ixy, iz)%XSD(ig) / hzfm(iz)
        MyPhi = GcPhi(ixy, iz, ig); NeighPhi = GcPhi(ixy, ineigh, ig)
        IF(ineigh .NE. 0 .AND. ineigh .NE. nzfm+1) THEN
          neighbeta = GcPinXS(ixy, ineigh)%XSD(ig) / hzfm(ineigh)
        ENDIF
        phisum = myphi + neighphi
        Dtil = 2._8 * neighbeta * mybeta / (mybeta + neighbeta)
        jfdm = Dtil * (myphi - neighphi)
        pdhat = GcPinXS(ixy, iz)%AxDhat(i, ig)
        dhat = (jfdm - jnet(ig)) / phisum
!        IF(abs(dhat-pdhat) .GT. 10._8*Dtil) THEN
!          ALPHA = DTIL/(DTIL-ABS(DHAT-pdhat))
!          DHAT = DHAT * ALPHA + Pdhat*(1.-ALPHA)
!        ENDIF
        GcPinXS(ixy, iz)%AxDtil(i, ig) = Dtil
        GcPinXS(ixy, iz)%AxDhat(i, ig) = Dhat  
      ENDDO
    ENDDO
  ENDDO
ENDDO

NULLIFY(SubPlaneMap, hzfm)

END SUBROUTINE

SUBROUTINE GcCmfdPsiUpdt(PhiC, PsiC)
USE PARAM
USE CMFD_MOD,    ONLY : nxy,         myzbf,       myzef, &
                        SubPlaneMap, PinVolFm
USE GcCMfd_mod,  ONLY : ngc,         GcPinXS
IMPLICIT NONE
REAL, POINTER :: PhiC(:, :, :), PsiC(:, :)
INTEGER :: ig, iz, iz0, ixy
DO iz = myzbf, myzef
  iz0 = SubPlaneMap(iz)
  DO ixy = 1, nxy
    PsiC(ixy, iz) = 0
    DO ig = 1, ngc
      PsiC(ixy, iz) = PsiC(ixy, iz) + GcPinXS(ixy, iz)%xsnf(ig) * PhiC(ixy, iz, ig)
    ENDDO
    PsiC(ixy, iz) = PsiC(ixy, iz) * PinVolFm(ixy, iz0)
  ENDDO
ENDDO
END SUBROUTINE


SUBROUTINE GcCmfdSrcUpdt(SRC, psi, phi, eigv, ig)
USE PARAM
USE TYPEDEF,    ONLY : scatmat
USE CMFD_MOD,   ONLY : nxy,                                   &
                       myzb,        myze,         myzbf,       myzef,      &
                       SubPlaneMap, PinVolFm
USE GcCMFD_MOD, ONLY : GcPinXS
IMPLICIT NONE
REAL, POINTER :: SRC(:, :), psi(:, :), phi(:, :, :)
REAL :: Eigv
TYPE(scatmat), POINTER :: XSS
REAL :: reigv, SS
INTEGER :: ig

INTEGER :: iz, iz0, ixy, ig2, igb, ige
reigv = one/eigv

DO iz = myzbf, myzef
  iz0 = SubPlaneMap(iz)
  DO ixy = 1, nxy
    SRC(ixy, iz) = psi(ixy, iz) * reigv * GcPinXs(ixy, iz)%Chi(ig)
    XSS => GcPinXS(ixy, iz)%Xss(ig) 
    igb = XSS%ib; ige = XSS%ie
    SS = 0
    DO ig2 = igb, ige
      SS = SS + phi(ixy, iz, ig2) * XSS%From(ig2)
    ENDDO
    SRC(ixy, iz) = SRC(ixy, iz) + SS * PinVolFm(ixy, iz0)
  ENDDO
ENDDO
NULLIFY(XSS)

END SUBROUTINE

SUBROUTINE GcCmfdEigUpdt(psi, psid, eigv, psierr, PE)
USE PARAM
USE TYPEDEF,  ONLY : PE_TYPE
USE CMFD_MOD, ONLY : nxy,           myzbf,       myzef,            &
                     SubPlaneMap
#ifdef MPI_ENV
USE MPIComm_mod, ONLY : REDUCE
#endif
IMPLICIT NONE
TYPE(PE_TYPE), OPTIONAL :: PE
REAL, POINTER :: psi(:, :), psid(:, :)
REAL :: psierr, eigv
REAL :: psipsi, psipsid
INTEGER :: ig, iz, iz0, ixy
#ifdef MPI_ENV
REAL :: buf(3), buf0(3) 
#endif
psipsi = zero; psipsid = zero; psierr = zero
DO iz = myzbf, myzef
  DO ixy = 1, nxy
    psipsi = psipsi + Psi(ixy, iz) * Psi(ixy, iz)
    psipsid = psipsid + Psi(ixy, iz) * PsiD(ixy, iz)
    psierr = psierr + (Psi(ixy, iz) - PsiD(ixy, iz))**2
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

FUNCTION GcResidualError(phi, psi, eigv,  PE, constsrc)
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
REAL :: GcResidualError
REAL, POINTER :: phi(:, :, :), psi(:, :)
REAL :: eigv
TYPE(PE_TYPE) :: PE
REAL, OPTIONAL :: ConstSrc
INTEGER :: ig, ig0, iz, iz0, ixy, ibd, ineigh, nbd
REAL :: vol, LMNT, tsrc, area
REAL :: myphi, neighphi, jsum, jnet, dtil, dhat
#ifdef MPI_ENV
REAL :: buf(2), buf0(2)
#endif
GcResidualError = 0; tsrc =0
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
    
    SRC(ixy, iz) = Psi(ixy, iz) / Eigv / PinVolFm(ixy, iz0)
    DO ig = 1, ngc
      DO ig0 = GcPinXS(ixy, iz)%XSS(ig)%ib, GcPinXS(ixy, iz)%XSS(ig)%ie
        Src(ixy, iz) = Src(ixy, iz) + phi(ixy, iz, ig0) * GcPinXS(ixy, iz)%XSS(ig)%From(ig0)
      ENDDO
    ENDDO
    IF(PRESENT(ConstSrc)) SRC(ixy, iz) = SRC(ixy, iz) + ConstSrc * ngc
    SRC(ixy, iz) = Src(ixy, iz) * PinVolFm(ixy, iz0)
    lmnt = lmnt - SRC(ixy, iz) 
    tsrc = tsrc +  src(ixy, iz) * src(ixy, iz)
    GcResidualError = GcResidualError + lmnt * lmnt    
  ENDDO
ENDDO
#else

#endif
#ifdef MPI_ENV 
buf0  = (/GcResidualError, tsrc/)
CALL REDUCE(buf0, buf, 2, PE%MPI_CMFD_COMM, .TRUE.)
GcResidualError = buf(1); tsrc = buf(2)
#endif
GcResidualError = SQRT(GcResidualError/tsrc)

END FUNCTION

SUBROUTINE MgCmfdSolUpdt(CmInfo, GroupInfo, w)
USE PARAM
USE TYPEDEF,   ONLY : CmInfo_Type, PE_TYPE,    GroupInfo_Type
USE CMFD_MOD,  ONLY : nxy,         myzbf,       myzef,               &
                      src,                                           &
                      SubPlaneMap, PinVol,      PinVolFm,            &
                      hzfm,        PinNeighIdx, ng
USE GcCmfd_Mod, ONLY : ngc,       GcPinXS,     GcCmfdLs
IMPLICIT NONE

TYPE(CmInfo_Type) :: CmInfo
TYPE(GroupInfo_Type) :: GroupInfo
REAL :: w 

INTEGER :: ixy, iz, iz0, ig, ig0, ig1, ig2
REAL :: fmult, temp, w2

REAL, POINTER :: GcPhiC(:, :, :), PhiFm(:, :, :)

GcPhiC => CMInfo%GcPhiC; PhiFm => CmInfo%PhiFm
w2=0.0
DO iz = myzbf, myzef
  DO ixy = 1, nxy
    DO ig = 1, ng
      ig0 = GroupInfo%InvGCStruct(ig)
!      ig1 = GroupInfo%GcStruct(1, ig);ig2 = GroupInfo%GcStruct(2, ig)
      fmult = GcPhiC(ixy, iz, ig0) / GcPinXs(ixy,iz)%phi(ig0)
      !temp= PhiFm(ixy, iz, ig) * fmult
      PhiFm(ixy, iz, ig) = PhiFm(ixy, iz, ig) * fmult !* w2 + temp * (1-w2)
    ENDDO
  ENDDO
ENDDO

NULLIFY(GcPhiC, PhiFm)

END SUBROUTINE