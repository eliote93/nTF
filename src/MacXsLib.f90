#include <defines.h>
#include <DefDBG.h>
MODULE MacXsLib_Mod
USE PARAM
USE TYPEDEF,         ONLY : XsMac_Type,     GroupInfo_Type,      Fxrinfo_type
USE XSLIB_MOD
USE XsUtil_mod,      ONLY : AllocXsMac,          AllocMacIsoXs,  &
                            LineIntPol, LineIntPol2,   CalcWgt,   CalcWgt2,   &!TempIntPol_RIA,
                            XsTempInterpolation, P1XsTempInterpolation,            &
                            P2XsTempInterpolation, P3XsTempInterpolation,     &
                            AllocPxsMac, LineIntPol2_2Var!,  AllocPxsSMac ! for OTF RIF CALCULATION...
USE BasicOperation,  ONLY : CP_CA, CP_VA,   MULTI_CA,            DotProduct
USE Allocs
IMPLICIT NONE
INTEGER :: nisomax = 1000
LOGICAL :: lalloc = .FALSE.
!INTEGER, PRIVATE :: i, j, k
!REAL, POINTER :: IsoXsMacA(:, :), ISoXsMacnf(:, :), IsoXsMacS0(:, :), IsoXsMacTr(:, :),   &
!                 IsoXsMacf(:, :), IsoXsMackf(:, :)
TYPE(XsMac_Type) :: XsMac0
INTERFACE MacXsBase
  MODULE PROCEDURE MacXsBase_Fxr
  MODULE PROCEDURE MacXsBase_CrCsp
  MODULE PROCEDURE MacXsBase_gen
END INTERFACE

INTERFACE BaseMacSTr
  MODULE PROCEDURE BaseMacSTr_Fxr
   MODULE PROCEDURE BaseMacSTr_Csp
  MODULE PROCEDURE BaseMacSTr_Gen
END INTERFACE

INTERFACE MacXsScatMatrix
  MODULE PROCEDURE MacXsScatMatrix_Fxr
  MODULE PROCEDURE MacXsScatMatrix_CSP
  MODULE PROCEDURE MacXsScatMatrix_Gen
END INTERFACE

INTERFACE MacP1XsScatMatrix
  MODULE PROCEDURE MacP1XsScatMatrix_Fxr
  MODULE PROCEDURE MacP1XsScatMatrix_Csp
  MODULE PROCEDURE MacP1XsScatMatrix_Gen
END INTERFACE

INTERFACE MacP2XsScatMatrix
  MODULE PROCEDURE MacP2XsScatMatrix_Fxr
  MODULE PROCEDURE MacP2XsScatMatrix_Csp
  MODULE PROCEDURE MacP2XsScatMatrix_Gen
END INTERFACE

INTERFACE MacP3XsScatMatrix
  MODULE PROCEDURE MacP3XsScatMatrix_Fxr
  MODULE PROCEDURE MacP3XsScatMatrix_Csp
  MODULE PROCEDURE MacP3XsScatMatrix_Gen
END INTERFACE

!INTERFACE MacP2XsScatMatrix
!  MODULE PROCEDURE MacP2XsScatMatrix_Fxr
!  MODULE PROCEDURE MacP2XsScatMatrix_Gen
!END INTERFACE

INTERFACE MacXsNf
  MODULE PROCEDURE MacXsNf_Fxr
  MODULE PROCEDURE MacXsNf_CrCsp
  MODULE PROCEDURE MacXsNf_Gen
END INTERFACE

INTERFACE MacXsKf
  MODULE PROCEDURE MacXsKf_Fxr
  MODULE PROCEDURE MacXsKf_CrCsp
  MODULE PROCEDURE MacXsKf_Gen
END INTERFACE



INTERFACE GetMacChi
  MODULE PROCEDURE GetMacChi_Fxr
  MODULE PROCEDURE GetMacChi_Gen
END INTERFACE


CONTAINS


!SUBROUTINE AllocIsotopicMacXs()
!IMPLICIT NONE
!INTEGER :: ng
!CALL Dmalloc(IsoXsMacA, nelthel, noghel)
!CALL Dmalloc(IsoXsMacnf, nelthel, noghel)
!CALL Dmalloc(IsoXsMacf, nelthel, noghel)
!CALL Dmalloc(IsoXsMackf, nelthel, noghel)
!CALL Dmalloc(IsoXsMacS0, nelthel, noghel)
!CALL Dmalloc(IsoXsMacTr, nelthel, noghel)
!lalloc = .TRUE.
!END SUBROUTINE

SUBROUTINE MacXsBase_Fxr(XsMac, FxrInfo, igb, ige, ng, eigv, kin, lIsoXsOut)
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
TYPE(Fxrinfo_type) :: FxrInfo
INTEGER :: igb, ige, ng
REAL :: eigv
LOGICAL :: kin, lIsoXsOut

REAL :: temp
INTEGER :: niso
INTEGER, POINTER :: idiso(:)
REAL, POINTER ::pnum(:)

niso = FxrInfo%niso; temp = FxrInfo%Temp
pnum => FxrInfo%pnum
idiso => FxrInfo%idiso
CALL MacXsBase_Gen(XsMac, temp, niso,  FxrInfo%idiso,  FxrInfo%pnum, igb, ige, ng, eigv, kin, lIsoXsOut)

NULLIFY(pnum, idiso)
END SUBROUTINE

SUBROUTINE MacXsBase_CrCsp(XsMac, FxrInfo, igb, ige, ng, eigv, kin, lIsoXsOut, lCrCspFtn)
USE XsUtil_mod,     ONLY : GetXsMacDat,        ReturnXsMacDat
USE CrCsp_Mod,      ONLY : MacXsBaseCsp
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
TYPE(Fxrinfo_type) :: FxrInfo
INTEGER :: igb, ige, ng
REAL :: eigv
LOGICAL :: kin, lIsoXsOut, lCrCspFtn

TYPE(XsMac_Type), POINTER :: XsMac1, XsMac2
REAL :: temp
INTEGER :: i
INTEGER :: niso
INTEGER, POINTER :: idiso(:)
REAL, POINTER ::pnum(:)

IF(lCrCspFtn .AND. FxrInfo%lCrCspFtn) THEN
  CALL GetXsMacDat(XsMac1, ng, lIsoXsout)
  Call GetXsMacDat(XsMac2, ng, lIsoXsOut)
  niso = FxrInfo%CspFxr%niso(1); Temp = FxrInfo%Temp
  Pnum => FxrInfo%CspFxr%pnum(:, 1);   idiso => FxrInfo%CspFxr%isolist(:,1)
  CALL MacXsBase_Gen(XsMac1, temp, niso, idiso, pnum, igb, ige, ng, eigv, kin, lIsoXsOut)
  niso = FxrInfo%CspFxr%niso(2)
  Pnum => FxrInfo%CspFxr%pnum(:, 2);   idiso => FxrInfo%CspFxr%isolist(:, 2)
  CALL MacXsBase_Gen(XsMac2, temp, niso, idiso, pnum, igb, ige, ng, eigv, kin, lIsoXsOut)
  !Interpolation
  CALL MacXsBaseCsp(XsMac, XsMac1, XsMac2, FxrInfo%CspFxr, FxrInfo%niso, igb, ige, lIsoXsOut)
  CALL ReturnXsMacDat(XsMac1); CALL ReturnXsMacDat(XsMac2)
ELSE
  niso = FxrInfo%niso; temp = FxrInfo%Temp
  pnum => FxrInfo%pnum; idiso => FxrInfo%idiso
  CALL MacXsBase_Gen(XsMac, temp, niso, idiso, pnum, igb, ige, ng, eigv, kin, lIsoXsOut)
ENDIF

NULLIFY(pnum, idiso)
END SUBROUTINE

SUBROUTINE MacXsBase_Gen(XsMac, temp, niso, idiso, pnum, igb, ige, ng, eigv, kin, lIsoXsOut)
IMPLICIT NONE
TYPE(XsMac_Type), INTENT(INOUT) :: XsMac
INTEGER, INTENT(IN) :: niso, igb, ige, ng
INTEGER, INTENT(IN) :: idiso(niso)
REAL, INTENT(IN) :: pnum(niso)
REAL, INTENT(IN) :: eigv, temp
LOGICAL, INTENT(IN) :: kin, lIsoXsOut

TYPE(LIBDATA), POINTER :: isodata
INTEGER :: ig, iso, id, it1, it2, p1idx, it1p1, it2p1
REAL :: wt1, wt2, wt1p1, wt2p1, reigv, ND                        !Temperature interpolation
REAL, POINTER :: xsmact(:), xsmaca(:), xsmacs(:), xsmactr(:), xsmacstr(:)
REAL, POINTER :: xsmacf(:), xsmacnf(:), xsmackf(:)
REAL :: IsoXsMacT(niso, ng), IsoXsMacA(niso, ng), IsoXsMacS0(niso, ng), IsoXsMacS1(niso, ng), IsoXsMacSS(niso, ng),   &
        IsoXsMacTr(niso, ng), IsoXsMacf(niso, ng), ISoXsMacnf(niso, ng), IsoXsMackf(niso, ng)
REAL :: IsoXsRadCap(niso, ng)  !-- JSU EDIT 20170727

IF(.NOT. XsMac%lalloc) THEN
  XsMac%ng = ng
  CALL AllocXsMac(XsMac)
ENDIF
IF(lIsoXsOut .AND. .NOT. XsMac%lIsoAlloc) THEN
  CALL AllocMacIsoXs(XsMac, ng, nelthel)
  !CALL AllocMacIsoXs(XsMac, ng, niso)      ! 16/12/08 why not niso?  >> 17/01/13 r560 rollback, : niso increasing
ENDIF
IsoXsMacf(1:niso, igb:ige) = 0.
IsoXsMacNf(1:niso, igb:ige) = 0.
IsoXsMackf(1:niso, igb:ige) = 0.

xsmact  => XsMac%xsmact;  xsmaca  => XsMac%xsmaca;   xsmacs  => XsMac%xsmacs
xsmactr => XsMac%xsmactr; xsmacstr => XsMac%xsmacstr
xsmacf => XsMac%xsmacf;   xsmacnf => XsMac%xsmacnf;  xsmackf => XsMac%XsMacKf

reigv = 1._8/eigv

DO iso = 1, niso
  id = mapnucl(idiso(iso));   isodata => ldiso(id)
  ND = pnum(iso)
  p1idx = mapnuclp13(idiso(iso))
  !Temperature Interpolation
  CALL XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2)
  if (p1idx.eq.0) then
    wt1p1=0._8; wt2p1=0._8
    it1p1=1; it2p1=1
  else
    CALL p1XsTempInterpolation(p1idx, isodata, temp, wt1p1, wt2p1, it1p1, it2p1)
  endif
  DO ig = igb, ige
    IsoXsMacA(iso, ig) = ND * (wt2 * isodata%siga(ig, it2) + wt1 * isodata%siga(ig, it1))
    IsoXsMacS0(iso, ig) =  ND * (wt2 * isodata%sigs(ig, it2) + wt1 * isodata%sigs(ig, it1))
    IsoXsMacS1(iso, ig) =  ND * (wt2p1 * isodata%sigsp1(ig, it2p1) + wt1p1 * isodata%sigsp1(ig, it1p1))
    IsoXsMacSS(iso, ig) =  ND * (wt2 * isodata%sigss(ig, it2) + wt1 * isodata%sigss(ig, it1))
    IsoXsMactr(iso, ig) =  ND * (wt2 * isodata%sigtr(ig, it2) + wt1 * isodata%sigtr(ig, it1))
    IsoXsMacT(iso, ig) = IsoXsMacA(iso, ig) + IsoXsMacS0(iso, ig) !--- BYS edit : total XS in isotope-wise
    IsoXsRadCap(iso, ig) = IsoXsMacA(iso, ig) - IsoXsMacf(iso, ig) !-- JSU EDIT 20170727
    IF (isodata%ifis.eq.0) CYCLE
    IsoXsMacf(iso, ig) = ND * (wt2 * isodata%sigf(ig, it2) + wt1 * isodata%sigf(ig, it1))
    IsoXsMacnf(iso, ig) = ND * (wt2 * isodata%signf(ig, it2) + wt1 * isodata%signf(ig, it1))
    IsoXsMackf(iso, ig) = IsoXsMacf(iso, ig) * isodata%kappa
  ENDDO
  NULLIFY(isodata) !Free the pointing variable
ENDDO !Isotopic sweeep

!Summation
DO ig = igb, ige
  xsmaca(ig) = SUM(IsoXsMacA(1:niso, ig)); xsmacs(ig) = SUM(IsoXsMacS0(1:niso, ig))
  xsmactr(ig) = SUM(IsoXsMactr(1:niso, ig))
  xsmacstr(ig) = xsmactr(ig) - xsmaca(ig); xsmact(ig) = xsmaca(ig) + xsmacs(ig)
  xsmacf(ig) = SUM(IsoXsMacF(1:niso, ig)); xsmacnf(ig) = SUM(IsoXsMacNf(1:niso, ig)); xsmackf(ig) = SUM(IsoXsMacKf(1:niso, ig))
ENDDO

IF(lIsoXsOut) THEN
  XsMac%IsoXsMacT(1:niso, igb:ige) = IsoXsMacT(1:niso, igb:ige)
  XsMac%IsoXsMacA(1:niso, igb:ige) = IsoXsMacA(1:niso, igb:ige)
  XsMac%IsoXsMacS0(1:niso, igb:ige) = IsoXsMacS0(1:niso, igb:ige)
  XsMac%IsoXsMacS1(1:niso, igb:ige) = IsoXsMacS1(1:niso, igb:ige)
  XsMac%IsoXsMacSS(1:niso, igb:ige) = IsoXsMacSS(1:niso, igb:ige)
  XsMac%IsoXsMacTr(1:niso, igb:ige) = IsoXsMacTr(1:niso, igb:ige)
  XsMac%IsoXsMacf(1:niso, igb:ige) = IsoXsMacf(1:niso, igb:ige)
  XsMac%ISoXsMacnf(1:niso, igb:ige) = IsoXsMacNf(1:niso, igb:ige)
  XsMac%IsoXsMackf(1:niso, igb:ige) = IsoXsMacKf(1:niso, igb:ige)
  XsMac%IsoXsRadCap(1:niso, igb:ige) = IsoXsRadCap(1:niso, igb:ige)
ENDIF

IF(kin) xsmacnf(igb:ige) = xsmacnf(igb:ige)*reigv

NULLIFY(xsmact);  NULLIFY(xsmaca);   NULLIFY(xsmacs)
NULLIFY(xsmactr); NULLIFY(xsmacstr);
NULLIFY(xsmacf);  NULLIFY(xsmacnf);  NULLIFY(xsmackf)
END SUBROUTINE

SUBROUTINE MacXsAF(IsoXsMacA, IsoXsMacF, temp, niso, idiso, pnum, igb, ige, ng)
IMPLICIT NONE
REAL :: IsoXsMacA(niso, ng), IsoXsMacF(niso, ng)
INTEGER, INTENT(IN) :: niso, igb, ige, ng
INTEGER, INTENT(IN) :: idiso(niso)
REAL, INTENT(IN) :: pnum(niso)
REAL, INTENT(IN) :: temp

TYPE(LIBDATA), POINTER :: isodata
INTEGER :: ig, iso, id, it1, it2
REAL :: wt1, wt2, ND                        !Temperature interpolation

IsoXsMacf(1:niso, igb:ige) = 0.

DO iso = 1, niso
  id = mapnucl(idiso(iso));   isodata => ldiso(id)
  ND = pnum(iso)
  CALL XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2)
  DO ig = igb, ige
    IsoXsMacA(iso, ig) = ND * (wt2 * isodata%siga(ig, it2) + wt1 * isodata%siga(ig, it1))
    IF (isodata%ifis.eq.0) CYCLE
    IsoXsMacf(iso, ig) = ND * (wt2 * isodata%sigf(ig, it2) + wt1 * isodata%sigf(ig, it1))
  ENDDO
ENDDO !Isotopic sweeep
END SUBROUTINE

SUBROUTINE MacXsANF(XsMacA, XsMacNF, FxrInfo, igb, ige, ng)
IMPLICIT NONE
REAL :: XsMacA(ng), XsMacNF(ng)
TYPE(Fxrinfo_type) :: FxrInfo
INTEGER, INTENT(IN) :: igb, ige, ng

TYPE(LIBDATA), POINTER :: isodata
INTEGER :: ig, iso, id, it1, it2
REAL :: wt1, wt2, ND                        !Temperature interpolation

REAL :: temp
INTEGER :: niso
INTEGER, POINTER :: idiso(:)
REAL, POINTER ::pnum(:)

niso = FxrInfo%niso; temp = FxrInfo%Temp
pnum => FxrInfo%pnum
idiso => FxrInfo%idiso

XsMacA(igb:ige) = 0.
XsMacNF(igb:ige) = 0.

DO iso = 1, niso
  id = mapnucl(idiso(iso));   isodata => ldiso(id)
  ND = pnum(iso)
  CALL XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2)
  DO ig = igb, ige
    XsMacA(ig) = XsMacA(ig) + ND * (wt2 * isodata%siga(ig, it2) + wt1 * isodata%siga(ig, it1))
    IF (isodata%ifis.eq.0) CYCLE
    XsMacNF(ig) = XsMacNF(ig) + ND * (wt2 * isodata%signf(ig, it2) + wt1 * isodata%signf(ig, it1))
  ENDDO
ENDDO !Isotopic sweeep
END SUBROUTINE

SUBROUTINE MacXsAF_CrCsp(IsoXsMacA, IsoXsMacF, FxrInfo, igb, ige, ng, lCrCspFtn)
IMPLICIT NONE
REAL, POINTER :: IsoXsMacA(:,:), IsoXsMacF(:,:)
TYPE(Fxrinfo_type) :: FxrInfo
INTEGER :: igb, ige, ng
LOGICAL :: lCrCspFtn

REAL :: temp
INTEGER :: i, ig
INTEGER :: niso
INTEGER, POINTER :: idiso(:)
REAL, POINTER ::pnum(:)
REAL, POINTER :: IsoA1(:,:), IsoF1(:,:)
REAL :: wt1, wt2

IF(lCrCspFtn .AND. FxrInfo%lCrCspFtn) THEN
  niso = FxrInfo%CspFxr%niso(1); Temp = FxrInfo%Temp
  Pnum => FxrInfo%CspFxr%pnum(:, 1);   idiso => FxrInfo%CspFxr%isolist(:,1)
  CALL MacXsAF(IsoXsMacA, IsoXsMacF, temp, niso, idiso, pnum, igb, ige, ng)
  niso = FxrInfo%CspFxr%niso(2)
  Pnum => FxrInfo%CspFxr%pnum(:, 2);   idiso => FxrInfo%CspFxr%isolist(:, 2)
  ALLOCATE(IsoA1(niso,ng),IsoF1(niso,ng))
  CALL MacXsAF(IsoA1,IsoF1,temp,niso,idiso,pnum,igb,ige,ng)
  !Interpolation
  DO ig = igb, ige
    wt1 = FxrInfo%CspFXR%fcsp(ig,1); wt2 = FxrInfo%CspFXR%fcsp(ig,2)
    IsoXsMacA(:,ig) = wt1*IsoXsMacA(:,ig)+wt2*IsoA1(:,ig)
    IsoXsMacF(:,ig) = wt1*IsoXsMacF(:,ig)+wt2*IsoF1(:,ig)
  END DO
  DEALLOCATE(IsoA1, IsoF1)
ELSE
  niso = FxrInfo%niso; temp = FxrInfo%Temp
  pnum => FxrInfo%pnum; idiso => FxrInfo%idiso
  CALL MacXsAF(IsoXsMacA, IsoXsMacF, temp, niso, idiso, pnum, igb, ige, ng)
ENDIF

NULLIFY(pnum, idiso)
END SUBROUTINE

SUBROUTINE MacXsBase_Cusping(XsMac, iso, FxrInfo, igb, ige, ng, eigv, kin, lIsoXsOut)
USE Material_Mod, ONLY : Mixture
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
TYPE(Fxrinfo_type) :: FxrInfo
INTEGER :: igb, ige, ng
REAL :: eigv
LOGICAL :: kin, lIsoXsOut

REAL :: temp
INTEGER :: niso
INTEGER, POINTER :: idiso(:)
REAL, POINTER ::pnum(:)
INTEGER :: iso

niso = Mixture(iso)%niso; temp = FxrInfo%Temp
pnum => Mixture(iso)%pnum
idiso => Mixture(iso)%idiso

CALL MacXsBase_Gen(XsMac, temp, niso, idiso, pnum, igb, ige, ng, eigv, kin, lIsoXsOut)

NULLIFY(pnum, idiso)
END SUBROUTINE

SUBROUTINE MacXS1gBase(xsa, xsf, xsn2n, xsn3n, temp, niso, idiso, pnum, phi, ng)
IMPLICIT NONE

INTEGER :: niso, ng
REAL :: xsa(niso), xsf(niso), xsn2n(niso), xsn3n(niso), pnum(niso)
REAL :: phi(ng)
REAL :: temp
INTEGER :: idiso(niso)

TYPE(LIBDATA), POINTER :: isodata
REAL :: IsoXsMacA(niso, ng), IsoXsMacS0(niso, ng), IsoXsMacTr(niso, ng), IsoXsMacf(niso, ng)
INTEGER :: ig, iso, id, idn2n, idn3n, it1, it2
REAL :: wt1, wt2, phisum                       !Temperature interpolation

REAL :: phitil;



!IF(.NOT. lAlloc) CALL AllocIsotopicMacXs()
xsa = 0; xsf = 0;

DO iso = 1, niso
  id = mapnucl(idiso(iso));   isodata => ldiso(id)
  idn2n=mapn2n(id);idn3n=mapn3n(id)
  IF (idn2n.NE.0) THEN
  DO ig = 1, ng
      IsoXsMacS0(iso, ig) = pnum(iso)*isodata%sign2n(ig)
    END DO
  ELSE
    IsoXsMacS0 = 0.;
  END IF
  IF (idn3n.NE.0) THEN
    DO ig = 1, ng
      IsoXsMacTr(iso, ig) = pnum(iso)*isodata%sign3n(ig)
  ENDDO
  ELSE
    IsoXsMacTr = 0.;
  END IF
ENDDO !Isotopic sweeep

!Group Condense
xsn2n(1:niso) = 0.
xsn3n(1:niso) = 0.


!OPT -- Edit by LHG : 2020.07.13
phisum = 0
DO ig= 1, ng
  phisum = phisum + phi(ig)
ENDDO
DO ig= 1, ng
  phitil = phi(ig)/phisum
DO iso = 1, niso
    xsn2n(iso) = xsn2n(iso) + IsoXsMacs0(iso, ig) * phitil
    xsn3n(iso) = xsn3n(iso) + IsoXsMactr(iso, ig) * phitil
ENDDO

ENDDO
END SUBROUTINE

SUBROUTINE MacXsNf_CrCsp(XsMac, Fxr, ig1, ig2, ng, eigv, kin, lCrCspFtn)
USE XsUtil_mod,    ONLY : GetXsMacDat, ReturnXsMacDat
USE CrCsp_Mod,     ONLY : MacXsNfCsp
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
TYPE(Fxrinfo_type) :: Fxr
REAL :: eigv
INTEGER :: ig1, ig2, ng
LOGICAL :: kin, lCrCspFtn

TYPE(XsMac_TYpe), POINTER :: XsMac1, XsMac2
REAL :: temp
INTEGER :: niso
INTEGER, POINTER :: idiso(:)
REAL, POINTER ::pnum(:)

IF(lCrCspFtn .AND. Fxr%lCrCspFtn) THEN
  CALL GetXsMacDat(XsMac1, ng, .FALSE.)
  Call GetXsMacDat(XsMac2, ng, .FALSE.)
  niso = Fxr%CspFxr%niso(1); Temp = Fxr%Temp
  Pnum => Fxr%CspFxr%pnum(:, 1);   idiso => Fxr%CspFxr%isolist(:,1)
  CALL MacXsNf_gen(XsMac1, niso, temp, idiso, pnum, ig1, ig2, ng, eigv, kin)

  niso = Fxr%CspFxr%niso(2); Temp = Fxr%Temp
  Pnum => Fxr%CspFxr%pnum(:, 2);   idiso => Fxr%CspFxr%isolist(:,2)
  CALL MacXsNf_gen(XsMac2, niso, temp, idiso, pnum, ig1, ig2, ng, eigv, kin)

  CALL MacXsNfCsp(XsMac, XsMac1, XsMac2, Fxr%CspFxr, ig1, ig2)
  NULLIFY(pnum); NULLIFY(idiso)
  CALL ReturnXsMacDat(XsMac1); CALL ReturnXsMacDat(XsMac2)
ELSE
  niso = Fxr%niso; temp = Fxr%Temp
  pnum => Fxr%pnum; idiso => Fxr%idiso
  CALL MacXsNf_gen(XsMac, niso, temp, idiso, pnum, ig1, ig2, ng, eigv, kin)
  NULLIFY(pnum); NULLIFY(idiso)
ENDIF

END SUBROUTINE

SUBROUTINE MacXsNf_Fxr(XsMac, Fxr, ig1, ig2, ng, eigv, kin)
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
TYPE(Fxrinfo_type) :: Fxr
REAL :: eigv
INTEGER :: ig1, ig2, ng
LOGICAL :: kin
REAL :: temp
INTEGER :: niso
INTEGER, POINTER :: idiso(:)
REAL, POINTER ::pnum(:)

niso = Fxr%niso; temp = Fxr%Temp
pnum => Fxr%pnum; idiso => Fxr%idiso
CALL MacXsNf_gen(XsMac, niso, temp, idiso, pnum, ig1, ig2, ng, eigv, kin)
NULLIFY(pnum); NULLIFY(idiso)

END SUBROUTINE

SUBROUTINE MacXsNf_Gen(XsMac, niso, temp, idiso, pnum, ig1, ig2, ng, eigv, kin)
USE XsUtil_mod,  ONLY : AllocXsMac,          AllocMacIsoXs,  FreeXsIsoMac
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
REAL :: temp, eigv
INTEGER :: niso
INTEGER :: idiso(niso)
REAL :: pnum(niso)
INTEGER :: ig1, ig2, ng
LOGICAL :: kin

REAL, POINTER :: xsmacnf(:), ISOMACNF(:,:)
TYPE(LIBDATA), POINTER :: isodata
INTEGER :: iso, id, ig
REAL :: reigv

INTEGER :: it1, it2
REAL :: wt1, wt2

IF(.NOT. XsMac%lalloc) THEN
  XsMac%ng = ng
  CALL AllocXsMac(XsMac)
ENDIF
!IF (.NOT. XsMac%lIsoAlloc .OR. XsMac%niso .NE. niso) THEN  ! for debugging... JSU 2019/08/12
!  CALL FreeXsIsoMac(XsMac)
!  CALL AllocMacIsoXs(XsMac,NG,NISO)
!END IF

XsMacNf => XsMac%XsMacNf;! ISOMACNF=> XsMac%ISOXsMacNf

XsMacNf(ig1:ig2) = 0.

!ISOMACNF = 0._8
reigv = 1._8 / eigv

DO iso = 1, niso
  id = MapNucl(idiso(iso)); isodata => ldiso(id)
  IF(mapfis(idiso(iso)) .EQ. 0) CYCLE
  CALL XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2)
  DO ig = ig1, ig2
    XsMacNf(ig) = XsMacNf(ig) + pnum(iso) * (wt1 * isodata%signf(ig, it1) + wt2 * isodata%signf(ig, it2))
!    ISOMACNF(iso, ig) = pnum(iso) * (wt1 * isodata%signf(ig, it1) + wt2 * isodata%signf(ig, it2))   ! for debugging... JSU 2019/08/12
  ENDDO
ENDDO

IF(kin) XsMacNf(ig1:ig2) = XsMacNf(ig1:ig2)*reigv

NULLIFY(XsMacNf, isodata)
END SUBROUTINE

!BYSedit
SUBROUTINE IsoMacXsNf(XsMac, niso, temp, idiso, pnum, ig1, ig2, ng, eigv, kin)
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
REAL :: temp, eigv
INTEGER :: niso
INTEGER :: idiso(niso)
REAL :: pnum(niso)
INTEGER :: ig1, ig2, ng
LOGICAL :: kin

REAL, POINTER :: xsmacnf(:), IsoXsMacKf(:, :), IsoXsMacnf(:, :)
TYPE(LIBDATA), POINTER :: isodata
INTEGER :: iso, id, ig
REAL :: reigv

INTEGER :: it1, it2
REAL :: wt1, wt2

IF(.NOT. XsMac%lalloc) THEN
  XsMac%ng = ng
  CALL AllocXsMac(XsMac)
ENDIF
IF(.NOT. XsMac%lIsoAlloc) THEN
  CALL AllocMacIsoXs(XsMac, ng, nelthel)
ENDIF

XsMacNf => XsMac%XsMacNf
IsoXsMacKf => XsMac%IsoXsMacKf
IsoXsMacnf => XsMac%IsoXsMacnf

CALL CP_CA(XsMacNf(ig1:ig2), 0._8, ig2 - ig1 + 1)

reigv = 1._8 / eigv

DO iso = 1, niso
  id = MapNucl(idiso(iso)); isodata => ldiso(id)
  IF(mapfis(idiso(iso)) .EQ. 0) CYCLE
  CALL XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2)
  DO ig = ig1, ig2
    XsMacNf(ig) = XsMacNf(ig) + pnum(iso) * (wt1 * isodata%signf(ig, it1) + wt2 * isodata%signf(ig, it2))
    IsoXsMacnf(iso, ig) = pnum(iso) * (wt1 * isodata%signf(ig, it1) + wt2 * isodata%signf(ig, it2))
    IsoXsMacKf(iso, ig) = pnum(iso) * (wt1 * isodata%sigf(ig, it1) + wt2 * isodata%sigf(ig, it2)) *isodata%kappa
    !--- BYS edit missing isodata%kappa term added/ signf-> sigf
  ENDDO
ENDDO

NULLIFY(XsMacNf, isodata)
END SUBROUTINE
!--- BYS edit end
SUBROUTINE MacXsScatMatrix_Csp(XsMac, FXR, igb, ige, ng, GroupInfo, lscat1, lCrCspFtn)
USE XsUtil_mod,        ONLY : GetXsMacDat,       ReturnXsMacDat
USE CrCsp_Mod,         ONLY : MacXsSmCsp
USE CNTL,              ONLY : nTracerCntl
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
TYPE(FxrInfo_Type) :: Fxr
INTEGER, INTENT(IN) :: igb, ige, ng
TYPE(GroupInfo_Type), INTENT(IN) :: GroupInfo
LOGICAL :: lscat1, lCrCspFtn

TYPE(XsMac_Type), POINTER :: XsMac1, XsMac2
REAL :: temp
INTEGER :: niso
INTEGER, POINTER :: idiso(:)
REAL, POINTER ::pnum(:)
REAL(4), POINTER :: fresoSIso(:,:),fresoSSIso(:,:),fresoS1Iso(:,:)

IF(lCrCspFtn .AND. Fxr%lCrCspFtn) THEN
  CALL GetXsMacDat(XsMac1, ng, .FALSE.)
  Call GetXsMacDat(XsMac2, ng, .FALSE.)
  niso = Fxr%CspFxr%niso(1); Temp = Fxr%Temp
  Pnum => Fxr%CspFxr%pnum(:, 1);   idiso => Fxr%CspFxr%isolist(:,1)
  if (nTracerCntl%lRst) then
    fresoSIso => Fxr%fresoSIso; fresoSSIso => Fxr%fresoSSIso; fresoS1Iso => Fxr%fresoS1Iso
    CALL MacXsScatMatrix_GenR(XsMac1, Temp, niso, idiso, pnum, igb, ige, ng, GroupInfo, lscat1, Fxr%lres, fresoSIso, fresoSSIso, fresoS1Iso)
  else
    CALL MacXsScatMatrix_Gen(XsMac1, Temp, niso, idiso, pnum, igb, ige, ng, GroupInfo, lscat1)
  endif
  niso = Fxr%CspFxr%niso(2); Temp = Fxr%Temp
  Pnum => Fxr%CspFxr%pnum(:, 2);   idiso => Fxr%CspFxr%isolist(:,2)
  if (nTracerCntl%lRst) then
    CALL MacXsScatMatrix_GenR(XsMac2, Temp, niso, idiso, pnum, igb, ige, ng, GroupInfo, lscat1, Fxr%lres, fresoSIso, fresoSSIso, fresoS1Iso)
  else
    CALL MacXsScatMatrix_Gen(XsMac1, Temp, niso, idiso, pnum, igb, ige, ng, GroupInfo, lscat1)
  endif
  CALL MacXsSmCsp(XsMac, XsMac1, XsMac2, Fxr%CspFxr, igb, ige, ng)
  !CALL MacXsNfCsp(XsMac, XsMac1, XsMac2, Fxr%CspFxr, ig1, ig2)
  NULLIFY(pnum); NULLIFY(idiso)
  CALL ReturnXsMacDat(XsMac1); CALL ReturnXsMacDat(XsMac2)
ELSE
  niso = Fxr%niso; temp = Fxr%Temp
  pnum => Fxr%pnum; idiso => Fxr%idiso
  if (nTracerCntl%lRst) then
    fresoSIso => Fxr%fresoSIso; fresoSSIso => Fxr%fresoSSIso; fresoS1Iso => Fxr%fresoS1Iso
    CALL MacXsScatMatrix_GenR(XsMac, Temp, niso, idiso, pnum, igb, ige, ng, GroupInfo, lscat1, Fxr%lres, fresoSIso, fresoSSIso, fresoS1Iso)
  else
    CALL MacXsScatMatrix_Gen(XsMac, Temp, niso, idiso, pnum, igb, ige, ng, GroupInfo, lscat1)
  endif
  NULLIFY(pnum); NULLIFY(idiso)
ENDIF
END SUBROUTINE

SUBROUTINE MacXsScatMatrix_Fxr(XsMac, FXR, igb, ige, ng, GroupInfo, lscat1)
USE CNTL,              ONLY : nTracerCntl
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
TYPE(FxrInfo_Type) :: Fxr
INTEGER, INTENT(IN) :: igb, ige, ng
TYPE(GroupInfo_Type), INTENT(IN) :: GroupInfo
LOGICAL :: lscat1

REAL :: temp
INTEGER :: niso
INTEGER, POINTER :: idiso(:)
REAL, POINTER ::pnum(:)
REAL(4), POINTER :: fresoSIso(:,:),fresoSSIso(:,:),fresoS1Iso(:,:)

niso = Fxr%niso; temp = Fxr%Temp
pnum => Fxr%pnum; idiso => Fxr%idiso
if (nTracerCntl%lRst) then
  fresoSIso => Fxr%fresoSIso; fresoSSIso => Fxr%fresoSSIso; fresoS1Iso => Fxr%fresoS1Iso
  CALL MacXsScatMatrix_GenR(XsMac, Temp, niso, idiso, pnum, igb, ige, ng, GroupInfo, lscat1, Fxr%lres, fresoSIso, fresoSSIso, fresoS1Iso)
else
  CALL MacXsScatMatrix_Gen(XsMac, Temp, niso, idiso, pnum, igb, ige, ng, GroupInfo, lscat1)
endif
NULLIFY(pnum); NULLIFY(idiso)
END SUBROUTINE

SUBROUTINE MacXsScatMatrix_GenR(XsMac, temp, niso, idiso, pnum, igb, ige, ng, GroupInfo, lscat1, fxrlres, fresoSIso, fresoSSIso, fresoS1Iso)
! Resonance Scattering IS CONSIDERED.
! XsMacSm, XsMacS, XsMacStr are updated.
USE CNTL, ONLY : nTracerCntl
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
INTEGER, INTENT(IN) :: niso, igb, ige, ng
INTEGER, INTENT(IN) :: idiso(niso)
REAL, INTENT(IN) :: pnum(niso)
REAL, INTENT(IN) :: temp
REAL(4), POINTER :: fresoSIso(:,:), fresoSSIso(:,:), fresoS1Iso(:,:)
TYPE(GroupInfo_Type), INTENT(IN) :: GroupInfo

LOGICAL :: lscat1, fxrlres, lres

TYPE(LIBDATA), POINTER :: isodata
INTEGER :: i, j, k
INTEGER :: iso, id, p1idx, it1, it2, it1p1, it2p1, ig, ig2, ig3, ib, ie, ind
INTEGER :: ScRange1(2), ScRange2(2) !Scattering Range
REAL :: wt1, wt2, wt1p1, wt2p1, a1, a2, ND
REAL :: sigs, sigss, sigsp1, sigsout, sigs1, sigs2, sigs_lib, sigss_lib, sigsout_lib, rat
REAL, POINTER :: Xsmacsm(:,:), XsMacS(:), XsMacStr(:)

IF(.NOT. XsMac%lalloc) THEN
  XsMac%ng = ng
  CALL AllocXsMac(XsMac)
ENDIF
!IF(.NOT. lAlloc) CALL AllocIsotopicMacXs()

XsMacSm => XsMac%XsMacSm
XsMacS => XsMac%XsMacS
XsMacStr => XsMac%XsMacStr

XsMacSm(1:ng,1:ng) = 0.
XsMacS(1:ng) = 0.
XsMacStr(1:ng) = 0.

IF ( fxrlres .AND. nTracerCntl%lSubGrpSweep ) THEN  !@@@@@ Resonance Treatment

  DO iso = 1, niso !%%%%% Isotope Sweep

    id = mapnucl(idiso(iso));   isodata => ldiso(id)
    p1idx = mapnuclp13(idiso(iso))
    ND = pnum(iso)
    CALL XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2)
    if (p1idx.ne.0) CALL p1XsTempInterpolation(p1idx, isodata, temp, wt1p1, wt2p1, it1p1, it2p1)
    DO ig = igb, ige  !!! Group Sweep

      IF (isodata%lreso) then !!*****!! Resonance isotope !!*****!!

        !Determine Scattering Range
        ScRange1 = (/isodata%sm(ig, it1)%ib,isodata%sm(ig, it1)%ie/)
        ScRange2 = (/isodata%sm(ig, it2)%ib,isodata%sm(ig, it2)%ie/)
        ib = min(ScRange1(1), ScRange2(1)); ie = max(ScRange1(2), ScRange2(2))

        DO ig2 = ib, ie  !!!&&!!! In Scattering Group Sweep

          IF ( (ig2 .GE. igresb) .AND. (ig2 .LE. igrese) ) THEN   !!! If ig2 is a resonance group

            IF ( ig2 .NE. ig ) THEN ! Non Self-scattering
              sigs_lib = wt1 * isodata%sigs(ig2, it1) + wt2 * isodata%sigs(ig2, it2)
              sigss_lib = wt1 * isodata%sigss(ig2, it1) + wt2 * isodata%sigss(ig2, it2)
              sigsout_lib = sigs_lib - sigss_lib

              sigs1 = 0; sigs2 = 0
              ind = (ig2 - ScRange1(1)) * (ig2 - ScRange1(2))
              IF(ind .LE. 0) sigs1 = isodata%sm(ig, it1)%from(ig2)
              ind = (ig2 - ScRange2(1)) * (ig2 - ScRange2(2))
              IF(ind .LE. 0) sigs2 = isodata%sm(ig, it2)%from(ig2)
              rat = (wt1 * sigs1 + wt2 * sigs2) / sigsout_lib

              sigs = sigs_lib * fresoSIso(iso,ig2) ! Resonance treated sigs
              sigss = sigss_lib * fresoSSIso(iso,ig2) ! Resonance treated sigss
              sigsout = sigs - sigss

              XsMacSm(ig2, ig) = XsMacSm(ig2, ig) + ND * rat * sigsout
            ELSE ! Self-scattering
              sigss = (wt1 * isodata%sigss(ig2, it1) + wt2 * isodata%sigss(ig2, it2)) * fresoSSIso(iso,ig2)
              if (p1idx.ne.0) then
                sigsp1 = (wt1p1 * isodata%sigsp1(ig2, it1p1) + wt2p1 * isodata%sigsp1(ig2, it2p1)) * fresoS1Iso(iso,ig2)
              else
                sigsp1 = (wt1 * (isodata%sigs(ig2, it1)-isodata%sigstr(ig2, it1)) + wt2 * (isodata%sigs(ig2, it2)-isodata%sigstr(ig2, it2))) * fresoSIso(iso,ig2)
              endif
              XsMacSm(ig2, ig) = XsMacSm(ig2, ig) + ND * (sigss - sigsp1)
            ENDIF

          ELSE   !!! If ig2 is not a resonance group

            sigs1 = 0; sigs2 = 0
            ind = (ig2 - ScRange1(1)) * (ig2 - ScRange1(2))
            IF(ind .LE. 0) sigs1 = isodata%sm(ig, it1)%from(ig2)
            ind = (ig2 - ScRange2(1)) * (ig2 - ScRange2(2))
            IF(ind .LE. 0) sigs2 = isodata%sm(ig, it2)%from(ig2)
            XsMacSm(ig2, ig) = XsMacSm(ig2, ig) + ND * (wt1 * sigs1 + wt2 * sigs2)

          ENDIF

        ENDDO     !!!&&!!! END In Scattering Group Sweep

        IF ( (ig .GE. igresb) .AND. (ig .LE. igrese) ) THEN
            sigs = ND * (wt1 * isodata%sigs(ig, it1) + wt2 * isodata%sigs(ig, it2)) *  fresoSIso(iso,ig)
            if (p1idx.ne.0) then
              sigsp1 = ND * (wt1p1 * isodata%sigsp1(ig, it1p1) + wt2p1 * isodata%sigsp1(ig, it2p1)) * fresoS1Iso(iso,ig)
            else
              sigsp1 = ND * (wt1 * (isodata%sigs(ig, it1)-isodata%sigstr(ig, it1)) + wt2 * (isodata%sigs(ig, it2)-isodata%sigstr(ig, it2))) * fresoSIso(iso,ig)
            endif
            XsMacS(ig) =  XsMacS(ig) + sigs
            XsMacStr(ig) =  XsMacStr(ig) + (sigs - sigsp1)
        ELSE
            XsMacS(ig) =  XsMacS(ig) + ND * (wt2 * isodata%sigs(ig, it2) + wt1 * isodata%sigs(ig, it1))
            XsMacStr(ig) =  XsMacStr(ig) + ND * (wt2 * isodata%sigstr(ig, it2) + wt1 * isodata%sigstr(ig, it1))
        ENDIF


      ELSE     !!*****!! Non Resonance isotope !!*****!!

        !Determine Scattering Range
        ScRange1 = (/isodata%sm(ig, it1)%ib,isodata%sm(ig, it1)%ie/)
        ScRange2 = (/isodata%sm(ig, it2)%ib,isodata%sm(ig, it2)%ie/)
        ib = min(ScRange1(1), ScRange2(1)); ie = max(ScRange1(2), ScRange2(2))
        DO ig2 = ib, ie  !In Scattering Group Sweep
          sigs1 = 0; sigs2 = 0
          ind = (ig2 - ScRange1(1)) * (ig2 - ScRange1(2))
          IF(ind .LE. 0) sigs1 = isodata%sm(ig, it1)%from(ig2)
          ind = (ig2 - ScRange2(1)) * (ig2 - ScRange2(2))
          IF(ind .LE. 0) sigs2 = isodata%sm(ig, it2)%from(ig2)
          XsMacSm(ig2, ig) = XsMacSm(ig2, ig) + ND * (wt1 * sigs1 + wt2 * sigs2)
        ENDDO
        XsMacS(ig) =  XsMacS(ig) + ND * (wt2 * isodata%sigs(ig, it2) + wt1 * isodata%sigs(ig, it1))
        XsMacStr(ig) =  XsMacStr(ig) + ND * (wt2 * isodata%sigstr(ig, it2) + wt1 * isodata%sigstr(ig, it1))

      ENDIF

    ENDDO !!! END Group Sweep

    NULLIFY(isodata)

  ENDDO   !%%%%% END Isotope Sweep

ELSE !@@@@@ Non Resonance Treatment

  DO iso = 1, niso
    id = mapnucl(idiso(iso));   isodata => ldiso(id)
    ND = pnum(iso)
    CALL XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2)
    DO ig = igb, ige
      !Determine Scattering Range
      ScRange1 = (/isodata%sm(ig, it1)%ib,isodata%sm(ig, it1)%ie/)
      ScRange2 = (/isodata%sm(ig, it2)%ib,isodata%sm(ig, it2)%ie/)
      ib = min(ScRange1(1), ScRange2(1)); ie = max(ScRange1(2), ScRange2(2))
      DO ig2 = ib, ie  !In Scattering Group Sweep
        sigs1 = 0; sigs2 = 0
        ind = (ig2 - ScRange1(1)) * (ig2 - ScRange1(2))
        IF(ind .LE. 0) sigs1 = isodata%sm(ig, it1)%from(ig2)
        ind = (ig2 - ScRange2(1)) * (ig2 - ScRange2(2))
        IF(ind .LE. 0) sigs2 = isodata%sm(ig, it2)%from(ig2)
        XsMacSm(ig2, ig) = XsMacSm(ig2, ig) + ND * (wt1 * sigs1 + wt2 * sigs2)
      ENDDO
      XsMacS(ig) =  XsMacS(ig) + ND * (wt2 * isodata%sigs(ig, it2) + wt1 * isodata%sigs(ig, it1))
      XsMacStr(ig) =  XsMacStr(ig) + ND * (wt2 * isodata%sigstr(ig, it2) + wt1 * isodata%sigstr(ig, it1))
    ENDDO
    NULLIFY(isodata)
  ENDDO

ENDIF

DO ig = igb, ige
  IF(lScat1) XsMacSm(ig, ig) = XsMacSm(ig, ig) + (XsMacS(ig) - XsMacStr(ig))
ENDDO

NULLIFY(XsMacSm, XsMacs, XsMacStr)
END SUBROUTINE
SUBROUTINE MacXsScatMatrix_Gen(XsMac, temp, niso, idiso, pnum, igb, ige, ng, GroupInfo, lscat1)
! Resonance Scattering IS NOT CONSIDERED.
! XsMacSm, XsMacS, XsMacStr are updated.
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
INTEGER, INTENT(IN) :: niso, igb, ige, ng
INTEGER, INTENT(IN) :: idiso(niso)
REAL, INTENT(IN) :: pnum(niso)
REAL, INTENT(IN) :: temp

TYPE(GroupInfo_Type), INTENT(IN) :: GroupInfo

LOGICAL :: lscat1

TYPE(LIBDATA), POINTER :: isodata
INTEGER :: iso, id, it1, it2, ig, ig2, ib, ie, ind
INTEGER :: ScRange1(2), ScRange2(2) !Scattering Range
REAL :: wt1, wt2, ND, ND1, ND2
REAL :: sigs1, sigs2
REAL, POINTER :: Xsmacsm(:,:), XsMacS(:), XsMacStr(:)

IF(.NOT. XsMac%lalloc) THEN
  XsMac%ng = ng
  CALL AllocXsMac(XsMac)
ENDIF
!IF(.NOT. lAlloc) CALL AllocIsotopicMacXs()

XsMacSm => XsMac%XsMacSm
XsMacS => XsMac%XsMacS
XsMacStr => XsMac%XsMacStr

XsMacSm(1:ng,1:ng) = 0.
XsMacS(1:ng) = 0.
XsMacStr(1:ng) = 0.

DO iso = 1, niso
  id = mapnucl(idiso(iso));   isodata => ldiso(id)
  ND = pnum(iso)
  CALL XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2)
  ND1 = ND*wt1; ND2 = ND*wt2
  DO ig = igb, ige
    !Determine Scattering Range
    ScRange1 = (/isodata%sm(ig, it1)%ib,isodata%sm(ig, it1)%ie/)
    ScRange2 = (/isodata%sm(ig, it2)%ib,isodata%sm(ig, it2)%ie/)
!    ib = min(ScRange1(1), ScRange2(1)); ie = max(ScRange1(2), ScRange2(2))
!    DO ig2 = ib, ie  !In Scattering Group Sweep
!      sigs1 = 0; sigs2 = 0
!      ind = (ig2 - ScRange1(1)) * (ig2 - ScRange1(2))
!      IF(ind .LE. 0) sigs1 = isodata%sm(ig, it1)%from(ig2)
!      ind = (ig2 - ScRange2(1)) * (ig2 - ScRange2(2))
!      IF(ind .LE. 0) sigs2 = isodata%sm(ig, it2)%from(ig2)
!      XsMacSm(ig2, ig) = XsMacSm(ig2, ig) + ND * (wt1 * sigs1 + wt2 * sigs2)
!    ENDDO
    DO ig2 = ScRange1(1),ScRange1(2)
      XsMacSm(ig2,ig) = XsMacSm(ig2,ig)+ND1*isodata%sm(ig,it1)%from(ig2)
    END DO
    DO ig2 = ScRange2(1),ScRange2(2)
      XsMacSm(ig2,ig) = XsMacSm(ig2,ig)+ND2*isodata%sm(ig,it2)%from(ig2)
    ENDDO
    !XsMacS(ig) =  XsMacS(ig) + ND * (wt2 * isodata%sigs(ig, it2) + wt1 * isodata%sigs(ig, it1))
    !XsMacStr(ig) =  XsMacStr(ig) + ND * (wt2 * isodata%sigstr(ig, it2) + wt1 * isodata%sigstr(ig, it1))
    XsMacS(ig) =  XsMacS(ig) + ND2 * isodata%sigs(ig, it2) + ND1 * isodata%sigs(ig, it1)
    XsMacStr(ig) =  XsMacStr(ig) + ND2 * isodata%sigstr(ig, it2) + ND1 * isodata%sigstr(ig, it1)
  ENDDO
  NULLIFY(isodata)
ENDDO


DO ig = igb, ige
  IF(lScat1) XsMacSm(ig, ig) = XsMacSm(ig, ig) + (XsMacS(ig) - XsMacStr(ig))
ENDDO

NULLIFY(XsMacSm, XsMacs, XsMacStr)
END SUBROUTINE

SUBROUTINE MacXsScatMatrix_Cusping(XsMac, iso, Fxr, igb, ige, ng, GroupInfo, lscat1)
USE Material_Mod,  ONLY: Mixture
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
TYPE(FxrInfo_Type) :: Fxr
INTEGER, INTENT(IN) :: igb, ige, ng
TYPE(GroupInfo_Type), INTENT(IN) :: GroupInfo
LOGICAL :: lscat1
INTEGER :: iso
REAL :: temp
INTEGER :: niso
INTEGER, POINTER :: idiso(:)
REAL, POINTER ::pnum(:)

niso = Mixture(iso)%niso; temp = Fxr%Temp
pnum => Mixture(iso)%pnum; idiso => Mixture(iso)%idiso
CALL MacXsScatMatrix_Gen(XsMac, Temp, niso, idiso, pnum, igb, ige, ng, GroupInfo, lscat1)
NULLIFY(pnum); NULLIFY(idiso)
END SUBROUTINE


SUBROUTINE IsoMacXsScatMatrix_GenR(IsoXsMacSm, temp, niso, idiso, pnum, igb, ige, ng, GroupInfo, lscat1, fxrlres, fresoSIso, fresoSSIso, fresoS1Iso)
! Resonance Scattering IS CONSIDERED.
USE CNTL, ONLY : nTracerCntl
IMPLICIT NONE
INTEGER, INTENT(IN) :: niso, igb, ige, ng
INTEGER, INTENT(IN) :: idiso(niso)
REAL, INTENT(IN) :: pnum(niso)
REAL, INTENT(IN) :: temp
REAL(4), POINTER :: fresoSIso(:,:), fresoSSIso(:,:), fresoS1Iso(:,:)
TYPE(GroupInfo_Type), INTENT(IN) :: GroupInfo
LOGICAL :: lscat1, fxrlres, lres

TYPE(LIBDATA), POINTER :: isodata
INTEGER :: i, j, k
INTEGER :: iso, id, p1idx, it1, it2, it1p1, it2p1, ig, ig2, ig3, ib, ie, ind
INTEGER :: ScRange1(2), ScRange2(2) !Scattering Range
REAL :: wt1, wt2, wt1p1, wt2p1, a1, a2, ND
REAL :: sigs, sigss, sigsp1, sigsout, sigs1, sigs2, sigs_lib, sigss_lib, sigsout_lib, rat
REAL :: IsoXsmacsm(niso,ng,ng),IsoXsMacS(niso,ng),IsoXsMacStr(niso,ng)

IsoXsMacSm(1:niso,1:ng,1:ng) = 0._8
IsoXsMacS(1:niso,1:ng) = 0._8
IsoXsMacStr(1:niso,1:ng) = 0._8
IF ( fxrlres .AND. nTracerCntl%lSubGrpSweep ) THEN  !@@@@@ Resonance Treatment

  DO iso = 1, niso !%%%%% Isotope Sweep

    id = mapnucl(idiso(iso));   isodata => ldiso(id)
    p1idx = mapnuclp13(idiso(iso))
    ND = pnum(iso)
    CALL XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2)
    if (p1idx.ne.0) CALL p1XsTempInterpolation(p1idx, isodata, temp, wt1p1, wt2p1, it1p1, it2p1)
    DO ig = igb, ige  !!! Group Sweep

      IF (isodata%lreso) then !!*****!! Resonance isotope !!*****!!

        !Determine Scattering Range
        ScRange1 = (/isodata%sm(ig, it1)%ib,isodata%sm(ig, it1)%ie/)
        ScRange2 = (/isodata%sm(ig, it2)%ib,isodata%sm(ig, it2)%ie/)
        ib = min(ScRange1(1), ScRange2(1)); ie = max(ScRange1(2), ScRange2(2))

        DO ig2 = ib, ie  !!!&&!!! In Scattering Group Sweep

          IF ( (ig2 .GE. igresb) .AND. (ig2 .LE. igrese) ) THEN   !!! If ig2 is a resonance group

            IF ( ig2 .NE. ig ) THEN ! Non Self-scattering
              sigs_lib = wt1 * isodata%sigs(ig2, it1) + wt2 * isodata%sigs(ig2, it2)
              sigss_lib = wt1 * isodata%sigss(ig2, it1) + wt2 * isodata%sigss(ig2, it2)
              sigsout_lib = sigs_lib - sigss_lib

              sigs1 = 0; sigs2 = 0
              ind = (ig2 - ScRange1(1)) * (ig2 - ScRange1(2))
              IF(ind .LE. 0) sigs1 = isodata%sm(ig, it1)%from(ig2)
              ind = (ig2 - ScRange2(1)) * (ig2 - ScRange2(2))
              IF(ind .LE. 0) sigs2 = isodata%sm(ig, it2)%from(ig2)
              rat = (wt1 * sigs1 + wt2 * sigs2) / sigsout_lib

              sigs = sigs_lib * fresoSIso(iso,ig2) ! Resonance treated sigs
              sigss = sigss_lib * fresoSSIso(iso,ig2) ! Resonance treated sigss
              sigsout = sigs - sigss

              IsoXsMacSm(iso, ig2, ig) = IsoXsMacSm(iso, ig2, ig) + ND * rat * sigsout
            ELSE ! Self-scattering
              sigss = (wt1 * isodata%sigss(ig2, it1) + wt2 * isodata%sigss(ig2, it2)) * fresoSSIso(iso,ig2)
              if (p1idx.ne.0) then
                sigsp1 = (wt1p1 * isodata%sigsp1(ig2, it1p1) + wt2p1 * isodata%sigsp1(ig2, it2p1)) * fresoS1Iso(iso,ig2)
              else
                sigsp1 = (wt1 * (isodata%sigs(ig2, it1)-isodata%sigstr(ig2, it1)) + wt2 * (isodata%sigs(ig2, it2)-isodata%sigstr(ig2, it2))) * fresoSIso(iso,ig2)
              endif
              IsoXsMacSm(iso, ig2, ig) = IsoXsMacSm(iso, ig2, ig) + ND * (sigss - sigsp1)
            ENDIF

          ELSE   !!! If ig2 is not a resonance group

            sigs1 = 0; sigs2 = 0
            ind = (ig2 - ScRange1(1)) * (ig2 - ScRange1(2))
            IF(ind .LE. 0) sigs1 = isodata%sm(ig, it1)%from(ig2)
            ind = (ig2 - ScRange2(1)) * (ig2 - ScRange2(2))
            IF(ind .LE. 0) sigs2 = isodata%sm(ig, it2)%from(ig2)
            IsoXsMacSm(iso, ig2, ig) = IsoXsMacSm(iso, ig2, ig) + ND * (wt1 * sigs1 + wt2 * sigs2)

          ENDIF

        ENDDO     !!!&&!!! END In Scattering Group Sweep

        IF ( (ig .GE. igresb) .AND. (ig .LE. igrese) ) THEN
            sigs = ND * (wt1 * isodata%sigs(ig, it1) + wt2 * isodata%sigs(ig, it2)) *  fresoSIso(iso,ig)
            if (p1idx.ne.0) then
              sigsp1 = ND * (wt1p1 * isodata%sigsp1(ig, it1p1) + wt2p1 * isodata%sigsp1(ig, it2p1)) * fresoS1Iso(iso,ig)
            else
              sigsp1 = ND * (wt1 * (isodata%sigs(ig, it1)-isodata%sigstr(ig, it1)) + wt2 * (isodata%sigs(ig, it2)-isodata%sigstr(ig, it2))) * fresoSIso(iso,ig)
            endif
            IsoXsMacS(iso, ig) =  IsoXsMacS(iso, ig) + sigs
            IsoXsMacStr(iso, ig) =  IsoXsMacStr(iso, ig) + (sigs - sigsp1)
        ELSE
            IsoXsMacS(iso, ig) =  IsoXsMacS(iso, ig) + ND * (wt2 * isodata%sigs(ig, it2) + wt1 * isodata%sigs(ig, it1))
            IsoXsMacStr(iso, ig) =  IsoXsMacStr(iso, ig) + ND * (wt2 * isodata%sigstr(ig, it2) + wt1 * isodata%sigstr(ig, it1))
        ENDIF


      ELSE     !!*****!! Non Resonance isotope !!*****!!

        !Determine Scattering Range
        ScRange1 = (/isodata%sm(ig, it1)%ib,isodata%sm(ig, it1)%ie/)
        ScRange2 = (/isodata%sm(ig, it2)%ib,isodata%sm(ig, it2)%ie/)
        ib = min(ScRange1(1), ScRange2(1)); ie = max(ScRange1(2), ScRange2(2))
        DO ig2 = ib, ie  !In Scattering Group Sweep
          sigs1 = 0; sigs2 = 0
          ind = (ig2 - ScRange1(1)) * (ig2 - ScRange1(2))
          IF(ind .LE. 0) sigs1 = isodata%sm(ig, it1)%from(ig2)
          ind = (ig2 - ScRange2(1)) * (ig2 - ScRange2(2))
          IF(ind .LE. 0) sigs2 = isodata%sm(ig, it2)%from(ig2)
          IsoXsMacSm(iso, ig2, ig) = IsoXsMacSm(iso, ig2, ig) + ND * (wt1 * sigs1 + wt2 * sigs2)
        ENDDO
        IsoXsMacS(iso, ig) =  IsoXsMacS(iso, ig) + ND * (wt2 * isodata%sigs(ig, it2) + wt1 * isodata%sigs(ig, it1))
        IsoXsMacStr(iso, ig) =  IsoXsMacStr(iso, ig) + ND * (wt2 * isodata%sigstr(ig, it2) + wt1 * isodata%sigstr(ig, it1))

      ENDIF

    ENDDO !!! END Group Sweep

    NULLIFY(isodata)

  ENDDO   !%%%%% END Isotope Sweep

ELSE !@@@@@ Non Resonance Treatment

  DO iso = 1, niso
    id = mapnucl(idiso(iso));   isodata => ldiso(id)
    ND = pnum(iso)
    CALL XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2)
    DO ig = igb, ige
      !Determine Scattering Range
      ScRange1 = (/isodata%sm(ig, it1)%ib,isodata%sm(ig, it1)%ie/)
      ScRange2 = (/isodata%sm(ig, it2)%ib,isodata%sm(ig, it2)%ie/)
      ib = min(ScRange1(1), ScRange2(1)); ie = max(ScRange1(2), ScRange2(2))
      DO ig2 = ib, ie  !In Scattering Group Sweep
        sigs1 = 0; sigs2 = 0
        ind = (ig2 - ScRange1(1)) * (ig2 - ScRange1(2))
        IF(ind .LE. 0) sigs1 = isodata%sm(ig, it1)%from(ig2)
        ind = (ig2 - ScRange2(1)) * (ig2 - ScRange2(2))
        IF(ind .LE. 0) sigs2 = isodata%sm(ig, it2)%from(ig2)
        IsoXsMacSm(iso, ig2, ig) = IsoXsMacSm(iso, ig2, ig) + ND * (wt1 * sigs1 + wt2 * sigs2)
      ENDDO
      IsoXsMacS(iso, ig) =  IsoXsMacS(iso, ig) + ND * (wt2 * isodata%sigs(ig, it2) + wt1 * isodata%sigs(ig, it1))
      IsoXsMacStr(iso, ig) =  IsoXsMacStr(iso, ig) + ND * (wt2 * isodata%sigstr(ig, it2) + wt1 * isodata%sigstr(ig, it1))
    ENDDO
    NULLIFY(isodata)
  ENDDO

ENDIF

IF(lScat1) THEN
  DO iso = 1, niso
    DO ig = igb, ige
      IsoXsMacSm(iso, ig, ig) = IsoXsMacSm(iso, ig, ig) + (IsoXsMacS(iso, ig) - IsoXsMacStr(iso, ig))
    ENDDO
  ENDDO
ENDIF

END SUBROUTINE
SUBROUTINE IsoMacXsScatMatrix_Gen(IsoXsMacSm, temp, niso, idiso, pnum, igb, ige, ng, GroupInfo, lscat1)
! Resonance Scattering IS NOT CONSIDERED.
IMPLICIT NONE
INTEGER, INTENT(IN) :: niso, igb, ige, ng
INTEGER, INTENT(IN) :: idiso(niso)
REAL, INTENT(IN) :: pnum(niso)
REAL, INTENT(IN) :: temp
TYPE(GroupInfo_Type), INTENT(IN) :: GroupInfo
LOGICAL :: lscat1

TYPE(LIBDATA), POINTER :: isodata
INTEGER :: i, j, k
INTEGER :: iso, id, it1, it2, ig, ig2, ig3, ib, ie, ind
INTEGER :: ScRange1(2), ScRange2(2) !Scattering Range
REAL :: wt1, wt2, ND
REAL :: sigs1, sigs2
REAL :: IsoXsmacsm(niso,ng,ng),IsoXsMacS(niso,ng),IsoXsMacStr(niso,ng)

IsoXsMacSm(1:niso,1:ng,1:ng) = 0._8
IsoXsMacS(1:niso,1:ng) = 0._8
IsoXsMacStr(1:niso,1:ng) = 0._8

DO iso = 1, niso
  id = mapnucl(idiso(iso));   isodata => ldiso(id)
  ND = pnum(iso)
  CALL XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2)
  DO ig = igb, ige
    !Determine Scattering Range
    ScRange1 = (/isodata%sm(ig, it1)%ib,isodata%sm(ig, it1)%ie/)
    ScRange2 = (/isodata%sm(ig, it2)%ib,isodata%sm(ig, it2)%ie/)
    ib = min(ScRange1(1), ScRange2(1)); ie = max(ScRange1(2), ScRange2(2))
    DO ig2 = ib, ie  !In Scattering Group Sweep
      sigs1 = 0; sigs2 = 0
      ind = (ig2 - ScRange1(1)) * (ig2 - ScRange1(2))
      IF(ind .LE. 0) sigs1 = isodata%sm(ig, it1)%from(ig2)
      ind = (ig2 - ScRange2(1)) * (ig2 - ScRange2(2))
      IF(ind .LE. 0) sigs2 = isodata%sm(ig, it2)%from(ig2)
      IsoXsMacSm(iso, ig2, ig) = IsoXsMacSm(iso, ig2, ig) + ND * (wt1 * sigs1 + wt2 * sigs2)
    ENDDO
    IsoXsMacS(iso, ig) =  IsoXsMacS(iso, ig) + ND * (wt2 * isodata%sigs(ig, it2) + wt1 * isodata%sigs(ig, it1))
    IsoXsMacStr(iso, ig) =  IsoXsMacStr(iso, ig) + ND * (wt2 * isodata%sigstr(ig, it2) + wt1 * isodata%sigstr(ig, it1))
  ENDDO
  NULLIFY(isodata)
ENDDO


IF(lScat1) THEN
  DO iso = 1, niso
    DO ig = igb, ige
      IsoXsMacSm(iso, ig, ig) = IsoXsMacSm(iso, ig, ig) + (IsoXsMacS(iso, ig) - IsoXsMacStr(iso, ig))
    ENDDO
  ENDDO
ENDIF

END SUBROUTINE
SUBROUTINE MacP1XsScatMatrix_Csp(XsMac, FXR, igb, ige, ng, GroupInfo, lscat1, lScatSum, lCrCspFtn)
USE XsUtil_mod,        ONLY : GetXsMacDat,       ReturnXsMacDat
USE CrCsp_Mod,         ONLY : MacXsP1SmCsp
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
TYPE(FxrInfo_Type) :: Fxr
INTEGER, INTENT(IN) :: igb, ige, ng
TYPE(GroupInfo_Type), INTENT(IN) :: GroupInfo
LOGICAL :: lscat1, lscatsum, lCrCspFtn

TYPE(XsMac_Type), POINTER :: XsMac1, XsMac2
REAL :: temp
INTEGER :: niso
INTEGER, POINTER :: idiso(:)
REAL, POINTER ::pnum(:)
IF(lCrCspFtn .AND. Fxr%lCrCspFtn) THEN
  CALL GetXsMacDat(XsMac1, ng, .FALSE.)
  Call GetXsMacDat(XsMac2, ng, .FALSE.)
  niso = Fxr%CspFxr%niso(1); Temp = Fxr%Temp
  Pnum => Fxr%CspFxr%pnum(:, 1);   idiso => Fxr%CspFxr%isolist(:,1)
  CALL MacP1XsScatMatrix_Gen(XsMac1, Temp, niso, idiso, pnum, igb, ige, ng, GroupInfo)

  niso = Fxr%CspFxr%niso(2); Temp = Fxr%Temp
  Pnum => Fxr%CspFxr%pnum(:, 2);   idiso => Fxr%CspFxr%isolist(:,2)
  CALL MacP1XsScatMatrix_Gen(XsMac2, Temp, niso, idiso, pnum, igb, ige, ng, GroupInfo)
  CALL MacXsP1SmCsp(XsMac, XsMac1, XsMac2, Fxr%CspFxr, igb, ige, ng)
  !CALL MacXsNfCsp(XsMac, XsMac1, XsMac2, Fxr%CspFxr, ig1, ig2)
  NULLIFY(pnum); NULLIFY(idiso)
  CALL ReturnXsMacDat(XsMac1); CALL ReturnXsMacDat(XsMac2)
ELSE
  niso = Fxr%niso; temp = Fxr%Temp
  pnum => Fxr%pnum; idiso => Fxr%idiso
  CALL MacP1XsScatMatrix_Gen(XsMac, Temp, niso, idiso, pnum, igb, ige, ng, GroupInfo)
  NULLIFY(pnum); NULLIFY(idiso)
ENDIF
END SUBROUTINE


SUBROUTINE MacP1XsScatMatrix_Fxr(XsMac, FXR, igb, ige, ng, GroupInfo)
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
TYPE(FxrInfo_Type) :: Fxr
INTEGER, INTENT(IN) :: igb, ige, ng
TYPE(GroupInfo_Type), INTENT(IN) :: GroupInfo

REAL :: temp
INTEGER :: niso
INTEGER, POINTER :: idiso(:)
REAL, POINTER ::pnum(:)

niso = Fxr%niso; temp = Fxr%Temp
pnum => Fxr%pnum; idiso => Fxr%idiso
CALL MacP1XsScatMatrix_Gen(XsMac, Temp, niso, idiso, pnum, igb, ige, ng, GroupInfo)
NULLIFY(pnum); NULLIFY(idiso)
END SUBROUTINE

SUBROUTINE MacP1XsScatMatrix_Gen(XsMac, temp, niso, idiso, pnum, igb, ige, ng, GroupInfo)
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
INTEGER, INTENT(IN) :: niso, igb, ige, ng
INTEGER, INTENT(IN) :: idiso(niso)
REAL, INTENT(IN) :: pnum(niso)
REAL, INTENT(IN) :: temp
TYPE(GroupInfo_Type), INTENT(IN) :: GroupInfo

LOGICAL :: lscat1

TYPE(LIBDATA), POINTER :: isodata
INTEGER :: i, j, k
INTEGER :: iso, p1idx, id, it1, it2, ig, ig2, ib, ie, ind
INTEGER :: ScRange1(2), ScRange2(2) !Scattering Range
REAL :: wt1, wt2
REAL :: sigs1, sigs2
REAL, POINTER :: XsMacP1Sm(:,:)

IF(.NOT. XsMac%lalloc) THEN
  XsMac%ng = ng
  CALL AllocXsMac(XsMac)
ENDIF
!IF(.NOT. lAlloc) CALL AllocIsotopicMacXs()

XsMacP1Sm => XsMac%XsMacP1Sm

XsMacP1Sm(1:ng,1:ng) = 0.

!CALL CP_CA(XsMacS, 0._8, ng)
DO iso = 1, niso
  p1idx = mapnuclp13(idiso(iso))
  IF(p1idx .eq. 0) CYCLE
  id = idnp1hel(p1idx);   isodata => ldiso(id)
  CALL P1XsTempInterpolation(p1idx, isodata, temp, wt1, wt2, it1, it2)
  DO ig = igb, ige
    ScRange1 = (/isodata%smp1(ig, it1)%ib,isodata%smp1(ig, it1)%ie/)
    ScRange2 = (/isodata%smp1(ig, it2)%ib,isodata%smp1(ig, it2)%ie/)
    ib = min(ScRange1(1), ScRange2(1)); ie = max(ScRange1(2), ScRange2(2))
    DO ig2= ib, ie  !In Scattering Group Sweep
      sigs1 = 0; sigs2 = 0
      ind = (ig2 - ScRange1(1)) * (ig2 - ScRange1(2))
      IF(ind .LE. 0) sigs1 = isodata%smp1(ig, it1)%from(ig2)
      ind = (ig2 - ScRange2(1)) * (ig2 - ScRange2(2))
      IF(ind .LE. 0) sigs2 = isodata%smp1(ig, it2)%from(ig2)
      XsMacP1Sm(ig2, ig) = XsMacP1Sm(ig2, ig) + pnum(iso) * (wt1 * sigs1 + wt2 * sigs2)
    ENDDO
  ENDDO
  NULLIFY(isodata)
ENDDO

NULLIFY(XsMacP1Sm)

END SUBROUTINE

SUBROUTINE MacP1XsScatMatrix_Cusping(XsMac, iso, FxrInfo, igb, ige, ng, GroupInfo)
USE Material_Mod, ONLY : Mixture
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
INTEGER :: iso
TYPE(Fxrinfo_type) :: FxrInfo
INTEGER :: igb, ige, ng
TYPE(GroupInfo_Type) :: GroupInfo

REAL, POINTER :: pnum(:)
INTEGER, POINTER :: idiso(:)
REAL :: temp
INTEGER :: niso

niso = Mixture(iso)%niso
temp = FxrInfo%temp
idiso => Mixture(iso)%idiso
pnum => Mixture(iso)%pnum

CALL MacP1XsScatMatrix_Gen(XsMac, temp, niso, idiso, pnum, igb, ige, ng, GroupInfo)
NULLIFY(pnum, idiso)

END SUBROUTINE
SUBROUTINE MacP2XsScatMatrix_Cusping(XsMac, iso, FxrInfo, igb, ige, ng, GroupInfo)
USE Material_Mod, ONLY : Mixture
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
INTEGER :: iso
TYPE(Fxrinfo_type) :: FxrInfo
INTEGER :: igb, ige, ng
TYPE(GroupInfo_Type) :: GroupInfo

REAL, POINTER :: pnum(:)
INTEGER, POINTER :: idiso(:)
REAL :: temp
INTEGER :: niso

niso = Mixture(iso)%niso
temp = FxrInfo%temp
idiso => Mixture(iso)%idiso
pnum => Mixture(iso)%pnum

CALL MacP2XsScatMatrix_Gen(XsMac, temp, niso, idiso, pnum, igb, ige, ng, GroupInfo)
NULLIFY(pnum, idiso)

END SUBROUTINE
SUBROUTINE MacP3XsScatMatrix_Cusping(XsMac, iso, FxrInfo, igb, ige, ng, GroupInfo)
USE Material_Mod, ONLY : Mixture
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
INTEGER :: iso
TYPE(Fxrinfo_type) :: FxrInfo
INTEGER :: igb, ige, ng
TYPE(GroupInfo_Type) :: GroupInfo

REAL, POINTER :: pnum(:)
INTEGER, POINTER :: idiso(:)
REAL :: temp
INTEGER :: niso

niso = Mixture(iso)%niso
temp = FxrInfo%temp
idiso => Mixture(iso)%idiso
pnum => Mixture(iso)%pnum

CALL MacP3XsScatMatrix_Gen(XsMac, temp, niso, idiso, pnum, igb, ige, ng, GroupInfo)
NULLIFY(pnum, idiso)

END SUBROUTINE

SUBROUTINE MacP2XsScatMatrix_Csp(XsMac, FXR, igb, ige, ng, GroupInfo, lscat1, lScatSum, lCrCspFtn)
USE XsUtil_mod,        ONLY : GetXsMacDat,       ReturnXsMacDat
USE CrCsp_Mod,         ONLY : MacXsP2SmCsp
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
TYPE(FxrInfo_Type) :: Fxr
INTEGER, INTENT(IN) :: igb, ige, ng
TYPE(GroupInfo_Type), INTENT(IN) :: GroupInfo
LOGICAL :: lscat1, lscatsum, lCrCspFtn

TYPE(XsMac_Type), POINTER :: XsMac1, XsMac2
REAL :: temp
INTEGER :: niso
INTEGER, POINTER :: idiso(:)
REAL, POINTER ::pnum(:)
IF(lCrCspFtn .AND. Fxr%lCrCspFtn) THEN
  CALL GetXsMacDat(XsMac1, ng, .FALSE.)
  Call GetXsMacDat(XsMac2, ng, .FALSE.)
  niso = Fxr%CspFxr%niso(1); Temp = Fxr%Temp
  Pnum => Fxr%CspFxr%pnum(:, 1);   idiso => Fxr%CspFxr%isolist(:,1)
  CALL MacP2XsScatMatrix_Gen(XsMac1, Temp, niso, idiso, pnum, igb, ige, ng, GroupInfo)

  niso = Fxr%CspFxr%niso(2); Temp = Fxr%Temp
  Pnum => Fxr%CspFxr%pnum(:, 2);   idiso => Fxr%CspFxr%isolist(:,2)
  CALL MacP2XsScatMatrix_Gen(XsMac2, Temp, niso, idiso, pnum, igb, ige, ng, GroupInfo)
  CALL MacXsP2SmCsp(XsMac, XsMac1, XsMac2, Fxr%CspFxr, igb, ige, ng)
  !CALL MacXsNfCsp(XsMac, XsMac1, XsMac2, Fxr%CspFxr, ig1, ig2)
  NULLIFY(pnum); NULLIFY(idiso)
  CALL ReturnXsMacDat(XsMac1); CALL ReturnXsMacDat(XsMac2)
ELSE
  niso = Fxr%niso; temp = Fxr%Temp
  pnum => Fxr%pnum; idiso => Fxr%idiso
  CALL MacP2XsScatMatrix_Gen(XsMac, Temp, niso, idiso, pnum, igb, ige, ng, GroupInfo)
  NULLIFY(pnum); NULLIFY(idiso)
ENDIF
END SUBROUTINE


SUBROUTINE MacP2XsScatMatrix_Fxr(XsMac, FXR, igb, ige, ng, GroupInfo)
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
TYPE(FxrInfo_Type) :: Fxr
INTEGER, INTENT(IN) :: igb, ige, ng
TYPE(GroupInfo_Type), INTENT(IN) :: GroupInfo

REAL :: temp
INTEGER :: niso
INTEGER, POINTER :: idiso(:)
REAL, POINTER ::pnum(:)

niso = Fxr%niso; temp = Fxr%Temp
pnum => Fxr%pnum; idiso => Fxr%idiso
CALL MacP2XsScatMatrix_Gen(XsMac, Temp, niso, idiso, pnum, igb, ige, ng, GroupInfo)
NULLIFY(pnum); NULLIFY(idiso)
END SUBROUTINE

SUBROUTINE MacP2XsScatMatrix_Gen(XsMac, temp, niso, idiso, pnum, igb, ige, ng, GroupInfo)
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
INTEGER, INTENT(IN) :: niso, igb, ige, ng
INTEGER, INTENT(IN) :: idiso(niso)
REAL, INTENT(IN) :: pnum(niso)
REAL, INTENT(IN) :: temp
TYPE(GroupInfo_Type), INTENT(IN) :: GroupInfo

LOGICAL :: lscat1

TYPE(LIBDATA), POINTER :: isodata
INTEGER :: i, j, k
INTEGER :: iso, p2idx, id, it1, it2, ig, ig2, ib, ie, ind
INTEGER :: ScRange1(2), ScRange2(2) !Scattering Range
REAL :: wt1, wt2
REAL :: sigs1, sigs2
REAL, POINTER :: XsMacP2Sm(:,:)


IF(.NOT. XsMac%lalloc) THEN
  XsMac%ng = ng
  CALL AllocXsMac(XsMac)
ENDIF
!IF(.NOT. lAlloc) CALL AllocIsotopicMacXs()

XsMacP2Sm => XsMac%XsMacP2Sm


XsMacP2Sm(1:ng,1:ng) = 0.

!CALL CP_CA(XsMacS, 0._8, ng)
DO iso = 1, niso
  p2idx = mapnuclp13(idiso(iso))
  IF(p2idx .eq. 0) CYCLE
  id = idnp1hel(p2idx);   isodata => ldiso(id)
  CALL P2XsTempInterpolation(p2idx, isodata, temp, wt1, wt2, it1, it2)
  DO ig = igb, ige
    ScRange1 = (/isodata%smp2(ig, it1)%ib,isodata%smp2(ig, it1)%ie/)
    ScRange2 = (/isodata%smp2(ig, it2)%ib,isodata%smp2(ig, it2)%ie/)
    ib = min(ScRange1(1), ScRange2(1)); ie = max(ScRange1(2), ScRange2(2))
    DO ig2= ib, ie  !In Scattering Group Sweep
      sigs1 = 0; sigs2 = 0
      ind = (ig2 - ScRange1(1)) * (ig2 - ScRange1(2))
      IF(ind .LE. 0) sigs1 = isodata%smp2(ig, it1)%from(ig2)
      ind = (ig2 - ScRange2(1)) * (ig2 - ScRange2(2))
      IF(ind .LE. 0) sigs2 = isodata%smp2(ig, it2)%from(ig2)
      XsMacP2Sm(ig2, ig) = XsMacP2Sm(ig2, ig) + pnum(iso) * (wt1 * sigs1 + wt2 * sigs2)
    ENDDO
  ENDDO
  NULLIFY(isodata)
ENDDO

NULLIFY(XsMacP2Sm)

END SUBROUTINE

SUBROUTINE MacP3XsScatMatrix_Csp(XsMac, FXR, igb, ige, ng, GroupInfo, lscat1, lScatSum, lCrCspFtn)
USE XsUtil_mod,        ONLY : GetXsMacDat,       ReturnXsMacDat
USE CrCsp_Mod,         ONLY : MacXsP3SmCsp
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
TYPE(FxrInfo_Type) :: Fxr
INTEGER, INTENT(IN) :: igb, ige, ng
TYPE(GroupInfo_Type), INTENT(IN) :: GroupInfo
LOGICAL :: lscat1, lscatsum, lCrCspFtn

TYPE(XsMac_Type), POINTER :: XsMac1, XsMac2
REAL :: temp
INTEGER :: niso
INTEGER, POINTER :: idiso(:)
REAL, POINTER ::pnum(:)
IF(lCrCspFtn .AND. Fxr%lCrCspFtn) THEN
  CALL GetXsMacDat(XsMac1, ng, .FALSE.)
  Call GetXsMacDat(XsMac2, ng, .FALSE.)
  niso = Fxr%CspFxr%niso(1); Temp = Fxr%Temp
  Pnum => Fxr%CspFxr%pnum(:, 1);   idiso => Fxr%CspFxr%isolist(:,1)
  CALL MacP3XsScatMatrix_Gen(XsMac1, Temp, niso, idiso, pnum, igb, ige, ng, GroupInfo)

  niso = Fxr%CspFxr%niso(2); Temp = Fxr%Temp
  Pnum => Fxr%CspFxr%pnum(:, 2);   idiso => Fxr%CspFxr%isolist(:,2)
  CALL MacP3XsScatMatrix_Gen(XsMac2, Temp, niso, idiso, pnum, igb, ige, ng, GroupInfo)
  CALL MacXsP3SmCsp(XsMac, XsMac1, XsMac2, Fxr%CspFxr, igb, ige, ng)
  !CALL MacXsNfCsp(XsMac, XsMac1, XsMac2, Fxr%CspFxr, ig1, ig2)
  NULLIFY(pnum); NULLIFY(idiso)
  CALL ReturnXsMacDat(XsMac1); CALL ReturnXsMacDat(XsMac2)
ELSE
  niso = Fxr%niso; temp = Fxr%Temp
  pnum => Fxr%pnum; idiso => Fxr%idiso
  CALL MacP3XsScatMatrix_Gen(XsMac, Temp, niso, idiso, pnum, igb, ige, ng, GroupInfo)
  NULLIFY(pnum); NULLIFY(idiso)
ENDIF
END SUBROUTINE

SUBROUTINE MacP3XsScatMatrix_Fxr(XsMac, FXR, igb, ige, ng, GroupInfo)
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
TYPE(FxrInfo_Type) :: Fxr
INTEGER, INTENT(IN) :: igb, ige, ng
TYPE(GroupInfo_Type), INTENT(IN) :: GroupInfo

REAL :: temp
INTEGER :: niso
INTEGER, POINTER :: idiso(:)
REAL, POINTER ::pnum(:)

niso = Fxr%niso; temp = Fxr%Temp
pnum => Fxr%pnum; idiso => Fxr%idiso
CALL MacP3XsScatMatrix_Gen(XsMac, Temp, niso, idiso, pnum, igb, ige, ng, GroupInfo)
NULLIFY(pnum); NULLIFY(idiso)
END SUBROUTINE

SUBROUTINE MacP3XsScatMatrix_Gen(XsMac, temp, niso, idiso, pnum, igb, ige, ng, GroupInfo)
USE XSLIB_MOD, ONLY : idnp1hel, Mapnuclp13, LIBDATA, ldiso
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
INTEGER, INTENT(IN) :: niso, igb, ige, ng
INTEGER, INTENT(IN) :: idiso(niso)
REAL, INTENT(IN) :: pnum(niso)
REAL, INTENT(IN) :: temp
TYPE(GroupInfo_Type), INTENT(IN) :: GroupInfo

LOGICAL :: lscat1

TYPE(LIBDATA), POINTER :: isodata
INTEGER :: i, j, k
INTEGER :: iso, p3idx, id, it1, it2, ig, ig2, ib, ie, ind
INTEGER :: ScRange1(2), ScRange2(2) !Scattering Range
REAL :: wt1, wt2
REAL :: sigs1, sigs2
REAL, POINTER :: XsMacP3Sm(:,:)


IF(.NOT. XsMac%lalloc) THEN
  XsMac%ng = ng
  CALL AllocXsMac(XsMac)
ENDIF
!IF(.NOT. lAlloc) CALL AllocIsotopicMacXs()

XsMacP3Sm => XsMac%XsMacP3Sm


XsMacP3Sm(1:ng,1:ng) = 0.

!CALL CP_CA(XsMacS, 0._8, ng)
DO iso = 1, niso
  p3idx = mapnuclp13(idiso(iso))
  IF(p3idx .eq. 0) CYCLE
  id = idnp1hel(p3idx);   isodata => ldiso(id)
  CALL P3XsTempInterpolation(p3idx, isodata, temp, wt1, wt2, it1, it2)
  DO ig = igb, ige
    ScRange1 = (/isodata%smp3(ig, it1)%ib,isodata%smp3(ig, it1)%ie/)
    ScRange2 = (/isodata%smp3(ig, it2)%ib,isodata%smp3(ig, it2)%ie/)
    ib = min(ScRange1(1), ScRange2(1)); ie = max(ScRange1(2), ScRange2(2))
    DO ig2= ib, ie  !In Scattering Group Sweep
      sigs1 = 0; sigs2 = 0
      ind = (ig2 - ScRange1(1)) * (ig2 - ScRange1(2))
      IF(ind .LE. 0) sigs1 = isodata%smp3(ig, it1)%from(ig2)
      ind = (ig2 - ScRange2(1)) * (ig2 - ScRange2(2))
      IF(ind .LE. 0) sigs2 = isodata%smp3(ig, it2)%from(ig2)
      XsMacP3Sm(ig2, ig) = XsMacP3Sm(ig2, ig) + pnum(iso) * (wt1 * sigs1 + wt2 * sigs2)
    ENDDO
  ENDDO
  NULLIFY(isodata)
ENDDO

NULLIFY(XsMacP3Sm)

END SUBROUTINE

SUBROUTINE BaseMacSTr_Csp(XsMac, Fxr, ig, ng, lCrCspFtn)
USE TypeDef,        ONLY : FxrInfo_Type
USE XsUtil_mod,     ONLY : GetXsMacDat,          ReturnXsMacDat
USE CrCsp_Mod,      ONLY : MacXsStrCsp
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
TYPE(FxrInfo_Type) :: Fxr
INTEGER :: ig, ng
LOGICAL :: lCrCspFtn

TYPE(XsMac_Type), POINTER :: XsMac1, XsMac2
INTEGER :: niso
REAL :: Temp
REAL, POINTER :: pnum(:)
INTEGER, POINTER :: idiso(:)

IF(lCrCspFtn .AND. Fxr%lCrCspFtn) THEN
  CALL GetXsMacDat(XsMac1, ng, .FALSE.); CALL GetXsMacDat(XsMac2, ng, .FALSE.)
  niso = Fxr%CspFxr%niso(1); Temp = Fxr%Temp
  Pnum => Fxr%CspFxr%pnum(:, 1);   idiso => Fxr%CspFxr%isolist(:,1)
  CAll BaseMacSTr_Gen(XsMac1, temp, niso, idiso, pnum, ig, ng )

  niso = Fxr%CspFxr%niso(2); Temp = Fxr%Temp
  Pnum => Fxr%CspFxr%pnum(:, 2);   idiso => Fxr%CspFxr%isolist(:, 2)
  CAll BaseMacSTr_Gen(XsMac2, temp, niso, idiso, pnum, ig, ng )

  CALL MacXsStrCsp(XsMac, XsMac1, XsMac2, Fxr%CspFxr, ig, ig)
  CALL ReturnXsMacDat(XsMac1); CALL ReturnXsMacDat(XsMac2)
  NULLIFY(pnum, idiso)
ELSE
  CAll BaseMacSTr_Gen(XsMac, Fxr%temp, Fxr%niso, Fxr%idiso, Fxr%pnum, ig, ng )
ENDIF
END SUBROUTINE


SUBROUTINE BaseMacSTr_Fxr(XsMac, Fxr, ig, ng)
USE TypeDef,   ONLY : FxrInfo_Type
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
TYPE(FxrInfo_Type) :: Fxr
INTEGER :: ig, ng
CAll BaseMacSTr_Gen(XsMac, Fxr%temp, Fxr%niso, Fxr%idiso, Fxr%pnum, ig, ng )
END SUBROUTINE


SUBROUTINE BaseMacSTr_Gen(XsMac, temp, niso, idiso, pnum, ig, ng)
! Transport Corrected Xs Lib
! Resonance Scattering IS NOT CONSIDERED.
! XsMacStr is updated.
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
REAL, INTENT(IN) :: Temp
INTEGER, INTENT(IN) :: niso, ig, ng
INTEGER, INTENT(IN) :: idiso(niso)
REAL, INTENT(IN) :: pnum(niso)

TYPE(LIBDATA), POINTER :: isodata
INTEGER :: ig2, id, iso, it1, it2
INTEGER :: ioutbeg, ioutend, ind, ind2
INTEGER :: i, j, k
INTEGER :: ScRange1(2), ScRange2(2)
REAL :: wt1, wt2, sigs1, sigs2, sigssum, ND
REAL, POINTER :: XsMacstr(:)

IF(.NOT. XsMac%lalloc) THEN
  XsMac%ng = ng
  CALL AllocXsMac(XsMac)
ENDIF
!IF(.NOT. lAlloc) CALL AllocIsotopicMacXs()

XsMacStr => XsMac%XsMacStr
!XsMacStr(ig) = 0._8
SigsSum = 0
DO iso = 1, niso
  id = mapnucl(idiso(iso)); isodata => ldiso(id)
  CALL XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2)
  ScRange1 = (/isodata%sm(ig, it1)%ioutsb, isodata%sm(ig, it1)%ioutse/)
  ScRange2 = (/isodata%sm(ig, it2)%ioutsb, isodata%sm(ig, it2)%ioutse/)
!  ioutbeg = MIN(ScRange1(1), ScRange2(1)); ioutend = MAX(ScRange1(2), ScRange2(2))
  !Out Scattering Range
!  DO ig2 = ioutbeg, ioutend
!  sigs1 =0; sigs2 = 0
!  ind = (ig2 - ScRange2(1)) * (ig2 - ScRange2(2))
!  ind2 = (ig - IsoData%sm(ig2, it2)%ib) * (ig - IsoData%sm(ig2, it2)%ie)
!  IF(ind .le. 0 .and. ind2 .le. 0) sigs2 = isodata%sm(ig2, it2)%from(ig)
!  ind = (ig2 - ScRange1(1)) * (ig2 - ScRange1(2))
!  ind2 = (ig - IsoData%sm(ig2, it1)%ib) * (ig - IsoData%sm(ig2, it1)%ie)
!  IF(ind .le. 0 .and. ind2 .le. 0) sigs1 = isodata%sm(ig2, it1)%from(ig)
!  sigssum = sigssum + pnum(iso) * (wt1 * sigs1 + wt2 * sigs2)
!  ENDDO
  ND = pnum(iso)*wt2
  DO ig2 = ScRange2(1), ScRange2(2)
  ind2 = (ig - IsoData%sm(ig2, it2)%ib) * (ig - IsoData%sm(ig2, it2)%ie)
    IF (ind2 .LE. 0) sigssum = sigssum+ND*isodata%sm(ig2,it2)%from(ig)
  END DO
  ND = pnum(iso)*wt1
  DO ig2 = ScRange1(1), ScRange1(2)
  ind2 = (ig - IsoData%sm(ig2, it1)%ib) * (ig - IsoData%sm(ig2, it1)%ie)
    IF (ind2 .LE. 0) sigssum = sigssum+ND*isodata%sm(ig2,it1)%from(ig)
  ENDDO
ENDDO

XsMacStr(ig) = SigsSum

NULLIFY(XsMacStr); NULLIFY(isodata)
END SUBROUTINE

SUBROUTINE BaseMacStr_Cusping(XsMac, iso, Fxr, ig, ng)
USE Material_Mod,  ONLY: Mixture
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
TYPE(FxrInfo_Type) :: Fxr
INTEGER, INTENT(IN) :: ig, ng
INTEGER :: iso
REAL :: temp
INTEGER :: niso
INTEGER, POINTER :: idiso(:)
REAL, POINTER ::pnum(:)

niso = Mixture(iso)%niso; temp = Fxr%Temp
pnum => Mixture(iso)%pnum; idiso => Mixture(iso)%idiso
CALL BaseMacStr_Gen(XsMac, Temp, niso, idiso, pnum, ig, ng)
NULLIFY(pnum); NULLIFY(idiso)
END SUBROUTINE

SUBROUTINE EffRIFPin(mypin, TempAvgsq, iz, lAIC, PE)
USE TYPEDEF, ONLY : Pin_Type, ResVarPin_Type, PE_Type
USE CNTL, ONLY : nTracerCntl
USE OMP_LIB
USE Core_mod, ONLY : GroupInfo
IMPLICIT NONE
TYPE(PE_Type) :: PE
TYPE(ResVarPin_Type) :: mypin
TYPE(libdata), POINTER :: isodata,jsodata,repisodata
TYPE(riftype), POINTER :: isrif

INTEGER,INTENT(IN) :: iz
REAL,INTENT(IN) :: TempAvgsq
LOGICAL,INTENT(IN) :: lAIC

INTEGER :: iso, jso, kso, idres, jdres, kdres, id, jd, icat, repnid, repid, rifid, nThread, niso, nmlv, nmlv1g
INTEGER :: nlv, nlvflx, ilv, ig, ig1, ig2, it, ir, idxtemp(2), idxsig(2,2), idxrat(2), nRiTemp, npot, igt
INTEGER, POINTER :: idiso(:)

REAL :: siglpiso, micsigb, xsa, xsf, phi, lvabs, wgtabs, lvfis, wgtfis, lograt, logxsa, f, repria, ria, var(200)
REAL :: rifa(2), rifs(2), riff(2), rifa_iso, rifs_iso, riff_iso, rifrat, wgtsig(2,2), wgttemp(2), wgtrat(2)
REAL :: escxs(200), ind, jnd, adjintmlg, maclv_log, lvabs_log, wgtlvabs, phimult, phihom, siglpiso_dum
REAL, POINTER :: pnum(:), mlgmaclv(:)

ig1 = GroupInfo%iresoGrp1; ig2 = GroupInfo%iresoGrp2

niso = myPin%niso
pnum => myPin%pnum
idiso => myPin%idiso

!nThread = PE%nThread
!CALL OMP_SET_NUM_THREADS(nThread)

nmlv = mlgdata(iz)%f_nmaclv
nmlv1g = mlgdata(iz)%f_nmaclv1G

if (.not.nTracerCntl%lED) then
  if (.not.lAIC) then
    mlgmaclv => mlgdata(iz)%f_maclv_log
  else
    mlgmaclv => mlgdata(iz)%f_maclv1G_log
  endif
ELSE  ! EDIT JSU 2019-0405
  igt = mypin%igt
  mlgmaclv => mlgdata(iz)%f_maclv_pin_log(:,igt)
endif

myPin%rifa = 1._8
myPin%riff = 1._8
IF (nTracerCntl%lRST) myPin%rifs = 1._8

DO iso = 1, niso
  id = mapnucl(idiso(iso)); isodata => ldiso(id)
  idres = mapnuclres(idiso(iso),iz);
  IF ( .NOT. isodata%lreso ) CYCLE
  IF ( isodata%lclad ) CYCLE
  IF ( pnum(iso).le.0._8 ) CYCLE
  CALL calcWgt( TempAvgsq, isodata%rtempsq, isodata%nrtemp, wgttemp, idxtemp, 'A' ) ! index and wgt for temperature
  ind = myPin%pnum_res(myPin%idx_Res(iso))
  if (ind.le.0._8) cycle
  DO ig = ig1, ig2
    !siglpiso = ind * isodata%lamsigp(ig) ! JSU ver.
    siglpiso = pnum(iso)*isodata%lamsigp(ig) ! original ver.
    DO jso = 1, niso
      IF ( iso .EQ. jso ) CYCLE ! Same isotope
      jd = mapnucl(idiso(jso)); jsodata => ldiso(jd)
      !if (isodata%rifid(jd).ne.0 .AND. jsodata%lreso .AND. ( .NOT. jsodata%lclad )) CYCLE ! isotope in RIFL
      if (isodata%rifid(jd).ne.0 .AND. (jsodata%lreso .OR. ( .NOT. jsodata%lclad )) ) CYCLE
      siglpiso = siglpiso + pnum(jso) * jsodata%lamsigp(ig)
    ENDDO ! DO jso = 1, niso

    nlv = isodata%nlv; nlvflx = isodata%nlvflx
    ! Make Pin-averaged ALONE XS FOR interpolation of RIFL
    if (nTracerCntl%lMLG) then
      xsa = 0._8; xsf = 0._8; phi = 1._8; escxs = 0._8
      if (.not.lAIC) then
      adjintmlg = ind / myPin%FnAdj(ig)
        do ilv = 1, nlv
          lvabs = isodata%lvabs(ilv,ig)
          maclv_log = dlog(lvabs * adjintmlg)
          var(1:nmlv) = mypin%avgxseq_mg(1:nmlv,ig)
          escxs(ilv) = LineIntPol2( maclv_log, nmlv, mlgmaclv, var(1:nmlv) )
          micsigb = (escxs(ilv) + siglpiso) / ind
          !wgtabs = LineIntPol( TempAvgsq, isodata%nrtemp, isodata%rtempsq, isodata%wgtabs(ilv,:,ig) )
          wgtabs = wgttemp(2)*isodata%wgtabs(ilv,idxtemp(2),ig)+wgttemp(1)*isodata%wgtabs(ilv,idxtemp(1),ig)
          wgtlvabs = wgtabs * lvabs; phimult = 1._8 / (lvabs + micsigb); phihom = micsigb * phimult
          phi = phi - wgtlvabs * phimult
          xsa = xsa + wgtlvabs * phihom
          if ( isodata%ityp .NE. 3 ) cycle
          lvfis = isodata%lvfis(ilv,ig)
          !wgtfis = LineIntPol( TempAvgsq, isodata%nrtemp, isodata%rtempsq, isodata%wgtfis(ilv,:,ig) )
          wgtfis = wgttemp(2)*isodata%wgtfis(ilv,idxtemp(2),ig)+wgttemp(1)*isodata%wgtfis(ilv,idxtemp(1),ig)
          xsf = xsf + wgtfis * lvfis * phihom
        enddo
      else
        do ilv = 1, nlv
          lvabs = isodata%lvabs(ilv,ig)
          maclv_log = dlog(lvabs*ind)
          var(1:nmlv1g) = mypin%avgxseq_1g(1:nmlv1g)
          escxs(ilv) = LineIntPol2( maclv_log, nmlv1g, mlgmaclv, var(1:nmlv1g) )
          micsigb = (escxs(ilv) + siglpiso) / ind
          !wgtabs = LineIntPol( TempAvgsq, isodata%nrtemp, isodata%rtempsq, isodata%wgtabs(ilv,:,ig) )
          wgtabs = wgttemp(2)*isodata%wgtabs(ilv,idxtemp(2),ig)+wgttemp(1)*isodata%wgtabs(ilv,idxtemp(1),ig)
          wgtlvabs = wgtabs * lvabs; phimult = 1._8 / (lvabs + micsigb); phihom = micsigb * phimult
          phi = phi - wgtlvabs * phimult
          xsa = xsa + wgtlvabs * phihom
          if ( isodata%ityp .NE. 3 ) cycle
          lvfis = isodata%lvfis(ilv,ig)
          !wgtfis = LineIntPol( TempAvgsq, isodata%nrtemp, isodata%rtempsq, isodata%wgtfis(ilv,:,ig) )
          wgtfis = wgttemp(2)*isodata%wgtfis(ilv,idxtemp(2),ig)+wgttemp(2)*isodata%wgtfis(ilv,idxtemp(1),ig)
          xsf = xsf + wgtfis * lvfis * phihom
        enddo
      endif

    else ! .NOT. MLG
      if (nTracerCntl%lCAT) then  ! CATEGORY

        icat = isodata%icat
        repnid = ResoCat(icat)%repid; repid = mapnucl(repnid); repisodata => ldiso(repid)
        nRiTemp = repisodata%nrtemp; npot = repisodata%nsig0
        var(1:nRiTemp) = repisodata%ri_a(npot,ig, 1:nRiTemp)
        repria = LineIntpol(TempAvgsq, nRiTemp, repisodata%rtempsq(1:nRiTemp), var(1:nRiTemp) )
        nRiTemp = isodata%nrtemp; npot = isodata%nsig0
        var(1:nRiTemp) = isodata%ri_a(npot,ig, 1:nRiTemp)
        ria = LineIntpol(TempAvgsq, nRiTemp, isodata%rtempsq(1:nRiTemp), var(1:nRiTemp) )
        f = repria/ria
        escxs = 0._8
        if (nTracerCntl%l4Lv) then
          do ilv = 1, nlv
            lvabs_log = dlog(isodata%lvabs(ilv,ig) * f)
            var(1:repisodata%nlvflx) = mypin%avgxseq(1:repisodata%nlvflx,icat,ig)
            escxs(ilv) = LineIntPol2( lvabs_log, repisodata%nlvflx, repisodata%lvflx_log(1:repisodata%nlvflx,ig), var(1:repisodata%nlvflx) )
          enddo
        else
          do ilv = 1, nlv
            lvabs_log = dlog(isodata%lvabs(ilv,ig) * f)
            var(1:repisodata%nlv) = mypin%avgxseq(1:repisodata%nlv,icat,ig)
            escxs(ilv) = LineIntPol2( lvabs_log, repisodata%nlv, repisodata%lvabs_log(1:repisodata%nlv,ig), var(1:repisodata%nlv) )
          enddo
        endif

      else ! ISOTOPE

        escxs = 0._8
        if (nTracerCntl%l4Lv) then
          do ilv = 1, nlv
            lvabs_log = isodata%lvabs_log(ilv,ig)
            var(1:nlvflx) = mypin%avgxseq(1:nlvflx,idres,ig)
            escxs(ilv) = LineIntPol2( lvabs_log, nlvflx, isodata%lvflx_log(1:nlvflx,ig), var(1:nlvflx) )
          enddo
        else
          do ilv = 1, nlv
            escxs(ilv) = mypin%avgxseq(ilv,idres,ig)
          enddo
        endif

      endif

      xsa = 0._8; xsf = 0._8; phi = 1._8
      do ilv = 1, nlv
        micsigb = (escxs(ilv) + siglpiso) / ind
        lvabs = isodata%lvabs(ilv,ig)
        wgtabs = LineIntPol( TempAvgsq, isodata%nrtemp, isodata%rtempsq, isodata%wgtabs(ilv,:,ig) )
        wgtlvabs = wgtabs * lvabs; phimult = 1._8 / (lvabs + micsigb); phihom = micsigb * phimult
        phi = phi - wgtlvabs * phimult
        xsa = xsa + wgtlvabs * phihom
        if ( isodata%ityp .NE. 3 ) cycle
        lvfis = isodata%lvfis(ilv,ig)
        wgtfis = LineIntPol( TempAvgsq, isodata%nrtemp, isodata%rtempsq, isodata%wgtfis(ilv,:,ig) )
        xsf = xsf + wgtfis * lvfis * phihom
      enddo

    endif

    xsa = xsa / phi ! pin average xs for RI interpolation
    if (xsa.le.0._8) cycle
    if ( isodata%ityp .EQ. 3 ) xsf = xsf / phi
    logxsa = dlog(xsa)

    DO it = 1, 2
      var(1:isodata%nsig0) = isodata%xsalog(1:isodata%nsig0,ig,idxtemp(it))
      CALL calcWgt2( logxsa, var(1:isodata%nsig0) , isodata%nsig0, wgtsig(:,it), idxsig(:,it), 'D' ) ! idx and wgt for xsa (RI)
    ENDDO
    rifa= 0._8; rifs = 0._8; riff = 0._8
    do jso = 1, niso
      IF (jso.EQ.iso) CYCLE
      jd = mapnucl(idiso(jso))
      rifid = isodata%rifid(jd)
      if ( rifid.eq.0 ) CYCLE
      jnd = mypin%pnum_Res(mypin%idx_Res(jso))
      IF (jnd.LE.0.) CYCLE
      lograt = dlog(jnd/ind)
    !  PRINT *, idiso(iso), myPin%idiso_Res(jso), myPin%pnumrat(jso,iso)
      isrif => isodata%rif(rifid)
      ! Get N.D. ratio index and weight
      CALL calcWgt( lograt, isrif%ratlog, isrif%nrat, wgtrat, idxrat, 'A' ) ! idx and wgt for pnum ratio
      do it = 1, 2
        ! absorption
        rifa_iso = 0._8
        do ir = 1, 2
          rifrat = isrif%abs(idxrat(ir),idxsig(1,it),ig,idxtemp(it)) * wgtsig(1,it) + &
                   isrif%abs(idxrat(ir),idxsig(2,it),ig,idxtemp(it)) * wgtsig(2,it)
          rifa_iso = rifa_iso + rifrat * wgtrat(ir)
        enddo ! ir: pnum ratio
        rifa(it) = rifa(it) + ( rifa_iso - 1._8 )
        ! scattering
        if (nTracerCntl%lRST) then
          rifs_iso = 0._8
          do ir = 1, 2
            rifrat = isrif%sct(idxrat(ir),idxsig(1,it),ig,idxtemp(it)) * wgtsig(1,it) + &
                     isrif%sct(idxrat(ir),idxsig(2,it),ig,idxtemp(it)) * wgtsig(2,it)
            rifs_iso = rifs_iso + rifrat * wgtrat(ir)
          enddo ! ir: pnum ratio
          rifs(it) = rifs(it) + ( rifs_iso - 1._8 )
        endif
        if ( isodata%ityp .NE. 3 ) CYCLE
        ! fission
        riff_iso = 0._8
        do ir = 1, 2
          rifrat = isrif%fis(idxrat(ir),idxsig(1,it),ig,idxtemp(it)) * wgtsig(1,it) + &
                   isrif%fis(idxrat(ir),idxsig(2,it),ig,idxtemp(it)) * wgtsig(2,it)
          riff_iso = riff_iso + rifrat * wgtrat(ir)
        enddo! ir: pnum ratio
        riff(it) = riff(it) + ( riff_iso - 1._8 )
      enddo ! do it=1,2
    enddo ! DO jso = 1, niso
    rifa = rifa + 1._8
    riff = riff + 1._8
    myPin%rifa(idres,ig) = rifa(1) * wgttemp(1) + rifa(2) * wgttemp(2)
    myPin%riff(idres,ig) = riff(1) * wgttemp(1) + riff(2) * wgttemp(2)
    if (nTracerCntl%lRST) then
      rifs = rifs + 1._8
      myPin%rifs(idres,ig) = rifs(1) * wgttemp(1) + rifs(2) * wgttemp(2)
    endif
  ENDDO ! DO ig = ig1, ig2
ENDDO ! DO iso = 1, niso

NULLIFY(isodata, jsodata, repisodata)
NULLIFY(isrif, pnum, idiso)
END SUBROUTINE

SUBROUTINE EffMacXS(XsMac, mypin, myFxr, tempref, niso, ig, ng, lIsoXsOut, iz, PE, ifxr)
! Inside of OpenMP
USE TYPEDEF, ONLY : FxrInfo_Type, Pin_Type, ResVarPin_Type, PE_Type
USE CNTL, ONLY : nTracerCntl
USE OMP_LIB
USE PointXSRT_MOD, ONLY : neg, TempIntIdxWgtPSM, icgb, icge, delu_pw
!USE XSLIB_MOD,     ONLY : pwxs_type, pwxs
#ifdef __PGI
USE IEEE_ARITHMETIC   !--- CNJ Edit : F2003 Standard
#endif
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac     !Microscopic XS
TYPE(PE_Type) :: PE
TYPE(ResVarPin_Type) :: mypin
TYPE(FxrInfo_Type) :: myFxr
TYPE(libdata), POINTER :: isodata, jsodata, repisodata

INTEGER,INTENT(IN) :: ig, niso, ng, iz, ifxr
REAL,INTENT(IN) :: tempref
LOGICAL,INTENT(IN) :: lIsoXsOut

INTEGER :: iso, jso, id, jd, idres, jdres, repnid, repid, icat, idxsig(2,2), ir, kso, idxrat(2), kdres
INTEGER :: nlv, nlvflx, ilv, it, it1, it2, idxtemp(2), idxtempavg(2), npot, nRiTemp, nmaclv, ibon, nThread
INTEGER :: rifid, igt, iufg, idpxs
INTEGER,PARAMETER :: nbonmax = 10000
INTEGER,POINTER :: idiso(:)

REAL :: temp, TempAvgsq, Tempsq, wt1, wt2, wgttemp(2), wgttempavg(2), invFnAdj
REAL :: siglp, siglp_noreso, bgxs, ind, lvabs, wgtabs, xsa, phi, logxsa, xss, xsss, xsf, lvfis, wgtfis
REAL :: miclv, miclv_log, lvabs_log, maclv_log, maclv, adjintmlg, Ft, dum2, f, repria, ria, rs, rifa(2), rifs(2), riff(2)
REAL :: wgtrat(2), jnd, rifa_iso, rifs_iso, riff_iso, rifrat, lograt
REAL :: SumSigA, SigaAnother, phimult, phihom, fl2, fl1, invfl1, wgtsig(2,2), s(2), wgtlvabs, var(200), flavg
REAL,DIMENSION(niso, ng) :: IsoXsMacA, IsoXsMacS, IsoXsMacSS, IsoXsMacS1, IsoXsMacF, IsoXSMacNF, IsoMacNu, IsoXsMacKF, IsoMacK
REAL,DIMENSION(niso) :: isoxsmaca_bon, isoxsmacf_bon, xsmaca0, phiiso, ratind, wphi
REAL,DIMENSION(nlvmax,niso) :: FtIso, micsigb, macsigb, wgtabs_bon, wgtfis_bon, maclvabs, maclvfis
REAL,POINTER :: xsmaca(:), xsmacf(:), xsmacs(:), xsmacstr(:), pnum(:), mlgmaclv(:), xsmacnf(:), xsmackf(:)
REAL(4),POINTER :: xseq(:)

LOGICAL :: lOTFRIFISO
INTEGER :: ipgs, ipge, idx1, idx2
REAL :: wgt1, wgt2, siglptot, sigptot, sigp_noreso, sigp, Ftiso_p, sumflx_ufgmix, sumflx_ufgiso, sumRR_ufgmix, sumRR_ufgiso, rif_temp, Presoesc, lambda, delMLG
REAL,POINTER, DIMENSION(:) :: PXST_MIX, PXSA_MIX, PXSF_MIX, PXSS_MIX, PXST_ISO, PXSA_ISO, PXSF_ISO, PXSS_ISO, UFGFLX_MIX, UFGFLX_ISO, PXSLS_MIX, PXSLS_ISO
TYPE(pwxs_type), pointer :: pxs
TYPE(riftype), POINTER :: isrif

IF(.NOT. XsMac%lalloc) THEN
  XsMac%ng = ng
  CALL AllocXsMac(XsMac)
ENDIF

IF( lIsoXsOut .AND. ( .NOT. XsMac%lIsoAlloc) ) THEN
  CALL AllocMacIsoXs(XsMac, ng, nelthel) ! 19/09/24 using nelthel considering increasable niso |--Edit by JSU
ENDIF

IF (nTRACERCntl%lOTFRIF .AND. .NOT.XsMac%lAllocpxs)   CALL AllocPXsMac(XsMac)
!IF (nTRACERCntl%lOTFRIF .AND. nTRACERCntl%lRST .AND. .NOT.XsMac%lAllocpxss) CALL AllocPXsSMac(XsMac)

TempAvgsq = tempref
temp = myFxr%temp
Tempsq = dsqrt(temp)
if (myFxr%lCLD) TempAvgsq=Tempsq

pnum => myFxr%pnum; idiso => myFxr%idiso

xsmaca => XsMac%xsmaca; xsmacf => XsMac%xsmacf;     xsmacnf => XsMac%xsMacnF
xsmacs => XsMac%xsmacs; xsmacstr => XsMac%xsmacstr; xsmackf => XsMac%xsmackf

IsoXsMacF(1:niso,ig) = 0.
IsoXsMacnF(1:niso,ig) = 0.
IsoMacNu(1:niso,ig) = 0.
IsoXsMacKF(1:niso,ig) = 0.
IsoMacK(1:niso,ig) = 0.

DO iso = 1, niso
  id = mapnucl(idiso(iso));   isodata => ldiso(id)
  CALL XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2)
  IsoXsMacA(iso, ig) = pnum(iso) * (wt2 * isodata%siga(ig, it2) + wt1 * isodata%siga(ig, it1))
  IsoXsMacS(iso, ig) = pnum(iso) * (wt2 * isodata%sigs(ig, it2) + wt1 * isodata%sigs(ig, it1))
  IsoXsMacS1(iso, ig) = IsoXsMacS(iso, ig) - pnum(iso) * (wt2 * isodata%sigstr(ig, it2) + wt1 * isodata%sigstr(ig, it1))
  IsoXsMacSS(iso, ig) = pnum(iso) * (wt2 * isodata%sigss(ig, it2) + wt1 * isodata%sigss(ig, it1))
  IF ( isodata%ifis .NE. 0 ) THEN
    IsoXsMacF(iso, ig) = pnum(iso) * (wt2 * isodata%sigf(ig, it2) + wt1 * isodata%sigf(ig, it1))
    IsoXsMacNF(iso, ig) = pnum(iso) * (wt2 * isodata%signf(ig, it2) + wt1 * isodata%signf(ig, it1))
    IsoMacK(iso, ig) = isodata%kappa
    IsoXsMacKF(iso, ig) = IsoXsMacF(iso, ig) * isodata%kappa
    IsoMacNu(iso, ig) = IsoXsMacNF(iso, ig) / IsoXsMacF(iso, ig)
  END IF
ENDDO

IF (nTracerCNTL%lRIF) THEN

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RIF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FtIso(1:nlvmax,1:niso) = 1._8
  micsigb(1:nlvmax,1:niso) = 0._8

  siglp_noreso = 0._8; sigp_noreso = 0._8
  DO jso = 1, niso
    jd = mapnucl(idiso(jso)); jsodata => ldiso(jd)
    IF ( .not.jsodata%lreso ) THEN
      siglp_noreso = siglp_noreso + pnum(jso) * jsodata%lamsigp(ig)
      sigp_noreso = sigp_noreso + pnum(jso) * jsodata%sigp
    END IF
  ENDDO
  ! SUBGROUP EFFECTIVE XS
  IF (nTracerCNTL%lMLG) THEN
    !******************************** MLG ********************************!
    if (.not.nTracerCNTL%lED) THEN
      if (myFxr%lCLD) then
        mlgmaclv => mlgdata0%c_maclv1G_log
        xseq => myFxr%XsEq_c_1g
        nmaclv = mlgdata0%c_nmaclv1G
        invFnAdj = 1._8
        delMLG = mlgdata0%del_c1g
      elseif (myFxr%lAIC) then
        mlgmaclv => mlgdata(iz)%f_maclv1G_log
        xseq => myFxr%XsEq_f_1G
        nmaclv = mlgdata(iz)%f_nmaclv1G
        invFnAdj = 1._8
        delMLG = mlgdata(iz)%del_f1g
      else
        mlgmaclv => mlgdata(iz)%f_maclv_log
        xseq => myFxr%XsEq_f_mg(:,ig)
        nmaclv = mlgdata(iz)%f_nmaclv
        invFnAdj = 1._8 / myFxr%FnAdj(ig)
        delMLG = mlgdata(iz)%del_fmg
      endif
    else ! equivalent Dancoff...
      if (myFxr%lCLD) then
        mlgmaclv => mlgdata0%c_maclv1G_log
        xseq => myFxr%XsEq_c_1g
        nmaclv = mlgdata0%c_nmaclv1G
        invFnAdj = 1._8
        delMLG = mlgdata0%del_c1g
      elseif (myFxr%lAIC) then
        igt = mypin%igt
        mlgmaclv => mlgdata(iz)%f_maclv_pin_log(:,igt)
        xseq => myFxr%XsEq_f_1G
        nmaclv = mlgdata(iz)%f_nmaclv1G
        invFnAdj = 1._8
        delMLG = mlgdata(iz)%del_f1g
      else
        igt = mypin%igt
        mlgmaclv => mlgdata(iz)%f_maclv_pin_log(:,igt)
        xseq => myFxr%XsEq_f_mg(:,ig)
        nmaclv = mlgdata(iz)%f_nmaclv
        invFnAdj = 1._8 / myFxr%FnAdj(ig)
        delMLG = mlgdata(iz)%del_fmg
      endif
    endif

    DO iso = 1, niso
      id = mapnucl(idiso(iso)); isodata => ldiso(id)
      IF (.NOT.isodata%lreso) CYCLE
      idres = mapnuclRes(idiso(iso),iz)
      ind = myFXR%pnum_Res(myFXR%idx_Res(iso));
      siglp = siglp_noreso + ind * isodata%lamsigp(ig)
      IF (ind.le.0._8) CYCLE
      nlv = isodata%nlv
      IF (.NOT.myFxr%lCLD.AND..not.myFxr%lAIC) then
        adjintmlg = ind * invFnAdj
        do ilv = 1, nlv
          maclv_log = dlog(isodata%lvabs(ilv,ig) * adjintmlg)
          var(1:nmaclv) = xseq(1:nmaclv)
          bgxs = LineIntPol2( maclv_log, nmaclv, mlgmaclv, var(1:nmaclv) )
          var(1:nmaclv) = myFxr%FtAdj(:,ig)
          FtIso(ilv,iso) = LineIntPol2( maclv_log, nmaclv, mlgmaclv, var(1:nmaclv) )
          micsigb(ilv,iso) = (bgxs + siglp) / ind / FtIso(ilv,iso)
        enddo
      ELSE
        do ilv = 1, nlv
          maclv_log = dlog(isodata%lvabs(ilv,ig) * ind)
          var(1:nmaclv) = xseq(1:nmaclv)
          bgxs = LineIntPol2( maclv_log, nmaclv, mlgmaclv, var(1:nmaclv) )
          FtIso(ilv,iso) = 1._8
          micsigb(ilv,iso) = (bgxs + siglp) / ind
        enddo
      ENDIF
    ENDDO ! DO iso = 1, niso
    ! OTF RIF MIXTURE XS AND FLUX
    IF (nTRACERCntl%lOTFRIF) THEN
      ipgs = icgb(ig)
      ipge = icge(ig)
      PXST_MIX => XsMac%PXST_MIX;  PXSA_MIX => XsMac%PXSA_MIX;  PXSF_MIX => XsMac%PXSF_MIX;  PXSS_MIX => XsMac%PXSS_MIX; PXSLS_MIX => XsMac%PXSLS_MIX
      PXST_ISO => XsMac%PXST_ISO;  PXSA_ISO => XsMac%PXSA_ISO;  PXSF_ISO => XsMac%PXSF_ISO;  PXSS_ISO => XsMac%PXSS_ISO; PXSLS_ISO => XsMac%PXSLS_ISO
      UFGFLX_MIX => XsMac%UFGFLX_MIX;  UFGFLX_ISO => XsMac%UFGFLX_ISO
      ! set mixture pxs
      !PXST_MIX = 0._8
      PXSA_MIX(ipgs:ipge) = 0._8;      PXSF_MIX(ipgs:ipge) = 0._8;      PXSS_MIX(ipgs:ipge) = 0._8; PXSLS_MIX = 0._8;
      siglptot = 0._8; sigptot = 0._8;
      DO iso = 1, niso
        id = mapnucl(idiso(iso));   isodata => ldiso(id)
        idpxs = isodata%pxsid
        sigptot = sigptot + pnum(iso)*isodata%sigp ! NR
        siglptot = siglptot + pnum(iso)*isodata%lamsigp(ig) ! IR
        IF (idpxs .EQ. 0) THEN ! NUCLIDE NOT IN PXS LIBRARY
          PXSS_MIX(ipgs:ipge) = PXSS_MIX(ipgs:ipge) + pnum(iso)*isodata%sigp
          PXSLS_MIX(ipgs:ipge) = PXSLS_MIX(ipgs:ipge) + pnum(iso)*isodata%lamsigp(ig)
        ELSE ! NUCLIDE IN PXS LIBRARY
          pxs => pwxs(idpxs)
          CALL TempIntIdxWgtPSM(pxs, temp, it1, it2, wgt1, wgt2)
          PXSA_MIX(ipgs:ipge) = PXSA_MIX(ipgs:ipge) + pnum(iso)*(pxs%xsa(ipgs:ipge,it1)*wgt1+pxs%xsa(ipgs:ipge,it2)*wgt2)
          PXSS_MIX(ipgs:ipge) = PXSS_MIX(ipgs:ipge) + pnum(iso)*(pxs%xss(ipgs:ipge,it1)*wgt1+pxs%xss(ipgs:ipge,it2)*wgt2)
          PXSLS_MIX(ipgs:ipge) = PXSLS_MIX(ipgs:ipge) + pnum(iso)*(pxs%xss(ipgs:ipge,it1)*wgt1+pxs%xss(ipgs:ipge,it2)*wgt2)*isodata%lamsigp(ig)/isodata%sigp
          IF (pxs%lfis) PXSF_MIX(ipgs:ipge) = PXSF_MIX(ipgs:ipge) + pnum(iso)*(pxs%xsf(ipgs:ipge,it1)*wgt1+pxs%xsf(ipgs:ipge,it2)*wgt2)
        END IF
      END DO
      !DO iufg = ipgs, ipge
      !  IF (PXSLS_MIX(iufg) .LT. 0._8) PXSLS_MIX(iufg) = 0._8
      !END DO
      PXST_MIX(ipgs:ipge) = PXSA_MIX(ipgs:ipge) + PXSS_MIX(ipgs:ipge)
      Presoesc = 1._8
      DO iufg = ipgs, ipge
        maclv_log = dlog((PXSA_MIX(iufg)+PXSLS_MIX(iufg)) * invFnAdj)
        !maclv_log = dlog((PXSA_MIX(iufg)) * invFnAdj)
        Presoesc = Presoesc * exp(-PXSA_MIX(iufg)/PXST_MIX(iufg)*delu_pw(iufg,ig)*0.5)
        var(1:nmaclv) = xseq(1:nmaclv)
        !bgxs = LineIntPol2( maclv_log, nmaclv, mlgmaclv, var(1:nmaclv) ) + siglptot
        bgxs = LineIntPol2( maclv_log, nmaclv, mlgmaclv, var(1:nmaclv) )
        !UFGFLX_MIX(iufg) = (bgxs  + siglptot) / (PXSA_MIX(iufg) + bgxs + PXSLS_MIX(iufg)) * Presoesc
        UFGFLX_MIX(iufg) = (bgxs  + siglptot) / (PXSA_MIX(iufg) + bgxs + siglptot) * Presoesc
        !UFGFLX_MIX(iufg) = bgxs / (PXSA_MIX(iufg) + bgxs) * Presoesc
        Presoesc = Presoesc * exp(-PXSA_MIX(iufg)/PXST_MIX(iufg)*delu_pw(iufg,ig)*0.5)
#ifdef DEBUG
        WRITE(98, '(I5, I8, A4, I7, 5ES14.6)') ifxr, iufg, ' MIX', 0, PXSA_MIX(iufg), PXST_MIX(iufg), Presoesc, UFGFLX_MIX(iufg), bgxs
#endif
      END DO
      sumflx_ufgmix = DOT_PRODUCT(UFGFLX_MIX(ipgs:ipge),delu_pw(ipgs:ipge,ig))
    END IF

  ELSE

    IF (nTracerCNTL%lCAT) THEN
      !******************************** CAT ********************************!
      DO iso = 1, niso
        id = mapnucl(idiso(iso)); isodata => ldiso(id)
        IF ( .NOT. isodata%lreso ) CYCLE
        idres = mapnuclRes(idiso(iso),iz)
        ind = myFXR%pnum_Res(myFXR%idx_Res(iso))
        siglp = siglp_noreso + ind * isodata%lamsigp(ig)
        IF (ind.le.0._8) CYCLE
        nlv = isodata%nlv; nlvflx = isodata%nlvflx
        ! Category information
        icat = isodata%icat
        repnid = ResoCat(icat)%repid; repid = mapnucl(repnid); repisodata => ldiso(repid)
        nRiTemp = repisodata%nrtemp; npot = repisodata%nsig0
        var(1:nRiTemp) = repisodata%ri_a(npot,ig, 1:nRiTemp)
        repria = LineIntpol(TempAvgsq, nRiTemp, repisodata%rtempsq(1:nRiTemp), var(1:nRiTemp) )
        nRiTemp = isodata%nrtemp; npot = isodata%nsig0
        var(1:nRiTemp) = isodata%ri_a(npot,ig, 1:nRiTemp)
        ria = LineIntpol(TempAvgsq, nRiTemp, isodata%rtempsq(1:nRiTemp), var(1:nRiTemp) )
        f = repria/ria
        do ilv = 1, nlv
          lvabs_log = dlog(isodata%lvabs(ilv,ig) * f)
          if (nTracerCNTL%l4Lv) then
            var(1:repisodata%nlvflx) = myFxr%xseq(1:repisodata%nlvflx,icat,ig)
            bgxs = LineIntPol2( lvabs_log, repisodata%nlvflx, repisodata%lvflx_log(1:repisodata%nlvflx,ig), var(1:repisodata%nlvflx) )
            var(1:repisodata%nlvflx) = myFxr%NDAF(1:repisodata%nlvflx,icat,ig)
            FtIso(ilv,iso) = LineIntPol2( lvabs_log, repisodata%nlvflx, repisodata%lvflx_log(1:repisodata%nlvflx,ig), var(1:repisodata%nlvflx) )
          else
            var(1:repisodata%nlv) = myFxr%xseq(1:repisodata%nlv,icat,ig)
            bgxs = LineIntPol2( lvabs_log, repisodata%nlv, repisodata%lvabs_log(1:repisodata%nlv,ig), var(1:repisodata%nlv) )
            var(1:repisodata%nlv) = myFxr%NDAF(1:repisodata%nlv,icat,ig)
            FtIso(ilv,iso) = LineIntPol2( lvabs_log, repisodata%nlv, repisodata%lvabs_log(1:repisodata%nlv,ig), var(1:repisodata%nlv) )
          endif
          micsigb(ilv,iso) = (bgxs + siglp) / ind / FtIso(ilv,iso)
        enddo
      ENDDO ! DO iso = 1, niso

    ELSE
      !******************************** ISO ********************************!
      DO iso = 1, niso
        id = mapnucl(idiso(iso)); isodata => ldiso(id)
        IF (.NOT.isodata%lreso) CYCLE
        idres = mapnuclres(idiso(iso),iz)
        ind = myFXR%pnum_Res(myFXR%idx_Res(iso))
        siglp = siglp_noreso + ind * isodata%lamsigp(ig)
        IF (ind.le.0._8) CYCLE
        nlv = isodata%nlv; nlvflx = isodata%nlvflx
        do ilv = 1, nlv
          miclv_log = isodata%lvabs_log(ilv,ig)
          if (nTracerCNTL%l4Lv) then
              var(1:nlvflx) = myFxr%xseq(1:nlvflx,idres,ig)
              bgxs = LineIntPol2( miclv_log, nlvflx, isodata%lvflx_log(1:nlvflx,ig), var(1:nlvflx) )
              var(1:nlvflx) = myFxr%NDAF(1:nlvflx,idres,ig)
              FtIso(ilv,iso) = LineIntPol2( miclv_log, nlvflx, isodata%lvflx_log(1:nlvflx,ig), var(1:nlvflx) )
          else
              bgxs = myFxr%xseq(ilv,idres,ig)
              FtIso(ilv,iso) = myFxr%NDAF(ilv,idres,ig)
          endif
          micsigb(ilv,iso) = (bgxs + siglp) / ind / FtIso(ilv,iso)
        enddo
      ENDDO

    ENDIF

  ENDIF

  DO iso = 1, niso
    id = mapnucl(idiso(iso)); isodata => ldiso(id)
    idpxs = isodata%pxsid
    IF (.NOT.isodata%lreso) CYCLE
    IF (pnum(iso).le.0._8) CYCLE
    lOTFRIFISO = (nTRACERCntl%lOTFRIF.AND..NOT.(idpxs.EQ.0))
    idres = mapnuclres(idiso(iso),iz)
    ind = myFxr%pnum_Res(myFxr%idx_Res(iso))
    IF (ind.le.0._8) CYCLE
    ! Calculation of Effective XS of nuclides alone...
    xsa = 0._8; xsf = 0._8; phi = 1._8
    nlv = isodata%nlv
    call calcWgt(TempAvgsq, isodata%rtempsq, isodata%nrtemp, wgttempavg, idxtempavg, 'A')
    do ilv = 1, nlv
      lvabs = isodata%lvabs(ilv,ig)
      !wgtabs = LineIntPol( TempAvgsq, isodata%nrtemp, isodata%rtempsq, isodata%wgtabs(ilv,:,ig) )
      wgtabs = wgttempavg(2)*isodata%wgtabs(ilv,idxtempavg(2),ig)+wgttempavg(1)*isodata%wgtabs(ilv,idxtempavg(1),ig)
!      IF (myFxr%FsrIdxSt.EQ.97 .AND. ig.EQ.10 .AND. iso.EQ.2) print*, 'W', wgtabs, lvabs
      wgtlvabs = wgtabs * lvabs; phimult = 1._8 / (lvabs + micsigb(ilv,iso)); phihom = FtIso(ilv,iso) * micsigb(ilv,iso) * phimult
      phi = phi - wgtlvabs * phimult
      xsa = xsa + wgtlvabs * phihom
      if ( isodata%ityp .NE. 3 ) cycle
      lvfis = isodata%lvfis(ilv,ig)
      !wgtfis = LineIntPol( TempAvgsq, isodata%nrtemp, isodata%rtempsq, isodata%wgtfis(ilv,:,ig) )
      wgtfis = wgttempavg(2)*isodata%wgtfis(ilv,idxtempavg(2),ig)+wgttempavg(1)*isodata%wgtfis(ilv,idxtempavg(1),ig)
      xsf = xsf + wgtfis * lvfis * phihom
    enddo
    xsa = xsa / phi; xsf = xsf / phi
    if (xsa.le.0._8) cycle

    call calcWgt( Tempsq, isodata%rtempsq, isodata%nrtemp, wgttemp, idxtemp, 'A')
    ! RIF CALCULATION...
    IF (.NOT.lOTFRIFISO.OR..NOT.nTRACERCntl%lMLG) THEN
    IF (.not.myFxr%lCLD) THEN
      IF (.NOT. nTRACERCntl%lRIFFXR) THEN ! Pinwise RIF
        IsoXsMacA(iso, ig) = xsa * pnum(iso) * myPin%rifa(idres,ig)
        if ( isodata%ityp .EQ. 3 ) then
          !xsf = xsf / phi
          IsoXsMacF(iso, ig) = xsf * pnum(iso) * myPin%riff(idres,ig)
          IsoXsMacNF(iso, ig) = IsoXsMacF(iso, ig) * IsoMacNu(iso, ig)
          IsoXsMacKF(iso, ig) = IsoXsMacF(iso, ig) * IsoMacK(iso, ig)
        endif
        IF (nTracerCNTL%lRST) THEN
          logxsa = dlog(xsa)
          DO it = 1, 2
            var(1:isodata%nsig0) = isodata%xsalog(1:isodata%nsig0,ig,idxtemp(it))
            CALL calcWgt2( logxsa, var(1:isodata%nsig0), isodata%nsig0, wgtsig(:,it), idxsig(:,it), 'D' )
            s(it) = ( isodata%xss(idxsig(1,it), ig, idxtemp(it)) ** wgtsig(1, it) ) * ( isodata%xss(idxsig(2,it), ig, idxtemp(it)) ** wgtsig(2, it) )
          ENDDO
          xss = s(1) * wgttemp(1) + s(2) * wgttemp(2)
#ifdef __PGI
          if (ieee_is_nan(xss)) cycle
#else
          if (isnan(xss)) cycle
#endif
          xsss = xss * isodata%ssf(ig)
          IsoXsMacS(iso, ig) = xss * pnum(iso) * myPin%rifs(idres,ig)
          IsoXsMacS1(iso, ig) = IsoXsMacS(iso, ig) * isodata%mu
          IsoXsMacSS(iso, ig) = xsss * pnum(iso) * myPin%rifs(idres,ig)
        ENDIF
      ELSE ! FXR-WISE RIF
        logxsa = dlog(xsa)
        DO it = 1, 2
          var(1:isodata%nsig0) = isodata%xsalog(1:isodata%nsig0,ig,idxtemp(it))
          CALL calcWgt2( logxsa, var(1:isodata%nsig0) , isodata%nsig0, wgtsig(:,it), idxsig(:,it), 'D' )
        END DO
        rifa= 0._8; rifs = 0._8; riff = 0._8
        do jso = 1, myFXR%niso
          IF (iso.EQ.jso) CYCLE
          jdres = mapnuclres(myFXR%idiso(jso),iz)
          jd = mapnucl(myFXR%idiso(jso))
          rifid = isodata%rifid(jd)
          if ( rifid.eq.0 ) CYCLE
          jnd = myFxr%pnum_Res(myFxr%idx_Res(jso))
          IF (jnd.LE.0.) CYCLE
          lograt = dlog(jnd/ind) ! log of pnum ratio
          isrif => isodata%rif(rifid)
          ! Get N.D. ratio index and weight
          CALL calcWgt( lograt, isrif%ratlog, isrif%nrat, wgtrat, idxrat, 'A' )
          do it = 1, 2
            ! absorption
            rifa_iso = 0._8
            do ir = 1, 2
              rifrat = isrif%abs(idxrat(ir),idxsig(1,it),ig,idxtemp(it)) * wgtsig(1,it) + &
                       isrif%abs(idxrat(ir),idxsig(2,it),ig,idxtemp(it)) * wgtsig(2,it)
              rifa_iso = rifa_iso + rifrat * wgtrat(ir)
            enddo    ! ir (pnum ratio)
            rifa(it) = rifa(it) + ( rifa_iso - 1._8 )
            ! scattering
            if (nTracerCntl%lRST) then
              rifs_iso = 0._8
              do ir = 1, 2
                rifrat = isrif%sct(idxrat(ir),idxsig(1,it),ig,idxtemp(it)) * wgtsig(1,it) + &
                         isrif%sct(idxrat(ir),idxsig(2,it),ig,idxtemp(it)) * wgtsig(2,it)
                rifs_iso = rifs_iso + rifrat * wgtrat(ir)
              enddo! ir (pnum ratio)
              rifs(it) = rifs(it) + ( rifs_iso - 1._8 )
            endif
            ! fission
            if (isodata%ityp .EQ. 3) THEN
            riff_iso = 0._8
            do ir = 1, 2
              rifrat = isrif%fis(idxrat(ir),idxsig(1,it),ig,idxtemp(it)) * wgtsig(1,it) + &
                       isrif%fis(idxrat(ir),idxsig(2,it),ig,idxtemp(it)) * wgtsig(2,it)
              riff_iso = riff_iso + rifrat * wgtrat(ir)
            enddo! ir (pnum ratio)
            riff(it) = riff(it) + ( riff_iso - 1._8 )
            END IF
          enddo ! do it=1,2
        enddo ! DO jso = 1, niso

        rifa = rifa + 1._8
        riff = riff + 1._8
        IsoXsMacA(iso, ig) = xsa * pnum(iso) * (rifa(1) * wgttemp(1) + rifa(2) * wgttemp(2))
#ifdef DEBUG
        WRITE(99, '(2I7, I3,A3,5ES14.6)') ifxr, idiso(iso), ig, 'A', (rifa(1) * wgttemp(1) + rifa(2) * wgttemp(2))
#endif
        if ( isodata%ityp .EQ. 3 ) then
          !xsf = xsf / phi
          IsoXsMacF(iso, ig) = xsf * pnum(iso) * (riff(1) * wgttemp(1) + riff(2) * wgttemp(2))
          IsoXsMacNF(iso, ig) = IsoXsMacF(iso, ig) * IsoMacNu(iso, ig)
          IsoXsMacKF(iso, ig) = IsoXsMacF(iso, ig) * IsoMacK(iso, ig)
#ifdef DEBUG
        WRITE(99, '(2I7, I3,A3,5ES14.6)') ifxr, idiso(iso), ig, 'F', (rifa(1) * wgttemp(1) + rifa(2) * wgttemp(2))
#endif
        endif
      END IF
    else ! cladding
      logxsa = dlog(xsa)
      DO it = 1, 2
        var(1:isodata%nsig0) = isodata%xsalog(1:isodata%nsig0,ig,idxtemp(it))
        CALL calcWgt2( logxsa, var(1:isodata%nsig0) , isodata%nsig0, wgtsig(:,it), idxsig(:,it), 'D' )
      ENDDO

      rifa= 0._8; rifs = 0._8; riff = 0._8
      do jso = 1, niso
        if ( jso .EQ. iso ) cycle ! exclude itself
        jd = mapnucl(idiso(jso)); jdres = mapnuclres(idiso(jso),iz)
        rifid = isodata%rifid(jd)
        if ( rifid.eq.0 ) CYCLE
        jnd = myFxr%pnum_Res(myFxr%idx_Res(jso))
        IF (jnd .LE. 0) CYCLE
        lograt = dlog(jnd/ind)
        isrif => isodata%rif(rifid)
        ! Get N.D. ratio index and weight
        CALL calcWgt( lograt, isrif%ratlog, isrif%nrat, wgtrat, idxrat, 'A' )
        do it = 1, 2
          ! absorption
          rifa_iso = 0._8
          do ir = 1, 2
            rifrat = isrif%abs(idxrat(ir),idxsig(1,it),ig,idxtemp(it)) * wgtsig(1,it) + &
                     isrif%abs(idxrat(ir),idxsig(2,it),ig,idxtemp(it)) * wgtsig(2,it)
            rifa_iso = rifa_iso + rifrat * wgtrat(ir)
          enddo
          rifa(it) = rifa(it) + ( rifa_iso - 1._8 )
          ! scattering
          if (nTracerCntl%lRST) then
            rifs_iso = 0._8
            do ir = 1, 2
              rifrat = isrif%sct(idxrat(ir),idxsig(1,it),ig,idxtemp(it)) * wgtsig(1,it) + &
                       isrif%sct(idxrat(ir),idxsig(2,it),ig,idxtemp(it)) * wgtsig(2,it)
              rifs_iso = rifs_iso + rifrat * wgtrat(ir)
            enddo
            rifs(it) = rifs(it) + ( rifs_iso - 1._8 )
          endif
          if ( isodata%ityp .NE. 3 ) CYCLE ! cladding..?
          ! fission
          riff_iso = 0._8
          do ir = 1, 2
            rifrat = isrif%fis(idxrat(ir),idxsig(1,it),ig,idxtemp(it)) * wgtsig(1,it) + &
                     isrif%fis(idxrat(ir),idxsig(2,it),ig,idxtemp(it)) * wgtsig(2,it)
            riff_iso = riff_iso + rifrat * wgtrat(ir)
          enddo
          riff(it) = riff(it) + ( riff_iso - 1._8 )
        enddo ! do it=1,2
      enddo ! DO jso = 1, niso
      rifa = rifa + 1._8
      riff = riff + 1._8
      IsoXsMacA(iso, ig) = xsa * pnum(iso) * (rifa(1) * wgttemp(1) + rifa(2) * wgttemp(2))
      if ( isodata%ityp .EQ. 3 ) then
        !xsf = xsf / phi
        IsoXsMacF(iso, ig) = xsf * pnum(iso) * (riff(1) * wgttemp(1) + riff(2) * wgttemp(2))
        IsoXsMacNF(iso, ig) = IsoXsMacF(iso, ig) * IsoMacNu(iso, ig)
        IsoXsMacKF(iso, ig) = IsoXsMacF(iso, ig) * IsoMacK(iso, ig)
      endif
      IF (nTracerCNTL%lRST) THEN
        DO it = 1, 2
          s(it) = ( isodata%xss(idxsig(1,it), ig, idxtemp(it)) ** wgtsig(1, it) ) * ( isodata%xss(idxsig(2,it), ig, idxtemp(it)) ** wgtsig(2, it) )
        ENDDO
        xss = s(1) * wgttemp(1) + s(2) * wgttemp(2)
#ifdef __PGI
        if (ieee_is_nan(xss)) cycle
#else
        if (isnan(xss)) cycle
#endif
        xsss = xss * isodata%ssf(ig)
        rifs = rifs + 1._8
        rs = rifs(1) * wgttemp(1) + rifs(2) * wgttemp(2)
        IsoXsMacS(iso, ig) = xss * pnum(iso) * rs
        IsoXsMacS1(iso, ig) = IsoXsMacS(iso, ig) * isodata%mu
        IsoXsMacSS(iso, ig) = xsss * pnum(iso) * rs
      ENDIF
    ENDIF ! RIFL (FUEL OR CLAD)
    ELSE ! OTFRIF
      pxs => pwxs(idpxs)
      !ind = myFXR%pnum_Res(myFXR%idx_Res(iso))
      ind = pnum(iso)
      if (IND .LT. 1E-15) CYCLE
      siglp = siglp_noreso + ind * isodata%lamsigp(ig)
      sigp  = sigp_noreso  + ind * isodata%sigp
      CALL TempIntIdxWgtPSM(pxs, temp, it1, it2, wgt1, wgt2)
      PXSA_ISO(ipgs:ipge) = (pxs%xsa(ipgs:ipge,it1)*wgt1+pxs%xsa(ipgs:ipge,it2)*wgt2)*pnum(iso)
      PXSS_ISO(ipgs:ipge) = (pxs%xss(ipgs:ipge,it1)*wgt1+pxs%xss(ipgs:ipge,it2)*wgt2)*pnum(iso)+siglp
      PXSLS_ISO(ipgs:ipge) = (pxs%xss(ipgs:ipge,it1)*wgt1+pxs%xss(ipgs:ipge,it2)*wgt2)*pnum(iso)*isodata%lamsigp(ig)/isodata%sigp
      !DO iufg = ipgs, ipge
      !  IF (PXSLS_ISO(iufg) .LT. 0._8) PXSLS_ISO(iufg) = 0._8
      !END DO
      IF (pxs%lfis) PXSF_ISO(ipgs:ipge)= (pxs%xsf(ipgs:ipge,it1)*wgt1+pxs%xsf(ipgs:ipge,it2)*wgt2)*pnum(iso)
      PXST_ISO(ipgs:ipge) = PXSA_ISO(ipgs:ipge) + PXSS_ISO(ipgs:ipge)
      Presoesc = 1._8
      IF (myFxr%lCLD) THEN
        DO iufg = ipgs, ipge
          Presoesc = Presoesc * exp(-PXSA_ISO(iufg)/PXST_ISO(iufg)*delu_pw(iufg,ig)*0.5)
          maclv_log = dlog((PXSA_ISO(iufg)+PXSLS_ISO(iufg)) * invFnAdj) ! *ind
          !maclv_log = dlog((PXSA_ISO(iufg)) * invFnAdj) ! *ind
          var(1:nmaclv) = xseq(1:nmaclv)
          bgxs = LineIntPol2( maclv_log, nmaclv, mlgmaclv, var(1:nmaclv) )
          Ftiso_p = 1._8
          !bgxs = (bgxs + siglp)! / ind
          !bgxs = (bgxs + sigp)! / ind
!          UFGFLX_ISO(iufg) = (bgxs + siglp) / (PXSA_ISO(iufg) + bgxs + PXSLS_ISO(iufg) + siglp_noreso) * Presoesc
          UFGFLX_ISO(iufg) = (bgxs + siglp) / (PXSA_ISO(iufg) + bgxs + siglp) * Presoesc
          !UFGFLX_ISO(iufg) = bgxs / (PXSA_ISO(iufg) + bgxs) * Presoesc
          Presoesc = Presoesc * exp(-PXSA_ISO(iufg)/PXST_ISO(iufg)*delu_pw(iufg,ig)*0.5)
#ifdef DEBUG
          WRITE(98, '(I5, I8, A4, i7, 5ES14.6)') ifxr, iufg, ' ISO', idiso(iso), PXSA_ISO(iufg), PXST_ISO(iufg), Presoesc, UFGFLX_ISO(iufg), bgxs
#endif
        END DO
      ELSE
        DO iufg = ipgs, ipge
          Presoesc = Presoesc * exp(-PXSA_ISO(iufg)/PXST_ISO(iufg)*delu_pw(iufg,ig)*0.5)
          maclv_log = dlog((PXSA_ISO(iufg)+PXSLS_ISO(iufg)) * invFnAdj) ! *ind
          !maclv_log = dlog((PXSA_ISO(iufg)) * invFnAdj) ! *ind
          !CALL LineIntPol2_2Var(maclv_log, bgxs, Ftiso_p, nmaclv, mlgmaclv, xseq(1:nmaclv), myFxr%FtAdj(:,ig))
          var(1:nmaclv) = xseq(1:nmaclv)
          bgxs = LineIntPol2( maclv_log, nmaclv, mlgmaclv, var(1:nmaclv) )
          Ftiso_p = 1._8
          !bgxs = (bgxs + siglp)! / Ftiso_p !  / ind
          !bgxs = (bgxs + sigp) / Ftiso_p !  / ind
          !UFGFLX_ISO(iufg) = (bgxs+siglp) / (PXSA_ISO(iufg) + bgxs + PXSLS_ISO(iufg) + siglp_noreso) * Presoesc
          UFGFLX_ISO(iufg) = (bgxs + siglp) / (PXSA_ISO(iufg) + bgxs + siglp) * Presoesc
          !UFGFLX_ISO(iufg) = bgxs / (PXSA_ISO(iufg) + bgxs) * Presoesc
          Presoesc = Presoesc * exp(-PXSA_ISO(iufg)/PXST_ISO(iufg)*delu_pw(iufg,ig)*0.5)
#ifdef DEBUG
          WRITE(98, '(I5, I8, A4, i7, 5ES14.6)') ifxr, iufg, ' ISO', idiso(iso), PXSA_ISO(iufg), PXST_ISO(iufg), Presoesc, UFGFLX_ISO(iufg), bgxs
#endif
        END DO
      END IF
      sumflx_ufgiso = DOT_PRODUCT(UFGFLX_ISO(ipgs:ipge),delu_pw(ipgs:ipge,ig))
      ! absorption
      sumRR_ufgmix = DOT_PRODUCT(PXSA_ISO(ipgs:ipge)*UFGFLX_MIX(ipgs:ipge),delu_pw(ipgs:ipge,ig))
      sumRR_ufgiso = DOT_PRODUCT(PXSA_ISO(ipgs:ipge)*UFGFLX_ISO(ipgs:ipge),delu_pw(ipgs:ipge,ig))
      rif_temp = sumRR_ufgmix / sumRR_ufgiso * sumflx_ufgiso / sumflx_ufgmix
#ifdef DirectUSERI
      IsoXsMacA(iso, ig) = sumRR_ufgmix / sumflx_ufgmix * ind
#else
      IsoXsMacA(iso, ig) = xsa * rif_temp * ind
#endif
#ifdef DEBUG
      WRITE(96, '(I5, I7, I3,A3,5ES14.6)') ifxr, idiso(iso), ig, 'A', rif_temp, sumRR_ufgmix, sumRR_ufgiso, sumflx_ufgmix, sumflx_ufgiso
#endif
      !IsoXsMacA(iso, ig) = sumRR_ufgmix / sumflx_ufgmix * ind
      ! fission
      if ( isodata%ityp .EQ. 3 ) then
        !xsf = xsf / phi
        IF (pxs%lfis) THEN
          sumRR_ufgmix = DOT_PRODUCT(PXSF_ISO(ipgs:ipge)*UFGFLX_MIX(ipgs:ipge),delu_pw(ipgs:ipge,ig))
          sumRR_ufgiso = DOT_PRODUCT(PXSF_ISO(ipgs:ipge)*UFGFLX_ISO(ipgs:ipge),delu_pw(ipgs:ipge,ig))
          rif_temp = sumRR_ufgmix / sumRR_ufgiso * sumflx_ufgiso / sumflx_ufgmix
#ifdef DirectUSERI
          IsoXsMacF(iso, ig) = sumRR_ufgmix / sumflx_ufgmix * ind
#else
          IsoXsMacF(iso, ig) = xsf * rif_temp * ind
#endif
          IsoXsMacNF(iso, ig) = IsoXsMacF(iso, ig) * IsoMacNu(iso, ig)
#ifdef DEBUG
          WRITE(96, '(I5, I7, I3,A3,5ES14.6)') ifxr, idiso(iso), ig, 'F', rif_temp, sumRR_ufgmix, sumRR_ufgiso, sumflx_ufgiso, sumflx_ufgmix
#endif
        !  IsoXsMacF(iso, ig) = sumRR_ufgmix / sumflx_ufgmix * ind
        !  IsoXsMacNF(iso, ig) = IsoXsMacF(iso, ig) * IsoMacNu(iso, ig)
        ELSE
          IsoXsMacF(iso, ig) = xsf * ind
          IsoXsMacNF(iso, ig) = IsoXsMacF(iso, ig) * IsoMacNu(iso, ig)
          !WRITE(96, '(I5, I7, I3,A3,ES14.6)') ifxr, idiso(iso), ig, 'F', 1.0, 0.0, 0.0, 0.0, 0.0
        END IF
      endif

    END IF ! RIF (RIFL OR OTF)
  ENDDO ! iso = 1, niso

ELSE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BON !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  isoxsmaca_bon(1:niso) = IsoXsMacA(1:niso, ig)
  isoxsmacf_bon(1:niso) = IsoXsMacF(1:niso, ig)
  macsigb(1:nlvmax,1:niso) = 0._8
  maclvabs(1:nlvmax,1:niso) = 0._8
  wgtabs_bon(1:nlvmax,1:niso) = 0._8
  maclvfis(1:nlvmax,1:niso) = 0._8
  wgtfis_bon(1:nlvmax,1:niso) = 0._8
  ratind(1:niso) = 1._8

  siglp = 0._8
  DO jso = 1, niso
    jd = mapnucl(idiso(jso)); jsodata => ldiso(jd)
    siglp = siglp + pnum(jso) * jsodata%lamsigp(ig)
  ENDDO ! DO jso = 1, niso

  SumSigA = 0._8
  DO iso = 1, niso
    id = mapnucl(idiso(iso)); isodata => ldiso(id)
    IF ( .NOT. isodata%lreso ) CYCLE
    SumSigA = SumSigA + isoxsmaca_bon(iso)
  enddo ! DO iso = 1, niso
  DO iso = 1, niso
    id = mapnucl(idiso(iso)); isodata => ldiso(id)
    IF ( .NOT. isodata%lreso ) CYCLE
    wphi(iso) = IsoXsMacA(iso, ig) / SumSigA
  END DO

  IF (nTracerCNTL%lMLG) THEN
    !******************************** MLG ********************************!
    if (.not.nTracerCNTL%lED) THEN
      if (myFxr%lCLD) then
        mlgmaclv => mlgdata0%c_maclv1G_log
        xseq => myFxr%XsEq_c_1g
        nmaclv = mlgdata0%c_nmaclv1G
      elseif (myFxr%lAIC) then
        mlgmaclv => mlgdata(iz)%f_maclv1G_log
        xseq => myFxr%XsEq_f_1G
        nmaclv = mlgdata(iz)%f_nmaclv1G
      else
        mlgmaclv => mlgdata(iz)%f_maclv_log
        xseq => myFxr%XsEq_f_mg(:,ig)
        nmaclv = mlgdata(iz)%f_nmaclv
        invFnAdj = 1._8 / myFxr%FnAdj(ig)
      endif
    else
      if (myFxr%lCLD) then
        mlgmaclv => mlgdata0%c_maclv1G_log
        xseq => myFxr%XsEq_c_1g
        nmaclv = mlgdata0%c_nmaclv1G
      elseif (myFxr%lAIC) then
        igt = mypin%igt
        mlgmaclv => mlgdata(iz)%f_maclv_pin_log(:,igt)
        xseq => myFxr%XsEq_f_1G
        nmaclv = mlgdata(iz)%f_nmaclv1G
      else
        igt = mypin%igt
        mlgmaclv => mlgdata(iz)%f_maclv_pin_log(:,igt)
        xseq => myFxr%XsEq_f_mg(:,ig)
        nmaclv = mlgdata(iz)%f_nmaclv
        invFnAdj = 1._8 / myFxr%FnAdj(ig)
      endif
    endif

    DO iso = 1, niso
      id = mapnucl(idiso(iso)); isodata => ldiso(id)
      IF ( .NOT. isodata%lreso ) CYCLE
      idres = mapnuclRes(idiso(iso),iz)
      ind = 0._8
      DO jso = 1, niso
        jdres = mapnuclres(idiso(jso),iz)
        IF (idres.eq.jdres) ind = ind + pnum(jso)
      ENDDO ! DO jso = 1, niso
      IF (ind.le.0._8) CYCLE
      ratind(iso)=pnum(iso)/ind
      nlv = isodata%nlv; nlvflx = isodata%nlvflx

      IF (.NOT.myFxr%lCLD.and..NOT.myFxr%lAIC) then
        adjintmlg = ind * invFnAdj
        do ilv = 1, nlv
          maclv_log = dlog(isodata%lvabs(ilv,ig) * adjintmlg)
          var(1:nmaclv) = xseq(1:nmaclv)
          bgxs = LineIntPol2( maclv_log, nmaclv, mlgmaclv, var(1:nmaclv) )
          macsigb(ilv,iso) = bgxs + siglp
          var(1:nmaclv) = myFxr%FtAdj(1:nmaclv,ig)
          Ft = LineIntPol2( maclv_log, nmaclv, mlgmaclv, var(1:nmaclv) )
          dum2 = Ft * ind
          maclvabs(ilv,iso) =  dum2 * isodata%lvabs(ilv,ig)
          wgtabs_bon(ilv,iso) = LineIntPol( TempAvgsq, isodata%nrtemp, isodata%rtempsq, isodata%wgtabs(ilv,:,ig) )
          if ( isodata%ityp .EQ. 3 ) then
            maclvfis(ilv,iso) = dum2 * isodata%lvfis(ilv,ig)
            wgtfis_bon(ilv,iso) = LineIntPol( TempAvgsq, isodata%nrtemp, isodata%rtempsq, isodata%wgtfis(ilv,:,ig) )
          endif
        enddo
      ELSE
        do ilv = 1, nlv
          maclv = isodata%lvabs(ilv,ig) * ind
          maclv_log = dlog(maclv)
          var(1:nmaclv) = xseq(1:nmaclv)
          bgxs = LineIntPol2( maclv_log, nmaclv, mlgmaclv, var(1:nmaclv) )
          macsigb(ilv,iso) = bgxs + siglp
          maclvabs(ilv,iso) = maclv
          wgtabs_bon(ilv,iso) = LineIntPol( TempAvgsq, isodata%nrtemp, isodata%rtempsq, isodata%wgtabs(ilv,:,ig) )
          if ( isodata%ityp .EQ. 3 ) then
            maclvfis(ilv,iso) = isodata%lvfis(ilv,ig) * ind
            wgtfis_bon(ilv,iso) = LineIntPol( TempAvgsq, isodata%nrtemp, isodata%rtempsq, isodata%wgtfis(ilv,:,ig) )
          endif
        enddo
      ENDIF
    ENDDO ! DO iso = 1, niso

  ELSE

    IF (nTracerCNTL%lCAT) THEN
      !******************************** CAT ********************************!
      DO iso = 1, niso
        id = mapnucl(idiso(iso)); isodata => ldiso(id)
        IF ( .NOT. isodata%lreso ) CYCLE
        idres = mapnuclRes(idiso(iso),iz)
        ind = 0._8
        DO jso = 1, niso
          jdres = mapnuclres(idiso(jso),iz)
          IF (idres.eq.jdres) ind = ind + pnum(jso)
        ENDDO ! DO jso = 1, niso
        IF (ind.le.0._8) CYCLE
        ratind(iso)=pnum(iso)/ind
        nlv = isodata%nlv; nlvflx = isodata%nlvflx
        icat = isodata%icat
        repnid = ResoCat(icat)%repid; repid = mapnucl(repnid); repisodata => ldiso(repid)
        nRiTemp = repisodata%nrtemp; npot = repisodata%nsig0
        var(1:nRiTemp) = repisodata%ri_a(npot,ig, 1:nRiTemp)
        repria = LineIntpol(TempAvgsq, nRiTemp, repisodata%rtempsq(1:nRiTemp), var(1:nRiTemp))
        nRiTemp = isodata%nrtemp; npot = isodata%nsig0
        var(1:nRiTemp) = isodata%ri_a(npot,ig, 1:nRiTemp)
        ria = LineIntpol(TempAvgsq, nRiTemp, isodata%rtempsq(1:nRiTemp), var(1:nRiTemp))
        f = repria/ria
        do ilv = 1, nlv
          lvabs_log = dlog(isodata%lvabs(ilv,ig) * f)
          if (nTracerCNTL%l4Lv) then
            var(1:repisodata%nlvflx) = myFxr%xseq(1:repisodata%nlvflx,icat,ig)
            bgxs = LineIntPol2( lvabs_log, repisodata%nlvflx, repisodata%lvflx_log(1:repisodata%nlvflx,ig), var(1:repisodata%nlvflx) )
            var(1:repisodata%nlvflx) = myFxr%NDAF(1:repisodata%nlvflx,icat,ig)
            Ft = LineIntPol2( lvabs_log, repisodata%nlvflx, repisodata%lvflx_log(1:repisodata%nlvflx,ig), var(1:repisodata%nlvflx) )
          else
            var(1:repisodata%nlv) = myFxr%xseq(1:repisodata%nlv,icat,ig)
            bgxs = LineIntPol2( lvabs_log, repisodata%nlv, repisodata%lvabs_log(1:repisodata%nlv,ig), var(1:repisodata%nlv) )
            var(1:repisodata%nlv) = myFxr%NDAF(1:repisodata%nlv,icat,ig)
            Ft = LineIntPol2( lvabs_log, repisodata%nlv, repisodata%lvabs_log(1:repisodata%nlv,ig), var(1:repisodata%nlv) )
          endif
          macsigb(ilv,iso) = bgxs + siglp
          dum2 = Ft * ind
          maclvabs(ilv,iso) = dum2 * isodata%lvabs(ilv,ig)
          wgtabs_bon(ilv,iso) = LineIntPol( TempAvgsq, isodata%nrtemp, isodata%rtempsq, isodata%wgtabs(ilv,:,ig) )
          if ( isodata%ityp .EQ. 3 ) then
            maclvfis(ilv,iso) = dum2 * isodata%lvfis(ilv,ig)
            wgtfis_bon(ilv,iso) = LineIntPol( TempAvgsq, isodata%nrtemp, isodata%rtempsq, isodata%wgtfis(ilv,:,ig) )
          endif
        enddo
      ENDDO ! DO iso = 1, niso
    ELSE
      !******************************** ISO ********************************!
      DO iso = 1, niso
        id = mapnucl(idiso(iso)); isodata => ldiso(id)
        IF ( .NOT. isodata%lreso ) CYCLE
        idres = mapnuclRes(idiso(iso),iz)
        ind = 0._8
        DO jso = 1, niso
          jdres = mapnuclres(idiso(jso),iz)
          IF (idres.eq.jdres) ind = ind + pnum(jso)
        ENDDO ! DO jso = 1, niso
        IF (ind.le.0._8) CYCLE
        ratind(iso)=pnum(iso)/ind
        nlv = isodata%nlv; nlvflx = isodata%nlvflx

        do ilv = 1, nlv
          miclv = isodata%lvabs(ilv,ig)
          miclv_log = dlog(miclv)
          if (nTracerCNTL%l4Lv) then
            var(1:nlvflx) = myFxr%xseq(1:nlvflx,idres,ig)
            bgxs = LineIntPol2( miclv_log, nlvflx, isodata%lvflx_log(1:nlvflx,ig), var(1:nlvflx))
            var(1:nlvflx) = myFxr%NDAF(1:nlvflx,idres,ig)
            Ft = LineIntPol2( miclv_log, nlvflx, isodata%lvflx_log(1:nlvflx,ig), var(1:nlvflx))
          else
            bgxs = myFxr%xseq(ilv,idres,ig)
            Ft = myFxr%NDAF(ilv,idres,ig)
          endif
          macsigb(ilv,iso) = bgxs + siglp
          dum2 = Ft * ind
          maclvabs(ilv,iso) = dum2 * miclv
          wgtabs_bon(ilv,iso) = LineIntPol( TempAvgsq, isodata%nrtemp, isodata%rtempsq, isodata%wgtabs(ilv,:,ig) )
          if ( isodata%ityp .EQ. 3 ) then
            maclvfis(ilv,iso) = dum2 * isodata%lvfis(ilv,ig)
            wgtfis_bon(ilv,iso) = LineIntPol( TempAvgsq, isodata%nrtemp, isodata%rtempsq, isodata%wgtfis(ilv,:,ig) )
          endif
        enddo
      ENDDO ! DO iso = 1, niso
    ENDIF

  ENDIF

    fl2 = 1._8
    DO ibon = 1, nbonmax

      DO iso = 1, niso
        xsmaca0(iso) = isoxsmaca_bon(iso)
      ENDDO

      fl1 = 1._8
      DO iso = 1, niso
        id = mapnucl(idiso(iso)); isodata => ldiso(id)
        IF ( .NOT. isodata%lreso ) CYCLE
        IF ( pnum(iso).le.0._8 ) CYCLE
        nlv = isodata%nlv

        phiiso(iso) = 0._8
        SigaAnother =  SumSigA - xsmaca0(iso)
        isoxsmaca_bon(iso) = 0._8
        if ( isodata%ityp .EQ. 3 ) isoxsmacf_bon(iso) = 0._8

        do ilv = 1, nlv
          phimult = 1._8 / (maclvabs(ilv,iso) + SigaAnother + macsigb(ilv,iso))
          phihom = macsigb(ilv,iso) * phimult
          wgtabs = wgtabs_bon(ilv, iso) * maclvabs(ilv,iso)
          phiiso(iso) = phiiso(iso) + wgtabs * phimult
          isoxsmaca_bon(iso) = isoxsmaca_bon(iso) + wgtabs * PhiHom
          if ( isodata%ityp .NE. 3 ) cycle
          isoxsmacf_bon(iso) = isoxsmacf_bon(iso) + wgtfis_bon(ilv, iso) * maclvfis(ilv,iso) * phihom
        enddo
        fl1 = fl1 - phiiso(iso)
      ENDDO
      flavg = 0._8
      DO iso = 1, niso
        id = mapnucl(idiso(iso)); isodata => ldiso(id)
        IF ( .NOT. isodata%lreso ) CYCLE
        flavg = flavg + phiiso(iso)*wphi(iso)
      END DO

      SumSigA = 0._8; invfl1 = 1._8 / fl1
      IF (fl1 .LE. 0._8)  invfl1 = 1._8 / flavg
      DO iso = 1, niso
        id = mapnucl(idiso(iso)); isodata => ldiso(id)
        IF ( .NOT. isodata%lreso ) CYCLE
        IF ( pnum(iso).le.0._8 ) CYCLE
        isoxsmaca_bon(iso) = isoxsmaca_bon(iso) * invfl1 * ratind(iso)
        IF(isodata%ityp .EQ. 3) isoxsmacf_bon(iso) = isoxsmacf_bon(iso) * invfl1 * ratind(iso)
        SumSigA = SumSigA + isoxsmaca_bon(iso)
      ENDDO

      IF(ibon .NE. 1) THEN
        IF(abs(1._8-fl1/fl2) .LE. epsm10) EXIT
      ENDIF

      fl2 = fl1
    ENDDO ! DO ibon = 1, nbonmax

  DO iso = 1, niso
    IsoXsMacA(iso, ig) = isoxsmaca_bon(iso)
    IsoXsMacF(iso, ig) = isoxsmacf_bon(iso)
    IsoXsMacNF(iso, ig) = IsoXsMacF(iso, ig) * IsoMacNu(iso, ig)
    IsoXsMacKF(iso, ig) = IsoXsMacF(iso, ig) * IsoMacK(iso, ig)
  ENDDO

  IF (nTracerCNTL%lRST) THEN
      DO iso = 1, niso
        id = mapnucl(idiso(iso)); isodata => ldiso(id)
        IF ( .NOT. isodata%lreso ) CYCLE
        logxsa = IsoXsMacA(iso, ig)/pnum(iso)
        if (logxsa.le.0._8) CYCLE
        logxsa = dlog(logxsa)
        CALL calcWgt( Tempsq, isodata%rtempsq, isodata%nrtemp, wgttemp, idxtemp, 'A' )
        DO it = 1, 2
          var(1:isodata%nsig0) = isodata%xsalog(1:isodata%nsig0,ig,idxtemp(it))
          CALL calcWgt2( logxsa, var(1:isodata%nsig0), isodata%nsig0, wgtsig(:,it), idxsig(:,it), 'D' )
          s(it) = ( isodata%xss(idxsig(1,it), ig, idxtemp(it)) ** wgtsig(1, it) ) * ( isodata%xss(idxsig(2,it), ig, idxtemp(it)) ** wgtsig(2, it) )
        ENDDO
        xss = s(1) * wgttemp(1) + s(2) * wgttemp(2)
#ifdef __PGI
        if (ieee_is_nan(xss)) cycle
#else
        if (isnan(xss)) cycle
#endif
        xsss = xss * isodata%ssf(ig)
        IsoXsMacS(iso, ig) = xss * pnum(iso)
        IsoXsMacS1(iso, ig) = IsoXsMacS(iso, ig) * isodata%mu
        IsoXsMacSS(iso, ig) = xsss * pnum(iso)
      ENDDO
  ENDIF

ENDIF

xsmaca(ig) = 0.; xsmacs(ig)  = 0.; xsmacf(ig) = 0.; xsmacNf(ig) = 0.; xsmackf(ig) = 0.;
DO iso = 1, niso
  xsmaca(ig) = xsmaca(ig) + IsoXsMacA(iso, ig)
  xsmacs(ig) = xsmacs(ig) + IsoXsMacS(iso, ig)
  xsmacf(ig) = xsmacf(ig) + IsoXsMacF(iso, ig)
  xsmacnf(ig) = xsmacnf(ig) + IsoXsMacNF(iso, ig)
  xsmackf(ig) = xsmackf(ig) + IsoXsMacKF(iso, ig)
ENDDO
xsmacstr(ig) = xsmacs(ig)
DO iso = 1, niso
  xsmacstr(ig) = xsmacstr(ig) - IsoXsMacS1(iso, ig)
ENDDO

IF(lIsoXsOut) THEN
  XsMac%IsoXsMacA(1:niso, ig) = IsoXsMacA(1:niso, ig)
  XsMac%IsoXsMacS0(1:niso, ig)= IsoXsMacS(1:niso, ig)
  XsMac%IsoXsMacS1(1:niso, ig)= IsoXsMacS1(1:niso, ig)
  XsMac%IsoXsMacF(1:niso, ig) = IsoXsMacF(1:niso, ig)
  XsMac%IsoXsMacNF(1:niso, ig)= IsoXsMacNF(1:niso, ig)
  XsMac%IsoXsMackf(1:niso, ig) = IsoXsMacKF(1:niso, ig)
  XsMac%IsoXsMacSS(1:niso, ig)= IsoXsMacSS(1:niso, ig)
ENDIF

NULLIFY(xsmaca, xsmacf, xsmacs, xsmacstr, xsmacnf, xsmackf)
NULLIFY(isodata, jsodata, repisodata)
NULLIFY(pnum, idiso, mlgmaclv, xseq)

END SUBROUTINE

SUBROUTINE GetMacChi_Fxr(Fxr, phi, ig1, ig2, nchi, ng, lchip)
IMPLICIT NONE
TYPE(FxrInfo_Type) :: Fxr
INTEGER :: ig1, ig2, nchi, ng
REAL :: phi(ng)
LOGICAL :: lchip
CALL GetMacChi_Gen(Fxr%chi, Fxr%lres, Fxr%fresoFIso, Fxr%Temp, Fxr%niso, Fxr%idiso, Fxr%pnum, phi, ig1, ig2, nchi, ng, lchip)
END SUBROUTINE

SUBROUTINE GetMacChi_Gen(chi, lres, fResoFIso, Temp, niso, idiso, pnum, phi, ig1, ig2, nchi, ng, lchip)
USE CNTL,   ONLY : nTracerCntl
USE Core_mod, ONLY : GroupInfo
IMPLICIT NONE
INTEGER :: niso, ig1, ig2, nchi, ng
INTEGER :: idiso(niso)
REAL :: pnum(niso)
REAL(4),POINTER :: fResoFIso(:,:)
REAL :: Temp
REAL :: CHI(nchi), phi(ng)
LOGICAL :: lres
LOGICAL :: lchip
TYPE(LIBDATA), POINTER :: isodata
INTEGER :: id, iso, ig
REAL :: PSI, localPsi
!Temp
INTEGER  :: it1, it2
REAL :: wt1, wt2

REAL :: IsoXsMacNf(niso, ng)
!LOGICAL :: lresChi, lres
INTEGER :: iResGrpBeg, iResGrpEnd
iResGrpBeg = GroupInfo%nofg + 1;
iResGrpEnd = GroupInfo%nofg + GroupInfo%norg
!lResChi=.TRUE.
!lResChi=.FALSE.
!lreschi=nTracerCntl%lreschi
!IF(.NOT. lAlloc) CALL AllocIsotopicMacXs()
IsoXsMacNf(1:niso, 1:ng) = 0.
!Initialize the chi
chi(1:nchi) = 0.
PSI = 0
DO iso = 1, niso
  id = mapnucl(idiso(iso))
  isodata => ldiso(id)
  IF(isodata%ityp .le. 1) CYCLE
  IF(isodata%chi(1) .le. epsm10) CYCLE
  !Temp Interpolation variables
  CALL XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2)
  !Obtain a Fission Source
  localpsi = 0._8
  DO ig = 1, ng
    IsoXsMacNf(iso, ig) = (wt1*isodata%signf(ig, it1) + wt2*isodata%signf(ig, it2))*pnum(iso)
    IF(lres.and.nTracerCntl%lrestrmt)THEN
        IF( ig .GE. iResGrpBeg .AND. ig .LE. iResGrpEnd )THEN
            IF(fResoFIso(iso, ig) .NE.1 )THEN
                IsoXsMacNf(iso, ig) = IsoXsMacNf(iso, ig) * fResoFIso(iso, ig)
            ENDIF
        ENDIF
    ENDIF
    localpsi = localpsi + IsoXsMacNf(iso, ig) * phi(ig)
  ENDDO
  psi = psi + localpsi
  !Calculate a Fission Neutron Dist.
  DO ig = ig1, nchihel
    IF(lchip) THEN
      chi(ig) = chi(ig) + localpsi * isodata%chip(ig)
    ELSE
      chi(ig) = chi(ig) + localpsi * isodata%chi(ig)
    END IF
  ENDDO
ENDDO

DO ig = ig1, min(nchihel, ig2)
  chi(ig) = chi(ig) / psi
ENDDO
IF(ASSOCIATED(isodata)) NULLIFY(isodata)

END SUBROUTINE

SUBROUTINE GetMacChid(Fxr, phi, chid, ig1, ig2, niso, ng)
USE CNTL,   ONLY : nTracerCntl
USE Core_mod, ONLY : GroupInfo
IMPLICIT NONE
TYPE(FxrInfo_Type) :: Fxr
REAL :: phi(ng), chid(ng)
INTEGER :: niso, ig1, ig2, ng

TYPE(LIBDATA), POINTER :: isodata
REAL(4),POINTER :: fResoFIso(:,:)
REAL :: pnum(niso)
REAL :: IsoXsMacNf(niso, ng)
REAL :: PSI, localPsi
REAL :: Temp
REAL :: wt1, wt2
INTEGER :: idiso(niso)
INTEGER :: iResGrpBeg, iResGrpEnd
INTEGER :: id, iso, ig
INTEGER  :: it1, it2
LOGICAL :: lres

fResoFIso => Fxr%FresoFIso
Temp = Fxr%Temp
pnum = Fxr%pnum(1:niso)
idiso = Fxr%idiso(1:niso)
lres = Fxr%lres

iResGrpBeg = GroupInfo%nofg + 1;
iResGrpEnd = GroupInfo%nofg + GroupInfo%norg
CALL CP_CA(IsoXsMacNf(1:niso, 1:ng), 0._8, niso, ng)
!Initialize the chi
CALL CP_CA(chid, 0._8, ng)
PSI = 0
DO iso = 1, niso
  id = mapnucl(idiso(iso))
  isodata => ldiso(id)
  IF(isodata%ityp .le. 1) CYCLE
  !IF(isodata%ichi .le. 1) CYCLE
  IF(isodata%chi(1) .le. epsm10) CYCLE
  !Temp Interpolation variables
  CALL XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2)
  !Obtain a Fission Source
  localpsi = 0._8
  DO ig = 1, ng
    IsoXsMacNf(iso, ig) = (wt1*isodata%signf(ig, it1) + wt2*isodata%signf(ig, it2))*pnum(iso)
    IF(lres.and.nTracerCntl%lrestrmt)THEN
        IF( ig .GE. iResGrpBeg .AND. ig .LE. iResGrpEnd )THEN
            IF(fResoFIso(iso, ig) .NE.1 )THEN
                IsoXsMacNf(iso, ig) = IsoXsMacNf(iso, ig) * fResoFIso(iso, ig)
            ENDIF
        ENDIF
    ENDIF
    localpsi = localpsi + IsoXsMacNf(iso, ig) * phi(ig)
  ENDDO
  localpsi = localpsi * isodata%beta(0)
  psi = psi + localpsi
  !Calculate a Fission Neutron Dist.
  DO ig = 1, ng
    chid(ig) = chid(ig) + localpsi * isodata%chid(ig)
  ENDDO
ENDDO

DO ig = 1, ng
  chid(ig) = chid(ig) / psi
ENDDO
IF(ASSOCIATED(isodata)) NULLIFY(isodata)
NULLIFY(fResoFIso)

END SUBROUTINE

SUBROUTINE GetMacChidk(Fxr, phi, chidk, ig1, ig2, niso, ng, nprec)
USE CNTL,   ONLY : nTracerCntl
USE Core_mod, ONLY : GroupInfo
IMPLICIT NONE
TYPE(FxrInfo_Type) :: Fxr
REAL :: phi(ng), chidk(ng, nprec)
INTEGER :: niso, ig1, ig2, ng, nprec

TYPE(LIBDATA), POINTER :: isodata
REAL(4),POINTER :: fResoFIso(:,:)
REAL :: pnum(niso)
REAL :: IsoXsMacNf(niso, ng)
REAL :: PSI(nprec), localPsi, localpsik(nprec)
REAL :: Temp
REAL :: wt1, wt2
INTEGER :: idiso(niso)
INTEGER :: iResGrpBeg, iResGrpEnd
INTEGER :: id, iso, ig, iprec
INTEGER  :: it1, it2
LOGICAL :: lres

fResoFIso => Fxr%FresoFIso
Temp = Fxr%Temp
pnum = Fxr%pnum(1:niso)
idiso = Fxr%idiso(1:niso)
lres = Fxr%lres

iResGrpBeg = GroupInfo%nofg + 1;
iResGrpEnd = GroupInfo%nofg + GroupInfo%norg
CALL CP_CA(IsoXsMacNf(1:niso, 1:ng), 0._8, niso, ng)
!Initialize the chi
CALL CP_CA(chidk, 0._8, ng, nprec)
PSI = 0
DO iso = 1, niso
  id = mapnucl(idiso(iso))
  isodata => ldiso(id)
  IF(isodata%ityp .le. 1) CYCLE
  !IF(isodata%ichi .le. 1) CYCLE
  IF(isodata%chi(1) .le. epsm10) CYCLE
  !Temp Interpolation variables
  CALL XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2)
  !Obtain a Fission Source
  localpsi = 0._8
  DO ig = 1, ng
    IsoXsMacNf(iso, ig) = (wt1*isodata%signf(ig, it1) + wt2*isodata%signf(ig, it2))*pnum(iso)
    IF(lres.and.nTracerCntl%lrestrmt)THEN
        IF( ig .GE. iResGrpBeg .AND. ig .LE. iResGrpEnd )THEN
            IF(fResoFIso(iso, ig) .NE.1 )THEN
                IsoXsMacNf(iso, ig) = IsoXsMacNf(iso, ig) * fResoFIso(iso, ig)
            ENDIF
        ENDIF
    ENDIF
    localpsi = localpsi + IsoXsMacNf(iso, ig) * phi(ig)
  ENDDO
  DO iprec = 1, nprec
    localpsik(iprec) = localpsi * isodata%beta(iprec)
    psi(iprec) = psi(iprec) + localpsik(iprec)
    !Calculate a Fission Neutron Dist.
    DO ig = 1, ng
      chidk(ig, iprec) = chidk(ig, iprec) + localpsik(iprec) * isodata%chidk(ig, iprec)
    ENDDO
  END DO
ENDDO

DO iprec = 1, nprec
  DO ig = 1, ng
    chidk(ig, iprec) = chidk(ig, iprec) / psi(iprec)
  ENDDO
END DO
IF(ASSOCIATED(isodata)) NULLIFY(isodata)
NULLIFY(fResoFIso)

END SUBROUTINE

SUBROUTINE MacXsKf_CrCsp(XsMac, Fxr, ig1, ig2, ng, eigv, kin, lCrCspFtn)
USE XsUtil_mod,   ONLY : GetXsMacDat,   ReturnXsMacDat
USE CrCsp_MOD,    ONLY : MacXsKfCsp
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
TYPE(Fxrinfo_type) :: Fxr
REAL :: eigv
INTEGER :: ig1, ig2, ng
LOGICAL :: kin, lCrCspFtn

TYPE(XsMac_TYPE), POINTER :: XsMac1, XsMac2

REAL :: temp
INTEGER :: niso
INTEGER, POINTER :: idiso(:)
REAL, POINTER ::pnum(:)

IF(lCrCspFtn .AND. Fxr%lCrCspFtn) THEN
  CALL GetXsMacDat(XsMac1, ng, .FALSE.)
  Call GetXsMacDat(XsMac2, ng, .FALSE.)
  niso = Fxr%CspFxr%niso(1); Temp = Fxr%Temp
  Pnum => Fxr%CspFxr%pnum(:, 1);   idiso => Fxr%CspFxr%isolist(:,1)
  CALL MacXskf_gen(XsMac1, niso, temp, idiso, pnum, ig1, ig2, ng, eigv, kin)

  niso = Fxr%CspFxr%niso(2); Temp = Fxr%Temp
  Pnum => Fxr%CspFxr%pnum(:, 2);   idiso => Fxr%CspFxr%isolist(:,2)
  CALL MacXskf_gen(XsMac2, niso, temp, idiso, pnum, ig1, ig2, ng, eigv, kin)

  CALL MacXskfCsp(XsMac, XsMac1, XsMac2, Fxr%CspFxr, ig1, ig2)
  NULLIFY(pnum); NULLIFY(idiso)
  CALL ReturnXsMacDat(XsMac1); CALL ReturnXsMacDat(XsMac2)
ELSE
  niso = Fxr%niso; temp = Fxr%Temp
  pnum => Fxr%pnum; idiso => Fxr%idiso
  CALL MacXsKf_gen(XsMac, niso, temp, idiso, pnum, ig1, ig2, ng, eigv, kin)
  NULLIFY(pnum); NULLIFY(idiso)
ENDIF
END SUBROUTINE


SUBROUTINE MacXsKf_Fxr(XsMac, Fxr, ig1, ig2, ng, eigv, kin)
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
TYPE(Fxrinfo_type) :: Fxr
REAL :: eigv
INTEGER :: ig1, ig2, ng
LOGICAL :: kin
REAL :: temp
INTEGER :: niso
INTEGER, POINTER :: idiso(:)
REAL, POINTER ::pnum(:)

niso = Fxr%niso; temp = Fxr%Temp
pnum => Fxr%pnum; idiso => Fxr%idiso

CALL MacXsKf_gen(XsMac, niso, temp, idiso, pnum, ig1, ig2, ng, eigv, kin)
NULLIFY(pnum); NULLIFY(idiso)
END SUBROUTINE

SUBROUTINE MacXskf_Gen(XsMac, niso, temp, idiso, pnum, ig1, ig2, ng, eigv, kin)
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
REAL :: temp, eigv
INTEGER :: niso
INTEGER :: idiso(niso)
REAL :: pnum(niso)
INTEGER :: ig1, ig2, ng
LOGICAL :: kin

REAL, POINTER :: xsmackf(:)
TYPE(LIBDATA), POINTER :: isodata
INTEGER :: iso, id, ig
REAL :: reigv

INTEGER :: it1, it2
REAL :: wt1, wt2

IF(.NOT. XsMac%lalloc) THEN
  XsMac%ng = ng
  CALL AllocXsMac(XsMac)
ENDIF

XsMacKf => XsMac%XsMacKf
xsmackf(ig1:ig2) = 0.

reigv = 1._8 / eigv

DO iso = 1, niso
  id = MapNucl(idiso(iso)); isodata => ldiso(id)
  IF(mapfis(idiso(iso)) .EQ. 0) CYCLE
  CALL XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2)
  DO ig = ig1, ig2
    XsMacKf(ig) = XsMacKf(ig) + pnum(iso) * (wt1 * isodata%sigf(ig, it1) + wt2 * isodata%sigf(ig, it2))*isodata%kappa
  ENDDO
ENDDO

!IF(kin) CALL  MULTI_CA(reigv, XsMacNf(ig1 : ig2), ig2 - ig1 + 1)
NULLIFY(XsMacKf, isodata)
END SUBROUTINE

SUBROUTINE SetCriticalNumberDensity(Temp)
! Criteria: SigTotal_Mactro < 1.E7
TYPE(LIBDATA), POINTER :: isodata
REAL :: Temp, wt1, wt2, xst_mic
INTEGER :: id, ig, it1, it2

ig = noghel
DO id = 1, nelthel
  isodata => ldiso(id)
  CALL XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2)
  IsoData%Crit_nd = epsm10
  xst_mic = wt1 * isodata%siga(ig, it1) + wt2 * isodata%siga(ig, it2)
  xst_mic = xst_mic + wt1 * isodata%sigs(ig, it1) + wt2 * isodata%sigs(ig, it2)
  IF(abs(xst_mic) .GT. epsm10) IsoData%Crit_nd = epsm10 / xst_mic
ENDDO
END SUBROUTINE

! Subroutine to generate effective XS with PSM
SUBROUTINE PSMEffXSGen(Core, Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
!Effective Xs Generation
USE PARAM
USE PointXSRT_MOD
USE TYPEDEF,        ONLY : CoreInfo_Type,       FxrInfo_Type,          GroupInfo_Type,     &
                           PE_Type,             THInfo_Type,           pin_type,           &
                           ResVarpin_type,      cell_type,             XsMac_Type
USE CNTL,           ONLY : nTracerCntl_Type
!USE MacXsLib_mod,   ONLY : MacXsBase
USE BasicOperation, ONLY : CP_VA
USE IOUTIL,         ONLY : message
USE FILES,          ONLY : io8
USE TH_Mod,         ONLY : GetPinFuelTemp
USE MPIComm_Mod,    ONLY : MPI_SYNC
USE OMP_LIB
USE TIMER,            ONLY : nTracer_dclock, TimeChk
!
IMPLICIT NONE
!
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(THInfo_Type) :: THInfo
REAL :: eigv
TYPE(GroupInfo_TYPE) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
!
TYPE(Pin_Type), POINTER :: PIN(:)
TYPE(ResVarPin_Type), POINTER :: ResVarPin(:,:),mypin
TYPE(Cell_Type), POINTER :: CellInfo(:), myCell
TYPE(FxrInfo_Type), POINTER :: myFxr
TYPE(XsMac_Type), SAVE :: XsMac(nThreadMax)
TYPE(PSMMacXs_Type), SAVE :: PSMXsMac(nThreadMax)
!
INTEGER :: nxy, myzb, myze, iz, ipin, xyb, xye, FxrIdxSt, nlocalFxr, j, ifxr, icel  !--- (xyb, xye)CNJ Edit : Domain Decomposition + MPI
INTEGER :: ng, ig, iresoGrp1, iresoGrp2, niso, iso, tid
REAL :: PinFuelTempAvgsq, XsMacAold, XsMacFold, XsMacNFold, XsMackFold, XsMacSold, XsMacStrold,ndsumcat(3),volsum
REAL :: isoXsMacAold(500),isoXsMacFold(500),isoXsMacNFold(500),isoXsMackFold(500),isoXsMacSold(500),isoXsMacS1old(500),isoXsMacSSold(500),isoxsmaccapold(500)
LOGICAL :: lxslib, lAIC, lmcxs
INTEGER :: nfueldiv, idxfuelsrt, iPSMcel, id, i, ifuelring, idpxs
TYPE(LIBDATA), POINTER :: isodata
REAL :: TimeEnd, TimeBeg
INTEGER, POINTER :: idiso(:)
!
IF (.not.nTracerCntl%lrestrmt) RETURN
lxslib  = nTracerCntl%lxslib
IF(.NOT. lxslib) RETURN
!
WRITE(mesg, '(A)') 'Update PSM Effective XSs...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)
TimeBeg = nTracer_dclock(FALSE,FALSE)
Pin => Core%Pin; CellInfo => Core%CellInfo
ResVarPin => Core%ResVarPin;
nxy = Core%nxy;  myzb = PE%myzb; myze = PE%myze

!--- CNJ Edit : Domain Decomposition + MPI
xyb = PE%myPinBeg; xye = PE%myPinEnd
IF (PE%RTMASTER) THEN
  xyb = 1; xye = nxy
ENDIF

ng = GroupInfo%ng
iresoGrp1 = GroupInfo%nofg + 1; iresoGrp2 = GroupInfo%nofg + GroupInfo%norg

CALL UpdtCoreIsoInfoPSM(Core, Fxr, PE)
!
CALL omp_set_num_threads(PE%nThread)
tid = 1
DO iz = myzb, myze
  IF (.NOT. Core%lFuelPlane(iz)) CYCLE
  !$OMP PARALLEL DEFAULT(SHARED)      &
  !$OMP PRIVATE(tid, ipin, mypin, FxrIdxSt, icel, myCell, iPSMcel, nLocalFxr, PinFuelTempAvgsq, ndsumcat, volsum, j, ifxr, myfxr, niso, &
  !$OMP         iso, id, isodata, nfueldiv, idxfuelsrt, ig, idiso, ifuelring, idpxs, &
  !$OMP         isoXsMacAold,isoXsMacFold,isoXsMacNFold,isoXsMackFold,isoXsMacSold,isoXsMacS1old,isoXsMacSSold,isoxsmaccapold)
  !$ tid = omp_get_thread_num()+1
  !$OMP DO
  DO ipin = xyb, xye
    mypin => ResVarPin(ipin,iz)
    FxrIdxSt = Pin(ipin)%FxrIdxSt
    icel = Pin(ipin)%Cell(iz)
    myCell => CellInfo(icel)
    IF (.NOT.myCell%lPSMcel) CYCLE
    iPSMcel = myCell%icelPSM
    nlocalFxr = CellInfo(icel)%nFxr
    PinFuelTempAvgsq = dsqrt(GetPinFuelTemp(Core, Fxr, iz, ipin))
    ! Set the Potential, total XS o each FXR  in case the moderator condition changes -> after debug sink to TH module
    mypin%siglpM = 0._8
    mypin%siglpMcat = 0._8
    mypin%mavgm = 0._8
    ndsumcat = 0._8
    volsum = 0._8
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      myFxr => Fxr(ifxr, iz)
      IF( myFxr%lfuel ) EXIT
      niso = myFxr%niso
      volsum = volsum + myCell%fxrvol(nLocalFxr-j+1)
      DO iso = 1, niso
        id = mapnucl(myFxr%idiso(iso));
        isodata => ldiso(id)
        ndsumcat(isodata%icatm) = ndsumcat(isodata%icatm) + myFXR%pnum(iso) * myCell%fxrvol(nLocalFxr-j+1)
        mypin%mavgm(isodata%icatm) = mypin%mavgm(isodata%icatm) + myFXR%pnum(iso) * isodata%aw * myCell%fxrvol(nLocalFxr-j+1)
        IF (isodata%sigp.EQ. 0._8) THEN
          mypin%siglpM = mypin%siglpM + myFXR%pnum(iso) * isodata%lamsigp1G * myCell%fxrvol(nLocalFxr-j+1)
          mypin%siglpMcat(isodata%icatm) = mypin%siglpMcat(isodata%icatm) + myFXR%pnum(iso) * isodata%lamsigp1G * myCell%fxrvol(nLocalFxr-j+1)
        ELSE
          mypin%siglpM = mypin%siglpM + myFXR%pnum(iso) * isodata%sigp * myCell%fxrvol(nLocalFxr-j+1)
          mypin%siglpMcat(isodata%icatm) = mypin%siglpMcat(isodata%icatm) + myFXR%pnum(iso) * isodata%sigp * myCell%fxrvol(nLocalFxr-j+1)
        END IF
      END DO
    END DO
    mypin%siglpM = mypin%siglpM / volsum
    myPin%lmsrc = .TRUE.
    DO i = 1, 3
      IF (ndsumcat(i) .LT. 1E-20) THEN
        myPin%lmsrc(i) = .FALSE.
        CYCLE
      END IF
      mypin%siglpMcat(i) = mypin%siglpMcat(i) / volsum
      mypin%mavgm(i) = mypin%mavgm(i) / ndsumcat(i)
      mypin%alpham(i) = ((mypin%mavgm(i)-1._8)/(mypin%mavgm(i)+1._8))**2
      mypin%invsubalpham(i) = 1._8/(1._8-mypin%alpham(i))
      mypin%sigpm_inval(i) = mypin%siglpMcat(i)*mypin%invsubalpham(i)
    END DO
    !
    nfueldiv = myCell%nfueldiv
    idxfuelsrt = FxrIdxSt+myCell%nnonfuel
    IF (.NOT.myPin%lfuel .OR. myCell%lrect) CYCLE
    CALL EffMacXsPSM(PSMXsMac(tid), Fxr(idxfuelsrt:FxrIdxSt+nLocalFXR-1, iz), mypin, PinFuelTempAvgsq, nfueldiv, iresoGrp1, iresoGrp2, ng, .TRUE.,iPSMcel)

    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      myFxr => Fxr(ifxr, iz)
      IF (idxfuelsrt.GT.ifxr) CYCLE ! beyond the fuel ring !.OR.ifxr.GT.FxrIdxSt+nLocalFXR-1
      IF( .NOT. myFxr%lfuel ) CYCLE
      niso = myFxr%niso
      myfxr%idiso_pastpsm(1:niso) = myfxr%idiso(1:niso)
      do ig = iResoGrp1, iResoGrp2
        do iso = 1, niso
          myFXR%fresoAIso(iso,ig) = 1._8
          myFXR%fresoFIso(iso,ig) = 1._8
        enddo
      enddo
      IF (nTRACERCntl%lgamma) THEN
        do ig = iResoGrp1, iResoGrp2
            do iso = 1, niso
              myFxr%fresocapIso(iso, ig) = 1._8
            enddo
        enddo
      END IF
      IF (nTracerCntl%lRST) THEN
        do ig = iResoGrp1, iResoGrp2
          do iso = 1, niso
            myFXR%fresoSIso(iso,ig) = 1._8
            myFXR%fresoSSIso(iso,ig) = 1._8
            myFXR%fresoS1Iso(iso,ig) = 1._8
          enddo
        enddo
      ENDIF
      idiso =>  myFXR%idiso
      ifuelring = nfueldiv - (ifxr - idxfuelsrt)
      CALL MacXsBase(XsMac(tid), myFxr, iResoGrp1, iResoGrp2, ng, eigv, FALSE, TRUE)
      !nuclides without pxs information (fresoXiso = 1._8)
      DO iso = 1, niso
        id = mapnucl(idiso(iso));    isodata => ldiso(id)
        idpxs = isodata%pxsid
        IF (idpxs .NE. 0) CYCLE
        !IF (idpxs .EQ. 0 .OR. idiso(iso) .LT. 90000) THEN
        DO ig = iResoGrp1, iResoGrp2
          PSMXsMac(tid)%IsoMGXsMacA(iso,ig,ifuelring) = XsMac(tid)%IsoXsMacA(iso,ig)
          IF (isoXsMacFold(iso) .GT. epsm8) THEN
            PSMXsMac(tid)%IsoMGXsMacF(iso,ig,ifuelring) = XsMac(tid)%IsoXsMacF(iso,ig)
          END IF
        END DO ! IG
        !END IF
      END DO ! ISO
      ! Calculate freso
      DO ig = iResoGrp1, iResoGrp2
        CALL CP_VA(isoXsMacFold(1:niso),XsMac(tid)%IsoXsMacF(1:niso,ig),niso)
        CALL CP_VA(isoXsMacNFold(1:niso),XsMac(tid)%IsoXsMacNF(1:niso,ig),niso)
        CALL CP_VA(isoXsMackFold(1:niso),XsMac(tid)%IsoXsMackF(1:niso,ig),niso)
        CALL CP_VA(isoXsMacAold(1:niso),XsMac(tid)%IsoXsMacA(1:niso,ig),niso)
        IF (nTracerCntl%lRST) THEN
          CALL CP_VA(isoXsMacSold(1:niso),XsMac(tid)%IsoXsMacS0(1:niso,ig),niso)
          CALL CP_VA(isoXsMacS1old(1:niso),XsMac(tid)%IsoXsMacS1(1:niso,ig),niso)
          CALL CP_VA(isoXsMacSSold(1:niso),XsMac(tid)%IsoXsMacSS(1:niso,ig),niso)
        ENDIF
        CALL CP_VA(IsoXsMacCapOld(1:niso),XsMac(tid)%IsoXsRadCap(1:niso,ig),niso)     !-- JSU EDIT 20170727
        !
        DO iso = 1, niso
          myFXR%fresoAIso(iso,ig) = PSMXsMac(tid)%IsoMGXsMacA(iso,ig,ifuelring) / isoXsMacAold(iso)
          IF (isoXsMacFold(iso) .GT. epsm8)   myFXR%fresoFIso(iso,ig) = PSMXsMac(tid)%IsoMGXsMacF(iso,ig,ifuelring) / isoXsMacFold(iso)
        END DO

        IF (nTRACERCntl%lGamma) THEN
          IF(isoXsMacFold(iso) .gt. epsm8) THEN ! Radioactive Capture resonance fraction...
            DO iso = 1, niso
              myFxr%fresocapIso(iso, ig) = (PSMXsMac(tid)%IsoMGXsMacA(iso,ig,ifuelring) - PSMXsMac(tid)%IsoMGXsMacF(iso,ig,ifuelring)) / (isoXsMacAold(iso) - isoXsMacFold(iso))
            END DO
          ELSE
            DO iso = 1, niso
              myFxr%fresocapIso(iso, ig) = PSMXsMac(tid)%IsoMGXsMacA(iso,ig,ifuelring) / isoXsMacAold(iso)
            END DO
          END IF
        END IF

        IF (nTracerCntl%lRST) THEN
          DO iso = 1, niso
            if (abs(isoXsMacSold(iso)).GT. epsm8) myFXR%fresoSIso(iso,ig) = PSMXsMac(tid)%IsoMGXsMacS(iso,ig,ifuelring)/isoXsMacSold(iso)
          ENDDO
        ENDIF

      ENDDO ! ig
    ENDDO ! fxr
  ENDDO ! ipin
!  STOP 'DEBUGGING ** SUBGROUPEFFXS LINE 231'
  !$OMP END DO
  !$OMP END PARALLEL
ENDDO ! iz
TimeEnd = nTracer_dclock(FALSE,FALSE)
TimeChk%PSMTime = TimeChk%PSMTime + TimeEnd - TimeBeg
CALL MPI_SYNC(PE%MPI_NTRACER_COMM)

NULLIFY(PIN, ResVarPin, CellINfo, myFxr, myPin, myCell)
END SUBROUTINE


SUBROUTINE GetfresoFXR(Core, Fxr, THInfo, eigv, GroupInfo, nTracerCntl, PE)
! Subroutine to get freso of each FXR's macroscopic XSs based on fresoISO
USE PARAM
USE PointXSRT_MOD
USE TYPEDEF,        ONLY : CoreInfo_Type,       FxrInfo_Type,          GroupInfo_Type,     &
                           PE_Type,             THInfo_Type,           pin_type,           &
                           ResVarpin_type,      cell_type,             XsMac_Type
USE CNTL,           ONLY : nTracerCntl_Type
!USE MacXsLib_mod,   ONLY : MacXsBase
USE BasicOperation, ONLY : CP_VA
USE IOUTIL,         ONLY : message
USE FILES,          ONLY : io8
USE TH_Mod,         ONLY : GetPinFuelTemp
USE XSLIB_MOD,      ONLY : libdata, ldiso, mapnucl
USE MPIComm_Mod,    ONLY : MPI_SYNC
USE OMP_LIB
USE TIMER,            ONLY : nTracer_dclock, TimeChk
!
IMPLICIT NONE
!
TYPE(CoreInfo_Type) :: Core
TYPE(FxrInfo_Type), POINTER :: Fxr(:, :)
TYPE(THInfo_Type) :: THInfo
REAL :: eigv
TYPE(GroupInfo_TYPE) :: GroupInfo
TYPE(nTracerCntl_Type) :: nTracerCntl
TYPE(PE_Type) :: PE
!
TYPE(Pin_Type), POINTER :: PIN(:)
TYPE(ResVarPin_Type), POINTER :: ResVarPin(:,:),mypin
TYPE(Cell_Type), POINTER :: CellInfo(:), myCell
TYPE(FxrInfo_Type), POINTER :: myFxr
TYPE(XsMac_Type), SAVE :: XsMac(nThreadMax)
TYPE(PSMMacXs_Type), SAVE :: PSMXsMac(nThreadMax)
!
INTEGER :: nxy, myzb, myze, iz, ipin, xyb, xye, FxrIdxSt, nlocalFxr, j, ifxr, icel  !--- (xyb, xye)CNJ Edit : Domain Decomposition + MPI
INTEGER :: ng, ig, iresoGrp1, iresoGrp2, niso, iso, tid
REAL :: PinFuelTempAvgsq, XsMacAold, XsMacFold, XsMacNFold, XsMackFold, XsMacSold, XsMacStrold
REAL :: XsMacAnew, XsMacFnew, XsMacNFnew, XsMackFnew, XsMacSnew,XsMacStrnew
REAL :: isoXsMacAold(500),isoXsMacFold(500),isoXsMacNFold(500),isoXsMackFold(500),isoXsMacSold(500),isoXsMacS1old(500),isoXsMacSSold(500),isoxsmaccapold(500)
LOGICAL :: lxslib, lAIC, lmcxs
INTEGER :: nfueldiv, iPSMcel, id, i, ifuelring, idpxs, idxfuelsrt
TYPE(LIBDATA), POINTER :: isodata
REAL :: TimeEnd, TimeBeg
INTEGER, POINTER :: idiso(:), idiso_psm(:)
INTEGER :: idxiso_psm(1000)
INTEGER :: jso
!
IF (.not.nTracerCntl%lrestrmt) RETURN
lxslib  = nTracerCntl%lxslib
IF(.NOT. lxslib) RETURN
!
!WRITE(mesg, '(A)') 'Update FXR freso by PSM...'
IF(PE%Master) CALL message(io8, TRUE, TRUE, mesg)
TimeBeg = nTracer_dclock(FALSE,FALSE)
Pin => Core%Pin; CellInfo => Core%CellInfo
ResVarPin => Core%ResVarPin;
nxy = Core%nxy;  myzb = PE%myzb; myze = PE%myze

!--- CNJ Edit : Domain Decomposition + MPI
xyb = PE%myPinBeg; xye = PE%myPinEnd
IF (PE%RTMASTER) THEN
  xyb = 1; xye = nxy
ENDIF

ng = GroupInfo%ng
iresoGrp1 = GroupInfo%nofg + 1; iresoGrp2 = GroupInfo%nofg + GroupInfo%norg

CALL UpdtCoreIsoInfoPSM(Core, Fxr, PE)
!
CALL omp_set_num_threads(PE%nThread)
tid = 1
DO iz = myzb, myze
  IF (.NOT. Core%lFuelPlane(iz)) CYCLE
  !$OMP PARALLEL DEFAULT(SHARED)      &
  !$OMP PRIVATE(tid, ipin, mypin, FxrIdxSt, icel, myCell, nLocalFxr, j, ifxr, myfxr, niso, &
  !$OMP         iso, id, isodata, idxfuelsrt, ig, idiso, idpxs, XsMacAold, XsMacFold, XsMacNFold, XsMackFold, XsMacSold, &
  !$OMP         isoXsMacAold,isoXsMacFold,isoXsMacNFold,isoXsMackFold,isoXsMacSold, &
  !$OMP         XsMacAnew, XsMacFnew, XsMacNFnew, XsMackFnew, XsMacSnew, idxiso_psm, idiso_psm, jso)
  !$ tid = omp_get_thread_num()+1
  !$OMP DO
  DO ipin = xyb, xye
    mypin => ResVarPin(ipin,iz)
    FxrIdxSt = Pin(ipin)%FxrIdxSt
    icel = Pin(ipin)%Cell(iz)
    myCell => CellInfo(icel)
    IF (.NOT.myCell%lPSMcel) CYCLE
    nlocalFxr = CellInfo(icel)%nFxr
    idxfuelsrt = FxrIdxSt+myCell%nnonfuel
    DO j = 1, nLocalFxr
      ifxr = FxrIdxSt + j -1
      myFxr => Fxr(ifxr, iz)
      IF( .NOT. myFxr%lfuel ) CYCLE
      IF (idxfuelsrt.GT.ifxr) CYCLE ! beyond the fuel ring !.OR.ifxr.GT.FxrIdxSt+nLocalFXR-1
      niso = myFxr%niso
      ! initialize the macroscopic freso of each FXR in pin
      do ig = iResoGrp1, iResoGrp2
        myFxr%fresoa(ig) = 1._8
        myFxr%fresof(ig) = 1._8
        myFxr%fresoNf(ig) = 1._8
        myFxr%fresokf(ig) = 1._8
      enddo
      IF (nTracerCntl%lRST) THEN
        do ig = iResoGrp1, iResoGrp2
          myFxr%fresos(ig) = 1._8
        enddo
      ENDIF
      idiso =>  myFXR%idiso
      idiso_psm => myFXR%idiso_pastpsm
      DO iso = 1, niso
        DO jso = 1, GroupInfo%ntiso
          IF (idiso_psm(jso).EQ.idiso(iso)) EXIT
        END DO
        idxiso_psm(iso) = jso
      END DO
      CALL MacXsBase(XsMac(tid), myFxr, iResoGrp1, iResoGrp2, ng, eigv, FALSE, TRUE)
      ! Calculate freso
      DO ig = iResoGrp1, iResoGrp2
        XsMacAold = XsMac(tid)%XsMacA(ig); XsMacFold = XsMac(tid)%XsMacF(ig);
        XsMacNFold = XsMac(tid)%XsMacNF(ig); XsMackFold = XsMac(tid)%XsMackF(ig)
        XsMacSold = XsMac(tid)%XsMacS(ig)
        XsMacAnew   = XsMacAold; XsMacFnew   = XsMacFold;
        XsMacNFnew  = XsMacNFold; XsMackFnew  = XsMackFold;
        XsMacSnew   = XsMacSold;
        CALL CP_VA(isoXsMacFold(1:niso),XsMac(tid)%IsoXsMacF(1:niso,ig),niso)
        CALL CP_VA(isoXsMacNFold(1:niso),XsMac(tid)%IsoXsMacNF(1:niso,ig),niso)
        CALL CP_VA(isoXsMackFold(1:niso),XsMac(tid)%IsoXsMackF(1:niso,ig),niso)
        CALL CP_VA(isoXsMacAold(1:niso),XsMac(tid)%IsoXsMacA(1:niso,ig),niso)
        IF (nTracerCntl%lRST) THEN
          CALL CP_VA(isoXsMacSold(1:niso),XsMac(tid)%IsoXsMacS0(1:niso,ig),niso)
        ENDIF
        !
        DO iso = 1, niso
          id = mapnucl(idiso(iso));    isodata => ldiso(id)
          idpxs = isodata%pxsid
          IF (idpxs .EQ. 0) CYCLE
          XsMacAnew   = XsMacAnew   + (myFxr%fresoAISO(idxiso_psm(iso),ig) - 1._8)*isoXsMacAold(iso)
          XsMacFnew   = XsMacFnew   + (myFxr%fresoFISO(idxiso_psm(iso),ig) - 1._8)*isoXsMacFold(iso)
          XsMacNFnew  = XsMacNFnew  + (myFxr%fresoFISO(idxiso_psm(iso),ig) - 1._8)*isoXsMacNFold(iso)
          XsMackFnew  = XsMackFnew  + (myFxr%fresoFISO(idxiso_psm(iso),ig) - 1._8)*isoXsMacKFold(iso)
        END DO
        IF (nTracerCntl%lRST) THEN
          DO iso = 1, niso
            XsMacSnew   = XsMacSnew   + (myFxr%fresoSISO(idxiso_psm(iso),ig) - 1._8)*isoXsMacSold(iso)
          END DO
        END IF
        !
        myFxr%fresoa(ig) = XsMacAnew / XsMacAold
        IF(XsMacFold .gt. epsm8) THEN
          myFxr%fresof(ig) = XsMacFnew / XsMacFold
          myFxr%fresonf(ig) = XsMacNFnew / XsMacNFold
          myFxr%fresokf(ig) = XsMackFnew / XsMacKFold
        END IF
        IF (nTracerCntl%lRST) myFxr%fresos(ig) = XsMacSnew / XsMacSold
      ENDDO ! ig
    ENDDO ! fxr
  ENDDO ! ipin
!  STOP 'DEBUGGING ** SUBGROUPEFFXS LINE 231'
  !$OMP END DO
  !$OMP END PARALLEL
ENDDO ! iz
TimeEnd = nTracer_dclock(FALSE,FALSE)
TimeChk%PSMfresoTime = TimeChk%PSMfresoTime + TimeEnd - TimeBeg
CALL MPI_SYNC(PE%MPI_NTRACER_COMM)
!DEALLOCATE(idxiso_psm)
NULLIFY(PIN, ResVarPin, CellINfo, myFxr, myPin, myCell)
END SUBROUTINE
END MODULE
