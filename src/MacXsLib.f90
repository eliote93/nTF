#include <defines.h>
MODULE MacXsLib_Mod
USE PARAM
USE TYPEDEF,         ONLY : XsMac_Type,     GroupInfo_Type,      Fxrinfo_type
USE XSLIB_MOD
USE XsUtil_mod,      ONLY : AllocXsMac,          AllocMacIsoXs,  &
                            LineIntPol, LineIntPol2,   CalcWgt,   CalcWgt2,   &!TempIntPol_RIA,      
                            XsTempInterpolation, P1XsTempInterpolation,            &
                            P2XsTempInterpolation, P3XsTempInterpolation
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

CALL CP_CA(XsMac%xsmact(igb:ige), 0._8 , ige - igb + 1)
CALL CP_CA(XsMac%xsmaca(igb:ige), 0._8 , ige - igb + 1)
CALL CP_CA(XsMac%xsmacs(igb:ige), 0._8 , ige - igb + 1)
CALL CP_CA(XsMac%xsmactr(igb:ige), 0._8 , ige - igb + 1)
CALL CP_CA(XsMac%xsmacstr(igb:ige), 0._8 , ige - igb + 1)
CALL CP_CA(XsMac%xsmacf(igb:ige), 0._8 , ige - igb + 1)
CALL CP_CA(XsMac%xsmacnf(igb:ige), 0._8 , ige - igb + 1)
CALL CP_CA(XsMac%xsmackf(igb:ige), 0._8 , ige - igb + 1)

CALL CP_CA(IsoXsMacT(1:niso, igb:ige), 0._8, niso, ige - igb + 1)
CALL CP_CA(IsoXsMacA(1:niso, igb:ige), 0._8, niso, ige - igb + 1)
CALL CP_CA(IsoXsMacS0(1:niso, igb:ige), 0._8, niso, ige - igb + 1)
CALL CP_CA(IsoXsMacS1(1:niso, igb:ige), 0._8, niso, ige - igb + 1)
CALL CP_CA(IsoXsMacSS(1:niso, igb:ige), 0._8, niso, ige - igb + 1)
CALL CP_CA(IsoXsMacTr(1:niso, igb:ige), 0._8, niso, ige - igb + 1)
CALL CP_CA(IsoXsMacf(1:niso, igb:ige), 0._8, niso, ige - igb + 1)
CALL CP_CA(IsoXsMacNf(1:niso, igb:ige), 0._8, niso, ige - igb + 1)
CALL CP_CA(IsoXsMackf(1:niso, igb:ige), 0._8, niso, ige - igb + 1)

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
  CALL CP_VA(XsMac%IsoXsMacT(1:niso, igb:ige), IsoXsMacT(1:niso, igb:ige), niso, ige-igb+1)
  CALL CP_VA(XsMac%IsoXsMacA(1:niso, igb:ige), IsoXsMacA(1:niso, igb:ige), niso, ige-igb+1)
  CALL CP_VA(XsMac%IsoXsMacS0(1:niso, igb:ige), IsoXsMacS0(1:niso, igb:ige), niso, ige-igb+1)
  CALL CP_VA(XsMac%IsoXsMacS1(1:niso, igb:ige), IsoXsMacS1(1:niso, igb:ige), niso, ige-igb+1)
  CALL CP_VA(XsMac%IsoXsMacSS(1:niso, igb:ige), IsoXsMacSS(1:niso, igb:ige), niso, ige-igb+1)
  CALL CP_VA(XsMac%IsoXsMactr(1:niso, igb:ige), IsoXsMactr(1:niso, igb:ige), niso, ige-igb+1)
  CALL CP_VA(XsMac%IsoXsMacf(1:niso, igb:ige), IsoXsMacf(1:niso, igb:ige), niso, ige-igb+1)
  CALL CP_VA(XsMac%IsoXsMacNf(1:niso, igb:ige), IsoXsMacNf(1:niso, igb:ige), niso, ige-igb+1)
  CALL CP_VA(XsMac%IsoXsMacKf(1:niso, igb:ige), IsoXsMacKf(1:niso, igb:ige), niso, ige-igb+1)
  CALL CP_VA(XsMac%IsoXsRadCap(1:niso, igb:ige), IsoXsRadCap(1:niso, igb:ige), niso, ige-igb+1) !-- JSU EDIT 20170727
ENDIF

IF(kin) CALL MULTI_CA(reigv, xsmacnf(igb:ige), ige - igb + 1)

NULLIFY(xsmact);  NULLIFY(xsmaca);   NULLIFY(xsmacs)
NULLIFY(xsmactr); NULLIFY(xsmacstr); 
NULLIFY(xsmacf);  NULLIFY(xsmacnf);  NULLIFY(xsmackf)
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

!IF(.NOT. lAlloc) CALL AllocIsotopicMacXs()



CALL CP_CA(IsoXsMacA(1:niso, 1:ng), 0._8, niso, ng)
CALL CP_CA(IsoXsMacf(1:niso, 1:ng), 0._8, niso, ng)
CALL CP_CA(IsoXsMacs0(1:niso, 1:ng), 0._8, niso, ng)  !n2n Reaction for Temp
CALL CP_CA(IsoXsMacTr(1:niso, 1:ng), 0._8, niso, ng)  !n3n Reaction for Temp

DO iso = 1, niso
  id = mapnucl(idiso(iso));   isodata => ldiso(id)
  idn2n=mapn2n(id);idn3n=mapn3n(id)
  !Temperature Interpolation
  CALL XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2) 
  DO ig = 1, ng
    IsoXsMacA(iso, ig) = pnum(iso) * (wt2 * isodata%siga(ig, it2) + wt1 * isodata%siga(ig, it1))
    if (isodata%ifis.ne.0) IsoXsMacf(iso, ig) = pnum(iso) * (wt2 * isodata%sigf(ig, it2) + wt1 * isodata%sigf(ig, it1))
    IF(idn2n .NE. 0) IsoXsMacs0(iso, ig) = pnum(iso) * isodata%sign2n(ig)
    IF(idn3n .NE. 0) IsoXsMactr(iso, ig) = pnum(iso) * isodata%sign3n(ig)
  ENDDO
  NULLIFY(isodata) !Free the pointing variable
ENDDO !Isotopic sweeep

!Group Condense
CALL CP_CA(xsa, 0._8, niso)
CALL CP_CA(xsf, 0._8, niso)
CALL CP_CA(xsn2n, 0._8, niso)
CALL CP_CA(xsn3n, 0._8, niso)

phisum = 0
DO ig= 1, ng
  DO iso = 1, niso
    xsa(iso) = xsa(iso) + IsoXsMacA(iso, ig) * phi(ig)
    xsf(iso) = xsf(iso) + IsoXsMacf(iso, ig) * phi(ig)
    xsn2n(iso) = xsn2n(iso) + IsoXsMacs0(iso, ig) * phi(ig)
    xsn3n(iso) = xsn3n(iso) + IsoXsMactr(iso, ig) * phi(ig)
  ENDDO
  phisum = phisum + phi(ig)
ENDDO
DO iso = 1, niso
  xsa(iso) = xsa(iso) / phisum; xsf(iso) = xsf(iso) / phisum
  xsn2n(iso) = xsn2n(iso) / phisum; xsn3n(iso) = xsn3n(iso) / phisum
ENDDO

NULLIFY(isodata)
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
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
REAL :: temp, eigv
INTEGER :: niso
INTEGER :: idiso(niso)
REAL :: pnum(niso)
INTEGER :: ig1, ig2, ng
LOGICAL :: kin

REAL, POINTER :: xsmacnf(:)
TYPE(LIBDATA), POINTER :: isodata
INTEGER :: iso, id, ig
REAL :: reigv

INTEGER :: it1, it2
REAL :: wt1, wt2

IF(.NOT. XsMac%lalloc) THEN
  XsMac%ng = ng
  CALL AllocXsMac(XsMac)
ENDIF  

XsMacNf => XsMac%XsMacNf
CALL CP_CA(XsMacNf(ig1:ig2), 0._8, ig2 - ig1 + 1)

reigv = 1._8 / eigv

DO iso = 1, niso
  id = MapNucl(idiso(iso)); isodata => ldiso(id)
  IF(mapfis(idiso(iso)) .EQ. 0) CYCLE
  CALL XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2)
  DO ig = ig1, ig2
    XsMacNf(ig) = XsMacNf(ig) + pnum(iso) * (wt1 * isodata%signf(ig, it1) + wt2 * isodata%signf(ig, it2))
  ENDDO
ENDDO

IF(kin) CALL  MULTI_CA(reigv, XsMacNf(ig1 : ig2), ig2 - ig1 + 1)
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
CALL CP_CA(XsMacSm, 0._8, ng, ng)
CALL CP_CA(XsMacS, 0._8, ng)
CALL CP_CA(XsMacStr, 0._8, ng)
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
REAL :: wt1, wt2, ND
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
CALL CP_CA(XsMacSm, 0._8, ng, ng)
CALL CP_CA(XsMacS, 0._8, ng)
CALL CP_CA(XsMacStr, 0._8, ng)
    
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
  

DO ig = igb, ige
  IF(lScat1) XsMacSm(ig, ig) = XsMacSm(ig, ig) + (XsMacS(ig) - XsMacStr(ig))
ENDDO

NULLIFY(XsMacSm, XsMacs, XsMacStr)
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

CALL CP_CA(XsMacP1Sm, 0._8, ng, ng)
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


CALL CP_CA(XsMacP2Sm, 0._8, ng, ng)

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


CALL CP_CA(XsMacP3Sm, 0._8, ng, ng)

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
REAL :: wt1, wt2, sigs1, sigs2, sigssum
REAL, POINTER :: XsMacstr(:)

IF(.NOT. XsMac%lalloc) THEN
  XsMac%ng = ng
  CALL AllocXsMac(XsMac)
ENDIF  
!IF(.NOT. lAlloc) CALL AllocIsotopicMacXs()

XsMacStr => XsMac%XsMacStr
XsMacStr(ig) = 0._8
SigsSum = 0
DO iso = 1, niso
  id = mapnucl(idiso(iso)); isodata => ldiso(id)
  CALL XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2)
  ScRange1 = (/isodata%sm(ig, it1)%ioutsb, isodata%sm(ig, it1)%ioutse/)
  ScRange2 = (/isodata%sm(ig, it2)%ioutsb, isodata%sm(ig, it2)%ioutse/)  
  ioutbeg = MIN(ScRange1(1), ScRange2(1)); ioutend = MAX(ScRange1(2), ScRange2(2))
  !Out Scattering Range
  DO ig2 = ioutbeg, ioutend
  sigs1 =0; sigs2 = 0
  ind = (ig2 - ScRange2(1)) * (ig2 - ScRange2(2))
  ind2 = (ig - IsoData%sm(ig2, it2)%ib) * (ig - IsoData%sm(ig2, it2)%ie)
  IF(ind .le. 0 .and. ind2 .le. 0) sigs2 = isodata%sm(ig2, it2)%from(ig)
  ind = (ig2 - ScRange1(1)) * (ig2 - ScRange1(2))
  ind2 = (ig - IsoData%sm(ig2, it1)%ib) * (ig - IsoData%sm(ig2, it1)%ie)
  IF(ind .le. 0 .and. ind2 .le. 0) sigs1 = isodata%sm(ig2, it1)%from(ig)
  sigssum = sigssum + pnum(iso) * (wt1 * sigs1 + wt2 * sigs2)
  ENDDO
ENDDO

XsMacStr(ig) = SigsSum

NULLIFY(XsMacStr); NULLIFY(isodata)
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
REAL :: escxs(20), ind, jnd, adjintmlg, maclv_log, lvabs_log, wgtlvabs, phimult, phihom
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
  
  CALL calcWgt( TempAvgsq, isodata%rtempsq, isodata%nrtemp, wgttemp, idxtemp, 'A' )
      
  DO ig = ig1, ig2
    ind = pnum(iso)
    siglpiso = ind * isodata%lamsigp(ig)    
    DO jso = 1, niso
      IF ( iso .EQ. jso ) CYCLE
      jd = mapnucl(idiso(jso)); jsodata => ldiso(jd)
      jdres = mapnuclres(idiso(jso),iz)
      if (idres.eq.jdres) ind = ind + pnum(jso)
      IF ( jsodata%lreso .AND. ( .NOT. jsodata%lclad ) ) THEN
          if (isodata%rifid(jd).ne.0) CYCLE
          siglpiso = siglpiso + pnum(jso) * jsodata%lamsigp(ig)
      ELSE
          siglpiso = siglpiso + pnum(jso) * jsodata%lamsigp(ig)  
      ENDIF
    ENDDO ! DO jso = 1, niso
    if (ind.le.0._8) cycle
    
    nlv = isodata%nlv; nlvflx = isodata%nlvflx
    
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
          wgtabs = LineIntPol( TempAvgsq, isodata%nrtemp, isodata%rtempsq, isodata%wgtabs(ilv,:,ig) )
          wgtlvabs = wgtabs * lvabs; phimult = 1._8 / (lvabs + micsigb); phihom = micsigb * phimult
          phi = phi - wgtlvabs * phimult
          xsa = xsa + wgtlvabs * phihom
          if ( isodata%ityp .NE. 3 ) cycle
          lvfis = isodata%lvfis(ilv,ig)
          wgtfis = LineIntPol( TempAvgsq, isodata%nrtemp, isodata%rtempsq, isodata%wgtfis(ilv,:,ig) ) 
          xsf = xsf + wgtfis * lvfis * phihom
        enddo
      else
        do ilv = 1, nlv
          lvabs = isodata%lvabs(ilv,ig)  
          maclv_log = dlog(lvabs*ind)
          var(1:nmlv1g) = mypin%avgxseq_1g(1:nmlv1g)
          escxs(ilv) = LineIntPol2( maclv_log, nmlv1g, mlgmaclv, var(1:nmlv1g) )
          micsigb = (escxs(ilv) + siglpiso) / ind
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
      
    else
        
      if (nTracerCntl%lCAT) then  
          
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
        
      else
          
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
        
    xsa = xsa / phi
    if (xsa.le.0._8) cycle
    if ( isodata%ityp .EQ. 3 ) xsf = xsf / phi
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
      jnd = 0._8
      DO kso = 1, niso
        kdres = mapnuclres(idiso(kso),iz)
        if (jdres.eq.kdres) jnd = jnd + pnum(kso)
      ENDDO 
      IF ( jnd.le.0._8 ) CYCLE
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
        
        if ( isodata%ityp .NE. 3 ) CYCLE
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

SUBROUTINE EffMacXS(XsMac, mypin, myFxr, tempref, niso, ig, ng, lIsoXsOut, iz, PE)
USE TYPEDEF, ONLY : FxrInfo_Type, Pin_Type, ResVarPin_Type, PE_Type
USE CNTL, ONLY : nTracerCntl
USE OMP_LIB
#ifdef __PGI
USE IEEE_ARITHMETIC   !--- CNJ Edit : F2003 Standard
#endif
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac     !Microscopic XS
TYPE(PE_Type) :: PE
TYPE(ResVarPin_Type) :: mypin
TYPE(FxrInfo_Type) :: myFxr
TYPE(libdata), POINTER :: isodata, jsodata, repisodata

INTEGER,INTENT(IN) :: ig, niso, ng, iz
REAL,INTENT(IN) :: tempref
LOGICAL,INTENT(IN) :: lIsoXsOut

INTEGER :: iso, jso, id, jd, idres, jdres, repnid, repid, icat, idxsig(2,2), ir, kso, idxrat(2), kdres
INTEGER :: nlv, nlvflx, ilv, it, it1, it2, idxtemp(2), npot, nRiTemp, nmaclv, ibon, nThread
INTEGER :: rifid, igt
INTEGER,PARAMETER :: nbonmax = 10000
INTEGER,POINTER :: idiso(:)

REAL :: temp, TempAvgsq, Tempsq, wt1, wt2, wgttemp(2), invFnAdj
REAL :: siglp, siglp_noreso, bgxs, ind, lvabs, wgtabs, xsa, phi, logxsa, xss, xsss, xsf, lvfis, wgtfis
REAL :: miclv, miclv_log, lvabs_log, maclv_log, maclv, adjintmlg, Ft, dum2, f, repria, ria, rs, rifa(2), rifs(2), riff(2)
REAL :: wgtrat(2), jnd, rifa_iso, rifs_iso, riff_iso, rifrat, lograt
REAL :: SumSigA, SigaAnother, phimult, phihom, fl2, fl1, invfl1, wgtsig(2,2), s(2), wgtlvabs, var(200)
REAL,DIMENSION(niso, ng) :: IsoXsMacA, IsoXsMacS, IsoXsMacSS, IsoXsMacS1, IsoXsMacF
REAL,DIMENSION(niso) :: isoxsmaca_bon, isoxsmacf_bon, xsmaca0, phiiso, ratind
REAL,DIMENSION(nlvmax,niso) :: FtIso, micsigb, macsigb, wgtabs_bon, wgtfis_bon, maclvabs, maclvfis
REAL,POINTER :: xsmaca(:), xsmacf(:), xsmacs(:), xsmacstr(:), pnum(:), mlgmaclv(:)
REAL(4),POINTER :: xseq(:)

TYPE(riftype), POINTER :: isrif

IF(.NOT. XsMac%lalloc) THEN
  XsMac%ng = ng
  CALL AllocXsMac(XsMac)
ENDIF  

IF( lIsoXsOut .AND. ( .NOT. XsMac%lIsoAlloc ) ) THEN
  CALL AllocMacIsoXs(XsMac, ng, nelthel)
ENDIF

TempAvgsq = tempref
temp = myFxr%temp
Tempsq = dsqrt(temp)
if (myFxr%lCLD) TempAvgsq=Tempsq

pnum => myFxr%pnum; idiso => myFxr%idiso

xsmaca => XsMac%xsmaca; xsmacf => XsMac%xsmacf
xsmacs => XsMac%xsmacs; xsmacstr => XsMac%xsmacstr

CALL CP_CA(IsoXsMacA(1:niso, ig), 0._8, niso)
CALL CP_CA(IsoXsMacS(1:niso, ig), 0._8, niso)
CALL CP_CA(IsoXsMacS1(1:niso, ig), 0._8, niso)
CALL CP_CA(IsoXsMacSS(1:niso, ig), 0._8, niso)
CALL CP_CA(IsoXsMacF(1:niso, ig), 0._8, niso)

DO iso = 1, niso
  id = mapnucl(idiso(iso));   isodata => ldiso(id)
  CALL XsTempInterpolation(id, isodata, temp, wt1, wt2, it1, it2) 
  IsoXsMacA(iso, ig) = pnum(iso) * (wt2 * isodata%siga(ig, it2) + wt1 * isodata%siga(ig, it1))
  IsoXsMacS(iso, ig) = pnum(iso) * (wt2 * isodata%sigs(ig, it2) + wt1 * isodata%sigs(ig, it1))
  IsoXsMacS1(iso, ig) = IsoXsMacS(iso, ig) - pnum(iso) * (wt2 * isodata%sigstr(ig, it2) + wt1 * isodata%sigstr(ig, it1))
  IsoXsMacSS(iso, ig) = pnum(iso) * (wt2 * isodata%sigss(ig, it2) + wt1 * isodata%sigss(ig, it1))
  IF ( isodata%ifis .NE. 0 ) IsoXsMacF(iso, ig) = pnum(iso) * (wt2 * isodata%sigf(ig, it2) + wt1 * isodata%sigf(ig, it1))
ENDDO
    
IF (nTracerCNTL%lRIF) THEN
    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RIF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  FtIso(1:nlvmax,1:niso) = 1._8
  micsigb(1:nlvmax,1:niso) = 0._8

  siglp_noreso = 0._8
  DO jso = 1, niso
    jd = mapnucl(idiso(jso)); jsodata => ldiso(jd)
    IF ( .not.jsodata%lreso ) siglp_noreso = siglp_noreso + pnum(jso) * jsodata%lamsigp(ig)
  ENDDO 
    
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
      IF (.NOT.isodata%lreso) CYCLE
      idres = mapnuclRes(idiso(iso),iz)
      ind = 0._8; siglp = siglp_noreso
      DO jso = 1, niso
        jd = mapnucl(idiso(jso)); jsodata => ldiso(jd)
        jdres = mapnuclres(idiso(jso),iz)
        IF (idres.eq.jdres) THEN
          ind = ind + pnum(jso)
          siglp = siglp + pnum(jso) * jsodata%lamsigp(ig)  
        ENDIF
      ENDDO ! DO jso = 1, niso   
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
    
  ELSE
      
    IF (nTracerCNTL%lCAT) THEN 
      !******************************** CAT ********************************!  
      DO iso = 1, niso
        id = mapnucl(idiso(iso)); isodata => ldiso(id)
        IF ( .NOT. isodata%lreso ) CYCLE       
        idres = mapnuclRes(idiso(iso),iz)
        ind = 0._8; siglp = siglp_noreso
        DO jso = 1, niso
          jd = mapnucl(idiso(jso)); jsodata => ldiso(jd)
          jdres = mapnuclres(idiso(jso),iz)
          IF (idres.eq.jdres) THEN
            ind = ind + pnum(jso)
            siglp = siglp + pnum(jso) * jsodata%lamsigp(ig)  
          ENDIF
        ENDDO ! DO jso = 1, niso            
        IF (ind.le.0._8) CYCLE       
        nlv = isodata%nlv; nlvflx = isodata%nlvflx          
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
        ind = 0._8; siglp = siglp_noreso
        DO jso = 1, niso
          jd = mapnucl(idiso(jso)); jsodata => ldiso(jd)
          jdres = mapnuclres(idiso(jso),iz)
          IF (idres.eq.jdres) THEN
            ind = ind + pnum(jso) 
            siglp = siglp + pnum(jso) * jsodata%lamsigp(ig)
          ENDIF
        ENDDO ! DO jso = 1, niso
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
    IF (.NOT.isodata%lreso) CYCLE        
    IF (pnum(iso).le.0._8) CYCLE
    idres = mapnuclres(idiso(iso),iz)  
    !------- JSU EDIT 2019-0405
    ind = 0._8
    DO jso = 1, niso
      jd = mapnucl(idiso(jso)); jsodata => ldiso(jd)
      jdres = mapnuclres(idiso(jso),iz)
      IF (idres.eq.jdres) THEN
        ind = ind + pnum(jso) 
      ENDIF
    ENDDO ! DO jso = 1, niso
    IF (ind.le.0._8) CYCLE
    !------- END OF JSU EDIT 2019-0405
    ! Calculation of Effective XS of nuclides alone...
    xsa = 0._8; xsf = 0._8; phi = 1._8
    nlv = isodata%nlv
    do ilv = 1, nlv
      lvabs = isodata%lvabs(ilv,ig)
      
      wgtabs = LineIntPol( TempAvgsq, isodata%nrtemp, isodata%rtempsq, isodata%wgtabs(ilv,:,ig) ) 
      wgtlvabs = wgtabs * lvabs; phimult = 1._8 / (lvabs + micsigb(ilv,iso)); phihom = FtIso(ilv,iso) * micsigb(ilv,iso) * phimult
      phi = phi - wgtlvabs * phimult
      xsa = xsa + wgtlvabs * phihom
      if ( isodata%ityp .NE. 3 ) cycle
      lvfis = isodata%lvfis(ilv,ig)
      wgtfis = LineIntPol( TempAvgsq, isodata%nrtemp, isodata%rtempsq, isodata%wgtfis(ilv,:,ig) )
      xsf = xsf + wgtfis * lvfis * phihom
    enddo
    xsa = xsa / phi
    if (xsa.le.0._8) cycle

    if (.not.myFxr%lCLD) then
        
      IsoXsMacA(iso, ig) = xsa * pnum(iso) * myPin%rifa(idres,ig) 
      if ( isodata%ityp .EQ. 3 ) then
        xsf = xsf / phi
        IsoXsMacF(iso, ig) = xsf * pnum(iso) * myPin%riff(idres,ig) 
      endif
    
      IF (nTracerCNTL%lRST) THEN
        logxsa = dlog(xsa)
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
        IsoXsMacS(iso, ig) = xss * pnum(iso) * myPin%rifs(idres,ig) 
        IsoXsMacS1(iso, ig) = IsoXsMacS(iso, ig) * isodata%mu
        IsoXsMacSS(iso, ig) = xsss * pnum(iso) * myPin%rifs(idres,ig) 
      ENDIF  
      
    else
        
      logxsa = dlog(xsa)    
      CALL calcWgt( Tempsq, isodata%rtempsq, isodata%nrtemp, wgttemp, idxtemp, 'A' )  
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
        jnd = 0._8
        DO kso = 1, niso
          kdres = mapnuclres(idiso(kso),iz)
          if (jdres.eq.kdres) jnd = jnd + pnum(kso)
        ENDDO 
        IF ( jnd.le.0._8 ) CYCLE
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
        
          if ( isodata%ityp .NE. 3 ) CYCLE
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
        xsf = xsf / phi
        IsoXsMacF(iso, ig) = xsf * pnum(iso) * (riff(1) * wgttemp(1) + riff(2) * wgttemp(2)) 
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
        
    endif
  ENDDO
  
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
      
      SumSigA = 0._8; invfl1 = 1._8 / fl1
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
      
xsmaca(ig) = 0.; xsmacs(ig)  = 0.; xsmacf(ig) = 0.
DO iso = 1, niso
  xsmaca(ig) = xsmaca(ig) + IsoXsMacA(iso, ig)
  xsmacs(ig) = xsmacs(ig) + IsoXsMacS(iso, ig)
  xsmacf(ig) = xsmacf(ig) + IsoXsMacF(iso, ig)
ENDDO
xsmacstr(ig) = xsmacs(ig)
DO iso = 1, niso
  xsmacstr(ig) = xsmacstr(ig) - IsoXsMacS1(iso, ig)
ENDDO

IF(lIsoXsOut) THEN
  CALL CP_VA(XsMac%IsoXsMacA(1:niso, ig), IsoXsMacA(1:niso, ig), niso)
  CALL CP_VA(XsMac%IsoXsMacS0(1:niso, ig), IsoXsMacS(1:niso, ig), niso)
  CALL CP_VA(XsMac%IsoXsMacS1(1:niso, ig), IsoXsMacS1(1:niso, ig), niso)
  CALL CP_VA(XsMac%IsoXsMacF(1:niso, ig), IsoXsMacF(1:niso, ig), niso)
  CALL CP_VA(XsMac%IsoXsMacSS(1:niso, ig), IsoXsMacSS(1:niso, ig), niso)
ENDIF

NULLIFY(xsmaca, xsmacf, xsmacs, xsmacstr)
NULLIFY(isodata, jsodata, repisodata)
NULLIFY(pnum, idiso, mlgmaclv, xseq)

END SUBROUTINE

SUBROUTINE GetMacChi_Fxr(Fxr, phi, ig1, ig2, nchi, ng)
IMPLICIT NONE
TYPE(FxrInfo_Type) :: Fxr
INTEGER :: ig1, ig2, nchi, ng
REAL :: phi(ng)
CALL GetMacChi_Gen(Fxr%chi, Fxr%lres, Fxr%fresoFIso, Fxr%Temp, Fxr%niso, Fxr%idiso, Fxr%pnum, phi, ig1, ig2, nchi, ng)  
END SUBROUTINE

SUBROUTINE GetMacChi_Gen(chi, lres, fResoFIso, Temp, niso, idiso, pnum, phi, ig1, ig2, nchi, ng) 
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
CALL CP_CA(IsoXsMacNf(1:niso, 1:ng), 0._8, niso, ng)
!Initialize the chi 
CALL CP_CA(chi, 0._8, nchi)
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
    chi(ig) = chi(ig) + localpsi * isodata%chi(ig)
  ENDDO
ENDDO

DO ig = ig1, min(nchihel, ig2)
  chi(ig) = chi(ig) / psi
ENDDO
IF(ASSOCIATED(isodata)) NULLIFY(isodata)

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
CALL CP_CA(XsMacKf(ig1:ig2), 0._8, ig2 - ig1 + 1)

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
END MODULE
