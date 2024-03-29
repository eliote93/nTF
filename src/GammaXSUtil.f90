#include <defines.h>
#ifdef __GAMMA_TRANSPORT
MODULE GamXSUtil
! XS UTILITIES MODULE
!       : ALLOCATION, INTERPOLATION, TEMPORARY DATA
USE PARAM
USE GammaTYPEDEF,    ONLY : GamMacXS_TYPE
USE XSLIB_MOD,       ONLY : GAMLIBDATA
USE allocs
USE TYPEDEF,         ONLY : XsMac_Type

IMPLICIT NONE

INTEGER, PARAMETER, PRIVATE ::ndatmax = 200
TYPE(GamMacXS_TYPE), POINTER, PRIVATE :: XsMacDat(:)
LOGICAL, PRIVATE :: luse(ndatmax) = .FALSE.
LOGICAL, PRIVATE :: lAlloc_XsMacDat = .FALSE.

CONTAINS

! TEMPORARY DATA TREATMENT
!
SUBROUTINE GetGamXsMacDat(XsMac, ngg, lIsoXsOut)
USE XSLib_MOD,      ONLY : nelthel
IMPLICIT NONE
TYPE(GamMacXS_TYPE), POINTER :: XsMac
LOGICAL :: lIsoXsOut
INTEGER :: ngg

INTEGER :: i
LOGICAL :: lget

!$OMP CRITICAL
IF(.NOT. lAlloc_XsMacDat) THEN
  ALLOCATE(XsMacDat(nDatMax))
  lAlloc_XsMacDat = .TRUE.
ENDIF
lget = .FALSE.
DO i = 1, nDatMax
  IF(lUse(i)) CYCLE
  XsMacDat(i)%id = i;
  lUse(i) = .TRUE.
  XsMac => XsMacDat(i)
  IF(XsMac%lalloc .AND. XsMac%ngg .NE. ngg) THEN
    CALL FreeGamXsMac(XsMac)
    IF(XsMac%lisoalloc) CALL FreeGamXsIsoMac(xsMac)
    IF(XsMac%lIsoSMAlloc) CALL FreeGamSMISO(XsMac)
  ENDIF
  IF(.NOT. XsMac%lAlloc) THEN
    XsMac%ngg = ngg; CALL AllocGamMacXs(XsMac)
  ENDIF
  IF(lIsoXsOut .AND. .NOT. XsMac%lIsoAlloc) THEN
    XsMac%ngg = ngg
    XsMac%niso = nelthel
    CALL AllocGamMacIsoXs(XsMac)
  END IF
  EXIT
ENDDO
!$OMP END CRITICAL
END SUBROUTINE

SUBROUTINE ReturnGamXsMacDat(XsMac)
IMPLICIT NONE
TYPE(GamMacXS_TYPE), POINTER :: XsMac
!$OMP CRITICAL
lUse(XsMac%id) = .FALSE.
NULLIFY(XsMac)
!$OMP END CRITICAL
END SUBROUTINE

SUBROUTINE FreeGamXsMac(XsMac)
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XsMac
IF(.NOT. ASSOCIATED(XsMac%XsMacT)) RETURN
DEALLOCATE(XsMac%MacKERMA, XsMac%XsMacT)
DEALLOCATE(XsMac%XsMacTR, XsMac%XsMacA)
DEALLOCATE(XsMac%XsMacS, XsMac%XsMacSTR)
DEALLOCATE(XsMac%XsMacSM)
IF(XsMac%lAllocSm)THEN
    DEALLOCATE(XsMac%MacGSM1, XsMac%MacGSM2, XsMac%MacGSM3)
    XsMac%lAlloCSM=.FALSE.
ENDIF
XsMac%lAlloc = .FALSE.
XsMac%ngg = 0
END SUBROUTINE

SUBROUTINE FreeGamSMISO(XsMac)
TYPE(GamMacXS_TYPE) :: XsMac
IF(.NOT. XsMac%lIsoSMAlloc) RETURN
DEALLOCATE(XsMac%IsoSM)
XsMac%lIsoSMAlloc = .FALSE.
END SUBROUTINE

SUBROUTINE FreeGamXsIsoMac(XsIsoMac)
TYPE(GamMacXS_TYPE) :: XsIsoMac
IF(.NOT. XsIsoMac%lIsoAlloc) RETURN
DEALLOCATE(XsIsoMac%IsoKERMA);  DEALLOCATE(XsIsoMac%IsoXsMacT)
DEALLOCATE(XsIsoMac%IsoXsMacA);    DEALLOCATE(XsIsoMac%IsoXsMacTR)
DEALLOCATE(XsIsoMac%IsoXsMacS);    DEALLOCATE(XsIsoMac%IsoXsMacSTR)
XsIsoMac%lIsoAlloc = .FALSE.
XsIsoMac%niso = 0
END SUBROUTINE
!***************************************************************************************************
!   Allocate Gamma Cross Section Type arrays and matirces
!
SUBROUTINE AllocGamMacXs(XsMac)
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XsMac
INTEGER :: ngg
ngg = XsMac%ngg
ALLOCATE(XsMac%MacKERMA(ngg), XsMac%XsMacT(ngg))
ALLOCATE(XsMac%XsMacTR(ngg), XsMac%XsMacA(ngg))
ALLOCATE(XsMac%XsMacS(ngg), XsMac%XsMacSTR(ngg))
ALLOCATE(XsMac%XsMacSM(ngg,ngg))
ALLOCATE(XsMac%MacGSM1(ngg,ngg), XsMac%MacGSM2(ngg,ngg), XsMac%MacGSM3(ngg,ngg))
XsMac%lallocsm = .TRUE.
XsMac%lalloc = .TRUE.
END SUBROUTINE

SUBROUTINE AllocGamMacIsoXs(XsMac)
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XsMac
INTEGER :: ngg, niso
ngg = XsMac%ngg
niso = XsMac%niso
ALLOCATE(XsMac%IsoKERMA(niso, ngg), XsMac%IsoXsMacT(niso, ngg), XsMac%IsoXsMacTR(niso, ngg))
ALLOCATE(XsMac%IsoXsMacA(niso, ngg), XsMac%IsoXsMacS(niso, ngg), XsMac%IsoXsMacSTR(niso, ngg))
XsMac%lIsoAlloc = .TRUE.
END SUBROUTINE

SUBROUTINE AllocGamIsoSM(XsMac)
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XsMac
INTEGER :: ngg, niso
ngg = XsMac%ngg
niso = XsMac%niso
ALLOCATE(XsMac%IsoSM(ngg, ngg, niso))
XsMac%lIsoSMAlloc = .TRUE.
END SUBROUTINE

SUBROUTINE AllocGamProdMat(XSMAT)
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XSMAT
INTEGER :: igg, ngg, ng, ig, niso
ngg = XSMAT%ngg
ng = XSMAT%ng
niso = XSMAT%niso
ALLOCATE(XSMAT%lexist(niso, 4))
ALLOCATE(XSMAT%GProdTot(ng,ngg))
XSMAT%lProdAlloc = .TRUE.
END SUBROUTINE

SUBROUTINE FreeGamProdMat(XSMAT)
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XSMAT
DEALLOCATE(XSMAT%lexist)
DEALLOCATE(XSMAT%GProdTot)
XSMAT%lProdAlloc = .FALSE.
END SUBROUTINE

SUBROUTINE AllocIsoGamProdMat(XSMAT)
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XSMAT
INTEGER :: igg, ngg, ng, ig, niso
ngg = XSMAT%ngg
ng = XSMAT%ng
niso = XSMAT%niso
ALLOCATE(XSMAT%IsoGProdFis(niso,ng,ngg))
ALLOCATE(XSMAT%IsoGProdRad(niso,ng,ngg))
ALLOCATE(XSMAT%IsoGProdInel(niso,ng,ngg))
ALLOCATE(XSMAT%IsoGProdNnel(niso,ng,ngg))
XSMAT%lIsoProdAlloc = .TRUE.
END SUBROUTINE

SUBROUTINE FreeIsoGamProdMat(XSMAT)
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XSMAT
DEALLOCATE(XSMAT%IsoGProdFis)
DEALLOCATE(XSMAT%IsoGProdRad)
DEALLOCATE(XSMAT%IsoGProdInel)
DEALLOCATE(XSMAT%IsoGProdNnel)
XSMAT%lIsoProdAlloc = .FALSE.
END SUBROUTINE

! ALLOCATIONE SUBROUTINE TO CALCULATE LOCAL KAPPA                         |-- JSU EDIT 2017.09.14. |
!SUBROUTINE AllocLocalKappa(XsMac) 
!IMPLICIT NONE
!TYPE(GamMacXs_TYPE) :: XsMac
!INTEGER :: igg, ngg, ig, ng, niso, iso
!ng = XsMac%ng
!ngg = XsMac%ngg
!niso = XsMac%niso
!ALLOCATE(XsMac%FisLocal(ng))     
!ALLOCATE(XsMac%QN2N(ng), XsMac%QN3N(ng))
!ALLOCATE(XsMac%InelLocal(ng))    
! ALLOCATE(XsMac%InelLosS(ng, ng))  
! ALLOCATE(XsMac%GProdInel(ng, ngg)) 
!ALLOCATE(XsMac%LocalQ(ng))
!XsMac%lKappaAlloc = .TRUE.
!END SUBROUTINE

! ALLOCATIONE SUBROUTINE TO CALCULATE LOCAL KERMA |-- JSU EDIT 2019.08.12. |
SUBROUTINE AllocMacKERMA(XsMac) 
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
INTEGER :: igg, ngg, ig, ng, niso, iso
ng = XsMac%ng
niso = XsMac%niso
ALLOCATE(XsMac%MacKERMA_t(ng));ALLOCATE(XsMac%MacKERMA_s(ng))
ALLOCATE(XsMac%MacKERMA_d(ng));ALLOCATE(XsMac%MacKERMA_p(ng))
ALLOCATE(XsMac%MacKERMA_f(ng))
ALLOCATE(XsMac%MacDelkf(ng))
XsMac%lKERMAAlloc = .TRUE.
END SUBROUTINE

SUBROUTINE FreeMacKERMA(XsMac)
TYPE(XsMac_Type) :: XsMac
IF(.NOT. XsMac%lKERMAAlloc) RETURN
DEALLOCATE(XsMac%MacKERMA_t);DEALLOCATE(XsMac%MacKERMA_s)
DEALLOCATE(XsMac%MacKERMA_d);DEALLOCATE(XsMac%MacKERMA_p)
DEALLOCATE(XsMac%MacKERMA_f)
DEALLOCATE(XsMac%MacDelkf)
XsMac%lKERMAAlloc = .FALSE.
END SUBROUTINE

SUBROUTINE AllocIsoMacKERMA(XsMac) 
IMPLICIT NONE
TYPE(XsMac_Type) :: XsMac
INTEGER :: igg, ngg, ig, ng, niso, iso
ng = XsMac%ng
niso = XsMac%niso
ALLOCATE(XsMac%IsoMacKERMA_t(niso,ng))
ALLOCATE(XsMac%IsoMacKERMA_s(niso,ng))
ALLOCATE(XsMac%IsoMacKERMA_d(niso,ng))
ALLOCATE(XsMac%IsoMacKERMA_p(niso,ng))
ALLOCATE(XsMac%IsoMacKERMA_f(niso,ng))
ALLOCATE(XsMac%IsoMacDelkf(niso,ng))
XsMac%lISOKERMAAlloc = .TRUE.
END SUBROUTINE

SUBROUTINE FreeIsoMacKERMA(XsMac)
TYPE(XsMac_Type) :: XsMac
IF(.NOT. XsMac%lISOKERMAAlloc) RETURN
DEALLOCATE(XsMac%IsoMacKERMA_t)
DEALLOCATE(XsMac%IsoMacKERMA_s)
DEALLOCATE(XsMac%IsoMacKERMA_d)
DEALLOCATE(XsMac%IsoMacKERMA_p)
DEALLOCATE(XsMac%IsoMacKERMA_f)
DEALLOCATE(XsMac%IsoMacDelkf)
XsMac%lISOKERMAAlloc = .FALSE.
END SUBROUTINE
!***************************************************************************************************
!   Cusping Pin Interpolation
!
SUBROUTINE IntMacBaseCsp(XsMac, XsMac1, XsMac2, CspFxr, niso0, igb, ige, lIsoXsOut)
USE TYPEDEF,           ONLY : CspFxr_Type
USE BasicOperation,    ONLY : CP_CA
USE GammaTYPEDEF,      ONLY : GamMacXS_TYPE
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XsMac, XsMac1, XsMac2
TYPE(CspFxr_Type) :: CspFxr

INTEGER :: igb, ige, niso0
LOGICAL :: lIsoXsOut

REAL :: wt1, wt2
INTEGER :: ig, iso, i, iloc

DO ig = igb, ige
  wt1 = CspFxr%fcsp(ig, 1); wt2 = CspFxr%fcsp(ig, 2)
  XsMac%MacKERMA(ig) = XsMac1%MacKERMA(ig) * wt1 + XsMac2%MacKERMA(ig) * wt2
  XsMac%XsMacT(ig)   =   XsMac1%XsMacT(ig) * wt1 + XsMac2%XsMacT(ig)   * wt2
  XsMac%XsMacTr(ig)  =  XsMac1%XsMacTr(ig) * wt1 + XsMac2%XsMacTr(ig)  * wt2
  XsMac%XsMacA(ig)   =   XsMac1%XsMacA(ig) * wt1 + XsMac2%XsMacA(ig)   * wt2
  XsMac%XsMacS(ig)   =   XsMac1%XsMacS(ig) * wt1 + XsMac2%XsMacS(ig)   * wt2
  XsMac%XsMacStr(ig) = XsMac1%XsMacStr(ig) * wt1 + XsMac2%XsMacStr(ig) * wt2
ENDDO
IF(.NOT. lIsoXsOut) RETURN
CALL CP_CA(XsMac%IsoKERMA(1:niso0, igb:ige), 0._8, niso0, ige-igb+1)
CALL CP_CA(XsMac%IsoXsMacT(1:niso0, igb:ige), 0._8, niso0, ige-igb+1)
CALL CP_CA(XsMac%IsoXsMacTR(1:niso0, igb:ige), 0._8, niso0, ige-igb+1)
CALL CP_CA(XsMac%IsoXsMacA(1:niso0, igb:ige), 0._8, niso0, ige-igb+1)
CALL CP_CA(XsMac%IsoXsMacS(1:niso0, igb:ige), 0._8, niso0, ige-igb+1)
CALL CP_CA(XsMac%IsoXsMacSTR(1:niso0, igb:ige), 0._8, niso0, ige-igb+1)

DO ig = igb, ige
  DO iso = 1, CspFxr%niso(1)
    iloc = CspFxr%isomap(iso, 1)
    XsMac%IsoKERMA(iloc, ig) = XsMac1%IsoKERMA(iso, ig) * CspFxr%fcsp(ig, 1)
    XsMac%IsoXsMacT(iloc, ig) = XsMac1%IsoXsMacT(iso, ig) * CspFxr%fcsp(ig, 1)
    XsMac%IsoXsMacTR(iloc, ig) = XsMac1%IsoXsMacTR(iso, ig) * CspFxr%fcsp(ig, 1)
    XsMac%IsoXsMacA(iloc, ig) = XsMac1%IsoXsMacA(iso, ig) * CspFxr%fcsp(ig, 1)
    XsMac%IsoXsMacS(iloc, ig) = XsMac1%IsoXsMacS(iso, ig) * CspFxr%fcsp(ig, 1)
    XsMac%IsoXsMacSTR(iloc, ig) = XsMac1%IsoXsMacSTR(iso, ig) * CspFxr%fcsp(ig, 1)
  ENDDO
  DO iso = 1, CspFxr%niso(2)
    iloc = CspFxr%isomap(iso, 2)
    XsMac%IsoKERMA(iloc, ig) = XsMac%IsoKERMA(iloc, ig) + XsMac2%IsoKERMA(iso, ig) * CspFxr%fcsp(ig, 2)
    XsMac%IsoXsMacT(iloc, ig) = XsMac%IsoXsMacT(iloc, ig) + XsMac2%IsoXsMacT(iso, ig) * CspFxr%fcsp(ig, 2)
    XsMac%IsoXsMacTR(iloc, ig) = XsMac%IsoXsMacTR(iloc, ig) + XsMac2%IsoXsMacTR(iso, ig) * CspFxr%fcsp(ig, 2)
    XsMac%IsoXsMacA(iloc, ig) = XsMac%IsoXsMacA(iloc, ig) + XsMac2%IsoXsMacA(iso, ig) * CspFxr%fcsp(ig, 2)
    XsMac%IsoXsMacS(iloc, ig) = XsMac%IsoXsMacS(iloc, ig) + XsMac2%IsoXsMacS(iso, ig) * CspFxr%fcsp(ig, 2)
    XsMac%IsoXsMacSTR(iloc, ig) = XsMac%IsoXsMacSTR(iloc, ig) + XsMac2%IsoXsMacSTR(iso, ig) * CspFxr%fcsp(ig, 2)
  ENDDO
ENDDO

END SUBROUTINE

SUBROUTINE IntScatMatCsp(XsMac, XsMac1, XsMac2, CspFxr, niso0, igb, ige, ngg)
USE TYPEDEF,           ONLY : CspFxr_Type
USE BasicOperation,    ONLY : CP_CA
USE GammaTYPEDEF,      ONLY : GamMacXS_TYPE
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XsMac, XsMac1, XsMac2
TYPE(CspFxr_Type) :: CspFxr
INTEGER :: niso0, igb, ige, ngg

REAL :: wt1, wt2
INTEGER :: igg, igg2

DO igg = igb, ige
  DO igg2 = 1, ngg
    wt1 = CspFxr%fcsp(igg2,1); wt2 = CspFxr%fcsp(igg2,2)
    XsMac%XsMacSm(igg2, igg) = XsMac1%XsMacSm(igg2, igg) * wt1 + XsMac2%XsMacSm(igg2, igg) * wt2
  END DO
  wt1 = CspFxr%fcsp(igg,1); wt2 = CspFxr%fcsp(igg,2)
  XsMac%XsMacS(igg) = XsMac1%XsMacS(igg) * wt1 + XsMac2%XsMacS(igg) * wt2
END DO

END SUBROUTINE

SUBROUTINE IntProdMatCsp(XsMac, XsMac1, XsMac2, CspFxr, niso0, igb, ige, ng)
USE TYPEDEF,           ONLY : CspFxr_Type
USE BasicOperation,    ONLY : CP_CA
USE GammaTYPEDEF,      ONLY : GamMacXS_TYPE
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XsMac, XsMac1, XsMac2
TYPE(CspFxr_Type) :: CspFxr
INTEGER :: niso0, igb, ige, ng

REAL :: wt1, wt2
INTEGER :: igg, ig

XsMac%GProdTot = 0.
DO igg = igb, ige
  DO ig = 1, ng
    wt1 = CspFxr%fcsp(ig,1); wt2 = CspFxr%fcsp(ig,2)
    XsMac%GProdTot(ig, igg) = XsMac1%GProdTot(ig, igg) * wt1 + XsMac2%GProdTot(ig, igg) * wt2
  END DO
END DO

END SUBROUTINE

SUBROUTINE mappingphatom(nele,ntiso)
!  -------------------------------------------------------------------------------------------------
!  Mapping for Photoatomic reaction data...
!  -------------------------------------------------------------------------------------------------
USE PARAM
USE XSlib_MOD,  ONLY : nuclidhel,   elementGAM,   mapEleGAM, ldiso, phatom, iso2elm
USE ALLOCS,     ONLY : dmalloc

IMPLICIT NONE

INTEGER :: nele, ntiso
INTEGER :: i,j
INTEGER :: idelem

! Nuclide List
CALL dmalloc(elementGAM,nele)
CALL dmalloc(iso2elm,ntiso)

DO i = 1, nele
  elementGAM(i) = phatom(i)%nid
END DO

! isotope -> element index matching
Do j = 1, ntiso
  idelem = (nuclidhel(j)/1000)*1000
  Do i = 1, nele
    IF (elementGAM(i).EQ.idelem) EXIT
  END DO
  iso2elm(j) = elementGAM(i)
  ldiso(j)%Phatom => phatom(i)
End Do

mapEleGAM = 0
DO i=1,nele
  mapEleGAM(elementGAM(i)) = i
END DO

RETURN
END SUBROUTINE
      
END MODULE
#endif
