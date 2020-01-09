#include <defines.h>
#ifdef __GAMMA_TRANSPORT    
Module GamXsLib_Mod
!-- JSU EDIT
! MODULE THAT CONTAINS SUBROUTINES USED IN MOC CROSS SECTION AND SOURCE GENERATION
USE PARAM
USE TYPEDEF,         ONLY : Fxrinfo_type  
USE GamXSUtil,       ONLY : AllocGamMacXs,  AllocGamMacIsoXs,  AllocGamProdMat,                     &
                            GamXsTempInterpolation
USE XsUtil_mod,      ONLY :  XsTempInterpolation
!USE Allocs
!
IMPLICIT NONE

! INTERFACE
INTERFACE GamXsBase
  MODULE PROCEDURE GamXsBase_Fxr
  MODULE PROCEDURE GamXsBase_Csp
  MODULE PROCEDURE GamXsBase_gen
END INTERFACE

INTERFACE GamTotScatXs
  MODULE PROCEDURE GamTotScatXs_Fxr
  MODULE PROCEDURE GamTotScatXs_Gen
END INTERFACE

INTERFACE GamScatMatrix
  MODULE PROCEDURE GamScatMatrix_Fxr
  MODULE PROCEDURE GamScatMatrix_Csp
  MODULE PROCEDURE GamScatMatrix_Gen
END INTERFACE

INTERFACE GamP1XsScatMatrix
  MODULE PROCEDURE GamP1XsScatMatrix_Fxr
  MODULE PROCEDURE GamP1XsScatMatrix_Gen
END INTERFACE

INTERFACE GamP2XsScatMatrix
  MODULE PROCEDURE GamP2XsScatMatrix_Fxr
  MODULE PROCEDURE GamP2XsScatMatrix_Gen
END INTERFACE

INTERFACE GamP3XsScatMatrix
  MODULE PROCEDURE GamP3XsScatMatrix_Fxr
  MODULE PROCEDURE GamP3XsScatMatrix_Gen
END INTERFACE

INTERFACE GamProdMatrix
  MODULE PROCEDURE GamProdMatrix_Fxr
  MODULE PROCEDURE GamProdMatrix_Csp
  MODULE PROCEDURE GamProdMatrix_Gen
END INTERFACE

! EXPLICIT KAPPA CALCULATION                                              |-- JSU EDIT 2017.09.13. |
!- FOR PHOTON
!    3.b. PHOTON ENERGY DUE TO INELASTIC SCATTERING 
!INTERFACE GamInelMatrix
!  MODULE PROCEDURE GamInelMatrix_Fxr
!  MODULE PROCEDURE GamInelMatrix_Csp
!  MODULE PROCEDURE GamInelMatrix_Gen
!END INTERFACE

!- FOR NEUTRON PART
!    1. LOCALLY DEPOSITED FISSION ENERGY
!    2. N2N, N3N REACTION ENERGY LOSS (TOTAL RECOVERABLE ENERGY DECREASE
!    3.a. INELASTIC SCATTERING LOCAL DEPOSIT ENERGY (NEUTRON ENERGY LOSS MATRIX)
INTERFACE GetLocalQMat
  MODULE PROCEDURE GetLocalQMat_Fxr
!  MODULE PROCEDURE GetLocalQMat_Csp
  MODULE PROCEDURE GetLocalQMat_Gen
END INTERFACE

  CONTAINS
  
!***************************************************************************************************
! BASIC PHOTO-ATOMIC REACTION DATA
!
!   GamXsBase SUBROUTINES
!             GENERATES MACROSCOPIC CROSS SECTION OF PHOTO-ATOMIC REACTIONS
SUBROUTINE GamXsBase_Fxr(XsMac, Fxr, igb, ige, ngg, lIsoXsOut)
USE GammaTYPEDEF,   ONLY : GamMacXS_TYPE
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XsMac
TYPE(Fxrinfo_type) :: Fxr
INTEGER :: igb, ige, ngg
REAL :: eigv
LOGICAL :: kin, lIsoXsOut

CALL GamXsBase_Gen(XsMac, Fxr%niso, Fxr%idiso, Fxr%pnum, igb, ige, ngg, lIsoXsOut)
END SUBROUTINE

SUBROUTINE GamXsBase_Csp(XsMac, Fxr, igb, ige, ngg, lIsoXsOut, lCrCspFtn)
USE GamXSUtil,       ONLY :  GetGamXsMacDat,  ReturnGamXsMacDat, IntMacBaseCsp
USE GammaTYPEDEF,    ONLY : GamMacXS_TYPE
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XsMac
TYPE(Fxrinfo_type) :: Fxr
INTEGER :: igb, ige, ngg
REAL :: eigv
LOGICAL :: kin, lIsoXsOut, lCrCspFtn

TYPE(GamMacXS_TYPE), POINTER :: XsMac1, XsMac2
INTEGER :: i
INTEGER :: niso
INTEGER, POINTER :: idiso(:)
REAL, POINTER :: pnum(:)

IF(lCrCspFtn .AND. Fxr%lCrCspFtn) THEN
  CALL GetGamXsMacDat(XsMac1, ngg, lIsoXsout)
  CALL GetGamXsMacDat(XsMac2, ngg, lIsoXsout)
  
  niso = Fxr%CspFxr%niso(1)
  Pnum => Fxr%CspFxr%pnum(:, 1);     idiso => Fxr%CspFxr%isolist(:, 1)
  CALL GamXsBase_Gen(XsMac1, niso, idiso, pnum, igb, ige, ngg, lIsoXsOut)
  
  niso = Fxr%CspFxr%niso(2)
  Pnum => Fxr%CspFxr%pnum(:, 2);     idiso => Fxr%CspFxr%isolist(:, 2)
  CALL GamXsBase_Gen(XsMac2, niso, idiso, pnum, igb, ige, ngg, lIsoXsOut)
  
  ! Interpolation of Two Csp
  CALL IntMacBaseCsp(XsMac, XsMac1, XsMac2, Fxr%CspFxr, Fxr%niso, igb, ige, lIsoXsOut)
  CALL ReturnGamXsMacDat(XsMac1);        CALL ReturnGamXsMacDat(XsMac2)
ELSE
  niso = Fxr%niso
  pnum => Fxr%pnum;                 idiso => Fxr%idiso
  CALL GamXsBase_Gen(XsMac, niso, idiso, pnum, igb, ige, ngg, lIsoXsOut)
ENDIF
END SUBROUTINE

SUBROUTINE GamXsBase_Gen(XsMac, niso, idiso, pnum, igb, ige, ngg, lIsoXsOut)
USE GammaLibdata,   ONLY : neltGAM, Gldiso, mapnuclELM
USE GammaTYPEDEF,   ONLY : GamMacXS_TYPE, GAMLIBDATA
IMPLICIT NONE
! Input variables
TYPE(GamMacXS_TYPE) :: XsMac
INTEGER :: igb, ige, ngg, niso
INTEGER :: idiso(niso)
REAL :: pnum(niso)
LOGICAL :: lIsoXsOut

! Pointers
TYPE(GAMLIBDATA), POINTER :: Gisodata

! Local Variable 
REAL :: IsoKERMA(niso, ngg), IsoXsMacT(niso, ngg), IsoXsMacTR(niso, ngg), IsoXsMacA(niso, ngg), &
        IsoXsMacS(niso, ngg), IsoXsMacSTR(niso, ngg)
INTEGER :: ig, igg, id, it1, it2
REAL :: reigv
INTEGER :: iso

! allocate macroscopic data
IF (.NOT.XsMac%lalloc) THEN
  XsMac%ngg = ngg
  CALL AllocGamMacXs(XsMac) 
END IF
IF (.NOT.XsMac%lIsoAlloc) THEN
  XsMac%ngg = ngg
  XsMac%niso = neltGAM
  CALL AllocGamMacIsoXs(XsMac)
END IF
! Initialization
XsMac%MacKERMA = 0._8
XsMac%XsMacSTR = 0._8
XsMac%XsMacTR = 0._8
XsMac%XsMacT = 0._8
XsMac%XsMacA = 0._8
XsMac%XsMacS = 0._8

IsoKERMA = 0._8
IsoXsMacSTR = 0._8
IsoXsMacTR = 0._8
IsoXsMacT = 0._8
IsoXsMacA = 0._8
IsoXsMacS = 0._8
! Isotope-wise macro scopic cross section
!    Due to independent of temp., simply multiplying pnum
DO iso = 1, niso
  id = mapnuclELM(idiso(iso))
  IF (id .EQ. 0) THEN
    PRINT *, 'ISOTOPE', idiso(iso), 'DOES NOT EXIST IN GAMMA LIBRARY(element)'
    CYCLE
  END IF
  Gisodata => Gldiso(id)
  DO igg = igb, ige
    IsoKERMA(iso, igg) = pnum(iso) * Gisodata%KERMA(igg)
    IsoXsMacTR(iso, igg) = pnum(iso) * Gisodata%GSIGTR(igg)
    IsoXsMacA(iso, igg) = pnum(iso) * Gisodata%GSIGA(igg)
    IsoXsMacS(iso, igg) = pnum(iso) * Gisodata%GSIGS(igg)
    IsoXsMacT(iso, igg) = IsoXsMacA(iso, igg) + IsoXsMacS(iso, igg)
  END DO ! group loop
  NULLIFY(Gisodata) 
END DO ! iso loop
! Macroscopic cross section of FXR (Summing iso-wise)
DO igg = igb, ige
  XsMac%MacKERMA(igg) = SUM(IsoKERMA(1:niso, igg))
  XsMac%XsMacTR(igg) = SUM(IsoXsMacTR(1:niso, igg))
  XsMac%XsMacA(igg) = SUM(IsoXsMacA(1:niso, igg))
  XsMac%XsMacS(igg) = SUM(IsoXsMacS(1:niso, igg))
  XsMac%XsMacSTR(igg) = XsMac%XsMacTR(igg) - XsMac%XsMacA(igg)
  XsMac%XsMacT(igg) = XsMac%XsMacA(igg) + XsMac%XsMacS(igg)
END DO
! Isotope-wise output Handling
IF (lIsoXsOut) THEN
  XsMac%IsoKERMA(1:niso, igb:ige) = IsoKERMA(1:niso, igb:ige)
  XsMac%IsoXsMacT(1:niso, igb:ige) = IsoXsMacT(1:niso, igb:ige)
  XsMac%IsoXsMacTR(1:niso, igb:ige) = IsoXsMacTR(1:niso, igb:ige)
  XsMac%IsoXsMacA(1:niso, igb:ige) = IsoXsMacA(1:niso, igb:ige)
  XsMac%IsoXsMacS(1:niso, igb:ige) = IsoXsMacS(1:niso, igb:ige)
END IF
END SUBROUTINE

!   GamTotScatXs SUBROUTINES
!               GENERATES MACROSCOPIC TRANSPORT CORRE
!           *** Meaningless in photo-atomic rection!
SUBROUTINE GamTotScatXs_Fxr(XsMac, Fxr, igg, ngg)
USE GammaTYPEDEF,   ONLY : GamMacXS_TYPE
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XsMac
TYPE(FXRINFO_TYPE) :: Fxr
INTEGER :: igg, ngg
CALL GamTotScatXs_Gen(XsMac, Fxr%niso, Fxr%idiso, Fxr%pnum, igg, ngg)
END SUBROUTINE

SUBROUTINE GamTotScatXs_Gen(XsMac, niso, idiso, pnum, igg, ngg)
USE GammaLibdata,   ONLY : neltGAM, Gldiso, mapnuclELM
USE GammaTYPEDEF,   ONLY : GamMacXS_TYPE, GAMLIBDATA
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XsMac
INTEGER :: niso, igg, ngg
INTEGER :: idiso(niso)
REAL :: pnum(niso)

TYPE(GAMLIBDATA), POINTER :: gisodata
INTEGER :: iso , idelem, ig2, ioutbeg, ioutend
INTEGER :: ScRange(2)
REAL :: SigsSum, sigs0sum
REAL, POINTER :: XsMacSTR(:)

IF(.NOT. XsMac%lalloc) THEN
  XsMac%ngg = ngg
  CALL AllocGamMacXs(XsMac)
ENDIF  

XsMacStr => XsMac%XsMacSTR
XsMacSTR(igg) = 0._8
SigsSum = 0
sigs0sum = 0
DO iso = 1, niso
  idelem = mapnuclELM(idiso(iso))
  IF (idelem .EQ. 0) CYCLE
  gisodata => Gldiso(idelem)
  ScRange = (/gisodata%GSM(igg)%ioutsb, gisodata%GSM(igg)%ioutse/)
  ioutbeg = ScRange(1); ioutend = ScRange(2)
  !Out Scattering Range
  DO ig2 = ioutbeg, ioutend
    sigssum = sigssum + pnum(iso) * gisodata%GSM(ig2)%from(igg)
  ENDDO
  NULLIFY(gisodata)
ENDDO

XsMacStr(igg) = SigsSum
!XsMac%XsMacS(ig) = Sigs0sum

NULLIFY(XsMacStr)
END SUBROUTINE

!   GamScatMatrix SUBROUTINES
!                 GENERATES MACROSCOPIC 0-TH ORDER SCATTERING MATRIX
SUBROUTINE GamScatMatrix_Fxr(XsMac, FXR, igb, ige, ngg, lscat1)
USE GammaTYPEDEF,   ONLY : GamMacXS_TYPE
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XsMac
TYPE(FXRInfo_TYPE) :: FXR
INTEGER :: igb, ige, ngg
LOGICAL :: lscat1
CALL GamScatMatrix_Gen(XsMac, FXR%niso, FXR%idiso, FXR%pnum, igb, ige, ngg, lscat1)
END SUBROUTINE

SUBROUTINE GamScatMatrix_Csp(XsMac, FXR, igb, ige, ngg, lscat1, lCrCspFtn)
USE GamXSUtil,      ONLY :  GetGamXsMacDat,  ReturnGamXsMacDat,  IntScatMatCsp
USE GammaTYPEDEF,   ONLY : GamMacXS_TYPE
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XsMac
TYPE(FXRInfo_TYPE) :: FXR
INTEGER :: igb, ige, ngg
LOGICAL :: lscat1, lCrCspFtn

TYPE(GamMacXS_TYPE), POINTER :: XsMac1, XsMac2
INTEGER :: i
INTEGER :: niso
REAL :: temp
INTEGER, POINTER :: idiso(:)
REAL, POINTER :: pnum(:)

IF(lCrCspFtn .AND. Fxr%lCrCspFtn) THEN
  CALL GetGamXsMacDat(XsMac1, ngg, .FALSE.)
  CALL GetGamXsMacDat(XsMac2, ngg, .FALSE.)
  
  niso = Fxr%CspFxr%niso(1)
  Pnum => Fxr%CspFxr%pnum(:, 1);     idiso => Fxr%CspFxr%isolist(:, 1)
  CALL GamScatMatrix_Gen(XsMac1, niso, idiso, pnum, igb, ige, ngg, lscat1)
  
  niso = Fxr%CspFxr%niso(2)
  Pnum => Fxr%CspFxr%pnum(:, 2);     idiso => Fxr%CspFxr%isolist(:, 2)
  CALL GamScatMatrix_Gen(XsMac2, niso, idiso, pnum, igb, ige, ngg, lscat1)
  
  ! Interpolation of Two Csp
  CALL IntScatMatCsp(XsMac, XsMac1, XsMac2, Fxr%CspFxr, Fxr%niso, igb, ige, ngg)
  CALL ReturnGamXsMacDat(XsMac1);        CALL ReturnGamXsMacDat(XsMac2)
  
ELSE
  niso = Fxr%niso;                  Temp = Fxr%Temp
  pnum => Fxr%pnum;                 idiso => Fxr%idiso
  CALL GamScatMatrix_Gen(XsMac, niso, idiso, pnum, igb, ige, ngg, lscat1)
ENDIF
END SUBROUTINE

SUBROUTINE GamScatMatrix_Gen(XsMac, niso, idiso, pnum, igb, ige, ngg, lscat1)
USE GammaLibdata,   ONLY : Gldiso, mapnuclELM
USE GammaTYPEDEF,   ONLY : GamMacXS_TYPE, GAMLIBDATA
IMPLICIT NONE
! Input Variables
TYPE(GamMacXS_TYPE) :: XsMac
INTEGER :: niso
INTEGER :: idiso(niso)
REAL :: pnum(niso)
INTEGER :: igb, ige, ngg
LOGICAL :: lscat1
! Pointers
REAL, POINTER :: Xsmacsm(:,:), XsMacS(:), XsMacStr(:)
TYPE(GAMLIBDATA), POINTER :: gisodata
! Local Variables
INTEGER :: id, igg, ib, ie, igg2, iso

! Allocate Macro data
IF(.NOT.XsMac%lalloc) THEN
  XsMac%ngg = ngg
  CALL AllocGamMacXs(XsMac)
END IF
! Pointing
XSmacsm => XsMac%XsMacSM
XSmacS => XsMac%XsMacS
XsMacStr => XsMac%XsMacStr
! Initialization
XsMacSm = 0._8
XsMacS = 0._8
XsMacStr = 0._8

DO iso = 1, niso
  id = mapnuclELM(idiso(iso))
  IF (id .EQ. 0) THEN
    PRINT *, 'ISOTOPE', idiso(iso), 'DOES NOT EXIST IN GAMMA LIBRARY(element)'
    CYCLE
  END IF
  gisodata => Gldiso(id)
  DO igg = igb, ige
    ib = gisodata%GSM(igg)%ib; ie = gisodata%GSM(igg)%ie
    DO igg2 = ib, ie
      XSMacSm(igg2,igg) = XsMacSm(igg2, igg) + pnum(iso) * gisodata%GSM(igg)%from(igg2) ! igg2 -> igg
    END DO
    XsMacS(igg) = XsMacS(igg) + pnum(iso) * gisodata%GSIGS(igg)
    XsMacSTR(igg) = XsMacSTR(igg) + pnum(iso) * gisodata%GSIGSTR(igg)
  END DO
  NULLIFY(gisodata)
END DO
IF(lScat1) THEN
  DO igg = igb, ige
    XsMacSm(igg, igg) = XsMacSm(igg, igg) + (XsMacS(igg) - XsMacStr(igg))
  ENDDO
END IF
! Nullifying Pointers
NULLIFY(XSmacsm, XSmacS, XsMacStr)
END SUBROUTINE


!***************************************************************************************************
! PN SCATTERING DATA
!
!   GamP1XsScatMatrix SUBROUTINES
!                     GENERATES P1 ORDER MACRO SCATTERING MATRIX
SUBROUTINE GamP1XsScatMatrix_Fxr(XsMac, FXR, igb, ige, ngg)
USE GammaTYPEDEF,   ONLY : GamMacXS_TYPE
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XsMac
TYPE(FXRInfo_TYPE) :: FXR
INTEGER :: igb, ige, ngg
CALL GamP1XsScatMatrix_Gen(XsMac, Fxr%niso, Fxr%idiso, Fxr%pnum, igb, ige, ngg)
END SUBROUTINE

SUBROUTINE GamP1XsScatMatrix_Gen(XsMac, niso, idiso, pnum, igb, ige, ngg)
USE GammaLibdata,   ONLY : Gldiso, mapnuclELM
USE GammaTYPEDEF,   ONLY : GamMacXS_TYPE, GAMLIBDATA
IMPLICIT NONE
! INPUT VARIABLES
TYPE(GamMacXS_TYPE) :: XsMac
INTEGER :: niso
INTEGER :: idiso(niso)
REAL :: pnum(niso)
INTEGER :: igb, ige, ngg
! POINTERS
REAL, POINTER :: XsMacP1Sm(:,:)
TYPE(GAMLIBDATA), POINTER :: gisodata
! LOCAL VARIABLES
INTEGER :: iso, id, igg, igg2, ib, ie

IF (.NOT. XsMac%lallocsm) THEN
  XsMac%ngg = ngg
  CALL AllocGamMacXs(XsMac)
END IF

XsMacP1Sm => XsMac%MacGSM1
XsMacP1Sm = 0._8

DO iso = 1, niso
  id = mapnuclELM(idiso(iso))
  IF (id .EQ. 0) CYCLE
  gisodata => Gldiso(id)
  DO igg = igb, ige
    ib = gisodata%GSM1(igg)%ib
    ie = gisodata%GSM1(igg)%ie
    DO igg2 = ib, ie
      XsMacP1Sm(igg2, igg) = XsMacP1Sm(igg2, igg) + pnum(iso) * gisodata%GSM1(igg)%from(igg2)
    END DO
  END DO
  NULLIFY(gisodata)
END DO

NULLIFY(XsMacP1Sm)
END SUBROUTINE

!   GamP2XsScatMatrix SUBROUTINES
!                     GENERATES P2 ORDER MACRO SCATTERING MATRIX
SUBROUTINE GamP2XsScatMatrix_Fxr(XsMac, FXR, igb, ige, ngg)
USE GammaTYPEDEF,   ONLY : GamMacXS_TYPE
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XsMac
TYPE(FXRInfo_TYPE) :: FXR
INTEGER :: igb, ige, ngg
CALL GamP2XsScatMatrix_Gen(XsMac, Fxr%niso, Fxr%idiso, Fxr%pnum, igb, ige, ngg)
END SUBROUTINE

SUBROUTINE GamP2XsScatMatrix_Gen(XsMac, niso, idiso, pnum, igb, ige, ngg)
USE GammaTYPEDEF,   ONLY : GamMacXS_TYPE, GAMLIBDATA
USE GammaLibdata,   ONLY : Gldiso, mapnuclELM
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XsMac
INTEGER :: niso
INTEGER :: idiso(niso)
REAL :: pnum(niso)
INTEGER :: igb, ige, ngg

REAL, POINTER :: XsMacP2Sm(:,:)
TYPE(GAMLIBDATA), POINTER :: gisodata

INTEGER :: iso, id, igg, igg2, ib, ie

IF (.NOT. XsMac%lallocsm) THEN
  XsMac%ngg = ngg
  CALL AllocGamMacXs(XsMac)
END IF

XsMacP2Sm => XsMac%MacGSM2
XsMacP2Sm = 0._8

DO iso = 1, niso
  id = mapnuclELM(idiso(iso))
  IF (id .EQ. 0) CYCLE
  gisodata => Gldiso(id)
  DO igg = igb, ige
    ib = gisodata%GSM2(igg)%ib
    ie = gisodata%GSM2(igg)%ie
    DO igg2 = ib, ie
      XsMacP2Sm(igg2, igg) = XsMacP2Sm(igg2, igg) + pnum(iso) * gisodata%GSM2(igg)%from(igg2)
    END DO
  END DO
  NULLIFY(gisodata)
END DO

NULLIFY(XsMacP2Sm)
END SUBROUTINE

!   GamP3XsScatMatrix SUBROUTINES
!                     GENERATES P3 ORDER MACRO SCATTERING MATRIX
SUBROUTINE GamP3XsScatMatrix_Fxr(XsMac, FXR, igb, ige, ngg)
USE GammaTYPEDEF,   ONLY : GamMacXS_TYPE
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XsMac
TYPE(FXRInfo_TYPE) :: FXR
INTEGER :: igb, ige, ngg
CALL GamP3XsScatMatrix_Gen(XsMac, Fxr%niso, Fxr%idiso, Fxr%pnum, igb, ige, ngg)
END SUBROUTINE

SUBROUTINE GamP3XsScatMatrix_Gen(XsMac, niso, idiso, pnum, igb, ige, ngg)
USE GammaTYPEDEF,   ONLY : GamMacXS_TYPE, GAMLIBDATA
USE GammaLibdata,   ONLY : Gldiso, mapnuclELM
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XsMac
INTEGER :: niso
INTEGER :: idiso(niso)
REAL :: pnum(niso)
INTEGER :: igb, ige, ngg

REAL, POINTER :: XsMacP3Sm(:,:)
TYPE(GAMLIBDATA), POINTER :: gisodata

INTEGER :: iso, id, igg, igg2, ib, ie

IF (.NOT. XsMac%lallocsm) THEN
  XsMac%ngg = ngg
  CALL AllocGamMacXs(XsMac)
END IF

XsMacP3Sm => XsMac%MacGSM3
XsMacP3Sm = 0._8

DO iso = 1, niso
  id = mapnuclELM(idiso(iso))
  IF (id .EQ. 0) CYCLE
  gisodata => Gldiso(id)
  DO igg = igb, ige
    ib = gisodata%GSM3(igg)%ib
    ie = gisodata%GSM3(igg)%ie
    DO igg2 = ib, ie
      XsMacP3Sm(igg2, igg) = XsMacP3Sm(igg2, igg) + pnum(iso) * gisodata%GSM3(igg)%from(igg2)
    END DO
  END DO
  NULLIFY(gisodata)
END DO

NULLIFY(XsMacP3Sm)
END SUBROUTINE


!***************************************************************************************************
! PHOTON PRODUCTION DATA
!
!   GamProdMatrix SUBROUTINE
!                 GenerateS Photon production matrix
SUBROUTINE GamProdMatrix_Fxr(XSMac, Fxr, igb, ige, ng, ngg, GroupInfo)
USE GamXSUtil,       ONLY :  GetGamXsMacDat,  ReturnGamXsMacDat,  IntProdMatCsp
USE GammaTYPEDEF,    ONLY :  GAMMAGROUPINFO_TYPE, GamMacXS_TYPE
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XSMac
TYPE(FxrInfo_TYPE) :: Fxr
INTEGER :: igb, ige, ng, ngg
TYPE(GAMMAGROUPINFO_TYPE) :: GroupInfo

CALL GamProdMatrix_Gen(XSMac, Fxr%temp, Fxr%niso, Fxr%idiso, Fxr%pnum,              &
  Fxr%resonfIso, Fxr%resocapIso, igb, ige, ng, ngg, GroupInfo, Fxr%lres)

END SUBROUTINE

SUBROUTINE GamProdMatrix_Csp(XSMac, Fxr, igb, ige, ng, ngg, GroupInfo, lCrCspFtn)
USE GamXSUtil,       ONLY :  GetGamXsMacDat,  ReturnGamXsMacDat,  IntProdMatCsp
USE GammaTYPEDEF,    ONLY :  GAMMAGROUPINFO_TYPE, GamMacXS_TYPE
IMPLICIT NONE
TYPE(GAMMAGROUPINFO_TYPE) :: GroupInfo
TYPE(FxrInfo_TYPE) :: Fxr
INTEGER :: ngg, ng, igb, ige
TYPE(GamMacXS_TYPE) :: XSMac
LOGICAL :: lCrCspFtn
LOGICAL :: lres

INTEGER :: niso
REAL :: temp
INTEGER, POINTER :: idiso(:)
REAL, POINTER :: pnum(:)
REAL, POINTER :: resonfIso(:,:), resocapIso(:,:)

TYPE(GamMacXS_TYPE), POINTER :: XsMac1, XsMac2

temp = Fxr%temp
resonfIso => Fxr%resonfIso
resocapIso => Fxr%resocapIso
lres = Fxr%lres
IF (lCrCspFtn.AND.Fxr%lCrCspFtn) THEN
  CALL GetGamXsMacDat(XsMac1, ngg, .FALSE.)
  CALL GetGamXsMacDat(XsMac2, ngg, .FALSE.)
  
  niso = Fxr%CspFxr%niso(1)
  Pnum => Fxr%CspFxr%pnum(:, 1);     idiso => Fxr%CspFxr%isolist(:, 1)
  CALL GamProdMatrix_Gen(XSMac, temp, niso, idiso, Pnum, resonfIso, resocapIso,                     &
               igb, ige, ng, ngg, GroupInfo, lres)
  niso = Fxr%CspFxr%niso(2)
  Pnum => Fxr%CspFxr%pnum(:, 2);     idiso => Fxr%CspFxr%isolist(:, 2)
  CALL GamProdMatrix_Gen(XSMac, temp, niso, idiso, Pnum, resonfIso, resocapIso,                     &
               igb, ige, ng, ngg, GroupInfo, lres)
  ! Interpolation of Two Csp
  CALL IntProdMatCsp(XsMac, XsMac1, XsMac2, Fxr%CspFxr, Fxr%niso, igb, ige, ng)
  CALL ReturnGamXsMacDat(XsMac1);        CALL ReturnGamXsMacDat(XsMac2)
  
ELSE
  niso = Fxr%niso
  Pnum => Fxr%pnum;     idiso => Fxr%idiso
  CALL GamProdMatrix_Gen(XSMac, temp, niso, idiso, Pnum, resonfIso, resocapIso,                     &
                            igb, ige, ng, ngg, GroupInfo, lres)
END IF

END SUBROUTINE
  
SUBROUTINE GamProdMatrix_Gen(XSMac, temp, niso, idiso, pnum, resonfIso, resocapIso,                 &
                    igb, ige, ng, ngg, GroupInfo, lres)
USE GamXSUtil,       ONLY : AllocGamProdMat
USE GammaLibdata,    ONLY : neltGAM, Gldiso, imtRAD, imtINEL, imtFIS, mapnuclGAM
USE GammaTYPEDEF,    ONLY : GAMMAGROUPINFO_TYPE, GamMacXS_TYPE, GAMLIBDATA
IMPLICIT NONE
! INPUT VARIABLES
TYPE(GAMMAGROUPINFO_TYPE) :: GroupInfo
INTEGER :: ngg, ng, igb, ige
TYPE(GamMacXS_TYPE) :: XSMac
INTEGER :: niso
INTEGER :: idiso(niso)
REAL :: pnum(niso)
REAL :: temp
REAL :: resonfIso(niso,ng), resocapIso(niso,ng)
LOGICAL :: lres

! POINTERS
TYPE(GAMLIBDATA), POINTER :: gisodata
REAL, POINTER, DIMENSION(:,:,:) :: IsoProdFis, IsoProdRad, IsoProdInel
REAL, POINTER, DIMENSION(:,:)  :: MatTot
LOGICAL, POINTER, DIMENSION(:) :: lfis, lrad, linel
INTEGER, POINTER, DIMENSION(:) :: ifisb, ifise, iradb, irade, iinelb, iinele

! LOCAL VARIABLES
INTEGER :: iso, id, ig, igg, imt, nx1, nx2
REAL :: wt1, wt2
INTEGER :: it1, it2
INTEGER :: nofg, norg
LOGICAL, allocatable :: lresogrp(:)
INTEGER :: tempgrpb, tempgrpe

! ALLOCATION FOR PRODUCTION MATRIX
IF(.NOT. XSMac%lProdAlloc) THEN
  XSMac%ng = ng
  XSMac%ngg = ngg
  XSMac%niso = neltGAM
  CALL AllocGamProdMat(XSMac)
ENDIF
XSMac%GProdTot = 0.
! DISTINGUISH RESONANCE GROUP
ALLOCATE(lresogrp(ng))
nofg = GroupInfo%nofg; norg = GroupInfo%norg
lresogrp = .FALSE.
IF (lres) THEN
  DO ig = 1, ng
    if(ig.gt.nofg.and.ig.le.(nofg+norg)) lresogrp(ig) = .TRUE.
  END DO
END IF
! INITIALIZATION
XSMac%IsoGProdFis = 0._8
XSMac%IsoGProdRad = 0._8
XSMac%IsoGProdInel = 0._8
XSMac%GProdTot = 0._8
XSMac%ifisb = 0; XSMac%ifise = 0; XSMac%iradb = 0; XSMac%irade = 0
XSMac%iinelb = 0; XSMac%iinele = 0
XSMac%lfis = .FALSE.; XSMac%lrad = .FALSE.; XSMac%linel = .FALSE.
! POINTING
IsoProdFis => XSMac%IsoGProdFis; IsoProdRad => XSMac%IsoGProdRad; IsoProdInel => XSMac%IsoGProdInel
MatTot => XSMac%GProdTot
lfis => XSMac%lfis; lrad => XSMac%lrad; linel => XSMac%linel
ifisb => XSMac%ifisb; iradb => XSMac%iradb; iinelb => XSMac%iinelb;
ifise => XSMac%ifise; irade => XSMac%irade; iinele => XSMac%iinele
! ISOTOPEWISE (TEMP. DEPENDENCE + PNUM)
DO iso = 1, niso
  id = mapnuclGAM(idiso(iso))
  IF (id .EQ. 0) THEN
    print *, 'SUBROUTINE GamProdMatrix -- No photon production data in library for', idiso(iso)
    CYCLE
  END IF
  gisodata => Gldiso(id)
  IF (gisodata%ityp .LT. 1) STOP 'SUBROUTINE GamProdMatrix -- GAMMA LIBRARY ERROR! TYPE INCORRECT'
  CALL GamXsTempInterpolation(id, Gisodata, temp, wt1, wt2, it1, it2) ! TEMPERATURE DEPENDENCE
  IF (gisodata%lphoton(imtFIS)) THEN  ! FISSION
    nx1 = gisodata%ppgrlow(imtFIS); nx2 = gisodata%ppgrup(imtFIS)
    ifisb(iso) = nx1; ifise(iso) = nx2
    IF(nx2.GE.igb .AND. nx1.LE.ige) THEN    ! Data existing region is beyond interest
      lfis(iso) = .TRUE.
      tempgrpb = MAX(nx1, igb)
      tempgrpe = MIN(nx2, ige)
      DO igg = tempgrpb, tempgrpe
        nx1 = gisodata%GPM(igg,it1,imtFIS)%ib; nx2 = gisodata%GPM(igg,it1,imtFIS)%ie
        DO ig = nx1, nx2
          IsoProdFis(iso,ig,igg) = pnum(iso) * (wt2 * gisodata%GPM(igg,it2,imtFIS)%from(ig)          &
                              + wt1 * gisodata%GPM(igg,it1,imtFIS)%from(ig))
        END DO
      END DO
    END IF
  END IF
  IF (gisodata%lphoton(imtRAD)) THEN  ! RADIOACTIVE CAPTURE
    nx1 = gisodata%ppgrlow(imtRAD); nx2 = gisodata%ppgrup(imtRAD)
    iradb(iso) = nx1; irade(iso) = nx2
    IF(nx2.GE.igb .AND. nx1.LE.ige) THEN    ! Data existing region is beyond interest
      lrad(iso) = .TRUE.
      tempgrpb = MAX(nx1, igb)
      tempgrpe = MIN(nx2, ige)
      DO igg = tempgrpb, tempgrpe
        nx1 = gisodata%GPM(igg,it1,imtRAD)%ib; nx2 = gisodata%GPM(igg,it1,imtRAD)%ie
        DO ig = nx1, nx2
          IsoProdRad(iso,ig,igg) = pnum(iso) * (wt2 * gisodata%GPM(igg,it2,imtRAD)%from(ig)          &
                              + wt1 * gisodata%GPM(igg,it1,imtRAD)%from(ig))
        END DO
      END DO
    END IF
  END IF
  IF (gisodata%lphoton(imtINEL)) THEN ! INELASTIC SCATTERING (no need for isotope wise data)
    nx1 = gisodata%ppgrlow(imtINEL); nx2 = gisodata%ppgrup(imtINEL)
    iinelb(iso) = nx1; iinele(iso) = nx2
    IF(nx2.GE.igb .AND. nx1.LE.ige) THEN    ! Data existing region is beyond interest
      linel(iso) = .TRUE.
      tempgrpb = MAX(nx1, igb)
      tempgrpe = MIN(nx2, ige)
      DO igg = tempgrpb, tempgrpe
        nx1 = gisodata%GPM(igg,it1,imtINEL)%ib; nx2 = gisodata%GPM(igg,it1,imtINEL)%ie
        DO ig = nx1, nx2
          IsoProdInel(iso,ig,igg) = pnum(iso) * (wt2 * gisodata%GPM(igg,it2,imtINEL)%from(ig)         &
                              + wt1 * gisodata%GPM(igg,it1,imtINEL)%from(ig))
        END DO
      END DO
    END IF
  END IF
  NULLIFY(gisodata)
END DO

IF (.NOT.any(lfis.OR.lrad.OR.linel)) RETURN

DO igg = igb, ige
  DO ig = 1, ng
    IF (lresogrp(ig)) THEN
      DO iso = 1, niso
        IF(lfis(iso)) MatTot(ig, igg) = MatTot(ig, igg) + IsoProdFis(iso,ig,igg) * resonfIso(iso, ig)
        IF(lrad(iso)) MatTot(ig, igg) = MatTot(ig, igg) + IsoProdRad(iso,ig,igg) * resocapIso(iso, ig)
        IF(linel(iso)) MatTot(ig, igg) = MatTot(ig, igg) + IsoProdInel(iso,ig,igg)
      END DO
    ELSE
      DO iso = 1, niso
        IF(lfis(iso)) MatTot(ig, igg) = MatTot(ig, igg) + IsoProdFis(iso,ig,igg)
        IF(lrad(iso)) MatTot(ig, igg) = MatTot(ig, igg) + IsoProdRad(iso,ig,igg)
        IF(linel(iso)) MatTot(ig, igg) = MatTot(ig, igg) + IsoProdInel(iso,ig,igg)
      END DO
    END IF
  END DO
END DO

DEALLOCATE(lresogrp)
NULLIFY(IsoProdFis, IsoProdRad, IsoProdInel, MatTot)
NULLIFY(lfis, lrad, linel)
NULLIFY(ifisb, iradb, iinelb, ifise, irade, iinele)
END SUBROUTINE


!***************************************************************************************************
! LOCAL ENERGY DEPOSIT DATA DUE TO NEUTRON                                |-- JSU EDIT 2017.09.13. |
!--  FOR EXPLICIT KAPPA CALCULATION
!
!   GetLocalQMat  SUBROUTINE
!                 GenerateS Locally Deposited Energy Matrix due to Neutron Reaction
!                   1. LOCALLY DEPOSITED FISSION ENERGY
!                   2. N2N, N3N REACTION ENERGY LOSS (TOTAL RECOVERABLE ENERGY DECREASE
!                   3.a. INELASTIC SCATTERING LOCAL DEPOSIT ENERGY (NEUTRON ENERGY LOSS MATRIX)
SUBROUTINE GetLocalQMat_FXR(XSMac, Fxr, igb, ige, ng, ngg, lfis)
USE TYPEDEF,         ONLY : FXRINFO_TYPE
USE GammaTYPEDEF,    ONLY : GamMacXS_TYPE
USE CORE_MOD,        ONLY : GroupInfo
IMPLICIT NONE
TYPE(FXRINFO_TYPE) :: FXR
TYPE(GamMacXS_TYPE) :: XSMac
INTEGER :: igb, ige, ng, ngg
LOGICAL :: lfis
INTEGER :: iResoGrpBeg, iResoGrpEnd, norg

iResoGrpBeg = GroupInfo%nofg + 1
iResoGrpEnd = GroupInfo%nofg + GroupInfo%norg
norg = GroupInfo%norg

CALL GetLocalQMat_GEN(XSMac, FXR%temp, FXR%niso, FXR%idiso, FXR%pnum, igb, ige, ng, ngg, FXR%fresonf, iResoGrpBeg, iResoGrpEnd, lfis)
                                                                                                      
END SUBROUTINE

!SUBROUTINE GetLocalQMat_CSP(XSMac, Fxr, igb, ige, ng, ngg, GroupInfo, lCrCspFtn)
!IMPLICIT NONE
!END SUBROUTINE
!
SUBROUTINE GetLocalQMat_GEN(XSMac, temp, niso, idiso, pnum, igb, ige, ng, ngg, fresonf, iResoGrpBeg, iResoGrpEnd, lfis)
! NEUTRON INFORMATION
USE XsTypeDef,       ONLY : LIBDATA
USE XSLIB_MOD,       ONLY : ldiso, nelthel, mapnucl, mapfis, itempmap, nelrhel
! GAMMA INFORMATION
USE GamXSUtil,       ONLY : AllocLocalKappa
USE GammaLibdata,    ONLY : neltGAM, Gldiso, imtINEL, mapnuclGAM
USE GammaTYPEDEF,    ONLY : GamMacXS_TYPE, GAMLIBDATA
USE GammaCore_mod,   ONLY : GamGroupInfo
USE CORE_MOD,        ONLY : GroupInfo

IMPLICIT NONE

! INPUT VARIABLES
INTEGER :: ngg, ng, igb, ige
TYPE(GamMacXS_TYPE) :: XSMac
INTEGER :: niso
INTEGER :: idiso(niso)
REAL :: pnum(niso)
REAL :: temp
INTEGER :: iResoGrpBeg, iResoGrpEnd, norg
REAL :: fresonf(iResoGrpBeg:iResoGrpEnd)
LOGICAL :: lfis

! POINTING VARIABLES
TYPE(LIBDATA), POINTER :: isodata
TYPE(GAMLIBDATA), POINTER :: gisodata
REAL, POINTER, DIMENSION(:) :: FisLocal, QN2N, QN3N, InelLocal
REAL, POINTER, DIMENSION(:, :) :: InelLoss, GProdInel
REAL, POINTER, DIMENSION(:) :: GamAvgE, NeuAvgE

! LOCAL VARIABLES
INTEGER :: nid, gid, inmn
INTEGER :: iso, ig, ig2, igg
INTEGER :: it1, it2
INTEGER :: ScRange1(2), ScRange2(2)
INTEGER :: ioutbeg, ioutend
REAL :: wt1, wt2
REAL :: isokappa, ison2n, ison3n
REAL :: sigs1, sigs2
INTEGER :: ind1, ind2
!
!iResoGrpBeg = GroupInfo%nofg + 1
!iResoGrpEnd = GroupInfo%nofg + GroupInfo%norg
!norg = GroupInfo%norg
iResoGrpBeg = MAX(iResoGrpBeg, igb)
iResoGrpEnd = MIN(iResoGrpEnd, ige)

! ALLOCATION
IF (.NOT. XsMac%lKappaAlloc) THEN
  XsMac%ng = ng
  XsMac%ngg = ngg
  XsMac%niso = nelthel
  CALL AllocLocalKappa(XsMac)
END IF
! INITIALIZATION
XsMac%FisLocal = 0.  
XsMac%QN2N = 0.
XsMac%QN3N = 0.
! XsMac%InelLosS = 0.
! XsMac%GProdInel = 0.
XsMac%InelLocal = 0.
XsMac%LocalQ = 0.

! POINTING 
FisLocal => XsMac%FisLocal
QN2N => XsMac%QN2N
QN3N => XsMac%QN3N
! InelLoss => XsMac%InelLoss
! GProdInel => XsMac%GProdInel
InelLocal => XsMac%InelLocal

GamAvgE => GamGroupInfo%GamAvgE
NeuAvgE => GamGroupInfo%NeuAvgE


! LOCAL DEPOSIT FISSION ENERGY and N2N, N3N REACTION ENERGY
DO iso = 1, niso
  isokappa = 0. ;ison2n = 0. ;ison3n = 0.
  nid = MapNucl(idiso(iso)); isodata => ldiso(nid)
  IF(nid .gt. nelrhel) CYCLE
  inmn = isodata%inmn
  gid = MapNuclGam(idiso(iso))
  IF (gid .NE. 0) THEN
    gisodata => Gldiso(gid)
    ! NMN REACTION
    IF (inmn .EQ. 1 .OR. inmn .EQ. 3) THEN ! N2N REACTION
      ison2n = gisodata%n2nQ
      DO ig = igb, ige
        QN2N(ig) = QN2N(ig) + pnum(iso) * isodata%sign2n(ig) * ison2n
      END DO
    END IF
    IF (inmn .EQ. 2 .OR. inmn .EQ. 3) THEN ! N3N REACTION
      ison3n = gisodata%n3nQ
      DO ig = igb, ige
        QN3N(ig) = QN3N(ig) + pnum(iso) * isodata%sign3n(ig) * ison3n
      END DO
    END IF
    ! INELASTIC SCATTERING
    DO ig = igb, ige
      ! ENERGY LOSS DUE TO INELASTIC SCATTERING
      CALL XsTempInterpolation(nid, isodata, temp, wt1, wt2, it1, it2)
      ScRange1 = (/isodata%sm(ig, it1)%ioutsb, isodata%sm(ig, it1)%ioutse/)
      ScRange2 = (/isodata%sm(ig, it2)%ioutsb, isodata%sm(ig, it2)%ioutse/)  
      ioutbeg = MIN(ScRange1(1), ScRange2(1)); ioutend = MAX(ScRange1(2), ScRange2(2))
      DO ig2 = ioutbeg, ioutend ! ig -> ig2
        sigs1 = 0.; sigs2 = 0.
        ind1 = (ig2 - ScRange1(1)) * (ig2 - ScRange1(2))
        ind2 = (ig - IsoData%sm(ig2, it1)%ib) * (ig - IsoData%sm(ig2, it1)%ie)
        IF(ind1 .le. 0 .and. ind2 .le. 0) sigs1 = isodata%sm(ig2, it1)%from(ig)
        ind1 = (ig2 - ScRange2(1)) * (ig2 - ScRange2(2))
        ind2 = (ig - IsoData%sm(ig2, it2)%ib) * (ig - IsoData%sm(ig2, it2)%ie)
        IF(ind1 .le. 0 .and. ind2 .le. 0) sigs2 = isodata%sm(ig2, it2)%from(ig)
        InelLocal(ig) = InelLocal(ig) + pnum(iso) * (wt1 * sigs1 + wt2 * sigs2) * (NeuAvgE(ig) - NeuAvgE(ig2))
      END DO
      
      ! ENERGY CONVERTED AS PHOTON
      CALL GamXsTempInterpolation(gid, Gisodata, temp, wt1, wt2, it1, it2)
      ScRange1 = (/gisodata%GPMOUTB(ig, it1, imtInel), gisodata%GPMOUTE(ig, it1, imtInel)/)
      ScRange2 = (/gisodata%GPMOUTB(ig, it2, imtInel), gisodata%GPMOUTE(ig, it2, imtInel)/)
      ioutbeg = MIN(ScRange1(1), ScRange2(1)); ioutend = MAX(ScRange1(2), ScRange2(2))
      DO igg = ioutbeg, ioutend
        sigs1 = 0.; sigs2 = 0.
        ind1 = (igg - ScRange1(1)) * (igg - ScRange1(2))
        ind2 = (ig - gisodata%GPM(igg, it1, imtInel)%ib) * (ig - gisodata%GPM(igg, it1, imtInel)%ie)
        IF(ind1 .le. 0 .and. ind2 .le. 0) sigs1 = gisodata%GPM(igg, it1, imtInel)%from(ig)
        ind1 = (igg - ScRange2(1)) * (igg - ScRange2(2))
        ind2 = (ig - gisodata%GPM(igg, it1, imtInel)%ib) * (ig - gisodata%GPM(igg, it1, imtInel)%ie)
        IF(ind1 .le. 0 .and. ind2 .le. 0) sigs2 = gisodata%GPM(igg, it2, imtInel)%from(ig)
        InelLocal(ig) = InelLocal(ig) - pnum(iso) * (wt1 * sigs1 + wt2 * sigs2) * GamAvgE(igg)
      END DO
    END DO
  END IF
  
  ! FISSION REACTION
  IF(mapfis(idiso(iso)) .NE. 0) THEN
    CALL XsTempInterpolation(nid, isodata, temp, wt1, wt2, it1, it2)
    IF (gid .NE. 0) THEN
      isokappa = gisodata%EXKAPPA(1) + gisodata%EXKAPPA(6)  ! FISSION PRODUCTS AND DELAYED BETA 
    ELSE
     isokappa = isodata%kappa
    END IF
    
    DO ig = igb, ige
      FisLocal(ig) = FisLocal(ig) + pnum(iso) * (wt1 * isodata%sigf(ig, it1) + wt2 * isodata%sigf(ig, it2)) * isokappa
    END DO
  END IF
END DO ! ISO LOOP

IF (lfis) THEN ! RESONANCE TREATMENT FOR FISSION
 FisLocal(iResoGrpBeg:iResoGrpEnd) = FisLocal(iResoGrpBeg:iResoGrpEnd) * fresonf(iResoGrpBeg:iResoGrpEnd)
END IF

XsMac%LocalQ = FisLocal + QN2N + QN3N + InelLocal

END SUBROUTINE


End Module
#endif