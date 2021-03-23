#include <defines.h>
#ifdef __GAMMA_TRANSPORT    
Module GamXsLib_Mod
!-- JSU EDIT
! MODULE THAT CONTAINS SUBROUTINES USED IN MOC CROSS SECTION AND SOURCE GENERATION
USE PARAM
USE GamXSUtil,       ONLY : AllocGamMacXs,  AllocGamMacIsoXs,  AllocGamProdMat,   & 
                            AllocGamIsoSM,  FreeGamSMISO
!                            GamXsTempInterpolation
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
USE TYPEDEF,         ONLY : Fxrinfo_type  
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
USE TYPEDEF,        ONLY : Fxrinfo_type
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
USE GammaTYPEDEF,   ONLY : GamMacXS_TYPE
USE XSLIB_MOD,      ONLY : GAMLIBDATA, nelmGAM, phatom, mapnucl, ldiso
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
IF ((.NOT.XsMac%lIsoAlloc).AND.lIsoXsOut) THEN
  XsMac%ngg = ngg
  XsMac%niso = nelmGAM
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
  id = mapnucl(idiso(iso))
  IF (id .EQ. 0) THEN
    PRINT *, 'ISOTOPE', idiso(iso), 'DOES NOT EXIST IN GAMMA LIBRARY(element)'
    CYCLE
  END IF
  Gisodata => ldiso(id)%phatom
  DO igg = igb, ige
    IsoKERMA(iso, igg) = pnum(iso) * Gisodata%KERMA(igg)
    IsoXsMacTR(iso, igg) = pnum(iso) * Gisodata%SIGTR(igg)
    IsoXsMacA(iso, igg) = pnum(iso) * Gisodata%SIGA(igg)
    IsoXsMacS(iso, igg) = pnum(iso) * Gisodata%SIGS(igg)
    IsoXsMacT(iso, igg) = IsoXsMacA(iso, igg) + IsoXsMacS(iso, igg)
    IsoXsMacSTR(iso, igg) = IsoXsMacTR(iso, igg) - IsoXsMacA(iso, igg)
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
  XsMac%IsoXsMacSTR(1:niso, igb:ige) = IsoXsMacSTR(1:niso, igb:ige)
END IF
END SUBROUTINE

!   GamScatMatrix SUBROUTINES
!                 GENERATES MACROSCOPIC 0-TH ORDER SCATTERING MATRIX
SUBROUTINE GamScatMatrix_Fxr(XsMac, FXR, igb, ige, ngg, lscat1,lISOscat)
USE GammaTYPEDEF,   ONLY : GamMacXS_TYPE
USE TYPEDEF,        ONLY : Fxrinfo_type
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XsMac
TYPE(FXRInfo_TYPE) :: FXR
INTEGER :: igb, ige, ngg
LOGICAL :: lscat1,lISOscat
CALL GamScatMatrix_Gen(XsMac, FXR%niso, FXR%idiso, FXR%pnum, igb, ige, ngg, lscat1, lISOscat)
END SUBROUTINE

SUBROUTINE GamScatMatrix_Csp(XsMac, FXR, igb, ige, ngg, lscat1, lCrCspFtn,lISOscat)
USE GamXSUtil,      ONLY : GetGamXsMacDat,  ReturnGamXsMacDat,  IntScatMatCsp
USE GammaTYPEDEF,   ONLY : GamMacXS_TYPE
USE TYPEDEF,         ONLY : Fxrinfo_type  
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XsMac
TYPE(FXRInfo_TYPE) :: FXR
INTEGER :: igb, ige, ngg
LOGICAL :: lscat1, lCrCspFtn, lISOscat

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
  CALL GamScatMatrix_Gen(XsMac1, niso, idiso, pnum, igb, ige, ngg, lscat1,lISOscat)
  
  niso = Fxr%CspFxr%niso(2)
  Pnum => Fxr%CspFxr%pnum(:, 2);     idiso => Fxr%CspFxr%isolist(:, 2)
  CALL GamScatMatrix_Gen(XsMac2, niso, idiso, pnum, igb, ige, ngg, lscat1,lISOscat)
  
  ! Interpolation of Two Csp
  CALL IntScatMatCsp(XsMac, XsMac1, XsMac2, Fxr%CspFxr, Fxr%niso, igb, ige, ngg)
  CALL ReturnGamXsMacDat(XsMac1);        CALL ReturnGamXsMacDat(XsMac2)
  
ELSE
  niso = Fxr%niso;                  Temp = Fxr%Temp
  pnum => Fxr%pnum;                 idiso => Fxr%idiso
  CALL GamScatMatrix_Gen(XsMac, niso, idiso, pnum, igb, ige, ngg, lscat1,lISOscat)
ENDIF
END SUBROUTINE

SUBROUTINE GamScatMatrix_Gen(XsMac, niso, idiso, pnum, igb, ige, ngg, lscat1, lISOscat)
USE XSLIB_MOD,   ONLY : phatom, GAMLIBDATA, mapnucl, ldiso
USE GammaTYPEDEF,   ONLY : GamMacXS_TYPE
IMPLICIT NONE
! Input Variables
TYPE(GamMacXS_TYPE) :: XsMac
INTEGER :: niso
INTEGER :: idiso(niso)
REAL :: pnum(niso)
INTEGER :: igb, ige, ngg
LOGICAL :: lscat1, lISOscat
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

IF (lISOscat) THEN
  CALL FreeGamSMISO(XsMac)
  XsMac%ngg = ngg
  XsMac%niso = niso
  CALL AllocGamIsoSM(XsMac)
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
  id = mapnucl(idiso(iso))
  IF (id .EQ. 0) THEN
    PRINT *, 'ELEMENT', idiso(iso), 'DOES NOT EXIST IN GAMMA LIBRARY(element) ** SCATMAT GEN'
    CYCLE
  END IF
  gisodata => ldiso(id)%phatom
  DO igg = igb, ige
    ib = gisodata%SM(igg)%ib; ie = gisodata%SM(igg)%ie
    DO igg2 = ib, ie
      XSMacSm(igg2,igg) = XsMacSm(igg2, igg) + pnum(iso) * gisodata%SM(igg)%from(igg2) ! igg2 -> igg
    END DO
    XsMacS(igg) = XsMacS(igg) + pnum(iso) * gisodata%SIGS(igg)
    XsMacSTR(igg) = XsMacSTR(igg) + pnum(iso) * gisodata%SIGSTR(igg)
  END DO
  NULLIFY(gisodata)
END DO
IF (lISOscat) THEN
  XsMac%IsoSM = 0.
  IF (lScat1) THEN
    DO iso = 1, niso
      ID = MAPNUCL(IDISO(ISO))
      if (ID.EQ.0) THEN
        CYCLE
      END IF
      gisodata => ldiso(id)%phatom
      DO igg = igb, ige
        ib = gisodata%SM(igg)%ib; ie = gisodata%SM(igg)%ie
        DO igg2 = ib, ie
          XsMac%IsoSM(igg2,igg,iso) = pnum(iso) * gisodata%SM(igg)%from(igg2) ! igg2 -> igg
        END DO
        XsMac%IsoSM(igg, igg, iso) =pnum(iso) * (gisodata%SM(igg)%from(igg) + gisodata%SIGS(igg) - gisodata%SIGSTR(igg))
      END DO
    END DO
  ELSE
    DO iso = 1, niso
      ID = MAPNUCL(IDISO(ISO))
      if (ID.EQ.0) THEN
        CYCLE
      END IF
      gisodata => ldiso(id)%phatom
      DO igg = igb, ige
        ib = gisodata%SM(igg)%ib; ie = gisodata%SM(igg)%ie
        DO igg2 = ib, ie
          XsMac%IsoSM(igg2,igg,iso) = pnum(iso) * gisodata%SM(igg)%from(igg2) ! igg2 -> igg
        END DO
      END DO
    END DO
  END IF
  NULLIFY(gisodata)
END IF

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
USE TYPEDEF,        ONLY : Fxrinfo_type
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XsMac
TYPE(FXRInfo_TYPE) :: FXR
INTEGER :: igb, ige, ngg
CALL GamP1XsScatMatrix_Gen(XsMac, Fxr%niso, Fxr%idiso, Fxr%pnum, igb, ige, ngg)
END SUBROUTINE

SUBROUTINE GamP1XsScatMatrix_Gen(XsMac, niso, idiso, pnum, igb, ige, ngg)
USE XSLIB_MOD,   ONLY : phatom, mapEleGAM, GAMLIBDATA
USE GammaTYPEDEF,   ONLY : GamMacXS_TYPE
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
  id = mapEleGAM(idiso(iso))
  IF (id .EQ. 0) CYCLE
  gisodata => phatom(id)
  DO igg = igb, ige
    ib = gisodata%SMp1(igg)%ib
    ie = gisodata%SMp1(igg)%ie
    DO igg2 = ib, ie
      XsMacP1Sm(igg2, igg) = XsMacP1Sm(igg2, igg) + pnum(iso) * gisodata%SMp1(igg)%from(igg2)
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
USE TYPEDEF,        ONLY : Fxrinfo_type
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XsMac
TYPE(FXRInfo_TYPE) :: FXR
INTEGER :: igb, ige, ngg
CALL GamP2XsScatMatrix_Gen(XsMac, Fxr%niso, Fxr%idiso, Fxr%pnum, igb, ige, ngg)
END SUBROUTINE

SUBROUTINE GamP2XsScatMatrix_Gen(XsMac, niso, idiso, pnum, igb, ige, ngg)
USE GammaTYPEDEF,   ONLY : GamMacXS_TYPE
USE XSLIB_MOD,      ONLY : phatom, mapEleGAM, GAMLIBDATA
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
  id = mapEleGAM(idiso(iso))
  IF (id .EQ. 0) CYCLE
  gisodata => phatom(id)
  DO igg = igb, ige
    ib = gisodata%SMp2(igg)%ib
    ie = gisodata%SMp2(igg)%ie
    DO igg2 = ib, ie
      XsMacP2Sm(igg2, igg) = XsMacP2Sm(igg2, igg) + pnum(iso) * gisodata%SMp2(igg)%from(igg2)
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
USE TYPEDEF,         ONLY : Fxrinfo_type  
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XsMac
TYPE(FXRInfo_TYPE) :: FXR
INTEGER :: igb, ige, ngg
CALL GamP3XsScatMatrix_Gen(XsMac, Fxr%niso, Fxr%idiso, Fxr%pnum, igb, ige, ngg)
END SUBROUTINE

SUBROUTINE GamP3XsScatMatrix_Gen(XsMac, niso, idiso, pnum, igb, ige, ngg)
USE GammaTYPEDEF,   ONLY : GamMacXS_TYPE
USE XSLIB_MOD,      ONLY : phatom, mapEleGAM, GAMLIBDATA
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
  id = mapEleGAM(idiso(iso))
  IF (id .EQ. 0) CYCLE
  gisodata => phatom(id)
  DO igg = igb, ige
    ib = gisodata%SMp3(igg)%ib
    ie = gisodata%SMp3(igg)%ie
    DO igg2 = ib, ie
      XsMacP3Sm(igg2, igg) = XsMacP3Sm(igg2, igg) + pnum(iso) * gisodata%SMp3(igg)%from(igg2)
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
SUBROUTINE GamProdMatrix_Fxr(XSMac, Fxr, igb, ige, ng, ngg, GroupInfo, lIsoOut)
USE GamXSUtil,       ONLY : GetGamXsMacDat,  ReturnGamXsMacDat,  IntProdMatCsp
USE GammaTYPEDEF,    ONLY : GamMacXS_TYPE
USE TYPEDEF,         ONLY : GROUPINFO_TYPE, Fxrinfo_type
IMPLICIT NONE
TYPE(GamMacXS_TYPE) :: XSMac
TYPE(FxrInfo_TYPE) :: Fxr
INTEGER :: igb, ige, ng, ngg
TYPE(GROUPINFO_TYPE) :: GroupInfo
LOGICAL :: lIsoOut

CALL GamProdMatrix_Gen(XSMac, Fxr%temp, Fxr%niso, Fxr%idiso, Fxr%pnum,              &
  Fxr%FresofIso, Fxr%fresocapIso, igb, ige, ng, ngg, GroupInfo, Fxr%lres, lIsoOut)

END SUBROUTINE

SUBROUTINE GamProdMatrix_Csp(XSMac, Fxr, igb, ige, ng, ngg, GroupInfo, lCrCspFtn, lIsoOut)
USE GamXSUtil,       ONLY :  GetGamXsMacDat,  ReturnGamXsMacDat,  IntProdMatCsp
USE GammaTYPEDEF,    ONLY :  GamMacXS_TYPE
USE TYPEDEF,         ONLY : GROUPINFO_TYPE,  FxrInfo_TYPE
IMPLICIT NONE
TYPE(GROUPINFO_TYPE) :: GroupInfo
TYPE(FxrInfo_TYPE) :: Fxr
INTEGER :: ngg, ng, igb, ige
TYPE(GamMacXS_TYPE) :: XSMac
LOGICAL :: lCrCspFtn, lIsoOut
LOGICAL :: lres

INTEGER :: niso
REAL :: temp
INTEGER, POINTER :: idiso(:)
REAL, POINTER :: pnum(:)
REAL(4), POINTER :: FresofIso(:,:), fresocapIso(:,:)

TYPE(GamMacXS_TYPE), POINTER :: XsMac1, XsMac2

temp = Fxr%temp
FresofIso => Fxr%FresofIso
fresocapIso => Fxr%fresocapIso
lres = Fxr%lres
IF (lCrCspFtn.AND.Fxr%lCrCspFtn) THEN
  CALL GetGamXsMacDat(XsMac1, ngg, .FALSE.)
  CALL GetGamXsMacDat(XsMac2, ngg, .FALSE.)
  
  niso = Fxr%CspFxr%niso(1)
  Pnum => Fxr%CspFxr%pnum(:, 1);     idiso => Fxr%CspFxr%isolist(:, 1)
  CALL GamProdMatrix_Gen(XSMac, temp, niso, idiso, Pnum, FresofIso, fresocapIso,                     &
               igb, ige, ng, ngg, GroupInfo, lres, lIsoOut)
  niso = Fxr%CspFxr%niso(2)
  Pnum => Fxr%CspFxr%pnum(:, 2);     idiso => Fxr%CspFxr%isolist(:, 2)
  CALL GamProdMatrix_Gen(XSMac, temp, niso, idiso, Pnum, FresofIso, fresocapIso,                     &
               igb, ige, ng, ngg, GroupInfo, lres, lIsoOut)
  ! Interpolation of Two Csp
  CALL IntProdMatCsp(XsMac, XsMac1, XsMac2, Fxr%CspFxr, Fxr%niso, igb, ige, ng)
  CALL ReturnGamXsMacDat(XsMac1);        CALL ReturnGamXsMacDat(XsMac2)
  
ELSE
  niso = Fxr%niso
  Pnum => Fxr%pnum;     idiso => Fxr%idiso
  CALL GamProdMatrix_Gen(XSMac, temp, niso, idiso, Pnum, FresofIso, fresocapIso,                     &
                            igb, ige, ng, ngg, GroupInfo, lres, lIsoOut)
END IF

END SUBROUTINE
  
SUBROUTINE GamProdMatrix_Gen(XSMac, temp, niso, idiso, pnum, FresofIso, fresocapIso,                 &
                    igb, ige, ng, ngg, GroupInfo, lres, lIsoOut)
USE GamXSUtil,       ONLY : AllocGamProdMat, AllocIsoGamPRodMat, FreeGamProdMat, FreeIsoGamPRodMat
USE XSLIB_mod,       ONLY : nelthel, ldiso, mapnucl, libdata
USE GammaTYPEDEF,    ONLY : GamMacXS_TYPE
USE TYPEDEF,         ONLY : GROUPINFO_TYPE
IMPLICIT NONE
! INPUT VARIABLES
TYPE(GROUPINFO_TYPE) :: GroupInfo
INTEGER :: ngg, ng, igb, ige
TYPE(GamMacXS_TYPE) :: XSMac
INTEGER :: niso
INTEGER :: idiso(niso)
REAL :: pnum(niso)
REAL :: temp
REAL(4),POINTER :: FresofIso(:,:), fresocapIso(:,:)
REAL(4), POINTER :: FresoReact(:, :)
LOGICAL :: lres, lIsoOut

! POINTERS
TYPE(libdata), POINTER :: lib
REAL, POINTER, DIMENSION(:,:,:) :: IsoProdFis, IsoProdRad, IsoProdInel, IsoProdNnel
REAL, POINTER, DIMENSION(:,:)  :: MatTot
LOGICAL, POINTER, DIMENSION(:,:) :: lexist
!INTEGER, POINTER, DIMENSION(:) :: ifisb, ifise, iradb, irade, iinelb, iinele

! LOCAL VARIABLES
INTEGER :: iso, id, ig, igg, imt, nx1, nx2
REAL :: wt1, wt2
INTEGER :: it1, it2
INTEGER :: nofg, norg
INTEGER :: ProdRange1(2), ProdRange2(2) !Produciton Range
INTEGER :: ind
REAL :: prod1, prod2
LOGICAL, allocatable :: lresogrp(:)
INTEGER :: tempgrpb, tempgrpe
REAL, POINTER :: IsoProd(:, :, :)

! ALLOCATION FOR PRODUCTION MATRIX
IF(.NOT. XSMac%lProdAlloc .OR. niso .NE. XSMac%niso) THEN
  IF (XSMac%lProdAlloc) CALL FreeGamProdMat(XsMac)
  XSMac%ng = ng
  XSMac%ngg = ngg
  XSMac%niso = niso
  CALL AllocGamProdMat(XSMac)
  IF (lIsoOut) THEN
    IF (XSMac%lIsoProdAlloc.AND.niso .NE. XSMac%niso) CALL FreeIsoGamProdMat(XsMac)
    CALL AllocIsoGamPRodMat(XsMAc)
  END IF
ENDIF
XSMac%GProdTot = 0.
XsMac%lexist = .FALSE.
MatTot => XSMac%GProdTot
lexist  => XsMac%lexist

! INITIALIZATION
IF(lIsoOut) THEN
  XSMac%IsoGProdFis = 0._8
  XSMac%IsoGProdRad = 0._8
  XSMac%IsoGProdInel = 0._8
  XSMac%IsoGProdNnel = 0._8
  XSMac%GProdTot = 0._8
  ! POINTING
  IsoProdFis => XSMac%IsoGProdFis
  IsoProdRad => XSMac%IsoGProdRad
  IsoProdInel => XSMac%IsoGProdInel
  IsoProdNnel => XSMac%IsoGProdNnel
ELSE
  ALLOCATE(IsoProdFis(niso, ng, igb:ige))
  ALLOCATE(IsoProdRad(niso, ng, igb:ige))
  ALLOCATE(IsoProdInel(niso, ng, igb:ige))
  ALLOCATE(IsoProdNnel(niso, ng, igb:ige))
  IsoProdFis = 0._8
  IsoProdRad = 0._8
  IsoProdInel= 0._8
  IsoProdNnel= 0._8
END IF
! ISOTOPEWISE (TEMP. DEPENDENCE + PNUM)
DO iso = 1, niso
  id = mapnucl(idiso(iso))
  IF (id .EQ. 0) THEN
    print *, 'SUBROUTINE GamProdMatrix -- No photon production data in library for', idiso(iso)
    CYCLE
  END IF
  lib => ldiso(id)
  IF (lib%ipp .EQ. 0)  CYCLE ! Nuclides whose photon production doesn't exist
  CALL XsTempInterpolation(id, lib, temp, wt1, wt2, it1, it2) ! TEMPERATURE DEPENDENCE
  DO imt = 1, 4
    IF (.NOT.lib%lphoton(imt)) CYCLE
    SELECT CASE(imt)
    CASE(1)
      IsoProd => IsoProdFis
    CASE(2)
      IsoProd => IsoProdRad
    CASE(3)
      IsoProd => IsoProdInel
    CASE(4)
      IsoProd => IsoProdNnel
    END SELECT
    nx1 = lib%ppm(imt)%iglow; nx2 = lib%ppm(imt)%igup
    IF(nx2.GE.igb .AND. nx1.LE.ige) THEN    ! Data existing region is beyond interest
      lexist(iso, imt) = .TRUE.
      tempgrpb = MAX(nx1, igb)
      tempgrpe = MIN(nx2, ige)
      DO igg = tempgrpb, tempgrpe
        ProdRange1 = (/lib%ppm(imt)%mat(igg,it1)%ib, lib%ppm(imt)%mat(igg,it1)%ie/)
        ProdRange2 = (/lib%ppm(imt)%mat(igg,it2)%ib, lib%ppm(imt)%mat(igg,it2)%ie/)
        nx1 = min(ProdRange1(1), ProdRange2(1)); nx2 = max(ProdRange1(2), ProdRange2(2))
        !nx1 = lib%ppm(imt)%mat(igg,it1)%ib; nx1 = lib%ppm(imt)%mat(igg,it2)%ie
        DO ig = nx1, nx2
          prod1 = 0; prod2 = 0;
          ind = (ig - ProdRange1(1)) * (ig - ProdRange1(2))
          IF(ind .LE. 0) prod1 = lib%ppm(imt)%mat(igg,it1)%from(ig)
          ind = (ig - ProdRange2(1)) * (ig - ProdRange2(2))
          IF(ind .LE. 0) prod2 = lib%ppm(imt)%mat(igg,it2)%from(ig)
          IsoProd(iso,ig,igg) = pnum(iso) * (wt1 * prod1 + wt2 * prod2)
        END DO
      END DO
    END IF
  END DO
END DO

IF (.NOT.any(lexist)) RETURN
IF (lres) THEN
  DO imt = 1, 2
    SELECT CASE(imt)
    CASE(1)
      IsoProd => IsoProdFis
      FresoReact => FresofIso
    CASE(2)
      IsoProd => IsoProdRad
      FresoReact => fresocapIso
!    CASE(3)
!      IsoProd => IsoProdInel
!      FresoReact => fresocapIso
    END SELECT
    DO iso = 1, niso
    IF (.NOT.lexist(iso, imt)) CYCLE
    DO igg = igb, ige
      DO ig = nofg+1, nofg+norg
          IsoProd(iso,ig,igg)  = IsoProd(iso,ig,igg) * FresoReact(iso, ig)
        END DO
      END DO
    END DO
  END DO
END IF
DO igg = igb, ige
  DO ig = 1, ng
    DO iso = 1, niso
      MatTot(ig, igg) = MatTot(ig, igg) + IsoProdFis(iso,ig,igg)
      MatTot(ig, igg) = MatTot(ig, igg) + IsoProdRad(iso,ig,igg)
      MatTot(ig, igg) = MatTot(ig, igg) + IsoProdInel(iso,ig,igg)
      MatTot(ig, igg) = MatTot(ig, igg) + IsoProdNnel(iso,ig,igg)
    END DO
  END DO
END DO

IF (.NOT.lIsoOut) THEN
  DEALLOCATE(IsoProdFis)
  DEALLOCATE(IsoProdRad)
  DEALLOCATE(IsoProdInel)
  DEALLOCATE(IsoProdNnel)
END IF
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
SUBROUTINE GetLocalQMat_FXR(XSMac, Fxr, igb, ige, ng, lisoout)
USE TYPEDEF,         ONLY : FXRINFO_TYPE, XsMac_Type
USE CORE_MOD,        ONLY : GroupInfo
IMPLICIT NONE
TYPE(FXRINFO_TYPE) :: FXR
TYPE(XsMac_Type) :: XSMac
INTEGER :: igb, ige, ng, ngg
INTEGER :: iResoGrpBeg, iResoGrpEnd, norg
LOGICAL :: lisoout

iResoGrpBeg = GroupInfo%nofg + 1
iResoGrpEnd = GroupInfo%nofg + GroupInfo%norg
norg = GroupInfo%norg

CALL GetLocalQMat_GEN(XSMac, FXR%temp, FXR%niso, FXR%idiso, FXR%pnum, igb, ige, ng, FXR%fresof, FXR%fresocapIso, FXR%fresoFIso, iResoGrpBeg, iResoGrpEnd, lisoout, FXR%lres)
                                                                                                      
END SUBROUTINE

!SUBROUTINE GetLocalQMat_CSP(XSMac, Fxr, igb, ige, ng, ngg, GroupInfo, lCrCspFtn)
!IMPLICIT NONE
!END SUBROUTINE
!
SUBROUTINE GetLocalQMat_GEN(XSMac, temp, niso, idiso, pnum, igb, ige, ng, fresof, fresocapIso, fresoFIso, iResoGrpBeg, iResoGrpEnd, lisoout, lres)
! NEUTRON INFORMATION
USE XSLIB_MOD,       ONLY : ldiso, nelthel, mapnucl, mapfis, itempmap, nelrhel, LIBDATA
! GAMMA INFORMATION
USE GamXSUtil,       ONLY : AllocIsoMacKERMA, AllocMacKERMA
USE CORE_MOD,        ONLY : GroupInfo
USE TYPEDEF,         ONLY : XsMac_Type

IMPLICIT NONE

! INPUT VARIABLES
INTEGER :: ng, igb, ige
TYPE(XsMac_Type) :: XSMac
INTEGER :: niso
INTEGER :: idiso(niso)
REAL :: pnum(niso)
REAL :: temp
INTEGER :: iResoGrpBeg, iResoGrpEnd, norg
REAL(4), POINTER :: fresof(:), fresocapIso(:, :), fresoFIso(:, :)
LOGICAL :: lisoout, lres

! POINTING VARIABLES
TYPE(LIBDATA), POINTER :: isodata

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
INTEGER :: imt
REAL, POINTER, DIMENSION(:, :) :: IsoMacKERMA_t, IsoMacKERMA_s, IsoMacKERMA_d, IsoMacKERMA_p, IsoMacKERMA_f
REAL, POINTER, DIMENSION(:, :) :: IsoMacDelkf

iResoGrpBeg = MAX(iResoGrpBeg, igb)
iResoGrpEnd = MIN(iResoGrpEnd, ige)

! ALLOCATION
IF (.NOT. XsMac%lKERMAAlloc) THEN
  XsMac%ng = ng
  XsMac%niso = nelthel
  CALL AllocMacKERMA(XsMac)
END IF
IF (lisoout .AND. .NOT. XsMac%lISOKERMAAlloc) THEN
  XsMac%ng = ng
  XsMac%niso = nelthel
  CALL AllocIsoMacKERMA(XsMac)
END IF

! Local Heat Deposition with KERMA calculation (NEUTRON INDUCED...)
IF (lisoout) THEN
  XsMac%IsoMacKERMA_t = 0.;XsMac%IsoMacKERMA_s = 0.;XsMac%IsoMacKERMA_d = 0.;XsMac%IsoMacKERMA_p = 0.;XsMac%IsoMacKERMA_f = 0.
  XsMac%IsoMacDelkf = 0.;
  IsoMacKERMA_t => XsMac%IsoMacKERMA_t
  IsoMacKERMA_s => XsMac%IsoMacKERMA_s
  IsoMacKERMA_d => XsMac%IsoMacKERMA_d
  IsoMacKERMA_p => XsMac%IsoMacKERMA_p
  IsoMacKERMA_f => XsMac%IsoMacKERMA_f
  IsoMacDelkf   => XsMac%IsoMacDelkf
ELSE
  ALLOCATE(IsoMacKERMA_t(niso, ng), IsoMacKERMA_s(niso, ng), IsoMacKERMA_d(niso, ng), IsoMacKERMA_p(niso, ng), IsoMacKERMA_f(niso, ng))
  ALLOCATE(IsoMacDelkf(niso, ng))
  IsoMacKERMA_t = 0.;IsoMacKERMA_s = 0.;IsoMacKERMA_d = 0.;IsoMacKERMA_p = 0.;IsoMacKERMA_f = 0.
  IsoMacDelkf = 0.
END IF


XsMac%MacKERMA_t = 0.;XsMac%MacKERMA_s = 0.;XsMac%MacKERMA_d = 0.;XsMac%MacKERMA_p = 0.;XsMac%MacKERMA_f = 0.
XsMac%MacDelkf = 0.
DO iso = 1,niso
  nid = MapNucl(idiso(iso)); isodata => ldiso(nid)
  CALL XsTempINterpolation(nid, isodata, temp, wt1, wt2, it1, it2)
  DO ig = igb, ige
    IsoMacKERMA_t(iso,ig) =  pnum(iso) * (wt1 * isodata%kerma_t(ig, it1) + wt2 * isodata%kerma_t(ig, it2))
    IsoMacKERMA_s(iso,ig) =  pnum(iso) * (wt1 * isodata%kerma_s(ig, it1) + wt2 * isodata%kerma_s(ig, it2))
    IsoMacKERMA_d(iso,ig) =  pnum(iso) * (wt1 * isodata%kerma_d(ig, it1) + wt2 * isodata%kerma_d(ig, it2))
    IsoMacKERMA_p(iso,ig) =  pnum(iso) * (wt1 * isodata%kerma_p(ig, it1) + wt2 * isodata%kerma_p(ig, it2))
    IsoMacKERMA_t(iso,ig) = IsoMacKERMA_t(iso,ig)-IsoMacKERMA_d(iso,ig)
  END DO
  IF(lres .AND. isodata%lreso) THEN
    DO ig = iResoGrpBeg, iResoGrpEnd
      IsoMacKERMA_d(iso,ig)=IsoMacKERMA_d(iso,ig)*fresocapIso(iso, ig)
    END DO
  END IF
  DO ig = igb, ige
    IsoMacKERMA_t(iso,ig) = IsoMacKERMA_t(iso,ig)+IsoMacKERMA_d(iso,ig)
  END DO
  IF(isodata%ifis.GT.0) THEN
    DO ig=igb, ige
      IsoMacKERMA_f(iso,ig) = pnum(iso) * (wt1 * isodata%kerma_f(ig, it1) + wt2 * isodata%kerma_f(ig, it2))
      IsoMacKERMA_t(iso,ig)= IsoMacKERMA_t(iso,ig) - IsoMacKERMA_f(iso,ig)
      IsoMacDelkf(iso,ig)   =  pnum(iso) * (wt1 * isodata%sigf(ig,it1) + wt2 * isodata%sigf(ig,it2))
    END DO
    IsoMacDelkf(iso, igb:ige) = IsoMacDelkf(iso, igb:ige) * (isodata%exkappa(3)+isodata%exkappa(5)+isodata%exkappa(6))
    IF(lres .AND. isodata%lreso) THEN
      DO ig = iResoGrpBeg, iResoGrpEnd
        IsoMacKERMA_f(iso,ig)=IsoMacKERMA_f(iso,ig)*fresofIso(iso, ig)
        IsoMacDelkf(iso,ig) = IsoMacDelkf(iso,ig)*fresofIso(iso, ig)
      END DO
    END IF
    DO ig = igb, ige
!#define delayedcontainedheat
#ifdef delayedcontainedheat
      IsoMacKERMA_t(iso,ig)= IsoMacKERMA_t(iso,ig) + IsoMacKERMA_f(iso,ig) + IsoMacDelkf(iso,ig)
#else
      IsoMacKERMA_t(iso,ig)= IsoMacKERMA_t(iso,ig) + IsoMacKERMA_f(iso,ig)! + IsoMacDelkf(iso,ig)
#endif
    END DO
    
  END IF
END DO

DO iso=1,niso
  nid = MapNucl(idiso(iso)); isodata => ldiso(nid)
  DO ig = igb, ige
!    XSMac%MacKERMA_t(ig)= XSMac%MacKERMA_t(ig) + IsoMacKERMA_t(iso,ig)
    XSMac%MacKERMA_s(ig)= XSMac%MacKERMA_s(ig) + IsoMacKERMA_s(iso,ig)
    XSMac%MacKERMA_d(ig)= XSMac%MacKERMA_d(ig) + IsoMacKERMA_d(iso,ig)
    XSMac%MacKERMA_p(ig)= XSMac%MacKERMA_p(ig) + IsoMacKERMA_p(iso,ig)
    XSMac%MacKERMA_T(ig) = XSMac%MacKERMA_T(ig) + IsoMacKERMA_T(iso,ig)
  END DO
  IF (isodata%ifis.EQ.0) CYCLE
    DO ig=igb, ige
      XsMac%MacKERMA_f(ig) = XsMac%MacKERMA_f(ig) + IsoMacKERMA_f(iso,ig)
      XsMac%MacDelkf(ig) = XsMac%MacDelkf(ig) + IsoMacDelkf(iso,ig)
    END DO
END DO

IF(.NOT.lisoout) THEN
  DEALLOCATE(IsoMacKERMA_t, IsoMacKERMA_s, IsoMacKERMA_d, IsoMacKERMA_p, IsoMacKERMA_f,IsoMacDelkf)
END IF


END SUBROUTINE

End Module
#endif