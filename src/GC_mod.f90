MODULE GC_mod
USE PARAM
USE TYPEDEF,         ONLY : CoreInfo_Type,FmInfo_Type, THInfo_Type, CMInfo_Type,Cell_Type,Pin_Type, &
                            FxrInfo_Type, PinXs_Type, XsMac_Type, GroupInfo_Type,PE_TYPE, Powerdist_TYPE, &
                            Asy_Type, Asyinfo_type
USE BasicOperation, ONLY : CP_VA, CP_CA, MULTI_VA, MULTI_CA
USE files,           ONLY : caseid
USE XSLIB_MOD
USE CNTL,         ONLY : nTracerCntl_Type
USE CritSpec_mod,     ONLY : GetDiffusionCoeff
USE Core_mod,         ONLY : eigv
IMPLICIT NONE

TYPE(AsyInfo_Type), POINTER :: AsyInfo(:)
TYPE(Asy_Type), POINTER :: Asy(:)
TYPE(Pin_Type), POINTER :: Pin(:)
TYPE(Cell_Type), POINTER :: Cell(:)
TYPE(LIBDATA), POINTER :: isodata
TYPE(FxrInfo_Type), POINTER :: myFXR
REAL :: TempRef   !THinfo%RefFuelTemp(:)

!--- Lists of Constants
INTEGER, parameter :: isosize=500 ! maximum # of isotope list size( depends on the memory)
INTEGER :: nxs = 6                ! # of XS kinds (!0 for total, 1: tr, 2: abs, 3: rmv, 4: fis, 5: nufis, 6: kfis)
INTEGER :: ngrp = 2
INTEGER :: nExtPol = 3
INTEGER :: nglist(5)
REAL(4) :: GrpLBList(16,5)
DATA nglist /2, 4, 8, 16, 1/
DATA GrpLBList / 6.2506E-01,        0._8,        0._8,        0._8,        0._8,        0._8,        0._8,        0._8, &
                       0._8,        0._8,        0._8,        0._8,        0._8,        0._8,        0._8,        0._8, &
                9.1188E+3_8, 3.9279E+0_8, 6.2506E-1_8,        0._8,        0._8,        0._8,        0._8,        0._8, &
                       0._8,        0._8,        0._8,        0._8,        0._8,        0._8,        0._8,        0._8, &
                2.2313E+6_8, 8.2085E+5_8, 9.1188E+3_8, 1.3007E+2_8, 3.9279E+0_8, 6.2506E-1_8, 1.4572E-1_8,        0._8, &
                       0._8,        0._8,        0._8,        0._8,        0._8,        0._8,        0._8,        0._8, &
                2.2313E+6_8, 8.2085E+5_8, 9.1188E+3_8, 1.3007E+2_8, 1.3710E+1_8, 6.4760E+0_8, 3.9279E+0_8, 1.4574E+0_8, &
                1.0722E+0_8, 6.2506E-1_8, 3.5767E-1_8, 1.4572E-1_8, 8.1968E-2_8, 4.2755E-2_8, 1.2396E-2_8,        0._8, &
                       0._8,        0._8,        0._8,        0._8,        0._8,        0._8,        0._8,        0._8, &
                       0._8,        0._8,        0._8,        0._8,        0._8,        0._8,        0._8,        0._8  /
!--- End of Lists of Constants

!--- Flux variables
REAL, POINTER :: Phi(:), PhiVol(:), PinPhiVol(:,:)          !Phi(ng), PhiVol(ng), PinPhiVol(ng,nxy)
REAL, POINTER :: GrpIdx(:), PhiG(:) , PinPhiVolG(:,:,:)       !GrpIdx(ng), PhiG(ngrp), PinPhiVolG(ngrp,nx,ny)
REAL, POINTER :: FsrPhi(:,:), FsrVol(:), FsrPhiVol(:,:)  !FsrPhi(nCoreFsr,ng), FsrVol(nCoreFsr), FsrPhiVol(nCoreFsr,ng)
REAL :: fphi,tphi, rphi
REAL, POINTER :: PhiC(:, :, :)
!--- end of Flux variables

!REAL, POINTER :: MacXsSm(:,:)
!REAL, POINTER :: MacXsSmP1(:,:)

!--- Isotope variables
INTEGER :: nisotot
REAL :: isoNumden(isosize), h2oNumden, uo2Numden ! O+H, O-H number density  ! 500 > 100  14/02/27
INTEGER :: isoName(isosize), isoList(isosize,2)
LOGICAL :: checkiso(isosize)
REAL :: myn2n, u238n2n, u238phi

!--- Isotope-wise XS variables
REAL,POINTER :: isoMicXs(:,:,:),isoMicSm(:,:,:), isoMicXs2g(:,:,:), isoMicSm2g(:,:,:) !tr, a, r, f, nu, k
REAL,POINTER :: isoMacXs(:,:,:),isoMacSm(:,:,:), isoMacXs2g(:,:,:), isoMacSm2g(:,:,:) !D, a, r, f, nf, kf
!isoMicXs(nisotot,0:nxs,ng), isoMicSm(nisotot,ng,ng),isoMicXs2g(nisotot,0:nxs,2),isoMicSm2g(nisotot,2,2)
!isoMacXs(0:nisotot,0:nxs,ng), isoMacSm(0:nisotot,ng,ng),isoMacXs2g(0:nisotot,0:nxs,2),isoMacSm2g(nisotot,2,2)
REAL, POINTER :: IsoMacXsSm(:,:,:)  !IsoMacXsSm(isosize,ng,ng)
!-P1 variable- ! BYS EDIT 15/12/15 P1SM
REAL, POINTER :: isoMicSmP1(:,:,:)  , IsoMacSmP1(:,:,:)  !isoMicSmP1(nisotot,ng,ng), isoMacSmP1(0:nisotot,ng,ng)
REAL, POINTER :: isoMicSm2gP1(:,:,:), IsoMacSm2gP1(:,:,:)  !isoMicSmP1(nisotot,2,2), isoMacSmP1(0:nisotot,2,2)
REAL, POINTER :: IsoMacXsSmP1(:,:,:) !IsoMacXsSmP1(isosize,ng,ng)

!--- FXR, FSR, Macro values
TYPE(PinXS_Type), POINTER :: PinXS(:, :)
!--- end of FXR, FSR variables


!--- GC values
REAL :: Bsq, Msq, kinf, keff
REAL, POINTER :: buexp(:,:), pwr(:,:)
TYPE(PowerDist_Type) :: PowerDist
REAL :: NormFac, PFac

!--- B1 critical correction
LOGICAL :: lB1
REAL, POINTER :: f_ratio(:), phicrit(:)
REAL, POINTER :: Dng(:), XstrB1(:), Xstr(:), Xst(:)


!--- ADF variables
REAL, POINTER :: ADF(:,:), ADF2g(:,:), avgADF2g(:), ADF2Gsurf(:,:) !ADF(4,ng), ADF2g(4,ngrp), avgADF2g(ngrp)
REAL, POINTER :: avgBPhiG(:), conPhiG(:,:) !avgBPhiG(ngrp)  !conPhiG(3, ngrp) center/ avg/ ig
REAL, POINTER :: SurfJG(:,:,:) ! SurfJG(inout2, nbd4, ngrp)
REAL, POINTER :: SurfJ(:,:,:) ! SurfJG(inout2, nbd4, ng)
REAL, POINTER :: SurfPhiG(:,:) ! SurfPhiG(nbd4, ngrp)
REAL, POINTER :: ChiG(:)


!--- local index
REAL :: hz
INTEGER :: nCoreFsr, nCoreFxr, nCorexy
INTEGER :: nx, ny, nz, nxy, nxya
INTEGER :: ipin, iasy, iz, iasytype
INTEGER :: npinlist
INTEGER :: xbg, xed, yst, yed

!--- Logical variables
LOGICAL :: lSA, lCB !lSingle Assembly, lCheckerBoard
LOGICAL :: lZdir
LOGICAL :: lPinwise, lPin
LOGICAL :: lGapHom

LOGICAL :: lPDQ_GCL = .FALSE.
LOGICAL :: lASY_GCL = .FALSE.
LOGICAL :: lASY_MIC = .TRUE.
LOGICAL :: lPIN_GCL = .FALSE.
LOGICAL :: lPIN_MIC = .FALSE.
!--- end of logical variables

INTEGER :: LFP_ID

END MODULE

MODULE GCpin_mod
IMPLICIT NONE
    !REAL, POINTER :: pinXSr(:,:), pinXSnf(:,:), pinXSss(:,:,:), pinChi(:,:)
    REAL, POINTER :: pinRmv(:,:), pinFis(:,:), pinSS(:,:,:), pinFlux(:,:), pinR(:,:), pinJ(:,:) !ipin,ng
    REAL, POINTER :: pinSout(:,:), pinSin(:,:), pinAbs(:,:)
    REAL, POINTER :: pinSoutG(:,:), pinSinG(:,:)
    REAL, POINTER :: pinJdir(:,:,:)

    !REAL, POINTER :: pinRmvR(:,:), pinFisR(:,:), pinJR(:,:)
    !REAL, POINTER :: Srcr(:,:), Lossr(:,:), pinres(:,:)
    !REAL, POINTER :: fsrphir(:,:), fsrxs(:,:,:)
END MODULE
