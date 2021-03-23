module CNTL
USE PARAM

LOGICAL :: prtscrn = .true.        !
LOGICAL :: lWriteRst = .false.     !Write Restart File

TYPE OutpCntl_Type
  SEQUENCE

  INTEGER :: IsoOutList(0:1000) = 0
  INTEGER :: FluxOutList(4, 0:100) = 0
  INTEGER :: VisMod = 0         !
  !PinXS Out
  LOGICAL :: lPinXsOut = .FALSE.
  INTEGER :: npinxsout = 0
  INTEGER, POINTER :: PinXsOutList(:)
  CHARACTER(256), POINTER :: PinXsFn(:)
  !Cusping Function Generation Output
  LOGICAL :: lCspGenOut = .FALSE.
  INTEGER :: CspGenBankId = 0
  !Effective regional Cross Section Out
  INTEGER :: nRegXsOut = 0
  INTEGER :: RegXsOutList(7, 1000) ! /[iz, ixa, iya, ix, iy, ir1, ir2], ireg /
  INTEGER :: IsoXsOutList(0:200, 1000) !
  LOGICAL :: RegXsOutASM(1000),RegXsOutPIN(1000),RegXsOutFXR(1000)

  !Effective scattering XS Out
  INTEGER :: nRegMAT0Out = 0
  INTEGER :: RegMAT0OutList(7, 1000) ! /[iz, ixa, iya, ix, iy, ir1, ir2], ireg /
  INTEGER :: IsoMAT0OutList(0:100, 1000) !  
  LOGICAL :: RegMAT0OutASM(1000),RegMAT0OutPIN(1000),RegMAT0OutFXR(1000)

  !Effective MG regional Cross Section Out
  INTEGER :: nRegMGXsOut = 0
  INTEGER :: RegMGXsOutList(7, 1000) ! /[iz, ixa, iya, ix, iy, ir1, ir2], ireg /
  INTEGER :: MGBdry(0:47, 1000) ! 
  INTEGER :: nMG(1000)
  
  !Condensed Flux Moments Out 
  INTEGER :: nRegPhimOut = 0
  INTEGER :: RegPhimOutList(7, 1000),PhimOrdOutList(0:21,1000)
  LOGICAL :: RegPhimOutASM(1000),RegPhimOutPIN(1000),RegPhimOutFXR(1000)

  !Effective regional Cross Section Out WHILE DEPLETION
  INTEGER :: nDeplRegXsOut = 0
  INTEGER :: DeplRegXsOutList(7, 1000) ! /[iz, ixa, iya, ix, iy, ir1, ir2], ireg /
  INTEGER :: DeplIsoXsOutList(0:100, 1000) ! 
  
  !spectral SPH factor out
  LOGICAL :: lSSPHout = .FALSE.

  !Binary Output
  LOGICAL :: lBOutp = .FALSE.
  LOGICAL :: lFlux_BOutp = .FALSE.
  LOGICAL :: lTh_BOutp = .FALSE.
  LOGICAL :: lDepl_BOutp = .FALSE.
  INTEGER :: FluxEdit = 2
  INTEGER :: DeplStepBOutp(0:1000) = -1
  INTEGER :: IsoBoutp(0:1000) = 0
END TYPE

type MCXS_type
  sequence
  integer :: iz,iasyx,iasyy,ix,iy,ifxr1,ifxr2,niso
  integer,pointer :: idiso(:)
  real,pointer :: isoxsa(:,:,:),isoxsnf(:,:,:)    
end type

TYPE nTracerCntl_Type
  SEQUENCE
  LOGICAL :: lXsLib = .TRUE.
  LOGICAL :: lBenchXs = .FALSE.
  LOGICAL :: lFXRLib = .FALSE.
  INTEGER :: XsOrder = -1
  LOGICAL :: lTrCorrection = .TRUE.
  LOGICAL :: lSSPH = .FALSE., lSSPHreg=.FALSE.
  LOGICAL :: lMLG=.TRUE., lRIF=.TRUE., lRST=.FALSE., lCAT=.TRUE., l4Lv=.TRUE.
  LOGICAL :: lED=.FALSE. ! Effective Dancoff
  LOGICAL :: lOldBon = .FALSE. ! Changhyun
  INTEGER :: nMLG=8, nMLGc=5

  INTEGER :: libTyp = 2
  LOGICAL :: lWriteRst = .FALSE.
  LOGICAL :: lWriteDeplRst = .FALSE.
  LOGICAL :: LpShf = .FALSE.
  LOGICAL :: lRstCal = .FALSE.
  LOGICAL :: l3Dim = .TRUE.
  INTEGER :: nDim = 3

  LOGICAL :: lrestrmt = .FALSE.

  INTEGER :: lProblem = lsseigv

  REAL :: PowerCore
  REAL :: PowerLevel = 1._8  ! Power Level
  REAL :: PowerFA      ! Nominal Assembly Power at the Full Power Condition (MW)
  REAL :: TempInlet    ! Inlet Coolant Temperature (C)
  REAL :: Pexit        ! Core Exit Coolant Pressure (MPa)
  REAL :: fMdotfa      ! Nominal Assembly Flow Rate (kg/sec)

  LOGICAL :: lCMFD = .TRUE.
  LOGICAL :: lGcCmfd = .FALSE.
  LOGICAL :: lGcCmfdGrp = .FALSE.
  INTEGER :: nGC
  INTEGER :: GcStructMod = 0  ! Predefined Group Structure : 0
                              ! Lower Energy Bound Value   : 1
                              ! Group Lower Bound Index    : 2
  REAL :: ELB(100)           ! Energy Lower Bound
  INTEGER :: GLB(100)        ! Group Lower Bound IndexT
  LOGICAL :: lSubPlane = .FALSE.  !SUbplane Scheme

  LOGICAL :: lLinSrc = .FALSE.    !linear Source
  LOGICAL :: lLinSrcCASMO = .FALSE.   !--- CNJ Edit : CASMO Linear Source
  LOGICAL :: lHybrid = .FALSE.        !--- CNJ Edit : CASMO Linear Source, Domain Decomposition
  LOGICAL :: lNodeMajor = .FALSE.     !--- CNJ Edit : Node Majors
  LOGICAL :: lDomainDcmp = .FALSE.    !--- CNJ Edit : Domain Decomposition
  !--- CNJ Edit : Angular Multigrid Ray Tracing
  LOGICAL :: lMultigrid = .FALSE.
  INTEGER :: MultigridLV = 4
  INTEGER :: gridNum, gridStr(1000)
  INTEGER :: AxSolver = lP3SENM   !Axial Sorver : P3 Default > r560d / 17/01/17
  LOGICAL :: lnonFuelpinFDM = .FALSE. !--- BYS Edit : non fuel pin FDM (skip Nodal) 
  LOGICAL :: lResetResErr = .FALSE.   !--- BYS Edit : Reset Residual Error check in CMFD
  LOGICAL :: lDhom = .FALSE.          !--- BYS Edit : true if hom. by D in axial
  
  INTEGER :: LkgSplitLv = 0
  
  LOGICAL :: lMOCUR = .FALSE.
  LOGICAL :: lOptUR = .TRUE.
  REAL :: UserUR = 0.75 
  
  LOGICAL :: lAxRefFDM = .FALSE.

  LOGICAL :: lFastMOC = .FALSE.
  INTEGER :: FastMOCLv = 0
  LOGICAL :: lAFSS = .FALSE.
  !Library Option

  LOGICAL :: lMATRA = .FALSE.
  LOGICAL :: lborontrack=.FALSE.
  LOGICAL :: lFeedBack = .FALSE.
  LOGICAL :: lSubGrpSweep = .FALSE.
  LOGICAL :: lSimpleTh = .FALSE.
  INTEGER :: ThCh_mod = 0                   ! T-H channel Mode 0 : Pin By Pin, 1 : Assembly Calculation, 2: 4 box,
  LOGICAL :: lFuelCondFDM = .FALSE.         ! Fuel Conduction FDM Solver
  LOGICAL :: LIHP = .FALSE.                 ! Internal Heat Profile
  LOGICAL :: lEffTemp = .FALSE.
  LOGICAL :: lUserDefTemp = .FALSE.

  LOGICAL :: lInitBoron  = .FALSE.          ! Boron Insertion
  LOGICAL :: lSearch = .FALSE.
  LOGICAL :: lBoronSearch = .FALSE.         ! Boron Search
  LOGICAL :: lCrSearch = .FALSE.            ! Control Rod Search
  !LOGICAL :: lB10Depl = .TRUE.             ! Boron-10 Depletion
  !INTEGER :: B10DeplMod = 1                ! Boron-10 Depletion   0 : TurnOff 1: Online  2 : Post Boron
  REAL :: BoronPPM = 0._8
  REAL :: Target_eigv = 1.0_8

  LOGICAL :: lDeplVarInit = .TRUE.   !Initialize Depletion Flux solution
  LOGICAL :: lXeDyn = .FALSE.
  LOGICAL :: lEqXe = .FALSE.          ! XE EQ
  LOGICAL :: lTrXe = .FALSE.           ! Tr Eq
  LOGICAL :: lReadXeEq = .FALSE.      ! Read Xe Eq
  LOGICAL :: lWriteXeEq = .FALSE.     ! Read Xe Eq

  LOGICAL :: lBsq = .FALSE.
  REAL :: Bsq = 0._8

  LOGICAL :: lBranch
  LOGICAL :: lCritSpec = .FALSE.
  LOGICAL :: lGC = .FALSE.    ! --- BYS edit : option for GroupConstantGeneration
  LOGICAL :: lPDQ_GCL = .FALSE.
  LOGICAL :: lASY_GCL = .FALSE.
  LOGICAL :: lASY_MIC = .TRUE.
  LOGICAL :: lPIN_GCL = .FALSE.
  LOGICAL :: lPIN_MIC = .FALSE.
  INTEGER :: GC_spec = 2      ! --- BYS edit : default as Solution Spectrum
  LOGICAL :: lPinwiseGC = .FALSE.
  LOGICAL :: lGCGapHom = .TRUE.
  LOGICAL :: lGCRST = .FALSE.
  INTEGER :: nGCgrp = 2       ! --- BYS edit : number of condensed group 2,4,8
  LOGICAL :: lEFTsearch = .FALSE.
  LOGICAL :: lEMTsearch = .FALSE.
  INTEGER :: nEFTpts = 15
  REAL    :: crit_EFT = 5 ! pcm
  REAL    :: P_EFT(50)
  INTEGER :: EFTUniCoolant = 0
  REAL    :: EFT_boron(100)
  INTEGER :: EFT_nboron = 1
  !--- CNJ Edit : Group-wise Albedo (BYS Request)
  LOGICAL :: lGroupAlbedo = .FALSE.
  LOGICAL :: lNodeTime = .FALSE.
  !
  LOGICAL :: lGap = .FALSE.
  LOGICAL :: lUsrGap = .FALSE.       !User Defined Gap Structure
  LOGICAL :: lhGap = .FALSE.  ! --- BYS edit : homogenized gap cell with fuel cell

  LOGICAL :: lRayGen =.TRUE.
  INTEGER :: RayGenOpt = 0
  
  LOGICAL :: lAutoBaffle = .FALSE.
  LOGICAL :: lAutoBarrel = .FALSE.
  LOGICAL :: lAutoRingStruct = .FALSE.
  LOGICAL :: lMiscStruct = .FALSE.
  !INTEGER :: BaffleMix(2)
  !INTEGER :: StructCell

  LOGICAL :: lDcpl = .FALSE.
  LOGICAL :: lDcplCal = .FALSE.
  LOGICAL :: lConstSrc = .TRUE.
  LOGICAL :: lDcplCmfdXsGen = .TRUE.
  LOGICAL :: lXsFtn = .TRUE.
  LOGICAL :: lXsFtnExt = .FALSE.
  INTEGER :: XsFtnMod = 1        ! 1 = Macro XS, 2 = Micro
  INTEGER :: DcplIterOpt(10)
  REAL :: ConstSrc = epsm4
  LOGICAL :: morecase = .TRUE.
  INTEGER :: TRSolver = 1 ! 1 = MOC, 2 = MC
  LOGICAL :: lAfSrc = .FALSE.
  ! RM EDIT
  LOGICAL :: lScatBd = .FALSE.  ! Scat boundary >> need to rename

!Control Rod Control Part
  LOGICAL :: lCrInfo = .FALSE.
  LOGICAL :: lCrIn = .FALSE.
  LOGICAL :: lCrCsp = .TRUE.
  LOGICAL :: lCrCspFtn = .TRUE.

  LOGICAL :: lFastScat1 = .TRUE.
  LOGICAL :: lScat1 = .FALSE.
  LOGICAL :: lGammaScat1 = .FALSE.   !--- CNJ Edit : Gamma Transport Pn Option
  LOGICAL :: lMacro = .FALSE.        !--- CNJ Edit : Pre-generated Macro Cross Section
  INTEGER :: ScatOd = 0
  INTEGER :: GammaScatOd = 0   !--- CNJ Edit : Gamma Transport Order

  INTEGER :: CaseNo = 0
  INTEGER :: CalNoId = 0
  !INTEGER :: IsoOutList(0:1000) = 0
  !INTEGER :: FluxOutList(4, 0:100) = 0
  TYPE(OutpCntl_Type) :: OutpCntl
  LOGICAL :: lFSRXS, lFSRPhi
  LOGICAL :: lBCRopt = .TRUE.
  LOGICAL :: lDetailOutput=.FALSE.
  
! PHOTON CALCULATION     !--- JSU EDIT 20170721
  LOGICAL :: lGamma = .false.
  ! Hexagonal calculation !--- KSC EDIT 20180521
  LOGICAL :: lHex = .FALSE.
  INTEGER :: nCP_er=0
  INTEGER :: nCP=0
  
  CHARACTER*512 :: MCXSfile
  LOGICAL :: lMCXS=.false.
  integer :: nmcxs
  type(MCXS_type),pointer :: mcxs(:)
  
  !--- JSR Edit : nTIG Restart
  LOGICAL :: lnTIGRst = .FALSE.
  LOGICAL :: lCooling = .FALSE. ! --- 180809
! Fast reactor...
  LOGICAL :: lisotxs = .FALSE.

END TYPE



TYPE(nTracerCntl_Type) nTracerCntl
TYPE(nTracerCntl_Type) DcplControl

#ifndef __GFORTRAN__
DATA nTracerCntl%OutpCntl%IsoOutList(1:18)/                                   &
                     90232,      92233,      92234,      92235,      92236,   &
                     92238,      93237,      93239,      94233,      94239,   &
                     94240,      94241,      94242,      95241,      95242,   &
                     95242,      96242,      96244                            /

DATA nTracerCntl%OutpCntl%IsoBoutp(1:18)/                                     &
                     90232,      92233,      92234,      92235,      92236,   &
                     92238,      93237,      93239,      94233,      94239,   &
                     94240,      94241,      94242,      95241,      95242,   &
                     95242,      96242,      96244                            /
#endif

CONTAINS
SUBROUTINE CopyCntlOption(RefCntl, TargetCntl)
USE PARAM
IMPLICIT NONE
TYPE(nTracerCntl_Type) RefCntl, TargetCntl
TargetCntl%lXsLib = RefCntl%lXsLib; TargetCntl%lBenchXs = RefCntl%lBenchXs
TargetCntl%lScat1 = RefCntl%lScat1; TargetCntl%libTyp = RefCntl%libTyp
TargetCntl%lWriteRst = RefCntl%lWriteRst; TargetCntl%l3Dim = RefCntl%l3Dim
TargetCntl%nDim = RefCntl%nDim
TargetCntl%lProblem = RefCntl%lProblem; TargetCntl%PowerFA = RefCntl%PowerFA
TargetCntl%TempInlet = RefCntl%TempInlet; TargetCntl%Pexit = RefCntl%Pexit
TargetCntl%fMdotfa = RefCntl%fMdotfa; TargetCntl%lCMFD = RefCntl%lCMFD
TargetCntl%lSubPlane = RefCntl%lSubPlane; TargetCntl%lDcpl = RefCntl%lDcpl
TargetCntl%lLinSrc = RefCntl%lLinSrc; TargetCntl%AxSolver = RefCntl%AxSolver
TargetCntl%lFeedBack = RefCntl%lFeedBack; TargetCntl%lSubGrpSweep = RefCntl%lSubGrpSweep
TargetCntl%lInitBoron = RefCntl%lInitBoron; TargetCntl%BoronPPM = RefCntl%BoronPPM
TargetCntl%lBsq = RefCntl%lBsq; TargetCntl%Bsq = RefCntl%Bsq
TargetCntl%lDcplCal = RefCntl%lDcplCal; TargetCntl%lxsFtn = RefCntl%lxsFtn
TargetCntl%XsFtnMod = RefCntl%XsFtnMod
END SUBROUTINE
END module